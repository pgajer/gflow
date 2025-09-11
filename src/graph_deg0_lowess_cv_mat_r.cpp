#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions

// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

#include <vector>                   // For std::vector

#include "graph_deg0_lowess_cv_mat.hpp" // For graph_deg0_lowess_cv_mat_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++

/**
 * @brief R interface for matrix version of graph degree 0 LOWESS with CV
 */
extern "C" SEXP S_graph_deg0_lowess_cv_mat(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_use_uniform_weights,
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose) {

    // Convert input parameters
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert Y matrix/list from R
    std::vector<std::vector<double>> Y;

    // Check if s_Y is a matrix or a list
    if (Rf_isMatrix(s_Y)) {
        // Handle matrix format (convert column-major R matrix to row-major C++ vector of vectors)
        SEXP s_Y_dim = Rf_getAttrib(s_Y, R_DimSymbol);
        int nrow = INTEGER(s_Y_dim)[0];
        int ncol = INTEGER(s_Y_dim)[1];

        // Each column of the matrix becomes a response variable vector
        Y.resize(ncol);
        double* Y_ptr = REAL(s_Y);

        for (int j = 0; j < ncol; j++) {
            Y[j].resize(nrow);
            for (int i = 0; i < nrow; i++) {
                // Access column-major matrix: Y_ptr[i + j*nrow]
                Y[j][i] = Y_ptr[i + j*nrow];
            }
        }
    } else if (Rf_isNewList(s_Y)) {
        // Handle list format
        int n_cols = LENGTH(s_Y);
        Y.resize(n_cols);

        for (int j = 0; j < n_cols; j++) {
            SEXP s_y_j = VECTOR_ELT(s_Y, j);
            if (!Rf_isReal(s_y_j)) {
                Rf_error("All elements of Y must be numeric vectors");
            }

            int n_rows = LENGTH(s_y_j);
            Y[j].resize(n_rows);
            double* y_ptr = REAL(s_y_j);

            for (int i = 0; i < n_rows; i++) {
                Y[j][i] = y_ptr[i];
            }
        }
    } else {
        Rf_error("Y must be a numeric matrix or a list of numeric vectors");
    }

    // Get scalar parameters
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];
    size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
    bool log_grid = (LOGICAL(s_log_grid)[0] == 1);
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool use_uniform_weights = (LOGICAL(s_use_uniform_weights)[0] == 1);
    size_t n_folds = (size_t)INTEGER(s_n_folds)[0];
    bool with_bw_predictions = (LOGICAL(s_with_bw_predictions)[0] == 1);
    double precision = REAL(s_precision)[0];
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    // Create graph and run algorithm
    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);
    graph_deg0_lowess_cv_mat_t result = graph.graph_deg0_lowess_cv_mat(
        Y,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        log_grid,
        kernel_type,
        dist_normalization_factor,
        use_uniform_weights,
        n_folds,
        with_bw_predictions,
        precision,
        verbose
    );

    // Create the return list
    const char* names[] = {
        "predictions",
        "bw_predictions",
        "bw_errors",
        "bws",
        "opt_bws",
        "opt_bw_idxs",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Track protection count
    int protect_count = 0;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    protect_count++;

    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    protect_count++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);

    // Helper function to convert vector to SEXP
    auto vec_to_sexp = [&protect_count](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        protect_count++;
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set predictions (list of numeric vectors, one per response variable)
    SEXP r_predictions = PROTECT(Rf_allocVector(VECSXP, result.predictions.size()));
    protect_count++;
    for (size_t j = 0; j < result.predictions.size(); j++) {
        SET_VECTOR_ELT(r_predictions, j, vec_to_sexp(result.predictions[j]));
    }
    SET_VECTOR_ELT(r_result, 0, r_predictions);

    // Set bw_predictions (if requested) - nested list structure
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        // Create outer list (one element per response variable)
        SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions.size()));
        protect_count++;
        
        for (size_t j = 0; j < result.bw_predictions.size(); j++) {
            // Create inner list (one element per bandwidth)
            SEXP r_bw_preds_j = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions[j].size()));
            protect_count++;
            
            for (size_t bw_idx = 0; bw_idx < result.bw_predictions[j].size(); bw_idx++) {
                SET_VECTOR_ELT(r_bw_preds_j, bw_idx, vec_to_sexp(result.bw_predictions[j][bw_idx]));
            }
            
            SET_VECTOR_ELT(r_bw_predictions, j, r_bw_preds_j);
        }
        
        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
    } else {
        SET_VECTOR_ELT(r_result, 1, Rf_allocVector(VECSXP, 0)); // Empty list
    }

    // Set bw_errors (list of numeric vectors, one per response variable)
    SEXP r_bw_errors = PROTECT(Rf_allocVector(VECSXP, result.bw_errors.size()));
    protect_count++;
    for (size_t j = 0; j < result.bw_errors.size(); j++) {
        SET_VECTOR_ELT(r_bw_errors, j, vec_to_sexp(result.bw_errors[j]));
    }
    SET_VECTOR_ELT(r_result, 2, r_bw_errors);

    // Set bws
    SET_VECTOR_ELT(r_result, 3, vec_to_sexp(result.bws));

    // Set opt_bws
    SET_VECTOR_ELT(r_result, 4, vec_to_sexp(result.opt_bws));

    // Set opt_bw_idxs (convert to 1-based for R)
    SEXP r_opt_bw_idxs = PROTECT(Rf_allocVector(INTSXP, result.opt_bw_idxs.size()));
    protect_count++;
    for (size_t j = 0; j < result.opt_bw_idxs.size(); j++) {
        INTEGER(r_opt_bw_idxs)[j] = result.opt_bw_idxs[j] + 1; // Convert to 1-based for R
    }
    SET_VECTOR_ELT(r_result, 5, r_opt_bw_idxs);

    UNPROTECT(protect_count);
    return r_result;
}
