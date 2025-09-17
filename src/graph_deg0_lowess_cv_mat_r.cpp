#include "graph_deg0_lowess_cv_mat.hpp" // For graph_deg0_lowess_cv_mat_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++

#include <vector>                   // For std::vector

#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions

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
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert Y matrix/list from R
    std::vector<std::vector<double>> Y;

    // Check if s_Y is a matrix or a list
    if (Rf_isMatrix(s_Y)) {
        // Handle matrix format (convert column-major R matrix to row-major C++ vector of vectors)
        SEXP s_dim = PROTECT(Rf_getAttrib(s_Y, R_DimSymbol));
        if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) < 1) {
            UNPROTECT(1);
            Rf_error("Y must be a matrix with a valid integer 'dim' attribute.");
        }
        const int nrow = INTEGER(s_dim)[0];
        const int ncol = INTEGER(s_dim)[1];
        UNPROTECT(1); // s_dim

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

    // Scalars (defensive extraction)
    const double min_bw_factor            = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor            = Rf_asReal(s_max_bw_factor);
    const size_t n_bws                    = (size_t)Rf_asInteger(s_n_bws);
    const bool   log_grid                 = (Rf_asLogical(s_log_grid) == TRUE);
    const size_t kernel_type              = (size_t)Rf_asInteger(s_kernel_type);
    const double dist_norm_factor         = Rf_asReal(s_dist_normalization_factor);
    const bool   use_uniform_weights      = (Rf_asLogical(s_use_uniform_weights) == TRUE);
    const size_t n_folds                  = (size_t)Rf_asInteger(s_n_folds);
    const bool   with_bw_predictions      = (Rf_asLogical(s_with_bw_predictions) == TRUE);
    const double precision                = Rf_asReal(s_precision);
    const bool   verbose                  = (Rf_asLogical(s_verbose) == TRUE);

    // Create graph and run algorithm
    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);
    graph_deg0_lowess_cv_mat_t result = graph.graph_deg0_lowess_cv_mat(
        Y,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        log_grid,
        kernel_type,
        dist_norm_factor,
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

    // Assemble R result
    const int n_elements = 6;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));

    // names
    {
        SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
        for (int i = 0; i < n_elements; i++) {
            SET_STRING_ELT(r_result_names, i, Rf_mkChar(names[i]));
        }
        Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
        UNPROTECT(1);
    }


    // Helper function to convert vector to SEXP - user has to unprotect the result
    auto vec_to_sexp = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set predictions (list of numeric vectors, one per response variable)
    {
        SEXP r_predictions = PROTECT(Rf_allocVector(VECSXP, result.predictions.size()));
        for (size_t j = 0; j < result.predictions.size(); j++) {
            SEXP s = vec_to_sexp(result.predictions[j]);
            SET_VECTOR_ELT(r_predictions, j, s);
            UNPROTECT(1); // s - protected in vec_to_sexp()
        }
        SET_VECTOR_ELT(r_result, 0, r_predictions);
        UNPROTECT(1);
    }

    // Set bw_predictions (if requested) - nested list structure
    {
        if (with_bw_predictions && !result.bw_predictions.empty()) {
            // Create outer list (one element per response variable)
            SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions.size()));

            for (size_t j = 0; j < result.bw_predictions.size(); j++) {
                // Create inner list (one element per bandwidth)
                SEXP r_bw_preds_j = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions[j].size()));

                for (size_t bw_idx = 0; bw_idx < result.bw_predictions[j].size(); bw_idx++) {
                    SEXP s = vec_to_sexp(result.bw_predictions[j][bw_idx]);
                    SET_VECTOR_ELT(r_bw_preds_j, bw_idx, s);
                    UNPROTECT(1); // s - protected in vec_to_sexp()
                }

                SET_VECTOR_ELT(r_bw_predictions, j, r_bw_preds_j);
                UNPROTECT(1); // r_bw_preds_j
            }

            SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
            UNPROTECT(1); // r_bw_predictions
        } else {
            SET_VECTOR_ELT(r_result, 1, Rf_allocVector(VECSXP, 0)); // Empty list
        }
    }

    // Set bw_errors (list of numeric vectors, one per response variable)
    {
        SEXP r_bw_errors = PROTECT(Rf_allocVector(VECSXP, result.bw_errors.size()));

        for (size_t j = 0; j < result.bw_errors.size(); j++) {
            SEXP s = vec_to_sexp(result.bw_errors[j]);
            SET_VECTOR_ELT(r_bw_errors, j, s);
            UNPROTECT(1); // s
        }
        SET_VECTOR_ELT(r_result, 2, r_bw_errors);
        UNPROTECT(1); // r_bw_errors
    }

    // Set bws
    {
        SEXP s = vec_to_sexp(result.bws);
        SET_VECTOR_ELT(r_result, 3, s);
        UNPROTECT(1); // s
    }

    // Set opt_bws
    {
        SEXP s = vec_to_sexp(result.opt_bws);
        SET_VECTOR_ELT(r_result, 4, s);
        UNPROTECT(1); // s
    }

    // Set opt_bw_idxs (convert to 1-based for R)
    {
        SEXP r_opt_bw_idxs = PROTECT(Rf_allocVector(INTSXP, result.opt_bw_idxs.size()));

        for (size_t j = 0; j < result.opt_bw_idxs.size(); j++) {
            INTEGER(r_opt_bw_idxs)[j] = result.opt_bw_idxs[j] + 1; // Convert to 1-based for R
        }
        SET_VECTOR_ELT(r_result, 5, r_opt_bw_idxs);
        UNPROTECT(1); // r_opt_bw_idxs
    }

    UNPROTECT(1); // r_result
    return r_result;
}
