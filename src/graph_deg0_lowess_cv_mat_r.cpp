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

    // predictions (slot 0): list<REALSXP vector>
    {
        const int L = (int) result.predictions.size();
        SEXP r_predictions = PROTECT(Rf_allocVector(VECSXP, L));
        for (int j = 0; j < L; ++j) {
            const int n = (int) result.predictions[(size_t)j].size();
            SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
            if (n > 0) {
                double* p = REAL(v);
                const auto& src = result.predictions[(size_t)j];
                std::copy(src.begin(), src.end(), p);
            }
            SET_VECTOR_ELT(r_predictions, j, v);
            UNPROTECT(1);  // v
        }
        SET_VECTOR_ELT(r_result, 0, r_predictions);
        UNPROTECT(1);  // r_predictions
    }

    // bw_predictions (slot 1): list<list<REALSXP vector>> if requested, else empty list
    {
        if (with_bw_predictions && !result.bw_predictions.empty()) {
            const int L = (int) result.bw_predictions.size();  // responses
            SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, L));

            for (int j = 0; j < L; ++j) {
                const int B = (int) result.bw_predictions[(size_t)j].size();  // bandwidths
                SEXP r_bw_preds_j = PROTECT(Rf_allocVector(VECSXP, B));

                for (int bw_idx = 0; bw_idx < B; ++bw_idx) {
                    const auto& src = result.bw_predictions[(size_t)j][(size_t)bw_idx];
                    const int n = (int) src.size();
                    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
                    if (n > 0) {
                        double* p = REAL(v);
                        std::copy(src.begin(), src.end(), p);
                    }
                    SET_VECTOR_ELT(r_bw_preds_j, bw_idx, v);
                    UNPROTECT(1);  // v
                }

                SET_VECTOR_ELT(r_bw_predictions, j, r_bw_preds_j);
                UNPROTECT(1);  // r_bw_preds_j
            }

            SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
            UNPROTECT(1);  // r_bw_predictions
        } else {
            SEXP empty = PROTECT(Rf_allocVector(VECSXP, 0));
            SET_VECTOR_ELT(r_result, 1, empty);
            UNPROTECT(1);  // empty
        }
    }

    // bw_errors (slot 2): list<REALSXP vector>
    {
        const int L = (int) result.bw_errors.size();
        SEXP r_bw_errors = PROTECT(Rf_allocVector(VECSXP, L));

        for (int j = 0; j < L; ++j) {
            const auto& src = result.bw_errors[(size_t)j];
            const int n = (int) src.size();
            SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
            if (n > 0) {
                double* p = REAL(v);
                std::copy(src.begin(), src.end(), p);
            }
            SET_VECTOR_ELT(r_bw_errors, j, v);
            UNPROTECT(1);  // v
        }

        SET_VECTOR_ELT(r_result, 2, r_bw_errors);
        UNPROTECT(1);  // r_bw_errors
    }

    // bws (slot 3): REALSXP vector
    {
        const int n = (int) result.bws.size();
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        if (n > 0) {
            double* p = REAL(s);
            std::copy(result.bws.begin(), result.bws.end(), p);
        }
        SET_VECTOR_ELT(r_result, 3, s);
        UNPROTECT(1);  // s
    }

    // opt_bws (slot 4): REALSXP vector
    {
        const int n = (int) result.opt_bws.size();
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        if (n > 0) {
            double* p = REAL(s);
            std::copy(result.opt_bws.begin(), result.opt_bws.end(), p);
        }
        SET_VECTOR_ELT(r_result, 4, s);
        UNPROTECT(1);  // s
    }

    // opt_bw_idxs (slot 5): INTSXP vector (1-based)
    {
        const int L = (int) result.opt_bw_idxs.size();
        SEXP r_opt_bw_idxs = PROTECT(Rf_allocVector(INTSXP, L));
        if (L > 0) {
            int* ip = INTEGER(r_opt_bw_idxs);
            for (int j = 0; j < L; ++j) ip[j] = (int) result.opt_bw_idxs[(size_t)j] + 1;
        }
        SET_VECTOR_ELT(r_result, 5, r_opt_bw_idxs);
        UNPROTECT(1);  // r_opt_bw_idxs
    }

    UNPROTECT(1);  // r_result
    return r_result;
}
