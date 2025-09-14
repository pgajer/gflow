/**
 * @brief Fixed version of graph_deg0_lowess_cv.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_graph_deg0_lowess_cv
 * 
 * Issues fixed:
 * 1. S_graph_deg0_lowess_cv (line 559): UNPROTECT(variable) - protect_count
 * 2. Lambda function captures and modifies protect_count
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constant
 * 2. Removed lambda that captures and modifies protect_count
 * 3. Used container-first pattern consistently
 */

#include <vector>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers
#include "set_wgraph.hpp"
#include "graph_deg0_lowess_cv.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations
struct graph_deg0_lowess_cv_t {
    std::vector<double> predictions;
    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> bw_errors;
    std::vector<double> bws;
    double opt_bw;
    size_t opt_bw_idx;
};

class set_wgraph_t {
public:
    set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                 const std::vector<std::vector<double>>& weight_list);
    
    graph_deg0_lowess_cv_t graph_deg0_lowess_cv(
        const std::vector<double>& y,
        double min_bw_factor,
        double max_bw_factor,
        size_t n_bws,
        bool log_grid,
        size_t kernel_type,
        double dist_normalization_factor,
        bool use_uniform_weights,
        size_t n_folds,
        bool with_bw_predictions,
        double precision,
        bool verbose);
};

extern "C" {

/**
 * Fixed version of S_graph_deg0_lowess_cv
 * Fixes: UNPROTECT(variable) at line 559, lambda capture issue
 * Solution: Use container-first pattern, remove lambda, use literal UNPROTECT
 */
SEXP S_graph_deg0_lowess_cv(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
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

    // Convert numeric vector
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

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
    graph_deg0_lowess_cv_t result = graph.graph_deg0_lowess_cv(
        y,
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

    // Create the return list (container-first pattern)
    const int n_elements = 6;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));

    // Set names
    {
        SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
        SET_STRING_ELT(r_result_names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(r_result_names, 1, Rf_mkChar("bw_predictions"));
        SET_STRING_ELT(r_result_names, 2, Rf_mkChar("bw_errors"));
        SET_STRING_ELT(r_result_names, 3, Rf_mkChar("bws"));
        SET_STRING_ELT(r_result_names, 4, Rf_mkChar("opt_bw"));
        SET_STRING_ELT(r_result_names, 5, Rf_mkChar("opt_bw_idx"));
        Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
        UNPROTECT(1); // r_result_names
    }

    // 0: predictions
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, result.predictions.size()));
        double* ptr = REAL(r_vec);
        std::copy(result.predictions.begin(), result.predictions.end(), ptr);
        SET_VECTOR_ELT(r_result, 0, r_vec);
        UNPROTECT(1);
    }

    // 1: bw_predictions (if requested)
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions.size()));
        
        for (size_t i = 0; i < result.bw_predictions.size(); i++) {
            SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, result.bw_predictions[i].size()));
            double* ptr = REAL(r_vec);
            std::copy(result.bw_predictions[i].begin(), result.bw_predictions[i].end(), ptr);
            SET_VECTOR_ELT(r_bw_predictions, i, r_vec);
            UNPROTECT(1); // r_vec
        }
        
        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
        UNPROTECT(1); // r_bw_predictions
    } else {
        SET_VECTOR_ELT(r_result, 1, Rf_allocVector(VECSXP, 0)); // Empty list
    }

    // 2: bw_errors
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, result.bw_errors.size()));
        double* ptr = REAL(r_vec);
        std::copy(result.bw_errors.begin(), result.bw_errors.end(), ptr);
        SET_VECTOR_ELT(r_result, 2, r_vec);
        UNPROTECT(1);
    }

    // 3: bws
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, result.bws.size()));
        double* ptr = REAL(r_vec);
        std::copy(result.bws.begin(), result.bws.end(), ptr);
        SET_VECTOR_ELT(r_result, 3, r_vec);
        UNPROTECT(1);
    }

    // 4: opt_bw
    {
        SEXP r_opt_bw = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(r_opt_bw)[0] = result.opt_bw;
        SET_VECTOR_ELT(r_result, 4, r_opt_bw);
        UNPROTECT(1);
    }

    // 5: opt_bw_idx
    {
        SEXP r_opt_bw_idx = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_opt_bw_idx)[0] = result.opt_bw_idx + 1; // Convert to 1-based for R
        SET_VECTOR_ELT(r_result, 5, r_opt_bw_idx);
        UNPROTECT(1);
    }

    UNPROTECT(1); // r_result
    return r_result;
}

} // extern "C"