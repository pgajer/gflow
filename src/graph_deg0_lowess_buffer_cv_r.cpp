#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length
#undef eval

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <mutex>                    // For std::mutex
#include <execution>                // For std::execution::par_unseq
#include <atomic>                   // For std::atomic
#include <thread>                   // For std::thread::hardware_concurrenyc

#include "graph_deg0_lowess_buffer_cv.hpp" // For graph_deg0_lowess_buffer_cv_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

extern "C" {
	SEXP S_graph_deg0_lowess_buffer_cv(
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
    SEXP s_buffer_hops,
    SEXP s_auto_buffer_hops,
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose);
}

/**
 * @brief R interface for graph degree 0 LOWESS with buffer zone cross-validation
 *
 * @details This function serves as the interface between R and C++ implementation
 * of the buffer zone cross-validation method for graph-based LOWESS. It handles
 * the conversion of R objects to C++ types, calls the C++ implementation, and
 * returns the results in R-compatible format.
 *
 * @param s_adj_list R list of integer vectors representing adjacency lists
 * @param s_weight_list R list of numeric vectors with edge weights
 * @param s_y R numeric vector of response values
 * @param s_min_bw_factor R numeric scalar for minimum bandwidth factor
 * @param s_max_bw_factor R numeric scalar for maximum bandwidth factor
 * @param s_n_bws R integer scalar for number of bandwidths to test
 * @param s_log_grid R logical scalar for log-spaced bandwidth grid
 * @param s_kernel_type R integer scalar for kernel function type
 * @param s_dist_normalization_factor R numeric scalar for distance normalization
 * @param s_buffer_hops R integer scalar for buffer zone hop distance
 * @param s_auto_buffer_hops R logical scalar to auto-determine buffer size
 * @param s_n_folds R integer scalar for number of cross-validation folds
 * @param s_with_bw_predictions R logical scalar to store bandwidth predictions
 * @param s_precision R numeric scalar for precision
 * @param s_verbose R logical scalar for verbose output
 *
 * @return SEXP R list containing cross-validation results
 */
SEXP S_graph_deg0_lowess_buffer_cv(
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
    SEXP s_buffer_hops,        // New parameter
    SEXP s_auto_buffer_hops,   // New parameter
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
    bool log_grid = LOGICAL(s_log_grid)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
	bool use_uniform_weights = LOGICAL(s_use_uniform_weights)[0];

    // New parameters for buffer zone method
    size_t buffer_hops = (size_t)INTEGER(s_buffer_hops)[0];
    bool auto_buffer_hops = LOGICAL(s_auto_buffer_hops)[0];

    size_t n_folds = (size_t)INTEGER(s_n_folds)[0];
    bool with_bw_predictions = LOGICAL(s_with_bw_predictions)[0];
    double precision = REAL(s_precision)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    // Create graph and run algorithm
    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);
    graph_deg0_lowess_buffer_cv_t result = graph.graph_deg0_lowess_buffer_cv(
        y,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        log_grid,
        kernel_type,
        dist_normalization_factor,
		use_uniform_weights,
        buffer_hops,
        auto_buffer_hops,
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
        "opt_bw",
        "opt_bw_idx",
        "buffer_hops_used",  // New return value
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Track protection count
    int protect_count = 0;
    SEXP r_result = PROTECT(allocVector(VECSXP, n_elements));
    protect_count++;

    SEXP r_result_names = PROTECT(allocVector(STRSXP, n_elements));
    protect_count++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_result_names, i, mkChar(names[i]));
    }
    setAttrib(r_result, R_NamesSymbol, r_result_names);

    // Helper function to convert vector to SEXP
    auto vec_to_sexp = [&protect_count](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(REALSXP, vec.size()));
        protect_count++;
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set predictions
    SET_VECTOR_ELT(r_result, 0, vec_to_sexp(result.predictions));

    // Set bw_predictions (if requested)
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        SEXP r_bw_predictions = PROTECT(allocVector(VECSXP, result.bw_predictions.size()));
        protect_count++;

        for (size_t i = 0; i < result.bw_predictions.size(); i++) {
            SET_VECTOR_ELT(r_bw_predictions, i, vec_to_sexp(result.bw_predictions[i]));
        }

        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
    } else {
        SET_VECTOR_ELT(r_result, 1, allocVector(VECSXP, 0)); // Empty list
    }

    // Set bw_errors
    SET_VECTOR_ELT(r_result, 2, vec_to_sexp(result.bw_errors));

    // Set bws
    SET_VECTOR_ELT(r_result, 3, vec_to_sexp(result.bws));

    // Set opt_bw
    SEXP r_opt_bw = PROTECT(allocVector(REALSXP, 1));
    protect_count++;
    REAL(r_opt_bw)[0] = result.opt_bw;
    SET_VECTOR_ELT(r_result, 4, r_opt_bw);

    // Set opt_bw_idx
    SEXP r_opt_bw_idx = PROTECT(allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(r_opt_bw_idx)[0] = result.opt_bw_idx + 1; // Convert to 1-based for R
    SET_VECTOR_ELT(r_result, 5, r_opt_bw_idx);

    // Set buffer_hops_used (new return value)
    SEXP r_buffer_hops_used = PROTECT(allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(r_buffer_hops_used)[0] = result.buffer_hops_used;
    SET_VECTOR_ELT(r_result, 6, r_buffer_hops_used);

    UNPROTECT(protect_count);
    return r_result;
}
