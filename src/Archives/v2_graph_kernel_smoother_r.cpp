#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length

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

#include "graph_kernel_smoother.hpp" // For graph_kernel_smoother_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

extern "C" {
    SEXP S_graph_kernel_smoother(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_vertex_hbhd_min_size,
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
 * @brief R interface for graph kernel smoother with buffer zone cross-validation
 *
 * @details This function serves as the interface between R and C++ implementation
 * of the buffer zone cross-validation method for graph-based kernel smoothing. It handles
 * the conversion of R objects to C++ types, calls the C++ implementation, and
 * returns the results in R-compatible format.
 *
 * The implementation computes kernel-weighted estimates of E(Y|G) using graph-based distances
 * and provides multiple evaluation metrics to assess performance at different bandwidths.
 *
 * @param s_adj_list R list of integer vectors representing adjacency lists
 * @param s_weight_list R list of numeric vectors with edge weights
 * @param s_y R numeric vector of response values
 * @param s_min_bw_factor R numeric scalar for minimum bandwidth factor
 * @param s_max_bw_factor R numeric scalar for maximum bandwidth factor
 * @param s_n_bws R integer scalar for number of bandwidths to test
 * @param s_log_grid R logical scalar for log-spaced bandwidth grid
 * @param s_vertex_hbhd_min_size R integer for minimum neighborhood size
 * @param s_kernel_type R integer scalar for kernel function type
 * @param s_dist_normalization_factor R numeric scalar for distance normalization
 * @param s_use_uniform_weights R logical for uniform vs kernel weights
 * @param s_buffer_hops R integer scalar for buffer zone hop distance
 * @param s_auto_buffer_hops R logical scalar to auto-determine buffer size
 * @param s_n_folds R integer scalar for number of cross-validation folds
 * @param s_with_bw_predictions R logical scalar to store bandwidth predictions
 * @param s_precision R numeric scalar for precision
 * @param s_verbose R logical scalar for verbose output
 *
 * @return SEXP R list containing cross-validation results and predictions
 */
SEXP S_graph_kernel_smoother(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_vertex_hbhd_min_size,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_use_uniform_weights,
    SEXP s_buffer_hops,
    SEXP s_auto_buffer_hops,
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
    size_t vertex_hbhd_min_size = (size_t)INTEGER(s_vertex_hbhd_min_size)[0];
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
    graph_kernel_smoother_t result = graph.graph_kernel_smoother(
        y,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        log_grid,
        vertex_hbhd_min_size,
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

    // Create the return list - one entry for each field in graph_kernel_smoother_t
    const char* names[] = {
        "predictions",                  // Optimal bandwidth predictions
        "bw_predictions",               // Predictions for each bandwidth
        "bw_mean_sq_errors",            // Mean squared errors for each bandwidth
        "bw_mean_abs_errors",           // Mean absolute errors for each bandwidth
        "vertex_min_bws",               // Minimum bandwidth for each vertex
        "opt_bw_idx",                   // Index of optimal bandwidth
        "buffer_hops_used",             // Buffer hops used
        "bw_lextr_count",               // Local extrema count for each bandwidth
        "bw_total_sq_curvature",        // Total squared curvature for each bandwidth
        "bw_total_sq_norma_curvature",  // Total squared normalized curvature
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

    auto size_t_vec_to_sexp = [&protect_count](const std::vector<size_t>& vec) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(INTSXP, vec.size()));
        protect_count++;
        int* ptr = INTEGER(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set predictions
    SET_VECTOR_ELT(r_result, 0, vec_to_sexp(result.predictions));

    // Set bw_predictions (if requested) as a matrix instead of list
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        // Determine matrix dimensions
        // Number of rows = length of each vector
        // Number of columns = number of vectors
        size_t ncols = result.bw_predictions.size();
        size_t nrows = 0;

        // Find the length of vectors to determine number of rows
        // Assuming all vectors have the same length
        if (ncols > 0) {
            nrows = result.bw_predictions[0].size();
        }

        // Create a numeric matrix
        SEXP r_bw_predictions = PROTECT(allocMatrix(REALSXP, nrows, ncols));
        protect_count++;

        // Fill the matrix with values
        double *matrix_ptr = REAL(r_bw_predictions);

        for (size_t j = 0; j < ncols; j++) {
            const auto& col = result.bw_predictions[j];
            // Check if the column has the expected number of rows
            if (col.size() != nrows) {
                error("Inconsistent vector lengths in bw_predictions");
            }

            // Copy data to the matrix (R matrices are column-major)
            // Each vector becomes one column
            for (size_t i = 0; i < nrows; i++) {
                matrix_ptr[i + j * nrows] = col[i];
            }
        }

        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
    } else {
        // Empty matrix 0x0
        SEXP r_bw_predictions = PROTECT(allocMatrix(REALSXP, 0, 0));
        protect_count++;
        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
    }

    // Set bw_mean_sq_errors
    SET_VECTOR_ELT(r_result, 2, vec_to_sexp(result.bw_mean_sq_errors));

    // Set bw_mean_abs_errors
    SET_VECTOR_ELT(r_result, 3, vec_to_sexp(result.bw_mean_abs_errors));

    // Set vertex_min_bws
    SET_VECTOR_ELT(r_result, 4, vec_to_sexp(result.vertex_min_bws));

    // Set opt_bw_idx
    SEXP r_opt_bw_idx = PROTECT(allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(r_opt_bw_idx)[0] = result.opt_bw_idx + 1; // Convert to 1-based for R
    SET_VECTOR_ELT(r_result, 5, r_opt_bw_idx);

    // Set buffer_hops_used
    SEXP r_buffer_hops_used = PROTECT(allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(r_buffer_hops_used)[0] = result.buffer_hops_used;
    SET_VECTOR_ELT(r_result, 6, r_buffer_hops_used);

    // Set bw_lextr_count
    SET_VECTOR_ELT(r_result, 7, size_t_vec_to_sexp(result.bw_lextr_count));

    // Set bw_total_sq_curvature
    SET_VECTOR_ELT(r_result, 8, vec_to_sexp(result.bw_total_sq_curvature));

    // Set bw_total_sq_norma_curvature
    SET_VECTOR_ELT(r_result, 9, vec_to_sexp(result.bw_total_sq_norma_curvature));

    UNPROTECT(protect_count);
    return r_result;
}
