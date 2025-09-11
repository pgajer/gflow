#include "graph_kernel_smoother.hpp" // For graph_kernel_smoother_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <mutex>                    // For std::mutex
#include <execution>                // For std::execution::par_unseq
#include <atomic>                   // For std::atomic

#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions

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
 * @brief  R interface for graph kernel smoother with buffer‐zone cross‐validation
 *
 * @details Converts R inputs into C++ types, invokes set_wgraph_t::graph_kernel_smoother,
 * and returns only the fields in graph_kernel_smoother_t:
 *   - predictions:          Optimal‐bandwidth predictions (length |V|)
 *   - bw_predictions:       Per‐bandwidth predictions (|V|×n_bws matrix)
 *   - bw_mean_abs_errors:   Mean absolute Rf_error for each bandwidth (length n_bws)
 *   - vertex_min_bws:       Per‐vertex minimum bandwidths (length |V|)
 *   - opt_bw_idx:           1‐based index of the chosen bandwidth in R
 *   - buffer_hops_used:     Number of buffer hops actually applied
 *
 * @param s_adj_list                R list of integer vectors for adjacency
 * @param s_weight_list             R list of numeric vectors for edge weights
 * @param s_y                       Numeric vector of responses
 * @param s_min_bw_factor           Numeric scalar: min bandwidth factor
 * @param s_max_bw_factor           Numeric scalar: max bandwidth factor
 * @param s_n_bws                   Integer scalar: number of bandwidths
 * @param s_log_grid                Logical scalar: log‐spaced grid?
 * @param s_vertex_hbhd_min_size    Integer: min neighborhood size per vertex
 * @param s_kernel_type             Integer: kernel type code
 * @param s_dist_normalization_factor Numeric: distance normalization factor
 * @param s_use_uniform_weights     Logical: uniform vs. kernel weights
 * @param s_buffer_hops             Integer: buffer‐zone hop distance
 * @param s_auto_buffer_hops        Logical: auto‐determine buffer hops?
 * @param s_n_folds                 Integer: number of CV folds
 * @param s_with_bw_predictions     Logical: store per‐bandwidth predictions?
 * @param s_precision               Numeric: precision for bandwidth grid
 * @param s_verbose                 Logical: verbose output?
 *
 * @return SEXP R list with components:
 *         "predictions", "bw_predictions", "bw_mean_abs_errors",
 *         "vertex_min_bws", "opt_bw_idx", "buffer_hops_used".
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
    size_t vertex_hbhd_min_size = (size_t)INTEGER(s_vertex_hbhd_min_size)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool use_uniform_weights = (LOGICAL(s_use_uniform_weights)[0] == 1);

    // New parameters for buffer zone method
    size_t buffer_hops = (size_t)INTEGER(s_buffer_hops)[0];
    bool auto_buffer_hops = (LOGICAL(s_auto_buffer_hops)[0] == 1);

    size_t n_folds = (size_t)INTEGER(s_n_folds)[0];
    bool with_bw_predictions = (LOGICAL(s_with_bw_predictions)[0] == 1);
    double precision = REAL(s_precision)[0];
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    // Run the C++ smoother
    set_wgraph_t graph(adj_list, weight_list);
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

    // Prepare return list names
    const char* names[] = {
        "predictions",
        "bw_predictions",
        "bw_mean_abs_errors",
        "vertex_min_bws",
        "opt_bw_idx",
        "buffer_hops_used",
        NULL
    };
    int n_elements = 6;

    // Allocate result list
    int protect_count = 0;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements)); protect_count++;
    SEXP r_names  = PROTECT(Rf_allocVector(STRSXP, n_elements)); protect_count++;
    for (int i = 0; i < n_elements; ++i) {
        SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    // Helper lambdas
    auto vec_to_sexp = [&protect_count](const std::vector<double>& v) {
        SEXP x = PROTECT(Rf_allocVector(REALSXP, v.size())); protect_count++;
        std::copy(v.begin(), v.end(), REAL(x));
        return x;
    };
    auto mat_to_sexp = [&protect_count](const std::vector<std::vector<double>>& M) {
        size_t ncol = M.size();
        size_t nrow = ncol ? M[0].size() : 0;
        SEXP x = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); protect_count++;
        double* ptr = REAL(x);
        for (size_t j = 0; j < ncol; ++j) {
            if (M[j].size() != nrow)
                Rf_error("Inconsistent bw_predictions dimensions");
            for (size_t i = 0; i < nrow; ++i)
                ptr[i + j * nrow] = M[j][i];
        }
        return x;
    };
    auto scalar_int = [&protect_count](size_t v, bool one_based=false) {
        SEXP x = PROTECT(Rf_allocVector(INTSXP, 1)); protect_count++;
        INTEGER(x)[0] = one_based ? (int)v + 1 : (int)v;
        return x;
    };

    // Fill in fields
    SET_VECTOR_ELT(r_result, 0, vec_to_sexp(result.predictions));
    SET_VECTOR_ELT(r_result, 1, mat_to_sexp(result.bw_predictions));
    SET_VECTOR_ELT(r_result, 2, vec_to_sexp(result.bw_mean_abs_errors));
    SET_VECTOR_ELT(r_result, 3, vec_to_sexp(result.vertex_min_bws));
    SET_VECTOR_ELT(r_result, 4, scalar_int(result.opt_bw_idx, true));
    SET_VECTOR_ELT(r_result, 5, scalar_int(result.buffer_hops_used));

    UNPROTECT(protect_count);
    return r_result;
}
