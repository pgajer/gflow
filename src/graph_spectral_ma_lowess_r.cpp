#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <vector>
#include <chrono>
#include <numeric> // for std::iota()

#include "graph_spectral_lowess.hpp"  // For graph_spectral_lowess_result_t
#include "bandwidth_utils.hpp"        // For get_candidate_bws
#include "error_utils.h"              // For REPORT_ERROR
#include "set_wgraph.hpp"             // For set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "progress_utils.hpp"         // For elapsed_time

extern "C" {
    SEXP S_graph_spectral_ma_lowess(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_n_evectors,
        // bw parameters
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        // kernel parameters
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        // model parameters
        SEXP s_blending_coef,
        // other
        SEXP s_precision,
        SEXP s_verbose
    );
}

/**
 * @brief R interface for graph_spectral_ma_lowess
 *
 * Handles conversion between R and C++ data structures for the model-averaged
 * spectral LOWESS algorithm.
 */
SEXP S_graph_spectral_ma_lowess(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_n_evectors,
    // bw parameters
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    // kernel parameters
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    // model parameters
    SEXP s_blending_coef,
    // other
    SEXP s_precision,
    SEXP s_verbose
) {
    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    size_t n_evectors = (size_t)INTEGER(s_n_evectors)[0];

    // bw parameters
    size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
    bool log_grid = LOGICAL(s_log_grid)[0];
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];

    // kernel parameters
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];

    // model parameters
    double blending_coef = REAL(s_blending_coef)[0];

    // other
    double precision = REAL(s_precision)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    // Create the graph
    set_wgraph_t graph(adj_list, weight_list);

    // Find diameter if not already computed
    if (graph.graph_diameter <= 0) {
        auto [end1, diam] = graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = graph.get_vertex_eccentricity(end1);
        graph.graph_diameter = diameter;
    }

    // Call the C++ function
    graph_spectral_lowess_t res = graph.graph_spectral_ma_lowess(
        y,
        n_evectors,
        // bw parameters
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        // kernel parameters
        dist_normalization_factor,
        kernel_type,
        // model parameters
        blending_coef,
        // other
        precision,
        verbose
    );

    // Create the return list using R's C API
    const char* names[] = {
        "predictions",
        "errors",
        "scale",
        "graph_diameter",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Create list and protect it
    int protect_count = 0;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements)); protect_count++;

    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements)); protect_count++;
    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // Helper function to convert vector to SEXP
    auto create_numeric_vector = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    SET_VECTOR_ELT(result, 0, create_numeric_vector(res.predictions)); protect_count++;
    SET_VECTOR_ELT(result, 1, create_numeric_vector(res.errors)); protect_count++;
    SET_VECTOR_ELT(result, 2, create_numeric_vector(res.scale)); protect_count++;

    SEXP diam = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;
    REAL(diam)[0] = graph.graph_diameter;
    SET_VECTOR_ELT(result, 3, diam);

    UNPROTECT(protect_count);
    return result;
}
