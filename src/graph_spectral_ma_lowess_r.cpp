#include "graph_spectral_lowess.hpp"  // For graph_spectral_lowess_result_t
#include "bandwidth_utils.hpp"        // For get_candidate_bws
#include "error_utils.h"              // For REPORT_ERROR
#include "set_wgraph.hpp"             // For set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "progress_utils.hpp"         // For elapsed_time

#include <vector>
#include <chrono>
#include <numeric> // for std::iota()

#include <R.h>
#include <Rinternals.h>

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
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

    // -------- Scalars / flags (defensive extraction) --------
    const size_t n_evectors = static_cast<size_t>(Rf_asInteger(s_n_evectors));
    // bw parameters
    const size_t n_bws      = static_cast<size_t>(Rf_asInteger(s_n_bws));
    const bool   log_grid   = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);
    // kernel parameters
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type = static_cast<size_t>(Rf_asInteger(s_kernel_type));
    // model parameters
    const double blending_coef = Rf_asReal(s_blending_coef);
    // other
    const double precision = Rf_asReal(s_precision);
    const bool   verbose   = (Rf_asLogical(s_verbose) == TRUE);

    // -------- Create the graph --------
    set_wgraph_t graph(adj_list, weight_list);

    // -------- Ensure diameter is computed --------
    if (graph.graph_diameter <= 0) {
        auto [end1, diam] = graph.get_vertex_eccentricity(0);  // Start from vertex 0
        (void)diam;
        auto [end2, diameter] = graph.get_vertex_eccentricity(end1);
        (void)end2;
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

    // -------- Build result (container-first; fixed UNPROTECT counts) --------
    const int n_elements = 4;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));

    // names
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    SET_STRING_ELT(result_names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(result_names, 1, Rf_mkChar("errors"));
    SET_STRING_ELT(result_names, 2, Rf_mkChar("scale"));
    SET_STRING_ELT(result_names, 3, Rf_mkChar("graph_diameter"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(1); // result_names

    // Helper function to convert vector to SEXP
    auto create_numeric_vector = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        std::copy(vec.begin(), vec.end(), REAL(r_vec));
        return r_vec;
    };

    // 0: predictions
    {
        SEXP s = create_numeric_vector(res.predictions);
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1); // s - protected in create_numeric_vector()
    }

    // 1: errors
    {
        SEXP s = create_numeric_vector(res.errors);
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1); // s - protected in create_numeric_vector()
    }

    // 2: scale
    {
        SEXP s = create_numeric_vector(res.scale);
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1); // s - protected in create_numeric_vector()
    }

    // 3: graph_diameter
    {
        SEXP diam = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(diam)[0] = graph.graph_diameter;
        SET_VECTOR_ELT(result, 3, diam);
        UNPROTECT(1);
    }

    UNPROTECT(1); // result
    return result;
}
