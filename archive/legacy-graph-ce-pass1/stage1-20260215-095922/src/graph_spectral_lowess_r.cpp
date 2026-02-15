#include "graph_spectral_lowess.hpp"   // For graph_spectral_lowess_t
#include "error_utils.h"               // For REPORT_ERROR
#include "SEXP_cpp_conversion_utils.hpp"
#include "uniform_grid_graph.hpp"

#include <vector>
#include <queue>
#include <chrono>
#include <numeric> // for std::iota()

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

extern "C" {
    SEXP S_graph_spectral_lowess(
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
        // other
        SEXP s_precision,
        SEXP s_n_cleveland_iterations,
        SEXP s_verbose
        );
}

/**
 * @brief R interface to the graph spectral LOWESS algorithm
 *
 * @details This function provides an interface between R and the C++ implementation of graph spectral LOWESS.
 * It performs local regression on graph data using spectral embeddings:
 * 1. Constructs a weighted graph from the provided adjacency and weight lists
 * 2. Computes a spectral embedding of vertices using Laplacian eigenvectors
 * 3. For each vertex, fits local linear models at multiple bandwidths
 * 4. Selects optimal bandwidth based on cross-validation Rf_error
 * 5. Returns predictions, errors, and local scales
 *
 * @param s_adj_list [SEXP] R list of integer vectors representing adjacency lists for graph vertices
 * @param s_weight_list [SEXP] R list of numeric vectors with corresponding edge weights
 * @param s_y [SEXP] R numeric vector of response values for graph vertices
 * @param s_n_evectors [SEXP] R integer specifying number of eigenvectors to use in spectral embedding
 * @param s_n_bws [SEXP] R integer specifying number of candidate bandwidths to evaluate
 * @param s_log_grid [SEXP] R logical indicating whether to use logarithmic spacing for bandwidths
 * @param s_min_bw_factor [SEXP] R numeric for minimum bandwidth as fraction of graph diameter
 * @param s_max_bw_factor [SEXP] R numeric for maximum bandwidth as fraction of graph diameter
 * @param s_dist_normalization_factor [SEXP] R numeric for distance normalization in kernel weights
 * @param s_kernel_type [SEXP] R integer specifying the kernel function for vertex weighting
 * @param s_precision [SEXP] R numeric specifying precision tolerance for optimization
 * @param s_verbose [SEXP] R logical indicating whether to display progress information
 *
 * @return [SEXP] R list with components:
 *   - predictions: Numeric vector of fitted values for each vertex
 *   - errors: Numeric vector of leave-one-out cross-validation errors
 *   - scale: Numeric vector of optimal bandwidths (local scales) for each vertex
 *   - graph_diameter: Numeric scalar with computed graph diameter
 *
 * @note The function uses spectral embedding to transform graph distances into
 * a Euclidean space suitable for local linear regression.
 *
 * @see graph_spectral_lowess (C++ implementation)
 * @see set_wgraph_t::find_minimum_radius_for_domain_min_size
 */
SEXP S_graph_spectral_lowess(
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
    // other
    SEXP s_precision,
    SEXP s_n_cleveland_iterations,
    SEXP s_verbose
    ) {
    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

    // Scalars / flags (defensive extraction)
    const size_t n_evectors              = (size_t) Rf_asInteger(s_n_evectors);
    const size_t n_bws                   = (size_t) Rf_asInteger(s_n_bws);
    const bool   log_grid                = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor           = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor           = Rf_asReal(s_max_bw_factor);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type             = (size_t) Rf_asInteger(s_kernel_type);
    const double precision               = Rf_asReal(s_precision);
    const size_t n_cleveland_iterations  = (size_t) Rf_asInteger(s_n_cleveland_iterations);
    const bool   verbose                 = (Rf_asLogical(s_verbose) == TRUE);

    // Build grid graph and compute diameter endpoints
    uniform_grid_graph_t grid_graph = uniform_grid_graph_t(adj_list, weight_list);

    // Find diameter endpoints
    auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);   // Start from vertex 0
    auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
    grid_graph.graph_diameter = diameter;

    // Run core computation
    graph_spectral_lowess_t res = grid_graph.graph_spectral_lowess(
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
        // other
        precision,
        n_cleveland_iterations,
        verbose
        );

    // Assemble result: list(predictions, errors, scale, graph_diameter)
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));

    // names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 1, Rf_mkChar("errors"));
        SET_STRING_ELT(names, 2, Rf_mkChar("scale"));
        SET_STRING_ELT(names, 3, Rf_mkChar("graph_diameter"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    // 0: predictions
    {
        const R_xlen_t n = (R_xlen_t) res.predictions.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.predictions.begin(), res.predictions.end(), REAL(v));
        SET_VECTOR_ELT(result, 0, v);
        UNPROTECT(1); // v -> [1]
    }

    // 1: errors
    {
        const R_xlen_t n = (R_xlen_t) res.errors.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.errors.begin(), res.errors.end(), REAL(v));
        SET_VECTOR_ELT(result, 1, v);
        UNPROTECT(1); // v -> [1]
    }

    // 2: scale
    {
        const R_xlen_t n = (R_xlen_t) res.scale.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.scale.begin(), res.scale.end(), REAL(v));
        SET_VECTOR_ELT(result, 2, v);
        UNPROTECT(1); // v -> [1]
    }

    // 3: graph_diameter (scalar)
    {
        SEXP v = PROTECT(Rf_allocVector(REALSXP, 1)); // [2]
        REAL(v)[0] = grid_graph.graph_diameter;
        SET_VECTOR_ELT(result, 3, v);
        UNPROTECT(1); // v -> [1]
    }

    UNPROTECT(1); // result
    return result;
}
