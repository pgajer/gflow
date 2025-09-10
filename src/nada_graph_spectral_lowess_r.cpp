#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

#include <vector>
#include <queue>
#include <chrono>
#include <numeric> // for std::iota()

#include "nada_graph_spectral_lowess.hpp"   // For graph_spectral_lowess_t
#include "error_utils.h"               // For REPORT_ERROR
#include "SEXP_cpp_conversion_utils.hpp"
#include "uniform_grid_graph.hpp"

extern "C" {
	SEXP S_nada_graph_spectral_lowess(
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
SEXP S_nada_graph_spectral_lowess(
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

    // other
    double precision = REAL(s_precision)[0];
	size_t n_cleveland_iterations = (size_t)INTEGER(s_n_cleveland_iterations)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);

    // Find diameter endpoints
    auto [end1, diam] = graph.get_vertex_eccentricity(0);  // Start from vertex 0
    auto [end2, diameter] = graph.get_vertex_eccentricity(end1);
    graph.graph_diameter = diameter;


    // Call the C++ function
    nada_graph_spectral_lowess_t res = graph.nada_graph_spectral_lowess(
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

    // Create the return list using R's C API
    const char* names[] = {
        "predictions",
        "bw_predictions",
        "errors",
        "global_bws",
        "opt_bw_idx",
        "opt_bw",
        "graph_diameter",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Create list and protect it
    int protect_count = 0;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    protect_count++;

    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    protect_count++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // Helper function to convert vector to SEXP
    auto vec_to_sexp = [&protect_count](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        protect_count++;
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set values
    SET_VECTOR_ELT(result, 0, vec_to_sexp(res.predictions));


    // Convert matrix of bw predictions
    SEXP bw_pred_r = PROTECT(Rf_allocMatrix(REALSXP,
                                            res.bw_predictions[0].size(),
                                            res.bw_predictions.size())); protect_count++;
    for(size_t i = 0; i < res.bw_predictions.size(); i++) {
        std::copy(res.bw_predictions[i].begin(),
                  res.bw_predictions[i].end(),
                  REAL(bw_pred_r) + i * res.bw_predictions[0].size());
    }
    SET_VECTOR_ELT(result, 1, bw_pred_r);

    SET_VECTOR_ELT(result, 2, vec_to_sexp(res.errors));
    SET_VECTOR_ELT(result, 3, vec_to_sexp(res.global_bws));

    SEXP s_opt_bw_idx = PROTECT(Rf_allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(s_opt_bw_idx)[0] = res.opt_bw_idx;
    SET_VECTOR_ELT(result, 4, s_opt_bw_idx);

    // result.opt_bw = global_bws[opt_bw_idx];
    SEXP s_opt_bw = PROTECT(Rf_allocVector(REALSXP, 1));
    protect_count++;
    REAL(s_opt_bw)[0] = res.opt_bw;
    SET_VECTOR_ELT(result, 5, s_opt_bw);

    SEXP s_graph_diameter = PROTECT(Rf_allocVector(REALSXP, 1));
    protect_count++;
    REAL(s_graph_diameter)[0] = graph.graph_diameter;
    SET_VECTOR_ELT(result, 6, s_graph_diameter);

    UNPROTECT(protect_count);
    return result;
}
