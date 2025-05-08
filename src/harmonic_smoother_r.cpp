#include <R.h>
#include <Rinternals.h>
 // Undefine conflicting macros after including R headers
#undef length

#include "harmonic_smoother.hpp"
#include "set_wgraph.hpp"  // for set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
    SEXP S_perform_harmonic_smoothing(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_harmonic_predictions,
        SEXP s_region_vertices,
        SEXP s_max_iterations,
        SEXP s_tolerance
        );

	SEXP S_harmonic_smoother(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_harmonic_predictions,
		SEXP s_region_vertices,
		SEXP s_max_iterations,
		SEXP s_tolerance,
		SEXP s_record_frequency,
		SEXP s_stability_window,
		SEXP s_stability_threshold
		);
}

/**
 * @brief R interface to perform harmonic smoothing on a weighted graph
 *
 * @details This function provides an R interface to the C++ implementation of harmonic
 *          smoothing on a weighted graph. It takes a graph structure (adjacency list and
 *          weights), function values to be smoothed, a set of region vertices defining
 *          the smoothing domain, and convergence parameters. The function preserves
 *          values at the boundary of the region while smoothing interior values.
 *
 * The function performs the following steps:
 * 1. Converts R data structures to C++ structures
 * 2. Constructs a weighted graph from the adjacency and weight lists
 * 3. Creates an unordered set of region vertices for efficient lookups
 * 4. Performs harmonic smoothing within the specified region
 * 5. Returns the updated function values as an R numeric vector
 *
 * @param s_adj_list R list of integer vectors representing the adjacency list
 * @param s_weight_list R list of numeric vectors representing edge weights
 * @param s_harmonic_predictions R numeric vector of function values to be smoothed
 * @param s_region_vertices R integer vector of vertex indices defining the region (1-based)
 * @param s_max_iterations R integer scalar for maximum number of iterations
 * @param s_tolerance R numeric scalar for convergence tolerance
 *
 * @return R numeric vector containing the smoothed function values
 *
 * @note This function expects 1-based vertex indices from R and converts them to 0-based
 *       indices for C++ processing. The returned values maintain the original vector size.
 * @note The region vertices must be valid indices within the range of the graph size.
 *
 * @see perform_harmonic_smoothing
 * @see set_wgraph_t
 */
SEXP S_perform_harmonic_smoothing(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_harmonic_predictions,
    SEXP s_region_vertices,
    SEXP s_max_iterations,
    SEXP s_tolerance
    ) {

    // Convert input from R format
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert function values to smooth
    std::vector<double> harmonic_predictions(REAL(s_harmonic_predictions),
											 REAL(s_harmonic_predictions) + LENGTH(s_harmonic_predictions));

    // Convert region vertices (using actual length of s_region_vertices)
    size_t n_region_vertices = LENGTH(s_region_vertices);

    // Create unordered set for efficient lookup
    std::unordered_set<size_t> region_vertices_set;
    for (size_t i = 0; i < n_region_vertices; ++i) {
        // Convert 1-based R indices to 0-based C++ indices
        region_vertices_set.insert(INTEGER(s_region_vertices)[i] - 1);
    }

    // Get smoothing parameters
    int max_iterations = INTEGER(s_max_iterations)[0];
    double tolerance = REAL(s_tolerance)[0];

    // Create and initialize the graph
    set_wgraph_t graph(adj_list, weight_list);
    graph.compute_graph_diameter();

    // Perform the harmonic smoothing
    graph.perform_harmonic_smoothing(
        harmonic_predictions,
        region_vertices_set,
        max_iterations,
        tolerance);

    // Convert back to R numeric vector
    SEXP r_harmonic_predictions = PROTECT(allocVector(REALSXP, harmonic_predictions.size()));
    std::copy(harmonic_predictions.begin(), harmonic_predictions.end(), REAL(r_harmonic_predictions));
    UNPROTECT(1);

    return r_harmonic_predictions;
}


/**
 * @brief R interface to advanced harmonic smoothing with prediction landscape tracking
 *
 * @details This function provides an R interface to the C++ implementation of harmonic
 *          smoothing with topology tracking. It allows monitoring how local extrema
 *          and their basins evolve during the smoothing process to identify the optimal
 *          amount of smoothing that preserves significant features.
 *
 * The function performs the following steps:
 * 1. Converts R data structures to C++ structures
 * 2. Constructs a weighted graph from the adjacency and weight lists
 * 3. Creates an unordered set of region vertices for efficient lookups
 * 4. Performs harmonic smoothing with topology tracking
 * 5. Converts the results back to R-compatible structures
 * 6. Returns a list containing smoothed values and topology evolution information
 *
 * @param s_adj_list R list of integer vectors representing the adjacency list
 * @param s_weight_list R list of numeric vectors representing edge weights
 * @param s_harmonic_predictions R numeric vector of function values to be smoothed
 * @param s_region_vertices R integer vector of vertex indices defining the region (1-based)
 * @param s_max_iterations R integer scalar for maximum number of iterations
 * @param s_tolerance R numeric scalar for convergence tolerance
 * @param s_record_frequency R integer scalar for frequency of recording states
 * @param s_stability_window R integer scalar for window size to check stability
 * @param s_stability_threshold R numeric scalar for topology stability threshold
 *
 * @return R list containing:
 *         - harmonic_predictions: Smoothed function values
 *         - i_harmonic_predictions: Matrix of function values at recorded iterations
 *         - i_basins: List of matrices representing extrema at each iteration
 *         - stable_iteration: Iteration at which topology stabilized
 *         - topology_differences: Differences between consecutive iterations
 *
 * @note This function expects 1-based vertex indices from R and converts them to 0-based
 *       indices for C++ processing. All returned indices maintain the 1-based convention.
 *
 * @see harmonic_smoother
 * @see set_wgraph_t
 * @see compute_basins
 */
SEXP S_harmonic_smoother(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_harmonic_predictions,
    SEXP s_region_vertices,
    SEXP s_max_iterations,
    SEXP s_tolerance,
    SEXP s_record_frequency,
    SEXP s_stability_window,
    SEXP s_stability_threshold
	) {

    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    std::vector<double> harmonic_predictions(REAL(s_harmonic_predictions),
                                             REAL(s_harmonic_predictions) + LENGTH(s_harmonic_predictions));

    size_t n_region_vertices = LENGTH(s_region_vertices);
    std::unordered_set<size_t> region_vertices_set;
    for (size_t i = 0; i < n_region_vertices; ++i) {
        region_vertices_set.insert(INTEGER(s_region_vertices)[i] - 1);
    }

    int max_iterations         = INTEGER(s_max_iterations)[0];
    double tolerance           = REAL(s_tolerance)[0];
    int record_frequency       = INTEGER(s_record_frequency)[0];
    size_t stability_window    = static_cast<size_t>(INTEGER(s_stability_window)[0]);
    double stability_threshold = REAL(s_stability_threshold)[0];

    set_wgraph_t graph(adj_list, weight_list);
    graph.compute_graph_diameter();

    harmonic_smoother_t result = graph.harmonic_smoother(
        harmonic_predictions,
        region_vertices_set,
        max_iterations,
        tolerance,
        record_frequency,
        stability_window,
        stability_threshold
    );

    const char* names[] = {
        "harmonic_predictions",
        "i_harmonic_predictions",
        "i_basins",
        "stable_iteration",
        "basin_cx_differences",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements]) ++n_elements;


    SEXP r_result = PROTECT(allocVector(VECSXP, n_elements));
    SEXP r_result_names = PROTECT(allocVector(STRSXP, n_elements));

    for (int i = 0; i < n_elements; ++i) {
        SET_STRING_ELT(r_result_names, i, mkChar(names[i]));
    }
    setAttrib(r_result, R_NamesSymbol, r_result_names);
	UNPROTECT(1); // r_result_names

	// r_harmonic_predictions
    SEXP r_harmonic_predictions = PROTECT(allocVector(REALSXP, harmonic_predictions.size()));
    std::copy(harmonic_predictions.begin(), harmonic_predictions.end(), REAL(r_harmonic_predictions));
    SET_VECTOR_ELT(r_result, 0, r_harmonic_predictions);
	UNPROTECT(1);

    int nvert = result.i_harmonic_predictions[0].size();
    int niters = result.i_harmonic_predictions.size();
    SEXP r_i_harmonic_predictions = PROTECT(allocMatrix(REALSXP, nvert, niters));
    for (int j = 0; j < niters; ++j) {
        for (int i = 0; i < nvert; ++i) {
            REAL(r_i_harmonic_predictions)[i + nvert * j] = result.i_harmonic_predictions[j][i];
        }
    }
    SET_VECTOR_ELT(r_result, 1, r_i_harmonic_predictions);
	UNPROTECT(1); // r_harmonic_predictions

	// r_i_basins
    SEXP r_i_basins = PROTECT(allocVector(VECSXP, result.i_basins.size()));
    for (size_t iter = 0; iter < result.i_basins.size(); ++iter) {
        const auto& basin_map = result.i_basins[iter];
        size_t n_extrema = basin_map.size();
        SEXP basin_matrix = PROTECT(allocMatrix(REALSXP, n_extrema, 2));

        size_t row = 0;
        for (const auto& [vertex, basin] : basin_map) {
            REAL(basin_matrix)[row] = vertex + 1;
            REAL(basin_matrix)[row + n_extrema] = basin.is_maximum;
            row++;
        }

        SEXP col_names = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(col_names, 0, mkChar("evertex"));
        SET_STRING_ELT(col_names, 1, mkChar("is_max"));
        setAttrib(basin_matrix, R_DimNamesSymbol, list2(R_NilValue, col_names));

        SET_VECTOR_ELT(r_i_basins, iter, basin_matrix);
        UNPROTECT(2); // basin_matrix, col_names
    }
    SET_VECTOR_ELT(r_result, 2, r_i_basins);
	UNPROTECT(1); // r_i_basins


    SEXP r_stable_iteration = PROTECT(ScalarInteger(result.stable_iteration));
    SET_VECTOR_ELT(r_result, 3, r_stable_iteration);
	UNPROTECT(1);

    SEXP r_basin_cx_differences = PROTECT(allocVector(REALSXP, result.basin_cx_differences.size()));
    std::copy(result.basin_cx_differences.begin(), result.basin_cx_differences.end(), REAL(r_basin_cx_differences));
    SET_VECTOR_ELT(r_result, 4, r_basin_cx_differences);
	UNPROTECT(1);

    UNPROTECT(1); // r_result_names, r_result, each SET_VECTOR_ELT target (5)
	return r_result;
}
