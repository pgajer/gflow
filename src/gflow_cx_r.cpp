#include <R.h>
#include <Rinternals.h>
   // Undefine conflicting macros after including R headers
#undef length
#undef eval

#include "set_wgraph.hpp"  // for set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"
#include "gflow_cx.hpp"
#include "cpp_utils.hpp"

extern "C" {
	SEXP S_compute_extrema_hop_nbhds(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y
		);

	SEXP S_create_gflow_cx(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_hop_idx_thld,
		SEXP s_smoother_type,
		SEXP s_max_outer_iterations,
		SEXP s_max_inner_iterations,
		SEXP s_smoothing_tolerance,
		SEXP s_sigma,
		SEXP s_process_in_order,
		SEXP s_verbose,
		SEXP s_detailed_recording
		);

	SEXP S_apply_harmonic_extension(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_region_vertices,
		SEXP s_boundary_values,
		SEXP s_smoother_type,
		SEXP s_max_iterations,
		SEXP s_tolerance,
		SEXP s_sigma,
		SEXP s_record_iterations,
		SEXP s_verbose
		);
}

/**
 * @brief R interface for computing extrema hop neighborhoods in a graph-based function
 *
 * @details This function analyzes a function defined on a graph to identify local extrema
 *          (both minima and maxima) and computes their hop neighborhoods. For each extremum,
 *          it determines the maximum hop distance at which the extremum property is maintained,
 *          along with information about the boundary where this property is first violated.
 *
 *          The function:
 *          1. Constructs a weighted graph from the provided adjacency and weight lists
 *          2. Computes the graph diameter for reference
 *          3. Calls the C++ implementation of `compute_extrema_hop_nbhds` to find all extrema
 *             and their hop neighborhoods
 *          4. Converts these results to R data structures
 *
 * @param s_adj_list SEXP representing the graph's adjacency list (from R)
 * @param s_weight_list SEXP representing the edge weights for the graph
 * @param s_y SEXP containing function values at each vertex as a numeric vector
 *
 * @return SEXP containing a named list with two components:
 *         - lmin_hop_nbhds: List of hop neighborhoods for local minima
 *         - lmax_hop_nbhds: List of hop neighborhoods for local maxima
 *         Each neighborhood contains:
 *         - vertex: The extremum vertex index (1-based for R)
 *         - hop_idx: The maximum hop radius where extremum condition holds
 *         - nbhd_df: Matrix of vertices and their hop distances within the neighborhood
 *         - nbhd_bd_df: Matrix of boundary vertices and their function values
 *
 * @note All vertex indices are converted to 1-based indexing for R compatibility
 */
SEXP S_compute_extrema_hop_nbhds(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y
	) {

	// ---- Input validation and conversion ----
	std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	// Convert numeric vector directly
	double* y_ptr = REAL(s_y);
	std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

	// ---- Call C++ implementation ----
	set_wgraph_t graph(adj_list, weight_list);

	// Ensure graph diameter is available
	graph.compute_graph_diameter();

	// Compute all extrema hop neighborhoods at once
	auto [lmin_hop_nbhd_map, lmax_hop_nbhd_map] = graph.compute_extrema_hop_nbhds(y);

	// ---- Build R objects ----
	SEXP r_lmin_hop_nbhds = PROTECT(allocVector(VECSXP, lmin_hop_nbhd_map.size()));
	SEXP r_lmax_hop_nbhds = PROTECT(allocVector(VECSXP, lmax_hop_nbhd_map.size()));

	// Helper lambda to populate one neighborhood
	auto fill_nbhd = [&](const hop_nbhd_t& nbhd, SEXP r_nbhd) {
		// Names: vertex, hop_idx, nbhd_df, nbhd_bd_df
		SEXP names = PROTECT(allocVector(STRSXP, 5));
		SET_STRING_ELT(names, 0, mkChar("vertex"));
		SET_STRING_ELT(names, 1, mkChar("value"));
		SET_STRING_ELT(names, 2, mkChar("hop_idx"));
		SET_STRING_ELT(names, 3, mkChar("nbhd_df"));
		SET_STRING_ELT(names, 4, mkChar("nbhd_bd_df"));
		setAttrib(r_nbhd, R_NamesSymbol, names);

		bool is_global_extremum = (nbhd.hop_idx == std::numeric_limits<size_t>::max());

		// vertex (1-based)
		SEXP r_vertex = PROTECT(ScalarInteger(
									static_cast<int>(nbhd.vertex) + 1
									));

		// value
		SEXP r_value  = PROTECT(ScalarReal(y[nbhd.vertex]));

		// hop index
		SEXP r_hop_idx;
		if (is_global_extremum) {
			r_hop_idx = PROTECT(ScalarReal(R_PosInf));
		} else {
			r_hop_idx = PROTECT(ScalarReal(static_cast<double>(nbhd.hop_idx)));
		}

		// nbhd matrix
		SEXP r_nbhd_df;
		if (is_global_extremum) {
			r_nbhd_df = PROTECT(allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = nbhd.hop_dist_map.size();
			r_nbhd_df = PROTECT(allocMatrix(REALSXP, m, 2));
			double *pr = REAL(r_nbhd_df);
			size_t i = 0;
			for (const auto& [v,d] : nbhd.hop_dist_map) {
				pr[i]     = v + 1;
				pr[i + m] = d;
				i++;
			}
		}

		// nbhd_bd matrix
		SEXP r_nbhd_bd_df;
		if (is_global_extremum) {
			r_nbhd_bd_df = PROTECT(allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = nbhd.y_nbhd_bd_map.size();
			r_nbhd_bd_df = PROTECT(allocMatrix(REALSXP, m, 2));
			double *pr_bd = REAL(r_nbhd_bd_df);
			size_t i = 0;
			for (const auto& [v, y_val] : nbhd.y_nbhd_bd_map) {
				pr_bd[i]     = v + 1;
				pr_bd[i + m] = y_val;
				i++;
			}
		}

		// Set list entries
		SET_VECTOR_ELT(r_nbhd, 0, r_vertex);
		SET_VECTOR_ELT(r_nbhd, 1, r_value);
		SET_VECTOR_ELT(r_nbhd, 2, r_hop_idx);
		SET_VECTOR_ELT(r_nbhd, 3, r_nbhd_df);
		SET_VECTOR_ELT(r_nbhd, 4, r_nbhd_bd_df);

		UNPROTECT(6); // names, r_vertex, r_value, r_hop_idx, r_nbhd_df, r_nbhd_bd_df
	};

	// Fill minima nbhd's
	size_t i = 0;
	for (const auto& [_, nbhd] : lmin_hop_nbhd_map) {
		SEXP r_nbhd = PROTECT(allocVector(VECSXP, 5));
		fill_nbhd(nbhd, r_nbhd);
		SET_VECTOR_ELT(r_lmin_hop_nbhds, i++, r_nbhd);
		UNPROTECT(1); // r_nbhd
	}

	// Fill maxima nbhd's
	i = 0;
	for (const auto& [_, nbhd] : lmax_hop_nbhd_map) {
		SEXP r_nbhd = PROTECT(allocVector(VECSXP, 5));
		fill_nbhd(nbhd, r_nbhd);
		SET_VECTOR_ELT(r_lmax_hop_nbhds, i++, r_nbhd);
		UNPROTECT(1); // r_nbhd
	}

	// Name the two top‐level components
	SEXP r_result = PROTECT(allocVector(VECSXP, 2));
	SEXP r_result_names = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(r_result_names, 0, mkChar("lmin_hop_nbhds"));
	SET_STRING_ELT(r_result_names, 1, mkChar("lmax_hop_nbhds"));
	setAttrib(r_result, R_NamesSymbol, r_result_names);
	SET_VECTOR_ELT(r_result, 0, r_lmin_hop_nbhds);
	SET_VECTOR_ELT(r_result, 1, r_lmax_hop_nbhds);

	UNPROTECT(4); // r_lmin_hop_nbhds, r_lmax_hop_nbhds, r_result, r_result_names

	return r_result;
}


/**
 * @brief R interface for creating a gradient flow complex with harmonic extension
 *
 * @details This function identifies local extrema in a function defined on a graph and
 *          applies harmonic extension to smooth out spurious extrema. It first computes
 *          the hop neighborhoods for all local extrema, identifying which are spurious
 *          based on the provided hop index threshold. Then it applies the specified
 *          harmonic extension method to smooth these spurious extrema.
 *
 *          The function:
 *          1. Constructs a weighted graph from the provided adjacency and weight lists
 *          2. Identifies local extrema (minima and maxima) and their hop neighborhoods
 *          3. Applies the specified smoothing method to eliminate spurious extrema
 *          4. Returns the smoothed function values and extrema information
 *
 * @param s_adj_list SEXP representing the graph's adjacency list (from R)
 * @param s_weight_list SEXP representing the edge weights for the graph
 * @param s_y SEXP containing function values at each vertex as a numeric vector
 * @param s_hop_idx_thld SEXP integer threshold for hop index; extrema with hop_idx <= this value
 *                      are considered spurious and will be smoothed
 * @param s_smoother_type SEXP integer selecting the smoothing method to use:
 *                      0: Weighted Mean
 *                      1: Harmonic Iterative
 *                      2: Harmonic Eigen
 *                      3: Hybrid Biharmonic-Harmonic
 *                      4: Boundary Smoothed Harmonic
 * @param s_max_outer_iterations SEXP integer maximum number of global smoothing iterations
 * @param s_max_inner_iterations SEXP integer maximum number of iterations for each smoothing operation
 * @param s_smoothing_tolerance SEXP numeric convergence tolerance for smoothing algorithms
 * @param s_sigma SEXP numeric parameter controlling the smoothing kernel width
 * @param s_process_in_order SEXP logical whether to process extrema in ascending order of hop index
 * @param s_verbose SEXP logical whether to print progress information
 * @param s_detailed_recording SEXP logical whether to record detailed information about each
 *                           smoothing step for visualization
 *
 * @return SEXP containing a named list with:
 *         - harmonic_predictions: Vector of smoothed function values
 *         - lmin_hop_nbhds: List of hop neighborhoods for local minima
 *         - lmax_hop_nbhds: List of hop neighborhoods for local maxima
 *         - smoothing_history: (If detailed_recording=TRUE) List of records detailing
 *           each smoothing step
 *
 * @note All vertex indices are converted to 1-based indexing in the returned R structures
 *       for compatibility with R conventions.
 *
 * @see compute_extrema_hop_nbhds(), perform_weighted_mean_hop_disk_extension(),
 *      harmonic_extender(), harmonic_extension_eigen(),
 *      hybrid_biharmonic_harmonic_extension(), boundary_smoothed_harmonic_extension()
 */
SEXP S_create_gflow_cx(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_hop_idx_thld,
	SEXP s_smoother_type,
	SEXP s_max_outer_iterations,
	SEXP s_max_inner_iterations,
	SEXP s_smoothing_tolerance,
	SEXP s_sigma,
	SEXP s_process_in_order,
	SEXP s_verbose,
	SEXP s_detailed_recording
	) {
	// Convert parameters
	std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	double* y_ptr = REAL(s_y);
	std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

	size_t hop_idx_thld = static_cast<size_t>(asInteger(s_hop_idx_thld));
	smoother_type_t smoother_type = static_cast<smoother_type_t>(asInteger(s_smoother_type));
	int max_outer_iterations = asInteger(s_max_outer_iterations);
	int max_inner_iterations = asInteger(s_max_inner_iterations);
	double smoothing_tolerance = asReal(s_smoothing_tolerance);
	double sigma = asReal(s_sigma);
	bool process_in_order = asLogical(s_process_in_order);
	bool verbose = asLogical(s_verbose);
	bool detailed_recording = asLogical(s_detailed_recording);

	// Create graph and compute
	set_wgraph_t graph(adj_list, weight_list);

	gflow_cx_t result = graph.create_gflow_cx(
		y, hop_idx_thld, smoother_type,
		max_outer_iterations, max_inner_iterations,
		smoothing_tolerance, sigma,
		process_in_order, verbose, detailed_recording
		);

	// ---- Build R objects ----
	SEXP r_result = PROTECT(allocVector(VECSXP, detailed_recording ? 4 : 3));
	SEXP r_names = PROTECT(allocVector(STRSXP, detailed_recording ? 4 : 3));

	// Add names to result list
	SET_STRING_ELT(r_names, 0, mkChar("harmonic_predictions"));
	SET_STRING_ELT(r_names, 1, mkChar("lmin_hop_nbhds"));
	SET_STRING_ELT(r_names, 2, mkChar("lmax_hop_nbhds"));
	if (detailed_recording) {
		SET_STRING_ELT(r_names, 3, mkChar("smoothing_history"));
	}
	setAttrib(r_result, R_NamesSymbol, r_names);

	// Add predictions
	SEXP r_predictions = PROTECT(allocVector(REALSXP, result.harmonic_predictions.size()));
	double* pred_ptr = REAL(r_predictions);
	for (size_t i = 0; i < result.harmonic_predictions.size(); i++) {
		pred_ptr[i] = result.harmonic_predictions[i];
	}
	SET_VECTOR_ELT(r_result, 0, r_predictions);

	SEXP r_lmin_hop_nbhds = PROTECT(allocVector(VECSXP, result.lmin_hop_nbhd_map.size()));
	SEXP r_lmax_hop_nbhds = PROTECT(allocVector(VECSXP, result.lmax_hop_nbhd_map.size()));

	// Helper lambda to populate one neighborhood
	auto fill_nbhd = [&](const hop_nbhd_t& nbhd, SEXP r_nbhd) {
		// Names: vertex, hop_idx, nbhd_df, nbhd_bd_df
		SEXP names = PROTECT(allocVector(STRSXP, 5));
		SET_STRING_ELT(names, 0, mkChar("vertex"));
		SET_STRING_ELT(names, 1, mkChar("value"));
		SET_STRING_ELT(names, 2, mkChar("hop_idx"));
		SET_STRING_ELT(names, 3, mkChar("nbhd_df"));
		SET_STRING_ELT(names, 4, mkChar("nbhd_bd_df"));
		setAttrib(r_nbhd, R_NamesSymbol, names);

		bool is_global_extremum = (nbhd.hop_idx == std::numeric_limits<size_t>::max());

		// vertex (1-based)
		SEXP r_vertex = PROTECT(ScalarInteger(
									static_cast<int>(nbhd.vertex) + 1
									));

		// value
		SEXP r_value  = PROTECT(ScalarReal(y[nbhd.vertex]));

		// hop index
		SEXP r_hop_idx;
		if (is_global_extremum) {
			r_hop_idx = PROTECT(ScalarReal(R_PosInf));
		} else {
			r_hop_idx = PROTECT(ScalarReal(static_cast<double>(nbhd.hop_idx)));
		}

		// nbhd matrix
		SEXP r_nbhd_df;
		if (is_global_extremum) {
			r_nbhd_df = PROTECT(allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = nbhd.hop_dist_map.size();
			r_nbhd_df = PROTECT(allocMatrix(REALSXP, m, 2));
			double *pr = REAL(r_nbhd_df);
			size_t i = 0;
			for (const auto& [v,d] : nbhd.hop_dist_map) {
				pr[i]     = v + 1;
				pr[i + m] = d;
				i++;
			}
		}

		// nbhd_bd matrix
		SEXP r_nbhd_bd_df;
		if (is_global_extremum) {
			r_nbhd_bd_df = PROTECT(allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = nbhd.y_nbhd_bd_map.size();
			r_nbhd_bd_df = PROTECT(allocMatrix(REALSXP, m, 2));
			double *pr_bd = REAL(r_nbhd_bd_df);
			size_t i = 0;
			for (const auto& [v, y_val] : nbhd.y_nbhd_bd_map) {
				pr_bd[i]     = v + 1;
				pr_bd[i + m] = y_val;
				i++;
			}
		}

		// Set list entries
		SET_VECTOR_ELT(r_nbhd, 0, r_vertex);
		SET_VECTOR_ELT(r_nbhd, 1, r_value);
		SET_VECTOR_ELT(r_nbhd, 2, r_hop_idx);
		SET_VECTOR_ELT(r_nbhd, 3, r_nbhd_df);
		SET_VECTOR_ELT(r_nbhd, 4, r_nbhd_bd_df);

		UNPROTECT(6); // names, r_vertex, r_value, r_hop_idx, r_nbhd_df, r_nbhd_bd_df
	};

	// Fill minima nbhd's
	size_t i = 0;
	for (const auto& [_, nbhd] : result.lmin_hop_nbhd_map) {
		SEXP r_nbhd = PROTECT(allocVector(VECSXP, 5));
		fill_nbhd(nbhd, r_nbhd);
		SET_VECTOR_ELT(r_lmin_hop_nbhds, i++, r_nbhd);
		UNPROTECT(1); // r_nbhd
	}

	// Fill maxima nbhd's
	i = 0;
	for (const auto& [_, nbhd] : result.lmax_hop_nbhd_map) {
		SEXP r_nbhd = PROTECT(allocVector(VECSXP, 5));
		fill_nbhd(nbhd, r_nbhd);
		SET_VECTOR_ELT(r_lmax_hop_nbhds, i++, r_nbhd);
		UNPROTECT(1); // r_nbhd
	}

	SET_VECTOR_ELT(r_result, 1, r_lmin_hop_nbhds);
	SET_VECTOR_ELT(r_result, 2, r_lmax_hop_nbhds);


	// Add detailed smoothing history if requested
	if (detailed_recording) {
		int n_steps = result.smoothing_history.size();
		SEXP r_history = PROTECT(allocVector(VECSXP, n_steps));

		for (int i = 0; i < n_steps; i++) {
			auto& step = result.smoothing_history[i];

			// Create list for this step
			SEXP r_step = PROTECT(allocVector(VECSXP, 7));
			SEXP r_step_names = PROTECT(allocVector(STRSXP, 7));

			SET_STRING_ELT(r_step_names, 0, mkChar("vertex"));
			SET_STRING_ELT(r_step_names, 1, mkChar("is_minimum"));
			SET_STRING_ELT(r_step_names, 2, mkChar("hop_idx"));
			SET_STRING_ELT(r_step_names, 3, mkChar("smoother"));
			SET_STRING_ELT(r_step_names, 4, mkChar("before"));
			SET_STRING_ELT(r_step_names, 5, mkChar("after"));
			SET_STRING_ELT(r_step_names, 6, mkChar("region"));

			setAttrib(r_step, R_NamesSymbol, r_step_names);

			// Fill in the values
			SET_VECTOR_ELT(r_step, 0, ScalarInteger(step.vertex + 1));  // Convert to 1-indexed
			SET_VECTOR_ELT(r_step, 1, ScalarLogical(step.is_minimum));
			SET_VECTOR_ELT(r_step, 2, ScalarInteger(step.hop_idx));
			SET_VECTOR_ELT(r_step, 3, ScalarInteger(static_cast<int>(step.smoother)));

			// Convert before and after vectors
			SEXP r_before = PROTECT(allocVector(REALSXP, step.before.size()));
			SEXP r_after = PROTECT(allocVector(REALSXP, step.after.size()));

			for (size_t j = 0; j < step.before.size(); j++) {
				REAL(r_before)[j] = step.before[j];
				REAL(r_after)[j] = step.after[j];
			}

			SET_VECTOR_ELT(r_step, 4, r_before);
			SET_VECTOR_ELT(r_step, 5, r_after);

			// Convert region set to vector
			SEXP r_region = PROTECT(allocVector(INTSXP, step.region.size()));
			int* region_ptr = INTEGER(r_region);
			int k = 0;
			for (size_t v : step.region) {
				region_ptr[k++] = v + 1;  // Convert to 1-indexed
			}

			SET_VECTOR_ELT(r_step, 6, r_region);

			// Add step to history
			SET_VECTOR_ELT(r_history, i, r_step);

			UNPROTECT(5);  // r_step, r_step_names, r_before, r_after, r_region
		}

		SET_VECTOR_ELT(r_result, 3, r_history);
		UNPROTECT(1);  // r_history
	}

	UNPROTECT(5);  // r_result, r_names, r_predictions, r_lmin_hop_nbhds, r_lmax_hop_nbhds

	return r_result;
}


/**
 * @brief Applies a specified harmonic extension method to a local subgraph
 *
 * @details This function applies one of several harmonic extension methods to
 *          smooth a function over a specified region of vertices, with known
 *          boundary values. It can handle various methods including weighted mean,
 *          iterative harmonic, eigen-based harmonic, and biharmonic/hybrid extensions.
 *
 * @param s_adj_list SEXP representing the graph's adjacency list (from R)
 * @param s_weight_list SEXP representing the edge weights for the graph
 * @param s_y SEXP containing function values at each vertex as a numeric vector
 * @param s_region_vertices SEXP containing indices of interior vertices (region to smooth)
 * @param s_boundary_values SEXP containing named list of boundary vertex indices and values
 * @param s_smoother_type SEXP integer selecting the smoothing method to use
 * @param s_max_iterations SEXP integer maximum number of iterations for iterative methods
 * @param s_tolerance SEXP numeric convergence tolerance
 * @param s_sigma SEXP numeric parameter for weighted mean method
 * @param s_record_iterations SEXP logical whether to record iteration progress
 * @param s_verbose SEXP logical whether to print progress information
 *
 * @return SEXP containing a named list with:
 *         - smoothed_values: Vector of smoothed function values
 *         - iterations: List of function values at each iteration (if requested)
 *         - n_iterations: Number of iterations performed (for iterative methods)
 *         - max_change: Maximum change in values at final iteration
 *         - converged: Whether the iterative method converged
 */
SEXP S_apply_harmonic_extension(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_region_vertices,
	SEXP s_boundary_values,
	SEXP s_smoother_type,
	SEXP s_max_iterations,
	SEXP s_tolerance,
	SEXP s_sigma,
	SEXP s_record_iterations,
	SEXP s_verbose) {

	// ---- Input validation and conversion ----
	std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	// Convert numeric y vector
	double* y_ptr = REAL(s_y);
	std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

	// Convert region vertices (convert from 1-indexed R to 0-indexed C++)
	int* region_ptr = INTEGER(s_region_vertices);
	std::unordered_set<size_t> region_vertices;
	for (int i = 0; i < LENGTH(s_region_vertices); i++) {
		region_vertices.insert(static_cast<size_t>(region_ptr[i]));
	}

	// Convert boundary values
	// For this we expect a named list where names are vertex indices (as characters)
	// and values are the boundary values
	std::unordered_map<size_t, double> boundary_values;

	if (TYPEOF(s_boundary_values) == VECSXP) {
		SEXP names = Rf_getAttrib(s_boundary_values, R_NamesSymbol);
		if (names != R_NilValue) {
			for (int i = 0; i < LENGTH(s_boundary_values); i++) {
				// Get the vertex index from name (assuming it's the vertex index)
				size_t vertex = static_cast<size_t>(std::stoi(CHAR(STRING_ELT(names, i))));
				double value = REAL(VECTOR_ELT(s_boundary_values, i))[0];
				boundary_values[vertex] = value;
			}
		} else {
			REPORT_ERROR("Boundary values must be a named list with vertex indices as names");
		}
	} else {
		REPORT_ERROR("Boundary values must be a list");
	}

	// Get other parameters
	smoother_type_t smoother_type = static_cast<smoother_type_t>(asInteger(s_smoother_type));
	int max_iterations = asInteger(s_max_iterations);
	double tolerance = asReal(s_tolerance);
	double sigma = asReal(s_sigma);
	bool record_iterations = asLogical(s_record_iterations);
	bool verbose = asLogical(s_verbose);

	// ---- Create the graph and prepare for computation ----
	set_wgraph_t graph(adj_list, weight_list);

	// Initialize result vector with the original values
	std::vector<double> smoothed_values = y;

	// ---- Apply the requested harmonic extension method ----
	// Vector to store iterations if requested
	std::vector<std::vector<double>> iterations;
	int iterations_performed = 0;
	double max_change = 0.0;
	bool converged = false;

	switch (smoother_type) {
	case smoother_type_t::WMEAN: {
		// For weighted mean, we need to create a hop_nbhd_t structure
		// This is a bit of a simplification since we're receiving the region and boundary
		// directly rather than computing them from a hop index

		// Find the extremum vertex (just pick first vertex in region for simplicity)
		size_t vertex = *region_vertices.begin();

		// Create hop_nbhd_t structure
		hop_nbhd_t hop_nbhd;
		hop_nbhd.vertex = vertex;
		hop_nbhd.hop_idx = 1;  // Doesn't matter for weighted mean

		// Add all region vertices to hop_dist_map with arbitrary distances
		// (distances aren't used by weighted mean)
		for (size_t v : region_vertices) {
			hop_nbhd.hop_dist_map[v] = 0;  // All at "distance 0" for simplicity
		}

		// Add boundary values
		hop_nbhd.y_nbhd_bd_map = boundary_values;

		// Apply weighted mean extension
		if (record_iterations) {
			iterations.push_back(smoothed_values);
		}

		graph.perform_weighted_mean_hop_disk_extension(
			smoothed_values,
			hop_nbhd,
			max_iterations,
			tolerance,
			sigma,
			verbose
			);

		// Since this method doesn't return convergence info, we'll just say it completed
		iterations_performed = 1;
		converged = true;
		max_change = 0.0;

		if (record_iterations) {
			iterations.push_back(smoothed_values);
		}

		break;
	}

	case smoother_type_t::HARMONIC_IT: {
		// Apply iterative harmonic extension
		int record_frequency = record_iterations ? 1 : max_iterations;

		harmonic_extender_t result = graph.harmonic_extender(
			boundary_values,
			region_vertices,
			max_iterations,
			tolerance,
			record_frequency,
			verbose
			);

		// Get results
		if (!result.iterations.empty()) {
			smoothed_values = result.iterations.back();

			// Copy iterations if requested
			if (record_iterations) {
				iterations = result.iterations;
			}

			iterations_performed = result.num_iterations;
			converged = result.converged;
			max_change = result.max_change_final;
		}

		break;
	}

	case smoother_type_t::HARMONIC_EIGEN: {
		// Apply eigen-based harmonic extension
		double regularization = 1e-10;

		std::vector<double> result = graph.harmonic_extension_eigen(
			boundary_values,
			region_vertices,
			regularization,
			verbose
			);

		// Update smoothed values
		smoothed_values = result;

		// Since this is a direct solver, just record before/after
		if (record_iterations) {
			iterations.push_back(y);
			iterations.push_back(smoothed_values);
		}

		// This is a direct method so it always "converges" in one step
		iterations_performed = 1;
		converged = true;
		max_change = 0.0;

		break;
	}

	case smoother_type_t::HYBRID_BIHARMONIC_HARMONIC: {
		// Apply hybrid biharmonic-harmonic extension
		int boundary_blend_distance = 2;

		std::vector<double> result = graph.hybrid_biharmonic_harmonic_extension(
			boundary_values,
			region_vertices,
			boundary_blend_distance,
			verbose
			);

		// Update smoothed values
		smoothed_values = result;

		// Since this is a direct solver, just record before/after
		if (record_iterations) {
			iterations.push_back(y);
			iterations.push_back(smoothed_values);
		}

		// This is a direct method so it always "converges" in one step
		iterations_performed = 1;
		converged = true;
		max_change = 0.0;

		break;
	}

	case smoother_type_t::BOUNDARY_SMOOTHED_HARMONIC: {
		// Apply boundary smoothed harmonic extension
		int boundary_blend_distance = 2;

		std::vector<double> result = graph.boundary_smoothed_harmonic_extension(
			boundary_values,
			region_vertices,
			boundary_blend_distance,
			verbose
			);

		// Update smoothed values
		smoothed_values = result;

		// Since this is a direct solver, just record before/after
		if (record_iterations) {
			iterations.push_back(y);
			iterations.push_back(smoothed_values);
		}

		// This is a direct method so it always "converges" in one step
		iterations_performed = 1;
		converged = true;
		max_change = 0.0;

		break;
	}

	default:
		REPORT_ERROR("Unknown smoother type: %d", static_cast<int>(smoother_type));
	}

	// ---- Build R return value ----
	SEXP r_result = PROTECT(allocVector(VECSXP, 5));
	SEXP r_names = PROTECT(allocVector(STRSXP, 5));

	// Set names for the result list
	SET_STRING_ELT(r_names, 0, mkChar("smoothed_values"));
	SET_STRING_ELT(r_names, 1, mkChar("iterations"));
	SET_STRING_ELT(r_names, 2, mkChar("n_iterations"));
	SET_STRING_ELT(r_names, 3, mkChar("max_change"));
	SET_STRING_ELT(r_names, 4, mkChar("converged"));
	setAttrib(r_result, R_NamesSymbol, r_names);

	// Create and set smoothed_values
	SEXP r_smoothed_values = PROTECT(allocVector(REALSXP, smoothed_values.size()));
	for (size_t i = 0; i < smoothed_values.size(); i++) {
		REAL(r_smoothed_values)[i] = smoothed_values[i];
	}
	SET_VECTOR_ELT(r_result, 0, r_smoothed_values);

	// Create and set iterations
	SEXP r_iterations = PROTECT(allocVector(VECSXP, iterations.size()));
	for (size_t i = 0; i < iterations.size(); i++) {
		SEXP iteration = PROTECT(allocVector(REALSXP, iterations[i].size()));
		for (size_t j = 0; j < iterations[i].size(); j++) {
			REAL(iteration)[j] = iterations[i][j];
		}
		SET_VECTOR_ELT(r_iterations, i, iteration);
		UNPROTECT(1);  // iteration
	}
	SET_VECTOR_ELT(r_result, 1, r_iterations);

	// Set other values
	SET_VECTOR_ELT(r_result, 2, ScalarInteger(iterations_performed));
	SET_VECTOR_ELT(r_result, 3, ScalarReal(max_change));
	SET_VECTOR_ELT(r_result, 4, ScalarLogical(converged));

	UNPROTECT(4);  // r_result, r_names, r_smoothed_values, r_iterations

	return r_result;
}
