#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
// #include "error_utils.h"

#include <vector>

#include <R.h>
#include <Rinternals.h>
/**
 * @brief SEXP interface for computing gradient basins of attraction for local extrema.
 *
 * @details This function provides an R interface to the gradient basin computation.
 *          It identifies all local extrema in the graph and computes their gradient basins,
 *          returning separate lists for minima and maxima basins with complete trajectory
 *          information including predecessors and terminal extrema.
 *
 *          For each local extremum, the basin contains all vertices reachable via
 *          monotone paths (ascending for minima, descending for maxima), along with
 *          predecessor information for trajectory reconstruction and identification
 *          of terminal extrema where trajectories terminate.
 *
 * @param s_adj_list R list of integer vectors representing the adjacency list (0-based in C++)
 * @param s_weight_list R list of numeric vectors representing edge weights
 * @param s_y R numeric vector of function values at each vertex
 *
 * @return R list with two components:
 *         - lmin_basins: list of basin structures for local minima
 *         - lmax_basins: list of basin structures for local maxima
 *         Each basin structure contains:
 *         - vertex: 1-based vertex index of the extremum
 *         - value: function value at the extremum
 *         - hop_idx: maximum hop distance within the basin
 *         - basin_df: matrix with columns (vertex, hop_distance) for all basin members
 *         - basin_bd_df: matrix with columns (vertex, y_value) for boundary vertices
 *         - predecessors: integer vector where predecessors[i] is predecessor of vertex i (0 = no predecessor)
 *         - terminal_extrema: integer vector of terminal extrema vertices (1-based)
 */
extern "C" SEXP S_compute_basins_of_attraction(
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

	// ---- Build graph ----
	set_wgraph_t graph(adj_list, weight_list);

	// ---- Identify local extrema ----
	size_t n = y.size();
	std::vector<size_t> local_minima;
	std::vector<size_t> local_maxima;

	for (size_t v = 0; v < n; ++v) {
		bool is_local_min = true;
		bool is_local_max = true;

		for (const auto& edge : graph.adjacency_list[v]) {
			size_t u = edge.vertex;
			if (y[u] <= y[v]) {
				is_local_min = false;
			}
			if (y[u] >= y[v]) {
				is_local_max = false;
			}
		}

		if (is_local_min) {
			local_minima.push_back(v);
		}
		if (is_local_max) {
			local_maxima.push_back(v);
		}
	}

	// ---- Compute basins for all extrema ----
	std::vector<gradient_basin_t> lmin_basins;
	std::vector<gradient_basin_t> lmax_basins;

	for (size_t vertex : local_minima) {
		lmin_basins.push_back(graph.compute_basin_of_attraction(vertex, y, false));
	}

	for (size_t vertex : local_maxima) {
		lmax_basins.push_back(graph.compute_basin_of_attraction(vertex, y, true));
	}

	// ---- Build R objects ----
	SEXP r_lmin_basins = PROTECT(Rf_allocVector(VECSXP, lmin_basins.size()));
	SEXP r_lmax_basins = PROTECT(Rf_allocVector(VECSXP, lmax_basins.size()));

	// Helper lambda to populate one gradient basin
	auto fill_basin = [&](const gradient_basin_t& basin, SEXP r_basin) {
		// Names: vertex, value, hop_idx, basin_df, basin_bd_df, predecessors, terminal_extrema
		SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
		SET_STRING_ELT(names, 0, Rf_mkChar("vertex"));
		SET_STRING_ELT(names, 1, Rf_mkChar("value"));
		SET_STRING_ELT(names, 2, Rf_mkChar("hop_idx"));
		SET_STRING_ELT(names, 3, Rf_mkChar("basin_df"));
		SET_STRING_ELT(names, 4, Rf_mkChar("basin_bd_df"));
		SET_STRING_ELT(names, 5, Rf_mkChar("predecessors"));
		SET_STRING_ELT(names, 6, Rf_mkChar("terminal_extrema"));
		Rf_setAttrib(r_basin, R_NamesSymbol, names);

		bool is_empty_basin = (basin.hop_idx == std::numeric_limits<size_t>::max());

		// vertex (1-based)
		SEXP r_vertex = PROTECT(Rf_ScalarInteger(
									static_cast<int>(basin.vertex) + 1
									));

		// value
		SEXP r_value  = PROTECT(Rf_ScalarReal(basin.value));

		// hop index
		SEXP r_hop_idx;
		if (is_empty_basin) {
			r_hop_idx = PROTECT(Rf_ScalarReal(R_NaN));
		} else {
			r_hop_idx = PROTECT(Rf_ScalarReal(static_cast<double>(basin.hop_idx)));
		}

		// basin matrix (vertex, hop_distance)
		SEXP r_basin_df;
		if (is_empty_basin) {
			r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = basin.hop_dist_map.size();
			r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
			double *pr = REAL(r_basin_df);
			size_t i = 0;
			for (const auto& [v, d] : basin.hop_dist_map) {
				pr[i]     = v + 1;  // convert to 1-based
				pr[i + m] = d;
				i++;
			}
		}

		// basin_bd matrix (vertex, y_value)
		SEXP r_basin_bd_df;
		if (is_empty_basin) {
			r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2)); // empty matrix
		} else {
			size_t m = basin.y_nbhd_bd_map.size();
			r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
			double *pr_bd = REAL(r_basin_bd_df);
			size_t i = 0;
			for (const auto& [v, y_val] : basin.y_nbhd_bd_map) {
				pr_bd[i]     = v + 1;  // convert to 1-based
				pr_bd[i + m] = y_val;
				i++;
			}
		}

		// predecessors vector
		// Create a vector of size n where predecessors[i] = predecessor of vertex i
		// Use 0 to indicate no predecessor (since R uses 1-based indexing)
		SEXP r_predecessors;
		if (is_empty_basin) {
			r_predecessors = PROTECT(Rf_allocVector(INTSXP, n));
			int *p_pred = INTEGER(r_predecessors);
			for (size_t i = 0; i < n; ++i) {
				p_pred[i] = 0;  // no predecessor
			}
		} else {
			r_predecessors = PROTECT(Rf_allocVector(INTSXP, n));
			int *p_pred = INTEGER(r_predecessors);

			// Initialize all to 0 (no predecessor)
			for (size_t i = 0; i < n; ++i) {
				p_pred[i] = 0;
			}

			// Fill in predecessors for basin vertices
			for (const auto& [v, pred_v] : basin.predecessors) {
				if (pred_v != std::numeric_limits<size_t>::max()) {
					p_pred[v] = static_cast<int>(pred_v) + 1;  // convert to 1-based
				}
				// if pred_v is max(), leave as 0 (no predecessor - this is the origin)
			}
		}

		// terminal_extrema vector
		SEXP r_terminal_extrema;
		if (is_empty_basin || basin.terminal_extrema.empty()) {
			r_terminal_extrema = PROTECT(Rf_allocVector(INTSXP, 0)); // empty vector
		} else {
			size_t m = basin.terminal_extrema.size();
			r_terminal_extrema = PROTECT(Rf_allocVector(INTSXP, m));
			int *p_term = INTEGER(r_terminal_extrema);
			for (size_t i = 0; i < m; ++i) {
				p_term[i] = static_cast<int>(basin.terminal_extrema[i]) + 1;  // convert to 1-based
			}
		}

		// Set list entries
		SET_VECTOR_ELT(r_basin, 0, r_vertex);
		SET_VECTOR_ELT(r_basin, 1, r_value);
		SET_VECTOR_ELT(r_basin, 2, r_hop_idx);
		SET_VECTOR_ELT(r_basin, 3, r_basin_df);
		SET_VECTOR_ELT(r_basin, 4, r_basin_bd_df);
		SET_VECTOR_ELT(r_basin, 5, r_predecessors);
		SET_VECTOR_ELT(r_basin, 6, r_terminal_extrema);

		UNPROTECT(8); // names, r_vertex, r_value, r_hop_idx, r_basin_df, r_basin_bd_df, r_predecessors, r_terminal_extrema
	};

	// Fill minima basins
	for (size_t i = 0; i < lmin_basins.size(); ++i) {
		SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 7));
		fill_basin(lmin_basins[i], r_basin);
		SET_VECTOR_ELT(r_lmin_basins, i, r_basin);
		UNPROTECT(1); // r_basin
	}

	// Fill maxima basins
	for (size_t i = 0; i < lmax_basins.size(); ++i) {
		SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 7));
		fill_basin(lmax_basins[i], r_basin);
		SET_VECTOR_ELT(r_lmax_basins, i, r_basin);
		UNPROTECT(1); // r_basin
	}

	// Name the two top-level components
	SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
	SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, 2));
	SET_STRING_ELT(r_result_names, 0, Rf_mkChar("lmin_basins"));
	SET_STRING_ELT(r_result_names, 1, Rf_mkChar("lmax_basins"));
	Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
	SET_VECTOR_ELT(r_result, 0, r_lmin_basins);
	SET_VECTOR_ELT(r_result, 1, r_lmax_basins);

	UNPROTECT(4); // r_lmin_basins, r_lmax_basins, r_result, r_result_names

	return r_result;
}
