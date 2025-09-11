#include "mst_completion_graphs.hpp"
#include "error_utils.h"                 // REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp" // convert_adj_list_from_R, convert_weight_list_from_R

#include <vector>
#include <algorithm>      // std::copy

#include <R.h>            // For Rprintf, etc.
#include <Rinternals.h>   // For SEXP macros

/**
 * @brief Extracts edge weights from a set_wgraph_t into a numeric R vector
 *
 * This function returns only the weights of unique edges (undirected),
 * avoiding duplicates from symmetric adjacency lists.
 */
SEXP get_graph_edge_weights(const set_wgraph_t& G) {
	auto& adj = G.adjacency_list;
	size_t n = adj.size();
	std::vector<double> weights;

	for (size_t u = 0; u < n; ++u) {
		for (const auto& nbr : adj[u]) {
			size_t v = nbr.vertex;
			if (u < v) { // avoid double-counting
				weights.push_back(nbr.weight);
			}
		}
	}

	SEXP res = PROTECT(Rf_allocVector(REALSXP, weights.size()));
	std::copy(weights.begin(), weights.end(), REAL(res));

	UNPROTECT(1);
	return res;
}

/**
 * Convert a set_wgraph_t to a 3-column R matrix [from, to, weight]
 */
SEXP graph_to_edge_matrix(const set_wgraph_t& G) {
	auto& adj = G.adjacency_list;
	size_t n = adj.size();

	std::vector<std::tuple<int, int, double>> edges;

	for (size_t u = 0; u < n; ++u) {
		for (const auto& nbr : adj[u]) {
			size_t v = nbr.vertex;
			double w = nbr.weight;
			if (u < v) {
				edges.emplace_back(u + 1, v + 1, w);  // 1-based indexing for R
			}
		}
	}

	size_t m = edges.size();
	SEXP res = PROTECT(Rf_allocMatrix(REALSXP, m, 3));
	double* ptr = REAL(res);
	for (size_t i = 0; i < m; ++i) {
		ptr[i] = std::get<0>(edges[i]);  // from
		ptr[i + m] = std::get<1>(edges[i]);  // to
		ptr[i + 2 * m] = std::get<2>(edges[i]);  // weight
	}
	UNPROTECT(1);
	return res;
}

#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <utility>  // for std::pair

/**
 * @brief Converts a set_wgraph_t to a pair of R lists:
 *        - adjacency list: integer vectors (1-based indices)
 *        - weight list: numeric vectors (edge weights)
 *
 * @param G A set_wgraph_t object representing the graph
 * @return std::pair<SEXP, SEXP> where first = adj list, second = weight list
 */
std::pair<SEXP, SEXP>
create_r_graph_from_set_wgraph(const set_wgraph_t& G) {
	const auto& adj_list = G.adjacency_list;
	size_t n_vertices = adj_list.size();

	// Prepare temporary C++ storage for adjacency and weights
	std::vector<std::vector<int>> adj_vect(n_vertices);
	std::vector<std::vector<double>> dist_vect(n_vertices);

	for (size_t i = 0; i < n_vertices; ++i) {
		for (const auto& nbr : adj_list[i]) {
			adj_vect[i].push_back(static_cast<int>(nbr.vertex));
			dist_vect[i].push_back(nbr.weight);
		}
	}

	SEXP r_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
	SEXP r_weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));

	for (size_t i = 0; i < n_vertices; ++i) {
		// Integer vector for neighbors
		SEXP RA = PROTECT(Rf_allocVector(INTSXP, adj_vect[i].size()));
		int* A = INTEGER(RA);
		for (size_t j = 0; j < adj_vect[i].size(); ++j) {
			A[j] = adj_vect[i][j] + 1; // 1-based indexing
		}
		SET_VECTOR_ELT(r_adj_list, i, RA);
		UNPROTECT(1); // RA

		// Numeric vector for weights
		SEXP RD = PROTECT(Rf_allocVector(REALSXP, dist_vect[i].size()));
		double* D = REAL(RD);
		for (size_t j = 0; j < dist_vect[i].size(); ++j) {
			D[j] = dist_vect[i][j];
		}
		SET_VECTOR_ELT(r_weight_list, i, RD);
		UNPROTECT(1); // RD
	}

	return std::make_pair(r_adj_list, r_weight_list);
}

/**
 * @brief R interface to compute MST and its completion graph using a quantile-based threshold.
 *
 * This function is a C-level interface callable from R via `.Call()`. It takes a numeric matrix
 * of data points, computes the minimal spanning tree (MST) using Prim's algorithm via `C_mstree`,
 * and adds extra edges between pairs of points whose Euclidean distance is below a specified quantile
 * threshold of MST edge lengths.
 *
 * @param s_X A numeric R matrix (type REALSXP), where each row is a data point (n_points x n_dims).
 * @param s_q_thld A numeric scalar (REALSXP) between 0 and 1 specifying the quantile threshold.
 * @param s_verbose A logical scalar (LGLSXP) controlling verbosity of output.
 *
 * @return An R list with the following components:
 *   - mst_adj_list: A list of integer vectors, where the i-th element contains the 1-based indices
 *                   of vertices adjacent to vertex i in the MST.
 *   - mst_weight_list: A list of numeric vectors, where the i-th element contains the weights of
 *                      edges from vertex i in the MST (in the same order as mst_adj_list).
 *   - cmst_adj_list: A list of integer vectors representing the adjacency list of the completed MST
 *                    (i.e., MST + additional short edges based on a quantile threshold).
 *   - cmst_weight_list: A list of numeric vectors representing edge weights for the completed MST
 *                       adjacency structure.
 *   - mst_edge_weights: A numeric vector of MST edge weights, in the order they were added during
 *                       Prim's algorithm.
 *
 * @throws Rf_error if input types or sizes are invalid.
 *
 * @note This function uses only the base R C API and does not depend on Rcpp.
 *
 * Example (in R):
 * \code{r}
 *   res <- .Call("S_create_mst_completion_graph", X, 0.9, TRUE)
 * \endcode
 */
extern "C" SEXP S_create_mst_completion_graph(
	SEXP s_X,
	SEXP s_q_thld,
	SEXP s_verbose
	) {
	// Input checking
	if (!Rf_isReal(s_X) || !Rf_isMatrix(s_X)) {
		REPORT_ERROR("X must be a numeric matrix");
	}
	if (!Rf_isReal(s_q_thld) || Rf_length(s_q_thld) != 1) {
		REPORT_ERROR("q_thld must be a numeric scalar");
	}
	if (!Rf_isLogical(s_verbose) || Rf_length(s_verbose) != 1) {
		REPORT_ERROR("verbose must be a logical scalar");
	}

	double q_thld = REAL(s_q_thld)[0];
	bool verbose = (LOGICAL(s_verbose)[0] == 1);

	std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));

	mst_completion_graph_t res = create_mst_completion_graph(X, q_thld, verbose);

	// Create R return list
	size_t nprot = 0;
	auto [r_mst_adj_list, r_mst_weights_list] = create_r_graph_from_set_wgraph(res.mstree); nprot += 2;
	auto [r_cmst_adj_list, r_cmst_weights_list] = create_r_graph_from_set_wgraph(res.completed_mstree); nprot += 2;

	size_t m = res.mstree_edge_weights.size();
	SEXP r_mst_weights = PROTECT(Rf_allocVector(REALSXP, m)); nprot++;
	std::copy(res.mstree_edge_weights.begin(), res.mstree_edge_weights.end(), REAL(r_mst_weights));

		// Output list components
	SEXP r_list = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
	SET_VECTOR_ELT(r_list, 0, r_mst_adj_list);
	SET_VECTOR_ELT(r_list, 1, r_mst_weights_list);
	SET_VECTOR_ELT(r_list, 2, r_cmst_adj_list);
	SET_VECTOR_ELT(r_list, 3, r_cmst_weights_list);
	SET_VECTOR_ELT(r_list, 4, r_mst_weights);

	SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
	SET_STRING_ELT(r_names, 0, Rf_mkChar("mst_adj_list"));
	SET_STRING_ELT(r_names, 1, Rf_mkChar("mst_weight_list"));
	SET_STRING_ELT(r_names, 2, Rf_mkChar("cmst_adj_list"));
	SET_STRING_ELT(r_names, 3, Rf_mkChar("cmst_weight_list"));
	SET_STRING_ELT(r_names, 4, Rf_mkChar("mst_edge_weights"));
	Rf_setAttrib(r_list, R_NamesSymbol, r_names);

	UNPROTECT(nprot);

	return r_list;
}

