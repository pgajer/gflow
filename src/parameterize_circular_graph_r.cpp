#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
	SEXP S_parameterize_circular_graph(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_use_edge_lengths
		);
}


/**
 * @brief Parameterizes a circular graph structure using spectral methods and returns the result to R
 *
 * This function takes adjacency list and weight list from R, constructs a weighted graph,
 * applies spectral parameterization, and returns the result to R as a named list.
 *
 * The parameterization is done using eigenvectors of the graph Laplacian to embed the
 * graph in a circle. The resulting angles represent the positions of vertices on the circle.
 *
 * @param s_adj_list SEXP containing an R list of integer vectors representing adjacency lists
 * @param s_weight_list SEXP containing an R list of numeric vectors representing edge weights
 * @param s_use_edge_lengths SEXP containing a logical value indicating whether to use edge weights
 *
 * @return SEXP containing a list with three elements:
 *         - angles: numeric vector of angles (in radians) for each vertex
 *         - eig_vec2: numeric vector containing the second eigenvector
 *         - eig_vec3: numeric vector containing the third eigenvector
 *
 * @note The function assumes that the input adjacency list and weight list have the same structure
 *       and that the graph is connected.
 */
SEXP S_parameterize_circular_graph(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_use_edge_lengths
	) {

	// Convert input parameters using R's C API
	std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
	bool use_edge_lengths = LOGICAL(s_use_edge_lengths)[0];

	set_wgraph_t graph(adj_list, weight_list);

	circular_param_result_t res = graph.parameterize_circular_graph(use_edge_lengths);

	// Create the return list using R's C API
	const char* names[] = {
		"angles",
		"eig_vec2",
		"eig_vec3",
		"eig_vec4",
		"eig_vec5",
		"eig_vec6",
		NULL
	};

	int n_elements = 0;
	while (names[n_elements] != NULL) n_elements++;

	// Create list and protect it
	int protect_count = 0;
	SEXP result = PROTECT(allocVector(VECSXP, n_elements)); protect_count++;

	// Set names
	SEXP result_names = PROTECT(allocVector(STRSXP, n_elements));
	for (int i = 0; i < n_elements; i++) {
		SET_STRING_ELT(result_names, i, mkChar(names[i]));
	}
	setAttrib(result, R_NamesSymbol, result_names);
	UNPROTECT(1); // for result_names

	// Helper function to convert vector to SEXP
	auto create_numeric_vector = [](const std::vector<double>& vec) -> SEXP {
		SEXP r_vec = PROTECT(allocVector(REALSXP, vec.size()));
		double* ptr = REAL(r_vec);
		std::copy(vec.begin(), vec.end(), ptr);
		return r_vec;
	};

	SET_VECTOR_ELT(result, 0, create_numeric_vector(res.angles)); protect_count++;
	SET_VECTOR_ELT(result, 1, create_numeric_vector(res.eig_vec2)); protect_count++;
	SET_VECTOR_ELT(result, 2, create_numeric_vector(res.eig_vec3)); protect_count++;
	SET_VECTOR_ELT(result, 3, create_numeric_vector(res.eig_vec4)); protect_count++;
	SET_VECTOR_ELT(result, 4, create_numeric_vector(res.eig_vec5)); protect_count++;
	SET_VECTOR_ELT(result, 5, create_numeric_vector(res.eig_vec6)); protect_count++;

	UNPROTECT(protect_count);

	return result;
}
