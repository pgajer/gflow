#include <R.h>
#include <Rinternals.h>

#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
	SEXP S_construct_function_aware_graph(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_function_values,
		SEXP s_weight_type,
		SEXP s_epsilon,
		SEXP s_lambda,
		SEXP s_alpha,
		SEXP s_beta,
		SEXP s_tau,
		SEXP s_p,
		SEXP s_q,
		SEXP s_r,
		SEXP s_normalize,
		SEXP s_weight_thld
		);

	SEXP S_analyze_function_aware_weights(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_function_values,
		SEXP s_weight_types,
		SEXP s_epsilon,
		SEXP s_lambda,
		SEXP s_alpha,
		SEXP s_beta,
		SEXP s_tau,
		SEXP s_p,
		SEXP s_q,
		SEXP s_r
		);
}

/**
 * R interface for construct_function_aware_graph
 */
SEXP S_construct_function_aware_graph(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_function_values,
	SEXP s_weight_type,
	SEXP s_epsilon,
	SEXP s_lambda,
	SEXP s_alpha,
	SEXP s_beta,
	SEXP s_tau,
	SEXP s_p,
	SEXP s_q,
	SEXP s_r,
	SEXP s_normalize,
	SEXP s_weight_thld
	) {

	// Convert inputs
	std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	// Convert function values with proper memory management
	std::unique_ptr<std::vector<double>> function_values_ptr = Rvect_to_CppVect_double(s_function_values);
	if (!function_values_ptr) {
		Rf_error("Failed to convert function values");
	}
	std::vector<double> function_values = std::move(*function_values_ptr);

	int weight_type = INTEGER(s_weight_type)[0];
	double epsilon = REAL(s_epsilon)[0];
	double lambda = REAL(s_lambda)[0];
	double alpha = REAL(s_alpha)[0];
	double beta = REAL(s_beta)[0];
	double tau = REAL(s_tau)[0];
	double p = REAL(s_p)[0];
	double q = REAL(s_q)[0];
	double r = REAL(s_r)[0];
	bool normalize = LOGICAL(s_normalize)[0];
	double weight_thld = REAL(s_weight_thld)[0];

	// Create graph and compute function-aware version
	set_wgraph_t graph(adj_list, weight_list);
	set_wgraph_t result_graph = graph.construct_function_aware_graph(
		function_values,
		weight_type,
		epsilon,
		lambda,
		alpha,
		beta,
		tau,
		p,
		q,
		r,
		normalize,
		weight_thld
		);

	// Convert result to R format
	return convert_wgraph_to_R(result_graph);
}

/**
 * @brief Analyzes weight distributions for function-aware graph construction
 *
 * @details
 * This function takes a graph and function values defined on its vertices, and analyzes
 * how different weighting schemes would modify the edge weights. It returns the complete
 * distribution of modified weights for several weighting schemes, which can be used to
 * determine appropriate thresholds for edge pruning.
 *
 * The function supports various weighting schemes:
 * - Inverse relationship: w_new = w_old / (|f(i) - f(j)| + epsilon)
 * - Direct relationship: w_new = w_old * |f(i) - f(j)|
 * - Exponential decay: w_new = w_old * exp(-lambda * |f(i) - f(j)|)
 * - Power law: w_new = w_old * |f(i) - f(j)|^(-alpha)
 * - Sigmoid: w_new = w_old * 1/(1 + exp(beta * (|f(i) - f(j)| - tau)))
 *
 * @param s_adj_list SEXP R list of adjacency lists
 * @param s_weight_list SEXP R list of edge weights
 * @param s_function_values SEXP R numeric vector of function values at each vertex
 * @param s_weight_type SEXP R integer specifying the weighting scheme to analyze (0-4)
 * @param s_params SEXP R numeric vector of parameters for the weighting function
 *
 * @return SEXP A named R list containing weight distributions for each weighting scheme
 *         The first element contains the original, unmodified weights.
 *
 * @note The caller is responsible for unprotecting the returned SEXP
 */
SEXP S_analyze_function_aware_weights(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_function_values,
	SEXP s_weight_types,
	SEXP s_epsilon,
	SEXP s_lambda,
	SEXP s_alpha,
	SEXP s_beta,
	SEXP s_tau,
	SEXP s_p,
	SEXP s_q,
	SEXP s_r
	) {

	// Validate inputs
	if (!Rf_isVector(s_weight_types) || !Rf_isInteger(s_weight_types)) {
		Rf_error("weight_types must be an integer vector");
	}

	// Convert inputs
	std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	// Convert function values with proper memory management
	std::unique_ptr<std::vector<double>> function_values_ptr = Rvect_to_CppVect_double(s_function_values);
	if (!function_values_ptr) {
		Rf_error("Failed to convert function values");
	}
	std::vector<double> function_values = std::move(*function_values_ptr);

	// Convert weight types to vector
	int* weight_types_ptr = INTEGER(s_weight_types);
	int weight_types_len = LENGTH(s_weight_types);
	std::vector<int> weight_types(weight_types_ptr, weight_types_ptr + weight_types_len);

	// Get other parameters
	double epsilon = REAL(s_epsilon)[0];
	double lambda = REAL(s_lambda)[0];
	double alpha = REAL(s_alpha)[0];
	double beta = REAL(s_beta)[0];
	double tau = REAL(s_tau)[0];
	double p = REAL(s_p)[0];
	double q = REAL(s_q)[0];
	double r = REAL(s_r)[0];

	// Create graph
	set_wgraph_t graph(adj_list, weight_list);

	// Check that function values match the number of vertices
	if (function_values.size() != graph.num_vertices()) {
		Rf_error("Length of function_values must match the number of vertices");
	}

	// Compute weight distributions
	std::vector<std::vector<double>> results = graph.analyze_function_aware_weights(
		function_values,
		weight_types,
		epsilon,
		lambda,
		alpha,
		beta,
		tau,
		p,
		q,
		r
		);

	// Create return list - now one element larger to include original weights
	int nprot = 0;
	SEXP r_list = PROTECT(Rf_allocVector(VECSXP, results.size())); nprot++;

	// Set names
	const char* type_names[] = {"original", "inverse", "direct", "exp", "power", "sigmoid", "lp_embedding"};
	SEXP r_list_names = PROTECT(Rf_allocVector(STRSXP, results.size())); nprot++;

	// Set first element as original weights
	SET_STRING_ELT(r_list_names, 0, Rf_mkChar("original"));

	// Convert original weights to R and add to list
	SEXP r_orig_weights = PROTECT(Rf_allocVector(REALSXP, results[0].size()));
	for (size_t j = 0; j < results[0].size(); j++) {
		REAL(r_orig_weights)[j] = results[0][j];
	}
	SET_VECTOR_ELT(r_list, 0, r_orig_weights);
	UNPROTECT(1); // r_orig_weights

	// Process the rest of the weight types
	for (size_t i = 0; i < weight_types.size(); i++) {
		int type = weight_types[i];
		if (type >= 0 && type <= 5) {
			SET_STRING_ELT(r_list_names, i + 1, Rf_mkChar(type_names[type + 1]));
		} else {
			SET_STRING_ELT(r_list_names, i + 1, Rf_mkChar("unknown"));
		}

		// Convert each weight distribution to R and add to list
		SEXP r_weights = PROTECT(Rf_allocVector(REALSXP, results[i + 1].size()));
		for (size_t j = 0; j < results[i + 1].size(); j++) {
			REAL(r_weights)[j] = results[i + 1][j];
		}
		SET_VECTOR_ELT(r_list, i + 1, r_weights);
		UNPROTECT(1); // r_weights
	}

	Rf_setAttrib(r_list, R_NamesSymbol, r_list_names);
	UNPROTECT(nprot); // r_list and r_list_names

	return r_list;
}
