#include "nerve_cx.hpp"
#include "simplex_weights.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "error_utils.h"   // For REPORT_ERROR

extern "C" {
#include "nerve_cx_r.h"
}

/**
 * @brief Finalizer for nerve complex objects
 */
void nerve_complex_finalizer(SEXP s_complex_ptr) {
	if (NULL == R_ExternalPtrAddr(s_complex_ptr)) return;
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	delete complex;
	R_ClearExternalPtr(s_complex_ptr);
}

/**
 * @brief Creates a nerve complex from point coordinates
 *
 * @param s_X SEXP matrix of point coordinates
 * @param s_k SEXP integer for k nearest neighbors
 * @param s_max_dim SEXP integer for maximum simplex dimension
 * @return SEXP external pointer to the nerve complex
 */
SEXP S_create_nerve_complex(SEXP s_X, SEXP s_k, SEXP s_max_dim) {
	// Convert inputs
	PROTECT(s_X = Rf_coerceVector(s_X, REALSXP));
	int* dimX = INTEGER(Rf_getAttrib(s_X, R_DimSymbol));
	int n_points = dimX[0];
	int n_dims = dimX[1];

	int k       = Rf_asInteger(s_k);
	int max_dim = Rf_asInteger(s_max_dim);

	// Convert coordinates to C++ format
	double* X = REAL(s_X);
	std::vector<std::vector<double>> coords(n_points);
	for (int i = 0; i < n_points; i++) {
		coords[i].resize(n_dims);
		for (int j = 0; j < n_dims; j++) {
			coords[i][j] = X[i + n_points * j]; // Column-major to row-major
		}
	}
	UNPROTECT(1); // s_X

	// Create return list
	SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 4));
	{
		// Set names
		SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
		SET_STRING_ELT(names, 0, Rf_mkChar("complex_ptr"));
		SET_STRING_ELT(names, 1, Rf_mkChar("n_vertices"));
		SET_STRING_ELT(names, 2, Rf_mkChar("simplex_counts"));
		SET_STRING_ELT(names, 3, Rf_mkChar("max_dimension"));
		Rf_setAttrib(r_result, R_NamesSymbol, names);
		UNPROTECT(1);
	}

	// Create the nerve complex
	{
		nerve_complex_t* complex = new nerve_complex_t(coords, k, max_dim);

		// Create external pointer
		SEXP complex_ptr = R_MakeExternalPtr(complex, R_NilValue, R_NilValue);

		// Register finalizer
		R_RegisterCFinalizerEx(complex_ptr, (R_CFinalizer_t) nerve_complex_finalizer, TRUE);

		SET_VECTOR_ELT(r_result, 0, complex_ptr);
		UNPROTECT(1);
	}

	SET_VECTOR_ELT(r_result, 1, Rf_ScalarInteger(n_points));

	// Get simplex counts
	{
		SEXP simplex_counts = PROTECT(Rf_allocVector(INTSXP, max_dim + 1));
		int* counts = INTEGER(simplex_counts);
		for (int d = 0; d <= max_dim; d++) {
			counts[d] = complex->num_simplices(d);
		}
		SET_VECTOR_ELT(r_result, 2, simplex_counts);
		UNPROTECT(1);
	}

	SET_VECTOR_ELT(r_result, 3, Rf_ScalarInteger(max_dim));

	UNPROTECT(1);
	return r_result;
}

/**
 * @brief Sets function values on vertices of the nerve complex
 */
SEXP S_set_function_values(SEXP s_complex_ptr, SEXP s_values) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Convert values
	PROTECT(s_values = Rf_coerceVector(s_values, REALSXP));
	int n_values = LENGTH(s_values);
	double* values = REAL(s_values);

	if (n_values != (int)complex->num_vertices()) {
		REPORT_ERROR("Number of values (%d) does not match number of vertices (%d)",
			  n_values, complex->num_vertices());
	}

	// Set function values
	std::vector<double> function_values(n_values);
	for (int i = 0; i < n_values; i++) {
		function_values[i] = values[i];
	}
	complex->set_function_values(function_values);

	UNPROTECT(1);
	return R_NilValue;
}

/**
 * @brief Sets weight scheme for the nerve complex
 */
SEXP S_set_weight_scheme(SEXP s_complex_ptr, SEXP s_weight_type, SEXP s_params) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Convert weight type
	const char* weight_type = CHAR(STRING_ELT(s_weight_type, 0));

	// Convert parameters
	PROTECT(s_params = Rf_coerceVector(s_params, REALSXP));
	int n_params = LENGTH(s_params);
	double* params = REAL(s_params);

	// Set weight functions
	size_t max_dim = complex->max_dimension;

	if (strcmp(weight_type, "uniform") == 0) {
		for (size_t d = 0; d <= max_dim; d++) {
			complex->set_weight_function(d, simplex_weights::create_weight_function(simplex_weights::uniform_weight));
		}
	}
	else if (strcmp(weight_type, "inverse_distance") == 0) {
		for (size_t d = 0; d <= max_dim; d++) {
			complex->set_weight_function(d, simplex_weights::create_weight_function(simplex_weights::inverse_distance_weight));
		}
	}
	else if (strcmp(weight_type, "gaussian") == 0) {
		if (n_params < 1) {
			REPORT_ERROR("Gaussian weight scheme requires sigma parameter");
		}
		double sigma = params[0];
		for (size_t d = 0; d <= max_dim; d++) {
			complex->set_weight_function(d, simplex_weights::create_gaussian_weight(sigma));
		}
	}
	else if (strcmp(weight_type, "volume") == 0) {
		if (n_params < 1) {
			REPORT_ERROR("Volume weight scheme requires alpha parameter");
		}
		double alpha = params[0];
		for (size_t d = 0; d <= max_dim; d++) {
			complex->set_weight_function(d, simplex_weights::create_volume_weight(alpha));
		}
	}
	else if (strcmp(weight_type, "gradient") == 0) {
		if (n_params < 1) {
			REPORT_ERROR("Gradient weight scheme requires gamma parameter");
		}
		double gamma = params[0];
		for (size_t d = 0; d <= max_dim; d++) {
			complex->set_weight_function(d, simplex_weights::create_gradient_weight(gamma));
		}
	}
	else {
		REPORT_ERROR("Unknown weight type: %s", weight_type);
	}

	// Update weights
	complex->update_weights();

	UNPROTECT(1);
	return R_NilValue;
}

/**
 * @brief Solves the full Laplacian problem
 */
SEXP S_solve_full_laplacian(SEXP s_complex_ptr, SEXP s_lambda, SEXP s_dim_weights) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Convert lambda
	double lambda = REAL(s_lambda)[0];

	// Convert dimension weights
	PROTECT(s_dim_weights = Rf_coerceVector(s_dim_weights, REALSXP));
	int n_weights = LENGTH(s_dim_weights);
	double* weights = REAL(s_dim_weights);
	UNPROTECT(1); // s_dim_weights

	if (n_weights < (int)complex->max_dimension + 1) {
		REPORT_ERROR("Not enough dimension weights provided");
	}

	// Convert weights to C++ vector
	std::vector<double> dim_weights(weights, weights + n_weights);

	// Solve the full Laplacian problem
	std::vector<double> result = complex->solve_full_laplacian(lambda, dim_weights);

	// Convert result to R vector
	SEXP r_result = PROTECT(Rf_allocVector(REALSXP, result.size()));
	double* r_result_ptr = REAL(r_result);
	for (size_t i = 0; i < result.size(); i++) {
		r_result_ptr[i] = result[i];
	}

	UNPROTECT(1); // r_result
	return r_result;
}

/**
 * @brief Gets simplex counts for each dimension
 */
SEXP S_get_simplex_counts(SEXP s_complex_ptr) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Get simplex counts
	SEXP counts = PROTECT(Rf_allocVector(INTSXP, complex->max_dimension + 1));
	int* counts_ptr = INTEGER(counts);
	for (size_t d = 0; d <= complex->max_dimension; d++) {
		counts_ptr[d] = complex->num_simplices(d);
	}

	UNPROTECT(1);
	return counts;
}

/**
 * @brief Extracts the 1-skeleton as a graph
 */
SEXP S_extract_skeleton_graph(SEXP s_complex_ptr) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Get the skeleton
	if (!complex->skeleton) {
		REPORT_ERROR("Nerve complex does not have a valid skeleton");
	}

	// Convert to R format
	return convert_wgraph_to_R(*(complex->skeleton));
}
