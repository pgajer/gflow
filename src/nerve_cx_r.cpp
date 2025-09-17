#include "nerve_cx.hpp"
#include "simplex_weights.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "error_utils.h"   // For REPORT_ERROR

#include <Rinternals.h>
#include <vector>
#include <stdexcept>

extern "C" {
	SEXP S_create_nerve_complex(SEXP s_X, SEXP s_k, SEXP s_max_dim);
    SEXP S_extract_skeleton_graph(SEXP s_complex_ptr);
	SEXP S_set_function_values(SEXP s_complex_ptr, SEXP s_values);
	SEXP S_set_weight_scheme(SEXP s_complex_ptr, SEXP s_weight_type, SEXP s_params);
	SEXP S_solve_full_laplacian(SEXP s_complex_ptr, SEXP s_lambda, SEXP s_dim_weights);
	SEXP S_get_simplex_counts(SEXP s_complex_ptr);
	SEXP S_extract_skeleton_graph(SEXP s_complex_ptr);
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
    // ---- Validate X: numeric double matrix with dim=c(n_points, n_dims)
    if (TYPEOF(s_X) != REALSXP)
        Rf_error("X must be a numeric (double) matrix.");
    SEXP s_dim = PROTECT(Rf_getAttrib(s_X, R_DimSymbol)); // +1
    if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) != 2) {
        UNPROTECT(1);
        Rf_error("X must be a numeric matrix with a valid 'dim' of length 2.");
    }
    const int n_points = INTEGER(s_dim)[0];
    const int n_dims   = INTEGER(s_dim)[1];
    if (n_points <= 0 || n_dims <= 0) {
        UNPROTECT(1);
        Rf_error("'X' must have positive dimensions.");
    }
    UNPROTECT(1); // s_dim

    // ---- Validate k and max_dim
    if (!Rf_isInteger(s_k) || !Rf_isInteger(s_max_dim))
        Rf_error("'k' and 'max_dim' must be integer.");
    const int k       = Rf_asInteger(s_k);
    const int max_dim = Rf_asInteger(s_max_dim);
    if (k < 2)        Rf_error("'k' must be at least 2.");
    if (max_dim < 1)  Rf_error("'max_dim' must be at least 1.");

    // ---- Copy coordinates (column-major R -> row-major C++)
    const double* X = REAL(s_X);
    std::vector<std::vector<double>> coords(n_points, std::vector<double>(n_dims));
    for (int j = 0; j < n_dims; ++j) {
        const double* col = X + (size_t)j * (size_t)n_points;
        for (int i = 0; i < n_points; ++i) {
            coords[i][j] = col[i];
        }
    }

    // ---- Allocate result list and set names
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 4));            // +1
    SEXP names    = PROTECT(Rf_allocVector(STRSXP, 4));            // +2
    SET_STRING_ELT(names, 0, Rf_mkChar("complex_ptr"));
    SET_STRING_ELT(names, 1, Rf_mkChar("n_vertices"));
    SET_STRING_ELT(names, 2, Rf_mkChar("simplex_counts"));
    SET_STRING_ELT(names, 3, Rf_mkChar("max_dimension"));
    Rf_setAttrib(r_result, R_NamesSymbol, names);
    UNPROTECT(1); // names -> +1

    // ---- Build complex (C++ heap), wrap in external pointer, register finalizer
    nerve_complex_t* nc = nullptr;
    try {
        nc = new nerve_complex_t(coords, k, max_dim);
    } catch (...) {
        UNPROTECT(1); // r_result
        Rf_error("Failed to construct nerve complex.");
    }

    SEXP complex_ptr = PROTECT(R_MakeExternalPtr(nc, R_NilValue, R_NilValue)); // +2
    R_RegisterCFinalizerEx(complex_ptr, (R_CFinalizer_t)nerve_complex_finalizer, TRUE);
    SET_VECTOR_ELT(r_result, 0, complex_ptr);
    UNPROTECT(1); // complex_ptr (now safely referenced by r_result) -> +1

    // ---- n_vertices
    SEXP n_vertices_sexp = PROTECT(Rf_ScalarInteger(n_points));    // +2
    SET_VECTOR_ELT(r_result, 1, n_vertices_sexp);
    UNPROTECT(1); // n_vertices_sexp -> +1

    // ---- simplex_counts: length max_dim + 1
    SEXP simplex_counts = PROTECT(Rf_allocVector(INTSXP, max_dim + 1)); // +2
    int* counts = INTEGER(simplex_counts);
    for (int d = 0; d <= max_dim; ++d) {
        counts[d] = nc->num_simplices(d);
    }
    SET_VECTOR_ELT(r_result, 2, simplex_counts);
    UNPROTECT(1); // simplex_counts -> +1

    // ---- max_dimension
    SEXP max_dim_sexp = PROTECT(Rf_ScalarInteger(max_dim));        // +2
    SET_VECTOR_ELT(r_result, 3, max_dim_sexp);
    UNPROTECT(1); // max_dim_sexp -> +1

    UNPROTECT(1); // r_result -> +0
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
