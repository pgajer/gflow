// R interface function

#include "nerve_cx.hpp"
#include "error_utils.h"   // For REPORT_ERROR

extern "C" {
#include "nerve_cx_r.h"
}

// Helper function to get an element from an R list by name
// This is a common utility function used in R C interfaces
static SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < LENGTH(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    }
    return elmt;
}

/**
 * @brief R interface for nerve_cx_spectral_filter function
 */
SEXP S_nerve_cx_spectral_filter(
	SEXP s_complex_ptr,
	SEXP s_y,
	SEXP s_laplacian_type,
	SEXP s_filter_type,
	SEXP s_laplacian_power,
	SEXP s_dim_weights,
	SEXP s_kernel_params,
	SEXP s_n_evectors,
	SEXP s_n_candidates,
	SEXP s_log_grid,
	SEXP s_with_t_predictions,
	SEXP s_verbose
	) {
	// Get the nerve complex
	nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
	if (complex == NULL) {
		REPORT_ERROR("Invalid nerve complex pointer");
	}

	// Convert input parameters
	PROTECT(s_y = Rf_coerceVector(s_y, REALSXP));
	int n_vertices = LENGTH(s_y);
	double* y_ptr = REAL(s_y);

	if (n_vertices != (int)complex->num_vertices()) {
		REPORT_ERROR("Length of y (%d) does not match number of vertices (%d)",
					 n_vertices, complex->num_vertices());
	}

	laplacian_type_t laplacian_type = static_cast<laplacian_type_t>(INTEGER(s_laplacian_type)[0]);
	filter_type_t filter_type = static_cast<filter_type_t>(INTEGER(s_filter_type)[0]);
	size_t laplacian_power = INTEGER(s_laplacian_power)[0];

	PROTECT(s_dim_weights = Rf_coerceVector(s_dim_weights, REALSXP));
	int n_dim_weights = LENGTH(s_dim_weights);
	double* dim_weights_ptr = REAL(s_dim_weights);

	std::vector<double> dim_weights(dim_weights_ptr, dim_weights_ptr + n_dim_weights);

	// Process kernel parameters
	kernel_params_t kernel_params;

	SEXP s_tau_factor = getListElement(s_kernel_params, "tau_factor");
	SEXP s_radius_factor = getListElement(s_kernel_params, "radius_factor");
	SEXP s_kernel_type = getListElement(s_kernel_params, "kernel_type");

	if (s_tau_factor != R_NilValue) {
		kernel_params.tau_factor = REAL(s_tau_factor)[0];
	}

	if (s_radius_factor != R_NilValue) {
		kernel_params.radius_factor = REAL(s_radius_factor)[0];
	}

	if (s_kernel_type != R_NilValue) {
		kernel_params.kernel_type = static_cast<kernel_type_t>(INTEGER(s_kernel_type)[0]);
	}

	size_t n_evectors = INTEGER(s_n_evectors)[0];
	size_t n_candidates = INTEGER(s_n_candidates)[0];
	bool log_grid = LOGICAL(s_log_grid)[0];
	bool with_t_predictions = LOGICAL(s_with_t_predictions)[0];
	bool verbose = LOGICAL(s_verbose)[0];

	// Convert y to std::vector
	std::vector<double> y(y_ptr, y_ptr + n_vertices);

	// Call the function
	nerve_cx_spectral_filter_t result = complex->nerve_cx_spectral_filter(
		y, laplacian_type, filter_type, laplacian_power, dim_weights,
		kernel_params, n_evectors, n_candidates, log_grid,
		with_t_predictions, verbose
		);

	// Prepare return list structure
	SEXP ret = PROTECT(Rf_allocVector(VECSXP, 7));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));

	// 1. Predictions
	SEXP predictions = PROTECT(Rf_allocVector(REALSXP, n_vertices));
	double* pred_ptr = REAL(predictions);
	for (size_t i = 0; i < n_vertices; ++i) {
		pred_ptr[i] = result.predictions[i];
	}
	SET_VECTOR_ELT(ret, 0, predictions);
	SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));

	// 2. Optimal parameter
	SEXP opt_param = PROTECT(Rf_ScalarReal(result.candidate_ts[result.opt_t_idx]));
	SET_VECTOR_ELT(ret, 1, opt_param);
	SET_STRING_ELT(names, 1, Rf_mkChar("optimal_parameter"));

	// 3. GCV score
	SEXP gcv_score = PROTECT(Rf_ScalarReal(result.gcv_min_score));
	SET_VECTOR_ELT(ret, 2, gcv_score);
	SET_STRING_ELT(names, 2, Rf_mkChar("gcv_score"));

	// 4. Computation time
	SEXP compute_time = PROTECT(Rf_ScalarReal(result.compute_time_ms));
	SET_VECTOR_ELT(ret, 3, compute_time);
	SET_STRING_ELT(names, 3, Rf_mkChar("compute_time_ms"));

	// 5. All parameters
	SEXP all_params = PROTECT(Rf_allocVector(REALSXP, result.candidate_ts.size()));
	double* params_ptr = REAL(all_params);
	for (size_t i = 0; i < result.candidate_ts.size(); ++i) {
		params_ptr[i] = result.candidate_ts[i];
	}
	SET_VECTOR_ELT(ret, 4, all_params);
	SET_STRING_ELT(names, 4, Rf_mkChar("all_parameters"));

	// 6. All GCV scores
	SEXP all_gcv = PROTECT(Rf_allocVector(REALSXP, result.gcv_scores.size()));
	double* gcv_ptr = REAL(all_gcv);
	for (size_t i = 0; i < result.gcv_scores.size(); ++i) {
		gcv_ptr[i] = result.gcv_scores[i];
	}
	SET_VECTOR_ELT(ret, 5, all_gcv);
	SET_STRING_ELT(names, 5, Rf_mkChar("all_gcv_scores"));

	// 7. All predictions (if requested)
	if (with_t_predictions) {
		SEXP all_predictions = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, result.t_predictions.size()));
		double* all_pred_ptr = REAL(all_predictions);

		for (size_t j = 0; j < result.t_predictions.size(); ++j) {
			for (size_t i = 0; i < n_vertices; ++i) {
				all_pred_ptr[i + j * n_vertices] = result.t_predictions[j][i];
			}
		}

		SET_VECTOR_ELT(ret, 6, all_predictions);
		SET_STRING_ELT(names, 6, Rf_mkChar("all_predictions"));
		UNPROTECT(1);
	}
	else {
		SET_VECTOR_ELT(ret, 6, R_NilValue);
		SET_STRING_ELT(names, 6, Rf_mkChar("all_predictions"));
	}

	// Set the list names
	Rf_setAttrib(ret, R_NamesSymbol, names);

	UNPROTECT(10);  // Unprotect everything
	return ret;
}
