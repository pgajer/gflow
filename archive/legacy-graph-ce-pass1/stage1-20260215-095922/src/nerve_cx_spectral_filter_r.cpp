// R interface function

#include "nerve_cx.hpp"
#include "error_utils.h"   // For REPORT_ERROR

// Helper function to get an element from an R list by name
static SEXP getListElement(SEXP list, const char* key) {
    if (TYPEOF(list) != VECSXP || key == NULL) return R_NilValue;

    SEXP names = PROTECT(Rf_getAttrib(list, R_NamesSymbol)); // +1
    if (names == R_NilValue || TYPEOF(names) != STRSXP || LENGTH(names) != LENGTH(list)) {
        UNPROTECT(1);
        return R_NilValue;
    }

    const int n = LENGTH(list);
    for (int i = 0; i < n; ++i) {
        SEXP nm = STRING_ELT(names, i);
        if (nm == R_NilValue) continue;
        const char* cname = CHAR(nm); // safe given names is STRSXP
        if (cname && strcmp(cname, key) == 0) {
            SEXP ans = VECTOR_ELT(list, i); // borrowing, no PROTECT required for return
            UNPROTECT(1); // names
            return ans;
        }
    }
    UNPROTECT(1); // names
    return R_NilValue;
}

/**
 * @brief R interface for nerve_cx_spectral_filter function
 */
extern "C" SEXP S_nerve_cx_spectral_filter(
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
    // ---- Validate external pointer
    if (TYPEOF(s_complex_ptr) != EXTPTRSXP)
        Rf_error("complex_ptr must be an external pointer.");
    nerve_complex_t* complex = (nerve_complex_t*) R_ExternalPtrAddr(s_complex_ptr);
    if (!complex) Rf_error("Invalid nerve complex pointer.");

    // ---- Validate and read y
    if (TYPEOF(s_y) != REALSXP)
        Rf_error("'y' must be a numeric (double) vector.");
    const int n_vertices = LENGTH(s_y);
    if (n_vertices <= 0) Rf_error("'y' must have positive length.");
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_vertices);

    if (n_vertices != (int) complex->num_vertices()) {
        Rf_error("Length of y (%d) does not match number of vertices (%d).",
                 n_vertices, (int)complex->num_vertices());
    }

    // ---- Scalar params (LENGTH-first; reject wrong types)
    if (!Rf_isInteger(s_laplacian_type) || !Rf_isInteger(s_filter_type) || !Rf_isInteger(s_laplacian_power))
        Rf_error("'laplacian_type', 'filter_type', and 'laplacian_power' must be integer.");
    if (!Rf_isInteger(s_n_evectors) || !Rf_isInteger(s_n_candidates))
        Rf_error("'n_evectors' and 'n_candidates' must be integer.");
    if (TYPEOF(s_log_grid) != LGLSXP || TYPEOF(s_with_t_predictions) != LGLSXP || TYPEOF(s_verbose) != LGLSXP)
        Rf_error("'log_grid', 'with_t_predictions', and 'verbose' must be logical.");

    const laplacian_type_t laplacian_type = (laplacian_type_t) Rf_asInteger(s_laplacian_type);
    const filter_type_t    filter_type    = (filter_type_t)    Rf_asInteger(s_filter_type);
    const int              lap_power_i    = Rf_asInteger(s_laplacian_power);
    if (lap_power_i < 0) Rf_error("'laplacian_power' must be nonnegative.");
    const size_t laplacian_power = (size_t) lap_power_i;

    const int n_evectors  = Rf_asInteger(s_n_evectors);
    const int n_candidates= Rf_asInteger(s_n_candidates);
    if (n_evectors <= 0 || n_candidates <= 0)
        Rf_error("'n_evectors' and 'n_candidates' must be positive.");

    const int log_grid          = LOGICAL(s_log_grid)[0];
    const int with_t_predictions= LOGICAL(s_with_t_predictions)[0];
    const int verbose           = LOGICAL(s_verbose)[0];

    // ---- dim_weights
    if (TYPEOF(s_dim_weights) != REALSXP)
        Rf_error("'dim_weights' must be a numeric (double) vector.");
    const int n_dim_weights = LENGTH(s_dim_weights);
    if (n_dim_weights <= 0) Rf_error("'dim_weights' must have positive length.");
    std::vector<double> dim_weights(REAL(s_dim_weights), REAL(s_dim_weights) + n_dim_weights);

    // ---- kernel params list
    if (TYPEOF(s_kernel_params) != VECSXP)
        Rf_error("'kernel_params' must be a named list.");
    // Safe fetches (return R_NilValue if missing)
    SEXP s_tau_factor    = getListElement(s_kernel_params, "tau_factor");
    SEXP s_radius_factor = getListElement(s_kernel_params, "radius_factor");
    SEXP s_kernel_type   = getListElement(s_kernel_params, "kernel_type");

    kernel_params_t kernel_params{};
    if (s_tau_factor != R_NilValue) {
        if (TYPEOF(s_tau_factor) != REALSXP || LENGTH(s_tau_factor) < 1)
            Rf_error("'tau_factor' must be numeric length>=1.");
        kernel_params.tau_factor = REAL(s_tau_factor)[0];
    }
    if (s_radius_factor != R_NilValue) {
        if (TYPEOF(s_radius_factor) != REALSXP || LENGTH(s_radius_factor) < 1)
            Rf_error("'radius_factor' must be numeric length>=1.");
        kernel_params.radius_factor = REAL(s_radius_factor)[0];
    }
    if (s_kernel_type != R_NilValue) {
        if (!Rf_isInteger(s_kernel_type) || LENGTH(s_kernel_type) < 1)
            Rf_error("'kernel_type' must be integer length>=1.");
        kernel_params.kernel_type = (kernel_type_t) INTEGER(s_kernel_type)[0];
    }

    // ---- Call core
    nerve_cx_spectral_filter_t result =
        complex->nerve_cx_spectral_filter(
            y, laplacian_type, filter_type, laplacian_power, dim_weights,
            kernel_params, (size_t)n_evectors, (size_t)n_candidates,
            log_grid != 0, with_t_predictions != 0, verbose != 0
        );

    // ---- Build return list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 7));         // +1
    SEXP names    = PROTECT(Rf_allocVector(STRSXP, 7));         // +2
    SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 1, Rf_mkChar("optimal_parameter"));
    SET_STRING_ELT(names, 2, Rf_mkChar("gcv_score"));
    SET_STRING_ELT(names, 3, Rf_mkChar("compute_time_ms"));
    SET_STRING_ELT(names, 4, Rf_mkChar("all_parameters"));
    SET_STRING_ELT(names, 5, Rf_mkChar("all_gcv_scores"));
    SET_STRING_ELT(names, 6, Rf_mkChar("all_predictions"));
    Rf_setAttrib(r_result, R_NamesSymbol, names);
    UNPROTECT(1); // names -> +1

    // 1) predictions
    {
        SEXP predictions = PROTECT(Rf_allocVector(REALSXP, n_vertices)); // +2
        double* p = REAL(predictions);
        for (int i = 0; i < n_vertices; ++i) p[i] = result.predictions[(size_t)i];
        SET_VECTOR_ELT(r_result, 0, predictions);
        UNPROTECT(1); // +1
    }

    // 2) optimal parameter
    {
        const int opt_idx = (int) result.opt_t_idx;
        if (opt_idx < 0 || opt_idx >= (int)result.candidate_ts.size())
            Rf_error("Internal error: opt_t_idx out of range.");
        SEXP opt_param = PROTECT(Rf_ScalarReal(result.candidate_ts[(size_t)opt_idx])); // +2
        SET_VECTOR_ELT(r_result, 1, opt_param);
        UNPROTECT(1); // +1
    }

    // 3) gcv score
    {
        SEXP gcv_score = PROTECT(Rf_ScalarReal(result.gcv_min_score)); // +2
        SET_VECTOR_ELT(r_result, 2, gcv_score);
        UNPROTECT(1); // +1
    }

    // 4) compute time
    {
        SEXP compute_time = PROTECT(Rf_ScalarReal(result.compute_time_ms)); // +2
        SET_VECTOR_ELT(r_result, 3, compute_time);
        UNPROTECT(1); // +1
    }

    // 5) all parameters
    {
        const int m = (int) result.candidate_ts.size();
        SEXP all_params = PROTECT(Rf_allocVector(REALSXP, m)); // +2
        double* p = REAL(all_params);
        for (int i = 0; i < m; ++i) p[i] = result.candidate_ts[(size_t)i];
        SET_VECTOR_ELT(r_result, 4, all_params);
        UNPROTECT(1); // +1
    }

    // 6) all gcv scores
    {
        const int m = (int) result.gcv_scores.size();
        SEXP all_gcv = PROTECT(Rf_allocVector(REALSXP, m)); // +2
        double* p = REAL(all_gcv);
        for (int i = 0; i < m; ++i) p[i] = result.gcv_scores[(size_t)i];
        SET_VECTOR_ELT(r_result, 5, all_gcv);
        UNPROTECT(1); // +1
    }

    // 7) all predictions (optional)
    if (with_t_predictions) {
        const int nT = (int) result.t_predictions.size();
        SEXP all_predictions = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, nT)); // +2
        double* A = REAL(all_predictions);
        for (int j = 0; j < nT; ++j) {
            const std::vector<double>& col = result.t_predictions[(size_t)j];
            if ((int)col.size() != n_vertices)
                Rf_error("Internal error: t_predictions[%d] has wrong length.", j);
            // column-major fill
            for (int i = 0; i < n_vertices; ++i)
                A[i + (size_t)j * (size_t)n_vertices] = col[(size_t)i];
        }
        SET_VECTOR_ELT(r_result, 6, all_predictions);
        UNPROTECT(1); // +1
    } else {
        SET_VECTOR_ELT(r_result, 6, R_NilValue);
    }

    UNPROTECT(1); // r_result -> +0
    return r_result;
}
