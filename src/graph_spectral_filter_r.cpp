#include "set_wgraph.hpp"               // set_wgraph_t
#include "graph_spectral_filter.hpp"    // graph_spectral_filter_t
#include "SEXP_cpp_conversion_utils.hpp"// convert_adj_list_from_R, convert_weight_list_from_R, Rvect_to_CppVect_double

#include <R.h>
#include <Rinternals.h>

extern "C" SEXP S_graph_spectral_filter(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_laplacian_type,
	SEXP s_filter_type,
	SEXP s_laplacian_power,
	SEXP s_kernel_tau_factor,
	SEXP s_kernel_radius_factor,
	SEXP s_kernel_type,
	SEXP s_kernel_adaptive,
	SEXP s_min_radius_factor,
	SEXP s_max_radius_factor,
	SEXP s_domain_min_size,
	SEXP s_precision,
	SEXP s_n_evectors_to_compute,
	SEXP s_n_candidates,
	SEXP s_log_grid,
	SEXP s_with_t_predictions,
	SEXP s_verbose
	) {

	// 1) Convert inputs
	std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

	// Convert enum types
	laplacian_type_t laplacian_type = static_cast<laplacian_type_t>(Rf_asInteger(s_laplacian_type));
	filter_type_t filter_type = static_cast<filter_type_t>(Rf_asInteger(s_filter_type));

	// Convert numeric parameters
	size_t laplacian_power = static_cast<size_t>(Rf_asInteger(s_laplacian_power));

	// Setup kernel parameters
	kernel_params_t kernel_params;
	kernel_params.tau_factor        = Rf_asReal(s_kernel_tau_factor);
	kernel_params.radius_factor     = Rf_asReal(s_kernel_radius_factor);
	kernel_params.kernel_type       = static_cast<kernel_type_t>(Rf_asInteger(s_kernel_type));
	kernel_params.adaptive          = (Rf_asLogical(s_kernel_adaptive) == TRUE);
	kernel_params.min_radius_factor = Rf_asReal(s_min_radius_factor);
	kernel_params.max_radius_factor = Rf_asReal(s_max_radius_factor);
	kernel_params.domain_min_size   = static_cast<size_t>(Rf_asInteger(s_domain_min_size));
	kernel_params.precision         = Rf_asReal(s_precision);

	// Convert remaining parameters
	size_t n_evectors_to_compute = static_cast<size_t>(Rf_asInteger(s_n_evectors_to_compute));
	size_t n_candidates          = static_cast<size_t>(Rf_asInteger(s_n_candidates));
	bool log_grid                = (Rf_asLogical(s_log_grid) == TRUE);
	bool with_t_predictions      = (Rf_asLogical(s_with_t_predictions) == TRUE);
	bool verbose                 = (Rf_asLogical(s_verbose) == TRUE);

	// 2) Call the C++ graph filter
	set_wgraph_t graph(adj_list, weight_list);
	graph.compute_graph_diameter();

	graph_spectral_filter_t result = graph.graph_spectral_filter(
		y,
		laplacian_type,
		filter_type,
		laplacian_power,
		kernel_params,
		n_evectors_to_compute,
		n_candidates,
		log_grid,
		with_t_predictions,
		verbose
		);

	// 3) Build the R list to return
	// Determine how many fields to return based on whether we include all predictions for each t value
	int n_fields = with_t_predictions ? 13 : 12;
	SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_fields));

	// Set names for the result list
	{
		SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
		SET_STRING_ELT(names, 0, Rf_mkChar("evalues"));
		SET_STRING_ELT(names, 1, Rf_mkChar("evectors"));
		SET_STRING_ELT(names, 2, Rf_mkChar("candidate_ts"));
		SET_STRING_ELT(names, 3, Rf_mkChar("gcv_scores"));
		SET_STRING_ELT(names, 4, Rf_mkChar("opt_t_idx"));
		SET_STRING_ELT(names, 5, Rf_mkChar("predictions"));
		SET_STRING_ELT(names, 6, Rf_mkChar("t_predictions"));
		SET_STRING_ELT(names, 7, Rf_mkChar("laplacian_type"));
		SET_STRING_ELT(names, 8, Rf_mkChar("filter_type"));
		SET_STRING_ELT(names, 9, Rf_mkChar("laplacian_power"));
		SET_STRING_ELT(names, 10, Rf_mkChar("kernel_params"));
		SET_STRING_ELT(names, 11, Rf_mkChar("compute_time_ms"));
		if (with_t_predictions) {
			SET_STRING_ELT(names, 12, Rf_mkChar("gcv_min_score"));
		}
		Rf_setAttrib(r_result, R_NamesSymbol, names);
		UNPROTECT(1); // names
	}


	// 3a evalues
	{
		SEXP Revalues = PROTECT(Rf_allocVector(REALSXP, result.evalues.size()));
		for (int i = 0; i < (int)result.evalues.size(); ++i)
			REAL(Revalues)[i] = result.evalues[i];
		SET_VECTOR_ELT(r_result, 0, Revalues);
		UNPROTECT(1); // Revalues
	}

	// 3b) evectors
	{
		int nrow = result.evectors.rows();
		int ncol = result.evectors.cols();
		SEXP Revectors = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
		for (int c = 0; c < ncol; ++c)
			for (int r = 0; r < nrow; ++r)
				REAL(Revectors)[r + nrow * c] = result.evectors(r, c);
		SET_VECTOR_ELT(r_result, 1, Revectors);
		UNPROTECT(1); // Revectors
	}

	// 3c) candidate_ts
	{
		SEXP Rts = PROTECT(Rf_allocVector(REALSXP, result.candidate_ts.size()));
		for (int i = 0; i < (int)result.candidate_ts.size(); ++i)
			REAL(Rts)[i] = result.candidate_ts[i];
		SET_VECTOR_ELT(r_result, 2, Rts);
		UNPROTECT(1);
	}

	// 3d) gcv_scores
	{
		SEXP Rgcv = PROTECT(Rf_allocVector(REALSXP, result.gcv_scores.size()));
		for (int i = 0; i < (int)result.gcv_scores.size(); ++i)
			REAL(Rgcv)[i] = result.gcv_scores[i];
		SET_VECTOR_ELT(r_result, 3, Rgcv);
		UNPROTECT(1);
	}

	// 3e) opt_t_idx
	{
		SEXP Ropt = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Ropt)[0] = static_cast<int>(result.opt_t_idx + 1);   // 1‐based for R
		SET_VECTOR_ELT(r_result, 4, Ropt);
		UNPROTECT(1);
	}

	// 3f) predictions
	{
		SEXP Rpred = PROTECT(Rf_allocVector(REALSXP, result.predictions.size()));
		for (int i = 0; i < (int)result.predictions.size(); ++i)
			REAL(Rpred)[i] = result.predictions[i];
		SET_VECTOR_ELT(r_result, 5, Rpred);
		UNPROTECT(1);
	}

	// 3g) t_predictions (optional)
	{
		if (with_t_predictions) {
			SEXP Rtp = PROTECT(Rf_allocMatrix(
								   REALSXP,
								   result.t_predictions[0].size(),
								   result.t_predictions.size()
								   ));
			int nvert = result.t_predictions[0].size();
			for (int j = 0; j < (int)result.t_predictions.size(); ++j) {
				for (int i = 0; i < nvert; ++i) {
					REAL(Rtp)[i + nvert * j] = result.t_predictions[j][i];
				}
			}
			SET_VECTOR_ELT(r_result, 6, Rtp);
			UNPROTECT(1);
		} else {
			SET_VECTOR_ELT(r_result, 6, R_NilValue);
		}

	}

	// 3h) laplacian_type
	{
		SEXP Rlap_type = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Rlap_type)[0] = static_cast<int>(result.laplacian_type);
		SET_VECTOR_ELT(r_result, 7, Rlap_type);
		UNPROTECT(1);
	}

	// 3i) filter_type
	{
		SEXP Rfilter_type = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Rfilter_type)[0] = static_cast<int>(result.filter_type);
		SET_VECTOR_ELT(r_result, 8, Rfilter_type);
		UNPROTECT(1);
	}

	// 3j) laplacian_power
	{
		SEXP Rlap_power = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Rlap_power)[0] = static_cast<int>(result.laplacian_power);
		SET_VECTOR_ELT(r_result, 9, Rlap_power);
		UNPROTECT(1);
	}

	// 3k) kernel_params
	{
		SEXP Rkernel_params = PROTECT(Rf_allocVector(VECSXP, 8));

		SEXP Rtau_factor = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rtau_factor)[0] = result.kernel_params.tau_factor;
		SET_VECTOR_ELT(Rkernel_params, 0, Rtau_factor);

		SEXP Rradius_factor = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rradius_factor)[0] = result.kernel_params.radius_factor;
		SET_VECTOR_ELT(Rkernel_params, 1, Rradius_factor);

		SEXP Rkernel_type = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Rkernel_type)[0] = static_cast<int>(result.kernel_params.kernel_type);
		SET_VECTOR_ELT(Rkernel_params, 2, Rkernel_type);

		SEXP Radaptive = PROTECT(Rf_allocVector(LGLSXP, 1));
		LOGICAL(Radaptive)[0] = result.kernel_params.adaptive;
		SET_VECTOR_ELT(Rkernel_params, 3, Radaptive);

		SEXP Rmin_radius_factor = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rmin_radius_factor)[0] = result.kernel_params.min_radius_factor;
		SET_VECTOR_ELT(Rkernel_params, 4, Rmin_radius_factor);

		SEXP Rmax_radius_factor = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rmax_radius_factor)[0] = result.kernel_params.max_radius_factor;
		SET_VECTOR_ELT(Rkernel_params, 5, Rmax_radius_factor);

		SEXP Rdomain_min_size = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(Rdomain_min_size)[0] = static_cast<int>(result.kernel_params.domain_min_size);
		SET_VECTOR_ELT(Rkernel_params, 6, Rdomain_min_size);

		SEXP Rprecision = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rprecision)[0] = result.kernel_params.precision;
		SET_VECTOR_ELT(Rkernel_params, 7, Rprecision);

		// Set names for kernel_params
		SEXP kernel_params_names = PROTECT(Rf_allocVector(STRSXP, 8));
		SET_STRING_ELT(kernel_params_names, 0, Rf_mkChar("tau_factor"));
		SET_STRING_ELT(kernel_params_names, 1, Rf_mkChar("radius_factor"));
		SET_STRING_ELT(kernel_params_names, 2, Rf_mkChar("kernel_type"));
		SET_STRING_ELT(kernel_params_names, 3, Rf_mkChar("adaptive"));
		SET_STRING_ELT(kernel_params_names, 4, Rf_mkChar("min_radius_factor"));
		SET_STRING_ELT(kernel_params_names, 5, Rf_mkChar("max_radius_factor"));
		SET_STRING_ELT(kernel_params_names, 6, Rf_mkChar("domain_min_size"));
		SET_STRING_ELT(kernel_params_names, 7, Rf_mkChar("precision"));
		Rf_setAttrib(Rkernel_params, R_NamesSymbol, kernel_params_names);

		SET_VECTOR_ELT(r_result, 10, Rkernel_params);
		UNPROTECT(1);
	}

	// 3l) compute_time_ms
	{
		SEXP Rcompute_time = PROTECT(Rf_allocVector(REALSXP, 1));
		REAL(Rcompute_time)[0] = result.compute_time_ms;
		SET_VECTOR_ELT(r_result, 11, Rcompute_time);
		UNPROTECT(1);
	}

	// 3m) gcv_min_score (if with_t_predictions, this is element 12)
	{
		if (with_t_predictions) {
			SEXP Rgcv_min = PROTECT(Rf_allocVector(REALSXP, 1));
			REAL(Rgcv_min)[0] = result.gcv_min_score;
			SET_VECTOR_ELT(r_result, 12, Rgcv_min);
			UNPROTECT(1);
		}
	}

	UNPROTECT(1);
	return r_result;
}
