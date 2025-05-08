#include <R.h>
#include <Rinternals.h>
#undef length             // avoid clash with R

#include "set_wgraph.hpp"                   // set_wgraph_t
#include "klaps_spectral_filter.hpp"   // klaps_spectral_filter_t
#include "SEXP_cpp_conversion_utils.hpp"    // convert_adj_list_from_R, convert_weight_list_from_R, Rvect_to_CppVect_double

extern "C" SEXP S_klaps_spectral_filter(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_filter_type,
	SEXP s_n_evectors_to_compute,
	SEXP s_tau_factor,
	SEXP s_radius_factor,
	SEXP s_laplacian_power,
	SEXP s_n_candidates,
	SEXP s_log_grid,
	SEXP s_with_t_predictions,
	SEXP s_verbose
	) {

	// 1) Convert inputs
	std::vector<std::vector<int>> adj_list =
		convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list =
		convert_weight_list_from_R(s_weight_list);

	const double* y_ptr = REAL(s_y);
	R_xlen_t n = XLENGTH(s_y);
	std::vector<double> y(y_ptr, y_ptr + n);

	filter_type_t filter_type    = static_cast<filter_type_t>(INTEGER(s_filter_type)[0]);
	size_t n_evectors_to_compute = static_cast<size_t>(INTEGER(s_n_evectors_to_compute)[0]);
	double tau_factor            = REAL(s_tau_factor)[0];
	double radius_factor         = REAL(s_radius_factor)[0];
	size_t laplacian_power       = static_cast<size_t>(INTEGER(s_laplacian_power)[0]);
	size_t n_candidates          = static_cast<size_t>(INTEGER(s_n_candidates)[0]);
	bool log_grid                = LOGICAL(s_log_grid)[0];
	bool with_t_predictions      = LOGICAL(s_with_t_predictions)[0];
	bool verbose                 = LOGICAL(s_verbose)[0];

	// 2) Call the C++ smoother
	set_wgraph_t graph(adj_list, weight_list);
	graph.compute_graph_diameter();

	klaps_spectral_filter_t result = graph.klaps_spectral_filter(
		y,
		filter_type,
		n_evectors_to_compute,
		tau_factor,
		radius_factor,
		laplacian_power,
		n_candidates,
		log_grid,
		with_t_predictions,
		verbose
		);

	// 3) Build the R list to return
	int nprot = 0;
	int n_fields = with_t_predictions ? 7 : 6;
	SEXP Rresult = PROTECT(allocVector(VECSXP, n_fields)); nprot++;

	// 3a) evalues
	SEXP Revalues = PROTECT(allocVector(REALSXP, result.evalues.size())); nprot++;
	for (int i = 0; i < (int)result.evalues.size(); ++i)
		REAL(Revalues)[i] = result.evalues[i];
	SET_VECTOR_ELT(Rresult, 0, Revalues);

	// 3b) evectors
	int nrow = result.evectors.rows();
	int ncol = result.evectors.cols();
	SEXP Revectors = PROTECT(allocMatrix(REALSXP, nrow, ncol)); nprot++;
	for (int c = 0; c < ncol; ++c)
		for (int r = 0; r < nrow; ++r)
			REAL(Revectors)[r + nrow * c] = result.evectors(r, c);
	SET_VECTOR_ELT(Rresult, 1, Revectors);

	// 3c) candidate_ts
	SEXP Rts = PROTECT(allocVector(REALSXP, result.candidate_ts.size())); nprot++;
	for (int i = 0; i < (int)result.candidate_ts.size(); ++i)
		REAL(Rts)[i] = result.candidate_ts[i];
	SET_VECTOR_ELT(Rresult, 2, Rts);

	// 3d) gcv_scores
	SEXP Rgcv = PROTECT(allocVector(REALSXP, result.gcv_scores.size())); nprot++;
	for (int i = 0; i < (int)result.gcv_scores.size(); ++i)
		REAL(Rgcv)[i] = result.gcv_scores[i];
	SET_VECTOR_ELT(Rresult, 3, Rgcv);

	// 3e) opt_t_idx
	SEXP Ropt = PROTECT(allocVector(INTSXP, 1)); nprot++;
	INTEGER(Ropt)[0] = static_cast<int>(result.opt_t_idx + 1);   // 1‐based for R
	SET_VECTOR_ELT(Rresult, 4, Ropt);

	// 3f) predictions
	SEXP Rpred = PROTECT(allocVector(REALSXP, result.predictions.size())); nprot++;
	for (int i = 0; i < (int)result.predictions.size(); ++i)
		REAL(Rpred)[i] = result.predictions[i];
	SET_VECTOR_ELT(Rresult, 5, Rpred);

	// 3g) t_predictions (optional)
	if (with_t_predictions) {
		SEXP Rtp = PROTECT(allocMatrix(
							   REALSXP,
							   result.t_predictions[0].size(),
							   result.t_predictions.size()
							   )); nprot++;
		int nvert = result.t_predictions[0].size();
		for (int j = 0; j < (int)result.t_predictions.size(); ++j) {
			for (int i = 0; i < nvert; ++i) {
				REAL(Rtp)[i + nvert * j] = result.t_predictions[j][i];
			}
		}
		SET_VECTOR_ELT(Rresult, 6, Rtp);
	}

	// 4) Set names
	SEXP names = PROTECT(allocVector(STRSXP, n_fields)); nprot++;
	SET_STRING_ELT(names, 0, mkChar("evalues"));
	SET_STRING_ELT(names, 1, mkChar("evectors"));
	SET_STRING_ELT(names, 2, mkChar("candidate_ts"));
	SET_STRING_ELT(names, 3, mkChar("gcv_scores"));
	SET_STRING_ELT(names, 4, mkChar("opt_t_idx"));
	SET_STRING_ELT(names, 5, mkChar("predictions"));
	if (with_t_predictions)
		SET_STRING_ELT(names, 6, mkChar("t_predictions"));
	setAttrib(Rresult, R_NamesSymbol, names);

	UNPROTECT(nprot);
	return Rresult;
}
