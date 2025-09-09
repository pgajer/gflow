#include <R.h>            // For Rprintf, etc.
#include <Rinternals.h>   // For SEXP macros
#undef length             // avoid clash with R
#undef eval

#include <vector>
#include <algorithm>      // std::copy

#include "set_wgraph.hpp"                // set_wgraph_t
#include "klaps_low_pass_smoother.hpp"   // klaps_low_pass_smoother_t declaration
#include "error_utils.h"                 // REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp" // convert_adj_list_from_R, convert_weight_list_from_R

extern "C" {
	SEXP S_klaps_low_pass_smoother(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_n_evectors_to_compute,
		SEXP s_min_num_eigenvectors,
		SEXP s_max_num_eigenvectors,
		SEXP s_tau_factor,
		SEXP s_radius_factor,
		SEXP s_laplacian_power,
		SEXP s_n_candidates,
		SEXP s_log_grid,
		SEXP s_energy_threshold,
		SEXP s_with_k_predictions,
		SEXP s_verbose
		) {
		// Convert inputs
		std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
		std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

		const double* y_ptr = REAL(s_y);
		R_xlen_t n = XLENGTH(s_y);
		std::vector<double> y(y_ptr, y_ptr + n);

		size_t n_evectors_to_compute  = static_cast<size_t>(INTEGER(s_n_evectors_to_compute)[0]);
		size_t min_num_eigenvectors   = static_cast<size_t>(INTEGER(s_min_num_eigenvectors)[0]);
		size_t max_num_eigenvectors   = static_cast<size_t>(INTEGER(s_max_num_eigenvectors)[0]);
		double tau_factor             = REAL(s_tau_factor)[0];
		double radius_factor          = REAL(s_radius_factor)[0];
		size_t laplacian_power        = static_cast<size_t>(INTEGER(s_laplacian_power)[0]);
		size_t n_candidates           = static_cast<size_t>(INTEGER(s_n_candidates)[0]);
		bool   log_grid               = LOGICAL(s_log_grid)[0];
		double energy_threshold       = REAL(s_energy_threshold)[0];
		bool   with_k_predictions     = LOGICAL(s_with_k_predictions)[0];
		bool   verbose                = LOGICAL(s_verbose)[0];

		// Call the C++ smoother
		set_wgraph_t graph(adj_list, weight_list);
		graph.compute_graph_diameter();

		klaps_low_pass_smoother_t result =
			graph.klaps_low_pass_smoother(
				y,
				n_evectors_to_compute,
				min_num_eigenvectors,
				max_num_eigenvectors,
				tau_factor,
				radius_factor,
				laplacian_power,
				n_candidates,
				log_grid,
				energy_threshold,
				with_k_predictions,
				verbose
				);

		// Prepare output list
		const char* names[] = {
			"evalues",
			"evectors",
			"candidate_ks",
			"eigengaps",
			"opt_k_eigengap",
			"gcv_scores",
			"opt_k_gcv",
			"spectral_energy",
			"opt_k_spectral_energy",
			"used_method",
			"predictions",
			"k_predictions",
			NULL
		};
		int n_elements = 12;

		int protect_count = 0;
		SEXP r_result = PROTECT(allocVector(VECSXP, n_elements)); protect_count++;
		SEXP r_names  = PROTECT(allocVector(STRSXP, n_elements)); protect_count++;
		for (int i = 0; i < n_elements; ++i) {
			SET_STRING_ELT(r_names, i, mkChar(names[i]));
		}
		setAttrib(r_result, R_NamesSymbol, r_names);

		// Helpers for conversion
		auto vec_to_real = [&](const std::vector<double>& v) {
			SEXP x = PROTECT(allocVector(REALSXP, v.size())); protect_count++;
			std::copy(v.begin(), v.end(), REAL(x));
			return x;
		};
		auto mat_to_real = [&](const Eigen::MatrixXd& M) {
			SEXP x = PROTECT(allocMatrix(REALSXP, M.rows(), M.cols())); protect_count++;
			double* ptr = REAL(x);
			for (int j = 0; j < M.cols(); ++j)
				for (int i = 0; i < M.rows(); ++i)
					ptr[i + j * M.rows()] = M(i,j);
			return x;
		};
		auto vec_to_int = [&](const std::vector<size_t>& v, bool one_based=false) {
			SEXP x = PROTECT(allocVector(INTSXP, v.size())); protect_count++;
			for (R_xlen_t i = 0; i < v.size(); ++i)
				INTEGER(x)[i] = one_based ? (int)v[i] + 1 : (int)v[i];
			return x;
		};
		auto scalar_int = [&](size_t v, bool one_based=false) {
			SEXP x = PROTECT(allocVector(INTSXP, 1)); protect_count++;
			INTEGER(x)[0] = one_based ? (int)v + 1 : (int)v;
			return x;
		};

		// Fill in slots
		SET_VECTOR_ELT(r_result, 0, vec_to_real(result.evalues));             // evalues
		SET_VECTOR_ELT(r_result, 1, mat_to_real(result.evectors));            // evectors
		SET_VECTOR_ELT(r_result, 2, vec_to_int(result.candidate_ks, true));   // candidate_ks
		SET_VECTOR_ELT(r_result, 3, vec_to_real(result.eigengaps));           // eigengaps
		SET_VECTOR_ELT(r_result, 4, scalar_int(result.opt_k_eigengap, true)); // opt_k_eigengap
		SET_VECTOR_ELT(r_result, 5, vec_to_real(result.gcv_scores));          // gcv_scores
		SET_VECTOR_ELT(r_result, 6, scalar_int(result.opt_k_gcv, true));      // opt_k_gcv
		SET_VECTOR_ELT(r_result, 7, vec_to_real(result.spectral_energy));     // spectral_energy
		SET_VECTOR_ELT(r_result, 8, scalar_int(result.opt_k_spectral_energy, true)); // opt_k_spectral_energy

		// used_method: 0=Eigengap,1=GCV,2=EnergyThreshold
		SEXP Rmethod = PROTECT(allocVector(INTSXP, 1)); protect_count++;
		INTEGER(Rmethod)[0] = static_cast<int>(result.used_method);
		SET_VECTOR_ELT(r_result, 9, Rmethod);

		SET_VECTOR_ELT(r_result,10, vec_to_real(result.predictions));  // predictions

		// k_predictions: n_vertices × n_candidates matrix (possibly empty)
		{
			size_t ncol = result.k_predictions.size();
			size_t nrow = ncol ? result.k_predictions[0].size() : 0;
			SEXP Kpred = PROTECT(allocMatrix(REALSXP, nrow, ncol)); protect_count++;
			double* ptr = REAL(Kpred);
			for (size_t j = 0; j < ncol; ++j) {
				for (size_t i = 0; i < nrow; ++i) {
					ptr[i + j * nrow] = result.k_predictions[j][i];
				}
			}
			SET_VECTOR_ELT(r_result,11, Kpred);
		}

		UNPROTECT(protect_count);
		return r_result;
	}
}  // extern "C"
