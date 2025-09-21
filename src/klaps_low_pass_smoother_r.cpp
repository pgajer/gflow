#include "set_wgraph.hpp"                // set_wgraph_t
#include "klaps_low_pass_smoother.hpp"   // klaps_low_pass_smoother_t declaration
#include "error_utils.h"                 // REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp" // convert_adj_list_from_R, convert_weight_list_from_R

#include <vector>
#include <algorithm>      // std::copy

#include <R.h>            // For Rprintf, etc.
#include <Rinternals.h>   // For SEXP macros

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
		bool   log_grid               = (LOGICAL(s_log_grid)[0] == 1);
		double energy_threshold       = REAL(s_energy_threshold)[0];
		bool   with_k_predictions     = (LOGICAL(s_with_k_predictions)[0] == 1);
		bool   verbose                = (LOGICAL(s_verbose)[0] == 1);

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
		SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
		{
			SEXP r_names  = PROTECT(Rf_allocVector(STRSXP, n_elements));
			for (int i = 0; i < n_elements; ++i) {
				SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
			}
			Rf_setAttrib(r_result, R_NamesSymbol, r_names);
			UNPROTECT(1);
		}

		// evalues (slot 0): REALSXP vector
		{
			const int n = (int) result.evalues.size();
			SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
			if (n > 0) {
				double* p = REAL(s);
				std::copy(result.evalues.begin(), result.evalues.end(), p);
			}
			SET_VECTOR_ELT(r_result, 0, s);
			UNPROTECT(1); // s
		}

		// evectors (slot 1): REALSXP matrix (rows=M.rows, cols=M.cols)
		{
			const int nrow = (int) result.evectors.rows();
			const int ncol = (int) result.evectors.cols();
			SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
			double* ptr = REAL(s);
			for (int j = 0; j < ncol; ++j)
				for (int i = 0; i < nrow; ++i)
					ptr[i + (size_t)j * nrow] = result.evectors(i, j);
			SET_VECTOR_ELT(r_result, 1, s);
			UNPROTECT(1);  // s
		}

		// candidate_ks (slot 2): INTSXP vector (1-based)
		{
			const int n = (int) result.candidate_ks.size();
			SEXP s = PROTECT(Rf_allocVector(INTSXP, n));
			if (n > 0) {
				int* ip = INTEGER(s);
				for (int i = 0; i < n; ++i) ip[i] = (int) result.candidate_ks[(size_t)i] + 1;
			}
			SET_VECTOR_ELT(r_result, 2, s);
			UNPROTECT(1);  // s
		}

		// eigengaps (slot 3): REALSXP vector
		{
			const int n = (int) result.eigengaps.size();
			SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
			if (n > 0) {
				double* p = REAL(s);
				std::copy(result.eigengaps.begin(), result.eigengaps.end(), p);
			}
			SET_VECTOR_ELT(r_result, 3, s);
			UNPROTECT(1);  // s
		}

		// opt_k_eigengap (slot 4): INTSXP scalar (1-based)
		{
			SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
			INTEGER(s)[0] = (int) result.opt_k_eigengap + 1;
			SET_VECTOR_ELT(r_result, 4, s);
			UNPROTECT(1);  // s
		}

		// gcv_scores (slot 5): REALSXP vector
		{
			const int n = (int) result.gcv_scores.size();
			SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
			if (n > 0) {
				double* p = REAL(s);
				std::copy(result.gcv_scores.begin(), result.gcv_scores.end(), p);
			}
			SET_VECTOR_ELT(r_result, 5, s);
			UNPROTECT(1);  // s
		}

		// opt_k_gcv (slot 6): INTSXP scalar (1-based)
		{
			SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
			INTEGER(s)[0] = (int) result.opt_k_gcv + 1;
			SET_VECTOR_ELT(r_result, 6, s);
			UNPROTECT(1);  // s
		}

		// spectral_energy (slot 7): REALSXP vector
		{
			const int n = (int) result.spectral_energy.size();
			SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
			if (n > 0) {
				double* p = REAL(s);
				std::copy(result.spectral_energy.begin(), result.spectral_energy.end(), p);
			}
			SET_VECTOR_ELT(r_result, 7, s);
			UNPROTECT(1);  // s
		}

		// opt_k_spectral_energy (slot 8): INTSXP scalar (1-based)
		{
			SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
			INTEGER(s)[0] = (int) result.opt_k_spectral_energy + 1;
			SET_VECTOR_ELT(r_result, 8, s);
			UNPROTECT(1);  // s
		}

		// used_method (slot 9): INTSXP scalar (0=Eigengap, 1=GCV, 2=EnergyThreshold)
		{
			SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
			INTEGER(s)[0] = (int) result.used_method;
			SET_VECTOR_ELT(r_result, 9, s);
			UNPROTECT(1);  // s
		}

		// predictions (slot 10): REALSXP vector
		{
			const int n = (int) result.predictions.size();
			SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
			if (n > 0) {
				double* p = REAL(s);
				std::copy(result.predictions.begin(), result.predictions.end(), p);
			}
			SET_VECTOR_ELT(r_result, 10, s);
			UNPROTECT(1);  // s
		}

		// k_predictions (slot 11): REALSXP matrix [nrow x ncol]
		{
			const size_t ncol_sz = result.k_predictions.size();
			const size_t nrow_sz = ncol_sz ? result.k_predictions[0].size() : 0;
			const int ncol = (int) ncol_sz;
			const int nrow = (int) nrow_sz;

			// (optional) shape check
			for (size_t j = 0; j < ncol_sz; ++j) {
				if (result.k_predictions[j].size() != nrow_sz)
					Rf_error("Inconsistent k_predictions dimensions");
			}

			SEXP Kpred = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
			double* ptr = REAL(Kpred);
			for (int j = 0; j < ncol; ++j)
				for (int i = 0; i < nrow; ++i)
					ptr[i + (size_t)j * nrow] = result.k_predictions[(size_t)j][(size_t)i];

			SET_VECTOR_ELT(r_result, 11, Kpred);
			UNPROTECT(1);  // Kpred
		}

		UNPROTECT(1);
		return r_result;
	}
}   // extern "C"
