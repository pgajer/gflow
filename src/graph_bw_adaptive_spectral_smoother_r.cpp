#include "graph_bw_adaptive_spectral_smoother.hpp" // For graph_bw_adaptive_spectral_smoother_t
#include "error_utils.h"                  // For REPORT_ERROR
#include "set_wgraph.hpp"                 // For set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"  // For converting R objects to C++

#include <vector>               // For std::vector
#include <string>               // For std::string
#include <algorithm>            // For std::copy

#include <R.h>                  // For Rprintf, Rf_error, etc.
#include <Rinternals.h>         // For SEXP and R macros

extern "C" {
	SEXP S_graph_bw_adaptive_spectral_smoother(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_n_evectors,
		SEXP s_min_bw_factor,
		SEXP s_max_bw_factor,
		SEXP s_n_bws,
		SEXP s_log_grid,
		SEXP s_kernel_type,
		SEXP s_dist_normalization_factor,
		SEXP s_precision,
		SEXP s_use_global_bw_grid,
		SEXP s_with_bw_predictions,
		SEXP s_with_vertex_bw_errors,
		SEXP s_verbose
		);
}

/**
 * @brief R interface for graph_bw_adaptive_spectral_smoother
 *
 * @details Converts R inputs into C++ types, invokes
 * set_wgraph_t::graph_bw_adaptive_spectral_smoother, and returns
 * the fields in graph_bw_adaptive_spectral_smoother_t:
 *   - predictions:           Optimal-bandwidth predictions (length |V|)
 *   - bw_predictions:        Per-bandwidth predictions (|V| x n_bws matrix)
 *   - bw_mean_abs_errors:    Mean absolute Rf_error for each bandwidth
 *   - vertex_min_bws:        Per-vertex minimum bandwidths
 *   - opt_bw_idx:            1-based index of chosen bandwidth
 *
 * @param s_adj_list                R list of integer vectors for adjacency
 * @param s_weight_list             R list of numeric vectors for edge weights
 * @param s_y                       Numeric vector of responses
 * @param s_n_evectors              Integer: number of Laplacian eigenvectors
 * @param s_min_bw_factor           Numeric: min bandwidth factor
 * @param s_max_bw_factor           Numeric: max bandwidth factor
 * @param s_n_bws                   Integer: number of bandwidths
 * @param s_log_grid                Logical: log-spaced grid?
 * @param s_kernel_type             Integer: kernel type code
 * @param s_dist_normalization_factor Numeric: distance normalization factor
 * @param s_precision               Numeric: precision for bandwidth grid
 * @param s_verbose                 Logical: verbose output?
 *
 * @return SEXP R list with fields: "predictions", "bw_predictions", "bw_mean_abs_errors",
 *         "vertex_min_bws", "opt_bw_idx"
 */
SEXP S_graph_bw_adaptive_spectral_smoother(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_n_evectors,
	SEXP s_min_bw_factor,
	SEXP s_max_bw_factor,
	SEXP s_n_bws,
	SEXP s_log_grid,
	SEXP s_kernel_type,
	SEXP s_dist_normalization_factor,
	SEXP s_precision,
	SEXP s_use_global_bw_grid,
	SEXP s_with_bw_predictions,
	SEXP s_with_vertex_bw_errors,
	SEXP s_verbose
	) {

	std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
	std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

	std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

	size_t n_evectors = (size_t)INTEGER(s_n_evectors)[0];
	double min_bw_factor = REAL(s_min_bw_factor)[0];
	double max_bw_factor = REAL(s_max_bw_factor)[0];
	size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
	bool log_grid = (LOGICAL(s_log_grid)[0] == 1);
	size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
	double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
	double precision = REAL(s_precision)[0];
	bool use_global_bw_grid = (LOGICAL(s_use_global_bw_grid)[0] == 1);
	bool with_bw_predictions = (LOGICAL(s_with_bw_predictions)[0] == 1);
	bool with_vertex_bw_errors = (LOGICAL(s_with_vertex_bw_errors)[0] == 1);
	bool verbose = (LOGICAL(s_verbose)[0] == 1);

	set_wgraph_t graph(adj_list, weight_list);
	graph_bw_adaptive_spectral_smoother_t result = graph.graph_bw_adaptive_spectral_smoother(
		y,
		n_evectors,
		min_bw_factor,
		max_bw_factor,
		n_bws,
		log_grid,
		kernel_type,
		dist_normalization_factor,
		precision,
		use_global_bw_grid,
		with_bw_predictions,
		with_vertex_bw_errors,
		verbose
		);

	const char* names[] = {
		"predictions",
		"bw_predictions",
		"bw_mean_abs_errors",
		"vertex_min_bws",
		"opt_bw_idx",
		NULL
	};

	int n_elements = 5;
	SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
	{
		SEXP r_names  = PROTECT(Rf_allocVector(STRSXP, n_elements));
		for (int i = 0; i < n_elements; ++i) {
			SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
		}
		Rf_setAttrib(r_result, R_NamesSymbol, r_names);
		UNPROTECT(1);
	}

	// predictions (slot 0): numeric vector
	{
		const int n = (int) result.predictions.size();
		SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
		if (n > 0) {
			double* ptr = REAL(s);
			std::copy(result.predictions.begin(), result.predictions.end(), ptr);
		}
		SET_VECTOR_ELT(r_result, 0, s);
		UNPROTECT(1); // s
	}

	// bw_predictions (slot 1): numeric matrix [nrow x ncol]
	{
		const size_t ncol_sz = result.bw_predictions.size();
		const size_t nrow_sz = ncol_sz ? result.bw_predictions[0].size() : 0;

		// validate shapes before allocating
		for (size_t j = 0; j < ncol_sz; ++j) {
			if (result.bw_predictions[j].size() != nrow_sz)
				Rf_error("Inconsistent bw_predictions dimensions");
		}

		const int ncol = (int) ncol_sz;
		const int nrow = (int) nrow_sz;

		SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
		double* ptr = REAL(s);
		for (int j = 0; j < ncol; ++j) {
			const auto& col = result.bw_predictions[(size_t)j];
			for (int i = 0; i < nrow; ++i)
				ptr[i + (size_t)j * nrow] = col[(size_t)i];
		}
		SET_VECTOR_ELT(r_result, 1, s);
		UNPROTECT(1);  // s
	}

	// bw_mean_abs_errors (slot 2): numeric vector
	{
		const int n = (int) result.bw_mean_abs_errors.size();
		SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
		if (n > 0) {
			double* ptr = REAL(s);
			std::copy(result.bw_mean_abs_errors.begin(),
					  result.bw_mean_abs_errors.end(), ptr);
		}
		SET_VECTOR_ELT(r_result, 2, s);
		UNPROTECT(1);  // s
	}

	// vertex_min_bws (slot 3): numeric vector
	{
		const int n = (int) result.vertex_min_bws.size();
		SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
		if (n > 0) {
			double* ptr = REAL(s);
			std::copy(result.vertex_min_bws.begin(),
					  result.vertex_min_bws.end(), ptr);
		}
		SET_VECTOR_ELT(r_result, 3, s);
		UNPROTECT(1);  // s
	}

	// opt_bw_idx (slot 4): integer scalar, 1-based
	{
		SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
		INTEGER(s)[0] = (int) result.opt_bw_idx + 1;
		SET_VECTOR_ELT(r_result, 4, s);
		UNPROTECT(1);  // s
	}

	UNPROTECT(1);
	return r_result;
}
