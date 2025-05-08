#ifndef UGGMALOG_R_H_
#define UGGMALOG_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_uggmalog(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_best_models_coverage_factor,
		SEXP s_min_bw_factor,
		SEXP s_max_bw_factor,
		SEXP s_n_bws,
		SEXP s_grid_size,
		SEXP s_start_vertex,
		SEXP s_snap_tolerance,
		SEXP s_dist_normalization_factor,
		SEXP s_min_path_size,
		SEXP s_diff_threshold,
		SEXP s_kernel_type,
		SEXP s_fit_quadratic,
		SEXP s_max_iterations,
		SEXP s_ridge_lambda,
		SEXP s_tolerance,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // UGGMALOG_R_H_
