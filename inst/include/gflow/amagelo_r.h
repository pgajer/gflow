#ifndef AMAGELO_R_H_
#define AMAGELO_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_amagelo(
		SEXP s_x,
		SEXP s_y,
		SEXP s_grid_size,
		SEXP s_min_bw_factor,
		SEXP s_max_bw_factor,
		SEXP s_n_bws,
		SEXP s_use_global_bw_grid,
		SEXP s_with_bw_predictions,
		SEXP s_log_grid,
		SEXP s_domain_min_size,
		SEXP s_kernel_type,
		SEXP s_dist_normalization_factor,
		SEXP s_n_cleveland_iterations,
		SEXP s_blending_coef,
		SEXP s_use_linear_blending,
		SEXP s_precision,
		SEXP s_small_depth_threshold,
		SEXP s_depth_similarity_tol,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // AMAGELO_R_H_
