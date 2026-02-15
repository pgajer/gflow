#ifndef ADAPTIVE_UGGMALO_R_H_
#define ADAPTIVE_UGGMALO_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_adaptive_uggmalo(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_min_path_size,
		SEXP s_n_grid_vertices,
		SEXP s_n_bws,
		SEXP s_min_bw_factor,
		SEXP s_max_bw_factor,
		SEXP s_max_iterations,
		SEXP s_precision,
		SEXP s_dist_normalization_factor,
		SEXP s_kernel_type,
		SEXP s_tolerance,
		SEXP s_n_bb,
		SEXP s_cri_probability,
		SEXP s_n_perms,
		SEXP s_blending_coef,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // ADAPTIVE_UGGMALO_R_H_
