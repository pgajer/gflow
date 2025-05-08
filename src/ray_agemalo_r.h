#ifndef RAY_AGEMALO_R_H_
#define RAY_AGEMALO_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_ray_agemalo(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_min_path_size,
		// packing parameters
		SEXP s_n_packing_vertices,
		SEXP s_max_packing_iterations,
		SEXP s_packing_precision,
		// bw parameters
		SEXP s_n_bws,
		SEXP s_log_grid,
		SEXP s_min_bw_factor,
		SEXP s_max_bw_factor,
		// kernel parameters
		SEXP s_dist_normalization_factor,
		SEXP s_kernel_type,
		// model
		SEXP s_model_tolerance,
		SEXP s_blending_coef,
		// Bayesian bootstrap parameters
		SEXP s_n_bb,
		SEXP s_cri_probability,
		// permutation parameters
		SEXP s_n_perms,
		// verbose
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // RAY_AGEMALO_R_H_
