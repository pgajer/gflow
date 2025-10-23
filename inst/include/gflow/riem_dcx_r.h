#ifndef RIEM_DCX_R_H_
#define RIEM_DCX_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_fit_rdgraph_regression(
		SEXP s_X,
		SEXP s_y,
		SEXP s_k,
		SEXP s_use_counting_measure,
		SEXP s_density_normalization,

		SEXP s_t_diffusion,
		SEXP s_beta_damping,
		SEXP s_gamma_modulation,
		SEXP s_t_scale_factor,
		SEXP s_beta_coefficient_factor,

		SEXP s_n_eigenpairs,
		SEXP s_filter_type,
		SEXP s_epsilon_y,
		SEXP s_epsilon_rho,
		SEXP s_max_iterations,

		SEXP s_max_ratio_threshold,
		SEXP s_threshold_percentile,
		SEXP s_density_alpha,
		SEXP s_density_epsilon,
		SEXP s_compute_extremality,

		SEXP s_p_threshold,
		SEXP s_max_hop,
		SEXP s_test_stage,
		SEXP s_verbose
		);

	SEXP S_compute_hop_extremp_radii_batch(
		SEXP s_adj_list,
		SEXP s_edge_densities,
		SEXP s_vertex_densities,
		SEXP s_candidates,
		SEXP s_y,
		SEXP s_p_threshold,
		SEXP s_detect_maxima,
		SEXP s_max_hop
		);

	SEXP S_compute_basins_of_attraction(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y
		);

#ifdef __cplusplus
}
#endif

#endif // RIEM_DCX_R_H_
