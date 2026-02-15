#ifndef GRAPH_DIFFUSION_SMOOTHER_R_H_
#define GRAPH_DIFFUSION_SMOOTHER_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_graph_diffusion_smoother(SEXP s_adj_list,
									SEXP s_weight_list,
									SEXP s_y,
									SEXP s_n_time_steps,
									SEXP s_step_factor,
									SEXP s_binary_threshold,
									SEXP s_ikernel,
									SEXP s_dist_normalization_factor,
									SEXP s_n_CVs,
									SEXP s_n_CV_folds,
									SEXP s_epsilon,
									SEXP s_verbose,
									SEXP s_seed);

	SEXP S_ext_graph_diffusion_smoother(SEXP Rgraph,
										SEXP Rd,
										SEXP Rweights,
										SEXP Ry,
										SEXP Rn_time_steps,
										SEXP Rstep_factor,
										SEXP Rnormalize,
										SEXP Rpreserve_local_maxima,
										SEXP Rlocal_maximum_weight_factor,
										SEXP Rpreserve_local_extrema,
										SEXP Rimputation_method,
										SEXP Rmax_iterations,
										SEXP Rconvergence_threshold,
										SEXP Rapply_binary_threshold,
										SEXP Rbinary_threshold,
										SEXP Rikernel,
										SEXP Rdist_normalization_factor,
										SEXP Rn_CVs,
										SEXP Rn_CV_folds,
										SEXP Repsilon,
										SEXP Rseed,
										SEXP Rn_cores,
										SEXP Rverbose);

	SEXP S_instrumented_gds(SEXP s_graph,
							SEXP s_edge_lengths,
							SEXP s_y,
							SEXP s_y_true,
							SEXP s_n_time_steps,
							SEXP s_base_step_factor,
							SEXP s_use_pure_laplacian,
							SEXP s_ikernel,
							SEXP s_kernel_scale,
							SEXP s_increase_factor,
							SEXP s_decrease_factor,
							SEXP s_oscillation_factor,
							SEXP s_min_step,
							SEXP s_max_step);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_DIFFUSION_SMOOTHER_R_H_
