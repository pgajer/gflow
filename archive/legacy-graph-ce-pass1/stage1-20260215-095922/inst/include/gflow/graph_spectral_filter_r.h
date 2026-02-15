#ifndef GRAPH_SPECTRAL_FILTER_R_H_
#define GRAPH_SPECTRAL_FILTER_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_graph_spectral_filter(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_laplacian_type,
		SEXP s_filter_type,
		SEXP s_laplacian_power,
		SEXP s_kernel_tau_factor,
		SEXP s_kernel_radius_factor,
		SEXP s_kernel_type,
		SEXP s_kernel_adaptive,
		SEXP s_min_radius_factor,
		SEXP s_max_radius_factor,
		SEXP s_domain_min_size,
		SEXP s_precision,
		SEXP s_n_evectors_to_compute,
		SEXP s_n_candidates,
		SEXP s_log_grid,
		SEXP s_with_t_predictions,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_SPECTRAL_FILTER_R_H_
