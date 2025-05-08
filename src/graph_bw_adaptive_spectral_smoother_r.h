#ifndef GRAPH_BW_ADAPTIVE_SPECTRAL_SMOOTHER_R_H_
#define GRAPH_BW_ADAPTIVE_SPECTRAL_SMOOTHER_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

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
	SEXP s_with_global_bws,
	SEXP s_with_bw_predictions,
	SEXP s_with_vertex_bw_errors,
	SEXP s_verbose
	);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_BW_ADAPTIVE_SPECTRAL_SMOOTHER_R_H_
