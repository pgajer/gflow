#ifndef GRAPH_SPECTRAL_LOWESS_MAT_R_H_
#define GRAPH_SPECTRAL_LOWESS_MAT_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_graph_spectral_lowess_mat(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_Y,
        SEXP s_n_evectors,
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_precision,
        SEXP s_n_cleveland_iterations,
        SEXP s_with_errors,
        SEXP s_with_scale,
        SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_SPECTRAL_LOWESS_MAT_R_H_
