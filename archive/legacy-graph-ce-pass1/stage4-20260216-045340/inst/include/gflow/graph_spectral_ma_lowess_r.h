#ifndef GRAPH_SPECTRAL_MA_LOWESS_R_H_
#define GRAPH_SPECTRAL_MA_LOWESS_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_graph_spectral_ma_lowess(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_n_evectors,
        // bw parameters
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        // kernel parameters
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        // model parameters
        SEXP s_blending_coef,
        // other
        SEXP s_precision,
        SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_SPECTRAL_MA_LOWESS_R_H_
