#ifndef GRAPH_DEG0_LOWESS_BUFFER_CV_R_H_
#define GRAPH_DEG0_LOWESS_BUFFER_CV_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_graph_deg0_lowess_buffer_cv(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
	SEXP s_use_uniform_weights,
    SEXP s_buffer_hops,        
    SEXP s_auto_buffer_hops,   
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_DEG0_LOWESS_BUFFER_CV_R_H_
