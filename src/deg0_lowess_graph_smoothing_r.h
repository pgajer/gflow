#ifndef DEG0_LOWESS_GRAPH_SMOOTHING_R_H_
#define DEG0_LOWESS_GRAPH_SMOOTHING_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

    SEXP S_deg0_lowess_graph_smoothing(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_X,
        SEXP s_max_iterations,
        SEXP s_convergence_threshold,
        SEXP s_convergence_type,
        SEXP s_k,
        SEXP s_pruning_thld,
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_n_folds,
        SEXP s_use_uniform_weights,
		SEXP s_outlier_thld,
        SEXP s_with_bw_predictions,
        SEXP s_switch_to_residuals_after,
        SEXP s_verbose
    );

#ifdef __cplusplus
}
#endif
#endif // DEG0_LOWESS_GRAPH_SMOOTHING_R_H_
