#ifndef GRAPH_DEG0_LOWESS_R_H_
#define GRAPH_DEG0_LOWESS_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

SEXP S_graph_deg0_lowess(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_bandwidth,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_verbose);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_DEG0_LOWESS_R_H_
