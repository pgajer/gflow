#ifndef METRIC_GRAPH_LOWPASS_R_H_
#define METRIC_GRAPH_LOWPASS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_metric_graph_lowpass_operator(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_conductance_rule,
    SEXP s_conductance_epsilon,
    SEXP s_conductance_alpha,
    SEXP s_conductance_sigma,
    SEXP s_sigma_rule,
    SEXP s_sigma_quantile,
    SEXP s_local_k,
    SEXP s_laplacian_type
);

SEXP S_metric_graph_lowpass_spectrum(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_conductance_rule,
    SEXP s_conductance_epsilon,
    SEXP s_conductance_alpha,
    SEXP s_conductance_sigma,
    SEXP s_sigma_rule,
    SEXP s_sigma_quantile,
    SEXP s_local_k,
    SEXP s_laplacian_type,
    SEXP s_n_eigenpairs,
    SEXP s_eigen_solver,
    SEXP s_dense_eigen_threshold,
    SEXP s_dense_fallback_threshold,
    SEXP s_dense_fallback,
    SEXP s_verbose
);

#ifdef __cplusplus
}
#endif

#endif
