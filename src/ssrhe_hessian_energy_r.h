#ifndef SSRHE_HESSIAN_ENERGY_R_H_
#define SSRHE_HESSIAN_ENERGY_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_ssrhe_hessian_operator(
    SEXP s_X,
    SEXP s_k,
    SEXP s_tangent_dim,
    SEXP s_nn_index,
    SEXP s_support_index,
    SEXP s_tangent_dim_rule,
    SEXP s_eigen_tolerance,
    SEXP s_derivative_order,
    SEXP s_stabilizer,
    SEXP s_pinv_tol,
    SEXP s_local_solver,
    SEXP s_normal_equations_max_condition,
    SEXP s_verbose
);

#ifdef __cplusplus
}
#endif

#endif
