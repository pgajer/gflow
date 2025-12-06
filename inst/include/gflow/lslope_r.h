#ifndef LSLOPE_R_H
#define LSLOPE_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief SEXP interface for gradient-restricted local slope (instrumented)
 *
 * Computes asymmetric association measures along gradient direction with
 * full diagnostic output.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector of directing function values
 * @param s_z R numeric vector of response function values
 * @param s_type R character: "slope", "normalized", or "sign"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
 * @param s_sigmoid_alpha R numeric: sigmoid scale (0 = auto-calibrate)
 * @param s_ascending R logical: use ascending (TRUE) or descending (FALSE) gradient
 *
 * @return R list with coefficients and diagnostics
 */
SEXP S_lslope_gradient_instrumented(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_sigmoid_alpha,
    SEXP s_ascending
);

/**
 * @brief SEXP interface for gradient-restricted local slope (production)
 *
 * Streamlined version returning only coefficient vector.
 */
SEXP S_lslope_gradient(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_sigmoid_alpha,
    SEXP s_ascending
);

/**
 * @brief SEXP interface for neighborhood local regression coefficient
 *
 * Computes local regression coefficient using all neighborhood edges.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector of directing function values
 * @param s_z R numeric vector of response function values
 * @param s_weight_type R character: "unit" or "derivative"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
 *
 * @return R list with coefficients and diagnostics
 */
SEXP S_lslope_neighborhood(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_weight_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon
);

#ifdef __cplusplus
}
#endif

#endif // LSLOPE_R_H
