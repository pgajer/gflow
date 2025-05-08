#ifndef MAGELOG_R_H_
#define MAGELOG_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

SEXP S_magelog(
    SEXP x_r,
    SEXP y_r,
    SEXP grid_size_r,  // Add this parameter
    SEXP fit_quadratic_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP min_points_r,
    SEXP cv_folds_r,
    SEXP n_bws_r,
    SEXP min_bw_factor_r,
    SEXP max_bw_factor_r,
    SEXP max_iterations_r,
    SEXP ridge_lambda_r,
    SEXP tolerance_r,
    SEXP with_bw_predictions_r);

#ifdef __cplusplus
}
#endif
#endif // MAGELOG_R_H_
