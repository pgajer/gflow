#ifndef ADAPTIVE_MAELOG_R_H_
#define ADAPTIVE_MAELOG_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

   SEXP S_adaptive_maelog(
	   SEXP x_r,
	   SEXP y_r,
	   SEXP max_adapt_iterations_r,
	   SEXP convergence_threshold_r,
	   SEXP c_min_r,
	   SEXP c_max_r,
	   SEXP power_r,
	   SEXP kernel_type_r,
	   SEXP min_points_r,
	   SEXP max_iterations_r,
	   SEXP ridge_lambda_r,
	   SEXP tolerance_r
	   );

#ifdef __cplusplus
}
#endif
#endif // ADAPTIVE_MAELOG_R_H_
