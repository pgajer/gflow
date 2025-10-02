#ifndef RIEM_DCX_R_H_
#define RIEM_DCX_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_build_nerve_from_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	SEXP S_riem_dcx_summary(SEXP);
	SEXP S_get_simplices(SEXP s_dcx_ptr, SEXP s_dim);
	SEXP S_get_metric_diagonal(SEXP s_dcx_ptr, SEXP s_dim);

#ifdef __cplusplus
}
#endif

#endif // RIEM_DCX_R_H_
