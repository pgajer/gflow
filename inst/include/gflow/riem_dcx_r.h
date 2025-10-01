#ifndef RIEM_DCX_R_H_
#define RIEM_DCX_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_build_nerve_from_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
	SEXP S_riem_dcx_summary(SEXP);

#ifdef __cplusplus
}
#endif

#endif // RIEM_DCX_R_H_
