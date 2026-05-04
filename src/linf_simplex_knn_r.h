#ifndef LINFINITY_SIMPLEX_KNN_R_H
#define LINFINITY_SIMPLEX_KNN_R_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_linf_simplex_knn(SEXP s_X, SEXP s_k, SEXP s_linf_tol);

#ifdef __cplusplus
}
#endif

#endif // LINFINITY_SIMPLEX_KNN_R_H
