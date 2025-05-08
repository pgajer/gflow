#ifndef KNN_R_H_
#define KNN_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

    SEXP S_kNN(SEXP RX, SEXP Rk);
    SEXP S_kNN_v2(SEXP RX, SEXP Rk);

#ifdef __cplusplus
}
#endif
#endif // KNN_R_H_
