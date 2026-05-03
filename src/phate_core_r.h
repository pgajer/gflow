#ifndef GFLOW_PHATE_CORE_R_H_
#define GFLOW_PHATE_CORE_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_phate_build_kernel(SEXP s_D,
                          SEXP s_k,
                          SEXP s_decay,
                          SEXP s_thresh,
                          SEXP s_bandwidth_scale,
                          SEXP s_symm_code,
                          SEXP s_diag_one);

#ifdef __cplusplus
}
#endif

#endif // GFLOW_PHATE_CORE_R_H_
