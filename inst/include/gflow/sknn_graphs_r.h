#ifndef GFLOW_SKNN_GRAPHS_R_H
#define GFLOW_SKNN_GRAPHS_R_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_create_sknn_graph(SEXP s_X,
                         SEXP s_k,
                         SEXP s_connect_components,
                         SEXP s_connect_method);

#ifdef __cplusplus
}
#endif

#endif
