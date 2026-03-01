#ifndef PATH_GRAPHS_R_H_
#define PATH_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_create_path_graph_series(SEXP s_adj_list,
                                SEXP s_weight_list,
                                SEXP s_h_values);

SEXP S_create_path_graph_plus(SEXP s_adj_list,
                              SEXP s_edge_length_list,
                              SEXP s_h);

SEXP S_create_path_graph_plm(SEXP s_adj_list,
                             SEXP s_edge_length_list,
                             SEXP s_h);

#ifdef __cplusplus
}
#endif

#endif // PATH_GRAPHS_R_H_
