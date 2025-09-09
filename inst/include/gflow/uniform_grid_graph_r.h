#ifndef UNIFORM_GRID_GRAPH_R_H_
#define UNIFORM_GRID_GRAPH_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_create_uniform_grid_graph(SEXP s_input_adj_list,
									 SEXP s_input_weight_list,
									 SEXP s_grid_size,
									 SEXP s_start_vertex,
									 SEXP s_snap_tolerance);
#ifdef __cplusplus
}
#endif
#endif // UNIFORM_GRID_GRAPH_R_H_
