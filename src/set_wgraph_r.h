#ifndef SET_WGRAPH_R_H_
#define SET_WGRAPH_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_compute_edge_weight_deviations(SEXP s_adj_list,
										  SEXP s_weight_list);

	SEXP S_compute_edge_weight_rel_deviations(SEXP s_adj_list,
											  SEXP s_weight_list);

	SEXP S_remove_redundant_edges(SEXP s_adj_list,
								  SEXP s_weight_list);

	SEXP S_find_graph_paths_within_radius(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_start,
		SEXP s_radius);

#ifdef __cplusplus
}
#endif
#endif // SET_WGRAPH_R_H_
