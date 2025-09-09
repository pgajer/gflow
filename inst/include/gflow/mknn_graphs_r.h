#ifndef MKNN_GRAPHS_R_H_
#define MKNN_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_create_mknn_graph(SEXP s_X, SEXP s_k);

	SEXP S_create_mknn_graphs(
		SEXP s_X,
		SEXP s_kmin,
		SEXP s_kmax,
		// pruning parameters
		SEXP s_max_path_edge_ratio_thld,
		SEXP s_path_edge_ratio_percentile,
		// other
		SEXP s_compute_full,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // MKNN_GRAPHS_R_H_
