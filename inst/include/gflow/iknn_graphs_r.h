#ifndef IKNN_GRAPHS_R_H_
#define IKNN_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_verify_pruning(SEXP s_X,
						  SEXP s_k,
						  SEXP s_max_alt_path_length);

	SEXP S_create_single_iknn_graph(
		SEXP s_X,
		SEXP s_k,
		SEXP s_pruning_thld,
		SEXP s_compute_full
		);

	SEXP S_create_iknn_graphs(
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
#endif // IKNN_GRAPHS_R_H_
