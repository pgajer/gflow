#ifndef MST_COMPLETION_GRAPH_R_H_
#define MST_COMPLETION_GRAPH_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_create_mst_completion_graph(
	SEXP s_X,
	SEXP s_q_thld,
	SEXP s_verbose
	);

#ifdef __cplusplus
}
#endif
#endif // MST_COMPLETION_GRAPH_R_H_
