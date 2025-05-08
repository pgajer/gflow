#ifndef KLAPS_SPECTRAL_FILTER_R_H_
#define KLAPS_SPECTRAL_FILTER_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_klaps_spectral_filter(
	SEXP s_graph,
	SEXP s_edge_length,
	SEXP s_y,
	SEXP s_filter_type,
	SEXP s_n_evectors_to_compute,
	SEXP s_laplacian_power,
	SEXP s_tau_factor,
	SEXP s_n_candidates,
	SEXP s_log_grid,
	SEXP s_with_t_predictions,
	SEXP s_verbose
	);

#ifdef __cplusplus
}
#endif
#endif // KLAPS_SPECTRAL_FILTER_R_H_
