#ifndef KLAPS_LOW_PASS_SMOOTHER_R_H_
#define KLAPS_LOW_PASS_SMOOTHER_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_klaps_low_pass_smoother(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_n_evectors_to_compute,
	SEXP s_min_num_eigenvectors,
	SEXP s_max_num_eigenvectors,
	SEXP s_tau_factor,
	SEXP s_radius_factor,
	SEXP s_laplacian_power,
	SEXP s_n_candidates,
	SEXP s_log_grid,
	SEXP s_energy_threshold,
	SEXP s_with_k_predictions,
	SEXP s_verbose
);

#ifdef __cplusplus
}
#endif
#endif // KLAPS_LOW_PASS_SMOOTHER_R_H_
