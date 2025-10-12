#ifndef RIEM_DCX_R_H_
#define RIEM_DCX_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_fit_knn_riem_graph_regression(
		SEXP s_X,
		SEXP s_y,
		SEXP s_k,
		SEXP s_use_counting_measure,
		SEXP s_density_normalization,
		SEXP s_t_diffusion,
		SEXP s_beta_damping,
		SEXP s_gamma_modulation,
		SEXP s_n_eigenpairs,
		SEXP s_filter_type,
		SEXP s_epsilon_y,
		SEXP s_epsilon_rho,
		SEXP s_max_iterations,
		SEXP s_max_ratio_threshold,
		SEXP s_threshold_percentile,
		SEXP s_test_stage,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif

#endif // RIEM_DCX_R_H_
