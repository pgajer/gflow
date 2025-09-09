#ifndef PGMALO_R_H_
#define PGMALO_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_pgmalo(
        SEXP neighbors_r,
        SEXP edge_lengths_r,
        SEXP y_r,
        SEXP y_true_r,
        SEXP use_median_r,
        SEXP h_min_r,
        SEXP h_max_r,
        SEXP p_r,
        SEXP n_bb_r,
        SEXP bb_max_distance_deviation_r,
        SEXP n_CVs_r,
        SEXP n_CV_folds_r,
        SEXP seed_r,
        SEXP kernel_type_r,
        SEXP dist_normalization_factor_r,
        SEXP epsilon_r,
        SEXP verbose_r);

	SEXP S_upgmalo(SEXP s_x,
				   SEXP s_y,
				   SEXP s_y_true,
				   SEXP s_use_median,
				   SEXP s_h_min,
				   SEXP s_h_max,
				   SEXP s_p,
				   SEXP s_n_bb,
				   SEXP s_bb_max_distance_deviation,
				   SEXP s_n_CVs,
				   SEXP s_n_CV_folds,
				   SEXP s_seed,
				   SEXP s_ikernel,
				   SEXP s_dist_normalization_factor,
				   SEXP s_epsilon,
				   SEXP s_verbose);

#ifdef __cplusplus
}
#endif
#endif // PGMALO_R_H_
