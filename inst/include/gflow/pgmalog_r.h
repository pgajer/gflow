#ifndef PGMALOG_R_H_
#define PGMALOG_R_H_
#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_upgmalog(SEXP s_x,
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
				   SEXP s_n_cores,
				   SEXP s_dist_normalization_factor,
				   SEXP s_epsilon,
				   SEXP s_verbose);

#ifdef __cplusplus
}
#endif
#endif // PGMALOG_R_H_
