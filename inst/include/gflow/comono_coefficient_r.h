#ifndef COMONO_COEFFICIENT_R_H
#define COMONO_COEFFICIENT_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_comono(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_type
		);

	SEXP S_comono_matrix(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_Z,
		SEXP s_type
		);

#ifdef __cplusplus
}
#endif

#endif // COMONO_COEFFICIENT_R_H
