/*!
 * @file magelo_external_bb.h
 * @brief Header for Bayesian bootstrap functions with external weights
 */

#ifndef MAGELO_EXTERNAL_BB_H
#define MAGELO_EXTERNAL_BB_H

#include <R.h>
#include <Rinternals.h>

/* C functions */
void C_llm_1D_fit_and_predict_global_BB_external(
    const int    *Tnn_i,
    const double *Tnn_w,
    const double *Tnn_x,
    const double *Tnn_y,
    const int    *rybinary,
    const int    *maxK,
    const int    *rnrTnn,
    const int    *rncTnn,
    const int    *rnx,
    const int    *rdeg,
    const int    *rnBB,
    const double *lambda,
    double       *gbbEy);

void C_get_BB_Eyg_external(
    const int    *rn_BB,
    const int    *Tnn_i,
    const double *Tnn_x,
    const double *Tnn_y,
    const int    *rybinary,
    const double *Tnn_w,
    const int    *rnx,
    const int    *rnrTnn,
    const int    *rncTnn,
    const int    *max_K,
    const int    *rdegree,
    const int    *Tgrid_nn_i,
    const double *Tgrid_nn_x,
    const double *Tgrid_nn_w,
    const int    *rnrTgrid_nn,
    const int    *rncTgrid_nn,
    const int    *grid_max_K,
    const double *lambda,
    double       *bb_dEyg,
    double       *bb_Eyg);

void C_generate_dirichlet_weights(const int *rn,
                                  const int *rn_BB,
                                  double    *lambda);

/* SEXP wrappers for R */
SEXP S_generate_dirichlet_weights(SEXP s_n, SEXP s_n_BB);

SEXP S_llm_1D_fit_and_predict_global_BB_external(
    SEXP s_Tnn_i,
    SEXP s_Tnn_w,
    SEXP s_Tnn_x,
    SEXP s_Tnn_y,
    SEXP s_ybinary,
    SEXP s_maxK,
    SEXP s_nrTnn,
    SEXP s_ncTnn,
    SEXP s_nx,
    SEXP s_deg,
    SEXP s_nBB,
    SEXP s_lambda);

SEXP S_get_BB_Eyg_external(
    SEXP s_n_BB,
    SEXP s_Tnn_i,
    SEXP s_Tnn_x,
    SEXP s_Tnn_y,
    SEXP s_ybinary,
    SEXP s_Tnn_w,
    SEXP s_nx,
    SEXP s_nrTnn,
    SEXP s_ncTnn,
    SEXP s_max_K,
    SEXP s_degree,
    SEXP s_Tgrid_nn_i,
    SEXP s_Tgrid_nn_x,
    SEXP s_Tgrid_nn_w,
    SEXP s_nrTgrid_nn,
    SEXP s_ncTgrid_nn,
    SEXP s_grid_max_K,
    SEXP s_lambda);

#endif /* MAGELO_EXTERNAL_BB_H */
