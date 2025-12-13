/*!
 * @file magelo_external_bb.c
 * @brief Bayesian bootstrap functions that accept external (pre-generated) weights
 *
 * These functions enable paired BB comparisons where the same weights are used
 * for both signal and null computations. This is essential for proper paired
 * hypothesis testing in fassoc1.test() and similar functions.
 *
 * @author Pawel Gajer
 * @date 2025
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Forward declarations of existing functions from magelo.c */
extern void C_llm_1D_beta(const double *Tnn_x,
                          const double *Tnn_y,
                          double *Tnn_w,
                          const int *maxK,
                          const int *rnrTnn,
                          const int *rncTnn,
                          const int *rdeg,
                          double *beta);

extern void C_wpredict_1D(const double *beta,
                          const int *Tnn_i,
                          const double *Tnn_w,
                          const double *Tnn_x,
                          const int *maxK,
                          const int *rnrTnn,
                          const int *rncTnn,
                          const int *rdeg,
                          const int *rnx,
                          const int *rybinary,
                          double *Ey);

extern void C_llm_1D_fit_and_predict(const int *Tnn_i,
                                     double *Tnn_w,
                                     const double *Tnn_x,
                                     const double *Tnn_y,
                                     const int *rybinary,
                                     const int *maxK,
                                     const int *nrTnn,
                                     const int *rncTnn,
                                     const int *rnx,
                                     const int *rdeg,
                                     double *Ey,
                                     double *beta);

#ifndef CHECK_PTR
#define CHECK_PTR(ptr) if ((ptr) == NULL) { Rf_error("Memory allocation failed"); }
#endif


/*!
 * @brief Creates BB predictions using externally provided weights (global reweighting)
 *
 * This function is identical to C_llm_1D_fit_and_predict_global_BB except that
 * it accepts pre-generated Dirichlet weights instead of generating them internally.
 * This enables paired comparisons where signal and null use the same weights.
 *
 * @param Tnn_i      Matrix of indices of K nearest neighbors (nrTnn x ncTnn)
 * @param Tnn_w      Matrix of base weights (nrTnn x ncTnn). Weights must sum to 1.
 * @param Tnn_x      Matrix of x values over neighbors (nrTnn x ncTnn)
 * @param Tnn_y      Matrix of y values over neighbors (nrTnn x ncTnn)
 * @param rybinary   Pointer to binary indicator (1 = restrict Ey to [0,1])
 * @param maxK       Array of indices where weights are non-zero
 * @param rnrTnn     Pointer to number of rows of Tnn_ matrices
 * @param rncTnn     Pointer to number of columns (grid size)
 * @param rnx        Pointer to number of data points
 * @param rdeg       Pointer to polynomial degree (1 or 2)
 * @param rnBB       Pointer to number of BB iterations
 * @param lambda     External weights matrix (nx x nBB), column-major
 * @param gbbEy      Output: BB estimates of Ey (nx x nBB)
 */
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
    double       *gbbEy)
{
    int nrTnn = rnrTnn[0];
    int ncTnn = rncTnn[0];
    int nTnn  = ncTnn * nrTnn;
    int nx    = rnx[0];
    int nBB   = rnBB[0];

    int deg   = rdeg[0];
    int ncoef = deg + 1;

    /* Allocate working arrays */
    double *beta = (double*)calloc(ncoef * ncTnn, sizeof(double));
    CHECK_PTR(beta);

    double *Tnn_BB_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_BB_w);

    double *bbEy = (double*)malloc(nx * sizeof(double));
    CHECK_PTR(bbEy);

    int inx;
    for (int iboot = 0; iboot < nBB; iboot++)
    {
        inx = iboot * nx;

        /* Use external weights from column iboot of lambda matrix */
        const double *lambda_b = lambda + iboot * nx;

        /* Apply lambda weights to neighbor structure */
        for (int j = 0; j < nTnn; j++)
            Tnn_BB_w[j] = lambda_b[Tnn_i[j]] * Tnn_w[j];

        /* Fit and predict */
        C_llm_1D_fit_and_predict(Tnn_i, Tnn_BB_w, Tnn_x, Tnn_y, rybinary,
                                 maxK, rnrTnn, rncTnn, rnx, rdeg, bbEy, beta);

        /* Store results */
        for (int i = 0; i < nx; i++)
            gbbEy[i + inx] = bbEy[i];
    }

    free(beta);
    free(Tnn_BB_w);
    free(bbEy);
}


/*!
 * @brief Creates BB gpredictions and derivatives using external weights
 *
 * This function is identical to C_get_BB_Eyg except that it accepts
 * pre-generated Dirichlet weights instead of generating them internally.
 *
 * @param rn_BB        Pointer to number of BB iterations
 * @param Tnn_i        Matrix of neighbor indices
 * @param Tnn_x        Matrix of x values over neighbors
 * @param Tnn_y        Matrix of y values over neighbors
 * @param rybinary     Pointer to binary indicator
 * @param Tnn_w        Matrix of base weights
 * @param rnx          Pointer to number of data points
 * @param rnrTnn       Pointer to number of rows of Tnn_ matrices
 * @param rncTnn       Pointer to number of columns (grid size)
 * @param max_K        Array of max neighbor indices
 * @param rdegree      Pointer to polynomial degree
 * @param Tgrid_nn_i   Grid neighbor indices
 * @param Tgrid_nn_x   Grid neighbor x values
 * @param Tgrid_nn_w   Grid neighbor weights
 * @param rnrTgrid_nn  Pointer to rows of grid neighbor matrices
 * @param rncTgrid_nn  Pointer to cols of grid neighbor matrices
 * @param grid_max_K   Grid max neighbor indices
 * @param lambda       External weights matrix (nx x n_BB), column-major
 * @param bb_dEyg      Output: BB derivative estimates (ncTnn x n_BB)
 * @param bb_Eyg       Output: BB gprediction estimates (ncTnn x n_BB)
 */
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
    double       *bb_Eyg)
{
    int n_BB   = rn_BB[0];
    int nx     = rnx[0];
    int ncTnn  = rncTnn[0];
    int nrTnn  = rnrTnn[0];
    int nTnn   = ncTnn * nrTnn;
    int deg    = rdegree[0];
    int ncoef  = deg + 1;

    /* Allocate working arrays */
    double *bb_Tnn_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(bb_Tnn_w);

    double *bb_beta = (double*)malloc(ncoef * ncTnn * sizeof(double));
    CHECK_PTR(bb_beta);

    for (int i_BB = 0; i_BB < n_BB; i_BB++)
    {
        /* Use external weights from column i_BB of lambda matrix */
        const double *lambda_b = lambda + i_BB * nx;

        /* Apply lambda weights to neighbor structure */
        for (int j = 0; j < nTnn; j++)
            bb_Tnn_w[j] = lambda_b[Tnn_i[j]] * Tnn_w[j];

        /* Fit local linear models */
        C_llm_1D_beta(Tnn_x, Tnn_y, bb_Tnn_w, max_K, rnrTnn, rncTnn, rdegree, bb_beta);

        /* Predict on grid */
        C_wpredict_1D(bb_beta, Tgrid_nn_i, Tgrid_nn_w, Tgrid_nn_x, grid_max_K,
                      rnrTgrid_nn, rncTgrid_nn, rdegree, rncTnn, rybinary,
                      bb_Eyg + i_BB * ncTnn);

        /* Extract derivatives (slope coefficients) */
        for (int j = 0; j < ncTnn; j++)
            bb_dEyg[j + i_BB * ncTnn] = bb_beta[1 + j * ncoef];
    }

    free(bb_Tnn_w);
    free(bb_beta);
}


/*!
 * @brief Generates Dirichlet(1,...,1) weights (uniform on simplex)
 *
 * This is a utility function that generates random weights summing to n
 * (not normalized to 1, as BB uses n * Dirichlet weights).
 *
 * @param n      Number of weights to generate
 * @param n_BB   Number of weight vectors to generate
 * @param lambda Output matrix (n x n_BB), column-major
 */
void C_generate_dirichlet_weights(const int *rn,
                                  const int *rn_BB,
                                  double    *lambda)
{
    int n    = rn[0];
    int n_BB = rn_BB[0];

    for (int b = 0; b < n_BB; b++)
    {
        double *lambda_b = lambda + b * n;
        double sum = 0.0;

        /* Generate Exp(1) random variables */
        for (int i = 0; i < n; i++)
        {
            lambda_b[i] = -log(unif_rand());
            sum += lambda_b[i];
        }

        /* Normalize to sum to n (not 1, as BB uses n * w) */
        double scale = (double)n / sum;
        for (int i = 0; i < n; i++)
            lambda_b[i] *= scale;
    }
}


/* R-callable wrapper for C_generate_dirichlet_weights */
SEXP S_generate_dirichlet_weights(SEXP s_n, SEXP s_n_BB)
{
    int n    = INTEGER(s_n)[0];
    int n_BB = INTEGER(s_n_BB)[0];

    SEXP s_lambda = PROTECT(Rf_allocMatrix(REALSXP, n, n_BB));
    double *lambda = REAL(s_lambda);

    GetRNGstate();

    int rn    = n;
    int rn_BB = n_BB;
    C_generate_dirichlet_weights(&rn, &rn_BB, lambda);

    PutRNGstate();

    UNPROTECT(1);
    return s_lambda;
}


/* R-callable wrapper for C_llm_1D_fit_and_predict_global_BB_external */
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
    SEXP s_lambda)
{
    int nrTnn = INTEGER(s_nrTnn)[0];
    int ncTnn = INTEGER(s_ncTnn)[0];
    int nx    = INTEGER(s_nx)[0];
    int nBB   = INTEGER(s_nBB)[0];

    SEXP s_gbbEy = PROTECT(Rf_allocMatrix(REALSXP, nx, nBB));

    C_llm_1D_fit_and_predict_global_BB_external(
        INTEGER(s_Tnn_i),
        REAL(s_Tnn_w),
        REAL(s_Tnn_x),
        REAL(s_Tnn_y),
        INTEGER(s_ybinary),
        INTEGER(s_maxK),
        INTEGER(s_nrTnn),
        INTEGER(s_ncTnn),
        INTEGER(s_nx),
        INTEGER(s_deg),
        INTEGER(s_nBB),
        REAL(s_lambda),
        REAL(s_gbbEy));

    UNPROTECT(1);
    return s_gbbEy;
}


/* R-callable wrapper for C_get_BB_Eyg_external */
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
    SEXP s_lambda)
{
    int n_BB  = INTEGER(s_n_BB)[0];
    int ncTnn = INTEGER(s_ncTnn)[0];

    /* Allocate output matrices */
    SEXP s_bb_Eyg  = PROTECT(Rf_allocMatrix(REALSXP, ncTnn, n_BB));
    SEXP s_bb_dEyg = PROTECT(Rf_allocMatrix(REALSXP, ncTnn, n_BB));

    C_get_BB_Eyg_external(
        INTEGER(s_n_BB),
        INTEGER(s_Tnn_i),
        REAL(s_Tnn_x),
        REAL(s_Tnn_y),
        INTEGER(s_ybinary),
        REAL(s_Tnn_w),
        INTEGER(s_nx),
        INTEGER(s_nrTnn),
        INTEGER(s_ncTnn),
        INTEGER(s_max_K),
        INTEGER(s_degree),
        INTEGER(s_Tgrid_nn_i),
        REAL(s_Tgrid_nn_x),
        REAL(s_Tgrid_nn_w),
        INTEGER(s_nrTgrid_nn),
        INTEGER(s_ncTgrid_nn),
        INTEGER(s_grid_max_K),
        REAL(s_lambda),
        REAL(s_bb_dEyg),
        REAL(s_bb_Eyg));

    /* Create result list */
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP s_names  = PROTECT(Rf_allocVector(STRSXP, 2));

    SET_VECTOR_ELT(s_result, 0, s_bb_Eyg);
    SET_VECTOR_ELT(s_result, 1, s_bb_dEyg);

    SET_STRING_ELT(s_names, 0, Rf_mkChar("bb.gpredictions"));
    SET_STRING_ELT(s_names, 1, Rf_mkChar("bb.dgpredictions"));

    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    UNPROTECT(4);
    return s_result;
}
