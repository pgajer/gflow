/*
  A collection of local linear regression (llm) and robust local linear regression (rllm) routines
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>

#include <stdlib.h>

#include "msr2.h"
#include "lm.h"
#include "kernels.h"  // for initialize_kernel()
#include "sampling.h" // for C_rsimplex() and C_runif_simplex()
#include "stats_utils.h" // for get_folds(), C_columnwise_eval() and iarray_t

static void (*kernel_with_stop_fn)(const double*, const int*, const double*, int*, double*);

/* get_folds() is now defined in stats_utils.c and declared in stats_utils.h */

/*
  \brief Gets bandwidths in the columns of d = t(nn.d).

  In each column, i, of the distances to K NN's, d = t(nn.d), there are least minK[i] elements with the distance to the NN's less than bw

  That is, bws[i] is defined as follows

  bws[i] = bw if d[minK[i]-1, i] = d[minK[i] - 1 + ir] < bw

  OR

  bws[i] = d[minK[i], i], if minK[i] < nr. Otherwise

  bws[i] = d[minK[i]-1, i] + eps


  \param d      A numeric 1D array corresponding to a nn.d matrix in R.
  \param rnr    A reference to the number of rows of d.
  \param rnc    A reference to the number of columns of d.
  \param minK   An array such that minK[i] is the required mininum __number__ of elements in the row of the i-th column of d.
  \param rbw    A reference to the value of a bandwidth parameter.

  \param bws    An output bandwidths array indicating the support of a given linear model.

*/
void C_get_bws_with_minK_a(const double *d,
                           const int    *rnr,
                           const int    *rnc,
                           const int    *minK,
                           const double *rbw,
                                 double *bws)
{
    int nr    = rnr[0]; // number of rows of d
    int nr1   = nr - 1; // index of the last element in each row of d
    int nc    = rnc[0]; // number of columns of d
    double bw = rbw[0];

    int iminK;         // holds minK[i] - 1 = the index of the minK[i]-th element in the i-th column of d
    double eps = 0.01; // a small constant added to z = d[iminK + ir] to make
                       // bws[i] = z + eps when iminK=nr1 or iminK < nr1 but
                       // z >= d[minK + ir]
    int ir;            // holds i * nr; corresponds to the index of the i-th column of d

    for ( int i = 0; i < nc; i++ )
    {
      ir = i * nr;
      iminK = minK[i] - 1;

      if ( iminK < 0 )
        Rf_error("ERROR in get_bws_with_minK_a(); file %s at line %d;  iminK=%d < 0!", __FILE__, __LINE__, iminK);

      if ( iminK >= nr )
        Rf_error("ERROR in get_bws_with_minK_a(); file %s at line %d;  iminK = (minK[%d] - 1) = %d >= nr=%d !", __FILE__, __LINE__, i, minK[i]-1, nr);

      if ( d[iminK + ir] < bw )
      {
        bws[i] = bw;

      } else if ( (iminK < nr1) && (d[iminK + 1 + ir] > d[iminK + ir]) ) { // it is not the last element and the next one is greater than that current one

        bws[i] = d[iminK + 1 + ir];

      } else {

        bws[i] = d[iminK + ir] + eps;
      }
    }
}

/*
  \brief Gets bandwidths

  Finds a bandwidth for each column of d, such that bws[i] = bw if d[iminK + ir]
  < bw and bws[i] = d[iminK + 1 + ir] = d[minK + ir], otherwise.

  \param d      A numeric 1D array corresponding to a nn.d matrix in R.
  \param rnr    A reference to the number of rows of the matrix associated with x.
  \param rnc    A reference to the number of columns of the matrix associated with x.
  \param rminK  A reference to the mininum __number_ of elements in the rows of d with weights > 0.
  \param rbw    A reference to a non-negative normalization constant.

  \param bws    An output bandwidths array indicating the support of a given linear model.

*/
void C_get_bws(const double *d,
             const int    *rnr,
             const int    *rnc,
             const int    *rminK,
             const double *rbw,
                   double *bws)
{
   int nr    = rnr[0];
   int nr1   = nr - 1;
   int nc    = rnc[0];
   int minK  = rminK[0];
   int iminK = minK - 1; // the index of the minK-th element
   double bw = rbw[0];

   double eps = 0.01; // a small constant added to z = d[iminK + ir] to make
                      // bws[i] = z + eps when iminK=nr1 or iminK < nr1 but
                      // z >= d[minK + ir]
   int ir;  // holds i * nr
   double z;
   for ( int i = 0; i < nc; i++ )
   {
     ir = i * nr;
     z = d[iminK + ir];

     if ( z < bw )
     {
       bws[i] = bw;

     } else if ( iminK < nr1 && z < d[minK + ir] ) {

       bws[i] = d[minK + ir]; // This makes only sense if d has at least (minK +
                              // 1) rows. Thus, nr >= minK + 1. The violation of
                              // this condition is nr < minK + 1 <=> nr - 1 <
                              // minK. In the case nr = minK, I could have
                              // defined bws[i] = d[iminK + ir] + 1, but I don't
                              // want to put here the 'if' clause.
     } else {

       bws[i] = z + eps;
     }
   }
}

/*!
    \fn void C_columnwise_weighting(double* x, int *rnr, int *rnc, int *rikernel, double* w)

    \brief Apply a kernel to the clumns of x, where x is a 1D array associated with an R matrix t(nn.d).

    This version of column-wise wieghting also populates maxK array, with
    maxK[i] = the index in x[,i] of the last non-zero element. maxK and w are
    output variables, both need to have memory allocated to them by the user.

    \param x         A 1D array of associated with distances to NN's (dimension nr-by-nc)
    \param rnr       A reference to the number of rows of nnDist
    \param rnc       A reference to the number of rows of nnDist (the number of nearest neighbors)
    \param rikernel  A reference to kernel index

    \param bws       An array of bandwidths.

    \param maxK      An output array of length nc with maxK[i] being the __index__ of
                     the last weight >0 element in the i-th column of x (=nn_w[,i]). If all
                     weights are 0 in the i-th column of x, then maxK[i] = -1.

    \param w         An output array of weights of the same size as x.

*/
void C_columnwise_weighting(const double *x,
                            const int    *rnr,
                            const int    *rnc,
                            const int    *rikernel,
                            const double *bws,
                                  int    *maxK,
                                  double *w) {

    int nr      = rnr[0];
    int nc      = rnc[0];
    int ikernel = rikernel[0];
    int nr1     = nr - 1;

    // Initializing maxK, setting all values to nr1.
    // This is done so that if a column of x has only non-zero weights, then
    // maxK is not updated and the value of nr1 is returned, which is the index
    // of the last element of each column of x.
    for ( int j = 0; j < nc; j++ )
      maxK[j] = nr1;

    // All elements of the array w need to be set to 0, in case maxK is changed
    // outside of this routine and elements that suppose to be 0 are 0.
    int n = nr * nc;
    for ( int i = 0; i < n; i++ )
      w[i] = 0;

    switch( ikernel )
    {
      case EPANECHNIKOV:
        kernel_with_stop_fn = C_epanechnikov_kernel_with_stop;
        break;
      case TRIANGULAR:
        kernel_with_stop_fn = C_triangular_kernel_with_stop;
        break;
      case TREXPONENTIAL:
        kernel_with_stop_fn = C_tr_exponential_kernel_with_stop;
        break;
      default:
        Rf_error("Unknown kernel in columnwise_weighting(): file %s line %d", __FILE__, __LINE__);
    }

    int jr; // holds j*nr = the index of the j-th column of x
    for (int j = 0; j < nc; j++)
    {
      jr = j*nr;
      kernel_with_stop_fn(x + jr, rnr, bws + j, maxK + j, w + jr);
    }
}


/*!
    \brief The main loop of the 1D local linear model that fits locally linear models. This version utilizes maxK parameter.

    \param Tnn_x      A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
    \param Tnn_y      A matrix of y values over K nearest neighbors of each element of the grid.
    \param Tnn_w      A matrix of weights over K nearest neighbors of each element of the grid.
    \param ma xK      An array of indices indicating the range where weights are not 0. Indices =< maxK[i] have weights > 0.
    \param rnrTnn     A reference to the number of rows of the above three matrices.
    \param rncTnn     A reference to the number of columns of the above three matrices.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.

    \param beta       An output array of the coefficients of the models.

    NOTE: Tnn_w is TS normalized within this function. Thus, weights do not have to sum up to 1.

*/
void C_llm_1D_beta(const double *Tnn_x,
                   const double *Tnn_y,
                         double *Tnn_w,
                   const int    *maxK,
                   const int    *rnrTnn,
                   const int    *rncTnn,
                   const int    *rdeg,
                         double *beta) {

    int ncTnn    = rncTnn[0];  // grid size = number of grid elements (break
                               // points); ncTnn-1 = number of grid intervals =
                               // numer of columns of Tnn_x, Tnn_y, Tnn_w
    int nrTnn    = rnrTnn[0];  // number of nearest neighbors = number of rows of Tnn_x ( as Tnn_x = t(nn.x) in R), Tnn_y, Tnn_w
    int deg      = rdeg[0];    // degree of x in the linear models
    int ncoef    = deg + 1;    // number of columns of beta
    int status   = 0;

    // Tnn_x, Tnn_y, Tnn_w are 1D arrays corresonding to matrices in R that
    // are passed to C in a row-major form (as they are transposed first).

    // beta was also transposed before passing it to this routine and in R it
    // was a matrix of dim n-by-(deg+1).

    // Normalizing Tnn_w so that the sum along each column is 1
    int G, jK;
    double s;
    for ( int j = 0; j < ncTnn; j++ )
    {
      G = maxK[j] + 1;
      jK = j * nrTnn;
      s = 0;
      for ( int k = 0; k < G; k++ )
        s += Tnn_w[k + jK];

      for ( int k = 0; k < G; k++ )
        Tnn_w[k + jK] /= s;
    }

    if ( deg == 1 )
    {
      int iK; // hold i * nrTnn
      int G;  // (index of the last non-zero element) + 1
      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        // Tnn_x + i*nrTnn      is the beginning of the array of x values over nrTnn NN's of the given grid element
        // Tnn_y + i*nrTnn      is the beginning of the array of y values over nrTnn NN's of the given grid element
        // Tnn_w + i*nrTnn is the beginning of the array of weights over nrTnn NN's of the given grid element
        // beta + i*ncoef   is the beginning of the array coefficients of the corresponding linear model
        iK = i*nrTnn;
        G = maxK[i] + 1;
        C_flmw_1D(Tnn_x + iK, Tnn_y + iK, Tnn_w + iK, &G, beta + i*ncoef, &status);

      } // end of for i

    } else { // deg > 1 case

      int nvar = 2; // number of columns of the new array
      // an array holding Tnn_x + i*nrTnn array of size nrTnn and it's squared values
      double *nnX = (double*)malloc(nvar * nrTnn * sizeof(double));
      CHECK_PTR( nnX );

      int iK;
      int G;  // (index of the last non-zero element) + 1
      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        iK = i*nrTnn;
        G = maxK[i] + 1;
        for ( int j = 0; j < G; j++ )
        {
          nnX[j] = Tnn_x[j + iK];
          nnX[j+G] = nnX[j]*nnX[j];
        }
        C_flmw(nnX, Tnn_y + iK, Tnn_w + iK, &G, &nvar, beta + i*ncoef, &status);
      }

      free(nnX);
    }
}

/*!
    \fn void C_columnwise_TS_norm(double *x, int *rnr, int *rnc, double *nx)

    \brief This routine divides each clumn of x by the sum of that colums
    elements, if the sum is not 0. If the sum is 0, nothing is done to the
    column. The result is updating x to it's normalized form.

    \param x         An array with nr rows and nc columns.
    \param rnr       A reference to the number of rows of x.
    \param rnc       A reference to the number of columns of x.

    \param nx        An output array normalized rows of x.

*/
void C_columnwise_TS_norm(const double *x,
                          const int    *rnr,
                          const int    *rnc,
                                double *nx)
{
    int nr = rnr[0];
    int nc = rnc[0];

    int ir;   // holds i * nr
    double s; // sum of column elements
    for ( int i = 0; i < nc; i++ )
    {
      s = 0;
      ir = i * nr;
      for ( int j = 0; j < nr; j++ )
      {
        s += x[j + ir];
        nx[j + ir] = x[j + ir];
      }

      if ( s != 0 )
      {
        for ( int j = 0; j < nr; j++ )
          nx[j + ir] /= s;
      }
    }
}

/* C_columnwise_eval() is now defined in stats_utils.c and declared in stats_utils.h */

/*!
    \brief Nearest neighbor weighted mean with maxK parameter.

    \param nn_i       nn_i[i,] indices of nc nearest neighbors of the i-th sample/point
    \param nn_Ey      nn_Ey[i,] mean values of y over nc nearest neighbors of the i-th sample/point
    \param nn_w       nn_w[i,] weights of y over nc nearest neighbors of the i-th sample/point
    \param maxK       An array of length ng with maxK[i] being the __index__ of the last non-zero weight in the i-th column nn_w[,i] of nn_w.
    \param rK         A reference of the number of rows of each of the above tables
    \param rng        A reference of the number of columns of each of the above tables
    \param rnx        A reference of the length of x (points over which Ey is estimated).

    \param Ey         An output variable of the nearest neighbor weighted means of y.

*/
void C_nn_wmean_maxK(const int    *nn_i,
                     const double *nn_Ey,
                     const double *nn_w,
                     const int    *maxK,
                     const int    *rK,  // number of rows of nn_'s
                     const int    *rng, // number of columns of nn_'s
                     const int    *rnx,
                           double *Ey) {
    int K  = rK[0];
    int ng = rng[0];
    int nx = rnx[0]; // nn_i[i]'s are in the range 0 .. (nx-1), as they index values of x, over which Ey is estimated

    // initializing Ey
    for ( int i = 0; i < nx; i++ )
      Ey[i] = 0;

    double* Ew = calloc(nx, sizeof(double));
    CHECK_PTR(Ew);

    int iK; // holds i * K
    int G;  // (index of the last non-zero element) + 1

    for ( int i = 0; i < ng; i++ ) {
      iK = i * K;
      G = maxK[i] + 1;
      for ( int j = 0; j < G; j++ )
      {
        Ey[nn_i[j + iK]] += nn_Ey[j + iK] * nn_w[j + iK];
        Ew[nn_i[j + iK]] += nn_w[j + iK];
      }
    }

    for ( int i = 0; i < nx; i++ )
    {
      if ( Ew[i] > 0.0 ){
        Ey[i] /= Ew[i];
      } else {
        //Ey[i] = R_NaN; // In prediction case we expect weights of some elements to be all 0 and Ey should be 0 in this case;
        Ey[i] = 0;
      }
    }

    free(Ew);
}


/*!
    \brief The main loop of the 1D local linear model that fits locally linear models. This version utilizes maxK parameter.

    \param Tnn_x      A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
    \param Tnn_y      A matrix of y values over K nearest neighbors of each element of the grid.
    \param Tnn_w      A matrix of weights over K nearest neighbors of each element of the grid.
    \param ma xK      An array of indices indicating the range where weights are not 0. Indices =< maxK[i] have weights > 0.
    \param rnrTnn     A reference to the number of rows of the above three matrices.
    \param rncTnn     A reference to the number of columns of the above three matrices.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.

    \param beta       An output array of the coefficients of the models.

    NOTE: Tnn_w is TS normalized within this function. Thus, weights do not have to sum up to 1.

*/
void C_llm_1D_beta_perms(const double *Tnn_x,
                         const int    *Tnn_i,
                         const double *y,
                         const int    *rny,
                         double       *Tnn_w,
                         const int    *maxK,
                         const int    *rnrTnn,
                         const int    *rncTnn,
                         const int    *rdeg,
                         const int    *rn_perms,
                         double *beta_perms)
{
    int n_perms = rn_perms[0];
    int ny      = rny[0];
    int ncTnn   = rncTnn[0];   // grid size = number of columns of Tnn_i, Tnn_d, Tnn_x and the number of elements of maxK
    int nrTnn   = rnrTnn[0];   // number of nearest neighbors = number of columns of Tnn_ matrices

    double *perm_y = (double *)malloc(ny * sizeof(double));
    CHECK_PTR(perm_y);

    // Copy y content to y_perm.
    memcpy(perm_y, y, ny * sizeof(double));

    double *Tnn_y = (double *)malloc(ncTnn * nrTnn * sizeof(double));
    CHECK_PTR(Tnn_y);

    for ( int i = 0; i < n_perms; i++ )
    {
      // in-place permutation of y
      // C_dpermute(perm_y, ny);  for testing purposes

      // row Rf_eval premuted y with nn.i
      C_columnwise_eval(Tnn_i,
                        rnrTnn,
                        rncTnn,
                        perm_y,
                        Tnn_y);

      // beta points to the i-th matrix of coefficients in the 3d array of permuted y coefficients
      double *beta = beta_perms + i * nrTnn * ncTnn;

      // call
      C_llm_1D_beta(Tnn_x,
                    Tnn_y,
                    Tnn_w,
                    maxK,
                    rnrTnn,
                    rncTnn,
                    rdeg,
                    beta);
    }

    free(perm_y);
    free(Tnn_y);
}


/*!

  Predicts the mean y over a vector of points x.

  \param beta       Coefficients of the local linear model.
  \param Tnn_i       Grid indices of NN's.
  \param Tnn_d       Distances to NN's within the grid.
  \param Tnn_x       The values of the predictor variable over indices of NN's.
  \param rnrTnn         A reference to the number of rows of Tnn_ matrices.
  \param rncTnn        A reference to the number of columns of the above four matrices, which is the same as the number of grid points.
  \param rdeg       A reference to the degree of x in the linear models.
  \param rikernel   The integeer index of a kernel used for generating weights.
  \param rnx        A reference of the length of x (points over which Ey is estimated).
  \param Ey         The predicted values of the mean of y over the grid points of x.

*/
void C_predict_1D(const double *beta,
                  const int    *Tnn_i,
                  const double *Tnn_d,
                  const double *Tnn_x,
                  const int    *maxK,
                  const int    *rnrTnn,
                  const int    *rncTnn,
                  const int    *rdeg,
                  const int    *rikernel,
                  const int    *rnx,
                        double *Ey)
{
    int ncTnn   = rncTnn[0];   // grid size = number of columns of Tnn_i, Tnn_d, Tnn_x and the number of elements of maxK
    int nrTnn   = rnrTnn[0];   // number of nearest neighbors = number of columns of Tnn_ matrices
    int deg     = rdeg[0];     // degree of x in the linear models
    int ncoef   = deg + 1;     // number of columns of beta
    int ikernel = rikernel[0]; // kernel index

    initialize_kernel(ikernel, 1.0);

    double *Tnn_w = (double *)malloc(nrTnn * ncTnn * sizeof(double));
    CHECK_PTR(Tnn_w);

    int iK;
    for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
    {
      iK = i * nrTnn;
      kernel_fn(Tnn_d + iK, nrTnn, Tnn_w + iK);
    }

    //
    // Mean y values over grid NN's
    //
    double *Tnn_Ey = (double *)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_Ey);

    int incoef;
    int G; // (index of the last non-zero element) + 1

    if ( deg == 1 )
    {
      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        iK = i * nrTnn;
        incoef = i * ncoef;
        G = maxK[i] + 1;
        for ( int j = 0; j < G; j++ )
          Tnn_Ey[j + iK] = beta[incoef] + Tnn_x[j + iK]* beta[1 + incoef];

      }
    } else {

      double x;
      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        iK = i * nrTnn;
        incoef = i * ncoef;
        G = maxK[i] + 1;
        for ( int j = 0; j < G; j++ )
        {
          x = Tnn_x[j + iK];
          Tnn_Ey[j + iK] = beta[incoef] + x*beta[1 + incoef] + x*x*beta[2 + incoef];
        }
      }
    }

    //
    // Computing Ey using Tnn_wmean()
    //
    //nn_wmean(Tnn_i, Tnn_Ey, Tnn_w, rnrTnn, rncTnn, rnx, Ey);
    C_nn_wmean_maxK(Tnn_i, Tnn_Ey, Tnn_w, maxK, rnrTnn, rncTnn, rnx, Ey);

    free(Tnn_Ey);
    free(Tnn_w);
}


/*!

  Predicts the mean y over NN's of each grid point using the coefficients of the
  corresponding weighted linear model. The predicted mean values form a matrix
  nn.Eyg that is then used to predict the weighted means of nn.Eyg at each grid
  point.

  \param beta       Coefficients of the local linear model.
  \param nn_i       Grid indices of NN's.
  \param nn_w       Weights to NN's within the grid.
  \param nn_x       The values of the predictor variable over indices of NN's.
  \param maxK       A vector of indices indicating the range where weights are not 0.
  \param rnrTnn     A reference to the number of rows of nn_ arrays.
  \param rncTnn     A reference to the number of columns of the above four matrices, which is the same as the number of grid points.
  \param rdeg       A reference to the degree of x in the linear models.
  \param rnx        A reference of the length of x (points over which Ey is estimated).
  \param rybinary   A reference to the integer indicator with 0 meaning y is not a binary variable.
  \param Ey         An output array of predicted values of the mean of y over x.

*/
void C_wpredict_1D(const double *beta,
                   const int    *Tnn_i,
                   const double *Tnn_w,
                   const double *Tnn_x,
                   const int    *maxK,
                   const int    *rnrTnn,
                   const int    *rncTnn,
                   const int    *rdeg,
                   const int    *rnx,
                   const int    *rybinary,
                         double *Ey) {

    int ncTnn   = rncTnn[0];  // grid size
    int nrTnn   = rnrTnn[0];  // number of nearest neighbors = number of columns of Tnn_ matrices
    int deg     = rdeg[0];    // degree of x in the linear models
    int ncoef   = deg + 1;    // number of columns of beta
    int ybinary = rybinary[0];

    // Predicted Ey over Tnn_x
    double *Tnn_Ey = (double *)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_Ey);

    int incoef, iK;
    int G; // (index of the last non-zero element) + 1

    if ( deg == 1 )
    {
      if ( ybinary == 0 )
      {
        for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
        {
          iK = i * nrTnn;
          incoef = i * ncoef;
          G = maxK[i] + 1;

          for ( int j = 0; j < G; j++ )
            Tnn_Ey[j + iK] = beta[incoef] + Tnn_x[j + iK]* beta[1 + incoef];
        }
      } else {

        double y;
        for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
        {
          iK = i * nrTnn;
          incoef = i * ncoef;
          G = maxK[i] + 1;

          for ( int j = 0; j < G; j++ )
          {
            y = beta[incoef] + Tnn_x[j + iK]* beta[1 + incoef];
            if ( y > 1 )
              y = 1;
            else if ( y < 0 )
              y = 0;
            Tnn_Ey[j + iK] = y;
          }
        }
      }

    } else if ( deg == 0 ) {

      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        iK = i * nrTnn;
        incoef = i * ncoef;
        G = maxK[i] + 1;

        for ( int j = 0; j < G; j++ )
          Tnn_Ey[j + iK] = beta[incoef];
      }

    } else {

      double x;
      for ( int i = 0; i < ncTnn; i++ )// iterating over elements of the grid
      {
        iK = i * nrTnn;
        incoef = i * ncoef;
        G = maxK[i] + 1;

        for ( int j = 0; j < G; j++ )
        {
          x = Tnn_x[j + iK];
          Tnn_Ey[j + iK] = beta[incoef] + x*beta[1 + incoef] + x*x*beta[2 + incoef];
        }
      }
    }

    // computing weighted means of Tnn_Ey with weights Tnn_w over indices 0, ... , (nx-1)
    C_nn_wmean_maxK(Tnn_i, Tnn_Ey, Tnn_w, maxK, rnrTnn, rncTnn, rnx, Ey);

    free(Tnn_Ey);
}

/*!

  Generates LOO llm_1D predictions of the mean y over x.

  \param x              1D array of predictor values.
  \param rnx            A reference to the number of elements of x.
  \param Tnn_i           1D array corresponding to a nrTnn-by-ncTnn matrix of x NN's indices to xg.
  \param Tnn_w           1D array corresponding to a nrTnn-by-ncTnn matrix of x NN's weights.
  \param Tnn_x           1D array corresponding to a nrTnn-by-ncTnn matrix of x values over NN's of xg.
  \param Tnn_y           1D array corresponding to a nrTnn-by-ncTnn matrix of y values over NN's of xg.
  \param rnrTnn             A reference to the number of rows of Tnn_ arrays.
  \param rncTnn            A reference to the number of columns of Tnn_ arrays, which is the same as the number of grid points.
  \param rdeg           A reference to the degree of x in the linear model.
  \param rmax_n_models  A reference to the maximal number of models assocaited with a single point.

  \param Ey             An output array of predicted values of the mean of y over x using models trained on LOO data.

*/
void C_loo_llm_1D(const double *x,
                  const int    *rnx,
                  const int    *Tnn_i,
                        double *Tnn_w,
                  const double *Tnn_x,
                  const double *Tnn_y,
                  const int    *rnrTnn,
                  const int    *rncTnn,
                  const int    *rdeg,
                  const int    *rmax_n_models,
                        double *Ey) {

    int nx           = rnx[0];           // number of elements of x
    int nrTnn        = rnrTnn[0];        // number of nearest neighbors = number of rows of Tnn_i, Tnn_w, Tnn_x, Tnn_y
    int ncTnn        = rncTnn[0];        // number of grid elements = number of columns of Tnn_i, Tnn_w, Tnn_x, Tnn_y
    int deg          = rdeg[0];          // degree of x in linear models
    int max_n_models = rmax_n_models[0]; // the maximal number of models assocaited with a single point.

    // Identifying grid (model) indices of models that cotain x[i] in their domains
    double *x_lm = (double*)calloc( max_n_models * nx, sizeof(double) ); // x_lm[,i] - xg indices of models associated with the i-th element of x
    CHECK_PTR(x_lm);

    int *nx_lm = (int*)calloc( nx, sizeof(int) ); // nx_lm[i] - The number of linear models associated with the i-th element of x
    CHECK_PTR(nx_lm);

    int *Tnn_ix = (int*)malloc( ncTnn * nx * sizeof(int) ); // Tnn_ix - index of x in Tnn_i
    CHECK_PTR(Tnn_ix);

    int ix;
    int iK;  // ig * nrTnn
    for ( int ig = 0; ig < ncTnn; ig++ )
    {
      iK = ig * nrTnn;
      for ( int j = 0; j < nrTnn; j++ )
      {
        ix = Tnn_i[j + iK]; // index of x in the i-th column and j-th row of Tnn_i = the j-th NN of the i-th xg
        Tnn_ix[ ig + ix*ncTnn ] = j;

        if ( Tnn_w[j + iK] > 0 )
        {
          x_lm[ nx_lm[ix] + ix*max_n_models] = ig;
          nx_lm[ix]++;
        }

      } // end of for j
    } // end of for ig


    // For each x[i] building linear models specified by x_lm[,i] using all
    // Tnn_x[,j], Tnn_y[,j], Tnn_w[,j] except x[i], y[i] and w[i], by setting w[i] to 0

    // single model Ey and w values
    double *sm_Ey = (double*)malloc( max_n_models * sizeof(double) );
    CHECK_PTR(sm_Ey);

    double *sm_w = (double*)malloc( max_n_models * sizeof(double) );
    CHECK_PTR(sm_w);

    // beta1 and beta2 carry coefficients of the linear models
    int ncoef = deg + 1;
    double *beta1 = (double*)malloc( ncoef * sizeof(double) );
    CHECK_PTR(beta1);

    double *beta2 = (double*)malloc( ncoef * sizeof(double) );
    CHECK_PTR(beta2);

    // After setting to 0 in xg[j] the weight of the i-th element of x, the
    // weights need to be rescaled so that they sum up to 1. This array holds
    // rescaled values.
    double *w = (double*)malloc( nrTnn * sizeof(double) );
    CHECK_PTR(w);

    int status = 0;
    int ixg;            // model index
    double cw;          // complement weight

    if ( deg == 1 )
    {
      for ( int i = 0; i < nx; i++ )
      {
        for ( int j = 0; j < nx_lm[i]; j++ )
        {
          ixg = x_lm[j + i*max_n_models];
          iK  = ixg*nrTnn;
          ix  = Tnn_ix[ ixg + i*ncTnn ]; // index of i in Tnn_i's ixg column - it is the index at which w needs to be set to 0

          sm_w[j] = *(Tnn_w + ix + iK); // storing the current value of *(Tnn_w + ix + iK), so we can reset it to its original value when we are done with the given i
          *(Tnn_w + ix + iK) = 0;

          cw = 1.0 - sm_w[j]; // this is what the sum of w's suppose to be after wetting *(Tnn_w + ix + iK) to
          for ( int k = 0; k < nrTnn; k++ ) {
            w[k] = *(Tnn_w + k + iK) / cw;
          }

          C_flmw_1D(Tnn_x + iK, Tnn_y + iK, w, &nrTnn, beta1, &status); // ideally we would want to replace nrTnn with k <= nrTnn so that w[k]>0 and w[k+1]=0
          sm_Ey[j] = beta1[0] + beta1[1] * x[i]; // prediction of Ey given the above LOO model

          // Resetting the value of *(Tnn_w + ix + iK) to what it was before we set it to 0
          *(Tnn_w + ix + iK) = sm_w[j] ;

        } // end of for j

        Ey[i] = wmean(sm_Ey, sm_w, nx_lm[i]);

      } // end of for i

    } else {

      int K2 = 2 * nrTnn;
      double *nnX;
      int nvar = 2;

      // an array holding x^2
      nnX = (double*)malloc( nvar * nrTnn * sizeof(double) );
      CHECK_PTR( nnX );

      for ( int i = 0; i < nx; i++ )
      {
        for ( int j = 0; j < nx_lm[i]; j++ )
        {
          ixg = x_lm[j + i*max_n_models];
          iK = ixg*nrTnn;
          ix = Tnn_ix[ ixg + i*ncTnn ]; // index of i in Tnn_i's ixg column - it is the index at which w needs to be set to 0

          sm_w[j] = *(Tnn_w + ix + iK);
          *(Tnn_w + ix + iK) = 0;

          cw = 1.0 - sm_w[j]; // this is what the sum of w's suppose to be after wetting *(Tnn_w + ix + iK) to
          for ( int k = 0; k < nrTnn; k++ ) {
            w[k] = *(Tnn_w + k + iK) / cw;
          }

          // copying the first nrTnn elements of Tnn_x + j + iK to nnX
          for ( int j = 0; j < nrTnn; j++ )
            nnX[j] = *(Tnn_x + j + iK);

          // Setting the remaining nrTnn elmenets of nnX to the square of the fist nrTnn
          for ( int j = nrTnn; j < K2; j++ )
            nnX[j] = nnX[j-nrTnn]*nnX[j-nrTnn];

          C_flmw(nnX, Tnn_y + iK, w, &nrTnn, &nvar, beta2, &status); // ideally we would want to replace nrTnn with k <= nrTnn so that w[k]>0 and w[k+1]=0

          sm_Ey[j] = beta2[0] + beta2[1]*x[i] + beta2[2]*x[i]*x[i];

          // Resetting the value of *(Tnn_w + ix + iK) to what it was before we set it to 0
          *(Tnn_w + ix + iK) = sm_w[j] ;

        } // end of for j

        Ey[i] = wmean(sm_Ey, sm_w, nx_lm[i]);

      } // end of for i

      free(nnX);
    }

    free(sm_Ey);
    free(sm_w);

    free(w);

    free(beta1);
    free(beta2);

    free(x_lm);
    free(nx_lm);
    free(Tnn_ix);
}


/*!

  Generates LOO llm_1D predictions of the mean y over x, based on degree 0 model.

  \param rnx            A reference to the number of elements of x.
  \param Tnn_i          A 1D array corresponding to a nrTnn-by-ncTnn matrix of x NN's indices to xg.
  \param Tnn_w          A 1D array corresponding to a nrTnn-by-ncTnn matrix of x NN's weights.
  \param Tnn_y          A 1D array corresponding to a nrTnn-by-ncTnn matrix of y values over NN's of xg.
  \param rnrTnn         A reference to the number of rows of Tnn_ arrays.
  \param rncTnn            A reference to the number of columns of Tnn_ arrays, which is the same as the number of grid points.
  \param Ey             The predicted values of the mean of y over x using models trained on LOO data.

*/
void C_deg0_loo_llm_1D(const int    *rnx,
                       const int    *Tnn_i,
                             double *Tnn_w,
                       const double *Tnn_y,
                       const int    *rnrTnn,
                       const int    *rncTnn,
                             double *Ey) {

    int nx    = rnx[0];     // number of elements of x
    int nrTnn = rnrTnn[0];  // number of nearest neighbors = number of rows of Tnn_i, Tnn_w, Tnn_x, Tnn_y
    int ncTnn = rncTnn[0];  // number of grid elements = number of columns of Tnn_i, Tnn_w, Tnn_x, Tnn_y
    int n     = ncTnn * nx;

    double *Eyg = (double*)malloc( ncTnn * sizeof(double) );
    CHECK_PTR(Eyg);

    int *Tnn_ix = (int*)malloc( n * sizeof(int) ); // Tnn_ix - is the inverse of Tnn_i: Tnn_ix[ ig + Tnn_i[j + ig*nrTnn]*ncTnn ] = j; in other words: Tnn_ix[ig, Tnn_i[j, jg]] = j
    CHECK_PTR(Tnn_ix);

    double *w = (double*)malloc( ncTnn * sizeof(double) );
    CHECK_PTR(w);

    // initializing Tnn_ix
    for ( int i = 0; i < n; i++ )
      Tnn_ix[i] = -1;

    int inrTnn;  // ig * nrTnn
    int ix;
    for ( int ig = 0; ig < ncTnn; ig++ )
    {
      inrTnn = ig * nrTnn;
      for ( int j = 0; j < nrTnn; j++ )
      {
        ix = Tnn_i[j + inrTnn]; // index of x in the i-th column and j-th row of Tnn_i = the j-th NN of the i-th xg
        Tnn_ix[ ig + ix*ncTnn ] = j;

      } // end of for j
    } // end of for i

    int j;
    for ( int ix = 0; ix < nx; ix++ )
    {
      for ( int ig = 0; ig < ncTnn; ig++ )
      {
        j = Tnn_ix[ ig + ix*ncTnn ]; // ix = Tnn_i[j + inrTnn]
        if ( j > -1 )
        {
          inrTnn  = ig * nrTnn;
          w[ig] = *(Tnn_w + j + inrTnn); // storing the current value of *(Tnn_w + j +
                                    // inrTnn), so we can reset it to its original
                                    // value when we are done with the given i
          *(Tnn_w + j + inrTnn) = 0;

        } else {

          w[ig] = 0;
        }

      } // end of for ig

      // Estimating Eyg using the data without the ix-th element of x
      //Tnn_wmean(Tnn_i, Tnn_y, Tnn_w, rnrTnn, rncTnn, rnx, Eyg);

      for ( int ig = 0; ig < ncTnn; ig++ )
      {
        inrTnn = ig * nrTnn;
        Eyg[ig] = wmean(Tnn_y + inrTnn, Tnn_w + inrTnn, nrTnn);
      }

      Ey[ix] = wmean(Eyg, w, ncTnn);

      // Resetting the value of *(Tnn_w + j + inrTnn) to what it was before we set it to 0
      for ( int ig = 0; ig < ncTnn; ig++ )
      {
        j = Tnn_ix[ ig + ix*ncTnn ]; // ix = Tnn_i[j + inrTnn]
        if ( j > -1 )
          *(Tnn_w + j + ig*nrTnn) = w[ig];
      }

    } // end of for ix


    free(Eyg);
    free(Tnn_ix);
    free(w);
}


/*!
    \brief A version of C_llm_1D_fit_and_predict() with nn.y replaced by a matrix with the number of columns equal to x used to construct nn.x.

    \param Y        A matrix with the number of columns that is the same as the length of x that was used to construct nn.* matrices. The main application of this routine is for the case when the rows of Y are permutations some y.
    \param nrY      A reference to the numer of rows of Y.
    \param ncY      A reference to the numer of columns of Y.
    \param Tnn_i     A matrix of indices of nrTnn nearest neighbors of each element of the grid, where nrTnn is determined in the parent routine.
    \param Tnn_w     A matrix of weights over nrTnn nearest neighbors of each element of the grid. Weights must sum up to 1!
    \param Tnn_x     A matrix of x values over nrTnn nearest neighbors of each element of the grid.
    \param nrTnn     A reference to the number of rows of the above Tnn_ matrices.
    \param ncTnn     A reference to the number of columns of the above Tnn_ matrices.
    \param maxK      An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
    \param deg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param EYg       An output array of predicted values of the mean of the rows of Y over X using local linear models.

*/
void C_mllm_1D_fit_and_predict(const double *Y,
                               const int    *nrY,
                               const int    *ncY,
                               const int    *rybinary,

                               const int    *Tnn_i,
                                     double *Tnn_w,
                               const double *Tnn_x,
                               const int    *nrTnn,
                               const int    *ncTnn,

                               const int    *maxK,
                               const int    *deg,
                                     double *EYg)
{
    // Allocating memory for beta
    double *beta = (double*)malloc( (*deg + 1) * (*ncTnn) * sizeof(double) );
    CHECK_PTR(beta);

    double *Tnn_y = (double*)malloc( (*nrTnn) * (*ncTnn) * sizeof(double) );
    CHECK_PTR(Tnn_y);

    int j, k;
    for ( int i = 0; i < *ncY; i++ )
    {
      j = i * (*nrY);
      k = i * (*ncTnn);

      C_columnwise_eval(Tnn_i, nrTnn, ncTnn, Y + j, Tnn_y);
      C_llm_1D_beta(Tnn_x, Tnn_y, Tnn_w, maxK, nrTnn, ncTnn, deg, beta);
      C_wpredict_1D(beta, Tnn_i, Tnn_w, Tnn_x, maxK, nrTnn, ncTnn, deg, ncY, rybinary, EYg + k);
    }

    free(Tnn_y);
    free(beta);
}

/*!
    \brief Combines llm_1D_beta() and wpredict_1D()

    \param Tnn_i      A matrix of indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
    \param Tnn_w      A matrix of weights over K nearest neighbors of each element of the grid. Weights must sum up to 1!
    \param Tnn_x      A matrix of x values over K nearest neighbors of each element of the grid.
    \param Tnn_y      A matrix of y values over K nearest neighbors of each element of the grid.
    \param maxK       An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
    \param nrTnn      A reference to the number of rows of the above Tnn_ matrices.
    \param rncTnn        A reference to the number of columns of the above Tnn_ matrices.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param rnx        A reference to the length of x.

    \param Ey         An output array of predicted values of the mean of y over X using X.grid local linear models.
    \param beta       An output array of the coefficients of the models.

*/
void C_llm_1D_fit_and_predict(const int    *Tnn_i,
                                    double *Tnn_w,
                              const double *Tnn_x,
                              const double *Tnn_y,
                              const int    *rybinary,
                              const int    *maxK,
                              const int    *nrTnn,
                              const int    *rncTnn,
                              const int    *rnx,
                              const int    *rdeg,
                                    double *Ey,
                                    double *beta)
{
    C_llm_1D_beta(Tnn_x, Tnn_y, Tnn_w, maxK, nrTnn, rncTnn, rdeg, beta);
    C_wpredict_1D(beta, Tnn_i, Tnn_w, Tnn_x, maxK, nrTnn, rncTnn, rdeg, rnx, rybinary, Ey);
}

/*!
    \brief Creates BB CI's of a 1D rllm model's Ey esimates.

    \param Tnn_i      A matrix of indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
    \param Tnn_w      A matrix of weights over K nearest neighbors of each element of the grid. Weights must sum up to 1!
    \param Tnn_x      A matrix of x values over K nearest neighbors of each element of the grid.
    \param Tnn_y      A matrix of y values over K nearest neighbors of each element of the grid.
    \param maxK       An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
    \param rnrTnn     A reference to the number of rows of the above Tnn_ matrices.
    \param rncTnn     A reference to the number of columns of the above Tnn_ matrices.
    \param rnx        A reference to the number of elements of x.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param rnBB       A reference to the number of Bayesian bootstrap iterations.
    \param Ey         Ey estimates.

    \param Ey_CI     An output array of credible interval radii estimates.

*/
void C_llm_1D_fit_and_predict_BB_CrI(const int    *Tnn_i,
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
                                     const double *Ey,
                                           double *Ey_CI) {
    int nrTnn     = rnrTnn[0];
    int ncTnn = rncTnn[0];
    int nx    = rnx[0];
    int nBB   = rnBB[0];

    int deg = rdeg[0];    // degree of x in the linear models
    int ncoef = deg + 1; // number of rows of beta; number of columns is ncTnn

    // array holding coefficients of the local linear models
    double *beta = (double*)calloc(ncoef * ncTnn, sizeof(double));
    CHECK_PTR(beta);

    // An array carrying Bayesian boostraped columns
    double *lambda = (double*)malloc( nrTnn * sizeof(double));
    CHECK_PTR(lambda);

    // This holds BB absolute errors for a given index of ncTnn
    double *bb = (double*)malloc( nBB * sizeof(double));
    CHECK_PTR(bb);

    // This holds BB absolute Rf_error estimates
    double *ae = (double*)malloc( nx * nBB * sizeof(double));
    CHECK_PTR(ae);

    // BB weights array
    double *Tnn_BB_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_BB_w);

    // BB Ey estimate
    double *bbEy = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(bbEy);

    int G;    // holds maxK[i] + 1
    int jnrTnn;
    int inx;

    double madC = 1.4826; //MAE scaling factor ensuring consistency between MAD
                          //and standard deviation when the data is normally
                          //distributed.
    madC *= 1.96;         // To get the boundary of credible interval (assuming normality).

    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      inx = iboot*nx;

      // Modifying Tnn_w by random samples from the (nrTnn-1)-dimensional simplex
      for ( int j = 0; j < ncTnn; j++ )
      {
        G = maxK[j] + 1;
        jnrTnn = j * nrTnn;

        //rCsimplex(Tnn_w + jK, &G, rC, lambda); // sampling from a Dirichlet distribution with the parameter w=(Tnn_w[0,i], ... , Tnn_w[nr-1, i]). The mean of that distribution is w/sum(w).
        C_rsimplex(Tnn_w + jnrTnn, &G, lambda); // sampling from a Dirichlet distribution with the parameter w=(Tnn_w[0,i], ... , Tnn_w[nr-1, i]). The mean of that distribution is w/sum(w).

        for ( int k = 0; k < G; k++ )
          Tnn_BB_w[k + jnrTnn] = lambda[k];
      }

      C_llm_1D_fit_and_predict(Tnn_i, Tnn_BB_w, Tnn_x, Tnn_y, rybinary, maxK, rnrTnn, rncTnn, rnx, rdeg, bbEy, beta);

      for ( int i = 0; i < nx; i++ )
        ae[i + inx] = fabs( bbEy[i] - Ey[i] );
    }

    for ( int i = 0; i < nx; i++ )
    {
      for ( int iboot = 0; iboot < nBB; iboot++ )
        bb[iboot] = ae[i + iboot*nx];

      Ey_CI[i] = madC * median(bb, nBB);
    }

    free(beta);
    free(lambda);
    free(bb);
    free(ae);
    free(Tnn_BB_w);
    free(bbEy);
}


/*!
    \brief Creates BB CI's of a 1D rllm model's Ey esimates using global reweighting of the elelements of x.

    \param Tnn_i      A matrix of indices of nrTnn nearest neighbors of each element of the grid, where nrTnn is determined in the parent routine.
    \param Tnn_w      A matrix of weights over nrTnn nearest neighbors of each element of the grid. Weights must sum up to 1!
    \param Tnn_x      A matrix of x values over nrTnn nearest neighbors of each element of the grid.
    \param Tnn_y      A matrix of y values over nrTnn nearest neighbors of each element of the grid.
    \param maxK       An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
    \param rnrTnn     A reference to the number of rows of the above Tnn_ matrices.
    \param rncTnn     A reference to the number of columns of the above Tnn_ matrices.
    \param rnx        A reference to the number of elements of x.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param rnBB       A reference to the number of Bayesian bootstrap iterations.
    \param Ey         Ey estimates.

    \param Ey_CI     An output array of credible interval radii estimates.

*/
void C_llm_1D_fit_and_predict_global_BB_CrI(const int    *Tnn_i,
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
                                            const double *Ey,
                                                  double *Ey_CI) {
    int nrTnn   = rnrTnn[0];
    int ncTnn  = rncTnn[0];
    int nTnn = ncTnn * nrTnn;
    int nx  = rnx[0];
    int nBB = rnBB[0];

    int deg = rdeg[0];    // degree of x in the linear models
    int ncoef = deg + 1; // number of rows of beta; number of columns is ncTnn

    // An array holding coefficients of the local linear models
    double *beta = (double*)calloc(ncoef * ncTnn, sizeof(double));
    CHECK_PTR(beta);

    // An array carrying Bayesian boostrap of x
    double *lambda = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(lambda);

    // Holds BB absolute errors for a given index of ncTnn
    double *bb = (double*)malloc( nBB * sizeof(double));
    CHECK_PTR(bb);

    // Holds BB absolute Rf_error estimates
    double *ae = (double*)malloc( nx * nBB * sizeof(double));
    CHECK_PTR(ae);

    // BB weights array
    double *Tnn_BB_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_BB_w);

    // BB Ey estimate
    double *bbEy = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(bbEy);


    double madC = 1.4826; //MAE scaling factor ensuring consistency between MAD
                          //and standard deviation when the data is normally
                          //distributed.
    madC *= 1.96;         // To get the boundary of credible interval (assuming normality).
    int inx;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      inx = iboot*nx;

      C_runif_simplex(rnx, lambda);

      for ( int j = 0; j < nTnn; j++ )
        Tnn_BB_w[j] = lambda[Tnn_i[j]] * Tnn_w[j];

      C_llm_1D_fit_and_predict(Tnn_i, Tnn_BB_w, Tnn_x, Tnn_y, rybinary, maxK, rnrTnn, rncTnn, rnx, rdeg, bbEy, beta);

      for ( int i = 0; i < nx; i++ )
        ae[i + inx] = fabs( bbEy[i] - Ey[i] );
    }

    for ( int i = 0; i < nx; i++ )
    {
      for ( int iboot = 0; iboot < nBB; iboot++ )
        bb[iboot] = ae[i + iboot*nx];

      Ey_CI[i] = madC * median(bb, nBB);
    }

    free(beta);
    free(lambda);
    free(bb);
    free(ae);
    free(Tnn_BB_w);
    free(bbEy);
}

/*!
  \brief Creates BB predictions of rllm.1D() model's Ey esimates using global reweighting of the elelements of x.

  \param Tnn_i      A matrix of indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
  \param Tnn_w      A matrix of weights over K nearest neighbors of each element of the grid. Weights must sum up to 1!
  \param Tnn_x      A matrix of x values over K nearest neighbors of each element of the grid.
  \param Tnn_y      A matrix of y values over K nearest neighbors of each element of the grid.
  \param maxK       An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
  \param rnrTnn     A reference to the number of rows of the above Tnn_ matrices.
  \param rncTnn     A reference to the number of columns of the above Tnn_ matrices.
  \param rnx        A reference to the number of elements of x.
  \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
  \param rnBB       A reference to the number of Bayesian bootstrap iterations.

  \param gbbEy      An output array of BB estimates of Ey.

*/
void C_llm_1D_fit_and_predict_global_BB(const int    *Tnn_i,
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
                                        double *gbbEy)
{
    int nrTnn   = rnrTnn[0];
    int ncTnn  = rncTnn[0];
    int nTnn = ncTnn * nrTnn;
    int nx  = rnx[0];
    int nBB = rnBB[0];

    int deg = rdeg[0];    // degree of x in the linear models
    int ncoef = deg + 1; // number of rows of beta; number of columns is ncTnn

    // array holding coefficients of the local linear models
    double *beta = (double*)calloc(ncoef * ncTnn, sizeof(double));
    CHECK_PTR(beta);

    // An array carrying Bayesian boostrap of x
    double *lambda = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(lambda);

    // BB weights array
    double *Tnn_BB_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_BB_w);

    // BB Ey estimate
    double *bbEy = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(bbEy);

    int inx;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      inx = iboot*nx;

      C_runif_simplex(rnx, lambda);

      for ( int j = 0; j < nTnn; j++ )
        Tnn_BB_w[j] = lambda[Tnn_i[j]] * Tnn_w[j];

      C_llm_1D_fit_and_predict(Tnn_i, Tnn_BB_w, Tnn_x, Tnn_y, rybinary, maxK, rnrTnn, rncTnn, rnx, rdeg, bbEy, beta);

      for ( int i = 0; i < nx; i++ )
        gbbEy[i + inx] = bbEy[i];
    }

    free(beta);
    free(lambda);
    free(Tnn_BB_w);
    free(bbEy);
}

/*!
    \brief Creates BB quantile-based CI's for rllm.1D() model's Ey esimates.

    \param rybinary  A reference to a logical variable. If TRUE, predicted Ey's will be trimmed to the closed interval \[0, 1\].
    \param Tnn_i      A matrix of indices of nrTnn nearest neighbors of each element of the grid, where nrTnn is determined in the parent routine.
    \param Tnn_w      A matrix of weights over nrTnn nearest neighbors of each element of the grid. Weights must sum up to 1!
    \param Tnn_x      A matrix of x values over nrTnn nearest neighbors of each element of the grid.
    \param Tnn_y      A matrix of y values over nrTnn nearest neighbors of each element of the grid.
    \param maxK       An array of indices indicating the range where weights are not 0. Indices <= maxK[i] have weights > 0.
    \param rnrTnn     A reference to the number of rows of the above Tnn_ matrices.
    \param rncTnn     A reference to the number of columns of the above Tnn_ matrices.
    \param rnx        A reference to the number of elements of x.
    \param rdeg       A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param rnBB       A reference to the number of Bayesian bootstrap iterations.
    \param ralpha     A reference to the confidence level.

    \param Ey_CI      An output array of BB estimates of Ey CI's.

*/
void C_llm_1D_fit_and_predict_global_BB_qCrI(const int    *rybinary,
                                             const int    *Tnn_i,
                                             const double *Tnn_w,
                                             const double *Tnn_x,
                                             const double *Tnn_y,
                                             const int    *maxK,
                                             const int    *rnrTnn,
                                             const int    *rncTnn,
                                             const int    *rnx,
                                             const int    *rdeg,
                                             const int    *rnBB,
                                             const double *ralpha,
                                                   double *Ey_CI)
{
    //
    // Generating BB Ey's
    //
    int ybinary  = rybinary[0];
    int nx       = rnx[0];
    int nBB      = rnBB[0];

    // BB Ey estimate
    double *bbEy = (double*)malloc( nx * nBB * sizeof(double));
    CHECK_PTR(bbEy);

    C_llm_1D_fit_and_predict_global_BB(Tnn_i, Tnn_w, Tnn_x, Tnn_y, rybinary, maxK, rnrTnn, rncTnn, rnx, rdeg, rnBB, bbEy);

    if ( ybinary )
    {
      int n = nx * nBB;
      for ( int j = 0; j < n; j++ )
      {
        if ( bbEy[j] < 0 ){
          bbEy[j] = 0;
        } else if ( bbEy[j] > 1 ){
          bbEy[j] = 1;
        }
      }
    }

    //
    // Computing (p, 1-p) quantiles of of each row of bbEy
    //
    double alpha = ralpha[0];
    double p = alpha / 2;
    double probs[] = {p, 1-p};
    int two = 2;

    double *bb = (double*)malloc( nBB * sizeof(double)); // Holds BB Ey estimates of each row of bbEy
    CHECK_PTR(bb);

    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < nBB; j++ )
        bb[j] = bbEy[i + nx*j];

      C_quantiles(bb, rnBB, probs, &two, Ey_CI + 2*i);
    }

    free(bbEy);
    free(bb);
}


/*!
    \brief Performs cross-validation of llm_1D() for degrees 1 and 2 returning MAE.

    This is a modifications of cv_mae_2D().

    \param rnfolds   A reference to the number of fold in cross-validation.
    \param rnreps    A reference to the number of repetitions of cross-validation.
    \param y         A response variable of length nx.
    \param Tnn_i     A matrix of indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
    \param Tnn_d     A matrix of distances over K nearest neighbors of each element of the grid.
    \param Tnn_x     A matrix of x values over Tnn_i.
    \param Tnn_y     A matrix of y values over Tnn_i.
    \param rnrTnn    A reference to the number of rows of the Tnn_ arrays.
    \param rncTnn    A reference to the number of columns of the Tnn_ arrays.
    \param nx        The length of x.
    \param rdeg      A reference to the degree of the polynomial of x in the linear regresion models. The only allowed values are 1 and 2.
    \param rbw       A reference to the value of bandwidth parameter.
    \param rminK     A reference to the mininum __number_ of elements in the rows of Tnn_i/Tnn_d that need to have weights > 0.
    \param rikernel  The integer index of a kernel used for generating weights.

    \param rMAE      A reference to estimated MAE.
*/
void C_cv_mae_1D(const int    *rnfolds,
                 const int    *rnreps,
                 const double *y,
                 const int    *Tnn_i,
                 const double *Tnn_d,
                 const double *Tnn_x,
                 const double *Tnn_y,
                 const int    *rybinary,
                 const int    *rnrTnn,
                 const int    *rncTnn,
                 const int    *rnx,
                 const int    *rdeg,
                 const double *rbw,
                 const int    *rminK,
                 const int    *rikernel,
                       double *rMAE) {

    int nfolds = rnfolds[0];  // number of folds
    int nreps  = rnreps[0];   // number of repetitions of CV
    int nx     = rnx[0];      // number of rows of X before trasposition
    int ncTnn  = rncTnn[0];   // grid size
    int K      = rnrTnn[0];   // number of nearest neighbors = number of rows of Tnn_ matrices
    int deg    = rdeg[0];     // degree of x in the linear models
    int minK   = rminK[0];    // minimal number of points of X present in each local neighbor of each X.grid point
    int ncoef  = deg + 1;     // number of rows of beta; number of columns is ncTnn

    // Mean values of y over elements of the grid.
    double *Ey = (double*)malloc(nx * nreps * sizeof(double));
    CHECK_PTR(Ey);

    // array holding coefficients of the local linear models
    double *beta = (double*)calloc(ncoef * ncTnn, sizeof(double));
    CHECK_PTR(beta);

    // weights array
    double *Tnn_w = (double*)calloc(K * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_w);

    // complementary weights array
    double *Tnn_cw = (double*)calloc(K * ncTnn, sizeof(double));
    CHECK_PTR(Tnn_cw);

    // Fold (F) (ind)ices indicator array
    int *Find = (int*)calloc(nx, sizeof(int));
    CHECK_PTR(Find);

    // Number of fold indices found in Tnn_i[,j] for each j
    int *nFinds = (int*)calloc(ncTnn, sizeof(int));
    CHECK_PTR(nFinds);

    // Turning minK into array with minK_a[j] = minK + nFinds[j] = the ninimal number of elements in the rows of Tnn_i/Tnn_d that need to have their weights > 0.
    int *minK_a = (int*)calloc(ncTnn, sizeof(int));
    CHECK_PTR(minK_a);

    // Allocating memory for the maxK parameter of columnwise_weighting; which
    // is an array of indices indicating the range where weights are > 0.
    // Indices <= maxK[i] have weights > 0.
    int *maxK = (int*)calloc(ncTnn, sizeof(int));
    CHECK_PTR(maxK);

    // Bandwidths
    double *bws = (double*)malloc(ncTnn * sizeof(double));
    CHECK_PTR(bws);

    // holds results of prediction over hold out set
    double *hoEy = (double*)calloc(nx, sizeof(double));
    CHECK_PTR(hoEy);

    iarray_t *folds; // array with nfolds elements each being array_t holding an array of size ~ 0.1*|X|
    int *F;          // array holding indices of the current fold
    int nF;          // number of elements of I
    int rep_nx;      // holds nx * rep
    int G;           // (index of the last non-zero element) + 1
    int jK;          // j * K

    for ( int rep = 0; rep < nreps; rep++ ) {

      rep_nx = nx * rep;
      folds = get_folds(nx, nfolds);

      for ( int ifold = 0; ifold < nfolds; ifold++ )
      {
        F  = folds[ifold].x;    // indices of the fold with values between 0 and nx-1
        nF = folds[ifold].size; // number of the elements in F

        for ( int j = 0; j < nF; j++ )
          Find[F[j]] = 1;

        // Finding number of fold indices in Tnn_i[,j]; the corresponding number of fold indices is nFinds[j]
        for ( int j = 0; j < ncTnn; j++ )
        {
          jK = j * K;
          for ( int k = 0; k < minK; k++ )
            nFinds[j] +=  Find[Tnn_i[k + jK]];
        }

        // Turning Kmin into array Kmin_a adding nFinds
        for ( int j = 0; j < ncTnn; j++ )
          minK_a[j] = minK + nFinds[j];

        // normalizing Tnn_d using minK_a
        C_get_bws_with_minK_a(Tnn_d, rnrTnn, rncTnn, minK_a, rbw, bws);

        // Appling kernel to the normalized distances
        C_columnwise_weighting(Tnn_d, rnrTnn, rncTnn, rikernel, bws, maxK, Tnn_w);

        // Modifying Tnn_w and creating Tnn_cw.
        // Tnn_w needs to be 0 at the indices of F
        for ( int j = 0; j < ncTnn; j++ )
        {
          jK = j * K;
          G = maxK[j] + 1;

          for ( int k = 0; k < G; k++ )
          {
            Tnn_cw[k + jK] = Tnn_w[k + jK] * Find[Tnn_i[k + jK]]; // This sets to 0 terms not from F
            Tnn_w[k + jK] *= 1 - Find[Tnn_i[k + jK]]; // Find[Tnn_i[k + jK]] = 1,
                                                    // when Tnn_i[k + jK] is from
                                                    // F, so we are multiplying
                                                    // by 0 terms corresponding
                                                    // to elements of F
          }
        }

        // Tnn_w and Tnn_cw are not normalized
        //
        // For model building weights MUST be normalized. Thus, we need to normalize Tnn_w.
        // For prediction, the weights do not have to be normalized.
        //
        // Building models using Tnn_w. Since elements of F have weights 0, they
        // do not participate in the model building.
        C_llm_1D_beta(Tnn_x, Tnn_y, Tnn_w, maxK, rnrTnn, rncTnn, rdeg, beta);

        // Generating predictions using Tnn_cw
        C_wpredict_1D(beta, Tnn_i, Tnn_cw, Tnn_x, maxK, rnrTnn, rncTnn, rdeg, rnx, rybinary, hoEy);

        for ( int j = 0; j < nF; j++ )
          Ey[F[j] + rep_nx] = hoEy[F[j]];

        // Resetting Find
        for ( int j = 0; j < nF; j++ )
          Find[F[j]] = 0; // reseting the fold indices to 0

        for ( int j = 0; j < ncTnn; j++ )
          nFinds[j] = 0;

      } // END OF for ( int ifold = 0; ifold < nfolds; ifold++ )

      // freeing memory allocated to folds[i].x and folds from within get_folds(nx, nfolds)
      for ( int i = 0; i < nfolds; i++ )
        free(folds[i].x);
      free(folds);

    } // END OF for ( int rep = 0; rep < nreps; rep++ )

    // MAE loop
    double mae = 0;
    for ( int rep = 0; rep < nreps; rep++ )
    {
      rep_nx = nx * rep;

      for ( int i = 0; i < nx; i++ )
        mae += fabs( y[i] - Ey[i + rep_nx] );
    }

    *rMAE = mae / ( nreps * nx );

    free(Ey);
    free(beta);
    free(Tnn_w);
    free(Tnn_cw);
    free(Find);
    free(nFinds);
    free(minK_a);
    free(maxK);
    free(bws);
    free(hoEy);
}

/*
  Given a vector, bws, of bandwidths this routine estimates a matrix of Eyg's
  each column corresonding for the Eyg estimate corresponding to different
  bandwidth.

  loop over bws values
  at each iterations
  1) Calculate nn.w
  2) fit local linear models at the grid points using the calculated nn.w's
  3) use the resulting model coefficients to predict Eyg

  \param bws         An array of bandwidth values.
  \param rnbws       A reference to the length of bws.
  \param Tnn_i       A matrix of K nearest neighbor indices of each element of the grid, where K is determined in the parent routine.
  \param Tnn_d       A matrix of distances to K nearest neighbors.
  \param Tnn_x       A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
  \param Tnn_y       A matrix of y values over K nearest neighbors of each element of the grid.
  \param rnrTnn      A reference to the number of rows of the above three matrices.
  \param rncTnn      A reference to the number of columns of the above three matrices.
  \param max_K       A vector of __indices__ indicating the range where weights are not 0.
  \param rdegree     A reference to the degree of the local models.
  \param rminK       A reference to the minimal number of points in each set of NNs.
  \param Tgrid_nn_i  A matrix of grid asociated NN indices.
  \param Tgrid_nn_d  A matrix of grid asociated distances to grid NNs.
  \param Tgrid_nn_x  A matrix of grid asociated x values over NNs.
  \param rnrTgrid_nn A reference to the number of rows of the above three matrices.
  \param rncTgrid_nn A reference to the number of columns of the above three matrices.
  \param Eygs        An output matrix of Eyg vectors (columns of Eygs) associated with bws values.

*/
void C_get_Eygs(const double *bws,
                const int    *rn_bws,
                const int    *Tnn_i,
                const double *Tnn_d,
                const double *Tnn_x,
                const double *Tnn_y,
                const int    *rybinary,
                const int    *rnrTnn,
                const int    *rncTnn,
                const int    *rdegree,
                const int    *rminK,
                const int    *Tgrid_nn_i,
                const double *Tgrid_nn_d,
                const double *Tgrid_nn_x,
                const int    *rnrTgrid_nn,
                const int    *rncTgrid_nn,
                      double *Eygs) {

    int ncTnn      = rncTnn[0];      // grid size
    int nrTnn      = rnrTnn[0];      // number of nearest neighbors = number of rows of Tnn_x ( as Tnn_x = t(nn.x) in R), Tnn_y, Tnn_w
    int deg        = rdegree[0];     // degree of local linear models
    int ncoef      = deg + 1;
    int ncTgrid_nn = rncTgrid_nn[0];
    int nrTgrid_nn = rnrTgrid_nn[0];  // number of nearest neighbors = number of rows of Tnn_x ( as Tnn_x = t(nn.x) in R), Tnn_y, Tnn_w

    // local radii
    double *Tnn_r = (double*)malloc(nrTnn * ncTnn * sizeof(double));
    CHECK_PTR(Tnn_r);

    // weights array
    double *Tnn_w = (double*)malloc(nrTnn * ncTnn * sizeof(double));
    CHECK_PTR(Tnn_w);

    // grid weights array
    double *Tgrid_nn_w = (double*)malloc(nrTgrid_nn * ncTgrid_nn * sizeof(double));
    CHECK_PTR(Tgrid_nn_w);

    // maxK
    int *max_K = (int*)malloc(ncTnn * sizeof(int));
    CHECK_PTR(max_K);

    // grid max_K
    int *grid_max_K = (int*)malloc(ncTnn * sizeof(int));
    CHECK_PTR(grid_max_K);

    // Allocating memory for beta
    double *beta = (double*)malloc( ncoef * ncTnn * sizeof(double) );
    CHECK_PTR(beta);

    int ikernel = 1;

    //for ( int i_bws = 0; i_bws < n_bws; i_bws++ )
    int i_bws = 0;
    {
      // Calculate Tnn_r
      C_get_bws(Tnn_d, rnrTnn, rncTnn, rminK, bws + i_bws, Tnn_r);

      // Calculate Tnn_w and max_K
      C_columnwise_weighting(Tnn_d, rnrTnn, rncTnn, &ikernel, Tnn_r, max_K, Tnn_w);

      // Fit local linear models at the grid points using the modified nn.w's
      C_llm_1D_beta(Tnn_x, Tnn_y, Tnn_w, max_K, rnrTnn, rncTnn, rdegree, beta);

      // Getting Tgrid_nn_w and grid_max_K
      C_columnwise_weighting(Tgrid_nn_d, rnrTgrid_nn, rncTgrid_nn, &ikernel, Tnn_r, grid_max_K, Tgrid_nn_w);

      // Use the resulting model coefficients to predict Eyg
      C_wpredict_1D(beta, Tgrid_nn_i, Tgrid_nn_w, Tgrid_nn_x, grid_max_K, rnrTgrid_nn, rncTgrid_nn, rdegree, rncTnn, rybinary, Eygs + i_bws*ncTnn);
    }

    free(max_K);
    free(grid_max_K);
    free(Tnn_r);
    free(Tnn_w);
    free(Tgrid_nn_w);
    free(beta);
}


/*
  Generates Bayesian bootstrap estimates of Eyg

  loop over n.BB iterations
  at each iterations
  1) select random sample, lambda, from a simple of dim nx-1
  2) change nn.w by lambda
  3) fit local linear models at the grid points using the modified nn.w's
  4) use the resulting model coefficients to predict bb.Eyg (BB version of Eyg)

  \param rn_BB       A reference to the number of Bayesian bootstrap iterations.
  \param Tnn_x       A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
  \param Tnn_y       A matrix of y values over K nearest neighbors of each element of the grid.
  \param Tnn_w       A matrix of NN weights.
  \param rxn         A reference to the number of elements of x.
  \param rnrTnn      A reference to the number of rows of the above three matrices.
  \param rncTnn      A reference to the number of columns of the above three matrices.
  \param max_K       A vector of __indices__ indicating the range where weights are not 0.
  \param rdegree     A reference to the degree of the local models.
  \param Tgrid_nn_i  A matrix of grid asociated NN indices.
  \param Tgrid_nn_w  A matrix of grid asociated NN weights.
  \param Tgrid_nn_x  A matrix of grid asociated x values over NNs.
  \param rnrTgrid_nn A reference to the number of rows of the above three matrices.
  \param rncTgrid_nn A reference to the number of columns of the above three matrices.
  \param grid_max_K  A vector of grid asociated max.K values.

*/
void C_get_BB_Eyg(const int    *rn_BB,
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
                        double *bb_dEyg,
                        double *bb_Eyg)
{

    int n_BB       = rn_BB[0];       // number of Bayesian bootstrap iterations
    int nx         = rnx[0];
    int ncTnn      = rncTnn[0];      // grid size
    int nrTnn      = rnrTnn[0];      // number of nearest neighbors = number of rows of Tnn_x ( as Tnn_x = t(nn.x) in R), Tnn_y, Tnn_w
    int nTnn       = ncTnn * nrTnn;
    int deg        = rdegree[0];     // degree of local linear models
    int ncoef      = deg + 1;
    //int ncTgrid_nn = rncTgrid_nn[0]; // grid size
    //int nrTgrid_nn = rnrTgrid_nn[0]; // number of nearest neighbors = number of rows of Tnn_x ( as Tnn_x = t(nn.x) in R), Tnn_y, Tnn_w

    // An array carrying Bayesian boostrap probabilities of x
    double *lambda = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(lambda);

    // BB weights array
    double *bb_Tnn_w = (double*)calloc(nrTnn * ncTnn, sizeof(double));
    CHECK_PTR(bb_Tnn_w);

    // Allocating memory for beta
    double *bb_beta = (double*)malloc( ncoef * ncTnn * sizeof(double) );
    CHECK_PTR(bb_beta);

    for ( int i_BB = 0; i_BB < n_BB; i_BB++ )
    {
      // Select random sample, lambda, from a simple of dim nx-1
      C_runif_simplex(rnx, lambda);

      // Change nn.w by lambda
      for ( int j = 0; j < nTnn; j++ )
        bb_Tnn_w[j] = lambda[Tnn_i[j]] * Tnn_w[j];

      // Fit local linear models at the grid points using the modified nn.w's
      C_llm_1D_beta(Tnn_x, Tnn_y, bb_Tnn_w, max_K, rnrTnn, rncTnn, rdegree, bb_beta);

      // Use the resulting model coefficients to predict bb.Eyg (BB version of Eyg)
      C_wpredict_1D(bb_beta, Tgrid_nn_i, Tgrid_nn_w, Tgrid_nn_x, grid_max_K, rnrTgrid_nn, rncTgrid_nn, rdegree, rncTnn, rybinary, bb_Eyg + i_BB*ncTnn);

      for ( int j = 0; j < ncTnn; j++ )
        bb_dEyg[j + i_BB*ncTnn] = bb_beta[1 + j*ncoef];
    }

    free(lambda);
    free(bb_Tnn_w);
    free(bb_beta);
}


/*
  Estimates credible interals of Eyg values using Bayesian bootstrap.

  \param rybinary    A reference to a logical variable. If TRUE, predicted Eyg's will be trimmed to the closed interval \[0, 1\].
  \param rn_BB       A reference to the number of Bayesian bootstrap iterations.
  \param Tnn_x       A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
  \param Tnn_y       A matrix of y values over K nearest neighbors of each element of the grid.
  \param Tnn_w       A matrix of NN weights.
  \param rxn         A reference to the number of elements of x.
  \param rnrTnn      A reference to the number of rows of the above three matrices.
  \param rncTnn      A reference to the number of columns of the above three matrices.
  \param max_K       A vector of __indices__ indicating the range where weights are not 0.
  \param rdegree     A reference to the degree of the local models.
  \param Tgrid_nn_i  A matrix of grid asociated NN indices.
  \param Tgrid_nn_w  A matrix of grid asociated NN weights.
  \param Tgrid_nn_x  A matrix of grid asociated x values over NNs.
  \param rnrTgrid_nn A reference to the number of rows of the above three matrices.
  \param rncTgrid_nn A reference to the number of columns of the above three matrices.
  \param grid_max_K  A vector of grid asociated max.K values.
  \param ralpha      A reference to the significance level.
*/
void C_get_Eyg_CrI(const int    *rybinary,
                   const int    *rn_BB,
                   const int    *Tnn_i,
                   const double *Tnn_x,
                   const double *Tnn_y,
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
                   const double *ralpha,
                         double *Eyg_CrI)
{
    int n_BB     = rn_BB[0];       // number of Bayesian bootstrap iterations
    int ncTnn    = rncTnn[0];      // grid size
    //int ybinary  = rybinary[0];
    int n        = ncTnn * n_BB;

    double *bb_Eyg = (double*)malloc( n * sizeof(double));
    CHECK_PTR(bb_Eyg);

    // Allocating memory for beta
    double *bb_dEyg = (double*)malloc( ncTnn * n_BB * sizeof(double) );
    CHECK_PTR(bb_dEyg);

    C_get_BB_Eyg(rn_BB,
                 Tnn_i,
                 Tnn_x,
                 Tnn_y,
                 rybinary,
                 Tnn_w,
                 rnx,
                 rnrTnn,
                 rncTnn,
                 max_K,
                 rdegree,
                 Tgrid_nn_i,
                 Tgrid_nn_x,
                 Tgrid_nn_w,
                 rnrTgrid_nn,
                 rncTgrid_nn,
                 grid_max_K,
                 bb_dEyg,
                 bb_Eyg);


    #if 0
    if ( ybinary )
    {
      for ( int j = 0; j < n; j++ )
      {
        if ( bb_Eyg[j] < 0 ){
          bb_Eyg[j] = 0;
        } else if ( bb_Eyg[j] > 1 ){
          bb_Eyg[j] = 1;
        }
      }
    }
    #endif

    //
    // Computing (p, 1-p) quantiles of of each row of bb_Eyg
    //
    double alpha = ralpha[0];
    double p = alpha / 2;
    double probs[] = {p, 1-p};
    int two = 2;

    double *bb = (double*)malloc( n_BB * sizeof(double)); // Holds BB Ey estimates of each row of bb_Eyg
    CHECK_PTR(bb);

    for ( int i = 0; i < ncTnn; i++ )
    {
      for ( int j = 0; j < n_BB; j++ )
        bb[j] = bb_Eyg[i + ncTnn*j];

      C_quantiles(bb, rn_BB, probs, &two, Eyg_CrI + 2*i);
    }

    free(bb_Eyg);
    free(bb_dEyg);
    free(bb);
}

// structure to hold a value of a function and the associated weight
typedef struct {
  double y;
  double w;
} yw_t;


// This allows sorting yw_t's in descending order
int cmp_desc_yw_t (const void *a, const void *b)
{
  yw_t *ya = (yw_t *)a;
  yw_t *yb = (yw_t *)b;

  if ( ya->w < yb->w ){
    return 1;
  } else if ( ya->w > yb->w ){
    return -1;
  } else {
    return 0;
  }
}

/*!

  \brief Returns the weighted mean of y values of yw_t array weighted by the
  corresponding w component of the struct.

  \param yw  A pointer to a yw_t array.
  \param n   The number of elements of yt.
*/
double yw_t_wmean( yw_t *yw, int n )
{
    double s = 0;
    double yws = 0;
    for ( int i = 0; i < n; i++ )
      s += yw[i].w;

    if ( s > 0 )
    {
      for ( int i = 0; i < n; i++ )
        yws += yw[i].y * yw[i].w;

      yws /= s;

    } else {

      for ( int i = 0; i < n; i++ )
        yws += yw[i].y;

      yws /= n;
    }

    return yws;
}

/*!
    \brief Interpolates the mean values, Eyg, of y defined over a grid using nn_w weights.

    \param Eyg        Mean values of y over grid elements.
    \param nn_i       An array associated with a K-by-ng matrix of nearest neighbors (NN) indices.
    \param nn_d       An array associated with a K-by-ng matrix of nearest neighbors (NN) distances.
    \param rK         A reference of the number of rows of the nn_ arrays.
    \param rng        A reference of the number of columns of the nn_ arrays.

    \param rnNN       A reference of the number of NN to use for finding interpolated
                      values. In principle, if we set nNN to the dimension of X, then it will be a
                      generalization of a linear interpolation. NOTE: nNN cannot be greater than ng!

    \param rnx        A reference of the length of x (points over which Ey is estimated).

    \param Ey         An output array of weighted mean values of y estimates using Eyg with nn_w weights.
*/
void C_interpolate_Eyg(const double *Eyg,
                       const int    *nn_i,
                       const double *nn_d,
                       const int    *rK,
                       const int    *rng,
                       const int    *rnNN,
                       const int    *rnx,
                             double *Ey)
{
    int K    = rK[0];
    int ng   = rng[0];
    int nNN  = rnNN[0];
    int nx = rnx[0]; // nn_i[i]'s are in the range 0 .. (nx-1), as they index values of x, over which Ey is estimated

    // initializing Ey
    for ( int i = 0; i < nx; i++ )
      Ey[i] = 0;

    int *nn_ix = (int*)malloc( ng * nx * sizeof(int) ); // nn_ix - index of x in nn_i
    CHECK_PTR(nn_ix);

    // initializing nn_ix
    int n = ng * nx;
    for ( int i = 0; i < n; i++ )
      nn_ix[i] = -1;

    int ix;
    int igx;
    int iK;  // ig * K
    for ( int ig = 0; ig < ng; ig++ )
    {
      iK  = ig * K;
      igx = ig * nx;
      for ( int j = 0; j < K; j++ )
      {
        ix = nn_i[j + iK]; // index of x in the i-th column and j-th row of nn_i = the j-th NN of the i-th xg
        nn_ix[ ix + igx ] = j;
      }
    }

    // Estimating NN weights
    double *bws = (double*)malloc(ng * sizeof(double));
    CHECK_PTR(bws);

    for ( int ig = 0; ig < ng; ig++ )
      bws[ig] = 1.0;

    double *nn_w = (double*)calloc(K * ng, sizeof(double));
    CHECK_PTR(nn_w);

    int *maxK = (int*)calloc(ng, sizeof(int));
    CHECK_PTR(maxK);

    int ikernel = (int)EPANECHNIKOV;
    C_columnwise_weighting(nn_d, rK, rng, &ikernel, bws, maxK, nn_w);

    yw_t *yw = (yw_t*)malloc(ng * sizeof(yw_t));
    CHECK_PTR(yw);

    int j;
    for ( int ix = 0; ix < nx; ix++ )
    {
      for ( int ig = 0; ig < ng; ig++ )
      {
        j = nn_ix[ix + ig*nx];
        if ( j > -1 )
        {
          yw[ig].y = Eyg[ig];
          yw[ig].w = nn_w[ j + ig*K];

        } else {

          yw[ig].y = 0;
          yw[ig].w = 0;
        }

      }

      qsort( yw, ng, sizeof(yw_t), cmp_desc_yw_t );

      Ey[ix] = yw_t_wmean( yw, nNN );
    }

    free(nn_ix);
    free(bws);
    free(nn_w);
    free(maxK);
    free(yw);
}


/*!
    \brief Performs cross-validation of llm_xD() for degree 0 model cvEy estimates.

    \param rnfolds   A reference to the number of fold in cross-validation.
    \param rnreps    A reference to the number of repetitions of cross-validation.
    \param rnNN      A reference to the number of nearest neighbors used in interpolate_Eyg() to estimate Ey given Ey.grid.
    \param rybinary  A reference to a binar indicator, ybinary, such that, ybinary=1 restricts the values of Eyg (and hence Ey) to [0,1].
    \param nn_i      An array corresponding to a K-by-ng matrix of indices of K nearest neighbors of each element of a grid.
    \param nn_d      An array corresponding to a K-by-ng matrix of distances over K nearest neighbors from the elements of the grid to X.
    \param y         A response variable.
    \param rK        A reference to the number of rows of the nn_ arrays.
    \param rng       A reference to the number of columns of the nn_ arrays.
    \param rnrX      A reference to the number of rows of X.
    \param rbw       A reference to the value of bandwidth parameter.
    \param rminK     A reference to the mininum __number_ of elements in the rows of nn_i/nn_d that need to have weights > 0.
    \param rikernel  The integeer index of a kernel used for generating weights.

    \param cvEy      The output variable of predicted values of the mean of y over
                     hold-out folds over nreps replications of CV; dim nrX-by-nreps.
*/
void C_cv_deg0_llm(const int    *rnfolds,
                   const int    *rnreps,
                   const int    *rnNN,
                   const int    *rybinary,
                   const int    *nn_i,
                   const double *nn_d,
                   const double *y,
                   const int    *rK,
                   const int    *rng,
                   const int    *rnrX,
                   const double *rbw,
                   const int    *rminK,
                   const int    *rikernel,
                         double *cvEy)
{
    int nfolds = rnfolds[0]; // number of folds
    int nreps  = rnreps[0];  // number of repetitions of CV
    int nrX    = rnrX[0];    // number of rows of X before trasposition
    int ng     = rng[0];     // grid size
    int K      = rK[0];      // number of nearest neighbors = number of rows of nn_ matrices
    int minK   = rminK[0];   // minimal number of points of X present in each local neighbor of each X.grid point
    int nNN    = rnNN[0];
    int ybinary = rybinary[0];

    // Mean values of y over elements of the grid.
    double *Eyg = (double*)malloc(ng * sizeof(double));
    CHECK_PTR(Eyg);

    // weights array
    double *nn_w = (double*)calloc(K * ng, sizeof(double));
    CHECK_PTR(nn_w);

    // complementary weights array
    double *nn_cw = (double*)calloc(K * ng, sizeof(double));
    CHECK_PTR(nn_cw);

    // Fold (F) (ind)ices indicator array
    int *Find = (int*)calloc(nrX, sizeof(int));
    CHECK_PTR(Find);

    // Number of fold indices found in nn_i[,j] for each j
    int *nFinds = (int*)calloc(ng, sizeof(int));
    CHECK_PTR(nFinds);

    // Turning minK into array with minK_a[j] = minK + nFinds[j] = the minimal number of elements of the rows of nn_i/nn_d with weights > 0.
    int *minK_a = (int*)calloc(ng, sizeof(int));
    CHECK_PTR(minK_a);

    // Allocating memory for the maxK parameter of columnwise_weighting
    int *maxK = (int*)calloc(ng, sizeof(int));
    CHECK_PTR(maxK);

    // Bandwidths
    double *bws = (double*)malloc(ng * sizeof(double));
    CHECK_PTR(bws);

    // holds results of prediction over hold out set
    double *hoEy = (double*)calloc(nrX, sizeof(double));
    CHECK_PTR(hoEy);

    // array of y values over NN's indices of nn_i
    double *nn_y = (double*)calloc(K * ng, sizeof(double));
    CHECK_PTR(nn_y);

    //
    // Creating nn_y
    //
    C_columnwise_eval(nn_i, rK, rng, y, nn_y);

    int nF;          // number of elements of I
    iarray_t *folds; // array with nfolds elements each being array_t holding an array of size ~ 0.1*|X|
    int *F;          // array holding indices of the current fold
    int rep_nrX;     // holds nrX * rep
    int G;           // (index of the last non-zero element) + 1
    int jK;          // holds j * K

    for ( int rep = 0; rep < nreps; rep++ )
    {
      folds   = get_folds(nrX, nfolds);
      rep_nrX = nrX * rep;

      for ( int ifold = 0; ifold < nfolds; ifold++ )
      {
        F  = folds[ifold].x;    // indices of the fold with values between 0 and nrX-1
        nF = folds[ifold].size; // number of the elements in F

        for ( int j = 0; j < nF; j++ )
          Find[F[j]] = 1;

        // Finding number of fold indices in nn_i[,j]; the corresponding number of fold indices is nFinds[j]
        for ( int j = 0; j < ng; j++ )
        {
          jK = j * K;
          for ( int k = 0; k < minK; k++ )
            nFinds[j] +=  Find[nn_i[k + jK]];
        }

        // Turning Kmin into array Kmin_a adding nFinds
        for ( int j = 0; j < ng; j++ )
          minK_a[j] = minK + nFinds[j];

        // getting bws' using minK_a
        C_get_bws_with_minK_a(nn_d, rK, rng, minK_a, rbw, bws);

        // Appling kernel to the normalized distances
        C_columnwise_weighting(nn_d, rK, rng, rikernel, bws, maxK, nn_w);

        // Modifying nn_w and creating nn_cw.
        // nn_w needs to be 0 at the indices of F
        for ( int j = 0; j < ng; j++ )
        {
          jK = j * K;
          G = maxK[j] + 1;

          for ( int k = 0; k < G; k++ )
          {
            nn_cw[k + jK] = nn_w[k + jK] * Find[nn_i[k + jK]]; // This sets to 0 terms not from F
            nn_w[k + jK] *= 1 - Find[nn_i[k + jK]]; // Find[nn_i[k + jK]] = 1,
                                                    // when nn_i[k + jK] is from
                                                    // F, so we are multiplying
                                                    // by 0 terms corresponding
                                                    // to elements of F
          }
        }

        // nn_w and nn_cw are not normalized
        //
        // For model building weights MUST be normalized. Thus, we need to normalize nn_w.
        // For prediction, the weights do not have to be normalized.

        // Normalizing nn_w and nn_cw
        double s; // sum of weights in each column
        for ( int j = 0; j < ng; j++ )
        {
          jK = j * K;
          G = maxK[j] + 1;

          // nn_w normalization
          s = 0.0;
          for ( int k = 0; k < G; k++ )
            s += nn_w[k + jK];

          for ( int k = 0; k < G; k++ )
            nn_w[k + jK] /= s;

          // nn_cw normalization
          s = 0.0;
          for ( int k = 0; k < G; k++ )
            s += nn_cw[k + jK];

          for ( int k = 0; k < G; k++ )
            nn_cw[k + jK] /= s;
        }

        // Estimating the mean values, Eyg, of y over elements of the grid.
        C_columnwise_wmean(nn_y, nn_w, maxK, rK, rng, Eyg);

        if ( ybinary == 1 )
        {
          for ( int j = 0; j < ng; j++ )
          {
            if ( Eyg[j] < 0 ){
              Eyg[j] = 0;
            } else if ( Eyg[j] > 1 ){
              Eyg[j] = 1;
            }
          }
        }

        // Generating predicted mean y values of the elements of the fold using nn_cw
        C_interpolate_Eyg(Eyg, nn_i, nn_d, rK, rng, &nNN, rnrX, hoEy);

        for ( int j = 0; j < nrX; j++ )
          if ( !R_FINITE(hoEy[j]) )
          {
            Rprintf("\n\nERROR in cv_deg0_llm() rep=%d ifold=%d: hoEy[%d] is not finite!\n", rep, ifold, j);
            Rf_error("STOP");
          }

        for ( int j = 0; j < nF; j++ )
          cvEy[F[j] + rep_nrX] = hoEy[F[j]];

        // Resetting Find
        for ( int j = 0; j < nF; j++ )
          Find[F[j]] = 0; // reseting the fold indices to 0

        for ( int j = 0; j < ng; j++ )
          nFinds[j] = 0;

      } // END OF for ( int ifold = 0; ifold < nfolds; ifold++ )

      // freeing memory allocated to folds[i].x and folds from within get_folds(nrX, nfolds)
      for ( int i = 0; i < nfolds; i++ )
        free(folds[i].x);
      free(folds);

    } // END OF for ( int rep = 0; rep < nreps; rep++ )

    free(Eyg);
    free(nn_w);
    free(nn_cw);
    free(Find);
    free(nFinds);
    free(minK_a);
    free(maxK);
    free(hoEy);
    free(bws);
}



/*!
    \brief Performs cross-validation of llm_xD() for degree 0 model returning mean absolute error (MAE).

    \param rnfolds   A reference to the number of fold in cross-validation.
    \param rnreps    A reference to the number of repetitions of cross-validation.
    \param rnNN      A reference to the number of nearest neighbors used in interpolate_Eyg() to estimate Ey given Ey.grid.
    \param rybinary  A reference to a binar indicator, ybinary, such that, ybinary=1 restricts the values of Eyg (and hence Ey) to [0,1].
    \param nn_i      An array corresponding to a K-by-ng matrix of indices of K nearest neighbors of each element of a grid.
    \param nn_d      An array corresponding to a K-by-ng matrix of distances over K nearest neighbors from the elements of the grid to X.
    \param y         A response variable.
    \param rK        A reference to the number of rows of the nn_ arrays.
    \param rng       A reference to the number of columns of the nn_ arrays.
    \param rnrX      A reference to the number of rows of X.
    \param rbw       A reference to the value of bandwidth parameter.
    \param rminK     A reference to the mininum __number_ of elements in the rows of nn_i/nn_d that need to have weights > 0.
    \param rikernel  The integeer index of a kernel used for generating weights.
    \param rMAE      A reference to estimated MAE.
*/
void C_cv_deg0_mae(const int    *rnfolds,
                   const int    *rnreps,
                   const int    *rnNN,
                   const int    *rybinary,
                   const int    *nn_i,
                   const double *nn_d,
                   const double *y,
                   const int    *rK,
                   const int    *rng,
                   const int    *rnrX,
                   const double *rbw,
                   const int    *rminK,
                   const int    *rikernel,
                         double *rMAE)
{
    int nreps  = rnreps[0];  // number of repetitions of CV
    int nrX    = rnrX[0];    // number of rows of X before trasposition

    // Mean values of y over elements of the grid.
    double *cvEy = (double*)malloc(nrX * nreps * sizeof(double));
    CHECK_PTR(cvEy);

    C_cv_deg0_llm(rnfolds,
                rnreps,
                rnNN,
                rybinary,
                nn_i,
                nn_d,
                y,
                rK,
                rng,
                rnrX,
                rbw,
                rminK,
                rikernel,
                cvEy);

    // MAE loop
    int rep_nrX; // holds nrX * rep
    double mae = 0;

    for ( int rep = 0; rep < nreps; rep++ )
    {
      rep_nrX = nrX * rep;

      for ( int i = 0; i < nrX; i++ )
        mae += fabs( y[i] - cvEy[i + rep_nrX] );
    }
    mae /= nreps * nrX;

    *rMAE = mae;

    free(cvEy);
}

/*!
    \brief Performs cross-validation of llm_xD() for degree 0 model returning binary loss function y(1-p) + (1-y)p.

    \param rnfolds   A reference to the number of fold in cross-validation.
    \param rnreps    A reference to the number of repetitions of cross-validation.
    \param rnNN      A reference to the number of nearest neighbors used in interpolate_Eyg() to estimate Ey given Ey.grid.
    \param rybinary  A reference to a binar indicator, ybinary, such that, ybinary=1 restricts the values of Eyg (and hence Ey) to [0,1].
    \param nn_i      An array corresponding to a K-by-ng matrix of indices of K nearest neighbors of each element of a grid.
    \param nn_d      An array corresponding to a K-by-ng matrix of distances over K nearest neighbors from the elements of the grid to X.
    \param y         A response variable.
    \param rK        A reference to the number of rows of the nn_ arrays.
    \param rng       A reference to the number of columns of the nn_ arrays.
    \param rnrX      A reference to the number of rows of X.
    \param rbw       A reference to the value of bandwidth parameter.
    \param rminK     A reference to the mininum __number_ of elements in the rows of nn_i/nn_d that need to have weights > 0.
    \param rikernel  The integeer index of a kernel used for generating weights.
    \param rbinloss  A reference to estimated binary loss.
*/
void C_cv_deg0_binloss(const int    *rnfolds,
                       const int    *rnreps,
                       const int    *rnNN,
                       const int    *rybinary,
                       const int    *nn_i,
                       const double *nn_d,
                       const double *y,
                       const int    *rK,
                       const int    *rng,
                       const int    *rnrX,
                       const double *rbw,
                       const int    *rminK,
                       const int    *rikernel,
                             double *rbinloss)
{
    int nreps  = rnreps[0];  // number of repetitions of CV
    int nrX    = rnrX[0];    // number of rows of X before trasposition

    // Mean values of y over elements of the grid.
    double *cvEy = (double*)malloc(nrX * nreps * sizeof(double));
    CHECK_PTR(cvEy);

    C_cv_deg0_llm(rnfolds,
                rnreps,
                rnNN,
                rybinary,
                nn_i,
                nn_d,
                y,
                rK,
                rng,
                rnrX,
                rbw,
                rminK,
                rikernel,
                cvEy);

    // binloss loop
    int rep_nrX; // holds nrX * rep
    double bl = 0;

    for ( int rep = 0; rep < nreps; rep++ )
    {
      rep_nrX = nrX * rep;

      for ( int i = 0; i < nrX; i++ )
        bl +=  y[i]*(1 - cvEy[i + rep_nrX]) + (1 - y[i])*cvEy[i + rep_nrX];
    }
    bl /= nreps * nrX;

    *rbinloss = bl;

    free(cvEy);
}
