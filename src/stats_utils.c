#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gflow_macros.h"
#include "stats_utils.h"

/*!
  Samples `n` points from a uniform distribution within a hypercube of dimension `dim`.

  \param n   The number of points to sample.
  \param dim The dimensionality of the hypercube.
  \param L   A pointer to an array containing the lower bounds for each dimension.
  \param w   The length of each edge in the hypercube.

  \return    A pointer to an array containing the sampled points. The caller is responsible for freeing this memory.
*/
double* runif_hcube(int n, int dim, double* L, double w) {
  // Allocate memory for the points
  double* X = malloc(n * dim * sizeof(double));
  if (X == NULL) {
    perror("runif_hcube(): Could not allocate memory for X");
    return NULL;
  }

  // Generate random points using R's RNG
  GetRNGstate();
  for (int i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      double random_val = unif_rand();
      X[d + i * dim] = L[d] + w * random_val;
    }
  }
  PutRNGstate();

  return X;
}


/*!
   Compares two doubles.
*/
int dcmp( void const *e1, void const *e2 )
{
    double v1 = *( ( double * ) e1 );
    double v2 = *( ( double * ) e2 );

    return ( 0 < v2 - v1) ? -1 : ( v1 - v2 > 0 ) ? 1 : 0;
}


/*!

  \fn median( double *a, int n )

  Median of double array of size n; done in place (array will be sorted)

  \param a  A double array of size n.
  \param n  The length of a.

*/
double median(double *a, int n )
{
    double median = 0;

    if ( n == 1 ){
        return a[0];
    } else if ( n == 0 )
    {
      error("ERROR: median cannot be computed for array of size 0!\n");
    }

    qsort( a, n, sizeof(double), dcmp );

    if( n % 2 == 0 )
        median = ( a[n/2-1] + a[n/2] ) / 2.0;
    else
        median = a[n/2];

    return median;
}

/*!
    \fn double mean( double *x, int n )

    \brief The mean of x

    \param x  - input 1D array
    \param n  - number of elements of x

*/
double mean(const  double *x, int n )
{
    double s = 0;
    for ( int i = 0; i < n; ++i )
      s += x[i];

    return s/n;
}

/*!
    \fn double wmean( double *x, double *w, int n )

    \brief weighted mean of x

    \param x  An input 1D array.
    \param w  Weights.
    \param n  The number of elements of x and w.

*/
double wmean(const double *x, const double *w, int n )
{
    double s = 0;
    double sw = 0;
    for ( int i = 0; i < n; ++i )
    {
      s += x[i] * w[i];
      sw += w[i];
    }

    return s/sw;
}

/* cmp_double is defined in c_utils.c */
extern int cmp_double(const void *a, const void *b);

/*!
  \brief Returns quantiles of an array.

  \param x       An input 1D array.
  \param rn      A reference to the number of elements of x.
  \param probs   An array of probabilities at which the quantiles are to be computed.
  \param rnprobs A reference to the number of elements of probs.
  \param quants  The output array of quantiles of x.

  see https://www-users.york.ac.uk/~mb55/intro/quantile.htm
*/
void C_quantiles(const double *x, const int *rn, const double *probs, const int *rnprobs, double *quants)
{
    int n = rn[0];
    int nprobs = rnprobs[0];

    // Creating a copy of x so that the sorting of x does not change the order of elements of x
    double *y = (double*)malloc(n * sizeof(double));
    CHECK_PTR(y);

    for ( int i = 0; i < n; i++ )
      y[i] = x[i];

    qsort(y, n, sizeof(double), cmp_double);

    double q; // = probs[i] * n
    int j;    // the integer part of p
    for ( int i = 0; i < nprobs; i++ )
    {
      //q = probs[i] * (n + 1);
      q = probs[i] * (n - 1);
      j = (int)floor(q);
      quants[i] = y[j] + (y[j+1] - y[j]) * (q-j);
    }

    free(y);
}

/*!
    \brief Weighted mean performed on the columns of Y.

    \param Y        An array of transformed Y matrix.
    \param nrY      A reference to the number of rows of Y.
    \param ncY      A reference to the number of columns of Y.

    \param Tnn_i     An array of NN's indices corresponding to a t(nn.i) matrix in R.
    \param Tnn_w     An array of weights on K NN's corresponding to a t(nn.w) matrix in R.
    \param nrTnn     A reference to the number of rows of Tnn_ matrices.
    \param ncTnn     A reference to the number of columns of Tnn_ matrices.

    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param EYg       A output array of length weighted mean arrays, one for each row of Y.
*/
void C_matrix_wmeans(const double *Y,
                     const int    *nrY,
                     const int    *ncY,
                     const int    *Tnn_i,
                     const double *Tnn_w,
                     const int    *nrTnn,
                     const int    *ncTnn,
                     const int    *maxK,
                           double *EYg) {

    double *Tnn_y = (double*)malloc( (*nrTnn) * (*ncTnn) * sizeof(double) );
    CHECK_PTR(Tnn_y);

    int j, k;
    for ( int i = 0; i < (*ncY); i++ )
    {
      j = i * (*nrY);
      k = i * (*ncTnn);

      C_columnwise_eval(Tnn_i, nrTnn, ncTnn, Y + j, Tnn_y);
      C_columnwise_wmean(Tnn_y, Tnn_w, maxK, nrTnn, ncTnn, EYg + k);
    }

    free(Tnn_y);
}


/*!
    \brief Computes column-wise weighted means.
    
    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param Ey        An output array of column-wise weighted means.
*/
void C_columnwise_wmean(const double *nn_y,
                        const double *nn_w,
                        const int *maxK,
                        const int *rnr,
                        const int *rnc,
                        double *Ey) {
    int nr = rnr[0];
    int nc = rnc[0];
    
    int inr;  // holds i * nr
    int G;    // holds maxK[i] + 1
    double ys, s;
    
    for (int i = 0; i < nc; i++) {
        G = maxK[i] + 1;
        inr = i * nr;
        ys = 0;
        s = 0;
        
        for (int j = 0; j < G; j++) {
            ys += nn_y[j + inr] * nn_w[j + inr];
            s += nn_w[j + inr];
        }
        
        if (s > 0) {
            Ey[i] = ys / s;
        } else {
            Ey[i] = 0;
        }
    }
}

/*!
    \fn void C_columnwise_wmean_BB(double *nn_y, double *nn_w, int *maxK, int *rnr, int *rnc, int *rnBB, double *Ey)

    \brief Bayesian bootstrap of column-wise weighted means.

    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y; The sum of weights does not have to be 1 along columns.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param rnBB      A reference to the number of Bayesian bootstrap iterations.

    \param Ey        An output nc-by-nBB array of Bayesian bootstrap of column-wise weighted means.

    ALERT: The column-wise sum of weights must be > 0 !!! But they don't have to sum to 1 column-wise.
*/
void C_modified_columnwise_wmean_BB(const double *nn_y,
                                    const double *nn_w,
                                    const int    *maxK,
                                    const int    *rnr,
                                    const int    *rnc,
                                    const int    *rnBB,
                                          double *Ey) {
    int nr = rnr[0];
    int nc = rnc[0];
    int nBB = rnBB[0];

    // An array carrying Bayesian boostraped columns
    double *lambda = (double*)malloc( nr * sizeof(double));
    CHECK_PTR(lambda);

    int inr;  // holds i * nr
    int inc;  // holds iboot * nc
    int G;    // holds maxK[i] + 1
    double s;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      for ( int i = 0; i < nc; i++ )
      {
        G = maxK[i]+1;
        inr = i * nr;
        inc = iboot*nc;

        C_rsimplex(nn_w + inr, &G, lambda); // sampling from a Dirichlet distribution with the parameter w=(nn_w[0,i], ... , nn_w[nr-1, i]). The mean of that distribution is w/sum(w).

        s = 0;
        for ( int j = 0; j < G; j++ )
          s += nn_y[j + inr] * lambda[j];

        Ey[i + inc] = s;
      }
    }

    free(lambda);
}

/*!
    \fn void C_columnwise_wmean_BB(double *nn_y, double *nn_w, int *maxK, int *rnr, int *rnc, int *rnBB, double *Ey)

    \brief Bayesian bootstrap of column-wise weighted means.

    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y; The sum of weights does not have to be 1 along columns.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param rnBB      A reference to the number of Bayesian bootstrap iterations.

    \param Ey        An output nc-by-nBB array of Bayesian bootstrap of column-wise weighted means.

    ALERT: The column-wise sum of weights must be > 0 !!! But they don't have to sum to 1 column-wise.
*/
void C_columnwise_wmean_BB(const double *nn_y,
                           const double *nn_w,
                           const int    *maxK,
                           const int    *rnr,
                           const int    *rnc,
                           const int    *rnBB,
                                 double *Ey) {
    int nr = rnr[0];
    int nc = rnc[0];
    int nBB = rnBB[0];

    // An array carrying Bayesian boostraped columns
    double *lambda = (double*)malloc( nr * sizeof(double));
    CHECK_PTR(lambda);

    int inr;  // holds i * nr
    int inc;  // holds iboot * nc
    int G;    // holds maxK[i] + 1
    double w, ys, s;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      for ( int i = 0; i < nc; i++ )
      {
        G = maxK[i]+1;
        C_runif_simplex(&G, lambda);

        inr = i * nr;
        inc = iboot*nc;
        ys = 0;
        s = 0;

        for ( int j = 0; j < G; j++ )
        {
          w = nn_w[j + inr] * lambda[j];
          ys += nn_y[j + inr] * w;
          s  += w;
        }

        Ey[i + inc] = ys / s;
      }
    }

    free(lambda);
}

/*!
    \brief Bayesian bootstrap estimate of credible intervals (CI).

    \param rybinary  A reference to a binary (0/1 valueed) variable. If 1, predicted Ey will be trimmed to the unit interval \[0, 1\].
    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y; The sum of weights does not have to be 1 along columns.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param rnBB      A reference to the number of Bayesian bootstrap iterations.
    \param ralpha    A reference to the confidence level.

    \param Eyg_CI    An output array of BB estimates of Eyg CI's.

    ALERT: The column-wise sum of weights must be > 0 !!! But they don't have to sum to 1 column-wise.
*/
void C_columnwise_wmean_BB_qCrI(const int    *rybinary,
                                const double *nn_y,
                                const double *nn_w,
                                const int    *maxK,
                                const int    *rnr,
                                const int    *rnc,
                                const int    *rnBB,
                                const double *ralpha,
                                      double *Eyg_CI) {
    //
    // Generating BB Eyg's
    //
    int ybinary  = rybinary[0];
    int nc       = rnc[0];
    int nBB      = rnBB[0];
    double alpha = ralpha[0];

    // BB Ey estimate
    double *bbEy = (double*)malloc( nc * nBB * sizeof(double));
    CHECK_PTR(bbEy);

    C_columnwise_wmean_BB(nn_y, nn_w, maxK, rnr, rnc, rnBB, bbEy);

    if ( ybinary )
    {
      int n = nc * nBB;
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
    double p = alpha / 2;
    double probs[] = {p, 1-p};
    int two = 2;

    double *bb = (double*)malloc( nBB * sizeof(double)); // Holds BB Ey estimates of each row of bbEy
    CHECK_PTR(bb);

    for ( int i = 0; i < nc; i++ )
    {
      for ( int j = 0; j < nBB; j++ )
        bb[j] = bbEy[i + nc*j];

      C_quantiles(bb, rnBB, probs, &two, Eyg_CI + 2*i);
    }

    free(bbEy);
    free(bb);
}


/*!
    \fn void C_columnwise_wmean_BB_CI_1(double *Ey, double *nn_y, double *nn_w, int *maxK, int *rnr, int *rnc, int *rnBB, double *Ey_CI)

    \brief Bayesian bootstrap Median Absolute Deviations (MAD) of column-wise weighted means, centered at Ey.

    \param Ey        Ey estimate of length nc.
    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y; The sum of weights does not have to be 1 along columns.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param rnBB      A reference to the number of Bayesian bootstrap iterations.

    \param Ey_CI     An output length nc array of CI's of Bayesian bootstrap of column-wise weighted means.

    ALERT: The column-wise sum of weights must be > 0 !!! But they don't have to sum to 1 column-wise.
*/
void C_columnwise_wmean_BB_CrI_1(const double *Ey,
                                const double *nn_y,
                                const double *nn_w,
                                const int *maxK,
                                const int *rnr,
                                const int *rnc,
                                const int *rnBB,
                                      double *Ey_CI) {
    int nr = rnr[0];
    int nc = rnc[0];
    int nBB = rnBB[0];

    // BB of Ey
    double *Ey_BB = (double*)malloc( nc * nBB * sizeof(double));
    CHECK_PTR(Ey_BB);

    // An array carrying Bayesian boostraped columns
    double *lambda = (double*)malloc( nr * sizeof(double));
    CHECK_PTR(lambda);

    int inr;  // holds i * nr
    int inc;  // holds iboot * nc
    int G;    // holds maxK[i] + 1
    double w, ys, s;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      for ( int i = 0; i < nc; i++ )
      {
        G = maxK[i]+1;
        C_runif_simplex(&G, lambda);

        inr = i * nr;
        inc = iboot*nc;
        ys = 0;
        s = 0;

        for ( int j = 0; j < G; j++ )
        {
          w = nn_w[j + inr] * lambda[j];
          ys += nn_y[j + inr] * w;
          s  += w;
        }

        Ey_BB[i + inc] = ys / s;
      }
    }

    // This holds results of BB for a given index of nc
    double *bb = (double*)malloc( nBB * sizeof(double));
    CHECK_PTR(bb);

    // Row-wise Medians of Ey_BB
    double madC = 1.4826; //MAE scaling factor ensuring consistency between MAD
                          //and standard deviation when the data is normally
                          //distributed.
    madC *= 1.9;          // To get the boundary of credible interval (assuming normality).
    for ( int i = 0; i < nc; i++ )
    {
      for ( int iboot = 0; iboot < nBB; iboot++ )
      {
        bb[iboot] = fabs( Ey_BB[i + iboot*nc] - Ey[i] );
      }

      Ey_CI[i] = madC * median(bb, nBB);
    }

    free(lambda);
    free(Ey_BB);
    free(bb);
}


/*!
    \fn void C_columnwise_wmean_BB_CrI_2(double *nn_y, double *nn_w, int *maxK, int *rnr, int *rnc, int *rnBB, double *Ey_CI)

    \brief Another version of Bayesian bootstrap Median Absolute Deviations
    (MAD) of column-wise weighted means, centered at Ey. This version is much
    more memory efficient, but pays for it by calling runif_simplex() nBB*nc
    instead of nc times.

    \param Ey        Ey estimate of length nc.
    \param nn_y      An array of y values with nr rows and nc columns.
    \param nn_w      An array of weights of the same dimension as nn_y; The sum of weights does not have to be 1 along columns.
    \param maxK      An array of length nc with maxK[i] being the __index__ of the last non-zero weight in the i-th column of nn_w.
    \param rnr       A reference to the number of rows of nn_ arrays.
    \param rnc       A reference to the number of columns of nn_ arrays.
    \param rnBB      A reference to the number of Bayesian bootstrap iterations.

    \param Ey_CI     An output length nc array of CI's of Bayesian bootstrap of column-wise weighted means.

    ALERT: The column-wise sum of weights must be > 0 !!! But they don't have to sum to 1 column-wise.
*/
void C_columnwise_wmean_BB_CrI_2(const double *Ey,
                                const double *nn_y,
                                const double *nn_w,
                                const int    *maxK,
                                const int    *rnr,
                                const int    *rnc,
                                const int    *rnBB,
                                      double *Ey_CI) {
    int nr = rnr[0];
    int nc = rnc[0];
    int nBB = rnBB[0];

    // An array carrying Bayesian boostraped columns
    double *lambda = (double*)malloc( nr * sizeof(double));
    CHECK_PTR(lambda);

    // This holds results of BB for a given index of nc
    double *bb = (double*)malloc( nBB * sizeof(double));
    CHECK_PTR(bb);

    int inr;  // holds i * nr
    int G;    // holds maxK[i] + 1
    double w, ys, s;

    double madC = 1.4826; //MAE scaling factor ensuring consistency between MAD
                          //and standard deviation when the data is normally
                          //distributed.
    madC *= 1.9;          // To get the boundary of credible interval (assuming normality).

    for ( int i = 0; i < nc; i++ )
    {
      G = maxK[i]+1;

      for ( int iboot = 0; iboot < nBB; iboot++ )
      {
        C_runif_simplex(&G, lambda);

        inr = i * nr;
        ys = 0;
        s = 0;

        for ( int j = 0; j < G; j++ )
        {
          w = nn_w[j + inr] * lambda[j];
          ys += nn_y[j + inr] * w;
          s  += w;
        }

        bb[iboot] = fabs( (ys / s) - Ey[i] );
      }

      Ey_CI[i] = madC * median(bb, nBB);
    }

    free(lambda);
    free(bb);
}

/*!
    \fn void C_columnwise_eval(int *nn_i, double *x, int *rK, int *rng, int *rnx, double *nn_x)

    \brief Evaluates x over an K-by-ng array nn_i of indices of get.knnx(x, xg, k)

    \param nn_i      An array of indices of K nearest neighbors of the i-th element of x.
    \param rK        A reference to the number of rows of nn_i.
    \param rng       A reference to the number of columns of nn_i.
    \param x         An array of nx elements.

    \param nn_x      A return array of the same dimensions as nn_i with values of x at the indices of nn_i.

*/
void C_columnwise_eval(const int    *nn_i,
                       const int    *rK,
                       const int    *rng,
                       const double *x,
                             double *nn_x) {
    int K  = rK[0];
    int ng = rng[0];

    int iK;   // holds i * K
    for ( int i = 0; i < ng; i++ )
    {
      iK = i * K;
      for ( int j = 0; j < K; j++ )
        nn_x[j + iK] = x[nn_i[j + iK]];
    }
}

/*!
    \fn void C_columnwise_eval(int *nn_i, double *x, int *rK, int *rng, int *rnx, double *nn_x)

    \brief Evaluates x over an K-by-ng array nn_i of indices of get.knnx(x, xg, k)

    \param TX        An array corresponding to the transpose of a numeric matrix, X, in R. The output matrix is the transpose of X / y (row-wise). The output overwrites TX.
    \param rnrTX     A reference to the number of rows of TX.
    \param rncTX     A reference to the number of columns of TX.
    \param y         An array of ncTX elements.
*/
void C_mat_columnwise_divide(      double *TX,
                             const int    *rnrTX,
                             const int    *rncTX,
                             const double *y) {
    int nrTX = rnrTX[0];
    int ncTX = rncTX[0];

    int inrTX;   // holds i * nrTX
    for ( int i = 0; i < ncTX; i++ )
    {
      if ( i % 10000 )
        Rprintf("\r%d", i);

      inrTX = i * nrTX;
      for ( int j = 0; j < nrTX; j++ )
        TX[j + inrTX] /= y[i];
    }

    Rprintf("\r");
}


/*
  \brief Normalizes a distance matrix

  normalize.dist(x, min.K, bw) divides the columns of the input distance matrix,
  x, by 'bw' if the distance of the min.K-st element of the row is <bw, if it is
  not, then all elements are divided by the distance of the min.K-st element.
  The radius, r, is bw in the first case and x[minK + ir] in the second.

  \param d      A numeric 1D array corresponding to a nn.d matrix in R.
  \param rnr    A reference to the number of rows of the matrix associated with x.
  \param rnc    A reference to the number of columns of the matrix associated with x.
  \param rminK     A reference to the mininum __number_ of elements in the rows of nn_i/nn_d that need to have weights > 0.
  \param rbw    A reference to a non-negative normalization constant.

  \param nd     An output normalized distance array.
  \param r      An output radius of support of the given model array.

  ALERT: This routine makes sense only in the situation when nr > minK. Thus,
  the user has to make sure this condition is satisfied before calling
  normalize_dist().

*/
void C_normalize_dist(const double *d,
                      const int    *rnr,
                      const int    *rnc,
                      const int    *rminK,
                      const double *rbw,
                            double *nd,
                            double *r) {
   int nr    = rnr[0];
   int nc    = rnc[0];
   int minK  = rminK[0];
   double bw = rbw[0];

   if ( nr - 1 < minK ) // In theory, we would need nr >= minK, but to simplify
                        // the way we enforce the minimal number of rows of nn_d
                        // to have weights > 0, requires that nr > minK, so if
                        // nr <= minK <=> nr - 1 < minK, we need to throw an
                        // error. See the big comment below.
   {
     error("ERROR in normalize_dist(); file %s  at line %d;  nr - 1 < minK !", __FILE__, __LINE__);
   }

   int ir;  // holds i * nr
   double z;
   for ( int i = 0; i < nc; i++ )
   {
     ir = i * nr;
     z = d[minK + ir]; // Formally, we should have written this as z =
                       // x[iminK + 1 + ir], where iminK = minK - 1
                       // is the __index__ of the minK-th element in the
                       // row of x. If z < bw, then x[iminK + 1 + ir] / bw
                       // < 1 and so the weight of the (iminK + 1)-th
                       // element will be > 0.
     if ( z < bw )
     {
       r[i] = bw;
       for ( int j = 0; j < nr; j++ )
         nd[j + ir] = d[j + ir] / bw;
     } else {
       z = d[minK + ir]; // The same comment as above !
       r[i] = z;
       for ( int j = 0; j < nr; j++ )
         nd[j + ir] = d[j + ir] / z;
     }
   }
}


/*
  \brief Normalizes a distance matrix

  A version of normalize_dist() whose minK parameter is an array of length nc.

  \param x      A numeric 1D array corresponding to an t(nn.d) matrix in R.
  \param rnr    A reference to the number of rows of the matrix associated with x.
  \param rnc    A reference to the number of columns of the matrix associated with x.
  \param minK   An array of length nr, such that, minK[i] is the mininum __number__ of elements in each row of x with weights > 0.
  \param rbw    A reference to a non-negative normalization constant.
  \param y      The output normalized array.

  ALERT: minK[j] has to be always < nr !!!

  This routine is used by cv_llm_2D()

*/
void C_normalize_dist_with_minK_a(const double *x,
                                  const int    *rnr,
                                  const int    *rnc,
                                  const int    *minK,
                                  const double *rbw,
                                        double *y) {
   int nr    = rnr[0];
   int nc    = rnc[0];
   double bw = rbw[0];

   int ir;  // holds i * nr
   double z;
   for ( int i = 0; i < nc; i++ )
   {
     ir = i * nr;
     z = x[minK[i] + ir]; // Formally, we should have written this as z =
                          // x[iminK[i] + 1 + ir], where iminK[i] = minK[i] - 1
                          // is the __index__ of the minK[i]-th element in the
                          // row of x. If z < bw, then x[iminK[i] + 1 + ir] / bw
                          // < 1 and so the weight of the (iminK[i] + 1)-th
                          // element will be > 0.
     if ( z < bw )
     {
       for ( int j = 0; j < nr; j++ )
         y[j + ir] = x[j + ir] / bw;
     } else {
       z = x[minK[i] + ir]; // The same comment as above !
       for ( int j = 0; j < nr; j++ )
         y[j + ir] = x[j + ir] / z;
     }
   }
}


/*
  \brief Samples with replacement from 0, ... ,(n-1)

  \param rn      A reference to the length of the sequence of indices the samples come from.

  \param mult    A pointer to an int array of length n with multiplicities of indices
                 selected. An index has multiplicity m if it was selected m times.
*/
void C_samplewr( const int *rn, int *mult) {
    int n = rn[0];

    for ( int i = 0; i < n; i++ )
      mult[i] = 0;

    GetRNGstate();
    double max = n - 0.01;
    int s; // selection index
    for ( int i = 0; i < n; i++ )
    {
      s = (int)runif(0.0, max);
      mult[s]++;
    }
    PutRNGstate();
}

/*!

  \brief Performs in place permutation of an int array using Fisher-Yates
  shuffle as implemented by Durstenfeld.

  \param x An int array.
  \param n The size of the array.

*/
void C_permute(int *x, int n)
{
  int k;
  int tmp;

  GetRNGstate();

  while ( n > 1 )
  {
    //k = rand() % n; // 0 <= k < n.
    k = (int)(n * runif(0.0,1.0));// 0 <= k < n.
    n--;            // n is now the last pertinent index;
    tmp  = x[n];    // swap x[n] with x[k] (does nothing if k == n).
    x[n] = x[k];
    x[k] = tmp;
  }

  PutRNGstate();
}

/*!

  \brief Performs in place permutation of a double array using Fisher-Yates
  shuffle as implemented by Durstenfeld.

  \param x An int array.
  \param n The size of the array.

*/
void C_dpermute(double *x, int n)
{
  int k;
  int tmp;

  GetRNGstate();

  while ( n > 1 )
  {
    //k = rand() % n; // 0 <= k < n.
    k = (int)(n * runif(0.0,1.0));// 0 <= k < n.
    n--;            // n is now the last pertinent index;
    tmp  = x[n];    // swap x[n] with x[k] (does nothing if k == n).
    x[n] = x[k];
    x[k] = tmp;
  }

  PutRNGstate();
}

/*!

  \brief Performs in place permutation of an int array using Fisher-Yates
  shuffle as implemented by Durstenfeld.

  \param x  An int array of indices.
  \param rn A reference to the size of the array.
  \param y  An array of permuted indices.

*/
void C_vpermute(int *x, int *rn, int *y)
{
    int n = rn[0];
    int k;
    int tmp;

    for ( int i = 0; i < n; i++ )
      y[i] = x[i];

    GetRNGstate();

    while ( n > 1 )
    {
      k = (int)(n * runif(0.0,1.0));// 0 <= k < n.
      n--;            // n is now the last pertinent index;
      tmp  = y[n];    // swap y[n] with y[k] (does nothing if k == n).
      y[n] = y[k];
      y[k] = tmp;
    }

    PutRNGstate();
}

/*!
    \brief Generates a split of integers from 0 to (n-1) into 'nfolds' folds of size ~ n/nfolds

    \param n      The number of integers from 0 to (n-1) that will be split into
                    'nfolds' groups of more or less equal size.

    \param nfolds The number of folds, that is, the number of groups of more or
                    less equal size into which the set of integers is going to be split.

    \return Returns array_t of size nfolds.


    ALERT: The user is responsible for freeing memory from the returned folds array_t !!!


    NOTE: in the future add double *y parameter, which of length nrow(X), so
    that the distribution of y on X and on X[-I] is more or less the same.
*/
iarray_t * get_folds(int n, int nfolds) {
    // Creating an array of consecutive integers 0, 1, ... , (n-1)
    int *x = (int*)malloc(n * sizeof(int));

    for ( int i = 0; i < n; i++ )
      x[i] = i;

    // Permuting it in-place
    C_permute(x, n);

    // Determining sizes of folds
    int fold_size = n / nfolds;
    int r = n % nfolds;

    int *fold_sizes = (int*)malloc(nfolds * sizeof(int));
    for ( int i = 0; i < nfolds; i++ )
      fold_sizes[i] = fold_size;

    int k = 0;
    while ( r > 0 && k < nfolds )
    {
      fold_sizes[k] += 1;
      r--;
      k++;
    }

    // Creating folds iarray_t
    iarray_t *folds = (iarray_t*)malloc( nfolds * sizeof(iarray_t));

    k = 0;
    for ( int i = 0; i < nfolds; i++ )
    {
      folds[i].x = (int*)malloc( fold_sizes[i] * sizeof(int) );
      folds[i].size = fold_sizes[i];

      for ( int j = 0; j < fold_sizes[i]; j++ )
        folds[i].x[j] = x[k++];
    }

    free(x);
    free(fold_sizes);

    return folds;
}


/*!
    \brief Generates a split of integers from 0 to (n-1) into 'nfolds' folds of size ~ n/nfolds

    \param rn      A reference to the number of integers from 0 to (n-1) that will be split into
                    'nfolds' groups of more or less equal size.

    \param rnfolds A reference to the number of folds, that is, the number of groups of more or
                    less equal size into which the set of integers is going to be split.

    \param folds   An output array of length n ( with memory allocated in the R interface to this function) holding folds.

*/
void C_v_get_folds(const int *rn,
                 const int *rnfolds,
                       int *folds)
{
    int n      = rn[0];
    int nfolds = rnfolds[0];

    iarray_t *F = get_folds(n, nfolds);

    for ( int ifold = 0; ifold < nfolds; ifold++ )
    {
      int *fold = F[ifold].x;    // indices of the fold with values between 0 and nrX-1
      int nfold = F[ifold].size; // number of the elements in fold

      for ( int i = 0; i < nfold; i++ )
        folds[fold[i]] = ifold;
    }

    // freeing memory allocated to folds[i].x and folds from within get_folds(nrX, nfolds)
    for ( int i = 0; i < nfolds; i++ )
      free(F[i].x);
    free(F);
}

/* C_quantiles is already defined earlier in the file at line 146 */

/*
  Winsorizes a double array

  \param y   A numeric vector of length n.
  \param rn  A reference to the length of y.
  \param rp  A reference to the proportion of the bottom and  top 100*p percent of y that will be replace by p and the 1-p quantile, respectively.
  \param wy  The output array of winsorized y values.
*/
void C_winsorize(const double *y, const int *rn, const double *rp, double *wy)
{
    int n = rn[0];
    double p = rp[0];

    double probs[2] = {p, 1-p};
    int nprobs = 2;
    double quants[2];

    C_quantiles(y, &n, probs, &nprobs, quants);

    for ( int i = 0; i < n; i++ )
      if ( y[i] < quants[0] ) wy[i] = quants[0];
      else if ( y[i] > quants[1] ) wy[i] = quants[1];
      else wy[i] = y[i];
}

/*!
    \brief Returns an estimate of the probability that the samples from the
    density function from which x is sampled have values larger than x0.

    \param x     An array of double values.
    \param rnx   A reference to the number of elements of x.
    \param rz    A reference to the selected value, z, for which we want to
                   know the proportion of elments of x greater than z.

    \param p     The proportion of x values greater than z.
*/
void C_pdistr(const double *x,
              const int    *rnx,
              const double *rz,
                    double *p)
{
    int nx = rnx[0];
    double z = rz[0];

    // Creating a copy of x so that the sorting of x does not change the order of elements of x
    double *y = (double*)malloc(nx * sizeof(double));
    CHECK_PTR(y);

    for ( int i = 0; i < nx; i++ )
      y[i] = x[i];

    qsort(y, nx, sizeof(double), cmp_double);

    int i = 0;    // the integer part of p
    while ( y[i] < z && i < nx ) {
      i++;
    }

    free(y);

    *p = (double)(nx - i) / nx;
}


/*!
    \brief Returns an estimate of the probability that the samples from the
    density function from which x is sampled have values larger than x0.

    It is an SEXP version of the above function.

    \param sx   A vector of real values.
    \param sz   A real value, z, for which we want to know the proportion of elments of x greater than z.

    \return     The proportion of x values greater than z.
*/
SEXP S_pdistr( SEXP sx, SEXP sz)

{
    int nx   = length(sx);
    double z = asReal(sz);
    double *x = REAL(sx);

    // Creating a copy of x so that the sorting of x does not change the order of elements of x
    double *y = (double*)malloc(nx * sizeof(double));
    CHECK_PTR(y);

    for ( int i = 0; i < nx; i++ )
      y[i] = x[i];

    qsort(y, nx, sizeof(double), cmp_double);

    int i = 0;    // the integer part of p
    while ( y[i] < z && i < nx ) {
      i++;
    }

    free(y);

    // Creating output variable
    SEXP out = PROTECT(allocVector(REALSXP, 1));
    double *pout = REAL(out);

    *pout = (double)(nx - i) / nx;

    UNPROTECT(1);

    return out;
}


/*!
   Pearson correlation coefficient

   \param x   An array of double values of length n.
   \param y   An array of double values of length n.
   \param rn  A reference to the length of x and y.
   \param rc  A reference to the value of the correlation coefficient between x and y.
*/
void C_pearson_cor(const double *x, const double *y, const int *rn, double *rc)
{
    int n = rn[0];

    double sxy =0, sx =0, sy = 0, s2x =0, s2y = 0;
    for ( int i = 0; i < n; i++ )
    {
      sxy += x[i] * y[i];
      sx  += x[i];
      sy  += y[i];
      s2x += x[i]*x[i];
      s2y += y[i]*y[i];
    }

    double c = n * sxy - sx * sy;
    c /= sqrt(n*s2x - sx*sx) * sqrt(n*s2y - sy*sy);

    *rc = c;
}


/*!
   Weighted covariance

   \param x    An array of double values of length n.
   \param y    An array of double values of length n.
   \param w    An array of non-negative weights of length n.
   \param rn   A reference to the length of x, y and w.
   \param rwc  A reference to the weighted covariance between x and y.

   NOTE: The weights do not have to sum up to 1.
*/
void C_wcov(const double *x, const double *y, const double *w, const int *rn, double *rwc)
{
    int n = rn[0];

    double wx = wmean(x, w, n);
    double wy = wmean(y, w, n);

    double wc = 0, sw = 0;
    for ( int i = 0; i < n; i++ )
    {
      wc += w[i] * (x[i] - wx) * (y[i] - wy);
      sw += w[i];
    }

    *rwc = wc / sw;
}

/*!
   Pearson weighted correlation

   \param x    An array of double values of length n.
   \param y    An array of double values of length n.
   \param w    An array of non-negative weights of length n.
   \param rn   A reference to the length of x, y and w.
   \param rwc  A reference to the weighted correlation between x and y.

   NOTE: The weights do not have to sum up to 1.
*/
void C_pearson_wcor(const double *x, const double *y, const double *w, const int *rn, double *rwc)
{
    C_wcov(x, y, w, rn, rwc);

    double wx = 0, wy = 0;
    C_wcov(x, x, w, rn, &wx);
    C_wcov(y, y, w, rn, &wy);

    if ( wx > 0 && wy > 0 ){
      *rwc /= sqrt(wx * wy);
    } else {
      *rwc = 0;
    }
}

/*!
   Bayesian bootstrap credible intervals estimates of Pearson weighted correlations

   \param nn_y1   The transpose of the matrix of y1 values of nn_i.
   \param nn_y2   The transpose of the matrix of y2 values of nn_i.
   \param nn_i    A matrix of NN indices.
   \param nn_w    A matrix of NN weights.
   \param rK      A reference to the number rows of nn_ arrays.
   \param rng     A reference to the number columns of nn_ arrays.
   \param rnx     A reference to the number of elements of x over which y1 and y2 are defined.
   \param rnBB    A reference to the number of Bayesian bootstrap iterations.
   \param ralpha  A reference to the level of significance.
   \param qCI     An array of BB credible intervals.

   NOTE: The weights do not have to sum up to 1.
*/
void C_pearson_wcor_BB_qCrI(const double *nn_y1,
                            const double *nn_y2,
                            const int    *nn_i,
                            const double *nn_w,
                            const int    *rK,
                            const int    *rng,
                            const int    *rnx,
                            const int    *rnBB,
                            const double *ralpha,
                                  double *qCI)
{
    int K   = rK[0];
    int ng  = rng[0];
    int nx  = rnx[0];
    int nBB = rnBB[0];

    // BB wcor estimate
    double *bbwcor = (double*)malloc( ng * nBB * sizeof(double));
    CHECK_PTR(bbwcor);

    // An array of Bayesian boostrap weights
    double *lambda = (double*)malloc( nx * sizeof(double));
    CHECK_PTR(lambda);

    // An array of modified by Bayesian boostrap weights
    double *w = (double*)malloc( K * sizeof(double));
    CHECK_PTR(w);

    int jK, ing;
    for ( int iboot = 0; iboot < nBB; iboot++ )
    {
      C_runif_simplex(rnx, lambda);

      ing = iboot * ng;
      for ( int j = 0; j < ng; j++ ) // bootstraping the j-th columns of nn_y1, nn_y2, nn_w
      {
        jK = j * K;

        for ( int i = 0; i < K; i++ )
          w[i] = nn_w[i + jK] * lambda[nn_i[i + jK]];

        C_pearson_wcor(nn_y1 + jK, nn_y2 + jK, w, rK, bbwcor + ing + j);
      }

    } // END OF for iboot

    //
    // Computing (p, 1-p) quantiles of of each row of bbwcor
    //
    double p = *ralpha / 2;
    double probs[] = {p, 1-p};
    int two = 2;

    double *bb = (double*)malloc( nBB * sizeof(double)); // Holds BB wcor() estimates of each row of bbwcor
    CHECK_PTR(bb);

    for ( int j = 0; j < ng; j++ )
    {
      for ( int iboot = 0; iboot < nBB; iboot++ )
        bb[iboot] = bbwcor[j + iboot*ng];

      C_quantiles(bb, rnBB, probs, &two, qCI + 2*j);
    }

    free(bb);
    free(w);
    free(bbwcor);
    free(lambda);
}


/*!
   Estimates local correlation on X.grid

   \param nn_i  A matrix of indicies of X with ng rows and K columns.
   \param nn_w  A matrix of weights with ng rows and K columns.
   \param Y     A matrix with nrX rows.

*/
SEXP S_lwcor( SEXP Snn_i, SEXP Snn_w , SEXP SY)
{
    #define DEBUG_S_LWCOR 0

    int nprot = 0;
    PROTECT( Snn_i = AS_NUMERIC(Snn_i)); nprot++;
    PROTECT( Snn_w = AS_NUMERIC(Snn_w)); nprot++;
    PROTECT( SY = AS_NUMERIC(SY)); nprot++;

    double *nn_i = REAL(Snn_i);
    double *nn_w = REAL(Snn_w);
    double *Y = REAL(SY);

    int *dim_nn = INTEGER( GET_DIM(Snn_i) );
    int *dimY   = INTEGER( GET_DIM(SY) );

    int ng = dim_nn[0];
    int K  = dim_nn[1];
    int nrY = dimY[0];
    int ncY = dimY[1];
    int npairs = ncY * (ncY - 1) / 2;

    #if DEBUG_S_LWCOR
    Rprintf("S_lwcor(): ng=%d  K=%d  nrY=%d  ncY=%d  npairs=%d\n", ng, K, nrY, ncY, npairs);
    #endif

    SEXP Sans; /* Create SEXP to hold result */
    PROTECT( Sans = allocMatrix(REALSXP, ng, npairs)); nprot++;
    double *ans = REAL( Sans );

    SEXP Sw; /* Create SEXP to hold weights */
    PROTECT( Sw = allocVector(REALSXP, K)); nprot++;
    double *w = REAL( Sw );

    SEXP Sy1;
    PROTECT( Sy1 = allocVector(REALSXP, nrY)); nprot++;
    double *y1 = REAL( Sy1 );

    SEXP Sy2;
    PROTECT( Sy2 = allocVector(REALSXP, nrY)); nprot++;
    double *y2 = REAL( Sy2 );

    double *yi, *yj;
    int jKng;
    int ipair = 0;
    int ii;
    int k;
    int ncY1 = ncY - 1;
    for ( int i = 0; i < ncY1; i++ )
    {
      for ( int j = i + 1; j < ncY; j++ )
      {
        for ( int ig = 0; ig < ng; ig++ )
        {
          yi = Y + i*nrY; // i-th column of Y
          yj = Y + j*nrY; // j-th column of Y
          k = 0;
          for ( int jK = 0; jK < K; jK++ )
          {
            jKng = jK * ng;
            ii = nn_i[ig + jKng];
            if ( R_FINITE(yi[ii]) && R_FINITE(yj[ii]) )
            {
              y1[k] = yi[ii]; // Y[nn_i[ig,],i]
              y2[k] = yj[ii]; // Y[nn_i[ig,],j]
              w[k] = nn_w[ig + jKng];
              k++;
            }
          }

          C_pearson_wcor(y1, y2, w, &k, ans + ig + ipair*ng);

          #if DEBUG_S_LWCOR
          Rprintf("i=%d j=%d ig=%d ipair=%d lcor=%.5f\n", i, j, ig, ipair, *(ans + ig + ipair*ng));

          char *file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/y1.csv";
          write_double_array(y1, k, file);
          Rprintf("S_lwcor() created %s\n", file);

          file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/y2.csv";
          write_double_array(y2, k, file);
          Rprintf("S_lwcor() created %s\n", file);

          file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/w.csv";
          write_double_array(w, k, file);
          Rprintf("S_lwcor() created %s\n", file);

          error("STOP");
          #endif

        } // END of for ig
        ipair++;
      }
    }

    UNPROTECT( nprot ); /* Wrap up; */

    return( Sans );
}


/*!
   Estimates local correlation between y and the columns of Y over X

   \param Snn_i  A matrix of indicies of X with ng rows and K columns.
   \param Snn_w  A matrix of weights with ng rows and K columns.
   \param SY     A matrix with nrX rows.
   \param Sy     An array with length(y)=nrow(X), such the the local

*/
SEXP S_lwcor_yY( SEXP Snn_i, SEXP Snn_w , SEXP SY, SEXP Sy )
{
    int nprot = 0;
    PROTECT( Snn_i = AS_NUMERIC(Snn_i)); nprot++;
    PROTECT( Snn_w = AS_NUMERIC(Snn_w)); nprot++;
    PROTECT( SY = AS_NUMERIC(SY)); nprot++;
    PROTECT( Sy = AS_NUMERIC(Sy)); nprot++;

    double *nn_i = REAL(Snn_i);
    double *nn_w = REAL(Snn_w);
    double *Y = REAL(SY);
    double *y = REAL(Sy);

    int *dim_nn = INTEGER( GET_DIM(Snn_i) );
    int *dimY   = INTEGER( GET_DIM(SY) );

    int ng = dim_nn[0];
    int K  = dim_nn[1];
    int nrY = dimY[0];
    int ncY = dimY[1];

    SEXP Sres; /* Create SEXP to hold result */
    PROTECT( Sres = allocMatrix(REALSXP, ng, ncY) ); nprot++;
    double *res = REAL( Sres );

    SEXP Sw; /* Create SEXP to hold weights */
    PROTECT( Sw = allocVector(REALSXP, K)); nprot++;
    double *w = REAL( Sw );

    SEXP Sy1;
    PROTECT( Sy1 = allocVector(REALSXP, nrY)); nprot++;
    double *y1 = REAL( Sy1 );

    SEXP Sy2;
    PROTECT( Sy2 = allocVector(REALSXP, nrY)); nprot++;
    double *y2 = REAL( Sy2 );

    double *yi;
    int jKng;
    int ii;
    int k;
    for ( int i = 0; i < ncY; i++ )
    {
      for ( int ig = 0; ig < ng; ig++ )
      {
        yi = Y + i*nrY; // i-th column of Y
        k = 0;
        for ( int jK = 0; jK < K; jK++ )
        {
          jKng = jK * ng;
          ii = nn_i[ig + jKng];
          if ( R_FINITE(yi[ii]) && R_FINITE(y[ii]) )
          {
            y1[k] = yi[ii]; // Y[nn_i[ig,],i]
            y2[k] = y[ii];  // y[nn_i[ig,]]
            w[k] = nn_w[ig + jKng];
            k++;
          }
        }

        C_pearson_wcor(y1, y2, w, &k, res + ig + i*ng);

        #if DEBUG_S_LWCOR_YY
        Rprintf("i=%d j=%d ig=%d lcor=%.5f\n", i, j, ig, *(res + ig + i*ng));

        char *file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/y1.csv";
        write_double_array(y1, k, file);
        Rprintf("S_lwcor() created %s\n", file);

        file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/y2.csv";
        write_double_array(y2, k, file);
        Rprintf("S_lwcor() created %s\n", file);

        file = "/Users/pgajer/projects/rllm/data-debugging/S_lwcor/w.csv";
        write_double_array(w, k, file);
        Rprintf("S_lwcor() created %s\n", file);

        error("STOP");
        #endif

      } // END OF for ig
    } // END OF for ( int i = 0; i < ncY; i++ )

    UNPROTECT( nprot ); /* Wrap up; */

    return( Sres );
}

/*!
  Computes density-adjusted distance between a set of points in n-dimensional space.

  \param X          The set of points in x-dim space for which the distance matrix is to be computed.
                    It comes from an n-by-x matrix in R that was transposed before
                    passing to this routine. Thus, it is now an x-by-n matrix passed to C in a
                    column-major format.
  \param rdim       A reference to the numer of columns of X (before the transpose operation).
  \param rnrX       A reference of the number of rows of X (before the transpose operation).
  \param density    An double array holding values of a density function. All its values need to be non-negative. This assumption can be weakened by assuming for any pair of points (p,q) the sum density(p) + density(q) is greater than zero.

  \param dist       A symmetric n-by-n array of distances between the points of X. It is assumed that the memory to dist has been allocated in the R interface function.
 */
void C_density_distance(const double *X,
                        const int    *rdim,
                        const int    *rnrX,
                        const double *density,
                        double       *dist)
{
    int dim = *rdim; // the dimensionality of the set of points
    int n   = *rnrX; // the number of points
    int i, j, k;
    double euclidean_dist, adjusted_dist;

    // Loop through each pair of points
    for (i = 0; i < n; ++i) {
        for (j = i; j < n; ++j) {
            euclidean_dist = 0.0;

            // Calculate Euclidean distance
            for (k = 0; k < dim; ++k) {
                euclidean_dist += pow(X[i * dim + k] - X[j * dim + k], 2);
            }
            euclidean_dist = sqrt(euclidean_dist);

            // Adjust the distance using density values
            double mean_density = (density[i] + density[j]) / 2.0;
            adjusted_dist = euclidean_dist / mean_density;

            // Fill in the symmetric distance matrix
            dist[i * n + j] = adjusted_dist;
            dist[j * n + i] = adjusted_dist;
        }
    }
}

/**
 * @brief Generates random query points based on the range of values in X.
 *
 * For each dimension in X, this function identifies the range (minimum to maximum)
 * of values and generates random query points within this range. The resulting points
 * are stored in the Q array.
 *
 * @param X    Pointer to the original data matrix. This matrix is represented as
 *             a one-dimensional double array in column-major order.
 * @param nX   Number of columns (points) in X.
 * @param dim  Dimensionality of the data (number of rows in X).
 * @param Q    Pointer to an array where the generated query points will be stored.
 *             The array should be pre-allocated with sufficient space to hold
 *             `dim * nQ` doubles.
 * @param nQ   Number of query points to be generated.
 *
 * @note The Q array should be pre-allocated before calling this function.
 * After the function call, Q will hold the generated query points, stored
 * in column-major order.
 */
void generate_Q(const double* X, int nX, int dim, double* Q, int nQ) {

    GetRNGstate();

    for (int d = 0; d < dim; d++) {
      // Find the minimum and maximum values in the d-th dimension of X
      double min_val = X[d];
      double max_val = X[d];
      for (int i = 1; i < nX; i++) {
        if (X[d + i*dim] < min_val) min_val = X[d + i*dim];
        if (X[d + i*dim] > max_val) max_val = X[d + i*dim];
      }

      // Generate the d-th coordinate values for Q points
      for (int i = 0; i < nQ; i++) {
        Q[d + i*dim] = runif(min_val, max_val); //rand_range(min_val, max_val);
      }
    }
}


/*!
  Generates a random matrix based on the range of values in X.

  For each dimension in X, this function identifies the range (minimum to maximum)
  of values and generates random query points within this range. The resulting points
  are stored in the Q array.

  \param X    A numerica matrix.
  \param nX   The number of rows in X.
  \param dim  The number of columns of X.
  \param Q    A random matrix to be produced from X.
  \param nQ   The rows of Q.

 */
void C_rmatrix(const double* X,
               const int *rnX,
               const int *rdim,
               double* Q,
               const int *rnQ) {
    generate_Q(X, *rnX, *rdim, Q, *rnQ);
}
