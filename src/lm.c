#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <stdlib.h>

#include "msr2.h"

/*!
    \fn void C_fsolve(double *a, double *b, int *nr, int *nc, int *status)

    \brief Solves a x = b with the solution x clobbering b (LAPACK).

    Solves a x = b with the solution x clobbering b (LAPACK) first call DPOTRF and then DPOTRS.

    \param a      An input 1D array corresponding to a symmetric positive defined matrix in column-major format (from R).
    \param b      An input 1D array in column-major format (from R).
    \param nr     The number of rows of a and b.
    \param nc     The number of columns of b.
    \param status The binary variable set to 0 if there are no errors (like EXIT_SUCCESS) and 1 if there were errors.

*/
void C_fsolve(double *a, double *b, int *nr, int *nc, int *status) {
    int info;
    *status = 0;
    F77_CALL(dpotrf)("L", nr, a, nr, &info FCONE);
    if (info != 0) {
      *status = 1;
      //Rprintf("C_fsolve(): Cholesky decomposition failed\n");
    } else {
      F77_CALL(dpotrs)("L", nr, nc, a, nr, b, nr, &info FCONE);
      if (info != 0) {
        *status = 1;
        //error("F77_CALL(dpotrs) crashed!");
      }
    }
}


/*!
    \fn void C_flm_1D(double *x, double *y, int *n, double *beta)

    \brief 1D linear model.

    One dimensional linean model using fsolve().

    \param x      A predictor 1D numeric array.
    \param y      An outcome variable 1D numeric array.
    \param n      The number of elements of x and y arrays.
    \param beta   An array of length 2 with the intercept and slop of the model

*/
void C_flm_1D(const double *x,
              const double *y,
              const int    *rn,
                    double *beta) {
    int n = rn[0];

    // A cross-product X^t * X, where X = cbind(1, x)
    double XtX[4];
    XtX[0] = n;

    double s  = 0;
    double s2 = 0;
    for ( int i = 0; i < n; i++ )
    {
      s += x[i];
      s2 += x[i]*x[i];
    }

    XtX[1] = s;
    XtX[2] = s;
    XtX[3] = s2;

    // A cross-product X^t * y, where X = cbind(1, x)
    s  = 0;
    s2 = 0;
    for ( int i = 0; i < n; i++ )
    {
      s += y[i];
      s2 += x[i]*y[i];
    }

    beta[0] = s;
    beta[1] = s2;

    // solving A*beta = B
    int nr = 2;
    int nc = 1;
    int status = 0;
    C_fsolve(XtX, beta, &nr, &nc, &status);
}


/*!
    \fn void C_flmw_1D(double *x, double *y, double *w, int *n, double *beta, int *status)

    \brief 1D linear model.

    One dimensional linean model using fsolve().

    \param x      A predictor 1D numeric array.
    \param y      An outcome variable 1D numeric array.
    \param w      Weigths; all must be non-negative, and they have to sum to 1!
    \param rn     The reference to the number of elements of x, y and w arrays.
    \param beta   An array of length 2 with the intercept and slop of the model
    \param status The binary variable set to 0 if there are no errors (like EXIT_SUCCESS) and 1 if there were errors.

*/
void C_flmw_1D(const double *x,
               const double *y,
               const double *w,
               const int    *rn,
                     double *beta,
                     int    *status) {
    int n = rn[0];

    // A cross-product X^t * W * X, where X = cbind(1, x) and W = diag(w)
    double XtWX[4];

    double sw   = 0; // sum of w's
    double swx  = 0; // sum of w[i]*x[i]'s
    double swx2 = 0; // sum of w[i]*x[i]*x[i]
    for ( int i = 0; i < n; i++ )
    {
      sw   += w[i];
      swx  += w[i]*x[i];
      swx2 += w[i]*x[i]*x[i];
    }

    XtWX[0] = sw;
    XtWX[1] = swx;
    XtWX[2] = swx;
    XtWX[3] = swx2;

    // A cross-product X^t * W * y, where X = cbind(1, x)
    double swy  = 0;
    double swxy = 0;
    for ( int i = 0; i < n; i++ )
    {
      swy  += w[i]*y[i];
      swxy += w[i]*x[i]*y[i];
    }

    beta[0] = swy;
    beta[1] = swxy;

    // solving A*beta = B
    int nr = 2;
    int nc = 1;
    *status = 0;
    C_fsolve(XtWX, beta, &nr, &nc, status);
}


/*!
    \fn void C_flm(double *x, double *y, int *rnr, int *rnc, double *beta)

    \brief General ordinary least squares linear model.

    General ordinary least squares linean model using fsolve().

    \param x    - predictor numeric 1D array derived from an n-by-p matrix.
    \param y    - outcome variable - 1D numeric array.
    \param rnr  - reference to the number of rows of the matrix associated with x.
    \param rnc  - reference to the number of columns of the matrix associated with x.
    \param beta - array of length p+1 with the coefficients of the model.

*/
void C_flm(const double *x,
           const double *y,
           const int    *rnr,
           const int    *rnc,
           double *beta) {
    int nr = rnr[0];
    int nc = rnc[0];
    int p  = nc + 1;

    //
    // Cross-product C = X^t * X, where X = cbind(1, x)
    //
    double *C = calloc( p*p, sizeof(double) ); // column-major array
    CHECK_PTR( C );

    // C's first row and first column
    C[0] = nr;
    int j1;
    int jr;
    for ( int j = 0; j < nc; j++ )
    {
      j1 = j + 1;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
        C[j1] += x[i + jr];

      C[j1*p] = C[j1];
    }

    // C's diagonal
    double z;
    int ii;
    for ( int j = 0; j < nc; j++ )
    {
      ii = j+1 + (j+1)*p;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
      {
        z = x[i + jr];
        C[ii] += z * z;
      }
    }

    // Above/Below diagnoal elements of C that are not in the first row/column
    int ir;
    for ( int i = 0; i < nc; i++ )
    {
      ir = i * nr;
      for ( int j = (i+1); j < nc; j++ )
      {
        ii = (i+1) + (j+1)*p;
        jr = j * nr;
        for ( int k =  0; k < nr; k++ )
          C[ii] += x[k + ir] * x[k + jr];

        C[(j+1) + (i+1)*p] = C[ii];
      }
    }


    //
    // Cross-product beta = X^t * y, where X = cbind(1, x)
    //
    beta[0] = 0;
    for ( int i = 0; i < nr; i++ )
      beta[0] += y[i];

    for ( int j = 0; j < nc; j++ )
    {
      j1 = j + 1;
      beta[j1] = 0;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
      {
        beta[j1] += y[i]*x[i + jr];
      }
    }

    // solving A*beta = B
    int n1 = 1;
    int status = 0;
    C_fsolve(C, beta, &p, &n1, &status);
    free(C);
}


/*!
    \fn void C_flmw(double *x, double *y, double *w, int *rnr, int *rnc, double *beta, int *status)

    \brief General ordinary least squares linear model with weights.

    General ordinary least squares linean model with weights using fsolve().

    \param x      A double array derived from an nr-by-nc matrix.
    \param y      An outcome variable - an array of length nr.
    \param w      Weights array of length nr.
    \param rnr    A reference to the number of rows of the matrix associated with x.
    \param rnc    A reference to the number of columns of the matrix associated with x.

    \param beta   An output array of length nc+1 with the coefficients of the model.
    \param status The binary variable set to 0 if there are no errors (like EXIT_SUCCESS) and 1 if there were errors.

    NOTE: A common source of error is w not being normalized, that is sum(w) != 1.
          Thus, the first step in debugging should be checking if the sum of w is 1.

*/
void C_flmw(const double *x,
            const double *y,
            const double *w,
            const int    *rnr,
            const int    *rnc,
                  double *beta,
                  int    *status) {
  int nr = rnr[0];
    int nc = rnc[0];
    int p  = nc + 1;

    //
    // Cross-product C = X^t * W * X, where X = cbind(1, x) and W = diag(w)
    //
    double *C = calloc( p*p, sizeof(double) ); // column-major p-by-p symmertic array
    CHECK_PTR( C );

    // C's first row and first column
    C[0] = 1; // sum of w[i]'s
    int j1;
    int jr;
    for ( int j = 0; j < nc; j++ )
    {
      j1 = j + 1;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
        C[j1] += w[i]*x[i + jr];

      C[j1*p] = C[j1];
    }

    // C's diagonal
    double z;
    int jj;
    for ( int j = 0; j < nc; j++ )
    {
      jj = j+1 + (j+1)*p;
      jr = j * nr;
      C[jj] = 0;
      for ( int i = 0; i < nr; i++ )
      {
        z = x[i + jr];
        C[jj] += z * w[i] * z;
      }
    }

    // Above/Below diagnoal elements of C that are not in the first row/column
    int ir;
    int ij;
    for ( int i = 0; i < nc; i++ )
    {
      ir = i * nr;
      for ( int j = (i+1); j < nc; j++ )
      {
        ij = (i+1) + (j+1)*p;
        jr = j * nr;
        for ( int k =  0; k < nr; k++ )
          C[ij] += x[k + ir] * w[k] * x[k + jr];

        C[(j+1) + (i+1)*p] = C[ij];
      }
    }

    //
    // Cross-product beta = X^t * W * y, where X = cbind(1, x)
    //
    beta[0] = 0;
    for ( int i = 0; i < nr; i++ )
      beta[0] += w[i] * y[i];

    for ( int j = 0; j < nc; j++ )
    {
      beta[j+1] = 0;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
      {
        beta[j+1] += w[i] * y[i] * x[i + jr];
      }
    }

    // solving A*beta = B. In the call to fsolve beta is B and then is overwritten with the solution beta
    int n1 = 1;
    *status = 0;
    C_fsolve(C, beta, &p, &n1, status);

    if ( *status > 0 ) {
      //Rprintf("C_fsolve() returned error status 1 in C_flmw(). Setting beta to (wmean(y,w),0,0,0)");
      beta[0] = wmean(y, w, nr);
      for ( int i = 1; i < p; i++ )
        beta[i] = 0;
    }

    free(C);
}



/*!
    \fn void C_flmw_blas(double *x, double *y, double *w, int *rnr, int *rnc, double *beta)

    \brief General ordinary least squares linear model with weights with X^t*W*X and X^t*W*y computed with BLAS.

    General ordinary least squares linean model with weights using fsolve().

    \param x      A double array derived from an nr-by-nc matrix.
    \param y      An outcome variable - an array of length nr.
    \param w      Weights array of length nr.
    \param rnr    A reference to the number of rows of the matrix associated with x.
    \param rnc    A reference to the number of columns of the matrix associated with x.

    \param beta   An output array of length p+1 with the coefficients of the model.
    \param status The binary variable set to 0 if there are no errors (like EXIT_SUCCESS) and 1 if there were errors.

    NOTE: A common source of error is w not being normalized, that is sum(w) != 1.
          Thus, the first step in debugging should be checking if the sum of w is 1.

*/
void C_flmw_blas(const double *x,
                 const double *y,
                 const double *w,
                 const int    *rnr,
                 const int    *rnc,
                 double *beta,
                 int    *status)
{
    int nr = rnr[0];
    int nc = rnc[0];
    int p  = nc + 1;

    //
    // Cross-product C = X^t * W * X, where X = cbind(1, x) and W = diag(w)
    //
    double *C = calloc( p*p, sizeof(double) ); // column-major p-by-p symmertic array
    CHECK_PTR( C );

    // C's first row and first column
    C[0] = 1; // sum of w[i]'s
    int j1;
    int jr;
    for ( int j = 0; j < nc; j++ )
    {
      j1 = j + 1;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
        C[j1] += w[i]*x[i + jr];

      C[j1*p] = C[j1];
    }

    // C's diagonal
    double z;
    int jj;
    for ( int j = 0; j < nc; j++ )
    {
      jj = j+1 + (j+1)*p;
      jr = j * nr;
      C[jj] = 0;
      for ( int i = 0; i < nr; i++ )
      {
        z = x[i + jr];
        C[jj] += z * w[i] * z;
      }
    }

    // Above/Below diagnoal elements of C that are not in the first row/column
    int ir;
    int ij;
    for ( int i = 0; i < nc; i++ )
    {
      ir = i * nr;
      for ( int j = (i+1); j < nc; j++ )
      {
        ij = (i+1) + (j+1)*p;
        jr = j * nr;
        for ( int k =  0; k < nr; k++ )
          C[ij] += x[k + ir] * w[k] * x[k + jr];

        C[(j+1) + (i+1)*p] = C[ij];
      }
    }

    //
    // Cross-product beta = X^t * W * y, where X = cbind(1, x)
    //
    beta[0] = 0;
    for ( int i = 0; i < nr; i++ )
      beta[0] += w[i] * y[i];

    for ( int j = 0; j < nc; j++ )
    {
      beta[j+1] = 0;
      jr = j * nr;
      for ( int i = 0; i < nr; i++ )
      {
        beta[j+1] += w[i] * y[i] * x[i + jr];
      }
    }

    // solving A*beta = B. In the call to fsolve beta is B and then is overwritten with the solution beta
    int n1 = 1;
    *status = 0;
    C_fsolve(C, beta, &p, &n1, status);
    free(C);
}
