#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <stdlib.h>
#include <stdio.h>

#include "msr2.h"


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

  // Generate random points
  for (int i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      double random_val = ((double) rand() / (RAND_MAX));
      X[d + i * dim] = L[d] + w * random_val;
    }
  }

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

int cmp_double(const void *a, const void *b);

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
