#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>

extern int cmp_double(const void *a, const void *b); // in c_utils.c

/**
 * @brief Samples from the uniform distribution over (K-1)-dimensional simplex.
 *
 * This function generates a random array of K non-negative doubles that sum up to 1,
 * effectively sampling a point from the (K-1)-dimensional simplex. It uses the
 * "ordered differences" method to ensure uniform distribution over the simplex.
 *
 * @details The algorithm works as follows:
 * 1. Generate K-1 uniform random numbers between 0 and 1.
 * 2. Sort these numbers in ascending order.
 * 3. Calculate the differences between adjacent pairs of these sorted numbers,
 *    including 0 as the lower bound and 1 as the upper bound.
 * 4. These differences form the K components of the sampled point on the simplex.
 *
 * This method ensures that the resulting point is uniformly distributed over the simplex.
 *
 * @param rK Pointer to an integer specifying the dimension of the simplex plus one.
 *           The value pointed to by rK should be greater than 1.
 * @param[out] lambda Pointer to a pre-allocated array of doubles where the output will be stored.
 *                    This array should have at least K elements, where K is the value pointed to by rK.
 *
 * @pre The caller must ensure that:
 *      - rK is not NULL and points to a valid integer greater than 1.
 *      - lambda is not NULL and points to a pre-allocated array of at least K doubles.
 *
 * @post Upon return:
 *       - lambda[0] to lambda[K-1] will contain K non-negative doubles that sum to 1.
 *       - The random number generator state will be updated.
 *
 * @note This function uses R's random number generator functions (GetRNGstate and PutRNGstate).
 *       It is designed to be called from R via .Call() or similar interfaces.
 *
 * @Rf_warning This function does not perform input validation. Incorrect inputs may lead to
 *          undefined behavior.
 *
 * @see For background on the algorithm:
 *      Smith, R. L., & Tierney, L. (1996). Exact Transition Probabilities for the
 *      Independence Metropolis Sampler. Technical Report.
 *
 * @todo Evaluate numerical stability for very large K values.
 */
void C_runif_simplex( const int *rK, double *lambda)
{
    int K = rK[0];
    int K1 = K - 1;

    GetRNGstate();

    for ( int i = 0; i < K1; i++ )
      lambda[i] = runif(0.0,1.0);

    PutRNGstate();

    qsort(lambda, K1, sizeof(double), cmp_double);

    lambda[K - 1] =  1.0 - lambda[K - 2];

    for ( int i = K - 2; i > 0; i-- )
      lambda[i] -= lambda[i-1];
}

/*
  \brief Samples from the Dirichlet distribution with parameter w.

  Generates a random array of K non-negative doubles that sum up to 1. The mean of that distribution is w.

  \param w       An array of n non-negative values constituting the parameters of the Dirichlet distribution from which this routine is doing random sampling.
  \param rn      A reference to the length of w and lambda.
  \param lambda  The output array.

*/
void C_rsimplex(const double *w, const int *rn, double *lambda)
{
    int n = rn[0];

    GetRNGstate();
    double s = 0;
    for ( int i = 0; i < n; i++ )
    {
      lambda[i] = rgamma(w[i],1.0);
      s += lambda[i];
    }
    PutRNGstate();

    for ( int i = 0; i < n; i++ )
      lambda[i] /= s;
}
