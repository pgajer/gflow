/*
  Wasserstein distance
*/

#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

#include "msr2.h"

extern int cmp_double(const void *a, const void *b);

/**
 * @brief Computes the Wasserstein distance between two one-dimensional empirical distributions using Kahan summation.
 *
 * This function calculates the 1D Wasserstein distance (also known as Earth Mover's Distance)
 * between two samples from one-dimensional distributions. It uses the Kahan summation algorithm
 * for improved numerical stability, especially for large samples or samples with widely varying magnitudes.
 *
 * @param x Pointer to an array of doubles, representing a sample from the first distribution.
 * @param y Pointer to an array of doubles, representing a sample from the second distribution.
 * @param rn Pointer to an integer, representing the number of points in both x and y.
 * @param d Pointer to a double, where the computed Wasserstein distance will be stored.
 *
 * @note The function assumes that both input samples have the same size, specified by *rn.
 * @note The function allocates temporary memory for sorting, which is freed before returning.
 *
 * @see https://en.wikipedia.org/wiki/Kahan_summation_algorithm for details on the Kahan summation algorithm.
 */
void C_wasserstein_distance_1D(const double *x,
                               const double *y,
                               const int    *rn,
                                     double *d)
{
    int n = *rn;
    double *x_sorted = (double *)malloc(n * sizeof(double));
    CHECK_PTR(x_sorted);

    double *y_sorted = (double *)malloc(n * sizeof(double));
    CHECK_PTR(y_sorted);

    // Copy and sort the samples
    memcpy(x_sorted, x, n * sizeof(double));
    memcpy(y_sorted, y, n * sizeof(double));
    qsort(x_sorted, n, sizeof(double), cmp_double);
    qsort(y_sorted, n, sizeof(double), cmp_double);

    // Compute the Wasserstein distance using Kahan summation
    double sum = 0.0;
    double c = 0.0;  // A running compensation for lost low-order bits
    for (int i = 0; i < n; i++) {
        double y = fabs(x_sorted[i] - y_sorted[i]) - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    *d = sum / n;

    // Free the allocated memory
    free(x_sorted);
    free(y_sorted);
}
