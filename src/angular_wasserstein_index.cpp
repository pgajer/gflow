#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "stats_utils.h"
#include "kNN.h"
#include "wasserstein_dist.h" // for C_wasserstein_distance_1D()

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void C_wasserstein_distance_1D(const double *x,
                               const double *y,
                               const int    *rn,
                                     double *d);

extern "C" {
    SEXP S_angular_wasserstein_index(SEXP s_X, SEXP s_Y, SEXP s_k);
}


/**
 * @brief Compute the Angular Wasserstein Index between two sets of points.
 *
 * This function calculates the Angular Wasserstein Index, a measure of dissimilarity
 * between two point sets X and Y, which is particularly useful for comparing denoised
 * versions of a dataset. The index is based on the angular distribution of k-nearest
 * neighbors for each point.
 *
 * The Angular Wasserstein Index is defined as:
 * I_W(X, Y) = (1 / (n * (k-1))) * sum(d_W(alpha_X(x), alpha_Y(x)))
 * where:
 * - n is the number of points in X (and Y)
 * - k is the number of nearest neighbors considered
 * - alpha_X(x) and alpha_Y(x) are the angular distributions of the k-1 nearest neighbors
 *   of point x in X and Y respectively
 * - d_W is the 1D Wasserstein distance between these angular distributions
 *
 * @param X The first set of points, represented as a vector of vectors where each inner
 *          vector is a point in n-dimensional space.
 * @param Y The second set of points, with the same structure and dimensions as X.
 * @param k The number of nearest neighbors to consider (including the point itself).
 *          Must be at least 2.
 *
 * @return The Angular Wasserstein Index between X and Y.
 *
 * @throws std::invalid_argument if X and Y have different dimensions or if k < 2.
 *
 * @note This function assumes that X and Y represent the same set of points before
 *       and after some transformation (e.g., denoising), and that they are in
 *       correspondence (i.e., X[i] corresponds to Y[i]).
 *
 * @note The function uses a k-nearest neighbors algorithm and a 1D Wasserstein
 *       distance calculation as subroutines. The efficiency of these implementations
 *       will significantly impact the overall performance of this function.
 *
 * @Rf_warning This function may be computationally expensive for large datasets or
 *          high values of k. Consider performance implications for your use case.
 *
 * Example usage:
 * @code
 * std::vector<std::vector<double>> X = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
 * std::vector<std::vector<double>> Y = {{1.1, 2.1}, {3.1, 4.1}, {5.1, 6.1}};
 * int k = 3;
 * double index = angular_wasserstein_index(X, Y, k);
 * @endcode
 */
#include <numeric>   // std::inner_product
#include <algorithm> // std::clamp
#include <cmath>

double angular_wasserstein_index(const std::vector<std::vector<double>>& X,
                                 const std::vector<std::vector<double>>& Y,
                                 int k) {
    // Preconditions are reportedly checked in the R wrapper, but add cheap asserts here too.
    const int n_X_points = static_cast<int>(X.size());
    const int dim        = static_cast<int>(X.empty() ? 0 : X[0].size());
    if (n_X_points == 0 || dim == 0 || k <= 1) return 0.0; // nothing to compare
    if (static_cast<int>(Y.size()) != n_X_points || static_cast<int>(Y[0].size()) != dim)
        Rf_error("X and Y must have identical shapes");

    const int k_minus_one = k - 1;

    // RAII: result holds flat vectors of size n*k
    const auto knn_res = kNN(X, k);  // indices and distances are std::vector<>

    double total_dist = 0.0;
    std::vector<double> angles_X(k_minus_one);
    std::vector<double> angles_Y(k_minus_one);

    // scratch allocated once per point
    std::vector<double> reference_X(dim), reference_Y(dim);
    std::vector<double> vec_X(dim), vec_Y(dim);

    for (int point = 0; point < n_X_points; ++point) {
        const int ref_idx = knn_res.indices[static_cast<size_t>(point) * k + (k - 1)];

        // reference direction (last neighbor)
        for (int d = 0; d < dim; ++d) {
            reference_X[d] = X[ref_idx][d] - X[point][d];
            reference_Y[d] = Y[ref_idx][d] - Y[point][d];
        }

        const double ref_norm_X = std::sqrt(std::inner_product(reference_X.begin(), reference_X.end(),
                                                               reference_X.begin(), 0.0));
        const double ref_norm_Y = std::sqrt(std::inner_product(reference_Y.begin(), reference_Y.end(),
                                                               reference_Y.begin(), 0.0));

        // Guard against zero norms (coincident points). If zero, angle defaults to 0.
        const bool ref_X_zero = (ref_norm_X == 0.0);
        const bool ref_Y_zero = (ref_norm_Y == 0.0);

        // j = 1..k-1 (skip j=0 which is the point itself)
        for (int j = 1; j < k; ++j) {
            const int neighbor_idx = knn_res.indices[static_cast<size_t>(point) * k + j];

            for (int d = 0; d < dim; ++d) {
                vec_X[d] = X[neighbor_idx][d] - X[point][d];
                vec_Y[d] = Y[neighbor_idx][d] - Y[point][d];
            }

            const double vec_norm_X = std::sqrt(std::inner_product(vec_X.begin(), vec_X.end(),
                                                                   vec_X.begin(), 0.0));
            const double vec_norm_Y = std::sqrt(std::inner_product(vec_Y.begin(), vec_Y.end(),
                                                                   vec_Y.begin(), 0.0));

            double cosX = 1.0, cosY = 1.0;
            if (!ref_X_zero && vec_norm_X != 0.0) {
                const double dotX = std::inner_product(vec_X.begin(), vec_X.end(), reference_X.begin(), 0.0);
                cosX = std::clamp(dotX / (vec_norm_X * ref_norm_X), -1.0, 1.0);
            }
            if (!ref_Y_zero && vec_norm_Y != 0.0) {
                const double dotY = std::inner_product(vec_Y.begin(), vec_Y.end(), reference_Y.begin(), 0.0);
                cosY = std::clamp(dotY / (vec_norm_Y * ref_norm_Y), -1.0, 1.0);
            }

            angles_X[j - 1] = std::acos(cosX);
            angles_Y[j - 1] = std::acos(cosY);
        }

        double dist = 0.0;
        C_wasserstein_distance_1D(angles_X.data(), angles_Y.data(),
                                  const_cast<int*>(&k_minus_one), &dist);
        total_dist += dist;
    }

    // No deletes: knn_res holds std::vector and frees itself
    return total_dist / (n_X_points * k_minus_one);
}

/**
 * @brief R interface for computing the Angular Wasserstein Index between two sets of points.
 *
 * This function serves as an interface between R and the C++ implementation of the
 * Angular Wasserstein Index calculation. It takes R objects as input, converts them
 * to C++ data structures, computes the index, and returns the result as an R object.
 *
 * The Angular Wasserstein Index is a measure of dissimilarity between two point sets,
 * particularly useful for comparing denoised versions of a dataset. It is based on
 * the angular distribution of k-nearest neighbors for each point.
 *
 * @param s_X SEXP object representing the first set of points as an R matrix.
 *            Each row should represent a point in n-dimensional space.
 * @param s_Y SEXP object representing the second set of points as an R matrix.
 *            It should have the same dimensions as s_X.
 * @param s_k SEXP object representing the number of nearest neighbors to consider.
 *            This should be an R integer with a single value >= 2.
 *
 * @return SEXP object containing a single numeric value representing the
 *         Angular Wasserstein Index between X and Y.
 *
 * @note This function is intended to be called from R using .Call(). It should
 *       not be called directly from C or C++ code.
 *
 * @Rf_warning Ensure that the input matrices s_X and s_Y are numeric and have the
 *          same dimensions. The function will raise an R Rf_error if these
 *          conditions are not met.
 *
 * R usage:
 * ```r
 * .Call("S_angular_wasserstein_index", X, Y, as.integer(k))
 * ```
 * where X and Y are numeric matrices and k is an integer >= 2.
 *
 * Implementation details:
 * 1. Converts R matrices to C++ vector<vector<double>> using Rmatrix_to_cpp().
 * 2. Extracts k value from R integer.
 * 3. Calls C++ angular_wasserstein_index() function.
 * 4. Returns result as R numeric vector of length 1.
 *
 * @see angular_wasserstein_index() for the underlying C++ implementation.
 * @see Rmatrix_to_cpp() for the R to C++ matrix conversion function.
 */
SEXP S_angular_wasserstein_index(SEXP s_X, SEXP s_Y, SEXP s_k) {

    // Converting R matrices to C++ vectors
    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));
    std::vector<std::vector<double>> Y = std::move(*Rmatrix_to_cpp(s_Y));

    // Extracting k value
    int k = INTEGER(s_k)[0];

    // Computing the Angular Wasserstein Index
    double d = angular_wasserstein_index(X, Y, k);

    // Creating an R numeric vector to hold the result
    SEXP result = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(result)[0] = d;

    UNPROTECT(1);

    return result;
}
