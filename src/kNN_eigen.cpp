#include "kNN.h"
#include <Eigen/Sparse>

/**
 * @brief Compute k-nearest neighbors from Eigen sparse matrix
 *
 * This function provides a C++ interface to kNN computation that works
 * directly with Eigen sparse matrices, avoiding SEXP manipulation.
 * It converts the sparse matrix to the dense format expected by the
 * ANN library, computes kNN, and returns results in C++ containers.
 *
 * @param X Feature matrix (n_points × n_features, sparse)
 * @param k Number of nearest neighbors
 * @return knn_result_t containing indices and distances
 *
 * @throws std::invalid_argument if matrix dimensions are invalid or k is out of range
 */
knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
) {
    const int n = static_cast<int>(X.rows());
    const int d = static_cast<int>(X.cols());

    if (n <= 0 || d <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    if (k <= 0 || k > n) {
        throw std::invalid_argument("k must be in range [1, n]");
    }

    // Convert Eigen sparse matrix to dense vector-of-vectors format
    // required by kNN() function
    std::vector<std::vector<double>> X_dense(n, std::vector<double>(d));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            X_dense[i][j] = X.coeff(i, j);
        }
    }

    // Call the existing kNN helper function
    return kNN(X_dense, k);
}
