#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <cmath>
#include <memory>

#include "msr2.h"
#include "cpp_utils.h"
#include "SEXP_cpp_conversion_utils.h"
#include "graph_diffussion_smoother.h"
#include "stats_utils.h"
#include "kernels.h"
#include "kNN.h"

// Forward declarations of helper functions
SEXP create_R_list(const std::vector<std::vector<std::vector<double>>>& X_traj,
                   const std::vector<double>& cv_errors,
                   int best_step);

extern void initialize_kernel(int ikernel);
extern void (*kernel_fn)(const double*, int, double*);
extern kNN_result_t kNN(const std::vector<std::vector<double>>& X, int k);

// Define the result struct
struct mean_shift_smoother_results_t {
    std::vector<std::vector<std::vector<double>>> X_traj;
    std::vector<double> cv_errors;
    int best_step;
};

extern "C" {
    SEXP S_mean_shift_matrix_smoother(SEXP s_X,
                                      SEXP s_k,
                                      SEXP s_n_steps,
                                      SEXP s_step_size,
                                      SEXP s_n_CVs,
                                      SEXP s_n_CV_folds,
                                      SEXP s_ikernel,
                                      SEXP s_dist_normalization_factor);
}

/**
 * @brief Performs mean shift smoothing on a dataset and computes cross-validation errors.
 *
 * This function applies the mean shift algorithm to smooth a given dataset. It iteratively
 * moves each point towards areas of higher density in the feature space. The function also
 * performs cross-validation to estimate the error at each smoothing step.
 *
 * @param X The input dataset, where each inner vector represents a point and its features.
 *          Type: const std::vector<std::vector<double>>&
 *          Each inner vector should have the same size, representing the number of features.
 *
 * @param k The number of nearest neighbors to consider for density estimation.
 *          Type: int
 *          Should be a positive integer, typically much smaller than the dataset size.
 *
 * @param n_steps The maximum number of smoothing steps to perform.
 *                Type: int
 *                Should be a positive integer.
 *
 * @param step_size The step size for updating point positions in each iteration.
 *                  Type: double
 *                  Should be a positive value, typically between 0 and 1.
 *
 * @param n_CVs The number of cross-validation runs to perform.
 *              Type: int
 *              Should be a positive integer.
 *
 * @param n_CV_folds The number of folds to use in each cross-validation run.
 *                   Type: int
 *                   Should be a positive integer, typically between 3 and 10.
 *
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 *                Type: int
 *                Default: 1 (Epanechnikov kernel)
 *
 * @param dist_normalization_factor A scaling factor applied to the maximum distance between
 *                                  a vertex and its neighbors. This ensures non-zero weights
 *                                  even when all distances are equal.
 *                                  Type: double
 *                                  Default: 1.01
 *                                  Should be greater than 1.0.
 *
 * @return std::unique_ptr<mean_shift_smoother_results_t>
 *         A unique pointer to a struct containing:
 *         - X_traj: A vector of smoothed datasets, one for each step.
 *         - cv_errors: A vector of cross-validation errors, one for each step.
 *         - best_step: The step with the lowest cross-validation error.
 *
 * @throws May throw exceptions related to memory allocation or invalid input parameters.
 *
 * @note This function assumes the existence of a kNN function with the signature:
 *       kNN_result_t kNN(const std::vector<std::vector<double>>& X, int k);
 *       where kNN_result_t is a struct containing int* indices and double* distances.
 *
 * @note The function uses a random number generator for cross-validation. The seed is fixed
 *       for reproducibility, but this may be adjusted if different random splits are desired.
 *
 * @warning This function may be computationally expensive for large datasets or high values of n_steps.
 *
 * Example usage:
 * @code
 * std::vector<std::vector<double>> X = {{1.0, 2.0}, {2.0, 3.0}, {5.0, 6.0}, {-1.0, -2.0}};
 * int k = 2;
 * int n_steps = 10;
 * double step_size = 0.1;
 * int n_CVs = 5;
 * int n_CV_folds = 3;
 *
 * auto results = mean_shift_smoother(X, k, n_steps, step_size, n_CVs, n_CV_folds);
 *
 * std::cout << "Best step: " << results->best_step << std::endl;
 * std::cout << "CV error at best step: " << results->cv_errors[results->best_step] << std::endl;
 * @endcode
 */
std::unique_ptr<mean_shift_smoother_results_t> mean_shift_matrix_smoother(const std::vector<std::vector<double>>& X,
                                                                          int k,
                                                                          int n_steps,
                                                                          double step_size,
                                                                          int n_CVs,
                                                                          int n_CV_folds,
                                                                          int ikernel = 1,
                                                                          double dist_normalization_factor = 1.01) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    initialize_kernel(ikernel);

    // Precompute KNN for all points
    auto knn_res = kNN(X, k);

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(k);
    std::vector<double> local_distances(k);

    // Initialize X_traj
    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimate gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);

            double max_dist = 0.0;
            for (int j = 0; j < k; ++j) {
                local_distances[j] = knn_res.distances[point * k + j];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k; ++j) {
                local_distances[j] /= max_dist;
            }

            kernel_fn(local_distances.data(), k, kernel_weights.data());

            double total_weight = 0.0;
            for (int i = 0; i < k; ++i) {
                int neighbor_idx = knn_res.indices[point * k + i];
                for (int j = 0; j < n_features; ++j) {
                    new_point[j] += kernel_weights[i] * X[neighbor_idx][j];
                }
                total_weight += kernel_weights[i];
            }

            // Compute gradient
            double gradient_norm = 0.0;
            for (int j = 0; j < n_features; ++j) {
                new_point[j] /= total_weight;
                gradient_vector_field[point][j] = new_point[j] - X_traj[step-1][point][j];
                gradient_norm += gradient_vector_field[point][j] * gradient_vector_field[point][j];
            }
            gradient_norm = std::sqrt(gradient_norm);

            // Normalize gradient
            if (gradient_norm > 0) {
                for (int j = 0; j < n_features; ++j) {
                    gradient_vector_field[point][j] /= gradient_norm;
                }
            }
        }

        // Update X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            double max_dist = 0.0;
            for (int j = 0; j < k; ++j) {
                local_distances[j] = knn_res.distances[point * k + j];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k; ++j) {
                local_distances[j] /= max_dist;
            }

            kernel_fn(local_distances.data(), k, kernel_weights.data());

            for (int i = 0; i < k; ++i) {
                int neighbor_idx = knn_res.indices[point * k + i];
                double proj = 0.0;
                for (int j = 0; j < n_features; ++j) {
                    proj += (X[neighbor_idx][j] - X[point][j]) * gradient_vector_field[point][j];
                }
                gradient_magnitude += kernel_weights[i] * std::abs(proj);
                total_weight += kernel_weights[i];
            }
            gradient_magnitude /= total_weight;

            // Update point position
            for (int j = 0; j < n_features; ++j) {
                X_traj[step][point][j] = X_traj[step-1][point][j] + step_size * gradient_magnitude * gradient_vector_field[point][j];
            }
        }
    }

    // Cross-validation error estimation
    std::vector<double> cv_errors(n_steps, 0.0);
    std::mt19937 gen(42);  // Random number generator

    for (int cv = 0; cv < n_CVs; ++cv) {
        std::vector<int> indices(n_points);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int fold = 0; fold < n_CV_folds; ++fold) {
            int fold_size = n_points / n_CV_folds;
            int start = fold * fold_size;
            int end = (fold == n_CV_folds - 1) ? n_points : (fold + 1) * fold_size;

            std::vector<int> test_indices(indices.begin() + start, indices.begin() + end);
            std::vector<int> train_indices;
            train_indices.reserve(n_points - test_indices.size());
            for (int i = 0; i < n_points; ++i) {
                if (std::find(test_indices.begin(), test_indices.end(), i) == test_indices.end()) {
                    train_indices.push_back(i);
                }
            }

            for (int step = 0; step < n_steps; ++step) {
                double error = 0.0;
                for (int test_point : test_indices) {
                    std::vector<double> estimate(n_features, 0.0);
                    double total_weight = 0.0;

                    double max_dist = 0.0;
                    for (int j = 0; j < k; ++j) {
                        local_distances[j] = knn_res.distances[test_point * k + j];
                        if (local_distances[j] > max_dist) max_dist = local_distances[j];
                    }

                    if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                    max_dist *= dist_normalization_factor;

                    for (int j = 0; j < k; ++j) {
                        local_distances[j] /= max_dist;
                    }

                    kernel_fn(local_distances.data(), k, kernel_weights.data());

                    for (int i = 0; i < k; ++i) {
                        int neighbor = knn_res.indices[test_point * k + i];
                        if (std::find(train_indices.begin(), train_indices.end(), neighbor) != train_indices.end()) {
                            for (int j = 0; j < n_features; ++j) {
                                estimate[j] += kernel_weights[i] * X_traj[step][neighbor][j];
                            }
                            total_weight += kernel_weights[i];
                        }
                    }

                    for (int j = 0; j < n_features; ++j) {
                        estimate[j] /= total_weight;
                        error += std::abs(X[test_point][j] - estimate[j]);
                    }
                }
                cv_errors[step] += error / (test_indices.size() * n_features);
            }
        }
    }

    for (double& error : cv_errors) {
        error /= (n_CVs * n_CV_folds);
    }

    // Find the step with minimum CV error
    int best_step = std::min_element(cv_errors.begin(), cv_errors.end()) - cv_errors.begin();

    // Create and return the results struct as a unique pointer
    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->cv_errors = std::move(cv_errors);
    results->best_step = best_step;

    return results;
}

#if 0
std::unique_ptr<mean_shift_smoother_results_t> old_mean_shift_smoother(const std::vector<std::vector<double>>& X,
                                                                   int k,
                                                                   int n_steps,
                                                                   double step_size,
                                                                   int n_CVs,
                                                                   int n_CV_folds,
                                                                   int ikernel = 1,
                                                                   double dist_normalization_factor = 1.01) {

    const double epsilon = 1e-10;  // Small constant to avoid division by zero
    const int n_points = X.size();
    const int n_features = X[0].size();

    initialize_kernel(ikernel);

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(k + 1);
    std::vector<double> local_distances(k + 1);

    // Initialize X_traj
    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimate gradient for each point
        for (int point = 0; point < n_points; ++point) {
            auto knn_res = kNN(X, k);
            std::vector<double> new_point(n_features, 0.0);

            local_distances[0] = 0;  // Distance to self is 0
            double max_dist = 0.0;
            for (int j = 0; j < k; ++j) {
                local_distances[j + 1] = knn_res.distances[point * k + j];
                if (local_distances[j + 1] > max_dist)
                    max_dist = local_distances[j + 1];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero

            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k; ++j)
                local_distances[j + 1] /= max_dist;

            kernel_fn(local_distances.data(), k + 1, kernel_weights.data());

            double total_weight = 0.0;
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < n_features; ++j) {
                    new_point[j] += kernel_weights[i] * X[knn_res.indices[point * k + i]][j];
                }
                total_weight += kernel_weights[i];
            }

            // Compute gradient
            double gradient_norm = 0.0;
            for (int j = 0; j < n_features; ++j) {
                new_point[j] /= total_weight;
                gradient_vector_field[point][j] = new_point[j] - X_traj[step-1][point][j];
                gradient_norm += gradient_vector_field[point][j] * gradient_vector_field[point][j];
            }
            gradient_norm = std::sqrt(gradient_norm);

            // Normalize gradient
            if (gradient_norm > epsilon) {
                for (int j = 0; j < n_features; ++j) {
                    gradient_vector_field[point][j] /= gradient_norm;
                }
            }
        }

        // Update X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            auto knn_res = kNN(X, k);
            double total_weight = 0.0;

            local_distances[0] = 0;  // Distance to self is 0
            double max_dist = 0.0;
            for (int j = 0; j < k; ++j) {
                local_distances[j + 1] = knn_res.distances[point * k + j];
                if (local_distances[j + 1] > max_dist)
                    max_dist = local_distances[j + 1];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero

            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k; ++j)
                local_distances[j + 1] /= max_dist;

            kernel_fn(local_distances.data(), k + 1, kernel_weights.data());

            for (int i = 0; i < k; ++i) {
                double proj = 0.0;
                for (int j = 0; j < n_features; ++j) {
                    proj += (X[knn_res.indices[point * k + i]][j] - X[point][j]) * gradient_vector_field[point][j];
                }
                gradient_magnitude += kernel_weights[i] * std::abs(proj);
                total_weight += kernel_weights[i];
            }
            gradient_magnitude /= total_weight;

            // Update point position
            for (int j = 0; j < n_features; ++j) {
                X_traj[step][point][j] = X_traj[step-1][point][j] + step_size * gradient_magnitude * gradient_vector_field[point][j];
            }
        }
    }

    // Cross-validation error estimation
    std::vector<double> cv_errors(n_steps, 0.0);
    std::mt19937 gen(42);  // Random number generator

    for (int cv = 0; cv < n_CVs; ++cv) {
        std::vector<int> indices(n_points);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int fold = 0; fold < n_CV_folds; ++fold) {
            int fold_size = n_points / n_CV_folds;
            int start = fold * fold_size;
            int end = (fold == n_CV_folds - 1) ? n_points : (fold + 1) * fold_size;

            std::vector<int> test_indices(indices.begin() + start, indices.begin() + end);
            std::vector<int> train_indices;
            train_indices.reserve(n_points - test_indices.size());
            for (int i = 0; i < n_points; ++i) {
                if (std::find(test_indices.begin(), test_indices.end(), i) == test_indices.end()) {
                    train_indices.push_back(i);
                }
            }

            for (int step = 0; step < n_steps; ++step) {
                double error = 0.0;
                for (int test_point : test_indices) {
                    std::vector<double> estimate(n_features, 0.0);
                    auto knn_res = kNN(X, k);
                    double total_weight = 0.0;

                    for (int i = 0; i < k; ++i) {
                        int neighbor = knn_res.indices[test_point * k + i];
                        if (std::find(train_indices.begin(), train_indices.end(), neighbor) != train_indices.end()) {
                            double weight = kernel(knn_res.distances[test_point * k + i]);
                            for (int j = 0; j < n_features; ++j) {
                                estimate[j] += weight * X_traj[step][neighbor][j];
                            }
                            total_weight += weight;
                        }
                    }

                    for (int j = 0; j < n_features; ++j) {
                        estimate[j] /= total_weight;
                        error += std::abs(X[test_point][j] - estimate[j]);
                    }
                }
                cv_errors[step] += error / (test_indices.size() * n_features);
            }
        }
    }

    for (double& error : cv_errors) {
        error /= (n_CVs * n_CV_folds);
    }

    // Find the step with minimum CV error
    int best_step = std::min_element(cv_errors.begin(), cv_errors.end()) - cv_errors.begin();

    // Create and return the results struct as a unique pointer
    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->cv_errors = std::move(cv_errors);
    results->best_step = best_step;

    return results;
}
#endif


/**
 * @brief R-compatible version of mean shift smoothing algorithm with cross-validation.
 *
 * This function performs mean shift smoothing on a dataset and computes cross-validation errors.
 * It is designed to be called from R using .Call() and uses SEXP types for all parameters and
 * the return value.
 *
 * @param s_X SEXP (list of numeric vectors)
 *        The input dataset, where each inner vector represents a point and its features.
 *        In R, this should be a list of numeric vectors or a matrix.
 *
 * @param s_k SEXP (integer)
 *        The number of nearest neighbors to consider for density estimation.
 *        In R, this should be a single integer value.
 *
 * @param s_n_steps SEXP (integer)
 *        The maximum number of smoothing steps to perform.
 *        In R, this should be a single integer value.
 *
 * @param s_step_size SEXP (numeric)
 *        The step size for updating point positions in each iteration.
 *        In R, this should be a single numeric value, typically between 0 and 1.
 *
 * @param s_n_CVs SEXP (integer)
 *        The number of cross-validation runs to perform.
 *        In R, this should be a single integer value.
 *
 * @param s_n_CV_folds SEXP (integer)
 *        The number of folds to use in each cross-validation run.
 *        In R, this should be a single integer value, typically between 3 and 10.
 *
 * @param s_ikernel SEXP (integer)
 *        An integer specifying the kernel function to use. Valid values are:
 *        1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 *        In R, this should be a single integer value.
 *
 * @param s_dist_normalization_factor SEXP (numeric)
 *        A scaling factor applied to the maximum distance between a vertex and its neighbors.
 *        This ensures non-zero weights even when all distances are equal.
 *        In R, this should be a single numeric value, typically slightly above 1.0.
 *
 * @return SEXP (list)
 *         Returns an R list containing:
 *         - X_traj: A list of matrices, each representing the smoothed dataset at each step.
 *         - cv_errors: A numeric vector of cross-validation errors, one for each step.
 *         - best_step: An integer indicating the step with the lowest cross-validation error.
 *
 * @note This function assumes the existence of helper functions:
 *       - R_list_of_dvectors_to_cpp_vector_of_dvectors: Converts R list/matrix to C++ vector of vectors.
 *       - initialize_kernel: Sets up the kernel function based on ikernel value.
 *       - kernel_fn: The kernel function used for weighting.
 *       - kNN: Computes k-nearest neighbors for the dataset.
 *
 * @warning This function may be computationally expensive for large datasets or high values of n_steps.
 *          Ensure sufficient memory is available, especially for large datasets.
 *
 * Example usage in R:
 * @code
 * # Assuming the shared object is compiled and loaded
 * dyn.load("path_to_your_compiled_shared_object.so")
 *
 * # Prepare input data
 * X <- matrix(rnorm(100*2), ncol=2)  # 100 points in 2D space
 * k <- 5
 * n_steps <- 10
 * step_size <- 0.1
 * n_CVs <- 5
 * n_CV_folds <- 3
 * ikernel <- 1  # Epanechnikov kernel
 * dist_normalization_factor <- 1.01
 *
 * # Call the function
 * result <- .Call("S_mean_shift_matrix_smoother",
 *                 X,
 *                 as.integer(k),
 *                 as.integer(n_steps),
 *                 as.double(step_size),
 *                 as.integer(n_CVs),
 *                 as.integer(n_CV_folds),
 *                 as.integer(ikernel),
 *                 as.double(dist_normalization_factor))
 *
 * # Access results
 * X_traj <- result$X_traj
 * cv_errors <- result$cv_errors
 * best_step <- result$best_step
 *
 * # Plot CV errors
 * plot(cv_errors, type='l', xlab='Step', ylab='CV Error')
 * points(best_step, cv_errors[best_step], col='red', pch=19)
 * @endcode
 */
SEXP S_mean_shift_matrix_smoother(SEXP s_X,
                                  SEXP s_k,
                                  SEXP s_n_steps,
                                  SEXP s_step_size,
                                  SEXP s_n_CVs,
                                  SEXP s_n_CV_folds,
                                  SEXP s_ikernel,
                                  SEXP s_dist_normalization_factor) {

    // Convert SEXP inputs to C++ types
    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));
    int k = INTEGER(s_k)[0];
    int n_steps = INTEGER(s_n_steps)[0];
    double step_size = REAL(s_step_size)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    // Call the C++ function
    auto results = mean_shift_matrix_smoother(X, k, n_steps, step_size, n_CVs, n_CV_folds, ikernel, dist_normalization_factor);

    // Check if results is valid
    if (!results) {
        Rf_error("mean_shift_smoother returned null pointer");
    }

    // Check if X_traj is valid
    if (results->X_traj.empty() || results->X_traj[0].empty() || results->X_traj[0][0].empty()) {
        Rf_error("mean_shift_smoother returned invalid X_traj");
    }

    // Check if cv_errors is valid
    if (results->cv_errors.empty()) {
        Rf_error("mean_shift_smoother returned empty cv_errors");
    }

    // Check if best_step is valid
    if (results->best_step < 0 || results->best_step >= static_cast<int>(results->cv_errors.size())) {
        Rf_error("mean_shift_smoother returned invalid best_step");
    }

    // Convert results to R list
    return create_R_list(results->X_traj, results->cv_errors, results->best_step);
}



/**
 * @brief Helper function to convert mean shift smoothing results to an R list.
 *
 * This function takes the results of the mean shift smoothing algorithm (trajectory of points,
 * cross-validation errors, and the best step) and converts them into an R list structure.
 * It is designed to be used internally by the S_mean_shift_matrix_smoother function to prepare
 * its return value for R.
 *
 * @param X_traj const std::vector<std::vector<std::vector<double>>>&
 *        A 3D vector representing the trajectory of points across all smoothing steps.
 *        - First dimension: steps
 *        - Second dimension: points
 *        - Third dimension: features of each point
 *
 * @param cv_errors const std::vector<double>&
 *        A vector of cross-validation errors, one for each smoothing step.
 *
 * @param best_step int
 *        The index of the step with the lowest cross-validation error.
 *
 * @return SEXP
 *         Returns an R list containing:
 *         - X_traj: A list of matrices, each representing the smoothed dataset at each step.
 *         - cv_errors: A numeric vector of cross-validation errors.
 *         - best_step: An integer indicating the step with the lowest cross-validation error.
 *
 * @note This function uses R's C API to create R objects. It handles the necessary memory
 *       protection (PROTECT/UNPROTECT) to prevent garbage collection issues.
 *
 * @warning This function assumes that the input vectors are not empty. It does not perform
 *          extensive error checking on the inputs.
 *
 * Example usage (within C++ code):
 * @code
 * std::vector<std::vector<std::vector<double>>> X_traj = // ... populated during smoothing
 * std::vector<double> cv_errors = // ... computed during cross-validation
 * int best_step = // ... determined from cv_errors
 *
 * SEXP result = create_R_list(X_traj, cv_errors, best_step);
 * // result is now ready to be returned to R
 * @endcode
 */
SEXP create_R_list(const std::vector<std::vector<std::vector<double>>>& X_traj,
                   const std::vector<double>& cv_errors,
                   int best_step) {

    int n_steps = X_traj.size();
    int n_points = X_traj[0].size();
    int n_features = X_traj[0][0].size();

    int nprot = 0;
    SEXP X_traj_r = PROTECT(X_traj_r = allocVector(VECSXP, n_steps)); nprot++;
    for (int i = 0; i < n_steps; ++i) {
        SEXP step_matrix = PROTECT(step_matrix = allocMatrix(REALSXP, n_points, n_features));
        double* ptr = REAL(step_matrix);
        for (int j = 0; j < n_points; ++j) {
            for (int k = 0; k < n_features; ++k) {
                ptr[j + k * n_points] = X_traj[i][j][k];
            }
        }
        SET_VECTOR_ELT(X_traj_r, i, step_matrix);
        UNPROTECT(1);
    }

    // cv_errors
    SEXP cv_errors_r = PROTECT(cv_errors_r = allocVector(REALSXP, cv_errors.size())); nprot++;
    std::copy(cv_errors.begin(), cv_errors.end(), REAL(cv_errors_r));

    // best_step
    SEXP best_step_r = PROTECT(best_step_r = allocVector(INTSXP, 1)); nprot++;
    INTEGER(best_step_r)[0] = best_step;

    // Set list elements
    SEXP result = PROTECT(result = allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(result, 0, X_traj_r);
    SET_VECTOR_ELT(result, 1, cv_errors_r);
    SET_VECTOR_ELT(result, 2, best_step_r);

    // Set names
    SEXP names = PROTECT(names = allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("X_traj"));
    SET_STRING_ELT(names, 1, mkChar("cv_errors"));
    SET_STRING_ELT(names, 2, mkChar("best_step"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result;
}
