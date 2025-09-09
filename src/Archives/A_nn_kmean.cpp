#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <memory>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <stack>
#include <numeric>
#include <unordered_map>
#include <random>
#include <chrono>

#include "msr2.h"
#include "msr2_cpp_utils.h"
#include "msr2_Cpp_to_R_utils.h"
#include "msr2_graph_diffussion_smoother.h"
#include "msr2_stats_utils.h"
#include "msr2_kernels.h"

std::unique_ptr<std::vector<std::vector<double>>> dist_to_weights(const std::vector<std::vector<int>>& graph,
                                                                  const std::vector<std::vector<double>>& d,
                                                                  int ikernel,
                                                                  double dist_normalization_factor);

std::unique_ptr<std::vector<double>> cv_imputation(const std::set<int>& test_set,
                                                   const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& edge_lengths,
                                                   const std::vector<double>& y,
                                                   bool y_binary,
                                                   imputation_method_t imputation_method,
                                                   iterative_imputation_params_t iterative_params,
                                                   bool apply_binary_threshold,
                                                   double binary_threshold,
                                                   int ikernel,
                                                   double dist_normalization_factor);

std::unique_ptr<std::vector<double>> prop_nbhrs_with_smaller_y(const std::vector<std::vector<int>>& graph,
                                                               const std::vector<double>& y );

extern "C" {
    SEXP S_nn_mean(SEXP Rgraph, SEXP Ry, SEXP Rn_itrs);
    SEXP S_nn_wmean(SEXP Rgraph, SEXP Rw, SEXP Ry, SEXP Rn_itrs, SEXP Repsilon);
    SEXP S_nn_kmean(SEXP Rgraph,
                    SEXP Rd,
                    SEXP Ry,
                    SEXP Rikernel,
                    SEXP Rn_itrs,
                    SEXP Rdist_normalization_factor);
    SEXP S_nn_kmean_cv(SEXP Rgraph,
                       SEXP Rd,
                       SEXP Ry,
                       SEXP Rimputation_method,
                       SEXP Rmax_iterations,
                       SEXP Rconvergence_threshold,
                       SEXP Rapply_binary_threshold,
                       SEXP Rbinary_threshold,
                       SEXP Rikernel,
                       SEXP Rn_itrs,
                       SEXP Rdist_normalization_factor,
                       SEXP Rn_CVs,
                       SEXP Rn_CV_folds,
                       SEXP Repsilon,
                       SEXP Ruse_low_pass_filter,
                       SEXP Rpreserve_local_extrema,
                       SEXP Rmin_plambda,
                       SEXP Rmax_plambda,
                       SEXP Rseed);

    SEXP S_cv_imputation(SEXP Rtest_set,
                         SEXP Rgraph,
                         SEXP Rd,
                         SEXP Ry,
                         SEXP Ry_binary,
                         SEXP Rimputation_method,
                         SEXP Rmax_iterations,
                         SEXP Rconvergence_threshold,
                         SEXP Rapply_binary_threshold,
                         SEXP Rbinary_threshold,
                         SEXP Rikernel,
                         SEXP Rdist_normalization_factor);
}

/**
 * @brief Computes the nearest neighbor mean of a vector over vertices of a graph.
 *
 * This function calculates the mean value of the nearest neighbors for each
 * vertex in the graph. The vertex is included in the mean calculation. The
 * graph is represented as an adjacency matrix, and the values are provided in a
 * vector.
 *
 * @param graph A vector of vectors representing the adjacency matrix of the graph. Each inner vector contains the indices of neighboring vertices.
 * @param y A vector of double values associated with each vertex in the graph.
 * @param n_itrs The number of NN mean iterations.
 *
 * @return A unique pointer to a vector of doubles containing the nearest neighbor mean for each vertex.
 *
 * Usage example:
 * @code
 * std::vector<std::vector<int>> graph = {{1, 2}, {0, 2}, {0, 1}};
 * std::vector<double> y = {1.0, 2.0, 3.0};
 * std::unique_ptr<std::vector<double>> y_nn_mean_uptr = nn_mean(graph, y);
 * std::vector<double> y_nn_mean = std::move(*y_nn_mean_uptr);
 * @endcode
 *
 * The result vector contains the mean value of the nearest neighbors for each vertex in the graph.
 * The memory management of the result vector is handled automatically using std::unique_ptr.
 */
std::unique_ptr<std::vector<double>> nn_mean(const std::vector<std::vector<int>>& graph,
                                             const std::vector<double>& y,
                                             int n_itrs) {
    auto result = std::make_unique<std::vector<double>>(y.size(), 0.0);

    double sum = 0;
    for (size_t i = 0; i < graph.size(); ++i) {
        sum = y[i];
        if (graph[i].empty()) {
            (*result)[i] = sum;
            continue;
        }

        for (int neighbor : graph[i]) {
            sum += y[neighbor];
        }
        (*result)[i] = sum / (graph[i].size() + 1);
    }

    if (n_itrs > 1) {
        for (int itr = 1; itr < n_itrs; itr++) {
            for (size_t i = 0; i < graph.size(); ++i) {
                sum = (*result)[i];
                if (graph[i].empty()) {
                    (*result)[i] = sum;
                    continue;
                }

                for (int neighbor : graph[i]) {
                    sum += (*result)[neighbor];
                }
                (*result)[i] = sum / (graph[i].size() + 1);
            }
        }
    }

    return result;
}


/**
 * @brief Computes the nearest neighbor mean of a numeric vector over vertices of a graph.
 *
 * This function calculates the mean value of the nearest neighbors for each vertex in the graph.
 * The graph is represented as a list of integer vectors (adjacency list), and the values are provided in a numeric vector.
 *
 * @param Rgraph A list of integer vectors representing the adjacency matrix of the graph. Each inner vector contains the indices of neighboring vertices.
 * @param Ry A numeric vector of values associated with each vertex in the graph.
 * @param Rn_itrs The number of NN mean iterations.
 *
 * @return A numeric vector containing the nearest neighbor mean for each vertex.
 *
 * Usage example in R:
 * @code
 * graph <- list(c(1, 2), c(0, 2), c(0, 1))
 * y <- c(1.0, 2.0, 3.0)
 * y_nn_mean <- .Call("S_nn_mean", graph, y)
 * @endcode
 *
 * The result vector contains the mean value of the nearest neighbors for each vertex in the graph.
 */
SEXP S_nn_mean(SEXP Rgraph, SEXP Ry, SEXP Rn_itrs) {

    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));

    int y_length = LENGTH(Ry);
    std::vector<double> y_vec(REAL(Ry), REAL(Ry) + y_length);

    int n_itrs = INTEGER(Rn_itrs)[0];

    // Call the nn_mean function
    std::unique_ptr<std::vector<double>> result = nn_mean(graph, y_vec, n_itrs);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return Rresult;
}


/**
 * Computes the weighted nearest neighbor mean of a vector over vertices of a
 * graph using normalized edge weights.
 *
 * This function calculates the weighted mean value of the nearest neighbors for each
 * vertex in the graph. The vertex itself is included in the mean calculation with the highest weight.
 * Edge weights are normalized to ensure the vertex always has the highest influence.
 *
 * @note 1) The function normalizes edge weights for each vertex such that the
 *       maximum neighbor weight is (1 - epsilon), and the vertex itself always
 *       has a weight of 1.0. This ensures that the vertex always has the
 *       highest influence in its own mean calculation.
 *
 * 2) The way vertext weights is artificially (separately from its neighbors)
 * defined is fixed in nn_kmean function, where the input involves distances of
 * the neighbors to the vertex and then a kernel function is used to compute
 * weights.
 *
 * @param graph A vector of vectors representing the adjacency matrix of the graph.
 *              Each inner vector contains the indices of neighboring vertices.
 * @param w A vector of vectors representing the edge weights.
 *          w[i][j] is the weight of the edge between vertex i and its j-th neighbor.
 * @param y A vector of double values associated with each vertex in the graph.
 * @param n_itrs The number of iterations for computing the nearest neighbor mean.
 * @param epsilon A parameter controlling the normalization of weights. It determines how much more
 *                weight the vertex has compared to its highest-weighted neighbor.
 *                Must be between 0 and 1. Default is 0.1.
 *
 * @return A unique pointer to a vector of doubles containing the nearest neighbor mean for each vertex.
 *
 *
 * Usage example:
 * @code
 * std::vector<std::vector<int>> graph = {{1, 2}, {0, 2}, {0, 1}};
 * std::vector<std::vector<double>> w = {{0.5, 0.3}, {0.5, 0.7}, {0.3, 0.7}};
 * std::vector<double> y = {1.0, 2.0, 3.0};
 * int n_itrs = 3;
 * double epsilon = 0.1;
 * std::unique_ptr<std::vector<double>> y_nn_mean_uptr = nn_mean(graph, w, y, n_itrs, epsilon);
 * std::vector<double> y_nn_mean = std::move(*y_nn_mean_uptr);
 * @endcode
 *
 * The result vector contains the weighted mean value of the nearest neighbors for each vertex in the graph.
 * The memory management of the result vector is handled automatically using std::unique_ptr.
 */
std::unique_ptr<std::vector<double>> nn_wmean(const std::vector<std::vector<int>>& graph,
                                              const std::vector<std::vector<double>>& w,
                                              const std::vector<double>& y,
                                              int n_itrs,
                                              double epsilon = 0.1) {
    auto result = std::make_unique<std::vector<double>>(y.size(), 0.0);
    const double tau = 1.0 - epsilon;

    for (size_t i = 0; i < graph.size(); ++i) {
        double weighted_sum = y[i];  // Include the vertex itself
        double weight_sum = 1.0;     // Weight for the vertex itself

        if (!graph[i].empty()) {
            // Find the maximum weight among neighbors
            double max_weight = *std::max_element(w[i].begin(), w[i].end());

            // Normalize factor
            double normalize_factor = tau / max_weight;

            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                double normalized_weight = w[i][j] * normalize_factor;
                weighted_sum += y[neighbor] * normalized_weight;
                weight_sum += normalized_weight;
            }
        }

        (*result)[i] = weighted_sum / weight_sum;
    }

    if (n_itrs > 1) {
        for (int itr = 1; itr < n_itrs; itr++) {
            auto temp = std::make_unique<std::vector<double>>(y.size(), 0.0);

            for (size_t i = 0; i < graph.size(); ++i) {
                double weighted_sum = (*result)[i];  // Include the vertex itself
                double weight_sum = 1.0;             // Weight for the vertex itself

                if (!graph[i].empty()) {
                    // Find the maximum weight among neighbors
                    double max_weight = *std::max_element(w[i].begin(), w[i].end());

                    // Normalize factor
                    double normalize_factor = tau / max_weight;

                    for (size_t j = 0; j < graph[i].size(); ++j) {
                        int neighbor = graph[i][j];
                        double normalized_weight = w[i][j] * normalize_factor;
                        weighted_sum += (*result)[neighbor] * normalized_weight;
                        weight_sum += normalized_weight;
                    }
                }

                (*temp)[i] = weighted_sum / weight_sum;
            }

            result = std::move(temp);
        }
    }

    return result;
}

/**
 * @brief R interface for the nn_wmean function.
 *
 * This function serves as an interface between R and the C++ nn_wmean function.
 * It converts R objects to C++ types, calls nn_wmean, and converts the result back to an R object.
 *
 * @param Rgraph An R list representing the graph structure.
 * @param Rw An R list of numeric vectors representing edge weights.
 * @param Ry An R numeric vector of values associated with each vertex.
 * @param Rn_itrs An R integer specifying the number of iterations.
 * @param Repsilon An R numeric value for the epsilon parameter.
 *
 * @return An R numeric vector containing the weighted nearest neighbor mean for each vertex.
 */
SEXP S_nn_wmean(SEXP Rgraph, SEXP Rw, SEXP Ry, SEXP Rn_itrs, SEXP Repsilon) {
    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    std::vector<std::vector<double>> weights = std::move(*Rdouble_vectvect_to_Cpp_double_vectvect(Rw));
    int y_length = LENGTH(Ry);
    std::vector<double> y_vec(REAL(Ry), REAL(Ry) + y_length);
    int n_itrs = INTEGER(Rn_itrs)[0];
    double epsilon = REAL(Repsilon)[0];

    // Call the nn_wmean function
    std::unique_ptr<std::vector<double>> result = nn_wmean(graph, weights, y_vec, n_itrs, epsilon);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return Rresult;
}


/**
 * @brief Computes the kernel-weighted nearest neighbor mean of a vector over vertices of a graph.
 *
 * This function calculates the kernel-weighted mean value of the nearest neighbors for
 * each vertex in the graph. The weighting is based on the distances between vertices,
 * transformed by a specified kernel function. The calculation can be iterated multiple times.
 *
 * @param graph   A vector of vectors representing the graph structure. Each inner vector
 *                contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors of edge lengths. edge_lengths[i][j] is the length of the edge between the i-th vertex and its j-th neighbor as listed in graph[i].
 * @param y       A vector of double values associated with each vertex in the graph.
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov),
 *                2 (Triangular),
 *                3 (Truncated Exponential),
 *                4 (Normal).
 * @param n_itrs  The number of iterations for computing the nearest neighbor mean.
 *                If n_itrs > 1, the function will use the results from the previous
 *                iteration as input for the next.
 *
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization
 * range. Default value is 1.01.
 *
 * @return A unique pointer to a vector of doubles containing the kernel-weighted
 *         nearest neighbor mean for each vertex.
 *
 * @note The function includes the vertex itself in the mean calculation with a
 *       distance of 0 and a weight of 1. The memory management of the result
 *       vector is handled automatically using std::unique_ptr.
 *
 * @example
 *   std::vector<std::vector<int>> graph = {{1, 2}, {0, 2}, {0, 1}};
 *   std::vector<std::vector<double>> edge_lengths = {{0.1, 0.2}, {0.1, 0.3}, {0.2, 0.3}};
 *   std::vector<double> y = {1.0, 2.0, 3.0};
 *   int ikernel = 1;  // Epanechnikov kernel
 *   int n_itrs = 3;
 *   auto result = nn_kmean(graph, edge_lengths, y, ikernel, n_itrs);
 */
std::unique_ptr<std::vector<double>> nn_kmean(const std::vector<std::vector<int>>& graph,
                                              const std::vector<std::vector<double>>& edge_lengths,
                                              const std::vector<double>& y,
                                              int ikernel,
                                              int n_itrs,
                                              double dist_normalization_factor = 1.01) {

#define DEBUG__nn_kmean 0

    auto result = std::make_unique<std::vector<double>>(y.size(), 0.0);

    initialize_kernel(ikernel);

     // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

    for (int itr = 0; itr < n_itrs; itr++) {
        auto temp = std::make_unique<std::vector<double>>(y.size(), 0.0);

#if DEBUG__nn_kmean
        Rprintf("itr: %d\n", itr);
#endif

        for (size_t i = 0; i < graph.size(); ++i) {
#if DEBUG__nn_kmean
            Rprintf("\n\n----------\ni: %d\n", (int)i);
            print_vect(graph[i], "graph[i]");
            int m = graph[i].size() + 1;
#endif
            if (!graph[i].empty()) {
                distances[0] = 0;  // Distance to self is 0
                double max_dist = 0.0;
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    distances[j + 1] = edge_lengths[i][j];
                    if (distances[j + 1] > max_dist)
                        max_dist = distances[j + 1];
                }

                if (max_dist == 0) max_dist = 1;  // Avoid division by zero

#if DEBUG__nn_kmean
                print_vect(distances, "distances", m);
                Rprintf("max_dist: %.2f\n", max_dist);
#endif
                max_dist *= dist_normalization_factor;

#if DEBUG__nn_kmean
                Rprintf("max_dist: %.2f\n", max_dist);
#endif

                for (size_t j = 0; j < graph[i].size(); ++j)
                    distances[j + 1] /= max_dist;

#if DEBUG__nn_kmean
                print_vect(distances, "norm distances", m);
#endif
                int n = graph[i].size() + 1;
                kernel_fn(distances.data(), n, kernel_weights.data());
#if DEBUG__nn_kmean
                print_vect(kernel_weights, "kernel_weights", m);
#endif
                double weighted_sum = (itr == 0) ? y[i] : (*result)[i];
                weighted_sum *= kernel_weights[0];
                double weight_sum = kernel_weights[0];
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    int neighbor = graph[i][j];
                    double neighbor_value = (itr == 0) ? y[neighbor] : (*result)[neighbor];
                    weighted_sum += kernel_weights[j + 1] * neighbor_value;
                    weight_sum += kernel_weights[j + 1];
                }

                (*temp)[i] = weighted_sum / weight_sum;
#if DEBUG__nn_kmean
                Rprintf("weighted_sum: %.4f\n", weighted_sum);
                Rprintf("weight_sum: %.4f\n", weight_sum);
                Rprintf("kmean: %.4f\n", (*temp)[i]);
                if (i == 2) error("DEBUGGING");
#endif
            } else {
                (*temp)[i] = (itr == 0) ? y[i] : (*result)[i];
            }
        }

        result = std::move(temp);
    }

    return result;
}

/**
 * @brief R interface for the nn_kmean function.
 *
 * This function serves as an interface between R and the C++ nn_kmean function.
 * It converts R objects to C++ types, calls nn_kmean, and converts the result back to an R object.
 *
 * @param Rgraph An R list representing the graph structure.
 * @param Rd An R list of numeric vectors representing distances between vertices.
 * @param Ry An R numeric vector of values associated with each vertex.
 * @param Rikernel An R integer specifying the kernel function to use.
 * @param Rn_itrs An R integer specifying the number of iterations.
 * @param Rdist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization range.
 *
 * @return An R numeric vector containing the kernel-weighted nearest neighbor mean for each vertex.
 */
SEXP S_nn_kmean(SEXP Rgraph,
                SEXP Rd,
                SEXP Ry,
                SEXP Rikernel,
                SEXP Rn_itrs,
                SEXP Rdist_normalization_factor) {
    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    std::vector<std::vector<double>> distances = std::move(*Rdouble_vectvect_to_Cpp_double_vectvect(Rd));

    int y_length = LENGTH(Ry);
    std::vector<double> y_vec(REAL(Ry), REAL(Ry) + y_length);

    int ikernel = INTEGER(Rikernel)[0];
    int n_itrs = INTEGER(Rn_itrs)[0];
    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];

    // Call the nn_kmean function
    std::unique_ptr<std::vector<double>> result = nn_kmean(graph, distances, y_vec, ikernel, n_itrs, dist_normalization_factor);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);

    return Rresult;
}

/**
 * @brief Computes the cross-validation error of kernel-weighted nearest neighbor mean of a function over vertices of a graph.
 *
 * This function calculates the cross-validation error of the kernel-weighted mean value of the nearest neighbors for each vertex in the graph.
 * The weighting is based on either hop counts or distances between vertices, transformed by a specified kernel function.
 * For binary outcomes, various imputation methods are available to determine the threshold for classification.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of edge lengths. edge_lengths[i][j] is the distance from vertex i to its j-th neighbor as listed in graph[i].
 *          If empty, the hop index strategy will be used for imputation. Otherwise, the kernel distance strategy will be used.
 * @param y A vector of double values associated with each vertex in the graph. If y contains only values 0 and 1, it will be treated as a binary outcome.
 * @param imputation_method Specifies the method for imputing binary outcomes. Options are:
 *        - LOCAL_MEAN_THRESHOLD: Uses the mean of y computed over the training vertices (default).
 *        - NEIGHBORHOOD_MATCHING: Uses a matching method based on local neighborhood statistics.
 *        - ITERATIVE_NEIGHBORHOOD_MATCHING: Uses an iterative version of the NEIGHBORHOOD_MATCHING method.
 *        - SUPPLIED_THRESHOLD: Uses a user-supplied threshold value.
 *        - GLOBAL_MEAN_THRESHOLD: Uses the global mean of y across all vertices.
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param binary_threshold The threshold value to use when imputation_method is SUPPLIED_THRESHOLD. Default is 0.5.
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal). Default is 1.
 * @param n_itrs The number of iterations for computing the nearest neighbor mean.
 *               If n_itrs > 1, the function will use the results from the previous iteration as input for the next. Default is 1.
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization
 * range. Default value is 1.01.
 * @param n_CVs The number of cross-validation rounds. Default is 0.
 * @param n_CV_folds The number of folds in each cross-validation round. Default is 10.
 * @param epsilon A small positive constant for clipping the estimated conditional expectation values of a binary variable. Default is 1e-10.
 * @param use_low_pass_filter Whether to apply a low-pass filter to smooth the results. Default is false.
 * @param preserve_local_extrema If true and using low-pass filter, adjusts weights to preserve local extrema in the smoothing process. Default is false.
 * @param min_plambda Lower bound on the proportion of eigenvectors to use for smoothing. Only used if use_low_pass_filter is true. Default is 0.01.
 * @param max_plambda Upper bound on the proportion of eigenvectors to use for smoothing. Only used if use_low_pass_filter is true. Default is 0.30.
 * @param seed A seed for the random number generator to ensure reproducibility. Default is 0.
 *
 * @return A unique pointer to a vector of cross-validation errors of the kernel-weighted nearest neighbor mean for each vertex.
 *
 * @note The choice of imputation_method affects how binary outcomes are classified:
 *       - LOCAL_MEAN_THRESHOLD adapts to local data distribution in each cross-validation fold.
 *       - NEIGHBORHOOD_MATCHING attempts to preserve local statistical properties of the graph.
 *       - ITERATIVE_NEIGHBORHOOD_MATCHING an iterative version of the NEIGHBORHOOD_MATCHING method.
 *       - SUPPLIED_THRESHOLD allows for domain-specific threshold selection.
 *       - GLOBAL_MEAN_THRESHOLD uses a single threshold based on the entire dataset.
 *
 * @example
 * // Using LOCAL_MEAN_THRESHOLD (default)
 * auto result = nn_kmean_cv(graph, distances, y);
 *
 * // Using SUPPLIED_THRESHOLD with a custom threshold
 * auto result = nn_kmean_cv(graph, distances, y, imputation_method_t::SUPPLIED_THRESHOLD, 0.7);
 *
 * // Using NEIGHBORHOOD_MATCHING
 * auto result = nn_kmean_cv(graph, distances, y, imputation_method_t::NEIGHBORHOOD_MATCHING);
 */
std::unique_ptr<std::vector<double>> nn_kmean_cv(const std::vector<std::vector<int>>& graph,
                                                 const std::vector<std::vector<double>>& edge_lengths,
                                                 const std::vector<double>& y,
                                                 imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                                                 iterative_imputation_params_t iterative_params = {},
                                                 bool apply_binary_threshold = true,
                                                 double binary_threshold = 0.5, // Only used if imputation_method is SUPPLIED_THRESHOLD
                                                 int ikernel = 1,
                                                 int n_itrs = 1,
                                                 double dist_normalization_factor = 1.01,
                                                 int n_CVs = 0,
                                                 int n_CV_folds = 10,
                                                 double epsilon = 1e-10,
                                                 bool use_low_pass_filter = false,
                                                 bool preserve_local_extrema = false,
                                                 double min_plambda = 0.01,
                                                 double max_plambda = 0.30,
                                                 unsigned int seed = 0) {
#define DEBUG__nn_kmean_cv 0

    int n_vertices = y.size();

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && imputation_method == imputation_method_t::GLOBAL_MEAN_THRESHOLD) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    auto cv_error = std::make_unique<std::vector<double>>(n_vertices, 0.0); // Mean Absolute Deviation error - for each vertex the mean over all folds and cross-validation iterations
    std::vector<int> cv_error_count(n_vertices, 0);

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    if (n_CVs > 0) {
        // Creating a set version of the adjacency matrix of the graph
        std::vector<std::set<int>> set_graph(n_vertices);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
        }

        int fold_size = n_vertices / n_CV_folds;

        //
        // The main cross-validation loop
        //
        for (int cv = 0; cv < n_CVs; ++cv) {
#if DEBUG__nn_kmean_cv
            Rprintf("cv: %d\n", cv);
#endif
            // Creating a test set
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    // int vertex = rand() % n_vertices;
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

#if DEBUG__nn_kmean_cv
            print_set(test_set, "test_set");
#endif
            std::vector<double> cv_y = std::move(*cv_imputation(test_set,
                                                                graph,
                                                                edge_lengths,
                                                                y,
                                                                y_binary,
                                                                imputation_method,
                                                                iterative_params,
                                                                apply_binary_threshold,
                                                                binary_threshold,
                                                                ikernel,
                                                                dist_normalization_factor));

            // Estimating the conditional expectation of cv_y
            std::vector<double> Ecv_y = std::move(*nn_kmean(graph, edge_lengths, cv_y, ikernel, n_itrs, dist_normalization_factor));

            if (use_low_pass_filter) {

                 // Initialize weights for each vertex
                std::vector<double> weights(n_vertices, 1.0);

                if (preserve_local_extrema) {
                    // Preserve local extrema by adjusting weights based on local topology
                    // This helps maintain important features in the data while still smoothing

                    // Compute the proportion of neighbors with smaller y values for each vertex
                    std::vector<double> p = std::move(*prop_nbhrs_with_smaller_y(graph, Ecv_y));

                    // Adjust weights to emphasize vertices that are local extrema
                    // Vertices with p close to 0 or 1 are likely local extrema
                    // The weight formula p * (1-p) gives lower weights to these vertices,
                    // preserving their values in the smoothing process
                    for (int vertex = 0; vertex < n_vertices; vertex++)
                        weights[vertex] = p[vertex] * (1 - p[vertex]);
                }

                // Apply spectral graph smoothing
                // This process uses the graph structure to smooth the Ecv_y values
                // It's particularly useful for reducing noise while preserving the underlying structure
                std::unique_ptr<graph_spectral_smoother_result_t> res = graph_spectral_smoother(graph,
                                                                                                edge_lengths,
                                                                                                weights,
                                                                                                Ecv_y,
                                                                                                imputation_method,
                                                                                                iterative_params,
                                                                                                apply_binary_threshold,
                                                                                                binary_threshold,
                                                                                                ikernel,
                                                                                                dist_normalization_factor,
                                                                                                n_CVs,
                                                                                                n_CV_folds,
                                                                                                epsilon,
                                                                                                min_plambda,
                                                                                                max_plambda,
                                                                                                seed);

                // Update Ecv_y with the smoothed values
                Ecv_y.assign(res->y_smoothed.data(), res->y_smoothed.data() + res->y_smoothed.size());
            }

            // Computing error over test vertices
            if (y_binary) {
                // For binary outcomes, compute cross-entropy loss
                // This is more appropriate than mean absolute error for binary classification
                for (const auto& vertex : test_set) {
                    // Clip Ecv_y to avoid log(0) or log(1) issues
                    double clipped_Ecv_y = std::max(epsilon, std::min(1.0 - epsilon, Ecv_y[vertex]));

                    // Compute negative log-likelihood (cross-entropy loss)
                    (*cv_error)[vertex] += -(y[vertex] * log(clipped_Ecv_y) + (1 - y[vertex]) * log(1 - clipped_Ecv_y));
                    cv_error_count[vertex]++;
                }
            } else {
                // For non-binary outcomes, compute mean absolute error
                for (const auto& vertex : test_set) {
                    (*cv_error)[vertex] += std::abs(Ecv_y[vertex] - y[vertex]);
                    cv_error_count[vertex]++;
                }
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex])
            (*cv_error)[vertex] /= cv_error_count[vertex];
    }

    return cv_error;
}


/**
 * @brief R interface for the nn_kmean_cv function.
 *
 * This function serves as an interface between R and the C++ nn_kmean_cv function.
 * It performs cross-validation for kernel-weighted nearest neighbor mean estimation
 * on a graph, with options for different imputation methods and low-pass filtering.
 *
 * @param Rgraph An R list representing the graph structure. Each element of the list
 *               should be an integer vector containing the indices of neighboring vertices.
 * @param Rd An R list of numeric vectors representing distances between vertices.
 *           The structure should match that of Rgraph. If empty, the hop index strategy
 *           will be used for imputation. Otherwise, the kernel distance strategy will be used.
 * @param Ry An R numeric vector of values associated with each vertex in the graph.
 * @param Rimputation_method An R integer specifying the imputation method to use.
 *                           Valid values correspond to the imputation_method_t enum.
 * @param Rmax_iterations The number of iterations in the iterative matching method.
 * @param Rconvergence_threshold The convergence threshold in the iterative matching method.
 * @param Rapply_binary_threshold An R logical indicating whether to apply binary thresholding.
 * @param Rbinary_threshold An R numeric value specifying the threshold for binary classification.
 * @param Rikernel An R integer specifying the kernel function to use. Valid values are:
 *                 1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 * @param Rn_itrs An R integer specifying the number of iterations for computing
 *                the nearest neighbor mean.
 * @param Rdist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization range.
 * @param Rn_CVs An R integer specifying the number of cross-validation rounds.
 * @param Rn_CV_folds An R integer specifying the number of folds in each cross-validation round.
 * @param Repsilon An R numeric value specifying a small positive constant for
 *                 clipping estimated conditional expectation values of a binary variable.
 * @param Ruse_low_pass_filter An R logical indicating whether to apply low-pass filtering.
 * @param Rpreserve_local_extrema An R logical indicating whether to preserve local extrema
 *                                during low-pass filtering.
 * @param Rmin_plambda An R numeric value specifying the minimum proportion of eigenvectors
 *                     to use in spectral smoothing.
 * @param Rmax_plambda An R numeric value specifying the maximum proportion of eigenvectors
 *                     to use in spectral smoothing.
 * @param Rseed An R integer used as a seed for the random number generator.
 *
 * @return An R numeric vector containing the cross-validation errors of the
 *         kernel-weighted nearest neighbor mean for each vertex.
 *
 * @note This function converts R objects to C++ types, calls the nn_kmean_cv function,
 *       and then converts the result back to an R object. It uses PROTECT/UNPROTECT
 *       for proper memory management in R. Ensure that the imputation_method_t enum
 *       is properly defined and accessible in the scope where this function is defined.
 *
 * @see nn_kmean_cv
 */
SEXP S_nn_kmean_cv(SEXP Rgraph,
                   SEXP Rd,
                   SEXP Ry,
                   SEXP Rimputation_method,
                   SEXP Rmax_iterations,
                   SEXP Rconvergence_threshold,
                   SEXP Rapply_binary_threshold,
                   SEXP Rbinary_threshold,
                   SEXP Rikernel,
                   SEXP Rn_itrs,
                   SEXP Rdist_normalization_factor,
                   SEXP Rn_CVs,
                   SEXP Rn_CV_folds,
                   SEXP Repsilon,
                   SEXP Ruse_low_pass_filter,
                   SEXP Rpreserve_local_extrema,
                   SEXP Rmin_plambda,
                   SEXP Rmax_plambda,
                   SEXP Rseed) {
    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    std::vector<std::vector<double>> distances = std::move(*Rdouble_vectvect_to_Cpp_double_vectvect(Rd));
    int y_length = LENGTH(Ry);
    std::vector<double> y_vec(REAL(Ry), REAL(Ry) + y_length);
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);
    int max_iterations = INTEGER(Rmax_iterations)[0];
    double convergence_threshold = REAL(Rconvergence_threshold)[0];
    bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];
    double binary_threshold = REAL(Rbinary_threshold)[0];
    int ikernel = INTEGER(Rikernel)[0];
    int n_itrs = INTEGER(Rn_itrs)[0];
    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];
    int n_CVs = INTEGER(Rn_CVs)[0];
    int n_CV_folds = INTEGER(Rn_CV_folds)[0];
    double epsilon = REAL(Repsilon)[0];
    bool use_low_pass_filter = LOGICAL(Ruse_low_pass_filter)[0];
    bool preserve_local_extrema = LOGICAL(Rpreserve_local_extrema)[0];
    double min_plambda = REAL(Rmin_plambda)[0];
    double max_plambda = REAL(Rmax_plambda)[0];
    unsigned int seed = (unsigned int)INTEGER(Rseed)[0];

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Call the nn_kmean_cv function
    std::unique_ptr<std::vector<double>> result = nn_kmean_cv(graph, distances, y_vec,
                                                              imputation_method,
                                                              iterative_params,
                                                              apply_binary_threshold,
                                                              binary_threshold,
                                                              ikernel,
                                                              n_itrs,
                                                              dist_normalization_factor,
                                                              n_CVs,
                                                              n_CV_folds,
                                                              epsilon,
                                                              use_low_pass_filter,
                                                              preserve_local_extrema,
                                                              min_plambda, max_plambda,
                                                              seed);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return Rresult;
}


/**
 * @brief Performs cross-validation imputation of graph vertex function values.
 *
 * This function imputes the values of a graph vertex function over a test set using only
 * training set vertices. If the distance object is empty hop indices are treated as distanced.
 *
 * @param test_set Set of vertices to be imputed (test vertices).
 * @param graph The graph structure represented as an adjacency list.
 * @param edge_lengths A vector of edge lengths.
 * @param y Vector of observed values for the graph vertex function.
 * @param y_binary Boolean indicating if y is binary (true) or continuous (false).
 * @param imputation_method Method used for imputation. Options include:
 *                          - LOCAL_MEAN_THRESHOLD
 *                          - GLOBAL_MEAN_THRESHOLD
 *                          - NEIGHBORHOOD_MATCHING
 *                          - ITERATIVE_NEIGHBORHOOD_MATCHING (New)
 *                          - SUPPLIED_THRESHOLD
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param apply_binary_threshold Whether to apply binary threshold for binary y (default: true).
 * @param binary_threshold Threshold for binary classification (default: 0.5).
 * @param ikernel Kernel type for distance weighting (default: 1).
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization range.
 *
 * @return std::unique_ptr<std::vector<double>> Pointer to vector of imputed y values.
 *
 * @details The function first constructs a test_graph that contains the
 * relationships between test vertices and their nearest training vertices. For
 * hop index strategy, it stores hop counts; for kernel distance strategy, it
 * stores distances. Hop counts are treated as distances.
 *
 * For binary y with NEIGHBORHOOD_MATCHING imputation method, it uses a specialized
 * algorithm that computes expected y values for training vertices and uses these to
 * impute test vertices based on their neighborhoods.
 *
 * For other cases, it performs weighted averaging of neighboring training vertices'
 * y values, using kernel-weighted distances.
 *
 * For binary y, it can apply a threshold to convert continuous imputed values to binary.
 *
 * @note The function assumes that the graph is connected and that each test vertex
 * has at least one path to a training vertex.
 */
std::unique_ptr<std::vector<double>> cv_imputation(const std::set<int>& test_set,
                                                   const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& edge_lengths,
                                                   const std::vector<double>& y,
                                                   bool y_binary,
                                                   imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                                                   iterative_imputation_params_t iterative_params = {},
                                                   bool apply_binary_threshold = true,
                                                   double binary_threshold = 0.5, // Only used if imputation_method is SUPPLIED_THRESHOLD
                                                   int ikernel = 1,
                                                   double dist_normalization_factor = 1.01) {
    int n_vertices = y.size();
    bool use_hop_index = edge_lengths.empty();
    initialize_kernel(ikernel);
    std::vector<double> cv_y(y); // Initializing cv_y to y

    // Creating a training set as the complement of the test set
    std::vector<int> all_vertices(n_vertices);
    std::iota(all_vertices.begin(), all_vertices.end(), 0);
    std::set<int> training_set;
    std::set_difference(all_vertices.begin(), all_vertices.end(), test_set.begin(), test_set.end(), std::inserter(training_set, training_set.begin()));

    //
    // Constructing test_graph
    //
    std::vector<std::map<int, double>> test_graph(n_vertices);

    if (use_hop_index) {
        // For hop index strategy, we store hop counts and treat it as edge length or distance
        for (const auto& test_vertex : test_set) {
            std::queue<std::pair<int, int>> bfs_queue;
            std::unordered_set<int> visited;

            bfs_queue.push({test_vertex, 0});
            visited.insert(test_vertex);

            while (!bfs_queue.empty()) {
                auto [curr_vertex, hop_count] = bfs_queue.front();
                bfs_queue.pop();

                if (training_set.find(curr_vertex) != training_set.end()) {
                    test_graph[test_vertex][curr_vertex] = hop_count;
                    continue;
                }

                for (const auto& neighbor : graph[curr_vertex]) {
                    if (visited.find(neighbor) == visited.end()) {
                        bfs_queue.push({neighbor, hop_count + 1});
                        visited.insert(neighbor);
                    }
                }
            }
        }
    } else {
        // For kernel distance strategy, we store distances between the given vertex and its neighbors
        for (const auto& test_vertex : test_set) {
            std::priority_queue<std::pair<double, int>,
                                std::vector<std::pair<double, int>>,
                                std::greater<std::pair<double, int>>> pq;
            std::vector<double> distances(n_vertices, std::numeric_limits<double>::max());

            pq.push({0.0, test_vertex});
            distances[test_vertex] = 0.0;

            while (!pq.empty()) {
                double dist = pq.top().first;
                int curr_vertex = pq.top().second;
                pq.pop();

                if (training_set.find(curr_vertex) != training_set.end()) {
                    test_graph[test_vertex][curr_vertex] = dist;
                    continue;  // We've found a training vertex, no need to explore further
                }

                if (dist > distances[curr_vertex]) continue;

                for (size_t i = 0; i < graph[curr_vertex].size(); ++i) {
                    int neighbor = graph[curr_vertex][i];
                    double new_dist = dist + edge_lengths[curr_vertex][i];

                    if (new_dist < distances[neighbor]) {
                        distances[neighbor] = new_dist;
                        pq.push({new_dist, neighbor});
                    }
                }
            }
        }
    }

    if (y_binary && imputation_method == imputation_method_t::NEIGHBORHOOD_MATCHING) {

        /**
         * @brief Lambda function to impute binary values using distance-weighted neighborhood matching.
         *
         * This function implements an enhanced neighborhood matching method for binary value imputation,
         * incorporating distance information to weight the influence of neighboring vertices. It first
         * computes the expected y values (Ey) for training vertices, then uses these along with
         * distance-based kernel weights to impute values for test vertices.
         *
         * @param graph The original graph structure.
         * @param test_graph A pre-computed graph structure containing relevant training neighbors
         *                   for each test vertex, including their distances.
         * @param test_set The set of vertices to be imputed (test vertices).
         * @param training_set The set of vertices with known values (training vertices).
         * @param kernel_function Pointer to the kernel function used for distance-based weighting.
         *
         * @note This function captures y, cv_y, and dist_normalization_factor by reference.
         *
         * @details The imputation process works as follows:
         *   1. For each training vertex, compute Ey as the mean y of its training neighbors.
         *   2. For each test vertex:
         *      a. Normalize distances to neighbors using max distance and dist_normalization_factor.
         *      b. Apply kernel function to normalized distances to get weights.
         *      c. Compute weighted averages of y and Ey for neighboring vertices.
         *      d. Calculate mean_y_with_0 and mean_y_with_1, incorporating the kernel weight at zero distance.
         *      e. Assign 1 if mean_y_with_1 is closer to mean_Ey, otherwise assign 0.
         *
         * @note This function assumes that test_graph is constructed such that each test vertex
         *       has at least one neighbor from the training set, possibly at hop distances > 1.
         * @note The function modifies the cv_y vector directly, which is captured by reference.
         */
        auto impute_binary_value_by_neighborhood_matching = [&y,&cv_y,&dist_normalization_factor](const std::vector<std::vector<int>>& graph,
                                                                                                  const std::vector<std::map<int, double>>& test_graph,
                                                                                                  const std::set<int>& test_set,
                                                                                                  const std::set<int>& training_set,
                                                                                                  void (*kernel_function)(const double*, int, double*)) {
            // Computing Ey for every vertex of the training set using only training vertices
            std::vector<double> Ey(y.size(), 0.0);
            for (const auto& vertex : training_set) {
                double sum = 0;
                int neighbor_counter = 0;
                for (const auto& neighbor : graph[vertex]) {
                    if (training_set.find(neighbor) != training_set.end()) {
                        sum += y[neighbor];
                        neighbor_counter++;
                    }
                }
                Ey[vertex] = (neighbor_counter > 0) ? sum / neighbor_counter : y[vertex];
            }
            // For every vertex of the test set, set cv_y to 0 or 1
            for (const auto& vertex : test_set) {
                double weighted_sum_y = 0.0, weighted_sum_Ey = 0.0;
                double total_weight  = 0.0;
                double max_dist      = 0.0;
                for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                    if (distance > max_dist)
                        max_dist = distance;
                }
                if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                max_dist *= dist_normalization_factor;

                for (const auto& [neighbor, distance] : test_graph[vertex]) {
                    double normalized_distance = distance / max_dist;
                    double kernel_weight;
                    kernel_function(&normalized_distance, 1, &kernel_weight);
                    weighted_sum_y += kernel_weight * y[neighbor];
                    weighted_sum_Ey += kernel_weight * Ey[neighbor];
                    total_weight += kernel_weight;
                }
                double mean_Ey = weighted_sum_Ey / total_weight;
                double zero = 0;
                double kernel_weight_at_zero;
                kernel_function(&zero, 1, &kernel_weight_at_zero);
                double mean_y_with_0 = weighted_sum_y / (total_weight + kernel_weight_at_zero);
                double mean_y_with_1 = (weighted_sum_y + kernel_weight_at_zero) / (total_weight + kernel_weight_at_zero);
                cv_y[vertex] = (std::abs(mean_y_with_1 - mean_Ey) < std::abs(mean_y_with_0 - mean_Ey)) ? 1.0 : 0.0;
            }
        };

        impute_binary_value_by_neighborhood_matching(graph, test_graph, test_set, training_set, kernel_fn);

    } else if (!y_binary && imputation_method == imputation_method_t::NEIGHBORHOOD_MATCHING) {

        /**
         * @brief Lambda function to impute continuous values using distance-weighted neighborhood matching.
         *
         * This function implements an enhanced neighborhood matching method for continuous value imputation,
         * incorporating distance information to weight the influence of neighboring vertices. It first
         * computes the expected y values (Ey) for training vertices, then uses these along with
         * distance-based kernel weights to impute values for test vertices.
         *
         * @param graph The original graph structure.
         * @param test_graph A pre-computed graph structure containing relevant training neighbors
         *                   for each test vertex, including their distances.
         * @param test_set The set of vertices to be imputed (test vertices).
         * @param training_set The set of vertices with known values (training vertices).
         * @param kernel_function Pointer to the kernel function used for distance-based weighting.
         *
         * @note This function captures y, cv_y, and dist_normalization_factor by reference.
         *
         * @details The imputation process works as follows:
         *   1. For each training vertex, compute Ey as the mean y of its training neighbors.
         *   2. For each test vertex:
         *      a. Normalize distances to neighbors using max distance and dist_normalization_factor.
         *      b. Apply kernel function to normalized distances to get weights.
         *      c. Compute weighted averages of y and Ey for neighboring vertices.
         *      d. Calculate an adjustment factor as half the difference between mean_Ey and mean_y.
         *      e. Impute the value as mean_y plus the adjustment factor.
         *      f. Optionally clip the imputed value to the range of original y values.
         *
         * @note This function assumes that test_graph is constructed such that each test vertex
         *       has at least one neighbor from the training set, possibly at hop distances > 1.
         * @note The function modifies the cv_y vector directly, which is captured by reference.
         * @note The imputation method attempts to balance preserving local characteristics with
         *       adjusting towards expected neighborhood characteristics.
         */
        auto impute_continuous_value_by_neighborhood_matching = [&y,&cv_y,&dist_normalization_factor](
            const std::vector<std::vector<int>>& graph,
            const std::vector<std::map<int, double>>& test_graph,
            const std::set<int>& test_set,
            const std::set<int>& training_set,
            void (*kernel_function)(const double*, int, double*)) {
            // Computing Ey for every vertex of the training set using only training vertices
            std::vector<double> Ey(y.size(), 0.0);
            for (const auto& vertex : training_set) {
                double sum = 0;
                int neighbor_counter = 0;
                for (const auto& neighbor : graph[vertex]) {
                    if (training_set.find(neighbor) != training_set.end()) {
                        sum += y[neighbor];
                        neighbor_counter++;
                    }
                }
                Ey[vertex] = (neighbor_counter > 0) ? sum / neighbor_counter : y[vertex];
            }
            // For every vertex of the test set, impute a value
            for (const auto& vertex : test_set) {
                double weighted_sum_y = 0.0, weighted_sum_Ey = 0.0;
                double total_weight = 0.0;
                double max_dist = 0.0;
                for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                    if (distance > max_dist)
                        max_dist = distance;
                }
                if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                max_dist *= dist_normalization_factor;

                for (const auto& [neighbor, distance] : test_graph[vertex]) {
                    double normalized_distance = distance / max_dist;
                    double kernel_weight;
                    kernel_function(&normalized_distance, 1, &kernel_weight);
                    weighted_sum_y += kernel_weight * y[neighbor];
                    weighted_sum_Ey += kernel_weight * Ey[neighbor];
                    total_weight += kernel_weight;
                }
                double mean_y = weighted_sum_y / total_weight;
                double mean_Ey = weighted_sum_Ey / total_weight;

                // Adjust the imputed value to minimize the difference between local mean and neighborhood mean
                double adjustment_factor = (mean_Ey - mean_y) / 2.0; // Use half the difference as adjustment
                cv_y[vertex] = mean_y + adjustment_factor;

                // Optional: Clip the imputed value to the range of original y values
                double min_y = *std::min_element(y.begin(), y.end());
                double max_y = *std::max_element(y.begin(), y.end());
                cv_y[vertex] = std::max(min_y, std::min(max_y, cv_y[vertex]));
            }
        };

        impute_continuous_value_by_neighborhood_matching(graph, test_graph, test_set, training_set, kernel_fn);

    } else if (y_binary && imputation_method == imputation_method_t::ITERATIVE_NEIGHBORHOOD_MATCHING) {

        /**
         * @brief Iteratively imputes binary values using distance-weighted neighborhood matching.
         *
         * This function implements an iterative refinement of the neighborhood matching method
         * for binary value imputation, incorporating distance information to weight the influence
         * of neighboring vertices. It repeatedly updates the estimates for test vertices
         * based on the current estimates of their neighbors, using distance-based kernel weights.
         *
         * @param graph The original graph structure.
         * @param test_graph A pre-computed graph structure containing relevant training neighbors,
         *                   including their distances.
         * @param test_set The set of vertices to be imputed (test vertices).
         * @param training_set The set of vertices with known values (training vertices).
         * @param kernel_function Pointer to the kernel function used for distance-based weighting.
         * @param max_iterations Maximum number of iterations for refinement.
         * @param convergence_threshold Threshold for considering the process converged.
         *
         * @note This function captures y, cv_y, and dist_normalization_factor by reference.
         *
         * @details The imputation process works as follows:
         *   1. For each iteration (up to max_iterations):
         *      a. Compute Ey for all vertices using current estimates.
         *      b. For each test vertex:
         *         i. Normalize distances to neighbors using max distance and dist_normalization_factor.
         *         ii. Apply kernel function to normalized distances to get weights.
         *         iii. Compute weighted averages of y and Ey for neighboring vertices.
         *         iv. Calculate mean_y_with_0 and mean_y_with_1, incorporating the kernel weight at zero distance.
         *         v. Assign 1 if mean_y_with_1 is closer to mean_Ey, otherwise assign 0.
         *      c. Check for convergence based on maximum change in estimates.
         *
         * @note This function assumes that test_graph is constructed such that each test vertex
         *       has at least one neighbor from the training set, possibly at hop distances > 1.
         * @note The function modifies the cv_y vector directly, which is captured by reference.
         */
        auto impute_binary_value_by_iterative_neighborhood_matching =
            [&y,&cv_y,&dist_normalization_factor](const std::vector<std::vector<int>>& graph,
                                                  const std::vector<std::map<int, double>>& test_graph,
                                                  const std::set<int>& test_set,
                                                  const std::set<int>& training_set,
                                                  void (*kernel_function)(const double*, int, double*),
                                                  int max_iterations = 10,
                                                  double convergence_threshold = 1e-6) {
                std::vector<double> prev_cv_y(y.size());
                for (int iteration = 0; iteration < max_iterations; ++iteration) {
                    prev_cv_y = cv_y;
                    // Compute Ey for all vertices using current estimates
                    std::vector<double> Ey(y.size(), 0.0);
                    for (int vertex = 0; vertex < static_cast<int>(y.size()); ++vertex) {
                        double sum = 0;
                        int neighbor_counter = 0;
                        for (const auto& neighbor : graph[vertex]) {
                            sum += cv_y[neighbor];
                            neighbor_counter++;
                        }
                        Ey[vertex] = (neighbor_counter > 0) ? sum / neighbor_counter : cv_y[vertex];
                    }
                    // Update estimates for test vertices
                    for (const auto& vertex : test_set) {
                        double weighted_sum_y = 0.0, weighted_sum_Ey = 0.0;
                        double total_weight = 0.0;
                        double max_dist = 0.0;
                        for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                            if (distance > max_dist)
                                max_dist = distance;
                        }
                        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                        max_dist *= dist_normalization_factor;

                        for (const auto& [neighbor, distance] : test_graph[vertex]) {
                            double normalized_distance = distance / max_dist;
                            double kernel_weight;
                            kernel_function(&normalized_distance, 1, &kernel_weight);
                            weighted_sum_y += kernel_weight * cv_y[neighbor];
                            weighted_sum_Ey += kernel_weight * Ey[neighbor];
                            total_weight += kernel_weight;
                        }
                        double mean_Ey = weighted_sum_Ey / total_weight;
                        double zero = 0;
                        double kernel_weight_at_zero;
                        kernel_function(&zero, 1, &kernel_weight_at_zero);
                        double mean_y_with_0 = weighted_sum_y / (total_weight + kernel_weight_at_zero);
                        double mean_y_with_1 = (weighted_sum_y + kernel_weight_at_zero) / (total_weight + kernel_weight_at_zero);
                        cv_y[vertex] = (std::abs(mean_y_with_1 - mean_Ey) < std::abs(mean_y_with_0 - mean_Ey)) ? 1.0 : 0.0;
                    }
                    // Check for convergence
                    double max_change = 0.0;
                    for (const auto& vertex : test_set) {
                        max_change = std::max(max_change, std::abs(cv_y[vertex] - prev_cv_y[vertex]));
                    }
                    if (max_change < convergence_threshold) {
                        break;
                    }
                }
            };

        impute_binary_value_by_iterative_neighborhood_matching(graph, test_graph, test_set, training_set,
                                                               kernel_fn,
                                                               iterative_params.max_iterations,
                                                               iterative_params.convergence_threshold);

    } else if (!y_binary && imputation_method == imputation_method_t::ITERATIVE_NEIGHBORHOOD_MATCHING) {

        /**
         * @brief Iteratively imputes continuous values using distance-weighted neighborhood matching.
         *
         * This function implements an iterative refinement of the neighborhood matching method
         * for continuous value imputation, incorporating distance information to weight the influence
         * of neighboring vertices. It repeatedly updates the estimates for test vertices
         * based on the current estimates of their neighbors, using distance-based kernel weights
         * and adjusting for local and neighborhood means.
         *
         * @param graph The original graph structure.
         * @param test_graph A pre-computed graph structure containing relevant training neighbors,
         *                   including their distances.
         * @param test_set The set of vertices to be imputed (test vertices).
         * @param training_set The set of vertices with known values (training vertices).
         * @param kernel_function Pointer to the kernel function used for distance-based weighting.
         * @param max_iterations Maximum number of iterations for refinement.
         * @param convergence_threshold Threshold for considering the process converged.
         *
         * @note This function captures y, cv_y, and dist_normalization_factor by reference.
         *
         * @details The imputation process works as follows:
         *   1. For each iteration (up to max_iterations):
         *      a. Compute Ey for all vertices using current estimates.
         *      b. For each test vertex:
         *         i. Normalize distances to neighbors using max distance and dist_normalization_factor.
         *         ii. Apply kernel function to normalized distances to get weights.
         *         iii. Compute weighted averages of y and Ey for neighboring vertices.
         *         iv. Calculate an adjustment factor as half the difference between mean_Ey and mean_y.
         *         v. Impute the value as mean_y plus the adjustment factor, clipped to the range of original y values.
         *      c. Check for convergence based on maximum change in estimates.
         *
         * @note This function assumes that test_graph is constructed such that each test vertex
         *       has at least one neighbor from the training set, possibly at hop distances > 1.
         * @note The function modifies the cv_y vector directly, which is captured by reference.
         * @note The imputation method attempts to balance preserving local characteristics with
         *       adjusting towards expected neighborhood characteristics, while considering distance-based influence.
         */
        auto impute_continuous_value_by_iterative_neighborhood_matching =
            [&y,&cv_y,&dist_normalization_factor](const std::vector<std::vector<int>>& graph,
                                                  const std::vector<std::map<int, double>>& test_graph,
                                                  const std::set<int>& test_set,
                                                  const std::set<int>& training_set,
                                                  void (*kernel_function)(const double*, int, double*),
                                                  int max_iterations = 10,
                                                  double convergence_threshold = 1e-6) {
                std::vector<double> prev_cv_y(y.size());
                double min_y = *std::min_element(y.begin(), y.end());
                double max_y = *std::max_element(y.begin(), y.end());
                for (int iteration = 0; iteration < max_iterations; ++iteration) {
                    prev_cv_y = cv_y;
                    // Compute Ey for all vertices using current estimates
                    std::vector<double> Ey(y.size(), 0.0);
                    for (int vertex = 0; vertex < static_cast<int>(y.size()); ++vertex) {
                        double sum = 0;
                        int neighbor_counter = 0;
                        for (const auto& neighbor : graph[vertex]) {
                            sum += cv_y[neighbor];
                            neighbor_counter++;
                        }
                        Ey[vertex] = (neighbor_counter > 0) ? sum / neighbor_counter : cv_y[vertex];
                    }
                    // Update estimates for test vertices
                    for (const auto& vertex : test_set) {
                        double weighted_sum_y = 0.0, weighted_sum_Ey = 0.0;
                        double total_weight = 0.0;
                        double max_dist = 0.0;
                        for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                            if (distance > max_dist)
                                max_dist = distance;
                        }
                        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                        max_dist *= dist_normalization_factor;

                        for (const auto& [neighbor, distance] : test_graph[vertex]) {
                            double normalized_distance = distance / max_dist;
                            double kernel_weight;
                            kernel_function(&normalized_distance, 1, &kernel_weight);
                            weighted_sum_y += kernel_weight * cv_y[neighbor];
                            weighted_sum_Ey += kernel_weight * Ey[neighbor];
                            total_weight += kernel_weight;
                        }
                        double mean_y = weighted_sum_y / total_weight;
                        double mean_Ey = weighted_sum_Ey / total_weight;
                        double adjustment_factor = (mean_Ey - mean_y) / 2.0;
                        cv_y[vertex] = std::max(min_y, std::min(max_y, mean_y + adjustment_factor));
                    }
                    // Check for convergence
                    double max_change = 0.0;
                    for (const auto& vertex : test_set) {
                        max_change = std::max(max_change, std::abs(cv_y[vertex] - prev_cv_y[vertex]));
                    }
                    if (max_change < convergence_threshold) {
                        break;
                    }
                }
            };

        impute_continuous_value_by_iterative_neighborhood_matching(graph, test_graph, test_set, training_set,
                                                                   kernel_fn,
                                                                   iterative_params.max_iterations,
                                                                   iterative_params.convergence_threshold);
    } else { // remaining imputation methods

        /**
         * @brief Imputes values for test vertices using kernel-based weighting.
         *
         * This function imputes values for vertices in the test set by computing
         * a weighted average of their neighbors' values. The weight is determined
         * by applying a kernel function to the normalized distances between vertices.
         *
         * @param test_set Set of vertices to impute values for.
         * @param test_graph Graph structure containing distances between vertices.
         * @param y Vector of original values for all vertices.
         * @param kernel_function Pointer to the kernel function to be used for weighting.
         *
         * @note This function modifies the cv_y vector captured by reference.
         */
        auto kernel_imputation = [&cv_y,&dist_normalization_factor](
            const std::set<int>& test_set,
            const std::vector<std::map<int, double>>& test_graph,
            const std::vector<double>& y,
            void (*kernel_function)(const double*, int, double*)) {
            for (const auto& vertex : test_set) {
                double weighted_sum  = 0.0;
                double total_weight  = 0.0;
                double max_dist      = 0.0;
                for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                    if (distance > max_dist)
                        max_dist = distance;
                }
                if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                max_dist *= dist_normalization_factor;
                for (const auto& [neighbor_index, distance] : test_graph[vertex]) {
                    double normalized_distance = distance / max_dist;
                    double kernel_weight;
                    kernel_function(&normalized_distance, 1, &kernel_weight);
                    weighted_sum += kernel_weight * y[neighbor_index];
                    total_weight += kernel_weight;
                }
                cv_y[vertex] = (total_weight > 0) ? weighted_sum / total_weight : y[vertex];
            }
        };

        kernel_imputation(test_set, test_graph, y, kernel_fn);  // Modifies cv_y using kernel-based imputation

        if (y_binary && apply_binary_threshold) {
            if (imputation_method == imputation_method_t::LOCAL_MEAN_THRESHOLD) {

                auto compute_local_mean_threshold = [&y](const std::set<int>& training_set) {
                    double sum = 0.0;
                    for (int idx : training_set) {
                        sum += y[idx];
                    }
                    return sum / training_set.size();
                };

                binary_threshold = compute_local_mean_threshold(training_set);
            }

            for (const auto& vertex : test_set)
                cv_y[vertex] = (cv_y[vertex] > binary_threshold) ? 1.0 : 0.0;
        }
    }

    std::unique_ptr<std::vector<double>> cv_y_uptr = std::make_unique<std::vector<double>>(std::move(cv_y));

    return cv_y_uptr;
}

/**
 * @brief R interface for cross-validation imputation on graphs
 *
 * This function provides an R interface to perform cross-validation imputation
 * on graph-structured data. It supports various imputation methods including
 * local mean, neighborhood matching, and iterative approaches for both binary
 * and continuous data.
 *
 * @param Rtest_set An integer vector of indices for the test set vertices (1-based)
 * @param Rgraph A list of integer vectors representing the graph structure. NOTE: It is assumed that the components of Rgraph are 0-based.
 * @param Rd A list of numeric vectors representing distances (empty list if using hop index)
 * @param Ry A numeric vector of original vertex values
 * @param Ry_binary A logical value indicating whether the data is binary
 * @param Rimputation_method An integer representing the imputation method:
 *        0: LOCAL_MEAN, 1: LOCAL_MEDIAN, 2: LOCAL_ABSOLUTE_MEDIAN,
 *        3: NEIGHBORHOOD_MATCHING, 4: ITERATIVE_NEIGHBORHOOD_MATCHING,
 *        5: LOCAL_MEAN_THRESHOLD, 6: SUPPLIED_THRESHOLD
 * @param Rmax_iterations An integer for the maximum number of iterations (for iterative methods)
 * @param Rconvergence_threshold A numeric value for the convergence threshold (for iterative methods)
 * @param Rapply_binary_threshold A logical value indicating whether to apply binary thresholding
 * @param Rbinary_threshold A numeric value for the binary threshold
 * @param Rikernel An integer representing the kernel type:
 *        0: Box, 1: Triangular, 2: Epanechnikov, 3: Gaussian
 * @param Rdist_normalization_factor A numeric value for distance normalization
 *
 * @return An R numeric vector containing the imputed values for the test set
 *
 * @details This function converts R objects to C++ types, calls the C++ implementation
 * of cv_imputation, and then converts the result back to an R vector. It handles
 * memory management for R objects and ensures proper conversion between R and C++ data types.
 *
 * The function supports various imputation strategies and can handle both binary and
 * continuous data. It allows for flexible graph representations and distance metrics,
 * making it suitable for a wide range of graph-structured imputation tasks.
 *
 * @note This function assumes that the input data is properly formatted and valid.
 * It does not perform extensive error checking on the inputs.
 *
 * @seealso cv_imputation for the underlying C++ implementation
 */
SEXP S_cv_imputation(SEXP Rtest_set,
                     SEXP Rgraph,
                     SEXP Rd,
                     SEXP Ry,
                     SEXP Ry_binary,
                     SEXP Rimputation_method,
                     SEXP Rmax_iterations,
                     SEXP Rconvergence_threshold,
                     SEXP Rapply_binary_threshold,
                     SEXP Rbinary_threshold,
                     SEXP Rikernel,
                     SEXP Rdist_normalization_factor) {

    // Protect R objects from garbage collection
    int nprot = 0;
    PROTECT(Rtest_set = coerceVector(Rtest_set, INTSXP)); nprot++;
    PROTECT(Rgraph = coerceVector(Rgraph, VECSXP)); nprot++;
    PROTECT(Rd = coerceVector(Rd, VECSXP)); nprot++;
    PROTECT(Ry = coerceVector(Ry, REALSXP));nprot++;
    PROTECT(Ry_binary = coerceVector(Ry_binary, LGLSXP)); nprot++;
    PROTECT(Rimputation_method = coerceVector(Rimputation_method, INTSXP)); nprot++;
    PROTECT(Rmax_iterations = coerceVector(Rmax_iterations, INTSXP)); nprot++;
    PROTECT(Rconvergence_threshold = coerceVector(Rconvergence_threshold, REALSXP)); nprot++;
    PROTECT(Rapply_binary_threshold = coerceVector(Rapply_binary_threshold, LGLSXP)); nprot++;
    PROTECT(Rbinary_threshold = coerceVector(Rbinary_threshold, REALSXP)); nprot++;
    PROTECT(Rikernel = coerceVector(Rikernel, INTSXP)); nprot++;
    PROTECT(Rdist_normalization_factor = coerceVector(Rdist_normalization_factor, REALSXP)); nprot++;

    // Convert R objects to C++ types
    std::set<int> test_set;
    for (int i = 0; i < LENGTH(Rtest_set); ++i) {
        test_set.insert(INTEGER(Rtest_set)[i] - 1);  // R indices are 1-based
    }

    std::vector<std::vector<int>> graph(LENGTH(Rgraph));
    for (int i = 0; i < LENGTH(Rgraph); ++i) {
        SEXP Rneighbors = VECTOR_ELT(Rgraph, i);
        for (int j = 0; j < LENGTH(Rneighbors); ++j) {
            graph[i].push_back(INTEGER(Rneighbors)[j]);
        }
    }

    std::vector<std::vector<double>> edge_lengths;
    if (LENGTH(Rd) > 0) {
        edge_lengths.resize(LENGTH(Rd));
        for (int i = 0; i < LENGTH(Rd); ++i) {
            SEXP Rdists = VECTOR_ELT(Rd, i);
            edge_lengths[i].assign(REAL(Rdists), REAL(Rdists) + LENGTH(Rdists));
        }
    }

    std::vector<double> y(REAL(Ry), REAL(Ry) + LENGTH(Ry));
    bool y_binary = LOGICAL(Ry_binary)[0];
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);
    iterative_imputation_params_t iterative_params = {
        INTEGER(Rmax_iterations)[0],
        REAL(Rconvergence_threshold)[0]
    };
    bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];
    double binary_threshold = REAL(Rbinary_threshold)[0];
    int ikernel = INTEGER(Rikernel)[0];
    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];

    // Call the C++ function
    auto result = cv_imputation(test_set,
                                graph,
                                edge_lengths,
                                y,
                                y_binary,
                                imputation_method,
                                iterative_params,
                                apply_binary_threshold,
                                binary_threshold,
                                ikernel,
                                dist_normalization_factor);

    // Convert the result to an R vector
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size())); nprot++;
    std::copy(result->begin(), result->end(), REAL(Rresult));

    // Unprotect R objects
    UNPROTECT(nprot);

    return Rresult;
}

/**
 * @brief Computes edge weights for a graph given edge lengths and a kernel function.
 *
 * This function calculates weights for each vertex and its edges in the graph
 * based on the provided edge lengths and a specified kernel function. It normalizes
 * the distances, applies the kernel function to compute initial weights, and then
 * normalizes these weights so that they sum to 1 for each vertex (including the
 * vertex's self-weight).
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector
 *              contains the indices of neighboring vertices for a given vertex.
 * @param d A vector of vectors containing the edge lengths. Each inner vector corresponds
 *          to a vertex and contains the lengths of edges to its neighbors.
 * @param ikernel An integer specifying the kernel function to use. Default is 1.
 * @param dist_normalization_factor A scaling factor applied to the maximum distance between
 *                                  a vertex and its neighbors. This ensures non-zero weights
 *                                  even when all distances are equal, by slightly increasing
 *                                  the normalization range. Default value is 1.01.
 *
 * @return A unique pointer to a vector of vectors containing the computed weights.
 *         The outer vector corresponds to vertices, and each inner vector contains
 *         weights for the vertex and its edges. The first weight in each inner
 *         vector is the self-weight of the vertex, followed by weights for its
 *         neighboring edges. All weights are normalized so that they sum to 1 for each vertex.
 *
 * @note This function initializes the kernel function based on the provided ikernel parameter.
 * @note If the maximum distance for a vertex is 0, it's set to 1 to avoid division by zero.
 * @note The function assumes that the graph and d vectors have the same number of vertices
 *       and that the edge information in both is consistent.
 *
 * @see initialize_kernel
 */
std::unique_ptr<std::vector<std::vector<double>>> dist_to_weights(const std::vector<std::vector<int>>& graph,
                                                                  const std::vector<std::vector<double>>& edge_lengths,
                                                                  int ikernel = 1,
                                                                  double dist_normalization_factor = 1.01) {
    int n_vertices = graph.size();
    initialize_kernel(ikernel);
    std::vector<std::vector<double>> weights(n_vertices);
    double zero = 0;

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        double max_dist = 0.0;
        for (const auto& edge_length : edge_lengths[vertex]) {
            if (edge_length > max_dist)
                max_dist = edge_length;
        }
        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
        max_dist *= dist_normalization_factor;

        std::vector<double> vertex_weights;
        vertex_weights.reserve(graph[vertex].size() + 1);
        double kernel_weight;
        kernel_fn(&zero, 1, &kernel_weight);
        vertex_weights.push_back(kernel_weight);

        for (const auto& edge_length : edge_lengths[vertex]) {
            double normalized_weight = edge_length / max_dist;
            kernel_fn(&normalized_weight, 1, &kernel_weight);
            vertex_weights.push_back(kernel_weight);
        }

        double total_weight = 0;
        for (const auto& weight : vertex_weights)
            total_weight += weight;
        for (size_t i = 0; i < vertex_weights.size(); i++)
            vertex_weights[i] /= total_weight;

        weights[vertex] = std::move(vertex_weights);
    }

    return std::make_unique<std::vector<std::vector<double>>>(std::move(weights));
}
