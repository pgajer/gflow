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
#include <cmath>

#include "msr2.h"
#include "cpp_utils.h"
#include "SEXP_cpp_conversion_utils.h"
#include "graph_diffusion_smoother.h"
#include "stats_utils.h"
#include "kernels.h"

void scale_to_range(std::vector<double>& x, double ymin, double ymax);
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

    SEXP S_nn_kmean(SEXP s_graph,
                    SEXP s_edge_lengths,
                    SEXP s_y,
                    SEXP s_ikernel,
                    SEXP s_n_itrs,
                    SEXP s_dist_normalization_factor,
                    SEXP s_range_adjust,
                    SEXP s_range_factor);

    SEXP S_nn_kmean_cv(SEXP s_graph,
                       SEXP s_edge_lengths,
                       SEXP s_y,
                       SEXP s_range_adjust,
                       SEXP s_range_factor,
                       SEXP s_imputation_method,
                       SEXP s_max_iterations,
                       SEXP s_convergence_threshold,
                       SEXP s_apply_binary_threshold,
                       SEXP s_binary_threshold,
                       SEXP s_ikernel,
                       SEXP s_n_itrs,
                       SEXP s_dist_normalization_factor,
                       SEXP s_n_CVs,
                       SEXP s_n_CV_folds,
                       SEXP s_epsilon,
                       SEXP s_use_low_pass_filter,
                       SEXP s_preserve_local_extrema,
                       SEXP s_min_plambda,
                       SEXP s_max_plambda,
                       SEXP s_seed);

    SEXP S_vanilla_nn_kmean_cv(SEXP s_graph,
                               SEXP s_edge_lengths,
                               SEXP s_y,
                               SEXP s_imputation_method,
                               SEXP s_max_iterations,
                               SEXP s_convergence_threshold,
                               SEXP s_apply_binary_threshold,
                               SEXP s_binary_threshold,
                               SEXP s_ikernel,
                               SEXP s_n_itrs,
                               SEXP s_dist_normalization_factor,
                               SEXP s_n_CVs,
                               SEXP s_n_CV_folds,
                               SEXP s_epsilon,
                               SEXP s_seed);
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
 * @param s_graph A list of integer vectors representing the adjacency matrix of the graph. Each inner vector contains the indices of neighboring vertices.
 * @param s_y A numeric vector of values associated with each vertex in the graph.
 * @param s_n_itrs The number of NN mean iterations.
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
SEXP S_nn_mean(SEXP s_graph, SEXP s_y, SEXP s_n_itrs) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);

    int y_length = LENGTH(s_y);
    std::vector<double> y_vec(REAL(s_y), REAL(s_y) + y_length);

    int n_itrs = INTEGER(s_n_itrs)[0];

    // Call the nn_mean function
    std::unique_ptr<std::vector<double>> result = nn_mean(graph, y_vec, n_itrs);

    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(s_result)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return s_result;
}


/**
 * Computes the weighted mean of a function over vertices of a graph using
 * normalized edge weights.
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
 * @param s_graph An R list representing the graph structure.
 * @param s_w An R list of numeric vectors representing edge weights.
 * @param s_y An R numeric vector of values associated with each vertex.
 * @param s_n_itrs An R integer specifying the number of iterations.
 * @param s_epsilon An R numeric value for the epsilon parameter.
 *
 * @return An R numeric vector containing the weighted nearest neighbor mean for each vertex.
 */
SEXP S_nn_wmean(SEXP s_graph, SEXP s_w, SEXP s_y, SEXP s_n_itrs, SEXP s_epsilon) {
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> weights = Rweights_to_vector(s_w);
    int y_length = LENGTH(s_y);
    std::vector<double> y_vec(REAL(s_y), REAL(s_y) + y_length);
    int n_itrs = INTEGER(s_n_itrs)[0];
    double epsilon = REAL(s_epsilon)[0];

    // Call the nn_wmean function
    std::unique_ptr<std::vector<double>> result = nn_wmean(graph, weights, y_vec, n_itrs, epsilon);

    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(s_result)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return s_result;
}


/**
 * @brief Computes the kernel-weighted nearest neighbor mean of a function over vertices of a graph.
 *
 * This function calculates the kernel-weighted mean of a function y over the vertices of a graph G.
 * For each vertex v in G, the kernel-weighted mean of y at v, kmean(y|v, G), is computed as follows:
 *
 * 1. Consider v and all its neighboring vertices (including v itself).
 * 2. Assign weights to these vertices based on their edge lengths to v, using a specified kernel function.
 * 3. Compute the weighted average of y values over these vertices.
 *
 * The edge length for v to itself is set to 0, ensuring it always receives the maximum kernel weight.
 * This process can be iterated multiple times, using the results from the previous iteration as input for the next.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector
 *              contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors of edge lengths. edge_lengths[i][j] is the length
 *                     of the edge between the i-th vertex and its j-th neighbor as listed in graph[i].
 * @param y A vector of double values associated with each vertex in the graph.
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 * @param n_itrs The number of iterations for computing the nearest neighbor mean.
 * @param dist_normalization_factor A scaling factor applied to the maximum distance between
 *                                  a vertex and its neighbors. This ensures non-zero weights
 *                                  even when all distances are equal. Default value is 1.01.
 * @param range_adjust A boolean flag indicating whether to adjust the range of the results
 *                     after each iteration (except the first). Default is false.
 * @param range_factor A factor used to expand the range of the results when range_adjust is true.
 *                     The range is expanded by range_factor * (max - min) in each direction.
 *                     Default value is 0.01.
 *
 * @return A unique pointer to a vector of doubles containing the kernel-weighted
 *         nearest neighbor mean for each vertex.
 *
 * @note The function's memory management is handled automatically using std::unique_ptr.
 * @note When range_adjust is true, the function uses the min and max values of the first
 *       kernel mean (itr = 0) to adjust the range of subsequent iterations.
 *
 * @example
 *   std::vector<std::vector<int>> graph = {{1, 2}, {0, 2}, {0, 1}};
 *   std::vector<std::vector<double>> edge_lengths = {{0.1, 0.2}, {0.1, 0.3}, {0.2, 0.3}};
 *   std::vector<double> y = {1.0, 2.0, 3.0};
 *   int ikernel = 1;  // Epanechnikov kernel
 *   int n_itrs = 3;
 *   bool range_adjust = true;
 *   double range_factor = 0.02;
 *   auto result = nn_kmean(graph, edge_lengths, y, ikernel, n_itrs, 1.01, range_adjust, range_factor);
 */
std::unique_ptr<std::vector<double>> nn_kmean(const std::vector<std::vector<int>>& graph,
                                              const std::vector<std::vector<double>>& edge_lengths,
                                              const std::vector<double>& y,
                                              int ikernel,
                                              int n_itrs,
                                              double dist_normalization_factor = 1.01,
                                              bool range_adjust = false,
                                              double range_factor = 0.01) {

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

        for (size_t i = 0; i < graph.size(); ++i) {
            if (!graph[i].empty()) {
                distances[0] = 0;  // Distance to self is 0
                double max_dist = 0.0;
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    distances[j + 1] = edge_lengths[i][j];
                    if (distances[j + 1] > max_dist)
                        max_dist = distances[j + 1];
                }

                if (max_dist == 0) max_dist = 1;  // Avoid division by zero

                max_dist *= dist_normalization_factor;

                for (size_t j = 0; j < graph[i].size(); ++j)
                    distances[j + 1] /= max_dist;

                int n = graph[i].size() + 1;
                kernel_fn(distances.data(), n, kernel_weights.data());

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
            } else {
                (*temp)[i] = (itr == 0) ? y[i] : (*result)[i];
            } // END OF if (!graph[i].empty())
        } // END OF for (size_t i = 0; i < graph.size(); ++i)

        if (itr > 0 && range_adjust) {
            double ymin = *std::min_element((*result).begin(), (*result).end());
            double ymax = *std::max_element((*result).begin(), (*result).end());
            if (range_factor != 0) {
                // expanding the range interval by range_factor * delta in each direction
                double delta = range_factor * (ymax - ymin);
                ymin -= delta;
                ymax += delta;
            }
            scale_to_range(*temp, ymin, ymax);
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
 * @param s_graph An R list representing the graph structure.
 * @param s_edge_lengths An R list of numeric vectors representing edge lengths.
 * @param s_y An R numeric vector of values associated with each vertex.
 * @param s_ikernel An R integer specifying the kernel function to use.
 * @param s_n_itrs An R integer specifying the number of iterations.
 * @param s_dist_normalization_factor An R numeric value. A scaling factor applied to the maximum
 *        distance between a vertex and its neighbors. This ensures non-zero weights
 *        even when all distances are equal, by slightly increasing the normalization range.
 * @param s_range_adjust An R logical value indicating whether to adjust the range of the results
 *        after each iteration (except the first).
 * @param s_range_factor An R numeric value used to expand the range of the results when range_adjust is true.
 *        The range is expanded by range_factor * (max - min) in each direction.
 *
 * @return An R numeric vector containing the kernel-weighted nearest neighbor mean for each vertex.
 *
 * @note Default values for s_range_adjust and s_range_factor should be set in R, not in this C++ interface.
 */
SEXP S_nn_kmean(SEXP s_graph,
                SEXP s_edge_lengths,
                SEXP s_y,
                SEXP s_ikernel,
                SEXP s_n_itrs,
                SEXP s_dist_normalization_factor,
                SEXP s_range_adjust,
                SEXP s_range_factor) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);

    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);

    int ikernel = INTEGER(s_ikernel)[0];
    int n_itrs = INTEGER(s_n_itrs)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool range_adjust = LOGICAL(s_range_adjust)[0];
    double range_factor = REAL(s_range_factor)[0];

    // Call the nn_kmean function
    std::unique_ptr<std::vector<double>> result = nn_kmean(graph,
                                                           edge_lengths,
                                                           y,
                                                           ikernel,
                                                           n_itrs,
                                                           dist_normalization_factor,
                                                           range_adjust,
                                                           range_factor);

    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(s_result)[i] = (*result)[i];
    }
    UNPROTECT(1);

    return s_result;
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
 * @param range_adjust A boolean flag indicating whether to adjust the range of the results
 *                     after each iteration (except the first). Default is false.
 * @param range_factor A factor used to expand the range of the results when range_adjust is true.
 *                     The range is expanded by range_factor * (max - min) in each direction.
 *                     Default value is 0.01.
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
                                                 bool range_adjust = false,
                                                 double range_factor = 0.01,
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
            std::vector<double> Ecv_y = std::move(*nn_kmean(graph, edge_lengths, cv_y, ikernel, n_itrs, dist_normalization_factor, range_adjust, range_factor));

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

            #if 0
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
            #endif

            // For non-binary outcomes, compute mean absolute error
            for (const auto& vertex : test_set) {
                (*cv_error)[vertex] += std::abs(Ecv_y[vertex] - y[vertex]);
                cv_error_count[vertex]++;
            }
            //}

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
 * @param s_graph An R list representing the graph structure. Each element of the list
 *               should be an integer vector containing the indices of neighboring vertices.
 * @param s_edge_lengths An R list of numeric vectors representing edge lengths.
 *           The structure should match that of s_graph. If empty, the hop index strategy
 *           will be used for imputation. Otherwise, the kernel distance strategy will be used.
 * @param s_y An R numeric vector of values associated with each vertex in the graph.
 * @param s_range_adjust An R logical value indicating whether to adjust the range of the results
 *        after each iteration (except the first).
 * @param s_range_factor An R numeric value used to expand the range of the results when range_adjust is true.
 *        The range is expanded by range_factor * (max - min) in each direction.
 * @param s_imputation_method An R integer specifying the imputation method to use.
 *                           Valid values correspond to the imputation_method_t enum.
 * @param s_max_iterations The number of iterations in the iterative matching method.
 * @param s_convergence_threshold The convergence threshold in the iterative matching method.
 * @param s_apply_binary_threshold An R logical indicating whether to apply binary thresholding.
 * @param s_binary_threshold An R numeric value specifying the threshold for binary classification.
 * @param s_ikernel An R integer specifying the kernel function to use. Valid values are:
 *                 1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 * @param s_n_itrs An R integer specifying the number of iterations for computing
 *                the nearest neighbor mean.
 * @param s_dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization range.
 * @param s_n_CVs An R integer specifying the number of cross-validation rounds.
 * @param s_n_CV_folds An R integer specifying the number of folds in each cross-validation round.
 * @param s_epsilon An R numeric value specifying a small positive constant for
 *                 clipping estimated conditional expectation values of a binary variable.
 * @param s_use_low_pass_filter An R logical indicating whether to apply low-pass filtering.
 * @param s_preserve_local_extrema An R logical indicating whether to preserve local extrema
 *                                during low-pass filtering.
 * @param s_min_plambda An R numeric value specifying the minimum proportion of eigenvectors
 *                     to use in spectral smoothing.
 * @param s_max_plambda An R numeric value specifying the maximum proportion of eigenvectors
 *                     to use in spectral smoothing.
 * @param s_seed An R integer used as a seed for the random number generator.
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
SEXP S_nn_kmean_cv(SEXP s_graph,
                   SEXP s_edge_lengths,
                   SEXP s_y,
                   SEXP s_range_adjust,
                   SEXP s_range_factor,
                   SEXP s_imputation_method,
                   SEXP s_max_iterations,
                   SEXP s_convergence_threshold,
                   SEXP s_apply_binary_threshold,
                   SEXP s_binary_threshold,
                   SEXP s_ikernel,
                   SEXP s_n_itrs,
                   SEXP s_dist_normalization_factor,
                   SEXP s_n_CVs,
                   SEXP s_n_CV_folds,
                   SEXP s_epsilon,
                   SEXP s_use_low_pass_filter,
                   SEXP s_preserve_local_extrema,
                   SEXP s_min_plambda,
                   SEXP s_max_plambda,
                   SEXP s_seed) {
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y_vec(REAL(s_y), REAL(s_y) + y_length);
    bool range_adjust = LOGICAL(s_range_adjust)[0];
    double range_factor = REAL(s_range_factor)[0];
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(s_imputation_method)[0]);
    int max_iterations = INTEGER(s_max_iterations)[0];
    double convergence_threshold = REAL(s_convergence_threshold)[0];
    bool apply_binary_threshold = LOGICAL(s_apply_binary_threshold)[0];
    double binary_threshold = REAL(s_binary_threshold)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_itrs = INTEGER(s_n_itrs)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool use_low_pass_filter = LOGICAL(s_use_low_pass_filter)[0];
    bool preserve_local_extrema = LOGICAL(s_preserve_local_extrema)[0];
    double min_plambda = REAL(s_min_plambda)[0];
    double max_plambda = REAL(s_max_plambda)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Call the nn_kmean_cv function
    std::unique_ptr<std::vector<double>> result = nn_kmean_cv(graph,
                                                              edge_lengths,
                                                              y_vec,
                                                              range_adjust,
                                                              range_factor,
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


// -------------------------------------------------------------------------------------
//
// vanilla nn_kmean and nn_kmean_cv
//
// -------------------------------------------------------------------------------------

std::unique_ptr<std::vector<double>> vanilla_nn_kmean(const std::vector<std::vector<int>>& graph,
                                                      const std::vector<std::vector<double>>& edge_lengths,
                                                      const std::vector<double>& y,
                                                      int ikernel,
                                                      int n_itrs,
                                                      double dist_normalization_factor = 1.01) {


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

        for (size_t i = 0; i < graph.size(); ++i) {
            if (!graph[i].empty()) {
                distances[0] = 0;  // Distance to self is 0
                double max_dist = 0.0;
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    distances[j + 1] = edge_lengths[i][j];
                    if (distances[j + 1] > max_dist)
                        max_dist = distances[j + 1];
                }

                if (max_dist == 0) max_dist = 1;  // Avoid division by zero

                max_dist *= dist_normalization_factor;

                for (size_t j = 0; j < graph[i].size(); ++j)
                    distances[j + 1] /= max_dist;

                int n = graph[i].size() + 1;
                kernel_fn(distances.data(), n, kernel_weights.data());

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
            } else {
                (*temp)[i] = (itr == 0) ? y[i] : (*result)[i];
            } // END OF if (!graph[i].empty())
        } // END OF for (size_t i = 0; i < graph.size(); ++i)

        result = std::move(temp);
    }

    return result;
}



std::unique_ptr<std::vector<double>> vanilla_nn_kmean_cv(const std::vector<std::vector<int>>& graph,
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
                                                         unsigned int seed = 0) {
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
            std::vector<double> Ecv_y = std::move(*vanilla_nn_kmean(graph, edge_lengths, cv_y, ikernel, n_itrs, dist_normalization_factor));

            // Computing the mean absolute deviation error
            for (const auto& vertex : test_set) {
                (*cv_error)[vertex] += std::abs(Ecv_y[vertex] - y[vertex]);
                cv_error_count[vertex]++;
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex])
            (*cv_error)[vertex] /= cv_error_count[vertex];
    }

    return cv_error;
}


SEXP S_vanilla_nn_kmean_cv(SEXP s_graph,
                           SEXP s_edge_lengths,
                           SEXP s_y,
                           SEXP s_imputation_method,
                           SEXP s_max_iterations,
                           SEXP s_convergence_threshold,
                           SEXP s_apply_binary_threshold,
                           SEXP s_binary_threshold,
                           SEXP s_ikernel,
                           SEXP s_n_itrs,
                           SEXP s_dist_normalization_factor,
                           SEXP s_n_CVs,
                           SEXP s_n_CV_folds,
                           SEXP s_epsilon,
                           SEXP s_seed) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y_vec(REAL(s_y), REAL(s_y) + y_length);
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(s_imputation_method)[0]);
    int max_iterations = INTEGER(s_max_iterations)[0];
    double convergence_threshold = REAL(s_convergence_threshold)[0];
    bool apply_binary_threshold = LOGICAL(s_apply_binary_threshold)[0];
    double binary_threshold = REAL(s_binary_threshold)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_itrs = INTEGER(s_n_itrs)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Call the nn_kmean_cv function
    std::unique_ptr<std::vector<double>> result = vanilla_nn_kmean_cv(graph,
                                                                      edge_lengths,
                                                                      y_vec,
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
                                                                      seed);
    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return Rresult;
}




// ----------------------------------------------------------------------------------------------------
//
// 2-skeleton based kmean estimates
//
// ----------------------------------------------------------------------------------------------------
#if 0

Here is an implementation of an algorithm that computes the kernel-weighted
nearest neighbor mean of a function over vertices of a graph. The graph is
constructed from a dataset using kNN method. More precisely, it is the
1-skeleton of the nerve of the covering of the dataset by kNN sets. I wonder if
adding 2-skeleton would allow for more refined estimate of the conditional
expectation of y over the 2-skeleton nerve of the kNN covering? To be precise a
three points p1, p2, p3 of the dataset X would be part of a 2-dim simplex
\Delta, if the corresponding kNN sets N(p1), N(p2), N(p3) all intersect.
Supplemting the existence of the simplex could be a data structure that holds
the area of the triangle span by p1, p2 and p3. This would be a 2-dim analogue
of the length of an edge in the graph case. Do you have any ideas on how this
extra structure be used to build a more refined model of the conditional
    expectation of y?


Thank you for sharing this interesting algorithm and your thoughts on extending it. Your idea of incorporating the 2-skeleton of the nerve complex is intriguing and could potentially lead to a more refined estimate of the conditional expectation. Let's explore how we might extend this concept and use the additional structure.

This extended version of the algorithm incorporates the 2-skeleton (triangles) into the estimation process. Here's a breakdown of the key changes and how they contribute to a more refined model:

1. Triangle Structure: We introduce a `Triangle` struct that holds the indices of the three vertices and the area of the triangle. This captures the 2-dimensional simplex information.

2. Vertex-Triangle Mapping: We create a `vertex_triangles` vector that maps each vertex to the triangles it's part of. This allows us to efficiently access the relevant 2-skeleton information for each vertex.

3. Triangle Contributions: In the main computation loop, we now include contributions from the triangles in addition to the edges. The weight of a triangle's contribution is based on its area, normalized by the square of the maximum distance (since area scales quadratically with distance).

4. Weighted Average: The final estimate for each vertex now includes:
   - The vertex's own value
   - Contributions from connected edges (1-skeleton)
   - Contributions from triangles it's part of (2-skeleton)

This approach allows for a more nuanced estimation of the conditional expectation by incorporating higher-dimensional relationships between data points. The triangle contributions can capture non-linear relationships that might not be apparent from edge connections alone.

To further refine this model, you could consider:

1. Adaptive Weighting: Adjust the relative importance of edge vs. triangle contributions based on local data characteristics.

2. Higher-order Simplices: Extend to 3-simplices (tetrahedra) or beyond, though the computational complexity would increase significantly.

3. Persistent Homology: Incorporate ideas from topological data analysis to focus on the most significant simplicial structures.

4. Non-Euclidean Geometries: Adapt the area and distance calculations for non-Euclidean spaces if your data lies on a manifold.

5. Machine Learning Integration: Use the simplex structure as features in a machine learning model to predict y values.



[PG] Thank you! it is very insightful. Don't you think that triangle_weight should also be computed using kernel function the same way edge weights are computed? In this way triangles with small area would contribute more to the final estimate.


You've made an excellent observation! You're absolutely right that applying a kernel function to the triangle weights, similar to how we handle edge weights, could provide a more nuanced and potentially more accurate contribution from the 2-skeleton structure. This approach would indeed give more weight to smaller triangles, which could be beneficial in capturing local structure more precisely.

Let's modify the algorithm to incorporate this insight. We'll use the kernel function for both edge and triangle weights.

Here are the key changes and their implications:

1. Triangle Distance Metric: We now use the square root of the triangle area as a "distance" metric for triangles. This choice allows us to map the 2D area to a 1D distance-like measure, which can be input into our kernel function.

2. Normalization: Both edge distances and triangle "distances" are normalized separately. This ensures that both types of structures contribute meaningfully to the final estimate, regardless of the scale differences between edge lengths and triangle areas.

3. Kernel Weights for Triangles: We now apply the kernel function to both edge distances and triangle "distances". This means that smaller triangles (those with smaller areas) will indeed contribute more to the final estimate, just as closer neighbors do in the edge-based calculations.

4. Balanced Contribution: By using the kernel function for both edges and triangles, we ensure that the contribution of each structure type (1-skeleton and 2-skeleton) is balanced and consistent with the overall approach of the algorithm.

This modification allows the algorithm to be more sensitive to local structure in both the 1-skeleton (graph) and 2-skeleton (triangles) representations of the data. Smaller triangles, which might represent areas of the data space where points are clustered more closely, will now have a stronger influence on the estimation.

Some potential benefits of this approach:

1. Improved Local Sensitivity: The algorithm can now better capture fine-grained local structure in the data.
2. Consistency: Using the same kernel function for both edges and triangles provides a consistent framework for weighting different geometric structures.
3. Adaptability: The relative influence of edges vs. triangles can be tuned by adjusting the normalization factors or potentially using different kernel functions for each.

To further refine this approach, you might consider:

1. Adaptive Kernel Selection: Use different kernel functions for edges and triangles based on the characteristics of your data.
2. Higher-order Simplices: Extend this approach to higher-order simplices (tetrahedra, etc.) using appropriate distance metrics.
3. Persistent Homology Integration: Incorporate ideas from topological data analysis to focus on the most significant simplicial structures across different scales.

                                                                                                                              Would you like to discuss any of these potential refinements or explore how this modified algorithm might perform on different types of datasets?
#endif


#if 0
struct Triangle {
    int v1, v2, v3;
    double area;
};

std::unique_ptr<std::vector<double>> extended_nn_kmean(
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<Triangle>& triangles,
    const std::vector<double>& y,
    int ikernel,
    int n_itrs,
    double dist_normalization_factor = 1.01,
    bool range_adjust = false,
    double range_factor = 0.01) {

    auto result = std::make_unique<std::vector<double>>(y.size(), 0.0);
    initialize_kernel(ikernel);

    // Determine the maximum number of neighbors and triangles for any vertex
    size_t max_neighbors = 0;
    size_t max_triangles = 0;
    std::vector<std::vector<const Triangle*>> vertex_triangles(graph.size());
    for (size_t i = 0; i < graph.size(); ++i) {
        max_neighbors = std::max(max_neighbors, graph[i].size());
        for (const auto& triangle : triangles) {
            if (triangle.v1 == i || triangle.v2 == i || triangle.v3 == i) {
                vertex_triangles[i].push_back(&triangle);
            }
        }
        max_triangles = std::max(max_triangles, vertex_triangles[i].size());
    }

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights_edges(max_neighbors + 1);
    std::vector<double> kernel_weights_triangles(max_triangles);
    std::vector<double> distances_edges(max_neighbors + 1);
    std::vector<double> distances_triangles(max_triangles);

    for (int itr = 0; itr < n_itrs; itr++) {
        auto temp = std::make_unique<std::vector<double>>(y.size(), 0.0);

        for (size_t i = 0; i < graph.size(); ++i) {
            if (!graph[i].empty()) {
                distances_edges[0] = 0;  // Distance to self is 0
                double max_dist = 0.0;
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    distances_edges[j + 1] = edge_lengths[i][j];
                    if (distances_edges[j + 1] > max_dist)
                        max_dist = distances_edges[j + 1];
                }
                if (max_dist == 0) max_dist = 1;  // Avoid division by zero
                max_dist *= dist_normalization_factor;

                // Normalize edge distances
                for (size_t j = 0; j < graph[i].size(); ++j)
                    distances_edges[j + 1] /= max_dist;

                // Compute triangle distances (using area as a proxy for distance)
                double max_area = 0.0;
                for (size_t j = 0; j < vertex_triangles[i].size(); ++j) {
                    distances_triangles[j] = std::sqrt(vertex_triangles[i][j]->area);  // Use square root of area as "distance"
                    if (distances_triangles[j] > max_area)
                        max_area = distances_triangles[j];
                }
                if (max_area == 0) max_area = 1;  // Avoid division by zero

                // Normalize triangle distances
                for (size_t j = 0; j < vertex_triangles[i].size(); ++j)
                    distances_triangles[j] /= max_area;

                // Compute kernel weights
                int n_edges = graph[i].size() + 1;
                int n_triangles = vertex_triangles[i].size();
                kernel_fn(distances_edges.data(), n_edges, kernel_weights_edges.data());
                kernel_fn(distances_triangles.data(), n_triangles, kernel_weights_triangles.data());

                double weighted_sum = (itr == 0) ? y[i] : (*result)[i];
                weighted_sum *= kernel_weights_edges[0];
                double weight_sum = kernel_weights_edges[0];

                // Edge contributions
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    int neighbor = graph[i][j];
                    double edge_weight = kernel_weights_edges[j + 1];
                    weighted_sum += edge_weight * ((itr == 0) ? y[neighbor] : (*result)[neighbor]);
                    weight_sum += edge_weight;
                }

                // Triangle contributions
                for (size_t j = 0; j < vertex_triangles[i].size(); ++j) {
                    const Triangle* triangle = vertex_triangles[i][j];
                    int v1 = triangle->v1, v2 = triangle->v2, v3 = triangle->v3;
                    double triangle_weight = kernel_weights_triangles[j];
                    double triangle_value = ((itr == 0) ? (y[v1] + y[v2] + y[v3]) / 3 :
                                             ((*result)[v1] + (*result)[v2] + (*result)[v3]) / 3);
                    weighted_sum += triangle_weight * triangle_value;
                    weight_sum += triangle_weight;
                }

                (*temp)[i] = weighted_sum / weight_sum;
            }
        }

        result = std::move(temp);
    }

    if (range_adjust) {
        // Range adjustment logic (unchanged)
    }

    return result;
}
#endif
