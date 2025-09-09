
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
#include "graph_diffussion_smoother.h"
#include "stats_utils.h"
#include "kernels.h"

extern "C" {
    SEXP S_graph_kernel_weighted_mean(SEXP s_graph,
                                      SEXP s_edge_lengths,
                                      SEXP s_y,
                                      SEXP s_ikernel,
                                      SEXP s_dist_normalization_factor);

    SEXP S_graph_kernel_weighted_mean_cv(SEXP s_graph,
                                     SEXP s_edge_lengths,
                                     SEXP s_y,
                                     SEXP s_ikernel,
                                     SEXP s_dist_normalization_factor,
                                     SEXP s_n_CVs,
                                     SEXP s_n_CV_folds,
                                     SEXP s_epsilon,
                                     SEXP s_seed);
}


/**
 * @brief Computes the kernel-weighted nearest neighbor mean of a function over vertices of a graph, with support for vertex weighting.
 *
 * This function calculates the kernel-weighted mean of values associated with vertices in a graph,
 * taking into account user-provided weights for each vertex. It supports excluding vertices with
 * zero total weight and handles cases where vertices have no neighbors.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains
 *              the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths (or distances) of edges in the graph.
 *                     It should have the same structure as the 'graph' parameter.
 * @param weights A vector of weights for each vertex in the graph. A weight of 0 will cause
 *                the vertex to be excluded from calculations.
 * @param y A vector of values associated with each vertex in the graph.
 * @param ikernel An integer specifying the kernel function to use for weighting.
 * @param dist_normalization_factor A factor used to normalize distances in the graph.
 *                                  Default value is 1.01.
 *
 * @return A std::pair containing:
 *         - First: A vector of doubles representing the kernel-weighted mean for each vertex.
 *         - Second: A vector of integers containing the indices of excluded vertices (those with zero total weight).
 *
 * @note Vertices with no neighbors (empty adjacency list in 'graph') are handled specially:
 *       - If the vertex's weight is 0, it's added to the excluded vertices list.
 *       - If the vertex's weight is non-zero, its original value from 'y' is used.
 *
 * @warning This function assumes that the input vectors (graph, edge_lengths, weights, y)
 *          all have the same number of elements, corresponding to the total number of vertices in the graph.
 *
 * @pre The function 'initialize_kernel' must be called before this function to set up the kernel function.
 * @pre The function 'kernel_fn' must be properly defined to calculate kernel weights based on distances.
 *
 * @throws May throw exceptions if vector accesses are out of bounds or if mathematical operations
 *         (like division by zero) fail. Proper input validation should be done before calling this function.
 *
 * @see initialize_kernel, kernel_fn
 */
std::pair<std::vector<double>, std::vector<int>> graph_kernel_weighted_mean_with_weights(const std::vector<std::vector<int>>& graph,
                                                                                         const std::vector<std::vector<double>>& edge_lengths,
                                                                                         const std::vector<double>& weights,
                                                                                         const std::vector<double>& y,
                                                                                         int ikernel,
                                                                                         double dist_normalization_factor = 1.01) {
    auto kmean = std::vector<double>(y.size(), 0.0);
    std::vector<int> excluded_vertices;
    initialize_kernel(ikernel);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

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

            double weight_sum = weights[i] * kernel_weights[0];
            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                weight_sum += weights[neighbor] * kernel_weights[j + 1];
            }

            if (weight_sum == 0) {
                excluded_vertices.push_back(i);
                kmean[i] = 0; // or any other appropriate value to indicate exclusion
            } else {
                double weighted_sum = y[i] * weights[i] * kernel_weights[0];
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    int neighbor = graph[i][j];
                    weighted_sum += weights[neighbor] * kernel_weights[j + 1] * y[neighbor];
                }
                kmean[i] = weighted_sum / weight_sum;
            }
        } else {
            if (weights[i] == 0) {
                excluded_vertices.push_back(i);
                kmean[i] = 0; // or any other appropriate value to indicate exclusion
            } else {
                kmean[i] = y[i];
            }
        }
    }

    return std::make_pair(kmean, excluded_vertices);
}

/**
* Computes Kernel-Weighted Mean on a Graph
*
* This function calculates the kernel-weighted mean of values associated with vertices in a graph.
* It uses the graph structure, edge lengths, and a specified kernel function to compute weighted
* averages for each vertex.
*
* @param graph A list of integer vectors. Each element represents a vertex, and the vector
*   contains indices of its neighboring vertices.
* @param edge_lengths A list of numeric vectors. Each element corresponds to a vertex, and the
*   vector contains the lengths (or distances) of edges to its neighbors. The structure should
*   match that of the `graph` parameter.
* @param y A numeric vector of values associated with each vertex in the graph.
* @param ikernel An integer specifying the kernel function to use for weighting.
* @param dist_normalization_factor A numeric value used to normalize distances in the graph.
*   Default is 1.01.
*
* @return A numeric vector containing the kernel-weighted mean for each vertex in the graph.
*
* @details
* The function performs the following steps for each vertex:
* 1. Normalizes the distances to its neighbors.
* 2. Applies the specified kernel function to these normalized distances.
* 3. Computes a weighted average of the vertex's value and its neighbors' values,
*    using the kernel weights.
*
* For vertices with no neighbors, their original value from `y` is retained.
*
*/
std::vector<double> graph_kernel_weighted_mean(const std::vector<std::vector<int>>& graph,
                                               const std::vector<std::vector<double>>& edge_lengths,
                                               const std::vector<double>& y,
                                               int ikernel,
                                               double dist_normalization_factor = 1.01) {
    auto kmean = std::vector<double>(y.size(), 0.0);
    initialize_kernel(ikernel);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

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

            double weighted_sum = y[i] * kernel_weights[0];
            double weight_sum = kernel_weights[0];

            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                weighted_sum += kernel_weights[j + 1] * y[neighbor];
                weight_sum += kernel_weights[j + 1];
            }

            kmean[i] = weighted_sum / weight_sum;
        } else {
            kmean[i] = y[i];
        }
    }
    return kmean;
}

/**
 * @brief R interface for computing kernel-weighted mean on a graph.
 *
 * This function serves as an interface between R and the C++ implementation of
 * the graph kernel-weighted mean algorithm. It takes R objects as input,
 * converts them to C++ data structures, calls the C++ function, and returns
 * the result as an R object.
 *
 * @param s_graph An R list of integer vectors representing the graph structure.
 *        Each element of the list corresponds to a vertex, and contains the
 *        indices of its neighboring vertices.
 * @param s_edge_lengths An R list of numeric vectors containing the lengths
 *        (or distances) of edges in the graph. The structure should match
 *        that of s_graph.
 * @param s_y An R numeric vector of values associated with each vertex in
 *        the graph.
 * @param s_ikernel An R integer specifying the kernel function to use for
 *        weighting.
 * @param s_dist_normalization_factor An R numeric value used to normalize
 *        distances in the graph.
 *
 * @return An R numeric vector containing the kernel-weighted mean for each
 *         vertex in the graph.
 *
 * @details
 * The function performs the following steps:
 * 1. Converts R data structures to C++ data structures.
 * 2. Calls the C++ implementation of graph_kernel_weighted_mean.
 * 3. Converts the C++ result back to an R numeric vector.
 *
 * @note
 * - This function assumes that the input R objects are of the correct type
 *   and structure. No extensive error checking is performed.
 * - The function uses PROTECT/UNPROTECT for memory management as per R's
 *   C interface guidelines.
 *
 * @see graph_kernel_weighted_mean, R_list_of_ivectors_to_cpp_vector_of_ivectors,
 *      R_list_of_dvectors_to_cpp_vector_of_dvectors
 */
SEXP S_graph_kernel_weighted_mean(SEXP s_graph,
                                  SEXP s_edge_lengths,
                                  SEXP s_y,
                                  SEXP s_ikernel,
                                  SEXP s_dist_normalization_factor) {

    std::vector<std::vector<int>> graph = std::move(*R_list_of_ivectors_to_cpp_vector_of_ivectors(s_graph));
    std::vector<std::vector<double>> edge_lengths = std::move(*R_list_of_dvectors_to_cpp_vector_of_dvectors(s_edge_lengths));
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    std::vector<double> result = graph_kernel_weighted_mean(graph,
                                                            edge_lengths,
                                                            y,
                                                            ikernel,
                                                            dist_normalization_factor);
    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result.size()));
    for (size_t i = 0; i < result.size(); ++i) {
        REAL(s_result)[i] = result[i];
    }
    UNPROTECT(1);

    return s_result;
}


/**
 * @brief Performs cross-validation for kernel-weighted mean estimation on a graph.
 *
 * This function implements a cross-validation procedure to evaluate the performance
 * of kernel-weighted mean estimation on a graph structure. It supports both binary
 * and continuous outcomes, using appropriate error metrics for each case.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector
 *        contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths (or distances) of edges
 *        in the graph. It should have the same structure as the 'graph' parameter.
 * @param y A vector of values associated with each vertex in the graph. These can be
 *        binary (0 or 1) or continuous values.
 * @param ikernel An integer specifying the kernel function to use for weighting.
 *        Default is 1.
 * @param dist_normalization_factor A factor used to normalize distances in the graph.
 *        Default is 1.01.
 * @param n_CVs The number of cross-validation iterations to perform. Must be greater than 0.
 * @param n_CV_folds The number of folds to use in each cross-validation iteration.
 *        Default is 10.
 * @param epsilon A small value used to avoid numerical issues in log calculations for
 *        binary outcomes. Default is 1e-10.
 * @param seed Seed for the random number generator. If 0, the current time is used.
 *        Default is 0.
 *
 * @return A vector of doubles containing the average cross-validation error for each vertex.
 *         Vertices for which no error could be computed (due to being excluded in all
 *         iterations) will have a NaN value.
 *
 * @throws std::invalid_argument if n_CVs is less than or equal to 0.
 *
 * @details
 * The function performs the following steps:
 * 1. Initializes data structures and random number generator.
 * 2. For each CV iteration:
 *    a. Selects a test set of vertices.
 *    b. Computes kernel-weighted means using the remaining vertices as training data.
 *    c. Calculates errors for the test set vertices.
 * 3. Averages the errors for each vertex across all CV iterations.
 *
 * For binary outcomes (y contains only 0 and 1), it computes cross-entropy loss.
 * For continuous outcomes, it computes mean absolute error.
 *
 * @note
 * - The function assumes that the graph_kernel_weighted_mean_with_weights function
 *   is available and correctly implemented.
 * - Vertices that are excluded in all CV iterations (due to graph structure or
 *   weighting) will have NaN as their final error value.
 * - The function uses C++11 random number generation facilities for reproducibility.
 *
 * @see graph_kernel_weighted_mean_with_weights
 */
std::vector<double> graph_kernel_weighted_mean_cv(const std::vector<std::vector<int>>& graph,
                                                  const std::vector<std::vector<double>>& edge_lengths,
                                                  const std::vector<double>& y,
                                                  int ikernel = 1,
                                                  double dist_normalization_factor = 1.01,
                                                  int n_CVs = 0,
                                                  int n_CV_folds = 10,
                                                  double epsilon = 1e-10,
                                                  unsigned int seed = 0) {
    if (n_CVs <= 0) {
        throw std::invalid_argument("n_CVs has to be greater than 0");
    }

    int n_vertices = y.size();

    //bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    bool y_binary = false;

    auto cv_error = std::vector<double>(n_vertices, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> cv_error_count(n_vertices, 0);

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    // Creating a set version of the adjacency matrix of the graph
    std::vector<std::set<int>> set_graph(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
    }

    int fold_size = n_vertices / n_CV_folds;

    std::vector<double> weights(n_vertices, 1.0);

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
                int vertex = uni(rng);
                test_set.insert(vertex);
            }
        }

        // Reseting all weights to 1 and set weights over test_set to 0
        std::fill(weights.begin(), weights.end(), 1.0);
        for (const auto& vertex : test_set) {
            weights[vertex] = 0.0;
        }

        // Estimating the conditional expectation of cv_y
        auto res = graph_kernel_weighted_mean_with_weights(graph,
                                                           edge_lengths,
                                                           weights,
                                                           y,
                                                           ikernel,
                                                           dist_normalization_factor);
        std::vector<double> Ecv_y = res.first;
        std::vector<int> excluded_vertices = res.second;

        // Computing a set difference between test_set and excluded_vertices
        std::set<int> valid_test_set;
        for (const auto& vertex : test_set) {
            if (std::find(excluded_vertices.begin(), excluded_vertices.end(), vertex) == excluded_vertices.end()) {
                valid_test_set.insert(vertex);
            }
        }

        // Checking if valid_test_set is empty
        if (valid_test_set.empty()) {
            continue;  // Skip this iteration if no valid test vertices
        }

        // Computing error over test vertices
        if (y_binary) {
            // For binary outcomes, computing cross-entropy loss
            for (const auto& vertex : valid_test_set) {
                // Clip Ecv_y to avoid log(0) or log(1) issues
                double clipped_Ecv_y = std::max(epsilon, std::min(1.0 - epsilon, Ecv_y[vertex]));
                // Compute negative log-likelihood (cross-entropy loss)
                cv_error[vertex] = (std::isnan(cv_error[vertex])) ?
                    -(y[vertex] * log(clipped_Ecv_y) + (1 - y[vertex]) * log(1 - clipped_Ecv_y)) :
                    cv_error[vertex] + -(y[vertex] * log(clipped_Ecv_y) + (1 - y[vertex]) * log(1 - clipped_Ecv_y));
                cv_error_count[vertex]++;
            }
        } else {
            // For non-binary outcomes, compute mean absolute error
            for (const auto& vertex : valid_test_set) {
                cv_error[vertex] = (std::isnan(cv_error[vertex])) ?
                    std::abs(Ecv_y[vertex] - y[vertex]) :
                    cv_error[vertex] + std::abs(Ecv_y[vertex] - y[vertex]);
                cv_error_count[vertex]++;
            }
        }

    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV error, leaving NaN for vertices with no estimates
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex] > 0) {
            cv_error[vertex] /= cv_error_count[vertex];
        }
    }

    return cv_error;
}

/**
 * @brief R interface for cross-validation of kernel-weighted mean on a graph.
 *
 * This function serves as an interface between R and the C++ implementation of
 * the cross-validation procedure for graph kernel-weighted mean algorithm.
 * It takes R objects as input, converts them to C++ data structures, calls
 * the C++ function, and returns the cross-validation errors as an R object.
 *
 * @param s_graph An R list of integer vectors representing the graph structure.
 *        Each element of the list corresponds to a vertex, and contains the
 *        indices of its neighboring vertices.
 * @param s_edge_lengths An R list of numeric vectors containing the lengths
 *        (or distances) of edges in the graph. The structure should match
 *        that of s_graph.
 * @param s_y An R numeric vector of values associated with each vertex in
 *        the graph.
 * @param s_ikernel An R integer specifying the kernel function to use for
 *        weighting.
 * @param s_dist_normalization_factor An R numeric value used to normalize
 *        distances in the graph.
 * @param s_n_CVs An R integer specifying the number of cross-validation
 *        iterations to perform.
 * @param s_n_CV_folds An R integer specifying the number of folds to use in
 *        each cross-validation iteration.
 * @param s_epsilon An R numeric value used as a small constant to avoid
 *        numerical issues in error calculations.
 * @param s_seed An R integer used as a seed for the random number generator.
 *
 * @return An R numeric vector containing the cross-validation errors for
 *         each vertex in the graph.
 *
 * @details
 * The function performs the following steps:
 * 1. Converts R data structures to C++ data structures.
 * 2. Calls the C++ implementation of graph_kernel_weighted_mean_cv.
 * 3. Converts the C++ result (cross-validation errors) back to an R numeric vector.
 *
 * @note
 * - This function assumes that the input R objects are of the correct type
 *   and structure. No extensive error checking is performed.
 * - The function uses PROTECT/UNPROTECT for memory management as per R's
 *   C interface guidelines.
 * - Vertices for which no cross-validation error could be computed (due to
 *   being excluded in all iterations) will have a NaN value in the result.
 *
 * @see graph_kernel_weighted_mean_cv, R_list_of_ivectors_to_cpp_vector_of_ivectors,
 *      R_list_of_dvectors_to_cpp_vector_of_dvectors
 */
SEXP S_graph_kernel_weighted_mean_cv(SEXP s_graph,
                                     SEXP s_edge_lengths,
                                     SEXP s_y,
                                     SEXP s_ikernel,
                                     SEXP s_dist_normalization_factor,
                                     SEXP s_n_CVs,
                                     SEXP s_n_CV_folds,
                                     SEXP s_epsilon,
                                     SEXP s_seed) {
    std::vector<std::vector<int>> graph = std::move(*R_list_of_ivectors_to_cpp_vector_of_ivectors(s_graph));
    std::vector<std::vector<double>> edge_lengths = std::move(*R_list_of_dvectors_to_cpp_vector_of_dvectors(s_edge_lengths));
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    std::vector<double> cv_errors =  graph_kernel_weighted_mean_cv(graph,
                                                                   edge_lengths,
                                                                   y,
                                                                   ikernel,
                                                                   dist_normalization_factor,
                                                                   n_CVs,
                                                                   n_CV_folds,
                                                                   epsilon,
                                                                   seed);

    // Convert result to SEXP and return
    SEXP s_cv_errors = PROTECT(allocVector(REALSXP, cv_errors.size()));
    for (size_t i = 0; i < cv_errors.size(); ++i) {
        REAL(s_cv_errors)[i] = cv_errors[i];
    }
    UNPROTECT(1);
    return s_cv_errors;
}
