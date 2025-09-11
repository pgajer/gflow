#include "sampling.h" // for C_runif_simplex()
#include "msr2.h"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "stats_utils.h"
#include "kernels.h"
#include "predictive_errors.hpp"
#include "adaptive_nbhd_size.hpp"

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

#include <R.h>
#include <Rinternals.h>

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>
create_chain_graph(const std::vector<double>& x);

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> create_hHN_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h);

extern "C" {
    SEXP S_graph_kmean(SEXP s_graph,
                       SEXP s_edge_lengths,
                       SEXP s_y,
                       SEXP s_ikernel,
                       SEXP s_dist_normalization_factor);

    SEXP S_graph_kmean_wmad_cv(SEXP s_graph,
                               SEXP s_edge_lengths,
                               SEXP s_y,
                               SEXP s_ikernel,
                               SEXP s_dist_normalization_factor,
                               SEXP s_n_CVs,
                               SEXP s_n_CV_folds,
                               SEXP s_seed,
                               SEXP s_use_weighted_MAD_error);

    SEXP S_graph_kmean_cv(SEXP s_graph,
                          SEXP s_edge_lengths,
                          SEXP s_y,
                          SEXP s_ikernel,
                          SEXP s_dist_normalization_factor,
                          SEXP s_n_CVs,
                          SEXP s_n_CV_folds,
                          SEXP s_seed);

    SEXP S_univariate_gkmm(SEXP s_x,
                           SEXP s_y,
                           SEXP s_y_true,
                           SEXP s_use_median,
                           SEXP s_h_min,
                           SEXP s_h_max,
                           SEXP s_n_CVs,
                           SEXP s_n_CV_folds,
                           SEXP s_p,
                           SEXP s_n_bb,
                           SEXP s_ikernel,
                           SEXP s_n_cores,
                           SEXP s_dist_normalization_factor,
                           SEXP s_epsilon,
                           SEXP s_seed);
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
std::vector<double> graph_kmean(const std::vector<std::vector<int>>& graph,
                                const std::vector<std::vector<double>>& edge_lengths,
                                const std::vector<double>& y,
                                int ikernel,
                                double dist_normalization_factor = 1.01) {
    auto kmean = std::vector<double>(y.size(), 0.0);
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

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
 * 2. Calls the C++ implementation of graph_kmean.
 * 3. Converts the C++ result back to an R numeric vector.
 *
 * @note
 * - This function assumes that the input R objects are of the correct type
 *   and structure. No extensive Rf_error checking is performed.
 * - The function uses PROTECT/UNPROTECT for memory management as per R's
 *   C interface guidelines.
 *
 * @see graph_kmean
 *      R_list_of_dvectors_to_cpp_vector_of_dvectors
 */
SEXP S_graph_kmean(SEXP s_graph,
                   SEXP s_edge_lengths,
                   SEXP s_y,
                   SEXP s_ikernel,
                   SEXP s_dist_normalization_factor) {

    std::vector<std::vector<int>> graph           = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    std::vector<double> result = graph_kmean(graph,
                                             edge_lengths,
                                             y,
                                             ikernel,
                                             dist_normalization_factor);
    // Convert result to SEXP and return
    SEXP s_result = PROTECT(Rf_allocVector(REALSXP, result.size()));
    for (size_t i = 0; i < result.size(); ++i) {
        REAL(s_result)[i] = result[i];
    }
    UNPROTECT(1);

    return s_result;
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
 * @Rf_warning This function assumes that the input vectors (graph, edge_lengths, weights, y)
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
std::pair<std::vector<double>, std::vector<int>> graph_kmean_with_weights(const std::vector<std::vector<int>>& graph,
                                                                          const std::vector<std::vector<double>>& edge_lengths,
                                                                          const std::vector<double>& weights,
                                                                          const std::vector<double>& y,
                                                                          int ikernel,
                                                                          double dist_normalization_factor = 1.01) {
    auto kmean = std::vector<double>(y.size(), 0.0);
    std::vector<int> excluded_vertices;
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

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
                kmean[i] = weights[i] * y[i];
            }
        }
    }

    return std::make_pair(kmean, excluded_vertices);
}

/**
 * @brief Performs cross-validation for kernel-weighted mean estimation on a graph.
 *
 * This function implements a cross-validation procedure to evaluate the performance
 * of kernel-weighted mean estimation on a graph structure. It supports both binary
 * and continuous outcomes, using appropriate Rf_error metrics for each case.
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
 * @param use_weighted_MAD_error A boolean flag to control the use of weighted Mean Absolute
 *        Deviation (MAD) Rf_error. Default is false.
 *        When true, the function computes a weighted MAD Rf_error to address class imbalance in of the binary variable y.
 *        The weights are calculated as follows:
 *        - For y_i = 0: weight = 1 / (1 - q)
 *        - For y_i = 1: weight = 1 / q
 *        where q is the proportion of samples where y_i = 1.
 *        This weighting scheme ensures that errors from both classes contribute equally
 *        to the final Rf_error metric, regardless of their relative frequencies in the dataset.
 *        When false, the standard (unweighted) MAD Rf_error is used.
 *
 * @return A vector of doubles containing the average cross-validation Rf_error for each vertex.
 *         Vertices for which no Rf_error could be computed (due to being excluded in all
 *         iterations) will have a NaN value.
 *
 * @throws std::invalid_argument if n_CVs is less than or equal to 0.
 *
 * @details
 * The function performs the following steps:
 * 1. Initializes data structures and random number generator.
 * 2. If use_weighted_MAD_error is true, computes the weights for each class.
 * 3. For each CV iteration:
 *    a. Selects a test set of vertices.
 *    b. Computes kernel-weighted means using the remaining vertices as training data.
 *    c. Calculates errors for the test set vertices, applying weights if use_weighted_MAD_error is true.
 * 4. Averages the errors for each vertex across all CV iterations.
 *
 * @note
 * - The function assumes that the graph_kmean_with_weights function
 *   is available and correctly implemented.
 * - Vertices that are excluded in all CV iterations (due to graph structure or
 *   weighting) will have NaN as their final Rf_error value.
 * - The function uses C++11 random number generation facilities for reproducibility.
 * - The weighted MAD Rf_error option is particularly useful for imbalanced binary variable y,
 *   as it gives equal importance to errors in both classes.
 *
 * @see graph_kmean_with_weights
 */
std::vector<double> graph_kmean_wmad_cv(const std::vector<std::vector<int>>& graph,
                                        const std::vector<std::vector<double>>& edge_lengths,
                                        const std::vector<double>& y,
                                        int ikernel = 1,
                                        double dist_normalization_factor = 1.01,
                                        int n_CVs = 0,
                                        int n_CV_folds = 10,
                                        unsigned int seed = 0,
                                        bool use_weighted_MAD_error = false) {
    if (n_CVs <= 0) {
        Rf_error("n_CVs has to be greater than 0");
    }

    int n_vertices = y.size();

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

    // Variables for weighted MAD Rf_error
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    use_weighted_MAD_error &= y_binary;
    double n1 = 0;
    if (use_weighted_MAD_error) {
        for (const auto& yi : y) {
            if (yi == 1) n1++;
        }
    }
    double q = n1 / n_vertices;
    double alpha0 = use_weighted_MAD_error ? 1.0 / (1.0 - q) : 1.0;
    double alpha1 = use_weighted_MAD_error ? 1.0 / q : 1.0;

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
        auto res = graph_kmean_with_weights(graph,
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

        // Computing cross-validation Rf_error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double Rf_error = std::abs(Ecv_y[vertex] - y[vertex]);
            double weight = (y[vertex] == 1) ? alpha1 : alpha0;

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = weight * Rf_error;
            } else {
                cv_error[vertex] += weight * Rf_error;
            }
            cv_error_count[vertex]++;
        }
    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV Rf_error, leaving NaN for vertices with no estimates
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
 *        numerical issues in Rf_error calculations.
 * @param s_seed An R integer used as a seed for the random number generator.
 * @param s_use_weighted_MAD_error An R logical value indicating whether to use
 *        weighted Mean Absolute Deviation (MAD) Rf_error calculation.
 *
 * @return An R numeric vector containing the cross-validation errors for
 *         each vertex in the graph.
 *
 * @details
 * The function performs the following steps:
 * 1. Converts R data structures to C++ data structures.
 * 2. Calls the C++ implementation of graph_kmean_cv.
 * 3. Converts the C++ result (cross-validation errors) back to an R numeric vector.
 *
 * If s_use_weighted_MAD_error is TRUE, the function computes a weighted MAD Rf_error
 * to address class imbalance in the binary variable y. The weights are:
 * - For y_i = 0: weight = 1 / (1 - q)
 * - For y_i = 1: weight = 1 / q
 * where q is the proportion of samples where y_i = 1.
 *
 * @note
 * - This function assumes that the input R objects are of the correct type
 *   and structure. No extensive Rf_error checking is performed.
 * - The function uses PROTECT/UNPROTECT for memory management as per R's
 *   C interface guidelines.
 * - Vertices for which no cross-validation Rf_error could be computed (due to
 *   being excluded in all iterations) will have a NaN value in the result.
 * - The weighted MAD Rf_error option is particularly useful for imbalanced datasets
 *   in binary classification problems, as it gives equal importance to errors
 *   in both classes.
 *
 * @see graph_kmean_wmad_cv
 *      R_list_of_dvectors_to_cpp_vector_of_dvectors
 */
SEXP S_graph_kmean_wmad_cv(SEXP s_graph,
                           SEXP s_edge_lengths,
                           SEXP s_y,
                           SEXP s_ikernel,
                           SEXP s_dist_normalization_factor,
                           SEXP s_n_CVs,
                           SEXP s_n_CV_folds,
                           SEXP s_seed,
                           SEXP s_use_weighted_MAD_error) {

    std::vector<std::vector<int>> graph           = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];
    bool use_weighted_MAD_error = (LOGICAL(s_use_weighted_MAD_error)[0] == 1);

    std::vector<double> cv_errors =  graph_kmean_wmad_cv(graph,
                                                         edge_lengths,
                                                         y,
                                                         ikernel,
                                                         dist_normalization_factor,
                                                         n_CVs,
                                                         n_CV_folds,
                                                         seed,
                                                         use_weighted_MAD_error);
    // Convert result to SEXP and return
    SEXP s_cv_errors = PROTECT(Rf_allocVector(REALSXP, cv_errors.size()));
    for (size_t i = 0; i < cv_errors.size(); ++i) {
        REAL(s_cv_errors)[i] = cv_errors[i];
    }
    UNPROTECT(1);
    return s_cv_errors;
}



/**
 * @brief Performs cross-validation for kernel-weighted mean estimation on a graph.
 *
 * This function implements a cross-validation procedure to evaluate the performance
 * of kernel-weighted mean estimation on a graph structure. It supports both binary
 * and continuous outcomes, using Mean Absolute Deviation (MAD) as the Rf_error metric.
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
 * @param epsilon A small value used to avoid numerical issues in calculations.
 *        Default is 1e-10.
 * @param seed Seed for the random number generator. If 0, the current time is used.
 *        Default is 0.
 *
 * @return A vector of doubles containing the average cross-validation Rf_error for each vertex.
 *         Vertices for which no Rf_error could be computed (due to being excluded in all
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
 *    c. Calculates MAD errors for the test set vertices.
 * 3. Averages the errors for each vertex across all CV iterations.
 *
 * @note
 * - The function assumes that the graph_kmean_with_weights function
 *   is available and correctly implemented.
 * - Vertices that are excluded in all CV iterations (due to graph structure or
 *   weighting) will have NaN as their final Rf_error value.
 * - The function uses C++11 random number generation facilities for reproducibility.
 * - This function uses standard (unweighted) Mean Absolute Deviation as the Rf_error metric.
 *   For a version that uses weighted MAD Rf_error, see graph_kmean_wmad_cv.
 *
 * @see graph_kmean_with_weights, graph_kmean_wmad_cv
 */
std::vector<double> graph_kmean_cv(const std::vector<std::vector<int>>& graph,
                                   const std::vector<std::vector<double>>& edge_lengths,
                                   const std::vector<double>& y,
                                   int ikernel = 1,
                                   double dist_normalization_factor = 1.01,
                                   int n_CVs = 0,
                                   int n_CV_folds = 10,
                                   unsigned int seed = 0) {

    int n_vertices = static_cast<int>(y.size());

    #if 0
    if (graph.empty()) {
        Rf_error("Input graph is empty");
    }
    if (graph.size() != n_vertices) {
        Rf_error("Number of vertices in graph (%d) doesn't match length of y (%d)",
                 (int)graph.size(), n_vertices);
    }
    // Check that each vertex's adjacency list exists
    for (int i = 0; i < n_vertices; ++i) {
        if (graph[i].empty()) {
            Rf_warning("Vertex %d has no neighbors", i);
        }
    }
    #endif

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
        auto res = graph_kmean_with_weights(graph,
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

        // Computing cross-validation Rf_error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double Rf_error = std::abs(Ecv_y[vertex] - y[vertex]);

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = Rf_error;
            } else {
                cv_error[vertex] += Rf_error;
            }
            cv_error_count[vertex]++;
        }
    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV Rf_error, leaving NaN for vertices with no estimates
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
 * @param s_ikernel An R integer specifying the kernel function to use for weighting.
 * @param s_dist_normalization_factor An R numeric value used to normalize
 *        distances in the graph.
 * @param s_n_CVs An R integer specifying the number of cross-validation
 *        iterations to perform.
 * @param s_n_CV_folds An R integer specifying the number of folds to use in
 *        each cross-validation iteration.
 * @param s_epsilon An R numeric value used as a small constant to avoid
 *        numerical issues in Rf_error calculations.
 * @param s_seed An R integer used as a seed for the random number generator.
 *
 * @return An R numeric vector containing the cross-validation errors for
 *         each vertex in the graph.
 *
 * @details
 * The function performs the following steps:
 * 1. Converts R data structures to C++ data structures.
 * 2. Calls the C++ implementation of graph_kmean_cv.
 * 3. Converts the C++ result (cross-validation errors) back to an R numeric vector.
 *
 * This function uses Mean Absolute Deviation (MAD) as the Rf_error metric for both
 * binary and continuous outcomes.
 *
 * @note
 * - This function assumes that the input R objects are of the correct type
 *   and structure. No extensive Rf_error checking is performed.
 * - The function uses PROTECT/UNPROTECT for memory management as per R's
 *   C interface guidelines.
 * - Vertices for which no cross-validation Rf_error could be computed (due to
 *   being excluded in all iterations) will have a NaN value in the result.
 * - This function uses standard (unweighted) Mean Absolute Deviation as the Rf_error metric.
 *   For a version that uses weighted MAD Rf_error, see S_graph_kmean_wmad_cv.
 *
 * @see graph_kmean_cv
 *      R_list_of_dvectors_to_cpp_vector_of_dvectors, S_graph_kmean_wmad_cv
 */
SEXP S_graph_kmean_cv(SEXP s_graph,
                      SEXP s_edge_lengths,
                      SEXP s_y,
                      SEXP s_ikernel,
                      SEXP s_dist_normalization_factor,
                      SEXP s_n_CVs,
                      SEXP s_n_CV_folds,
                      SEXP s_seed) {

    std::vector<std::vector<int>> graph           = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    std::vector<double> cv_errors =  graph_kmean_cv(graph,
                                                    edge_lengths,
                                                    y,
                                                    ikernel,
                                                    dist_normalization_factor,
                                                    n_CVs,
                                                    n_CV_folds,
                                                    seed);
    // Convert result to SEXP and return
    SEXP s_cv_errors = PROTECT(Rf_allocVector(REALSXP, cv_errors.size()));
    for (size_t i = 0; i < cv_errors.size(); ++i) {
        REAL(s_cv_errors)[i] = cv_errors[i];
    }
    UNPROTECT(1);
    return s_cv_errors;
}

/**
 * @brief Computes the kernel-weighted nearest neighbor mean of a function over vertices of a graph, with support for Bayesian Bootstrap weights.
 *
 * This function calculates the kernel-weighted mean of a function y over the vertices of a graph,
 * using edge lengths to determine the kernel weights and incorporating vertex-specific weights
 * for Bayesian Bootstrap applications. The weights are derived from Bayesian bootstrap, ensuring
 * they are strictly positive and sum to 1.
 *
 * For each vertex i, the weighted mean is computed as:
 * kmean[i] = sum(weights[j] * kernel_weights[j] * y[j]) / sum(weights[j] * kernel_weights[j])
 * where j runs over i and its neighbors, and kernel_weights are determined by the normalized distances.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains
 *              the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths of edges corresponding to the
 *                    neighbors in the graph structure.
 * @param weights A vector of strictly positive weights derived from Bayesian bootstrap that sum to 1.
 * @param y A vector of function values at each vertex.
 * @param ikernel An integer specifying the kernel function to use.
 * @param dist_normalization_factor A factor used to normalize distances. Default is 1.01.
 * @param epsilon Threshold for checking effectively zero weight sums. Default is 1e-15.
 *
 * @return A vector of kernel-weighted means for each vertex in the graph.
 *
 * @note Input validation is performed at the R level before calling this function.
 * @note The function assumes that each vertex has at least one neighbor.
 * @note For each vertex, distances are normalized by max_distance * dist_normalization_factor.
 *       If all distances for a vertex are 0, normalization uses max_dist = 1.
 * @note The kernel function should be initialized before calling this function using initialize_kernel().
 *
 * @see initialize_kernel(), kernel_fn(), C_runif_simplex()
 */
std::vector<double> graph_kmean_with_bb_weigths(const std::vector<std::vector<int>>& graph,
                                                const std::vector<std::vector<double>>& edge_lengths,
                                                const std::vector<double>& weights,
                                                const std::vector<double>& y,
                                                int ikernel,
                                                double dist_normalization_factor = 1.01,
                                                double epsilon = 1e-15) {

    auto kmean = std::vector<double>(y.size(), 0.0);
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and distances vectors
    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

    for (size_t i = 0; i < graph.size(); ++i) {
        // Assume each vertex has at least one neighbor
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

        double weighted_sum = y[i] * weights[i] * kernel_weights[0];
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            weighted_sum += weights[neighbor] * kernel_weights[j + 1] * y[neighbor];
        }

        if (weight_sum < epsilon) {
            kmean[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        kmean[i] = weighted_sum / weight_sum;
    }

    return kmean;
}

/**
 * @brief Computes Bayesian bootstraps of graph kernel weighted means.
 *
 * This function performs Bayesian bootstrap sampling to estimate the distribution
 * of graph kernel weighted means. It uses the graph structure, edge lengths, and
 * observed values to compute kernel-weighted means for each bootstrap iteration.
 *
 * The weights for each bootstrap iteration are generated using C_runif_simplex,
 * which implements the "ordered differences" method for sampling from the
 * (K-1)-dimensional simplex. While theoretically the weights could be exactly
 * zero, in practice this is extremely unlikely due to the floating-point
 * implementation. The weights always sum to 1, and the epsilon parameter
 * provides additional protection against numerical instability in the unlikely
 * case of effectively zero weights over a given set of k-nearest neighbors.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector
 *              contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths of edges corresponding
 *                     to the neighbors in the graph structure.
 * @param y A vector of observed values at each vertex of the graph.
 * @param n_bb The number of Bayesian bootstrap iterations to perform.
 * @param n_cores Number of cores to use for parallel computation:
 *                - n_cores = 1: serial execution
 *                - n_cores > 1: parallel execution with specified number of cores
 * @param ikernel An integer specifying the kernel function to use.
 * @param dist_normalization_factor A factor used to normalize distances in the kernel
 *                                  computation (default: 1.01).
 * @param epsilon Small positive number to check for effectively zero weights sum (default: 1e-15).
 *                While weights are theoretically positive with probability 1, this parameter
 *                provides protection against numerical underflow, using the original vertex
 *                value if the sum of kernel-weighted weights falls below epsilon.
 *
 * @return A vector of n_bb vectors, where each inner vector has length equal to
 *         the number of vertices in the graph, containing the bootstrap estimates
 *         of kernel-weighted means for each vertex.
 *
 * @see graph_kmean_with_bb_weigths for the underlying weighted mean computation.
 * @see C_runif_simplex for the implementation of simplex sampling.
 *
 * @throws May throw exceptions if memory allocation fails or if input sizes are inconsistent.
 *
 * @pre The sizes of graph, edge_lengths, and y must be consistent.
 * @pre n_bb must be positive.
 * @pre ikernel must be a valid kernel function identifier.
 * @pre dist_normalization_factor must be positive.
 * @pre epsilon must be positive.
 *
 * @Rf_warning This function may be computationally intensive for large graphs or high
 *          numbers of bootstrap iterations.
 */
std::vector<std::vector<double>> graph_kmean_bb(
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    int n_bb,
    int n_cores = 1,
    int ikernel = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15)
{
    const int n_vertices = (int)y.size();
    std::vector<std::vector<double>> bb_Ey(n_bb);

#ifdef _OPENMP
    if (n_cores > 1) omp_set_num_threads(n_cores);
#endif

    // If C_runif_simplex uses R RNG, precompute weights sequentially:
    std::vector<std::vector<double>> all_w(n_bb, std::vector<double>(n_vertices));
    for (int i = 0; i < n_bb; ++i) {
        C_runif_simplex(&n_vertices, all_w[i].data());   // sequential, safe
    }

    #ifdef _OPENMP
#pragma omp parallel for if(n_cores > 1) schedule(static)
    #endif
    for (int iboot = 0; iboot < n_bb; ++iboot) {
        bb_Ey[iboot] = graph_kmean_with_bb_weigths(
            graph, edge_lengths, all_w[iboot], y, ikernel, dist_normalization_factor, epsilon);
    }

    return bb_Ey;
}

#if 0
std::vector<std::vector<double>> graph_kmean_bb(const std::vector<std::vector<int>>& graph,
                                                const std::vector<std::vector<double>>& edge_lengths,
                                                const std::vector<double>& y,
                                                int n_bb,
                                                int ikernel,
double dist_normalization_factor = 1.01,
                                                double epsilon = 1e-15) {
    int n_points = y.size();
    std::vector<double> weights(n_points);
    std::vector<std::vector<double>> bb_Ey(n_bb);

    for (int iboot = 0; iboot < n_bb; iboot++) {
        C_runif_simplex(&n_points, weights.data());
        bb_Ey[iboot] = graph_kmean_with_bb_weigths(graph,
                                                   edge_lengths,
                                                   weights,
                                                   y,
                                                   ikernel,
                                                   dist_normalization_factor,
                                                   epsilon);
    }

    return bb_Ey;
}
#endif

/**
 * @brief Performs graph kernel weighted mean model Bayesian bootstrap credible intervals estimation.
 *
 * @details This function combines two steps:
 *          1. Performs graph kernel weighted mean model Bayesian bootstrap sampling
 *          2. Calculates credible intervals from the bootstrap estimates
 *
 *          The central location (mean/median) and confidence level of the intervals
 *          can be configured via parameters.
 *
 * @param graph Adjacency matrix of the graph [n_vertices x n_vertices]
 * @param edge_lengths Edge weights/distances matrix [n_vertices x n_vertices]
 * @param y Observed values at each vertex [n_vertices]
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default 500)
 * @param use_median Use median instead of mean for central location (default: false)
 * @param n_cores Number of cores for parallel computation (default: 1)
 *                - 1: serial execution
 *                - >1: parallel execution
 * @param ikernel Kernel function selector (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 *
 * @return bb_cri_t struct containing:
 *         - bb_Ey: Central location (mean/median) for each vertex
 *         - cri_L: Lower credible interval bounds
 *         - cri_U: Upper credible interval bounds
 *
 * @throws std::invalid_argument If:
 *         - Input matrices have inconsistent dimensions
 *         - n_bb <= 0
 *         - n_cores <= 0
 *         - p is not in (0,1)
 *
 * @see graph_kmean_bb() for bootstrap computation details
 * @see bb_cri() for credible interval computation details
 */
bb_cri_t graph_kmean_bb_cri(const std::vector<std::vector<int>>& graph,
                            const std::vector<std::vector<double>>& edge_lengths,
                            const std::vector<double>& y,
                            double p = 0.95,
                            int n_bb = 500,
                            bool use_median = false,
                            int n_cores = 1,
                            int ikernel = 1,
                            double dist_normalization_factor = 1.01,
                            double epsilon = 1e-15) {

    // Perform bootstrap iterations
    std::vector<std::vector<double>> bb_Eys = graph_kmean_bb(graph,
                                                             edge_lengths,
                                                             y,
                                                             n_bb,
                                                             n_cores,
                                                             ikernel,
                                                             dist_normalization_factor,
                                                             epsilon);
    // Calculate credible intervals
    bb_cri_t results = bb_cri(bb_Eys, use_median, p);
    return results;
}

/**
 * @brief Optimizes graph k-means estimation by finding the optimal neighborhood size (h)
 *        and computing credible intervals.
 *
 * @details This function performs several steps:
 *          1. For each h in [h_min, h_max]:
 *             - Constructs h-hop neighborhood (hHN) graphs
 *             - Computes cross-validation errors
 *          2. Finds optimal h that minimizes CV Rf_error
 *          3. Computes conditional expectations using optimal h
 *          4. Optional: Computes Bayesian bootstrap credible intervals
 *
 * @param graph        A sparce graph adjacency list [n_vertices]
 * @param edge_lengths An edge lengths list [n_vertices]
 * @param y Observed values at each vertex [n_vertices]
 * @param y_true True values for Rf_error calculation [n_vertices]
 * @param h_min Minimum neighborhood size to consider (default: 2)
 * @param h_max Maximum neighborhood size to consider (default: 30)
 * @param n_CVs Number of cross-validation iterations (default: 1000)
 * @param n_CV_folds Number of folds for cross-validation (default: 10)
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 500, 0 to skip)
 * @param use_median Use median instead of mean for central tendency (default: false)
 * @param ikernel Kernel function selector (default: 1)
 * @param n_cores Number of cores for parallel computation (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param seed Random seed for reproducibility (default: 0)
 *
 * @return adaptive_nbhd_size_t struct containing:
 *         - hhn_graphs: Vector of h-hop neighborhood graphs for each h
 *         - cv_errors: Cross-validation errors for each h
 *         - opt_h: Optimal h value that minimizes CV Rf_error
 *         - opt_h_graph: Graph with optimal h
 *         - condEy: Conditional expectations using optimal h
 *         - bb_condEy: Bootstrap central tendency (if n_bb > 0)
 *         - cri_L: Lower credible interval bounds (if n_bb > 0)
 *         - cri_U: Upper credible interval bounds (if n_bb > 0)
 *
 * @throws std::invalid_argument If:
 *         - Input matrices have inconsistent dimensions
 *         - h_min > h_max
 *         - Invalid parameter values (n_CVs, n_CV_folds, p, etc.)
 */
adaptive_nbhd_size_t gkmm(const std::vector<std::vector<int>>& graph,
                          const std::vector<std::vector<double>>& edge_lengths,
                          const std::vector<double>& y,
                          const std::vector<double>& y_true,
                          bool use_median = false,
                          int h_min = 2,
                          int h_max = 30,
                          int n_CVs = 1000,
                          int n_CV_folds = 10,
                          double p = 0.95,
                          int n_bb = 500,
                          int ikernel = 1,
                          int n_cores = 1,
                          double dist_normalization_factor = 1.01,
                          double epsilon = 1e-15,
                          unsigned int seed = 0) {

    #define DEBUG__gkmm 0
    #if DEBUG__gkmm
    Rprintf("In gkmm()\n");
    print_vect_vect(graph,"graph");
    print_vect_vect(edge_lengths,"edge_lengths");
    #endif

    int n_vertices = static_cast<int>(y.size());

    adaptive_nbhd_size_t results;
    results.graphs.resize(h_max - h_min + 1);
    results.cv_errors.resize(h_max - h_min + 1);
    results.h_values.resize(h_max - h_min + 1);

    for (int i = 0, h = h_min; h <= h_max; h++, i++) {
        results.h_values[i] = h;

        // Creating an hHN graph for the given value of h
        auto hhn_graph = create_hHN_graph(graph, edge_lengths, h);

        #if DEBUG__gkmm
        Rprintf("h: %d\n", h);
        Rprintf("hhn_graph.first.size(): %d\n", (int)hhn_graph.first.size());
        Rprintf("hhn_graph.second.size(): %d\n", (int)hhn_graph.second.size());

        print_vect_vect(hhn_graph.first,"hhn_graph.first");
        print_vect_vect(hhn_graph.second,"hhn_graph.second");
        Rf_error("Debugging");
        #endif

        // Computing CV errors
        auto errors = graph_kmean_cv(hhn_graph.first,
                                     hhn_graph.second,
                                     y,
                                     ikernel,
                                     dist_normalization_factor,
                                     n_CVs,
                                     n_CV_folds,
                                     seed);

        // Calculate mean Rf_error across vertices
        double total_error = std::accumulate(errors.begin(), errors.end(), 0.0);
        results.cv_errors[i] = total_error / n_vertices;
        results.graphs[i] = std::move(hhn_graph);
    }

    // Find the optimal h (minimum CV Rf_error)
    auto min_it = std::min_element(results.cv_errors.begin(), results.cv_errors.end());
    int opt_h_idx = std::distance(results.cv_errors.begin(), min_it);
    results.opt_h = h_min + opt_h_idx;

    // Store optimal graph
    results.opt_h_graph = results.graphs[opt_h_idx];

    // Compute conditional expectations using optimal graph
    results.condEy = graph_kmean(results.opt_h_graph.first,
                                 results.opt_h_graph.second,
                                 y,
                                 ikernel,
                                 dist_normalization_factor);

    // Compute true errors
    if (!y_true.empty() && (int)y_true.size() == n_vertices) {
        results.true_errors.resize(n_vertices);
        for (int i = 0; i < n_vertices; i++) {
            results.true_errors[i] = std::abs(y_true[i] - results.condEy[i]);
        }
    } else {
        results.true_errors.clear();  // Ensure empty if no true values
    }

    // Optional: Compute bootstrap credible intervals
    if (n_bb > 0) {
        bb_cri_t bb_cri_results = graph_kmean_bb_cri(results.opt_h_graph.first,
                                                     results.opt_h_graph.second,
                                                     y,
                                                     p,
                                                     n_bb,
                                                     use_median,
                                                     n_cores,
                                                     ikernel,
                                                     dist_normalization_factor,
                                                     epsilon);

        results.bb_condEy = std::move(bb_cri_results.bb_Ey);
        results.cri_L = std::move(bb_cri_results.cri_L);
        results.cri_U = std::move(bb_cri_results.cri_U);
    }

    #if DEBUG__gkmm
    // Before returning results, verify the data is stored correctly
    Rprintf("Verification before return:\n");
    for (size_t i = 0; i < results.graphs.size(); ++i) {
        Rprintf("h = %d:\n", h_min + static_cast<int>(i));
        Rprintf("  Graph size: %zu\n", results.graphs[i].first.size());
        Rprintf("  CV Rf_error: %f\n", results.cv_errors[i]);
    }
    #endif

    return results;
}

/**
 * @brief Performs adaptive neighborhood size k-means estimation for univariate data.
 *
 * @details This function performs nonparametric regression using graph k-means with
 *          adaptive neighborhood size selection for univariate data (x,y). It:
 *          1. Constructs a chain graph from sorted x values as adjacency lists
 *          2. Computes corresponding edge lengths lists
 *          3. Applies adaptive neighborhood size graph k-means to estimate E(y|x)
 *
 * @pre x must be sorted in ascending order
 *
 * @param x Sorted predictor values [n_points]
 * @param y Response values corresponding to x [n_points]
 * @param y_true True response values for Rf_error calculation [n_points]
 * @param use_median Use median instead of mean for central tendency (default: false)
 * @param h_min Minimum neighborhood size to consider (default: 2)
 * @param h_max Maximum neighborhood size to consider (default: 30)
 * @param n_CVs Number of cross-validation iterations (default: 1000)
 * @param n_CV_folds Number of folds for cross-validation (default: 10)
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 500)
 * @param ikernel Kernel function selector (default: 1)
 * @param n_cores Number of cores for parallel computation (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param seed Random seed for reproducibility (default: 0)
 *
 * @return adaptive_nbhd_size_t struct containing results from graph k-means analysis
 *
 * @throws std::invalid_argument If:
 *         - Input vectors have different sizes
 *         - x is not sorted
 *         - Less than 2 vertices
 *         - Invalid parameter values
 */
adaptive_nbhd_size_t
univariate_gkmm(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     const std::vector<double>& y_true,
                                     bool use_median = false,
                                     int h_min = 2,
                                     int h_max = 30,
                                     int n_CVs = 1000,
                                     int n_CV_folds = 10,
                                     double p = 0.95,
                                     int n_bb = 500,
                                     int ikernel = 1,
                                     int n_cores = 1,
                                     double dist_normalization_factor = 1.01,
                                     double epsilon = 1e-15,
                                     unsigned int seed = 0) {

    auto [x_graph, x_edge_lengths] = create_chain_graph(x);

    // Apply adaptive neighborhood size graph k-means
    return gkmm(x_graph,
                x_edge_lengths,
                y,
                y_true,
                use_median,
                h_min,
                h_max,
                n_CVs,
                n_CV_folds,
                p,
                n_bb,
                ikernel,
                n_cores,
                dist_normalization_factor,
                epsilon,
                seed);
}

/**
 * @brief R interface for adaptive neighborhood size k-means on univariate data
 *
 * @details Performs k-means clustering with adaptive neighborhood sizes on univariate data.
 *          The function:
 *          1. Constructs h-hop neighborhood graphs for a range of h values
 *          2. Performs cross-validation to select optimal neighborhood size
 *          3. Computes conditional expectations and credible intervals via bootstrap
 *          4. Evaluates true errors if ground truth is provided
 *
 * @param s_x SEXP containing sorted numeric vector of data points
 * @param s_y SEXP containing numeric vector of response values
 * @param s_y_true SEXP containing numeric vector of true response values (optional)
 * @param s_use_median SEXP logical; TRUE to use median, FALSE to use mean
 * @param s_h_min SEXP integer; minimum neighborhood size (h) to consider
 * @param s_h_max SEXP integer; maximum neighborhood size (h) to consider
 * @param s_n_CVs SEXP integer; number of cross-validation iterations
 * @param s_n_CV_folds SEXP integer; number of cross-validation folds
 * @param s_p SEXP numeric; confidence level for credible intervals (0,1)
 * @param s_n_bb SEXP integer; number of bootstrap iterations
 * @param s_ikernel SEXP integer; kernel function selector (1-4)
 * @param s_n_cores SEXP integer; number of cores for parallel processing
 * @param s_dist_normalization_factor SEXP numeric; scaling factor for distances
 * @param s_epsilon SEXP numeric; numerical stability parameter
 * @param s_seed SEXP integer; random seed for reproducibility
 *
 * @return SEXP containing a named list with components:
 *   - h_values: Vector of tested neighborhood sizes
 *   - hHN_graphs: List of h-hop neighborhood graphs, each containing:
 *     * adj_list: Adjacency list representation
 *     * edge_lengths: Corresponding edge weights
 *   - cv_errors: Cross-validation errors for each h
 *   - opt_h: Optimal neighborhood size selected by CV
 *   - opt_h_graph_adj_list: Adjacency list for optimal h
 *   - opt_h_graph_edge_lengths: Edge weights for optimal h
 *   - condEy: Conditional expectations at optimal h
 *   - bb_condEy: Bootstrap samples of conditional expectations
 *   - cri_L: Lower bounds of credible intervals
 *   - cri_U: Upper bounds of credible intervals
 *   - true_error: True Rf_error at optimal h (if y_true provided)
 *
 * @note Input vectors s_x and s_y must have the same length.
 *       The s_x vector must be sorted in ascending order.
 *       If provided, s_y_true must match the length of s_x.
 *
 * @throws R Rf_error if input validation fails or memory allocation errors occur
 */
SEXP S_univariate_gkmm(SEXP s_x,
                       SEXP s_y,
                       SEXP s_y_true,
                       SEXP s_use_median,
                       SEXP s_h_min,
                       SEXP s_h_max,
                       SEXP s_n_CVs,
                       SEXP s_n_CV_folds,
                       SEXP s_p,
                       SEXP s_n_bb,
                       SEXP s_ikernel,
                       SEXP s_n_cores,
                       SEXP s_dist_normalization_factor,
                       SEXP s_epsilon,
                       SEXP s_seed) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(s_y_true) == n_points) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    bool use_median = (LOGICAL(s_use_median)[0] == 1);
    int h_min = INTEGER(s_h_min)[0];
    int h_max = INTEGER(s_h_max)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double p = REAL(s_p)[0];
    int n_bb = INTEGER(s_n_bb)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    auto cpp_results = univariate_gkmm(x,
                                       y,
                                       y_true,
                                       use_median,
                                       h_min,
                                       h_max,
                                       n_CVs,
                                       n_CV_folds,
                                       p,
                                       n_bb,
                                       ikernel,
                                       n_cores,
                                       dist_normalization_factor,
                                       epsilon,
                                       seed);

    // Construct s_hHN_graphs
    SEXP s_hHN_graphs = PROTECT(Rf_allocVector(VECSXP, h_max - h_min + 1)); n_protected++;
    SEXP h_names = PROTECT(Rf_allocVector(STRSXP, h_max - h_min + 1)); n_protected++;

    for (int i = 0; i < h_max - h_min + 1; i++) {
        SEXP h_graph = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP h_graph_names = PROTECT(Rf_allocVector(STRSXP, 2));

        SEXP adj_list = convert_vector_vector_int_to_R(cpp_results.graphs[i].first); UNPROTECT(1);
        SEXP edge_lengths = convert_vector_vector_double_to_R(cpp_results.graphs[i].second); UNPROTECT(1);

        SET_VECTOR_ELT(h_graph, 0, adj_list);
        SET_VECTOR_ELT(h_graph, 1, edge_lengths);

        SET_STRING_ELT(h_graph_names, 0, Rf_mkChar("adj_list"));
        SET_STRING_ELT(h_graph_names, 1, Rf_mkChar("edge_lengths"));
        Rf_setAttrib(h_graph, R_NamesSymbol, h_graph_names);

        SET_VECTOR_ELT(s_hHN_graphs, i, h_graph);

        char h_label[32];
        snprintf(h_label, sizeof(h_label), "h_%d", h_min + i);
        SET_STRING_ELT(h_names, i, Rf_mkChar(h_label));

        UNPROTECT(2);  // h_graph and h_graph_names
    }
    Rf_setAttrib(s_hHN_graphs, R_NamesSymbol, h_names);

    // Creating return list
    const int N_COMPONENTS = 11;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    SEXP s_h_values = convert_vector_int_to_R(cpp_results.h_values); n_protected++;
    SET_VECTOR_ELT(result, 0, s_h_values);
    SET_VECTOR_ELT(result, 1, s_hHN_graphs);

    if (!cpp_results.cv_errors.empty() && n_points > 0) {
        SEXP s_cv_errors = convert_vector_double_to_R(cpp_results.cv_errors); n_protected++;
        SET_VECTOR_ELT(result, 2, s_cv_errors);
    } else {
        SET_VECTOR_ELT(result, 2, R_NilValue);
    }

    SEXP s_opt_h = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h)[0] = cpp_results.opt_h;
    SET_VECTOR_ELT(result, 3, s_opt_h);

    SEXP s_opt_h_graph_adj_list = convert_vector_vector_int_to_R(cpp_results.opt_h_graph.first); UNPROTECT(1);
    SET_VECTOR_ELT(result, 4, s_opt_h_graph_adj_list);

    SEXP s_opt_h_graph_edge_lengths = convert_vector_vector_double_to_R(cpp_results.opt_h_graph.second); UNPROTECT(1);
    SET_VECTOR_ELT(result, 5, s_opt_h_graph_edge_lengths);

    SEXP s_condEy = convert_vector_double_to_R(cpp_results.condEy); n_protected++;
    SET_VECTOR_ELT(result, 6, s_condEy);

    if (cpp_results.bb_condEy.size() > 0) {
        SEXP s_bb_condEy = convert_vector_double_to_R(cpp_results.bb_condEy); n_protected++;
        SET_VECTOR_ELT(result, 7, s_bb_condEy);

        SEXP s_cri_L = convert_vector_double_to_R(cpp_results.cri_L); n_protected++;
        SET_VECTOR_ELT(result, 8, s_cri_L);

        SEXP s_cri_U = convert_vector_double_to_R(cpp_results.cri_U); n_protected++;
        SET_VECTOR_ELT(result, 9, s_cri_U);
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }

    if (cpp_results.true_errors.size() > 0) {
        SEXP s_true_error = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        double mean_true_error = std::accumulate(cpp_results.true_errors.begin(),
                                                 cpp_results.true_errors.end(), 0.0) /  cpp_results.true_errors.size();
        REAL(s_true_error)[0] = mean_true_error;
        SET_VECTOR_ELT(result, 10, s_true_error);
    } else {
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    // Setting names for return list
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, Rf_mkChar("h_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("graphs"));
    SET_STRING_ELT(names, 2, Rf_mkChar("h_cv_errors"));
    SET_STRING_ELT(names, 3, Rf_mkChar("opt_h"));
    SET_STRING_ELT(names, 4, Rf_mkChar("opt_graph_adj_list"));
    SET_STRING_ELT(names, 5, Rf_mkChar("opt_graph_edge_lengths"));
    SET_STRING_ELT(names, 6, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 7, Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 8, Rf_mkChar("opt_ci_lower"));
    SET_STRING_ELT(names, 9, Rf_mkChar("opt_ci_upper"));
    SET_STRING_ELT(names, 10, Rf_mkChar("true_error"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}
