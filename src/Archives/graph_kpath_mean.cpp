#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <algorithm> // for std::max
#include <random>

#include "cpp_utils.h"
#include "SEXP_cpp_conversion_utils.h"
#include "stats_utils.h"
#include "kernels.h"

std::vector<int> shortest_path_by_hops(const std::vector<std::vector<int>>& graph,
                                       int start, int end);

extern "C" {
    SEXP S_graph_kpath_mean(SEXP s_graph,
                            SEXP s_edge_lengths,
                            SEXP s_core_graph,
                            SEXP s_y,
                            SEXP s_ikernel,
                            SEXP s_dist_normalization_factor);
}


/**
 * @brief Computes a kernel-weighted mean of vertex values using path-averaged output values
 *
 * @details This function implements a modified kernel smoothing algorithm for graph-structured data
 * where the neighbor values are computed as means along shortest paths. For each vertex,
 * it computes a weighted average of:
 * 1. The vertex's own value
 * 2. The mean values along shortest paths to each of its neighbors (excluding the initial vertex)
 *
 * The weights are determined by a kernel function applied to the normalized edge distances.
 * The algorithm proceeds as follows:
 * 1. For each vertex i:
 *    - Normalizes the distances to its neighbors
 *    - Computes kernel weights based on these distances
 *    - For each neighbor j:
 *      * Finds the shortest path from i to j
 *      * Computes the mean of values along this path (excluding vertex i)
 *      * Weights this path mean by the kernel weight
 *    - Combines all weighted values to produce the final smoothed value
 *
 * @param graph Vector of vectors representing the graph adjacency list
 * @param edge_lengths Vector of vectors containing the corresponding edge lengths
 * @param core_graph Graph structure used for shortest path computation (typically unweighted)
 * @param y Vector of vertex values to be smoothed
 * @param ikernel Integer specifying which kernel function to use
 * @param dist_normalization_factor Factor used to normalize distances (default: 1.01)
 *
 * @return Vector of smoothed values for each vertex
 *
 * @throws std::invalid_argument If the shortest path computation encounters invalid vertices
 *
 * @note
 * - Vertices with no neighbors retain their original values
 * - The function assumes the existence of an initialized kernel function (initialize_kernel)
 * - Path means exclude the initial vertex to avoid self-bias
 * - Edge distances are normalized by the maximum distance times the normalization factor
 *
 * @pre
 * - graph, edge_lengths, and y must have consistent dimensions
 * - core_graph must represent a valid graph structure for path finding
 * - ikernel must be a valid kernel function identifier
 *
 * @see graph_kernel_weighted_mean, shortest_path_by_hops, initialize_kernel
 *
 * Time Complexity: O(V * (E + (V + E) * P)) where:
 * - V is the number of vertices
 * - E is the maximum number of edges per vertex
 * - P is the average path length
 */
std::vector<double> graph_kpath_mean(const std::vector<std::vector<int>>& graph,
                                                    const std::vector<std::vector<double>>& edge_lengths,
                                                    const std::vector<std::vector<int>>& core_graph,
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

                std::vector<int> path = shortest_path_by_hops(core_graph, i, neighbor);
                if (path.empty()) continue;

                double y_mean_along_path = 0;
                int path_len = path.size();
                for (int l = 1; l < path_len; ++l) {
                    y_mean_along_path += y[path[l]];
                }
                y_mean_along_path /= path_len - 1;

                weighted_sum += kernel_weights[j + 1] * y_mean_along_path;
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
 * @brief R interface wrapper for graph_kpath_mean function
 *
 * @details This function provides an R-callable wrapper for the C++ implementation
 * of the graph kernel weighted path mean algorithm. It handles the conversion between
 * R data structures (SEXP) and C++ types, manages R memory protection, and ensures
 * proper data transfer between R and C++.
 *
 * The function performs the following steps:
 * 1. Converts R objects to C++ data structures
 * 2. Calls the core C++ implementation
 * 3. Converts the result back to an R object
 *
 * @param s_graph SEXP containing the graph structure as an R list of integer vectors
 *                representing adjacency lists (1-based indices in R converted to 0-based in C++)
 * @param s_edge_lengths SEXP containing edge lengths as an R list of numeric vectors
 *                      corresponding to the edges in s_graph
 * @param s_core_graph SEXP containing the core graph structure used for path finding
 *                     as an R list of integer vectors
 * @param s_y SEXP containing the numeric vector of vertex values to be smoothed
 * @param s_ikernel SEXP containing a single integer specifying the kernel function to use
 * @param s_dist_normalization_factor SEXP containing a single numeric value for distance normalization
 *
 * @return SEXP containing a numeric vector of smoothed values for each vertex
 *
 * @note
 * - The function assumes R indexing (1-based) and converts to C++ indexing (0-based) internally
 * - Memory protection is handled using PROTECT/UNPROTECT macros
 * - The returned SEXP is properly protected before return
 *
 * @pre
 * - s_graph must be a list of integer vectors
 * - s_edge_lengths must be a list of numeric vectors
 * - s_core_graph must be a list of integer vectors
 * - s_y must be a numeric vector
 * - s_ikernel must be a single integer
 * - s_dist_normalization_factor must be a single numeric value
 * - The dimensions of s_graph and s_edge_lengths must match
 * - The length of s_y must match the number of vertices in the graph
 *
 * @throws R error via Rf_error() if input validation fails
 *
 * @see graph_kpath_mean, Rgraph_to_vector, Rweights_to_vector
 *
 * Example R usage:
 * ```r
 * result <- .Call("S_graph_kpath_mean",
 *                 graph,
 *                 edge_lengths,
 *                 core_graph,
 *                 y,
 *                 as.integer(kernel_type),
 *                 as.numeric(dist_norm_factor))
 * ```
 */
SEXP S_graph_kpath_mean(SEXP s_graph,
                        SEXP s_edge_lengths,
                        SEXP s_core_graph,
                        SEXP s_y,
                        SEXP s_ikernel,
                        SEXP s_dist_normalization_factor) {

    if (TYPEOF(s_y) != REALSXP) {
        Rf_error("y must be a numeric vector");
    }
    if (TYPEOF(s_ikernel) != INTSXP) {
        Rf_error("ikernel must be an integer");
    }
    if (TYPEOF(s_dist_normalization_factor) != REALSXP) {
        Rf_error("dist_normalization_factor must be numeric");
    }

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    std::vector<double> result = graph_kpath_mean(graph,
                                                  edge_lengths,
                                                  core_graph,
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
 * @brief Computes the kernel-weighted path-based mean of a function over vertices of a graph, with support for vertex weighting.
 *
 * This function calculates the kernel-weighted mean of values associated with vertices in a graph,
 * incorporating path information from a core graph structure. For each vertex pair, it considers the
 * mean value along the shortest path between them rather than just direct neighbor values.
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains
 *              the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths (or distances) of edges in the graph.
 *                     It should have the same structure as the 'graph' parameter.
 * @param core_graph A vector of vectors representing the core graph structure used for finding
 *                   shortest paths between vertices.
 * @param weights A vector of weights for each vertex in the graph. A weight of 0 will cause
 *                the vertex to be excluded from calculations.
 * @param y A vector of values associated with each vertex in the graph.
 * @param ikernel An integer specifying the kernel function to use for weighting.
 * @param dist_normalization_factor A factor used to normalize distances in the graph.
 *                                  Default value is 1.01.
 *
 * @return A std::pair containing:
 *         - First: A vector of doubles representing the path-based kernel-weighted mean for each vertex.
 *         - Second: A vector of integers containing the indices of excluded vertices (those with zero total weight).
 *
 * @note Vertices with no neighbors (empty adjacency list in 'graph') are handled specially:
 *       - If the vertex's weight is 0, it's added to the excluded vertices list.
 *       - If the vertex's weight is non-zero, its original value from 'y' is used.
 *
 * @warning This function assumes that the input vectors (graph, edge_lengths, core_graph, weights, y)
 *          all have appropriate sizes, and that the core_graph structure is valid for path finding.
 *
 * @pre The function 'initialize_kernel' must be called before this function to set up the kernel function.
 * @pre The function 'kernel_fn' must be properly defined to calculate kernel weights based on distances.
 * @pre The function 'shortest_path_by_hops' must be properly defined to find paths in the core graph.
 *
 * @throws May throw exceptions if vector accesses are out of bounds or if mathematical operations
 *         (like division by zero) fail. Proper input validation should be done before calling this function.
 *
 * @see initialize_kernel, kernel_fn, shortest_path_by_hops
 */
std::pair<std::vector<double>, std::vector<int>> graph_kpath_mean_with_weights(
        const std::vector<std::vector<int>>& graph,
        const std::vector<std::vector<double>>& edge_lengths,
        const std::vector<std::vector<int>>& core_graph,
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

            // Calculate total weight first to check for exclusion
            double weight_sum = weights[i] * kernel_weights[0];
            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                std::vector<int> path = shortest_path_by_hops(core_graph, i, neighbor);
                if (!path.empty()) {
                    weight_sum += weights[neighbor] * kernel_weights[j + 1];
                }
            }

            if (weight_sum == 0) {
                excluded_vertices.push_back(i);
                kmean[i] = 0;
            } else {
                double weighted_sum = weights[i] * y[i] * kernel_weights[0];
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    int neighbor = graph[i][j];
                    std::vector<int> path = shortest_path_by_hops(core_graph, i, neighbor);
                    if (path.empty()) continue;

                    double y_mean_along_path = 0;
                    int path_len = path.size();
                    for (int l = 1; l < path_len; ++l) {
                        y_mean_along_path += y[path[l]];
                    }
                    y_mean_along_path /= path_len - 1;

                    weighted_sum += weights[neighbor] * kernel_weights[j + 1] * y_mean_along_path;
                }
                kmean[i] = weighted_sum / weight_sum;
            }
        } else {
            if (weights[i] == 0) {
                excluded_vertices.push_back(i);
                kmean[i] = 0;
            } else {
                kmean[i] = weights[i] * y[i];
            }
        }
    }

    return std::make_pair(kmean, excluded_vertices);
}


/**
 * @brief Performs cross-validation for path-based kernel-weighted means on a graph.
 *
 * This function estimates the prediction error of path-based kernel-weighted means through
 * cross-validation. For each round of cross-validation, it randomly selects test vertices,
 * computes path-based means excluding these vertices, and calculates the absolute deviation
 * between predicted and actual values.
 *
 * @param graph A vector of vectors representing the graph structure where graph[i] contains
 *              indices of vertices adjacent to vertex i.
 * @param edge_lengths A vector of vectors containing the lengths/weights of edges where
 *                    edge_lengths[i][j] is the length of edge from vertex i to its jth neighbor.
 * @param core_graph A vector of vectors representing the core graph structure used for finding
 *                   shortest paths between vertices.
 * @param y A vector of values associated with each vertex that will be predicted.
 * @param ikernel Integer specifying which kernel function to use for weighting (default: 1).
 * @param dist_normalization_factor Factor used to normalize distances (default: 1.01).
 * @param n_CVs Number of cross-validation rounds to perform (must be > 0).
 * @param n_CV_folds Number of folds for cross-validation (default: 10).
 * @param epsilon Small value to prevent numerical instability (default: 1e-10).
 * @param seed Random seed for test set selection (default: 0, uses system time).
 *
 * @return Vector of doubles containing the average cross-validation error for each vertex.
 *         Vertices that were never in a test set or had no valid predictions will have NaN values.
 *
 * @note The function uses absolute deviation as the error metric.
 * @note If seed=0, the current system time is used as the random seed.
 * @note For vertices with no neighbors or those that were always excluded, the result will be NaN.
 *
 * @throws std::invalid_argument if n_CVs <= 0
 *
 * @pre The graph, edge_lengths, core_graph, and y vectors must have consistent sizes
 * @pre The function assumes initialize_kernel and kernel_fn are properly set up
 * @pre shortest_path_by_hops must be properly defined for path finding in core_graph
 *
 * @see graph_kpath_mean_with_weights, shortest_path_by_hops
 */
std::vector<double> graph_kpath_mean_cv(const std::vector<std::vector<int>>& graph,
                                        const std::vector<std::vector<double>>& edge_lengths,
                                        const std::vector<std::vector<int>>& core_graph,
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

        // Estimating the conditional expectation of cv_y using path-based means
        auto res = graph_kpath_mean_with_weights(graph,
                                               edge_lengths,
                                               core_graph,
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

        // Computing cross-validation error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double error = std::abs(Ecv_y[vertex] - y[vertex]);

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = error;
            } else {
                cv_error[vertex] += error;
            }
            cv_error_count[vertex]++;
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


SEXP S_graph_kpath_mean_cv(SEXP s_graph,
                           SEXP s_edge_lengths,
                           SEXP s_core_graph,
                           SEXP s_y,
                           SEXP s_ikernel,
                           SEXP s_dist_normalization_factor,
                           SEXP s_n_CVs,
                           SEXP s_n_CV_folds,
                           SEXP s_epsilon,
                           SEXP s_seed) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    std::vector<double> cv_errors =  graph_kpath_mean_cv(graph,
                                                         edge_lengths,
                                                         core_graph,
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
