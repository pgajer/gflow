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

extern "C" {
    SEXP S_graph_cached_kpath_mean(SEXP s_graph,
                                   SEXP s_edge_lengths,
                                   SEXP s_core_graph,
                                   SEXP s_shortest_paths,
                                   SEXP s_y,
                                   SEXP s_ikernel,
                                   SEXP s_dist_normalization_factor);

    SEXP S_graph_cached_kpath_mean_cv(SEXP s_graph,
                                      SEXP s_edge_lengths,
                                      SEXP s_core_graph,
                                      SEXP s_shortest_paths,
                                      SEXP s_y,
                                      SEXP s_ikernel,
                                      SEXP s_dist_normalization_factor,
                                      SEXP s_n_CVs,
                                      SEXP s_n_CV_folds,
                                      SEXP s_epsilon,
                                      SEXP s_seed);
}

/**
 * @brief Computes kernel-weighted means of vertex values using pre-computed shortest paths
 *
 * @details This function implements a modified kernel smoothing algorithm for graph-structured data
 * where neighbor values are computed as means along pre-computed shortest paths. For each vertex,
 * it computes a weighted average of:
 * 1. The vertex's own value
 * 2. The mean values along shortest paths to each of its neighbors (excluding the initial vertex)
 *
 * The weights are determined by a kernel function applied to the normalized edge distances.
 * The algorithm uses pre-computed shortest paths stored in a map for efficiency.
 *
 * The algorithm proceeds as follows:
 * 1. For each vertex i:
 *    - Normalizes the distances to its neighbors using the maximum distance
 *    - Computes kernel weights based on these normalized distances
 *    - For each neighbor j:
 *      * Retrieves the pre-computed shortest path from i to j
 *      * Computes the mean of values along this path (excluding vertex i)
 *      * Weights this path mean by the kernel weight
 *    - Combines all weighted values to produce the final smoothed value
 *
 * @param graph Vector of vectors representing the graph adjacency list
 * @param edge_lengths Vector of vectors containing the corresponding edge lengths
 * @param core_graph Graph structure used for path validation
 * @param shortest_paths Map containing pre-computed shortest paths between vertex pairs.
 *                      Keys are pairs (source, target) where source < target
 * @param y Vector of vertex values to be smoothed
 * @param ikernel Integer specifying which kernel function to use
 * @param dist_normalization_factor Factor used to normalize distances (default: 1.01)
 *
 * @return Vector of smoothed values for each vertex
 *
 * @throws R error if a required shortest path is not found in the shortest_paths map
 * @throws R error if a retrieved path is empty
 *
 * @note
 * - Vertices with no neighbors retain their original values
 * - The function assumes the existence of an initialized kernel function (initialize_kernel)
 * - Path means exclude the initial vertex to avoid self-bias
 * - Edge distances are normalized by the maximum distance times the normalization factor
 * - The shortest_paths map must contain paths for all connected vertex pairs
 * - For any vertex pair (i,j), the key in shortest_paths must be (min(i,j), max(i,j))
 *
 * @pre
 * - graph, edge_lengths, and y must have consistent dimensions
 * - core_graph must represent a valid graph structure
 * - shortest_paths must contain all necessary paths
 * - ikernel must be a valid kernel function identifier
 *
 * @see graph_kpath_mean, initialize_kernel
 *
 * Time Complexity: O(V * E) where:
 * - V is the number of vertices
 * - E is the maximum number of edges per vertex
 * Note: Path lookup is O(1) due to pre-computation
 *
 * Space Complexity: O(V^2) for the shortest_paths map in worst case
 *
 * @example
 * ```cpp
 * std::vector<std::vector<int>> graph = {{1,2}, {0,2}, {0,1}};
 * std::vector<std::vector<double>> lengths = {{1.0,2.0}, {1.0,1.5}, {2.0,1.5}};
 * std::vector<std::vector<int>> core_graph = graph;
 * std::map<std::pair<int,int>, std::vector<int>> paths;
 * paths[{0,1}] = {0,1};
 * paths[{0,2}] = {0,2};
 * paths[{1,2}] = {1,2};
 * std::vector<double> y = {1.0, 2.0, 3.0};
 * auto result = graph_cached_kpath_mean(graph, lengths, core_graph,
 *                                                     paths, y, 1);
 * ```
 */
std::vector<double> graph_cached_kpath_mean(const std::vector<std::vector<int>>& graph,
                                            const std::vector<std::vector<double>>& edge_lengths,
                                            const std::vector<std::vector<int>>& core_graph,
                                            const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
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

                auto path_key = std::make_pair(
                    std::min<size_t>(i, neighbor),
                    std::max<size_t>(i, neighbor)
                    );

                auto path_it = shortest_paths.find(path_key);
                if (path_it == shortest_paths.end()) {
                    Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                            static_cast<int>(i), neighbor);
                    error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                }

                const auto& path = path_it->second;
                if (path.empty()) {
                    Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                            static_cast<int>(i), neighbor);
                    error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                }

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
 * @brief R interface for computing kernel-weighted means of vertex values using pre-computed shortest paths
 *
 * @details This function provides an R interface to graph_cached_kpath_mean. Input validation
 * is handled by the calling R function graph.cached.kpath.mean.
 *
 * @param s_graph SEXP representing an R list of integer vectors defining the graph adjacency list
 * @param s_edge_lengths SEXP representing an R list of numeric vectors containing edge lengths
 * @param s_core_graph SEXP representing an R list of integer vectors defining the core graph structure
 * @param s_shortest_paths SEXP representing an R list of integer vectors containing shortest paths
 * @param s_y SEXP representing an R numeric vector of vertex values to be smoothed
 * @param s_ikernel SEXP representing an R integer specifying the kernel function to use
 * @param s_dist_normalization_factor SEXP representing an R numeric scalar for distance normalization
 *
 * @return SEXP A numeric vector containing the smoothed values for each vertex
 *
 * @note
 * - All vertex indices in R inputs are assumed to be 0-based (conversion done in R)
 * - The function handles R object protection internally
 *
 * @see graph_cached_kpath_mean, graph.cached.kpath.mean
 */
SEXP S_graph_cached_kpath_mean(SEXP s_graph,
                               SEXP s_edge_lengths,
                               SEXP s_core_graph,
                               SEXP s_shortest_paths,
                               SEXP s_y,
                               SEXP s_ikernel,
                               SEXP s_dist_normalization_factor) {
    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths = shortest_paths_Rlist_to_cpp_map(s_shortest_paths);
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    // Compute result
    std::vector<double> result = graph_cached_kpath_mean(graph,
                                                         edge_lengths,
                                                         core_graph,
                                                         shortest_paths,
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
 * @brief Computes the kernel-weighted path-based mean using pre-computed shortest paths, with support for vertex weighting.
 *
 * This function calculates the kernel-weighted mean of values associated with vertices in a graph,
 * using pre-computed shortest paths between vertices for efficiency. It supports vertex weighting
 * and handles excluded vertices (those with zero total weight).
 *
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains
 *              the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors containing the lengths (or distances) of edges in the graph.
 *                     It should have the same structure as the 'graph' parameter.
 * @param core_graph A vector of vectors representing the core graph structure.
 * @param shortest_paths A map containing pre-computed shortest paths between vertex pairs.
 *                      Keys are pairs of vertices (i,j) where i < j.
 * @param weights A vector of weights for each vertex in the graph.
 * @param y A vector of values associated with each vertex in the graph.
 * @param ikernel An integer specifying the kernel function to use for weighting.
 * @param dist_normalization_factor A factor used to normalize distances in the graph.
 *                                  Default value is 1.01.
 *
 * @return A std::pair containing:
 *         - First: A vector of doubles representing the path-based kernel-weighted mean for each vertex.
 *         - Second: A vector of integers containing the indices of excluded vertices (those with zero total weight).
 *
 * @note This version uses cached shortest paths for efficiency compared to computing them on-the-fly.
 * @throws Errors if required paths are missing from the shortest_paths map.
 *
 * @pre All input vectors must have appropriate sizes
 * @pre shortest_paths must contain all required vertex pairs
 * @pre initialize_kernel must be called before this function
 */
std::pair<std::vector<double>, std::vector<int>> graph_cached_kpath_mean_with_weights(
        const std::vector<std::vector<int>>& graph,
        const std::vector<std::vector<double>>& edge_lengths,
        const std::vector<std::vector<int>>& core_graph,
        const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
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

            // Calculate weighted sum and check for exclusion
            double weight_sum = weights[i] * kernel_weights[0];
            double weighted_sum = weights[i] * y[i] * kernel_weights[0];

            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                auto path_key = std::make_pair(
                    std::min<size_t>(i, neighbor),
                    std::max<size_t>(i, neighbor)
                );

                auto path_it = shortest_paths.find(path_key);
                if (path_it == shortest_paths.end()) {
                    Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                            static_cast<int>(i), neighbor);
                    error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                }

                const auto& path = path_it->second;
                if (path.empty()) {
                    Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                            static_cast<int>(i), neighbor);
                    error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                }

                // Calculate weighted mean along path
                double path_weighted_sum = 0.0;
                double path_weight_sum = 0.0;
                for (int l = 1; l < path.size(); ++l) {
                    path_weighted_sum += weights[path[l]] * y[path[l]];
                    path_weight_sum += weights[path[l]];
                }

                if (path_weight_sum > 0) {
                    double y_mean_along_path = path_weighted_sum / path_weight_sum;
                    weighted_sum += kernel_weights[j + 1] * y_mean_along_path;
                    weight_sum += kernel_weights[j + 1];
                }
            }

            if (weight_sum == 0) {
                excluded_vertices.push_back(i);
                kmean[i] = 0;
            } else {
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
 * @brief Performs cross-validation for path-based kernel-weighted means on a graph using pre-computed shortest paths.
 *
 * This function estimates the prediction error of path-based kernel-weighted means through
 * cross-validation, utilizing pre-computed shortest paths for efficiency. For each round of
 * cross-validation, it randomly selects test vertices, computes path-based means excluding
 * these vertices, and calculates the absolute deviation between predicted and actual values.
 *
 * @param graph A vector of vectors representing the graph structure where graph[i] contains
 *              indices of vertices adjacent to vertex i.
 * @param edge_lengths A vector of vectors containing the lengths/weights of edges where
 *                    edge_lengths[i][j] is the length of edge from vertex i to its jth neighbor.
 * @param core_graph A vector of vectors representing the core graph structure used for finding
 *                   shortest paths between vertices.
 * @param shortest_paths A map containing pre-computed shortest paths between vertex pairs.
 *                      Keys are pairs of vertices (i,j) where i < j.
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
 * @note This version uses cached shortest paths for efficiency compared to computing them on-the-fly.
 *
 * @throws std::invalid_argument if n_CVs <= 0
 *
 * @pre The graph, edge_lengths, core_graph, and y vectors must have consistent sizes
 * @pre shortest_paths must contain all required vertex pairs
 * @pre The function assumes initialize_kernel and kernel_fn are properly set up
 *
 * @see graph_cached_kpath_mean_with_weights
 */
std::vector<double> graph_cached_kpath_mean_cv(
        const std::vector<std::vector<int>>& graph,
        const std::vector<std::vector<double>>& edge_lengths,
        const std::vector<std::vector<int>>& core_graph,
        const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
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

        // Estimating the conditional expectation of cv_y using cached path-based means
        auto res = graph_cached_kpath_mean_with_weights(
            graph,
            edge_lengths,
            core_graph,
            shortest_paths,
            weights,
            y,
            ikernel,
            dist_normalization_factor
        );
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

/**
 * @brief R interface for cross-validation of path-based kernel-weighted means using pre-computed shortest paths.
 *
 * This function serves as an interface between R and the C++ implementation of graph_cached_kpath_mean_cv.
 * It converts R objects to C++ data structures, performs the cross-validation, and returns the results
 * back to R.
 *
 * @param s_graph SEXP containing the graph structure (list of integer vectors)
 * @param s_edge_lengths SEXP containing edge lengths (list of numeric vectors)
 * @param s_core_graph SEXP containing the core graph structure (list of integer vectors)
 * @param s_shortest_paths SEXP containing pre-computed shortest paths (list with names as "i_j")
 * @param s_y SEXP containing vertex values (numeric vector)
 * @param s_ikernel SEXP containing kernel type (single integer)
 * @param s_dist_normalization_factor SEXP containing normalization factor (single numeric)
 * @param s_n_CVs SEXP containing number of CV iterations (single integer)
 * @param s_n_CV_folds SEXP containing number of CV folds (single integer)
 * @param s_epsilon SEXP containing epsilon value for numerical stability (single numeric)
 * @param s_seed SEXP containing random seed (single integer)
 *
 * @return SEXP containing a numeric vector of cross-validation errors for each vertex
 *
 * @note The function assumes all input SEXPs are of the correct type and properly structured
 * @note Memory management is handled using PROTECT/UNPROTECT for R compatibility
 *
 * @pre s_graph must be a list of integer vectors representing the graph adjacency list
 * @pre s_edge_lengths must be a list of numeric vectors matching s_graph structure
 * @pre s_core_graph must be a list of integer vectors representing the core graph
 * @pre s_shortest_paths must be a list with names in "i_j" format containing integer vectors
 * @pre s_y must be a numeric vector with length matching the number of vertices
 * @pre s_ikernel must be a single integer
 * @pre s_dist_normalization_factor must be a single numeric value
 * @pre s_n_CVs must be a positive integer
 * @pre s_n_CV_folds must be a positive integer
 * @pre s_epsilon must be a positive numeric value
 * @pre s_seed must be a non-negative integer
 *
 * @see graph_cached_kpath_mean_cv, Rgraph_to_vector, Rweights_to_vector, shortest_paths_Rlist_to_cpp_map
 */
SEXP S_graph_cached_kpath_mean_cv(SEXP s_graph,
                                  SEXP s_edge_lengths,
                                  SEXP s_core_graph,
                                  SEXP s_shortest_paths,
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
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths = shortest_paths_Rlist_to_cpp_map(s_shortest_paths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    std::vector<double> cv_errors =  graph_cached_kpath_mean_cv(graph,
                                                                edge_lengths,
                                                                core_graph,
                                                                shortest_paths,
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
