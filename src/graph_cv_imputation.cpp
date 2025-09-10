#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

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
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_diffusion_smoother.hpp"
#include "stats_utils.h"
#include "kernels.h"

std::unique_ptr<std::vector<std::vector<double>>> dist_to_weights(const std::vector<std::vector<int>>& graph,
                                                                  const std::vector<std::vector<double>>& d,
                                                                  int ikernel,
                                                                  double dist_normalization_factor);

extern "C" {
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
    double scale = 1.0;
    initialize_kernel(ikernel, scale);
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

        impute_binary_value_by_iterative_neighborhood_matching(graph, test_graph, test_set,
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

        impute_continuous_value_by_iterative_neighborhood_matching(graph, test_graph, test_set,
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
 * It does not perform extensive Rf_error checking on the inputs.
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
    PROTECT(Rtest_set = Rf_coerceVector(Rtest_set, INTSXP)); nprot++;
    PROTECT(Rgraph = Rf_coerceVector(Rgraph, VECSXP)); nprot++;
    PROTECT(Rd = Rf_coerceVector(Rd, VECSXP)); nprot++;
    PROTECT(Ry = Rf_coerceVector(Ry, REALSXP));nprot++;
    PROTECT(Ry_binary = Rf_coerceVector(Ry_binary, LGLSXP)); nprot++;
    PROTECT(Rimputation_method = Rf_coerceVector(Rimputation_method, INTSXP)); nprot++;
    PROTECT(Rmax_iterations = Rf_coerceVector(Rmax_iterations, INTSXP)); nprot++;
    PROTECT(Rconvergence_threshold = Rf_coerceVector(Rconvergence_threshold, REALSXP)); nprot++;
    PROTECT(Rapply_binary_threshold = Rf_coerceVector(Rapply_binary_threshold, LGLSXP)); nprot++;
    PROTECT(Rbinary_threshold = Rf_coerceVector(Rbinary_threshold, REALSXP)); nprot++;
    PROTECT(Rikernel = Rf_coerceVector(Rikernel, INTSXP)); nprot++;
    PROTECT(Rdist_normalization_factor = Rf_coerceVector(Rdist_normalization_factor, REALSXP)); nprot++;

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
    SEXP Rresult = PROTECT(Rf_allocVector(REALSXP, result->size())); nprot++;
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
    double scale = 1.0;
    initialize_kernel(ikernel, scale);
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
