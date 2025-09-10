/**
 * @file wasserstein_perm_test.cpp
 * @brief Implementation of Wasserstein permutation test for vertex-wise statistical analysis
 * @author Pawel Gajer
 * @date Feb 16, 2025
 *
 * @details
 * This file implements a permutation test using Wasserstein distance metrics to analyze
 * differences between backbone and null model predictions at the vertex level. The test
 * performs bootstrap sampling to generate a null distribution and computes effect sizes
 * and p-values for each vertex.
 *
 * The implementation supports parallel execution through OpenMP and uses thread-local
 * storage to ensure thread safety during computations.
 *
 * @note Requires OpenMP support for parallel execution
 *
 * @dependencies
 * - Standard Template Library (STL)
 * - OpenMP
 * - Custom Wasserstein distance calculation (C_wasserstein_distance_1D)
 *
 * @see vertex_wasserstein_perm_test_results_t
 * @see C_wasserstein_distance_1D
 */

#include <vector>          // for std::vector
#include <random>          // for std::random_device, std::mt19937, std::uniform_int_distribution
#include <utility>         // for std::pair
#include <algorithm>       // for std::for_each, std::count_if
#include <numeric>         // for std::iota
#include <execution>       // for std::execution::seq, std::execution::par_unseq
#include <thread>          // for std::thread::hardware_concurrency()
#include <cstddef>         // for size_t

#include "wasserstein_perm_test.hpp"  // Contains vertex_wasserstein_perm_test_results_t definition and C_wasserstein_distance_1D
#include "exec_policy.hpp"
#include "omp_compat.h"

/**
 * @brief Performs a vertex-wise Wasserstein permutation test
 *
 * @param bb_predictions Vector of vectors containing backbone model predictions for each vertex
 *                      Outer vector: samples, Inner vector: vertices
 * @param null_predictions Vector of vectors containing null model predictions for each vertex
 *                        Outer vector: samples, Inner vector: vertices
 * @param n_bootstraps Number of bootstrap samples to generate for null distribution [default: 1000]
 * @param alpha Significance level for hypothesis testing [default: 0.05]
 *
 * @return vertex_wasserstein_perm_test_results_t Structure containing:
 *         - p_values: Vector of p-values for each vertex
 *         - effect_sizes: Vector of Wasserstein distances (effect sizes) for each vertex
 *         - significant_vertices: Boolean vector indicating significance at alpha level
 *         - null_distances: Vector of Wasserstein distances from bootstrap sampling
 *
 * @throws No explicit exceptions, but STL containers may throw std::bad_alloc
 *
 * @note The function uses thread-local storage and OpenMP for parallel computation
 * @note The number of tests performed is min(bb_predictions.size(), null_predictions.size())
 *
 * @Rf_warning Input vectors must be non-empty and have consistent dimensions
 * @Rf_warning Requires sufficient memory to store bootstrap samples and results
 */
vertex_wasserstein_perm_test_results_t vertex_wasserstein_perm_test(
    const std::vector<std::vector<double>>& bb_predictions,
    const std::vector<std::vector<double>>& null_predictions,
    size_t n_bootstraps,
    double alpha
    ) {
    const size_t n_vertices = bb_predictions[0].size();
    const size_t n_bb = bb_predictions.size();
    const size_t n_perms = null_predictions.size();
    const size_t n_tests = std::min(n_bb, n_perms);

    vertex_wasserstein_perm_test_results_t results;
    results.p_values.resize(n_vertices);
    results.effect_sizes.resize(n_vertices);
    results.significant_vertices.resize(n_vertices);
    results.null_distances.resize(n_bootstraps);  // Changed from reserve to resize

    // Extract values for each vertex from null predictions
    std::vector<std::vector<double>> null_values(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        null_values[i].resize(n_tests);
        for (size_t j = 0; j < n_tests; ++j) {
            null_values[i][j] = null_predictions[j][i];
        }
    }

    // Generate bootstrap pairs for null distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::pair<size_t, size_t>> pairs;
    pairs.reserve(n_bootstraps);

    // Generate unique pairs for bootstrap sampling
    for (size_t k = 0; k < n_bootstraps; ++k) {
        size_t i = 0, j = 0;
        do {
            i = std::uniform_int_distribution<size_t>(0, n_vertices - 1)(gen);
            j = std::uniform_int_distribution<size_t>(0, n_vertices - 1)(gen);
        } while (i == j);
        pairs.emplace_back(i, j);
    }

    // Compute null distribution of Wasserstein distances
    std::vector<size_t> bootstrap_indices(n_bootstraps);
    std::iota(bootstrap_indices.begin(), bootstrap_indices.end(), 0);

    gflow::for_each(gflow::seq,
                  bootstrap_indices.begin(),
                  bootstrap_indices.end(),
                  [&](size_t k) {  // Changed from int to size_t to match index type
                      const auto& [i, j] = pairs[k];
                      int n_test_int = static_cast<int>(n_tests);
                      C_wasserstein_distance_1D(
                          null_values[i].data(),
                          null_values[j].data(),
                          &n_test_int,
                          &results.null_distances[k]  // Fixed: use k instead of undefined variable
                          );
                  });

    // Extract bb values and compute effect sizes for each vertex
    std::vector<size_t> vertex_indices(n_vertices);
    std::iota(vertex_indices.begin(), vertex_indices.end(), 0);

    // Create thread-local storage for vertex_bb_values
    const size_t n_threads = std::thread::hardware_concurrency();
    std::vector<std::vector<double>> thread_local_bb_values(n_threads, std::vector<double>(n_tests));

    gflow::for_each(gflow::seq,
                  vertex_indices.begin(),
                  vertex_indices.end(),
                  [&](size_t i) {  // Changed from int to size_t
                      // Get thread-local storage
                      const int thread_id = omp_get_thread_num() % n_threads;
                      auto& vertex_bb_values = thread_local_bb_values[thread_id];

                      // Get bootstrap distribution for this vertex
                      for (size_t j = 0; j < n_tests; ++j) {
                          vertex_bb_values[j] = bb_predictions[j][i];
                      }

                      // Compute Wasserstein distance between bb and null distributions
                      int n_test_int = static_cast<int>(n_tests);
                      C_wasserstein_distance_1D(
                          vertex_bb_values.data(),
                          null_values[i].data(),
                          &n_test_int,
                          &results.effect_sizes[i]
                          );

                      // Compute p-value
                      size_t more_extreme = std::count_if(
                          results.null_distances.begin(),
                          results.null_distances.end(),
                          [es = results.effect_sizes[i]](double d) { return d >= es; }
                          );

                      results.p_values[i] = static_cast<double>(more_extreme) / n_bootstraps;
                      results.significant_vertices[i] = results.p_values[i] < alpha;
                  });

    return results;
}
