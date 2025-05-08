#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <thread>                   // For std::thread::hardware_concurrency
#include <mutex>                    // For std::mutex

#include "graph_deg0_lowess_cv_mat.hpp" // For graph_deg0_lowess_cv_mat_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

/**
 * @brief Matrix version of degree 0 LOWESS with cross-validation for bandwidth selection
 *
 * @details This function implements a matrix version of graph-based extension of LOWESS degree 0
 * (locally weighted average) with spatially-stratified cross-validation for
 * bandwidth selection. The algorithm:
 * 1. Creates a maximal packing of vertices to serve as fold seed points
 * 2. Assigns all vertices to the nearest seed point to form spatially coherent folds
 * 3. For each candidate bandwidth, performs cross-validation across the folds for each response variable
 * 4. Selects the bandwidth with the lowest cross-validation error for each response variable
 * 5. Fits the final model with the optimal bandwidth for each response variable
 *
 * @param Y Matrix of response values where Y[j][i] is the j-th response variable at the i-th vertex
 * @param min_bw_factor Minimum bandwidth as a factor of graph diameter
 * @param max_bw_factor Maximum bandwidth as a factor of graph diameter
 * @param n_bws Number of bandwidths to test
 * @param log_grid Whether to use logarithmic spacing for bandwidth grid
 * @param kernel_type Type of kernel function for weighting
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 * @param n_folds Number of cross-validation folds
 * @param with_bw_predictions Whether to compute and store predictions for all bandwidths
 * @param precision Precision for bandwidth grid computation
 * @param verbose Whether to print progress information
 *
 * @return graph_deg0_lowess_cv_mat_t Structure containing optimal model results and CV diagnostics for each response variable
 */
graph_deg0_lowess_cv_mat_t set_wgraph_t::graph_deg0_lowess_cv_mat(
    const std::vector<std::vector<double>>& Y,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool log_grid,
    size_t kernel_type,
    double dist_normalization_factor,
    bool use_uniform_weights,
    size_t n_folds,
    bool with_bw_predictions,
    double precision,
    bool verbose) {

    auto start_time = std::chrono::steady_clock::now();

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    // Validate input
    if (Y.empty()) {
        REPORT_ERROR("Input matrix Y is empty");
    }
    
    size_t n_response_vars = Y.size();
    if (n_response_vars == 0) {
        REPORT_ERROR("No response variables provided");
    }
    
    size_t n_vertices = Y[0].size();
    for (size_t j = 1; j < n_response_vars; ++j) {
        if (Y[j].size() != n_vertices) {
            REPORT_ERROR("Inconsistent number of vertices across response variables");
        }
    }
    
    if (n_vertices != adjacency_list.size()) {
        REPORT_ERROR("Number of vertices in Y (%zu) does not match graph size (%zu)",
                   n_vertices, adjacency_list.size());
    }

    // Set number of threads for parallel computation
    unsigned int available_threads = std::thread::hardware_concurrency();
    if (available_threads == 0) available_threads = 4; // Fallback if detection fails

    if (verbose) {
        Rprintf("Starting graph_deg0_lowess_cv_mat() algorithm\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Number of response variables: %zu\n", n_response_vars);
        Rprintf("Number of CV folds: %zu\n", n_folds);
        Rprintf("Using %u threads for parallel operations\n", available_threads);
        Rprintf("Using uniform weights: %s\n", use_uniform_weights ? "TRUE" : "FALSE");
    }

    // 1. Create spatially stratified CV folds
    if (verbose) {
        Rprintf("Creating spatially stratified CV folds...\n");
    }

    // Determine appropriate seed density based on graph size
    double cv_seed_factor;
    if (n_vertices < 500) {
        // For small graphs: 15-20% of vertices as seeds
        cv_seed_factor = 0.2;  // Using upper end of range
    } else if (n_vertices < 5000) {
        // For medium graphs: 5-10% of vertices as seeds
        cv_seed_factor = 0.1;  // Using upper end of range
    } else {
        // For large graphs: 1-5% of vertices as seeds
        cv_seed_factor = 0.05;  // Using upper end of range
    }

    // Calculate target number of seeds based on graph size and desired number of folds
    size_t target_seeds = std::max(size_t(n_vertices * cv_seed_factor), n_folds * 3);
    target_seeds = std::min(target_seeds, n_vertices / 10);  // Cap to avoid too many seeds

    if (verbose) {
        Rprintf("Creating maximal packing with target of %zu seeds for %zu folds\n",
                target_seeds, n_folds);
    }

    // Select seed points using maximal packing with specified number of seeds
    size_t max_iterations = 100;  // Reasonable default
    std::vector<size_t> fold_seeds = create_maximal_packing(
        target_seeds,
        max_iterations,
        precision);

    // Ensure we have enough seeds
    size_t actual_seeds = fold_seeds.size();
    if (verbose) {
        Rprintf("Created maximal packing with %zu seed points\n", actual_seeds);
    }

    // If we have too few seeds, we can't do the requested number of folds
    size_t actual_folds = std::min(n_folds, actual_seeds);
    if (actual_folds < n_folds && verbose) {
        Rprintf("Warning: Could only create %zu folds (requested %zu)\n",
                actual_folds, n_folds);
    }

    // Create CV folds by grouping seeds if necessary
    std::vector<std::vector<size_t>> seed_groups;
    if (actual_seeds <= n_folds) {
        // Case 1: Not enough seeds, each seed becomes its own fold
        seed_groups.resize(actual_seeds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i].push_back(fold_seeds[i]);
        }
    } else {
        // Case 2: More seeds than folds, distribute seeds among folds
        seed_groups.resize(n_folds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i % n_folds].push_back(fold_seeds[i]);
        }
    }

    // Assign vertices to folds based on closest seed group
    std::vector<std::vector<size_t>> folds(actual_folds);
    std::vector<size_t> fold_assignment(n_vertices, SIZE_MAX);

    // For each vertex, find the closest seed group
    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
        double min_distance = std::numeric_limits<double>::infinity();
        size_t closest_fold = 0;

        // Check distance to each fold (via its seed groups)
        for (size_t fold_idx = 0; fold_idx < seed_groups.size(); ++fold_idx) {
            // Find minimum distance to any seed in this group
            for (size_t seed : seed_groups[fold_idx]) {
                double distance = compute_shortest_path_distance(vertex, seed);
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_fold = fold_idx;
                }
            }
        }

        // Assign vertex to closest fold
        folds[closest_fold].push_back(vertex);
        fold_assignment[vertex] = closest_fold;
    }

    if (verbose) {
        Rprintf("Created %zu folds with sizes: ", actual_folds);
        for (size_t i = 0; i < actual_folds; ++i) {
            Rprintf("%zu ", folds[i].size());
        }
        Rprintf("\n");
    }

    // 2. Define minimum and maximum bandwidth based on graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
    }

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    if (verbose) {
        Rprintf("Graph diameter: %f\n", graph_diameter);
        Rprintf("Bandwidth range: [%.4f, %.4f]\n", min_bw, max_bw);
    }

    // Generate bandwidth grid
    std::vector<double> bw_grid = get_candidate_bws(
        min_bw, max_bw, n_bws, log_grid, precision);

    if (verbose) {
        Rprintf("Testing %zu candidate bandwidths\n", bw_grid.size());
    }

    // Initialize result structure
    graph_deg0_lowess_cv_mat_t result;
    result.bws = bw_grid;
    
    // Initialize per-response variable structures
    result.bw_errors.resize(n_response_vars);
    result.opt_bws.resize(n_response_vars);
    result.opt_bw_idxs.resize(n_response_vars);
    result.predictions.resize(n_response_vars, std::vector<double>(n_vertices, 0.0));
    
    for (size_t j = 0; j < n_response_vars; ++j) {
        result.bw_errors[j].resize(bw_grid.size(), 0.0);
    }
    
    if (with_bw_predictions) {
        result.bw_predictions.resize(n_response_vars);
        for (size_t j = 0; j < n_response_vars; ++j) {
            result.bw_predictions[j].resize(bw_grid.size(), 
                                          std::vector<double>(n_vertices, 0.0));
        }
    }

    // 3. Define the prediction function that can handle custom weights
    auto predict_with_weights = [&](double bandwidth, const std::vector<double>& weights, const std::vector<double>& y) 
        -> std::pair<std::vector<double>, std::vector<bool>> {
        
        std::vector<double> predictions(n_vertices);
        std::vector<bool> excluded(n_vertices, false);

        // Process each vertex
        for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
            // Find vertices within bandwidth radius
            auto vertex_map = find_vertices_within_radius(vertex, bandwidth);
            size_t nbhd_size = vertex_map.size();

            // Prepare arrays for kernel calculation
            std::vector<double> normalized_dists(nbhd_size);
            std::vector<size_t> neighbors(nbhd_size);
            std::vector<double> neighbor_weights(nbhd_size);
            
            size_t counter = 0;
            for (const auto& [neighbor, distance] : vertex_map) {
                normalized_dists[counter] = distance;
                neighbors[counter] = neighbor;
                neighbor_weights[counter] = weights[neighbor];
                counter++;
            }

            // Compute kernel weights
            std::vector<double> kernel_weights(nbhd_size);

            // Check if using uniform weights
            if (use_uniform_weights) {
                // Set all weights to 1/nbhd_size for uniform weighting
                //double uniform_weight = 1.0 / nbhd_size;
                double uniform_weight = 1.0;
                std::fill(kernel_weights.begin(), kernel_weights.end(), uniform_weight);
            } else {

                // Distance normalization strategy
                double max_dist = 0.0;
                for (size_t k = 0; k < nbhd_size; ++k) {
                    max_dist = std::max(max_dist, normalized_dists[k]);
                }
                if (max_dist == 0) max_dist = 1;
                max_dist *= dist_normalization_factor;

                for (size_t k = 0; k < nbhd_size; ++k) {
                    normalized_dists[k] /= max_dist;
                }

                // Use regular kernel-based weighting
                kernel_fn(normalized_dists.data(), nbhd_size, kernel_weights.data());
            }

            // Apply user-specified weights (for CV)
            for (size_t k = 0; k < nbhd_size; ++k) {
                kernel_weights[k] *= neighbor_weights[k];
            }

            // Compute weighted average
            double weighted_sum = 0.0;
            double weight_sum = 0.0;

            for (size_t k = 0; k < nbhd_size; ++k) {
                weighted_sum += kernel_weights[k] * y[neighbors[k]];
                weight_sum += kernel_weights[k];
            }

            // Check if we need to exclude this vertex
            if (weight_sum <= 0) {
                excluded[vertex] = true;
                predictions[vertex] = std::numeric_limits<double>::quiet_NaN();
            } else {
                predictions[vertex] = weighted_sum / weight_sum;
            }
        }

        return {predictions, excluded};
    };

    // 4. Perform cross-validation for each bandwidth and each response variable
    if (verbose) {
        Rprintf("Starting cross-validation across %zu bandwidths and %zu response variables...\n", 
                bw_grid.size(), n_response_vars);
    }

    // Setup for multi-threaded processing
    std::mutex mutex;
    progress_tracker_t progress(actual_folds * n_response_vars, "Processing folds", 1);
    std::atomic<size_t> completed_tasks(0);

    // Process each response variable and fold
    for (size_t j = 0; j < n_response_vars; ++j) {
        std::vector<size_t> fold_test_counts(actual_folds, 0);
        std::vector<double> fold_total_errors(actual_folds, 0.0);

        // Process each fold in parallel
        std::vector<size_t> fold_indices(actual_folds);
        std::iota(fold_indices.begin(), fold_indices.end(), 0);

        std::for_each(std::execution::par_unseq, fold_indices.begin(), fold_indices.end(),
                     [&](size_t fold) {
            try {
                if (verbose) {
                    std::lock_guard<std::mutex> lock(mutex);
                    Rprintf("Processing fold %zu/%zu for response variable %zu/%zu\n", 
                            fold+1, actual_folds, j+1, n_response_vars);
                }

                // Generate weights for this fold (1.0 for training, 0.0 for testing)
                std::vector<double> fold_weights(n_vertices, 1.0);
                for (size_t vertex : folds[fold]) {
                    fold_weights[vertex] = 0.0;  // Test vertex
                }

                // Process each bandwidth
                for (size_t bw_idx = 0; bw_idx < bw_grid.size(); ++bw_idx) {
                    // Get predictions with this bandwidth
                    auto [predictions, excluded] = 
                        predict_with_weights(bw_grid[bw_idx], fold_weights, Y[j]);

                    // Calculate error on test vertices
                    double fold_error = 0.0;
                    size_t test_count = 0;

                    for (size_t vertex : folds[fold]) {
                        if (!excluded[vertex]) {  // Only count non-excluded vertices
                            fold_error += std::pow(predictions[vertex] - Y[j][vertex], 2);
                            test_count++;
                        }
                    }

                    // Store fold results (thread-safe)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        result.bw_errors[j][bw_idx] += fold_error;
                        fold_test_counts[fold] = test_count;
                        fold_total_errors[fold] += fold_error;
                    }

                    // Store predictions if requested
                    if (with_bw_predictions) {
                        std::lock_guard<std::mutex> lock(mutex);
                        // Only update predictions for test vertices from this fold
                        for (size_t vertex : folds[fold]) {
                            result.bw_predictions[j][bw_idx][vertex] = predictions[vertex];
                        }
                    }
                }

                // Update progress
                if (verbose) {
                    size_t tasks_completed = completed_tasks.fetch_add(1, std::memory_order_relaxed) + 1;
                    std::lock_guard<std::mutex> lock(mutex);
                    progress.update(tasks_completed);
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(mutex);
                REPORT_ERROR("Error processing fold %zu for response variable %zu: %s", 
                             fold, j, e.what());
            }
        });

        // Normalize errors by total test points
        size_t total_test_points = 0;
        for (auto count : fold_test_counts) {
            total_test_points += count;
        }

        if (total_test_points > 0) {
            for (auto& error : result.bw_errors[j]) {
                error /= total_test_points;
            }
        }

        // 4. Find the optimal bandwidth for this response variable
        auto min_it = std::min_element(result.bw_errors[j].begin(), result.bw_errors[j].end());
        result.opt_bw_idxs[j] = min_it - result.bw_errors[j].begin();
        result.opt_bws[j] = bw_grid[result.opt_bw_idxs[j]];

        if (verbose) {
            Rprintf("Response variable %zu - Optimal bandwidth index: %zu\n", j, result.opt_bw_idxs[j]);
            Rprintf("Response variable %zu - Optimal bandwidth: %.4f\n", j, result.opt_bws[j]);
            Rprintf("Response variable %zu - CV error: %.6f\n", j, result.bw_errors[j][result.opt_bw_idxs[j]]);
        }
    }

    if (verbose) {
        progress.finish();
    }

    // 5. Final prediction with optimal bandwidth for each response variable
    if (verbose) {
        Rprintf("Computing final predictions using optimal bandwidths...\n");
    }

    // Process each response variable in parallel
    std::vector<size_t> response_indices(n_response_vars);
    std::iota(response_indices.begin(), response_indices.end(), 0);

    std::for_each(std::execution::par_unseq, response_indices.begin(), response_indices.end(),
                 [&](size_t j) {
        try {
            std::vector<double> uniform_weights(n_vertices, 1.0);  // Use all vertices
            auto [final_predictions, _] = predict_with_weights(result.opt_bws[j], uniform_weights, Y[j]);
            result.predictions[j] = final_predictions;
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(mutex);
            REPORT_ERROR("Error computing final predictions for response variable %zu: %s", j, e.what());
        }
    });

    if (verbose) {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        Rprintf("Completed graph_deg0_lowess_cv_mat in %.2f seconds\n", elapsed);
    }

    return result;
}
