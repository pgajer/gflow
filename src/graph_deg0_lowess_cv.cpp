#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <mutex>                    // For std::mutex
#include <execution>                // For std::execution::par_unseq
#include <atomic>                   // For std::atomic
#include <thread>                   // For std::thread::hardware_concurrenyc

#include "exec_policy.hpp"
#include "graph_deg0_lowess_cv.hpp" // For graph_deg0_lowess_cv_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

extern "C" {
    SEXP S_graph_deg0_lowess_cv(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_use_uniform_weights,
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose);
}

/**
 * @brief Performs degree 0 LOWESS with cross-validation for bandwidth selection
 *
 * @details This function implements a graph-based extension of LOWESS degree 0
 * (locally weighted average) with spatially-stratified cross-validation for
 * bandwidth selection. The algorithm:
 * 1. Creates a maximal packing of vertices to serve as fold seed points
 * 2. Assigns all vertices to the nearest seed point to form spatially coherent folds
 * 3. For each candidate bandwidth, performs cross-validation across the folds
 * 4. Selects the bandwidth with the lowest cross-validation Rf_error
 * 5. Fits the final model with the optimal bandwidth
 *
 * @param y Response values at each vertex in the graph
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
 * @return graph_deg0_lowess_cv_t Structure containing optimal model results and CV diagnostics
 */
graph_deg0_lowess_cv_t set_wgraph_t::graph_deg0_lowess_cv(
    const std::vector<double>& y,
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

    unsigned int available_threads = std::thread::hardware_concurrency();
    if (available_threads == 0) available_threads = 4; // Fallback if detection fails

    size_t n_vertices = adjacency_list.size();

    if (verbose) {
        Rprintf("Starting graph_deg0_lowess_cv() algorithm\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
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


    // Define minimum and maximum bandwidth based on graph diameter
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
    graph_deg0_lowess_cv_t result;
    result.bws = bw_grid;
    result.bw_errors.resize(bw_grid.size(), 0.0);
    if (with_bw_predictions) {
        result.bw_predictions.resize(bw_grid.size(),
                                     std::vector<double>(n_vertices, 0.0));
    }

    // 2. Define the prediction function that can handle custom weights
    auto predict_with_weights = [&](double bandwidth, const std::vector<double>& weights) 
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
                // double uniform_weight = 1.0 / nbhd_size;
                double uniform_weight = 1.0;
                std::fill(kernel_weights.begin(), kernel_weights.end(), uniform_weight);
            } else {
                // Distance normalization strategy:
                // 1. Find the maximum distance (max_dist) in the vertex neighborhood
                // 2. Scale max_dist by dist_normalization_factor (typically 1.1)
                // 3. Divide all distances by this scaled maximum
                // This approach ensures all normalized distances fall within [0, 1/dist_normalization_factor],
                // which is approximately [0, 0.91] when dist_normalization_factor = 1.1
                // This keeps all distances within the effective support of the kernel functions,
                // as most kernels in this implementation have support on [-1, 1]

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

    // 3. Perform cross-validation for each bandwidth
    if (verbose) {
        Rprintf("Starting cross-validation across %zu bandwidths...\n", bw_grid.size());
    }

    // Setup for multi-threaded processing
    std::mutex mutex;
    progress_tracker_t progress(actual_folds, "Processing folds", 1);
    std::atomic<size_t> completed_tasks(0);

    // Process each fold in parallel
    std::vector<size_t> fold_indices(actual_folds);
    std::iota(fold_indices.begin(), fold_indices.end(), 0);

    std::vector<size_t> fold_test_counts(actual_folds, 0);
    std::vector<double> fold_total_errors(actual_folds, 0.0);

    //std::for_each(std::execution::par_unseq, fold_indices.begin(), fold_indices.end(),
    std::for_each(GFLOW_EXEC_POLICY, fold_indices.begin(), fold_indices.end(),
                  [&](size_t fold) {
                      try {
                          if (verbose) {
                              std::lock_guard<std::mutex> lock(mutex);
                              Rprintf("Processing fold %zu/%zu\n", fold+1, actual_folds);
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
                                  predict_with_weights(bw_grid[bw_idx], fold_weights);

                              // Calculate Rf_error on test vertices
                              double fold_error = 0.0;
                              size_t test_count = 0;

                              for (size_t vertex : folds[fold]) {
                                  if (!excluded[vertex]) {  // Only count non-excluded vertices
                                      fold_error += std::pow(predictions[vertex] - y[vertex], 2);
                                      test_count++;
                                  }
                              }

                              // Store fold results
                              if (test_count > 0) {
                                  std::lock_guard<std::mutex> lock(mutex);
                                  result.bw_errors[bw_idx] += fold_error;
                                  fold_test_counts[fold] = test_count;
                                  fold_total_errors[fold] += fold_error;
                              }

                              // Store predictions if requested
                              if (with_bw_predictions) {
                                  std::lock_guard<std::mutex> lock(mutex);
                                  // Only update predictions for test vertices from this fold
                                  for (size_t vertex : folds[fold]) {
                                      result.bw_predictions[bw_idx][vertex] = predictions[vertex];
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
                          REPORT_ERROR("Error processing fold %zu: %s", fold, e.what());
                      }
                  });

    // Normalize errors by total test points
    size_t total_test_points = 0;
    for (auto count : fold_test_counts) {
        total_test_points += count;
    }

    if (total_test_points > 0) {
        for (auto& Rf_error : result.bw_errors) {
            Rf_error /= total_test_points;
        }
    }

    // 4. Find the optimal bandwidth
    auto min_it = std::min_element(result.bw_errors.begin(), result.bw_errors.end());
    result.opt_bw_idx = min_it - result.bw_errors.begin();
    result.opt_bw = bw_grid[result.opt_bw_idx];

    if (verbose) {
        Rprintf("Optimal bandwidth index: %zu\n", result.opt_bw_idx);
        Rprintf("Optimal bandwidth: %.4f\n", result.opt_bw);
        Rprintf("CV Rf_error: %.6f\n", result.bw_errors[result.opt_bw_idx]);
    }

    // 5. Final prediction with optimal bandwidth
    std::vector<double> uniform_weights(n_vertices, 1.0);  // Use all vertices
    auto [final_predictions, _] = predict_with_weights(result.opt_bw, uniform_weights);
    result.predictions = final_predictions;

    if (verbose) {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        Rprintf("Completed graph_deg0_lowess_cv in %.2f seconds\n", elapsed);
    }

    return result;
}

/**
 * @brief R interface for graph degree 0 LOWESS with CV
 */
SEXP S_graph_deg0_lowess_cv(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_use_uniform_weights,
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose) {

    // Convert input parameters
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Get scalar parameters
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];
    size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
    bool log_grid = LOGICAL(s_log_grid)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool use_uniform_weights = LOGICAL(s_use_uniform_weights)[0];
    size_t n_folds = (size_t)INTEGER(s_n_folds)[0];
    bool with_bw_predictions = LOGICAL(s_with_bw_predictions)[0];
    double precision = REAL(s_precision)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    // Create graph and run algorithm
    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);
    graph_deg0_lowess_cv_t result = graph.graph_deg0_lowess_cv(
        y,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        log_grid,
        kernel_type,
        dist_normalization_factor,
        use_uniform_weights,
        n_folds,
        with_bw_predictions,
        precision,
        verbose
    );

    // Create the return list
    const char* names[] = {
        "predictions",
        "bw_predictions",
        "bw_errors",
        "bws",
        "opt_bw",
        "opt_bw_idx",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Track protection count
    int protect_count = 0;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    protect_count++;

    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    protect_count++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);

    // Helper function to convert vector to SEXP
    auto vec_to_sexp = [&protect_count](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        protect_count++;
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Set predictions
    SET_VECTOR_ELT(r_result, 0, vec_to_sexp(result.predictions));

    // Set bw_predictions (if requested)
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, result.bw_predictions.size()));
        protect_count++;
        
        for (size_t i = 0; i < result.bw_predictions.size(); i++) {
            SET_VECTOR_ELT(r_bw_predictions, i, vec_to_sexp(result.bw_predictions[i]));
        }
        
        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
    } else {
        SET_VECTOR_ELT(r_result, 1, Rf_allocVector(VECSXP, 0)); // Empty list
    }

    // Set bw_errors
    SET_VECTOR_ELT(r_result, 2, vec_to_sexp(result.bw_errors));

    // Set bws
    SET_VECTOR_ELT(r_result, 3, vec_to_sexp(result.bws));

    // Set opt_bw
    SEXP r_opt_bw = PROTECT(Rf_allocVector(REALSXP, 1));
    protect_count++;
    REAL(r_opt_bw)[0] = result.opt_bw;
    SET_VECTOR_ELT(r_result, 4, r_opt_bw);

    // Set opt_bw_idx
    SEXP r_opt_bw_idx = PROTECT(Rf_allocVector(INTSXP, 1));
    protect_count++;
    INTEGER(r_opt_bw_idx)[0] = result.opt_bw_idx + 1; // Convert to 1-based for R
    SET_VECTOR_ELT(r_result, 5, r_opt_bw_idx);

    UNPROTECT(protect_count);
    return r_result;
}
