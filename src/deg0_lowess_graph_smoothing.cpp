#include <Rinternals.h>  // For Rprintf
#include <R_ext/Print.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <chrono>
#include <memory>

#include "set_wgraph.hpp"
#include "iknn_graphs_r.h"
#include "iknn_graphs.hpp"
#include "graph_deg0_lowess_cv.hpp"
#include "progress_utils.hpp"
#include "deg0_lowess_graph_smoothing.hpp"

deg0_lowess_graph_smoothing_t deg0_lowess_graph_smoothing(
    set_wgraph_t& initial_graph,
    const std::vector<std::vector<double>>& X,
    // Iteration control parameters
    size_t max_iterations,
    double convergence_threshold,
    ConvergenceCriteriaType convergence_type,
    // Graph construction parameters
    size_t k,
    double pruning_thld,
    // Bandwidth parameters
    size_t n_bws,
    bool log_grid,
    double min_bw_factor,
    double max_bw_factor,
    // Kernel parameters
    double dist_normalization_factor,
    size_t kernel_type,
    // Cross-validation parameters
    size_t n_folds,
    bool use_uniform_weights,
    // Other parameters
	double outlier_thld,
    bool with_bw_predictions,
    // Boosting parameters
    size_t switch_to_residuals_after,
    bool verbose) {

    // Start timing for the entire process
    auto total_start_time = std::chrono::steady_clock::now();

    // Validation
    if (X.empty()) {
        REPORT_ERROR("Empty data matrix provided");
    }

    if (max_iterations == 0) {
        REPORT_ERROR("Maximum iterations must be greater than 0");
    }

    // Initialize result structure
    deg0_lowess_graph_smoothing_t result;
    result.iterations_performed = 0;
    result.used_boosting = (switch_to_residuals_after < max_iterations);

    // Add initial graph and data matrix
    result.smoothed_graphs.push_back(initial_graph);
    result.smoothed_X.push_back(X);

    // Initialize vertex mapping (at first, each vertex maps to itself)
    std::vector<size_t> initial_mapping(initial_graph.num_vertices());
    for (size_t i = 0; i < initial_mapping.size(); ++i) {
        initial_mapping[i] = i;
    }
    result.vertex_mappings.push_back(initial_mapping);

    // Create a copy of the original data to use for calculating residuals
    const std::vector<std::vector<double>> original_X = X;

    // Initialize a cumulative smoothed version for boosting mode
    std::vector<std::vector<double>> cumulative_smoothed_X;

    // Progress tracking
    std::unique_ptr<progress_tracker_t> progress;
    if (verbose) {
        progress = std::make_unique<progress_tracker_t>(max_iterations, "Degree 0 LOWESS Graph Smoothing", 1);
        Rprintf("Starting degree 0 LOWESS graph smoothing with max %zu iterations\n", max_iterations);
        Rprintf("Convergence threshold: %.6f, Type: %d\n", convergence_threshold,
                static_cast<int>(convergence_type));
        if (switch_to_residuals_after < max_iterations) {
            Rprintf("Using boosting mode: switching to residual smoothing after %zu iterations\n",
                    switch_to_residuals_after);
        }
    }

    // Main iteration loop
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        auto iteration_start_time = std::chrono::steady_clock::now();

        if (verbose) {
            Rprintf("\n--- Iteration %zu/%zu ---\n", iter + 1, max_iterations);
        }

        // 1. Get current graph and data
        set_wgraph_t& current_graph = result.smoothed_graphs.back();

        // Determine the input data for this iteration
        std::vector<std::vector<double>> current_input_X;

        // Check if we should use direct smoothing or residual smoothing
        bool use_residual_smoothing = (iter >= switch_to_residuals_after);

        if (use_residual_smoothing) {
            if (verbose) {
                Rprintf("Using residual smoothing (boosting mode)\n");
            }

            // In first residual iteration, initialize cumulative smoothed X
            if (iter == switch_to_residuals_after) {
                cumulative_smoothed_X = result.smoothed_X.back();
            }

            // Calculate residuals between original data and current cumulative estimate
            current_input_X.resize(original_X.size(), std::vector<double>(original_X[0].size()));
            for (size_t j = 0; j < original_X.size(); ++j) {
                for (size_t i = 0; i < original_X[j].size(); ++i) {
                    current_input_X[j][i] = original_X[j][i] - cumulative_smoothed_X[j][i];
                }
            }

            if (verbose) {
                Rprintf("Calculated residuals for boosting\n");
            }
        } else {
            // Use direct smoothing - use the latest smoothed data
            current_input_X = result.smoothed_X.back();

            if (verbose && iter > 0) {
                Rprintf("Using direct smoothing\n");
            }
        }

        // 2. Apply degree 0 LOWESS to estimate conditional expectations
        if (verbose) {
            Rprintf("Computing degree 0 LOWESS smoothing with cross-validation...\n");
        }

        auto lowess_start_time = std::chrono::steady_clock::now();

        // Apply degree 0 LOWESS to each feature
        std::vector<std::vector<double>> smoothed_features;
        smoothed_features.resize(current_input_X.size());

        // Process each feature sequentially
        for (size_t j = 0; j < current_input_X.size(); ++j) {
            if (verbose) {
                Rprintf("Processing feature %zu/%zu\n", j + 1, current_input_X.size());
            }

            // Apply degree 0 LOWESS with cross-validation to this feature
            graph_deg0_lowess_cv_t lowess_result = current_graph.graph_deg0_lowess_cv(
                current_input_X[j],
                min_bw_factor,
                max_bw_factor,
                n_bws,
                log_grid,
                kernel_type,
                dist_normalization_factor,
                use_uniform_weights,
                n_folds,
                with_bw_predictions,
                1e-6,  // precision
                false // verbose
            );

            // Store the predictions for this feature
            smoothed_features[j] = lowess_result.predictions;

            if (verbose) {
                Rprintf("Feature %zu: Optimal bandwidth = %.6f (index %zu)\n",
                       j + 1, lowess_result.opt_bw, lowess_result.opt_bw_idx);
            }
        }

        if (verbose) {
            elapsed_time(lowess_start_time, "Degree 0 LOWESS complete", true);
        }

        // 3. Create new data matrix from LOWESS predictions
        std::vector<std::vector<double>> new_smoothed_component = smoothed_features;

        // Calculate the new smoothed X based on mode
        std::vector<std::vector<double>> new_X;

        if (use_residual_smoothing) {
            // In boosting mode, add the smoothed residuals to the cumulative result
            new_X.resize(original_X.size(), std::vector<double>(original_X[0].size()));
            for (size_t j = 0; j < original_X.size(); ++j) {
                for (size_t i = 0; i < original_X[j].size(); ++i) {
                    // Add smoothed residuals to the current cumulative estimate
                    new_X[j][i] = cumulative_smoothed_X[j][i] + new_smoothed_component[j][i];
                }
            }

            // Update the cumulative estimate for next iteration
            cumulative_smoothed_X = new_X;

            if (verbose) {
                Rprintf("Updated cumulative smoothed data with smoothed residuals\n");
            }
        } else {
            // In direct smoothing mode, just use the LOWESS output
            new_X = new_smoothed_component;
        }

        // 4. Create new graph from smoothed data
        set_wgraph_t new_graph;
        size_t ncc = 0;
        if (iter < max_iterations - 1) {  // Skip for last iteration to save time
            new_graph = create_iknn_graph_from_matrix(new_X, k, pruning_thld, verbose);

            // Find the number of connected components of the new graph
            ncc = new_graph.count_connected_components();
            if (ncc > 1) {

                if (verbose) {
                    Rprintf("Found %zu connected components in iteration %zu.\n", ncc, iter + 1);
                }

                // Identify the components
                std::vector<std::vector<size_t>> components = new_graph.get_connected_components();

                // Find the largest component (assumed to be the main cluster)
                size_t largest_component_idx = 0;
                size_t largest_component_size = 0;

                for (size_t i = 0; i < components.size(); ++i) {
                    if (components[i].size() > largest_component_size) {
                        largest_component_size = components[i].size();
                        largest_component_idx = i;
                    }
                }

                // Add outlier vertices to the outlier list (using the current mapping to get original indices)
                size_t n_outliers = 0;
                std::vector<size_t>& current_mapping = result.vertex_mappings.back();
                for (size_t i = 0; i < components.size(); ++i) {
                    if (i != largest_component_idx) {  // Not the largest component
                        for (size_t vertex : components[i]) {
                            // Map back to original index
                            result.outlier_indices.push_back(current_mapping[vertex]);
                            n_outliers++;
                        }
                    }
                }

                // Check if the number of outliers does not exceed the outlier_thld limit
                size_t n_vertices_of_new_graph = new_graph.num_vertices();
                double prop_of_outliers = (double)(n_outliers) / (double)(n_vertices_of_new_graph);

                if (prop_of_outliers > outlier_thld) {
                    Rprintf("The proportion of outliers: %.3f (%zu out of %zu) exceeds the outlier threshold limit %.3f\n. Aborting smoothing iterations\n",
                            prop_of_outliers, n_outliers, n_vertices_of_new_graph, outlier_thld);
                    break;
                }

                // Create a new mapping for the remaining vertices
                std::vector<size_t> new_mapping;
                new_mapping.reserve(components[largest_component_idx].size());

                // Create a mapping from old indices to new indices
                std::unordered_map<size_t, size_t> old_to_new_idx;
                for (size_t i = 0; i < components[largest_component_idx].size(); ++i) {
                    size_t old_idx = components[largest_component_idx][i];
                    old_to_new_idx[old_idx] = i;
                    // Map to original indices
                    new_mapping.push_back(current_mapping[old_idx]);
                }

                // Create a new subgraph with only the largest component
                set_wgraph_t subgraph = new_graph.create_subgraph(components[largest_component_idx]);

                // Create a new data matrix with only the vertices in the largest component
                std::vector<std::vector<double>> subgraph_X;
                subgraph_X.resize(new_X.size(), std::vector<double>(components[largest_component_idx].size()));

                for (size_t j = 0; j < new_X.size(); ++j) {
                    for (size_t i = 0; i < components[largest_component_idx].size(); ++i) {
                        size_t old_idx = components[largest_component_idx][i];
                        subgraph_X[j][i] = new_X[j][old_idx];
                    }
                }

                // Replace the current graph and data with the subgraph
                new_graph = std::move(subgraph);
                new_X = std::move(subgraph_X);

                // Store the new mapping
                result.vertex_mappings.push_back(new_mapping);

                if (verbose) {
                    Rprintf("Continuing with largest component (%zu vertices, %.1f%% of total)\n",
                            largest_component_size,
                            100.0 * largest_component_size / current_mapping.size());
                    Rprintf("Identified %zu outlier vertices in this iteration\n",
                            components.size() - 1);
                }

            } else {
                // No new components found, so the mapping stays the same
                result.vertex_mappings.push_back(result.vertex_mappings.back());
            }

            // Find diameter endpoints
            auto [end1, diam] = new_graph.get_vertex_eccentricity(0);  // Start from vertex 0
            auto [end2, diameter] = new_graph.get_vertex_eccentricity(end1);
            new_graph.graph_diameter = diameter;
        }

        // 5. Compute convergence metric against previous iteration's result
        const std::vector<std::vector<double>>& previous_X = result.smoothed_X.back();
        const std::vector<size_t>& previous_mapping = result.vertex_mappings.back();
        const std::vector<size_t>& current_mapping = ncc > 1 ?
            result.vertex_mappings.back() : result.vertex_mappings[result.vertex_mappings.size() - 2];

        double diff;
        // Use the appropriate difference computation based on matrix sizes
        if (previous_X[0].size() != new_X[0].size()) {
            if (verbose) {
                Rprintf("Using mapping-aware difference computation for differently sized matrices\n");
            }
            diff = compute_difference_with_mapping(previous_X, new_X, convergence_type,
                                                   current_mapping, previous_mapping);
        } else {
            // Use the original compute_difference function when matrices are the same size
            diff = compute_difference(previous_X, new_X, convergence_type);
        }

        result.convergence_metrics.push_back(diff);

        if (verbose) {
            Rprintf("Iteration %zu convergence metric: %.6f\n", iter + 1, diff);
        }

        // 6. Store results
        result.smoothed_X.push_back(new_X);
        if (iter < max_iterations - 1) {  // Skip for last iteration
            result.smoothed_graphs.push_back(new_graph);
        }

        // 7. Update iterations performed
        result.iterations_performed = iter + 1;

        // 8. Update progress
        if (verbose) {
            elapsed_time(iteration_start_time, "Iteration complete", true);
            progress->update(iter + 1);
        }

        // 9. Check for convergence
        if (diff <= convergence_threshold) {
            if (verbose) {
                Rprintf("Convergence reached at iteration %zu: %.6f <= %.6f\n",
                        iter + 1, diff, convergence_threshold);
            }
            break;
        }
    }

    // Final progress and timing
    if (verbose) {
        if (progress) {
            progress->finish();
        }
        elapsed_time(total_start_time, "Total processing time", true);
        Rprintf("Completed %d iterations\n", result.iterations_performed);

        // Print final convergence value
        if (!result.convergence_metrics.empty()) {
            Rprintf("Final convergence value: %.6f\n", result.convergence_metrics.back());
        }

        // Print whether boosting was used
        if (result.used_boosting) {
            Rprintf("Used boosting mode: switched to residual smoothing after %zu iterations\n",
                    switch_to_residuals_after);
        } else {
            Rprintf("Used direct smoothing mode for all iterations\n");
        }
    }

    return result;
}
