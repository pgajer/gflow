#include <Rinternals.h>  // For Rprintf
#include <R_ext/Print.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

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
#include "graph_spectral_lowess_mat.hpp"
#include "progress_utils.hpp"
#include "spectral_lowess_graph_smoothing.hpp"

/**
 * @brief Implements iterative spectral LOWESS graph smoothing with optional boosting
 * 
 * This function iteratively applies spectral LOWESS smoothing to a graph and its associated data matrix.
 * It supports two modes of operation:
 * 1. Direct smoothing: Features from the data matrix are smoothed directly
 * 2. Residual smoothing (boosting): After a specified number of iterations, switches to 
 *    smoothing the residuals between the original data and the current estimates
 * 
 * @param initial_graph Initial graph structure
 * @param X Initial data matrix [n_features][n_samples]
 * @param max_iterations Maximum number of iterations to perform
 * @param convergence_threshold Threshold for convergence
 * @param convergence_type Type of convergence criteria to use
 * @param k Number of nearest neighbors for kNN graph construction
 * @param pruning_thld Threshold for pruning edges
 * @param n_evectors Number of eigenvectors for spectral embedding
 * @param n_bws Number of candidate bandwidths
 * @param log_grid Whether to use logarithmic spacing for bandwidth grid
 * @param min_bw_factor Factor for minimum bandwidth
 * @param max_bw_factor Factor for maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances
 * @param kernel_type Type of kernel function to use
 * @param n_cleveland_iterations Number of robustness iterations
 * @param compute_errors Whether to compute prediction errors
 * @param compute_scales Whether to compute bandwidth/scale info
 * @param switch_to_residuals_after Number of direct smoothing iterations before switching to residual smoothing.
 *                                 Set to 0 to always use residual smoothing (boosting).
 *                                 Set to a value >= max_iterations to never switch (traditional smoothing).
 * @param verbose Whether to print progress information
 * @return spectral_lowess_graph_smoothing_t Structure containing results of the smoothing process
 */
spectral_lowess_graph_smoothing_t spectral_lowess_graph_smoothing(
    const set_wgraph_t& initial_graph,
    const std::vector<std::vector<double>>& X,
    // Iteration control parameters
    size_t max_iterations,
    double convergence_threshold,
    ConvergenceCriteriaType convergence_type,
    // Graph construction parameters
    size_t k,
    double pruning_thld,
    // graph_spectral_lowess_mat parameters
    size_t n_evectors,
    // Bandwidth parameters
    size_t n_bws,
    bool log_grid,
    double min_bw_factor,
    double max_bw_factor,
    // Kernel parameters
    double dist_normalization_factor,
    size_t kernel_type,
    // Other
    size_t n_cleveland_iterations,
    bool compute_errors,
    bool compute_scales,
    // Boosting parameters
    size_t switch_to_residuals_after,
    bool verbose
    ) {
    
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
    spectral_lowess_graph_smoothing_t result;
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
        progress = std::make_unique<progress_tracker_t>(max_iterations, "Graph Smoothing", 1);
        Rprintf("Starting spectral LOWESS graph smoothing with max %zu iterations\n", max_iterations);
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
        const set_wgraph_t& current_graph = result.smoothed_graphs.back();
        
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
        
        // 2. Apply spectral LOWESS to estimate conditional expectations
        if (verbose) {
            Rprintf("Computing spectral LOWESS smoothing...\n");
        }
        
        auto lowess_start_time = std::chrono::steady_clock::now();
        graph_spectral_lowess_mat_t lowess_result = current_graph.graph_spectral_lowess_mat(
            current_input_X,
            n_evectors,
            n_bws,
            log_grid,
            min_bw_factor,
            max_bw_factor,
            dist_normalization_factor,
            kernel_type,
            1e-6,  // precision
            n_cleveland_iterations,
            compute_errors,
            compute_scales,
            verbose
        );
        
        if (verbose) {
            elapsed_time(lowess_start_time, "Spectral LOWESS complete", true);
        }
        
        // 3. Create new data matrix from LOWESS predictions
        std::vector<std::vector<double>> new_smoothed_component = lowess_result.predictions;
        
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

                // Check if the number of outliers does not exceed the follwing outlier_thld limit
                size_t n_vertices_of_new_graph = new_graph.num_vertices();
                double outlier_thld = 0.05;
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
        // Use the new compute_difference_with_mapping function when we have different sized matrices
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

/**
 * @brief Computes the difference between two data matrices according to the specified criteria
 *
 * @param X1 First data matrix [n_features][n_samples]
 * @param X2 Second data matrix [n_features][n_samples]
 * @param convergence_type Type of convergence criteria to use
 * @return double The computed difference metric
 */
double compute_difference(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type) {
    
    if (X1.empty() || X2.empty() || X1.size() != X2.size() || 
        X1[0].size() != X2[0].size()) {
        REPORT_ERROR("Invalid input matrices for difference computation");
    }
    
    size_t n_features = X1.size();
    size_t n_samples = X1[0].size();
    
    switch (convergence_type) {
        case ConvergenceCriteriaType::MAX_ABSOLUTE_DIFF: {
            double max_diff = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (size_t i = 0; i < n_samples; ++i) {
                    double diff = std::abs(X1[j][i] - X2[j][i]);
                    max_diff = std::max(max_diff, diff);
                }
            }
            return max_diff;
        }
        
        case ConvergenceCriteriaType::MEAN_ABSOLUTE_DIFF: {
            double sum_diff = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (size_t i = 0; i < n_samples; ++i) {
                    sum_diff += std::abs(X1[j][i] - X2[j][i]);
                }
            }
            return sum_diff / (n_features * n_samples);
        }
        
        case ConvergenceCriteriaType::RELATIVE_CHANGE: {
            double max_rel_change = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (size_t i = 0; i < n_samples; ++i) {
                    double abs_val = std::abs(X1[j][i]);
                    if (abs_val > std::numeric_limits<double>::epsilon()) {
                        double rel_change = std::abs(X1[j][i] - X2[j][i]) / abs_val;
                        max_rel_change = std::max(max_rel_change, rel_change);
                    } else if (std::abs(X2[j][i]) > std::numeric_limits<double>::epsilon()) {
                        // If X1 is zero but X2 isn't, this is a large relative change
                        max_rel_change = 1.0;  // Consider it 100% change
                    }
                }
            }
            return max_rel_change;
        }
        
        default:
            REPORT_ERROR("Unknown convergence criteria type");
            return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * @brief Compute the difference between two data matrices according to specified criteria
 *
 * This version handles matrices of different sizes by using the vertex mapping to
 * compare only the shared vertices between both matrices.
 *
 * @param X1 First data matrix
 * @param X2 Second data matrix (possibly with different dimensions)
 * @param convergence_type Type of convergence criteria to use
 * @param current_mapping Mapping of vertices in X2 to their original indices
 * @param previous_mapping Mapping of vertices in X1 to their original indices
 * @return double The computed difference according to the specified criteria
 */
double compute_difference_with_mapping(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type,
    const std::vector<size_t>& current_mapping,
    const std::vector<size_t>& previous_mapping) {

    if (X1.empty() || X2.empty() || X1.size() != X2.size()) {
        REPORT_ERROR("Invalid input matrices for difference computation");
    }

    // Find common vertices between the two mappings
    std::unordered_set<size_t> current_original_indices(current_mapping.begin(), current_mapping.end());

    // Create mappings for lookup: original_index -> position in each matrix
    std::unordered_map<size_t, size_t> prev_idx_map;
    for (size_t i = 0; i < previous_mapping.size(); ++i) {
        prev_idx_map[previous_mapping[i]] = i;
    }

    std::unordered_map<size_t, size_t> curr_idx_map;
    for (size_t i = 0; i < current_mapping.size(); ++i) {
        curr_idx_map[current_mapping[i]] = i;
    }

    // Create pairs of indices to compare (positions in each matrix)
    std::vector<std::pair<size_t, size_t>> comparison_indices;
    for (size_t orig_idx : current_original_indices) {
        if (prev_idx_map.find(orig_idx) != prev_idx_map.end()) {
            // This vertex exists in both matrices
            comparison_indices.push_back({prev_idx_map[orig_idx], curr_idx_map[orig_idx]});
        }
    }

    if (comparison_indices.empty()) {
        // No common vertices - can't compute difference
        REPORT_ERROR("No common vertices between iterations for difference computation");
    }

    size_t n_features = X1.size();
    size_t n_common_vertices = comparison_indices.size();

    switch (convergence_type) {
        case ConvergenceCriteriaType::MAX_ABSOLUTE_DIFF: {
            double max_diff = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (const auto& [idx1, idx2] : comparison_indices) {
                    double diff = std::abs(X1[j][idx1] - X2[j][idx2]);
                    max_diff = std::max(max_diff, diff);
                }
            }
            return max_diff;
        }

        case ConvergenceCriteriaType::MEAN_ABSOLUTE_DIFF: {
            double sum_diff = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (const auto& [idx1, idx2] : comparison_indices) {
                    sum_diff += std::abs(X1[j][idx1] - X2[j][idx2]);
                }
            }
            return sum_diff / (n_features * n_common_vertices);
        }

        case ConvergenceCriteriaType::RELATIVE_CHANGE: {
            double max_rel_change = 0.0;
            for (size_t j = 0; j < n_features; ++j) {
                for (const auto& [idx1, idx2] : comparison_indices) {
                    double abs_val = std::abs(X1[j][idx1]);
                    if (abs_val > std::numeric_limits<double>::epsilon()) {
                        double rel_change = std::abs(X1[j][idx1] - X2[j][idx2]) / abs_val;
                        max_rel_change = std::max(max_rel_change, rel_change);
                    } else if (std::abs(X2[j][idx2]) > std::numeric_limits<double>::epsilon()) {
                        // If X1 is zero but X2 isn't, this is a large relative change
                        max_rel_change = 1.0;  // Consider it 100% change
                    }
                }
            }
            return max_rel_change;
        }

        default:
            REPORT_ERROR("Unknown convergence criteria type");
            return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * @brief Helper function to create an ikNN graph from a data matrix
 *
 * @param X Data matrix [n_features][n_samples]
 * @param k Number of nearest neighbors
 * @param pruning_thld Threshold for pruning edges
 * @param verbose Whether to print progress information
 * @return set_wgraph_t The constructed ikNN graph
 */
set_wgraph_t create_iknn_graph_from_matrix(
    const std::vector<std::vector<double>>& X,
    size_t k,
    double pruning_thld,
    bool verbose) {

    if (X.empty()) {
        REPORT_ERROR("Empty data matrix provided for graph construction");
    }

    size_t n_features = X.size();
    size_t n_samples = X[0].size();

    // Create R matrix from X (transpose to match R's column-major format)
    SEXP r_matrix = PROTECT(allocMatrix(REALSXP, n_samples, n_features));
    double* data = REAL(r_matrix);

    for (size_t j = 0; j < n_features; ++j) {
        if (X[j].size() != n_samples) {
            PROTECT(r_matrix); // Ensure protection before error
            REPORT_ERROR("Inconsistent sample size across features");
        }
        for (size_t i = 0; i < n_samples; ++i) {
            data[i + j * n_samples] = X[j][i];
        }
    }

    // Create parameters for S_create_single_iknn_graph
    SEXP r_k = PROTECT(ScalarInteger(k));
    SEXP r_pruning_thld = PROTECT(ScalarReal(pruning_thld));
    SEXP r_compute_full = PROTECT(ScalarLogical(0)); // Don't need full components

    // Call S_create_single_iknn_graph
    if (verbose) {
        Rprintf("Creating ikNN graph with k=%zu, pruning_thld=%.4f\n", k, pruning_thld);
    }

    auto start_time = std::chrono::steady_clock::now();
    SEXP r_graph = PROTECT(S_create_single_iknn_graph(r_matrix, r_k, r_pruning_thld, r_compute_full));

    if (verbose) {
        elapsed_time(start_time, "Graph construction complete", true);
    }

    // Extract adjacency list and weights
    SEXP pruned_adj_list = VECTOR_ELT(r_graph, 4);  // pruned_adj_list is at index 4
    SEXP pruned_weight_list = VECTOR_ELT(r_graph, 5);  // pruned_weight_list is at index 5

    // Convert to set_wgraph_t format
    std::vector<std::vector<int>> adj_list(n_samples);
    std::vector<std::vector<double>> weight_list(n_samples);

    for (size_t i = 0; i < n_samples; ++i) {
        SEXP r_neighbors = VECTOR_ELT(pruned_adj_list, i);
        SEXP r_weights = VECTOR_ELT(pruned_weight_list, i);

        int* neighbors = INTEGER(r_neighbors);
        double* weights = REAL(r_weights);

        int n_neighbors = LENGTH(r_neighbors);

        adj_list[i].resize(n_neighbors);
        weight_list[i].resize(n_neighbors);

        for (int j = 0; j < n_neighbors; ++j) {
            // Convert from 1-based indexing (R) to 0-based indexing (C++)
            adj_list[i][j] = neighbors[j] - 1;
            weight_list[i][j] = weights[j];
        }
    }

    set_wgraph_t result(adj_list, weight_list);

    UNPROTECT(5);  // r_matrix, r_k, r_pruning_thld, r_compute_full, r_graph
    return result;
}


/**
 * @brief Expand a reduced data matrix back to original size using a vertex mapping
 *
 * @param reduced_X The reduced data matrix
 * @param mapping Mapping from reduced indices to original indices
 * @param original_size The size of the original matrix
 * @return std::vector<std::vector<double>> The expanded matrix
 */
std::vector<std::vector<double>> expand_to_original_size(
    const std::vector<std::vector<double>>& reduced_X,
    const std::vector<size_t>& mapping,
    size_t original_size) {

    std::vector<std::vector<double>> expanded_X;
    if (reduced_X.empty()) return expanded_X;

    // Initialize with appropriate size
    expanded_X.resize(reduced_X.size(), std::vector<double>(original_size, 0.0));

    // Fill the expanded matrix using the mapping
    for (size_t j = 0; j < reduced_X.size(); ++j) {
        for (size_t i = 0; i < mapping.size(); ++i) {
            size_t original_idx = mapping[i];
            expanded_X[j][original_idx] = reduced_X[j][i];
        }
    }

    return expanded_X;
}
