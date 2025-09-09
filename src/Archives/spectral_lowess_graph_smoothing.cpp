#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <chrono>
#include <memory>
#include <Rinternals.h>  // For Rprintf
#include <R_ext/Print.h>

#include "set_wgraph.hpp"
#include "iknn_graphs_r.h"
#include "graph_spectral_lowess_mat.hpp"
#include "progress_utils.hpp"
#include "iknn_graphs.hpp"
#include "spectral_lowess_graph_smoothing.hpp"

/**
 * @brief Computes the absolute difference between two data matrices
 * 
 * @param X1 First data matrix [n_features][n_samples]
 * @param X2 Second data matrix [n_features][n_samples]
 * @param convergence_type Type of convergence criteria to use
 * @return double The computed difference metric according to the specified convergence criteria
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
 * @brief Creates an ikNN graph from a data matrix using the specified parameters
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
 * @brief Implements iterative spectral LOWESS graph smoothing
 * 
 * This function iteratively applies spectral LOWESS smoothing to a graph and its associated data matrix.
 * At each iteration, conditional expectations of features are estimated using graph_spectral_lowess_mat,
 * and a new graph is constructed from the smoothed data.
 * 
 * The process continues until either the maximum number of iterations is reached or
 * the change in the data matrix falls below the specified convergence threshold.
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
    spectral_lowess_graph_smoothing_t result;
    result.iterations_performed = 0;
    
    // Add initial graph and data matrix
    result.smoothed_graphs.push_back(initial_graph);
    result.smoothed_X.push_back(X);
    
    // Progress tracking
    std::unique_ptr<progress_tracker_t> progress;
    if (verbose) {
        progress = std::make_unique<progress_tracker_t>(max_iterations, "Graph Smoothing", 1);
        Rprintf("Starting spectral LOWESS graph smoothing with max %zu iterations\n", max_iterations);
        Rprintf("Convergence threshold: %.6f, Type: %d\n", convergence_threshold, 
                static_cast<int>(convergence_type));
    }
    
    // Main iteration loop
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        auto iteration_start_time = std::chrono::steady_clock::now();
        
        if (verbose) {
            Rprintf("\n--- Iteration %zu/%zu ---\n", iter + 1, max_iterations);
        }
        
        // 1. Get current graph and data
        const set_wgraph_t& current_graph = result.smoothed_graphs.back();
        const std::vector<std::vector<double>>& current_X = result.smoothed_X.back();
        
        // 2. Apply spectral LOWESS to estimate conditional expectations
        if (verbose) {
            Rprintf("Computing spectral LOWESS smoothing...\n");
        }
        
        auto lowess_start_time = std::chrono::steady_clock::now();
        graph_spectral_lowess_mat_t lowess_result = current_graph.graph_spectral_lowess_mat(
            current_X,
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
        std::vector<std::vector<double>> new_X = lowess_result.predictions;
        
        // 4. Create new graph from smoothed data
        set_wgraph_t new_graph;
        if (iter < max_iterations - 1) {  // Skip for last iteration to save time
            new_graph = create_iknn_graph_from_matrix(new_X, k, pruning_thld, verbose);

            // Find the number of connected components of the new graph
            size_t ncc = new_graph.count_connected_components();
            if (ncc > 1) {
                Rprintf("A new graph has %zu connected components. Aborting the smoothing process.\n", ncc);
                break;
            }

            // Find diameter endpoints
            auto [end1, diam] = new_graph.get_vertex_eccentricity(0);  // Start from vertex 0
            auto [end2, diameter] = new_graph.get_vertex_eccentricity(end1);
            new_graph.graph_diameter = diameter;
        }
        
        // 5. Compute convergence metric
        double diff = compute_difference(current_X, new_X, convergence_type);
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
    }
    
    return result;
}
