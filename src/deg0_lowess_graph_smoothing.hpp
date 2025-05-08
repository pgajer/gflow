#ifndef DEG0_LOWESS_GRAPH_SMOOTHING_HPP
#define DEG0_LOWESS_GRAPH_SMOOTHING_HPP

#include <vector>
#include "set_wgraph.hpp"
#include "spectral_lowess_graph_smoothing.hpp"  // For ConvergenceCriteriaType

/**
 * @brief Structure containing the results of the degree 0 LOWESS graph smoothing process
 */
struct deg0_lowess_graph_smoothing_t {
    std::vector<set_wgraph_t> smoothed_graphs;                // Sequence of smoothed graphs at each iteration
    std::vector<std::vector<std::vector<double>>> smoothed_X; // Sequence of smoothed data matrices at each iteration
    std::vector<std::vector<size_t>> vertex_mappings;         // For each iteration, maps current vertices to original indices
    std::vector<size_t> outlier_indices;                      // Stores indices of vertices identified as outliers
    std::vector<double> convergence_metrics;                  // Convergence metric at each iteration
    int iterations_performed;                                 // Actual number of iterations performed
    bool used_boosting;                                       // Whether boosting was used
};

/**
 * @brief Implements iterative degree 0 LOWESS graph smoothing
 *
 * This function iteratively applies degree 0 LOWESS smoothing to a graph and its associated data matrix.
 * At each iteration, conditional expectations of features are estimated using local weighted averages,
 * and a new graph is constructed from the smoothed data.
 *
 * Supports two modes of operation:
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
 * @param n_bws Number of candidate bandwidths
 * @param log_grid Whether to use logarithmic spacing for bandwidth grid
 * @param min_bw_factor Factor for minimum bandwidth
 * @param max_bw_factor Factor for maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances
 * @param kernel_type Type of kernel function to use
 * @param n_folds Number of cross-validation folds
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 * @param n_cleveland_iterations Number of robustness iterations
 * @param with_bw_predictions Whether to store predictions for all bandwidths
 * @param switch_to_residuals_after Number of direct smoothing iterations before switching to residual smoothing.
 *                                 Set to 0 to always use residual smoothing (boosting).
 *                                 Set to a value >= max_iterations to never switch (traditional smoothing).
 * @param verbose Whether to print progress information
 * @return deg0_lowess_graph_smoothing_t Structure containing results of the smoothing process
 */
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
    bool verbose);

// Declaration of helper functions
set_wgraph_t create_iknn_graph_from_matrix(
    const std::vector<std::vector<double>>& X,
    size_t k,
    double pruning_thld,
    bool verbose);

double compute_difference(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type);

double compute_difference_with_mapping(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type,
    const std::vector<size_t>& current_mapping,
    const std::vector<size_t>& previous_mapping);

std::vector<std::vector<double>> expand_to_original_size(
    const std::vector<std::vector<double>>& reduced_X,
    const std::vector<size_t>& mapping,
    size_t original_size);

#endif // DEG0_LOWESS_GRAPH_SMOOTHING_HPP
