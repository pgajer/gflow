#ifndef SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP
#define SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP

#include "set_wgraph.hpp"

#include <vector>
#include <cstddef>
using std::size_t;

/**
 * @brief Defines the type of convergence criteria to use for the iterative process
 */
enum class ConvergenceCriteriaType {
    MAX_ABSOLUTE_DIFF,   // Maximum absolute difference across all features
    MEAN_ABSOLUTE_DIFF,  // Mean absolute difference across all features
    RELATIVE_CHANGE      // Maximum relative change across all features
};

/**
 * @brief Structure containing the results of the spectral LOWESS graph smoothing process
 */
struct spectral_lowess_graph_smoothing_t {
    std::vector<set_wgraph_t> smoothed_graphs;                // Sequence of smoothed graphs at each iteration
    std::vector<std::vector<std::vector<double>>> smoothed_X; // Sequence of smoothed data matrices at each iteration

    std::vector<std::vector<size_t>> vertex_mappings;         // For each iteration, maps current vertices to original indices
    std::vector<size_t> outlier_indices;                      // Stores indices of vertices identified as outliers

    std::vector<double> convergence_metrics;                  // Convergence metric at each iteration
    int iterations_performed;                                 // Actual number of iterations performed
    bool used_boosting;                                       // Whether boosting was used
};

/**
 * @brief Implements iterative spectral LOWESS graph smoothing
 * 
 * This function iteratively applies spectral LOWESS smoothing to a graph and its associated data matrix.
 * At each iteration, conditional expectations of features are estimated using graph_spectral_lowess_mat,
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
    bool verbose);

/**
 * @brief Helper function to create an ikNN graph from a data matrix
 */
set_wgraph_t create_iknn_graph_from_matrix(
    const std::vector<std::vector<double>>& X,
    size_t k,
    double max_path_edge_ratio_thld,
    double path_edge_ratio_percentile,
    bool verbose);

/**
 * @brief Computes the difference between two data matrices according to the specified criteria
 */
double compute_difference(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type);


/**
 * @brief Compute the difference between two data matrices according to specified criteria
 */
double compute_difference_with_mapping(
    const std::vector<std::vector<double>>& X1,
    const std::vector<std::vector<double>>& X2,
    ConvergenceCriteriaType convergence_type,
    const std::vector<size_t>& current_mapping,
    const std::vector<size_t>& previous_mapping);

/**
 * @brief Expand a reduced data matrix back to original size using a vertex mapping
 */
std::vector<std::vector<double>> expand_to_original_size(
    const std::vector<std::vector<double>>& reduced_X,
    const std::vector<size_t>& mapping,
    size_t original_size);

#endif // SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP
