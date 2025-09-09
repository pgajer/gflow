#ifndef SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP
#define SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP

#include <vector>
#include "set_wgraph.hpp"

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
    std::vector<set_wgraph_t> smoothed_graphs;              // Sequence of smoothed graphs at each iteration
    std::vector<std::vector<std::vector<double>>> smoothed_X; // Sequence of smoothed data matrices at each iteration
    std::vector<double> convergence_metrics;                // Convergence metric at each iteration
    int iterations_performed;                              // Actual number of iterations performed
};

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
    bool verbose);

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
    bool verbose);

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
    ConvergenceCriteriaType convergence_type);

#endif // SPECTRAL_LOWESS_GRAPH_SMOOTHING_HPP
