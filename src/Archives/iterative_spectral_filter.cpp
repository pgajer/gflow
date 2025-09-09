#include "set_wgraph.hpp"


struct iterative_spectral_filter_t {
	std::vector<double> predictions;
	std::vector<std::vector<double>> iterations_history;
	std::vector<double> convergence_history;
	bool converged;
	int iterations_performed;
	set_wgraph_t final_graph;
};

/**
 * @brief Performs iterative spectral smoothing with function-aware graph adaptation
 *
 * @details
 * This function implements an iterative smoothing algorithm that alternates between:
 * 1. Computing a spectral estimate of the function on the current graph
 * 2. Updating the graph weights based on the current function estimate
 * The process continues until convergence or a maximum number of iterations.
 *
 * @param y Original function values
 * @param weight_type Type of weight modification to use (see construct_function_aware_graph)
 * @param weight_params Parameters for the weight modification function
 * @param n_iterations Maximum number of iterations
 * @param tol Convergence tolerance (algorithm stops when the change in function estimate is below this)
 * @param alpha Interpolation factor between iterations (0-1, where 1 means use only new estimate)
 * @param normalize Whether to normalize weights after modification
 * @param spectral_params Parameters for spectral estimation:
 *   spectral_params[0] = laplacian_type (see graph_spectral_filter)
 *   spectral_params[1] = filter_type
 *   spectral_params[2] = laplacian_power
 *   spectral_params[3] = kernel_tau_factor
 *   spectral_params[4] = n_evectors_to_compute
 * @param verbose Whether to print progress information
 * @return A struct containing:
 *   - predictions: Final smoothed function values
 *   - iterations_history: Function estimates at each iteration
 *   - convergence_history: Change in function estimates between iterations
 *   - converged: Whether the algorithm converged
 *   - iterations_performed: Number of iterations performed
 */
iterative_spectral_filter_t set_wgraph_t::iterative_spectral_filter(
	const std::vector<double>& y,
	int weight_type,
	const std::vector<double>& weight_params,
	int n_iterations,
	double tol,
	double alpha,
	bool normalize,
	const std::vector<double>& spectral_params,
	bool verbose) const {

	// Initialize result structure
	iterative_spectral_filter_t result;
	result.converged = false;
	result.iterations_performed = 0;

	// Initialize history if requested
	result.iterations_history.push_back(y);

	// Initialize current estimate and graph
	std::vector<double> current_estimate = y;
	set_wgraph_t current_graph = *this;

	// Parse spectral parameters or use defaults
	laplacian_type_t laplacian_type = spectral_params.size() > 0
		? static_cast<laplacian_type_t>(static_cast<int>(spectral_params[0]))
		: laplacian_type_t::NORMALIZED;

	filter_type_t filter_type = spectral_params.size() > 1
		? static_cast<filter_type_t>(static_cast<int>(spectral_params[1]))
		: filter_type_t::HEAT;

	size_t laplacian_power = spectral_params.size() > 2
		? static_cast<size_t>(spectral_params[2])
		: 1;

	double kernel_tau_factor = spectral_params.size() > 3
		? spectral_params[3]
		: 0.05;

	size_t n_evectors_to_compute = spectral_params.size() > 4
		? static_cast<size_t>(spectral_params[4])
		: std::min<size_t>(50, adjacency_list.size() / 2);

	// Set up kernel parameters
	kernel_params_t kernel_params;
	kernel_params.tau_factor = kernel_tau_factor;

	// Main iteration loop
	for (int iter = 0; iter < n_iterations; ++iter) {
		if (verbose) {
			Rprintf("Iteration %d/%d\n", iter + 1, n_iterations);
		}

		// Step 1: Compute spectral estimate on current graph
		graph_spectral_filter_t filter_result = current_graph.graph_spectral_filter(
			current_estimate,
			laplacian_type,
			filter_type,
			laplacian_power,
			kernel_params,
			n_evectors_to_compute,
			40,  // n_candidates
			true, // log_grid
			false, // with_t_predictions
			verbose
			);

		// Get the predictions from the filter
		std::vector<double> new_estimate = filter_result.predictions;

		// Step 2: Interpolate between current and new estimate if alpha < 1
		if (alpha < 1.0) {
			for (size_t i = 0; i < new_estimate.size(); ++i) {
				new_estimate[i] = (1.0 - alpha) * current_estimate[i] + alpha * new_estimate[i];
			}
		}

		// Step 3: Compute change for convergence check
		double max_change = 0.0;
		for (size_t i = 0; i < new_estimate.size(); ++i) {
			max_change = std::max(max_change, std::abs(new_estimate[i] - current_estimate[i]));
		}

		result.convergence_history.push_back(max_change);

		if (verbose) {
			Rprintf("  Max change: %g\n", max_change);
		}

		// Step 4: Check for convergence
		if (max_change < tol) {
			if (verbose) {
				Rprintf("Converged after %d iterations\n", iter + 1);
			}
			result.converged = true;
			result.iterations_performed = iter + 1;
			current_estimate = new_estimate;
			break;
		}

		// Step 5: Update graph based on new estimate
		current_graph = current_graph.construct_function_aware_graph(
			new_estimate,
			weight_type,
			weight_params,
			normalize
			);

		// Step 6: Update current estimate
		current_estimate = new_estimate;

		// Record history if requested
		result.iterations_history.push_back(current_estimate);

		result.iterations_performed = iter + 1;
	}

	// Set final predictions
	result.predictions = current_estimate;
	result.final_graph = current_graph;

	return result;
}
