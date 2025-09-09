#include <R.h>                     // Rprintf
// Undefine conflicting macros from R headers
#undef length
#undef eval

#include "graph_spectral_filter.hpp"
#include "set_wgraph.hpp"                   // for set_wgraph_t, adjacency_list
#include "klaps_low_pass_smoother.hpp"      // for compute_graph_laplacian_spectrum()
#include "error_utils.h"                    // REPORT_ERROR()
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <chrono>

/**
 * @brief Smooths a signal on a graph using spectral filtering with various Laplacian operators and filter types.
 *
 * @details
 * This function implements a comprehensive framework for spectral filtering of signals defined on graph vertices.
 * It supports multiple types of graph Laplacians, spectral filters, and parameter configurations to achieve
 * different smoothing characteristics.
 *
 * The general spectral filtering process follows these steps:
 * 1. Constructs a specific graph Laplacian operator based on the selected `laplacian_type`
 * 2. Computes the eigendecomposition of this Laplacian (or its transformation)
 * 3. Projects the signal onto the eigenbasis (Graph Fourier Transform)
 * 4. Applies filter weights based on the eigenvalues and selected `filter_type`
 * 5. Reconstructs smoothed signals for a grid of filter parameter values
 * 6. Selects the optimal parameter using Generalized Cross-Validation (GCV)
 *
 * The function provides extensive flexibility through:
 * - Different Laplacian constructions (standard, normalized, kernel-based, etc.)
 * - Various spectral filter types (heat, Gaussian, cubic spline, etc.)
 * - Parameter controls for localization and smoothness
 *
 * @param y Vector of signal values defined on graph vertices (size = n_vertices)
 *
 * @param laplacian_type Type of graph Laplacian to use:
 *     - STANDARD: L = D - A, the combinatorial Laplacian. Provides basic smoothing that minimizes
 *       first differences across edges. Best for general-purpose smoothing on regularly structured graphs.
 *
 *     - NORMALIZED: L_norm = D^(-1/2) L D^(-1/2), the normalized Laplacian. Accounts for varying
 *       vertex degrees, giving more balanced smoothing on irregular graphs. Recommended for graphs with
 *       highly variable connectivity.
 *
 *     - RANDOM_WALK: L_rw = D^(-1) L, the random walk Laplacian. Similar to normalized Laplacian
 *       but with asymmetric normalization. Useful when modeling diffusion processes.
 *
 *     - KERNEL: L_kernel = D_kernel - W_kernel, a Laplacian constructed using distance-based kernel
 *       weights rather than adjacency weights. Provides smoother spectral response, especially on
 *       irregular or noisy graphs. Good for capturing fine geometric structures.
 *
 *     - NORMALIZED_KERNEL: Normalized version of the kernel Laplacian. Combines benefits of
 *       normalization and kernel-based edge weighting. Excellent for irregular graphs where
 *       geometric relationships matter.
 *
 *     - ADAPTIVE_KERNEL: Kernel Laplacian with locally adaptive bandwidth. Adjusts smoothing
 *       automatically based on local graph density. Best for graphs with highly variable density.
 *
 *     - SHIFTED: I - L, the shifted standard Laplacian. Inverts the spectrum to emphasize smooth
 *       components. Particularly effective for chain graphs and 1D signals.
 *
 *     - SHIFTED_KERNEL: I - L_kernel, the shifted kernel Laplacian. Combines kernel-based edge
 *       weighting with spectral inversion. Excellent for chain graphs and when geometric
 *       relationships are important.
 *
 *     - REGULARIZED: L + ε*I, adds small regularization to ensure positive definiteness.
 *       Helps with numerical stability in filtering operations.
 *
 *     - REGULARIZED_KERNEL: L_kernel + ε*I, regularized version of kernel Laplacian.
 *
 *     - MULTI_SCALE: Weighted combination of kernel Laplacians at different scales.
 *       Captures both fine and coarse features simultaneously. Good for signals with
 *       features at multiple scales.
 *
 * @param filter_type Type of spectral filter to apply:
 *     - HEAT: exp(-t*λ), classic heat kernel filter. Provides smooth decay across frequencies
 *       with more pronounced filtering at higher frequencies. General-purpose choice for
 *       most applications.
 *
 *     - GAUSSIAN: exp(-t*λ²), Gaussian spectral filter. More aggressive decay at higher
 *       frequencies than heat kernel. Produces very smooth results with minimal ringing.
 *
 *     - NON_NEGATIVE: exp(-t*max(λ,0)), truncated heat kernel that only attenuates
 *       non-negative eigenvalues. Useful when dealing with signed Laplacians.
 *
 *     - CUBIC_SPLINE: 1/(1+t*λ²), filter that mimics cubic spline behavior. Minimizes second
 *       derivatives, producing the smoothest results while preserving linear trends.
 *       Excellent for chain graphs and when spline-like smoothness is desired.
 *
 *     - EXPONENTIAL: exp(-t*sqrt(λ)), less aggressive decay than heat kernel.
 *       Maintains more mid-frequency components.
 *
 *     - MEXICAN_HAT: λ*exp(-t*λ²), band-pass filter that enhances mid-frequencies
 *       while attenuating both low and high frequencies. Useful for edge detection.
 *
 *     - IDEAL_LOW_PASS: 1 for λ < t, 0 otherwise. Sharp cutoff filter that retains
 *       all frequencies below threshold t. Can introduce ringing but provides strict
 *       frequency separation.
 *
 *     - BUTTERWORTH: 1/(1+(λ/t)^(2*n)), smoother cutoff than ideal filter.
 *       Produces less ringing while maintaining good frequency separation.
 *
 *     - TIKHONOV: 1/(1+t*λ), first-order smoothing filter. Less aggressive than
 *       cubic spline but more than heat kernel at preserving trends.
 *
 *     - POLYNOMIAL: (1-λ/λ_max)^p for λ < λ_max, polynomial decay filter.
 *       Control decay rate with power p. Useful for custom filtering needs.
 *
 *     - INVERSE_COSINE: cos(π*λ/(2*λ_max)), smooth filter with cosine profile.
 *       Provides gentle transition from low to high frequencies.
 *
 *     - ADAPTIVE: Data-driven filter that adapts to signal properties.
 *       Useful when optimal filtering is not known a priori.
 *
 * @param laplacian_power Power to which the Laplacian (or shifted Laplacian) is raised:
 *     - 1: Standard filtering with direct Laplacian. Suitable for most applications.
 *     - 2: Squared Laplacian (L²) minimizes second derivatives, providing smoother results
 *          similar to cubic splines. Preserves linear trends better than power 1.
 *     - Higher values: Increasingly smooth results but may oversmooth detail.
 *       Odd powers maintain sign characteristics of eigenvalues, while even powers
 *       make all eigenvalues positive.
 *
 * @param kernel_params Parameters for kernel-based Laplacian construction:
 *     - tau_factor: Determines kernel bandwidth as a fraction of graph diameter.
 *       Smaller values (0.001-0.01) create highly localized neighborhoods,
 *       larger values (0.1-0.5) create more global influence.
 *
 *     - radius_factor: Multiplier for search radius when finding vertices within
 *       kernel range. Larger values include more vertices in kernel computation
 *       but increase computational cost. Typically 2-5 is sufficient.
 *
 *     - kernel_type: Shape of kernel function (GAUSSIAN, EXPONENTIAL, HEAT, etc.).
 *       Different kernels have different decay characteristics:
 *         * GAUSSIAN: Classic bell curve, smooth decay from center
 *         * EXPONENTIAL: Sharper peak, heavier tails than Gaussian
 *         * HEAT: Similar to Gaussian but with different scaling
 *         * TRICUBE: Compact support, smooth transition to zero
 *         * EPANECHNIKOV: Parabolic shape, optimal in statistical sense
 *         * UNIFORM: Equal weighting within radius, sharp cutoff
 *         * Others: Variations with different smoothness properties
 *
 *     - adaptive: If true, kernel bandwidth adapts to local graph density,
 *       providing consistent behavior across variable-density regions.
 *
 * @param n_evectors Number of Laplacian eigenpairs to compute:
 *     - Higher values provide more accurate spectral representation but increase
 *       computational cost. Typically 100-200 is sufficient for most graphs.
 *     - Should be much smaller than graph size for large graphs (n_vertices >> n_evectors).
 *     - For very precise filtering, may need to increase this value.
 *
 * @param n_candidates Number of filter parameter values (diffusion times t) to evaluate:
 *     - Controls granularity of parameter search. Higher values give more precise
 *       optimal parameter selection but increase computation time.
 *     - Typically 20-50 is sufficient for smooth parameter landscapes.
 *
 * @param log_grid If true, use logarithmically-spaced filter parameter values:
 *     - Logarithmic spacing (true) is usually better as filter response often
 *       varies more dramatically at small parameter values.
 *     - Linear spacing (false) may be better when parameter range is narrow
 *       or when focusing on a specific region of parameter space.
 *
 * @param with_t_predictions If true, store all smoothed signals for each filter parameter:
 *     - Useful for visualization, parameter sensitivity analysis, or creating
 *       animations of smoothing effect as parameter varies.
 *     - Increases memory usage proportionally to n_candidates.
 *
 * @param verbose If true, print progress information during computation:
 *     - Helpful for monitoring long-running computations on large graphs.
 *     - Provides insight into eigendecomposition convergence and parameter search.
 *
 * @return A graph_spectral_filter_t structure containing:
 *   - evalues: Eigenvalues of the Laplacian operator
 *   - evectors: Corresponding eigenvectors
 *   - candidate_ts: Grid of filter parameter values tested
 *   - gcv_scores: Generalized Cross-Validation score for each parameter
 *   - opt_t_idx: Index of the optimal parameter value
 *   - predictions: Smoothed signal at the optimal parameter
 *   - t_predictions: (Optional) All smoothed signals at each parameter value
 *   - laplacian_type: Type of Laplacian used
 *   - filter_type: Type of spectral filter applied
 *   - laplacian_power: Power to which Laplacian was raised
 *   - kernel_params: Kernel parameters used
 *   - compute_time_ms: Computation time in milliseconds
 *   - gcv_min_score: Minimum GCV score achieved
 *
 * @note Recommended parameter combinations for specific use cases:
 *   - For general graph smoothing: STANDARD Laplacian with HEAT filter, power=1
 *   - For spline-like smoothing on chain graphs: STANDARD Laplacian with CUBIC_SPLINE filter, power=2
 *   - For preserving geometric structure: KERNEL Laplacian with GAUSSIAN filter, power=1
 *   - For irregularly connected graphs: NORMALIZED_KERNEL Laplacian with HEAT filter, power=1
 *   - For multi-scale features: MULTI_SCALE Laplacian with HEAT filter, power=1
 *
 * @see compute_graph_laplacian_spectrum_generic for details on Laplacian spectral computation
 * @see kernel_type_t for available kernel functions
 * @see filter_type_t for available spectral filters
 */
/**
 * @brief Implements spectral filtering for graph signals with various Laplacian operators and filter types.
 *
 * This implementation handles all Laplacian types and filter types defined in the documentation,
 * with proper parameter tracking and result storage.
 */
graph_spectral_filter_t
set_wgraph_t::graph_spectral_filter(
	const std::vector<double>& y,
	laplacian_type_t laplacian_type,
	filter_type_t filter_type,
	size_t laplacian_power,
	kernel_params_t kernel_params,
	size_t n_evectors_to_compute,
	size_t n_candidates,
	bool   log_grid,
	bool   with_t_predictions,
	bool   verbose
	) const {
	// Initialize result structure
	graph_spectral_filter_t result;
	size_t n_vertices = adjacency_list.size();

	// Store parameters in result structure
	result.laplacian_type = laplacian_type;
	result.filter_type = filter_type;
	result.laplacian_power = laplacian_power;
	result.kernel_params = kernel_params;

	// Record start time for performance tracking
	auto start_time = std::chrono::steady_clock::now();

	if (verbose) {
		Rprintf("Starting graph_spectral_filter with:\n");
		Rprintf("  - Laplacian type: %s\n", laplacian_type_to_string(laplacian_type).c_str());
		Rprintf("  - Filter type: %d\n", static_cast<int>(filter_type));
		Rprintf("  - Laplacian power: %zu\n", laplacian_power);
		Rprintf("  - Number of vertices: %zu\n", n_vertices);
		Rprintf("  - Number of eigenvectors: %zu\n", n_evectors_to_compute);
	}

	// Step 1: Compute the Laplacian spectrum using the generic function
	// This function handles all types of Laplacians with proper eigenvalue sorting
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> spectrum = compute_graph_laplacian_spectrum_generic(
		y,
		n_evectors_to_compute,
		laplacian_type,
		kernel_params,
		laplacian_power,
		verbose
		);

	// Extract eigenvalues and eigenvectors
	size_t m = spectrum.first.size();
	result.evalues.assign(spectrum.first.data(), spectrum.first.data() + m);
	result.evectors = std::move(spectrum.second);

	// Step 2: Compute Graph Fourier Transform
	Eigen::VectorXd y_ev = Eigen::Map<const Eigen::VectorXd>(y.data(), n_vertices);
	Eigen::VectorXd gft = result.evectors.transpose() * y_ev;

	// Step 3: Generate candidate diffusion times
	double eps = 1e-11;
	double t_max;

	// Determine appropriate t_max based on eigenvalue distribution
	double eval_max = *std::max_element(result.evalues.begin(), result.evalues.end(),
										[](double a, double b) { return std::abs(a) < std::abs(b); });
	// double eval_min = *std::min_element(result.evalues.begin(), result.evalues.end(),
	// 									[](double a, double b) { return std::abs(a) < std::abs(b); });

	// Different filter types need different t_max ranges
	if (filter_type == filter_type_t::HEAT ||
		filter_type == filter_type_t::GAUSSIAN ||
		filter_type == filter_type_t::NON_NEGATIVE ||
		filter_type == filter_type_t::EXPONENTIAL ||
		filter_type == filter_type_t::MEXICAN_HAT) {
		// For exponential decay filters, use -log(eps)/max eigenvalue
		t_max = (eval_max > 0) ? -std::log(eps) / eval_max : 1.0;
	} else if (filter_type == filter_type_t::CUBIC_SPLINE ||
			   filter_type == filter_type_t::TIKHONOV ||
			   filter_type == filter_type_t::BUTTERWORTH) {
		// For rational filters, use 1/epsilon ~ 10^11
		t_max = 1.0 / eps;
	} else if (filter_type == filter_type_t::IDEAL_LOW_PASS) {
		// For ideal lowpass, use the max eigenvalue
		t_max = eval_max;
	} else if (filter_type == filter_type_t::POLYNOMIAL ||
			   filter_type == filter_type_t::INVERSE_COSINE) {
		// For polynomial/cosine filters, just use 1.0 as they work with normalized values
		t_max = 1.0;
	} else if (filter_type == filter_type_t::ADAPTIVE) {
		// For adaptive, compute based on signal and eigenvalue statistics
		// This is a simple heuristic - could be improved
		double signal_variance = (y_ev.array() - y_ev.mean()).square().sum() / n_vertices;
		t_max = signal_variance / (eval_max * eval_max);
	} else {
		// Default case
		t_max = (eval_max > 0) ? -std::log(eps) / eval_max : 1.0;
	}

	if (verbose) {
		Rprintf("  - Maximum eigenvalue: %.6e\n", eval_max);
		Rprintf("  - Maximum filter parameter: %.6e\n", t_max);
	}

	// Generate diffusion time grid
	result.candidate_ts.reserve(n_candidates);
	for (size_t j = 0; j < n_candidates; ++j) {
		// Skip t=0 case by starting at j+1
		double frac = double(j + 1) / double(n_candidates);
		double t = log_grid
			? std::exp(std::log(eps) + std::log(t_max / eps) * frac)
			: frac * t_max;
		result.candidate_ts.push_back(t);
	}

	// Step 4: Loop over candidate parameters, compute GCV scores and predictions
	result.gcv_scores.resize(n_candidates);
	if (with_t_predictions) {
		result.t_predictions.assign(n_candidates, std::vector<double>(n_vertices));
	}

	for (size_t idx = 0; idx < n_candidates; ++idx) {
		double t = result.candidate_ts[idx];

		// Compute spectral weights based on filter type
		Eigen::ArrayXd w;
		Eigen::ArrayXd lambda = Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m);

		// Apply filter function based on selected type
		switch (filter_type) {
		case filter_type_t::HEAT:
			// Heat kernel filter
			w = (-t * lambda).exp();
			break;

		case filter_type_t::GAUSSIAN:
			// Gaussian filter
			w = (-t * lambda.square()).exp();
			break;

		case filter_type_t::NON_NEGATIVE:
			// Non-negative truncated filter
			lambda = lambda.max(0); // Replace negative values with 0
			w = (-t * lambda).exp();
			break;

		case filter_type_t::CUBIC_SPLINE:
			// Cubic spline-like filter
			w = 1.0 / (1.0 + t * lambda.square());
			break;

		case filter_type_t::EXPONENTIAL:
			// Exponential filter
			w = (-t * lambda.abs().sqrt()).exp();
			break;

		case filter_type_t::MEXICAN_HAT:
			// Mexican hat (wavelet) filter
			w = lambda * (-t * lambda.square()).exp();
			break;

		case filter_type_t::IDEAL_LOW_PASS:
			// Ideal low-pass filter
			w = (lambda < t).cast<double>();
			break;

		case filter_type_t::BUTTERWORTH:
			// Butterworth filter (order 2)
		{
			double n = 2.0; // Filter order
			w = 1.0 / (1.0 + Eigen::pow(lambda / t, 2 * n));
		}
		break;

		case filter_type_t::TIKHONOV:
			// Tikhonov filter
			w = 1.0 / (1.0 + t * lambda);
			break;

		case filter_type_t::POLYNOMIAL:
			// Polynomial filter
		{
			double p = 3.0; // Polynomial degree
			w = ((lambda < eval_max).cast<double>()).cwiseProduct(
				Eigen::pow(1.0 - lambda / eval_max, p)
				);
		}
		break;

 		case filter_type_t::INVERSE_COSINE:
			// Inverse cosine filter
			w = ((lambda < eval_max).cast<double>()).cwiseProduct(
				Eigen::cos(M_PI * lambda / (2 * eval_max))
				);
			break;

		case filter_type_t::ADAPTIVE:
			// Adaptive filter (simple version - could be more sophisticated)
		{
			double threshold = eval_max * 0.1; // 10% threshold
			Eigen::ArrayXd exp_term = (-t * (lambda - threshold).square()).exp();
			w = exp_term + ((lambda < threshold).cast<double>()).cwiseProduct(
				(1.0 - exp_term)
				);
		}
		break;

		default:
			REPORT_ERROR("Unsupported filter_type");
		}

		// Compute smoothed signal
		Eigen::VectorXd y_t = result.evectors * (w * gft.array()).matrix();

		// Calculate GCV score: ||y - y_t||^2 / (n - trace(S_t))^2
		double trS = w.sum(); // Trace of smoothing matrix
		double norm2 = (y_ev - y_t).squaredNorm();
		double denom = double(n_vertices) - trS;

		// Handle edge case where trace is close to n_vertices
		if (std::abs(denom) < 1e-10) {
			denom = 1e-10; // Avoid division by near-zero
		}

		result.gcv_scores[idx] = norm2 / (denom * denom);

		// Store predictions for this parameter value
		if (with_t_predictions) {
			auto &buf = result.t_predictions[idx];
			for (size_t v = 0; v < n_vertices; ++v) {
				buf[v] = y_t[v];
			}
		}
	}

	// Step 5: Find optimal parameter
	result.opt_t_idx = std::min_element(
		result.gcv_scores.begin(),
		result.gcv_scores.end()
		) - result.gcv_scores.begin();

	result.gcv_min_score = result.gcv_scores[result.opt_t_idx];

	if (verbose) {
		Rprintf("  - Optimal parameter index: %zu\n", result.opt_t_idx);
		Rprintf("  - Optimal parameter value: %.6e\n", result.candidate_ts[result.opt_t_idx]);
		Rprintf("  - Minimum GCV score: %.6e\n", result.gcv_min_score);
	}

	// Step 6: Generate final predictions at optimal parameter
	if (with_t_predictions) {
		// If we stored all parameter predictions, just copy the optimal one
		result.predictions = result.t_predictions[result.opt_t_idx];
	} else {
		// Otherwise, compute predictions at optimal parameter
		double t_opt = result.candidate_ts[result.opt_t_idx];

		// Recompute filter weights at optimal parameter
		Eigen::ArrayXd w;
		Eigen::ArrayXd lambda = Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m);

		// Apply same filter function as above (could refactor to avoid duplication)
		switch (filter_type) {
		case filter_type_t::HEAT:
			w = (-t_opt * lambda).exp();
			break;
		case filter_type_t::GAUSSIAN:
			w = (-t_opt * lambda.square()).exp();
			break;
		case filter_type_t::NON_NEGATIVE:
			lambda = lambda.max(0);
			w = (-t_opt * lambda).exp();
			break;
		case filter_type_t::CUBIC_SPLINE:
			w = 1.0 / (1.0 + t_opt * lambda.square());
			break;
		case filter_type_t::EXPONENTIAL:
			w = (-t_opt * lambda.abs().sqrt()).exp();
			break;
		case filter_type_t::MEXICAN_HAT:
			w = lambda * (-t_opt * lambda.square()).exp();
			break;
		case filter_type_t::IDEAL_LOW_PASS:
			w = (lambda < t_opt).cast<double>();
			break;
		case filter_type_t::BUTTERWORTH:
		{
			double n = 2.0; // Filter order
			w = 1.0 / (1.0 + Eigen::pow(lambda / t_opt, 2 * n));
		}
		break;
		case filter_type_t::TIKHONOV:
			w = 1.0 / (1.0 + t_opt * lambda);
			break;
		case filter_type_t::POLYNOMIAL:
		{
			double p = 3.0; // Polynomial degree
			w = ((lambda < eval_max).cast<double>()).cwiseProduct(
				Eigen::pow(1.0 - lambda / eval_max, p)
				);
		}
		break;
		case filter_type_t::INVERSE_COSINE:
			w = ((lambda < eval_max).cast<double>()).cwiseProduct(
				Eigen::cos(M_PI * lambda / (2 * eval_max))
				);
			break;
		case filter_type_t::ADAPTIVE:
		{
			double threshold = eval_max * 0.1;
			Eigen::ArrayXd exp_term = (-t_opt * (lambda - threshold).square()).exp();
			w = exp_term + ((lambda < threshold).cast<double>()).cwiseProduct(
				(1.0 - exp_term)
				);
		}
		break;
		default:
			REPORT_ERROR("Unsupported filter_type");
		}

		// Compute final smoothed signal
		Eigen::VectorXd y_opt = result.evectors * (w * gft.array()).matrix();

		// Store in result
		result.predictions.resize(n_vertices);
		for (size_t v = 0; v < n_vertices; ++v) {
			result.predictions[v] = y_opt[v];
		}
	}

	// Record computation time
	auto end_time = std::chrono::steady_clock::now();
	result.compute_time_ms = std::chrono::duration<double, std::milli>(
		end_time - start_time).count();

	if (verbose) {
		Rprintf("  - Computation time: %.2f ms\n", result.compute_time_ms);
	}

	return result;
}
