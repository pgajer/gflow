#include <R.h>                     // Rprintf
// Undefine conflicting macros from R headers
#undef length

#include "klaps_spectral_filter.hpp"
#include "set_wgraph.hpp"                   // for set_wgraph_t, adjacency_list
#include "klaps_low_pass_smoother.hpp"      // for compute_graph_laplacian_spectrum()
#include "error_utils.h"                    // REPORT_ERROR()
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

/**
 * @brief Smooths a signal on a graph using a modified heat kernel filter with the shifted Laplacian.
 *
 * This function implements a spectral filtering approach using a heat kernel-inspired operation
 * adapted for the shifted Laplacian operator (I-L)^k. For a given vertex signal y, it computes
 * a family of smoothed signals by applying exponential weights to the eigencomponents:
 *   y_t = \sum_{i=0}^{m-1} exp(-t * |λ_i|) * <y, v_i> * v_i
 *
 * The method:
 * 1. Computes the eigen-decomposition of the shifted kernel Laplacian (I-L)^k
 * 2. Projects the signal onto the eigenbasis (Graph Fourier Transform)
 * 3. Applies exponential decay weights based on diffusion time t and eigenvalue magnitudes
 * 4. Reconstructs smoothed signals for a grid of candidate diffusion times
 * 5. Selects the optimal t using Generalized Cross-Validation (GCV)
 *
 * Note: We use |λ_i| (absolute values of eigenvalues) to ensure numerical stability,
 * as the shifted Laplacian can have negative eigenvalues.
 *
 * @param y                      Vertex signal to be smoothed (size = n_vertices)
 *
 * @param n_evectors_to_compute  Number of Laplacian eigenpairs to compute (m)
 *
 * @param tau_factor
 *     Positive real scaling factor that determines the kernel bandwidth \(\tau\) as a fraction of the graph diameter.
 *     The kernel bandwidth is computed as:
 *         \[
 *         \tau = \text{tau\_factor} \times \text{graph diameter}
 *         \]
 *     Smaller `tau_factor` values result in more localized smoothing; larger values make smoothing more global.
 *
 * @param radius_factor A real scaling factor of tau radius that is not less than 1.
 *
 * @param laplacian_power
 *     Positive odd integer specifying the power to which (I - L) is raised.
 *     Higher values apply stronger smoothing by repeatedly reinforcing low-pass filtering.
 *     Typically, larger `laplacian_power` requires smaller `tau_factor` to maintain the smoothing scale.
 *
 * @param n_candidates           Number of diffusion times t to evaluate
 * @param log_grid               If true, use logarithmically-spaced t values; otherwise linear spacing
 * @param with_t_predictions     If true, store all smoothed signals y_t for each diffusion time
 * @param verbose                If true, print progress information
 *
 * @return A klaps_spectral_filter_t containing:
 *   - evalues: Eigenvalues of the shifted Laplacian
 *   - evectors: Corresponding eigenvectors
 *   - candidate_ts: Grid of diffusion times tested
 *   - gcv_scores: GCV score for each diffusion time
 *   - opt_t_idx: Index of the optimal diffusion time
 *   - predictions: Smoothed signal at the optimal diffusion time
 *   - t_predictions: (Optional) All smoothed signals at each diffusion time
 */
klaps_spectral_filter_t
set_wgraph_t::klaps_spectral_filter(
	const std::vector<double>& y,
	filter_type_t filter_type,
	size_t n_evectors_to_compute,
	double tau_factor,
	double radius_factor,
	size_t laplacian_power,
	size_t n_candidates,
	bool   log_grid,
	bool   with_t_predictions,
	bool   verbose
	) const {

	klaps_spectral_filter_t result;
	size_t n_vertices = adjacency_list.size();

	// 1) Eigen‐decomposition
	double tau = tau_factor * graph_diameter;
	auto [evals, V] = compute_graph_shifted_kernel_laplacian_spectrum(
		n_evectors_to_compute,
		tau,
		radius_factor,
		laplacian_power,
		verbose
		);

	size_t m = evals.size();

	// 2) Sort ascending
	std::vector<std::pair<double, Eigen::VectorXd>> pairs;
	pairs.reserve(m);
	for (size_t i = 0; i < m; ++i)
		pairs.emplace_back(evals[i], V.col(i));
	std::sort(pairs.begin(), pairs.end(),
			  [](auto &a, auto &b){ return a.first < b.first; });
	for (size_t i = 0; i < m; ++i) {
		evals[i]= pairs[i].first;
		V.col(i)= pairs[i].second;
	}
	result.evalues.assign(evals.data(), evals.data() + m);
	result.evectors = std::move(V);

	// 3) Graph Fourier coefficients
	Eigen::VectorXd y_ev = Eigen::Map<const Eigen::VectorXd>(y.data(), n_vertices);
	Eigen::VectorXd gft  = result.evectors.transpose() * y_ev;

	// 4) Build diffusion‐time grid in [0, t_max]
	double eps = 1e-11;
	double eval_max = result.evalues.back();
	double t_max = (eval_max > 0)
		? -std::log(eps) / eval_max
		: 1.0;

	if (verbose) {
		// Rprintf("t_max: %.3f\n", t_max);
	}

	result.candidate_ts.reserve(n_candidates);
	for (size_t j = 0; j < n_candidates; ++j) {
		// Skip t=0 case
		double frac = double(j + 1) / double(n_candidates);
		double t = log_grid
			? std::exp(std::log(eps) + std::log(t_max / eps) * frac)
			: frac * t_max;
		result.candidate_ts.push_back(t);
	}

	// 5) Loop over t_j, compute GCV and (optionally) store y_{t_j}
	result.gcv_scores.resize(n_candidates);
	if (with_t_predictions)
		result.t_predictions.assign(n_candidates, std::vector<double>(n_vertices));

	for (size_t idx = 0; idx < n_candidates; ++idx) {
		double t = result.candidate_ts[idx];

		// Compute spectral weights based on filter type
		Eigen::ArrayXd w;
		if (filter_type == filter_type_t::HEAT) {
			// Option 1: Standard heat kernel filter
			w = (-t * Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m)).exp();
		} else if (filter_type == filter_type_t::GAUSSIAN) {
			// Option 2: Gaussian filter (squared eigenvalues)
			w = (-t * Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m).square()).exp();
		} else if (filter_type == filter_type_t::NON_NEGATIVE) {
			// Option 3: Non-negative truncated heat kernel
			Eigen::ArrayXd lambda = Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m);
			lambda = lambda.max(0); // Replace negative values with 0
			w = (-t * lambda).exp();
		} else {
			REPORT_ERROR("Unsupported filter_type. Must be 'HEAT', 'GAUSSIAN', or 'NON_NEGATIVE'");
		}

		// y_t = V * (w * gft)
		Eigen::VectorXd y_t = result.evectors * (w * gft.array()).matrix();

		// GCV: ||y - y_t||^2 / (n - trace(S_t))^2, with trace(S_t) = sum_i w_i
		double trS = w.sum();
		double norm2 = (y_ev - y_t).squaredNorm();
		double denom = double(n_vertices) - trS;
		result.gcv_scores[idx] = norm2 / (denom * denom);

		if (with_t_predictions) {
			auto &buf = result.t_predictions[idx];
			for (size_t v = 0; v < n_vertices; ++v)
				buf[v] = y_t[v];
		} else if (idx == 0) {
			// pre‐allocate predictions once
			result.predictions.resize(n_vertices);
		}
	}

	// 6) Pick optimal t
	result.opt_t_idx = std::min_element(
		result.gcv_scores.begin(),
		result.gcv_scores.end()
		) - result.gcv_scores.begin();

	// 7) Fill in result.predictions if not already stored
	if (with_t_predictions) {
		result.predictions = result.t_predictions[result.opt_t_idx];
	} else {
		double t = result.candidate_ts[result.opt_t_idx];
		// Use the same fix option as above
		Eigen::ArrayXd w = (-t * Eigen::Map<Eigen::ArrayXd>(result.evalues.data(), m).abs()).exp();
		Eigen::VectorXd y_t = result.evectors * (w * gft.array()).matrix();
		for (size_t v = 0; v < n_vertices; ++v)
			result.predictions[v] = y_t[v];
	}

	bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
	if (y_binary) {
		for (size_t i = 0; i < n_vertices; ++i) {
			result.predictions[i] = std::clamp(result.predictions[i], 0.0, 1.0);
		}

		if (with_t_predictions) {
			//result.t_predictions.assign(n_candidates, std::vector<double>(n_vertices));
			for (size_t j = 0; j < n_candidates; ++j) {
				for (size_t i = 0; i < n_vertices; ++i) {
					result.t_predictions[j][i] = std::clamp(result.t_predictions[j][i], 0.0, 1.0);
				}
			}
		}
	}

	return result;
}
