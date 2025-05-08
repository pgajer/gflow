// Core C++ STL
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <execution>

// I/O
#include <cstdio>
#include <iostream>

// Eigen for linear algebra
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// Spectra for eigen decomposition
#include <SymEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseSymMatProd.h>

// Project-specific headers
#include "graph_bw_adaptive_spectral_smoother.hpp" // For graph_bw_adaptive_spectral_smoother_t
#include "bandwidth_utils.hpp"          // For get_candidate_bws
#include "kernels.h"                    // For kernel functions
#include "error_utils.h"                // For REPORT_ERROR
#include "mlm.hpp"                      // For lm_t structure
#include "set_wgraph.hpp"               // For set_wgraph_t
#include "mabilo.hpp"                   // For uwmabilo()
#include "cpp_stats_utils.hpp"          // For running_window_average()

Eigen::MatrixXd create_spectral_embedding(
	const std::map<size_t, double>& vertex_map,
	const Eigen::MatrixXd& eigenvectors,
	size_t n_evectors);

/**
 * @brief Spectral-based graph smoother with global bandwidth selection and optional diagnostics
 *
 * @details This function estimates the conditional expectation \f$ \mathbb{E}[Y \mid G] \f$ over a graph
 * using locally weighted linear regression models in a spectral embedding space.
 * It follows a global bandwidth selection strategy by evaluating prediction error over a shared
 * set of bandwidths across all vertices and selecting the bandwidth index with the lowest average error.
 *
 * Each vertex is embedded into a lower-dimensional Euclidean space using Laplacian eigenvectors, and
 * kernel-weighted linear models are fit over local neighborhoods defined by increasing bandwidths.
 *
 * @param y Vector of response values, one per vertex
 * @param n_evectors Number of Laplacian eigenvectors used for spectral embedding
 * @param min_bw_factor Minimum bandwidth factor (relative to graph diameter)
 * @param max_bw_factor Maximum bandwidth factor (relative to graph diameter)
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, bandwidth grid is log-spaced; otherwise, linearly spaced
 * @param kernel_type Type of kernel function used for computing weights
 * @param dist_normalization_factor Factor used to scale maximum neighborhood distances
 * @param precision Numerical precision used in binary search and grid generation
 * @param use_global_bw_grid If true, use global candidate bandwidths schema
 * @param with_bw_predictions If true, returns predictions for all bandwidths (bw_predictions)
 * @param with_vertex_bw_errors If true, returns LOOCV error per vertex and bandwidth index
 * @param verbose If true, prints progress and debugging information
 *
 * @return graph_spectral_smoother_t structure containing:
 * - predictions: Final prediction vector using the globally optimal bandwidth index
 * - bw_predictions: Matrix of predictions for all bandwidths (only if with_bw_predictions = true)
 * - vertex_bw_errors: Matrix of LOOCV errors per vertex and bandwidth (only if with_vertex_bw_errors = true)
 * - bw_mean_abs_errors: Mean absolute error per bandwidth across all vertices
 * - vertex_min_bws: Minimum usable bandwidth per vertex
 * - opt_bw_idx: Index of the globally selected optimal bandwidth
 */
graph_bw_adaptive_spectral_smoother_t set_wgraph_t::graph_bw_adaptive_spectral_smoother(
	const std::vector<double>& y,
	size_t n_evectors,
	double min_bw_factor,
	double max_bw_factor,
	size_t n_bws,
	bool log_grid,
	size_t kernel_type,
	double dist_normalization_factor,
	double precision,
	bool use_global_bw_grid,
	bool with_bw_predictions,
	bool with_vertex_bw_errors,
	bool verbose
	) {

	size_t n_vertices = adjacency_list.size();
	graph_bw_adaptive_spectral_smoother_t result;
	result.predictions.resize(n_vertices);
	result.vertex_min_bws.resize(n_vertices);
	result.bw_mean_abs_errors.resize(n_bws, 0.0);
	result.opt_bw_idx = 0;
	result.bw_predictions.resize(n_bws, std::vector<double>(n_vertices, 0.0));

	std::vector<std::vector<double>> vertex_errors;
	if (with_vertex_bw_errors)
		vertex_errors.resize(n_vertices, std::vector<double>(n_bws, std::numeric_limits<double>::infinity()));

	Eigen::MatrixXd eigenvectors = compute_graph_laplacian_eigenvectors(n_evectors, verbose);

	compute_graph_diameter();
	double min_bw = min_bw_factor * graph_diameter;
	double max_bw = max_bw_factor * graph_diameter;

	std::vector<double> candidate_bws;
	if (use_global_bw_grid) {
		candidate_bws = get_candidate_bws(
			min_bw,
			max_bw,
			n_bws,
			log_grid,
			precision);
	}

	initialize_kernel(kernel_type, 1.0);

	// Set Eigen to use available threads for parallel computation
	unsigned int available_threads = std::thread::hardware_concurrency();
	if (available_threads == 0) available_threads = 4;  // Fallback if detection fails
	Eigen::setNbThreads(available_threads);

	// Minimum number of vertices needed for a robust n_evectors-dimensional model
	size_t domain_min_size = n_evectors + 5;  // Adding extra points for stability

	if (verbose) {
		Rprintf("Starting graph_bw_adaptive_spectral_smoother()\n");
		Rprintf("Number of vertices: %zu\n", n_vertices);
		Rprintf("Graph diameter: %f\n", graph_diameter);
		Rprintf("min_bw: %.4f\n", min_bw);
		Rprintf("max_bw: %.4f\n", max_bw);

		Rprintf("\nNumber of eigenvectors: %zu\n", n_evectors);
		Rprintf("domain_min_size: %zu\n\n", domain_min_size);
		Rprintf("Using %u threads for Eigen operations\n", available_threads);
		Rprintf("\n");
	}

	for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

		double vertex_min_bw = find_minimum_radius_for_domain_min_size(
			vertex,
			min_bw,
			max_bw,
			domain_min_size,
			precision
			);

		std::vector<double> vertex_candidate_bws;
		if (!use_global_bw_grid) {
			if (vertex_min_bw < min_bw) vertex_min_bw = min_bw;

			result.vertex_min_bws[vertex] = vertex_min_bw;

			vertex_candidate_bws = get_candidate_bws(
				vertex_min_bw, max_bw, n_bws, log_grid, precision);

		} else {
			result.vertex_min_bws[vertex] = min_bw;

			vertex_candidate_bws = candidate_bws;
			for (size_t bw_idx = 0; bw_idx < vertex_candidate_bws.size(); ++bw_idx) {
				if (vertex_candidate_bws[bw_idx] < vertex_min_bw) {
					vertex_candidate_bws[bw_idx] = vertex_min_bw;
				}
			}
		}

		auto local_map = find_vertices_within_radius(vertex, max_bw);

		for (size_t bw_idx = 0; bw_idx < vertex_candidate_bws.size(); ++bw_idx) {

			double current_bw = vertex_candidate_bws[bw_idx];

			std::map<size_t, double> subdomain;
			for (const auto& [v, dist] : local_map) {
				if (dist <= current_bw)
					subdomain[v] = dist;
			}

			if (subdomain.size() < n_evectors + 5)
				continue;

			Eigen::MatrixXd embedding = create_spectral_embedding(
				subdomain,
				eigenvectors,
				n_evectors);

			lm_t model = cleveland_fit_linear_model(
				embedding, y, subdomain, dist_normalization_factor);

			auto it = std::find(model.vertices.begin(), model.vertices.end(), vertex);
			if (it == model.vertices.end()) continue;
			size_t idx = std::distance(model.vertices.begin(), it);
			double prediction = model.predictions[idx];
			double error = model.errors[idx];

			result.bw_predictions[bw_idx][vertex] = prediction;

			if (with_vertex_bw_errors)
				vertex_errors[vertex][bw_idx] = error;

			result.bw_mean_abs_errors[bw_idx] += std::abs(prediction - y[vertex]);
		}
	}

	for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
		result.bw_mean_abs_errors[bw_idx] /= n_vertices;
	}

	result.opt_bw_idx = std::min_element(result.bw_mean_abs_errors.begin(),
										 result.bw_mean_abs_errors.end()) - result.bw_mean_abs_errors.begin();

	if (!with_bw_predictions) {
		result.predictions = std::move(result.bw_predictions[result.opt_bw_idx]);

		// Clear all other bandwidth predictions to free memory
		for (size_t i = 0; i < n_bws; ++i) {
			if (i != result.opt_bw_idx) {
				std::vector<double>().swap(result.bw_predictions[i]);  // Swap idiom for capacity release
			}
		}

		// Optionally shrink the bw_predictions outer vector
		// but only if you want to aggressively trim to just one element:
		// result.bw_predictions.resize(1); // <--- only if you're okay changing its shape
	} else {
		result.predictions = result.bw_predictions[result.opt_bw_idx];
	}


	return result;
}


