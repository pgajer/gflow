#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>              // For eigenvalue computation
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <algorithm>                    // For std::min_element
#include <numeric>                      // For std::accumulate
#include <map>                          // For std::map
#include <cmath>                        // For std::log2

// for debugging
//#include <filesystem>
//#include <fstream>
#include "omp_compat.h"
#include "cpp_utils.hpp"                // For debugging and elapsed.time
#include "graph_spectral_lowess.hpp"    // For graph_spectral_lowess_t
#include "bandwidth_utils.hpp"          // For get_candidate_bws
#include "kernels.h"                    // For kernel functions
#include "error_utils.h"                // For REPORT_ERROR
#include "mlm.hpp"                      // For lm_t structure
#include "set_wgraph.hpp"               // For set_wgraph_t
#include "cpp_stats_utils.hpp"          // For running_window_average()

Eigen::MatrixXd create_spectral_embedding(
	const std::map<size_t, double>& vertex_map,
	const Eigen::MatrixXd& eigenvectors,
	size_t n_evectors);

/**
 * @brief Spectral-based locally weighted regression on weighted graphs
 *
 * MOTIVATION
 *
 * Estimating smooth functions on graphs presents a fundamental challenge: how
 * do we adapt regression bandwidth to local geometric structure while accounting
 * for intrinsic graph topology? Classical LOWESS on Euclidean domains uses
 * ambient distance to define neighborhoods, but graphs possess intrinsic geometry
 * that may differ substantially from any embedding. This function addresses the
 * challenge by combining spectral methods for dimension reduction with local
 * weighted regression, allowing bandwidth adaptation while respecting graph structure.
 *
 * ALGORITHM OVERVIEW
 *
 * The algorithm proceeds in three main phases. First, we compute a global spectral
 * embedding using Laplacian eigenvectors. These eigenvectors provide coordinates
 * that encode graph structure through the diffusion geometry of random walks. The
 * eigenvectors corresponding to smallest eigenvalues capture large-scale structure,
 * while larger eigenvalues encode local variation.
 *
 * Second, for each vertex we identify candidate bandwidths spanning from a
 * minimum ensuring sufficient local sample size to a maximum encompassing
 * substantial portions of the graph. For each bandwidth, we extract the induced
 * neighborhood and project it into the spectral embedding space. This creates
 * a Euclidean domain where local regression can be performed while respecting
 * the graph's intrinsic geometry.
 *
 * Third, we fit weighted linear models in the spectral coordinates for each
 * bandwidth, using kernel weights that decay with graph distance. Cross-validation
 * error guides bandwidth selection, with optional smoothing of the error profile
 * to improve stability. The final prediction at each vertex comes from the
 * locally optimal model.
 *
 * SPECTRAL EMBEDDING CONSTRUCTION
 *
 * The Laplacian matrix L = D - A encodes the graph structure, where D is the
 * degree matrix and A the adjacency matrix. For weighted graphs, entries A[i,j]
 * represent edge lengths or costs. The Laplacian's eigendecomposition L = VΛV^T
 * provides coordinates: the i-th eigenvector V[:,i] assigns each vertex a value
 * reflecting its position in the graph's diffusion geometry.
 *
 * We compute the n_evectors smallest eigenvectors. The smallest (trivial)
 * eigenvector is constant, corresponding to eigenvalue zero. Subsequent
 * eigenvectors partition the graph at progressively finer scales. Using these
 * eigenvectors as coordinates provides an n_evectors-dimensional Euclidean
 * representation where Euclidean distance approximates diffusion distance on
 * the original graph.
 *
 * LOCAL NEIGHBORHOOD REGRESSION
 *
 * At each vertex i, we identify all vertices within graph distance bw, where
 * bw is a candidate bandwidth. The spectral embedding restricted to this
 * neighborhood gives coordinates X[j,1:n_evectors] for each neighbor j, and
 * we have response values y[j]. We fit the weighted linear model:
 *
 *   ŷ[j] = β₀ + β₁X[j,1] + ... + βₙ X[j,n_evectors]
 *
 * with weights w[j] = K(d[i,j] / h), where d[i,j] is graph distance from i
 * to j, h = dist_normalization_factor controls decay rate, and K is the kernel
 * function specified by kernel_type.
 *
 * The weighted least squares problem minimizes:
 *
 *   ∑ⱼ w[j](y[j] - ŷ[j])²
 *
 * Leave-one-out cross-validation estimates prediction error without requiring
 * a held-out test set. For each vertex j in the neighborhood, we refit the
 * model excluding j and predict ŷ₋ⱼ[j]. The LOOCV error averages squared
 * prediction errors across all vertices.
 *
 * BANDWIDTH SELECTION
 *
 * Bandwidth selection balances bias and variance. Small bandwidths yield
 * unbiased estimates but high variance due to limited data. Large bandwidths
 * reduce variance through extensive smoothing but introduce bias by averaging
 * over regions where the response varies substantially.
 *
 * We evaluate candidate bandwidths on a grid between vertex-specific minimums
 * (ensuring sufficient sample size) and a global maximum (capturing substantial
 * graph structure). For each bandwidth, we compute LOOCV error and select the
 * bandwidth minimizing this error.
 *
 * The algorithm optionally smooths the error profile across bandwidths using
 * a secondary smoothing procedure. This removes noise in error estimates that
 * can arise from sampling variability, particularly for vertices with irregular
 * neighborhoods. The smoothed profile often yields more stable bandwidth selection.
 *
 * PARAMETERS
 *
 * @param y Response values at graph vertices (length = number of vertices)
 * @param n_evectors Dimension of spectral embedding. Controls trade-off between
 *        computational cost and embedding fidelity. Typical values: 5-20.
 *        Larger values capture finer geometric detail but increase cost.
 *
 * @param n_bws Number of candidate bandwidths evaluated at each vertex.
 *        More candidates improve optimization but increase computation.
 *        Typical values: 10-30.
 *
 * @param log_grid If true, candidate bandwidths space logarithmically between
 *        minimum and maximum. Logarithmic spacing concentrates candidates at
 *        smaller bandwidths where error profiles often vary most rapidly.
 *        If false, use linear spacing.
 *
 * @param min_bw_factor Minimum bandwidth = min_bw_factor × graph_diameter.
 *        Controls smallest neighborhood size. Values 0.05-0.15 ensure sufficient
 *        local sample size while maintaining locality.
 *
 * @param max_bw_factor Maximum bandwidth = max_bw_factor × graph_diameter.
 *        Controls largest neighborhood size. Values 0.3-0.8 balance global
 *        information with local adaptation.
 *
 * @param dist_normalization_factor Scaling factor h for kernel weights.
 *        Larger values produce slower decay, effectively widening the kernel.
 *        Interacts with bandwidth: both control smoothing but at different scales.
 *
 * @param kernel_type Specifies kernel function K for distance-based weighting.
 *        Common choices: Gaussian, Epanechnikov, triangular, uniform.
 *        Kernel shape affects local influence structure.
 *
 * @param precision Numerical tolerance for distance comparisons and binary search
 *        in bandwidth determination. Typical value: 1e-6.
 *
 * @param n_cleveland_iterations Number of robustifying iterations for weighted
 *        least squares. Cleveland's original LOWESS uses 2-3 iterations to
 *        downweight outliers. Set to 1 for standard weighted regression.
 *
 * @param verbose If true, print progress information including eigendecomposition
 *        status and percentage of vertices processed.
 *
 * @return graph_spectral_lowess_t structure containing:
 *   - predictions: Smoothed response values at each vertex
 *   - errors: Estimated prediction error (from LOOCV) at each vertex
 *   - scale: Optimal bandwidth selected for each vertex
 *
 * LIMITATIONS AND EXTENSIONS
 *
 * The current implementation uses a global Laplacian eigendecomposition.
 * This means eigenvectors encode global graph structure, which may not
 * optimally represent local geometry in heterogeneous graphs. A natural
 * extension computes neighborhood-specific Laplacians for each vertex,
 * providing spectral coordinates adapted to local structure.
 *
 * The Laplacian construction uses edge weights from the adjacency matrix,
 * which represent ambient geometry. An alternative constructs a Riemannian
 * Laplacian from a density-derived metric that accounts for data distribution
 * and response coherence. This would provide intrinsic geometry more suitable
 * for statistical inference.
 *
 * @see set_wgraph_t For the weighted graph structure and distance computations
 * @see cleveland_fit_linear_model For the weighted regression implementation
 * @see get_candidate_bws For bandwidth grid generation
 */
graph_spectral_lowess_t set_wgraph_t::graph_spectral_lowess(
	const std::vector<double>& y,
	size_t n_evectors,
	// bw parameters
	size_t n_bws,
	bool log_grid,
	double min_bw_factor,
	double max_bw_factor,
	// kernel parameters
	double dist_normalization_factor,
	size_t kernel_type,
	// other
	double precision,
	size_t n_cleveland_iterations,
	bool verbose
	) const {

#define DEBUG__graph_spectral_lowess 0

	auto start_time = std::chrono::steady_clock::now();

	// Initialize kernel function
	initialize_kernel(kernel_type, 1.0);

	// Set Eigen to use available threads for parallel computation
	unsigned int available_threads = gflow_get_max_threads();
	if (available_threads == 0) available_threads = 4;  // Fallback if detection fails
	Eigen::setNbThreads(available_threads);

	// Define minimum and maximum bandwidth based on graph diameter
	if (graph_diameter <= 0) {
		REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
	}

	double min_bw = min_bw_factor * graph_diameter;
	double max_bw = max_bw_factor * graph_diameter;

	// Minimum number of vertices needed for a robust n_evectors-dimensional model
	size_t domain_min_size = n_evectors + 5;  // Adding extra points for stability

	if (verbose) {
		Rprintf("Starting graph_spectral_lowess() algorithm\n");
		Rprintf("Number of vertices: %zu\n", adjacency_list.size());
		Rprintf("Graph diameter: %f\n", graph_diameter);
		Rprintf("min_bw: %.4f\n", min_bw);
		Rprintf("max_bw: %.4f\n", max_bw);

		Rprintf("\nNumber of eigenvectors: %zu\n", n_evectors);
		Rprintf("domain_min_size: %zu\n\n", domain_min_size);
		Rprintf("Using %u threads for Eigen operations\n", available_threads);
		Rprintf("\n");
	}

	// Initialize result structure
	size_t n_vertices = adjacency_list.size();
	graph_spectral_lowess_t result;
	result.predictions.resize(n_vertices);
	result.errors.resize(n_vertices);
	result.scale.resize(n_vertices);

	// Construct the Laplacian matrix
	if (verbose) {
		Rprintf("Constructing Laplacian matrix...\n");
	}

	// Create adjacency matrix as a sparse matrix
	Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
	std::vector<Eigen::Triplet<double>> triples;
	triples.reserve(adjacency_list.size() * 2);  // Estimate for undirected graph

	// Use edge weights for Laplacian construction
	//bool use_edge_lengths = true;  // Default to using actual edge lengths
	for (size_t i = 0; i < n_vertices; ++i) {
		for (const auto& edge : adjacency_list[i]) {
			size_t j = edge.vertex;
			//double weight = use_edge_lengths ? edge.weight : 1.0;
			if (i < j) {   // Ensure each edge is added only once
				triples.push_back(Eigen::Triplet<double>(i, j, edge.weight));
				triples.push_back(Eigen::Triplet<double>(j, i, edge.weight));
			}
		}
	}
	A.setFromTriplets(triples.begin(), triples.end());

	// Compute degree matrix
	Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
	for (int k = 0; k < A.outerSize(); ++k) {
		double sum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
			sum += it.value();
		}
		D.insert(k, k) = sum;
	}

	// Compute Laplacian: L = D - A
	Eigen::SparseMatrix<double> L = D - A;

	// Add small regularization for numerical stability
	for (size_t i = 0; i < n_vertices; i++) {
		L.coeffRef(i, i) += 1e-8;
	}

	// Compute eigenvalues and eigenvectors of the Laplacian
	if (verbose) {
		Rprintf("Computing Laplacian eigenvectors...\n");
	}

	Spectra::SparseSymMatProd<double> op(L);

	// Ensure we request enough eigenvectors
	int nev = std::min(n_evectors + 5, n_vertices);  // Request a few extra
	int ncv = std::min(2 * nev, (int)n_vertices);  // Control parameter for the algorithm

	// Ensure nev < ncv
	if (nev >= ncv) {
		ncv = nev + 1;
	}

	// Construct eigen solver to find eigenvalues closest to 0
	// ncv rule-of-thumb: 3*nev is often safer for hard problems
	int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
	ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n

	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv_default);
	eigs.init();
	int maxit = 1000;  // Increase maximum iterations
	double tol = 1e-10;  // Tolerance for convergence
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);
	Eigen::MatrixXd eigenvectors;

	if (eigs.info() != Spectra::CompInfo::Successful) {
		// if (verbose) {
		// 	Rprintf("Initial eigenvalue computation failed. Attempting with adjusted parameters...\n");
		// }
		// Define fallback parameters to try
		std::vector<std::pair<int, double>> attempts = {
			{2000, 1e-8},    // More iterations, slightly relaxed tolerance
			{3000, 1e-6},    // Even more iterations, more relaxed tolerance
			{5000, 1e-4} // Final attempt with very relaxed parameters
		};
		// Try with original ncv first
		bool success = false;
		for (const auto& params : attempts) {
			int adjusted_maxit = params.first;
			double adjusted_tol = params.second;
			eigs.init();
			eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
			if (eigs.info() == Spectra::CompInfo::Successful) {
				if (verbose) {
					Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
							ncv, adjusted_maxit, adjusted_tol);
				}
				eigenvectors = eigs.eigenvectors();  // Add this line to extract eigenvectors
				success = true;
				break;
			}
		}
		// If still not successful, try with increased ncv values
		if (!success) {
			int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
			std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

			for (const auto& multiplier : ncv_multipliers) {
				int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
				// Create a new solver with adjusted ncv
				Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
				for (const auto& params : attempts) {
					int adjusted_maxit = params.first;
					double adjusted_tol = params.second;
					adjusted_eigs.init();
					adjusted_eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
					if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
						if (verbose) {
							Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
									adjusted_ncv, adjusted_maxit, adjusted_tol);
						}
						eigenvectors = adjusted_eigs.eigenvectors();
						success = true;
						break;
					}
				}
				if (success) {
					break;
				}
			}
		}
		// If all attempts failed, report an Rf_error
		if (!success) {
			REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
		}
	} else {
		// Get eigenvectors
		eigenvectors = eigs.eigenvectors();
	}

	// Progress counter for serial processing
	size_t progress_counter = 0;
	const size_t progress_step = std::max<size_t>(1, n_vertices / 20);

	// Process each vertex (serial)
	for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
		try {
			// Find minimum bandwidth that ensures enough vertices for modeling
			double vertex_min_bw = find_minimum_radius_for_domain_min_size(
				vertex, min_bw, max_bw, domain_min_size, precision
				);

			// Check bandwidth constraints
			if (vertex_min_bw >= max_bw) {
				REPORT_ERROR("Required minimum bandwidth (%.4f) exceeds maximum bandwidth (%.4f) for vertex %zu",
							 vertex_min_bw, max_bw, vertex);
			}

			// Generate candidate bandwidths
			std::vector<double> candidate_bws = get_candidate_bws(
				vertex_min_bw, max_bw, n_bws, log_grid, precision
				);

			// Find all vertices within maximum radius
			auto reachable_vertices_map = find_vertices_within_radius(vertex, max_bw);

			// Create storage for bandwidth evaluations
			std::vector<double> bandwidth_errors(candidate_bws.size(),
												 std::numeric_limits<double>::infinity());
			std::vector<lm_t> bandwidth_models(candidate_bws.size());

			// Evaluate each bandwidth
			for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {
				double current_bw = candidate_bws[bw_idx];

				// Filter vertices by current bandwidth
				std::map<size_t, double> local_vertex_map;
				for (const auto& [v, dist] : reachable_vertices_map) {
					if (dist <= current_bw) local_vertex_map[v] = dist;
				}

				// Skip if we don't have enough vertices
				if (local_vertex_map.size() < domain_min_size) continue;

				// Create embedding using eigenvectors
				Eigen::MatrixXd embedding = create_spectral_embedding(
					local_vertex_map, eigenvectors, n_evectors
					);

				// Fit weighted linear model
				lm_t model = cleveland_fit_linear_model(
					embedding, y, local_vertex_map,
					dist_normalization_factor, n_cleveland_iterations
					);

				// Store model and error
				bandwidth_errors[bw_idx] = model.mean_error;
				bandwidth_models[bw_idx] = std::move(model);
			}

			// Find best bandwidth (minimum error)
			std::vector<double>::iterator min_error_it;
			bool smooth_bandwidth_errors = true;

			if (smooth_bandwidth_errors) {
				double q_thld = 1.0 / 3.0;
				int k = std::max(2, static_cast<int>(q_thld * candidate_bws.size()));

				std::vector<double> finite_filled_errors = bandwidth_errors;
				double max_finite_error = -std::numeric_limits<double>::infinity();
				for (double e : finite_filled_errors) {
					if (std::isfinite(e) && e > max_finite_error) {
						max_finite_error = e;
					}
				}
				if (std::isfinite(max_finite_error)) {
					for (double& e : finite_filled_errors) {
						if (!std::isfinite(e)) {
							e = max_finite_error;
						}
					}
					bandwidth_errors = running_window_average(finite_filled_errors, k);
				}
			}
			min_error_it = std::min_element(bandwidth_errors.begin(), bandwidth_errors.end());

			if (min_error_it != bandwidth_errors.end() && std::isfinite(*min_error_it)) {
				// Index of best bandwidth
				size_t best_bw_idx = static_cast<size_t>(min_error_it - bandwidth_errors.begin());

				// Store results for this vertex
				const auto& best_model = bandwidth_models[best_bw_idx];

				if (best_model.vertices.empty()) {
					Rprintf("vertex: %zu\n", vertex);
					Rprintf("best_bw_idx: %zu\n", best_bw_idx);
					REPORT_ERROR("ERROR\n");
				}

				// Find this vertex's prediction in the model
				auto it = std::find(best_model.vertices.begin(),
									best_model.vertices.end(), vertex);

				// Store prediction, error, and scale
				if (it != best_model.vertices.end()) {
					size_t idx = static_cast<size_t>(it - best_model.vertices.begin());
					result.predictions[vertex] = best_model.predictions[idx];
					result.errors[vertex]      = best_model.errors[idx];
				} else {
					Rprintf("vertex: %zu\n", vertex);
					Rprintf("best_bw_idx: %zu\n", best_bw_idx);
					Rprintf("best_model.vertices.size(): %zu\n", best_model.vertices.size());
					REPORT_ERROR("ERROR: vertex %zu not found in best_model.vertices\n", vertex);
				}

				result.scale[vertex] = candidate_bws[best_bw_idx];

			} else {
				// No valid models found
				result.predictions[vertex] = std::numeric_limits<double>::quiet_NaN();
				result.errors[vertex]      = std::numeric_limits<double>::infinity();
				result.scale[vertex]       = max_bw;
				REPORT_ERROR("ERROR: No valid models found for vertex: %zu\n", vertex);
			}

			// Progress
			if (verbose) {
				++progress_counter;
				if (progress_counter % progress_step == 0 || progress_counter == n_vertices) {
					double percentage = 100.0 * progress_counter / static_cast<double>(n_vertices);
					Rprintf("\rProcessing vertices: %.1f%% complete (%zu/%zu)",
							percentage, progress_counter, n_vertices);
					R_FlushConsole();
				}
			}
		} catch (const std::exception& e) {
			REPORT_ERROR("Error processing vertex %zu: %s", vertex, e.what());
		}
	}

	auto is_binary01 = [](const std::vector<double>& yy, double tol = 1e-12) -> bool {
        for (double v : yy) {
            if (!(std::fabs(v) <= tol || std::fabs(v - 1.0) <= tol)) {
                return false;
            }
        }
        return true;
    };

    const bool y_binary = is_binary01(y);

	if (y_binary) {
		for (auto& p : result.predictions) {
			p = std::clamp(p, 0.0, 1.0);
		}
	}

	if (verbose) {
		Rprintf("\nCompleted graph_spectral_lowess in %.2f seconds\n",
				std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count());
	}

	return result;
}
