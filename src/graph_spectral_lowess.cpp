#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>              // For eigenvalue computation
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <execution>                    // For std::execution::seq/par
#include <thread>                       // For std::thread::hardware_concurrency
#include <algorithm>                    // For std::min_element
#include <numeric>                      // For std::accumulate
#include <map>                          // For std::map
#include <cmath>                        // For std::log2

// for debugging
//#include <filesystem>
//#include <fstream>
#include "cpp_utils.hpp"                // For debugging and elapsed.time
#include "exec_policy.hpp"
#include "graph_spectral_lowess.hpp"    // For graph_spectral_lowess_t
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
 * @brief Implements spectral-based locally weighted regression for graph data
 *
 * @details This function computes a local linear approximation of the response variable
 * at each vertex using spectral embedding of local neighborhoods. The algorithm:
 * 1. Computes the Laplacian eigenvectors for dimension reduction
 * 2. For each vertex, identifies neighboring vertices within varying bandwidths
 * 3. Creates spectral embeddings of these local neighborhoods
 * 4. Fits weighted linear models and selects optimal bandwidth based on LOOCV Rf_error
 * 5. Produces smoothed predictions, Rf_error estimates, and local scale information
 *
 * @param y Response values at each vertex in the graph
 * @param n_evectors Number of eigenvectors to use for the spectral embedding
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing for bandwidth grid; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances in kernel weight calculations
 * @param kernel_type Type of kernel function to use for weighting (e.g., Gaussian, triangular)
 * @param precision Precision threshold for binary search and numerical comparisons
 * @param verbose Whether to print progress information
 *
 * @return graph_spectral_lowess_t Structure containing:
 * - predictions: Smoothed values at each vertex
 * - errors: Estimated prediction errors
 * - scale: Local bandwidth/scale parameter for each vertex
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
	unsigned int available_threads = std::thread::hardware_concurrency();
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
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
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

	// Progress counter for parallel processing
	std::atomic<size_t> progress_counter{0};
	const size_t progress_step = std::max(size_t(1), n_vertices / 20);

	// Process each vertex
	std::vector<size_t> vertices(n_vertices);
	std::iota(vertices.begin(), vertices.end(), 0);

	// Use std::execution::par_unseq for parallel processing if available
	gflow::for_each(gflow::seq, vertices.begin(), vertices.end(),
				  [&](size_t vertex) {
					  try {
						  // Find minimum bandwidth that ensures enough vertices for modeling
						  double vertex_min_bw = find_minimum_radius_for_domain_min_size(
							  vertex,
							  min_bw,
							  max_bw,
							  domain_min_size,
							  precision
							  );

						  // Check bandwidth constraints
						  if (vertex_min_bw >= max_bw) {
							  REPORT_ERROR("Required minimum bandwidth (%.4f) exceeds maximum bandwidth (%.4f) for vertex %zu",
										   vertex_min_bw, max_bw, vertex);
						  }

						  // Generate candidate bandwidths
						  std::vector<double> candidate_bws = get_candidate_bws(
							  vertex_min_bw,
							  max_bw,
							  n_bws,
							  log_grid,
							  precision
							  );

						  // Find all vertices within maximum radius
						  auto reachable_vertices_map = find_vertices_within_radius(vertex, max_bw);

						  // Create storage for bandwidth evaluations
						  std::vector<double> bandwidth_errors(candidate_bws.size(), std::numeric_limits<double>::infinity());
						  std::vector<lm_t> bandwidth_models(candidate_bws.size());

						  // Evaluate each bandwidth
						  for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {
							  double current_bw = candidate_bws[bw_idx];

							  // Filter vertices by current bandwidth
							  std::map<size_t, double> local_vertex_map;
							  for (const auto& [v, dist] : reachable_vertices_map) {
								  if (dist <= current_bw) {
									  local_vertex_map[v] = dist;
								  }
							  }

							  // Skip if we don't have enough vertices
							  if (local_vertex_map.size() < domain_min_size) {
								  continue;
							  }

							  // Create embedding using eigenvectors
							  Eigen::MatrixXd embedding = create_spectral_embedding(
								  local_vertex_map,
								  eigenvectors,
								  n_evectors
								  );

							  // Fit weighted linear model
							  lm_t model = cleveland_fit_linear_model(
								  embedding,
								  y,
								  local_vertex_map,
								  dist_normalization_factor,
								  n_cleveland_iterations);

							  // Store model and Rf_error
							  bandwidth_errors[bw_idx] = model.mean_error;
							  bandwidth_models[bw_idx] = std::move(model);
						  }

						  // Find best bandwidth (minimum Rf_error)
						  std::vector<double>::iterator min_error_it;
						  bool smooth_bandwidth_errors = true;

						  if (smooth_bandwidth_errors) {
#if 1
							  double q_thld = 1.0 / 3.0;
							  int k = std::max(2, (int)(q_thld * n_bws));
							  double epsilon = 1e-10;
							  std::vector<double> null_vector;
							  auto bandwidth_errors_fit = uwmabilo(candidate_bws,
																   bandwidth_errors,
																   null_vector,
																   k,
																   k,
																   kernel_type,
																   dist_normalization_factor,
																   epsilon,
																   false);
							  bandwidth_errors = std::move(bandwidth_errors_fit.predictions);
#else

							  bandwidth_errors = running_window_average(bandwidth_errors, 2);
							  // print_vect(bandwidth_errors, "bandwidth_errors");
#endif

							  min_error_it = std::min_element(bandwidth_errors.begin(), bandwidth_errors.end());

						  } else {
							  min_error_it = std::min_element(bandwidth_errors.begin(), bandwidth_errors.end());
						  }

						  if (min_error_it != bandwidth_errors.end() && std::isfinite(*min_error_it)) {

							  // Get index of best bandwidth
							  size_t best_bw_idx = min_error_it - bandwidth_errors.begin();

							  // Store results for this vertex
							  const auto& best_model = bandwidth_models[best_bw_idx];

							  if (best_model.vertices.size() == 0) {
								  Rprintf("vertex: %zu\n", vertex);
								  Rprintf("best_bw_idx: %zu\n", best_bw_idx);
								  REPORT_ERROR("ERROR\n");
							  }

							  // Find this vertex's prediction in the model
							  auto it = std::find_if(best_model.vertices.begin(), best_model.vertices.end(),
													 [&](auto& v) { return v == vertex; });

							  // Store prediction, Rf_error, and scale
							  if (it != best_model.vertices.end()) {
								  size_t idx = it - best_model.vertices.begin();
								  result.predictions[vertex] = best_model.predictions[idx];
								  result.errors[vertex] = best_model.errors[idx];
							  } else {
								  // Fallback if vertex not found in predictions
								  Rprintf("vertex: %zu\n", vertex);
								  Rprintf("best_bw_idx: %zu\n", best_bw_idx);
								  Rprintf("best_model.vertices.size(): %zu\n", best_model.vertices.size());
								  //print_vect(best_model.vertices, "best_model.vertices");

								  REPORT_ERROR("ERROR: vertex %zu not found in best_model.vertices\n", vertex);
							  }

							  result.scale[vertex] = candidate_bws[best_bw_idx];

						  } else {
							  // No valid models found
							  result.predictions[vertex] = std::numeric_limits<double>::quiet_NaN();
							  result.errors[vertex] = std::numeric_limits<double>::infinity();
							  result.scale[vertex] = max_bw;
							  REPORT_ERROR("ERROR: No valid models found for vertex: %zu\n", vertex);
						  }

						  // Update progress counter
						  if (verbose) {
							  size_t current = ++progress_counter;
							  if (current % progress_step == 0 || current == n_vertices) {
								  double percentage = 100.0 * current / n_vertices;
								  Rprintf("\rProcessing vertices: %.1f%% complete (%zu/%zu)",
										  percentage, current, n_vertices);
								  R_FlushConsole();
							  }
						  }
					  } catch (const std::exception& e) {
						  REPORT_ERROR("Error processing vertex %zu: %s", vertex, e.what());
					  }
				  }
		);

	bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
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
