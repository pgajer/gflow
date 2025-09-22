#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>             // For eigenvalue computation
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <algorithm>                   // For std::min_element
#include <numeric>                     // For std::accumulate
#include <map>                         // For std::map
#include <cmath>                       // For std::log2

// for debugging
// #include <filesystem>
// #include <fstream>
#include "cpp_utils.hpp"               // For debugging and elapsed.time
#include "omp_compat.h"
#include "exec_policy.hpp"
#include "nada_graph_spectral_lowess.hpp" // For nada_graph_spectral_lowess_t
#include "bandwidth_utils.hpp"         // For get_candidate_bws
#include "kernels.h"                   // For kernel functions
#include "error_utils.h"               // For REPORT_ERROR
#include "mlm.hpp"                     // For lm_t structure
#include "set_wgraph.hpp"              // For set_wgraph_t
#include "mabilo.hpp"                  // For uwmabilo()
#include "cpp_stats_utils.hpp"         // For running_window_average()

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
 *         - predictions: Smoothed values at each vertex
 *         - errors: Estimated prediction errors
 *         - scale: Local bandwidth/scale parameter for each vertex
 */
nada_graph_spectral_lowess_t set_wgraph_t::nada_graph_spectral_lowess(
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


// #define DEBUG__nada_graph_spectral_lowess 0

    auto start_time = std::chrono::steady_clock::now();

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    // Set Eigen to use available threads for parallel computation
    unsigned int available_threads = gflow_get_max_threads();
    if (available_threads == 0) available_threads = 4; // Fallback if detection fails
    Eigen::setNbThreads(available_threads);

    // Define minimum and maximum bandwidth based on graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
    }

	size_t n_vertices = adjacency_list.size();

	double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    // Minimum number of vertices needed for a robust n_evectors-dimensional model
    size_t domain_min_size = n_evectors + 5; // Adding extra points for stability

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

	#if 0
	// gather vertex-wise bw_min info and write to file
	{
		// debugging/testing - saving to file
		std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data/";
		if (!std::filesystem::exists(debug_dir)) {
			if (!std::filesystem::create_directories(debug_dir)) {
				REPORT_ERROR("ERROR: Failed to create debug directory: %s\n", debug_dir.c_str());
			}
		}
		std::string bw_mins_file_name = debug_dir + "bw_mins.csv";
		std::ofstream bw_mins_file(bw_mins_file_name);
		if (!bw_mins_file.is_open()) {
			REPORT_ERROR("ERROR: Failed to open file for writing bw mins: %s\n",
						 bw_mins_file_name.c_str());
		}

		// Process each vertex
		std::vector<size_t> vertices(n_vertices);
		std::iota(vertices.begin(), vertices.end(), 0);

		gflow::for_each(gflow::seq, vertices.begin(), vertices.end(),
					  [&](size_t vertex) {

						  // Find minimum bandwidth that ensures enough vertices for modeling
						  // Ensure we have enough vertices for the model
						  double vertex_min_bw = find_minimum_radius_for_domain_min_size(
							  vertex,
							  min_bw,
							  max_bw,
							  domain_min_size,
							  precision
							  );

						  bw_mins_file << vertex_min_bw << "\n";
					  }
			);

		bw_mins_file.close();
		Rprintf("\nbw mins were written to file: %s\n", bw_mins_file_name.c_str());
	}
	#endif

    // Initialize result structure
    nada_graph_spectral_lowess_t result;
    result.predictions.resize(n_vertices);
    result.errors.resize(n_vertices);

    // Construct the Laplacian matrix
    if (verbose) {
        Rprintf("Constructing Laplacian matrix...\n");
    }

    // Create adjacency matrix as a sparse matrix
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;
    triples.reserve(adjacency_list.size() * 2); // Estimate for undirected graph

    // Use edge weights for Laplacian construction
    //bool use_edge_lengths = true; // Default to using actual edge lengths
    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            //double weight = use_edge_lengths ? edge.weight : 1.0;
            if (i < j) {  // Ensure each edge is added only once
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
    int nev = std::min(n_evectors + 5, n_vertices); // Request a few extra
    int ncv = std::min(2 * nev, (int)n_vertices); // Control parameter for the algorithm

    // Ensure nev < ncv
    if (nev >= ncv) {
        ncv = nev + 1;
    }

	// Construct eigen solver to find eigenvalues closest to 0
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
	eigs.init();
	int maxit = 1000; // Increase maximum iterations
	double tol = 1e-10; // Tolerance for convergence
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);
	Eigen::MatrixXd eigenvectors;

	if (eigs.info() != Spectra::CompInfo::Successful) {
		// if (verbose) {
		// 	Rprintf("Initial eigenvalue computation failed. Attempting with adjusted parameters...\n");
		// }
		// Define fallback parameters to try
		std::vector<std::pair<int, double>> attempts = {
			{2000, 1e-8},   // More iterations, slightly relaxed tolerance
			{3000, 1e-6},   // Even more iterations, more relaxed tolerance
			{5000, 1e-4}    // Final attempt with very relaxed parameters
		};
		// Try with original ncv first
		bool success = false;
		for (const auto& params : attempts) {
			int adjusted_maxit = params.first;
			double adjusted_tol = params.second;
			// if (verbose) {
			// 	Rprintf("Trying with original ncv=%d, maxit=%d, tol=%g\n",
			// 			ncv, adjusted_maxit, adjusted_tol);
			// }
			eigs.init();
			eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
			if (eigs.info() == Spectra::CompInfo::Successful) {
				if (verbose) {
					Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
							ncv, adjusted_maxit, adjusted_tol);
				}
				eigenvectors = eigs.eigenvectors(); // Add this line to extract eigenvectors
				success = true;
				break;
			}
		}
		// If still not successful, try with increased ncv values
		if (!success) {

			// Define multipliers for ncv
			// std::vector<int> ncv_multipliers = {2, 4, 8};
			// int max_ncv = std::min(16 * nev, (int)n_vertices);

			int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
			std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

			for (const auto& multiplier : ncv_multipliers) {
				int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
				// if (verbose) {
				// 	Rprintf("Trying with increased ncv=%d\n", adjusted_ncv);
				// }
				// Create a new solver with adjusted ncv
				Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
				for (const auto& params : attempts) {
					int adjusted_maxit = params.first;
					double adjusted_tol = params.second;
					// if (verbose) {
					// 	Rprintf("Trying with ncv=%d, maxit=%d, tol=%g\n",
					// 			adjusted_ncv, adjusted_maxit, adjusted_tol);
					// }
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

	// 1. Collect minimum bandwidths for all vertices
	std::vector<double> vertex_min_bws(n_vertices);
	for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
		vertex_min_bws[vertex] = find_minimum_radius_for_domain_min_size(
			vertex, min_bw, max_bw, domain_min_size, precision);
	}

	// 2. Compute the median minimum bandwidth
	auto compute_median = [](std::vector<double> values) -> double {
		size_t n = values.size();
		if (n == 0) return 0.0;
		std::sort(values.begin(), values.end());
		if (n % 2 == 0) {
			return (values[n/2 - 1] + values[n/2]) / 2.0;
		}
		return values[n/2];
	};
	double median_min_bw = compute_median(vertex_min_bws);

	// 3. Generate global bandwidth grid
	std::vector<double> global_bws = get_candidate_bws(
		median_min_bw, max_bw, n_bws, log_grid, precision);

	// 4. Create storage for bandwidth evaluations across all vertices
	std::vector<std::vector<lm_t>> all_bandwidth_models(global_bws.size(),
														std::vector<lm_t>(n_vertices));
	std::vector<double> mean_bandwidth_errors(global_bws.size(), 0.0);

	// 5. Process each bandwidth
	for (size_t bw_idx = 0; bw_idx < global_bws.size(); ++bw_idx) {
		double current_bw = global_bws[bw_idx];

		// Process each vertex with this bandwidth
		for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

			double vertex_current_bw = (current_bw < vertex_min_bws[vertex]) ? vertex_min_bws[vertex] : current_bw;

			// Find vertices within this bandwidth
			auto unordered_local_vertex_map = find_vertices_within_radius(vertex, vertex_current_bw);
			std::map<size_t, double> local_vertex_map(unordered_local_vertex_map.begin(), unordered_local_vertex_map.end());

			if (local_vertex_map.size() < domain_min_size) {
				REPORT_ERROR("local_vertex_map.size() < domain_min_size");
			}

			// Create embedding using eigenvectors
			Eigen::MatrixXd embedding = create_spectral_embedding(
				local_vertex_map,
				eigenvectors,
				n_evectors
				);

			lm_t model = cleveland_fit_linear_model(
				embedding,
				y,
				local_vertex_map,
				dist_normalization_factor,
				n_cleveland_iterations);

			// Store model and Rf_error
			mean_bandwidth_errors[bw_idx] += model.mean_error;
			all_bandwidth_models[bw_idx][vertex] = std::move(model);
		}

		mean_bandwidth_errors[bw_idx] /= n_vertices;
	}

	// 6. Find bandwidth with the smallest mean prediction Rf_error
	auto min_error_it = std::min_element(mean_bandwidth_errors.begin(), mean_bandwidth_errors.end());

	if (!std::isfinite(*min_error_it)) {
		REPORT_ERROR("ERROR: min_error_it is not finite\n");
	}

	size_t best_bw_idx;
	if (min_error_it != mean_bandwidth_errors.end()) {
		// Get index of best bandwidth
		best_bw_idx = min_error_it - mean_bandwidth_errors.begin();

		for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

			// Store results for this vertex
			const auto& best_model = all_bandwidth_models[best_bw_idx][vertex];

			if (best_model.vertices.size() == 0) {
				Rprintf("vertex: %zu\n", vertex);
				Rprintf("best_bw_idx: %zu\n", best_bw_idx);
				REPORT_ERROR("ERROR: best_model.vertices.size(): %zu\n", best_model.vertices.size());
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
				Rprintf("vertex: %zu\n", vertex);
				Rprintf("best_bw_idx: %zu\n", best_bw_idx);
				Rprintf("best_model.vertices.size(): %zu\n", best_model.vertices.size());
				//print_vect(best_model.vertices, "best_model.vertices");
				REPORT_ERROR("ERROR: vertex %zu not found in best_model.vertices\n", vertex);
			}
		}

	} else {
		REPORT_ERROR("ERROR: No minimal mean_bandwidth_error was found\n");
	}

	// populating bw_predictions
	std::vector<std::vector<double>> bw_predictions(global_bws.size(),
													std::vector<double>(n_vertices));

	for (size_t bw_idx = 0; bw_idx < global_bws.size(); ++bw_idx) {

		for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

			// Store results for this vertex
			const auto& model = all_bandwidth_models[bw_idx][vertex];

			if (model.vertices.size() == 0) {
				Rprintf("vertex: %zu\n", vertex);
				Rprintf("bw_idx: %zu\n", bw_idx);
				REPORT_ERROR("ERROR: model.vertices.size(): %zu\n", model.vertices.size());
			}

			// Find this vertex's prediction in the model
			auto it = std::find_if(model.vertices.begin(), model.vertices.end(),
								   [&](auto& v) { return v == vertex; });

			// Store prediction, Rf_error, and scale
			if (it != model.vertices.end()) {
				size_t idx = it - model.vertices.begin();
				bw_predictions[bw_idx][vertex] = model.predictions[idx];

			} else {
				Rprintf("vertex: %zu\n", vertex);
				Rprintf("bw_idx: %zu\n", bw_idx);
				Rprintf("model.vertices.size(): %zu\n", model.vertices.size());
				//print_vect(model.vertices, "model.vertices");
				REPORT_ERROR("ERROR: vertex %zu not found in model.vertices\n", vertex);
			}
		}
	}

    if (verbose) {
        Rprintf("\nCompleted graph_spectral_lowess in %.2f seconds\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count());
    }

	result.opt_bw_idx     = best_bw_idx;
	result.opt_bw         = global_bws[best_bw_idx];
	result.global_bws     = std::move(global_bws);
	result.bw_predictions = std::move(bw_predictions);

    return result;
}
