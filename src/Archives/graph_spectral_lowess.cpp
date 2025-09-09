#include <Spectra/SymEigsSolver.h>     // For eigenvalue computation
#include <execution>                   // For std::execution::seq/par
#include <thread>                      // For std::thread::hardware_concurrency
#include <algorithm>                   // For std::min_element
#include <numeric>                     // For std::accumulate

#include "bandwidth_utils.hpp"         // For get_candidate_bws
#include "kernels.h"                   // For kernel functions
#include "error_utils.h"               // For REPORT_ERROR
#include "lm.hpp"                      // For lm_t structure

lm_t fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    size_t kernel_type,
    double dist_normalization_factor);

Eigen::MatrixXd create_spectral_embedding(
    const std::unordered_map<size_t, double>& vertices,
    const Eigen::MatrixXd& eigenvectors,
    size_t n_evectors);

struct graph_spectral_lowess_result_t {
	std::vector<double> predictions;      ///< predictions[i] is an estimate of E(Y|G) at the i-th vertex
    std::vector<double> errors;           ///< errors[i] is an estimate of the prediction error at the i-th vertex
    std::vector<double> scale;            ///< scale[i] is a local scale (radius/bandwidth) at the i-th vertex
};

/**
 *
 * @param y Response values at each vertex in the original graph
 * @param n_evectors Number of eigenvectors to use for the construction of a local embedding of the set of reachable vertices around each vertex
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances in the graph
 * @param kernel_type Type of kernel function used for weight calculation (e.g., Gaussian, triangular)
 * @param verbose Whether to print progress information
 */
graph_spectral_lowess_result_t set_wgraph_t::graph_spectral_lowess(
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
	bool verbose
	) const {

	if (verbose) {
        Rprintf("Starting graph_spectral_lowess()\n");
        Rprintf("n_vertices(grid_graph): %zu\n", adjacency_list.size());
		Rprintf("graph_diameter: %f\n", graph_diameter);
	}

	initialize_kernel(kernel_type, 1.0);

	// Set Eigen to use 12 threads for parallel computation
	Eigen::setNbThreads(12);

	// Determined graph diameter
	auto grid_graph = uniform_grid_graph_t(adjacency_list);
	auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
	auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
	graph_diameter = diameter;

	// Define minimum and maximum bandwidth as a fraction of graph diameter
	double max_bw = max_bw_factor * graph_diameter;
	double min_bw = min_bw_factor * graph_diameter;

	size_t domain_min_size = n_evectors + 1; // minimal number of elements in a general position(!) for the linear model of n_evectors variable to be fittable

	int n_vertices = adjacency_list.size(); // Number of vertices
	graph_spectral_lowess_result_t result;
    result.predictions.resize(n_vertices, INFINITY);
    result.errors.resize(n_vertices, INFINITY);
    result.scale.resize(n_vertices, INFINITY);

	// Step 1: Construct the adjacency matrix and degree matrix (unchanged)
	// Adjacency matrix as a sparse matrix
	Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
	std::vector<Eigen::Triplet<double>> triples;
	for (int i = 0; i < n_vertices; ++i) {
		for (const auto& edge : adjacency_list[i]) {
			int j = edge.vertex;
			// Use edge length as weight if specified, otherwise use binary weights
			double weight = use_edge_lengths ? edge.weight : 1.0;
			if (i < j) {  // Ensure each edge is added only once
				triples.push_back(Eigen::Triplet<double>(i, j, weight));
				triples.push_back(Eigen::Triplet<double>(j, i, weight));
			}
		}
	}
	A.setFromTriplets(triples.begin(), triples.end());

	// Compute Laplacian matrix as a sparse matrix
	Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
	for (int k = 0; k < A.outerSize(); ++k) {
		double sum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
			sum += it.value();
		}
		D.insert(k, k) = sum;
	}
	Eigen::SparseMatrix<double> L = D - A;

	// Add a small regularization term to improve numerical stability
	for (int i = 0; i < n_vertices; i++) {
		L.coeffRef(i, i) += 1e-8;
	}

	// Step 3: Compute eigenvalues and eigenvectors of the Laplacian
	// Eigenvalue decomposition
	Spectra::SparseSymMatProd<double> op(L);

	// Ensure ncv is within bounds: nev < ncv <= n
	int nev = 5 * n_evectors;
	int ncv = std::min(4 * nev, n_vertices); // Adjust ncv to be within bounds
	// Ensure nev < ncv
	if (nev >= ncv) {
		nev = ncv - 1;
	}

	if (nev < 1) {
		nev = 1;
	}

	// Construct eigen solver object to find eigenvalues closest to 0
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
	eigs.init();
	int maxit = 1000; // Increase maximum iterations (default is often too low)
	double tol = 1e-16; // Adjust tolerance if needed
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);


	if (eigs.info() != Spectra::CompInfo::Successful)
		REPORT_ERROR("Eigenvalue estimation with Spectra failed.");

	{
		// for testing/debugging only
		Eigen::VectorXd eigenvalues = eigs.eigenvalues();
		Rprintf("Eigenvalue 1: %.8f\n", eigenvalues(0));
		Rprintf("Eigenvalue 2: %.8f\n", eigenvalues(1));
		Rprintf("Eigenvalue 3: %.8f\n", eigenvalues(2));
		Rprintf("Eigenvalue 4: %.8f\n", eigenvalues(3));
		Rprintf("Eigenvalue 5: %.8f\n", eigenvalues(4));
		Rprintf("Eigenvalue 6: %.8f\n", eigenvalues(5));

		// Eigen::MatrixXd evectors = eigs.eigenvectors();
		// Eigen::VectorXd eig_vec2 = evectors.col(1);
		// Eigen::VectorXd eig_vec3 = evectors.col(2);
		// Eigen::VectorXd eig_vec4 = evectors.col(3);
		// Eigen::VectorXd eig_vec5 = evectors.col(4);
		// Eigen::VectorXd eig_vec6 = evectors.col(5);
	}

	// Step 4: For each vertex compute ...

	// Using for_each with parallel execution policy
    std::vector<size_t> vertices(n_vertices);
    std::iota(vertices.begin(), vertices.end(), 0);

    std::for_each(std::execution::seq, vertices.begin(), vertices.end(),
				  [&](size_t vertex) {

					  // Ensure minimum bandwidth requirements
					  double min_min_bw = find_minimum_radius_for_domain_min_size(
						  vertex,
						  min_bw,
						  max_bw,
						  domain_min_size,
						  precision
						  );

					  if (min_bw < min_min_bw) {
						  min_bw = min_min_bw;
					  }
					  if (min_bw >= max_bw) {
						  REPORT_ERROR("Minimum bandwidth (%f) >= maximum bandwidth (%f)\n",
									   min_bw, max_bw_factor * graph_diameter);
					  }

					  // Generate candidate bandwidths
					  std::vector<double> candidate_bws = get_candidate_bws(
						  min_bw,
						  max_bw,
						  n_bws,
						  log_grid,
						  precision
						  );

					  // Find vertices within the given radius
					  std::unordered_map<size_t, double> reachable_vertices = find_vertices_within_radius(vertex, max_bw);

					  // Evaluate each bandwidth
					  std::vector<double> local_predictions(candidate_bws.size());
					  std::vector<double> local_error(candidate_bws.size(), 0.0);

					  std::unordered_map<size_t, lm_t> model_map;

					  // For each bandwidth, fit models and compute errors
					  for (size_t bw_index = 0; bw_index < candidate_bws.size(); ++bw_index) {
						  double current_bw = candidate_bws[bw_index];

						  // Filter vertices by current bandwidth
						  std::map<size_t, double> local_vertex_map;
						  for (const auto& [v, dist] : reachable_vertices) {
							  if (dist <= current_bw) {
								  local_vertex_map[v] = dist;
							  }
						  }

						  // Create embedding using n_evectors eigenvectors
						  Eigen::MatrixXd embedding = create_spectral_embedding(local_vertex_map, n_evectors);

						  // Fit weighted linear model and compute LOOCV errors
						  lm_t lm_fit = fit_linear_model(embedding, y, local_vertex_map,
														 vertex,
														 dist_normalization_factor,
														 local_predictions[bw_index]);

						  local_error[bw_index] = lm_fit.mean_error;
						  model_map[bw_index]   = std::move(lm_fit);
					  }

					  // Find optimal bandwidth
					  size_t opt_bw_index = std::min_element(local_error.begin(), local_error.end()) - local_error.begin();
					  double opt_bw = candidate_bws[opt_bw_index];

					  // Store results
					  result.predictions[vertex] = std::move(model_map[opt_bw_index].predictions);
					  result.errors[vertex]      = std::move(model_map[opt_bw_index].errors);
					  result.scale[vertex]       = opt_bw;
		}
	);

	return result;
}
