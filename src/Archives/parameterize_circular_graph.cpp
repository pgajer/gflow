#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseSymMatProd.h>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

#include "set_wgraph.hpp"
#include "error_utils.h" // for REPORT_ERROR()

/**
 * Parameterize a circular graph structure using spectral methods.
 *
 * @param use_edge_lengths Whether to use edge lengths as weights in the graph Laplacian
 * @return A struct containing angles and eigenvectors for the circular parameterization
 */
circular_param_result_t set_wgraph_t::parameterize_circular_graph(
	bool use_edge_lengths) const {

	// Set Eigen to use 12 threads for parallel computation
	Eigen::setNbThreads(12);

	int n_vertices = adjacency_list.size(); // Number of vertices

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
	int nev = 3;
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
	//eigs.compute(Spectra::SortRule::SmallestAlge);
	int maxit = 1000; // Increase maximum iterations (default is often too low)
	double tol = 1e-16; // Adjust tolerance if needed
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);


	if (eigs.info() != Spectra::CompInfo::Successful)
		REPORT_ERROR("Eigenvalue estimation with Spectra failed.");

	Eigen::MatrixXd evectors = eigs.eigenvectors();

	// The first eigenvector corresponds to the smallest eigenvalue
	// For a connected graph, this is approximately constant
	// We want the 2nd and 3rd eigenvectors
	Eigen::VectorXd eig_vec2 = evectors.col(1);
	Eigen::VectorXd eig_vec3 = evectors.col(2);

	// Step 5: Calculate circular coordinates (angles) from the eigenvectors
	std::vector<double> angles(n_vertices);
	for (int i = 0; i < n_vertices; i++) {
		angles[i] = std::atan2(eig_vec3(i), eig_vec2(i));
		if (angles[i] < 0) {
			angles[i] += 2 * M_PI;
		}
	}

	// Return the struct with all the results
	circular_param_result_t result;
	result.angles = angles;
	result.eig_vec2.resize(n_vertices);
	result.eig_vec3.resize(n_vertices);

	// Copy Eigen vectors to std::vectors
	for (int i = 0; i < n_vertices; ++i) {
		result.eig_vec2[i] = eig_vec2(i);
		result.eig_vec3[i] = eig_vec3(i);
	}

	return result;
}

// this version is very slow
#if 0
circular_param_result_t set_wgraph_t::parameterize_circular_graph(
	bool use_edge_lengths) const {

	// Set Eigen to use multiple threads for parallel computation
	Eigen::setNbThreads(12);

	int n_vertices = adjacency_list.size();  // Number of vertices

	// Step 1: Construct the adjacency matrix and degree matrix
	// Adjacency matrix as a sparse matrix
	Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
	std::vector<Eigen::Triplet<double>> triples;
	for (int i = 0; i < n_vertices; ++i) {
		for (const auto& edge : adjacency_list[i]) {
			int j = edge.vertex;
			// Use edge length as weight if specified, otherwise use binary weights
			double weight = use_edge_lengths ? edge.weight : 1.0;
			if (i < j) {   // Ensure each edge is added only once
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

	// Step 3: Compute eigenvalues and eigenvectors of the Laplacian using Eigen
	// Convert to dense for eigendecomposition - for large graphs, consider iterative methods
	Eigen::MatrixXd L_dense = Eigen::MatrixXd(L);

	// Perform eigendecomposition
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(L_dense);

	if (eigen_solver.info() != Eigen::Success) {
		REPORT_ERROR("Eigenvalue computation with Eigen failed.");
	}

	// The eigenvectors are sorted by eigenvalue (ascending)
	// The first eigenvector (constant) is not useful, so we take the 2nd and 3rd
	Eigen::VectorXd eig_vec2 = eigen_solver.eigenvectors().col(1);
	Eigen::VectorXd eig_vec3 = eigen_solver.eigenvectors().col(2);

	// Step 5: Calculate circular coordinates (angles) from the eigenvectors
	std::vector<double> angles(n_vertices);
	for (int i = 0; i < n_vertices; i++) {
		angles[i] = std::atan2(eig_vec3(i), eig_vec2(i));
		if (angles[i] < 0) {
			angles[i] += 2 * M_PI;
		}
	}

	// Return the struct with all the results
	circular_param_result_t result;
	result.angles = angles;
	result.eig_vec2.resize(n_vertices);
	result.eig_vec3.resize(n_vertices);

	// Copy Eigen vectors to std::vectors
	for (int i = 0; i < n_vertices; ++i) {
		result.eig_vec2[i] = eig_vec2(i);
		result.eig_vec3[i] = eig_vec3(i);
	}

	return result;
}
#endif


#if 0
circular_param_result_t set_wgraph_t::parameterize_circular_graph(
	bool use_edge_lengths) const {

	size_t n_vertices = adjacency_list.size(); // Number of vertices

	// Step 1: Construct the adjacency matrix and degree matrix
	Eigen::SparseMatrix<double> adj_matrix(n_vertices, n_vertices);// Adjacency matrix
	Eigen::VectorXd degree_vec = Eigen::VectorXd::Zero(n_vertices);// Diagonal of degree matrix

	// Fill in the adjacency matrix
	std::vector<Eigen::Triplet<double>> triplets;
	for (size_t i = 0; i < n_vertices; i++) {
		for (const auto& edge : adjacency_list[i]) {
			int j = edge.vertex;
			// Use edge length as weight if specified, otherwise use binary weights
			double weight = use_edge_lengths ? edge.weight : 1.0;

			triplets.push_back(Eigen::Triplet<double>(i, j, weight));
			degree_vec(i) += weight; // Add to degree
		}
	}
	adj_matrix.setFromTriplets(triplets.begin(), triplets.end());

	// Step 2: Construct the normalized Laplacian L = I - D^(-1/2) A D^(-1/2)
	Eigen::SparseMatrix<double> laplacian(n_vertices, n_vertices);
	Eigen::SparseMatrix<double> identity(n_vertices, n_vertices);
	identity.setIdentity();

	// Calculate D^(-1/2)
	Eigen::VectorXd degree_inv_sqrt(n_vertices);
	for (size_t i = 0; i < n_vertices; i++) {
		degree_inv_sqrt(i) = 1.0 / std::sqrt(std::max(degree_vec(i), 1e-10)); // Avoid division by zero
	}

	// Apply D^(-1/2) on both sides of A
	Eigen::SparseMatrix<double> degree_inv_sqrt_mat(n_vertices, n_vertices);
	for (size_t i = 0; i < n_vertices; i++) {
		degree_inv_sqrt_mat.insert(i, i) = degree_inv_sqrt(i);
	}

	Eigen::SparseMatrix<double> normalized_adj = degree_inv_sqrt_mat * adj_matrix * degree_inv_sqrt_mat;
	laplacian = identity - normalized_adj;

	// Step 3: Compute eigenvalues and eigenvectors of the Laplacian
	Eigen::SparseSymMatProd<double> op(laplacian);
	Eigen::SpectraSymShiftSolve<double, Eigen::Ascending> eigen_solver;
	eigen_solver.compute(op, 3, 10); // Get 3 smallest eigenvalues, with 10 Lanczos vectors

#if 0
	// Convert sparse to dense for eigenvalue computation
	Eigen::MatrixXd laplacian_dense = Eigen::MatrixXd(laplacian);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(laplacian_dense);

	if (eigen_solver.info() != Eigen::Success) {
		// Handle eigenvalue computation error
		REPORT_ERROR("Eigenvalue computation failed");
	}
#endif


	// Step 4: Extract the second and third eigenvectors (sorted by eigenvalue)
	// The first eigenvector (corresponding to smallest eigenvalue) is typically constant
	// and not useful for embedding
	Eigen::VectorXd eig_vec2 = eigen_solver.eigenvectors().col(1); // Second eigenvector
	Eigen::VectorXd eig_vec3 = eigen_solver.eigenvectors().col(2); // Third eigenvector

	// Step 5: Calculate circular coordinates (angles) from the eigenvectors
	std::vector<double> angles(n_vertices);
	for (size_t i = 0; i < n_vertices; i++) {
		// Convert to polar coordinates - atan2 gives angle in range [-π, π]
		angles[i] = std::atan2(eig_vec3(i), eig_vec2(i));

		// Normalize to [0, 2π) range if desired
		if (angles[i] < 0) {
			angles[i] += 2 * M_PI;
		}
	}

	// Return the struct with all the results
	circular_param_result_t result;
	result.angles = angles;
	result.eig_vec2.resize(eig_vec2.size());
	result.eig_vec3.resize(eig_vec3.size());

	// Copy Eigen vectors to std::vectors
	for (size_t i = 0; i < eig_vec2.size(); ++i) {
		result.eig_vec2[i] = eig_vec2(i);
		result.eig_vec3[i] = eig_vec3(i);
	}

	return result;
}
#endif

/**
 * Alternative version that allows specifying a reference vertex to set as angle 0
 *
 * @param reference_vertex The vertex to use as reference (angle 0)
 * @param use_edge_lengths Whether to use edge lengths as weights in the graph Laplacian
 * @return A struct containing angles and eigenvectors for the circular parameterization
 */
circular_param_result_t set_wgraph_t::parameterize_circular_graph_with_reference(
	size_t reference_vertex,
	bool use_edge_lengths) const {

	// Get the results using spectral method
	circular_param_result_t result = parameterize_circular_graph(use_edge_lengths);

	if (reference_vertex < adjacency_list.size()) {
		// Get the angle of the reference vertex
		double ref_angle = result.angles[reference_vertex];

		// Shift all angles so that reference_vertex has angle 0
		for (size_t i = 0; i < result.angles.size(); i++) {
			result.angles[i] = fmod(result.angles[i] - ref_angle + 2 * M_PI, 2 * M_PI);
		}
	}

	return result;
}

// Example usage:
/*
  void example_usage() {
  set_wgraph_t my_graph;
  // ... populate my_graph ...

  // Get circular parameterization
  circular_param_result_t result = my_graph.parameterize_circular_graph(true);

  // Access the results
  for (size_t i = 0; i < result.angles.size(); i++) {
  std::cout << "Vertex " << i << ": angle = " << result.angles[i] << std::endl;
  }

  // If you want to access the eigenvectors for further analysis
  std::cout << "2nd eigenvector: " << result.eig_vec2.transpose() << std::endl;
  std::cout << "3rd eigenvector: " << result.eig_vec3.transpose() << std::endl;

  // Or use the reference vertex version
  circular_param_result_t ref_result = my_graph.parameterize_circular_graph_with_reference(0, true);
  // ... use ref_result ...
}
*/
