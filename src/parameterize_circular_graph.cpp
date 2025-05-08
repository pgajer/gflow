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

// normalized laplacian
#if 0
circular_param_result_t set_wgraph_t::parameterize_circular_graph(
	bool use_edge_lengths) const {
	// Set Eigen to use 12 threads for parallel computation
	Eigen::setNbThreads(12);
	int n_vertices = adjacency_list.size();  // Number of vertices

	// Step 1: Construct the adjacency matrix and degree matrix (unchanged)
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

	// Compute degree matrix D and D^(-1/2) directly
	Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
	Eigen::SparseMatrix<double> D_inv_sqrt(n_vertices, n_vertices);

	for (int k = 0; k < A.outerSize(); ++k) {
		double sum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
			sum += it.value();
		}
		// Add small epsilon to avoid division by zero
		double d_value = sum + 1e-10;
		D.insert(k, k) = d_value;
		D_inv_sqrt.insert(k, k) = 1.0 / std::sqrt(d_value);
	}

	// Compute normalized Laplacian L = I - D^(-1/2) A D^(-1/2)
	Eigen::SparseMatrix<double> I(n_vertices, n_vertices);
	I.setIdentity();

	// Efficient computation of D^(-1/2) A D^(-1/2)
	Eigen::SparseMatrix<double> temp = D_inv_sqrt * A * D_inv_sqrt;
	Eigen::SparseMatrix<double> L = I - temp;

	// Add a small regularization term to improve numerical stability
	for (int i = 0; i < n_vertices; i++) {
		L.coeffRef(i, i) += 1e-8;
	}

	// Step 3: Compute eigenvalues and eigenvectors of the Laplacian
	// Eigenvalue decomposition
	Spectra::SparseSymMatProd<double> op(L);

	// Ensure ncv is within bounds: nev < ncv <= n
	int nev = 6;
	int ncv = std::min(4 * nev, n_vertices);  // Adjust ncv to be within bounds

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

	int maxit = 1000;  // Increase maximum iterations
	double tol = 1e-16;  // Adjust tolerance
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

	if (eigs.info() != Spectra::CompInfo::Successful)
		REPORT_ERROR("Eigenvalue estimation with Spectra failed.");


	// debugging
	Eigen::VectorXd eigenvalues = eigs.eigenvalues();
	Rprintf("Eigenvalue 1: %.8f\n", eigenvalues(0));
	Rprintf("Eigenvalue 2: %.8f\n", eigenvalues(1));
	Rprintf("Eigenvalue 3: %.8f\n", eigenvalues(2));
	Rprintf("Eigenvalue 4: %.8f\n", eigenvalues(3));
	Rprintf("Eigenvalue 5: %.8f\n", eigenvalues(4));
	Rprintf("Eigenvalue 6: %.8f\n", eigenvalues(5));


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
#endif

// this function uses standard laplacian and produces figure eight instead of circular graph when using the second and third eigenvectors
#if 1
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
	int nev = 10;
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

	// debugging
	Eigen::VectorXd eigenvalues = eigs.eigenvalues();
	Rprintf("Eigenvalue 1: %.8f\n", eigenvalues(0));
	Rprintf("Eigenvalue 2: %.8f\n", eigenvalues(1));
	Rprintf("Eigenvalue 3: %.8f\n", eigenvalues(2));
	Rprintf("Eigenvalue 4: %.8f\n", eigenvalues(3));
	Rprintf("Eigenvalue 5: %.8f\n", eigenvalues(4));
	Rprintf("Eigenvalue 6: %.8f\n", eigenvalues(5));

	Eigen::MatrixXd evectors = eigs.eigenvectors();

	// The first eigenvector corresponds to the smallest eigenvalue
	// For a connected graph, this is approximately constant
	// We want the 2nd and 3rd eigenvectors
	Eigen::VectorXd eig_vec2 = evectors.col(1);
	Eigen::VectorXd eig_vec3 = evectors.col(2);
	Eigen::VectorXd eig_vec4 = evectors.col(3);
	Eigen::VectorXd eig_vec5 = evectors.col(4);
	Eigen::VectorXd eig_vec6 = evectors.col(5);

	// Step 5: Calculate circular coordinates (angles) from the eigenvectors
	std::vector<double> angles(n_vertices);
	// Normalize eigenvectors first to ensure uniform scaling
	double norm2 = eig_vec2.norm();
	double norm3 = eig_vec3.norm();
	eig_vec2 /= norm2;
	eig_vec3 /= norm3;

	// Calculate angles differently
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
	result.eig_vec4.resize(n_vertices);
	result.eig_vec5.resize(n_vertices);
	result.eig_vec6.resize(n_vertices);


	// Copy Eigen vectors to std::vectors
	for (int i = 0; i < n_vertices; ++i) {
		result.eig_vec2[i] = eig_vec2(i);
		result.eig_vec3[i] = eig_vec3(i);
		result.eig_vec4[i] = eig_vec4(i);
		result.eig_vec5[i] = eig_vec5(i);
		result.eig_vec6[i] = eig_vec6(i);
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
