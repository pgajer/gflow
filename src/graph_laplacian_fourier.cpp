#include <R.h>
#include <Rinternals.h>

#include <Eigen/Dense>
#include <vector>
#include <optional>

#include "graph_laplacian_fourier.hpp"

/**
 * @brief Creates a weighted chain graph from sorted x values
 *
 * Constructs a chain graph where consecutive points are connected with edges
 * weighted by their Euclidean distance.
 *
 * @param x Sorted vector of x coordinates
 * @return graph_data_t Structure containing the adjacency and Laplacian matrices
 *
 * @pre x must be sorted in ascending order
 * @pre x must contain at least 2 points
 *
 * @note The function assumes x is sorted. If x is not sorted, the results will be incorrect.
 */
graph_data_t create_chain_graph(const Eigen::VectorXd& x) {
	graph_data_t graph;
	int n = x.size();

	// Initialize adjacency matrix
	graph.adjacency_matrix = Eigen::MatrixXd::Zero(n, n);

	// Construct weighted edges based on distances between consecutive points
	for (int i = 0; i < n - 1; ++i) {
		double weight = std::abs(x(i + 1) - x(i));
		graph.adjacency_matrix(i, i + 1) = weight;
		graph.adjacency_matrix(i + 1, i) = weight;
	}

	// Compute Laplacian matrix: L = D - A
	Eigen::VectorXd degrees = graph.adjacency_matrix.rowwise().sum();
	// Convert diagonal matrix to regular matrix before subtraction
	graph.laplacian_matrix = degrees.asDiagonal().toDenseMatrix() - graph.adjacency_matrix;

	return graph;
}

/**
 * @brief Computes the eigendecomposition of the graph Laplacian matrix
 *
 * Updates the graph_data structure with computed eigenvectors and eigenvalues
 * of the Laplacian matrix.
 *
 * @param[in,out] graph Structure containing the Laplacian matrix to decompose
 *
 * @pre graph.laplacian_matrix must be symmetric
 * @post graph.eigenvectors and graph.eigenvalues will be populated
 */
void compute_laplacian_eigendecomposition(graph_data_t& graph) {
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(graph.laplacian_matrix);
	graph.eigenvalues = solver.eigenvalues();
	graph.eigenvectors = solver.eigenvectors();
}

/**
 * @brief Computes Fourier coefficients of a signal in the graph spectral domain
 *
 * Projects the signal onto the eigenvectors of the graph Laplacian to obtain
 * its Fourier coefficients.
 *
 * @param y Signal values vector
 * @param graph Graph structure containing the eigenvectors
 * @return Eigen::VectorXd Vector of Fourier coefficients
 *
 * @pre y.size() must match the number of graph vertices
 * @pre graph.eigenvectors must be computed (via compute_laplacian_eigendecomposition)
 */
Eigen::VectorXd compute_fourier_coefficients(
	const Eigen::VectorXd& y,
	const graph_data_t& graph
	) {
	return graph.eigenvectors.transpose() * y;
}

/**
 * @brief Main analysis function that performs spectral analysis of a signal on a graph
 *
 * Creates or reuses a chain graph structure and computes the Fourier coefficients
 * of the input signal with respect to the graph Laplacian eigenbasis.
 *
 * @param x Sorted vector of x coordinates
 * @param y Signal values vector
 * @param existing_graph Optional existing graph structure to reuse
 * @return spectral_result_t Complete analysis results
 *
 * @pre x and y must have the same size
 * @pre If x is provided, it must be sorted in ascending order
 * @pre If existing_graph is provided, its size must match y.size()
 *
 * @note When existing_graph is provided, x is ignored
 */
spectral_result_t analyze_signal(
	const Eigen::VectorXd& x,
	const Eigen::VectorXd& y,
	const std::optional<graph_data_t>& existing_graph) {
	spectral_result_t result;

	// Create new graph or use existing one
	if (!existing_graph) {
		result.graph_data = create_chain_graph(x);
		compute_laplacian_eigendecomposition(result.graph_data);
	} else {
		result.graph_data = *existing_graph;
	}

	// Compute Fourier coefficients
	result.fourier_coeffs = compute_fourier_coefficients(y, result.graph_data);

return result;
}
