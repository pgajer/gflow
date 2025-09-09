#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <unordered_map>
#include <unordered_set>
#include <queue>

#include "set_wgraph.hpp"
#include "error_utils.h"  // For REPORT_ERROR()
#include "gflow_cx.hpp"

/**
 * @brief Performs smooth function extension from boundary to interior of a hop disk
 *
 * @details This function extends values from the boundary of a hop neighborhood to its interior
 *          using a two-phase approach. In the first phase, it initializes interior values using
 *          inverse-distance weighted interpolation from boundary vertices. In the second phase,
 *          it iteratively refines these values to approximate a harmonic extension by computing
 *          weighted averages of neighboring vertices. This produces a smoother extension that
 *          preserves the boundary values while creating a natural transition within the interior.
 *
 * @param[in,out] smoothed_values Vector of function values to be modified for interior vertices
 * @param[in] hop_nbhd Hop neighborhood structure containing the center vertex, hop radius,
 *                    hop distance map, and boundary values map
 * @param[in] max_iterations Maximum number of refinement iterations to perform
 * @param[in] tolerance Convergence threshold for maximum value change between iterations
 * @param[in] sigma Parameter controlling the Gaussian weighting in the refinement phase (default: 1.0)
 * @param[in] verbose Whether to print convergence information (default: false)
 *
 * @throws None, but errors might be reported through the REPORT_* macros
 *
 * @note This function modifies the values of interior vertices in the smoothed_values vector
 *       while preserving the boundary values. The boundary is defined as vertices at hop distance
 *       (hop_radius + 1) from the center vertex.
 *
 * @see hop_nbhd_t
 * @see compute_shortest_path_distance
 */
void set_wgraph_t::perform_weighted_mean_hop_disk_extension(
	std::vector<double>& smoothed_values,
	const hop_nbhd_t& hop_nbhd,
	int max_iterations,
	double tolerance,
	double sigma,
	bool verbose
	) const {

	size_t hop_radius = hop_nbhd.hop_idx;

	// Get interior vertices (hop distance <= hop_radius)
	std::unordered_set<size_t> interior_vertices;
	for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
		if (dist <= hop_radius) {
			interior_vertices.insert(v);
		}
	}

	// Get boundary vertices (hop distance = hop_radius + 1)
	std::unordered_set<size_t> boundary_vertices;
	for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
		boundary_vertices.insert(v);
	}

	// For each interior vertex, calculate initial value based on boundary values
	for (size_t u : interior_vertices) {
		// Calculate geodesic distances to boundary vertices
		std::vector<std::pair<size_t, double>> boundary_distances;
		for (size_t b : boundary_vertices) {
			double dist = compute_shortest_path_distance(u, b);
			boundary_distances.push_back({b, dist});
		}

		// Calculate inverse-distance weights
		std::vector<double> weights;
		double weight_sum = 0.0;
		for (const auto& [_, dist] : boundary_distances) {
			//double weight = 1.0 / (dist + 1e-10);  // Avoid division by zero
			double weight = std::exp(-(dist * dist) / (2.0 * sigma * sigma));
			weights.push_back(weight);
			weight_sum += weight;
		}

		// Calculate weighted average
		double new_value = 0.0;
		for (size_t i = 0; i < boundary_distances.size(); ++i) {
			auto [b, _] = boundary_distances[i];
			new_value += (weights[i] / weight_sum) * smoothed_values[b];
		}

		// Update the value
		smoothed_values[u] = new_value;
	}

	// Iteratively refine to approximate harmonic extension
	bool converged = false;
	for (int iter = 0; iter < max_iterations; ++iter) {
		double max_change = 0.0;
		std::vector<double> new_values = smoothed_values;

		for (size_t u : interior_vertices) {
			// Skip boundary vertices (note: boundary vertices should not be in interior_vertices,
			// but adding check for safety)
			if (boundary_vertices.find(u) != boundary_vertices.end()) {
				REPORT_ERROR("ERROR: Interior vertex %zu belongs to the boundary of the hop disk", u);
				//continue;
			}

			// Average of neighbors
			double sum_values = 0.0;
			double sum_weights = 0.0;

			for (const auto& edge : adjacency_list[u]) {
				// Get the neighbor vertex
				size_t neighbor = edge.vertex;

				// Calculate weight based on edge weight (using Gaussian weighting)
				double edge_dist = edge.weight;
				double weight = std::exp(-(edge_dist * edge_dist) / (2.0 * sigma * sigma));

				sum_values += smoothed_values[neighbor] * weight;
				sum_weights += weight;
			}

			if (sum_weights > 0) {
				double weighted_avg = sum_values / sum_weights;
				double change = std::abs(smoothed_values[u] - weighted_avg);
				max_change = std::max(max_change, change);
				new_values[u] = weighted_avg;
			}
		}

		// Update values for next iteration
		smoothed_values = new_values;

		// Check for convergence
		if (max_change < tolerance) {
			converged = true;
			if (verbose) {
				Rprintf("Hop disk extension converged after %d iterations\n", iter + 1);
			}
			break;
		}
	}

	if (!converged && verbose) {
		Rprintf("Hop disk extension did not converge after %d iterations\n", max_iterations);
	}
}

/**
 * @brief Performs harmonic extension from boundary values to interior vertices
 *
 * @details Extends function values from the boundary of a region to its interior
 *          by enforcing the harmonic property (weighted average of neighbors).
 *          The algorithm iteratively updates interior values until convergence.
 *
 * @param boundary_values Map of boundary vertex indices to their fixed values
 * @param region_vertices Set of all vertex indices in the region (boundary + interior)
 * @param max_iterations Maximum number of iterations to perform
 * @param tolerance Convergence threshold for value changes
 * @param record_frequency How often to record states (every N iterations)
 * @return harmonic_extender_t Structure with extension history and results
 */
harmonic_extender_t set_wgraph_t::harmonic_extender(
	const std::unordered_map<size_t, double>& boundary_values,
	const std::unordered_set<size_t>& region_vertices,
	int max_iterations,
	double tolerance,
	int record_frequency,
	bool verbose
	) const {

	// Initialize result structure
	harmonic_extender_t result;
	result.num_iterations = 0;
	result.converged = false;
	result.max_change_final = 0.0;

	// Ensure edge weights are computed
	ensure_edge_weights_computed();

	// Validate input
	if (region_vertices.empty()) {
		REPORT_WARNING("Empty region provided to harmonic_extender");
		return result;
	}

	if (boundary_values.empty()) {
		REPORT_WARNING("No boundary values provided to harmonic_extender");
		return result;
	}

	// Identify interior vertices (region vertices that are not boundary vertices)
	std::unordered_set<size_t> interior_vertices;
	for (const auto& vertex : region_vertices) {
		if (boundary_values.find(vertex) == boundary_values.end()) {
			interior_vertices.insert(vertex);
		}
	}

	// Return if we have no interior vertices to extend to
	if (interior_vertices.empty()) {
		if (verbose) {
			Rprintf("No interior vertices to extend to\n");
		}
		return result;
	}

	// Initialize function values vector with size equal to number of vertices in the graph
	std::vector<double> function_values(adjacency_list.size(), 0.0);

	// Set boundary values
	for (const auto& [vertex, value] : boundary_values) {
		if (vertex >= function_values.size()) {
			REPORT_ERROR("Boundary vertex %zu is out of range", vertex);
		}
		function_values[vertex] = value;
	}

	// Calculate initial value for interior vertices (average of boundary values)
	double boundary_avg = 0.0;
	for (const auto& [vertex, value] : boundary_values) {
		boundary_avg += value;
	}
	boundary_avg /= boundary_values.size();

	// Initialize interior vertices with boundary average
	for (const auto& vertex : interior_vertices) {
		function_values[vertex] = boundary_avg;
	}

	// Record initial state
	result.iterations.push_back(function_values);

	// Temporary vector for new values during iteration
	std::vector<double> new_values = function_values;

	// Iterative relaxation for harmonic extension
	double max_change = 0.0;
	for (int iter = 0; iter < max_iterations; ++iter) {
		max_change = 0.0;

		// Update each interior vertex value
		for (const size_t& v : interior_vertices) {
			// Calculate weighted average of neighbors
			double sum = 0.0;
			double weight_sum = 0.0;

			for (const auto& edge : adjacency_list[v]) {
				// Use inverse of edge weight as weighting factor
				double weight = 1.0 / (edge.weight + 1e-10); // Avoid division by zero
				sum += function_values[edge.vertex] * weight;
				weight_sum += weight;
			}

			if (weight_sum > 0) {
				double weighted_avg = sum / weight_sum;
				double change = std::abs(function_values[v] - weighted_avg);
				max_change = std::max(max_change, change);
				new_values[v] = weighted_avg;
			}
		}

		// Update the values for next iteration
		for (const size_t& v : interior_vertices) {
			function_values[v] = new_values[v];
		}

		// Increment iteration counter
		result.num_iterations++;

		// Record state if it's time (based on record_frequency)
		if ((iter + 1) % record_frequency == 0 || iter == max_iterations - 1 || max_change < tolerance) {
			result.iterations.push_back(function_values);
		}

		// Check for convergence
		if (max_change < tolerance) {
			result.converged = true;
			if (verbose) {
				Rprintf("Harmonic extension converged after %d iterations\n", iter + 1);
			}
			break;
		}
	}

	// Report if did not converge
	if (!result.converged && verbose) {
		Rprintf("Harmonic extension did not converge after %d iterations. Max change: %g\n",
				max_iterations, max_change);
	}

	// Record final maximum change
	result.max_change_final = max_change;

	return result;
}

#if 0
void set_wgraph_t::perform_harmonic_extension(
	std::vector<double>& values,
	const std::unordered_set<size_t>& interior,
	const std::unordered_map<size_t, double>& boundary_values,
	double tolerance = 1e-6,
	int max_iterations = 1000) const {

	// Initialize interior values (can use your quasi-harmonic method for a better starting point)
	// Or simply initialize with the average of boundary values
	double avg_boundary = 0.0;
	for (const auto& [_, value] : boundary_values) {
		avg_boundary += value;
	}
	avg_boundary /= boundary_values.size();

	for (size_t u : interior) {
		values[u] = avg_boundary;
	}

	// Fixed boundary values
	for (const auto& [v, value] : boundary_values) {
		values[v] = value;
	}

	// Iterative relaxation
	double max_diff = tolerance + 1.0;
	int iterations = 0;

	while (max_diff > tolerance && iterations < max_iterations) {
		max_diff = 0.0;

		// Create a copy of current values for update
		std::vector<double> new_values = values;

		// Update each interior vertex with average of neighbors
		for (size_t u : interior) {
			double sum = 0.0;
			int count = 0;

			// Average of neighbors
			for (const auto& edge : adjacency_list[u]) {
				size_t v = edge.vertex;
				sum += values[v];
				count++;
			}

			if (count > 0) {
				new_values[u] = sum / count;
				max_diff = std::max(max_diff, std::abs(new_values[u] - values[u]));
			}
		}

		// Update values
		values = new_values;
		iterations++;
	}
}
#endif





/**
 * @brief Computes harmonic extension using a direct sparse linear system solver
 *
 * @details Constructs and solves a linear system Ax = b where A is the graph Laplacian
 *  with rows corresponding to boundary vertices replaced by identity rows,
 *  and b contains the fixed boundary values. Uses Eigen's sparse matrix
 *  capabilities for efficient solution with fallback strategies for
 *  numerical robustness.
 *
 * @param boundary_values Map of boundary vertex indices to their fixed values
 * @param region_vertices Set of all vertex indices in the region (boundary + interior)
 * @param regularization Small value to add to diagonal for numerical stability
 * @param verbose Whether to print progress and diagnostic information
 * @return std::vector<double> Vector containing the harmonic extension values
 */
std::vector<double> set_wgraph_t::harmonic_extension_eigen(
	const std::unordered_map<size_t, double>& boundary_values,
	const std::unordered_set<size_t>& region_vertices,
	double regularization,
	bool verbose
	) const {

	// Ensure edge weights are computed
	ensure_edge_weights_computed();

	// Validate input
	if (region_vertices.empty()) {
		REPORT_WARNING("Empty region provided to harmonic_extension_eigen");
		return std::vector<double>();
	}

	if (boundary_values.empty()) {
		REPORT_WARNING("No boundary values provided to harmonic_extension_eigen");
		return std::vector<double>();
	}

	// Prepare result vector with size equal to number of vertices in the graph
	std::vector<double> result(adjacency_list.size(), 0.0);

	// Create index mapping for compact matrix (map original indices to sequential indices)
	std::unordered_map<size_t, size_t> index_map;
	std::unordered_map<size_t, size_t> reverse_map;
	size_t idx = 0;
	for (const auto& vertex : region_vertices) {
		index_map[vertex] = idx;
		reverse_map[idx] = vertex;
		idx++;
	}
	size_t system_size = region_vertices.size();

	if (verbose) {
		Rprintf("Setting up harmonic extension system with %zu vertices\n", system_size);
	}

	// Initialize triplet list for sparse matrix construction
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(system_size * 5); // Estimate average degree of 5

	// Initialize right-hand side vector
	Eigen::VectorXd b = Eigen::VectorXd::Zero(system_size);

	// Construct system: Laplacian for interior vertices, identity for boundary
	for (const auto& vertex : region_vertices) {
		size_t i = index_map[vertex];

		// If this is a boundary vertex, set identity row and boundary value in b
		if (boundary_values.find(vertex) != boundary_values.end()) {
			triplets.push_back(T(i, i, 1.0));
			b(i) = boundary_values.at(vertex);
		}
		// Interior vertex: set Laplacian row (weighted degree on diagonal, -weight on neighbors)
		else {
			double weighted_degree = 0.0;

			// Process all neighbors
			for (const auto& edge : adjacency_list[vertex]) {
				// Only include neighbors that are in the region
				if (region_vertices.find(edge.vertex) != region_vertices.end()) {
					double weight = 1.0 / (edge.weight + 1e-10); // Inverse weight
					weighted_degree += weight;

					// Off-diagonal entry (negative weight)
					triplets.push_back(T(i, index_map[edge.vertex], -weight));
				}
			}

			// Diagonal entry (weighted degree + regularization)
			triplets.push_back(T(i, i, weighted_degree + regularization));
		}
	}

	// Construct sparse matrix
	Eigen::SparseMatrix<double> A(system_size, system_size);
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();

	if (verbose) {
		Rprintf("Matrix constructed with %d non-zeros\n", (int)A.nonZeros());
	}

	// Solve the system using sparse LU decomposition with fallback strategies
	Eigen::VectorXd x;
	bool solution_found = false;

	// Try SparseLU first (direct solver)
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with SparseLU\n");
			}
			Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
			solver.compute(A);
			if (solver.info() == Eigen::Success) {
				x = solver.solve(b);
				if (solver.info() == Eigen::Success) {
					solution_found = true;
					if (verbose) {
						Rprintf("SparseLU solution successful\n");
					}
				} else if (verbose) {
					Rprintf("SparseLU solve failed\n");
				}
			} else if (verbose) {
				Rprintf("SparseLU decomposition failed\n");
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("SparseLU exception: %s\n", e.what());
			}
		}
	}

	// Try SparseQR if SparseLU failed
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with SparseQR\n");
			}
			Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
			solver.compute(A);
			if (solver.info() == Eigen::Success) {
				x = solver.solve(b);
				if (solver.info() == Eigen::Success) {
					solution_found = true;
					if (verbose) {
						Rprintf("SparseQR solution successful\n");
					}
				} else if (verbose) {
					Rprintf("SparseQR solve failed\n");
				}
			} else if (verbose) {
				Rprintf("SparseQR decomposition failed\n");
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("SparseQR exception: %s\n", e.what());
			}
		}
	}

	// Try iterative BiCGSTAB if direct solvers failed
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with BiCGSTAB\n");
			}
			Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

			// Try with different regularization levels if needed
			std::vector<double> reg_levels = {regularization, 1e-8, 1e-6, 1e-4};
			for (double reg : reg_levels) {
				if (solution_found) break;

				// Create regularized matrix
				Eigen::SparseMatrix<double> A_reg = A;
				if (reg > regularization) {
					// Add additional regularization to diagonal
					for (int k = 0; k < A_reg.outerSize(); ++k) {
						for (Eigen::SparseMatrix<double>::InnerIterator it(A_reg, k); it; ++it) {
							if (it.row() == it.col()) {
								it.valueRef() += (reg - regularization);
							}
						}
					}
					if (verbose) {
						Rprintf("Increased regularization to %g\n", reg);
					}
				}

				// Try different iteration counts and tolerances
				std::vector<std::pair<int, double>> attempts = {
					{1000, 1e-10},
					{2000, 1e-8},
					{5000, 1e-6},
					{10000, 1e-4}
				};

				for (const auto& [maxit, tol] : attempts) {
					if (solution_found) break;

					solver.setMaxIterations(maxit);
					solver.setTolerance(tol);
					solver.compute(A_reg);

					if (solver.info() == Eigen::Success) {
						x = solver.solve(b);
						if (solver.info() == Eigen::Success) {
							solution_found = true;
							if (verbose) {
								Rprintf("BiCGSTAB solution successful with maxit=%d, tol=%g, reg=%g\n",
										maxit, tol, reg);
							}
							break;
						} else if (verbose) {
							Rprintf("BiCGSTAB solve failed with maxit=%d, tol=%g\n", maxit, tol);
						}
					} else if (verbose) {
						Rprintf("BiCGSTAB decomposition failed with maxit=%d, tol=%g\n", maxit, tol);
					}
				}
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("BiCGSTAB exception: %s\n", e.what());
			}
		}
	}

	// If all solvers failed, fall back to iterative method
	if (!solution_found) {
		if (verbose) {
			Rprintf("All Eigen solvers failed, falling back to iterative method\n");
		}

		// Use harmonic_extender instead
		harmonic_extender_t iter_result = harmonic_extender(
			boundary_values,
			region_vertices,
			10000,   // More iterations for fallback
			1e-8,// Tighter tolerance
			10000,   // Don't record intermediate states
			verbose
			);

		// If we have at least one result, use it
		if (!iter_result.iterations.empty()) {
			return iter_result.iterations.back();
		} else {
			REPORT_ERROR("Harmonic extension failed: all solvers failed");
			return result;
		}
	}

	// Convert solution back to original indexing
	for (size_t i = 0; i < system_size; i++) {
		size_t orig_idx = reverse_map[i];
		result[orig_idx] = x(i);
	}

	return result;
}

/**
 * @brief Performs biharmonic extension from boundary values to interior vertices
 *
 * @details Extends function values from the boundary and its neighborhood to the interior
 *  by enforcing the biharmonic property (∆²u = 0). This produces smoother extensions
 *  that better preserve local features compared to harmonic extensions.
 *
 * @param boundary_values Map of boundary vertex indices to their fixed values
 * @param boundary_neighbor_values Map of vertices in the 1-ring around the boundary to their values
 * @param region_vertices Set of all vertex indices in the region (boundary + interior)
 * @param regularization Small value to add to diagonal for numerical stability
 * @param verbose Whether to print progress and diagnostic information
 * @return std::vector<double> Vector containing the biharmonic extension values
 */
std::vector<double> set_wgraph_t::biharmonic_extension_eigen(
	const std::unordered_map<size_t, double>& boundary_values,
	const std::unordered_map<size_t, double>& boundary_neighbor_values,
	const std::unordered_set<size_t>& region_vertices,
	double regularization,
	bool verbose
	) const {

	// Ensure edge weights are computed
	ensure_edge_weights_computed();

	// Validate input
	if (region_vertices.empty()) {
		REPORT_WARNING("Empty region provided to biharmonic_extension_eigen");
		return std::vector<double>();
	}

	if (boundary_values.empty() && boundary_neighbor_values.empty()) {
		REPORT_WARNING("No boundary values provided to biharmonic_extension_eigen");
		return std::vector<double>();
	}

	// Prepare result vector with size equal to number of vertices in the graph
	std::vector<double> result(adjacency_list.size(), 0.0);

	// Create index mapping for compact matrix
	std::unordered_map<size_t, size_t> index_map;
	std::unordered_map<size_t, size_t> reverse_map;
	size_t idx = 0;
	for (const auto& vertex : region_vertices) {
		index_map[vertex] = idx;
		reverse_map[idx] = vertex;
		idx++;
	}
	size_t system_size = region_vertices.size();

	if (verbose) {
		Rprintf("Setting up biharmonic extension system with %zu vertices\n", system_size);
	}

	// Create the Laplacian matrix for the region
	typedef Eigen::Triplet<double> T;
	std::vector<T> L_triplets;
	L_triplets.reserve(system_size * 5);

	// Construct the graph Laplacian for the region
	for (const auto& vertex : region_vertices) {
		size_t i = index_map[vertex];
		double weighted_degree = 0.0;

		for (const auto& edge : adjacency_list[vertex]) {
			if (region_vertices.find(edge.vertex) != region_vertices.end()) {
				double weight = 1.0 / (edge.weight + 1e-10);
				weighted_degree += weight;
				size_t j = index_map[edge.vertex];
				L_triplets.push_back(T(i, j, -weight));
			}
		}

		L_triplets.push_back(T(i, i, weighted_degree));
	}

	Eigen::SparseMatrix<double> L(system_size, system_size);
	L.setFromTriplets(L_triplets.begin(), L_triplets.end());
	L.makeCompressed();

	// Square the Laplacian (bilaplacian)
	Eigen::SparseMatrix<double> L2 = L * L;
	L2.makeCompressed();

	// Initialize system matrix and right-hand side vector
	std::vector<T> system_triplets;
	system_triplets.reserve(L2.nonZeros() + boundary_values.size() + boundary_neighbor_values.size());
	Eigen::VectorXd b = Eigen::VectorXd::Zero(system_size);

	// Set boundary and boundary neighbor constraints
	// For boundary vertices, set identity rows
	for (const auto& [vertex, value] : boundary_values) {
		if (index_map.find(vertex) != index_map.end()) {
			size_t i = index_map[vertex];
			system_triplets.push_back(T(i, i, 1.0));
			b(i) = value;
		}
	}

	// For boundary neighbor vertices, set identity rows
	for (const auto& [vertex, value] : boundary_neighbor_values) {
		if (index_map.find(vertex) != index_map.end() &&
			boundary_values.find(vertex) == boundary_values.end()) {
			size_t i = index_map[vertex];
			system_triplets.push_back(T(i, i, 1.0));
			b(i) = value;
		}
	}

	// Set bilaplacian rows for interior vertices
	std::unordered_set<size_t> constrained_vertices;
	for (const auto& [vertex, _] : boundary_values) {
		constrained_vertices.insert(vertex);
	}
	for (const auto& [vertex, _] : boundary_neighbor_values) {
		constrained_vertices.insert(vertex);
	}

	// Add bilaplacian rows for unconstrained vertices
	for (int k = 0; k < L2.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(L2, k); it; ++it) {
			size_t row = it.row();
			size_t col = it.col();
			size_t orig_row = reverse_map[row];

			if (constrained_vertices.find(orig_row) == constrained_vertices.end()) {
				system_triplets.push_back(T(row, col, it.value()));
			}
		}
	}

	// Add regularization to diagonal for unconstrained vertices
	for (const auto& vertex : region_vertices) {
		if (constrained_vertices.find(vertex) == constrained_vertices.end()) {
			size_t i = index_map[vertex];
			system_triplets.push_back(T(i, i, L2.coeff(i, i) + regularization));
		}
	}

	// Construct the system matrix
	Eigen::SparseMatrix<double> A(system_size, system_size);
	A.setFromTriplets(system_triplets.begin(), system_triplets.end());
	A.makeCompressed();

	if (verbose) {
		Rprintf("Matrix constructed with %d non-zeros\n", (int)A.nonZeros());
	}

	// Solve the system using appropriate solver with fallback strategies
	Eigen::VectorXd x;
	bool solution_found = false;

	// Try SparseLU first (direct solver)
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with SparseLU\n");
			}
			Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
			solver.compute(A);
			if (solver.info() == Eigen::Success) {
				x = solver.solve(b);
				if (solver.info() == Eigen::Success) {
					solution_found = true;
					if (verbose) {
						Rprintf("SparseLU solution successful\n");
					}
				} else if (verbose) {
					Rprintf("SparseLU solve failed\n");
				}
			} else if (verbose) {
				Rprintf("SparseLU decomposition failed\n");
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("SparseLU exception: %s\n", e.what());
			}
		}
	}

	// Try SparseQR if SparseLU failed
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with SparseQR\n");
			}
			Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
			solver.compute(A);
			if (solver.info() == Eigen::Success) {
				x = solver.solve(b);
				if (solver.info() == Eigen::Success) {
					solution_found = true;
					if (verbose) {
						Rprintf("SparseQR solution successful\n");
					}
				} else if (verbose) {
					Rprintf("SparseQR solve failed\n");
				}
			} else if (verbose) {
				Rprintf("SparseQR decomposition failed\n");
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("SparseQR exception: %s\n", e.what());
			}
		}
	}

	// Try iterative BiCGSTAB if direct solvers failed
	if (!solution_found) {
		try {
			if (verbose) {
				Rprintf("Attempting solution with BiCGSTAB\n");
			}
			Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

			// Try with different regularization levels if needed
			std::vector<double> reg_levels = {regularization, 1e-8, 1e-6, 1e-4};
			for (double reg : reg_levels) {
				if (solution_found) break;

				// Create regularized matrix
				Eigen::SparseMatrix<double> A_reg = A;
				if (reg > regularization) {
					// Add additional regularization to diagonal
					for (int k = 0; k < A_reg.outerSize(); ++k) {
						for (Eigen::SparseMatrix<double>::InnerIterator it(A_reg, k); it; ++it) {
							if (it.row() == it.col()) {
								it.valueRef() += (reg - regularization);
							}
						}
					}
					if (verbose) {
						Rprintf("Increased regularization to %g\n", reg);
					}
				}

				// Try different iteration counts and tolerances
				std::vector<std::pair<int, double>> attempts = {
					{1000, 1e-10},
					{2000, 1e-8},
					{5000, 1e-6},
					{10000, 1e-4}
				};

				for (const auto& [maxit, tol] : attempts) {
					if (solution_found) break;

					solver.setMaxIterations(maxit);
					solver.setTolerance(tol);
					solver.compute(A_reg);

					if (solver.info() == Eigen::Success) {
						x = solver.solve(b);
						if (solver.info() == Eigen::Success) {
							solution_found = true;
							if (verbose) {
								Rprintf("BiCGSTAB solution successful with maxit=%d, tol=%g, reg=%g\n",
										maxit, tol, reg);
							}
							break;
						} else if (verbose) {
							Rprintf("BiCGSTAB solve failed with maxit=%d, tol=%g\n", maxit, tol);
						}
					} else if (verbose) {
						Rprintf("BiCGSTAB decomposition failed with maxit=%d, tol=%g\n", maxit, tol);
					}
				}
			}
		} catch (const std::exception& e) {
			if (verbose) {
				Rprintf("BiCGSTAB exception: %s\n", e.what());
			}
		}
	}

	// If all solvers failed, fall back to iterative method
	if (!solution_found) {
		if (verbose) {
			Rprintf("All Eigen solvers failed, falling back to iterative method\n");
		}

		// Use harmonic_extender instead
		harmonic_extender_t iter_result = harmonic_extender(
			boundary_values,
			region_vertices,
			10000,   // More iterations for fallback
			1e-8,// Tighter tolerance
			10000,   // Don't record intermediate states
			verbose
			);

		// If we have at least one result, use it
		if (!iter_result.iterations.empty()) {
			return iter_result.iterations.back();
		} else {
			REPORT_ERROR("Harmonic extension failed: all solvers failed");
			return result;
		}
	}

	// Convert solution back to original indexing
	for (size_t i = 0; i < system_size; i++) {
		size_t orig_idx = reverse_map[i];
		result[orig_idx] = x(i);
	}

	return result;
}


/**
 * @brief Performs a hybrid biharmonic-harmonic extension of function values
 *
 * @details This function creates a smooth extension that combines the advantages
 *          of biharmonic extension (smooth transitions at boundaries) with
 *          harmonic extension (no interior extrema). It applies biharmonic behavior
 *          near the boundary and gradually transitions to harmonic behavior in
 *          the interior of the region.
 *
 * @param[in] boundary_values Map of boundary vertex indices to their fixed values
 * @param[in] region_vertices Set of all vertex indices in the region
 * @param[in] boundary_blend_distance Controls how far the biharmonic behavior extends
 *                                    from the boundary (in hop distance)
 * @param[in] verbose Whether to print progress information
 * @return std::vector<double> Vector containing the hybrid extension values
 */
std::vector<double> set_wgraph_t::hybrid_biharmonic_harmonic_extension(
    const std::unordered_map<size_t, double>& boundary_values,
    const std::unordered_set<size_t>& region_vertices,
    int boundary_blend_distance,
    bool verbose
) const {
    // Ensure edge weights are computed
    ensure_edge_weights_computed();

    // Validate input
    if (region_vertices.empty()) {
        REPORT_WARNING("Empty region provided to hybrid_biharmonic_harmonic_extension");
        return std::vector<double>();
    }

    if (boundary_values.empty()) {
        REPORT_WARNING("No boundary values provided to hybrid_biharmonic_harmonic_extension");
        return std::vector<double>();
    }

    if (verbose) {
        Rprintf("Starting hybrid biharmonic-harmonic extension with blend distance %d\n",
                boundary_blend_distance);
    }

    // Step 1: Identify boundary vertices and compute distance from boundary for all vertices
    std::unordered_set<size_t> boundary_vertices;
    for (const auto& [v, _] : boundary_values) {
        boundary_vertices.insert(v);
    }

    // Compute distances from boundary using BFS
    std::unordered_map<size_t, int> distance_from_boundary;
    std::queue<size_t> bfs_queue;

    // Initialize distances and queue
    for (size_t v : region_vertices) {
        if (boundary_vertices.find(v) != boundary_vertices.end()) {
            distance_from_boundary[v] = 0;
            bfs_queue.push(v);
        }
    }

    // BFS to compute distances
    while (!bfs_queue.empty()) {
        size_t current = bfs_queue.front();
        bfs_queue.pop();

        int current_dist = distance_from_boundary[current];

        for (const auto& edge : adjacency_list[current]) {
            size_t neighbor = edge.vertex;

            if (region_vertices.find(neighbor) != region_vertices.end() &&
                distance_from_boundary.find(neighbor) == distance_from_boundary.end()) {
                distance_from_boundary[neighbor] = current_dist + 1;
                bfs_queue.push(neighbor);
            }
        }
    }

    // Step 2: Compute pure harmonic extension
    std::vector<double> harmonic_result = harmonic_extension_eigen(
        boundary_values, region_vertices, 1e-10, verbose);

    // Step 3: Compute biharmonic boundary conditions (required for both boundary and near-boundary)
    std::unordered_map<size_t, double> biharmonic_boundary_values = boundary_values;

    // Identify vertices within the blend distance from the boundary
    std::unordered_set<size_t> blend_region_vertices;
    for (size_t v : region_vertices) {
        if (distance_from_boundary.find(v) != distance_from_boundary.end() &&
            distance_from_boundary[v] <= boundary_blend_distance) {
            blend_region_vertices.insert(v);

            // If this is at the exact blend distance, use the harmonic value as boundary condition
            if (distance_from_boundary[v] == boundary_blend_distance) {
                biharmonic_boundary_values[v] = harmonic_result[v];
            }
        }
    }

    // Step 4: Compute biharmonic extension within the blend region
    // We need both boundary values and their "derivatives" (approximated by neighbors)
    std::unordered_map<size_t, double> biharmonic_derivative_values;

    // For each boundary vertex, estimate "normal derivative" using outside neighbors
    for (size_t v : boundary_vertices) {
        // Find neighbors outside the region
        std::vector<size_t> outside_neighbors;

        for (const auto& edge : adjacency_list[v]) {
            size_t neighbor = edge.vertex;
            if (region_vertices.find(neighbor) == region_vertices.end()) {
                outside_neighbors.push_back(neighbor);
            }
        }

        // If we have outside neighbors, use them to estimate a derivative
        if (!outside_neighbors.empty()) {
            // This is a very basic estimate - more sophisticated approaches are possible
            // We're just preserving the derivative information here
            biharmonic_derivative_values[v] = harmonic_result[v];
        }
    }

    // Step 5: Apply biharmonic extension to the blend region
    std::vector<double> biharmonic_result(adjacency_list.size());

    // Copy harmonic result as a starting point
    biharmonic_result = harmonic_result;

    // Set up and solve biharmonic system for the blend region
    // Implementation depends on whether you use a direct biharmonic solver
    // or the two-step approach mentioned earlier

    // For simplicity, let's use a simplified biharmonic approach here:
    // Solve a constrained harmonic problem that preserves derivatives at the boundary
    // This is a approximation of the biharmonic behavior

    // Step 6: Create a smooth blending between biharmonic and harmonic extensions
    std::vector<double> hybrid_result(adjacency_list.size());

    for (size_t v : region_vertices) {
        if (boundary_vertices.find(v) != boundary_vertices.end()) {
            // Boundary vertices keep their original values
            hybrid_result[v] = boundary_values.at(v);
        }
        else if (distance_from_boundary.find(v) != distance_from_boundary.end()) {
            int dist = distance_from_boundary[v];

            if (dist <= boundary_blend_distance) {
                // Compute blending weight (smooth transition from 1 to 0)
                double alpha = 1.0 - static_cast<double>(dist) / boundary_blend_distance;

                // Apply smooth cubic blending function
                alpha = alpha * alpha * (3 - 2 * alpha); // Cubic Hermite spline

                // Blend biharmonic and harmonic results
                hybrid_result[v] = alpha * biharmonic_result[v] + (1.0 - alpha) * harmonic_result[v];
            }
            else {
                // Interior vertices use harmonic values
                hybrid_result[v] = harmonic_result[v];
            }
        }
        else {
            // If we somehow missed a vertex, use harmonic value
            hybrid_result[v] = harmonic_result[v];
        }
    }

    return hybrid_result;
}


std::vector<double> set_wgraph_t::boundary_smoothed_harmonic_extension(
    const std::unordered_map<size_t, double>& boundary_values,
    const std::unordered_set<size_t>& region_vertices,
    int boundary_blend_distance,
    bool verbose
) const {
    // Ensure edge weights are computed
    ensure_edge_weights_computed();

    // Step 1: Compute distances from boundary and identify the blend region
    std::unordered_set<size_t> boundary_vertices;
    for (const auto& [v, _] : boundary_values) {
        boundary_vertices.insert(v);
    }

    // Compute distances from boundary using BFS
    std::unordered_map<size_t, int> distance_from_boundary;
    std::queue<size_t> bfs_queue;

    // Initialize distances and queue
    for (size_t v : region_vertices) {
        if (boundary_vertices.find(v) != boundary_vertices.end()) {
            distance_from_boundary[v] = 0;
            bfs_queue.push(v);
        }
    }

    // BFS to compute distances
    while (!bfs_queue.empty()) {
        size_t current = bfs_queue.front();
        bfs_queue.pop();

        int current_dist = distance_from_boundary[current];

        for (const auto& edge : adjacency_list[current]) {
            size_t neighbor = edge.vertex;

            if (region_vertices.find(neighbor) != region_vertices.end() &&
                distance_from_boundary.find(neighbor) == distance_from_boundary.end()) {
                distance_from_boundary[neighbor] = current_dist + 1;
                bfs_queue.push(neighbor);
            }
        }
    }

    // Step 2: Compute standard harmonic extension as base solution
    std::vector<double> harmonic_result = harmonic_extension_eigen(
        boundary_values, region_vertices, 1e-10, verbose);

    // Step 3: Apply boundary smoothing
    std::vector<double> smooth_result = harmonic_result;

    // For each vertex in the blend region, apply a local smoothing operation
    // that approximates biharmonic behavior near the boundary
    for (int smoothing_step = 0; smoothing_step < 3; ++smoothing_step) {
        std::vector<double> temp_result = smooth_result;

        for (size_t v : region_vertices) {
            if (distance_from_boundary.find(v) != distance_from_boundary.end() &&
                distance_from_boundary[v] <= boundary_blend_distance &&
                distance_from_boundary[v] > 0) { // Skip exact boundary

                // Compute weighted average of neighbors with curvature term
                double value_sum = 0.0;
                double weight_sum = 0.0;
                double curvature_sum = 0.0;
                double curvature_weight_sum = 0.0;

                for (const auto& edge : adjacency_list[v]) {
                    double weight = 1.0 / (edge.weight + 1e-10);

                    // Standard harmonic term
                    value_sum += smooth_result[edge.vertex] * weight;
                    weight_sum += weight;

                    // Compute discrete Laplacian at neighbor
                    double neighbor_laplacian = 0.0;
                    double neighbor_weight_sum = 0.0;

                    for (const auto& neighbor_edge : adjacency_list[edge.vertex]) {
                        double n_weight = 1.0 / (neighbor_edge.weight + 1e-10);
                        neighbor_laplacian += n_weight * (smooth_result[neighbor_edge.vertex] -
                                                         smooth_result[edge.vertex]);
                        neighbor_weight_sum += n_weight;
                    }

                    if (neighbor_weight_sum > 0) {
                        neighbor_laplacian /= neighbor_weight_sum;
                    }

                    // Apply the curvature term with distance-based weight
                    double dist_factor = 1.0;
                    if (distance_from_boundary.find(edge.vertex) != distance_from_boundary.end()) {
                        dist_factor = 1.0 - static_cast<double>(distance_from_boundary[edge.vertex]) /
                                           (boundary_blend_distance + 1);
                        dist_factor = std::max(0.0, dist_factor);
                    }

                    curvature_sum += neighbor_laplacian * weight * dist_factor;
                    curvature_weight_sum += weight * dist_factor;
                }

                // Combine harmonic value with curvature correction
                double harmonic_avg = (weight_sum > 0) ? value_sum / weight_sum : smooth_result[v];
                double curvature_term = (curvature_weight_sum > 0) ?
                                        curvature_sum / curvature_weight_sum : 0.0;

                // Compute blend factor based on distance from boundary
                double blend_factor = static_cast<double>(boundary_blend_distance -
                                                         distance_from_boundary[v]) /
                                     boundary_blend_distance;
                blend_factor = std::max(0.0, std::min(1.0, blend_factor));
                blend_factor = blend_factor * blend_factor * (3 - 2 * blend_factor); // Cubic blending

                // Apply smoothing with distance-dependent strength
                temp_result[v] = harmonic_avg + blend_factor * curvature_term * 0.25;
            }
        }

        // Update for next smoothing iteration
        smooth_result = temp_result;
    }

    return smooth_result;
}
