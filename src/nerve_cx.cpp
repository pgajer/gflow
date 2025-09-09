
#include <ANN/ANN.h>
#include <queue>
#include <algorithm>

#include "nerve_cx.hpp"
#include "error_utils.h"   // REPORT_ERROR()
#include "cpp_utils.hpp"

void nerve_complex_t::initialize_from_knn(size_t k) {
    // Compute k-neighborhoods for each vertex
    compute_k_neighborhoods(k);

    // Create 1-skeleton (graph)
    skeleton = std::make_unique<set_wgraph_t>(coordinates.size());

    // First add 0-simplices (vertices)
    for (size_t i = 0; i < coordinates.size(); ++i) {
        std::vector<size_t> vertex = {i};
        // Vertices have weight 1.0
        add_simplex(vertex, 1.0);
    }

	// Add 1-simplices (edges)
	for (size_t i = 0; i < coordinates.size(); ++i) {
		for (size_t j = i + 1; j < coordinates.size(); ++j) {
			// Check if neighborhoods intersect
			bool have_common_vertex = false;
			double min_total_distance = std::numeric_limits<double>::max();

			// Find the common neighbor that minimizes the sum of distances
			for (size_t v : knn_idx[i]) {
				if (knn_idx[j].find(v) != knn_idx[j].end()) {
					have_common_vertex = true;

					// Get distances to the common vertex
					double dist_iv = std::sqrt(knn_dist[i][v]);
					double dist_jv = std::sqrt(knn_dist[j][v]);
					double total_distance = dist_iv + dist_jv;

					// Keep track of minimum total distance
					if (total_distance < min_total_distance) {
						min_total_distance = total_distance;
					}
				}
			}

			if (have_common_vertex) {
				// Use the minimum total distance as the edge weight
				std::vector<size_t> edge = {i, j};
				add_simplex(edge, min_total_distance);

				// Add edge to 1-skeleton
				skeleton->add_edge(i, j, min_total_distance);
			}
		}
	}

	// Add 2-simplices (triangles) if max_dimension >= 2
	if (max_dimension >= 2) {

		auto find_edge = [&](size_t a, size_t b) -> double {
			std::vector<size_t> edge = {std::min(a, b), std::max(a, b)};
			auto it = simplices[1].find(edge);
			return it != simplices[1].end() ? it->second.weight : 0.0;
		};

		for (size_t i = 0; i < coordinates.size(); ++i) {
			for (size_t j = i + 1; j < coordinates.size(); ++j) {
				// Skip if edge {i,j} doesn't exist
				if (simplices[1].find({i, j}) == simplices[1].end()) continue;

				for (size_t l = j + 1; l < coordinates.size(); ++l) {
					// Skip if edges {i,l} or {j,l} don't exist
					if (simplices[1].find({i, l}) == simplices[1].end() ||
						simplices[1].find({j, l}) == simplices[1].end()) continue;

					// For 2-simplices (triangles)
					bool have_common_vertex = false;

					// Check if three neighborhoods have a common vertex
					for (size_t v : knn_idx[i]) {
						if (knn_idx[j].find(v) != knn_idx[j].end() &&
							knn_idx[l].find(v) != knn_idx[l].end()) {
							have_common_vertex = true;
							break;  // We only need to know if they intersect
						}
					}

					if (have_common_vertex) {

						// Get edge lengths from 1-simplex weights
						double d_ij = find_edge(i, j);
						double d_il = find_edge(i, l);
						double d_jl = find_edge(j, l);

						// Compute area using Heron's formula
						double s = (d_ij + d_il + d_jl) / 2.0;
						double area = std::sqrt(s * (s - d_ij) * (s - d_il) * (s - d_jl));

						// Use area as the triangle weight
						std::vector<size_t> triangle = {i, j, l};
						add_simplex(triangle, area);
					}
				}
			}
		}
	}

	// Add 3-simplices (tetrahedra) if max_dimension >= 3
	if (max_dimension >= 3) {

		// Helper function to find triangle area
		auto find_triangle_area = [&](const std::vector<size_t>& triangle) -> double {
			std::vector<size_t> sorted_triangle = triangle;
			std::sort(sorted_triangle.begin(), sorted_triangle.end());
			auto it = simplices[2].find(sorted_triangle);
			return it != simplices[2].end() ? it->second.weight : 0.0;
		};

		// Helper to find edge length
		auto find_edge_length = [&](size_t a, size_t b) -> double {
			std::vector<size_t> edge = {std::min(a, b), std::max(a, b)};
			auto it = simplices[1].find(edge);
			return it != simplices[1].end() ? it->second.weight : 0.0;
		};

		for (size_t i = 0; i < coordinates.size(); ++i) {
			for (size_t j = i + 1; j < coordinates.size(); ++j) {
				// Skip if edge {i,j} doesn't exist
				if (simplices[1].find({i, j}) == simplices[1].end()) continue;

				for (size_t l = j + 1; l < coordinates.size(); ++l) {
					// Skip if required edges don't exist
					if (simplices[1].find({i, l}) == simplices[1].end() ||
						simplices[1].find({j, l}) == simplices[1].end()) continue;

					// Skip if triangle {i,j,l} doesn't exist
					std::vector<size_t> tri_ijl = {i, j, l};
					std::sort(tri_ijl.begin(), tri_ijl.end());
					if (simplices[2].find(tri_ijl) == simplices[2].end()) continue;

					for (size_t m = l + 1; m < coordinates.size(); ++m) {
						// Skip if required edges don't exist
						if (simplices[1].find({i, m}) == simplices[1].end() ||
							simplices[1].find({j, m}) == simplices[1].end() ||
							simplices[1].find({l, m}) == simplices[1].end()) continue;

						// Skip if required triangles don't exist
						std::vector<size_t> tri_ijm = {i, j, m};
						std::vector<size_t> tri_ilm = {i, l, m};
						std::vector<size_t> tri_jlm = {j, l, m};
						std::sort(tri_ijm.begin(), tri_ijm.end());
						std::sort(tri_ilm.begin(), tri_ilm.end());
						std::sort(tri_jlm.begin(), tri_jlm.end());

						if (simplices[2].find(tri_ijm) == simplices[2].end() ||
							simplices[2].find(tri_ilm) == simplices[2].end() ||
							simplices[2].find(tri_jlm) == simplices[2].end()) continue;

						// Check if all four neighborhoods have a common vertex
						bool have_common_vertex = false;

						for (size_t v : knn_idx[i]) {
							if (knn_idx[j].find(v) != knn_idx[j].end() &&
								knn_idx[l].find(v) != knn_idx[l].end() &&
								knn_idx[m].find(v) != knn_idx[m].end()) {
								have_common_vertex = true;
								break;
							}
						}

						if (have_common_vertex) {

							// Get areas of the four faces
							double area_ijl = find_triangle_area(tri_ijl);
							double area_ijm = find_triangle_area(tri_ijm);
							double area_ilm = find_triangle_area(tri_ilm);
							double area_jlm = find_triangle_area(tri_jlm);

							// Calculate volume
							// We'll use a combination of the face areas and edge lengths
							// to approximate the volume of the tetrahedron

							// Calculate average side length
							double avg_edge = (
								find_edge_length(i, j) + find_edge_length(i, l) +
								find_edge_length(i, m) + find_edge_length(j, l) +
								find_edge_length(j, m) + find_edge_length(l, m)
								) / 6.0;

							// Calculate average face area
							double avg_area = (area_ijl + area_ijm + area_ilm + area_jlm) / 4.0;

							// Approximate height (for a regular tetrahedron)
							double height = avg_edge * std::sqrt(6.0) / 3.0;

							// Volume of a tetrahedron = (1/3) * base area * height
							double volume = avg_area * height / 3.0;

							// Use volume as the tetrahedron weight
							std::vector<size_t> tetrahedron = {i, j, l, m};
							add_simplex(tetrahedron, volume);
						}
					}
				}
			}
		}
	}

	if (max_dimension > 3) {
		REPORT_ERROR("Currently we can only build nerve comples for max_dimension not greater than 3\n");
	}
}

#if 0
bool nerve_complex_t::check_nerve_condition(
	const std::vector<size_t>& simplex, size_t k
	) const {

	if (simplex.empty()) return false;

	// Create a set to hold the intersection
	std::unordered_set<size_t> intersection;

	// Initialize with the vertices from the first neighborhood
	for (const auto& edge : knn_neighbors[simplex[0]]) {
		intersection.insert(edge.vertex);
	}

	// Find intersection with other neighborhoods
	for (size_t i = 1; i < simplex.size(); ++i) {
		std::unordered_set<size_t> temp;
		std::unordered_set<size_t> neighbor_vertices;

		// Convert edge_info_t to vertex indices
		for (const auto& edge : knn_neighbors[simplex[i]]) {
			neighbor_vertices.insert(edge.vertex);
		}

		// Find intersection
		for (size_t v : intersection) {
			if (neighbor_vertices.find(v) != neighbor_vertices.end()) {
				temp.insert(v);
			}
		}

		intersection = std::move(temp);

		if (intersection.empty()) {
			return false;
		}
	}

	return !intersection.empty();
}
#endif

void nerve_complex_t::add_simplex(const std::vector<size_t>& simplex, double weight) {
	if (simplex.empty()) return;

	size_t dim = simplex.size() - 1;
	if (dim >= simplices.size()) {
		simplices.resize(dim + 1);
	}

	// Create sorted copy of simplex vertices
	std::vector<size_t> sorted_simplex = simplex;
	std::sort(sorted_simplex.begin(), sorted_simplex.end());

	// Add the simplex to the complex
	simplex_info_t info(sorted_simplex, weight);
	simplices[dim][sorted_simplex] = info;
}

void nerve_complex_t::set_weight_function(size_t dim, weight_function_t weight_func) {
	if (dim >= weight_calculators.size()) {
		weight_calculators.resize(dim + 1);
	}
	weight_calculators[dim] = weight_func;
}

void nerve_complex_t::update_weights() {
	for (size_t dim = 0; dim < simplices.size(); ++dim) {
		if (dim >= weight_calculators.size() || !weight_calculators[dim]) continue;

		for (auto& [simplex, info] : simplices[dim]) {
			info.weight = weight_calculators[dim](simplex, coordinates, function_values);
		}
	}
}

Eigen::SparseMatrix<double> nerve_complex_t::boundary_operator(size_t p) const {
	if (p == 0 || p >= simplices.size()) {
		// Boundary of 0-simplices is empty, and p must be valid
		return Eigen::SparseMatrix<double>(0, 0);
	}

	const auto& p_simplices = simplices[p];
	const auto& p_minus_1_simplices = simplices[p-1];

	size_t rows = p_minus_1_simplices.size();
	size_t cols = p_simplices.size();

	Eigen::SparseMatrix<double> boundary(rows, cols);
	std::vector<Eigen::Triplet<double>> triplets;

	// Create maps from simplices to their indices
	std::unordered_map<std::vector<size_t>, size_t, simplex_hash_t> p_simplex_to_idx;
	std::unordered_map<std::vector<size_t>, size_t, simplex_hash_t> p_minus_1_simplex_to_idx;

	size_t idx = 0;
	for (const auto& [simplex, _] : p_simplices) {
		p_simplex_to_idx[simplex] = idx++;
	}

	idx = 0;
	for (const auto& [simplex, _] : p_minus_1_simplices) {
		p_minus_1_simplex_to_idx[simplex] = idx++;
	}

	// Compute boundary matrix entries
	for (const auto& [simplex, _] : p_simplices) {
		size_t col = p_simplex_to_idx[simplex];

		for (size_t i = 0; i < simplex.size(); ++i) {
			// Create (p-1)-simplex by removing i-th vertex
			std::vector<size_t> face;
			for (size_t j = 0; j < simplex.size(); ++j) {
				if (j != i) face.push_back(simplex[j]);
			}

			// Find the index of this face
			if (p_minus_1_simplex_to_idx.find(face) != p_minus_1_simplex_to_idx.end()) {
				size_t row = p_minus_1_simplex_to_idx[face];
				double sign = (i % 2 == 0) ? 1.0 : -1.0;

				triplets.emplace_back(row, col, sign);
			}
		}
	}

	boundary.setFromTriplets(triplets.begin(), triplets.end());
	return boundary;
}

Eigen::SparseMatrix<double> nerve_complex_t::hodge_laplacian(size_t p) const {
	if (p >= simplices.size()) {
		return Eigen::SparseMatrix<double>(0, 0);
	}

	const auto& p_simplices = simplices[p];
	size_t dim = p_simplices.size();

	Eigen::SparseMatrix<double> laplacian(dim, dim);

	if (p > 0) {
		// Down Laplacian: ∂_p^* ∂_p
		Eigen::SparseMatrix<double> boundary_p = boundary_operator(p);
		Eigen::SparseMatrix<double> down_laplacian = boundary_p.transpose() * boundary_p;

		// Apply weights
		std::vector<Eigen::Triplet<double>> down_triplets;
		for (int k = 0; k < down_laplacian.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(down_laplacian, k); it; ++it) {
				down_triplets.emplace_back(it.row(), it.col(), it.value());
			}
		}

		laplacian.setFromTriplets(down_triplets.begin(), down_triplets.end());
	}

	if (p < simplices.size() - 1) {
		// Up Laplacian: ∂_{p+1} ∂_{p+1}^*
		Eigen::SparseMatrix<double> boundary_p_plus_1 = boundary_operator(p + 1);
		Eigen::SparseMatrix<double> up_laplacian = boundary_p_plus_1 * boundary_p_plus_1.transpose();

		// Apply weights
		std::vector<Eigen::Triplet<double>> up_triplets;
		for (int k = 0; k < up_laplacian.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(up_laplacian, k); it; ++it) {
				up_triplets.emplace_back(it.row(), it.col(), it.value());
			}
		}

		// Add to existing laplacian
		if (laplacian.nonZeros() == 0) {
			laplacian.setFromTriplets(up_triplets.begin(), up_triplets.end());
		} else {
			Eigen::SparseMatrix<double> up_matrix(dim, dim);
			up_matrix.setFromTriplets(up_triplets.begin(), up_triplets.end());
			laplacian += up_matrix;
		}
	}

	return laplacian;
}

Eigen::SparseMatrix<double> nerve_complex_t::full_laplacian(
	const std::vector<double>& dim_weights
	) const {
	size_t n_vertices = num_vertices();
	Eigen::SparseMatrix<double> full_lap(n_vertices, n_vertices);

	// Start with the 0-dimensional Hodge Laplacian (graph Laplacian)
	Eigen::SparseMatrix<double> lap_0 = hodge_laplacian(0);
	full_lap = lap_0;

	// Add weighted contributions from higher-dimensional Laplacians
	for (size_t p = 1; p < std::min(dim_weights.size(), simplices.size()); ++p) {
		if (dim_weights[p] <= 0.0) continue;

		// Construct the operator B_p to incorporate p-dimensional information
		Eigen::SparseMatrix<double> B_p(n_vertices, n_vertices);
		std::vector<Eigen::Triplet<double>> triplets;

		// For each p-simplex, add its contribution
		for (const auto& [simplex, info] : simplices[p]) {
			double weight = info.weight * dim_weights[p];

			// For each vertex in the simplex
			for (size_t v_i : simplex) {
				// Add diagonal contribution
				triplets.emplace_back(v_i, v_i, weight);

				// Add off-diagonal contributions
				double off_diag = -weight / simplex.size();
				for (size_t v_j : simplex) {
					if (v_i != v_j) {
						triplets.emplace_back(v_i, v_j, off_diag);
					}
				}
			}
		}

		B_p.setFromTriplets(triplets.begin(), triplets.end());
		full_lap += B_p;
	}

	return full_lap;
}

std::vector<double> nerve_complex_t::solve_full_laplacian(
	double lambda,
	const std::vector<double>& dim_weights
	) const {
	size_t n_vertices = num_vertices();

	// Construct the full Laplacian
	Eigen::SparseMatrix<double> L_full = full_laplacian(dim_weights);

	// Add regularization
	for (int i = 0; i < L_full.outerSize(); ++i) {
		L_full.coeffRef(i, i) += lambda;
	}

	// Convert function values to Eigen vector
	Eigen::VectorXd y_vec(n_vertices);
	for (size_t i = 0; i < n_vertices; ++i) {
		y_vec(i) = function_values[i];
	}

	// Solve the system
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(L_full);

	if (solver.info() != Eigen::Success) {
		// Fall back to a more robust but slower solver
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;
		cg_solver.compute(L_full);
		Eigen::VectorXd f_vec = cg_solver.solve(lambda * y_vec);

		// Convert back to std::vector
		std::vector<double> result(n_vertices);
		for (size_t i = 0; i < n_vertices; ++i) {
			result[i] = f_vec(i);
		}
		return result;
	}

	Eigen::VectorXd f_vec = solver.solve(lambda * y_vec);

	// Convert back to std::vector
	std::vector<double> result(n_vertices);
	for (size_t i = 0; i < n_vertices; ++i) {
		result[i] = f_vec(i);
	}

	return result;
}

/**
 * @brief Compute k-neighborhoods for all vertices using ANN library
 * @param k Number of nearest neighbors (including the point itself)
 * @return Vector of neighborhoods, where each neighborhood is a vector of vertex indices
 */
void nerve_complex_t::compute_k_neighborhoods(size_t k) {
    size_t n_points = coordinates.size();
    size_t n_dims = coordinates[0].size();

    // Resize the storage for neighborhoods
    knn_idx.resize(n_points);
    knn_dist.resize(n_points);

    // Create ANNpointArray to hold coordinates
    ANNpointArray data_pts = annAllocPts(n_points, n_dims);

    // Copy coordinates into ANN data structure
    for (size_t i = 0; i < n_points; i++) {
        for (size_t j = 0; j < n_dims; j++) {
            data_pts[i][j] = coordinates[i][j];
        }
    }

    // Build kd-tree
    ANNkd_tree* kdTree = new ANNkd_tree(data_pts, n_points, n_dims);

    // Allocate arrays for kNN search results
    ANNidxArray nn_idx = new ANNidx[k];
    ANNdistArray nn_dist = new ANNdist[k];

    // Perform k-NN search for each point
    for (size_t i = 0; i < n_points; i++) {
        // Clear any existing data
        knn_idx[i].clear();
        knn_dist[i].clear();

        // ANN will include the point itself in the results
        kdTree->annkSearch(data_pts[i], k, nn_idx, nn_dist, 0.0);

        // Store the results
        for (size_t j = 0; j < k; j++) {
            size_t neighbor_idx = nn_idx[j];
            double distance = nn_dist[j];

            knn_idx[i].insert(neighbor_idx);
            knn_dist[i][neighbor_idx] = distance;
        }
    }

    // Clean up ANN resources
    delete [] nn_idx;
    delete [] nn_dist;
    delete kdTree;
    annDeallocPts(data_pts);
    annClose(); // Close ANN
}

bool nerve_complex_t::can_form_higher_simplex(
	const std::vector<size_t>& simplex1,
	const std::vector<size_t>& simplex2,
	size_t target_dim
	) const {
	if (simplex1.size() != target_dim - 1 || simplex2.size() != target_dim - 1) {
		return false;
	}

	// Check if simplices share exactly target_dim-2 vertices
	// (one vertex less than their dimension)
	std::vector<size_t> common;
	std::set_intersection(
		simplex1.begin(), simplex1.end(),
		simplex2.begin(), simplex2.end(),
		std::back_inserter(common)
		);

	return common.size() == target_dim - 2;
}

std::vector<size_t> nerve_complex_t::merge_simplices(
	const std::vector<size_t>& simplex1,
	const std::vector<size_t>& simplex2
	) const {
	std::vector<size_t> merged;
	std::set_union(
		simplex1.begin(), simplex1.end(),
		simplex2.begin(), simplex2.end(),
		std::back_inserter(merged)
		);
	return merged;
}
