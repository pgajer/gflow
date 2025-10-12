/**
 * @file riem_dcx_initialization.cpp
 * @brief Initialization routines for Riemannian simplicial complex regression
 *
 * This file implements the geometric initialization phase that constructs
 * the simplicial complex from k-NN data, computes initial densities, and
 * assembles the initial Hodge Laplacian operators.
 */

#include "riem_dcx.hpp"
#include "set_wgraph.hpp"
#include "kNN.h"
#include <algorithm>
#include <limits>
#include <unordered_set>

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
    );

/**
 * @brief Initialize Riemannian simplicial complex from k-NN structure
 *
 * This function constructs the complete geometric structure needed for
 * regression from k-nearest neighbor data. The initialization proceeds
 * through several phases that build the simplicial complex, compute
 * probability densities, construct the Riemannian metric, and assemble
 * the Hodge Laplacian operators.
 *
 * The construction follows van der Waerden's pedagogical approach: we begin
 * with the motivating observation that the k-NN graph captures local
 * geometric structure in the data. From this graph, we build a simplicial
 * complex whose topology and geometry encode both the intrinsic structure
 * of the feature space and the distribution of data points.
 *
 * PHASE 1: GEOMETRIC CONSTRUCTION
 *
 * We start by computing k-nearest neighborhoods for each data point. These
 * neighborhoods define a covering of the data space, with each neighborhood
 * serving as a local coordinate chart. The nerve theorem ensures that the
 * topology of this covering is captured by the simplicial complex whose
 * vertices correspond to data points and whose edges connect points with
 * overlapping neighborhoods.
 *
 * The edge construction uses an efficient O(nk) algorithm that examines
 * neighborhood intersections only between neighbors. For each pair of
 * vertices i,j that share neighborhood overlap, we create an edge weighted
 * by the intersection size and the minimum path length through common
 * neighbors. This weighting encodes both combinatorial (intersection size)
 * and metric (path length) information.
 *
 * Geometric edge pruning removes spurious connections that arise from
 * ambient space geometry but do not reflect intrinsic structure. We apply
 * ratio-based pruning that compares each edge length to the typical scale
 * of edges in the local neighborhood, removing edges whose length ratio
 * exceeds a threshold. This cleaning step produces a 1-skeleton that better
 * represents the intrinsic geometry.
 *
 * PHASE 2: DENSITY INITIALIZATION
 *
 * We compute initial probability densities from a reference measure on the
 * data points. The reference measure can be either uniform (counting measure)
 * or weighted by local density estimates. From this reference measure, we
 * derive vertex densities by summing measure over neighborhoods, and edge
 * densities by summing measure over neighborhood intersections.
 *
 * These densities provide the initial probability distribution that will
 * evolve during iterative refinement. The normalization ensures that vertex
 * densities sum to n and edge densities sum to n_edges, maintaining
 * interpretability as probability masses.
 *
 * PHASE 3: METRIC CONSTRUCTION
 *
 * The Riemannian metric is constructed from the density distributions. At
 * vertices (dimension 0), the metric is always diagonal with M_0 = diag(ρ_0).
 * At edges (dimension 1), the full mass matrix M_1 captures geometric
 * interactions through triple neighborhood intersections.
 *
 * For two edges sharing a vertex, their inner product is determined by the
 * density mass in the triple intersection of their endpoints' neighborhoods.
 * This construction ensures positive semidefiniteness while encoding the
 * geometric relationship between edges through their shared vertex structure.
 *
 * PHASE 4: LAPLACIAN ASSEMBLY
 *
 * The Hodge Laplacian L_0 is assembled from the boundary operator B_1 and
 * the metric tensors M_0, M_1 according to the formula L_0 = B_1 M_1^{-1} B_1^T M_0.
 * This operator encodes diffusion dynamics on the complex, with diffusion
 * rates determined by the Riemannian structure.
 *
 * The normalized Laplacian L_sym = M_0^{-1/2} L_0 M_0^{-1/2} is also computed
 * for spectral analysis and filtering operations.
 *
 * @param X Data matrix (n × d) where rows are observations and columns are features
 * @param k Number of nearest neighbors for graph construction
 * @param use_counting_measure If true, use uniform reference measure; if false,
 *                             use density-weighted measure based on local scales
 * @param density_normalization Power for density weighting (when not using counting measure).
 *                              Typical value: 1.0 gives 1/d_local weighting
 * @param max_ratio_threshold Maximum edge length ratio for geometric pruning.
 *                            Edges longer than this ratio times the local scale are removed.
 *                            Typical value: 2.0 to 3.0
 * @param threshold_percentile Percentile for computing local scale in geometric pruning.
 *                             Typical value: 0.5 (median)
 *
 * @post Simplicial complex S[0], S[1] populated with vertices and edges
 * @post Densities rho.rho[0], rho.rho[1] computed and normalized
 * @post Metric tensors g.M[0], g.M[1] constructed
 * @post Hodge Laplacian L.L[0] and normalized Laplacian L.L0_sym assembled
 * @post neighbor_sets member variable populated for subsequent use
 * @post reference_measure member variable populated
 *
 * @note This function performs substantial computation: O(nk^2) for edge
 *       construction, O(nk^3) for metric assembly. For large n or k, this
 *       initialization can take significant time.
 */
/**
 * @file riem_dcx_initialization.cpp
 * @brief Initialization routines for Riemannian simplicial complex regression
 */

#include "riem_dcx.hpp"
#include "set_wgraph.hpp"
#include "kNN.h"
#include "debug_serialization.hpp"
#include <algorithm>
#include <limits>
#include <unordered_set>

// ============================================================
// DEBUG CONTROL
// ============================================================
// Set to true to enable debug serialization
// Creates files in /tmp/gflow_debug/initialize_from_knn/
#define DEBUG_INITIALIZE_FROM_KNN true

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
);

void riem_dcx_t::initialize_from_knn(
    const spmat_t& X,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double max_ratio_threshold,
    double threshold_percentile
) {
    const size_t n_points = static_cast<size_t>(X.rows());

    // Setup debug directory if enabled
    std::string debug_dir;
    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_dir = "/tmp/gflow_debug/initialize_from_knn/";
        // Note: In production code, you'd want to ensure this directory exists
        // For debugging, create it manually: mkdir -p /tmp/gflow_debug/initialize_from_knn/
    }

    // ================================================================
    // PHASE 0: INITIALIZE DIMENSION STRUCTURE
    // ================================================================

    pmax = -1;
    pmax = 0;
    g.M.resize(1);
    g.M_solver.resize(1);
    L.B.resize(1);
    L.L.resize(1);
    g.M[0] = spmat_t();
    g.M_solver[0].reset();
    L.B[0] = spmat_t();
    L.L[0] = spmat_t();

    // ================================================================
    // PHASE 1A: K-NN COMPUTATION AND NEIGHBORHOOD CONSTRUCTION
    // ================================================================

    knn_result_t knn_result = compute_knn_from_eigen(X, k);

    std::vector<std::vector<index_t>> knn_indices(n_points, std::vector<index_t>(k));
    std::vector<std::vector<double>> knn_distances(n_points, std::vector<double>(k));

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < k; ++j) {
            const size_t offset = i * k + j;
            knn_indices[i][j] = static_cast<index_t>(knn_result.indices[offset]);
            knn_distances[i][j] = knn_result.distances[offset];
        }
    }

    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_serialization::save_knn_result(
            debug_dir + "phase_1a_knn_result.bin",
            knn_indices, knn_distances, n_points, k
        );
    }

    neighbor_sets.resize(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_serialization::save_neighbor_sets(
            debug_dir + "phase_1a_neighbor_sets.bin",
            neighbor_sets
        );
    }

    // ================================================================
    // PHASE 1B: EDGE CONSTRUCTION VIA NEIGHBORHOOD INTERSECTIONS
    // ================================================================

	// Initialize vertex_cofaces with self-loops
	vertex_cofaces.resize(n_points);
	for (index_t i = 0; i < n_points; ++i) {
		vertex_cofaces[i].reserve(k);

		// Add self-loop at position [0]
		// Convention: vertex_cofaces[i][0] is always the vertex itself
		vertex_cofaces[i].push_back(neighbor_info_t{
				i,  // vertex_index = i (self-reference)
				i,  // simplex_index = i (vertex index)
				neighbor_sets[i].size(),
				0.0,
				0.0
			});
	}

	// Build edge topology by checking ALL pairs for neighborhood overlap
	std::vector<index_t> intersection;
	intersection.reserve(k);

	for (size_t i = 0; i < n_points - 1; ++i) {  // ← Check ALL pairs
		const std::vector<index_t>& neighbors_i = knn_indices[i];

		for (size_t j = i + 1; j < n_points; ++j) {  // ← ALL j > i, not just neighbors

			// Compute neighborhood intersection
			intersection.clear();
			for (index_t v : neighbor_sets[i]) {
				if (neighbor_sets[j].find(v) != neighbor_sets[j].end()) {
					intersection.push_back(v);
				}
			}

			size_t common_count = intersection.size();
			if (common_count == 0) continue;  // No edge if neighborhoods don't overlap

			// Compute minimum path length through common neighbors
			double min_dist = std::numeric_limits<double>::max();
			for (index_t x_k : intersection) {
				auto it_i = std::find(neighbors_i.begin(), neighbors_i.end(), x_k);
				if (it_i == neighbors_i.end()) continue;
				size_t idx_i = it_i - neighbors_i.begin();
				double dist_i_k = knn_distances[i][idx_i];

				const std::vector<index_t>& neighbors_j = knn_indices[j];
				auto it_j = std::find(neighbors_j.begin(), neighbors_j.end(), x_k);
				if (it_j == neighbors_j.end()) continue;
				size_t idx_j = it_j - neighbors_j.begin();
				double dist_j_k = knn_distances[j][idx_j];

				min_dist = std::min(min_dist, dist_i_k + dist_j_k);
			}

			// Add bidirectional edges
			vertex_cofaces[i].push_back(neighbor_info_t{
					static_cast<index_t>(j),
					0,
					common_count,
					min_dist,
					0.0
				});

			vertex_cofaces[j].push_back(neighbor_info_t{
					static_cast<index_t>(i),
					0,
					common_count,
					min_dist,
					0.0
				});
		}
	}


    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_serialization::save_vertex_cofaces(
            debug_dir + "phase_1b_vertex_cofaces_pre_pruning.bin",
            vertex_cofaces
        );

        // Extract edge list for easier comparison
        std::vector<std::pair<index_t, index_t>> edges_pre_pruning;
        std::vector<double> weights_pre_pruning;
        for (size_t i = 0; i < n_points; ++i) {
            for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
                index_t j = vertex_cofaces[i][k_idx].vertex_index;
                if (j > i) {
                    edges_pre_pruning.push_back({static_cast<index_t>(i), j});
                    weights_pre_pruning.push_back(vertex_cofaces[i][k_idx].dist);
                }
            }
        }
        debug_serialization::save_edge_list(
            debug_dir + "phase_1b_edges_pre_pruning.bin",
            edges_pre_pruning, weights_pre_pruning, "PHASE_1B_PRE_PRUNING"
        );
    }

    // ================================================================
    // PHASE 1C: GEOMETRIC EDGE PRUNING
    // ================================================================

    set_wgraph_t temp_graph(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
            const auto& edge_info = vertex_cofaces[i][k_idx];
            index_t j = edge_info.vertex_index;
            if (j > i) {
                temp_graph.add_edge(i, j, edge_info.dist);
            }
        }
    }

    size_t n_edges_before_pruning = 0;
    if (DEBUG_INITIALIZE_FROM_KNN) {
        for (size_t i = 0; i < n_points; ++i) {
            n_edges_before_pruning += vertex_cofaces[i].size() - 1;
        }
        n_edges_before_pruning /= 2;

		Rprintf("=== Before geometric pruning ===\n");
		Rprintf("Edges before: %zu\n", n_edges_before_pruning);
    }

    set_wgraph_t pruned_graph = temp_graph.prune_edges_geometrically(
        max_ratio_threshold, threshold_percentile
    );

    if (DEBUG_INITIALIZE_FROM_KNN) {
		size_t n_edges_after_pruning = 0;
		for (const auto& nbrs : pruned_graph.adjacency_list) {
			n_edges_after_pruning += nbrs.size();
		}
		n_edges_after_pruning /= 2;

		Rprintf("=== After geometric pruning ===\n");
		Rprintf("Edges after: %zu\n", n_edges_after_pruning);
		Rprintf("Edges removed: %zu\n", n_edges_before_pruning - n_edges_after_pruning);

        std::vector<std::pair<index_t, index_t>> edges_post_pruning;
        std::vector<double> weights_post_pruning;
        for (size_t i = 0; i < n_points; ++i) {
            for (const auto& edge_info : pruned_graph.adjacency_list[i]) {
                index_t j = edge_info.vertex;
                if (j > i) {
                    edges_post_pruning.push_back({static_cast<index_t>(i), j});
                    weights_post_pruning.push_back(edge_info.weight);
                }
            }
        }

        debug_serialization::save_edge_list(
            debug_dir + "phase_1c_edges_post_pruning.bin",
            edges_post_pruning, weights_post_pruning, "PHASE_1C_POST_PRUNING"
        );

        debug_serialization::save_pruning_params(
            debug_dir + "phase_1c_pruning_params.bin",
            max_ratio_threshold, threshold_percentile,
            n_edges_before_pruning, edges_post_pruning.size()
        );
    }

    // ================================================================
    // PHASE 1D: FILTER vertex_cofaces IN PLACE
    // ================================================================

    for (size_t i = 0; i < n_points; ++i) {
        const auto& pruned_neighbors = pruned_graph.adjacency_list[i];

        std::unordered_set<index_t> survivors;
        for (const auto& edge_info : pruned_neighbors) {
            survivors.insert(edge_info.vertex);
        }

        auto new_end = std::remove_if(
            vertex_cofaces[i].begin() + 1,
            vertex_cofaces[i].end(),
            [&survivors](const neighbor_info_t& info) {
                return !survivors.count(info.vertex_index);
            }
        );
        vertex_cofaces[i].erase(new_end, vertex_cofaces[i].end());
    }

    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_serialization::save_vertex_cofaces(
            debug_dir + "phase_1d_vertex_cofaces_post_pruning.bin",
            vertex_cofaces
        );
    }

    // ================================================================
    // PHASE 1E: BUILD edge_registry AND ASSIGN FINAL INDICES
    // ================================================================

    std::vector<std::pair<index_t, index_t>> edge_list;
    edge_list.reserve(pruned_graph.adjacency_list.size() * k / 2);

    for (size_t i = 0; i < n_points; ++i) {
        for (const auto& edge_info : pruned_graph.adjacency_list[i]) {
            index_t j = edge_info.vertex;
            if (j > i) {
                edge_list.push_back({static_cast<index_t>(i), j});
            }
        }
    }

    const index_t n_edges = edge_list.size();
    edge_registry.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        edge_registry[e] = {edge_list[e].first, edge_list[e].second};
    }

    extend_by_one_dim(n_edges);

    for (index_t e = 0; e < n_edges; ++e) {
        const auto [i, j] = edge_registry[e];

        for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
            if (vertex_cofaces[i][k_idx].vertex_index == j) {
                vertex_cofaces[i][k_idx].simplex_index = e;
                break;
            }
        }

        for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
            if (vertex_cofaces[j][k_idx].vertex_index == i) {
                vertex_cofaces[j][k_idx].simplex_index = e;
                break;
            }
        }
    }

    if (DEBUG_INITIALIZE_FROM_KNN) {
        debug_serialization::save_vertex_cofaces(
            debug_dir + "phase_1e_final_vertex_cofaces.bin",
            vertex_cofaces
        );

        std::vector<double> final_weights;
        for (const auto& [i, j] : edge_list) {
            for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
                if (vertex_cofaces[i][k_idx].vertex_index == j) {
                    final_weights.push_back(vertex_cofaces[i][k_idx].dist);
                    break;
                }
            }
        }
        debug_serialization::save_edge_list(
            debug_dir + "phase_1e_final_edge_registry.bin",
            edge_list, final_weights, "PHASE_1E_FINAL"
        );

        // Compute and save connectivity
        std::vector<int> component_ids;
        size_t n_components = debug_serialization::compute_connected_components_from_vertex_cofaces(
            vertex_cofaces, component_ids
        );
        debug_serialization::save_connectivity(
            debug_dir + "phase_1e_connectivity.bin",
            component_ids, n_components
        );
    }

		// ================================================================
	// PHASE 1F: BUILD TRIANGLES AND POPULATE edge_cofaces
	// ================================================================

	// Initialize edge_cofaces with self-loops
	edge_cofaces.resize(n_edges);

	for (size_t e = 0; e < n_edges; ++e) {
		edge_cofaces[e].reserve(k);  // self + ~k triangles

		// Add edge self-loop at position [0]
		const auto [i, j] = edge_registry[e];

		// Find edge length from vertex_cofaces
		double edge_length = 0.0;
		for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
			if (vertex_cofaces[i][k_idx].vertex_index == j) {
				edge_length = vertex_cofaces[i][k_idx].dist;
				break;
			}
		}

		edge_cofaces[e].push_back(neighbor_info_t{
				static_cast<index_t>(i),  // Store first vertex as reference
				static_cast<index_t>(e),  // simplex_index = edge index
				0,                        // isize not used for edge self-loop
				edge_length,              // Edge length
				0.0                       // density (will be set by compute_initial_densities)
			});
	}

	// Build triangles and populate edge_cofaces simultaneously
	std::vector<std::array<index_t, 3>> triangle_list;
	triangle_list.reserve(n_points * k);

	for (size_t i = 0; i < n_points; ++i) {
		// Iterate over pairs of incident edges
		for (size_t a = 1; a < vertex_cofaces[i].size(); ++a) {
			index_t j = vertex_cofaces[i][a].vertex_index;
			// REMOVED: index_t edge_ij_idx = vertex_cofaces[i][a].simplex_index;

			for (size_t b = a + 1; b < vertex_cofaces[i].size(); ++b) {
				index_t s = vertex_cofaces[i][b].vertex_index;
				// REMOVED: index_t edge_is_idx = vertex_cofaces[i][b].simplex_index;

				// Check if edge [j,s] exists
				index_t edge_js_idx = NO_EDGE;
				for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
					if (vertex_cofaces[j][k_idx].vertex_index == s) {
						edge_js_idx = vertex_cofaces[j][k_idx].simplex_index;
						break;
					}
				}

				if (edge_js_idx == NO_EDGE) continue;  // No edge [j,s], so no triangle

				// Verify non-empty triple intersection
				bool has_intersection = false;
				for (index_t v : neighbor_sets[i]) {
					if (neighbor_sets[j].find(v) != neighbor_sets[j].end() &&
						neighbor_sets[s].find(v) != neighbor_sets[s].end()) {
						has_intersection = true;
						break;
					}
				}

				if (!has_intersection) continue;

				// Triangle exists: [i, j, s]
				std::array<index_t, 3> tri = {i, j, s};
				std::sort(tri.begin(), tri.end());
				triangle_list.push_back(tri);
			}
		}
	}

	// Remove duplicate triangles
	std::sort(triangle_list.begin(), triangle_list.end());
	triangle_list.erase(std::unique(triangle_list.begin(), triangle_list.end()),
						triangle_list.end());

	const index_t n_triangles = triangle_list.size();

	// ================================================================
	// PHASE 1G: POPULATE edge_cofaces WITH TRIANGLES
	// ================================================================

	// Now populate edge_cofaces with triangle information
	for (index_t t = 0; t < n_triangles; ++t) {
		const auto& tri = triangle_list[t];
		const index_t v0 = tri[0];
		const index_t v1 = tri[1];
		const index_t v2 = tri[2];

		// Find the three edge indices
		// Edge [v0, v1]
		index_t e01 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
			if (vertex_cofaces[v0][k].vertex_index == v1) {
				e01 = vertex_cofaces[v0][k].simplex_index;
				break;
			}
		}

		// Edge [v0, v2]
		index_t e02 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
			if (vertex_cofaces[v0][k].vertex_index == v2) {
				e02 = vertex_cofaces[v0][k].simplex_index;
				break;
			}
		}

		// Edge [v1, v2]
		index_t e12 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
			if (vertex_cofaces[v1][k].vertex_index == v2) {
				e12 = vertex_cofaces[v1][k].simplex_index;
				break;
			}
		}

		// Compute triple intersection size (for isize)
		size_t triple_isize = 0;
		for (index_t v : neighbor_sets[v0]) {
			if (neighbor_sets[v1].find(v) != neighbor_sets[v1].end() &&
				neighbor_sets[v2].find(v) != neighbor_sets[v2].end()) {
				triple_isize++;
			}
		}

		// Add triangle to edge_cofaces for each of the three edges
		// For edge [v0,v1], the third vertex is v2
		edge_cofaces[e01].push_back(neighbor_info_t{
				v2,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,    // Intersection size
				0.0,             // dist (not used for triangles)
				0.0              // density (will be set if we compute triangle densities)
			});

		// For edge [v0,v2], the third vertex is v1
		edge_cofaces[e02].push_back(neighbor_info_t{
				v1,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,
				0.0,
				0.0
			});

		// For edge [v1,v2], the third vertex is v0
		edge_cofaces[e12].push_back(neighbor_info_t{
				v0,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,
				0.0,
				0.0
			});
	}


	Rprintf("\nIn riem_dcx_t::initialize_from_knn() DEBUGGING early exit at line 609\n");
    return;


    // ================================================================
    // PHASE 2: DENSITY INITIALIZATION
    // ================================================================

    // Initialize reference measure
    initialize_reference_measure(
        knn_indices,
        knn_distances,
        use_counting_measure,
        density_normalization
    );

    // Compute initial densities from reference measure
    compute_initial_densities();

    // ================================================================
    // PHASE 3: METRIC CONSTRUCTION
    // ================================================================

    initialize_metric_from_density();

	// ================================================================
    // PHASE 4: BUILD BOUNDARY OPERATORS
    // ================================================================

	build_boundary_operator_from_edges();      // B[1]: edges → vertices

	if (edge_cofaces.size() > 0) {
		build_boundary_operator_from_triangles();  // B[2]: triangles → edges
	}

    // ================================================================
    // PHASE 5: LAPLACIAN ASSEMBLY
    // ================================================================

    assemble_operators();
}

// -------------------
#if 0
void riem_dcx_t::initialize_from_knn(
    const spmat_t& X,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double max_ratio_threshold,
    double threshold_percentile
) {
    const size_t n_points = static_cast<size_t>(X.rows());

    // ================================================================
    // PHASE 0: INITIALIZE DIMENSION STRUCTURE
    // ================================================================

    // Initialize pmax and resize dimension-indexed containers
    pmax = -1;  // Start with no dimensions

    // Add vertex dimension (dimension 0)
    pmax = 0;
    g.M.resize(1);
    g.M_solver.resize(1);
    L.B.resize(1);
    L.L.resize(1);

    // Initialize M[0] as empty (will be set by initialize_metric_from_density())
    g.M[0] = spmat_t();
    g.M_solver[0].reset();
    L.B[0] = spmat_t();  // No boundary for dimension 0
    L.L[0] = spmat_t();  // Will be assembled later

    // ================================================================
    // PHASE 1A: K-NN COMPUTATION AND NEIGHBORHOOD CONSTRUCTION
    // ================================================================

    // Compute k-nearest neighbors for all points
    knn_result_t knn_result = compute_knn_from_eigen(X, k);

    // Extract into structured format for efficient access
    std::vector<std::vector<index_t>> knn_indices(n_points, std::vector<index_t>(k));
    std::vector<std::vector<double>> knn_distances(n_points, std::vector<double>(k));

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < k; ++j) {
            const size_t offset = i * k + j;
            knn_indices[i][j] = static_cast<index_t>(knn_result.indices[offset]);
            knn_distances[i][j] = knn_result.distances[offset];
        }
    }

    // Build neighborhood sets for efficient intersection queries
	neighbor_sets.resize(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

	// ================================================================
	// PHASE 1B: EDGE CONSTRUCTION VIA NEIGHBORHOOD INTERSECTIONS
	// ================================================================

	// Initialize vertex_cofaces with self-loops
	vertex_cofaces.resize(n_points);
	for (index_t i = 0; i < n_points; ++i) {
		vertex_cofaces[i].reserve(k);

		// Add self-loop at position [0]
		// Convention: vertex_cofaces[i][0] is always the vertex itself
		vertex_cofaces[i].push_back(neighbor_info_t{
				i,  // vertex_index = i (self-reference)
				i,  // simplex_index = i (vertex index)
				neighbor_sets[i].size(),
				0.0,
				0.0
			});
	}

	// Build edge topology by checking ALL pairs for neighborhood overlap
	std::vector<index_t> intersection;
	intersection.reserve(k);

	for (size_t i = 0; i < n_points - 1; ++i) {  // ← Check ALL pairs
		const std::vector<index_t>& neighbors_i = knn_indices[i];

		for (size_t j = i + 1; j < n_points; ++j) {  // ← ALL j > i, not just neighbors

			// Compute neighborhood intersection
			intersection.clear();
			for (index_t v : neighbor_sets[i]) {
				if (neighbor_sets[j].find(v) != neighbor_sets[j].end()) {
					intersection.push_back(v);
				}
			}

			size_t common_count = intersection.size();
			if (common_count == 0) continue;  // No edge if neighborhoods don't overlap

			// Compute minimum path length through common neighbors
			double min_dist = std::numeric_limits<double>::max();
			for (index_t x_k : intersection) {
				auto it_i = std::find(neighbors_i.begin(), neighbors_i.end(), x_k);
				if (it_i == neighbors_i.end()) continue;
				size_t idx_i = it_i - neighbors_i.begin();
				double dist_i_k = knn_distances[i][idx_i];

				const std::vector<index_t>& neighbors_j = knn_indices[j];
				auto it_j = std::find(neighbors_j.begin(), neighbors_j.end(), x_k);
				if (it_j == neighbors_j.end()) continue;
				size_t idx_j = it_j - neighbors_j.begin();
				double dist_j_k = knn_distances[j][idx_j];

				min_dist = std::min(min_dist, dist_i_k + dist_j_k);
			}

			// Add bidirectional edges
			vertex_cofaces[i].push_back(neighbor_info_t{
					static_cast<index_t>(j),
					0,
					common_count,
					min_dist,
					0.0
				});

			vertex_cofaces[j].push_back(neighbor_info_t{
					static_cast<index_t>(i),
					0,
					common_count,
					min_dist,
					0.0
				});
		}
	}

	// ================================================================
	// PHASE 1C: GEOMETRIC EDGE PRUNING
	// ================================================================

	// Build temporary graph directly from vertex_cofaces
	set_wgraph_t temp_graph(n_points);
	for (size_t i = 0; i < n_points; ++i) {
		// Skip self-loop at [0], iterate over edges [1:]
		for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
			const auto& edge_info = vertex_cofaces[i][k];
			index_t j = edge_info.vertex_index;

			if (j > i) {   // Only add each edge once
				temp_graph.add_edge(i, j, edge_info.dist);
			}
		}
	}

	// Apply geometric pruning
	set_wgraph_t pruned_graph = temp_graph.prune_edges_geometrically(
		max_ratio_threshold, threshold_percentile
		);

	// ================================================================
	// PHASE 1D: FILTER vertex_cofaces IN PLACE
	// ================================================================

	for (size_t i = 0; i < n_points; ++i) {
		const auto& pruned_neighbors = pruned_graph.adjacency_list[i];

		// Build set of surviving neighbors
		std::unordered_set<index_t> survivors;
		for (const auto& edge_info : pruned_neighbors) {
			survivors.insert(edge_info.vertex);
		}

		// Remove-erase idiom: filter out non-surviving edges
		auto new_end = std::remove_if(
			vertex_cofaces[i].begin() + 1,  // Skip self-loop at [0]
			vertex_cofaces[i].end(),
			[&survivors](const neighbor_info_t& info) {
				return !survivors.count(info.vertex_index);
			}
			);
		vertex_cofaces[i].erase(new_end, vertex_cofaces[i].end());
	}

	// ================================================================
	// PHASE 1E: BUILD edge_registry AND ASSIGN FINAL INDICES
	// ================================================================

	// Build edge_registry from pruned_graph
	std::vector<std::pair<index_t, index_t>> edge_list;
	edge_list.reserve(pruned_graph.adjacency_list.size() * k / 2);   // Rough estimate

	for (size_t i = 0; i < n_points; ++i) {
		for (const auto& edge_info : pruned_graph.adjacency_list[i]) {
			index_t j = edge_info.vertex;
			if (j > i) {   // Only add each edge once
				edge_list.push_back({static_cast<index_t>(i), j});
			}
		}
	}

	const index_t n_edges = edge_list.size();
	edge_registry.resize(n_edges);
	for (index_t e = 0; e < n_edges; ++e) {
		edge_registry[e] = {edge_list[e].first, edge_list[e].second};
	}

	extend_by_one_dim(n_edges);  // Now pmax = 1, g.M resized to include M[1]

	// Update simplex_index in vertex_cofaces with final edge indices
	for (index_t e = 0; e < n_edges; ++e) {
		const auto [i, j] = edge_registry[e];

		// Update in vertex_cofaces[i]
		for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
			if (vertex_cofaces[i][k].vertex_index == j) {
				vertex_cofaces[i][k].simplex_index = e;
				break;
			}
		}

		// Update in vertex_cofaces[j]
		for (size_t k = 1; k < vertex_cofaces[j].size(); ++k) {
			if (vertex_cofaces[j][k].vertex_index == i) {
				vertex_cofaces[j][k].simplex_index = e;
				break;
			}
		}
	}

	// ================================================================
	// PHASE 1F: BUILD TRIANGLES AND POPULATE edge_cofaces
	// ================================================================

	// Initialize edge_cofaces with self-loops
	edge_cofaces.resize(n_edges);

	for (size_t e = 0; e < n_edges; ++e) {
		edge_cofaces[e].reserve(k);  // self + ~k triangles

		// Add edge self-loop at position [0]
		const auto [i, j] = edge_registry[e];

		// Find edge length from vertex_cofaces
		double edge_length = 0.0;
		for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
			if (vertex_cofaces[i][k_idx].vertex_index == j) {
				edge_length = vertex_cofaces[i][k_idx].dist;
				break;
			}
		}

		edge_cofaces[e].push_back(neighbor_info_t{
				static_cast<index_t>(i),  // Store first vertex as reference
				static_cast<index_t>(e),  // simplex_index = edge index
				0,                        // isize not used for edge self-loop
				edge_length,              // Edge length
				0.0                       // density (will be set by compute_initial_densities)
			});
	}

	// Build triangles and populate edge_cofaces simultaneously
	std::vector<std::array<index_t, 3>> triangle_list;
	triangle_list.reserve(n_points * k);

	for (size_t i = 0; i < n_points; ++i) {
		// Iterate over pairs of incident edges
		for (size_t a = 1; a < vertex_cofaces[i].size(); ++a) {
			index_t j = vertex_cofaces[i][a].vertex_index;
			// REMOVED: index_t edge_ij_idx = vertex_cofaces[i][a].simplex_index;

			for (size_t b = a + 1; b < vertex_cofaces[i].size(); ++b) {
				index_t s = vertex_cofaces[i][b].vertex_index;
				// REMOVED: index_t edge_is_idx = vertex_cofaces[i][b].simplex_index;

				// Check if edge [j,s] exists
				index_t edge_js_idx = NO_EDGE;
				for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
					if (vertex_cofaces[j][k_idx].vertex_index == s) {
						edge_js_idx = vertex_cofaces[j][k_idx].simplex_index;
						break;
					}
				}

				if (edge_js_idx == NO_EDGE) continue;  // No edge [j,s], so no triangle

				// Verify non-empty triple intersection
				bool has_intersection = false;
				for (index_t v : neighbor_sets[i]) {
					if (neighbor_sets[j].find(v) != neighbor_sets[j].end() &&
						neighbor_sets[s].find(v) != neighbor_sets[s].end()) {
						has_intersection = true;
						break;
					}
				}

				if (!has_intersection) continue;

				// Triangle exists: [i, j, s]
				std::array<index_t, 3> tri = {i, j, s};
				std::sort(tri.begin(), tri.end());
				triangle_list.push_back(tri);
			}
		}
	}

	// Remove duplicate triangles
	std::sort(triangle_list.begin(), triangle_list.end());
	triangle_list.erase(std::unique(triangle_list.begin(), triangle_list.end()),
						triangle_list.end());

	const index_t n_triangles = triangle_list.size();

	// ================================================================
	// PHASE 1G: POPULATE edge_cofaces WITH TRIANGLES
	// ================================================================

	// Now populate edge_cofaces with triangle information
	for (index_t t = 0; t < n_triangles; ++t) {
		const auto& tri = triangle_list[t];
		const index_t v0 = tri[0];
		const index_t v1 = tri[1];
		const index_t v2 = tri[2];

		// Find the three edge indices
		// Edge [v0, v1]
		index_t e01 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
			if (vertex_cofaces[v0][k].vertex_index == v1) {
				e01 = vertex_cofaces[v0][k].simplex_index;
				break;
			}
		}

		// Edge [v0, v2]
		index_t e02 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
			if (vertex_cofaces[v0][k].vertex_index == v2) {
				e02 = vertex_cofaces[v0][k].simplex_index;
				break;
			}
		}

		// Edge [v1, v2]
		index_t e12 = NO_EDGE;
		for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
			if (vertex_cofaces[v1][k].vertex_index == v2) {
				e12 = vertex_cofaces[v1][k].simplex_index;
				break;
			}
		}

		// Compute triple intersection size (for isize)
		size_t triple_isize = 0;
		for (index_t v : neighbor_sets[v0]) {
			if (neighbor_sets[v1].find(v) != neighbor_sets[v1].end() &&
				neighbor_sets[v2].find(v) != neighbor_sets[v2].end()) {
				triple_isize++;
			}
		}

		// Add triangle to edge_cofaces for each of the three edges
		// For edge [v0,v1], the third vertex is v2
		edge_cofaces[e01].push_back(neighbor_info_t{
				v2,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,    // Intersection size
				0.0,             // dist (not used for triangles)
				0.0              // density (will be set if we compute triangle densities)
			});

		// For edge [v0,v2], the third vertex is v1
		edge_cofaces[e02].push_back(neighbor_info_t{
				v1,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,
				0.0,
				0.0
			});

		// For edge [v1,v2], the third vertex is v0
		edge_cofaces[e12].push_back(neighbor_info_t{
				v0,              // vertex_index = third vertex
				t,               // simplex_index = triangle index
				triple_isize,
				0.0,
				0.0
			});
	}

    // ================================================================
    // PHASE 2: DENSITY INITIALIZATION
    // ================================================================

    // Initialize reference measure
    initialize_reference_measure(
        knn_indices,
        knn_distances,
        use_counting_measure,
        density_normalization
    );

    // Compute initial densities from reference measure
    compute_initial_densities();

    // ================================================================
    // PHASE 3: METRIC CONSTRUCTION
    // ================================================================

    initialize_metric_from_density();

	// ================================================================
    // PHASE 4: BUILD BOUNDARY OPERATORS
    // ================================================================

	build_boundary_operator_from_edges();      // B[1]: edges → vertices

	if (edge_cofaces.size() > 0) {
		build_boundary_operator_from_triangles();  // B[2]: triangles → edges
	}

    // ================================================================
    // PHASE 5: LAPLACIAN ASSEMBLY
    // ================================================================

    assemble_operators();
}
#endif
