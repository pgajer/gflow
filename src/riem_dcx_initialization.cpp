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
 * @pre pmax and dimensional structure must be initialized via init_dims()
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
    std::vector<std::unordered_set<index_t>> neighbor_sets(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

    // ================================================================
    // PHASE 1B: EDGE CONSTRUCTION VIA NEIGHBORHOOD INTERSECTIONS
    // ================================================================

    // Build adjacency list representation
    std::vector<std::vector<iknn_vertex_t>> adjacency_list(n_points);
    std::vector<index_t> intersection;
    intersection.reserve(k);

    for (size_t i = 0; i < n_points; ++i) {
        const std::vector<index_t>& neighbors_i = knn_indices[i];

        // Examine only neighbors j > i to avoid duplicates
        for (index_t j : neighbor_sets[i]) {
            if (j <= i) continue;

            // Compute neighborhood intersection
            intersection.clear();
            for (index_t v : neighbor_sets[i]) {
                if (neighbor_sets[j].find(v) != neighbor_sets[j].end()) {
                    intersection.push_back(v);
                }
            }

            size_t common_count = intersection.size();
            if (common_count == 0) continue;

            // Compute minimum path length through common neighbors
            double min_dist = std::numeric_limits<double>::max();
            for (index_t x_k : intersection) {
                // Distance from i to x_k
                auto it_i = std::find(neighbors_i.begin(), neighbors_i.end(), x_k);
                size_t idx_i = it_i - neighbors_i.begin();
                double dist_i_k = knn_distances[i][idx_i];

                // Distance from j to x_k
                const std::vector<index_t>& neighbors_j = knn_indices[j];
                auto it_j = std::find(neighbors_j.begin(), neighbors_j.end(), x_k);
                size_t idx_j = it_j - neighbors_j.begin();
                double dist_j_k = knn_distances[j][idx_j];

                min_dist = std::min(min_dist, dist_i_k + dist_j_k);
            }

            // Add bidirectional edges
            adjacency_list[i].emplace_back(iknn_vertex_t{j, common_count, min_dist});
            adjacency_list[j].emplace_back(iknn_vertex_t{i, common_count, min_dist});
        }
    }

    // Store neighborhood sets for later use in density computation
    this->neighbor_sets = std::move(neighbor_sets);

    // ================================================================
    // PHASE 1C: GEOMETRIC EDGE PRUNING
    // ================================================================

    // Extract edges into temporary structure
    struct temp_edge {
        index_t i, j;
        size_t isize;
        double dist;
    };
    std::vector<temp_edge> all_edges;

    for (size_t i = 0; i < n_points; ++i) {
        for (const auto& neighbor : adjacency_list[i]) {
            size_t j = neighbor.index;
            if (j > i) {
                all_edges.push_back({
                    static_cast<index_t>(i),
                    static_cast<index_t>(j),
                    neighbor.isize,
                    neighbor.dist
                });
            }
        }
    }

    // Apply geometric pruning via temporary graph
    set_wgraph_t temp_graph(n_points);
    for (const auto& edge : all_edges) {
        temp_graph.add_edge(edge.i, edge.j, edge.dist);
    }

    set_wgraph_t pruned_graph = temp_graph.prune_edges_geometrically(
        max_ratio_threshold, threshold_percentile
    );

    // Rebuild edge list with only retained edges
    std::vector<temp_edge> pruned_edges;
    for (const auto& edge : all_edges) {
        // Check if edge survives pruning
        bool exists = false;
        for (const auto& nbr : pruned_graph.adjacency_list[edge.i]) {
            if (nbr.vertex == edge.j) {
                exists = true;
                break;
            }
        }
        if (exists) {
            pruned_edges.push_back(edge);
        }
    }

    all_edges = std::move(pruned_edges);

    // ================================================================
    // PHASE 1D: POPULATE SIMPLEX TABLES
    // ================================================================

    // Build edge simplex table
    const index_t n_edges = all_edges.size();
    extend_by_one_dim(n_edges);

    S[1].simplex_verts.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        std::vector<index_t> verts = {all_edges[e].i, all_edges[e].j};
        S[1].simplex_verts[e] = verts;
        S[1].id_of[verts] = e;
    }

    // Build vertex simplex table
    S[0].simplex_verts.resize(n_points);
    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        S[0].simplex_verts[i] = {i};
        S[0].id_of[{i}] = i;
    }

	// ================================================================
	// PHASE 1E: POPULATE STAR TABLES
	// ================================================================

	// Build stars[0]: for each vertex, list incident edges
	stars[0].resize(n_points);
	for (index_t e = 0; e < n_edges; ++e) {
		const std::vector<index_t>& edge_verts = S[1].simplex_verts[e];
		stars[0].star_over[edge_verts[0]].push_back(e);
		stars[0].star_over[edge_verts[1]].push_back(e);
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
    // PHASE 4: LAPLACIAN ASSEMBLY
    // ================================================================

    assemble_operators();

	// ================================================================
	// PHASE 5: BUILD 2-SIMPLICES (TRIANGLES)
	// ================================================================

	std::vector<std::array<index_t, 3>> triangle_list;
	triangle_list.reserve(n_points * k);

	for (size_t i = 0; i < n_points; ++i) {
		const std::vector<index_t>& incident = stars[0].star_over[i];

		for (size_t a = 0; a < incident.size(); ++a) {
			const std::vector<index_t>& e_ij_verts = S[1].simplex_verts[incident[a]];
			index_t j = (e_ij_verts[0] == i) ? e_ij_verts[1] : e_ij_verts[0];

			for (size_t b = a + 1; b < incident.size(); ++b) {
				const std::vector<index_t>& e_is_verts = S[1].simplex_verts[incident[b]];
				index_t s = (e_is_verts[0] == i) ? e_is_verts[1] : e_is_verts[0];

				// Check if edge [j,s] exists
				std::vector<index_t> edge_js = {std::min(j,s), std::max(j,s)};
				if (S[1].id_of.find(edge_js) == S[1].id_of.end()) {
					continue;
				}

				// Verify non-empty triple intersection
				bool has_intersection = false;
				for (index_t v : this->neighbor_sets[i]) {
					if (this->neighbor_sets[j].find(v) != this->neighbor_sets[j].end() &&
						this->neighbor_sets[s].find(v) != this->neighbor_sets[s].end()) {
						has_intersection = true;
						break;
					}
				}

				if (has_intersection) {
					std::array<index_t, 3> tri = {i, j, s};
					std::sort(tri.begin(), tri.end());
					triangle_list.push_back(tri);
				}
			}
		}
	}

	// Remove duplicates
	std::sort(triangle_list.begin(), triangle_list.end());
	triangle_list.erase(std::unique(triangle_list.begin(), triangle_list.end()),
						triangle_list.end());

	const index_t n_triangles = triangle_list.size();

	if (n_triangles > 0) {
		extend_by_one_dim(n_triangles);

		S[2].simplex_verts.resize(n_triangles);
		for (index_t t = 0; t < n_triangles; ++t) {
			std::vector<index_t> tri_verts = {
				triangle_list[t][0],
				triangle_list[t][1],
				triangle_list[t][2]
			};
			S[2].simplex_verts[t] = tri_verts;
			S[2].id_of[tri_verts] = t;
		}
	}
}
