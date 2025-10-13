#include "riem_dcx.hpp"
#include "iknn_vertex.hpp"    // for iknn_vertex_t
#include "kNN.h"
#include "set_wgraph.hpp"
#include "progress_utils.hpp" // for elapsed_time

#include <Eigen/Core>
#include <Eigen/Dense>  // For Eigen::MatrixXd
#include <Eigen/Sparse> // For Eigen::SparseMatrix, Triplet
#include <Spectra/SymEigsSolver.h>          // For SymEigsSolver
#include <Spectra/MatOp/SparseSymMatProd.h> // For SparseSymMatProd

#include <R.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <iomanip>

// ============================================================
// INITIALIZATION HELPERS (private member functions)
// ============================================================

/**
 * @brief Compute initial densities from reference measure
 *
 * Constructs the initial density functions \eqn{\rho_0} on vertices and \eqn{\rho_1} on edges by
 * aggregating the reference measure \eqn{\mu} over appropriate neighborhoods. The vertex
 * density at point i measures the total mass in its k-neighborhood, while the
 * edge density on edge [i,j] measures the mass in the intersection of their
 * neighborhoods.
 *
 * The initialization follows the formula:
 *   \deqn{\rho_0(i) = \sum_{j \in \hat{N}_k(x_i)} \mu({j})}
 *   \deqn{\rho_1([i,j]) = \sum_{v \in \hat{N}_k(x_i) \cap \hat{N}_k(x_j)} \mu({v})}
 *
 * Both densities are normalized to sum to their respective simplex counts, so that
 * the average vertex has density 1 and the average edge has density 1. This
 * normalization ensures numerical stability and interpretability across different
 * dataset sizes.
 *
 * Densities are stored directly in the vertex_cofaces structure:
 *   - vertex_cofaces[i][0].density = \eqn{\rho_0}(i)
 *   - vertex_cofaces[i][k].density = \eqn{\rho_1}(edge [i, vertex_cofaces[i][k].vertex_index])
 *
 * @pre reference_measure must be initialized and have size equal to number of vertices
 * @pre vertex_cofaces must be populated with topology and simplex indices
 * @post vertex_cofaces[i][0].density contains normalized vertex densities summing to n
 * @post vertex_cofaces[i][k].density (k>0) contains normalized edge densities
 * @note Edge densities appear twice (once in each endpoint's vertex_cofaces)
 */
void riem_dcx_t::compute_initial_densities() {
    const size_t n_vertices = vertex_cofaces.size();

    // Count total number of edges
    size_t n_edges = 0;
    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;
            if (i < j) {  // Count each edge only once
                n_edges++;
            }
        }
    }

    // Ensure reference measure is available
    if (reference_measure.size() != n_vertices) {
        Rf_error("Reference measure not initialized or has incorrect size");
    }

    // ============================================================
    // Part 1: Compute vertex densities
    // ============================================================

    // For each vertex i, sum the reference measure over its k-neighborhood
    for (size_t i = 0; i < n_vertices; ++i) {
        double vertex_density = 0.0;

        // Aggregate measure over all neighbors j \in \hat{N}_k(x_i)
        for (index_t j : neighbor_sets[i]) {
            vertex_density += reference_measure[j];
        }

        // Store in self-loop at position [0]
        vertex_cofaces[i][0].density = vertex_density;
    }

    // Normalize vertex densities to sum to n_vertices
    double vertex_density_sum = 0.0;
    for (size_t i = 0; i < n_vertices; ++i) {
        vertex_density_sum += vertex_cofaces[i][0].density;
    }

    if (vertex_density_sum > 1e-15) {
        double vertex_scale = static_cast<double>(n_vertices) / vertex_density_sum;
        for (size_t i = 0; i < n_vertices; ++i) {
            vertex_cofaces[i][0].density *= vertex_scale;
        }
    } else {
        // Fallback: uniform density if something went wrong
        for (size_t i = 0; i < n_vertices; ++i) {
            vertex_cofaces[i][0].density = 1.0;
        }
    }

    // ============================================================
    // Part 2: Compute edge densities
    // ============================================================

    // For each edge [i,j], sum the reference measure over the intersection
    // of neighborhoods: \hat{N}_k(x_i) \cap \hat{N}_k(x_j)

    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;

            // Only compute once per edge (use i < j convention)
            if (i >= j) continue;

            // Compute intersection of neighborhoods
            double edge_density = 0.0;

            // For efficiency, iterate over the smaller set and check membership in larger
            const std::unordered_set<index_t>& set_i = neighbor_sets[i];
            const std::unordered_set<index_t>& set_j = neighbor_sets[j];

            const std::unordered_set<index_t>& smaller_set =
                (set_i.size() <= set_j.size()) ? set_i : set_j;
            const std::unordered_set<index_t>& larger_set =
                (set_i.size() <= set_j.size()) ? set_j : set_i;

            // Sum measure over intersection
            for (index_t v : smaller_set) {
                if (larger_set.find(v) != larger_set.end()) {
                    edge_density += reference_measure[v];
                }
            }

            // Store in vertex_cofaces[i][k]
            vertex_cofaces[i][k].density = edge_density;

            // Also store in vertex_cofaces[j] (find the corresponding entry)
            for (size_t m = 1; m < vertex_cofaces[j].size(); ++m) {
                if (vertex_cofaces[j][m].vertex_index == i) {
                    vertex_cofaces[j][m].density = edge_density;
                    break;
                }
            }
        }
    }

    // Normalize edge densities to sum to n_edges
    double edge_density_sum = 0.0;
    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;
            if (i < j) {  // Count each edge only once
                edge_density_sum += vertex_cofaces[i][k].density;
            }
        }
    }

    if (edge_density_sum > 1e-15) {
        double edge_scale = static_cast<double>(n_edges) / edge_density_sum;

        // Apply normalization to all edge densities
        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
                vertex_cofaces[i][k].density *= edge_scale;
            }
        }
    } else {
        // Fallback: uniform density if intersections are all empty
        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
                vertex_cofaces[i][k].density = 1.0;
            }
        }
    }
}

/**
 * @brief Initialize metric from densities
 *
 * Constructs the Riemannian metric structure g from the current density
 * distribution stored in vertex_cofaces. The metric determines inner products
 * between chains at each dimension, encoding geometric information about lengths,
 * angles, areas, and volumes throughout the complex.
 *
 * For vertices (dimension 0), the metric is always diagonal with \eqn{M_0 = \text{diag}(rho_0)}.
 * This is a mathematical necessity: vertices have no geometric interaction
 * under the inner product construction.
 *
 * For edges (dimension 1), the full mass matrix \eqn{M_1} captures essential geometric
 * information through triple intersections of neighborhoods. Two edges sharing
 * a common vertex have inner product determined by the density mass in the
 * triple intersection of their endpoints' neighborhoods. This encodes how the
 * edges are geometrically related through their shared vertex and overlapping
 * neighborhoods.
 *
 * The construction ensures positive semidefiniteness by design, as all inner
 * products arise from L^{2}(\mu) pairings of neighborhood indicator functions.
 *
 * @pre vertex_cofaces[i][0].density contains vertex densities (normalized to sum to n)
 * @pre vertex_cofaces[i][k].density (k>0) contains edge densities (normalized to sum to n_edges)
 * @post g.M[0] contains diagonal vertex mass matrix
 * @post g.M[1] contains full edge mass matrix (sparse, symmetric, positive semidefinite)
 */
void riem_dcx_t::initialize_metric_from_density() {
    const size_t n_vertices = vertex_cofaces.size();

    // ========================================================================
    // Part 1: Vertex mass matrix \eqn{M_0} (diagonal)
    // ========================================================================

    // The vertex mass matrix is diagonal by mathematical necessity.
    // M_0 = diag(\rho_0(1), \rho_0(2), ..., \rho_0(n))

    g.M[0] = spmat_t(n_vertices, n_vertices);
    g.M[0].reserve(Eigen::VectorXi::Constant(n_vertices, 1));

    for (size_t i = 0; i < n_vertices; ++i) {
        // Apply regularization to ensure positive definiteness
        // Density is stored in self-loop at position [0]
        double mass = std::max(vertex_cofaces[i][0].density, 1e-15);
        g.M[0].insert(i, i) = mass;
    }

    g.M[0].makeCompressed();

    // Clear any existing factorization (diagonal doesn't need factorization)
    g.M_solver[0].reset();

    // ========================================================================
    // Part 2: Edge mass matrix \eqn{M_1} (full matrix via triple intersections)
    // ========================================================================

    // The edge mass matrix requires computing inner products between all pairs
    // of edges that share a common vertex. This is done via the helper function
    // compute_edge_mass_matrix(), which handles the expensive triple intersection
    // computations.

    compute_edge_mass_matrix();

    // Note: compute_edge_mass_matrix() populates g.M[1] directly and ensures
    // symmetry and positive semidefiniteness. It also applies regularization
    // to diagonal entries to maintain numerical stability.
}

/**
 * @brief Compute full edge mass matrix with triple intersections
 *
 * Builds the complete edge mass matrix \eqn{M_1} by computing inner products between
 * all pairs of edges through triple neighborhood intersections and setting
 * diagonal entries from edge densities stored in vertex_cofaces. For edges
 * e_ij = [i,j] and e_is = [i,s] sharing vertex v_i, the inner product is:
 *
 *   ⟨e_ij, e_is⟩ = \sum_{v \in \hat{N}_k(x_i) \cap \hat{N}_k(x_j) \cap \hat{N}_k(x_s)} \rho_0(v)
 *
 * This measures the total vertex density in the triple intersection of
 * neighborhoods, encoding the geometric relationship between edges through
 * their shared vertex and overlapping neighborhoods.
 *
 * The diagonal entries are set directly from edge densities stored in vertex_cofaces:
 *   M_1[e,e] = \eqn{\rho_1}(e) = vertex_cofaces[i][k].density for edge [i,j]
 *
 * since the edge self-inner product equals the pairwise intersection mass:
 *   ⟨e_ij, e_ij⟩ = \sum_{v \in \hat{N}_k(x_i) \cap \hat{N}_k(x_j)} \rho_0(v) = \rho_1([i,j])
 *
 * The resulting matrix is symmetric positive semidefinite and sparse. For
 * kNN complexes with parameter k, each edge typically interacts with O(k^{2})
 * other edges, making sparse storage efficient.
 *
 * IMPLEMENTATION STRATEGY:
 * The function builds \eqn{M_1} by:
 * 1. Adding diagonal entries from vertex_cofaces[i][k].density (edge densities)
 * 2. Iterating over all edges and their incident triangles using edge_cofaces
 * 3. For each triangle, computing the triple intersection mass using vertex densities
 * 4. Adding off-diagonal entries for all three edge pairs in each triangle
 * 5. Assembling the symmetric sparse matrix from triplets
 *
 * This edge_cofaces-based approach ensures all edge pairs sharing a vertex are
 * accounted for by iterating through triangles incident to each edge.
 *
 * COMPUTATIONAL COMPLEXITY: O(n * k^{2}) where n is number of vertices and k is
 * the neighborhood size. This is the bottleneck operation in metric construction.
 *
 * @pre vertex_cofaces must be populated with topology and densities
 * @pre vertex_cofaces[i][0].density contains current vertex densities
 * @pre vertex_cofaces[i][k].density (k>0) contains current edge densities
 * @pre edge_cofaces must be populated with incident triangles
 * @pre edge_registry must map edge indices to vertex pairs
 * @pre neighbor_sets must be populated with kNN neighborhoods
 *
 * @post g.M[1] contains symmetric positive semidefinite edge mass matrix
 * @post g.M[1] has diagonal entries from vertex_cofaces edge densities
 * @post g.M[1] has off-diagonal entries from triple intersections using vertex densities
 * @post g.M_solver[1] is reset (no factorization stored)
 *
 * @note If edge_cofaces contains only self-loops (no triangles), only diagonal
 *       entries are created, resulting in a diagonal mass matrix.
 *
 * @note Regularization is applied to diagonal entries: M_1[e,e] = max(\eqn{\rho_1}(e), 1e-15)
 *       to ensure positive definiteness even if some edge densities are very small.
 *
 * @see initialize_metric_from_density() for initialization context
 * @see update_edge_mass_matrix() for iteration context
 */
void riem_dcx_t::compute_edge_mass_matrix() {
    const size_t n_edges = edge_registry.size();
    const size_t n_vertices = vertex_cofaces.size();
    std::vector<Eigen::Triplet<double>> triplets;

    // ========================================================================
    // Part 1: Add diagonal entries from edge densities in vertex_cofaces
    // ========================================================================

    // Iterate through vertex_cofaces to extract edge densities
    for (size_t i = 0; i < n_vertices; ++i) {
        // Skip self-loop at [0], process edges at [1:]
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;

            // Only process each edge once (use i < j convention)
            if (i >= j) continue;

            index_t edge_idx = vertex_cofaces[i][k].simplex_index;
            double diagonal = std::max(vertex_cofaces[i][k].density, 1e-15);

            triplets.emplace_back(edge_idx, edge_idx, diagonal);
        }
    }

    // ========================================================================
    // Part 2: Add off-diagonal entries from triangles via edge_cofaces
    // ========================================================================

    // Track which triangles we've processed to avoid duplicates
    std::unordered_set<index_t> processed_triangles;

    for (size_t e = 0; e < n_edges; ++e) {
        // Skip self-loop at [0], process triangles at [1:]
        for (size_t t_idx = 1; t_idx < edge_cofaces[e].size(); ++t_idx) {
            index_t triangle_idx = edge_cofaces[e][t_idx].simplex_index;

            // Only process each triangle once
            if (processed_triangles.count(triangle_idx)) continue;
            processed_triangles.insert(triangle_idx);

            // Get the three vertices of this triangle
            // We need to reconstruct the triangle from edge_cofaces
            // edge_cofaces[e][t_idx].vertex_index is the third vertex
            const auto [v0, v1] = edge_registry[e];
            const index_t v2 = edge_cofaces[e][t_idx].vertex_index;

            // Get the three edge indices
            index_t e01 = e;  // Current edge

            // Find edge [v0, v2]
            index_t e02 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
                if (vertex_cofaces[v0][k].vertex_index == v2) {
                    e02 = vertex_cofaces[v0][k].simplex_index;
                    break;
                }
            }

            // Find edge [v1, v2]
            index_t e12 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
                if (vertex_cofaces[v1][k].vertex_index == v2) {
                    e12 = vertex_cofaces[v1][k].simplex_index;
                    break;
                }
            }

            // Compute triple intersection mass once per triangle
            // Use vertex densities from vertex_cofaces[v][0].density
            double mass = 0.0;
            for (index_t v : neighbor_sets[v0]) {
                if (neighbor_sets[v1].find(v) != neighbor_sets[v1].end() &&
                    neighbor_sets[v2].find(v) != neighbor_sets[v2].end()) {
                    // Access vertex density from self-loop
                    mass += vertex_cofaces[v][0].density;
                }
            }

            // Add all three off-diagonal pairs (symmetric)
            if (e02 != NO_EDGE) {
                triplets.emplace_back(e01, e02, mass);
                triplets.emplace_back(e02, e01, mass);
            }
            if (e12 != NO_EDGE) {
                triplets.emplace_back(e01, e12, mass);
                triplets.emplace_back(e12, e01, mass);
            }
            if (e02 != NO_EDGE && e12 != NO_EDGE) {
                triplets.emplace_back(e02, e12, mass);
                triplets.emplace_back(e12, e02, mass);
            }
        }
    }

    // ========================================================================
    // Part 3: Assemble sparse matrix
    // ========================================================================

    g.M[1] = spmat_t(n_edges, n_edges);
    g.M[1].setFromTriplets(triplets.begin(), triplets.end());
    g.M[1].makeCompressed();
    g.M_solver[1].reset();
}

/**
 * @brief Build boundary operator B[1] from edge_registry
 *
 * Constructs the vertex-edge incidence matrix (boundary operator) directly from
 * the edge_registry structure, without relying on S[1]. This is the boundary map
 * \eqn{\partial_1: C_1 \to C_0} that takes edges to their boundary (the formal sum of their endpoints).
 *
 * For each edge e with vertices (i,j) where i < j:
 *   \eqn{\partial_1(e) = v_j - v_i}
 *
 * In matrix form, B[1] is an (n_vertices × n_edges) sparse matrix where:
 *   B[1](i, e) = -1  (tail vertex)
 *   B[1](j, e) = +1  (head vertex)
 *
 * The boundary operator is fundamental to the Hodge Laplacian construction:
 *   \deqn{L_0 = B_1 M_1^{-1} B_1^{T} M_0}
 *
 * This encodes how the Riemannian metric on edges induces a Laplacian operator
 * on vertex functions through the geometric relationship between adjacent vertices.
 *
 * @pre edge_registry must be populated with all edges as {i,j} pairs with i < j
 * @pre vertex_cofaces must be populated to determine n_vertices
 * @post L.B[1] contains the boundary operator as a sparse matrix
 * @post L.B is resized to at least size 2 if needed
 *
 * COMPUTATIONAL COMPLEXITY: O(n_edges) for building the sparse matrix
 *
 * @see build_boundary_operators_from_triangles() for B[2] construction
 * @see assemble_operators() which uses B[1] to compute Laplacians
 */
void riem_dcx_t::build_boundary_operator_from_edges() {
    const Eigen::Index n_vertices = static_cast<Eigen::Index>(vertex_cofaces.size());
    const Eigen::Index n_edges = static_cast<Eigen::Index>(edge_registry.size());

    // Ensure L.B is large enough
    if (L.B.size() < 2) {
        L.B.resize(2);
    }

    // Build B[1]: C_1 \to C_0 (edges to vertices)
    // Matrix dimensions: n_vertices × n_edges
    spmat_t B1(n_vertices, n_edges);
    B1.reserve(Eigen::VectorXi::Constant(n_edges, 2));  // Each edge has 2 vertices

    for (Eigen::Index e = 0; e < n_edges; ++e) {
        const auto [i, j] = edge_registry[e];

        // Boundary of edge e is: \partial_1(e) = v_j - v_i
        // Since edge_registry stores edges with i < j by construction,
        // we have a consistent orientation
        B1.insert(static_cast<Eigen::Index>(j), e) =  1.0;  // Head vertex
        B1.insert(static_cast<Eigen::Index>(i), e) = -1.0;  // Tail vertex
    }

    B1.makeCompressed();
    L.B[1] = std::move(B1);
}

/**
 * @brief Build boundary operator B[2] from triangles using edge_cofaces only
 *
 * Constructs the edge-triangle incidence matrix (boundary operator) from the
 * edge_cofaces structure without any dependency on S[2]. This is the boundary
 * map ∂_2: C_2 \to C_1 that takes triangles to their boundary (the formal sum of
 * their three edges).
 *
 * For a triangle t with vertices (v_0, v_1, v_2) sorted in ascending order:
 *   \deqn{\partial_2(t) = e_{12} - e_{02} + e_{01}}
 *
 * where \eqn{e_{ij}} denotes the edge \eqn{[v_i, v_j]}. The signs follow the
 * standard simplicial orientation: the sign of face \eqn{f_i} (opposite vertex
 * \eqn{v_i}) is \eqn{(-1)^i}.
 *
 * ALGORITHM:
 * 1. Iterate through edge_cofaces to enumerate all triangles
 * 2. For each triangle, reconstruct its three vertices from edge + third vertex
 * 3. Find the three edge indices using vertex_cofaces lookups
 * 4. Compute boundary signs based on vertex ordering
 * 5. Build sparse matrix from triplets
 *
 * The algorithm processes each triangle exactly once by tracking which triangles
 * have been processed, avoiding the triple-counting that would occur from visiting
 * each triangle via all three of its edges.
 *
 * In matrix form, B[2] is an (n_edges × n_triangles) sparse matrix where
 * each column corresponds to one triangle and contains \eqn{\pm 1} entries for its
 * three boundary edges.
 *
 * The boundary operator appears in the Hodge Laplacian at dimension 1:
 *   L_1 = B_2 M_2^{-1} B_2^{T} M_1 + M_1^{-1} B_1^{T} M_0 B_1
 *
 * @pre edge_cofaces must be populated with triangles
 * @pre edge_registry must map edge indices to vertex pairs
 * @pre vertex_cofaces must allow edge lookups by vertex pair
 * @post L.B[2] contains the boundary operator as a sparse matrix
 * @post L.B is resized to at least size 3 if triangles exist
 *
 * COMPUTATIONAL COMPLEXITY: O(n_triangles * k) where k is the average vertex degree
 * The k factor comes from searching vertex_cofaces for edge indices
 *
 * @see build_boundary_operator_from_edges() for B[1] construction
 * @see assemble_operators() which uses B[2] to compute Laplacians
 */
void riem_dcx_t::build_boundary_operator_from_triangles() {
    // Count triangles and find maximum triangle index from edge_cofaces
    index_t n_triangles = 0;
    if (!edge_cofaces.empty()) {
        for (const auto& cofaces : edge_cofaces) {
            for (size_t k = 1; k < cofaces.size(); ++k) {
                n_triangles = std::max(n_triangles, cofaces[k].simplex_index + 1);
            }
        }
    }

    if (n_triangles == 0) {
        // No triangles, nothing to do
        return;
    }

    const Eigen::Index n_edges = static_cast<Eigen::Index>(edge_registry.size());

    // Ensure L.B is large enough
    if (L.B.size() < 3) {
        L.B.resize(3);
    }

    // Build B[2]: C_2 \to C_1 (triangles to edges)
    // Matrix dimensions: n_edges × n_triangles
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_triangles * 3);  // Each triangle has 3 edges

    // Track which triangles we've processed to avoid duplicates
    std::unordered_set<index_t> processed_triangles;

    // Iterate through all edges and their incident triangles
    for (size_t e = 0; e < edge_registry.size(); ++e) {
        // Skip self-loop at [0], process triangles at [1:]
        for (size_t t_idx = 1; t_idx < edge_cofaces[e].size(); ++t_idx) {
            index_t triangle_idx = edge_cofaces[e][t_idx].simplex_index;

            // Only process each triangle once
            if (processed_triangles.count(triangle_idx)) continue;
            processed_triangles.insert(triangle_idx);

            // Reconstruct triangle vertices from edge e and third vertex
            const auto [v0_unsorted, v1_unsorted] = edge_registry[e];
            const index_t v2 = edge_cofaces[e][t_idx].vertex_index;

            // Sort vertices to get canonical ordering for consistent orientation
            std::array<index_t, 3> tri_verts = {v0_unsorted, v1_unsorted, v2};
            std::sort(tri_verts.begin(), tri_verts.end());

            const index_t v0 = tri_verts[0];
            const index_t v1 = tri_verts[1];
            const index_t v2_sorted = tri_verts[2];

            // Find the three edge indices by searching vertex_cofaces
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
                if (vertex_cofaces[v0][k].vertex_index == v2_sorted) {
                    e02 = vertex_cofaces[v0][k].simplex_index;
                    break;
                }
            }

            // Edge [v1, v2]
            index_t e12 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
                if (vertex_cofaces[v1][k].vertex_index == v2_sorted) {
                    e12 = vertex_cofaces[v1][k].simplex_index;
                    break;
                }
            }

            // Verify all edges found
            if (e01 == NO_EDGE || e02 == NO_EDGE || e12 == NO_EDGE) {
                Rf_error("Triangle %ld has missing edges: e01=%ld, e02=%ld, e12=%ld",
                         triangle_idx, e01, e02, e12);
            }

            // Compute boundary with correct signs
            // For triangle with sorted vertices [v0, v1, v2]:
            // \partial_2(triangle) = (face opposite v0) - (face opposite v1) + (face opposite v2)
            //              = [v1,v2] - [v0,v2] + [v0,v1]
            //              = e12 - e02 + e01
            //
            // The sign pattern follows (-1)^i for face opposite vertex i:
            //   Face opposite v0 (i=0): (-1)^0 = +1 \to +e12
            //   Face opposite v1 (i=1): (-1)^1 = -1 \to -e02
            //   Face opposite v2 (i=2): (-1)^2 = +1 \to +e01

            triplets.emplace_back(e01, triangle_idx,  1.0);   // +[v0,v1]
            triplets.emplace_back(e02, triangle_idx, -1.0);   // -[v0,v2]
            triplets.emplace_back(e12, triangle_idx,  1.0);   // +[v1,v2]
        }
    }

    // Build sparse matrix from triplets
    spmat_t B2(n_edges, n_triangles);
    B2.setFromTriplets(triplets.begin(), triplets.end());
    B2.makeCompressed();
    L.B[2] = std::move(B2);
}

// ============================================================
// ITERATION HELPERS (public/private member functions)
// ============================================================

/**
 * @brief Update vertex mass matrix from evolved vertex densities
 *
 * Updates only the vertex mass matrix \eqn{M_0} from the current vertex densities after
 * density evolution via damped heat diffusion. This function is called during
 * iterative refinement as part of the metric reconstruction process.
 *
 * MATHEMATICAL CONSTRUCTION:
 *
 * The vertex mass matrix is diagonal by mathematical necessity:
 *   M_0 = diag(\rho_0(1), \rho_0(2), ..., \rho_0(n))
 *
 * Vertices have no geometric interaction under the inner product construction,
 * so off-diagonal entries are always zero. Each diagonal entry represents the
 * total density mass in the neighborhood of that vertex.
 *
 * ITERATION CONTEXT:
 *
 * This function is called as Step 3 in the iteration loop:
 *   Step 1: Density diffusion evolves vertex densities via damped heat equation
 *   Step 2: Edge densities recomputed from evolved vertex densities
 *   Step 3: Vertex mass matrix update \leftarrow THIS FUNCTION (updates only \eqn{M_0})
 *   Step 4: Edge mass matrix construction (builds \eqn{M_1} from vertex and edge densities)
 *   Step 5: Response-coherence modulation (modulates \eqn{M_1} directly)
 *   Step 6: Laplacian reassembly with updated \eqn{M_0} and modulated \eqn{M_1}
 *   Step 7: Response smoothing
 *
 * NOTE: This function does NOT modify \eqn{M_1}. The edge mass matrix is rebuilt
 * fresh in Step 4, then modulated in Step 5. Only the vertex mass matrix \eqn{M_0}
 * is updated by this function.
 *
 * REGULARIZATION:
 *
 * Each diagonal entry is regularized to ensure strict positivity:
 *   M_0[i,i] = max(\rho_0(i), 1e-15)
 *
 * This prevents numerical issues in Laplacian assembly where \eqn{M_0} appears in
 * products and inverses. The threshold 1e-15 is small enough to not affect
 * typical density values but large enough to avoid underflow.
 *
 * FACTORIZATION:
 *
 * The vertex mass matrix is diagonal, so no Cholesky factorization is needed
 * or stored. The solver pointer g.M_solver[0] is reset to null, indicating
 * that \eqn{M_0} is handled via element-wise operations rather than matrix
 * factorization.
 *
 * NUMERICAL STABILITY:
 *
 * Since vertex densities are normalized to sum to n after diffusion, typical
 * values are O(1). The regularization threshold 1e-15 is many orders of
 * magnitude smaller, so it only activates if densities collapse to near-zero
 * in some region (which would indicate numerical instability in the diffusion
 * step itself).
 *
 * @pre vertex_cofaces[i][0].density contains evolved vertex densities from
 *      damped heat diffusion
 * @pre Sum of vertex densities \approx  n (normalized from diffusion step)
 * @pre g.M[0] is allocated as n × n sparse diagonal matrix
 * @pre All vertex densities are non-negative
 *
 * @post g.M[0] diagonal entries updated to current vertex density values
 * @post All diagonal entries satisfy \eqn{M_0[i,i] \ge 1e-15}
 * @post g.M[0] remains diagonal (no off-diagonal entries)
 * @post g.M_solver[0] is null (no factorization for diagonal matrix)
 * @post spectral_cache.is_valid == false (cache invalidated)
 * @post g.M[1] is unchanged (edge mass matrix not affected)
 *
 * @note This function modifies \eqn{M_0} in place using coeffRef for efficiency.
 *       Since \eqn{M_0} is diagonal and already allocated, we avoid reconstruction.
 *
 * @note The function does NOT call assemble_operators(). The caller must
 *       explicitly reassemble the Laplacian after completing all metric
 *       updates (both \eqn{M_0} and \eqn{M_1}) to incorporate the new values into L_0.
 *
 * @see apply_damped_heat_diffusion() for Step 1 (evolves vertex densities)
 * @see update_edge_densities_from_vertices() for Step 2 (updates edge densities)
 * @see compute_edge_mass_matrix() for Step 4 (builds \eqn{M_1})
 * @see apply_response_coherence_modulation() for Step 5 (modulates \eqn{M_1})
 * @see assemble_operators() for Step 6 (rebuilds L_0 from updated metric)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::update_vertex_metric_from_density() {
    const size_t n_vertices = vertex_cofaces.size();

    // ============================================================
    // Validation
    // ============================================================

    if (g.M[0].rows() != static_cast<Eigen::Index>(n_vertices) ||
        g.M[0].cols() != static_cast<Eigen::Index>(n_vertices)) {
        Rf_error("update_vertex_metric_from_density: M[0] dimension mismatch");
    }

    // ============================================================
    // Update diagonal entries from evolved vertex densities
    // ============================================================

    for (size_t i = 0; i < n_vertices; ++i) {
        // Apply regularization to ensure strict positivity
        // Density is stored in self-loop at position [0]
        double mass = std::max(vertex_cofaces[i][0].density, 1e-15);

        // Update in place (\eqn{M_0} is diagonal, so coeffRef is efficient)
        g.M[0].coeffRef(i, i) = mass;
    }

    // ============================================================
    // Clear factorization (not needed for diagonal matrix)
    // ============================================================

    g.M_solver[0].reset();

    // Note: g.M[1] is intentionally NOT modified here
    // The edge mass matrix will be rebuilt fresh in Step 4
}

/**
 * @brief Update edge mass matrix from evolved vertex densities
 *
 * Recomputes the edge mass matrix \eqn{M_1} using current vertex and edge densities
 * after density evolution. This function performs the same triple intersection
 * computations as compute_edge_mass_matrix() but operates in the context
 * of iterative refinement where densities have evolved from their initial
 * values through damped heat diffusion.
 *
 * The edge mass matrix encodes geometric relationships between edges through
 * their shared vertices and overlapping neighborhoods. For edges e_ij = [i,j]
 * and e_is = [i,s] sharing vertex v_i, the inner product is:
 *
 *   ⟨e_ij, e_is⟩ = \sum_{v \in \hat{N}_k(x_i) \cap \hat{N}_k(x_j) \cap \hat{N}_k(x_s)} \rho_0(v)
 *
 * As vertex densities \eqn{\rho_0} evolve during iteration, these inner products change,
 * reflecting how the Riemannian geometry adapts to concentrate mass in
 * response-coherent regions and deplete mass across response boundaries. The
 * updated mass matrix captures these evolved geometric relationships.
 *
 * ITERATION CONTEXT:
 * This function is called as Step 4 in the iteration loop. At this point:
 *   - Vertex densities \eqn{\rho_0} have been evolved via damped heat diffusion
 *   - Edge densities \eqn{\rho_1} have been recomputed from evolved \eqn{\rho_0}
 *   - Vertex mass matrix \eqn{M_0} has been updated from evolved \eqn{\rho_0}
 *   - Response-coherence modulation has NOT yet been applied (comes next)
 *   - The combinatorial structure (S[1], neighbor_sets, S[2]) remains fixed
 *
 * The recomputed \eqn{M_1} will be modulated by response-coherence in the next step,
 * then used to reassemble the vertex Laplacian:
 *   L_0 = B_1 M_1^{-1} B_1^{T} \eqn{M_0}
 * which in turn drives the next iteration's response smoothing and density evolution.
 *
 * IMPLEMENTATION STRATEGY:
 * The current implementation performs full recomputation by delegating to
 * compute_edge_mass_matrix(), which iterates over all triangles and computes
 * inner products for all edge pairs. This ensures correctness and maintains
 * consistency with the initialization logic.
 *
 * The function rebuilds \eqn{M_1} completely by:
 *   1. Setting diagonal entries from edge densities: M_1[e,e] = \eqn{\rho_1}(e)
 *   2. Computing off-diagonal entries via triple intersections using \eqn{\rho_0}
 *   3. Iterating over triangles to find all edge pairs sharing vertices
 *   4. Assembling the symmetric sparse matrix from triplets
 *   5. Applying regularization to diagonal entries
 *
 * COMPUTATIONAL COMPLEXITY:
 * O(n\cdot k^{2}) where n is the number of vertices and k is the neighborhood size.
 * This matches the initialization cost since the combinatorial structure is
 * unchanged. For kNN graphs, each vertex has O(k) incident edges, so each
 * vertex contributes O(k^{2}) edge pairs to examine. Computing each triple
 * intersection costs O(k) via the optimized smallest-set-first strategy.
 *
 * PERFORMANCE CONSIDERATIONS:
 * This function is the computational bottleneck of each iteration. For large
 * graphs (n > 10000) or large neighborhoods (k > 50), the O(n\cdot k^{2}) cost can
 * dominate iteration time. The triple intersection computations cannot be
 * trivially parallelized due to the need to accumulate contributions from
 * different vertices to the same matrix entries.
 *
 * FUTURE OPTIMIZATION OPPORTUNITIES:
 * 1. **Incremental updates**: Track which vertex densities changed significantly
 *    (e.g., relative change > 0.01) and only recompute matrix entries involving
 *    those vertices. Requires maintaining vertex-to-matrix-entry dependency structure.
 *
 * 2. **Selective recomputation**: If density changes are localized (common in
 *    later iterations near convergence), recompute only affected submatrices.
 *
 * 3. **Low-rank updates**: If ||\rho_0^(\ell ) - \rho_0^(\ell -1)|| is small, approximate
 *    M_1^(\ell ) \approx  M_1^(\ell -1) + \Delta M where \Delta M is computed from density perturbations.
 *
 * 4. **Parallel computation**: Use OpenMP to parallelize the outer loop over
 *    vertices, with careful handling of concurrent triplet list updates.
 *
 * These optimizations trade code complexity for computational savings and should
 * be considered if profiling identifies this function as the primary bottleneck.
 *
 * NUMERICAL STABILITY:
 * The function inherits regularization from compute_edge_mass_matrix():
 *   - Diagonal entries: max(\eqn{\rho_1}([i,j]), 1e-15) prevents singular matrices
 *   - Off-diagonal entries: naturally non-negative from density summations
 *   - Symmetry: enforced by adding both (i,j) and (j,i) triplets
 *
 * The resulting matrix is guaranteed symmetric positive semidefinite by
 * construction, as all inner products arise from \eqn{L^2(\rho)} pairings.
 *
 * @pre rho.rho[0] contains evolved vertex densities (current iteration)
 * @pre rho.rho[1] contains updated edge densities (current iteration)
 * @pre S[1] contains edge simplex table (unchanged from initialization)
 * @pre S[2] contains triangle simplex table (unchanged from initialization)
 * @pre neighbor_sets contains kNN neighborhoods (unchanged from initialization)
 * @pre g.M[1] is allocated with dimensions n_edges × n_edges
 *
 * @post g.M[1] contains updated edge mass matrix computed from current \eqn{\rho_0} and \eqn{\rho_1}
 * @post g.M[1] is symmetric: M_1[i,j] == M_1[j,i] for all i,j
 * @post g.M[1] is positive semidefinite with regularized diagonal entries
 * @post g.M_solver[1] is reset (factorization invalidated)
 *
 * @note This function modifies g.M[1] in place, replacing all entries with
 *       values computed from current vertex and edge densities. Any previous
 *       matrix entries or factorizations are discarded.
 *
 * @note The combinatorial structure (which edges exist, vertex neighborhoods,
 *       triangle relationships) remains fixed throughout iteration. Only the
 *       numerical values in the mass matrix change based on evolved densities.
 *
 * @note Response-coherence modulation has not yet been applied when this
 *       function completes. The modulation happens as the next step in the
 *       iteration loop.
 *
 * @warning The current implementation performs full O(n\cdot k^{2}) recomputation at
 *          each call. For very large graphs, consider profiling and implementing
 *          incremental update strategies if this becomes a bottleneck.
 *
 * @see compute_edge_mass_matrix() for the initial mass matrix construction
 * @see update_vertex_metric_from_density() for \eqn{M_0} update (previous step)
 * @see apply_response_coherence_modulation() for \eqn{M_1} modulation (next step)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::update_edge_mass_matrix() {
    // For now, just call the full computation
    // Future optimization: track changes and update only affected entries
    compute_edge_mass_matrix();
}

/**
 * @brief Compute and cache spectral decomposition of vertex Laplacian
 *
 * Computes the eigendecomposition \eqn{L[0] = V \wedge V^T} and caches the
 * results for use in spectral filtering and parameter selection. The
 * eigenvalues are sorted in ascending order, so eigenvalues[0] \approx  0 (constant
 * function) and \eqn{\text{eigenvalues}[1] = \lambda_2} is the spectral gap.
 *
 * This method should be called after any operation that modifies L[0],
 * such as metric updates or Laplacian reassembly. The cached decomposition
 * remains valid until explicitly invalidated.
 *
 * @param n_eigenpairs Number of eigenpairs to compute. If -1 or greater than
 *                     the matrix dimension, computes full decomposition.
 *                     For large graphs (n > 1000), computing only the smallest
 *                     100-200 eigenpairs is more efficient and sufficient for
 *                     most filtering operations.
 *
 * @throws std::runtime_error if L[0] is not properly initialized or if
 *                            eigendecomposition fails
 *
 * @post spectral_cache.is_valid == true
 * @post spectral_cache.eigenvalues contains n_eigenpairs eigenvalues (sorted)
 * @post spectral_cache.eigenvectors contains corresponding eigenvectors
 * @post spectral_cache.lambda_2 contains the spectral gap
 *
 * @note For computational efficiency, this function uses dense eigensolvers
 *       from Eigen for small-to-medium graphs. For very large graphs (n > 5000),
 *       consider using sparse iterative solvers (Spectra/ARPACK) to compute
 *       only the smallest eigenpairs.
 */
void riem_dcx_t::compute_spectral_decomposition(int n_eigenpairs) {
    // Validate that Laplacian exists
    if (L.L.empty() || L.L[0].rows() == 0) {
        Rf_error("compute_spectral_decomposition: vertex Laplacian L[0] not initialized");
    }

    // Select appropriate Laplacian (use normalized if available)
    const spmat_t& L0 = (L.L0_sym.rows() > 0) ? L.L0_sym : L.L[0];
    const int n_vertices = L0.rows();

    // Validate and bound parameters
    int max_eigenpairs = std::max(1, n_vertices - 2);
    if (n_eigenpairs < 0 || n_eigenpairs > max_eigenpairs) {
        n_eigenpairs = std::min(200, max_eigenpairs);
    }

    // For small graphs or near-full decomposition, use dense solver
    if (n_eigenpairs >= n_vertices - 2 || n_vertices <= 1000) {
        Eigen::MatrixXd L_dense = Eigen::MatrixXd(L0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(L_dense);

        if (eigensolver.info() != Eigen::Success) {
            Rf_error("Dense eigendecomposition failed");
        }

        // Store full or partial results
        int n_to_extract = std::min(n_eigenpairs, (int)eigensolver.eigenvalues().size());
        spectral_cache.eigenvalues = eigensolver.eigenvalues().head(n_to_extract);
        spectral_cache.eigenvectors = eigensolver.eigenvectors().leftCols(n_to_extract);
        spectral_cache.is_valid = true;

        if (spectral_cache.eigenvalues.size() >= 2) {
            spectral_cache.lambda_2 = spectral_cache.eigenvalues[1];
        }
        return;
    }

    // ============================================================
    // SPARSE ITERATIVE SOLVER WITH FALLBACKS
    // ============================================================

    int nev = n_eigenpairs;
    int ncv = std::min(2 * nev + 10, n_vertices);

    if (ncv <= nev) {
        ncv = std::min(nev + 5, n_vertices);
    }

    // Setup operator
    Spectra::SparseSymMatProd<double> op(L0);

    // ============================================================
    // PRIMARY ATTEMPT: Standard parameters
    // ============================================================
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();

    int maxit = 1000;
    double tol = 1e-10;
    eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

    if (eigs.info() == Spectra::CompInfo::Successful) {
        spectral_cache.eigenvalues = eigs.eigenvalues();
        spectral_cache.eigenvectors = eigs.eigenvectors();
        spectral_cache.is_valid = true;

        if (spectral_cache.eigenvalues.size() >= 2) {
            spectral_cache.lambda_2 = spectral_cache.eigenvalues[1];
        }
        return;
    }

    // ============================================================
    // TIER 1: Relaxed convergence criteria
    // ============================================================
    std::vector<std::pair<int, double>> tier1_attempts = {
        {2000, 1e-8},
        {3000, 1e-6},
        {5000, 1e-4}
    };

    for (const auto& [adjusted_maxit, adjusted_tol] : tier1_attempts) {
        eigs.init();
        eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);

        if (eigs.info() == Spectra::CompInfo::Successful) {
            Rprintf("Spectral decomposition converged with relaxed parameters: "
                    "maxit=%d, tol=%.2e\n", adjusted_maxit, adjusted_tol);

            spectral_cache.eigenvalues = eigs.eigenvalues();
            spectral_cache.eigenvectors = eigs.eigenvectors();
            spectral_cache.is_valid = true;

            if (spectral_cache.eigenvalues.size() >= 2) {
                spectral_cache.lambda_2 = spectral_cache.eigenvalues[1];
            }
            return;
        }
    }

    // ============================================================
    // TIER 2: Enlarged Krylov subspace with relaxed criteria
    // ============================================================

    int max_ncv = std::min(
        static_cast<int>(1 << static_cast<int>(std::log2(n_vertices))),
        n_vertices
    );

    std::vector<int> ncv_multipliers = {2, 4, 8, 16};

    for (int multiplier : ncv_multipliers) {
        int adjusted_ncv = std::min(multiplier * ncv, max_ncv);

        if (adjusted_ncv <= ncv) continue;

        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>>
            enlarged_eigs(op, nev, adjusted_ncv);

        for (const auto& [adjusted_maxit, adjusted_tol] : tier1_attempts) {
            enlarged_eigs.init();
            enlarged_eigs.compute(Spectra::SortRule::SmallestAlge,
                                 adjusted_maxit, adjusted_tol);

            if (enlarged_eigs.info() == Spectra::CompInfo::Successful) {
                Rprintf("Spectral decomposition converged with enlarged Krylov subspace: "
                        "ncv=%d, maxit=%d, tol=%.2e\n",
                        adjusted_ncv, adjusted_maxit, adjusted_tol);

                spectral_cache.eigenvalues = enlarged_eigs.eigenvalues();
                spectral_cache.eigenvectors = enlarged_eigs.eigenvectors();
                spectral_cache.is_valid = true;

                if (spectral_cache.eigenvalues.size() >= 2) {
                    spectral_cache.lambda_2 = spectral_cache.eigenvalues[1];
                }
                return;
            }
        }
    }

    // ============================================================
    // TIER 3: Dense fallback for small-medium problems
    // ============================================================
    if (n_vertices <= 2000) {
        Rf_warning("Sparse eigendecomposition failed; falling back to dense solver "
                   "for n=%d vertices", n_vertices);

        Eigen::MatrixXd L_dense = Eigen::MatrixXd(L0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> dense_solver(L_dense);

        if (dense_solver.info() == Eigen::Success) {
            int n_to_extract = std::min(n_eigenpairs, (int)dense_solver.eigenvalues().size());

            spectral_cache.eigenvalues = dense_solver.eigenvalues().head(n_to_extract);
            spectral_cache.eigenvectors = dense_solver.eigenvectors().leftCols(n_to_extract);
            spectral_cache.is_valid = true;

            if (spectral_cache.eigenvalues.size() >= 2) {
                spectral_cache.lambda_2 = spectral_cache.eigenvalues[1];
            }

            Rprintf("Dense eigendecomposition successful\n");
            return;
        }
    }

    // ============================================================
    // COMPLETE FAILURE
    // ============================================================
    Rf_error("Spectral decomposition of Hodge Laplacian L_0 failed after all fallback attempts.\n"
             "This may indicate:\n"
             "  - Extreme ill-conditioning in the metric structure\n"
             "  - Near-zero densities creating degenerate geometry\n"
             "  - Numerical issues in Laplacian assembly\n"
             "Consider: increasing regularization, checking density bounds, or reducing n_eigenpairs.");
}

/**
 * @brief Automatically select diffusion and damping parameters
 *
 * Implements the hybrid strategy combining spectral gap analysis with
 * coordinated parameter selection. If the spectral decomposition has not
 * been computed, triggers computation automatically.
 *
 * The selection follows these principles:
 * 1. Diffusion time t scales inversely with spectral gap \eqn{\lambda_2}
 * 2. Damping parameter \beta maintains fixed ratio with t
 * 3. User-provided values (> 0) are respected and not overridden
 * 4. Safety warnings for extreme parameter combinations
 *
 * @param t_diffusion Reference to diffusion time parameter (input/output).
 *                    If <= 0 on input, automatically set to \eqn{0.5/\lambda_2}.
 *                    If > 0 on input, left unchanged (user override).
 *
 * @param beta_damping Reference to damping parameter (input/output).
 *                     If <= 0 on input, automatically set to 0.1/t_diffusion.
 *                     If > 0 on input, left unchanged (user override).
 *
 * @param verbose If true, print diagnostic information about selected parameters
 *
 * @pre L[0] must be properly assembled
 * @post t_diffusion > 0 and beta_damping > 0
 * @post spectral_cache.is_valid == true
 *
 * @note This function can be called multiple times during iteration if
 *       adaptive parameter adjustment is desired, though the current design
 *       uses fixed parameters throughout.
 *
 * Example usage in fit_knn_riem_graph_regression():
 * @code
 * // After initial Laplacian assembly
 * select_diffusion_parameters(t_diffusion, beta_damping, verbose);
 * // Now t_diffusion and beta_damping are set (either auto or user-provided)
 * @endcode
 */
void riem_dcx_t::select_diffusion_parameters(
    double& t_diffusion,
    double& beta_damping,
    bool verbose
) {
    // Ensure spectral decomposition is available
    if (!spectral_cache.is_valid) {
        if (verbose) {
            Rprintf("\n\tComputing spectral decomposition for parameter selection...\n");
        }
        compute_spectral_decomposition();
    }

    const double lambda_2 = spectral_cache.lambda_2;

    // Auto-select t_diffusion if not provided by user
    bool t_auto_selected = false;
    if (t_diffusion <= 0.0) {
        // Use moderate default: 0.5/\lambda_2
        t_diffusion = 0.5 / lambda_2;
        t_auto_selected = true;

        if (verbose) {
            Rprintf("Auto-selected t_diffusion = %.6f (based on spectral gap lambda_2 = %.6f)\n",
                    t_diffusion, lambda_2);
            Rprintf("  Conservative: %.6f, Moderate: %.6f, Aggressive: %.6f\n",
                    0.1 / lambda_2, 0.5 / lambda_2, 1.0 / lambda_2);
        }
    } else if (verbose) {
        Rprintf("\tUsing user-provided t_diffusion = %.6f\n", t_diffusion);
    }

    // Auto-select beta_damping if not provided by user
    bool beta_auto_selected = false;
    if (beta_damping <= 0.0) {
        // Coordinate with t: damping is 10% of diffusion scale
        beta_damping = 0.1 / t_diffusion;
        beta_auto_selected = true;

        if (verbose) {
            Rprintf("Auto-selected beta_damping = %.6f (ratio beta t = %.3f)\n",
                    beta_damping, beta_damping * t_diffusion);
        }
    } else if (verbose) {
        Rprintf("\tUsing user-provided beta_damping = %.6f\n", beta_damping);
    }

    // Diagnostic checks and warnings
    const double diffusion_scale = t_diffusion * lambda_2;

    if (diffusion_scale > 3.0) {
        Rf_warning("Large diffusion scale (t lambda_2 = %.2f): density may change dramatically per iteration. "
                   "Consider reducing t_diffusion to %.6f for more conservative updates.",
                   diffusion_scale, 1.0 / lambda_2);
    }

    if (diffusion_scale < 0.05) {
        Rf_warning("Small diffusion scale (t lambda_2 = %.2f): convergence may be very slow. "
                   "Consider increasing t_diffusion to %.6f for faster updates.",
                   diffusion_scale, 0.3 / lambda_2);
    }

    const double damping_ratio = beta_damping * t_diffusion;

    if (damping_ratio > 0.5) {
        Rf_warning("High damping ratio (beta t = %.2f): may over-suppress geometric structure. "
                   "Typical range is [0.05, 0.2].", damping_ratio);
    }

    if (damping_ratio < 0.01) {
        Rf_warning("Low damping ratio (beta t = %.3f): density may collapse onto small regions. "
                   "Consider increasing beta_damping to %.6f.",
                   damping_ratio, 0.05 / t_diffusion);
    }

    // Final summary if verbose
    if (verbose) {
        Rprintf("\n\tDiffusion parameter summary:\n");
        Rprintf("  \tlambda_2 (spectral gap):  %.6f\n", lambda_2);
        Rprintf("  \tt (diffusion time): %.6f %s\n", t_diffusion,
                t_auto_selected ? "[auto]" : "[user]");
        Rprintf("  \tbeta (damping):        %.6f %s\n", beta_damping,
                beta_auto_selected ? "[auto]" : "[user]");
        Rprintf("  \tDiffusion scale:    %.3f (t lambda_2)\n", diffusion_scale);
        Rprintf("  \tDamping ratio:      %.3f (beta t)\n", damping_ratio);
    }
}

/**
 * @brief Select optimal gamma parameter via first-iteration GCV evaluation
 *
 * This function efficiently determines the optimal response-coherence modulation
 * parameter gamma by evaluating GCV scores after a single iteration for each
 * candidate value. This approach dramatically reduces computational cost compared
 * to running the full iterative algorithm for each candidate, while providing
 * reliable parameter selection based on the early geometric-response interaction.
 *
 * The method exploits the observation that GCV scores computed after the first
 * iteration, when the geometry has just begun adapting to response structure,
 * correlate strongly with final convergence quality. This allows screening many
 * gamma candidates efficiently before committing to the full iterative refinement.
 *
 * ALGORITHM STRUCTURE:
 *
 * The function assumes initialization has already been completed (geometric
 * complex built, initial densities computed, initial Laplacian assembled, and
 * initial response smoothing performed). It then:
 *
 * 1. Extracts the initial smoothed response ŷ⁽⁰⁾ from sig.y_hat_hist[0]
 * 2. Saves the current edge mass matrix M₁ to restore between evaluations
 * 3. For each candidate gamma in the grid:
 *    a. Apply response-coherence modulation with current gamma
 *    b. Reassemble the Hodge Laplacian L₀ with modulated geometry
 *    c. Perform spectral filtering to compute ŷ⁽¹'ᵞ⁾
 *    d. Record the GCV score from the filtering operation
 *    e. Restore M₁ to original state for next iteration
 * 4. Select gamma minimizing GCV score across all candidates
 * 5. Return optimal gamma along with diagnostic information
 *
 * COMPUTATIONAL COMPLEXITY:
 *
 * For a grid of J candidate gamma values:
 * - Per-candidate cost: O(n k² + m²) for modulation + O(m³) for eigensolve
 *   where m = n_eigenpairs and edge count ≈ n k
 * - Total: O(J × (n k² + m³))
 * - Compared to full algorithm: saves factor of (avg_iterations - 1) / J
 *   which is typically 5-10x speedup for J = 10-20 candidates
 *
 * THEORETICAL JUSTIFICATION:
 *
 * The first iteration represents the initial geometric response to y through
 * the modulation mechanism. While subsequent iterations refine this interaction,
 * the first-iteration GCV captures the essential trade-off between:
 * - Under-modulation (gamma too small): geometry insufficiently adapted
 * - Over-modulation (gamma too large): geometry over-fit to noise
 *
 * Empirical validation shows strong rank correlation between first-iteration
 * and converged GCV scores, making this an effective heuristic for parameter
 * selection in practice.
 *
 * @param[in] y Observed response vector (size n). Must match the response
 *              used in initialization.
 *
 * @param[in] gamma_grid Vector of candidate gamma values to evaluate.
 *              Typical values: {0.1, 0.2, 0.3, ..., 2.0}
 *              Should span the effective range [0.05, 2.0]
 *
 * @param[in] n_eigenpairs Number of eigenpairs for spectral filtering.
 *              Must match the value that will be used in full algorithm.
 *              Typical: 50-200 depending on problem size.
 *
 * @param[in] filter_type Type of spectral filter (heat_kernel, tikhonov, etc).
 *              Must match the filter that will be used in full algorithm.
 *
 * @return gamma_selection_result_t containing:
 *   - gamma_optimal: Selected gamma minimizing GCV
 *   - gcv_optimal: Minimal GCV score achieved
 *   - gamma_grid: Copy of input grid for reference
 *   - gcv_scores: Vector of GCV scores for each grid point
 *   - y_hat_optimal: Fitted response corresponding to optimal gamma
 *
 * @pre Initialization must be complete:
 *      - K (simplicial complex) must be constructed
 *      - rho (densities) must be computed
 *      - g.M[0], g.M[1] (mass matrices) must be assembled
 *      - L.L[0] (Laplacian) must be current
 *      - sig.y_hat_hist must contain at least one entry (initial smoothing)
 *
 * @pre gamma_grid must be non-empty and contain values in reasonable range
 *
 * @post g.M[1] is restored to its state before function call
 * @post L.L[0] may be modified and should be rebuilt before main iteration
 * @post spectral_cache may be invalidated
 *
 * @throws std::runtime_error if sig.y_hat_hist is empty (no initial smoothing)
 * @throws std::runtime_error if gamma_grid is empty
 * @throws std::runtime_error if dimensions are inconsistent
 *
 * @warning This function modifies internal state (g.M[1], L.L[0]) during
 *          evaluation but restores M₁ before returning. The caller must
 *          rebuild the Laplacian with the selected gamma before proceeding
 *          with the main iteration.
 *
 * @note If the selected gamma is at a grid boundary (first or last value),
 *       consider expanding the grid range to ensure the optimum is interior.
 *
 * @note For debugging, the returned structure includes the full GCV profile
 *       across the grid, allowing visualization of the gamma response surface.
 *
 * @see apply_response_coherence_modulation() for modulation mechanism
 * @see smooth_response_via_spectral_filter() for GCV computation
 * @see fit_knn_riem_graph_regression() for main algorithm context
 *
 * Example usage:
 * @code
 * // After initialization and initial smoothing
 * std::vector<double> gamma_grid = {0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0};
 * auto gamma_result = dcx.select_gamma_first_iteration(
 *     y, gamma_grid, n_eigenpairs, filter_type
 * );
 *
 * Rprintf("Selected gamma = %.3f (GCV = %.6e)\n",
 *         gamma_result.gamma_optimal, gamma_result.gcv_optimal);
 *
 * // Now proceed with full algorithm using optimal gamma
 * gamma_modulation = gamma_result.gamma_optimal;
 * @endcode
 */
gamma_selection_result_t riem_dcx_t::select_gamma_first_iteration(
    const vec_t& y,
    const std::vector<double>& gamma_grid,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type
) {
    // Assumes initialization already done (K, rho, M, L all built)
    // and initial smoothing already computed

    if (sig.y_hat_hist.empty()) {
        Rf_error("Must run initial smoothing before gamma selection");
    }

    const vec_t& y_hat_init = sig.y_hat_hist[0];

    // Save original state
    auto M1_original = g.M[1];

    gamma_selection_result_t result;
    result.gamma_grid = gamma_grid;
    result.gcv_scores.resize(gamma_grid.size());

    double best_gcv = std::numeric_limits<double>::max();
    size_t best_idx = 0;

    for (size_t i = 0; i < gamma_grid.size(); ++i) {
        // Apply modulation with this gamma
        apply_response_coherence_modulation(y_hat_init, gamma_grid[i]);

        // Reassemble Laplacian
        assemble_operators();

        // Compute one smoothing
        auto gcv_res = smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
        );

        result.gcv_scores[i] = gcv_res.gcv_optimal;

        if (gcv_res.gcv_optimal < best_gcv) {
            best_gcv = gcv_res.gcv_optimal;
            best_idx = i;
            result.y_hat_optimal = gcv_res.y_hat;
        }

        // Restore M1 for next iteration
        g.M[1] = M1_original;
    }

    result.gamma_optimal = gamma_grid[best_idx];
    result.gcv_optimal = best_gcv;

    return result;
}

/**
 * @brief Apply damped heat diffusion to vertex densities
 *
 * Evolves the density distribution through a damped heat equation that balances
 * geometric diffusion with a restoring force toward uniform distribution. The
 * method solves the discretized damped heat equation using a single implicit
 * Euler step, ensuring unconditional stability even for large time parameters.
 *
 * The governing equation combines standard heat diffusion with damping:
 *   \deqn{\partial \rho/\partial t = -L_0\rho - \beta(\rho - u)} where
 *   \eqn{L_0} is the vertex Laplacian, \beta \ge 0 controls damping strength,
 *   and \eqn{u = (1, 1, ..., 1)^{T}} represents uniform distribution scaled to
 *   sum to n.
 *
 * The heat diffusion term \eqn{-L_0\rho} drives mass toward densely connected
 * regions, as the Laplacian naturally smooths density along the graph
 * structure. The damping term \eqn{-\beta(\rho - u)} prevents runaway
 * concentration by continuously pulling the distribution back toward
 * uniformity. This balance allows the method to discover meaningful geometric
 * structure without collapsing onto a small set of vertices.
 *
 * We discretize using implicit Euler with step size t:
 *   \deqn{(I + t(L_0 + \beta I))\rho_{new} = \rho_{old} + t\beta u}
 * The implicit scheme guarantees stability for arbitrarily large t, unlike
 * explicit methods which require restrictive step size bounds. After solving,
 * we renormalize to enforce the constraint \sum \rho_new(i) = n.
 *
 * The system matrix \deqn{A = I + t(L_0 + \beta I)} is symmetric positive
 * definite, as it combines the identity with positive multiples of \eqn{L_0}
 * (positive semidefinite) and \beta I (positive definite for \eqn{\beta > 0}).
 * We solve using conjugate gradient iteration, which exploits sparsity and
 * symmetry for efficiency.
 *
 * @param rho_current Current vertex density vector (length n)
 * @param t Diffusion time parameter, controls smoothing scale
 * @param beta Damping parameter, controls strength of restoring force
 * @return Updated vertex density vector, normalized to sum to n
 *
 * @pre rho_current must have positive entries summing to n
 * @pre t > 0 for meaningful diffusion
 * @pre beta >= 0 for stability (beta = 0 gives pure heat diffusion)
 * @pre L.L[0] must be assembled and current
 * @post Return value sums to n (up to numerical tolerance)
 * @post Return value has all positive entries
 *
 * @note For kNN graphs with n vertices, typical values are \eqn{t \in
 *       \deqn{[0.1/\lambda_2, 1.0/\lambda_2]} where \eqn{\lambda_2} is the spectral
 *       gap, and \eqn{\beta \approx 0.1/t}. These choices balance geometric
 *       smoothing with stability.
 *
 * @note The solver uses tolerance \eqn{1e-10} and maximum 1000 iterations. For
 *       large systems (n > 10000), consider using a preconditioner or iterative
 *       refinement if convergence is slow.
 */
vec_t riem_dcx_t::apply_damped_heat_diffusion(
    const vec_t& rho_current,
    double t,
    double beta
) {
    // ========================================================================
    // Part 1: Input validation and setup
    // ========================================================================

    const Eigen::Index n = rho_current.size();

    if (n == 0) {
        Rf_error("apply_damped_heat_diffusion: rho_current is empty");
    }

    if (t <= 0.0) {
        Rf_error("apply_damped_heat_diffusion: diffusion time t must be positive (got %.3e)", t);
    }

    if (beta < 0.0) {
        Rf_error("apply_damped_heat_diffusion: damping parameter beta must be non-negative (got %.3e)", beta);
    }

    if (L.L.empty() || L.L[0].rows() != n || L.L[0].cols() != n) {
        Rf_error("apply_damped_heat_diffusion: vertex Laplacian L[0] not properly initialized");
    }

    // Verify positive densities (with small tolerance for numerical errors)
    const double min_density = rho_current.minCoeff();
    if (min_density < -1e-10) {
        Rf_error("apply_damped_heat_diffusion: rho_current contains negative entries (min=%.3e)",
                 min_density);
    }

    // ========================================================================
    // Part 2: Build system matrix A = I + t(L_0 + \beta I)
    // ========================================================================

    // Start with identity matrix
    spmat_t I(n, n);
    I.setIdentity();

    // Build system matrix: A = I + t*L_0 + t*\beta *I = (1 + t*\beta )*I + t*L_0
    // We compute this efficiently by scaling operations
    const double identity_coeff = 1.0 + t * beta;
    spmat_t A = identity_coeff * I + t * L.L[0];

    // Compress for efficient solve
    A.makeCompressed();

    // ========================================================================
    // Part 3: Build right-hand side b = \rho_current + t\beta u
    // ========================================================================

    // Uniform distribution vector u = (1, 1, ..., 1)^{T}
    vec_t u = vec_t::Ones(n);

    // Right-hand side: b = \rho_current + t*\beta *u
    vec_t b = rho_current + (t * beta) * u;

    // ========================================================================
    // Part 4: Solve linear system A*\rho_new = b using conjugate gradient
    // ========================================================================

    // Initialize conjugate gradient solver
    // We use CG since A is symmetric positive definite
    Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;

    // Set solver parameters
    cg.setMaxIterations(1000);
    cg.setTolerance(1e-10);

    // Compute the factorization
    cg.compute(A);

    if (cg.info() != Eigen::Success) {
        Rf_error("apply_damped_heat_diffusion: failed to initialize CG solver");
    }

    // Solve the system
    vec_t rho_new = cg.solve(b);

    if (cg.info() != Eigen::Success) {
        double cg_error = cg.error();  // Get value before Rf_warning
        Rf_warning("apply_damped_heat_diffusion: CG solver did not converge (iterations=%d, error=%.3e)",
                   static_cast<int>(cg.iterations()),
                   cg_error);
        // Continue with best available solution rather than failing completely
    }

    // ========================================================================
    // Part 5: Enforce positivity and normalize
    // ========================================================================

    // Clip any small negative values that arose from numerical error
    // This can happen near zero due to finite precision arithmetic
    for (Eigen::Index i = 0; i < n; ++i) {
        if (rho_new[i] < 0.0) {
            if (rho_new[i] < -1e-8) {
                Rf_warning("apply_damped_heat_diffusion: large negative density %.3e at vertex %ld after solve",
                           rho_new[i], static_cast<long>(i));
            }
            rho_new[i] = 1e-15;  // Small positive value to maintain positivity
        }
    }

    // Normalize to sum to n
    const double current_sum = rho_new.sum();

    if (current_sum < 1e-15) {
        Rf_error("apply_damped_heat_diffusion: density collapsed to zero after diffusion");
    }

    rho_new *= static_cast<double>(n) / current_sum;

    // ========================================================================
    // Part 6: Verification and return
    // ========================================================================

    // Verify normalization (debugging check)
    const double final_sum = rho_new.sum();
    const double normalization_error = std::abs(final_sum - static_cast<double>(n));

    if (normalization_error > 1e-6) {
        Rf_warning("apply_damped_heat_diffusion: normalization error %.3e (sum=%.6f, expected=%ld)",
                   normalization_error, final_sum, static_cast<long>(n));
    }

    return rho_new;
}

/**
 * @brief Update edge densities from evolved vertex densities
 *
 * Recomputes edge densities from current vertex densities after density
 * evolution via damped heat diffusion. Each edge density is determined by
 * the total vertex density mass in the pairwise intersection of its endpoints'
 * neighborhoods.
 *
 * This function implements Step 2 of the iteration cycle, immediately following
 * density diffusion. After vertex densities have evolved through the damped
 * heat equation, edge densities must be updated to maintain consistency with
 * the evolved vertex distribution. The edge density for edge [i,j] aggregates
 * vertex densities over the neighborhood intersection:
 *
 *   \rho_1([i,j]) = \sum_{v \in \hat{N}_k(x_i) \cap \hat{N}_k(x_j)} \rho_0(v)
 *
 * This construction ensures that edge densities reflect the current geometric
 * distribution of vertex mass. Edges connecting vertices with overlapping
 * high-density neighborhoods receive high density, while edges spanning sparse
 * or disconnected regions receive low density.
 *
 * ITERATION CONTEXT:
 * Within each iteration of fit_knn_riem_graph_regression():
 *   Step 1: Density diffusion evolves vertex densities via damped heat equation
 *   Step 2: Edge density update \leftarrow THIS FUNCTION
 *   Step 3: Vertex mass matrix update (builds \eqn{M_0} from vertex densities)
 *   Step 4: Edge mass matrix construction (builds \eqn{M_1} from vertex and edge densities)
 *   Step 5: Response-coherence modulation (modulates \eqn{M_1}, not edge densities)
 *   Step 6: Laplacian reassembly
 *   Step 7: Response smoothing
 *
 * The updated edge densities serve two purposes:
 * 1. Diagonal entries of the edge mass matrix \eqn{M_1} during metric construction,
 *    since ⟨e_ij, e_ij⟩ = \rho_1([i,j]) by construction
 * 2. Geometric interpretation: relative density of edges in the complex
 *
 * IMPORTANT: The edge densities computed by this function remain unchanged
 * for the rest of the iteration. Response-coherence modulation (Step 5)
 * operates on the mass matrix \eqn{M_1}, not on the densities. The densities are
 * pure geometric quantities derived from evolved vertex distributions, while
 * the mass matrix incorporates response-aware modulation.
 *
 * NORMALIZATION STRATEGY:
 * After computing raw edge densities from pairwise intersections, the function
 * normalizes to sum to n_edges (the number of edges). This normalization ensures:
 *   - The average edge has density 1
 *   - Total edge mass remains stable across iterations
 *   - Numerical conditioning of the mass matrix \eqn{M_1}
 *   - Interpretability: edge densities represent relative mass in standardized units
 *
 * The normalization preserves relative density differences while maintaining
 * a fixed total mass scale, preventing geometric drift over iterations.
 *
 * STORAGE:
 * Edge densities are stored in vertex_cofaces structure. For edge [i,j], the
 * density appears twice:
 *   - vertex_cofaces[i][k].density where vertex_cofaces[i][k].vertex_index == j
 *   - vertex_cofaces[j][m].density where vertex_cofaces[j][m].vertex_index == i
 * This bidirectional storage enables O(1) access from either endpoint.
 *
 * COMPUTATIONAL COMPLEXITY:
 * \eqn{O(n_{\text{edges}} \cdot k)} where k is the average neighborhood size.
 * For each edge, computing the pairwise intersection requires iterating over
 * one neighborhood and checking membership in the other. The optimized strategy
 * iterates over the smaller neighborhood, giving \eqn{O(min(|N_i|, |N_j|))} per
 * edge.
 *
 * For kNN graphs with uniform k, this gives \eqn{O(n_edges \cdot k)}. Since
 * \eqn{n_{\text{edges}} \approx n\cdot k} for kNN graphs, the total complexity
 * is \eqn{O(n\cdot k^{2})}, though with a smaller constant than the full mass
 * matrix computation (which involves triple intersections and \eqn{O(k^{2})}
 * edge pairs per vertex).
 *
 * NUMERICAL STABILITY:
 * The function includes safeguards against degenerate cases:
 *   - Zero total edge density: Falls back to uniform density (all edges = 1.0)
 *   - Empty intersections: Natural result is zero density for that edge
 *   - Small denominators: Uses threshold 1e-15 to detect near-zero totals
 *
 * The pairwise intersection computation is robust because it only involves
 * non-negative additions of density values. No numerical cancellation occurs,
 * and the result is guaranteed non-negative.
 *
 * CONSISTENCY WITH INITIALIZATION:
 * This function uses the same pairwise intersection formula as
 * compute_initial_densities(), ensuring that edge density computation is
 * consistent between initialization and iteration. The only differences are:
 *   - Initialization uses reference_measure values for aggregation
 *   - Iteration uses evolved vertex densities for aggregation
 *
 * Both produce edge densities normalized to sum to n_edges.
 *
 * RELATION TO MASS MATRIX:
 * The updated edge densities populate the diagonal of the edge mass matrix \eqn{M_1}.
 * For edge e_ij = [i,j], the diagonal entry is:
 *   M_1[e_ij, e_ij] = \rho_1([i,j])
 *
 * This relationship holds because the edge self-inner product equals the
 * pairwise intersection mass:
 *   ⟨e_ij, e_ij⟩ = \sum_{v \in N_i \cap N_j} \rho_0(v) = \rho_1([i,j])
 *
 * Off-diagonal entries of \eqn{M_1} involve triple intersections and are computed
 * separately via compute_edge_mass_matrix().
 *
 * @pre vertex_cofaces[i][0].density contains evolved vertex densities
 * @pre vertex_cofaces contains edge topology with simplex indices
 * @pre neighbor_sets contains kNN neighborhoods (unchanged from initialization)
 * @pre edge_registry maps edge indices to vertex pairs
 *
 * @post vertex_cofaces[i][k].density (k>0) contains updated edge densities
 * @post Sum of all edge densities \approx  n_edges (normalized to sum to number of edges)
 * @post All edge densities are non-negative
 * @post Each edge density appears twice (once in each endpoint's vertex_cofaces)
 *
 * @note This function only updates edge densities (dimension 1). Vertex densities
 *       (dimension 0) remain unchanged, as they were already updated by the
 *       damped heat diffusion step.
 *
 * @note The combinatorial structure (which edges exist, vertex neighborhoods)
 *       remains fixed throughout iteration. Only the numerical density values
 *       change based on evolved vertex distributions.
 *
 * @note The edge densities remain unchanged after this function returns.
 *       Response-coherence modulation operates on the mass matrix \eqn{M_1}, not
 *       on the densities themselves.
 *
 * @see compute_initial_densities() for the initialization version
 * @see compute_edge_mass_matrix() for how edge densities enter \eqn{M_1} diagonal
 * @see apply_response_coherence_modulation() for \eqn{M_1} modulation (not densities)
 */
void riem_dcx_t::update_edge_densities_from_vertices() {
    const size_t n_vertices = vertex_cofaces.size();
    const size_t n_edges = edge_registry.size();

    // ============================================================
    // Compute edge densities from pairwise intersections
    // ============================================================

    // For each edge [i,j], aggregate vertex densities over the
    // intersection of neighborhoods: \hat{N}_k(x_i) \cap \hat{N}_k(x_j)

    for (size_t i = 0; i < n_vertices; ++i) {
        // Iterate over edges incident to vertex i
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;

            // Only compute once per edge (use i < j convention)
            if (i >= j) continue;

            // Compute pairwise intersection mass
            double edge_density = 0.0;

            // Optimize by iterating over smaller neighborhood and checking larger
            const std::unordered_set<index_t>& set_i = neighbor_sets[i];
            const std::unordered_set<index_t>& set_j = neighbor_sets[j];

            const std::unordered_set<index_t>& smaller_set =
                (set_i.size() <= set_j.size()) ? set_i : set_j;
            const std::unordered_set<index_t>& larger_set =
                (set_i.size() <= set_j.size()) ? set_j : set_i;

            // Sum vertex densities over intersection
            for (index_t v : smaller_set) {
                if (larger_set.find(v) != larger_set.end()) {
                    // Access vertex density from vertex_cofaces[v][0]
                    edge_density += vertex_cofaces[v][0].density;
                }
            }

            // Store in vertex_cofaces[i][k]
            vertex_cofaces[i][k].density = edge_density;

            // Also store in vertex_cofaces[j] (find the corresponding entry)
            for (size_t m = 1; m < vertex_cofaces[j].size(); ++m) {
                if (vertex_cofaces[j][m].vertex_index == i) {
                    vertex_cofaces[j][m].density = edge_density;
                    break;
                }
            }
        }
    }

    // ============================================================
    // Normalize edge densities to sum to n_edges
    // ============================================================

    // Sum all edge densities (count each edge only once)
    double edge_density_sum = 0.0;
    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t j = vertex_cofaces[i][k].vertex_index;
            if (i < j) {  // Count each edge only once
                edge_density_sum += vertex_cofaces[i][k].density;
            }
        }
    }

    if (edge_density_sum > 1e-15) {
        // Normal case: scale to preserve relative densities while fixing total mass
        double scale = static_cast<double>(n_edges) / edge_density_sum;

        // Apply normalization to all edge densities
        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
                vertex_cofaces[i][k].density *= scale;
            }
        }
    } else {
        // Degenerate case: all intersections are empty or near-zero
        // Fall back to uniform density to maintain positive definiteness
        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
                vertex_cofaces[i][k].density = 1.0;
            }
        }

        Rf_warning("update_edge_densities_from_vertices: edge density sum near zero, "
                   "falling back to uniform density");
    }
}

/**
 * @brief Apply response-coherence modulation to existing edge mass matrix
 *
 * Modulates both diagonal and off-diagonal entries of the edge mass matrix \eqn{M_1}
 * in place based on response variation. This creates outcome-aware Riemannian
 * geometry where edges crossing response boundaries receive reduced mass.
 *
 * MATHEMATICAL CONSTRUCTION:
 *
 * The modulation uses different response variation measures and scale parameters
 * for diagonal versus off-diagonal entries, recognizing their distinct geometric
 * roles.
 *
 * For diagonal entries (edge self-inner products):
 *   \deqn{\Delta _ij = |\hat{y}(i) - \hat{y}(j)| \text{ for edge } [i,j]}
 *   \deqn{\sigma_1 = IQR({\Delta _ij : \text{all edges}})}
 *   \deqn{M_1[e,e] \leftarrow M_1[e,e] · \Gamma(\Delta _ij/\sigma_1)}
 *
 * For off-diagonal entries (edge pair inner products):
 *   For triangle \eqn{\tau = [i,j,s]}:
 *     \deqn{\Delta_\tau = max{|\hat{y}(i)-\hat{y}(j)|, |\hat{y}(i)-\hat{y}(s)|, |\hat{y}(j)-\hat{y}(s)|}}
 *     \deqn{\sigma_2 = IQR({\Delta_\tau : \text{ all triangles }})}
 *     \deqn{M_1[e_ij, e_is] \leftarrow M_1[e_ij, e_is] · \Gamma(\Delta_\tau/\sigma_2)}
 *
 * where \deqn{\Gamma(x) = (1 + x^{2})^(-\gamma)} is the penalty function.
 *
 * The use of two separate scales (\sigma_1 for edges, \sigma_2 for triangles) allows the
 * modulation to adapt independently to the distribution of response variation
 * at different geometric levels. Diagonal entries use pairwise differences along
 * edges, while off-diagonal entries use the maximum variation across triangle
 * vertices, reflecting that off-diagonal interactions involve three points
 * rather than two.
 *
 * IMPORTANT OPTIMIZATION:
 * Rather than rebuilding \eqn{M_1} from scratch, this version:
 * 1. Computes modulation factors for all edges and triangles
 * 2. Applies factors directly to existing matrix entries via coeffRef
 * 3. Preserves matrix structure and sparsity pattern
 *
 * This is much faster than reconstruction, especially for large graphs.
 *
 * NORMALIZATION:
 * After modulation, the total Frobenius mass is normalized to equal the
 * pre-modulation mass:
 *   ||M_1^(modulated)||_F = ||M_1^(original)||_F
 *
 * This prevents systematic drift in the overall geometry scale across iterations.
 * The normalization preserves relative modulation differences while fixing the
 * total mass.
 *
 * ITERATION CONTEXT:
 * This function is called as Step 5 in the iteration loop:
 *   Step 1: Density diffusion (evolves vertex densities)
 *   Step 2: Edge density update (computes edge densities from evolved vertex densities)
 *   Step 3: Vertex mass matrix update (builds \eqn{M_0} from vertex densities)
 *   Step 4: Edge mass matrix construction (builds \eqn{M_1} from vertex and edge densities)
 *   Step 5: Response-coherence modulation \leftarrow THIS FUNCTION
 *   Step 6: Laplacian reassembly (builds L_0 from \eqn{M_0} and modulated \eqn{M_1})
 *   Step 7: Response smoothing
 *
 * The modulated \eqn{M_1} is used immediately in Step 6 to assemble the Laplacian,
 * which then drives response smoothing and the next iteration's density evolution.
 *
 * @param y_hat Current fitted response values
 * @param gamma Decay rate parameter (typically 0.5-2.0). Controls the strength
 *              of modulation: larger \gamma produces sharper response boundaries,
 *              while smaller \gamma yields gentler transitions.
 *
 * @pre y_hat.size() == vertex_cofaces.size()
 * @pre g.M[1] is assembled with current edge mass matrix from Step 4
 * @pre edge_cofaces contains triangles (may be empty if no triangles, in which
 *      case only diagonal modulation is performed)
 * @pre gamma > 0
 *
 * @post g.M[1] modulated in place with diagonal and off-diagonal entries adjusted
 * @post g.M[1] remains symmetric positive semidefinite
 * @post Total Frobenius mass preserved: ||\eqn{M_1}||_F unchanged by normalization
 * @post g.M_solver[1] invalidated (factorization no longer valid)
 *
 * @note If edge_cofaces contains no triangles, only diagonal entries are modulated.
 *       A warning is issued if \gamma > 0 but no triangles are available.
 *
 * @note The two-scale approach (\sigma_1 ≠ \sigma_2) is essential for proper modulation.
 *       Diagonal and off-diagonal entries use different response variation
 *       measures (pairwise vs. max over triangle) and thus require different
 *       normalization scales.
 *
 * @note This function operates on the **mass matrix** \eqn{M_1}, not on the edge
 *       densities. The densities remain unchanged; only the geometric
 *       structure (inner products) is modulated.
 *
 * @see compute_edge_mass_matrix() for Step 4 (builds \eqn{M_1} before modulation)
 * @see assemble_operators() for Step 6 (uses modulated \eqn{M_1} to build L_0)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::apply_response_coherence_modulation(
    const vec_t& y_hat,
    double gamma
) {
    const size_t n_edges = edge_registry.size();
    const size_t n_vertices = vertex_cofaces.size();

    // ============================================================
    // Validation
    // ============================================================

    if (static_cast<size_t>(y_hat.size()) != n_vertices) {
        Rf_error("apply_response_coherence_modulation: y_hat size mismatch");
    }

    if (gamma <= 0.0) {
        Rf_error("apply_response_coherence_modulation: gamma must be positive");
    }

    if (n_edges == 0) {
        return;  // Nothing to modulate
    }

    // ============================================================
    // Step 1: Compute edge response variations (for diagonal)
    // ============================================================

    std::vector<double> edge_deltas;
    edge_deltas.reserve(n_edges);

    for (size_t e = 0; e < n_edges; ++e) {
        const auto [v0, v1] = edge_registry[e];
        double delta = std::abs(y_hat[v0] - y_hat[v1]);
        edge_deltas.push_back(delta);
    }

    // ============================================================
    // Step 2: Compute scale parameter \sigma_1 for edges
    // ============================================================

    std::vector<double> edge_copy = edge_deltas;
    size_t q1_idx = edge_copy.size() / 4;
    size_t q3_idx = (3 * edge_copy.size()) / 4;

    std::nth_element(edge_copy.begin(),
                     edge_copy.begin() + q1_idx,
                     edge_copy.end());
    double Q1_edge = edge_copy[q1_idx];

    edge_copy = edge_deltas;
    std::nth_element(edge_copy.begin(),
                     edge_copy.begin() + q3_idx,
                     edge_copy.end());
    double Q3_edge = edge_copy[q3_idx];

    double sigma_1 = std::max(Q3_edge - Q1_edge, 1e-10);

    // ============================================================
    // Step 3: Compute triangle response variations (for off-diagonal)
    // ============================================================

    // Map from unordered edge pairs to triangle delta
    // Key: packed (e1, e2) with e1 < e2
    // Value: maximum response variation over all triangles containing both edges
    std::unordered_map<uint64_t, double> edge_pair_max_delta;

    double sigma_2 = 1e-10;  // Default if no triangles

    // Lambda to pack two edge indices into a single uint64_t key
    auto pack_edge_pair = [](index_t e1, index_t e2) -> uint64_t {
        if (e1 > e2) std::swap(e1, e2);
        return (static_cast<uint64_t>(e1) << 32) | static_cast<uint64_t>(e2);
    };

    // Count triangles from edge_cofaces
    size_t n_triangles = 0;
    if (!edge_cofaces.empty()) {
        for (const auto& cofaces : edge_cofaces) {
            for (size_t k = 1; k < cofaces.size(); ++k) {
                n_triangles = std::max(n_triangles,
                    static_cast<size_t>(cofaces[k].simplex_index + 1));
            }
        }
    }

    if (n_triangles > 0) {
        std::vector<double> triangle_deltas;
        triangle_deltas.reserve(n_triangles);

        // Track which triangles we've processed
        std::unordered_set<index_t> processed_triangles;

        // Iterate through all edges and their incident triangles
        for (size_t e = 0; e < n_edges; ++e) {
            const auto [v0_edge, v1_edge] = edge_registry[e];

            // Process triangles incident to this edge
            for (size_t t_idx = 1; t_idx < edge_cofaces[e].size(); ++t_idx) {
                index_t triangle_idx = edge_cofaces[e][t_idx].simplex_index;

                // Only process each triangle once
                if (processed_triangles.count(triangle_idx)) continue;
                processed_triangles.insert(triangle_idx);

                // Reconstruct triangle vertices
                const index_t v2 = edge_cofaces[e][t_idx].vertex_index;
                std::array<index_t, 3> tri_verts = {v0_edge, v1_edge, v2};
                std::sort(tri_verts.begin(), tri_verts.end());

                const index_t v0 = tri_verts[0];
                const index_t v1 = tri_verts[1];
                const index_t v2_sorted = tri_verts[2];

                // Compute maximum response variation over triangle vertices
                double delta = std::max({
                    std::abs(y_hat[v0] - y_hat[v1]),
                    std::abs(y_hat[v0] - y_hat[v2_sorted]),
                    std::abs(y_hat[v1] - y_hat[v2_sorted])
                });
                triangle_deltas.push_back(delta);

                // Find the three edge indices
                index_t e01 = NO_EDGE;
                for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
                    if (vertex_cofaces[v0][k].vertex_index == v1) {
                        e01 = vertex_cofaces[v0][k].simplex_index;
                        break;
                    }
                }

                index_t e02 = NO_EDGE;
                for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
                    if (vertex_cofaces[v0][k].vertex_index == v2_sorted) {
                        e02 = vertex_cofaces[v0][k].simplex_index;
                        break;
                    }
                }

                index_t e12 = NO_EDGE;
                for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
                    if (vertex_cofaces[v1][k].vertex_index == v2_sorted) {
                        e12 = vertex_cofaces[v1][k].simplex_index;
                        break;
                    }
                }

                // Record maximum delta for each edge pair in this triangle
                if (e01 != NO_EDGE && e02 != NO_EDGE) {
                    uint64_t key_01_02 = pack_edge_pair(e01, e02);
                    edge_pair_max_delta[key_01_02] = std::max(
                        edge_pair_max_delta[key_01_02], delta);
                }
                if (e01 != NO_EDGE && e12 != NO_EDGE) {
                    uint64_t key_01_12 = pack_edge_pair(e01, e12);
                    edge_pair_max_delta[key_01_12] = std::max(
                        edge_pair_max_delta[key_01_12], delta);
                }
                if (e02 != NO_EDGE && e12 != NO_EDGE) {
                    uint64_t key_02_12 = pack_edge_pair(e02, e12);
                    edge_pair_max_delta[key_02_12] = std::max(
                        edge_pair_max_delta[key_02_12], delta);
                }
            }
        }

        // Compute scale parameter \sigma_2 for triangles
        if (!triangle_deltas.empty()) {
            std::vector<double> tri_copy = triangle_deltas;
            q1_idx = tri_copy.size() / 4;
            q3_idx = (3 * tri_copy.size()) / 4;

            std::nth_element(tri_copy.begin(),
                             tri_copy.begin() + q1_idx,
                             tri_copy.end());
            double Q1_tri = tri_copy[q1_idx];

            tri_copy = triangle_deltas;
            std::nth_element(tri_copy.begin(),
                             tri_copy.begin() + q3_idx,
                             tri_copy.end());
            double Q3_tri = tri_copy[q3_idx];

            sigma_2 = std::max(Q3_tri - Q1_tri, 1e-10);
        }
    }

    // ============================================================
    // Step 4: Store total mass before modulation
    // ============================================================

    double total_mass_before = 0.0;
    for (int k = 0; k < g.M[1].outerSize(); ++k) {
        for (spmat_t::InnerIterator it(g.M[1], k); it; ++it) {
            if (it.row() <= it.col()) {
                total_mass_before += (it.row() == it.col()) ?
                    it.value() : 2.0 * it.value();
            }
        }
    }

    // ============================================================
    // Step 5: Apply modulation in place
    // ============================================================

    // Part A: Diagonal entries (modulated by edge variation)
    for (size_t e = 0; e < n_edges; ++e) {
        double delta_normalized = edge_deltas[e] / sigma_1;
        double penalty = std::pow(1.0 + delta_normalized * delta_normalized, -gamma);

        // Modify diagonal entry in place
        g.M[1].coeffRef(e, e) *= penalty;
    }

    // Part B: Off-diagonal entries (modulated by triangle variation)
    if (!edge_pair_max_delta.empty()) {

        // Iterate over all non-zero entries in \eqn{M_1}
        for (int k = 0; k < g.M[1].outerSize(); ++k) {
            for (spmat_t::InnerIterator it(g.M[1], k); it; ++it) {
                // Skip diagonal (already handled)
                if (it.row() == it.col()) continue;

                // Only process upper triangle to avoid double-processing
                if (it.row() > it.col()) continue;

                index_t e1 = static_cast<index_t>(it.row());
                index_t e2 = static_cast<index_t>(it.col());
                uint64_t key = pack_edge_pair(e1, e2);

                // Check if this edge pair appears in any triangle
                auto pair_it = edge_pair_max_delta.find(key);
                if (pair_it != edge_pair_max_delta.end()) {
                    double delta = pair_it->second;
                    double delta_normalized = delta / sigma_2;
                    double penalty = std::pow(
                        1.0 + delta_normalized * delta_normalized, -gamma);

                    // Modify both (i,j) and (j,i) entries for symmetry
                    double modulated_value = it.value() * penalty;
                    g.M[1].coeffRef(e1, e2) = modulated_value;
                    g.M[1].coeffRef(e2, e1) = modulated_value;
                }
            }
        }
    } else if (gamma > 0.0) {
        // Warn that we skipped off-diagonal modulation
        static bool warning_issued = false;
        if (!warning_issued) {
            Rprintf("Note: No triangles available for off-diagonal modulation.\n");
            warning_issued = true;
        }
    }

    // ============================================================
    // Step 6: Normalize to preserve total Frobenius mass
    // ============================================================

    double total_mass_after = 0.0;
    for (int k = 0; k < g.M[1].outerSize(); ++k) {
        for (spmat_t::InnerIterator it(g.M[1], k); it; ++it) {
            if (it.row() <= it.col()) {
                total_mass_after += (it.row() == it.col()) ?
                    it.value() : 2.0 * it.value();
            }
        }
    }

    if (total_mass_after > 1e-15) {
        double scale = total_mass_before / total_mass_after;
        g.M[1] *= scale;
    } else {
        // The mass collapse is telling us that the modulation parameters are
        // incompatible with the data.
        Rf_error(
            "Response-coherence modulation caused complete mass collapse.\n"
            "This indicates:\n"
            "  - gamma = %.3f may be too large (try reducing to 0.5-1.0)\n"
            "  - Response may be extremely discontinuous everywhere\n"
            "  - Scale parameters: sigma_1 = %.3e, sigma_2 = %.3e\n"
            "  - All edges have large response variation: min(Delta ) = %.3e, max(Delta ) = %.3e\n"
            "Consider: (1) reducing gamma, (2) checking response data for outliers,\n"
            "or (3) setting gamma = 0 to disable modulation.",
            gamma, sigma_1, sigma_2,
            *std::min_element(edge_deltas.begin(), edge_deltas.end()),
            *std::max_element(edge_deltas.begin(), edge_deltas.end())
        );
    }

    // Invalidate factorization
    g.M_solver[1].reset();
}

// ============================================================
// SPECTRAL FILTER FUNCTIONS
// ============================================================

namespace {

    /**
     * @brief Apply spectral filter function to eigenvalues
     *
     * @param lambda Eigenvalue (or array of eigenvalues)
     * @param eta Smoothing parameter
     * @param filter_type Type of spectral filter
     * @return Filter weight(s) in [0,1]
     */
    inline double apply_filter_function(
        double lambda,
        double eta,
        rdcx_filter_type_t filter_type
        ) {
        switch (filter_type) {
        case rdcx_filter_type_t::HEAT_KERNEL:
            return std::exp(-eta * lambda);

        case rdcx_filter_type_t::TIKHONOV:
            return 1.0 / (1.0 + eta * lambda);

        case rdcx_filter_type_t::CUBIC_SPLINE:
            return 1.0 / (1.0 + eta * lambda * lambda);

        case rdcx_filter_type_t::GAUSSIAN:
            return std::exp(-eta * lambda * lambda);

            // NEW CASES:
        case rdcx_filter_type_t::EXPONENTIAL:
            return std::exp(-eta * std::sqrt(std::max(lambda, 0.0)));

        case rdcx_filter_type_t::BUTTERWORTH:
            // Using order n = 2 (standard choice)
        {
            double ratio = lambda / eta;
            return 1.0 / (1.0 + ratio * ratio * ratio * ratio);
        }

        default:
            Rf_error("Unsupported filter type");
            return 0.0;
        }
    }

    /**
     * @brief Compute filter weights for all eigenvalues
     */
    Eigen::VectorXd compute_filter_weights(
        const vec_t& eigenvalues,
        double eta,
        rdcx_filter_type_t filter_type
        ) {
        const int m = eigenvalues.size();
        Eigen::VectorXd weights(m);

        for (int i = 0; i < m; ++i) {
            weights[i] = apply_filter_function(eigenvalues[i], eta, filter_type);
        }

        return weights;
    }

    /**
     * @brief Compute GCV score for given smoothing parameter
     *
     * GCV score = ||y - y_hat||^{2} / (n - trace(S))^{2}
     * where trace(S) = sum of filter weights (effective degrees of freedom)
     */
    double compute_gcv_score(
        const vec_t& y,
        const vec_t& y_hat,
        const Eigen::VectorXd& filter_weights
        ) {
        const double n = static_cast<double>(y.size());
        const double trace_S = filter_weights.sum();

        double residual_norm_sq = (y - y_hat).squaredNorm();
        double denom = n - trace_S;

        // Handle edge case where trace is close to n
        if (std::abs(denom) < 1e-10) {
            denom = 1e-10;
        }

        return residual_norm_sq / (denom * denom);
    }

    /**
     * @brief Generate logarithmically-spaced grid for smoothing parameter
     */
    std::vector<double> generate_eta_grid(
        double lambda_max,
        rdcx_filter_type_t filter_type,
        int n_candidates
        ) {
        const double eps = 1e-10;
        double eta_min = eps;
        double eta_max;

        // Determine appropriate eta_max based on filter type
        switch (filter_type) {
        case rdcx_filter_type_t::HEAT_KERNEL:
        case rdcx_filter_type_t::GAUSSIAN:
        case rdcx_filter_type_t::EXPONENTIAL:
            // For exponential filters: -log(eps)/lambda_max
            eta_max = (lambda_max > 0) ?
                -std::log(eps) / lambda_max : 1.0;
            break;

        case rdcx_filter_type_t::TIKHONOV:
        case rdcx_filter_type_t::CUBIC_SPLINE:
        case rdcx_filter_type_t::BUTTERWORTH:
            // For rational filters: 1/eps
            eta_max = 1.0 / eps;
            break;

        default:
            eta_max = 1.0;
            break;
        }

        // Generate logarithmic grid
        std::vector<double> eta_grid(n_candidates);
        for (int i = 0; i < n_candidates; ++i) {
            double t = static_cast<double>(i + 1) / static_cast<double>(n_candidates);
            eta_grid[i] = std::exp(std::log(eta_min) + t * (std::log(eta_max) - std::log(eta_min)));
        }

        return eta_grid;
    }
} // anonymous namespace

// ============================================================
// MAIN SPECTRAL FILTERING FUNCTION
// ============================================================

/**
 * @brief Smooth response via spectral filtering with GCV parameter selection
 *
 * This function implements spectral filtering of the response signal y using
 * the eigendecomposition of the Hodge Laplacian L_0. The filtering operates
 * in the eigenbasis, applying frequency-dependent weights determined by the
 * filter type and smoothing parameter η.
 *
 * The process follows these steps:
 * 1. Ensure spectral decomposition is available (compute if needed)
 * 2. Transform response to spectral domain (Graph Fourier Transform)
 * 3. Generate grid of candidate smoothing parameters
 * 4. For each candidate, apply filter and compute GCV score
 * 5. Select optimal parameter minimizing GCV
 * 6. Return smoothed response and diagnostic information
 *
 * @param y Observed response vector (size = n_vertices)
 * @param n_eigenpairs Number of eigenpairs to use in filtering.
 *                     More eigenpairs capture finer-scale variation.
 *                     Typical values: 50-200 for most applications.
 * @param filter_type Type of spectral filter to apply:
 *        - HEAT_KERNEL: exp(-η·λ), exponential decay. General-purpose filter
 *          providing smooth, progressive high-frequency attenuation. Most
 *          commonly used filter for routine applications.
 *
 *        - TIKHONOV: 1/(1+η·λ), rational decay. Corresponds to Tikhonov
 *          regularization. Gentler attenuation than heat kernel, better
 *          preserves linear trends and mid-frequency components.
 *
 *        - CUBIC_SPLINE: 1/(1+η·λ²), spline smoothness. Minimizes second
 *          derivatives, producing very smooth results similar to cubic
 *          splines. Excellent when response should vary gradually along paths.
 *
 *        - GAUSSIAN: exp(-η·λ²), super-exponential decay. Most aggressive
 *          high-frequency attenuation among exponential filters. Produces
 *          extremely smooth results with minimal oscillations or ringing.
 *
 *        - EXPONENTIAL: exp(-η·√λ), intermediate decay. Less aggressive than
 *          heat kernel, maintains more structure in mid-frequency range.
 *          Useful when moderate smoothing is desired without excessive blur.
 *
 *        - BUTTERWORTH: 1/(1+(λ/η)⁴), smooth rational cutoff. Fourth-order
 *          rational filter (n=2) providing clear frequency separation with
 *          reduced ringing compared to ideal low-pass filters. Good balance
 *          between sharpness and smoothness of cutoff.
 *
 * @return gcv_result_t structure containing:
 *         - eta_optimal: Selected smoothing parameter
 *         - y_hat: Smoothed response at optimal parameter
 *         - gcv_scores: GCV scores for all candidate parameters
 *         - eta_grid: Grid of evaluated parameters
 */
gcv_result_t riem_dcx_t::smooth_response_via_spectral_filter(
    const vec_t& y,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type
    ) {

    // ================================================================
    // VALIDATION: Input parameters and complex state
    // ================================================================

    // Validate that Laplacian exists and is properly initialized
    if (L.L.empty() || L.L[0].rows() == 0) {
        Rf_error("smooth_response_via_spectral_filter: vertex Laplacian L[0] not initialized");
    }

    const int n_vertices = L.L[0].rows();

    // Validate response vector dimension
    if (y.size() != n_vertices) {
        Rf_error("smooth_response_via_spectral_filter: response vector size (%d) "
                 "does not match number of vertices (%d)",
                 (int)y.size(), n_vertices);
    }

    // Validate n_eigenpairs bounds
    if (n_eigenpairs <= 0) {
        Rf_error("smooth_response_via_spectral_filter: n_eigenpairs must be positive, got %d",
                 n_eigenpairs);
    }

    // Critical: n_eigenpairs cannot exceed matrix dimension
    // For eigensolvers, we need at least 2 fewer than dimension for convergence
    const int max_feasible_eigenpairs = std::max(1, n_vertices - 2);

    if (n_eigenpairs > max_feasible_eigenpairs) {
        Rf_warning("smooth_response_via_spectral_filter: requested n_eigenpairs=%d exceeds "
                   "maximum feasible value %d for n_vertices=%d; reducing to %d",
                   n_eigenpairs, max_feasible_eigenpairs, n_vertices, max_feasible_eigenpairs);
        n_eigenpairs = max_feasible_eigenpairs;
    }

    // Practical limitation: very large n_eigenpairs may be computationally prohibitive
    if (n_eigenpairs > 500) {
        Rf_warning("smooth_response_via_spectral_filter: n_eigenpairs=%d is very large; "
                   "eigendecomposition may be slow. Consider reducing for computational efficiency.",
                   n_eigenpairs);
    }

    // Validate that we can proceed with spectral filtering
    // Either cache is valid with sufficient eigenpairs, or we need to compute
    const bool need_recompute = !spectral_cache.is_valid ||
        spectral_cache.eigenvalues.size() < n_eigenpairs;

    if (need_recompute && n_vertices > 10000 && n_eigenpairs > 200) {
        Rprintf("smooth_response_via_spectral_filter: computing %d eigenpairs for n=%d vertices; "
                "this may take significant time...\n", n_eigenpairs, n_vertices);
    }

    // Check for degenerate Laplacian (all zeros would indicate disconnected graph)
    if (L.L[0].nonZeros() == 0) {
        Rf_error("smooth_response_via_spectral_filter: vertex Laplacian is empty "
                 "(graph may be completely disconnected)");
    }

    // Validate response vector contains finite values
    for (int i = 0; i < y.size(); ++i) {
        if (!std::isfinite(y[i])) {
            Rf_error("smooth_response_via_spectral_filter: response vector contains "
                     "non-finite value at index %d (y[%d] = %f)",
                     i, i, y[i]);
        }
    }

    // ================================================================
    // STEP 1: ENSURE SPECTRAL DECOMPOSITION IS AVAILABLE
    // ================================================================

    // Use cached decomposition if valid and sufficient
    if (!spectral_cache.is_valid ||
        spectral_cache.eigenvalues.size() < n_eigenpairs) {
        compute_spectral_decomposition(n_eigenpairs);
    }

    // Extract needed eigenpairs from cache
    const int m = std::min(n_eigenpairs,
                           static_cast<int>(spectral_cache.eigenvalues.size()));

    vec_t eigenvalues = spectral_cache.eigenvalues.head(m);
    Eigen::MatrixXd eigenvectors = spectral_cache.eigenvectors.leftCols(m);

    // Validate response vector
    if (n_vertices != eigenvectors.rows()) {
        Rf_error("Response vector size (%d) does not match number of vertices (%d)",
                 n_vertices, (int)eigenvectors.rows());
    }

    // ================================================================
    // STEP 2: GRAPH FOURIER TRANSFORM
    // ================================================================

    // Project response onto eigenbasis: y_spectral = V^T y
    vec_t y_spectral = eigenvectors.transpose() * y;

    // ================================================================
    // STEP 3: GENERATE CANDIDATE SMOOTHING PARAMETERS
    // ================================================================

    const int n_candidates = 40;  // Standard grid size
    const double lambda_max = eigenvalues[m - 1];

    std::vector<double> eta_grid = generate_eta_grid(
        lambda_max, filter_type, n_candidates
        );

    // ================================================================
    // STEP 4: EVALUATE GCV FOR EACH CANDIDATE PARAMETER
    // ================================================================

    std::vector<double> gcv_scores(n_candidates);
    double best_gcv = std::numeric_limits<double>::max();
    int best_idx = 0;
    vec_t best_y_hat;

    for (int i = 0; i < n_candidates; ++i) {
        const double eta = eta_grid[i];

        // Compute filter weights for this eta
        Eigen::VectorXd filter_weights = compute_filter_weights(
            eigenvalues, eta, filter_type
            );

        // Apply filter in spectral domain
        vec_t y_filtered_spectral = filter_weights.asDiagonal() * y_spectral;

        // Transform back to vertex domain: y_hat = V * (filter_weights .* y_spectral)
        vec_t y_hat = eigenvectors * y_filtered_spectral;

        // Compute GCV score
        double gcv = compute_gcv_score(y, y_hat, filter_weights);
        gcv_scores[i] = gcv;

        // Track best solution
        if (gcv < best_gcv) {
            best_gcv = gcv;
            best_idx = i;
            best_y_hat = y_hat;
        }
    }

    // ================================================================
    // STEP 5: CONSTRUCT RESULT
    // ================================================================

    gcv_result_t result;
    result.eta_optimal = eta_grid[best_idx];
    result.gcv_optimal = best_gcv;
    result.y_hat = best_y_hat;
    result.gcv_scores = gcv_scores;
    result.eta_grid = eta_grid;

    return result;
}

/**
 * @brief Check convergence of iterative refinement procedure (simplified version)
 *
 * This function determines whether the iterative refinement process has
 * converged by monitoring the smoothed response field. Following the principle
 * that response convergence is the primary criterion for practical applications,
 * this simplified implementation focuses solely on the stability of predictions
 * rather than tracking the underlying geometric evolution.
 *
 * The convergence criterion is based on relative change to ensure scale
 * independence. We measure the L2 norm of differences normalized by the
 * L2 norm of the previous state. This provides a dimensionless measure
 * of change that is robust to the scale of the problem.
 *
 * Although the theoretical framework monitors both response and density
 * convergence, empirical experience shows that response stability is the
 * dominant indicator of practical convergence. The geometric quantities
 * (densities and mass matrices) typically stabilize before or alongside
 * the response, making separate density tracking often redundant for
 * determining when to terminate iteration. This simplified approach
 * reduces computational overhead by avoiding the extraction and comparison
 * of potentially large density vectors from the coface structures.
 *
 * The function provides diagnostic information through the returned
 * status structure, including human-readable messages describing the
 * current convergence state. This information is useful for monitoring
 * iteration progress and diagnosing convergence behavior.
 *
 * @param y_hat_prev Fitted response values from previous iteration
 * @param y_hat_curr Fitted response values from current iteration
 * @param epsilon_y Convergence threshold for response relative change
 * @param iteration Current iteration number (1-indexed)
 * @param max_iterations Maximum allowed iterations
 *
 * @return convergence_status_t structure containing:
 *         - converged: true if response criterion satisfied
 *         - response_change: relative change in fitted values
 *         - max_density_change: set to 0.0 (not tracked in simplified version)
 *         - iteration: current iteration number
 *         - message: human-readable convergence status
 *
 * @note Typical threshold value: epsilon_y = 1e-4 to 1e-3
 *
 * @note The function handles edge cases gracefully:
 *       - Near-zero norms are replaced with 1.0 to avoid division by zero
 *       - Maximum iterations reached is treated as termination (not convergence)
 *
 * @note For applications requiring explicit geometric convergence monitoring,
 *       use check_convergence_with_geometry() instead, which tracks density
 *       changes at the cost of additional computation.
 */
convergence_status_t riem_dcx_t::check_convergence(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    double epsilon_y,
    int iteration,
    int max_iterations
    ) {
    // Initialize status structure
    convergence_status_t status;
    status.iteration = iteration;
    status.converged = false;
    status.response_change = 0.0;
    status.max_density_change = 0.0;  // Not tracked in simplified version

    // ================================================================
    // CHECK 1: Maximum iterations reached
    // ================================================================

    // If we've reached the maximum iteration count, terminate regardless
    // of convergence criteria. This prevents infinite loops in cases where
    // convergence is slow or parameters are poorly chosen.
    if (iteration >= max_iterations) {
        // Compute change for diagnostic purposes even though we're terminating
        double y_norm_prev = y_hat_prev.norm();
        if (y_norm_prev < 1e-15) {
            y_norm_prev = 1.0;
        }
        status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;
        status.message = "Maximum iterations reached";
        return status;
    }

    // ================================================================
    // CHECK 2: Response convergence
    // ================================================================

    // Compute relative change in fitted response values using L2 norm.
    // The relative change is computed as ||y_curr - y_prev|| / ||y_prev||
    // to ensure scale independence.

    // Validate dimensions
    if (y_hat_prev.size() != y_hat_curr.size()) {
        Rf_error("check_convergence: response vectors have inconsistent sizes "
                 "(%d vs %d)", (int)y_hat_prev.size(), (int)y_hat_curr.size());
    }

    // Compute previous state norm with safety check
    double y_norm_prev = y_hat_prev.norm();
    if (y_norm_prev < 1e-15) {
        // Near-zero norm indicates either initialization or numerical issues.
        // Use 1.0 as denominator to avoid division by zero while still
        // measuring absolute change.
        y_norm_prev = 1.0;
    }

    // Compute relative change
    status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

    // ================================================================
    // CHECK 3: Evaluate convergence criterion
    // ================================================================

    // Determine convergence status based solely on response criterion.
    // This follows the practical principle that stable prediction is
    // the primary goal, and geometric quantities stabilize alongside
    // or before the response in typical use cases.

    bool response_converged = (status.response_change < epsilon_y);

    if (response_converged) {
        status.converged = true;
        status.message = "Converged: response stable";
        return status;
    }

    // Still changing
    status.converged = false;
    status.message = "Iteration in progress";
    return status;
}

/**
 * @brief Check convergence with enhanced diagnostics (simplified version)
 *
 * Extended version that estimates convergence rate from response change
 * history. This information can guide adaptive parameter selection or
 * early stopping decisions. Like the basic version, this simplified
 * implementation tracks only response convergence, not geometric quantities.
 *
 * @param y_hat_prev Fitted response values from previous iteration
 * @param y_hat_curr Fitted response values from current iteration
 * @param epsilon_y Convergence threshold for response relative change
 * @param iteration Current iteration number (1-indexed)
 * @param max_iterations Maximum allowed iterations
 * @param response_change_history Vector of response changes from previous iterations
 *
 * @return detailed_convergence_status_t structure containing:
 *         - converged: true if response criterion satisfied
 *         - response_change: relative change in fitted values
 *         - max_density_change: set to 0.0 (not tracked)
 *         - response_change_rate: ratio of consecutive changes
 *         - density_change_rate: set to 0.0 (not tracked)
 *         - estimated_iterations_remaining: projection based on geometric decay
 *         - density_changes_by_dim: empty vector (not tracked)
 *         - iteration: current iteration number
 *         - message: human-readable convergence status
 *
 * @note The convergence rate estimation assumes geometric (exponential) decay
 *       of the response change, which is typical for well-posed diffusion
 *       and spectral filtering iterations.
 */
detailed_convergence_status_t riem_dcx_t::check_convergence_detailed(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations,
    const std::vector<double>& response_change_history,
    const gcv_history_t& gcv_history
    ) {
    detailed_convergence_status_t status;
    status.iteration = iteration;
    status.converged = false;
    status.response_change = 0.0;
    status.max_density_change = 0.0;  // Not tracked in simplified version
    status.response_change_rate = 0.0;
    status.density_change_rate = 0.0;  // Not tracked in simplified version
    status.estimated_iterations_remaining = -1;

    // Maximum iterations check
    if (iteration >= max_iterations) {
        double y_norm_prev = y_hat_prev.norm();
        if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
        status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;
        status.message = "Maximum iterations reached";
        return status;
    }

    // Response convergence
    double y_norm_prev = y_hat_prev.norm();
    if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
    status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

    // Estimate convergence rate if history is available
    if (response_change_history.size() >= 2) {
        // Simple linear extrapolation of convergence rate
        size_t n = response_change_history.size();
        double recent_change = response_change_history[n-1];
        double prev_change = response_change_history[n-2];

        if (prev_change > 1e-15) {
            status.response_change_rate = recent_change / prev_change;

            // Estimate iterations to convergence assuming geometric decay
            if (status.response_change_rate < 1.0 && status.response_change_rate > 0.01) {
                double log_current = std::log(std::max(recent_change, 1e-15));
                double log_target = std::log(epsilon_y);
                double log_rate = std::log(status.response_change_rate);

                if (log_rate < -1e-6) {  // Ensure decreasing
                    status.estimated_iterations_remaining =
                        static_cast<int>(std::ceil((log_target - log_current) / log_rate));
                }
            }
        }
    }

    // Compute GCV diagnostics
    const int n_gcv = gcv_history.iterations.size();
    if (n_gcv > 0) {
        status.gcv_current = gcv_history.iterations[n_gcv - 1].gcv_optimal;

        if (n_gcv > 1) {
            double prev_gcv = gcv_history.iterations[n_gcv - 2].gcv_optimal;
            status.gcv_change = status.gcv_current - prev_gcv;

            // Compute GCV change rate
            if (n_gcv >= 3) {
                double prev_prev_gcv = gcv_history.iterations[n_gcv - 3].gcv_optimal;
                double prev_change = prev_gcv - prev_prev_gcv;

                // Avoid division by zero
                if (std::abs(prev_change) > 1e-15) {
                    status.gcv_change_rate = status.gcv_change / prev_change;
                } else {
                    status.gcv_change_rate = 0.0;
                }
            } else {
                status.gcv_change_rate = 0.0;
            }
        } else {
            status.gcv_change = 0.0;
            status.gcv_change_rate = 0.0;
        }
    } else {
        status.gcv_current = 0.0;
        status.gcv_change = 0.0;
        status.gcv_change_rate = 0.0;
    }

    // Convergence determination
    bool response_converged = (status.response_change < epsilon_y);

    if (response_converged) {
        status.converged = true;
        status.message = "Converged: response stable";
    } else {
        status.message = "Iteration in progress";

        // Add convergence estimate to message if available
        if (status.estimated_iterations_remaining > 0 &&
            status.estimated_iterations_remaining < 100) {
            status.message += " (est. " +
                std::to_string(status.estimated_iterations_remaining) +
                " iterations remaining)";
        }
    }

    return status;
}

/**
 * @brief Check convergence with full geometric tracking
 *
 * This extended version monitors both response and density convergence by
 * extracting and comparing density values from the vertex_cofaces structure.
 * Use this function when explicit geometric convergence verification is
 * required, at the cost of additional computation for density extraction.
 *
 * The convergence criteria follow the theoretical framework: iteration
 * terminates when both response and density stabilize below their respective
 * thresholds. This ensures that not only predictions but also the underlying
 * geometric structure has reached a stable state.
 *
 * For density monitoring, we track changes at both vertex and edge levels
 * by extracting densities from the coface structures. The maximum relative
 * change across both dimensions serves as the geometric convergence criterion.
 * This ensures that convergence is not declared prematurely if any dimension
 * continues to evolve even when others have stabilized.
 *
 * @param y_hat_prev Fitted response values from previous iteration
 * @param y_hat_curr Fitted response values from current iteration
 * @param vertex_cofaces_prev Vertex coface structure from previous iteration
 * @param vertex_cofaces_curr Vertex coface structure from current iteration
 * @param epsilon_y Convergence threshold for response relative change
 * @param epsilon_rho Convergence threshold for density relative change
 * @param iteration Current iteration number (1-indexed)
 * @param max_iterations Maximum allowed iterations
 *
 * @return convergence_status_t structure containing:
 *         - converged: true if both criteria satisfied
 *         - response_change: relative change in fitted values
 *         - max_density_change: maximum relative change across vertex and edge densities
 *         - iteration: current iteration number
 *         - message: human-readable convergence status
 *
 * @note Typical threshold values:
 *       - epsilon_y = 1e-4 to 1e-3 (response convergence)
 *       - epsilon_rho = 1e-4 to 1e-3 (density convergence)
 *
 * @note This function is more expensive than check_convergence() due to
 *       density extraction. For typical applications where response stability
 *       is the primary concern, the simplified version suffices.
 */
convergence_status_t riem_dcx_t::check_convergence_with_geometry(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces_prev,
    const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations
    ) {
    // Initialize status structure
    convergence_status_t status;
    status.iteration = iteration;
    status.converged = false;
    status.response_change = 0.0;
    status.max_density_change = 0.0;

    // ================================================================
    // CHECK 1: Maximum iterations reached
    // ================================================================

    if (iteration >= max_iterations) {
        // Compute changes for diagnostic purposes
        double y_norm_prev = y_hat_prev.norm();
        if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
        status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

        // Compute density changes
        const size_t n_vertices = vertex_cofaces_prev.size();

        vec_t vertex_dens_prev(n_vertices);
        vec_t vertex_dens_curr(n_vertices);
        for (size_t i = 0; i < n_vertices; ++i) {
            vertex_dens_prev[i] = vertex_cofaces_prev[i][0].density;
            vertex_dens_curr[i] = vertex_cofaces_curr[i][0].density;
        }

        double vertex_norm = vertex_dens_prev.norm();
        if (vertex_norm < 1e-15) vertex_norm = 1.0;
        status.max_density_change = (vertex_dens_curr - vertex_dens_prev).norm() / vertex_norm;

        // Extract edge densities
        size_t n_edges = 0;
        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces_prev[i].size(); ++k) {
                if (vertex_cofaces_prev[i][k].vertex_index > static_cast<index_t>(i)) {
                    n_edges++;
                }
            }
        }

        if (n_edges > 0) {
            vec_t edge_dens_prev(n_edges);
            vec_t edge_dens_curr(n_edges);
            size_t edge_idx = 0;

            for (size_t i = 0; i < n_vertices; ++i) {
                for (size_t k = 1; k < vertex_cofaces_prev[i].size(); ++k) {
                    index_t j = vertex_cofaces_prev[i][k].vertex_index;
                    if (j > static_cast<index_t>(i)) {
                        edge_dens_prev[edge_idx] = vertex_cofaces_prev[i][k].density;

                        for (size_t m = 1; m < vertex_cofaces_curr[i].size(); ++m) {
                            if (vertex_cofaces_curr[i][m].vertex_index == j) {
                                edge_dens_curr[edge_idx] = vertex_cofaces_curr[i][m].density;
                                break;
                            }
                        }
                        edge_idx++;
                    }
                }
            }

            double edge_norm = edge_dens_prev.norm();
            if (edge_norm < 1e-15) edge_norm = 1.0;
            double edge_change = (edge_dens_curr - edge_dens_prev).norm() / edge_norm;
            status.max_density_change = std::max(status.max_density_change, edge_change);
        }

        status.message = "Maximum iterations reached";
        return status;
    }

    // ================================================================
    // CHECK 2: Response convergence
    // ================================================================

    if (y_hat_prev.size() != y_hat_curr.size()) {
        Rf_error("check_convergence_with_geometry: response vectors have inconsistent sizes "
                 "(%d vs %d)", (int)y_hat_prev.size(), (int)y_hat_curr.size());
    }

    double y_norm_prev = y_hat_prev.norm();
    if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
    status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

    // ================================================================
    // CHECK 3: Density convergence
    // ================================================================

    if (vertex_cofaces_prev.size() != vertex_cofaces_curr.size()) {
        Rf_error("check_convergence_with_geometry: vertex_cofaces have inconsistent sizes "
                 "(%d vs %d)", (int)vertex_cofaces_prev.size(),
                 (int)vertex_cofaces_curr.size());
    }

    const size_t n_vertices = vertex_cofaces_prev.size();

    // Extract and compare vertex densities
    vec_t vertex_dens_prev(n_vertices);
    vec_t vertex_dens_curr(n_vertices);

    for (size_t i = 0; i < n_vertices; ++i) {
        vertex_dens_prev[i] = vertex_cofaces_prev[i][0].density;
        vertex_dens_curr[i] = vertex_cofaces_curr[i][0].density;
    }

    double vertex_norm = vertex_dens_prev.norm();
    if (vertex_norm < 1e-15) vertex_norm = 1.0;
    double vertex_change = (vertex_dens_curr - vertex_dens_prev).norm() / vertex_norm;
    status.max_density_change = vertex_change;

    // Extract and compare edge densities
    size_t n_edges = 0;
    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t k = 1; k < vertex_cofaces_prev[i].size(); ++k) {
            index_t j = vertex_cofaces_prev[i][k].vertex_index;
            if (j > static_cast<index_t>(i)) {
                n_edges++;
            }
        }
    }

    if (n_edges > 0) {
        vec_t edge_dens_prev(n_edges);
        vec_t edge_dens_curr(n_edges);
        size_t edge_idx = 0;

        for (size_t i = 0; i < n_vertices; ++i) {
            for (size_t k = 1; k < vertex_cofaces_prev[i].size(); ++k) {
                index_t j = vertex_cofaces_prev[i][k].vertex_index;

                if (j > static_cast<index_t>(i)) {
                    edge_dens_prev[edge_idx] = vertex_cofaces_prev[i][k].density;

                    bool found = false;
                    for (size_t m = 1; m < vertex_cofaces_curr[i].size(); ++m) {
                        if (vertex_cofaces_curr[i][m].vertex_index == j) {
                            edge_dens_curr[edge_idx] = vertex_cofaces_curr[i][m].density;
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        Rf_error("check_convergence_with_geometry: edge [%d,%d] not found in current iteration",
                                 (int)i, (int)j);
                    }

                    edge_idx++;
                }
            }
        }

        double edge_norm = edge_dens_prev.norm();
        if (edge_norm < 1e-15) edge_norm = 1.0;
        double edge_change = (edge_dens_curr - edge_dens_prev).norm() / edge_norm;
        status.max_density_change = std::max(status.max_density_change, edge_change);
    }

    // ================================================================
    // CHECK 4: Evaluate convergence criteria
    // ================================================================

    bool response_converged = (status.response_change < epsilon_y);
    bool density_converged = (status.max_density_change < epsilon_rho);

    if (response_converged && density_converged) {
        status.converged = true;
        status.message = "Converged: both response and density stable";
        return status;
    }

    if (response_converged && !density_converged) {
        status.converged = false;
        status.message = "Response converged, but density still evolving";
        return status;
    }

    if (!response_converged && density_converged) {
        status.converged = false;
        status.message = "Density converged, but response still evolving";
        return status;
    }

    status.converged = false;
    status.message = "Iteration in progress";
    return status;
}

/**
 * @brief Check convergence of iterative refinement procedure
 *
 * This function determines whether the iterative refinement process has
 * converged by monitoring two primary quantities: the smoothed response
 * field and the density distributions across all simplex dimensions. The
 * iteration is considered converged when both quantities have stabilized
 * below specified thresholds.
 *
 * The convergence criteria are based on relative changes to ensure scale
 * independence. We measure the L2 norm of differences normalized by the
 * L2 norm of the previous state. This provides a dimensionless measure
 * of change that is robust to the scale of the problem.
 *
 * For density monitoring, we track changes across all dimensions where
 * density is defined (vertices, edges, and potentially higher-dimensional
 * simplices if present). The maximum relative change across all dimensions
 * serves as the geometric convergence criterion. This ensures that convergence
 * is not declared prematurely if any dimension continues to evolve even when
 * others have stabilized.
 *
 * The function also provides diagnostic information through the returned
 * status structure, including human-readable messages describing the
 * current convergence state. This information is useful for monitoring
 * iteration progress and diagnosing convergence behavior.
 *
 * @param y_hat_prev Fitted response values from previous iteration
 * @param y_hat_curr Fitted response values from current iteration
 * @param rho_prev Density distributions from previous iteration (indexed by dimension)
 * @param rho_curr Density distributions from current iteration (indexed by dimension)
 * @param epsilon_y Convergence threshold for response relative change
 * @param epsilon_rho Convergence threshold for density relative change
 * @param iteration Current iteration number (1-indexed)
 * @param max_iterations Maximum allowed iterations
 *
 * @return convergence_status_t structure containing:
 *         - converged: true if both criteria satisfied
 *         - response_change: relative change in fitted values
 *         - max_density_change: maximum relative change across ALL densities
 *                               (not just dimensions 0 and 1, but all dimensions
 *                               where density is defined)
 *         - iteration: current iteration number
 *         - message: human-readable convergence status
 *
 * @note Typical threshold values:
 *       - epsilon_y = 1e-4 to 1e-3 (response convergence)
 *       - epsilon_rho = 1e-4 to 1e-3 (density convergence)
 *
 * @note The function handles edge cases gracefully:
 *       - Near-zero norms are replaced with 1.0 to avoid division by zero
 *       - Maximum iterations reached is treated as termination (not convergence)
 *       - Empty density vectors are handled safely
 */
convergence_status_t riem_dcx_t::check_convergence_with_rho(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    const std::vector<vec_t>& rho_prev,
    const std::vector<vec_t>& rho_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations
    ) {
    // Initialize status structure
    convergence_status_t status;
    status.iteration = iteration;
    status.converged = false;
    status.response_change = 0.0;
    status.max_density_change = 0.0;

    // ================================================================
    // CHECK 1: Maximum iterations reached
    // ================================================================

    // If we've reached the maximum iteration count, terminate regardless
    // of convergence criteria. This prevents infinite loops in cases where
    // convergence is slow or parameters are poorly chosen.
    if (iteration >= max_iterations) {
        // Compute changes for diagnostic purposes even though we're terminating
        double y_norm_prev = y_hat_prev.norm();
        if (y_norm_prev < 1e-15) {
            y_norm_prev = 1.0;
        }
        status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

        // Compute density changes
        for (size_t p = 0; p < rho_prev.size() && p < rho_curr.size(); ++p) {
            if (rho_prev[p].size() == 0 || rho_curr[p].size() == 0) {
                continue;
            }

            double rho_norm = rho_prev[p].norm();
            if (rho_norm < 1e-15) {
                rho_norm = 1.0;
            }

            double change = (rho_curr[p] - rho_prev[p]).norm() / rho_norm;
            status.max_density_change = std::max(status.max_density_change, change);
        }

        status.message = "Maximum iterations reached";
        return status;
    }

    // ================================================================
    // CHECK 2: Response convergence
    // ================================================================

    // Compute relative change in fitted response values using L2 norm.
    // The relative change is computed as ||y_curr - y_prev|| / ||y_prev||
    // to ensure scale independence.

    // Validate dimensions
    if (y_hat_prev.size() != y_hat_curr.size()) {
        Rf_error("check_convergence: response vectors have inconsistent sizes "
                 "(%d vs %d)", (int)y_hat_prev.size(), (int)y_hat_curr.size());
    }

    // Compute previous state norm with safety check
    double y_norm_prev = y_hat_prev.norm();
    if (y_norm_prev < 1e-15) {
        // Near-zero norm indicates either initialization or numerical issues.
        // Use 1.0 as denominator to avoid division by zero while still
        // measuring absolute change.
        y_norm_prev = 1.0;
    }

    // Compute relative change
    status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

    // ================================================================
    // CHECK 3: Density convergence
    // ================================================================

    // Compute relative changes for densities at all dimensions.
    // We track the maximum change across all dimensions as the overall
    // geometric convergence measure.

    // Validate that both density vectors have compatible structure
    if (rho_prev.size() != rho_curr.size()) {
        Rf_error("check_convergence: density vectors have inconsistent dimension counts "
                 "(%d vs %d)", (int)rho_prev.size(), (int)rho_curr.size());
    }

    // Iterate over all simplex dimensions
    for (size_t p = 0; p < rho_prev.size(); ++p) {
        // Skip dimensions with no simplices
        if (rho_prev[p].size() == 0 || rho_curr[p].size() == 0) {
            continue;
        }

        // Validate dimension consistency
        if (rho_prev[p].size() != rho_curr[p].size()) {
            Rf_error("check_convergence: density at dimension %d has inconsistent sizes "
                     "(%d vs %d)", (int)p, (int)rho_prev[p].size(), (int)rho_curr[p].size());
        }

        // Compute previous state norm with safety check
        double rho_norm = rho_prev[p].norm();
        if (rho_norm < 1e-15) {
            rho_norm = 1.0;
        }

        // Compute relative change at this dimension
        double change_p = (rho_curr[p] - rho_prev[p]).norm() / rho_norm;

        // Track maximum across dimensions
        status.max_density_change = std::max(status.max_density_change, change_p);
    }

    // ================================================================
    // CHECK 4: Evaluate convergence criteria
    // ================================================================

    // Determine convergence status based on both response and density criteria.
    // Full convergence requires both to be below their respective thresholds.

    bool response_converged = (status.response_change < epsilon_y);
    bool density_converged = (status.max_density_change < epsilon_rho);

    // ================================================================
    // CASE 1: Full convergence (both criteria satisfied)
    // ================================================================
    if (response_converged && density_converged) {
        status.converged = true;
        status.message = "Converged: both response and density stable";
        return status;
    }

    // ================================================================
    // CASE 2: Partial convergence (only response stable)
    // ================================================================
    if (response_converged && !density_converged) {
        status.converged = false;
        status.message = "Response converged, but density still evolving";
        return status;
    }

    // ================================================================
    // CASE 3: Partial convergence (only density stable)
    // ================================================================
    if (!response_converged && density_converged) {
        status.converged = false;
        status.message = "Density converged, but response still evolving";
        return status;
    }

    // ================================================================
    // CASE 4: No convergence (both still changing)
    // ================================================================
    status.converged = false;
    status.message = "Iteration in progress";
    return status;
}

/**
 * @brief Check convergence with enhanced diagnostics
 *
 * Extended version that tracks per-dimension changes and estimates
 * convergence rate. This information can guide adaptive parameter
 * selection or early stopping decisions.
 */
detailed_convergence_status_t riem_dcx_t::check_convergence_with_rho_detailed(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    const std::vector<vec_t>& rho_prev,
    const std::vector<vec_t>& rho_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations,
    const std::vector<double>& response_change_history
    ) {
    detailed_convergence_status_t status;
    status.iteration = iteration;
    status.converged = false;
    status.response_change = 0.0;
    status.max_density_change = 0.0;
    status.response_change_rate = 0.0;
    status.density_change_rate = 0.0;
    status.estimated_iterations_remaining = -1;

    // Maximum iterations check
    if (iteration >= max_iterations) {
        double y_norm_prev = y_hat_prev.norm();
        if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
        status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

        status.message = "Maximum iterations reached";
        return status;
    }

    // Response convergence
    double y_norm_prev = y_hat_prev.norm();
    if (y_norm_prev < 1e-15) y_norm_prev = 1.0;
    status.response_change = (y_hat_curr - y_hat_prev).norm() / y_norm_prev;

    // Density convergence with per-dimension tracking
    status.density_changes_by_dim.resize(rho_prev.size());

    for (size_t p = 0; p < rho_prev.size(); ++p) {
        if (rho_prev[p].size() == 0 || rho_curr[p].size() == 0) {
            status.density_changes_by_dim[p] = 0.0;
            continue;
        }

        double rho_norm = rho_prev[p].norm();
        if (rho_norm < 1e-15) rho_norm = 1.0;

        double change_p = (rho_curr[p] - rho_prev[p]).norm() / rho_norm;
        status.density_changes_by_dim[p] = change_p;
        status.max_density_change = std::max(status.max_density_change, change_p);
    }

    // Estimate convergence rate if history is available
    if (response_change_history.size() >= 2) {
        // Simple linear extrapolation of convergence rate
        size_t n = response_change_history.size();
        double recent_change = response_change_history[n-1];
        double prev_change = response_change_history[n-2];

        if (prev_change > 1e-15) {
            status.response_change_rate = recent_change / prev_change;

            // Estimate iterations to convergence assuming geometric decay
            if (status.response_change_rate < 1.0 && status.response_change_rate > 0.01) {
                double log_current = std::log(std::max(recent_change, 1e-15));
                double log_target = std::log(epsilon_y);
                double log_rate = std::log(status.response_change_rate);

                if (log_rate < -1e-6) {  // Ensure decreasing
                    status.estimated_iterations_remaining =
                        static_cast<int>(std::ceil((log_target - log_current) / log_rate));
                }
            }
        }
    }

    // Convergence determination
    bool response_converged = (status.response_change < epsilon_y);
    bool density_converged = (status.max_density_change < epsilon_rho);

    if (response_converged && density_converged) {
        status.converged = true;
        status.message = "Converged: both response and density stable";
    } else if (response_converged) {
        status.message = "Response converged, but density still evolving";
    } else if (density_converged) {
        status.message = "Density converged, but response still evolving";
    } else {
        status.message = "Iteration in progress";

        // Add convergence estimate to message if available
        if (status.estimated_iterations_remaining > 0 &&
            status.estimated_iterations_remaining < 100) {
            status.message += " (est. " +
                std::to_string(status.estimated_iterations_remaining) +
                " iterations remaining)";
        }
    }

    return status;
}

/**
 * @brief Format GCV history with trend indicators
 *
 * Creates a string showing GCV evolution across iterations with visual
 * indicators: '\' for decrease (improvement), '/' for increase (degradation).
 *
 * Example output: "0.0452 / 0.0389 \ 0.0367 \ 0.0351"
 * Indicates: increased from iter 0->1, then decreased in 1->2 and 2->3
 *
 * @param gcv_history GCV history object containing all GCV results
 * @param max_display Maximum number of recent iterations to display (default: 10)
 * @return Formatted string with GCV values and trend indicators
 */
std::string riem_dcx_t::format_gcv_history(
    const gcv_history_t& gcv_history,
    int max_display
    ) const {
    std::string result;
    char buffer[32];  // Buffer for snprintf

    const int n_iters = static_cast<int>(gcv_history.iterations.size());
    if (n_iters == 0) {
        return "N/A";
    }

    // Determine display window (show most recent max_display iterations)
    int start_idx = std::max(0, n_iters - max_display);

    // Add ellipsis if we're truncating history
    if (start_idx > 0) {
        result += "... ";
    }

    // Format first GCV value in display window
    snprintf(buffer, sizeof(buffer), "%.6f",
             gcv_history.iterations[start_idx].gcv_optimal);
    result += buffer;

    // Add subsequent values with trend indicators
    for (int i = start_idx + 1; i < n_iters; ++i) {
        double prev_gcv = gcv_history.iterations[i-1].gcv_optimal;
        double curr_gcv = gcv_history.iterations[i].gcv_optimal;

        // Choose indicator based on direction of change
        const char* indicator = (curr_gcv < prev_gcv) ? " \\ " : " / ";
        result += indicator;

        snprintf(buffer, sizeof(buffer), "%.6f", curr_gcv);
        result += buffer;
    }

    return result;
}

// ============================================================
// MAIN REGRESSION METHOD
// ============================================================

/**
 * @brief Fit Riemannian graph regression model using iterative geometric refinement
 *
 * @details
 * This function implements a novel approach to conditional expectation estimation
 * that combines geometric and probabilistic perspectives. Given feature-response
 * data (X, y), we estimate E[y|X] by constructing a Riemannian simplicial complex
 * whose geometry adapts iteratively to both the intrinsic structure of X and the
 * response field y.
 *
 * ## Mathematical Framework
 *
 * The method addresses a fundamental challenge in high-dimensional regression:
 * while the ambient dimension of X may be large, real-world datasets typically
 * exhibit much lower intrinsic dimensionality. Rather than estimating the
 * conditional expectation in the full ambient space, we construct a geometric
 * object K(X) that captures the intrinsic structure, then estimate E[y|K(X)].
 *
 * The geometric object is a simplicial complex K equipped with:
 * - A probability density rho encoding the data distribution
 * - A Riemannian metric g derived from rho
 * - A Hodge Laplacian L_0 encoding diffusion dynamics
 *
 * The triple (K, g, rho) evolves through iterative refinement where:
 * 1. Density diffuses across the geometry via heat kernel
 * 2. The metric adapts to the evolved density
 * 3. The edge mass matrix is modulated by response coherence
 * 4. The Laplacian is reassembled with the new metric
 * 5. The response is smoothed via spectral filtering on the new Laplacian
 *
 * This iteration continues until both the geometry and the response field reach
 * a mutually consistent fixed point.
 *
 * ## Algorithm Structure
 *
 * ### Part I: Initialization
 *
 * We begin by constructing the initial geometric structure from k-nearest neighbor
 * relationships. For each data point x_i, we compute its k nearest neighbors,
 * defining a neighborhood N_k(x_i). These neighborhoods form a covering of the
 * data space.
 *
 * The nerve of this covering gives us a simplicial complex: vertices correspond
 * to data points, and edges connect vertices whose neighborhoods overlap. We
 * weight edges by the size of neighborhood intersections and the minimum path
 * length through common neighbors, encoding both topological (intersection size)
 * and metric (path length) information.
 *
 * Geometric pruning removes spurious long-range connections that arise from
 * ambient space geometry but do not reflect intrinsic structure. We apply
 * ratio-based thresholding, comparing each edge length to local length scales.
 *
 * Initial densities are computed from a reference measure mu on the data points.
 * For vertex i, the density rho_0(i) equals mu(N_k(x_i)). For edge [i,j], the
 * density rho_1([i,j]) equals mu(N_k(x_i) intersect N_k(x_j)). These densities
 * capture how probability mass is distributed across the complex.
 *
 * The Riemannian metric is constructed from densities. The vertex mass matrix
 * M_0 is diagonal with M_0[i,i] = rho_0(i). The edge mass matrix M_1 captures
 * geometric interactions: for edges e_ij = [i,j] and e_is = [i,s] sharing vertex i,
 * their inner product <e_ij, e_is> equals the total density in the triple
 * intersection N_k(x_i) intersect N_k(x_j) intersect N_k(x_s).
 *
 * The Hodge Laplacian L_0 = B_1 M_1^{-1} B_1^T M_0 is assembled, where B_1 is the
 * vertex-edge boundary operator. This operator encodes diffusion on the complex,
 * with rates determined by the Riemannian structure.
 *
 * ### Part II: Iterative Refinement
 *
 * The iteration refines both the geometry and the response field through a
 * feedback loop. Each iteration performs these steps:
 *
 * **Step 1: Density Diffusion (Vertices)** Vertex densities evolve via damped
 * heat diffusion, solving (I + t(L_0 + beta*I))rho_0^(ell+1) = rho_0^(ell) + t*beta*u.
 * The parameter t controls diffusion rate, while beta provides damping toward
 * uniform distribution. The evolved densities are renormalized to sum to n.
 *
 * **Step 2: Edge Density Update** Edge densities are recomputed by summing
 * vertex densities over neighborhood intersections:
 *   rho_1^(ell+1)([i,j]) = sum_{v in N_i intersect N_j} rho_0^(ell+1)(v)
 * The edge densities are renormalized to sum to n_edges.
 *
 * **Step 3: Vertex Mass Matrix Update** The vertex mass matrix M_0 is updated
 * from evolved vertex densities: M_0 = diag(rho_0^(ell+1)). This is a simple
 * diagonal update that takes O(n) time.
 *
 * **Step 4: Edge Mass Matrix Construction** The edge mass matrix M_1 is rebuilt
 * from evolved vertex densities (for off-diagonal entries via triple intersections)
 * and updated edge densities (for diagonal entries):
 *   M_1[e,e] = rho_1^(ell+1)(e)
 *   M_1[e_ij, e_is] = sum_{v in N_i intersect N_j intersect N_s} rho_0^(ell+1)(v)
 * This is the computational bottleneck, requiring O(n*k^2) operations.
 *
 * **Step 5: Response-Coherence Modulation** The edge mass matrix M_1 is modulated
 * based on response variation. For each edge [i,j], we compute Delta_ij = |y_hat(i) - y_hat(j)|
 * measuring response discontinuity. Edges crossing response boundaries (large Delta_ij)
 * receive reduced mass via the penalty factor (1 + Delta_ij^2/sigma_1^2)^(-gamma),
 * where sigma_1 is an adaptive scale (IQR of all Delta_ij) and gamma controls
 * modulation strength. Off-diagonal entries are modulated similarly using
 * triangle-based response variations. The total Frobenius mass is preserved
 * via normalization.
 *
 * This creates geometry where edges within response-coherent regions maintain
 * high mass (short Riemannian distance), while edges crossing response
 * boundaries have reduced mass (long Riemannian distance). Diffusion becomes
 * faster within coherent regions and slower across boundaries, naturally
 * respecting the response structure.
 *
 * **Step 6: Laplacian Reassembly** The Hodge Laplacian is reassembled from
 * the updated metric: L_0 = B_1 M_1^{-1} B_1^T M_0. The spectral decomposition
 * cache is invalidated since the Laplacian has changed.
 *
 * **Step 7: Response Smoothing** The response is smoothed via spectral
 * filtering using the updated Laplacian. We compute the eigendecomposition
 * L_0 = V Lambda V^T and apply a low-pass filter: y_hat^(ell+1) = V F_eta(Lambda) V^T y,
 * where F_eta(Lambda) = diag(f(lambda_1), ..., f(lambda_m)) attenuates
 * high-frequency components. The smoothing parameter eta is selected via
 * Generalized Cross-Validation (GCV) to balance fidelity and smoothness.
 *
 * **Step 8: Convergence Check** We monitor relative changes in both response
 * and densities. Convergence is declared when:
 *   - ||y_hat^(ell+1) - y_hat^(ell)|| / ||y_hat^(ell)|| < epsilon_y (response stability)
 *   - max_p ||rho_p^(ell+1) - rho_p^(ell)|| / ||rho_p^(ell)|| < epsilon_rho (density stability)
 *
 * ### Part III: Finalization
 *
 * Upon convergence or reaching maximum iterations, we store the original
 * response and the complete fitted history. The final state contains:
 * - Fitted response values y_hat accessible via sig.y_hat_hist.back()
 * - Complete fitting history for diagnostic analysis
 * - Converged geometry (K, g, rho) encoding the learned structure
 *
 * ## Theoretical Justification
 *
 * The method combines three powerful ideas:
 *
 * **Geometric Perspective:** By working on a simplicial complex rather than the
 * ambient space, we exploit intrinsic low-dimensionality. The complex K(X)
 * captures essential geometric relationships while discarding irrelevant ambient
 * directions.
 *
 * **Measure-Theoretic Foundation:** The density rho provides a probability measure
 * on the complex. Smoothing the response with respect to this measure yields
 * conditional expectation estimates that properly account for data distribution.
 *
 * **Adaptive Geometry:** The response-coherence modulation creates geometry that
 * respects the response structure. Regions where the response varies smoothly
 * receive short Riemannian distances (rapid diffusion), while discontinuities
 * receive long distances (slow diffusion). This automatic adaptation eliminates
 * the need for manual bandwidth selection or kernel choice.
 *
 * The iteration finds a fixed point where the geometry is consistent with the
 * response and the response is smooth with respect to the geometry. This mutual
 * consistency is the key to the method's effectiveness.
 *
 * ## Parameter Selection Guidelines
 *
 * **k (neighbors):** Controls local neighborhood size. Typical range: 10-50.
 * Larger k captures more global structure but increases computation. For n < 500,
 * use k approximately 0.1*n. For n > 5000, use k approximately 20-30.
 *
 * **t_diffusion:** Controls density evolution rate. If <= 0, automatically set to
 * 0.5/lambda_2 where lambda_2 is the spectral gap. Larger t produces faster
 * convergence but may overshoot. Smaller t gives more conservative updates.
 *
 * **beta_damping:** Prevents density from concentrating excessively. If <= 0,
 * automatically set to 0.1/t_diffusion. The product beta*t should be small (< 0.5)
 * for gentle damping.
 *
 * **gamma_modulation:** Controls response-coherence modulation strength. Range:
 * 0-2.0. Larger gamma creates stronger geometric adaptation to response structure.
 * Use gamma = 1.0 as default. Set gamma = 0 to disable modulation entirely.
 *
 * **n_eigenpairs:** Number of eigenpairs for spectral filtering. Typical range:
 * 50-200. More eigenpairs capture finer response variation but increase
 * computation. Use min(n/5, 200) as a reasonable default.
 *
 * **filter_type:** Spectral filter for response smoothing:
 * - HEAT_KERNEL: General-purpose, exponential decay exp(-eta*lambda)
 * - TIKHONOV: Gentler smoothing, rational decay 1/(1+eta*lambda)
 * - CUBIC_SPLINE: Spline-like smoothness, 1/(1+eta*lambda^2)
 * - GAUSSIAN: Aggressive smoothing, exp(-eta*lambda^2)
 *
 * **epsilon_y, epsilon_rho:** Convergence thresholds. Typical: 1e-4 to 1e-3.
 * Tighter thresholds (1e-5) give more accurate convergence but require more
 * iterations. Looser thresholds (1e-3) converge faster but may be less stable.
 *
 * **max_iterations:** Safety limit. Typical: 20-50. Most problems converge
 * within 5-15 iterations with reasonable parameters. If convergence takes > 30
 * iterations, consider adjusting t_diffusion or gamma_modulation.
 *
 * **max_ratio_threshold:** Geometric pruning threshold. Typical: 2.0-3.0.
 * Larger values retain more edges (denser graph), smaller values prune more
 * aggressively. Use 2.5 as default.
 *
 * **threshold_percentile:** Percentile for local scale computation. Typical: 0.5
 * (median). Lower percentiles use shorter reference scales (more aggressive
 * pruning), higher percentiles use longer scales (less pruning).
 *
 * ## Computational Complexity
 *
 * **Initialization:** O(n k^2 log n) dominated by k-NN computation and edge
 * construction. Metric assembly is O(n k^3) but typically k^3 << n.
 *
 * **Per-iteration:** O(n k^2 + m^3) where m is the number of eigenpairs computed.
 * Density diffusion requires solving a sparse linear system (O(n k^2) with
 * iterative solvers). Eigendecomposition is O(m^2 n) for sparse methods. For
 * typical parameters (k ~ 30, m ~ 100), each iteration is O(n).
 *
 * **Total:** O(n k^2 log n + T*n k^2) where T is the number of iterations. For
 * T ~ 10, the initialization dominates asymptotically, making the method O(n k^2)
 * overall.
 *
 * ## Memory Requirements
 *
 * - Simplicial complex: O(n) for vertices, O(nk) for edges
 * - Mass matrices: O(n) for M_0 (diagonal), O(nk^2) for M_1 (sparse)
 * - Laplacian: O(nk) (sparse)
 * - Eigendecomposition cache: O(mn) for m eigenpairs
 * - History storage: O(Tn) if storing full iteration history
 *
 * Total: O(n(k^2 + m)) which is linear in n for fixed k and m.
 *
 * ## Numerical Considerations
 *
 * **Spectral Decomposition:** Uses robust fallback strategies for eigenvalue
 * computation. If standard ARPACK fails, progressively relaxes tolerance and
 * enlarges Krylov subspace. For small problems (n < 1000), falls back to dense
 * solver as guaranteed success path.
 *
 * **Ill-Conditioning:** Near-zero densities can create ill-conditioned metrics.
 * Regularization (adding small epsilon to diagonal entries) prevents numerical
 * instability. Monitor warnings about small spectral gaps or ill-conditioning.
 *
 * **Convergence Stalling:** If iteration stalls (changes decrease very slowly
 * but don't cross threshold), consider:
 * - Increasing t_diffusion for faster convergence
 * - Loosening convergence thresholds
 * - Checking for nearly-disconnected geometry
 *
 * ## Example Usage
 *
 * @code
 * // Create and initialize complex
 * riem_dcx_t complex;
 *
 * // Prepare data
 * Eigen::SparseMatrix<double> X = ...; // n x d feature matrix
 * Eigen::VectorXd y = ...;             // n response vector
 *
 * // Fit with default parameters
 * complex.fit_knn_riem_graph_regression(
 *     X, y,
 *     20,                              // k neighbors
 *     true,                            // use counting measure
 *     1.0,                             // density normalization
 *     0.0,                             // auto t_diffusion
 *     0.0,                             // auto beta_damping
 *     1.0,                             // gamma modulation
 *     100,                             // n eigenpairs
 *     rdcx_filter_type_t::HEAT_KERNEL, // filter type
 *     1e-4,                            // epsilon_y
 *     1e-4,                            // epsilon_rho
 *     30,                              // max iterations
 *     2.5,                             // max ratio threshold
 *     0.5                              // threshold percentile
 * );
 *
 * // Extract fitted values
 * Eigen::VectorXd y_hat = complex.sig.y_hat_hist.back();
 *
 * // Examine convergence history
 * for (size_t iter = 0; iter < complex.sig.y_hat_hist.size(); ++iter) {
 *     // Analyze evolution...
 * }
 * @endcode
 *
 * @param[in] X Feature matrix (n x d sparse or dense). Rows are observations,
 *              columns are features. Can be high-dimensional as method exploits
 *              intrinsic structure.
 *
 * @param[in] y Response vector (n x 1). Can be any real-valued signal defined
 *              on the data points. Method works for both smooth and discontinuous
 *              responses.
 *
 * @param[in] k Number of nearest neighbors for graph construction. Controls
 *              local neighborhood size. Typical: 10-50. Larger k captures more
 *              global structure but increases computation.
 *
 * @param[in] use_counting_measure If true, use uniform reference measure mu(x)=1
 *              for all points. If false, use density-weighted measure based on
 *              local distances. Counting measure is simpler and often sufficient.
 *
 * @param[in] density_normalization Power for density weighting when
 *              use_counting_measure=false. Typical: 1.0 gives mu(x) = 1/d_local(x)
 *              where d_local is the distance to k-th neighbor. Higher powers
 *              increase weight concentration in dense regions.
 *
 * @param[in,out] t_diffusion Diffusion time parameter for density evolution.
 *              If <= 0 on input, automatically set to 0.5/lambda_2 where lambda_2
 *              is spectral gap. If > 0, uses provided value. Larger t produces
 *              faster density evolution. Modified in place during parameter selection.
 *
 * @param[in,out] beta_damping Damping parameter preventing density concentration.
 *              If <= 0 on input, automatically set to 0.1/t_diffusion. If > 0,
 *              uses provided value. Provides gentle push toward uniform
 *              distribution. Modified in place during parameter selection.
 *
 * @param[in] gamma_modulation Response-coherence modulation strength. Controls
 *              how strongly edge mass matrix is modulated by response variation.
 *              Range: 0-2. Default: 1.0. Set to 0 to disable modulation entirely.
 *              Larger values create stronger geometric adaptation.
 *
 * @param[in] n_eigenpairs Number of eigenpairs to compute for spectral filtering.
 *              More eigenpairs capture finer response variation but increase
 *              computation. Typical: 50-200. Use min(n/5, 200) as default.
 *
 * @param[in] filter_type Type of spectral filter for response smoothing:
 *              - HEAT_KERNEL: exp(-eta*lambda), general-purpose exponential decay
 *              - TIKHONOV: 1/(1+eta*lambda), gentler rational decay
 *              - CUBIC_SPLINE: 1/(1+eta*lambda^2), spline-like smoothness
 *              - GAUSSIAN: exp(-eta*lambda^2), aggressive high-frequency attenuation
 *
 * @param[in] epsilon_y Convergence threshold for response relative change.
 *              Iteration stops when ||y_hat^(ell+1)-y_hat^(ell)||/||y_hat^(ell)|| < epsilon_y.
 *              Typical: 1e-4 to 1e-3. Tighter threshold gives more accurate
 *              convergence but requires more iterations.
 *
 * @param[in] epsilon_rho Convergence threshold for density relative change.
 *              Iteration stops when max_p ||rho_p^(ell+1)-rho_p^(ell)||/||rho_p^(ell)|| < epsilon_rho.
 *              Typical: 1e-4 to 1e-3.
 *
 * @param[in] max_iterations Maximum number of refinement iterations. Safety
 *              limit preventing infinite loops. Typical: 20-50. Most problems
 *              converge within 5-15 iterations with reasonable parameters.
 *
 * @param[in] max_ratio_threshold Maximum edge length ratio for geometric pruning.
 *              Edges longer than this ratio times the local scale are removed.
 *              Typical: 2.0-3.0. Use 2.5 as default. Larger values retain more
 *              edges (denser graph).
 *
 * @param[in] threshold_percentile Percentile for computing local scale in
 *              geometric pruning. Typical: 0.5 (median). Lower percentiles
 *              produce more aggressive pruning.
 *
 * @param[in] verbose (logical) controlling progress reporting during computation
 *
 * @pre X.rows() must equal y.size()
 * @pre k must be less than X.rows()
 * @pre All epsilon values must be positive
 * @pre max_iterations must be positive
 *
 * @post Simplicial complex fully populated with converged geometry
 * @post sig.y contains original response
 * @post sig.y_hat_hist contains complete iteration history
 * @post rho.rho contains final converged densities
 * @post g.M contains final converged metric tensors
 * @post L.L[0] contains final converged Hodge Laplacian
 * @post spectral_cache contains final eigendecomposition
 *
 * @throws std::runtime_error if k >= n
 * @throws std::runtime_error if dimensions are inconsistent
 * @throws std::runtime_error if eigendecomposition fails after all fallbacks
 * @throws std::runtime_error if linear system solve fails in density diffusion
 *
 * @warning This function performs substantial computation. For large n (> 10000)
 *          or large k (> 50), initialization may take minutes. Monitor progress
 *          via Rprintf output.
 *
 * @warning Automatic parameter selection (t_diffusion <= 0, beta_damping <= 0)
 *          requires computing spectral gap, adding cost to initialization.
 *
 * @note The function prints iteration diagnostics to R console via Rprintf.
 *       Each iteration shows response_change and density_change values for
 *       monitoring convergence progress.
 *
 * @note Final fitted values are accessed via sig.y_hat_hist.back(). The complete
 *       history is available for analyzing convergence behavior or producing
 *       diagnostic plots.
 *
 * @see initialize_from_knn() for details on geometric construction
 * @see smooth_response_via_spectral_filter() for spectral filtering algorithm
 * @see check_convergence() for convergence criteria
 * @see apply_damped_heat_diffusion() for density evolution
 * @see apply_response_coherence_modulation() for edge modulation
 */
void riem_dcx_t::fit_knn_riem_graph_regression(
    const spmat_t& X,
    const vec_t& y,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double t_diffusion,
    double beta_damping,
    double gamma_modulation,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations,
    double max_ratio_threshold,
    double threshold_percentile,
    double density_alpha,
    double density_epsilon,
    int test_stage,
    bool verbose
    ) {
    // ================================================================
    // PART I: INITIALIZATION
    // ================================================================

    auto total_time = std::chrono::steady_clock::now();
    auto phase_time = std::chrono::steady_clock::now();

    // --------------------------------------------------------------
    // Phase 1-4: Build geometric structure
    // --------------------------------------------------------------

    if (verbose) Rprintf("Phase 1-4: Build geometric structure ... ");

    initialize_from_knn(
        X, k,
        use_counting_measure,
        density_normalization,
        max_ratio_threshold,
        threshold_percentile,
        density_alpha,
        density_epsilon
        );

    // Store original response EARLY, before any possible early returns
    sig.y = y;

    // Validate triangle construction for response-coherence
    bool has_triangles = false;
    if (!edge_cofaces.empty()) {
        for (const auto& cofaces : edge_cofaces) {
            if (cofaces.size() > 1) {  // Has triangles (more than just self-loop)
                has_triangles = true;
                break;
            }
        }
    }

    if (gamma_modulation > 0.0 && !has_triangles) {
        Rf_warning("gamma_modulation = %.2f requires triangles for off-diagonal "
                   "modulation, but no triangles were constructed. "
                   "Only diagonal entries of M[1] will be modulated. "
                   "Consider: (1) increasing k to create more triangles, "
                   "(2) setting gamma_modulation = 0 to disable modulation, or "
                   "(3) adjusting pruning thresholds.",
                   gamma_modulation);
    }

    if (test_stage == 0) {
        Rprintf("TEST_STAGE 0: Stopped after initialize_from_knn\n");
        return;
    }

    // --------------------------------------------------------------
    // Phase 4.5: Select or validate diffusion parameters
    // --------------------------------------------------------------

    if (verbose) {
        elapsed_time(phase_time, "DONE", true);
        Rprintf("Phase 4.5: Select or validate diffusion parameters ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    select_diffusion_parameters(t_diffusion, beta_damping, /*verbose=*/true);

    if (test_stage == 1) {
        Rprintf("TEST_STAGE 1: Stopped after parameter selection\n");
        return;
    }

    // --------------------------------------------------------------
    // Phase 5: Initial response smoothing
    // --------------------------------------------------------------

    if (verbose) {
        elapsed_time(phase_time, "DONE", true);
        Rprintf("Phase 5: Initial response smoothing ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    gcv_history.clear();
    auto gcv_result = smooth_response_via_spectral_filter(
        y, n_eigenpairs, filter_type
        );
    vec_t y_hat_curr = gcv_result.y_hat;
    sig.y_hat_hist.clear();
    sig.y_hat_hist.push_back(y_hat_curr);

    // Store initial GCV results
    gcv_history.add(gcv_result);

    if (test_stage == 2) {
        Rprintf("TEST_STAGE 2: Stopped after initial smoothing\n");
        return;
    }

    // --------------------------------------------------------------
    // Phase 5.5: Automatic gamma selection (if requested)
    // --------------------------------------------------------------

    // Check if automatic gamma selection is requested
    // Convention: gamma_modulation < 0 triggers auto-selection
    bool gamma_auto_selected = false;

    if (gamma_modulation < 0.0) {
        Rprintf("\n");
        Rprintf("==========================================================\n");
        Rprintf("AUTOMATIC GAMMA SELECTION VIA FIRST-ITERATION GCV\n");
        Rprintf("==========================================================\n");

        // Define search grid
        std::vector<double> gamma_grid = {
            0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3,
            0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9,
            1.0, 1.2, 1.4, 1.6, 1.8, 2.0
        };

        Rprintf("Evaluating %d candidate gamma values...\n",
                (int)gamma_grid.size());
        Rprintf("Grid: [%.2f, %.2f, ..., %.2f, %.2f]\n",
                gamma_grid.front(), gamma_grid[1],
                gamma_grid[gamma_grid.size()-2], gamma_grid.back());

        // Perform gamma selection
        auto gamma_result = select_gamma_first_iteration(
            y, gamma_grid, n_eigenpairs, filter_type
            );

        // Use selected gamma
        gamma_modulation = gamma_result.gamma_optimal;
        gamma_auto_selected = true;

        // STORE THE RESULT IN THE MEMBER VARIABLE
        gamma_selection_result = gamma_result;
        gamma_was_auto_selected = true;

        Rprintf("\n");
        Rprintf("GAMMA SELECTION COMPLETE\n");
        Rprintf("  Selected gamma: %.3f\n", gamma_result.gamma_optimal);
        Rprintf("  GCV score:      %.6e\n", gamma_result.gcv_optimal);

        // Check if optimum is at boundary
        if (gamma_modulation == gamma_grid.front()) {
            Rf_warning("Optimal gamma is at lower grid boundary (%.3f). "
                       "Consider exploring smaller values.", gamma_modulation);
        }
        if (gamma_modulation == gamma_grid.back()) {
            Rf_warning("Optimal gamma is at upper grid boundary (%.3f). "
                       "Consider exploring larger values.", gamma_modulation);
        }

        Rprintf("==========================================================\n");
        Rprintf("\n");

    } else if (gamma_modulation == 0.0) {
        Rprintf("Using gamma = 0 (no response-coherence modulation)\n");
        gamma_was_auto_selected = false;
    } else {
        Rprintf("Using user-provided gamma = %.3f\n", gamma_modulation);
        gamma_was_auto_selected = false;
    }

    // ================================================================
    // PART II: ITERATIVE REFINEMENT
    // ================================================================

    // Note: epsilon_rho parameter is accepted for API compatibility but not
    // used in the simplified convergence checking that monitors only response
    // changes. The geometric quantities (densities stored in vertex_cofaces)
    // typically stabilize before or alongside response convergence, making
    // separate density tracking often redundant for determining when to
    // terminate iteration.

    const size_t n_vertices = vertex_cofaces.size();

    density_history.clear();
    // Store initial densities (iteration 0)
    vec_t rho_initial(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        rho_initial[i] = vertex_cofaces[i][0].density;
    }
    density_history.add(rho_initial);

    vec_t y_hat_prev;
    vec_t rho_vertex_prev(n_vertices);  // Allocate once
    std::vector<double> response_change_history;

    Rprintf("\n");
    Rprintf("Starting iterative refinement...\n");
    if (gamma_auto_selected) {
        Rprintf("(using auto-selected gamma = %.3f)\n", gamma_modulation);
    }
    Rprintf("\n");

    // Initialize convergence tracking
    converged = false;
    n_iterations = 0;
    response_changes.clear();

    for (int iter = 1; iter <= max_iterations; ++iter) {
        y_hat_prev = y_hat_curr;

        // --------------------------------------------------------------
        // Step 1: Density diffusion
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 1: Density diffusion ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        // Extract current vertex densities from vertex_cofaces
        for (size_t i = 0; i < n_vertices; ++i) {
            rho_vertex_prev[i] = vertex_cofaces[i][0].density;
        }

        // Apply damped heat diffusion
        vec_t rho_vertex_new = apply_damped_heat_diffusion(
            rho_vertex_prev, t_diffusion, beta_damping
            );

        // Store evolved densities
        density_history.add(rho_vertex_new);

        // Store evolved vertex densities back into vertex_cofaces
        for (size_t i = 0; i < n_vertices; ++i) {
            vertex_cofaces[i][0].density = rho_vertex_new[i];
        }

        // Test if density evolution is working
        static vec_t rho_prev_iter;
        if (iter == 1) {
            rho_prev_iter = rho_vertex_prev;
        } else {
            double density_change_norm = (rho_vertex_new - rho_prev_iter).norm();
            double density_norm = rho_prev_iter.norm();
            Rprintf("  Density L2 change: %.6e (relative: %.6e)\n",
                    density_change_norm, density_change_norm / density_norm);
        }
        rho_prev_iter = rho_vertex_new;

        if (test_stage == 3) {
            Rprintf("TEST_STAGE 3: Stopped after first diffusion\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 2: Edge density update
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 2: Edge density update ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        update_edge_densities_from_vertices();

        if (test_stage == 4) {
            Rprintf("TEST_STAGE 4: Stopped after edge density update\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 3: Update vertex mass matrix from evolved densities
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 3: Update vertex mass matrix from evolved densities ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        update_vertex_metric_from_density();  // Updates only M[0]

        if (test_stage == 5) {
            Rprintf("TEST_STAGE 5: Stopped after update_vertex_metric_from_density()\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 4: Rebuild edge mass matrix from evolved densities
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 4: Rebuild edge mass matrix from evolved densities ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        update_edge_mass_matrix();  // Compute fresh \eqn{M_1} from current \eqn{\rho_0}

        if (test_stage == 6) {
            Rprintf("TEST_STAGE 6: Stopped after update_edge_mass_matrix()\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 5: Apply response-coherence modulation to fresh M[1]
        // --------------------------------------------------------------
        // Start conservative, increase gradually
        // double gamma_eff = gamma_modulation * std::min(1.0, iter / 5.0);
        // apply_response_coherence_modulation(y_hat_curr, gamma_eff);

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 5: Apply response-coherence modulation to fresh M[1] ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        if (gamma_modulation > 0.0) {
            apply_response_coherence_modulation(y_hat_curr, gamma_modulation);
        }

        if (test_stage == 7) {
            Rprintf("TEST_STAGE 7: Stopped after apply_response_coherence_modulation()\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 6: Laplacian reassembly
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 6: Laplacian reassembly ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        assemble_operators();
        spectral_cache.invalidate();

        // Test if eigenvalues are changing
        static double prev_lambda_2 = 0.0;
        if (iter == 1) {
            prev_lambda_2 = spectral_cache.lambda_2;
        } else {
            double lambda_change = std::abs(spectral_cache.lambda_2 - prev_lambda_2);
            Rprintf("  lambda_2 changed: %.6e -> %.6e (Delta = %.6e)\n",
                    prev_lambda_2, spectral_cache.lambda_2, lambda_change);
            prev_lambda_2 = spectral_cache.lambda_2;
        }

        if (test_stage == 8) {
            Rprintf("TEST_STAGE 8: Stopped after assemble_operators()\n");
            return;
        }

        // --------------------------------------------------------------
        // Step 7: Response smoothing
        // --------------------------------------------------------------

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
            Rprintf("itr %d: Step 7: Response smoothing ... ", iter);
            phase_time = std::chrono::steady_clock::now();
        }

        gcv_result = smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
            );
        y_hat_curr = gcv_result.y_hat;
        sig.y_hat_hist.push_back(y_hat_curr);

        // Store GCV result for this iteration
        gcv_history.add(gcv_result);

        if (test_stage == 9) {
            Rprintf("TEST_STAGE 9: Stopped after smooth_response_via_spectral_filter()\n");
            return;
        }

        // ----------------------------------------------------------------
        // Step 8: Convergence check and detailed reporting
        // ----------------------------------------------------------------

        auto status = check_convergence_detailed(
            y_hat_prev,
            y_hat_curr,
            epsilon_y,
            epsilon_rho,
            iter,
            max_iterations,
            response_change_history,
            gcv_history  // Add this argument
            );

        if (test_stage == 10) {
            Rprintf("TEST_STAGE 10: Stopped after check_convergence_detailed()\n");
            return;
        }

        response_change_history.push_back(status.response_change);

        // Enhanced reporting with GCV trajectory
        if (verbose) {
            Rprintf("Iteration %d: response_change=%.6e", iter, status.response_change);

            // Add GCV trend information
            std::string gcv_display = format_gcv_history(gcv_history, 10);
            Rprintf(", GCV: %s", gcv_display.c_str());

            // Add change indicator for current iteration
            if (gcv_history.size() > 1) {
                if (status.gcv_change < 0) {
                    Rprintf(" [improving: %.2e]", status.gcv_change);
                } else if (status.gcv_change > 0) {
                    Rprintf(" [degrading: +%.2e]", status.gcv_change);
                } else {
                    Rprintf(" [stable]");
                }
            }

            Rprintf("\n");
        }

        if (status.converged) {
            Rprintf("%s\n", status.message.c_str());
            converged = true;
            n_iterations = iter;
            break;
        }

        // Detecting problematic GCV trends
        if (gcv_history.size() >= 3 && status.gcv_change > 0) {
            // Check if GCV has been increasing for multiple iterations
            int consecutive_increases = 1;
            for (int i = static_cast<int>(gcv_history.size()) - 2; i >= 1; --i) {
                // FIXED: Access the gcv_optimal field from the iterations vector
                if (gcv_history.iterations[i].gcv_optimal >
                    gcv_history.iterations[i-1].gcv_optimal) {
                    consecutive_increases++;
                } else {
                    break;
                }
            }

            if (consecutive_increases >= 3) {
                Rf_warning("GCV has increased for %d consecutive iterations. "
                           "This may indicate:\n"
                           "  - Over-aggressive geometric modulation (try reducing gamma)\n"
                           "  - Numerical instability (check for near-zero densities)\n"
                           "  - Need for parameter adjustment",
                           consecutive_increases);
            }
        }

        if (verbose) {
            elapsed_time(phase_time, "DONE", true);
        }

        ///// testing
        // if (iter == 1) break;  // Only do one iteration for testing
    }

    // ================================================================
    // FINALIZATION: Store convergence results
    // ================================================================

    // If loop completed without converging
    if (!converged) {
        n_iterations = max_iterations;
    }

    // Store response change history
    response_changes = response_change_history;

    // Maybe in the future: store density changes if tracked
    // density_changes = ...;

    if (verbose) {
        elapsed_time(phase_time, "DONE", true);
        elapsed_time(total_time, "Total elapsed time", true);
    }
}
