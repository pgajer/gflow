#include "riem_dcx.hpp"
#include "iknn_vertex.hpp" // for iknn_vertex_t
#include "kNN.h"
#include "set_wgraph.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>  // For Eigen::MatrixXd
#include <Eigen/Sparse> // For Eigen::SparseMatrix, Triplet
#include <Spectra/SymEigsSolver.h>          // For SymEigsSolver
#include <Spectra/MatOp/SparseSymMatProd.h> // For SparseSymMatProd

#include <R.h>

#include <cmath>
#include <limits>
#include <algorithm>

// ============================================================
// INITIALIZATION HELPERS (private member functions)
// ============================================================

/**
 * @brief Compute initial densities from reference measure
 *
 * Constructs the initial density functions ρ₀ on vertices and ρ₁ on edges by
 * aggregating the reference measure μ over appropriate neighborhoods. The vertex
 * density at point i measures the total mass in its k-neighborhood, while the
 * edge density on edge [i,j] measures the mass in the intersection of their
 * neighborhoods.
 *
 * The initialization follows the formula:
 *   ρ₀(i) = Σ_{j ∈ Ň_k(x_i)} μ({j})
 *   ρ₁([i,j]) = Σ_{v ∈ Ň_k(x_i) ∩ Ň_k(x_j)} μ({v})
 *
 * Both densities are normalized to sum to their respective simplex counts, so that
 * the average vertex has density 1 and the average edge has density 1. This
 * normalization ensures numerical stability and interpretability across different
 * dataset sizes.
 *
 * Densities are stored directly in the vertex_cofaces structure:
 *   - vertex_cofaces[i][0].density = ρ₀(i)
 *   - vertex_cofaces[i][k].density = ρ₁(edge [i, vertex_cofaces[i][k].vertex_index])
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

        // Aggregate measure over all neighbors j ∈ Ň_k(x_i)
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
    // of neighborhoods: Ň_k(x_i) ∩ Ň_k(x_j)

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
 * distribution ρ. The metric determines inner products between chains at
 * each dimension, encoding geometric information about lengths, angles,
 * areas, and volumes throughout the complex.
 *
 * For vertices (dimension 0), the metric is always diagonal with M₀ = diag(ρ₀).
 * This is a mathematical necessity: vertices have no geometric interaction
 * under the inner product construction.
 *
 * For edges (dimension 1), the full mass matrix M₁ captures essential geometric
 * information through triple intersections of neighborhoods. Two edges sharing
 * a common vertex have inner product determined by the density mass in the
 * triple intersection of their endpoints' neighborhoods. This encodes how the
 * edges are geometrically related through their shared vertex and overlapping
 * neighborhoods.
 *
 * The construction ensures positive semidefiniteness by design, as all inner
 * products arise from L²(μ) pairings of neighborhood indicator functions.
 *
 * @pre rho.rho[0] must contain vertex densities (normalized to sum to n)
 * @pre rho.rho[1] must contain edge densities (normalized to sum to n_edges)
 * @post g.M[0] contains diagonal vertex mass matrix
 * @post g.M[1] contains full edge mass matrix (sparse, symmetric, positive semidefinite)
 */
void riem_dcx_t::initialize_metric_from_density() {
    const size_t n_vertices = S[0].size();

    // ========================================================================
    // Part 1: Vertex mass matrix M₀ (diagonal)
    // ========================================================================

    // The vertex mass matrix is diagonal by mathematical necessity.
    // M₀ = diag(ρ₀(1), ρ₀(2), ..., ρ₀(n))

    g.M[0] = spmat_t(n_vertices, n_vertices);
    g.M[0].reserve(Eigen::VectorXi::Constant(n_vertices, 1));

    for (size_t i = 0; i < n_vertices; ++i) {
        // Apply regularization to ensure positive definiteness
        double mass = std::max(rho.rho[0][i], 1e-15);
        g.M[0].insert(i, i) = mass;
    }

    g.M[0].makeCompressed();

    // Clear any existing factorization (diagonal doesn't need factorization)
    g.M_solver[0].reset();

    // ========================================================================
    // Part 2: Edge mass matrix M₁ (full matrix via triple intersections)
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
 * Builds the complete edge mass matrix M₁ by computing inner products between
 * all pairs of edges through triple neighborhood intersections and setting
 * diagonal entries from edge densities. For edges e_ij = [i,j] and e_is = [i,s]
 * sharing vertex v_i, the inner product is:
 *
 *   ⟨e_ij, e_is⟩ = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j) ∩ N̂_k(x_s)} ρ₀(v)
 *
 * This measures the total vertex density in the triple intersection of
 * neighborhoods, encoding the geometric relationship between edges through
 * their shared vertex and overlapping neighborhoods.
 *
 * The diagonal entries are set directly from edge densities:
 *   M₁[e,e] = ρ₁(e)
 *
 * since the edge self-inner product equals the pairwise intersection mass:
 *   ⟨e_ij, e_ij⟩ = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j)} ρ₀(v) = ρ₁([i,j])
 *
 * The resulting matrix is symmetric positive semidefinite and sparse. For
 * kNN complexes with parameter k, each edge typically interacts with O(k²)
 * other edges, making sparse storage efficient.
 *
 * IMPLEMENTATION STRATEGY:
 * The function builds M₁ by:
 * 1. Adding diagonal entries from rho.rho[1] (edge densities)
 * 2. Iterating over all triangles in S[2]
 * 3. For each triangle [i,j,s], computing the triple intersection mass
 * 4. Adding off-diagonal entries for all three edge pairs in the triangle
 * 5. Assembling the symmetric sparse matrix from triplets
 *
 * This triangle-based approach ensures all edge pairs sharing a vertex are
 * accounted for, as every such pair appears in at least one triangle (if
 * the complex includes 2-simplices).
 *
 * COMPUTATIONAL COMPLEXITY: O(n * k²) where n is number of vertices and k is
 * the neighborhood size. This is the bottleneck operation in metric construction.
 *
 * The complexity arises from:
 * - Iterating over O(n·k²) triangles
 * - Computing each triple intersection in O(k) time
 * - Building sparse matrix from O(n·k²) triplets
 *
 * @pre S[0] and S[1] must be populated with vertices and edges
 * @pre S[2] must be populated with triangles (for off-diagonal entries)
 * @pre rho.rho[0] must contain current vertex densities
 * @pre rho.rho[1] must contain current edge densities
 * @pre neighbor_sets must be populated with kNN neighborhoods
 *
 * @post g.M[1] contains symmetric positive semidefinite edge mass matrix
 * @post g.M[1] has diagonal entries from rho.rho[1]
 * @post g.M[1] has off-diagonal entries from triple intersections using rho.rho[0]
 * @post g.M_solver[1] is reset (no factorization stored)
 *
 * @note If S[2] is empty (no triangles), only diagonal entries are created,
 *       resulting in a diagonal mass matrix. This loses the full Riemannian
 *       structure and reduces the method to a simpler graph-based approach.
 *
 * @note Regularization is applied to diagonal entries: M₁[e,e] = max(ρ₁(e), 1e-15)
 *       to ensure positive definiteness even if some edge densities are very small.
 *
 * @see initialize_metric_from_density() for initialization context
 * @see update_edge_mass_matrix() for iteration context
 * @see compute_edge_inner_product() for the triple intersection computation
 */
void riem_dcx_t::compute_edge_mass_matrix() {
    const size_t n_edges = S[1].size();
    std::vector<Eigen::Triplet<double>> triplets;

    // Add diagonal entries
    for (size_t e = 0; e < n_edges; ++e) {
        double diagonal = std::max(rho.rho[1][e], 1e-15);
        triplets.emplace_back(e, e, diagonal);
    }

    // Add off-diagonal entries from triangles
    if (S.size() > 2 && S[2].size() > 0) {
        for (size_t t = 0; t < S[2].size(); ++t) {
            const auto& tri_verts = S[2].simplex_verts[t];

            // Get the three edge IDs
            index_t e01 = S[1].get_id({tri_verts[0], tri_verts[1]});
            index_t e02 = S[1].get_id({tri_verts[0], tri_verts[2]});
            index_t e12 = S[1].get_id({tri_verts[1], tri_verts[2]});

            // Compute triple intersection mass once per triangle
            double mass = 0.0;
            for (index_t v : neighbor_sets[tri_verts[0]]) {
                if (neighbor_sets[tri_verts[1]].find(v) != neighbor_sets[tri_verts[1]].end() &&
                    neighbor_sets[tri_verts[2]].find(v) != neighbor_sets[tri_verts[2]].end()) {
                    mass += rho.rho[0][v];
                }
            }

            // Add all three off-diagonal pairs
            triplets.emplace_back(e01, e02, mass);
            triplets.emplace_back(e02, e01, mass);
            triplets.emplace_back(e01, e12, mass);
            triplets.emplace_back(e12, e01, mass);
            triplets.emplace_back(e02, e12, mass);
            triplets.emplace_back(e12, e02, mass);
        }
    }

    g.M[1] = spmat_t(n_edges, n_edges);
    g.M[1].setFromTriplets(triplets.begin(), triplets.end());
    g.M[1].makeCompressed();
    g.M_solver[1].reset();
}

/**
 * @brief Compute inner product between two edges sharing a vertex
 *
 * For edges e_ij = [i,j] and e_is = [i,s] sharing vertex v_i, computes:
 *   ⟨e_ij, e_is⟩ = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j) ∩ N̂_k(x_s)} ρ₀(v)
 *
 * This is the total vertex density in the triple intersection of neighborhoods,
 * representing the geometric relationship between the two edges through their
 * shared vertex and overlapping neighborhoods.
 *
 * COMPUTATIONAL STRATEGY:
 * To compute the triple intersection efficiently, we iterate over the smallest
 * set and check membership in the other two sets. For typical kNN graphs with
 * moderate k, this gives O(k) complexity per inner product computation.
 *
 * @param e1 Index of first edge
 * @param e2 Index of second edge
 * @param vertex_i Index of shared vertex
 * @return Inner product ⟨e1, e2⟩ at vertex_i
 */
double riem_dcx_t::compute_edge_inner_product(
    index_t e1,
    index_t e2,
    index_t vertex_i
) const {
    // Get the vertices of both edges
    const std::vector<index_t>& edge1_verts = S[1].simplex_verts[e1];
    const std::vector<index_t>& edge2_verts = S[1].simplex_verts[e2];

    // Find the other endpoints (not vertex_i)
    const index_t j = (edge1_verts[0] == vertex_i) ? edge1_verts[1] : edge1_verts[0];
    const index_t s = (edge2_verts[0] == vertex_i) ? edge2_verts[1] : edge2_verts[0];

    // Get the three neighborhood sets
    const std::unordered_set<index_t>& N_i = neighbor_sets[vertex_i];
    const std::unordered_set<index_t>& N_j = neighbor_sets[j];
    const std::unordered_set<index_t>& N_s = neighbor_sets[s];

    // Find the smallest set for efficiency
    const std::unordered_set<index_t>* smallest_set = &N_i;
    size_t min_size = N_i.size();

    if (N_j.size() < min_size) {
        smallest_set = &N_j;
        min_size = N_j.size();
    }

    if (N_s.size() < min_size) {
        smallest_set = &N_s;
    }

    // Compute triple intersection mass by iterating over smallest set
    // and checking membership in the other two
    double triple_intersection_mass = 0.0;

    for (index_t v : *smallest_set) {
        // Check if v is in all three neighborhoods
        bool in_N_i = (smallest_set == &N_i) || (N_i.find(v) != N_i.end());
        bool in_N_j = (smallest_set == &N_j) || (N_j.find(v) != N_j.end());
        bool in_N_s = (smallest_set == &N_s) || (N_s.find(v) != N_s.end());

        if (in_N_i && in_N_j && in_N_s) {
            // v is in the triple intersection, add its density
            triple_intersection_mass += rho.rho[0][v];
        }
    }

    return triple_intersection_mass;
}

// ============================================================
// ITERATION HELPERS (public/private member functions)
// ============================================================

/**
 * @brief Update vertex mass matrix from evolved vertex densities
 *
 * Updates only the vertex mass matrix M₀ from the current vertex densities after
 * density evolution via damped heat diffusion. This function is called during
 * iterative refinement as part of the metric reconstruction process.
 *
 * MATHEMATICAL CONSTRUCTION:
 *
 * The vertex mass matrix is diagonal by mathematical necessity:
 *   M₀ = diag(ρ₀(1), ρ₀(2), ..., ρ₀(n))
 *
 * Vertices have no geometric interaction under the inner product construction,
 * so off-diagonal entries are always zero. Each diagonal entry represents the
 * total density mass in the neighborhood of that vertex.
 *
 * ITERATION CONTEXT:
 *
 * This function is called as Step 3 in the iteration loop:
 *   Step 1: Density diffusion evolves ρ₀ via damped heat equation
 *   Step 2: Edge densities ρ₁ recomputed from evolved ρ₀
 *   Step 3: Vertex mass matrix update ← THIS FUNCTION (updates only M₀)
 *   Step 4: Edge mass matrix construction (builds M₁ from ρ₀ and ρ₁)
 *   Step 5: Response-coherence modulation (modulates M₁ directly)
 *   Step 6: Laplacian reassembly with updated M₀ and modulated M₁
 *   Step 7: Response smoothing
 *
 * NOTE: This function does NOT modify M₁. The edge mass matrix is rebuilt
 * fresh in Step 4, then modulated in Step 5. Only the vertex mass matrix M₀
 * is updated by this function.
 *
 * REGULARIZATION:
 *
 * Each diagonal entry is regularized to ensure strict positivity:
 *   M₀[i,i] = max(ρ₀(i), 1e-15)
 *
 * This prevents numerical issues in Laplacian assembly where M₀ appears in
 * products and inverses. The threshold 1e-15 is small enough to not affect
 * typical density values but large enough to avoid underflow.
 *
 * FACTORIZATION:
 *
 * The vertex mass matrix is diagonal, so no Cholesky factorization is needed
 * or stored. The solver pointer g.M_solver[0] is reset to null, indicating
 * that M₀ is handled via element-wise operations rather than matrix
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
 * @pre rho.rho[0] contains evolved vertex densities from damped heat diffusion
 * @pre rho.rho[0].sum() ≈ n (normalized from diffusion step)
 * @pre g.M[0] is allocated as n × n sparse diagonal matrix
 * @pre All entries of rho.rho[0] are non-negative
 *
 * @post g.M[0] diagonal entries updated to current ρ₀ values
 * @post All diagonal entries satisfy M₀[i,i] ≥ 1e-15
 * @post g.M[0] remains diagonal (no off-diagonal entries)
 * @post g.M_solver[0] is null (no factorization for diagonal matrix)
 * @post spectral_cache.is_valid == false (cache invalidated)
 * @post g.M[1] is unchanged (edge mass matrix not affected)
 *
 * @note This function modifies M₀ in place using coeffRef for efficiency.
 *       Since M₀ is diagonal and already allocated, we avoid reconstruction.
 *
 * @note The function does NOT call assemble_operators(). The caller must
 *       explicitly reassemble the Laplacian after completing all metric
 *       updates (both M₀ and M₁) to incorporate the new values into L₀.
 *
 * @see apply_damped_heat_diffusion() for Step 1 (evolves ρ₀)
 * @see update_edge_mass_matrix() for Step 4 (builds M₁)
 * @see apply_response_coherence_modulation() for Step 5 (modulates M₁)
 * @see assemble_operators() for Step 6 (rebuilds L₀ from updated metric)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::update_vertex_metric_from_density() {
    const size_t n_vertices = S[0].size();

    // ============================================================
    // Validation
    // ============================================================

    if (static_cast<size_t>(rho.rho[0].size()) != n_vertices) {
        Rf_error("update_vertex_metric_from_density: vertex density size mismatch");
    }

    if (g.M[0].rows() != static_cast<Eigen::Index>(n_vertices) ||
        g.M[0].cols() != static_cast<Eigen::Index>(n_vertices)) {
        Rf_error("update_vertex_metric_from_density: M[0] dimension mismatch");
    }

    // ============================================================
    // Update diagonal entries from evolved vertex densities
    // ============================================================

    for (size_t i = 0; i < n_vertices; ++i) {
        // Apply regularization to ensure strict positivity
        double mass = std::max(rho.rho[0][i], 1e-15);

        // Update in place (M₀ is diagonal, so coeffRef is efficient)
        g.M[0].coeffRef(i, i) = mass;
    }

    // ============================================================
    // Clear factorization (not needed for diagonal matrix)
    // ============================================================

    g.M_solver[0].reset();

    // Note: g.M[1] is intentionally NOT modified here
    // The edge mass matrix was already modulated by response-coherence
    // in Step 3 and should not be overwritten
}

/**
 * @brief Update edge mass matrix from evolved vertex densities
 *
 * Recomputes the edge mass matrix M₁ using current vertex and edge densities
 * after density evolution. This function performs the same triple intersection
 * computations as compute_edge_mass_matrix() but operates in the context
 * of iterative refinement where densities have evolved from their initial
 * values through damped heat diffusion.
 *
 * The edge mass matrix encodes geometric relationships between edges through
 * their shared vertices and overlapping neighborhoods. For edges e_ij = [i,j]
 * and e_is = [i,s] sharing vertex v_i, the inner product is:
 *
 *   ⟨e_ij, e_is⟩ = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j) ∩ N̂_k(x_s)} ρ₀(v)
 *
 * As vertex densities ρ₀ evolve during iteration, these inner products change,
 * reflecting how the Riemannian geometry adapts to concentrate mass in
 * response-coherent regions and deplete mass across response boundaries. The
 * updated mass matrix captures these evolved geometric relationships.
 *
 * ITERATION CONTEXT:
 * This function is called as Step 4 in the iteration loop. At this point:
 *   - Vertex densities ρ₀ have been evolved via damped heat diffusion
 *   - Edge densities ρ₁ have been recomputed from evolved ρ₀
 *   - Vertex mass matrix M₀ has been updated from evolved ρ₀
 *   - Response-coherence modulation has NOT yet been applied (comes next)
 *   - The combinatorial structure (S[1], neighbor_sets, S[2]) remains fixed
 *
 * The recomputed M₁ will be modulated by response-coherence in the next step,
 * then used to reassemble the vertex Laplacian:
 *   L₀ = B₁ M₁⁻¹ B₁ᵀ M₀
 * which in turn drives the next iteration's response smoothing and density evolution.
 *
 * IMPLEMENTATION STRATEGY:
 * The current implementation performs full recomputation by delegating to
 * compute_edge_mass_matrix(), which iterates over all triangles and computes
 * inner products for all edge pairs. This ensures correctness and maintains
 * consistency with the initialization logic.
 *
 * The function rebuilds M₁ completely by:
 *   1. Setting diagonal entries from edge densities: M₁[e,e] = ρ₁(e)
 *   2. Computing off-diagonal entries via triple intersections using ρ₀
 *   3. Iterating over triangles to find all edge pairs sharing vertices
 *   4. Assembling the symmetric sparse matrix from triplets
 *   5. Applying regularization to diagonal entries
 *
 * COMPUTATIONAL COMPLEXITY:
 * O(n·k²) where n is the number of vertices and k is the neighborhood size.
 * This matches the initialization cost since the combinatorial structure is
 * unchanged. For kNN graphs, each vertex has O(k) incident edges, so each
 * vertex contributes O(k²) edge pairs to examine. Computing each triple
 * intersection costs O(k) via the optimized smallest-set-first strategy.
 *
 * PERFORMANCE CONSIDERATIONS:
 * This function is the computational bottleneck of each iteration. For large
 * graphs (n > 10000) or large neighborhoods (k > 50), the O(n·k²) cost can
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
 * 3. **Low-rank updates**: If ||ρ₀^(ℓ) - ρ₀^(ℓ-1)|| is small, approximate
 *    M₁^(ℓ) ≈ M₁^(ℓ-1) + ΔM where ΔM is computed from density perturbations.
 *
 * 4. **Parallel computation**: Use OpenMP to parallelize the outer loop over
 *    vertices, with careful handling of concurrent triplet list updates.
 *
 * These optimizations trade code complexity for computational savings and should
 * be considered if profiling identifies this function as the primary bottleneck.
 *
 * NUMERICAL STABILITY:
 * The function inherits regularization from compute_edge_mass_matrix():
 *   - Diagonal entries: max(ρ₁([i,j]), 1e-15) prevents singular matrices
 *   - Off-diagonal entries: naturally non-negative from density summations
 *   - Symmetry: enforced by adding both (i,j) and (j,i) triplets
 *
 * The resulting matrix is guaranteed symmetric positive semidefinite by
 * construction, as all inner products arise from L²(ρ) pairings.
 *
 * @pre rho.rho[0] contains evolved vertex densities (current iteration)
 * @pre rho.rho[1] contains updated edge densities (current iteration)
 * @pre S[1] contains edge simplex table (unchanged from initialization)
 * @pre S[2] contains triangle simplex table (unchanged from initialization)
 * @pre neighbor_sets contains kNN neighborhoods (unchanged from initialization)
 * @pre g.M[1] is allocated with dimensions n_edges × n_edges
 *
 * @post g.M[1] contains updated edge mass matrix computed from current ρ₀ and ρ₁
 * @post g.M[1] is symmetric: M₁[i,j] == M₁[j,i] for all i,j
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
 * @warning The current implementation performs full O(n·k²) recomputation at
 *          each call. For very large graphs, consider profiling and implementing
 *          incremental update strategies if this becomes a bottleneck.
 *
 * @see compute_edge_mass_matrix() for the initial mass matrix construction
 * @see update_vertex_metric_from_density() for M₀ update (previous step)
 * @see apply_response_coherence_modulation() for M₁ modulation (next step)
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
 * Computes the eigendecomposition L[0] = V Λ V^T and caches the results
 * for use in spectral filtering and parameter selection. The eigenvalues
 * are sorted in ascending order, so eigenvalues[0] ≈ 0 (constant function)
 * and eigenvalues[1] = λ₂ is the spectral gap.
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
 * 1. Diffusion time t scales inversely with spectral gap λ₂
 * 2. Damping parameter β maintains fixed ratio with t
 * 3. User-provided values (> 0) are respected and not overridden
 * 4. Safety warnings for extreme parameter combinations
 *
 * @param t_diffusion Reference to diffusion time parameter (input/output).
 *                    If <= 0 on input, automatically set to 0.5/λ₂.
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
            Rprintf("Computing spectral decomposition for parameter selection...\n");
        }
        compute_spectral_decomposition();
    }

    const double lambda_2 = spectral_cache.lambda_2;

    // Auto-select t_diffusion if not provided by user
    bool t_auto_selected = false;
    if (t_diffusion <= 0.0) {
        // Use moderate default: 0.5/λ₂
        t_diffusion = 0.5 / lambda_2;
        t_auto_selected = true;

        if (verbose) {
            Rprintf("Auto-selected t_diffusion = %.6f (based on spectral gap λ₂ = %.6f)\n",
                    t_diffusion, lambda_2);
            Rprintf("  Conservative: %.6f, Moderate: %.6f, Aggressive: %.6f\n",
                    0.1 / lambda_2, 0.5 / lambda_2, 1.0 / lambda_2);
        }
    } else if (verbose) {
        Rprintf("Using user-provided t_diffusion = %.6f\n", t_diffusion);
    }

    // Auto-select beta_damping if not provided by user
    bool beta_auto_selected = false;
    if (beta_damping <= 0.0) {
        // Coordinate with t: damping is 10% of diffusion scale
        beta_damping = 0.1 / t_diffusion;
        beta_auto_selected = true;

        if (verbose) {
            Rprintf("Auto-selected beta_damping = %.6f (ratio β·t = %.3f)\n",
                    beta_damping, beta_damping * t_diffusion);
        }
    } else if (verbose) {
        Rprintf("Using user-provided beta_damping = %.6f\n", beta_damping);
    }

    // Diagnostic checks and warnings
    const double diffusion_scale = t_diffusion * lambda_2;

    if (diffusion_scale > 3.0) {
        Rf_warning("Large diffusion scale (t·λ₂ = %.2f): density may change dramatically per iteration. "
                   "Consider reducing t_diffusion to %.6f for more conservative updates.",
                   diffusion_scale, 1.0 / lambda_2);
    }

    if (diffusion_scale < 0.05) {
        Rf_warning("Small diffusion scale (t·λ₂ = %.2f): convergence may be very slow. "
                   "Consider increasing t_diffusion to %.6f for faster updates.",
                   diffusion_scale, 0.3 / lambda_2);
    }

    const double damping_ratio = beta_damping * t_diffusion;

    if (damping_ratio > 0.5) {
        Rf_warning("High damping ratio (β·t = %.2f): may over-suppress geometric structure. "
                   "Typical range is [0.05, 0.2].", damping_ratio);
    }

    if (damping_ratio < 0.01) {
        Rf_warning("Low damping ratio (β·t = %.3f): density may collapse onto small regions. "
                   "Consider increasing beta_damping to %.6f.",
                   damping_ratio, 0.05 / t_diffusion);
    }

    // Final summary if verbose
    if (verbose) {
        Rprintf("\nDiffusion parameter summary:\n");
        Rprintf("  λ₂ (spectral gap):  %.6f\n", lambda_2);
        Rprintf("  t (diffusion time): %.6f %s\n", t_diffusion,
                t_auto_selected ? "[auto]" : "[user]");
        Rprintf("  β (damping):        %.6f %s\n", beta_damping,
                beta_auto_selected ? "[auto]" : "[user]");
        Rprintf("  Diffusion scale:    %.3f (t·λ₂)\n", diffusion_scale);
        Rprintf("  Damping ratio:      %.3f (β·t)\n", damping_ratio);
    }
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
 *   ∂ρ/∂t = -L₀ρ - β(ρ - u)
 * where L₀ is the vertex Laplacian, β ≥ 0 controls damping strength, and
 * u = (1, 1, ..., 1)ᵀ represents uniform distribution scaled to sum to n.
 *
 * The heat diffusion term -L₀ρ drives mass toward densely connected regions,
 * as the Laplacian naturally smooths density along the graph structure. The
 * damping term -β(ρ - u) prevents runaway concentration by continuously pulling
 * the distribution back toward uniformity. This balance allows the method to
 * discover meaningful geometric structure without collapsing onto a small set
 * of vertices.
 *
 * We discretize using implicit Euler with step size t:
 *   (I + t(L₀ + βI))ρ_new = ρ_old + tβu
 * The implicit scheme guarantees stability for arbitrarily large t, unlike
 * explicit methods which require restrictive step size bounds. After solving,
 * we renormalize to enforce the constraint Σρ_new(i) = n.
 *
 * The system matrix A = I + t(L₀ + βI) is symmetric positive definite, as it
 * combines the identity with positive multiples of L₀ (positive semidefinite)
 * and βI (positive definite for β > 0). We solve using conjugate gradient
 * iteration, which exploits sparsity and symmetry for efficiency.
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
 * @note For kNN graphs with n vertices, typical values are t ∈ [0.1/λ₂, 1.0/λ₂]
 *       where λ₂ is the spectral gap, and β ≈ 0.1/t. These choices balance
 *       geometric smoothing with stability.
 *
 * @note The solver uses tolerance 1e-10 and maximum 1000 iterations. For large
 *       systems (n > 10000), consider using a preconditioner or iterative
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
    // Part 2: Build system matrix A = I + t(L₀ + βI)
    // ========================================================================

    // Start with identity matrix
    spmat_t I(n, n);
    I.setIdentity();

    // Build system matrix: A = I + t*L₀ + t*β*I = (1 + t*β)*I + t*L₀
    // We compute this efficiently by scaling operations
    const double identity_coeff = 1.0 + t * beta;
    spmat_t A = identity_coeff * I + t * L.L[0];

    // Compress for efficient solve
    A.makeCompressed();

    // ========================================================================
    // Part 3: Build right-hand side b = ρ_current + tβu
    // ========================================================================

    // Uniform distribution vector u = (1, 1, ..., 1)ᵀ
    vec_t u = vec_t::Ones(n);

    // Right-hand side: b = ρ_current + t*β*u
    vec_t b = rho_current + (t * beta) * u;

    // ========================================================================
    // Part 4: Solve linear system A*ρ_new = b using conjugate gradient
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
 * density diffusion. After vertex densities ρ₀ have evolved through the damped
 * heat equation, edge densities ρ₁ must be updated to maintain consistency with
 * the evolved vertex distribution. The edge density for edge [i,j] aggregates
 * vertex densities over the neighborhood intersection:
 *
 *   ρ₁([i,j]) = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j)} ρ₀(v)
 *
 * This construction ensures that edge densities reflect the current geometric
 * distribution of vertex mass. Edges connecting vertices with overlapping
 * high-density neighborhoods receive high density, while edges spanning sparse
 * or disconnected regions receive low density.
 *
 * ITERATION CONTEXT:
 * Within each iteration of fit_knn_riem_graph_regression():
 *   Step 1: Density diffusion evolves ρ₀ via damped heat equation
 *   Step 2: Edge density update ← THIS FUNCTION
 *   Step 3: Vertex mass matrix update (builds M₀ from ρ₀)
 *   Step 4: Edge mass matrix construction (builds M₁ from ρ₀ and ρ₁)
 *   Step 5: Response-coherence modulation (modulates M₁, not ρ₁)
 *   Step 6: Laplacian reassembly
 *   Step 7: Response smoothing
 *
 * The updated edge densities serve two purposes:
 * 1. Diagonal entries of the edge mass matrix M₁ during metric construction,
 *    since ⟨e_ij, e_ij⟩ = ρ₁([i,j]) by construction
 * 2. Geometric interpretation: relative density of edges in the complex
 *
 * IMPORTANT: The edge densities ρ₁ computed by this function remain unchanged
 * for the rest of the iteration. Response-coherence modulation (Step 5)
 * operates on the mass matrix M₁, not on the densities ρ₁. The densities are
 * pure geometric quantities derived from evolved vertex distributions, while
 * the mass matrix incorporates response-aware modulation.
 *
 * NORMALIZATION STRATEGY:
 * After computing raw edge densities from pairwise intersections, the function
 * normalizes to sum to n_edges (the number of edges). This normalization ensures:
 *   - The average edge has density 1
 *   - Total edge mass remains stable across iterations
 *   - Numerical conditioning of the mass matrix M₁
 *   - Interpretability: edge densities represent relative mass in standardized units
 *
 * The normalization preserves relative density differences while maintaining
 * a fixed total mass scale, preventing geometric drift over iterations.
 *
 * COMPUTATIONAL COMPLEXITY:
 * O(n_edges · k) where k is the average neighborhood size. For each edge,
 * computing the pairwise intersection requires iterating over one neighborhood
 * and checking membership in the other. The optimized strategy iterates over
 * the smaller neighborhood, giving O(min(|N_i|, |N_j|)) per edge.
 *
 * For kNN graphs with uniform k, this gives O(n_edges · k). Since n_edges ≈ n·k
 * for kNN graphs, the total complexity is O(n·k²), though with a smaller
 * constant than the full mass matrix computation (which involves triple
 * intersections and O(k²) edge pairs per vertex).
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
 *   - Iteration uses evolved rho.rho[0] values for aggregation
 *
 * Both produce edge densities normalized to sum to n_edges.
 *
 * RELATION TO MASS MATRIX:
 * The updated edge densities populate the diagonal of the edge mass matrix M₁.
 * For edge e_ij = [i,j], the diagonal entry is:
 *   M₁[e_ij, e_ij] = ρ₁([i,j])
 *
 * This relationship holds because the edge self-inner product equals the
 * pairwise intersection mass:
 *   ⟨e_ij, e_ij⟩ = Σ_{v ∈ N_i ∩ N_j} ρ₀(v) = ρ₁([i,j])
 *
 * Off-diagonal entries of M₁ involve triple intersections and are computed
 * separately via compute_edge_mass_matrix().
 *
 * @pre rho.rho[0] contains evolved vertex densities from current iteration
 * @pre S[1] contains edge simplex table (unchanged from initialization)
 * @pre neighbor_sets contains kNN neighborhoods (unchanged from initialization)
 * @pre rho.rho[1] is allocated with size equal to number of edges
 *
 * @post rho.rho[1] contains updated edge densities computed from current ρ₀
 * @post rho.rho[1].sum() ≈ n_edges (normalized to sum to number of edges)
 * @post All entries of rho.rho[1] are non-negative
 *
 * @note This function only updates edge densities (dimension 1). Vertex densities
 *       (dimension 0) remain unchanged, as they were already updated by the
 *       damped heat diffusion step.
 *
 * @note The combinatorial structure (which edges exist, vertex neighborhoods)
 *       remains fixed throughout iteration. Only the numerical density values
 *       change based on evolved vertex distributions.
 *
 * @note The edge densities ρ₁ remain unchanged after this function returns.
 *       Response-coherence modulation operates on the mass matrix M₁, not
 *       on the densities themselves.
 *
 * @see compute_initial_densities() for the initialization version
 * @see update_edge_mass_matrix() for how ρ₁ enters M₁ diagonal
 * @see apply_response_coherence_modulation() for M₁ modulation (not ρ₁)
 */
void riem_dcx_t::update_edge_densities_from_vertices() {
    const size_t n_edges = S[1].size();

    // Ensure vertex densities are available
    if (static_cast<size_t>(rho.rho[0].size()) != S[0].size()) {
        Rf_error("update_edge_densities_from_vertices: vertex densities not properly initialized");
    }

    // ============================================================
    // Compute edge densities from pairwise intersections
    // ============================================================

    // For each edge [i,j], aggregate vertex densities over the
    // intersection of neighborhoods: N̂_k(x_i) ∩ N̂_k(x_j)

    for (size_t e = 0; e < n_edges; ++e) {
        // Get the two vertices of this edge
        const std::vector<index_t>& edge_verts = S[1].simplex_verts[e];
        const index_t i = edge_verts[0];
        const index_t j = edge_verts[1];

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
                edge_density += rho.rho[0][v];
            }
        }

        rho.rho[1][e] = edge_density;
    }

    // ============================================================
    // Normalize edge densities to sum to n_edges
    // ============================================================

    double edge_density_sum = rho.rho[1].sum();

    if (edge_density_sum > 1e-15) {
        // Normal case: scale to preserve relative densities while fixing total mass
        rho.rho[1] *= static_cast<double>(n_edges) / edge_density_sum;
    } else {
        // Degenerate case: all intersections are empty or near-zero
        // Fall back to uniform density to maintain positive definiteness
        rho.rho[1].setConstant(1.0);

        Rf_warning("update_edge_densities_from_vertices: edge density sum near zero, "
                   "falling back to uniform density");
    }
}

/**
 * @brief Apply response-coherence modulation to existing edge mass matrix
 *
 * Modulates both diagonal and off-diagonal entries of the edge mass matrix M₁
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
 *   Δ_ij = |ŷ(i) - ŷ(j)| for edge [i,j]
 *   σ₁ = IQR({Δ_ij : all edges})
 *   M₁[e,e] ← M₁[e,e] · Γ(Δ_ij/σ₁)
 *
 * For off-diagonal entries (edge pair inner products):
 *   For triangle τ = [i,j,s] ∈ S[2]:
 *     Δ_τ = max{|ŷ(i)-ŷ(j)|, |ŷ(i)-ŷ(s)|, |ŷ(j)-ŷ(s)|}
 *     σ₂ = IQR({Δ_τ : all triangles})
 *     M₁[e_ij, e_is] ← M₁[e_ij, e_is] · Γ(Δ_τ/σ₂)
 *
 * where Γ(x) = (1 + x²)^(-γ) is the penalty function.
 *
 * The use of two separate scales (σ₁ for edges, σ₂ for triangles) allows the
 * modulation to adapt independently to the distribution of response variation
 * at different geometric levels. Diagonal entries use pairwise differences along
 * edges, while off-diagonal entries use the maximum variation across triangle
 * vertices, reflecting that off-diagonal interactions involve three points
 * rather than two.
 *
 * IMPORTANT OPTIMIZATION:
 * Rather than rebuilding M₁ from scratch, this version:
 * 1. Computes modulation factors for all edges and triangles
 * 2. Applies factors directly to existing matrix entries via coeffRef
 * 3. Preserves matrix structure and sparsity pattern
 *
 * This is much faster than reconstruction, especially for large graphs.
 *
 * NORMALIZATION:
 * After modulation, the total Frobenius mass is normalized to equal the
 * pre-modulation mass:
 *   ||M₁^(modulated)||_F = ||M₁^(original)||_F
 *
 * This prevents systematic drift in the overall geometry scale across iterations.
 * The normalization preserves relative modulation differences while fixing the
 * total mass.
 *
 * ITERATION CONTEXT:
 * This function is called as Step 5 in the iteration loop:
 *   Step 1: Density diffusion (evolves ρ₀)
 *   Step 2: Edge density update (computes ρ₁ from evolved ρ₀)
 *   Step 3: Vertex mass matrix update (builds M₀ from ρ₀)
 *   Step 4: Edge mass matrix construction (builds M₁ from ρ₀ and ρ₁)
 *   Step 5: Response-coherence modulation ← THIS FUNCTION
 *   Step 6: Laplacian reassembly (builds L₀ from M₀ and modulated M₁)
 *   Step 7: Response smoothing
 *
 * The modulated M₁ is used immediately in Step 6 to assemble the Laplacian,
 * which then drives response smoothing and the next iteration's density evolution.
 *
 * @param y_hat Current fitted response values
 * @param gamma Decay rate parameter (typically 0.5-2.0). Controls the strength
 *              of modulation: larger γ produces sharper response boundaries,
 *              while smaller γ yields gentler transitions.
 *
 * @pre y_hat.size() == S[0].size()
 * @pre g.M[1] is assembled with current edge mass matrix from Step 4
 * @pre S[2] contains triangles (may be empty if pmax < 2, in which case only
 *      diagonal modulation is performed)
 * @pre gamma > 0
 *
 * @post g.M[1] modulated in place with diagonal and off-diagonal entries adjusted
 * @post g.M[1] remains symmetric positive semidefinite
 * @post Total Frobenius mass preserved: ||M₁||_F unchanged by normalization
 * @post g.M_solver[1] invalidated (factorization no longer valid)
 *
 * @note If S[2] is empty (no triangles), only diagonal entries are modulated.
 *       A warning is issued if γ > 0 but no triangles are available.
 *
 * @note The two-scale approach (σ₁ ≠ σ₂) is essential for proper modulation.
 *       Diagonal and off-diagonal entries use different response variation
 *       measures (pairwise vs. max over triangle) and thus require different
 *       normalization scales.
 *
 * @note This function operates on the **mass matrix** M₁, not on the edge
 *       densities ρ₁. The densities remain unchanged; only the geometric
 *       structure (inner products) is modulated.
 *
 * @see update_edge_mass_matrix() for Step 4 (builds M₁ before modulation)
 * @see assemble_operators() for Step 6 (uses modulated M₁ to build L₀)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::apply_response_coherence_modulation(
    const vec_t& y_hat,
    double gamma
) {
    const size_t n_edges = S[1].size();
    const size_t n_vertices = S[0].size();

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
        const std::vector<index_t>& verts = S[1].simplex_verts[e];
        double delta = std::abs(y_hat[verts[0]] - y_hat[verts[1]]);
        edge_deltas.push_back(delta);
    }

    // ============================================================
    // Step 2: Compute scale parameter σ₁ for edges
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

    if (S.size() > 2 && S[2].size() > 0) {
        std::vector<double> triangle_deltas;
        triangle_deltas.reserve(S[2].size());

        for (size_t t = 0; t < S[2].size(); ++t) {
            const auto& verts = S[2].simplex_verts[t];

            // Compute maximum response variation over triangle vertices
            double delta = std::max({
                std::abs(y_hat[verts[0]] - y_hat[verts[1]]),
                std::abs(y_hat[verts[0]] - y_hat[verts[2]]),
                std::abs(y_hat[verts[1]] - y_hat[verts[2]])
            });
            triangle_deltas.push_back(delta);

            // Get the three edge IDs for this triangle
            std::vector<index_t> e01_verts = {verts[0], verts[1]};
            std::vector<index_t> e02_verts = {verts[0], verts[2]};
            std::vector<index_t> e12_verts = {verts[1], verts[2]};

            index_t e01 = S[1].get_id(e01_verts);
            index_t e02 = S[1].get_id(e02_verts);
            index_t e12 = S[1].get_id(e12_verts);

            // Record maximum delta for each edge pair in this triangle
            uint64_t key_01_02 = pack_edge_pair(e01, e02);
            uint64_t key_01_12 = pack_edge_pair(e01, e12);
            uint64_t key_02_12 = pack_edge_pair(e02, e12);

            edge_pair_max_delta[key_01_02] = std::max(
                edge_pair_max_delta[key_01_02], delta);
            edge_pair_max_delta[key_01_12] = std::max(
                edge_pair_max_delta[key_01_12], delta);
            edge_pair_max_delta[key_02_12] = std::max(
                edge_pair_max_delta[key_02_12], delta);
        }

        // Compute scale parameter σ₂ for triangles
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

        // Iterate over all non-zero entries in M₁
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
            Rprintf("Note: No triangles available...\n");
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
            "  - Scale parameters: σ₁ = %.3e, σ₂ = %.3e\n"
            "  - All edges have large response variation: min(Δ) = %.3e, max(Δ) = %.3e\n"
            "Consider: (1) reducing gamma, (2) checking response data for outliers,\n"
            "or (3) setting gamma = 0 to disable modulation.",
            gamma, sigma_1, sigma_2,
            *std::min_element(edge_deltas.begin(), edge_deltas.end()),
            *std::max_element(edge_deltas.begin(), edge_deltas.end())
            );

        // Rf_warning(
        //     "Response-coherence modulation collapsed all mass (gamma = %.3f too large). "
        //     "Disabling modulation for this iteration. Reduce gamma parameter.",
        //     gamma
        // );
        // // Restore to pre-modulation state by rebuilding from densities
        // compute_edge_mass_matrix();
        // return;  // Skip normalization, use unmodulated matrix
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
 * GCV score = ||y - y_hat||² / (n - trace(S))²
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
        // For exponential filters: -log(eps)/lambda_max
        eta_max = (lambda_max > 0) ? -std::log(eps) / lambda_max : 1.0;
        break;

    case rdcx_filter_type_t::TIKHONOV:
    case rdcx_filter_type_t::CUBIC_SPLINE:
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
 *        - HEAT_KERNEL: exp(-ηλ), exponential decay
 *        - TIKHONOV: 1/(1+ηλ), rational decay
 *        - CUBIC_SPLINE: 1/(1+ηλ²), spline smoothness
 *        - GAUSSIAN: exp(-ηλ²), aggressive high-frequency attenuation
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
    result.y_hat = best_y_hat;
    result.gcv_scores = gcv_scores;
    result.eta_grid = eta_grid;

    return result;
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
convergence_status_t riem_dcx_t::check_convergence(
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
detailed_convergence_status_t riem_dcx_t::check_convergence_detailed(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    const std::vector<vec_t>& rho_prev,
    const std::vector<vec_t>& rho_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations,
    const std::vector<double>& response_change_history = {}
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
 * 0.5-2.0. Larger gamma creates stronger geometric adaptation to response structure.
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
 * complex.init_dims(1, {n_points, 0});  // Max dimension 1, n points
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
 * @pre Complex must be initialized via init_dims() before calling
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
    int test_stage
) {
    // ================================================================
    // PART I: INITIALIZATION
    // ================================================================

    // Phase 1-4: Build geometric structure
    initialize_from_knn(
        X, k,
        use_counting_measure,
        density_normalization,
        max_ratio_threshold,
        threshold_percentile
    );

    if (test_stage == 0) {
        Rprintf("TEST_STAGE 0: Stopped after initialize_from_knn\n");
        return;
    }

    // Validate triangle construction for response-coherence
    if (gamma_modulation > 0.0 && (S.size() <= 2 || S[2].size() == 0)) {
        Rf_warning("gamma_modulation = %.2f requires triangles for off-diagonal "
                   "modulation, but no triangles were constructed. "
                   "Only diagonal entries of M₁ will be modulated. "
                   "Consider: (1) increasing k to create more triangles, "
                   "(2) setting gamma_modulation = 0 to disable modulation, or "
                   "(3) adjusting pruning thresholds.",
                   gamma_modulation);
    }

    // Phase 4.5: Select or validate diffusion parameters
    select_diffusion_parameters(t_diffusion, beta_damping, /*verbose=*/true);

    if (test_stage == 1) {
        Rprintf("TEST_STAGE 1: Stopped after parameter selection\n");
        return;
    }

    // Phase 5: Initial response smoothing
    auto gcv_result = smooth_response_via_spectral_filter(
        y, n_eigenpairs, filter_type
    );
    vec_t y_hat_curr = gcv_result.y_hat;
    sig.y_hat_hist.clear();
    sig.y_hat_hist.push_back(y_hat_curr);

    if (test_stage == 2) {
        Rprintf("TEST_STAGE 2: Stopped after initial smoothing\n");
        return;
    }

    // ================================================================
    // PART II: ITERATIVE REFINEMENT
    // ================================================================

    vec_t y_hat_prev;
    std::vector<vec_t> rho_prev;
    std::vector<double> response_change_history;

    for (int iter = 1; iter <= max_iterations; ++iter) {
        y_hat_prev = y_hat_curr;
        rho_prev = rho.rho;

        // Step 1: Density diffusion
        rho.rho[0] = apply_damped_heat_diffusion(
            rho.rho[0], t_diffusion, beta_damping
        );

        if (test_stage == 3) {
            Rprintf("TEST_STAGE 3: Stopped after first diffusion\n");
            return;
        }

        // Step 2: Edge density update
        update_edge_densities_from_vertices();

        if (test_stage == 4) {
            Rprintf("TEST_STAGE 4: Stopped after edge density update\n");
            return;
        }

        // Step 3: Update vertex mass matrix from evolved densities
        update_vertex_metric_from_density();  // Updates only M[0]

        if (test_stage == 5) {
            Rprintf("TEST_STAGE 5: Stopped after update_vertex_metric_from_density()\n");
            return;
        }

        // Step 4: Rebuild edge mass matrix from evolved densities
        update_edge_mass_matrix();  // Compute fresh M₁ from current ρ₀

        if (test_stage == 6) {
            Rprintf("TEST_STAGE 6: Stopped after update_edge_mass_matrix()\n");
            return;
        }

        // Step 5: Apply response-coherence modulation to fresh M₁
        apply_response_coherence_modulation(y_hat_curr, gamma_modulation);

        if (test_stage == 7) {
            Rprintf("TEST_STAGE 7: Stopped after apply_response_coherence_modulation()\n");
            return;
        }

        // Step 6: Laplacian reassembly
        assemble_operators();
        spectral_cache.invalidate();

        if (test_stage == 8) {
            Rprintf("TEST_STAGE 8: Stopped after assemble_operators()\n");
            return;
        }

        // Step 7: Response smoothing
        gcv_result = smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
        );
        y_hat_curr = gcv_result.y_hat;
        sig.y_hat_hist.push_back(y_hat_curr);

        if (test_stage == 9) {
            Rprintf("TEST_STAGE 9: Stopped after smooth_response_via_spectral_filter()\n");
            return;
        }

        // Step 8: Convergence check
        auto status = check_convergence_detailed(
            y_hat_prev, y_hat_curr,
            rho_prev, rho.rho,
            epsilon_y, epsilon_rho,
            iter, max_iterations,
            response_change_history
        );

        if (test_stage == 10) {
            Rprintf("TEST_STAGE 10: Stopped after check_convergence_detailed()\n");
            return;
        }

        response_change_history.push_back(status.response_change);

        Rprintf("Iteration %d: response_change=%.6f, density_change=%.6f\n",
                iter, status.response_change, status.max_density_change);

        for (size_t p = 0; p < status.density_changes_by_dim.size(); ++p) {
            Rprintf("  Density[%d] change: %.6f\n",
                    (int)p, status.density_changes_by_dim[p]);
        }
        if (status.estimated_iterations_remaining > 0) {
            Rprintf("  Estimated iterations remaining: %d\n",
                    status.estimated_iterations_remaining);
        }

        if (status.converged) {
            Rprintf("%s\n", status.message.c_str());
            break;
        }

        ///// testing
        if (iter == 1) break;  // Only do one iteration for testing
    }

    // ================================================================
    // PART III: FINALIZATION
    // ================================================================

    sig.y = y;
}
