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
// INTERNAL HELPERS (anonymous namespace)
// ============================================================

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
    );

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
 * @pre reference_measure must be initialized and have size equal to number of vertices
 * @post rho.rho[0] contains normalized vertex densities summing to n
 * @post rho.rho[1] contains normalized edge densities summing to n_edges
 */
void riem_dcx_t::compute_initial_densities() {
    const size_t n_vertices = S[0].size();
    const size_t n_edges = S[1].size();

    // Ensure reference measure is available
    if (reference_measure.size() != n_vertices) {
        Rf_error("Reference measure not initialized or has incorrect size");
    }

    // ============================================================
    // Part 1: Compute vertex densities
    // ============================================================

    // For each vertex i, sum the reference measure over its k-neighborhood
    rho.rho[0] = vec_t::Zero(n_vertices);

    for (size_t i = 0; i < n_vertices; ++i) {
        double vertex_density = 0.0;

        // Aggregate measure over all neighbors j ∈ Ň_k(x_i)
        for (index_t j : neighbor_sets[i]) {
            vertex_density += reference_measure[j];
        }

        rho.rho[0][i] = vertex_density;
    }

    // Normalize vertex densities to sum to n
    double vertex_density_sum = rho.rho[0].sum();

    if (vertex_density_sum > 1e-15) {
        rho.rho[0] *= static_cast<double>(n_vertices) / vertex_density_sum;
    } else {
        // Fallback: uniform density if something went wrong
        rho.rho[0].setConstant(1.0);
    }

    // ============================================================
    // Part 2: Compute edge densities
    // ============================================================

    // For each edge [i,j], sum the reference measure over the intersection
    // of neighborhoods: Ň_k(x_i) ∩ Ň_k(x_j)
    rho.rho[1] = vec_t::Zero(n_edges);

    for (size_t e = 0; e < n_edges; ++e) {
        // Get the two vertices of this edge
        const std::vector<index_t>& edge_verts = S[1].simplex_verts[e];
        const index_t i = edge_verts[0];
        const index_t j = edge_verts[1];

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

        rho.rho[1][e] = edge_density;
    }

    // Normalize edge densities to sum to n_edges
    double edge_density_sum = rho.rho[1].sum();

    if (edge_density_sum > 1e-15) {
        rho.rho[1] *= static_cast<double>(n_edges) / edge_density_sum;
    } else {
        // Fallback: uniform density if intersections are all empty
        rho.rho[1].setConstant(1.0);
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
 * all pairs of edges through triple neighborhood intersections. For edges
 * e_ij = [i,j] and e_is = [i,s] sharing vertex v_i, the inner product is:
 *
 *   ⟨e_ij, e_is⟩ = Σ_{v ∈ N̂_k(x_i) ∩ N̂_k(x_j) ∩ N̂_k(x_s)} ρ₀(v)
 *
 * This measures the total vertex density in the triple intersection of
 * neighborhoods, encoding the geometric relationship between edges through
 * their shared vertex and overlapping neighborhoods.
 *
 * The resulting matrix is symmetric positive semidefinite and sparse. For
 * kNN complexes with parameter k, each edge typically interacts with O(k²)
 * other edges, making sparse storage efficient.
 *
 * COMPUTATIONAL COMPLEXITY: O(n * k²) where n is number of vertices and k is
 * the neighborhood size. This is the bottleneck operation in metric construction.
 *
 * @pre S[0] and S[1] must be populated with vertices and edges
 * @pre rho.rho[0] must contain current vertex densities
 * @pre stars[0] must be populated (maps vertices to their incident edges)
 * @post g.M[1] contains symmetric positive semidefinite edge mass matrix
 */
void riem_dcx_t::compute_edge_mass_matrix() {
    const size_t n_edges = S[1].size();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_edges * 10);

    const size_t n_vertices = S[0].size();

    for (size_t i = 0; i < n_vertices; ++i) {
        const std::vector<index_t>& incident_edges = stars[0].star_over[i];
        const size_t n_incident = incident_edges.size();

        if (n_incident == 0) continue;

        for (size_t a = 0; a < n_incident; ++a) {
            const index_t e_ij_idx = incident_edges[a];
            const std::vector<index_t>& e_ij_verts = S[1].simplex_verts[e_ij_idx];

            // Only process diagonal when we're at the SMALLER endpoint
            // Since edges are stored sorted: e_ij_verts[0] < e_ij_verts[1]
            if (i == e_ij_verts[0]) {
                // Add diagonal entry (only once per edge)
                double diagonal = std::max(rho.rho[1][e_ij_idx], 1e-15);
                triplets.emplace_back(e_ij_idx, e_ij_idx, diagonal);
            }

            // Off-diagonal entries (process all pairs as before)
            for (size_t b = a + 1; b < n_incident; ++b) {
                const index_t e_is_idx = incident_edges[b];

                double inner_product = compute_edge_inner_product(e_ij_idx, e_is_idx, i);

                triplets.emplace_back(e_ij_idx, e_is_idx, inner_product);
                triplets.emplace_back(e_is_idx, e_ij_idx, inner_product);
            }
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
 * @brief Update metric from evolved densities during iterative refinement
 *
 * Reconstructs the Riemannian metric structure g from the current density
 * distribution ρ after density evolution steps. This function is called
 * during each iteration of the regression algorithm after vertex densities
 * have been updated via damped heat diffusion and edge densities have been
 * recomputed from the evolved vertex distribution.
 *
 * The metric update follows the same mathematical construction as initialization
 * but operates on evolved densities rather than initial values. The vertex mass
 * matrix M₀ is updated by replacing diagonal entries with current vertex densities.
 * The edge mass matrix M₁ is recomputed via triple neighborhood intersections
 * using the evolved vertex densities, capturing how the geometric relationships
 * between edges change as density concentrates in response-coherent regions and
 * depletes across response boundaries.
 *
 * A critical side effect of metric updating is invalidation of the spectral cache.
 * Since the Laplacian L₀ depends on the mass matrices through L₀ = B₁ M₁⁻¹ B₁ᵀ M₀,
 * any change to the metric renders previously computed eigendecompositions invalid.
 * The function explicitly invalidates the spectral cache to ensure that subsequent
 * response smoothing operations trigger fresh eigendecomposition with the updated
 * geometry.
 *
 * ITERATION CONTEXT:
 * This function is called as Step 4 in the iteration loop of fit_knn_riem_graph_regression():
 *   Step 1: Density diffusion (evolves ρ₀)
 *   Step 2: Edge density update (derives ρ₁ from evolved ρ₀)
 *   Step 3: Response-coherence modulation (adjusts ρ₁ based on response variation)
 *   Step 4: Metric update ← THIS FUNCTION
 *   Step 5: Laplacian reassembly (builds L₀ from updated M₀, M₁)
 *   Step 6: Response smoothing (solves for ŷ using updated L₀)
 *   Step 7: Convergence check
 *
 * COMPUTATIONAL COST:
 * The dominant cost is recomputing M₁ via update_edge_mass_matrix(), which
 * performs O(n·k²) triple intersection computations. Updating M₀ is O(n) and
 * negligible by comparison. For large graphs (n > 10000), this step can become
 * the primary computational bottleneck of each iteration.
 *
 * FUTURE OPTIMIZATION:
 * The current implementation performs full recomputation of M₁ at each iteration.
 * Potential optimizations include:
 *   - Incremental updates tracking which vertex densities changed significantly
 *   - Lazy evaluation deferring M₁ update until Laplacian assembly
 *   - Low-rank approximations when density changes are small
 * These optimizations would require tracking density change history and accepting
 * increased code complexity in exchange for reduced iteration cost.
 *
 * @pre rho.rho[0] contains evolved vertex densities (normalized to sum to n)
 * @pre rho.rho[1] contains updated edge densities (normalized to sum to n_edges)
 * @pre S[0], S[1], neighbor_sets, and stars[0] remain unchanged from initialization
 *
 * @post g.M[0] diagonal entries updated to current vertex densities
 * @post g.M[1] recomputed with current vertex densities via triple intersections
 * @post spectral_cache.is_valid == false (cache invalidated)
 * @post Laplacian L.L[0] remains unchanged (caller must invoke assemble_operators())
 *
 * @note Unlike initialize_metric_from_density(), this function assumes the
 *       sparse matrix structures g.M[0] and g.M[1] are already allocated with
 *       correct dimensions. It updates entries in place rather than reconstructing
 *       from scratch, though the current implementation of update_edge_mass_matrix()
 *       does perform full reconstruction of M₁.
 *
 * @note This function does NOT automatically reassemble the Laplacian. The caller
 *       must explicitly call assemble_operators() after metric update to rebuild
 *       L₀ from the updated mass matrices.
 *
 * @see initialize_metric_from_density() for the initial metric construction
 * @see update_edge_mass_matrix() for the edge mass matrix recomputation
 * @see fit_knn_riem_graph_regression() for the complete iteration context
 */
void riem_dcx_t::update_metric_from_density() {
    // Rebuild M₀ from updated vertex densities
    const size_t n_vertices = S[0].size();

    for (size_t i = 0; i < n_vertices; ++i) {
        double mass = std::max(rho.rho[0][i], 1e-15);
        g.M[0].coeffRef(i, i) = mass;
    }

    // Rebuild M₁ from updated vertex and edge densities
    update_edge_mass_matrix();

    // Invalidate spectral cache since Laplacian will change
    spectral_cache.invalidate();
}

/**
 * @brief Update edge mass matrix from evolved vertex densities
 *
 * Recomputes the edge mass matrix M₁ using current vertex densities after
 * density evolution. This function performs the same triple intersection
 * computations as compute_edge_mass_matrix() but operates in the context
 * of iterative refinement where densities have evolved from their initial
 * values through damped heat diffusion and response-coherence modulation.
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
 * This function is called by update_metric_from_density() as part of Step 4
 * in the iteration loop. At this point in each iteration:
 *   - Vertex densities ρ₀ have been evolved via damped heat diffusion
 *   - Edge densities ρ₁ have been recomputed from evolved ρ₀
 *   - Response-coherence modulation has adjusted ρ₁ based on fitted values
 *   - The combinatorial structure (S[1], neighbor_sets, stars[0]) remains fixed
 *
 * The recomputed M₁ will be used to reassemble the vertex Laplacian:
 *   L₀ = B₁ M₁⁻¹ B₁ᵀ M₀
 * which in turn drives the next iteration's response smoothing and density evolution.
 *
 * IMPLEMENTATION STRATEGY:
 * The current implementation performs full recomputation by delegating to
 * compute_edge_mass_matrix(), which iterates over all vertices and computes
 * inner products for all edge pairs in each vertex's star. This ensures
 * correctness and maintains consistency with the initialization logic.
 *
 * The function rebuilds M₁ completely by:
 *   1. Iterating over all vertices i
 *   2. For each vertex, examining all incident edge pairs (e_ij, e_is)
 *   3. Computing inner products via triple intersections using current ρ₀
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
 * @pre neighbor_sets contains kNN neighborhoods (unchanged from initialization)
 * @pre stars[0] contains vertex-to-incident-edges mapping (unchanged from initialization)
 * @pre g.M[1] is allocated with dimensions n_edges × n_edges
 *
 * @post g.M[1] contains updated edge mass matrix computed from current ρ₀
 * @post g.M[1] is symmetric: M₁[i,j] == M₁[j,i] for all i,j
 * @post g.M[1] is positive semidefinite with regularized diagonal entries
 * @post g.M_solver[1] is reset (factorization invalidated)
 *
 * @note This function modifies g.M[1] in place, replacing all entries with
 *       values computed from current vertex densities. Any previous matrix
 *       entries or factorizations are discarded.
 *
 * @note The combinatorial structure (which edges exist, vertex neighborhoods,
 *       star relationships) remains fixed throughout iteration. Only the
 *       numerical values in the mass matrix change based on evolved densities.
 *
 * @warning The current implementation performs full O(n·k²) recomputation at
 *          each call. For very large graphs, consider profiling and implementing
 *          incremental update strategies if this becomes a bottleneck.
 *
 * @see compute_edge_mass_matrix() for the initial mass matrix construction
 * @see update_metric_from_density() for the calling context
 * @see compute_edge_inner_product() for the triple intersection computation
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
        Rf_warning("apply_damped_heat_diffusion: CG solver did not converge (iterations=%d, error=%.3e)",
                   static_cast<int>(cg.iterations()),
                   cg.error());
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
 *   Step 3: Response-coherence modulation adjusts ρ₁ based on fitted values
 *   Step 4: Metric update rebuilds M₀ and M₁ from evolved densities
 *
 * The updated edge densities serve dual purposes:
 * 1. Input to response-coherence modulation (Step 3), which further adjusts
 *    edge densities based on response variation
 * 2. Diagonal entries of the edge mass matrix M₁ during metric update (Step 4),
 *    since ⟨e_ij, e_ij⟩ = ρ₁([i,j]) by construction
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
 * @note After this function returns, edge densities may be further modified by
 *       response-coherence modulation (Step 3) before being used in metric
 *       construction (Step 4).
 *
 * @see compute_initial_densities() for the initialization version
 * @see apply_response_coherence_modulation() for subsequent edge density adjustment
 * @see update_metric_from_density() for how updated densities enter the metric
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
 * @brief Apply response-coherence modulation to edge densities
 *
 * Adjusts edge densities based on response variation across each edge, creating
 * outcome-aware geometry where edges within response-coherent regions retain high
 * density (short distances in the Riemannian metric) while edges crossing response
 * boundaries receive low density (long distances, acting as diffusion barriers).
 *
 * This function implements Step 3 of the iteration cycle, following edge density
 * update and preceding metric reconstruction. The modulation mechanism is central
 * to the algorithm's ability to learn geometry that respects response structure:
 * regions where the response varies smoothly become well-connected through
 * high-mass edges, while response discontinuities create geometric barriers
 * that prevent diffusion across boundaries.
 *
 * MATHEMATICAL CONSTRUCTION:
 * For each edge [i,j], we measure response variation as the absolute difference
 * between endpoint fitted values:
 *   Δ_ij = |ŷ(i) - ŷ(j)|
 *
 * We compute an adaptive scale parameter σ₁ as the interquartile range (IQR) of
 * all edge variations {Δ_ij}. The IQR provides robust scale estimation that
 * resists outliers and adapts to the current response structure.
 *
 * Each edge density is then multiplied by a penalty function:
 *   Γ(Δ_ij) = (1 + Δ²_ij/σ²₁)^(-γ)
 *   ρ₁([i,j]) ← ρ₁([i,j]) · Γ(Δ_ij)
 *
 * The penalty function Γ takes values in (0,1]:
 *   - Γ(0) = 1: Edges with no response variation retain full density
 *   - Γ(Δ) → 0 as Δ → ∞: Edges with large variation lose density
 *   - γ controls decay rate: larger γ creates sharper boundaries
 *
 * After applying penalties, we renormalize edge densities to preserve total mass,
 * ensuring geometric scale stability across iterations.
 *
 * ITERATION CONTEXT:
 * Within each iteration of fit_knn_riem_graph_regression():
 *   Step 1: Density diffusion evolves vertex densities ρ₀
 *   Step 2: Edge densities ρ₁ recomputed from evolved ρ₀
 *   Step 3: Response-coherence modulation ← THIS FUNCTION
 *   Step 4: Metric update uses modulated ρ₁ to rebuild M₁
 *   Step 5: Laplacian reassembly with updated metric
 *   Step 6: Response smoothing with updated Laplacian
 *
 * The modulated edge densities directly influence the edge mass matrix M₁ diagonal
 * entries and thus the vertex Laplacian structure. High-penalty edges effectively
 * disconnect regions with different response values, allowing the response
 * smoothing (Step 6) to maintain sharp transitions while remaining smooth within
 * coherent regions.
 *
 * ADAPTIVE SCALE PARAMETER:
 * The IQR-based scale σ₁ = IQR({Δ_ij}) adapts to the current response distribution:
 *   - Early iterations: Large σ₁ if response variation is widespread
 *   - Late iterations: Smaller σ₁ as geometry focuses on remaining boundaries
 *   - Robustness: IQR ignores outliers, focusing on typical variation
 *
 * We use std::nth_element for efficient O(n) quartile computation without full
 * sorting. The first and third quartiles (Q1, Q3) define IQR = Q3 - Q1.
 *
 * DEGENERATE CASES:
 * Several edge cases require careful handling:
 *
 * 1. **Near-constant response** (σ₁ < 1e-10): If all edges have nearly identical
 *    response variation, the scale becomes ill-defined. We set σ₁ = 1e-10 to
 *    avoid division by zero, effectively disabling modulation when response
 *    variation is negligible.
 *
 * 2. **Uniform edge densities**: If all Δ_ij are equal, IQR = 0. The regularization
 *    σ₁ = max(σ₁, 1e-10) handles this case.
 *
 * 3. **Zero total mass after modulation**: If all edges receive strong penalties,
 *    the sum could approach zero. We check for this and issue a warning while
 *    proceeding with the degenerate normalization.
 *
 * GAMMA PARAMETER INTERPRETATION:
 * The exponent γ controls the strength of response-based geometric adaptation:
 *
 *   γ = 0.5: Gentle modulation, gradual boundaries
 *     Example: Δ = σ gives Γ = (2)^(-0.5) ≈ 0.71 (29% reduction)
 *
 *   γ = 1.0: Moderate modulation (recommended default)
 *     Example: Δ = σ gives Γ = (2)^(-1) = 0.50 (50% reduction)
 *
 *   γ = 2.0: Strong modulation, sharp boundaries
 *     Example: Δ = σ gives Γ = (2)^(-2) = 0.25 (75% reduction)
 *
 * For boundary-crossing edges with Δ = 3σ:
 *   γ = 0.5: Γ ≈ 0.32 (68% reduction)
 *   γ = 1.0: Γ = 0.10 (90% reduction)
 *   γ = 2.0: Γ = 0.01 (99% reduction)
 *
 * Typical applications use γ ∈ [0.5, 2], with γ = 1 providing good balance.
 *
 * NORMALIZATION STRATEGY:
 * After applying penalties, we normalize to preserve total edge mass:
 *   ρ₁ ← ρ₁ · (Σ ρ₁_old) / (Σ ρ₁_new)
 *
 * This ensures:
 *   - Total mass constant: Σ ρ₁ = n_edges before and after modulation
 *   - Only relative distribution changes, not overall scale
 *   - Geometric stability across iterations
 *   - Numerical conditioning of mass matrix
 *
 * The normalization compensates for the global mass loss from penalties, rescaling
 * the distribution while preserving the relative penalty effects.
 *
 * VERTEX DENSITIES UNCHANGED:
 * Critically, this function only modulates edge densities (dimension p=1). Vertex
 * densities ρ₀ remain unchanged, as they represent the underlying data distribution
 * and should not depend on response structure. This separation ensures:
 *   - Vertex masses reflect data density, not response variation
 *   - Response-aware geometry emerges only through edge connectivity
 *   - Interpretability: vertices = observations, edges = relationships
 *
 * COMPUTATIONAL COMPLEXITY:
 * O(n_edges log n_edges) for the IQR computation via nth_element, plus O(n_edges)
 * for computing variations and applying penalties. The dominant cost is the
 * nth_element calls for quartile computation.
 *
 * For very large graphs (n_edges > 10^6), consider approximating IQR by sampling
 * a random subset of edges, though this is rarely necessary in practice.
 *
 * NUMERICAL STABILITY:
 * The penalty function is numerically stable:
 *   - No cancellation: only additions and multiplications
 *   - Bounded output: Γ ∈ (0, 1] always
 *   - Smooth behavior: continuous and differentiable
 *   - Protected division: σ₁ regularized to avoid zero denominator
 *
 * The power function (1 + x)^(-γ) is well-conditioned for x ≥ 0 and moderate γ.
 *
 * @param y_hat Current fitted response values (from previous iteration's smoothing)
 * @param gamma Decay rate parameter controlling boundary sharpness (typically 0.5-2.0)
 *
 * @pre y_hat.size() == S[0].size() (fitted values for all vertices)
 * @pre rho.rho[1] contains current edge densities (from Step 2 of iteration)
 * @pre S[1] contains edge simplex table with vertex endpoints
 * @pre gamma > 0 (positive decay rate)
 *
 * @post rho.rho[1] contains modulated edge densities
 * @post rho.rho[1].sum() ≈ n_edges (preserved total mass)
 * @post All entries of rho.rho[1] are positive (penalties reduce but never zero)
 * @post rho.rho[0] is unchanged (vertex densities not modulated)
 *
 * @note This function modifies edge densities in place. The original edge density
 *       distribution is not preserved, as it's replaced by the modulated version.
 *
 * @note If gamma = 0, all edges receive Γ = 1 (no modulation). If gamma is very
 *       large (> 10), even small response variations create strong penalties,
 *       potentially fragmenting the geometry excessively.
 *
 * @warning Very small gamma (< 0.1) provides minimal geometric adaptation and may
 *          not effectively create response boundaries. Very large gamma (> 5) may
 *          over-fragment the geometry, preventing diffusion even within coherent
 *          regions.
 *
 * @see update_edge_densities_from_vertices() for Step 2 (edge density computation)
 * @see update_metric_from_density() for Step 4 (using modulated densities)
 * @see fit_knn_riem_graph_regression() for complete iteration context
 */
void riem_dcx_t::apply_response_coherence_modulation(
    const vec_t& y_hat,
    double gamma
) {
    const size_t n_edges = S[1].size();

    // Validate inputs
    if (static_cast<size_t>(y_hat.size()) != S[0].size()) {
        Rf_error("apply_response_coherence_modulation: y_hat size does not match number of vertices");
    }

    if (gamma <= 0.0) {
        Rf_error("apply_response_coherence_modulation: gamma must be positive");
    }

    // ============================================================
    // Step 1: Compute response variations for all edges
    // ============================================================

    std::vector<double> delta_values;
    delta_values.reserve(n_edges);

    for (size_t e = 0; e < n_edges; ++e) {
        // Get edge endpoints
        const std::vector<index_t>& edge_verts = S[1].simplex_verts[e];
        const index_t i = edge_verts[0];
        const index_t j = edge_verts[1];

        // Compute response variation: Δ_ij = |ŷ(i) - ŷ(j)|
        double delta = std::abs(y_hat[i] - y_hat[j]);
        delta_values.push_back(delta);
    }

    // ============================================================
    // Step 2: Compute adaptive scale parameter (IQR)
    // ============================================================

    // Compute first quartile (Q1) at 25th percentile
    std::vector<double> delta_copy1 = delta_values;
    size_t q1_index = delta_copy1.size() / 4;
    std::nth_element(delta_copy1.begin(),
                     delta_copy1.begin() + q1_index,
                     delta_copy1.end());
    double Q1 = delta_copy1[q1_index];

    // Compute third quartile (Q3) at 75th percentile
    std::vector<double> delta_copy2 = delta_values;
    size_t q3_index = (3 * delta_copy2.size()) / 4;
    std::nth_element(delta_copy2.begin(),
                     delta_copy2.begin() + q3_index,
                     delta_copy2.end());
    double Q3 = delta_copy2[q3_index];

    // Compute IQR = Q3 - Q1
    double sigma_1 = Q3 - Q1;

    // Regularize to avoid division by zero when response is nearly constant
    if (sigma_1 < 1e-10) {
        sigma_1 = 1e-10;

        // Note: Not issuing warning here as near-constant response is valid
        // and simply means no modulation is needed
    }

    // ============================================================
    // Step 3: Store total mass before modulation for normalization
    // ============================================================

    double total_mass_before = rho.rho[1].sum();

    // ============================================================
    // Step 4: Apply penalty function to each edge
    // ============================================================

    for (size_t e = 0; e < n_edges; ++e) {
        // Compute normalized variation
        double delta = delta_values[e];
        double delta_normalized = delta / sigma_1;

        // Compute penalty: Γ(Δ) = (1 + Δ²/σ²)^(-γ)
        //                       = (1 + (Δ/σ)²)^(-γ)
        double penalty = std::pow(1.0 + delta_normalized * delta_normalized, -gamma);

        // Apply penalty to edge density
        rho.rho[1][e] *= penalty;
    }

    // ============================================================
    // Step 5: Renormalize to preserve total edge mass
    // ============================================================

    double total_mass_after = rho.rho[1].sum();

    if (total_mass_after > 1e-15) {
        // Normal case: rescale to restore original total mass
        rho.rho[1] *= total_mass_before / total_mass_after;
    } else {
        // Degenerate case: all edges received extreme penalties
        // Fall back to uniform density and warn
        rho.rho[1].setConstant(1.0);

        Rf_warning("apply_response_coherence_modulation: edge density sum near zero "
                   "after modulation, falling back to uniform density. "
                   "Consider reducing gamma parameter.");
    }

    // Note: Vertex densities rho.rho[0] are intentionally NOT modified
    // They represent data distribution and should not depend on response
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
    const int n_vertices = y.size();
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
 * field and the density distributions across simplex dimensions. The
 * iteration is considered converged when both quantities have stabilized
 * below specified thresholds.
 *
 * The convergence criteria are based on relative changes to ensure scale
 * independence. We measure the L2 norm of differences normalized by the
 * L2 norm of the previous state. This provides a dimensionless measure
 * of change that is robust to the scale of the problem.
 *
 * For density monitoring, we track changes across all dimensions where
 * density is defined. In the standard implementation, this includes vertex
 * densities (dimension 0) and edge densities (dimension 1). The maximum
 * relative change across all dimensions serves as the geometric convergence
 * criterion.
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
 *         - max_density_change: maximum relative change across all densities
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
    double threshold_percentile
    ) {

#define DEBUG_FIT_KNN_RIEM_GRAPH_REGRESSION 0

    // ================================================================
    // PART I: INITIALIZATION
    // ================================================================

    // ----------------------------------------------------------------
    // Phase 1a: Build 1-skeleton geometry
    // ----------------------------------------------------------------

#if DEBUG_FIT_KNN_RIEM_GRAPH_REGRESSION
    Rprintf("Phase 1: Computing k-NN neighborhoods...\n");
#endif

    const size_t n_points = static_cast<size_t>(X.rows());

    // Step 1: Compute k-NN using Eigen-compatible wrapper
    knn_result_t knn_result = compute_knn_from_eigen(X, k);

    // Extract results into convenient format
    std::vector<std::vector<index_t>> knn_indices(n_points, std::vector<index_t>(k));
    std::vector<std::vector<double>> knn_distances(n_points, std::vector<double>(k));
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < k; ++j) {
            const size_t offset = i * k + j;
            knn_indices[i][j] = static_cast<index_t>(knn_result.indices[offset]);
            knn_distances[i][j] = knn_result.distances[offset];
        }
    }

    // Step 2: Build neighborhood sets
    std::vector<std::unordered_set<index_t>> neighbor_sets(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

    // Step 3: Build edges via O(nk) neighborhood intersection algorithm
    std::vector<std::vector<iknn_vertex_t>> adjacency_list(n_points);
    std::vector<index_t> intersection;
    intersection.reserve(k);  // Pre-allocate for efficiency

    for (size_t i = 0; i < n_points; ++i) {
        const std::vector<index_t>& neighbors_i = knn_indices[i];

        // Only examine neighbors j where j > i to avoid duplicate edges
        for (index_t j : neighbor_sets[i]) {
            if (j <= i) continue;  // Skip self and already-processed pairs

            // Compute intersection of neighborhoods
            intersection.clear();
            for (index_t v : neighbor_sets[i]) {
                if (neighbor_sets[j].find(v) != neighbor_sets[j].end()) {
                    intersection.push_back(v);
                }
            }

            size_t common_count = intersection.size();

            if (common_count > 0) {
                // Compute minimum distance through common neighbors
                double min_dist = std::numeric_limits<double>::max();

                for (index_t x_k : intersection) {
                    // Find x_k in i's neighbor list to get distance
                    auto it_i = std::find(neighbors_i.begin(), neighbors_i.end(), x_k);
                    size_t idx_i = it_i - neighbors_i.begin();
                    double dist_i_k = knn_distances[i][idx_i];

                    // Find x_k in j's neighbor list to get distance
                    const std::vector<index_t>& neighbors_j = knn_indices[j];
                    auto it_j = std::find(neighbors_j.begin(), neighbors_j.end(), x_k);
                    size_t idx_j = it_j - neighbors_j.begin();
                    double dist_j_k = knn_distances[j][idx_j];

                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                // Add edges (bidirectional)
                adjacency_list[i].emplace_back(iknn_vertex_t{j, common_count, min_dist});
                adjacency_list[j].emplace_back(iknn_vertex_t{i, common_count, min_dist});
            }
        }
    }


    // Step 4: Store neighbor sets as member variable for later use
    this->neighbor_sets = std::move(neighbor_sets);

    // ----------------------------------------------------------------
    // Phase 1b: Edge construction with geometric edge pruning
    // ----------------------------------------------------------------
    // Extract edges from adjacency list
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

    // Apply geometric pruning
    set_wgraph_t temp_graph(n_points);
    for (const auto& edge : all_edges) {
        temp_graph.add_edge(edge.i, edge.j, edge.dist);
    }

    set_wgraph_t pruned_graph = temp_graph.prune_edges_geometrically(
        max_ratio_threshold, threshold_percentile);

    // Rebuild edge list with only retained edges
    std::vector<temp_edge> pruned_edges;
    for (const auto& edge : all_edges) {
        // Check if edge still exists in pruned graph
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

    // Now build simplex structures from all_edges
    const index_t n_edges = all_edges.size();
    extend_by_one_dim(n_edges);

    S[1].simplex_verts.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        std::vector<index_t> verts = {all_edges[e].i, all_edges[e].j};
        S[1].simplex_verts[e] = verts;
        S[1].id_of[verts] = e;
    }

    // ----------------------------------------------------------------
    // Phase 1c: Build 0-Simplices (Vertices)
    // ----------------------------------------------------------------
    S[0].simplex_verts.resize(n_points);
    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        S[0].simplex_verts[i] = {i};
        S[0].id_of[{i}] = i;
    }

    // ----------------------------------------------------------------
    // Phase 1d: Initialize reference measure
    // ----------------------------------------------------------------
    // Initilizting reference_measure
    initialize_reference_measure(
        knn_indices,
        knn_distances,
        use_counting_measure,
        density_normalization
        );

    // ----------------------------------------------------------------
    // Phase 1e: Initialize vertex masses: m_i = μ(N̂_k(x_i))
    // ----------------------------------------------------------------
    vec_t vertex_masses = vec_t::Zero(n_points);
    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        double mass = 0.0;
        for (index_t j : neighbor_sets[i]) {
            mass += reference_measure[j];
        }
        vertex_masses[i] = mass;
    }

    g.M[0] = spmat_t(n_points, n_points);
    g.M[0].reserve(Eigen::VectorXi::Constant(n_points, 1));
    for (size_t i = 0; i < n_points; ++i) {
        g.M[0].insert(i, i) = std::max(vertex_masses[i], 1e-15);
    }
    g.M[0].makeCompressed();


#if DEBUG_FIT_KNN_RIEM_GRAPH_REGRESSION
    // Count total edges
    size_t total_edges = 0;
    for (const auto& adj : adjacency_list) {
        total_edges += adj.size();
    }
    Rprintf("Phase 1 complete: %zu vertices, %zu directed edges (%zu undirected)\n",
            n_points, total_edges, total_edges / 2);
#endif

    // ----------------------------------------------------------------
    // Phase 2: Compute initial densities
    // ----------------------------------------------------------------
    compute_initial_densities();

    // ----------------------------------------------------------------
    // Phase 3: Build initial metric
    // ----------------------------------------------------------------
    initialize_metric_from_density();

    // ----------------------------------------------------------------
    // Phase 4: Assemble initial Laplacian
    // ----------------------------------------------------------------
    assemble_operators();

    // ----------------------------------------------------------------
    // Phase 4.5: Select diffusion parameters (auto or validate user input)
    // ----------------------------------------------------------------
    select_diffusion_parameters(t_diffusion, beta_damping, /*verbose=*/true);

    // ----------------------------------------------------------------
    // Phase 5: Initial response smoothing
    // ----------------------------------------------------------------
    auto gcv_result = this->smooth_response_via_spectral_filter(
        y, n_eigenpairs, filter_type
    );
    vec_t y_hat_curr = gcv_result.y_hat;
    sig.y_hat_hist.clear();
    sig.y_hat_hist.push_back(y_hat_curr);

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
        rho.rho[0] = apply_damped_heat_diffusion(rho.rho[0], t_diffusion, beta_damping);

        // Step 2: Edge density update
        update_edge_densities_from_vertices();

        // Step 3: Response-coherence modulation
        apply_response_coherence_modulation(y_hat_curr, gamma_modulation);

        // Step 4: Metric update
        update_metric_from_density();

        // Step 5: Laplacian reassembly
        assemble_operators();

        // Invalidate spectral cache since L[0] changed
        spectral_cache.invalidate();

        // Step 6: Response smoothing
        gcv_result = smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
        );
        y_hat_curr = gcv_result.y_hat;
        sig.y_hat_hist.push_back(y_hat_curr);

        // Step 7: Convergence check
        const bool research_mode = true;
        if (research_mode) {
            auto status = check_convergence_detailed(
                y_hat_prev, y_hat_curr,
                rho_prev, rho.rho,
                epsilon_y, epsilon_rho,
                iter, max_iterations,
                response_change_history
                );

            // Store history for next iteration
            response_change_history.push_back(status.response_change);

            // Log detailed diagnostics
            Rprintf("Iteration %d: response_change=%.6f, density_change=%.6f\n",
                    iter, status.response_change, status.max_density_change);

            // Optional: print detailed per-dimension info
            if (1) {
                for (size_t p = 0; p < status.density_changes_by_dim.size(); ++p) {
                    Rprintf("  Density[%d] change: %.6f\n",
                            (int)p, status.density_changes_by_dim[p]);
                }
                if (status.estimated_iterations_remaining > 0) {
                    Rprintf("  Estimated iterations remaining: %d\n",
                            status.estimated_iterations_remaining);
                }
            }

            if (status.converged) {
                Rprintf("%s\n", status.message.c_str());
                break;
            }
        } else {
            auto status = check_convergence(
                y_hat_prev, y_hat_curr,
                rho_prev, rho.rho,
                epsilon_y, epsilon_rho,
                iter, max_iterations
                );

            Rprintf("Iteration %d: response_change=%.6f, density_change=%.6f\n",
                    iter, status.response_change, status.max_density_change);

            if (status.converged) {
                Rprintf("%s\n", status.message.c_str());
                break;
            }
        }
    }

    // ================================================================
    // PART III: FINALIZATION
    // ================================================================

    // Store original response
    sig.y = y;

    // Final state already in place
}
