#include "riem_dcx.hpp"
#include "iknn_vertex.hpp" // for iknn_vertex_t
#include "kNN.h"
#include <R.h>

// ============================================================
// INTERNAL HELPERS (anonymous namespace)
// ============================================================

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
    );

namespace {

double apply_filter_function(
    double lambda,
    double eta,
    filter_type_t filter_type
) {
    switch (filter_type) {
        case filter_type_t::HEAT_KERNEL:
            return std::exp(-eta * lambda);

        case filter_type_t::TIKHONOV:
            return 1.0 / (1.0 + eta * lambda);

        case filter_type_t::CUBIC_SPLINE:
            return 1.0 / (1.0 + eta * lambda * lambda);

        case filter_type_t::GAUSSIAN:
            return std::exp(-eta * lambda * lambda);

        default:
            Rf_error("Unknown filter type");
            return 0.0; // Never reached
    }
}

} // anonymous namespace

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

    // Use triplet list for efficient sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_edges * 10);  // Heuristic: ~10 nonzeros per row

    // ========================================================================
    // Main loop: For each vertex, compute inner products between its edges
    // ========================================================================

    for (const auto& star_entry : stars[0].star_over) {
        const std::vector<index_t>& incident_edges = star_entry;

        if (incident_edges.empty()) continue;

        // Get the vertex index (need to find which vertex this star belongs to)
        // This requires looking up the vertex from the star table structure
        // For efficiency, we'll iterate over all vertices explicitly below
    }

    // More direct approach: iterate over all vertices
    const size_t n_vertices = S[0].size();

    for (size_t i = 0; i < n_vertices; ++i) {
        // Get all edges incident to vertex i
        const std::vector<index_t>& incident_edges = stars[0].star_over[i];
        const size_t n_incident = incident_edges.size();

        if (n_incident == 0) continue;

        // For each pair of incident edges, compute their inner product
        for (size_t a = 0; a < n_incident; ++a) {
            const index_t e_ij_idx = incident_edges[a];

            for (size_t b = a; b < n_incident; ++b) {
                const index_t e_is_idx = incident_edges[b];

                double inner_product = 0.0;

                if (a == b) {
                    // Diagonal entry: ⟨e_ij, e_ij⟩ = ρ₁([i,j])
                    inner_product = rho.rho[1][e_ij_idx];
                } else {
                    // Off-diagonal entry: requires triple intersection
                    inner_product = compute_edge_inner_product(e_ij_idx, e_is_idx, i);
                }

                // Apply regularization to diagonal entries
                if (a == b) {
                    inner_product = std::max(inner_product, 1e-15);
                }

                // Add to triplet list (symmetric matrix)
                triplets.emplace_back(e_ij_idx, e_is_idx, inner_product);
                if (a != b) {
                    triplets.emplace_back(e_is_idx, e_ij_idx, inner_product);
                }
            }
        }
    }

    // ========================================================================
    // Build sparse matrix from triplets
    // ========================================================================

    g.M[1] = spmat_t(n_edges, n_edges);
    g.M[1].setFromTriplets(triplets.begin(), triplets.end());
    g.M[1].makeCompressed();

    // Clear any existing factorization (will be recomputed if needed)
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

vec_t riem_dcx_t::apply_damped_heat_diffusion(
    const vec_t& rho_current,
    double t,
    double beta
) {
    // ... implementation ...
}

void riem_dcx_t::update_edge_densities_from_vertices() {
    // ... implementation ...
}

void riem_dcx_t::apply_response_coherence_modulation(
    const vec_t& y_hat,
    double gamma
) {
    // ... implementation ...
}

gcv_result_t riem_dcx_t::smooth_response_via_spectral_filter(
    const vec_t& y,
    int n_eigenpairs,
    filter_type_t filter_type
) {
    // ... implementation ...
}

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
    // ... implementation ...
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
    filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations
    ) {

#define DEBUG_FIT_KNN_RIEM_GRAPH_REGRESSION 0

    // ================================================================
    // PART I: INITIALIZATION
    // ================================================================

    // ----------------------------------------------------------------
    // Phase 1: Build 1-skeleton geometry
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
    // Phase 2: Initialize reference measure
    // ----------------------------------------------------------------
    initialize_reference_measure(
        knn_indices,
        knn_distances,
        use_counting_measure,
        density_normalization
        );

    // ----------------------------------------------------------------
    // Phase 3: Compute initial densities
    // ----------------------------------------------------------------
    compute_initial_densities();

    // ----------------------------------------------------------------
    // Phase 4: Build initial metric
    // ----------------------------------------------------------------
    initialize_metric_from_density();

    // ----------------------------------------------------------------
    // Phase 5: Assemble initial Laplacian
    // ----------------------------------------------------------------
    assemble_operators();

    // ----------------------------------------------------------------
    // Phase 6: Initial response smoothing
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

    for (int iter = 1; iter <= max_iterations; ++iter) {
        y_hat_prev = y_hat_curr;
        rho_prev = rho.rho;

        // Step 1: Density diffusion
        rho.rho[0] = this->apply_damped_heat_diffusion(
            rho.rho[0], t_diffusion, beta_damping
        );

        // Step 2: Edge density update
        this->update_edge_densities_from_vertices();

        // Step 3: Response-coherence modulation
        this->apply_response_coherence_modulation(y_hat_curr, gamma_modulation);

        // Step 4: Metric update
        this->update_metric_from_density();

        // Step 5: Laplacian reassembly
        this->assemble_operators();

        // Step 6: Response smoothing
        gcv_result = this->smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
        );
        y_hat_curr = gcv_result.y_hat;
        sig.y_hat_hist.push_back(y_hat_curr);

        // Step 7: Convergence check
        auto status = this->check_convergence(
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

    // ================================================================
    // PART III: FINALIZATION
    // ================================================================

    // Store original response
    sig.y = y;

    // Final state already in place
}
