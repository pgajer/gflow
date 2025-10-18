/**
 * @file riem_dcx_coboundary_full.cpp
 * @brief Complete implementation of coboundary operator with non-diagonal metric
 *
 * This file provides the full implementation of the coboundary operator ∂₁*
 * for non-diagonal Riemannian metrics, where edges incident to the same vertex
 * have non-orthogonal inner products encoded in local Gram matrices.
 *
 */

#include "riem_dcx.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <limits>

/**
 * @brief Compute coboundary operator ∂₁* with full non-diagonal metric
 *
 * Computes ∂₁*(f) for a vertex function f using the complete Riemannian
 * structure where edges incident to the same vertex may have non-orthogonal
 * inner products encoded in the edge mass matrix M[1].
 *
 * MATHEMATICAL FORMULATION:
 * Solves M₁ x = B₁ᵀ M₀ f where:
 * - M₁ is the full edge mass matrix (potentially non-diagonal)
 * - B₁ is the boundary operator
 * - M₀ is the vertex mass matrix (diagonal)
 * - x gives the coboundary ∂₁*(f)
 *
 * COMPUTATIONAL METHODS:
 * - If M₁ is diagonal: Direct solve O(m), where m is the number of edges
 * - If M₁ prefactored: Cholesky solve O(m)
 * - If use_iterative: Conjugate gradient O(k·nnz(M₁))
 * - Otherwise: Factor then solve O(m^1.5) first call, O(m) subsequent
 *
 * This method may cache the Cholesky factorization of M₁ for efficiency,
 * making subsequent calls much faster.
 *
 * @param f Vertex function (0-chain) with length n_vertices
 * @param use_iterative Force iterative CG solve even if factorization available
 * @param cg_tol Conjugate gradient convergence tolerance
 * @param cg_maxit Maximum CG iterations
 *
 * @return Edge function (1-chain) with length n_edges giving ∂₁*(f)
 *
 * @complexity O(m^1.5) first call with factorization, O(m) subsequent calls
 *
 * @post May populate g.M_solver[1] with Cholesky factorization
 *
 * @see compute_edge_mass_matrix() constructs the non-diagonal M₁
 * @see compute_all_extremality_scores_full() uses this for extrema detection
 */
vec_t riem_dcx_t::compute_coboundary_del1star_full(
    const vec_t& f,
    bool use_iterative,
    double cg_tol,
    int cg_maxit
    ) const {

    // ========================================================================
    // VALIDATION AND DIMENSION EXTRACTION
    // ========================================================================

    const size_t n_vertices = vertex_cofaces.size();
    const size_t n_edges = edge_registry.size();
    const Eigen::Index n = static_cast<Eigen::Index>(n_vertices);
    const Eigen::Index m = static_cast<Eigen::Index>(n_edges);

    // Validate input dimension
    if (f.size() != n) {
        Rf_error("compute_coboundary_del1star_full: f.size()=%lld must equal n_vertices=%zu",
                 static_cast<long long>(f.size()), n_vertices);
    }

    // Validate metric structures
    if (g.M.size() < 2) {
        Rf_error("compute_coboundary_del1star_full: metric not initialized (g.M.size()=%zu < 2)",
                 g.M.size());
    }

    if (g.M[0].rows() != n) {
        Rf_error("compute_coboundary_del1star_full: vertex metric dimension mismatch");
    }

    if (g.M[1].rows() != m) {
        Rf_error("compute_coboundary_del1star_full: edge metric dimension mismatch");
    }

    // Validate boundary operator exists
    if (L.B.size() < 2 || L.B[1].rows() != n || L.B[1].cols() != m) {
        Rf_error("compute_coboundary_del1star_full: boundary operator B[1] not properly initialized");
    }

    // ========================================================================
    // COMPUTE RIGHT-HAND SIDE: b = B₁ᵀ M₀ f
    // ========================================================================

    // For diagonal M₀, we can compute B₁ᵀ M₀ f efficiently.
    // Since B₁(i, e) = -1 for tail, B₁(j, e) = +1 for head of edge e = [i,j]:
    //   (B₁ᵀ M₀ f)[e] = w_j f_j - w_i f_i

    vec_t b = vec_t::Zero(m);

    for (Eigen::Index e = 0; e < m; ++e) {
        const auto [i, j] = edge_registry[e];
        const Eigen::Index ei = static_cast<Eigen::Index>(i);
        const Eigen::Index ej = static_cast<Eigen::Index>(j);

        // Extract vertex weights (diagonal of M₀)
        const double w_i = std::max(g.M[0].coeff(ei, ei), 1e-15);
        const double w_j = std::max(g.M[0].coeff(ej, ej), 1e-15);

        // Compute weighted boundary value
        b[e] = w_j * f[ej] - w_i * f[ei];
    }

    // ========================================================================
    // SOLVE M₁ x = b FOR DIFFERENT METRIC STRUCTURES
    // ========================================================================

    vec_t x;  // Solution vector

    // ------------------------------------------------------------------------
    // CASE 1: DIAGONAL METRIC (Fast path)
    // ------------------------------------------------------------------------

    if (g.is_diagonal(1)) {
        // Diagonal solve: x[e] = b[e] / M₁(e,e)
        x = vec_t(m);
        for (Eigen::Index e = 0; e < m; ++e) {
            const double diag = std::max(g.M[1].coeff(e, e), 1e-15);
            x[e] = b[e] / diag;
        }
        return x;
    }

    // ------------------------------------------------------------------------
    // CASE 2: NON-DIAGONAL METRIC WITH ITERATIVE SOLVE
    // ------------------------------------------------------------------------

    if (use_iterative) {
        // Use conjugate gradient with diagonal preconditioning

        // Extract diagonal for preconditioning
        vec_t diag_M1 = g.get_diagonal(1);
        for (Eigen::Index i = 0; i < diag_M1.size(); ++i) {
            diag_M1[i] = std::max(diag_M1[i], 1e-15);
        }

        // Setup conjugate gradient solver
        Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;
        cg.setTolerance(cg_tol);
        cg.setMaxIterations(cg_maxit);

        cg.compute(g.M[1]);

        if (cg.info() != Eigen::Success) {
            Rf_error("compute_coboundary_del1star_full: CG setup failed");
        }

        // Solve
        x = cg.solve(b);

        if (cg.info() != Eigen::Success) {
            Rf_warning("compute_coboundary_del1star_full: CG did not fully converge "
                      "(iterations=%d, error=%.3e)",
                      static_cast<int>(cg.iterations()), cg.error());
            // Continue with best available solution
        }

        return x;
    }

    // ------------------------------------------------------------------------
    // CASE 3: NON-DIAGONAL METRIC WITH DIRECT SOLVE
    // ------------------------------------------------------------------------

    // Check if factorization already exists
    if (g.M_solver[1] && g.M_solver[1]->info() == Eigen::Success) {
        // Use existing factorization
        x = g.M_solver[1]->solve(b);

        if (g.M_solver[1]->info() != Eigen::Success) {
            Rf_error("compute_coboundary_del1star_full: Cholesky solve failed");
        }

        return x;
    }

    // Need to factor M₁ first
    // Note: This modifies g.M_solver[1], but method is marked const
    // We cast away const for factorization caching (doesn't change mathematical state)

    auto* mutable_g = const_cast<metric_family_t*>(&g);

    auto solver = std::make_unique<Eigen::SimplicialLLT<spmat_t>>();
    solver->compute(g.M[1]);

    if (solver->info() != Eigen::Success) {
        const char* reason;
        switch (solver->info()) {
        case Eigen::NumericalIssue:
            reason = "numerical issue (matrix may not be positive definite)";
            break;
        case Eigen::NoConvergence:
            reason = "no convergence";
            break;
        case Eigen::InvalidInput:
            reason = "invalid input";
            break;
        default:
            reason = "unknown error";
            break;
        }

        Rf_error("compute_coboundary_del1star_full: M[1] factorization failed: %s "
                 "(size=%ld, nnz=%ld)",
                 reason, static_cast<long>(g.M[1].rows()),
                 static_cast<long>(g.M[1].nonZeros()));
    }

    // Solve with new factorization
    x = solver->solve(b);

    if (solver->info() != Eigen::Success) {
        Rf_error("compute_coboundary_del1star_full: Cholesky solve failed after factorization");
    }

    // Cache the factorization for future use
    mutable_g->M_solver[1] = std::move(solver);

    return x;
}

/**
 * @brief Compute extremality scores using full non-diagonal metric
 *
 * Computes extr(v; f, g) for all vertices using the complete Riemannian
 * structure including non-orthogonal edge inner products.
 *
 * The extremality score measures directional coherence of the discrete
 * gradient field ∂₁*(f), computed using the full non-diagonal metric:
 *
 *   extr(v; f, g) = -∑_{e∋v} ρ₁(e)·∂₁*(f)(e) / ∑_{e∋v} ρ₁(e)·|∂₁*(f)(e)|
 *
 * where ∂₁*(f) is obtained by solving M₁ x = B₁ᵀ M₀ f with non-diagonal M₁.
 *
 * For weakly correlated edges, both versions give similar results. For
 * strongly correlated edges (dense triple intersections, high-degree vertices),
 * the full version can produce significantly different and more accurate scores.
 *
 * @param f Vertex function (0-chain)
 * @param use_iterative Use iterative solver for coboundary (default: false)
 * @param cg_tol CG tolerance (default: 1e-10)
 * @param cg_maxit CG max iterations (default: 1000)
 *
 * @return Vector of extremality scores in [-1, 1] or NaN for isolated vertices
 *
 * @complexity O(m^1.5 + n) first call, O(m + n) subsequent calls
 *             (factorization cost amortized over multiple calls)
 *
 * @note This method internally calls compute_coboundary_del1star_full()
 *       which may cache the factorization of M₁ for efficiency.
 *
 * @see compute_coboundary_del1star_full() for coboundary computation details
 * @see compute_all_extremality_scores() for faster diagonal approximation
 */
vec_t riem_dcx_t::compute_all_extremality_scores_full(
    const vec_t& f,
    bool use_iterative,
    double cg_tol,
    int cg_maxit
    ) const {

    const size_t n_vertices = vertex_cofaces.size();
    const Eigen::Index n = static_cast<Eigen::Index>(n_vertices);

    // ========================================================================
    // PHASE 1: COMPUTE COBOUNDARY WITH FULL METRIC
    // ========================================================================

    const vec_t del1star_f = compute_coboundary_del1star_full(
        f, use_iterative, cg_tol, cg_maxit
    );

    // ========================================================================
    // PHASE 2: AGGREGATE AT EACH VERTEX
    // ========================================================================

    vec_t extremality = vec_t::Constant(n, std::numeric_limits<double>::quiet_NaN());

    for (size_t v = 0; v < n_vertices; ++v) {
        // Skip self-loop at position 0
        if (vertex_cofaces[v].size() <= 1) {
            // Isolated vertex (no neighbors)
            continue;
        }

        double numerator = 0.0;
        double denominator = 0.0;

        // Iterate through incident edges
        for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
            const auto& neighbor_info = vertex_cofaces[v][k];

            const index_t edge_idx = neighbor_info.simplex_index;

            // Edge density is the diagonal entry of M₁
            // (stored in neighbor_info.density for convenience)
            const double rho_edge = std::max(neighbor_info.density, 1e-15);

            // Get coboundary value for this edge
            double del1star_val = del1star_f[static_cast<Eigen::Index>(edge_idx)];

            // ORIENTATION CORRECTION
            // Edge registry stores [i,j] with i < j
            // ∂₁*(f)[e] is oriented from i → j
            // If v is the head (j), reverse sign for v's perspective
            const auto [edge_tail, edge_head] = edge_registry[edge_idx];

            if (v == edge_head) {
                del1star_val = -del1star_val;
            }

            // Accumulate weighted contributions
            numerator += rho_edge * del1star_val;
            denominator += rho_edge * std::abs(del1star_val);
        }

        // Compute extremality score with sign reversal
        // (negative because inward gradients indicate maximum)
        if (denominator > 1e-15) {
            extremality[static_cast<Eigen::Index>(v)] = -numerator / denominator;
        }
        // else: leave as NaN (all incident gradients are zero)
    }

    return extremality;
}
