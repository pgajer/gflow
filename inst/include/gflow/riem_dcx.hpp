#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include <unordered_set>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <R.h>

// Type aliases
using index_t = size_t;
using vec_t = Eigen::VectorXd;
using spmat_t = Eigen::SparseMatrix<double>;

// ================================================================
// SUPPORTING ENUMERATIONS AND STRUCTURES
// ================================================================

/**
 * @brief Filter type enumeration for spectral filtering
 */
enum class rdcx_filter_type_t {
    HEAT_KERNEL,   ///< f(λ) = exp(-ηλ), exponential decay
    TIKHONOV,      ///< f(λ) = 1/(1 + ηλ), rational decay
    CUBIC_SPLINE,  ///< f(λ) = 1/(1 + ηλ²), spline smoothness
    GAUSSIAN       ///< f(λ) = exp(-ηλ^2)
};

/**
 * @brief Result structure for GCV-based spectral filtering
 */
struct gcv_result_t {
    double eta_optimal;              ///< Optimal smoothing parameter
    vec_t y_hat;                     ///< Smoothed response vector
    std::vector<double> gcv_scores;  ///< GCV scores for each eta candidate
    std::vector<double> eta_grid;    ///< Grid of eta values evaluated
};

/**
 * @brief Convergence status for iterative refinement
 */
struct convergence_status_t {
    bool converged;              ///< True if convergence criteria met
    double response_change;      ///< Relative change in fitted values
    double max_density_change;   ///< Maximum relative change in densities
    int iteration;               ///< Current iteration number
    std::string message;         ///< Human-readable convergence message
};

/**
 * @brief Complete result structure for regression fitting
 */
struct regression_result_t {
    vec_t y_hat_final;                      ///< Final fitted values
    int n_iterations;                       ///< Number of iterations performed
    bool converged;                         ///< Convergence status
    std::vector<double> response_changes;   ///< History of response changes
    std::vector<double> density_changes;    ///< History of density changes
    std::vector<vec_t> y_hat_history;       ///< Full fitted values history
};

// ================================================================
// COMPONENT STRUCTURES (definitions follow)
// ================================================================

// ========================= Utility: hash for simplex keys =========================
struct vec_hash_t {
    std::size_t operator()(const std::vector<index_t>& v) const noexcept {
        std::size_t h = 1469598103934665603ull; // FNV-like
        for (auto x : v) {
            h ^= x + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        }
        return h;
    }
};

/**
 * @brief Lookup table for simplices of fixed dimension
 *
 * Maintains a bidirectional mapping between simplices (represented as sorted
 * vertex lists) and their integer identifiers. All simplices in a table must
 * have the same dimension (number of vertices). Used in discrete Morse theory
 * computations to efficiently reference and retrieve simplices.
 *
 * The vertex indices in each simplex are maintained in sorted order to ensure
 * consistent lookup and comparison.
 */
struct simplex_table_t {
    /// Storage for simplex vertex lists, indexed by simplex ID
    std::vector<std::vector<index_t>> simplex_verts;

    /// Map from sorted vertex list to simplex ID for fast lookup
    std::unordered_map<std::vector<index_t>, index_t, vec_hash_t> id_of;

    /**
     * @brief Get the number of simplices in the table
     * @return The count of stored simplices
     */
    index_t size() const { return static_cast<index_t>(simplex_verts.size()); }

    /**
     * @brief Add a simplex to the table or retrieve its existing ID
     *
     * @param v Vector of vertex indices defining the simplex
     * @return The simplex ID (newly assigned or existing)
     *
     * If the simplex already exists, returns its ID. Otherwise, assigns a new
     * ID and stores the simplex. The input vertices are sorted internally.
     * All simplices must have the same dimension; adding a simplex with a
     * different number of vertices than existing simplices triggers an error.
     */
    index_t add_simplex(std::vector<index_t> v) {
        std::sort(v.begin(), v.end());

        // Validate dimension consistency
        if (!simplex_verts.empty()) {
            size_t expected_dim = simplex_verts[0].size();
            if (v.size() != expected_dim) {
                Rf_error("Dimension mismatch: expected %zu vertices, got %zu",
                         expected_dim, v.size());
            }
        }

        auto it = id_of.find(v);
        if (it != id_of.end()) return it->second;

        index_t id = size();
        simplex_verts.push_back(v);
        id_of.emplace(std::move(v), id);
        return id;
    }

    /**
     * @brief Retrieve the ID of an existing simplex
     *
     * @param v Vector of vertex indices (will be sorted for lookup)
     * @return The simplex ID
     *
     * @note Triggers an R error if the simplex is not found in the table
     */
    index_t get_id(const std::vector<index_t>& v) const {
        auto key = v;
        std::sort(key.begin(), key.end());
        auto it = id_of.find(key);
        if (it == id_of.end()) Rf_error("Simplex not found");
        return it->second;
    }
};

/**
 * @brief Star of simplices in the complex
 *
 * For each p-simplex, maintains the list of (p+1)-simplices in its star.
 * The star of a simplex consists of all higher-dimensional simplices that
 * contain it as a face. Used in discrete Morse theory to traverse the
 * complex upward through dimensions.
 */
struct star_table_t {
    /// For each p-simplex (indexed by ID), stores IDs of (p+1)-simplices in its star
    std::vector<std::vector<index_t>> star_over;

    /**
     * @brief Initialize storage for a given number of p-simplices
     * @param n_from Number of p-simplices to allocate storage for
     */
    void resize(index_t n_from) { star_over.assign(n_from, {}); }
};

/**
 * @brief Riemannian metric on chain spaces
 *
 * Encodes the geometry of the simplicial complex through mass (metric) matrices
 * at each dimension. The metric determines inner products between chains and
 * influences all geometric computations including Laplacians and heat flows.
 * Supports both diagonal and general sparse symmetric positive definite matrices.
 */
struct metric_family_t {
    /// Sparse symmetric positive definite mass matrices for each dimension
    std::vector<spmat_t> M;

    /// Cholesky factorizations of M for efficient solving (null if diagonal)
    std::vector<std::unique_ptr<Eigen::SimplicialLLT<spmat_t>>> M_solver;

    void init_dims(const std::vector<index_t>& n_by_dim) {
        const int P = static_cast<int>(n_by_dim.size()) - 1;
        M.resize(P+1); M_solver.resize(P+1);
        for (int p=0; p<=P; ++p) {
            const Eigen::Index n = static_cast<Eigen::Index>(n_by_dim[p]);
            M[p] = spmat_t(n, n);
            M[p].setIdentity();
            M_solver[p].reset();
        }
    }

    bool is_diagonal(int p) const {
        return (M[p].nonZeros() == M[p].rows());
    }

    /**
     * @brief Extract diagonal entries of M[p]
     * @param p Dimension index
     * @return Vector of diagonal values
     */
    vec_t get_diagonal(int p) const {
        const Eigen::Index n = M[p].rows();
        vec_t d(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            d[i] = M[p].coeff(i, i);
        }
        return d;
    }

    /**
     * @brief Set M[p] as a diagonal matrix from vector
     * @param p Dimension index
     * @param diag Diagonal values
     */
    void set_diagonal(int p, const vec_t& diag) {
        const Eigen::Index n = M[p].rows();
        M[p].setZero();
        M[p].reserve(Eigen::VectorXi::Constant(n, 1));
        for (Eigen::Index i = 0; i < n; ++i) {
            M[p].insert(i, i) = std::max(diag[i], 1e-15);
        }
        M[p].makeCompressed();
        M_solver[p].reset();
    }

    void normalize() {
        if (!M.empty() && M[0].rows() > 0) {
            vec_t m0 = get_diagonal(0);
            double s = m0.sum();
            if (s > 0) {
                m0 /= s;
                set_diagonal(0, m0);
            }
        }
        if (M.size() > 1 && M[1].rows() > 0) {
            vec_t m1 = get_diagonal(1);
            double meanw = m1.mean();
            if (meanw > 0) {
                m1 /= meanw;
                set_diagonal(1, m1);
            }
        }
    }

    vec_t apply_Minv(int p, const vec_t& x) const {
        if (is_diagonal(p)) {
            vec_t y = x;
            for (Eigen::Index i = 0; i < y.size(); ++i) {
                double d = std::max(M[p].coeff(i, i), 1e-15);
                y[i] /= d;
            }
            return y;
        } else {
            if (!M_solver[p]) Rf_error("M[p] solver not initialized");
            return M_solver[p]->solve(x);
        }
    }

    /**
     * @brief Compute and store Cholesky factorization of M[p]
     * @param p Dimension index
     *
     * Factorizes the metric matrix for efficient repeated solves. Only needed
     * for non-diagonal metrics. Triggers an error if factorization fails (which
     * indicates the metric is not positive definite).
     */
    void refactor(int p) {
        if (!is_diagonal(p)) {
            auto solver = std::make_unique<Eigen::SimplicialLLT<spmat_t>>();
            solver->compute(M[p]);

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
                char buf[256];
                std::snprintf(buf, sizeof(buf),
                              "M[%d] factorization failed: %s (size=%ld, nnz=%ld)",
                              p, reason,
                              (long)M[p].rows(),
                              (long)M[p].nonZeros());
                Rf_error("%s", buf);
            }

            M_solver[p] = std::move(solver);
        } else {
            M_solver[p].reset();
        }
    }
};

/**
 * @brief Hodge Laplacians and related operators on the chain complex
 *
 * Manages the discrete exterior calculus operators: boundary maps, Hodge
 * Laplacians, and derived operators for heat diffusion and regularization.
 * The Laplacians encode the geometric structure of the complex and enable
 * solving partial differential equations discretized on the simplicial domain.
 */
struct laplacian_bundle_t {
    /// Boundary maps B[p]: C_p -> C_{p-1} for each dimension
    std::vector<spmat_t> B;

    /// Hodge Laplacians L[p] = B[p+1]^T M[p+1]^{-1} B[p+1] + M[p]^{-1} B[p] M[p-1] B[p]^T
    std::vector<spmat_t> L;

    /// Symmetrized vertex Laplacian L0_sym = M[0]^{-1/2} L[0] M[0]^{-1/2}
    spmat_t L0_sym;

    /// Edge conductances (optional weights) for dimension 0 Laplacian
    vec_t c1;

    /**
     * @brief Apply inverse metric from the right to boundary map
     * @param g Metric family
     * @param p Target dimension
     * @param Bp Boundary map B[p]
     * @return B[p] M[p]^{-1} (as a sparse matrix)
     *
     * Helper for assembling Laplacians. Assumes diagonal metric for efficiency.
     */
    spmat_t right_apply_Minv(const metric_family_t& g, int p, const spmat_t& Bp) const {
        spmat_t A(Bp.rows(), Bp.cols());
        A.reserve(Bp.nonZeros());
        for (int k=0; k<Bp.outerSize(); ++k) {
            for (spmat_t::InnerIterator it(Bp,k); it; ++it) {
                double d = std::max(g.M[p].coeff(it.row(), it.row()), 1e-15);
                A.insert(it.row(), it.col()) = it.value() / d;
            }
        }
        A.makeCompressed();
        return A;
    }

    /**
     * @brief Apply inverse metric from the left to a matrix
     * @param g Metric family
     * @param p Dimension index
     * @param T Input matrix
     * @return M[p]^{-1} T (as a sparse matrix)
     *
     * Helper for assembling Laplacians. Assumes diagonal metric for efficiency.
     */
    spmat_t left_apply_Minv(const metric_family_t& g, int p, const spmat_t& T) const {
        spmat_t D(T.rows(), T.cols());
        D.reserve(T.nonZeros());
        for (int k=0; k<T.outerSize(); ++k) {
            for (spmat_t::InnerIterator it(T,k); it; ++it) {
                double d = std::max(g.M[p].coeff(it.row(), it.row()), 1e-15);
                D.insert(it.row(), it.col()) = it.value() / d;
            }
        }
        D.makeCompressed();
        return D;
    }

    /**
     * @brief Build symmetrized (normalized) vertex Laplacian for spectral analysis
     *
     * @param g Metric family containing vertex masses M[0] and edge information
     *
     * @details
     * This function constructs the normalized graph Laplacian L0_sym, which is the
     * symmetrized version of the standard graph Laplacian. This form is particularly
     * useful for spectral analysis because it has eigenvalues in [0, 2] and satisfies
     * nice properties for diffusion processes.
     *
     * MATHEMATICAL FORMULATION:
     *
     * Given a weighted graph with:
     *   - Vertices V = {1, ..., n}
     *   - Edges E with conductances c_ij > 0
     *   - Vertex masses m_i > 0 (from metric M[0])
     *
     * We construct:
     * 1. Standard graph Laplacian: L_div
     * 2. Mass matrix diagonal: D = diag(m_1, ..., m_n)
     * 3. Normalized Laplacian: L_sym = D^{-1/2} L_div D^{-1/2}
     *
     * The normalized Laplacian has the property that for functions f on vertices:
     *   <f, L_sym f>_{ℓ²} = (1/2) Σ_{edges (i,j)} c_ij (f_i/√m_i - f_j/√m_j)²
     *
     * This measures the "smoothness" of f with respect to the geometry.
     */
    void build_L0_sym_if_needed(const metric_family_t& g) {
        // ========== VALIDATION: Check if we can build L0_sym ==========

        // Need at least 2 dimensions (vertices and edges)
        if (B.size() < 2) {
            L0_sym.resize(0, 0);
            return;
        }

        // Need a non-empty boundary operator B[1] (vertex-edge incidence)
        if (B[1].nonZeros() == 0) {
            L0_sym.resize(0, 0);
            return;
        }

        // Need vertex metric M[0]
        if (g.M.empty() || g.M[0].rows() == 0) {
            L0_sym.resize(0, 0);
            return;
        }

        // Need conductances c1 for each edge
        // c1[e] represents the "strength" or "weight" of edge e
        if (c1.size() != B[1].cols()) {
            L0_sym.resize(0, 0);
            return;
        }

        // ========== EXTRACT DIMENSIONS ==========

        const spmat_t& B1 = B[1];           // Boundary operator ∂₁: C₁ → C₀
        const Eigen::Index n = B1.rows();   // Number of vertices
        const Eigen::Index m = B1.cols();   // Number of edges

        // ========== BUILD STANDARD GRAPH LAPLACIAN L_div ==========

        // The graph Laplacian has the form:
        //   (L_div)_ii = Σ_{j∼i} c_ij  (sum of conductances of incident edges)
        //   (L_div)_ij = -c_ij          (negative conductance if i,j are adjacent)
        //   (L_div)_ij = 0              (zero if i,j not adjacent)
        //
        // This can be written as: L_div = B₁ C B₁ᵀ
        // where C = diag(c₁, ..., c_m) is the diagonal matrix of edge conductances.

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(4 * m);  // Each edge contributes 4 entries to L_div

        // Iterate over each edge
        for (int e = 0; e < m; ++e) {
            // ===== DECODE EDGE ENDPOINTS FROM BOUNDARY OPERATOR =====

            // The boundary operator B₁ has exactly 2 non-zero entries per column:
            //   B₁[i, e] = -1  (tail vertex of edge e)
            //   B₁[j, e] = +1  (head vertex of edge e)
            //
            // This encodes the oriented edge as: e = [i → j]

            int i = -1;  // Tail vertex (receives -1 in B₁)
            int j = -1;  // Head vertex (receives +1 in B₁)

            for (spmat_t::InnerIterator it(B1, e); it; ++it) {
                if (it.value() > 0) {
                    j = it.row();  // Positive coefficient → head vertex
                } else if (it.value() < 0) {
                    i = it.row();  // Negative coefficient → tail vertex
                }
            }

            // Skip malformed edges (should not happen in valid complex)
            if (i < 0 || j < 0) continue;

            // ===== GET EDGE CONDUCTANCE =====

            double c = c1[e];

            // Skip edges with non-positive conductance (degenerate or removed)
            if (c <= 0) continue;

            // ===== ADD CONTRIBUTIONS TO LAPLACIAN =====

            // For an edge (i,j) with conductance c, the Laplacian contributions are:
            //
            //   L[i,i] ← L[i,i] + c    (degree term at vertex i)
            //   L[j,j] ← L[j,j] + c    (degree term at vertex j)
            //   L[i,j] ← L[i,j] - c    (off-diagonal coupling)
            //   L[j,i] ← L[j,i] - c    (symmetric)
            //
            // This ensures L_div is symmetric and positive semi-definite.

            triplets.emplace_back(i, i,  c);   // Diagonal: vertex i degree
            triplets.emplace_back(j, j,  c);   // Diagonal: vertex j degree
            triplets.emplace_back(i, j, -c);   // Off-diagonal: coupling
            triplets.emplace_back(j, i, -c);   // Off-diagonal: coupling (symmetric)
        }

        // Assemble sparse matrix from triplets
        // Note: setFromTriplets automatically sums contributions to same (row, col)
        spmat_t L0_div(n, n);
        L0_div.setFromTriplets(triplets.begin(), triplets.end());

        // ========== BUILD NORMALIZATION MATRIX D^{-1/2} ==========

        // The vertex masses m_i are stored on the diagonal of M[0].
        // We construct the diagonal matrix:
        //   D^{-1/2} = diag(1/√m₁, ..., 1/√m_n)
        //
        // This will be used to normalize the Laplacian.

        spmat_t Dm05(n, n);
        Dm05.reserve(Eigen::VectorXi::Constant(n, 1));  // Reserve space for diagonal

        for (Eigen::Index i = 0; i < n; ++i) {
            // Extract vertex mass (with safety threshold to avoid division by zero)
            double m = std::max(g.M[0].coeff(i, i), 1e-15);

            // Insert 1/√m_i on the diagonal
            Dm05.insert(i, i) = 1.0 / std::sqrt(m);
        }
        Dm05.makeCompressed();

        // ========== COMPUTE NORMALIZED LAPLACIAN ==========

        // The normalized (symmetrized) Laplacian is:
        //   L_sym = D^{-1/2} L_div D^{-1/2}
        //
        // This has several nice properties:
        // 1. Symmetric: L_sym = L_symᵀ
        // 2. Eigenvalues in [0, 2]
        // 3. Smallest eigenvalue is 0 (with eigenvector ∝ √masses)
        // 4. Number of 0 eigenvalues = number of connected components
        // 5. Well-conditioned for iterative solvers
        //
        // The normalization makes vertices with large mass contribute less
        // to the energy, balancing the influence of high-degree and low-degree vertices.

        L0_sym = Dm05 * L0_div * Dm05;
    }

    /**
     * @brief Assemble all Hodge Laplacians from boundary maps and metric
     * @param g Metric family defining the inner products
     *
     * Computes L[p] for all dimensions using the Hodge decomposition formula.
     * Each Laplacian combines contributions from boundary maps above and below
     * dimension p, weighted by the metric.
     */
    void assemble_all(const metric_family_t& g) {
        const int P = static_cast<int>(g.M.size()) - 1;
        L.assign(P + 1, spmat_t());

        for (int p = 0; p <= P; ++p) {
            const Eigen::Index n = g.M[p].rows();

            spmat_t up, down;
            bool has_up = false, has_down = false;

            // Up term: ∂_{p+1} ∂_{p+1}^* = B[p+1] M[p+1]^{-1} B[p+1]^T M[p]
            if (p + 1 <= P && B.size() > static_cast<size_t>(p + 1) && B[p+1].nonZeros() > 0) {
                spmat_t temp = B[p+1] * left_apply_Minv(g, p + 1, B[p+1].transpose());
                up = temp * g.M[p];
                has_up = true;
            }

            // Down term: ∂_p^* ∂_p = M[p]^{-1} B[p]^T M[p-1] B[p]
            if (p >= 1 && B.size() > static_cast<size_t>(p) && B[p].nonZeros() > 0) {
                spmat_t temp = B[p].transpose() * g.M[p - 1] * B[p];
                down = left_apply_Minv(g, p, temp);
                has_down = true;
            }

            // Assemble Laplacian
            if (has_up && has_down) {
                L[p] = up + down;
            } else if (has_up) {
                L[p] = up;
            } else if (has_down) {
                L[p] = down;
            } else {
                L[p] = spmat_t(n, n);
            }
        }

        build_L0_sym_if_needed(g);
    }

    /**
     * @brief Apply heat kernel exp(-t L[p]) to a vector
     * @param p Dimension index
     * @param x Input vector
     * @param t Diffusion time
     * @param m_steps Number of backward Euler steps for approximation
     * @return Approximation to exp(-t L[p]) x
     *
     * Uses m-step backward Euler time stepping to approximate the heat
     * semigroup. Larger m_steps give better approximations at higher cost.
     */
    vec_t heat_apply(int p, const vec_t& x, double t, int m_steps=8) const {
        if (p<0 || p>=static_cast<int>(L.size())) Rf_error("heat_apply: bad p");
        if (t<=0) return x;
        spmat_t A = L[p];
        spmat_t I(A.rows(), A.cols()); I.setIdentity();
        double h = t / std::max(1, m_steps);
        A = I + h * A;
        Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;
        cg.setTolerance(1e-8); cg.setMaxIterations(1000); cg.compute(A);
        vec_t y = x;
        for (int s=0; s<m_steps; ++s) y = cg.solve(y);
        return y;
    }

    /**
     * @brief Solve Tikhonov regularization problem (I + eta L[p]) x = b
     * @param p Dimension index
     * @param b Right-hand side vector
     * @param eta Regularization parameter
     * @param tol Conjugate gradient tolerance
     * @param maxit Maximum CG iterations
     * @return Solution vector x
     *
     * Standard ridge regression or Tikhonov regularization against the
     * Laplacian. Used for signal smoothing and fitting with geometric
     * regularization on the complex.
     */
    vec_t tikhonov_solve(int p, const vec_t& b, double eta,
                         double tol=1e-8, int maxit=1000) const {
        if (p<0 || p>=static_cast<int>(L.size())) Rf_error("tikhonov: bad p");
        spmat_t A = L[p];
        spmat_t I(A.rows(), A.cols()); I.setIdentity();
        A = I + eta * A;
        Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;
        cg.setTolerance(tol); cg.setMaxIterations(maxit); cg.compute(A);
        return cg.solve(b);
    }
};

/**
 * @brief Density functions on chains at each dimension
 *
 * Maintains density (or mass distribution) functions as coefficient vectors
 * for chains at each dimension of the complex. In the Morse-Smale regression
 * framework, these densities evolve according to geometric flows and encode
 * the probability distribution of data across the simplicial structure.
 */
struct density_family_t {
    /// Density vector for each dimension p (length n_p)
    std::vector<vec_t> rho;

    /**
     * @brief Initialize density vectors for all dimensions
     * @param n_by_dim Vector where n_by_dim[p] is the number of p-simplices
     *
     * Allocates zero-initialized density vectors with appropriate dimensions
     * for the entire chain complex.
     */
    void init_dims(const std::vector<index_t>& n_by_dim) {
        rho.resize(n_by_dim.size());
        for (index_t p = 0; p < n_by_dim.size(); ++p)
            rho[p] = vec_t::Zero((Eigen::Index)n_by_dim[p]);
    }
};

/**
 * @brief Signal observations and fitted values
 *
 * Manages the observed signal and maintains history of fitted signal values
 * during iterative regression. The signal y typically represents scalar
 * observations at vertices (0-simplices), while y_hat_hist tracks the
 * evolution of fitted values through the optimization procedure.
 */
struct signal_state_t {
    /// Observed signal vector (typically at vertices)
    vec_t y; ///< Original response

    /// History of fitted signal vectors across iterations
    std::vector<vec_t> y_hat_hist; ///< Fitted values history

    signal_state_t() = default;

    /**
     * @brief Clear the fitted value history
     *
     * Useful for resetting state between independent runs or freeing memory
     * when history is no longer needed.
     */
    void clear_history() { y_hat_hist.clear(); }
};

// ================================================================
// MAIN DATA STRUCTURE
// ================================================================

/**
 * @brief Riemannian simplicial complex for geometric regression
 *
 * Central data structure representing a simplicial complex equipped with
 * a Riemannian metric, Hodge Laplacian operators, and signal processing
 * capabilities for conditional expectation estimation.
 */
struct riem_dcx_t {
    // ----------------------------------------------------------------
    // Core Complex Structure
    // ----------------------------------------------------------------

    int pmax = 1;                        ///< Maximum simplex dimension
    std::vector<simplex_table_t> S;      ///< Simplex tables by dimension
    std::vector<star_table_t> stars;     ///< Star neighborhoods

    // ----------------------------------------------------------------
    // Geometric Structure
    // ----------------------------------------------------------------

    metric_family_t g;                   ///< Mass matrices (Riemannian metric)
    laplacian_bundle_t L;                ///< Hodge Laplacian operators
    density_family_t rho;                ///< Probability densities

    // ----------------------------------------------------------------
    // Signal State
    // ----------------------------------------------------------------

    signal_state_t sig;                  ///< Response and fitted values

    // ----------------------------------------------------------------
    // Geometric Data (for regression)
    // ----------------------------------------------------------------

    std::vector<double> edge_lengths;           ///< Edge lengths (1-skeleton)
    std::vector<double> reference_measure;      ///< Reference probability measure

    // ----------------------------------------------------------------
    // Convergence Tracking (for regression)
    // ----------------------------------------------------------------

    bool converged = false;                      ///< Convergence flag
    int n_iterations = 0;                        ///< Number of iterations performed
    double final_response_change = 0.0;          ///< Final response convergence metric
    double final_density_change = 0.0;           ///< Final density convergence metric
    std::vector<double> response_changes;        ///< Response change history
    std::vector<double> density_changes;         ///< Density change history

    // ================================================================
    // CONSTRUCTION METHODS
    // ================================================================

    /**
     * @brief Build nerve complex from kNN covering (general purpose)
     */
    void build_nerve_from_knn(
        const spmat_t& X,
        index_t k,
        index_t max_p,
        bool use_counting_measure,
        double density_normalization
    );

    /**
     * @brief Fit kNN Riemannian graph regression model
     */
    void fit_knn_riem_graph_regression(
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
    );

    /**
     * @brief Legacy wrapper for backward compatibility
     */
    void build_knn_riem_dcx(
        const spmat_t& X,
        const vec_t& y,
        index_t k,
        index_t max_p,
        bool use_counting_measure,
        double density_normalization
    );

    // ================================================================
    // ITERATION HELPER METHODS
    // ================================================================

    /**
     * @brief Apply damped heat diffusion to vertex density
     */
    vec_t apply_damped_heat_diffusion(
        const vec_t& rho_current,
        double t,
        double beta
    );

    /**
     * @brief Update edge densities from evolved vertex densities
     */
    void update_edge_densities_from_vertices();

    /**
     * @brief Modulate edge densities based on response variation
     */
    void apply_response_coherence_modulation(
        const vec_t& y_hat,
        double gamma
    );

    /**
     * @brief Smooth response via spectral filtering with GCV
     */
    gcv_result_t smooth_response_via_spectral_filter(
        const vec_t& y,
        int n_eigenpairs,
        rdcx_filter_type_t filter_type
    );

    /**
     * @brief Check convergence criteria
     */
    convergence_status_t check_convergence(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        const std::vector<vec_t>& rho_prev,
        const std::vector<vec_t>& rho_curr,
        double epsilon_y,
        double epsilon_rho,
        int iteration,
        int max_iterations
    );

    // ================================================================
    // EXISTING METHODS
    // ================================================================

    /**
     * @brief Initialize all structures with given dimensions
     * @param pmax_ Maximum chain dimension
     * @param n_by_dim Number of simplices at each dimension
     *
     * Allocates and initializes all component structures. Must be called
     * before adding simplices or performing computations.
     */
    void init_dims(int pmax_, const std::vector<index_t>& n_by_dim) {
        pmax = pmax_;
        S.resize(pmax+1);
        stars.resize(pmax+1);
        //inn.resize(pmax+1);
        rho.init_dims(n_by_dim);
        g.init_dims(n_by_dim);
        L.B.assign(pmax+1, spmat_t());
        L.L.assign(pmax+1, spmat_t());
    }

    /**
     * @brief Build boundary operator B[1] from stored edge list
     *
     * Constructs the vertex-edge incidence matrix from the edge simplex table.
     * Each edge (i,j) with i < j contributes (+1) at vertex j and (-1) at
     * vertex i in the corresponding column of B[1].
     */
    void build_incidence_from_edges() {
        if (S.size()<2) return;
        const auto& edges = S[1].simplex_verts;
        const Eigen::Index n = static_cast<Eigen::Index>(S[0].size());
        const Eigen::Index m = static_cast<Eigen::Index>(edges.size());
        spmat_t B1(n, m);
        B1.reserve(Eigen::VectorXi::Constant(m,2));
        for (Eigen::Index e=0; e<m; ++e) {
            auto v = edges[e];
            if (v.size()!=2)
                Rf_error("Edge must have exactly 2 vertices");
            Eigen::Index i = static_cast<Eigen::Index>(v[0]);
            Eigen::Index j = static_cast<Eigen::Index>(v[1]);
            if (i<j) { B1.insert(j,e) =  1.0; B1.insert(i,e) = -1.0; }
            else     { B1.insert(i,e) =  1.0; B1.insert(j,e) = -1.0; }
        }
        B1.makeCompressed();
        if (L.B.size()<2) L.B.resize(2);
        L.B[1] = std::move(B1);
    }

    /**
     * @brief Assemble all Laplacian operators from current metric
     *
     * Updates all Hodge Laplacians based on the current state of the metric
     * and boundary maps. Should be called after any change to the metric.
     */
    void assemble_operators() { L.assemble_all(g); }

    /**
     * @brief Perform one iteration of density evolution and signal fitting
     * @param t_vertex Heat diffusion time for vertex density
     * @param eta_y Tikhonov regularization parameter for signal fitting
     * @param with_edge_density Whether to compute edge densities (unused)
     * @param rebuild_knn Whether to rebuild nearest neighbor structure (unused)
     *
     * Executes one step of the iterative algorithm: evolves vertex density
     * via heat diffusion from empirical distribution, fits signal with
     * regularization if available, updates metric from density, and
     * reassembles operators. This implements the core computational loop
     * of the Morse-Smale regression method.
     */
    void iteration_step(double t_vertex, double eta_y,
                        bool with_edge_density=false,
                        bool rebuild_knn=false) {
        (void)with_edge_density; (void)rebuild_knn;
        const Eigen::Index n = static_cast<Eigen::Index>(S[0].size());
        vec_t emp = vec_t::Ones(n);
        emp /= std::max<double>(emp.sum(), 1e-15);
        rho.rho[0] = L.heat_apply(0, emp, t_vertex, 8);

        if (sig.y.size()==n) {
            vec_t yhat = L.tikhonov_solve(0, sig.y, eta_y);
            sig.y_hat_hist.push_back(yhat);
        }

        if (!g.M.empty() && g.M[0].rows()==n) {
            const double alpha = 0.2;
            rho.rho[0] = rho.rho[0].cwiseMax(1e-15);
            vec_t m0 = g.get_diagonal(0);
            m0 = (1.0 - alpha) * m0 + alpha * rho.rho[0];
            g.set_diagonal(0, m0);
            g.normalize();
        }

        assemble_operators();
    }

    /**
     * @brief Extend the complex by one dimension, preserving all existing structure
     *
     * @details
     * Appends a new chain degree to the discrete chain complex, incrementing pmax
     * by one and initializing all data structures for the new dimension. This
     * operation preserves all existing simplex tables, metric structures, boundary
     * operators, and Laplacians at lower dimensions. The function is designed for
     * incremental complex construction where simplices are built dimension by
     * dimension.
     *
     * MOTIVATION:
     * During nerve complex construction from k-NN neighborhoods, we discover
     * simplices dimension by dimension: first vertices, then edges (by checking
     * pairwise neighborhood intersections), then triangles (by checking triple
     * intersections), and so forth. Each dimension depends on having the previous
     * dimension fully constructed to check face existence. This function enables
     * that incremental workflow without destroying lower-dimensional data.
     *
     * OPERATION:
     * Let p_new = pmax + 1 denote the new dimension being added. The function:
     *
     * 1. Increments pmax to p_new
     * 2. Resizes all dimension-indexed containers to accommodate p_new:
     *    - S (simplex tables)
     *    - stars (upward incidence)
     *    - inn (local inner products)
     *    - rho.rho (density vectors)
     *    - g.M (metric diagonal and matrices)
     *    - g.M_solver (Cholesky factorizations, if maintained)
     *    - L.B (boundary operators)
     *    - L.L (Hodge Laplacians)
     *
     * 3. Initializes the new dimension p_new with default values:
     *    - S[p_new]: Empty simplex table (to be populated by caller)
     *    - rho.rho[p_new]: Zero density vector of length n_new
     *    - g.M[p_new]: Sparse diagonal matrix with ones on diagonal
     *    - g.M_solver[p_new]: Null pointer (no factorization yet)
     *    - stars[p_new]: Star table with one bucket per (p_new-1)-simplex
     *    - inn[p_new]: Empty local inner product blocks
     *    - L.B[p_new]: Empty boundary operator matrix (shape: |C_{p_new-1}| × n_new)
     *    - L.L[p_new]: Empty Laplacian (to be assembled after faces are defined)
     *
     * USAGE PATTERN:
     * The typical usage in complex construction is:
     *
     * @code
     * // Start with vertices only (dimension 0)
     * riem_dcx_t dcx;
     * std::vector<index_t> n_by_dim = {n_vertices};
     * dcx.init_dims(0, n_by_dim);
     * // ... populate S[0] with vertices ...
     *
     * // Discover edges exist
     * std::vector<std::vector<index_t>> edge_list = discover_edges();
     * dcx.extend_by_one_dim(edge_list.size());  // Now pmax = 1
     * // ... populate S[1] with edges ...
     * // ... build L.B[1] (vertex-edge incidence) ...
     *
     * // Discover triangles exist
     * std::vector<std::vector<index_t>> triangle_list = discover_triangles();
     * dcx.extend_by_one_dim(triangle_list.size());  // Now pmax = 2
     * // ... populate S[2] with triangles ...
     * // ... build L.B[2] (edge-triangle incidence) ...
     * @endcode
     *
     * STAR TABLE SIZING:
     * The star table stars[p_new] is used to record which (p_new+1)-simplices
     * contain each p_new-simplex. However, at the time of creating dimension p_new,
     * no (p_new+1)-simplices exist yet. The star table is sized based on the
     * number of (p_new-1)-simplices, which is what stars[p_new-1] indexes over.
     *
     * Specifically:
     *   - stars[p_new].star_over has length S[p_new-1].size()
     *   - Each entry stars[p_new].star_over[sigma_id] will eventually list the
     *     p_new-simplices that have (p_new-1)-simplex sigma_id in their boundary
     *   - For p_new = 0, stars[0] has length 0 (no (-1)-simplices)
     *
     * METRIC INITIALIZATION:
     * The new dimension starts with identity metric (diagonal ones). This is a
     * safe default that can be overwritten by the caller based on geometric
     * measurements (e.g., intersection measures for nerve complexes). The metric
     * determines inner products on chain spaces and must be positive definite.
     *
     * BOUNDARY OPERATOR:
     * The boundary map B[p_new]: C_{p_new} → C_{p_new-1} is allocated with the
     * correct dimensions (rows = |C_{p_new-1}|, cols = |C_{p_new}| = n_new) but
     * contains no entries initially. The caller must populate this after defining
     * the faces of each p_new-simplex. For a p_new-simplex with faces
     * {f_0, ..., f_p_new}, the boundary operator has entries:
     *     B[p_new](f_i, simplex_id) = (-1)^i
     *
     * LAPLACIAN ASSEMBLY:
     * The Hodge Laplacian L[p_new] is left empty and must be computed via
     * assemble_operators() after:
     *   1. All boundary operators B[p_new+1] and B[p_new] are populated
     *   2. The metric g.M[p_new] reflects the actual geometry
     *
     * The Laplacian at dimension p is given by the Hodge decomposition:
     *     L[p] = B[p+1]^T M[p+1]^{-1} B[p+1] + M[p]^{-1} B[p] M[p-1] B[p]^T
     * where the first term represents coboundary (curl) contribution and the
     * second term represents boundary (divergence) contribution.
     *
     * INVARIANTS MAINTAINED:
     * After calling extend_by_one_dim(n_new):
     *   - pmax has increased by 1
     *   - S.size() == stars.size() == inn.size() == pmax + 1
     *   - g.M.size() == L.B.size() == L.L.size() == pmax + 1
     *   - rho.rho.size() == pmax + 1
     *   - S[pmax].size() == 0 (empty, to be populated)
     *   - g.M[pmax] is diagonal with n_new ones
     *   - L.B[pmax] is empty sparse matrix of shape (|C_{pmax-1}|, n_new)
     *   - All structures at dimensions 0, ..., pmax-1 are unchanged
     *
     * EFFICIENCY:
     * This operation is O(n_new) for initializing the new dimension's data
     * structures, plus O(pmax) for resizing dimension-indexed containers.
     * Critically, it does NOT rebuild or copy data from lower dimensions,
     * making it suitable for use in tight loops during complex construction.
     *
     * ALTERNATIVE TO init_dims():
     * The older init_dims(pmax, n_by_dim) function resizes all dimensions at once
     * and may require rebuilding lower-dimensional structures if called repeatedly.
     * In contrast, extend_by_one_dim() is designed for incremental growth:
     *   - init_dims(): Use once at the start to establish dimension 0
     *   - extend_by_one_dim(): Use repeatedly to grow the complex upward
     *
     * THREAD SAFETY:
     * This function modifies the riem_dcx_t object and is NOT thread-safe. Do not
     * call concurrently from multiple threads on the same object.
     *
     * @param n_new Number of simplices at the new dimension p_new = pmax + 1.
     *              Must be non-negative. If zero, the dimension is added but
     *              contains no simplices (useful for maintaining consistent
     *              dimension indexing even when higher simplices don't exist).
     *
     * @post pmax is incremented by 1
     * @post All dimension-indexed containers have size pmax + 1
     * @post S[pmax] is empty (size 0)
     * @post g.M[pmax] is diagonal identity matrix of size n_new × n_new
     * @post rho.rho[pmax] is zero vector of length n_new
     * @post L.B[pmax] is empty sparse matrix of shape (|S[pmax-1]|, n_new)
     * @post L.L[pmax] is empty (to be assembled later)
     * @post All data at dimensions 0, ..., pmax-1 (pre-increment) are unchanged
     *
     * @note The caller must subsequently:
     *       1. Populate S[pmax].simplex_verts with actual simplices
     *       2. Update S[pmax].id_of for simplex lookup
     *       3. Set g.M[pmax] diagonal entries based on geometric measurements
     *       4. Build L.B[pmax] by iterating over simplices and their faces
     *       5. Call assemble_operators() to compute L.L[pmax]
     *
     * @note If n_new = 0, the dimension exists but is empty. This is valid and
     *       sometimes useful when the construction algorithm discovers that no
     *       simplices of a given dimension exist (e.g., no triangles in a sparse
     *       graph). Having the dimension allocated prevents index errors in
     *       downstream code that assumes dimension p_new exists.
     *
     * @see init_dims() - Initialize all dimensions at once (use for dimension 0)
     * @see assemble_operators() - Compute Laplacians after boundary maps are built
     * @see build_incidence_from_edges() - Example of building B[1] from S[1]
     *
     * @example
     * // Incremental construction of a 2-dimensional complex
     * riem_dcx_t dcx;
     *
     * // Start with 5 vertices
     * std::vector<index_t> n_by_dim = {5};
     * dcx.init_dims(0, n_by_dim);
     * for (index_t i = 0; i < 5; ++i) {
     *     dcx.S[0].add_simplex({i});
     * }
     *
     * // Add 8 edges
     * dcx.extend_by_one_dim(8);  // pmax is now 1
     * dcx.S[1].add_simplex({0, 1});
     * dcx.S[1].add_simplex({1, 2});
     * // ... add remaining edges ...
     * dcx.build_incidence_from_edges();
     *
     * // Add 3 triangles
     * dcx.extend_by_one_dim(3);  // pmax is now 2
     * dcx.S[2].add_simplex({0, 1, 2});
     * dcx.S[2].add_simplex({1, 2, 3});
     * dcx.S[2].add_simplex({2, 3, 4});
     * // ... build B[2] from triangle faces ...
     *
     * // Compute all Laplacians
     * dcx.assemble_operators();
     */
    void extend_by_one_dim(index_t n_new) {
        const int p_new = pmax + 1;

        // Increment maximum dimension and resize all containers
        pmax = p_new;
        S.resize(pmax + 1);
        stars.resize(pmax + 1);
        //inn.resize(pmax + 1);

        // Resize density family
        rho.rho.resize(pmax + 1);

        // Resize metric family
        g.M.resize(pmax + 1);
        g.M_solver.resize(pmax + 1);

        // Resize operator family
        L.B.resize(pmax + 1);
        L.L.resize(pmax + 1);

        // Initialize combinatorics at new dimension (empty for now)
        S[p_new].simplex_verts.clear();
        S[p_new].id_of.clear();

        // Initialize density at new dimension (zeros)
        rho.rho[p_new] = vec_t::Zero(static_cast<Eigen::Index>(n_new));

        // Initialize metric at new dimension (identity)
        g.M[p_new] = spmat_t(static_cast<Eigen::Index>(n_new),
                             static_cast<Eigen::Index>(n_new));
        g.M[p_new].reserve(Eigen::VectorXi::Constant(static_cast<int>(n_new), 1));
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n_new); ++i) {
            g.M[p_new].insert(i, i) = 1.0;
        }
        g.M[p_new].makeCompressed();
        g.M_solver[p_new].reset();

        // Initialize star table at new dimension
        if (p_new >= 1) {
            stars[p_new].resize(S[p_new - 1].size());
        } else {
            stars[p_new].resize(0);
        }

        // Initialize local inner products
        // inn[p_new].resize((p_new >= 1) ? S[p_new - 1].size() : 0);

        // Initialize boundary operator (empty matrix with correct shape)
        const Eigen::Index rows = (p_new >= 1) ?
            static_cast<Eigen::Index>(S[p_new - 1].size()) : 0;
        const Eigen::Index cols = static_cast<Eigen::Index>(n_new);
        L.B[p_new] = spmat_t(rows, cols);

        // Initialize Laplacian (empty, to be assembled later)
        L.L[p_new] = spmat_t();
    }

private:

    std::vector<std::unordered_set<index_t>> neighbor_sets;

    /**
     * @brief Cached spectral decomposition of vertex Laplacian
     *
     * Stores eigenvalues and eigenvectors of L[0] to avoid recomputation.
     * Updated whenever the Laplacian is reassembled during iteration.
     */
    struct spectral_cache_t {
        vec_t eigenvalues;           ///< Eigenvalues in ascending order
        Eigen::MatrixXd eigenvectors; ///< Corresponding eigenvectors (columns)
        bool is_valid;               ///< True if cache contains current L[0] spectrum
        double lambda_2;             ///< Spectral gap (second smallest eigenvalue)

        spectral_cache_t() : is_valid(false), lambda_2(0.0) {}

        void invalidate() { is_valid = false; }
    } spectral_cache;

    // ================================================================
    // INTERNAL HELPERS
    // ================================================================

    /**
     * @brief Compute and cache spectral decomposition of vertex Laplacian
     */
    void compute_spectral_decomposition(int n_eigenpairs = -1);

    /**
     * @brief Automatically select diffusion and damping parameters
     */
    void select_diffusion_parameters(
        double& t_diffusion,
        double& beta_damping,
        bool verbose = false
        );

    /**
     * @brief Initialize reference measure for density computation
     */
    void initialize_reference_measure(
        const std::vector<std::vector<index_t>>& knn_neighbors,
        const std::vector<std::vector<double>>& knn_distances,
        bool use_counting_measure,
        double density_normalization
    );

    /**
     * @brief Compute initial densities from reference measure
     */
    void compute_initial_densities();

    /**
     * @brief Initialize metric from densities
     */
    void initialize_metric_from_density();

    /**
     * @brief Update metric from current densities
     */
    void update_metric_from_density();

    /**
     * @brief Compute full edge mass matrix with triple intersections
     */
    void compute_edge_mass_matrix();

    /**
     * @brief Update edge mass matrix from current vertex densities
     */
    void update_edge_mass_matrix();

    /**
     * @brief Compute inner product between edges in a star
     */
    double compute_edge_inner_product(
        index_t e1,
        index_t e2,
        index_t vertex_i
    ) const;

    /**
     * @brief Compute simplex volume via Cayley-Menger determinant
     */
    double compute_simplex_volume(
        const std::vector<index_t>& vertices
    ) const;
};
