#pragma once

#include "iknn_vertex.hpp"

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

/// Sentinel value indicating no vertex (used in self-loops)
constexpr index_t NO_VERTEX = std::numeric_limits<index_t>::max();

/// Sentinel value indicating no edge (used during construction)
constexpr index_t NO_EDGE = std::numeric_limits<index_t>::max();

// ================================================================
// SUPPORTING ENUMERATIONS AND STRUCTURES
// ================================================================

/**
 * @brief Filter type enumeration for spectral filtering
 */
enum class rdcx_filter_type_t {
    HEAT_KERNEL,        ///< f(λ) = exp(-ηλ)
    TIKHONOV,           ///< f(λ) = 1/(1 + ηλ)
    CUBIC_SPLINE,       ///< f(λ) = 1/(1 + ηλ²)
    GAUSSIAN,           ///< f(λ) = exp(-ηλ²)
    EXPONENTIAL,        ///< f(λ) = exp(-η·√λ)
    BUTTERWORTH         ///< f(λ) = 1/(1+(λ/η)^(2n))
};

/**
 * @brief Result structure for GCV-based spectral filtering
 */
struct gcv_result_t {
    double eta_optimal;              ///< Optimal smoothing parameter
    double gcv_optimal;              ///< Optimal GCV
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
 * @brief Enhanced convergence check with additional diagnostics
 *
 * This extended version provides more detailed diagnostic information,
 * including per-dimension density changes and convergence predictions.
 * Useful for research and debugging, but not necessary for production use.
 *
 * @note This function is optional and not currently used in the main
 *       iteration loop. It demonstrates how more detailed convergence
 *       monitoring could be implemented if needed.
 */
struct detailed_convergence_status_t {
    bool converged;
    double response_change;
    std::vector<double> density_changes_by_dim;
    double max_density_change;
    int iteration;
    std::string message;

    // Rate diagnostics
    double response_change_rate;
    double density_change_rate;
    int estimated_iterations_remaining;

    // GCV diagnostics
    double gcv_current;
    double gcv_change;        // Current - previous (negative means improvement)
    double gcv_change_rate;   // Rate of GCV change
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
     */
    spmat_t right_apply_Minv(const metric_family_t& g, int p, const spmat_t& Bp) const;

    /**
     * @brief Apply inverse metric from the left to a matrix
     */
    spmat_t left_apply_Minv(const metric_family_t& g, int p, const spmat_t& T) const;

    /**
     * @brief Build symmetrized (normalized) vertex Laplacian for spectral analysis
     */
    void build_L0_sym_if_needed(const metric_family_t& g);

    /**
     * @brief Assemble all Hodge Laplacians from boundary maps and metric
     */
    void assemble_all(const metric_family_t& g);

    /**
     * @brief Apply heat kernel exp(-t L[p]) to a vector
     */
    vec_t heat_apply(
        int p,
        const vec_t& x,
        double t,
        int m_steps
        ) const;

    /**
     * @brief Solve Tikhonov regularization problem (I + eta L[p]) x = b
     */
    vec_t tikhonov_solve(
        int p,
        const vec_t& b,
        double eta,
        double tol,
        int maxit
        ) const;
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

struct neighbor_info_t {
    index_t vertex_index;    ///< Vertex index
    index_t simplex_index;   ///< Simplex index in S[p]
    size_t isize;            ///< Neighborhood intersection size
    double dist;             ///< Geometric distance
    double density;          ///< Density from ρ[p]
};

struct gamma_selection_result_t {
    double gamma_optimal;
    double gcv_optimal;
    std::vector<double> gamma_grid;
    std::vector<double> gcv_scores;
    vec_t y_hat_optimal;
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

    /// Cofaces of vertices: [i][0] is vertex i, [i][j>0] are incident edge
    std::vector<std::vector<neighbor_info_t>> vertex_cofaces;

    // vertex_cofaces[i][0] = {
    //     .index = i,                        // The vertex itself
    //     .simplex_index = i,
    //     .isize = |N̂_k(i)|,
    //     .dist = 0.0,
    //     .density = ρ_0[i]
    // };

    // vertex_cofaces[i][j] = {  // j > 0
    //     .index = neighbor_vertex_idx,     // The other endpoint of edge [i, neighbor_vertex_idx]
    //     .simplex_index = edge_idx,        // For density lookup
    //     .isize = |N̂_k(i) ∩ N̂_k(neighbor)|,
    //     .dist = d(i, neighbor),
    //     .density = ρ_1[edge_idx]          // Edge density
    // };

    // Edge registry for indexing
    std::vector<std::array<index_t, 2>> edge_registry;  ///< edge_registry[e] = {i, j}

    /// Cofaces of edges: [e][0] is edge e, [e][j>0] are incident triangles
    std::vector<std::vector<neighbor_info_t>> edge_cofaces;

    // edge_cofaces[e][0] = {
    //     .vertex_index = -1,                  // Special marker: -1
    //     .simplex_index = e,                  // The edge index
    //     .isize = 0,                          // Not applicable
    //     .dist = edge_lengths[e],             // Edge length (useful!)
    //     .density = ρ_1[e]                    // Edge density
    // };

    // edge_cofaces[e][j] = {  // j > 0
    //     .vertex_index = third_vertex_idx,    // The third vertex forming triangle
    //     .simplex_index = triangle_idx,       // Triangle index in S[2]
    //     .isize = triangle_intersection_size,
    //     .dist = 0.0,                         // Or some measure
    //     .density = ρ_2[triangle_idx]         // Triangle density
    // };

    // ----------------------------------------------------------------
    // Geometric Structure
    // ----------------------------------------------------------------

    metric_family_t g;                   ///< Mass matrices (Riemannian metric)
    laplacian_bundle_t L;                ///< Hodge Laplacian operators


    // ----------------------------------------------------------------
    // Signal State
    // ----------------------------------------------------------------

    signal_state_t sig;                  ///< Response and fitted values

    // GCV tracking across iterations
    struct gcv_history_t {
        std::vector<gcv_result_t> iterations;

        void clear() {
            iterations.clear();
        }

        void add(const gcv_result_t& result) {
            iterations.push_back(result);
        }

        size_t size() const {
            return iterations.size();
        }

        bool empty() const {
            return iterations.empty();
        }
    } gcv_history;

    struct density_history_t {
        std::vector<vec_t> rho_vertex;  // One vector per iteration

        void clear() {
            rho_vertex.clear();
        }

        void add(const vec_t& rho) {
            rho_vertex.push_back(rho);
        }
    } density_history;

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

    // Gamma selection result (empty if not performed)
    gamma_selection_result_t gamma_selection_result;
    bool gamma_was_auto_selected = false;

    // ================================================================
    // CONSTRUCTION METHODS
    // ================================================================

    /**
     * @brief Fit Riemannian graph regression model using iterative geometric refinement
     *
     * @param X Feature matrix (sparse, n x d)
     * @param y Response vector (length n)
     * @param k Number of nearest neighbors
     * @param use_counting_measure If true, use counting measure; else distance-based
     * @param density_normalization Normalization mode for densities
     * @param t_diffusion Diffusion time parameter (if <= 0, auto-select using t_scale_factor)
     * @param beta_damping Damping parameter (if <= 0, auto-select using beta_coefficient_factor)
     * @param gamma_modulation Response-coherence modulation strength
     * @param t_scale_factor Scale factor for auto-selecting t (controls diffusion scale t*lambda_2)
     * @param beta_coefficient_factor Factor for auto-selecting beta (controls damping coefficient beta*t)
     * @param n_eigenpairs Number of eigenpairs for spectral filtering
     * @param filter_type Type of spectral filter to apply
     * @param epsilon_y Convergence threshold for response
     * @param epsilon_rho Convergence threshold for density
     * @param max_iterations Maximum number of iterations
     * @param max_ratio_threshold Maximum ratio threshold for edge pruning
     * @param threshold_percentile Percentile for adaptive threshold
     * @param density_alpha Power for distance-based density (default 1.5)
     * @param density_epsilon Regularization for distance-based density (default 1e-10)
     * @param test_stage For debugging: stop at specific stage (default -1, run all)
     * @param verbose Print progress information
     */
    void fit_rdgraph_regression(
        const spmat_t& X,
        const vec_t& y,
        index_t k,
        bool use_counting_measure,
        double density_normalization,
        double t_diffusion,
        double beta_damping,
        double gamma_modulation,
        double t_scale_factor,
        double beta_coefficient_factor,
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
        );

    /**
     * @brief Build boundary operator B[1] from edge_registry
     */
    void build_boundary_operator_from_edges();

    /**
     * @brief Build boundary operator B[2] from triangles using edge_cofaces only
     */
    void build_boundary_operator_from_triangles();

    /**
     * @brief Compute the coboundary operator ∂₁* with full non-diagonal metric
     *
     * This implements the solution of: M₁ x = B₁ᵀ M₀ f where M₁ may be
     *   non-diagonal (full edge mass matrix with correlations).
     */
    vec_t compute_coboundary_del1star_full(
        const vec_t& f,
        bool use_iterative,
        double cg_tol = 1e-10,
        int cg_maxit = 1000
        ) const;

    /**
     * @brief Compute extremality scores using full non-diagonal metric
     */
    vec_t compute_all_extremality_scores_full(
        const vec_t& f,
        bool use_iterative = false,
        double cg_tol = 1e-10,
        int cg_maxit = 1000
        ) const;

    // ================================================================
    // ITERATION HELPER METHODS
    // ================================================================

    /**
     * @brief Select optimal gamma parameter via first-iteration GCV evaluation
     */
    gamma_selection_result_t select_gamma_first_iteration(
        const vec_t& y,
        const std::vector<double>& gamma_grid,
        int n_eigenpairs,
        rdcx_filter_type_t filter_type,
        bool verbose
        );

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
        rdcx_filter_type_t filter_type,
        bool verbose
        );

    /**
     * @brief Check convergence of iterative refinement procedure (simplified version)
     */
    convergence_status_t check_convergence(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        double epsilon_y,
        int iteration,
        int max_iterations
        );

    /**
     * @brief Check convergence with enhanced diagnostics (simplified version)
     */
    detailed_convergence_status_t check_convergence_detailed(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        double epsilon_y,
        double epsilon_rho,
        int iteration,
        int max_iterations,
        const std::vector<double>& response_change_history,
        const gcv_history_t& gcv_history
        );

    /**
     * @brief Check convergence with full geometric tracking
     */
    convergence_status_t check_convergence_with_geometry(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces_prev,
        const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces_curr,
        double epsilon_y,
        double epsilon_rho,
        int iteration,
        int max_iterations
        );

    // old convergence helper functions with rho
    convergence_status_t check_convergence_with_rho(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        const std::vector<vec_t>& rho_prev,
        const std::vector<vec_t>& rho_curr,
        double epsilon_y,
        double epsilon_rho,
        int iteration,
        int max_iterations
        );

    detailed_convergence_status_t check_convergence_with_rho_detailed(
        const vec_t& y_hat_prev,
        const vec_t& y_hat_curr,
        const std::vector<vec_t>& rho_prev,
        const std::vector<vec_t>& rho_curr,
        double epsilon_y,
        double epsilon_rho,
        int iteration,
        int max_iterations,
        const std::vector<double>& response_change_history = {}
        );

    // ================================================================
    // METHODS
    // ================================================================

    /**
     * @brief Initialize Riemannian simplicial complex from k-NN structure
     */
    void initialize_from_knn(
        const spmat_t& X,
        index_t k,
        bool use_counting_measure,
        double density_normalization,
        double max_ratio_threshold,
        double threshold_percentile,
        double density_alpha,
        double density_epsilon
        );

    /**
     * @brief Compute number of connected components in graph
     * @return Number of connected components (1 = connected)
     */
    int compute_connected_components();

    /**
     * @brief Assemble all Laplacian operators from current metric
     */
    void assemble_operators();

    /**
     * @brief Extend the complex by one dimension, preserving all existing structure
     */
    void extend_by_one_dim(index_t n_new);

    /**
     * @brief Batch computation of hop-extremp radii
     */
    std::vector<size_t> compute_hop_extremp_radii_batch(
        const std::vector<size_t>& vertices,
        const vec_t& y,
        double p_threshold = 0.90,
        bool detect_maxima = true,
        size_t max_hop = 20
        ) const;

private:

    std::vector<std::unordered_set<index_t>> neighbor_sets;

    // ================================================================
    // INTERNAL HELPERS
    // ================================================================

    /**
     * @brief Compute effective degrees for all vertices
     *
     * Uses initial edge densities from vertex_cofaces to compute connectivity strength.
     */
    vec_t compute_effective_degrees() const;

    /**
     * @brief Compute hop-extremp radius for a single vertex
     *
     * @param vertex Vertex index
     * @param y Function values (typically from sig.y_hat_hist.back())
     * @param p_threshold Extremp threshold (default 0.90)
     * @param detect_maxima True for maxp, false for minp
     * @param max_hop Maximum hop distance to explore
     */
    size_t compute_hop_extremp_radius(
        size_t vertex,
        const vec_t& y,
        double p_threshold = 0.90,
        bool detect_maxima = true,
        size_t max_hop = 20
        ) const;

    /**
     * @brief Compute and cache spectral decomposition of vertex Laplacian
     */
    void compute_spectral_decomposition(
        int n_eigenpairs,
        bool verbose
        );

    /**
     * @brief Automatically select diffusion and damping parameters
     */
    void select_diffusion_parameters(
        double& t_diffusion,
        double& beta_damping,
        double t_scale_factor,
        double beta_coefficient_factor,
        int n_eigenpairs,
        bool verbose
        );

    /**
     * @brief Initialize reference measure for density computation
     */
    void initialize_reference_measure(
        const std::vector<std::vector<index_t>>& knn_neighbors,
        const std::vector<std::vector<double>>& knn_distances,
        bool use_counting_measure,
        double density_normalization,
        double alpha,
        double epsilon
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
     * @brief Update vertex mass matrix from evolved vertex densities
     */
    void update_vertex_metric_from_density();

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

    std::string format_gcv_history(
        const gcv_history_t& gcv_history,
        int max_display
        ) const;
};
