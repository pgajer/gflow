#include "riem_dcx.hpp"
#include "kNN_r.h"            // for S_kNN()

#include <unordered_set>
#include <numeric>    // For std::accumulate

#include <ANN/ANN.h>  // ANN library header

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/**
 * @brief Initialize reference measure for density computation
 *
 * Constructs the reference measure μ that assigns base weights to vertices
 * before any geometric evolution. The reference measure provides the starting
 * point for density estimation.
 *
 * Two approaches serve different data characteristics:
 * - Counting measure treats all points equally (appropriate for uniform sampling)
 * - Distance-based measure down-weights points in dense regions (appropriate
 *   when sampling density varies)
 *
 * The formula μ({x}) = (ε + d_k(x))^(-α) with α ∈ [1,2] provides a robust
 * local density surrogate without bandwidth selection, addressing kernel density
 * estimation instability in moderate-to-high dimensions.
 *
 * @param knn_neighbors k-nearest neighbor indices for each point
 * @param knn_distances k-nearest neighbor distances for each point
 * @param use_counting_measure If true, use uniform weights; if false, use distance-based
 * @param density_normalization Target sum for normalized weights (typically n)
 */
void riem_dcx_t::initialize_reference_measure(
    const std::vector<std::vector<index_t>>& knn_neighbors,
    const std::vector<std::vector<double>>& knn_distances,
    bool use_counting_measure,
    double density_normalization
) {
    const size_t n = knn_neighbors.size();
    std::vector<double> vertex_weights(n);

    if (use_counting_measure) {
		// Counting measure: each vertex has unit weight
        std::fill(vertex_weights.begin(), vertex_weights.end(), 1.0);
    } else {
		// Distance-based measure using d_k distances
        // Formula: w(x) = (ε + d_k(x))^(-α)
        // where d_k(x) is the distance to the kth nearest neighbor

		const double epsilon = 1e-10;  // Regularization to prevent division by zero
        const double alpha = 1.5;      // Exponent in [1, 2] for density estimation

        // For each vertex, extract distance to kth (last) neighbor
        for (size_t i = 0; i < n; ++i) {
            // The last element in knn_distances[i] is d_k(x_i)
            const double d_k = knn_distances[i].back();

            // Compute weight: points in dense regions (small d_k) get higher weights
            vertex_weights[i] = std::pow(epsilon + d_k, -alpha);
        }
    }

    // Normalize to target sum
    double total = std::accumulate(vertex_weights.begin(),
                                   vertex_weights.end(), 0.0);
    double target = (density_normalization > 0.0) ?
                    density_normalization : static_cast<double>(n);
    if (total > 1e-15) {
        double scale = target / total;
        for (double& w : vertex_weights) w *= scale;
    }

    // Store in member variable for use by other modules
    this->reference_measure = std::move(vertex_weights);
}

/**
 * @brief Extend the complex by one dimension, preserving all existing structure
 *
 * @details
 * Appends a new chain degree to the discrete chain complex, incrementing pmax
 * by one and initializing all data structures for the new dimension. This
 * operation preserves all existing metric structures, boundary operators, and
 * Laplacians at lower dimensions. The function is designed for incremental
 * complex construction where simplices are built dimension by dimension.
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
 *    - g.M (metric diagonal and matrices)
 *    - g.M_solver (Cholesky factorizations, if maintained)
 *    - L.B (boundary operators)
 *    - L.L (Hodge Laplacians)
 *
 * 3. Initializes the new dimension p_new with default values:
 *    - g.M[p_new]: Sparse diagonal matrix with ones on diagonal
 *    - g.M_solver[p_new]: Null pointer (no factorization yet)
 *    - L.B[p_new]: Empty boundary operator matrix (shape: n_{p-1} × n_new)
 *    - L.L[p_new]: Empty Laplacian (to be assembled after boundary operator is built)
 *
 * USAGE PATTERN:
 * The typical usage in complex construction is:
 *
 * @code
 * // Start with vertices only (dimension 0)
 * riem_dcx_t dcx;
 * dcx.pmax = 0;
 * dcx.vertex_cofaces.resize(n_vertices);
 * // ... populate vertex_cofaces with self-loops and edges ...
 *
 * // Add edge dimension
 * dcx.extend_by_one_dim(n_edges);  // Now pmax = 1
 * dcx.edge_registry.resize(n_edges);
 * // ... populate edge_registry ...
 * dcx.build_boundary_operator_from_edges();
 *
 * // Add triangle dimension
 * dcx.extend_by_one_dim(n_triangles);  // Now pmax = 2
 * dcx.edge_cofaces.resize(n_edges);
 * // ... populate edge_cofaces with triangles ...
 * dcx.build_boundary_operator_from_triangles();
 *
 * // Compute all Laplacians
 * dcx.assemble_operators();
 * @endcode
 *
 * METRIC INITIALIZATION:
 * The new dimension starts with identity metric (diagonal ones). This is a
 * safe default that will be overwritten by the caller based on geometric
 * measurements:
 *   - For vertices: M₀ = diag(ρ₀) set by initialize_metric_from_density()
 *   - For edges: M₁ computed by compute_edge_mass_matrix()
 *   - For triangles: M₂ would follow similar pattern if implemented
 *
 * BOUNDARY OPERATOR:
 * The boundary map B[p_new]: C_{p_new} → C_{p_new-1} is allocated with the
 * correct dimensions but contains no entries initially. The caller must
 * populate this by calling the appropriate boundary construction function:
 *   - B[1]: build_boundary_operator_from_edges()
 *   - B[2]: build_boundary_operator_from_triangles()
 *
 * The boundary operator is built from coface structures (vertex_cofaces,
 * edge_cofaces) rather than simplex tables, making this function independent
 * of the deprecated S and stars structures.
 *
 * DIMENSION COUNTING:
 * Unlike the old version that used S[p-1].size() to determine the number of
 * (p-1)-simplices, this version determines dimensions from coface structures:
 *   - n_0 (vertices) = vertex_cofaces.size()
 *   - n_1 (edges) = edge_registry.size()
 *   - n_2 (triangles) = maximum triangle index in edge_cofaces + 1
 *
 * LAPLACIAN ASSEMBLY:
 * The Hodge Laplacian L[p_new] is left empty and must be computed via
 * assemble_operators() after:
 *   1. All boundary operators B[p_new+1] and B[p_new] are populated
 *   2. The metric g.M[p_new] reflects the actual geometry
 *
 * The Laplacian at dimension p is given by the Hodge decomposition:
 *     L[p] = B[p+1]^T M[p+1]^{-1} B[p+1] + M[p]^{-1} B[p] M[p-1] B[p]^T
 *
 * INVARIANTS MAINTAINED:
 * After calling extend_by_one_dim(n_new):
 *   - pmax has increased by 1
 *   - g.M.size() == L.B.size() == L.L.size() == pmax + 1
 *   - g.M[pmax] is diagonal with n_new ones
 *   - L.B[pmax] is empty sparse matrix of shape (n_{pmax-1}, n_new)
 *   - All structures at dimensions 0, ..., pmax-1 are unchanged
 *
 * EFFICIENCY:
 * This operation is O(n_new) for initializing the new dimension's data
 * structures, plus O(pmax) for resizing dimension-indexed containers.
 * Critically, it does NOT rebuild or copy data from lower dimensions,
 * making it suitable for use during complex construction.
 *
 * THREAD SAFETY:
 * This function modifies the riem_dcx_t object and is NOT thread-safe. Do not
 * call concurrently from multiple threads on the same object.
 *
 * @param n_new Number of simplices at the new dimension p_new = pmax + 1.
 *              Must be non-negative. If zero, the dimension is added but
 *              contains no simplices (useful for maintaining consistent
 *              dimension indexing).
 *
 * @post pmax is incremented by 1
 * @post All dimension-indexed containers have size pmax + 1
 * @post g.M[pmax] is diagonal identity matrix of size n_new × n_new
 * @post L.B[pmax] is empty sparse matrix of appropriate shape
 * @post L.L[pmax] is empty (to be assembled later)
 * @post All data at dimensions 0, ..., pmax-1 (pre-increment) are unchanged
 *
 * @note The caller must subsequently:
 *       1. Populate appropriate coface structures (edge_registry, edge_cofaces, etc.)
 *       2. Set g.M[pmax] entries based on geometric measurements (via metric functions)
 *       3. Build L.B[pmax] by calling boundary construction functions
 *       4. Call assemble_operators() to compute L.L[pmax]
 *
 * @note If n_new = 0, the dimension exists but is empty. This is valid when
 *       no simplices of a given dimension exist (e.g., no triangles in a
 *       sparse graph).
 *
 * @see build_boundary_operator_from_edges() - Builds B[1] from edge_registry
 * @see build_boundary_operator_from_triangles() - Builds B[2] from edge_cofaces
 * @see assemble_operators() - Computes Laplacians after boundary maps are built
 * @see initialize_metric_from_density() - Sets metric from densities in cofaces
 */
void riem_dcx_t::extend_by_one_dim(index_t n_new) {
	const int p_new = pmax + 1;

	// Increment maximum dimension
	pmax = p_new;

	// Resize metric family
	g.M.resize(pmax + 1);
	g.M_solver.resize(pmax + 1);

	// Resize operator family
	L.B.resize(pmax + 1);
	L.L.resize(pmax + 1);

	// Initialize metric at new dimension (identity)
	g.M[p_new] = spmat_t(static_cast<Eigen::Index>(n_new),
						 static_cast<Eigen::Index>(n_new));
	g.M[p_new].reserve(Eigen::VectorXi::Constant(static_cast<int>(n_new), 1));
	for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n_new); ++i) {
		g.M[p_new].insert(i, i) = 1.0;
	}
	g.M[p_new].makeCompressed();
	g.M_solver[p_new].reset();

	// Initialize boundary operator (empty matrix with correct shape)
	// Determine number of (p-1)-simplices from coface structures
	Eigen::Index n_prev = 0;
	if (p_new == 1) {
		// Boundary of edges: n_0 = number of vertices
		n_prev = static_cast<Eigen::Index>(vertex_cofaces.size());
	} else if (p_new == 2) {
		// Boundary of triangles: n_1 = number of edges
		n_prev = static_cast<Eigen::Index>(edge_registry.size());
	} else if (p_new >= 3) {
		// For higher dimensions, would need additional coface structures
		// For now, this is not implemented as we only handle up to dimension 2
		n_prev = 0;
	}

	const Eigen::Index cols = static_cast<Eigen::Index>(n_new);
	L.B[p_new] = spmat_t(n_prev, cols);

	// Initialize Laplacian (empty, to be assembled later)
	L.L[p_new] = spmat_t();
}
