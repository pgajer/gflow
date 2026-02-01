/**
 * @file riem_dcx_initialization.cpp
 *
 * @brief Initialization routines for Riemannian simplicial complex regression
 *
 * This file implements the geometric initialization phase that constructs
 * the simplicial complex from k-NN data, computes initial densities, and
 * assembles the initial Hodge Laplacian operators.
 *
 * DEBUG INSTRUMENTATION
 * =====================
 * This file contains debug serialization code controlled by DEBUG_INITIALIZE_FROM_KNN.
 *
 * To enable debugging:
 * 1. Set DEBUG_INITIALIZE_FROM_KNN to 1 in this file
 * 2. Create directory: mkdir -p /tmp/gflow_debug/initialize_from_knn
 * 3. Recompile and run
 * 4. Debug files will be written to /tmp/gflow_debug/initialize_from_knn/
 *
 * Debug files include:
 * - phase_1a_knn_result.bin: k-NN computation results
 * - phase_1b_edges_pre_pruning.bin: Edge list before geometric pruning
 * - phase_1c_edges_post_pruning.bin: Edge list after geometric pruning
 * - phase_1e_connectivity.bin: Final connectivity information
 *
 * Use R/debug_comparison.R::compare_debug_outputs() to compare with create_iknn_graph
 *
 * IMPORTANT: Debug code must remain disabled (flag = 0) for CRAN submission.
 */

// ============================================================
// DEBUG CONTROL
// ============================================================
// Creates files in /tmp/gflow_debug/initialize_from_knn/
#define DEBUG_INITIALIZE_FROM_KNN 0

#include "riem_dcx.hpp"
#include "set_wgraph.hpp"
#include "kNN.h"
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <numeric>
#include <cmath>

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
    );

static inline double median_of_copy(std::vector<double> x, bool average_even = false) {
    const size_t n = x.size();
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();

    auto mid = x.begin() + (n / 2);
    std::nth_element(x.begin(), mid, x.end());
    double m_hi = *mid;

    if (!average_even || (n % 2 == 1)) {
        return m_hi; // odd n -> exact median; even n -> upper median
    }

    auto mid_lo = x.begin() + (n / 2 - 1);
    std::nth_element(x.begin(), mid_lo, x.end());
    double m_lo = *mid_lo;

    return 0.5 * (m_lo + m_hi);
}

static inline double positive_floor(double scale_hint) {
    const double base = 1e-12;
    if (!std::isfinite(scale_hint) || scale_hint <= 0.0) scale_hint = 1.0;
    return base * std::max(1.0, scale_hint);
}

void validate_and_log_compress_weights_in_place(
    std::vector<double>& vertex_weights,
    bool average_even_median = false,
    bool clamp_nonpositive_to_zero = false
) {
    const size_t n = vertex_weights.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(vertex_weights[i])) {
            Rf_error("vertex_weights contains non-finite values.");
        }
        if (vertex_weights[i] < 0.0) {
            if (clamp_nonpositive_to_zero) vertex_weights[i] = 0.0;
            else Rf_error("vertex_weights contains negative values.");
        }
    }

    double w_median = median_of_copy(vertex_weights, average_even_median);

    if (!std::isfinite(w_median) || w_median <= 0.0) {
        w_median = positive_floor(1.0);
    } else {
        w_median = std::max(w_median, positive_floor(w_median));
    }

    for (size_t i = 0; i < n; ++i) {
        vertex_weights[i] = std::log1p(vertex_weights[i] / w_median);
    }
}

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

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
);

/**
 * @brief Initialize geometric complex from k-nearest neighbor structure
 *
 * Constructs the 1-skeleton simplicial complex from feature matrix X by building
 * a k-nearest neighbor graph and initializing all geometric structures required
 * for the regression framework. This function establishes the topological foundation,
 * validates graph connectivity, initializes density measures, and prepares the
 * geometric operators for subsequent iterative refinement.
 *
 * The initialization proceeds through five phases:
 *
 * **Phase 1A: k-NN Construction** computes the k nearest neighbors for each point
 * using Euclidean distance in the feature space. The neighborhood information
 * defines local covering sets that will determine edge connectivity.
 *
 * **Phase 1B: Edge Construction** forms edges between vertices i and j when their
 * k-neighborhoods overlap (share at least one common neighbor). For each edge,
 * the function records the intersection size and computes the minimum path length
 * through common neighbors as an initial edge weight estimate.
 *
 * **Connectivity Validation** verifies that the resulting graph forms a single
 * connected component. Disconnected graphs cannot support the geometric diffusion
 * and spectral methods underlying the regression framework. If connectivity fails,
 * the function terminates with an error message suggesting parameter adjustments.
 *
 * **Phase 1C: Geometric Edge Pruning** (optional) removes edges with anomalously
 * large weights relative to the local weight distribution. This pruning step
 * eliminates spurious long-range connections that could distort the geometry,
 * typically arising from the k-NN construction in regions of varying density.
 *
 * **Phase 2: Density Initialization** constructs the reference measure μ and
 * computes initial vertex and edge densities by aggregating μ over appropriate
 * neighborhoods. The reference measure can be uniform (counting measure) or
 * density-weighted based on k-nearest neighbor distances.
 *
 * **Phase 3: Metric Construction** builds the initial Riemannian metric from
 * the initialized densities, constructing mass matrices M_0 (vertices) and M_1
 * (edges) that encode the geometric structure.
 *
 * **Phase 4: Build Boundary Operators** constructs the discrete differential
 * operators B_1 and B_2 encoding the simplicial complex topology.
 *
 * **Phase 5: Laplacian Assembly** assembles the Hodge Laplacian operators L_0
 * and L_1 from the metric and boundary operators, preparing them for use in
 * spectral filtering and density evolution.
 *
 * @param X Feature matrix (n × d) in sparse format, where n is the number of
 *          observations and d is the number of features. Each row represents
 *          one observation in feature space. Features should be appropriately
 *          scaled if they have different units or ranges.
 *
 * @param k Number of nearest neighbors to use for graph construction. Must
 *          satisfy 2 ≤ k < n. Larger k creates more edges and makes connectivity
 *          more likely but may include geometrically irrelevant long-range
 *          connections. Typical values: 10-50. The minimum k for connectivity
 *          depends on data geometry but is often around log(n) + 5.
 *
 * @param use_counting_measure If true, use uniform reference measure μ(x) = 1
 *                             for all vertices (counting measure). If false,
 *                             use distance-based reference measure
 *                             μ(x) = (ε + d_k(x))^(-α) where d_k(x) is the
 *                             distance to the k-th nearest neighbor. Counting
 *                             measure is appropriate for uniformly sampled data;
 *                             distance-based measure adapts to varying sampling
 *                             density.
 *
 * @param density_normalization Target sum for normalized density measures.
 *                              If 0 (default), densities are normalized to sum
 *                              to n, giving average vertex density = 1. If positive,
 *                              densities sum to the specified value. Most users
 *                              should use the default (0 for automatic n-normalization).
 *
 * @param max_ratio_threshold Maximum allowed ratio between edge weight and local
 *                           median edge weight during geometric pruning. Edges
 *                           with weights exceeding this ratio times the local
 *                           median are removed as anomalous long-range connections.
 *                           Set to Inf to disable pruning. Typical value: 3-10.
 *                           Only used when pruning is enabled.
 *
 * @param threshold_percentile Percentile (in [0,1]) of local edge weight
 *                            distribution used as reference for pruning threshold.
 *                            Value of 0.5 uses local median; 0.75 uses 75th
 *                            percentile. Higher percentiles create more aggressive
 *                            pruning. Only used when pruning is enabled.
 *
 * @param density_alpha Exponent α in [1,2] for distance-based reference measure
 *                     formula μ(x) = (ε + d_k(x))^(-α). Larger values (near 2)
 *                     create stronger density-adaptive weighting, concentrating
 *                     measure on densely sampled regions. Smaller values (near 1)
 *                     produce more uniform weighting. Only used when
 *                     use_counting_measure = false. Default: 1.5.
 *
 * @param density_epsilon Regularization parameter ε > 0 in reference measure
 *                       formula μ(x) = (ε + d_k(x))^(-α). Prevents numerical
 *                       issues when nearest neighbor distances are very small.
 *                       Should be small relative to typical k-NN distances but
 *                       large enough for stability. Only used when
 *                       use_counting_measure = false. Default: 1e-10.
 *
 * @throws Rf_error If the constructed k-NN graph is disconnected (has more than
 *                  one connected component). Error message suggests increasing k
 *                  or removing outlier samples.
 *
 * @pre The feature matrix X must have at least 2 rows (n ≥ 2)
 * @pre Parameter k must satisfy 1 ≤ k < n
 * @pre If use_counting_measure = false, density_alpha must be in [1, 2]
 * @pre If use_counting_measure = false, density_epsilon must be positive
 *
 * @post vertex_cofaces is populated with graph topology
 * @post reference_measure contains initialized vertex weights
 * @post Initial densities rho.rho[0] and rho.rho[1] are computed and normalized
 * @post Metric matrices g.M[0] and g.M[1] are constructed
 * @post Boundary operators B[1] and B[2] are assembled
 * @post Hodge Laplacians L.L[0] and L.L[1] are assembled
 * @post The object is ready for fit_rdgraph_regression() to proceed
 *       with iterative refinement
 *
 * @note To explore graph connectivity before calling this function, use
 *       create.iknn.graphs() or create.single.iknn.graph() from R to examine
 *       connectivity structure across different k values.
 *
 * @note Edge pruning (controlled by max_ratio_threshold) is optional but
 *       recommended for data with highly variable local density, where k-NN
 *       construction may create spurious long-range edges in sparse regions.
 *
 * @see fit_rdgraph_regression() for the main regression function
 * @see compute_connected_components() for connectivity validation
 * @see initialize_reference_measure() for density weighting details
 */
void riem_dcx_t::initialize_from_knn(
    const spmat_t& X,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double max_ratio_threshold,
	double path_edge_ratio_percentile,
    double threshold_percentile,
	double density_alpha,
	double density_epsilon,
	verbose_level_t verbose_level
) {
    const size_t n_points = static_cast<size_t>(X.rows());

    // Setup debug directory if enabled
    std::string debug_dir;
    #if DEBUG_INITIALIZE_FROM_KNN
        debug_dir = "/tmp/gflow_debug/initialize_from_knn/";
        // Note: In production code, you'd want to ensure this directory exists
        // For debugging, create it manually: mkdir -p /tmp/gflow_debug/initialize_from_knn/
    #endif

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

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_knn_result(
            debug_dir + "phase_1a_knn_result.bin",
            knn_indices, knn_distances, n_points, k
        );
    #endif

    neighbor_sets.resize(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_neighbor_sets(
            debug_dir + "phase_1a_neighbor_sets.bin",
            neighbor_sets
        );
    #endif

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


    #if DEBUG_INITIALIZE_FROM_KNN
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
    #endif


	// ================================================================
	// CONNECTIVITY VALIDATION
	// ================================================================

	/**
	 * Verify that the constructed k-nearest neighbor graph is connected.
	 *
	 * The graph must form a single connected component for the geometric
	 * diffusion and spectral methods to function properly. If the graph has
	 * multiple components, geometric information cannot propagate between
	 * components, and the Laplacian eigendecomposition will have multiple
	 * zero eigenvalues corresponding to the disconnected pieces.
	 *
	 * Common causes of disconnected graphs include:
	 * - k value too small relative to data geometry
	 * - Isolated outlier points far from main data cloud
	 * - Data naturally divided into separate clusters
	 *
	 * Solutions:
	 * - Increase k to create more connections
	 * - Remove outlier observations before fitting
	 * - Process each connected component separately if appropriate
	 */
	int n_conn_comps = compute_connected_components();

	if (n_conn_comps > 1) {
		Rf_error("The constructed k-NN graph (k=%d) has %d connected components.\n"
				 "The regression requires a fully connected graph.\n"
				 "Solutions:\n"
				 "  1. Increase k to create more edges\n"
				 "  2. Remove outlier samples that are isolated\n"
				 "  3. Use create.iknn.graphs() to explore connectivity across k values",
				 (int)k, n_conn_comps);
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


    #if DEBUG_INITIALIZE_FROM_KNN
	size_t n_edges_before_pruning = 0;
    for (size_t i = 0; i < n_points; ++i) {
		n_edges_before_pruning += vertex_cofaces[i].size() - 1;
	}
	n_edges_before_pruning /= 2;

	Rprintf("=== Before geometric pruning ===\n");
	Rprintf("Edges before: %zu\n", n_edges_before_pruning);
    #endif

    set_wgraph_t pruned_graph = temp_graph.prune_edges_geometrically(
        max_ratio_threshold,
		path_edge_ratio_percentile
    );

    #if DEBUG_INITIALIZE_FROM_KNN
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
    #endif

	// ---- Quantile-based edge length pruning
	if (threshold_percentile > 0.0) {
        #if DEBUG_INITIALIZE_FROM_KNN
		Rprintf("=== Applying quantile-based edge length pruning ===\n");
		Rprintf("threshold_percentile=%.3f (removing top %.1f%% of edges by length)\n",
				threshold_percentile, 100.0 * (1.0 - threshold_percentile));
        #endif

		pruned_graph = pruned_graph.prune_long_edges(threshold_percentile);

		// Recount edges after quantile pruning
		size_t n_edges_in_pruned_graph_sz = 0;
		for (const auto& nbrs : pruned_graph.adjacency_list) {
			n_edges_in_pruned_graph_sz += nbrs.size();
		}
		n_edges_in_pruned_graph_sz /= 2;

        #if DEBUG_INITIALIZE_FROM_KNN
		Rprintf("=== After quantile pruning ===\n");
		Rprintf("Edges after quantile pruning: %zu\n", n_edges_in_pruned_graph_sz);
		Rprintf("Additional edges removed: %zu\n",
				n_edges_after_geometric_sz - n_edges_in_pruned_graph_sz);
        #endif
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

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_vertex_cofaces(
            debug_dir + "phase_1d_vertex_cofaces_post_pruning.bin",
            vertex_cofaces
        );
    #endif

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

    #if DEBUG_INITIALIZE_FROM_KNN
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
    #endif

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
        density_normalization,
		density_alpha,
		density_epsilon,
		verbose_level
    );

    // Compute initial densities from reference measure
    compute_initial_densities(verbose_level);

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

/**
 * Initialize reference measure μ that provides base weights for vertices.
 *
 * The reference measure serves as the starting point for density estimation
 * before any geometric evolution. We support two approaches:
 *
 * When use_counting_measure = true, all points receive uniform weight,
 * appropriate for uniformly sampled data or when sampling density should not
 * influence the geometry.
 *
 * When use_counting_measure = false, we use distance-based weighting:
 *   μ({x_i}) = (ε + d_k(x_i))^(-α)
 * where d_k(x_i) is the distance to the k-th nearest neighbor. This formula
 * provides a robust local density surrogate without requiring bandwidth
 * selection, addressing kernel density estimation instability in moderate
 * to high dimensions.
 *
 * The exponent α controls the strength of density adaptation. Larger α values
 * (near 2) create stronger weighting toward densely sampled regions, while
 * smaller values (near 1) produce more uniform weighting. The regularization
 * parameter ε prevents numerical issues when distances are very small.
 *
 * After computing raw weights, the function normalizes them to sum to the
 * target value (typically n), ensuring numerical stability and consistent
 * interpretation across datasets of different sizes.
 *
 * @param knn_neighbors k-nearest neighbor indices for each point
 * @param knn_distances k-nearest neighbor distances for each point
 * @param use_counting_measure If true, use uniform weights; if false, use distance-based
 * @param density_normalization Target sum for normalized weights (typically n)
 * @param density_alpha Exponent in [1, 2] for density estimation.
 * @param density_epsilon Regularization to prevent division by zero
 */
void riem_dcx_t::initialize_reference_measure(
    const std::vector<std::vector<index_t>>& knn_neighbors,
    const std::vector<std::vector<double>>& knn_distances,
    bool use_counting_measure,
    double density_normalization,
    double density_alpha,
    double density_epsilon,
	verbose_level_t verbose_level
) {
    const size_t n = knn_neighbors.size();
    std::vector<double> vertex_weights(n);

    if (use_counting_measure) {
        // Counting measure: uniform weights
        std::fill(vertex_weights.begin(), vertex_weights.end(), 1.0);
    } else {
        // ================================================================
        // ROBUST DISTANCE-BASED MEASURE
        // ================================================================

		dk_raw.clear();
		dk_clamped.clear();
		dk_clamped_low.clear();
		dk_clamped_high.clear();
		dk_lower = NA_REAL;
		dk_upper = NA_REAL;

        // Step 1: Collect all k-th neighbor distances
        std::vector<double> all_dk(n);
        for (size_t i = 0; i < n; ++i) {
            all_dk[i] = knn_distances[i].back();  // Last neighbor distance
        }
		dk_raw = all_dk;
		dk_clamped = all_dk;

        // Step 2: Compute robust statistics for scaling
        std::vector<double> sorted_dk = all_dk;
        std::sort(sorted_dk.begin(), sorted_dk.end());

        double d_median = sorted_dk[n / 2];
        // double d_q1 = sorted_dk[n / 4];
        // double d_q3 = sorted_dk[(3 * n) / 4];
        // double d_iqr = d_q3 - d_q1;

        // Step 3: Define reasonable bounds to prevent extreme ratios
        // Allow at most 100:1 ratio in raw d_k values
        double d_min_allowed = d_median / 10.0;
        double d_max_allowed = d_median * 10.0;
		dk_lower = d_min_allowed;
		dk_upper = d_max_allowed;

		if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
			Rprintf("\n\tReference measure: d_k range [%.3e, %.3e], median=%.3e\n",
					sorted_dk[0], sorted_dk[n-1], d_median);
		}

        int n_clamped_low = 0;
        int n_clamped_high = 0;

        // Step 4: Compute weights with bounded d_k
        for (size_t i = 0; i < n; ++i) {
            double d_k = all_dk[i];

            // Clamp to reasonable bounds
            if (d_k < d_min_allowed) {
                d_k = d_min_allowed;
				dk_clamped_low.push_back((index_t)i);
                n_clamped_low++;
            }
            if (d_k > d_max_allowed) {
                d_k = d_max_allowed;
				dk_clamped_high.push_back((index_t)i);
                n_clamped_high++;
            }

			dk_clamped[i] = d_k;

            // Apply power law with regularization
            // w(x) = (ε + d_k)^(-α)
            vertex_weights[i] = std::pow(density_epsilon + d_k, -density_alpha);
        }

        if (n_clamped_low > 0 || n_clamped_high > 0) {
			if (vl_at_least(verbose_level, verbose_level_t::TRACE))
            Rprintf("\tClamped %d vertices to lower bound, %d to upper bound (prevents extreme weights)\n",
                    n_clamped_low, n_clamped_high);
        }

        // Step 5: Additional safety - limit final weight ratio to 1000:1
        double w_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
        double w_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());
        double initial_ratio = w_max / w_min;

        if (initial_ratio > 1000.0) {
			if (vl_at_least(verbose_level, verbose_level_t::TRACE))
				Rprintf("\tInitial weight ratio %.2e exceeds 1000:1. Applying compression...\n",
						initial_ratio);

            // Use log-compression to reduce dynamic range
            // w_compressed = log(1 + w/w_median)
			validate_and_log_compress_weights_in_place(vertex_weights);

            w_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
            w_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());
            double compressed_ratio = w_max / w_min;

			if (vl_at_least(verbose_level, verbose_level_t::TRACE))
				Rprintf("\tAfter compression: weight ratio = %.2e\n", compressed_ratio);
        }
    }

    // Normalize to target sum
    double total = std::accumulate(vertex_weights.begin(), vertex_weights.end(), 0.0);
    double target = (density_normalization > 0.0) ?
                    density_normalization : static_cast<double>(n);

    if (total > 1e-15) {
        double scale = target / total;
        for (double& w : vertex_weights) {
            w *= scale;
        }
    }

    // Final safety check
    double final_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
    double final_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());
    double final_ratio = final_max / final_min;

	if (vl_at_least(verbose_level, verbose_level_t::TRACE))
		Rprintf("\tFinal reference measure: range [%.3e, %.3e], ratio=%.2e\n",
				final_min, final_max, final_ratio);

    if (final_ratio > 1e6) {
        Rf_warning("Reference measure has extreme ratio %.2e. "
                   "Consider using use_counting_measure=TRUE for more uniform initialization.",
                   final_ratio);
    }

    this->reference_measure = std::move(vertex_weights);
}
