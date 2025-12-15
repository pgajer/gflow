/**
 * @file gfc.hpp
 * @brief Gradient Flow Complex (GFC) computation with refinement pipeline
 *
 * This header defines data structures and functions for computing refined
 * gradient flow complexes from scalar functions on weighted graphs. The
 * refinement pipeline includes relative value filtering, overlap-based
 * clustering, geometric filtering, and basin expansion.
 *
 * The GFC computation serves as the foundation for association analysis
 * via gradient flow, enabling robust measurement of relationships between
 * functions defined on the same graph structure.
 */

#ifndef GFC_HPP
#define GFC_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <queue>
#include <cmath>

#include <Eigen/Core>

// Forward declaration
struct set_wgraph_t;

using std::size_t;

// ============================================================================
// Parameter Structures
// ============================================================================

/**
 * @brief Parameters controlling the GFC refinement pipeline
 *
 * This structure encapsulates all parameters that control the basin
 * refinement process, from initial computation through filtering,
 * clustering, and expansion stages.
 */
struct gfc_params_t {
    // Edge length threshold for basin construction (quantile)
    double edge_length_quantile_thld = 0.9;

    // Relative value filtering parameters
    bool apply_relvalue_filter = true;
    double min_rel_value_max = 1.1;   ///< Minimum relative value for maxima
    double max_rel_value_min = 0.9;   ///< Maximum relative value for minima

    // Clustering parameters
    bool apply_maxima_clustering = true;
    bool apply_minima_clustering = true;
    double max_overlap_threshold = 0.15;  ///< Overlap distance threshold for maxima
    double min_overlap_threshold = 0.15;  ///< Overlap distance threshold for minima

    // Geometric filtering parameters
    bool apply_geometric_filter = true;
    double p_mean_nbrs_dist_threshold = 0.9;  ///< Percentile threshold for mean neighbor distance (maxima only)
    double p_mean_hopk_dist_threshold = 0.9;  ///< Percentile threshold for hop-k distance
    double p_deg_threshold = 0.9;             ///< Percentile threshold for degree
    int min_basin_size = 10;                  ///< Minimum basin size to retain

    // Basin expansion
    bool expand_basins = true;  ///< Whether to expand basins to cover all vertices

    // Auxiliary parameters
    int hop_k = 2;  ///< Hop distance for summary statistics

    // Default constructor with sensible defaults
    gfc_params_t() = default;
};

// ============================================================================
// Basin Data Structures
// ============================================================================

/**
 * @brief Compact representation of a single basin for output
 *
 * This structure provides a lightweight representation of a basin
 * suitable for returning to R, containing only the essential information
 * needed for downstream analysis.
 */
struct basin_compact_t {
    size_t extremum_vertex;            ///< 0-based vertex index of extremum
    double extremum_value;             ///< Function value at extremum
    bool is_maximum;                   ///< true = maximum basin, false = minimum
    std::vector<size_t> vertices;      ///< Basin member vertices (0-based)
    std::vector<int> hop_distances;    ///< Hop distance from extremum for each vertex
    int max_hop_distance;              ///< Maximum hop distance in basin

    basin_compact_t()
        : extremum_vertex(0), extremum_value(0.0), is_maximum(true), max_hop_distance(0) {}
};

/**
 * @brief Summary statistics for a single extremum
 *
 * Contains computed statistics about each extremum and its basin,
 * used for filtering decisions and downstream reporting.
 */
struct extremum_summary_t {
    size_t vertex;              ///< 0-based vertex index
    double value;               ///< Function value at extremum
    double rel_value;           ///< Relative to median (value / median)
    bool is_maximum;            ///< Type of extremum
    int basin_size;             ///< Number of vertices in basin
    int hop_index;              ///< Maximum hop distance in basin (hop_idx)
    double p_mean_nbrs_dist;    ///< Percentile of mean distance to neighbors
    double p_mean_hopk_dist;    ///< Percentile of mean hop-k distance
    double deg_percentile;      ///< Degree percentile of extremum vertex
    int degree;                 ///< Degree of extremum vertex

    extremum_summary_t()
        : vertex(0), value(0.0), rel_value(1.0), is_maximum(true),
          basin_size(0), hop_index(0), p_mean_nbrs_dist(0.0),
          p_mean_hopk_dist(0.0), deg_percentile(0.0), degree(0) {}
};

/**
 * @brief Record of counts at each refinement stage
 */
struct stage_counts_t {
    std::string stage_name;
    int n_max_before;
    int n_max_after;
    int n_min_before;
    int n_min_after;

    stage_counts_t()
        : stage_name(""), n_max_before(0), n_max_after(0),
          n_min_before(0), n_min_after(0) {}

    stage_counts_t(const std::string& name, int max_b, int max_a, int min_b, int min_a)
        : stage_name(name), n_max_before(max_b), n_max_after(max_a),
          n_min_before(min_b), n_min_after(min_a) {}
};

// ============================================================================
// Result Structures
// ============================================================================

/**
 * @brief Complete GFC result for a single function
 *
 * Contains all output from the GFC computation including refined basins,
 * summary statistics, membership information, and stage history.
 */
struct gfc_result_t {
    // Refined basins
    std::vector<basin_compact_t> max_basins;  ///< Maximum basins after refinement
    std::vector<basin_compact_t> min_basins;  ///< Minimum basins after refinement

    // Summary statistics for each retained extremum
    std::vector<extremum_summary_t> max_summaries;
    std::vector<extremum_summary_t> min_summaries;

    // Vertex-to-basin membership (for multiplicity handling)
    // max_membership[v] = indices of max basins containing vertex v
    std::vector<std::vector<int>> max_membership;
    std::vector<std::vector<int>> min_membership;

    // Expanded basin assignments (if expand_basins = true)
    // expanded_max_assignment[v] = index of assigned max basin (-1 if none)
    std::vector<int> expanded_max_assignment;
    std::vector<int> expanded_min_assignment;

    // Stage history for reporting
    std::vector<stage_counts_t> stage_history;

    // Metadata
    size_t n_vertices;
    double y_median;

    // Parameters used (for reproducibility)
    gfc_params_t params;

    gfc_result_t() : n_vertices(0), y_median(0.0) {}
};

// ============================================================================
// Helper Function Declarations
// ============================================================================

/**
 * @brief Compute overlap distance matrix using Szymkiewicz-Simpson coefficient
 *
 * The overlap distance between basins A and B is defined as:
 *   d(A, B) = 1 - |A ∩ B| / min(|A|, |B|)
 *
 * This measure equals 0 when one basin is completely contained in the other,
 * making it ideal for detecting nested basin structures.
 *
 * @param basin_vertices Vector of vectors, each containing vertex indices for a basin
 * @return Symmetric matrix of pairwise overlap distances
 */
Eigen::MatrixXd compute_overlap_distance_matrix(
    const std::vector<std::vector<size_t>>& basin_vertices
);

/**
 * @brief Create adjacency list for threshold graph from distance matrix
 *
 * Creates a graph where basins are connected if their overlap distance
 * is strictly less than the threshold.
 *
 * @param dist_matrix Symmetric distance matrix
 * @param threshold Distance threshold (edges created where dist < threshold)
 * @return Adjacency list (0-based indices)
 */
std::vector<std::vector<int>> create_threshold_graph(
    const Eigen::MatrixXd& dist_matrix,
    double threshold
);

/**
 * @brief Find connected components of a graph
 *
 * @param adj_list Adjacency list (0-based indices)
 * @return Vector of component assignments (0-based component indices)
 */
std::vector<int> find_connected_components(
    const std::vector<std::vector<int>>& adj_list
);

/**
 * @brief Compute edge length threshold from quantile
 *
 * @param graph The weighted graph
 * @param quantile Quantile value in (0, 1]
 * @return Edge length threshold
 */
double compute_edge_length_threshold(
    const set_wgraph_t& graph,
    double quantile
);

/**
 * @brief Expand basins to cover all vertices using shortest-path assignment
 *
 * For each vertex not covered by any basin, assigns it to the nearest
 * basin based on weighted shortest path distance.
 *
 * @param graph The weighted graph
 * @param basin_vertices Vector of basin vertex sets
 * @param n_vertices Total number of vertices
 * @return Assignment vector: assignment[v] = basin index for vertex v
 */
std::vector<int> expand_basins_to_cover(
    const set_wgraph_t& graph,
    const std::vector<std::vector<size_t>>& basin_vertices,
    size_t n_vertices
);

// ============================================================================
// Main GFC Computation Functions
// ============================================================================

/**
 * @brief Compute refined gradient flow complex for a single function
 *
 * This is the main entry point for GFC computation. It performs:
 * 1. Initial basin computation using BFS from each local extremum
 * 2. Relative value filtering (optional)
 * 3. Overlap-based clustering and merging (optional)
 * 4. Geometric filtering by hop distance and degree (optional)
 * 5. Basin expansion to cover all vertices (optional)
 *
 * @param graph The weighted graph structure
 * @param y Function values at each vertex
 * @param params Refinement parameters
 * @param verbose Print progress messages
 * @return Complete GFC result
 */
gfc_result_t compute_gfc(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const gfc_params_t& params,
    bool verbose = false
);

/**
 * @brief Compute GFC for multiple functions over the same graph
 *
 * Efficiently computes GFC for each column of a matrix, reusing
 * graph structure across all computations. Supports OpenMP
 * parallelization over columns.
 *
 * @param graph The weighted graph structure
 * @param Y Matrix of function values (n_vertices x n_functions)
 * @param params Refinement parameters (applied to all functions)
 * @param n_cores Number of OpenMP threads (1 = sequential)
 * @param verbose Print progress messages
 * @return Vector of GFC results, one per function
 */
std::vector<gfc_result_t> compute_gfc_matrix(
    const set_wgraph_t& graph,
    const Eigen::MatrixXd& Y,
    const gfc_params_t& params,
    int n_cores = 1,
    bool verbose = false
);

// ============================================================================
// Internal Pipeline Functions (exposed for testing)
// ============================================================================

namespace gfc_internal {

/**
 * @brief Filter extrema by relative value
 *
 * Removes maxima with relative value < min_rel_value_max
 * and minima with relative value > max_rel_value_min.
 *
 * @param summaries Vector of extremum summaries (modified in place)
 * @param basins Vector of basins (modified in place)
 * @param min_rel_value_max Threshold for maxima
 * @param max_rel_value_min Threshold for minima
 * @param is_maximum Whether these are maximum basins
 */
void filter_by_relvalue(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double min_rel_value_max,
    double max_rel_value_min,
    bool is_maximum
);

/**
 * @brief Cluster basins by overlap and merge within clusters
 *
 * Performs single-linkage clustering based on overlap distance,
 * then merges basins within each cluster by taking the union
 * of vertices and keeping the most extreme representative.
 *
 * @param summaries Vector of extremum summaries (modified in place)
 * @param basins Vector of basins (modified in place)
 * @param overlap_threshold Distance threshold for clustering
 * @param is_maximum Whether these are maximum basins
 */
void cluster_and_merge_basins(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double overlap_threshold,
    bool is_maximum
);

/**
 * @brief Filter basins by geometric characteristics
 *
 * Removes basins whose extrema have:
 * - Mean neighbor distance percentile >= threshold (maxima only)
 * - Mean hop-k distance percentile >= threshold
 * - Degree percentile >= threshold
 * - Basin size < minimum
 *
 * @param summaries Vector of extremum summaries (modified in place)
 * @param basins Vector of basins (modified in place)
 * @param p_mean_nbrs_dist_threshold Threshold for neighbor distance percentile (maxima only)
 * @param p_mean_hopk_dist_threshold Threshold for hop-k distance percentile
 * @param p_deg_threshold Threshold for degree percentile
 * @param min_basin_size Minimum basin size
 * @param is_maximum Whether these are maximum basins
 */
void filter_by_geometry(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double p_mean_nbrs_dist_threshold,
    double p_mean_hopk_dist_threshold,
    double p_deg_threshold,
    int min_basin_size,
    bool is_maximum
);

/**
 * @brief Compute summary statistics for basins
 *
 * Computes relative values, hop distances, and geometric measures
 * for each basin.
 *
 * @param graph The weighted graph
 * @param basins Vector of basins
 * @param y Function values
 * @param y_median Median of y values
 * @param hop_k Hop distance for statistics
 * @return Vector of summary statistics
 */
std::vector<extremum_summary_t> compute_basin_summaries(
    const set_wgraph_t& graph,
    const std::vector<basin_compact_t>& basins,
    const std::vector<double>& y,
    double y_median,
    int hop_k
);

/**
 * @brief Build membership vectors from basins
 *
 * For each vertex, determines which basins contain it.
 *
 * @param basins Vector of basins
 * @param n_vertices Total number of vertices
 * @return membership[v] = vector of basin indices containing v
 */
std::vector<std::vector<int>> build_membership_vectors(
    const std::vector<basin_compact_t>& basins,
    size_t n_vertices
);

} // namespace gfc_internal

#endif // GFC_HPP
