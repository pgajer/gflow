/**
 * @file gfc_flow.hpp
 * @brief Trajectory-based Gradient Flow Complex (GFC) computation
 *
 * This header defines data structures and functions for computing gradient
 * flow complexes using a trajectory-first approach. Rather than growing
 * basins outward from extrema (as in compute_gfc), this approach traces
 * gradient flow trajectories and derives basins as byproducts.
 *
 * The trajectory-first approach naturally covers all vertices without
 * requiring a separate expansion step, and treats gradient flow lines
 * as the fundamental geometric primitives.
 *
 * Modulation options allow the gradient direction to be influenced by
 * local density and/or edge length, enabling more geometrically-aware
 * flow computation.
 */

#ifndef GFC_FLOW_HPP
#define GFC_FLOW_HPP

#include "gfc.hpp"              // For gfc_params_t, basin_compact_t, etc.
#include "gflow_modulation.hpp" // For gflow_modulation_t

#include <vector>
#include <map>
#include <set>
#include <cstddef>
#include <utility>
#include <string>

// Forward declaration
struct set_wgraph_t;

using std::size_t;

// ============================================================================
// Trajectory Data Structures
// ============================================================================

/**
 * @brief A single gradient flow trajectory
 *
 * Represents a path through the graph following gradient flow from
 * a local minimum to a local maximum (or from an arbitrary starting
 * point to both extrema via joining).
 */
struct gflow_trajectory_t {
    std::vector<size_t> vertices;  ///< Ordered vertices along trajectory (0-based)
    size_t start_vertex;           ///< Starting vertex (local min for full trajectories)
    size_t end_vertex;             ///< Ending vertex (local max for full trajectories)
    bool starts_at_lmin;           ///< True if trajectory starts at local minimum
    bool ends_at_lmax;             ///< True if trajectory ends at local maximum
    double total_change;           ///< Total function change: y[end] - y[start]
    int trajectory_id;             ///< Unique identifier for this trajectory

    gflow_trajectory_t()
        : start_vertex(0), end_vertex(0),
          starts_at_lmin(false), ends_at_lmax(false),
          total_change(0.0), trajectory_id(-1) {}
};

// ============================================================================
// Parameter Structure
// ============================================================================

/**
 * @brief Parameters for trajectory-based GFC computation
 *
 * Extends gfc_params_t with trajectory-specific parameters including
 * modulation type and optional density weights.
 */
struct gfc_flow_params_t : public gfc_params_t {
    /// Modulation strategy for gradient direction selection
    gflow_modulation_t modulation = gflow_modulation_t::NONE;

    /// Whether to store full trajectory information in result
    bool store_trajectories = true;

    /// Maximum trajectory length (vertices) before giving up
    /// Prevents infinite loops in degenerate cases
    size_t max_trajectory_length = 10000;

    gfc_flow_params_t() = default;

    /// Construct from base parameters
    explicit gfc_flow_params_t(const gfc_params_t& base)
        : gfc_params_t(base),
          modulation(gflow_modulation_t::NONE),
          store_trajectories(true),
          max_trajectory_length(10000) {}
};

// ============================================================================
// Result Structure
// ============================================================================

/**
 * @brief Complete result from trajectory-based GFC computation
 *
 * Extends gfc_result_t conceptually with trajectory information and
 * additional statistics about the flow computation.
 */
struct gfc_flow_result_t {
    // ---- Basins (same structure as gfc_result_t) ----
    std::vector<basin_compact_t> max_basins;  ///< Maximum basins after refinement
    std::vector<basin_compact_t> min_basins;  ///< Minimum basins after refinement

    // Summary statistics for each retained extremum
    std::vector<extremum_summary_t> max_summaries;
    std::vector<extremum_summary_t> min_summaries;

    // Vertex-to-basin membership
    std::vector<std::vector<int>> max_membership;
    std::vector<std::vector<int>> min_membership;

    // Basin assignments (always computed in flow approach, no expansion needed)
    std::vector<int> max_assignment;  ///< max_assignment[v] = max basin index for vertex v
    std::vector<int> min_assignment;  ///< min_assignment[v] = min basin index for vertex v

    // Stage history for reporting
    std::vector<stage_counts_t> stage_history;

    // ---- Trajectory-specific outputs ----

    /// All computed trajectories (if store_trajectories = true)
    std::vector<gflow_trajectory_t> trajectories;

    /// Vertex to trajectory assignment: vertex_trajectory[v] = trajectory index
    /// containing vertex v (first trajectory if multiple)
    std::vector<int> vertex_trajectory;

    /// Number of trajectories started from local minima
    int n_lmin_trajectories;

    /// Number of trajectories started from non-extremal vertices (joins)
    int n_join_trajectories;

    // ---- Metadata ----
    size_t n_vertices;
    double y_median;

    // Parameters used (for reproducibility)
    gfc_flow_params_t params;

    gfc_flow_result_t()
        : n_lmin_trajectories(0), n_join_trajectories(0),
          n_vertices(0), y_median(0.0) {}
};

// ============================================================================
// Main Computation Functions (Free Functions)
// ============================================================================

/**
 * @brief Compute gradient flow complex using trajectory-first approach
 *
 * This function computes basins of attraction by tracing gradient flow
 * trajectories rather than growing basins outward from extrema. The
 * algorithm proceeds as follows:
 *
 * 1. Identify all local minima using neighborhood comparison
 * 2. From each local minimum, trace an ascending trajectory following
 *    the steepest gradient (with optional modulation) until reaching
 *    a local maximum
 * 3. Assign all trajectory vertices to both the starting min-basin
 *    and ending max-basin
 * 4. For any unvisited vertices, trace both ascending and descending
 *    trajectories, joining them at that vertex
 * 5. Apply the same filtering/merging pipeline as compute_gfc():
 *    - Relative value filtering
 *    - Overlap-based clustering
 *    - Geometric filtering
 *
 * The trajectory-first approach naturally covers all vertices without
 * requiring a separate expansion step.
 *
 * @param graph The weighted graph structure
 * @param y Function values at each vertex
 * @param params Flow computation and refinement parameters
 * @param density Optional vertex density weights for modulation
 *                (required if modulation uses DENSITY)
 * @param verbose Print progress messages
 * @return Complete GFC flow result including trajectories
 */
gfc_flow_result_t compute_gfc_flow(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density = {},
    bool verbose = false
);

/**
 * @brief Compute GFC flow for multiple functions over the same graph
 *
 * Efficiently computes trajectory-based GFC for each column of a matrix,
 * reusing graph structure across all computations. Supports OpenMP
 * parallelization over columns.
 *
 * @param graph The weighted graph structure
 * @param Y Matrix of function values (n_vertices x n_functions)
 * @param params Flow computation parameters (applied to all functions)
 * @param density Optional vertex density weights
 * @param n_cores Number of OpenMP threads (1 = sequential)
 * @param verbose Print progress messages
 * @return Vector of GFC flow results, one per function
 */
std::vector<gfc_flow_result_t> compute_gfc_flow_matrix(
    const set_wgraph_t& graph,
    const Eigen::MatrixXd& Y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density = {},
    int n_cores = 1,
    bool verbose = false
);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Compute modulated gradient score for an edge
 *
 * Computes the score used to select the "steepest" direction based
 * on the modulation strategy.
 *
 * @param delta_y Function difference: y[target] - y[source]
 * @param edge_length Length of the edge
 * @param target_density Density at target vertex (used if modulation includes DENSITY)
 * @param modulation Modulation strategy
 * @return Computed score (higher = steeper ascent)
 */
inline double compute_modulated_score(
    double delta_y,
    double edge_length,
    double target_density,
    gflow_modulation_t modulation
) {
    // Avoid division by zero
    const double eps = 1e-12;

    switch (modulation) {
        case gflow_modulation_t::NONE:
            return delta_y;

        case gflow_modulation_t::DENSITY:
            return target_density * delta_y;

        case gflow_modulation_t::EDGELEN:
            return delta_y / (edge_length + eps);

        case gflow_modulation_t::DENSITY_EDGELEN:
            return (target_density * delta_y) / (edge_length + eps);

        default:
            return delta_y;
    }
}

#endif // GFC_FLOW_HPP
