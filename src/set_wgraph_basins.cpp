#include <queue>
#include "set_wgraph.hpp"
#include "error_utils.h"
#include <R_ext/Print.h>

/**
 * @brief Recursively enumerate all paths from a starting vertex to the origin extremum.
 *
 * This function uses the all_predecessors map to recursively build all possible paths
 * from start_vertex back to the origin of the basin. Due to convergence in gradient flow,
 * the number of paths can grow exponentially with basin size.
 *
 * @param start_vertex Vertex to start path enumeration from
 * @param origin_vertex The basin origin (extremum)
 * @param all_predecessors Map of vertex to all its valid predecessors
 * @param current_path Current path being constructed (passed by value for recursion)
 * @param all_paths Output vector accumulating all discovered paths
 * @param max_paths Safety limit to prevent memory exhaustion (default 1000000)
 *
 * @note Paths are stored in order from start_vertex to origin_vertex
 */
void enumerate_all_paths(
    size_t start_vertex,
    size_t origin_vertex,
    const std::unordered_map<size_t, std::vector<size_t>>& all_predecessors,
    std::vector<size_t> current_path,
    std::vector<std::vector<size_t>>& all_paths,
    size_t max_paths = 1000000
) {
    // Safety check to prevent memory exhaustion
    if (all_paths.size() >= max_paths) {
        return;
    }

    current_path.push_back(start_vertex);

    // Base case: reached the origin
    if (start_vertex == origin_vertex) {
        all_paths.push_back(current_path);
        return;
    }

    // Recursive case: explore all predecessors
    auto it = all_predecessors.find(start_vertex);
    if (it != all_predecessors.end() && !it->second.empty()) {
        for (size_t pred : it->second) {
            enumerate_all_paths(pred, origin_vertex, all_predecessors,
                              current_path, all_paths, max_paths);
        }
    } else {
		// If no predecessors found and we haven't reached origin, path is incomplete
		// (shouldn't happen in well-formed basin, but we handle it gracefully)
		REPORT_ERROR("No predecessors found and we haven't reached origin, path is incomplete!");
	}
}

/**
 * @brief Computes the gradient basin of attraction for a given local extremum with terminal extrema
 *
 * @details This function determines the complete basin of attraction by following all gradient
 *          flow trajectories where the extremum property is maintained. The basin consists of
 *          all vertices reachable from the given extremum via monotone paths, along with
 *          predecessor information that enables reconstruction of individual gradient trajectories
 *          and identification of terminal extrema where trajectories terminate.
 *
 *          The function performs a breadth-first search starting from the given vertex and
 *          includes all vertices that can be reached via monotone paths. For a local minimum,
 *          monotonicity requires that function values strictly increase along each edge of the
 *          path (y[w] > y[prev] for each step). For a local maximum, function values must
 *          strictly decrease along each edge (y[w] < y[prev] for each step). The search
 *          continues along each path as long as the monotonicity condition is satisfied.
 *
 *          In addition to identifying basin membership, the function detects terminal extrema,
 *          which are local extrema of the opposite type reachable within the basin. For a
 *          descending basin from a local maximum M_j, terminal extrema are local minima
 *          {m_1, m_2, ...} where gradient trajectories terminate. Similarly, for an ascending
 *          basin from a local minimum m_i, terminal extrema are local maxima {M_1, M_2, ...}.
 *          These terminal extrema define the valid gradient flow cells of the form (m_i, M_j),
 *          which represent unions of all gradient trajectories connecting the extremum pair.
 *
 *          The predecessor map enables complete trajectory reconstruction. For any vertex v in
 *          the basin, following the predecessor chain backwards (v -> predecessors[v] ->
 *          predecessors[predecessors[v]] -> ...) traces the unique gradient trajectory from
 *          the origin extremum to v. This supports downstream analyses such as monotonicity
 *          testing for biomarker discovery along gradient flow cells.
 *
 *          Global extrema are processed identically to other vertices, with their basin
 *          potentially encompassing the entire graph if no monotonicity violations exist.
 *
 * Basin membership criteria:
 * - For a local minimum at vertex: A vertex u belongs to the basin if it can be reached from
 *   vertex via a path where y strictly increases along each edge
 * - For a local maximum at vertex: A vertex u belongs to the basin if it can be reached from
 *   vertex via a path where y strictly decreases along each edge
 *
 * Terminal extrema detection:
 * - For a descending basin (from maximum): A basin vertex u is a terminal minimum if all its
 *   neighbors w satisfy y[w] >= y[u]
 * - For an ascending basin (from minimum): A basin vertex u is a terminal maximum if all its
 *   neighbors w satisfy y[w] <= y[u]
 *
 * The hop_idx field represents:
 * - The maximum hop distance (graph distance) from vertex to any basin member
 * - std::numeric_limits<size_t>::max() if the basin is empty (vertex is not a 1-hop extremum)
 * - Otherwise, the actual maximum distance within the basin
 *
 * @param vertex Index of the vertex (local extremum) whose basin to compute
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, computes basin for a maximum; otherwise, for a minimum
 *
 * @return gradient_basin_t Structure containing:
 *         - vertex: The origin extremum
 *         - value: Function value at origin
 *         - is_maximum: Basin type indicator
 *         - hop_idx: Maximum hop distance in basin
 *         - hop_dist_map: All basin vertices mapped to their hop distances
 *         - predecessors: Predecessor map for trajectory reconstruction
 *         - y_nbhd_bd_map: Boundary vertices (adjacent to but outside basin) with their values
 *         - terminal_extrema: Vector of terminal extremum indices (access values via y vector)
 *
 * @note The function returns a basin with empty hop_dist_map if vertex is not a local extremum
 *       in its 1-hop neighborhood, indicated by hop_idx = std::numeric_limits<size_t>::max()
 *
 * @note Terminal extrema values should be accessed from the original y vector as y[terminal_extrema[i]]
 *       to avoid redundant storage and improve iteration performance
 */
/**
 * @brief Compute gradient basin with optional trajectory enumeration.
 *
 * Performs breadth-first search from a local extremum to identify all vertices
 * reachable via monotone paths. Terminals are identified structurally as vertices
 * with no successors in the basin (s-terminals). When requested, enumerates all
 * gradient flow trajectories from terminals to the extremum.
 *
 * @param vertex Index of the extremum vertex
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basins (from maximum), false for ascending (from minimum)
 * @param with_trajectories If true, enumerate all trajectories (computationally expensive)
 *
 * @return gradient_basin_t Complete basin structure with optional trajectories
 */
/**
 * @brief Compute gradient basin with optional trajectory enumeration.
 *
 * Performs breadth-first search from a local extremum to identify all vertices
 * reachable via monotone paths. After BFS, reconstructs all valid monotone edges
 * within the basin (not just BFS tree edges) to correctly identify gradient flow
 * structure. Terminals are identified structurally as vertices with no successors.
 *
 * @param vertex Index of the extremum vertex
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basins (from maximum), false for ascending (from minimum)
 * @param with_trajectories If true, enumerate all trajectories (computationally expensive)
 *
 * @return gradient_basin_t Complete basin structure with optional trajectories
 */
/**
 * @brief Compute gradient basin with optional trajectory enumeration using shortest paths.
 *
 * Performs breadth-first search from a local extremum to identify all vertices
 * reachable via monotone paths. Builds all_predecessors during BFS, including only
 * edges that maintain optimal hop distance (shortest monotone paths). Terminals are
 * identified structurally as vertices with no successors. This approach ensures
 * perfect trajectory coverage while maintaining computational tractability.
 *
 * @param vertex Index of the extremum vertex
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basins (from maximum), false for ascending (from minimum)
 * @param with_trajectories If true, enumerate all shortest monotone paths
 *
 * @return gradient_basin_t Complete basin structure with optional trajectories
 */
// Strategy: Full Monotone DAG with m-Terminals Only
gradient_basin_t set_wgraph_t::compute_basin_of_attraction(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima,
    bool with_trajectories
    ) const {

    gradient_basin_t result;
    result.vertex = vertex;
    result.value = y[vertex];
    result.is_maximum = detect_maxima;

    // Check if vertex is a local extremum in its 1-hop neighborhood
    bool is_local_extremum = true;
    for (const auto& edge : adjacency_list[vertex]) {
        size_t u = edge.vertex;
        if (detect_maxima) {
            if (y[u] >= y[vertex]) {
                is_local_extremum = false;
                break;
            }
        } else {
            if (y[u] <= y[vertex]) {
                is_local_extremum = false;
                break;
            }
        }
    }

    if (!is_local_extremum) {
        result.hop_idx = std::numeric_limits<size_t>::max();
        return result;
    }

    // BFS to determine basin membership ONLY
    // (hop distance used only for diagnostics, not for predecessor selection)
    std::vector<int> hop_distance(adjacency_list.size(), -1);
    std::queue<size_t> q;

    result.hop_dist_map[vertex] = 0;
    hop_distance[vertex] = 0;
    q.push(vertex);

    size_t max_hop_distance = 0;

    // Simple BFS to find all reachable vertices via monotone edges
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();

        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (hop_distance[v] != -1) {
                continue;  // Already visited
            }

            // Check monotonicity
            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) {
                continue;
            }

            // Add to basin
            hop_distance[v] = hop_distance[u] + 1;
            result.hop_dist_map[v] = static_cast<size_t>(hop_distance[v]);
            q.push(v);

            max_hop_distance = std::max(max_hop_distance,
                                       static_cast<size_t>(hop_distance[v]));
        }
    }

    result.hop_idx = max_hop_distance;

    // Now build ALL_PREDECESSORS with ALL monotone edges in the basin
    // (no hop distance restriction)
    for (const auto& [v, _] : result.hop_dist_map) {
        result.all_predecessors[v] = std::vector<size_t>();
    }

    for (const auto& [u, hop_u] : result.hop_dist_map) {
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            // Is v in the basin?
            if (result.hop_dist_map.find(v) == result.hop_dist_map.end()) {
                continue;
            }

            // Check monotonicity for edge u -> v
            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (monotone) {
                result.all_predecessors[v].push_back(u);
            }
        }
    }

    // Identify m-terminals ONLY (local minima within basin)
    std::vector<size_t> m_terminals;
    for (const auto& [u, hop] : result.hop_dist_map) {
        if (u == vertex) continue;  // Skip origin

        bool is_local_min_in_basin = true;
        for (const auto& edge : adjacency_list[u]) {
            size_t w = edge.vertex;

            // Only consider neighbors in the basin
            if (result.hop_dist_map.find(w) == result.hop_dist_map.end()) {
                continue;
            }

            if (detect_maxima) {
                if (y[w] < y[u]) {
                    is_local_min_in_basin = false;
                    break;
                }
            } else {
                if (y[w] > y[u]) {
                    is_local_min_in_basin = false;
                    break;
                }
            }
        }

        if (is_local_min_in_basin) {
            m_terminals.push_back(u);
        }
    }

    // Use ONLY m-terminals
    result.terminal_extrema = m_terminals;

    Rprintf("\n=== Basin Structure Diagnostic ===\n");
    Rprintf("Basin vertex: %zu (y=%.6f)\n", vertex, y[vertex]);
    Rprintf("Basin size: %zu vertices\n", result.hop_dist_map.size());
    Rprintf("Max hop distance: %zu\n", max_hop_distance);
    Rprintf("m-terminals (local min in basin): %zu\n", m_terminals.size());
    Rprintf("Using FULL monotone DAG (all monotone edges)\n");
    Rprintf("===================================\n\n");

    // Compute trajectories if requested
    if (with_trajectories) {
        // Set a per-terminal path limit to prevent explosion
        const size_t MAX_PATHS_PER_TERMINAL = 50000;

        // Enumerate all monotone paths from m-terminals
        for (size_t terminal : result.terminal_extrema) {
            trajectory_set_t traj_set;
            traj_set.terminal_vertex = terminal;

            std::vector<size_t> empty_path;
            enumerate_all_paths(terminal, vertex, result.all_predecessors,
                              empty_path, traj_set.trajectories);

            // Check if we hit path limit
            if (traj_set.trajectories.size() >= MAX_PATHS_PER_TERMINAL) {
                Rprintf("WARNING: Terminal %zu hit path limit (%zu paths)\n",
                       terminal, traj_set.trajectories.size());
            }

            result.trajectory_sets.push_back(traj_set);
        }

        // Coverage diagnostic
        std::unordered_set<size_t> vertices_in_trajectories;
        for (const auto& traj_set : result.trajectory_sets) {
            for (const auto& traj : traj_set.trajectories) {
                for (size_t v : traj) {
                    vertices_in_trajectories.insert(v);
                }
            }
        }

        size_t basin_size = result.hop_dist_map.size();
        size_t covered = vertices_in_trajectories.size();

        Rprintf("\n=== Trajectory Coverage Diagnostic ===\n");
        Rprintf("Basin vertex: %zu, Basin size: %zu\n", vertex, basin_size);
        Rprintf("Vertices in trajectories: %zu\n", covered);
        Rprintf("Uncovered vertices: %zu\n", basin_size - covered);

        if (covered == basin_size) {
            Rprintf("SUCCESS: Perfect trajectory coverage (100%%)!\n");
        } else {
            Rprintf("WARNING: Incomplete coverage (%.1f%%)\n",
                   100.0 * covered / basin_size);
        }

        // Path statistics
        size_t total_paths = 0;
        size_t min_paths = std::numeric_limits<size_t>::max();
        size_t max_paths = 0;
        std::vector<size_t> path_lengths;

        for (const auto& traj_set : result.trajectory_sets) {
            size_t n = traj_set.trajectories.size();
            total_paths += n;
            if (n > 0) {
                min_paths = std::min(min_paths, n);
                max_paths = std::max(max_paths, n);
            }

            // Sample path lengths
            for (const auto& traj : traj_set.trajectories) {
                path_lengths.push_back(traj.size());
            }
        }

        if (total_paths > 0) {
            // Compute path length statistics
            std::sort(path_lengths.begin(), path_lengths.end());
            size_t median_length = path_lengths[path_lengths.size() / 2];
            size_t min_length = path_lengths.front();
            size_t max_length = path_lengths.back();

            Rprintf("\nPath enumeration statistics:\n");
            Rprintf("  Total paths: %zu\n", total_paths);
            Rprintf("  Paths per terminal: min=%zu, avg=%.1f, max=%zu\n",
                   min_paths,
                   static_cast<double>(total_paths) / result.terminal_extrema.size(),
                   max_paths);
            Rprintf("  Path lengths: min=%zu, median=%zu, max=%zu\n",
                   min_length, median_length, max_length);
        }

        Rprintf("======================================\n\n");
    }

    // Collect boundary vertices
    std::unordered_set<size_t> basin_vertices;
    for (const auto& [v, _] : result.hop_dist_map) {
        basin_vertices.insert(v);
    }

    std::unordered_set<size_t> boundary_found;
    for (const auto& [v, _] : result.hop_dist_map) {
        for (const auto& edge : adjacency_list[v]) {
            size_t u = edge.vertex;
            if (basin_vertices.find(u) == basin_vertices.end() &&
                boundary_found.find(u) == boundary_found.end()) {
                boundary_found.insert(u);
                result.y_nbhd_bd_map[u] = y[u];
            }
        }
    }

    return result;
}

