#include <queue>
#include "set_wgraph.hpp"
#include "error_utils.h"
#include <R_ext/Print.h>

/**
 * @brief Enumerate paths with early stopping after K paths found
 *
 * Uses depth-first recursive enumeration but stops after finding K complete paths.
 * Simpler and more efficient than BFS over path space for finding K paths.
 */
void enumerate_k_paths_recursive(
    size_t current,
    size_t target,
    const std::unordered_map<size_t, std::vector<size_t>>& all_predecessors,
    std::vector<size_t>& current_path,
    std::vector<std::vector<size_t>>& result_paths,
    size_t max_paths
) {
    // Early termination if we have enough paths
    if (result_paths.size() >= max_paths) {
        return;
    }

    // Base case: reached target
    if (current == target) {
        result_paths.push_back(current_path);
        return;
    }

    // Recursive case: try each predecessor
    auto it = all_predecessors.find(current);
    if (it != all_predecessors.end()) {
        for (size_t pred : it->second) {
            // Check we haven't already visited this vertex (cycle detection)
            if (std::find(current_path.begin(), current_path.end(), pred)
                != current_path.end()) {
                continue;
            }

            current_path.push_back(pred);
            enumerate_k_paths_recursive(pred, target, all_predecessors,
                                       current_path, result_paths, max_paths);
            current_path.pop_back();

            // Early exit if we have enough paths
            if (result_paths.size() >= max_paths) {
                return;
            }
        }
    }
}

std::vector<std::vector<size_t>> enumerate_k_shortest_paths(
    size_t terminal,
    size_t origin,
    const std::unordered_map<size_t, std::vector<size_t>>& all_predecessors,
    size_t max_paths
) {
    std::vector<std::vector<size_t>> result_paths;
    std::vector<size_t> current_path = {terminal};

    enumerate_k_paths_recursive(terminal, origin, all_predecessors,
                                current_path, result_paths, max_paths);

    return result_paths;
}

/**
 * @brief Compute gradient basin with K-shortest path trajectory enumeration.
 *
 * Performs breadth-first search from a local extremum to identify all vertices
 * reachable via monotone paths (basin membership). Then constructs the complete
 * monotone directed acyclic graph including all monotone edges within the basin.
 * Terminals are identified as local extrema within the basin (m-terminals), ensuring
 * all trajectories represent complete gradient descent to true local minima.
 *
 * For computational tractability, enumerates only the K shortest paths from each
 * terminal (default K=20), providing sufficient statistical power while maintaining
 * efficiency. Path length is measured by number of edges (hop distance).
 *
 * @param vertex Index of the extremum vertex (basin origin)
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basins (from maximum), false for ascending
 * @param with_trajectories If true, enumerate K-shortest paths from m-terminals
 * @param k_paths Number of shortest paths to enumerate per terminal (default: 20)
 *
 * @return gradient_basin_t Complete basin structure with K-shortest trajectories
 */
gradient_basin_t set_wgraph_t::compute_basin_of_attraction(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima,
    bool with_trajectories,
    size_t k_paths
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
            // For maximum: all neighbors must have strictly smaller values
            if (y[u] >= y[vertex]) {
                is_local_extremum = false;
                break;
            }
        } else {
            // For minimum: all neighbors must have strictly larger values
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

    // BFS to determine basin membership
    // Hop distance used only for diagnostics, not for predecessor selection
    std::vector<int> hop_distance(adjacency_list.size(), -1);
    std::queue<size_t> q;

    result.hop_dist_map[vertex] = 0;
    hop_distance[vertex] = 0;
    q.push(vertex);

    size_t max_hop_distance = 0;

    // BFS: identify all vertices reachable via monotone edges
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();

        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (hop_distance[v] != -1) {
                continue;  // Already visited
            }

            // Check monotonicity condition for edge u -> v
            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) {
                continue;  // Skip non-monotone edges
            }

            // Add v to basin
            hop_distance[v] = hop_distance[u] + 1;
            result.hop_dist_map[v] = static_cast<size_t>(hop_distance[v]);
            q.push(v);

            max_hop_distance = std::max(max_hop_distance,
                                       static_cast<size_t>(hop_distance[v]));
        }
    }

    result.hop_idx = max_hop_distance;

    // Build complete monotone DAG: include ALL monotone edges within basin
    // This captures the full gradient flow structure, not just BFS tree edges
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

    // Identify m-terminals: local minima/maxima within the basin
    // These represent true gradient descent endpoints
    std::vector<size_t> m_terminals;
    for (const auto& [u, hop] : result.hop_dist_map) {
        if (u == vertex) continue;  // Skip origin

        bool is_local_extremum_in_basin = true;
        for (const auto& edge : adjacency_list[u]) {
            size_t w = edge.vertex;

            // Only consider neighbors within the basin
            if (result.hop_dist_map.find(w) == result.hop_dist_map.end()) {
                continue;
            }

            if (detect_maxima) {
                // For descending basin: terminal if no downhill neighbors in basin
                if (y[w] < y[u]) {
                    is_local_extremum_in_basin = false;
                    break;
                }
            } else {
                // For ascending basin: terminal if no uphill neighbors in basin
                if (y[w] > y[u]) {
                    is_local_extremum_in_basin = false;
                    break;
                }
            }
        }

        if (is_local_extremum_in_basin) {
            m_terminals.push_back(u);
        }
    }

    // Use m-terminals for cell structure
    result.terminal_extrema = m_terminals;

    // Diagnostic output
    if (with_trajectories) {
        Rprintf("\n=== Basin Structure Diagnostic ===\n");
        Rprintf("Basin vertex: %zu (y=%.6f)\n", vertex, y[vertex]);
        Rprintf("Basin size: %zu vertices\n", result.hop_dist_map.size());
        Rprintf("Max hop distance: %zu\n", max_hop_distance);
        Rprintf("m-terminals: %zu\n", m_terminals.size());
        Rprintf("Strategy: K-shortest paths (K=%zu per terminal)\n", k_paths);
        Rprintf("Using FULL monotone DAG (all monotone edges)\n");
        Rprintf("===================================\n\n");
    }

    // Enumerate K-shortest paths from each m-terminal
    if (with_trajectories) {
        size_t total_paths = 0;
        size_t terminals_at_limit = 0;

        for (size_t terminal : result.terminal_extrema) {
            trajectory_set_t traj_set;
            traj_set.terminal_vertex = terminal;

            // Enumerate K shortest paths from this terminal
            traj_set.trajectories = enumerate_k_shortest_paths(
                terminal,
                vertex,
                result.all_predecessors,
                k_paths
            );

            total_paths += traj_set.trajectories.size();
            if (traj_set.trajectories.size() >= k_paths) {
                terminals_at_limit++;
            }

            result.trajectory_sets.push_back(traj_set);
        }

        // Verify trajectory coverage
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
        size_t uncovered = basin_size - covered;

        Rprintf("\n=== Trajectory Coverage Diagnostic ===\n");
        Rprintf("Basin vertex: %zu, Basin size: %zu\n", vertex, basin_size);
        Rprintf("Vertices in trajectories: %zu\n", covered);
        Rprintf("Uncovered vertices: %zu\n", uncovered);

        if (uncovered == 0) {
            Rprintf("SUCCESS: Perfect trajectory coverage (100%%)!\n");
        } else {
            Rprintf("WARNING: Incomplete coverage (%.1f%%)\n",
                   100.0 * covered / basin_size);
            Rprintf("Consider increasing K (currently %zu) for complete coverage\n", k_paths);
        }

        // Path enumeration statistics
        size_t min_paths = std::numeric_limits<size_t>::max();
        size_t max_paths_found = 0;
        std::vector<size_t> path_lengths;

        for (const auto& traj_set : result.trajectory_sets) {
            size_t n = traj_set.trajectories.size();
            if (n > 0) {
                min_paths = std::min(min_paths, n);
                max_paths_found = std::max(max_paths_found, n);
            }

            // Collect path lengths for statistics
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

            double avg_length = 0.0;
            for (size_t len : path_lengths) {
                avg_length += len;
            }
            avg_length /= path_lengths.size();

            Rprintf("\nPath enumeration statistics:\n");
            Rprintf("  Total paths: %zu\n", total_paths);
            Rprintf("  Terminals at K-limit: %zu / %zu\n",
                   terminals_at_limit, result.terminal_extrema.size());
            Rprintf("  Paths per terminal: min=%zu, avg=%.1f, max=%zu\n",
                   min_paths,
                   static_cast<double>(total_paths) / result.terminal_extrema.size(),
                   max_paths_found);
            Rprintf("  Path lengths: min=%zu, median=%zu, mean=%.1f, max=%zu\n",
                   min_length, median_length, avg_length, max_length);
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
