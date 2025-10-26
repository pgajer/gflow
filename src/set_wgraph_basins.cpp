#include <queue>
#include <random>
#include "set_wgraph.hpp"
#include "error_utils.h"
#include <R_ext/Print.h>

/**
 * @brief Trace path from vertex upward to origin following predecessors
 *
 * @param start Starting vertex
 * @param origin Target vertex (basin maximum/minimum)
 * @param all_predecessors Map from vertex to predecessors
 * @return Path from start to origin (inclusive)
 */
std::vector<size_t> trace_path_to_origin(
    size_t start,
    size_t origin,
    const std::unordered_map<size_t, std::vector<size_t>>& all_predecessors
) {
    std::vector<size_t> path;
    size_t current = start;
    std::unordered_set<size_t> visited;

    while (current != origin) {
        path.push_back(current);
        visited.insert(current);

        auto it = all_predecessors.find(current);
        if (it == all_predecessors.end() || it->second.empty()) {
            // Shouldn't happen in a valid basin, but be safe
            Rprintf("WARNING: trace_path_to_origin stuck at vertex %zu\n", current);
            break;
        }

        // Pick first predecessor arbitrarily
        size_t next = it->second[0];

        // Cycle detection
        if (visited.find(next) != visited.end()) {
            Rprintf("WARNING: cycle detected at vertex %zu\n", next);
            break;
        }

        current = next;
    }

    path.push_back(origin);
    return path;
}

/**
 * @brief Trace path from vertex downward to a terminal following successors
 *
 * @param start Starting vertex
 * @param successors Map from vertex to successors
 * @param m_terminals Set of m-terminal vertices
 * @return Path from start to terminal (inclusive), with terminal at end
 */
std::vector<size_t> trace_path_to_terminal(
    size_t start,
    const std::unordered_map<size_t, std::vector<size_t>>& successors,
    const std::unordered_set<size_t>& m_terminals
) {
    std::vector<size_t> path;
    size_t current = start;
    std::unordered_set<size_t> visited;

    // If start is already a terminal, return single-vertex path
    if (m_terminals.find(start) != m_terminals.end()) {
        path.push_back(start);
        return path;
    }

    while (m_terminals.find(current) == m_terminals.end()) {
        path.push_back(current);
        visited.insert(current);

        auto it = successors.find(current);
        if (it == successors.end() || it->second.empty()) {
            // Reached a leaf that's not an m-terminal (shouldn't happen)
            Rprintf("WARNING: trace_path_to_terminal stuck at non-terminal vertex %zu\n",
                   current);
            break;
        }

        // Pick first successor arbitrarily
        size_t next = it->second[0];

        // Cycle detection
        if (visited.find(next) != visited.end()) {
            Rprintf("WARNING: cycle detected at vertex %zu\n", next);
            break;
        }

        current = next;
    }

    path.push_back(current);  // Add the terminal
    return path;
}

/**
 * @brief Enumerate paths with early stopping after K paths found
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
            // Cycle detection
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
 * @brief Compute gradient basin with K-shortest paths plus gap-filling for complete coverage
 */
gradient_basin_t set_wgraph_t::compute_basin_of_attraction(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima,
    bool with_trajectories,
    size_t k_paths
    ) const {

    #define DEBUG__compute_basin_of_attraction 0

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

    // BFS to determine basin membership
    std::vector<int> hop_distance(adjacency_list.size(), -1);
    std::queue<size_t> q;

    result.hop_dist_map[vertex] = 0;
    hop_distance[vertex] = 0;
    q.push(vertex);

    size_t max_hop_distance = 0;

    while (!q.empty()) {
        size_t u = q.front();
        q.pop();

        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (hop_distance[v] != -1) {
                continue;
            }

            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) {
                continue;
            }

            hop_distance[v] = hop_distance[u] + 1;
            result.hop_dist_map[v] = static_cast<size_t>(hop_distance[v]);
            q.push(v);

            max_hop_distance = std::max(max_hop_distance,
                                       static_cast<size_t>(hop_distance[v]));
        }
    }

    result.hop_idx = max_hop_distance;

    // Build complete monotone DAG
    for (const auto& [v, _] : result.hop_dist_map) {
        result.all_predecessors[v] = std::vector<size_t>();
    }

    for (const auto& [u, hop_u] : result.hop_dist_map) {
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (result.hop_dist_map.find(v) == result.hop_dist_map.end()) {
                continue;
            }

            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (monotone) {
                result.all_predecessors[v].push_back(u);
            }
        }
    }

    // Build successor map from predecessors
    std::unordered_map<size_t, std::vector<size_t>> successors;
    for (const auto& [v, _] : result.hop_dist_map) {
        successors[v] = std::vector<size_t>();
    }

    for (const auto& [v, preds] : result.all_predecessors) {
        for (size_t pred : preds) {
            successors[pred].push_back(v);
        }
    }

    // Identify m-terminals
    std::vector<size_t> m_terminals;
    for (const auto& [u, hop] : result.hop_dist_map) {
        if (u == vertex) continue;

        bool is_local_extremum_in_basin = true;
        for (const auto& edge : adjacency_list[u]) {
            size_t w = edge.vertex;

            if (result.hop_dist_map.find(w) == result.hop_dist_map.end()) {
                continue;
            }

            if (detect_maxima) {
                if (y[w] < y[u]) {
                    is_local_extremum_in_basin = false;
                    break;
                }
            } else {
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

    result.terminal_extrema = m_terminals;

	#if DEBUG__compute_basin_of_attraction
    if (with_trajectories) {
        Rprintf("\n=== Basin Structure Diagnostic ===\n");
        Rprintf("Basin vertex: %zu (y=%.6f)\n", vertex, y[vertex]);
        Rprintf("Basin size: %zu vertices\n", result.hop_dist_map.size());
        Rprintf("Max hop distance: %zu\n", max_hop_distance);
        Rprintf("m-terminals: %zu\n", m_terminals.size());
        Rprintf("Strategy: K-shortest paths (K=%zu) + gap-filling\n", k_paths);
        Rprintf("===================================\n\n");
    }
	#endif

    // STEP 1: Enumerate K paths from each m-terminal
    if (with_trajectories) {
        size_t initial_paths = 0;

        for (size_t terminal : result.terminal_extrema) {
            trajectory_set_t traj_set;
            traj_set.terminal_vertex = terminal;

            traj_set.trajectories = enumerate_k_shortest_paths(
                terminal,
                vertex,
                result.all_predecessors,
                k_paths
            );

            initial_paths += traj_set.trajectories.size();
            result.trajectory_sets.push_back(traj_set);
        }

		#if DEBUG__compute_basin_of_attraction
        Rprintf("Phase 1: Enumerated %zu paths from %zu m-terminals\n",
				initial_paths, m_terminals.size());
		#endif

        // STEP 2: Identify uncovered vertices
        std::unordered_set<size_t> covered_vertices;
        for (const auto& traj_set : result.trajectory_sets) {
            for (const auto& traj : traj_set.trajectories) {
                for (size_t v : traj) {
                    covered_vertices.insert(v);
                }
            }
        }

        std::vector<size_t> uncovered;
        for (const auto& [v, _] : result.hop_dist_map) {
            if (covered_vertices.find(v) == covered_vertices.end()) {
                uncovered.push_back(v);
            }
        }

        size_t initial_uncovered = uncovered.size();

		#if DEBUG__compute_basin_of_attraction
        Rprintf("Initial coverage: %zu / %zu vertices (%.1f%%)\n",
               covered_vertices.size(), result.hop_dist_map.size(),
               100.0 * covered_vertices.size() / result.hop_dist_map.size());
        #endif

        if (initial_uncovered > 0) {
			#if DEBUG__compute_basin_of_attraction
            Rprintf("Phase 2: Gap-filling for %zu uncovered vertices\n", initial_uncovered);
			#endif

            // Randomize order for unbiased selection
			std::random_device rd;
			std::mt19937 rng(rd());
			std::shuffle(uncovered.begin(), uncovered.end(), rng);

            // Create set of m-terminals for quick lookup
            std::unordered_set<size_t> m_terminal_set(m_terminals.begin(), m_terminals.end());

            // Map from terminal to its trajectory_set index
            std::unordered_map<size_t, size_t> terminal_to_idx;
            for (size_t i = 0; i < result.terminal_extrema.size(); ++i) {
                terminal_to_idx[result.terminal_extrema[i]] = i;
            }

            size_t gap_fill_count = 0;

            for (size_t v : uncovered) {
                // Skip if now covered by a previous gap-fill
                if (covered_vertices.find(v) != covered_vertices.end()) {
                    continue;
                }

                // Trace upward to origin
                std::vector<size_t> ascending_path = trace_path_to_origin(
                    v, vertex, result.all_predecessors);

                // Trace downward to terminal
                std::vector<size_t> descending_path = trace_path_to_terminal(
                    v, successors, m_terminal_set);

                if (descending_path.empty()) {
                    Rprintf("WARNING: Could not trace vertex %zu to terminal\n", v);
                    continue;
                }

                // Get terminal at end of descending path
                size_t terminal = descending_path.back();

                // Construct complete path: [terminal, ..., v, ..., origin]
                // descending_path goes from v to terminal, so reverse it
                std::vector<size_t> complete_path;
                for (auto it = descending_path.rbegin(); it != descending_path.rend(); ++it) {
                    complete_path.push_back(*it);
                }
                // ascending_path goes from v to origin, skip v (already included)
                for (size_t i = 1; i < ascending_path.size(); ++i) {
                    complete_path.push_back(ascending_path[i]);
                }

                // Add to appropriate terminal's trajectory set
                auto idx_it = terminal_to_idx.find(terminal);
                if (idx_it != terminal_to_idx.end()) {
                    result.trajectory_sets[idx_it->second].trajectories.push_back(complete_path);
                    gap_fill_count++;

                    // Mark all vertices in this path as covered
                    for (size_t vertex_in_path : complete_path) {
                        covered_vertices.insert(vertex_in_path);
                    }
                }
            }

            #if DEBUG__compute_basin_of_attraction
            Rprintf("Gap-filling: Added %zu paths, covered %zu additional vertices\n",
					gap_fill_count, covered_vertices.size() - (result.hop_dist_map.size() - initial_uncovered));
			#endif
        }

        // Final coverage check
        size_t final_uncovered = 0;
        for (const auto& [v, _] : result.hop_dist_map) {
            if (covered_vertices.find(v) == covered_vertices.end()) {
                final_uncovered++;
            }
        }

		#if DEBUG__compute_basin_of_attraction
        Rprintf("\n=== Final Coverage Summary ===\n");
        Rprintf("Total vertices: %zu\n", result.hop_dist_map.size());
        Rprintf("Covered vertices: %zu\n", covered_vertices.size());
        Rprintf("Uncovered vertices: %zu\n", final_uncovered);

        if (final_uncovered == 0) {
            Rprintf("SUCCESS: Perfect coverage (100%%)!\n");
        } else {
            Rprintf("WARNING: Incomplete coverage (%.1f%%)\n",
                   100.0 * covered_vertices.size() / result.hop_dist_map.size());
        }
        #endif

        // Path statistics
        size_t total_paths = 0;
        std::vector<size_t> path_lengths;

        for (const auto& traj_set : result.trajectory_sets) {
            total_paths += traj_set.trajectories.size();
            for (const auto& traj : traj_set.trajectories) {
                path_lengths.push_back(traj.size());
            }
        }

		#if DEBUG__compute_basin_of_attraction
        if (!path_lengths.empty()) {
            std::sort(path_lengths.begin(), path_lengths.end());
            size_t median_length = path_lengths[path_lengths.size() / 2];
            size_t min_length = path_lengths.front();
            size_t max_length = path_lengths.back();

            double avg_length = 0.0;
            for (size_t len : path_lengths) {
                avg_length += len;
            }
            avg_length /= path_lengths.size();

            Rprintf("\nPath statistics:\n");
            Rprintf("  Total paths: %zu\n", total_paths);
            Rprintf("  Avg paths per terminal: %.1f\n",
                   static_cast<double>(total_paths) / result.terminal_extrema.size());
            Rprintf("  Path lengths: min=%zu, median=%zu, mean=%.1f, max=%zu\n",
					min_length, median_length, avg_length, max_length);
        }
		Rprintf("================================\n\n");
		#endif
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
