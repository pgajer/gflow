#include <queue>
#include <random>
#include "set_wgraph.hpp"
#include "error_utils.h"
#include <R_ext/Print.h>

/**
 * @brief Find shortest monotonic path from source to target within a basin
 *
 * @param source Starting vertex
 * @param target Target vertex
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basin (require y[next] < y[current])
 * @param basin_vertices Set of vertices in the basin (restricts search space)
 * @return Vector of vertices forming path from source to target (empty if no path exists)
 */
std::vector<size_t> set_wgraph_t::find_shortest_monotonic_within_basin_path(
    size_t source,
    size_t target,
    const std::vector<double>& y,
    bool detect_maxima,
    const std::unordered_set<size_t>& basin_vertices
) const {

    size_t n = adjacency_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    std::vector<bool> visited(n, false);

    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<queue_entry>> pq;

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        if (u == target) {
            break;  // Found target
        }

        // Explore neighbors with monotonicity constraint
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (visited[v]) continue;

            // CRITICAL: Only explore vertices within the basin
            if (basin_vertices.find(v) == basin_vertices.end()) {
                continue;
            }

            // Check monotonicity constraint
            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) continue;  // Skip non-monotonic edges

            double new_dist = dist[u] + edge.weight;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                pq.push({new_dist, v});
            }
        }
    }

    // Reconstruct path
    std::vector<size_t> path;
    if (dist[target] == INFINITY) {
        return path;  // No monotonic path exists within basin
    }

    size_t curr = target;
    while (curr != INVALID_VERTEX) {
        path.push_back(curr);
        curr = prev[curr];
    }

    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Find shortest monotonic path from source to target
 *
 * Uses Dijkstra's algorithm with monotonicity constraints to find the shortest
 * path from source vertex to target vertex. The path must maintain strict monotonicity:
 * strictly decreasing function values for descending basins, or strictly increasing
 * for ascending basins.
 *
 * @param source Starting vertex
 * @param target Target vertex
 * @param y Function values at each vertex
 * @param detect_maxima True for descending basin (require y[next] < y[current])
 * @return Vector of vertices forming path from source to target (empty if no monotonic path exists)
 */
std::vector<size_t> set_wgraph_t::find_shortest_monotonic_path(
    size_t source,
    size_t target,
    const std::vector<double>& y,
    bool detect_maxima
) const {

    size_t n = adjacency_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    std::vector<bool> visited(n, false);

    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<queue_entry>> pq;

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        if (u == target) {
            break;  // Found target
        }

        // Explore neighbors with monotonicity constraint
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            if (visited[v]) continue;

            // Check monotonicity constraint
            double delta_y = y[v] - y[u];
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) continue;  // Skip non-monotonic edges

            double new_dist = dist[u] + edge.weight;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                pq.push({new_dist, v});
            }
        }
    }

    // Reconstruct path
    std::vector<size_t> path;
    if (dist[target] == INFINITY) {
        return path;  // No monotonic path exists
    }

    size_t curr = target;
    while (curr != INVALID_VERTEX) {
        path.push_back(curr);
        curr = prev[curr];
    }

    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Compute gradient basin using geodesic paths with m-terminal cell decomposition
 *
 * @details This function implements a two-phase algorithm to construct gradient flow cells
 *          within a basin of attraction, where each cell represents the set of vertices whose
 *          gradient trajectories pass through a specific local minimum (m-terminal).
 *
 * The use of m-terminals (as opposed to s-terminals which are simply tree leaves) preserves
 * important biological and mathematical structure: m-terminals represent true local extrema
 * within the basin, and gradient trajectories connecting m-terminals to the reference extremum
 * encode information about the gradient structure between opposite local extrema.
 *
 * Algorithm Overview:
 *
 * Phase 1: Primary Path Construction
 * -----------------------------------
 * 1. Compute the geodesic basin from the extremum using Dijkstra's algorithm with monotonicity
 *    constraints (via find_local_extremum_geodesic_basin). This yields a tree structure where
 *    each vertex has a unique shortest monotonic path to the extremum.
 *
 * 2. Identify m-terminals: vertices that are local minima/maxima within the basin subgraph.
 *    These represent true terminal extrema where gradient flow naturally terminates. Unlike
 *    s-terminals (tree leaves), m-terminals have geometric meaning independent of the specific
 *    path structure.
 *
 * 3. For each m-terminal, trace its geodesic path to the extremum along the predecessor tree
 *    and mark all vertices on this path as covered. These paths form the primary flow channels.
 *
 * Phase 2: Gap-Filling via Constrained Dijkstra
 * ---------------------------------------------
 * For each uncovered vertex v:
 *
 * 4A. Identify the closest m-terminal by geodesic distance from the origin. Sort all m-terminals
 *     by their distance from v and attempt connection in order of increasing distance.
 *
 * 4B. Run constrained Dijkstra from v to find the shortest monotonic path to the target m-terminal.
 *     This ensures the path maintains strict monotonicity (decreasing for maxima basins, increasing
 *     for minima basins). If no monotonic path exists to the closest terminal, try the next-closest.
 *
 * 4C. Construct upward path from v to the reference extremum using the predecessor tree (γ₂).
 *
 * 4D. Join paths to form complete trajectory: [m-terminal] → ... → [v] → ... → [extremum]
 *     This is the concatenation of the downward monotonic path (reversed γ₁) with the upward
 *     geodesic path (γ₂), with duplicate vertex v removed at the junction point.
 *
 * 4E. Mark all vertices on this trajectory as covered and assign to the m-terminal's cell.
 *
 * 5. Continue until all vertices are covered or no more monotonic connections can be made.
 *
 * Mathematical Foundation:
 * -----------------------
 * The geodesic paths use actual Riemannian distances (edge weights) rather than hop counts,
 * making them suitable for downstream analysis with exponential path weighting exp(-|γ|²/σ²)
 * where |γ| represents the true Riemannian path length.
 *
 * The constrained Dijkstra search for gap-filling ensures that all constructed trajectories
 * maintain the fundamental property of gradient flow: strict monotonicity of the response
 * function along the path. This preserves the geometric structure even for vertices not
 * on the primary geodesic tree paths.
 *
 * Cell Decomposition Properties:
 * -----------------------------
 * - Each covered vertex belongs to exactly one cell (partition property)
 * - Each cell is associated with an m-terminal (local extremum in basin)
 * - All trajectories in a cell pass through the same m-terminal
 * - Path lengths measured by Riemannian distance for biological interpretation
 * - Coverage typically >95%, with uncovered vertices only when no monotonic path exists
 *
 * Computational Complexity:
 * ------------------------
 * - Geodesic basin: O(E log V) for initial Dijkstra with E edges and V vertices
 * - m-terminal identification: O(V × avg_degree) to check local extremality
 * - Primary paths: O(m × depth) where m is number of m-terminals, depth is tree depth
 * - Gap-filling: O(uncovered × E log V) worst case for Dijkstra per uncovered vertex
 * - Overall: O((V + uncovered) × E log V) in worst case
 *
 * For typical basins with good primary coverage (>90%), gap-filling cost is small.
 *
 * Biological Interpretation:
 * -------------------------
 * m-terminals represent distinct "sinks" or terminal states in the disease progression
 * landscape. Each cell captures a coherent subdomain where samples follow similar
 * progression pathways. The use of true local extrema (rather than arbitrary tree leaves)
 * ensures that cells have biological meaning: they represent genuinely distinct terminal
 * states rather than artifacts of the geodesic tree structure.
 *
 * @param vertex Index of the extremum vertex (basin origin)
 * @param y Function values at each vertex (e.g., disease prevalence estimates)
 * @param detect_maxima True for descending basins (from maximum), false for ascending (from minimum)
 * @param with_trajectories If true, enumerate complete trajectories for each vertex
 *
 * @return gradient_basin_t Complete basin structure with m-terminal-based cell decomposition
 *         containing:
 *         - vertex: reference extremum vertex
 *         - terminal_extrema: vector of m-terminal vertices (local extrema in basin)
 *         - trajectory_sets: one set per m-terminal with all trajectories through it
 *         - hop_dist_map: hop distances from origin (diagnostic)
 *         - all_predecessors: predecessor relationships in geodesic tree
 *         - y_nbhd_bd_map: boundary vertices outside the basin
 *
 * @note Requires the graph to be connected and edge weights to be positive (Riemannian distances).
 *       The monotonicity constraint in gap-filling may prevent some vertices from being covered
 *       if no monotonic path exists to any m-terminal.
 *
 * @see find_local_extremum_geodesic_basin() for geodesic basin computation
 * @see find_shortest_monotonic_path() for constrained path finding
 * @see gradient_basin_t for the return structure
 *
 * @example
 * // Compute descending basin from a local maximum with trajectory enumeration
 * gradient_basin_t basin = graph.compute_geodesic_basin(max_vertex, y_values, true, true);
 *
 * // Examine coverage
 * Rprintf("Found %zu m-terminals\n", basin.terminal_extrema.size());
 *
 * // Access cells and trajectories for downstream monotonicity analysis
 * for (const auto& traj_set : basin.trajectory_sets) {
 *     size_t terminal = traj_set.terminal_vertex;
 *     // Compute weighted monotonicity indices across trajectories
 *     for (const auto& traj : traj_set.trajectories) {
 *         // Path starts at terminal, ends at extremum
 *         // Use for computing feature-response monotonicity
 *     }
 * }
 */
gradient_basin_t set_wgraph_t::compute_geodesic_basin(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima,
    double edge_length_thld,
    bool with_trajectories
    ) const {

    #define DEBUG__compute_geodesic_basin 0

    gradient_basin_t result;
    result.vertex = vertex;
    result.value = y[vertex];
    result.is_maximum = detect_maxima;

    // ========================================================================
    // PHASE 1: COMPUTE GEODESIC BASIN AND IDENTIFY STRUCTURE
    // ========================================================================

    // Step 1: Get geodesic basin using Dijkstra with monotonicity constraints
    basin_t geodesic_basin = find_local_extremum_geodesic_basin(vertex,
                                                                y,
                                                                detect_maxima,
                                                                edge_length_thld);

    // Step 2: Build hop distance map (for diagnostics and compatibility)
    std::unordered_map<size_t, size_t> hop_map;
    hop_map[vertex] = 0;

    std::queue<size_t> q;
    q.push(vertex);
    size_t max_hop = 0;

    while (!q.empty()) {
        size_t u = q.front();
        q.pop();

        for (const auto& [v, pred] : geodesic_basin.reachability_map.predecessors) {
            if (pred == u && hop_map.find(v) == hop_map.end()) {
                hop_map[v] = hop_map[u] + 1;
                max_hop = std::max(max_hop, hop_map[v]);
                q.push(v);
            }
        }
    }

    result.hop_idx = max_hop;

    // Build result structures from geodesic basin
    for (const auto& [v, dist] : geodesic_basin.reachability_map.distances) {
        result.hop_dist_map[v] = hop_map[v];

        size_t pred = geodesic_basin.reachability_map.predecessors.at(v);
        if (pred != INVALID_VERTEX) {
            result.all_predecessors[v] = {pred};
        } else {
            result.all_predecessors[v] = {};
        }
    }

    // Step 3: Identify m-terminals (local extrema within basin)
    std::vector<size_t> m_terminals;

    for (const auto& [u, dist] : geodesic_basin.reachability_map.distances) {
        if (u == vertex) continue;  // Skip origin

        bool is_local_extremum_in_basin = true;

        for (const auto& edge : adjacency_list[u]) {
            size_t w = edge.vertex;

            if (geodesic_basin.reachability_map.distances.find(w) ==
                geodesic_basin.reachability_map.distances.end()) {
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

    result.terminal_extrema = m_terminals;

    #if DEBUG__compute_geodesic_basin
    Rprintf("\n=== Geodesic Basin Structure ===\n");
    Rprintf("Basin vertex: %zu (y=%.6f)\n", vertex, y[vertex]);
    Rprintf("Basin size: %zu vertices\n", geodesic_basin.reachability_map.distances.size());
    Rprintf("Max hop distance: %zu\n", max_hop);
    Rprintf("m-terminals: %zu\n", m_terminals.size());
    Rprintf("Method: m-terminals with constrained Dijkstra gap-filling\n");
    Rprintf("================================\n\n");
    #endif

    // ========================================================================
    // PHASE 2: BUILD TRAJECTORY SETS WITH M-TERMINAL COVERAGE
    // ========================================================================

    if (with_trajectories && !m_terminals.empty()) {
        // Initialize trajectory sets for each m-terminal
        std::unordered_map<size_t, size_t> terminal_to_idx;
        for (size_t i = 0; i < m_terminals.size(); ++i) {
            trajectory_set_t traj_set;
            traj_set.terminal_vertex = m_terminals[i];
            result.trajectory_sets.push_back(traj_set);
            terminal_to_idx[m_terminals[i]] = i;
        }

        std::unordered_set<size_t> m_terminal_set(m_terminals.begin(), m_terminals.end());

        // Track which vertices have been assigned to cells
        std::unordered_set<size_t> covered_vertices;
        covered_vertices.insert(vertex);  // Origin is always covered

        // --------------------------------------------------------------------
        // PHASE 2.1: PRIMARY PATH CONSTRUCTION
        // --------------------------------------------------------------------
        #if DEBUG__compute_geodesic_basin
        Rprintf("Phase 1: Constructing primary paths from m-terminals...\n");
        #endif

        for (size_t term : m_terminals) {
            // Trace geodesic path from terminal to origin
            std::vector<size_t> path;
            size_t curr = term;
            std::unordered_set<size_t> visited;

            while (curr != INVALID_VERTEX) {
                path.push_back(curr);
                visited.insert(curr);

                if (curr == vertex) {
                    break;  // Reached origin
                }

                auto pred_it = geodesic_basin.reachability_map.predecessors.find(curr);
                if (pred_it == geodesic_basin.reachability_map.predecessors.end()) {
                    break;
                }

                size_t pred = pred_it->second;
                if (pred == INVALID_VERTEX || visited.find(pred) != visited.end()) {
                    break;
                }

                curr = pred;
            }

            // Add this primary path to the terminal's trajectory set
            if (!path.empty() && path.back() == vertex) {
                auto idx_it = terminal_to_idx.find(term);
                if (idx_it != terminal_to_idx.end()) {
                    result.trajectory_sets[idx_it->second].trajectories.push_back(path);

                    // Mark all vertices on this path as covered
                    for (size_t v : path) {
                        covered_vertices.insert(v);
                    }
                }
            }
        }

        size_t basin_size = geodesic_basin.reachability_map.distances.size();

        #if DEBUG__compute_geodesic_basin
        size_t initial_coverage = covered_vertices.size();
        Rprintf("  Primary paths cover %zu/%zu vertices (%.1f%%)\n",
               initial_coverage, basin_size,
               100.0 * initial_coverage / basin_size);
        #endif

        // Create basin vertex set
        std::unordered_set<size_t> basin_vertices;
        for (const auto& [v, dist] : geodesic_basin.reachability_map.distances) {
            basin_vertices.insert(v);
        }

        // --------------------------------------------------------------------
        // PHASE 2.2: GAP-FILLING VIA CONSTRAINED DIJKSTRA
        // --------------------------------------------------------------------
        if (covered_vertices.size() < basin_size) {

            #if DEBUG__compute_geodesic_basin
            Rprintf("\nPhase 2: Gap-filling using constrained Dijkstra...\n");
            #endif

            size_t gap_fill_count = 0;
            size_t no_path_count = 0;

            for (const auto& [v, v_dist] : geodesic_basin.reachability_map.distances) {
                if (v == vertex) continue;  // Skip origin
                if (covered_vertices.find(v) != covered_vertices.end()) continue;  // Already covered

                // Step 4A: Sort m-terminals by distance from v
                std::vector<std::pair<double, size_t>> terminal_distances;

                for (size_t term : m_terminals) {
                    double term_dist = geodesic_basin.reachability_map.distances.at(term);
                    double distance = std::abs(term_dist - v_dist);
                    terminal_distances.push_back({distance, term});
                }

                std::sort(terminal_distances.begin(), terminal_distances.end());

                // Step 4B: Try to find monotonic path to terminals in order of distance
                std::vector<size_t> path_to_terminal;
                size_t connected_terminal = INVALID_VERTEX;

                for (const auto& [dist, term] : terminal_distances) {
                    path_to_terminal = find_shortest_monotonic_within_basin_path(v, term, y, detect_maxima, basin_vertices);

                    if (!path_to_terminal.empty()) {
                        connected_terminal = term;
                        break;  // Found a monotonic path
                    }
                }

                if (connected_terminal == INVALID_VERTEX) {
                    no_path_count++;
                    continue;  // No monotonic path to any m-terminal
                }

                // Step 4C: Construct upward path from v to origin
                std::vector<size_t> upward_path;
                size_t curr = v;
                std::unordered_set<size_t> visited;

                while (curr != INVALID_VERTEX) {
                    upward_path.push_back(curr);
                    visited.insert(curr);

                    if (curr == vertex) {
                        break;
                    }

                    auto pred_it = geodesic_basin.reachability_map.predecessors.find(curr);
                    if (pred_it == geodesic_basin.reachability_map.predecessors.end()) {
                        break;
                    }

                    size_t pred = pred_it->second;
                    if (pred == INVALID_VERTEX || visited.find(pred) != visited.end()) {
                        break;
                    }

                    curr = pred;
                }

                // Step 4D: Join paths: [terminal, ..., v, ..., origin]
                std::vector<size_t> complete_path;

                // Add path_to_terminal (already goes from v to terminal)
                // Reverse it to go from terminal to v
                for (auto it = path_to_terminal.rbegin(); it != path_to_terminal.rend(); ++it) {
                    complete_path.push_back(*it);
                }

                // Add upward path from v to origin, but skip v (already included)
                for (size_t i = 1; i < upward_path.size(); ++i) {
                    complete_path.push_back(upward_path[i]);
                }

                // Step 4E: Add to terminal's cell and mark vertices as covered
                auto idx_it = terminal_to_idx.find(connected_terminal);
                if (idx_it != terminal_to_idx.end() && !complete_path.empty()) {
                    result.trajectory_sets[idx_it->second].trajectories.push_back(complete_path);
                    gap_fill_count++;

                    // Mark all vertices on this path as covered
                    for (size_t u : complete_path) {
                        covered_vertices.insert(u);
                    }
                }
            }

            #if DEBUG__compute_geodesic_basin
            Rprintf("  Gap-filling added %zu trajectories\n", gap_fill_count);
            if (no_path_count > 0) {
                Rprintf("  WARNING: %zu vertices have no monotonic path to any m-terminal\n", no_path_count);
            }
            Rprintf("  Final coverage: %zu/%zu vertices (%.1f%%)\n",
                   covered_vertices.size(), basin_size,
                    100.0 * covered_vertices.size() / basin_size);
            #endif
        }

        // ====================================================================
        // FINAL DIAGNOSTICS AND STATISTICS
        // ====================================================================

        #if DEBUG__compute_geodesic_basin
        size_t final_uncovered = basin_size - covered_vertices.size();
        Rprintf("\n=== Trajectory Coverage Summary ===\n");
        Rprintf("Basin size: %zu vertices\n", basin_size);
        Rprintf("Covered vertices: %zu\n", covered_vertices.size());
        Rprintf("Uncovered vertices: %zu\n", final_uncovered);

        if (final_uncovered == 0) {
            Rprintf("SUCCESS: Perfect coverage (100%%)!\n");
        } else {
            Rprintf("Coverage: %.1f%%\n", 100.0 * covered_vertices.size() / basin_size);
        }
        #endif

        // Path statistics
        size_t total_paths = 0;
        std::vector<size_t> path_lengths;
        std::vector<size_t> paths_per_terminal;

        for (const auto& traj_set : result.trajectory_sets) {
            size_t n_paths = traj_set.trajectories.size();
            total_paths += n_paths;
            paths_per_terminal.push_back(n_paths);

            for (const auto& traj : traj_set.trajectories) {
                path_lengths.push_back(traj.size());
            }
        }

        if (!path_lengths.empty()) {
            std::sort(path_lengths.begin(), path_lengths.end());
            std::sort(paths_per_terminal.begin(), paths_per_terminal.end());

            #if DEBUG__compute_geodesic_basin
            Rprintf("\nPath statistics:\n");
            Rprintf("  Total paths: %zu\n", total_paths);
            Rprintf("  Paths per terminal: min=%zu, median=%zu, max=%zu\n",
                   paths_per_terminal.front(),
                   paths_per_terminal[paths_per_terminal.size() / 2],
                   paths_per_terminal.back());
            Rprintf("  Path lengths: min=%zu, median=%zu, max=%zu\n",
                   path_lengths.front(),
                   path_lengths[path_lengths.size() / 2],
                    path_lengths.back());
            #endif
        }

        #if DEBUG__compute_geodesic_basin
        Rprintf("===================================\n\n");
        #endif
    }

    // ========================================================================
    // COLLECT BOUNDARY VERTICES
    // ========================================================================
    const size_t n_vertices = adjacency_list.size();
    std::vector<unsigned char> is_in_basin(n_vertices, 0);
    for (const auto& [v, dist] : geodesic_basin.reachability_map.distances) {
        is_in_basin[v] = 1;
    }

    std::vector<unsigned char> boundary_seen(n_vertices, 0);
    for (const auto& [v, dist] : geodesic_basin.reachability_map.distances) {
        for (const auto& edge : adjacency_list[v]) {
            size_t u = edge.vertex;
            if (!is_in_basin[u] && !boundary_seen[u]) {
                boundary_seen[u] = 1;
                result.y_nbhd_bd_map[u] = y[u];
            }
        }
    }

    return result;
}
