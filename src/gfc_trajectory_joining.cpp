/**
 * @file gfc_trajectory_joining.cpp
 * @brief Implementation of trajectory joining functions for GFC
 *
 * This file contains the implementations of functions for joining gradient
 * trajectories across spurious extrema to create complete min-to-max paths.
 *
 */

#include "gfc.hpp"
#include "set_wgraph.hpp"

// ============================================================================
// Trajectory Helper Functions
// ============================================================================

/**
 * @brief Find shortest path from a vertex to the nearest target vertex.
 *
 * Uses Dijkstra's algorithm to find the closest target and the path to it.
 *
 * @param graph The weighted graph
 * @param source Starting vertex
 * @param targets Set of target vertices
 * @return Pair of (target_vertex, path from source to target), or (MAX, empty) if unreachable
 */
std::pair<size_t, std::vector<size_t>> find_shortest_path_to_nearest(
    const set_wgraph_t& graph,
    size_t source,
    const std::unordered_set<size_t>& targets
) {
    const size_t n = graph.num_vertices();
    const size_t INF_VERTEX = std::numeric_limits<size_t>::max();
    const double INF_DIST = std::numeric_limits<double>::infinity();

    if (targets.empty()) {
        return {INF_VERTEX, {}};
    }

    // Check if source is already a target
    if (targets.count(source)) {
        return {source, {source}};
    }

    // Dijkstra's algorithm
    std::vector<double> dist(n, INF_DIST);
    std::vector<size_t> parent(n, INF_VERTEX);
    dist[source] = 0.0;

    // Min-heap: (distance, vertex)
    using pq_entry = std::pair<double, size_t>;
    std::priority_queue<pq_entry, std::vector<pq_entry>, std::greater<pq_entry>> pq;
    pq.push({0.0, source});

    size_t nearest_target = INF_VERTEX;
    double nearest_dist = INF_DIST;

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        // Skip if we've already found a better path
        if (d > dist[u]) continue;

        // Stop early if we've reached the nearest target
        if (d >= nearest_dist) break;

        // Check if u is a target
        if (targets.count(u) && d < nearest_dist) {
            nearest_target = u;
            nearest_dist = d;
            continue;  // Don't expand further from targets
        }

        // Expand neighbors
        for (const auto& edge : graph.adjacency_list[u]) {
            size_t v = edge.vertex;
            double new_dist = d + edge.weight;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                parent[v] = u;
                pq.push({new_dist, v});
            }
        }
    }

    if (nearest_target == INF_VERTEX) {
        return {INF_VERTEX, {}};
    }

    // Reconstruct path from source to nearest_target
    std::vector<size_t> path;
    size_t current = nearest_target;
    while (current != INF_VERTEX) {
        path.push_back(current);
        current = parent[current];
    }
    std::reverse(path.begin(), path.end());

    return {nearest_target, path};
}

/**
 * @brief Reverse a path (create a copy with reversed order)
 */
std::vector<size_t> reverse_path(const std::vector<size_t>& path) {
    std::vector<size_t> reversed(path.rbegin(), path.rend());
    return reversed;
}

/**
 * @brief Join two paths at their connection point
 *
 * Assumes path1 ends at the same vertex where path2 starts.
 * The connection vertex appears only once in the result.
 *
 * @param path1 First path (ends at connection point)
 * @param path2 Second path (starts at connection point)
 * @return Joined path
 */
std::vector<size_t> join_paths(
    const std::vector<size_t>& path1,
    const std::vector<size_t>& path2
) {
    if (path1.empty()) return path2;
    if (path2.empty()) return path1;

    std::vector<size_t> joined;
    joined.reserve(path1.size() + path2.size() - 1);

    // Add all of path1
    joined.insert(joined.end(), path1.begin(), path1.end());

    // Add path2 starting from index 1 (skip connection point)
    if (path2.size() > 1) {
        joined.insert(joined.end(), path2.begin() + 1, path2.end());
    }

    return joined;
}

/**
 * @brief Compute sum of edge weights along a path
 */
double compute_path_length(
    const set_wgraph_t& graph,
    const std::vector<size_t>& path
) {
    if (path.size() < 2) return 0.0;

    double total = 0.0;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        total += graph.get_edge_weight(path[i], path[i + 1]);
    }

    return total;
}

/**
 * @brief Build cell map from joined trajectories
 *
 * Groups trajectory indices by their (min_vertex, max_vertex) pairs.
 */
std::map<std::pair<size_t, size_t>, std::vector<size_t>> build_cell_map(
    const std::vector<joined_trajectory_t>& joined_trajectories
) {
    std::map<std::pair<size_t, size_t>, std::vector<size_t>> cell_map;

    for (size_t i = 0; i < joined_trajectories.size(); ++i) {
        const auto& jt = joined_trajectories[i];
        auto key = std::make_pair(jt.min_vertex, jt.max_vertex);
        cell_map[key].push_back(i);
    }

    return cell_map;
}

// ============================================================================
// Extension Chain Finding (Recursive)
// ============================================================================

/**
 * @brief Find shortest extension chain to each reachable non-spurious target.
 *
 * Returns at most one chain per non-spurious target (the shortest found).
 * This avoids combinatorial explosion by not enumerating all possible paths.
 */
std::unordered_map<size_t, extension_chain_t> find_shortest_extensions(
    size_t spurious_vertex,
    bool is_min,
    const std::unordered_map<size_t, gradient_basin_t>& all_max_basins,
    const std::unordered_map<size_t, gradient_basin_t>& all_min_basins,
    const std::unordered_set<size_t>& non_spurious_targets,
    int depth_remaining,
    std::unordered_set<size_t> visited = {}
) {
    std::unordered_map<size_t, extension_chain_t> best_per_target;

    if (depth_remaining <= 0) return best_per_target;
    if (visited.count(spurious_vertex)) return best_per_target;

    visited.insert(spurious_vertex);

    const auto& basin_map = is_min ? all_min_basins : all_max_basins;

    auto it = basin_map.find(spurious_vertex);
    if (it == basin_map.end()) {
        return best_per_target;
    }

    const gradient_basin_t& basin = it->second;

    for (const auto& traj_set : basin.trajectory_sets) {
        size_t terminal = traj_set.terminal_vertex;

        // Find shortest trajectory in this set (by vertex count)
        const std::vector<size_t>* shortest_path = nullptr;
        for (const auto& path : traj_set.trajectories) {
            if (!shortest_path || path.size() < shortest_path->size()) {
                shortest_path = &path;
            }
        }

        if (!shortest_path || shortest_path->empty()) continue;

        if (non_spurious_targets.count(terminal)) {
            // Found non-spurious target directly
            extension_chain_t chain;
            chain.target_vertex = terminal;
            chain.extension_path = *shortest_path;
            chain.intermediates.push_back(spurious_vertex);

            auto existing = best_per_target.find(terminal);
            if (existing == best_per_target.end() ||
                chain.extension_path.size() < existing->second.extension_path.size()) {
                best_per_target[terminal] = std::move(chain);
            }
        } else {
            // Terminal is also spurious - recurse
            bool terminal_is_min = !is_min;

            auto sub_results = find_shortest_extensions(
                terminal,
                terminal_is_min,
                all_max_basins,
                all_min_basins,
                non_spurious_targets,
                depth_remaining - 1,
                visited
            );

            // Extend each sub-result through current spurious vertex
            for (auto& [target, sub_chain] : sub_results) {
                std::vector<size_t> full_path = join_paths(sub_chain.extension_path, *shortest_path);

                extension_chain_t extended;
                extended.target_vertex = target;
                extended.extension_path = std::move(full_path);
                extended.intermediates = sub_chain.intermediates;
                extended.intermediates.push_back(spurious_vertex);

                auto existing = best_per_target.find(target);
                if (existing == best_per_target.end() ||
                    extended.extension_path.size() < existing->second.extension_path.size()) {
                    best_per_target[target] = std::move(extended);
                }
            }
        }
    }

    return best_per_target;
}

// ============================================================================
// Main Joined Trajectory Computation
// ============================================================================

/**
 * @brief Compute one representative (shortest) trajectory per Morse-Smale cell.
 *
 * For each (non-spurious min, non-spurious max) pair that has connecting
 * gradient trajectories, stores the shortest trajectory found.
 *
 * Extensions are restricted to non-spurious minima within each maximum's basin,
 * ensuring trajectories stay within geometrically connected regions.
 */
std::vector<joined_trajectory_t> compute_joined_trajectories(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const std::unordered_map<size_t, gradient_basin_t>& all_max_basins_full,
    const std::unordered_map<size_t, gradient_basin_t>& all_min_basins_full,
    const std::unordered_set<size_t>& non_spurious_max,
    const std::unordered_set<size_t>& non_spurious_min,
    int max_chain_depth,
    bool verbose
) {
    // Track best (shortest) trajectory per cell
    std::map<std::pair<size_t, size_t>, joined_trajectory_t> best_per_cell;

    size_t direct_count = 0;
    size_t chain_extended_count = 0;
    size_t geodesic_extended_count = 0;
    size_t geodesic_fallback_count = 0;
    size_t unreachable_count = 0;

    if (verbose) {
        Rprintf("  Computing representative trajectories (one per cell)...\n");
        Rprintf("  Non-spurious: %zu maxima, %zu minima\n",
                non_spurious_max.size(), non_spurious_min.size());
    }

    size_t max_idx = 0;
    for (size_t max_v : non_spurious_max) {
        ++max_idx;

        auto max_it = all_max_basins_full.find(max_v);
        if (max_it == all_max_basins_full.end()) continue;

        const gradient_basin_t& max_basin = max_it->second;

        // ====================================================================
        // Build set of vertices in this basin
        // ====================================================================
        std::unordered_set<size_t> basin_vertex_set;
        for (const auto& [v, hop] : max_basin.hop_dist_map) {
            basin_vertex_set.insert(v);
        }

        // ====================================================================
        // Filter non-spurious minima to only those within this basin
        // ====================================================================
        std::unordered_set<size_t> valid_target_min;
        for (size_t min_v : non_spurious_min) {
            if (basin_vertex_set.count(min_v)) {
                valid_target_min.insert(min_v);
            }
        }

        if (verbose) {
            Rprintf("  Processing max %zu/%zu (vertex %zu): %zu terminal sets, %zu valid target minima\n",
                    max_idx, non_spurious_max.size(), max_v + 1,
                    max_basin.trajectory_sets.size(), valid_target_min.size());
        }

        // Skip if no valid targets in this basin
        if (valid_target_min.empty()) {
            if (verbose) {
                Rprintf("    No valid target minima in basin, skipping\n");
            }
            continue;
        }

        // ====================================================================
        // Process each trajectory set (grouped by terminal)
        // ====================================================================
        for (const auto& traj_set : max_basin.trajectory_sets) {
            size_t terminal_min = traj_set.terminal_vertex;

            // Find shortest trajectory in this set
            const std::vector<size_t>* shortest_orig = nullptr;
            for (const auto& path : traj_set.trajectories) {
                if (!shortest_orig || path.size() < shortest_orig->size()) {
                    shortest_orig = &path;
                }
            }

            if (!shortest_orig || shortest_orig->empty()) continue;

            // ----------------------------------------------------------------
            // Case 1: Direct connection to valid non-spurious minimum
            // ----------------------------------------------------------------
            if (valid_target_min.count(terminal_min)) {
                auto key = std::make_pair(terminal_min, max_v);
                auto existing = best_per_cell.find(key);

                if (existing == best_per_cell.end() ||
                    shortest_orig->size() < existing->second.path.size()) {

                    joined_trajectory_t jt;
                    jt.min_vertex = terminal_min;
                    jt.max_vertex = max_v;
                    jt.path = *shortest_orig;
                    jt.intermediate_extrema = {};
                    jt.total_change = y[max_v] - y[terminal_min];
                    jt.path_length = compute_path_length(graph, jt.path);

                    best_per_cell[key] = std::move(jt);
                    ++direct_count;
                }
                continue;
            }

            // ----------------------------------------------------------------
            // Case 2: Boundary vertex (not a local minimum, no basin)
            // Use geodesic to nearest valid target
            // ----------------------------------------------------------------
            if (all_min_basins_full.find(terminal_min) == all_min_basins_full.end()) {
                auto [nearest_min, extension_path] = find_shortest_path_to_nearest(
                    graph, terminal_min, valid_target_min
                );

                if (nearest_min == std::numeric_limits<size_t>::max()) {
                    ++unreachable_count;
                    continue;
                }

                auto key = std::make_pair(nearest_min, max_v);
                std::vector<size_t> full_path = join_paths(extension_path, *shortest_orig);

                auto existing = best_per_cell.find(key);
                if (existing == best_per_cell.end() ||
                    full_path.size() < existing->second.path.size()) {

                    joined_trajectory_t jt;
                    jt.min_vertex = nearest_min;
                    jt.max_vertex = max_v;
                    jt.path = std::move(full_path);
                    jt.intermediate_extrema = {terminal_min};
                    jt.total_change = y[max_v] - y[nearest_min];
                    jt.path_length = compute_path_length(graph, jt.path);

                    best_per_cell[key] = std::move(jt);
                    ++geodesic_extended_count;
                }
                continue;
            }

            // ----------------------------------------------------------------
            // Case 3: Spurious local minimum - try chain extension first
            // ----------------------------------------------------------------
            auto extensions = find_shortest_extensions(
                terminal_min,
                true,  // terminal_min is a minimum
                all_max_basins_full,
                all_min_basins_full,
                valid_target_min,  // Only targets within this basin
                max_chain_depth,
                {}
            );

            if (!extensions.empty()) {
                // Successfully found chain extensions
                for (const auto& [target_min, ext_chain] : extensions) {
                    std::vector<size_t> full_path = join_paths(ext_chain.extension_path, *shortest_orig);

                    auto key = std::make_pair(target_min, max_v);
                    auto existing = best_per_cell.find(key);

                    if (existing == best_per_cell.end() ||
                        full_path.size() < existing->second.path.size()) {

                        joined_trajectory_t jt;
                        jt.min_vertex = target_min;
                        jt.max_vertex = max_v;
                        jt.path = std::move(full_path);
                        jt.intermediate_extrema = ext_chain.intermediates;
                        jt.total_change = y[max_v] - y[target_min];
                        jt.path_length = compute_path_length(graph, jt.path);

                        best_per_cell[key] = std::move(jt);
                        ++chain_extended_count;
                    }
                }
            } else {
                // Chain extension failed - fallback to geodesic
                if (verbose) {
                    Rprintf("    Warning: spurious min %zu could not reach any valid min via chain "
                            "(depth limit = %d), trying geodesic...\n",
                            terminal_min + 1, max_chain_depth);
                }

                auto [nearest_min, extension_path] = find_shortest_path_to_nearest(
                    graph, terminal_min, valid_target_min
                );

                if (nearest_min == std::numeric_limits<size_t>::max()) {
                    ++unreachable_count;
                    continue;
                }

                auto key = std::make_pair(nearest_min, max_v);
                std::vector<size_t> full_path = join_paths(extension_path, *shortest_orig);

                auto existing = best_per_cell.find(key);
                if (existing == best_per_cell.end() ||
                    full_path.size() < existing->second.path.size()) {

                    joined_trajectory_t jt;
                    jt.min_vertex = nearest_min;
                    jt.max_vertex = max_v;
                    jt.path = std::move(full_path);
                    jt.intermediate_extrema = {terminal_min};
                    jt.total_change = y[max_v] - y[nearest_min];
                    jt.path_length = compute_path_length(graph, jt.path);

                    best_per_cell[key] = std::move(jt);
                    ++geodesic_fallback_count;
                }
            }
        }
    }

    // Convert map to vector
    std::vector<joined_trajectory_t> result;
    result.reserve(best_per_cell.size());
    for (auto& [key, jt] : best_per_cell) {
        result.push_back(std::move(jt));
    }

    if (verbose) {
        Rprintf("  Trajectory joining complete: %zu cells\n", result.size());
        Rprintf("    Direct: %zu, Chain-extended: %zu, Geodesic-extended: %zu\n",
                direct_count, chain_extended_count, geodesic_extended_count);
        Rprintf("    Geodesic-fallback: %zu, Unreachable: %zu\n",
                geodesic_fallback_count, unreachable_count);
    }

    return result;
}

/**
 * @brief Diagnostic function to trace extension search for a specific vertex
 */
void debug_extension_search(
    size_t query_vertex,
    bool query_is_min,
    const std::unordered_map<size_t, gradient_basin_t>& all_max_basins,
    const std::unordered_map<size_t, gradient_basin_t>& all_min_basins,
    const std::unordered_set<size_t>& non_spurious_max,
    const std::unordered_set<size_t>& non_spurious_min,
    int max_depth
) {
    Rprintf("\n=== DEBUG: Extension search for vertex %zu (%s) ===\n",
            query_vertex + 1, query_is_min ? "min" : "max");

    // Check if vertex exists in appropriate basin map
    const auto& basin_map = query_is_min ? all_min_basins : all_max_basins;
    const auto& target_set = query_is_min ? non_spurious_min : non_spurious_max;

    auto it = basin_map.find(query_vertex);
    if (it == basin_map.end()) {
        Rprintf("  ERROR: Vertex %zu has NO BASIN in all_%s_basins!\n",
                query_vertex + 1, query_is_min ? "min" : "max");
        Rprintf("  This means it was not identified as a local %s.\n",
                query_is_min ? "minimum" : "maximum");

        // Check the other basin map
        const auto& other_map = query_is_min ? all_max_basins : all_min_basins;
        if (other_map.count(query_vertex)) {
            Rprintf("  NOTE: Vertex %zu IS in all_%s_basins (opposite type).\n",
                    query_vertex + 1, query_is_min ? "max" : "min");
        }
        return;
    }

    const gradient_basin_t& basin = it->second;

    Rprintf("  Basin found: %zu trajectory sets, %zu terminal extrema\n",
            basin.trajectory_sets.size(), basin.terminal_extrema.size());

    // Check if query vertex is in non-spurious set
    if (target_set.count(query_vertex)) {
        Rprintf("  NOTE: Vertex %zu is already NON-SPURIOUS!\n", query_vertex + 1);
    }

    // Examine each terminal
    Rprintf("\n  Terminal analysis:\n");
    for (size_t i = 0; i < basin.trajectory_sets.size(); ++i) {
        const auto& traj_set = basin.trajectory_sets[i];
        size_t terminal = traj_set.terminal_vertex;

        // Terminals in a MIN basin are maxima (trajectories go FROM max TO min)
        bool terminal_is_max = query_is_min;
        const auto& terminal_target_set = terminal_is_max ? non_spurious_max : non_spurious_min;

        bool is_non_spurious = terminal_target_set.count(terminal) > 0;
        const char* status = is_non_spurious ? "NON-SPURIOUS" : "spurious";

        Rprintf("    Terminal %zu: vertex %zu (%s %s), %zu trajectories\n",
                i + 1, terminal + 1,
                status,
                terminal_is_max ? "MAX" : "MIN",
                traj_set.trajectories.size());

        // If spurious, check if it has a basin
        if (!is_non_spurious) {
            const auto& next_basin_map = terminal_is_max ? all_max_basins : all_min_basins;
            if (next_basin_map.count(terminal)) {
                const auto& next_basin = next_basin_map.at(terminal);
                Rprintf("      -> Has basin with %zu trajectory sets\n",
                        next_basin.trajectory_sets.size());
            } else {
                Rprintf("      -> NO BASIN (cannot extend further)\n");
            }
        }
    }

    // Now do a traced recursive search
    Rprintf("\n  Recursive trace (depth limit = %d):\n", max_depth);

    std::function<void(size_t, bool, int, std::string)> trace_search;
    trace_search = [&](size_t vertex, bool is_min, int depth, std::string indent) {
        if (depth <= 0) {
            Rprintf("%s[DEPTH LIMIT REACHED]\n", indent.c_str());
            return;
        }

        const auto& bmap = is_min ? all_min_basins : all_max_basins;
        auto bit = bmap.find(vertex);
        if (bit == bmap.end()) {
            Rprintf("%s[NO BASIN for %zu]\n", indent.c_str(), vertex + 1);
            return;
        }

        const auto& b = bit->second;
        for (const auto& ts : b.trajectory_sets) {
            size_t term = ts.terminal_vertex;
            bool term_is_max = is_min;
            const auto& term_targets = term_is_max ? non_spurious_max : non_spurious_min;

            if (term_targets.count(term)) {
                Rprintf("%sTerminal %zu: NON-SPURIOUS %s -> SUCCESS PATH FOUND!\n",
                        indent.c_str(), term + 1, term_is_max ? "MAX" : "MIN");
            } else {
                Rprintf("%sTerminal %zu: spurious %s -> recursing...\n",
                        indent.c_str(), term + 1, term_is_max ? "max" : "min");
                trace_search(term, !is_min, depth - 1, indent + "  ");
            }
        }
    };

    trace_search(query_vertex, query_is_min, max_depth, "    ");

    Rprintf("\n=== END DEBUG ===\n\n");
}
