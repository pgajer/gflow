/**
 * @file gfc_flow.cpp
 * @brief Implementation of trajectory-based Gradient Flow Complex computation
 *
 * This file implements the trajectory-first approach to GFC computation,
 * including member functions of set_wgraph_t for trajectory following
 * and free functions for the main computation pipeline.
 *
 * The trajectory-first approach traces gradient flow lines as the primary
 * operation, with basins emerging as the collection of vertices along
 * trajectories sharing the same endpoints.
 */

#include "gfc_flow.hpp"
#include "gfc.hpp"
#include "set_wgraph.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>

// ============================================================================
// Member Function Implementations for set_wgraph_t
// ============================================================================

int set_wgraph_t::check_nbr_extremum_type(
    size_t vertex,
    const std::vector<double>& y
) const {
    if (vertex >= adjacency_list.size()) {
        return 0;
    }

    const auto& neighbors = adjacency_list[vertex];
    if (neighbors.empty()) {
        // Isolated vertex: could be considered both min and max
        // Return 0 to indicate not a meaningful extremum
        return 0;
    }

    bool is_max = true;
    bool is_min = true;
    double y_v = y[vertex];

    for (const auto& edge : neighbors) {
        double y_u = y[edge.vertex];
        if (y_u >= y_v) {
            is_max = false;
        }
        if (y_u <= y_v) {
            is_min = false;
        }
        // Early exit if neither
        if (!is_max && !is_min) {
            return 0;
        }
    }

    if (is_max) return 1;
    if (is_min) return -1;
    return 0;
}

std::vector<size_t> set_wgraph_t::follow_gradient_trajectory(
    size_t start_vertex,
    const std::vector<double>& y,
    bool ascending,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    double edge_length_thld,
    const std::vector<bool>& visited,
    size_t max_length
) const {
    std::vector<size_t> trajectory;
    const size_t n = adjacency_list.size();

    // Validate inputs
    if (start_vertex >= n || y.size() != n) {
        return trajectory;
    }

    // Check if density is required but not provided
    bool needs_density = (modulation == gflow_modulation_t::DENSITY ||
                          modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_density && density.size() != n) {
        // Fall back to NONE or EDGELEN
        if (modulation == gflow_modulation_t::DENSITY) {
            modulation = gflow_modulation_t::NONE;
        } else {
            modulation = gflow_modulation_t::EDGELEN;
        }
    }

    // Track vertices in this trajectory to avoid self-loops
    std::unordered_set<size_t> trajectory_set;

    trajectory.push_back(start_vertex);
    trajectory_set.insert(start_vertex);

    size_t current = start_vertex;

    while (trajectory.size() < max_length) {
        // Find the best neighbor to move to
        size_t best_next = INVALID_VERTEX;
        double best_score = ascending ? -std::numeric_limits<double>::infinity()
                                      : std::numeric_limits<double>::infinity();

        for (const auto& edge : adjacency_list[current]) {
            size_t u = edge.vertex;
            double edge_len = edge.weight;

            // Skip if edge is too long
            if (edge_len > edge_length_thld) {
                continue;
            }

            // Skip if already in this trajectory (avoid loops)
            if (trajectory_set.count(u) > 0) {
                continue;
            }

            double delta_y = y[u] - y[current];

            // Skip if not moving in desired direction
            if (ascending && delta_y <= 0) {
                continue;
            }
            if (!ascending && delta_y >= 0) {
                continue;
            }

            // Get density for modulation
            double target_dens = needs_density ? density[u] : 1.0;

            // Compute modulated score
            double score = compute_modulated_score(delta_y, edge_len, target_dens, modulation);

            // Check if this is better
            bool is_better = ascending ? (score > best_score) : (score < best_score);
            if (is_better) {
                best_score = score;
                best_next = u;
            }
        }

        // If no valid neighbor found, we've reached a local extremum
        if (best_next == INVALID_VERTEX) {
            break;
        }

        // If best neighbor is already visited (by another trajectory), stop here
        // but include it as the endpoint
        if (visited[best_next]) {
            trajectory.push_back(best_next);
            break;
        }

        // Move to best neighbor
        trajectory.push_back(best_next);
        trajectory_set.insert(best_next);
        current = best_next;
    }

    return trajectory;
}

gflow_trajectory_t set_wgraph_t::join_trajectories_at_vertex(
    size_t vertex,
    const std::vector<double>& y,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    double edge_length_thld,
    const std::vector<bool>& visited,
    size_t max_length
) const {
    gflow_trajectory_t result;
    const size_t n = adjacency_list.size();

    if (vertex >= n || visited[vertex]) {
        return result;  // Empty trajectory
    }

    // Trace descending trajectory (toward local minimum)
    std::vector<size_t> desc_traj = follow_gradient_trajectory(
        vertex, y, false, modulation, density, edge_length_thld, visited, max_length
    );

    // Trace ascending trajectory (toward local maximum)
    std::vector<size_t> asc_traj = follow_gradient_trajectory(
        vertex, y, true, modulation, density, edge_length_thld, visited, max_length
    );

    // Join: desc_traj is ordered [vertex, ..., min_endpoint]
    //       asc_traj is ordered [vertex, ..., max_endpoint]
    // We want [min_endpoint, ..., vertex, ..., max_endpoint]

    // Reverse descending trajectory to get [min_endpoint, ..., vertex]
    std::reverse(desc_traj.begin(), desc_traj.end());

    // Build joined trajectory
    // desc_traj ends with vertex, asc_traj starts with vertex
    // So we concatenate desc_traj + asc_traj[1:]
    result.vertices = std::move(desc_traj);
    for (size_t i = 1; i < asc_traj.size(); ++i) {
        result.vertices.push_back(asc_traj[i]);
    }

    if (result.vertices.empty()) {
        return result;
    }

    // Set trajectory properties
    result.start_vertex = result.vertices.front();
    result.end_vertex = result.vertices.back();

    // Check if endpoints are actual local extrema
    result.starts_at_lmin = (check_nbr_extremum_type(result.start_vertex, y) == -1);
    result.ends_at_lmax = (check_nbr_extremum_type(result.end_vertex, y) == 1);

    // Compute total change
    result.total_change = y[result.end_vertex] - y[result.start_vertex];

    return result;
}

// ============================================================================
// Free Function Implementations
// ============================================================================

gfc_flow_result_t compute_gfc_flow(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density,
    bool verbose
) {
    gfc_flow_result_t result;
    result.n_vertices = graph.num_vertices();
    result.params = params;

    const size_t n = result.n_vertices;

    if (n == 0 || y.size() != n) {
        return result;
    }

    // Compute median of y
    std::vector<double> y_sorted = y;
    std::sort(y_sorted.begin(), y_sorted.end());
    result.y_median = y_sorted[n / 2];

    if (verbose) {
        Rprintf("GFC Flow computation: %zu vertices, median = %.4f\n", n, result.y_median);
        Rprintf("  Modulation: %s\n", gflow_modulation_to_string(params.modulation).c_str());
    }

    // ========================================================================
    // Step 1: Find local extrema and compute edge length threshold
    // ========================================================================

    if (verbose) {
        Rprintf("Step 1: Finding local extrema and preparing trajectory computation...\n");
    }

    auto [lmin_vertices, lmax_vertices] = graph.find_nbr_extrema(y);

    // Create sets for fast lookup
    std::unordered_set<size_t> lmin_set(lmin_vertices.begin(), lmin_vertices.end());
    std::unordered_set<size_t> lmax_set(lmax_vertices.begin(), lmax_vertices.end());

    // Compute edge length threshold
    double edge_length_thld = graph.compute_quantile_edge_length(params.edge_length_quantile_thld);

    if (verbose) {
        Rprintf("  Found %zu local minima, %zu local maxima\n",
                lmin_vertices.size(), lmax_vertices.size());
        Rprintf("  Edge length threshold (%.2f quantile): %.4f\n",
                params.edge_length_quantile_thld, edge_length_thld);
    }

    // ========================================================================
    // Step 2: Trace trajectories from all local minima
    // ========================================================================

    if (verbose) {
        Rprintf("Step 2: Tracing ascending trajectories from local minima...\n");
    }

    std::vector<bool> visited(n, false);
    std::vector<int> vertex_to_traj(n, -1);

    // Maps to accumulate basin memberships
    // Key: extremum vertex, Value: set of vertices in its basin
    std::map<size_t, std::set<size_t>> ascending_basins;   // lmin -> vertices
    std::map<size_t, std::set<size_t>> descending_basins;  // lmax -> vertices

    int trajectory_id = 0;
    result.n_lmin_trajectories = 0;

    // Process local minima in order of increasing y value for determinism
    std::vector<size_t> lmin_sorted = lmin_vertices;
    std::sort(lmin_sorted.begin(), lmin_sorted.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });

    for (size_t lmin : lmin_sorted) {
        // Trace ascending trajectory from this local minimum
        std::vector<size_t> traj = graph.follow_gradient_trajectory(
            lmin, y, true,
            params.modulation, density,
            edge_length_thld, visited,
            params.max_trajectory_length
        );

        if (traj.empty()) {
            continue;
        }

        // The endpoint should be a local maximum (or a previously visited vertex)
        size_t lmax = traj.back();

        // Assign all trajectory vertices to basins
        for (size_t v : traj) {
            ascending_basins[lmin].insert(v);
            descending_basins[lmax].insert(v);

            if (!visited[v]) {
                visited[v] = true;
                vertex_to_traj[v] = trajectory_id;
            }
        }

        // Store trajectory if requested
        if (params.store_trajectories) {
            gflow_trajectory_t gft;
            gft.vertices = std::move(traj);
            gft.start_vertex = lmin;
            gft.end_vertex = lmax;
            gft.starts_at_lmin = true;
            gft.ends_at_lmax = (lmax_set.count(lmax) > 0);
            gft.total_change = y[lmax] - y[lmin];
            gft.trajectory_id = trajectory_id;
            result.trajectories.push_back(std::move(gft));
        }

        ++trajectory_id;
        ++result.n_lmin_trajectories;
    }

    if (verbose) {
        size_t n_visited = std::count(visited.begin(), visited.end(), true);
        Rprintf("  Traced %d trajectories, covering %zu/%zu vertices (%.1f%%)\n",
                result.n_lmin_trajectories, n_visited, n,
                100.0 * n_visited / n);
    }

    // ========================================================================
    // Step 3: Handle unvisited vertices by joining trajectories
    // ========================================================================

    if (verbose) {
        Rprintf("Step 3: Processing unvisited vertices via trajectory joining...\n");
    }

    result.n_join_trajectories = 0;

    // Find unvisited vertices
    std::vector<size_t> unvisited;
    for (size_t v = 0; v < n; ++v) {
        if (!visited[v]) {
            unvisited.push_back(v);
        }
    }

    // Process unvisited vertices
    // Sort by y value for deterministic processing
    std::sort(unvisited.begin(), unvisited.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });

    for (size_t v : unvisited) {
        if (visited[v]) {
            // May have been visited by a previous join
            continue;
        }

        gflow_trajectory_t joined = graph.join_trajectories_at_vertex(
            v, y, params.modulation, density,
            edge_length_thld, visited,
            params.max_trajectory_length
        );

        if (joined.vertices.empty()) {
            // Isolated vertex or error: mark as visited anyway
            visited[v] = true;
            continue;
        }

        // Assign trajectory vertices to basins
        size_t lmin = joined.start_vertex;
        size_t lmax = joined.end_vertex;

        for (size_t u : joined.vertices) {
            ascending_basins[lmin].insert(u);
            descending_basins[lmax].insert(u);

            if (!visited[u]) {
                visited[u] = true;
                vertex_to_traj[u] = trajectory_id;
            }
        }

        // Store trajectory if requested
        if (params.store_trajectories) {
            joined.trajectory_id = trajectory_id;
            result.trajectories.push_back(std::move(joined));
        }

        ++trajectory_id;
        ++result.n_join_trajectories;
    }

    if (verbose) {
        size_t n_visited = std::count(visited.begin(), visited.end(), true);
        Rprintf("  Joined %d additional trajectories, total coverage: %zu/%zu (%.1f%%)\n",
                result.n_join_trajectories, n_visited, n,
                100.0 * n_visited / n);
    }

    result.vertex_trajectory = std::move(vertex_to_traj);

    // ========================================================================
    // Step 4: Convert basin maps to basin_compact_t structures
    // ========================================================================

    if (verbose) {
        Rprintf("Step 4: Building basin structures...\n");
    }

    std::vector<basin_compact_t> min_basins;
    std::vector<basin_compact_t> max_basins;

    // Build minimum basins (ascending basins keyed by local minima)
    for (const auto& [lmin, vertices] : ascending_basins) {
        basin_compact_t basin;
        basin.extremum_vertex = lmin;
        basin.extremum_value = y[lmin];
        basin.is_maximum = false;

        basin.vertices.assign(vertices.begin(), vertices.end());

        // Compute hop distances from extremum (simplified: using BFS)
        basin.hop_distances.resize(basin.vertices.size(), 0);
        std::unordered_map<size_t, int> hop_map;
        hop_map[lmin] = 0;

        std::queue<size_t> bfs_queue;
        bfs_queue.push(lmin);
        int max_hop = 0;

        while (!bfs_queue.empty()) {
            size_t u = bfs_queue.front();
            bfs_queue.pop();
            int u_hop = hop_map[u];

            for (const auto& edge : graph.adjacency_list[u]) {
                size_t v = edge.vertex;
                if (vertices.count(v) && hop_map.find(v) == hop_map.end()) {
                    hop_map[v] = u_hop + 1;
                    max_hop = std::max(max_hop, u_hop + 1);
                    bfs_queue.push(v);
                }
            }
        }

        // Fill hop distances
        for (size_t i = 0; i < basin.vertices.size(); ++i) {
            auto it = hop_map.find(basin.vertices[i]);
            basin.hop_distances[i] = (it != hop_map.end()) ? it->second : -1;
        }
        basin.max_hop_distance = max_hop;

        min_basins.push_back(std::move(basin));
    }

    // Build maximum basins (descending basins keyed by local maxima)
    for (const auto& [lmax, vertices] : descending_basins) {
        basin_compact_t basin;
        basin.extremum_vertex = lmax;
        basin.extremum_value = y[lmax];
        basin.is_maximum = true;

        basin.vertices.assign(vertices.begin(), vertices.end());

        // Compute hop distances from extremum
        basin.hop_distances.resize(basin.vertices.size(), 0);
        std::unordered_map<size_t, int> hop_map;
        hop_map[lmax] = 0;

        std::queue<size_t> bfs_queue;
        bfs_queue.push(lmax);
        int max_hop = 0;

        while (!bfs_queue.empty()) {
            size_t u = bfs_queue.front();
            bfs_queue.pop();
            int u_hop = hop_map[u];

            for (const auto& edge : graph.adjacency_list[u]) {
                size_t v = edge.vertex;
                if (vertices.count(v) && hop_map.find(v) == hop_map.end()) {
                    hop_map[v] = u_hop + 1;
                    max_hop = std::max(max_hop, u_hop + 1);
                    bfs_queue.push(v);
                }
            }
        }

        for (size_t i = 0; i < basin.vertices.size(); ++i) {
            auto it = hop_map.find(basin.vertices[i]);
            basin.hop_distances[i] = (it != hop_map.end()) ? it->second : -1;
        }
        basin.max_hop_distance = max_hop;

        max_basins.push_back(std::move(basin));
    }

    // Compute initial summaries
    auto max_summaries = gfc_internal::compute_basin_summaries(
        graph, max_basins, y, result.y_median, params.hop_k
    );
    auto min_summaries = gfc_internal::compute_basin_summaries(
        graph, min_basins, y, result.y_median, params.hop_k
    );

    // Record initial counts
    result.stage_history.emplace_back(
        "initial",
        static_cast<int>(max_basins.size()),
        static_cast<int>(max_basins.size()),
        static_cast<int>(min_basins.size()),
        static_cast<int>(min_basins.size())
    );

    if (verbose) {
        Rprintf("  Built %zu min basins, %zu max basins\n",
                min_basins.size(), max_basins.size());
    }

    // ========================================================================
    // Step 5: Apply filtering pipeline (reuse from gfc_internal)
    // ========================================================================

    // Step 5a: Filter by relative values
    if (params.apply_relvalue_filter) {
        if (verbose) {
            Rprintf("Step 5a: Filtering by relative values (max >= %.2f, min <= %.2f)...\n",
                    params.min_rel_value_max, params.max_rel_value_min);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::filter_by_relvalue(
            max_summaries, max_basins,
            params.min_rel_value_max, params.max_rel_value_min, true
        );
        gfc_internal::filter_by_relvalue(
            min_summaries, min_basins,
            params.min_rel_value_max, params.max_rel_value_min, false
        );

        result.stage_history.emplace_back(
            "relvalue",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Retained %zu maxima and %zu minima\n",
                    max_basins.size(), min_basins.size());
        }
    }

    // Step 5b: Cluster and merge maxima
    if (params.apply_maxima_clustering && max_basins.size() > 1) {
        if (verbose) {
            Rprintf("Step 5b: Clustering maxima (overlap threshold = %.2f)...\n",
                    params.max_overlap_threshold);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::cluster_and_merge_basins(
            max_summaries, max_basins,
            params.max_overlap_threshold, true
        );

        result.stage_history.emplace_back(
            "merge.max",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Result: %zu maxima after merging\n", max_basins.size());
        }
    }

    // Step 5c: Cluster and merge minima
    if (params.apply_minima_clustering && min_basins.size() > 1) {
        if (verbose) {
            Rprintf("Step 5c: Clustering minima (overlap threshold = %.2f)...\n",
                    params.min_overlap_threshold);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::cluster_and_merge_basins(
            min_summaries, min_basins,
            params.min_overlap_threshold, false
        );

        result.stage_history.emplace_back(
            "merge.min",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Result: %zu minima after merging\n", min_basins.size());
        }
    }

    // Step 5d: Filter by geometric characteristics
    if (params.apply_geometric_filter) {
        if (verbose) {
            Rprintf("Step 5d: Geometric filtering (nbrs < %.2f, hopk < %.2f, deg < %.2f, size >= %d)...\n",
                    params.p_mean_nbrs_dist_threshold,
                    params.p_mean_hopk_dist_threshold,
                    params.p_deg_threshold,
                    params.min_basin_size);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        // Recompute summaries after merging
        max_summaries = gfc_internal::compute_basin_summaries(
            graph, max_basins, y, result.y_median, params.hop_k
        );
        min_summaries = gfc_internal::compute_basin_summaries(
            graph, min_basins, y, result.y_median, params.hop_k
        );

        gfc_internal::filter_by_geometry(
            max_summaries,
            max_basins,
            params.p_mean_nbrs_dist_threshold,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size
        );

        gfc_internal::filter_by_geometry(
            min_summaries,
            min_basins,
            params.p_mean_nbrs_dist_threshold,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size
        );

        result.stage_history.emplace_back(
            "geometric",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Retained %zu maxima and %zu minima\n",
                    max_basins.size(), min_basins.size());
        }
    }

    // ========================================================================
    // Step 6: Build membership vectors and assignments
    // ========================================================================

    result.max_membership = gfc_internal::build_membership_vectors(max_basins, n);
    result.min_membership = gfc_internal::build_membership_vectors(min_basins, n);

    // Build assignment vectors (first basin for each vertex)
    result.max_assignment.resize(n, -1);
    result.min_assignment.resize(n, -1);

    for (size_t v = 0; v < n; ++v) {
        if (!result.max_membership[v].empty()) {
            result.max_assignment[v] = result.max_membership[v][0];
        }
        if (!result.min_membership[v].empty()) {
            result.min_assignment[v] = result.min_membership[v][0];
        }
    }

    // ========================================================================
    // Finalize result
    // ========================================================================

    result.max_basins = std::move(max_basins);
    result.min_basins = std::move(min_basins);
    result.max_summaries = std::move(max_summaries);
    result.min_summaries = std::move(min_summaries);

    if (verbose) {
        Rprintf("GFC Flow computation complete: %zu maxima, %zu minima, %zu trajectories\n",
                result.max_basins.size(), result.min_basins.size(),
                result.trajectories.size());
    }

    return result;
}

std::vector<gfc_flow_result_t> compute_gfc_flow_matrix(
    const set_wgraph_t& graph,
    const Eigen::MatrixXd& Y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density,
    int n_cores,
    bool verbose
) {
    const int n = static_cast<int>(Y.rows());
    const int p = static_cast<int>(Y.cols());

    if (verbose) {
        Rprintf("GFC Flow matrix computation: %d functions, %d vertices\n", p, n);
#ifdef _OPENMP
        Rprintf("  OpenMP threads: %d\n", n_cores);
#else
        Rprintf("  OpenMP: not available (sequential execution)\n");
#endif
    }

    std::vector<gfc_flow_result_t> results(p);

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    int progress_interval = std::max(1, p / 20);

    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
    for (int j = 0; j < p; ++j) {
        std::vector<double> y_j(n);
        for (int i = 0; i < n; ++i) {
            y_j[i] = Y(i, j);
        }

        results[j] = compute_gfc_flow(graph, y_j, params, density, false);

        #pragma omp critical
        {
            if (verbose && n_cores == 1 && (j + 1) % progress_interval == 0) {
                Rprintf("  Processed %d/%d functions (%.0f%%)\n",
                        j + 1, p, 100.0 * (j + 1) / p);
            }
        }
    }

    if (verbose) {
        Rprintf("  Complete.\n");
    }

    return results;
}
