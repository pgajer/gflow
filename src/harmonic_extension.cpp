/**
 * @file harmonic_extension.cpp
 * @brief Implementation of harmonic extension for trajectory coordinates
 */

#include "harmonic_extension.hpp"
#include "set_wgraph.hpp"
#include <queue>
#include <cmath>
#include <algorithm>
#include <limits>
#include <R_ext/Print.h>

/**
 * @brief Compute arc-length coordinates for trajectory vertices
 *
 * @param trajectory Ordered vertices from minimum to maximum
 * @return Pair of (coordinates in [0,1], total trajectory length)
 */
std::pair<std::vector<double>, double>
set_wgraph_t::compute_arc_length_coords(
    const std::vector<size_t>& trajectory
) const {

    const size_t n = trajectory.size();
    std::vector<double> coords(n, 0.0);

    if (n <= 1) {
        return {coords, 0.0};
    }

    // Compute cumulative distances
    std::vector<double> cumulative(n, 0.0);

    for (size_t i = 1; i < n; ++i) {
        double edge_len = get_edge_weight(trajectory[i-1], trajectory[i]);
        cumulative[i] = cumulative[i-1] + edge_len;
    }

    double total_length = cumulative[n-1];

    // Normalize to [0, 1]
    if (total_length > 0) {
        for (size_t i = 0; i < n; ++i) {
            coords[i] = cumulative[i] / total_length;
        }
    }

    return {coords, total_length};
}

/**
 * @brief Solve harmonic extension using Gauss-Seidel iteration
 *
 * @param tubular_vertices Vertices in tubular neighborhood
 * @param trajectory_set Set of trajectory vertices (Dirichlet boundary)
 * @param trajectory_coords Map from trajectory vertex to its coordinate
 * @param initial_coords Initial coordinates for all tubular vertices
 *                       (should be trajectory coord of nearest trajectory vertex)
 * @param use_edge_weights Use inverse edge length weights
 * @param max_iterations Maximum number of iterations
 * @param tolerance Convergence tolerance
 * @param[out] n_iterations Number of iterations performed
 * @param[out] final_max_change Final maximum change
 * @return Extended coordinates for all tubular vertices
 */
std::vector<double> set_wgraph_t::solve_harmonic_extension(
    const std::vector<size_t>& tubular_vertices,
    const std::unordered_set<size_t>& trajectory_set,
    const std::unordered_map<size_t, double>& trajectory_coords,
    const std::vector<double>& initial_coords,
    bool use_edge_weights,
    int max_iterations,
    double tolerance,
    int& n_iterations,
    double& final_max_change
) const {

    const size_t n = tubular_vertices.size();

    // Build vertex to index map
    std::unordered_map<size_t, size_t> vertex_to_idx;
    for (size_t i = 0; i < n; ++i) {
        vertex_to_idx[tubular_vertices[i]] = i;
    }

    // Create set of tubular vertices for fast lookup
    std::unordered_set<size_t> tubular_set(
        tubular_vertices.begin(), tubular_vertices.end()
    );

    // Initialize coordinates from initial_coords or trajectory values
    std::vector<double> coords(n);
    std::vector<bool> is_fixed(n, false);

    for (size_t i = 0; i < n; ++i) {
        size_t v = tubular_vertices[i];

        if (trajectory_set.count(v) > 0) {
            // Fixed Dirichlet boundary - use exact trajectory coordinate
            coords[i] = trajectory_coords.at(v);
            is_fixed[i] = true;
        } else {
            // Use provided initial coordinate (nearest trajectory vertex)
            coords[i] = initial_coords[i];
        }
    }

    // Precompute neighbor weights for each vertex
    std::vector<std::vector<std::pair<size_t, double>>> neighbor_weights(n);

    for (size_t i = 0; i < n; ++i) {
        size_t v = tubular_vertices[i];

        for (const auto& edge : adjacency_list[v]) {
            size_t w = edge.vertex;

            // Only include neighbors in tubular neighborhood
            if (tubular_set.count(w) == 0) continue;

            size_t j = vertex_to_idx.at(w);
            double weight = use_edge_weights ? (1.0 / edge.weight) : 1.0;

            neighbor_weights[i].push_back({j, weight});
        }
    }

    // Gauss-Seidel iteration
    n_iterations = 0;
    final_max_change = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_change = 0.0;

        for (size_t i = 0; i < n; ++i) {
            if (is_fixed[i]) continue;

            // Compute weighted average of neighbors
            double weighted_sum = 0.0;
            double weight_total = 0.0;

            for (const auto& [j, w] : neighbor_weights[i]) {
                weighted_sum += w * coords[j];
                weight_total += w;
            }

            if (weight_total > 0) {
                double new_coord = weighted_sum / weight_total;
                double change = std::abs(new_coord - coords[i]);
                max_change = std::max(max_change, change);
                coords[i] = new_coord;
            }
        }

        n_iterations = iter + 1;
        final_max_change = max_change;

        if (max_change < tolerance) {
            break;
        }
    }

    return coords;
}

/**
 * @brief Compute tubular neighborhood using hop distance (BFS)
 */
tubular_neighborhood_t set_wgraph_t::compute_tubular_neighborhood_hop(
    const std::vector<size_t>& trajectory,
    int hop_radius,
    const std::unordered_set<size_t>& basin_restriction
) const {

    tubular_neighborhood_t result;

    // Track visited vertices: vertex -> (hop dist, geodesic dist, traj idx)
    struct vertex_info {
        int hop_dist;
        double geo_dist;
        size_t traj_idx;
    };
    std::unordered_map<size_t, vertex_info> visited;

    // Build trajectory vertex to index map
    std::unordered_map<size_t, size_t> traj_vertex_to_idx;
    for (size_t i = 0; i < trajectory.size(); ++i) {
        traj_vertex_to_idx[trajectory[i]] = i;
    }

    // BFS queue: (vertex, hop distance, geodesic distance, source traj index)
    std::queue<std::tuple<size_t, int, double, size_t>> bfs_queue;

    // Initialize with trajectory vertices
    for (size_t i = 0; i < trajectory.size(); ++i) {
        size_t v = trajectory[i];
        if (basin_restriction.empty() || basin_restriction.count(v) > 0) {
            visited[v] = {0, 0.0, i};
            bfs_queue.push({v, 0, 0.0, i});
        }
    }

    // BFS expansion
    while (!bfs_queue.empty()) {
        auto [u, hop_dist, geo_dist, traj_idx] = bfs_queue.front();
        bfs_queue.pop();

        if (hop_dist < hop_radius) {
            for (const auto& edge : adjacency_list[u]) {
                size_t w = edge.vertex;

                if (!basin_restriction.empty() && basin_restriction.count(w) == 0) {
                    continue;
                }

                if (visited.count(w) > 0) {
                    continue;
                }

                // Check if w is on trajectory
                auto traj_it = traj_vertex_to_idx.find(w);
                if (traj_it != traj_vertex_to_idx.end()) {
                    visited[w] = {0, 0.0, traj_it->second};
                } else {
                    double new_geo_dist = geo_dist + edge.weight;
                    visited[w] = {hop_dist + 1, new_geo_dist, traj_idx};
                    bfs_queue.push({w, hop_dist + 1, new_geo_dist, traj_idx});
                }
            }
        }
    }

    // Collect results
    result.vertices.reserve(visited.size());
    result.hop_distances.reserve(visited.size());
    result.geodesic_distances.reserve(visited.size());
    result.nearest_traj_idx.reserve(visited.size());

    for (const auto& [v, info] : visited) {
        result.vertices.push_back(v);
        result.hop_distances.push_back(info.hop_dist);
        result.geodesic_distances.push_back(info.geo_dist);
        result.nearest_traj_idx.push_back(info.traj_idx);
    }

    return result;
}

/**
 * @brief Compute tubular neighborhood using geodesic distance (multi-source Dijkstra)
 */
tubular_neighborhood_t set_wgraph_t::compute_tubular_neighborhood_geodesic(
    const std::vector<size_t>& trajectory,
    double geodesic_radius,
    const std::unordered_set<size_t>& basin_restriction
) const {

    tubular_neighborhood_t result;

    const size_t n = adjacency_list.size();

    // Distance and predecessor tracking
    std::vector<double> geo_dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> hop_dist(n, std::numeric_limits<int>::max());
    std::vector<size_t> nearest_traj(n, INVALID_VERTEX);
    std::vector<bool> finalized(n, false);

    // Build trajectory vertex to index map
    std::unordered_map<size_t, size_t> traj_vertex_to_idx;
    for (size_t i = 0; i < trajectory.size(); ++i) {
        traj_vertex_to_idx[trajectory[i]] = i;
    }

    // Priority queue: (geodesic distance, hop distance, vertex, traj index)
    using pq_entry = std::tuple<double, int, size_t, size_t>;
    std::priority_queue<pq_entry, std::vector<pq_entry>, std::greater<pq_entry>> pq;

    // Initialize with trajectory vertices
    for (size_t i = 0; i < trajectory.size(); ++i) {
        size_t v = trajectory[i];
        if (basin_restriction.empty() || basin_restriction.count(v) > 0) {
            geo_dist[v] = 0.0;
            hop_dist[v] = 0;
            nearest_traj[v] = i;
            pq.push({0.0, 0, v, i});
        }
    }

    // Multi-source Dijkstra
    while (!pq.empty()) {
        auto [d, h, u, traj_idx] = pq.top();
        pq.pop();

        if (finalized[u]) {
            continue;
        }
        finalized[u] = true;

        // Stop if beyond geodesic radius (for this vertex)
        // But continue processing queue for other paths
        if (d > geodesic_radius) {
            continue;
        }

        // Explore neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t w = edge.vertex;

            if (finalized[w]) {
                continue;
            }

            if (!basin_restriction.empty() && basin_restriction.count(w) == 0) {
                continue;
            }

            double new_geo_dist = d + edge.weight;

            // Only add if within radius
            if (new_geo_dist > geodesic_radius) {
                continue;
            }

            // Check if w is on trajectory - use its own index
            auto traj_it = traj_vertex_to_idx.find(w);
            size_t w_traj_idx = (traj_it != traj_vertex_to_idx.end())
                                 ? traj_it->second : traj_idx;
            double w_geo_dist = (traj_it != traj_vertex_to_idx.end())
                                 ? 0.0 : new_geo_dist;
            int w_hop_dist = (traj_it != traj_vertex_to_idx.end())
                              ? 0 : h + 1;

            if (w_geo_dist < geo_dist[w]) {
                geo_dist[w] = w_geo_dist;
                hop_dist[w] = w_hop_dist;
                nearest_traj[w] = w_traj_idx;
                pq.push({w_geo_dist, w_hop_dist, w, w_traj_idx});
            }
        }
    }

    // Collect results - only vertices that were reached within radius
    for (size_t v = 0; v < n; ++v) {
        if (finalized[v] && geo_dist[v] <= geodesic_radius) {
            result.vertices.push_back(v);
            result.hop_distances.push_back(hop_dist[v]);
            result.geodesic_distances.push_back(geo_dist[v]);
            result.nearest_traj_idx.push_back(nearest_traj[v]);
        }
    }

    return result;
}

/**
 * @brief Compute tubular neighborhood (dispatcher)
 */
tubular_neighborhood_t set_wgraph_t::compute_tubular_neighborhood(
    const std::vector<size_t>& trajectory,
    double radius,
    tube_radius_type_t radius_type,
    const std::unordered_set<size_t>& basin_restriction
) const {

    if (radius_type == tube_radius_type_t::HOP) {
        return compute_tubular_neighborhood_hop(
            trajectory, static_cast<int>(radius), basin_restriction
        );
    } else {
        return compute_tubular_neighborhood_geodesic(
            trajectory, radius, basin_restriction
        );
    }
}

/**
 * @brief Compute harmonic extension of trajectory coordinates
 */
harmonic_extension_result_t set_wgraph_t::compute_harmonic_extension(
    const std::vector<size_t>& trajectory,
    const harmonic_extension_params_t& params,
    bool verbose
) const {

    harmonic_extension_result_t result;
    result.trajectory = trajectory;
    result.tube_type = params.tube_type;
    result.tube_radius = params.tube_radius;

    if (trajectory.empty()) {
        if (verbose) {
            Rprintf("Warning: Empty trajectory provided\n");
        }
        return result;
    }

    const char* type_str = (params.tube_type == tube_radius_type_t::HOP)
                            ? "hop" : "geodesic";

    if (verbose) {
        Rprintf("Computing harmonic extension for trajectory of %zu vertices\n",
                trajectory.size());
        Rprintf("  Tube radius: %.2f (%s)\n", params.tube_radius, type_str);
        Rprintf("  Laplacian weights: %s\n",
                params.use_edge_weights ? "inverse length" : "unit");
    }

    // Step 1: Compute arc-length coordinates for trajectory
    auto [traj_coords, total_length] = compute_arc_length_coords(trajectory);
    result.trajectory_coords = traj_coords;
    result.trajectory_length = total_length;

    if (verbose) {
        Rprintf("  Trajectory length: %.4f\n", total_length);
    }

    // Step 2: Compute tubular neighborhood
    auto tube_nbhd = compute_tubular_neighborhood(
        trajectory,
        params.tube_radius,
        params.tube_type,
        params.basin_restriction
    );

    result.tubular_vertices = tube_nbhd.vertices;
    result.hop_distances = tube_nbhd.hop_distances;
    result.geodesic_distances = tube_nbhd.geodesic_distances;
    result.nearest_traj_idx = tube_nbhd.nearest_traj_idx;

    if (verbose) {
        Rprintf("  Tubular neighborhood: %zu vertices\n", tube_nbhd.vertices.size());

        // Summary statistics
        if (!tube_nbhd.geodesic_distances.empty()) {
            double max_geo = *std::max_element(
                tube_nbhd.geodesic_distances.begin(),
                tube_nbhd.geodesic_distances.end()
            );
            int max_hop = *std::max_element(
                tube_nbhd.hop_distances.begin(),
                tube_nbhd.hop_distances.end()
            );
            Rprintf("    Max geodesic distance: %.4f\n", max_geo);
            Rprintf("    Max hop distance: %d\n", max_hop);
        }

        // Count by hop distance
        std::map<int, int> hop_counts;
        for (int h : tube_nbhd.hop_distances) {
            hop_counts[h]++;
        }
        for (const auto& [h, count] : hop_counts) {
            Rprintf("    Hop %d: %d vertices\n", h, count);
        }
    }

    // Build trajectory coordinate map for solver
    std::unordered_set<size_t> trajectory_set(trajectory.begin(), trajectory.end());
    std::unordered_map<size_t, double> trajectory_coord_map;

    for (size_t i = 0; i < trajectory.size(); ++i) {
        trajectory_coord_map[trajectory[i]] = traj_coords[i];
    }

    // Step 3: Build initial coordinates from nearest trajectory vertex
    std::vector<double> initial_coords(tube_nbhd.vertices.size());
    for (size_t i = 0; i < tube_nbhd.vertices.size(); ++i) {
        initial_coords[i] = traj_coords[tube_nbhd.nearest_traj_idx[i]];
    }

    if (verbose) {
        Rprintf("  Initial coords from nearest-neighbor projection\n");
    }

    // Step 4: Solve harmonic extension
    int n_iter;
    double max_change;

    result.extended_coords = solve_harmonic_extension(
        tube_nbhd.vertices,
        trajectory_set,
        trajectory_coord_map,
        initial_coords,
        params.use_edge_weights,
        params.max_iterations,
        params.tolerance,
        n_iter,
        max_change
    );

    result.n_iterations = n_iter;
    result.final_max_change = max_change;

    if (verbose) {
        Rprintf("  Solver converged in %d iterations (max change: %.2e)\n",
                n_iter, max_change);
    }

    // Build vertex to index map
    for (size_t i = 0; i < tube_nbhd.vertices.size(); ++i) {
        result.vertex_to_idx[tube_nbhd.vertices[i]] = i;
    }

    // Validate coordinates are in [0, 1]
    double min_coord = *std::min_element(
        result.extended_coords.begin(), result.extended_coords.end()
    );
    double max_coord = *std::max_element(
        result.extended_coords.begin(), result.extended_coords.end()
    );

    if (verbose) {
        Rprintf("  Coordinate range: [%.4f, %.4f]\n", min_coord, max_coord);
    }

    return result;
}

/**
 * @brief Extract the trajectory with maximal mean density from a cell
 *
 * Given all trajectories in a cell (connecting minimum to maximum),
 * selects the one with the highest mean vertex density.
 *
 * @param trajectories Vector of trajectories (each a vector of vertex indices)
 * @param density Density values for each vertex
 * @return Index of the trajectory with maximal mean density
 */
size_t select_max_density_trajectory(
    const std::vector<std::vector<size_t>>& trajectories,
    const std::vector<double>& density
) {
    if (trajectories.empty()) {
        return 0;
    }

    size_t best_idx = 0;
    double best_mean_density = -std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < trajectories.size(); ++i) {
        const auto& traj = trajectories[i];

        if (traj.empty()) continue;

        double sum_density = 0.0;
        for (size_t v : traj) {
            if (v < density.size()) {
                sum_density += density[v];
            }
        }

        double mean_density = sum_density / static_cast<double>(traj.size());

        if (mean_density > best_mean_density) {
            best_mean_density = mean_density;
            best_idx = i;
        }
    }

    return best_idx;
}
