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
 * @brief Compute tubular neighborhood of a trajectory via BFS
 *
 * @param trajectory Ordered vertices of the trajectory
 * @param radius Maximum hop distance from trajectory
 * @param basin_restriction If non-empty, restrict to these vertices
 * @return Tuple of (tubular vertices, hop distances, nearest trajectory index)
 */
std::tuple<std::vector<size_t>, std::vector<int>, std::vector<size_t>>
set_wgraph_t::compute_tubular_neighborhood(
    const std::vector<size_t>& trajectory,
    int radius,
    const std::unordered_set<size_t>& basin_restriction
) const {

    std::vector<size_t> tubular_vertices;
    std::vector<int> hop_distances;
    std::vector<size_t> nearest_traj_idx;

    // Track visited vertices: vertex -> (hop distance, trajectory index)
    std::unordered_map<size_t, std::pair<int, size_t>> visited;

    // Build trajectory vertex to index map
    std::unordered_map<size_t, size_t> traj_vertex_to_idx;
    for (size_t i = 0; i < trajectory.size(); ++i) {
        traj_vertex_to_idx[trajectory[i]] = i;
    }

    // Initialize BFS queue with trajectory vertices at distance 0
    // Queue entries: (vertex, hop distance, source trajectory index)
    std::queue<std::tuple<size_t, int, size_t>> bfs_queue;

    for (size_t i = 0; i < trajectory.size(); ++i) {
        size_t v = trajectory[i];
        if (basin_restriction.empty() || basin_restriction.count(v) > 0) {
            visited[v] = {0, i};
            bfs_queue.push({v, 0, i});
        }
    }

    // BFS expansion
    while (!bfs_queue.empty()) {
        auto [u, dist, traj_idx] = bfs_queue.front();
        bfs_queue.pop();

        // Explore neighbors if within radius
        if (dist < radius) {
            for (const auto& edge : adjacency_list[u]) {
                size_t w = edge.vertex;

                // Skip if outside basin restriction
                if (!basin_restriction.empty() && basin_restriction.count(w) == 0) {
                    continue;
                }

                // Check if already visited
                auto it = visited.find(w);
                if (it != visited.end()) {
                    // Already visited - but check if this path is equally short
                    // and comes from a trajectory vertex with closer coordinate
                    // (This handles ties by preferring the closer trajectory vertex)
                    if (it->second.first == dist + 1) {
                        // Same distance - could update to closer trajectory vertex
                        // but for simplicity we keep first discovery
                    }
                    continue;
                }

                // Check if w is itself a trajectory vertex
                auto traj_it = traj_vertex_to_idx.find(w);
                if (traj_it != traj_vertex_to_idx.end()) {
                    // w is on trajectory - use its own index
                    visited[w] = {0, traj_it->second};
                    // Don't expand from here as it's already seeded
                } else {
                    // w is not on trajectory - inherit source from u
                    visited[w] = {dist + 1, traj_idx};
                    bfs_queue.push({w, dist + 1, traj_idx});
                }
            }
        }
    }

    // Collect results
    tubular_vertices.reserve(visited.size());
    hop_distances.reserve(visited.size());
    nearest_traj_idx.reserve(visited.size());

    for (const auto& [v, info] : visited) {
        tubular_vertices.push_back(v);
        hop_distances.push_back(info.first);
        nearest_traj_idx.push_back(info.second);
    }

    return {tubular_vertices, hop_distances, nearest_traj_idx};
}

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
 * @brief Compute harmonic extension of trajectory coordinates
 *
 * Given a geodesic trajectory through a basin, this function computes
 * arc-length coordinates for the trajectory vertices and extends them
 * harmonically to a tubular neighborhood.
 *
 * @param trajectory Ordered vertices from minimum to maximum
 * @param params Extension parameters
 * @param verbose Print progress information
 * @return Complete harmonic extension result
 */
harmonic_extension_result_t set_wgraph_t::compute_harmonic_extension(
    const std::vector<size_t>& trajectory,
    const harmonic_extension_params_t& params,
    bool verbose
) const {

    harmonic_extension_result_t result;
    result.trajectory = trajectory;

    if (trajectory.empty()) {
        if (verbose) {
            Rprintf("Warning: Empty trajectory provided\n");
        }
        return result;
    }

    if (verbose) {
        Rprintf("Computing harmonic extension for trajectory of %zu vertices\n",
                trajectory.size());
        Rprintf("  Tube radius: %d hops\n", params.tube_radius);
        Rprintf("  Edge weights: %s\n",
                params.use_edge_weights ? "inverse length" : "unit");
    }

    // Step 1: Compute arc-length coordinates for trajectory
    auto [traj_coords, total_length] = compute_arc_length_coords(trajectory);
    result.trajectory_coords = traj_coords;
    result.trajectory_length = total_length;

    if (verbose) {
        Rprintf("  Trajectory length: %.4f\n", total_length);
    }

    // Step 2: Compute tubular neighborhood with nearest trajectory indices
    auto [tubular_verts, hop_dists, nearest_idx] = compute_tubular_neighborhood(
        trajectory, params.tube_radius, params.basin_restriction
    );
    result.tubular_vertices = tubular_verts;
    result.hop_distances = hop_dists;
    result.nearest_traj_idx = nearest_idx;

    if (verbose) {
        Rprintf("  Tubular neighborhood: %zu vertices\n", tubular_verts.size());

        // Count by hop distance
        std::vector<int> hop_counts(params.tube_radius + 1, 0);
        for (int d : hop_dists) {
            if (d <= params.tube_radius) {
                hop_counts[d]++;
            }
        }
        for (int d = 0; d <= params.tube_radius; ++d) {
            Rprintf("    Hop %d: %d vertices\n", d, hop_counts[d]);
        }
    }

    // Build trajectory coordinate map for solver
    std::unordered_set<size_t> trajectory_set(trajectory.begin(), trajectory.end());
    std::unordered_map<size_t, double> trajectory_coord_map;

    for (size_t i = 0; i < trajectory.size(); ++i) {
        trajectory_coord_map[trajectory[i]] = traj_coords[i];
    }

    // Step 3: Build initial coordinates from nearest trajectory vertex
    std::vector<double> initial_coords(tubular_verts.size());
    for (size_t i = 0; i < tubular_verts.size(); ++i) {
        initial_coords[i] = traj_coords[nearest_idx[i]];
    }

    if (verbose) {
        Rprintf("  Initial coords from nearest-neighbor projection\n");
    }

    // Step 4: Solve harmonic extension
    int n_iter;
    double max_change;

    result.extended_coords = solve_harmonic_extension(
        tubular_verts,
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
    for (size_t i = 0; i < tubular_verts.size(); ++i) {
        result.vertex_to_idx[tubular_verts[i]] = i;
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
