#include <queue>
#include "set_wgraph.hpp"
#include <R_ext/Print.h>

/**
 * @brief Computes the basin of attraction for a given local extremum.
 *
 * @details This function determines the complete basin of attraction by following all paths
 *          where the extremum property is maintained. Unlike the hop neighborhood version,
 *          this function does not stop at the first violation but continues exploring all
 *          non-violating paths to their full extent.
 *
 *          The function performs a breadth-first search (BFS) starting from the given vertex
 *          and includes all vertices that can be reached via monotone paths (strictly
 *          increasing for minima, strictly decreasing for maxima). The search continues
 *          along each path as long as the extremum condition relative to the origin vertex
 *          is satisfied.
 *
 *          Global extrema are processed identically to other vertices, with their basin
 *          potentially encompassing the entire graph if no violations exist.
 *
 * Basin membership criteria:
 * - For a local minimum at v: vertex u is in the basin if it can be reached from v
 *   via a path where y strictly increases along each edge (y[w] > y[prev] for each step)
 * - For a local maximum at v: vertex u is in the basin if it can be reached from v
 *   via a path where y strictly decreases along each edge (y[w] < y[prev] for each step)
 *
 * The hop_idx field in the result now represents:
 * - The maximum hop distance within the basin
 * - std::numeric_limits<size_t>::max() if the basin is empty (vertex is not a 1-hop extremum)
 * - Otherwise, the actual maximum distance from vertex to any basin member
 *
 * @param vertex Index of the vertex (local extremum) whose basin to compute
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, computes basin for a maximum; otherwise, for a minimum
 * @return hop_nbhd_t Structure containing vertex, max distance, hop distance map, and boundary map
 */
hop_nbhd_t set_wgraph_t::compute_basin_of_attraction(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
    ) const {

    hop_nbhd_t result;
    result.vertex  = vertex;

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

    // BFS to find the complete basin of attraction
    std::vector<int> hop_distance(adjacency_list.size(), -1);
    std::vector<size_t> predecessor(adjacency_list.size(), std::numeric_limits<size_t>::max());
    std::queue<size_t> q;

    // Initialize with the start vertex
    result.hop_dist_map[vertex] = 0;
    hop_distance[vertex] = 0;
    q.push(vertex);

    size_t max_hop_distance = 0;

    // BFS with monotonicity constraints along edges
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();

        // For vertices other than the extremum, verify monotonicity along incoming edge
        if (u != vertex) {
            size_t prev = predecessor[u];
            double delta_y = y[u] - y[prev];

            // For maxima: y must decrease along trajectory (delta < 0)
            // For minima: y must increase along trajectory (delta > 0)
            bool monotone = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!monotone) {
                // Monotonicity violation, do not include this vertex in basin
                continue;
            }

            // Add vertex to basin
            result.hop_dist_map[u] = static_cast<size_t>(hop_distance[u]);
            max_hop_distance = std::max(max_hop_distance,
                                       static_cast<size_t>(hop_distance[u]));
        }

        // Explore neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            // If not yet visited
            if (hop_distance[v] == -1) {
                hop_distance[v] = hop_distance[u] + 1;
                predecessor[v] = u;
                q.push(v);
            }
        }
    }

    result.hop_idx = max_hop_distance;

    // Collect boundary vertices (neighbors of basin vertices that are not in the basin)
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
