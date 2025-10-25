#include <queue>
#include "set_wgraph.hpp"
#include <R_ext/Print.h>

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
gradient_basin_t set_wgraph_t::compute_basin_of_attraction(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
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

    // BFS to find the complete basin of attraction
    std::vector<int> hop_distance(adjacency_list.size(), -1);
    std::vector<size_t> predecessor(adjacency_list.size(), std::numeric_limits<size_t>::max());
    std::queue<size_t> q;

    // Initialize with the start vertex
    result.hop_dist_map[vertex] = 0;
    hop_distance[vertex] = 0;
    predecessor[vertex] = std::numeric_limits<size_t>::max();
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
            result.predecessors[u] = prev;
            max_hop_distance = std::max(max_hop_distance,
                                       static_cast<size_t>(hop_distance[u]));

            // Check if u is a terminal extremum (opposite type from origin)
            bool is_terminal_extremum = true;
            for (const auto& edge : adjacency_list[u]) {
                size_t w = edge.vertex;
                // For descending basin (from max), terminals are minima: neighbors satisfy y[w] >= y[u]
                // For ascending basin (from min), terminals are maxima: neighbors satisfy y[w] <= y[u]
                if (detect_maxima) {
                    if (y[w] < y[u]) {
                        is_terminal_extremum = false;
                        break;
                    }
                } else {
                    if (y[w] > y[u]) {
                        is_terminal_extremum = false;
                        break;
                    }
                }
            }

            if (is_terminal_extremum) {
                result.terminal_extrema.push_back(u);
            }
        }

        // Explore neighbors for BFS expansion
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
