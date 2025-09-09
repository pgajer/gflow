
/**
 * @brief Finds all shortest paths within a specified radius from a start vertex
 *
 * This function implements a modified version of Dijkstra's algorithm to find all shortest
 * paths from a start vertex to other vertices in the graph, limited by a maximum distance (radius).
 * The algorithm operates in two phases:
 * 1. A bounded Dijkstra's algorithm to find distances and predecessors
 * 2. Path reconstruction to generate actual paths from the collected information
 *
 * The function caches results for subsequent calls with the same start vertex.
 *
 * @param start The index of the starting vertex
 * @param radius The maximum allowed distance for paths from the start vertex
 *
 * @return shortest_paths_t A structure containing:
 *         - paths: Vector of path_t objects, each containing:
 *           * vertices: Sequence of vertices in the path
 *           * total_weight: Total distance of the path
 *         - reachable_vertices: Set of all vertices reachable within the radius
 *
 * @note Paths are returned in descending order of their total weights
 * @note The function caches results in paths_cache for efficiency in subsequent calls
 *
 * @pre The start vertex index must be valid (0 <= start < adjacency_list.size())
 * @pre The radius must be non-negative
 *
 * @complexity Time: O((V + E) * log V) where V is the number of vertices and E is the number of edges
 *            Space: O(V + E) for storing distances, predecessors, and the result paths
 */
shortest_paths_t set_wgraph_t::find_shortest_paths_within_radius(size_t start, double radius) const {
    // Check cache first
    auto cache_it = paths_cache.find(start);
    if (cache_it != paths_cache.end()) {
        return cache_it->second;
    }

    // Phase 1: a bounded version of Dijkstra's algorithm producing
    // std::unordered_map<size_t, std::pair<double, size_t>> map from vertex to its
    // (distance, predecessor) pair, containing only vertices within radius
    //
    // 1. Initializes distances to INFINITY except start vertex (distance 0)
    // 2. Uses priority queue to process vertices in order of increasing distance
    // 3. For each vertex u:
    //    - If distance to u exceeds radius, terminates that branch
    //    - Updates distances to adjacent vertices if they can be reached
    //      with shorter path through u and total distance <= radius
    // 4. Records (distance, predecessor) pairs for all reached vertices

    std::unordered_map<size_t, std::pair<double, int>> dp_map; // dp_map[vertex] = <distance, prev>
    size_t n = adjacency_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<int> prev(n, -1);
    dist[start] = 0;

    std::priority_queue<std::pair<double, int>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto top = pq.top();
        double d = -top.first;
        int u = top.second;
        pq.pop();

        if (d > radius) break;
        if (d > dist[u]) continue;

        dp_map.emplace(u, std::pair<double, int>(d, prev[u]));

        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double w = edge.weight;

            if (dist[u] + w < dist[v] && dist[u] + w <= radius) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});
            }
        }
    }

    // Phase 2: Uses dp_map to coonstruct a vector of corresponding paths
    shortest_paths_t result;

    struct vertex_info_t {
        size_t vertex;
        double distance;
    };

    // Create and sort vertex info
    std::vector<vertex_info_t> vertex_info;
    vertex_info.reserve(dp_map.size());
    for (const auto& [vertex, info] : dp_map) {
        if (vertex != start) {  // exclude target vertex itself
            vertex_info.emplace_back(vertex, info.first);
        }
    }

    // Sort by distance in descending order
    std::sort(vertex_info.begin(), vertex_info.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    // result.reachable_vertices will be used to keep track of used vertices to avoid duplicate paths
    result.reachable_vertices.insert(start);  // target vertex is always used

    // Process vertices in order of decreasing distance
    for (const auto& info : vertex_info) {
        // Skip if this vertex was already used in another path
        if (result.reachable_vertices.count(info.vertex) > 0) {
            continue;
        }

        path_t new_path;
        new_path.total_weight = info.distance;  // distance from start to target

        // Reconstruct the path from this vertex to target
        size_t curr = info.vertex;
        while (curr != -1) {
            new_path.vertices.push_back(curr);
            result.reachable_vertices.insert(curr);  // mark vertex as used

            auto it = dp_map.find(curr);
            if (it == dp_map.end()) {
                REPORT_ERROR("Path reconstruction failed: vertex not found in path info");
            }

            curr = it->second.second;  // move to predecessor
        }

        // Reverse the path to get correct order (from target to end)
        std::reverse(new_path.vertices.begin(), new_path.vertices.end());

        result.paths.push_back(new_path);
    }

    // Sort paths in the descending order of their total weights
    std::sort(result.paths.begin(), result.paths.end());

    paths_cache[start] = result;

    return result;
}
