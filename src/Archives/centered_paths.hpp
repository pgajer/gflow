#ifndef CENTERED_PATHS_H_
#define CENTERED_PATHS_H_

#include <vector>
#include <queue>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>  // for std::reverse

/**
 * @brief Structure representing a path in a graph
 *
 * Stores information about a path including its vertices, total weight (length),
 * and how far a specified vertex is from the center of the path.
 */
struct path_t {
    std::vector<int> vertices;    ///< Sequence of vertices in the path
    double total_weight;          ///< Total length of the path
    double center_offset;         ///< Distance of target vertex from path center (0 = perfectly centered)

    path_t() : total_weight(0), center_offset(INFINITY) {}

    // For priority queue ordering
    bool operator<(const path_t& other) const {
        if (std::abs(total_weight - other.total_weight) > 1e-10) {
            return total_weight < other.total_weight;  // prefer paths closer to target length
        }
        return center_offset > other.center_offset;    // prefer more centered paths
    }
};

/**
 * @brief Class for finding centered paths in a weighted graph
 *
 * Finds paths of specified length that contain a given vertex,
 * attempting to position that vertex as close to the center as possible.
 */
struct path_finder_t {

    const std::vector<std::vector<int>>& adj_list;
    const std::vector<std::vector<double>>& weight_list;
    const double target_length;  // 2 * bandwidth
    const double length_tolerance;  // acceptable deviation from target_length

    /**
     * @brief Constructor for path finder
     *
     * @param adj Adjacency list representation of the graph where adj[i] contains indices of vertices adjacent to vertex i
     * @param weights Weight list where weights[i] contains weights of edges from vertex i to vertices in adj[i]
     * @param bandwidth Half of the desired path length
     * @param tolerance Acceptable deviation from target path length as a fraction (default: 0.1)
     */
    path_finder_t(
        const std::vector<std::vector<int>>& adj,
        const std::vector<std::vector<double>>& weights,
        double bandwidth,
        double tolerance = 0.1
    ) : adj_list(adj),
        weight_list(weights),
        target_length(2.0 * bandwidth),
        length_tolerance(tolerance) {}

    /**
     * @brief Find multiple centered paths containing a specified vertex
     *
     * Searches for paths that:
     * 1. Have length approximately 2 * bandwidth
     * 2. Contain the specified vertex
     * 3. Have the specified vertex as close to center as possible
     * 4. Share no vertices except the specified vertex
     *
     * @param v Target vertex to be centered
     * @param max_paths Maximum number of paths to return
     * @return Vector of path_t objects, ordered by quality (considering length match and centeredness)
     */
    std::vector<path_t> find_centered_paths(
        int v,
        int max_paths = 5
    ) {
        std::vector<path_t> result_paths;
        std::set<int> used_vertices = {v};  // track vertices used in accepted paths
        std::vector<double> radius_multipliers = {1.2, 1.5, 1.8, 2.0};

        for (double radius_mult : radius_multipliers) {
            auto candidates = find_paths_with_radius(v, target_length/2 * radius_mult);

            // Sort candidates by quality (length match and centeredness)
            std::priority_queue<path_t> pq;
            for (const auto& path : candidates) {
                pq.push(path);
            }

            // Try to add disjoint paths
            while (!pq.empty() && result_paths.size() < max_paths) {
                path_t current = pq.top();
                pq.pop();

                bool is_disjoint = true;
                for (int vertex : current.vertices) {
                    if (vertex != v && used_vertices.count(vertex) > 0) {
                        is_disjoint = false;
                        break;
                    }
                }

                if (is_disjoint) {
                    // Add all vertices except v to used_vertices
                    for (int vertex : current.vertices) {
                        if (vertex != v) {
                            used_vertices.insert(vertex);
                        }
                    }
                    result_paths.push_back(current);
                }
            }

            if (!result_paths.empty()) break;  // Stop if we found any valid paths
        }

        return result_paths;
    }

private:
    std::vector<path_t> find_paths_with_radius(int v, double search_radius) {
        std::vector<path_t> candidates;
        int n = adj_list.size();

        // Get distances from v
        auto [dist_from_v, prev_from_v] = run_dijkstra(v);

        // Find endpoints within search_radius
        std::vector<int> endpoints;
        for (int i = 0; i < n; i++) {
            if (dist_from_v[i] <= search_radius) {
                endpoints.push_back(i);
            }
        }

        // Try all endpoint pairs
        for (int start : endpoints) {
            auto [dist_from_start, prev_from_start] = run_dijkstra(start);

            for (int end : endpoints) {
                if (start >= end) continue;

                // Check if path length is close to target
                double path_length = dist_from_start[end];
                if (std::abs(path_length - target_length) >
                    target_length * length_tolerance) continue;

                // Reconstruct path and verify v is on it
                auto path = reconstruct_path(start, end, prev_from_start);
                if (!is_vertex_on_path(v, path)) continue;

                // Calculate path properties
                path_t candidate;
                candidate.vertices = path;
                candidate.total_weight = path_length;
                candidate.center_offset = calculate_center_offset(v, path);

                candidates.push_back(candidate);
            }
        }

        return candidates;
    }

    std::pair<std::vector<double>, std::vector<int>> run_dijkstra(int start) {
        int n = adj_list.size();
        std::vector<double> dist(n, INFINITY);
        std::vector<int> prev(n, -1);
        dist[start] = 0;

        std::priority_queue<std::pair<double, int>> pq;
        pq.push({0, start});

        while (!pq.empty()) {
            double d = -pq.top().first;
            int u = pq.top().second;
            pq.pop();

            if (d > dist[u]) continue;

            for (size_t i = 0; i < adj_list[u].size(); i++) {
                int v = adj_list[u][i];
                double w = weight_list[u][i];

                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    prev[v] = u;
                    pq.push({-dist[v], v});
                }
            }
        }

        return {dist, prev};
    }

    std::vector<int> reconstruct_path(
        int start,
        int end,
        const std::vector<int>& prev
        ) {
        std::vector<int> path;
        for (int curr = end; curr != -1; curr = prev[curr]) {
            path.push_back(curr);
        }
        std::reverse(path.begin(), path.end());

        // Validate that path starts at start vertex
        if (path.front() != start) {
            error("Path reconstruction failed: path does not start at the specified vertex");
        }

        return path;
    }



    bool is_vertex_on_path(int v, const std::vector<int>& path) {
        return std::find(path.begin(), path.end(), v) != path.end();
    }

    double calculate_center_offset(
        int v,
        const std::vector<int>& path
    ) {
        double path_position = 0;
        for (size_t i = 0; i < path.size(); i++) {
            if (path[i] == v) {
                path_position = i;
                break;
            }
        }
        return std::abs(path_position / (path.size() - 1) - 0.5);
    }
};

#endif // CENTERED_PATHS_H_
