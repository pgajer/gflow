#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

#include <vector>
#include <queue>
#include <limits>
#include <memory>
#include <algorithm>
#include <unordered_set>

#include "error_utils.h" // for REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
    SEXP S_shortest_path(SEXP s_graph, SEXP s_edge_lengths, SEXP s_vertices);
}

/**
 * @brief Computes shortest path distances between specified vertices in a weighted graph.
 *
 * This function calculates the shortest path distances between all pairs of vertices
 * specified in the input vector. It uses Dijkstra's algorithm to compute these distances
 * efficiently. The function returns a flattened distance matrix where the element at
 * index (i * n + j) represents the shortest distance from vertices[i] to vertices[j],
 * where n is the number of specified vertices.
 *
 * @param graph A vector of vectors representing the adjacency list of the graph.
 *              graph[i] contains the indices of vertices adjacent to vertex i.
 * @param edge_lengths A vector of vectors representing the lengths of edges in the graph.
 *                     edge_lengths[i][j] is the length of the edge from vertex i to its j-th neighbor.
 * @param vertices A vector of vertex indices for which to compute the shortest path distances.
 *
 * @return A unique pointer to a vector of doubles representing the flattened distance matrix.
 *         The matrix is stored in row-major order, where the element at index (i * n + j)
 *         is the shortest distance from vertices[i] to vertices[j]. If there is no path
 *         between two vertices, the corresponding distance will be set to infinity
 *         (std::numeric_limits<double>::infinity()).
 *
 * @throws std::invalid_argument If the graph and edge_lengths structures are inconsistent,
 *                               or if any vertex in 'vertices' is out of range.
 *
 * @note The time complexity of this function is O(k * (V + E) * log(V)), where k is the
 *       number of vertices in the 'vertices' vector, V is the total number of vertices in
 *       the graph, and E is the total number of edges. The space complexity is O(V + k^2).
 *
 * @warning This function assumes that the graph is well-formed and that all vertex
 *          indices are valid. It does not perform extensive input validation.
 *
 * Example usage:
 * @code
 * std::vector<std::vector<int>> graph = {{1, 2}, {0, 2, 3}, {0, 1, 3, 4}, {1, 2, 4}, {2, 3}};
 * std::vector<std::vector<double>> edge_lengths = {{1.0, 4.0}, {1.0, 2.0, 5.0}, {4.0, 2.0, 1.0, 3.0}, {5.0, 1.0, 2.0}, {3.0, 2.0}};
 * std::vector<int> vertices = {0, 2, 4};
 *
 * auto distance_matrix_ptr = shortest_path(graph, edge_lengths, vertices);
 *
 * // To access the distance from vertices[1] to vertices[2]:
 * double distance = (*distance_matrix_ptr)[1 * vertices.size() + 2];
 * @endcode
 */
std::unique_ptr<std::vector<double>> shortest_path(const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& edge_lengths,
                                                   std::vector<int> vertices) {
    const int n = vertices.size();
    auto distance_matrix = std::make_unique<std::vector<double>>(n * n, std::numeric_limits<double>::infinity());

    // Lambda function to get index in the flattened matrix
    auto get_index = [n](int i, int j) { return i * n + j; };

    // Custom comparator for priority queue
    auto cmp = [](std::pair<int, double> left, std::pair<int, double> right) {
        return left.second > right.second;
    };

    for (int i = 0; i < n; ++i) {
        int start = vertices[i];
        std::vector<double> dist(graph.size(), std::numeric_limits<double>::infinity());
        std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, decltype(cmp)> pq(cmp);

        dist[start] = 0;
        pq.push({start, 0});

        while (!pq.empty()) {
            int u = pq.top().first;
            double d = pq.top().second;
            pq.pop();

            if (d > dist[u]) continue;

            // Update distance in the matrix if the vertex is in our target set
            auto it = std::find(vertices.begin(), vertices.end(), u);
            if (it != vertices.end()) {
                int j = std::distance(vertices.begin(), it);
                (*distance_matrix)[get_index(i, j)] = d;
            }

            for (size_t v = 0; v < graph[u].size(); ++v) {
                int to = graph[u][v];
                double len = edge_lengths[u][v];

                if (dist[u] + len < dist[to]) {
                    dist[to] = dist[u] + len;
                    pq.push({to, dist[to]});
                }
            }
        }
    }

    return distance_matrix;
}

/**
 * @brief R interface for the C++ shortest_path function.
 *
 * This function serves as an interface between R and the C++ implementation of the shortest_path
 * algorithm. It takes R objects as input, converts them to appropriate C++ data structures,
 * calls the C++ shortest_path function, and then converts the result back to an R object.
 *
 * @param s_graph SEXP object representing the graph structure. Expected to be a list of integer
 *                vectors, where each vector represents the neighbors of a vertex.
 * @param s_edge_lengths SEXP object representing the edge lengths. Expected to be a list of
 *                       double vectors, corresponding to the lengths of edges in s_graph.
 * @param s_vertices SEXP object representing the vertices to compute distances between.
 *                   Can be either an integer or numeric vector in R.
 *
 * @return SEXP A PROTECT'd R matrix (REALSXP) containing the pairwise shortest distances
 *              between the specified vertices. The matrix is square with dimensions
 *              length(s_vertices) x length(s_vertices). The caller (R) is responsible for
 *              garbage collecting this object.
 *
 * @throws Rf_error If s_vertices is neither an integer nor a numeric vector.
 *
 * @warning Ensure that s_graph and s_edge_lengths have consistent structures. The function
 *          assumes that the input is well-formed and does not perform extensive input validation.
 *
 * @see shortest_path The underlying C++ function that performs the actual computation.
 * @see flat_vector_to_R_matrix Used internally to convert the C++ result to an R matrix.
 *
 * Example R usage:
 * @code
 * # In R:
 * result <- .Call("S_shortest_path", graph, edge_lengths, vertices)
 * @endcode
 *
 * Where:
 * - graph is a list of integer vectors
 * - edge_lengths is a list of numeric vectors
 * - vertices is an integer or numeric vector
 */
SEXP S_shortest_path(SEXP s_graph, SEXP s_edge_lengths, SEXP s_vertices) {

    // Convert R objects to C++ objects
    std::vector<std::vector<int>> graph           = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(s_edge_lengths);

    // Convert s_vertices to std::vector<int>
    int vertices_length = Rf_length(s_vertices);
    std::vector<int> vertices(vertices_length);

    // Check if s_vertices is integer or real
    switch(TYPEOF(s_vertices)) {
        case INTSXP: {
            int* vertices_ptr = INTEGER(s_vertices);
            for (int i = 0; i < vertices_length; ++i) {
                vertices[i] = vertices_ptr[i];
            }
            break;
        }
        case REALSXP: {
            double* vertices_ptr = REAL(s_vertices);
            for (int i = 0; i < vertices_length; ++i) {
                vertices[i] = static_cast<int>(vertices_ptr[i]);
            }
            break;
        }
        default:
            error("s_vertices must be a vector of integers or reals");
    }

    // Call the C++ shortest_path function
    std::unique_ptr<std::vector<double>> result = shortest_path(graph, edge_lengths, vertices);

    // Convert the result to an R matrix
    SEXP r_matrix = PROTECT(flat_vector_to_R_matrix(*result, vertices_length, vertices_length));

    UNPROTECT(1);
    return r_matrix;
}




/**
 * @brief Computes the shortest path between two vertices in a weighted graph.
 *
 * This function calculates the shortest path between two specified vertices
 * using Dijkstra's algorithm and returns the sequence of vertices along that path.
 *
 * @param graph A vector of vectors representing the adjacency list of the graph.
 *              graph[i] contains the indices of vertices adjacent to vertex i.
 * @param edge_lengths A vector of vectors representing the lengths of edges in the graph.
 *                     edge_lengths[i][j] is the length of the edge from vertex i to its j-th neighbor.
 * @param start The index of the starting vertex.
 * @param end The index of the ending vertex.
 *
 * @return A vector of integers representing the sequence of vertices along the shortest path
 *         from 'start' to 'end', including both endpoints. If there is no path between the two
 *         vertices, the function returns an empty vector.
 *
 * @throws std::invalid_argument If the graph and edge_lengths structures are inconsistent,
 *                               or if 'start' or 'end' is out of range.
 *
 * @note The time complexity of this function is O((V + E) * log(V)), where V is the total
 *       number of vertices in the graph and E is the total number of edges.
 *
 * @warning This function assumes that the graph is well-formed and that all vertex
 *          indices are valid. It does not perform extensive input validation.
 *
 * Example usage:
 * @code
 * std::vector<std::vector<int>> graph = {{1, 2}, {0, 2, 3}, {0, 1, 3, 4}, {1, 2, 4}, {2, 3}};
 * std::vector<std::vector<double>> edge_lengths = {
 *     {1.0, 4.0},           // Edges from vertex 0
 *     {1.0, 2.0, 5.0},      // Edges from vertex 1
 *     {4.0, 2.0, 1.0, 3.0}, // Edges from vertex 2
 *     {5.0, 1.0, 2.0},      // Edges from vertex 3
 *     {3.0, 2.0}            // Edges from vertex 4
 * };
 * int start = 0;
 * int end = 4;
 *
 * std::vector<int> path = shortest_path_bw_two_vertices(graph, edge_lengths, start, end);
 *
 * @endcode
 */
std::vector<int> shortest_path_bw_two_vertices(const std::vector<std::vector<int>>& graph,
                                               const std::vector<std::vector<double>>& edge_lengths,
                                               int start, int end) {
    const int V = graph.size();

    if (start < 0 || start >= V || end < 0 || end >= V) {
        Rf_error("Start or end vertex is out of range.");
    }

    // Distance vector and predecessor vector
    std::vector<double> dist(V, std::numeric_limits<double>::infinity());
    std::vector<int> prev(V, -1);

    // Custom comparator for priority queue
    auto cmp = [](const std::pair<int, double>& left, const std::pair<int, double>& right) {
        return left.second > right.second;
    };

    // Min-priority queue
    std::priority_queue<std::pair<int, double>,
                        std::vector<std::pair<int, double>>,
                        decltype(cmp)> pq(cmp);

    dist[start] = 0.0;
    pq.push({start, 0.0});

    while (!pq.empty()) {
        int u = pq.top().first;
        double d = pq.top().second;
        pq.pop();

        if (u == end) {
            break;  // Found the shortest path to the destination
        }

        if (d > dist[u]) {
            continue;  // Skip outdated entries
        }

        for (size_t i = 0; i < graph[u].size(); ++i) {
            int v = graph[u][i];
            double len = edge_lengths[u][i];

            if (dist[u] + len < dist[v]) {
                dist[v] = dist[u] + len;
                prev[v] = u;
                pq.push({v, dist[v]});
            }
        }
    }

    // Reconstruct the path from start to end
    std::vector<int> path;
    if (dist[end] == std::numeric_limits<double>::infinity()) {
        // No path exists
        return path;
    }

    for (int at = end; at != -1; at = prev[at]) {
        path.push_back(at);
    }

    std::reverse(path.begin(), path.end());

    return path;
}


/**
 * @brief Finds the path with minimum number of edges between two vertices in an unweighted graph
 *
 * @param graph Adjacency list representation of the graph
 * @param start Starting vertex index
 * @param end Target vertex index
 * @return std::vector<int> Vector of vertex indices representing the shortest path from start to end.
 *         Empty vector if no path exists.
 * @throws std::invalid_argument if start or end vertex indices are out of range
 */
std::vector<int> shortest_path_by_hops(const std::vector<std::vector<int>>& graph,
                                       int start, int end) {
    const int n_vertices = graph.size();
    if (start < 0 || start >= n_vertices || end < 0 || end >= n_vertices) {
        Rf_error("Start or end vertex is out of range.");
    }

    // Check if vertices are the same
    // if (start == end) {
    //     return {start};
    // }

    // Check if vertices are neighbors
    if (std::find(graph[start].begin(), graph[start].end(), end) != graph[start].end()) {
        return {start, end};
    }

    // If not neighbors, proceed with BFS
    std::vector<int> dist(n_vertices, std::numeric_limits<int>::max());
    std::vector<int> prev(n_vertices, -1);
    std::queue<int> q;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == end) {
            break;  // Found the shortest path to the destination
        }

        for (int v : graph[u]) {
            if (dist[v] == std::numeric_limits<int>::max()) {
                dist[v] = dist[u] + 1;
                prev[v] = u;
                q.push(v);
            }
        }
    }

    // Reconstruct the path from start to end
    std::vector<int> path;
    if (dist[end] == std::numeric_limits<int>::max()) {
        // No path exists
        return path;
    }

    for (int at = end; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());
    return path;
}


/**
 * @brief Find shortest paths from source vertex to all vertices using Dijkstra's algorithm
 *
 * @param adj_list Adjacency list representation of the graph
 * @param weight_list Weight list corresponding to adjacency list edges
 * @param start Source vertex from which to compute shortest paths
 * @return Pair containing:
 *         - Vector of shortest distances to each vertex
 *         - Vector of previous vertices in shortest paths
 */
std::pair<std::vector<double>, std::vector<int>> find_all_shortest_paths_from_vertex(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start
) {
    int n = adj_list.size();
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

/**
 * @brief Reconstructs a path from start to end vertex using a previous vertex array
 *
 * @param start The starting vertex ID
 * @param end The ending vertex ID
 * @param prev Vector where prev[i] contains the previous vertex in the shortest path to vertex i.
 *             A value of -1 indicates no previous vertex (i.e., the start of the path)
 *
 * @return std::vector<int> The reconstructed path from start to end vertex (inclusive)
 *
 * @throws std::runtime_error if the reconstructed path doesn't begin at the specified start vertex
 *
 * @details This function reconstructs a path through a graph using a "previous vertex" array
 *          typically generated by shortest path algorithms like Dijkstra's or BFS.
 *          It works backwards from the end vertex, following the prev array until reaching
 *          the start of the path (indicated by -1), then reverses the result to get
 *          the path in the correct order.
 *
 * @note The function assumes that prev[i] contains valid vertex IDs or -1
 *       The resulting path includes both the start and end vertices
 */
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

    if (path.front() != start) {
        error("Path reconstruction failed: path does not start at the specified vertex");
    }

    return path;
}

/**
 * @brief Computes shortest path distances from a start vertex to specified target vertices using Dijkstra's algorithm
 *
 * @param adj_list Adjacency list representing the graph where adj_list[u] contains vertices adjacent to u
 * @param weight_list Weight list where weight_list[u][i] is the weight of edge from u to adj_list[u][i]
 * @param start_vertex Starting vertex for path computation
 * @param target_vertices Vector of vertices to compute shortest paths to
 * @return std::vector<double> Vector containing shortest path distances to target vertices (INFINITY if no path exists)
 *
 * @throws std::invalid_argument if adjacency and weight lists have inconsistent sizes
 * @throws std::out_of_range if start_vertex or any target vertex is out of range
 *
 * @note Time complexity: O((V' + E') log V') where V' and E' are the number of vertices and edges
 *       explored before finding all targets (can be much less than total V and E)
 * @note Space complexity: O(V)
 */
std::vector<double> compute_shortest_path_distances(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    size_t start_vertex,
    const std::vector<size_t>& target_vertices
) {
    // Validate input
    if (adj_list.size() != weight_list.size()) {
        REPORT_ERROR("Adjacency and weight lists must have the same size");
    }

    size_t n_vertices = adj_list.size();
    if (start_vertex >= n_vertices) {
        REPORT_ERROR("Start vertex out of range");
    }

    for (size_t v : target_vertices) {
        if (v >= n_vertices) {
            REPORT_ERROR("Target vertex out of range");
        }
    }

    // Create set of remaining target vertices to find
    std::unordered_set<size_t> remaining_targets(target_vertices.begin(), target_vertices.end());
    if (remaining_targets.empty()) {
        return {};
    }

    // If start vertex is one of the targets, remove it immediately
    if (remaining_targets.count(start_vertex)) {
        remaining_targets.erase(start_vertex);
    }

    // Initialize distances only for vertices we've seen
    std::unordered_map<size_t, double> dist;
    dist[start_vertex] = 0;

    // Min-heap priority queue: pairs of (distance, vertex)
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, start_vertex});

    while (!pq.empty() && !remaining_targets.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        d = -d;  // Convert back to positive distance

        // Skip if we already found a better path
        if (dist.count(u) && d > dist[u]) continue;

        // Process all neighbors
        for (size_t i = 0; i < adj_list[u].size(); i++) {
            size_t v = static_cast<size_t>(adj_list[u][i]);
            double w = weight_list[u][i];

            double new_dist = dist[u] + w;
            if (!dist.count(v) || new_dist < dist[v]) {
                dist[v] = new_dist;
                pq.push({-new_dist, v});

                // If this is a target vertex, remove it from remaining targets
                if (remaining_targets.count(v)) {
                    remaining_targets.erase(v);
                    if (remaining_targets.empty()) break;  // Exit early if all targets found
                }
            }
        }
    }

    // Extract distances for target vertices
    std::vector<double> result;
    result.reserve(target_vertices.size());
    for (size_t v : target_vertices) {
        // If we haven't reached a target, it's not reachable
        result.push_back(dist.count(v) ? dist[v] : std::numeric_limits<double>::infinity());
    }

    return result;
}
