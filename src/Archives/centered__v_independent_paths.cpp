#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>  // for std::reverse
#include <sstream>
#include <string>

#include "centered_paths.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_shortest_path.hpp"
#include "cpp_utils.hpp"

extern "C" {
    SEXP S_find_centered_v_independent_paths(
    SEXP adj_list_s,
    SEXP weight_list_s,
    SEXP vertex_s,
    SEXP target_length_s,
    SEXP length_tolerance_s,
    SEXP relative_center_offset_tolerance_s);
}

#if 0
/**
 * Prints a vector of bounded path information to stdout using Rprintf.
 * This function is designed to display path information in a clear, tabular format
 * that makes it easy to trace the path from any vertex back to its start point.
 *
 * The output is formatted as a table with columns for:
 * - Vertex ID: The current vertex in the path
 * - Distance: The accumulated distance from the start to this vertex
 * - Previous: The previous vertex in the path (useful for path reconstruction)
 *
 * @param paths Vector of bounded_path_info_t structures containing path information
 * @param name Optional name to be printed before the path information
 */
void print_bounded_paths(const std::vector<bounded_path_info_t>& paths,
                         const std::string& name = "") {

    // First, we'll check if there's anything to print
    if (paths.empty()) {
        Rprintf("No path information available.\n");
        return;
    }

    // If a name was provided, print it on a separate line
    if (!name.empty()) {
        Rprintf("%s:\n", name.c_str());
    }

    // Print the number of vertices in the path information
    Rprintf("Path information for %d vertices:\n\n", static_cast<int>(paths.size()));

    // Create a header for our table with appropriate spacing
    // We use consistent column widths to ensure alignment
    Rprintf("%-10s %-15s %-10s\n", "Vertex", "Distance", "Previous");

    // Add a separator line to visually distinguish the header
    // The length matches the total width of our columns
    Rprintf("----------------------------------------\n");

    // Iterate through each path entry and print its information
    for (const auto& path_info : paths) {
        // For vertices with no previous vertex (like the start vertex),
        // we'll display "-" instead of -1 for better readability
        std::string prev_vertex = (path_info.prev == -1) ?
                                 "-" :
                                 std::to_string(path_info.prev);

        // Print each row with consistent formatting:
        // - Vertex ID is left-aligned in a 10-character field
        // - Distance is displayed with 4 decimal places in a 15-character field
        // - Previous vertex is left-aligned in a 10-character field
        Rprintf("%-10d %-15.4f %-10s\n",
                path_info.vertex,
                path_info.distance,
                prev_vertex.c_str());
    }

    // Add a final newline for better separation from subsequent output
    Rprintf("\n");
}

/**
 * Prints the contents of a vector of paths to stdout using Rprintf
 * Each path is formatted to show its components in a clear, readable manner
 *
 * @param paths Vector of path_t structures to be printed
 */
void print_paths(const std::vector<path_t>& paths) {
    // Print header for the output
    Rprintf("Number of paths: %d\n", static_cast<int>(paths.size()));

    // Iterate through each path
    for (size_t i = 0; i < paths.size(); ++i) {
        const path_t& path = paths[i];

        // Convert path vertices to string representation
        std::stringstream vertices_str;
        vertices_str << "[";
        for (size_t j = 0; j < path.vertices.size(); ++j) {
            vertices_str << path.vertices[j];
            if (j < path.vertices.size() - 1) {
                vertices_str << ", ";
            }
        }
        vertices_str << "]";

        // Print path information with proper formatting
        Rprintf("\nPath %d:\n", static_cast<int>(i + 1));
        Rprintf("  Vertices:       %s\n", vertices_str.str().c_str());
        Rprintf("  Target Vertex:  %d\n", path.target_vertex);
        Rprintf("  Total Weight:   %.4f\n", path.total_weight);

        // Handle special case for infinite center offset
        if (std::isinf(path.center_offset)) {
            Rprintf("  Center Offset:  Inf\n");
        } else {
            Rprintf("  Center Offset:  %.4f\n", path.center_offset);
        }

        // Add separator between paths for better readability
        if (i < paths.size() - 1) {
            Rprintf("----------------------------------------\n");
        }
    }
    Rprintf("\n");  // Final newline for clean separation from subsequent output
}
#endif

/**
 * @brief Calculates how centered a vertex is in a path using relative center offset
 *
 * @param v The vertex ID to check for centrality
 * @param path Vector of vertex IDs representing the path
 *
 * @return double The center offset value:
 *         - 0.0 means the vertex is exactly in the middle of the path
 *         - 0.5 means the vertex is at either end of the path
 *         - Values between 0.0 and 0.5 indicate how far from center the vertex is
 *         - INFINITY if the vertex is not in the path
 *
 * @details Calculates the normalized position of vertex v in the path, where:
 *          - The first vertex in the path has position 0.0
 *          - The last vertex in the path has position 1.0
 *          Then computes the absolute difference from 0.5 to determine
 *          how far the vertex is from the center position.
 */
static double calculate_relative_center_offset(int v, const std::vector<int>& path) {
    auto it = std::find(path.begin(), path.end(), v);
    if (it == path.end()) return INFINITY;

    double path_position = std::distance(path.begin(), it);
    return std::abs(path_position / (path.size() - 1) - 0.5);
}

/**
 * @brief Reconstructs a path from start to end vertex using bounded path information
 *
 * @param start The starting vertex ID
 * @param end The ending vertex ID
 * @param path_info Vector of bounded_path_info_t structures containing path information
 *                  for vertices reachable within a certain radius
 *
 * @return std::vector<int> The reconstructed path from start to end vertex (inclusive)
 *
 * @throws std::runtime_error if:
 *         - Any vertex in the path is not found in path_info
 *         - The reconstructed path doesn't begin at the specified start vertex
 *
 * @details This function reconstructs a path through a graph using bounded_path_info_t
 *          structures that contain the previous vertex information. It works backwards
 *          from the end vertex, following the prev pointers until reaching the start
 *          vertex, then reverses the result to get the path in the correct order.
 *
 * @note The function expects that if a vertex is referenced in prev, its information
 *       must exist in path_info
 * @note The resulting path includes both the start and end vertices
 *
 * @see bounded_path_info_t for the structure containing path information
 */
static std::vector<int> reconstruct_bounded_path(
    int start,
    int end,
    const std::unordered_map<int, std::pair<double, int>>& vertex_map
) {
    std::vector<int> path;
    int curr = end;
    while (curr != -1) {
        auto it = vertex_map.find(curr);
        if (it == vertex_map.end()) {
            error("Path reconstruction failed: vertex not found in path info");
        }
        path.push_back(curr);
        curr = it->second.second;  // because it->second is a std::pair<double, int>, and we want the second element (prev vertex) of this pair.
    }

    std::reverse(path.begin(), path.end());
    if (path.front() != start) {
        error("Path reconstruction failed: path does not start at the specified vertex");
    }
    return path;
}


/**
 * @brief Performs a bounded Dijkstra search to find shortest paths within a given radius
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param weight_list Weight list where weight_list[i] contains the weights of edges
 *                    from vertex i to its adjacent vertices in adj_list[i]
 * @param start Starting vertex for the search
 * @param radius Maximum allowed distance from the start vertex
 *
 * @return std::unordered_map<int, std::pair<double, int>> Map from vertex to its
 *         (distance, predecessor) pair, containing only vertices within radius
 *
 * @details The function implements a bounded version of Dijkstra's algorithm:
 *          1. Initializes distances to INFINITY except start vertex (distance 0)
 *          2. Uses priority queue to process vertices in order of increasing distance
 *          3. For each vertex u:
 *             - If distance to u exceeds radius, terminates that branch
 *             - Updates distances to adjacent vertices if they can be reached
 *               with shorter path through u and total distance <= radius
 *          4. Records (distance, predecessor) pairs for all reached vertices
 *
 * @note The search terminates early for any path exceeding the radius
 * @note Only vertices reachable within the radius appear in the result map
 *
 * @pre Graph must be weighted and directed. For undirected graphs, each edge should
 *      appear in both directions in the adjacency list
 * @pre Edge weights must be non-negative
 * @pre radius must be non-negative
 * @pre weight_list[i].size() must equal adj_list[i].size() for all i
 *
 * @see reconstruct_bounded_path for path reconstruction using the returned map
 */
static std::unordered_map<int, std::pair<double, int>> find_shortest_paths_within_radius(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start,
    double radius
) {
    std::unordered_map<int, std::pair<double, int>> result; // result[vertex] = <distance, prev>
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

        if (d > radius) break;
        if (d > dist[u]) continue;

        result.emplace(u, std::pair<double, int>(d, prev[u]));

        for (size_t i = 0; i < adj_list[u].size(); i++) {
            int v = adj_list[u][i];
            double w = weight_list[u][i];
            if (dist[u] + w < dist[v] && dist[u] + w <= radius) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});
            }
        }
    }

    return result;
}

//
// Finding v-independent paths from BFS results algorithm
//

// 1. Initialize paths:
//    - Create sorted vertex_info as above
//    - Initialize empty vector<PathInfo> paths

// 2. Generate initial paths:
//    For each VertexInfo in vertex_info (already sorted by distance):
//    - If !available, continue
//    - Generate path using reconstruct_bounded_path
//    - Mark all vertices in path as !available
//    - Add to paths with its vertex_set

// 3. Filter for v-independence:
//    - Start with longest path (paths[0])
//    - For each subsequent path i:
//      - Check intersection of vertex_sets (excluding v)
//      - If intersection empty, keep path
//      - Otherwise, discard
//    - Repeat with next longest remaining path

struct vertex_info_t {
    int vertex;
    double distance;
    bool available;
};

struct path_info_t {
    std::vector<int> vertices;
    std::unordered_set<int> vertex_set;
    double length;

    path_info_t(const std::vector<int>& verts, double len)
        : vertices(verts), length(len) {
        vertex_set.insert(vertices.begin(), vertices.end());
    }
};

/**
 * @brief Finds v-independent paths from bounded Dijkstra search results
 *
 * @param bfs_map Map from vertex to (distance, predecessor) pair from find_shortest_paths_within_radius
 * @param v Target vertex that must be in all paths
 * @param max_length Maximum allowed path length
 *
 * @return Vector of path_t structures representing v-independent paths
 *
 * @details Implements efficient path finding using:
 *          1. Single initial sort of vertices by distance from target
 *          2. O(1) vertex availability tracking
 *          3. O(1) path intersection testing via hash sets
 *          4. Efficient path reconstruction using predecessor information
 *          5. Progressive filtering to ensure v-independence of paths
 *
 * @note Paths are returned in order of decreasing length
 * @note All paths contain vertex v
 * @note Any two paths in the result share only vertex v
 *
 * @see find_shortest_paths_within_radius for generation of the bfs_map
 * @see reconstruct_bounded_path for path reconstruction details
 * @see path_t for path structure definition
 */
static std::vector<path_t> find_v_independent_paths(
    const std::unordered_map<int, std::pair<double, int>>& bfs_map,
    int v,
    double max_length
) {
    // Create and sort vertex info
    std::vector<vertex_info_t> vertex_info;
    vertex_info.reserve(bfs_map.size());
    for (const auto& [vertex, info] : bfs_map) {
        if (vertex != v) {  // exclude v itself
            vertex_info.push_back({vertex, info.first, true});
        }
    }

    // Sort by distance in descending order
    std::sort(vertex_info.begin(), vertex_info.end(),
        [](const vertex_info_t& a, const vertex_info_t& b) {
            return a.distance > b.distance;
        });

    // Generate initial paths
    std::vector<path_t> paths;
    for (const auto& info : vertex_info) {
        if (!info.available) continue;

        // Reconstruct path from v to current vertex
        std::vector<int> vertices = reconstruct_bounded_path(v, info.vertex, bfs_map);

        // Check path length
        if (vertices.size() < 2) continue;  // Skip trivial paths

        // Get path length from endpoint's distance
        double path_length = info.distance;
        if (path_length > max_length) continue;

        // Create path_t structure
        path_t path;
        path.vertices = vertices;
        path.target_vertex = v;
        path.total_weight = path_length;
        path.center_offset = calculate_relative_center_offset(v, vertices);

        // Calculate distances to target for each vertex
        path.dist_to_target.reserve(vertices.size());
        for (int vertex : vertices) {
            auto it = bfs_map.find(vertex);
            path.dist_to_target.push_back(
                it != bfs_map.end() ? it->second.first : INFINITY
            );
        }

        // Mark vertices in path as unavailable
        for (int vertex : vertices) {
            for (auto& v_info : vertex_info) {
                if (v_info.vertex == vertex) {
                    v_info.available = false;
                    break;
                }
            }
        }

        paths.push_back(path);
    }

    // Filter for v-independence
    std::vector<path_t> independent_paths;
    for (const auto& path : paths) {
        bool is_independent = true;

        // Create vertex set excluding v
        std::unordered_set<int> path_vertices(
            path.vertices.begin(), path.vertices.end());
        path_vertices.erase(v);

        // Check against all previously accepted independent paths
        for (const auto& ind_path : independent_paths) {
            std::unordered_set<int> ind_vertices(
                ind_path.vertices.begin(), ind_path.vertices.end());
            ind_vertices.erase(v);

            // Check for intersection
            for (int vertex : path_vertices) {
                if (ind_vertices.count(vertex) > 0) {
                    is_independent = false;
                    break;
                }
            }
            if (!is_independent) break;
        }

        if (is_independent) {
            independent_paths.push_back(path);
        }
    }

    return independent_paths;
}

/**
 * @brief Finds v-independent paths with v positioned close to center
 *
 * @param initial_v_independent_paths Input vector of path_t structures from find_v_independent_paths
 * @param v Target vertex that must be in all paths
 * @param target_length Desired path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 * @param bfs_map Map from vertex to (distance, predecessor) pair from find_shortest_paths_within_radius
 * @param relative_center_offset_tolerance Tolerance for selecting paths with near-optimal center positioning
 *
 * @return Vector of path_t structures representing v-independent paths with v positioned near center
 *
 * @details The function implements a two-phase algorithm:
 *          1. For pairs of paths from initial set whose combined length is near target_length:
 *             - Finds shortest path between their endpoints
 *             - Verifies v is included and length meets criteria
 *          2. From valid paths found:
 *             - Uses existing center offset calculations
 *             - Selects paths where v is optimally centered within tolerance
 *
 * @note If input contains only one path, returns that path
 * @note Returns empty vector if no paths meeting criteria are found
 * @note All returned paths contain v and are v-independent
 *
 * @see path_t for the structure containing path and distance information
 */
static std::vector<path_t> find_centered_v_independent_paths(
    const std::vector<path_t>& initial_v_independent_paths,
    int v,
    double target_length,
    double length_tolerance,
    const std::unordered_map<int, std::pair<double, int>>& bfs_map,
    double relative_center_offset_tolerance
) {
    if (initial_v_independent_paths.empty()) return std::vector<path_t>();
    if (initial_v_independent_paths.size() == 1) return initial_v_independent_paths;

    std::vector<path_t> centered_paths;
    double length_eps = target_length * length_tolerance;

    // For each pair of paths
    for (size_t i = 0; i < initial_v_independent_paths.size(); i++) {
        for (size_t j = i + 1; j < initial_v_independent_paths.size(); j++) {
            const auto& path_i = initial_v_independent_paths[i];
            const auto& path_j = initial_v_independent_paths[j];

            // Get endpoints different from v
            int e_i = (path_i.vertices.front() == v) ?
                     path_i.vertices.back() : path_i.vertices.front();
            int e_j = (path_j.vertices.front() == v) ?
                     path_j.vertices.back() : path_j.vertices.front();

            // Reconstruct path between endpoints
            std::vector<int> vertices = reconstruct_bounded_path(e_i, e_j, bfs_map);

            // Verify v is in path
            if (std::find(vertices.begin(), vertices.end(), v) == vertices.end()) continue;

            // Get path length
            auto it = bfs_map.find(e_j);
            if (it == bfs_map.end()) continue;
            double path_length = it->second.first;

            if (std::abs(path_length - target_length) <= length_eps) {
                path_t new_path;
                new_path.vertices = vertices;
                new_path.target_vertex = v;
                new_path.total_weight = path_length;
                new_path.center_offset = calculate_relative_center_offset(v, vertices);

                // Calculate distances to target
                new_path.dist_to_target.reserve(vertices.size());
                for (int vertex : vertices) {
                    auto dist_it = bfs_map.find(vertex);
                    new_path.dist_to_target.push_back(
                        dist_it != bfs_map.end() ? dist_it->second.first : INFINITY
                    );
                }

                centered_paths.push_back(new_path);
            }
        }
    }

    if (centered_paths.empty()) return std::vector<path_t>();

    // Find minimum center offset
    double min_offset = std::numeric_limits<double>::infinity();
    for (const auto& path : centered_paths) {
        min_offset = std::min(min_offset, path.center_offset);
    }

    // Select paths with offset within tolerance
    std::vector<path_t> result;
    for (const auto& path : centered_paths) {
        if (path.center_offset <= min_offset + relative_center_offset_tolerance) {
            result.push_back(path);
        }
    }

    return result;
}

/**
 * @brief Finds v-independent paths with v positioned close to center (convenience overload)
 *
 * @param adj_list Adjacency list representation of the graph
 * @param weight_list Edge weights corresponding to adjacency list
 * @param v Target vertex that must be included in paths
 * @param target_length Desired total path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 * @param relative_center_offset_tolerance Tolerance for selecting paths with near-optimal center positioning
 *
 * @return Vector of path_t structures representing v-independent paths with v positioned near center
 *
 * @details Convenience function that combines three steps:
 *          1. Performs bounded Dijkstra search from vertex v
 *          2. Finds initial set of v-independent paths
 *          3. Selects paths with v optimally centered
 *
 * @note This overload handles the complete pipeline from graph to centered paths
 * @note The maximum search radius is automatically set to target_length + tolerance
 *
 * @see find_shortest_paths_within_radius for the bounded Dijkstra implementation
 * @see find_v_independent_paths for path independence criteria
 * @see find_centered_v_independent_paths for center positioning details
 * @see path_t for path structure definition
 */
static std::vector<path_t> find_centered_v_independent_paths(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int v,
    double target_length,
    double length_tolerance,
    double relative_center_offset_tolerance
) {
    double max_radius = target_length + length_tolerance * target_length;

    auto bfs_map = find_shortest_paths_within_radius(
        adj_list,
        weight_list,
        v,
        max_radius);

    auto initial_v_independent_paths = find_v_independent_paths(
        bfs_map,
        v,
        max_radius);

    return find_centered_v_independent_paths(
        initial_v_independent_paths,
        v,
        target_length,
        length_tolerance,
        bfs_map,
        relative_center_offset_tolerance);
}

/**
 * @brief R interface for finding v-independent paths centered around a target vertex
 *
 * @param adj_list_s SEXP List of integer vectors representing graph adjacency (1-based indices)
 * @param weight_list_s SEXP List of numeric vectors containing edge weights
 * @param vertex_s SEXP Integer scalar specifying target vertex v (1-based index)
 * @param target_length_s SEXP Numeric scalar specifying desired total path length
 * @param length_tolerance_s SEXP Numeric scalar specifying allowed deviation from target length as fraction
 * @param relative_center_offset_tolerance_s SEXP Numeric scalar specifying tolerance for path centering
 *
 * @return SEXP (list) where each element is a list containing:
 *         - vertices: integer vector of path vertices (1-based indices)
 *         - total_weight: numeric scalar of path length
 *         - center_offset: numeric scalar measuring target vertex position (0 = center)
 *         - dist_to_target: numeric vector of distances from each vertex to target v
 *
 * @details Converts R inputs to C++ types, calls find_centered_v_independent_paths(),
 *          and converts the resulting path_t structures back to R lists.
 *          Handles conversion between 0-based (C++) and 1-based (R) indexing.
 *
 * @note All vertex indices in input must be positive integers
 * @note All edge weights must be non-negative
 * @note Tolerances must be positive
 */
SEXP S_find_centered_v_independent_paths(
    SEXP adj_list_s,
    SEXP weight_list_s,
    SEXP vertex_s,
    SEXP target_length_s,
    SEXP length_tolerance_s,
    SEXP relative_center_offset_tolerance_s
) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list = Rgraph_to_vector(adj_list_s);
    std::vector<std::vector<double>> weight_list = Rweights_to_vector(weight_list_s);
    int vertex = INTEGER(vertex_s)[0] - 1;  // Convert to 0-based indexing
    double target_length = REAL(target_length_s)[0];
    double length_tolerance = REAL(length_tolerance_s)[0];
    double relative_center_offset_tolerance = REAL(relative_center_offset_tolerance_s)[0];

    // Find centered v-independent paths
    std::vector<path_t> paths = find_centered_v_independent_paths(
        adj_list,
        weight_list,
        vertex,
        target_length,
        length_tolerance,
        relative_center_offset_tolerance
    );

    // Convert results to R list
    SEXP result = PROTECT(allocVector(VECSXP, paths.size()));

    for(size_t i = 0; i < paths.size(); i++) {
        // Create list for current path (4 elements: vertices, weight, offset, distances)
        SEXP path = PROTECT(allocVector(VECSXP, 4));

        // Convert vertices to R (adding 1 for 1-based indexing)
        SEXP vertices = PROTECT(allocVector(INTSXP, paths[i].vertices.size()));
        for(size_t j = 0; j < paths[i].vertices.size(); j++) {
            INTEGER(vertices)[j] = paths[i].vertices[j] + 1;
        }

        // Create weight and offset elements
        SEXP total_weight = PROTECT(ScalarReal(paths[i].total_weight));
        SEXP center_offset = PROTECT(ScalarReal(paths[i].center_offset));

        // Convert distances to target
        SEXP dist_to_target = PROTECT(allocVector(REALSXP, paths[i].dist_to_target.size()));
        for(size_t j = 0; j < paths[i].dist_to_target.size(); j++) {
            REAL(dist_to_target)[j] = paths[i].dist_to_target[j];
        }

        // Set list elements
        SET_VECTOR_ELT(path, 0, vertices);
        SET_VECTOR_ELT(path, 1, total_weight);
        SET_VECTOR_ELT(path, 2, center_offset);
        SET_VECTOR_ELT(path, 3, dist_to_target);

        // Add path to result
        SET_VECTOR_ELT(result, i, path);

        UNPROTECT(4);  // vertices, total_weight, center_offset, dist_to_target
        UNPROTECT(1);  // path
    }

    UNPROTECT(1);  // result
    return result;
}
