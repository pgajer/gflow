#include "centered_paths.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_shortest_path.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "kernels.h"

#include <vector>
#include <queue>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>  // for std::reverse
#include <sstream>
#include <string>
#include <limits>
#include <map>
#include <unordered_set>
#include <utility>
#include <numeric> // for std::accumulate

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_get_path_data(
        SEXP adj_list_s,
        SEXP weight_list_s,
        SEXP y_s,
        SEXP ref_vertex_s,
        SEXP bandwidth_s,
        SEXP dist_normalization_factor_s,
        SEXP min_path_size_s,
        SEXP diff_threshold_s,
        SEXP kernel_type_s,
        SEXP verbose_s
        );

    SEXP S_ugg_get_path_data(
        SEXP adj_list_s,
        SEXP weight_list_s,
        SEXP grid_size_s,
        SEXP y_s,
        SEXP ref_vertex_s,
        SEXP bandwidth_s,
        SEXP dist_normalization_factor_s,
        SEXP min_path_size_s,
        SEXP diff_threshold_s,
        SEXP kernel_type_s,
        SEXP verbose_s);
}

const double EPSILON = 1e-10;  // Threshold for floating-point comparisons

/**
 * @brief Prints formatted information about vertices and their distances
 *
 * @details Displays a table of vertices and their associated distances.
 *          The output is formatted in columns for better readability.
 *          Special handling is included for infinite distances.
 *
 * @param vertex_info Vector of vertex_info_t structures containing vertex numbers
 *                    and their associated distances
 * @param name Optional name to label the output (defaults to empty string)
 *
 * @note Special cases:
 *       - Infinite distances are printed as "Inf"
 *       - If no vertices are present, prints a message indicating empty data
 *       - If name is provided, it's printed as a header
 *
 * Example output format:
 * @code
 * Vertex | Distance
 * -------|----------
 *      1 |   5.2000
 *      2 |   0.0000
 *      3 |      Inf
 * @endcode
 */
void print_vertex_info(const std::vector<vertex_info_t>& vertex_info,
                       const std::string& name = "") {
    // If empty, print message and return
    if (vertex_info.empty()) {
        Rprintf("No vertex information available.\n");
        return;
    }

    // Print name if provided
    if (!name.empty()) {
        Rprintf("%s:\n", name.c_str());
    }

    // Print header
    Rprintf("Vertex | Distance\n");
    Rprintf("-------|----------\n");

    // Print each vertex's information
    for (const auto& info : vertex_info) {
        if (std::isinf(info.distance)) {
            Rprintf("%6zu |      Inf\n", info.vertex);
        } else {
            Rprintf("%6zu | %8.4f\n", info.vertex, info.distance);
        }
    }

    Rprintf("\n");  // Add final newline for clean separation
}


/**
 * @brief Prints the contents of a vertex-to-path map showing distances and predecessors
 * @param vmap The unordered map containing vertex -> (distance, predecessor) pairs
 * @param name Optional name to label the output (defaults to empty string)
 */
void print_vertex_path_map(
    const std::unordered_map<size_t, std::pair<double, size_t>>& vmap,
    const std::string& name = ""
) {
    // First, we'll check if there's anything to print
    if (vmap.empty()) {
        Rprintf("No path information available.\n");
        return;
    }

    // If a name was provided, print it on a separate line
    if (!name.empty()) {
        Rprintf("%s:\n", name.c_str());
    }

    // Print header
    Rprintf("Vertex | Distance | Predecessor\n");
    Rprintf("-------|----------|------------\n");

    // Print each vertex's information
    for (const auto& [vertex, info] : vmap) {
        const auto& [distance, predecessor] = info;

        // Handle special case for infinity
        if (distance == std::numeric_limits<double>::infinity()) {
            Rprintf("%6zu | INF      | ", vertex);
        } else {
            Rprintf("%6zu | %8.2f | ", vertex, distance);
        }

        // Handle special case for no predecessor (e.g., start vertex)
        if (predecessor == INVALID_VERTEX) {
            Rprintf("START\n");
        } else {
            Rprintf("%zu\n", predecessor);
        }
    }
    Rprintf("\n");  // Add a blank line at the end for readability
}


/**
 * Prints the contents of a vector of paths to stdout using Rprintf
 * Each path is formatted to show its components in a clear, readable manner
 *
 * @param paths Vector of ref_vertex_path_t structures to be printed
 */
void print_vect_paths(const std::vector<ref_vertex_path_t>& paths,
                      const std::string& name = "") {

    // If a name was provided, print it on a separate line
    if (!name.empty()) {
        Rprintf("%s:\n", name.c_str());
    }

    // Print header for the output
    Rprintf("Number of paths: %d\n", static_cast<int>(paths.size()));

    // Iterate through each path
    for (size_t i = 0; i < paths.size(); ++i) {
        const ref_vertex_path_t& path = paths[i];

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
        //Rprintf("  Target Vertex:  %zu\n", path.target_vertex);
        //Rprintf("  Total Weight:   %.4f\n", path.total_weight);

        // Handle special case for infinite center offset
        // if (std::isinf(path.center_offset)) {
        //     Rprintf("  Center Offset:  Inf\n");
        // } else {
        //     Rprintf("  Center Offset:  %.4f\n", path.center_offset);
        // }

        // Add separator between paths for better readability
        if (i < paths.size() - 1) {
            Rprintf("----------------------------------------\n");
        }
    }
    Rprintf("\n");  // Final newline for clean separation from subsequent output
}

/**
 * @brief Prints a formatted representation of paths stored in a map
 *
 * @details For each path in the map, prints:
 *          - The sequence of vertices in the path
 *          - The target vertex
 *          - Total path weight
 *          - Center offset (distance of target vertex from path center)
 *          - Distances from each vertex to the target vertex (if available)
 *
 * @param paths A map where keys are vertex indices and values are ref_vertex_path_t structures
 *              containing path information including vertices, weights, and distances
 * @param name Optional name to label the output (defaults to empty string)
 *
 * @note Special cases:
 *       - Infinite values are printed as "Inf"
 *       - Paths are separated by dashed lines for readability
 *       - If name is provided, it's printed as a header
 *
 * Example output format:
 * @code
 * Path for vertex 1:
 *   Vertices:       [1, 2, 3]
 *   Target Vertex:  2
 *   Total Weight:   10.5000
 *   Center Offset:  0.0000
 *   Distances to target:
 *     vertex 1 -> target: 5.2000
 *     vertex 2 -> target: 0.0000
 *     vertex 3 -> target: 5.3000
 * @endcode
 */
void print_map_paths(const std::map<size_t, ref_vertex_path_t>& paths,
                     const std::string& name = "") {
    // If a name was provided, print it on a separate line
    if (!name.empty()) {
        Rprintf("%s:\n", name.c_str());
    }

    // Print header for the output
    Rprintf("Number of paths: %zu\n", paths.size());

    // Iterate through each path in the map
    for (const auto& [vertex_idx, path] : paths) {
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
        Rprintf("\nPath for vertex %zu:\n", vertex_idx);
        Rprintf("  Vertices:       %s\n", vertices_str.str().c_str());
        // Rprintf("  Target Vertex:  %zu\n", path.target_vertex);
        // Rprintf("  Total Weight:   %.4f\n", path.total_weight);

        // // Handle special case for infinite center offset
        // if (std::isinf(path.center_offset)) {
        //     Rprintf("  Center Offset:  Inf\n");
        // } else {
        //     Rprintf("  Center Offset:  %.4f\n", path.center_offset);
        // }

        // Print distances to target if available
        if (!path.dist_to_ref_vertex.empty()) {
            Rprintf("  Distances to target:\n");
            for (size_t i = 0; i < path.dist_to_ref_vertex.size(); ++i) {
                if (std::isinf(path.dist_to_ref_vertex[i])) {
                    Rprintf("    vertex %zu -> target: Inf\n", path.vertices[i]);
                } else {
                    Rprintf("    vertex %zu -> target: %.4f\n",
                           path.vertices[i], path.dist_to_ref_vertex[i]);
                }
            }
        }

        // Add separator between paths for better readability
        if (std::next(paths.find(vertex_idx)) != paths.end()) {
            Rprintf("----------------------------------------\n");
        }
    }
    Rprintf("\n");  // Final newline for clean separation from subsequent output
}

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
static double calculate_relative_center_offset(size_t v, const std::vector<size_t>& path) {
    auto it = std::find(path.begin(), path.end(), v);
    if (it == path.end()) return INFINITY;

    double path_position = std::distance(path.begin(), it);
    return std::abs(path_position / (path.size() - 1) - 0.5);
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
 * @return std::unordered_map<size_t, std::pair<double, size_t>> Map from vertex to its
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
static std::unordered_map<size_t, std::pair<double, size_t>> find_shortest_paths_within_radius(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    size_t start,
    double radius
    ) {

    std::unordered_map<size_t, std::pair<double, size_t>> result; // result[vertex] = <distance, prev>
    size_t n = adj_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    dist[start] = 0;

    // debugging ------
    #if 0
    Rprintf("In find_shortest_paths_within_radius()\n");
    print_vect_vect(adj_list, "adj_list");
    print_vect_vect(weight_list, "weight_list");
    Rprintf("start: %d\tradius: %f\n", start, radius);
    #endif

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto top = pq.top();
        double d = -top.first;
        size_t u = top.second;
        pq.pop();

        if (d > radius) break;
        if (d > dist[u]) continue;

        result.emplace(u, std::pair<double, size_t>(d, prev[u]));

        for (size_t i = 0; i < adj_list[u].size(); i++) {

            // debugging
            // Rprintf("i: %zu\tu: %zu\n", i, u);
            size_t v = adj_list[u][i];
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

/**
 * @brief Reconstructs unique paths from vertices to a target vertex using BFS results
 *
 * @details Processes BFS map results to create ref_vertex_path_t objects containing:
 *          - Complete vertex sequences for each path
 *          - Total path weights
 *          - Center offsets
 *          - Distances to target vertex
 *
 * @param ref_vertex Target vertex that all paths lead to
 * @param bfs_map Map of vertex -> (distance, predecessor) pairs from BFS results
 * @return std::vector<ref_vertex_path_t> Vector of reconstructed unique paths
 *
 * @note Algorithm:
 *       1. Orders vertices by distance from target
 *       2. Reconstructs paths starting from furthest vertices
 *       3. Skips vertices already included in previous paths
 */
static std::vector<ref_vertex_path_t> reconstruct_ref_vertex_paths(
    size_t ref_vertex,
    const std::unordered_map<size_t, std::pair<double, size_t>>& bfs_map
) {
    std::vector<ref_vertex_path_t> result;

    // Create and sort vertex info
    std::vector<vertex_info_t> vertex_info;
    vertex_info.reserve(bfs_map.size());
    for (const auto& [vertex, info] : bfs_map) {
        if (vertex != ref_vertex) {  // exclude target vertex itself
            //vertex_info.emplace_back(vertex, info.first);
            vertex_info.push_back({vertex, info.first});
        }
    }

    // Sort by distance in descending order
    std::sort(vertex_info.begin(), vertex_info.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    // Keep track of used vertices to avoid Rf_duplicate paths
    std::unordered_set<size_t> used_vertices;
    used_vertices.insert(ref_vertex);  // target vertex is always used

    // Process vertices in order of decreasing distance
    for (const auto& info : vertex_info) {
        // Skip if this vertex was already used in another path
        if (used_vertices.count(info.vertex) > 0) {
            continue;
        }

        ref_vertex_path_t new_path;
        // new_path.target_vertex = ref_vertex;
        // new_path.total_weight = info.distance;  // distance from start to target

        // Reconstruct the path from this vertex to target
        size_t curr = info.vertex;
        while (curr != INVALID_VERTEX) {
            new_path.vertices.push_back(curr);
            used_vertices.insert(curr);  // mark vertex as used

            auto it = bfs_map.find(curr);
            if (it == bfs_map.end()) {
                Rf_error("Path reconstruction failed: vertex not found in path info");
            }

            // Store distance to target for this vertex
            new_path.dist_to_ref_vertex.push_back(it->second.first);
            curr = it->second.second;  // move to predecessor
        }

        // Reverse the path to get correct order (from target to end)
        std::reverse(new_path.vertices.begin(), new_path.vertices.end());
        std::reverse(new_path.dist_to_ref_vertex.begin(), new_path.dist_to_ref_vertex.end());

        // Calculate center offset
        if (new_path.vertices.size() > 1) {
            new_path.center_offset = calculate_relative_center_offset(ref_vertex, new_path.vertices);
        }

        result.push_back(new_path);
    }

    // Sort paths by total weight and center offset (using ref_vertex_path_t's operator<)
    std::sort(result.begin(), result.end());

    return result;
}

// ------------------------------------------------------------------------------
//
// implementation of find_paths_meeting_size_requirement()
//
// ------------------------------------------------------------------------------

/**
 * @brief Computes and populates all metrics for a path including distances, weights, and kernel values
 *
 * @details Calculates and sets the following metrics in the input path_data:
 *          - x_path: distances along the path from initial vertex
 *          - w_path: kernel-weighted values based on distance from reference vertex
 *          - y_path: y-values for each vertex in the path
 *          - total_weight: sum of edge weights along the path
 *          - rel_center_offset: relative position of reference vertex in path
 *
 * @param[in,out] path_data Path data structure to be populated with metrics
 * @param[in] adj_list Graph adjacency list
 * @param[in] weight_list Edge weights corresponding to adjacency list
 * @param[in] y Vector of y-values for each vertex
 * @param[in] dist_normalization_factor Factor for normalizing distances in kernel calculations
 *
 * @pre path_data.vertices and path_data.ref_vertex must be properly initialized
 * @pre adj_list and weight_list must represent a valid graph structure
 * @pre y must contain values for all vertices in the path
 *
 * @throws std::runtime_error if reference vertex is not found in path vertices
 */
static void compute_path_metrics(
    path_data_t& path_data,
    const std::vector<double>& y,
    double dist_normalization_factor
    ) {

    if (path_data.x_path.size() != path_data.vertices.size()) {
        REPORT_ERROR("x_path size mismatch with vertices size");
    }

    size_t path_n_vertices = path_data.vertices.size();
    path_data.w_path.resize(path_n_vertices);
    path_data.y_path.resize(path_n_vertices);
    std::vector<double> d_path(path_n_vertices);

    // Calculate y_path
    for (size_t k = 0; k < path_n_vertices; ++k) {
        path_data.y_path[k] = y[path_data.vertices[k]];
    }

    // Calculate weights using kernel
    //if (path_data.ref_vertex_index < 0) {
    // path_data.ref_vertex_index = get_ref_vertex_index(path_data.ref_vertex, path_data.vertices);
    //}
    double x_ref = 0; //path_data.x_path[path_data.ref_vertex_index];
    double max_dist = 0.0;

    for (size_t k = 0; k < path_n_vertices; ++k) {
        d_path[k] = std::abs(path_data.x_path[k] - x_ref);
        max_dist = std::max(max_dist, d_path[k]);
    }

    if (max_dist == 0) max_dist = 1;
    max_dist *= dist_normalization_factor;

    for (size_t k = 0; k < path_n_vertices; ++k) {
        d_path[k] /= max_dist;
    }

    kernel_fn(d_path.data(), path_n_vertices, path_data.w_path.data());

    // Normalize weights
    double total_w_path = std::accumulate(path_data.w_path.begin(), path_data.w_path.end(), 0.0);
    for (size_t k = 0; k < path_n_vertices; ++k) {
        path_data.w_path[k] /= total_w_path;
    }

    path_data.rel_center_offset = calculate_relative_center_offset(path_data.ref_vertex, path_data.vertices); // may be easier to compute it using  path_data.ref_vertex_index and path_data.vertices.size() <<---
}


/**
 * @brief Reconstructs complete paths from current vertex information gathered during BFS
 *
 * @details For each vertex discovered during BFS that hasn't been used in a previous path,
 *          reconstructs a complete path from that vertex back to the start vertex using
 *          predecessor information. The paths are constructed starting from vertices
 *          farthest from the start vertex.
 *
 * @param[in] vertex_info Map containing for each vertex:
 *                        - distance from start (get<0>)
 *                        - predecessor vertex (get<1>)
 *                        - path length to this vertex (get<2>)
 * @param[in] start Starting vertex that all paths lead to
 *
 * @return Vector of ref_vertex_path_t objects, each containing:
 *         - vertices: sequence of vertices in the path
 *         - dist_to_ref_vertex: distances from each vertex to the target
 *         - total_weight: total path weight
 *         - target_vertex: the start vertex
 *
 * @pre vertex_info must contain valid predecessor information that leads back to start vertex
 * @pre vertex_info must contain accurate distance information for path sorting
 *
 * @note Vertices are processed in order of decreasing distance from start to ensure
 *       longer paths are constructed first
 * @note Each vertex (except start) is used in at most one path to avoid duplicates
 *
 * @throws std::runtime_error if path reconstruction fails due to missing vertex information
 */
static std::vector<ref_vertex_path_t> reconstruct_current_paths(
    const std::unordered_map<size_t, std::tuple<double, size_t, size_t>>& vertex_info,
    size_t start
) {
    std::vector<ref_vertex_path_t> paths;
    std::unordered_set<size_t> used_vertices;
    used_vertices.insert(start);

    // Create ordered vertex info for processing
    std::vector<std::pair<size_t, double>> vertices;
    for (const auto& [vertex, info] : vertex_info) {
        if (vertex != start) {
            vertices.emplace_back(vertex, std::get<0>(info));
        }
    }

    // Sort by distance
    std::sort(vertices.begin(), vertices.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Reconstruct paths from farthest vertices
    for (const auto& [vertex, distance] : vertices) {
        if (used_vertices.count(vertex) > 0) continue;

        ref_vertex_path_t new_path;
        new_path.target_vertex = start;
        new_path.total_weight = distance;

        // Build path from this vertex back to start
        size_t curr = vertex;
        while (curr != INVALID_VERTEX) {
            new_path.vertices.push_back(curr);
            used_vertices.insert(curr);

            auto it = vertex_info.find(curr);
            if (it == vertex_info.end()) {
                REPORT_ERROR("Vertex not found in path info during reconstruction");
            }

            new_path.dist_to_ref_vertex.push_back(std::get<0>(it->second));
            curr = std::get<1>(it->second);  // Get predecessor
        }

        std::reverse(new_path.vertices.begin(), new_path.vertices.end());
        std::reverse(new_path.dist_to_ref_vertex.begin(), new_path.dist_to_ref_vertex.end());

        paths.push_back(std::move(new_path));
    }

    return paths;
}


/**
 * @brief Constructs a composite path by joining two paths through their common reference vertex
 *
 * @details Modifies the input composite_path by:
 *          1. Taking the reverse of path2
 *          2. Joining it with path1 (excluding Rf_duplicate reference vertex)
 *          3. Computing all necessary path metrics
 *          4. Adjusting total weight to account for shared edge at reference vertex
 *
 * @param[out] composite_path Output structure to store the constructed path
 * @param[in] path1 First path to join
 * @param[in] path2 Second path to join (will be reversed)
 * @param[in] adj_list Graph adjacency list
 * @param[in] weight_list Edge weights corresponding to adjacency list
 * @param[in] y Vector of y-values for each vertex
 * @param[in] dist_normalization_factor Factor for normalizing distances in kernel calculations
 *
 * @pre path1 and path2 must share composite_path.ref_vertex
 * @pre adj_list and weight_list must represent a valid graph structure
 * @pre y must contain values for all vertices in both paths
 *
 * @note The function will compute all path metrics (x_path, w_path, y_path) for the
 *       constructed composite path
 */
static void construct_composite_path(
    path_data_t& composite_path,
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    const std::vector<double>& y,
    double dist_normalization_factor
) {
    // Reserve space for the composite path
    size_t combined_length = path1.vertices.size() + path2.vertices.size();
    if (combined_length < path1.vertices.size()) { // Check for overflow
        REPORT_ERROR("Path length overflow detected");
    }
    combined_length -= 1; // Now safe to subtract

    composite_path.vertices.reserve(combined_length);


    #if 0
    Rprintf("\nIn construct_composite_path()\n");
    auto print_ref_vertex_path_t = [](const ref_vertex_path_t& path, const std::string& name = "") {
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
        Rprintf("\nPath %s:\n", name.c_str());  // Fixed: using name instead of names
        Rprintf("  Vertices:       %s\n", vertices_str.str().c_str());
        Rprintf("  Target Vertex:  %d\n", path.target_vertex);
        Rprintf("  Total Weight:   %.4f\n", path.total_weight);
    };  // Note the semicolon at the end of the lambda definition

    print_ref_vertex_path_t(path1, "path1");
    print_ref_vertex_path_t(path2, "path2");
    #endif

    // Add reversed path2
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path2.vertices.rbegin(),
        path2.vertices.rend()
    );

    // Add path1 (excluding the reference vertex)
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path1.vertices.begin() + 1,
        path1.vertices.end()
    );

    // print_vect(composite_path.vertices,"composite_path.vertices");

    // Calculate total weight
    composite_path.total_weight = path1.total_weight - path1.dist_to_ref_vertex[1] + path2.total_weight;

    // Given distances from the start vertex of path1 and path2 as recored in
    // dist_to_ref_vertex component of ref_vertex_path_t, we can reconstruct x_path from these
    composite_path.x_path.resize(combined_length);
    composite_path.x_path[0] = 0.0;
    size_t n_path2 = path2.vertices.size();
    for (size_t i = 1; i < n_path2; i++) {
        composite_path.x_path[i] = path2.total_weight - path2.dist_to_ref_vertex[n_path2 - i];
    }

    double x2 = composite_path.x_path[n_path2 - 1];
    for (size_t i = n_path2; i < combined_length; i++) {
        composite_path.x_path[i] = x2 + path1.dist_to_ref_vertex[n_path2 - i + 1];
    }

    // Compute all path metrics (w_path, y_path, etc.)
    compute_path_metrics(
        composite_path,
        y,
        dist_normalization_factor
    );
}

/**
 * @brief Finds paths that meet minimum size requirement, either as single paths or composites
 *
 * @param adj_list Adjacency list representation of the graph
 * @param weight_list Corresponding edge weights
 * @param start Starting vertex
 * @param min_path_size Minimum required number of vertices in resulting path(s)
 * @param y Vector of y-values for each vertex
 * @param dist_normalization_factor Factor for normalizing distances in kernel calculations
 * @return std::vector<path_data_t> Vector of path data objects meeting size requirement
 */
static std::vector<path_data_t> find_paths_meeting_size_requirement(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    size_t ref_vertex,
    size_t min_path_size,
    size_t diff_threshold,
    const std::vector<double>& y,
    double dist_normalization_factor
) {
    std::vector<path_data_t> result;
    size_t n = adj_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    dist[ref_vertex] = 0;

    auto paths_explore_different_directions = [diff_threshold](const ref_vertex_path_t& path1, const ref_vertex_path_t& path2) -> bool {
        size_t remaining1 = path1.vertices.size() - 1;
        size_t remaining2 = path2.vertices.size() - 1;
        size_t check_length = std::min({diff_threshold, remaining1, remaining2});

        // Compare initial segments of both paths
        for (size_t i = 1; i <= check_length; i++) {
            // Look for path1's current vertex in path2's initial segment
            auto it = std::find(path2.vertices.begin() + 1,
                                path2.vertices.begin() + 1 + check_length,
                                path1.vertices[i]);

            // If we found this vertex in path2, paths are not different
            if (it != (path2.vertices.begin() + 1 + check_length)) {
                return false;
            }
        }

        // If we get here, no shared vertices were found
        return true;
    };

    struct queue_entry_t {
        double distance;
        size_t vertex;
        size_t path_length;

        bool operator>(const queue_entry_t& other) const {
            return distance > other.distance;
        }
    };

    std::priority_queue<queue_entry_t, std::vector<queue_entry_t>, std::greater<>> pq;
    pq.push({0, ref_vertex, 1});

    std::unordered_map<size_t, std::tuple<double, size_t, size_t>> vertex_info;

    while (!pq.empty()) {
        auto [d, u, path_len] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        vertex_info[u] = std::make_tuple(d, prev[u], path_len);

        // In find_paths_meeting_size_requirement, in the single path case:
        if (path_len >= min_path_size) {
            path_data_t path_data;
            path_data.ref_vertex = ref_vertex;
            //path_data.ref_vertex_index = 0;

            // Reconstruct path and distances
            std::vector<double> distances;
            size_t curr = u;
            while (curr != INVALID_VERTEX) {
                path_data.vertices.push_back(curr);
                distances.push_back(dist[curr]);  // Store distances
                curr = std::get<1>(vertex_info[curr]);
            }

            // Reverse both vectors since we built them backwards
            std::reverse(path_data.vertices.begin(), path_data.vertices.end());
            std::reverse(distances.begin(), distances.end());

            // Set x_path - these are distances from start
            path_data.x_path = std::move(distances);
            path_data.total_weight = dist[u];

            compute_path_metrics(path_data,
                                 y,
                                 dist_normalization_factor);
            result.push_back(std::move(path_data));
            return result;
        }

        // Explore neighbors
        for (size_t i = 0; i < adj_list[u].size(); i++) {
            size_t v = adj_list[u][i];
            double w = weight_list[u][i];
            double new_dist = dist[u] + w;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                pq.push({new_dist, v, path_len + 1});
            }
        }

        // Check for composite paths
        if (!vertex_info.empty() && path_len > 1) {
            auto paths = reconstruct_current_paths(vertex_info, ref_vertex);

            for (size_t i = 0; i < paths.size(); ++i) {
                for (size_t j = i + 1; j < paths.size(); ++j) {

                    if (!paths_explore_different_directions(paths[i], paths[j])) continue;

                    size_t combined_length = paths[i].vertices.size() + paths[j].vertices.size() - 1;
                    if (combined_length >= min_path_size) {
                        path_data_t composite_path;
                        composite_path.ref_vertex = ref_vertex;
                        //composite_path.ref_vertex_index = paths[j].vertices.size() - 1; // after reversing paths[j].vertices; ref_vertex is the last vertex of the path

                        // Construct composite path
                        construct_composite_path(
                            composite_path,
                            paths[i],
                            paths[j],
                            y,
                            dist_normalization_factor
                        );

                        result.push_back(std::move(composite_path));
                        return result;
                    }
                }
            }
        }
    }

    return result;
}

/**
 * @brief Finds paths through a reference vertex and computes their metrics
 *
 * @details First attempts to find paths within a given bandwidth radius. If no
 *          suitable paths are found (either single path or composite paths meeting
 *          minimum size requirement), expands search without radius constraint.
 *          For composite paths, ensures that the initial segments (up to diff_threshold
 *          vertices after ref_vertex) are completely different, thus capturing different
 *          local directions in the graph.
 *
 * @param[in] adj_list Graph adjacency list
 * @param[in] weight_list Edge weights corresponding to adjacency list
 * @param[in] y Vector of y-values for each vertex
 * @param[in] ref_vertex Reference vertex all paths must contain
 * @param[in] bandwidth Initial radius for bounded path search
 * @param[in] dist_normalization_factor Factor for normalizing distances in kernel calculations
 * @param[in] min_path_size Minimum required number of vertices in path
 * @param[in] diff_threshold Number of vertices after ref_vertex that must be different
 *                          for paths to be considered exploring different directions
 * @param[in] kernel_type Type of kernel function to use
 * @param[in] verbose Whether to print debug information
 *
 * @return Vector of path_data_t objects meeting size requirements
 *
 * @note If no paths meeting size requirements are found within bandwidth,
 *       automatically expands search without radius constraint.
 *       Composite paths are only created from pairs of paths that explore different
 *       directions, as determined by having no common vertices in their first
 *       diff_threshold vertices after ref_vertex.
 */
static std::vector<path_data_t> get_path_data(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    bool verbose
) {
    initialize_kernel(kernel_type, 1.0);

    // First try with bandwidth constraint
    auto shortest_paths_map = find_shortest_paths_within_radius(
        adj_list,
        weight_list,
        ref_vertex,
        bandwidth
    );

    std::vector<ref_vertex_path_t> paths = reconstruct_ref_vertex_paths(ref_vertex, shortest_paths_map);
    size_t n_paths = paths.size();

    if (verbose) {
        // debugging
        Rprintf("n_paths: %zu\n", n_paths);
        print_vect_paths(paths, "reconstructed paths from shortest_paths_map");
    }

    auto paths_explore_different_directions = [diff_threshold](const ref_vertex_path_t& path1, const ref_vertex_path_t& path2) -> bool {
        size_t remaining1 = path1.vertices.size() - 1;
        size_t remaining2 = path2.vertices.size() - 1;
        size_t check_length = std::min({diff_threshold, remaining1, remaining2});

        // Compare initial segments of both paths
        for (size_t i = 1; i <= check_length; i++) {
            // Look for path1's current vertex in path2's initial segment
            auto it = std::find(path2.vertices.begin() + 1,
                                path2.vertices.begin() + 1 + check_length,
                                path1.vertices[i]);

            // If we found this vertex in path2, paths are not different
            if (it != (path2.vertices.begin() + 1 + check_length)) {
                return false;
            }
        }

        // If we get here, no shared vertices were found
        return true;
    };

    // Check if paths within bandwidth are sufficient
    bool needs_expansion = false;

    if (n_paths < 1) {
        needs_expansion = true;
        if (needs_expansion && verbose) {
            REPORT_WARNING(
                "Warning: No paths were found in the local nbhd with radius %f\n",
                bandwidth
            );
        }
    } else if (n_paths == 1) {
        needs_expansion = paths[0].vertices.size() < min_path_size;
        if (needs_expansion && verbose) {
            REPORT_WARNING(
                "Warning: Local nbhd has only one path and its length is %zu - that is less than min_path_size: %zu\n"
                "Expanded local neighborhood so there is at least composite path with min_path_size vertices.\n",
                paths[0].vertices.size(), min_path_size
            );
        }
    } else if (n_paths > 1) {
        size_t max_composite_size = 0;
        for (size_t i = 0; i < n_paths; ++i) {
            for (size_t j = i + 1; j < n_paths; ++j) {
                if (paths_explore_different_directions(paths[i], paths[j])) {
                    max_composite_size = std::max(
                        max_composite_size,
                        paths[i].vertices.size() + paths[j].vertices.size() - 1
                        );
                }
            }
        }
        needs_expansion = max_composite_size < min_path_size;
        if (needs_expansion && verbose) {
            REPORT_WARNING(
                "Warning: None of the composed paths has min_path_size: %zu vertices\n"
                "Expanded local neighborhood so there is at least composite path with min_path_size vertices.\n",
                min_path_size
            );
        }
    }

    // If current paths don't meet requirements, expand search
    if (needs_expansion) {
        return find_paths_meeting_size_requirement(
            adj_list,
            weight_list,
            ref_vertex,
            min_path_size,
            diff_threshold,
            y,
            dist_normalization_factor
        );
    }

    // Process paths within bandwidth that meet requirements
    std::vector<path_data_t> result;

    if (n_paths == 1) {
        path_data_t path_data;
        path_data.ref_vertex = ref_vertex;
        //path_data.ref_vertex_index = 0;
        path_data.vertices = std::move(paths[0].vertices);
        path_data.x_path   = std::move(paths[0].dist_to_ref_vertex);
        path_data.total_weight = paths[0].total_weight;

        compute_path_metrics(
            path_data,
            y,
            dist_normalization_factor
            );

        result.push_back(std::move(path_data));
    } else {
        for (size_t i = 0; i < n_paths; ++i) {
            for (size_t j = i + 1; j < n_paths; ++j) {
                size_t comb_path_size = paths[i].vertices.size() + paths[j].vertices.size() - 1;
                if (comb_path_size < min_path_size) continue;

                if (!paths_explore_different_directions(paths[i], paths[j])) continue;

                if (verbose) {
                    Rprintf("i: %zu  j: %zu\n", i, j);
                    print_vect(paths[i].vertices, "paths[i].vertices");
                    print_vect(paths[j].vertices, "paths[j].vertices");
                }

                path_data_t path_data;
                path_data.ref_vertex = ref_vertex;

                construct_composite_path(
                    path_data,
                    paths[i],
                    paths[j],
                    y,
                    dist_normalization_factor
                );

                result.push_back(std::move(path_data));
            }
        }
    }

    return result;
}


/**
 * @brief R interface for computing path data centered around a reference vertex
 *        with associated distances, weights, and values along the paths
 *
 * @param adj_list_s SEXP List of integer vectors representing graph adjacency (1-based indices)
 *                   Must be a valid graph representation with matching dimensions
 * @param weight_list_s SEXP List of numeric vectors containing edge weights
 *                     Must match adj_list_s structure exactly
 * @param y_s SEXP Numeric vector containing values associated with each vertex
 * @param ref_vertex_s SEXP Integer scalar specifying reference vertex (1-based index)
 *                    Must be within valid range for the graph
 * @param bandwidth_s SEXP Numeric scalar specifying maximum path distance from reference vertex
 * @param dist_normalization_factor_s SEXP Numeric scalar for distance normalization in kernel computations
 * @param min_path_size_s SEXP Integer scalar specifying minimum required path length
 * @param kernel_type_s SEXP Integer scalar specifying kernel function type (1-10)
 *
 * @return SEXP (list) where each element is a list with named components:
 *         - "vertices": integer vector of path vertices (1-based indices)
 *         - "ref_vertex": integer scalar of reference vertex (1-based index)
 *         - "rel_center_offset": numeric scalar measuring reference vertex position (0 = center, 0.5 = endpoint)
 *         - "total_weight": numeric scalar of path length
 *         - "x_path": numeric vector of cumulative distances along path from ref_vertex
 *         - "w_path": numeric vector of kernel weights for each vertex
 *         - "y_path": numeric vector of y-values for path vertices
 *
 * @note All vertex indices in input must be positive integers
 * @note All edge weights must be non-negative
 * @note Kernel type must be between 1 and 10 (see MSR2_KERNELS.H for definitions)
 */
SEXP S_get_path_data(
    SEXP adj_list_s,
    SEXP weight_list_s,
    SEXP y_s,
    SEXP ref_vertex_s,
    SEXP bandwidth_s,
    SEXP dist_normalization_factor_s,
    SEXP min_path_size_s,
    SEXP diff_threshold_s,
    SEXP kernel_type_s,
    SEXP verbose_s
    ) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(adj_list_s);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_s);

    size_t n_vertices = LENGTH(y_s);
    std::vector<double> y(REAL(y_s), REAL(y_s) + n_vertices);

    size_t ref_vertex = INTEGER(ref_vertex_s)[0];
    double bandwidth = REAL(bandwidth_s)[0];
    double dist_normalization_factor = REAL(dist_normalization_factor_s)[0];
    size_t min_path_size = INTEGER(min_path_size_s)[0];
    size_t diff_threshold = INTEGER(diff_threshold_s)[0];
    size_t kernel_type = INTEGER(kernel_type_s)[0];
    bool verbose = (LOGICAL(verbose_s)[0] == 1);

    // print_vect_vect(adj_list,"adj_list");

    std::vector<path_data_t> paths = get_path_data(
        adj_list,
        weight_list,
        y,
        ref_vertex,
        bandwidth,
        dist_normalization_factor,
        min_path_size,
        diff_threshold,
        kernel_type,
        verbose
        );

    size_t n_paths = paths.size();

    // Convert results to R list
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_paths));

    // Define component names matching path_data_t structure
    const std::vector<std::string> path_comps_names = {
        "vertices", "ref_vertex", "rel_center_offset", "total_weight",
        "x_path", "w_path", "y_path"
    };

    // Create SEXP for names once
    SEXP names = PROTECT(Rf_allocVector(STRSXP, path_comps_names.size()));
    for(size_t i = 0; i < path_comps_names.size(); i++) {
        SET_STRING_ELT(names, i, Rf_mkChar(path_comps_names[i].c_str()));
    }

    for(size_t i = 0; i < n_paths; i++) {
        // Create list for current path (7 elements matching path_data_t)
        SEXP path = PROTECT(Rf_allocVector(VECSXP, 7));

        // Convert vertices to R (adding 1 for 1-based indexing)
        SEXP vertices = PROTECT(Rf_allocVector(INTSXP, paths[i].vertices.size()));
        for(size_t j = 0; j < paths[i].vertices.size(); j++) {
            INTEGER(vertices)[j] = paths[i].vertices[j] + 1;
        }

        // Create scalar elements
        SEXP ref_vertex_r = PROTECT(Rf_ScalarInteger(paths[i].ref_vertex + 1));
        SEXP rel_center_offset = PROTECT(Rf_ScalarReal(paths[i].rel_center_offset));
        SEXP total_weight = PROTECT(Rf_ScalarReal(paths[i].total_weight));

        // Convert vector elements
        SEXP x_path = PROTECT(Rf_allocVector(REALSXP, paths[i].x_path.size()));
        SEXP w_path = PROTECT(Rf_allocVector(REALSXP, paths[i].w_path.size()));
        SEXP y_path = PROTECT(Rf_allocVector(REALSXP, paths[i].y_path.size()));

        // Copy vector data
        std::copy(paths[i].x_path.begin(), paths[i].x_path.end(), REAL(x_path));
        std::copy(paths[i].w_path.begin(), paths[i].w_path.end(), REAL(w_path));
        std::copy(paths[i].y_path.begin(), paths[i].y_path.end(), REAL(y_path));

        // Set list elements
        SET_VECTOR_ELT(path, 0, vertices);
        SET_VECTOR_ELT(path, 1, ref_vertex_r);
        SET_VECTOR_ELT(path, 2, rel_center_offset);
        SET_VECTOR_ELT(path, 3, total_weight);
        SET_VECTOR_ELT(path, 4, x_path);
        SET_VECTOR_ELT(path, 5, w_path);
        SET_VECTOR_ELT(path, 6, y_path);

        // Set names for the current path list
        Rf_setAttrib(path, R_NamesSymbol, names);

        // Add path to result
        SET_VECTOR_ELT(result, i, path);

        // Unprotect all elements created in this iteration
        UNPROTECT(7);  // vertices, ref_vertex_r, rel_center_offset, total_weight, x_path, w_path, y_path
        UNPROTECT(1);  // path
    }

    UNPROTECT(2);  // result, names
    return result;
}



// --------------------------------------------------------------------------------------------------------
//
// Implementation of ugg_get_path_data()
//
// --------------------------------------------------------------------------------------------------------

// Helper to check if a vertex is from the original graph (has y-value)
bool is_original_vertex(size_t vertex, const std::vector<double>& y) {
    return vertex < y.size();
}

bool is_original_vertex(size_t vertex, size_t n_original_vertices) {
    return vertex < n_original_vertices;
}

// Get neighbors and weights from uniform grid graph
std::vector<std::pair<size_t, double>> get_neighbors(
    const uniform_grid_graph_t& graph,
    size_t vertex
) {
    std::vector<std::pair<size_t, double>> neighbors;
    for (const auto& [v, w] : graph.adjacency_list[vertex]) {
        neighbors.push_back({v, w});
    }
    return neighbors;
}

/**
 * @brief Calculates the relative center offset of a grid vertex with respect to a path,
 *        accounting for its precise position between path vertices
 *
 * @details
 * This function determines how far a grid vertex is from the center of a path by finding
 * its exact position relative to the path vertices. The function implements a continuous
 * measure of centrality that works both when the grid vertex coincides with a path vertex
 * and when it lies between path vertices.
 *
 * The calculation follows these steps:
 * 1. Compute shortest distances from the grid vertex to all reachable vertices
 * 2. Find the consecutive pair of path vertices that minimize total distance to the grid vertex
 * 3. Calculate the grid vertex's relative position between these vertices
 * 4. Convert this position to a centrality measure
 *
 * The centrality measure is calculated as:
 * @code
 * abs((idx_before + offset)/(path_size - 1) - 0.5)
 * @endcode
 * where:
 * - idx_before is the index of the path vertex before the grid vertex
 * - offset = d_before/(d_before + d_after) represents the relative position between vertices
 * - d_before is the distance to the path vertex before the grid vertex
 * - d_after is the distance to the next path vertex
 *
 * Special cases:
 * - When the grid vertex coincides with a path vertex, offset = 0
 * - When the grid vertex is at an endpoint, the function returns 0.5
 * - When the grid vertex lies on the first edge, idx_before = 0
 * - When the grid vertex lies on the last edge, idx_before = path_size - 2
 *
 * Example:
 * For a path [A,B,C,D,E] and a grid vertex G between B and C:
 * @code
 * If: distances[B] = 1.2, distances[C] = 0.8
 * Then: idx_before = 1 (B's index)
 *       offset = 1.2/(1.2 + 0.8) = 0.6
 *       result = |(1 + 0.6)/4 - 0.5| = |0.4 - 0.5| = 0.1
 * @endcode
 *
 * @param grid_vertex The reference vertex from the uniform grid graph
 * @param path_vertices Vector of vertices forming the path (contains only original graph vertices)
 * @param graph The uniform grid graph structure containing connectivity information
 *
 * @return
 * - A value between 0 and 0.5, where:
 *   - 0.0 indicates the grid vertex is at the exact center of the path
 *   - 0.5 indicates the grid vertex is at one of the endpoints
 *   - Values between 0 and 0.5 indicate relative distance from center
 * - Returns INFINITY if no valid position can be found
 *
 * @pre
 * - path_vertices must contain at least two vertices
 * - grid_vertex must be a valid vertex in the uniform grid graph
 * - path_vertices must contain only vertices from the original graph
 *
 * @note
 * The function uses a breadth-first search to compute distances, ensuring accurate
 * results even when the grid vertex is not directly connected to path vertices.
 * Floating-point comparisons use an EPSILON constant to handle numerical precision issues.
 *
 * @Rf_warning
 * Performance depends on the size of the grid graph, as BFS must explore all vertices
 * reachable from the grid vertex to compute accurate distances.
 */
double calculate_grid_relative_center_offset(
    size_t grid_vertex,
    const std::vector<size_t>& path_vertices,
    const uniform_grid_graph_t& graph
) {
    // Compute distances as before
    std::unordered_map<size_t, double> distances;
    std::queue<std::pair<size_t, double>> q;
    q.push({grid_vertex, 0.0});

    // fprintf(stderr,"\nIn calculate_grid_relative_center_offset()\n");
    // fprintf(stderr,"grid_vertex: %zu\n", grid_vertex);
    // print_vect(path_vertices, "path_vertices");

    while (!q.empty()) {
        auto [curr, dist] = q.front();
        q.pop();
        if (distances.count(curr)) continue;

        // fprintf(stderr,"curr: %zu\n", curr);
        //distances[curr] = dist;
        distances.try_emplace(curr, dist);

        for (const auto& [next, weight] : graph.adjacency_list[curr]) {
            if (!distances.count(next)) {
                q.push({next, dist + weight});
            }
        }
    }

    // Find two consecutive vertices with minimum total distance to grid_vertex
    double min_total_dist = INFINITY;
    size_t idx_before = INVALID_VERTEX;
    double d_before, d_after;

    for (size_t i = 0; i < path_vertices.size() - 1; i++) {
        if (!distances.count(path_vertices[i]) ||
            !distances.count(path_vertices[i+1])) continue;

        double dist1 = distances[path_vertices[i]];
        double dist2 = distances[path_vertices[i+1]];
        double total_dist = dist1 + dist2;

        if (total_dist < min_total_dist) {
            min_total_dist = total_dist;
            idx_before = i;
            d_before = dist1;
            d_after = dist2;
        }
    }

    // Special cases for endpoints
    if (idx_before == INVALID_VERTEX) {
        // Check if grid_vertex coincides with an endpoint
        if (distances.count(path_vertices.front()) &&
            distances[path_vertices.front()] < EPSILON) {
            return 0.5;  // At first vertex
        }
        if (distances.count(path_vertices.back()) &&
            distances[path_vertices.back()] < EPSILON) {
            return 0.5;  // At last vertex
        }
        return INFINITY;  // No valid position found
    }

    // Calculate offset
    double offset = d_before / (d_before + d_after);
    return std::abs((idx_before + offset) / (path_vertices.size() - 1) - 0.5);
}

/**
 * @brief Computes metrics for a path in a uniform grid graph
 *
 * @param path_data Path data structure to be populated with metrics
 * @param graph The uniform grid graph structure
 * @param y Vector of values at the original vertices
 * @param dist_normalization_factor Factor used to normalize distances for kernel weight calculation
 *
 * This function computes several path metrics:
 * - y_path: y-values along the path vertices
 * - x_path: cumulative distances along the path from the start
 * - w_path: kernel weights based on distances from reference vertex
 * - total_weight: total path length
 * - rel_center_offset: measure of how centered the reference vertex is
 */
void compute_grid_path_metrics(
    path_data_t& path_data,
    const uniform_grid_graph_t& graph,
    const std::vector<double>& y,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights
    ) {

    // fprintf(stderr, "\nEntering compute_grid_path_metrics\n");
    // fprintf(stderr, "Path data state:\n");
    // fprintf(stderr, "  ref_vertex: %zu\n", path_data.ref_vertex);
    // //fprintf(stderr, "  ref_vertex_index: %zu\n", path_data.ref_vertex_index);
    // fprintf(stderr, "  vertices size: %zu\n", path_data.vertices.size());
    // fprintf(stderr, "  vertices: ");
    // for (size_t v : path_data.vertices) {
    //     fprintf(stderr, "%zu ", v);
    // }
    // fprintf(stderr, "\n");


    size_t path_n_vertices = path_data.vertices.size();

    path_data.w_path.resize(path_n_vertices);
    path_data.y_path.resize(path_n_vertices);

    // Calculate actual distances along the path using grid graph
    path_data.total_weight = 0;

    for (size_t k = 0; k < path_n_vertices; ++k) {
        path_data.y_path[k] = y[path_data.vertices[k]];
    }
    if (path_data.x_path.size() != path_n_vertices) {
        path_data.x_path.resize(path_n_vertices);
        path_data.x_path[0] = 0;
        path_data.total_weight = 0;  // Need to initialize this
        for (size_t k = 1; k < path_n_vertices; ++k) {
            double weight = graph.get_parent_graph_edge_weight(
                edge_weights,
                path_data.vertices[k-1],
                path_data.vertices[k]
                );

            size_t v1 = path_data.vertices[k-1];
            size_t v2 = path_data.vertices[k];

            auto v_min = std::min(v1, v2);
            auto v_max = std::max(v1, v2);

            auto it = edge_weights.find({v_min, v_max});
            if (it == edge_weights.end()) {
                Rprintf("\nIn compute_grid_path_metrics()\n");
                Rprintf("path_data.vertices[%zu]: %zu\n", k - 1, v1);
                Rprintf("path_data.vertices[%zu]: %zu\n", k, v2);
                print_vect(path_data.vertices, "path_data.vertices");

                REPORT_ERROR("Edge (%zu,%zu) not found in parent graph weights, computing grid distance\n",
                             v_min, v_max);
            }

            path_data.x_path[k] = path_data.x_path[k-1] + weight;
            path_data.total_weight += weight;  // Add weight, not dist
        }
    }

    // Validate initialization
    if (path_data.ref_vertex == INVALID_VERTEX) {
        REPORT_ERROR("path_data ref_vertex not properly initialized - has value INVALID_VERTEX");
    }

    // // Validate bounds
    // if (path_data.ref_vertex_index >= path_data.x_path.size()) {
    //         REPORT_ERROR("Found path_data.ref_vertex_index: %zu is greater than path_data.x_path.size(): %zu\n",
    //                  path_data.ref_vertex_index, path_data.x_path.size());

    // }

    double x_ref = 0; //path_data.x_path[path_data.ref_vertex_index];
    std::vector<double> d_path(path_n_vertices);
    double max_dist = 0.0;

    for (size_t k = 0; k < path_n_vertices; ++k) {
        d_path[k] = std::abs(path_data.x_path[k] - x_ref);
        max_dist = std::max(max_dist, d_path[k]);
    }

    if (max_dist == 0) max_dist = 1;
    max_dist *= dist_normalization_factor;

    for (size_t k = 0; k < path_n_vertices; ++k) {
        d_path[k] /= max_dist;
    }

    kernel_fn(d_path.data(), path_n_vertices, path_data.w_path.data());

    // Normalize weights
    double total_w_path = std::accumulate(
        path_data.w_path.begin(),
        path_data.w_path.end(),
        0.0
    );

    for (size_t k = 0; k < path_n_vertices; ++k) {
        path_data.w_path[k] /= total_w_path;
    }

    path_data.rel_center_offset = calculate_grid_relative_center_offset(
        path_data.ref_vertex,
        path_data.vertices,
        graph
    );
}

/**
 * @brief Constructs a composite path from two paths in a uniform grid graph
 *
 * @param composite_path Output path data structure for the composite path
 * @param path1 First input path
 * @param path2 Second input path
 * @param graph The uniform grid graph structure
 * @param y Vector of values at the original vertices
 * @param dist_normalization_factor Factor used to normalize distances for kernel weights
 *
 * Creates a composite path by:
 * 1. Reversing path2 and concatenating it with path1 (excluding shared reference vertex)
 * 2. Computing path metrics for the composite path
 * The reference vertex becomes the junction point between the two paths.
 */
void construct_grid_composite_path(
    path_data_t& composite_path,
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    const uniform_grid_graph_t& graph,
    const std::vector<double>& y,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights
    ) {

    // fprintf(stderr, "\nEntering construct_grid_composite_path\n");
    // fprintf(stderr, "Path1 size: %zu, Path2 size: %zu\n",
    //         path1.vertices.size(), path2.vertices.size());
    // fprintf(stderr, "Reference vertex: %zu\n", composite_path.ref_vertex);

#if 0
    // Validate input paths
    if (path1.vertices.empty() || path2.vertices.empty()) {
        REPORT_ERROR("Empty path provided for composite path construction");
    }
    if (path1.vertices.size() < 2 || path2.vertices.size() < 2) {
        REPORT_ERROR("Paths must have at least 2 vertices for composite path construction");
    }
#endif

    // Reserve space for the composite path
    size_t combined_length = path1.vertices.size() + path2.vertices.size();
    if (combined_length < path1.vertices.size()) { // Check for overflow
        REPORT_ERROR("Path length overflow detected");
    }
    combined_length -= 1; // Now safe to subtract

    composite_path.vertices.reserve(combined_length);

    // Add reversed path2
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path2.vertices.rbegin(),
        path2.vertices.rend()
    );

    // Add path1 (excluding the first vertex if it is the same as in path2)
    size_t offset = 0;
    if (path1.vertices[0] == path2.vertices[0]) {
        offset = 1;
    }
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path1.vertices.begin() + offset,
        path1.vertices.end()
        );

    // Calculate total weight using grid graph distances
    composite_path.total_weight = path1.total_weight + path2.total_weight;

    double ref_vertex_offset = path1.dist_to_ref_vertex[1] / (path1.dist_to_ref_vertex[1] + path2.dist_to_ref_vertex[1]);
    // The centrality measure is calculated as:
    // abs((idx_before + ref_vertex_offset)/(path_size - 1) - 0.5)
    double idx_before =  path1.dist_to_ref_vertex.size() - 2;
    double path_size_minus_one = composite_path.vertices.size() - 1.0;
    composite_path.rel_center_offset = fabs((idx_before + ref_vertex_offset) / path_size_minus_one - 0.5);

#if 1
    // Compute all path metrics
    Rprintf("\nBefore compute_grid_path_metrics() in construct_grid_composite_path()\n");

    // After constructing composite path
    // fprintf(stderr, "Composite path vertices: ");
    // for (size_t v : composite_path.vertices) {
    //     fprintf(stderr, "%zu ", v);
    // }
    // //fprintf(stderr, "\nRef vertex index: %zu\n", composite_path.ref_vertex_index);


    Rprintf("\nIn construct_grid_composite_path()\n");
    auto print_ref_vertex_path_t = [](const ref_vertex_path_t& path, const std::string& name = "") {
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
        Rprintf("\nPath %s:\n", name.c_str());  // Fixed: using name instead of names
        Rprintf("  Vertices:       %s\n", vertices_str.str().c_str());
        Rprintf("  Target Vertex:  %zu\n", path.target_vertex);
        Rprintf("  Total Weight:   %.4f\n", path.total_weight);
    };

    print_ref_vertex_path_t(path1, "path1");
    print_ref_vertex_path_t(path2, "path2");

    print_vect(composite_path.vertices,"composite_path.vertices");

#endif

    compute_grid_path_metrics(
        composite_path,
        graph,
        y,
        dist_normalization_factor,
        edge_weights
    );
}

/**
 * @brief Reconstructs a path in a grid graph keeping only original vertices
 *
 * @param shortest_paths_map Map of vertices to their distance and predecessor from BFS
 * @param vertex Starting vertex for path reconstruction
 * @param n_original_vertices Number of original vertices in the graph
 * @return std::vector<int> Path containing only original vertices
 *
 * Reconstructs a path from the BFS results, filtering out grid vertices and
 * keeping only the original vertices that contain data values.
 */
std::vector<size_t> reconstruct_path_with_grid(
    const std::unordered_map<size_t, std::pair<double, size_t>>& shortest_paths_map,
    size_t vertex,
    size_t n_original_vertices
) {
    std::vector<size_t> path;

    // Only add to path if it's an original vertex
    if (is_original_vertex(vertex, n_original_vertices)) {
        path.push_back(vertex);
    }

    while (true) {
        auto it = shortest_paths_map.find(vertex);
        if (it == shortest_paths_map.end()) break;

        vertex = it->second.second;  // Get prev vertex
        if (vertex == INVALID_VERTEX) break;

        // Only add to path if it's an original vertex
        if (is_original_vertex(vertex, n_original_vertices)) {
            path.push_back(vertex);
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Finds shortest paths within a given radius in a uniform grid graph
 *
 * @param graph The uniform grid graph structure
 * @param start Starting vertex for the search
 * @param radius Maximum distance to search
 * @return std::unordered_map<size_t, std::pair<double, size_t>> Map of vertices to their (distance, predecessor)
 *
 * Implements Dijkstra's algorithm to find shortest paths in the grid graph,
 * stopping when path distances exceed the given radius. Keeps track of all vertices
 * (both original and grid) for accurate path reconstruction.
 */
std::unordered_map<size_t, std::pair<double, size_t>> find_grid_shortest_paths_within_radius(
    const uniform_grid_graph_t& graph,
    size_t start,
    double radius
) {
    std::unordered_map<size_t, std::pair<double, size_t>> result;           // vertex -> (distance from start_vertex, predecessor along the shortest path from start_vertex to vertex)
    std::vector<double> dist(graph.adjacency_list.size(), INFINITY);
    std::vector<size_t> prev(graph.adjacency_list.size(), INVALID_VERTEX);
    dist[start] = 0;

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        d = -d;
        pq.pop();

        if (d > radius) break;
        if (d > dist[u]) continue;

        // Store ALL vertices we visit, in order to reconstruct the path - we can check in the reconstructed path which vertices are original (their indices are less than  graph.n_original_vertices)
        result.emplace(u, std::pair<double, size_t>(d, prev[u]));

        for (const auto& [v, w] : graph.adjacency_list[u]) {
            if (dist[u] + w < dist[v] && dist[u] + w <= radius) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});
            }
        }
    }

    return result;
}

/**
 * @brief Reconstructs paths from current search state in a grid graph
 *
 * @param vertex_info Map containing vertex information (distance, predecessor, path length)
 * @param start Starting vertex (reference vertex)
 * @param n_original_vertices Number of original vertices in the graph
 * @return std::vector<ref_vertex_path_t> Vector of reconstructed paths
 *
 * Reconstructs paths from the current search state, filtering to keep only
 * original vertices. Orders paths by distance and ensures no vertex is used
 * in multiple paths.
 */
std::vector<ref_vertex_path_t> reconstruct_current_grid_paths(
    const std::unordered_map<size_t, std::tuple<double, size_t, size_t>>& vertex_info,
    size_t start,
    size_t n_original_vertices
) {
    std::vector<ref_vertex_path_t> paths;
    std::unordered_set<size_t> used_vertices;
    used_vertices.insert(start);

    // Create ordered vertex info for processing
    std::vector<std::pair<size_t, double>> vertices;
    for (const auto& [vertex, info] : vertex_info) {
        if (vertex != start && is_original_vertex(vertex, n_original_vertices)) {
            vertices.emplace_back(vertex, std::get<0>(info));  // vertex and its distance
        }
    }

    // Sort by distance
    std::sort(vertices.begin(), vertices.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Reconstruct paths from farthest vertices
    for (const auto& [vertex, distance] : vertices) {
        if (used_vertices.count(vertex) > 0) continue;

        ref_vertex_path_t new_path;
        new_path.target_vertex = start;
        new_path.total_weight = distance;

        // Reconstruct path, keeping only original vertices
        size_t curr = vertex;
        while (curr != INVALID_VERTEX) {
            if (is_original_vertex(curr, n_original_vertices)) {
                new_path.vertices.push_back(curr);
                used_vertices.insert(curr);
                new_path.dist_to_ref_vertex.push_back(std::get<0>(vertex_info.at(curr)));
            }
            curr = std::get<1>(vertex_info.at(curr));  // get predecessor
        }

        std::reverse(new_path.vertices.begin(), new_path.vertices.end());
        std::reverse(new_path.dist_to_ref_vertex.begin(), new_path.dist_to_ref_vertex.end());

        paths.push_back(std::move(new_path));
    }

    return paths;
}

/**
 * @brief Determines if two paths explore different directions by comparing their
 *        original vertices after the reference vertex
 *
 * This function examines whether paths take meaningfully different routes through
 * the graph by looking at their original vertices. It ignores grid vertices since
 * they don't carry data values and are just intermediary points for distance
 * calculations.
 *
 * @param path1 First path to compare
 * @param path2 Second path to compare
 * @param diff_threshold Maximum number of original vertices to compare after ref_vertex
 * @return true if paths explore different directions (no shared original vertices)
 */
bool paths_explore_different_directions(
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    size_t diff_threshold
) {
    // Create sets to store the original vertices we want to compare
    std::set<size_t> set1, set2;

    // Skip index 0 as it's the reference vertex shared by both paths
    // Insert up to diff_threshold vertices from path1
    size_t count1 = 0;
    for (size_t i = 1; i < path1.vertices.size() && count1 < diff_threshold; i++) {
        set1.insert(path1.vertices[i]);
        count1++;
    }

    // Insert up to diff_threshold vertices from path2
    size_t count2 = 0;
    for (size_t i = 1; i < path2.vertices.size() && count2 < diff_threshold; i++) {
        set2.insert(path2.vertices[i]);
        count2++;
    }

    // Find any vertices that appear in both sets
    std::vector<size_t> intersection;
    std::set_intersection(
        set1.begin(), set1.end(),
        set2.begin(), set2.end(),
        std::back_inserter(intersection)
    );

    // Paths explore different directions if they share no vertices
    return intersection.empty();
}

/**
 * @brief Finds paths meeting minimum size requirement in a uniform grid graph
 *
 * @param graph The uniform grid graph structure
 * @param y Vector of values at original vertices
 * @param ref_vertex Reference vertex for path search
 * @param min_path_size Minimum number of vertices required in paths
 * @param diff_threshold Threshold for determining if paths explore different directions
 * @param dist_normalization_factor Factor for normalizing distances in kernel calculations
 * @return std::vector<path_data_t> Vector of paths meeting requirements
 *
 * Implements a modified shortest path search that tracks path lengths and tries to find:
 * 1. A single path meeting the size requirement
 * 2. If no single path exists, looks for composite paths meeting the requirement
 * Only counts original vertices (non-grid vertices) when checking path sizes.
 */
std::vector<path_data_t> find_grid_paths_meeting_size_requirement(
    const uniform_grid_graph_t& graph,
    const std::vector<double>& y,
    size_t ref_vertex,
    size_t min_path_size,
    size_t diff_threshold,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights
) {
    std::vector<path_data_t> result;
    std::vector<double> dist(graph.adjacency_list.size(), INFINITY);
    std::vector<size_t> prev(graph.adjacency_list.size(), INVALID_VERTEX);
    std::vector<size_t> path_length(graph.adjacency_list.size(), 1);
    dist[ref_vertex] = 0;

    struct queue_entry_t {
        double distance;
        size_t vertex;
        size_t path_length;

        bool operator>(const queue_entry_t& other) const {
            return distance > other.distance;
        }
    };

    std::priority_queue<queue_entry_t, std::vector<queue_entry_t>, std::greater<>> pq;
    pq.push({0, ref_vertex, 1});

    // Store complete vertex information for path reconstruction
    std::unordered_map<size_t, std::tuple<double, size_t, size_t>> vertex_info;  // vertex -> (distance, prev, path_length)

    while (!pq.empty()) {
        auto [d, u, curr_path_len] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        // Store information for ALL vertices we visit
        vertex_info[u] = std::make_tuple(d, prev[u], curr_path_len);

        // First try to find a single path meeting the requirement
        if (is_original_vertex(u, graph.n_original_vertices) && curr_path_len >= min_path_size) {
            path_data_t path_data;
            path_data.ref_vertex = ref_vertex;

            // Reconstruct path keeping only original vertices
            std::vector<size_t> full_path;
            size_t curr = u;
            while (curr != INVALID_VERTEX) {
                if (is_original_vertex(curr, graph.n_original_vertices)) {
                    full_path.push_back(curr);
                }
                curr = std::get<1>(vertex_info[curr]);
            }
            std::reverse(full_path.begin(), full_path.end());

            path_data.vertices = full_path;

            path_data.vertices = full_path;

            // Find position of reference vertex in the path
            auto it = std::find(path_data.vertices.begin(),
                                path_data.vertices.end(),
                                path_data.ref_vertex);
            if (it == path_data.vertices.end()) {
                REPORT_ERROR("Reference vertex not found in path vertices");
            }
            //path_data.ref_vertex_index = it - path_data.vertices.begin();

            // debugging
            // Rprintf("\nBefore compute_grid_path_metrics() in find_grid_paths_meeting_size_requirement()\n"); // debugging
            compute_grid_path_metrics(path_data,
                                      graph,
                                      y,
                                      dist_normalization_factor,
                                      edge_weights);
            result.push_back(std::move(path_data));
            return result;
        }

        // Explore neighbors
        for (const auto& [v, w] : graph.adjacency_list[u]) {
            double new_dist = dist[u] + w;
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;

                size_t new_path_len = curr_path_len;
                if (is_original_vertex(v, graph.n_original_vertices)) {
                    new_path_len++;
                }

                path_length[v] = new_path_len;
                pq.push({new_dist, v, new_path_len});
            }
        }

        // Check for composite paths after exploring enough vertices
        if (!vertex_info.empty() && curr_path_len > 1) {
            auto paths = reconstruct_current_grid_paths(vertex_info, ref_vertex, graph.n_original_vertices);

            // Look for pairs of paths that can form valid composite paths
            for (size_t i = 0; i < paths.size(); ++i) {
                for (size_t j = i + 1; j < paths.size(); ++j) {
                    // Check combined Rf_length(remember these paths now contain only original vertices)
                    size_t combined_length = paths[i].vertices.size() + paths[j].vertices.size() - 1;
                    if (combined_length < min_path_size) continue;

                    // Check if paths explore different directions
                    if (!paths_explore_different_directions(paths[i], paths[j], diff_threshold)) continue;

                    // Create composite path
                    path_data_t composite_path;
                    composite_path.ref_vertex = ref_vertex;
                    construct_grid_composite_path(
                        composite_path,
                        paths[i],
                        paths[j],
                        graph,
                        y,
                        dist_normalization_factor,
                        edge_weights
                    );

                    result.push_back(std::move(composite_path));
                    return result;
                }
            }
        }
    }

    return result;
}

/**
 * @brief Precomputes shortest path information from a reference vertex within a maximum bandwidth radius
 *        for efficient path-based computations.
 *
 * This function computes shortest paths using Dijkstra's algorithm from the reference vertex
 * to all reachable vertices within the specified bandwidth radius. The paths are stored in a
 * format optimized for later use, including:
 * - Distances from the reference vertex to all reachable vertices
 * - Predecessor information for path reconstruction
 * - A list of original vertices (excluding grid vertices) sorted by distance in descending order
 *
 * @param uniform_grid_graph The uniform grid graph containing both original vertices and
 *                          additional grid vertices used for path refinement
 * @param ref_vertex The source vertex from which to compute all shortest paths
 * @param max_bandwidth The maximum path Rf_length(radius) to consider when computing paths.
 *                     Paths longer than this value are not computed or stored
 * @return reachability_map_t A structure containing:
 *         - distances: Map of vertex ID to its distance from ref_vertex
 *         - predecessors: Map of vertex ID to its predecessor in the shortest path
 *         - sorted_vertices: Vector of original vertices (excluding grid vertices)
 *                          sorted by distance from ref_vertex in descending order
 *         - ref_vertex: The source vertex used for path computation
 */
reachability_map_t precompute_max_bandwidth_paths(
    const uniform_grid_graph_t& uniform_grid_graph,
    size_t ref_vertex,
    double max_bandwidth) {

    reachability_map_t result;
    result.ref_vertex = ref_vertex;

    // Compute shortest paths once for maximum bandwidth
    auto shortest_paths_map = find_grid_shortest_paths_within_radius(
        uniform_grid_graph,
        ref_vertex,
        max_bandwidth
    );

    // Store distances, predecessors and sorted vertices (in descending order by distance)
    for (const auto& [vertex, dist_parent_info] : shortest_paths_map) {
        result.distances[vertex] = dist_parent_info.first;
        result.predecessors[vertex] = dist_parent_info.second;
        if (vertex != ref_vertex &&
            is_original_vertex(vertex, uniform_grid_graph.n_original_vertices)) {
            result.sorted_vertices.push_back({vertex, dist_parent_info.first});
        }
    }

    std::sort(result.sorted_vertices.begin(), result.sorted_vertices.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    return result;
}

/**
 * @brief Processes paths to generate path_data_t objects
 *
 * @param paths Vector of ref_vertex_path_t objects to process
 * @param uniform_grid_graph The uniform grid graph
 * @param y Vector of y-values associated with vertices
 * @param min_path_size Minimum required path size
 * @param diff_threshold Threshold for determining different directions
 * @param dist_normalization_factor Factor for normalizing distances
 * @return Vector of processed path_data_t objects
 */
std::vector<path_data_t> process_paths(
    const std::vector<ref_vertex_path_t>& paths,
    const uniform_grid_graph_t& uniform_grid_graph,
    const std::vector<double>& y,
    size_t min_path_size,
    size_t diff_threshold,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights) {

    // fprintf(stderr, "\nEntering process_paths()\n");
    // fprintf(stderr, "Number of input paths: %zu\n", paths.size());

    if (paths.empty()) {
        REPORT_ERROR("paths vector is empty");
        return {};  // Return empty result for empty input
    }

    std::vector<path_data_t> result;

    if (paths.size() == 1 && paths[0].vertices.size() >= min_path_size) {

        // fprintf(stderr, "Processing single path case\n");
        // fprintf(stderr, "Path vertices: ");
        // for (size_t v : paths[0].vertices) {
        //     fprintf(stderr, "%zu ", v);
        // }
        // fprintf(stderr, "\nTarget vertex: %zu\n", paths[0].target_vertex);

        // Single path case
        path_data_t path_data;
        path_data.vertices   = std::move(paths[0].vertices);
        path_data.x_path     = std::move(paths[0].dist_to_ref_vertex);
        path_data.ref_vertex = paths[0].target_vertex;
        path_data.rel_center_offset = 0.5; // the ref vertex is at the boundary of the path

        // Find position of target vertex in the path

        // WARNING !!!! <<---

        // If path_data.ref_vertex is a grid vertex that is not an original
        // vertex, then it will not be a part of path_data.vertices and so it
        // does not make sense to talk about its index within path_data.vertices


        path_data.total_weight = paths[0].total_weight;

        compute_grid_path_metrics(
            path_data,
            uniform_grid_graph,
            y,
            dist_normalization_factor,
            edge_weights
        );

        result.push_back(std::move(path_data));
    }
    else if (paths.size() > 1) {
        // Look for composite paths
        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = i + 1; j < paths.size(); ++j) {
                size_t combined_size = paths[i].vertices.size() +
                                  paths[j].vertices.size() - 1;

                if (combined_size < min_path_size) continue;

                // Check if paths explore different directions
                if (!paths_explore_different_directions(paths[i], paths[j], diff_threshold)) continue;

                path_data_t composite_path;
                composite_path.ref_vertex = paths[i].target_vertex;

                construct_grid_composite_path( // this is where rel_center_offset would need to be defined/derived
                    composite_path,
                    paths[i],
                    paths[j],
                    uniform_grid_graph,
                    y,
                    dist_normalization_factor,
                    edge_weights
                );

                result.push_back(std::move(composite_path));
            }
        }
    }

    return result;
}

/**
 * @brief Generates path data for a specific bandwidth using precomputed information
 */
std::vector<path_data_t> generate_paths_for_bandwidth(
    const uniform_grid_graph_t& uniform_grid_graph,
    const reachability_map_t& reachability_map,
    const std::vector<double>& y,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    const edge_weights_t& edge_weights) {

    // fprintf(stderr, "\nEntering generate_paths_for_bandwidth()\n");
    // fprintf(stderr, "Reference vertex: %zu\n", reachability_map.ref_vertex);
    // fprintf(stderr, "Number of sorted vertices: %zu\n", reachability_map.sorted_vertices.size());
    //
    // fprintf(stderr, "\nbandwidth: %f\n\n", bandwidth);

    //print_vect(reachability_map.sorted_vertices,"reachability_map.sorted_vertices");
    // print_umap(reachability_map.distances, "reachability_map.distances");

    std::vector<ref_vertex_path_t> paths;
    std::unordered_set<size_t> used_vertices;

    // Use only vertices within current bandwidth
    for (const auto& vertex_dist_info : reachability_map.sorted_vertices) {

        if (vertex_dist_info.distance > bandwidth) continue;
        if (used_vertices.count(vertex_dist_info.vertex) > 0) continue;

        ref_vertex_path_t new_path;
        new_path.target_vertex = reachability_map.ref_vertex;
        new_path.total_weight = vertex_dist_info.distance;

        // Reconstruct path using precomputed predecessors
        std::vector<size_t> full_path;
        size_t curr = vertex_dist_info.vertex;

        // fprintf(stderr, "\nReconstructing path from vertex %zu\n", curr);

        while (curr != INVALID_VERTEX) {
            // fprintf(stderr, "Processing vertex %zu (is_original: %d)\n",
            //         curr, (int)is_original_vertex(curr, y));
            if (is_original_vertex(curr, y)) {
                new_path.vertices.push_back(curr);
                used_vertices.insert(curr);

                auto dist_it = reachability_map.distances.find(curr);
                if (dist_it != reachability_map.distances.end()) {
                    new_path.dist_to_ref_vertex.push_back(dist_it->second);
                }
            }

            auto pred_it = reachability_map.predecessors.find(curr);
            if (pred_it == reachability_map.predecessors.end()) break;
            curr = pred_it->second;
        }

        std::reverse(new_path.vertices.begin(), new_path.vertices.end());
        std::reverse(new_path.dist_to_ref_vertex.begin(), new_path.dist_to_ref_vertex.end());


        print_vect(new_path.vertices, "new_path.vertices");
        print_vect(new_path.dist_to_ref_vertex, "new_path.dist_to_ref_vertex");


        if (new_path.vertices.size() >= 2) {
            #if 0
            // relative center offset needs to be calculated for the final paths constructed from these outgoing from the ref vertex paths
            new_path.center_offset = calculate_grid_relative_center_offset(
                reachability_map.ref_vertex,
                new_path.vertices,
                uniform_grid_graph
                );
            #endif
            paths.push_back(std::move(new_path));
        }
    }

    return process_paths(paths,
                         uniform_grid_graph,
                         y,
                         min_path_size,
                         diff_threshold,
                         dist_normalization_factor,
                         edge_weights);
}


/**
 * @brief Modified version of ugg_get_path_data that uses precomputed paths
 *
 * @param uniform_grid_graph The uniform grid graph
 * @param y Vector of y-values associated with vertices
 * @param ref_vertex Reference vertex for path computation
 * @param bandwidth Current bandwidth to consider
 * @param reachability_map Precomputed path information for maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances
 * @param min_path_size Minimum required path size
 * @param diff_threshold Threshold for determining different directions
 * @param kernel_type Type of kernel to use for weighting
 * @return Vector of processed path_data_t objects
 */
std::vector<path_data_t> ugg_get_path_data_efficient(
    const uniform_grid_graph_t& uniform_grid_graph,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    const reachability_map_t& reachability_map,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    const edge_weights_t& edge_weights
    ) {

    initialize_kernel(kernel_type, 1.0);

    std::vector<path_data_t> result = generate_paths_for_bandwidth(
        uniform_grid_graph,
        reachability_map,
        y,
        bandwidth,
        dist_normalization_factor,
        min_path_size,
        diff_threshold,
        edge_weights
        );

    // fprintf(stderr, "In ugg_get_path_data_efficient()\nNumber path_data_t objects found: %zu\n", result.size());

    if (result.empty()) {
        return find_grid_paths_meeting_size_requirement(
            uniform_grid_graph,
            y,
            ref_vertex,
            min_path_size,
            diff_threshold,
            dist_normalization_factor,
            edge_weights
        );
    }

    return result;
}


/**
 * @brief Computes path data for paths centered around a reference vertex in a uniform grid graph
 *
 * This function analyzes a uniform grid graph to find and process paths that include
 * a specified reference vertex. It computes various metrics and properties for these paths,
 * applying kernel-based weightings and considering distance constraints.
 *
 * @param uniform_grid_graph The uniform grid graph to analyze
 * @param y Vector of y-values associated with vertices
 * @param ref_vertex Index of the reference vertex to center paths around
 * @param bandwidth Maximum allowed distance for path construction
 * @param dist_normalization_factor Factor for normalizing distance calculations
 * @param min_path_size Minimum required number of vertices in a valid path
 * @param diff_threshold Threshold for determining if paths explore different directions
 * @param kernel_type Type of kernel to use for weighting calculations
 * @param verbose Whether to print detailed progress information
 *
 * @return Vector of path_data_t structures containing analyzed path information
 *
 * @throws Reports Rf_error if y.size() doesn't match number of original vertices
 *
 * @details The function performs the following steps:
 * 1. Verifies input data consistency
 * 2. Finds shortest paths within specified bandwidth
 * 3. Reconstructs paths keeping only original vertices
 * 4. Processes single paths meeting size requirements
 * 5. Combines compatible paths into composite paths if needed
 * 6. Falls back to unrestricted path search if no suitable paths found
 */
std::vector<path_data_t> ugg_get_path_data(
    const uniform_grid_graph_t& uniform_grid_graph,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    const edge_weights_t& edge_weights,
    bool verbose
    ) {

    // We can verify at the start that y.size() matches our expectation
    if (y.size() != uniform_grid_graph.n_original_vertices) {
        REPORT_ERROR("Mismatch between y size and number of original vertices");
    }

    initialize_kernel(kernel_type, 1.0);

    // First try with bandwidth constraint
    auto shortest_paths_map = find_grid_shortest_paths_within_radius(
        uniform_grid_graph,
        ref_vertex,
        bandwidth
        );

    // Reconstruct paths similarly to get_path_data() function, but only keeping original vertices within the uniform_grid_graph
    std::vector<ref_vertex_path_t> paths;
    std::unordered_set<size_t> used_vertices;

    // Create ordered vertex info for processing
    std::vector<vertex_info_t> vertex_info;
    for (const auto& [vertex, info] : shortest_paths_map) {
        if (vertex != ref_vertex &&
            is_original_vertex(vertex, uniform_grid_graph.n_original_vertices)) {
            vertex_info.push_back({vertex, info.first});
        }
    }

    // Sort by distance
    std::sort(vertex_info.begin(), vertex_info.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    // Reconstruct paths
    for (const auto& info : vertex_info) {
        if (used_vertices.count(info.vertex) > 0) continue;

        ref_vertex_path_t new_path;
        new_path.target_vertex = ref_vertex;
        new_path.total_weight = info.distance;

        // Reconstruct path, keeping only original vertices
        std::vector<size_t> full_path;
        size_t curr = info.vertex;
        while (curr != INVALID_VERTEX) {
            if (is_original_vertex(curr, y)) {
                new_path.vertices.push_back(curr);
                used_vertices.insert(curr);
            }

            auto it = shortest_paths_map.find(curr);
            if (it == shortest_paths_map.end()) break;

            new_path.dist_to_ref_vertex.push_back(it->second.first);
            curr = it->second.second;
        }

        // Reverse the path to get correct order (from target to end)
        std::reverse(new_path.vertices.begin(), new_path.vertices.end());
        std::reverse(new_path.dist_to_ref_vertex.begin(), new_path.dist_to_ref_vertex.end());

        if (new_path.vertices.size() >= 2) {
            new_path.center_offset = calculate_grid_relative_center_offset(
                ref_vertex,
                new_path.vertices,
                uniform_grid_graph
            );
            paths.push_back(std::move(new_path));
        }
    }

    if (verbose) {
        size_t n_paths = paths.size();
        Rprintf("n_paths: %zu\n", n_paths);
        print_vect_paths(paths, "reconstructed paths from shortest_paths_map");
    }

    // Process paths similarly to original version
    std::vector<path_data_t> result;

    if (paths.size() == 1 && paths[0].vertices.size() >= min_path_size) {
        path_data_t path_data;
        path_data.ref_vertex = ref_vertex;
        path_data.vertices = std::move(paths[0].vertices);
        path_data.x_path   = std::move(paths[0].dist_to_ref_vertex);
        path_data.total_weight = paths[0].total_weight;

        // Compute metrics using uniform grid graph
        Rprintf("\nBefore compute_grid_path_metrics() in ugg_get_path_data()\n"); // debugging
        compute_grid_path_metrics(
            path_data,
            uniform_grid_graph,
            y,
            dist_normalization_factor,
            edge_weights
        );

        result.push_back(std::move(path_data));
    } else if (paths.size() > 1) {
        // Look for composite paths
        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = i + 1; j < paths.size(); ++j) {
                size_t combined_size = paths[i].vertices.size() +
                                  paths[j].vertices.size() - 1;

                if (combined_size < min_path_size) continue;

                // Check if paths explore different directions
                if (!paths_explore_different_directions(paths[i], paths[j], diff_threshold)) continue;

                path_data_t composite_path;
                composite_path.ref_vertex = ref_vertex;

                // Construct composite path keeping only original vertices
                construct_grid_composite_path(
                    composite_path,
                    paths[i],
                    paths[j],
                    uniform_grid_graph,
                    y,
                    dist_normalization_factor,
                    edge_weights
                );

                result.push_back(std::move(composite_path));
            }
        }
    }

    // If no suitable paths found, expand search without bandwidth constraint
    if (result.empty()) {
        return find_grid_paths_meeting_size_requirement(
            uniform_grid_graph,
            y,
            ref_vertex,
            min_path_size,
            diff_threshold,
            dist_normalization_factor,
            edge_weights
        );
    }

    return result;
}

/**
 * @brief R interface for computing path data centered around a reference vertex
 *        with associated distances, weights, and values along the paths
 *
 * @param adj_list_s SEXP List of integer vectors representing graph adjacency (1-based indices)
 *                   Must be a valid graph representation with matching dimensions
 * @param weight_list_s SEXP List of numeric vectors containing edge weights
 *                     Must match adj_list_s structure exactly
 * @param y_s SEXP Numeric vector containing values associated with each vertex
 * @param ref_vertex_s SEXP Integer scalar specifying reference vertex (1-based index)
 *                    Must be within valid range for the graph
 * @param bandwidth_s SEXP Numeric scalar specifying maximum path distance from reference vertex
 * @param dist_normalization_factor_s SEXP Numeric scalar for distance normalization in kernel computations
 * @param min_path_size_s SEXP Integer scalar specifying minimum required path length
 * @param kernel_type_s SEXP Integer scalar specifying kernel function type (1-10)
 *
 * @return SEXP (list) where each element is a list with named components:
 *         - "vertices": integer vector of path vertices (1-based indices)
 *         - "ref_vertex": integer scalar of reference vertex (1-based index)
 *         - "rel_center_offset": numeric scalar measuring reference vertex position (0 = center, 0.5 = endpoint)
 *         - "total_weight": numeric scalar of path length
 *         - "x_path": numeric vector of cumulative distances along path from ref_vertex
 *         - "w_path": numeric vector of kernel weights for each vertex
 *         - "y_path": numeric vector of y-values for path vertices
 *
 * @note All vertex indices in input must be positive integers
 * @note All edge weights must be non-negative
 * @note Kernel type must be between 1 and 10 (see MSR2_KERNELS.H for definitions)
 */
SEXP S_ugg_get_path_data(
    SEXP adj_list_s,
    SEXP weight_list_s,
    SEXP grid_size_s,
    SEXP y_s,
    SEXP ref_vertex_s,
    SEXP bandwidth_s,
    SEXP dist_normalization_factor_s,
    SEXP min_path_size_s,
    SEXP diff_threshold_s,
    SEXP kernel_type_s,
    SEXP verbose_s
    ) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(adj_list_s);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_s);

    size_t grid_size   = (size_t) Rf_asInteger(grid_size_s);
    size_t ref_vertex  = (size_t) Rf_asInteger(ref_vertex_s);

    std::vector<double> y(REAL(y_s), REAL(y_s) + LENGTH(y_s));

    double bandwidth                 = Rf_asReal(bandwidth_s);
    double dist_normalization_factor = Rf_asReal(dist_normalization_factor_s);

    size_t min_path_size  = (size_t) Rf_asInteger(min_path_size_s);
    size_t diff_threshold = (size_t) Rf_asInteger(diff_threshold_s);
    size_t kernel_type    = (size_t) Rf_asInteger(kernel_type_s);

    bool verbose = (Rf_asLogical(verbose_s) == TRUE);

    double snap_tolerance = 0.1;
    uniform_grid_graph_t uniform_grid_graph = create_uniform_grid_graph(
        adj_list,
        weight_list,
        grid_size,
        ref_vertex,
        snap_tolerance);

    edge_weights_t edge_weights = precompute_edge_weights(adj_list, weight_list);

    std::vector<path_data_t> paths = ugg_get_path_data(
        uniform_grid_graph,
        y,
        ref_vertex,
        bandwidth,
        dist_normalization_factor,
        min_path_size,
        diff_threshold,
        kernel_type,
        edge_weights,
        verbose
        );

    size_t n_paths = paths.size();

    // Convert results to R list
    const size_t RESULT_LIST_SIZE = n_paths + 3;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, RESULT_LIST_SIZE));

    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, RESULT_LIST_SIZE));
    for(size_t i = 0; i < n_paths; i++) {
        std::string path_id = "path" + std::to_string(i);
        SET_STRING_ELT(result_names, i, Rf_mkChar(path_id.c_str()));
    }

    SET_STRING_ELT(result_names, n_paths + 0, Rf_mkChar("ugg_adj_list"));
    SET_STRING_ELT(result_names, n_paths + 1, Rf_mkChar("ugg_weight_list"));
    SET_STRING_ELT(result_names, n_paths + 2, Rf_mkChar("ugg_grid_vertices"));

    Rf_setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(1);  // result_names

    // Define component names matching path_data_t structure
    const std::vector<std::string> path_comps_names = {
        "vertices", "ref_vertex", "rel_center_offset", "total_weight",
        "x_path", "w_path", "y_path"
    };

    // Create SEXP for names once
    SEXP names = PROTECT(Rf_allocVector(STRSXP, path_comps_names.size()));
    for(size_t i = 0; i < path_comps_names.size(); i++) {
        SET_STRING_ELT(names, i, Rf_mkChar(path_comps_names[i].c_str()));
    }

    for(size_t i = 0; i < n_paths; i++) {
        // Create list for current path (7 elements matching path_data_t)
        SEXP path = PROTECT(Rf_allocVector(VECSXP, 7));

        // Convert vertices to R (adding 1 for 1-based indexing)
        SEXP vertices = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)paths[i].vertices.size()));
        for(size_t j = 0; j < paths[i].vertices.size(); j++) {
            INTEGER(vertices)[j] = paths[i].vertices[j] + 1;
        }

        // Create scalar elements
        SEXP ref_vertex_r = PROTECT(Rf_ScalarInteger(paths[i].ref_vertex + 1));
        SEXP rel_center_offset = PROTECT(Rf_ScalarReal(paths[i].rel_center_offset));
        SEXP total_weight = PROTECT(Rf_ScalarReal(paths[i].total_weight));

        // Convert vector elements
        SEXP x_path = PROTECT(Rf_allocVector(REALSXP, paths[i].x_path.size()));
        SEXP w_path = PROTECT(Rf_allocVector(REALSXP, paths[i].w_path.size()));
        SEXP y_path = PROTECT(Rf_allocVector(REALSXP, paths[i].y_path.size()));

        // Copy vector data
        std::copy(paths[i].x_path.begin(), paths[i].x_path.end(), REAL(x_path));
        std::copy(paths[i].w_path.begin(), paths[i].w_path.end(), REAL(w_path));
        std::copy(paths[i].y_path.begin(), paths[i].y_path.end(), REAL(y_path));

        // Set list elements
        SET_VECTOR_ELT(path, 0, vertices);
        SET_VECTOR_ELT(path, 1, ref_vertex_r);
        SET_VECTOR_ELT(path, 2, rel_center_offset);
        SET_VECTOR_ELT(path, 3, total_weight);
        SET_VECTOR_ELT(path, 4, x_path);
        SET_VECTOR_ELT(path, 5, w_path);
        SET_VECTOR_ELT(path, 6, y_path);

        // Set names for the current path list
        Rf_setAttrib(path, R_NamesSymbol, names);

        // Add path to result
        SET_VECTOR_ELT(result, i, path);

        // Unprotect all elements created in this iteration
        UNPROTECT(8);  // vertices, ref_vertex_r, rel_center_offset, total_weight, x_path, w_path, y_path, path
    }
    UNPROTECT(1);  // names

    // Extract adjacency and weight lists from result
    size_t n_total_vertices = uniform_grid_graph.adjacency_list.size();
    SEXP r_adj_list = PROTECT(Rf_allocVector(VECSXP, n_total_vertices));
    SEXP r_weight_list = PROTECT(Rf_allocVector(VECSXP, n_total_vertices));

    // Convert the set-based representation back to R lists
    for (size_t i = 0; i < n_total_vertices; ++i) {
        const auto& neighbors = uniform_grid_graph.adjacency_list[i];

        // Create vectors for this vertex's adjacency list and weights
        SEXP r_adj = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)neighbors.size()));
        SEXP r_weights = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)neighbors.size()));

        // Fill the vectors
        size_t idx = 0;
        for (const auto& [neighbor, weight] : neighbors) {
            // Convert to 1-based indices for R
            INTEGER(r_adj)[idx] = neighbor + 1;
            REAL(r_weights)[idx] = weight;
            ++idx;
        }

        SET_VECTOR_ELT(r_adj_list, i, r_adj);
        SET_VECTOR_ELT(r_weight_list, i, r_weights);
        UNPROTECT(2); // for r_adj and r_weights
    }
    SET_VECTOR_ELT(result, n_paths + 0, r_adj_list);
    SET_VECTOR_ELT(result, n_paths + 1, r_weight_list);
    UNPROTECT(2); // for r_adj_list and r_weights_list

    // Create grid vertices vector (1-based indices)
    SEXP r_grid_vertices = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)uniform_grid_graph.grid_vertices.size()));

    size_t counter = 0;
    for (const auto& i : uniform_grid_graph.grid_vertices) {
        // Convert to 1-based indices for R
        INTEGER(r_grid_vertices)[counter++] = i + 1;
    }
    SET_VECTOR_ELT(result, n_paths + 2, r_grid_vertices);
    UNPROTECT(1); // r_grid_vertices

    UNPROTECT(1); // result

    return result;
}
