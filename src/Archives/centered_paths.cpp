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

struct vpath_info_t {
    double distance;
    std::vector<int> path;
    bool contains_v;
    std::vector<double> dist_to_v;  // Added to track distances to v for each vertex in path
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

extern "C" {
    SEXP S_find_centered_paths(SEXP adj_list_s,
                               SEXP weight_list_s,
                               SEXP vertex_s,
                               SEXP search_radius_s,
                               SEXP target_length_s,
                               SEXP length_tolerance_s,
                               SEXP min_path_length_s);
}



// Helper structure to store bounded path information
//
// It is meant to be used in std::vector<bounded_path_info_t> that replaces
// std::vector<int> prev - storing information about previous vertex for all
// paths starting at some target vertex; in case of large graphs and only a
// handful of paths in a vicinity of the target vertex that we are considering
// std::vector<bounded_path_info_t> is much more efficient
struct bounded_path_info_t {
    int vertex;           // The vertex id
    double distance;      // Distance from start to this vertex
    int prev;             // Previous vertex in path
};

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
#endif

#if 0
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
#endif

/**
 * @brief Find paths containing a target vertex that meet specific criteria using bounded searches
 *        Optimized to process vertices from furthest to closest to enable early termination
 *
 * @details
 *
 *
 *  The function processes potential start vertices in order of decreasing
 *          distance from the target vertex v. If no valid paths are found starting
 *          from a vertex at distance d, all vertices closer than d are skipped
 *          as they cannot form valid paths either.
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param weight_list Weight list where weight_list[i] contains the weights of edges
 *                    from vertex i to its adjacent vertices in adj_list[i]
 * @param v Target vertex that each shortest paths must include to be considered for selection
 * @param target_length Desired total path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 *                        (e.g., 0.1 means ±10% of target_length)
 * @param min_path_length Minimum number of vertices required in valid paths
 *
 * @return std::vector<path_t> Vector of path candidates that meet all criteria
 *
 */
static std::vector<path_t> find_paths_with_radius(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int v,
    double target_length,
    double length_tolerance,
    int min_path_length
) {
    std::vector<path_t> candidates;


    // Find all vertices within search_radius of v
    auto vertices_near_v = find_shortest_paths_within_radius(
        adj_list, weight_list, v, search_radius + length_tolerance);


    // debugging ------
    Rprintf("In find_paths_with_radius()\n");
    print_vect_vect(adj_list, "adj_list");
    print_vect_vect(weight_list, "weight_list");
    print_bounded_paths(vertices_near_v, "vertices_near_v");

    // Sort by distance from v, furthest first
    std::sort(vertices_near_v.begin(), vertices_near_v.end(),
        [](const bounded_path_info_t& a, const bounded_path_info_t& b) {
            return a.distance > b.distance;
        });

    // debugging ------
    print_bounded_paths(vertices_near_v, "vertices_near_v after sorting");


    // Keep track of minimum distance that could possibly work
    double min_working_distance = 0.0;

    for (const auto& start_info : vertices_near_v) {
        // If this vertex is closer than our minimum working distance,
        // all remaining vertices will be too close to form valid paths
        if (start_info.distance < min_working_distance) {
            break;  // Early termination
        }

        bool found_valid_path = false;

        // Find all vertices reachable from this start point within (target_length + length_tolerance) radius
        auto vertices_from_start = find_shortest_paths_within_radius(
            adj_list, weight_list, start_info.vertex, target_length + length_tolerance);

        // Sort end vertices by distance too for consistent ordering
        std::sort(vertices_from_start.begin(), vertices_from_start.end(),
            [](const bounded_path_info_t& a, const bounded_path_info_t& b) {
                return a.distance > b.distance;
            });

        for (const auto& end_info : vertices_from_start) {
            // Skip if start >= end to avoid duplicate paths
            if (start_info.vertex >= end_info.vertex) continue;

            double path_length = end_info.distance;

            // Skip if path length is not within tolerance of target length
            if (std::abs(path_length - target_length) > target_length * length_tolerance) {
                continue;
            }

            // Reconstruct the path and verify it includes vertex v
            try {
                auto path = reconstruct_bounded_path(
                    start_info.vertex, end_info.vertex, vertices_from_start);

                // Skip if path is too short or doesn't include target vertex
                if (path.size() <= min_path_length ||
                    std::find(path.begin(), path.end(), v) == path.end()) {
                    continue;
                }

                // Create and store the candidate path
                path_t candidate;
                candidate.vertices = path;
                candidate.target_vertex = v;
                candidate.total_weight = path_length;
                candidate.center_offset = calculate_relative_center_offset(v, path);
                candidates.push_back(candidate);

                found_valid_path = true;
            } catch (const std::runtime_error& e) {
                // Skip paths that can't be reconstructed
                continue;
            }
        }

        // If we didn't find any valid paths from this distance,
        // update min_working_distance as no closer vertices can work
        if (!found_valid_path) {
            min_working_distance = start_info.distance;
        }
    }

    return candidates;
}





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

/**
 * @brief Calculate how centered a vertex is in a path
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
static double calculate_center_offset(int v, const std::vector<int>& path) {
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
    const std::vector<bounded_path_info_t>& path_info
) {
    std::vector<int> path;
    int curr = end;

    // Find end vertex in path_info
    auto find_vertex = [&path_info](int v) -> const bounded_path_info_t* {
        auto it = std::find_if(path_info.begin(), path_info.end(),
            [v](const bounded_path_info_t& info) { return info.vertex == v; });
        return it != path_info.end() ? &(*it) : nullptr;
    };

    while (curr != -1) {
        const bounded_path_info_t* curr_info = find_vertex(curr);
        if (!curr_info) {
            error("Path reconstruction failed: vertex not found in path info");
        }
        path.push_back(curr);
        curr = curr_info->prev;
    }

    std::reverse(path.begin(), path.end());
    if (path.front() != start) {
        error("Path reconstruction failed: path does not start at the specified vertex");
    }
    return path;
}

/**
 * @brief Find paths containing a target vertex that meet specific criteria using bounded searches
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param weight_list Weight list where weight_list[i] contains the weights of edges
 *                    from vertex i to its adjacent vertices in adj_list[i]
 * @param v Target vertex that must be included in the paths
 * @param search_radius Maximum allowed distance from vertex v to path endpoints
 * @param target_length Desired total path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 *                        (e.g., 0.1 means ±10% of target_length)
 * @param min_path_length Minimum number of vertices required in valid paths
 *
 * @return std::vector<path_t> Vector of path candidates that meet all criteria
 *
 * @details The function efficiently finds paths using bounded searches:
 *          1. Finds all vertices within search_radius of target vertex v
 *          2. For each such vertex (potential start point):
 *             - Finds all vertices within search_radius of this start point
 *             - For each reachable vertex (potential end point):
 *               * Verifies path length is within tolerance of target_length
 *               * Reconstructs the path and checks it includes vertex v
 *               * Verifies minimum path length requirement
 *               * Calculates centrality of vertex v in the path
 *          3. Returns all valid paths as path_t structures
 *
 * @note The function only considers paths where start < end to avoid duplicate paths
 *       in opposite directions
 * @note Uses bounded Dijkstra searches for efficiency, only exploring vertices
 *       within search_radius
 *
 * @pre Graph must be weighted and directed. For undirected graphs, each edge should
 *      appear in both directions in the adjacency list
 * @pre Edge weights must be non-negative
 * @pre search_radius, target_length, and length_tolerance must be non-negative
 * @pre min_path_length must be positive
 * @pre weight_list[i].size() must equal adj_list[i].size() for all i
 *
 * @see bounded_path_info_t for the structure containing path information
 * @see find_shortest_paths_within_radius for the underlying bounded path search
 * @see reconstruct_bounded_path for path reconstruction details
 * @see path_t for the structure containing path information
 * @see calculate_center_offset for how path centrality is computed
 */
static std::vector<bounded_path_info_t> find_shortest_paths_within_radius(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start,
    double radius
) {
    std::vector<bounded_path_info_t> result;
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

        result.push_back({u, d, prev[u]});

        for (size_t i = 0; i < adj_list[u].size(); i++) {
            int v = adj_list[u][i]; // crash linen <<---
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


#if 0
/**
 * @brief Find paths that contain a target vertex and meet specific criteria
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param weight_list Weight list where weight_list[i] contains the weights of edges
 *                    from vertex i to its adjacent vertices in adj_list[i]
 * @param v Target vertex that must be included in the paths
 * @param search_radius Maximum allowed distance from vertex v to path endpoints
 * @param target_length Desired total path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 *                        (e.g., 0.1 means ±10% of target_length)
 * @param min_path_length Minimum number of vertices required in valid paths
 *
 * @return std::vector<path_t> Vector of path candidates that meet all criteria
 *
 * @details The function finds paths by:
 *          1. Computing shortest paths from vertex v to all other vertices
 *          2. Identifying potential endpoints within search_radius of v
 *          3. For each pair of endpoints (start, end):
 *             - Computes the shortest path between them
 *             - Checks if the path length is within tolerance of target_length
 *             - Verifies the path includes vertex v and meets minimum length
 *             - Calculates how centered v is in the path
 *          4. Returns all paths meeting the criteria as path_t structures
 *
 * @note The function only considers paths where start < end to avoid duplicate paths
 *       in opposite directions
 *
 * @see path_t for the structure containing path information
 * @see calculate_center_offset for how path centrality is computed
 * @see find_all_shortest_paths_from_vertex for the underlying shortest path algorithm
 */
static std::vector<path_t> find_paths_with_radius_not_optimized(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int v,
    double search_radius,
    double target_length,
    double length_tolerance,
    int min_path_length
) {
    std::vector<path_t> candidates;

    // Find all vertices within search_radius of v
    auto vertices_near_v = find_shortest_paths_within_radius(
        adj_list, weight_list, v, search_radius);

    // Try all pairs of vertices within radius as potential endpoints
    for (const auto& start_info : vertices_near_v) {
        // Find all vertices reachable from this start point within radius
        auto vertices_from_start = find_shortest_paths_within_radius(
            adj_list, weight_list, start_info.vertex, search_radius);

        for (const auto& end_info : vertices_from_start) {
            // Skip if start >= end to avoid duplicate paths
            if (start_info.vertex >= end_info.vertex) continue;

            double path_length = end_info.distance;

            // Skip if path length is not within tolerance of target length
            if (std::abs(path_length - target_length) > target_length * length_tolerance) {
                continue;
            }

            // Reconstruct the path and verify it includes vertex v
            auto path = reconstruct_bounded_path(
                start_info.vertex, end_info.vertex, vertices_from_start);

            // Skip if path is too short or doesn't include target vertex
            if (path.size() <= min_path_length ||
                std::find(path.begin(), path.end(), v) == path.end()) {
                continue;
            }

            // Create and store the candidate path
            path_t candidate;
            candidate.vertices = path;
            candidate.target_vertex = v;
            candidate.total_weight = path_length;
            candidate.center_offset = calculate_center_offset(v, path);
            candidates.push_back(candidate);
        }
    }

    return candidates;
}
#endif


/**
 * @brief Find paths containing a target vertex that meet specific criteria using bounded searches
 *        Optimized to process vertices from furthest to closest to enable early termination
 *
 * @details
 *
 *
 *  The function processes potential start vertices in order of decreasing
 *          distance from the target vertex v. If no valid paths are found starting
 *          from a vertex at distance d, all vertices closer than d are skipped
 *          as they cannot form valid paths either.
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param weight_list Weight list where weight_list[i] contains the weights of edges
 *                    from vertex i to its adjacent vertices in adj_list[i]
 * @param v Target vertex that each shortest paths must include to be considered for selection
 * @param search_radius Maximum allowed distance from vertex v to path endpoints
 * @param target_length Desired total path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 *                        (e.g., 0.1 means ±10% of target_length)
 * @param min_path_length Minimum number of vertices required in valid paths
 *
 * @return std::vector<path_t> Vector of path candidates that meet all criteria
 *
 */
static std::vector<path_t> find_paths_with_radius(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int v,
    double search_radius,
    double target_length,
    double length_tolerance,
    int min_path_length
) {
    std::vector<path_t> candidates;

    // Find all vertices within search_radius of v
    auto vertices_near_v = find_shortest_paths_within_radius(
        adj_list, weight_list, v, search_radius + length_tolerance);


    // debugging ------
    Rprintf("In find_paths_with_radius()\n");
    print_vect_vect(adj_list, "adj_list");
    print_vect_vect(weight_list, "weight_list");
    print_bounded_paths(vertices_near_v, "vertices_near_v");

    // Sort by distance from v, furthest first
    std::sort(vertices_near_v.begin(), vertices_near_v.end(),
        [](const bounded_path_info_t& a, const bounded_path_info_t& b) {
            return a.distance > b.distance;
        });

    // debugging ------
    print_bounded_paths(vertices_near_v, "vertices_near_v after sorting");


    // Keep track of minimum distance that could possibly work
    double min_working_distance = 0.0;

    for (const auto& start_info : vertices_near_v) {
        // If this vertex is closer than our minimum working distance,
        // all remaining vertices will be too close to form valid paths
        if (start_info.distance < min_working_distance) {
            break;  // Early termination
        }

        bool found_valid_path = false;

        // Find all vertices reachable from this start point within (target_length + length_tolerance) radius
        auto vertices_from_start = find_shortest_paths_within_radius(
            adj_list, weight_list, start_info.vertex, target_length + length_tolerance);

        // Sort end vertices by distance too for consistent ordering
        std::sort(vertices_from_start.begin(), vertices_from_start.end(),
            [](const bounded_path_info_t& a, const bounded_path_info_t& b) {
                return a.distance > b.distance;
            });

        for (const auto& end_info : vertices_from_start) {
            // Skip if start >= end to avoid duplicate paths
            if (start_info.vertex >= end_info.vertex) continue;

            double path_length = end_info.distance;

            // Skip if path length is not within tolerance of target length
            if (std::abs(path_length - target_length) > target_length * length_tolerance) {
                continue;
            }

            // Reconstruct the path and verify it includes vertex v
            try {
                auto path = reconstruct_bounded_path(
                    start_info.vertex, end_info.vertex, vertices_from_start);

                // Skip if path is too short or doesn't include target vertex
                if (path.size() <= min_path_length ||
                    std::find(path.begin(), path.end(), v) == path.end()) {
                    continue;
                }

                // Create and store the candidate path
                path_t candidate;
                candidate.vertices = path;
                candidate.target_vertex = v;
                candidate.total_weight = path_length;
                candidate.center_offset = calculate_center_offset(v, path);
                candidates.push_back(candidate);

                found_valid_path = true;
            } catch (const std::runtime_error& e) {
                // Skip paths that can't be reconstructed
                continue;
            }
        }

        // If we didn't find any valid paths from this distance,
        // update min_working_distance as no closer vertices can work
        if (!found_valid_path) {
            min_working_distance = start_info.distance;
        }
    }

    return candidates;
}

/**
 * @brief R interface for finding centered paths in a graph
 *
 * @param adj_list_R List of integer vectors representing graph adjacency (1-based indices)
 * @param weight_list_R List of numeric vectors containing edge weights
 * @param vertex_R Integer scalar specifying target vertex (1-based index)
 * @param bandwidth_R Numeric scalar specifying desired half-length of paths
 * @param tolerance_R Numeric scalar specifying acceptable deviation from target length
 * @param max_paths_R Integer scalar specifying maximum number of paths to return
 *
 * @return SEXP (list) where each element is a list containing:
 *         - vertices: integer vector of path vertices (1-based indices)
 *         - total_weight: numeric scalar of path length
 *         - center_offset: numeric scalar indicating how centered the target vertex is
 */
SEXP S_find_centered_paths(SEXP adj_list_s,
                           SEXP weight_list_s,
                           SEXP vertex_s,
                           SEXP search_radius_s,
                           SEXP target_length_s,
                           SEXP length_tolerance_s,
                           SEXP min_path_length_s) {

    std::vector<std::vector<int>> adj_list        = Rgraph_to_vector(adj_list_s);
    std::vector<std::vector<double>> weight_list  = Rweights_to_vector(weight_list_s);

    int vertex = INTEGER(vertex_s)[0];
    double search_radius = REAL(search_radius_s)[0];
    double target_length = REAL(target_length_s)[0];
    double length_tolerance = REAL(length_tolerance_s)[0];
    int min_path_length = INTEGER(min_path_length_s)[0];

    // Create path finder and get results
    std::vector<path_t> paths = find_paths_with_radius(
        adj_list,
        weight_list,
        vertex,
        search_radius,
        target_length,
        length_tolerance,
        min_path_length);

    // Convert results to R list
    SEXP result = PROTECT(allocVector(VECSXP, paths.size()));

    // debugging
    Rprintf("In S_find_centered_paths(): ");
    Rprintf("vertex: %d\n", vertex);
    Rprintf("Number of paths found: %d\n", (int)paths.size());


    for(size_t i = 0; i < paths.size(); i++) {
        SEXP path = PROTECT(allocVector(VECSXP, 3));

        // Convert vertices to R (adding 1 to convert to 1-based indexing)
        SEXP vertices = PROTECT(allocVector(INTSXP, paths[i].vertices.size()));
        for(size_t j = 0; j < paths[i].vertices.size(); j++) {
            INTEGER(vertices)[j] = paths[i].vertices[j] + 1; // changing to 1-based integers
        }

        // Create weight and offset elements
        SEXP total_weight = PROTECT(ScalarReal(paths[i].total_weight));
        SEXP center_offset = PROTECT(ScalarReal(paths[i].center_offset));

        // debugging
        Rprintf("\n\n--------\ni: %d\n",(int)i);
        print_vect(paths[i].vertices,"paths[i].vertices");
        Rprintf("total_weight: %f\n", paths[i].total_weight);
        Rprintf("center_offset: %f\n", paths[i].center_offset);

        // Combine into path list
        SET_VECTOR_ELT(path, 0, vertices);
        SET_VECTOR_ELT(path, 1, total_weight);
        SET_VECTOR_ELT(path, 2, center_offset);

        // Add to result
        SET_VECTOR_ELT(result, i, path);
        UNPROTECT(4);
    }

    UNPROTECT(1);

    return result;
}


///
///
///

/**
 * @brief Finds v-independent paths with v positioned close to center
 *
 * @param initial_v_independent_paths Input vector of v-independent paths from find_v_independent_paths. Those are v-independent paths all starting at v that have length bound by max_length
 * @param v Target vertex that must be in all paths
 * @param target_length Desired path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 * @param bfs_map Map from vertex to (distance, predecessor) pair from find_shortest_paths_within_radius
 * @param relative_center_offset_tolerance Tolerance for selecting paths with near-optimal center positioning
 *
 * @return Vector of vectors of integers, where each inner vector represents a v-independent path
 *         with v positioned near the center
 *
 * @details The function implements a two-phase algorithm:
 *          1. For pairs of paths from initial_v_independent_paths whose combined length is near target_length:
 *             - Finds shortest path between their endpoints
 *             - Verifies v is included and length meets criteria
 *          2. From valid paths found:
 *             - Calculates center offset for v in each path
 *             - Selects paths where v is optimally centered within tolerance
 *
 * @note If initial_v_independent_paths contains only one path, returns that path
 * @note Returns empty vector if no paths meeting criteria are found
 * @note All returned paths contain v and are v-independent
 *
 * @see calculate_relative_center_offset for center position calculation
 * @see find_v_independent_paths for generation of input paths
 */
static std::vector<std::vector<int>> find_centered_v_independent_paths(
    const std::vector<std::vector<int>>& initial_v_independent_paths,
    int v,
    double target_length,
    double length_tolerance,
    const std::unordered_map<int, std::pair<double, int>>& bfs_map,
    double relative_center_offset_tolerance
) {
    if (initial_v_independent_paths.empty()) return std::vector<std::vector<int>>();
    if (initial_v_independent_paths.size() == 1) return initial_v_independent_paths;

    // Structure to store path info with center offset
    struct centered_path_t {
        std::vector<int> vertices;
        double center_offset;

        centered_path_t(const std::vector<int>& v, double c)
            : vertices(v), center_offset(c) {}
    };

    std::vector<centered_path_t> centered_paths;
    double length_eps = target_length * length_tolerance;

    // For each pair of paths in initial_v_independent_paths
    for (size_t i = 0; i < initial_v_independent_paths.size(); i++) {
        for (size_t j = i + 1; j < initial_v_independent_paths.size(); j++) {
            const auto& path_i = initial_v_independent_paths[i];
            const auto& path_j = initial_v_independent_paths[j];

            // Get endpoints different from v
            int e_i = (path_i.front() == v) ? path_i.back() : path_i.front();
            int e_j = (path_j.front() == v) ? path_j.back() : path_j.front();

            // Check if path between endpoints exists in bfs_map
            auto it_j = bfs_map.find(e_j);
            if (it_j == bfs_map.end()) continue;

            // Reconstruct path between endpoints
            std::vector<int> path_ij = reconstruct_bounded_path(e_i, e_j, bfs_map);

            // Check if path contains v
            auto v_pos = std::find(path_ij.begin(), path_ij.end(), v);
            if (v_pos == path_ij.end()) continue;

            // Check path length
            auto it_length = bfs_map.find(e_j);
            if (it_length == bfs_map.end()) continue;
            double path_length = it_length->second.first;

            if (std::abs(path_length - target_length) <= length_eps) {
                double offset = calculate_relative_center_offset(v, path_ij);
                centered_paths.emplace_back(path_ij, offset);
            }
        }
    }

    if (centered_paths.empty()) return std::vector<std::vector<int>>();

    // Find minimum center offset
    double min_offset = std::numeric_limits<double>::infinity();
    for (const auto& path : centered_paths) {
        min_offset = std::min(min_offset, path.center_offset);
    }

    // Select paths with offset within tolerance
    std::vector<std::vector<int>> result;
    for (const auto& path : centered_paths) {
        if (path.center_offset <= min_offset + relative_center_offset_tolerance) {
            result.push_back(path.vertices);
        }
    }

    return result;
}



/**
 * @brief Finds v-independent paths from bounded Dijkstra search results
 *
 * @param bfs_map Map from vertex to (distance, predecessor) pair from find_shortest_paths_within_radius
 * @param v Target vertex that must be in all paths
 * @param max_length Maximum allowed path length
 *
 * @return Vector of vectors of integers, where each inner vector represents a v-independent path
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
 */
static std::vector<std::vector<int>> find_v_independent_paths(
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
    std::vector<path_info_t> paths;
    for (const auto& info : vertex_info) {
        if (!info.available) continue;

        // Reconstruct path from v to current vertex
        std::vector<int> path = reconstruct_bounded_path(v, info.vertex, bfs_map);

        // Check path length
        if (path.size() < 2) continue;  // Skip trivial paths

        // Get path length from endpoint's distance in bfs_map
        auto it = bfs_map.find(info.vertex);
        if (it == bfs_map.end()) {
            error("Path vertex not found in distance map");
        }
        double path_length = it->second.first;  // distance is stored in first

        if (path_length > max_length) continue;

        // Mark vertices in path as unavailable
        for (int vertex : path) {
            for (auto& v_info : vertex_info) {
                if (v_info.vertex == vertex) {
                    v_info.available = false;
                    break;
                }
            }
        }

        // Add to paths
        paths.emplace_back(path, path_length);
    }

    // Filter for v-independence
    std::vector<path_info_t> independent_paths;
    for (size_t i = 0; i < paths.size(); i++) {
        bool is_independent = true;
        // Remove v from vertex set for intersection check
        paths[i].vertex_set.erase(v);

        // Check against all previously accepted independent paths
        for (const auto& ind_path : independent_paths) {
            // Check intersection
            for (int vertex : paths[i].vertex_set) {
                if (ind_path.vertex_set.count(vertex) > 0) {
                    is_independent = false;
                    break;
                }
            }
            if (!is_independent) break;
        }

        if (is_independent) {
            independent_paths.push_back(paths[i]);
        }

        // Restore v to vertex set
        paths[i].vertex_set.insert(v);
    }

    // Extract just the vertex sequences for return
    std::vector<std::vector<int>> result;
    result.reserve(independent_paths.size());
    for (const auto& path : independent_paths) {
        result.push_back(path.vertices);
    }

    return result;
}


/**
 * @brief Finds v-independent paths with v positioned close to center
 *
 * @param initial_v_independent_paths Input vector of v-independent paths from find_v_independent_paths
 * @param v Target vertex that must be in all paths
 * @param target_length Desired path length
 * @param length_tolerance Allowed deviation from target_length as a fraction
 * @param bfs_map Map from vertex to (distance, predecessor) pair from find_shortest_paths_within_radius
 * @param relative_center_offset_tolerance Tolerance for selecting paths with near-optimal center positioning
 *
 * @return Vector of path_t structures representing v-independent paths with v positioned near center
 *         Each path_t includes distances from each vertex to the target vertex v
 *
 * @details The function implements a two-phase algorithm:
 *          1. For pairs of paths from initial set whose combined length is near target_length:
 *             - Finds shortest path between their endpoints
 *             - Verifies v is included and length meets criteria
 *          2. From valid paths found:
 *             - Calculates center offset for v in each path
 *             - Computes distances from each vertex to target v
 *             - Selects paths where v is optimally centered within tolerance
 *
 * @note If input contains only one path, returns that path as path_t
 * @note Returns empty vector if no paths meeting criteria are found
 * @note All returned paths contain v and are v-independent
 *
 * @see path_t for the structure containing path and distance information
 * @see calculate_relative_center_offset for center position calculation
 * @see find_v_independent_paths for generation of input paths
 */
static std::vector<path_t> find_centered_v_independent_paths(
    const std::vector<std::vector<int>>& initial_v_independent_paths,
    int v,
    double target_length,
    double length_tolerance,
    const std::unordered_map<int, std::pair<double, int>>& bfs_map,
    double relative_center_offset_tolerance
) {
    if (initial_v_independent_paths.empty()) return std::vector<path_t>();

    if (initial_v_independent_paths.size() == 1) {
        path_t single_path;
        single_path.vertices = initial_v_independent_paths[0];
        single_path.target_vertex = v;
        single_path.center_offset = calculate_relative_center_offset(v, single_path.vertices);

        // Calculate total weight and distances for single path
        auto end_it = bfs_map.find(single_path.vertices.back());
        if (end_it != bfs_map.end()) {
            single_path.total_weight = end_it->second.first;
        }

        // Calculate distances to target for each vertex
        single_path.dist_to_target.reserve(single_path.vertices.size());
        for (int vertex : single_path.vertices) {
            auto it = bfs_map.find(vertex);
            single_path.dist_to_target.push_back(
                it != bfs_map.end() ? it->second.first : INFINITY
            );
        }

        return std::vector<path_t>{single_path};
    }

    // Temporary structure for path construction
    struct path_candidate_t {
        std::vector<int> vertices;
        double total_weight;
        double center_offset;

        path_candidate_t(const std::vector<int>& v, double w, double c)
            : vertices(v), total_weight(w), center_offset(c) {}
    };

    std::vector<path_candidate_t> path_candidates;
    double length_eps = target_length * length_tolerance;

    // For each pair of paths
    for (size_t i = 0; i < initial_v_independent_paths.size(); i++) {
        for (size_t j = i + 1; j < initial_v_independent_paths.size(); j++) {
            const auto& path_i = initial_v_independent_paths[i];
            const auto& path_j = initial_v_independent_paths[j];

            // Get endpoints different from v
            int e_i = (path_i.front() == v) ? path_i.back() : path_i.front();
            int e_j = (path_j.front() == v) ? path_j.back() : path_j.front();

            auto it_j = bfs_map.find(e_j);
            if (it_j == bfs_map.end()) continue;

            // Reconstruct path between endpoints
            std::vector<int> path_ij = reconstruct_bounded_path(e_i, e_j, bfs_map);

            // Verify v is in path and get path length
            if (std::find(path_ij.begin(), path_ij.end(), v) == path_ij.end()) continue;

            double path_length = it_j->second.first;
            if (std::abs(path_length - target_length) <= length_eps) {
                double offset = calculate_relative_center_offset(v, path_ij);
                path_candidates.emplace_back(path_ij, path_length, offset);
            }
        }
    }

    if (path_candidates.empty()) return std::vector<path_t>();

    // Find minimum center offset
    double min_offset = std::numeric_limits<double>::infinity();
    for (const auto& path : path_candidates) {
        min_offset = std::min(min_offset, path.center_offset);
    }

    // Convert qualifying candidates to path_t
    std::vector<path_t> result;
    for (const auto& candidate : path_candidates) {
        if (candidate.center_offset <= min_offset + relative_center_offset_tolerance) {
            path_t path;
            path.vertices = candidate.vertices;
            path.target_vertex = v;
            path.total_weight = candidate.total_weight;
            path.center_offset = candidate.center_offset;

            // Calculate distances to target for each vertex
            path.dist_to_target.reserve(path.vertices.size());
            for (int vertex : path.vertices) {
                auto it = bfs_map.find(vertex);
                path.dist_to_target.push_back(
                    it != bfs_map.end() ? it->second.first : INFINITY
                );
            }

            result.push_back(path);
        }
    }

    return result;
}



/**
 * @brief R interface for finding centered paths in a graph
 *
 * @param adj_list_R List of integer vectors representing graph adjacency (1-based indices)
 * @param weight_list_R List of numeric vectors containing edge weights
 * @param vertex_R Integer scalar specifying target vertex (1-based index)
 * @param bandwidth_R Numeric scalar specifying desired half-length of paths
 * @param tolerance_R Numeric scalar specifying acceptable deviation from target length
 * @param max_paths_R Integer scalar specifying maximum number of paths to return
 *
 * @return SEXP (list) where each element is a list containing:
 *         - vertices: integer vector of path vertices (1-based indices)
 *         - total_weight: numeric scalar of path length
 *         - center_offset: numeric scalar indicating how centered the target vertex is
 */
SEXP S_find_centered_paths(SEXP adj_list_s,
                           SEXP weight_list_s,
                           SEXP vertex_s,
                           SEXP search_radius_s,
                           SEXP target_length_s,
                           SEXP length_tolerance_s,
                           SEXP min_path_length_s) {

    std::vector<std::vector<int>> adj_list        = Rgraph_to_vector(adj_list_s);
    std::vector<std::vector<double>> weight_list  = Rweights_to_vector(weight_list_s);

    int vertex = INTEGER(vertex_s)[0];
    double search_radius = REAL(search_radius_s)[0];
    double target_length = REAL(target_length_s)[0];
    double length_tolerance = REAL(length_tolerance_s)[0];
    int min_path_length = INTEGER(min_path_length_s)[0];

    // Create path finder and get results
    std::vector<path_t> paths = find_paths_with_radius(
        adj_list,
        weight_list,
        vertex,
        search_radius,
        target_length,
        length_tolerance,
        min_path_length);

    // Convert results to R list
    SEXP result = PROTECT(allocVector(VECSXP, paths.size()));

    // debugging
    Rprintf("In S_find_centered_paths(): ");
    Rprintf("vertex: %d\n", vertex);
    Rprintf("Number of paths found: %d\n", (int)paths.size());


    for(size_t i = 0; i < paths.size(); i++) {
        SEXP path = PROTECT(allocVector(VECSXP, 3));

        // Convert vertices to R (adding 1 to convert to 1-based indexing)
        SEXP vertices = PROTECT(allocVector(INTSXP, paths[i].vertices.size()));
        for(size_t j = 0; j < paths[i].vertices.size(); j++) {
            INTEGER(vertices)[j] = paths[i].vertices[j] + 1; // changing to 1-based integers
        }

        // Create weight and offset elements
        SEXP total_weight = PROTECT(ScalarReal(paths[i].total_weight));
        SEXP center_offset = PROTECT(ScalarReal(paths[i].center_offset));

        // debugging
        Rprintf("\n\n--------\ni: %d\n",(int)i);
        print_vect(paths[i].vertices,"paths[i].vertices");
        Rprintf("total_weight: %f\n", paths[i].total_weight);
        Rprintf("center_offset: %f\n", paths[i].center_offset);

        // Combine into path list
        SET_VECTOR_ELT(path, 0, vertices);
        SET_VECTOR_ELT(path, 1, total_weight);
        SET_VECTOR_ELT(path, 2, center_offset);

        // Add to result
        SET_VECTOR_ELT(result, i, path);
        UNPROTECT(4);
    }

    UNPROTECT(1);

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



static std::vector<path_t> find_v_centered_paths(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int v,
    double target_length,
    double length_tolerance,
    double relative_center_offset_tolerance
    ) {

    double max_radius = target_length + length_tolerance * target_length;

    // debugging
    Rprintf("\n\nIn find_v_centered_paths()\n");
    Rprintf("max_radius: %f\n", max_radius);

    auto bfs_map = find_shortest_paths_within_radius(
        adj_list,
        weight_list,
        v,
        max_radius);

    // debugging
    print_vertex_path_map(bfs_map, "result of find_shortest_paths_within_radius()");

    double min_target_length = target_length * (1 - length_tolerance);
    double max_target_length = target_length * (1 + length_tolerance);

    Rprintf("min_target_length: %f\n", min_target_length);
    Rprintf("max_target_length: %f\n", max_target_length);

    std::vector<path_t> paths = reconstruct_paths(v, bfs_map);
    int n_paths = paths.size();

    print_vect_paths(paths, "reconstructed paths from bfs_map");

    std::vector<path_t> centered_paths;

    if (n_paths == 1) {
        return paths;
    } else if (n_paths > 1) {

        std::unordered_set<std::pair<int,int>> candidate_vertex_pairs; // a set of BFS vertices from vertex_info so that d_i + d_j is within the range [min_target_length, max_target_length]

        int i = 0;
        while (i < n_paths) {
            int e_i = paths[i].vertices[paths[i].vertices.size() - 1]; // end vertex of the i-th path
            int d_i = paths[i].total_weight;
            int j = i + 1;

            candidate_vertex_pairs.clear();

            while (j < n_paths) {
                int e_j = paths[j].vertices[paths[j].vertices.size() - 1];
                int d_j = paths[j].total_weight;

                // The shortest path between e_i and e_j may be shorter than the
                // union of paths joining e_i with v and v with e_j. Thus, the
                // length of this shortest path may be less than d_i + d_j. The
                // shortest path between e_i and e_j has the length less than d_i +
                // d_j if and only if v is not a part of that path. Since we are
                // looking for shortest paths between vertices from the BFS nbhd of
                // v, we could have computed the shortest path between them and
                // check if its length is d_i + d_j. But since we are looking for
                // paths passing through v of length in the range
                //
                // [target_length * (1 - length_tolerance), target_length * (1 + length_tolerance)]
                //
                // it would be waistful to do that calculation if d_i + d_j was
                // outside of this range.
                //
                // Moreover, since vertex_info is sorted in the descending order of
                // distances from v
                //
                // if d_i + d_j < target_length * (1 - length_tolerance)
                //
                // then this condition will be true for all j' > j and so we break
                // the inner while loop if this is the case

                if (d_i + d_j < min_target_length) {
                    i++;
                    break;
                } else if (d_i + d_j > max_target_length) {
                    j++;
                    continue;
                } else {
                    candidate_vertices.insert(e_j);
                }
                j++;
            }

            Rprintf("\n\ni: %d\n", i);
            print_uset(candidate_vertices, "candidate_vertices");

            // Find shortest paths paths between e_i and candidate_vertices
            // filtering out paths that do not contain v
            std::map<int, path_t> good_shortest_paths = find_v_shortest_paths(
                adj_list,
                weight_list,
                e_i,
                candidate_vertices,
                v,
                d_i);

            print_map_paths(good_shortest_paths, "good_shortest_paths");

            for (const auto& [w, path] : good_shortest_paths) {
                centered_paths.emplace_back(path);
            }

            print_vect_paths(centered_paths, "centered_paths");

            i++;
        }


    } else {
        error("No paths were found");
    }



    if (centered_paths.empty()) return std::vector<path_t>();

    // Find minimum center offset
    double min_offset = std::numeric_limits<double>::infinity();
    for (const auto& path : centered_paths) {
        min_offset = std::min(min_offset, path.center_offset);
    }

    // First, sort paths by length in descending order
    std::sort(centered_paths.begin(), centered_paths.end(),
              [](const path_t& a, const path_t& b) {
                  return a.total_weight > b.total_weight;
              });

    // Create a vector of flags to mark paths that should be kept
    std::vector<bool> keep_path(centered_paths.size(), true);

    // Check each path against longer paths
    for (size_t i = 1; i < centered_paths.size(); ++i) {
        if (!keep_path[i]) continue;  // Skip if already marked for removal

        // Check against all longer paths
        for (size_t j = 0; j < i; ++j) {
            if (!keep_path[j]) continue;  // Skip if the longer path was marked for removal

            if (is_subpath(centered_paths[i], centered_paths[j])) {
                keep_path[i] = false;
                break;
            }
        }
    }

    // Construct final result with non-subpath paths within offset tolerance
    std::vector<path_t> result;
    for (size_t i = 0; i < centered_paths.size(); ++i) {
        if (keep_path[i] &&
            centered_paths[i].center_offset <= min_offset + relative_center_offset_tolerance) {
            result.push_back(centered_paths[i]);
        }
    }

    print_vect_paths(result, "result of find_v_centered_paths()");

    return result;
}


#if 0
/**
 * @brief Checks if one path is a sub-path of another path.
 *
 * A path is considered a sub-path if all its vertices appear in the same order
 * within the potential super-path. The function uses a sliding window approach
 * to efficiently check for sequence matching.
 *
 * @param potential_subpath The path being checked as a potential sub-path
 * @param potential_superpath The longer path being checked as a potential super-path
 *
 * @return true if potential_subpath is a sub-path of potential_superpath,
 *         false otherwise
 *
 * @note The function assumes both paths are valid and contain at least one vertex.
 * @note If potential_subpath is longer than or equal in length to potential_superpath,
 *       the function returns false without performing sequence matching.
 *
 * Example:
 * @code
 * path_t path1 = {vertices: {1, 2, 3}, ...};
 * path_t path2 = {vertices: {1, 2, 3, 4, 5}, ...};
 * bool is_sub = is_subpath(path1, path2);  // Returns true
 * @endcode
 */
static bool is_subpath(const path_t& potential_subpath, const path_t& potential_superpath) {
    if (potential_subpath.vertices.size() >= potential_superpath.vertices.size()) {
        return false;
    }

    // Use sliding window approach to check if smaller path's vertices appear in sequence
    const auto& sub_vertices = potential_subpath.vertices;
    const auto& super_vertices = potential_superpath.vertices;

    for (size_t i = 0; i <= super_vertices.size() - sub_vertices.size(); ++i) {
        bool match = true;
        for (size_t j = 0; j < sub_vertices.size(); ++j) {
            if (super_vertices[i + j] != sub_vertices[j]) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }
    return false;
}
#endif


/**
 * @brief Finds paths of specified length with a target vertex positioned near the center.
 *
 * This function searches for paths in a graph that:
 * 1. Have a total length within the specified tolerance of the target length
 * 2. Contain the specified vertex (v) positioned close to the center
 * 3. Are not sub-paths of other valid paths in the result set
 *
 * @param adj_list Adjacency list representation of the graph where adj_list[i]
 *                 contains the indices of vertices adjacent to vertex i
 * @param weight_list Weight list corresponding to adj_list where weight_list[i][j]
 *                    is the weight of the edge from vertex i to adj_list[i][j]
 * @param v The vertex that should be positioned near the center of the paths
 * @param target_length The desired total length of the paths
 * @param length_tolerance Relative tolerance for path length (actual_length must be
 *                        within target_length ± (target_length * length_tolerance))
 * @param relative_center_offset_tolerance Maximum allowed deviation from the most
 *                                        centered position of v in the path
 *
 * @return Vector of path_t objects representing the found paths. Empty vector if no
 *         paths meeting the criteria are found.
 *
 * @note The function first finds all paths meeting the length and centering criteria,
 *       then filters out any paths that are sub-paths of longer valid paths.
 * @note Paths are considered valid only if they contain the target vertex v.
 * @note The relative center offset is calculated as the absolute difference between
 *       the actual position of v and the middle position, divided by path length.
 *
 * Example:
 * @code
 * std::vector<std::vector<int>> adj_list = {{1, 2}, {0, 2}, {0, 1}};
 * std::vector<std::vector<double>> weight_list = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
 * auto paths = find_v_centered_paths(adj_list, weight_list, 1, 2.0, 0.1, 0.1);
 * @endcode
 */


// Key improvements and modifications:

// Added is_subpath helper function that efficiently checks if one path is contained within another using a sliding window approach.
// Modified the path filtering process to:

// First sort paths by total weight in descending order
// Use a boolean vector to mark paths for keeping/removal
// Check each path against only longer paths (optimization)
// Combine sub-path elimination with center offset filtering in one pass


// Potential areas for further improvement:

// a) Performance optimizations:

// Could use parallel processing for sub-path checking on large path sets
// Could implement early termination in is_subpath checks based on path properties
// Could use hash-based path comparison for faster matching

// b) Memory optimizations:

// Could implement in-place filtering instead of creating new result vector
// Could use bit vector instead of bool vector for keep_path

// c) Robustness improvements:

// Could add tolerance comparison for path weights to handle floating-point imprecision
// Could add validation for path consistency (ensure v is present, check for cycles)

// Would you like me to implement any of these additional improvements?


static std::vector<path_data_t> get_path_data(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    int ref_vertex,
    double bandwidth,
    double dist_normalization_factor = 1.01,
    int min_path_size = 5,
    int kernel_type = 7
    ) {

    initialize_kernel(kernel_type, 1.0);

    auto bfs_map = find_shortest_paths_within_radius(
        adj_list,
        weight_list,
        ref_vertex,
        bandwidth);

    std::vector<path_t> paths = reconstruct_paths(ref_vertex, bfs_map);
    int n_paths = paths.size();

    // debugging
    print_vect_paths(paths, "reconstructed paths from bfs_map");

    std::vector<path_data_t> result;
    std::vector<double> d_path; // used to compute distances to ref_vertex along the given path

    if (n_paths == 1 && paths[0].vertices.size() < min_path_size) {
        REPORT_WARNING("Warning: Local nbhd has only one path and its length is %d - that is less than min_path_size: %d\n", (int)paths[0].vertices.size(), min_path_size);
        return std::vector<path_data_t>();
    } else if (n_paths == 1) {
        path_data_t path_data;
        path_data.ref_vertex = ref_vertex;
        path_data.vertices = std::move(paths[0].vertices);
        path_data.total_weight = paths[0].total_weight;

        path_data.rel_center_offset = calculate_relative_center_offset(ref_vertex, path_data.vertices);

        // Calculating x_path := distance along the path from the initial vertex of the path
        int path_n_vertices = path_data.vertices.size();
        path_data.x_path.resize(path_n_vertices);
        path_data.w_path.resize(path_n_vertices);
        path_data.y_path.resize(path_n_vertices);
        d_path.resize(path_n_vertices);

        // Extract y values, weights, and compute distances
        path_data.x_path[0] = 0;
        for (int k = 0; k < path_n_vertices; ++k) {
            path_data.y_path[k] = y[path_data.vertices[k]];
            if (k > 0) {
                auto neighbors = adj_list[path_data.vertices[k - 1]];
                auto neighbor_weights = weight_list[path_data.vertices[k - 1]];
                auto it = std::find(neighbors.begin(), neighbors.end(), path_data.vertices[k]);
                if (it == neighbors.end()) {
                    REPORT_ERROR("Edge not found in adjacency list between vertices %d and %d",
                                 path_data.vertices[k-1], path_data.vertices[k]);
                }
                int neighbor_k = it - neighbors.begin();
                path_data.x_path[k] = path_data.x_path[k - 1] + neighbor_weights[neighbor_k];
            }
        }

        double x_ref = path_data.x_path[0];
        double max_dist = 0.0;
        for (int k = 0; k < path_n_vertices; ++k) {
            d_path[k] = std::abs(path_data.x_path[k] - x_ref);
            max_dist = std::max(max_dist, d_path[k]);
        }

        if (max_dist == 0) max_dist = 1;
        max_dist *= dist_normalization_factor;

        for (int k = 0; k < path_n_vertices; ++k) {
            d_path[k] /= max_dist;
        }

        kernel_fn(d_path.data(), path_n_vertices, path_data.w_path.data());

        double total_w_path = std::accumulate(path_data.w_path.begin(), path_data.w_path.end(), 0.0);
        for (int k = 0; k < path_n_vertices; ++k)
            path_data.w_path[k] /= total_w_path;

        result.emplace_back(path_data);

    } else if (n_paths > 1) {
        // check if there is a compsite path with at least min_path_size vertices
        int max_size = 0;
        for (int i = 0; i < n_paths; ++i) {
            for (int j = i + 1; j < n_paths; ++j) {
                int comb_path_size = paths[i].vertices.size() + paths[j].vertices.size();
                if (comb_path_size > max_size)
                    max_size = comb_path_size;
            }
        }

        if (max_size < min_path_size) {
            Rprintf("Warning: None of the composed paths has min_path_size: %d vertices\n", min_path_size);
            return std::vector<path_data_t>();
        }

        for (int i = 0; i < n_paths; ++i) {
            for (int j = i + 1; j < n_paths; ++j) {
                int comb_path_n_vertices = paths[i].vertices.size() + paths[j].vertices.size() - 1; // removing ref_vertex from one of the paths
                if (comb_path_n_vertices < min_path_size) continue;

                Rprintf("i: %d  j: %d\n", i, j);
                print_vect(paths[i].vertices, "paths[i].vertices");
                print_vect(paths[j].vertices, "paths[j].vertices");

                path_data_t path_data;
                path_data.ref_vertex = ref_vertex;
                // Combining paths[i].vertices and paths[j].vertices removing ref_vertex from one of them
                // Reserve the total space needed (single allocation)
                path_data.vertices.reserve(
                    paths[i].vertices.size() + paths[j].vertices.size() - 1
                    );
                // Insert j's vertices in reverse order using a reverse iterator
                path_data.vertices.insert(
                    path_data.vertices.end(),
                    paths[j].vertices.rbegin(),  // Reverse iterator start
                    paths[j].vertices.rend()     // Reverse iterator end
                    );
                // Insert i's vertices normally but starting from the second vertex
                path_data.vertices.insert(
                    path_data.vertices.end(),
                    paths[i].vertices.begin() + 1, // skip the first element
                    paths[i].vertices.end()
                    );

                // Compute correctly the total weight by subtracting the weight of the first edge of path_data[i].vertices
                // Find the weight of the first edge in paths[i] (between ref_vertex and the second vertex in the path)
                auto it = std::find(adj_list[ref_vertex].begin(), adj_list[ref_vertex].end(), paths[i].vertices[1]);
                if (it != adj_list[ref_vertex].end()) {  // Add safety check
                    int first_edge_weight_index = std::distance(adj_list[ref_vertex].begin(), it);
                    path_data.total_weight = paths[j].total_weight + paths[i].total_weight - weight_list[ref_vertex][first_edge_weight_index];
                } else {
                    REPORT_ERROR("Edge (%d,%d) not found in adjacency list", ref_vertex, paths[i].vertices[1]);
                }
                // Calculating relative center offset of ref_vertex in path_data.vertices
                path_data.rel_center_offset = calculate_relative_center_offset(ref_vertex, path_data.vertices);
                // Calculating x_path := distance along the path from the initial vertex of the path
                path_data.x_path.resize(comb_path_n_vertices);
                path_data.w_path.resize(comb_path_n_vertices);
                path_data.y_path.resize(comb_path_n_vertices);
                d_path.resize(comb_path_n_vertices);

                // Extract y values, weights, and compute distances
                path_data.x_path[0] = 0;
                for (int k = 0; k < comb_path_n_vertices; ++k) {
                    path_data.y_path[k] = y[path_data.vertices[k]];
                    if (k > 0) {
                        auto neighbors = adj_list[path_data.vertices[k - 1]];
                        auto neighbor_weights = weight_list[path_data.vertices[k - 1]];
                        auto it = std::find(neighbors.begin(), neighbors.end(), path_data.vertices[k]);
                        int neighbor_k = it - neighbors.begin();
                        path_data.x_path[k] = path_data.x_path[k - 1] + neighbor_weights[neighbor_k];
                    }
                }

                int ref_idx = paths[j].vertices.size() - 1; // last position in the reversed paths[j].vertices
                double x_ref = path_data.x_path[ref_idx];
                double max_dist = 0.0;
                for (int k = 0; k < comb_path_n_vertices; ++k) {
                    d_path[k] = std::abs(path_data.x_path[k] - x_ref);
                    max_dist = std::max(max_dist, d_path[k]);
                }

                if (max_dist == 0) max_dist = 1;
                max_dist *= dist_normalization_factor;

                for (int k = 0; k < comb_path_n_vertices; ++k) {
                    d_path[k] /= max_dist;
                }

                kernel_fn(d_path.data(), comb_path_n_vertices, path_data.w_path.data());

                double total_w_path = std::accumulate(path_data.w_path.begin(), path_data.w_path.end(), 0.0);
                for (int k = 0; k < comb_path_n_vertices; ++k)
                    path_data.w_path[k] /= total_w_path;

                result.emplace_back(path_data);
            }
        }
    }

    return result;
}


/**
 * @brief Finds shortest paths from a start vertex to candidate vertices that must pass through a specific vertex.
 *
 * @details This function implements a modified Dijkstra's algorithm to find shortest paths from
 *          a start vertex to a set of candidate vertices. Only paths that pass through vertex 'v'
 *          at distance dist_start_to_v are included in the result. The algorithm includes early
 *          pruning optimization: paths where vertex 'v' hasn't been discovered and the distance
 *          already exceeds dist_start_to_v are discarded.
 *
 * @param adj_list      Adjacency list representing the graph structure. adj_list[i] contains
 *                      the indices of vertices adjacent to vertex i.
 * @param weight_list   List of edge weights corresponding to adj_list. weight_list[i][j] is
 *                      the weight of the edge from vertex i to its j-th neighbor in adj_list[i].
 * @param start         Index of the starting vertex.
 * @param candidate_vertices Set of vertex indices for which we want to find shortest paths.
 *                      Only paths to these vertices that include vertex 'v' will be returned.
 * @param v             Index of the vertex that must be included in each path (target_vertex).
 * @param dist_start_to_v The required distance from 'start' to 'v' in valid paths.
 *
 * @return A map where keys are vertices from candidate_vertices and values are path_t structures
 *         containing detailed information about the paths, including:
 *         - Sequence of vertices in the path
 *         - Target vertex (v)
 *         - Total path weight
 *         - How well v is centered in the path (center_offset)
 *         - Distances from each vertex in the path to v
 *
 * @pre adj_list and weight_list must have the same size (number of vertices).
 * @pre For each vertex i, adj_list[i] and weight_list[i] must have the same size.
 * @pre All indices in adj_list must be valid vertex indices (0 to n-1, where n is the number of vertices).
 * @pre All weights must be non-negative.
 * @pre start, v, and all vertices in candidate_vertices must be valid vertex indices.
 * @pre dist_start_to_v must be the exact distance from start to v in the graph.
 *
 * @throw std::invalid_argument if preconditions are not met.
 *
 * @note The function uses floating-point arithmetic with epsilon comparisons to handle
 *       numerical precision issues.
 *
 * @complexity Time: O((E + V) log V) where V is the number of vertices and E is the number of edges.
 *             Space: O(V) for storing path information.
 */
std::map<int, path_t> find_v_shortest_paths(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start,
    std::unordered_set<int>& candidate_vertices,
    int v,
    double dist_start_to_v) {

    const int n = adj_list.size();
    const double INF = std::numeric_limits<double>::infinity();
    const double EPSILON = 1e-9;
    const double DIST_START_TO_V_PLUS_EPSILON = dist_start_to_v + EPSILON;

    auto cmp = [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        return a.first > b.first;
    };

    std::vector<vpath_info_t> vertex_info(n, {INF, {}, false, {}});
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        decltype(cmp)> pq(cmp);

    vertex_info[start].distance = 0;
    vertex_info[start].path = {start};
    vertex_info[start].dist_to_v = {INF};  // Initialize distance to v
    pq.push({0, start});

    while (!pq.empty()) {
        auto [curr_dist, curr] = pq.top();
        pq.pop();

        if (curr_dist > vertex_info[curr].distance) continue;

        if (!vertex_info[curr].contains_v && curr_dist > DIST_START_TO_V_PLUS_EPSILON) {
            continue;
        }

        if (curr == v) {
            vertex_info[curr].contains_v = true;
            // Update distances to v for all vertices in the current path
            double v_dist = vertex_info[curr].distance;
            for (size_t i = 0; i < vertex_info[curr].path.size(); ++i) {
                int path_vertex = vertex_info[curr].path[i];
                vertex_info[curr].dist_to_v[i] = v_dist - vertex_info[path_vertex].distance;
            }
        }

        for (size_t i = 0; i < adj_list[curr].size(); ++i) {
            int next = adj_list[curr][i];
            double weight = weight_list[curr][i];
            double new_dist = curr_dist + weight;

            if (new_dist < vertex_info[next].distance - EPSILON) {
                vertex_info[next].distance = new_dist;
                vertex_info[next].path = vertex_info[curr].path;
                vertex_info[next].path.push_back(next);
                vertex_info[next].contains_v = vertex_info[curr].contains_v;

                // Copy distances to v and add distance for the new vertex
                vertex_info[next].dist_to_v = vertex_info[curr].dist_to_v;
                if (vertex_info[next].contains_v) {
                    // If v is already in the path, calculate actual distance to v
                    vertex_info[next].dist_to_v.push_back(new_dist - vertex_info[v].distance);
                } else {
                    // If v not found yet, mark as INF
                    vertex_info[next].dist_to_v.push_back(INF);
                }

                pq.push({new_dist, next});
            }
        }
    }

    std::map<int, path_t> result;
    for (int candidate : candidate_vertices) {
        if (vertex_info[candidate].contains_v) {
            path_t path;
            path.target_vertex = v;
            path.vertices = std::move(vertex_info[candidate].path);
            path.total_weight = vertex_info[candidate].distance;
            path.center_offset = calculate_relative_center_offset(v, path.vertices);
            path.dist_to_target = std::move(vertex_info[candidate].dist_to_v);
            result[candidate] = std::move(path);
        }
    }

    print_map_paths(result, "result of find_v_shortest_paths()");

    return result;
}
