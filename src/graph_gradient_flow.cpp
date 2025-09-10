#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <queue>
#include <random>     // for std::mt19937
#include <sstream>
#include <unordered_set>
#include <functional>
#include <map>

#include <filesystem> // For debugging

#include "reachability_map.hpp"
#include "cpp_utils.hpp"
#include "cpp_stats_utils.hpp"
#include "set_wgraph.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp"
#include "gradient_flow.hpp"

extern "C" {
    SEXP S_construct_graph_gradient_flow(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_scale,
        SEXP s_quantile_scale_thld,
        SEXP s_with_trajectories
        );
}

bool break_duplicate_values(std::vector<double>& y, double noise_magnitude = 1e-10);

/**
 * @brief Checks if a vertex is a local minimum or maximum within its radius-neighborhood
 *
 * Determines whether the given vertex is a local extremum by comparing its value
 * with all vertices in its radius-neighborhood (defined by paths_result).
 *
 * @param vertex The vertex to check for extremum property
 * @param paths_result Contains reachable vertices within the radius-neighborhood
 * @param y Vector of values associated with vertices
 *
 * @return std::optional<std::pair<int, bool>> where:
 *         - If vertex is an extremum: pair of (vertex_index, is_maximum)
 *           * is_maximum = true for local maximum
 *           * is_maximum = false for local minimum
 *         - If not an extremum: std::nullopt
 *
 * @pre y must not contain Rf_duplicate values to ensure unique extrema
 * @pre vertex index must be valid (0 <= vertex < y.size())
 * @pre y.size() must match the graph size
 */
std::optional<std::pair<size_t, bool>> set_wgraph_t::check_local_extremum(
    size_t vertex,
    const shortest_paths_t& paths_result,
    const std::vector<double>& y) const {

    bool is_minimum = true;
    bool is_maximum = true;

    #if 0
    // debugging
    size_t sel_vertex = 379;
    if (vertex == sel_vertex)
    {
        Rprintf("In check_local_extremum()\n");
        Rprintf("sel_vertex: %zu\ty[sel_vertex]: %f\n", sel_vertex, y[sel_vertex]);
        Rprintf("neighbor\tdelta\n");
        for (size_t neighbor : paths_result.reachable_vertices) {
            Rprintf("%zu\t%f\n",
                    neighbor, y[neighbor] - y[sel_vertex]);
        }
        Rprintf("\n");
    }
    #endif

    for (size_t neighbor : paths_result.reachable_vertices) {
        if (y[neighbor] < y[vertex]) {
            is_minimum = false;
            break;
        }
    }
    for (size_t neighbor : paths_result.reachable_vertices) {
        if (y[neighbor] > y[vertex]) {
            is_maximum = false;
            break;
        }
    }

    if (is_minimum) return std::make_pair(vertex, false);  // false for minimum
    if (is_maximum) return std::make_pair(vertex, true);   // true for maximum
    return std::nullopt;
}

/**
 * @brief Reconstructs non-overlapping paths from a reference vertex to reachable vertices.
 *
 * This function creates paths from the reference vertex to other vertices using a greedy approach:
 * 1. Processes vertices in order of decreasing distance from the reference vertex
 * 2. For each vertex, reconstructs a path back to the reference vertex
 * 3. Ensures no vertex appears in more than one path
 * 4. Returns paths starting from the reference vertex
 *
 * @param reachability_map A reachability map containing distances, predecessors, and sorted vertices
 * @return A vector of vertex_path_t structures, each representing a non-overlapping path
 */
std::vector<vertex_path_t> set_wgraph_t::reconstruct_graph_paths(
    const reachability_map_t& reachability_map
    ) const {

    std::vector<vertex_path_t> paths;

    // Create a set of vertices to process, ordered in decreasing order by distance
    struct VertexInfoComparator {
        bool operator()(const vertex_info_t& a, const vertex_info_t& b) const {
            return a.distance > b.distance;
        }
    };

    std::set<vertex_info_t, VertexInfoComparator> to_process(
        reachability_map.sorted_vertices.begin(),
        reachability_map.sorted_vertices.end());

    // CXX20 version
    // std::set<vertex_info_t, decltype([](const vertex_info_t& a, const vertex_info_t& b) {
    //     return a.distance > b.distance;
    // })> to_process(reachability_map.sorted_vertices.begin(),
    //               reachability_map.sorted_vertices.end());

    // Track discovered vertices to ensure non-overlapping paths
    std::vector<bool> discovered(adjacency_list.size(), false);

    // Mark reference vertex as discovered
    discovered[reachability_map.ref_vertex] = true;

    while (!to_process.empty()) {
        // Get the furthest unprocessed vertex
        auto current = *to_process.begin();
        to_process.erase(to_process.begin());

        // Skip if already discovered in a previous path
        if (discovered[current.vertex]) {
            continue;
        }

        // Initialize new path
        vertex_path_t path;

        size_t vertex = current.vertex;

        // Reconstruct path from current vertex back to reference vertex
        while (vertex != reachability_map.ref_vertex &&
               vertex != INVALID_VERTEX) {
            // Skip if we hit a vertex that's already in another path
            if (vertex != current.vertex && discovered[vertex]) {
                path.vertices.clear();
                path.dist_to_ref_vertex.clear();
                break;
            }

            path.vertices.push_back(vertex);
            path.dist_to_ref_vertex.push_back(reachability_map.distances.at(vertex));

            // Mark vertex as discovered
            discovered[vertex] = true;

            // Move to predecessor
            vertex = reachability_map.predecessors.at(vertex);
        }

        // If path reconstruction was successful
        if (!path.vertices.empty()) {
            // Add reference vertex to the path
            path.vertices.push_back(reachability_map.ref_vertex);
            path.dist_to_ref_vertex.push_back(0.0);

            // Reverse the path to start from reference vertex
            std::reverse(path.vertices.begin(), path.vertices.end());
            std::reverse(path.dist_to_ref_vertex.begin(), path.dist_to_ref_vertex.end());

            // Remove vertices from to_process that are now part of this path
            auto it = to_process.begin();
            while (it != to_process.end()) {
                if (discovered[it->vertex]) {
                    it = to_process.erase(it);
                } else {
                    ++it;
                }
            }

            paths.push_back(std::move(path));
        }
    }

    return paths;
}

/**
 * @brief Identifies the endpoints of non-overlapping paths from a reference vertex,
 *        including their complete shortest paths.
 *
 * This function finds endpoints of non-overlapping paths from the reference vertex by:
 * 1. Using reconstruct_paths to create non-overlapping paths
 * 2. Extracting the furthest vertex (endpoint) from each path
 * 3. Returning the endpoints with their distances and complete paths from the reference vertex
 *
 * @param reachability_map A reachability map containing distances, predecessors, and sorted vertices
 * @return A vector of vertex_shortest_path_info_t structures containing endpoint vertices,
 *         their distances, and complete paths from the reference vertex
 */
std::vector<vertex_shortest_path_info_t> set_wgraph_t::get_vertex_shortest_paths(
    const reachability_map_t& reachability_map) const {

    // Reconstruct non-overlapping paths
    // std::vector<vertex_path_t> paths = reconstruct_graph_paths(
    //     const_cast<reachability_map_t&>(reachability_map));

    std::vector<vertex_path_t> paths = reconstruct_graph_paths(reachability_map);

    // Extract endpoints with their complete paths
    std::vector<vertex_shortest_path_info_t> endpoints;
    endpoints.reserve(paths.size());  // Pre-allocate for efficiency

    for (const auto& path : paths) {
        // Skip empty paths
        if (path.vertices.empty()) {
            continue;
        }

        // Get the endpoint (furthest vertex in the path)
        size_t endpoint_vertex = path.vertices.back();
        double endpoint_dist = path.dist_to_ref_vertex.back();

        // Add to endpoints list with the complete path
        endpoints.push_back({endpoint_vertex, endpoint_dist, path.vertices});
    }

    return endpoints;
}


/**
 * @brief Computes a reachability map from a reference vertex within a specified radius.
 *
 * This function performs a bounded Dijkstra's algorithm to find all vertices that are
 * reachable from the specified reference vertex within the given radius. For each
 * reachable vertex, it calculates the shortest distance from the reference vertex and
 * records the predecessor in this shortest path.
 *
 * The function maintains three key pieces of information:
 * 1. A mapping of vertices to their shortest distances from the reference vertex
 * 2. A mapping of vertices to their predecessors in the shortest path
 * 3. A sorted list of vertices in descending order of their distances
 *
 * The time complexity is O((V + E) log V) where V is the number of vertices within
 * the radius and E is the number of edges connecting them.
 *
 * @param ref_vertex The reference vertex from which to compute distances
 * @param radius The maximum distance to search (excluding paths longer than this value)
 * @return A reachability_map_t structure containing the computed information
 *
 * @see reachability_map_t
 * @see find_graph_paths_within_radius
 */
reachability_map_t set_wgraph_t::compute_graph_reachability_map(
    size_t ref_vertex,
    double radius) const {

    reachability_map_t result;
    result.ref_vertex = ref_vertex;

    // Implement a bounded Dijkstra's algorithm
    size_t n = adjacency_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);  // Using INVALID_VERTEX as sentinel
    dist[ref_vertex] = 0;

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, ref_vertex});

    while (!pq.empty()) {
        double d = -pq.top().first;  // Negate to get actual distance (priority queue is max-heap)
        size_t u = pq.top().second;
        pq.pop();

        if (d > radius) break;  // Stop if we exceed the radius
        if (d > dist[u]) continue;  // Skip if we've found a better path already

        // Add vertex information to result
        result.distances[u] = d;
        result.predecessors[u] = prev[u];
        if (u != ref_vertex) {  // Don't include reference vertex in sorted_vertices
            result.sorted_vertices.push_back({u, d});
        }

        // Explore neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double w = edge.weight;

            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});  // Negate for priority queue (smaller distances first)
            }
        }
    }

    // Sort vertices by distance in descending order
    std::sort(result.sorted_vertices.begin(), result.sorted_vertices.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    return result;
}

/**
 * @brief Checks for Rf_duplicate values in a vector of doubles, considering numerical precision
 *
 * This function can detect both exact duplicates and near-duplicates within a specified
 * tolerance. It returns information about the duplicates found to help with debugging
 * and potential resolution.
 *
 * @param y Vector of values to check for duplicates
 * @param tolerance Maximum allowed difference between values to consider them duplicates
 *                 (default: 1e-10)
 * @return std::vector<std::pair<size_t, size_t>> Pairs of indices where duplicates were found
 *
 * @note The tolerance parameter helps handle floating-point precision issues
 * @note Complexity: O(n log n) where n is the size of y
 */
std::vector<std::pair<size_t, size_t>> find_duplicate_values(
    const std::vector<double>& y,
    double tolerance = 1e-10
    ) {

    std::vector<std::pair<size_t, size_t>> duplicates;

    // Create pairs of (value, index) for sorting
    std::vector<std::pair<double, size_t>> indexed_values;
    indexed_values.reserve(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        indexed_values.emplace_back(y[i], i);
    }

    // Sort by value
    std::sort(indexed_values.begin(), indexed_values.end());

    // Check consecutive elements for duplicates
    for (size_t i = 1; i < indexed_values.size(); ++i) {
        double diff = std::abs(indexed_values[i].first - indexed_values[i-1].first);
        if (diff <= tolerance) {
            duplicates.emplace_back(
                indexed_values[i-1].second,
                indexed_values[i].second
            );
        }
    }

    return duplicates;
}


/**
 * @brief Validation function that throws an exception if duplicates are found
 *
 * @param y Vector of values to validate
 * @param tolerance Maximum allowed difference between values
 * @throws std::invalid_argument if duplicates are found
 */
void validate_no_duplicates(const std::vector<double>& y, double tolerance = 1e-10) {
    auto duplicates = find_duplicate_values(y, tolerance);
    if (!duplicates.empty()) {
        std::ostringstream oss;
        oss << "Duplicate values found in y at indices: ";
        for (const auto& [i, j] : duplicates) {
            oss << "(" << i << "," << j << ": " << y[i] << "," << y[j] << ") ";
        }
        REPORT_ERROR(oss.str().c_str());
    }
}

/**
 * @brief Adds small random noise to break ties in case of Rf_duplicate values
 *
 * @param y Vector of values to modify
 * @param noise_magnitude Maximum magnitude of noise to add (default: 1e-10)
 * @return true if modifications were made, false if no duplicates found
 */
bool break_duplicate_values(std::vector<double>& y, double noise_magnitude) {
    auto duplicates = find_duplicate_values(y);
    if (duplicates.empty()) {
        return false;
    }

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-noise_magnitude, noise_magnitude);

    // Track indices that need noise
    std::unordered_set<size_t> indices_to_modify;
    for (const auto& [i, j] : duplicates) {
        indices_to_modify.insert(i);
        indices_to_modify.insert(j);
    }

    // Add noise to all duplicates
    for (size_t idx : indices_to_modify) {
        y[idx] += dis(gen);
    }

    return true;
}

/**
 * @brief Adds small random noise only to non-unique values in a vector
 *
 * This is especially useful for extrema detection algorithms that require
 * distinct values to properly identify local minima and maxima
 *
 * @param y Vector of values to selectively apply jitter to
 * @return A copy of the input vector with jitter applied only to duplicates
 */
std::vector<double> selective_jitter(const std::vector<double>& y) {
    // Return empty vector if input is empty
    if (y.empty()) {
        return {};
    }

    // Make a copy of the input vector
    std::vector<double> result = y;

    // Calculate appropriate jitter magnitude (z/5 where z is the smallest difference)
    double jitter_magnitude = 0.0;
    {
        // Make a copy for sorting
        std::vector<double> sorted = y;
        std::sort(sorted.begin(), sorted.end());

        // Find smallest non-zero difference between consecutive values
        double min_diff = std::numeric_limits<double>::max();
        for (size_t i = 1; i < sorted.size(); ++i) {
            double diff = sorted[i] - sorted[i-1];
            if (diff > 0 && diff < min_diff) {
                min_diff = diff;
            }
        }

        // If all values are the same or there's only one value
        if (min_diff == std::numeric_limits<double>::max()) {
            // Use a small fraction of the value itself (or 1e-10 if value is 0)
            min_diff = std::abs(sorted[0]) > 0 ? std::abs(sorted[0]) * 1e-5 : 1e-10;
        }

        // Calculate jitter magnitude as z/5
        jitter_magnitude = min_diff / 5.0;
    }

    // Use the break_duplicate_values function with calculated jitter magnitude
    break_duplicate_values(result, jitter_magnitude);

    return result;
}


/**
 * @brief Adds small random noise to values similar to R's jitter function
 *
 * @param y Vector of values to add jitter to
 * @return A copy of the input vector with jitter applied
 */
std::vector<double> jitter(const std::vector<double>& y) {
    // Return empty vector if input is empty
    if (y.empty()) {
        return {};
    }

    // Make a copy of the input vector for sorting
    std::vector<double> sorted = y;
    std::sort(sorted.begin(), sorted.end());

    // Find smallest non-zero difference between consecutive values
    double min_diff = std::numeric_limits<double>::max();
    for (size_t i = 1; i < sorted.size(); ++i) {
        double diff = sorted[i] - sorted[i-1];
        if (diff > 0 && diff < min_diff) {
            min_diff = diff;
        }
    }

    // If all values are the same or there's only one value
    if (min_diff == std::numeric_limits<double>::max()) {
        // Use a small fraction of the value itself (or 1e-10 if value is 0)
        min_diff = std::abs(sorted[0]) > 0 ? std::abs(sorted[0]) * 1e-5 : 1e-10;
    }

    // Calculate jitter magnitude as z/5
    double jitter_magnitude = min_diff / 5.0;

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-jitter_magnitude, jitter_magnitude);

    // Apply jitter to a copy of the input vector
    std::vector<double> result = y;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] += dis(gen);
    }

    return result;
}

/**
 * @brief Conditionally adds small random noise only to non-unique values
 *
 * This function checks if duplicates exist, and only applies jitter when
 * necessary. It's optimized for extrema detection algorithms that require
 * distinct values to properly identify local minima and maxima.
 *
 * @param y Vector of values to potentially modify (modified in-place)
 * @param tolerance Maximum allowed difference between values to consider them duplicates
 *
 * @return true if duplicates were found and modifications were made, false otherwise
 */
bool conditional_selective_jitter(
    std::vector<double>& y,
    double tolerance = 1e-10
    ) {
    // Return false if input is empty
    if (y.empty()) {
        return false;
    }

    // First check if there are any duplicates
    auto duplicates = find_duplicate_values(y, tolerance);
    if (duplicates.empty()) {
        // No duplicates found, no need to apply jitter

        // debugging
        // Rprintf("In conditional_selective_jitter(): No duplicates found, no need to apply jitter\n");

        return false;
    }

    // Calculate appropriate jitter magnitude (z/5 where z is the smallest difference)
    double jitter_magnitude = 0.0;
    {
        // Make a copy for sorting
        std::vector<double> sorted = y;
        std::sort(sorted.begin(), sorted.end());

        // Find smallest non-zero difference between consecutive values
        double min_diff = std::numeric_limits<double>::max();
        for (size_t i = 1; i < sorted.size(); ++i) {
            double diff = sorted[i] - sorted[i-1];
            if (diff > 0 && diff < min_diff) {
                min_diff = diff;
            }
        }

        // Rprintf("In conditional_selective_jitter(): min_diff: %.8f\n", min_diff);


        // If all values are the same or there's only one value
        if (min_diff == std::numeric_limits<double>::max()) {
            // Use a small fraction of the value itself (or 1e-10 if value is 0)
            min_diff = std::abs(sorted[0]) > 0 ? std::abs(sorted[0]) * 1e-5 : 1e-10;
        }

        // Calculate jitter magnitude as z/5
        jitter_magnitude = min_diff / 5.0;
    }

    // debugging
    // Rprintf("In conditional_selective_jitter(): jitter_magnitude: %.8f\n", jitter_magnitude);


    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-jitter_magnitude, jitter_magnitude);

    // Track indices that need noise
    std::unordered_set<size_t> indices_to_modify;
    for (const auto& [i, j] : duplicates) {
        indices_to_modify.insert(i);
        indices_to_modify.insert(j);
    }

    // debugging
    // print_uset(indices_to_modify, "indices_to_modify");

    // Add noise to all duplicates
    for (size_t idx : indices_to_modify) {
        y[idx] += dis(gen);
    }

    // After applying jitter, verify that we've actually resolved all duplicates
    // In rare cases, the jitter might not be enough to create uniqueness
    auto remaining_duplicates = find_duplicate_values(y, tolerance);
    if (!remaining_duplicates.empty()) {
        // If duplicates still exist, apply a slightly larger jitter to those values
        std::unordered_set<size_t> indices_to_fix;
        for (const auto& [i, j] : remaining_duplicates) {
            indices_to_fix.insert(i);
            indices_to_fix.insert(j);
        }

        // Use double the jitter magnitude for stubborn duplicates
        std::uniform_real_distribution<> stronger_dis(-jitter_magnitude * 2.0, jitter_magnitude * 2.0);
        for (size_t idx : indices_to_fix) {
            y[idx] += stronger_dis(gen);
        }
    }

    return true;
}

/**
 * @brief Constructs a gradient trajectory from a starting vertex, handling local extrema.
 *
 * This function builds a trajectory that follows the gradient of a function (y) over the graph,
 * either in ascending or descending direction. The trajectory is constructed step by step,
 * finding the path with optimal rate of change at each iteration. It handles local extrema by:
 *
 * 1. Checking if any extrema exist in the constructed path segment
 * 2. Truncating the path at the first encountered extremum when appropriate
 * 3. Potentially switching direction when encountering an opposing extremum
 *    (e.g., a minimum during ascent or a maximum during descent)
 *
 * The function uses a local scale parameter at each vertex to determine the search radius
 * for finding the next segment of the trajectory. This radius represents an approximate
 * region where the function y is well-approximated by a linear model.
 *
 * @param start The vertex from which to start the trajectory construction
 * @param ascending_init Initial direction flag - true for ascending trajectory (following increasing y values),
 *                        false for descending trajectory (following decreasing y values)
 * @param scale Vector of local scale values for each vertex (search radius for gradient computation)
 * @param y Vector of function values at each vertex (the function whose gradient we're following)
 * @param extrema_map Map of vertices that are local extrema of y, with flags indicating
 *                      whether each is a maximum (true) or minimum (false)
 *
 * @return A gradient_trajectory_t structure containing:
 *         - path: The sequence of vertices forming the gradient trajectory
 *         - ascending: The final direction of the trajectory, which may differ from ascending_init
 *                     if the trajectory encountered an opposing extremum
 *
 * @note The trajectory may be truncated if it encounters a "terminal" extremum
 *       (maximum during ascent or minimum during descent). The function will also
 *       automatically expand the search radius if no suitable path endpoints are found
 *       within the initial scale radius.
 *
 * @see compute_graph_reachability_map
 * @see get_vertex_shortest_paths
 * @see compute_gradient_flow
 */
gradient_trajectory_t set_wgraph_t::construct_trajectory(
    size_t start,
    bool ascending,
    const std::vector<double>& scale,
    const std::vector<double>& y,
    const std::unordered_map<size_t, bool>& extrema_map,
    double long_edge_lower_thld,
    double long_edge_upper_thld
    ) const {

#define DEBUG__construct_trajectory 1

    // Initialize result structure
    gradient_trajectory_t result;
    result.path.push_back(start);
    result.ends_at_critical = false;
    result.ends_at_lmax = false;

#if DEBUG__construct_trajectory
    Rprintf("Starting trajectory construction from vertex %zu (ascending = %d)\n", start, ascending);
    Rprintf("Value at start vertex: %f\n", y[start]);
#endif

    size_t current = start;
    bool should_continue = true;

    // Track vertices in the trajectory to prevent loops
    std::unordered_set<size_t> trajectory_vertices;
    trajectory_vertices.insert(start);


#if DEBUG__construct_trajectory
    // debugging/testing - saving to file
    std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data/";
    if (!std::filesystem::exists(debug_dir)) {
        if (!std::filesystem::create_directories(debug_dir)) {
            REPORT_ERROR("ERROR: Failed to create debug directory: %s\n", debug_dir.c_str());
        }
    }
#endif

    // Main trajectory building loop
    while (should_continue) {

#if DEBUG__construct_trajectory
        Rprintf("Current vertex: %zu, Value: %f, Scale: %f\n", current, y[current], scale[current]);
#endif

        // Compute reachability map from current vertex
        reachability_map_t reachability_map = compute_graph_reachability_map(current, scale[current]);
        std::vector<vertex_path_t> paths = reconstruct_graph_paths(reachability_map);

#if DEBUG__construct_trajectory
        Rprintf("Number of paths found: %zu\n", paths.size());
#endif

        // If no paths found, gradually increase radius until we find some
        if (paths.empty()) {
            double multiplier = 2.0;
            double radius = scale[current];
            while (paths.empty()) {
                radius *= multiplier;
                reachability_map = compute_graph_reachability_map(current, radius);
                paths = reconstruct_graph_paths(reachability_map);

                // Safety check to avoid infinite loops
                if (radius > 1000 * scale[current]) {
                    REPORT_WARNING("Failed to find paths even with greatly enlarged radius");
                    return result;
                }
            }
        }

        // Filter out paths that would create loops in the trajectory
        std::vector<vertex_path_t> filtered_paths;
        for (const auto& path : paths) {
            bool creates_loop = false;
            for (size_t i = 1; i < path.vertices.size(); ++i) {
                // Check if this vertex (other than the current one) is already in our trajectory
                if (trajectory_vertices.count(path.vertices[i]) > 0) {
                    creates_loop = true;
                    break;
                }
            }
            if (!creates_loop) {
                filtered_paths.push_back(path);
            }
        }
        paths = filtered_paths;

#if DEBUG__construct_trajectory
        Rprintf("Number of filtered paths found: %zu\n", paths.size());

        for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
            std::string path_file_name = debug_dir + "vertex_" + std::to_string(current + 1) + "_path_" + std::to_string(path_idx) + ".csv";
            std::ofstream path_file(path_file_name);
            if (!path_file.is_open()) {
                REPORT_ERROR("WARNING: Could not open file %s for writing\n", path_file_name.c_str());
            }
            const auto& path = paths[path_idx];
            // Write header
            path_file << "vertex,dist,value\n";
            for (size_t i = 0; i < path.vertices.size(); ++i) {
                size_t vertex = path.vertices[i];
                path_file << (vertex + 1) << "," << path.dist_to_ref_vertex[i] << "," << y[vertex] << "\n"; // 1-indexed vertex for R compatibility
            }
            path_file.close();
        }

        Rprintf("Vertex %zu paths written to directory: %s\n", current, debug_dir.c_str());
#endif


        // If all paths would create loops, we're stuck - terminate
        if (paths.empty()) {
            REPORT_WARNING("WARNING: in construct_trajectory() with start vertex: %zu all paths would create loops. Terminating search!!!", start);
            break;
        }

        // Store path evaluation results
        struct evaluated_path_t {
            size_t path_index;
            double quality_metric;
            double total_change;
            size_t best_vertex_index;
            bool has_extremum;
            bool is_maximum;
            size_t extremum_index;
        };
        std::vector<evaluated_path_t> evaluated_paths;

        // Evaluate each path
        for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {

#if DEBUG__construct_trajectory
            std::string path_file_name = debug_dir + "vertex_" + std::to_string(current + 1) + "_path_" + std::to_string(path_idx) + "_quality_metrics.csv";
            std::ofstream path_file(path_file_name);
            if (!path_file.is_open()) {
                REPORT_ERROR("WARNING: Could not open file %s for writing\n", path_file_name.c_str());
            }
            path_file << "vertex,tot_chg,MI,pf,path_dist,adj_rate,QM\n";
#endif

            const auto& path = paths[path_idx];
            evaluated_path_t Rf_eval;
            Rf_eval.path_index = path_idx;
            Rf_eval.has_extremum = false;
            Rf_eval.best_vertex_index = 0;
            Rf_eval.quality_metric = ascending ? -std::numeric_limits<double>::max() :
                                             std::numeric_limits<double>::max();

            // Check each vertex along the path
            for (size_t i = 1; i < path.vertices.size(); ++i) {
                size_t vertex = path.vertices[i];

                // Calculate total function change
                double total_change = y[vertex] - y[current];

                // Calculate cumulative absolute changes
                double cumulative_absolute_changes = 0.0;
                double prev_y = y[current];
                for (size_t j = 1; j <= i; ++j) {
                    size_t v = path.vertices[j];
                    cumulative_absolute_changes += std::abs(y[v] - prev_y);
                    prev_y = y[v];
                }

                // Calculate monotonicity index
                double mono_index = (cumulative_absolute_changes > 0.0) ?
                    std::abs(total_change) / cumulative_absolute_changes : 1.0;

                // Extract edge lengths
                std::vector<double> edge_lengths;
                edge_lengths.reserve(i);
                for (size_t j = 1; j <= i; ++j) {
                    edge_lengths.push_back(path.dist_to_ref_vertex[j] - path.dist_to_ref_vertex[j-1]);
                }

                // Calculate maximum thresholded edge length
                double max_thldd_edge_length = calculate_path_max_threshold(
                    edge_lengths, long_edge_lower_thld, long_edge_upper_thld);

                // Calculate penalty factor based on edge length
                double penalty_factor = 1.0 + max_thldd_edge_length; // Ranges from 1.0 to 2.0

                // Calculate adjusted rate
                double path_distance = path.dist_to_ref_vertex[i];
                double adjusted_rate = 0.0;
                if (path_distance > 0.0) {
                    adjusted_rate = total_change / std::pow(path_distance, penalty_factor);
                }

                // Calculate quality metric
                double quality_metric = mono_index * adjusted_rate;

#if DEBUG__construct_trajectory
                path_file << (vertex + 1) << ","
                          << total_change << ","
                          << mono_index << ","
                          << penalty_factor << ","
                          << path_distance  << ","
                          << adjusted_rate  << ","
                          << quality_metric << "\n";
#endif
                // Check if vertex is an extremum
                auto extremum_it = extrema_map.find(vertex);
                if (extremum_it != extrema_map.end()) {
                    Rf_eval.has_extremum = true;
                    Rf_eval.is_maximum = extremum_it->second;
                    Rf_eval.extremum_index = i;
                    Rf_eval.quality_metric = quality_metric;
                    Rf_eval.total_change = total_change;
                    break;
                }

                // Determine if this is the best vertex so far based on direction
                bool is_better = false;
                if (ascending) {
                    is_better = (quality_metric > 0 && quality_metric > Rf_eval.quality_metric);
                } else {
                    is_better = (quality_metric < 0 && quality_metric < Rf_eval.quality_metric);
                }

                if (is_better) {
                    Rf_eval.quality_metric = quality_metric;
                    Rf_eval.total_change = total_change;
                    Rf_eval.best_vertex_index = i;
                }
            } // END OF for (size_t i = 1; i < path.vertices.size(); ++i)

#if DEBUG__construct_trajectory
            path_file.close();

            bool cond = Rf_eval.has_extremum ||
                (ascending && Rf_eval.quality_metric > 0) ||  // Only positive quality for ascending
                (!ascending && Rf_eval.quality_metric < 0);   // Only negative quality for descending

            Rprintf("Rf_eval.has_extremum: %d\nascending: %d\nRf_eval.quality_metric: %.4f\ncond: %d\n",
                    (int)Rf_eval.has_extremum, (int)ascending, Rf_eval.quality_metric, (int)cond);
#endif

            // Only include paths that have valid metrics
            if (Rf_eval.has_extremum ||
                (ascending && Rf_eval.quality_metric > 0) ||  // Only positive quality for ascending
                (!ascending && Rf_eval.quality_metric < 0)    // Only negative quality for descending
                ) {
                evaluated_paths.push_back(Rf_eval);
            }
        } // END OF for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx)

#if DEBUG__construct_trajectory
        Rprintf("Found %zu valid paths from vertex %zu\n", evaluated_paths.size(), current);
        Rf_error("DEBUGGING");
#endif

        // If we have no valid paths, terminate
        if (evaluated_paths.empty()) {
            REPORT_WARNING("WARNING: in construct_trajectory() with start vertex: %zu evaluated_paths is empty. Breaking search!!!", start);
            break;
        }

        // Find the best path based on extrema and quality metrics
        size_t best_path_idx = 0;
        bool found_best = false;

        // Find the best path with an appropriate extremum
        double best_extremum_quality = ascending ? -std::numeric_limits<double>::max() :
            std::numeric_limits<double>::max();
        bool found_appropriate_extremum = false;

        for (size_t i = 0; i < evaluated_paths.size(); ++i) {
            const auto& Rf_eval = evaluated_paths[i];
            if (Rf_eval.has_extremum &&
                ((ascending && Rf_eval.is_maximum) || (!ascending && !Rf_eval.is_maximum))) {
                // For ascending, select highest quality metric among appropriate extrema
                // For descending, select lowest quality metric among appropriate extrema
                bool is_better_extremum = ascending ?
                    (Rf_eval.quality_metric > best_extremum_quality) :
                    (Rf_eval.quality_metric < best_extremum_quality);

                if (!found_appropriate_extremum || is_better_extremum) {
                    found_appropriate_extremum = true;
                    best_extremum_quality = Rf_eval.quality_metric;
                    best_path_idx = i;
                }
            }
        }

        if (found_appropriate_extremum) {
            found_best = true;

            // Process the selected path with an extremum...
            const auto& Rf_eval = evaluated_paths[best_path_idx];

            // Mark that we're ending at a critical point
            result.quality_metric = Rf_eval.quality_metric;
            result.total_change = Rf_eval.total_change;
            result.ends_at_critical = true;
            result.ends_at_lmax = Rf_eval.is_maximum;

            // Add path up to and including the extremum
            const auto& best_path = paths[Rf_eval.path_index];
            result.path.insert(
                result.path.end(),
                best_path.vertices.begin() + 1,
                best_path.vertices.begin() + Rf_eval.extremum_index + 1
                );

            // We found a critical point, so terminate
            should_continue = false;
        }

        // If no appropriate extremum found, use quality metrics
        if (!found_best) {
            // Find path with best quality metric
            if (ascending) {
                double best_quality = -std::numeric_limits<double>::max();
                for (size_t i = 0; i < evaluated_paths.size(); ++i) {
                    if (evaluated_paths[i].quality_metric > best_quality) {
                        best_quality = evaluated_paths[i].quality_metric;
                        best_path_idx = i;
                        found_best = true;
                    }
                }
            } else {
                double best_quality = std::numeric_limits<double>::max();
                for (size_t i = 0; i < evaluated_paths.size(); ++i) {
                    if (evaluated_paths[i].quality_metric < best_quality) {
                        best_quality = evaluated_paths[i].quality_metric;
                        best_path_idx = i;
                        found_best = true;
                    }
                }
            }

            if (found_best) {
                const auto& Rf_eval = evaluated_paths[best_path_idx];
                result.quality_metric = Rf_eval.quality_metric;
                result.total_change = Rf_eval.total_change;

                const auto& best_path = paths[Rf_eval.path_index];

                // Determine the appropriate index up to which we should add vertices
                size_t end_index = Rf_eval.has_extremum ? Rf_eval.extremum_index : Rf_eval.best_vertex_index;

                // Add vertices to our trajectory set
                for (size_t j = 1; j <= end_index; ++j) {
                    trajectory_vertices.insert(best_path.vertices[j]);
                }

                // Add vertices to result path
                result.path.insert(
                    result.path.end(),
                    best_path.vertices.begin() + 1,
                    best_path.vertices.begin() + end_index + 1);

                // Update current vertex to the end of the selected path segment
                current = best_path.vertices[end_index];

                // If we found an extremum, mark that we're ending at a critical point
                if (Rf_eval.has_extremum) {
                    result.ends_at_critical = true;
                    result.ends_at_lmax = Rf_eval.is_maximum;
                    should_continue = false;  // Terminate the trajectory
                }
            }
        }

        // If we couldn't find a best path, terminate
        if (!found_best) {
            REPORT_WARNING("WARNING: in construct_trajectory() with start vertex: %zu we couldn't find a best path. Terminating search!!!", start);
            break;
        }
    }

    return result;
}


/**
 * @brief Computes the gradient flow structure including gradient flow cells
 *
 * This method performs a comprehensive analysis of the function y defined over
 * the graph vertices, computing gradient trajectories and decomposing the graph
 * into gradient flow cells. The computation proceeds in several stages:
 *
 * 1. Local Extrema Detection:
 *    - Identifies tau-local minima and maxima using the radius parameter
 *    - Breaks ties in Rf_duplicate values if necessary
 *
 * 2. Gradient Trajectory Computation:
 *    - For each unprocessed vertex, constructs ascending and descending trajectories
 *    - Combines trajectories to form complete paths between local minimum-maximum pairs
 *    - Each trajectory follows paths of steepest ascent/descent within tau-neighborhoods
 *
 * 3. Gradient Flow complex Construction:
 *    - Groups trajectories by their extremal pairs to form basins and pro-cells
 *
 * The resulting gradient_flow_t structure contains:
 * - trajectories: All gradient trajectories in the graph
 * - local_extrema: Pairs of (vertex index, is_maximum) for all local extrema
 * - cells: Gradient flow cells, each defined by a minimum-maximum pair
 *
 *
 * @param y Vector of function values defined on graph vertices
 * @param scale Vertex radii used for local extrema detection and gradient flow computation
 * @param quantile_scale_thld Scale values are truncated to the interval (0, scale_quantile), where scale_quantile = quantile(scale, prob = quantile_scale_thld);

 * @return gradient_flow_t Structure containing the complete gradient flow complex information
 *
 * @throws std::runtime_error If y's size doesn't match the number of vertices
 *
 * @Rf_warning The function modifies input vector y if Rf_duplicate values are found,
 *          adding small noise to break ties
 *
 * @see gradient_flow_t For the complete structure definition
 * @see find_connected_components() For the cell decomposition algorithm
 * @see check_local_extremum() For local extrema detection details
 * @see construct_trajectory() For gradient trajectory construction details
 */
gradient_flow_t set_wgraph_t::compute_gradient_flow(
    std::vector<double>& y,
    std::vector<double>& scale,
    double quantile_scale_thld
    ) const {

    #define DEBUG__compute_gradient_flow 0

    gradient_flow_t result;

    // Check if there are Rf_duplicate values, and if there are, add to them random
    // noise runif(min = -z/5, max = z/5), where z is the smallest difference of
    // sorted y values
    double tolerance = 1e-10;
    if (conditional_selective_jitter(y, tolerance)) {
        result.messages = "Small noise added to break ties in Rf_duplicate values";
        REPORT_WARNING("Small noise added to break ties in Rf_duplicate values\n");
    }

    // Get 90th and 95th percentiles of edge weights those are unusually long edges that may cross basin boundaries
    std::vector<double> edge_weight_percentiles = compute_weight_percentiles({0.9, 0.95});
    double q90 = edge_weight_percentiles[0];
    double q95 = edge_weight_percentiles[1];

    // Clear cache and initialize unprocessed vertices
    paths_cache.clear();
    unprocessed_vertices.clear();
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        unprocessed_vertices.insert(i);
    }

    // debug
    // Rprintf("In compute_gradient_flow(): scale.size(): %zu\n", scale.size());

    // Finding small scale
    double small_scale = quantile(scale, quantile_scale_thld);

    // Truncate all scale values to at most small_scale
    for (auto& s : scale) {
        s = std::min(s, small_scale);
    }

    // Find local extrema
    for (size_t v : unprocessed_vertices) {
        auto shortest_paths = find_graph_paths_within_radius(v, scale[v]);
        auto extremum = check_local_extremum(v, shortest_paths, y);
        if (extremum) {
            result.local_extrema[extremum->first] = extremum->second;
        }
    }

    // debugging
    // print_umap(result.local_extrema, "result.local_extrema");
    // Rf_error("compute_gradient_flow() DEBUG");

#if 0
    // After computing local extrema
    Rprintf("Found %zu local extrema:\n", result.local_extrema.size());
    for (const auto& [vertex, is_maximum] : result.local_extrema) {
        Rprintf("  Vertex %zu: %s, Value: %f\n",
                vertex, is_maximum ? "Maximum" : "Minimum", y[vertex]);
    }

    // After scale calculation
    Rprintf("Scale values used: min=%f, max=%f, quantile=%f\n",
            *std::min_element(scale.begin(), scale.end()),
            *std::max_element(scale.begin(), scale.end()),
            small_scale);
#endif

    while (!unprocessed_vertices.empty()) {
        size_t v = *unprocessed_vertices.begin();

        auto ascending_traj = construct_trajectory(
            v,
            true,
            scale,
            y,
            result.local_extrema,
            q90,
            q95
            );

        auto descending_traj = construct_trajectory(
            v,
            false,
            scale,
            y,
            result.local_extrema,
            q90,
            q95
            );

        // Combine trajectories as before
        gradient_flow_t::trajectory_t trajectory;
        trajectory.vertices.insert(
            trajectory.vertices.end(),
            descending_traj.path.rbegin(),
            descending_traj.path.rend() - 1  // Skip the starting vertex to avoid duplication
        );
        trajectory.vertices.insert(
            trajectory.vertices.end(),
            ascending_traj.path.begin(),
            ascending_traj.path.end()
        );

        #if DEBUG__compute_gradient_flow
        // debuggging
        print_vect(ascending_traj.path, "ascending_traj.path", 1);
        print_vect(descending_traj.path, "descending_traj.path", 1);
        print_vect(trajectory.vertices, "trajectory.vertices", 1);
        #endif

        // Determine trajectory type using the same logic but with the new trajectory objects
        if ((descending_traj.ends_at_critical && !descending_traj.ends_at_lmax) &&
            (ascending_traj.ends_at_critical && ascending_traj.ends_at_lmax)) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_LMAX;
        } else if (
            (descending_traj.ends_at_critical && !descending_traj.ends_at_lmax && !ascending_traj.ends_at_critical) ||
            (ascending_traj.ends_at_critical && !ascending_traj.ends_at_lmax && !descending_traj.ends_at_critical)
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_ONLY;
        } else if (
            (descending_traj.ends_at_critical && descending_traj.ends_at_lmax && !ascending_traj.ends_at_critical) ||
            (ascending_traj.ends_at_critical && ascending_traj.ends_at_lmax && !descending_traj.ends_at_critical)
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMAX_ONLY;
        } else if (
            descending_traj.ends_at_critical && !descending_traj.ends_at_lmax &&
            ascending_traj.ends_at_critical && !ascending_traj.ends_at_lmax
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_LMIN;
        } else if (
            descending_traj.ends_at_critical && descending_traj.ends_at_lmax &&
            ascending_traj.ends_at_critical && ascending_traj.ends_at_lmax
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMAX_LMAX;
        } else {
            trajectory.trajectory_type = gradient_flow_t::UNKNOWN;
        }

        // Transfer quality metrics
        // Average of both trajectories <<--- we need to recompute this
        trajectory.quality_metric = (ascending_traj.quality_metric + descending_traj.quality_metric) / 2.0;
        // Sum of absolute changes <<--- we need to recompute this
        trajectory.total_change = std::abs(ascending_traj.total_change) + std::abs(descending_traj.total_change);

        result.trajectories.push_back(trajectory);

        // Remove processed vertices
        for (size_t vertex : trajectory.vertices) {
            unprocessed_vertices.erase(vertex);
        }

        #if DEBUG__compute_gradient_flow
        //REPORT_ERROR("DEBUGGING compute_gradient_flow\n");
        #endif
    }

    // Computing ascending and descending basins
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];

        // Handle each trajectory based on its type
        switch (traj.trajectory_type) {
        case gradient_flow_t::LMIN_LMAX:
            // Add all vertices to both basins using range insertion
            result.ascending_basin_map[traj.vertices.front()].insert(
                traj.vertices.begin(), traj.vertices.end());
            result.descending_basin_map[traj.vertices.back()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::LMIN_ONLY:
            // Only min is a critical point
            result.ascending_basin_map[traj.vertices.front()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::LMAX_ONLY:
            // Only max is a critical point
            result.descending_basin_map[traj.vertices.back()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::UNKNOWN:
            REPORT_WARNING("Found trajectory of UNKNOWN type");
            break;

        default:
            REPORT_WARNING("Found trajectory with no type assignment");
            break;
        }
    }

    // Create map of (min, max) pairs to the set = the union of the correspoing trajectory vertices
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_LMAX) {
            result.cell_map[{traj.vertices.front(), traj.vertices.back()}].insert(
                traj.vertices.begin(), traj.vertices.end());
        }
    }
    // After all proper pro-cells are populated we add to them corresponding dangling trajectories
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_ONLY) {
            // How to insert traj.vertices into all pro-cells that have traj.vertices.front() as the local minimum?
            // loop over all pro-cell with traj.vertices.front() as the lmin
            // cell_map[{traj.vertices.front(), ???}].insert(
            //     traj.vertices.begin(), traj.vertices.end());
        }
    }

    // Create map of (min, max) pairs to the set = the union of the correspoing trajectory vertices
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_LMAX) {
            result.cell_map[{traj.vertices.front(), traj.vertices.back()}].insert(
                traj.vertices.begin(), traj.vertices.end());
        }
    }

    // After all proper cells are populated we add to them corresponding dangling trajectories
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];

        if (traj.trajectory_type == gradient_flow_t::LMIN_ONLY) {
            // Add vertices to all pro-cells that have the same local minimum
            size_t lmin = traj.vertices.front();

            for (auto& [key, vertices] : result.cell_map) {
                if (key.first == lmin) {
                    // This pro-cell has the same local minimum
                    vertices.insert(traj.vertices.begin(), traj.vertices.end());
                }
            }
        }
        else if (traj.trajectory_type == gradient_flow_t::LMAX_ONLY) {
            // Add vertices to all pro-cells that have the same local maximum
            size_t lmax = traj.vertices.back();

            for (auto& [key, vertices] : result.cell_map) {
                if (key.second == lmax) {
                    // This pro-cell has the same local maximum
                    vertices.insert(traj.vertices.begin(), traj.vertices.end());
                }
            }
        }
    }

    result.scale = std::move(scale);

    return result;
}

/**
 * @brief Standalone function to construct gradient flow for a graph
 *
 * Wrapper function that constructs the gradient flow structure for a given graph
 * using specified vertex values and bandwidth (radius).
 *
 * @param graph The weighted graph to analyze
 * @param y Vector of values associated with vertices
 * @param bandwidth Maximum distance to consider for neighboring vertices
 *
 * @return gradient_flow_t The computed gradient flow structure
 *
 * @pre y.size() must match the number of vertices in graph
 * @pre bandwidth must be positive
 * @pre y must not contain Rf_duplicate values
 *
 * @see set_wgraph_t::compute_gradient_flow
 */
gradient_flow_t construct_graph_gradient_flow(
    const set_wgraph_t& graph,
    std::vector<double>& y,
    double bandwidth,
    double quantile_scale_thld
    ) {

    std::vector<double> scale(y.size(), bandwidth);

    return graph.compute_gradient_flow(y, scale, quantile_scale_thld);
}


/**
 * @brief Computes gradient flow structure for a function over a weighted graph (R interface)
 *
 * This function provides an R interface to compute gradient flow trajectories, basins, and
 * pro-cells for a given weighted graph and function over its vertices. It identifies local extrema
 * (minima and maxima) and constructs trajectories that follow the gradient of the function.
 *
 * @param s_adj_list SEXP (list) Adjacency lists for each vertex where s_adj_list[[i]]
 *                   contains indices of vertices adjacent to vertex i (1-based indices in R)
 * @param s_weight_list SEXP (list) Edge weights where s_weight_list[[i]][j] is the weight
 *                      of edge from vertex i to s_adj_list[[i]][j]
 * @param s_y SEXP (numeric) Vector of function values at graph vertices
 * @param s_scale SEXP (numeric) Vector of scale values for each vertex, controlling the
 *                local neighborhood size for gradient computation and extrema detection
 * @param s_with_trajectories SEXP (logical) Whether to include trajectory information in output
 *
 * @return SEXP (list) with components:
 *   - local_extrema: (matrix) Nx2 matrix where N is number of extrema. Each row contains
 *                    (vertex_index, is_maximum)
 *   - trajectories: (list or NULL) If s_with_trajectories is TRUE, list of trajectory data where
 *                  each element contains:
 *                    - vertices: Integer vector of vertex indices
 *                    - type: Type of trajectory (LMIN_LMAX, LMIN_ONLY, LMAX_ONLY, or UNKNOWN)
 *   - basins: (list) Contains ascending and descending basins:
 *             - ascending: List of basins around local minima
 *               - local_min: Index of local minimum vertex
 *               - vertices: Vertices in the ascending basin
 *             - descending: List of basins around local maxima
 *               - local_max: Index of local maximum vertex
 *               - vertices: Vertices in the descending basin
 *   - cells: (list) gradient flow cells, each containing:
 *               - local_min: Index of local minimum vertex
 *               - local_max: Index of local maximum vertex
 *               - vertices: Vertices in the pro-cell
 *   - messages: (character) Any messages or warnings generated during computation
 *
 * @note All vertex indices in the return value are 1-based (R convention)
 * @note Trajectories are included in output only if s_with_trajectories is TRUE
 *
 * @see gradient_flow_t For the underlying C++ structure
 * @see compute_gradient_flow For the core computation function
 */
SEXP S_construct_graph_gradient_flow(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_scale,
    SEXP s_quantile_scale_thld,
    SEXP s_with_trajectories
    ) {
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    size_t n_vertices = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_vertices);
    std::vector<double> scale(REAL(s_scale), REAL(s_scale) + n_vertices);
    double quantile_scale_thld = REAL(s_quantile_scale_thld)[0];
    bool with_trajectories = LOGICAL(s_with_trajectories)[0];

    //print_vect(scale, "scale in S_construct_graph_gradient_flow");

    set_wgraph_t graph(adj_list, weight_list);
    gradient_flow_t result = graph.compute_gradient_flow(
        y,
        scale,
        quantile_scale_thld
        );

    // Create names for the list elements
    const char* names[] = {
        "local_extrema",
        "trajectories",
        "basins",
        "cells",
        "messages",
        NULL};

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Convert results to R list
    size_t n_protected = 0;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements)); n_protected++;
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    for (size_t i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    UNPROTECT(1); // for r_names

    // Set local extrema
    size_t n_extrema = result.local_extrema.size();
    SEXP r_extrema = PROTECT(Rf_allocMatrix(INTSXP, n_extrema, 2));
    int* extrema_ptr = INTEGER(r_extrema);
    size_t i = 0;
    for (const auto& [vertex, is_maximum] : result.local_extrema) {
        extrema_ptr[i] = vertex + 1;  // Convert to 1-based indexing
        extrema_ptr[i + n_extrema] = is_maximum;  // is_maximum flag
        i++;
    }
    SET_VECTOR_ELT(r_result, 0, r_extrema);
    UNPROTECT(1);

    // Set trajectories
    if (with_trajectories) {
        size_t n_trajectories = result.trajectories.size();
        SEXP r_trajectories = PROTECT(Rf_allocVector(VECSXP, n_trajectories));
        for (size_t i = 0; i < n_trajectories; i++) {
            const auto& traj = result.trajectories[i];

            // Vertices
            SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, traj.vertices.size()));
            int* traj_ptr = INTEGER(r_vertices);
            for (size_t j = 0; j < traj.vertices.size(); j++) {
                traj_ptr[j] = traj.vertices[j] + 1;  // Convert to 1-based indexing
            }

            #if 0
            if (i == 0) {
                Rprintf("\ni: %zu\n",i);
                print_vect(traj.vertices, "traj.vertices", 1);
                Rf_error("S_construct_graph_gradient_flow(): DEBUGGING\n");
            }
            #endif

            // Trajectory type as string
            SEXP r_type = PROTECT(Rf_allocVector(STRSXP, 1));
            const char* type_str;
            switch (traj.trajectory_type) {
            case gradient_flow_t::LMIN_LMAX:
                type_str = "LMIN_LMAX";
                break;
            case gradient_flow_t::LMIN_ONLY:
                type_str = "LMIN_ONLY";
                break;
            case gradient_flow_t::LMAX_ONLY:
                type_str = "LMAX_ONLY";
                break;
            default:
                type_str = "UNKNOWN";
            }
            SET_STRING_ELT(r_type, 0, Rf_mkChar(type_str));

            // Create a list with vertices and trajectory type
            SEXP r_traj = PROTECT(Rf_allocVector(VECSXP, 4));

            SET_VECTOR_ELT(r_traj, 0, r_vertices);
            SET_VECTOR_ELT(r_traj, 1, r_type);
            UNPROTECT(2); // for r_vertices and r_type

            // Add quality metric
            SEXP r_quality = PROTECT(Rf_allocVector(REALSXP, 1));
            REAL(r_quality)[0] = result.trajectories[i].quality_metric;
            SET_VECTOR_ELT(r_traj, 2, r_quality);
            UNPROTECT(1);

            // Add total change
            SEXP r_change = PROTECT(Rf_allocVector(REALSXP, 1));
            REAL(r_change)[0] = result.trajectories[i].total_change;
            SET_VECTOR_ELT(r_traj, 3, r_change);
            UNPROTECT(1);

            // Update names
            SEXP r_traj_names = PROTECT(Rf_allocVector(STRSXP, 4));
            SET_STRING_ELT(r_traj_names, 0, Rf_mkChar("vertices"));
            SET_STRING_ELT(r_traj_names, 1, Rf_mkChar("type"));
            SET_STRING_ELT(r_traj_names, 2, Rf_mkChar("quality_metric"));
            SET_STRING_ELT(r_traj_names, 3, Rf_mkChar("total_change"));
            Rf_setAttrib(r_traj, R_NamesSymbol, r_traj_names);
            UNPROTECT(1);

            SET_VECTOR_ELT(r_trajectories, i, r_traj);
            UNPROTECT(1); // for r_traj
        }
        SET_VECTOR_ELT(r_result, 1, r_trajectories);
        UNPROTECT(1); // for r_trajectories
    } else {
        SET_VECTOR_ELT(r_result, 1, R_NilValue);
    }

    // Set basins
    SEXP r_basins = PROTECT(Rf_allocVector(VECSXP, 2));

    // Ascending basins
    size_t n_ascending = result.ascending_basin_map.size();
    SEXP r_ascending = PROTECT(Rf_allocVector(VECSXP, n_ascending));
    i = 0;
    for (const auto& [min_vertex, basin] : result.ascending_basin_map) {
        // Create a list with local_min and vertices
        SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 2));

        // Local minimum
        SEXP r_min = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_min)[0] = min_vertex + 1;  // Convert to 1-based indexing

        // Vertices in the basin
        SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, basin.size()));
        int* vertices_ptr = INTEGER(r_vertices);
        size_t j = 0;
        for (const auto& vertex : basin) {
            vertices_ptr[j++] = vertex + 1;  // Convert to 1-based indexing
        }

        SET_VECTOR_ELT(r_basin, 0, r_min);
        SET_VECTOR_ELT(r_basin, 1, r_vertices);
        UNPROTECT(2);  // for r_min, r_vertices

        // Set names for basin components
        SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(r_basin_names, 0, Rf_mkChar("local_min"));
        SET_STRING_ELT(r_basin_names, 1, Rf_mkChar("vertices"));
        Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);
        UNPROTECT(1);  // for r_basin_names

        SET_VECTOR_ELT(r_ascending, i++, r_basin);
        UNPROTECT(1);  // for r_basin
    }

    // Descending basins
    size_t n_descending = result.descending_basin_map.size();
    SEXP r_descending = PROTECT(Rf_allocVector(VECSXP, n_descending));
    i = 0;
    for (const auto& [max_vertex, basin] : result.descending_basin_map) {
        // Create a list with local_max and vertices
        SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 2));

        // Local maximum
        SEXP r_max = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_max)[0] = max_vertex + 1;  // Convert to 1-based indexing

        // Vertices in the basin
        SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, basin.size()));
        int* vertices_ptr = INTEGER(r_vertices);
        size_t j = 0;
        for (const auto& vertex : basin) {
            vertices_ptr[j++] = vertex + 1;  // Convert to 1-based indexing
        }

        SET_VECTOR_ELT(r_basin, 0, r_max);
        SET_VECTOR_ELT(r_basin, 1, r_vertices);
        UNPROTECT(2);  // for r_max, r_vertices

        // Set names for basin components
        SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(r_basin_names, 0, Rf_mkChar("local_max"));
        SET_STRING_ELT(r_basin_names, 1, Rf_mkChar("vertices"));
        Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);
        UNPROTECT(1);  // for r_basin_names

        SET_VECTOR_ELT(r_descending, i++, r_basin);
        UNPROTECT(1);  // for r_basin
    }

    SET_VECTOR_ELT(r_basins, 0, r_ascending);
    SET_VECTOR_ELT(r_basins, 1, r_descending);
    UNPROTECT(2);  // for r_ascending, r_descending

    // Set names for basins components
    SEXP r_basins_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_basins_names, 0, Rf_mkChar("ascending"));
    SET_STRING_ELT(r_basins_names, 1, Rf_mkChar("descending"));
    Rf_setAttrib(r_basins, R_NamesSymbol, r_basins_names);
    UNPROTECT(1);  // for r_basins_names

    SET_VECTOR_ELT(r_result, 2, r_basins);
    UNPROTECT(1);  // for r_basins

    // Set cells
    size_t n_cells = result.cell_map.size();
    SEXP r_cells = PROTECT(Rf_allocVector(VECSXP, n_cells));
    i = 0;
    for (const auto& [key_pair, vertices] : result.cell_map) {
        // Create a list with local_min, local_max, and vertices
        SEXP r_cell = PROTECT(Rf_allocVector(VECSXP, 3));

        // Local minimum
        SEXP r_min = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_min)[0] = key_pair.first + 1;  // Convert to 1-based indexing

        // Local maximum
        SEXP r_max = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_max)[0] = key_pair.second + 1;  // Convert to 1-based indexing

        // Vertices in the cell
        SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, vertices.size()));
        int* vertices_ptr = INTEGER(r_vertices);
        size_t j = 0;
        for (const auto& vertex : vertices) {
            vertices_ptr[j++] = vertex + 1;  // Convert to 1-based indexing
        }

        SET_VECTOR_ELT(r_cell, 0, r_min);
        SET_VECTOR_ELT(r_cell, 1, r_max);
        SET_VECTOR_ELT(r_cell, 2, r_vertices);
        UNPROTECT(3);  // for r_min, r_max, r_vertices

        // Set names for cell components
        SEXP r_cell_names = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(r_cell_names, 0, Rf_mkChar("local_min"));
        SET_STRING_ELT(r_cell_names, 1, Rf_mkChar("local_max"));
        SET_STRING_ELT(r_cell_names, 2, Rf_mkChar("vertices"));
        Rf_setAttrib(r_cell, R_NamesSymbol, r_cell_names);
        UNPROTECT(1);  // for r_cell_names

        SET_VECTOR_ELT(r_cells, i++, r_cell);
        UNPROTECT(1);  // for r_cell
    }
    SET_VECTOR_ELT(r_result, 3, r_cells);
    UNPROTECT(1);  // for r_cells

    // Set messages
    SEXP r_messages = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(r_messages, 0, Rf_mkChar(result.messages.c_str()));
    SET_VECTOR_ELT(r_result, 4, r_messages);
    UNPROTECT(1);  // for r_messages

    UNPROTECT(n_protected);
    return r_result;
}

