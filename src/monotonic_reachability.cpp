#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

#include <vector>         // std::vector
#include <queue>          // std::priority_queue
#include <random>         // std::mt19937
#include <sstream>        // optional: debugging/logging
#include <unordered_set>  // used in other modules?
#include <functional>     // used in std::pair / comparator
#include <map>            // for monotonic_reachability_map_t::info
#include <utility>        // for std::pair
#include <cmath>          // std::abs
#include <algorithm>      // std::reverse, std::sort

#include "reachability_map.hpp"
#include "cpp_utils.hpp"
#include "cpp_stats_utils.hpp"
#include "set_wgraph.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp"
#include "gradient_flow.hpp"


/**
 * @brief Computes a monotonic reachability map optimizing for path monotonicity
 *
 * This function performs a modified Dijkstra's algorithm to find paths from a
 * reference vertex to all reachable vertices within a specified radius,
 * optimizing for the monotonicity index rather than shortest distance. The
 * monotonicity index measures how consistently the function changes in one
 * direction (monotonically increasing or decreasing) along a path.
 *
 * @param ref_vertex The reference vertex from which to compute monotonic paths
 * @param y Vector of function values at each vertex
 * @param radius The maximum distance to search from the reference vertex
 * @param ascending If true, prioritize paths with increasing y values; if false, decreasing
 *
 * @return A monotonic_reachability_map_t structure containing the computed information
 */
monotonic_reachability_map_t set_wgraph_t::compute_monotonic_reachability_map(
    size_t ref_vertex,
    const std::vector<double>& y,
    double radius,
    bool ascending
    ) const {

    monotonic_reachability_map_t result;
    result.ref_vertex = ref_vertex;

    // Validate input
    if (ref_vertex >= adjacency_list.size()) {
        REPORT_ERROR("Reference vertex %zu is out of bounds (graph has %zu vertices)",
                     ref_vertex, adjacency_list.size());
    }

    if (y.size() != adjacency_list.size()) {
        REPORT_ERROR("Function value vector y size (%zu) doesn't match graph size (%zu)",
                     y.size(), adjacency_list.size());
    }

    // Initialize data structures
    size_t n = adjacency_list.size();
    std::vector<bool> visited(n, false);

    // Custom comparator for priority queue based on monotonicity index
    // Since monotonicity index is always non-negative and we want to maximize it in both cases,
    // we use the same comparison regardless of direction
    auto comparator = [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
        return a.second < b.second; // We want max monotonicity index regardless of direction
    };

    // Priority queue: pair of (vertex, monotonicity_index)
    std::priority_queue<std::pair<size_t, double>,
                        std::vector<std::pair<size_t, double>>,
                        decltype(comparator)> pq(comparator);

    // Initialize reference vertex
    monotonic_path_info_t ref_info;
    ref_info.total_change = 0.0;
    ref_info.cum_abs_change = 0.0;
    ref_info.monotonicity_index = 1.0; // By definition
    ref_info.predecessor = INVALID_VERTEX;
    ref_info.distance = 0.0;

    result.info[ref_vertex] = ref_info;

    // Add reference vertex to priority queue
    pq.push({ref_vertex, ascending ? -1.0 : 1.0}); // Adjust sign based on direction

    // Process vertices in order of monotonicity index
    while (!pq.empty()) {
        // Get vertex with best monotonicity index
        size_t u = pq.top().first;
        pq.pop();

        // Skip if already visited
        if (visited[u]) {
            continue;
        }

        // Mark as visited
        visited[u] = true;

        // Skip if beyond radius
        if (result.info[u].distance > radius) {
            continue;
        }

        // Store vertex info for result
        if (u != ref_vertex) {
            result.sorted_vertices.push_back({u, result.info[u].monotonicity_index});
        }

        // Explore neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double edge_weight = edge.weight;

            // Skip if already visited
            if (visited[v]) {
                continue;
            }

            // Calculate new path metrics
            double new_distance = result.info[u].distance + edge_weight;

            // Skip if beyond radius
            if (new_distance > radius) {
                continue;
            }

            double y_change = y[v] - y[u];
            double new_total_change = result.info[u].total_change + y_change;
            double new_cum_abs_change = result.info[u].cum_abs_change + std::abs(y_change);

            // Calculate monotonicity index
            double new_mono_index;
            if (new_cum_abs_change > 0.0) {
                new_mono_index = std::abs(new_total_change) / new_cum_abs_change;
            } else {
                new_mono_index = 1.0; // Perfect monotonicity when no change
            }

            // Direction preferences are handled via final path selection
            // rather than during path discovery
            // We accept all paths and rely on the monotonicity index
            // and final filtering to find the best ones

            // Update path info
            monotonic_path_info_t path_info;
            path_info.total_change = new_total_change;
            path_info.cum_abs_change = new_cum_abs_change;
            path_info.monotonicity_index = new_mono_index;
            path_info.predecessor = u;
            path_info.distance = new_distance;

            result.info[v] = path_info;

            // Add to priority queue - use negative because priority queue is min-heap
            // but we want to process vertices with highest monotonicity index first
            pq.push({v, -new_mono_index});
        }
    }

    // Sort vertices by monotonicity index (highest first)
    std::sort(result.sorted_vertices.begin(), result.sorted_vertices.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance; // Use distance field for monotonicity
              });

    return result;
}

/**
 * @brief Reconstructs the monotonic path from reference vertex to a target vertex
 *
 * @param map The monotonic reachability map from compute_monotonic_reachability_map
 * @param target_vertex The vertex to which the path should be reconstructed
 *
 * @return path_t The reconstructed path from reference to target
 */
path_t set_wgraph_t::reconstruct_monotonic_path(
    const monotonic_reachability_map_t& map,
    size_t target_vertex
    ) const {

    path_t result;

    // Check if target exists in the map
    auto it = map.info.find(target_vertex);
    if (it == map.info.end()) {
        // No path exists
        return result;
    }

    // Start from target and follow predecessors to source
    std::vector<size_t> rev_path;
    std::vector<double> rev_distances;

    size_t current = target_vertex;
    double current_distance = it->second.distance;

    // Build path in reverse
    while (current != INVALID_VERTEX) {
        rev_path.push_back(current);
        rev_distances.push_back(current_distance);

        if (current == map.ref_vertex) {
            break;
        }

        auto info_it = map.info.find(current);
        if (info_it == map.info.end()) {
            // This shouldn't happen if the map is valid
            REPORT_ERROR("Broken path: vertex %zu not found in map", current);
        }

        current = info_it->second.predecessor;
        current_distance = (current == map.ref_vertex) ? 0.0 : map.info.at(current).distance;
    }

    // Reverse path to get source to target
    std::reverse(rev_path.begin(), rev_path.end());
    std::reverse(rev_distances.begin(), rev_distances.end());

    // Fill result
    result.vertices = rev_path;
    result.distances = rev_distances;
    result.ref_vertex_index = 0; // Reference vertex is always at index 0
    result.total_weight = rev_distances.back();

    return result;
}

/**
 * @brief Finds the vertex with the best monotonicity index for gradient estimation
 *
 * This function identifies the vertex with the highest monotonicity index that
 * has a change in y-value in the desired direction from the reference vertex.
 * It can be used to find the direction of steepest consistent ascent or descent.
 *
 * @param map The monotonic reachability map from compute_monotonic_reachability_map
 * @param min_distance Minimum distance from reference vertex (to avoid local noise)
 * @param min_path_length Minimum number of vertices in the path (including reference)
 * @param ascending If true, find best path with positive total change; if false, negative total change
 *
 * @return std::pair<size_t, double> The best vertex and its monotonicity index,
 *         or {INVALID_VERTEX, 0.0} if no suitable vertex found
 */
std::pair<size_t, double> set_wgraph_t::find_best_gradient_vertex(
    const monotonic_reachability_map_t& map,
    double min_distance,
    size_t min_path_length,
    bool ascending
    ) const {

    size_t best_vertex = INVALID_VERTEX;
    double best_mono_index = -1.0;
    double best_total_change = 0.0;

    // Traverse all vertices in the map
    for (const auto& [vertex, info] : map.info) {
        // Skip reference vertex
        if (vertex == map.ref_vertex) {
            continue;
        }

        // Check direction constraint - now applied during final selection
        if ((ascending && info.total_change <= 0.0) ||
            (!ascending && info.total_change >= 0.0)) {
            continue;
        }

        // Check distance constraint
        if (info.distance < min_distance) {
            continue;
        }

        // Check path length constraint (would require reconstructing path)
        path_t path = reconstruct_monotonic_path(map, vertex);
        if (path.vertices.size() < min_path_length) {
            continue;
        }

        // Find vertex with best monotonicity index
        if (info.monotonicity_index > best_mono_index ||
            (info.monotonicity_index == best_mono_index &&
             std::abs(info.total_change) > std::abs(best_total_change))) {

            best_vertex = vertex;
            best_mono_index = info.monotonicity_index;
            best_total_change = info.total_change;
        }
    }

    return {best_vertex, best_mono_index};
}

/**
 * @brief Finds the minimum radius required to include at least domain_min_size vertices
 *
 * @param vertex The center vertex from which to measure distances
 * @param lower_bound Initial lower bound for binary search
 * @param upper_bound Initial upper bound for binary search
 * @param domain_min_size Minimum number of vertices required in the neighborhood
 * @param precision Precision threshold for terminating binary search
 *
 * @return double The minimum radius that ensures at least domain_min_size vertices
 *                are within the neighborhood of the vertex
 *
 * @details This function performs a binary search to find the smallest radius value
 * for which the vertex neighborhood (set of all vertices within the radius) contains
 * at least domain_min_size elements. This is required to ensure sufficient data
 * points for linear model fitting.
 */
double set_wgraph_t::find_minimum_radius_for_domain_min_size(
    size_t vertex,
    double lower_bound,
    double upper_bound,
    size_t domain_min_size,
    double precision
    ) const {

    // Check if the lower bound already satisfies the condition
    std::unordered_map<size_t, double> neighborhood = find_vertices_within_radius(vertex, lower_bound);
    if (neighborhood.size() >= domain_min_size) {
        return lower_bound;
    }

    // Check if the upper bound fails to satisfy the condition
    neighborhood = find_vertices_within_radius(vertex, upper_bound);
    if (neighborhood.size() < domain_min_size) {
        Rprintf("\n---------------------\nERROR\nvertex: %zu\n"
                "Not enough vertices (found: %zu, needed: %zu) within upper_bound: %.4f\n",
                vertex + 1, neighborhood.size(), domain_min_size, upper_bound);
        REPORT_ERROR("ERROR: Insufficient vertices in maximum radius\n");
    }

    // Binary search for the minimum radius
    double thld = precision * graph_diameter;
    while ((upper_bound - lower_bound) > thld) {
        double mid = (lower_bound + upper_bound) / 2.0;
        neighborhood = find_vertices_within_radius(vertex, mid);

        if (neighborhood.size() == domain_min_size) {
            upper_bound = mid;
            break;
        } else if (neighborhood.size() > domain_min_size) {
            // If condition is satisfied, try a smaller radius
            upper_bound = mid;
        } else {
            // If condition is not satisfied, try a larger radius
            lower_bound = mid;
        }
    }

    // Return the upper bound as the minimum radius that satisfies the condition
    return upper_bound;
}

/**
 * @brief Detects local extrema (maxima or minima) within a specified radius
 *
 * This function identifies vertices that are local extrema by checking if there exists
 * a disk neighborhood around each vertex where all neighbors have lower (for maxima) or higher
 * (for minima) function values than the vertex itself.
 *
 * @param y Vector of function values at each vertex
 * @param max_radius Maximum radius to search for neighborhoods
 * @param min_neighborhood_size Minimum number of vertices required in a neighborhood
 * @param detect_maxima If true, detect local maxima; if false, detect local minima
 *
 * @return Vector of local_extremum_t structures for each detected extremum
 */
std::vector<local_extremum_t> set_wgraph_t::detect_local_extrema(
    const std::vector<double>& y,
    double max_radius,
    size_t min_neighborhood_size,
    bool detect_maxima
    ) const {

    std::vector<local_extremum_t> extrema;

    // Validate input
    if (y.size() != adjacency_list.size()) {
        REPORT_ERROR("Function value vector y size (%zu) doesn't match graph size (%zu)",
                    y.size(), adjacency_list.size());
    }

    // Process each vertex as a potential extremum
    for (size_t vertex = 0; vertex < adjacency_list.size(); ++vertex) {
        // Run bounded Dijkstra's algorithm from this vertex
        std::vector<double> distances(adjacency_list.size(), std::numeric_limits<double>::infinity());
        std::vector<bool> processed(adjacency_list.size(), false);

        // Use a priority queue for Dijkstra - pair of (distance, vertex)
        std::priority_queue<std::pair<double, size_t>,
                           std::vector<std::pair<double, size_t>>,
                           std::greater<std::pair<double, size_t>>> pq;

        // Initialize start vertex
        distances[vertex] = 0.0;
        pq.push({0.0, vertex});

        // Store vertices in order of discovery for later processing
        std::vector<std::pair<size_t, double>> discovered_vertices;
        std::vector<size_t> valid_vertices;

        // Dijkstra's algorithm
        while (!pq.empty()) {
            double dist = pq.top().first;
            size_t curr = pq.top().second;
            pq.pop();

            // Skip if already processed
            if (processed[curr]) {
                continue;
            }

            // Stop if beyond radius
            if (dist > max_radius) {
                break;
            }

            // Mark as processed
            processed[curr] = true;

            // Store vertex with its distance
            if (curr != vertex) {  // Don't include the center vertex itself
                discovered_vertices.push_back({curr, dist});
            }

            // Explore neighbors
            for (const auto& edge : adjacency_list[curr]) {
                size_t neighbor = edge.vertex;
                double weight = edge.weight;

                double new_dist = dist + weight;

                // Update distance if shorter path found
                if (new_dist < distances[neighbor]) {
                    distances[neighbor] = new_dist;
                    pq.push({new_dist, neighbor});
                }
            }
        }

        // If we didn't discover enough vertices, skip
        if (discovered_vertices.size() < min_neighborhood_size) {
            continue;
        }

        // Sort vertices by increasing distance
        std::sort(discovered_vertices.begin(), discovered_vertices.end(),
                 [](const auto& a, const auto& b) {
                     return a.second < b.second;
                 });

        // Check for extremum property at different neighborhood sizes
        size_t valid_neighbors = 0;
        double valid_radius = 0.0;
        bool is_extremum = false;

        for (const auto& [neighbor, dist] : discovered_vertices) {
            // Calculate delta_y with appropriate sign
            double delta_y = y[neighbor] - y[vertex];

            // Check extremum condition based on detect_maxima flag
            bool condition_met = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (condition_met) {
                // This neighbor supports the extremum property
                valid_neighbors++;
                valid_radius = dist;  // Update the valid radius
                valid_vertices.push_back(neighbor);
            } else {
                // This neighbor violates the extremum property
                // Check if we have enough valid neighbors before this violation
                if (valid_neighbors >= min_neighborhood_size) {
                    is_extremum = true;
                }
                break;  // Stop checking further neighbors
            }
        }

        // Handle case where all neighbors support the extremum property
        if (valid_neighbors == discovered_vertices.size() && valid_neighbors >= min_neighborhood_size) {
            is_extremum = true;
            valid_radius = discovered_vertices.back().second;  // Use the furthest valid distance
        }

        // If we found an extremum, add it to the results
        if (is_extremum) {
            local_extremum_t extremum;
            extremum.vertex = vertex;
            extremum.value = y[vertex];
            extremum.radius = valid_radius;
            extremum.neighborhood_size = valid_neighbors;
            extremum.is_maximum = detect_maxima;

            extremum.vertices.reserve(valid_neighbors + 1);
            extremum.vertices = std::move(valid_vertices);

            #if 0
            extremum.vertices.push_back(vertex);
            for (const auto& [neighbor, dist] : discovered_vertices) {
                // Calculate delta_y with appropriate sign
                double delta_y = y[neighbor] - y[vertex];
                // Check extremum condition based on detect_maxima flag
                bool condition_met = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);
                if (condition_met) {
                    extremum.vertices.push_back(neighbor);
                }
            }
            #endif

            extrema.push_back(extremum);
        }
    }

    return extrema;
}

/**
 * @brief Detects local maxima within a specified radius
 *
 * Convenience wrapper around detect_local_extrema for finding local maxima.
 *
 * @param y Vector of function values at each vertex
 * @param max_radius Maximum radius to search for neighborhoods
 * @param min_neighborhood_size Minimum number of vertices required in a neighborhood
 *
 * @return Vector of local_extremum_t structures for each detected maximum
 */
std::vector<local_extremum_t> set_wgraph_t::detect_local_maxima(
    const std::vector<double>& y,
    double max_radius,
    size_t min_neighborhood_size
    ) const {

    return detect_local_extrema(y, max_radius, min_neighborhood_size, true);
}

/**
 * @brief Detects local minima within a specified radius
 *
 * Convenience wrapper around detect_local_extrema for finding local minima.
 *
 * @param y Vector of function values at each vertex
 * @param max_radius Maximum radius to search for neighborhoods
 * @param min_neighborhood_size Minimum number of vertices required in a neighborhood
 *
 * @return Vector of local_extremum_t structures for each detected minimum
 */
std::vector<local_extremum_t> set_wgraph_t::detect_local_minima(
    const std::vector<double>& y,
    double max_radius,
    size_t min_neighborhood_size
    ) const {

    return detect_local_extrema(y, max_radius, min_neighborhood_size, false);
}
