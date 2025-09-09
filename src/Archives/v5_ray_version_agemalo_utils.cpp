#include <vector>             // For std::vector
#include <queue>              // For std::priority_queue
#include <utility>            // For std::pair
#include <functional>         // For std::greater
#include <limits>             // For std::numeric_limits
#include <cmath>              // For std::fabs
#include <numeric>            // For std::accumulate

#include "set_wgraph.hpp"     // For set_wgraph_t class definition
#include "error_utils.h"      // For REPORT_ERROR()
#include "weighted_correlation.hpp" // For calculate_weighted_correlation()
#include "kernels.h"
#include "progress_utils.hpp" // for elapsed.time
#include "opt_bw.hpp"

/**
 * @brief Finds all shortest paths within a specified radius from a start vertex
 *
 * This function implements a modified version of Dijkstra's algorithm to find all shortest
 * paths from a start vertex to other vertices in the graph, limited by a maximum distance (radius).
 * The algorithm operates in two phases:
 * 1. A bounded Dijkstra's algorithm to find distances and predecessors
 * 2. Path reconstruction to generate actual paths from the collected information
 *
 * Each vertex in the reachable set is mapped to its position within the shortest path that contains it,
 * allowing efficient lookup and subpath extraction for any reachable vertex.
 *
 * @note The algorithm terminates exploration once vertex distances exceed the radius,
 *       but ensures all possible shorter paths within radius are discovered first
 *
 * The function caches results for subsequent calls with the same start vertex. (disabled!!!)
 *
 * @param start The index of the starting vertex
 * @param radius The maximum allowed distance for paths from the start vertex
 *
 * @return shortest_paths_t A structure containing:
 *         - paths: Vector of path_t objects, each containing:
 *           * vertices: Sequence of vertices in the path
 *           * total_weight: Total distance of the path
 *         - reachable_vertices: Set of all vertices reachable within the radius
 *         - vertex_to_path_map: Maps each vertex to its subpath information:
 *           * path_idx: Index of the path containing the vertex
 *           * vertex_idx: Position of the vertex within that path
 *
 * @note Paths are returned in descending order of their total weights
 * @note The vertex_to_path_map enables O(1) lookup of any vertex's position within a path,
 *       allowing efficient extraction of subpaths between any reachable vertices
 * @note Each vertex is mapped to the shortest path that contains it when multiple paths are possible
 * @note The function caches results in paths_cache for efficiency in subsequent calls
 *
 * @pre The start vertex index must be valid (0 <= start < adjacency_list.size())
 * @pre The radius must be non-negative
 *
 * @complexity Time: O((V + E) * log V) where V is the number of vertices and E is the number of edges
 *            Space: O(V + E) for storing distances, predecessors, and the result paths
 *
 * @see subpath_t Structure containing path index and vertex position information
 */
shortest_paths_t set_wgraph_t::find_graph_paths_within_radius(size_t start, double radius) const {
    // Phase 1: a bounded version of Dijkstra's algorithm
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

            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});
            }
        }
    }

	if (dp_map.empty() || dp_map.size() == 1) {  // Only contains start vertex
		// Handle case where no vertices are reachable within radius
		shortest_paths_t result;
		result.reachable_vertices.insert(start);
		REPORT_WARNING("No vertices are reachable from vertex %zu within radius %f\n",
					   start, radius);
		return result;
	}

    // Phase 2: Uses dp_map to construct a vector of corresponding paths
    shortest_paths_t result;

    struct vertex_info_t {
        size_t vertex;
        double distance;
    };

    // Create and sort vertex info
    std::vector<vertex_info_t> vertex_info;
    vertex_info.reserve(dp_map.size());
    for (const auto& [vertex, info] : dp_map) {
        if (vertex != start) {  // exclude start vertex itself
            vertex_info.emplace_back(vertex, info.first);
        }
    }

    // Sort by distance in descending order
    std::sort(vertex_info.begin(), vertex_info.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    // result.reachable_vertices will be used to keep track of used vertices to avoid duplicate paths
    result.reachable_vertices.insert(start);  // start vertex is always used

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
        std::vector<size_t> temp_vertices;
        std::vector<double> temp_distances;

        double curr_distance = info.distance;
        temp_distances.push_back(curr_distance);

        temp_vertices.push_back(curr);
        result.reachable_vertices.insert(curr);  // mark vertex as used

        while (true) {
            auto it = dp_map.find(curr);
            if (it == dp_map.end()) {
                REPORT_ERROR("Path reconstruction failed: vertex not found in path info");
            }

            curr = it->second.second;  // move to predecessor
            if (curr == -1) break;  // reached the start vertex

            double prev_distance = it->second.first;  // distance of predecessor from start
            curr_distance = prev_distance;
            temp_distances.push_back(curr_distance);

            temp_vertices.push_back(curr);
            result.reachable_vertices.insert(curr);  // mark vertex as used
        }

        // Reverse the path to get correct order (from start to end)
        std::reverse(temp_vertices.begin(), temp_vertices.end());
        std::reverse(temp_distances.begin(), temp_distances.end());

        new_path.vertices = std::move(temp_vertices);
        new_path.distances = std::move(temp_distances);

        result.paths.push_back(new_path);
    }

    // Sort paths in the descending order of their total weights
    std::sort(result.paths.begin(), result.paths.end());

    // Now populate the vertex_to_path_map with detailed subpath information
    // We'll process paths from shortest to longest (reverse of how they're sorted)
    // to ensure vertices are mapped to the shortest paths that contain them
    for (int path_idx = result.paths.size() - 1; path_idx >= 0; path_idx--) {
        const auto& path = result.paths[path_idx];

        for (size_t vertex_idx = 1; vertex_idx < path.vertices.size(); vertex_idx++) { // starting from 1 as we are not interested in paths with one vertex
            size_t vertex = path.vertices[vertex_idx];

            // Only add mapping if this vertex doesn't already have one
            // This ensures each vertex is mapped to the shortest path containing it
            if (result.vertex_to_path_map.find(vertex) == result.vertex_to_path_map.end()) {
                subpath_t subpath;
                subpath.path_idx = path_idx;
                subpath.vertex_idx = vertex_idx;

                result.vertex_to_path_map[vertex] = subpath;
            }
        }
    }

    return result;
}

/**
 * @brief Extracts x/y/w values along a geodesic ray path for use in weighted linear models
 *
 * @details For a given path in the graph, this function extracts the distance values,
 * corresponding response values, and calculates kernel weights based on distance from
 * the reference vertex. The weights are normalized so they sum to 1. This function is
 * designed to prepare data for fitting weighted linear models along geodesic rays.
 *
 * The function performs these key steps:
 * 1. Copies the vertices from the path
 * 2. Uses the pre-computed distances from the reference vertex
 * 3. Normalizes distances using the specified factor
 * 4. Applies the globally defined kernel function to calculate weights
 * 5. Normalizes weights to sum to 1
 * 6. Collects response (y) values for each vertex
 *
 * @param[in] y Vector of response values for all vertices in the original graph
 * @param[in] path The path object containing vertex indices, distances, and total weight
 * @param[in] dist_normalization_factor Factor to adjust the normalization of distances
 *
 * @return A gray_xyw_t object containing vertices, distances, response values and weights
 *
 * @note This function assumes that initialize_kernel() has been called to set up the
 * global kernel_fn pointer before this function is called
 *
 * @see initialize_kernel()
 * @see gray_xyw_t
 */
gray_xyw_t set_wgraph_t::get_xyw_along_path(
    const std::vector<double>& y,
    path_t& path,
    double dist_normalization_factor
) const {
    gray_xyw_t result;

    // Get number of vertices in the path
    size_t n_vertices = path.vertices.size();
    if (n_vertices == 0) {
        return result; // Return empty result for empty path
    }


	// Check that distances and vertices match in length
	if (path.distances.size() != n_vertices) {
		// Handle inconsistent path data
		REPORT_ERROR("Path has inconsistent data: %zu vertices but %zu distances",
					 n_vertices, path.distances.size());
		return result;
	}

    // Use the pre-computed distances directly
    result.x_path = path.distances;

	#if 0
	// weithts centered at the init vertex of each path
	{
		// Calculate distances to use for kernel weighting
		std::vector<double> normalized_distances = result.x_path;

		// Normalize distances for kernel function
		double max_dist = normalized_distances.back(); // Maximum distance in the path
		if (max_dist <= 0.0) max_dist = 1.0;          // Safety check

		max_dist *= dist_normalization_factor;        // Apply normalization factor

		for (size_t i = 0; i < n_vertices; ++i) {
			normalized_distances[i] /= max_dist;
		}

		// Calculate kernel weights
		result.w_path.resize(n_vertices);
		kernel_fn(normalized_distances.data(), n_vertices, result.w_path.data());
	}


	// uniform weights
	{
		double dn = (double)result.vertices.size();
		result.w_path.resize(n_vertices);
		for (size_t i = 0; i < result.vertices.size(); i++) {
			result.w_path[i] = 1.0 / dn;
		}
	}
    #endif

	// weithts centered at the center of the path
	{
		size_t n_path_vertices = result.x_path.size();

		// finding the mid vertex
		size_t ref_index;
		if ((n_path_vertices - 1) % 2 == 0) {
			ref_index = (n_path_vertices - 1) / 2;
		} else {
			ref_index = n_path_vertices / 2;
		}

		// Computing distances from ref_index to all other vertices
		double ref_dist = result.x_path[ref_index];
		std::vector<double> ref_vertex_distances(n_path_vertices);
		for (size_t i = 0; i < n_path_vertices; i++) {
			ref_vertex_distances[i] = std::abs(result.x_path[i] - ref_dist);
		}

		// Normalize distances for kernel function
		double max_dist = std::max(ref_vertex_distances.front(), ref_vertex_distances.back());
		if (max_dist <= 0.0) max_dist = 1.0;          // Safety check

		max_dist *= dist_normalization_factor;        // Apply normalization factor

		for (size_t i = 0; i < n_path_vertices; ++i) {
			ref_vertex_distances[i] /= max_dist;
		}

		// Calculate kernel weights
		result.w_path.resize(n_vertices);
		kernel_fn(ref_vertex_distances.data(), n_path_vertices, result.w_path.data());
	}

    // Normalize weights to sum to 1
    double total_weight = std::accumulate(result.w_path.begin(), result.w_path.end(), 0.0);
    if (total_weight > 0.0) {
        for (size_t i = 0; i < n_vertices; ++i) {
            result.w_path[i] /= total_weight;
        }
    }

    // Extract response values for each vertex in the path
    result.y_path.resize(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        result.y_path[i] = y[path.vertices[i]];
    }

	// move vertices
    result.vertices = std::move(path.vertices);

    return result;
}


/**
 * @brief Find the minimum bandwidth that ensures at least one geodesic ray has the minimum required size
 *
 * This function performs a binary search to find the smallest bandwidth value for which
 * at least one geodesic ray starting from the specified grid vertex contains at least
 * min_path_size vertices. This is crucial for ensuring sufficient data points for model fitting.
 *
 * The function takes into account the graph's packing radius, ensuring the minimum bandwidth
 * is at least the maximum packing radius to guarantee coverage of all vertices.
 *
 * @param grid_vertex The starting vertex for geodesic rays
 * @param min_path_size Minimum number of vertices required in a geodesic ray
 * @param max_bw_factor Maximum bandwidth as a fraction of graph diameter
 * @param lower_bound_factor Lower bound factor for binary search (default: 0.001)
 * @param precision Precision threshold for terminating binary search (default: 0.001)
 *
 * @return The minimum bandwidth value that satisfies the path size requirement
 *
 * @note The function uses find_graph_paths_within_radius() to discover paths at each bandwidth
 */
double set_wgraph_t::find_minimum_bandwidth(
    size_t grid_vertex,
    size_t min_path_size,
    double max_bw_factor,
    double lower_bound_factor,
    double precision) const {

    // Initialize search bounds
    double graph_diameter = this->graph_diameter;
    if (graph_diameter <= 0.0) {
        // Find diameter if not already calculated
        auto [end1, diam] = this->get_vertex_eccentricity(0);
        auto [end2, diameter] = this->get_vertex_eccentricity(end1);
        graph_diameter = diameter;
    }

    // Initial search bounds
    double lower_bound = lower_bound_factor * graph_diameter;
    double upper_bound = max_bw_factor * graph_diameter;

    // Ensure the lower bound is at least the maximum packing radius
    // This guarantees all vertices are covered by at least one ball
    if (this->max_packing_radius > 0) { // I don't think we need this if
        lower_bound = std::max(lower_bound, this->max_packing_radius);
    }

    // Check if the lower bound already satisfies the condition
    shortest_paths_t paths = this->find_graph_paths_within_radius(grid_vertex, lower_bound);
    if (has_sufficient_path_size(paths, min_path_size)) {
        return lower_bound;
    }

    // Check if the upper bound fails to satisfy the condition
    paths = this->find_graph_paths_within_radius(grid_vertex, upper_bound);
    if (!has_sufficient_path_size(paths, min_path_size)) {
		Rprintf("max_bw_factor: %.4f\ngraph_diameter: %.4f\nupper_bound: %.4f\n",
				max_bw_factor, graph_diameter, upper_bound);
		REPORT_ERROR("Even the maximum bandwidth we could not find a single geodesic ray containing at least min_path_size: %zu vertices", min_path_size);
        // return upper_bound;
    }

    // Binary search for the minimum bandwidth
    while ((upper_bound - lower_bound) > precision * graph_diameter) {
        double mid = (lower_bound + upper_bound) / 2.0;
        paths = this->find_graph_paths_within_radius(grid_vertex, mid);

        if (has_sufficient_path_size(paths, min_path_size)) {
            // If condition is satisfied, try a smaller bandwidth
            upper_bound = mid;
        } else {
            // If condition is not satisfied, try a larger bandwidth
            lower_bound = mid;
        }
    }

    // Return the upper bound as the minimum bandwidth that satisfies the condition
    return upper_bound;
}

/**
 * @brief Check if any path in the shortest_paths structure has the minimum required size
 *
 * @param paths The shortest_paths_t structure containing paths
 * @param min_path_size The minimum required path size
 * @return true if at least one path has the minimum size, false otherwise
 */
bool set_wgraph_t::has_sufficient_path_size(
    const shortest_paths_t& paths,
    size_t min_path_size) const {

    for (const auto& path : paths.paths) {
        if (path.vertices.size() >= min_path_size) {
            return true;
        }
    }

    return false;
}


/**
 * @brief Calculate directional threshold for clustering geodesic rays
 *
 * This function determines an adaptive directional threshold that controls how
 * geodesic rays are clustered based on their similarity in direction. The threshold
 * varies with bandwidth to provide appropriate clustering at different scales:
 *   - At small bandwidths (bw ≤ bw_min), we use alpha = 1.0, preserving more distinct rays
 *   - At large bandwidths (bw ≥ dir_thld_factor * graph_diameter), we use alpha = 0.25,
 *     leading to more aggressive clustering
 *   - For intermediate bandwidths, we linearly interpolate alpha
 *
 * The directional threshold is then calculated as alpha * bw.
 *
 * @param bw Current bandwidth being evaluated
 * @param bw_min Minimum valid bandwidth for the current grid vertex
 * @param graph_diameter Diameter of the graph
 * @param min_path_size Minimum path size requirement (used for minimum threshold)
 * @param dir_thld_factor Factor determining when maximum clustering begins (default: 1/3)
 * @param alpha_min Minimum alpha value for large bandwidths.
 *
 * @return The directional threshold value to be used for ray clustering
 */
double calculate_directional_threshold(
    double bw,
    double bw_min,
    double graph_diameter,
    size_t min_path_size,
    double dir_thld_factor = 1.0/3.0,
	double alpha_min = 0.25) {

    // Define bandwidth transition point where we reach maximum clustering
    double bw_transition = dir_thld_factor * graph_diameter;

    // Calculate the alpha value using piecewise linear function
    double alpha;
    if (bw <= bw_min) {
        // For small bandwidths, use alpha = 1.0 to preserve distinct rays
        alpha = 1.0;
    } else if (bw >= bw_transition) {
        // For large bandwidths, use alpha = 0.25 for aggressive clustering
        alpha = alpha_min;
    } else {
        // Linear interpolation for intermediate bandwidths
        // As bw increases from bw_min to bw_transition, alpha decreases from 1.0 to alpha_min
        double fraction = (bw - bw_min) / (bw_transition - bw_min);
        alpha = 1.0 - fraction * (1.0 - alpha_min);
    }

    // Calculate directional threshold as alpha * bw
    double directional_threshold = alpha * bw;

    // Ensure the threshold allows for minimum path coverage
    // This is a safety mechanism to ensure we can find paths with at least min_path_size vertices
    double min_threshold = min_path_size * (graph_diameter / 100.0);

    // Return the maximum of the calculated threshold and the minimum required threshold
    return std::max(directional_threshold, min_threshold);
}

/**
 * @brief Calculate Jaccard similarity threshold for path clustering
 *
 * This function determines how strict to be when clustering paths based on
 * Jaccard similarity of their vertices. Similar to directional threshold,
 * it uses a piecewise linear function to transition from strict clustering
 * at small bandwidths to more lenient clustering at large bandwidths.
 *
 * @param bw Current bandwidth being evaluated
 * @param bw_min Minimum valid bandwidth for the current grid vertex
 * @param graph_diameter Diameter of the graph
 * @param min_path_size Minimum path size requirement (used for minimum threshold)
 * @param jaccard_thld_factor Factor determining when maximum leniency begins (default: 1/3)
 * @param beta_min Minimum Jaccard threshold for large bandwidths (default: 0.3)
 *
 * @return The Jaccard similarity threshold value to use for path clustering
 */
double calculate_jaccard_threshold(
    double bw,
    double bw_min,
    double graph_diameter,
    double jaccard_thld_factor = 1.0/3.0,
    double beta_min = 0.3) {

    // Define bandwidth transition point where we reach maximum leniency
    double bw_transition = jaccard_thld_factor * graph_diameter;

    // Calculate the threshold value using piecewise linear function
    double beta;
    if (bw <= bw_min) {
        // For small bandwidths, require near-perfect similarity (beta = 1.0)
        beta = 1.0;
    } else if (bw >= bw_transition) {
        // For large bandwidths, allow more lenient clustering (beta = beta_min)
        beta = beta_min;
    } else {
        // Linear interpolation for intermediate bandwidths
        // As bw increases from bw_min to bw_transition, beta decreases from 1.0 to beta_min
        double fraction = (bw - bw_min) / (bw_transition - bw_min);
        beta = 1.0 - fraction * (1.0 - beta_min);
    }

    return beta;
}

/**
 * @brief Cluster similar geodesic rays and select representatives
 *
 * This function reduces the number of geodesic rays by clustering those that explore
 * similar directions from a reference vertex. For each cluster, it selects the ray
 * with the highest weighted correlation R-squared value as the representative.
 *
 * The clustering is controlled by a directional threshold which determines how
 * similar rays must be to belong to the same cluster. Rays are considered similar
 * if they share vertices within the directional threshold distance from the reference.
 *
 * @param grid_vertex Reference vertex from which geodesic rays originate
 * @param shortest_paths Collection of geodesic rays found within the bandwidth
 * @param y Response values at each vertex
 * @param directional_threshold Distance threshold for clustering similar rays
 * @param dist_normalization_factor Factor for normalizing distances when computing weights
 * @param bw Current bandwidth value (used for computing kernel weights)
 *
 * @return Vector of selected geodesic rays (one per direction cluster)
 */
std::vector<path_t> set_wgraph_t::cluster_and_select_rays(
    size_t grid_vertex,
    const shortest_paths_t& shortest_paths,
    const std::vector<double>& y,
    double directional_threshold,
	double jaccard_thld,
	size_t min_path_size,
    double dist_normalization_factor) const {

	if (shortest_paths.paths.empty()) {
		REPORT_WARNING("shortest_paths.paths is empty for grid_vertex: %zu\n", grid_vertex);
		return {}; // Return empty result if no paths
	}

    // Create clusters based on shared vertices within directional_threshold
    std::vector<std::vector<size_t>> clusters;
    std::vector<bool> processed(shortest_paths.paths.size(), false);

    // Process each path
    for (size_t i = 0; i < shortest_paths.paths.size(); ++i) {
        if (processed[i]) continue;

        // Start a new cluster with this path
        std::vector<size_t> cluster;
        cluster.push_back(i);
        processed[i] = true;

        const auto& path_i = shortest_paths.paths[i];

        // Find similar paths to add to this cluster
        for (size_t j = i + 1; j < shortest_paths.paths.size(); ++j) {
            if (processed[j]) continue;

            const auto& path_j = shortest_paths.paths[j];

            // Check if paths explore similar directions
            if (paths_share_direction(path_i, path_j, grid_vertex, directional_threshold, jaccard_thld)) {
                cluster.push_back(j);
                processed[j] = true;
            }
        }

        clusters.push_back(std::move(cluster));
    }

    // Select best path from each cluster based on weighted correlation
    std::vector<path_t> selected_paths;
    for (const auto& cluster : clusters) {
        size_t best_path_idx = select_best_path(cluster, shortest_paths, y, min_path_size, dist_normalization_factor);
        selected_paths.push_back(shortest_paths.paths[best_path_idx]);
    }

    return selected_paths;
}

/**
 * @brief Determine if two paths explore similar directions from a reference vertex
 *
 * Two paths are considered to explore similar directions if they share vertices within
 * the directional threshold distance from the reference vertex. This function checks
 * the early segments of both paths to see if they overlap sufficiently.
 *
 * @param path_i First path to compare
 * @param path_j Second path to compare
 * @param ref_vertex Reference vertex from which both paths originate
 * @param directional_threshold Distance threshold for similarity check
 *
 * @return true if paths explore similar directions, false otherwise
 */
bool set_wgraph_t::paths_share_direction(
    const path_t& path_i,
    const path_t& path_j,
    size_t ref_vertex,
    double directional_threshold,
	double jaccard_thld) const {

    // Find vertices within directional_threshold from the reference vertex
    std::unordered_set<size_t> path_i_vertices;
    std::unordered_set<size_t> path_j_vertices;

    // For first path, collect vertices within threshold
    for (size_t idx = 0; idx < path_i.vertices.size(); ++idx) {
        if (path_i.distances[idx] > directional_threshold) break;

        // Skip the reference vertex itself
        if (path_i.vertices[idx] != ref_vertex) {
            path_i_vertices.insert(path_i.vertices[idx]);
        }
    }

    // For second path, collect vertices within threshold
    for (size_t idx = 0; idx < path_j.vertices.size(); ++idx) {
        if (path_j.distances[idx] > directional_threshold) break;

        // Skip the reference vertex itself
        if (path_j.vertices[idx] != ref_vertex) {
            path_j_vertices.insert(path_j.vertices[idx]);
        }
    }

    // If either set is empty after excluding the reference vertex, they don't share direction
    if (path_i_vertices.empty() || path_j_vertices.empty()) {
        return false;
    }

    // Count shared vertices
    size_t intersection_size = 0;
    for (size_t vertex : path_i_vertices) {
        if (path_j_vertices.count(vertex) > 0) {
            intersection_size++;
        }
    }

    // Compute Jaccard similarity (intersection over union)
    size_t union_size = path_i_vertices.size() + path_j_vertices.size() - intersection_size;
    double jaccard = static_cast<double>(intersection_size) / union_size;

    // Paths share direction if Jaccard similarity exceeds a threshold (e.g., 0.3)
    return jaccard >= jaccard_thld;
}

/**
 * @brief Select the best path from a cluster based on weighted correlation R-squared
 *
 * For each path in the cluster, this function:
 * 1. Extracts distance (d), response (y), and calculated weights (w) along the path
 * 2. Computes weighted correlation between d and y
 * 3. Squares the correlation to get R-squared
 * 4. Selects the path with the highest R-squared value
 *
 * @param cluster Indices of paths in a directional cluster
 * @param shortest_paths Collection of all geodesic rays
 * @param y Response values at each vertex
 * @param bw Current bandwidth value (used for computing weights)
 * @param dist_normalization_factor Factor for normalizing distances
 *
 * @return Index of the best path in the cluster
 */
size_t set_wgraph_t::select_best_path(
    const std::vector<size_t>& cluster,
    const shortest_paths_t& shortest_paths,
    const std::vector<double>& y,
	size_t min_path_size,
    double dist_normalization_factor) const {

	double best_rsquared = -1.0;
	size_t best_path_idx = INVALID_VERTEX; // Initialize to invalid value
	const double MIN_ACCEPTABLE_RSQUARED = 0.1; // Or other threshold

    for (size_t path_idx : cluster) {
        path_t path = shortest_paths.paths[path_idx];

        // Skip paths that are too short
        if (path.vertices.size() < min_path_size) {
            continue;
        }

        // Get x/y/w data along the path
        gray_xyw_t xyw_data = get_xyw_along_path(y, path, dist_normalization_factor);

        // Compute weighted correlation between distance and response
        double correlation = calculate_weighted_correlation(
            xyw_data.x_path,
            xyw_data.y_path,
            xyw_data.w_path
        );

        // Compute R-squared
        double rsquared = correlation * correlation;

		// Update best path if this one has higher R-squared and meets minimum threshold
		if (rsquared > best_rsquared && rsquared >= MIN_ACCEPTABLE_RSQUARED) {
			best_rsquared = rsquared;
			best_path_idx = path_idx;
		}
    }

	// After loop, validate we found a good path
	if (best_path_idx == INVALID_VERTEX) {
		// No good paths found, choose longest path or path with most data points
		size_t longest_path_idx = cluster[0];
		size_t max_length = shortest_paths.paths[cluster[0]].vertices.size();

		for (size_t path_idx : cluster) {
			if (shortest_paths.paths[path_idx].vertices.size() > max_length) {
				max_length = shortest_paths.paths[path_idx].vertices.size();
				longest_path_idx = path_idx;
			}
		}

		return longest_path_idx;
	}

    return best_path_idx;
}


/**
 * @brief Find the optimal bandwidth that minimizes mean prediction error
 *
 * This function performs a binary search to find the bandwidth value that
 * minimizes the mean prediction error across all ray models. For each
 * evaluated bandwidth, it:
 * 1. Finds all geodesic rays within that radius
 * 2. Clusters similar rays and selects representatives
 * 3. Fits weighted linear models to each selected ray
 * 4. Computes the average prediction error across all models
 *
 * @param grid_vertex The grid vertex from which to find geodesic rays
 * @param y Response values at each vertex
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param min_bw Minimum bandwidth to consider (must allow at least one valid ray)
 * @param max_bw Maximum bandwidth to consider (typically a fraction of graph diameter)
 * @param dist_normalization_factor Factor for normalizing distances in kernel calculations
 * @param y_binary Whether the response variable is binary (affects model fitting)
 * @param tolerance Convergence tolerance for model fitting
 * @param precision Precision for bandwidth binary search
 * @param weights Optional weights for vertices (for bootstrapping)
 *
 * @return A pair containing the optimal bandwidth and the fitted models at that bandwidth
 */
std::pair<double, std::vector<ext_ulm_t>> set_wgraph_t::find_optimal_bandwidth(
    size_t grid_vertex,
    const std::vector<double>& y,
    size_t min_path_size,
    double min_bw,
    double max_bw,
    double dist_normalization_factor,
    bool y_binary,
    double tolerance,
    double precision,
	bool verbose,
    const std::optional<std::vector<double>>& weights) const {

	auto start_time = std::chrono::steady_clock::now();
    bool is_time_consuming = false;

    // Check if min_bw already provides good models
    auto [lower_error, lower_models] = evaluate_bandwidth(
        grid_vertex, y, min_path_size, min_bw, min_bw,
        dist_normalization_factor, y_binary, tolerance, weights
    );

    // Check if max_bw provides better models
    auto [upper_error, upper_models] = evaluate_bandwidth(
        grid_vertex, y, min_path_size, max_bw, min_bw,
        dist_normalization_factor, y_binary, tolerance, weights
    );

    // If min_bw gives better results, use it
    if (lower_error <= upper_error) {
        return {min_bw, lower_models};
    }

    // Initialize binary search
    double lower_bound = min_bw;
    double upper_bound = max_bw;
    double best_bw = min_bw;
    double best_error = lower_error;
    std::vector<ext_ulm_t> best_models = lower_models;

	// Binary search for optimal bandwidth
	int iteration_count = 0;
	const int MAX_ITERATIONS = 50; // Reasonable upper limit

	while ((upper_bound - lower_bound) > precision * max_bw && iteration_count < MAX_ITERATIONS) {
		iteration_count++;

		// If taking many iterations, mark as time-consuming
        if (iteration_count > 10) {
            is_time_consuming = true;
        }

        // Try two points in the interval to determine direction
        double mid1 = lower_bound + (upper_bound - lower_bound) / 3;
        double mid2 = upper_bound - (upper_bound - lower_bound) / 3;

        auto [error1, models1] = evaluate_bandwidth(
            grid_vertex, y, min_path_size, mid1, min_bw,
            dist_normalization_factor, y_binary, tolerance, weights
        );

        auto [error2, models2] = evaluate_bandwidth(
            grid_vertex, y, min_path_size, mid2, min_bw,
            dist_normalization_factor, y_binary, tolerance, weights
        );

        // Update best bandwidth if found
        if (error1 < best_error) {
            best_error = error1;
            best_bw = mid1;
            best_models = models1;
        }

        if (error2 < best_error) {
            best_error = error2;
            best_bw = mid2;
            best_models = models2;
        }

        // Narrow search range based on which third has the better value
        if (error1 < error2) {
            upper_bound = mid2;
        } else {
            lower_bound = mid1;
        }
    }

	if (iteration_count >= MAX_ITERATIONS) {
		REPORT_WARNING("Binary search for optimal bandwidth did not converge in %d iterations.\n", MAX_ITERATIONS);
		is_time_consuming = true;
	}

	if (verbose && is_time_consuming) {
        elapsed_time(start_time, "Bandwidth optimization completed", true);
    }

    return {best_bw, best_models};
}

/**
 * @brief Find optimal bandwidth using a grid search over candidate bandwidths
 *
 * This function evaluates prediction error at a set of candidate bandwidths
 * evenly distributed over a range, either on a linear or logarithmic scale.
 *
 * @param grid_vertex The grid vertex from which to find geodesic rays
 * @param y Response values at each vertex
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param min_bw Minimum bandwidth to consider
 * @param max_bw Maximum bandwidth to consider
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing; if false, use linear spacing
 * @param dist_normalization_factor Factor for normalizing distances
 * @param y_binary Whether response is binary (0/1)
 * @param tolerance Convergence tolerance for model fitting
 * @param precision Precision for numerical comparisons
 * @param verbose Whether to print progress information
 * @param weights Optional vertex weights (for bootstrapping)
 *
 * @return opt_bw_t structure containing candidate bandwidths, errors, optimal bandwidth, and models
 */
opt_bw_t set_wgraph_t::find_optimal_bandwidth_over_grid(
    size_t grid_vertex,
    const std::vector<double>& y,
    size_t min_path_size,
    double min_bw,
    double max_bw,
    size_t n_bws,
    bool log_grid,
    double dist_normalization_factor,
    bool y_binary,
    double tolerance,
    double precision,
    bool verbose,
    const std::optional<std::vector<double>>& weights) const {

    auto start_time = std::chrono::steady_clock::now();
    bool is_time_consuming = false;

    // Initialize result structure
    opt_bw_t result;

    // Handle the case where min_bw and max_bw are very close
    double bw_range = max_bw - min_bw;
    double min_spacing = precision;  // Minimum meaningful difference between bandwidths

    // If range is too small for requested number of bandwidths
    if (bw_range < (n_bws - 1) * min_spacing) {
        // Reduce the number of bandwidths based on the available range
        size_t effective_n_bws = std::max(size_t(2), size_t(bw_range / min_spacing) + 1);
        if (effective_n_bws < n_bws) {
            REPORT_WARNING("Warning: Bandwidth range (%.6f) is too small for %zu distinct values. "
                           "Reducing to %zu bandwidths.\n",
                           bw_range, n_bws, effective_n_bws);
            n_bws = effective_n_bws;
        }
    }

    // Generate candidate bandwidths using specified spacing strategy
    result.bws.reserve(n_bws);

    // Always include min_bw
    result.bws.push_back(min_bw);

    if (n_bws > 2) {
        // Generate intermediate bandwidths
        if (log_grid && min_bw > 0) {
            // Logarithmic spacing
            double log_min = std::log(min_bw);
            double log_max = std::log(max_bw);
            double log_step = (log_max - log_min) / (n_bws - 1);

            for (size_t i = 1; i < n_bws - 1; ++i) {
                double log_bw = log_min + i * log_step;
                result.bws.push_back(std::exp(log_bw));
            }
        } else {
            // Linear spacing
            double step = bw_range / (n_bws - 1);

            for (size_t i = 1; i < n_bws - 1; ++i) {
                result.bws.push_back(min_bw + i * step);
            }
        }
    }

    // Always include max_bw
    if (n_bws > 1) {
        result.bws.push_back(max_bw);
    }

    // Ensure we have the right number of bandwidths
    if (result.bws.size() != n_bws) {
        REPORT_WARNING("Warning: Generated %zu bandwidths instead of requested %zu\n",
                       result.bws.size(), n_bws);
    }

    // Initialize storage for errors
    result.errors.resize(result.bws.size(), std::numeric_limits<double>::infinity());

    // Track best bandwidth and models
    double best_error = std::numeric_limits<double>::infinity();
    size_t best_idx = 0;

    // Evaluate each candidate bandwidth
    for (size_t i = 0; i < result.bws.size(); ++i) {
        double bw = result.bws[i];

        // Evaluate this bandwidth
        auto [error, models] = evaluate_bandwidth(
            grid_vertex, y, min_path_size, bw, min_bw,
            dist_normalization_factor, y_binary, tolerance, weights
        );

        // Store error
        result.errors[i] = error;

        // Update best if this is better
        if (error < best_error) {
            best_error = error;
            best_idx = i;
            result.models = models;
        }

        // Check for time-consuming operation
        if (i > 0 && i % 5 == 0) {
            is_time_consuming = true;
        }
    }

    // Set optimal bandwidth
    result.opt_bw = result.bws[best_idx];

    if (verbose && is_time_consuming) {
        elapsed_time(start_time, "Grid bandwidth search completed", true);

        // Print bandwidth search results
        Rprintf("Bandwidth search for vertex %zu:\n", grid_vertex);
        Rprintf("  Best bandwidth: %.6f (index %zu of %zu)\n",
                result.opt_bw, best_idx, result.bws.size());
        Rprintf("  Error range: %.6f to %.6f\n",
                *std::min_element(result.errors.begin(), result.errors.end()),
                *std::max_element(result.errors.begin(), result.errors.end()));
    }

    return result;
}


/**
 * @brief Evaluate the prediction error for models at a specific bandwidth
 *
 * This function evaluates the mean prediction error for a given bandwidth by:
 * 1. Finding all geodesic rays within the bandwidth radius
 * 2. Clustering similar rays and selecting representatives
 * 3. Fitting weighted linear models to each selected ray
 * 4. Computing the average prediction error across all models
 *
 * @param grid_vertex The grid vertex from which to find geodesic rays
 * @param y Response values at each vertex
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param bw Bandwidth to evaluate
 * @param dist_normalization_factor Factor for normalizing distances in kernel calculations
 * @param y_binary Whether the response variable is binary (affects model fitting)
 * @param tolerance Convergence tolerance for model fitting
 * @param weights Optional weights for vertices (for bootstrapping)
 *
 * @return A pair containing the mean prediction error and the fitted models
 */
std::pair<double, std::vector<ext_ulm_t>> set_wgraph_t::evaluate_bandwidth(
    size_t grid_vertex,
    const std::vector<double>& y,
    size_t min_path_size,
    double bw,
	double min_bw,
    double dist_normalization_factor,
    bool y_binary,
    double tolerance,
    const std::optional<std::vector<double>>& weights) const {

    // Find all geodesic rays within bandwidth radius
    shortest_paths_t shortest_paths = find_graph_paths_within_radius(grid_vertex, bw);

	if (shortest_paths.paths.empty()) {
		// No paths found at this bandwidth
		return {std::numeric_limits<double>::infinity(), {}};
	}

	// Calculate directional threshold for this bandwidth
    double directional_threshold = calculate_directional_threshold(
        bw,
        min_bw,
        graph_diameter,
        min_path_size
		);

	// Calculate Jaccard threshold for this bandwidth
	double jaccard_thld_factor = 1.0/3.0;
	double beta_min = 0.3;
	double jaccard_thld = calculate_jaccard_threshold(
		bw,
		min_bw,
		graph_diameter,
		jaccard_thld_factor,
		beta_min
		);

    // Cluster similar rays and select representatives
    std::vector<path_t> selected_paths = cluster_and_select_rays(
        grid_vertex,
        shortest_paths,
        y,
        directional_threshold,
        jaccard_thld,
		min_path_size,
        dist_normalization_factor
    );

	if (selected_paths.empty()) {
		// No representative paths after clustering
		return {std::numeric_limits<double>::infinity(), {}};
	}

    // Fit weighted linear models to selected paths
    std::vector<ext_ulm_t> models;
    models.reserve(selected_paths.size());

    double total_error = 0.0;
    size_t valid_model_count = 0;

    for (auto& path : selected_paths) {
        // Skip paths that are too short
        if (path.vertices.size() < min_path_size) {
            continue;
        }

        // Extract data along the path for model fitting
        gray_xyw_t xyw_path = get_xyw_along_path(y, path, dist_normalization_factor);

        // Apply additional weights if provided (for bootstrap samples)
        if (weights) {
            for (size_t i = 0; i < xyw_path.vertices.size(); i++) {
                xyw_path.w_path[i] *= (*weights)[xyw_path.vertices[i]];
            }
        }

        // Fit weighted linear model to the data
		int n_iter = 0;
		double robust_scale = 6.0;
        ulm_t fit_result = cleveland_ulm(
            xyw_path.x_path.data(),
            xyw_path.y_path.data(),
            xyw_path.w_path,
            y_binary,
            tolerance,
			n_iter,
			robust_scale
        );

        // Ensure we have error estimates
        if (fit_result.errors.empty()) {
            REPORT_ERROR("fit_result.errors is empty\n");
        }

        // Create extended model object with all metadata needed for future use
        ext_ulm_t extended_result(fit_result);
        extended_result.bw = bw;
        extended_result.vertices = xyw_path.vertices;
        extended_result.x_path = xyw_path.x_path;
        extended_result.y_path = xyw_path.y_path;
        extended_result.w_path = xyw_path.w_path;
        extended_result.predictions = fit_result.predictions;
        extended_result.errors = fit_result.errors;

        // Calculate mean error for this model
        double mean_error = std::accumulate(fit_result.errors.begin(),
                                         fit_result.errors.end(), 0.0) / fit_result.errors.size();
        extended_result.mean_error = mean_error;

        // Add this model's error to the total for bandwidth evaluation
        total_error += mean_error;
        valid_model_count++;

        // Add this model to our collection
        models.push_back(std::move(extended_result));
    }

	if (valid_model_count == 0) {
		Rprintf("grid_vertex: %zu bw: %.4f\n",
				grid_vertex, bw);
		REPORT_ERROR("valid_model_count is 0\n");
	}

    // Compute average error across all models for this bandwidth
    double mean_error = valid_model_count > 0 ?
                      total_error / valid_model_count :
                      std::numeric_limits<double>::infinity();

    // double mean_error = models.empty() ? std::numeric_limits<double>::infinity()
	// 	: total_error / models.size();

    return {mean_error, models};
}
