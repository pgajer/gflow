#include "find_and_filter_composite_paths.h"

/**
 * @brief Combines search and filtering operations to find qualified paths with adaptive radius adjustment
 *
 * @details This function provides a complete solution for finding paths that meet specific size
 * requirements by:
 *   1. Computing initial reachability map using the provided radius
 *   2. Reconstructing paths and checking if they meet size criteria
 *   3. Adaptively increasing the search radius if needed until suitable paths are found
 *   4. Filtering the final paths to include only those meeting the requirements
 *
 * The function effectively combines multiple path-finding steps into a single call, making
 * client code simpler and more maintainable.
 *
 * @param ref_vertex Reference vertex from which to start the search
 * @param initial_radius Initial search radius to use
 * @param min_path_size Minimum number of original graph vertices required in paths
 * @param min_num_grid_vertices Minimum number of grid vertices required in paths
 * @param graph_diameter Diameter of the original graph (used to determine maximum allowed radius)
 * @param final_radius Optional pointer to store the final radius used for successful path discovery
 *
 * @return Reachability map containing only paths that satisfy the size requirements
 */
reachability_map_t uniform_grid_graph_t::find_and_filter_composite_paths(
	size_t ref_vertex,
	double initial_radius,
	size_t min_path_size,
	size_t min_num_grid_vertices,
	double graph_diameter,
	double* final_radius
	) const {
	// Start with the initial radius
	double current_radius = initial_radius;

	// Maximum allowed radius based on graph diameter (to prevent excessive searching)
	double max_allowed_radius = graph_diameter * 0.9;
	// if (initial_radius >= max_allowed_radius) {
	// 	max_allowed_radius = graph_diameter * 0.75;
	// }

	// Keep track of our best reachability map
	reachability_map_t best_map;
	bool found_suitable_paths = false;

	// Try increasingly larger radii until we find suitable paths or reach the maximum
	while (current_radius <= max_allowed_radius) {
		// Compute reachability map with current radius
		reachability_map_t current_map = compute_reachability_map(ref_vertex, current_radius);

		// Reconstruct paths from the reachability map
		auto grid_vertex_paths = reconstruct_paths(current_map);

		// Create composite paths that combine grid and original vertices
		auto composite_paths = create_composite_paths(grid_vertex_paths);

		// Check if we have at least one path meeting our size requirements
		if (has_min_size_path(composite_paths, min_path_size, min_num_grid_vertices)) {
			// Filter composite paths to include only those meeting requirements
			auto filtered_paths = filter_composite_paths(
				composite_paths,
				min_path_size,
				min_num_grid_vertices
				);

			// Now we need to filter the reachability map to include only vertices in filtered paths
			// First, collect all vertices from filtered paths
			std::unordered_set<size_t> vertices_to_keep;
			for (const auto& path : filtered_paths) {
				// Add original vertices
				vertices_to_keep.insert(path.vertices.begin(), path.vertices.end());
				// Add grid vertices
				vertices_to_keep.insert(path.grid_vertices.begin(), path.grid_vertices.end());
			}

			// Create a filtered reachability map
			reachability_map_t filtered_map;
			filtered_map.ref_vertex = current_map.ref_vertex;

			// Include only vertices that are in our filtered paths
			for (const auto& [vertex, distance] : current_map.distances) {
				if (vertices_to_keep.count(vertex) > 0 || vertex == ref_vertex) {
					filtered_map.distances[vertex] = distance;
					filtered_map.predecessors[vertex] = current_map.predecessors.at(vertex);

					// Add to sorted vertices if it's an original vertex (not the reference)
					if (vertex != ref_vertex && is_original_vertex(vertex)) {
						filtered_map.sorted_vertices.push_back({vertex, distance});
					}
				}
			}

			// Re-sort the sorted_vertices by distance
			std::sort(filtered_map.sorted_vertices.begin(), filtered_map.sorted_vertices.end(),
					  [](const vertex_info_t& a, const vertex_info_t& b) {
						  return a.distance > b.distance;
					  });

			// Store this as our best map
			best_map = std::move(filtered_map);
			found_suitable_paths = true;
			break;
		}

		// Increase the radius for the next iteration
		current_radius *= 1.5;
	}

	// Store the final radius if the caller requested it
	if (final_radius) {
		*final_radius = current_radius;
	}

	// If we couldn't find suitable paths, log a warning
	if (!found_suitable_paths) {
		REPORT_WARNING("Warning: Could not find paths of size %zu even with extended search radius.\n", min_path_size);
	}

	return best_map;
}
