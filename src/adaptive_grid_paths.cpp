#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h"
#include "kernels.h"

#include <vector>
#include <queue>
#include <set>
#include <unordered_set>
#include <cmath>
#include <utility>
#include <tuple>
#include <numeric>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/**
 * @brief Discovers paths that meet minimum size requirements by adaptively adjusting search radius
 *
 * @details This method attempts to find paths that contain at least min_path_size original graph
 * vertices and min_num_grid_vertices grid vertices. If such paths are not found with the initial
 * bandwidth (current_max_bw), the method incrementally increases the search radius up to a maximum
 * determined by the graph_diameter until suitable paths are found or the maximum allowed bandwidth
 * is reached.
 *
 * @param grid_vertex Starting grid vertex for the search
 * @param min_path_size Minimum number of original graph vertices required in paths
 * @param min_num_grid_vertices Minimum number of grid vertices required in paths
 * @param current_max_bw Initial maximum bandwidth to use for search
 * @param graph_diameter Diameter of the original graph (used to determine maximum allowed bandwidth)
 * @param final_max_bw Optional pointer to store the final bandwidth used for successful path discovery
 *
 * @return Vector of composite paths that meet the size requirements, or empty vector if no such paths found
 *
 * @note This method only returns paths that satisfy the size requirements, filtering out any that don't.
 */
std::vector<compose_path_t> uniform_grid_graph_t::find_min_size_composite_paths(
	size_t grid_vertex,
	size_t min_path_size,
	size_t min_num_grid_vertices,
	double current_max_bw,
	double graph_diameter,
	double* final_max_bw
	) const {

	double max_allowed_bw = graph_diameter * 0.9;
	// if (current_max_bw >= max_allowed_bw) {
	// 	max_allowed_bw = graph_diameter * 0.75;
	// }

	// Start with 1.5 times the current bandwidth
	double new_max_bw = current_max_bw * 1.5;
	std::vector<compose_path_t> result;

	while (new_max_bw <= max_allowed_bw) {
		// Compute new reachability map with increased radius
		auto reachability_map  = compute_reachability_map(grid_vertex, new_max_bw);
		auto grid_vertex_paths = reconstruct_paths(reachability_map);
		auto composite_paths   = create_composite_paths(grid_vertex_paths);

		// Filter paths that meet criteria
		result = filter_composite_paths(composite_paths, min_path_size, min_num_grid_vertices);

		// If we found at least one path of sufficient size
		if (!result.empty()) {
			if (final_max_bw) *final_max_bw = new_max_bw;
			return result;
		}

		// Increase radius for next iteration
		new_max_bw *= 1.5;
	}

	if (final_max_bw) *final_max_bw = new_max_bw;
	REPORT_WARNING("Warning: Could not find paths of size %zu even with extended search radius.\n", min_path_size);

	return result; // Empty vector
}

/**
 * @brief Filters a collection of paths to include only those meeting minimum size requirements
 *
 * @details This function examines each path in the input collection and selects only those
 * that contain at least min_path_size original graph vertices and min_num_grid_vertices grid
 * vertices. This filtering is useful for ensuring that paths have sufficient coverage of both
 * the original graph and the grid representation for accurate modeling.
 *
 * @param paths Vector of composite paths to filter
 * @param min_path_size Minimum number of original graph vertices required in a valid path
 * @param min_num_grid_vertices Minimum number of grid vertices required in a valid path
 *
 * @return Vector containing only the paths that satisfy both size requirements
 */
std::vector<compose_path_t> uniform_grid_graph_t::filter_composite_paths(
	const std::vector<compose_path_t>& paths,
	size_t min_path_size,
	size_t min_num_grid_vertices
	) const {
	// Create a result vector to store the filtered paths
	std::vector<compose_path_t> filtered_paths;

	// Reserve space in the filtered_paths vector to avoid multiple reallocations
	// This is an optimization that can significantly improve performance when
	// we expect a substantial number of paths to meet our criteria
	filtered_paths.reserve(paths.size());

	// Iterate through each path in the input collection
	for (const auto& path : paths) {
		// Check if the current path meets both size requirements:
		// 1. At least min_path_size original graph vertices
		// 2. At least min_num_grid_vertices grid vertices
		if (path.vertices.size() >= min_path_size &&
			path.grid_vertices.size() >= min_num_grid_vertices) {
			// This path meets our criteria, so add it to the filtered collection
			filtered_paths.push_back(path);
		}
	}

	// If we allocated more space than needed, we can shrink the container
	// to free up unused memory (optional optimization)
	filtered_paths.shrink_to_fit();

	// Return the collection of paths that passed our filtering criteria
	return filtered_paths;
}

/**
 * @brief Combines search and filtering operations to find qualified paths with adaptive radius adjustment
 *
 * @details This function provides a complete solution for finding paths that meet specific size
 * requirements by:
 *   1. Computing initial reachability map using the provided radius
 *   2. Reconstructing paths and checking if they meet size criteria
 *   3. Adaptively increasing the search radius if needed until suitable paths are found
 *   4. Filtering the final composite paths to include only those meeting the requirements
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
 * @return Vector of composite paths that satisfy the size requirements
 */
std::vector<compose_path_t> uniform_grid_graph_t::find_and_filter_composite_paths(
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
	double max_allowed_radius = graph_diameter;
	// if (initial_radius >= max_allowed_radius) {
	// 	max_allowed_radius = graph_diameter * 0.75;
	// }

	if (ref_vertex == 93) {
		Rprintf("In uniform_grid_graph_t::find_and_filter_composite_paths()\n");
		Rprintf("initial_radius: %.3f\n",initial_radius);
		Rprintf("graph_diameter: %.3f\n", graph_diameter);
		Rprintf("max_allowed_radius: %.3f\n", max_allowed_radius);
	}

	// Keep track of our best paths
	std::vector<compose_path_t> filtered_paths;

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
			filtered_paths = filter_composite_paths(
				composite_paths,
				min_path_size,
				min_num_grid_vertices
				);

			// We found suitable paths, so we're done
			break;
		}

		// Increase the radius for the next iteration
		current_radius *= 1.5;
	}

	// Store the final radius if the caller requested it
	if (final_radius) {
		*final_radius = current_radius;
	}

	// If we couldn't find suitable paths, log a Rf_warning
	if (filtered_paths.empty()) {
		REPORT_WARNING("In uniform_grid_graph_t::find_and_filter_composite_paths() Warning: Could not find paths of size %zu even with extended search radius.\n",
					   min_path_size);
	}

	return filtered_paths;
}
