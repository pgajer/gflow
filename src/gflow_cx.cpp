#include <queue>
#include "set_wgraph.hpp"
#include <R_ext/Print.h>

/**
 * @brief Compute the local extremum hop index for a given vertex.
 *
 * @details This function determines the maximum number of hops for which a vertex remains
 *          a strict local extremum. If the vertex is not a neighborhood local extremum,
 *          the returned index is 0. Otherwise, the function performs a breadth-first search
 *          (BFS) starting from the given vertex and identifies the minimal hop h at which a
 *          violation occurs (i.e., a vertex with a more extreme value is found).
 *
 * The hop index is defined as:
 * - For a local minimum: the smallest hop distance at which y[w] < y[v]
 * - For a local maximum: the smallest hop distance at which y[w] > y[v]
 *
 * If no such vertex is found, the hop index is the graph diameter.
 *
 * @param vertex Index of the vertex to evaluate
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, checks for maxima; otherwise, for minima
 * @return size_t The local extremum hop index (0 if not a local extremum)
 */
size_t set_wgraph_t::compute_extremum_hop_index(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
) const {
    // Check if vertex is a local extremum in its 1-hop neighborhood
    for (const auto& edge : adjacency_list[vertex]) {
        size_t u = edge.vertex;
        if (detect_maxima) {
            if (y[u] >= y[vertex]) return 0;
        } else {
            if (y[u] <= y[vertex]) return 0;
        }
    }

    // BFS to find first violating vertex
    std::vector<bool> visited(adjacency_list.size(), false);
    std::queue<std::pair<size_t, size_t>> q; // (vertex, hop count)
    visited[vertex] = true;
    q.push({vertex, 0});

    while (!q.empty()) {
        auto [v, hop] = q.front();
        q.pop();

        for (const auto& edge : adjacency_list[v]) {
            size_t u = edge.vertex;
            if (visited[u]) continue;
            visited[u] = true;

            if ((detect_maxima && y[u] > y[vertex]) || (!detect_maxima && y[u] < y[vertex])) {
                return hop + 1; // violation found at this hop
            }

            q.push({u, hop + 1});
        }
    }

    // No violation found: return max possible distance (graph diameter)
    return graph_diameter;
}

/**
 * @brief Computes the hop neighborhood where a given vertex ceases to be a local extremum.
 *
 * @details This function determines the maximum hop distance for which a vertex remains
 *          a strict local extremum and returns the corresponding hop neighborhood structure.
 *
 *          The function first checks if the input vertex is a **global extremum** (i.e., the global
 *          minimum or maximum of the function y). If so, the function terminates early and sets
 *          `hop_idx` to `std::numeric_limits<size_t>::max()` to indicate this special case.
 *
 *          Otherwise, the function performs a breadth-first search (BFS) starting from the given vertex
 *          and identifies the minimal hop h at which a violation occurs (i.e., a vertex with
 *          a more extreme value is found). It processes vertices level by level to ensure that:
 *          1. All vertices at each hop level are fully explored before moving to the next level
 *          2. When violations are found at a level, all vertices at that level are still processed
 *          3. The search stops after processing the hop level where violations first occur
 *
 * The hop index is defined as:
 * - For a local minimum: the smallest hop distance at which y[w] < y[v]
 * - For a local maximum: the smallest hop distance at which y[w] > y[v]
 *
 * Special cases:
 * - If the vertex is not a 1-hop local extremum, hop_idx is 0.
 * - If the vertex is a global extremum, hop_idx is set to std::numeric_limits<size_t>::max().
 * - If no violation is found in the entire graph, hop_idx is the graph diameter.
 *
 * @param vertex Index of the vertex to evaluate
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, checks for maxima; otherwise, for minima
 * @return hop_nbhd_t Structure containing the vertex, hop index, hop distance map, and boundary map
 */
hop_nbhd_t set_wgraph_t::compute_extremum_hop_nbhd(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
) const {

    hop_nbhd_t result;
    result.vertex  = vertex;

	// Checking if 'vertex' is a global minimum or maximum of y
	size_t min_idx = 0, max_idx = 0;
	for (size_t i = 1; i < y.size(); ++i) {
		if (y[i] < y[min_idx]) min_idx = i;
		if (y[i] > y[max_idx]) max_idx = i;
	}

	if ((detect_maxima && vertex == max_idx) || (!detect_maxima && vertex == min_idx)) {
		// Rprintf("vertex: %zu is a global %s\n", vertex, detect_maxima ? "maximum" : "minimum");
		result.hop_idx = std::numeric_limits<size_t>::max();
		return result;
	}

    // Check if vertex is a local extremum in its 1-hop neighborhood
    bool is_local_extremum = true;
    for (const auto& edge : adjacency_list[vertex]) {
        size_t u = edge.vertex;
        if (detect_maxima) {
            if (y[u] >= y[vertex]) {
                is_local_extremum = false;
                break;
            }
        } else {
            if (y[u] <= y[vertex]) {
                is_local_extremum = false;
                break;
            }
        }
    }

    if (!is_local_extremum) {
        result.hop_idx = 0;
        return result;
    }

    // BFS to find the extremum hop neighborhood
    std::vector<bool> visited(adjacency_list.size(), false);
    std::vector<size_t> current_level;  // Vertices at current level
    std::vector<size_t> next_level;     // Vertices at next level

    // Initialize with the start vertex
    result.hop_dist_map[vertex] = 0;
    visited[vertex] = true;
    current_level.push_back(vertex);

    size_t current_hop = 0;
    bool violation_found = false;

    // Process BFS level by level
    while (!current_level.empty()) {
        // Process all vertices at the current level
        for (size_t v : current_level) {
            for (const auto& edge : adjacency_list[v]) {
                size_t u = edge.vertex;
                if (visited[u]) continue;

                visited[u] = true;

                // Check if this vertex violates the extremum condition
                bool is_violation = (detect_maxima && y[u] > y[vertex]) ||
                                   (!detect_maxima && y[u] < y[vertex]);

                if (is_violation) {
                    violation_found = true;
                    // Don't add to boundary map yet, we'll do that in a separate pass
                } else {
                    // This vertex is part of the next level to explore
                    next_level.push_back(u);

                    // Add to hop_dist_map as it's at current_hop + 1
                    result.hop_dist_map[u] = current_hop + 1;
                }
            }
        }

        // If we found violations at this level, we've found our hop index
        if (violation_found) {
            result.hop_idx = current_hop;
			// Rprintf("Found violation at current_hop: %zu\n", current_hop);
            break;
        }

        // Move to the next level
        current_hop++;
        current_level.swap(next_level);
        next_level.clear();
    }

	// Rprintf("current_hop: %zu\n", current_hop);

    // If we've exited the loop with violations found, we need to:
    // 1. Remove vertices from hop_dist_map with hop distance > result.hop_idx
    // 2. Ensure we have all boundary vertices at hop_idx + 1
    if (violation_found) {
        // Remove vertices with hop distance > result.hop_idx
        std::vector<size_t> vertices_to_remove;
        for (const auto& [v, dist] : result.hop_dist_map) {
            if (dist > result.hop_idx) {
                vertices_to_remove.push_back(v);
            }
        }

        for (size_t v : vertices_to_remove) {
            result.hop_dist_map.erase(v);
        }

        // We need another BFS to find all vertices at exactly hop_idx + 1
        // Reset visited flags except for vertices we've already processed
        std::fill(visited.begin(), visited.end(), false);
        for (const auto& [v, _] : result.hop_dist_map) {
            visited[v] = true;
        }

        // Start a new BFS from vertices at hop_idx
        std::vector<size_t> boundary_search;
        for (const auto& [v, dist] : result.hop_dist_map) {
            if (dist == result.hop_idx) {
                boundary_search.push_back(v);
            }
        }

        // Find all vertices at hop_idx + 1
        for (size_t v : boundary_search) {
            for (const auto& edge : adjacency_list[v]) {
                size_t u = edge.vertex;
                if (visited[u]) continue;

                visited[u] = true;

                // Add all vertices at hop_idx + 1 to the boundary map, regardless of whether
                // they violate the extremum condition (we want the complete boundary)
                result.y_nbhd_bd_map[u] = y[u];
            }
        }
    } else {
		result.hop_idx = current_hop;
	}

    return result;
}


/**
 * @brief Analyze extrema and their hop neighborhoods in a function defined on a graph
 *
 * @param y Vector of function values at each vertex
 * @return A pair of maps from extrema vertices to their hop neighborhoods for minima and maxima
 */
std::pair<std::unordered_map<size_t, hop_nbhd_t>, std::unordered_map<size_t, hop_nbhd_t>>
set_wgraph_t::compute_extrema_hop_nbhds(
	const std::vector<double>& y
	) const {

	// Find all local extrema
	auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

	// Initialize maps for results
	std::unordered_map<size_t, hop_nbhd_t> lmin_hop_nbhd_map;
	std::unordered_map<size_t, hop_nbhd_t> lmax_hop_nbhd_map;

	// Compute hop neighborhoods for all local minima
	bool detect_maxima = false;
	for (const auto& vertex : nbr_lmin) {
		hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, y, detect_maxima);
		lmin_hop_nbhd_map[vertex] = hop_nbhd;
	}

	// Compute hop neighborhoods for all local maxima
	detect_maxima = true;
	for (const auto& vertex : nbr_lmax) {
		hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, y, detect_maxima);
		lmax_hop_nbhd_map[vertex] = hop_nbhd;
	}

	return {lmin_hop_nbhd_map, lmax_hop_nbhd_map};
}

std::vector<size_t> set_wgraph_t::compute_ascending_gradient_trajectory(
	const std::vector<double>& values,
	size_t start_vertex) const {

	std::vector<size_t> trajectory = {start_vertex};
	size_t current = start_vertex;
	bool is_local_max = false;

	while (!is_local_max) {
		double max_increase = 0.0;
		size_t next_vertex = current;  // Default to staying at current (local max)

		for (const auto& edge : adjacency_list[current]) {
			size_t neighbor = edge.vertex;
			double increase = values[neighbor] - values[current];

			if (increase > max_increase) {
				max_increase = increase;
				next_vertex = neighbor;
			}
		}

		// If we can't find a higher neighbor, we've reached a local maximum
		if (next_vertex == current) {
			is_local_max = true;
		} else {
			trajectory.push_back(next_vertex);
			current = next_vertex;
		}
	}

	return trajectory;
}

/**
 * @brief Creates a gradient flow complex with iterative harmonic extension
 *
 * @details This function identifies local extrema in a function defined on a graph and
 *          iteratively applies harmonic extension to smooth out spurious extrema. The process
 *          continues until either all spurious extrema (those with hop index <= threshold)
 *          are eliminated or the maximum number of iterations is reached.
 *
 *          The algorithm follows these key steps:
 *          1. Identify local extrema (minima and maxima) and compute their hop neighborhoods
 *          2. Classify extrema as either spurious (hop_idx <= hop_idx_thld) or non-spurious
 *          3. Apply the selected smoothing method to each spurious extremum
 *          4. Recompute extrema after smoothing to check if:
 *             - No new spurious extrema were introduced
 *             - All non-spurious extrema were preserved
 *             - Any remaining spurious extrema need further processing
 *          5. Repeat until convergence criteria are met or max iterations reached
 *
 * @param[in] y Vector of function values at each vertex
 * @param[in] hop_idx_thld Threshold for hop index; extrema with hop_idx <= this value
 *                         are considered spurious and will be smoothed
 * @param[in] smoother_type The type of smoothing method to use (weighted mean, harmonic, etc.)
 * @param[in] max_outer_iterations Maximum number of global smoothing iterations
 * @param[in] max_inner_iterations Maximum number of iterations for each smoothing operation
 * @param[in] smoothing_tolerance Convergence tolerance for smoothing algorithms
 * @param[in] sigma Parameter controlling the smoothing kernel width for weighted mean method
 * @param[in] process_in_order Whether to process extrema in ascending order of hop index
 * @param[in] verbose Whether to print progress information
 * @param[in] detailed_recording Whether to record detailed information about each smoothing step
 *
 * @return gflow_cx_t A structure containing:
 *         - harmonic_predictions: Smoothed function values
 *         - lmin_hop_nbhd_map: Map of local minima to their hop neighborhoods
 *         - lmax_hop_nbhd_map: Map of local maxima to their hop neighborhoods
 *         - smoothing_history: (If detailed_recording=TRUE) Records of each smoothing step
 *
 * @note The iterative nature of this algorithm helps ensure that smoothing operations don't
 *       introduce new unwanted extrema or eliminate important (non-spurious) ones. At each
 *       iteration, all current spurious extrema are processed in a single pass before
 *       validation and the next iteration.
 *
 * @see compute_extrema_hop_nbhds(), perform_weighted_mean_hop_disk_extension(),
 *      harmonic_extender(), harmonic_extension_eigen(),
 *      hybrid_biharmonic_harmonic_extension(), boundary_smoothed_harmonic_extension()
 */
gflow_cx_t set_wgraph_t::create_gflow_cx(
    const std::vector<double>& y,
    size_t hop_idx_thld,
    smoother_type_t smoother_type,
    int max_outer_iterations,
    int max_inner_iterations,
    double smoothing_tolerance,
    double sigma,
    bool process_in_order,
    bool verbose,
    bool detailed_recording
) const {

    // Initialize the smoothed function with original values
    std::vector<double> smoothed_y = y;

    struct extremum_info_t {
        size_t vertex;
        size_t hop_idx;
        bool is_minimum;
    };

    // Initialize gflow_cx with initial data
    gflow_cx_t gflow_cx;
    std::vector<smoothing_step_t> smoothing_history;

    //------------------------------------------------------------------
    // PHASE 1: Initial Analysis - Find extrema
    //------------------------------------------------------------------
    auto [lmin_hop_nbhd_map, lmax_hop_nbhd_map] = compute_extrema_hop_nbhds(y);

    // Store the initial extrema maps
    gflow_cx.lmin_hop_nbhd_map = lmin_hop_nbhd_map;
    gflow_cx.lmax_hop_nbhd_map = lmax_hop_nbhd_map;

    //------------------------------------------------------------------
    // PHASE 2: Iterative Smoothing
    //------------------------------------------------------------------

    // Track extrema that should be preserved (hop_idx > hop_idx_thld)
    std::unordered_set<size_t> non_spurious_lmin;
    std::unordered_set<size_t> non_spurious_lmax;

    for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            non_spurious_lmin.insert(vertex);
        }
    }

    for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            non_spurious_lmax.insert(vertex);
        }
    }

    // Iterative smoothing process
    bool smoothing_complete = false;
    int outer_iteration = 0;

    while (!smoothing_complete && outer_iteration < max_outer_iterations) {
        if (verbose) {
            Rprintf("Outer smoothing iteration %d/%d\n", outer_iteration + 1, max_outer_iterations);
        }

        // Collect extrema to process in this iteration
        std::vector<extremum_info_t> extrema_to_process;

        // Collect hypothetically spurious minima
        for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
            if (hop_nbhd.hop_idx <= hop_idx_thld) {
                extrema_to_process.push_back({vertex, hop_nbhd.hop_idx, true});
            }
        }

        // Collect hypothetically spurious maxima
        for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
            if (hop_nbhd.hop_idx <= hop_idx_thld) {
                extrema_to_process.push_back({vertex, hop_nbhd.hop_idx, false});
            }
        }

        // If no spurious extrema remain, we're done
        if (extrema_to_process.empty()) {
            smoothing_complete = true;
            if (verbose) {
                Rprintf("No spurious extrema remain after %d iterations, smoothing complete\n",
                        outer_iteration);
            }
            break;
        }

        // Sort extrema by hop_idx if requested
        if (process_in_order) {
            std::sort(extrema_to_process.begin(), extrema_to_process.end(),
                     [](const extremum_info_t& a, const extremum_info_t& b) {
                         return a.hop_idx < b.hop_idx;
                     });
        }

        // Process each extremum
        for (const auto& extremum : extrema_to_process) {
            size_t vertex = extremum.vertex;
            bool is_minimum = extremum.is_minimum;

            // Get the appropriate hop neighborhood
            const hop_nbhd_t& hop_nbhd = is_minimum ?
                lmin_hop_nbhd_map[vertex] : lmax_hop_nbhd_map[vertex];

            // Record pre-smoothing state if detailed recording is enabled
            smoothing_step_t step;
            if (detailed_recording) {
                step.vertex     = vertex;
                step.is_minimum = is_minimum;
                step.hop_idx    = hop_nbhd.hop_idx;
                step.smoother   = smoother_type;
                step.before     = smoothed_y;
            }

            // Get interior vertices (hop distance <= hop_radius)
            std::unordered_set<size_t> region_vertices;
            for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
                if (dist <= hop_nbhd.hop_idx) {
                    region_vertices.insert(v);
                }
            }

            // Get boundary vertices (hop distance = hop_radius + 1)
            std::unordered_map<size_t, double> boundary_values;
            for (const auto& [v, val] : hop_nbhd.y_nbhd_bd_map) {
                boundary_values[v] = smoothed_y[v];
            }

            if (detailed_recording) {
                step.region          = region_vertices;
                step.boundary_values = boundary_values;
            }

            // Apply the appropriate smoothing method
            switch (smoother_type) {
                case smoother_type_t::WMEAN:
                    perform_weighted_mean_hop_disk_extension(
                        smoothed_y,
                        hop_nbhd,
                        max_inner_iterations,
                        smoothing_tolerance,
                        sigma,
                        verbose
                    );
                    break;

                case smoother_type_t::HARMONIC_IT:
                    {
                        int record_frequency = 10;
                        harmonic_extender_t res = harmonic_extender(
                            boundary_values,
                            region_vertices,
                            max_inner_iterations,
                            smoothing_tolerance,
                            record_frequency,
                            verbose
                        );

                        // Update smoothed_y with the final iteration results
                        if (!res.iterations.empty()) {
                            auto& final_iter = res.iterations.back();
                            for (size_t v : region_vertices) {
                                smoothed_y[v] = final_iter[v];
                            }
                        }
                    }
                    break;

                case smoother_type_t::HARMONIC_EIGEN:
                    {
                        double regularization = 1e-10;
                        auto harmonic_result = harmonic_extension_eigen(
                            boundary_values,
                            region_vertices,
                            regularization,
                            verbose
                        );

                        // Update smoothed_y with harmonic extension result
                        for (size_t v : region_vertices) {
                            smoothed_y[v] = harmonic_result[v];
                        }
                    }
                    break;

                case smoother_type_t::HYBRID_BIHARMONIC_HARMONIC:
                    {
                        int boundary_blend_distance = 2;
                        auto hybrid_result = hybrid_biharmonic_harmonic_extension(
                            boundary_values,
                            region_vertices,
                            boundary_blend_distance,
                            verbose
                        );

                        // Update smoothed_y with hybrid extension result
                        for (size_t v : region_vertices) {
                            smoothed_y[v] = hybrid_result[v];
                        }
                    }
                    break;

                case smoother_type_t::BOUNDARY_SMOOTHED_HARMONIC:
                    {
                        int boundary_blend_distance = 2;
                        auto boundary_smoothed_result = boundary_smoothed_harmonic_extension(
                            boundary_values,
                            region_vertices,
                            boundary_blend_distance,
                            verbose
                        );

                        // Update smoothed_y with boundary smoothed extension result
                        for (size_t v : region_vertices) {
                            smoothed_y[v] = boundary_smoothed_result[v];
                        }
                    }
                    break;
            }

            // Record post-smoothing state if detailed recording is enabled
            if (detailed_recording) {
                step.after = smoothed_y;
                smoothing_history.push_back(step);
            }
        }

        // Recompute extrema for validation
        auto [new_lmin_hop_nbhd_map, new_lmax_hop_nbhd_map] = compute_extrema_hop_nbhds(smoothed_y);

        // Check if validation criteria are met:
        // 1. No new spurious extrema
        // 2. All non-spurious extrema preserved
        bool validation_passed = true;

        // Check for new spurious minima
        for (const auto& [vertex, hop_nbhd] : new_lmin_hop_nbhd_map) {
            if (hop_nbhd.hop_idx <= hop_idx_thld) {
                // If this is a new spurious minimum, or if it was previously non-spurious
                if (lmin_hop_nbhd_map.find(vertex) == lmin_hop_nbhd_map.end() ||
                    non_spurious_lmin.find(vertex) != non_spurious_lmin.end()) {
                    if (verbose) {
                        Rprintf("Validation failed: Found spurious minimum at vertex %zu\n", vertex);
                    }
                    validation_passed = false;
                    break;
                }
            }
        }

        // Check that all non-spurious minima are preserved
        if (validation_passed) {
            for (size_t vertex : non_spurious_lmin) {
                if (new_lmin_hop_nbhd_map.find(vertex) == new_lmin_hop_nbhd_map.end()) {
                    if (verbose) {
                        Rprintf("Validation failed: Lost non-spurious minimum at vertex %zu\n", vertex);
                    }
                    validation_passed = false;
                    break;
                }
            }
        }

        // Analogous checks for maxima
        if (validation_passed) {
            for (const auto& [vertex, hop_nbhd] : new_lmax_hop_nbhd_map) {
                if (hop_nbhd.hop_idx <= hop_idx_thld) {
                    if (lmax_hop_nbhd_map.find(vertex) == lmax_hop_nbhd_map.end() ||
                        non_spurious_lmax.find(vertex) != non_spurious_lmax.end()) {
                        if (verbose) {
                            Rprintf("Validation failed: Found spurious maximum at vertex %zu\n", vertex);
                        }
                        validation_passed = false;
                        break;
                    }
                }
            }
        }

        if (validation_passed) {
            for (size_t vertex : non_spurious_lmax) {
                if (new_lmax_hop_nbhd_map.find(vertex) == new_lmax_hop_nbhd_map.end()) {
                    if (verbose) {
                        Rprintf("Validation failed: Lost non-spurious maximum at vertex %zu\n", vertex);
                    }
                    validation_passed = false;
                    break;
                }
            }
        }

        // Update extrema maps for next iteration
        lmin_hop_nbhd_map = new_lmin_hop_nbhd_map;
        lmax_hop_nbhd_map = new_lmax_hop_nbhd_map;

        // Check if we're done
        if (validation_passed) {
            // No new spurious extrema and all non-spurious extrema preserved
            bool any_spurious_extrema = false;

            for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
                if (hop_nbhd.hop_idx <= hop_idx_thld &&
                    non_spurious_lmin.find(vertex) == non_spurious_lmin.end()) {
                    any_spurious_extrema = true;
                    break;
                }
            }

            if (!any_spurious_extrema) {
                for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
                    if (hop_nbhd.hop_idx <= hop_idx_thld &&
                        non_spurious_lmax.find(vertex) == non_spurious_lmax.end()) {
                        any_spurious_extrema = true;
                        break;
                    }
                }
            }

            if (!any_spurious_extrema) {
                smoothing_complete = true;
                if (verbose) {
                    Rprintf("All spurious extrema have been eliminated after %d iterations\n",
                            outer_iteration + 1);
                }
            }
        }

        outer_iteration++;
    }

    // Store final results
    gflow_cx.lmin_hop_nbhd_map = lmin_hop_nbhd_map;
    gflow_cx.lmax_hop_nbhd_map = lmax_hop_nbhd_map;
    gflow_cx.harmonic_predictions = smoothed_y;

    // If detailed recording was requested, store the history
    if (detailed_recording) {
        gflow_cx.smoothing_history = smoothing_history;
    }

    return gflow_cx;
}
