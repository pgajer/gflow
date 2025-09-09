	#if 0
	gradient_flow_t::trajectory_t construct_trajectory_with_distance_sorting(
		size_t start_vertex,
		const std::vector<double>& scale,
		const std::vector<double>& y) const;
	#endif

	std::pair<std::vector<size_t>, std::vector<size_t>> construct_both_trajectories(
		size_t start_vertex,
		const std::vector<double>& scale,
		const std::vector<double>& y) const;
#if 0
	std::vector<size_t> construct_trajectory(
		size_t start,
		bool ascending,
		double radius,
		const std::vector<double>& y) const;

	std::vector<size_t> construct_trajectory(
		size_t start,
		bool ascending,
		const std::vector<double>& scale,
		const std::vector<double>& y) const;
	#endif

/**
 * @brief Simultaneously constructs ascending and descending trajectories from a starting vertex
 *
 * This function efficiently computes both ascending and descending gradient flow trajectories
 * from a given vertex in a single pass. By computing the neighborhood just once, it reduces
 * the computational overhead compared to separate trajectory calculations.
 *
 * The function handles flat regions (where function values are equal) by using topological
 * information from the graph structure. For chain graphs, it follows the natural ordering of
 * vertices when function values are equal.
 *
 * @param start_vertex Starting vertex for both trajectories
 * @param scale Vector of scale parameters for each vertex, determining neighborhood radius
 * @param y Vector of function values at each vertex
 * @return A pair of vectors containing {ascending_trajectory, descending_trajectory}
 * @see find_shortest_paths_within_radius
 * @see compute_gradient_flow
 */
std::pair<std::vector<size_t>, std::vector<size_t>> set_wgraph_t::construct_both_trajectories(
    size_t start_vertex,
    const std::vector<double>& scale,
    const std::vector<double>& y) const {

    std::vector<size_t> ascending = {start_vertex};
    std::vector<size_t> descending = {start_vertex};
    std::set<size_t> visited_asc = {start_vertex};
    std::set<size_t> visited_desc = {start_vertex};

    // Compute neighborhood just once
    shortest_paths_t shortest_paths = find_graph_paths_within_radius(start_vertex, scale[start_vertex]);
    //auto shortest_paths = find_shortest_paths_within_radius(start_vertex, scale[start_vertex]);

    // Process ascending trajectory
    size_t current_asc = start_vertex;
    while (true) {
        double max_value = y[current_asc];
        size_t next_asc = INVALID_VERTEX;

        for (const auto& [neighbor, _] : adjacency_list[current_asc]) {
            if (!visited_asc.count(neighbor) &&
                shortest_paths.reachable_vertices.count(neighbor) &&
                y[neighbor] > max_value) {
                max_value = y[neighbor];
                next_asc = neighbor;
            }
        }

        if (next_asc == INVALID_VERTEX) break;

        ascending.push_back(next_asc);
        visited_asc.insert(next_asc);
        current_asc = next_asc;
    }

    // Process descending trajectory (similar logic)
    size_t current_desc = start_vertex;
    while (true) {
        double min_value = y[current_desc];
        size_t next_desc = INVALID_VERTEX;

        for (const auto& [neighbor, _] : adjacency_list[current_desc]) {
            if (!visited_desc.count(neighbor) &&
                shortest_paths.reachable_vertices.count(neighbor) &&
                y[neighbor] < min_value) {
                min_value = y[neighbor];
                next_desc = neighbor;
            }
        }

        if (next_desc == INVALID_VERTEX) break;

        descending.push_back(next_desc);
        visited_desc.insert(next_desc);
        current_desc = next_desc;
    }

    return {ascending, descending};
}

/**
 * @brief Constructs a gradient flow trajectory with distance-based sorting for flat regions
 *
 * This function creates a complete gradient flow trajectory that addresses the challenge of
 * flat regions (equal function values) by incorporating graph distance information. The
 * algorithm works in several steps:
 *
 * 1. Constructs preliminary ascending and descending trajectories
 * 2. Identifies the furthest vertex from the starting point
 * 3. Determines whether this furthest vertex is a local minimum or maximum
 * 4. Computes distances from this furthest vertex to all trajectory vertices
 * 5. Sorts vertices based on both function values and distances to create a consistent flow
 * 6. Calculates the total weight of the resulting trajectory
 *
 * This approach is particularly effective for chain graphs with regions of equal values,
 * as it creates a natural flow that follows the graph structure while respecting the
 * function gradient where it exists.
 *
 * @param start_vertex Starting vertex for the trajectory
 * @param scale Vector of scale parameters for each vertex, determining neighborhood radius
 * @param y Vector of function values at each vertex
 * @return A gradient_flow_t::trajectory_t structure containing ordered vertices and total path weight
 * @see construct_both_trajectories
 * @see compute_distances_from_vertex
 * @see compute_graph_distance
 */
#if 0
gradient_flow_t::trajectory_t set_wgraph_t::construct_trajectory_with_distance_sorting(
    size_t start_vertex,
    const std::vector<double>& scale,
    const std::vector<double>& y) const {

    // Step 1: Construct preliminary trajectories
    auto [ascending, descending] = construct_both_trajectories(start_vertex, scale, y);

    // Step 2: Find furthest vertex (could be from ascending or descending)
    size_t furthest_vertex;
    double max_distance = 0;

    // Check ascending path for furthest vertex
    for (size_t v : ascending) {
        double dist = compute_graph_distance(start_vertex, v);
        if (dist > max_distance) {
            max_distance = dist;
            furthest_vertex = v;
        }
    }

    // Check descending path for furthest vertex
    for (size_t v : descending) {
        double dist = compute_graph_distance(start_vertex, v);
        if (dist > max_distance) {
            max_distance = dist;
            furthest_vertex = v;
        }
    }

    // Step 3: Determine if furthest vertex is minimum or maximum
    bool furthest_is_minimum = y[furthest_vertex] < y[start_vertex];

    // Step 4: Compute distances from furthest vertex to all vertices
    auto distances_from_furthest = compute_distances_from_vertex(furthest_vertex);

    // Step 5: Create combined trajectory with all vertices
    std::vector<size_t> all_vertices;
    all_vertices.insert(all_vertices.end(), ascending.begin(), ascending.end());

    // Add descending vertices (exclude start to avoid duplication)
    for (size_t i = 1; i < descending.size(); ++i) {
        all_vertices.push_back(descending[i]);
    }

    // Step 6: Sort vertices based on distance from furthest vertex
    // and the nature of furthest vertex (min/max)
    std::sort(all_vertices.begin(), all_vertices.end(),
        [&](size_t a, size_t b) {
            if (std::abs(y[a] - y[b]) < 1e-10) {
                // If function values are equal, use distance as tiebreaker
                if (furthest_is_minimum) {
                    // Ascending order from minimum
                    return distances_from_furthest[a] < distances_from_furthest[b];
                } else {
                    // Descending order from maximum
                    return distances_from_furthest[a] > distances_from_furthest[b];
                }
            } else {
                // Normal case - use function values
                return y[a] < y[b];
            }
        });

    // Create trajectory
    gradient_flow_t::trajectory_t trajectory;
    trajectory.vertices = all_vertices;

    // Calculate total weight of the trajectory
    double total_weight = 0.0;
    for (size_t i = 0; i < all_vertices.size() - 1; ++i) {
        size_t v1 = all_vertices[i];
        size_t v2 = all_vertices[i + 1];

        // Find edge weight in adjacency list
        bool edge_found = false;
        for (const auto& edge : adjacency_list[v1]) {
            if (edge.vertex == v2) {
                total_weight += edge.weight;
                edge_found = true;
                break;
            }
        }

        // If vertices are not directly connected (which might happen after sorting),
        // we need to find the shortest path between them
        if (!edge_found) {
            // Get shortest path between these vertices
            auto path_result = find_shortest_path(v1, v2);

            // Add weights of all edges in the shortest path
            for (size_t j = 0; j < path_result.path.size() - 1; ++j) {
                size_t path_v1 = path_result.path[j];
                size_t path_v2 = path_result.path[j + 1];

                for (const auto& edge : adjacency_list[path_v1]) {
                    if (edge.vertex == path_v2) {
                        total_weight += edge.weight;
                        break;
                    }
                }
            }
        }
    }

    trajectory.total_weight = total_weight;

    return trajectory;
}
#endif

/**
 * @brief Constructs a trajectory following the gradient from a starting vertex
 *
 * Builds a path that follows either ascending or descending gradient of the values
 * in vector y, constrained to neighbors within the specified radius. The trajectory
 * continues until reaching a local extremum.
 *
 * @param start Starting vertex for the trajectory
 * @param ascending If true, follows increasing values; if false, follows decreasing values
 * @param radius Maximum distance to consider for neighboring vertices
 * @param y Vector of values associated with vertices
 *
 * @return std::vector<size_t> Sequence of vertices forming the trajectory
 *
 * @pre start must be a valid vertex index
 * @pre radius must be positive
 * @pre y must not contain duplicate values
 * @pre y.size() must match the graph size
 */
std::vector<size_t> set_wgraph_t::construct_trajectory(
    size_t start,
    bool ascending,
    double radius,
    const std::vector<double>& y) const {

    std::vector<size_t> trajectory{start};
    size_t current = start;


    Rprintf("In construct_trajectory(double radius) start: %zu\tascending: %d\n", start, (int)ascending);

    while (true) {
        // Compute reachability map from current vertex with radius based on scale
        reachability_map_t reachability_map = compute_graph_reachability_map(current, radius);
        std::vector<vertex_shortest_path_info_t> path_endpoints = find_graph_path_endpoints(reachability_map);

        {
            Rprintf("In construct_trajectory() start: %zu\tascending: %d\n", start, (int)ascending);
            Rprintf("path_endpoints\n");
            for (const auto& p : path_endpoints) {
                Rprintf("\nvertex: %zu\tdistance: %.3f\n", p.vertex, p.distance);
                print_vect(p.path,"path");
                for (const auto& v : p.path) {
                    Rprintf("v: %zu\ty[v]: %.3f\n", v, y[v]);
                }
                Rprintf("\n");
            }
        }

        // If no endpoints found, gradually increase radius until we find some
        if (path_endpoints.empty()) {
            REPORT_WARNING("path_endpoints.size() == 0: enlarging the radius\n");
            double multiplier = 2.0;
            double radius = radius;
            while(path_endpoints.empty()) {
                radius *= multiplier;
                reachability_map = compute_graph_reachability_map(current, radius);
                path_endpoints = find_graph_path_endpoints(reachability_map);
            }
        }

        // Find the path endpoint with the optimal rate of change
        size_t next = current;  // Default to current if no better vertex found
        std::vector<size_t> best_path;  // Store the path to the best next vertex

        if (ascending) {

            Rprintf("\nIn ascending\n");

            // For ascending trajectory, find endpoint with highest positive rate of change
            double max_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                // Calculate rate of change = change in y value / distance
                double rate = (y[v.vertex] - y[current]) / v.distance;


                Rprintf("v.vertex: %zu\tdrate: %.3f\t(rate > max_rate_of_change): %d\n", v.vertex, rate, (int)(rate > max_rate_of_change));

                if (rate > max_rate_of_change) {
                    max_rate_of_change = rate;
                    next = v.vertex;
                    best_path = std::move(v.path);
                }
            }
        } else {

            Rprintf("\nIn descending\n");

            // For descending trajectory, find endpoint with lowest (most negative) rate of change
            double min_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                double rate = (y[v.vertex] - y[current]) / v.distance;

                    Rprintf("v.vertex: %zu\tdrate: %.3f\t(rate < min_rate_of_change): %d\n", v.vertex, rate, (int)(rate < min_rate_of_change));

                if (rate < min_rate_of_change) {
                    min_rate_of_change = rate;
                    next = v.vertex;
                    best_path = std::move(v.path);
                }
            }
        }

        // If no better vertex found, we've reached an extremum
        if (next == current) break;

        // Since current is the first vertex in the path by construction,
        // we can directly add vertices from index 1 onward to avoid duplication
        if (best_path.size() > 1) {
            trajectory.insert(trajectory.end(), best_path.begin() + 1, best_path.end());
        } else {
            // This would only happen if the path contains only the current vertex
            // (unlikely but possible edge case)
            REPORT_WARNING("Path has only the current vertex, no progress possible");
        }

        current = next;  // Move to the next vertex for the next iteration
    }



    #if 0 // old
    while (true) {
        shortest_paths_t shortest_paths = find_graph_paths_within_radius(current, radius);
        std::unordered_set<int> neighbors = shortest_paths.reachable_vertices;

        // Find next vertex (maximum for ascending, minimum for descending)
        size_t next = current;  // Default to current if no better vertex found
        for (size_t neighbor : neighbors) {
            if (ascending && y[neighbor] > y[next]) {
                next = neighbor;
            } else if (!ascending && y[neighbor] < y[next]) {
                next = neighbor;
            }
        }

        // If no better vertex found, we've reached an extremum
        if (next == current) break;

        current = next;
        trajectory.push_back(current);
    }
    #endif

    return trajectory;
}

/**
 * @brief Constructs a trajectory by following function gradients
 *
 * This function builds a path through the graph by repeatedly moving to adjacent vertices
 * that increase (or decrease) the value of a function defined by the vector y.
 * It uses locally adaptive search radii based on the scale vector to determine
 * the neighborhood for each step.
 *
 * The function continues until it reaches a local extremum - a vertex where no
 * neighbors have a better (higher for ascending, lower for descending) function value.
 *
 * @param start The index of the starting vertex for the trajectory
 * @param ascending If true, move toward increasing y values; if false, toward decreasing values
 * @param scale Vector of adaptive radius values for each vertex
 * @param y Vector of function values at each vertex
 *
 * @return A vector of vertex indices representing the trajectory through the graph
 *
 * @note If no neighbors are found within the initial radius (scale[current]),
 *       the function automatically enlarges the search radius until neighbors are found
 *
 * @note The trajectory follows the steepest ascent/descent path according to the
 *       function values in y
 *
 * @note The trajectory is constructed by connecting shortest paths between successive
 *       chosen vertices
 *
 * @see find_shortest_paths_within_radius
 */
#if 0
std::vector<size_t> set_wgraph_t::construct_trajectory(
    size_t start,
    bool ascending,
    const std::vector<double>& scale,
    const std::vector<double>& y) const {

    std::vector<size_t> trajectory{start};
    size_t current = start;

    Rprintf("\n----------------------------------------\n"
            "In construct_trajectory(std::vector<double>& scale) start: %zu\tascending: %d\n", start, (int)ascending);

    while (true) {
        // Compute reachability map from current vertex with radius based on scale
        reachability_map_t reachability_map = compute_graph_reachability_map(current, scale[current]);
        std::vector<vertex_shortest_path_info_t> path_endpoints = find_graph_path_endpoints(reachability_map);

        {
            Rprintf("in while(true)\tstart: %zu\tascending: %d\n", start, (int)ascending);
            Rprintf("path_endpoints\n");
            for (const auto& p : path_endpoints) {
                Rprintf("\nendpt: %zu\tdistance: %.3f\n", p.vertex, p.distance);
                print_vect(p.path,"path");
                for (const auto& v : p.path) {
                    Rprintf("y[%zu]: %.5f\n", v, y[v]);
                }
                Rprintf("\n");
            }
        }

        // If no endpoints found, gradually increase radius until we find some
        if (path_endpoints.empty()) {
            REPORT_WARNING("path_endpoints.size() == 0: enlarging the radius\n");
            double multiplier = 2.0;
            double radius = scale[current];
            while(path_endpoints.empty()) {
                radius *= multiplier;
                reachability_map = compute_graph_reachability_map(current, radius);
                path_endpoints = find_graph_path_endpoints(reachability_map);
            }
        }

        // Find the path endpoint with the optimal rate of change
        size_t next = current;  // Default to current if no better vertex found
        std::vector<size_t> best_path;  // Store the path to the best next vertex

        if (ascending) {

            Rprintf("\nIn if(ascending) block\n");

            // For ascending trajectory, find endpoint with highest positive rate of change
            double max_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                // Calculate rate of change = change in y value / distance
                double rate = (y[v.vertex] - y[current]) / v.distance;

                    Rprintf("v.vertex: %zu\trate: %.3f\t(rate > max_rate_of_change): %d\n", v.vertex, rate, (int)(rate > max_rate_of_change));

                if (rate > max_rate_of_change) {
                    max_rate_of_change = rate;
                    next = v.vertex;
                    best_path = std::move(v.path);
                }
            }
        } else {

            Rprintf("\nIn if(descending) block\n");

            // For descending trajectory, find endpoint with lowest (most negative) rate of change
            double min_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                double rate = (y[v.vertex] - y[current]) / v.distance;


                Rprintf("v.vertex: %zu\trate: %.3f\t(rate < min_rate_of_change): %d\n", v.vertex, rate, (int)(rate < min_rate_of_change));

                if (rate < min_rate_of_change) {
                    min_rate_of_change = rate;
                    next = v.vertex;
                    best_path = std::move(v.path);
                }
            }
        }

        // If no better vertex found, we've reached an extremum
        if (next == current) break;

        // Since current is the first vertex in the path by construction,
        // we can directly add vertices from index 1 onward to avoid duplication
        if (best_path.size() > 1) {
            trajectory.insert(trajectory.end(), best_path.begin() + 1, best_path.end());
        } else {
            // This would only happen if the path contains only the current vertex
            // (unlikely but possible edge case)
            REPORT_WARNING("Path has only the current vertex, no progress possible");
        }

        current = next;  // Move to the next vertex for the next iteration
    }

    return trajectory;
}
#endif


/**
 * @brief Constructs a trajectory by following function gradients
 *
 * This function builds a path through the graph by repeatedly moving to adjacent vertices
 * that increase (or decrease) the value of a function defined by the vector y.
 * It uses locally adaptive search radii based on the scale vector to determine
 * the neighborhood for each step.
 *
 * The function continues until it reaches a local extremum - a vertex where no
 * neighbors have a better (higher for ascending, lower for descending) function value.
 *
 * @param start The index of the starting vertex for the trajectory
 * @param ascending If true, move toward increasing y values; if false, toward decreasing values
 * @param scale Vector of adaptive radius values for each vertex
 * @param y Vector of function values at each vertex
 *
 * @return A vector of vertex indices representing the trajectory through the graph
 *
 * @note If no neighbors are found within the initial radius (scale[current]),
 *       the function automatically enlarges the search radius until neighbors are found
 *
 * @note The trajectory follows the steepest ascent/descent path according to the
 *       function values in y
 *
 * @note The trajectory is constructed by connecting shortest paths between successive
 *       chosen vertices
 *
 * @see find_shortest_paths_within_radius
 */
std::vector<size_t> set_wgraph_t::construct_trajectory(
    size_t start,
    bool ascending,
    const std::vector<double>& scale,
    const std::vector<double>& y) const {

    std::vector<size_t> trajectory{start};
    size_t current = start;

    while (true) {

        // shortest_paths_t shortest_paths = find_graph_paths_within_radius(current, scale[current]);
        // std::unordered_set<int> neighbors = shortest_paths.reachable_vertices;

        reachability_map_t reachability_map = compute_graph_reachability_map(current, scale[current]);
        std::vector<vertex_info_t> path_endpoints = find_graph_path_endpoints(reachability_map);

        if (neighbors.size() == 0) {
            REPORT_WARNING("neighbors.size() == 0: enlarging the radius\n");
            double multiplier = 2.0;
            double radius = scale[current];
            while(neighbors.size() == 0) {
                radius *= multiplier;
                shortest_paths = find_graph_paths_within_radius(current, radius);
                neighbors = shortest_paths.reachable_vertices;
            }
        }

        // Find next vertex (maximum for ascending, minimum for descending)
        size_t next = current;  // Default to current if no better vertex found
        for (size_t neighbor : neighbors) {
            if (ascending && y[neighbor] > y[next]) {
                next = neighbor;
            } else if (!ascending && y[neighbor] < y[next]) {
                next = neighbor;
            }
        }

        // If no better vertex found, we've reached an extremum
        if (next == current) break;

        // Extract the shortest path from 'current' to 'next'
        auto subpath_info = shortest_paths.vertex_to_path_map[next];
        std::vector<size_t> subpath(
            shortest_paths.paths[subpath_info.path_idx].vertices.begin() + 1, // current is already in trajectory so we need to start from the next vertex
            shortest_paths.paths[subpath_info.path_idx].vertices.begin() + subpath_info.vertex_idx + 1
            );

        trajectory.insert(trajectory.end(), subpath.begin(), subpath.end());
        current = next;
    }

    return trajectory;
}
