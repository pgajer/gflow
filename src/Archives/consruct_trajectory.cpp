

gradient_trajectory_t set_wgraph_t::construct_trajectory(
    size_t start,
    bool ascending,
    const std::vector<double>& scale,
    const std::vector<double>& y,
    const std::unordered_map<size_t, bool>& extrema_map,
    double long_edge_lower_thld,
    double long_edge_upper_thld
    ) const {

#define DEBUG__construct_trajectory 0
#if DEBUG__construct_trajectory
    Rprintf("\n----------------------------------------\n"
            "In construct_trajectory(std::vector<double>& scale) start: %zu\n", start + 1);
    Rprintf("extrema_map\nvertex\tis_lmax\n");
    for (const auto& [v, is_lmax] : extrema_map) {
        Rprintf("%zu\t%d\n", v + 1, (int)is_lmax);
    }
    Rprintf("\n");
#endif

    // Initialize result structure
    gradient_trajectory_t result;
    result.path.push_back(start);
    result.ends_at_critical = false;  // Whether the trajectory ends at local maximum or minimum
    result.ends_at_lmax     = false;  // Whether the trajectory ends at local maximum

    size_t current = start;
    bool should_continue = true;

    // Lambda to filter out paths that go through visited vertices
    auto filter_out_visited_vertices = [&visited_vertices](const std::vector<vertex_path_t>& paths) {
        std::vector<vertex_path_t> filtered_paths;
        for (const auto& path : paths) {
            // Check if any intermediate vertex in the path is already visited
            bool has_visited = false;
            for (size_t i = 1; i < path.path.size(); ++i) {
                if (visited_vertices.count(path.vertices[i]) > 0) {
                    has_visited = true;
                    break;
                }
            }
            if (!has_visited) {
                filtered_paths.push_back(path);
            }
        }

        return filtered_paths;
    };

    // Main trajectory building loop
    while (should_continue) {
        // Compute reachability map from current vertex with radius based on scale
        reachability_map_t reachability_map = compute_graph_reachability_map(current, scale[current]);
        // std::vector<vertex_shortest_path_info_t> path_endpoints = get_vertex_shortest_paths(reachability_map);  // <<--- this to be replaced by reconstruct_graph_paths()
        std::vector<vertex_path_t> paths = reconstruct_graph_paths(reachability_map);

        // If no endpoints found, gradually increase radius until we find some
        if (paths.empty()) {
            double multiplier = 2.0;
            double radius = scale[current];
            while (path_endpoints.empty()) {
                radius *= multiplier;
                reachability_map = compute_graph_reachability_map(current, radius);
                // path_endpoints = get_vertex_shortest_paths(reachability_map);
                paths = reconstruct_graph_paths(reachability_map);

                // Safety check to avoid infinite loops
                if (radius > 1000 * scale[current]) {
                    REPORT_WARNING("Failed to find path endpoints even with greatly enlarged radius");
                    return result;  // Return what we have so far
                }
            }
        }

        // Skip visited vertices to prevent cycles
        paths = filter_out_visited_vertices(paths);

        for (const auto& path : paths) {
            // For each path
            bool has_extremum = false;
            bool is_maximum = false;
            size_t extremum_index = 0;

            double best_rate = ascending ? 0.0 : 0.0; // Initialize to 0 for ascending, 0 for descending
            size_t best_rate_index = 0;
            size_t prev_vertex = current_vertex;
            double prev_y = y[current_vertex];
            double cumulative_distance = 0.0;

            // Check each vertex along the path
            for (size_t i = 1; i < path.vertices.size(); ++i) {
                size_t vertex = path.vertices[i];

                // Check if vertex is an extremum
                auto extremum_it = extrema_map.find(vertex);
                if (extremum_it != extrema_map.end()) {
                    has_extremum = true;
                    is_maximum = extremum_it->second;
                    extremum_index = i;
                    break;
                }

                // Calculate total function change
                double total_change = y[vertex] - y[current_vertex];

                // Calculate cumulative absolute changes
                double cumulative_absolute_changes = 0.0;
                double prev_y = y[current_vertex];
                for (size_t j = 1; j <= i; ++j) {
                    size_t w = path.vertices[j];
                    cumulative_absolute_changes += std::abs(y[w] - prev_y);
                    prev_y = y[w];
                }

                // Calculate monotonicity index
                //double mono_index = calculate_monotonicity_index(total_change, cumulative_absolute_changes);
                double mono_index = (cumulative_absolute_changes > 0) ? total_change / cumulative_absolute_changes : 0;

                // Determine the maximum value among all edge lengths after
                // applying a threshold function to each. This represents an L^∞
                // (L-infinity) norm analogue to path evenness, capturing the
                // worst-case scenario in path segment distribution.
                std::vector<double> edge_lengths(i);
                for (size_t j = 1; j <= i; ++j) {
                    edge_lengths[j-1] = path.dist_to_ref_vertex[j] - path.dist_to_ref_vertex[j-1];
                }
                // max_thldd_edge_length is between 0 and 1, with value 0
                // indicating no long edges along the path and 1, indicating
                // existence of a long edge within the path whose length is
                // above long_edge_upper_thld
                double max_thldd_edge_length = calculate_path_max_threshold(edge_lengths, long_edge_lower_thld, long_edge_upper_thld);

                // Calculate adjusted penalty factor based on evenness
                double penalty_factor = 2.0 - max_thldd_edge_length;  // Ranges from 1.0 (no long edges along the path) to 2.0 (there is at least one long edge)

                // Calculate penalized rate
                double adjusted_abs_rate = 0.0;
                if (path.dist_to_ref_vertex[i] > 0.0) {
                    adjusted_abs_rate = std::abs(total_change) / std::pow(path.dist_to_ref_vertex[i], penalty_factor);
                }

                // Calculate final quality metric
                double monotonicity_adjusted_rate = mono_index * adjusted_abs_rate;

                evaluated_paths.emplace_back(path, monotonicity_adjusted_rate, total_change);

                // For ascending direction, look for highest positive monotonicity_adjusted_rate
                if (ascending) {
                    if (monotonicity_adjusted_rate > best_rate) {
                        best_rate = monotonicity_adjusted_rate;
                        best_rate_index = i;
                    }
                }
                // For descending direction, look for lowest negative rate
                else {
                    if (monotonicity_adjusted_rate < best_rate) {
                        best_rate = monotonicity_adjusted_rate;
                        best_rate_index = i;
                    }
                }
            }

            if (!has_extremum) {
                paths_with_no_extrema.push_back(std::make_tuple(path, rate));
            } else if (is_maximum && !ascending_done) {
                paths_with_maxima.push_back(std::make_tuple(path, rate, extremum_index));
            } else if (!descending_done) {  // extremum is minimum
                paths_with_minima.push_back(std::make_tuple(path, rate, extremum_index));
            }


#if 0
        // old version

        // Search through endpoints for the optimal rate of change
        if (ascending) {
            // For ascending trajectory, find endpoint with highest positive rate of change
            double max_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                // Calculate rate of change = change in y value / distance
                double rate = (y[v.vertex] - y[current]) / v.distance;
                if (rate > max_rate_of_change) {
                    max_rate_of_change = rate;
                    next = v.vertex;
                    best_path = v.path;  // We don't use std::move since we may need to check multiple paths
                }
            }
        } else {
            // For descending trajectory, find endpoint with lowest (most negative) rate of change
            double min_rate_of_change = 0.0;
            for (const auto& v : path_endpoints) {
                double rate = (y[v.vertex] - y[current]) / v.distance;
                if (rate < min_rate_of_change) {
                    min_rate_of_change = rate;
                    next = v.vertex;
                    best_path = v.path;
                }
            }
        }

        // If no better vertex found, we've reached an extremum
        if (next == current || best_path.empty()) {
            break;
        }

        // Check if the best path contains any local extrema
        // We need to find the first extremum along the path
        size_t first_extremum_idx = best_path.size();  // Initialize to out-of-bounds

        for (size_t idx = 1; idx < best_path.size(); ++idx) {  // Start from 1 to skip the current vertex
            size_t vertex = best_path[idx];
            auto extremum_it = extrema_map.find(vertex);

            if (extremum_it != extrema_map.end()) {
                // Mark that we're ending at a critical point
                result.ends_at_critical = true;

                // Truncate path at this critical point
                first_extremum_idx = idx;
                result.ends_at_lmax = extremum_it->second;

                // Always terminate upon finding a critical point
                should_continue = false;
                break;
            }
        }

         // Truncate the path if needed
        if (first_extremum_idx < best_path.size()) {
            // Add the truncated path (up to and including the extremum)
            result.path.insert(result.path.end(), best_path.begin() + 1, best_path.begin() + first_extremum_idx + 1);

            // Update current to the extremum for next iteration
            current = best_path[first_extremum_idx];

        } else {
            // No extrema on the path, add the entire path segment
            result.path.insert(result.path.end(), best_path.begin() + 1, best_path.end());
            current = next;
        }
#endif

    } // EDN OF while (should_continue) {

    return result;
}
