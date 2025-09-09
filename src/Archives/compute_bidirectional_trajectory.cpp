/**
 * @brief Computes bidirectional gradient trajectories from a vertex, finding both ascending and descending paths
 *
 * This function simultaneously computes both ascending and descending gradient trajectories
 * starting from a given vertex. It categorizes paths based on whether they contain extrema
 * (local minima or maxima) and prioritizes finding complete trajectories that connect
 * local minima to local maxima.
 *
 * The algorithm handles several cases:
 * 1. When both local minima and maxima are found, it completes both trajectories
 * 2. When only one type of extrema is found, it completes that trajectory and continues searching for the other
 * 3. When no extrema are found, it follows the steepest gradient direction
 *
 * Cycles are prevented by tracking visited vertices, and paths are selected based on
 * the steepness of their gradient (rate of change).
 *
 * @param start The starting vertex index
 * @param scale Vector of scale values for each vertex, used to determine search radius
 * @param y Vector of function values at each vertex
 * @param extrema_map Map of vertex indices to boolean indicating whether vertex is a maximum (true) or minimum (false)
 * @return A pair of gradient trajectories {ascending_trajectory, descending_trajectory}
 */
std::pair<gradient_trajectory_t, gradient_trajectory_t>
set_wgraph_t::compute_bidirectional_trajectory(
    size_t start,
    const std::vector<double>& scale,
    const std::vector<double>& y,
    const std::unordered_map<size_t, bool>& extrema_map,
    double long_edge_lower_thld,
    double long_edge_upper_thld
    ) const {

#define DEBUG__compute_bidirectional_trajectory 0

    // Initialize trajectories
    gradient_trajectory_t ascending_traj;
    gradient_trajectory_t descending_traj;
    ascending_traj.path.push_back(start);
    descending_traj.path.push_back(start);

    // Initialize state variables
    size_t current_vertex = start;
    std::unordered_set<size_t> visited_vertices;
    visited_vertices.insert(start);
    bool ascending_done = false;
    bool descending_done = false;

    // Lambda to filter out paths that go through visited vertices
    auto filter_out_visited_vertices = [&visited_vertices](const std::vector<vertex_shortest_path_info_t>& paths) {
        std::vector<vertex_shortest_path_info_t> filtered_paths;
        for (const auto& path : paths) {
            if (visited_vertices.count(path.vertex) == 0) {
                // Check if any intermediate vertex in the path is already visited
                bool has_visited = false;
                for (size_t i = 1; i < path.path.size() - 1; ++i) {
                    if (visited_vertices.count(path.path[i]) > 0) {
                        has_visited = true;
                        break;
                    }
                }
                if (!has_visited) {
                    filtered_paths.push_back(path);
                }
            }
        }
        return filtered_paths;
    };


#if DEBUG__compute_bidirectional_trajectory
    // ToDo <<----

    // How do we deal with sub-paths???

    // 5. Replace the filter_by_positive_rate and filter_by_negative_rate functions
    // with functions that filter by quality metric sign:

    auto filter_by_positive_quality = [](
        const std::vector<std::tuple<vertex_shortest_path_info_t, double, double>>& paths) {

        std::vector<std::tuple<vertex_shortest_path_info_t, double, double>> filtered_paths;
        for (const auto& [path, quality, change] : paths) {
            if (quality > 0.0 && change > 0.0) {  // Ensure both quality and change are positive
                filtered_paths.push_back(std::make_tuple(path, quality, change));
            }
        }
        return filtered_paths;
    };
#endif

    // Lambda to filter paths by positive rate of change
    auto filter_by_positive_rate = [&y, &current_vertex](
        const std::vector<std::tuple<vertex_shortest_path_info_t, double>>& paths) {

        std::vector<std::tuple<vertex_shortest_path_info_t, double>> filtered_paths;
        for (const auto& [path, rate] : paths) {
            if (rate > 0) {
                filtered_paths.push_back(std::make_tuple(path, rate));
            }
        }
        return filtered_paths;
    };

    // Lambda to filter paths by negative rate of change
    auto filter_by_negative_rate = [&y, &current_vertex](
        const std::vector<std::tuple<vertex_shortest_path_info_t, double>>& paths) {

        std::vector<std::tuple<vertex_shortest_path_info_t, double>> filtered_paths;
        for (const auto& [path, rate] : paths) {
            if (rate < 0) {
                filtered_paths.push_back(std::make_tuple(path, rate));
            }
        }
        return filtered_paths;
    };

    // Main trajectory building loop - continues until both directions are complete
    while (!ascending_done || !descending_done) {
        reachability_map_t reachability_map = compute_graph_reachability_map(current_vertex, scale[current_vertex]);
        // std::vector<vertex_shortest_path_info_t> paths = get_vertex_shortest_paths(reachability_map);  // <<--- this to be removed
        std::vector<vertex_path_t> paths = reconstruct_graph_paths(reachability_map);

        // Skip visited vertices to prevent cycles
        paths = filter_out_visited_vertices(paths);

        // If no paths remain, terminate search
        if (paths.empty()) {
            break;
        }

#if DEBUG__compute_bidirectional_trajectory
        Rprintf("in while(...)\tcurrent_vertex: %zu\n", current_vertex + 1);

        // error("DEBUGGING\n");
        for (const auto& p : paths) {
            Rprintf("%zu:: ", p.vertex);
            print_vect(p.path,"path");
            #if 0
            for (const auto& v : p.path) {
                Rprintf("y[%zu]: %.5f\n", v, y[v]);
            }
            Rprintf("\n");
            #endif
        }
#endif

        // Categorize paths
        std::vector<std::tuple<vertex_path_t, double>> paths_with_no_extrema;
        std::vector<std::tuple<vertex_path_t, double, size_t>> paths_with_minima;
        std::vector<std::tuple<vertex_path_t, double, size_t>> paths_with_maxima;

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
        } // END OF for (size_t i = 1; i < path.vertices.size(); ++i)



#if DEBUG__compute_bidirectional_trajectory
        // old for loop
        for (const auto& path : paths) {
            // Calculate rate of change
            double rate = (y[path.vertex] - y[current_vertex]) / path.distance;

            // Check if path contains extrema
            bool has_extremum = false;
            bool is_maximum = false;
            size_t extremum_index = 0;

            for (size_t i = 1; i < path.vertices.size(); ++i) {
                size_t vertex = path.path[i];
                auto extremum_it = extrema_map.find(vertex);

                if (extremum_it != extrema_map.end()) {
                    has_extremum = true;
                    is_maximum = extremum_it->second;
                    extremum_index = i;
                    break;
                }
            }

            if (!has_extremum) {
                paths_with_no_extrema.push_back(std::make_tuple(path, rate));
            } else if (is_maximum && !ascending_done) {
                paths_with_maxima.push_back(std::make_tuple(path, rate, extremum_index));
            } else if (!descending_done) {  // extremum is minimum
                paths_with_minima.push_back(std::make_tuple(path, rate, extremum_index));
            }
        }
#endif

        // CASE 1: Found both types of extrema - optimal scenario
        if (!paths_with_minima.empty() && !paths_with_maxima.empty()) {

#if DEBUG__compute_bidirectional_trajectory
            Rprintf("Found both types of extrema\n");
#endif

            // Select best minimum (steepest negative rate)
            auto best_min_it = std::min_element(
                paths_with_minima.begin(),
                paths_with_minima.end(),
                [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
            );

            // Select best maximum (steepest positive rate)
            auto best_max_it = std::max_element(
                paths_with_maxima.begin(),
                paths_with_maxima.end(),
                [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
            );

            // Complete both trajectories
            if (!descending_done) {
                const auto& [min_path, min_rate, min_extremum_idx] = *best_min_it;

                // Add subpath from current vertex to the local minimum
                descending_traj.path.insert(
                    descending_traj.path.end(),
                    min_path.path.begin() + 1,  // Skip current vertex
                    min_path.path.begin() + min_extremum_idx + 1  // Include extremum
                );

#if DEBUG__compute_bidirectional_trajectory
                print_vect(descending_traj.path, "descending_traj.path", 1);
#endif

                descending_traj.ends_at_critical = true;
                descending_traj.ends_at_lmax = false;  // It's a minimum
                descending_done = true;

                // Mark these vertices as visited
                for (size_t i = 1; i <= min_extremum_idx; ++i) {
                    visited_vertices.insert(min_path.path[i]);
                }
            }

            if (!ascending_done) {
                const auto& [max_path, max_rate, max_extremum_idx] = *best_max_it;

                // Add subpath from current vertex to the local maximum
                ascending_traj.path.insert(
                    ascending_traj.path.end(),
                    max_path.path.begin() + 1,  // Skip current vertex
                    max_path.path.begin() + max_extremum_idx + 1  // Include extremum
                );

#if DEBUG__compute_bidirectional_trajectory
                print_vect(ascending_traj.path, "ascending_traj.path", 1);
#endif

                ascending_traj.ends_at_critical = true;
                ascending_traj.ends_at_lmax = true;  // It's a maximum
                ascending_done = true;

                // Mark these vertices as visited
                for (size_t i = 1; i <= max_extremum_idx; ++i) {
                    visited_vertices.insert(max_path.path[i]);
                }
            }

            break;  // Both trajectories complete

        // CASE 2: Found only minima
        } else if (!paths_with_minima.empty() && !descending_done) {

#if DEBUG__compute_bidirectional_trajectory
            Rprintf("Found only minima processing descending\n");
            // Rprintf("paths_with_minima.size(): %zu\n", paths_with_minima.size());
#endif

            // Select best minimum (steepest negative rate)
            auto best_min_it = std::min_element(
                paths_with_minima.begin(),
                paths_with_minima.end(),
                [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
            );

            const auto& [min_path, min_rate, min_extremum_idx] = *best_min_it;

            // Complete descending trajectory
            descending_traj.path.insert(
                descending_traj.path.end(),
                min_path.path.begin() + 1,  // Skip current vertex
                min_path.path.begin() + min_extremum_idx + 1  // Include extremum
            );

#if DEBUG__compute_bidirectional_trajectory
            print_vect(descending_traj.path, "descending_traj.path", 1);
#endif

            descending_traj.ends_at_critical = true;
            descending_traj.ends_at_lmax = false;  // It's a minimum
            descending_done = true;

            // Mark these vertices as visited
            for (size_t i = 1; i <= min_extremum_idx; ++i) {
                visited_vertices.insert(min_path.path[i]);
            }

            if (!ascending_done) {
                // Always reset to the last vertex of the ascending trajectory
                current_vertex = ascending_traj.path.back();
                continue; // Critical: restart the loop to recalculate paths for the new position
            }

        // CASE 3: Found only maxima
        } else if (!paths_with_maxima.empty() && !ascending_done) {

#if DEBUG__compute_bidirectional_trajectory
            Rprintf("Found only maxima processing ascending\n");
            //Rprintf("paths_with_maxima.size(): %zu\n", paths_with_maxima.size());
#endif

            // Select best maximum (steepest positive rate)
            auto best_max_it = std::max_element(
                paths_with_maxima.begin(),
                paths_with_maxima.end(),
                [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
            );

            const auto& [max_path, max_rate, max_extremum_idx] = *best_max_it;


            #if DEBUG__compute_bidirectional_trajectory
            print_vect(max_path.path,"max_path.path",1);
            #endif

            // Complete ascending trajectory
            ascending_traj.path.insert(
                ascending_traj.path.end(),
                max_path.path.begin() + 1,  // Skip current vertex
                max_path.path.begin() + max_extremum_idx + 1  // Include extremum
            );

            #if DEBUG__compute_bidirectional_trajectory
            print_vect(ascending_traj.path, "ascending_traj.path", 1);
            #endif

            ascending_traj.ends_at_critical = true;
            ascending_traj.ends_at_lmax = true;  // It's a maximum
            ascending_done = true;

            // Mark these vertices as visited
            for (size_t i = 1; i <= max_extremum_idx; ++i) {
                visited_vertices.insert(max_path.path[i]);
            }

            if (!descending_done) {
                // Always reset to the last vertex of the descending trajectory
                current_vertex = descending_traj.path.back();
                continue; // Critical: restart the loop to recalculate paths for the new position
            }

        // CASE 4: No extrema found - follow gradient
        } else if (!paths_with_no_extrema.empty()) {


            // For ascending direction (if not done)
            if (!ascending_done) {

#if DEBUG__compute_bidirectional_trajectory
                Rprintf("No extrema found: in if(!ascending_done)\n");
#endif

                auto positive_paths = filter_by_positive_rate(paths_with_no_extrema);

                if (!positive_paths.empty()) {
                    // Find path with steepest positive rate
                    auto best_asc_it = std::max_element(
                        positive_paths.begin(),
                        positive_paths.end(),
                        [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
                    );

                    const auto& [best_asc_path, best_asc_rate] = *best_asc_it;

                    // Update current vertex and continue search from there
                    current_vertex = best_asc_path.vertex;
                    visited_vertices.insert(current_vertex);

                    #if DEBUG__compute_bidirectional_trajectory
                    print_vect(best_asc_path.path,"best_asc_path.path",1);
                    #endif

                    // Add path segment to ascending trajectory
                    ascending_traj.path.insert(
                        ascending_traj.path.end(),
                        best_asc_path.path.begin() + 1,  // Skip current vertex
                        best_asc_path.path.end()
                    );

                    #if DEBUG__compute_bidirectional_trajectory
                    print_vect(ascending_traj.path, "ascending_traj.path", 1);
                    #endif

                    // Mark path vertices as visited
                    for (size_t i = 1; i < best_asc_path.path.size(); ++i) {
                        visited_vertices.insert(best_asc_path.path[i]);
                    }

                    continue;  // Continue search from new vertex
                } else {
                    ascending_done = true;  // No positive paths available
                }
            }

            // For descending direction (if not done)
            if (!descending_done) {

#if DEBUG__compute_bidirectional_trajectory
                Rprintf("No extrema found: in if(!descending_done)\n");
#endif


                auto negative_paths = filter_by_negative_rate(paths_with_no_extrema);

                if (!negative_paths.empty()) {
                    // Find path with steepest negative rate
                    auto best_desc_it = std::min_element(
                        negative_paths.begin(),
                        negative_paths.end(),
                        [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); }
                    );

                    const auto& [best_desc_path, best_desc_rate] = *best_desc_it;

                    // Update current vertex and continue search from there
                    current_vertex = best_desc_path.vertex;
                    visited_vertices.insert(current_vertex);

                    // Add path segment to descending trajectory
                    descending_traj.path.insert(
                        descending_traj.path.end(),
                        best_desc_path.path.begin() + 1,  // Skip current vertex
                        best_desc_path.path.end()
                    );

                    #if DEBUG__compute_bidirectional_trajectory
                    print_vect(descending_traj.path, "descending_traj.path", 1);
                    #endif

                    // Mark path vertices as visited
                    for (size_t i = 1; i < best_desc_path.path.size(); ++i) {
                        visited_vertices.insert(best_desc_path.path[i]);
                    }

                    continue;  // Continue search from new vertex
                } else {
                    descending_done = true;  // No negative paths available
                }
            }

        // CASE 5: No valid paths found
        } else {
            break;  // Terminate search
        }
    }

    #if DEBUG__compute_bidirectional_trajectory
    Rprintf("\nFinal results:\n");
    print_vect(ascending_traj.path, "ascending_traj.path", 1);
    print_vect(descending_traj.path, "descending_traj.path", 1);
    #endif

    return {ascending_traj, descending_traj};
}
