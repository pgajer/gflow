
/**
 * @brief Computes candidate bandwidth values from composite paths using a quantile-based strategy
 *
 * Computes min_bw (smallest distance that ensures at least one path has min_path_size original vertices
 * and min_num_grid_vertices grid vertices), max_bw (largest distance in any path), and generates n_bws
 * candidate bandwidth values using quantile-based subdivision of the [min_bw, max_bw] interval.
 *
 * @param composite_paths Vector of composite paths
 * @param min_path_size Minimum required number of original graph vertices
 * @param min_num_grid_vertices Minimum required number of grid vertices
 * @param n_bws Number of bandwidth values to generate
 *
 * @return std::tuple<double, double, std::vector<double>> Returns (min_bw, max_bw, candidate_bws)
 */
std::tuple<double, double, std::vector<double>> uniform_grid_graph_t::get_candidate_bws(
    const std::vector<compose_path_t>& composite_paths,
    size_t min_path_size,
    size_t min_num_grid_vertices,
    size_t n_bws) const {
    // Collect all unique distances from composite paths
    std::set<double> all_distances;
    for (const auto& path : composite_paths) {
        all_distances.insert(path.dist_to_ref_vertex.begin(),
                           path.dist_to_ref_vertex.end());
        // Also include grid distances
        // all_distances.insert(path.grid_dist_to_ref_vertex.begin(),
        //                    path.grid_dist_to_ref_vertex.end());
    }

    // Convert to vector for easier manipulation
    std::vector<double> sorted_distances(all_distances.begin(), all_distances.end());
    std::sort(sorted_distances.begin(), sorted_distances.end());

    // Find max_bw - it's the largest distance in any path
    double max_bw;
    if (sorted_distances.empty()) {

        Rprintf("composite_paths.size(): %zu\n",composite_paths.size());

        REPORT_ERROR("Error: No distances available to calculate bandwidths.");
        return {0.0, 0.0, std::vector<double>{}};
    } else {
        max_bw = sorted_distances.back();
    }

    // Find min_bw - smallest distance that ensures at least one path has min_path_size original vertices
    // and min_num_grid_vertices grid vertices
    double min_bw = max_bw;  // Start with max_bw and work down

    for (const auto& path : composite_paths) {
        // Examine each distance in the path
        for (size_t i = 0; i < path.vertices.size(); ++i) {
            double dist = path.dist_to_ref_vertex[i];

            // Count original vertices within this distance
            size_t orig_vertices_within_dist = std::count_if(
                path.dist_to_ref_vertex.begin(),
                path.dist_to_ref_vertex.end(),
                [dist](double d) { return d <= dist; }
            );

            // Count grid vertices within this distance
            size_t grid_vertices_within_dist = std::count_if(
                path.grid_dist_to_ref_vertex.begin(),
                path.grid_dist_to_ref_vertex.end(),
                [dist](double d) { return d <= dist; }
            );

            // Check if both conditions are satisfied and this distance is smaller than current min_bw
            if (orig_vertices_within_dist >= min_path_size &&
                grid_vertices_within_dist >= min_num_grid_vertices &&
                dist < min_bw) {
                min_bw = dist;
            }
        }

        #if 0
        // Also check distances from grid vertices, since those might have different distributions
        for (size_t i = 0; i < path.grid_vertices.size(); ++i) {
            double dist = path.grid_dist_to_ref_vertex[i];

            // Count original vertices within this distance
            size_t orig_vertices_within_dist = std::count_if(
                path.dist_to_ref_vertex.begin(),
                path.dist_to_ref_vertex.end(),
                [dist](double d) { return d <= dist; }
            );

            // Count grid vertices within this distance
            size_t grid_vertices_within_dist = std::count_if(
                path.grid_dist_to_ref_vertex.begin(),
                path.grid_dist_to_ref_vertex.end(),
                [dist](double d) { return d <= dist; }
            );

            // Check if both conditions are satisfied and this distance is smaller than current min_bw
            if (orig_vertices_within_dist >= min_path_size &&
                grid_vertices_within_dist >= min_num_grid_vertices &&
                dist < min_bw) {
                min_bw = dist;
            }
        }
        #endif
    }

    // Handle the case where min_bw and max_bw are very close
    double bw_range = max_bw - min_bw;
    double min_spacing = 1e-6;  // Minimum meaningful difference between bandwidths

    // If range is too small for requested number of bandwidths
    if (bw_range < (n_bws - 1) * min_spacing) {
        // Option 1: Reduce the number of bandwidths based on the available range
        size_t effective_n_bws = std::max(size_t(2), size_t(bw_range / min_spacing) + 1);

        if (effective_n_bws < n_bws) {
            //if (verbose) {
            REPORT_WARNING("Warning: Bandwidth range (%.6f) is too small for %zu distinct values. "
                           "Reducing to %zu bandwidths.\n",
                           bw_range, n_bws, effective_n_bws);
            //}
            n_bws = effective_n_bws;
        }

        #if 0
        // Option 2: Consider artificial expansion of the range
        // Only do this if we'd end up with very few bandwidths otherwise
        if (n_bws < 5 && max_bw > 0) {
            // Expand the range by adding a percentage to max_bw
            double expansion_factor = 0.1;  // 10% expansion
            double old_max_bw = max_bw;
            max_bw = max_bw * (1.0 + expansion_factor);

            if (verbose) {
                Rprintf("Note: Artificially expanding bandwidth range from [%.6f, %.6f] to [%.6f, %.6f] "
                        "to ensure adequate model diversity.\n",
                        min_bw, old_max_bw, min_bw, max_bw);
            }
        }
        #endif
    }

    // Filter distances to be within [min_bw, max_bw]
    auto it_min = std::lower_bound(sorted_distances.begin(), sorted_distances.end(), min_bw);
    auto it_max = std::upper_bound(sorted_distances.begin(), sorted_distances.end(), max_bw);
    std::vector<double> filtered_distances(it_min, it_max);

    // Generate candidate bandwidths using quantile strategy
    std::vector<double> candidate_bws;
    candidate_bws.reserve(n_bws);

    // Always include min_bw
    candidate_bws.push_back(min_bw);

    // Calculate indices for internal quantiles
    if (n_bws > 2) {
        size_t n_distances = filtered_distances.size();

        if (n_distances <= n_bws) {
            // If we have fewer unique distances than requested bandwidths,
            // just use all available distances
            candidate_bws = filtered_distances;
        } else {
            // Otherwise, use quantile-based selection
            for (size_t i = 1; i < n_bws - 1; ++i) {
                // Calculate index for this quantile
                size_t index = (i * (n_distances - 1)) / (n_bws - 1);
                candidate_bws.push_back(filtered_distances[index]);
            }
        }
    }

    // Add max_bw if we haven't already and if there are multiple values
    if (!candidate_bws.empty() && candidate_bws.back() != max_bw) {
        candidate_bws.push_back(max_bw);
    }

    // Ensure we have no duplicates
    std::sort(candidate_bws.begin(), candidate_bws.end());
    auto last = std::unique(candidate_bws.begin(), candidate_bws.end());
    candidate_bws.erase(last, candidate_bws.end());

    return {min_bw, max_bw, candidate_bws};
}
