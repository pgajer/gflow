#include <unordered_map>
#include <queue>
#include <algorithm>      // For std::find()

#include "harmonic_smoother.hpp"
#include "reachability_map.hpp"
#include "set_wgraph.hpp"
#include "basin.hpp"
#include "error_utils.h"  // For REPORT_ERROR()

/**
 * @brief Calculates a basin complex difference metric between two sets of basins
 *
 * @details This function quantifies how much the basin complex structure has changed
 *          between two iterations by examining changes in extrema and their basins.
 *          The metric incorporates:
 *          1. Appearance/disappearance of extrema
 *          2. Changes in basin sizes
 *          3. Shifts in extrema positions
 *
 * @param basins1 First basin map
 * @param basins2 Second basin map
 * @return Double value representing the basin complex difference (0.0 = identical)
 */
double set_wgraph_t::basin_cx_difference(
    const std::unordered_map<size_t, basin_t>& basins1,
    const std::unordered_map<size_t, basin_t>& basins2
    ) const {
    // If both are empty, they're identical
    if (basins1.empty() && basins2.empty()) {
        return 0.0;
    }

     // If one is empty and the other isn't, they're completely different
    if (basins1.empty() || basins2.empty()) {
        return 1.0;
    }

     // Count extrema that exist in both maps
    size_t common_extrema = 0;
    double size_difference_sum = 0.0;
    double total_basins = basins1.size() + basins2.size();

     // Maps to track which extrema have been matched
    std::unordered_set<size_t> matched_in_basins1;
    std::unordered_set<size_t> matched_in_basins2;

     // First pass: exact vertex matches
    for (const auto& [vertex1, basin1] : basins1) {
        if (basins2.find(vertex1) != basins2.end()) {
            const basin_t& basin2 = basins2.at(vertex1);

             // Found a matching extremum
            common_extrema++;
            matched_in_basins1.insert(vertex1);
            matched_in_basins2.insert(vertex1);

             // Calculate relative size difference
            double size1 = basin1.reachability_map.sorted_vertices.size();
            double size2 = basin2.reachability_map.sorted_vertices.size();
            double max_size = std::max(size1, size2);
            if (max_size > 0) {
                size_difference_sum += std::abs(size1 - size2) / max_size;
            }
        }
    }

     // Second pass: extrema that moved slightly (within a small radius)
    double proximity_threshold = 0.05 * graph_diameter; // 5% of graph diameter

    for (const auto& [vertex1, basin1] : basins1) {
        if (matched_in_basins1.find(vertex1) != matched_in_basins1.end()) {
            continue; // Already matched in first pass
        }

         // Look for nearby extrema in basins2
        for (const auto& [vertex2, basin2] : basins2) {
            if (matched_in_basins2.find(vertex2) != matched_in_basins2.end()) {
                continue; // Already matched
            }

             // Only match extrema of the same type (both maxima or both minima)
            if (basin1.is_maximum != basin2.is_maximum) {
                continue;
            }

             // Check if they're close enough
            double distance = compute_shortest_path_distance(vertex1, vertex2);
            if (distance <= proximity_threshold) {
                 // Found a nearby extremum
                common_extrema++;
                matched_in_basins1.insert(vertex1);
                matched_in_basins2.insert(vertex2);

                 // Calculate relative size difference
                double size1 = basin1.reachability_map.sorted_vertices.size();
                double size2 = basin2.reachability_map.sorted_vertices.size();
                double max_size = std::max(size1, size2);
                if (max_size > 0) {
                    size_difference_sum += std::abs(size1 - size2) / max_size;
                }

                 // Add a penalty for the distance
                size_difference_sum += distance / proximity_threshold;

                break; // Match found, move to next vertex1
            }
        }
    }

     // Calculate Dice-Sørensen similarity
    double dice = (2.0 * common_extrema) / total_basins;

     // Calculate average size difference for matched basins
    double avg_size_diff = (common_extrema > 0) ? size_difference_sum / common_extrema : 0.0;

     // Combine metrics: more weight to Dice, less to size differences
    double difference = 0.7 * (1.0 - dice) + 0.3 * avg_size_diff;

    return difference;
}

/**
 * @brief Identifies local extrema and their basins for a function on the graph
 *
 * @details This streamlined version of basin complex computation focuses on
 *          identifying local extrema and their basins without performing
 *          topological simplification. It's used for tracking how extrema
 *          evolve during the harmonic smoothing process.
 *
 * @param y Vector of function values at each vertex
 * @return Map from extrema vertex indices to their basin structures
 */
std::unordered_map<size_t, basin_t> set_wgraph_t::compute_basins(
    const std::vector<double>& y
) const {
    // 1. Find neighborhood-local extrema
    auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);
    std::unordered_map<size_t, basin_t> basins_map;

    // 2. Find basins for local minima
    bool detect_maxima = false;
    for (const auto& vertex : nbr_lmin) {
        auto basin = find_local_extremum_bfs_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
            basin.extremum_vertex = vertex;
            basins_map[vertex] = std::move(basin);
        }
    }

    // 3. Find basins for local maxima
    detect_maxima = true;
    for (const auto& vertex : nbr_lmax) {
        auto basin = find_local_extremum_bfs_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
            basin.extremum_vertex = vertex;
            basins_map[vertex] = std::move(basin);
        }
    }

    return basins_map;
}

/**
 * @brief Performs harmonic smoothing while tracking prediction landscape evolution
 *
 * @details Enhanced version of harmonic smoothing that monitors how local extrema
 *          and their basins evolve during the smoothing process. This allows for
 *          identifying the optimal amount of smoothing that reduces noise while
 *          preserving significant prediction landscape features.
 *
 * @param harmonic_predictions Vector of function values to be smoothed
 * @param region_vertices Set of vertex indices defining the region to smooth
 * @param max_iterations Maximum number of iterations to perform
 * @param tolerance Convergence threshold for value changes
 * @param record_frequency How often to record states (every N iterations)
 * @param stability_window Number of consecutive iterations to check for stability
 * @param stability_threshold Threshold for considering prediction landscape stable
 * @return harmonic_smoother_t Structure with smoothing history and results
 */
harmonic_smoother_t set_wgraph_t::harmonic_smoother(
    std::vector<double>& harmonic_predictions,
    std::unordered_set<size_t>& region_vertices,
    int max_iterations,
    double tolerance,
    int record_frequency,
    size_t stability_window,
    double stability_threshold
    ) const {
    // Initialize result structure
    harmonic_smoother_t result;

    // Return if region is empty or too small
    if (region_vertices.size() <= 5) {
        REPORT_WARNING("WARNING: region_vertices size %zu is too small for effective harmonic smoothing\n",
                       region_vertices.size());
        return result;
    }

    // 1. Identify boundary vertices (vertices in the region with degree 1 or neighbors outside)
    std::unordered_set<size_t> boundary_vertices;
    for (const size_t& v : region_vertices) {
        // Check if vertex is of degree 1 or has neighbors outside the region
        // Modifying the definition of the boundary to exclude deg 1 vertices as for harmonic smoothing we need to focus on vertices that have neighbors not in the region
        for (const auto& edge : adjacency_list[v]) {
            if (region_vertices.count(edge.vertex) == 0) {
                // This is a boundary vertex
                boundary_vertices.insert(v);
                break;
            }
        }
    }

    // 2. Create interior vertices set (region vertices that are not boundary vertices)
    std::unordered_set<size_t> interior_vertices;
    for (const auto& vertex : region_vertices) {
        if (boundary_vertices.find(vertex) == boundary_vertices.end()) {
            interior_vertices.insert(vertex);
        }
    }

    // Return if we have no interior vertices to smooth
    if (interior_vertices.empty()) {
        REPORT_WARNING("WARNING: No interior vertices to smooth\n");
        return result;
    }

    // 3. Record initial state
    result.i_harmonic_predictions.push_back(harmonic_predictions);
    result.i_basins.push_back(compute_basins(harmonic_predictions));

    // Temporary vector for new values during iteration
    std::vector<double> new_values = harmonic_predictions;

    // 4. Iterative relaxation for harmonic interpolation with prediction shape tracking
    bool converged = false;
    bool basin_cx_stable = false;

    for (int iter = 0; iter < max_iterations && !converged; ++iter) {
        // Reset convergence flag for this iteration
        converged = true;

        // Update each interior vertex value
        for (const size_t& v : interior_vertices) {
            // Average the values of neighbors
            double sum = 0.0;
            double weight_sum = 0.0;

            for (const auto& edge : adjacency_list[v]) {
                // Use inverse of edge weight as weighting factor
                double weight = 1.0 / (edge.weight + 1e-10); // Avoid division by zero
                sum += harmonic_predictions[edge.vertex] * weight;
                weight_sum += weight;
            }

            if (weight_sum > 0) {
                double weighted_avg = sum / weight_sum;
                if (std::abs(new_values[v] - weighted_avg) > tolerance) {
                    converged = false;
                }
                new_values[v] = weighted_avg;
            }
        }

        // Update the values for next iteration
        for (const size_t& v : interior_vertices) {
            harmonic_predictions[v] = new_values[v];
        }

        // Record state if it's time (based on record_frequency)
        if ((iter + 1) % record_frequency == 0 || converged) {
            // Save current harmonic predictions
            result.i_harmonic_predictions.push_back(harmonic_predictions);

            // Compute and save current basins
            auto current_basins = compute_basins(harmonic_predictions);
            result.i_basins.push_back(current_basins);

            // Calculate basin complex difference if we have at least two recorded states
            if (result.i_basins.size() >= 2) {
                size_t last_idx = result.i_basins.size() - 1;
                double diff = basin_cx_difference(
                    result.i_basins[last_idx - 1],
                    result.i_basins[last_idx]
                    );
                result.basin_cx_differences.push_back(diff);

                // Check if basin complex structure has stabilized
                if (!basin_cx_stable &&
                    result.is_basin_cx_stable(stability_window, stability_threshold)) {
                    basin_cx_stable = true;
                    result.stable_iteration = iter + 1;

                    // If we're only interested in finding the stable point, we could break here
                    // but we'll continue to track the full evolution
                }
            }
        }
    }

    // If we never detected stability, set stable_iteration to max_iterations
    if (result.stable_iteration == 0) {
        result.stable_iteration = max_iterations;
    }

    return result;
}



/**
 * @brief Performs harmonic repair of function values when one basin absorbs another
 *
 * @details This function implements the critical harmonic interpolation step during
 *          topological simplification of a basin complex. When a less significant
 *          basin is absorbed into a more significant one during cancellation of
 *          topological features, this function ensures that the function values
 *          in the absorbed region are smoothly repaired.
 *
 * The algorithm proceeds through these key steps:
 * 1. Identify all vertices in the absorbed basin
 * 2. Determine boundary vertices (vertices with neighbors outside the basin or
 *    vertices shared with the absorbing basin)
 * 3. Fix the values at boundary vertices as constraints
 * 4. Iteratively update interior vertex values using a weighted average of neighbor values
 * 5. Make special adjustments for the absorbed extremum to ensure smooth blending
 *
 * The harmonic repair preserves the main topological features of the absorbing
 * basin while eliminating the local extremum of the absorbed basin. This creates
 * a smoothly varying function with fewer spurious extrema.
 *
 * @param[in,out] harmonic_predictions  Vector of function values to be repaired in place
 * @param[in] absorbing_basin  The basin that will absorb another basin
 * @param[in] absorbed_basin   The basin being absorbed
 *
 * @pre Both basins must have valid vertices with indices less than the graph size
 * @pre The size of harmonic_predictions must match the number of vertices in the graph
 *
 * @post The values in harmonic_predictions are harmonically repaired within the absorbed
 *       basin, preserving boundary values and maintaining a smooth transition
 *
 * @note This function is a key component of the topological simplification process
 *       in gradient flow basin complex construction
 * @note The harmonic repair respects the geometry of the graph through edge-weighted averaging
 * @note Special care is taken for the absorbed extremum point to ensure proper blending
 *
 * @see basin_t
 * @see create_basin_cx
 * @see absorb_basin
 */
void set_wgraph_t::perform_harmonic_repair(
    std::vector<double>& harmonic_predictions,
    const basin_t& absorbing_basin,
    const basin_t& absorbed_basin
) const {
     // Collect all vertices in the absorbed basin
    std::unordered_set<size_t> basin_vertices;
    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        basin_vertices.insert(v_info.vertex);
    }

     // Identify boundary vertices (vertices in the basin with neighbors outside)
    std::unordered_map<size_t, double> boundary_values;

     // Add the absorbing basin's extremum as a boundary vertex with fixed value
    size_t absorbing_extremum = absorbing_basin.extremum_vertex;
    if (absorbing_extremum < harmonic_predictions.size()) {
        boundary_values[absorbing_extremum] = harmonic_predictions[absorbing_extremum];
    }

     // Add other boundary vertices from the absorbed basin
    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        size_t v = v_info.vertex;

         // Skip the extremum itself - it will be recalculated
        if (v == absorbed_basin.extremum_vertex) {
            continue;
        }

         // Check if vertex index is valid
        if (v >= adjacency_list.size()) {
            continue;
        }

         // Check if this vertex has neighbors outside the basin
        for (const auto& edge : adjacency_list[v]) {
            if (basin_vertices.count(edge.vertex) == 0) {
                 // This is a boundary vertex - keep its original value
                if (v < harmonic_predictions.size()) {
                    boundary_values[v] = harmonic_predictions[v];
                }
                break;
            }
        }
    }

     // Include common vertices between absorbing and absorbed basins as boundary values
    for (const auto& v_info : absorbing_basin.reachability_map.sorted_vertices) {
        size_t v = v_info.vertex;
        if (basin_vertices.count(v) > 0 && v != absorbed_basin.extremum_vertex) {
            if (v < harmonic_predictions.size()) {
                boundary_values[v] = harmonic_predictions[v];
            }
        }
    }

     // Simple Jacobi iteration for harmonic interpolation
    const int MAX_ITERATIONS = 100;
    const double TOLERANCE = 1e-6;

    std::vector<double> new_values = harmonic_predictions;

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        bool converged = true;

        for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
            size_t v = v_info.vertex;

             // Skip boundary vertices - their values are fixed
            if (boundary_values.count(v) > 0) {
                continue;
            }

             // Skip invalid vertex indices
            if (v >= adjacency_list.size() || v >= harmonic_predictions.size()) {
                continue;
            }

             // Average the values of neighbors
            double sum = 0.0;
            double weight_sum = 0.0;

            for (const auto& edge : adjacency_list[v]) {
                 // Check that neighbor index is valid
                if (edge.vertex >= harmonic_predictions.size()) {
                    continue;
                }

                 // Use inverse of edge weight as weighting factor
                double weight = 1.0 / (edge.weight + 1e-10); // Avoid division by zero
                sum += harmonic_predictions[edge.vertex] * weight;
                weight_sum += weight;
            }

            if (weight_sum > 0) {
                double weighted_avg = sum / weight_sum;
                if (std::abs(new_values[v] - weighted_avg) > TOLERANCE) {
                    converged = false;
                }
                new_values[v] = weighted_avg;
            }
        }

         // Update the values for next iteration
        for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
            size_t v = v_info.vertex;
            if (boundary_values.count(v) == 0 && v < harmonic_predictions.size()) {
                harmonic_predictions[v] = new_values[v];
            }
        }

        if (converged) {
            break;
        }
    }

     // Make a special adjustment for the absorbed extremum
    size_t absorbed_extremum = absorbed_basin.extremum_vertex;

     // Check validity of absorbed extremum
    if (absorbed_extremum >= adjacency_list.size() || absorbed_extremum >= harmonic_predictions.size()) {
        return;
    }

    double neighbor_avg = 0.0;
    double weight_sum = 0.0;

    for (const auto& edge : adjacency_list[absorbed_extremum]) {
        if (edge.vertex >= harmonic_predictions.size()) {
            continue;
        }

        double weight = 1.0 / (edge.weight + 1e-10);
        neighbor_avg += harmonic_predictions[edge.vertex] * weight;
        weight_sum += weight;
    }

    if (weight_sum > 0) {
        neighbor_avg /= weight_sum;

         // Calculate distance between extrema
        double extrema_distance = 1.0;
        if (absorbing_basin.reachability_map.distances.count(absorbed_extremum) > 0) {
            extrema_distance = absorbing_basin.reachability_map.distances.at(absorbed_extremum);
        } else {
             // Compute approximate distance if not directly available
            double d = compute_shortest_path_distance(
                absorbing_basin.extremum_vertex,
                absorbed_extremum
            );
            extrema_distance = std::max(1.0, d);
        }

         // Blend based on distance
        double blend_factor = 1.0 / (1.0 + extrema_distance);
        harmonic_predictions[absorbed_extremum] =
            blend_factor * harmonic_predictions[absorbing_basin.extremum_vertex] +
            (1.0 - blend_factor) * neighbor_avg;
    }
}


/**
 * @brief Harmonic (Dirichlet) smoothing within a region, keeping boundary values fixed.
 *
 * @details
 * This function performs harmonic smoothing (discrete Laplace equation) on the induced
 * subgraph G[R] of a weighted graph, where R is the region vertex set.
 *
 * Boundary vertices B are defined as vertices in R that have at least one neighbor
 * outside R. Boundary values are held fixed, and interior vertices I = R \ B are
 * iteratively updated so that for each v in I:
 *
 *   f(v) = sum_{u in N(v) ∩ R} c(v,u) f(u) / sum_{u in N(v) ∩ R} c(v,u),
 *
 * where c(v,u) are conductances. If edge weights represent distances, conductances
 * are taken as c = 1/(d + eps). If edge weights are already conductances/affinities,
 * use c = w directly.
 *
 * The function also computes and reports the residual measure:
 *
 *   max_{v in I} | f(v) - avg_w( f(N(v) ∩ R) ) |
 *
 * which quantifies how close the final solution is to being harmonic on the interior.
 *
 * @param[in,out] harmonic_predictions Full graph vector of values, modified in place on interior vertices.
 * @param[in] region_vertices Vertex set R to smooth within.
 * @param[in] max_iterations Maximum number of iterations.
 * @param[in] tolerance Convergence threshold applied to both max_change and max_residual.
 * @param[in] edge_weight_is_distance If true, treat edge.weight as distance (use 1/(d+eps));
 *                                    otherwise treat edge.weight as conductance (use w).
 * @param[in] verbose If true, prints progress and final diagnostics.
 *
 * @return harmonic_smoothing_stats_t Diagnostics for the run.
 */
harmonic_smoothing_stats_t set_wgraph_t::perform_harmonic_smoothing(
    std::vector<double>& harmonic_predictions,
    const std::unordered_set<size_t>& region_vertices,
    int max_iterations,
    double tolerance,
    bool edge_weight_is_distance,
    bool verbose
) const {

    harmonic_smoothing_stats_t stats;
    stats.num_region = region_vertices.size();

    ensure_edge_weights_computed();

    // ---- Basic validation
    if (region_vertices.empty()) {
        REPORT_WARNING("perform_harmonic_smoothing: region_vertices is empty");
        return stats;
    }
    if (harmonic_predictions.size() != adjacency_list.size()) {
        REPORT_ERROR("perform_harmonic_smoothing: harmonic_predictions size (%zu) != graph size (%zu)",
                     harmonic_predictions.size(), adjacency_list.size());
    }
    if (max_iterations <= 0) {
        REPORT_WARNING("perform_harmonic_smoothing: max_iterations <= 0 (got %d)", max_iterations);
        return stats;
    }
    if (!(tolerance > 0.0)) {
        REPORT_WARNING("perform_harmonic_smoothing: tolerance must be > 0 (got %g)", tolerance);
        return stats;
    }

    for (const size_t v : region_vertices) {
        if (v >= adjacency_list.size()) {
            REPORT_ERROR("perform_harmonic_smoothing: region vertex %zu out of range (graph size %zu)",
                         v, adjacency_list.size());
        }
    }

    const double eps = 1e-10;

    auto conductance_of = [&](double edge_weight) -> double {
        if (edge_weight_is_distance) {
            return 1.0 / (edge_weight + eps);
        }
        return (edge_weight > 0.0) ? edge_weight : 0.0;
    };

    // ---- Identify boundary B: v in R with at least one neighbor outside R
    std::unordered_set<size_t> boundary_vertices;
    boundary_vertices.reserve(region_vertices.size());

    for (const size_t v : region_vertices) {
        for (const auto& edge : adjacency_list[v]) {
            if (region_vertices.count(edge.vertex) == 0) {
                boundary_vertices.insert(v);
                break;
            }
        }
    }

    stats.num_boundary = boundary_vertices.size();

    if (boundary_vertices.empty()) {
        REPORT_WARNING("perform_harmonic_smoothing: boundary is empty; Dirichlet problem ill-posed on a closed region. No changes made.");
        return stats;
    }

    // ---- Interior I = R \ B
    std::unordered_set<size_t> interior_vertices;
    interior_vertices.reserve(region_vertices.size() - boundary_vertices.size());

    for (const size_t v : region_vertices) {
        if (boundary_vertices.count(v) == 0) {
            interior_vertices.insert(v);
        }
    }

    stats.num_interior = interior_vertices.size();

    if (interior_vertices.empty()) {
        if (verbose) {
            Rprintf("perform_harmonic_smoothing: no interior vertices; nothing to smooth.\n");
        }
        return stats;
    }

    auto weighted_avg_in_region = [&](size_t v, const std::vector<double>& values, double& out_wsum) -> double {
        double sum = 0.0;
        double wsum = 0.0;

        for (const auto& edge : adjacency_list[v]) {
            const size_t u = edge.vertex;
            if (region_vertices.count(u) == 0) {
                continue; // restrict to N(v) ∩ R
            }
            const double c = conductance_of(edge.weight);
            if (c <= 0.0) {
                continue;
            }
            sum += c * values[u];
            wsum += c;
        }

        out_wsum = wsum;
        return (wsum > 0.0) ? (sum / wsum) : values[v];
    };

    // Boundary hull for optional warning
    double boundary_min = std::numeric_limits<double>::infinity();
    double boundary_max = -std::numeric_limits<double>::infinity();
    for (const size_t b : boundary_vertices) {
        boundary_min = std::min(boundary_min, harmonic_predictions[b]);
        boundary_max = std::max(boundary_max, harmonic_predictions[b]);
    }

    // Reserve per-iteration diagnostics
    stats.max_change.reserve(static_cast<size_t>(max_iterations));
    stats.max_residual.reserve(static_cast<size_t>(max_iterations));

    // ---- Iterative relaxation (Jacobi)
    std::vector<double> new_values = harmonic_predictions;

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_change_iter = 0.0;

        // Propose updates
        for (const size_t v : interior_vertices) {
            double wsum = 0.0;
            const double avg = weighted_avg_in_region(v, harmonic_predictions, wsum);

            const double proposed = (wsum > 0.0) ? avg : harmonic_predictions[v];

            const double change = std::abs(proposed - harmonic_predictions[v]);
            if (change > max_change_iter) {
                max_change_iter = change;
            }
            new_values[v] = proposed;
        }

        // Commit updates
        for (const size_t v : interior_vertices) {
            harmonic_predictions[v] = new_values[v];
        }

        // Residual after update:
        // max_{v in I} | f(v) - avg_w( f(N(v) ∩ R) ) |
        double max_residual_iter = 0.0;
        for (const size_t v : interior_vertices) {
            double wsum = 0.0;
            const double avg = weighted_avg_in_region(v, harmonic_predictions, wsum);
            if (wsum <= 0.0) {
                continue;
            }
            const double residual = std::abs(harmonic_predictions[v] - avg);
            if (residual > max_residual_iter) {
                max_residual_iter = residual;
            }
        }

        // Record diagnostics
        stats.max_change.push_back(max_change_iter);
        stats.max_residual.push_back(max_residual_iter);
        stats.num_iterations = iter + 1;

        // Convergence: require both to be small
        if (max_change_iter < tolerance && max_residual_iter < tolerance) {
            stats.converged = true;
            break;
        }
    }

    // Optional convex-hull warning
    {
        const double hull_eps = 10.0 * tolerance;
        bool violated = false;
        for (const size_t v : interior_vertices) {
            const double x = harmonic_predictions[v];
            if (x < boundary_min - hull_eps || x > boundary_max + hull_eps) {
                violated = true;
                break;
            }
        }
        if (violated) {
            REPORT_WARNING("perform_harmonic_smoothing: interior values exceed boundary convex hull (possible negative/zero weights, leakage, or ill-conditioned conductances).");
        }
    }

    if (verbose) {
        const double last_change = stats.max_change.empty() ? 0.0 : stats.max_change.back();
        const double last_resid  = stats.max_residual.empty() ? 0.0 : stats.max_residual.back();
        Rprintf("perform_harmonic_smoothing: iters=%d converged=%s last_max_change=%g last_max_residual=%g |R|=%zu |B|=%zu |I|=%zu\n",
                stats.num_iterations,
                stats.converged ? "true" : "false",
                last_change,
                last_resid,
                stats.num_region,
                stats.num_boundary,
                stats.num_interior);
    }

    return stats;
}
