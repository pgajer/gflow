#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <mutex>                    // For std::mutex
#include <execution>                // For std::execution::par_unseq
#include <atomic>                   // For std::atomic
#include <thread>                   // For std::thread::hardware_concurrenyc
#include <queue>                    // Foe std::queue

#include "cpp_utils.hpp"            // For debugging

#include "graph_deg0_lowess_buffer_cv.hpp" // For graph_deg0_lowess_buffer_cv_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

/**
 * @brief Perform degree-0 LOWESS with buffer zone cross-validation for bandwidth selection
 *
 * @details This function implements a graph-based extension of LOWESS degree-0
 * (locally weighted average) with spatially-stratified cross-validation using
 * buffer zones to prevent spatial autocorrelation from biasing the performance estimates.
 * The algorithm:
 * 1. Creates a maximal packing of vertices to serve as fold seed points
 * 2. Assigns all vertices to the nearest seed point to form spatially coherent folds
 * 3. For each fold, creates a buffer zone around test vertices to ensure spatial independence
 * 4. For each candidate bandwidth, performs cross-validation across the folds
 * 5. Selects the bandwidth with the lowest cross-validation error
 * 6. Fits the final model with the optimal bandwidth
 *
 * @param y Response values at each vertex in the graph
 * @param min_bw_factor Minimum bandwidth as a factor of graph diameter
 * @param max_bw_factor Maximum bandwidth as a factor of graph diameter
 * @param n_bws Number of bandwidths to test
 * @param log_grid Whether to use logarithmic spacing for bandwidth grid
 * @param kernel_type Type of kernel function for weighting
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 * @param buffer_hops Number of hops for buffer zone around test vertices
 * @param auto_buffer_hops Whether to automatically determine optimal buffer size
 * @param n_folds Number of cross-validation folds
 * @param with_bw_predictions Whether to compute and store predictions for all bandwidths
 * @param precision Precision for bandwidth grid computation
 * @param verbose Whether to print progress information
 *
 * @return graph_deg0_lowess_cv_t Structure containing optimal model results and CV diagnostics
 *
 * @note If auto_buffer_hops is true, the buffer_hops parameter is ignored and
 * an optimal buffer size is determined based on spatial autocorrelation analysis.
 */
graph_deg0_lowess_buffer_cv_t set_wgraph_t::graph_deg0_lowess_buffer_cv(
    const std::vector<double>& y,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool log_grid,
    size_t kernel_type,
    double dist_normalization_factor,
    bool use_uniform_weights,
    size_t buffer_hops,
    bool auto_buffer_hops,
    size_t n_folds,
    bool with_bw_predictions,
    double precision,
    bool verbose) {

    size_t domain_min_size = 3;

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    // Initialize result structure
    graph_deg0_lowess_buffer_cv_t result;

    // Auto-determine buffer hops if requested
    if (auto_buffer_hops) {
        buffer_hops = determine_optimal_buffer_hops(y, verbose);
    }

    // Always set the buffer_hops_used field
    result.buffer_hops_used = buffer_hops;

    size_t n_vertices = adjacency_list.size();

    std::vector<std::vector<size_t>> folds;
    if (n_folds >= n_vertices) {  // LOO CV case

        n_folds = n_vertices;
        folds.resize(n_folds);
        for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
            folds[vertex] = {vertex};
        }
        // Estimate the diameter of the graph
        auto [end1, diam] = get_vertex_eccentricity(0);   // Start from vertex 0
        auto [end2, diameter] = get_vertex_eccentricity(end1);
        graph_diameter = diameter;
    } else {
        // Create spatially stratified folds
        folds = create_spatially_stratified_folds(n_folds);
    }

    if (graph_diameter <= 0) {
        REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
    }

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    // Generate bandwidth grid
    std::vector<double> bw_grid = get_candidate_bws(
        min_bw,
        max_bw,
        n_bws,
        log_grid,
        precision);

    result.bws = bw_grid;
    result.bw_errors.resize(bw_grid.size(), std::numeric_limits<double>::quiet_NaN());

    if (with_bw_predictions) {
        result.bw_predictions.resize(bw_grid.size(),
                                     std::vector<double>(y.size(), 0.0));
    }

    if (verbose) {
        Rprintf("Starting graph_deg0_lowess_buffer_cv() with Buffer Zone Method\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Number of CV folds: %zu\n", n_folds);
        Rprintf("Buffer zone hops: %zu\n", buffer_hops);

        Rprintf("min_bw_factor: %f\n", min_bw_factor);
        Rprintf("min_bw: %f\n", min_bw);

        Rprintf("max_bw_factor: %f\n", max_bw_factor);
        Rprintf("max_bw: %f\n", max_bw);

        Rprintf("graph_diameter: %f\n", graph_diameter);

        Rprintf("n_bws: %zu\n", n_bws);
        print_vect(bw_grid, "bw_grid");

        Rprintf("folds:\n");
        for (size_t fold = 0; fold < n_folds; ++fold) {
            Rprintf("%zu: ", fold);
            print_vect(folds[fold]);
        }
    }

    // For each bandwidth, perform cross-validation
    for (size_t bw_idx = 0; bw_idx < bw_grid.size(); ++bw_idx) {

        double bandwidth = bw_grid[bw_idx];
        double total_error = 0.0;
        size_t total_test_points = 0;

        size_t debugging_counter = 0;

        // Process each fold
        for (size_t fold = 0; fold < folds.size(); ++fold) {

            const std::vector<size_t>& test_vertices = folds[fold];

            // 1. Create buffer zone arond test vertices
            std::unordered_set<size_t> buffer_zone = create_buffer_zone(
                test_vertices
                buffer_hops
                );

            // 2. Training vertices are the vertices of the graph that do not belong to buffer zone
            std::vector<size_t> training_vertices;
            training_vertices.reserve(n_vertices - buffer_zone.size());   // Pre-allocate for efficiency

            for (size_t i = 0; i < n_vertices; ++i) {
                if (buffer_zone.find(i) == buffer_zone.end()) {
                    training_vertices.push_back(i);
                }
            }

            // 3. Generate predictions fot the training set vertices only given kernel mean given the current bandwidth
            std::vector<double> trainig_predictions = get_trainig_predictions(y, training_vertices, bandwidth);

            // Process each test vertex
            for (size_t vertex : test_vertices) {
                // Predict value for this test vertex using training vertices outside buffer zone
                auto prediction = predict_test_vertex_with_buffer(
                    vertex,
                    buffer_hops,
                    y,
                    folds[fold],
                    dist_normalization_factor,
                    use_uniform_weights);

                // if (debugging_counter >= 5) error("DEBUGGING");

                Rprintf("bw_idx: %zu fold: %zu v: %zu y[v]: %.3f y.pred[v]: %.3f", bw_idx, fold, vertex, y[vertex], prediction);

                // If prediction was possible
                if (!std::isnan(prediction)) {
                    // Calculate squared error
                    double error = std::pow(prediction - y[vertex], 2);
                    total_error += error;
                    total_test_points++;

                    Rprintf(" error: %.3f total_error: %.3f\n", error, total_error);

                    // Store prediction if requested
                    if (with_bw_predictions) {
                        result.bw_predictions[bw_idx][vertex] = prediction;
                    }
                } else {
                    REPORT_ERROR("ERROR processing bw %zu fold %zu vertex %zu. NaN prediction value for the vertex.",
                                 bw_idx, fold, vertex);
                }
            }
        }
        
        // Calculate average error for this bandwidth
        if (total_test_points > 0) {
            result.bw_errors[bw_idx] = total_error / total_test_points;
        } else {
            result.bw_errors[bw_idx] = std::numeric_limits<double>::infinity();
            REPORT_ERROR("ERROR processing bw %zu. No predictions found!", bw_idx);
        }
    }

    // Find optimal bandwidth
    auto min_it = std::min_element(result.bw_errors.begin(), result.bw_errors.end());
    result.opt_bw_idx = min_it - result.bw_errors.begin();
    result.opt_bw = bw_grid[result.opt_bw_idx];

    // Final prediction with optimal bandwidth
    result.predictions = predict_all_with_optimal_bandwidth(
        result.opt_bw,
        y,
        dist_normalization_factor,
        use_uniform_weights,
        min_bw,
        max_bw,
        domain_min_size
        );

    return result;
}

/**
 * @brief Predicts the value at a test vertex using training vertices outside a buffer zone
 *
 * @details
 * This function implements a spatial cross-validation approach for graph-based LOWESS
 * that preserves independence between test and training sets by using buffer zones.
 * For a given test vertex, it:
 * 1. Creates a buffer zone of specified hop distance around the test vertex
 * 2. Finds training vertices outside this buffer zone (starting at buffer_hops+1)
 * 3. Applies kernel weighting based on geodesic distances to compute the prediction
 *
 * The approach ensures that predictions for test vertices are made using only
 * training vertices that are spatially separated, which produces more realistic
 * cross-validation error estimates by preventing information leakage through
 * the graph structure.
 *
 * @param test_vertex Index of the vertex for which to make a prediction
 * @param bandwidth Current bandwidth value being evaluated in cross-validation
 * @param buffer_hops Number of hops to use for the buffer zone
 * @param y Vector of response values at each vertex in the graph
 * @param current_fold Vector of vertex indices that belong to the current test fold
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 *
 * @return Predicted value for the test vertex
 *
 * @throws std::runtime_error If no valid training vertices can be found outside the buffer zone
 * @throws std::runtime_error If the sum of weights is zero (should not occur with proper normalization)
 *
 * @note Training vertices are selected based on hop distance from the test vertex,
 *       not the bandwidth parameter. The bandwidth is applied in the weighting stage.
 * @note If use_uniform_weights is true, all training vertices have equal weight;
 *       otherwise, weights decrease with distance according to the kernel function.
 * @note The function shifts distances so the closest training vertex has distance 0,
 *       ensuring it receives the maximum kernel weight.
 */
double set_wgraph_t::predict_test_vertex_with_buffer(
    size_t test_vertex,
    size_t buffer_hops,
    const std::vector<double>& y,
    const std::vector<size_t>& current_fold,
    double dist_normalization_factor,
    bool use_uniform_weights) {

    // Lambda function to find vertices at specific hop distance
    auto find_vertices_at_hop = [this](size_t start_vertex, size_t target_hop) -> std::unordered_map<size_t, size_t> {
        std::unordered_map<size_t, size_t> hop_map;  // maps vertex -> hop distance
        std::queue<std::pair<size_t, size_t>> queue;  // (vertex, hop_level)
        std::unordered_set<size_t> visited;

        // Start BFS from the test vertex
        queue.push({start_vertex, 0});
        visited.insert(start_vertex);

        while (!queue.empty()) {
            auto [current, hop] = queue.front();
            queue.pop();

            // If we've reached our target hop, add to result
            if (hop == target_hop) {
                hop_map[current] = hop;
                continue;  // Don't explore beyond this hop
            }

            // If we're still within the buffer zone, keep exploring
            if (hop < target_hop) {
                for (const auto& edge : adjacency_list[current]) {
                    size_t neighbor = edge.vertex;
                    if (visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        queue.push({neighbor, hop + 1});
                    }
                }
            }
        }

        return hop_map;
    };

    // Create buffer zone around test vertex
    std::unordered_set<size_t> buffer_zone = create_buffer_zone({test_vertex}, buffer_hops);

    // Create set of current fold vertices for efficient lookup
    std::unordered_set<size_t> current_fold_set(current_fold.begin(), current_fold.end());

    // Start with buffer hop distance and increase until we find training vertices
    std::vector<size_t> training_vertices;
    std::vector<double> distances;
    size_t current_hop = buffer_hops + 1;
    bool found_training_vertices = false;

    while (!found_training_vertices && current_hop <= buffer_hops + 5) {  // Limit search depth
        // Find vertices at exactly this hop distance
        auto hop_vertices = find_vertices_at_hop(test_vertex, current_hop);

        // Filter to include only valid training vertices
        for (const auto& [vertex, hop] : hop_vertices) {
            // Skip if vertex is in current fold
            if (current_fold_set.find(vertex) != current_fold_set.end()) {
                continue;
            }

            // Calculate actual graph distance for weighting
            double distance = compute_shortest_path_distance(test_vertex, vertex);

            // Include all valid training vertices outside buffer zone
            training_vertices.push_back(vertex);
            distances.push_back(distance);
        }

        // Check if we found any valid training vertices
        if (!training_vertices.empty()) {
            found_training_vertices = true;
        } else {
            // Increase hop distance and try again
            current_hop++;
        }
    }

    if (training_vertices.empty()) {
        // we should never find ourselves here !!!
        REPORT_ERROR("No valid training vertices were found for test vertex: %zu\n", test_vertex);
    }

    // Calculate prediction based on training vertices
    std::vector<double> weights(training_vertices.size());

    if (use_uniform_weights) {
        // Use uniform weights
        std::fill(weights.begin(), weights.end(), 1.0);
    } else {
        // First, find the minimum distance for shifting
        double min_dist = *std::min_element(distances.begin(), distances.end());

        // Prepare normalized distances vector
        std::vector<double> normalized_dists(training_vertices.size());

        // The closest training vertices should receive the highest weight
        // Shift distances so the smallest is 0, then normalize
        for (size_t i = 0; i < distances.size(); ++i) {
            normalized_dists[i] = distances[i] - min_dist;   // Shift by minimum
        }

        // Find the new maximum after shifting
        double max_dist = *std::max_element(normalized_dists.begin(), normalized_dists.end());
        if (max_dist == 0) max_dist = 1;
        max_dist *= dist_normalization_factor;

        // Normalize by the scaled maximum
        for (size_t i = 0; i < normalized_dists.size(); ++i) {
            normalized_dists[i] /= max_dist;
        }

        // Use the kernel function to calculate weights
        kernel_fn(normalized_dists.data(), static_cast<int>(training_vertices.size()), weights.data());
    }

    // Compute weighted average
    double weighted_sum = 0.0;
    double weight_sum = 0.0;

    for (size_t i = 0; i < training_vertices.size(); ++i) {
        weighted_sum += weights[i] * y[training_vertices[i]];
        weight_sum += weights[i];
    }

    if (weight_sum <= 0) {
        // Based on how we defined normalized distances this should never happen
        REPORT_ERROR("Undefined prediction value for test vertex: %zu\n", test_vertex);
    }

    return weighted_sum / weight_sum;
}

/**
 * @brief Creates spatially stratified cross-validation folds for graph-based data
 *
 * @details
 * This function creates cross-validation folds that preserve spatial structure by:
 * 1. Generating a maximal packing of seed vertices distributed across the graph
 * 2. Assigning each vertex to the fold corresponding to its nearest seed
 *
 * This approach ensures that each fold contains vertices from different regions
 * of the graph, creating spatially coherent partitions. Unlike random partitioning,
 * spatially stratified folds maintain the underlying spatial structure of the data,
 * which is essential for evaluating methods on graph-structured data.
 *
 * The algorithm uses geodesic distances (shortest paths) on the graph to:
 * - Select seed points that are well-separated
 * - Assign vertices to folds based on proximity
 *
 * @param n_folds Desired number of cross-validation folds
 *
 * @return std::vector<std::vector<size_t>> A vector containing n_folds vectors,
 *         where each inner vector contains the vertex indices for that fold
 *
 * @note The actual number of folds may be less than requested if the graph structure
 *       cannot support the desired number of well-separated regions
 * @note The function uses maximal packing to ensure seed points are well-distributed
 *       across the graph, which helps create balanced folds
 * @note For small graphs, the function automatically adjusts the target number of
 *       seed points to ensure reasonable fold creation
 */
std::vector<std::vector<size_t>> set_wgraph_t::create_spatially_stratified_folds(
    size_t n_folds
    ) {
    // Use maximal packing to create spatially coherent folds
    size_t n_vertices = adjacency_list.size();

    // Create seed points using maximal packing
    size_t target_seeds = std::max(n_folds * 3, n_vertices / 20);
    size_t max_iterations = 100;
    double precision = 1e-6;
    std::vector<size_t> seed_vertices = create_maximal_packing(
        target_seeds, max_iterations, precision);

    size_t actual_seeds = seed_vertices.size();
    size_t actual_folds = std::min(n_folds, actual_seeds);

    // Group seeds into folds
    std::vector<std::vector<size_t>> seed_groups;
    if (actual_seeds <= actual_folds) {
        seed_groups.resize(actual_seeds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i].push_back(seed_vertices[i]);
        }
    } else {
        seed_groups.resize(actual_folds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i % actual_folds].push_back(seed_vertices[i]);
        }
    }

    // Assign vertices to nearest seed group
    std::vector<std::vector<size_t>> folds(actual_folds);

    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
        double min_distance = std::numeric_limits<double>::infinity();
        size_t closest_fold = 0;

        for (size_t fold_idx = 0; fold_idx < seed_groups.size(); ++fold_idx) {
            for (size_t seed : seed_groups[fold_idx]) {
                double distance = compute_shortest_path_distance(vertex, seed);
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_fold = fold_idx;
                }
            }
        }

        folds[closest_fold].push_back(vertex);
    }

    return folds;
}

/**
 * @brief Generate predictions for all vertices using the optimal bandwidth with adaptive domain size adjustment
 *
 * @details
 * This function produces the final predictions for all vertices in the graph
 * using the optimal bandwidth determined through buffer zone cross-validation.
 * It implements an adaptive approach that ensures robust predictions by:
 *
 * 1. Initially attempting to use the optimal bandwidth for each vertex
 * 2. Dynamically expanding the bandwidth if insufficient neighbors are found
 * 3. Applying consistent distance normalization and kernel weighting
 *
 * For each vertex, the algorithm:
 * - Searches for neighbors within the optimal bandwidth radius
 * - If too few vertices are found, increases the bandwidth to ensure at least domain_min_size neighbors
 * - Shifts distances so the closest neighbor has distance 0
 * - Applies kernel weighting based on normalized distances
 * - Computes a weighted average of response values
 *
 * This adaptive approach balances global bandwidth optimization (from cross-validation)
 * with local adaptivity to ensure stable predictions across the entire graph, including
 * in sparse regions where the optimal bandwidth might be insufficient.
 *
 * @param bandwidth Optimal bandwidth value determined by cross-validation
 * @param y Vector of response values at each vertex in the graph
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 * @param lower_bound Minimum value for bandwidth expansion (typically a small fraction of graph diameter)
 * @param upper_bound Maximum value for bandwidth expansion (typically a large fraction of graph diameter)
 * @param domain_min_size Minimum number of neighbors required for a stable prediction
 *
 * @return std::vector<double> Vector of predictions for all vertices in the graph
 *
 * @throws std::runtime_error If no valid neighbors can be found even after bandwidth expansion
 *
 * @note The function first attempts to use the optimal bandwidth determined by cross-validation,
 *       only expanding it when necessary to ensure prediction stability
 * @note The bandwidth expansion is bounded by lower_bound and upper_bound to prevent excessive search
 * @note The distance normalization approach (shifting by minimum distance and normalizing)
 *       ensures consistency with the cross-validation procedure
 * @note For uniform weighting (use_uniform_weights = true), all neighbors contribute equally
 *       regardless of distance; otherwise, kernel weights decrease with distance
 *
 * @see find_minimum_radius_for_domain_min_size For the algorithm that determines the
 *      minimum bandwidth needed to include domain_min_size neighbors
 * @see predict_test_vertex_with_buffer For the compatible prediction method used during cross-validation
 */
std::vector<double> set_wgraph_t::predict_all_with_optimal_bandwidth(
    double bandwidth,
    const std::vector<double>& y,
    double dist_normalization_factor,
    bool use_uniform_weights,
    double lower_bound,
    double upper_bound,
    size_t domain_min_size
    ) {

    size_t n_vertices = adjacency_list.size();
    std::vector<double> predictions(n_vertices);


    double precision = 1e-6;
    // Process each vertex
    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

        // Find vertices within bandwidth radius
        auto vertex_map = find_vertices_within_radius(vertex, bandwidth);

        if (vertex_map.size() < domain_min_size) {
            double min_bw = find_minimum_radius_for_domain_min_size(
                vertex,
                lower_bound,
                upper_bound,
                domain_min_size,
                precision);

            vertex_map = find_vertices_within_radius(vertex, min_bw);
            // Note that vertex_map.empty() = true can never happen as 'vertex' is a part of the map so vertex_map.size() >= 1
        }

        // Convert map to vectors for easier processing
        std::vector<size_t> neighbors;
        std::vector<double> distances;

        for (const auto& [neighbor, distance] : vertex_map) {
            neighbors.push_back(neighbor);
            distances.push_back(distance);
        }

        // Calculate weights
        std::vector<double> weights(neighbors.size());

        if (use_uniform_weights) {
            // Use uniform weights
            std::fill(weights.begin(), weights.end(), 1.0);
        } else {
            // First, find the minimum distance for shifting
            double min_dist = *std::min_element(distances.begin(), distances.end());

            // Prepare normalized distances vector
            std::vector<double> normalized_dists(neighbors.size());

            // Shift distances so the smallest is 0, then normalize
            for (size_t i = 0; i < distances.size(); ++i) {
                normalized_dists[i] = distances[i] - min_dist;   // Shift by minimum
            }

            // Find the new maximum after shifting
            double max_dist = *std::max_element(normalized_dists.begin(), normalized_dists.end());
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            // Normalize by the scaled maximum
            for (size_t i = 0; i < normalized_dists.size(); ++i) {
                normalized_dists[i] /= max_dist;
            }

            // Use the kernel function to calculate weights
            kernel_fn(normalized_dists.data(), static_cast<int>(neighbors.size()), weights.data());
        }

        // Compute weighted average
        double weighted_sum = 0.0;
        double weight_sum = 0.0;

        for (size_t i = 0; i < neighbors.size(); ++i) {
            weighted_sum += weights[i] * y[neighbors[i]];
            weight_sum += weights[i];
        }

        // Set prediction
        if (weight_sum <= 0) {
            predictions[vertex] = y[vertex];  // Default to original value
        } else {
            predictions[vertex] = weighted_sum / weight_sum;
        }
    }

    return predictions;
}


/**
 * @brief Creates a buffer zone around specified vertices up to a given hop distance
 *
 * @details
 * This function identifies all vertices within a specified hop distance of any
 * input vertex, creating a "buffer zone" around these vertices. Buffer zones
 * are essential for spatial cross-validation to ensure independence between
 * training and test sets in graph-structured data.
 *
 * The function uses breadth-first search (BFS) to identify vertices at increasing
 * hop distances, up to the specified maximum. The buffer zone includes:
 * - The input vertices themselves
 * - All vertices that can be reached within buffer_hops steps from any input vertex
 *
 * @param vertices Vector of vertex indices around which to create the buffer zone
 * @param buffer_hops Maximum number of hops to include in the buffer zone
 *
 * @return std::unordered_set<size_t> Set of vertex indices that form the buffer zone
 *
 * @note Hop distance is defined as the minimum number of edges in the shortest path
 *       between two vertices
 * @note The function handles potentially overlapping buffer zones when multiple
 *       input vertices are provided
 * @note If buffer_hops is 0, only the input vertices themselves are included
 *       in the buffer zone
 * @note The size of the buffer zone grows exponentially with buffer_hops in
 *       most graph structures
 *
 * @see predict_test_vertex_with_buffer() For application in buffer zone cross-validation
 */
std::unordered_set<size_t> set_wgraph_t::create_buffer_zone(
    const std::vector<size_t>& vertices,
    size_t buffer_hops) {

    std::unordered_set<size_t> buffer_zone;

    // Add all input vertices to the buffer zone
    for (size_t vertex : vertices) {
        buffer_zone.insert(vertex);
    }

    // Use breadth-first search to expand buffer zone
    std::queue<std::pair<size_t, size_t>> bfs_queue;  // (vertex, hop_level)
    std::unordered_set<size_t> visited = buffer_zone;

    // Initialize queue with starting vertices
    for (size_t vertex : vertices) {
        bfs_queue.push({vertex, 0});
    }

    // BFS traversal
    while (!bfs_queue.empty()) {
        auto [current_vertex, hop_level] = bfs_queue.front();
        bfs_queue.pop();

        // If we've reached the maximum hop level, stop expanding
        if (hop_level >= buffer_hops) {
            continue;
        }

        // Process neighbors
        for (const auto& edge : adjacency_list[current_vertex]) {
            size_t neighbor = edge.vertex;

            // If not visited yet
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                buffer_zone.insert(neighbor);
                bfs_queue.push({neighbor, hop_level + 1});
            }
        }
    }

    return buffer_zone;
}

/**
 * @brief Predict values for test vertices while respecting buffer zones
 *
 * @details This function generates predictions for test vertices by:
 * 1. Finding all vertices within the bandwidth radius
 * 2. Removing any vertices that are in the buffer zone
 * 3. Using only the remaining vertices (outside the buffer zone) to make predictions
 * 4. If no suitable vertices remain for prediction, marking the prediction as excluded
 *
 * This approach ensures proper separation between training and test data by
 * preventing any information leakage through the graph structure.
 *
 * @param bandwidth Current bandwidth value to use for predictions
 * @param weights Vector of weights for all vertices (0.0 for test/buffer vertices, 1.0 for training)
 * @param y Response values at each vertex in the graph
 * @param test_vertices Vector of vertex indices that are in the test set
 * @param buffer_zone Set of vertex indices that form the buffer zone
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param use_uniform_weights Whether to use uniform weights instead of kernel weights
 *
 * @return std::vector<prediction_result_t> Vector of prediction results for test vertices
 *
 * @note If no valid neighbors outside the buffer zone are found for a test vertex,
 * the prediction is marked as excluded (is_excluded = true)
 */
std::vector<prediction_result_t> set_wgraph_t::predict_with_buffer_zone(
    double bandwidth,
    const std::vector<double>& weights,
    const std::vector<double>& y,
    const std::vector<size_t>& test_vertices,
    const std::unordered_set<size_t>& buffer_zone,
    double dist_normalization_factor,
    bool use_uniform_weights) {

    std::vector<prediction_result_t> predictions(test_vertices.size());

    // Process each test vertex
    for (size_t i = 0; i < test_vertices.size(); ++i) {
        size_t test_vertex = test_vertices[i];

        // Find vertices within bandwidth radius (excluding buffer zone)
        auto vertex_map = find_vertices_within_radius(test_vertex, bandwidth);

        // Remove buffer zone vertices from consideration
        std::vector<size_t> neighbors;
        std::vector<double> distances;

        for (const auto& [neighbor, distance] : vertex_map) {
            // Skip if in buffer zone
            if (buffer_zone.find(neighbor) != buffer_zone.end()) {
                continue;
            }

            // Only use vertices with non-zero weight (training vertices)
            if (weights[neighbor] > 0) {
                neighbors.push_back(neighbor);
                distances.push_back(distance);
            }
        }

        // Check if we have enough neighbors for prediction
        size_t nbhd_size = neighbors.size();
        if (nbhd_size == 0) {
            // No valid neighbors outside buffer zone for prediction
            predictions[i].value = std::numeric_limits<double>::quiet_NaN();
            predictions[i].is_excluded = true;
            continue;
        }

        std::vector<double> kernel_weights(nbhd_size);

        // Check if using uniform weights
        if (use_uniform_weights) {
            // Set all weights to 1/nbhd_size for uniform weighting
            double uniform_weight = 1.0;
            std::fill(kernel_weights.begin(), kernel_weights.end(), uniform_weight);
        } else {
            // Distance normalization strategy:
            // 1. Find the maximum distance (max_dist) in the vertex neighborhood
            // 2. Scale max_dist by dist_normalization_factor (typically 1.1)
            // 3. Divide all distances by this scaled maximum
            // This approach ensures all normalized distances fall within [0, 1/dist_normalization_factor],
            // which is approximately [0, 0.91] when dist_normalization_factor = 1.1
            // This keeps all distances within the effective support of the kernel functions,
            // as most kernels in this implementation have support on [-1, 1]

            double max_dist = 0.0;
            for (size_t k = 0; k < nbhd_size; ++k) {
                max_dist = std::max(max_dist, distances[k]);
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            // Normalize distances and calculate kernel weights
            std::vector<double> normalized_dists(nbhd_size);

            for (size_t k = 0; k < nbhd_size; ++k) {
                normalized_dists[k] = distances[k] / max_dist;
            }

            kernel_fn(normalized_dists.data(), nbhd_size, kernel_weights.data());
        }

        // Compute weighted average
        double weighted_sum = 0.0;
        double weight_sum = 0.0;

        for (size_t k = 0; k < nbhd_size; ++k) {
            weighted_sum += kernel_weights[k] * y[neighbors[k]];
            weight_sum += kernel_weights[k];
        }

        // Set prediction
        if (weight_sum <= 0) {
            predictions[i].value = std::numeric_limits<double>::quiet_NaN();
            predictions[i].is_excluded = true;
        } else {
            predictions[i].value = weighted_sum / weight_sum;
            predictions[i].is_excluded = false;
        }
    }

    return predictions;
}

/**
 * @brief Determine the optimal buffer hop distance based on spatial autocorrelation
 *
 * @details This function analyzes the spatial autocorrelation in the data at
 * different hop distances using Moran's I statistic. It determines the optimal
 * buffer size as either:
 * 1. The hop distance where autocorrelation drops below a significance threshold, or
 * 2. The hop distance where autocorrelation shows the largest decrease
 *
 * This provides a data-driven approach to setting the buffer size based on the
 * actual spatial dependence structure in the data.
 *
 * @param y Response values at each vertex in the graph
 * @param verbose Whether to print progress information
 *
 * @return size_t The optimal number of hops for the buffer zone
 *
 * @note The function checks up to a maximum of 10 hop distances and uses
 * a threshold of 0.1 for determining significant autocorrelation.
 */
size_t set_wgraph_t::determine_optimal_buffer_hops(
    const std::vector<double>& y,
    bool verbose) {

    if (verbose) {
        Rprintf("Analyzing spatial autocorrelation to determine optimal buffer size...\n");
    }

    // Calculate Moran's I at different hop distances
    std::vector<double> morans_i_values;
    const size_t max_hops_to_check = 10;   // Reasonable upper limit

    for (size_t hop = 1; hop <= max_hops_to_check; ++hop) {
        double morans_i = calculate_morans_i(y, hop);
        morans_i_values.push_back(morans_i);

        if (verbose) {
            Rprintf("  Hop distance %zu: Moran's I = %.4f\n", hop, morans_i);
        }
    }

    // Find the hop distance where autocorrelation drops below significance threshold
    // A common threshold is 0.1 or when the value becomes statistically non-significant
    const double autocorr_threshold = 0.1;
    size_t optimal_hops = 1;   // Default to 1 if no suitable value found

    for (size_t hop = 0; hop < morans_i_values.size(); ++hop) {
        if (std::abs(morans_i_values[hop]) < autocorr_threshold) {
            optimal_hops = hop + 1;   // +1 because hop is 0-indexed
            break;
        }
    }

    // If no value drops below threshold, take the hop distance where
    // autocorrelation has reduced the most significantly
    if (optimal_hops == 1 && morans_i_values.size() > 2) {
        std::vector<double> diff_values(morans_i_values.size() - 1);
        for (size_t i = 0; i < diff_values.size(); ++i) {
            diff_values[i] = std::abs(morans_i_values[i] - morans_i_values[i+1]);
        }

        size_t max_diff_idx = std::max_element(diff_values.begin(), diff_values.end()) - diff_values.begin();
        optimal_hops = max_diff_idx + 2;   // +2 because we want the hop after the big drop and diff_idx is 0-indexed
    }

    return optimal_hops;
}

/**
 * @brief Calculate Moran's I spatial autocorrelation statistic
 *
 * @details Moran's I is a measure of spatial autocorrelation that indicates
 * how similar observations are based on their spatial proximity. This function:
 * 1. Constructs a spatial weights matrix (W) based on the specified hop distance
 * 2. Calculates the mean of the response variable
 * 3. Computes Moran's I statistic using the formula:
 *    I = (N/W_sum) * (Σᵢ Σⱼ wᵢⱼ(yᵢ-ȳ)(yⱼ-ȳ)) / (Σᵢ(yᵢ-ȳ)²)
 *
 * Values of I range from -1 (perfect dispersion) to 1 (perfect correlation),
 * with 0 indicating no spatial autocorrelation.
 *
 * @param y Response values at each vertex in the graph
 * @param hop_distance The hop distance to use for defining neighborhood relationships
 *
 * @return double Moran's I statistic value
 *
 * @note The function uses a breadth-first search approach to identify vertices
 * at exactly hop_distance away from each vertex.
 */
double set_wgraph_t::calculate_morans_i(
    const std::vector<double>& y,
    size_t hop_distance) {

    size_t n_vertices = adjacency_list.size();

    // Calculate mean of y
    double y_mean = 0.0;
    for (size_t i = 0; i < n_vertices; ++i) {
        y_mean += y[i];
    }
    y_mean /= n_vertices;

    // Create adjacency matrix for specified hop distance
    std::vector<std::vector<double>> W(n_vertices, std::vector<double>(n_vertices, 0.0));

    // For each vertex, find all vertices exactly hop_distance away
    for (size_t i = 0; i < n_vertices; ++i) {
        std::unordered_set<size_t> current_level = {i};
        std::unordered_set<size_t> next_level;
        std::unordered_set<size_t> visited = {i};

        // BFS to find vertices at the specified hop distance
        for (size_t hop = 1; hop <= hop_distance; ++hop) {
            next_level.clear();

            for (size_t vertex : current_level) {
                for (const auto& edge : adjacency_list[vertex]) {
                    size_t neighbor = edge.vertex;
                    if (visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        next_level.insert(neighbor);
                    }
                }
            }

            current_level = next_level;

            // If we've reached the target hop distance, set adjacency values
            if (hop == hop_distance) {
                for (size_t neighbor : current_level) {
                    W[i][neighbor] = 1.0;
                }
            }
        }
    }

    // Calculate Moran's I
    double numerator = 0.0;
    double denominator = 0.0;
    double W_sum = 0.0;

    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t j = 0; j < n_vertices; ++j) {
            if (i != j) {
                numerator += W[i][j] * (y[i] - y_mean) * (y[j] - y_mean);
                W_sum += W[i][j];
            }
        }
        denominator += (y[i] - y_mean) * (y[i] - y_mean);
    }

    if (W_sum == 0.0 || denominator == 0.0) {
        return 0.0;   // No spatial relationship or no variance
    }

    return (n_vertices / W_sum) * (numerator / denominator);
}


std::vector<double> set_wgraph_t::get_trainig_predictions(
    const std::vector<double>& y,
    const std::vector<double>& training_vertices,
    double bandwidth) const {

    std::vector<double> trainig_predictions(training_vertices.size());
    for (size_t vertex : trainig_vertices) {

        // 1. Find vertices within bandwidth radius
        auto vertex_map = find_vertices_within_radius(vertex, bandwidth);

        // 2. If the bandwidth is too small, increase it to find
        if (vertex_map.empty()) {
            double min_bw = find_minimum_radius_for_domain_min_size(
                vertex,
                lower_bound,
                upper_bound,
                domain_min_size,
                precision);

            vertex_map = find_vertices_within_radius(vertex, min_bw);
            if (vertex_map.empty()) {
                // we should never get here !!!
                REPORT_ERROR("ERROR: Even after expansion of the bandwidth we still have no vertices within the new bandwidth radius\n");
            }
        }

}

}
