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
#include <filesystem>
#include <experimental/filesystem>

#include "graph_kernel_smoother.hpp" // For graph_kernel_smoother_t
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
 * @param vertex_hbhd_min_size The minimal number of vertices that the smallest neighborhood of the given vertex has to have
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
graph_kernel_smoother_t set_wgraph_t::graph_kernel_smoother(
    const std::vector<double>& y,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool log_grid,
    size_t vertex_hbhd_min_size,
    size_t kernel_type,
    double dist_normalization_factor,
    bool use_uniform_weights,
    size_t buffer_hops,
    bool auto_buffer_hops,
    size_t n_folds,
    bool with_bw_predictions,
    double precision,
    bool verbose) {

    #define DEBUG__graph_kernel_smoother 0

    initialize_kernel(kernel_type, 1.0);

    size_t ncc = count_connected_components();
    if (ncc > 1) {
        REPORT_ERROR("The graph has to have only one connected component!!!\n");
    }

    size_t n_vertices = adjacency_list.size();

    std::vector<std::vector<size_t>> folds;
    if (n_folds >= n_vertices) { // LOO CV case

        n_folds = n_vertices;
        folds.resize(n_folds);
        for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
            folds[vertex] = {vertex};
        }
        // Estimate the diameter of the graph
        auto [end1, diam] = get_vertex_eccentricity(0);  // Start from vertex 0
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

    // Initialize result structure
    graph_kernel_smoother_t result;
    result.predictions.resize(n_vertices);          // kernel smoothing predictions for each vertex
    result.bw_mean_sq_errors.resize(n_bws, 0.0);    // mean_bw_sq_errors[bw_idx] is the mean prediction squared error for the bw_idx-th bandwidth
    result.bw_mean_abs_errors.resize(n_bws, 0.0);   // mean_bw_abs_errors[bw_idx] is the mean prediction absolute error for the bw_idx-th bandwidth
    result.vertex_min_bws.resize(n_vertices);       // vertex_min_bws[i] is the the i-th vertex's minimum bandwidth
    result.bw_lextr_count.resize(n_bws);
    result.bw_total_sq_curvature.resize(n_bws, 0.0);
    result.bw_total_sq_norma_curvature.resize(n_bws, 0.0);

    if (with_bw_predictions) {
        result.bw_predictions.resize(n_bws,
                                     std::vector<double>(n_vertices, 0.0));
    }

    // Auto-determine buffer hops if requested
    if (auto_buffer_hops) {
        buffer_hops = determine_optimal_buffer_hops(y, verbose);
    }
    result.buffer_hops_used = buffer_hops;

    if (verbose) {
        Rprintf("Starting graph_kernel_smoother() with Buffer Zone Method\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Number of CV folds: %zu\n", n_folds);
        Rprintf("Buffer zone hops: %zu\n", buffer_hops);
        Rprintf("min_bw_factor: %f\n", min_bw_factor);
        Rprintf("max_bw_factor: %f\n", max_bw_factor);
        Rprintf("min_bw: %f\n", min_bw);
        Rprintf("max_bw: %f\n", max_bw);
        Rprintf("graph_diameter: %f\n", graph_diameter);
    }

    std::vector<size_t> neighbors(n_vertices);
    std::vector<double> distances(n_vertices);

    // -----------------------------------------------------------------------------
    // Cross-validation section
    // -----------------------------------------------------------------------------

    // debugging
    std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data/";
    if (!std::filesystem::exists(debug_dir)) {
        if (!std::filesystem::create_directories(debug_dir)) {
            REPORT_ERROR("ERROR: Failed to create debug directory: %s\n", debug_dir.c_str());
        }
    }

    for (size_t fold = 0; fold < folds.size(); ++fold) {

        const std::vector<size_t>& test_vertices = folds[fold];
        std::unordered_set<size_t> test_set(
            test_vertices.begin(),
            test_vertices.end()
            );

        // 1. Create buffer zone arond test vertices
        std::unordered_set<size_t> buffer_zone = create_buffer_zone(
            test_vertices,
            buffer_hops
            );

        // 2. Training vertices are the vertices of the graph that do not belong to buffer zone
        std::unordered_set<size_t> training_set;
        training_set.reserve(n_vertices - buffer_zone.size());  // Pre-allocate for efficiency

        for (size_t i = 0; i < n_vertices; ++i) {
            if (buffer_zone.find(i) == buffer_zone.end()) {
                training_set.insert(i);
            }
        }

#if DEBUG__graph_kernel_smoother
        Rprintf("fold: %zu\n", fold);
        print_vect(test_vertices, "test_vertices");
        print_uset(buffer_zone, "buffer_zone");
        //print_uset(training_set, "training_set");
        // error("DEBUGGING\n");
#endif

        // 3. Generate predictions fot the test set vertices only given kernel mean given the current bandwidth
        for (const auto& vertex : test_set) {

            // 1. compute neighbor distance map of the given vertex for the maximum bw
            std::unordered_map<size_t, double> ngbr_dist_map = find_vertices_within_radius(vertex, max_bw);

            // 2. Use ngbr_dist_map to find vertex specific min_bw
            auto [sorted_vertices, vertex_min_bw] = get_sorted_vertices_and_min_radius(
                ngbr_dist_map,
                vertex_hbhd_min_size,
                training_set
                );

            if (vertex_min_bw < min_bw) {
                vertex_min_bw = min_bw;
            }

            // 3. Generate candidate bandwidths
            std::vector<double> candidate_bws = get_candidate_bws(
                vertex_min_bw,
                max_bw,
                n_bws,
                log_grid,
                precision
                );

            // 4. For each bw generate kernel-weighted mean prediction at the given vertex
            for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

                double current_bw = candidate_bws[bw_idx];

                size_t nbhd_size = 0;
                for (const auto& [v, dist] : sorted_vertices) {
                    if (dist >= current_bw)
                        break;
                    neighbors[nbhd_size] = v;
                    distances[nbhd_size] = dist;
                    nbhd_size++;
                }

                // Distance normalization strategy:
                // 1. Find the maximum distance (max_dist) in the vertex neighborhood
                // 2. Scale max_dist by dist_normalization_factor (typically 1.1)
                // 3. Divide all distances by this scaled maximum
                // This approach ensures all normalized distances fall within [0, 1/dist_normalization_factor],
                // which is approximately [0, 0.91] when dist_normalization_factor = 1.1
                // This keeps all distances within the effective support of the kernel functions,
                // as most kernels in this implementation have support on [-1, 1]

                double max_dist = distances[nbhd_size - 1];
                if (max_dist == 0) max_dist = 1.0;

                // Calculate prediction based on training vertices
                std::vector<double> weights(nbhd_size, 1.0);

                if (!use_uniform_weights) {

                    max_dist *= dist_normalization_factor;

                    // Normalize by the scaled maximum
                    for (size_t i = 0; i < nbhd_size; ++i) {
                        distances[i] /= max_dist;
                    }

                    // Use the kernel function to calculate weights
                    kernel_fn(distances.data(), static_cast<int>(nbhd_size), weights.data());
                }

                // Compute weighted average
                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                for (size_t i = 0; i < nbhd_size; ++i) {
                    weighted_sum += weights[i] * y[neighbors[i]];
                    weight_sum += weights[i];
                }

                if (weight_sum <= 0) {
                    // Based on how we defined normalized distances this should never happen
                    Rprintf("weight_sum: %.4f\n", weight_sum);
                    REPORT_ERROR("Undefined prediction value for vertex: %zu\n", vertex);
                }

                //test_bw_predictions[bw_idx][vertex] = weighted_sum / weight_sum;
                double prediction = weighted_sum / weight_sum;

                result.bw_mean_sq_errors[bw_idx]  += std::pow(std::abs(prediction - y[vertex]), 2.0);
                result.bw_mean_abs_errors[bw_idx] += std::abs(prediction - y[vertex]);

            } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)

        } // END OF for (const auto& vertex : training_set)

#if DEBUG__graph_kernel_smoother
        // debugging
        {
            // Dump training predictions to CSV
            std::filesystem::path file_path = debug_dir;
            file_path /= "training_predictions_fold_" + std::to_string(fold) + ".csv";

            std::ofstream ofs(file_path);
            if (!ofs.is_open()) {
                REPORT_ERROR("ERROR: Failed to open file for writing: %s\n", file_path.string().c_str());
            }

            // Optional: write a header row
            ofs << "vertex";
            for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
                ofs << ",bw_" << bw_idx;
            }
            ofs << "\n";

            // Write one row per vertex
            for (size_t vertex_idx = 0; vertex_idx < n_vertices; ++vertex_idx) {
                ofs << vertex_idx;
                for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
                    ofs << "," << training_bw_predictions[bw_idx][vertex_idx];
                }
                ofs << "\n";
            }

            ofs.close();

        }

        #if 0
        // Generate predictions over test vertices and estimate prediction errors
        for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
            for (const auto& test_vertex : test_vertices) {
                double prediction = generate_test_predictions_from_training(
                    test_vertex,
                    buffer_hops,
                    training_bw_predictions[bw_idx],
                    training_set,
                    dist_normalization_factor,
                    use_uniform_weights
                    );
                result.bw_mean_sq_errors[bw_idx] += std::pow(prediction - y[test_vertex], 2.0);
                result.bw_mean_abs_errors[bw_idx] += std::abs(prediction - y[test_vertex]);
            }
        }
        #endif

        // Collect test predictions in a map for debugging
        // -----------------------------------------------
        // key = vertex index, value = vector of length n_bws holding its per‑bw predictions
        std::unordered_map<size_t, std::vector<double>> test_bw_predictions;
        test_bw_predictions.reserve(test_vertices.size());

        // Initialize each test vertex’s vector
        for (auto v : test_vertices) {
            test_bw_predictions[v].assign(n_bws, 0.0);
        }

        // Generate & store predictions (and accumulate errors)
        for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
            for (auto test_vertex : test_vertices) {
                double prediction = generate_test_predictions_from_training(
                    test_vertex,
                    buffer_hops,
                    training_bw_predictions[bw_idx],
                    training_set,
                    dist_normalization_factor,
                    use_uniform_weights
                    );
                // store for later dump
                test_bw_predictions[test_vertex][bw_idx] = prediction;

                // accumulate your error metrics
                result.bw_mean_sq_errors[bw_idx]  += std::pow(std::abs(prediction - y[test_vertex]), 0.5);
                result.bw_mean_abs_errors[bw_idx] += std::abs(prediction - y[test_vertex]);
            }
        }

        // Dump test predictions to CSV
        // -----------------------------------------------
        std::filesystem::path test_file = debug_dir;
        test_file /= "test_predictions_fold_" + std::to_string(fold) + ".csv";

        std::ofstream ofs_test(test_file);
        if (!ofs_test.is_open()) {
            REPORT_ERROR("ERROR: Failed to open file for writing: %s\n",
                         test_file.string().c_str());
        }

        // header
        ofs_test << "vertex";
        for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
            ofs_test << ",bw_" << bw_idx;
        }
        ofs_test << "\n";

        // one row per test vertex (in original order)
        for (auto test_vertex : test_vertices) {
            ofs_test << test_vertex;
            for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
                ofs_test << "," << test_bw_predictions[test_vertex][bw_idx];
            }
            ofs_test << "\n";
        }

        ofs_test.close();
        //error("DEBUGGING\n");

#endif

    } // END OF for (size_t fold = 0; fold < folds.size(); ++fold)

    // Assuming the total number of test vertices was n_vertices
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        result.bw_mean_sq_errors[bw_idx]  /= n_vertices;
        result.bw_mean_abs_errors[bw_idx] /= n_vertices;
    }

    // -----------------------------------------------------------------------------
    // All data prediction section
    // -----------------------------------------------------------------------------
    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

        // 1. compute neighbor distance map of the given vertex for the maximum bw
        std::unordered_map<size_t, double> ngbr_dist_map = find_vertices_within_radius(vertex, max_bw);

        // 2. Use ngbr_dist_map to find vertex specific min_bw
        auto [sorted_vertices, vertex_min_bw] = get_sorted_vertices_and_min_radius(
            ngbr_dist_map,
            vertex_hbhd_min_size
            );

        if (vertex_min_bw < min_bw) {
            vertex_min_bw = min_bw;
        }

        result.vertex_min_bws[vertex] = vertex_min_bw;

        // 3. Generate candidate bandwidths
        std::vector<double> candidate_bws = get_candidate_bws(
            vertex_min_bw,
            max_bw,
            n_bws,
            log_grid,
            precision
            );

        // 4. For each bw generate kernel-weighted mean prediction at the given vertex
        for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {

            double current_bw = candidate_bws[bw_idx];

            size_t nbhd_size = 0;
            for (const auto& [v, dist] : sorted_vertices) {
                if (dist >= current_bw)
                    break;
                neighbors[nbhd_size] = v;
                distances[nbhd_size] = dist;
                nbhd_size++;
            }

            // Distance normalization strategy:
            // 1. Find the maximum distance (max_dist) in the vertex neighborhood
            // 2. Scale max_dist by dist_normalization_factor (typically 1.1)
            // 3. Divide all distances by this scaled maximum
            // This approach ensures all normalized distances fall within [0, 1/dist_normalization_factor],
            // which is approximately [0, 0.91] when dist_normalization_factor = 1.1
            // This keeps all distances within the effective support of the kernel functions,
            // as most kernels in this implementation have support on [-1, 1]

            double max_dist = distances[nbhd_size - 1];
            if (max_dist == 0) max_dist = 1.0;

            // Calculate prediction based on training vertices
            std::vector<double> weights(nbhd_size);

            if (use_uniform_weights) {
                // Use uniform weights
                std::fill(weights.begin(), weights.end(), 1.0);
            } else {
                max_dist *= dist_normalization_factor;

                // Normalize by the scaled maximum
                for (size_t i = 0; i < nbhd_size; ++i) {
                    distances[i] /= max_dist;
                }

                // Use the kernel function to calculate weights
                kernel_fn(distances.data(), static_cast<int>(nbhd_size), weights.data());
            }

            // Compute weighted average
            double weighted_sum = 0.0;
            double weight_sum = 0.0;
            for (size_t i = 0; i < nbhd_size; ++i) {
                weighted_sum += weights[i] * y[neighbors[i]];
                weight_sum += weights[i];
            }

            if (weight_sum <= 0) {
                // Based on how we defined normalized distances this should never happen
                Rprintf("weight_sum: %.4f\n", weight_sum);
                REPORT_ERROR("Undefined prediction value for vertex: %zu\n", vertex);
            }

            result.bw_predictions[bw_idx][vertex] = weighted_sum / weight_sum;

        } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)

    } // END OF for (size_t vertex = 0; vertex < n_vertices; ++vertex)

    // Compute the number of local extrema and different estimates of total curvature
    double delta = 1e-6 * compute_median_edge_length();
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        result.bw_lextr_count[bw_idx] = get_lextr_count(result.bw_predictions[bw_idx]);
        result.bw_total_sq_curvature[bw_idx] = get_total_sq_curvature(result.bw_predictions[bw_idx]);
        result.bw_total_sq_norma_curvature[bw_idx] = get_total_sq_normalized_curvature(result.bw_predictions[bw_idx], delta);
    }

    // Find optimal bandwidth using result.bw_mean_sq_errors
    auto min_it = std::min_element(result.bw_mean_abs_errors.begin(), result.bw_mean_abs_errors.end());
    result.opt_bw_idx = min_it - result.bw_mean_abs_errors.begin();
    result.predictions = result.bw_predictions[result.opt_bw_idx];

    return result;
}

/**
 * @brief Find the minimum radius needed to include a specified number of neighbor vertices
 *
 * @details
 * This function determines the minimum radius (bandwidth) required to include at least
 * vertex_hbhd_min_size neighbors around a vertex. It sorts the neighbors by distance
 * and returns the distance to the (vertex_hbhd_min_size + 1)th vertex. If there aren't
 * enough neighbors, it returns the maximum available distance.
 *
 * @param ngbr_dist_map Map containing vertex indices and their distances from a reference vertex
 * @param vertex_hbhd_min_size Minimum number of neighbors required in the vertex neighborhood
 *
 * @return double The minimum radius (bandwidth) needed to include at least vertex_hbhd_min_size neighbors.
 *         Returns the maximum available distance if there aren't enough neighbors.
 */
double set_wgraph_t::find_min_radius_for_neighbor_count(
    std::unordered_map<size_t, double>& ngbr_dist_map,
    size_t vertex_hbhd_min_size
    ) const {
    std::vector<std::pair<size_t, double>> sorted_vertices;
    sorted_vertices.reserve(ngbr_dist_map.size());

    // Convert map to vector for sorting
    for (const auto& [v, dist] : ngbr_dist_map) {
        sorted_vertices.emplace_back(v, dist);
    }

    // Sort by distance
    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    // Find minimum bandwidth that includes at least vertex_hbhd_min_size neighbors
    if (sorted_vertices.size() <= vertex_hbhd_min_size) {
        // Not enough vertices, return the maximum available distance
        return sorted_vertices.empty() ? 0.0 : sorted_vertices.back().second;
    }

    // Return the distance to the (vertex_hbhd_min_size + 1)th vertex
    return sorted_vertices[vertex_hbhd_min_size].second;
}

/**
 * @brief Get vertices sorted by distance and find minimum radius for neighborhood size
 *
 * @details
 * This function sorts vertices by their distance from a reference vertex and
 * determines the minimum radius needed to include at least vertex_hbhd_min_size neighbors.
 * It returns both the sorted vertices and the minimum radius value.
 *
 * @param ngbr_dist_map Map containing vertex indices and their distances from a reference vertex
 * @param vertex_hbhd_min_size Minimum number of neighbors required in the vertex neighborhood
 *
 * @return std::pair<std::vector<std::pair<size_t, double>>, double>
 *         First element: Vector of vertex-distance pairs sorted by distance
 *         Second element: Minimum radius needed to include at least vertex_hbhd_min_size neighbors
 */
std::pair<std::vector<std::pair<size_t, double>>, double>
set_wgraph_t::get_sorted_vertices_and_min_radius(
    const std::unordered_map<size_t, double>& ngbr_dist_map,
    size_t vertex_hbhd_min_size
    ) const {
    std::vector<std::pair<size_t, double>> sorted_vertices;
    sorted_vertices.reserve(ngbr_dist_map.size());

    // Convert map to vector for sorting
    for (const auto& [v, dist] : ngbr_dist_map) {
        sorted_vertices.emplace_back(v, dist);
    }

    // Sort by distance
    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    // Find minimum bandwidth that includes at least vertex_hbhd_min_size neighbors
    double min_radius = 0.0;
    if (sorted_vertices.size() <= vertex_hbhd_min_size) {
        // Not enough vertices, return the maximum available distance
        REPORT_ERROR("Not enough vertices found in ngbr_dist_map\n");
        //min_radius = sorted_vertices.empty() ? 0.0 : sorted_vertices.back().second;
    } else {
        // Return the distance to the (vertex_hbhd_min_size + 1)th vertex
        min_radius = sorted_vertices[vertex_hbhd_min_size].second;
    }

    return {sorted_vertices, min_radius};
}

/**
 * @brief Get training vertices sorted by distance and find minimum radius for neighborhood size
 *
 * @details
 * This function filters vertices from ngbr_dist_map that exist in the training_set,
 * sorts them by their distance from a reference vertex, and determines the minimum
 * radius needed to include at least vertex_hbhd_min_size neighbors.
 * It returns both the filtered sorted vertices and the minimum radius value.
 *
 * @param ngbr_dist_map Map containing vertex indices and their distances from a reference vertex
 * @param vertex_hbhd_min_size Minimum number of neighbors required in the vertex neighborhood
 * @param training_set Set of vertices that should be considered for inclusion
 *
 * @return std::pair<std::vector<std::pair<size_t, double>>, double>
 *         First element: Vector of vertex-distance pairs filtered to training set and sorted by distance
 *         Second element: Minimum radius needed to include at least vertex_hbhd_min_size neighbors
 *
 * @throws Reports an error if the number of filtered vertices is insufficient to meet vertex_hbhd_min_size
 */
std::pair<std::vector<std::pair<size_t, double>>, double>
set_wgraph_t::get_sorted_vertices_and_min_radius(
    const std::unordered_map<size_t, double>& ngbr_dist_map,
    size_t vertex_hbhd_min_size,
    std::unordered_set<size_t>& training_set
    ) const {

    std::vector<std::pair<size_t, double>> sorted_vertices;
    sorted_vertices.reserve(ngbr_dist_map.size());

    // Convert map to vector for sorting, only including vertices in training_set
    for (const auto& [v, dist] : ngbr_dist_map) {
        if (training_set.find(v) != training_set.end()) {
            sorted_vertices.emplace_back(v, dist);
        }
    }

    // Sort by distance
    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    // Find minimum bandwidth that includes at least vertex_hbhd_min_size neighbors
    double min_radius = 0.0;
    if (sorted_vertices.size() <= vertex_hbhd_min_size) {
        // Not enough vertices, return the maximum available distance
        REPORT_ERROR("Not enough vertices found in ngbr_dist_map\n");
    } else {
        // Return the distance to the (vertex_hbhd_min_size + 1)th vertex
        min_radius = sorted_vertices[vertex_hbhd_min_size].second;
    }

    return {sorted_vertices, min_radius};
}


/**
 * @brief Count the number of local extrema (maxima or minima) of a function over the graph
 *
 * @details
 * This function identifies and counts vertices that are local extrema of a given function
 * represented by values at each vertex. A vertex is considered:
 * - A local maximum if its value is greater than or equal to all its neighbors' values
 * - A local minimum if its value is less than or equal to all its neighbors' values
 *
 * Each vertex is compared with its immediate neighbors according to the graph's adjacency list.
 * Isolated vertices (with no neighbors) are always considered both local maxima and minima.
 *
 * @param y Vector of function values at each vertex of the graph
 *
 * @return size_t The total count of vertices that are either local maxima or local minima
 *
 * @note If a vertex has the same value as all its neighbors, it will be counted as both
 *       a local maximum and a local minimum, but will only contribute once to the total count.
 *
 * @see detect_local_extrema() For a function that returns the actual extrema vertices
 */
size_t set_wgraph_t::get_lextr_count(
    std::vector<double>& y
    ) const {
    size_t n_vertices = adjacency_list.size();
    size_t lextr_count = 0;
    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
        const auto& nbhrs = adjacency_list[vertex];
        bool is_lmax = true;
        for (const auto& edge_info : nbhrs) {
            size_t nbhr = edge_info.vertex;
            if (y[vertex] < y[nbhr]) {
                is_lmax = false;
                break;
            }
        }
        bool is_lmin = true;
        for (const auto& edge_info : nbhrs) {
            size_t nbhr = edge_info.vertex;
            if (y[vertex] > y[nbhr]) {
                is_lmin = false;
                break;
            }
        }
        if (is_lmax || is_lmin) {
            lextr_count++;
        }
    }
    return lextr_count;
}

/**
 * @brief  Multi‐target Dijkstra: shortest‐path distances from `from` to each v∈toSet.
 * @param  from    Source vertex index
 * @param  to_set   Set of target vertices
 * @return Map v -> distance (INFINITY if unreachable)
 */
std::unordered_map<size_t,double>
set_wgraph_t::compute_shortest_path_distances(
    size_t from,
    const std::unordered_set<size_t>& to_set
) const {
    // 1. Standard Dijkstra distance array
    std::vector<double> dist(adjacency_list.size(),
                             std::numeric_limits<double>::infinity());
    dist[from] = 0.0;

    // 2. Prepare result map and "remaining" copy of targets
    std::unordered_map<size_t,double> result;
    result.reserve(to_set.size());
    for (auto v : to_set)
        result[v] = std::numeric_limits<double>::infinity();

    // Remove source if it’s a target
    std::unordered_set<size_t> remaining = to_set;
    if (remaining.erase(from))
        result[from] = 0.0;

    // 3. Min‐heap of (distance, vertex)
    using pq_item = std::pair<double, size_t>;
    std::priority_queue<pq_item,
        std::vector<pq_item>,
        std::greater<pq_item>>
      pq;
    pq.push({0.0, from});

    // 4. Main loop: stop when all targets found
    while (!pq.empty() && !remaining.empty()) {
        auto [d,u] = pq.top();
        pq.pop();

        // Skip stale entries
        if (d > dist[u])
            continue;

        // If u is one of our targets, record & remove it
        if (remaining.erase(u)) {
            result[u] = d;
        }

        // Relax neighbors
        for (auto const& edge : adjacency_list[u]) {
            size_t v    = edge.vertex;
            double alt  = d + edge.weight;
            if (alt < dist[v]) {
                dist[v] = alt;
                pq.push({alt, v});
            }
        }
    }

    return result;
}

/**
 * @brief  Predicts y[test_vertex] from its nearest training vertices beyond a buffer zone.
 *
 * 1) Finds the set S of training vertices at the smallest hop‐distance > buffer_hops.
 * 2) Computes true graph‐distances d(v) for v∈S via multi‐target Dijkstra.
 * 3) Computes weights w(v) either uniformly or by applying kernel_fn to normalized d(v).
 * 4) Returns ∑_{v∈S} w(v)·y[v] / ∑_{v∈S} w(v).
 *
 * @param test_vertex             Index of the test vertex.
 * @param buffer_hops             Minimum hop‐distance to exclude “nearby” vertices.
 * @param predictions             Response vector of training vertices predictions (size = |V|) that will be use to estimate the prediction at the given test vertex.
 * @param training_set            Set of allowed training vertices.
 * @param dist_normalization_factor  Factor to scale distances before kernel.
 * @param use_uniform_weights     If true, all w(v)=1.
 * @return                        Predicted value for y[test_vertex].
 */
double set_wgraph_t::generate_test_predictions_from_training(
    size_t                                   test_vertex,
    size_t                                   buffer_hops,
    const std::vector<double>&               predictions,
    const std::unordered_set<size_t>&        training_set,
    double                                   dist_normalization_factor,
    bool                                     use_uniform_weights
) const {

    #define DEBUG__generate_test_predictions_from_training 0

    // --- 1. BFS to find the minimal hop‐layer > buffer_hops that contains training vertices
    auto find_hop_layer = [this, &training_set](
        size_t start, size_t target_hop
    ) {
        std::unordered_set<size_t> visited{start};
        std::queue<std::pair<size_t,size_t>> q;
        q.push({start,0});

        std::unordered_set<size_t> layer;
        while (!q.empty()) {
            auto [u, hop] = q.front(); q.pop();
            if (hop == target_hop) {
                if (training_set.count(u))
                    layer.insert(u);
                continue;  // do not enqueue neighbors of this layer
            }
            for (auto const& e : adjacency_list[u]) {
                if (visited.insert(e.vertex).second) {
                    q.push({e.vertex, hop+1});
                }
            }
        }
        return layer;
    };

    size_t n_vertices = adjacency_list.size();
    size_t hop  = buffer_hops + 1;
    size_t max_hop = n_vertices;   // worst‑case diameter ≤ n_vertices-1 :contentReference[oaicite:0]{index=0}
    std::unordered_set<size_t> candidates;

    // fail fast if no trainers at all
    if (training_set.empty()) {
        REPORT_ERROR("Empty training_set for test vertex %zu\n", test_vertex);
    }

    // climb hops until we find some or hit max_hop
    while (hop <= max_hop) {
        candidates = find_hop_layer(test_vertex, hop);
        if (!candidates.empty()) break;
        ++hop;
    }

    #if DEBUG__generate_test_predictions_from_training
        Rprintf("\nIn generate_test_predictions_from_training()\ntest_vertex: %zu\n", test_vertex);
        Rprintf("buffer_hops: %zu\n", buffer_hops);
        Rprintf("hop: %zu\n", hop);
        print_uset(candidates, "candidates");
    #endif

    if (candidates.empty()) {
        REPORT_ERROR(
            "No reachable training vertices found within %zu hops "
            "for test vertex %zu\n",
            max_hop, test_vertex
            );
    }

    // --- 2. Compute true graph‐distances only to these candidates
    std::unordered_map<size_t,double> dist_map =
        compute_shortest_path_distances(test_vertex, candidates);

    // --- 3. Pack into vectors for weighting
    size_t n_candidates = candidates.size();
    std::vector<size_t> verts; verts.reserve(n_candidates);
    std::vector<double> dists; dists.reserve(n_candidates);
    for (auto v : candidates) {
        verts.push_back(v);
        dists.push_back(dist_map[v]);      // ∞ if unreachable (shouldn’t happen)
    }

    // --- 4. Build weights
    std::vector<double> weights(n_candidates, 1.0);

    #if DEBUG__generate_test_predictions_from_training
    Rprintf("use_uniform_weights: %s\n",
            use_uniform_weights ? "TRUE" : "FALSE");
    #endif

    if (!use_uniform_weights && n_candidates > 1) {
        // shift so min=0
        // double min_d = *std::min_element(dists.begin(), dists.end());
        // for (auto& dd : dists) dd -= min_d;
        // // scale by factor
        // double mx = *std::max_element(dists.begin(), dists.end());
        // if (mx == 0) mx = 1.0;
        // mx *= dist_normalization_factor;
        // for (auto& dd : dists) dd /= mx;

        double total_dist = std::accumulate(dists.begin(), dists.end(), 0.0);
        for (auto& dist : dists) dist /= total_dist;

        #if DEBUG__generate_test_predictions_from_training
        Rprintf("total_dist: %.3f\n", total_dist);
        print_vect(dists, "dists");
        #endif

        // apply kernel
        kernel_fn(dists.data(), static_cast<int>(n_candidates), weights.data());
    }

    // --- 5. Weighted average
    double num=0.0, den=0.0;
    for (size_t i = 0; i < n_candidates; ++i) {
        num += weights[i] * predictions[ verts[i] ];
        den += weights[i];
    }
    if (den <= 0) {
        REPORT_ERROR("Zero total weight for test vertex %zu\n", test_vertex);
    }

    #if DEBUG__generate_test_predictions_from_training
    print_umap(dist_map, "dist_map");
    print_vect(verts, "verts");
    print_vect(dists, "dists");
    print_vect(weights, "weights");
    Rprintf("num: %.3f  den: %.3f  prediction: %.3f\n", num, den, num / den);
    #endif

    #if DEBUG__generate_test_predictions_from_training
    if (test_vertex == 4) {
        error("DEBUGGING\n");
    }
    #endif

    return num / den;
}

/**
 * @brief Compute total squared curvature using the random‑walk normalized Laplacian on an unweighted graph.
 *
 * For each vertex i, defines the Laplacian
 *   Δf(i) = f(i) - (1/deg(i)) ∑_{j ∈ N(i)} f(j),
 * where deg(i) = |N(i)| and N(i) is the neighbor set of i.
 * Isolated vertices (deg = 0) contribute zero curvature.
 *
 * @param predictions  Vector of length |V| containing f(i) for each vertex.
 * @return             Sum over all vertices of [Δf(i)]².
 */
double set_wgraph_t::get_total_sq_curvature(
    const std::vector<double>& predictions
) const {
    double total_sq = 0.0;
    size_t n = adjacency_list.size();

    for (size_t i = 0; i < n; ++i) {
        const auto& nbrs = adjacency_list[i];
        size_t deg = nbrs.size();
        if (deg == 0) {
            continue;  // no neighbors → zero curvature
        }

        // Compute average prediction over neighbors
        double sum = 0.0;
        for (const auto& edge : nbrs) {
            sum += predictions[edge.vertex];
        }
        double avg = sum / static_cast<double>(deg);

        // Laplacian value and accumulate its square
        double lap = predictions[i] - avg;
        total_sq += lap * lap;
    }

    return total_sq;
}


/**
 * @brief Compute sum of squared normalized curvature (Laplacian) with additive regularization.
 *
 * Uses:
 *   w_{ij} = 1 / (edge.length + delta)
 *   d_i = \sum_j w_{ij}
 *   \Delta f(i) = f(i) - (1 / d_i) \sum_j w_{ij} f(j)
 *
 * @param predictions  Vector of size |V| with function values (f(i) for each vertex).
 * @param delta        Additive regularization parameter (>= 0) to avoid division by zero.
 * @return             Sum of squared normalized Laplacian values: \sum_i [\Delta f(i)]^2.
 */
double set_wgraph_t::get_total_sq_normalized_curvature(
    const std::vector<double>& predictions,
    double delta
) const {
    double total_sq = 0.0;
    size_t n_vertices = adjacency_list.size();

    for (size_t i = 0; i < n_vertices; ++i) {
        const auto& neighbors = adjacency_list[i];
        if (neighbors.empty()) {
            continue;  // Isolated vertex contributes zero curvature
        }

        // 1) compute weighted degree d_i = sum_j 1/(length_ij + delta)
        double degree_i = 0.0;
        for (const auto& edge : neighbors) {
            degree_i += 1.0 / (edge.weight + delta);
        }

        // 2) compute weighted neighbor average: sum_j w_ij * f(j) / d_i
        double weighted_sum = 0.0;
        for (const auto& edge : neighbors) {
            double w = 1.0 / (edge.weight + delta);
            weighted_sum += w * predictions[edge.vertex];
        }
        double neighbor_avg = weighted_sum / degree_i;

        // 3) Laplacian at i: f(i) - neighbor_avg
        double lap_val = predictions[i] - neighbor_avg;
        total_sq += lap_val * lap_val;
    }

    return total_sq;
}

/**
 * @brief Compute the median edge length in the graph.
 *
 * Collects all edge weights from adjacency_list and returns their median value.
 * Utilizes std::nth_element for O(E) average‐case performance.
 *
 * @return Median of all edge lengths, or 0.0 if the graph contains no edges.
 */
double set_wgraph_t::compute_median_edge_length() const {

    // 1) Gather all edge lengths into a single vector
    size_t total_degree = std::accumulate(
        adjacency_list.begin(), adjacency_list.end(), size_t(0),
        [](size_t sum, const auto& nbrs) { return sum + nbrs.size(); }
    );
    size_t n_edges = total_degree / 2;  // integer division

    // 2) Collect each undirected edge exactly once
    std::vector<double> lengths;
    lengths.reserve(n_edges);
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        for (auto const& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            if (j > i) {
                lengths.push_back(edge.weight);
            }
        }
    }

    // 2) Handle empty graph
    if (lengths.empty()) {
        return 0.0;
    }

    size_t n = lengths.size();
    size_t mid = n / 2;

    // 3) Partition around the middle
    std::nth_element(lengths.begin(), lengths.begin() + mid, lengths.end());
    double median = lengths[mid];

    // 4) If even count, average with max of lower half
    if (n % 2 == 0) {
        double lowerMax = *std::max_element(lengths.begin(), lengths.begin() + mid);
        median = (lowerMax + median) / 2.0;
    }

    return median;
}
