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
    result.predictions.resize(n_vertices);        // kernel smoothing predictions for each vertex
    result.mean_bw_errors.resize(n_bws, std::numeric_limits<double>::quiet_NaN()); // mean_bw_errors[bw_idx] is the mean prediction error for the bw_idx-th bandwidth
    result.vertex_min_bws.resize(n_vertices);     // vertex_min_bws[i] is the the i-th vertex's minimum bandwidth
    result.vertex_opt_bws.resize(n_vertices);     // optimal bandwidth value for each vertex
    result.bw_lextr_count.resize(n_bws);

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

        #if 0
        Rprintf("n_bws: %zu\n", n_bws);
        print_vect(bw_grid, "bw_grid");

        Rprintf("folds:\n");
        for (size_t fold = 0; fold < n_folds; ++fold) {
            Rprintf("%zu: ", fold);
            print_vect(folds[fold]);
        }
        #endif
    }

    std::vector<size_t> neighbors(n_vertices);
    std::vector<double> distances(n_vertices);
    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {

#if DEBUG__graph_kernel_smoother
        Rprintf("\n\nvertex: %zu\n", vertex);
#endif

        // 1. compute neighbor distance map of the given vertex for the maximum bw
        std::unordered_map<size_t, double> ngbr_dist_map = find_vertices_within_radius(vertex, max_bw);

#if DEBUG__graph_kernel_smoother
        print_umap(ngbr_dist_map, "ngbr_dist_map");
#endif

        // 2. Use ngbr_dist_map to find vertex specific min_bw
        auto [sorted_vertices, vertex_min_bw] = get_sorted_vertices_and_min_radius(
            ngbr_dist_map,
            vertex_hbhd_min_size
            );

#if DEBUG__graph_kernel_smoother
        Rprintf("vertex_min_bw: %.4f\n", vertex_min_bw);
#endif

        if (vertex_min_bw < min_bw) {
            vertex_min_bw = min_bw;
        }

#if DEBUG__graph_kernel_smoother
        Rprintf("vertex_min_bw: %.4f\n", vertex_min_bw);
        print_vect_pair(sorted_vertices, "sorted_vertices");
#endif
     
        // 3. Generate candidate bandwidths
        std::vector<double> candidate_bws = get_candidate_bws(
            vertex_min_bw,
            max_bw,
            n_bws,
            log_grid,
            precision
            );

#if DEBUG__graph_kernel_smoother
        print_vect(candidate_bws, "candidate_bws");
#endif

        // 4. For each bw generate kernel-weighted mean prediction at the given vertex
        for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {

            double current_bw = candidate_bws[bw_idx];

#if DEBUG__graph_kernel_smoother
            Rprintf("bw_idx: %zu current_bw: %.3f\n", bw_idx, current_bw);
#endif

            size_t nbhd_size = 0;
            for (const auto& [v, dist] : sorted_vertices) {
                if (dist >= current_bw)
                    break;
                neighbors[nbhd_size] = v;
                distances[nbhd_size] = dist;
                nbhd_size++;
            }

#if DEBUG__graph_kernel_smoother
            Rprintf("nbhd_size: %zu\n", nbhd_size);
            {
                // Make sure nbhd_size ≤ neighbors.size()
                size_t k = std::min(nbhd_size, neighbors.size());

                // This allocates a new vector and copies elements [0, k)
                std::vector<size_t> neighbors_trimmed(
                    neighbors.begin(),
                    neighbors.begin() + k
                    );

                std::vector<double> distances_trimmed(
                    distances.begin(),
                    distances.begin() + k
                    );

                print_vect(neighbors_trimmed, "neighbors");
                print_vect(distances_trimmed, "distances");
            }
#endif

            double max_dist = distances[nbhd_size - 1];
            if (max_dist == 0) max_dist = 1.0;

#if DEBUG__graph_kernel_smoother
            Rprintf("max_dist: %.4f\n", max_dist);
#endif

            // Calculate prediction based on training vertices
            std::vector<double> weights(nbhd_size);

            if (use_uniform_weights) {
                // Use uniform weights
                std::fill(weights.begin(), weights.end(), 1.0);
            } else {
                max_dist *= dist_normalization_factor;

                // Normalize by the scaled maximum
                for (size_t i = 0; i < distances.size(); ++i) {
                    distances[i] /= max_dist;
                }

                // Use the kernel function to calculate weights
                kernel_fn(distances.data(), static_cast<int>(nbhd_size), weights.data());
            }

#if DEBUG__graph_kernel_smoother
            print_vect(weights, "weights");
#endif

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

#if DEBUG__graph_kernel_smoother
            Rprintf("weighted_sum: %.4f\n", weighted_sum);
            Rprintf("weight_sum: %.4f\n", weight_sum);
            Rprintf("result.bw_predictions[bw_idx][vertex]: %.4f y[vertex]: %.4f\n", result.bw_predictions[bw_idx][vertex], y[vertex]);
            error("DEBUGGING");
#endif

        } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)

    } // END OF for (size_t vertex = 0; vertex < n_vertices; ++vertex)

    // Compute the number of local extrema
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        result.bw_lextr_count[bw_idx] = get_lextr_count(result.bw_predictions[bw_idx]);
    }

    // Perform CV - TO BE DONE LATER

    // Find optimal bandwidth - TO BE DONE LATER
    // Compute the mean of each bw prediction errors
    // for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {
    //     double mean_error = std::accumulate(result.bw_errors[bw_idx].begin(),
    //                                         result.bw_errors[bw_idx].end(), 0.0) / n_vertices;
    //     result.mean_bw_errors = mean_error;
    // }
    // auto min_it = std::min_element(result.mean_bw_errors.begin(), result.mean_bw_errors.end());
    // result.opt_bw_idx = min_it - result.mean_bw_errors.begin();
    // result.opt_bw = bw_grid[result.opt_bw_idx];


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
        min_radius = sorted_vertices.empty() ? 0.0 : sorted_vertices.back().second;
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
