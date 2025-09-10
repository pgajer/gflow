#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

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
 * @brief  Perform kernel‐based smoothing on graph‐structured data with CV‐based bandwidth selection.
 *
 * This method computes at each vertex an estimate of E[Y|G] by applying a
 * Nadaraya–Watson style smoother on graph‐distances.  It supports either
 * leave‑one‑out or stratified K‑fold cross‐validation with a buffer‑zone
 * exclusion to avoid local leakage.  After evaluating each candidate
 * bandwidth on the test folds (mean absolute Rf_error), it selects the
 * optimal bandwidth and then recomputes final predictions on all vertices.
 *
 * @param y                        Vector of observed responses of length |V|.
 * @param min_bw_factor            Minimum bandwidth = min_bw_factor × graph diameter.
 * @param max_bw_factor            Maximum bandwidth = max_bw_factor × graph diameter.
 * @param n_bws                    Number of bandwidth values to evaluate.
 * @param log_grid                 If true, use logarithmic spacing of bandwidths.
 * @param vertex_hbhd_min_size     Minimum number of neighbors per vertex to define its lower‐bound bandwidth.
 * @param kernel_type              Enum / integer code of the kernel (Normal, Laplace, etc.).
 * @param dist_normalization_factor  Factor to scale distances before kernel evaluation.
 * @param use_uniform_weights      If true, all kernel weights = 1 (no distance weighting).
 * @param buffer_hops              Number of graph‐hop layers to exclude around each test vertex.
 * @param auto_buffer_hops         If true, override buffer_hops by an automatic routine.
 * @param n_folds                  Number of CV folds (if ≥|V|, performs leave‐one‐out).
 * @param with_bw_predictions      If true, store per‐bandwidth prediction vectors.
 * @param precision                Minimum spacing between successive bandwidths.
 * @param verbose                  If true, emit Rprintf diagnostics during execution.
 *
 * @return A graph_kernel_smoother_t containing:
 *   - predictions: final smoothed values at the optimal bandwidth,
 *   - bw_predictions: per‐bandwidth predictions (if requested),
 *   - bw_mean_abs_errors: CV‐mean absolute Rf_error for each bandwidth,
 *   - vertex_min_bws: per‐vertex lower‐bound bandwidths,
 *   - opt_bw_idx: index of the chosen bandwidth,
 *   - buffer_hops_used: the buffer_hops actually applied.
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

    compute_graph_diameter();

    size_t n_vertices = adjacency_list.size();

    std::vector<std::vector<size_t>> folds;
    if (n_folds >= n_vertices) {
        n_folds = n_vertices;
        folds.resize(n_folds);
        for (size_t v = 0; v < n_vertices; ++v)
            folds[v] = {v};
    } else {
        folds = create_spatially_stratified_folds(n_folds);
    }

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    // Initialize result structure
    graph_kernel_smoother_t result;
    result.predictions.resize(n_vertices);          // kernel smoothing predictions for each vertex
    result.bw_mean_abs_errors.resize(n_bws, 0.0);   // mean_bw_abs_errors[bw_idx] is the mean prediction absolute Rf_error for the bw_idx-th bandwidth
    result.vertex_min_bws.resize(n_vertices);       // vertex_min_bws[i] is the the i-th vertex's minimum bandwidth

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
        // Rf_error("DEBUGGING\n");
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

                if (weight_sum <= 0) { // this can happen only when all neighbors of the vertex are in the buffer zone; in this case we skip the vertex
                    continue;
                }

                double prediction = weighted_sum / weight_sum;
                result.bw_mean_abs_errors[bw_idx] += std::abs(prediction - y[vertex]);

            } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)
        } // END OF for (const auto& vertex : training_set)
    } // END OF for (size_t fold = 0; fold < folds.size(); ++fold)

    // Assuming the total number of test vertices was n_vertices
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        result.bw_mean_abs_errors[bw_idx] /= n_vertices;
    }

    // Find optimal bandwidth using result.bw_mean_sq_errors
    auto min_it = std::min_element(result.bw_mean_abs_errors.begin(), result.bw_mean_abs_errors.end());
    result.opt_bw_idx = min_it - result.bw_mean_abs_errors.begin();

    // -----------------------------------------------------------------------------
    // All data prediction section
    // -----------------------------------------------------------------------------
    if (with_bw_predictions) {
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

                // Distance normalization
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

        result.predictions = result.bw_predictions[result.opt_bw_idx];

    } else {

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
            {
                size_t bw_idx = result.opt_bw_idx;
                double current_bw = candidate_bws[bw_idx];

                size_t nbhd_size = 0;
                for (const auto& [v, dist] : sorted_vertices) {
                    if (dist >= current_bw)
                        break;
                    neighbors[nbhd_size] = v;
                    distances[nbhd_size] = dist;
                    nbhd_size++;
                }

                // Distance normalization
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

                result.predictions[vertex] = weighted_sum / weight_sum;

            } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)
        } // END OF for (size_t vertex = 0; vertex < n_vertices; ++vertex)
    }

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
 * @throws Reports an Rf_error if the number of filtered vertices is insufficient to meet vertex_hbhd_min_size
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
