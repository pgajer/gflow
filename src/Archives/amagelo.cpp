#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
    // Undefine conflicting macros from R headers
#undef length

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <queue>                    // Foe std::queue

#include "cpp_utils.hpp"            // For debugging
#include <filesystem>
#include <experimental/filesystem>

#include "amagelo.hpp" // For amagelo_t
#include "uniform_grid_graph.hpp"
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t
#include "sort_utils.hpp"           // For sort_by_x_keep_y_and_order()
#include "graph_utils.hpp"          // For create_chain_graph()
#include "kernel_utils.hpp"         // For get_weights()

/**
 * @brief  Perform model averaged uniform grid local linear smoothing of 1d data
 *
 * @param x                        Vector of predictor variable.
 * @param y                        Vector of observed responses.
 * @param grid_size                Number of elements of the uniform grid over [min(x), max(x)].
 * @param min_bw_factor            Minimum bandwidth = min_bw_factor × (max(x) - min(x)).
 * @param max_bw_factor            Maximum bandwidth = max_bw_factor × (max(x) - min(x)).
 * @param n_bws                    Number of bandwidth values to evaluate.
 * @param log_grid                 If true, use logarithmic spacing of bandwidths.
 * @param hbhd_min_size     Minimum number of neighbors per vertex to define its lower‐bound bandwidth.
 * @param kernel_type              Enum / integer code of the kernel (Normal, Laplace, etc.).
 * @param dist_normalization_factor  Factor to scale distances before kernel evaluation.
 * @param with_bw_predictions      If true, store per‐bandwidth prediction vectors.
 * @param precision                Minimum spacing between successive bandwidths.
 * @param verbose                  If true, emit Rprintf diagnostics during execution.
 *
 * @return A amagelo_t containing:
 *   - predictions: final smoothed values at the optimal bandwidth,
 *   - bw_predictions: per‐bandwidth predictions (if requested),
 *   - bw_errors: mean LOOCV errors for each bandwidth,
 *   - opt_bw_idx: index of the chosen bandwidth,
 *   - min_bw,
 *   - max_bw,
 *   - bws
 */
amagelo_t amagelo(
    const std::vector<double>& x,
    const std::vector<double>& y,
    size_t grid_size,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool use_global_bw_grid,
    bool with_bw_predictions,
    bool log_grid,
    size_t domain_min_size,
    size_t kernel_type,
    double dist_normalization_factor,
    size_t n_cleveland_iterations,
    double blending_coef,
    bool use_linear_blending,
    double precision,
    bool verbose
    ) {

    amagelo_t result;

    size_t n_original_vertices = x.size();
    if (domain_min_size > n_original_vertices) {
        REPORT_ERROR("domain_min_size (%zu) exceeds total points (%zu)\n",
                     domain_min_size, n_original_vertices);
    }

    sort_by_x_keep_y_and_order(
        x,
        y,
        result.x_sorted,
        result.y_sorted,
        result.order
        );

    // Create a uniform grid graph from result.x_sorted
    auto [adj_list, weight_list] = create_chain_graph(result.x_sorted);

    double min_diff = INFINITY;
    double d;
    for (size_t i = 1; i < n_original_vertices; i++) {
        if ((d = x[i] - x[i-1]) < min_diff) {
            min_diff = d;
        }
    }
    double snap_tolerance = 1e-2 * min_diff;
    size_t start_vertex = 0;
    uniform_grid_graph_t x_graph = create_uniform_grid_graph(
        adj_list,
        weight_list,
        grid_size,
        start_vertex,
        snap_tolerance);

    x_graph.compute_graph_diameter();

    std::unordered_map<size_t, double> grid_coords_map =
        x_graph.compute_shortest_path_distances(0, x_graph.grid_vertices);

    initialize_kernel(kernel_type, 1.0);

    double graph_diameter = x_graph.graph_diameter; // = result.x_sorted[n_original_vertices - 1] - result.x_sorted[0];

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    result.min_bw = min_bw;
    result.max_bw = max_bw;

    std::vector<double> candidate_bws;
    if (use_global_bw_grid) {
        candidate_bws = get_candidate_bws(
            min_bw,
            max_bw,
            n_bws,
            log_grid,
            precision);
        result.bws = candidate_bws;
    }

    // Initialize result structure
    result.predictions.resize(n_original_vertices);
    result.bw_predictions.resize(
        n_bws,
        std::vector<double>(n_original_vertices, 0.0)
        );
    result.bw_errors.resize(n_bws, 0.0);

#if 0
    result.grid_predictions.resize(grid_size);
    result.bw_grid_predictions.resize(
        n_bws,
        std::vector<double>(grid_size, 0.0)
        );
#endif

    if (verbose) {
        Rprintf("Starting amagelo() with Buffer Zone Method\n");
        Rprintf("Number of original vertices: %zu\n", n_original_vertices);
        Rprintf("Number of grid vertices: %zu\n", x_graph.grid_vertices.size());
        Rprintf("min_bw_factor: %f\n", min_bw_factor);
        Rprintf("max_bw_factor: %f\n", max_bw_factor);
        Rprintf("min_bw: %f\n", min_bw);
        Rprintf("max_bw: %f\n", max_bw);
        Rprintf("graph_diameter: %f\n", graph_diameter);
    }

    // weight/prediction/error/mean_error/be struct needed for mode averaging and local scale estimation; we use it to record weight/prediction/error/bw of the given vertex in each model where the vertext is in the support of the model
    struct wpeme_t {
        double weight;
        double prediction;
        double error;
        double mean_error;

        // Constructor needed for emplace_back(w,p,e,me)
        wpeme_t(double w, double p, double e, double me)
            : weight(w), prediction(p), error(e), mean_error(me) {}
    };

    std::vector<std::vector<std::vector<wpeme_t>>> bw_vertex_wpeme(n_bws); // wpe[bw_idx][i] stores a vector of {weight, prediction, error, mean_error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        bw_vertex_wpeme[bw_idx].resize(n_original_vertices);
    }

// disabling grid predictions
#if 0
    struct grid_wpme_t {
        double weight;
        double prediction;
        double mean_error;

        // Constructor needed for emplace_back(w,p,me)
        grid_wpme_t(double w, double p, double me)
            : weight(w), prediction(p), mean_error(me) {}
    };


    std::vector<
        std::unordered_map<size_t,std::vector<grid_wpme_t>>
        > bw_grid_vertex_wpme(n_bws);

    // Reserving all buckets so we never rehash mid‑flight
    for (auto &mp : bw_grid_vertex_wpme) {
        mp.reserve( x_graph.grid_vertices.size() );
    }

    // Seeding every key with an empty vector
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        auto &mp = bw_grid_vertex_wpme[bw_idx];
        for (auto gv : x_graph.grid_vertices) {
            mp.try_emplace(gv, std::vector<grid_wpme_t>{});
        }
    }
#endif

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    for (const auto& grid_vertex : x_graph.grid_vertices) {

        double grid_vertex_min_bw = x_graph.find_grid_minimum_radius_for_domain_min_size(
            grid_vertex,
            min_bw,
            max_bw,
            domain_min_size,
            precision
            );

        std::vector<double> grid_vertex_candidate_bws;
        if (!use_global_bw_grid) {
            if (grid_vertex_min_bw < min_bw) grid_vertex_min_bw = min_bw;

            grid_vertex_candidate_bws = get_candidate_bws(
                grid_vertex_min_bw, max_bw, n_bws, log_grid, precision);

        } else {

            grid_vertex_candidate_bws = candidate_bws;
            for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
                if (grid_vertex_candidate_bws[bw_idx] < grid_vertex_min_bw) {
                    grid_vertex_candidate_bws[bw_idx] = grid_vertex_min_bw;
                }
            }
        }

        auto [original_vertex_map, grid_vertex_map] = x_graph.find_original_vertices_within_radius(grid_vertex, max_bw);

        for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

            double current_bw = grid_vertex_candidate_bws[bw_idx];

            // Rprintf("\n------\nbw_idx: %zu  current_bw: %.3f\n", bw_idx, current_bw);

            std::vector<double> local_y;
            std::vector<double> local_x;
            std::vector<double> local_d;
            std::vector<size_t> local_vertices;
            for (const auto& [v, dist] : original_vertex_map) {
                if (dist <= current_bw) {
                    local_y.push_back(result.y_sorted[v]);
                    local_x.push_back(result.x_sorted[v]);
                    local_d.push_back(dist);
                    local_vertices.push_back(v);
                }
            }

            size_t nbhd_size = local_x.size();
            if (nbhd_size < domain_min_size) {
                REPORT_ERROR("local_x.size() < domain_min_size\n");
            }

            std::vector<double> local_w = get_weights(local_d, dist_normalization_factor);

#if 0
            // debugging
            print_vect(local_x, "local_x");
            print_vect(local_y, "local_y");
            print_vect(local_w, "local_w");
            print_vect(local_vertices, "local_vertices");
#endif

            double tolerance = 1e-6;
            double robust_scale = 6.0;
            ulm_t model = cleveland_ulm(
                local_x.data(),
                local_y.data(),
                local_w,
                y_binary,
                tolerance,
                n_cleveland_iterations,
                robust_scale);

            double model_mean_error = std::accumulate(model.errors.begin(), model.errors.end(), 0.0) / model.errors.size();

            for (size_t i = 0; i < nbhd_size; ++i) {
                bw_vertex_wpeme[bw_idx][ local_vertices[i] ].emplace_back(
                    local_w[i],
                    model.predictions[i],
                    model.errors[i],
                    model_mean_error
                    );
            }


// disabling grid predictions !!!
#if 0
            // Predict model values at local grid vertices
            std::vector<size_t> local_grid_vertices;
            std::vector<double> local_grid_d;
            for (const auto& [gv, dist] : grid_vertex_map) {
                if (dist <= current_bw) {
                    local_grid_vertices.push_back(gv);
                    local_grid_d.push_back(dist);
                }
            }

            // Use the actual x‐coordinate of each grid vertex to predict model values
            std::vector<double> grid_coords;
            grid_coords.reserve(local_grid_vertices.size());
            for (auto gv : local_grid_vertices)
                grid_coords.push_back(result.x_sorted[0] + grid_coords_map[gv]); // grid_coords_map[gv] is the distance from the first vertex (index 0) which has the x_sorted value of x_sorted[0]. Since for model building we are using x_sorted coordinates as predictors, we need to use the same coordinates for prediction

            auto local_grid_predictions = model.predict(grid_coords);

            // Generate local grid weights for model averaging
            std::vector<double> local_grid_w = get_weights(local_grid_d, dist_normalization_factor);

            if (local_grid_vertices.size() != local_grid_w.size()) {
                REPORT_ERROR("ERROR: local_grid_vertices.size(): %zu\n"
                             "local_grid_w.size(): %zu\n"
                             "They suppose to be the same\n",
                    local_grid_vertices.size(), local_grid_w.size());
            }

            if (local_grid_vertices.size() != local_grid_predictions.size()) {
                REPORT_ERROR("ERROR: local_grid_vertices.size(): %zu\n"
                             "local_grid_predictions.size(): %zu\n"
                             "They suppose to be the same\n",
                             local_grid_vertices.size(), local_grid_predictions.size());
            }

            for (size_t i = 0; i < local_grid_vertices.size(); ++i) {
                size_t gv = local_grid_vertices[i];
                double w  = local_grid_w[i];
                double p  = local_grid_predictions[i];
                try {
                    // 1) bounds‑checked access to the bw_idx slot
                    auto &grid_map = bw_grid_vertex_wpme.at(bw_idx);
                    // 2) bounds‑checked access to the gv key
                    auto &wp_vec   = grid_map.at(gv);
                    wp_vec.emplace_back(w, p, model_mean_error);
                }
                catch (const std::out_of_range &e) {
                    // You hit this if bw_idx or gv was invalid
                    Rprintf("ERROR in amagelo(): %s\n", e.what());
                    REPORT_ERROR("amagelo internal error: invalid map access (bw_idx=%zu, grid_vertex=%zu)\n",
                                 bw_idx, gv);
                }
            }
#endif
// disabling grid predictions !!!

        }
    }

    //
    // Model Averaging over the original vertices
    //
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

        for (size_t vertex = 0; vertex < n_original_vertices; ++vertex) {

            auto vertex_wpes = bw_vertex_wpeme[bw_idx][vertex];
            if (vertex_wpes.size() == 0) {
                REPORT_ERROR("ERROR: No prediction data found for vertex: %zu\n", vertex);
            }

            double prediction_sum = 0.0;
            double weight_sum     = 0.0;
            double error_sum      = 0.0;
            double effective_weight;

            for (const auto& x : vertex_wpes) {

                if (blending_coef == 0.0) {
                    // Pure position weight
                    effective_weight = x.weight;
                }
                else if (blending_coef == 1.0) {
                    // Full mean error influence
                    effective_weight = x.mean_error * x.weight;
                }
                else {
                    // Smooth interpolation between the two approaches
                    if (use_linear_blending) {
                        // Method 1: Linear interpolation of weights
                        effective_weight = (1.0 - blending_coef) * x.weight + blending_coef * (x.mean_error * x.weight);
                    } else {
                        // Method 2: Power-based scaling
                        effective_weight = x.weight * pow(x.mean_error, blending_coef);
                    }
                }

                prediction_sum += effective_weight * x.prediction;
                error_sum      += effective_weight * x.mean_error;
                weight_sum     += effective_weight;
            }

            if (weight_sum > 0) {
                if (weight_sum > 1e-10) {
                    result.bw_predictions[bw_idx][vertex] = prediction_sum / weight_sum;
                    result.bw_errors[bw_idx] += error_sum / weight_sum;
                } else {
                    REPORT_ERROR("ERROR: Very small weight sum encountered for vertex %zu. Setting predictions[i] to NaN\n", vertex);
                    //result.bw_predictions[bw_idx][vertex] = std::numeric_limits<double>::quiet_NaN();
                }
            } else {
                REPORT_ERROR("Weight sum = 0 for vertex %zu in predictions\n", vertex);
            }
        }
    }


#if 0
// disabling grid predictions !!!
    //
    // Model Averaging over grid vertices
    //
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

        for (const auto& grid_vertex : x_graph.grid_vertices) {

            if (!bw_grid_vertex_wpme[bw_idx].contains(grid_vertex)) {
                Rprintf("ERROR: grid_vertex: %zu  not found in bw_grid_vertex_wpme[bw_idx] for bw_idx: %zu\n",
                        grid_vertex, bw_idx);

                print_uset(x_graph.grid_vertices, "x_graph.grid_vertices");

                Rprintf("bw_grid_vertex_wpme[bw_idx].size(): %zu\n", bw_grid_vertex_wpme[bw_idx].size());

                REPORT_ERROR("ERROR: grid_vertex: %zu  not found in bw_grid_vertex_wpme[bw_idx] for bw_idx: %zu\n",
                    grid_vertex, bw_idx);
            }

            auto grid_vertex_wps = bw_grid_vertex_wpme[bw_idx][grid_vertex];
            if (grid_vertex_wps.size() == 0) {
                REPORT_ERROR("ERROR: No prediction data found for grid_vertex: %zu\n", grid_vertex);
            }

            double prediction_sum = 0.0;
            double weight_sum     = 0.0;
            double effective_weight;

            for (const auto& x : grid_vertex_wps) {

                if (blending_coef == 0.0) {
                    // Pure position weight
                    effective_weight = x.weight;
                }
                else if (blending_coef == 1.0) {
                    // Full mean error influence
                    effective_weight = x.mean_error * x.weight;
                }
                else {
                    // Smooth interpolation between the two approaches
                    if (use_linear_blending) {
                        // Method 1: Linear interpolation of weights
                        effective_weight = (1.0 - blending_coef) * x.weight + blending_coef * (x.mean_error * x.weight);
                    } else {
                        // Method 2: Power-based scaling
                        effective_weight = x.weight * pow(x.mean_error, blending_coef);
                    }
                }

                prediction_sum += effective_weight * x.prediction;
                weight_sum     += effective_weight;
            }

            if (weight_sum > 0) {
                if (weight_sum > 1e-10) {
                    result.bw_grid_predictions[bw_idx][grid_vertex] = prediction_sum / weight_sum;
                } else {
                    REPORT_ERROR("ERROR: Very small weight sum encountered for grid_vertex %zu. Setting predictions[i] to NaN\n", grid_vertex);
                    // predictions[grid_vertex] = std::numeric_limits<double>::quiet_NaN();
                }
            } else {
                REPORT_ERROR("Weight sum = 0 for grid_vertex %zu in predictions\n", grid_vertex);
            }
        }
    }
#endif

    //
    // Find optimal bw
    //
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        result.bw_errors[bw_idx] /= n_original_vertices;
    }

    result.opt_bw_idx = std::min_element(result.bw_errors.begin(),
                                         result.bw_errors.end()) - result.bw_errors.begin();

    if (!with_bw_predictions) {
        result.predictions = std::move(result.bw_predictions[result.opt_bw_idx]);

        // Clear all other bandwidth predictions to free memory
        for (size_t i = 0; i < n_bws; ++i) {
            if (i != result.opt_bw_idx) {
                std::vector<double>().swap(result.bw_predictions[i]);  // Swap idiom for capacity release
            }
        }

#if 0
// disabling grid predictions !!!
        result.grid_predictions = std::move(result.bw_grid_predictions[result.opt_bw_idx]);

        for (size_t i = 0; i < n_bws; ++i) {
            if (i != result.opt_bw_idx) {
                std::vector<double>().swap(result.bw_grid_predictions[i]);  // Swap idiom for capacity release
            }
        }
#endif

    } else {
        result.predictions = result.bw_predictions[result.opt_bw_idx];
        //result.grid_predictions = result.bw_grid_predictions[result.opt_bw_idx];

        if (y_binary) {
            for (size_t i = 0; i < result.predictions.size(); ++i) {
                result.predictions[i] = std::clamp(result.predictions[i], 0.0, 1.0);
            }
            //   for (size_t i = 0; i < result.grid_predictions.size(); ++i) {
            //         result.grid_predictions[i] = std::clamp(result.grid_predictions[i], 0.0, 1.0);
            //     }
        }
    }

    return result;
}

/**
 * @brief Finds original vertices within a specified radius of a reference vertex
 *
 * @param vertex The reference vertex from which distances are measured
 * @param radius Maximum distance threshold for inclusion
 * @return std::unordered_map<size_t, double> Map of vertices to their distances from the reference vertex
 *
 * @details This function computes the shortest path distance from the reference vertex
 * to all other vertices in the graph, and returns those vertices whose distance
 * is less than or equal to the specified radius along with their distances.
 * The reference vertex is always included in the result with distance 0.
 */
std::pair<std::unordered_map<size_t, double>,
          std::unordered_map<size_t, double>>
uniform_grid_graph_t::find_original_vertices_within_radius(
    size_t grid_vertex,
    double radius) const {

    // Rprintf("Entering find_original_vertices_within_radius(): grid_vertex: %zu\n", grid_vertex);

    std::unordered_map<size_t, double> original_vertex_map; // maps original vertices to their distances from the reference grid_vertex
    std::unordered_map<size_t, double> grid_vertex_map;     // maps grid vertices to their distances from the reference grid_vertex

    // Initialize priority queue for Dijkstra's algorithm
    // Using min-heap with pairs of (distance, grid_vertex)
    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<>> pq;

    // Track distances
    std::vector<double> distances(adjacency_list.size(), INFINITY);
    distances[grid_vertex] = 0.0;

    // Include the reference grid_vertex with distance 0 if it is an original vertex
    if (is_original_vertex(grid_vertex)) {
        original_vertex_map[grid_vertex] = 0.0;
    } else {
        grid_vertex_map[grid_vertex] = 0.0;
    }

    // Initialize with reference grid_vertex
    pq.push({0.0, grid_vertex});

    // Run Dijkstra's algorithm
    while (!pq.empty()) {
        auto [current_dist, current_vertex] = pq.top();
        pq.pop();

        // If we've exceeded the radius, we can terminate early
        if (current_dist > radius) {
            break;
        }

        // Skip if we've already found a shorter path
        if (current_dist > distances[current_vertex]) {
            continue;
        }

        // Add this vertex to the original_vertex_map with its distance
        if (is_original_vertex(current_vertex)) {
            original_vertex_map[current_vertex] = current_dist;
        }

        if (grid_vertices.contains(current_vertex)) {
            grid_vertex_map[current_vertex] = current_dist;
        }

        // Explore neighbors
        for (const auto& edge : adjacency_list[current_vertex]) {
            size_t neighbor = edge.vertex;
            double new_dist = current_dist + edge.weight;

            // Update distance if we found a shorter path
            if (new_dist < distances[neighbor]) {
                distances[neighbor] = new_dist;
                pq.push({new_dist, neighbor});
            }
        }
    }

    return std::make_pair(original_vertex_map, grid_vertex_map);
}


/**
 * @brief Finds the minimum radius required to include at least domain_min_size original vertices
 *
 * @param vertex The center vertex from which to measure distances
 * @param lower_bound Initial lower bound for binary search
 * @param upper_bound Initial upper bound for binary search
 * @param domain_min_size Minimum number of vertices required in the neighborhood
 * @param precision Precision threshold for terminating binary search
 *
 * @return double The minimum radius that ensures at least domain_min_size vertices
 *                are within the neighborhood of the vertex
 *
 * @details This function performs a binary search to find the smallest radius value
 * for which the vertex neighborhood (set of all vertices within the radius) contains
 * at least domain_min_size elements. This is required to ensure sufficient data
 * points for linear model fitting.
 */
double uniform_grid_graph_t::find_grid_minimum_radius_for_domain_min_size(
    size_t grid_vertex,
    double lower_bound,
    double upper_bound,
    size_t domain_min_size,
    double precision) const {

#define DEBUG__find_grid_minimum_radius_for_domain_min_size 0

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
    Rprintf("Entering find_grid_minimum_radius_for_domain_min_size()\ngrid_vertex: %zu\n"
            "lower_bound: %.3f\n"
            "upper_bound: %.3f\n"
            "domain_min_size: %zu\n",
            grid_vertex, lower_bound, upper_bound, domain_min_size);
#endif

    {
        // Check if the lower bound already satisfies the condition
        auto [orig_low, grid_low] = find_original_vertices_within_radius(grid_vertex, lower_bound);

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
        Rprintf("orig_low.size(): %zu\n", orig_low.size());
#endif

        if (orig_low.size() >= domain_min_size)
            return lower_bound;
    }

    {
        // Check if the upper bound fails to satisfy the condition
        auto [orig_high, grid_high] = find_original_vertices_within_radius(grid_vertex, upper_bound);

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
        Rprintf("orig_high.size(): %zu\n", orig_high.size());
#endif

        if (orig_high.size() < domain_min_size) {
            Rprintf("\n---------------------\nERROR\ngrid_vertex: %zu\n"
                    "Not enough vertices (found: %zu, needed: %zu) within upper_bound: %.4f\n",
                    grid_vertex + 1, orig_high.size(), domain_min_size, upper_bound);
            REPORT_ERROR("ERROR: Insufficient vertices in maximum radius\n");
        }
    }

    // Binary search for the minimum radius
#if DEBUG__find_grid_minimum_radius_for_domain_min_size
    Rprintf("(upper_bound - lower_bound): %.3f\n"
            "precision: %.3f\n"
            "graph_diameter: %.3f\n"
            "precision * graph_diameter: %.3f\n",
            upper_bound - lower_bound,
            precision,
            graph_diameter,
            precision * graph_diameter
        );
    size_t debug_counter = 0;
#endif

    double thld = precision * graph_diameter;
    while ((upper_bound - lower_bound) > thld) {

        double mid = 0.5 * (lower_bound + upper_bound);

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
        Rprintf("debug_counter: %zu\t"
                "lower_bound: %.3f\t"
                "upper_bound: %.3f\t"
                "(upper_bound - lower_bound): %.3f\t"
                "thld: %.3f\t"
                "mid: %.3f\n",
                debug_counter,
                lower_bound,
                upper_bound,
                upper_bound - lower_bound,
                thld,
                mid);
#endif

        auto [orig_mid, grid_mid] = find_original_vertices_within_radius(grid_vertex, mid);

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
        Rprintf("orig_mid.size(): %zu\n", orig_mid.size());
#endif

        if (orig_mid.size() == domain_min_size) {
            upper_bound = mid;
            break;
        } else if (orig_mid.size() > domain_min_size) {
            // If condition is satisfied, try a smaller radius
            upper_bound = mid;
            //Rprintf("Setting upper_bound to mid\n");

        } else {
            // If condition is not satisfied, try a larger radius
            lower_bound = mid;
            //Rprintf("Setting lower_bound to mid\n");
        }

#if DEBUG__find_grid_minimum_radius_for_domain_min_size
        if (debug_counter == 20) {
            error("Reached debug_counter: %zu\nAborting!\n", debug_counter);
        }

        debug_counter++;
#endif
    }

    // Return the upper bound as the minimum radius that satisfies the condition
    return upper_bound;
}
