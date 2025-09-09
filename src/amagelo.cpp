#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
    // Undefine conflicting macros from R headers
#undef length
#undef eval

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::sort and std::for_each
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
 * @brief Perform Adaptive Model Averaged Uniform Grid LOcal linear smoothing (AMAGELO) of 1D data
 *
 * AMAGELO is a nonparametric smoothing technique that uses:
 * - A uniform grid over the data domain
 * - Local linear models at each grid point
 * - Model averaging with adaptive weighting
 * - Bandwidth optimization using cross-validation
 * - Detection of local extrema with significance measures
 * - Optional triplet harmonic smoothing for small wiggles
 *
 * @param x                        Vector of predictor variable
 * @param y                        Vector of observed responses
 * @param grid_size                Number of elements in the uniform grid over [min(x), max(x)]
 * @param min_bw_factor            Minimum bandwidth = min_bw_factor × (max(x) - min(x))
 * @param max_bw_factor            Maximum bandwidth = max_bw_factor × (max(x) - min(x))
 * @param n_bws                    Number of bandwidth values to evaluate
 * @param use_global_bw_grid       If true, use the same bandwidth grid for all vertices
 * @param with_bw_predictions      If true, store per-bandwidth prediction vectors
 * @param log_grid                 If true, use logarithmic spacing of bandwidths
 * @param domain_min_size          Minimum number of points required in each local model's domain
 * @param kernel_type              Integer code for kernel function (e.g., Normal, Laplace)
 * @param dist_normalization_factor Factor to scale distances before kernel evaluation
 * @param n_cleveland_iterations   Number of robustness iterations in local linear fitting
 * @param blending_coef            Control parameter for model averaging (0=pure position weights, 1=full error influence)
 * @param use_linear_blending      If true, use linear blending instead of power-based scaling
 * @param precision                Minimum spacing between successive bandwidths in optimization
 * @param small_depth_threshold    Threshold for identifying small wiggles in triplet smoothing
 * @param depth_similarity_tol     Tolerance for depth similarity in triplet smoothing
 * @param verbose                  If true, emit diagnostic messages during execution
 *
 * @return An amagelo_t structure containing:
 *   - x_sorted, y_sorted: Data sorted by x values
 *   - order: Original indices of sorted data
 *   - grid_coords: x-coordinates of grid points
 *   - predictions: Final smoothed values at optimal bandwidth
 *   - bw_predictions: Per-bandwidth predictions (if requested)
 *   - grid_predictions: Predictions at grid points
 *   - harmonic_predictions: Predictions after triplet harmonic smoothing
 *   - local_extrema: Information about detected local extrema
 *   - monotonic_interval_proportions: Relative lengths of monotonicity intervals
 *   - change_scaled_monotonicity_index: Weighted signed average of directional changes, quantifying monotonicity strength and directionality; values close to \(+1\) or \(-1\) indicate strong global monotonic trends.
 *   - bw_errors: Mean cross-validation errors for each bandwidth
 *   - opt_bw_idx: Index of optimal bandwidth
 *   - min_bw, max_bw: Bandwidth range
 *   - bws: Vector of evaluated bandwidths
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
    double small_depth_threshold,
    double depth_similarity_tol,
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

    // Creating vector container of x_graph.grid_vertices
    std::vector<size_t> grid_vertices(x_graph.grid_vertices.begin(), x_graph.grid_vertices.end());
    std::sort(
        grid_vertices.begin(),
        grid_vertices.end(),
        [&](size_t lhs, size_t rhs) {
            return grid_coords_map[lhs] < grid_coords_map[rhs];
        }
        );

    std::map<size_t, size_t> grid_idx_map;
    for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {
        grid_idx_map[ grid_vertices[grid_idx] ] = grid_idx;
    }

    result.grid_coords.resize(grid_size);
    for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {
        result.grid_coords[grid_idx] = result.x_sorted[0] + grid_coords_map[ grid_vertices[grid_idx] ];
    }


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
    result.grid_predictions.resize(grid_size);
    result.bw_predictions.resize(
        n_bws,
        std::vector<double>(n_original_vertices, 0.0)
        );
    result.bw_errors.resize(n_bws, 0.0);

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

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    std::vector<std::vector<ulm_t>> bw_grid_models(n_bws);  // bw_grid_models[bw_idx][grid_idx] ulm_t model at grid_idx at bw_idx
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        bw_grid_models[bw_idx].resize(grid_size);
    }

    std::vector<std::vector<double>> bw_grid_bws(n_bws); //  bw_grid_bws[bw_idx][grid_idx] = current_bw
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {
        bw_grid_bws[bw_idx].resize(grid_size);
    }

    std::vector<std::unordered_map<size_t, double>> grid_vertex_maps(grid_size);

    //
    // Creating Local Models
    //
    for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {
        size_t grid_vertex = grid_vertices[grid_idx];

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
        grid_vertex_maps[grid_idx] = std::move(grid_vertex_map);

        for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

            double current_bw = grid_vertex_candidate_bws[bw_idx];
            bw_grid_bws[bw_idx][grid_idx] = current_bw;

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

            bw_grid_models[bw_idx][grid_idx] = std::move(model);
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

    } else {
        result.predictions = result.bw_predictions[result.opt_bw_idx];
        if (y_binary) {
            for (size_t i = 0; i < result.predictions.size(); ++i) {
                result.predictions[i] = std::clamp(result.predictions[i], 0.0, 1.0);
            }
        }
    }

    //
    // Creating optimal bw grid_predictions
    //
    struct grid_wpme_t {
        double weight;
        double prediction;
        double mean_error;

        // Constructor needed for emplace_back(w,p,me)
        grid_wpme_t(double w, double p, double me)
            : weight(w), prediction(p), mean_error(me) {}
    };

    std::vector<std::vector<grid_wpme_t>> grid_vertex_wpme(grid_size); // grid_vertex_wpme[grid_idx] a vector of wpme_t's IMPORTANT: grid_vertex_wpme is indexed by grid_idx no grid_vertex

    for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {

        const auto& model = bw_grid_models[result.opt_bw_idx][grid_idx];
        double model_mean_error = std::accumulate(model.errors.begin(), model.errors.end(), 0.0) / model.errors.size();

        // Predict model values at local grid vertices
        double current_bw = bw_grid_bws[result.opt_bw_idx][grid_idx];
        std::vector<size_t> local_grid_vertices;
        std::vector<double> local_grid_d;
        for (const auto& [gv, dist] : grid_vertex_maps[grid_idx]) {
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

        for (size_t i = 0; i < local_grid_vertices.size(); ++i) {
            grid_vertex_wpme[ grid_idx_map[ local_grid_vertices[i] ] ].emplace_back(
                local_grid_w[i],
                local_grid_predictions[i],
                model_mean_error
                );
        }
    }

    for (size_t grid_idx = 0; grid_idx < grid_size; ++grid_idx) {

        auto gv_wpmes = grid_vertex_wpme[grid_idx];
        if (gv_wpmes.size() == 0) {
            REPORT_ERROR("ERROR: No prediction data found for grid_idx: %zu\n", grid_idx);
        }

        double prediction_sum = 0.0;
        double weight_sum     = 0.0;
        double error_sum      = 0.0;
        double effective_weight;

        for (const auto& x : gv_wpmes) {

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

                result.grid_predictions[grid_idx] = prediction_sum / weight_sum;
                if (y_binary)
                    result.grid_predictions[grid_idx] = std::clamp(result.grid_predictions[grid_idx], 0.0, 1.0);

            } else {
                REPORT_ERROR("ERROR: Very small weight sum encountered for grid_idx %zu. Setting predictions[i] to NaN\n", grid_idx);
            }
        } else {
            REPORT_ERROR("Weight sum = 0 for grid_idx %zu in predictions\n", grid_idx);
        }
    }

#if 0
    // Compute Total‐Variation Monotonicity Index (TVMI)
    {
        const auto& gp = result.grid_predictions;
        const auto& gc = result.grid_coords;
        size_t G = gp.size();
        // 1. sum of absolute grid diffs
        double sum_abs = 0.0;
        for (size_t i = 1; i < G; ++i) {
            sum_abs += std::abs(gp[i] - gp[i-1]);
        }
        // 2. grid span (b - a)
        double span = gc.back() - gc.front();
        // avoid division by zero
        const double eps = 1e-12;
        if (span <= eps || sum_abs <= eps) {
            result.tvmi = 0.0;  // no variation ⇒ define TVMI=0
        } else {
            // ∆1 = (ŷ(b) - ŷ(a)) / span
            // double delta1 = (gp.back() - gp.front()) / span;
            // δ1 = (1/span) * sum_abs  ⇒ so delta1/δ1 = (gp.back()-gp.front())/sum_abs
            // result.tvmi = delta1 / (sum_abs / span);
            // simplifies to
            result.tvmi = (gp.back() - gp.front()) / sum_abs;
        }
    }

    // Compute relative monotonic interval lengths and the Simpson index
    {
        const auto& yhat  = result.grid_predictions; // holds ŷ values on uniform grid
        const auto& xgrid = result.grid_coords;      // holds corresponding x-coordinates
        size_t n = yhat.size();
        std::vector<double> monotonic_interval_lengths;
        if (n < 2) {
            monotonic_interval_lengths.push_back(xgrid.back() - xgrid.front());
        } else {
            size_t start_idx = 0;
            // Determine initial direction: +1 for increasing, -1 for decreasing, 0 for flat
            int direction = 0;
            for (size_t i = 1; i < n; ++i) {
                double diff = yhat[i] - yhat[i-1];
                int curr = (diff > 0) ? 1 : (diff < 0 ? -1 : 0);
                if (curr == 0) continue;  // Skip flat segments
                direction = curr;
                break;
            }
            // Iterate through grid, splitting intervals when direction changes
            for (size_t i = 1; i < n; ++i) {
                double diff = yhat[i] - yhat[i-1];
                int curr = (diff > 0) ? 1 : (diff < 0 ? -1 : 0);
                if (curr == 0 || curr == direction) {
                    // Continue current monotonic interval
                    continue;
                }
                // Direction change: close previous interval
                double length = xgrid[i-1] - xgrid[start_idx];
                monotonic_interval_lengths.push_back(length);
                // Start new interval at previous point
                start_idx = i - 1;
                direction = curr;
            }
            // Close last interval
            double last_length = xgrid[n-1] - xgrid[start_idx];
            monotonic_interval_lengths.push_back(last_length);
        }

        // relative monotonic interval
        double total = std::accumulate(monotonic_interval_lengths.begin(), monotonic_interval_lengths.end(), 0.0);
        std::vector<double>& p = result.monotonic_interval_proportions;
        p.resize(monotonic_interval_lengths.size());
        for (size_t i = 0; i < monotonic_interval_lengths.size(); ++i) {
            p[i] = monotonic_interval_lengths[i] / total;
        }

        // Simpson index = ∑ p_i^2
        result.simpson_index = std::inner_product(p.begin(), p.end(), p.begin(), 0.0);
    }
#endif

    // --- Local Extrema & Depths ---
    {
        result.local_extrema.clear();  // vector<extremum_t> field in amagelo_t

        const auto& yhat = result.predictions;
        const auto& x = result.x_sorted;
        size_t n = yhat.size();

        // Step 1: Detect local extrema (including endpoints)
        for (size_t i = 0; i < n; ++i) {
            double curr = yhat[i];
            bool is_extremum = false;
            bool is_max = false;

            if (i == 0 && n > 1) {
                // Left endpoint
                if (curr > yhat[i + 1]) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < yhat[i + 1]) {
                    is_extremum = true;
                    is_max = false;
                }
            } else if (i == n - 1 && n > 1) {
                // Right endpoint
                if (curr > yhat[i - 1]) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < yhat[i - 1]) {
                    is_extremum = true;
                    is_max = false;
                }
            } else if (i > 0 && i + 1 < n) {
                // Interior points
                double prev = yhat[i - 1], next = yhat[i + 1];
                if (curr > prev && curr > next) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < prev && curr < next) {
                    is_extremum = true;
                    is_max = false;
                }
            }

            if (is_extremum) {
                result.local_extrema.emplace_back(i, x[i], curr, is_max, 0.0);
            }
        }

        // Step 2: Compute depths of local extrema
        for (auto& e : result.local_extrema) {
            size_t i = e.idx;

            if (e.is_max) {
                // Local maximum: look for descending points on both sides
                size_t j_left = i, j_right = i;
                while (j_left > 0 && yhat[j_left] > yhat[j_left - 1]) --j_left;
                while (j_right + 1 < n && yhat[j_right] > yhat[j_right + 1]) ++j_right;

                double d_left  = yhat[i] - yhat[j_left];
                double d_right = yhat[i] - yhat[j_right];

                if (j_left == i) {
                    // No descent on the left — use right
                    e.depth = d_right;
                    e.depth_idx = j_right;
                } else if (j_right == i) {
                    // No descent on the right — use left
                    e.depth = d_left;
                    e.depth_idx = j_left;
                } else {
                    // Both sides descend — take shallower
                    if (d_left < d_right) {
                        e.depth = d_left;
                        e.depth_idx = j_left;
                    } else {
                        e.depth = d_right;
                        e.depth_idx = j_right;
                    }
                }
            } else {
                // Local minimum: look for ascending points on both sides
                size_t j_left = i, j_right = i;
                while (j_left > 0 && yhat[j_left] < yhat[j_left - 1]) --j_left;
                while (j_right + 1 < n && yhat[j_right] < yhat[j_right + 1]) ++j_right;

                double d_left  = yhat[j_left] - yhat[i];
                double d_right = yhat[j_right] - yhat[i];

                if (j_left == i) {
                    // No ascent on the left — use right
                    e.depth = d_right;
                    e.depth_idx = j_right;
                } else if (j_right == i) {
                    // No ascent on the right — use left
                    e.depth = d_left;
                    e.depth_idx = j_left;
                } else {
                    // Both sides ascend — take shallower
                    if (d_left < d_right) {
                        e.depth = d_left;
                        e.depth_idx = j_left;
                    } else {
                        e.depth = d_right;
                        e.depth_idx = j_right;
                    }
                }
            }
        }

        // Compute total depth and ŷ range
        double total_depth = 0.0;
        double yhat_min = *std::min_element(yhat.begin(), yhat.end());
        double yhat_max = *std::max_element(yhat.begin(), yhat.end());
        double yhat_range = yhat_max - yhat_min;

        for (const auto& e : result.local_extrema) {
            total_depth += e.depth;
        }

        for (auto& e : result.local_extrema) {
            e.rel_depth = (total_depth > 0.0) ? e.depth / total_depth : 0.0;
            e.range_rel_depth = (yhat_range > 0.0) ? e.depth / yhat_range : 0.0;
        }

        // Collect depth ratios
        std::vector<double> rel_depths, range_rel_depths;
        for (const auto& e : result.local_extrema) {
            rel_depths.push_back(e.rel_depth);
            range_rel_depths.push_back(e.range_rel_depth);
        }

#if 0
        // Sort for ECDF-based quantiles
        std::vector<double> sorted_rel = rel_depths;
        std::vector<double> sorted_range = range_rel_depths;
        std::sort(sorted_rel.begin(), sorted_rel.end());
        std::sort(sorted_range.begin(), sorted_range.end());

        // Assign empirical CDF quantiles
        for (auto& e : result.local_extrema) {
            e.q_rel_depth       = std::lower_bound(sorted_rel.begin(), sorted_rel.end(), e.rel_depth) - sorted_rel.begin();
            e.q_range_rel_depth = std::lower_bound(sorted_range.begin(), sorted_range.end(), e.range_rel_depth) - sorted_range.begin();

            e.q_rel_depth       /= static_cast<double>(sorted_rel.size());
            e.q_range_rel_depth /= static_cast<double>(sorted_range.size());
        }
#endif
    }

#if 0
    // --- Harmonic Smoothing ---
    result.harmonic_predictions = result.predictions;

    if (result.local_extrema.size() > 1) {

        // Step 1: Sort extrema by x-coordinate
        std::sort(result.local_extrema.begin(), result.local_extrema.end(),
                  [](const extremum_t& a, const extremum_t& b) { return a.x < b.x; });

        size_t m = result.local_extrema.size();

        // Step 2: Build idx → position map for fast lookup
        std::unordered_map<size_t, size_t> idx_to_pos;
        for (size_t k = 0; k < m; ++k) {
            idx_to_pos[result.local_extrema[k].idx] = k;
        }

        // Step 3: Process cancellable extremum pairs
        std::vector<bool> used(m, false); // track which extrema have been already repaired

        for (size_t k = 0; k < m; ++k) {
            if (used[k]) continue; // skip if already used in a cancellation

            const auto& e1 = result.local_extrema[k];

            auto it = idx_to_pos.find(e1.depth_idx);
            if (it == idx_to_pos.end()) continue; // no matching extremum

            size_t l = it->second;
            const auto& e2 = result.local_extrema[l];

            // Check mutual matching
            if (e2.depth_idx != e1.idx || used[l]) continue;

            // Found cancellable pair (e1, e2)
            used[k] = true;
            used[l] = true;

            size_t idx_start = std::min(e1.idx, e2.idx);
            size_t idx_end   = std::max(e1.idx, e2.idx);

            // Safety: if indices are too close, skip
            if (idx_end <= idx_start + 1) continue;

            double x_start = result.x_sorted[idx_start];
            double x_end   = result.x_sorted[idx_end];
            double y_start = result.predictions[idx_start];
            double y_end   = result.predictions[idx_end];

            if (x_end == x_start)
                continue; // avoid division by zero

            // Harmonic repair: linearly interpolate between (x_start, y_start) and (x_end, y_end)
            double slope = (y_end - y_start) / (x_end - x_start);

            for (size_t i = idx_start; i <= idx_end; ++i) {
                double xi = result.x_sorted[i];
                result.harmonic_predictions[i] = y_start + slope * (xi - x_start);
            }
        }
    }
#endif


    // --- Harmonic Smoothing ---
    result.harmonic_predictions = result.predictions;
    if (result.local_extrema.size() > 2) {

        // Step 1: Sort extrema by x-coordinate
        std::sort(result.local_extrema.begin(), result.local_extrema.end(),
                  [](const extremum_t& a, const extremum_t& b) { return a.x < b.x; });

        // Step 2: Parameters
        // const double small_depth_threshold = 0.05;  // you can tune this threshold
        // const double depth_similarity_tol  = 0.001;  // allow 1% difference to call depths "similar"

#if 0
        std::unordered_map<size_t, size_t> idx_to_pos;
        for (size_t k = 0; k < m; ++k) {
            idx_to_pos[result.local_extrema[k].idx] = k;
        }

        auto it = idx_to_pos.find(e1.depth_idx);
        if (it != idx_to_pos.end()) {
            const auto& e2 = result.local_extrema[it->second];
            if (e2.depth_idx == e1.idx) {
                // (e1, e2) cancellable pair
            }
        }
#endif

        // Step 3: Scan triplets
        size_t m = result.local_extrema.size();
        for (size_t k = 0; k + 2 < m; ++k) {
            const auto& e1 = result.local_extrema[k];
            const auto& e2 = result.local_extrema[k + 1];
            const auto& e3 = result.local_extrema[k + 2];

            bool similar_e1e2 = std::abs(e1.range_rel_depth - e2.range_rel_depth) <= depth_similarity_tol * std::max(e1.range_rel_depth, e2.range_rel_depth);
            bool similar_e2e3 = std::abs(e2.range_rel_depth - e3.range_rel_depth) <= depth_similarity_tol * std::max(e2.range_rel_depth, e3.range_rel_depth);

            if (
                !(similar_e1e2 || similar_e2e3) ||
                (similar_e1e2 && (e1.range_rel_depth > small_depth_threshold || e2.range_rel_depth > small_depth_threshold)) ||
                (similar_e2e3 && (e2.range_rel_depth > small_depth_threshold || e3.range_rel_depth > small_depth_threshold))
                ) {
                continue;
            }

            // Default choice: use (e1, e2, e3)
            size_t idx_start = e1.idx;
            size_t idx_mid   = e2.idx;
            size_t idx_end   = e3.idx;

            // Try to refine to (e2, e3, e4) if possible
            if (k + 3 < m) {
                const auto& e4 = result.local_extrema[k + 3];

                double x2 = result.x_sorted[e2.idx];
                double x3 = result.x_sorted[e3.idx];
                double x1 = result.x_sorted[e1.idx];
                double x4 = result.x_sorted[e4.idx];
                double x23_mean = 0.5 * (x2 + x3);

                double d1 = std::abs(x23_mean - x1);
                double d4 = std::abs(x23_mean - x4);

                if (d4 < d1) {
                    // (e2, e3, e4) is a better triplet
                    idx_start = e2.idx;
                    idx_mid   = e3.idx;
                    idx_end   = e4.idx;
                    ++k; // advance extra step since we "shifted" forward
                }
            }

            // Repair between idx_start and idx_end, with optional dip at idx_mid
            double x_start = result.x_sorted[idx_start];
            double x_end   = result.x_sorted[idx_end];
            double y_start = result.predictions[idx_start];
            double y_end   = result.predictions[idx_end];

            if (x_end == x_start)
                continue;  // safety

            if (std::abs(y_start - y_end) < 1e-6) {
                if (idx_mid > idx_start && idx_mid < idx_end) {
                    double x_mid = result.x_sorted[idx_mid];
                    double dip_frac = 0.01;
                    double y_dip = y_start - dip_frac * (std::abs(y_start) + 1.0);

                    double slope1 = (y_dip - y_start) / (x_mid - x_start);
                    for (size_t i = idx_start; i <= idx_mid; ++i) {
                        double xi = result.x_sorted[i];
                        result.harmonic_predictions[i] = y_start + slope1 * (xi - x_start);
                    }

                    double slope2 = (y_end - y_dip) / (x_end - x_mid);
                    for (size_t i = idx_mid + 1; i <= idx_end; ++i) {
                        double xi = result.x_sorted[i];
                        result.harmonic_predictions[i] = y_dip + slope2 * (xi - x_mid);
                    }
                } else {
                    // fallback linear
                    double slope = (y_end - y_start) / (x_end - x_start);
                    for (size_t i = idx_start; i <= idx_end; ++i) {
                        double xi = result.x_sorted[i];
                        result.harmonic_predictions[i] = y_start + slope * (xi - x_start);
                    }
                }
            } else {
                // Normal linear interpolation
                double slope = (y_end - y_start) / (x_end - x_start);
                for (size_t i = idx_start; i <= idx_end; ++i) {
                    double xi = result.x_sorted[i];
                    result.harmonic_predictions[i] = y_start + slope * (xi - x_start);
                }
            }

            // Skip next two extrema to avoid overlapping triplet repairs
            k += 1;
        }
    }

    // --- Harmonic Predictions Local Extrema & Depths ---
    {
        result.harmonic_predictions_local_extrema.clear();  // vector<extremum_t> field in amagelo_t

        const auto& yhat = result.harmonic_predictions;
        const auto& x = result.x_sorted;
        size_t n = yhat.size();

        // Step 1: Detect local extrema (including endpoints)
        for (size_t i = 0; i < n; ++i) {
            double curr = yhat[i];
            bool is_extremum = false;
            bool is_max = false;

            if (i == 0 && n > 1) {
                // Left endpoint
                if (curr > yhat[i + 1]) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < yhat[i + 1]) {
                    is_extremum = true;
                    is_max = false;
                }
            } else if (i == n - 1 && n > 1) {
                // Right endpoint
                if (curr > yhat[i - 1]) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < yhat[i - 1]) {
                    is_extremum = true;
                    is_max = false;
                }
            } else if (i > 0 && i + 1 < n) {
                // Interior points
                double prev = yhat[i - 1], next = yhat[i + 1];
                if (curr > prev && curr > next) {
                    is_extremum = true;
                    is_max = true;
                } else if (curr < prev && curr < next) {
                    is_extremum = true;
                    is_max = false;
                }
            }

            if (is_extremum) {
                result.harmonic_predictions_local_extrema.emplace_back(i, x[i], curr, is_max, 0.0);
            }
        }

        // Step 2: Compute depths of local extrema
        for (auto& e : result.harmonic_predictions_local_extrema) {
            size_t i = e.idx;

            if (e.is_max) {
                // Local maximum: look for descending points on both sides
                size_t j_left = i, j_right = i;
                while (j_left > 0 && yhat[j_left] > yhat[j_left - 1]) --j_left;
                while (j_right + 1 < n && yhat[j_right] > yhat[j_right + 1]) ++j_right;

                double d_left  = yhat[i] - yhat[j_left];
                double d_right = yhat[i] - yhat[j_right];

                if (j_left == i) {
                    // No descent on the left — use right
                    e.depth = d_right;
                    e.depth_idx = j_right;
                } else if (j_right == i) {
                    // No descent on the right — use left
                    e.depth = d_left;
                    e.depth_idx = j_left;
                } else {
                    // Both sides descend — take shallower
                    if (d_left < d_right) {
                        e.depth = d_left;
                        e.depth_idx = j_left;
                    } else {
                        e.depth = d_right;
                        e.depth_idx = j_right;
                    }
                }
            } else {
                // Local minimum: look for ascending points on both sides
                size_t j_left = i, j_right = i;
                while (j_left > 0 && yhat[j_left] < yhat[j_left - 1]) --j_left;
                while (j_right + 1 < n && yhat[j_right] < yhat[j_right + 1]) ++j_right;

                double d_left  = yhat[j_left] - yhat[i];
                double d_right = yhat[j_right] - yhat[i];

                if (j_left == i) {
                    // No ascent on the left — use right
                    e.depth = d_right;
                    e.depth_idx = j_right;
                } else if (j_right == i) {
                    // No ascent on the right — use left
                    e.depth = d_left;
                    e.depth_idx = j_left;
                } else {
                    // Both sides ascend — take shallower
                    if (d_left < d_right) {
                        e.depth = d_left;
                        e.depth_idx = j_left;
                    } else {
                        e.depth = d_right;
                        e.depth_idx = j_right;
                    }
                }
            }
        }

        // Compute total depth and ŷ range
        double total_depth = 0.0;
        double yhat_min = *std::min_element(yhat.begin(), yhat.end());
        double yhat_max = *std::max_element(yhat.begin(), yhat.end());
        double yhat_range = yhat_max - yhat_min;

        for (const auto& e : result.harmonic_predictions_local_extrema) {
            total_depth += e.depth;
        }

        for (auto& e : result.harmonic_predictions_local_extrema) {
            e.rel_depth = (total_depth > 0.0) ? e.depth / total_depth : 0.0;
            e.range_rel_depth = (yhat_range > 0.0) ? e.depth / yhat_range : 0.0;
        }

        // Collect depth ratios
        std::vector<double> rel_depths, range_rel_depths;
        for (const auto& e : result.harmonic_predictions_local_extrema) {
            rel_depths.push_back(e.rel_depth);
            range_rel_depths.push_back(e.range_rel_depth);
        }
    }

    //
    // Compute Change-Scaled Monotonicity Index
    //
    {
        const auto& y = result.harmonic_predictions;
        const auto& x = result.x_sorted;
        size_t n = y.size();

        if (n < 2) {
            result.change_scaled_monotonicity_index = 0.0;
        } else {
            std::vector<double> delta_y(n-1);
            std::vector<double> delta_x(n-1);
            std::vector<int> orientation(n-1);
            std::vector<double> weights(n-1);

            for (size_t i = 1; i < n; ++i) {
                delta_y[i-1] = y[i] - y[i-1];
                delta_x[i-1] = x[i] - x[i-1];
                orientation[i-1] = (delta_y[i-1] > 0) ? 1 : (delta_y[i-1] < 0) ? -1 : 0;
                weights[i-1] = std::abs(delta_y[i-1]) * delta_x[i-1];
            }

            double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
            if (total_weight < 1e-12) {
                result.change_scaled_monotonicity_index = 0.0;  // If essentially no change, define CSMI = 0
            } else {
                double csmi = 0.0;
                for (size_t i = 0; i < n-1; ++i) {
                    csmi += orientation[i] * (weights[i] / total_weight);
                }
                result.change_scaled_monotonicity_index = csmi;
            }
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

        if (grid_vertices.find(current_vertex) != grid_vertices.end()) {
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

    {
        // Check if the lower bound already satisfies the condition
        auto [orig_low, grid_low] = find_original_vertices_within_radius(grid_vertex, lower_bound);
        if (orig_low.size() >= domain_min_size)
            return lower_bound;
    }

    {
        // Check if the upper bound fails to satisfy the condition
        auto [orig_high, grid_high] = find_original_vertices_within_radius(grid_vertex, upper_bound);

        if (orig_high.size() < domain_min_size) {
            Rprintf("\n---------------------\nERROR\ngrid_vertex: %zu\n"
                    "Not enough vertices (found: %zu, needed: %zu) within upper_bound: %.4f\n",
                    grid_vertex + 1, orig_high.size(), domain_min_size, upper_bound);
            REPORT_ERROR("ERROR: Insufficient vertices in maximum radius\n");
        }
    }

    // Binary search for the minimum radius
    double thld = precision * graph_diameter;
    while ((upper_bound - lower_bound) > thld) {

        double mid = 0.5 * (lower_bound + upper_bound);

        auto [orig_mid, grid_mid] = find_original_vertices_within_radius(grid_vertex, mid);

        if (orig_mid.size() == domain_min_size) {
            upper_bound = mid;
            break;
        } else if (orig_mid.size() > domain_min_size) {
            // If condition is satisfied, try a smaller radius
            upper_bound = mid;

        } else {
            // If condition is not satisfied, try a larger radius
            lower_bound = mid;
        }
    }

    // Return the upper bound as the minimum radius that satisfies the condition
    return upper_bound;
}
