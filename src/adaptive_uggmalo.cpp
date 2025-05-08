/**
 * @brief Implementation of Adaptive UGGMALO algorithm
 * Provides model-averaged local linear regression on graph structures with
 * adaptive bandwidth selection and sophisticated model selection strategies.
 */

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <optional>
#include <queue>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>  // For std::reverse
#include <sstream>
#include <string>
#include <limits>
#include <map>
#include <utility>
#include <numeric>    // For std::accumulate
#include <execution>
#include <atomic>
#include <mutex>
#include <random>     // For std::mt19937
#include <chrono>
#include <omp.h>

#include "uniform_grid_graph.hpp"
#include "ulm.hpp"
#include "graph_utils.hpp"        // For get_grid_diameter()
#include "centered_paths.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_shortest_path.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h"          // For REPORT_ERROR()
#include "kernels.h"
#include "sampling.h"             // For C_runif_simplex()
#include "predictive_errors.hpp"
#include "wasserstein_perm_test.hpp"
#include "ext_ulm_priority_queue.hpp"
#include "graph_maximal_packing.hpp"
#include "bandwidth_utils.hpp"    // For get_candidate_bws()

extern "C" {
    SEXP S_adaptive_uggmalo(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_min_path_size,
        SEXP s_n_grid_vertices,
        SEXP s_n_bws,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_max_iterations,
        SEXP s_precision,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_tolerance,
        SEXP s_n_bb,
        SEXP s_cri_probability,
        SEXP s_n_perms,
        SEXP s_blending_coef,
        SEXP s_verbose
        );
}

/**
 * @brief Implements the Adaptive Uniform Grid Graph Model-Averaged LOcal Linear regression (Adaptive-UGGMALO) algorithm
 *
 * @details This function implements an advanced algorithm for analyzing weighted graphs using local path
 * linear models with model averaging and adaptive bandwidth selection. The algorithm combines uniform
 * grid representation with local linear modeling to capture both global and local structure in graph data.
 *
 * The algorithm proceeds through several key phases:
 * 1. Graph Analysis and Grid Construction
 *    - Computes graph diameter to determine maximum bandwidth
 *    - Creates uniform grid representation of input graph with specified number of vertices
 *    - Maps original vertices to grid vertices using snap tolerance
 *
 * 2. Model Generation and Selection
 *    - For each grid vertex:
 *      - Computes reachability maps within maximum bandwidth
 *      - Constructs paths through grid and original vertices
 *      - Ensures minimum path size requirements are met
 *      - Generates candidate bandwidths based on path distributions
 *
 * 3. Local Linear Modeling
 *    - For each bandwidth and valid path:
 *      - Fits weighted linear models to path data
 *      - Computes leave-one-out cross-validation errors
 *      - Maintains queue of models sorted by prediction error
 *
 * 4. Model Averaging
 *    - Combines predictions from multiple models using weighted averaging
 *    - Supports both standard and mean-error-weighted averaging schemes
 *    - Generates predictions for both original and grid vertices
 *
 * 5. Optional Statistical Analysis
 *    - Bayesian bootstrap estimation for confidence intervals
 *    - Permutation testing for statistical significance
 *    - Wasserstein distance-based effect size calculation
 *
 * @param adj_list Adjacency list representation of the graph, where adj_list[i] contains indices of vertices adjacent to vertex i
 * @param weight_list Corresponding edge weights, where weight_list[i][j] is the weight of edge from vertex i to adj_list[i][j]
 * @param y Response values at each vertex in the original graph
 * @param start_vertex Index of the vertex to start graph traversal
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param n_grid_vertices Number of vertices in the uniform grid representation
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param snap_tolerance Maximum distance for snapping original vertices to grid vertices
 * @param dist_normalization_factor Factor for normalizing distances in the graph
 * @param kernel_type Type of kernel function used for weight calculation (e.g., Gaussian, triangular)
 * @param tolerance Convergence tolerance for linear model fitting
 * @param n_bb Number of Bayesian bootstrap iterations (0 for no bootstrap)
 * @param cri_probability Probability level for confidence intervals in bootstrap
 * @param n_perms Number of permutation test iterations (0 for no permutation testing)
 * @param blending_coef Number between 0 and 1 allowing for smooth interpolation between position-based weights and mean-error times position-based weights using the formula effective_weight = x.weight * pow(x.mean_error, blending_coef);
 * @param verbose Whether to print progress information
 *
 * @return adaptive_uggmalo_result_t structure containing:
 *         - graph_diameter: Computed diameter of the input graph
 *         - grid_opt_bw: Optimal bandwidth for each grid vertex
 *         - predictions: Model-averaged predictions for original graph vertices
 *         - grid_predictions: Model-averaged predictions for grid vertices
 *         - Bootstrap results (if n_bb > 0):
 *           - bb_predictions: Bootstrap replication predictions
 *           - cri_lower: Lower confidence interval bounds
 *           - cri_upper: Upper confidence interval bounds
 *         - Permutation test results (if n_perms > 0):
 *           - null_predictions: Predictions under null hypothesis
 *           - permutation_tests: Statistical test results including p-values and effect sizes
 *
 * @note Performance Characteristics:
 * - Time Complexity: O(V * P * M * B) where:
 *   - V = number of grid vertices
 *   - P = average number of paths per vertex
 *   - M = maximum iterations in linear regression
 *   - B = number of bandwidths evaluated
 * - Space Complexity: O(V * N + B * P) where:
 *   - N = number of original vertices
 *   - B = number of bandwidths
 *   - P = average path size
 *
 * @note Implementation Details:
 * - Parallelization: Uses parallel execution for bootstrap and permutation iterations
 * - Thread Safety: Implements thread-local random number generation
 * - Memory Management: Preallocates vectors to minimize reallocations
 * - Error Handling:
 *   - Validates grid vertex count matches specification
 *   - Handles non-convergent fits by skipping affected models
 *   - Uses original response values as fallback for failed predictions
 *   - Implements error reporting for critical failures
 *
 * @warning
 * - Input graph must be connected for proper diameter calculation
 * - Edge weights must be non-negative
 * - n_grid_vertices should be chosen considering computational resources
 * - Bootstrap and permutation testing can be computationally intensive
 *
 * @see vertex_wasserstein_perm_test For details on permutation testing
 * @see bb_cri For details on bootstrap confidence interval calculation
 * @see create_uniform_grid_graph For grid construction algorithm
 */
struct adaptive_uggmalo_result_t {
    double graph_diameter;

    std::vector<double> predictions;      ///< predictions[i] is an averaged-model estimate of E(Y|G) at the i-th vertex of the original (before uniform grid construction) graph G
    std::vector<double> errors;           ///< errors[i] is an averaged-model estimate of the prediction error at the i-th vertex of the original (before uniform grid construction) graph
    std::vector<double> scale;            ///< scale[i] is a local scale at the i-th vertex of the original (before uniform grid construction) graph G, which is an approximate radius of a disk in G where predictions is well approximated by a linear model

    std::vector<double> grid_opt_bw;      ///< grid_opt_bw[grid_vertex] = the optimal bandwidth over all models at this vertex; this gives a local scale at each grid vertex
    std::unordered_map<size_t, double> grid_predictions_map; ///< model-averaged predictions at the grid vertices; Note that for the Morse-Smale cells construction we only need estimate of E(y|G) over grid points; this would require reformulation of the notion of tau-gradient flow

    // Bootstrap fields (if n_bb > 0)
    std::vector<std::vector<double>> bb_predictions;
    std::vector<double> cri_lower;
    std::vector<double> cri_upper;

    // Permutation fields (if n_perms > 0)
    std::vector<std::vector<double>> null_predictions;
    std::vector<double> null_predictions_cri_lower;
    std::vector<double> null_predictions_cri_upper;

    // Permutation Tests
    std::optional<vertex_wasserstein_perm_test_results_t> permutation_tests;
};

adaptive_uggmalo_result_t adaptive_uggmalo(
    const uniform_grid_graph_t& grid_graph,
    const std::vector<double>& y,
    size_t min_path_size,
    size_t n_bws,
    double min_bw_factor,
    double max_bw_factor,
    double dist_normalization_factor,
    size_t kernel_type,
    double tolerance,
    size_t n_bb,
    double cri_probability,
    size_t n_perms,
    double blending_coef,
    bool verbose
    ) {

    // std::for_each flag
    auto exec = std::execution::seq;
    // auto exec = std::execution::par_unseq;

    initialize_kernel(kernel_type, 1.0);

    // Initialize result structure
    adaptive_uggmalo_result_t result;

    if (grid_graph.graph_diameter < 0) {
        REPORT_ERROR("grid_graph.graph_diameter: %.4f\n", grid_graph.graph_diameter);
    }

    result.graph_diameter = grid_graph.graph_diameter;

    	// Define minimum and maximum bandwidth as a fraction of graph diameter
    double max_bw = max_bw_factor * grid_graph.graph_diameter;
    double min_bw = min_bw_factor * grid_graph.graph_diameter;
    // Ensure the min_bw is at least the maximum packing radius
    // This guarantees all vertices are covered by at least one ball
    min_bw = std::max(min_bw, grid_graph.max_packing_radius);

    if (verbose) {
        Rprintf("In adaptive_uggmalo()\n");

        Rprintf("min_bw_factor: %f\n", min_bw_factor);
        Rprintf("min_bw: %f\n", min_bw);

        Rprintf("max_bw_factor: %f\n", max_bw_factor);
        Rprintf("max_bw: %f\n", max_bw);

        Rprintf("graph_diameter: %f\n", result.graph_diameter);
        Rprintf("n_bws: %zu\n", n_bws);
        Rprintf("n_vertices(grid_graph): %zu\n", grid_graph.adjacency_list.size());
        Rprintf("grid_graph.grid_vertices.size(): %zu\n", grid_graph.grid_vertices.size());
    }

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    size_t n_vertices = grid_graph.n_original_vertices;
    result.predictions.resize(n_vertices, INFINITY);
    result.errors.resize(n_vertices, INFINITY);
    result.scale.resize(n_vertices, INFINITY);

    // Create models for a single grid vertex.
    // Returns: Min Heap: Queue of models sorted in the ascending order by mean error
    auto create_grid_vertex_models = [&](
        size_t grid_vertex,
        const std::vector<double>& y,
        const std::optional<std::vector<double>>& weights
        ) {

        // Fit local models over paths of different size
        ext_ulm_priority_queue model_queue; // a priority queue of models sorted in the increasing order of their mean predition error
        // the insertion operation into model_queue:
        // If there are no models whose path vertices intersect query model's path vertices - insert it
        // If there are models with intersecting vertices AND the new model has lower mean_error than ALL of them - remove them and insert the new model
        // If there are models with intersecting vertices AND at least one of them has lower mean_error than the new model - do not insert the new model

        // Run Bijkstra's algorithm with the restriction on the path length and
        // returns the corresponding data structures for the shortest paths
        // staring at the given grid vertex. Note that the search is done in the
        // grid graph so both the original and grid vertices are part of the
        // search and results.

        // Reconstruct the paths discovered in reachability_map stratifying the
        // findings by the original graph vertices and grid vertices.

        // Find the minimum bandwidth that ensures at least one geodesic ray has the minimum required size
        // NOTE that this is a conservative estimate as for the composite paths this value may be even smaller
        double precision = 0.001;
        double min_min_bw = grid_graph.find_minimum_bandwidth(
            grid_vertex,
            min_bw,
            max_bw,
            min_path_size,
            precision
            );
        if (min_bw < min_min_bw) {
            min_bw = min_min_bw;
        }
        if (min_bw >= max_bw) {
            REPORT_ERROR("Minimum bandwidth (%f) >= maximum bandwidth (%f)\n",
                         min_bw, max_bw_factor * grid_graph.graph_diameter);
        }

        bool log_grid = false;

        std::vector<double> candidate_bws = get_candidate_bws(
            min_bw,
            max_bw,
            n_bws,
            log_grid,
            precision
            );

        double final_radius_value = 0.0;
        size_t min_num_grid_vertices = 0;
        double initial_radius = max_bw;

        std::vector<compose_path_t> composite_paths = grid_graph.find_and_filter_composite_paths(
            grid_vertex,
            initial_radius,
            min_path_size,
            min_num_grid_vertices,
            result.graph_diameter,
            &final_radius_value
            );

        if (composite_paths.size() == 0) {
            REPORT_ERROR("ERROR: No composite paths detected in create_grid_vertex_models lambda grid_vertex: %zu\tcomposite_paths.size(): %zu\n", // <<----
                         grid_vertex, composite_paths.size());
        }

        if (final_radius_value > initial_radius) {
            Rprintf("Radius was increased from %.3f to %.3f\n", initial_radius, final_radius_value);
        }


        for (const auto& bw : candidate_bws) {

            auto bw_composite_paths = grid_graph.get_bw_composite_paths(bw,
                                                                        composite_paths,
                                                                        min_path_size,
                                                                        min_num_grid_vertices);

            if (bw_composite_paths.size() == 0) {
                REPORT_WARNING("bw: %.3f\tbw_composite_paths.size() is 0\n",bw);
            }

            for (auto& path : bw_composite_paths) {

                auto xyw_path = grid_graph.get_xyw_data(y, path, dist_normalization_factor);

                // Fit weighted linear model to the current path
                if (weights) {
                    // multiply path.w_path by the corresponding values of weights
                    for (size_t i = 0; i < path.vertices.size(); i++) {
                        xyw_path.w_path[i] *= (*weights)[path.vertices[i]];
                    }
                }

                ulm_t fit_result = ulm(
                    xyw_path.x_path.data(),
                    xyw_path.y_path.data(),
                    xyw_path.w_path,
                    y_binary,
                    tolerance
                    );

                // Create the extended object from fit_result
                ext_ulm_t extended_result(fit_result);
                extended_result.bw            = bw;

                #if 0
                extended_result.vertices      = path.vertices;
                extended_result.w_path        = xyw_path.w_path;
                extended_result.grid_w_path   = xyw_path.grid_w_path;
                extended_result.grid_vertices = path.grid_vertices;
                #else
                extended_result.vertices = std::move(path.vertices);
                path.vertices.clear(); // Explicitly reset after move
                extended_result.w_path        = std::move(xyw_path.w_path);
                xyw_path.w_path.clear();
                extended_result.grid_w_path   = std::move(xyw_path.grid_w_path);
                xyw_path.grid_w_path.clear();
                extended_result.grid_vertices = std::move(path.grid_vertices);
                path.grid_vertices.clear(); // Explicitly reset after move
                #endif

                // Calculate the mean of loocv_brier_errors vector
                // We'll use std::accumulate to sum all elements and divide by the vector size
                double sum = std::accumulate(fit_result.errors.begin(),
                                             fit_result.errors.end(),
                                             0.0);  // Use 0.0 to ensure floating-point arithmetic
                extended_result.mean_error = sum / fit_result.errors.size();

                extended_result.grid_predictions = fit_result.predict(xyw_path.grid_x_path);

                model_queue.custom_insert(extended_result);
            }
        }

        // auto model = model_queue.top();
        // Rprintf("grid_vertex: %zu n(models): %zu  bw: %.3f\n", grid_vertex, model_queue.size(), model.bw);

        return model_queue;
    }; // END OF create_grid_vertex_models()


    // weight/prediction/error/mean_error/be struct needed for mode averaging and local scale estimation; we use it to record weight/prediction/error/bw of the given vertex in each model where the vertext is in the support of the model
    struct wpe_t {
        double weight;
        double prediction;
        double error;
        double mean_error;
        double bw;

        // Constructor needed for emplace_back(x,y,z)
        wpe_t(double w, double p, double e, double me, double bw)
            : weight(w), prediction(p), error(e), mean_error(me), bw(bw) {}
    };

    struct grid_wpe_t {
        double weight;
        double prediction;
        double mean_error;
        double bw;

        // Constructor needed for emplace_back(x,y,z)
        grid_wpe_t(double w, double p, double me, double bw)
            : weight(w), prediction(p), mean_error(me), bw(bw) {}
    };


    auto average_models = [&blending_coef,&n_vertices,&verbose,&grid_graph](
        std::unordered_map<size_t, ext_ulm_priority_queue>& grid_vertex_models_map,
        const std::vector<double>& y,
        std::vector<double>& predictions,
        std::vector<double>& errors,
        std::vector<double>& scale,
        std::optional<std::unordered_map<size_t, double>>& grid_predictions_map
        ) {

        bool process_errors = errors.size() == predictions.size();
        bool process_scale  = scale.size() == predictions.size();

        std::vector<std::vector<wpe_t>> wpe(n_vertices); // wpe[i] stores a vector of {weight, prediction, error, mean_error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions
        std::unordered_map<size_t, std::vector<grid_wpe_t>> grid_wpe_map;

        // Rprintf("In average_models lambda grid_vertex_models_map.size(): %zu\n", grid_vertex_models_map.size());

        for (auto& [i, models_queue] : grid_vertex_models_map) {

            // auto model = models_queue.top();
            // Rprintf("model.grid_vertices.size(): %zu\n", model.grid_vertices.size());
            // if (models_queue.empty()) {
            //     Rprintf("models_queue is empty for i: %zu\n",i);
            // }

            while (!models_queue.empty()) {
                auto model = models_queue.top(); // ext_ulm_t object
                // Store weight, prediction and error for the given vertex from the given 'model' in wpe[ model.vertices[i] ]
                for (size_t i = 0; i < model.vertices.size(); ++i) {
                    wpe[ model.vertices[i] ].emplace_back(model.w_path[i], model.predictions[i], model.errors[i], model.mean_error, model.bw);
                }

                if (grid_predictions_map) {
                    for (size_t i = 0; i < model.grid_vertices.size(); ++i) {
                        size_t grid_idx = model.grid_vertices[i];
                        grid_wpe_map[grid_idx].emplace_back(model.grid_w_path[i], model.grid_predictions[i], model.mean_error, model.bw);
                    }
                }

                models_queue.pop();
            }
        }

        // Model averaging over the vertices of the original graph
        for (size_t i = 0; i < n_vertices; i++ ) {

            double prediction_sum = 0.0;
            double weight_sum     = 0.0;
            double error_sum      = 0.0;
            double local_scale    = 0.0; //INFINITY;
            double effective_weight;

            for (const auto& x : wpe[i]) {
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
                    // Method 1: Linear interpolation of weights
                    effective_weight = (1.0 - blending_coef) * x.weight + blending_coef * (x.mean_error * x.weight);
                    // Method 2: Power-based scaling
                    //effective_weight = x.weight * pow(x.mean_error, blending_coef);
                }

                prediction_sum += effective_weight * x.prediction;
                weight_sum     += effective_weight;

                if (process_errors) error_sum += effective_weight * x.error;
                // if (process_scale && x.bw < local_scale) local_scale = x.bw; // min bw method
                if (process_scale) local_scale += effective_weight * x.bw; // weighted mean bw method
            }

            if (weight_sum > 0) {
                if (weight_sum > 1e-10) {
                    predictions[i] = prediction_sum / weight_sum;
                    if (process_errors)
                        errors[i] = error_sum / weight_sum;
                    if (process_scale)
                        scale[i] = local_scale / weight_sum;
                } else {
                    if (verbose) {
                        REPORT_WARNING("Very small weight sum encountered for vertex %zu. Setting predictions[i] to y[i]\n", i);
                    }
                    predictions[i] = std::numeric_limits<double>::quiet_NaN();
                }
            } else {
                REPORT_WARNING("Weight sum = 0 for vertex %zu in predictions\n", i);
                predictions[i] = y[i];
            }
        }

        if (grid_predictions_map) {

            //Rprintf("grid_wpe_map.size(): %zu\n", grid_wpe_map.size());

            // Model averaging over the grid vertices
            for (auto& [i, wpe_entries] : grid_wpe_map) {

                //Rprintf("i: %zu\twpe_entries.size(): %zu\n", i, wpe_entries.size());

                auto x = wpe_entries[0];
                (*grid_predictions_map)[i] = x.prediction;

                if (wpe_entries.size() < 2) {
                    continue;
                }

                double prediction_sum = 0.0;
                double weight_sum = 0.0;
                double effective_weight;

                for (const auto& x : wpe_entries) {

                    if (std::isnan(x.weight) || (blending_coef > 0 && std::isnan(x.mean_error))) {
                        continue;
                    }

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
                        // Method 1: Linear interpolation of weights
                        effective_weight = (1.0 - blending_coef) * x.weight + blending_coef * (x.mean_error * x.weight);
                        // Method 2: Power-based scaling
                        //effective_weight = x.weight * pow(x.mean_error, blending_coef);
                    }

                    prediction_sum += effective_weight * x.prediction;
                    weight_sum += effective_weight;
                }

                if (weight_sum > 0) {
                    if (weight_sum > 1e-10) {
                        (*grid_predictions_map)[i] = prediction_sum / weight_sum;
                    } else {
                        REPORT_ERROR("Weight sum < 1e-10 for vertex %zu in grid predictions\n", i);
                    }
                } else {
                    Rprintf("\n\nERROR: weight_sum: %f for vertex %zu in grid predictions\n", weight_sum, i);
                    Rprintf("grid_wpe_map[%zu].size(): %zu\n", i, wpe_entries.size());

                    size_t counter = 0;
                    for (const auto& x : wpe_entries) {
                        Rprintf("%zu\tx.weight: %f\tx.prediction: %f\n", counter++, x.weight, x.prediction);
                    }

                    REPORT_ERROR("Weight sum = 0 for vertex %zu in grid predictions\n", i);
                }

                // Rprintf("grid_predictions_map.size(): %zu\n", (*grid_predictions_map).size());
            }
        } // END OF if (grid_predictions_map)
    }; // END OF auto average_models

    // Calculates predictions using grid vertex models. Returns void. Accepts memory allocated vector of predictions as a 'return value'
    auto grid_models_predictions = [&](
        const std::vector<double>& y,
        const std::optional<std::vector<double>>& weights,
        std::vector<double>& predictions,
        std::vector<double>& errors,
        std::vector<double>& scale,
        std::optional<std::unordered_map<size_t, double>>& grid_predictions_map
        ) {

        std::unordered_map<size_t, ext_ulm_priority_queue> grid_vertex_models_map;

        for (const auto& grid_vertex : grid_graph.grid_vertices) {
            ext_ulm_priority_queue models = create_grid_vertex_models(
                grid_vertex,
                y,
                weights
                );
            grid_vertex_models_map[grid_vertex] = models;
        }

         average_models(grid_vertex_models_map, y, predictions, errors, scale, grid_predictions_map);
    };

    // Create thread-local random number generators
    const unsigned int num_threads = std::thread::hardware_concurrency();
    std::vector<std::mt19937> thread_rngs(num_threads);
    for (unsigned int i = 0; i < num_threads; ++i) {
        std::random_device rd;
        thread_rngs[i].seed(rd());
    }

    std::mutex rng_mutex; // Mutex for thread-safe random number generation
    std::vector<double> empty_errors;
    std::vector<double> empty_scale;

    if (n_bb > 0) {
        // Perform Bayesian bootstrap estimation

        // Initialize storage for predictions and errors
        result.bb_predictions.resize(n_bb);
        for (auto& predictions : result.bb_predictions) {
            predictions.resize(n_vertices);
        }

        // Create indices for parallel iteration
        std::vector<int> bb_indices(n_bb);
        std::iota(bb_indices.begin(), bb_indices.end(), 0);

        // Progress tracking
        std::atomic<int> bootstrap_counter{0};

        std::for_each(exec,
                      bb_indices.begin(),
                      bb_indices.end(),
                      [&](int iboot) {
                          // Get thread-local RNG
                          const int thread_id = omp_get_thread_num() % num_threads;
                          auto& local_rng = thread_rngs[thread_id];

                          // Generate weights using thread-local RNG
                          std::vector<double> weights(n_vertices);
                          {
                              // Use thread-local RNG to generate weights
                              std::gamma_distribution<double> gamma(1.0, 1.0);
                              double sum = 0.0;
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] = gamma(local_rng);
                                  sum += weights[i];
                              }
                              // Normalize weights
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] /= sum;
                              }
                          }

                          std::optional<std::vector<double>> opt_weights(weights);
                          std::optional<std::unordered_map<size_t, double>> opt_nullptr_grid_pred;
                          grid_models_predictions(
                              y,
                              opt_weights,
                              result.bb_predictions[iboot],
                              empty_errors,
                              empty_scale,
                              opt_nullptr_grid_pred
                              );

                          // Thread-safe progress update without mutex
                          if (verbose) {
                              int current_count = ++bootstrap_counter;
                              //if (current_count % progress_chunk == 0) {
                              // Use \r to move cursor to start of line
                              REprintf("\rBootstrap progress: %d%%",
                                       static_cast<int>((100.0 * current_count) / n_bb));
                              //}
                          }
                      });

        bool use_median = true;
        bb_cri_t bb_cri_res = bb_cri(result.bb_predictions, use_median, cri_probability);

        result.predictions = std::move(bb_cri_res.bb_Ey);
        result.cri_lower   = std::move(bb_cri_res.cri_L);
        result.cri_upper   = std::move(bb_cri_res.cri_U);

    } else {
        // Perform standard prediction without bootstrap
        std::optional<std::vector<double>> opt_nullptr;
        std::optional<std::unordered_map<size_t, double>> opt_grid_predictions_map(result.grid_predictions_map);

        grid_models_predictions(
            y,
            opt_nullptr,
            result.predictions,
            result.errors,
            result.scale,
            opt_grid_predictions_map
            );

        if (opt_grid_predictions_map) {
            result.grid_predictions_map = *opt_grid_predictions_map;

            //print_umap(result.grid_predictions_map, "result.grid_predictions_map");
            // Rprintf("result.grid_predictions_map.size(): %zu\n", result.grid_predictions_map.size());
            // Rprintf("n_grid_vertices: %zu\ngrid_graph.grid_vertices.size(): %zu\n",
            //         n_grid_vertices, grid_graph.grid_vertices.size());
        }
    }

    if (n_perms > 0) {
        // Perform permutation testing if requested

        // Initialize permutation test results
        result.null_predictions.resize(n_perms);
        for (auto& pred : result.null_predictions) {
            pred.resize(n_vertices);
        }

        // Create indices for parallel iteration
        std::vector<int> perm_indices(n_perms);
        std::iota(perm_indices.begin(), perm_indices.end(), 0);

        // Create indices for permutation
        std::vector<size_t> indices(n_vertices);
        std::iota(indices.begin(), indices.end(), 0);

        // Atomic counter for tracking progress
        std::atomic<int> permutation_counter{0};

        // Create a random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        // Parallel execution of permutation iterations
        std::for_each(exec,
                      perm_indices.begin(),
                      perm_indices.end(),
                      [&](int iperm) {
                          // Create permuted y vector
                          std::vector<double> y_perm(y);  // Copy original y

                          // Generate permutation in thread-safe manner
                          {
                              std::lock_guard<std::mutex> lock(rng_mutex);
                              std::shuffle(indices.begin(), indices.end(), gen);
                          }

                          // Apply permutation
                          std::vector<double> y_shuffled(n_vertices);
                          for (size_t i = 0; i < n_vertices; ++i) {
                              y_shuffled[i] = y[indices[i]];
                          }

                          // Get thread-local RNG
                          const int thread_id = omp_get_thread_num() % num_threads;
                          auto& local_rng = thread_rngs[thread_id];
                          // Generate weights using thread-local RNG
                          std::vector<double> weights(n_vertices);
                          {
                              // Use thread-local RNG to generate weights
                              std::gamma_distribution<double> gamma(1.0, 1.0);
                              double sum = 0.0;
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] = gamma(local_rng);
                                  sum += weights[i];
                              }
                              // Normalize weights
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] /= sum;
                              }
                          }

                          // Process permuted data using optimal bandwidth
                          std::optional<std::vector<double>> opt_weights(weights);
                          std::optional<std::unordered_map<size_t, double>> opt_nullptr_grid_pred;
                          grid_models_predictions(
                              y_shuffled,
                              opt_weights,
                              result.null_predictions[iperm],
                              empty_errors,
                              empty_scale,
                              opt_nullptr_grid_pred
                              );

                          // Increment counter and update progress
                          if (verbose) {
                              int current_count = ++permutation_counter;
                              // Use \r to move cursor to start of line
                              REprintf("\rPermutation Test Progress: %d%%",
                                       static_cast<int>((100.0 * current_count) / n_perms));
                          }
                      });

        bool use_median = true;
        bb_cri_t null_cri_res = bb_cri(result.null_predictions, use_median, cri_probability);

        result.null_predictions_cri_lower = std::move(null_cri_res.cri_L);
        result.null_predictions_cri_upper = std::move(null_cri_res.cri_U);

        #if 0
        std::optional<std::vector<double>> opt_nullptr;
        std::optional<std::unordered_map<size_t, double>> opt_grid_predictions_map(result.grid_predictions_map);
        result.null_predictions_cri_lower.resize(n_vertices);
        result.null_predictions_cri_upper.resize(n_vertices);

        grid_models_predictions(
            null_cri_res.cri_L,
            opt_nullptr,
            result.null_predictions_cri_lower,
            opt_grid_predictions_map
            );

        grid_models_predictions(
            null_cri_res.cri_U,
            opt_nullptr,
            result.null_predictions_cri_upper,
            opt_grid_predictions_map
            );
        #endif

        if (n_bb > 0) {
            size_t n_bootstraps = 1000;
            double alpha = 0.05;
            result.permutation_tests = vertex_wasserstein_perm_test(
                result.bb_predictions,
                result.null_predictions,
                n_bootstraps,
                alpha
                );
        }
        if (verbose) {
            Rprintf("\nPermutation Test Completed.\n");
        }
    }

    return result;
}

/**
 * @brief R interface function for the Adaptive Uniform Grid Graph Model-Averaged LOcal Linear regression algorithm
 *
 * @details This function serves as the bridge between R and the C++ implementation of the Adaptive-UGGMALO
 * algorithm. It handles the conversion of R data structures to C++ types, manages memory protection for R
 * objects, and transforms the C++ results back into an R-compatible format. The function is designed to be
 * called from R using .Call() interface.
 *
 * The function performs several key operations:
 * 1. Data Structure Conversion
 *    - Converts R lists representing the graph structure to C++ vectors
 *    - Transforms R numeric vectors to C++ vectors while preserving data integrity
 *    - Handles R's 1-based indexing to C++'s 0-based indexing conversion
 *
 * 2. Memory Management
 *    - Implements proper PROTECT/UNPROTECT mechanisms for R objects
 *    - Ensures memory safety during R-to-C++ and C++-to-R conversions
 *    - Manages temporary storage for intermediate calculations
 *
 * 3. Result Processing
 *    - Creates an R list with named components for return values
 *    - Handles optional bootstrap and permutation test results
 *    - Manages conversion of C++ data structures back to R format
 *
 * @param s_adj_list [SEXP] R list of integer vectors representing graph adjacency lists
 * @param s_weight_list [SEXP] R list of numeric vectors containing edge weights
 * @param s_y [SEXP] R numeric vector of response values at each vertex
 * @param s_min_path_size [SEXP] R integer for minimum required path size
 * @param s_n_grid_vertices [SEXP] R integer specifying number of grid vertices
 * @param s_n_bws [SEXP] R integer for number of bandwidth values to evaluate
 * @param s_max_bw_factor [SEXP] R numeric scaling factor for maximum bandwidth
 * @param s_snap_tolerance [SEXP] R numeric tolerance for vertex snapping
 * @param s_dist_normalization_factor [SEXP] R numeric factor for distance normalization
 * @param s_kernel_type [SEXP] R integer specifying kernel function type
 * @param s_tolerance [SEXP] R numeric convergence tolerance for model fitting
 * @param s_n_bb [SEXP] R integer for number of bootstrap iterations
 * @param s_cri_probability [SEXP] R numeric confidence level for intervals
 * @param s_n_perms [SEXP] R integer for number of permutation test iterations
 * @param s_use_mean_error [SEXP] R logical for using mean error in model averaging
 * @param s_verbose [SEXP] R logical controlling progress output
 *
 * @return SEXP R list containing:
 *         - graph_diameter: Numeric value of computed graph diameter
 *         - grid_opt_bw: Numeric vector of optimal bandwidths for grid vertices
 *         - predictions: Numeric vector of predictions for original vertices
 *         - grid_predictions: Numeric vector of predictions for grid vertices
 *         - bb_predictions: Matrix of bootstrap predictions (if n_bb > 0)
 *         - cri_lower: Numeric vector of lower confidence bounds (if n_bb > 0)
 *         - cri_upper: Numeric vector of upper confidence bounds (if n_bb > 0)
 *         - null_predictions: Matrix of permutation test predictions (if n_perms > 0)
 *         - p_values: Numeric vector of vertex-wise p-values (if n_perms > 0)
 *         - effect_sizes: Numeric vector of effect sizes (if n_perms > 0)
 *         - significant_vertices: Logical vector indicating significance (if n_perms > 0)
 *
 * @note Implementation Details:
 * - Uses helper functions for converting between R and C++ data structures
 * - Implements careful error checking for R object types
 * - Maintains proper memory protection counting
 * - Handles matrix transposition for R's column-major format
 *
 * @warning
 * - All input parameters must be of the correct R type (numeric, integer, logical)
 * - Vertex indices in s_adj_list must be 1-based (R convention)
 * - The function assumes input validation has been performed in R
 * - Large graphs may require significant memory for bootstrap/permutation results
 *
 * @note Memory Management:
 * - Each created R object is protected using PROTECT()
 * - Protection stack is balanced using UNPROTECT() before return
 * - Temporary objects are protected during matrix/vector conversions
 *
 * @see adaptive_uggmalo For the underlying C++ implementation
 * @see convert_R_graph_to_vector For adjacency list conversion details
 * @see convert_R_weights_to_vector For weight list conversion details
 */
SEXP S_adaptive_uggmalo(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_path_size,
    SEXP s_n_grid_vertices,
    SEXP s_n_bws,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_max_iterations,
    SEXP s_precision,
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    SEXP s_tolerance,
    SEXP s_n_bb,
    SEXP s_cri_probability,
    SEXP s_n_perms,
    SEXP s_blending_coef,
    SEXP s_verbose
    ) {

    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Extract scalar parameters using R's C API
    size_t min_path_size = (size_t)INTEGER(s_min_path_size)[0];
    size_t n_grid_vertices = (size_t)INTEGER(s_n_grid_vertices)[0];
    size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];
    size_t max_iterations = (size_t)INTEGER(s_max_iterations)[0];
    double precision = REAL(s_precision)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double tolerance = REAL(s_tolerance)[0];
    size_t n_bb = (size_t)INTEGER(s_n_bb)[0];
    double cri_probability = REAL(s_cri_probability)[0];
    size_t n_perms = (size_t)INTEGER(s_n_perms)[0];
    double blending_coef = REAL(s_blending_coef)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    size_t n_vertices = adj_list.size();
    uniform_grid_graph_t grid_graph;

    if (n_grid_vertices < n_vertices) {

        // The returned uniform_grid_graph_t object inherits both the graph_diameter
        // and max_packing_radius values from the intermediate set_wgraph_t object,
        // making these calculated values available for further analysis.
        grid_graph = create_maximal_packing(adj_list,
                                            weight_list,
                                            n_grid_vertices,
                                            max_iterations,
                                            precision);
        Rprintf("grid_graph.graph_diameter: %.4f\n", grid_graph.graph_diameter);

    } else {
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), 0);

        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);

        // Find diameter endpoints
        auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
        grid_graph.graph_diameter = diameter;
    }

    // Call the C++ function
    adaptive_uggmalo_result_t res = adaptive_uggmalo(
        grid_graph,
        y,
        min_path_size,
        n_bws,
        min_bw_factor,
        max_bw_factor,
        dist_normalization_factor,
        kernel_type,
        tolerance,
        n_bb,
        cri_probability,
        n_perms,
        blending_coef,
        verbose
    );

    // Create the return list using R's C API
    const char* names[] = {
        "graph_diameter",
        "grid_opt_bw",
        "predictions",
        "errors",
        "scale",
        "grid_predictions",
        "bb_predictions",
        "cri_lower",
        "cri_upper",
        "null_predictions",
        "null_predictions_cri_lower",
        "null_predictions_cri_upper",
        "p_values",
        "effect_sizes",
        "significant_vertices",
        "null_distances",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Create list and protect it
    SEXP result = PROTECT(allocVector(VECSXP, n_elements));
    SEXP result_names = PROTECT(allocVector(STRSXP, n_elements));

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, mkChar(names[i]));
    }
    setAttrib(result, R_NamesSymbol, result_names);

    // Helper function to convert vector to SEXP
    auto create_numeric_vector = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(REALSXP, vec.size()));
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Helper function to convert matrix to SEXP
    auto create_numeric_matrix = [](const std::vector<std::vector<double>>& mat) -> SEXP {
        if (mat.empty()) return R_NilValue;
        size_t nrow = mat.size();
        size_t ncol = mat[0].size();

        SEXP r_mat = PROTECT(allocMatrix(REALSXP, nrow, ncol));
        double* ptr = REAL(r_mat);

        for (size_t i = 0; i < nrow; ++i) {
            for (size_t j = 0; j < ncol; ++j) {
                ptr[i + j * nrow] = mat[i][j];  // Column-major order for R
            }
        }
        return r_mat;
    };

    // Helper function to convert boolean vector to SEXP
    auto create_logical_vector = [](const std::vector<bool>& vec) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(LGLSXP, vec.size()));
        int* ptr = LOGICAL(r_vec);
        for (size_t i = 0; i < vec.size(); ++i) {
            ptr[i] = vec[i];
        }
        return r_vec;
    };

    // Set graph diameter
    SEXP diam = PROTECT(allocVector(REALSXP, 1));
    REAL(diam)[0] = res.graph_diameter;
    SET_VECTOR_ELT(result, 0, diam);

    // Set other numeric vectors
    SET_VECTOR_ELT(result, 1, create_numeric_vector(res.grid_opt_bw));
    SET_VECTOR_ELT(result, 2, create_numeric_vector(res.predictions));
    SET_VECTOR_ELT(result, 3, create_numeric_vector(res.errors));
    SET_VECTOR_ELT(result, 4, create_numeric_vector(res.scale));
    SET_VECTOR_ELT(result, 5, create_numeric_vector(map_to_vector(res.grid_predictions_map)));

    // Set bootstrap fields if available
    if (!res.bb_predictions.empty()) {
        SET_VECTOR_ELT(result, 6, create_numeric_matrix(res.bb_predictions));
        SET_VECTOR_ELT(result, 7, create_numeric_vector(res.cri_lower));
        SET_VECTOR_ELT(result, 8, create_numeric_vector(res.cri_upper));
    } else {
        // Set to NULL if not available
        SET_VECTOR_ELT(result, 6, R_NilValue);
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
    }

    // Set permutation fields if available
    if (!res.null_predictions.empty()) {
        SET_VECTOR_ELT(result, 9, create_numeric_matrix(res.null_predictions));
        SET_VECTOR_ELT(result, 10, create_numeric_vector(res.null_predictions_cri_lower));
        SET_VECTOR_ELT(result, 11, create_numeric_vector(res.null_predictions_cri_upper));
    } else {
        SET_VECTOR_ELT(result, 9, R_NilValue);
        SET_VECTOR_ELT(result, 10, R_NilValue);
        SET_VECTOR_ELT(result, 11, R_NilValue);
    }

    // Set permutation test results if available
    if (res.permutation_tests) {
        SET_VECTOR_ELT(result, 12, create_numeric_vector(res.permutation_tests->p_values));
        SET_VECTOR_ELT(result, 13, create_numeric_vector(res.permutation_tests->effect_sizes));
        SET_VECTOR_ELT(result, 14, create_logical_vector(res.permutation_tests->significant_vertices));
        SET_VECTOR_ELT(result, 15, create_numeric_vector(res.permutation_tests->null_distances));
    } else {
        for (int i = 12; i <= 15; i++) {
            SET_VECTOR_ELT(result, i, R_NilValue);
        }
    }

    // Count protected objects: result, result_names, diam, and all non-NULL vectors/matrices
    int total_protected = 3;  // Start with result, result_names, and diam
    for (int i = 1; i < n_elements; i++) {
        if (VECTOR_ELT(result, i) != R_NilValue) total_protected++;
    }

    UNPROTECT(total_protected);
    return result;
}

