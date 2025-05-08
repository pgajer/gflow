#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>  // for std::reverse
#include <sstream>
#include <string>
#include <limits>
#include <map>
#include <utility>
#include <numeric> // for std::accumulate
#include <execution>
#include <atomic>
#include <mutex>
#include <random>     // for std::mt19937
#include <chrono>
#include <omp.h>

#include "edge_weights.hpp"
#include "ulm.hpp"
#include "graph_utils.hpp"
#include "uniform_grid_graph.hpp"
#include "centered_paths.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_shortest_path.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "kernels.h"
#include "sampling.h" // for C_runif_simplex()
#include "predictive_errors.hpp"
#include "wasserstein_perm_test.hpp"

extern "C" {
    SEXP S_uggmalo(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_best_models_coverage_factor,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_n_bws,
        SEXP s_grid_size,
        SEXP s_start_vertex,
        SEXP s_snap_tolerance,
        SEXP s_dist_normalization_factor,
        SEXP s_min_path_size,
        SEXP s_diff_threshold,
        SEXP s_kernel_type,
        SEXP s_fit_quadratic,
        SEXP s_tolerance,
        SEXP s_n_bb,
        SEXP s_p,
        SEXP s_n_perms,
        SEXP s_verbose
        );
}

/**
 * @brief Implements the Uniform Grid Graph Model-Averaged LOcal Linear regression (UGGMALO) algorithm
 *
 * @details This function implements a sophisticated algorithm for analyzing weighted graphs using
 * local path linear models with model averaging. The algorithm performs the following main steps:
 * 1. Computes graph diameter and determines bandwidth range
 * 2. Creates a uniform grid representation of the input graph
 * 3. For each candidate bandwidth:
 *    - Processes paths through grid vertices
 *    - Fits local linear models to path data
 *    - Computes weighted predictions and errors
 * 4. Determines optimal bandwidth based on cross-validation errors
 *
 * The algorithm uses weighted linear regression on paths through the graph to create
 * local models, which are then combined using weighted averaging. Model evaluation
 * is performed using leave-one-out cross-validation with Brier score errors.
 *
 * @details Algorithm Performance Characteristics:
 * - Time Complexity: O(V * P * M) where:
 *   V = number of grid vertices
 *   P = average number of paths per vertex
 *   M = maximum iterations in linear regression
 * - Space Complexity: O(V * N) where N = number of original vertices
 *
 * @details Error Handling:
 * - Non-convergent fits are skipped with warnings
 * - Invalid predictions default to original y values
 * - Infinite errors indicate no valid predictions
 *
 * @details Model Selection Strategy:
 * - Models are ranked by mean LOOCV Brier score error
 * - Selection continues until coverage threshold is met
 * - Minimum path size constraint ensures model stability
 *
 * @param adj_list Vector of vectors containing adjacency lists for each vertex
 * @param weight_list Vector of vectors containing edge weights corresponding to adj_list
 * @param y Vector of response values (0 or 1) for original vertices
 * @param best_models_coverage_factor Proportion of the number of vertices of the given neighborhood (given by the value of a bandwidth) of the given grid vertex that the models with the smallest mean error have to cover to stop the best models selection process.
 * @param min_bw_factor Factor for minimum bandwidth (must be > 0)
 * @param max_bw_factor Factor for maximum bandwidth (must be < 1)
 * @param n_bws Number of bandwidth values to test
 * @param grid_size Number of grid points to create in uniform grid graph
 * @param start_vertex Starting vertex for grid construction
 * @param snap_tolerance Tolerance for snapping vertices to grid positions
 * @param dist_normalization_factor Factor for normalizing distances in kernel computations
 * @param min_path_size Minimum required path length in original vertices
 * @param diff_threshold Threshold for determining path direction differences
 * @param kernel_type Type of kernel function to use (1-10)
 * @param fit_quadratic Whether to include quadratic term in linear regression (default: false)
 * @param tolerance Convergence tolerance for linear regression (default: 1e-8)
 * @param with_errors Whether to compute LOOCV errors (default: true)
 * @param verbose Enable progress reporting and warnings (default: false)
 *
 * @return uggmalo_t structure containing:
 *         - candidate_bws: Vector of tested bandwidth values
 *         - bw_predictions: Predictions for each bandwidth
 *         - mean_errors: Mean errors for each bandwidth
 *         - opt_bw: Optimal bandwidth value
 *         - opt_predictions: Predictions for optimal bandwidth
 *         - opt_bw_index: Index of optimal bandwidth
 *         - graph_diameter: Computed diameter of the input graph
 *
 * @note Memory Usage: Scales primarily with grid_size and number of paths
 * @note Thread Safety: Function is not thread-safe due to shared state
 * @note Performance: Consider reducing grid_size for large graphs
 * @note Input validation is performed in the parent R function
 * @note Progress is reported every 5 bandwidths when verbose is true
 * @note Vertices with no valid predictions retain their original y values
 * @note Model averaging uses kernel-weighted means of predictions and errors
 * @note The algorithm handles non-convergent model fits by skipping those paths
 *
 * @warning May return infinite mean error if no valid predictions for a bandwidth
 * @warning Linear regression may not converge for some paths
 *
 * @pre All input vectors must not be empty
 * @pre best_models_coverage_factor must be in (0,1]
 * @pre weight_list values must be non-negative
 * @pre y values must be 0 or 1 only
 * @pre All vectors in adj_list and weight_list must have corresponding sizes
 * @pre y must have length equal to the number of vertices
 * @pre 0 < min_bw_factor < max_bw_factor < 1
 * @pre n_bws > 0
 * @pre grid_size > 0
 * @pre start_vertex must be a valid vertex index
 * @pre snap_tolerance > 0
 * @pre dist_normalization_factor > 0
 * @pre min_path_size > 0
 * @pre diff_threshold > 0
 * @pre 1 <= kernel_type <= 10
 *
 * @see create_uniform_grid_graph()
 * @see ugg_get_path_data()
 * @see eigen_ulogit_fit()
 */
struct uggmalo_t {
    std::vector<double> candidate_bws;
    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> mean_errors;
    double opt_bw;
    size_t opt_bw_index;
    std::vector<double> opt_predictions;
    double graph_diameter; // relative to start_vertex; ideally we should allow start_vertex to be undefined and then the true graph diameter would be computed; in some context this can be calculated by looking only at the distances between vertices of degree 1 or even better the evenness end vertices

    // Bootstrap fields
    std::vector<std::vector<double>> bb_predictions;  // All bootstrap predictions [n_bb][n_points]
    std::vector<double> bb_central;                   // Central location estimates
    std::vector<double> cri_lower;                    // Lower credible interval bounds
    std::vector<double> cri_upper;                    // Upper credible interval bounds
    // Permutation test fields
    std::vector<std::vector<double>> null_predictions;  // Predictions for each permutation [n_perms][n_points]
    // Vertex-level test results (only populated if both n_bb > 0 and n_perms > 0)
    std::optional<vertex_wasserstein_perm_test_results_t> permutation_tests; // New field
};

uggmalo_t uggmalo(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    double best_models_coverage_factor,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    size_t grid_size,
    size_t start_vertex,
    double snap_tolerance,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    bool fit_quadratic,
    double tolerance,
    size_t n_bb,          // Number of Bayesian bootstrap iterations
    double p,          // Probability level for credible intervals
    size_t n_perms,
    bool verbose
    ) {

    const bool with_errors = true; // we want path linear models to always return LOOCV error estimate

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Validation of input parameters will be done in the parent R function
    uggmalo_t result;

    // Step 1: Compute graph diameter and bandwidth range
    result.graph_diameter = get_vertex_eccentricity(adj_list, weight_list, start_vertex);

    // debugging
    Rprintf("In uggmalo()\n");
    Rprintf("start_vertex: %zu\n", start_vertex);

    // Generate candidate bandwidths
    result.candidate_bws.resize(n_bws);
    double min_bw = min_bw_factor * result.graph_diameter;
    double max_bw = max_bw_factor * result.graph_diameter;

    if (min_bw <= 0.0 || std::isnan(min_bw) || std::isinf(min_bw)) {
        REPORT_ERROR("Invalid minimum bandwidth calculated. Check graph_diameter and min_bw_factor. min_bw: %f\n", min_bw);
    }

    if (n_bws == 1) {
        result.candidate_bws[0] = min_bw;
    } else {
        double bw_step = (max_bw - min_bw) / (n_bws - 1);

        // Defining minimum acceptable step size relative to the bandwidth range
        double min_step_threshold = (max_bw - min_bw) * 1e-6;

        if (bw_step < min_step_threshold) {
            // Adjusting n_bws automatically
            size_t adjusted_n_bws = static_cast<size_t>((max_bw - min_bw) / min_step_threshold) + 1;
            if (verbose) {
                REPORT_WARNING("Adjusting n_bws from %zu to %zu to ensure adequate step size\n",
                               n_bws, adjusted_n_bws);
            }
            n_bws = adjusted_n_bws;
            bw_step = (max_bw - min_bw) / (n_bws - 1);
            result.candidate_bws.resize(n_bws);
        }

        for(size_t i = 0; i < n_bws; i++) {
            result.candidate_bws[i] = min_bw + i * bw_step;
        }
    }

    Rprintf("min_bw_factor: %f\n", min_bw_factor);
    Rprintf("max_bw_factor: %f\n", max_bw_factor);
    Rprintf("graph_diameter: %f\n", result.graph_diameter);
    Rprintf("min_bw: %f\n", min_bw);
    Rprintf("max_bw: %f\n", max_bw);
    Rprintf("n_bws: %zu\n", n_bws);
    print_vect(result.candidate_bws, "result.candidate_bws");

    // Step 2: Create uniform grid graph
    uniform_grid_graph_t uniform_grid_graph = create_uniform_grid_graph(
        adj_list,
        weight_list,
        grid_size,
        start_vertex,
        snap_tolerance
        );

    // debugging
    // std::vector<double> distances_from_start_vertex_to_grid_vertices = uniform_grid_graph.compute_shortest_path_distances(start_vertex);
    // Rprintf("\n\n");
    // uniform_grid_graph.print_grid_vertices();
    // print_vect(distances_from_start_vertex_to_grid_vertices, "distances_from_start_vertex_to_grid_vertices");

    // print_vect_vect(adj_list, "adj_list");
    // print_vect_vect(weight_list, "weight_list");
    // uniform_grid_graph.print("uniform_grid_graph");
    // Rprintf("\n\n");
    // error("uggmalo() debugging");


    edge_weights_t edge_weights = precompute_edge_weights(adj_list, weight_list);

    // Initialize storage for predictions and errors
    const size_t n_vertices = uniform_grid_graph.n_original_vertices;
    result.bw_predictions.resize(n_bws);
    for (auto& predictions : result.bw_predictions) {
        predictions.resize(n_vertices);
    }
    result.mean_errors.resize(n_bws);

    // Creating an extended eigen_ulogit_t structure that contains besides the eigen_ulogit_t components also: mean_error, vertices and w_path components. The models will be sorted w/r to mean_error
    struct ext_ulm_t : public ulm_t {
        double mean_error;            // mean LOOCV error of the model
        std::vector<size_t> vertices; // vertices of the path along which the model is estimated
        std::vector<double> w_path;   // weights at the vertices of the path

        // Add a constructor that takes the base class
        explicit ext_ulm_t(const ulm_t& base)
            : ulm_t(base),      // Initialize base class members
              mean_error(0.0),  // Initialize new members
              vertices(),       // Empty vector
              w_path()          // Empty vector
            {}

        // This operator enables automatic sorting in set/multiset in the ascending order from the smallest to largest model's mean_error - models with the smallest mean_error are the most desirable
        bool operator<(const ext_ulm_t& other) const {
            return mean_error < other.mean_error;
        }
    };

    // weight/prediction/error struct needed for mode averaging; we use it to record weight/prediction/error of the given vertex in each model where the vertext is in the support of the model
    struct wpe_t {
        double weight;
        double prediction;
        double error;

        // Constructor needed for emplace_back(x,y,z)
        wpe_t(double w, double p, double e)
            : weight(w), prediction(p), error(e) {}
    };

    // Precomputing paths for each grid vertex
    std::unordered_map<size_t, reachability_map_t> vertex_paths;
    double max_bandwidth = result.candidate_bws.back();  // Largest bandwidth

    for (const auto& grid_vertex : uniform_grid_graph.grid_vertices) {
        vertex_paths[grid_vertex] = precompute_max_bandwidth_paths(
            uniform_grid_graph,
            grid_vertex,
            max_bandwidth
            );
    }

    // Lambda function to process a single bandwidth value
    auto process_bw = [
        // Algorithm parameters
        &edge_weights,
        &vertex_paths,
        best_models_coverage_factor,
        dist_normalization_factor,
        min_path_size,
        diff_threshold,
        kernel_type,
        // Model fitting parameters
        y_binary,
        fit_quadratic,
        tolerance,
        with_errors,
        // Debug parameter
        &verbose
        ](
            double bandwidth,
            const uniform_grid_graph_t& uniform_grid_graph,
            const std::vector<double>& y,
            std::vector<double>& predictions,
            double& mean_error,
            std::optional<std::vector<double>> weights
            ) {

        if (weights && (*weights).size() != y.size()) {
            REPORT_ERROR("weights vector size does not match y size");
        }

        // This lambda processes a single bandwidth value to:
        // 1. Fit local linear models for paths through each grid vertex
        // 2. Select best models based on error and coverage criteria
        // 3. Compute weighted average predictions and errors
        // Returns: Updates predictions vector and mean_error for this bandwidth

        // Get the number of vertices in the original graph
        const size_t n_vertices = uniform_grid_graph.n_original_vertices;
        std::vector<std::vector<wpe_t>> wpe(n_vertices); // wpe[i] stores a vector of {weight, prediction, error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions

        // Phase 1: fit linear regression models to the path data associated with each grid vertex
        for (const auto& grid_vertex : uniform_grid_graph.grid_vertices) {
            // Each grid vertex serves as a reference point for finding paths
            // and fitting local linear models

            // debugging
            fprintf(stderr, "\nIn uggmalo process_bw() processing grid_vertex: %zu\n", grid_vertex);

            // Collect paths that pass through this grid vertex
            std::vector<path_data_t> paths = ugg_get_path_data_efficient(
                uniform_grid_graph,
                y,
                grid_vertex,
                bandwidth,
                vertex_paths[grid_vertex],
                dist_normalization_factor,
                min_path_size,
                diff_threshold,
                kernel_type,
                edge_weights
                );

            if (paths.empty()) {
                if (verbose) {
                    REPORT_WARNING("No valid paths found through vertex %d at bandwidth %f\n",
                                   grid_vertex, bandwidth);
                }
                continue;
            }

            fprintf(stderr, "After ugg_get_path_data_efficient()\nNumber path_data_t objects found: %zu\n", paths.size());


            // Fitting linear models to valid paths through the current grid vertex
            //
            // From these models we are going to select the smallest number of the best
            // models so that the union of their supports covers the set of original
            // vertices of the given neighborhood (given by the size of the given bandwidth)

            // Identifying the support set of all paths
            std::unordered_set<int> all_vertices;
            for (const auto& path : paths) {
                all_vertices.insert(path.vertices.begin(), path.vertices.end());
            }

            // Creating a container for models so that they can be sorted w/r mean_error
            std::multiset<ext_ulm_t> all_models;
            for (auto& path : paths) {
                // Fit weighted linear model to the current path

                if (weights) {
                    // multiply path.w_path by the corresponding values of weights
                    for (size_t i = 0; i < path.vertices.size(); i++) {
                        path.w_path[i] *= (*weights)[path.vertices[i]];
                    }
                }

                ulm_t fit_result = ulm(
                    path.x_path.data(),
                    path.y_path.data(),
                    path.w_path,
                    y_binary,
                    tolerance
                    );

                // Create the extended object from fit_result
                ext_ulm_t extended_result(fit_result);

                // Calculate the mean of loocv_brier_errors vector
                // We'll use std::accumulate to sum all elements and divide by the vector size
                double sum = std::accumulate(fit_result.errors.begin(),
                                             fit_result.errors.end(),
                                             0.0);  // Use 0.0 to ensure floating-point arithmetic

                extended_result.mean_error = sum / fit_result.errors.size();
                extended_result.vertices   = std::move(path.vertices);
                extended_result.w_path     = std::move(path.w_path);

                all_models.insert(extended_result);

            }

            // Iterating through all_models to select the minimal number of the
            // models with the smallest mean error that cover certain fraction
            // of the set of all vertices
            size_t coverage_thld = static_cast<size_t>(best_models_coverage_factor * all_vertices.size());
            if (coverage_thld < min_path_size) {
                coverage_thld = min_path_size;
            }
            if (coverage_thld > all_vertices.size()) {
                coverage_thld = all_vertices.size();
            }
            std::unordered_set<int> vertices_of_selected_models;
            for (const auto& model : all_models) {

                vertices_of_selected_models.insert(model.vertices.begin(), model.vertices.end());

                // Store weight, prediction and error for the given vertex from the given 'model' in wpe[ model.vertices[i] ]
                for (size_t i = 0; i < model.vertices.size(); ++i) {
                    wpe[ model.vertices[i] ].emplace_back(model.w_path[i], model.predictions[i], model.errors[i]);
                }

                if (vertices_of_selected_models.size() > coverage_thld) {
                    break;
                }
            }
        } // END OF for (const auto& grid_vertex : uniform_grid_graph.grid_vertices)

        // Phase 2: Model averaging
        double total_error = 0.0;
        size_t valid_vertex_count = 0;
        for (size_t i = 0; i < n_vertices; i++ ) {

            double prediction_sum = 0.0;
            double weight_sum   = 0.0;
            double error_sum  = 0.0;

            for (const auto& x : wpe[i]) {
                prediction_sum += x.weight * x.prediction;
                weight_sum += x.weight;
                error_sum +=  x.weight * x.error;
            }

            if (weight_sum > 0) {

                if (weight_sum < 1e-10) {
                    if (verbose) {
                        REPORT_WARNING("Very small weight sum encountered for vertex %zu. Setting predictions[i] to y[i]\n", i);
                    }
                    predictions[i] = std::numeric_limits<double>::quiet_NaN();
                } else {
                    total_error   += error_sum / weight_sum;
                    predictions[i] = prediction_sum / weight_sum;
                    valid_vertex_count++;
                }
            } else {
                predictions[i] = y[i];
            }
        }

        mean_error = valid_vertex_count > 0 ? total_error / valid_vertex_count :
            std::numeric_limits<double>::infinity();
    };

    // Step 3: Process each bandwidth in parallel
    std::vector<int> bw_indices(n_bws);
    std::iota(bw_indices.begin(), bw_indices.end(), 0);

    // Progress tracking
    std::atomic<int> bandwidth_counter{0};
    const size_t progress_chunk = std::max<size_t>(1, n_bws / 10);  // Report every 10% progress

    // Parallel execution of bandwidth processing
    //auto exec = std::execution::par_unseq,
    auto exec = std::execution::seq;
    std::for_each(exec,
                  bw_indices.begin(),
                  bw_indices.end(),
                  [&](int bw_idx) {
                      if (verbose && (bw_idx % 5 == 0)) {
                          // Thread-safe progress update
                          int current_count = ++bandwidth_counter;
                          if (current_count % progress_chunk == 0) {
                              REprintf("\rProcessing bandwidth %d%%",
                                       static_cast<int>((100.0 * current_count) / n_bws));
                          }
                      }

                      process_bw(
                          result.candidate_bws[bw_idx],
                          uniform_grid_graph,
                          y,
                          result.bw_predictions[bw_idx],
                          result.mean_errors[bw_idx],
                          std::nullopt
                          );

                      if (result.mean_errors[bw_idx] == std::numeric_limits<double>::infinity()) {
                          if (verbose) {
                              REprintf("Warning: No valid mean errors for bandwidth index %d and bandwidth %f\n",
                                       bw_idx, result.candidate_bws[bw_idx]);
                          }
                      }
                  });

    if (verbose) {
        Rprintf("\nBandwidth processing completed.\n");
    }

    // Step 4: Find optimal bandwidth
    if (std::all_of(result.mean_errors.begin(), result.mean_errors.end(),
                    [](double error) { return std::isinf(error); })) {
        if (verbose) {
            REPORT_ERROR("All bandwidths resulted in infinite errors");
        }
        // result.opt_bw_index = n_bws / 2;
    } else {
        result.opt_bw_index = std::distance(
            result.mean_errors.begin(),
            std::min_element(result.mean_errors.begin(), result.mean_errors.end())
            );
        result.opt_bw = result.candidate_bws[result.opt_bw_index];
        result.opt_predictions = result.bw_predictions[result.opt_bw_index];
    }

    // Step 5: Generate Bayesian bootstraps
    if (n_bb > 0) {
        if (p <= 0.0 || p >= 1.0) {
            REPORT_ERROR("Probability level p must be in (0,1)");
        }

        const int n_points = y.size();

        // Initialize results vectors
        result.bb_predictions.resize(n_bb);
        for (auto& pred : result.bb_predictions) {
            pred.resize(n_points);
        }

        // Create thread-local random number generators
        const unsigned int num_threads = std::thread::hardware_concurrency();
        std::vector<std::mt19937> thread_rngs(num_threads);
        for (unsigned int i = 0; i < num_threads; ++i) {
            std::random_device rd;
            thread_rngs[i].seed(rd());
        }

        // Create indices for parallel iteration
        std::vector<int> bb_indices(n_bb);
        std::iota(bb_indices.begin(), bb_indices.end(), 0);

        // Progress tracking
        std::atomic<int> bootstrap_counter{0};
        const size_t progress_chunk = std::max<size_t>(1, n_bb / 100);  // Report every 1% progress
        //auto start_time = std::chrono::steady_clock::now();

        // Parallel execution with thread-local RNG
        // std::execution::par_unseq,
        std::for_each(exec,
                      bb_indices.begin(),
                      bb_indices.end(),
                      [&](int iboot) {
                          // Get thread-local RNG
                          const int thread_id = omp_get_thread_num() % num_threads;
                          auto& local_rng = thread_rngs[thread_id];

                          // Generate weights using thread-local RNG
                          std::vector<double> weights(n_points);
                          {
                              // Use thread-local RNG to generate weights
                              std::gamma_distribution<double> gamma(1.0, 1.0);
                              double sum = 0.0;
                              for (int i = 0; i < n_points; ++i) {
                                  weights[i] = gamma(local_rng);
                                  sum += weights[i];
                              }
                              // Normalize weights
                              for (int i = 0; i < n_points; ++i) {
                                  weights[i] /= sum;
                              }
                          }

                          // Process bootstrap with thread-local weights
                          double bootstrap_error;
                          process_bw(
                              result.opt_bw,
                              uniform_grid_graph,
                              y,
                              result.bb_predictions[iboot],
                              bootstrap_error,
                              std::optional<std::vector<double>>(weights)
                              );

                          // Thread-safe progress update without mutex
                          if (verbose) {
                              int current_count = ++bootstrap_counter;
                              if (current_count % progress_chunk == 0) {
                                  // Use \r to move cursor to start of line
                                  REprintf("\rBootstrap progress: %d%%",
                                           static_cast<int>((100.0 * current_count) / n_bb));
                              }
                          }

                          #if 0
                          // Thread-safe progress update
                          int current_count = ++bootstrap_counter;
                          if (verbose && (current_count % progress_chunk == 0)) {
                              auto current_time = std::chrono::steady_clock::now();
                              auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                                  current_time - start_time).count();
                              double progress = (100.0 * current_count) / n_bb;

                              std::lock_guard<std::mutex> lock(console_mutex); // line (uggmalo.cpp:704) <<---
                              Rprintf("\rBootstrap progress: %.1f%% (%d/%zu) - Elapsed: %llds",
                                      progress, current_count, n_bb, elapsed);
                              //R_FlushConsole();
                          }
                          #endif
                      });


        bool use_median = true;
        bb_cri_t bb_cri_res = bb_cri(result.bb_predictions, use_median, p);

        result.bb_central = std::move(bb_cri_res.bb_Ey);
        result.cri_lower  = std::move(bb_cri_res.cri_L);
        result.cri_upper  = std::move(bb_cri_res.cri_U);

        if (verbose) {
            Rprintf("\nBootstrap completed.\n");
        }
    }

    // Step 6: Permutation Test
    if (n_perms > 0) {
        const size_t n_points = y.size();

        // Initialize permutation test results
        result.null_predictions.resize(n_perms);
        for (auto& pred : result.null_predictions) {
            pred.resize(n_points);
        }

        // Create indices for parallel iteration
        std::vector<int> perm_indices(n_perms);
        std::iota(perm_indices.begin(), perm_indices.end(), 0);

        // Create indices for permutation
        std::vector<size_t> indices(n_points);
        std::iota(indices.begin(), indices.end(), 0);

        // Mutex for thread-safe random number generation
        std::mutex rng_mutex;

        // Atomic counter for tracking progress
        std::atomic<int> permutation_counter{0};
        const size_t progress_chunk = std::max<size_t>(1, n_perms / 100);  // Report every 1% progress

        // Track time for progress updates
        // auto ptm = std::chrono::steady_clock::now();

        // Create a random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        // Parallel execution of permutation iterations
        // std::execution::par_unseq,
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
                          std::vector<double> y_shuffled(n_points);
                          for (size_t i = 0; i < n_points; ++i) {
                              y_shuffled[i] = y[indices[i]];
                          }

                          // Process permuted data using optimal bandwidth
                          double dummy_error;
                          process_bw(
                              result.opt_bw,
                              uniform_grid_graph,
                              y_shuffled,
                              result.null_predictions[iperm],
                              dummy_error,    // Dummy variable for unused error
                              std::nullopt   // No weights for permutation test
                              );

                          // Increment counter and update progress
                          if (verbose) {
                              int current_count = ++permutation_counter;
                              if (current_count % progress_chunk == 0) {
                                  // Use \r to move cursor to start of line
                                  REprintf("\rPermutation Test Progress: %d%%",
                                           static_cast<int>((100.0 * current_count) / n_perms));
                              }
                          }

                          #if 0
                          if (verbose) {
                              int current_count = ++permutation_counter;
                              {
                                  std::lock_guard<std::mutex> lock(console_mutex);
                                  char msg[100];
                                  snprintf(msg, sizeof(msg), "\rProcessed %d permutations", current_count);
                                  elapsed_time(ptm, msg, true);
                                  ptm = std::chrono::steady_clock::now();
                              }
                          }
                          #endif
                      });

        if (verbose) {
            Rprintf("\nPermutation Test Completed.\n");
        }
    }

    if (n_bb > 0 && n_perms > 0) {
        size_t n_bootstraps = 1000;
        double alpha = 0.05;
        result.permutation_tests = vertex_wasserstein_perm_test(
            result.bb_predictions,
            result.null_predictions,
            n_bootstraps,
            alpha
            );
    }

    return result;
}

/**
 * @brief R interface for the Uniform Grid Graph Model-Averaged LOcal Linear regression (UGGMALO) algorithm
 *
 * @details This function serves as the interface between R and the C++ implementation of the UGGMALO
 * algorithm. It converts R objects to C++ types, calls the core algorithm, and converts the results
 * back to R objects.
 *
 * @param s_adj_list SEXP (list) Adjacency list representation of the graph
 * @param s_weight_list SEXP (list) Edge weights corresponding to adj_list
 * @param s_y SEXP (numeric) Response values (0 or 1) for vertices
 * @param s_best_models_coverage_factor (numeric) Proportion of the set of given neighborhood of a grid vertex that the models with the smallest mean error have to cover to stop the best models selection process.
 * @param s_min_bw_factor SEXP (numeric) Minimum bandwidth factor
 * @param s_max_bw_factor SEXP (numeric) Maximum bandwidth factor
 * @param s_n_bws SEXP (integer) Number of bandwidth values to test
 * @param s_grid_size SEXP (integer) Number of grid points
 * @param s_start_vertex SEXP (integer) Starting vertex index (1-based in R)
 * @param s_snap_tolerance SEXP (numeric) Tolerance for grid point snapping
 * @param s_dist_normalization_factor SEXP (numeric) Distance normalization factor
 * @param s_min_path_size SEXP (integer) Minimum path size
 * @param s_diff_threshold SEXP (integer) Path direction difference threshold
 * @param s_kernel_type SEXP (integer) Kernel function type (1-10)
 * @param s_fit_quadratic SEXP (logical) Whether to fit quadratic term
 * @param s_tolerance SEXP (numeric) Convergence tolerance
 * @param s_verbose SEXP (logical) Enable verbose output
 *
 * @return SEXP (list) with components:
 *   - x_grid: Numeric vector of grid points
 *   - predictions: Numeric vector of predictions at original vertices
 *   - bw_grid_predictions: Matrix of predictions for each bandwidth
 *   - bw_grid_errors: Matrix of errors for each bandwidth
 *   - mean_brier_errors: Vector of mean Brier errors
 *   - opt_brier_bw_idx: Optimal bandwidth index (1-based)
 *   - bws: Vector of candidate bandwidths
 *   - fit_info: List containing fitting parameters:
 *     * fit_quadratic: Logical
 *     * pilot_bandwidth: Numeric
 *     * kernel_type: Integer
 *     * cv_folds: Integer
 *     * min_bw_factor: Numeric
 *     * max_bw_factor: Numeric
 *     * tolerance: Numeric
 *
 * @note All R indices are 1-based and converted to 0-based for C++
 * @note The function handles memory protection for all SEXP objects
 * @note Returns errors to R using the error() mechanism
 */
SEXP S_uggmalo(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_best_models_coverage_factor,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_grid_size,
    SEXP s_start_vertex,
    SEXP s_snap_tolerance,
    SEXP s_dist_normalization_factor,
    SEXP s_min_path_size,
    SEXP s_diff_threshold,
    SEXP s_kernel_type,
    SEXP s_fit_quadratic,
    SEXP s_tolerance,
    SEXP s_n_bb,
    SEXP s_p,
    SEXP s_n_perms,
    SEXP s_verbose
    ) {

    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    size_t n_vertices = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_vertices);

    double best_models_coverage_factor = REAL(s_best_models_coverage_factor)[0];
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];
    size_t n_bws = INTEGER(s_n_bws)[0];
    size_t grid_size = INTEGER(s_grid_size)[0];
    size_t start_vertex = INTEGER(s_start_vertex)[0];
    double snap_tolerance = REAL(s_snap_tolerance)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    size_t min_path_size = INTEGER(s_min_path_size)[0];
    size_t diff_threshold = INTEGER(s_diff_threshold)[0];
    size_t kernel_type = INTEGER(s_kernel_type)[0];
    bool fit_quadratic = LOGICAL(s_fit_quadratic)[0];
    double tolerance = REAL(s_tolerance)[0];
    size_t n_bb = INTEGER(s_n_bb)[0];
    double p = REAL(s_p)[0];
    size_t n_perms = INTEGER(s_n_perms)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    uggmalo_t res = uggmalo(
        adj_list,
        weight_list,
        y,
        best_models_coverage_factor,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        grid_size,
        start_vertex,
        snap_tolerance,
        dist_normalization_factor,
        min_path_size,
        diff_threshold,
        kernel_type,
        fit_quadratic,
        tolerance,
        n_bb,
        p,
        n_perms,
        verbose
        );


    // Create the return list with the CORRECT number of elements
    const char* names[] = {
        "candidate_bws",
        "bw_predictions",
        "mean_errors",
        "opt_bw",
        "opt_bw_idx",
        "opt_predictions",
        "graph_diameter",
        "bb_predictions",
        "bb_central",
        "cri_lower",
        "cri_upper",
        "null_predictions",
        "p_values",
        "effect_sizes",
        "significant_vertices",
        NULL  // Add NULL terminator
    };

    // Convert results to R list
    size_t n_protected = 0;
    size_t n_elements = sizeof(names)/sizeof(names[0]) - 1;  // subtract 1 for NULL terminator
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements)); n_protected++;

    // Create and set the names
    SEXP names_vec = PROTECT(Rf_allocVector(STRSXP, n_elements)); n_protected++;
    for(size_t i = 0; i < n_elements; i++) {
        SET_STRING_ELT(names_vec, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, names_vec);

    // Helper function to convert std::vector<double> to SEXP
    auto vec_to_sexp = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        std::copy(vec.begin(), vec.end(), REAL(r_vec));
        UNPROTECT(1);
        return r_vec;
    };

    // Helper function to convert std::vector<std::vector<double>> to SEXP
    auto matrix_to_sexp = [](const std::vector<std::vector<double>>& mat) -> SEXP {
        if (mat.empty()) return R_NilValue;
        size_t nrow = mat.size();
        size_t ncol = mat[0].size();
        SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, ncol, nrow));
        double* ptr = REAL(r_mat);
        for (size_t i = 0; i < nrow; ++i) {
            for (size_t j = 0; j < ncol; ++j) {
                ptr[j + i * ncol] = mat[i][j];
            }
        }
        UNPROTECT(1);
        return r_mat;
    };

    // Helper function to convert std::vector<bool> to SEXP
    auto bool_vec_to_sexp = [](const std::vector<bool>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(LGLSXP, vec.size()));
        for (size_t i = 0; i < vec.size(); ++i) {
            LOGICAL(r_vec)[i] = vec[i];
        }
        UNPROTECT(1);
        return r_vec;
    };

    // Initialize all elements to R_NilValue first
    for(size_t i = 0; i < n_elements; i++) {
        SET_VECTOR_ELT(result, i, R_NilValue);
    }

    // Populate existing fields
    SET_VECTOR_ELT(result, 0, vec_to_sexp(res.candidate_bws));
    SET_VECTOR_ELT(result, 1, matrix_to_sexp(res.bw_predictions));
    SET_VECTOR_ELT(result, 2, vec_to_sexp(res.mean_errors));
    SEXP opt_bw = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
    REAL(opt_bw)[0] = res.opt_bw;
    SET_VECTOR_ELT(result, 3, opt_bw);
    SEXP opt_idx = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(opt_idx)[0] = res.opt_bw_index + 1;  // Convert to 1-based indexing for R
    SET_VECTOR_ELT(result, 4, opt_idx);
    SET_VECTOR_ELT(result, 5, vec_to_sexp(res.opt_predictions));
    SEXP diam = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
    REAL(diam)[0] = res.graph_diameter;
    SET_VECTOR_ELT(result, 6, diam);

    // Populate bootstrap fields (if available)
    if (!res.bb_predictions.empty()) {
        SET_VECTOR_ELT(result, 7, matrix_to_sexp(res.bb_predictions));
        SET_VECTOR_ELT(result, 8, vec_to_sexp(res.bb_central));
        SET_VECTOR_ELT(result, 9, vec_to_sexp(res.cri_lower));
        SET_VECTOR_ELT(result, 10, vec_to_sexp(res.cri_upper));
    }

    // Populate permutation fields (if available)
    if (!res.null_predictions.empty()) {
        SET_VECTOR_ELT(result, 11, matrix_to_sexp(res.null_predictions));
    }

    // Populate vertex test results (if available)
    if (res.permutation_tests) {
        SET_VECTOR_ELT(result, 12, vec_to_sexp(res.permutation_tests->p_values));
        SET_VECTOR_ELT(result, 13, vec_to_sexp(res.permutation_tests->effect_sizes));
        SET_VECTOR_ELT(result, 14, bool_vec_to_sexp(res.permutation_tests->significant_vertices));
    }

    UNPROTECT(n_protected);
    return result;
}
