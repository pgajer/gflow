#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

// #include <omp.h>
#include "omp_compat.h"

// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <cmath>
#include <limits>
// #include <iostream>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <map>
#include <numeric>
#include <random>
#include <chrono>

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "Eigen_utils.h"
#include "msr2.h"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_diffusion_smoother.hpp"
#include "kernels.h"
#include "error_utils.h"

void scale_to_range(std::vector<double>& x, double ymin, double ymax);
double calculate_smallest_difference(const std::vector<double>& y);
std::unique_ptr<std::pair<std::vector<int>, std::vector<int>>> find_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y);
std::unique_ptr<std::vector<double>> prop_nbhrs_with_smaller_y(const std::vector<std::vector<int>>& graph,
                                                               const std::vector<double>& y );
std::unique_ptr<std::vector<double>> cv_imputation(const std::set<int>& test_set,
                                                   const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& d,
                                                   const std::vector<double>& y,
                                                   bool y_binary,
                                                   imputation_method_t imputation_method,
                                                   iterative_imputation_params_t iterative_params,
                                                   bool apply_binary_threshold,
                                                   double binary_threshold,
                                                   int ikernel,
                                                   double dist_normalization_factor);

extern "C" {
    SEXP S_graph_diffusion_smoother(SEXP s_adj_list,
                                    SEXP s_weight_list,
                                    SEXP s_y,
                                    SEXP s_n_time_steps,
                                    SEXP s_step_factor,
                                    SEXP s_binary_threshold,
                                    SEXP s_ikernel,
                                    SEXP s_dist_normalization_factor,
                                    SEXP s_n_CVs,
                                    SEXP s_n_CV_folds,
                                    SEXP s_epsilon,
                                    SEXP s_verbose,
                                    SEXP s_seed);

    SEXP S_ext_graph_diffusion_smoother(SEXP Rgraph,
                                    SEXP Rd,
                                    SEXP Rweights,
                                    SEXP Ry,
                                    SEXP Rn_time_steps,
                                    SEXP Rstep_factor,
                                    SEXP Rnormalize,
                                    SEXP Rpreserve_local_maxima,
                                    SEXP Rlocal_maximum_weight_factor,
                                    SEXP Rpreserve_local_extrema,
                                    SEXP Rimputation_method,
                                    SEXP Rmax_iterations,
                                    SEXP Rconvergence_threshold,
                                    SEXP Rapply_binary_threshold,
                                    SEXP Rbinary_threshold,
                                    SEXP Rikernel,
                                    SEXP Rdist_normalization_factor,
                                    SEXP Rn_CVs,
                                    SEXP Rn_CV_folds,
                                    SEXP Repsilon,
                                    SEXP Rseed,
                                    SEXP Rn_cores,
                                    SEXP Rverbose);

    SEXP S_instrumented_gds(SEXP s_graph,
                            SEXP s_edge_lengths,
                            SEXP s_y,
                            SEXP s_y_true,
                            SEXP s_n_time_steps,
                            SEXP s_base_step_factor,
                            SEXP s_use_pure_laplacian,
                            SEXP s_ikernel,
                            SEXP s_kernel_scale,
                            SEXP s_increase_factor,
                            SEXP s_decrease_factor,
                            SEXP s_oscillation_factor,
                            SEXP s_min_step,
                            SEXP s_max_step);
}

// version parallelizing cross-validation part
/**
 * @brief Applies a graph diffusion smoothing algorithm with cross-validation, supporting parallelization.
 *
 * This function smooths a function `y` defined over the vertices of a graph using
 * a diffusion process. It can handle both continuous and binary data, with various
 * imputation methods and normalization options. The cross-validation process can be
 * parallelized for improved performance on multi-core systems.
 *
 * @param graph The adjacency list of the graph. Each inner vector contains the indices of neighboring vertices.
 * @param edge_lengths A vector of vectors of edge lengths.
 * @param weights A vector of weights associated with each vertex for the diffusion process.
 * @param y The initial values of the function defined on the vertices of the graph.
 * @param n_time_steps The number of time steps for the diffusion process.
 * @param step_factor The step size factor for the diffusion process.
 * @param normalize Normalization option: 0 for no normalization, 1 for range adjustment, 2 for mean adjustment.
 *
 * @param preserve_local_maxima If true, adjusts the diffusion process to preserve local maxima in the input signal.
 *        This is done by modifying the weights based on the proportion of neighbors with smaller values.
 *        Vertices that are local maxima will have their weights set to zero, preventing them from being smoothed.
 *        Default is false.
 *
 * @param local_maximum_weight_factor Scale factor applied to weights when preserving local maxima or extrema.
 *        Values less than 1.0 reduce the influence of preserved features. Default is 1.0.
 *
 * @param preserve_local_extrema If true, adjusts the diffusion process to preserve both local maxima and minima
 *        in the input signal. This is achieved by modifying the weights based on how extreme a vertex's value is
 *        compared to its neighbors. Vertices that are local extrema (either maxima or minima) will have their
 *        weights reduced, limiting the amount of smoothing applied to them. Default is false.
 *
 * @note preserve_local_maxima and preserve_local_extrema are mutually exclusive. Only one can be set to true at a time.
 *       If both are set to true, the function will throw an exception.
 *
 * @param imputation_method The method used for imputing values during cross-validation.
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param apply_binary_threshold Whether to apply binary thresholding for binary data.
 * @param binary_threshold The threshold value for binary classification.
 * @param ikernel An integer specifying the kernel function to use for distance-based calculations.
 *        Set ikernel to 0 to use unweighted means.
 * @param dist_normalization_factor A scaling factor applied to the maximum distance between a vertex
 *        and its neighbors. This ensures non-zero weights even when all distances are equal,
 *        by slightly increasing the normalization range. Default value is 1.01.
 * @param n_CVs The number of cross-validation rounds.
 * @param n_CV_folds The number of folds in each cross-validation round.
 * @param epsilon A small positive constant for numerical stability in binary Rf_error calculations.
 * @param seed A seed for the random number generator used in cross-validation.
 * @param n_cores The number of CPU cores to use for parallel computation. If set to 1 (default),
 *        the function runs sequentially. For values greater than 1, OpenMP is used to parallelize
 *        the cross-validation process. If set to 0 or a value larger than the available cores,
 *        uses all available cores.
 * @param verbose If true, prints progress information and diagnostic messages during computation.
 *        Default is false.
 *
 * @return A unique pointer to a graph_diffusion_smoother_result_t structure containing:
 *         - y_traj: A vector of vectors representing the trajectory of y values over time.
 *                   Each inner vector contains the smoothed y values for all vertices at that step.
 *         - cv_errors: An (n_time_steps)-by-(n_CVs) matrix stored as a vector, containing
 *                     cross-validation errors for each time step and CV round.
 *         - mean_cv_errors: A vector of length n_time_steps containing the mean CV errors
 *                          across all CV rounds for each time step.
 *         - y_optimal: The smoothed y values at the optimal time step (where mean CV Rf_error is minimal).
 *                     If no CV is performed (n_CVs = 0), contains the final smoothed values.
 *         - n_time_steps: The number of time steps used in the diffusion process.
 *         - n_CVs: The number of cross-validation rounds performed.
 *         - optimal_time_step: The time step where the mean CV Rf_error is minimal.
 *                             Set to -1 if no CV is performed or if an Rf_error occurs.
 *         - min_cv_error: The minimum mean CV Rf_error value across all time steps.
 *                        Set to infinity if no CV is performed or if an Rf_error occurs.
 *
 * @throws std::runtime_error if input parameters are invalid or inconsistent.
 *
 * @note The function handles errors gracefully by:
 *       - Resetting CV data and continuing without CV if cross-validation fails
 *       - Providing fallback values for optimal solutions in Rf_error cases
 *       - Maintaining thread safety in parallel computations
 *       - Reporting detailed Rf_error messages when verbose is true
 *
 * @note The function supports both hop-index and distance-based neighbor calculations,
 *       determined by whether the edge_lengths parameter contains distances.
 *       It also handles binary data differently, using log-likelihood for Rf_error calculation.
 *       When n_cores > 1, the cross-validation process is parallelized using OpenMP.
 */
std::unique_ptr<graph_diffusion_smoother_result_t>
graph_diffusion_smoother_mp(const std::vector<std::vector<int>>& graph,
                            const std::vector<std::vector<double>>& edge_lengths,
                            std::vector<double>& weights,
                            const std::vector<double>& y,
                            int n_time_steps,
                            double step_factor,
                            int normalize,
                            bool preserve_local_maxima = false,
                            double local_maximum_weight_factor = 1.0,
                            bool preserve_local_extrema = false,
                            imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                            iterative_imputation_params_t iterative_params = {},
                            bool apply_binary_threshold = true,
                            double binary_threshold = 0.5,
                            int ikernel = 1,
                            double dist_normalization_factor = 1.1,
                            int n_CVs = 0,
                            int n_CV_folds = 10,
                            double epsilon = 1e-10,
                            unsigned int seed = 0,
                            int n_cores = 1,
                            bool verbose = false) {

    if (verbose) {
        Rprintf("In graph_diffusion_smoother_mp()\n");
    }

    // All parameter tests are done in the R function calling this one
    int n_vertices = y.size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && imputation_method == imputation_method_t::GLOBAL_MEAN_THRESHOLD) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    // Set up OpenMP
    int max_threads = omp_get_max_threads();
    int num_threads = (n_cores > 0 && n_cores <= max_threads) ? n_cores : max_threads;
    omp_set_num_threads(num_threads);

    // Initialize thread-local kernel weights and vertex edge lengths vectors
    std::vector<std::vector<double>> thread_kernel_weights(num_threads, std::vector<double>(max_neighbors));
    std::vector<std::vector<double>> thread_vertex_edge_lengths(num_threads, std::vector<double>(max_neighbors));

    bool mean_adjust = false;
    bool range_adjust = false;
    if (normalize == 1) {
        range_adjust = true;
    } else if (normalize == 2) {
        mean_adjust = true;
    }

    // computing the mean of y
    double Ey = 0;
    if (mean_adjust) {
        Ey = std::accumulate(y.begin(), y.end(), 0.0) / n_vertices;
    }

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    // computing min and max values of y
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());

    // Define thread-safe kernel diffusion loop
    auto kernel_diffusion_loop = [&](std::vector<double>& y_current,
                                   std::vector<double>& y_next,
                                   std::vector<double>& local_kernel_weights,
                                   std::vector<double>& local_vertex_edge_lengths) {

        for (int vertex = 0; vertex < n_vertices; ++vertex) {

            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                local_vertex_edge_lengths[j] = edge_lengths[vertex][j];
                if (local_vertex_edge_lengths[j] > max_dist)
                    max_dist = local_vertex_edge_lengths[j];
            }

            //max_dist += 1e-6;
            max_dist *= dist_normalization_factor;

            // Normalizing vertex edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j)
                local_vertex_edge_lengths[j] /= max_dist;

            kernel_fn(local_vertex_edge_lengths.data(), n_vertex_neighbors, local_kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += local_kernel_weights[j] * y_current[neighbor];
                total_weight += local_kernel_weights[j];
            }

            average_neighbor_value /= total_weight;

            y_next[vertex] = y_current[vertex] + weights[vertex] * step_factor *
                            (average_neighbor_value - y_current[vertex]);
        }
    };

    // The main loop for y_traj computation
    auto y_traj = std::vector<std::vector<double>>();
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);
    int time_step = 0;

    // Use thread 0's buffers for the main trajectory computation
    std::vector<double>& main_kernel_weights = thread_kernel_weights[0];
    std::vector<double>& main_vertex_edge_lengths = thread_vertex_edge_lengths[0];

    while (time_step < n_time_steps) {
        if (verbose) {
            Rprintf("\rtime_step: %d", time_step);
        }

        y_traj.push_back(y_current);

        if (preserve_local_maxima) {
            weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));
            if (local_maximum_weight_factor < 1) {
                for (int vertex = 0; vertex < n_vertices; vertex++)
                    weights[vertex] = local_maximum_weight_factor * weights[vertex];
            }
            for (int vertex = 0; vertex < n_vertices; vertex++)
                weights[vertex] = 1 - weights[vertex];
        } else if (preserve_local_extrema) {
            weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));
            if (local_maximum_weight_factor < 1) {
                for (int vertex = 0; vertex < n_vertices; vertex++)
                    weights[vertex] = local_maximum_weight_factor * weights[vertex];
            }
            for (int vertex = 0; vertex < n_vertices; vertex++)
                weights[vertex] = weights[vertex] * (1 - weights[vertex]);
        }

        kernel_diffusion_loop(y_current, y_next, main_kernel_weights, main_vertex_edge_lengths);

        if (range_adjust) {
            scale_to_range(y_next, ymin, ymax);
        } else if (mean_adjust) {
            double Ey_next = std::accumulate(y_next.begin(), y_next.end(), 0.0) / n_vertices;
            double Delta = Ey - Ey_next;
            for (auto& value : y_next) {
                value += Delta;
            }
        }

        std::swap(y_current, y_next);
        time_step++;
    }

    y_traj.push_back(y_current);

    std::vector<double> cv_errors(n_CVs * n_time_steps, 0.0);

    if (n_CVs > 0) {
        try {
            // Creating a set version of the adjacency matrix of the graph
            std::vector<std::set<int>> set_graph(n_vertices);
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
            }

            int fold_size = n_vertices / n_CV_folds;

            // The main parallel cross-validation loop
            #ifdef _OPENMP
            #pragma omp parallel
            #endif
            {
                int thread_id = omp_get_thread_num();
                std::vector<double>& local_kernel_weights = thread_kernel_weights[thread_id];
                std::vector<double>& local_vertex_edge_lengths = thread_vertex_edge_lengths[thread_id];

                // Thread-local variables
                std::vector<double> local_y_current = y;
                std::vector<double> local_y_next(n_vertices);
                std::vector<double> local_weights = weights;

                // Use consistent RNG seeding
                std::mt19937 local_rng(seed);
                std::uniform_int_distribution<int> local_uni(0, n_vertices - 1);

                #ifdef _OPENMP
#pragma omp for schedule(dynamic)
                #endif
                for (int cv = 0; cv < n_CVs; ++cv) {
                    if (verbose) {
                        #ifdef _OPENMP
#pragma omp critical
                        #endif
                        {
                            Rprintf("\rcv step: %d", cv);
                        }
                    }

                    // Creating a test set
                    std::set<int> test_set;
                    if (fold_size == 1 && n_vertices == n_CVs) {
                        test_set.insert(cv);
                    } else {
                        while ((int)test_set.size() < fold_size) {
                            int vertex = local_uni(local_rng);
                            test_set.insert(vertex);
                        }
                    }

                    int cv_imputation_ikernel = ikernel;
                    if (ikernel == 0)
                        cv_imputation_ikernel = 1;

                    local_y_current = std::move(*cv_imputation(test_set,
                                                             graph,
                                                             edge_lengths,
                                                             y,
                                                             y_binary,
                                                             imputation_method,
                                                             iterative_params,
                                                             apply_binary_threshold,
                                                             binary_threshold,
                                                             cv_imputation_ikernel,
                                                             dist_normalization_factor));

                    // Perform diffusion smoothing on the CV graph
                    for (int local_time_step = 0; local_time_step < n_time_steps; ++local_time_step) {
                        if (preserve_local_maxima) {
                            local_weights = std::move(*prop_nbhrs_with_smaller_y(graph, local_y_current));
                            if (local_maximum_weight_factor < 1) {
                                for (int vertex = 0; vertex < n_vertices; vertex++)
                                    local_weights[vertex] = local_maximum_weight_factor * local_weights[vertex];
                            }
                            for (int vertex = 0; vertex < n_vertices; vertex++)
                                local_weights[vertex] = 1 - local_weights[vertex];
                        } else if (preserve_local_extrema) {
                            local_weights = std::move(*prop_nbhrs_with_smaller_y(graph, local_y_current));
                            if (local_maximum_weight_factor < 1) {
                                for (int vertex = 0; vertex < n_vertices; vertex++)
                                    local_weights[vertex] = local_maximum_weight_factor * local_weights[vertex];
                            }
                            for (int vertex = 0; vertex < n_vertices; vertex++)
                                local_weights[vertex] = local_weights[vertex] * (1 - local_weights[vertex]);
                        }

                        kernel_diffusion_loop(local_y_current, local_y_next,
                                           local_kernel_weights, local_vertex_edge_lengths);

                        if (range_adjust) {
                            scale_to_range(local_y_next, ymin, ymax);
                        } else if (mean_adjust) {
                            double Ey_next = std::accumulate(local_y_next.begin(), local_y_next.end(), 0.0) / n_vertices;
                            double Delta = Ey - Ey_next;
                            for (auto& value : local_y_next) {
                                value += Delta;
                            }
                        }

                        // Computing CV Rf_error at the test vertices
                        double cv_error = 0.0;
                        if (y_binary) {
                            for (const auto& vertex : test_set) {
                                double clipped_y_next = std::clamp(local_y_next[vertex], epsilon, 1.0 - epsilon);
                                cv_error += y[vertex] * log(clipped_y_next) +
                                          (1 - y[vertex]) * log(1 - clipped_y_next);
                            }
                            cv_error *= -1;
                        } else {
                            for (const auto& vertex : test_set)
                                cv_error += std::abs(local_y_next[vertex] - y[vertex]);
                            cv_error /= test_set.size();
                        }

                        cv_errors[local_time_step + cv * n_time_steps] = cv_error;

                        std::swap(local_y_current, local_y_next);
                    }
                }
            }
        } catch (const std::exception& e) {
            Rprintf("\nError in cross-validation: %s\n", e.what());
            // Reset CV data on Rf_error
            std::fill(cv_errors.begin(), cv_errors.end(), 0.0);
            n_CVs = 0;  // Indicate that CV failed
        }
    }

    // Create and populate result
    auto result = std::make_unique<graph_diffusion_smoother_result_t>();
    result->n_time_steps = n_time_steps;
    result->n_CVs = n_CVs;
    result->y_traj = std::move(y_traj);
    result->cv_errors = std::move(cv_errors);

    // Initialize optimal values
    result->optimal_time_step = -1;
    result->min_cv_error = std::numeric_limits<double>::infinity();

    if (n_CVs > 0) {
        try {
            // Compute mean CV errors
            result->mean_cv_errors.resize(n_time_steps, 0.0);
            for (int time_step = 0; time_step < n_time_steps; ++time_step) {
                double sum_error = 0.0;
                for (int cv = 0; cv < n_CVs; ++cv) {
                    sum_error += result->cv_errors[time_step + cv * n_time_steps];
                }
                result->mean_cv_errors[time_step] = sum_error / n_CVs;
            }

            // Find optimal time step
            if (!result->mean_cv_errors.empty()) {
                auto min_iter = std::min_element(result->mean_cv_errors.begin(),
                                                 result->mean_cv_errors.end());
                result->optimal_time_step = std::distance(result->mean_cv_errors.begin(),
                                                          min_iter);
                result->min_cv_error = *min_iter;

                // Store optimal solution
                if (result->optimal_time_step >= 0 &&
                    result->optimal_time_step < static_cast<int>(result->y_traj.size())) {
                    result->y_optimal = result->y_traj[result->optimal_time_step];
                } else {
                    if (verbose) {
                        Rprintf("\nWarning: optimal_time_step out of bounds for y_traj\n");
                    }
                    result->y_optimal = y;
                }
            } else {
                if (verbose) {
                    Rprintf("\nWarning: mean_cv_errors vector is empty\n");
                }
                result->y_optimal = y;
            }
        } catch (const std::exception& e) {
            if (verbose) {
                Rprintf("\nError computing optimal time step: %s\n", e.what());
            }
            result->y_optimal = y;
        }
    } else {
        // When no CV is performed, use the final smoothed values
        result->y_optimal = !result->y_traj.empty() ? result->y_traj.back() : y;
    }

    if (verbose) {
        Rprintf("\nOptimal time step: %d (CV Rf_error: %f)\n",
                result->optimal_time_step,
                result->min_cv_error);
    }

    return result;
}


// serial extended version
/**
 * @brief Applies a graph diffusion smoothing algorithm with cross-validation.
 *
 * This function smooths a function `y` defined over the vertices of a graph using
 * a diffusion process. It can handle both continuous and binary data, with various
 * imputation methods and normalization options.
 *
 * @param graph The adjacency list of the graph. Each inner vector contains the indices of neighboring vertices.
 * @param edge_lengths A vector of vectors of edge lengths.
 * @param weights A vector of weights associated with each vertex for the diffusion process.
 * @param y The initial values of the function defined on the vertices of the graph.
 * @param n_time_steps The number of time steps for the diffusion process.
 * @param step_factor The step size factor for the diffusion process.
 * @param normalize Normalization option: 0 for no normalization, 1 for range adjustment, 2 for mean adjustment.
 *
 * @param preserve_local_maxima If true, adjusts the diffusion process to preserve local maxima in the input signal.
 *        This is done by modifying the weights based on the proportion of neighbors with smaller values.
 *        Vertices that are local maxima will have their weights set to zero, preventing them from being smoothed.
 *        Default is false.
 *
 * @param preserve_local_extrema If true, adjusts the diffusion process to preserve both local maxima and minima
 *        in the input signal. This is achieved by modifying the weights based on how extreme a vertex's value is
 *        compared to its neighbors. Vertices that are local extrema (either maxima or minima) will have their
 *        weights reduced, limiting the amount of smoothing applied to them. Default is false.
 *
 * @note preserve_local_maxima and preserve_local_extrema are mutually exclusive. Only one can be set to true at a time.
 *       If both are set to true, the function will throw an exception.
 *
 * @param imputation_method The method used for imputing values during cross-validation.
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param apply_binary_threshold Whether to apply binary thresholding for binary data.
 * @param binary_threshold The threshold value for binary classification.
 * @param ikernel An integer specifying the kernel function to use for distance-based calculations. Set ikernel to 0 to use unweighted means.
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization
 * range. Default value is 1.01.
 * @param n_CVs The number of cross-validation rounds.
 * @param n_CV_folds The number of folds in each cross-validation round.
 * @param epsilon A small positive constant for numerical stability in binary Rf_error calculations.
 * @param seed A seed for the random number generator used in cross-validation.
 *
 * @return A unique pointer to a graph_diffusion_smoother_result_t structure containing:
 *         - y_traj: A vector of vectors representing the trajectory of y values at each time step.
 *                   Each inner vector contains the smoothed y values for all vertices at that step.
 *         - cv_errors: An (n_time_steps)-by-(n_CVs) matrix stored as a vector, containing
 *                     cross-validation errors for each time step and CV round.
 *         - mean_cv_errors: A vector of length n_time_steps containing the mean CV errors
 *                          across all CV rounds for each time step.
 *         - y_optimal: The smoothed y values at the optimal time step (where mean CV Rf_error is minimal).
 *                     If no CV is performed (n_CVs = 0), contains the final smoothed values.
 *         - n_time_steps: The number of time steps used in the diffusion process.
 *         - n_CVs: The number of cross-validation rounds performed.
 *         - optimal_time_step: The time step where the mean CV Rf_error is minimal.
 *                             Set to -1 if no CV is performed or if an Rf_error occurs.
 *         - min_cv_error: The minimum mean CV Rf_error value across all time steps.
 *                        Set to infinity if no CV is performed or if an Rf_error occurs.
 *
 * @throws std::runtime_error if input parameters are invalid or inconsistent.
 *
 * @note The function supports both hop-index and distance-based neighbor calculations,
 *       determined by whether the 'd' parameter is empty or not.
 *       It also handles binary data differently, using log-likelihood for Rf_error calculation.
 *
 */
std::unique_ptr<graph_diffusion_smoother_result_t>
ext_graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                         const std::vector<std::vector<double>>& edge_lengths,
                         std::vector<double>& weights,
                         const std::vector<double>& y,
                         int n_time_steps,
                         double step_factor,
                         int normalize,
                         bool preserve_local_maxima = false,
                         double local_maximum_weight_factor = 1.0,
                         bool preserve_local_extrema = false,
                         imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                         iterative_imputation_params_t iterative_params = {},
                         bool apply_binary_threshold = true,
                         double binary_threshold = 0.5,
                         int ikernel = 1,
                         double dist_normalization_factor = 1.1,
                         int n_CVs = 0,
                         int n_CV_folds = 10,
                         double epsilon = 1e-10,
                         unsigned int seed = 0) {

    // All parameter tests are done in the R function calling this one
    int n_vertices = y.size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Determine the maximum number of neighbors any vertex has; NOTE: Only needed if ikernel > 0
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and vertex edge lengths vectors; NOTE: Only needed if ikernel > 0
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && imputation_method == imputation_method_t::GLOBAL_MEAN_THRESHOLD) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    bool mean_adjust = false;
    bool range_adjust = false;
    if (normalize == 1) {
        range_adjust = true;
    } else if (normalize == 2) {
        mean_adjust = true;
    }

    // computing the mean of y
    double Ey = 0;
    if (mean_adjust) {
        Ey = std::accumulate(y.begin(), y.end(), 0.0) / n_vertices;
    }

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    // computing min and max values of y
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());

    /**
     * @brief Performs one step of kernel-based diffusion smoothing on a graph.
     *
     * This lambda function applies a single iteration of diffusion smoothing using a kernel function.
     * It updates the values of all vertices based on their neighbors' values, weighted by a kernel function.
     *
     * @param y_current [in] The current values of all vertices before this diffusion step.
     * @param y_next [out] The updated values of all vertices after this diffusion step.
     *
     * @details
     * For each vertex:
     * 1. Computes the maximum distance to its neighbors.
     * 2. Normalizes the edge lengths by this maximum distance.
     * 3. Applies the kernel function to these normalized distances.
     * 4. Computes a weighted average of the neighbors' values using the kernel weights.
     * 5. Updates the vertex's value using this weighted average and the step factor.
     *
     * The function uses the following captured variables from its enclosing scope:
     * - graph: The adjacency list representation of the graph.
     * - edge_lengths: The lengths of edges between vertices.
     * - weights: The individual weights for each vertex.
     * - n_neighbors: The number of neighbors for each vertex.
     * - step_factor: The factor controlling the magnitude of each diffusion step.
     * - kernel_fn: The kernel function used for weighting.
     * - vertex_edge_lengths: A pre-allocated vector for storing edge lengths.
     * - kernel_weights: A pre-allocated vector for storing kernel weights.
     *
     * @note This function assumes that all necessary variables are properly initialized in the enclosing scope.
     * @note Isolated vertices (those with no neighbors) are skipped in the diffusion process.
     */
    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {

        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue; // Ensure no isolated vertices

            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                if (vertex_edge_lengths[j] > max_dist)
                    max_dist = vertex_edge_lengths[j];
            }

            //max_dist += 1e-6;  // Add a small value to prevent division by zero and ensure non-zero weights
            max_dist *= dist_normalization_factor;

            // Normalizing vertex edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j)
                vertex_edge_lengths[j] /= max_dist;

            kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }

            average_neighbor_value /= total_weight;

            y_next[vertex] = y_current[vertex] + weights[vertex] * step_factor * (average_neighbor_value - y_current[vertex]);
        }
    };

    // The main loop:
    auto y_traj = std::vector<std::vector<double>>();
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);
    int time_step = 0;
    while (time_step < n_time_steps) {

        Rprintf("\rtime_step: %d",time_step);

        y_traj.push_back(y_current);

        if (preserve_local_maxima) {
            // Preserve local maxima by adjusting weights based on local topology
            // This helps maintain the local maxima of y while still smoothing it

            // Compute the proportion of neighbors with smaller y values for each vertex
            weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));

            if (local_maximum_weight_factor < 1) {
                for (int vertex = 0; vertex < n_vertices; vertex++)
                    weights[vertex] = local_maximum_weight_factor * weights[vertex];
            }

            // Redefining weights so that the weights are zero at the local maxima
            for (int vertex = 0; vertex < n_vertices; vertex++)
                weights[vertex] = 1 - weights[vertex];

        } else if (preserve_local_extrema) {
            // Preserve local extrema by adjusting weights based on local topology
            // This helps maintain important features of y while still smoothing it

            // Compute the proportion of neighbors with smaller y values for each vertex
            weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));

            if (local_maximum_weight_factor < 1) {
                for (int vertex = 0; vertex < n_vertices; vertex++)
                    weights[vertex] = local_maximum_weight_factor * weights[vertex];
            }

            // Adjust weights to emphasize vertices that are local extrema
            // Vertices with p close to 0 or 1 are likely local extrema
            // The weight formula p * (1-p) gives lower weights to these vertices,
            // preserving their values in the smoothing process
            for (int vertex = 0; vertex < n_vertices; vertex++)
                weights[vertex] = weights[vertex] * (1 - weights[vertex]);
        }

        kernel_diffusion_loop(y_current, y_next);

        // Apply normalization if needed
        if (range_adjust) {
            scale_to_range(y_next, ymin, ymax);
        } else if (mean_adjust) {
            double Ey_next = std::accumulate(y_next.begin(), y_next.end(), 0.0) / n_vertices;
            double Delta = Ey - Ey_next;
            for (auto& value : y_next) {
                value += Delta;
            }
        }

        // Swap y_current and y_next
        std::swap(y_current, y_next);

        time_step++;
    }

    y_traj.push_back(y_current);

    std::vector<double> cv_errors(n_CVs * n_time_steps, 0.0); // Initialize with 0.0 when n_CVs = 0

    if (n_CVs > 0) {

        // Creating a set version of the adjacency matrix of the graph
        std::vector<std::set<int>> set_graph(n_vertices);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
        }

        int fold_size = n_vertices / n_CV_folds;

        //
        // The main cross-validation loop
        //
        for (int cv = 0; cv < n_CVs; ++cv) {

            Rprintf("\rcv step: %d",cv);

            // Creating a test set
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

            int cv_imputation_ikernel = ikernel;
            if (ikernel == 0)
                cv_imputation_ikernel = 1;

            y_current = std::move(*cv_imputation(test_set,
                                                 graph,
                                                 edge_lengths,
                                                 y,
                                                 y_binary,
                                                 imputation_method,
                                                 iterative_params,
                                                 apply_binary_threshold,
                                                 binary_threshold,
                                                 cv_imputation_ikernel,
                                                 dist_normalization_factor));

            // Perform diffusion smoothing on the CV graph at each step estimating the MAD Rf_error
            time_step = 0;
            while (time_step < n_time_steps) {

                if (preserve_local_maxima) {
                    weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));
                    if (local_maximum_weight_factor < 1) {
                        for (int vertex = 0; vertex < n_vertices; vertex++)
                            weights[vertex] = local_maximum_weight_factor * weights[vertex];
                    }
                    for (int vertex = 0; vertex < n_vertices; vertex++)
                        weights[vertex] = 1 - weights[vertex];

                } else if (preserve_local_extrema) {
                    weights = std::move(*prop_nbhrs_with_smaller_y(graph, y_current));
                    if (local_maximum_weight_factor < 1) {
                        for (int vertex = 0; vertex < n_vertices; vertex++)
                            weights[vertex] = local_maximum_weight_factor * weights[vertex];
                    }
                    for (int vertex = 0; vertex < n_vertices; vertex++)
                        weights[vertex] = weights[vertex] * (1 - weights[vertex]);
                }

                kernel_diffusion_loop(y_current, y_next);

                // Apply normalization if needed
                if (range_adjust) {
                    scale_to_range(y_next, ymin, ymax);
                } else if (mean_adjust) {
                    double Ey_next = std::accumulate(y_next.begin(), y_next.end(), 0.0) / n_vertices;
                    double Delta = Ey - Ey_next;
                    for (auto& value : y_next) {
                        value += Delta;
                    }
                }

                // Computing CV Rf_error at the test vertices
                double cv_error = 0.0;
                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        double clipped_y_next = std::clamp(y_next[vertex], epsilon, 1.0 - epsilon);
                        cv_error += y[vertex] * log(clipped_y_next) + (1 - y[vertex]) * log(1 - clipped_y_next);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set)
                        cv_error += std::abs(y_next[vertex] - y[vertex]);
                    cv_error /= test_set.size();
                }
                cv_errors[time_step + cv * n_time_steps] = cv_error;

                // Swap y_current and y_next
                std::swap(y_current, y_next);

                time_step++;
            }
        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    auto result = std::make_unique<graph_diffusion_smoother_result_t>();
    result->n_time_steps = n_time_steps;
    result->n_CVs = n_CVs;
    result->y_traj = std::move(y_traj);
    result->cv_errors = std::move(cv_errors);

    // Initialize optimal time step and minimum CV Rf_error with default values
    result->optimal_time_step = -1;  // -1 indicates no optimal time step found
    result->min_cv_error = std::numeric_limits<double>::infinity();

    // Compute mean CV errors and find optimal time step
    if (n_CVs > 0) {
        try {
            // Resize mean_cv_errors vector
            result->mean_cv_errors.resize(n_time_steps, 0.0);

            // Compute mean CV errors for each time step
            for (int time_step = 0; time_step < n_time_steps; ++time_step) {
                double sum_error = 0.0;
                for (int cv = 0; cv < n_CVs; ++cv) {
                    sum_error += result->cv_errors[time_step + cv * n_time_steps];
                }
                result->mean_cv_errors[time_step] = sum_error / n_CVs;
            }

            // Find the optimal time step (index of global minimum)
            if (!result->mean_cv_errors.empty()) {
                auto min_iter = std::min_element(result->mean_cv_errors.begin(),
                                                 result->mean_cv_errors.end());

                result->optimal_time_step = std::distance(result->mean_cv_errors.begin(),
                                                          min_iter);
                result->min_cv_error = *min_iter;

                // Store the optimal solution
                if (result->optimal_time_step >= 0 &&
                    result->optimal_time_step < static_cast<int>(result->y_traj.size())) {
                    result->y_optimal = result->y_traj[result->optimal_time_step];
                } else {
                    Rprintf("\nWarning: optimal_time_step out of bounds for y_traj\n");
                    // Initialize y_optimal with original y values as fallback
                    result->y_optimal = y;
                }

                // Rprintf("\nOptimal time step: %d (CV Rf_error: %f)\n",
                //         result->optimal_time_step,
                //         result->min_cv_error);
            } else {
                Rprintf("\nWarning: mean_cv_errors vector is empty\n");
                // Initialize y_optimal with original y values as fallback
                result->y_optimal = y;
            }
        } catch (const std::exception& e) {
            Rprintf("\nError computing optimal time step: %s\n", e.what());
            // Keep default values and initialize y_optimal with original y values
            result->y_optimal = y;
        }
    } else {
        Rprintf("\nNo cross-validation performed (n_CVs = 0)\n");
        // When no CV is performed, use the final smoothed values
        result->y_optimal = !result->y_traj.empty() ? result->y_traj.back() : y;
    }
    return result;
}



/**
 * @brief R interface for Graph Diffusion Smoothing with cross-validation
 *
 * This function provides an R interface to the C++ graph diffusion smoothing algorithm.
 * It processes input data from R, performs graph diffusion smoothing with optional
 * cross-validation, and returns results in a format compatible with R.
 *
 * @param Rgraph SEXP: List of integer vectors representing the graph adjacency list.
 *               Each vector contains 1-based indices of neighboring vertices.
 * @param Redge_length SEXP: List of numeric vectors containing edge lengths corresponding to neighbors.
 * @param Rweights SEXP: Numeric vector of weights for each vertex.
 * @param Ry SEXP: Numeric vector of values to be smoothed.
 * @param Rn_time_steps SEXP: Integer specifying the number of diffusion steps.
 * @param Rstep_factor SEXP: Numeric value controlling the magnitude of each diffusion step.
 * @param Rnormalize SEXP: Integer specifying normalization method (0: none, 1: range, 2: mean).
 * @param Rpreserve_local_maxima SEXP: Logical indicating whether to preserve local maxima.
 * @param Rlocal_maximum_weight_factor SEXP: Numeric value between 0 and 1 for weight adjustment.
 * @param Rpreserve_local_extrema SEXP: Logical indicating whether to preserve local extrema.
 * @param Rimputation_method SEXP: Integer specifying the imputation method for cross-validation.
 * @param Rmax_iterations SEXP: Integer specifying maximum iterations for iterative imputation.
 * @param Rconvergence_threshold SEXP: Numeric threshold for convergence in iterative methods.
 * @param Rapply_binary_threshold SEXP: Logical indicating whether to apply binary thresholding.
 * @param Rbinary_threshold SEXP: Numeric threshold value for binary classification.
 * @param Rikernel SEXP: Integer specifying the kernel function index.
 * @param Rdist_normalization_factor SEXP: Numeric factor for distance normalization (> 1.0).
 * @param Rn_CVs SEXP: Integer specifying number of cross-validation rounds.
 * @param Rn_CV_folds SEXP: Integer specifying number of folds in each CV round.
 * @param Repsilon SEXP: Numeric value for numerical stability in calculations.
 * @param Rseed SEXP: Integer seed for random number generation.
 * @param Rn_cores SEXP: Integer specifying number of CPU cores to use (>1 for parallel).
 * @param Rverbose SEXP: Logical controlling progress output.
 *
 * @return SEXP: A named list containing:
 *         - y_traj: List of numeric vectors representing the smoothing trajectory.
 *                   Each vector contains smoothed values at a time step.
 *         - cv_errors: Matrix (n_time_steps × n_CVs) of cross-validation errors.
 *                     Empty matrix if n_CVs = 0.
 *         - mean_cv_errors: Numeric vector of mean CV errors per time step.
 *         - y_optimal: Numeric vector of optimally smoothed values based on CV.
 *                     Final smoothed values if no CV performed.
 *         - n_time_steps: Integer indicating number of time steps used.
 *         - n_CVs: Integer indicating number of CV rounds performed.
 *         - optimal_time_step: Integer indicating time step with minimal CV Rf_error.
 *                             -1 if no CV performed.
 *         - min_cv_error: Numeric value of minimum mean CV Rf_error.
 *                        Inf if no CV performed.
 *
 * @throws Rf_error with descriptive message for:
 *         - Invalid input dimensions or parameter values
 *         - Memory allocation failures
 *         - Computation errors
 *         - Invalid result states
 *         - Inconsistent dimensions in results
 *
 * @note This function performs extensive parameter validation and Rf_error checking.
 *       It handles memory protection for R objects and ensures proper cleanup.
 *       Cross-validation can be parallelized when n_cores > 1.
 *       Progress information is printed when verbose = TRUE.
 *
 * @note The function supports both sequential and parallel computation modes:
 *       - Sequential mode (n_cores = 1): Uses graph_diffusion_smoother
 *       - Parallel mode (n_cores > 1): Uses graph_diffusion_smoother_mp
 *
 * @see graph_diffusion_smoother
 * @see graph_diffusion_smoother_mp
 */
SEXP S_ext_graph_diffusion_smoother(SEXP Rgraph,
                                    SEXP Redge_length,
                                    SEXP Rweights,
                                    SEXP Ry,
                                    SEXP Rn_time_steps,
                                    SEXP Rstep_factor,
                                    SEXP Rnormalize,
                                    SEXP Rpreserve_local_maxima,
                                    SEXP Rlocal_maximum_weight_factor,
                                    SEXP Rpreserve_local_extrema,
                                    SEXP Rimputation_method,
                                    SEXP Rmax_iterations,
                                    SEXP Rconvergence_threshold,
                                    SEXP Rapply_binary_threshold,
                                    SEXP Rbinary_threshold,
                                    SEXP Rikernel,
                                    SEXP Rdist_normalization_factor,
                                    SEXP Rn_CVs,
                                    SEXP Rn_CV_folds,
                                    SEXP Repsilon,
                                    SEXP Rseed,
                                    SEXP Rn_cores,
                                    SEXP Rverbose) {
    try {
        // Initial validation
        int n_vertices = LENGTH(Ry);
        if (n_vertices <= 0) {
            REPORT_ERROR("Input vector y must have positive length");
        }

        std::vector<std::vector<int>> graph = convert_adj_list_from_R(Rgraph);
        if (graph.size() != static_cast<size_t>(n_vertices)) {
            REPORT_ERROR("Graph size must match input vector length");
        }

        std::vector<std::vector<double>> edge_length = convert_weight_list_from_R(Redge_length);
        if (edge_length.size() != static_cast<size_t>(n_vertices)) {
            REPORT_ERROR("Edge lengths size must match input vector length");
        }

        int nprot = 0;

        // Parameter conversion with Rf_error checking
        PROTECT(Ry = Rf_coerceVector(Ry, REALSXP)); nprot++;
        double *y = REAL(Ry);
        if (!y) {
            UNPROTECT(nprot);
            Rf_error("Failed to convert input vector y");
        }

        PROTECT(Rweights = Rf_coerceVector(Rweights, REALSXP)); nprot++;
        double *weights = REAL(Rweights);
        if (!weights) {
            UNPROTECT(nprot);
            Rf_error("Failed to convert weights vector");
        }

        PROTECT(Rn_time_steps = Rf_coerceVector(Rn_time_steps, INTSXP)); nprot++;
        int n_time_steps = INTEGER(Rn_time_steps)[0];
        if (n_time_steps <= 0) {
            UNPROTECT(nprot);
            Rf_error("Number of time steps must be positive");
        }

        PROTECT(Rstep_factor = Rf_coerceVector(Rstep_factor, REALSXP)); nprot++;
        double step_factor = REAL(Rstep_factor)[0];
        if (step_factor <= 0) {
            UNPROTECT(nprot);
            Rf_error("Step factor must be positive");
        }

        PROTECT(Rnormalize = Rf_coerceVector(Rnormalize, INTSXP)); nprot++;
        int normalize = INTEGER(Rnormalize)[0];
        if (normalize < 0 || normalize > 2) {
            UNPROTECT(nprot);
            Rf_error("Normalize must be 0, 1, or 2");
        }

        PROTECT(Rpreserve_local_maxima = Rf_coerceVector(Rpreserve_local_maxima, LGLSXP)); nprot++;
        bool preserve_local_maxima = LOGICAL(Rpreserve_local_maxima)[0];

        PROTECT(Rlocal_maximum_weight_factor = Rf_coerceVector(Rlocal_maximum_weight_factor, REALSXP)); nprot++;
        double local_maximum_weight_factor = REAL(Rlocal_maximum_weight_factor)[0];
        if (local_maximum_weight_factor < 0 || local_maximum_weight_factor > 1) {
            UNPROTECT(nprot);
            Rf_error("Local maximum weight factor must be between 0 and 1");
        }

        PROTECT(Rpreserve_local_extrema = Rf_coerceVector(Rpreserve_local_extrema, LGLSXP)); nprot++;
        bool preserve_local_extrema = LOGICAL(Rpreserve_local_extrema)[0];

        if (preserve_local_maxima && preserve_local_extrema) {
            UNPROTECT(nprot);
            Rf_error("Cannot preserve both local maxima and local extrema simultaneously");
        }

        PROTECT(Rimputation_method = Rf_coerceVector(Rimputation_method, INTSXP)); nprot++;
        imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);

        PROTECT(Rmax_iterations = Rf_coerceVector(Rmax_iterations, INTSXP)); nprot++;
        int max_iterations = INTEGER(Rmax_iterations)[0];
        if (max_iterations <= 0) {
            UNPROTECT(nprot);
            Rf_error("Maximum iterations must be positive");
        }

        PROTECT(Rconvergence_threshold = Rf_coerceVector(Rconvergence_threshold, REALSXP)); nprot++;
        double convergence_threshold = REAL(Rconvergence_threshold)[0];
        if (convergence_threshold <= 0) {
            UNPROTECT(nprot);
            Rf_error("Convergence threshold must be positive");
        }

        PROTECT(Rapply_binary_threshold = Rf_coerceVector(Rapply_binary_threshold, LGLSXP)); nprot++;
        bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];

        PROTECT(Rbinary_threshold = Rf_coerceVector(Rbinary_threshold, REALSXP)); nprot++;
        double binary_threshold = REAL(Rbinary_threshold)[0];
        if (binary_threshold < 0 || binary_threshold > 1) {
            UNPROTECT(nprot);
            Rf_error("Binary threshold must be between 0 and 1");
        }

        PROTECT(Rikernel = Rf_coerceVector(Rikernel, INTSXP)); nprot++;
        int ikernel = INTEGER(Rikernel)[0];
        if (ikernel < 0) {
            UNPROTECT(nprot);
            Rf_error("Invalid kernel index");
        }

        PROTECT(Rdist_normalization_factor = Rf_coerceVector(Rdist_normalization_factor, REALSXP)); nprot++;
        double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];
        if (dist_normalization_factor <= 1.0) {
            UNPROTECT(nprot);
            Rf_error("Distance normalization factor must be greater than 1.0");
        }

        PROTECT(Rn_CVs = Rf_coerceVector(Rn_CVs, INTSXP)); nprot++;
        int n_CVs = INTEGER(Rn_CVs)[0];
        if (n_CVs < 0) {
            UNPROTECT(nprot);
            Rf_error("Number of CV rounds must be non-negative");
        }

        PROTECT(Rn_CV_folds = Rf_coerceVector(Rn_CV_folds, INTSXP)); nprot++;
        int n_CV_folds = INTEGER(Rn_CV_folds)[0];
        if (n_CVs > 0 && n_CV_folds <= 1) {
            UNPROTECT(nprot);
            Rf_error("Number of CV folds must be greater than 1 when performing cross-validation");
        }

        PROTECT(Repsilon = Rf_coerceVector(Repsilon, REALSXP)); nprot++;
        double epsilon = REAL(Repsilon)[0];
        if (epsilon <= 0) {
            UNPROTECT(nprot);
            Rf_error("Epsilon must be positive");
        }

        PROTECT(Rseed = Rf_coerceVector(Rseed, INTSXP)); nprot++;
        unsigned int seed = static_cast<unsigned int>(INTEGER(Rseed)[0]);

        PROTECT(Rn_cores = Rf_coerceVector(Rn_cores, INTSXP)); nprot++;
        int n_cores = INTEGER(Rn_cores)[0];
        if (n_cores < 0) {
            UNPROTECT(nprot);
            Rf_error("Number of cores must be non-negative");
        }

        PROTECT(Rverbose = Rf_coerceVector(Rverbose, LGLSXP)); nprot++;
        bool verbose = LOGICAL(Rverbose)[0];

        // Initialize parameters structure
        iterative_imputation_params_t iterative_params;
        iterative_params.max_iterations = max_iterations;
        iterative_params.convergence_threshold = convergence_threshold;

        // Create local copy of weights
        std::vector<double> local_weights(weights, weights + n_vertices);

        // Call appropriate smoother function
        std::unique_ptr<graph_diffusion_smoother_result_t> result;

        if (n_cores == 1) {
            result = ext_graph_diffusion_smoother(graph,
                                              edge_length,
                                              local_weights,
                                              std::vector<double>(y, y + n_vertices),
                                              n_time_steps,
                                              step_factor,
                                              normalize,
                                              preserve_local_maxima,
                                              local_maximum_weight_factor,
                                              preserve_local_extrema,
                                              imputation_method,
                                              iterative_params,
                                              apply_binary_threshold,
                                              binary_threshold,
                                              ikernel,
                                              dist_normalization_factor,
                                              n_CVs,
                                              n_CV_folds,
                                              epsilon,
                                              seed);
        } else {
            result = graph_diffusion_smoother_mp(graph,
                                                 edge_length,
                                                 local_weights,
                                                 std::vector<double>(y, y + n_vertices),
                                                 n_time_steps,
                                                 step_factor,
                                                 normalize,
                                                 preserve_local_maxima,
                                                 local_maximum_weight_factor,
                                                 preserve_local_extrema,
                                                 imputation_method,
                                                 iterative_params,
                                                 apply_binary_threshold,
                                                 binary_threshold,
                                                 ikernel,
                                                 dist_normalization_factor,
                                                 n_CVs,
                                                 n_CV_folds,
                                                 epsilon,
                                                 seed,
                                                 n_cores,
                                                 verbose);

        }

        // Check result
        if (!result) {
            UNPROTECT(nprot);
            Rf_error("Failed to compute graph diffusion smoothing");
        }

        // Validate result components
        if (result->y_traj.empty()) {
            UNPROTECT(nprot);
            Rf_error("Empty trajectory returned from graph diffusion smoothing");
        }

        if (result->y_optimal.empty()) {
            UNPROTECT(nprot);
            Rf_error("Empty optimal solution returned from graph diffusion smoothing");
        }

        // Add check for consistent dimensions
        if (!result->y_traj.empty() && result->y_traj[0].size() != static_cast<size_t>(n_vertices)) {
            UNPROTECT(nprot);
            Rf_error("Inconsistent dimensions in result trajectory");
        }

        if (result->y_optimal.size() != static_cast<size_t>(n_vertices)) {
            UNPROTECT(nprot);
            Rf_error("Inconsistent dimensions in optimal solution");
        }

        // Convert results to R objects
        SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 8)); nprot++;

        // Convert y_traj to R list of vectors
        SEXP r_y_traj = PROTECT(Rf_allocVector(VECSXP, result->y_traj.size())); nprot++;
        for (size_t t = 0; t < result->y_traj.size(); ++t) {
            SEXP r_vector = PROTECT(Rf_allocVector(REALSXP, n_vertices));
            std::copy(result->y_traj[t].begin(), result->y_traj[t].end(), REAL(r_vector));
            SET_VECTOR_ELT(r_y_traj, t, r_vector);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, 0, r_y_traj);

        // Convert cv_errors to R matrix
        if (result->n_CVs > 0) {
            if (result->cv_errors.size() != static_cast<size_t>(result->n_time_steps * result->n_CVs)) {
                UNPROTECT(nprot);
                Rf_error("Invalid cv_errors size");
            }
        }
        SEXP r_cv_errors = PROTECT(Rf_allocMatrix(REALSXP, result->n_time_steps, result->n_CVs)); nprot++;
        double* cv_errors_ptr = REAL(r_cv_errors);
        for (int t = 0; t < result->n_time_steps; ++t) {
            for (int cv = 0; cv < result->n_CVs; ++cv) {
                cv_errors_ptr[t + cv * result->n_time_steps] = result->cv_errors[t + cv * result->n_time_steps];
            }
        }
        SET_VECTOR_ELT(r_result, 1, r_cv_errors);

        // Convert mean_cv_errors to R vector
        SEXP r_mean_cv_errors = PROTECT(Rf_allocVector(REALSXP, result->mean_cv_errors.size())); nprot++;
        std::copy(result->mean_cv_errors.begin(), result->mean_cv_errors.end(), REAL(r_mean_cv_errors));
        SET_VECTOR_ELT(r_result, 2, r_mean_cv_errors);

        // Convert y_optimal to R vector
        SEXP r_y_optimal = PROTECT(Rf_allocVector(REALSXP, result->y_optimal.size())); nprot++;
        std::copy(result->y_optimal.begin(), result->y_optimal.end(), REAL(r_y_optimal));
        SET_VECTOR_ELT(r_result, 3, r_y_optimal);

        // Convert scalar values
        SEXP r_n_time_steps = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
        INTEGER(r_n_time_steps)[0] = result->n_time_steps;
        SET_VECTOR_ELT(r_result, 4, r_n_time_steps);

        SEXP r_n_CVs = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
        INTEGER(r_n_CVs)[0] = result->n_CVs;
        SET_VECTOR_ELT(r_result, 5, r_n_CVs);

        SEXP r_optimal_time_step = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
        INTEGER(r_optimal_time_step)[0] = result->optimal_time_step;
        SET_VECTOR_ELT(r_result, 6, r_optimal_time_step);

        SEXP r_min_cv_error = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
        REAL(r_min_cv_error)[0] = result->min_cv_error;
        SET_VECTOR_ELT(r_result, 7, r_min_cv_error);

        // Set names for the result list
        SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 8)); nprot++;
        SET_STRING_ELT(r_names, 0, Rf_mkChar("y_traj"));
        SET_STRING_ELT(r_names, 1, Rf_mkChar("cv_errors"));
        SET_STRING_ELT(r_names, 2, Rf_mkChar("mean_cv_errors"));
        SET_STRING_ELT(r_names, 3, Rf_mkChar("y_optimal"));
        SET_STRING_ELT(r_names, 4, Rf_mkChar("n_time_steps"));
        SET_STRING_ELT(r_names, 5, Rf_mkChar("n_CVs"));
        SET_STRING_ELT(r_names, 6, Rf_mkChar("optimal_time_step"));
        SET_STRING_ELT(r_names, 7, Rf_mkChar("min_cv_error"));
        Rf_setAttrib(r_result, R_NamesSymbol, r_names);

        // Add checks for allocation failure
        if (!r_result || !r_y_traj || !r_cv_errors || !r_mean_cv_errors ||
            !r_y_optimal || !r_n_time_steps || !r_n_CVs ||
            !r_optimal_time_step || !r_min_cv_error || !r_names) {
            UNPROTECT(nprot);
            Rf_error("Failed to allocate memory for results");
        }

        UNPROTECT(nprot);
        return r_result;

    } catch (const std::bad_alloc& e) {
        Rf_error("Memory allocation failed: %s", e.what());
    } catch (const std::exception& e) {
        Rf_error("Error during computation: %s", e.what());
    } catch (...) {
        Rf_error("Unknown Rf_error occurred");
    }

    // This return is never reached but silences compiler warnings
    return R_NilValue;
}


/*!
 * @brief Performs basic graph diffusion smoothing on vertex values.
 *
 * This function implements a kernel-based diffusion process on a graph where each vertex
 * value is updated based on weighted averages of its neighbors' values. The weights
 * are determined by applying a kernel function to normalized edge lengths.
 *
 * @param graph Vector of vectors representing the graph's adjacency list.
 *             Each inner vector contains the indices of neighboring vertices.
 * @param edge_lengths Vector of vectors containing the lengths/distances of edges.
 *                    edge_lengths[i][j] corresponds to the length of the edge between
 *                    vertex i and its j-th neighbor in graph[i].
 * @param y Initial values at each vertex.
 * @param n_time_steps Number of diffusion steps to perform.
 * @param step_factor Factor controlling the magnitude of each diffusion step (typically in (0,1]).
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels: 0-Constant, 1-Epanechnikov, 2-Triangular,
 *               3-TrExponential, 4-Laplace, 5-Normal, 6-Biweight, 7-Tricube, 8-Cosine
 * @param kernel_scale Scale parameter for Normal and Laplace kernels (default: 1.0).
 * @param dist_normalization_factor Factor used in edge length normalization (default: 1.01).
 *
 * @return Vector of vectors containing the trajectory of vertex values at each time step,
 *         including the initial and final states. Size is (n_time_steps + 1) × n_vertices.
 *
 * @note The function normalizes edge lengths by the maximum edge length for each vertex
 *       before applying the kernel function.
 *
 * @Rf_warning Vertices with no neighbors (isolated vertices) maintain their original values.
 */
std::vector<std::vector<double>> basic_graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                                                                const std::vector<std::vector<double>>& edge_lengths,
                                                                const std::vector<double>& y,
                                                                int n_time_steps,
                                                                double step_factor,
                                                                int ikernel = 1,
                                                                double kernel_scale = 1.0) {
    int n_vertices = y.size();
    initialize_kernel(ikernel, kernel_scale);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and vertex edge lengths vectors
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);

    // Calculate number of neighbors for each vertex
    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        n_neighbors[vertex] = graph[vertex].size();
    }

    // Define the diffusion step function
    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            // Find maximum distance to neighbors
            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                max_dist = std::max(max_dist, vertex_edge_lengths[j]);
            }

            max_dist += 1e-6;  // Prevent division by zero

            // Normalize edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] /= max_dist;
            }

            // Calculate kernel weights
            kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

            // Calculate weighted average of neighbor values
            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }

            average_neighbor_value /= total_weight;

            // Update vertex value
            y_next[vertex] = y_current[vertex] + step_factor * (average_neighbor_value - y_current[vertex]);
        }
    };

    // Main diffusion loop
    std::vector<std::vector<double>> y_traj;
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);

    int time_step = 0;
    while (time_step < n_time_steps) {
        y_traj.push_back(y_current);

        kernel_diffusion_loop(y_current, y_next);

        // Swap current and next states
        std::swap(y_current, y_next);

        time_step++;
    }

    // Add final state
    y_traj.push_back(y_current);

    return y_traj;
}

/*!
 * @brief Performs adaptive graph diffusion smoothing with momentum and dynamic step size adjustment.
 *
 * This function implements an advanced version of kernel-based diffusion on a graph where vertex
 * values are updated using adaptive step sizes and momentum terms. The adaptation mechanism
 * allows for faster convergence while avoiding oscillations.
 *
 * @details
 * The adaptive algorithm operates as follows:
 *
 * 1. Step Size Adaptation:
 *    - Each vertex maintains its own step size that adapts based on local behavior
 *    - Step sizes increase when changes are consistent (have the same sign)
 *    - Step sizes decrease when oscillations or overshooting is detected
 *    - Changes are measured through delta ratio: current_delta / previous_delta
 *    - Step sizes are bounded by [min_step, max_step] to ensure stability
 *
 * 2. Momentum Implementation:
 *    - Momentum terms track the history of changes at each vertex
 *    - New momentum = momentum * old_momentum + (1-momentum) * current_delta
 *    - Final update combines immediate delta with momentum term
 *    - Helps maintain consistent progress and overcome local barriers
 *
 * 3. Update Formula:
 *    For each vertex v:
 *    - delta = average_neighbor_value - current_value
 *    - momentum_term = momentum * old_momentum + (1-momentum) * delta
 *    - new_value = current_value + step_size * ((1-momentum) * delta + momentum * momentum_term)
 *
 * @param graph Vector of vectors representing the graph's adjacency list.
 *             Each inner vector contains the indices of neighboring vertices.
 * @param edge_lengths Vector of vectors containing the lengths/distances of edges.
 *                    edge_lengths[i][j] corresponds to the length of the edge between
 *                    vertex i and its j-th neighbor in graph[i].
 * @param y Initial values at each vertex.
 * @param n_time_steps Number of diffusion steps to perform.
 * @param base_step_factor Initial step size for all vertices (typically in (0,1]).
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels: 0-Constant, 1-Epanechnikov, 2-Triangular,
 *               3-TrExponential, 4-Laplace, 5-Normal, 6-Biweight, 7-Tricube, 8-Cosine
 * @param kernel_scale Scale parameter for Normal and Laplace kernels (default: 1.0).
 * @param dist_normalization_factor Factor used in edge length normalization (default: 1.01).
 * @param momentum Momentum coefficient controlling the influence of previous updates (default: 0.9).
 *                Higher values (closer to 1) give more weight to previous changes.
 * @param increase_factor Factor to multiply step size when changes are productive (default: 1.1).
 *                       Should be greater than 1.0.
 * @param decrease_factor Factor to multiply step size when oscillations occur (default: 0.8).
 *                       Should be less than 1.0.
 * @param min_step Minimum allowed step size to maintain progress (default: 0.01).
 * @param max_step Maximum allowed step size to prevent instability (default: 2.0).
 *
 * @return Vector of vectors containing the trajectory of vertex values at each time step,
 *         including the initial and final states. Size is (n_time_steps + 1) × n_vertices.
 *
 * @note
 * - The function normalizes edge lengths by the maximum edge length for each vertex
 *   before applying the kernel function.
 * - Step size adaptation is performed independently for each vertex.
 * - Momentum terms are maintained separately for each vertex.
 * - Isolated vertices (those with no neighbors) maintain their original values.
 *
 * @Rf_warning
 * - Very high increase_factor values may lead to instability
 * - Very low decrease_factor values may slow convergence unnecessarily
 * - High momentum values may cause overshooting in regions requiring precise adjustment
 * - The algorithm may require more memory than the basic version due to additional
 *   per-vertex storage for step sizes and momentum terms
 *
 * @see basic_graph_diffusion_smoother for the non-adaptive version
 */
std::vector<std::vector<double>> momentum_adaptive_graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                                                                            const std::vector<std::vector<double>>& edge_lengths,
                                                                            const std::vector<double>& y,
                                                                            int n_time_steps,
                                                                            double base_step_factor,
                                                                            int ikernel = 1,
                                                                            double kernel_scale = 1.0,
                                                                            double momentum = 0.9,
                                                                            double increase_factor = 1.1,
                                                                            double decrease_factor = 0.8,
                                                                            double min_step = 0.01,
                                                                            double max_step = 2.0) {

    int n_vertices = y.size();
    initialize_kernel(ikernel, kernel_scale);

    // Initialize step sizes and momentum for each vertex
    std::vector<double> step_sizes(n_vertices, base_step_factor);
    std::vector<double> prev_deltas(n_vertices, 0.0);
    std::vector<double> momentum_terms(n_vertices, 0.0);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and vertex edge lengths vectors
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);

    // Calculate number of neighbors for each vertex
    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        n_neighbors[vertex] = graph[vertex].size();
    }

    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            // Calculating average_neighbor_value
            // Find maximum distance to neighbors
            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                max_dist = std::max(max_dist, vertex_edge_lengths[j]);
            }

            max_dist += 1e-6;  // Prevent division by zero

            // Normalize edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] /= max_dist;
            }

            // Calculate kernel weights
            kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

            // Calculate weighted average of neighbor values
            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }

            average_neighbor_value /= total_weight;

            // Calculate current delta
            double delta = average_neighbor_value - y_current[vertex];

            // Update momentum term
            momentum_terms[vertex] = momentum * momentum_terms[vertex] + (1 - momentum) * delta;

            // Adapt step size based on delta history
            double oscillation_decrease_factor = 0.5;
            if (prev_deltas[vertex] != 0.0) {
                double delta_ratio = delta / prev_deltas[vertex];

                if (delta_ratio < 0) {
                    // Oscillation detected (sign change)
                    step_sizes[vertex] *= oscillation_decrease_factor;  // Aggressive decrease
                } else {
                    // No oscillation, check for overshooting
                    if (std::abs(delta) < std::abs(prev_deltas[vertex])) {
                        // Delta is decreasing - good progress
                        step_sizes[vertex] *= increase_factor;
                    } else {
                        // Overshooting detected
                        step_sizes[vertex] *= decrease_factor;  // Modest decrease
                    }
                }

                // Clamp step size to bounds
                step_sizes[vertex] = std::clamp(step_sizes[vertex], min_step, max_step);
            }

            // Update vertex value using momentum and adaptive step size
            y_next[vertex] = y_current[vertex] + step_sizes[vertex] * ((1 - momentum) * delta + momentum * momentum_terms[vertex]);

            // Store current delta for next iteration
            prev_deltas[vertex] = delta;
        }
    };

    // Main diffusion loop
    std::vector<std::vector<double>> y_traj;
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);

    int time_step = 0;
    while (time_step < n_time_steps) {
        y_traj.push_back(y_current);

        kernel_diffusion_loop(y_current, y_next);

        // Swap current and next states
        std::swap(y_current, y_next);

        time_step++;
    }

    // Add final state
    y_traj.push_back(y_current);

    return y_traj;
}


/*!
 * @brief Performs graph diffusion smoothing with adaptive step sizes.
 *
 * This function implements kernel-based diffusion on a graph where each vertex
 * maintains its own step size that adapts based on local behavior. The adaptation
 * mechanism responds to two types of non-optimal behavior:
 * 1. Oscillations (sign changes in delta) - aggressive step size reduction
 * 2. Overshooting (increasing magnitude without sign change) - modest step size reduction
 *
 * @param graph Vector of vectors representing the graph's adjacency list
 * @param edge_lengths Vector of vectors containing the lengths of edges
 * @param y Initial values at each vertex
 * @param n_time_steps Number of diffusion steps to perform
 * @param base_step_factor Initial step size for all vertices
 * @param ikernel Type of kernel function to use (default: 1)
 * @param kernel_scale Scale parameter for Normal and Laplace kernels (default: 1.0)
 * @param dist_normalization_factor Factor for edge length normalization (default: 1.01)
 * @param increase_factor Factor to multiply step size when changes are productive (default: 1.1)
 * @param decrease_factor Factor to multiply step size when overshooting (default: 0.8)
 * @param oscillation_factor Factor to multiply step size when oscillating (default: 0.5)
 * @param min_step Minimum allowed step size (default: 0.01)
 * @param max_step Maximum allowed step size (default: 2.0)
 */
std::vector<std::vector<double>> adaptive_graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                                                                   const std::vector<std::vector<double>>& edge_lengths,
                                                                   const std::vector<double>& y,
                                                                   int n_time_steps,
                                                                   double base_step_factor,
                                                                   int ikernel = 1,
                                                                   double kernel_scale = 1.0,
                                                                   double increase_factor = 1.1,
                                                                   double decrease_factor = 0.8,
                                                                   double oscillation_factor = 0.5,
                                                                   double min_step = 0.01,
                                                                   double max_step = 2.0) {

    int n_vertices = y.size();
    initialize_kernel(ikernel, kernel_scale);

    // Initialize per-vertex step sizes and previous deltas
    std::vector<double> step_sizes(n_vertices, base_step_factor);
    std::vector<double> prev_deltas(n_vertices, 0.0);

    // Initialize kernel-related variables
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);
    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        n_neighbors[vertex] = graph[vertex].size();
    }

    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            // Calculate weighted average of neighbors
            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                max_dist = std::max(max_dist, vertex_edge_lengths[j]);
            }
            max_dist += 1e-6;

            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] /= max_dist;
            }
            kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }
            average_neighbor_value /= total_weight;

            // Calculate current delta
            double delta = average_neighbor_value - y_current[vertex];

            // Adapt step size based on delta history
            if (prev_deltas[vertex] != 0.0) {
                double delta_ratio = delta / prev_deltas[vertex];
                if (delta_ratio < 0) {
                    // Oscillation detected (sign change)
                    step_sizes[vertex] *= oscillation_factor;
                } else {
                    // No oscillation, check for overshooting
                    if (std::abs(delta) < std::abs(prev_deltas[vertex])) {
                        // Delta is decreasing - good progress
                        step_sizes[vertex] *= increase_factor;
                    } else {
                        // Overshooting detected
                        step_sizes[vertex] *= decrease_factor;
                    }
                }
                step_sizes[vertex] = std::clamp(step_sizes[vertex], min_step, max_step);
            } // Note, that when prev_deltas[vertex] = 0.0, we do not adjust step_sizes[vertex], so it stays at what it was - is this an optimal strategy?

            // Simple update without momentum
            y_next[vertex] = y_current[vertex] + step_sizes[vertex] * delta;

            // Store current delta for next iteration
            prev_deltas[vertex] = delta;
        }
    };

    // Main diffusion loop
    std::vector<std::vector<double>> y_traj;
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);

    int time_step = 0;
    while (time_step < n_time_steps) {
        y_traj.push_back(y_current);
        kernel_diffusion_loop(y_current, y_next);
        std::swap(y_current, y_next);
        time_step++;
    }
    y_traj.push_back(y_current);

    return y_traj;
}



/*!
 * @brief Performs graph diffusion smoothing with adaptive step sizes and early stopping.
 *
 * This function implements kernel-based diffusion on a graph where each vertex
 * maintains its own step size that adapts based on local behavior. The process can
 * terminate early if convergence is detected.
 *
 * @details
 * Adaptation mechanism responds to two types of non-optimal behavior:
 * 1. Oscillations (sign changes in delta) - aggressive step size reduction
 * 2. Overshooting (increasing magnitude without sign change) - modest step size reduction
 *
 * Convergence is monitored through either:
 * - Average change: mean absolute difference between consecutive iterations
 * - Maximum change: maximum absolute difference between consecutive iterations
 *
 * The process stops when either:
 * 1. The selected change metric falls below the tolerance threshold
 * 2. The maximum number of time steps is reached
 *
 * @param graph Vector of vectors representing the graph's adjacency list
 * @param edge_lengths Vector of vectors containing the lengths of edges
 * @param y Initial values at each vertex
 * @param n_time_steps Maximum number of diffusion steps
 * @param base_step_factor Initial step size for all vertices
 * @param ikernel Type of kernel function to use (default: 1)
 * @param kernel_scale Scale parameter for Normal and Laplace kernels (default: 1.0)
 * @param dist_normalization_factor Factor for edge length normalization (default: 1.01)
 * @param increase_factor Factor to multiply step size when changes are productive (default: 1.1)
 * @param decrease_factor Factor to multiply step size when overshooting (default: 0.8)
 * @param oscillation_factor Factor to multiply step size when oscillating (default: 0.5)
 * @param min_step Minimum allowed step size (default: 0.01)
 * @param max_step Maximum allowed step size (default: 2.0)
 * @param convergence_criterion Type of convergence check to use (default: AVERAGE_CHANGE)
 * @param tolerance Threshold for convergence detection (default: 1e-6)
 * @param min_steps Minimum number of steps before checking convergence (default: 5)
 *
 * @return Vector of vectors containing the trajectory of vertex values up to convergence
 *         or maximum steps, including initial and final states.
 *
 * @note The function may return before n_time_steps iterations if convergence is detected,
 *       but will always perform at least min_steps iterations.
 */
std::vector<std::vector<double>> adaptive_graph_diffusion_smoother_with_stop(
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    int n_time_steps,
    double base_step_factor,
    int ikernel = 1,
    double kernel_scale = 1.0,
    double increase_factor = 1.1,
    double decrease_factor = 0.8,
    double oscillation_factor = 0.5,
    double min_step = 0.01,
    double max_step = 2.0,
    diffn_sm_convergence_criterion_t convergence_criterion = diffn_sm_convergence_criterion_t::AVERAGE_CHANGE,
    double tolerance = 1e-6,
    int min_steps = 5) {

    int n_vertices = y.size();
    initialize_kernel(ikernel, kernel_scale);

    // Helper functions for convergence checking
    auto compute_average_change = [](const std::vector<double>& y_current,
                                   const std::vector<double>& y_prev) -> double {
        double total_change = 0.0;
        for (size_t i = 0; i < y_current.size(); ++i) {
            total_change += std::abs(y_current[i] - y_prev[i]);
        }
        return total_change / y_current.size();
    };

    auto compute_max_change = [](const std::vector<double>& y_current,
                                const std::vector<double>& y_prev) -> double {
        double max_change = 0.0;
        for (size_t i = 0; i < y_current.size(); ++i) {
            max_change = std::max(max_change, std::abs(y_current[i] - y_prev[i]));
        }
        return max_change;
    };

    // Initialize per-vertex step sizes and previous deltas
    std::vector<double> step_sizes(n_vertices, base_step_factor);
    std::vector<double> prev_deltas(n_vertices, 0.0);

    // Initialize kernel-related variables
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);
    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        n_neighbors[vertex] = graph[vertex].size();
    }

    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            // Calculate weighted average of neighbors
            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                max_dist = std::max(max_dist, vertex_edge_lengths[j]);
            }
            max_dist += 1e-6;

            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] /= max_dist;
            }
            kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = graph[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }
            average_neighbor_value /= total_weight;

            // Calculate current delta
            double delta = average_neighbor_value - y_current[vertex];

            // Adapt step size based on delta history
            if (prev_deltas[vertex] != 0.0) {
                double delta_ratio = delta / prev_deltas[vertex];
                if (delta_ratio < 0) {
                    // Oscillation detected (sign change)
                    step_sizes[vertex] *= oscillation_factor;
                } else {
                    // No oscillation, check for overshooting
                    if (std::abs(delta) < std::abs(prev_deltas[vertex])) {
                        // Delta is decreasing - good progress
                        step_sizes[vertex] *= increase_factor;
                    } else {
                        // Overshooting detected
                        step_sizes[vertex] *= decrease_factor;
                    }
                }
                step_sizes[vertex] = std::clamp(step_sizes[vertex], min_step, max_step);
            } // Note, that when prev_deltas[vertex] = 0.0, we do not adjust step_sizes[vertex], so it stays at what it was - is this an optimal strategy?

            // Simple update without momentum
            y_next[vertex] = y_current[vertex] + step_sizes[vertex] * delta;

            // Store current delta for next iteration
            prev_deltas[vertex] = delta;
        }
    };

    // Main diffusion loop
    std::vector<std::vector<double>> y_traj;
    std::vector<double> y_current = y;
    std::vector<double> y_prev = y;
    std::vector<double> y_next(n_vertices);

    int time_step = 0;
    while (time_step < n_time_steps) {
        y_traj.push_back(y_current);
        y_prev = y_current;

        kernel_diffusion_loop(y_current, y_next);
        std::swap(y_current, y_next);

        // Check for convergence after minimum number of steps
        if (time_step >= min_steps) {
            double change = (convergence_criterion == diffn_sm_convergence_criterion_t::AVERAGE_CHANGE)
                          ? compute_average_change(y_current, y_prev)
                          : compute_max_change(y_current, y_prev);

            if (change < tolerance) {
                y_traj.push_back(y_current);  // Add final state
                return y_traj;
            }
        }

        time_step++;
    }

    y_traj.push_back(y_current);  // Add final state if max steps reached
    return y_traj;
}


// Other possible convergence criteria that could be added:

// Relative change: (current_change / initial_change) < relative_tolerance
// Moving average change: Check if change has been consistently small over several steps
// Energy-based criteria: Monitor the total "energy" of the system
// Gradient-based criteria: Check if the maximum gradient magnitude is below a threshold

// Additional features that could be useful:

// Return the number of steps taken
// Return the final change magnitude
// Option to track and return the convergence history
// Add a maximum runtime parameter


/*!
 * @brief Performs and instruments graph diffusion smoothing with adaptive step sizes.
 *
 * @details This function implements and monitors kernel-based or pure Laplacian diffusion on a graph
 * where each vertex maintains its own step size that adapts based on local behavior. The adaptation
 * mechanism responds to:
 * 1. Oscillations (sign changes in delta) - aggressive step size reduction
 * 2. Overshooting (increasing magnitude without sign change) - modest step size reduction
 *
 * The function tracks comprehensive performance metrics including:
 * - Evolution of vertex values over time
 * - Pre- and post-update deltas
 * - Step size adaptation history
 * - Global residual norms (L1)
 * - Energy metrics (smoothness, fidelity, Laplacian)
 * - Signal-to-noise ratio trajectory
 * - Curvature preservation metrics
 *
 * @param graph Vector of vectors representing the graph's adjacency list where graph[i] contains
 *        the indices of vertices adjacent to vertex i
 * @param edge_lengths Vector of vectors containing the lengths of edges where edge_lengths[i][j]
 *        is the length of the edge between vertex i and its j-th neighbor
 * @param y Initial values at each vertex (noisy data)
 * @param y_true Ground truth values at each vertex, used for performance evaluation
 * @param n_time_steps Number of diffusion steps to perform
 * @param base_step_factor Initial step size for all vertices
 * @param use_pure_laplacian If true, uses uniform weights (w_ij = 1) instead of kernel weights
 * @param ikernel Type of kernel function to use (ignored if use_pure_laplacian is true):
 *        - 1: Gaussian kernel
 *        - 2: Laplace kernel
 *        - 3: Uniform kernel
 * @param kernel_scale Scale parameter for Normal and Laplace kernels (ignored if use_pure_laplacian is true)
 * @param dist_normalization_factor Factor for edge length normalization (default: 1.01)
 * @param increase_factor Factor to multiply step size when changes are productive (default: 1.1)
 * @param decrease_factor Factor to multiply step size when overshooting (default: 0.8)
 * @param oscillation_factor Factor to multiply step size when oscillating (default: 0.5)
 * @param min_step Minimum allowed step size (default: 0.01)
 * @param max_step Maximum allowed step size (default: 2.0)
 *
 * @return A pair containing:
 *         - Vector of vectors containing vertex values at each time step
 *         - Performance metrics struct with the following members:
 *           - y_trajectory: Evolution of vertex values
 *           - pre_update_deltas: Delta values before updates
 *           - post_update_deltas: Delta values after updates
 *           - step_size_history: Step sizes for each vertex over time
 *           - global_residual_norm: L1 norm of deltas per iteration
 *           - max_absolute_delta: Maximum absolute delta per iteration
 *           - oscillation_events: Binary flags for oscillation detection
 *           - smoothness_energy: \f$ E_{smoothness} = \sum_{i} \sum_{j \in N(i)} w_{ij} (y_i - y_j)^2 \f$
 *           - fidelity_energy: \f$ E_{fidelity} = \sum_{i} (y_i - y^0_i)^2 \f$
 *           - laplacian_energy: \f$ E_{\Delta} = \sum_{i} |\sum_{j \in N(i)} w_{ij} (y_i - y_j)| \f$
 *           - energy_ratio: Ratio of smoothness to fidelity energy
 *           - initial_snr: Initial signal-to-noise ratio
 *           - snr_trajectory: SNR evolution over time
 *           - mean_absolute_deviation: MAD from ground truth over time
 *           - pointwise_curvature_error: Average absolute difference in Laplacian
 *           - integrated_curvature_error: Root mean square of Laplacian differences
 *
 * @note The function assumes that the graph is connected and that edge_lengths[i][j]
 *       corresponds to the edge between vertex i and the j-th vertex in graph[i].
 *
 * @Rf_warning Performance metrics involving y_true should only be used for algorithm
 *          tuning with synthetic data where ground truth is known. They are not
 *          available in real applications where true values are unknown.
 *
 * @see adaptive_graph_diffusion_smoother For the non-instrumented version
 *      graph_diffusion_smoother_performance_t For the performance metrics structure
 *
 * @example
 * ```cpp
 * // Create a simple chain graph
 * std::vector<std::vector<int>> graph = {{1}, {0, 2}, {1}};
 * std::vector<std::vector<double>> edge_lengths = {{1.0}, {1.0, 1.0}, {1.0}};
 *
 * // Generate synthetic data
 * std::vector<double> y_true = {0.0, 1.0, 0.0};
 * std::vector<double> y_noisy = {0.1, 0.9, 0.2};
 *
 * // Run instrumented diffusion
 * auto [trajectory, performance] = instrumented_gds(
 *     graph, edge_lengths, y_noisy, y_true,
 *     100,    // n_time_steps
 *     0.5,    // base_step_factor
 *     true    // use_pure_laplacian
 * );
 * ```
 */
graph_diffusion_smoother_performance_t
instrumented_gds(const std::vector<std::vector<int>>& graph,
                 const std::vector<std::vector<double>>& edge_lengths,
                 const std::vector<double>& y,
                 const std::vector<double>& y_true,
                 int n_time_steps,
                 double base_step_factor,
                 bool use_pure_laplacian = false,
                 int ikernel = 1,
                 double kernel_scale = 1.0,
                 double increase_factor = 1.1,
                 double decrease_factor = 0.8,
                 double oscillation_factor = 0.5,
                 double min_step = 0.01,
                 double max_step = 2.0) {

    // Input validation
    if (y.empty() || y.size() != y_true.size() || graph.size() != y.size() ||
        edge_lengths.size() != y.size()) {
        Rf_error("Invalid input dimensions");
    }

    int n_vertices = y.size();
    if (!use_pure_laplacian) {
        initialize_kernel(ikernel, kernel_scale);
    }

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and vertex edge lengths vectors
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);

    // Calculate number of neighbors for each vertex
    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        n_neighbors[vertex] = graph[vertex].size();
    }

    // Initialize performance tracking
    graph_diffusion_smoother_performance_t perf(n_vertices, n_time_steps);

    auto compute_smoothness_energy = [&](const std::vector<double>& y_current) {
        double energy = 0.0;
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            for (size_t j = 0; j < graph[vertex].size(); ++j) {
                int neighbor = graph[vertex][j];
                double weight = use_pure_laplacian ? 1.0 : kernel_weights[j];
                double diff = y_current[vertex] - y_current[neighbor];
                energy += weight * diff * diff;
            }
        }
        return energy;
    };

    auto compute_fidelity_energy = [&](const std::vector<double>& y_current) {
        double energy = 0.0;
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            double diff = y_current[vertex] - y[vertex];
            energy += diff * diff;
        }
        return energy;
    };

    // Initialize per-vertex step sizes and previous deltas
    std::vector<double> step_sizes(n_vertices, base_step_factor);
    std::vector<double> prev_deltas(n_vertices, 0.0);

    // Compute initial metrics
    perf.initial_snr = perf.compute_snr(y, y_true);

    // Helper functions for energy computations
    auto compute_laplacian = [&](const std::vector<double>& y_current, int vertex) {
        double laplacian = 0.0;
        for (size_t j = 0; j < graph[vertex].size(); ++j) {
            int neighbor = graph[vertex][j];
            double weight = use_pure_laplacian ? 1.0 : kernel_weights[j];
            laplacian += weight * (y_current[vertex] - y_current[neighbor]);
        }
        return laplacian;
    };

    auto compute_laplacian_energy = [&](const std::vector<double>& y_current) {
        double energy = 0.0;
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            energy += std::abs(compute_laplacian(y_current, vertex));
        }
        return energy;
    };

    auto compute_curvature_error = [&](const std::vector<double>& y_current) {
        double pointwise_sum = 0.0;
        double integrated_sum = 0.0;
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            double true_laplacian = compute_laplacian(y_true, vertex);
            double current_laplacian = compute_laplacian(y_current, vertex);
            double diff = current_laplacian - true_laplacian;
            pointwise_sum += std::abs(diff);
            integrated_sum += diff * diff;
        }
        return std::make_pair(pointwise_sum / n_vertices, std::sqrt(integrated_sum / n_vertices));
    };

    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        // Initialize delta vectors for this time step
        // Initialize delta vectors and oscillation events for this time step
        perf.pre_update_deltas.emplace_back(n_vertices);
        perf.post_update_deltas.emplace_back(n_vertices);
        perf.oscillation_events.emplace_back(n_vertices, false);
        std::vector<double> current_deltas(n_vertices);

        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            // Compute weights once if using kernel
            if (!use_pure_laplacian) {
                double max_dist = 0.0;
                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    vertex_edge_lengths[j] = edge_lengths[vertex][j];
                    max_dist = std::max(max_dist, vertex_edge_lengths[j]);
                }
                max_dist += 1e-6;

                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    vertex_edge_lengths[j] /= max_dist;
                }
                kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());
            }

            // Calculate weighted average using either pure or kernel-weighted Laplacian
            double weighted_sum = 0.0;
            double total_weight = 0.0;

            for (size_t j = 0; j < graph[vertex].size(); ++j) {
                int neighbor = graph[vertex][j];
                double weight = use_pure_laplacian ? 1.0 : kernel_weights[j];
                weighted_sum += weight * y_current[neighbor];
                total_weight += weight;
            }

            // Check for zero total weight
            double average_neighbor_value;
            if (total_weight > std::numeric_limits<double>::epsilon()) {
                average_neighbor_value = weighted_sum / total_weight;
            } else {
                // If total weight is effectively zero, keep current value
                average_neighbor_value = y_current[vertex];
            }

            double delta = average_neighbor_value - y_current[vertex];
            current_deltas[vertex] = delta;

            // Adapt step size based on delta history
            if (prev_deltas[vertex] != 0.0) {
                double delta_ratio = delta / prev_deltas[vertex];
                if (delta_ratio < 0) {
                    // Oscillation detected (sign change)
                    step_sizes[vertex] *= oscillation_factor;
                } else {
                    // No oscillation, check for overshooting
                    if (std::abs(delta) < std::abs(prev_deltas[vertex])) {
                        // Delta is decreasing - good progress
                        step_sizes[vertex] *= increase_factor;
                    } else {
                        // Overshooting detected
                        step_sizes[vertex] *= decrease_factor;
                    }
                }
                step_sizes[vertex] = std::clamp(step_sizes[vertex], min_step, max_step);
            } // Note, that when prev_deltas[vertex] = 0.0, we do not adjust step_sizes[vertex], so it stays at what it was - is this an optimal strategy?

            // Simple update without momentum
            y_next[vertex] = y_current[vertex] + step_sizes[vertex] * delta;

            // Store current delta for next iteration
            prev_deltas[vertex] = delta;
        }

        // Update performance metrics
        double l1_norm = 0.0;
        double max_delta = 0.0;
        for (double delta : current_deltas) {
            l1_norm += std::abs(delta);
            max_delta = std::max(max_delta, std::abs(delta));
        }

        // Update all performance metrics
        perf.global_residual_norm.push_back(l1_norm);
        perf.max_absolute_delta.push_back(max_delta);
        perf.step_size_history.push_back(step_sizes);

        double current_smoothness = compute_smoothness_energy(y_current);
        double current_fidelity = compute_fidelity_energy(y_current);

        // For energy calculations, add checks for NaN/Inf:
        auto is_finite = [](double x) {
            return !std::isnan(x) && !std::isinf(x);
        };

        auto safe_energy_update = [&](double energy_value, std::vector<double>& energy_vector) {
            if (is_finite(energy_value)) {
                energy_vector.push_back(energy_value);
            } else {
                // If energy is NaN/Inf, use previous value or zero
                double safe_value = energy_vector.empty() ? 0.0 : energy_vector.back();
                energy_vector.push_back(safe_value);
            }
        };

        safe_energy_update(current_smoothness, perf.smoothness_energy);
        safe_energy_update(current_fidelity, perf.fidelity_energy);

        perf.smoothness_energy.push_back(current_smoothness);
        perf.fidelity_energy.push_back(current_fidelity);

        // For energy ratio, add extra check for zero fidelity
        if (current_fidelity > std::numeric_limits<double>::epsilon()) {
            safe_energy_update(current_smoothness / current_fidelity, perf.energy_ratio);
        } else {
            safe_energy_update(perf.energy_ratio.empty() ? 0.0 : perf.energy_ratio.back(),
                               perf.energy_ratio);
        }

        safe_energy_update(compute_laplacian_energy(y_current), perf.laplacian_energy);

        auto [pointwise_curv, integrated_curv] = compute_curvature_error(y_current);
        perf.pointwise_curvature_error.push_back(pointwise_curv);
        perf.integrated_curvature_error.push_back(integrated_curv);

        perf.snr_trajectory.push_back(perf.compute_snr(y_current, y_true));
        perf.mean_absolute_deviation.push_back(perf.compute_mad(y_current, y_true));

        // Store pre-update deltas
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            perf.pre_update_deltas.back()[vertex] = current_deltas[vertex];
        }

        // After updates, store post-update deltas
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            perf.post_update_deltas.back()[vertex] = y_next[vertex] - y_current[vertex];
        }
    };

    // Main diffusion loop
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);

    int time_step = 0;
    while (time_step < n_time_steps) {
        perf.y_trajectory.push_back(y_current);
        kernel_diffusion_loop(y_current, y_next);
        std::swap(y_current, y_next);
        time_step++;
    }
    perf.y_trajectory.push_back(y_current);

    return perf;
}

/*!
 * @brief R interface for instrumented graph diffusion smoothing
 *
 * @details Converts R objects to C++ types, calls instrumented_gds, and converts results back to R
 *
 * @param s_graph SEXP containing adjacency list representation of graph
 * @param s_edge_lengths SEXP containing edge lengths for each vertex
 * @param s_y SEXP containing initial values at vertices
 * @param s_y_true SEXP containing ground truth values for performance evaluation
 * @param s_n_time_steps SEXP containing number of diffusion steps
 * @param s_base_step_factor SEXP containing initial step size
 * @param s_use_pure_laplacian SEXP logical for using uniform weights
 * @param s_ikernel SEXP containing kernel type (1-3)
 * @param s_kernel_scale SEXP containing kernel scale parameter
 * @param s_dist_normalization_factor SEXP containing edge length normalization factor
 * @param s_increase_factor SEXP containing step size increase factor
 * @param s_decrease_factor SEXP containing step size decrease factor
 * @param s_oscillation_factor SEXP containing oscillation reduction factor
 * @param s_min_step SEXP containing minimum step size
 * @param s_max_step SEXP containing maximum step size
 *
 * @return SEXP containing R list with performance metrics
 */
SEXP S_instrumented_gds(SEXP s_graph,
                        SEXP s_edge_lengths,
                        SEXP s_y,
                        SEXP s_y_true,
                        SEXP s_n_time_steps,
                        SEXP s_base_step_factor,
                        SEXP s_use_pure_laplacian,
                        SEXP s_ikernel,
                        SEXP s_kernel_scale,
                        SEXP s_increase_factor,
                        SEXP s_decrease_factor,
                        SEXP s_oscillation_factor,
                        SEXP s_min_step,
                        SEXP s_max_step) {
    int nprot = 0;

    // Convert input parameters
    std::vector<std::vector<int>> graph           = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(s_edge_lengths);

    PROTECT(s_y = Rf_coerceVector(s_y, REALSXP)); nprot++;
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

    PROTECT(s_y_true = Rf_coerceVector(s_y_true, REALSXP)); nprot++;
    std::vector<double> y_true(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));

    int n_time_steps = Rf_asInteger(s_n_time_steps);
    double base_step_factor = Rf_asReal(s_base_step_factor);
    bool use_pure_laplacian = Rf_asLogical(s_use_pure_laplacian);
    int ikernel = Rf_asInteger(s_ikernel);
    double kernel_scale = Rf_asReal(s_kernel_scale);
    double increase_factor = Rf_asReal(s_increase_factor);
    double decrease_factor = Rf_asReal(s_decrease_factor);
    double oscillation_factor = Rf_asReal(s_oscillation_factor);
    double min_step = Rf_asReal(s_min_step);
    double max_step = Rf_asReal(s_max_step);

    // Call C++ implementation
    auto perf = instrumented_gds(graph,
                                 edge_lengths,
                                 y,
                                 y_true,
                                 n_time_steps,
                                 base_step_factor,
                                 use_pure_laplacian,
                                 ikernel,
                                 kernel_scale,
                                 increase_factor,
                                 decrease_factor,
                                 oscillation_factor,
                                 min_step,
                                 max_step);

    // Convert performance metrics to R list
    const int N_RESULTS = 20;
    SEXP results = PROTECT(Rf_allocVector(VECSXP, N_RESULTS)); nprot++;

    // Trajectory and deltas
    SET_VECTOR_ELT(results, 0, convert_vector_vector_double_to_R(perf.y_trajectory));       nprot++ ;
    SET_VECTOR_ELT(results, 1, convert_vector_vector_double_to_R(perf.pre_update_deltas));  nprot++;
    SET_VECTOR_ELT(results, 2, convert_vector_vector_double_to_R(perf.post_update_deltas)); nprot++;
    SET_VECTOR_ELT(results, 3, convert_vector_vector_double_to_R(perf.step_size_history));  nprot++;

    // Scalar metrics
    SET_VECTOR_ELT(results, 4, convert_vector_double_to_R(perf.global_residual_norm)); nprot++;
    SET_VECTOR_ELT(results, 5, convert_vector_double_to_R(perf.max_absolute_delta));   nprot++;
    SET_VECTOR_ELT(results, 6, convert_vector_vector_bool_to_R(perf.oscillation_events)); nprot++;

    // Per-vertex event counts
    SET_VECTOR_ELT(results, 7, convert_vector_int_to_R(perf.oscillation_count_per_vertex)); nprot++;
    SET_VECTOR_ELT(results, 8, convert_vector_int_to_R(perf.increase_events_per_vertex));   nprot++;
    SET_VECTOR_ELT(results, 9, convert_vector_int_to_R(perf.decrease_events_per_vertex));   nprot++;
    SET_VECTOR_ELT(results, 10, convert_vector_int_to_R(perf.oscillation_reductions_per_vertex)); nprot++;

    // Energy metrics
    SET_VECTOR_ELT(results, 11, convert_vector_double_to_R(perf.smoothness_energy)); nprot++;
    SET_VECTOR_ELT(results, 12, convert_vector_double_to_R(perf.fidelity_energy));   nprot++;
    SET_VECTOR_ELT(results, 13, convert_vector_double_to_R(perf.laplacian_energy));  nprot++;
    SET_VECTOR_ELT(results, 14, convert_vector_double_to_R(perf.energy_ratio));      nprot++;

    // Truth-based metrics
    SEXP s_initial_snr;
    PROTECT(s_initial_snr = Rf_allocVector(REALSXP, 1)); nprot++;
    REAL(s_initial_snr)[0] = perf.initial_snr;
    SET_VECTOR_ELT(results, 15, s_initial_snr);

    SET_VECTOR_ELT(results, 16, convert_vector_double_to_R(perf.snr_trajectory));          nprot++;
    SET_VECTOR_ELT(results, 17, convert_vector_double_to_R(perf.mean_absolute_deviation)); nprot++;

    // Curvature errors
    SET_VECTOR_ELT(results, 18, convert_vector_double_to_R(perf.pointwise_curvature_error));  nprot++;
    SET_VECTOR_ELT(results, 19, convert_vector_double_to_R(perf.integrated_curvature_error)); nprot++;

    // Set names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_RESULTS)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("y_trajectory"));
    SET_STRING_ELT(names, 1, Rf_mkChar("pre_update_deltas"));
    SET_STRING_ELT(names, 2, Rf_mkChar("post_update_deltas"));
    SET_STRING_ELT(names, 3, Rf_mkChar("step_size_history"));
    SET_STRING_ELT(names, 4, Rf_mkChar("global_residual_norm"));
    SET_STRING_ELT(names, 5, Rf_mkChar("max_absolute_delta"));
    SET_STRING_ELT(names, 6, Rf_mkChar("oscillation_events"));
    SET_STRING_ELT(names, 7, Rf_mkChar("oscillation_count_per_vertex"));
    SET_STRING_ELT(names, 8, Rf_mkChar("increase_events_per_vertex"));
    SET_STRING_ELT(names, 9, Rf_mkChar("decrease_events_per_vertex"));
    SET_STRING_ELT(names, 10, Rf_mkChar("oscillation_reductions_per_vertex"));
    SET_STRING_ELT(names, 11, Rf_mkChar("smoothness_energy"));
    SET_STRING_ELT(names, 12, Rf_mkChar("fidelity_energy"));
    SET_STRING_ELT(names, 13, Rf_mkChar("laplacian_energy"));
    SET_STRING_ELT(names, 14, Rf_mkChar("energy_ratio"));
    SET_STRING_ELT(names, 15, Rf_mkChar("initial_snr"));
    SET_STRING_ELT(names, 16, Rf_mkChar("snr_trajectory"));
    SET_STRING_ELT(names, 17, Rf_mkChar("mean_absolute_deviation"));
    SET_STRING_ELT(names, 18, Rf_mkChar("pointwise_curvature_error"));
    SET_STRING_ELT(names, 19, Rf_mkChar("integrated_curvature_error"));

    Rf_setAttrib(results, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return results;
}


// serial streamlined version
/**
 * @brief Performs graph diffusion smoothing with cross-validation to find optimal diffusion time
 *
 * @details This function smooths vertex function values over a graph using diffusion with
 * cross-validation to identify the optimal diffusion time. The implementation uses a masked
 * diffusion approach where test vertices are excluded from the diffusion process by assigning
 * them zero weights, while still updating their values for Rf_error calculation.
 *
 * The algorithm works in three phases:
 * 1. Computes a complete diffusion trajectory for the full dataset
 * 2. Performs cross-validation by running multiple diffusion processes with different test sets
 * 3. Identifies the optimal diffusion time by finding the time step with minimum average Rf_error
 *
 * During cross-validation, test vertices are excluded from influencing other vertices (assigned
 * weight 0), but their values are still updated based on training vertices (weight 1) to evaluate
 * prediction accuracy at each time step.
 *
 * @param adj_list Adjacency list representation of the graph, where adj_list[i] contains indices of vertices
 *              adjacent to vertex i
 * @param weight_list Lengths of edges between vertices, where weight_list[i][j] is the length of the edge
 *                    between vertex i and its j-th neighbor in adj_list[i]
 * @param y Initial vertex function values to be smoothed
 * @param n_time_steps Number of diffusion time steps to perform
 * @param step_factor Factor controlling the magnitude of each diffusion step
 * @param binary_threshold Threshold for binary classification (if y is binary). If negative, uses mean of y.
 * @param ikernel Kernel function identifier for distance weighting (default: 1)
 * @param dist_normalization_factor Scaling factor applied to maximum neighbor distances (default: 1.1)
 * @param n_CVs Number of cross-validation runs to perform (default: 0, no cross-validation)
 * @param n_CV_folds Number of folds in each cross-validation run (default: 10)
 * @param epsilon Numerical tolerance for zero comparisons (default: 1e-10)
 * @param seed Random seed for test set selection (default: 0, uses current time)
 *
 * @return A graph_diffusion_smoother_result_t struct containing:
 *         - y_traj: Complete diffusion trajectory for all time steps
 *         - cv_errors: Individual CV errors for each time step and CV run
 *         - mean_cv_errors: Mean CV errors across all runs for each time step
 *         - y_optimal: Smoothed values at the optimal diffusion time
 *         - n_time_steps: Number of time steps used
 *         - n_CVs: Number of CV runs performed
 *         - optimal_time_step: Time step with minimum mean CV Rf_error
 *         - min_cv_error: Minimum mean CV Rf_error value
 *
 * @note The function detects if y is binary (containing only values 0.0 and 1.0) and
 *       uses log-likelihood for Rf_error calculation in that case. Otherwise, it uses
 *       mean absolute deviation (MAD) for continuous values.
 *
 * @note When n_CVs is 0, no cross-validation is performed and the final diffusion
 *       result (at time step n_time_steps) is returned as y_optimal.
 *
 * @note The optimal_time_step is set to -1 if no cross-validation is performed or if
 *       an Rf_error occurs during the computation of the optimal time step.
 *
 * @Rf_warning The function does not perform parameter validation. This should be done
 *          by the calling R function.
 *
 * @see graph_diffusion_smoother_result_t
 * @see initialize_kernel
 */
graph_diffusion_smoother_result_t
graph_diffusion_smoother(const std::vector<std::vector<int>>& adj_list,
                         const std::vector<std::vector<double>>& weight_list,
                         const std::vector<double>& y,
                         int n_time_steps,
                         double step_factor,
                         double binary_threshold = 0.5,
                         int ikernel = 1,
                         double dist_normalization_factor = 1.1,
                         int n_CVs = 0,
                         int n_CV_folds = 10,
                         double epsilon = 1e-10,
                         bool verbose = true,
                         unsigned int seed = 0) {

    // Add this near the beginning of the function
    auto total_start_time = std::chrono::steady_clock::now();

    // All parameter tests are done in the R function calling this one
    int n_vertices = y.size();

    if (verbose) {
        Rprintf("Starting graph diffusion smoothing with %d vertices, %d time steps, and %d CV runs\n",
                n_vertices, n_time_steps, n_CVs);
    }

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Determine the maximum number of neighbors any vertex has
    size_t max_neighbors = 0;
    for (const auto& neighbors : adj_list) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize kernel weights and vertex edge lengths vectors
    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_weight_list(max_neighbors);

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && binary_threshold < 0) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = adj_list[vertex].size();

    /**
     * @brief Performs one step of weighted kernel-based diffusion smoothing on a graph.
     *
     * This lambda function applies a single iteration of diffusion smoothing using a kernel function,
     * respecting the provided vertex weights (0 for test vertices, 1 for training vertices).
     *
     * @param y_current [in] The current values of all vertices before this diffusion step.
     * @param y_next [out] The updated values of all vertices after this diffusion step.
     * @param vertex_weights [in] The weights for each vertex (0 for test vertices, 1 for training vertices).
     */
    auto weighted_kernel_diffusion_loop = [&](std::vector<double>& y_current,
                                             std::vector<double>& y_next,
                                             const std::vector<double>& vertex_weights) {

        for (int vertex = 0; vertex < n_vertices; ++vertex) {

            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue; // Ensure no isolated vertices

            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_weight_list[j] = weight_list[vertex][j];
                if (vertex_weight_list[j] > max_dist)
                    max_dist = vertex_weight_list[j];
            }

            max_dist *= dist_normalization_factor;

            // Normalizing vertex edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j)
                vertex_weight_list[j] /= max_dist;

            kernel_fn(vertex_weight_list.data(), n_vertex_neighbors, kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;

            // Only consider neighbors with non-zero weights (training vertices)
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = adj_list[vertex][j];

                // Only include training vertices in the average calculation
                if (vertex_weights[neighbor] > epsilon) {
                    average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                    total_weight += kernel_weights[j];
                }
            }

            // If all neighbors are test vertices, maintain the current value
            if (total_weight < epsilon) {
                y_next[vertex] = y_current[vertex];
            } else {
                average_neighbor_value /= total_weight;

                // Update both training and test vertices based on diffusion
                y_next[vertex] = y_current[vertex] + step_factor * (average_neighbor_value - y_current[vertex]);
            }
        }
    };

    // Regular diffusion loop without weights for the standard trajectory
    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue; // Ensure no isolated vertices

            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_weight_list[j] = weight_list[vertex][j];
                if (vertex_weight_list[j] > max_dist)
                    max_dist = vertex_weight_list[j];
            }

            max_dist *= dist_normalization_factor;

            // Normalizing vertex edge lengths
            for (int j = 0; j < n_vertex_neighbors; ++j)
                vertex_weight_list[j] /= max_dist;

            kernel_fn(vertex_weight_list.data(), n_vertex_neighbors, kernel_weights.data());

            double average_neighbor_value = 0;
            double total_weight = 0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                int neighbor = adj_list[vertex][j];
                average_neighbor_value += kernel_weights[j] * y_current[neighbor];
                total_weight += kernel_weights[j];
            }

            average_neighbor_value /= total_weight;

            y_next[vertex] = y_current[vertex] + step_factor * (average_neighbor_value - y_current[vertex]);
        }
    };

    // The main loop: Calculate full trajectory without CV
    auto y_traj = std::vector<std::vector<double>>();
    std::vector<double> y_current = y;
    std::vector<double> y_next(n_vertices);
    int time_step = 1;

    auto trajectory_start_time = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("Computing full diffusion trajectory for %d time steps...\n", n_time_steps);
    }
    while (time_step < n_time_steps) {
        // Only print every few steps to avoid excessive output
        if (verbose && (time_step % 10 == 0 || time_step == n_time_steps - 1)) {
            Rprintf("\rStep %d of %d", time_step + 1, n_time_steps);
        }

        kernel_diffusion_loop(y_current, y_next);
        y_traj.push_back(y_next);
        std::swap(y_current, y_next);
        time_step++;
    }
    if (verbose) {
        Rprintf("\r");
        elapsed_time(trajectory_start_time, "Full diffusion trajectory computed in", true);
    }

    y_traj.push_back(y_current);

    std::vector<double> cv_errors(n_CVs * n_time_steps, 0.0); // Initialize with 0.0 when n_CVs = 0

    if (n_CVs > 0) {

        auto cv_start_time = std::chrono::steady_clock::now();
        if (verbose) {
            Rprintf("Starting cross-validation with %d runs and %d folds each\n", n_CVs, n_CV_folds);
        }

        int fold_size = n_vertices / n_CV_folds;

        // Cross-validation loop
        for (int cv = 0; cv < n_CVs; ++cv) {
            // auto cv_run_start_time = std::chrono::steady_clock::now();
            if (verbose) {
                Rprintf("\rCV run %d of %d: testing on %d vertices",
                        cv + 1, n_CVs, fold_size);
            }

            // Creating a test set
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

            // Create vertex weights: 0 for test vertices, 1 for training vertices
            std::vector<double> vertex_weights(n_vertices, 1.0);
            for (const auto& vertex : test_set) {
                vertex_weights[vertex] = 0.0;
            }

            // Start with original y values
            y_current = y;

            // Perform diffusion smoothing on the weighted graph at each step
            time_step = 0;
            while (time_step < n_time_steps) {
                // if (verbose && (time_step == 0 || time_step == n_time_steps - 1)) {
                //     Rprintf("  CV run %d: diffusion step %d of %d\n",
                //             cv + 1, time_step + 1, n_time_steps);
                // }

                // Run one step of weighted diffusion
                weighted_kernel_diffusion_loop(y_current, y_next, vertex_weights);

                // Computing CV Rf_error at the test vertices
                double cv_error = 0.0;
                if (y_binary) {
                    // Log-likelihood Rf_error for binary values
                    for (const auto& vertex : test_set) {
                        double clipped_y_next = std::clamp(y_next[vertex], epsilon, 1.0 - epsilon);
                        cv_error += y[vertex] * log(clipped_y_next) + (1 - y[vertex]) * log(1 - clipped_y_next);
                    }
                    cv_error *= -1;
                } else {
                    // MAD Rf_error for continuous values
                    for (const auto& vertex : test_set) {
                        cv_error += std::abs(y_next[vertex] - y[vertex]);
                    }
                    cv_error /= test_set.size();
                }

                cv_errors[time_step + cv * n_time_steps] = cv_error;

                // Swap y_current and y_next
                std::swap(y_current, y_next);

                time_step++;

                #if 0
                if (verbose) {
                    char message[100];
                    snprintf(message, sizeof(message), "CV run %d of %d completed in", cv + 1, n_CVs);
                    elapsed_time(cv_run_start_time, message, true);
                }
                #endif
            }
        }

        if (verbose) {
            Rprintf("\r");
            elapsed_time(cv_start_time, "All cross-validation runs completed in", true);
        }
    }

    graph_diffusion_smoother_result_t result;
    result.n_time_steps = n_time_steps;
    result.n_CVs = n_CVs;
    result.y_traj = std::move(y_traj);
    result.cv_errors = std::move(cv_errors);

    // Initialize optimal time step and minimum CV Rf_error with default values
    result.optimal_time_step = -1;  // -1 indicates no optimal time step found
    result.min_cv_error = std::numeric_limits<double>::infinity();

    // Compute mean CV errors and find optimal time step
    if (n_CVs > 0) {

        auto opt_start_time = std::chrono::steady_clock::now();
        if (verbose) {
            Rprintf("Computing optimal diffusion time from CV results\n");
        }

        try {
            // Resize mean_cv_errors vector
            result.mean_cv_errors.resize(n_time_steps, 0.0);

            // Compute mean CV errors for each time step
            for (int time_step = 0; time_step < n_time_steps; ++time_step) {
                double sum_error = 0.0;
                for (int cv = 0; cv < n_CVs; ++cv) {
                    sum_error += result.cv_errors[time_step + cv * n_time_steps];
                }
                result.mean_cv_errors[time_step] = sum_error / n_CVs;
            }

            // Find the optimal time step (index of global minimum)
            if (!result.mean_cv_errors.empty()) {
                auto min_iter = std::min_element(result.mean_cv_errors.begin(),
                                                 result.mean_cv_errors.end());

                result.optimal_time_step = std::distance(result.mean_cv_errors.begin(),
                                                          min_iter);
                result.min_cv_error = *min_iter;

                // Store the optimal solution
                if (result.optimal_time_step >= 0 &&
                    result.optimal_time_step < static_cast<int>(result.y_traj.size())) {
                    result.y_optimal = result.y_traj[result.optimal_time_step];
                } else {
                    REPORT_WARNING("\nWarning: optimal_time_step out of bounds for y_traj\n");
                    // Initialize y_optimal with original y values as fallback
                    result.y_optimal = y;
                }

                if (verbose) {
                    Rprintf("Found optimal diffusion time: step %d (Rf_error: %.6f)\n",
                            result.optimal_time_step, result.min_cv_error);
                }

            } else {
                REPORT_WARNING("\nWarning: mean_cv_errors vector is empty\n");
                // Initialize y_optimal with original y values as fallback
                result.y_optimal = y;
            }

            if (verbose) {
                elapsed_time(opt_start_time, "Optimal diffusion time computation completed in", true);
            }

        } catch (const std::exception& e) {
            REPORT_ERROR("\nComputing optimal time step: %s\n", e.what());
            // Keep default values and initialize y_optimal with original y values
            result.y_optimal = y;
        }
    } else {
        if (verbose) {
            Rprintf("No cross-validation performed, using final diffusion state\n");
        }
        // When no CV is performed, use the final smoothed values
        result.y_optimal = !result.y_traj.empty() ? result.y_traj.back() : y;
    }

    if (verbose) {
        elapsed_time(total_start_time, "Total graph diffusion processing completed in", true);
    }

    return result;
}

/**
 * @brief R interface to graph diffusion smoothing with cross-validation
 *
 * @details This function provides an R interface to the C++ implementation of graph
 * diffusion smoothing. It handles the conversion between R and C++ data structures,
 * calls the graph_diffusion_smoother function, and converts the results back to R objects.
 *
 * The function performs graph diffusion smoothing with optional cross-validation to
 * identify the optimal diffusion time. It uses a masked diffusion approach where
 * test vertices are excluded from the diffusion process by assigning them zero weights,
 * while still updating their values for Rf_error calculation.
 *
 * @param s_adj_list R list containing the adjacency list representation of the graph
 * @param s_weight_list R list containing the edge weights for the graph
 * @param s_y R numeric vector of values to be smoothed
 * @param s_n_time_steps R integer for number of diffusion time steps
 * @param s_step_factor R numeric for diffusion step factor
 * @param s_binary_threshold R numeric for binary classification threshold (if y is binary)
 * @param s_ikernel R integer specifying the kernel function for distance weighting
 * @param s_dist_normalization_factor R numeric for distance normalization factor
 * @param s_n_CVs R integer for number of cross-validation runs
 * @param s_n_CV_folds R integer for number of folds in each CV run
 * @param s_epsilon R numeric for numerical tolerance
 * @param s_seed R integer for random seed
 * @param s_verbose R logical for verbose output
 *
 * @return An R list containing:
 *  - y_traj: List of numeric vectors representing the trajectory of y values
 *  - cv_errors: Matrix of CV errors (rows=time steps, cols=CV runs)
 *  - mean_cv_errors: Numeric vector of mean CV errors for each time step
 *  - y_optimal: Numeric vector of smoothed values at the optimal time step
 *  - n_time_steps: Integer number of time steps used
 *  - n_CVs: Integer number of CV runs performed
 *  - optimal_time_step: Integer time step with minimum mean CV Rf_error
 *  - min_cv_error: Numeric minimum mean CV Rf_error value
 *
 * @note This function performs memory management using R's PROTECT/UNPROTECT mechanism
 * to prevent garbage collection of intermediate results.
 *
 * @see graph_diffusion_smoother
 */
SEXP S_graph_diffusion_smoother(SEXP s_adj_list,
                                SEXP s_weight_list,
                                SEXP s_y,
                                SEXP s_n_time_steps,
                                SEXP s_step_factor,
                                SEXP s_binary_threshold,
                                SEXP s_ikernel,
                                SEXP s_dist_normalization_factor,
                                SEXP s_n_CVs,
                                SEXP s_n_CV_folds,
                                SEXP s_epsilon,
                                SEXP s_verbose,
                                SEXP s_seed) {

    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Extract scalar parameters using R's C API
    int n_time_steps = INTEGER(s_n_time_steps)[0];
    double step_factor = REAL(s_step_factor)[0];
    double binary_threshold = REAL(s_binary_threshold)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = INTEGER(s_seed)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    int nprot = 0;

    // Call the graph diffusion smoother function
    graph_diffusion_smoother_result_t result = graph_diffusion_smoother(
        adj_list,
        weight_list,
        y,
        n_time_steps,
        step_factor,
        binary_threshold,
        ikernel,
        dist_normalization_factor,
        n_CVs,
        n_CV_folds,
        epsilon,
        verbose,
        seed
        );

    // Create the return list using R's C API
    const char* names[] = {
        "y_traj",
        "cv_errors",
        "mean_cv_errors",
        "y_optimal",
        "n_time_steps",
        "n_CVs",
        "optimal_time_step",
        "min_cv_error",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Create list and protect it
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements)); nprot++;
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_elements)); nprot++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    // Helper function to convert vector to SEXP
    auto create_numeric_vector = [&nprot](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size())); nprot++;
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Helper function to convert matrix to SEXP (column-major for R)
    auto create_numeric_matrix = [&nprot](const std::vector<double>& vec, int nrow, int ncol) -> SEXP {
        if (vec.empty() || nrow == 0 || ncol == 0) return R_NilValue;

        SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); nprot++;
        double* ptr = REAL(r_mat);

        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                ptr[i + j * nrow] = vec[i + j * nrow];  // Column-major order for R
            }
        }
        return r_mat;
    };

    // Convert y_traj to an R list of vectors
    SEXP r_y_traj = PROTECT(Rf_allocVector(VECSXP, result.y_traj.size())); nprot++;
    for (size_t t = 0; t < result.y_traj.size(); ++t) {
        SET_VECTOR_ELT(r_y_traj, t, create_numeric_vector(result.y_traj[t]));
    }
    SET_VECTOR_ELT(r_result, 0, r_y_traj);

    // Convert cv_errors to R matrix
    if (result.n_CVs > 0 && result.n_time_steps > 0) {
        SEXP r_cv_errors = create_numeric_matrix(result.cv_errors, result.n_time_steps, result.n_CVs);
        SET_VECTOR_ELT(r_result, 1, r_cv_errors);
    } else {
        SET_VECTOR_ELT(r_result, 1, R_NilValue);
    }

    // Convert mean_cv_errors to R vector
    SET_VECTOR_ELT(r_result, 2, create_numeric_vector(result.mean_cv_errors));

    // Convert y_optimal to R vector
    SET_VECTOR_ELT(r_result, 3, create_numeric_vector(result.y_optimal));

    // Convert scalar values
    SEXP r_n_time_steps = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
    INTEGER(r_n_time_steps)[0] = result.n_time_steps;
    SET_VECTOR_ELT(r_result, 4, r_n_time_steps);

    SEXP r_n_CVs = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
    INTEGER(r_n_CVs)[0] = result.n_CVs;
    SET_VECTOR_ELT(r_result, 5, r_n_CVs);

    SEXP r_optimal_time_step = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
    INTEGER(r_optimal_time_step)[0] = result.optimal_time_step;
    SET_VECTOR_ELT(r_result, 6, r_optimal_time_step);

    SEXP r_min_cv_error = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
    REAL(r_min_cv_error)[0] = result.min_cv_error;
    SET_VECTOR_ELT(r_result, 7, r_min_cv_error);

    UNPROTECT(nprot);
    return r_result;
}
