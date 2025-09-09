#include <R.h>
#include <Rinternals.h>

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

#include <omp.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <map>
#include <numeric>

#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseSymMatProd.h>

#include "msr2.h"
#include "msr2_cpp_utils.h"
#include "msr2_Cpp_to_R_utils.h"
#include "msr2_Eigen_utils.h"
#include "msr2_graph_diffussion_smoother.h"
#include "msr2_kernels.h"

void scale_to_range(std::vector<double>& x, double ymin, double ymax);
double calculate_smallest_difference(const std::vector<double>& y);
std::unique_ptr<std::pair<std::vector<int>, std::vector<int>>> find_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y);
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

SEXP S_graph_diffusion_smoother(SEXP Rgraph,
                                SEXP Rd,
                                SEXP Rweights,
                                SEXP Ry,
                                SEXP Rn_time_steps,
                                SEXP Rstep_factor,
                                SEXP Rnormalize,
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
                                SEXP Rseed);

SEXP S_graph_spectral_smoother(SEXP Rgraph,
                               SEXP Rd,
                               SEXP Rweights,
                               SEXP Ry,
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
                               SEXP Rmin_plambda,
                               SEXP Rmax_plambda,
                               SEXP Rseed);
}

// ------------------------------------------------------------------------------------------------------
//
// set version of diffusion smoother
//
//
// Used to be stop_criterion 2 of weighted_heat_smoother_trajectory
// with added CV
//
// ------------------------------------------------------------------------------------------


/**
 * @brief Estimates the function values for excluded vertices using their neighbors' values and computes the mean absolute deviation error.
 *
 * This function uses the weighted mean of the neighbors' smoothed values to estimate the function values for the excluded vertices
 * (test set). It then computes the mean absolute deviation (MAD) error between the estimated and original function values.
 *
 * @param y         The original function values.
 * @param y_t       The smoothed function values.
 * @param test_set  The set of indices of excluded vertices.
 * @param test_graph The test graph representing the neighborhood structure of excluded vertices.
 *
 * @return The mean absolute deviation prediction error for the excluded vertices.
 */
double graph_mean_nb_predict_error(const std::vector<double>& y,
                                   std::vector<double>& y_t,
                                   const std::set<int>& test_set,
                                   const std::vector<std::set<std::pair<int, int>>>& test_graph) {

    for (int excluded_vertex : test_set) {
        double sum_values = 0.0;
        for (const auto& pair : test_graph[excluded_vertex]) {
            int neighbor = pair.first;
            double weight = 1.0 / pair.second;
            sum_values += weight * y_t[neighbor];
        }
        y_t[excluded_vertex] = sum_values / test_graph[excluded_vertex].size();
    }

    // Calculate mean absolute deviation error for excluded vertices
    double error = 0.0;
    for (int excluded_vertex : test_set) {
        error += std::abs(y_t[excluded_vertex] - y[excluded_vertex]);
    }

    return error / test_set.size();
}

/**
 * @brief Creates a cross-validation (CV) graph by removing test vertices from the original graph.
 *
 * This function creates a CV graph from a set-based adjacency list of the original graph by removing the test vertices
 * and their corresponding edges. It returns a map of the original neighbor sets of the test vertices, which can be used
 * to reconstruct the original graph later.
 *
 * @param set_graph A reference to the graph represented as a vector of sets, where each set contains the neighbors of a vertex.
 * @param test_set  The set of vertices to be used as the test set for cross-validation.
 *
 * @return A unique pointer to a map, where the keys are the test vertices and the values are their corresponding neighbor sets from the original graph.
 *
 * @note The input graph `set_graph` is modified in-place by removing the test vertices and their edges. The returned map contains
 *       the original neighbor sets of the test vertices.
 */
std::unique_ptr<std::map<int, std::set<int>>>
create_cv_graph(std::vector<std::set<int>>& set_graph,
                const std::set<int>& test_set) {

    auto result = std::make_unique<std::map<int, std::set<int>>>();

    for (const auto& i : test_set)
        (*result)[i] = set_graph[i];

    for (const auto& i : test_set) {
        for (const auto& neighbor : set_graph[i]) {
            set_graph[neighbor].erase(i);
        }
        // set_graph[i].clear(); // Since when performing diffusion on a CV graph we iterate over training vertices set_graph[i] is never visited for 'i' an element of test set, so there is no need to clear it
    }

    return result;
}

/**
 * @brief Reconstructs the original graph by reinserting test vertices and their edges.
 *
 * This function reconstructs the original graph from a CV graph by reinserting the test vertices and their corresponding
 * edges using the stored neighbor sets from the `create_cv_graph` function.
 *
 * @param set_graph A reference to the CV graph represented as a vector of sets, which will be modified to represent the reconstructed original graph.
 * @param res       A unique pointer to a map, where the keys are the test vertices and the values are their corresponding neighbor sets from the original graph.
 *
 * @note The `set_graph` parameter is modified in-place to represent the reconstructed original graph. The `res` parameter
 *       is the result obtained from the `create_cv_graph` function.
 */
void reconstruct_graph(std::vector<std::set<int>>& set_graph,
                       std::unique_ptr<std::map<int, std::set<int>>>& res) {

    for (const auto& [i, i_neighbors] : (*res)) {
        // set_graph[i] = i_neighbors; // set_graph[i] was not cleared, so we don't have to reset it

        for (const auto& neighbor : i_neighbors) {
            std::set<int> union_set;
            std::set<int> i_set = {i};
            std::set_union(set_graph[neighbor].begin(), set_graph[neighbor].end(),
                           i_set.begin(), i_set.end(),
                           std::inserter(union_set, union_set.begin()));

            set_graph[neighbor] = union_set;
        }
    }
}

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
 * @param epsilon A small positive constant for numerical stability in binary error calculations.
 * @param seed A seed for the random number generator used in cross-validation.
 *
 * @return A unique pointer to a graph_diffusion_smoother_result_t structure containing:
 *         - y_traj: A vector of vectors representing the trajectory of y values over time.
 *         - cv_errors: A vector of cross-validation errors for each time step and CV round.
 *         - n_time_steps: The number of time steps used in the diffusion process.
 *         - n_CVs: The number of cross-validation rounds performed.
 *
 * @throws std::runtime_error if input parameters are invalid or inconsistent.
 *
 * @note The function supports both hop-index and distance-based neighbor calculations,
 *       determined by whether the 'd' parameter is empty or not.
 *       It also handles binary data differently, using log-likelihood for error calculation.
 *
 */

// July 26, 2024 version befor fixing the issue with updating the values sequentially and implementing the diffusiong smoothing using a lambda function
std::unique_ptr<graph_diffusion_smoother_result_t>
graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                         const std::vector<std::vector<double>>& edge_lengths,
                         const std::vector<double>& weights,
                         const std::vector<double>& y,
                         int n_time_steps,
                         double step_factor,
                         int normalize,
                         imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                         iterative_imputation_params_t iterative_params = {},
                         bool apply_binary_threshold = true,
                         double binary_threshold = 0.5,
                         int ikernel = 1,
                         double dist_normalization_factor = 1.01,
                         int n_CVs = 0,
                         int n_CV_folds = 10,
                         double epsilon = 1e-10,
                         unsigned int seed = 0) {

    // All parameter tests are done in the R function calling this one
    int n_vertices = y.size();

    if (ikernel)
        initialize_kernel(ikernel);

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
        for (int vertex = 0; vertex < n_vertices; ++vertex)
            Ey += y[vertex];
        Ey /= n_vertices;
    }

    auto y_traj = std::vector<std::vector<double>>();
    std::vector<double> y_t(y); // Initialize y_t to y

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    // computing min and max values of y
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());

    int time_step = 0;
    while (time_step < n_time_steps) {

        y_traj.push_back(y_t);

        if (ikernel == 0) {
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                double average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += y_t[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
            }
        } else {
            for (int vertex = 0; vertex < n_vertices; ++vertex) {

                int n_vertex_neighbors = n_neighbors[vertex];
                if (n_vertex_neighbors == 0) continue; // Ensure no isolated vertices

                double max_dist = 0.0;
                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    vertex_edge_lengths[j] = edge_lengths[vertex][j];
                    if (vertex_edge_lengths[j] > max_dist)
                        max_dist = vertex_edge_lengths[j];
                }

                // Normalizing vertex edge lengths
                for (int j = 0; j < n_vertex_neighbors; ++j)
                    vertex_edge_lengths[j] /= max_dist;

                kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

                double average_neighbor_value = 0;
                double total_weight = 0;
                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    int neighbor = graph[vertex][j];
                    average_neighbor_value += kernel_weights[j] * y_t[neighbor];
                    total_weight += kernel_weights[j];
                }

                average_neighbor_value /= total_weight;

                y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
            }
        }

        if (range_adjust) {
            scale_to_range(y_t, ymin, ymax);
        } else if (mean_adjust) {
            // Adjusting y_t so that its mean = mean(y)
            double Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += y_t[vertex];
            Ey_t /= n_vertices;

            double Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                y_t[vertex] += Delta;
        }

        time_step++;
    }

    y_traj.push_back(y_t);

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

            std::vector<double> y_t = std::move(*cv_imputation(test_set,
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

            // Perform diffusion smoothing on the CV graph at each step estimating the MAD error
            // NOTE: there is no range or mean adjustment in the CV loop
            int time_step = 0;
            while (time_step < n_time_steps) {

                if (ikernel == 0) {
                    for (int vertex = 0; vertex < n_vertices; ++vertex) {
                        if (n_neighbors[vertex] == 0) continue; // Ensure no isolated vertices
                        double average_neighbor_value = 0;
                        for (int neighbor : graph[vertex])
                            average_neighbor_value += y_t[neighbor];
                        average_neighbor_value /= n_neighbors[vertex];

                        y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
                    }
                } else {
                    for (int vertex = 0; vertex < n_vertices; ++vertex) {

                        if (n_neighbors[vertex] == 0) continue; // Ensure no isolated vertices

                        double max_dist = 0.0;
                        for (int j = 0; j < n_neighbors[vertex]; ++j) {
                            vertex_edge_lengths[j] = edge_lengths[vertex][j];
                            if (vertex_edge_lengths[j] > max_dist)
                                max_dist = vertex_edge_lengths[j];
                        }

                        // Normalizing vertex edge lengths
                        for (int j = 0; j < n_neighbors[vertex]; ++j)
                            vertex_edge_lengths[j] /= max_dist;

                        kernel_fn(vertex_edge_lengths.data(), n_neighbors[vertex], kernel_weights.data());

                        double average_neighbor_value = 0;
                        double total_weight = 0;
                        for (int j = 0; j < n_neighbors[vertex]; ++j) {
                            int neighbor = graph[vertex][j];
                            average_neighbor_value += kernel_weights[j] * y_t[neighbor];
                            total_weight += kernel_weights[j];
                        }

                        average_neighbor_value /= total_weight;

                        y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
                    }
                }

                // Computing CV error at the test vertices
                double cv_error = 0.0;
                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        double clipped_y_t = std::max(epsilon, std::min(1.0 - epsilon, y_t[vertex]));
                        cv_error += y[vertex] * log(clipped_y_t) + (1 - y[vertex]) * log(1 - clipped_y_t);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set)
                        cv_error += std::abs(y_t[vertex] - y[vertex]);
                    cv_error /= test_set.size();
                }
                cv_errors[time_step + cv * n_time_steps] = cv_error;

                time_step++;
            }
        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    auto result = std::make_unique<graph_diffusion_smoother_result_t>();
    result->n_time_steps = n_time_steps;
    result->n_CVs = n_CVs;
    result->y_traj = std::move(y_traj);
    result->cv_errors = std::move(cv_errors);

    return result;
}





std::unique_ptr<graph_diffusion_smoother_result_t>
graph_diffusion_smoother(const std::vector<std::vector<int>>& graph,
                         const std::vector<std::vector<double>>& edge_lengths,
                         const std::vector<double>& weights,
                         const std::vector<double>& y,
                         int n_time_steps,
                         double step_factor,
                         int normalize,
                         imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                         iterative_imputation_params_t iterative_params = {},
                         bool apply_binary_threshold = true,
                         double binary_threshold = 0.5,
                         int ikernel = 1,
                         double dist_normalization_factor = 1.01,
                         int n_CVs = 0,
                         int n_CV_folds = 10,
                         double epsilon = 1e-10,
                         unsigned int seed = 0) {

    #define DEBUG__graph_diffusion_smoother 0

    // All parameter tests are done in the R function calling this one
    int n_vertices = y.size();

    if (ikernel)
        initialize_kernel(ikernel);

#if DEBUG__graph_diffusion_smoother
    Rprintf("ikernel: %d\n", ikernel);
#endif

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
        for (int vertex = 0; vertex < n_vertices; ++vertex)
            Ey += y[vertex];
        Ey /= n_vertices;
    }

    auto y_traj = std::vector<std::vector<double>>();
    std::vector<double> y_t(y); // Initialize y_t to y

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    // computing min and max values of y
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());

    int time_step = 0;
    while (time_step < n_time_steps) {

#if DEBUG__graph_diffusion_smoother
        Rprintf("time_step: %d\n", time_step);
#endif
        y_traj.push_back(y_t);

        if (ikernel == 0) {
#if DEBUG__graph_diffusion_smoother
            Rprintf("In if (ikernel == 0)\n");
#endif
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                double average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += y_t[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
            }
#if DEBUG__graph_diffusion_smoother
            print_vect(y_t, "y_t");
            print_vect(weights, "weights");
#endif
        } else {
#if DEBUG__graph_diffusion_smoother
            Rprintf("In if (ikernel != 0)\n");
#endif
            for (int vertex = 0; vertex < n_vertices; ++vertex) {

                int n_vertex_neighbors = n_neighbors[vertex];
                if (n_vertex_neighbors == 0) continue; // Ensure no isolated vertices

                double max_dist = 0.0;
                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    vertex_edge_lengths[j] = edge_lengths[vertex][j];
                    if (vertex_edge_lengths[j] > max_dist)
                        max_dist = vertex_edge_lengths[j];
                }

                // Normalizing vertex edge lengths
                for (int j = 0; j < n_vertex_neighbors; ++j)
                    vertex_edge_lengths[j] /= max_dist;

                kernel_fn(vertex_edge_lengths.data(), n_vertex_neighbors, kernel_weights.data());

                double average_neighbor_value = 0;
                double total_weight = 0;
                for (int j = 0; j < n_vertex_neighbors; ++j) {
                    int neighbor = graph[vertex][j];
                    average_neighbor_value += kernel_weights[j] * y_t[neighbor];
                    total_weight += kernel_weights[j];
                }

                average_neighbor_value /= total_weight;

                y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
            }
        }

        if (range_adjust) {
            scale_to_range(y_t, ymin, ymax);
        } else if (mean_adjust) {
            // Adjusting y_t so that its mean = mean(y)
            double Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += y_t[vertex];
            Ey_t /= n_vertices;

            double Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                y_t[vertex] += Delta;
        }

        time_step++;
    }

    y_traj.push_back(y_t);

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

            std::vector<double> y_t = std::move(*cv_imputation(test_set,
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

            // Perform diffusion smoothing on the CV graph at each step estimating the MAD error
            // NOTE: there is no range or mean adjustment in the CV loop
            int time_step = 0;
            while (time_step < n_time_steps) {

                if (ikernel == 0) {
                    for (int vertex = 0; vertex < n_vertices; ++vertex) {
                        if (n_neighbors[vertex] == 0) continue; // Ensure no isolated vertices
                        double average_neighbor_value = 0;
                        for (int neighbor : graph[vertex])
                            average_neighbor_value += y_t[neighbor];
                        average_neighbor_value /= n_neighbors[vertex];

                        y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
                    }
                } else {
                    for (int vertex = 0; vertex < n_vertices; ++vertex) {

                        if (n_neighbors[vertex] == 0) continue; // Ensure no isolated vertices

                        double max_dist = 0.0;
                        for (int j = 0; j < n_neighbors[vertex]; ++j) {
                            vertex_edge_lengths[j] = edge_lengths[vertex][j];
                            if (vertex_edge_lengths[j] > max_dist)
                                max_dist = vertex_edge_lengths[j];
                        }

                        // Normalizing vertex edge lengths
                        for (int j = 0; j < n_neighbors[vertex]; ++j)
                            vertex_edge_lengths[j] /= max_dist;

                        kernel_fn(vertex_edge_lengths.data(), n_neighbors[vertex], kernel_weights.data());

                        double average_neighbor_value = 0;
                        double total_weight = 0;
                        for (int j = 0; j < n_neighbors[vertex]; ++j) {
                            int neighbor = graph[vertex][j];
                            average_neighbor_value += kernel_weights[j] * y_t[neighbor];
                            total_weight += kernel_weights[j];
                        }

                        average_neighbor_value /= total_weight;

                        y_t[vertex] += weights[vertex] * step_factor * (average_neighbor_value - y_t[vertex]);
                    }
                }

                // Computing CV error at the test vertices
                double cv_error = 0.0;
                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        double clipped_y_t = std::max(epsilon, std::min(1.0 - epsilon, y_t[vertex]));
                        cv_error += y[vertex] * log(clipped_y_t) + (1 - y[vertex]) * log(1 - clipped_y_t);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set)
                        cv_error += std::abs(y_t[vertex] - y[vertex]);
                    cv_error /= test_set.size();
                }
                cv_errors[time_step + cv * n_time_steps] = cv_error;

                time_step++;
            }
        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    auto result = std::make_unique<graph_diffusion_smoother_result_t>();
    result->n_time_steps = n_time_steps;
    result->n_CVs = n_CVs;
    result->y_traj = std::move(y_traj);
    result->cv_errors = std::move(cv_errors);

    return result;
}

/**
 * @brief R interface for the graph diffusion smoothing algorithm with cross-validation.
 *
 * This function serves as an interface between R and the C++ graph_diffusion_smoother function.
 * It applies a diffusion process to smooth a function defined on the vertices of a graph,
 * with options for various imputation methods, normalization, and cross-validation.
 *
 * @param Rgraph An R list representing the adjacency list of the graph (1-based indexing).
 * @param Rd An R list of numeric vectors representing distances between vertices.
 *           If NULL or empty, a hop index strategy will be used.
 * @param Rweights An R numeric vector of weights associated with each vertex for the diffusion process.
 * @param Ry An R numeric vector of initial function values defined on the vertices of the graph.
 * @param Rn_time_steps An R integer specifying the number of time steps for the diffusion process.
 * @param Rstep_factor An R numeric value specifying the step size factor for the diffusion process.
 * @param Rnormalize An R integer specifying the normalization option:
 *                   0 for no normalization, 1 for range adjustment, 2 for mean adjustment.
 * @param Rimputation_method An R integer specifying the imputation method to use during cross-validation.
 * @param Rmax_iterations The number of iterations in the iterative matching method.
 * @param Rconvergence_threshold The convergence threshold in the iterative matching method.
 * @param Rapply_binary_threshold An R logical indicating whether to apply binary thresholding for binary data.
 * @param Rbinary_threshold An R numeric value specifying the threshold for binary classification.
 * @param Rikernel An R integer specifying the kernel function to use for distance-based calculations.
 * @param Rn_CVs An R integer specifying the number of cross-validation rounds.
 * @param Rn_CV_folds An R integer specifying the number of folds in each cross-validation round.
 * @param Repsilon An R numeric value specifying a small positive constant for numerical stability.
 * @param Rseed An R integer used as a seed for the random number generator in cross-validation.
 *
 * @return An R list with the following components:
 *         - y_traj_matrix: A numeric matrix where each column represents the state of the
 *           function values at a particular time step.
 *         - cv_errors: A numeric matrix of cross-validation errors, where rows correspond
 *           to time steps and columns to CV rounds.
 *
 * @note This function converts R objects to C++ types, calls the graph_diffusion_smoother function,
 *       and then converts the results back to R objects. It uses PROTECT/UNPROTECT for proper
 *       memory management in R.
 *
 * @seealso graph_diffusion_smoother
 */
SEXP S_graph_diffusion_smoother(SEXP Rgraph,
                                SEXP Rd,
                                SEXP Rweights,
                                SEXP Ry,
                                SEXP Rn_time_steps,
                                SEXP Rstep_factor,
                                SEXP Rnormalize,
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
                                SEXP Rseed) {

    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    std::vector<std::vector<double>> d = std::move(*Rdouble_vectvect_to_Cpp_double_vectvect(Rd));
    int n_vertices = LENGTH(Ry);

    int nprot = 0;

    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    double *y = REAL(Ry);

    PROTECT(Rweights = coerceVector(Rweights, REALSXP)); nprot++;
    double *weights = REAL(Rweights);

    PROTECT(Rn_time_steps = coerceVector(Rn_time_steps, INTSXP)); nprot++;
    int n_time_steps = INTEGER(Rn_time_steps)[0];

    PROTECT(Rstep_factor = coerceVector(Rstep_factor, REALSXP)); nprot++;
    double step_factor = REAL(Rstep_factor)[0];

    PROTECT(Rnormalize = coerceVector(Rnormalize, INTSXP)); nprot++;
    int normalize = INTEGER(Rnormalize)[0];

    PROTECT(Rimputation_method = coerceVector(Rimputation_method, INTSXP)); nprot++;
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);

    int max_iterations = INTEGER(Rmax_iterations)[0];
    double convergence_threshold = REAL(Rconvergence_threshold)[0];

    PROTECT(Rapply_binary_threshold = coerceVector(Rapply_binary_threshold, LGLSXP)); nprot++;
    bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];

    PROTECT(Rbinary_threshold = coerceVector(Rbinary_threshold, REALSXP)); nprot++;
    double binary_threshold = REAL(Rbinary_threshold)[0];

    PROTECT(Rikernel = coerceVector(Rikernel, INTSXP)); nprot++;
    int ikernel = INTEGER(Rikernel)[0];

    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];

    PROTECT(Rn_CVs = coerceVector(Rn_CVs, INTSXP)); nprot++;
    int n_CVs = INTEGER(Rn_CVs)[0];

    PROTECT(Rn_CV_folds = coerceVector(Rn_CV_folds, INTSXP)); nprot++;
    int n_CV_folds = INTEGER(Rn_CV_folds)[0];

    PROTECT(Repsilon = coerceVector(Repsilon, REALSXP)); nprot++;
    double epsilon = REAL(Repsilon)[0];

    PROTECT(Rseed = coerceVector(Rseed, INTSXP)); nprot++;
    unsigned int seed = static_cast<unsigned int>(INTEGER(Rseed)[0]);

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    std::unique_ptr<graph_diffusion_smoother_result_t> res = graph_diffusion_smoother(graph,
                                                                                      d,
                                                                                      std::vector<double>(weights, weights + n_vertices),
                                                                                      std::vector<double>(y, y + n_vertices),
                                                                                      n_time_steps,
                                                                                      step_factor,
                                                                                      normalize,
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

    const auto& y_traj = res->y_traj;
    const auto& cv_errors = res->cv_errors;

    int n_iterations = y_traj.size();
    if (n_iterations != n_time_steps + 1) {
        Rprintf("n_iterations: %d\n", n_iterations);
        Rprintf("n_time_steps: %d\n", n_time_steps);
    }

    SEXP y_traj_matrix = PROTECT(allocMatrix(REALSXP, n_vertices, n_iterations)); nprot++;
    double *y_traj_ptr = REAL(y_traj_matrix);

    for (int i = 0; i < n_iterations; i++) {
        std::copy(y_traj[i].begin(), y_traj[i].end(), y_traj_ptr + i * n_vertices);
    }

    SEXP cv_errors_matrix = PROTECT(allocMatrix(REALSXP, n_time_steps, n_CVs)); nprot++;
    double *cv_errors_matrix_ptr = REAL(cv_errors_matrix);
    std::copy(cv_errors.begin(), cv_errors.end(), cv_errors_matrix_ptr);

    SEXP result_list = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(result_list, 0, y_traj_matrix);
    SET_VECTOR_ELT(result_list, 1, cv_errors_matrix);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); ++nprot;
    SET_STRING_ELT(names, 0, mkChar("y_traj_matrix"));
    SET_STRING_ELT(names, 1, mkChar("cv_errors"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result_list;
}



//
// Spectral diffusion smoothing
//


double mean_absolute_deviation(const Eigen::VectorXd& y_t,
                               const Eigen::VectorXd& y) {
    return (y_t - y).cwiseAbs().mean();
}

Eigen::VectorXd compute_y_t(const Eigen::VectorXd& evalues,
                            const Eigen::MatrixXd& evectors,
                            const Eigen::VectorXd& gft,
                            double t,
                            int n_eigenvectors) {
    int n_vertices = gft.size();
    Eigen::VectorXd y_t = Eigen::VectorXd::Zero(n_vertices);
    for (int i = 0; i < n_eigenvectors; ++i) {
        y_t += std::exp(-evalues[i] * t) * gft[i] * evectors.col(i);
    }
    return y_t;
}

double optimize_t(const Eigen::VectorXd& evalues,
                  const Eigen::MatrixXd& evectors,
                  const Eigen::VectorXd& gft,
                  const Eigen::VectorXd& y,
                  int n_eigenvectors) {
    auto objective_function = [&](double t) -> double {
        Eigen::VectorXd y_t = compute_y_t(evalues, evectors, gft, t, n_eigenvectors);
        return mean_absolute_deviation(y_t, y);
    };

    // Nelder-Mead simplex algorithm
    double t0 = 1.0;  // Initial guess for t
    double best_t = t0;
    double best_value = objective_function(t0);

    double step = 0.1;
    double tol = 1e-6;
    int max_iter = 100;

    for (int iter = 0; iter < max_iter; ++iter) {
        double t1 = best_t + step;
        double value1 = objective_function(t1);

        if (value1 < best_value) {
            best_t = t1;
            best_value = value1;
        } else {
            step /= 2.0;
        }

        if (step < tol) {
            break;
        }
    }

    return best_t;
}

/**
 * @brief Applies a graph spectral smoothing algorithm using a low-pass filter
 *        approach with cross-validation to determine the optimal number of eigenvectors.
 *
 * This function smooths a function `y` defined over the vertices of a graph by
 * projecting it onto the subspace spanned by the eigenvectors of the graph Laplacian.
 * The number of eigenvectors used is determined by cross-validation, aiming to minimize
 * the mean absolute deviation (MAD) between the original and smoothed function values.
 *
 * The algorithm can be summarized as follows:
 *
 * 1. **Input Validation**:
 *    - Ensure the lengths of `graph`, `y`, and `weights` are the same.
 *    - Ensure the proportion limits `min_plambda` and `max_plambda` are valid.
 *
 * 2. **Graph Construction**:
 *    - Construct the adjacency matrix `A` from the input graph.
 *    - Compute the degree matrix `D` and the Laplacian matrix `L = D - A`.
 *
 * 3. **Eigenvalue Decomposition**:
 *    - Use the Spectra library to compute the smallest eigenvalues and their
 *      corresponding eigenvectors of the Laplacian matrix.
 *
 * 4. **Cross-Validation**:
 *    - Split the vertices into training and test sets for `n_CVs` cross-validation rounds.
 *    - Impute initial values for the test vertices based on the training vertices using BFS.
 *    - Compute the Graph Fourier Transform (GFT) of the imputed function values.
 *    - For each number of eigenvectors from `min_num_eigenvectors` to `max_num_eigenvectors`,
 *      compute the low-pass filtered version of the function.
 *    - Calculate the MAD error at the test vertices and store the results.
 *
 * 5. **Mean CV Error Calculation**:
 *    - Compute the mean CV error for each number of eigenvectors.
 *    - Identify the number of eigenvectors that minimizes the mean CV error.
 *
 * 6. **Final Smoothing**:
 *    - Compute the low-pass filtered version of `y` using the optimal number of eigenvectors.
 *
 * 7. **Result Storage**:
 *    - Store the optimal number of eigenvectors, the smoothed function values,
 *      and the CV errors in the result structure.
 *
 * @param graph A 2D vector representing the adjacency list of the graph.
 * @param y A vector of function values defined on the vertices of the graph.
 * @param weights A vector of weights associated with the vertices.
 * @param min_plambda Lower bound on the number of eigenvectors to use for smoothing, expressed as a proportion.
 * @param max_plambda Upper bound on the number of eigenvectors to use for smoothing, expressed as a proportion.
 * @param n_CVs Number of cross-validation rounds.
 * @param n_CV_folds Number of folds in each cross-validation round.
 * @return A unique pointer to a `graph_spectral_smoother_result_t` structure containing the results of the smoothing.
 */
// This version uses dense Laplacian matrix representation
std::unique_ptr<graph_spectral_smoother_result_t>
dense_graph_spectral_smoother(std::vector<std::vector<int>>& graph,
                              const std::vector<double>& y,
                              const std::vector<double>& weights,
                              double min_plambda = 0.01, // lower bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.01 that corresponds to 1% of eigenvectors with the smallest eigenvalues
                              double max_plambda = 0.20, // upper bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.20 that corresponds to 20% of eigenvectors with the smallest eigenvalues
                              int n_CVs = 0,
                              int n_CV_folds = 10) {

#define DEBUG__graph_spectral_smoother 0

#if DEBUG__graph_spectral_smoother
    Rprintf("In graph_spectral_smoother()\n");
#endif

    int n_vertices = y.size();

    // Adjacency matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_vertices, n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        for (int neighbor : graph[vertex])
            A(vertex, neighbor) = A(neighbor, vertex) = 1;
    }

    // Compute Laplacian matrix
    Eigen::MatrixXd D = A.rowwise().sum().asDiagonal();
    Eigen::MatrixXd L = D - A;

    // Eigenvalue decomposition
    Spectra::DenseSymMatProd<double> op(L);

    int min_num_eigenvectors = static_cast<int>(min_plambda * n_vertices);
    int max_num_eigenvectors = static_cast<int>(max_plambda * n_vertices);

    if (min_num_eigenvectors <= 0)
        min_num_eigenvectors = 1;

    if (max_num_eigenvectors < min_num_eigenvectors)
        max_num_eigenvectors = min_num_eigenvectors;

    // Ensure ncv is within bounds: nev < ncv <= n
    int nev = std::min(max_num_eigenvectors + 1, n_vertices);
    int ncv = std::min(2 * nev, n_vertices); // Adjust ncv to be within bounds
    // Ensure nev < ncv
    if (nev >= ncv) {
        nev = ncv - 1;
    }

    if (nev < 1) {
        nev = 1;
        //throw std::runtime_error("Invalid value for nev. It must be at least 1.");
    }

    if (max_num_eigenvectors > nev) max_num_eigenvectors = nev;

#if DEBUG__graph_spectral_smoother
    Rprintf("Checking nev < ncv <= n_vertices\n");
    Rprintf("n_vertices: %d\n", n_vertices);
    Rprintf("nev: %d\n", nev);
    Rprintf("ncv: %d\n", ncv);
    Rprintf("min_num_eigenvectors: %d\n", min_num_eigenvectors);
    Rprintf("max_num_eigenvectors: %d\n", max_num_eigenvectors);
    //error("DEBUGGING");
#endif

    // Construct eigensolver object to find eigenvalues closest to 0
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    eigs.compute(Spectra::SortRule::SmallestAlge);

    if (eigs.info() != Spectra::CompInfo::Successful)
        throw std::runtime_error("Eigenvalue estimation with Spectra failed.");

    Eigen::VectorXd evalues = eigs.eigenvalues();
    Eigen::MatrixXd evectors = eigs.eigenvectors();

    int n_filters = max_num_eigenvectors - min_num_eigenvectors + 1;
    Eigen::MatrixXd cv_errors(n_filters, n_CVs);
    Eigen::MatrixXd low_pass_ys(n_vertices, n_filters);

#if DEBUG__graph_spectral_smoother
    print_Eigen_VectorXd(evalues, "evalues");
    Rprintf("\nn_filters: %d\n", n_filters);
    Rprintf("n_CVs: %d\n", n_CVs);
    //Rprintf("n_CVs * n_filters: %d\n\n", n_CVs * n_filters);
    Rprintf("nrows(evectors): %d\n", (int)evectors.rows());
    Rprintf("ncols(evectors): %d\n", (int)evectors.cols());
    Rprintf("nrows(cv_errors): %d\n", (int)cv_errors.rows());
    Rprintf("ncols(cv_errors): %d\n", (int)cv_errors.cols());
    Rprintf("nrows(low_pass_ys): %d\n", (int)low_pass_ys.rows());
    Rprintf("ncols(low_pass_ys): %d\n", (int)low_pass_ys.cols());
    Rprintf("length(evalues): %d\n", (int)evalues.size());
    // error("DEBUGGING");
#endif

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (n_CVs > 0) {

        double epsilon = 1e-15;

        // Creating a set version of the adjacency matrix of the graph
        std::vector<std::set<int>> set_graph(n_vertices);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
        }

        int fold_size = n_vertices / n_CV_folds;
        std::vector<int> all_vertices(n_vertices);
        std::iota(all_vertices.begin(), all_vertices.end(), 0);

        for (int cv = 0; cv < n_CVs; ++cv) {

#if DEBUG__graph_spectral_smoother
            Rprintf("cv: %d\n", cv);
#endif
            // Creating a test set
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    int vertex = rand() % n_vertices;
                    test_set.insert(vertex);
                }
            }

#if DEBUG__graph_spectral_smoother
            print_set(test_set, "test_set");
#endif

            // Creating a training set as the complement of the test set
            std::set<int> training_set;
            std::set_difference(all_vertices.begin(), all_vertices.end(), test_set.begin(), test_set.end(), std::inserter(training_set, training_set.begin()));

            // Creating a test_graph for initial imputation of y values over
            // test vertices given y values at the training vertices
            std::vector<std::set<std::pair<int, int>>> test_graph(n_vertices);

            for (const auto& test_vertex : test_set) {

                std::queue<std::pair<int, int>> bfs_queue;
                std::unordered_set<int> visited;

                bfs_queue.push({test_vertex, 1});
                visited.insert(test_vertex);

                while (!bfs_queue.empty()) {
                    auto [curr_vertex, hop_index] = bfs_queue.front();
                    bfs_queue.pop();

                    bool found_training_neighbor = false;
                    for (const auto& j : graph[curr_vertex]) {
                        if (training_set.find(j) != training_set.end()) {
                            test_graph[test_vertex].insert({j, hop_index});
                            found_training_neighbor = true;

                        } else if (visited.find(j) == visited.end()) {
                            bfs_queue.push({j, hop_index + 1});
                            visited.insert(j);
                        }
                    }

                    if (found_training_neighbor) {
                        break;
                    }
                }
            }

            // computing a CV version of y where the value of y over the test
            // vertices is derived from the values of y over the training vertices
            std::vector<double> cv_y(y); // Initializing cv_y to y

            // Impute y values at the test vertices using the means of the
            // function values at the training vertices
            for (const auto& vertex : test_set) {
                double mean_y_value = 0.0;
                double total_weight = 0.0;
                for (const auto& neighbor_pair : test_graph[vertex]) {
                    int neighbor_index = neighbor_pair.first;
                    double neighbor_weight = 1.0 / neighbor_pair.second;
                    mean_y_value += neighbor_weight * cv_y[neighbor_index];
                    total_weight += neighbor_weight;
                }
                mean_y_value /= total_weight;
                if (y_binary) {
                    cv_y[vertex] = (mean_y_value > 0.5) ? 1 : 0;
                } else {
                    cv_y[vertex] = mean_y_value;
                }
            }

            // Turn cv_y into an Eigen vector
            Eigen::VectorXd cv_y_evect = Eigen::Map<const Eigen::VectorXd>(cv_y.data(), cv_y.size());

            // Compute the Graph Fourier Transform (GFT) of cv_y
            Eigen::VectorXd gft_cv_y = evectors.transpose() * cv_y_evect;

            // print_Eigen_VectorXd(gft_cv_y, "gft_cv_y");

            // for each n_eigenvectors construct a low-pass filtered version of cv_y using the first 'n_eigenvectors' eigenvectors
            for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors;
                 n_eigenvectors <= max_num_eigenvectors;
                 ++n_eigenvectors, ++filter_index) {

                // assert(filter_index < n_filters);


                Eigen::VectorXd low_pass_cv_y = Eigen::VectorXd::Zero(n_vertices);
                //for (int i = 0; i < n_eigenvectors; ++i)
                for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
                    low_pass_cv_y += gft_cv_y[i] * evectors.col(i);

                // Computing MAD error at the test vertices
                double cv_error = 0.0;

                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        // cv_error += y[vertex] * log(low_pass_cv_y[vertex]) + (1 - y[vertex]) * log(1 - low_pass_cv_y[vertex]);
                        double clipped_low_pass_cv_y = std::max(epsilon, std::min(1.0 - epsilon, low_pass_cv_y[vertex]));
                        cv_error += y[vertex] * log(clipped_low_pass_cv_y) + (1 - y[vertex]) * log(1 - clipped_low_pass_cv_y);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set) {
                        cv_error += std::abs(low_pass_cv_y[vertex] - y[vertex]);
                        //cv_error += std::abs(low_pass_cv_y[vertex] - y[vertex]);
                    }
                    cv_error /= test_set.size();
                }
                //cv_errors[filter_index + cv * n_filters] = cv_error;
                cv_errors(filter_index, cv) = cv_error;
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    // computing row-wise means of the matrix cv_errors
    std::vector<double> mean_cv_errors(n_filters, 0.0);
    for (int filter_index = 0; filter_index < n_filters; ++filter_index) {
        for (int cv = 0; cv < n_CVs; ++cv) {
            mean_cv_errors[filter_index] += cv_errors(filter_index, cv); //cv_errors[filter_index + cv * n_filters];
        }
        mean_cv_errors[filter_index] /= n_CVs;
    }

    // finding the low-pass filter index with the smallest mean CV error
    double min_cv_error = mean_cv_errors[0];
    int opt_filter_index = 0;
    for (int filter_index = 1; filter_index < n_filters; ++filter_index) {
        if (mean_cv_errors[filter_index] < min_cv_error) {
            min_cv_error = mean_cv_errors[filter_index];
            opt_filter_index = filter_index;
        }
    }

    // computing the low-pass filtered version of y for the given opt_filter_index
    // Turning y into an Eigen vector
    Eigen::VectorXd y_evect = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());

    // Computing the Graph Fourier Transform (GFT) of y
    Eigen::VectorXd gft_y = evectors.transpose() * y_evect;

    Eigen::VectorXd opt_low_pass_y = Eigen::VectorXd::Zero(n_vertices);
    int n_eigenvectors = min_num_eigenvectors + opt_filter_index;
    for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
        opt_low_pass_y += gft_y[i] * evectors.col(i);

    for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors; n_eigenvectors < max_num_eigenvectors; ++n_eigenvectors, ++filter_index) {
        Eigen::VectorXd low_pass_y = Eigen::VectorXd::Zero(n_vertices);
        for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
            low_pass_y += gft_y[i] * evectors.col(i);

        for (int vertex = 0; vertex < n_vertices; ++vertex)
            low_pass_ys(vertex, filter_index) = low_pass_y[vertex];
    }

    auto result = std::make_unique<graph_spectral_smoother_result_t>();

    result->optimal_num_eigenvectors = n_eigenvectors;
    result->y_smoothed.assign(opt_low_pass_y.data(), opt_low_pass_y.data() + opt_low_pass_y.size());
    result->cv_errors = cv_errors;
    result->mean_cv_errors = mean_cv_errors;
    result->low_pass_ys = low_pass_ys;
    result->n_filters = n_filters;
    result->evalues = evalues;
    result->evectors = evectors;
    result->min_num_eigenvectors = min_num_eigenvectors;
    result->max_num_eigenvectors = max_num_eigenvectors;

    return result;
}

/**
 * @brief Applies a graph spectral smoothing algorithm using a low-pass filter
 *        approach with cross-validation to determine the optimal number of eigenvectors.
 *
 * This function smooths a function `y` defined over the vertices of a graph by
 * projecting it onto the subspace spanned by the eigenvectors of the graph Laplacian.
 * It uses sparse Laplacian matrix representation and supports various imputation methods
 * for cross-validation.
 *
 * The algorithm can be summarized as follows:
 *
 * 1. **Input Preparation**:
 *    - Determine if hop index or kernel distance strategy should be used based on `d`.
 *    - Check if `y` is binary and set the binary threshold if necessary.
 *
 * 2. **Graph Construction**:
 *    - Construct the adjacency matrix `A` from the input graph.
 *    - Compute the degree matrix `D` and the Laplacian matrix `L = D - A`.
 *
 * 3. **Eigenvalue Decomposition**:
 *    - Use the Spectra library to compute the smallest eigenvalues and their
 *      corresponding eigenvectors of the Laplacian matrix.
 *
 * 4. **Cross-Validation**:
 *    - Split the vertices into training and test sets for `n_CVs` cross-validation rounds.
 *    - Create a test graph for imputation based on the chosen strategy (hop index or kernel distance).
 *    - Impute values for test vertices using the specified imputation method.
 *    - Apply binary thresholding if necessary.
 *    - Compute the Graph Fourier Transform (GFT) of the imputed function values.
 *    - For each number of eigenvectors, compute the low-pass filtered version and calculate errors.
 *
 * 5. **Optimal Filter Selection**:
 *    - Compute the mean CV error for each number of eigenvectors.
 *    - Identify the number of eigenvectors that minimizes the mean CV error.
 *
 * 6. **Final Smoothing**:
 *    - Compute the low-pass filtered version of `y` using the optimal number of eigenvectors.
 *
 * 7. **Result Storage**:
 *    - Store the optimal number of eigenvectors, the smoothed function values,
 *      CV errors, and other relevant information in the result structure.
 *
 * @param graph The adjacency list of the graph.
 * @param d Distances between vertices. If empty, hop index strategy is used.
 * @param weights Weights associated with the vertices (currently unused).
 * @param y Function values defined on the vertices of the graph.
 * @param imputation_method Method for imputing values during cross-validation.
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param apply_binary_threshold Whether to apply binary thresholding.
 * @param binary_threshold Threshold for binary classification.
 * @param ikernel Integer specifying the kernel function to use.
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization
 * range. Default value is 1.01.
 * @param n_CVs Number of cross-validation rounds.
 * @param n_CV_folds Number of folds in each cross-validation round.
 * @param epsilon Small constant for numerical stability in binary case.
 * @param use_low_pass_filter Whether to use low-pass filtering (currently unused).
 * @param preserve_local_extrema Whether to preserve local extrema (currently unused).
 * @param min_plambda Lower bound on the proportion of eigenvectors to use.
 * @param max_plambda Upper bound on the proportion of eigenvectors to use.
 * @param seed Seed for random number generation.
 *
 * @return A unique pointer to a `graph_spectral_smoother_result_t` structure
 *         containing the results of the smoothing.
 *
 * @note The function supports both binary and continuous data, with special
 *       handling for binary cases including different error calculations.
 */
std::unique_ptr<graph_spectral_smoother_result_t>
graph_spectral_smoother(const std::vector<std::vector<int>>& graph,
                        const std::vector<std::vector<double>>& d,
                        const std::vector<double>& weights,
                        const std::vector<double>& y,
                        imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                        iterative_imputation_params_t iterative_params = {},
                        bool apply_binary_threshold = true,
                        double binary_threshold = 0.5,
                        int ikernel = 1,
                        double dist_normalization_factor = 1.01,
                        int n_CVs = 0,
                        int n_CV_folds = 10,
                        double epsilon = 1e-10,
                        double min_plambda = 0.01, // lower bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.01 that corresponds to 1% of eigenvectors with the smallest eigenvalues
                        double max_plambda = 0.30, // upper bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.20 that corresponds to 20% of eigenvectors with the smallest eigenvalues
                        unsigned int seed = 0) {

#define DEBUG__graph_spectral_smoother 0

#if DEBUG__graph_spectral_smoother
    Rprintf("In graph_spectral_smoother()\n");
#endif

    int n_vertices = y.size();

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

    // Adjacency matrix as a sparse matrix
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> tripletList;
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        for (int neighbor : graph[vertex]) {
            if (vertex < neighbor) {  // Ensure each edge is added only once
                tripletList.push_back(Eigen::Triplet<double>(vertex, neighbor, 1.0));
                tripletList.push_back(Eigen::Triplet<double>(neighbor, vertex, 1.0));
            }
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // Compute Laplacian matrix as a sparse matrix
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            sum += it.value();
        }
        D.insert(k, k) = sum;
    }
    Eigen::SparseMatrix<double> L = D - A;

    // Eigenvalue decomposition
    Spectra::SparseSymMatProd<double> op(L);

    int min_num_eigenvectors = static_cast<int>(min_plambda * n_vertices);
    int max_num_eigenvectors = static_cast<int>(max_plambda * n_vertices);

    if (min_num_eigenvectors <= 0)
        min_num_eigenvectors = 1;

    if (max_num_eigenvectors < min_num_eigenvectors)
        max_num_eigenvectors = min_num_eigenvectors;

    // Ensure ncv is within bounds: nev < ncv <= n
    int nev = std::min(max_num_eigenvectors + 1, n_vertices);
    int ncv = std::min(2 * nev, n_vertices); // Adjust ncv to be within bounds
    // Ensure nev < ncv
    if (nev >= ncv) {
        nev = ncv - 1;
    }

    if (nev < 1) {
        nev = 1;
    }

#if DEBUG__graph_spectral_smoother
    Rprintf("n_vertices: %d\n", n_vertices);
    Rprintf("nev: %d\n", nev);
    Rprintf("ncv: %d\n", ncv);
#endif

    // Construct eigen solver object to find eigenvalues closest to 0
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    eigs.compute(Spectra::SortRule::SmallestAlge);

    if (eigs.info() != Spectra::CompInfo::Successful)
        throw std::runtime_error("Eigenvalue estimation with Spectra failed.");

    Eigen::VectorXd evalues = eigs.eigenvalues();
    Eigen::MatrixXd evectors = eigs.eigenvectors();

    int n_filters = max_num_eigenvectors - min_num_eigenvectors + 1;
    Eigen::MatrixXd cv_errors(n_filters, n_CVs);
    Eigen::MatrixXd low_pass_ys(n_vertices, n_filters);

#if DEBUG__graph_spectral_smoother
    print_Eigen_VectorXd(evalues, "evalues");
    Rprintf("\nn_filters: %d\n", n_filters);
    Rprintf("n_CVs: %d\n", n_CVs);
    //Rprintf("n_CVs * n_filters: %d\n\n", n_CVs * n_filters);
    Rprintf("nrows(evectors): %d\n", (int)evectors.rows());
    Rprintf("ncols(evectors): %d\n", (int)evectors.cols());
    Rprintf("nrows(cv_errors): %d\n", (int)cv_errors.rows());
    Rprintf("ncols(cv_errors): %d\n", (int)cv_errors.cols());
    Rprintf("nrows(low_pass_ys): %d\n", (int)low_pass_ys.rows());
    Rprintf("ncols(low_pass_ys): %d\n", (int)low_pass_ys.cols());
#endif

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

#if DEBUG__graph_spectral_smoother
            Rprintf("cv: %d\n", cv);
#endif
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

#if DEBUG__graph_spectral_smoother
            print_set(test_set, "test_set");
#endif
            std::vector<double> cv_y = std::move(*cv_imputation(test_set,
                                                                graph,
                                                                d,
                                                                y,
                                                                y_binary,
                                                                imputation_method,
                                                                iterative_params,
                                                                apply_binary_threshold,
                                                                binary_threshold,
                                                                ikernel,
                                                                dist_normalization_factor));

            // Turn cv_y into an Eigen vector
            Eigen::VectorXd cv_y_evect = Eigen::Map<const Eigen::VectorXd>(cv_y.data(), cv_y.size());

            // Compute the Graph Fourier Transform (GFT) of cv_y
            Eigen::VectorXd gft_cv_y = evectors.transpose() * cv_y_evect;

            // print_Eigen_VectorXd(gft_cv_y, "gft_cv_y");

            // for each n_eigenvectors construct a low-pass filtered version of cv_y using the first 'n_eigenvectors' eigenvectors
            for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors;
                 n_eigenvectors <= max_num_eigenvectors;
                 ++n_eigenvectors, ++filter_index) {

                Eigen::VectorXd low_pass_cv_y = Eigen::VectorXd::Zero(n_vertices);
                for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
                    low_pass_cv_y += gft_cv_y[i] * evectors.col(i);

                // Computing MAD error at the test vertices
                double cv_error = 0.0;

                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        double clipped_low_pass_cv_y = std::max(epsilon, std::min(1.0 - epsilon, low_pass_cv_y[vertex]));
                        cv_error += y[vertex] * log(clipped_low_pass_cv_y) + (1 - y[vertex]) * log(1 - clipped_low_pass_cv_y);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set) {
                        cv_error += std::abs(low_pass_cv_y[vertex] - y[vertex]);
                    }
                    cv_error /= test_set.size();
                }
                cv_errors(filter_index, cv) = cv_error;
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    // computing row-wise means of the matrix cv_errors
    std::vector<double> mean_cv_errors(n_filters, 0.0);
    for (int filter_index = 0; filter_index < n_filters; ++filter_index) {
        for (int cv = 0; cv < n_CVs; ++cv) {
            mean_cv_errors[filter_index] += cv_errors(filter_index, cv); //cv_errors[filter_index + cv * n_filters];
        }
        mean_cv_errors[filter_index] /= n_CVs;
    }

    // finding the low-pass filter index with the smallest mean CV error
    double min_cv_error = mean_cv_errors[0];
    int opt_filter_index = 0;
    for (int filter_index = 1; filter_index < n_filters; ++filter_index) {
        if (mean_cv_errors[filter_index] < min_cv_error) {
            min_cv_error = mean_cv_errors[filter_index];
            opt_filter_index = filter_index;
        }
    }

    // computing the low-pass filtered version of y for the given opt_filter_index
    // Turning y into an Eigen vector
    Eigen::VectorXd y_evect = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());

    // Computing the Graph Fourier Transform (GFT) of y
    Eigen::VectorXd gft_y = evectors.transpose() * y_evect;

    Eigen::VectorXd opt_low_pass_y = Eigen::VectorXd::Zero(n_vertices);
    int n_eigenvectors = min_num_eigenvectors + opt_filter_index;
    for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
        opt_low_pass_y += gft_y[i] * evectors.col(i);

    for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors; n_eigenvectors <= max_num_eigenvectors; ++n_eigenvectors, ++filter_index) {
        Eigen::VectorXd low_pass_y = Eigen::VectorXd::Zero(n_vertices);
        for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
            low_pass_y += gft_y[i] * evectors.col(i);

        for (int vertex = 0; vertex < n_vertices; ++vertex)
            low_pass_ys(vertex, filter_index) = low_pass_y[vertex];
    }

    auto result = std::make_unique<graph_spectral_smoother_result_t>();

    result->optimal_num_eigenvectors = n_eigenvectors;
    result->y_smoothed.assign(opt_low_pass_y.data(), opt_low_pass_y.data() + opt_low_pass_y.size());
    result->cv_errors = cv_errors;
    result->mean_cv_errors = mean_cv_errors;
    result->low_pass_ys = low_pass_ys;
    result->n_filters = n_filters;
    result->evalues = evalues;
    result->evectors = evectors;
    result->min_num_eigenvectors = min_num_eigenvectors;
    result->max_num_eigenvectors = max_num_eigenvectors;

    return result;
}

/**
 * @brief R interface for the graph spectral smoothing algorithm using a low-pass filter approach with cross-validation.
 *
 * This function serves as an interface between R and the C++ graph_spectral_smoother function.
 * It applies a spectral smoothing algorithm to a function defined on the vertices of a graph,
 * using cross-validation to determine the optimal number of eigenvectors for the low-pass filter.
 *
 * The algorithm includes the following main steps:
 * 1. Graph construction and Laplacian computation
 * 2. Eigenvalue decomposition of the Laplacian
 * 3. Cross-validation for determining optimal smoothing parameters
 * 4. Final smoothing of the input function
 *
 * @param Rgraph An R list of integer vectors representing the graph structure (1-based indexing).
 * @param Rd An R list of numeric vectors representing distances between vertices.
 *           If empty, a hop index strategy will be used for imputation.
 * @param Rweights An R numeric vector of weights associated with the vertices (currently unused).
 * @param Ry An R numeric vector of function values defined on the vertices of the graph.
 * @param Rimputation_method An R integer specifying the imputation method to use.
 * @param Rmax_iterations The number of iterations in the iterative matching method.
 * @param Rconvergence_threshold The convergence threshold in the iterative matching method.
 * @param Rapply_binary_threshold An R logical indicating whether to apply binary thresholding.
 * @param Rbinary_threshold An R numeric value specifying the threshold for binary classification.
 * @param Rikernel An R integer specifying the kernel function to use.
 * @param Rn_CVs An R integer specifying the number of cross-validation rounds.
 * @param Rn_CV_folds An R integer specifying the number of folds in each cross-validation round.
 * @param Repsilon An R numeric value specifying a small constant for numerical stability.
 * @param Rmin_plambda An R numeric value representing the lower bound on the proportion of eigenvectors to use.
 * @param Rmax_plambda An R numeric value representing the upper bound on the proportion of eigenvectors to use.
 * @param Rseed An R integer used as a seed for the random number generator.
 *
 * @return An R list with the following components:
 *         - evalues: Eigenvalues of the graph Laplacian
 *         - evectors: Eigenvectors of the graph Laplacian
 *         - optimal_num_eigenvectors: The optimal number of eigenvectors determined by cross-validation
 *         - y_smoothed: The smoothed function values
 *         - cv_errors: Matrix of cross-validation errors for each number of eigenvectors and CV round
 *         - mean_cv_errors: Vector of mean cross-validation errors for each number of eigenvectors
 *         - low_pass_ys: Matrix of low-pass filtered versions of y for different numbers of eigenvectors
 *         - min_num_eigenvectors: Minimum number of eigenvectors used
 *         - max_num_eigenvectors: Maximum number of eigenvectors used
 *
 * @note This function converts R objects to C++ types, calls the graph_spectral_smoother function,
 *       and then converts the results back to R objects. It uses PROTECT/UNPROTECT for proper
 *       memory management in R.
 */
SEXP S_graph_spectral_smoother(SEXP Rgraph,
                               SEXP Rd,
                               SEXP Rweights,
                               SEXP Ry,
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
                               SEXP Rmin_plambda,
                               SEXP Rmax_plambda,
                               SEXP Rseed) {

    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    std::vector<std::vector<double>> d = std::move(*Rdouble_vectvect_to_Cpp_double_vectvect(Rd));
    std::vector<double> weights = std::move(*Rvect_to_CppVect_double(Rweights));
    std::vector<double> y = std::move(*Rvect_to_CppVect_double(Ry));

    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);
    int max_iterations = INTEGER(Rmax_iterations)[0];
    double convergence_threshold = REAL(Rconvergence_threshold)[0];
    bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];
    double binary_threshold = REAL(Rbinary_threshold)[0];
    int ikernel = INTEGER(Rikernel)[0];
    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];
    int n_CVs = INTEGER(Rn_CVs)[0];
    int n_CV_folds = INTEGER(Rn_CV_folds)[0];
    double epsilon = REAL(Repsilon)[0];
    double min_plambda = REAL(Rmin_plambda)[0];
    double max_plambda = REAL(Rmax_plambda)[0];
    unsigned int seed = static_cast<unsigned int>(INTEGER(Rseed)[0]);

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Call the graph_spectral_smoother function
    std::unique_ptr<graph_spectral_smoother_result_t> result = graph_spectral_smoother(graph,
                                                                                       d,
                                                                                       weights,
                                                                                       y,
                                                                                       imputation_method,
                                                                                       iterative_params,
                                                                                       apply_binary_threshold,
                                                                                       binary_threshold,
                                                                                       ikernel,
                                                                                       dist_normalization_factor,
                                                                                       n_CVs,
                                                                                       n_CV_folds,
                                                                                       epsilon,
                                                                                       min_plambda,
                                                                                       max_plambda,
                                                                                       seed);

    // Convert the results back to R types (this part remains largely the same)
    int nprot = 0;
    SEXP Rresult = PROTECT(allocVector(VECSXP, 9));  nprot++;

        // evalues
    int evalues_length = result->evalues.size();
    SEXP Revalues = PROTECT(allocVector(REALSXP, evalues_length)); nprot++;
    for (int i = 0; i < evalues_length; ++i) {
        REAL(Revalues)[i] = result->evalues[i];
    }
    SET_VECTOR_ELT(Rresult, 0, Revalues);

    // evectors
    int evectors_rows = result->evectors.rows();
    int evectors_cols = result->evectors.cols();
    SEXP Revectors = PROTECT(allocMatrix(REALSXP, evectors_rows, evectors_cols)); nprot++;
    for (int i = 0; i < evectors_rows; ++i) {
        for (int j = 0; j < evectors_cols; ++j) {
            REAL(Revectors)[i + evectors_rows * j] = result->evectors(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 1, Revectors);

    // optimal_num_eigenvectors
    SEXP Roptimal_num_eigenvectors = PROTECT(allocVector(INTSXP, 1)); nprot++;
    INTEGER(Roptimal_num_eigenvectors)[0] = result->optimal_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 2, Roptimal_num_eigenvectors);

    // y_smoothed
    int y_smoothed_length = result->y_smoothed.size();
    SEXP Ry_smoothed = PROTECT(allocVector(REALSXP, y_smoothed_length)); nprot++;
    for (int i = 0; i < y_smoothed_length; ++i) {
        REAL(Ry_smoothed)[i] = result->y_smoothed[i];
    }
    SET_VECTOR_ELT(Rresult, 3, Ry_smoothed);

    // cv_errors matrix
    int cv_errors_rows = result->cv_errors.rows();
    int cv_errors_cols = result->cv_errors.cols();
    SEXP Rcv_errors = PROTECT(allocMatrix(REALSXP, cv_errors_rows, cv_errors_cols)); nprot++;
    for (int i = 0; i < cv_errors_rows; ++i) {
        for (int j = 0; j < cv_errors_cols; ++j) {
            REAL(Rcv_errors)[i + cv_errors_rows * j] = result->cv_errors(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 4, Rcv_errors);

    // mean_cv_errors
    int mean_cv_errors_length = result->mean_cv_errors.size();
    SEXP Rmean_cv_errors = PROTECT(allocVector(REALSXP, mean_cv_errors_length)); nprot++;
    for (int i = 0; i < mean_cv_errors_length; ++i) {
        REAL(Rmean_cv_errors)[i] = result->mean_cv_errors[i];
    }
    SET_VECTOR_ELT(Rresult, 5, Rmean_cv_errors);

    // low_pass_ys matrix
    int low_pass_ys_rows = result->low_pass_ys.rows();
    int low_pass_ys_cols = result->low_pass_ys.cols();
    SEXP Rlow_pass_ys = PROTECT(allocMatrix(REALSXP, low_pass_ys_rows, low_pass_ys_cols)); nprot++;
    for (int i = 0; i < low_pass_ys_rows; ++i) {
        for (int j = 0; j < low_pass_ys_cols; ++j) {
            REAL(Rlow_pass_ys)[i + low_pass_ys_rows * j] = result->low_pass_ys(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 6, Rlow_pass_ys);

    SEXP Rmin_num_eigenvectors = PROTECT(allocVector(REALSXP, 1)); nprot++;
    REAL(Rmin_num_eigenvectors)[0] = result->min_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 7, Rmin_num_eigenvectors);

    SEXP Rmax_num_eigenvectors = PROTECT(allocVector(REALSXP, 1)); nprot++;
    REAL(Rmax_num_eigenvectors)[0] = result->max_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 8, Rmax_num_eigenvectors);

    // Set the names of the list components
    SEXP names = PROTECT(allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("evalues"));
    SET_STRING_ELT(names, 1, mkChar("evectors"));
    SET_STRING_ELT(names, 2, mkChar("optimal_num_eigenvectors"));
    SET_STRING_ELT(names, 3, mkChar("y_smoothed"));
    SET_STRING_ELT(names, 4, mkChar("cv_errors"));
    SET_STRING_ELT(names, 5, mkChar("mean_cv_errors"));
    SET_STRING_ELT(names, 6, mkChar("low_pass_ys"));
    SET_STRING_ELT(names, 7, mkChar("min_num_eigenvectors"));
    SET_STRING_ELT(names, 8, mkChar("max_num_eigenvectors"));
    setAttrib(Rresult, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return Rresult;
}
