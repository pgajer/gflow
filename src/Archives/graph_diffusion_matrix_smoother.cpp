#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

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
#include <random>
#include <chrono>

#include "msr2.h"
#include "cpp_utils.h"
#include "SEXP_cpp_conversion_utils.h"
#include "graph_diffusion_smoother.h"
#include "kernels.h"

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
    SEXP S_graph_diffusion_matrix_smoother(SEXP r_X,
                                           SEXP r_graph,
                                           SEXP r_edge_lengths,
                                           SEXP r_weights,
                                           SEXP r_n_time_steps,
                                           SEXP r_step_factor,
                                           SEXP r_normalize,
                                           SEXP r_imputation_method,
                                           SEXP r_max_iterations,
                                           SEXP r_convergence_threshold,
                                           SEXP r_ikernel,
                                           SEXP r_dist_normalization_factor,
                                           SEXP r_n_CVs,
                                           SEXP r_n_CV_folds,
                                           SEXP r_epsilon,
                                           SEXP r_seed);
}


/**
 * @brief Performs graph diffusion smoothing on a matrix of features.
 *
 * This function applies graph diffusion smoothing to each column of the input matrix X independently.
 * It creates a trajectory of smoothed versions of X, one for each time step, and computes
 * cross-validation errors at each step.
 *
 * @param X [in] The input matrix to be smoothed, where rows represent vertices and columns represent features.
 * @param graph [in] The adjacency list representation of the graph.
 * @param edge_lengths [in] The lengths of edges between vertices.
 * @param weights [in] The individual weights for each vertex.
 * @param n_time_steps [in] The number of diffusion steps to perform.
 * @param step_factor [in] The factor controlling the magnitude of each diffusion step.
 * @param normalize [in] The normalization method to use: 0 for none, 1 for range adjustment, 2 for mean adjustment.
 * @param imputation_method [in] The method to use for imputing missing values during cross-validation.
 * @param iterative_params [in] Parameters for iterative imputation methods.
 * @param ikernel [in] The index of the kernel function to use.
 * @param dist_normalization_factor [in] The factor used for normalizing distances in the kernel function.
 * @param n_CVs [in] The number of cross-validation runs to perform.
 * @param n_CV_folds [in] The number of folds to use in each cross-validation run.
 * @param epsilon [in] A small value used to prevent numerical instability.
 * @param seed [in] The seed for the random number generator used in cross-validation.
 *
 * @return A unique pointer to a graph_diffusion_matrix_smoother_result_t struct containing:
 *         - X_traj: A trajectory of diffusion smoothing versions of X.
 *         - mean_cv_error: A vector of mean cross-validation errors, one per time step.
 *
 * @note This function assumes that all input parameters are valid and properly initialized.
 * @note The function removes all code related to binary data handling from the original graph_diffusion_smoother.
 * @note Cross-validation errors are computed as the mean of errors across all features.
 *
 * @see graph_diffusion_smoother_result_t
 * @see imputation_method_t
 * @see iterative_imputation_params_t
 */
std::unique_ptr<graph_diffusion_matrix_smoother_result_t>
graph_diffusion_matrix_smoother(
    const std::vector<std::vector<double>>& X,
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& weights,
    int n_time_steps,
    double step_factor,
    int normalize,
    imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
    iterative_imputation_params_t iterative_params = {},
    int ikernel = 1,
    double dist_normalization_factor = 1.01,
    int n_CVs = 0,
    int n_CV_folds = 10,
    double epsilon = 1e-10,
    unsigned int seed = 0) {

    int n_vertices = X.size();
    int n_features = X[0].size();

    if (ikernel)
        initialize_kernel(ikernel);

    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    std::vector<double> kernel_weights(max_neighbors);
    std::vector<double> vertex_edge_lengths(max_neighbors);

    if (seed == 0) {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    bool mean_adjust = (normalize == 2);
    bool range_adjust = (normalize == 1);

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    auto kernel_diffusion_loop = [&](std::vector<double>& y_current, std::vector<double>& y_next) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            int n_vertex_neighbors = n_neighbors[vertex];
            if (n_vertex_neighbors == 0) continue;

            double max_dist = 0.0;
            for (int j = 0; j < n_vertex_neighbors; ++j) {
                vertex_edge_lengths[j] = edge_lengths[vertex][j];
                if (vertex_edge_lengths[j] > max_dist)
                    max_dist = vertex_edge_lengths[j];
            }

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

    auto X_traj = std::vector<std::vector<std::vector<double>>>(n_time_steps + 1, std::vector<std::vector<double>>(n_vertices, std::vector<double>(n_features)));
    X_traj[0] = X;

    std::vector<double> y_current(n_vertices);
    std::vector<double> y_next(n_vertices);

    for (int feature = 0; feature < n_features; ++feature) {
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            y_current[vertex] = X[vertex][feature];
        }

        double ymin = *std::min_element(y_current.begin(), y_current.end());
        double ymax = *std::max_element(y_current.begin(), y_current.end());
        double Ey = mean_adjust ? std::accumulate(y_current.begin(), y_current.end(), 0.0) / n_vertices : 0.0;

        for (int time_step = 1; time_step <= n_time_steps; ++time_step) {
            kernel_diffusion_loop(y_current, y_next);

            if (range_adjust) {
                scale_to_range(y_next, ymin, ymax);
            } else if (mean_adjust) {
                double Ey_next = std::accumulate(y_next.begin(), y_next.end(), 0.0) / n_vertices;
                double Delta = Ey - Ey_next;
                for (auto& value : y_next) {
                    value += Delta;
                }
            }

            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                X_traj[time_step][vertex][feature] = y_next[vertex];
            }

            std::swap(y_current, y_next);
        }
    }

    std::vector<double> mean_cv_errors(n_time_steps, 0.0);

    if (n_CVs > 0) {
        std::vector<std::vector<double>> cv_errors(n_CVs, std::vector<double>(n_time_steps, 0.0));

        std::vector<std::set<int>> set_graph(n_vertices);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
        }

        int fold_size = n_vertices / n_CV_folds;

        for (int cv = 0; cv < n_CVs; ++cv) {
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

            int cv_imputation_ikernel = ikernel == 0 ? 1 : ikernel;

            for (int feature = 0; feature < n_features; ++feature) {
                for (int vertex = 0; vertex < n_vertices; ++vertex) {
                    y_current[vertex] = X[vertex][feature];
                }

                auto imputed_y = cv_imputation(test_set, graph, edge_lengths, y_current,
                                               false, imputation_method, iterative_params,
                                               false, 0.0, cv_imputation_ikernel,
                                               dist_normalization_factor);

                y_current = std::move(*imputed_y);

                double ymin = *std::min_element(y_current.begin(), y_current.end());
                double ymax = *std::max_element(y_current.begin(), y_current.end());
                double Ey = mean_adjust ? std::accumulate(y_current.begin(), y_current.end(), 0.0) / n_vertices : 0.0;

                for (int time_step = 0; time_step < n_time_steps; ++time_step) {
                    kernel_diffusion_loop(y_current, y_next);

                    if (range_adjust) {
                        scale_to_range(y_next, ymin, ymax);
                    } else if (mean_adjust) {
                        double Ey_next = std::accumulate(y_next.begin(), y_next.end(), 0.0) / n_vertices;
                        double Delta = Ey - Ey_next;
                        for (auto& value : y_next) {
                            value += Delta;
                        }
                    }

                    double cv_error = std::accumulate(test_set.begin(), test_set.end(), 0.0,
                        [&](double sum, int vertex) { return sum + std::abs(y_next[vertex] - X[vertex][feature]); }
                    ) / test_set.size();

                    cv_errors[cv][time_step] += cv_error;

                    std::swap(y_current, y_next);
                }
            }

            // Normalize CV errors by the number of features
            for (int time_step = 0; time_step < n_time_steps; ++time_step) {
                cv_errors[cv][time_step] /= n_features;
            }
        }

        // Compute mean CV errors across all CVs
        for (int time_step = 0; time_step < n_time_steps; ++time_step) {
            mean_cv_errors[time_step] = std::accumulate(cv_errors.begin(), cv_errors.end(), 0.0,
                [time_step](double sum, const std::vector<double>& errors) { return sum + errors[time_step]; }
            ) / n_CVs;
        }
    }

    auto result = std::make_unique<graph_diffusion_matrix_smoother_result_t>();
    result->X_traj = std::move(X_traj);
    result->mean_cv_error = std::move(mean_cv_errors);

    return result;
}


/**
 * @brief Graph Diffusion Matrix Smoother (C++ implementation)
 *
 * This function performs graph diffusion smoothing on a matrix of features.
 * It applies the diffusion process to each column of the input matrix independently,
 * using the provided graph structure.
 *
 * @param r_X SEXP: A numeric matrix where rows represent vertices and columns represent features.
 * @param r_graph SEXP: A list of integer vectors, each containing indices of neighbors for each vertex.
 * @param r_edge_lengths SEXP: A list of numeric vectors containing edge lengths corresponding to neighbors.
 * @param r_weights SEXP: A numeric vector of weights for each vertex.
 * @param r_n_time_steps SEXP: An integer specifying the number of diffusion steps to perform.
 * @param r_step_factor SEXP: A numeric value controlling the magnitude of each diffusion step.
 * @param r_normalize SEXP: An integer specifying the normalization method (0, 1, or 2).
 * @param r_imputation_method SEXP: An integer specifying the imputation method for cross-validation.
 * @param r_iterative_params SEXP: A list of parameters for iterative imputation.
 * @param r_ikernel SEXP: An integer specifying the index of the kernel function to use.
 * @param r_dist_normalization_factor SEXP: A numeric value for normalizing distances in the kernel function.
 * @param r_n_CVs SEXP: An integer specifying the number of cross-validation runs to perform.
 * @param r_n_CV_folds SEXP: An integer specifying the number of folds in each cross-validation run.
 * @param r_epsilon SEXP: A small numeric value to prevent numerical instability.
 * @param r_seed SEXP: An integer seed for the random number generator used in cross-validation.
 *
 * @return SEXP: A list containing two elements:
 *         - X_traj: A list of matrices representing the trajectory of diffusion smoothing versions of X.
 *         - mean_cv_error: A numeric vector of mean cross-validation errors, one per time step.
 *
 * @note This function is intended to be called from R using the .Call interface.
 *       It should not be called directly from C++ code.
 */
SEXP S_graph_diffusion_matrix_smoother(SEXP r_X,
                                       SEXP r_graph,
                                       SEXP r_edge_lengths,
                                       SEXP r_weights,
                                       SEXP r_n_time_steps,
                                       SEXP r_step_factor,
                                       SEXP r_normalize,
                                       SEXP r_imputation_method,
                                       SEXP r_max_iterations,
                                       SEXP r_convergence_threshold,
                                       SEXP r_ikernel,
                                       SEXP r_dist_normalization_factor,
                                       SEXP r_n_CVs,
                                       SEXP r_n_CV_folds,
                                       SEXP r_epsilon,
                                       SEXP r_seed) {
    // Convert R matrix to C++ vector of vectors
    SEXP r_X_dim = Rf_getAttrib(r_X, R_DimSymbol);
    int n_rows = INTEGER(r_X_dim)[0];
    int n_cols = INTEGER(r_X_dim)[1];
    std::vector<std::vector<double>> X(n_rows, std::vector<double>(n_cols));
    double* r_X_ptr = REAL(r_X);
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            X[i][j] = r_X_ptr[i + j * n_rows];
        }
    }

    // Convert R list to C++ vector of vectors for graph
    std::vector<std::vector<int>> graph(Rf_length(r_graph));
    for (int i = 0; i < Rf_length(r_graph); ++i) {
        SEXP r_neighbors = VECTOR_ELT(r_graph, i);
        int* neighbors_ptr = INTEGER(r_neighbors);
        graph[i].assign(neighbors_ptr, neighbors_ptr + Rf_length(r_neighbors));
        // R uses 1-based indexing, convert to 0-based
        for (auto& neighbor : graph[i]) --neighbor;
    }

    // Convert R list to C++ vector of vectors for edge_lengths
    std::vector<std::vector<double>> edge_lengths(Rf_length(r_edge_lengths));
    for (int i = 0; i < Rf_length(r_edge_lengths); ++i) {
        SEXP r_lengths = VECTOR_ELT(r_edge_lengths, i);
        double* lengths_ptr = REAL(r_lengths);
        edge_lengths[i].assign(lengths_ptr, lengths_ptr + Rf_length(r_lengths));
    }

    // Convert R vector to C++ vector for weights
    std::vector<double> weights(REAL(r_weights), REAL(r_weights) + Rf_length(r_weights));

    // Convert other R parameters to C++ types
    int n_time_steps = Rf_asInteger(r_n_time_steps);
    double step_factor = Rf_asReal(r_step_factor);
    int normalize = Rf_asInteger(r_normalize);
    imputation_method_t imputation_method = static_cast<imputation_method_t>(Rf_asInteger(r_imputation_method));
    int max_iterations = INTEGER(r_max_iterations)[0];
    double convergence_threshold = REAL(r_convergence_threshold)[0];

    int ikernel = Rf_asInteger(r_ikernel);
    double dist_normalization_factor = Rf_asReal(r_dist_normalization_factor);
    int n_CVs = Rf_asInteger(r_n_CVs);
    int n_CV_folds = Rf_asInteger(r_n_CV_folds);
    double epsilon = Rf_asReal(r_epsilon);
    unsigned int seed = static_cast<unsigned int>(Rf_asInteger(r_seed));

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    auto result = graph_diffusion_matrix_smoother(X,
                                                  graph,
                                                  edge_lengths,
                                                  weights,
                                                  n_time_steps,
                                                  step_factor,
                                                  normalize,
                                                  imputation_method,
                                                  iterative_params,
                                                  ikernel,
                                                  dist_normalization_factor,
                                                  n_CVs,
                                                  n_CV_folds,
                                                  epsilon,
                                                  seed);
    // Convert the result back to R objects
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));

    // Convert X_traj to R list of matrices
    SEXP r_X_traj = PROTECT(Rf_allocVector(VECSXP, result->X_traj.size()));
    for (size_t t = 0; t < result->X_traj.size(); ++t) {
        SEXP r_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
        double* matrix_ptr = REAL(r_matrix);
        for (int i = 0; i < n_rows; ++i) {
            for (int j = 0; j < n_cols; ++j) {
                matrix_ptr[i + j * n_rows] = result->X_traj[t][i][j];
            }
        }
        SET_VECTOR_ELT(r_X_traj, t, r_matrix);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(r_result, 0, r_X_traj);

    // Convert mean_cv_error to R vector
    SEXP r_mean_cv_error = PROTECT(Rf_allocVector(REALSXP, result->mean_cv_error.size()));
    std::copy(result->mean_cv_error.begin(), result->mean_cv_error.end(), REAL(r_mean_cv_error));
    SET_VECTOR_ELT(r_result, 1, r_mean_cv_error);

    // Set names for the result list
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_names, 0, Rf_mkChar("X_traj"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("mean_cv_error"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(4);
    return r_result;
}
