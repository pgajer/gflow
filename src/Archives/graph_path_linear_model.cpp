#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>

#include "msr2.h"
#include "graph_diffusion_smoother.h"
#include "SEXP_cpp_conversion_utils.h"
#include "kernels.h"

std::unique_ptr<std::vector<double>> cv_imputation(const std::set<int>& test_set,
                                                   const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& edge_lengths,
                                                   const std::vector<double>& y,
                                                   bool y_binary,
                                                   imputation_method_t imputation_method,
                                                   iterative_imputation_params_t iterative_params,
                                                   bool apply_binary_threshold,
                                                   double binary_threshold,
                                                   int ikernel,
                                                   double dist_normalization_factor);

extern "C" {
    SEXP S_graph_kernel_local_linear_1D_model(SEXP s_X,
                                              SEXP s_graph,
                                              SEXP s_edge_lengths,
                                              SEXP s_y,
                                              SEXP s_ikernel,
                                              SEXP s_dist_normalization_factor);

    SEXP S_graph_kernel_local_linear_1D_model_cv(SEXP s_X,
                                                 SEXP s_graph,
                                                 SEXP s_edge_lengths,
                                                 SEXP s_y,
                                                 SEXP s_imputation_method,
                                                 SEXP s_max_iterations,
                                                 SEXP s_convergence_threshold,
                                                 SEXP s_apply_binary_threshold,
                                                 SEXP s_binary_threshold,
                                                 SEXP s_ikernel,
                                                 SEXP s_dist_normalization_factor,
                                                 SEXP s_n_CVs,
                                                 SEXP s_n_CV_folds,
                                                 SEXP s_seed);
}


/**
 * @brief Performs weighted least squares regression.
 *
 * This function calculates the coefficients (β₀ and β₁) for a linear regression
 * model y = β₀ + β₁x using weighted least squares method.
 *
 * @param x A vector of double values representing the independent variable.
 * @param y A vector of double values representing the dependent variable.
 * @param w A vector of double values representing the weights for each observation.
 *
 * @return A pair of doubles where:
 *         - first: β₀ (y-intercept)
 *         - second: β₁ (slope)
 *
 * @note The function assumes that x, y, and w are of equal length and non-empty.
 *       No input validation is performed within the function.
 *
 * @example
 *   std::vector<double> x = {1, 2, 3, 4, 5};
 *   std::vector<double> y = {2, 4, 5, 4, 5};
 *   std::vector<double> w = {1, 1, 1, 0.5, 0.5};
 *   auto [intercept, slope] = weighted_least_squares(x, y, w);
 */
std::pair<double, double> weighted_least_squares(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 const std::vector<double>& w) {
    double sum_w = 0, sum_wx = 0, sum_wy = 0, sum_wxx = 0, sum_wxy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum_w += w[i];
        sum_wx += w[i] * x[i];
        sum_wy += w[i] * y[i];
        sum_wxx += w[i] * x[i] * x[i];
        sum_wxy += w[i] * x[i] * y[i];
    }

    double denominator = sum_w * sum_wxx - sum_wx * sum_wx;
    double beta_0 = (sum_wxx * sum_wy - sum_wx * sum_wxy) / denominator;
    double beta_1 = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator;

    return {beta_0, beta_1};
}

/**
 * @brief Computes a graph kernel-weighted 1D local linear model.
 *
 * This function calculates a local linear model for each point in a one-dimensional
 * dataset, using a graph structure to define neighborhoods and kernel-weighted
 * least squares regression.
 *
 * @param X A vector of doubles representing points on a real line.
 * @param graph A vector of vectors representing the graph structure. Each inner
 *              vector contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors of edge lengths. edge_lengths[i][j] is
 *                     the length of the edge between the i-th vertex and its j-th
 *                     neighbor as listed in graph[i].
 * @param y A vector of double values representing the response variable for each point.
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 * @param dist_normalization_factor A scaling factor applied to the maximum distance
 *                                  between a vertex and its neighbors. Default is 1.01.
 *
 * @return A unique pointer to a vector of doubles containing the fitted y-intercept
 *         values (β₀) for each point, representing the local linear model prediction.
 *
 * @note If the input y values are all either 0 or 1, the output values are clamped
 *       to this range, suitable for binary response variables.
 * @note The function assumes that all input vectors are of appropriate and
 *       consistent sizes. No extensive input validation is performed.
 * @note Memory management of the result is handled automatically using std::unique_ptr.
 * @note The function uses a kernel function defined in kernels.h, which is set based on the ikernel parameter.
 *
 * @example
 *   std::vector<double> X = {1, 2, 3, 4, 5};
 *   std::vector<std::vector<int>> graph = {{1}, {0, 2}, {1, 3}, {2, 4}, {3}};
 *   std::vector<std::vector<double>> edge_lengths = {{1}, {1, 1}, {1, 1}, {1, 1}, {1}};
 *   std::vector<double> y = {0, 0.5, 1, 0.5, 1};
 *   int ikernel = 1;  // Epanechnikov kernel
 *   auto result = graph_kernel_local_linear_1D_model(X, graph, edge_lengths, y, ikernel);
 */
std::unique_ptr<std::vector<double>> graph_kernel_local_linear_1D_model(
    const std::vector<double>& X,
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    int ikernel,
    double dist_normalization_factor = 1.01) {

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    auto result = std::make_unique<std::vector<double>>(X.size());

    // Select kernel function
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    for (size_t i = 0; i < X.size(); ++i) {

        int local_nbhd_size = graph[i].size();
        int x_size = local_nbhd_size + 1;
        std::vector<double> local_x, local_y, local_weights;
        local_x.reserve(x_size);
        local_y.reserve(x_size);
        local_weights.reserve(x_size);

        // Add the point itself
        local_x.push_back(0);
        local_y.push_back(y[i]);
        local_weights.push_back(1.0);

        double max_dist = 0;
        for (size_t j = 0; j < local_nbhd_size; ++j) {
            int neighbor = graph[i][j];
            local_x.push_back(X[neighbor] - X[i]);
            local_y.push_back(y[neighbor]);
            max_dist = std::max(max_dist, edge_lengths[i][j]);
        }

        // Normalize distances and compute local_weights
        max_dist *= dist_normalization_factor;
        std::vector<double> normalized_dist(x_size);
        for (size_t j = 0; j < x_size; ++j) {
            normalized_dist[j] = std::abs(local_x[j]) / max_dist;
        }

        kernel_fn(normalized_dist.data(), x_size, local_weights.data());

        // Fit weighted least squares
        auto [beta_0, beta_1] = weighted_least_squares(local_x, local_y, local_weights);

        // Store the result, clamping to [0, 1] if y is binary
        if (y_binary) {
            (*result)[i] = std::clamp(beta_0, 0.0, 1.0);
        } else {
            (*result)[i] = beta_0;
        }
    }

    return result;
}

/**
 * @brief R/C++ interface function for graph kernel local linear model computation.
 *
 * This function serves as an interface between R and C++, allowing the graph kernel
 * local linear model to be computed from R while leveraging C++ performance.
 *
 * @param s_X SEXP representing a list of numeric vectors, each vector representing
 *            a point in one-dimensional space.
 * @param s_graph SEXP representing a list of integer vectors, each vector containing
 *                the indices of neighboring vertices for a given vertex in the graph.
 * @param s_edge_lengths SEXP representing a list of numeric vectors, each vector
 *                       containing the lengths of edges corresponding to the neighbors
 *                       in s_graph.
 * @param s_y SEXP representing a numeric vector of response variables corresponding
 *            to each point in s_X.
 * @param s_ikernel SEXP representing an integer specifying the kernel function to use.
 *                  Valid values are:
 *                  1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 * @param s_dist_normalization_factor SEXP representing a numeric value for the
 *                                    distance normalization factor.
 *
 * @return SEXP representing a numeric vector containing the fitted y-intercept
 *         values (β₀) for each point, representing the local linear model prediction.
 *
 * @note This function assumes that R_list_of_dvectors_to_cpp_vector_of_dvectors and
 *       Rgraph_to_vector are correctly implemented for SEXP to C++ conversions.
 * @note The function uses std::move for efficient transfer of large data structures,
 *       which is safe as it operates on temporary C++ objects, not the original R objects.
 * @note Memory management for C++ objects is handled automatically through RAII.
 * @note The function uses PROTECT/UNPROTECT for proper R garbage collection management.
 *
 * @warning This function does not perform extensive error checking on input parameters.
 *          Ensure that inputs are of the correct type and size before calling.
 *
 * @see graph_kernel_local_linear_model for the underlying C++ implementation.
 */
SEXP S_graph_kernel_local_linear_1D_model(SEXP s_X,
                                          SEXP s_graph,
                                          SEXP s_edge_lengths,
                                          SEXP s_y,
                                          SEXP s_ikernel,
                                          SEXP s_dist_normalization_factor) {

    std::vector<double> X = std::move(*Rvect_to_CppVect_double(s_X));
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);

    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);

    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    std::unique_ptr<std::vector<double>> result = graph_kernel_local_linear_1D_model(X, graph, edge_lengths, y, ikernel, dist_normalization_factor);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);

    return Rresult;
}


/**
 * @brief Computes the cross-validation error of kernel-weighted local linear models of a function over vertices of a graph.
 *
 * This function calculates the cross-validation error of the kernel-weighted
 * local linear models for each vertex in the graph. The weighting is based on
 * edge lengths, transformed by a specified kernel function. For binary
 * outcomes, cross-entropy loss is used, while mean absolute error is used for non-binary outcomes.
 *
 * @param X A vector of doubles representing points on a real line.
 * @param graph A vector of vectors representing the graph structure. Each inner vector contains the indices of neighboring vertices for a given vertex.
 * @param edge_lengths A vector of vectors of edge lengths. edge_lengths[i][j] is the length of the edge between the i-th vertex and its j-th neighbor as listed in graph[i].
 * @param y A vector of double values associated with each vertex in the graph. If y contains only values 0 and 1, it will be treated as a binary outcome.
 * @param imputation_method Specifies the method for imputing binary outcomes. Options are:
 *        - LOCAL_MEAN_THRESHOLD: Uses the mean of y computed over the training vertices (default).
 *        - NEIGHBORHOOD_MATCHING: Uses a matching method based on local neighborhood statistics.
 *        - ITERATIVE_NEIGHBORHOOD_MATCHING: Uses an iterative version of the NEIGHBORHOOD_MATCHING method.
 *        - SUPPLIED_THRESHOLD: Uses a user-supplied threshold value.
 *        - GLOBAL_MEAN_THRESHOLD: Uses the global mean of y across all vertices.
 * @param iterative_params Parameters for iterative imputation methods (only used with certain imputation methods).
 * @param apply_binary_threshold Whether to apply a threshold for binary classification.
 * @param binary_threshold The threshold value to use when imputation_method is SUPPLIED_THRESHOLD. Default is 0.5.
 * @param ikernel An integer specifying the kernel function to use. Valid values are:
 *                1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal). Default is 1.
 * @param dist_normalization_factor A scaling factor applied to the maximum distance between a vertex and its neighbors. Default value is 1.01.
 * @param n_CVs The number of cross-validation rounds. Default is 0.
 * @param n_CV_folds The number of folds in each cross-validation round. Default is 10.
 * @param epsilon A small positive constant for clipping the estimated conditional expectation values of a binary variable. Default is 1e-10.
 * @param seed A seed for the random number generator to ensure reproducibility. Default is 0.
 *
 * @return A unique pointer to a vector of cross-validation errors of the kernel-weighted local linear model for each vertex.
 *
 * @note For binary outcomes, cross-entropy loss is used as the error metric. For non-binary outcomes, mean absolute error is used.
 * @note If n_CVs is 0, no cross-validation is performed, and the returned vector will contain zeros.
 * @note The function uses C++11 random number generation facilities for reproducibility when a seed is provided.
 * @note The choice of imputation_method affects how binary outcomes are classified:
 *       - LOCAL_MEAN_THRESHOLD adapts to local data distribution in each cross-validation fold.
 *       - NEIGHBORHOOD_MATCHING attempts to preserve local statistical properties of the graph.
 *       - ITERATIVE_NEIGHBORHOOD_MATCHING an iterative version of the NEIGHBORHOOD_MATCHING method.
 *       - SUPPLIED_THRESHOLD allows for domain-specific threshold selection.
 *       - GLOBAL_MEAN_THRESHOLD uses a single threshold based on the entire dataset.
 *
 * @example
 * std::vector<double> X = {1, 2, 3, 4, 5};
 * std::vector<std::vector<int>> graph = {{1}, {0, 2}, {1, 3}, {2, 4}, {3}};
 * std::vector<std::vector<double>> edge_lengths = {{1}, {1, 1}, {1, 1}, {1, 1}, {1}};
 * std::vector<double> y = {0, 0.5, 1, 0.5, 1};
 * auto result = graph_kernel_local_linear_1D_model_cv(X, graph, edge_lengths, y, imputation_method_t::LOCAL_MEAN_THRESHOLD, {}, true, 0.5, 1, 1.01, 10, 5, 1e-10, 42);
 */
std::unique_ptr<std::vector<double>> graph_kernel_local_linear_1D_model_cv(
    const std::vector<double>& X,
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
    iterative_imputation_params_t iterative_params = {},
    bool apply_binary_threshold = true,
    double binary_threshold = 0.5, // Only used if imputation_method is SUPPLIED_THRESHOLD
    int ikernel = 1,
    double dist_normalization_factor = 1.01,
    int n_CVs = 0,
    int n_CV_folds = 10,
    unsigned int seed = 0) {

    int n_vertices = y.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && imputation_method == imputation_method_t::GLOBAL_MEAN_THRESHOLD) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    auto cv_error = std::make_unique<std::vector<double>>(n_vertices, 0.0); // Mean Absolute Deviation error - for each vertex the mean over all folds and cross-validation iterations
    std::vector<int> cv_error_count(n_vertices, 0);

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

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
                    // int vertex = rand() % n_vertices;
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

            std::vector<double> cv_y = std::move(*cv_imputation(test_set,
                                                                graph,
                                                                edge_lengths,
                                                                y,
                                                                y_binary,
                                                                imputation_method,
                                                                iterative_params,
                                                                apply_binary_threshold,
                                                                binary_threshold,
                                                                ikernel,
                                                                dist_normalization_factor));

            // Estimating the conditional expectation of cv_y
            std::vector<double> Ecv_y = std::move(*graph_kernel_local_linear_1D_model(X, graph, edge_lengths, cv_y, ikernel, dist_normalization_factor));

            #if 0
            // Computing error over test vertices
            if (y_binary) {
                // For binary outcomes, compute cross-entropy loss
                // This is more appropriate than mean absolute error for binary classification
                for (const auto& vertex : test_set) {
                    // Clip Ecv_y to avoid log(0) or log(1) issues
                    double clipped_Ecv_y = std::max(epsilon, std::min(1.0 - epsilon, Ecv_y[vertex]));

                    // Compute negative log-likelihood (cross-entropy loss)
                    (*cv_error)[vertex] += -(y[vertex] * log(clipped_Ecv_y) + (1 - y[vertex]) * log(1 - clipped_Ecv_y));
                    cv_error_count[vertex]++;
                }
            } else {
            #endif


            for (const auto& vertex : test_set) {
                (*cv_error)[vertex] += std::abs(Ecv_y[vertex] - y[vertex]);
                cv_error_count[vertex]++;
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex])
            (*cv_error)[vertex] /= cv_error_count[vertex];
    }

    return cv_error;
}


/**
 * @brief R interface for computing cross-validation error of kernel-weighted local linear models on a graph.
 *
 * This function serves as an interface between R and C++, allowing the cross-validation
 * error of kernel-weighted local linear models to be computed from R while leveraging
 * C++ performance. It converts R objects to C++ types, calls the C++ implementation,
 * and returns the results back to R.
 *
 * @param s_X SEXP representing a numeric vector of points on a real line.
 * @param s_graph SEXP representing a list of integer vectors, each containing the indices
 *                of neighboring vertices for a given vertex in the graph (1-based indexing).
 * @param s_edge_lengths SEXP representing a list of numeric vectors, each containing the
 *                       lengths of edges corresponding to the neighbors in s_graph.
 * @param s_y SEXP representing a numeric vector of response variables for each vertex.
 * @param s_imputation_method SEXP representing a character or integer specifying the imputation method. Options are:
 *        - local_mean_threshold: Uses the mean of y computed over the training vertices (default).
 *        - neighborhood_matching: Uses a matching method based on local neighborhood statistics.
 *        - iterative_neighborhood_matching: Uses an iterative version of the NEIGHBORHOOD_MATCHING method.
 *        - supplied_threshold: Uses a user-supplied threshold value.
 *        - global_mean_threshold: Uses the global mean of y across all vertices.
 * @param s_max_iterations The number of iterations in the iterative matching method.
 * @param s_convergence_threshold The convergence threshold in the iterative matching method.
 * @param s_apply_binary_threshold SEXP representing a logical value indicating whether
 *                                 to apply binary threshold.
 * @param s_binary_threshold SEXP representing a numeric value for the binary threshold.
 * @param s_ikernel SEXP representing an integer specifying the kernel function to use.
 * @param s_dist_normalization_factor SEXP representing a numeric value for distance normalization.
 * @param s_n_CVs SEXP representing an integer for the number of cross-validation rounds.
 * @param s_n_CV_folds SEXP representing an integer for the number of folds in each CV round.
 * @param s_epsilon SEXP representing a numeric value for the small constant in error calculation.
 * @param s_seed SEXP representing an integer seed for random number generation.
 *
 * @return SEXP representing a numeric vector of cross-validation errors for each vertex.
 *
 * @note This function assumes that all input SEXPs are of the correct type and structure.
 *       It's recommended to perform type checking and error handling in the R wrapper.
 * @note The graph indices in s_graph are assumed to be 1-based (R-style) and are
 *       converted to 0-based indices for the C++ function.
 * @note Memory management for R objects is handled using PROTECT/UNPROTECT.
 *
 * @see graph_kernel_local_linear_1D_model_cv for the underlying C++ implementation.
 */
SEXP S_graph_kernel_local_linear_1D_model_cv(SEXP s_X,
                                             SEXP s_graph,
                                             SEXP s_edge_lengths,
                                             SEXP s_y,
                                             SEXP s_imputation_method,
                                             SEXP s_max_iterations,
                                             SEXP s_convergence_threshold,
                                             SEXP s_apply_binary_threshold,
                                             SEXP s_binary_threshold,
                                             SEXP s_ikernel,
                                             SEXP s_dist_normalization_factor,
                                             SEXP s_n_CVs,
                                             SEXP s_n_CV_folds,
                                             SEXP s_seed) {
    // Convert X to std::vector<double>
    std::vector<double> X(REAL(s_X), REAL(s_X) + LENGTH(s_X));

    // Convert graph to std::vector<std::vector<int>>
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);

    // Convert edge_lengths to std::vector<std::vector<double>>
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);

    // Convert y to std::vector<double>
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));

    // Convert imputation_method
    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(s_imputation_method)[0]);

    // Convert iterative_params
    int max_iterations = INTEGER(s_max_iterations)[0];
    double convergence_threshold = REAL(s_convergence_threshold)[0];
    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Convert other parameters
    bool apply_binary_threshold = LOGICAL(s_apply_binary_threshold)[0];
    double binary_threshold = REAL(s_binary_threshold)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = static_cast<unsigned int>(INTEGER(s_seed)[0]);

    // Call the C++ function
    std::unique_ptr<std::vector<double>> result = graph_kernel_local_linear_1D_model_cv(X,
                                                                                        graph,
                                                                                        edge_lengths,
                                                                                        y,
                                                                                        imputation_method,
                                                                                        iterative_params,
                                                                                        apply_binary_threshold,
                                                                                        binary_threshold,
                                                                                        ikernel,
                                                                                        dist_normalization_factor,
                                                                                        n_CVs,
                                                                                        n_CV_folds,
                                                                                        seed);
    // Convert the result back to an R numeric vector
    SEXP r_result = PROTECT(allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(r_result)[i] = (*result)[i];
    }

    UNPROTECT(1);
    return r_result;
}


/**
 * @brief Optimizes graph path linear model estimation by finding the optimal neighborhood size (h)
 *        and computing credible intervals.
 *
 * @details This function performs several steps:
 *          1. For each h in [h_min, h_max]:
 *             - Constructs h-hop neighborhood (hHN) graphs
 *             - Computes cross-validation errors
 *          2. Finds optimal h that minimizes CV error
 *          3. Computes conditional expectations using optimal h
 *          4. Optional: Computes Bayesian bootstrap credible intervals
 *
 * @param graph Original adjacency matrix [n_vertices x n_vertices]
 * @param edge_lengths Edge weights/distances matrix [n_vertices x n_vertices]
 * @param y Observed values at each vertex [n_vertices]
 * @param y_true True values for error calculation [n_vertices]
 * @param h_min Minimum neighborhood size to consider (default: 2)
 * @param h_max Maximum neighborhood size to consider (default: 30)
 * @param n_CVs Number of cross-validation iterations (default: 1000)
 * @param n_CV_folds Number of folds for cross-validation (default: 10)
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 500, 0 to skip)
 * @param use_median Use median instead of mean for central tendency (default: false)
 * @param ikernel Kernel function selector (default: 1)
 * @param n_cores Number of cores for parallel computation (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param seed Random seed for reproducibility (default: 0)
 *
 * @return adaptive_nbhd_size_t struct containing:
 *         - graphs: Vector of h-hop neighborhood graphs for each h
 *         - cv_errors: Cross-validation errors for each h
 *         - opt_h: Optimal h value that minimizes CV error
 *         - opt_h_graph: Graph with optimal h
 *         - condEy: Conditional expectations using optimal h
 *         - bb_condEy: Bootstrap central tendency (if n_bb > 0)
 *         - cri_L: Lower credible interval bounds (if n_bb > 0)
 *         - cri_U: Upper credible interval bounds (if n_bb > 0)
 *
 * @throws std::invalid_argument If:
 *         - Input matrices have inconsistent dimensions
 *         - h_min > h_max
 *         - Invalid parameter values (n_CVs, n_CV_folds, p, etc.)
 */
#if 0
adaptive_nbhd_size_t
    adaptive_nbhd_size_graph_path_linear_model(const std::vector<std::vector<int>>& graph,
                                               const std::vector<std::vector<double>>& edge_lengths,
                                               const std::vector<double>& y,
                                               const std::vector<double>& y_true,
                                               bool use_median = false,
                                               int h_min = 2,
                                               int h_max = 30,
                                               int n_CVs = 1000,
                                               int n_CV_folds = 10,
                                               double p = 0.95,
                                               int n_bb = 500,
                                               int ikernel = 1,
                                               int n_cores = 1,
                                               double dist_normalization_factor = 1.01,
                                               double epsilon = 1e-15,
                                               unsigned int seed = 0) {

    #define DEBUG__adaptive_nbhd_size_graph_path_linear_model 0
    #if DEBUG__adaptive_nbhd_size_graph_path_linear_model
    Rprintf("In adaptive_nbhd_size_graph_path_linear_model()\n");
    print_vect_vect(graph,"graph");
    print_vect_vect(edge_lengths,"edge_lengths");
    #endif

    int n_vertices = static_cast<int>(y.size());

    adaptive_nbhd_size_t results;
    results.graphs.resize(h_max - h_min + 1);
    results.cv_errors.resize(h_max - h_min + 1);
    results.h_values.resize(h_max - h_min + 1);

    for (int i = 0, h = h_min; h <= h_max; h++, i++) {
        results.h_values[i] = h;

        // Creating an hHN graph for the given value of h
        auto hhn_graph = create_hHN_graph(graph, edge_lengths, h);

        #if DEBUG__adaptive_nbhd_size_graph_path_linear_model
        Rprintf("h: %d\n", h);
        Rprintf("hhn_graph.first.size(): %d\n", (int)hhn_graph.first.size());
        Rprintf("hhn_graph.second.size(): %d\n", (int)hhn_graph.second.size());

        print_vect_vect(hhn_graph.first,"hhn_graph.first");
        print_vect_vect(hhn_graph.second,"hhn_graph.second");
        Rf_error("Debugging");
        #endif

        // Computing CV errors
        auto errors = graph_kmean_cv(hhn_graph.first,
                                     hhn_graph.second,
                                     y,
                                     ikernel,
                                     dist_normalization_factor,
                                     n_CVs,
                                     n_CV_folds,
                                     epsilon,
                                     seed);

        // Calculate mean error across vertices
        double total_error = std::accumulate(errors.begin(), errors.end(), 0.0);
        results.cv_errors[i] = total_error / n_vertices;
        results.graphs[i] = std::move(hhn_graph);
    }

    // Find the optimal h (minimum CV error)
    auto min_it = std::min_element(results.cv_errors.begin(), results.cv_errors.end());
    int opt_h_idx = std::distance(results.cv_errors.begin(), min_it);
    results.opt_h = h_min + opt_h_idx;

    // Store optimal graph
    results.opt_h_graph = results.graphs[opt_h_idx];

    // Compute conditional expectations using optimal graph
    results.condEy = graph_kmean(results.opt_h_graph.first,
                                 results.opt_h_graph.second,
                                 y,
                                 ikernel,
                                 dist_normalization_factor);

    // Compute true errors
    if (!y_true.empty() && y_true.size() == n_vertices) {
        results.true_errors.resize(n_vertices);
        for (size_t i = 0; i < n_vertices; i++) {
            results.true_errors[i] = std::abs(y_true[i] - results.condEy[i]);
        }
    } else {
        results.true_errors.clear();  // Ensure empty if no true values
    }

    // Optional: Compute bootstrap credible intervals
    if (n_bb > 0) {
        bb_cri_t bb_cri_results = graph_kmean_bb_cri(results.opt_h_graph.first,
                                                     results.opt_h_graph.second,
                                                     y,
                                                     p,
                                                     n_bb,
                                                     use_median,
                                                     n_cores,
                                                     ikernel,
                                                     dist_normalization_factor,
                                                     epsilon);

        results.bb_condEy = std::move(bb_cri_results.bb_Ey);
        results.cri_L = std::move(bb_cri_results.cri_L);
        results.cri_U = std::move(bb_cri_results.cri_U);
    }

    #if DEBUG__adaptive_nbhd_size_graph_path_linear_model
    // Before returning results, verify the data is stored correctly
    Rprintf("Verification before return:\n");
    for (size_t i = 0; i < results.graphs.size(); ++i) {
        Rprintf("h = %d:\n", h_min + static_cast<int>(i));
        Rprintf("  Graph size: %zu\n", results.graphs[i].first.size());
        Rprintf("  CV error: %f\n", results.cv_errors[i]);
    }
    #endif

    return results;
}
#endif
