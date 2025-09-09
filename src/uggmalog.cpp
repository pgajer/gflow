#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

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

#include "ulogit.hpp"
#include "graph_utils.hpp"
#include "uniform_grid_graph.hpp"
#include "centered_paths.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_shortest_path.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "kernels.h"

extern "C" {
    SEXP S_uggmalog(
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
        SEXP s_max_iterations,
        SEXP s_ridge_lambda,
        SEXP s_tolerance,
        SEXP s_verbose
        );
}

/**
 * @brief Implements the Uniform Grid Graph Model-Averaged LOGistic regression (UGGMALOG) algorithm
 *
 * @details This function implements a sophisticated algorithm for analyzing weighted graphs using
 * local path logistic models with model averaging. The algorithm performs the following main steps:
 * 1. Computes graph diameter and determines bandwidth range
 * 2. Creates a uniform grid representation of the input graph
 * 3. For each candidate bandwidth:
 *    - Processes paths through grid vertices
 *    - Fits local logistic models to path data
 *    - Computes weighted predictions and errors
 * 4. Determines optimal bandwidth based on cross-validation errors
 *
 * The algorithm uses weighted logistic regression on paths through the graph to create
 * local models, which are then combined using weighted averaging. Model evaluation
 * is performed using leave-one-out cross-validation with Brier score errors.
 *
 * @details Algorithm Performance Characteristics:
 * - Time Complexity: O(V * P * M) where:
 *   V = number of grid vertices
 *   P = average number of paths per vertex
 *   M = maximum iterations in logistic regression
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
 * @param fit_quadratic Whether to include quadratic term in logistic regression (default: false)
 * @param max_iterations Maximum number of iterations for logistic regression fitting (default: 100)
 * @param ridge_lambda Ridge penalty parameter for logistic regression (default: 0.0)
 * @param tolerance Convergence tolerance for logistic regression (default: 1e-8)
 * @param with_errors Whether to compute LOOCV errors (default: true)
 * @param verbose Enable progress reporting and warnings (default: false)
 *
 * @return uggmalog_t structure containing:
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
 * @warning Logistic regression may not converge for some paths
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
struct uggmalog_t {
    std::vector<double> candidate_bws;
    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> mean_errors;
    double opt_bw;
    size_t opt_bw_index;
    std::vector<double> opt_predictions;
    double graph_diameter;
};

uggmalog_t uggmalog(
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
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool verbose
    ) {

    const bool with_errors = true;

    // Validation of input parameters will be done in the parent R function
    uggmalog_t result;

    // Step 1: Compute graph diameter and bandwidth range
    result.graph_diameter = get_vertex_eccentricity(adj_list, weight_list, start_vertex);

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

    // Step 2: Create uniform grid graph
    uniform_grid_graph_t uniform_grid_graph = create_uniform_grid_graph(
        adj_list,
        weight_list,
        grid_size,
        start_vertex,
        snap_tolerance
    );

    edge_weights_t edge_weights = precompute_edge_weights(adj_list, weight_list);

    // Initialize storage for predictions and errors
    const size_t n_vertices = uniform_grid_graph.n_original_vertices;
    result.bw_predictions.resize(n_bws);
    for (auto& predictions : result.bw_predictions) {
        predictions.resize(n_vertices);  // This actually creates the elements
    }
    result.mean_errors.resize(n_bws);

    // Creating an extended eigen_ulogit_t structure that contains mean_error component. The models will be sorted w/r to mean_error
    struct ext_eigen_ulogit_t : public eigen_ulogit_t {
        double mean_error;
        std::vector<size_t> vertices;
        std::vector<double> w_path;

        // Add a constructor that takes the base class
        explicit ext_eigen_ulogit_t(const eigen_ulogit_t& base)
            : eigen_ulogit_t(base),  // Initialize base class members
              mean_error(0.0),       // Initialize new members
              vertices(),            // Empty vector
              w_path()              // Empty vector
            {}

        // This operator enables automatic sorting in set/multiset
        bool operator<(const ext_eigen_ulogit_t& other) const {
            return mean_error < other.mean_error;
        }
    };

    // weight/prediction/error struct
    struct wpe_t {
        double weight;
        double prediction;
        double error;

        // Constructor needed for emplace_back(x,y,z)
        wpe_t(double w, double p, double e)
            : weight(w), prediction(p), error(e) {}
    };

    // Precomputing paths for each grid vertex
    std::unordered_map<int, reachability_map_t> vertex_paths;
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
        fit_quadratic,
        max_iterations,
        ridge_lambda,
        tolerance,
        with_errors,
        // Debug parameter
        &verbose
        ](
            double bandwidth,
            const uniform_grid_graph_t& uniform_grid_graph,
            const std::vector<double>& y,
            std::vector<double>& predictions,
            double& mean_error
            ) {

        // This lambda processes a single bandwidth value to:
        // 1. Fit local logistic models for paths through each grid vertex
        // 2. Select best models based on error and coverage criteria
        // 3. Compute weighted average predictions and errors
        // Returns: Updates predictions vector and mean_error for this bandwidth

        // Get the number of vertices in the original graph
        const size_t n_vertices = uniform_grid_graph.n_original_vertices;
        std::vector<std::vector<wpe_t>> wpe(n_vertices); // wpe[i] stores a vector of {weight, prediction, error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions

        // Phase 1: fit logistic regression models to the path data associated with each grid vertex
        for (const auto& grid_vertex : uniform_grid_graph.grid_vertices) {
            // Each grid vertex serves as a reference point for finding paths
            // and fitting local logistic models

            // Collect paths that pass through this grid vertex
            // std::vector<path_data_t> paths = ugg_get_path_data(
            //     uniform_grid_graph,
            //     y,
            //     grid_vertex,
            //     bandwidth,
            //     dist_normalization_factor,
            //     min_path_size,
            //     diff_threshold,
            //     kernel_type,
            //     verbose
            //     );
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

            // Fitting logistic models to valid paths through the current grid vertex
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
            std::multiset<ext_eigen_ulogit_t> all_models;
            for (const auto& path : paths) {
                // Fit weighted logistic model to the current path
                eigen_ulogit_t fit_result = eigen_ulogit_fit(
                    path.x_path.data(),
                    path.y_path.data(),
                    path.w_path,
                    fit_quadratic,
                    max_iterations,
                    ridge_lambda,
                    tolerance,
                    with_errors
                    );

                if (!fit_result.converged) {
                    if (verbose) {
                        REPORT_WARNING("Warning: Model fitting did not converge for path through vertex %d at bandwidth %f\n",
                                       grid_vertex, bandwidth);
                    }
                    continue;
                } else {
                    // Create the extended object from fit_result
                    ext_eigen_ulogit_t extended_result(fit_result);

                    // Calculate the mean of loocv_brier_errors vector
                    // We'll use std::accumulate to sum all elements and divide by the vector size
                    double sum = std::accumulate(fit_result.loocv_brier_errors.begin(),
                                                 fit_result.loocv_brier_errors.end(),
                                                 0.0);  // Use 0.0 to ensure floating-point arithmetic

                    extended_result.mean_error = sum / fit_result.loocv_brier_errors.size();
                    extended_result.vertices   = std::move(path.vertices);
                    extended_result.w_path     = std::move(path.w_path);

                    all_models.insert(extended_result);
                }
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
                    wpe[ model.vertices[i] ].emplace_back(model.w_path[i], model.predictions[i], model.loocv_brier_errors[i]);
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

    // Step 3: Process each bandwidth
    for (size_t bw_idx = 0; bw_idx < n_bws; ++bw_idx) {

        if (verbose && (bw_idx % 5 == 0)) {
            Rprintf("Processing bandwidth %zu of %zu\n", (size_t)bw_idx + 1, (size_t)n_bws);
        }

        process_bw(
            result.candidate_bws[bw_idx],
            uniform_grid_graph,
            y,
            result.bw_predictions[bw_idx],
            result.mean_errors[bw_idx]
            );

        if (result.mean_errors[bw_idx] == std::numeric_limits<double>::infinity()) {
            if (verbose) {
                REPORT_WARNING("Warning: No valid result.mean_errors[bw_idx] predictions for bandwidth index %d and bandwidth %f\n",
                               (int)bw_idx, result.candidate_bws[bw_idx]);
            }
        }
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

    return result;
}

/**
 * @brief R interface for the Uniform Grid Graph Model-Averaged LOGistic regression (UGGMALOG) algorithm
 *
 * @details This function serves as the interface between R and the C++ implementation of the UGGMALOG
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
 * @param s_max_iterations SEXP (integer) Maximum fitting iterations
 * @param s_ridge_lambda SEXP (numeric) Ridge penalty parameter
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
 *     * max_iterations: Integer
 *     * ridge_lambda: Numeric
 *     * tolerance: Numeric
 *
 * @note All R indices are 1-based and converted to 0-based for C++
 * @note The function handles memory protection for all SEXP objects
 * @note Returns errors to R using the error() mechanism
 */
SEXP S_uggmalog(
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
    SEXP s_max_iterations,
    SEXP s_ridge_lambda,
    SEXP s_tolerance,
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
    size_t max_iterations = INTEGER(s_max_iterations)[0];
    double ridge_lambda = REAL(s_ridge_lambda)[0];
    double tolerance = REAL(s_tolerance)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    uggmalog_t result = uggmalog(
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
        max_iterations,
        ridge_lambda,
        tolerance,
        verbose
        );

    // Convert results to R list
    size_t n_protected = 0;
    const size_t RESULT_LIST_SIZE = 6;
    SEXP r_result = PROTECT(allocVector(VECSXP, RESULT_LIST_SIZE)); n_protected++;

    // Updated names including x_grid
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, RESULT_LIST_SIZE)); n_protected++;
    SET_STRING_ELT(r_names, 0, Rf_mkChar("candidate_bws"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("bw_predictions"));
    SET_STRING_ELT(r_names, 2, Rf_mkChar("mean_errors"));
    SET_STRING_ELT(r_names, 3, Rf_mkChar("opt_bw_idx"));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("predictions")); // opt_predictions
    SET_STRING_ELT(r_names, 5, Rf_mkChar("graph_diameter"));

    // Candidate bandwidths
    SEXP r_candidate_bws = PROTECT(Rf_allocVector(REALSXP, result.candidate_bws.size())); n_protected++;
    std::copy(result.candidate_bws.begin(), result.candidate_bws.end(),
              REAL(r_candidate_bws));
    SET_VECTOR_ELT(r_result, 0, r_candidate_bws);

    // Bandwidth predictions
    if (n_bws != result.bw_predictions.size()) { // result.bw_predictions should have n_bws elements
        REPORT_ERROR("n_bws: %ld is not equal to result.bw_predictions.size(): %ld\n",
            n_bws, result.bw_predictions.size());
    }
    SEXP r_bw_predictions;
    r_bw_predictions = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_bws)); n_protected++; // the i-th column of r_bw_predictions is result.bw_predictions[i]
    for (size_t i = 0; i < n_bws; ++i) {
        if (!result.bw_predictions[i].empty()) {

            if (n_vertices != result.bw_predictions[i].size()) {
                REPORT_ERROR("n_vertices: %ld is not equal to result.bw_predictions.size(): %ld\n",
                             n_vertices, result.bw_predictions.size());
            }

            // Copy column by column for column-major order
            std::copy(result.bw_predictions[i].begin(),
                      result.bw_predictions[i].end(),
                      REAL(r_bw_predictions) + i * n_vertices);
        }
    }
    SET_VECTOR_ELT(r_result, 1, r_bw_predictions);

    // Mean errors
    SEXP r_mean_errors = PROTECT(Rf_allocVector(REALSXP, result.mean_errors.size())); n_protected++;
    std::copy(result.mean_errors.begin(), result.mean_errors.end(), REAL(r_mean_errors));
    SET_VECTOR_ELT(r_result, 2, r_mean_errors);

    // Optimal indices
    SEXP r_opt_bw_idx = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(r_opt_bw_idx)[0] = result.opt_bw_index + 1; // Convert to 1-based indexing
    SET_VECTOR_ELT(r_result, 3, r_opt_bw_idx);

    // Optimal predictions
    SEXP r_opt_predictions = PROTECT(Rf_allocVector(REALSXP, n_vertices)); n_protected++;
    std::copy(result.opt_predictions.begin(), result.opt_predictions.end(), REAL(r_opt_predictions));
    SET_VECTOR_ELT(r_result, 4, r_opt_predictions);

    // Graph diameter
    SEXP r_graph_diameter = PROTECT(ScalarReal(result.graph_diameter)); n_protected++;
    SET_VECTOR_ELT(r_result, 5, r_graph_diameter);

    // Set names for the main list
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(n_protected);

    return r_result;
}
