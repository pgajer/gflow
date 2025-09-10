#ifdef _OPENMP
//#include <omp.h>
#include "omp_compat.h"
#endif

// Undefine R's match macro if it exists
#ifdef match
    #undef match
#endif

#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <execution>
#include <atomic>
#include <mutex>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <unordered_map>

#include "exec_policy.hpp"
#include "sampling.h" // for C_runif_simplex()
#include "error_utils.h"
#include "pglm.h"
#include "msr2.h"
#include "path_graphs.hpp"
#include "cpp_utils.hpp"                 // for elapsed_time
#include "SEXP_cpp_conversion_utils.hpp"
#include "kernels.h"
#include "1D_linear_models.h"
#include "predictive_errors.hpp"

lm_loocv_t predict_lm_1d_loocv(const std::vector<double>& y,
                               const std::vector<double>& x,
                               const std::vector<double>& w,
                               const std::vector<int>& vertex_indices,
                               int ref_index,
                               bool y_binary,
                               double epsilon = 1e-8);

std::vector<lm_loocv_t> predict_lms_1d_loocv(const std::vector<double>& y,
                                             const std::vector<double>& x,
                                             const std::vector<int>& vertex_indices,
                                             const std::vector<std::vector<double>>& w_list,
                                             bool y_binary,
                                             double epsilon = 1e-8);

lm_mae_t predict_lm_1d_with_mae(const std::vector<double>& y,
                                std::vector<double>& x,
                                const std::vector<double>& w,
                                int ref_index,
                                double epsilon = 1e-8);

double predict_lm_1d(const std::vector<double>& y,
                     std::vector<double>& x,
                     const std::vector<double>& w,
                     int ref_index,
                     double epsilon = 1e-8);

std::pair<std::vector<double>, std::vector<int>> pgmalo_with_cv_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon);


/**
 * @brief Performs k-fold cross-validation for Path Graph Linear Models (PGMALO)
 *
 * @details This function evaluates PGMALO prediction performance using k-fold cross-validation.
 * For each iteration, it:
 * 1. Randomly selects test vertices
 * 2. Sets their weights to 0 (effectively excluding them from model fitting)
 * 3. Fits PGMALO using remaining vertices as training data
 * 4. Computes prediction errors for test vertices
 * 5. Averages errors across all iterations where predictions were possible
 *
 * Special handling is implemented when fold_size=1 and n_vertices=n_CVs, where each vertex
 * becomes its own test set exactly once.
 *
 * @param path_graph Pre-computed path graph structure containing topology and path information
 * @param y Response values for each vertex (can be binary or continuous)
 * @param n_CVs Number of cross-validation iterations to perform
 * @param n_CV_folds Number of folds for each cross-validation iteration
 * @param seed Random seed for test set selection (0 uses current time)
 * @param ikernel Kernel function identifier for weight computation
 * @param max_distance_deviation Maximum allowed deviation from optimal position in path
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical tolerance for zero comparisons (default: 1e-8)
 *
 * @return Vector of mean absolute deviation (MAD) errors for each vertex
 *         - For vertices that could be predicted in at least one iteration:
 *           returns average absolute Rf_error across all successful predictions
 *         - For vertices that could never be predicted: returns NaN
 *
 * @note
 * - Uses pgmalo_with_cv_weights() for model fitting with test set exclusion
 * - Vertices may be unpredictable in some iterations if they lack sufficient
 *   neighboring vertices with non-zero weights
 * - Error for each vertex is averaged only over iterations where prediction
 *   was possible
 * - Test sets are randomly selected unless fold_size=1 and n_vertices=n_CVs
 *
 * Example usage:
 * @code
 * path_graph_plm_t path_graph = create_path_graph_plm(adj_list, weights, h);
 * std::vector<double> y = {0.0, 1.0, 0.0, 1.0};  // Response values
 * int n_CVs = 10;        // Number of CV iterations
 * int n_CV_folds = 4;    // 4-fold cross-validation
 * unsigned int seed = 42; // Random seed
 *
 * auto cv_errors = pgmalo_cv(path_graph, y, n_CVs, n_CV_folds, seed,
 *                         1, 2, 1.01, 1e-8);
 * @endcode
 *
 * @pre path_graph must be properly initialized
 * @pre y.size() must equal number of vertices in path_graph
 * @pre n_CVs > 0
 * @pre n_CV_folds > 0 and n_CV_folds <= number of vertices
 *
 * @see path_graph_plm_t
 * @see pgmalo_with_cv_weights
 * @see create_path_graph_plm
 *
 * Time Complexity: O(n_CVs * V * P * N), where:
 * - V is number of vertices
 * - P is average number of paths per vertex
 * - N is average path length
 */
std::vector<double> pgmalo_cv(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_CVs,
    int n_CV_folds,
    unsigned int seed,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    int n_vertices = static_cast<int>(y.size());

    auto cv_error = std::vector<double>(n_vertices, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> cv_error_count(n_vertices, 0);

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    // Creating a set version of the adjacency matrix of the graph
    std::vector<std::set<int>> set_graph(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        set_graph[vertex].insert(path_graph.adj_list[vertex].begin(), path_graph.adj_list[vertex].end());
    }

    int fold_size = n_vertices / n_CV_folds;

    std::vector<double> weights(n_vertices, 1.0);

    // cross-validation loop
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

        // Reseting all weights to 1 and set weights over test_set to 0
        std::fill(weights.begin(), weights.end(), 1.0);
        for (const auto& vertex : test_set) {
            weights[vertex] = 0.0;
        }

        // Estimating the conditional expectation of cv_y
        auto res = pgmalo_with_cv_weights(path_graph,
                                        y,
                                        weights,
                                        ikernel,
                                        max_distance_deviation,
                                        dist_normalization_factor,
                                        epsilon);

        std::vector<double> Ecv_y = res.first;
        std::vector<int> excluded_vertices = res.second;

        // Computing a set difference between test_set and excluded_vertices
        std::set<int> valid_test_set;
        for (const auto& vertex : test_set) {
            if (std::find(excluded_vertices.begin(), excluded_vertices.end(), vertex) == excluded_vertices.end()) {
                valid_test_set.insert(vertex);
            }
        }

        // Checking if valid_test_set is empty
        if (valid_test_set.empty()) {
            continue;  // Skip this iteration if no valid test vertices
        }

        // Computing cross-validation Rf_error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double Rf_error = std::abs(Ecv_y[vertex] - y[vertex]);

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = Rf_error;
            } else {
                cv_error[vertex] += Rf_error;
            }
            cv_error_count[vertex]++;
        }
    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV Rf_error, leaving NaN for vertices with no estimates
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex] > 0) {
            cv_error[vertex] /= cv_error_count[vertex];
        }
    }

    return cv_error;
}

/**
 * @brief Performs parallel cross-validation for Gaussian Process Latent Model (PGMALO)
 *
 * @details This function implements a parallel cross-validation procedure for PGMALO
 * using TBB (Threading Building Blocks). It divides the data into folds and
 * estimates the model's predictive performance by computing cross-validation errors
 * for each vertex in the graph.
 *
 * The function uses parallel execution for the main cross-validation iterations while
 * ensuring thread safety through appropriate synchronization mechanisms. It handles
 * the creation of test sets, weight calculations, and Rf_error averaging.
 *
 * @param path_graph Reference to the path graph structure containing adjacency lists
 *                   and other graph-related data
 * @param y Vector of observed values at each vertex
 * @param n_CVs Number of cross-validation iterations to perform
 * @param n_CV_folds Number of folds for cross-validation
 * @param seed Random seed for reproducibility (0 for time-based seed)
 * @param ikernel Kernel function identifier for the PGMALO
 * @param max_distance_deviation Maximum allowed deviation in distance calculations
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical tolerance parameter (default: 1e-8)
 *
 * @return std::vector<double> Vector containing average cross-validation errors
 *         for each vertex. Vertices that were never in a test set or were excluded
 *         will have NaN values.
 *
 * @note The function uses atomic operations and mutex locks to ensure thread safety
 *       when updating shared data structures.
 *
 * @Rf_warning The size of input vector y must match the number of vertices in path_graph.
 *          The function assumes that path_graph contains valid adjacency lists for
 *          all vertices.
 *
 * @throws std::invalid_argument If input dimensions are inconsistent or invalid
 *
 * @see pgmalo_with_cv_weights
 * @see path_graph_plm_t
 *
 * @example
 * ```cpp
 * path_graph_plm_t graph;  // Initialize your graph
 * std::vector<double> observations{1.0, 2.0, 3.0, 4.0, 5.0};
 * int n_cvs = 10;
 * int n_folds = 5;
 * unsigned int seed = 42;
 *
 * auto cv_errors = pgmalo_cv_parallel(graph, observations, n_cvs, n_folds,
 *                                  seed, 1, 2, 1.01, 1e-8);
 * ```
 */
std::vector<double> pgmalo_cv_parallel(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_CVs,
    int n_CV_folds,
    unsigned int seed,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    int n_vertices = static_cast<int>(y.size());

    auto cv_error = std::vector<double>(n_vertices, std::numeric_limits<double>::quiet_NaN());
    std::vector<std::atomic<int>> cv_error_count(n_vertices);
    std::vector<std::mutex> vertex_mutexes(n_vertices);

    // Check if a seed was provided
    if (seed == 0) {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Creating a set version of the adjacency matrix of the graph
    std::vector<std::set<int>> set_graph(n_vertices);

    // First, populate the sets sequentially to avoid any race conditions
    for (int i = 0; i < n_vertices; ++i) {
        set_graph[i].insert(path_graph.adj_list[i].begin(),
                          path_graph.adj_list[i].end());
    }

    int fold_size = n_vertices / n_CV_folds;

    // Generate all CV iterations in advance
    std::vector<std::vector<int>> all_test_sets(n_CVs);
    {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> uni(0, n_vertices - 1);

        for (int cv = 0; cv < n_CVs; ++cv) {
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    test_set.insert(uni(rng));
                }
            }
            all_test_sets[cv].insert(all_test_sets[cv].end(),
                                   test_set.begin(),
                                   test_set.end());
        }
    }

    // Parallel execution of CV iterations
    std::vector<int> cv_indices(n_CVs);
    std::iota(cv_indices.begin(), cv_indices.end(), 0);

    std::for_each(GFLOW_EXEC_POLICY,
                  cv_indices.begin(),
                  cv_indices.end(),
                  [&](int cv) {
        std::vector<double> weights(n_vertices, 1.0);
        const auto& test_set = all_test_sets[cv];

        // Set weights for test set vertices to 0
        for (const auto& vertex : test_set) {
            weights[vertex] = 0.0;
        }

        // Estimate conditional expectation
        auto res = pgmalo_with_cv_weights(path_graph,
                                      y,
                                      weights,
                                      ikernel,
                                      max_distance_deviation,
                                      dist_normalization_factor,
                                      epsilon);

        std::vector<double> Ecv_y = res.first;
        std::vector<int> excluded_vertices = res.second;

        // Compute valid test set
        std::vector<int> valid_test_vertices;
        std::copy_if(test_set.begin(),
                    test_set.end(),
                    std::back_inserter(valid_test_vertices),
                    [&](int vertex) {
                        return std::find(excluded_vertices.begin(),
                                       excluded_vertices.end(),
                                       vertex) == excluded_vertices.end();
                    });

        // Update cross-validation errors
        for (const auto& vertex : valid_test_vertices) {
            double Rf_error = std::abs(Ecv_y[vertex] - y[vertex]);

            std::lock_guard<std::mutex> lock(vertex_mutexes[vertex]);
            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = Rf_error;
            } else {
                cv_error[vertex] += Rf_error;
            }
            cv_error_count[vertex]++;
        }
    });

    // Compute average CV Rf_error
    for (int i = 0; i < n_vertices; ++i) {
        if (cv_error_count[i] > 0) {
            cv_error[i] /= cv_error_count[i];
        }
    }

    return cv_error;
}


/**
 * @brief Computes Graph Path Linear Models (PGMALO) with cross-validation weights for vertex predictions
 *
 * @details This function implements a variant of PGMALO specifically designed for cross-validation,
 * where test set vertices are excluded from model fitting by setting their weights to 0.
 * The function processes paths in the graph to fit local linear models while respecting
 * the weight-based exclusions and tracking vertices that cannot be modeled due to
 * insufficient data points.
 *
 * The algorithm operates in three phases:
 * 1. Precomputes models for all valid paths (those with sufficient non-zero weight points)
 * 2. Processes each vertex to find its best predictive model among valid paths
 * 3. Performs model averaging for final predictions
 *
 * A vertex is considered "excluded" if no valid model can be fit for it due to:
 * - No paths containing the vertex
 * - Insufficient non-zero weight points in all available paths
 * - No valid models found in the model cache
 *
 * @param path_graph The path graph structure containing graph topology and path information
 * @param y Vector of response variables for each vertex
 * @param weights Binary weight vector (0.0 for test set, 1.0 for training set)
 * @param ikernel Kernel function identifier for weight computation
 * @param max_distance_deviation Maximum allowed deviation from optimal position in path
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical tolerance for zero comparisons (default: 1e-8)
 *
 * @return std::pair containing:
 *         - First: Vector of predicted values for each vertex
 *         - Second: Vector of indices of excluded vertices (those without valid models)
 *
 * @note For excluded vertices, the original y values are retained in the prediction vector
 *
 * @pre weights must contain only values 0.0 (test set) or 1.0 (training set)
 * @pre path_graph must be properly initialized with valid paths
 * @pre y and weights vectors must have size equal to number of vertices in path_graph
 *
 * Key parameters:
 * - Minimum number of non-zero weight points required for model fitting: 3
 * - Weights are treated as binary: values > epsilon are set to 1.0, others to 0.0
 *
 * @throws Rf_error If no valid paths are found for a vertex
 *
 * Example usage:
 * @code
 * // Create path graph and prepare data
 * path_graph_plm_t path_graph = create_path_graph_plm(adj_list, weight_list, h);
 * std::vector<double> y = {0.0, 1.0, 0.0, 1.0};  // Binary response
 * std::vector<double> weights = {1.0, 0.0, 1.0, 1.0};  // One vertex in test set
 *
 * // Compute PGMALO with cross-validation
 * auto [predictions, excluded] = pgmalo_with_cv_weights(
 *     path_graph, y, weights, 1, 2, 1.01, 1e-8
 * );
 * @endcode
 *
 * @see path_graph_plm_t
 * @see create_path_graph_plm
 *
 * @Rf_warning This function modifies the provided weights to be strictly binary (0.0 or 1.0)
 * @Rf_warning Performance depends on the number of valid paths and non-zero weight points
 *
 * Time Complexity: O(V * P * N), where:
 * - V is the number of vertices
 * - P is the average number of paths per vertex
 * - N is the average path length
 *
 * Space Complexity: O(V * P + M), where M is the size of the model cache
 */
std::pair<std::vector<double>, std::vector<int>> pgmalo_with_cv_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    const int min_n_non_zero_points = 3; // Minimum number of points needed for fitting
    int h = path_graph.h;
    int desired_path_length = h + 1;
    int n_vertices = path_graph.vertex_paths.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Initialize cache and collect unique paths
    std::map<itriplet_t, lm_loocv_t> model_cache;

    // Initialize working vectors for model computation
    std::vector<double> y_path;
    std::vector<double> x_path;
    std::vector<double> bb_weights_path;
    std::vector<double> w_path;
    std::vector<double> d_path;
    std::vector<std::vector<double>> w_list;

    #ifdef _OPENMP
#pragma omp critical(kernel_init)
    #endif
    {
        initialize_kernel(ikernel, 1.0);
    }

    // Phase 1: Precompute models for all paths
    for (const auto& path_endpoints : path_graph.longest_paths) {
        const std::vector<int>& path = path_graph.shortest_paths.at(path_endpoints);
        int path_n_vertices = path.size();

        // Count non-zero weight points in this path
        int n_non_zero_points = 0;
        for (int i = 0; i < path_n_vertices; ++i) {
            if (weights[path[i]] > epsilon) {
                n_non_zero_points++;
            }
        }

        // Skip paths with insufficient non-zero weight points
        if (n_non_zero_points < min_n_non_zero_points) {
            continue;
        }

        // Resize working vectors for current path
        y_path.resize(path_n_vertices);
        x_path.resize(path_n_vertices);
        bb_weights_path.resize(path_n_vertices);
        w_path.resize(path_n_vertices);
        d_path.resize(path_n_vertices);
        w_list.resize(path_n_vertices);

        // Initialize path coordinates
        x_path[0] = 0;

        // Extract values and compute path distances
        for (int i = 0; i < path_n_vertices; ++i) {
            y_path[i] = y[path[i]];
            // Convert weights to binary (0.0 or 1.0)
            bb_weights_path[i] = (weights[path[i]] > epsilon) ? 1.0 : 0.0;

            if (i > 0) {
                auto neighbors = path_graph.adj_list[path[i - 1]];
                auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                int neighbor_i = it - neighbors.begin();
                x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i];
            }
        }

        // Compute kernel weights for each position
        for (int pos = 0; pos < path_n_vertices; ++pos) {
            // Calculate distances to reference point
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[pos]);
                max_dist = std::max(max_dist, d_path[i]);
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            // Normalize distances and compute kernel weights
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] /= max_dist;
            }

            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            // Normalize kernel weights and apply binary weights
            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i) {
                w_path[i] = (w_path[i] / total_w_path) * bb_weights_path[i];
            }

            w_list[pos].assign(w_path.begin(), w_path.end());
        }

        // Compute models for all positions
        auto models = predict_lms_1d_loocv(y_path, x_path, path, w_list, y_binary, epsilon);

        // Cache computed models
        for (int pos = 0; pos < path_n_vertices; ++pos) {
            model_cache.emplace(itriplet_t{path_endpoints.first, path_endpoints.second, pos},
                                models[pos]);
        }
    }

    // Phase 2: Process each vertex using precomputed models
    std::vector<double> Ey(n_vertices);
    std::vector<int> excluded_vertices;

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get all paths containing the vertex
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(desired_path_length, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            excluded_vertices.push_back(vertex_i);
            Ey[vertex_i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Find valid paths based on distance to middle point
        std::vector<size_t> valid_path_indices;
        valid_path_indices.reserve(vertex_paths.size());

        int mid_pt = (desired_path_length - 1) / 2;
        int min_dist_to_mid_pt = h;

        // Find minimum distance to middle point
        for (const auto& path : vertex_paths) {
            min_dist_to_mid_pt = std::min(min_dist_to_mid_pt, std::abs(path.second - mid_pt));
        }

        // Collect paths within allowed deviation and with sufficient non-zero weights
        for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            const auto& path = vertex_paths[path_i].first;

            // Count non-zero weight points in this path
            int n_non_zero_points = 0;
            for (int i = 0; i < path.size(); ++i) {
                if (weights[path[i]] > epsilon) {
                    n_non_zero_points++;
                }
            }

            if (n_non_zero_points >= min_n_non_zero_points &&
                std::abs(vertex_paths[path_i].second - mid_pt) <=
                min_dist_to_mid_pt + max_distance_deviation) {
                valid_path_indices.push_back(path_i);
            }
        }

        if (valid_path_indices.empty()) {
            excluded_vertices.push_back(vertex_i);
            Ey[vertex_i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Find best model among valid paths
        double best_loocv_error = std::numeric_limits<double>::infinity();
        lm_loocv_t best_model;
        bool found_valid_model = false;

        for (const auto& path_i : valid_path_indices) {
            const auto& path = vertex_paths[path_i].first;
            int position_in_path = vertex_paths[path_i].second;

            auto it = model_cache.find(itriplet_t{path[0], path.back(), position_in_path});
            if (it != model_cache.end()) {
                found_valid_model = true;
                if (it->second.loocv_at_ref_vertex < best_loocv_error) {
                    best_loocv_error = it->second.loocv_at_ref_vertex;
                    best_model = it->second;
                }
            }
        }

        if (!found_valid_model) {
            excluded_vertices.push_back(vertex_i);
            Ey[vertex_i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Use the best model for prediction
        Ey[vertex_i] = best_model.predicted_value;
        if (y_binary) {
            Ey[vertex_i] = std::clamp(Ey[vertex_i], 0.0, 1.0);
        }
    }

    // Phase 3: Perform model averaging
    std::vector<double> final_Ey = Ey;  // Create copy for averaging
    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        // Skip excluded vertices
        if (std::find(excluded_vertices.begin(), excluded_vertices.end(), vertex_i) != excluded_vertices.end()) {
            continue;
        }

        const auto& vertex_path_graph_info = path_graph.vertex_paths[vertex_i];
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(desired_path_length, path_graph.shortest_paths);

        double total_weight = 0.0;
        double weighted_sum = 0.0;

        for (const auto& path_pair : vertex_paths) {
            const auto& path = path_pair.first;
            int pos = path_pair.second;

            // Check if this path has a valid model
            auto it = model_cache.find(itriplet_t{path[0], path.back(), pos});
            if (it != model_cache.end()) {
                const auto& model = it->second;
                double weight = model.w_values[model.ref_index];
                if (weight > epsilon) {
                    total_weight += weight;
                    weighted_sum += weight * Ey[vertex_i];
                }
            }
        }

        if (total_weight > epsilon) {
            final_Ey[vertex_i] = weighted_sum / total_weight;
            if (y_binary) {
                final_Ey[vertex_i] = std::clamp(final_Ey[vertex_i], 0.0, 1.0);
            }
        }
    }

    return std::make_pair(final_Ey, excluded_vertices);
}
