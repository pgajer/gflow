#ifdef _OPENMP
    #include <omp.h>
#endif

// Undefine R's match macro if it exists
#ifdef match
    #undef match
#endif

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length
#undef eval

#include <execution>
#include <atomic>
#include <mutex>
#include <omp.h>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <unordered_map>
#include <map>

#include "pgmalog.hpp"
#include "pgmalo.hpp"
#include "ulogit.hpp"
#include "ulm.hpp"
#include "memory_utils.hpp"
#include "progress_utils.hpp"
#include "error_utils.h"
#include "sampling.h" // for C_runif_simplex()
#include "pglm.h"
#include "msr2.h"
#include "path_graphs.hpp"
#include "cpp_utils.hpp"                 // for elapsed_time
#include "SEXP_cpp_conversion_utils.hpp"
#include "kernels.h"
#include "1D_linear_models.h"
#include "predictive_errors.hpp"

bb_cri_t pgmalog_bb_cri(const path_graph_plm_t& path_graph,
                     const std::vector<double>& y,
                     double p,
                     int n_bb,
                     int max_distance_deviation,
                     bool use_median,  // Added missing parameter
                     int ikernel,
                     double dist_normalization_factor,
                     double epsilon);

std::vector<double> pgmalog_cv(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_CVs,
    int n_CV_folds,
    unsigned int seed,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon);

std::vector<double> pgmalog_cv_parallel(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_CVs,
    int n_CV_folds,
    unsigned int seed,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon);

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>
create_chain_graph(const std::vector<double>& x);

path_graph_plm_t create_path_graph_plm(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h);

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
path_graph_plm_t sexp_to_path_graph_plm(SEXP s_path_graph);

extern "C" {
    SEXP S_upgmalog(SEXP s_x,
                   SEXP s_y,
                   SEXP s_y_true,
                   SEXP s_use_median,
                   SEXP s_h_min,
                   SEXP s_h_max,
                   SEXP s_p,
                   SEXP s_n_bb,
                   SEXP s_bb_max_distance_deviation,
                   SEXP s_n_CVs,
                   SEXP s_n_CV_folds,
                   SEXP s_seed,
                   SEXP s_ikernel,
                   SEXP s_n_cores,
                   SEXP s_dist_normalization_factor,
                   SEXP s_epsilon,
                   SEXP s_verbose);
}



/**
 * @brief Computes Graph Path Linear Models (PGMALOG) with cross-validation weights for vertex predictions
 *
 * @details This function implements a variant of PGMALOG specifically designed for cross-validation,
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
 * // Compute PGMALOG with cross-validation
 * auto [predictions, excluded] = pgmalog_with_cv_weights(
 *     path_graph, y, weights, 1, 2, 1.01, 1e-8
 * );
 * @endcode
 *
 * @see path_graph_plm_t
 * @see create_path_graph_plm
 *
 * @warning This function modifies the provided weights to be strictly binary (0.0 or 1.0)
 * @warning Performance depends on the number of valid paths and non-zero weight points
 *
 * Time Complexity: O(V * P * N), where:
 * - V is the number of vertices
 * - P is the average number of paths per vertex
 * - N is the average path length
 *
 * Space Complexity: O(V * P + M), where M is the size of the model cache
 */
std::pair<std::vector<double>, std::vector<int>> pgmalog_with_cv_weights(
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

    #pragma omp critical(kernel_init)
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



/**
* @brief Performs Bayesian bootstrap sampling for graph path linear models with path-specific weights
*
* This function implements a modified version of graph path linear models (PGMALOG) with Bayesian
* bootstrap sampling. For each path in the graph, it generates path-specific weights from a
* Dirichlet distribution parameterized by the path's kernel weights. This approach ensures that
* the Bayesian bootstrap estimates are centered around the unweighted graph path linear model.
*
* The algorithm follows these steps for each vertex:
* 1. Identifies all valid paths containing the vertex
* 2. For each valid path:
*    - Computes kernel weights based on distances
*    - Generates path-specific weights from Dirichlet(kernel_weights)
*    - Fits a weighted linear model using the sampled weights
*    - Computes LOOCV error for model selection
* 3. Selects the best model based on LOOCV error
* 4. Computes predictions using either direct model prediction or weighted average
*    of predictions from neighboring vertices
*
* @param path_graph The path graph structure containing adjacency lists, weights, and path information
* @param y Vector of response variables for each vertex
* @param ikernel Integer specifying the kernel type for weight computation
* @param max_distance_deviation Maximum allowed deviation from the minimum distance to mid-point
*                              when selecting valid paths. If > 0, enables weighted averaging
*                              of predictions from neighboring vertices
* @param min_required_paths Minimum number of valid paths required for each vertex
* @param dist_normalization_factor Factor used to normalize distances (default: 1.01)
* @param epsilon Small positive number for numerical stability (default: 1e-8)
*
* @return std::pair containing:
*         - First element: Vector of predicted values for each vertex
*         - Second element: Vector of LOOCV errors for each vertex
*
* @throws Rf_error if:
*         - max_distance_deviation is negative
*         - min_required_paths is less than 1
*         - Length of y doesn't match number of vertices
*         - No paths found for a vertex
*         - No valid paths found within deviation limits
*
* @note The Bayesian bootstrap is implemented using path-specific Dirichlet sampling,
*       where the parameters of the Dirichlet distribution are the kernel weights
*       computed for each path. This ensures that:
*       1. The weights are properly normalized within each path
*       2. The bootstrap estimates are centered around the original weighted model
*       3. The uncertainty quantification is path-specific
*
* @see C_rsimplex For the Dirichlet sampling implementation
* @see predict_lm_1d_loocv For the weighted linear model fitting
* @see path_graph_plm_t For the path graph data structure
*/
std::pair<std::vector<double>, std::vector<double>> pgmalog_with_bb_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int ikernel,
    int max_distance_deviation,
    int min_required_paths,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    // Input validation
    if (max_distance_deviation < 0) {
        Rf_error("max_distance_deviation must be non-negative");
    }
    if (min_required_paths < 1) {
        Rf_error("min_required_paths must be at least 1");
    }

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    std::vector<double> Ey(n_vertices);
    std::vector<double> errors(n_vertices);
    std::vector<lm_loocv_t> best_models(n_vertices);

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get all paths containing the vertex
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            Rf_error("No paths found for vertex %d", vertex_i);
        }

        // Find paths within allowed deviation
        std::vector<size_t> valid_path_indices;
        int mid_pt = (h - 1) / 2;
        int min_dist_to_mid_pt = h;

        // First pass: find minimum distance
        for (const auto& path : vertex_paths) {
            int dist_to_mid_pt = std::abs(path.second - mid_pt);
            min_dist_to_mid_pt = std::min(min_dist_to_mid_pt, dist_to_mid_pt);
        }

        // Second pass: collect valid paths
        for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            int dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt);
            if (dist_to_mid_pt <= min_dist_to_mid_pt + max_distance_deviation) {
                valid_path_indices.push_back(path_i);
            }
        }

        if (valid_path_indices.empty()) {
            Rf_error("No valid paths found for vertex %d within deviation limits", vertex_i);
        }

        // For each valid path, compute the LOOCV model
        double best_loocv_error = std::numeric_limits<double>::infinity();
        lm_loocv_t best_model;

        for (const auto& path_i : valid_path_indices) {
            std::vector<int> path = vertex_paths[path_i].first;
            int path_n_vertices = path.size();
            int position_in_path = vertex_paths[path_i].second;

            // Prepare data for linear model
            std::vector<double> y_path(path_n_vertices);
            std::vector<double> x_path(path_n_vertices);
            std::vector<double> d_path(path_n_vertices);
            x_path[0] = 0;

            // Extract y values and compute distances
            for (int i = 0; i < path_n_vertices; ++i) {
                y_path[i] = y[path[i]];
                if (i > 0) {
                    auto neighbors = path_graph.adj_list[path[i - 1]];
                    auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                    auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                    int neighbor_i = it - neighbors.begin();
                    x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i];
                }
            }

            // Compute kernel weights
            std::vector<double> w_path(path_n_vertices);
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[position_in_path]);
                max_dist = std::max(max_dist, d_path[i]);
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] /= max_dist;
            }

            #pragma omp critical(kernel_init)
            {
                initialize_kernel(ikernel, 1.0);
            }
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] /= total_w_path;

            // Generate Dirichlet samples using kernel weights as parameters
            std::vector<double> bb_weights_path(path_n_vertices);
            C_rsimplex(w_path.data(), &path_n_vertices, bb_weights_path.data());

            // Fit LOOCV model with sampled weights
            auto model = predict_lm_1d_loocv(y_path, x_path, bb_weights_path, path, position_in_path, y_binary, epsilon);

            // Update best model if this one has lower LOOCV error
            if (model.loocv_at_ref_vertex < best_loocv_error) {
                best_loocv_error = model.loocv_at_ref_vertex;
                best_model = model;
            }
        }

        best_models[vertex_i] = best_model;
        Ey[vertex_i] = best_model.predicted_value;
        errors[vertex_i] = best_model.loocv_at_ref_vertex;
    }

    // Handle prediction based on max_distance_deviation
    if (max_distance_deviation > 0) {
        for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
            lm_loocv_t best_model = best_models[vertex_i];
            // Get vertices of the best model path
            std::vector<int> path_vertex_indices = best_model.vertex_indices;
            std::vector<double> path_vertex_weights = best_model.w_values;
            int path_length = path_vertex_weights.size();

            // for each vertex from path_vertex_indices use its best model to predict the value of that model at 'vertex_i'
            // Of course for 'vertex_i' we have already computed the predicted value of its best model: best_model.predicted_value = Ey[vertex_i]
            double weight = path_vertex_weights[best_model.ref_index];
            double total_weight = weight;
            double vertex_Ey = weight * Ey[vertex_i];

            for (int i = 0; i < path_length; ++i) {
                int index =  path_vertex_indices[i];
                if (index != vertex_i) {
                    lm_loocv_t i_best_model = best_models[index];
                    // Check if vertex_i is in the model's training set
                    auto it = std::find(i_best_model.vertex_indices.begin(),
                                        i_best_model.vertex_indices.end(),
                                        vertex_i);
                    if (it != i_best_model.vertex_indices.end()) {
                        double i_Ey = i_best_model.predict(vertex_i);
                        weight = path_vertex_weights[i];
                        total_weight += weight;
                        vertex_Ey += weight * i_Ey;
                    }
                }
            }
            vertex_Ey /= total_weight;
            Ey[vertex_i] = vertex_Ey;
        }
    }

    return std::make_pair(Ey, errors);
}

/**
 * @brief Path Graph Linear Model with adaptive neighborhood size selection
 *
 * @details Fits a PGMALOG by testing different neighborhood sizes (h values) and selecting
 * the optimal size based on cross-validation errors. For each h value, constructs a path
 * graph and computes both global and local predictions. Optionally computes bootstrap
 * credible intervals.
 *
 *  This function performs several steps:
 *    1. For each odd h in [h_min, h_max]:
 *       - Constructs h-hop neighborhood Path Linear Model (PLM) graphs
 *       - Computes leave-one-out cross-validation (LOOCV) errors
 *   2. Global optimization:
 *       - Finds optimal h that minimizes mean CV error
 *       - Computes conditional expectations using optimal h
 *   3. Local optimization:
 *       - For each vertex, finds optimal h minimizing its CV error
 *       - Computes vertex-specific conditional expectations
 *   4. Optional: Computes Bayesian bootstrap credible intervals
 *       - Estimates central tendency of conditional expectations
 *       - Computes credible interval bounds at specified probability level
 *
 * @param neighbors Adjacency lists for each vertex
 * @param edge_lengths Corresponding edge lengths for each adjacency
 * @param y Response variables for each vertex
 * @param y_true True values for error computation (optional)
 * @param use_median Use median instead of mean for bootstrap intervals
 * @param h_min Minimum neighborhood size (must be odd)
 * @param h_max Maximum neighborhood size (must be odd)
 * @param p Confidence level for bootstrap intervals (0 < p < 1)
 * @param n_bb Number of bootstrap iterations (0 for no bootstrap)
 * @param bb_max_distance_deviation The maximal distance deviation from the min distance from the reference point applied to the bootstrap samples of y for the optimal h. To pick the optimal h max_distance_deviation is set to (h-1)/2. If bb_max_distance_deviation = -1, we set it to (h-1)/2 in the bootstrap loop
 * @param ikernel Kernel function identifier
 * @param n_cores Number of cores for parallel computation
 * @param dist_normalization_factor Distance normalization factor
 * @param epsilon Numerical stability threshold
 * @param seed Random seed for reproducibility
 * @param verbose Enable progress reporting
 *
 * @return pgmalo_t structure containing:
 *         - Optimal h value and corresponding graph
 *         - Global and local predictions
 *         - Cross-validation errors
 *         - Bootstrap credible intervals (if requested)
 */
pgmalo_t pgsmalog(const std::vector<std::vector<int>>& neighbors,
                                 const std::vector<std::vector<double>>& edge_lengths,
                                 const std::vector<double>& y,
                                 const std::vector<double>& y_true,
                                 bool use_median = false,
                                 int h_min = 3,
                                 int h_max = 31,
                                 double p = 0.95,
                                 int n_bb = 500,
                                 int bb_max_distance_deviation = 2,
                                 int n_CVs = 0,
                                 int n_CV_folds = 10,
                                 unsigned int seed = 0,
                                 int ikernel = 1,
                                 int n_cores = 1,
                                 double dist_normalization_factor = 1.01,
                                 double epsilon = 1e-15,
                                 bool verbose = true) {

    auto total_ptm = std::chrono::steady_clock::now();
    auto ptm = std::chrono::steady_clock::now();  // Start timing

    const int n_vertices = static_cast<int>(y.size());
    const int n_h_values = (h_max - h_min) / 2 + 1;

    #define DEBUG__pgsmalo 0
    #if DEBUG__pgsmalo
    Rprintf("\nIn pgsmalo()\n");
    Rprintf("Number of vertices: %d\n", n_vertices);
    print_vect_vect(neighbors, "neighbors");
    print_vect_vect(edge_lengths, "edge_lengths");
    #endif

    // Initialize results structure
    pgmalo_t results;
    results.graphs.resize(n_h_values);
    results.h_cv_errors.resize(n_h_values);
    results.h_values.resize(n_h_values);
    results.h_predictions.resize(n_h_values);

    // Storage for intermediate results
    //std::unordered_map<int, std::vector<double>> errors_map;
    //errors_map.reserve(n_h_values);
    std::unordered_map<int, std::vector<double>> predictions_map;
    predictions_map.reserve(n_h_values);

    // Initialize uniform weights
    std::vector<double> weights(n_vertices, 1.0);

    // Process each h value
    for (int i = 0, h = h_min; h <= h_max; h += 2, i++) {

        results.h_values[i] = h;

        //int local_max_distance_deviation = (h - 1) / 2;
        int local_max_distance_deviation = 0;

        if (verbose) {
            Rprintf("\nProcessing h = %d (%d/%d)\n", h, i + 1, n_h_values);
        }

        // Create path graph
        if (verbose) {
            Rprintf("\tCreating path graph ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto path_graph = create_path_graph_plm(neighbors, edge_lengths, h);

        if (verbose) elapsed_time(ptm, "DONE");

        // Compute predictions and errors
        if (verbose) {
            Rprintf("\tComputing predictions and errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto [predictions, errors] = pgmalog(path_graph,
                                            y,
                                            weights,
                                            ikernel,
                                            local_max_distance_deviation,
                                            dist_normalization_factor,
                                            epsilon,
                                            verbose);

        if (n_CVs) { // overwriting errors with CV errors

            errors = pgmalog_cv_parallel(path_graph,
                                        y,
                                        n_CVs,
                                        n_CV_folds,
                                        seed,
                                        ikernel,
                                        local_max_distance_deviation,
                                        dist_normalization_factor,
                                        epsilon);

            // Removing NaN elements from errors
            errors.erase(std::remove_if(errors.begin(), errors.end(),
                                        [](double x) { return std::isnan(x); }),
                         errors.end());
        }

        if (verbose) elapsed_time(ptm, "DONE");

        #if DEBUG__pgmalog
        print_vect(predictions,"predictions");
        print_vect(errors,"errors");
        #endif

        // Calculate mean error
        if (verbose) {
            Rprintf("\tCalculating mean error ... ");
            ptm = std::chrono::steady_clock::now();
        }

        const double mean_error = std::accumulate(errors.begin(), errors.end(), 0.0) / n_vertices;
        results.h_cv_errors[i] = mean_error;
        results.graphs[i] = std::move(path_graph);
        //errors_map[h] = std::move(errors);

        results.h_predictions[i] = predictions;

        predictions_map[h] = std::move(predictions);

        if (verbose) elapsed_time(ptm, "DONE");
    }

    // Find optimal h
    auto min_it = std::min_element(results.h_cv_errors.begin(), results.h_cv_errors.end());
    const int opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * opt_h_idx;

    #if 0
    // Compute local predictions
    if (verbose) {
        Rprintf("Computing local predictions ... ");
        ptm = std::chrono::steady_clock::now();
    }
    results.opt_local_predictions.resize(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        double min_error = std::numeric_limits<double>::infinity();
        int local_opt_h = h_min;

        for (int h_idx = 0; h_idx < n_h_values; h_idx++) {
            const int h = results.h_values[h_idx];
            if (errors_map[h][i] < min_error) {
                min_error = errors_map[h][i];
                local_opt_h = h;
            }
        }
        results.opt_local_predictions[i] = predictions_map[local_opt_h][i];
    }
    if (verbose) elapsed_time(ptm, "DONE");
    #endif

    // Store optimal results
    results.graph = std::move(results.graphs[opt_h_idx]);
    results.predictions = std::move(predictions_map[results.opt_h]);

    // Compute true errors if available
    if (!y_true.empty()) {
        if (y_true.size() != n_vertices) {
            Rf_error("y_true size (%zu) does not match number of vertices (%d)",
                    y_true.size(), n_vertices);
        }

        if (verbose) {
            Rprintf("Computing true errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        results.true_errors.resize(n_vertices);
        #pragma omp parallel for if(n_cores > 1) schedule(dynamic)
        for (int i = 0; i < n_vertices; i++) {
            results.true_errors[i] = std::abs(y_true[i] - results.predictions[i]);
        }

        if (verbose) elapsed_time(ptm, "DONE");
    }

    // Compute bootstrap intervals if requested
    if (n_bb > 0) {
        if (verbose) {
            Rprintf("Computing bootstrap intervals (n=%d) ... ", n_bb);
            ptm = std::chrono::steady_clock::now();
        }

        int local_max_distance_deviation = 0;
        if (bb_max_distance_deviation == -1) {
            local_max_distance_deviation = (results.opt_h - 1) / 2;
        } else {
            local_max_distance_deviation = bb_max_distance_deviation;
        }

        bb_cri_t bb_cri_results = pgmalog_bb_cri(results.graph,
                                              y,
                                              p,
                                              n_bb,
                                              local_max_distance_deviation,
                                              use_median,
                                              ikernel,
                                              dist_normalization_factor,
                                              epsilon);

        results.bb_predictions = std::move(bb_cri_results.bb_Ey);
        results.ci_lower = std::move(bb_cri_results.cri_L);
        results.ci_upper = std::move(bb_cri_results.cri_U);

        if (verbose) elapsed_time(ptm, "DONE");
    }

    if (verbose) {
        Rprintf("\nOptimal h: %d\n", results.opt_h);
        elapsed_time(total_ptm, "Total time");
    }

    return results;
}


/**
 * @brief Optimizes univariate path linear model estimation using path locally weighted linear models
 *
 * @details This is a specialized version of pgmalog for univariate data.
 *          It constructs a chain graph from the ordered x values and applies path linear model estimation.
 *          The function performs several steps:
 *          1. Creates a chain graph from ordered x values
 *          2. For each odd h in [h_min, h_max]:
 *             - Constructs h-hop neighborhood Path Linear Model (PLM) graphs
 *             - Computes leave-one-out cross-validation (LOOCV) errors
 *          3. Global optimization:
 *             - Finds optimal h that minimizes mean CV error
 *             - Computes conditional expectations using optimal h
 *          4. Local optimization:
 *             - For each vertex, finds optimal h minimizing its CV error
 *             - Computes vertex-specific conditional expectations
 *          5. Optional: Computes Bayesian bootstrap credible intervals
 *
 * @param x Vector of ordered x values
 * @param y Observed y values corresponding to x
 * @param y_true True y values for error calculation (optional)
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param use_median Use median instead of mean for central tendency (default: false)
 * @param h_min Minimum neighborhood size to consider (default: 3, must be odd)
 * @param h_max Maximum neighborhood size to consider (default: 31, must be odd)
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 500, 0 to skip)
 * @param bb_max_distance_deviation The maximal distance deviation from the min distance from the reference point applied to the bootstrap samples of y for the optimal h. To pick the optimal h max_distance_deviation is set to (h-1)/2. If bb_max_distance_deviation = -1, we set it to (h-1)/2 in the bootstrap loop
 * @param ikernel Kernel function selector (default: 1)
 * @param n_cores Number of cores for parallel computation (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages (default: true)
 *
 * @pre x and y must have the same size
 * @pre x and y cannot be empty
 * @pre h_min and h_max must be odd numbers
 * @pre y_true, if provided, must have the same size as y
 * @pre p must be in range (0,1)
 *
 * @return pgmalo_t struct containing:
 *   - graphs: Vector of PLM graphs for each h value
 *   - h_values: Vector of used h values (odd numbers from h_min to h_max)
 *   - cv_errors: Mean cross-validation errors for each h
 *   - true_errors: Mean absolute deviation from true values (if provided)
 *   - opt_h: Optimal h value minimizing mean CV error
 *   - graph: Graph structure for optimal h
 *   - predictions: Global conditional expectations using optimal h
 *   - local_predictions: Vertex-specific conditional expectations using locally optimal h
 *   - bb_predictions: Bootstrap central tendency estimates (if n_bb > 0)
 *   - opt_ci_lower, opt_ci_upper: Credible interval bounds (if n_bb > 0)
 *
 * @see pgmalog
 * @see create_chain_graph
 *
 * @note This function creates a chain graph where each vertex is connected to its immediate neighbors
 *       in the ordered sequence of x values.
 */
pgmalo_t upgmalog(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& y_true,
                  bool use_median,
                  int h_min,
                  int h_max,
                  double p,
                  int n_bb,
                  int bb_max_distance_deviation,
                  int n_CVs,
                  int n_CV_folds,
                  unsigned int seed,
                  int ikernel,
                  int n_cores,
                  double dist_normalization_factor,
                  double epsilon,
                  bool verbose) {

    // Add input validation
    if (x.empty() || y.empty()) {
        Rf_error("Input vectors x and y cannot be empty");
    }
    if (x.size() != y.size()) {
        Rf_error("Input vectors x and y must have the same size");
    }

    auto [x_graph, x_edge_lengths] = create_chain_graph(x);

    return pgsmalog(x_graph,
                    x_edge_lengths,
                    y,
                    y_true,
                    use_median,
                    h_min,
                    h_max,
                    p,
                    n_bb,
                    bb_max_distance_deviation,
                    n_CVs,
                    n_CV_folds,
                    seed,
                    ikernel,
                    n_cores,
                    dist_normalization_factor,
                    epsilon,
                    verbose);
}

/**
 * @brief R interface for univariate path linear model optimization using k-path locally weighted models
 *
 * @details This function provides an R interface to upgmalog.
 *          It converts R inputs to C++ types, calls the core implementation, and converts results
 *          back to R objects. The function performs:
 *          1. Input validation and conversion from R to C++
 *          2. Execution of univariate path linear model optimization
 *          3. Conversion of results to an R list
 *
 * @param s_x R numeric vector of ordered x values
 * @param s_y R numeric vector of observed y values
 * @param s_y_true R numeric vector of true y values (optional)
 * @param s_use_median R logical for using median instead of mean
 * @param s_h_min R integer for minimum neighborhood size (must be odd)
 * @param s_h_max R integer for maximum neighborhood size (must be odd)
 * @param s_p R numeric for probability level of credible intervals
 * @param s_n_bb R integer for number of bootstrap iterations
 * @param s_ikernel R integer for kernel function selection
 * @param s_n_cores R integer for number of cores to use
 * @param s_dist_normalization_factor R numeric for distance normalization
 * @param s_epsilon R numeric for numerical stability
 * @param s_verbose R logical for progress messages
 *
 * @return An R list containing:
 *   - h_values: Integer vector of used h values
 *   - opt_h: Optimal h value (numeric)
 *   - predictions: Numeric vector of conditional expectation estimates
 *   - opt_local_predictions: Numeric vector of conditional expectation estimates using locally (vertex-wise) optimal h values
 *   - opt_bb_predictions: Numeric vector of bootstrap estimates (if n_bb > 0)
 *   - opt_ci_lower: Numeric vector of lower credible bounds (if n_bb > 0)
 *   - opt_ci_upper: Numeric vector of upper credible bounds (if n_bb > 0)
 *   - h_cv_errors: Numeric vector of cross-validation errors
 *   - true_error: Mean true error (if y_true provided)
 *   - opt_graph_adj_list: List of adjacency lists for optimal graph
 *   - opt_graph_edge_lengths: List of edge lengths for optimal graph
 *
 * @note This function uses R's memory protection mechanism via PROTECT/UNPROTECT
 *
 * @seealso upgmalog
 *
 * @example
 * \dontrun{
 * # R example:
 * result <- .Call("S_upgmalog",
 *                 x = as.numeric(1:100),
 *                 y = rnorm(100),
 *                 y_true = NULL,
 *                 use_median = FALSE,
 *                 h_min = 3L,
 *                 h_max = 31L,
 *                 p = 0.95,
 *                 n_bb = 500L,
 *                 ikernel = 1L,
 *                 n_cores = 1L,
 *                 dist_normalization_factor = 1.01,
 *                 epsilon = 1e-15,
 *                 verbose = TRUE)
 * }
 */
SEXP S_upgmalog(SEXP s_x,
               SEXP s_y,
               SEXP s_y_true,
               SEXP s_use_median,
               SEXP s_h_min,
               SEXP s_h_max,
               SEXP s_p,
               SEXP s_n_bb,
               SEXP s_bb_max_distance_deviation,
               SEXP s_n_CVs,
               SEXP s_n_CV_folds,
               SEXP s_seed,
               SEXP s_ikernel,
               SEXP s_n_cores,
               SEXP s_dist_normalization_factor,
               SEXP s_epsilon,
               SEXP s_verbose) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(s_y_true) == n_points) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    //int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    bool use_median = LOGICAL(s_use_median)[0];
    int h_min = INTEGER(s_h_min)[0];
    int h_max = INTEGER(s_h_max)[0];
    double p = REAL(s_p)[0];
    int n_bb = INTEGER(s_n_bb)[0];
    int bb_max_distance_deviation = INTEGER(s_bb_max_distance_deviation)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    auto cpp_results = upgmalog(x,
                                y,
                                y_true,
                                use_median,
                                h_min,
                                h_max,
                                p,
                                n_bb,
                                bb_max_distance_deviation,
                                n_CVs,
                                n_CV_folds,
                                seed,
                                ikernel,
                                n_cores,
                                dist_normalization_factor,
                                epsilon,
                                verbose);
    // Creating return list
    const int N_COMPONENTS = 13;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(cpp_results.h_values)); n_protected++;

    if (!cpp_results.h_cv_errors.empty() && n_points > 0) {
        SET_VECTOR_ELT(result, 1, convert_vector_double_to_R(cpp_results.h_cv_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    SEXP s_opt_h_idx = PROTECT(allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h_idx)[0] = cpp_results.opt_h_idx + 1;
    SET_VECTOR_ELT(result, 2, s_opt_h_idx);

    SEXP s_opt_h = PROTECT(allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h)[0] = cpp_results.opt_h;
    SET_VECTOR_ELT(result, 3, s_opt_h);

    SET_VECTOR_ELT(result, 4, convert_vector_vector_int_to_R(cpp_results.graph.adj_list)); n_protected++;
    SET_VECTOR_ELT(result, 5, convert_vector_vector_double_to_R(cpp_results.graph.weight_list)); n_protected++;
    SET_VECTOR_ELT(result, 6, convert_vector_double_to_R(cpp_results.predictions)); n_protected++;
    SET_VECTOR_ELT(result, 7, convert_vector_double_to_R(cpp_results.local_predictions)); n_protected++;

    if (cpp_results.bb_predictions.size() > 0) {
        SET_VECTOR_ELT(result, 8, convert_vector_double_to_R(cpp_results.bb_predictions)); n_protected++;
        SET_VECTOR_ELT(result, 9, convert_vector_double_to_R(cpp_results.ci_lower)); n_protected++;
        SET_VECTOR_ELT(result, 10,convert_vector_double_to_R(cpp_results.ci_upper)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    if (cpp_results.true_errors.size() > 0) {
        SEXP s_true_error = PROTECT(allocVector(REALSXP, 1)); n_protected++;
        double mean_true_error = std::accumulate(cpp_results.true_errors.begin(),
                                                 cpp_results.true_errors.end(), 0.0) /  cpp_results.true_errors.size();
        REAL(s_true_error)[0] = mean_true_error;
        SET_VECTOR_ELT(result, 11, s_true_error);
    } else {
        SET_VECTOR_ELT(result, 11, R_NilValue);
    }

    SET_VECTOR_ELT(result, 12,
                   convert_vector_vector_double_to_R(cpp_results.h_predictions)); n_protected++;


    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("h_values"));
    SET_STRING_ELT(names, 1, mkChar("h_errors"));
    SET_STRING_ELT(names, 2, mkChar("opt_h_idx"));
    SET_STRING_ELT(names, 3, mkChar("opt_h"));
    SET_STRING_ELT(names, 4, mkChar("graph_adj_list"));
    SET_STRING_ELT(names, 5, mkChar("graph_edge_lengths"));
    SET_STRING_ELT(names, 6, mkChar("predictions"));
    SET_STRING_ELT(names, 7, mkChar("local_predictions"));
    SET_STRING_ELT(names, 8, mkChar("bb_predictions"));
    SET_STRING_ELT(names, 9, mkChar("ci_lower"));
    SET_STRING_ELT(names, 10, mkChar("ci_upper"));
    SET_STRING_ELT(names, 11, mkChar("true_error"));
    SET_STRING_ELT(names, 12, mkChar("h_predictions"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}


/**
 * @brief Implements Path Graph Linear Model fitting with precomputed path models and proper model averaging
 *
 * @details This function implements an optimized version of the Path Graph Linear Model algorithm
 * that precomputes and caches linear models for all paths, then performs proper model averaging
 * over all valid models containing each vertex. The algorithm consists of two main phases:
 *
 * Phase 1 - Model Precomputation:
 * - For each path of length h+1:
 *   - Computes path coordinates based on edge weights
 *   - Calculates kernel weights for each position in the path
 *   - Fits linear models for all positions
 *   - Associates models with vertices they contain
 *   - Filters models based on position deviation from path midpoint
 *
 * Phase 2 - Model Averaging:
 * - For each vertex:
 *   - Collects all valid models containing the vertex
 *   - Performs weighted averaging using original kernel weights
 *   - Computes weighted average predictions and errors
 *
 * The function uses kernel weighting to give more weight to nearby points when fitting local
 * linear models. The kernel weights are normalized and combined with provided sample weights.
 *
 * @param path_graph The path graph structure containing:
 *    - h: Maximum hop distance
 *    - adj_list: Adjacency lists for h-hop neighborhoods
 *    - weight_list: Accumulated weights to h-hop neighbors
 *    - shortest_paths: Map of vertex pairs to shortest paths between them
 *    - vertex_paths: For each vertex, lists of paths containing it
 *    - longest_paths: Vector of path endpoints for paths of length h
 * @param y Vector of response values for each vertex
 * @param weights Vector of sample weights for each vertex (e.g., from Bayesian bootstrap)
 * @param ikernel Integer specifying the kernel type for distance weighting
 * @param max_distance_deviation Maximum allowed deviation from path midpoint when selecting valid models
 * @param dist_normalization_factor Factor for normalizing distances in kernel weight computation (default: 1.01)
 * @param epsilon Small constant for numerical stability in linear model fitting (default: 1e-8)
 *
 * @return std::pair containing:
 *    - first: Vector of model-averaged predictions for each vertex
 *    - second: Vector of model-averaged LOOCV errors for each vertex
 *
 * @throws Rf_error if no valid models are found for any vertex
 *
 * @note The function handles binary response variables (y ∈ {0,1}) by clamping predictions
 * to [0,1]. For regular regression, no clamping is performed.
 *
 * @pre
 * - path_graph must be properly initialized with valid paths
 * - y.size() == weights.size() == number of vertices in path_graph
 * - All vectors in path_graph must be properly sized
 * - ikernel must specify a valid kernel type
 * - max_distance_deviation must be non-negative
 * - dist_normalization_factor must be positive
 * - epsilon must be positive
 *
 * @warning
 * - The function modifies the kernel state using initialize_kernel()
 * - Large graphs with many paths may require significant memory
 * - Performance depends heavily on path structure and kernel type
 *
 * @see path_graph_plm_t
 * @see vertex_model_info_t
 * @see predict_lms_1d_loocv
 * @see initialize_kernel
 * @see kernel_fn
 *
 * Example usage:
 * @code
 * path_graph_plm_t graph;  // Initialize graph structure
 * std::vector<double> y = {1.0, 2.0, 3.0};  // Response values
 * std::vector<double> w = {1.0, 1.0, 1.0};  // Sample weights
 * int kernel_type = 1;  // Gaussian kernel
 * int max_dev = 2;     // Maximum deviation from midpoint
 *
 * auto [predictions, errors] = pgmalog(
 *     graph, y, w, kernel_type, max_dev);
 * @endcode
 */


// path graph model averaged local linear model
// * @return std::pair containing:
// *    - First element: Vector of fitted values (estimated conditional expectations) for each vertex
// *    - Second element: Vector of weighted LOOCV errors for each vertex
std::pair<std::vector<double>, std::vector<double>> pgmalog(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int kernel_type,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon,
    bool verbose) {

    int h = path_graph.h; // must be even and such that h >= 4 and h <= n_vertices - 2, where n_vertices is the number of vertices of the graph
    int mid_vertex = h / 2;
    int path_n_vertices = h + 1;
    int n_vertices = path_graph.vertex_paths.size();

    std::vector<double> predictions(n_vertices);
    std::vector<double> errors(n_vertices);

    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("PGMALOG");

    if (verbose) {
        Rprintf("Starting PGMALOG computation\n");
        Rprintf("Path graph size: %d vertices\n", n_vertices);
        Rprintf("Path graph h: %d\n", h);
    }

    initialize_kernel(kernel_type, 1.0);

    // Prepare data for linear model
    std::vector<double> y_path(path_n_vertices);
    std::vector<double> x_path(path_n_vertices);
    std::vector<double> d_path(path_n_vertices);
    std::vector<double> w_path(path_n_vertices);

    // Constructing a mapping assigning to each edge - represented by a pair
    // (start, end) with start < end - its weight
    std::map<std::pair<int,int>, double> edge_weight_map;
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        auto neighbors = path_graph.adj_list[vertex];
        auto weights = path_graph.weight_list[vertex];
        for (size_t j = 0; j < neighbors.size(); j++) {
            if (vertex < neighbors[j]) {
                edge_weight_map[std::make_pair(vertex, neighbors[j])] = weights[j];
            }
        }
    }

    // Phase 1:
    //
    // 1) For each reference vertex, ref_vertex, find all path_n_vertices paths
    // containing the given vertex.
    //
    // 2) For each path compute the distance from the reference vertex to the
    // mid point of the path.
    //
    // 3) Let min_dist be the minimum of these distances
    //
    // 4) Select all paths whose distance to the min point is no more than min_dist + max_distance_deviation
    //
    // 5) For each of these paths
    //
    // a) fit a weighted linear model to y restricted to the path with x being
    // the distance from the initial vertex of the path. The wieghts are
    // computed using the distance from the reference vertex to other vertices
    // of the path.
    //
    // b) store that model in the vector vertex_models[ref_vertex] contatining
    // all models containing ref_vertex in their support
    //

    //auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("  Phase 1: Computing single-model predictions ... ");
    }
    auto phase1_ptm = std::chrono::steady_clock::now();

    struct pred_w_err_t {
        double prediction;
        double weight;
        double error;
    } pred_w_err;

    std::vector<std::vector<pred_w_err_t>> vertex_pred_w_err(n_vertices);

    for (int ref_vertex = 0; ref_vertex < n_vertices; ++ref_vertex) {

        auto vertex_path_graph_info = path_graph.vertex_paths[ref_vertex]; // a struct with two components:
        // containing_paths: List of (start,end) pairs for paths that contain the given vertex
        // position_in_path: Position of vertex in each path

        // Get all paths containing the vertex together with the positions within the path of the given vertex 'ref_vertex'
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(path_n_vertices, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            Rf_error("No paths found for vertex %d", ref_vertex);
        }

        // Find paths within allowed deviation
        std::vector<size_t> valid_path_indices;
        int min_dist_to_mid_vertex = h;

        // First pass: find minimum distance
        for (const auto& path : vertex_paths) {
            int dist_to_mid_vertex = std::abs(path.second - mid_vertex);
            min_dist_to_mid_vertex = std::min(min_dist_to_mid_vertex, dist_to_mid_vertex);
        }

        // Second pass: collect valid paths
        for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            int dist_to_mid_vertex = std::abs(vertex_paths[path_i].second - mid_vertex);
            if (dist_to_mid_vertex <= min_dist_to_mid_vertex + max_distance_deviation) {
                valid_path_indices.push_back(path_i);
            }
        }

        if (valid_path_indices.empty()) {
            Rf_error("No valid paths found for vertex %d within deviation limits", ref_vertex);
        }

        // For each valid path fit a weighted linear model of y restricted to the path
        for (const auto& path_i : valid_path_indices) {

            std::vector<int> path = vertex_paths[path_i].first;
            int position_in_path = vertex_paths[path_i].second;

            // Clear data vectors
            y_path.clear();
            x_path.clear();
            d_path.clear();
            x_path[0] = 0;
            y_path[0] = y[path[0]];

            // Extract y values and compute distances
            for (int i = 1; i < path_n_vertices; ++i) {
                y_path[i] = y[path[i]];
                x_path[i] = x_path[i - 1] + edge_weight_map[{std::min(path[i - 1], path[i]), std::max(path[i - 1], path[i])}];
            }

            // Compute kernel weights
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[position_in_path]);
                max_dist = std::max(max_dist, d_path[i]);
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] /= max_dist;
            }

            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] = w_path[i] / total_w_path * weights[path[i]];

            bool fit_quadratic = false;
            int max_iterations = 100;
            double ridge_lambda = 1e-6;
            bool with_errors = true;
            eigen_ulogit_t fit = eigen_ulogit_fit(
                x_path.data(),
                y_path.data(),
                w_path,
                fit_quadratic,
                max_iterations,
                ridge_lambda,
                epsilon,
                with_errors
                );

            // For each vertex of the path record predicted value, the weight, and the models LOOCV at that vertex
            for (int i = 0; i < path_n_vertices; ++i) {
                pred_w_err.prediction = fit.predictions[i];
                pred_w_err.weight     = w_path[i];
                pred_w_err.error      = fit.loocv_brier_errors[i];
                vertex_pred_w_err[path[i]].push_back(pred_w_err);
            }
        } // END OF for (const auto& path_i : valid_path_indices)
    } // END OF for (int ref_vertex = 0; ref_vertex < n_vertices; ++ref_vertex)

    if (verbose) {
        elapsed_time(phase1_ptm, "Done");
        mem_tracker.report();
    }

    //
    // Phase 2: Model averaging
    //
    if (verbose) {
        Rprintf("  Phase 2: Computing model-averaged predictions ... ");
    }
    auto phase2_ptm = std::chrono::steady_clock::now();

    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    double wmean_error = 0.0;

    for (int i = 0; i < n_vertices; i++) {
        weighted_sum = 0.0;
        weight_sum = 0.0;
        wmean_error = 0.0;
        for (const auto& v : vertex_pred_w_err[i]) {
            weighted_sum += v.weight * v.prediction;
            weight_sum   += v.weight;
            wmean_error  += v.weight * v.error;
        }
        predictions[i] = weighted_sum / weight_sum;
        errors[i] = wmean_error / weight_sum;
    }

    if (verbose) {
        elapsed_time(phase2_ptm, "Done");
        mem_tracker.report();
    }

    if (verbose) {
        elapsed_time(total_ptm, "\nTotal PGMALOG computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return std::make_pair(predictions, errors);
}

