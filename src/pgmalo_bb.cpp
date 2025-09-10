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
#include "pgmalo.hpp"
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
/**
 * @brief Performs parallel bootstrap calculations for Graph Path Linear Model using STL parallel algorithms
 *
 * @param path_graph The path graph structure containing adjacency lists
 * @param y Vector of observed values
 * @param n_bb Number of bootstrap iterations
 * @param max_distance_deviation Maximum allowed deviation in distance calculations
 * @param ikernel Kernel function identifier
 * @param dist_normalization_factor Factor for normalizing distances
 * @param epsilon Numerical tolerance parameter
 *
 * @return Vector of vectors containing bootstrap predictions for each iteration
 */
std::vector<std::vector<double>> pgmalo_bb(const path_graph_plm_t& path_graph,
                                         const std::vector<double>& y,
                                         int n_bb,
                                         int max_distance_deviation,
                                         int ikernel,
                                         double dist_normalization_factor = 1.01,
                                         double epsilon = 1e-8) {

    int n_vertices = static_cast<int>(y.size());

    // Initialize results vector
    std::vector<std::vector<double>> bb_Ey(n_bb);
    for (auto& Ey : bb_Ey) {
        Ey.resize(n_vertices);
    }

    // Create indices for parallel iteration
    std::vector<int> bb_indices(n_bb);
    std::iota(bb_indices.begin(), bb_indices.end(), 0);

    // Mutex for thread-safe random number generation
    std::mutex rng_mutex;

    // Parallel execution of bootstrap iterations
    gflow::for_each(GFLOW_EXEC_POLICY,
                    bb_indices.begin(),
                    bb_indices.end(),
                    [&](int iboot) {
                        // Thread-local weight vector
                        std::vector<double> weights(n_vertices);

                        // Generate weights in a thread-safe manner
                        {
                            std::lock_guard<std::mutex> lock(rng_mutex);
                            C_runif_simplex(&n_vertices, weights.data());
                        }

                        // Compute predictions for this bootstrap iteration
                        auto [predictions, errors] = spgmalo(path_graph,
                                                             y,
                                                             weights,
                                                             ikernel,
                                                             max_distance_deviation,
                                                             dist_normalization_factor,
                                                             epsilon);

                        // Store results - no need for mutex as each thread writes to its own index
                        bb_Ey[iboot] = std::move(predictions);
                    });

    return bb_Ey;
}

/**
 * @brief Performs Bayesian bootstrap estimation with credible intervals for graph path linear models
 *
 * @details This function implements a complete Bayesian bootstrap analysis in three steps:
 * 1. Performs multiple iterations of weighted graph path linear model estimation
 * 2. Computes central location estimates (mean or median) for each vertex
 * 3. Calculates credible intervals at the specified probability level
 *
 * The implementation supports parallel processing for large graphs and provides
 * robust uncertainty quantification through bootstrap resampling.
 *
 * @param path_graph PLM graph structure containing adjacency lists and path information
 * @param y Response variable values at each vertex
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 100)
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param use_median If true, uses median instead of mean for central location (default: false)
 * @param ikernel Kernel function selector (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 *
 * @return bb_cri_t struct containing:
 *         - bb_Ey: Central location estimates (mean/median) for each vertex
 *         - cri_L: Lower bounds of credible intervals
 *         - cri_U: Upper bounds of credible intervals
 *
 * @throws Rf_error if:
 *         - Input dimensions are inconsistent
 *         - Parameters are invalid (p ∉ (0,1), n_bb ≤ 0, n_cores ≤ 0)
 *         - Numerical instability occurs during computation
 *         - Memory allocation fails
 */
bb_cri_t pgmalo_bb_cri(const path_graph_plm_t& path_graph,
                     const std::vector<double>& y,
                     double p,
                     int n_bb,
                     int max_distance_deviation,
                     bool use_median,
                     int ikernel,
                     double dist_normalization_factor,
                     double epsilon) {

    // Perform bootstrap iterations
    std::vector<std::vector<double>> bb_Eys = pgmalo_bb(path_graph,
                                                      y,
                                                      n_bb,
                                                      max_distance_deviation,
                                                      ikernel,
                                                      dist_normalization_factor,
                                                      epsilon);

    // Calculate credible intervals
    return bb_cri(bb_Eys, use_median, p);
}


#if 0

/**
 * @brief Computes double-weighted local linear models with Bayesian bootstrap weights
 *
 * @details This function extends pgmalo_without_bb_weights by incorporating Bayesian bootstrap weights
 * into the kernel-weighted local linear models. The weights modify the contribution of each
 * vertex in both the model fitting and prediction stages.
 *
 * Key modifications from pgmalo_without_bb_weights:
 * 1. Kernel weights are multiplied by corresponding Bayesian bootstrap weights
 * 2. LOOCV computations account for the modified weighting scheme
 * 3. Path-based model combinations incorporate both kernel and bootstrap weights
 *
 * For each vertex in the graph, the function:
 * 1. Identifies valid paths based on position and deviation constraints
 * 2. For each path:
 *    - Computes kernel weights based on distances
 *    - Modifies weights by multiplying with Bayesian bootstrap weights
 *    - Fits local linear models using the combined weights
 *    - Computes LOOCV errors
 * 3. Selects or combines models based on max_distance_deviation parameter
 *
 * @param path_graph PLM graph structure containing:
 *    - Adjacency lists
 *    - Edge weights
 *    - Path information
 *    - Path length h (must be positive and odd)
 * @param y Vector of response variables (one per vertex)
 * @param weights Bayesian bootstrap weights with properties:
 *    - Same length as y
 *    - Sum to 1.0
 *    - All components are strictly positive
 * @param ikernel Integer specifying the kernel type for distance-based weighting
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 *        - If 0: Uses only the best model for each vertex
 *        - If > 0: Combines predictions from multiple models using path-based weights
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return std::pair containing:
 *    - First element: Vector of fitted values (estimated conditional expectations) for each vertex
 *    - Second element: Vector of weighted LOOCV errors for each vertex
 *
 * @throws Rf_error if:
 *    - Path length h is not positive or not odd
 *    - Invalid kernel type is specified
 *    - Length of y or weights doesn't match number of vertices
 *    - Weights don't sum to 1 (within epsilon)
 *    - Any weight is zero or negative
 *    - Distance normalization factor is not positive
 *    - max_distance_deviation is negative
 *    - min_required_paths is less than 1
 *    - No valid paths are found for a vertex within deviation limits
 *    - Numerical instability is detected during model fitting
 *
 * @note
 * 1. The function assumes weights are properly normalized (sum to 1)
 * 2. The modification of kernel weights by bootstrap weights affects:
 *    - Model fitting through weighted least squares
 *    - LOOCV Rf_error computation
 *    - Model averaging when max_distance_deviation > 0
 * 3. The returned LOOCV errors incorporate both kernel and bootstrap weights
 *
 * @see predict_lm_1d_loocv for details on the weighted LOOCV implementation
 * @see path_graph_plm_t for the expected graph structure
 *
 * Example usage:
 * @code
 * path_graph_plm_t graph = // ... initialize graph
 * std::vector<double> y = // ... response variables
 * std::vector<double> weights = // ... Bayesian bootstrap weights
 * int kernel_type = 1;  // e.g., Gaussian kernel
 * int max_dev = 2;      // allow up to 2 positions deviation from center
 * int min_paths = 3;    // require at least 3 valid paths
 *
 * auto [predictions, errors] = pgmalo_without_bb_weights_with_bb_weights(
 *     graph, y, weights, kernel_type, max_dev, min_paths);
 * @endcode
 */
std::pair<std::vector<double>, std::vector<double>> pgmalo_with_global_bb_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    // Input validation
    if (max_distance_deviation < 0) {
        Rf_error("max_distance_deviation must be non-negative");
    }

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }
    if (weights.size() != n_vertices) {
        Rf_error("Length of weights must match number of vertices");
    }

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Validate weights
    double weight_sum = 0.0;
    for (const auto& w : weights) {
        #if 0
        if (w <= 0) {
            Rf_error("All weights must be strictly positive");
        }
        #endif
        weight_sum += w;
    }
    if (std::abs(weight_sum - 1.0) > epsilon) {
        Rf_error("Weights must sum to 1");
    }

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
            std::vector<double> bb_weights_path(path_n_vertices);  // Bayesian bootstrap weights for this path
            x_path[0] = 0;

            // Extract y values, weights, and compute distances
            for (int i = 0; i < path_n_vertices; ++i) {
                y_path[i] = y[path[i]];
                bb_weights_path[i] = weights[path[i]];  // Get BB weight for this vertex
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

            double scale = 1.0;
            initialize_kernel(ikernel, scale);
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] /= total_w_path;

            // Multiply kernel weights by Bayesian bootstrap weights
            for (int i = 0; i < path_n_vertices; ++i) {
                w_path[i] *= bb_weights_path[i];
            }

            // Fit LOOCV model with combined weights
            auto model = predict_lm_1d_loocv(y_path, x_path, w_path, path, position_in_path, y_binary, epsilon);

            // Update best model if this one has lower LOOCV Rf_error
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
* @brief Performs Bayesian bootstrap sampling for graph path linear models with path-specific weights
*
* This function implements a modified version of graph path linear models (PGMALO) with Bayesian
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
*    - Computes LOOCV Rf_error for model selection
* 3. Selects the best model based on LOOCV Rf_error
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
std::pair<std::vector<double>, std::vector<double>> pgmalo_with_bb_weights(
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

            #ifdef _OPENMP
#pragma omp critical(kernel_init)
            #endif
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

            // Update best model if this one has lower LOOCV Rf_error
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
 * @brief Computes conditional expectations using path-based weighted local linear models
 *
 * @details This function estimates conditional expectations at each vertex using weighted
 * local linear models computed along paths in the graph. The algorithm uses a two-level
 * weighting scheme:
 *
 * 1. Distance-based kernel weights for individual model fitting:
 *    - For each path, distances between vertices are computed using edge weights
 *    - Kernel weights are calculated based on normalized distances from the target vertex
 *
 * 2. Path-based weights for combining predictions (when max_distance_deviation > 0):
 *    - For each vertex i, identifies its best model (lowest LOOCV Rf_error)
 *    - The final estimate at vertex i is a weighted average:
 *      \[
 *      \hat{y}(i) = \frac{\sum_j w_j \hat{y}_j}{\sum_j w_j}
 *      \]
 *      where j ranges over vertices in i's best path, w_j is the kernel weight of
 *      vertex j in that path, and \hat{y}_j is vertex j's best model prediction at vertex i
 *
 * Algorithm steps for each vertex:
 * 1. Find all paths of length h containing the vertex
 * 2. Filter valid paths based on position:
 *    - Compute minimum deviation from center position across all paths
 *    - Select paths where vertex position deviates at most max_distance_deviation from the minimum deviation
 * 3. For each valid path:
 *    - Compute distances and kernel weights
 *    - Fit weighted linear model and compute LOOCV errors
 *    - Track the model with lowest LOOCV Rf_error at the vertex
 * 4. Compute final predictions:
 *    - If max_distance_deviation = 0: Use best model's prediction directly
 *    - If max_distance_deviation > 0: Use weighted average of predictions from best models
 *
 * @param path_graph Graph structure containing adjacency lists, edge weights, and path info
 * @param y Response variables (one per vertex)
 * @param ikernel Kernel type for distance-based weighting
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 *        - 0: Use only best model predictions
 *        - >0: Use weighted average of predictions from multiple models
 * @param dist_normalization_factor Factor for scaling distances (default: 1.01)
 * @param epsilon Numerical stability threshold (default: 1e-8)
 *
 * @return std::pair containing:
 *    - Vector of estimated conditional expectations for each vertex
 *    - Vector of LOOCV errors for each vertex's best model
 *
 * @throws Rf_error if:
 *    - No paths found for a vertex
 *    - No valid paths within deviation limits
 *
 * @note
 * - Assumes graph is valid and connected
 * - Path length h must be positive and odd
 * - Input validation performed by calling R function
 */
std::pair<std::vector<double>, std::vector<double>> pgmalo_without_bb_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    double scale = 1.0;
    initialize_kernel(ikernel, scale); // the function sets up the global kernel function pointer based on the specified kernel type

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    #define DEBUG_pgmalo_without_bb_weights 0
    #if DEBUG_pgmalo_without_bb_weights
    Rprintf("\nIn pgmalo_without_bb_weights()\n");
    Rprintf("Number of vertices: %d\n", n_vertices);
    Rprintf("h: %d\n", h);
    //print_vect_vect(neighbors, "neighbors");
    //print_vect_vect(edge_lengths, "edge_lengths");
    #endif

    std::vector<double> Ey(n_vertices);
    std::vector<double> errors(n_vertices);
    std::vector<lm_loocv_t> best_models(n_vertices); // stores the best model for each vertex

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {

        #if DEBUG_pgmalo_without_bb_weights
        Rprintf("\nvertex_i: %d\n", vertex_i);
        #endif

        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i]; // a struct with two components:
        // containing_paths: List of (start,end) pairs for paths that contain the given vertex
        // position_in_path: Position of vertex in each path

        // Get all paths containing the vertex together with the positions within the path of the given vertex 'vertex_i'
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

        // For each valid path, compute the LOOCV Rf_error
        double best_loocv_error = std::numeric_limits<double>::infinity();
        lm_loocv_t best_model;

        for (const auto& path_i : valid_path_indices) {

            #if DEBUG_pgmalo_without_bb_weights
            Rprintf("path_i: %d\n", (int)path_i);
            #endif

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
                if (i > 0) { // if we had an (edge -> edge weight) map - this part could have been much faster
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

            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] /= total_w_path;

            //#if DEBUG_pgmalo_without_bb_weights
            #if 0
            print_vect(path, "path");
            Rprintf("position_in_path: %d\n", position_in_path);
            print_vect(y_path, "y_path");
            print_vect(x_path, "x_path");
            print_vect(w_path, "w_path");
            #endif

            // Fit a weighted linear modela and return as well model's LOOCV errors
            auto model = predict_lm_1d_loocv(y_path, x_path, w_path, path, position_in_path, y_binary, epsilon);

            // Update best model if this one has lower LOOCV Rf_error
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

        #if DEBUG_pgmalo_without_bb_weights
        Rprintf("\nIn the if (max_distance_deviation > 0) block\n");
        #endif

        for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
            lm_loocv_t best_model = best_models[vertex_i];
            // Get vertices of the best model path
            std::vector<int> path_vertex_indices = best_model.vertex_indices;
            std::vector<double> path_vertex_weights = best_model.w_values;
            int path_length = path_vertex_weights.size();

            #if DEBUG_pgmalo_without_bb_weights
            Rprintf("\nvertex_i: %d\n", vertex_i);
            #endif

            // for each vertex from path_vertex_indices use its best model to predict the value of that model at 'vertex_i'
            // Of course for 'vertex_i' we have already computed the predicted value of its best model: best_model.predicted_value = Ey[vertex_i]
            double weight = path_vertex_weights[best_model.ref_index];
            double total_weight = weight;
            double vertex_Ey = weight * Ey[vertex_i];

            for (int i = 0; i < path_length; ++i) {
                int index =  path_vertex_indices[i];
                if (index != vertex_i) {
                    lm_loocv_t i_best_model = best_models[index];

                    //#if DEBUG_pgmalo_without_bb_weights
                    #if 0
                    Rprintf("i: %d\tindex: %d\tvertex_i: %d\n", i, index, vertex_i);
                    print_vect(i_best_model.vertex_indices, "i_best_model.vertex_indices");
                    #endif

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

#endif
