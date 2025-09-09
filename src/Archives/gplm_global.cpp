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

#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <omp.h>
#include <unordered_map>

#include "msr2.h"
#include "path_graphs.h"
#include "cpp_utils.h"                 // for elapsed_time
#include "SEXP_cpp_conversion_utils.h"
#include "kernels.h"
#include "1D_linear_models.h"
#include "predictive_errors.h"
#include "adaptive_nbhd_size.h"

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
                               bool y_binary = false,
                               double epsilon = 1e-8);

extern "C" {
    SEXP S_univariate_gplm_global(SEXP s_x,
                                  SEXP s_y,
                                  SEXP s_y_true,
                                  SEXP s_max_distance_deviation,
                                  SEXP s_h_min,
                                  SEXP s_h_max,
                                  SEXP s_p,
                                  SEXP s_n_bb,
                                  SEXP s_ikernel,
                                  SEXP s_n_cores,
                                  SEXP s_dist_normalization_factor,
                                  SEXP s_epsilon,
                                  SEXP s_seed,
                                  SEXP s_verbose);
}

/**
 * @brief Performs Bayesian bootstrap estimation of generalized path-based local models (GPLM)
 *
 * This function computes multiple bootstrap replicates of GPLM estimates by:
 * 1. Generating Bayesian bootstrap weights for each vertex
 * 2. Finding optimal paths for each vertex
 * 3. Computing kernel-weighted local linear models along these paths
 * 4. Aggregating predictions when multiple valid paths exist
 *
 * The function supports parallel execution through OpenMP when n_cores > 1.
 *
 * @param path_graph Graph structure containing:
 *    - Vertex paths
 *    - Adjacency lists
 *    - Edge weights
 *    - Path length parameter h
 * @param y Vector of response values for each vertex
 * @param n_bb Number of Bayesian bootstrap replicates to compute
 * @param n_cores Number of cores to use for parallel computation:
 *    - If n_cores > 1: Uses OpenMP parallel processing
 *    - If n_cores = 1: Runs in serial mode
 * @param ikernel Integer specifying the kernel type for local weighting
 * @param max_distance_deviation Maximum allowed deviation from the optimal path distance:
 *    - Controls path selection flexibility
 *    - Must be non-negative
 *    - If > 0, enables model averaging across multiple paths
 * @param dist_normalization_factor Factor for normalizing path distances (default: 1.01)
 * @param epsilon Numerical tolerance for various computations (default: 1e-8):
 *    - Used in weight validation
 *    - Used in linear model fitting
 *
 * @return std::pair containing:
 *    - First: Vector of bootstrap estimates (n_bb × n_vertices matrix)
 *         Each row contains GPLM predictions for all vertices for one bootstrap replicate
 *    - Second: Vector of bootstrap errors (n_bb × n_vertices matrix)
 *         Each row contains leave-one-out cross-validation errors for all vertices
 *
 * @throws Rf_error if:
 *    - max_distance_deviation is negative
 *    - y.size() doesn't match number of vertices
 *    - No paths found for a vertex
 *    - No valid paths within deviation limits
 *
 * @note Performance Considerations:
 *    - Pre-allocates memory for all major data structures
 *    - Pre-computes path information before bootstrap loop
 *    - Uses thread-local storage for parallel execution
 *    - Employs move semantics for efficient return of large data structures
 *
 * @see gplm_with_global_bb_weights For the core computation on a single set of weights
 * @see path_graph_plm_t For the graph structure documentation
 *
 * @example
 * ```cpp
 * path_graph_plm_t graph = // ... initialize graph
 * std::vector<double> y = // ... response values
 * int n_bb = 500;        // number of bootstrap replicates
 * int n_cores = 4;       // number of parallel threads
 * auto result = gplm_global_bb(graph, y, n_bb, n_cores, 1, 2);
 * // result.first contains bootstrap estimates
 * // result.second contains bootstrap errors
 * ```
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> gplm_global_bb(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_bb,
    int n_cores,
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
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Pre-allocate all results
    std::vector<std::vector<double>> all_Ey(n_bb, std::vector<double>(n_vertices));
    std::vector<std::vector<double>> all_errors(n_bb, std::vector<double>(n_vertices));
    std::vector<std::vector<lm_loocv_t>> all_best_models(n_bb, std::vector<lm_loocv_t>(n_vertices));

    // Pre-allocate thread-local vectors
    std::vector<std::vector<double>> thread_weights(n_cores, std::vector<double>(n_vertices));
    std::vector<std::vector<double>> thread_w_path(n_cores, std::vector<double>(h));
    std::vector<std::vector<size_t>> thread_valid_path_indices(n_cores);
    for(auto& indices : thread_valid_path_indices) {
        indices.reserve(h);  // Reserve space for typical case
    }

    // Pre-compute path information for each vertex
    struct vertex_path_info_t {
        std::vector<std::pair<std::vector<int>, int>> paths;
        std::vector<size_t> valid_indices;
        int min_dist_to_mid_pt;
    };

    std::vector<vertex_path_info_t> vertex_path_info(n_vertices);
    int mid_pt = (h - 1) / 2;

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto& info = vertex_path_info[vertex_i];
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get all paths containing the vertex
        info.paths = vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);
        if (info.paths.empty()) {
            Rf_error("No paths found for vertex %d", vertex_i);
        }

        // Find minimum distance
        info.min_dist_to_mid_pt = h;
        for (const auto& path : info.paths) {
            int dist_to_mid_pt = std::abs(path.second - mid_pt);
            info.min_dist_to_mid_pt = std::min(info.min_dist_to_mid_pt, dist_to_mid_pt);
        }

        // Find valid paths
        for (size_t path_i = 0; path_i < info.paths.size(); ++path_i) {
            int dist_to_mid_pt = std::abs(info.paths[path_i].second - mid_pt);
            if (dist_to_mid_pt <= info.min_dist_to_mid_pt + max_distance_deviation) {
                info.valid_indices.push_back(path_i);
            }
        }

        if (info.valid_indices.empty()) {
            Rf_error("No valid paths found for vertex %d within deviation limits", vertex_i);
        }
    }

    // Pre-allocate vectors for linear models (per bootstrap replicate)
    std::vector<std::vector<double>> all_y_path(n_bb, std::vector<double>(h));
    std::vector<std::vector<double>> all_x_path(n_bb, std::vector<double>(h));
    std::vector<std::vector<double>> all_d_path(n_bb, std::vector<double>(h));
    std::vector<std::vector<double>> all_bb_weights_path(n_bb, std::vector<double>(h));

    auto process_bootstrap = [&](int iboot, int thread_id) {
        auto& weights = thread_weights[thread_id];
        auto& w_path = thread_w_path[thread_id];
        auto& y_path = all_y_path[iboot];
        auto& x_path = all_x_path[iboot];
        auto& d_path = all_d_path[iboot];
        auto& bb_weights_path = all_bb_weights_path[iboot];
        auto& best_models = all_best_models[iboot];

        // Generate bootstrap weights
        C_runif_simplex(&n_vertices, weights.data());

        for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
            const auto& vertex_info = vertex_path_info[vertex_i];

            // For each valid path, compute the LOOCV model
            double best_loocv_error = std::numeric_limits<double>::infinity();
            lm_loocv_t best_model;

            for (const auto& path_i : vertex_info.valid_indices) {
                const auto& path_pair = vertex_info.paths[path_i];
                const auto& path = path_pair.first;
                int position_in_path = path_pair.second;
                int path_n_vertices = path.size();

                // Prepare data for linear model
                x_path[0] = 0;

                // Extract y values, weights, and compute distances
                for (int i = 0; i < path_n_vertices; ++i) {
                    y_path[i] = y[path[i]];
                    bb_weights_path[i] = weights[path[i]];
                    if (i > 0) {
                        auto neighbors = path_graph.adj_list[path[i - 1]];
                        auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                        auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                        int neighbor_i = it - neighbors.begin();
                        x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i];
                    }
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

                double scale = 1.0;
                initialize_kernel(ikernel, scale);
                kernel_fn(d_path.data(), path_n_vertices, w_path.data());

                double total_w_path = std::accumulate(w_path.begin(),
                                                    w_path.begin() + path_n_vertices, 0.0);
                for (int i = 0; i < path_n_vertices; ++i)
                    w_path[i] /= total_w_path;

                // Multiply kernel weights by Bayesian bootstrap weights
                for (int i = 0; i < path_n_vertices; ++i) {
                    w_path[i] *= bb_weights_path[i];
                }

                // Fit LOOCV model with combined weights
                auto model = predict_lm_1d_loocv(y_path, x_path, w_path, path,
                                                 position_in_path, y_binary, epsilon);

                // Update best model if this one has lower LOOCV error
                if (model.loocv_at_ref_vertex < best_loocv_error) {
                    best_loocv_error = model.loocv_at_ref_vertex;
                    best_model = model;
                }
            }

            best_models[vertex_i] = best_model;
            all_Ey[iboot][vertex_i] = best_model.predicted_value;
            all_errors[iboot][vertex_i] = best_model.loocv_at_ref_vertex;
        }

        // Handle prediction based on max_distance_deviation
        if (max_distance_deviation > 0) {
            for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
                const auto& best_model = best_models[vertex_i];
                const auto& path_vertex_indices = best_model.vertex_indices;
                const auto& path_vertex_weights = best_model.w_values;
                int path_length = path_vertex_weights.size();

                double weight = path_vertex_weights[best_model.ref_index];
                double total_weight = weight;
                double vertex_Ey = weight * all_Ey[iboot][vertex_i];

                for (int i = 0; i < path_length; ++i) {
                    int index = path_vertex_indices[i];
                    if (index != vertex_i) {
                        const auto& i_best_model = best_models[index];
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
                all_Ey[iboot][vertex_i] = vertex_Ey;
            }
        }
    };

    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
        #pragma omp parallel for schedule(dynamic)
        for (int iboot = 0; iboot < n_bb; ++iboot) {
            int thread_id = omp_get_thread_num();
            process_bootstrap(iboot, thread_id);
        }
    } else {
        for (int iboot = 0; iboot < n_bb; ++iboot) {
            process_bootstrap(iboot, 0);
        }
    }

    return std::make_pair(std::move(all_Ey), std::move(all_errors));
}



/**
 * @brief Computes median predictions across bootstrap samples for each vertex
 *
 * @param predictions Matrix of predictions [bootstrap_idx][vertex_idx]
 *
 * @return Vector of median predictions for each vertex [vertex_idx]
 *
 * @note Uses nth_element for efficient median computation without full sorting
 */
std::vector<double> compute_median_predictions(const std::vector<std::vector<double>>& predictions) {
    std::vector<double> medians(predictions[0].size());
    std::vector<double> values;
    values.reserve(predictions.size());

    for (size_t v = 0; v < predictions[0].size(); ++v) {
        values.clear();
        for (const auto& pred : predictions) {
            values.push_back(pred[v]);
        }
        std::nth_element(values.begin(), values.begin() + values.size()/2, values.end());
        medians[v] = values[values.size()/2];
    }
    return medians;
}


/**
 * @brief Computes credible intervals from vertex-organized predictions
 *
 * @param vertex_predictions Matrix of predictions organized by vertex [vertex_idx][bootstrap_idx]
 * @param p Confidence level (e.g., 0.95 for 95% credible intervals)
 *
 * @return std::pair containing:
 *         - first: Vector of lower bounds [vertex_idx]
 *         - second: Vector of upper bounds [vertex_idx]
 *
 * @note Sorts predictions for each vertex to compute quantile-based intervals
 */
std::pair<std::vector<double>, std::vector<double>> compute_credible_intervals_from_predictions(
    const std::vector<std::vector<double>>& vertex_predictions,
    double p) {

    std::vector<double> lower(vertex_predictions.size());
    std::vector<double> upper(vertex_predictions.size());
    std::vector<double> sorted_values;

    for (size_t v = 0; v < vertex_predictions.size(); ++v) {
        sorted_values = vertex_predictions[v];
        std::sort(sorted_values.begin(), sorted_values.end());

        int lower_idx = static_cast<int>((1 - p) / 2 * sorted_values.size());
        int upper_idx = static_cast<int>((1 + p) / 2 * sorted_values.size());

        lower[v] = sorted_values[lower_idx];
        upper[v] = sorted_values[upper_idx];
    }

    return {lower, upper};
}




/**
 * @brief Compute median predictions and corresponding errors from bootstrap results
 */
std::pair<std::vector<double>, std::vector<double>> compute_median_predictions(
    const std::vector<std::vector<double>>& bb_predictions,
    const std::vector<std::vector<double>>& bb_errors) {

    int n_vertices = bb_predictions[0].size();
    int n_bb = bb_predictions.size();

    std::vector<double> median_predictions(n_vertices);
    std::vector<double> corresponding_errors(n_vertices);

    for (int v = 0; v < n_vertices; ++v) {
        // Create pairs of predictions and errors for this vertex
        std::vector<std::pair<double, double>> pred_error_pairs(n_bb);
        for (int b = 0; b < n_bb; ++b) {
            pred_error_pairs[b] = {bb_predictions[b][v], bb_errors[b][v]};
        }

        // Sort by predictions
        std::sort(pred_error_pairs.begin(), pred_error_pairs.end());

        // Get median prediction and its corresponding error
        int median_idx = n_bb / 2;
        median_predictions[v] = pred_error_pairs[median_idx].first;
        corresponding_errors[v] = pred_error_pairs[median_idx].second;
    }

    return {median_predictions, corresponding_errors};
}


/**
 * @brief Compute bootstrap credible intervals
 */
std::pair<std::vector<double>, std::vector<double>> compute_credible_intervals(
    const std::vector<std::vector<double>>& bb_predictions,
    double p) {

    int n_vertices = bb_predictions[0].size();
    int n_bb = bb_predictions.size();

    std::vector<double> lower(n_vertices);
    std::vector<double> upper(n_vertices);

    int lower_idx = static_cast<int>((1 - p) / 2 * n_bb);
    int upper_idx = static_cast<int>((1 + p) / 2 * n_bb);

    for (int v = 0; v < n_vertices; ++v) {
        std::vector<double> vertex_predictions(n_bb);
        for (int b = 0; b < n_bb; ++b) {
            vertex_predictions[b] = bb_predictions[b][v];
        }

        std::sort(vertex_predictions.begin(), vertex_predictions.end());

        lower[v] = vertex_predictions[lower_idx];
        upper[v] = vertex_predictions[upper_idx];
    }

    return {lower, upper};
}

/**
 * @brief Performs global path-based local mean analysis with adaptive neighborhood size selection
 *
 * This function analyzes a graph using different neighborhood sizes (h values) to find the optimal
 * neighborhood size for prediction. It uses bootstrap resampling to assess prediction accuracy
 * and compute credible intervals.
 *
 * @param neighbors Vector of neighbor indices for each vertex [vertex_idx][neighbor_idx]
 * @param edge_lengths Vector of edge lengths corresponding to neighbors [vertex_idx][neighbor_idx]
 * @param y Vector of observed values at each vertex [vertex_idx]
 * @param y_true Vector of true values for computing prediction error (optional) [vertex_idx]
 * @param max_distance_deviation Maximum allowed deviation in path distances for local computation
 * @param use_median Whether to use median instead of mean for local estimation (default: false)
 * @param h_min Minimum neighborhood size to test (must be odd, default: 3)
 * @param h_max Maximum neighborhood size to test (must be odd, default: 31)
 * @param p Confidence level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap samples (default: 500)
 * @param ikernel Kernel function selection (default: 1)
 * @param n_cores Number of CPU cores to use for parallel computation (default: 1)
 * @param dist_normalization_factor Factor for distance normalization (default: 1.01)
 * @param epsilon Small value to prevent division by zero (default: 1e-15)
 * @param seed Random seed for reproducibility (default: 0)
 * @param verbose Whether to print progress information (default: true)
 *
 * @return adaptive_nbhd_size_global_t containing:
 *         - Optimal neighborhood size and corresponding path graph
 *         - Per-h results including median predictions and CV errors
 *         - Bootstrap-based statistics for the optimal h value
 *         - Credible intervals for predictions
 *
 * @throws Rf_error if:
 *         - h_min > h_max
 *         - h_min or h_max is not odd
 *         - Input dimensions are inconsistent
 *         - p is not between 0 and 1
 *
 * @note The function organizes bootstrap results by vertex first, then bootstrap sample
 *       for efficient per-vertex analysis and credible interval computation.
 */
// new
adaptive_nbhd_size_global_t gplm_global(
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    const std::vector<double>& y_true,
    int max_distance_deviation,
    int h_min = 3,
    int h_max = 31,
    double p = 0.95,
    int n_bb = 500,
    int ikernel = 1,
    int n_cores = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15,
    unsigned int seed = 0,
    bool verbose = true) {

    auto total_ptm = std::chrono::steady_clock::now();
    auto ptm = std::chrono::steady_clock::now();  // Start timing

    seed = seed + 0;

    // Input validation
    if (h_min > h_max) {
        Rf_error("h_min must be less than or equal to h_max");
    }
    if (h_min % 2 == 0 || h_max % 2 == 0) {
        Rf_error("h_min and h_max must be odd numbers");
    }
    if (neighbors.size() != edge_lengths.size() || neighbors.size() != y.size()) {
        Rf_error("Inconsistent input dimensions");
    }
    if (p <= 0 || p >= 1) {
        Rf_error("p must be between 0 and 1");
    }

    int n_vertices = static_cast<int>(y.size());
    adaptive_nbhd_size_global_t results;

    int n_h_values = (h_max - h_min) / 2 + 1;
    std::vector<path_graph_plm_t> graphs(n_h_values);
    results.h_cv_errors.resize(n_h_values);
    results.h_values.resize(n_h_values);
    results.h_median_predictions.resize(n_h_values);
    results.h_vertex_cv_errors.resize(n_h_values);

    for (int i = 0, h = h_min; h <= h_max; h += 2, i++) {
        results.h_values[i] = h;

        if (verbose) {
            Rprintf("\nProcessing h = %d\n", h);
            Rprintf("\tCreating path graph ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto path_graph = create_path_graph_plm(neighbors, edge_lengths, h);
        if (verbose) elapsed_time(ptm, "DONE");

        if (verbose) {
            Rprintf("\tComputing bootstrap predictions and errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto [bb_predictions, bb_errors] = gplm_global_bb(path_graph,
                                                        y,
                                                        n_bb,
                                                        n_cores,
                                                        ikernel,
                                                        max_distance_deviation,
                                                        dist_normalization_factor,
                                                        epsilon);

        if (verbose) elapsed_time(ptm, "DONE");

        if (verbose) {
            Rprintf("\tReorganizing predictions and errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        // Compute median predictions for this h
        results.h_median_predictions[i] = compute_median_predictions(bb_predictions);

        // Reorganize CV errors by vertex
        results.h_vertex_cv_errors[i].resize(n_vertices);
        for (int v = 0; v < n_vertices; ++v) {
            results.h_vertex_cv_errors[i][v].resize(n_bb);
            for (int b = 0; b < n_bb; ++b) {
                results.h_vertex_cv_errors[i][v][b] = bb_errors[b][v];
            }
        }

        // Calculate mean error across vertices for this h
        results.h_cv_errors[i] = std::accumulate(bb_errors.begin(),
                                               bb_errors.end(),
                                               0.0,
                                               [](double sum, const std::vector<double>& errors) {
                                                   return sum + std::accumulate(errors.begin(),
                                                                             errors.end(),
                                                                             0.0);
                                               }) / (n_vertices * n_bb);

        graphs[i] = std::move(path_graph);

        if (verbose) elapsed_time(ptm, "DONE");
    }

    // Find optimal h
    auto min_it = std::min_element(results.h_cv_errors.begin(), results.h_cv_errors.end());
    int opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * opt_h_idx;

    // Store optimal graph
    results.opt_graph = std::move(graphs[opt_h_idx]);

    // Store optimal results
    results.opt_predictions = results.h_median_predictions[opt_h_idx];

    // Compute mean CV errors per vertex for optimal h
    results.opt_cv_errors.resize(n_vertices);
    for (int v = 0; v < n_vertices; ++v) {
        results.opt_cv_errors[v] = std::accumulate(results.h_vertex_cv_errors[opt_h_idx][v].begin(),
                                                 results.h_vertex_cv_errors[opt_h_idx][v].end(),
                                                 0.0) / n_bb;
    }

    // Compute credible intervals for optimal h
    auto [ci_lower, ci_upper] = compute_credible_intervals_from_predictions(
        results.h_vertex_cv_errors[opt_h_idx], p);
    results.opt_ci_lower = std::move(ci_lower);
    results.opt_ci_upper = std::move(ci_upper);

    // Computing true errors if available
    if (!y_true.empty() && y_true.size() == n_vertices) {
        if (verbose) {
            Rprintf("Computing true errors ... ");
            ptm = std::chrono::steady_clock::now();
        }
        results.opt_true_errors.resize(n_vertices);
        for (size_t i = 0; i < n_vertices; i++) {
            results.opt_true_errors[i] = std::abs(y_true[i] - results.opt_predictions[i]);
        }
        if (verbose) elapsed_time(ptm, "DONE");
    } else {
        results.opt_true_errors.clear();
    }

    if (verbose) elapsed_time(total_ptm, "Total elapsed time:");

    return results;
}

// old gplm_global
#if 0
/**
 * @brief Performs global adaptive neighborhood selection for path-based local models with bootstrap
 *
 * This function implements a comprehensive analysis pipeline that:
 * 1. Tests multiple neighborhood sizes (h values)
 * 2. Performs bootstrap replication for uncertainty quantification
 * 3. Computes cross-validation errors for model selection
 * 4. Determines the optimal neighborhood size
 * 5. Calculates credible intervals for predictions
 *
 * The algorithm uses path-based local models, where predictions at each vertex
 * are computed using weighted local linear regression along paths through the graph.
 * Bootstrap replication provides uncertainty quantification and robust model selection.
 *
 * @param neighbors Vector of adjacency lists representing the graph structure
 *                 Each inner vector contains indices of neighboring vertices
 *
 * @param edge_lengths Vector of edge weights corresponding to neighbors
 *                    Must match the structure of neighbors vector
 *
 * @param y Response values at each vertex
 *
 * @param y_true True response values for error calculation (optional)
 *               If empty, true errors won't be computed
 *
 * @param max_distance_deviation Maximum allowed deviation from optimal path distance
 *                             Controls path selection flexibility
 *
 * @param use_median Whether to use median instead of mean for aggregating predictions
 *                   Default: false
 *
 * @param h_min Minimum neighborhood size to test
 *              Must be odd and ≥ 3
 *              Default: 3
 *
 * @param h_max Maximum neighborhood size to test
 *              Must be odd and > h_min
 *              Default: 31
 *
 * @param p Credible interval probability level
 *          Must be between 0 and 1
 *          Default: 0.95 (95% credible intervals)
 *
 * @param n_bb Number of bootstrap replicates
 *             Default: 500
 *
 * @param ikernel Kernel function type for local weighting
 *                Default: 1
 *
 * @param n_cores Number of cores for parallel processing
 *                If > 1, enables OpenMP parallelization
 *                Default: 1
 *
 * @param dist_normalization_factor Factor for normalizing path distances
 *                                 Default: 1.01
 *
 * @param epsilon Numerical tolerance for computations
 *                Default: 1e-15
 *
 * @param seed Random seed for reproducibility
 *             Default: 0
 *
 * @param verbose Whether to print progress information
 *                Default: true
 *
 * @return adaptive_nbhd_size_global_t structure containing:
 *    - opt_graph: Path graph for the optimal h value
 *    - h_values: Vector of tested neighborhood sizes
 *    - opt_h: Optimal neighborhood size based on CV errors
 *    - h_cv_errors: Vector of mean cross-validation errors for each h
 *    - opt_predictions: Vector of median bootstrap predictions for optimal h
 *    - opt_cv_errors: Vector of per-vertex CV errors for optimal h
 *    - opt_true_errors: Vector of true prediction errors for optimal h (if y_true provided)
 *    - opt_ci_lower: Vector of lower bounds for credible intervals
 *    - opt_ci_upper: Vector of upper bounds for credible intervals
 *
 * @throws Rf_error if:
 *    - h_min > h_max
 *    - h_min or h_max is even
 *    - Input dimensions are inconsistent
 *    - p is not in (0,1)
 *    - No valid paths found for any vertex
 *
 * @note Performance Considerations:
 *    - Memory usage scales with number of vertices, bootstrap replicates,
 *      and number of h values tested
 *    - Parallel processing enabled when n_cores > 1
 *    - Progress information provided when verbose = true
 *
 * @see compute_median_predictions For bootstrap aggregation details
 * @see compute_credible_intervals For credible interval computation
 * @see path_graph_plm_t For underlying graph structure
 *
 * @example
 * ```cpp
 * // Create graph structure
 * std::vector<std::vector<int>> neighbors = {{1,2}, {0,2}, {0,1}};
 * std::vector<std::vector<double>> edge_lengths = {{1.0,1.0}, {1.0,1.0}, {1.0,1.0}};
 * std::vector<double> y = {1.0, 2.0, 3.0};
 *
 * // Run analysis
 * auto results = gplm_global(neighbors, edge_lengths, y, {}, 2);
 *
 * // Access results
 * int best_h = results.opt_h;
 * auto predictions = results.opt_predictions;
 * auto lower_ci = results.opt_ci_lower;
 * auto upper_ci = results.opt_ci_upper;
 * ```
 */
adaptive_nbhd_size_global_t gplm_global(
    const std::vector<std::vector<int>>& neighbors,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& y,
    const std::vector<double>& y_true,
    int max_distance_deviation,
    bool use_median = false,
    int h_min = 3,
    int h_max = 31,
    double p = 0.95,
    int n_bb = 500,
    int ikernel = 1,
    int n_cores = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15,
    unsigned int seed = 0,
    bool verbose = true) {

    double total_ptm = std::chrono::steady_clock::now();

    // Input validation
    if (h_min > h_max) {
        Rf_error("h_min must be less than or equal to h_max");
    }
    if (h_min % 2 == 0 || h_max % 2 == 0) {
        Rf_error("h_min and h_max must be odd numbers");
    }
    if (neighbors.size() != edge_lengths.size() || neighbors.size() != y.size()) {
        Rf_error("Inconsistent input dimensions");
    }
    if (p <= 0 || p >= 1) {
        Rf_error("p must be between 0 and 1");
    }

    int n_vertices = static_cast<int>(y.size());
    adaptive_nbhd_size_global_t results;

    int n_h_values = (h_max - h_min) / 2 + 1;
    std::vector<path_graph_plm_t> graphs(n_h_values);
    results.h_cv_errors.resize(n_h_values);    // Changed from cv_errors
    results.h_values.resize(n_h_values);

    std::vector<std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>> bootstrap_results(n_h_values);

    double ptm;
    for (int i = 0, h = h_min; h <= h_max; h += 2, i++) {
        results.h_values[i] = h;

        if (verbose) {
            Rprintf("\nProcessing h = %d\n", h);
            Rprintf("\tCreating path graph ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto path_graph = create_path_graph_plm(neighbors, edge_lengths, h);
        if (verbose) elapsed_time(ptm, "DONE");

        if (verbose) {
            Rprintf("\tComputing bootstrap predictions and errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        bootstrap_results[i] = gplm_global_bb(path_graph,
                                            y,
                                            n_bb,
                                            n_cores,
                                            ikernel,
                                            max_distance_deviation,
                                            dist_normalization_factor,
                                            epsilon);

        if (verbose) elapsed_time(ptm, "DONE");

        if (verbose) {
            Rprintf("\tComputing median predictions and errors ... ");
            ptm = std::chrono::steady_clock::now();
        }

        auto [median_preds, median_errors] = compute_median_predictions(
            bootstrap_results[i].first,
            bootstrap_results[i].second
        );

        // Calculate mean error across vertices
        results.h_cv_errors[i] = std::accumulate(median_errors.begin(),
                                               median_errors.end(), 0.0) / n_vertices;

        graphs[i] = std::move(path_graph);

        if (verbose) elapsed_time(ptm, "DONE");
    }

    // Find optimal h
    auto min_it = std::min_element(results.h_cv_errors.begin(), results.h_cv_errors.end());
    int opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * opt_h_idx;

    // Store optimal graph and results
    results.opt_graph = graphs[opt_h_idx];     // Changed from opt_h_graph

    // Store bootstrap results for optimal h
    results.bb_predictions = std::move(bootstrap_results[opt_h_idx].first);
    results.bb_errors = std::move(bootstrap_results[opt_h_idx].second);

    // Compute and store median predictions and errors for optimal h
    auto [opt_preds, opt_errs] = compute_median_predictions(
        results.bb_predictions,
        results.bb_errors
    );
    results.opt_predictions = std::move(opt_preds);    // Changed from predictions
    results.opt_cv_errors = std::move(opt_errs);      // Changed from vertex_cv_errors

    // Compute credible intervals
    auto [ci_lower, ci_upper] = compute_credible_intervals(results.bb_predictions, p);
    results.opt_ci_lower = std::move(ci_lower);    // Changed from opt_ci_lower
    results.opt_ci_upper = std::move(ci_upper);    // Changed from opt_ci_upper

    // Computing true errors if available
    if (!y_true.empty() && y_true.size() == n_vertices) {
        if (verbose) {
            Rprintf("Computing true errors ... ");
            ptm = std::chrono::steady_clock::now();
        }
        results.opt_true_errors.resize(n_vertices);    // Changed from true_errors
        for (size_t i = 0; i < n_vertices; i++) {
            results.opt_true_errors[i] = std::abs(y_true[i] - results.opt_predictions[i]);
        }
        if (verbose) elapsed_time(ptm, "DONE");
    } else {
        results.opt_true_errors.clear();  // Ensure empty if no true values
    }

    if (verbose) {
        double total_time = ((double)clock() / CLOCKS_PER_SEC) - total_ptm;
        Rprintf("\nTotal elapsed time: %.2f seconds\n", total_time);
    }

    return results;
}
#endif

/**
*  Univariate Generalized Path-Based Local Models with Global Adaptive Neighborhood Size
#'
*  @description
*  Fits a generalized path-based local model (GPLM) to univariate data with automatic neighborhood
*  size selection and uncertainty quantification via bootstrap. This function is a specialized
*  version of GPLM for one-dimensional predictor variables, automatically creating a chain graph
*  structure based on the ordered x values.
#'
*  @param x Numeric vector of predictor values, must be sorted in ascending order
*  @param y Numeric vector of response values
*  @param y_true Optional numeric vector of true response values for error calculation
*  @param max_distance_deviation Integer specifying maximum allowed deviation from optimal path distance
*  @param use_median Logical indicating whether to use median instead of mean for aggregation (default: FALSE)
*  @param h_min Minimum neighborhood size to consider, must be odd and >= 3 (default: 3)
*  @param h_max Maximum neighborhood size to consider, must be odd and > h_min (default: 31)
*  @param p Probability level for credible intervals, must be between 0 and 1 (default: 0.95)
*  @param n_bb Number of bootstrap replicates (default: 500)
*  @param ikernel Integer specifying kernel type for local weighting (default: 1)
*  @param n_cores Number of cores for parallel processing (default: 1)
*  @param dist_normalization_factor Factor for normalizing path distances (default: 1.01)
*  @param epsilon Numerical tolerance for computations (default: 1e-15)
*  @param seed Random seed for reproducibility (default: 0)
*  @param verbose Logical indicating whether to print progress information (default: TRUE)
*
*  @return An object of class 'adaptive_nbhd_size_global_t' containing:
*    - opt_graph: Path graph for the optimal h value
*    - h_values: Vector of tested neighborhood sizes
*    - opt_h: Optimal neighborhood size based on CV errors
*    - h_cv_errors: Vector of mean cross-validation errors for each h
*    - opt_predictions: Vector of median bootstrap predictions for optimal h
*    - opt_cv_errors: Vector of per-vertex CV errors for optimal h
*    - opt_true_errors: Vector of true prediction errors for optimal h (if y_true provided)
*    - opt_ci_lower: Vector of lower bounds for credible intervals
*    - opt_ci_upper: Vector of upper bounds for credible intervals
*
*  @details
*  The function implements the following steps:
*  1. Creates a chain graph from sorted x values
*  2. Tests multiple neighborhood sizes from h_min to h_max
*  3. For each h:
*     - Fits local linear models along paths
*     - Computes bootstrap replicates for uncertainty quantification
*     - Calculates cross-validation errors
*  4. Selects optimal h based on CV errors
*  5. Computes credible intervals from bootstrap replicates
#'
*  @section Parallel Processing:
*  When n_cores > 1, the function uses OpenMP for parallel bootstrap computations.
*  Set n_cores to the number of available CPU cores for maximum performance.
#'
*  @examples
*  \dontrun{
*  # Generate sample data
*  x <- sort(runif(100))
*  y <- sin(2*pi*x) + rnorm(100, 0, 0.1)
#'
*  # Fit model with default parameters
*  result <- univariate_gplm_global(x, y, y_true = numeric(0),
*                                   max_distance_deviation = 2)
#'
*  # Access results
*  optimal_h <- result$opt_h
*  predictions <- result$opt_predictions
*  lower_ci <- result$opt_ci_lower
*  upper_ci <- result$opt_ci_upper
*  }
#'
*  @seealso
*  \code{\link{gplm_global}} for the general graph version
#'
*  @references
*  Add relevant references to the methodology
#'
*  @note
*  - Input x values must be sorted in ascending order
*  - The number of observations must be at least h_min
*  - Parallel processing requires OpenMP support
#'
*  @export
*/
adaptive_nbhd_size_global_t univariate_gplm_global(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& y_true,
    int max_distance_deviation,
    int h_min = 3,
    int h_max = 31,
    double p = 0.95,
    int n_bb = 500,
    int ikernel = 1,
    int n_cores = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15,
    unsigned int seed = 0,
    bool verbose = true) {

    // Input validation
    if (x.empty() || y.empty()) {
        Rf_error("Input vectors x and y cannot be empty");
    }
    if (x.size() != y.size()) {
        Rf_error("Input vectors x and y must have the same size");
    }
    // Add validation for y_true if provided
    if (!y_true.empty() && y_true.size() != y.size()) {
        Rf_error("y_true must have the same size as y");
    }
    // Add validation for x values
    if (x.size() < h_min) {
        Rf_error("Number of observations must be at least h_min");
    }
    // Verify x is sorted
    if (!std::is_sorted(x.begin(), x.end())) {
        Rf_error("x values must be sorted in ascending order");
    }

    auto [x_graph, x_edge_lengths] = create_chain_graph(x);
    return gplm_global(x_graph,
                      x_edge_lengths,
                      y,
                      y_true,
                      max_distance_deviation,
                      h_min,
                      h_max,
                      p,
                      n_bb,
                      ikernel,
                      n_cores,
                      dist_normalization_factor,
                      epsilon,
                      seed,
                      verbose);
}


/**
 * @brief R/C++ interface function for univariate generalized path-based linear models
 *
 * This function provides the interface between R and C++ implementations, handling:
 * 1. Conversion of R objects to C++ types
 * 2. Memory protection for R objects
 * 3. Call to C++ implementation
 * 4. Conversion of results back to R objects
 * 5. Creation of named R list for return values
 *
 * @param s_x SEXP (NumericVector) Sorted predictor values
 * @param s_y SEXP (NumericVector) Response values
 * @param s_y_true SEXP (NumericVector) True response values (optional)
 * @param s_max_distance_deviation SEXP (IntegerVector) Maximum allowed path deviation
 * @param s_use_median SEXP (LogicalVector) Whether to use median for aggregation
 * @param s_h_min SEXP (IntegerVector) Minimum neighborhood size
 * @param s_h_max SEXP (IntegerVector) Maximum neighborhood size
 * @param s_p SEXP (NumericVector) Credible interval probability
 * @param s_n_bb SEXP (IntegerVector) Number of bootstrap replicates
 * @param s_ikernel SEXP (IntegerVector) Kernel type for weighting
 * @param s_n_cores SEXP (IntegerVector) Number of parallel threads
 * @param s_dist_normalization_factor SEXP (NumericVector) Distance normalization
 * @param s_epsilon SEXP (NumericVector) Computational tolerance
 * @param s_seed SEXP (IntegerVector) Random seed
 * @param s_verbose SEXP (LogicalVector) Progress output control
 *
 * @return SEXP (List) A named list containing:
 *    - h_values: Integer vector of tested neighborhood sizes
 *    - h_cv_errors: Numeric vector of cross-validation errors for each h
 *    - h_median_predictions: List of numeric vectors, median predictions for each h
 *    - h_vertex_cv_errors: List of matrices, CV errors for each h organized by vertex and bootstrap
 *    - opt_h: Optimal neighborhood size
 *    - opt_predictions: Numeric vector of predictions at optimal h
 *    - opt_cv_errors: Numeric vector of CV errors at optimal h
 *    - opt_ci_lower: Numeric vector of lower credible bounds
 *    - opt_ci_upper: Numeric vector of upper credible bounds
 *    - opt_true_error: Mean absolute true error (if y_true provided)
 *
 * @note Input Requirements:
 *    - s_x must contain sorted values
 *    - s_x and s_y must have same length
 *    - s_y_true (if provided) must match length of s_y
 *    - All scalar parameters must be length 1
 *
 * @throws Rf_error if:
 *    - Input validation fails in C++ implementation
 *    - Memory allocation fails
 *    - Type conversion errors occur
 *
 * @see univariate_gplm_global() For the underlying C++ implementation
 * @see convert_vector_int_to_R() For integer vector conversion
 * @see convert_vector_double_to_R() For double vector conversion
 *
 * @warning
 *    - Does not include opt_graph component in return value
 *    - R objects must be properly protected during all operations
 *    - Memory protection count must match number of PROTECT calls
 */
SEXP S_univariate_gplm_global(SEXP s_x,
                              SEXP s_y,
                              SEXP s_y_true,
                              SEXP s_max_distance_deviation,
                              SEXP s_h_min,
                              SEXP s_h_max,
                              SEXP s_p,
                              SEXP s_n_bb,
                              SEXP s_ikernel,
                              SEXP s_n_cores,
                              SEXP s_dist_normalization_factor,
                              SEXP s_epsilon,
                              SEXP s_seed,
                              SEXP s_verbose) {

    // Convert R inputs to C++
    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(s_y_true) == n_points) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    // Convert scalar parameters
    int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    int h_min = INTEGER(s_h_min)[0];
    int h_max = INTEGER(s_h_max)[0];
    double p = REAL(s_p)[0];
    int n_bb = INTEGER(s_n_bb)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    // Call C++ implementation
    auto cpp_results = univariate_gplm_global(x,
                                            y,
                                            y_true,
                                            max_distance_deviation,
                                            h_min,
                                            h_max,
                                            p,
                                            n_bb,
                                            ikernel,
                                            n_cores,
                                            dist_normalization_factor,
                                            epsilon,
                                            seed,
                                            verbose);

    // Create return list with updated number of components
    int n_protected = 0;
    const int N_COMPONENTS = 10;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Updated component indices
    const int IDX_H_VALUES = 0;
    const int IDX_H_CV_ERRORS = 1;
    const int IDX_H_MEDIAN_PREDICTIONS = 2;
    const int IDX_H_VERTEX_CV_ERRORS = 3;
    const int IDX_OPT_H = 4;
    const int IDX_OPT_PREDICTIONS = 5;
    const int IDX_OPT_CV_ERRORS = 6;
    const int IDX_OPT_CI_LOWER = 7;
    const int IDX_OPT_CI_UPPER = 8;
    const int IDX_OPT_TRUE_ERROR = 9;

    // Convert and set components
    SET_VECTOR_ELT(result, IDX_H_VALUES,
        PROTECT(convert_vector_int_to_R(cpp_results.h_values))); n_protected++;

    SET_VECTOR_ELT(result, IDX_H_CV_ERRORS,
        PROTECT(convert_vector_double_to_R(cpp_results.h_cv_errors))); n_protected++;

    SET_VECTOR_ELT(result, IDX_H_MEDIAN_PREDICTIONS,
        PROTECT(convert_vector_vector_double_to_R(cpp_results.h_median_predictions))); n_protected++;

    SET_VECTOR_ELT(result, IDX_H_VERTEX_CV_ERRORS,
        PROTECT(convert_vector_vector_vector_double_to_R(cpp_results.h_vertex_cv_errors))); n_protected++;

    SEXP s_opt_h = PROTECT(allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h)[0] = cpp_results.opt_h;
    SET_VECTOR_ELT(result, IDX_OPT_H, s_opt_h);

    SET_VECTOR_ELT(result, IDX_OPT_PREDICTIONS,
        PROTECT(convert_vector_double_to_R(cpp_results.opt_predictions))); n_protected++;

    SET_VECTOR_ELT(result, IDX_OPT_CV_ERRORS,
        PROTECT(convert_vector_double_to_R(cpp_results.opt_cv_errors))); n_protected++;

    SET_VECTOR_ELT(result, IDX_OPT_CI_LOWER,
        PROTECT(convert_vector_double_to_R(cpp_results.opt_ci_lower))); n_protected++;

    SET_VECTOR_ELT(result, IDX_OPT_CI_UPPER,
        PROTECT(convert_vector_double_to_R(cpp_results.opt_ci_upper))); n_protected++;

    if (!cpp_results.opt_true_errors.empty()) {
        SEXP s_true_error = PROTECT(allocVector(REALSXP, 1)); n_protected++;
        REAL(s_true_error)[0] = std::accumulate(cpp_results.opt_true_errors.begin(),
                                               cpp_results.opt_true_errors.end(), 0.0) /
                                cpp_results.opt_true_errors.size();
        SET_VECTOR_ELT(result, IDX_OPT_TRUE_ERROR, s_true_error);
    } else {
        SET_VECTOR_ELT(result, IDX_OPT_TRUE_ERROR, R_NilValue);
    }

    // Set names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, IDX_H_VALUES, mkChar("h_values"));
    SET_STRING_ELT(names, IDX_H_CV_ERRORS, mkChar("h_cv_errors"));
    SET_STRING_ELT(names, IDX_H_MEDIAN_PREDICTIONS, mkChar("h_median_predictions"));
    SET_STRING_ELT(names, IDX_H_VERTEX_CV_ERRORS, mkChar("h_vertex_cv_errors"));
    SET_STRING_ELT(names, IDX_OPT_H, mkChar("opt_h"));
    SET_STRING_ELT(names, IDX_OPT_PREDICTIONS, mkChar("opt_predictions"));
    SET_STRING_ELT(names, IDX_OPT_CV_ERRORS, mkChar("opt_cv_errors"));
    SET_STRING_ELT(names, IDX_OPT_CI_LOWER, mkChar("opt_ci_lower"));
    SET_STRING_ELT(names, IDX_OPT_CI_UPPER, mkChar("opt_ci_upper"));
    SET_STRING_ELT(names, IDX_OPT_TRUE_ERROR, mkChar("opt_true_error"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return result;
}
