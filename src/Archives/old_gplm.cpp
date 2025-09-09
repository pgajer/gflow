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

#include <chrono>
#include <execution>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <omp.h>
#include <unordered_map>

#include "error_utils.h"
#include "gplm.h"
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
    SEXP S_graph_kpath_lm(SEXP s_path_graph,
                          SEXP s_y,
                          SEXP s_ikernel,
                          SEXP s_dist_normalization_factor);

    SEXP S_graph_kpath_lm_cv(SEXP s_path_graph,
                             SEXP s_y,
                             SEXP s_ikernel,
                             SEXP s_dist_normalization_factor,
                             SEXP s_n_CVs,
                             SEXP s_n_CV_folds,
                             SEXP s_epsilon,
                             SEXP s_seed);

    SEXP S_graph_kpath_lm_flexible(SEXP s_path_graph,
                               SEXP s_y,
                               SEXP s_ikernel,
                               SEXP s_max_distance_deviation,
                               SEXP s_min_required_paths,
                               SEXP s_dist_normalization_factor,
                                   SEXP s_epsilon);

    SEXP S_graph_kpath_lm_flexible_cv(SEXP s_path_graph,
                                      SEXP s_y,
                                      SEXP s_ikernel,
                                      SEXP s_max_distance_deviation,
                                      SEXP s_min_required_paths,
                                      SEXP s_dist_normalization_factor,
                                      SEXP s_epsilon,
                                      SEXP s_n_CVs,
                                      SEXP s_n_CV_folds,
                                      SEXP s_seed);
}


/**
 * @brief Computes local linear models along paths in a graph using kernel weights
 *
 * For each vertex in the graph, this function:
 * 1. Finds all paths of length h containing the vertex
 * 2. Among these paths, selects those where the vertex is closest to the middle
 * 3. For each selected path:
 *    - Maps vertices to positions along the path using edge weights as distances
 *    - Computes kernel weights based on distances from the target vertex
 *    - Fits a local linear model using the kernel weights
 * 4. Averages predictions from all selected paths for the vertex
 *
 * @param path_graph PLM graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param ikernel Integer specifying the kernel type for distance-based weighting
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 *
 * @return Vector of fitted values for each vertex in the graph
 *
 * @throws Rf_error if:
 *         - Path length h is not positive or not odd
 *         - Invalid kernel type is specified
 *         - Length of y doesn't match number of vertices
 *         - Distance normalization factor is not positive
 *         - No valid paths are found for any vertex
 *
 * @note The model uses kernel weights based on the distance between vertices along each path.
 *       These distances are computed using the edge weights stored in the path_graph structure.
 *
 * @see predict_lm_1d for details on the local linear model fitting
 * @see path_graph_plm_t for the graph structure definition
 * @see graph_kpath_lm_with_weights for a version that supports additional vertex weights (e.g., for cross-validation)
 */
std::vector<double> graph_kpath_lm(const path_graph_plm_t& path_graph,
                                   const std::vector<double>& y,
                                   int ikernel,
                                   double dist_normalization_factor = 1.01) {

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }

    std::vector<double> Ey(n_vertices); // output vector of the mean local linear models

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Extracting full paths of length h, with the position of the given vertex within the path
        std::vector<std::pair<std::vector<int>, int>> vertex_paths = vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        int min_dist_to_mid_pt = h;
        std::unordered_set<int> indices_of_min_dist_to_mid_pt;
        int mid_pt = (h - 1) / 2;
        int dist_to_mid_pt;
        int n_vertex_paths = vertex_paths.size();

        for (int path_i = 0; path_i < n_vertex_paths; ++path_i) {
            if ((dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt)) < min_dist_to_mid_pt) {
                min_dist_to_mid_pt = dist_to_mid_pt;
                indices_of_min_dist_to_mid_pt.clear();
                indices_of_min_dist_to_mid_pt.insert(path_i);
            } else if (dist_to_mid_pt == min_dist_to_mid_pt) {
                indices_of_min_dist_to_mid_pt.insert(path_i);
            }
        }

        if (indices_of_min_dist_to_mid_pt.empty()) {
            Rf_error("No valid paths were found");
        }

        double vertex_Ey = 0;
        for (const auto& path_i : indices_of_min_dist_to_mid_pt) {
            std::vector<int> path = vertex_paths[path_i].first;
            int path_n_vertices   = path.size();
            int position_in_path  = vertex_paths[path_i].second;

            // computing x, y and w for a weighted linear 1d model
            std::vector<double> y_path(path_n_vertices); // y restricted to the path
            std::vector<double> x_path(path_n_vertices); // distance from the initial vertex of the path to each vertex within the path
            std::vector<double> d_path(path_n_vertices); // distance from the reference vertex to any other vertex - used in computing kernel weights
            x_path[0] = 0;
            for (int i = 0; i < path_n_vertices; ++i) {
                y_path[i] = y[path[i]];
                if (i > 0) {
                    auto neighbors = path_graph.adj_list[path[i - 1]]; // neighbors of path[i - 1]
                    auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                    // finding path[i] in neighbors
                    auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                    int neighbor_i = it - neighbors.begin(); // index of path[i] in neighbors
                    x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i]; // neighbor_weights[neighbor_i] is the weight/distance from path[i-1] to path[i]
                }
            }

            // Compute and apply kernel weights
            std::vector<double> w_path(path_n_vertices); // kernel weights
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[position_in_path]);
                if (d_path[i] > max_dist)
                    max_dist = d_path[i];
            }
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;
            for (int i = 0; i < path_n_vertices; ++i)
                d_path[i] /= max_dist;

            // creating kernel weights
            double scale = 1.0;
            initialize_kernel(ikernel, scale);
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            vertex_Ey += predict_lm_1d(y_path, x_path, w_path, position_in_path);
        }

        Ey[vertex_i] = vertex_Ey / indices_of_min_dist_to_mid_pt.size();
    }

    return Ey;
}

/**
 * @brief R interface for computing local linear models along paths in a graph
 *
 * @details This function serves as the R interface wrapper for the graph_kpath_lm algorithm.
 * It processes R objects (SEXP) into C++ data structures, calls the core implementation,
 * and converts the results back into an R object. The function implements local linear
 * model fitting along paths in a graph structure, where kernel weights are used based
 * on distances between vertices.
 *
 * @param s_path_graph SEXP containing the path graph structure with the following components:
 *        - Adjacency lists for vertices
 *        - Edge weights
 *        - Path information including length h
 * @param s_y SEXP containing numeric vector of response variables (one per vertex)
 * @param s_ikernel SEXP containing integer specifying the kernel type for weighting
 * @param s_dist_normalization_factor SEXP containing numeric value for distance normalization
 *
 * @return SEXP (REALSXP) containing vector of fitted values for each vertex
 *
 * @throws R error (via Rf_error) if:
 *         - Input conversions fail
 *         - Underlying graph_kpath_lm function throws an error
 *         - Memory allocation fails
 *
 * @note This function handles all necessary R object protection and memory management
 *       required for interfacing with R's garbage collector.
 *
 * @see graph_kpath_lm for the core implementation details
 * @see sexp_to_path_graph_plm for the conversion of R objects to C++ graph structures
 */
SEXP S_graph_kpath_lm(SEXP s_path_graph,
                      SEXP s_y,
                      SEXP s_ikernel,
                      SEXP s_dist_normalization_factor) {

    path_graph_plm_t path_graph = sexp_to_path_graph_plm(s_path_graph);
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    // Compute result
    std::vector<double> result = graph_kpath_lm(path_graph,
                                                y,
                                                ikernel,
                                                dist_normalization_factor);

    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result.size()));
    for (size_t i = 0; i < result.size(); ++i) {
        REAL(s_result)[i] = result[i];
    }
    UNPROTECT(1);

    return s_result;
}

/**
 * @brief Computes local linear models along paths in a graph with optional weights
 *
 * For each vertex in the graph, this function:
 * 1. Finds all paths of length h containing the vertex
 * 2. Among these paths, selects those where the vertex is closest to the middle
 * 3. For each selected path:
 *    - Computes distances between vertices along the path
 *    - Applies kernel weights based on these distances
 *    - Combines kernel weights with input weights (if provided)
 *    - Fits a local linear model
 * 4. Averages predictions from all valid paths for the vertex
 *
 * @param path_graph PLM graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param weights Optional vector of weights for each vertex. If empty, uniform weights are used.
 *                Must have same length as y if non-empty.
 * @param ikernel Integer specifying the kernel type for distance-based weighting
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 *
 * @return A pair containing:
 *         - Vector of fitted values for each vertex
 *         - Vector of indices of vertices where model fitting failed
 *
 * @throws Rf_error if:
 *         - Path length h is not positive or not odd
 *         - Invalid kernel type is specified
 *         - Input vector lengths don't match
 *         - Distance normalization factor is not positive
 *
 * @note Vertices are excluded from the model (added to excluded_vertices) if:
 *       - No valid paths of length h contain the vertex
 *       - All weights (kernel * input) are zero for all paths
 *       - Local linear model fitting fails for all paths
 *
 * @see predict_lm_1d for details on the local linear model fitting
 * @see path_graph_plm_t for the graph structure definition
 */
std::pair<std::vector<double>, std::vector<int>> graph_kpath_lm_with_weights(const path_graph_plm_t& path_graph,
                                                                             const std::vector<double>& y,
                                                                             const std::vector<double>& weights,
                                                                             int ikernel,
                                                                             double dist_normalization_factor = 1.01) {

    #define DEBUG__graph_kpath_lm_with_weights 0

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }
    if (!weights.empty() && weights.size() != n_vertices) {
        Rf_error("Length of weights must be 0 or match number of vertices");
    }

    std::vector<double> Ey(n_vertices, 0.0); // output vector of the path linear models
    std::vector<int> excluded_vertices; // vertices where model couldn't be fit

    #if DEBUG__graph_kpath_lm_with_weights
    Rprintf("\nIn graph_kpath_lm_with_weights()\n");
    Rprintf("h: %d\n", h);
    print_vect(weights, "weights");
    #endif

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            excluded_vertices.push_back(vertex_i);
            continue;
        }

        int min_dist_to_mid_pt = h;
        std::unordered_set<int> indices_of_min_dist_to_mid_pt;
        int mid_pt = (h - 1) / 2;
        int dist_to_mid_pt;
        int n_vertex_paths = vertex_paths.size();

        #if DEBUG__graph_kpath_lm_with_weights
        Rprintf("vertex_i: %d\n", vertex_i);
        Rprintf("vertex_paths\n");
        for (int path_i = 0; path_i < n_vertex_paths; ++path_i) {
            print_vect(vertex_paths[path_i].first, "vertex_path");
            Rprintf("vertex_index: %d\n", vertex_paths[path_i].second);
        }
        #endif

        for (int path_i = 0; path_i < n_vertex_paths; ++path_i) {
            if ((dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt)) < min_dist_to_mid_pt) {
                min_dist_to_mid_pt = dist_to_mid_pt;
                indices_of_min_dist_to_mid_pt.clear();
                indices_of_min_dist_to_mid_pt.insert(path_i);
            } else if (dist_to_mid_pt == min_dist_to_mid_pt) {
                indices_of_min_dist_to_mid_pt.insert(path_i);
            }
        }

        #if DEBUG__graph_kpath_lm_with_weights
        print_uset(indices_of_min_dist_to_mid_pt, "indices_of_min_dist_to_mid_pt");
        #endif

        if (indices_of_min_dist_to_mid_pt.empty()) {
            excluded_vertices.push_back(vertex_i);
            continue;
        }

        bool any_valid_path = false;
        double vertex_Ey = 0;
        int valid_paths_count = 0;

        for (const auto& path_i : indices_of_min_dist_to_mid_pt) {
            try {
                std::vector<int> path = vertex_paths[path_i].first;
                int path_n_vertices   = path.size();
                int position_in_path  = vertex_paths[path_i].second;

                #if DEBUG__graph_kpath_lm_with_weights
                Rprintf("path_i: %d\n", path_i);
                print_vect(path, "path");
                Rprintf("path_n_vertices: %d\n", path_n_vertices);
                Rprintf("position_in_path: %d\n", position_in_path);
                #endif

                std::vector<double> y_path(path_n_vertices);
                std::vector<double> x_path(path_n_vertices);
                std::vector<double> d_path(path_n_vertices);
                std::vector<double> w_path(path_n_vertices, 1.0);

                x_path[0] = 0;
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

                #if DEBUG__graph_kpath_lm_with_weights
                print_vect(x_path, "x_path");
                #endif

                if (weights.size() == n_vertices) {
                    for (int i = 0; i < path_n_vertices; ++i)
                        w_path[i] = weights[path[i]];
                }

                #if DEBUG__graph_kpath_lm_with_weights
                print_vect(w_path, "w_path");
                #endif

                std::vector<double> kernel_w_path(path_n_vertices);
                double max_dist = 0.0;
                for (int i = 0; i < path_n_vertices; ++i) {
                    d_path[i] = std::abs(x_path[i] - x_path[position_in_path]);
                    if (d_path[i] > max_dist)
                        max_dist = d_path[i];
                }
                if (max_dist == 0) max_dist = 1;
                max_dist *= dist_normalization_factor;
                for (int i = 0; i < path_n_vertices; ++i)
                    d_path[i] /= max_dist;

                #if DEBUG__graph_kpath_lm_with_weights
                print_vect(d_path, "d_path");
                #endif

                double scale = 1.0;
                initialize_kernel(ikernel, scale);
                kernel_fn(d_path.data(), path_n_vertices, kernel_w_path.data());

                #if DEBUG__graph_kpath_lm_with_weights
                print_vect(kernel_w_path, "kernel_w_path");
                #endif

                for (int i = 0; i < path_n_vertices; ++i)
                    w_path[i] *= kernel_w_path[i];

                // Check if all weights are zero
                if (std::all_of(w_path.begin(), w_path.end(),
                    [](double w) { return w == 0.0; })) {
                    continue;
                }

                // Check if all weights non-negative
                for (int i = 0; i < path_n_vertices; ++i) {
                    if (w_path[i] < 0) {
                        Rprintf("path_i: %d\n", path_i);
                        print_vect(w_path, "w_path");
                        print_vect(path, "path");
                        Rprintf("position_in_path: %d\n", position_in_path);
                        Rf_error("Negative weights are not allowed");
                    }
                }

                vertex_Ey += predict_lm_1d(y_path, x_path, w_path, position_in_path);
                valid_paths_count++;
                any_valid_path = true;
            }
            catch (const std::runtime_error&) {
                // Continue to next path if this one fails
                continue;
            }
        }

        if (!any_valid_path) {
            excluded_vertices.push_back(vertex_i);
        } else {
            Ey[vertex_i] = vertex_Ey / valid_paths_count;
        }
    }

    return std::make_pair(Ey, excluded_vertices);
}

/**
 * @brief Performs cross-validation for local linear models on graph paths
 *
 * @details Implements k-fold cross-validation for graph-based local linear models.
 * For each iteration, randomly selects test vertices, fits models using remaining vertices
 * as training data, and computes prediction errors for test vertices.
 *
 * @param path_graph Graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param ikernel Integer specifying kernel type for distance-based weighting (default: 1)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param n_CVs Number of cross-validation iterations to perform
 * @param n_CV_folds Number of folds for cross-validation (default: 10)
 * @param epsilon Small value to prevent numerical issues (default: 1e-10)
 * @param seed Random seed for reproducibility (default: 0 uses system time)
 *
 * @return Vector of cross-validation errors for each vertex (NaN for vertices never in test set)
 *
 * @throws std::invalid_argument if n_CVs <= 0
 *
 * @note
 * - If seed = 0, system time is used to initialize random number generator
 * - For n_vertices == n_CVs and fold_size == 1, implements leave-one-out CV
 * - Vertices that cannot be predicted (e.g., due to insufficient paths) get NaN errors
 *
 * @see graph_kpath_lm_with_weights for the underlying model fitting
 */
std::vector<double> graph_kpath_lm_cv(const path_graph_plm_t& path_graph,
                                      const std::vector<double>& y,
                                      int ikernel = 1,
                                      double dist_normalization_factor = 1.01,
                                      int n_CVs = 0,
                                      int n_CV_folds = 10,
                                      unsigned int seed = 0) {

    #define DEBUG__graph_kpath_lm_cv 0

    if (n_CVs <= 0) {
        Rf_error("n_CVs has to be greater than 0");
    }

    int n_vertices = y.size();

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

    #if DEBUG__graph_kpath_lm_cv
    Rprintf("\nIn graph_kpath_lm_cv()\n");
    #endif

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

        #if DEBUG__graph_kpath_lm_cv
        print_set(test_set, "test_set");
        #endif

        // Reseting all weights to 1 and set weights over test_set to 0
        std::fill(weights.begin(), weights.end(), 1.0);
        for (const auto& vertex : test_set) {
            weights[vertex] = 0.0;
        }

        #if DEBUG__graph_kpath_lm_cv
        print_vect(weights, "weights");
        #endif

        // Estimating the conditional expectation of cv_y using linear model with path means
        auto res = graph_kpath_lm_with_weights(path_graph,
                                               y,
                                               weights,
                                               ikernel,
                                               dist_normalization_factor);

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

        // Computing cross-validation error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double error = std::abs(Ecv_y[vertex] - y[vertex]);

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = error;
            } else {
                cv_error[vertex] += error;
            }
            cv_error_count[vertex]++;
        }
    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV error, leaving NaN for vertices with no estimates
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex] > 0) {
            cv_error[vertex] /= cv_error_count[vertex];
        }
    }

    return cv_error;
}


/**
 * @brief R interface for graph path local linear model cross-validation
 *
 * @details Converts R objects to C++ types and calls graph_kpath_lm_cv to perform
 * cross-validation of local linear models on graph paths.
 *
 * @param s_path_graph SEXP containing path graph structure
 * @param s_y SEXP containing response variables
 * @param s_ikernel SEXP containing kernel type (integer)
 * @param s_dist_normalization_factor SEXP containing distance normalization factor
 * @param s_n_CVs SEXP containing number of CV iterations
 * @param s_n_CV_folds SEXP containing number of CV folds
 * @param s_epsilon SEXP containing numerical stability parameter
 * @param s_seed SEXP containing random seed
 *
 * @return SEXP (REALSXP) containing vector of cross-validation errors
 *
 * @throws R error via Rf_error for invalid inputs
 *
 * @note
 * - Handles all necessary R object protection
 * - Returns NaN for vertices that couldn't be predicted
 * - Uses R's random number generator when seed = 0
 *
 * @see graph_kpath_lm_cv for implementation details
 */
SEXP S_graph_kpath_lm_cv(SEXP s_path_graph,
                         SEXP s_y,
                         SEXP s_ikernel,
                         SEXP s_dist_normalization_factor,
                         SEXP s_n_CVs,
                         SEXP s_n_CV_folds,
                         SEXP s_epsilon,
                         SEXP s_seed) {

    path_graph_plm_t path_graph = sexp_to_path_graph_plm(s_path_graph);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    double epsilon = REAL(s_epsilon)[0];
    epsilon = epsilon + 0.0;
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    std::vector<double> cv_errors =  graph_kpath_lm_cv(path_graph,
                                                       y,
                                                       ikernel,
                                                       dist_normalization_factor,
                                                       n_CVs,
                                                       n_CV_folds,
                                                       seed);
    // Convert result to SEXP and return
    SEXP s_cv_errors = PROTECT(allocVector(REALSXP, cv_errors.size()));
    for (size_t i = 0; i < cv_errors.size(); ++i) {
        REAL(s_cv_errors)[i] = cv_errors[i];
    }
    UNPROTECT(1);
    return s_cv_errors;
}


/**
 * @brief Computes MAE-weighted local linear models along paths in a graph
 *
 * This function extends graph_kpath_lm by weighting predictions from different paths
 * based on their Mean Absolute Error (MAE). Paths with lower MAE get higher weights
 * in the final prediction, improving the robustness of the estimates.
 *
 * For each vertex in the graph, this function:
 * 1. Finds all paths of length h containing the vertex
 * 2. Among these paths, selects those where the vertex is closest to the middle
 * 3. For each selected path:
 *    - Maps vertices to positions along the path using edge weights as distances
 *    - Computes kernel weights based on distances from the target vertex
 *    - Fits a local linear model using the kernel weights
 *    - Computes the MAE for the fitted model
 * 4. Combines predictions using weights w_i = (1/mae_i) / sum(1/mae_j)
 *
 * @param path_graph PLM graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param ikernel Integer specifying the kernel type for distance-based weighting
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return Vector of fitted values for each vertex in the graph
 *
 * @throws Rf_error if:
 *         - Path length h is not positive or not odd
 *         - Invalid kernel type is specified
 *         - Length of y doesn't match number of vertices
 *         - Distance normalization factor is not positive
 *         - No valid paths are found for any vertex
 *         - All MAE values for a vertex's paths are effectively zero
 *
 * @note The model uses two types of weights:
 *       1. Kernel weights based on the distance between vertices along each path
 *       2. Path weights based on the inverse of each path's MAE
 *
 * @see predict_lm_1d_with_mae for details on the local linear model fitting and MAE computation
 * @see path_graph_plm_t for the graph structure definition
 */
std::vector<double> graph_kpath_lm_weighted_mae(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int ikernel,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8) {

    int h = path_graph.h;
    int n_vertices = path_graph.vertex_paths.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }

    std::vector<double> Ey(n_vertices);

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get paths containing the vertex
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        // Find paths where vertex is closest to middle
        int min_dist_to_mid_pt = h;
        std::unordered_set<int> indices_of_min_dist_to_mid_pt;
        int mid_pt = (h - 1) / 2;

        for (int path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            int dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt);
            if (dist_to_mid_pt < min_dist_to_mid_pt) {
                min_dist_to_mid_pt = dist_to_mid_pt;
                indices_of_min_dist_to_mid_pt.clear();
                indices_of_min_dist_to_mid_pt.insert(path_i);
            } else if (dist_to_mid_pt == min_dist_to_mid_pt) {
                indices_of_min_dist_to_mid_pt.insert(path_i);
            }
        }

        if (indices_of_min_dist_to_mid_pt.empty()) {
            Rf_error("No valid paths were found");
        }

        // Store predictions and MAEs for each path
        std::vector<double> path_predictions;
        std::vector<double> path_maes;
        path_predictions.reserve(indices_of_min_dist_to_mid_pt.size());
        path_maes.reserve(indices_of_min_dist_to_mid_pt.size());

        // Compute predictions and MAEs for each path
        for (const auto& path_i : indices_of_min_dist_to_mid_pt) {
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

            double scale = 1.0;
            initialize_kernel(ikernel, scale);
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            // Fit model and get prediction with MAE
            auto fit_results = predict_lm_1d_with_mae(y_path, x_path, w_path, position_in_path, epsilon);
            path_predictions.push_back(fit_results.predicted_value);
            path_maes.push_back(fit_results.mae);
        }

        bool use_unweighted = false;
        double unweighted_prediction = 0.0;
        double sum_inverse_mae = 0.0;
        std::vector<double> mae_weights(path_maes.size());

        for (size_t i = 0; i < path_maes.size(); ++i) {
            if (path_maes[i] <= epsilon) {
                use_unweighted = true;
                unweighted_prediction = path_predictions[i];
                break;
            }
            mae_weights[i] = 1.0 / path_maes[i];
            sum_inverse_mae += mae_weights[i];
        }

        if (use_unweighted) {
            Ey[vertex_i] = unweighted_prediction;
            continue;  // Skip to next vertex
        }

        if (sum_inverse_mae <= epsilon) {
            Rf_error("Sum of inverse MAE weights is effectively zero for vertex %d", vertex_i);
        }

        // Compute weighted average prediction
        double weighted_prediction = 0.0;
        for (size_t i = 0; i < path_predictions.size(); ++i) {
            weighted_prediction += (mae_weights[i] / sum_inverse_mae) * path_predictions[i];
        }
        Ey[vertex_i] = weighted_prediction;
    }

    return Ey;
}


/**
 * @brief Computes MAE-weighted local linear models along paths with flexible vertex positioning
 *
 * This function extends graph_kpath_lm_weighted_mae by allowing paths where the vertex
 * position deviates from the optimal center position by a controlled amount. This increases
 * the number of paths used in prediction, potentially leading to more robust estimates.
 *
 * For each vertex in the graph, this function:
 * 1. Finds all paths of length h containing the vertex
 * 2. Among these paths:
 *    - Identifies the minimum distance to middle position
 *    - Selects paths where the vertex position deviates from optimal by at most max_distance_deviation
 *    - Ensures at least min_required_paths are selected (if available)
 * 3. For each selected path:
 *    - Maps vertices to positions using edge weights as distances
 *    - Computes kernel weights based on distances from the target vertex
 *    - Fits a local linear model using the kernel weights
 *    - Computes the MAE for the fitted model
 * 4. Combines predictions using weights w_i = (1/mae_i) / sum(1/mae_j)
 *
 * @param path_graph PLM graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param ikernel Integer specifying the kernel type for distance-based weighting
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param min_required_paths Minimum number of paths required for prediction
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return Vector of fitted values for each vertex in the graph
 *
 * @throws Rf_error if:
 *         - Path length h is not positive or not odd
 *         - Invalid kernel type is specified
 *         - Length of y doesn't match number of vertices
 *         - Distance normalization factor is not positive
 *         - max_distance_deviation is negative
 *         - min_required_paths is less than 1
 *         - Not enough valid paths are found for a vertex
 *         - All MAE values for a vertex's paths are effectively zero
 *
 * @note If fewer than min_required_paths are available even with max_distance_deviation,
 *       the function uses all available paths and issues a warning.
 */
std::vector<double> graph_kpath_lm_flexible(const path_graph_plm_t& path_graph,
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
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }

    std::vector<double> Ey(n_vertices);

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get all paths containing the vertex
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            Rf_error("No paths found for vertex %d", vertex_i);
        }

        // Find minimum distance to middle point
        int mid_pt = (h - 1) / 2;
        int min_dist_to_mid_pt = h;

        for (const auto& path : vertex_paths) {
            int dist_to_mid_pt = std::abs(path.second - mid_pt);
            min_dist_to_mid_pt = std::min(min_dist_to_mid_pt, dist_to_mid_pt);
        }

        // Select paths within allowed deviation
        std::vector<size_t> selected_path_indices;
        for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            int dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt);
            if (dist_to_mid_pt <= min_dist_to_mid_pt + max_distance_deviation) {
                selected_path_indices.push_back(path_i);
            }
        }

        if (selected_path_indices.empty()) {
            Rf_error("No valid paths found for vertex %d within deviation limits", vertex_i);
        }

        if (selected_path_indices.size() < min_required_paths) {
            Rprintf("Warning: Only %d paths found for vertex %d (minimum %d requested)\n",
                   (int)selected_path_indices.size(), vertex_i, min_required_paths);
        }

        // Store predictions and MAEs for selected paths
        std::vector<double> path_predictions;
        std::vector<double> path_maes;
        path_predictions.reserve(selected_path_indices.size());
        path_maes.reserve(selected_path_indices.size());

        // Compute predictions and MAEs for each selected path
        for (const auto& path_i : selected_path_indices) {
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

            double scale = 1.0;
            initialize_kernel(ikernel, scale);
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            // Fit model and get prediction with MAE
            auto fit_results = predict_lm_1d_with_mae(y_path, x_path, w_path, position_in_path, epsilon);
            path_predictions.push_back(fit_results.predicted_value);
            path_maes.push_back(fit_results.mae);
        }

        // Compute weights based on inverse MAE
        double sum_inverse_mae = 0.0;
        std::vector<double> mae_weights(path_maes.size());
        bool use_unweighted = false;
        double unweighted_prediction = 0.0;
        for (size_t i = 0; i < path_maes.size(); ++i) {
            if (path_maes[i] <= epsilon) {
                use_unweighted = true;
                unweighted_prediction = path_predictions[i];
                break;
            }
            mae_weights[i] = 1.0 / path_maes[i];
            sum_inverse_mae += mae_weights[i];
        }

        if (use_unweighted) {
            Ey[vertex_i] = unweighted_prediction;
            continue;  // Skip to next vertex
        }

        if (sum_inverse_mae <= epsilon) {
            Rf_error("Sum of inverse MAE weights is effectively zero for vertex %d", vertex_i);
        }

        // Compute weighted average prediction
        double weighted_prediction = 0.0;
        for (size_t i = 0; i < path_predictions.size(); ++i) {
            weighted_prediction += (mae_weights[i] / sum_inverse_mae) * path_predictions[i];
        }
        Ey[vertex_i] = weighted_prediction;
    }

    return Ey;
}


/**
 * @brief R interface for graph_kpath_lm_flexible function
 *
 * This function provides an R interface to the C++ implementation of the flexible path
 * local linear model algorithm. It converts R objects (SEXP) to C++ types, calls the
 * implementation function, and converts the results back to R format.
 *
 * @param s_path_graph SEXP containing the path graph structure (converted using sexp_to_path_graph_plm)
 * @param s_y SEXP containing numeric vector of response variables
 * @param s_ikernel SEXP containing integer specifying the kernel type
 * @param s_max_distance_deviation SEXP containing integer for maximum allowed deviation from optimal center position
 * @param s_min_required_paths SEXP containing integer for minimum number of paths required
 * @param s_dist_normalization_factor SEXP containing numeric value for distance normalization (default: 1.01)
 * @param s_epsilon SEXP containing numeric value for computational precision (default: 1e-8)
 *
 * @return SEXP containing numeric vector of fitted values for each vertex
 *
 * @throws R error if:
 *         - Any input SEXP is NULL or of incorrect type
 *         - Length of input vectors is incorrect
 *         - Parameters are out of valid ranges
 *         - Memory allocation fails
 *
 * @note This function handles R garbage collection through PROTECT/UNPROTECT
 */
SEXP S_graph_kpath_lm_flexible(SEXP s_path_graph,
                               SEXP s_y,
                               SEXP s_ikernel,
                               SEXP s_max_distance_deviation,
                               SEXP s_min_required_paths,
                               SEXP s_dist_normalization_factor,
                               SEXP s_epsilon) {

    if (!isReal(s_y) || !isInteger(s_ikernel) || !isInteger(s_max_distance_deviation) ||
        !isInteger(s_min_required_paths) || !isReal(s_dist_normalization_factor) ||
        !isReal(s_epsilon)) {
        error("Invalid argument type");
    }

    // Convert inputs
    path_graph_plm_t path_graph = sexp_to_path_graph_plm(s_path_graph);
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    int ikernel = INTEGER(s_ikernel)[0];
    int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    int min_required_paths = INTEGER(s_min_required_paths)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];

    // Compute result
    std::vector<double> result = graph_kpath_lm_flexible(path_graph,
                                                        y,
                                                        ikernel,
                                                        max_distance_deviation,
                                                        min_required_paths,
                                                        dist_normalization_factor,
                                                        epsilon);

    // Convert result to SEXP and return
    SEXP s_result = PROTECT(allocVector(REALSXP, result.size()));
    if (s_result == R_NilValue) {
        UNPROTECT(1);
        error("Memory allocation failed");
    }

    for (size_t i = 0; i < result.size(); ++i) {
        REAL(s_result)[i] = result[i];
    }

    UNPROTECT(1);
    return s_result;
}


/**
 * @brief Computes weighted local linear models along paths with flexible vertex positioning
 *
 * This function extends graph_kpath_lm_flexible by incorporating vertex-specific weights
 * in the local linear model fitting. It combines both kernel weights based on path distances
 * and user-provided vertex weights.
 *
 * @param path_graph Graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param weights Vector of non-negative weights for each vertex
 * @param ikernel Integer specifying the kernel type
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param min_required_paths Minimum number of paths required for prediction
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return std::pair containing:
 *         - Vector of fitted values for each vertex
 *         - Vector of indices of vertices where model couldn't be fit
 *
 * @throws Rf_error if:
 *         - Input vectors have inconsistent lengths
 *         - Parameters are out of valid ranges
 *         - No valid paths found for a vertex
 *         - Weights are negative
 *         - Numerical instability detected
 *
 * @note Final weights for each vertex in a path are products of:
 *       - User-provided vertex weights
 *       - Kernel weights based on path distances
 * @note Vertices are excluded from prediction if:
 *       - No valid paths contain the vertex
 *       - All path weights become zero
 *       - Numerical instability in weighted predictions
 */
std::pair<std::vector<double>, std::vector<int>>
graph_kpath_lm_flexible_with_weights(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
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
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }
    if (weights.size() != n_vertices) {
        Rf_error("Length of weights must match the number of vertices");
    }

    std::vector<double> Ey(n_vertices); // output vector of the predictions of the path linear models
    std::vector<int> excluded_vertices; // vertices where model couldn't be fit

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Get all paths containing the vertex
        std::vector<std::pair<std::vector<int>, int>> vertex_paths =
            vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        if (vertex_paths.empty()) {
            excluded_vertices.push_back(vertex_i);
            // Rf_error("No paths found for vertex %d", vertex_i);
        }

        // Find minimum distance to middle point
        int mid_pt = (h - 1) / 2;
        int min_dist_to_mid_pt = h;

        for (const auto& path : vertex_paths) {
            int dist_to_mid_pt = std::abs(path.second - mid_pt);
            min_dist_to_mid_pt = std::min(min_dist_to_mid_pt, dist_to_mid_pt);
        }

        // Select paths within allowed deviation
        std::vector<size_t> selected_path_indices;
        for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
            int dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt);
            if (dist_to_mid_pt <= min_dist_to_mid_pt + max_distance_deviation) {
                selected_path_indices.push_back(path_i);
            }
        }

        if (selected_path_indices.empty()) {
            excluded_vertices.push_back(vertex_i);
            //Rf_error("No valid paths found for vertex %d within deviation limits", vertex_i);
        }

        if (selected_path_indices.size() < min_required_paths) {
            Rprintf("Warning: Only %d paths found for vertex %d (minimum %d requested)\n",
                   (int)selected_path_indices.size(), vertex_i, min_required_paths);
        }

        // Store predictions and MAEs for selected paths
        std::vector<double> path_predictions;
        std::vector<double> path_maes;
        path_predictions.reserve(selected_path_indices.size());
        path_maes.reserve(selected_path_indices.size());

        // Compute predictions and MAEs for each selected path
        for (const auto& path_i : selected_path_indices) {
            std::vector<int> path = vertex_paths[path_i].first;
            int path_n_vertices = path.size();
            int position_in_path = vertex_paths[path_i].second;

            // Prepare data for linear model
            std::vector<double> y_path(path_n_vertices);
            std::vector<double> x_path(path_n_vertices);
            std::vector<double> d_path(path_n_vertices);
            std::vector<double> w_path(path_n_vertices, 1.0);
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

            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] = weights[path[i]];

            // Compute kernel weights
            std::vector<double> kernel_w_path(path_n_vertices);
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
            kernel_fn(d_path.data(), path_n_vertices, kernel_w_path.data());

            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] *= kernel_w_path[i];

            // Check if all weights are zero
            if (std::all_of(w_path.begin(), w_path.end(),
                            [](double w) { return w == 0.0; })) {
                excluded_vertices.push_back(vertex_i);
                continue;
            }

            // Fit model and get prediction with MAE
            auto fit_results = predict_lm_1d_with_mae(y_path, x_path, w_path, position_in_path, epsilon);
            path_predictions.push_back(fit_results.predicted_value);
            path_maes.push_back(fit_results.mae);
        }

        // Compute weights based on inverse MAE
        double sum_inverse_mae = 0.0;
        std::vector<double> mae_weights(path_maes.size());
        bool use_unweighted = false;
        double unweighted_prediction = 0.0;
        for (size_t i = 0; i < path_maes.size(); ++i) {
            if (path_maes[i] <= epsilon) {
                use_unweighted = true;
                unweighted_prediction = path_predictions[i];
                break;
            }
            mae_weights[i] = 1.0 / path_maes[i];
            sum_inverse_mae += mae_weights[i];
        }

        if (use_unweighted) {
            Ey[vertex_i] = unweighted_prediction;
            continue;  // Skip to next vertex
        }

        if (sum_inverse_mae <= epsilon) {
            Rf_error("Sum of inverse MAE weights is effectively zero for vertex %d", vertex_i);
        }

        // Compute weighted average prediction
        double weighted_prediction = 0.0;
        for (size_t i = 0; i < path_predictions.size(); ++i) {
            weighted_prediction += (mae_weights[i] / sum_inverse_mae) * path_predictions[i];
        }
        Ey[vertex_i] = weighted_prediction;
    }

    return std::make_pair(Ey, excluded_vertices);
}


/**
 * @brief Performs cross-validation for the flexible path local linear model
 *
 * This function implements k-fold cross-validation for graph_kpath_lm_flexible by repeatedly:
 * 1. Randomly selecting test vertices
 * 2. Setting their weights to zero (effectively excluding them from model fitting)
 * 3. Fitting the model using the remaining vertices
 * 4. Computing prediction errors for test vertices
 *
 * The cross-validation error is computed using absolute deviation loss function:
 * error = |predicted - observed|
 *
 * @param path_graph Graph structure containing adjacency lists, weights, and path information
 * @param y Vector of response variables (one per vertex)
 * @param ikernel Integer specifying kernel type (default: 1)
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param min_required_paths Minimum number of paths required for prediction
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 * @param n_CVs Number of cross-validation iterations to perform
 * @param n_CV_folds Number of folds for cross-validation (default: 10)
 * @param seed Random seed for reproducibility (default: 0, uses system time)
 *
 * @return Vector of average cross-validation errors for each vertex.
 *         NaN indicates vertices that were never successfully predicted.
 *
 * @note Special case: If fold_size=1 and n_vertices=n_CVs, performs leave-one-out CV
 *       systematically rather than random sampling
 *
 * @throws std::invalid_argument if:
 *         - n_CVs <= 0
 *         - Input vectors have inconsistent sizes
 *         - n_CV_folds <= 0 or > n_vertices
 *
 * @see graph_kpath_lm_flexible_with_weights for the underlying model fitting
 */
std::vector<double> graph_kpath_lm_flexible_cv(const path_graph_plm_t& path_graph,
                                              const std::vector<double>& y,
                                              int ikernel = 1,
                                              int max_distance_deviation = 1,
                                              int min_required_paths = 1,
                                              double dist_normalization_factor = 1.01,
                                              double epsilon = 1e-8,
                                              int n_CVs = 0,
                                              int n_CV_folds = 10,
                                              unsigned int seed = 0) {

    int n_vertices = y.size();
    if (n_vertices != path_graph.vertex_paths.size()) {
        Rf_error("Inconsistent sizes between y and path_graph");
    }

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

    // The main cross-validation loop
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

        // Estimating the conditional expectation of cv_y using linear model with path means
        auto res = graph_kpath_lm_flexible_with_weights(path_graph,
                                                        y,
                                                        weights,
                                                        ikernel,
                                                        max_distance_deviation,
                                                        min_required_paths,
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

        // Computing cross-validation error over test vertices using absolute deviation loss function
        for (const auto& vertex : valid_test_set) {
            double error = std::abs(Ecv_y[vertex] - y[vertex]);

            if (std::isnan(cv_error[vertex])) {
                cv_error[vertex] = error;
            } else {
                cv_error[vertex] += error;
            }
            cv_error_count[vertex]++;
        }
    } // END OF for (int cv = 0; cv < n_CVs; ++cv)

    // Compute average CV error, leaving NaN for vertices with no estimates
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        if (cv_error_count[vertex] > 0) {
            cv_error[vertex] /= cv_error_count[vertex];
        }
    }

    return cv_error;
}

/**
 * @brief R interface for cross-validation of flexible path local linear model
 *
 * This function provides an R interface to the C++ implementation of cross-validation
 * for the flexible path local linear model. It handles conversion between R and C++
 * data types and performs input validation.
 *
 * @param s_path_graph SEXP containing the path graph structure
 * @param s_y SEXP containing numeric vector of response variables
 * @param s_ikernel SEXP containing integer specifying kernel type
 * @param s_max_distance_deviation SEXP containing integer for maximum allowed deviation
 * @param s_min_required_paths SEXP containing integer for minimum required paths
 * @param s_dist_normalization_factor SEXP containing numeric normalization factor
 * @param s_epsilon SEXP containing numeric value for computational precision
 * @param s_n_CVs SEXP containing integer number of cross-validation iterations
 * @param s_n_CV_folds SEXP containing integer number of CV folds
 * @param s_seed SEXP containing integer random seed
 *
 * @return SEXP containing numeric vector of cross-validation errors
 *
 * @throws R error if:
 *         - Any input SEXP is NULL
 *         - Inputs have incorrect types
 *         - Memory allocation fails
 */
SEXP S_graph_kpath_lm_flexible_cv(SEXP s_path_graph,
                                 SEXP s_y,
                                 SEXP s_ikernel,
                                 SEXP s_max_distance_deviation,
                                 SEXP s_min_required_paths,
                                 SEXP s_dist_normalization_factor,
                                 SEXP s_epsilon,
                                 SEXP s_n_CVs,
                                 SEXP s_n_CV_folds,
                                 SEXP s_seed) {
    // Input validation
    if (!isReal(s_y) || !isInteger(s_ikernel) || !isInteger(s_max_distance_deviation) ||
        !isInteger(s_min_required_paths) || !isReal(s_dist_normalization_factor) ||
        !isReal(s_epsilon) || !isInteger(s_n_CVs) || !isInteger(s_n_CV_folds) ||
        !isInteger(s_seed)) {
        error("Invalid argument type");
    }
    if (s_path_graph == R_NilValue) {
        error("path_graph cannot be NULL");
    }

    // Convert inputs
    path_graph_plm_t path_graph = sexp_to_path_graph_plm(s_path_graph);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int ikernel = INTEGER(s_ikernel)[0];
    int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    int min_required_paths = INTEGER(s_min_required_paths)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    // Compute result
    std::vector<double> cv_errors = graph_kpath_lm_flexible_cv(path_graph,
                                                              y,
                                                              ikernel,
                                                              max_distance_deviation,
                                                              min_required_paths,
                                                              dist_normalization_factor,
                                                              epsilon,
                                                              n_CVs,
                                                              n_CV_folds,
                                                              seed);

    // Convert result to SEXP and return
    SEXP s_cv_errors = PROTECT(allocVector(REALSXP, cv_errors.size()));
    if (s_cv_errors == R_NilValue) {
        UNPROTECT(1);
        error("Memory allocation failed");
    }

    for (size_t i = 0; i < cv_errors.size(); ++i) {
        REAL(s_cv_errors)[i] = cv_errors[i];
    }

    UNPROTECT(1);
    return s_cv_errors;
}
