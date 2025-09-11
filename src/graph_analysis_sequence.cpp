#include "graph_diffusion_smoother.hpp"
#include "path_graphs.hpp"               // for path_graph_t and path_graph_plus_t
#include "MS_complex.h"                  // for MS_complex_t
#include "graph_diffusion_smoother.hpp"  // for graph_diffusion_smoother and related types
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"                 // for elapsed_time

#include <vector>
#include <map>
#include <utility>   // for std::pair
#include <memory>    // for std::unique_ptr
#include <algorithm> // for std::max_element
#include <stdexcept> // for std::invalid_argument

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {
    SEXP S_compute_graph_analysis_sequence(SEXP s_adj_list,
                                           SEXP s_weight_list,
                                           SEXP s_y,
                                           SEXP s_Ey,
                                           SEXP s_h_values,
                                           SEXP s_diffusion_params);
}


path_graph_plus_t create_path_graph_plus(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h);

path_graph_plm_t create_path_graph_plm(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h);

std::vector<double> graph_kpath_lm(const path_graph_plm_t& path_graph,
                                  const std::vector<double>& y,
                                  const std::vector<double>& weights, // this vector can either be empty or has to have the same number of elements as y
                                  int ikernel,
                                  double dist_normalization_factor = 1.01);

std::vector<double> graph_kmean(const std::vector<std::vector<int>>& graph,
                                               const std::vector<std::vector<double>>& edge_lengths,
                                               const std::vector<double>& y,
                                               int ikernel,
                                               double dist_normalization_factor = 1.01);

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
                             double dist_normalization_factor = 1.01,
                             int n_CVs = 0,
                             int n_CV_folds = 10,
                             double epsilon = 1e-10,
                             unsigned int seed = 0);

MS_complex_t graph_MS_cx(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<int>>& hop_list,
    const std::vector<std::vector<int>>& core_adj_list,
    const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
    const std::vector<double>& Ey);

/**
 * @brief Creates a sub path graph (plus version) for a smaller h value from an existing path graph with larger h
 *
 * @details Given P_{\bullet \leq h_max}(G) and h ≤ h_max, constructs P_{\bullet \leq h}(G)
 *          retaining the hop information
 *
 * @param source_graph The source path graph computed with a larger h value
 * @param h The desired (smaller) number of hops
 * @return path_graph_plus_t The path graph for the smaller h value, including hop counts
 *
 * @pre h must be less than or equal to the h value used to generate source_graph
 */
path_graph_plus_t create_sub_path_graph_plus(const path_graph_plus_t& source_graph, int h) {
    path_graph_plus_t result;
    const int n_vertices = source_graph.adj_list.size();
    result.adj_list.resize(n_vertices);
    result.weight_list.resize(n_vertices);
    result.hop_list.resize(n_vertices);

    for (int v = 0; v < n_vertices; ++v) {
        const auto& src_adj = source_graph.adj_list[v];
        const auto& src_weights = source_graph.weight_list[v];
        const auto& src_hops = source_graph.hop_list[v];

        // Pre-allocate with conservative estimate
        result.adj_list[v].reserve(src_adj.size());
        result.weight_list[v].reserve(src_adj.size());
        result.hop_list[v].reserve(src_adj.size());

        // Filter edges based on hop count
        for (size_t i = 0; i < src_adj.size(); ++i) {
            if (src_hops[i] <= h) {
                result.adj_list[v].push_back(src_adj[i]);
                result.weight_list[v].push_back(src_weights[i]);
                result.hop_list[v].push_back(src_hops[i]);

                // Only store path for smaller vertex index to larger
                if (v < src_adj[i]) {
                    result.shortest_paths[{v, src_adj[i]}] =
                        source_graph.shortest_paths.at({v, src_adj[i]});
                }
            }
        }
    }

    return result;
}

/**
 * @brief Generates a series of path graphs (plus version) for different h values efficiently
 *
 * @param adj_list The adjacency list of the input graph
 * @param weight_list The weight list of the input graph
 * @param h_values Vector of h values to generate graphs for
 * @return std::vector<path_graph_plus_t> Vector of path graphs for each requested h value
 *
 * @pre h_values must not be empty
 * @pre h_values must be sorted in ascending order for optimal performance
 */
std::vector<path_graph_plus_t> create_path_graph_series_plus(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<int>& h_values) {

    if (h_values.empty()) {
        Rf_error("h_values must not be empty");
    }

    // Use the maximum value in h_values as h_max
    const int h_max = *std::max_element(h_values.begin(), h_values.end());

    // First, create the graph for maximum h
    path_graph_plus_t max_graph = create_path_graph_plus(adj_list, weight_list, h_max);

    std::vector<path_graph_plus_t> result;
    result.reserve(h_values.size());

    // Generate graphs for each requested h value
    for (int h : h_values) {
        if (h == h_max) {
            result.push_back(max_graph);
        } else {
            result.push_back(create_sub_path_graph_plus(max_graph, h));
        }
    }

    return result;
}


/**
 * @brief Performs a sequence of graph analyses including path graph creation, kernel smoothing,
 *        diffusion smoothing, and Morse-Smale complex computation at multiple scales
 *
 * @param adj_list The adjacency list of the input graph
 * @param weight_list The weight list of the input graph
 * @param y Optional vector of vertex function values (can be empty if Ey is provided)
 * @param Ey Optional vector of pre-computed expected values (can be empty if y is provided)
 * @param h_values Vector of h values to use for path graph creation
 * @param diffusion_params Parameters for graph diffusion smoothing
 * @return std::vector<std::pair<MS_complex_t, std::vector<trajectory_t>>> Results for each h value
 *
 * @pre Exactly one of y or Ey must be non-empty
 * @pre h_values must not be empty
 */
struct diffusion_parameters_t {
    int n_time_steps;
    double step_factor;
    int normalize;
    bool preserve_local_maxima;
    double local_maximum_weight_factor;
    bool preserve_local_extrema;
    imputation_method_t imputation_method;
    iterative_imputation_params_t iterative_params;
    bool apply_binary_threshold;
    double binary_threshold;
    int ikernel;
    double dist_normalization_factor;
    int n_CVs;
    int n_CV_folds;
    double epsilon;
    unsigned int seed;
};

enum class cond_exp_estimator_t {
    MEAN,            // kernel mean over closed neighbor set (including the given vertex)
    PATH_MEAN,       // kernel mean over neighbor path means
    LM,              // linear model
};


std::vector<MS_complex_plus_t>
compute_graph_analysis_sequence(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    const std::vector<double>& Ey,
    const std::vector<int>& h_values,
    const diffusion_parameters_t& diffusion_params,
    cond_exp_estimator_t cond_exp_estimator = cond_exp_estimator_t::LM,
    bool use_diffn_smoother = false,
    bool verbose = false) {

    // Input validation
    if (y.empty() && Ey.empty()) {
        Rf_error("Either y or Ey must be provided");
    }
    if (!y.empty() && !Ey.empty()) {
        Rf_error("Only one of y or Ey should be provided");
    }
    if (h_values.empty()) {
        Rf_error("h_values cannot be empty");
    }

    // How to check if cond_exp_estimator has one of the allowed three values?


    auto ptm = std::chrono::steady_clock::now(); // for elapsed time
    std::vector<MS_complex_plus_t> results;
    std::vector<path_graph_plus_t> path_graphs;

    Rprintf("In compute_graph_analysis_sequence()\n");

    // Create path graphs based on h_values size
    if (h_values.size() == 1) { // Single h value case
        if (verbose) {
            Rprintf("Creating path graph ... ");

            ptm = std::chrono::steady_clock::now();
        }
        path_graph_plus_t pg_plus = create_path_graph_plus(adj_list, weight_list, h_values[0]);
        path_graphs.push_back(pg_plus);
        if (verbose) elapsed_time(ptm,"DONE",false);

    } else {
        // Multiple h values case
        if (verbose) {
            Rprintf("Creating path graphs ... ");
            ptm = std::chrono::steady_clock::now();
        }
        path_graphs = create_path_graph_series_plus(adj_list, weight_list, h_values);
        if (verbose) elapsed_time(ptm,"DONE",false);
    }

    std::vector<double> weights(adj_list.size(), 1.0); // Initialize weights - needed by graph_kpath_lm and graph_diffusion_smoother

    // Process each path graph
    for (size_t i = 0; i < path_graphs.size(); ++i) {
        const auto& pg = path_graphs[i];
        std::vector<double> current_Ey;
        std::vector<double> conditional_expectation; // estimated conditional expectation

        if (verbose) Rprintf("path graph: %d\n",(int)i);

        // Compute Ey if not provided
        if (Ey.empty()) {

            if (cond_exp_estimator == cond_exp_estimator_t::MEAN) {
                // Computing kernel weighted mean
                if (verbose) {
                    Rprintf("Computing kernel means ... ");
                    ptm = std::chrono::steady_clock::now();
                }
                conditional_expectation = graph_kmean(pg.adj_list,
                                                      pg.weight_list,
                                                      y,
                                                      diffusion_params.ikernel,
                                                      diffusion_params.dist_normalization_factor);
                if (verbose) elapsed_time(ptm,"DONE",false);

            } else if (cond_exp_estimator == cond_exp_estimator_t::PATH_MEAN) {

                if (verbose) {
                    Rprintf("Computing kernel path means ... ");
                    ptm = std::chrono::steady_clock::now();
                }
                // conditional_expectation = graph_kpath_mean(pg.adj_list,
                //                                            pg.weight_list,
                //                                            adj_list,
                //                                            y,
                //                                            diffusion_params.ikernel,
                //                                            diffusion_params.dist_normalization_factor);
                if (verbose) elapsed_time(ptm,"DONE",false);

            } else if (cond_exp_estimator == cond_exp_estimator_t::LM) {

                if (verbose) {
                    Rprintf("Computing local linear models ... ");
                    ptm = std::chrono::steady_clock::now();
                }
                #if 0
                conditional_expectation = graph_kpath_lm(llm_path_graph,
                                                         y,
                                                         weights,
                                                         diffusion_params.ikernel,
                                                         diffusion_params.dist_normalization_factor);
                #endif
                if (verbose) elapsed_time(ptm,"DONE",false);
            }


            if (use_diffn_smoother) {
                // Apply diffusion smoothing
                if (verbose) {
                    Rprintf("Computing diffusion smoothing of the kernel means ... ");
                    ptm = std::chrono::steady_clock::now();
                }

                auto diffusion_result = ext_graph_diffusion_smoother(
                    pg.adj_list,
                    pg.weight_list,
                    weights,
                    conditional_expectation,
                    diffusion_params.n_time_steps,
                    diffusion_params.step_factor,
                    diffusion_params.normalize,
                    diffusion_params.preserve_local_maxima,
                    diffusion_params.local_maximum_weight_factor,
                    diffusion_params.preserve_local_extrema,
                    diffusion_params.imputation_method,
                    diffusion_params.iterative_params,
                    diffusion_params.apply_binary_threshold,
                    diffusion_params.binary_threshold,
                    diffusion_params.ikernel,
                    diffusion_params.dist_normalization_factor,
                    diffusion_params.n_CVs,
                    diffusion_params.n_CV_folds,
                    diffusion_params.epsilon,
                    diffusion_params.seed
                    );
                if (verbose) elapsed_time(ptm,"DONE",false);

                // Check if CV was performed and we found an optimal solution
                if (!diffusion_result->mean_cv_errors.empty() && diffusion_result->optimal_time_step >= 0) {
                    current_Ey = diffusion_result->y_optimal; // The optimal smoothed values are directly available in y_optimal
                    // The corresponding Rf_error and time step
                    // double optimal_cv_error = result->min_cv_error;
                    // int optimal_time = result->optimal_time_step;
                } else {
                    current_Ey = diffusion_result->y_traj.back(); // Use the final smoothed values
                }
            } else {
                current_Ey = conditional_expectation;
            }

        } else {
            current_Ey = Ey;
        }

        // Compute MS complex
        Rprintf("Computing gradient trajectories complex ... ");
        auto ms_result = graph_MS_cx(pg.adj_list,
                                     pg.hop_list,
                                     adj_list, // core_adj_list is the input adj_list
                                     pg.shortest_paths,
                                     current_Ey);
        Rprintf("DONE\n");

        MS_complex_plus_t ms_cx_plus;
        // copy ms_result to ms_cx_plus
        ms_cx_plus.lmax_to_lmin        = ms_result.lmax_to_lmin;
        ms_cx_plus.lmin_to_lmax        = ms_result.lmin_to_lmax;
        ms_cx_plus.local_maxima        = ms_result.local_maxima;
        ms_cx_plus.local_minima        = ms_result.local_minima;
        ms_cx_plus.procells            = ms_result.procells;
        ms_cx_plus.cells               = ms_result.cells;
        ms_cx_plus.unique_trajectories = ms_result.unique_trajectories;
        ms_cx_plus.cell_trajectories   = ms_result.cell_trajectories;

        if (Ey.empty()) {
            ms_cx_plus.Ey = current_Ey;
        } else {
            ms_cx_plus.Ey = std::vector<double>(); // Initialize as empty vector
        }

        ms_cx_plus.h_value = h_values[i]; // Use the correct h-value from h_values vector
        ms_cx_plus.path_graph_adj_list    = pg.adj_list;
        ms_cx_plus.path_graph_weight_list = pg.weight_list; // Fixed member access
        ms_cx_plus.shortest_paths         = pg.shortest_paths; // Fixed member access

        if (Ey.empty()) {
            ms_cx_plus.Ey = current_Ey;  // Store computed Ey
        } else {
            ms_cx_plus.Ey = Ey;
        }

        results.push_back(ms_cx_plus);
    }

    return results;
}


/**
 * @brief R interface function for computing a sequence of graph analyses including path graph creation,
 *        kernel smoothing, diffusion smoothing, and Morse-Smale complex computation at multiple scales
 *
 * @param s_adj_list An R list representing the adjacency list of the input graph where each element
 *                   is an integer vector containing the indices of adjacent vertices
 *
 * @param s_weight_list An R list representing the weights of edges in the graph where each element
 *                      is a numeric vector containing the weights of edges to adjacent vertices
 *
 * @param s_y An R numeric vector representing function values over vertices, or R_NilValue if
 *            pre-computed expectation values are provided in s_Ey. Contains raw function values
 *            that will be smoothed using kernel weighted mean and diffusion smoothing
 *
 * @param s_Ey An R numeric vector representing pre-computed expected values over vertices,
 *             or R_NilValue if raw function values are provided in s_y. These values are
 *             used directly without additional smoothing
 *
 * @param s_h_values An R integer vector containing the h-values (hop distances) for which
 *                   to compute path graphs and subsequent analyses
 *
 * @param s_diffusion_params An R list containing parameters for diffusion smoothing with the following elements:
 *                          - n_time_steps: Integer, number of time steps for diffusion
 *                          - step_factor: Numeric, factor controlling step size
 *                          - normalize: Integer, normalization method
 *                          - preserve_local_maxima: Logical, whether to preserve local maxima
 *                          - local_maximum_weight_factor: Numeric, weight factor for local maxima
 *                          - preserve_local_extrema: Logical, whether to preserve local extrema
 *                          - imputation_method: Integer, method for imputing missing values
 *                          - iterative_params: List containing:
 *                              - max_iterations: Integer, maximum number of iterations
 *                              - convergence_threshold: Numeric, threshold for convergence
 *                          - apply_binary_threshold: Logical, whether to apply binary threshold
 *                          - binary_threshold: Numeric, threshold value
 *                          - ikernel: Integer, kernel type
 *                          - dist_normalization_factor: Numeric, factor for distance normalization
 *                          - n_CVs: Integer, number of cross-validations
 *                          - n_CV_folds: Integer, number of CV folds
 *                          - epsilon: Numeric, small value for numerical stability
 *                          - seed: Integer, random seed for reproducibility
 *
 * @return SEXP An R list of Morse-Smale complexes, one for each h-value. Each complex is represented
 *              as a list containing:
 *              - lmax_to_lmin: Named list mapping local maxima to sets of connected local minima
 *              - lmin_to_lmax: Named list mapping local minima to sets of connected local maxima
 *              - local_maxima: Integer vector of local maxima vertices
 *              - local_minima: Integer vector of local minima vertices
 *              - procells: Named list mapping (max,min) pairs to proto-cells
 *              - cells: Named list mapping (max,min) pairs to decomposed cells
 *              - unique_trajectories: List of integer vectors representing unique paths
 *              - cell_trajectories: Named list mapping cells to trajectory indices
 *              - path_graph_adj_list: Adjacency list of the path graph
 *              - path_graph_weight_list: Weight list of the path graph
 *              - shortest_paths: Named list mapping vertex pairs to shortest paths
 *              - Ey: Numeric vector of smoothed function values (empty if Ey was provided)
 *              Each complex also includes the h-value used for its computation
 *
 * @throws Rf_error if both s_y and s_Ey are R_NilValue
 *
 * @note Memory management is handled through PROTECT/UNPROTECT mechanism
 * @note The function performs extensive data structure conversions between R and C++
 * @note For large graphs or multiple h-values, this function can be computationally intensive
 * @note The path graphs generated for each h-value are retained in the output to allow reuse
 *
 * @example
 * # R usage example:
 * result <- compute_graph_analysis_sequence(
 *   adj_list = my_graph$adj_list,
 *   weight_list = my_graph$weights,
 *   y = vertex_function,
 *   Ey = NULL,
 *   h_values = c(1, 2, 3),
 *   diffusion_params = list(
 *     n_time_steps = 100,
 *     step_factor = 0.1,
 *     normalize = 1,
 *     # ... other parameters ...
 *   )
 * )
 */
SEXP S_compute_graph_analysis_sequence(SEXP s_adj_list,
                                       SEXP s_weight_list,
                                       SEXP s_y,
                                       SEXP s_Ey,
                                       SEXP s_h_values,
                                       SEXP s_diffusion_params) {
    // Converting R inputs to C++ types
    std::vector<std::vector<int>> adj_vect       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

    std::vector<double> y, Ey; // Changed to double from int
    if (s_y != R_NilValue) {
        PROTECT(s_y = Rf_coerceVector(s_y, REALSXP)); // Changed to REALSXP
        double* y_array = REAL(s_y);
        int y_len = LENGTH(s_y);
        y.assign(y_array, y_array + y_len);
        UNPROTECT(1);
    } else if (s_Ey != R_NilValue) {
        PROTECT(s_Ey = Rf_coerceVector(s_Ey, REALSXP)); // Changed to REALSXP
        double* Ey_array = REAL(s_Ey);
        int Ey_len = LENGTH(s_Ey);
        Ey.assign(Ey_array, Ey_array + Ey_len);
        UNPROTECT(1);
    } else {
        Rf_error("s_y and s_Ey cannot be both null.");
    }

    PROTECT(s_h_values = Rf_coerceVector(s_h_values, INTSXP));
    int* h_values_array = INTEGER(s_h_values);
    int h_values_len = LENGTH(s_h_values);
    std::vector<int> h_values(h_values_array, h_values_array + h_values_len);
    UNPROTECT(1);

    // Unpack diffusion parameters from R list
    diffusion_parameters_t diffusion_params;
    PROTECT(s_diffusion_params = Rf_coerceVector(s_diffusion_params, VECSXP));

    // Extract each parameter from the R list
    diffusion_params.n_time_steps = INTEGER(VECTOR_ELT(s_diffusion_params, 0))[0];
    diffusion_params.step_factor = REAL(VECTOR_ELT(s_diffusion_params, 1))[0];
    diffusion_params.normalize = INTEGER(VECTOR_ELT(s_diffusion_params, 2))[0];
    diffusion_params.preserve_local_maxima = LOGICAL(VECTOR_ELT(s_diffusion_params, 3))[0];
    diffusion_params.local_maximum_weight_factor = REAL(VECTOR_ELT(s_diffusion_params, 4))[0];
    diffusion_params.preserve_local_extrema = LOGICAL(VECTOR_ELT(s_diffusion_params, 5))[0];
    diffusion_params.imputation_method = static_cast<imputation_method_t>(INTEGER(VECTOR_ELT(s_diffusion_params, 6))[0]);

    // Get iterative params
    SEXP s_iterative_params = VECTOR_ELT(s_diffusion_params, 7);
    diffusion_params.iterative_params.max_iterations = INTEGER(VECTOR_ELT(s_iterative_params, 0))[0];
    diffusion_params.iterative_params.convergence_threshold = REAL(VECTOR_ELT(s_iterative_params, 1))[0];

    // Continue with remaining parameters
    diffusion_params.apply_binary_threshold = LOGICAL(VECTOR_ELT(s_diffusion_params, 8))[0];
    diffusion_params.binary_threshold = REAL(VECTOR_ELT(s_diffusion_params, 9))[0];
    diffusion_params.ikernel = INTEGER(VECTOR_ELT(s_diffusion_params, 10))[0];
    diffusion_params.dist_normalization_factor = REAL(VECTOR_ELT(s_diffusion_params, 11))[0];
    diffusion_params.n_CVs = INTEGER(VECTOR_ELT(s_diffusion_params, 12))[0];
    diffusion_params.n_CV_folds = INTEGER(VECTOR_ELT(s_diffusion_params, 13))[0];
    diffusion_params.epsilon = REAL(VECTOR_ELT(s_diffusion_params, 14))[0];
    diffusion_params.seed = INTEGER(VECTOR_ELT(s_diffusion_params, 15))[0];

    UNPROTECT(1);

    // Call the C++ function
    std::vector<MS_complex_plus_t> results = compute_graph_analysis_sequence(adj_vect,
                                                                             weight_vect,
                                                                             y,
                                                                             Ey,
                                                                             h_values,
                                                                             diffusion_params);

    // Convert results to R list
    SEXP r_results_list;
    PROTECT(r_results_list = Rf_allocVector(VECSXP, results.size()));

    for (size_t i = 0; i < results.size(); i++) {
        const auto& ms = results[i];

        // Create list for single MS complex
        SEXP ms_list;
        PROTECT(ms_list = Rf_allocVector(VECSXP, 12)); // Number of fields in MS_complex_plus_t

        // Convert each component
        SET_VECTOR_ELT(ms_list, 0, convert_set_to_R(ms.local_maxima));
        SET_VECTOR_ELT(ms_list, 1, convert_set_to_R(ms.local_minima));
        SET_VECTOR_ELT(ms_list, 2, convert_map_set_to_R(ms.lmax_to_lmin));
        SET_VECTOR_ELT(ms_list, 3, convert_map_set_to_R(ms.lmin_to_lmax));
        SET_VECTOR_ELT(ms_list, 4, convert_procells_to_R(ms.procells));
        SET_VECTOR_ELT(ms_list, 5, convert_map_vector_set_to_R(ms.cells));
        SET_VECTOR_ELT(ms_list, 6, convert_vector_vector_int_to_R(ms.unique_trajectories)); UNPROTECT(1);
        SET_VECTOR_ELT(ms_list, 7, convert_cell_trajectories_to_R(ms.cell_trajectories));
        SET_VECTOR_ELT(ms_list, 8, convert_vector_vector_int_to_R(ms.path_graph_adj_list)); UNPROTECT(1);
        SET_VECTOR_ELT(ms_list, 9, convert_vector_vector_double_to_R(ms.path_graph_weight_list)); UNPROTECT(1);
        SET_VECTOR_ELT(ms_list, 10, convert_map_vector_to_R(ms.shortest_paths));
        SET_VECTOR_ELT(ms_list, 11, convert_vector_double_to_R(ms.Ey)); UNPROTECT(1); // convert_vector_double_to_R(ms.Ey)

        // Set names
        SEXP names;
        PROTECT(names = Rf_allocVector(STRSXP, 12));
        SET_STRING_ELT(names, 0, Rf_mkChar("local_maxima"));
        SET_STRING_ELT(names, 1, Rf_mkChar("local_minima"));
        SET_STRING_ELT(names, 2, Rf_mkChar("lmax_to_lmin"));
        SET_STRING_ELT(names, 3, Rf_mkChar("lmin_to_lmax"));
        SET_STRING_ELT(names, 4, Rf_mkChar("procells"));
        SET_STRING_ELT(names, 5, Rf_mkChar("cells"));
        SET_STRING_ELT(names, 6, Rf_mkChar("unique_trajectories"));
        SET_STRING_ELT(names, 7, Rf_mkChar("cell_trajectories"));
        SET_STRING_ELT(names, 8, Rf_mkChar("path_graph_adj_list"));
        SET_STRING_ELT(names, 9, Rf_mkChar("path_graph_weight_list"));
        SET_STRING_ELT(names, 10, Rf_mkChar("shortest_paths"));
        SET_STRING_ELT(names, 11, Rf_mkChar("Ey"));

        Rf_setAttrib(ms_list, R_NamesSymbol, names);
        UNPROTECT(1); // names

        SET_VECTOR_ELT(r_results_list, i, ms_list);
        UNPROTECT(1); // ms_list
    }

    UNPROTECT(1); // r_results_list
    return r_results_list;
}
