#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

// #include <omp.h>
#include "omp_compat.h"

#include <future>
#include <mutex>
#include <execution>
#include <atomic>
#include <mutex>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <unordered_map>
#include <map>

#include "exec_policy.hpp"
#include "pgmalo.hpp"
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


extern "C" {

    SEXP S_pgmalo(
        SEXP neighbors_r,
        SEXP edge_lengths_r,
        SEXP y_r,
        SEXP y_true_r,
        SEXP use_median_r,
        SEXP h_min_r,
        SEXP h_max_r,
        SEXP p_r,
        SEXP n_bb_r,
        SEXP bb_max_distance_deviation_r,
        SEXP n_CVs_r,
        SEXP n_CV_folds_r,
        SEXP seed_r,
        SEXP kernel_type_r,
        SEXP dist_normalization_factor_r,
        SEXP epsilon_r,
        SEXP verbose_r);

    SEXP S_upgmalo(SEXP s_x,
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
                   SEXP s_dist_normalization_factor,
                   SEXP s_epsilon,
                   SEXP s_verbose);
}

bool validate_vertex_paths(const path_graph_plm_t& path_graph);

pgmalo_t pgmalo_mp(const std::vector<std::vector<int>>& neighbors,
                   const std::vector<std::vector<double>>& edge_lengths,
                   const std::vector<double>& y,
                   const std::vector<double>& y_true,
                   bool use_median = false,
                   int h_min = 3,
                   int h_max = 31,
                   double p = 0.95,
                   int n_bb = 500,
                   int bb_max_distance_deviation = 0,
                   int n_CVs = 0,
                   int n_CV_folds = 10,
                   unsigned int seed = 0,
                   int kernel_type = 7L,
                   double dist_normalization_factor = 1.01,
                   double epsilon = 1e-15,
                   bool verbose = false);

bb_cri_t pgmalo_bb_cri(const path_graph_plm_t& path_graph,
                     const std::vector<double>& y,
                     double p,
                     int n_bb,
                     int max_distance_deviation,
                     bool use_median,  // Added missing parameter
                     int ikernel,
                     double dist_normalization_factor,
                     double epsilon);

std::vector<double> pgmalo_cv(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    int n_CVs,
    int n_CV_folds,
    unsigned int seed,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon);

std::vector<double> pgmalo_cv_parallel(
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

path_graph_plm_t sexp_to_path_graph_plm(SEXP s_path_graph);

/**
 * @brief Path Graph Linear Model
 *
 * @details Fits a PGMALO by testing different neighborhood sizes (h values) and selecting
 * the optimal size based on cross-validation errors. For each h value, constructs a path
 * graph and computes both global and local predictions. Optionally computes bootstrap
 * credible intervals.
 *
 *  This function performs several steps:
 *    1. For each odd h in [h_min, h_max]:
 *       - Constructs h-hop neighborhood Path Linear Model (PLM) graphs
 *       - Computes leave-one-out cross-validation (LOOCV) errors
 *   2. Global optimization:
 *       - Finds optimal h that minimizes mean CV Rf_error
 *       - Computes conditional expectations using optimal h
 *   3. Local optimization:
 *       - For each vertex, finds optimal h minimizing its CV Rf_error
 *       - Computes vertex-specific conditional expectations
 *   4. Optional: Computes Bayesian bootstrap credible intervals
 *       - Estimates central tendency of conditional expectations
 *       - Computes credible interval bounds at specified probability level
 *
 * @param neighbors Adjacency lists for each vertex
 * @param edge_lengths Corresponding edge lengths for each adjacency
 * @param y Response variables for each vertex
 * @param y_true True values for Rf_error computation (optional)
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
pgmalo_t pgmalo(const std::vector<std::vector<int>>& neighbors,
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
                int ikernel = 7L,
                double dist_normalization_factor = 1.01,
                double epsilon = 1e-15,
                bool verbose = false) {

    auto total_ptm = std::chrono::steady_clock::now();
    auto ptm = std::chrono::steady_clock::now();  // Start timing

    const int n_vertices = static_cast<int>(y.size());
    const int n_h_values = (h_max - h_min) / 2 + 1;

    #define DEBUG__pgmalo 0
    #if DEBUG__pgmalo
    Rprintf("\nIn pgmalo()\n");
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
    std::unordered_map<int, std::vector<double>> predictions_map;
    predictions_map.reserve(n_h_values);

    // Initialize uniform weights
    std::vector<double> weights(n_vertices, 1.0);
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Process each h value
    for (int i = 0, h = h_min; h <= h_max; h += 2, i++) {
        results.h_values[i] = h;
        //int local_max_distance_deviation = (h - 1) / 2;
        int local_max_distance_deviation = bb_max_distance_deviation;
        //if (y_binary) local_max_distance_deviation = 0;

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

        auto [predictions, errors] = spgmalo(path_graph,
                                            y,
                                            weights,
                                            ikernel,
                                            local_max_distance_deviation,
                                            dist_normalization_factor,
                                            epsilon,
                                            verbose);

        if (n_CVs && y_binary) { // overwriting errors with CV errors

            errors = pgmalo_cv_parallel(path_graph,
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

        #if DEBUG__pgmalo
        print_vect(predictions,"predictions");
        print_vect(errors,"errors");
        #endif

        // Calculate mean Rf_error
        if (verbose) {
            Rprintf("\tCalculating mean Rf_error ... ");
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
    results.opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * results.opt_h_idx;

    #if 0
    // Compute local predictions
    if (verbose) {
        Rprintf("Computing local predictions ... ");
        ptm = std::chrono::steady_clock::now();
    }
    results.local_predictions.resize(n_vertices);
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
        results.local_predictions[i] = predictions_map[local_opt_h][i];
    }
    if (verbose) elapsed_time(ptm, "DONE");
    #endif

    // Store optimal results
    results.graph = std::move(results.graphs[results.opt_h_idx]);
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

        bb_cri_t bb_cri_results = pgmalo_bb_cri(results.graph,
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
 * @details This is a specialized version of pgmalo for univariate data.
 *          It constructs a chain graph from the ordered x values and applies path linear model estimation.
 *          The function performs several steps:
 *          1. Creates a chain graph from ordered x values
 *          2. For each odd h in [h_min, h_max]:
 *             - Constructs h-hop neighborhood Path Linear Model (PLM) graphs
 *             - Computes leave-one-out cross-validation (LOOCV) errors
 *          3. Global optimization:
 *             - Finds optimal h that minimizes mean CV Rf_error
 *             - Computes conditional expectations using optimal h
 *          4. Local optimization:
 *             - For each vertex, finds optimal h minimizing its CV Rf_error
 *             - Computes vertex-specific conditional expectations
 *          5. Optional: Computes Bayesian bootstrap credible intervals
 *
 * @param x Vector of ordered x values
 * @param y Observed y values corresponding to x
 * @param y_true True y values for Rf_error calculation (optional)
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
 *   - opt_h: Optimal h value minimizing mean CV Rf_error
 *   - opt_h_graph: Graph structure for optimal h
 *   - predictions: Global conditional expectations using optimal h
 *   - local_predictions: Vertex-specific conditional expectations using locally optimal h
 *   - bb_predictions: Bootstrap central tendency estimates (if n_bb > 0)
 *   - opt_ci_lower, opt_ci_upper: Credible interval bounds (if n_bb > 0)
 *
 * @see pgmalo
 * @see create_chain_graph
 *
 * @note This function creates a chain graph where each vertex is connected to its immediate neighbors
 *       in the ordered sequence of x values.
 */
pgmalo_t upgmalo(const std::vector<double>& x,
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

    return pgmalo(x_graph,
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
                  dist_normalization_factor,
                  epsilon,
                  verbose);
}

/**
 * @brief R interface for univariate path linear model optimization using k-path locally weighted models
 *
 * @details This function provides an R interface to upgmalo.
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
 * @param s_dist_normalization_factor R numeric for distance normalization
 * @param s_epsilon R numeric for numerical stability
 * @param s_verbose R logical for progress messages
 *
 * @return An R list containing:
 *   - h_values: Integer vector of used h values
 *   - opt_h: Optimal h value (numeric)
 *   - opt_predictions: Numeric vector of conditional expectation estimates
 *   - opt_local_predictions: Numeric vector of conditional expectation estimates using locally (vertex-wise) optimal h values
 *   - opt_bb_predictions: Numeric vector of bootstrap estimates (if n_bb > 0)
 *   - opt_ci_lower: Numeric vector of lower credible bounds (if n_bb > 0)
 *   - opt_ci_upper: Numeric vector of upper credible bounds (if n_bb > 0)
 *   - h_cv_errors: Numeric vector of cross-validation errors
 *   - true_error: Mean true Rf_error (if y_true provided)
 *   - opt_graph_adj_list: List of adjacency lists for optimal graph
 *   - opt_graph_edge_lengths: List of edge lengths for optimal graph
 *
 * @note This function uses R's memory protection mechanism via PROTECT/UNPROTECT
 *
 * @seealso upgmalo
 *
 * @example
 * \dontrun{
 * # R example:
 * result <- .Call("S_upgmalo",
 *                 x = as.numeric(1:100),
 *                 y = rnorm(100),
 *                 y_true = NULL,
 *                 use_median = FALSE,
 *                 h_min = 3L,
 *                 h_max = 31L,
 *                 p = 0.95,
 *                 n_bb = 500L,
 *                 ikernel = 1L,
 *                 dist_normalization_factor = 1.01,
 *                 epsilon = 1e-15,
 *                 verbose = TRUE)
 * }
 */
SEXP S_upgmalo(SEXP s_x,
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
               SEXP s_dist_normalization_factor,
               SEXP s_epsilon,
               SEXP s_verbose) {

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
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    auto cpp_results = upgmalo(x,
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
                               dist_normalization_factor,
                               epsilon,
                                       verbose);
    // Creating return list
    int n_protected = 0;  // Track number of PROTECT calls
    const int N_COMPONENTS = 13;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(cpp_results.h_values)); n_protected++;

    if (!cpp_results.h_cv_errors.empty() && n_points > 0) {
        SET_VECTOR_ELT(result, 1, convert_vector_double_to_R(cpp_results.h_cv_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    SEXP s_opt_h_idx = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h_idx)[0] = cpp_results.opt_h_idx + 1;
    SET_VECTOR_ELT(result, 2, s_opt_h_idx);

    SEXP s_opt_h = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
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
        SEXP s_true_error = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
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
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, Rf_mkChar("h_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("h_errors"));
    SET_STRING_ELT(names, 2, Rf_mkChar("opt_h_idx"));
    SET_STRING_ELT(names, 3, Rf_mkChar("opt_h"));
    SET_STRING_ELT(names, 4, Rf_mkChar("graph_adj_list"));
    SET_STRING_ELT(names, 5, Rf_mkChar("graph_edge_lengths"));
    SET_STRING_ELT(names, 6, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 7, Rf_mkChar("local_predictions"));
    SET_STRING_ELT(names, 8, Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 9, Rf_mkChar("ci_lower"));
    SET_STRING_ELT(names, 10, Rf_mkChar("ci_upper"));
    SET_STRING_ELT(names, 11, Rf_mkChar("true_error"));
    SET_STRING_ELT(names, 12, Rf_mkChar("h_predictions"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}


/**
 * @brief Implements Single Path Graph Linear Model
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
 * @Rf_warning
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
 * auto [predictions, errors] = spgmalo(
 *     graph, y, w, kernel_type, max_dev);
 * @endcode
 */
std::pair<std::vector<double>, std::vector<double>> spgmalo(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int kernel_type,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon,
    bool verbose) {

    #define DEBUG_single_path_graph_malm 0

    int h = path_graph.h; // must be even and such that h >= 4 and h <= n_vertices - 2, where n_vertices is the number of vertices of the graph
    int mid_vertex = h / 2;
    int path_n_vertices = h + 1;
    int n_vertices = path_graph.vertex_paths.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    std::vector<double> predictions(n_vertices);
    std::vector<double> errors(n_vertices);

    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("SPGMALO");

    if (verbose) {
        Rprintf("Starting SPGMALO computation\n");
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

    //------------------------------------------------------------------------------
    // Phase 1: Single Model Computation
    //------------------------------------------------------------------------------

    //     2. Compute kernel weights based on distances
    //     3. Fit local linear model using weighted least squares
    //     4. Store predictions, weights, and errors for each point in window


    //auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("  Phase 1: Computing single-model predictions ... ");
    }
    auto phase1_ptm = std::chrono::steady_clock::now();

    struct pred_w_err_t {
        double prediction;
        double weight;
        double Rf_error;
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

#if DEBUG_single_path_graph_malm
            Rprintf("path_i: %d\n", (int)path_i);
#endif

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

            // Fitting a weighted linear model and return as well model's LOOCV prediction errors
            // ulm_t constains predictions and errors
            ulm_t fit = ulm(x_path.data(),
                            y_path.data(),
                            w_path,
                            y_binary,
                            epsilon);

            // For each vertex of the path record predicted value, the weight, and the models LOOCV at that vertex
            for (int i = 0; i < path_n_vertices; ++i) {
                pred_w_err.prediction = fit.predictions[i];
                pred_w_err.weight     = w_path[i];
                pred_w_err.Rf_error      = fit.errors[i];
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
            wmean_error  += v.weight * v.Rf_error;
        }
        predictions[i] = weighted_sum / weight_sum;
        errors[i] = wmean_error / weight_sum;
    }

    if (verbose) {
        elapsed_time(phase2_ptm, "Done");
        mem_tracker.report();
    }

    if (verbose) {
        elapsed_time(total_ptm, "\nTotal SPGMALO computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return std::make_pair(predictions, errors);
}


/**
 * @brief R interface for Path Graph Linear Model Analysis
 *
 * @details An R interface wrapper for the pgmalo() function that handles data conversion between
 * R and C++ and returns results as an R list. This function processes graph data and fits a PGMALO
 * model by testing different neighborhood sizes, performing cross-validation, and optionally
 * computing bootstrap credible intervals.
 *
 * The function converts R input objects to C++ data structures, calls the core pgmalo() implementation,
 * and converts the results back to R objects. All memory management follows R's protection mechanisms
 * to ensure proper garbage collection.
 *
 * @param neighbors_r SEXP containing adjacency lists for each vertex
 * @param edge_lengths_r SEXP containing corresponding edge lengths for each adjacency
 * @param y_r SEXP containing response variables for each vertex
 * @param y_true_r SEXP containing true values for Rf_error computation (optional)
 * @param use_median_r SEXP logical for using median instead of mean for bootstrap intervals
 * @param h_min_r SEXP minimum neighborhood size (must be odd)
 * @param h_max_r SEXP maximum neighborhood size (must be odd)
 * @param p_r SEXP confidence level for bootstrap intervals (0 < p < 1)
 * @param n_bb_r SEXP number of bootstrap iterations
 * @param bb_max_distance_deviation_r SEXP maximal distance deviation for bootstrap samples
 * @param n_CVs_r SEXP number of cross-validations
 * @param n_CV_folds_r SEXP number of cross-validation folds
 * @param seed_r SEXP random seed for reproducibility
 * @param kernel_type_r SEXP kernel function identifier
 * @param n_cores_r SEXP number of cores for parallel computation
 * @param dist_normalization_factor_r SEXP distance normalization factor
 * @param epsilon_r SEXP numerical stability threshold
 * @param verbose_r SEXP enable progress reporting
 *
 * @return SEXP (R list) containing:
 *         - h_values: Integer vector of tested neighborhood sizes
 *         - opt_h: Optimal neighborhood size
 *         - opt_h_idx: Index of optimal neighborhood size
 *         - h_cv_errors: Vector of cross-validation errors for each h
 *         - true_errors: Vector of true errors (if y_true provided)
 *         - predictions: Vector of global predictions
 *         - local_predictions: Vector of local predictions
 *         - h_predictions: List of predictions for each h value
 *         - bb_predictions: Bootstrap predictions
 *         - ci_lower: Lower confidence interval bounds
 *         - ci_upper: Upper confidence interval bounds
 *         - has_bootstrap: Logical indicating if bootstrap results exist
 */
SEXP S_pgmalo(
    SEXP neighbors_r,
    SEXP edge_lengths_r,
    SEXP y_r,
    SEXP y_true_r,
    SEXP use_median_r,
    SEXP h_min_r,
    SEXP h_max_r,
    SEXP p_r,
    SEXP n_bb_r,
    SEXP bb_max_distance_deviation_r,
    SEXP n_CVs_r,
    SEXP n_CV_folds_r,
    SEXP seed_r,
    SEXP kernel_type_r,
    SEXP dist_normalization_factor_r,
    SEXP epsilon_r,
    SEXP verbose_r
    ) {

    std::vector<std::vector<int>> neighbors       = convert_adj_list_from_R(neighbors_r);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(edge_lengths_r);

    int n_vertices = LENGTH(y_r);
    std::vector<double> y(REAL(y_r), REAL(y_r) + n_vertices);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(y_true_r) == n_vertices) {
        y_true.assign(REAL(y_true_r), REAL(y_true_r) + n_vertices);
    }

    //int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    bool use_median = LOGICAL(use_median_r)[0];
    int h_min = INTEGER(h_min_r)[0];
    int h_max = INTEGER(h_max_r)[0];
    double p = REAL(p_r)[0];
    int n_bb = INTEGER(n_bb_r)[0];
    int bb_max_distance_deviation = INTEGER(bb_max_distance_deviation_r)[0];
    int n_CVs = INTEGER(n_CVs_r)[0];
    int n_CV_folds = INTEGER(n_CV_folds_r)[0];
    unsigned int seed = (unsigned int)INTEGER(seed_r)[0];
    int kernel_type = INTEGER(kernel_type_r)[0];
    double dist_normalization_factor = REAL(dist_normalization_factor_r)[0];
    double epsilon = REAL(epsilon_r)[0];
    bool verbose = LOGICAL(verbose_r)[0];

#if 0
    auto results = pgmalo_mp(neighbors,
                             edge_lengths,
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
                             kernel_type,
                             dist_normalization_factor,
                             epsilon,
                             verbose);
#else
    auto results = pgmalo(neighbors,
                          edge_lengths,
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
                          kernel_type,
                          dist_normalization_factor,
                          epsilon,
                          verbose);
#endif



    // Creating return list
    const int N_COMPONENTS = 11;
    int n_protected = 0;  // Track number of PROTECT calls
    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Convert and set h_values
    SET_VECTOR_ELT(result_r, 0, convert_vector_int_to_R(results.h_values)); n_protected++;

    // Set opt_h
    SEXP opt_h_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(opt_h_r)[0] = results.opt_h;
    SET_VECTOR_ELT(result_r, 1, opt_h_r);

    // Set opt_h_idx
    SEXP opt_h_idx_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(opt_h_idx_r)[0] = results.opt_h_idx;
    SET_VECTOR_ELT(result_r, 2, opt_h_idx_r);

    // Convert and set h_cv_errors
    SET_VECTOR_ELT(result_r, 3, convert_vector_double_to_R(results.h_cv_errors)); n_protected++;

    // Convert and set true_errors if they exist
    if (results.has_true_errors()) {
        SET_VECTOR_ELT(result_r, 4, convert_vector_double_to_R(results.true_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result_r, 4, R_NilValue);
    }

    // Convert and set predictions
    SET_VECTOR_ELT(result_r, 5, convert_vector_double_to_R(results.predictions)); n_protected++;

    // Convert and set local_predictions
    SET_VECTOR_ELT(result_r, 6, convert_vector_double_to_R(results.local_predictions)); n_protected++;

    // Convert and set h_predictions
    SEXP h_predictions_r = PROTECT(Rf_allocVector(VECSXP, results.h_predictions.size())); n_protected++;
    for (size_t i = 0; i < results.h_predictions.size(); i++) {
        SET_VECTOR_ELT(h_predictions_r, i,
                      convert_vector_double_to_R(results.h_predictions[i]));
        n_protected++;
    }
    SET_VECTOR_ELT(result_r, 7, h_predictions_r);

    // Set bootstrap results if they exist
    if (results.has_bootstrap_results()) {
        SET_VECTOR_ELT(result_r, 8, convert_vector_double_to_R(results.bb_predictions)); n_protected++;
        SET_VECTOR_ELT(result_r, 9, convert_vector_double_to_R(results.ci_lower)); n_protected++;
        SET_VECTOR_ELT(result_r, 10,convert_vector_double_to_R(results.ci_upper)); n_protected++;
    } else {
        SET_VECTOR_ELT(result_r, 8, R_NilValue);
        SET_VECTOR_ELT(result_r, 9, R_NilValue);
        SET_VECTOR_ELT(result_r, 10, R_NilValue);
    }

    // Set names for the list elements
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    const char* names_str[] = {
        "h_values",
        "opt_h",
        "opt_h_idx",
        "h_errors",
        "true_errors",
        "predictions",
        "local_predictions",
        "h_predictions",
        "bb_predictions",
        "ci_lower",
        "ci_upper"
    };
    for (int i = 0; i < N_COMPONENTS; i++) {
        SET_STRING_ELT(names, i, Rf_mkChar(names_str[i]));
    }
    Rf_setAttrib(result_r, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return result_r;
}

// Mutex for protecting console output
std::mutex console_mutex;

pgmalo_t pgmalo_mp(const std::vector<std::vector<int>>& neighbors,
                   const std::vector<std::vector<double>>& edge_lengths,
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
                   int kernel_type,
                   double dist_normalization_factor,
                   double epsilon,
                   bool verbose) {

    auto total_ptm = std::chrono::steady_clock::now();
    auto ptm = std::chrono::steady_clock::now();

    const int n_vertices = static_cast<int>(y.size());
    const int n_h_values = (h_max - h_min) / 2 + 1;

    // Initialize results structure
    pgmalo_t results;
    results.graphs.resize(n_h_values);
    results.h_cv_errors.resize(n_h_values);
    results.h_values.resize(n_h_values);
    results.h_predictions.resize(n_h_values);

    // Create vector of h values for parallel processing
    std::vector<int> h_values(n_h_values);
    std::generate(h_values.begin(), h_values.end(),
                 [h = h_min]() mutable { auto current = h; h += 2; return current; });

    // Thread-safe storage for intermediate results
    struct thread_safe_results_t {
        std::mutex mutex;
        std::unordered_map<int, std::vector<double>> predictions_map;

        void store_predictions(int h, std::vector<double>&& predictions) {
            std::lock_guard<std::mutex> lock(mutex);
            predictions_map[h] = std::move(predictions);
        }
    };

    thread_safe_results_t thread_safe_results;
    thread_safe_results.predictions_map.reserve(n_h_values);

    // Initialize uniform weights
    std::vector<double> weights(n_vertices, 1.0);
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Parallel processing of h values

    struct termination_state_t {
        std::mutex mutex;
        int failed_h = std::numeric_limits<int>::max();
        bool has_failed = false;
        std::vector<int> successful_h_values;

        bool should_terminate(int current_h) {
            std::lock_guard<std::mutex> lock(mutex);
            return has_failed && current_h >= failed_h;
        }

        void signal_failure(int h) {
            std::lock_guard<std::mutex> lock(mutex);
            if (!has_failed || h < failed_h) {
                has_failed = true;
                failed_h = h;
            }
        }

        void record_success(int h) {
            std::lock_guard<std::mutex> lock(mutex);
            successful_h_values.push_back(h);
        }
    } termination_state;

    std::for_each(GFLOW_EXEC_POLICY, h_values.begin(), h_values.end(),
        [&](int h) {
            int i = (h - h_min) / 2;  // Calculate index based on h value

            if (verbose) {
                {
                    std::lock_guard<std::mutex> lock(console_mutex);
                    Rprintf("\nProcessing h = %d (%d/%d)\n\tCreating path graph ... ",
                            h, i + 1, n_h_values);
                    ptm = std::chrono::steady_clock::now();
                }
            }

            auto path_graph = create_path_graph_plm(neighbors, edge_lengths, h);

            if (!validate_vertex_paths(path_graph)) {
                termination_state.signal_failure(h);
                return;
            }

            if (verbose) {
                {
                    std::lock_guard<std::mutex> lock(console_mutex);
                    elapsed_time(ptm, "DONE");
                }
            }

            // Compute predictions and errors
            if (verbose) {
                {
                    std::lock_guard<std::mutex> lock(console_mutex);
                    Rprintf("\tComputing predictions and errors ... ");
                    ptm = std::chrono::steady_clock::now();
                }
            }

            auto [predictions, errors] = spgmalo(path_graph,
                                               y,
                                               weights,
                                               kernel_type,
                                               bb_max_distance_deviation,
                                               dist_normalization_factor,
                                               epsilon,
                                               verbose);

            if (n_CVs && y_binary) {
                errors = pgmalo_cv_parallel(path_graph,
                                          y,
                                          n_CVs,
                                          n_CV_folds,
                                          seed,
                                          kernel_type,
                                          bb_max_distance_deviation,
                                          dist_normalization_factor,
                                          epsilon);

                errors.erase(std::remove_if(errors.begin(), errors.end(),
                                          [](double x) { return std::isnan(x); }),
                           errors.end());
            }

            if (verbose) {
                {
                    std::lock_guard<std::mutex> lock(console_mutex);
                    elapsed_time(ptm, "DONE");
                    Rprintf("\tCalculating mean Rf_error ... ");
                    ptm = std::chrono::steady_clock::now();
                }
            }

            const double mean_error = std::accumulate(errors.begin(), errors.end(), 0.0) / n_vertices;

            // Thread-safe updates to shared data structures
            {
                std::lock_guard<std::mutex> lock(thread_safe_results.mutex);
                results.h_cv_errors[i] = mean_error;
                results.graphs[i] = std::move(path_graph);
                results.h_values[i] = h;
                results.h_predictions[i] = predictions;
                thread_safe_results.predictions_map[h] = std::move(predictions);
            }

            if (verbose) {
                {
                    std::lock_guard<std::mutex> lock(console_mutex);
                    elapsed_time(ptm, "DONE");
                }
            }

            termination_state.record_success(h);
        });


    if (termination_state.has_failed) {
        // Sort successful h values to ensure they're in order
        std::sort(termination_state.successful_h_values.begin(),
                  termination_state.successful_h_values.end());

        // Resize result vectors to match successful h values
        size_t new_size = termination_state.successful_h_values.size();
        std::vector<path_graph_plm_t> new_graphs;
        std::vector<double> new_cv_errors;
        std::vector<int> new_h_values;
        std::vector<std::vector<double>> new_predictions;

        new_graphs.reserve(new_size);
        new_cv_errors.reserve(new_size);
        new_h_values.reserve(new_size);
        new_predictions.reserve(new_size);

        // Copy only successful results
        for (int h : termination_state.successful_h_values) {
            int i = (h - h_min) / 2;
            new_graphs.push_back(std::move(results.graphs[i]));
            new_cv_errors.push_back(results.h_cv_errors[i]);
            new_h_values.push_back(results.h_values[i]);
            new_predictions.push_back(std::move(results.h_predictions[i]));
        }

        // Update results structure with only successful data
        results.graphs = std::move(new_graphs);
        results.h_cv_errors = std::move(new_cv_errors);
        results.h_values = std::move(new_h_values);
        results.h_predictions = std::move(new_predictions);

        // Print Rf_warning about adjusted h range
        Rprintf("\nWarning: Path validation failed at h = %d. Analysis proceeding with ",
                termination_state.failed_h);
        if (results.h_values.empty()) {
            Rprintf("no valid h values.\n");
            // Handle the case where no h values were successful
            REPORT_ERROR("No valid h values found.");
        } else {
            Rprintf("h values in range [%d, %d].\n",
                    results.h_values.front(),
                    results.h_values.back());
        }
    }

    // Find optimal h
    auto min_it = std::min_element(results.h_cv_errors.begin(), results.h_cv_errors.end());
    results.opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * results.opt_h_idx;

    // Store optimal results
    results.graph = std::move(results.graphs[results.opt_h_idx]);
    results.predictions = std::move(thread_safe_results.predictions_map[results.opt_h]);

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

        bb_cri_t bb_cri_results = pgmalo_bb_cri(results.graph,
                                                y,
                                                p,
                                                n_bb,
                                                local_max_distance_deviation,
                                                use_median,
                                                kernel_type,
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
