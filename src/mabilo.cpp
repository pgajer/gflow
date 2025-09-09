#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

#include <execution>
#include <atomic>
#include <mutex>
#include <numeric>
#include <vector>
#include <random>
#include <algorithm>     // For std::max
#include <cmath>         // For std::fabs()

#include "exec_policy.hpp"
#include "mabilo.hpp"    // For mabilo_t
#include "sampling.h"    // For C_runif_simplex()
#include "cpp_utils.hpp" // For elapsed_time
#include "error_utils.h" // For REPORT_ERROR()
#include "kernels.h"     // For initialize_kernel()
#include "ulm.hpp"       // For ulm() utilities
#include "memory_utils.hpp"
#include "progress_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include "predictive_errors.hpp"

extern "C" {
    SEXP S_wmabilo(SEXP s_x,
                   SEXP s_y,
                   SEXP s_y_true,
                   SEXP s_w,
                   SEXP s_k_min,
                   SEXP s_k_max,
                   SEXP s_distance_kernel,
                   SEXP s_dist_normalization_factor,
                   SEXP s_epsilon,
                   SEXP s_verbose);

    SEXP S_mabilo(SEXP s_x,
                  SEXP s_y,
                  SEXP s_y_true,
                  SEXP s_k_min,
                  SEXP s_k_max,
                  SEXP s_n_bb,
                  SEXP s_p,
                  SEXP s_distance_kernel,
                  SEXP s_dist_normalization_factor,
                  SEXP s_epsilon,
                  SEXP s_verbose);

    SEXP S_mabilo_with_smoothed_errors(SEXP s_x,
                                       SEXP s_y,
                                       SEXP s_y_true,
                                       SEXP s_k_min,
                                       SEXP s_k_max,
                                       SEXP s_error_window_factor,
                                       SEXP s_distance_kernel,
                                       SEXP s_dist_normalization_factor,
                                       SEXP s_epsilon,
                                       SEXP s_verbose);
}

/**
 * @brief R interface for MABILO (Model-Averaged Locally Weighted Scatterplot Smoothing)
 *
 * @param s_x Vector of x coordinates (must be sorted)
 * @param s_y Vector of y coordinates (response values)
 * @param s_y_true Optional vector of true y values for error calculation
 * @param s_w Vector of weights for each point
 * @param s_k_min Minimum number of neighbors (must be positive)
 * @param s_k_max Maximum number of neighbors (must be greater than k_min)
 * @param s_distance_kernel Kernel type for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param s_dist_normalization_factor Factor for normalizing distances
 * @param s_epsilon Small number for numerical stability
 * @param s_verbose Whether to print progress information
 *
 * @return A list containing:
 * - k_values: Vector of k values tested
 * - opt_k: Optimal k value for model-averaged predictions
 * - opt_k_idx: Index of optimal k value
 * - k_mean_errors: Mean LOOCV errors for each k
 * - k_mean_true_errors: Mean true errors if y_true provided
 * - predictions: Model-averaged predictions using optimal k
 * - k_predictions: Model-averaged predictions for all k values
 *
 * @throws error if input vectors have inconsistent lengths or invalid parameters
 */
SEXP S_wmabilo(SEXP s_x,
                 SEXP s_y,
                 SEXP s_y_true,
                 SEXP s_w,
                 SEXP s_k_min,
                 SEXP s_k_max,
                 SEXP s_distance_kernel,
                 SEXP s_dist_normalization_factor,
                 SEXP s_epsilon,
                 SEXP s_verbose) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    bool y_true_exists = LENGTH(s_y_true) == n_points;
    if (y_true_exists) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    std::vector<double> w(REAL(s_w), REAL(s_w) + n_points);

    //int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    int k_min = INTEGER(s_k_min)[0];
    int k_max = INTEGER(s_k_max)[0];

    int distance_kernel = INTEGER(s_distance_kernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    mabilo_t wmabilo_results = wmabilo(x,
                                       y,
                                       y_true,
                                       w,
                                       k_min,
                                       k_max,
                                       distance_kernel,
                                       dist_normalization_factor,
                                       epsilon,
                                       verbose);

    // Creating return list
    const int N_COMPONENTS = 7;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    std::vector<int> k_values(wmabilo_results.k_mean_errors.size());
    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
        k_values[k_index] = k;
    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(k_values)); n_protected++;

    SEXP s_opt_k = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k)[0] = wmabilo_results.opt_k;
    SET_VECTOR_ELT(result, 1, s_opt_k);

    SEXP s_opt_k_idx = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k_idx)[0] = wmabilo_results.opt_k_idx;
    SET_VECTOR_ELT(result, 2, s_opt_k_idx);

    SET_VECTOR_ELT(result, 3, convert_vector_double_to_R(wmabilo_results.k_mean_errors)); n_protected++;

    // true errors
    if (y_true_exists) {
        SET_VECTOR_ELT(result, 4, convert_vector_double_to_R(wmabilo_results.k_mean_true_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 4, R_NilValue);
    }

    SET_VECTOR_ELT(result, 5, convert_vector_double_to_R(wmabilo_results.predictions)); n_protected++;
    SET_VECTOR_ELT(result, 6, convert_vector_vector_double_to_R(wmabilo_results.k_predictions)); n_protected++;

    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("k_values"));
    SET_STRING_ELT(names, 1, mkChar("opt_k"));
    SET_STRING_ELT(names, 2, mkChar("opt_k_idx"));
    SET_STRING_ELT(names, 3, mkChar("k_mean_errors"));
    SET_STRING_ELT(names, 4, mkChar("k_mean_true_errors"));
    SET_STRING_ELT(names, 5, mkChar("predictions"));
    SET_STRING_ELT(names, 6, mkChar("k_predictions"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}

/**
 * @brief Performs parallel Bayesian bootstrap calculations for MABILO
 *
 * @details Implements parallel bootstrap sampling for uncertainty quantification:
 * 1. Generates bootstrap weights using Bayesian bootstrap
 * 2. Performs MABILO fitting for each bootstrap sample
 * 3. Aggregates predictions across bootstrap iterations
 *
 * Thread safety is ensured through mutex-protected random number generation.
 *
 * @param x Vector of sorted x coordinates
 * @param y Vector of observed y values
 * @param k Number of neighbors for local fitting
 * @param n_bb Number of bootstrap iterations (must be positive)
 * @param distance_kernel Kernel function type (0: Tricube, 1: Epanechnikov, 2: Exponential)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-8)
 * @param verbose Enable progress messages (default: false)
 *
 * @return Vector of vectors containing bootstrap predictions for each iteration
 *
 * @throws Rf_error if parameters are invalid or computation fails
 * @note Thread-safe implementation using parallel execution
 */
std::vector<std::vector<double>> mabilo_bb(const std::vector<double>& x,
                                             const std::vector<double>& y,
                                             int k,
                                             int n_bb,
                                             int distance_kernel,
                                             double dist_normalization_factor = 1.01,
                                             double epsilon = 1e-8,
                                             bool verbose = false) {

    int n_points = static_cast<int>(y.size());

    // Initialize results vector
    std::vector<std::vector<double>> bb_predictions(n_bb);
    for (auto& Ey : bb_predictions) {
        Ey.resize(n_points);
    }

    // Create indices for parallel iteration
    std::vector<int> bb_indices(n_bb);
    std::iota(bb_indices.begin(), bb_indices.end(), 0);

    // Mutex for thread-safe random number generation
    std::mutex rng_mutex;

    // Parallel execution of bootstrap iterations
    //std::for_each(std::execution::par_unseq,
    std::for_each(GFLOW_EXEC_POLICY,
                  bb_indices.begin(),
                  bb_indices.end(),
                  [&](int iboot) {
        // Thread-local weight vector
        std::vector<double> weights(n_points);

        // Generate weights in a thread-safe manner
        {
            std::lock_guard<std::mutex> lock(rng_mutex);
            C_runif_simplex(&n_points, weights.data());
        }

        // Compute predictions for this bootstrap iteration
        std::vector<double> y_true;
        auto wmabilo_results = wmabilo(x,
                                       y,
                                       y_true,
                                       weights,
                                       k,
                                       k,
                                       distance_kernel,
                                       dist_normalization_factor,
                                       epsilon,
                                       verbose);

        // Store results - no need for mutex as each thread writes to its own index
        bb_predictions[iboot] = std::move(wmabilo_results.predictions);
    });

    return bb_predictions;
}

/**
 * @brief Computes Bayesian bootstrap predictions with credible intervals for MABILO
 *
 * @details Implements a complete Bayesian uncertainty analysis:
 * 1. Performs parallel bootstrap iterations using mabilo_bb
 * 2. Computes point estimates using median
 * 3. Calculates credible intervals at specified probability level
 *
 * @param x Vector of sorted x coordinates
 * @param y Vector of observed y values
 * @param k Number of neighbors for local fitting
 * @param n_bb Number of bootstrap iterations (must be positive)
 * @param p Probability level for credible intervals (must be in (0,1))
 * @param distance_kernel Kernel function type
 * @param dist_normalization_factor Distance normalization factor
 * @param epsilon Numerical stability parameter
 *
 * @return bb_cri_t structure containing:
 *         - bb_predictions: Median predictions across bootstrap iterations
 *         - cri_L: Lower bounds of credible intervals
 *         - cri_U: Upper bounds of credible intervals
 *
 * @throws Rf_error for invalid parameters or failed computation
 */
bb_cri_t mabilo_bb_cri(const std::vector<double>& x,
                         const std::vector<double>& y,
                         int k,
                         int n_bb,
                         double p,
                         int distance_kernel,
                         double dist_normalization_factor,
                         double epsilon) {

    // Perform bootstrap iterations
    std::vector<std::vector<double>> bb_predictionss = mabilo_bb(x,
                                                                   y,
                                                                   k,
                                                                   n_bb,
                                                                   distance_kernel,
                                                                   dist_normalization_factor,
                                                                   epsilon);

    // Calculate credible intervals
    bool use_median = true;
    return bb_cri(bb_predictionss, use_median, p);
}


/**
 * @brief Main interface for MABILO (Model-Averaged Locally Weighted Scatterplot Smoothing) with optional Bayesian bootstrap
 *
 * @details This function provides a complete implementation of MABILO with two main components:
 * 1. Core MABILO algorithm:
 *    - Fits local linear models using k-hop neighborhoods
 *    - Performs kernel-weighted model averaging
 *    - Finds optimal window size k through LOOCV
 *
 * 2. Optional Bayesian bootstrap analysis:
 *    - Computes bootstrap predictions using optimal k
 *    - Calculates credible intervals for uncertainty quantification
 *    - Provides central location estimates
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Vector of y values (response variable)
 * @param y_true Optional vector of true y values for error calculation
 * @param k_min Minimum number of neighbors on each side
 * @param k_max Maximum number of neighbors on each side
 * @param n_bb Number of Bayesian bootstrap iterations (0 to skip bootstrap)
 * @param p Probability level for credible intervals (used only if n_bb > 0)
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilo_t structure containing:
 *         - opt_k: Optimal k value for model-averaged predictions
 *         - predictions: Model-averaged predictions using optimal k
 *         - k_mean_errors: Mean LOOCV errors for each k
 *         - k_mean_true_errors: Mean true errors if y_true provided
 *         - k_predictions: Predictions for all k values
 *         If n_bb > 0, also includes:
 *         - bb_predictions: Central location of bootstrap estimates
 *         - cri_L: Lower bounds of credible intervals
 *         - cri_U: Upper bounds of credible intervals
 *
 * @throws std::invalid_argument if:
 *         - Input vectors have inconsistent lengths
 *         - k_min or k_max values are invalid
 *         - n_bb is negative
 *         - p is not in (0,1) when n_bb > 0
 *
 * @note
 * - Input x values must be sorted in ascending order
 * - Window size for each k is 2k + 1 (k points on each side plus center)
 * - Uses equal weights (1.0) for all observations in core algorithm
 * - Bootstrap analysis uses optimal k from core algorithm
 */
mabilo_t mabilo(const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& y_true,
                    int k_min,
                    int k_max,
                    int n_bb,
                    double p,
                    int distance_kernel,
                    double dist_normalization_factor,
                    double epsilon,
                    bool verbose) {

    mabilo_t uwmabilo_results = uwmabilo(x,
                                         y,
                                         y_true,
                                         k_min,
                                         k_max,
                                         distance_kernel,
                                         dist_normalization_factor,
                                         epsilon,
                                         verbose);

    if (n_bb) {
        bb_cri_t bb_res = mabilo_bb_cri(x,
                                        y,
                                        uwmabilo_results.opt_k,
                                        n_bb,
                                        p,
                                        distance_kernel,
                                        dist_normalization_factor,
                                        epsilon);
        uwmabilo_results.bb_predictions = std::move(bb_res.bb_Ey);
        uwmabilo_results.cri_L          = std::move(bb_res.cri_L);
        uwmabilo_results.cri_U          = std::move(bb_res.cri_U);
    }

    return uwmabilo_results;
}

/**
 * @brief R interface for MABILO with Bayesian bootstrap capability
 *
 * @details Provides an R interface to the MABILO algorithm with optional Bayesian bootstrap analysis.
 * Converts R objects to C++ types, calls the core implementation, and returns results in an R list.
 * Handles memory protection and type conversion following R's C interface guidelines.
 *
 * @param s_x R vector of x coordinates (numeric)
 * @param s_y R vector of y values (numeric)
 * @param s_y_true R vector of true y values, or NULL (numeric)
 * @param s_k_min Minimum number of neighbors (integer)
 * @param s_k_max Maximum number of neighbors (integer)
 * @param s_n_bb Number of bootstrap iterations (integer)
 * @param s_p Probability level for credible intervals (numeric)
 * @param s_distance_kernel Kernel function type (integer):
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param s_dist_normalization_factor Distance normalization factor (numeric)
 * @param s_epsilon Numerical stability parameter (numeric)
 * @param s_verbose Enable progress messages (logical)
 *
 * @return An R list containing:
 * - k_values: Integer vector of k values tested
 * - opt_k: Optimal k value
 * - opt_k_idx: Index of optimal k
 * - k_mean_errors: Vector of mean LOOCV errors for each k
 * - k_mean_true_errors: Vector of true errors if y_true provided, NULL otherwise
 * - predictions: Vector of model-averaged predictions
 * - k_predictions: Matrix of predictions for all k values
 * If bootstrap performed (n_bb > 0):
 * - bb_predictions: Vector of bootstrap central estimates
 * - cri_L: Vector of lower credible interval bounds
 * - cri_U: Vector of upper credible interval bounds
 *
 * @throws Rf_error if:
 * - Input vectors have inconsistent lengths
 * - Memory allocation fails
 * - Invalid parameter values provided
 *
 * @note
 * - All input vectors must have the same length
 * - x values must be sorted in ascending order
 * - Bootstrap results are NULL if n_bb = 0
 * - Uses R's protection stack for memory management
 *
 * @see mabilo() for core implementation details
 */
SEXP S_mabilo(SEXP s_x,
              SEXP s_y,
              SEXP s_y_true,
              SEXP s_k_min,
              SEXP s_k_max,
              SEXP s_n_bb,
              SEXP s_p,
              SEXP s_distance_kernel,
              SEXP s_dist_normalization_factor,
              SEXP s_epsilon,
              SEXP s_verbose) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    bool y_true_exists = LENGTH(s_y_true) == n_points;
    if (y_true_exists) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    int k_min = INTEGER(s_k_min)[0];
    int k_max = INTEGER(s_k_max)[0];

    int n_bb = INTEGER(s_n_bb)[0];
    double p = REAL(s_p)[0];

    int distance_kernel = INTEGER(s_distance_kernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    mabilo_t wmabilo_results = mabilo(x,
                                      y,
                                      y_true,
                                      k_min,
                                      k_max,
                                      n_bb,
                                      p,
                                      distance_kernel,
                                      dist_normalization_factor,
                                      epsilon,
                                      verbose);

    // Creating return list
    const int N_COMPONENTS = 10;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    std::vector<int> k_values(wmabilo_results.k_mean_errors.size());
    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
        k_values[k_index] = k;
    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(k_values)); n_protected++;

    SEXP s_opt_k = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k)[0] = wmabilo_results.opt_k;
    SET_VECTOR_ELT(result, 1, s_opt_k);

    SEXP s_opt_k_idx = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k_idx)[0] = wmabilo_results.opt_k_idx + 1;
    SET_VECTOR_ELT(result, 2, s_opt_k_idx);

    SET_VECTOR_ELT(result, 3, convert_vector_double_to_R(wmabilo_results.k_mean_errors)); n_protected++;

    // true errors
    if (y_true_exists) {
        SET_VECTOR_ELT(result, 4, convert_vector_double_to_R(wmabilo_results.k_mean_true_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 4, R_NilValue);
    }

    SET_VECTOR_ELT(result, 5, convert_vector_double_to_R(wmabilo_results.predictions)); n_protected++;
    SET_VECTOR_ELT(result, 6, convert_vector_vector_double_to_R(wmabilo_results.k_predictions)); n_protected++;

    if (n_bb > 0) {
        SET_VECTOR_ELT(result, 7, convert_vector_double_to_R(wmabilo_results.bb_predictions)); n_protected++;
        SET_VECTOR_ELT(result, 8, convert_vector_double_to_R(wmabilo_results.cri_L)); n_protected++;
        SET_VECTOR_ELT(result, 9, convert_vector_double_to_R(wmabilo_results.cri_U)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }

    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("k_values"));
    SET_STRING_ELT(names, 1, mkChar("opt_k"));
    SET_STRING_ELT(names, 2, mkChar("opt_k_idx"));
    SET_STRING_ELT(names, 3, mkChar("k_mean_errors"));
    SET_STRING_ELT(names, 4, mkChar("k_mean_true_errors"));
    SET_STRING_ELT(names, 5, mkChar("predictions"));
    SET_STRING_ELT(names, 6, mkChar("k_predictions"));
    SET_STRING_ELT(names, 7, mkChar("bb_predictions"));
    SET_STRING_ELT(names, 8, mkChar("cri_L"));
    SET_STRING_ELT(names, 9, mkChar("cri_U"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}

/**
 * @brief Smoothed error version of Model-Averaged LOWESS (MABILO) for robust local regression
 *
 * @details Similar to uwmabilo, but applies additional smoothing to the LOOCV errors
 * to reduce noise in k selection. The error curve is smoothed using uwmabilo with
 * a fixed window size of 0.25 * n_points.
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Observed y values corresponding to x (response variable)
 * @param y_true Optional true y values for error calculation
 * @param k_min Minimum number of neighbors on each side
 * @param k_max Maximum number of neighbors on each side
 * @param error_window_factor Factor to determine window size for error curve smoothing (default: 0.25).
 *                           The window size will be error_window_factor * n_points neighbors on each side.
 *                           Larger values create smoother error curves but may miss local structure.
 * @param distance_kernel Kernel function (0: Tricube, 1: Epanechnikov, 2: Exponential)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilo_t structure with additional smoothed_k_mean_errors field
 *
 * @throws std::invalid_argument for invalid input parameters
 * @throws Rf_error for numerical instability
 *
 * @note The smoothing helps prevent selecting suboptimal k values due to noise
 * in the error measurements
 */
mabilo_t mabilo_with_smoothed_errors(const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         const std::vector<double>& y_true,
                                         int k_min,
                                         int k_max,
                                         double error_window_factor,
                                         int distance_kernel,
                                         double dist_normalization_factor,
                                         double epsilon,
                                         bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILO");

    if (verbose) {
        Rprintf("Starting MABILO computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    mabilo_t results;

    auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 1: Computing models for different k values\n");
    }
    progress_tracker_t k_progress(k_max - k_min + 1, "Model computation");

    // std::vector<double> dists; // a vector of distances from the ref_pt within the window
    auto window_weights = [&x, &dist_normalization_factor](int start, int end, int ref_pt) {

        int window_size = end - start + 1;
        std::vector<double> dists(window_size);
        std::vector<double> weights(window_size);

        // Calculate distances to reference point
        double max_dist = 0.0;
        for (int i = 0; i < window_size; ++i) {
            dists[i] = std::abs(x[i + start] - x[ref_pt]);
            max_dist = std::max(max_dist, dists[i]);
        }

        if (max_dist) {
            max_dist *= dist_normalization_factor;

            // Normalize distances and compute kernel weights
            for (int i = 0; i < window_size; ++i) {
                dists[i] /= max_dist;
            }
        }

        kernel_fn(dists.data(), window_size, weights.data());

        // Normalize and rescale kernel weights by w
        double total_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (int i = 0; i < window_size; ++i)
            weights[i] = (weights[i] / total_weights);

        return weights;
    };

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    bool y_true_exists = !y_true.empty();

    int x_min_index = 0;
    int x_max_index = 0;
    int n_points_minus_one = n_points - 1;
    std::vector<double> w_window;
    int n_k_values = k_max - k_min + 1;

    // Storage for predictions across all k values
    std::vector<std::vector<double>> k_sm_predictions(n_k_values, std::vector<double>(n_points));
    std::vector<std::vector<double>> k_predictions(n_k_values, std::vector<double>(n_points));

    // Vectors for errors during single k iteration
    std::vector<double> k_errors(n_points);
    std::vector<double> k_true_errors(n_points);

    // Vectors in results struct to store mean errors for each k
    results.k_mean_errors.resize(n_k_values);
    results.k_mean_true_errors.resize(n_k_values);

    // Pre-allocate vectors outside the k loop with maximum possible sizes
    std::vector<std::pair<double, const ulm_plus_t*>> filtered_models;
    filtered_models.reserve(2 * k_max + 1);  // Maximum window size for any k

    std::vector<double> local_errors;
    local_errors.reserve(2 * k_max + 1);  // Maximum number of models for a point

    std::vector<double> all_errors;

    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++) {
        auto k_ptm = std::chrono::steady_clock::now();
        if (verbose) {
            Rprintf("\nProcessing k=%d (%d/%d) ... ",
                   k, k_index + 1, k_max - k_min + 1);
        }

        std::vector<std::vector<ulm_plus_t>> pt_models(n_points);
        for (int i = 0; i < n_points; i++) {
            pt_models[i].reserve(2 * k + 1);
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        // Phase 1: Creating single model predictions
        //
        // For each x value (reference point) create a window of size 2k + 1
        // around that value and fit a weighted linear model to x and y
        // restricted to the window with weights defined by a kernel on the
        // distances of restricted x values from the reference point. The value
        // of the weight vector at the boundary end of the window that is the
        // furthest from the reference point is close to 0 (this is where we use
        // dist_normalization_factor). The weights are symmetric around the
        // reference point.
        //
        // For each point (i - the index of that point) in the support of the given
        // model we insert that model in the vector pt_models[i] of all models that
        // have 'i' in their support.
        //
        // Here we may restrict ourselves to include only the model in pt_models[i] if
        // 'i' is not too far from the referece point of the model.

        if (verbose) {
            Rprintf("  Phase 1: Computing single-model predictions ... ");
        }
        auto phase1_ptm = std::chrono::steady_clock::now();

        initialize_kernel(distance_kernel, 1.0);

        for (int i = 0; i < n_points; i++) {

            // find the start and the end indices of the window around a ref_pt (x value) so that ref_pt is as much as possible in the middle of the window
            if (i > k_minus_one && i < n_points_minus_k) {
                x_min_index = i - k; // the first condition implies that x_min_index >= 0
                x_max_index = i + k; // the second condition implies that x_min_index < n_points
            } else if (i < k) {
                x_min_index = 0;
                x_max_index = two_k;
            } else if (i > n_points_minus_k_minus_one) {
                x_min_index = n_points_minus_one_minus_two_k;
                x_max_index = n_points_minus_one;
            }

            // Computing window weights
            w_window = window_weights(x_min_index, x_max_index, i);

            // Fitting a weighted linear model
            ulm_t fit = ulm(x.data() + x_min_index,
                            y.data() + x_min_index,
                            w_window,
                            y_binary,
                            epsilon);
            ulm_plus_t wlm_res;
            wlm_res.predictions = std::move(fit.predictions);
            wlm_res.errors      = std::move(fit.errors);
            wlm_res.x_min_index = x_min_index;
            wlm_res.x_max_index = x_max_index;

            // Store single-model prediction and its LOOCV squared errors
            k_errors[i] = wlm_res.errors[i - x_min_index];
            double y_prediction_at_ref_pt = wlm_res.predictions[i - x_min_index];
            k_sm_predictions[k_index][i] = y_prediction_at_ref_pt;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - y_prediction_at_ref_pt);
            }

            // For x indices around i insert that model into pt_models[i]
            for (int j = x_min_index; j <= x_max_index; j++) {
                wlm_res.w = w_window;  // Store the weights for later model averaging
                pt_models[j].push_back(wlm_res);
            }
        }

        if (verbose) {
            elapsed_time(phase1_ptm, "Done");
            mem_tracker.report();
        }

        // Phase 2: Model averaging
        // k_errors and k_true_errors are reused for model-averaged predictions
        if (verbose) {
            Rprintf("  Phase 2: Computing model-averaged predictions ... ");
        }
        auto phase2_ptm = std::chrono::steady_clock::now();


        for (int i = 0; i < n_points; i++) {

            double weighted_sum = 0.0;
            double weight_sum = 0.0;
            //int model_counter = 0;
            //local_errors.resize(pt_models[i].size());
            double wmean_error = 0.0;

            for (const auto& model : pt_models[i]) {
                int local_index = i - model.x_min_index;
                double weight = model.w[local_index];
                weighted_sum += weight * model.predictions[local_index];
                weight_sum += weight;
                wmean_error +=  weight * model.errors[local_index];
                //local_errors[model_counter++] = model.errors[local_index];
            }

            k_errors[i] = wmean_error / weight_sum;

            // Store model-averaged prediction and its errors
            k_predictions[k_index][i] = weighted_sum / weight_sum;

            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - k_predictions[k_index][i]);
            }
        }

        // Compute mean errors for model-averaged predictions at current k
        results.k_mean_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
        }

        if (verbose) {
            elapsed_time(phase2_ptm, "Done");
            mem_tracker.report();
        }

        if (verbose) {
            char message[100];  // Buffer large enough for the message
            snprintf(message, sizeof(message), "\nTotal time for k=%d: ", k);
            elapsed_time(k_ptm, message);
            k_progress.update(k_index + 1);
        }
    }

    if (verbose) {
        elapsed_time(models_ptm, "\nTotal model computation time: ");
    }

    // -------------------------------------------------------------------------------------------
    //
    // Phase 3: Find optimal k for model-averaged predictions
    // using approximated LOOCV squared errros for these predictions
    //
    // -------------------------------------------------------------------------------------------
    auto opt_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 3: Finding optimal  model averaged predictions over all k's ... ");
    }

    if (k_max > k_min) {

        // smoothing results.k_mean_errors using uwmabilo() with k = 0.25 * n_points
        int n_k_values = k_max - k_min + 1;
        int window_size = static_cast<int>(error_window_factor * n_k_values);
        // 2k + 1 = window_size => k = (window_size - 1) / 2
        int error_k = (window_size - 1) / 2;
        if (error_k < 1) error_k = 1;

        std::vector<double> k_values(n_k_values);
        for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
            k_values[k_index] = k;

        std::vector<double> true_errors;
        auto errors_fit = uwmabilo(k_values,
                                   results.k_mean_errors,
                                   true_errors,
                                   error_k,
                                   error_k,
                                   distance_kernel,
                                   dist_normalization_factor,
                                   epsilon,
                                   verbose);

        results.smoothed_k_mean_errors = std::move(errors_fit.predictions);

        auto min_it = std::min_element(results.smoothed_k_mean_errors.begin(), results.smoothed_k_mean_errors.end());
        results.opt_k_idx = std::distance(results.smoothed_k_mean_errors.begin(), min_it);

    } else {
        results.opt_k_idx = 0;
    }
    results.opt_k = k_min + results.opt_k_idx;
    results.predictions = k_predictions[results.opt_k_idx];
    results.k_predictions = std::move(k_predictions);
    if (verbose) {
        elapsed_time(opt_k_ptm, "Done");
        mem_tracker.report();
    }


    if (verbose) {
        elapsed_time(total_ptm, "\nTotal MABILO computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}

/**
 * @brief R interface to the smoothed error MABILO implementation
 *
 * @details Converts R objects to C++ types, calls mabilo_with_smoothed_errors(), and returns results
 * as an R list. Handles memory management and R object protection.
 *
 * @param s_x R vector of x values
 * @param s_y R vector of y values
 * @param s_y_true R vector of true y values (can be NULL)
 * @param s_k_min R integer for minimum k
 * @param s_k_max R integer for maximum k
 * @param s_error_window_factor Factor to determine window size for error curve smoothing (default: 0.25).
 *                           The window size will be error_window_factor * n_points neighbors on each side.
 *                           Larger values create smoother error curves but may miss local structure.
 * @param s_distance_kernel R integer for kernel type
 * @param s_dist_normalization_factor R numeric for distance normalization
 * @param s_epsilon R numeric for numerical stability
 * @param s_verbose R logical for progress messages
 *
 * @return R list containing:
 *         - k_values: Vector of k values used
 *         - opt_k: Optimal k value
 *         - opt_k_idx: Index of optimal k
 *         - k_mean_errors: Raw error values
 *         - smoothed_k_mean_errors: Smoothed error values
 *         - k_mean_true_errors: True errors (if y_true provided)
 *         - predictions: Final predictions
 *         - k_predictions: Predictions for all k values
 */
SEXP S_mabilo_with_smoothed_errors(SEXP s_x,
                                     SEXP s_y,
                                     SEXP s_y_true,
                                     SEXP s_k_min,
                                     SEXP s_k_max,
                                     SEXP s_error_window_factor,
                                     SEXP s_distance_kernel,
                                     SEXP s_dist_normalization_factor,
                                     SEXP s_epsilon,
                                     SEXP s_verbose) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    bool y_true_exists = LENGTH(s_y_true) == n_points;
    if (y_true_exists) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    int k_min = INTEGER(s_k_min)[0];
    int k_max = INTEGER(s_k_max)[0];

    double error_window_factor = REAL(s_error_window_factor)[0];

    int distance_kernel = INTEGER(s_distance_kernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    mabilo_t mabilo_with_smoothed_errors_results = mabilo_with_smoothed_errors(x,
                                                                               y,
                                                                               y_true,
                                                                               k_min,
                                                                               k_max,
                                                                               error_window_factor,
                                                                               distance_kernel,
                                                                               dist_normalization_factor,
                                                                               epsilon,
                                                                               verbose);

    // Creating return list
    const int N_COMPONENTS = 8;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    std::vector<int> k_values(mabilo_with_smoothed_errors_results.k_mean_errors.size());
    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
        k_values[k_index] = k;
    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(k_values)); n_protected++;

    SEXP s_opt_k = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k)[0] = mabilo_with_smoothed_errors_results.opt_k + 1;
    SET_VECTOR_ELT(result, 1, s_opt_k);

    SEXP s_opt_k_idx = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_k_idx)[0] = mabilo_with_smoothed_errors_results.opt_k_idx + 1;
    SET_VECTOR_ELT(result, 2, s_opt_k_idx);

    SET_VECTOR_ELT(result, 3, convert_vector_double_to_R(mabilo_with_smoothed_errors_results.k_mean_errors)); n_protected++;

    SET_VECTOR_ELT(result, 4, convert_vector_double_to_R(mabilo_with_smoothed_errors_results.smoothed_k_mean_errors)); n_protected++;

    // true errors
    if (y_true_exists) {
        SET_VECTOR_ELT(result, 5, convert_vector_double_to_R(mabilo_with_smoothed_errors_results.k_mean_true_errors)); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 5, R_NilValue);
    }

    SET_VECTOR_ELT(result, 6, convert_vector_double_to_R(mabilo_with_smoothed_errors_results.predictions)); n_protected++;
    SET_VECTOR_ELT(result, 7, convert_vector_vector_double_to_R(mabilo_with_smoothed_errors_results.k_predictions)); n_protected++;

    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("k_values"));
    SET_STRING_ELT(names, 1, mkChar("opt_k"));
    SET_STRING_ELT(names, 2, mkChar("opt_k_idx"));
    SET_STRING_ELT(names, 3, mkChar("k_mean_errors"));
    SET_STRING_ELT(names, 4, mkChar("smoothed_k_mean_errors"));
    SET_STRING_ELT(names, 5, mkChar("k_mean_true_errors"));
    SET_STRING_ELT(names, 6, mkChar("predictions"));
    SET_STRING_ELT(names, 7, mkChar("k_predictions"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}



/**
 * @brief Implements the Model-Averaged Bi-kNN LOcal linear model (MABILO) algorithm
 *
 * @details MABILO extends traditional LOWESS by incorporating model averaging with
 * bi-k nearest neighbor structure. The algorithm consists of three phases:
 *
 * Phase 1 - Single Model Computation:
 * - For each k in [k_min, k_max]:
 *   - Fit local linear models using k-hop neighborhoods
 *   - Compute predictions and their LOOCV errors
 *   - Store model information for each point in support
 *
 * Phase 2 - Model Averaging:
 * - For each point, compute weighted average predictions using:
 *   - All models containing the point in their support
 *   - Original kernel weights from model fitting
 *   - Weighted average of LOOCV errors
 *
 * Phase 3 - Optimal k Selection:
 * - Find k with minimum mean LOOCV error
 * - Return corresponding predictions and errors
 *
 * The algorithm uses k-hop neighbors instead of k-nearest neighbors, providing
 * more symmetric neighborhoods in 1D data.
 *
 * @param x Vector of predictor values (sorted in ascending order)
 * @param y Vector of response values corresponding to x
 * @param y_true Optional vector of true values for error calculation
 * @param w Vector of sample weights
 * @param k_min Minimum number of neighbors to consider
 * @param k_max Maximum number of neighbors to consider
 * @param distance_kernel Integer specifying kernel type for distance weighting
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small constant for numerical stability in model fitting
 * @param verbose Flag for detailed progress output
 *
 * @return mabilo_t structure containing:
 *   - opt_k: Optimal k value
 *   - opt_k_idx: Index of optimal k
 *   - predictions: Model-averaged predictions using optimal k
 *   - k_mean_errors: Mean LOOCV errors for each k
 *   - k_mean_true_errors: Mean absolute errors vs y_true (if provided)
 *   - k_predictions: Model-averaged predictions for all k values
 *
 * @pre
 * - x must be sorted in ascending order
 * - All input vectors must have the same size
 * - k_min <= k_max
 * - k_max <= (n_points - 1)/2
 *
 * @note
 * - The algorithm handles binary response variables (y ∈ {0,1})
 * - For boundary points, windows are adjusted to maintain constant width
 * - Memory usage scales with k_max and number of points
 *
 * @see ulm_t
 * @see kernel_fn
 */
mabilo_t wmabilo(const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& y_true,
                 const std::vector<double>& w,
                 int k_min,
                 int k_max,
                 int distance_kernel,
                 double dist_normalization_factor,
                 double epsilon,
                 bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILO");

    if (verbose) {
        Rprintf("Starting MABILO computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    initialize_kernel(distance_kernel, 1.0);

    auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 1: Computing models for different k values\n");
    }
    progress_tracker_t k_progress(k_max - k_min + 1, "Model computation");

    //------------------------------------------------------------------------------
    // Algorithm Overview
    //------------------------------------------------------------------------------
    // MABILO (Model-Averaged Bi-kNN LOcal linear model) implements local linear
    // regression with model averaging over k-hop neighborhoods. For each point,
    // it fits models using symmetric windows and averages them using kernel weights.

    //------------------------------------------------------------------------------
    // Window Weight Computation (Lambda Function)
    //------------------------------------------------------------------------------
    auto window_weights = [&x, &dist_normalization_factor, &w](int start, int end, int ref_pt) {

        // Computes kernel weights for a window of points:
        // 1. Calculates distances from reference point
        // 2. Normalizes distances using max distance
        // 3. Applies kernel function
        // 4. Normalizes weights and combines with sample weights

        int window_size = end - start + 1;
        std::vector<double> dists(window_size);
        std::vector<double> weights(window_size);

        // Calculate distances to reference point
        double max_dist = 0.0;
        for (int i = 0; i < window_size; ++i) {
            dists[i] = std::abs(x[i + start] - x[ref_pt]);
            max_dist = std::max(max_dist, dists[i]);
        }

        if (max_dist) {
            max_dist *= dist_normalization_factor;

            // Normalize distances and compute kernel weights
            for (int i = 0; i < window_size; ++i) {
                dists[i] /= max_dist;
            }
        }

        kernel_fn(dists.data(), window_size, weights.data());

        // Normalize and rescale kernel weights by w
        double total_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (int i = 0; i < window_size; ++i)
            weights[i] = (weights[i] / total_weights) * w[i + start];

        return weights;
    };

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    bool y_true_exists = !y_true.empty();

    int x_min_index = 0;
    int x_max_index = 0;
    int n_points_minus_one = n_points - 1;
    std::vector<double> w_window;
    int n_k_values = k_max - k_min + 1;

    // Storage for predictions across all k values
    std::vector<std::vector<double>> k_predictions(n_k_values, std::vector<double>(n_points));

    // Vectors for errors during single k iteration
    std::vector<double> k_errors(n_points);
    std::vector<double> k_true_errors(n_points);

    // Vectors in results struct to store mean errors for each k
    mabilo_t results;
    results.k_mean_errors.resize(n_k_values);
    results.k_mean_true_errors.resize(n_k_values);

    // Pre-allocate vectors outside the k loop with maximum possible sizes
    std::vector<std::pair<double, const ulm_plus_t*>> filtered_models;
    filtered_models.reserve(2 * k_max + 1);  // Maximum window size for any k

    std::vector<double> local_errors;
    local_errors.reserve(2 * k_max + 1);  // Maximum number of models for a point

    std::vector<double> all_errors;

    struct pred_w_err_t {
        double prediction;
        double weight;
        double error;
    } pred_w_err;

    std::vector<std::vector<pred_w_err_t>> pt_pred_w_err(n_points);
    for (int i = 0; i < n_points; i++) {
        pt_pred_w_err[i].reserve(2 * k_max + 1);
    }

    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++) {
        auto k_ptm = std::chrono::steady_clock::now();
        if (verbose) {
            Rprintf("\nProcessing k=%d (%d/%d) ... ",
                    k, k_index + 1, k_max - k_min + 1);
        }

        for (int i = 0; i < n_points; i++) {
            pt_pred_w_err[i].clear();
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        //------------------------------------------------------------------------------
        // Phase 1: Single Model Computation
        //------------------------------------------------------------------------------
        // For each k from k_min to k_max:
        //   For each point x[i]:
        //     1. Define window:
        //        - Interior points: symmetric k-hop neighborhood
        //        - Boundary points: adjusted window size maintaining total width
        //     2. Compute kernel weights based on distances
        //     3. Fit local linear model using weighted least squares
        //     4. Store predictions, weights, and errors for each point in window
        if (verbose) {
            Rprintf("  Phase 1: Computing single-model predictions ... ");
        }
        auto phase1_ptm = std::chrono::steady_clock::now();

        for (int i = 0; i < n_points; i++) {

            // find the start and the end indices of the window around a ref_pt (x value) so that ref_pt is as much as possible in the middle of the window
            if (i > k_minus_one && i < n_points_minus_k) {
                x_min_index = i - k; // the first condition implies that x_min_index >= 0
                x_max_index = i + k; // the second condition implies that x_min_index < n_points
            } else if (i < k) {
                x_min_index = 0;
                x_max_index = two_k;
            } else if (i > n_points_minus_k_minus_one) {
                x_min_index = n_points_minus_one_minus_two_k;
                x_max_index = n_points_minus_one;
            }

            // Computing window weights
            w_window = window_weights(x_min_index, x_max_index, i);

            // Fitting a weighted linear model
            ulm_t wlm_fit = ulm(x.data() + x_min_index,
                                y.data() + x_min_index,
                                w_window,
                                y_binary,
                                epsilon);

            // For each point of the window record predicted value, the weight, and the models LOOCV at that point
            int x_max_index_plus_one = x_max_index + 1;
            for (int s = 0, j = x_min_index; j < x_max_index_plus_one; s++, j++) {
                pred_w_err.prediction = wlm_fit.predictions[s];
                pred_w_err.weight     = w_window[s];
                pred_w_err.error      = wlm_fit.errors[s];
                pt_pred_w_err[j].push_back(pred_w_err);
            }
        }

        if (verbose) {
            elapsed_time(phase1_ptm, "Done");
            mem_tracker.report();
        }

        //------------------------------------------------------------------------------
        // Phase 2: Model Averaging
        //------------------------------------------------------------------------------
        // For each point x[i]:
        //   1. Collect all models containing the point
        //   2. Compute weighted average of predictions using kernel weights
        //   3. Compute weighted average of LOOCV errors
        //   4. If true values provided, compute absolute prediction errors
        if (verbose) {
            Rprintf("  Phase 2: Computing model-averaged predictions ... ");
        }
        auto phase2_ptm = std::chrono::steady_clock::now();

        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        double wmean_error = 0.0;

        for (int i = 0; i < n_points; i++) {
            weighted_sum = 0.0;
            weight_sum = 0.0;
            wmean_error = 0.0;
            for (const auto& v : pt_pred_w_err[i]) {
                weighted_sum += v.weight * v.prediction;
                weight_sum   += v.weight;
                wmean_error  += v.weight * v.error;
            }

            k_errors[i] = wmean_error / weight_sum;
            k_predictions[k_index][i] = weighted_sum / weight_sum;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - k_predictions[k_index][i]);
            }
        }

        // Compute mean errors for model-averaged predictions at current k
        results.k_mean_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
        }

        if (verbose) {
            elapsed_time(phase2_ptm, "Done");
            mem_tracker.report();
        }

        if (verbose) {
            char message[100];  // Buffer large enough for the message
            snprintf(message, sizeof(message), "\nTotal time for k=%d: ", k);
            elapsed_time(k_ptm, message);
            k_progress.update(k_index + 1);
        }
    }

    if (verbose) {
        elapsed_time(models_ptm, "\nTotal model computation time: ");
    }

    //------------------------------------------------------------------------------
    // Phase 3: Optimal k Selection
    //------------------------------------------------------------------------------
    // 1. Compare mean LOOCV errors across different k values
    // 2. Select k with minimum mean error
    // 3. Store corresponding predictions and errors
    auto opt_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 3: Finding optimal  model averaged predictions over all k's ... ");
    }

    if (k_max > k_min) {
        auto min_it = std::min_element(results.k_mean_errors.begin(), results.k_mean_errors.end());
        results.opt_k_idx = std::distance(results.k_mean_errors.begin(), min_it);
    } else {
        results.opt_k_idx = 0;
    }
    results.opt_k = k_min + results.opt_k_idx;
    results.predictions = k_predictions[results.opt_k_idx];
    results.k_predictions = std::move(k_predictions);
    if (verbose) {
        elapsed_time(opt_k_ptm, "Done");
        mem_tracker.report();
    }

    if (verbose) {
        elapsed_time(total_ptm, "\nTotal MABILO computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}


/**
 * @brief Implements the Model-Averaged Bi-kNN LOcal linear model (MABILO) algorithm
 *
 * @details MABILO extends traditional LOWESS by incorporating model averaging with
 * bi-k nearest neighbor structure. The algorithm consists of three phases:
 *
 * Phase 1 - Single Model Computation:
 * - For each k in [k_min, k_max]:
 *   - Fit local linear models using k-hop neighborhoods
 *   - Compute predictions and their LOOCV errors
 *   - Store model information for each point in support
 *
 * Phase 2 - Model Averaging:
 * - For each point, compute weighted average predictions using:
 *   - All models containing the point in their support
 *   - Original kernel weights from model fitting
 *   - Weighted average of LOOCV errors
 *
 * Phase 3 - Optimal k Selection:
 * - Find k with minimum mean LOOCV error
 * - Return corresponding predictions and errors
 *
 * The algorithm uses k-hop neighbors instead of k-nearest neighbors, providing
 * more symmetric neighborhoods in 1D data.
 *
 * @param x Vector of predictor values (sorted in ascending order)
 * @param y Vector of response values corresponding to x
 * @param y_true Optional vector of true values for error calculation
 * @param k_min Minimum number of neighbors to consider
 * @param k_max Maximum number of neighbors to consider
 * @param distance_kernel Integer specifying kernel type for distance weighting
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small constant for numerical stability in model fitting
 * @param verbose Flag for detailed progress output
 *
 * @return mabilo_t structure containing:
 *   - opt_k: Optimal k value
 *   - opt_k_idx: Index of optimal k
 *   - predictions: Model-averaged predictions using optimal k
 *   - k_mean_errors: Mean LOOCV errors for each k
 *   - k_mean_true_errors: Mean absolute errors vs y_true (if provided)
 *   - k_predictions: Model-averaged predictions for all k values
 *
 * @pre
 * - x must be sorted in ascending order
 * - All input vectors must have the same size
 * - k_min <= k_max
 * - k_max <= (n_points - 1)/2
 *
 * @note
 * - The algorithm handles binary response variables (y ∈ {0,1})
 * - For boundary points, windows are adjusted to maintain constant width
 * - Memory usage scales with k_max and number of points
 *
 * @see ulm_t
 * @see kernel_fn
 */
mabilo_t uwmabilo(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& y_true,
                  int k_min,
                  int k_max,
                  int distance_kernel,
                  double dist_normalization_factor,
                  double epsilon,
                  bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILO");

    if (verbose) {
        Rprintf("Starting MABILO computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    initialize_kernel(distance_kernel, 1.0);

    auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 1: Computing models for different k values\n");
    }
    progress_tracker_t k_progress(k_max - k_min + 1, "Model computation");

    //------------------------------------------------------------------------------
    // Algorithm Overview
    //------------------------------------------------------------------------------
    // MABILO (Model-Averaged Bi-kNN LOcal linear model) implements local linear
    // regression with model averaging over k-hop neighborhoods. For each point,
    // it fits models using symmetric windows and averages them using kernel weights.

    //------------------------------------------------------------------------------
    // Window Weight Computation (Lambda Function)
    //------------------------------------------------------------------------------
    auto window_weights = [&x, &dist_normalization_factor](int start, int end, int ref_pt) {

        // Computes kernel weights for a window of points:
        // 1. Calculates distances from reference point
        // 2. Normalizes distances using max distance
        // 3. Applies kernel function
        // 4. Normalizes weights and combines with sample weights

        int window_size = end - start + 1;
        std::vector<double> dists(window_size);
        std::vector<double> weights(window_size);

        // Calculate distances to reference point
        double max_dist = 0.0;
        for (int i = 0; i < window_size; ++i) {
            dists[i] = std::abs(x[i + start] - x[ref_pt]);
            max_dist = std::max(max_dist, dists[i]);
        }

        if (max_dist) {
            max_dist *= dist_normalization_factor;

            // Normalize distances and compute kernel weights
            for (int i = 0; i < window_size; ++i) {
                dists[i] /= max_dist;
            }
        }

        kernel_fn(dists.data(), window_size, weights.data());

        // Normalize and rescale kernel weights by w
        double total_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (int i = 0; i < window_size; ++i)
            weights[i] = (weights[i] / total_weights);

        return weights;
    };

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    bool y_true_exists = !y_true.empty();

    int x_min_index = 0;
    int x_max_index = 0;
    int n_points_minus_one = n_points - 1;
    std::vector<double> w_window;
    int n_k_values = k_max - k_min + 1;

    // Storage for predictions across all k values
    std::vector<std::vector<double>> k_predictions(n_k_values, std::vector<double>(n_points));

    // Vectors for errors during single k iteration
    std::vector<double> k_errors(n_points);
    std::vector<double> k_true_errors(n_points);

    // Vectors in results struct to store mean errors for each k
    mabilo_t results;
    results.k_mean_errors.resize(n_k_values);
    results.k_mean_true_errors.resize(n_k_values);

    // Pre-allocate vectors outside the k loop with maximum possible sizes
    std::vector<std::pair<double, const ulm_plus_t*>> filtered_models;
    filtered_models.reserve(2 * k_max + 1);  // Maximum window size for any k

    std::vector<double> local_errors;
    local_errors.reserve(2 * k_max + 1);  // Maximum number of models for a point

    std::vector<double> all_errors;

    struct pred_w_err_t {
        double prediction;
        double weight;
        double error;
    } pred_w_err;

    std::vector<std::vector<pred_w_err_t>> pt_pred_w_err(n_points);
    for (int i = 0; i < n_points; i++) {
        pt_pred_w_err[i].reserve(10 * k_max + 1);
    }

    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++) {
        auto k_ptm = std::chrono::steady_clock::now();
        if (verbose) {
            Rprintf("\nProcessing k=%d (%d/%d) ... ",
                    k, k_index + 1, k_max - k_min + 1);
        }

        for (int i = 0; i < n_points; i++) {
            pt_pred_w_err[i].clear();
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        //------------------------------------------------------------------------------
        // Phase 1: Single Model Computation
        //------------------------------------------------------------------------------
        // For each k from k_min to k_max:
        //   For each point x[i]:
        //     1. Define window:
        //        - Interior points: symmetric k-hop neighborhood
        //        - Boundary points: adjusted window size maintaining total width
        //     2. Compute kernel weights based on distances
        //     3. Fit local linear model using weighted least squares
        //     4. Store predictions, weights, and errors for each point in window
        if (verbose) {
            Rprintf("  Phase 1: Computing single-model predictions ... ");
        }
        auto phase1_ptm = std::chrono::steady_clock::now();

        for (int i = 0; i < n_points; i++) {

            // find the start and the end indices of the window around a ref_pt (x value) so that ref_pt is as much as possible in the middle of the window
            if (i > k_minus_one && i < n_points_minus_k) {
                x_min_index = i - k; // the first condition implies that x_min_index >= 0
                x_max_index = i + k; // the second condition implies that x_min_index < n_points
            } else if (i < k) {
                x_min_index = 0;
                x_max_index = two_k;
            } else if (i > n_points_minus_k_minus_one) {
                x_min_index = n_points_minus_one_minus_two_k;
                x_max_index = n_points_minus_one;
            }

            // Computing window weights
            w_window = window_weights(x_min_index, x_max_index, i);

            // Fitting a weighted linear model
            ulm_t wlm_fit = ulm(x.data() + x_min_index,
                                y.data() + x_min_index,
                                w_window,
                                y_binary,
                                epsilon);

            // For each point of the window record predicted value, the weight, and the models LOOCV at that point
            int x_max_index_plus_one = x_max_index + 1;
            for (int s = 0, j = x_min_index; j < x_max_index_plus_one; s++, j++) {
                pred_w_err.prediction = wlm_fit.predictions[s];
                pred_w_err.weight     = w_window[s];
                pred_w_err.error      = wlm_fit.errors[s];
                pt_pred_w_err[j].push_back(pred_w_err);
            }
        }

        if (verbose) {
            elapsed_time(phase1_ptm, "Done");
            mem_tracker.report();
        }

        //------------------------------------------------------------------------------
        // Phase 2: Model Averaging
        //------------------------------------------------------------------------------
        // For each point x[i]:
        //   1. Collect all models containing the point
        //   2. Compute weighted average of predictions using kernel weights
        //   3. Compute weighted average of LOOCV errors
        //   4. If true values provided, compute absolute prediction errors
        if (verbose) {
            Rprintf("  Phase 2: Computing model-averaged predictions ... ");
        }
        auto phase2_ptm = std::chrono::steady_clock::now();

        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        double wmean_error = 0.0;

        for (int i = 0; i < n_points; i++) {
            weighted_sum = 0.0;
            weight_sum = 0.0;
            wmean_error = 0.0;
            for (const auto& v : pt_pred_w_err[i]) {
                weighted_sum += v.weight * v.prediction;
                weight_sum   += v.weight;
                wmean_error  += v.weight * v.error;
            }

            k_errors[i] = wmean_error / weight_sum;
            k_predictions[k_index][i] = weighted_sum / weight_sum;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - k_predictions[k_index][i]);
            }
        }

        // Compute mean errors for model-averaged predictions at current k
        results.k_mean_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
        }

        if (verbose) {
            elapsed_time(phase2_ptm, "Done");
            mem_tracker.report();
        }

        if (verbose) {
            char message[100];  // Buffer large enough for the message
            snprintf(message, sizeof(message), "\nTotal time for k=%d: ", k);
            elapsed_time(k_ptm, message);
            k_progress.update(k_index + 1);
        }
    }

    if (verbose) {
        elapsed_time(models_ptm, "\nTotal model computation time: ");
    }

    //------------------------------------------------------------------------------
    // Phase 3: Optimal k Selection
    //------------------------------------------------------------------------------
    // 1. Compare mean LOOCV errors across different k values
    // 2. Select k with minimum mean error
    // 3. Store corresponding predictions and errors
    auto opt_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 3: Finding optimal  model averaged predictions over all k's ... ");
    }

    if (k_max > k_min) {
        auto min_it = std::min_element(results.k_mean_errors.begin(), results.k_mean_errors.end());
        results.opt_k_idx = std::distance(results.k_mean_errors.begin(), min_it);
    } else {
        results.opt_k_idx = 0;
    }
    results.opt_k = k_min + results.opt_k_idx;
    results.predictions = k_predictions[results.opt_k_idx];
    results.k_predictions = std::move(k_predictions);
    if (verbose) {
        elapsed_time(opt_k_ptm, "Done");
        mem_tracker.report();
    }

    if (verbose) {
        elapsed_time(total_ptm, "\nTotal MABILO computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}
