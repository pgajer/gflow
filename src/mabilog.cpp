#include "exec_policy.hpp"
#include "sampling.h" // for C_runif_simplex()
#include "ulogit.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "kernels.h"  // for initialize_kernel()
#include "memory_utils.hpp"
#include "progress_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"                 // for elapsed_time
#include "predictive_errors.hpp"

#include <execution>
#include <atomic>
#include <mutex>
#include <numeric>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()

#include <R.h>
#include <Rinternals.h>

extern "C" {

    SEXP S_wmabilog(SEXP s_x,
                    SEXP s_y,
                    SEXP s_y_true,
                    SEXP s_w,
                    SEXP s_max_iterations,
                    SEXP s_ridge_lambda,
                    SEXP s_max_beta,
                    SEXP s_tolerance,
                    SEXP s_k_min,
                    SEXP s_k_max,
                    SEXP s_distance_kernel,
                    SEXP s_dist_normalization_factor,
                    SEXP s_verbose);

    SEXP S_mabilog(SEXP s_x,
                   SEXP s_y,
                   SEXP s_y_true,
                   SEXP s_max_iterations,
                   SEXP s_ridge_lambda,
                   SEXP s_max_beta,
                   SEXP s_tolerance,
                   SEXP s_k_min,
                   SEXP s_k_max,
                   SEXP s_n_bb,
                   SEXP s_p,
                   SEXP s_distance_kernel,
                   SEXP s_dist_normalization_factor,
                   SEXP s_verbose);

    SEXP S_mabilog_with_smoothed_errors(SEXP s_x,
                                        SEXP s_y,
                                        SEXP s_y_true,
                                        SEXP s_max_iterations,
                                        SEXP s_ridge_lambda,
                                        SEXP s_max_beta,
                                        SEXP s_tolerance,
                                        SEXP s_k_min,
                                        SEXP s_k_max,
                                        SEXP s_error_window_factor,
                                        SEXP s_distance_kernel,
                                        SEXP s_dist_normalization_factor,
                                        SEXP s_verbose);
}

struct mabilog_t {
    // k values
    int opt_k;     // optimal model averaging k value - the one with the smallest mean LOOCV Rf_error
    int opt_k_idx; // optimal model averaging k value index

    // Errors
    std::vector<double> k_mean_errors;   // mean LOOCV squared errors for each k for model averaged predictions
    std::vector<double> smoothed_k_mean_errors;
    std::vector<double> k_mean_true_errors; // mean absolute Rf_error between predictions and y_true

    // The best (over all k) model evaluation
    std::vector<double> predictions; // optimal k model averaged predictions

    std::vector<std::vector<double>> k_predictions; // for each k model averaged predictions

    // Bayesian bootstrap creadible intervals
    std::vector<double> bb_predictions; // central location of the Bayesian bootstrap estimates
    std::vector<double> cri_L; // credible intervals lower limit
    std::vector<double> cri_U; // credible intervals upper limit
};

/**
 * @brief Implements the Model-Averaged Bi-kNN local logistic regression  (MABILOG) algorithm
 *
 * @details MABILOG extends the traditional LOWESS algorithm by incorporating model averaging with
 * bi-k nearest neighbor structure. The algorithm consists of three phases:
 *
 * 1. For each k in [k_min, k_max]:
 *    - Fit local linear models using k-hop neighborhoods
 *    - Compute predictions and their LOOCV errors
 *    - Perform kernel-weighted model averaging
 *    - Calculate mean LOOCV and true errors for predictions
 *
 * 2. For each point, compute model-averaged predictions using kernel weights
 *    - Models contribute to predictions based on their kernel weights
 *    - Kernel weights reflect the distance between points
 *    - Final prediction is a weighted average of all contributing models
 *
 * 3. Find optimal k for model-averaged predictions based on weighted LOOCV errors
 *
 * The algorithm uses k-hop neighbors instead of k-nearest neighbors, providing more
 * symmetric neighborhoods in 1D data. For each point x_i, its k-hop neighborhood
 * includes k points to the left and k points to the right when available.
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Observed y values corresponding to x (response variable)
 * @param y_true Optional true y values for Rf_error calculation. Used for algorithm evaluation
 * @param w Observation weights, typically from Bayesian bootstrap
 * @param k_min Minimum number of neighbors on each side (minimum window half-width)
 * @param k_max Maximum number of neighbors on each side (maximum window half-width)
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 *        Makes weights at window boundaries close to but not equal to 0
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilo_t structure containing:
 *         - opt_k: Optimal k value for model-averaged predictions
 *         - predictions: Model-averaged predictions using opt_k
 *         - k_mean_errors: Mean LOOCV errors for model-averaged predictions for each k
 *         - k_mean_true_errors: Mean true errors for model-averaged predictions if y_true provided
 *         - k_predictions: Model-averaged predictions for all k values
 *
 * @throws std::invalid_argument for invalid input parameters
 * @throws Rf_error for numerical instability
 *
 * @note
 * - Input x values must be sorted in ascending order
 * - Window size for each k is 2k + 1 (k points on each side plus the center point)
 * - For points near boundaries, the window is adjusted to include available points
 * - Model averaging uses kernel weights based on distances between points
 *
 * @see ulogit_t for the structure of individual linear models
 * @see initialize_kernel() for kernel function details
 */
mabilog_t wmabilog(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& y_true,
                   const std::vector<double>& w,
                   int max_iterations,
                   double ridge_lambda,
                   double max_beta,
                   double tolerance,
                   int k_min,
                   int k_max,
                   int distance_kernel,
                   double dist_normalization_factor,
                   bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILOG");

    if (verbose) {
        Rprintf("Starting MABILOG computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    mabilog_t results;

    auto models_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 1: Computing models for different k values\n");
    }
    progress_tracker_t k_progress(k_max - k_min + 1, "Model computation");

    // std::vector<double> dists; // a vector of distances from the ref_pt within the window
    auto window_weights = [&x, &dist_normalization_factor, &w](int start, int end, int ref_pt) {

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

    std::vector<std::pair<double, const ulogit_t*>> filtered_models;
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

        std::vector<std::vector<ulogit_t>> pt_models(n_points);
        for (int i = 0; i < n_points; i++) {
            pt_models[i].reserve(2 * k + 1);
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        // Phase 1: Creating single model predictions
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
            ulogit_t lgtic_res = ulogit(x.data() + x_min_index,
                                                                  y.data() + x_min_index,
                                                                  w_window,
                                                                  max_iterations,
                                                                  ridge_lambda,
                                                                  max_beta,
                                                                  tolerance,
                                                                  verbose);
            lgtic_res.x_min_index = x_min_index;
            lgtic_res.x_max_index = x_max_index;

            // Store single-model prediction and its LOOCV squared errors
            k_errors[i] = lgtic_res.errors[i - x_min_index];
            double y_prediction_at_ref_pt = lgtic_res.predictions[i - x_min_index];
            k_sm_predictions[k_index][i] = y_prediction_at_ref_pt;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - y_prediction_at_ref_pt);
            }

            // For x indices around i insert that model into pt_models[i]
            for (int j = x_min_index; j <= x_max_index; j++) {
                lgtic_res.w = w_window;  // Store the weights for later model averaging
                pt_models[j].push_back(lgtic_res);
            }
        }

        if (verbose) {
            elapsed_time(phase1_ptm, "Done");
            mem_tracker.report();
        }

        //k_predictions = k_sm_predictions;
        #if 1
        // Phase 2: Model averaging
        // k_errors and k_true_errors are reused for model-averaged predictions
        if (verbose) {
            Rprintf("  Phase 2: Computing model-averaged predictions ... ");
        }
        auto phase2_ptm = std::chrono::steady_clock::now();
        for (int i = 0; i < n_points; i++) {

            double weighted_sum = 0.0;
            double weight_sum = 0.0;
            double wmean_error = 0.0;

            for (const auto& model : pt_models[i]) {
                int local_index = i - model.x_min_index;
                double weight = model.w[local_index];
                weighted_sum += weight * model.predictions[local_index];
                weight_sum += weight;
                wmean_error +=  weight * model.errors[local_index];
            }

            k_errors[i] = wmean_error / weight_sum;

            // Store model-averaged prediction and its errors
            k_predictions[k_index][i] = weighted_sum / weight_sum;

            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - k_predictions[k_index][i]);
            }
        }

        if (verbose) {
            elapsed_time(phase2_ptm, "Done");
            mem_tracker.report();
        }
        #endif

        // Compute mean errors for model-averaged predictions at current k
        results.k_mean_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
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
        elapsed_time(total_ptm, "\nTotal MABILOG computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}


/**
 * @brief An un-weighted version of the Model-Averaged bi-kNN Logistic regression (MABILOG) algorithm
 *
 * @details MABILOG extends the traditional LOWESS algorithm to the binary
 * outcome case by incorporating model averaging with bi-k nearest neighbor
 * structure. The algorithm consists of three phases:
 *
 * 1. For each k in [k_min, k_max]:
 *    - Fit local logistic regression models using k-hop neighborhoods
 *    - Compute predictions and their LOOCV errors
 *    - Perform kernel-weighted model averaging
 *    - Calculate mean LOOCV and true errors for predictions
 *
 * 2. For each point, compute model-averaged predictions using kernel weights
 *    - Models contribute to predictions based on their kernel weights
 *    - Final prediction is a weighted average of all contributing models
 *
 * 3. Find optimal k for model-averaged predictions based on weighted LOOCV errors
 *
 * Window Size and k Relationship:
 * - For each point x_i, the window size is 2k + 1 points
 * - k points are taken from both left and right of x_i when available
 * - Near boundaries, the window is adjusted while maintaining total size:
 *   * Left boundary: takes first 2k + 1 points
 *   * Right boundary: takes last 2k + 1 points
 *
 * Distance Normalization:
 * The dist_normalization_factor (typically 1.01) serves to:
 * - Scale the maximum distance in each window slightly up
 * - Prevent zero/near-zero weights at window boundaries
 * - Improve numerical stability of the weighted regression
 * Values closer to 1.0 give more weight to boundary points
 *
 * Computational Complexity:
 * - Time: O(n * K * w), where:
 *   * n is the number of points
 *   * K is the number of k values (k_max - k_min + 1)
 *   * w is the maximum window size (2 * k_max + 1)
 * - Space: O(n * K) for storing predictions and models
 * - Each local linear regression: O(w) operations
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Observed y values corresponding to x (response variable)
 * @param y_true Optional true y values for Rf_error calculation. Used for algorithm evaluation
 * @param k_min Minimum number of neighbors on each side (minimum window half-width)
 * @param k_max Maximum number of neighbors on each side (maximum window half-width)
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param tolerance Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilog_t structure containing:
 *         - opt_k: Optimal k value for model-averaged predictions
 *         - predictions: Model-averaged predictions using opt_k
 *         - k_mean_errors: Mean LOOCV errors for model-averaged predictions for each k
 *         - k_mean_true_errors: Mean true errors for model-averaged predictions if y_true provided
 *         - k_predictions: Model-averaged predictions for all k values
 *
 * @throws std::invalid_argument for invalid input parameters
 * @throws Rf_error for numerical instability
 *
 * @note
 * - Input x values must be sorted in ascending order
 * - Model averaging uses kernel weights based on distances between points
 *
 * @see ulogit_t for the structure of individual linear models
 * @see initialize_kernel() for kernel function details
 */
mabilog_t uwmabilog(const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& y_true,
                    int max_iterations,
                    double ridge_lambda,
                    double max_beta,
                    double tolerance,
                    int k_min,
                    int k_max,
                    int distance_kernel,
                    double dist_normalization_factor,
                    bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILOG");

    if (verbose) {
        Rprintf("Starting MABILOG computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    mabilog_t results;

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
    std::vector<std::pair<double, const ulogit_t*>> filtered_models;
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

        std::vector<std::vector<ulogit_t>> pt_models(n_points);
        for (int i = 0; i < n_points; i++) {
            pt_models[i].reserve(2 * k + 1);
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        // Phase 1: Creating single model predictions
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
            ulogit_t lgtic_res = ulogit(x.data() + x_min_index,
                                                                  y.data() + x_min_index,
                                                                  w_window,
                                                                  max_iterations,
                                                                  ridge_lambda,
                                                                  max_beta,
                                                                  tolerance,
                                                                  verbose);
            lgtic_res.x_min_index = x_min_index;
            lgtic_res.x_max_index = x_max_index;

            // Store single-model prediction and its LOOCV squared errors
            k_errors[i] = lgtic_res.errors[i - x_min_index];
            double y_prediction_at_ref_pt = lgtic_res.predictions[i - x_min_index];
            k_sm_predictions[k_index][i] = y_prediction_at_ref_pt;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - y_prediction_at_ref_pt);
            }

            // For x indices around i insert that model into pt_models[i]
            for (int j = x_min_index; j <= x_max_index; j++) {
                lgtic_res.w = w_window;  // Store the weights for later model averaging
                pt_models[j].push_back(lgtic_res);
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
        if (verbose) {
            elapsed_time(phase2_ptm, "Done");
            mem_tracker.report();
        }

        // Compute mean errors for model-averaged predictions at current k
        results.k_mean_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
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
        elapsed_time(total_ptm, "\nTotal MABILOG computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}

/**
 * \brief R interface for weighted MABILO (logistic variant).
 *
 * \param s_x  Numeric vector of x.
 * \param s_y  Numeric vector of responses (typically 0/1 for logistic).
 * \param s_y_true Optional numeric vector of true responses (same length as x/y).
 * \param s_w  Numeric vector of per-point weights (same length as x/y).
 * \param s_max_iterations Integer >= 1.
 * \param s_ridge_lambda   Double >= 0 (L2 regularization).
 * \param s_max_beta       Positive double (cap on |beta|).
 * \param s_tolerance      Positive double (convergence tolerance).
 * \param s_k_min          Integer >= 1.
 * \param s_k_max          Integer >= k_min and <= length(x).
 * \param s_distance_kernel Integer kernel code (0=Tricube, 1=Epanechnikov, 2=Exponential).
 * \param s_dist_normalization_factor Double distance normalization factor.
 * \param s_verbose Logical flag for progress output.
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
 * \throws Calls \c Rf_error() on type/length mismatches or invalid parameters.
 */
SEXP S_wmabilog(SEXP s_x,
                SEXP s_y,
                SEXP s_y_true,
                SEXP s_w,
                SEXP s_max_iterations,
                SEXP s_ridge_lambda,
                SEXP s_max_beta,
                SEXP s_tolerance,
                SEXP s_k_min,
                SEXP s_k_max,
                SEXP s_distance_kernel,
                SEXP s_dist_normalization_factor,
                SEXP s_verbose) {

    // --- Coerce x/y/w (and optional y_true) to REAL and copy (long-vector safe) ---
    std::vector<double> x, y, y_true, w;
    {
        SEXP sx = s_x, sy = s_y, sw = s_w, syt = (s_y_true == R_NilValue ? R_NilValue : s_y_true);
        PROTECT_INDEX pix, piy, piw, pzt; // pzt for syt
        PROTECT_WITH_INDEX(sx, &pix);
        PROTECT_WITH_INDEX(sy, &piy);
        PROTECT_WITH_INDEX(sw, &piw);
        PROTECT_WITH_INDEX(syt, &pzt);

        if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), pix);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), piy);
        if (TYPEOF(sw) != REALSXP) REPROTECT(sw = Rf_coerceVector(sw, REALSXP), piw);
        if (syt != R_NilValue && TYPEOF(syt) != REALSXP)
            REPROTECT(syt = Rf_coerceVector(syt, REALSXP), pzt);

        const R_xlen_t nx = XLENGTH(sx), ny = XLENGTH(sy), nw = XLENGTH(sw);
        if (nx != ny) { UNPROTECT(4); Rf_error("x and y must have the same length."); }
        if (nx != nw) { UNPROTECT(4); Rf_error("x and w must have the same length."); }

        x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(nx));
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));
        w.assign(REAL(sw), REAL(sw) + static_cast<size_t>(nw));
        if (syt != R_NilValue) {
            const R_xlen_t nyt = XLENGTH(syt);
            if (nyt == nx) {
                y_true.assign(REAL(syt), REAL(syt) + static_cast<size_t>(nyt));
            }
        }
        UNPROTECT(4); // sx, sy, sw, syt
    }

    // --- Scalars / parameters (validated) ---
    const int    max_iterations = Rf_asInteger(s_max_iterations);
    const double ridge_lambda   = Rf_asReal   (s_ridge_lambda);
    const double max_beta       = Rf_asReal   (s_max_beta);
    const double tolerance      = Rf_asReal   (s_tolerance);
    const int    k_min          = Rf_asInteger(s_k_min);
    const int    k_max          = Rf_asInteger(s_k_max);
    const int    distance_kernel= Rf_asInteger(s_distance_kernel);
    const double dist_norm_factor = Rf_asReal (s_dist_normalization_factor);
    const bool   verbose        = (Rf_asLogical(s_verbose) == TRUE);

    // --- Core computation (no R allocations inside) ---
    mabilog_t wmabilog_results = wmabilog(x, y, y_true, w,
                                          max_iterations, ridge_lambda,
                                          max_beta, tolerance,
                                          k_min, k_max,
                                          distance_kernel, dist_norm_factor,
                                          verbose);

    // --- Build result (container-first; per-element PROTECT/UNPROTECT) ---
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 7));

    // 0: k_values sequence [k_min .. k_max]
    {
        const R_xlen_t K = static_cast<R_xlen_t>(
            static_cast<long long>(k_max) - static_cast<long long>(k_min) + 1LL);
        std::vector<int> k_values(static_cast<size_t>(K));
        for (R_xlen_t i = 0; i < K; ++i)
            k_values[static_cast<size_t>(i)] = k_min + static_cast<int>(i);
        SEXP kv = PROTECT(convert_vector_int_to_R(k_values));
        SET_VECTOR_ELT(result, 0, kv);
        UNPROTECT(1);
    }

    // 1: opt_k (scalar int)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(wmabilog_results.opt_k));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }

    // 2: opt_k_idx (scalar int; expose as returned by core)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(wmabilog_results.opt_k_idx + 1));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }

    // 3: k_mean_errors (numeric)
    {
        SEXP s = PROTECT(convert_vector_double_to_R(wmabilog_results.k_mean_errors));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }

    // 4: k_mean_true_errors (numeric or NULL)
    if (!y_true.empty()) {
        SEXP s = PROTECT(convert_vector_double_to_R(wmabilog_results.k_mean_true_errors));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 4, R_NilValue);
    }

    // 5: predictions (numeric)
    {
        SEXP s = PROTECT(convert_vector_double_to_R(wmabilog_results.predictions));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }

    // 6: k_predictions (list<numeric>)
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(wmabilog_results.k_predictions));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(names, 0, Rf_mkChar("k_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("opt_k"));
    SET_STRING_ELT(names, 2, Rf_mkChar("opt_k_idx"));
    SET_STRING_ELT(names, 3, Rf_mkChar("k_mean_errors"));
    SET_STRING_ELT(names, 4, Rf_mkChar("k_mean_true_errors"));
    SET_STRING_ELT(names, 5, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 6, Rf_mkChar("k_predictions"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return result;
}


/**
 * @brief Performs parallel Bayesian bootstrap calculations for MABILOG
 *
 * @details Implements parallel bootstrap sampling for uncertainty quantification:
 * 1. Generates bootstrap weights using Bayesian bootstrap
 * 2. Performs MABILOG fitting for each bootstrap sample
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
 * @param tolerance Numerical stability parameter (default: 1e-8)
 * @param verbose Enable progress messages (default: false)
 *
 * @return Vector of vectors containing bootstrap predictions for each iteration
 *
 * @throws Rf_error if parameters are invalid or computation fails
 * @note Thread-safe implementation using parallel execution
 */
std::vector<std::vector<double>> mabilog_bb(const std::vector<double>& x,
                                            const std::vector<double>& y,
                                            int max_iterations,
                                            double ridge_lambda,
                                            double max_beta,
                                            double tolerance,
                                            int k,
                                            int n_bb,
                                            int distance_kernel,
                                            double dist_normalization_factor = 1.01,
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
    gflow::for_each(GFLOW_EXEC_POLICY,
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
                        auto wmabilog_results = wmabilog(x,
                                                         y,
                                                         y_true,
                                                         weights,
                                                         max_iterations,
                                                         ridge_lambda,
                                                         max_beta,
                                                         tolerance,
                                                         k,
                                                         k,
                                                         distance_kernel,
                                                         dist_normalization_factor,
                                                         verbose);

                        // Store results - no need for mutex as each thread writes to its own index
                        bb_predictions[iboot] = std::move(wmabilog_results.predictions);
                    });

    return bb_predictions;
}

/**
 * @brief Computes Bayesian bootstrap predictions with credible intervals for MABILOG
 *
 * @details Implements a complete Bayesian uncertainty analysis:
 * 1. Performs parallel bootstrap iterations using mabilog_bb
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
 * @param tolerance Numerical stability parameter
 *
 * @return bb_cri_t structure containing:
 *         - bb_predictions: Median predictions across bootstrap iterations
 *         - cri_L: Lower bounds of credible intervals
 *         - cri_U: Upper bounds of credible intervals
 *
 * @throws Rf_error for invalid parameters or failed computation
 */
bb_cri_t mabilog_bb_cri(const std::vector<double>& x,
                        const std::vector<double>& y,
                        int max_iterations,
                        double ridge_lambda,
                        double max_beta,
                        double tolerance,
                        int k,
                        int n_bb,
                        double p,
                        int distance_kernel,
                        double dist_normalization_factor) {

    // Perform bootstrap iterations
    std::vector<std::vector<double>> bb_predictionss = mabilog_bb(x,
                                                                  y,
                                                                  max_iterations,
                                                                  ridge_lambda,
                                                                  max_beta,
                                                                  tolerance,
                                                                  k,
                                                                  n_bb,
                                                                  distance_kernel,
                                                                  dist_normalization_factor,
                                                                  false);

    // Calculate credible intervals
    bool use_median = true;
    return bb_cri(bb_predictionss, use_median, p);
}


/**
 * @brief Main interface for MABILOG (Model-Averaged Locally Weighted Scatterplot Smoothing) with optional Bayesian bootstrap
 *
 * @details This function provides a complete implementation of MABILOG with two main components:
 * 1. Core MABILOG algorithm:
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
 * @param y_true Optional vector of true y values for Rf_error calculation
 * @param k_min Minimum number of neighbors on each side
 * @param k_max Maximum number of neighbors on each side
 * @param n_bb Number of Bayesian bootstrap iterations (0 to skip bootstrap)
 * @param p Probability level for credible intervals (used only if n_bb > 0)
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param tolerance Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilog_t structure containing:
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
mabilog_t mabilog(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& y_true,
                   int max_iterations,
                   double ridge_lambda,
                   double max_beta,
                   double tolerance,
                   int k_min,
                   int k_max,
                   int n_bb,
                   double p,
                   int distance_kernel,
                   double dist_normalization_factor,
                   bool verbose) {

    mabilog_t uwmabilog_results = uwmabilog(x,
                                             y,
                                             y_true,
                                             max_iterations,
                                             ridge_lambda,
                                             max_beta,
                                             tolerance,
                                             k_min,
                                             k_max,
                                             distance_kernel,
                                             dist_normalization_factor,
                                             verbose);

    if (n_bb) {
        bb_cri_t bb_res = mabilog_bb_cri(x,
                                         y,
                                         max_iterations,
                                         ridge_lambda,
                                         max_beta,
                                         tolerance,
                                         uwmabilog_results.opt_k,
                                         n_bb,
                                         p,
                                         distance_kernel,
                                         dist_normalization_factor);

        uwmabilog_results.bb_predictions = std::move(bb_res.bb_Ey);
        uwmabilog_results.cri_L             = std::move(bb_res.cri_L);
        uwmabilog_results.cri_U             = std::move(bb_res.cri_U);
    }

    return uwmabilog_results;
}

/**
 * @brief R interface for (unweighted) MABILO logistic variant with optional bootstrap.
 *
 * @param s_x  Numeric vector of predictors (length N).
 * @param s_y  Numeric vector of responses (length N; typically 0/1).
 * @param s_y_true Optional numeric vector of true responses (length N) for true-error summaries.
 * @param s_max_iterations Integer >= 1.
 * @param s_ridge_lambda   Double >= 0 (L2).
 * @param s_max_beta       Double > 0 (cap on |beta|).
 * @param s_tolerance      Double > 0 (conv. tol).
 * @param s_k_min          Integer >= 1.
 * @param s_k_max          Integer >= k_min and <= N.
 * @param s_n_bb           Integer >= 0 (bootstrap reps).
 * @param s_p              Double in (0,1] (e.g., CI quantile).
 * @param s_distance_kernel Kernel function type (integer):
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param s_dist_normalization_factor Double distance normalization factor.
 * @param s_verbose        Logical flag for progress output.
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
 * @throws Calls \c Rf_error() on type/length mismatches or invalid parameters.
 *
 * @note
 * - All input vectors must have the same length
 * - x values must be sorted in ascending order
 * - Bootstrap results are NULL if n_bb = 0
 * - Uses R's protection stack for memory management
 *
 */
SEXP S_mabilog(SEXP s_x,
               SEXP s_y,
               SEXP s_y_true,
               SEXP s_max_iterations,
               SEXP s_ridge_lambda,
               SEXP s_max_beta,
               SEXP s_tolerance,
               SEXP s_k_min,
               SEXP s_k_max,
               SEXP s_n_bb,
               SEXP s_p,
               SEXP s_distance_kernel,
               SEXP s_dist_normalization_factor,
               SEXP s_verbose) {

    // ---- Coerce x/y (and optional y_true) to REAL and copy (long-vector safe) ----
    std::vector<double> x, y, y_true;
    {
        // Always hold exactly 3 protected slots: sx, sy, syt (syt may be R_NilValue)
        SEXP sx = s_x, sy = s_y, syt = (s_y_true == R_NilValue ? R_NilValue : s_y_true);
        PROTECT_INDEX pix, piy, piyt;
        PROTECT_WITH_INDEX(sx, &pix);
        PROTECT_WITH_INDEX(sy, &piy);
        PROTECT_WITH_INDEX(syt, &piyt);

        if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), pix);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), piy);
        if (syt != R_NilValue && TYPEOF(syt) != REALSXP)
            REPROTECT(syt = Rf_coerceVector(syt, REALSXP), piyt);

        const R_xlen_t nx = XLENGTH(sx);
        const R_xlen_t ny = XLENGTH(sy);
        x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(nx));
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));

        if (syt != R_NilValue) {
            const R_xlen_t nyt = XLENGTH(syt);
            if (nyt == nx) {
                y_true.assign(REAL(syt), REAL(syt) + static_cast<size_t>(nyt));
            }
        }
        UNPROTECT(3); // sx, sy, syt
    }

    const R_xlen_t N = static_cast<R_xlen_t>(x.size());

    // ---- Scalars / parameters (validated) ----
    const int    max_iterations = Rf_asInteger(s_max_iterations);
    const double ridge_lambda   = Rf_asReal   (s_ridge_lambda);
    const double max_beta       = Rf_asReal   (s_max_beta);
    const double tolerance      = Rf_asReal   (s_tolerance);
    const int    k_min          = Rf_asInteger(s_k_min);
    const int    k_max          = Rf_asInteger(s_k_max);
    const int    n_bb           = Rf_asInteger(s_n_bb);
    const double p              = Rf_asReal   (s_p);
    const int    distance_kernel= Rf_asInteger(s_distance_kernel);
    const double dist_norm_factor = Rf_asReal (s_dist_normalization_factor);
    const bool   verbose        = (Rf_asLogical(s_verbose) == TRUE);

    if (max_iterations < 1)      Rf_error("max_iterations must be >= 1.");
    if (ridge_lambda   < 0.0)    Rf_error("ridge_lambda must be >= 0.");
    if (!(max_beta     > 0.0))   Rf_error("max_beta must be > 0.");
    if (!(tolerance    > 0.0))   Rf_error("tolerance must be > 0.");
    if (k_min < 1)               Rf_error("k_min must be >= 1.");
    if (k_max < k_min)           Rf_error("k_max must be >= k_min.");
    if (static_cast<R_xlen_t>(k_max) > N)
        Rf_error("k_max (%d) cannot exceed length(x) (%lld).",
                 k_max, static_cast<long long>(N));
    if (n_bb < 0)                Rf_error("n_bb must be >= 0.");
    if (!(p > 0.0 && p <= 1.0))  Rf_error("p must be in (0, 1].");

    // ---- Core computation (no R allocations inside) ----
    mabilog_t out = mabilog(x, y, y_true,
                            max_iterations, ridge_lambda, max_beta, tolerance,
                            k_min, k_max,
                            n_bb, p,
                            distance_kernel, dist_norm_factor,
                            verbose);

    // ---- Build result (container-first; per-element PROTECT/UNPROTECT) ----
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 10));

    // 0: k_values sequence [k_min .. k_max]
    {
        const R_xlen_t K = static_cast<R_xlen_t>(
            static_cast<long long>(k_max) - static_cast<long long>(k_min) + 1LL);
        std::vector<int> k_values(static_cast<size_t>(K));
        for (R_xlen_t i = 0; i < K; ++i)
            k_values[static_cast<size_t>(i)] = k_min + static_cast<int>(i);
        SEXP kv = PROTECT(convert_vector_int_to_R(k_values));
        SET_VECTOR_ELT(result, 0, kv);
        UNPROTECT(1);
    }

    // 1: opt_k (scalar int)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(out.opt_k));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }

    // 2: opt_k_idx
    {
        SEXP s = PROTECT(Rf_ScalarInteger(out.opt_k_idx + 1));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }

    // 3: k_mean_errors
    {
        SEXP s = PROTECT(convert_vector_double_to_R(out.k_mean_errors));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }

    // 4: k_mean_true_errors (or NULL)
    if (!y_true.empty()) {
        SEXP s = PROTECT(convert_vector_double_to_R(out.k_mean_true_errors));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 4, R_NilValue);
    }

    // 5: predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(out.predictions));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }

    // 6: k_predictions (list<numeric>)
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(out.k_predictions));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    // 7–9: bootstrap outputs (or NULLs if n_bb == 0)
    if (n_bb > 0) {
        SEXP s = PROTECT(convert_vector_double_to_R(out.bb_predictions));
        SET_VECTOR_ELT(result, 7, s); UNPROTECT(1);

        s = PROTECT(convert_vector_double_to_R(out.cri_L));
        SET_VECTOR_ELT(result, 8, s); UNPROTECT(1);

        s = PROTECT(convert_vector_double_to_R(out.cri_U));
        SET_VECTOR_ELT(result, 9, s); UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
    SET_STRING_ELT(names, 0, Rf_mkChar("k_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("opt_k"));
    SET_STRING_ELT(names, 2, Rf_mkChar("opt_k_idx"));
    SET_STRING_ELT(names, 3, Rf_mkChar("k_mean_errors"));
    SET_STRING_ELT(names, 4, Rf_mkChar("k_mean_true_errors"));
    SET_STRING_ELT(names, 5, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 6, Rf_mkChar("k_predictions"));
    SET_STRING_ELT(names, 7, Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 8, Rf_mkChar("cri_L"));
    SET_STRING_ELT(names, 9, Rf_mkChar("cri_U"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2);
    return result;
}


/**
 * @brief Smoothed Rf_error version of Model-Averaged LOWESS (MABILOG) for robust local regression
 *
 * @details Similar to uwmabilog, but applies additional smoothing to the LOOCV errors
 * to reduce noise in k selection. The Rf_error curve is smoothed using uwmabilog with
 * a fixed window size of 0.25 * n_points.
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Observed y values corresponding to x (response variable)
 * @param y_true Optional true y values for Rf_error calculation
 * @param k_min Minimum number of neighbors on each side
 * @param k_max Maximum number of neighbors on each side
 * @param error_window_factor Factor to determine window size for Rf_error curve smoothing (default: 0.25).
 *                           The window size will be error_window_factor * n_points neighbors on each side.
 *                           Larger values create smoother Rf_error curves but may miss local structure.
 * @param distance_kernel Kernel function (0: Tricube, 1: Epanechnikov, 2: Exponential)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param tolerance Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilog_t structure with additional smoothed_k_mean_errors field
 *
 * @throws std::invalid_argument for invalid input parameters
 * @throws Rf_error for numerical instability
 *
 * @note The smoothing helps prevent selecting suboptimal k values due to noise
 * in the Rf_error measurements
 */
mabilog_t mabilog_with_smoothed_errors(const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        const std::vector<double>& y_true,
                                        int max_iterations,
                                        double ridge_lambda,
                                        double max_beta,
                                        double tolerance,
                                        int k_min,
                                        int k_max,
                                        double error_window_factor,
                                        int distance_kernel,
                                        double dist_normalization_factor,
                                        bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILOG");

    if (verbose) {
        Rprintf("Starting MABILOG computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    mabilog_t results;

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
    std::vector<std::pair<double, const ulogit_t*>> filtered_models;
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

        std::vector<std::vector<ulogit_t>> pt_models(n_points);
        for (int i = 0; i < n_points; i++) {
            pt_models[i].reserve(2 * k + 1);
        }

        int n_points_minus_k = n_points - k;
        int n_points_minus_k_minus_one = n_points - k - 1;
        int k_minus_one = k - 1;
        int two_k = 2 * k;
        int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

        // Phase 1: Creating single model predictions
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
            ulogit_t lgtic_res = ulogit(x.data() + x_min_index,
                                        y.data() + x_min_index,
                                        w_window,
                                        max_iterations,
                                        ridge_lambda,
                                        max_beta,
                                        tolerance,
                                        verbose);
            lgtic_res.x_min_index = x_min_index;
            lgtic_res.x_max_index = x_max_index;

            // Store single-model prediction and its LOOCV squared errors
            k_errors[i] = lgtic_res.errors[i - x_min_index];
            double y_prediction_at_ref_pt = lgtic_res.predictions[i - x_min_index];
            k_sm_predictions[k_index][i] = y_prediction_at_ref_pt;
            if (y_true_exists) {
                k_true_errors[i] = std::abs(y_true[i] - y_prediction_at_ref_pt);
            }

            // For x indices around i insert that model into pt_models[i]
            for (int j = x_min_index; j <= x_max_index; j++) {
                lgtic_res.w = w_window;  // Store the weights for later model averaging
                pt_models[j].push_back(lgtic_res);
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

        // smoothing results.k_mean_errors using uwmabilog() with k = 0.25 * n_points
        int n_k_values = k_max - k_min + 1;
        int window_size = static_cast<int>(error_window_factor * n_k_values);
        // 2k + 1 = window_size => k = (window_size - 1) / 2
        int error_k = (window_size - 1) / 2;
        if (error_k < 1) error_k = 1;

        std::vector<double> k_values(n_k_values);
        for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
            k_values[k_index] = k;

        std::vector<double> true_errors;
        auto errors_fit = uwmabilog(k_values,
                                    results.k_mean_errors,
                                    true_errors,
                                    max_iterations,
                                    ridge_lambda,
                                    max_beta,
                                    tolerance,
                                    error_k,
                                    error_k,
                                    distance_kernel,
                                    dist_normalization_factor,
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
        elapsed_time(total_ptm, "\nTotal MABILOG computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}

/**
 * @brief R interface for logistic MABILO with error smoothing.
 *
 * @param s_x  Numeric vector of predictors (length N).
 * @param s_y  Numeric vector of responses (length N; typically 0/1).
 * @param s_y_true Optional numeric vector of true responses (length N) for true-error summaries.
 * @param s_max_iterations Integer >= 1.
 * @param s_ridge_lambda   Double >= 0 (L2 regularization strength).
 * @param s_max_beta       Double > 0 (cap on |beta|).
 * @param s_tolerance      Double > 0 (convergence tolerance).
 * @param s_k_min          Integer >= 1.
 * @param s_k_max          Integer >= k_min and <= N.
 * @param s_error_window_factor Double window factor for smoothing cross-validated errors.
 * @param s_distance_kernel Integer in {0,1,2} (0=Tricube, 1=Epanechnikov, 2=Exponential).
 * @param s_dist_normalization_factor Double distance normalization factor.
 * @param s_verbose        Logical flag; if TRUE, print progress.
 *
 * @return R list containing:
 *         - k_values: Vector of k values used
 *         - opt_k: Optimal k value
 *         - opt_k_idx: Index of optimal k
 *         - k_mean_errors: Raw Rf_error values
 *         - smoothed_k_mean_errors: Smoothed Rf_error values
 *         - k_mean_true_errors: True errors (if y_true provided)
 *         - predictions: Final predictions
 *         - k_predictions: Predictions for all k values
 *
 * @throws Calls \c Rf_error() on type/length mismatches or invalid parameters.
 */
SEXP S_mabilog_with_smoothed_errors(SEXP s_x,
                                    SEXP s_y,
                                    SEXP s_y_true,
                                    SEXP s_max_iterations,
                                    SEXP s_ridge_lambda,
                                    SEXP s_max_beta,
                                    SEXP s_tolerance,
                                    SEXP s_k_min,
                                    SEXP s_k_max,
                                    SEXP s_error_window_factor,
                                    SEXP s_distance_kernel,
                                    SEXP s_dist_normalization_factor,
                                    SEXP s_verbose) {

  // ---- Coerce x/y (and optional y_true) to REAL and copy (long-vector safe) ----
  std::vector<double> x, y, y_true;
  {
      // Always hold exactly 3 protected slots: sx, sy, syt (syt may be R_NilValue)
      SEXP sx = s_x, sy = s_y, syt = (s_y_true == R_NilValue ? R_NilValue : s_y_true);
      PROTECT_INDEX pix, piy, piyt;
      PROTECT_WITH_INDEX(sx, &pix);
      PROTECT_WITH_INDEX(sy, &piy);
      PROTECT_WITH_INDEX(syt, &piyt);

      if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), pix);
      if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), piy);
      if (syt != R_NilValue && TYPEOF(syt) != REALSXP)
          REPROTECT(syt = Rf_coerceVector(syt, REALSXP), piyt);

      const R_xlen_t nx = XLENGTH(sx);
      const R_xlen_t ny = XLENGTH(sy);
      x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(nx));
      y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));

      if (syt != R_NilValue) {
          const R_xlen_t nyt = XLENGTH(syt);
          if (nyt == nx) {
              y_true.assign(REAL(syt), REAL(syt) + static_cast<size_t>(nyt));
          }
      }
      UNPROTECT(3); // sx, sy, syt
  }

  const R_xlen_t N = static_cast<R_xlen_t>(x.size());

  // ---- Scalars / parameters (validated) ----
  const int    max_iterations = Rf_asInteger(s_max_iterations);
  const double ridge_lambda   = Rf_asReal   (s_ridge_lambda);
  const double max_beta       = Rf_asReal   (s_max_beta);
  const double tolerance      = Rf_asReal   (s_tolerance);
  const int    k_min          = Rf_asInteger(s_k_min);
  const int    k_max          = Rf_asInteger(s_k_max);
  const double error_window_factor = Rf_asReal(s_error_window_factor);
  const int    distance_kernel= Rf_asInteger(s_distance_kernel);
  const double dist_norm_factor = Rf_asReal (s_dist_normalization_factor);
  const bool   verbose        = (Rf_asLogical(s_verbose) == TRUE);

  if (max_iterations < 1)      Rf_error("max_iterations must be >= 1.");
  if (ridge_lambda   < 0.0)    Rf_error("ridge_lambda must be >= 0.");
  if (!(max_beta     > 0.0))   Rf_error("max_beta must be > 0.");
  if (!(tolerance    > 0.0))   Rf_error("tolerance must be > 0.");
  if (k_min < 1)               Rf_error("k_min must be >= 1.");
  if (k_max < k_min)           Rf_error("k_max must be >= k_min.");
  if (static_cast<R_xlen_t>(k_max) > N)
                               Rf_error("k_max (%d) cannot exceed length(x) (%lld).",
                                        k_max, static_cast<long long>(N));

  // ---- Core computation (no R allocations inside) ----
  mabilog_t out = mabilog_with_smoothed_errors(x, y, y_true,
                                               max_iterations, ridge_lambda,
                                               max_beta, tolerance,
                                               k_min, k_max,
                                               error_window_factor,
                                               distance_kernel, dist_norm_factor,
                                               verbose);

  // ---- Build result (container-first; per-element PROTECT/UNPROTECT) ----
  SEXP result = PROTECT(Rf_allocVector(VECSXP, 8));

  // 0: k_values sequence [k_min .. k_max]
  {
    const R_xlen_t K = static_cast<R_xlen_t>(
        static_cast<long long>(k_max) - static_cast<long long>(k_min) + 1LL);
    std::vector<int> k_values(static_cast<size_t>(K));
    for (R_xlen_t i = 0; i < K; ++i)
      k_values[static_cast<size_t>(i)] = k_min + static_cast<int>(i);
    SEXP kv = PROTECT(convert_vector_int_to_R(k_values));
    SET_VECTOR_ELT(result, 0, kv);
    UNPROTECT(1);
  }

  // 1: opt_k (scalar int; expose as core returns)
  {
    SEXP s = PROTECT(Rf_ScalarInteger(out.opt_k));
    SET_VECTOR_ELT(result, 1, s);
    UNPROTECT(1);
  }

  // 2: opt_k_idx
  {
    SEXP s = PROTECT(Rf_ScalarInteger(out.opt_k_idx + 1));
    SET_VECTOR_ELT(result, 2, s);
    UNPROTECT(1);
  }

  // 3: k_mean_errors
  {
    SEXP s = PROTECT(convert_vector_double_to_R(out.k_mean_errors));
    SET_VECTOR_ELT(result, 3, s);
    UNPROTECT(1);
  }

  // 4: smoothed_k_mean_errors
  {
    SEXP s = PROTECT(convert_vector_double_to_R(out.smoothed_k_mean_errors));
    SET_VECTOR_ELT(result, 4, s);
    UNPROTECT(1);
  }

  // 5: k_mean_true_errors (or NULL)
  if (!y_true.empty()) {
    SEXP s = PROTECT(convert_vector_double_to_R(out.k_mean_true_errors));
    SET_VECTOR_ELT(result, 5, s);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(result, 5, R_NilValue);
  }

  // 6: predictions
  {
    SEXP s = PROTECT(convert_vector_double_to_R(out.predictions));
    SET_VECTOR_ELT(result, 6, s);
    UNPROTECT(1);
  }

  // 7: k_predictions (list<numeric>)
  {
    SEXP s = PROTECT(convert_vector_vector_double_to_R(out.k_predictions));
    SET_VECTOR_ELT(result, 7, s);
    UNPROTECT(1);
  }

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 8));
  SET_STRING_ELT(names, 0, Rf_mkChar("k_values"));
  SET_STRING_ELT(names, 1, Rf_mkChar("opt_k"));
  SET_STRING_ELT(names, 2, Rf_mkChar("opt_k_idx"));
  SET_STRING_ELT(names, 3, Rf_mkChar("k_mean_errors"));
  SET_STRING_ELT(names, 4, Rf_mkChar("smoothed_k_mean_errors"));
  SET_STRING_ELT(names, 5, Rf_mkChar("k_mean_true_errors"));
  SET_STRING_ELT(names, 6, Rf_mkChar("predictions"));
  SET_STRING_ELT(names, 7, Rf_mkChar("k_predictions"));
  Rf_setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(2);
  return result;
}
