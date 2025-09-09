#include "A_mabilo.h"

/**
 * @brief Implements the Model-Averaged Bi-kNN LOcal linear model (MABILO) algorithm for robust local regression
 *
 * @details MABILO extends the traditional LOWESS algorithm by incorporating model averaging with
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
 * @param y_true Optional true y values for error calculation. Used for algorithm evaluation
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
 * @see ulm_plus_t for the structure of individual linear models
 * @see initialize_kernel() for kernel function details
 */
mabilo_t old_wmabilo(const std::vector<double>& x,
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

    mabilo_t results;

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

#if 0
        struct pred_w_err_t {
        double prediction;
        double weight;
        double error;
    };

        std::vector<std::vector<pred_w_err_t>> pt_pred_w_err(n_points);
#endif

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
            ulm_plus_t wlm_res = ulm(x.data() + x_min_index,
                                y.data() + x_min_index,
                                w_window,
                                y_binary,
                                epsilon);
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
 * @brief An un-weighted version of the Model-Averaged LOWESS (MABILO) algorithm for robust local regression
 *
 * @details MABILO extends the traditional LOWESS algorithm by incorporating model averaging with
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
 * @param y_true Optional true y values for error calculation. Used for algorithm evaluation
 * @param k_min Minimum number of neighbors on each side (minimum window half-width)
 * @param k_max Maximum number of neighbors on each side (maximum window half-width)
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
 * @see ulm_t for the structure of individual linear models
 * @see initialize_kernel() for kernel function details
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
    mabilo_t results;
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
