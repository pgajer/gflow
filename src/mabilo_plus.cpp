// Experimental version of MABILO with different model kernels and Rf_error filtering strategies

#include "ulm.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "pglm.h"     // for itriplet_t
#include "kernels.h"  // for initialize_kernel()
#include "memory_utils.hpp"
#include "progress_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

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
    SEXP S_mabilo_plus(SEXP s_x,
                       SEXP s_y,
                       SEXP s_y_true,
                       SEXP s_w,
                       SEXP s_k_min,
                       SEXP s_k_max,
                       SEXP s_model_averaging_strategy,
                       SEXP s_error_filtering_strategy,
                       SEXP s_distance_kernel,
                       SEXP s_model_kernel,
                       //SEXP s_run_parallel,
                       SEXP s_dist_normalization_factor,
                       SEXP s_epsilon,
                       SEXP s_verbose);
}

struct mabilo_plus_t {
    // k values
    int opt_sm_k;     // optimal single model k value - the one with the smallest mean LOOCV Rf_error
    int opt_sm_k_idx; // optimal single model k value index
    int opt_ma_k;     // optimal model averaging k value - the one with the smallest mean LOOCV Rf_error
    int opt_ma_k_idx; // optimal model averaging k value index

    // Errors
    std::vector<double> k_mean_sm_errors;   // mean LOOCV squared errors for each k for single model predictions
    std::vector<double> k_mean_ma_errors;   // mean LOOCV squared errors for each k for model averaged predictions
    std::vector<double> k_mean_sm_true_errors; // mean absolute Rf_error between sm_predictions and y_true
    std::vector<double> k_mean_ma_true_errors; // mean absolute Rf_error between ma_predictions and y_true

    // The best (over all k) model evaluation
    std::vector<double> sm_predictions; // optimal k single model (before model averaging) predictions
    std::vector<double> ma_predictions; // optimal k model averaged predictions
    //std::vector<double> sm_errors;      // optimal k single model (before model averaging) Leave-One-Out Cross-Validation squared Rf_error
    //std::vector<double> ma_errors;      // optimal k model averaged approximate LOOCV squared errors

    std::vector<std::vector<double>> k_sm_predictions; // for each k single model predictions
    std::vector<std::vector<double>> k_ma_predictions; // for each k model averaged predictions
};


enum model_averaging_strategy_t {
    KERNEL_WEIGHTS_ONLY,
    ERROR_WEIGHTS_ONLY,
    KERNEL_AND_ERROR_WEIGHTS,
    KERNEL_WEIGHTS_WITH_FILTERING
};

enum error_filtering_strategy_t {
    GLOBAL_PERCENTILE,
    LOCAL_PERCENTILE,
    COMBINED_PERCENTILE,
    BEST_K_MODELS
};

#if 0
// outline of mabilo
mabilo_t mabilo(...) {
    // Initialize results vectors
    std::vector<std::vector<double>> k_predictions(n_k_values, std::vector<double>(n_points));
    results.k_mean_sm_errors.resize(n_k_values);
    results.k_mean_ma_errors.resize(n_k_values);
    // ... other initializations ...

    // Combined Phase 1 & 3: Single loop over k values
    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++) {
        // Create temporary storage for current k
        std::vector<std::vector<ulm_plus_t>> pt_models(n_points);
        for (int i = 0; i < n_points; i++) {
            pt_models[i].reserve(2 * k + 1);
        }

        // Phase 1: Compute models and single-model predictions
        for (int i = 0; i < n_points; i++) {
            // ... fit models and store in pt_models ...
            // ... compute single-model predictions and errors ...
        }

        // Phase 3: Model averaging for current k
        for (int i = 0; i < n_points; i++) {
            // ... perform model averaging using pt_models ...
            // ... store results in k_predictions ...
        }

        // Compute mean errors for current k
        results.k_mean_sm_errors[k_index] = ...;
        results.k_mean_ma_errors[k_index] = ...;

        // pt_models goes out of scope here, freeing memory
    }

    // Phase 2: Find optimal k for single models
    // ... find opt_sm_k using k_mean_sm_errors ...

    // Phase 4 (now 3): Find optimal k for model averaged predictions
    // ... find opt_ma_k using k_mean_ma_errors ...

    return results;
}
#endif

/**
 * @brief Implements the Model-Averaged LOWESS (MABILO) algorithm for robust local regression
 *
 * @details MABILO extends the traditional LOWESS algorithm by incorporating model averaging and
 * bi-k nearest neighbor structure. The algorithm consists of three phases:
 *
 * 1. For each k in [k_min, k_max]:
 *    - Fit local linear models using k-hop neighborhoods
 *    - Compute single-model predictions and their LOOCV errors
 *    - Perform model averaging using filtered models and chosen weighting strategy
 *    - Calculate mean LOOCV and true errors for both single-model and averaged predictions
 *
 * 2. Find optimal k for single-model predictions based on mean LOOCV errors
 *
 * 3. Find optimal k for model-averaged predictions based on surrogate LOOCV errors
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
 * @param model_averaging_strategy Strategy for combining model predictions:
 *        - KERNEL_WEIGHTS_ONLY: Uses only distance-based kernel weights
 *        - ERROR_WEIGHTS_ONLY: Uses only Rf_error-based weights
 *        - KERNEL_AND_ERROR_WEIGHTS: Combines both weight types
 * @param error_filtering_strategy Method for filtering models based on their errors:
 *        - GLOBAL_PERCENTILE: Uses global Rf_error distribution
 *        - LOCAL_PERCENTILE: Uses Rf_error distribution at each point
 *        - COMBINED_PERCENTILE: Applies both global and local filtering
 *        - BEST_K_MODELS: Selects k best models by Rf_error
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param model_kernel Kernel function for Rf_error-based weights (same options as distance_kernel)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 *        Makes weights at window boundaries close to but not equal to 0
 * @param epsilon Numerical stability parameter (default: 1e-15)
 * @param verbose Enable progress messages
 *
 * @return mabilo_t structure containing:
 *         - opt_sm_k: Optimal k value for single-model predictions
 *         - opt_ma_k: Optimal k value for model-averaged predictions
 *         - sm_predictions: Single-model predictions using opt_sm_k
 *         - ma_predictions: Model-averaged predictions using opt_ma_k
 *         - k_mean_sm_errors: Mean LOOCV errors for single-model predictions for each k
 *         - k_mean_ma_errors: Mean surrogate LOOCV errors for model-averaged predictions for each k
 *         - k_mean_sm_true_errors: Mean true errors for single-model predictions if y_true provided
 *         - k_mean_ma_true_errors: Mean true errors for model-averaged predictions if y_true provided
 *         - k_sm_predictions: Single-model predictions for all k values
 *         - k_ma_predictions: Model-averaged predictions for all k values
 *
 * @throws Rf_error for invalid input parameters or numerical instability
 *
 * @note
 * - Input x values must be ordered
 * - Window size for each k is 2k + 1 (k points on each side plus the center point)
 * - For points near boundaries, the window is adjusted to include available points
 * - Surrogate LOOCV errors for model-averaged predictions are derived from the
 *   minimal LOOCV errors of contributing models
 *
 *
 * @see ulm_plus_t for the structure of individual linear models
 * @see initialize_kernel() for kernel function details
 */
mabilo_plus_t mabilo_plus(const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& y_true,
                    const std::vector<double>& w,
                    int k_min,
                    int k_max,
                    model_averaging_strategy_t model_averaging_strategy,
                    error_filtering_strategy_t error_filtering_strategy,
                    int distance_kernel,
                    int model_kernel,
                    double dist_normalization_factor,
                    double epsilon,
                    bool verbose) {

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILO_PLUS");

    if (verbose) {
        Rprintf("Starting MABILO_PLUS computation\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    mabilo_plus_t results;

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
    std::vector<std::vector<double>> k_ma_predictions(n_k_values, std::vector<double>(n_points));

    // Vectors for errors during single k iteration
    std::vector<double> k_errors(n_points);
    std::vector<double> k_true_errors(n_points);

    // Vectors in results struct to store mean errors for each k
    results.k_mean_sm_errors.resize(n_k_values);
    results.k_mean_ma_errors.resize(n_k_values);
    results.k_mean_sm_true_errors.resize(n_k_values);
    results.k_mean_ma_true_errors.resize(n_k_values);

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

        // Compute mean errors for single models at current k
        results.k_mean_sm_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
        if (y_true_exists) {
            results.k_mean_sm_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
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

        initialize_kernel(model_kernel, 1.0);

        // First compute global Rf_error statistics if needed
        all_errors.reserve(n_points * (2 * k + 1));  // Approximate max size
        all_errors.clear();

        if (error_filtering_strategy == GLOBAL_PERCENTILE ||
            error_filtering_strategy == COMBINED_PERCENTILE) {
            for (int i = 0; i < n_points; i++) {
                for (const auto& model : pt_models[i]) {
                    int local_index = i - model.x_min_index;
                    all_errors.push_back(model.errors[local_index]);
                }
            }
            std::sort(all_errors.begin(), all_errors.end());
        }

        switch (model_averaging_strategy) {
        case KERNEL_WEIGHTS_ONLY: {
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
                k_ma_predictions[k_index][i] = weighted_sum / weight_sum;

                if (y_true_exists) {
                    k_true_errors[i] = std::abs(y_true[i] - k_ma_predictions[k_index][i]);
                }


                #if 0
                // Using the model with the minimal LOOCV Rf_error at the given
                // point as an estimate of the LOOCV squared Rf_error for the model
                // averaged predictions
                k_errors[i] = *std::min_element(local_errors.begin(), local_errors.end());
                #endif
            }

            // Compute mean errors for model-averaged predictions at current k
            results.k_mean_ma_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
            if (y_true_exists) {
                results.k_mean_ma_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
            }

            break;
        }

        case KERNEL_WEIGHTS_WITH_FILTERING: {

            std::vector<double> local_ws;
            std::vector<double> normalize_ws;

            for (int i = 0; i < n_points; i++) {

                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                int model_counter = 0;
                local_ws.resize(pt_models[i].size());
                normalize_ws.resize(pt_models[i].size());
                local_errors.resize(pt_models[i].size());
                double min_w = std::numeric_limits<double>::max();
                double max_w = 0.0;

                for (const auto& model : pt_models[i]) {
                    int local_index = i - model.x_min_index;
                    double weight = model.w[local_index];
                    weighted_sum += weight * model.predictions[local_index];
                    weight_sum += weight;
                    min_w = std::min(min_w, weight);
                    max_w = std::max(max_w, weight);
                    local_ws[model_counter] = weight;
                    local_errors[model_counter++] = model.errors[local_index];
                }

                // Store model-averaged prediction and its errors
                k_ma_predictions[k_index][i] = weighted_sum / weight_sum;

                if (y_true_exists) {
                    k_true_errors[i] = std::abs(y_true[i] - k_ma_predictions[k_index][i]);
                }

                double w_range = max_w - min_w;
                for (auto& w: local_ws) {
                    w = 1.0 - (w - min_w) / w_range;
                }

                std::vector<double> w_weights(local_ws.size());
                kernel_fn(local_ws.data(), local_ws.size(), w_weights.data());

                weighted_sum = 0.0;
                weight_sum = 0.0;
                for (size_t j = 0; j < w_weights.size(); j++) {
                    weighted_sum += local_errors[j] * w_weights[j];
                    weight_sum += w_weights[j];
                }

                k_errors[i] = weighted_sum / weight_sum;
            }

            // Compute mean errors for model-averaged predictions at current k
            results.k_mean_ma_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
            if (y_true_exists) {
                results.k_mean_ma_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
            }

            break;
        }

        case ERROR_WEIGHTS_ONLY: {
            for (int i = 0; i < n_points; i++) {
                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                filtered_models.clear();

                // Filter models based on chosen strategy
                switch (error_filtering_strategy) {
                case GLOBAL_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.25 * all_errors.size())];
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case LOCAL_PERCENTILE: {
                    local_errors.clear();
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        local_errors.push_back(model.errors[local_index]);
                    }
                    std::sort(local_errors.begin(), local_errors.end());
                    double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= local_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case COMBINED_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.5 * all_errors.size())];
                    local_errors.clear();

                    // First apply global filter
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            local_errors.push_back(model.errors[local_index]);
                        }
                    }

                    if (!local_errors.empty()) {
                        std::sort(local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case BEST_K_MODELS: {
                    const int k = 3;  // Number of best models to use
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        filtered_models.push_back({model.errors[local_index], &model});
                    }

                    if (filtered_models.size() > k) {
                        std::partial_sort(filtered_models.begin(),
                                          filtered_models.begin() + k,
                                          filtered_models.end(),
                                          [](const auto& a, const auto& b) {
                                              return a.first < b.first;
                                          });
                        filtered_models.resize(k);
                    }
                    break;
                }
                }

                // If we have filtered models, compute their min-max normalized weights
                if (!filtered_models.empty()) {
                    // Find min and max errors among filtered models
                    double min_error = filtered_models[0].first;
                    double max_error = filtered_models[0].first;
                    for (const auto& [Rf_error, _] : filtered_models) {
                        min_error = std::min(min_error, Rf_error);
                        max_error = std::max(max_error, Rf_error);
                    }

                    double error_range = max_error - min_error;
                    if (error_range > epsilon) {
                        for (const auto& [Rf_error, model] : filtered_models) {
                            int local_index = i - model->x_min_index;
                            double normalized_error = (Rf_error - min_error) / error_range;

                            // Use model_kernel for Rf_error weighting
                            std::vector<double> error_weight(1);
                            kernel_fn(&normalized_error, 1, error_weight.data());

                            weighted_sum += error_weight[0] * model->predictions[local_index];
                            weight_sum += error_weight[0];
                        }
                    } else {
                        // All errors are essentially the same, use simple average
                        for (const auto& [_, model] : filtered_models) {
                            int local_index = i - model->x_min_index;
                            weighted_sum += model->predictions[local_index];
                            weight_sum += 1.0;
                        }
                    }

                    k_ma_predictions[k_index][i] = weighted_sum / weight_sum;
                    k_errors[i] = min_error;
                } else {
                    // If no models passed filtering, use the best available model
                    auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                       [i](const auto& a, const auto& b) {
                                                           return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                       });
                    k_ma_predictions[k_index][i] = best_model->predictions[i - best_model->x_min_index];
                    k_errors[i] = best_model->errors[i - best_model->x_min_index];
                    if (y_true_exists) {
                        k_true_errors[i] = std::abs(y_true[i] - k_ma_predictions[k_index][i]);
                    }
                }
            }

            results.k_mean_ma_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
            if (y_true_exists) {
                results.k_mean_ma_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
            }

            break;
        }

        case KERNEL_AND_ERROR_WEIGHTS: {
            for (int i = 0; i < n_points; i++) {
                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                filtered_models.clear();

                // Filter models based on chosen strategy
                switch (error_filtering_strategy) {
                case GLOBAL_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.25 * all_errors.size())];
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case LOCAL_PERCENTILE: {
                    local_errors.clear();
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        local_errors.push_back(model.errors[local_index]);
                    }
                    std::sort(local_errors.begin(), local_errors.end());
                    double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= local_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case COMBINED_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.5 * all_errors.size())];
                    local_errors.clear();

                    // First apply global filter
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            local_errors.push_back(model.errors[local_index]);
                        }
                    }

                    if (!local_errors.empty()) {
                        std::sort(local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case BEST_K_MODELS: {
                    const int k = 3;  // Number of best models to use
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        filtered_models.push_back({model.errors[local_index], &model});
                    }

                    if (filtered_models.size() > k) {
                        std::partial_sort(filtered_models.begin(),
                                          filtered_models.begin() + k,
                                          filtered_models.end(),
                                          [](const auto& a, const auto& b) {
                                              return a.first < b.first;
                                          });
                        filtered_models.resize(k);
                    }
                    break;
                }
                }

                // If we have filtered models, compute combined weights
                if (!filtered_models.empty()) {
                    // Find min and max errors among filtered models
                    double min_error = filtered_models[0].first;
                    double max_error = filtered_models[0].first;
                    for (const auto& [Rf_error, _] : filtered_models) {
                        min_error = std::min(min_error, Rf_error);
                        max_error = std::max(max_error, Rf_error);
                    }

                    double error_range = max_error - min_error;
                    if (error_range > epsilon) {
                        for (const auto& [Rf_error, model] : filtered_models) {
                            int local_index = i - model->x_min_index;

                            // Get distance-based kernel weight (already normalized)
                            double kernel_weight = model->w[local_index];

                            // Compute Rf_error-based weight using model_kernel
                            double normalized_error = (Rf_error - min_error) / error_range;
                            std::vector<double> error_weight(1);
                            kernel_fn(&normalized_error, 1, error_weight.data());

                            // Combine weights by multiplication
                            double combined_weight = kernel_weight * error_weight[0];
                            weighted_sum += combined_weight * model->predictions[local_index];
                            weight_sum += combined_weight;
                        }
                    } else {
                        // All errors are essentially the same, use only kernel weights
                        for (const auto& [_, model] : filtered_models) {
                            int local_index = i - model->x_min_index;
                            double kernel_weight = model->w[local_index];
                            weighted_sum += kernel_weight * model->predictions[local_index];
                            weight_sum += kernel_weight;
                        }
                    }

                    k_ma_predictions[k_index][i] = weighted_sum / weight_sum;
                    k_errors[i] = min_error;
                } else {
                    // If no models passed filtering, use the best available model
                    auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                       [i](const auto& a, const auto& b) {
                                                           return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                       });
                    k_ma_predictions[k_index][i] = best_model->predictions[i - best_model->x_min_index];
                    k_errors[i] = best_model->errors[i - best_model->x_min_index];
                    if (y_true_exists) {
                        k_true_errors[i] = std::abs(y_true[i] - k_ma_predictions[k_index][i]);
                    }
                }
            }

            results.k_mean_ma_errors[k_index] = std::accumulate(k_errors.begin(), k_errors.end(), 0.0) / n_points;
            if (y_true_exists) {
                results.k_mean_ma_true_errors[k_index] = std::accumulate(k_true_errors.begin(), k_true_errors.end(), 0.0) / n_points;
            }

            break;
        }
        default:
            REPORT_ERROR("Unknown model averaging strategy");
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
    // Phase 3: Find optimal k for single-model predictions
    //
    // -------------------------------------------------------------------------------------------
    auto opt_sm_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 3: Finding optimal k ... ");
    }

    auto min_it = std::min_element(results.k_mean_sm_errors.begin(), results.k_mean_sm_errors.end());
    results.opt_sm_k_idx = std::distance(results.k_mean_sm_errors.begin(), min_it);
    results.opt_sm_k = k_min + results.opt_sm_k_idx;

    // Save optimal single-model predictions before they're overwritten
    results.sm_predictions = k_sm_predictions[results.opt_sm_k_idx];
    results.k_sm_predictions = k_sm_predictions;

    if (verbose) {
        elapsed_time(opt_sm_k_ptm, "Done");
        mem_tracker.report();
    }

    // -------------------------------------------------------------------------------------------
    //
    // Phase 4: Find optimal k for model-averaged predictions
    // using approximated LOOCV squared errros for these predictions
    //
    // -------------------------------------------------------------------------------------------
    auto opt_ma_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 4: Finding optimal model averaged predictions over all k's ... ");
    }
    min_it = std::min_element(results.k_mean_ma_errors.begin(), results.k_mean_ma_errors.end());
    results.opt_ma_k_idx = std::distance(results.k_mean_ma_errors.begin(), min_it);
    results.opt_ma_k = k_min + results.opt_ma_k_idx;
    results.ma_predictions = k_ma_predictions[results.opt_ma_k_idx];
    results.k_ma_predictions = std::move(k_ma_predictions);
    if (verbose) {
        elapsed_time(opt_ma_k_ptm, "Done");
        mem_tracker.report();
    }


    if (verbose) {
        elapsed_time(total_ptm, "\nTotal MABILO_PLUS computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}

/**
 * @brief R interface for MABILO (Model-Averaged Locally Weighted Scatterplot Smoothing)
 *
 * @param s_x Vector of x coordinates
 * @param s_y Vector of y coordinates (response values)
 * @param s_y_true Optional vector of true y values for Rf_error calculation
 * @param s_w Vector of weights for each point
 * @param s_k_min Minimum number of neighbors (must be positive)
 * @param s_k_max Maximum number of neighbors (must be greater than k_min and not exceed input size)
 * @param s_model_averaging_strategy Strategy for averaging models:
 *        - "kernel_weights_only": Uses only distance-based kernel weights
 *        - "error_weights_only": Uses only Rf_error-based weights
 *        - "kernel_and_error_weights": Combines both weight types
 * @param s_error_filtering_strategy Strategy for filtering models:
 *        - "global_percentile": Uses global Rf_error distribution
 *        - "local_percentile": Uses Rf_error distribution at each point
 *        - "combined_percentile": Applies both global and local filtering
 *        - "best_k_models": Selects k best models by Rf_error
 * @param s_distance_kernel Kernel type for distance-based weights (0: Tricube, 1: Epanechnikov, 2: Exponential)
 * @param s_model_kernel Kernel type for Rf_error-based weights (same options as distance_kernel)
 * @param s_dist_normalization_factor Factor for normalizing distances
 * @param s_epsilon Small number to prevent division by zero
 * @param s_verbose Whether to print progress information
 *
 * @return A list containing:
 * - k_values: Vector of k values tested
 * - opt_sm_k: Optimal k value for single-model predictions
 * - opt_sm_k_idx: Index of optimal k value for single-model predictions
 * - opt_ma_k: Optimal k value for model-averaged predictions
 * - opt_ma_k_idx: Index of optimal k value for model-averaged predictions
 * - k_mean_sm_errors: Mean LOOCV errors for single-model predictions for each k
 * - k_mean_ma_errors: Mean surrogate LOOCV errors for model-averaged predictions for each k
 * - k_mean_sm_true_errors: Mean true errors for single-model predictions if y_true provided
 * - k_mean_ma_true_errors: Mean true errors for model-averaged predictions if y_true provided
 * - sm_predictions: Single-model predictions using optimal k
 * - ma_predictions: Model-averaged predictions using optimal k
 * - k_sm_predictions: Single-model predictions for all k values
 * - k_ma_predictions: Model-averaged predictions for all k values
 *
 * @throws Rf_error if input vectors have inconsistent lengths
 * @throws Rf_error if k_min or k_max values are invalid
 * @throws Rf_error if unknown model averaging or Rf_error filtering strategy is provided
 */
SEXP S_mabilo_plus(SEXP s_x,
                SEXP s_y,
                SEXP s_y_true,
                SEXP s_w,
                SEXP s_k_min,
                SEXP s_k_max,
                SEXP s_model_averaging_strategy,
                SEXP s_error_filtering_strategy,
                SEXP s_distance_kernel,
                SEXP s_model_kernel,
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

    model_averaging_strategy_t model_averaging_strategy;
    error_filtering_strategy_t error_filtering_strategy;

    const char* strat = CHAR(STRING_ELT(s_model_averaging_strategy, 0));
    if (strcmp(strat, "kernel_weights_only") == 0) {
        model_averaging_strategy = KERNEL_WEIGHTS_ONLY;
    } else if (strcmp(strat, "error_weights_only") == 0) {
        model_averaging_strategy = ERROR_WEIGHTS_ONLY;
    } else if (strcmp(strat, "kernel_and_error_weights") == 0) {
        model_averaging_strategy = KERNEL_AND_ERROR_WEIGHTS;
    } else if (strcmp(strat, "kernel_weights_with_filtering") == 0) {
        model_averaging_strategy = KERNEL_WEIGHTS_WITH_FILTERING;
    } else {
        Rf_error("Unknown model averaging strategy");
    }

    const char* filter = CHAR(STRING_ELT(s_error_filtering_strategy, 0));
    if (strcmp(filter, "global_percentile") == 0) {
        error_filtering_strategy = GLOBAL_PERCENTILE;
    } else if (strcmp(filter, "local_percentile") == 0) {
        error_filtering_strategy = LOCAL_PERCENTILE;
    } else if (strcmp(filter, "combined_percentile") == 0) {
        error_filtering_strategy = COMBINED_PERCENTILE;
    } else if (strcmp(filter, "best_k_models") == 0) {
        error_filtering_strategy = BEST_K_MODELS;
    } else {
        Rf_error("Unknown Rf_error filtering strategy");
    }

    int distance_kernel = INTEGER(s_distance_kernel)[0];
    int model_kernel = INTEGER(s_model_kernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    mabilo_plus_t mabilo_plus_results = mabilo_plus(x,
                                                    y,
                                                    y_true,
                                                    w,
                                                    k_min,
                                                    k_max,
                                                    model_averaging_strategy,
                                                    error_filtering_strategy,
                                                    distance_kernel,
                                                    model_kernel,
                                                    dist_normalization_factor,
                                                    epsilon,
                                                    verbose);

    // Creating return list
    const int N_COMPONENTS = 13;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    std::vector<int> k_values(mabilo_plus_results.k_mean_sm_errors.size());
    for (int k_index = 0, k = k_min; k <= k_max; k++, k_index++)
        k_values[k_index] = k;
    SET_VECTOR_ELT(result, 0, convert_vector_int_to_R(k_values)); n_protected++;

    // opt_sm_k
    SEXP s_opt_sm_k = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_sm_k)[0] = mabilo_plus_results.opt_sm_k; // 1-base
    SET_VECTOR_ELT(result, 1, s_opt_sm_k);

    SEXP s_opt_sm_k_idx = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_sm_k_idx)[0] = mabilo_plus_results.opt_sm_k_idx + 1; // 1-base
    SET_VECTOR_ELT(result, 2, s_opt_sm_k_idx);

    // opt_ma_k
    SEXP s_opt_ma_k = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_ma_k)[0] = mabilo_plus_results.opt_ma_k; // 1-base
    SET_VECTOR_ELT(result, 3, s_opt_ma_k);

    SEXP s_opt_ma_k_idx = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
    INTEGER(s_opt_ma_k_idx)[0] = mabilo_plus_results.opt_ma_k_idx + 1; // 1-base
    SET_VECTOR_ELT(result, 4, s_opt_ma_k_idx);

    // k_mean_sm_errors
    SET_VECTOR_ELT(result, 5, PROTECT(convert_vector_double_to_R(mabilo_plus_results.k_mean_sm_errors))); n_protected++;
    // k_mean_ma_errors
    SET_VECTOR_ELT(result, 6, PROTECT(convert_vector_double_to_R(mabilo_plus_results.k_mean_ma_errors))); n_protected++;

    // true errors
    if (y_true_exists) {
        SET_VECTOR_ELT(result, 7, PROTECT(convert_vector_double_to_R(mabilo_plus_results.k_mean_sm_true_errors))); n_protected++;
        SET_VECTOR_ELT(result, 8, PROTECT(convert_vector_double_to_R(mabilo_plus_results.k_mean_ma_true_errors))); n_protected++;
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
    }

    // predictions
    SET_VECTOR_ELT(result, 9, PROTECT(convert_vector_double_to_R(mabilo_plus_results.sm_predictions))); n_protected++;
    SET_VECTOR_ELT(result, 10, PROTECT(convert_vector_double_to_R(mabilo_plus_results.ma_predictions))); n_protected++;

    // k_predictions
    SET_VECTOR_ELT(result, 11, convert_vector_vector_double_to_R(mabilo_plus_results.k_sm_predictions)); n_protected++;
    SET_VECTOR_ELT(result, 12, convert_vector_vector_double_to_R(mabilo_plus_results.k_ma_predictions)); n_protected++;

    // Setting names for return list
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, Rf_mkChar("k_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("opt_sm_k"));
    SET_STRING_ELT(names, 2, Rf_mkChar("opt_sm_k_idx"));
    SET_STRING_ELT(names, 3, Rf_mkChar("opt_ma_k"));
    SET_STRING_ELT(names, 4, Rf_mkChar("opt_ma_k_idx"));
    SET_STRING_ELT(names, 5, Rf_mkChar("k_mean_sm_errors"));
    SET_STRING_ELT(names, 6, Rf_mkChar("k_mean_ma_errors"));
    SET_STRING_ELT(names, 7, Rf_mkChar("k_mean_sm_true_errors"));
    SET_STRING_ELT(names, 8, Rf_mkChar("k_mean_ma_true_errors"));
    SET_STRING_ELT(names, 9, Rf_mkChar("sm_predictions"));
    SET_STRING_ELT(names, 10, Rf_mkChar("ma_predictions"));
    SET_STRING_ELT(names, 11, Rf_mkChar("k_sm_predictions"));
    SET_STRING_ELT(names, 12, Rf_mkChar("k_ma_predictions"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}
