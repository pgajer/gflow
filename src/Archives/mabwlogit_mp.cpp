#if 0

#include <R.h>
#include <Rinternals.h>

#undef length  // to resolve naming conflict between the R macro length defined in Rinternals.h and a member function in the C++ standard library's codecvt class

#include <execution>
#include <mutex>
#include <ANN/ANN.h>
#include <Eigen/Dense>
#include <vector>
#include <numeric>    // for std::iota
#include <random>     // for std::mt19937
#include <algorithm>  // for std::shuffle

#include "ulogit.hpp"
#include "kernels.h"
#include "mabwlogit.hpp"

#include "cpp_utils.h" // for debugging and elapsed.time
#include <fstream>     // for debugging and elapsed.time



/**
 * @brief Performs parallel model-averaged local logistic regression with bandwidth selection
 *
 * This function implements a parallel version of the model-averaged local logistic
 * regression algorithm using C++17's parallel algorithms. It parallelizes the
 * local model fitting process while maintaining thread safety for shared resources.
 * It supports the same three modes of operation as the sequential version:
 * 1. Fixed bandwidth mode (pilot_bandwidth > 0)
 * 2. LOOCV-based bandwidth selection (cv_folds = 0)
 * 3. K-fold cross-validation bandwidth selection (cv_folds > 0)
 *
 * For bandwidth selection, the function creates a grid of candidate bandwidths
 * around a rule-of-thumb value and selects the optimal bandwidth that minimizes
 * the cross-validation error.
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param fit_quadratic If true, fits quadratic term in addition to linear
 * @param pilot_bandwidth Initial bandwidth for local fitting (if > 0, used as fixed bandwidth)
 * @param kernel_type Integer specifying the kernel function type:
 *                    1: Epanechnikov, 2: Triangular, 4: Laplace,
 *                    5: Normal, 6: Biweight, 7: Tricube
 * @param min_points Minimal number of points in each local window
 * @param cv_folds Number of cross-validation folds (0 for LOOCV approximation)
 * @param n_bws Number of bandwidths in the grid for bandwidth selection
 * @param min_bw_factor Lower bound factor for bandwidth grid relative to rule-of-thumb
 * @param max_bw_factor Upper bound factor for bandwidth grid relative to rule-of-thumb
 * @param max_iterations Maximum number of iterations for logistic regression fitting
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for logistic regression fitting
 * @param with_errors If true, compute and return prediction errors
 * @param with_bw_preditions If true, store predictions for each bandwidth
 *
 * @return mabwlogit_t structure containing:
 *         - candidate_bandwidths: Vector of tested bandwidths
 *         - bw_predictions: Predictions for each bandwidth (if with_bw_preditions = true)
 *         - mean_errors: Mean CV error for each bandwidth
 *         - opt_bw: Optimal bandwidth
 *         - predictions: Final predictions using optimal bandwidth
 *         - errors: Prediction errors (if with_errors = true)
 *         - beta1s: First-order coefficients
 *         - beta2s: Second-order coefficients (if fit_quadratic = true)
 *         - Input parameters used in the fit
 *
 * @throw std::invalid_argument if input vectors are empty or of unequal length,
 *        or if cv_folds is negative
 * @throw std::runtime_error if no successful model fits are achieved
 *
 * @note The function requires at least 6 points in each local neighborhood for
 *       stable estimation.
 * @note This is the parallel version of the algorithm using C++17 parallel algorithms.
 * @note Thread-safe access to shared resources (e.g., k-d tree) is managed internally
 *       using mutexes.
 * @note Performance improvement over the sequential version depends on the size of
 *       the input data and the available hardware concurrency.
 */
mabwlogit_t mabwlogit_mp(
    const std::vector<double>& x,
    const std::vector<double>& y,
    bool fit_quadratic,
    double pilot_bandwidth,
    int kernel_type,
    int min_points,
    int cv_folds,
    int n_bws,
    double min_bw_factor,
    double max_bw_factor,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_errors,
    bool with_bw_preditions,
    bool verbose) {

    mabwlogit_t result;
    result.fit_quadratic   = fit_quadratic;
    result.pilot_bandwidth = pilot_bandwidth;
    result.kernel_type     = kernel_type;
    result.cv_folds        = cv_folds;
    result.min_bw_factor   = min_bw_factor;
    result.max_bw_factor   = max_bw_factor;
    result.max_iterations = max_iterations;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    // Create data points array for ANN
    int n_pts = x.size();

    // If verbose, display initial message
    if (verbose) {
        REprintf("Starting local model fitting for %d points\n", n_pts);
    }

    ANNpointArray data_points = annAllocPts(n_pts, 1);  // 1 dimension
    for(int i = 0; i < n_pts; i++) {
        data_points[i][0] = x[i];
    }

    // Build kd-tree
    ANNkd_tree* kdtree = new ANNkd_tree(data_points, n_pts, 1);

    // Lambda function for fitting local logistic regression with a specific bandwidth
    auto fit_local_logistic = [kdtree,&cv_folds,&verbose](
        const std::vector<double>& x,
        const std::vector<double>& y,
        double bandwidth,
        bool fit_quadratic,
        int kernel_type,
        int min_points,
        int max_iterations,
        double ridge_lambda,
        double tolerance,
        bool with_errors) {

        struct thread_local_resources_t {
            ANNpoint query_point;
            ANNidxArray nn_idx;
            ANNdistArray nn_dists;
            std::vector<double> local_x;
            std::vector<double> local_y;
            std::vector<double> local_w;
            std::vector<double> local_d;
            std::vector<int> local_indices;

            thread_local_resources_t(int min_points, int n_pts) {
                query_point = annAllocPt(1);
                nn_idx = new ANNidx[min_points];
                nn_dists = new ANNdist[min_points];
                local_x.reserve(n_pts);
                local_y.reserve(n_pts);
                local_w.reserve(n_pts);
                local_d.reserve(n_pts);
                local_indices.reserve(n_pts);
            }

            ~thread_local_resources_t() {
                annDeallocPt(query_point);
                delete[] nn_idx;
                delete[] nn_dists;
            }
        };

        int n_pts = x.size();
        // Initialize progress counter and mutex for progress tracking
        int progress = 0;
        std::mutex progress_mutex;
        std::vector<std::mutex> pt_pred_weights_mutexes(n_pts);  // One mutex per point
        static std::mutex eigen_mutex;

        std::vector<std::vector<point_t>> pt_pred_weights(n_pts);
        std::vector<double> predictions(n_pts);
        std::vector<double> errors;
        if (with_errors) errors.resize(n_pts);
        std::vector<double> beta1s(n_pts);
        std::vector<double> beta2s;  // either mean beta[2] coefficients (over different quadratic models) in the fit_quadratic case OR difference of consecutive beta[1] coefficients (over sorted x) assuming that x is approximately uniform grid
        if (fit_quadratic) beta2s.resize(n_pts);

        initialize_kernel(kernel_type, 1.0);

        // Pre-allocate arrays for ANN search
        ANNpoint query_point = annAllocPt(1);
        ANNidxArray nn_idx = new ANNidx[min_points];
        ANNdistArray nn_dists = new ANNdist[min_points];

        // Pre-allocate vectors
        std::vector<double> local_x, local_y, local_w, local_d;
        local_x.reserve(n_pts);
        local_y.reserve(n_pts);
        local_w.reserve(n_pts);
        local_d.reserve(n_pts);

        // Store indices along with local data for accurate mapping
        std::vector<int> local_indices;
        local_indices.reserve(n_pts);

        // Mutex for thread-safe kdtree access
        std::mutex kdtree_mutex;
        std::for_each(
            std::execution::par_unseq,
            std::begin(std::vector<int>(n_pts)),
            std::end(std::vector<int>(n_pts)),
            [&progress, &progress_mutex, &pt_pred_weights_mutexes, &bandwidth, &min_points, &x, &y, &pt_pred_weights, &fit_quadratic, &with_errors, &kdtree_mutex, &max_iterations, &ridge_lambda, &tolerance, kdtree, verbose, n_pts](int pti) {

                // Thread-local storage
                thread_local thread_local_resources_t resources(min_points, n_pts);

                double center = x[pti];

                // Clear vectors for reuse
                resources.local_x.clear();
                resources.local_y.clear();
                resources.local_w.clear();
                resources.local_d.clear();
                resources.local_indices.clear();

                int p_pts = n_pts / 20;

                // Update progress periodically (e.g., every 5%)
                if (verbose) {
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    progress++;
                    if (progress % p_pts == 0 || progress == n_pts) {  // Report every 5%
                        double percent = (100.0 * progress) / n_pts;
                        REprintf("\rProgress: %.1f%% (%d/%d points processed)",
                                 percent, progress, n_pts);
                        if (progress == n_pts) {
                            REprintf("\n");
                        }
                    }
                }

                // First try bandwidth-based neighborhood
                for (int i = 0; i < n_pts; ++i) {
                    double dist = std::abs(x[i] - center) / bandwidth;
                    if (dist < 1.0) {
                        resources.local_indices.push_back(i);
                    }
                }

                int n_window_pts = resources.local_indices.size();
                if (n_window_pts < min_points) {
                    // Protect kdtree access with mutex
                    {
                        std::lock_guard<std::mutex> lock(kdtree_mutex);
                        resources.query_point[0] = center;
                        kdtree->annkSearch(
                            resources.query_point,
                            min_points,
                            resources.nn_idx,
                            resources.nn_dists
                            );
                    }

                    // Process k-NN results
                    for (int i = 0; i < min_points; ++i) {
                        double shifted_x = x[resources.nn_idx[i]] - center;
                        resources.local_x.push_back(shifted_x);
                        resources.local_y.push_back(y[resources.nn_idx[i]]);
                        resources.local_d.push_back(std::abs(shifted_x));
                    }

                    double max_dist = *std::max_element(
                        resources.local_d.begin(),
                        resources.local_d.end()
                        );
                    if (max_dist < std::numeric_limits<double>::epsilon()) {
                        max_dist = std::numeric_limits<double>::epsilon();
                    }

                    for (auto& d : resources.local_d) {
                        d /= max_dist;
                    }

                    n_window_pts = min_points;
                } else {
                    for (int i = 0; i < n_window_pts; ++i) {
                        double shifted_x = x[resources.local_indices[i]] - center;
                        resources.local_x.push_back(shifted_x);
                        resources.local_y.push_back(y[resources.local_indices[i]]);
                        resources.local_d.push_back(std::abs(shifted_x) / bandwidth);
                    }
                }

                // Compute weights
                resources.local_w.resize(n_window_pts);
                kernel_fn(resources.local_d.data(), n_window_pts, resources.local_w.data());

                // Process weights
                double sum_weights = std::accumulate(
                    resources.local_w.begin(),
                    resources.local_w.end(),
                    0.0
                    );
                if (sum_weights < std::numeric_limits<double>::epsilon()) {
                    sum_weights = std::numeric_limits<double>::epsilon();
                }
                for (auto& w : resources.local_w) {
                    w /= sum_weights;
                }

                eigen_ulogit_t fit_result;
                { // Fit local model
                    std::lock_guard<std::mutex> lock(eigen_mutex);
                    fit_result = eigen_ulogit_fit(
                        resources.local_x.data(),
                        resources.local_y.data(),
                        resources.local_w,
                        fit_quadratic,
                        max_iterations,
                        ridge_lambda,
                        tolerance,
                        with_errors
                        );
                }


                // Store results for this point
                for (size_t i = 0; i < resources.local_x.size(); ++i) {
                    int orig_idx = resources.local_indices[i];
                    point_t pt;
                    pt.w = resources.local_w[i];
                    pt.p = fit_result.predictions[i];
                    pt.beta1 = fit_result.beta[1];
                    if (fit_quadratic) pt.beta2 = fit_result.beta[2];
                    if (with_errors) pt.e = fit_result.errors[i];

                    // Protect the push_back operation
                    {
                        std::lock_guard<std::mutex> lock(pt_pred_weights_mutexes[orig_idx]);
                        pt_pred_weights[orig_idx].push_back(pt);
                    }
                }
            });

        // Compute weighted averages
        for (int pti = 0; pti < n_pts; pti++) {
            const auto& v = pt_pred_weights[pti];
            if (v.empty()) {
                predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double total_weight = 0.0;
            double weighted_prediction = 0.0;
            double weighted_error = 0.0;
            double weighted_beta1 = 0.0;
            double weighted_beta2 = 0.0;
            for (const auto& pt : v) {
                weighted_prediction += pt.w * pt.p;
                weighted_beta1 += pt.w * pt.beta1;
                if (with_errors) weighted_error += pt.w * pt.e;
                if (fit_quadratic) weighted_beta2 += pt.w * pt.beta2;
                total_weight += pt.w;
            }

            predictions[pti] = total_weight > 0.0 ? weighted_prediction / total_weight :
                std::numeric_limits<double>::quiet_NaN();

            beta1s[pti] = total_weight > 0.0 ? weighted_beta1 / total_weight :
                std::numeric_limits<double>::quiet_NaN();


            if (with_errors) errors[pti] = total_weight > 0.0 ? weighted_error / total_weight :
                                 std::numeric_limits<double>::quiet_NaN();

            if (fit_quadratic) beta2s[pti] = total_weight > 0.0 ? weighted_beta2 / total_weight :
                                   std::numeric_limits<double>::quiet_NaN();

        }

        // Clean up ANN allocations
        annDeallocPt(query_point);
        delete[] nn_idx;
        delete[] nn_dists;

        local_logit_t ll;
        ll.preds = std::move(predictions);
        if (with_errors) ll.errs = std::move(errors);
        ll.beta1s = std::move(beta1s);
        if (fit_quadratic) ll.beta2s = std::move(beta2s);

        return ll;
    };

    // If pilot_bandwidth > 0, use it directly
    if (pilot_bandwidth > 0) {
        auto ll = fit_local_logistic(x,
                                     y,
                                     pilot_bandwidth,
                                     fit_quadratic,
                                     kernel_type,
                                     min_points,
                                     max_iterations,
                                     ridge_lambda,
                                     tolerance,
                                     with_errors);

        result.predictions = std::move(ll.preds);
        result.beta1s = std::move(ll.beta1s);

        if (with_errors) result.errors = std::move(ll.errs);

        if (fit_quadratic) {
            result.beta2s = std::move(ll.beta2s);
        } else {
            // First, create sorted indices of x values
            std::vector<size_t> sort_idx(x.size());
            std::iota(sort_idx.begin(), sort_idx.end(), 0);
            std::sort(sort_idx.begin(), sort_idx.end(),
                      [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });

            // Create sorted arrays of x and beta1s
            std::vector<double> x_sorted(x.size());
            std::vector<double> beta1s_sorted(x.size());
            for (size_t i = 0; i < x.size(); i++) {
                x_sorted[i] = x[sort_idx[i]];
                beta1s_sorted[i] = result.beta1s[sort_idx[i]];
            }

            // Initialize beta2s with the same size as x
            result.beta2s.resize(x.size(), std::numeric_limits<double>::quiet_NaN());

            // Compute numerical derivatives using central differences where possible
            // and forward/backward differences at the endpoints
            for (size_t i = 0; i < x_sorted.size(); i++) {
                if (i == 0) {
                    // Forward difference for first point
                    if (i + 1 < x_sorted.size()) {
                        result.beta2s[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i]) /
                            (x_sorted[i + 1] - x_sorted[i]);
                    }
                } else if (i == x_sorted.size() - 1) {
                    // Backward difference for last point
                    result.beta2s[sort_idx[i]] = (beta1s_sorted[i] - beta1s_sorted[i - 1]) /
                        (x_sorted[i] - x_sorted[i - 1]);
                } else {
                    // Central difference for interior points
                    result.beta2s[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i - 1]) /
                        (x_sorted[i + 1] - x_sorted[i - 1]);
                }
            }
        }

        result.opt_bw = std::numeric_limits<double>::quiet_NaN();

        return result;
    }

    // Create bandwidth grid spanning min_bw_factor to max_bw_factor times data width
    result.candidate_bandwidths.resize(n_bws);
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    double x_range = x_max - x_min;

    if (x_range < std::numeric_limits<double>::epsilon()) {
        error("Input x values are effectively constant");
    }

    double min_bw = min_bw_factor * x_range;
    double max_bw = max_bw_factor * x_range;
    double dx = (max_bw - min_bw) / (n_bws - 1);
    for(int i = 0; i < n_bws; i++) {
        result.candidate_bandwidths[i] = min_bw + i * dx;
    }

    if (with_bw_preditions) {
        result.bw_predictions.resize(n_bws);
        result.bw_errors.resize(n_bws);
    }
    result.mean_errors.resize(n_bws);

    // If cv_folds = 0, use LOOCV approximation
    if (cv_folds == 0) {

        result.bw_isq_beta2.resize(n_bws);

        for (int i = 0; i < n_bws; i++) {
            auto ll = fit_local_logistic(x,
                                         y,
                                         result.candidate_bandwidths[i],
                                         fit_quadratic,
                                         kernel_type,
                                         min_points,
                                         max_iterations,
                                         ridge_lambda,
                                         tolerance,
                                         true);
            // Compute mean error
            double total_error = 0.0;
            int valid_count = 0;
            for (const auto& err : ll.errs) {
                if (!std::isnan(err)) {
                    total_error += err;
                    valid_count++;
                }
            }
            result.mean_errors[i] = valid_count > 0 ? total_error / valid_count : std::numeric_limits<double>::infinity();

            if (fit_quadratic) {
                double total_sq_beta2 = 0.0;
                int valid_count = 0;
                for (const auto& b2 : ll.beta2s) {
                    if (!std::isnan(b2)) {
                        total_sq_beta2 += b2 * b2;
                        valid_count++;
                    }
                }
                result.bw_isq_beta2[i] = valid_count > 0 ? total_sq_beta2 / valid_count : std::numeric_limits<double>::infinity();
            } else {
                // First compute numerical derivatives as above
                std::vector<double> beta2s_current(x.size());

                // Create sorted arrays
                std::vector<size_t> sort_idx(x.size());
                std::iota(sort_idx.begin(), sort_idx.end(), 0);
                std::sort(sort_idx.begin(), sort_idx.end(),
                          [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });

                std::vector<double> x_sorted(x.size());
                std::vector<double> beta1s_sorted(x.size());
                for (size_t i = 0; i < x.size(); i++) {
                    x_sorted[i] = x[sort_idx[i]];
                    beta1s_sorted[i] = ll.beta1s[sort_idx[i]];
                }

                // Compute derivatives
                for (size_t i = 0; i < x_sorted.size(); i++) {
                    if (i == 0) {
                        beta2s_current[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i]) /
                            (x_sorted[i + 1] - x_sorted[i]);
                    } else if (i == x_sorted.size() - 1) {
                        beta2s_current[sort_idx[i]] = (beta1s_sorted[i] - beta1s_sorted[i - 1]) /
                            (x_sorted[i] - x_sorted[i - 1]);
                    } else {
                        beta2s_current[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i - 1]) /
                            (x_sorted[i + 1] - x_sorted[i - 1]);
                    }
                }

                // Compute integral of squared beta2s
                double total_sq_beta2 = 0.0;
                int valid_count = 0;
                for (const auto& b2 : beta2s_current) {
                    if (!std::isnan(b2)) {
                        total_sq_beta2 += b2 * b2;
                        valid_count++;
                    }
                }
                result.bw_isq_beta2[i] = valid_count > 0 ? total_sq_beta2 / valid_count : std::numeric_limits<double>::infinity();
            }

            if (with_bw_preditions) {
                result.bw_predictions[i] = std::move(ll.preds);
                result.bw_errors[i] = std::move(ll.errs);
            }
        }
    } else {
        // K-fold cross-validation
        int n = x.size();
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(indices.begin(), indices.end(), gen);

        int fold_size = n / cv_folds;
        for (int i = 0; i < n_bws; i++) {
            double total_error = 0.0;
            int valid_count = 0;

            for (int fold = 0; fold < cv_folds; fold++) {
                // Create training and validation sets
                std::vector<double> train_x, train_y, valid_x, valid_y;
                int start_idx = fold * fold_size;
                int end_idx = (fold == cv_folds - 1) ? n : (fold + 1) * fold_size;

                for (int j = 0; j < n; j++) {
                    if (j >= start_idx && j < end_idx) {
                        valid_x.push_back(x[indices[j]]);
                        valid_y.push_back(y[indices[j]]);
                    } else {
                        train_x.push_back(x[indices[j]]);
                        train_y.push_back(y[indices[j]]);
                    }
                }

                // Sort training data for efficient interpolation
                std::vector<size_t> train_sort_idx(train_x.size());
                std::iota(train_sort_idx.begin(), train_sort_idx.end(), 0);
                std::sort(train_sort_idx.begin(), train_sort_idx.end(),
                          [&train_x](size_t i1, size_t i2) { return train_x[i1] < train_x[i2]; });

                std::vector<double> sorted_train_x, sorted_train_y;
                for (size_t idx : train_sort_idx) {
                    sorted_train_x.push_back(train_x[idx]);
                    sorted_train_y.push_back(train_y[idx]);
                }

                // Fit on training set
                //auto [train_preds, _] = fit_local_logistic(sorted_train_x,
                auto ll = fit_local_logistic(sorted_train_x,
                                                           sorted_train_y,
                                                           result.candidate_bandwidths[i],
                                                           fit_quadratic,
                                                           kernel_type,
                                                           min_points,
                                                           max_iterations,
                                                           ridge_lambda,
                                                           tolerance,
                                                           false);

                auto train_preds = std::move(ll.preds);

                // Compute validation error using linear interpolation
                for (size_t j = 0; j < valid_x.size(); j++) {
                    // Find the bracketing interval in the sorted training data
                    auto upper = std::upper_bound(sorted_train_x.begin(), sorted_train_x.end(), valid_x[j]);

                    // Handle edge cases
                    if (upper == sorted_train_x.begin()) {
                        // Validation point is before first training point - use nearest neighbor
                        if (!std::isnan(train_preds[0])) {
                            total_error += std::abs(valid_y[j] - train_preds[0]);
                            valid_count++;
                        }
                        continue;
                    }

                    if (upper == sorted_train_x.end()) {
                        // Validation point is after last training point - use nearest neighbor
                        size_t last_idx = sorted_train_x.size() - 1;
                        if (!std::isnan(train_preds[last_idx])) {
                            total_error += std::abs(valid_y[j] - train_preds[last_idx]);
                            valid_count++;
                        }
                        continue;
                    }

                    // Get indices for the bracketing interval
                    size_t upper_idx = std::distance(sorted_train_x.begin(), upper);
                    size_t lower_idx = upper_idx - 1;

                    // Calculate interpolation weight (lambda)
                    double x_lower = sorted_train_x[lower_idx];
                    double x_upper = sorted_train_x[upper_idx];
                    double lambda = (valid_x[j] - x_lower) / (x_upper - x_lower);

                    // Perform linear interpolation of predictions
                    if (!std::isnan(train_preds[lower_idx]) && !std::isnan(train_preds[upper_idx])) {
                        double pred_j = (1.0 - lambda) * train_preds[lower_idx] +
                            lambda * train_preds[upper_idx];
                        total_error += std::abs(valid_y[j] - pred_j);
                        valid_count++;
                    }
                }
            }

            result.mean_errors[i] = valid_count > 0 ? total_error / valid_count :
                std::numeric_limits<double>::infinity();
        }
    }

    // Find optimal bandwidth
    auto min_it = std::min_element(result.mean_errors.begin(), result.mean_errors.end());
    int opt_idx = std::distance(result.mean_errors.begin(), min_it);
    result.opt_bw = result.candidate_bandwidths[opt_idx];

    // NOTE: if bw_predictions is used, the optimal model is already estimated

    // Compute final predictions using optimal bandwidth
    auto ll = fit_local_logistic(x,
                                 y,
                                 result.opt_bw,
                                 fit_quadratic,
                                 kernel_type,
                                 min_points,
                                 max_iterations,
                                 ridge_lambda,
                                 tolerance,
                                 with_errors);

    result.predictions = std::move(ll.preds);
    result.errors = std::move(ll.errs);

    if (fit_quadratic) {
        result.beta2s = std::move(ll.beta2s);
    } else {
        result.beta1s = std::move(ll.beta1s);

        // Compute numerical derivatives for final model
        std::vector<size_t> sort_idx(x.size());
        std::iota(sort_idx.begin(), sort_idx.end(), 0);
        std::sort(sort_idx.begin(), sort_idx.end(),
                  [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });

        std::vector<double> x_sorted(x.size());
        std::vector<double> beta1s_sorted(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            x_sorted[i] = x[sort_idx[i]];
            beta1s_sorted[i] = result.beta1s[sort_idx[i]];
        }

        result.beta2s.resize(x.size(), std::numeric_limits<double>::quiet_NaN());

        for (size_t i = 0; i < x_sorted.size(); i++) {
            if (i == 0) {
                if (i + 1 < x_sorted.size()) {
                    result.beta2s[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i]) /
                        (x_sorted[i + 1] - x_sorted[i]);
                }
            } else if (i == x_sorted.size() - 1) {
                result.beta2s[sort_idx[i]] = (beta1s_sorted[i] - beta1s_sorted[i - 1]) /
                    (x_sorted[i] - x_sorted[i - 1]);
            } else {
                result.beta2s[sort_idx[i]] = (beta1s_sorted[i + 1] - beta1s_sorted[i - 1]) /
                    (x_sorted[i + 1] - x_sorted[i - 1]);
            }
        }
    }

    // Clean up ANN allocations before returning
    annDeallocPts(data_points);
    delete kdtree;
    annClose();

    return result;
}
#endif
