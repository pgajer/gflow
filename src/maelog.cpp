#include <R.h>
#include <Rinternals.h>
#undef length  // to resolve naming conflict between the R macro length defined in Rinternals.h and a member function in the C++ standard library's codecvt class
#undef eval

#include <execution>
#include <mutex>
#include <vector>
#include <numeric>    // for std::iota
#include <random>     // for std::mt19937
#include <algorithm>  // for std::shuffle
#include <ANN/ANN.h>

#include <Eigen/Dense>

#include "ulogit.hpp"
#include "kernels.h"
#include "maelog.hpp"

#include "cpp_utils.hpp" // for debugging and elapsed.time
// #include <fstream>       // for debugging and elapsed.time

maelog_t maelog_mp(
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
    bool verbose);

double compute_std_dev(const std::vector<double>& x);
double compute_iqr(const std::vector<double>& x);

extern "C" {
    SEXP S_maelog(
        SEXP x_r,
        SEXP y_r,
        SEXP fit_quadratic_r,
        SEXP pilot_bandwidth_r,
        SEXP kernel_type_r,
        SEXP min_points_r,
        SEXP cv_folds_r,
        SEXP n_bws_r,
        SEXP min_bw_factor_r,
        SEXP max_bw_factor_r,
        SEXP max_iterations_r,
        SEXP ridge_lambda_r,
        SEXP tolerance_r,
        SEXP with_errors_r,
        SEXP with_bw_preditions_r
        //SEXP parallel_r,
        //SEXP verbose_r
        );
}

/**
 * @brief Selects bandwidth for local logistic regression
 *
 * Provides a rule-of-thumb method for bandwidth selection
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param kernel_type Integer specifying kernel function
 * @return Selected bandwidth value
 */
double rule_of_thumb_logit_bw(
    const std::vector<double>& x,
    int kernel_type) {

    // Modified constants for logistic regression
    const double NORMAL_LOGIT_CONST = 1.4;      // Increased from 0.9
    const double EPANECHNIKOV_LOGIT_CONST = 3.2; // Increased from 2.34
    const double BIWEIGHT_LOGIT_CONST = 3.8;    // Increased from 2.78
    const double TRIANGULAR_LOGIT_CONST = 3.5;  // Increased from 2.58
    const double TRICUBE_LOGIT_CONST = 4.2;     // Increased from 3.12
    const double LAPLACE_LOGIT_CONST = 2.0;     // Increased from 1.30

    double kernel_const;
    switch(kernel_type) {
        case 1: kernel_const = EPANECHNIKOV_LOGIT_CONST; break;
        case 2: kernel_const = TRIANGULAR_LOGIT_CONST; break;
        case 4: kernel_const = LAPLACE_LOGIT_CONST; break;
        case 5: kernel_const = NORMAL_LOGIT_CONST; break;
        case 6: kernel_const = BIWEIGHT_LOGIT_CONST; break;
        case 7: kernel_const = TRICUBE_LOGIT_CONST; break;
        default: kernel_const = NORMAL_LOGIT_CONST;
    }

    int n = x.size();

    // Rule of thumb bandwidth (modified for logistic regression)
    double sd = compute_std_dev(x);
    double iqr = compute_iqr(x);
    double iqr_based = (iqr < std::numeric_limits<double>::epsilon()) ? sd : iqr / 1.34;
    double min_spread = std::min(sd, iqr_based);
    double h_rot = kernel_const * min_spread * std::pow(n, -2.0/7.0);  // Changed exponent

    return h_rot;
}

/**
 * @brief Performs local linear or quadratic logistic regression
 *
 * @details This function implements a local logistic regression algorithm with CV-based global bandwidth
 * selection using Brier score. The algorithm fits local linear or quadratic logistic
 * models in the neighborhood of each data point, with the neighborhood size determined by either
 * a fixed bandwidth or one selected through cross-validation.
 *
 * The function supports three different error measures for bandwidth selection:
 * 1. Deviance-based criterion: $$-[y_i\log(\hat{p}_{(-i)}) + (1-y_i)\log(1-\hat{p}_{(-i)})]$$
 * 2. Brier score: $$(\hat{p}_{(-i)} - y_i)^2$$
 * 3. Absolute error: $$|y_i - \hat{p}_{(-i)}|$$
 *
 * where $$\hat{p}_{(-i)}$$ represents the leave-one-out prediction for observation i.
 *
 * The algorithm proceeds as follows:
 * 1. If pilot_bandwidth > 0, uses this fixed bandwidth for all local models
 * 2. Otherwise:
 *    a. Creates a grid of candidate bandwidths
 *    b. For each bandwidth:
 *       - Fits local models using either LOOCV or k-fold CV
 *       - Computes all three error measures
 *    c. Selects optimal bandwidths minimizing each error measure
 * 3. Returns predictions and model parameters for optimal or specified bandwidth(s)
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param fit_quadratic If true, includes quadratic terms in local models
 * @param pilot_bandwidth If > 0, uses this fixed bandwidth instead of selection
 * @param kernel_type Integer specifying kernel function (1=Gaussian, 2=Epanechnikov)
 * @param min_points Minimum number of points required in local neighborhood
 * @param cv_folds Number of cross-validation folds (0 for LOOCV approximation)
 * @param n_bws Number of bandwidths to test in grid search
 * @param min_bw_factor Lower bound factor for bandwidth grid relative to data range
 * @param max_bw_factor Upper bound factor for bandwidth grid relative to data range
 * @param max_iterations Maximum iterations for logistic regression fitting
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for model fitting
 * @param with_errors If true, computes and returns error measures
 * @param with_bw_predictions If true, retains predictions for all bandwidths
 *
 * @return maelog_t struct containing:
 * - candidate_bandwidths: Vector of tested bandwidths
 * - bw_predictions: Vector of prediction vectors for each bandwidth
 * - mean_deviance_errors: Mean deviance errors for each bandwidth
 * - mean_brier_errors: Mean Brier scores for each bandwidth
 * - mean_abs_errors: Mean absolute errors for each bandwidth
 * - opt_deviance_bw_idx: Index of bandwidth minimizing deviance
 * - opt_brier_bw_idx: Index of bandwidth minimizing Brier score
 * - opt_abs_bw_idx: Index of bandwidth minimizing absolute error
 * - beta1s: Linear coefficients from local models
 * - beta2s: Quadratic coefficients (if fit_quadratic=true)
 *
 * Additional members store input parameters and configuration:
 * - fit_quadratic: Whether quadratic terms were included
 * - pilot_bandwidth: Fixed bandwidth if specified
 * - kernel_type: Type of kernel function used
 * - cv_folds: Number of CV folds used
 * - min_bw_factor: Lower bound factor for bandwidth grid
 * - max_bw_factor: Upper bound factor for bandwidth grid
 * - max_iterations: Maximum fitting iterations
 * - ridge_lambda: Ridge regularization parameter
 * - tolerance: Convergence tolerance
 *
 * @note If with_bw_predictions=false, only predictions for unique optimal
 * bandwidths are retained in bw_predictions to conserve memory.
 *
 * @warning Input vectors x and y must be the same length and y must contain
 * only binary values (0 or 1). The function uses the ANN library for efficient
 * nearest neighbor searches, which must be properly initialized before calling.
 *
 * @see eigen_ulogit_fit For details on the local logistic regression fitting
 * @see initialize_kernel For kernel function initialization
 */
maelog_t maelog(
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
    bool with_bw_predictions) {

    maelog_t result;
    result.fit_quadratic = fit_quadratic;
    result.pilot_bandwidth = pilot_bandwidth;
    result.kernel_type = kernel_type;
    result.cv_folds = cv_folds;
    result.min_bw_factor = min_bw_factor;
    result.max_bw_factor = max_bw_factor;
    result.max_iterations = max_iterations;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    // Lambda function for fitting local logistic regression with a specific bandwidth
    auto fit_local_logistic = [&cv_folds](
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

        int n_pts = x.size();
        // Create data points array for ANN
        ANNpointArray data_points = annAllocPts(n_pts, 1);  // 1 dimension
        for(int i = 0; i < n_pts; i++) {
            data_points[i][0] = x[i];
        }
        // Build kd-tree
        ANNkd_tree* kdtree = new ANNkd_tree(data_points, n_pts, 1);

        std::vector<std::vector<point_t>> pt_pred_weights(n_pts);
        std::vector<double> predictions(n_pts);
        std::vector<double> brier_errors;
        if (with_errors) {
            brier_errors.resize(n_pts);
        }
        std::vector<double> beta1s(n_pts);
        std::vector<double> beta2s;
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
        std::vector<int> local_indices; // Store indices along with local data for accurate mapping
        local_indices.reserve(n_pts);

        // fitting individual models in a neighborhood of bandwidth 'bandwidth' of each point
        for (int pti = 0; pti < n_pts; pti++) {
            double center = x[pti];

            // Clear vectors for reuse
            local_x.clear();
            local_y.clear();
            local_w.clear();
            local_d.clear();
            local_indices.clear();

            // First try bandwidth-based neighborhood
            for (int i = 0; i < n_pts; ++i) {
                double dist = std::abs(x[i] - center) / bandwidth;
                if (dist < 1.0) {
                    local_indices.push_back(i);
                }
            }

            int n_window_pts = local_indices.size();
            if (n_window_pts < min_points) { // If not enough points, switch to k-NN

                query_point[0] = center;
                kdtree->annkSearch(query_point, min_points, nn_idx, nn_dists);

                //Rprintf("pti: %d\tn_window_pts: %d\n", pti, n_window_pts);
                //print_vect(local_indices, "local_indices");

                local_indices.clear();  // Clear existing indices
                // Calculate shifted x values and distances
                for (int i = 0; i < min_points; ++i) {
                    local_indices.push_back(nn_idx[i]);  // Store the new indices
                    double shifted_x = x[nn_idx[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[nn_idx[i]]);
                    local_d.push_back(std::abs(shifted_x));
                }

                //print_vect(local_indices, "After Upadate local_indices");

                // Find maximum distance and normalize
                double max_dist = *std::max_element(local_d.begin(), local_d.end());
                if (max_dist < std::numeric_limits<double>::epsilon()) {
                    max_dist = std::numeric_limits<double>::epsilon();
                }

                for (auto& d : local_d) {
                    d /= max_dist;
                }

                n_window_pts = min_points;

            } else {
                for (int i = 0; i < n_window_pts; ++i) {
                    double shifted_x = x[local_indices[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[local_indices[i]]);
                    local_d.push_back(std::abs(shifted_x) / bandwidth);
                }
            }

            local_w.resize(n_window_pts); // Resize local_w before using it
            kernel_fn(local_d.data(), n_window_pts, local_w.data());

            // Normalize weights with protection against zero sum
            double sum_weights = std::accumulate(local_w.begin(), local_w.end(), 0.0);
            if (sum_weights < std::numeric_limits<double>::epsilon()) {
                sum_weights = std::numeric_limits<double>::epsilon();
            }
            for (auto& w : local_w) {
                w /= sum_weights;
            }

            eigen_ulogit_t fit_result = eigen_ulogit_fit(
                local_x.data(),
                local_y.data(),
                local_w,
                fit_quadratic,
                max_iterations,
                ridge_lambda,
                tolerance,
                with_errors
            );

            // Store predictions and errors
            point_t pt;
            for (size_t i = 0; i < local_x.size(); ++i) {
                int orig_idx = local_indices[i];

                if (orig_idx < 0 || orig_idx > (int)(pt_pred_weights.size() - 1)) {
                    Rprintf("i: %d  orig_idx: %d  is out of range  pt_pred_weights.size(): %d\n", (int)i, orig_idx, (int)pt_pred_weights.size());
                }

                pt.w = local_w[i];
                pt.p = fit_result.predictions[i];
                pt.beta1 = fit_result.beta[1];
                if (fit_quadratic) pt.beta2 = fit_result.beta[2];
                if (with_errors) {
                    pt.brier_error = fit_result.loocv_brier_errors[i];
                }
                pt_pred_weights[orig_idx].push_back(pt);
            }
        }

        // Compute weighted averages
        for (int pti = 0; pti < n_pts; pti++) {
            const auto& v = pt_pred_weights[pti];
            if (v.empty()) {
                predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double total_weight = 0.0;
            double weighted_prediction = 0.0;
            double weighted_brier_error = 0.0;
            double weighted_beta1 = 0.0;
            double weighted_beta2 = 0.0;

            for (const auto& pt : v) {
                weighted_prediction += pt.w * pt.p;
                weighted_beta1 += pt.w * pt.beta1;
                if (with_errors) {
                    weighted_brier_error += pt.w * pt.brier_error;
                }
                if (fit_quadratic) weighted_beta2 += pt.w * pt.beta2;
                total_weight += pt.w;
            }

            if (total_weight > 0.0) {
                predictions[pti] = weighted_prediction / total_weight;
                beta1s[pti] = weighted_beta1 / total_weight;
                if (with_errors) {
                    brier_errors[pti] = weighted_brier_error / total_weight;
                }
                if (fit_quadratic) beta2s[pti] = weighted_beta2 / total_weight;
            } else {
                predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                beta1s[pti] = std::numeric_limits<double>::quiet_NaN();
                if (with_errors) {
                    brier_errors[pti] = std::numeric_limits<double>::quiet_NaN();
                }
                if (fit_quadratic) beta2s[pti] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // Clean up ANN allocations
        annDeallocPt(query_point);
        delete[] nn_idx;
        delete[] nn_dists;
        // Clean up ANN allocations
        annDeallocPts(data_points);
        delete kdtree;
        annClose();

        local_logit_t ll;
        ll.preds = std::move(predictions);
        ll.beta1s = std::move(beta1s);
        if (fit_quadratic) ll.beta2s = std::move(beta2s);
        if (with_errors) {
            ll.brier_errors = std::move(brier_errors);
        }

        return ll;
    };

    // If pilot_bandwidth > 0, use it directly
    if (pilot_bandwidth > 0) {
        auto ll = fit_local_logistic(x, y, pilot_bandwidth, fit_quadratic, kernel_type,
                                   min_points, max_iterations, ridge_lambda, tolerance, with_errors);

        result.beta1s = std::move(ll.beta1s);
        if (fit_quadratic) {
            result.beta2s = std::move(ll.beta2s);
        }

        result.bw_predictions.resize(1);
        result.bw_predictions[0] = std::move(ll.preds);

        result.opt_brier_bw_idx = -1;

        return result;
    }

    // Create bandwidth grid
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

    // Initialize error vectors
    result.mean_brier_errors.resize(n_bws);
    result.bw_predictions.resize(n_bws); // this vector is always initialized, even if with_bw_predictions = false, in which case only optimal bw predictions will be non-empty

    // Perform bandwidth selection
    if (cv_folds == 0) {
        // LOOCV approximation
        for (int i = 0; i < n_bws; i++) {
            auto ll = fit_local_logistic(x, y, result.candidate_bandwidths[i], fit_quadratic,
                                       kernel_type, min_points, max_iterations, ridge_lambda,
                                       tolerance, true);

            // Compute mean errors
            auto compute_mean = [](const std::vector<double>& errors) {
                double total = 0.0;
                int valid_count = 0;
                for (const auto& err : errors) {
                    if (!std::isnan(err)) {
                        total += err;
                        valid_count++;
                    }
                }
                return valid_count > 0 ? total / valid_count :
                                       std::numeric_limits<double>::infinity();
            };

            result.mean_brier_errors[i] = compute_mean(ll.brier_errors);

            result.bw_predictions[i] = std::move(ll.preds); // when with_bw_predictions = false, predictions corresponding to non-optimal bw's are going to be set to empty vectors
        }
    } else {
        // K-fold cross-validation implementation
        int n = x.size();
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(indices.begin(), indices.end(), gen);

        int fold_size = n / cv_folds;
        for (int i = 0; i < n_bws; i++) {
            double total_brier_error    = 0.0;
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
                            total_brier_error    += std::pow(train_preds[0] - valid_y[j], 2);
                            valid_count++;
                        }
                        continue;
                    }

                    if (upper == sorted_train_x.end()) {
                        // Validation point is after last training point - use nearest neighbor
                        size_t last_idx = sorted_train_x.size() - 1;
                        if (!std::isnan(train_preds[last_idx])) {
                            total_brier_error    += std::pow(train_preds[last_idx] - valid_y[j], 2);
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
                        total_brier_error    += std::pow(pred_j - valid_y[j], 2);
                        valid_count++;
                    }
                }
            }

            result.mean_brier_errors[i] = valid_count > 0 ? total_brier_error / valid_count :
                std::numeric_limits<double>::infinity();
        }
    }

    // Find optimal bandwidths for each error measure
    auto find_opt_idx = [](const std::vector<double>& errors) {
        return std::distance(errors.begin(),
                           std::min_element(errors.begin(), errors.end()));
    };

    result.opt_brier_bw_idx = find_opt_idx(result.mean_brier_errors);

    std::set<int> unique_opt_indices = {
        result.opt_brier_bw_idx
    };

    if (cv_folds == 0 && !with_bw_predictions) {
        // Clear non-optimal predictions
        for (size_t i = 0; i < result.bw_predictions.size(); i++) {
            if (unique_opt_indices.find(i) == unique_opt_indices.end()) {
                result.bw_predictions[i].clear();
            }
        }
    } else if (cv_folds > 0 && !with_bw_predictions) {
        // Store coefficients for each optimal bandwidth
        for (int idx : unique_opt_indices) {

            auto ll = fit_local_logistic(x, y, result.candidate_bandwidths[idx],
                                         fit_quadratic, kernel_type, min_points,
                                         max_iterations, ridge_lambda, tolerance, false);
            result.bw_predictions[idx] = std::move(ll.preds);
            // I am leaving it for the future
            if (ll.beta1s.size() > 0) result.beta1s = std::move(ll.beta1s);
            if (ll.beta2s.size() > 0) result.beta2s = std::move(ll.beta2s);
        }
    } else if (cv_folds > 0 && with_bw_predictions) {
        for (size_t i = 0; i < result.candidate_bandwidths.size(); i++) {
            auto ll = fit_local_logistic(x, y, result.candidate_bandwidths[i],
                                         fit_quadratic, kernel_type, min_points,
                                         max_iterations, ridge_lambda, tolerance, false);
            result.bw_predictions[i] = std::move(ll.preds);
            if (unique_opt_indices.find(i) != unique_opt_indices.end()) {
                if (ll.beta1s.size() > 0) result.beta1s = std::move(ll.beta1s);
                if (ll.beta2s.size() > 0) result.beta2s = std::move(ll.beta2s);
            }
        }
    }

    annClose();

    return result;
}


/**
 * @brief R interface for model averaged bandwidth logistic regression
 *
 * This function provides an R interface to the C++ implementation of model averaged
 * bandwidth logistic regression (maelog). It handles conversion between R and C++
 * data types, performs input validation, and manages R's garbage collection protection.
 *
 * @param x_r SEXP (NumericVector) Predictor variable values
 * @param y_r SEXP (NumericVector) Binary response variable values (0 or 1)
 * @param fit_quadratic_r SEXP (LogicalVector) Whether to include quadratic terms
 * @param pilot_bandwidth_r SEXP (NumericVector) Initial bandwidth value (if > 0, fixed bandwidth is used)
 * @param kernel_type_r SEXP (IntegerVector) Kernel function type (1-7)
 * @param min_points_r SEXP (IntegerVector) Minimum number of points required for local fitting
 * @param cv_folds_r SEXP (IntegerVector) Number of cross-validation folds (0 for LOOCV)
 * @param n_bws_r SEXP (IntegerVector) Number of bandwidth values to evaluate
 * @param min_bw_factor_r SEXP (NumericVector) Minimum bandwidth factor relative to data range
 * @param max_bw_factor_r SEXP (NumericVector) Maximum bandwidth factor relative to data range
 * @param max_iterations_r SEXP (IntegerVector) Maximum number of iterations for model fitting
 * @param ridge_lambda_r SEXP (NumericVector) Ridge regularization parameter
 * @param tolerance_r SEXP (NumericVector) Convergence tolerance for model fitting
 * @param with_errors_r SEXP (LogicalVector) Whether to compute and return prediction errors
 * @param with_bw_preditions_r SEXP (LogicalVector) Whether to return predictions for each bandwidth
 * @param verbose_r SEXP (LogicalVector) Whether to display progress information
 *
 * @return SEXP (List) A named list containing:
 *   - bw_predictions: Vector of prediction vectors, one per bandwidth. If with_bw_predictions=FALSE,
 *     only predictions for optimal bandwidths are kept (non-empty vectors).
 *   - mean_deviance_errors: Mean deviance errors for each candidate bandwidth
 *   - mean_brier_errors: Mean Brier errors for each candidate bandwidth
 *   - mean_abs_errors: Mean absolute errors for each candidate bandwidth
 *   - opt_deviance_bw_idx: Index of bandwidth minimizing deviance error
 *   - opt_brier_bw_idx: Index of bandwidth minimizing Brier error
 *   - opt_abs_bw_idx: Index of bandwidth minimizing absolute error
 *   - beta1s: Linear coefficients (only when pilot_bandwidth > 0)
 *   - beta2s: Quadratic coefficients (only when pilot_bandwidth > 0 and fit_quadratic=TRUE)
 *   - fit_info: List containing model parameters:
 *     - fit_quadratic: Whether quadratic terms were included
 *     - pilot_bandwidth: Fixed bandwidth if specified
 *     - kernel_type: Type of kernel function used
 *     - cv_folds: Number of CV folds used
 *     - min_bw_factor: Lower bound factor for bandwidth grid
 *     - max_bw_factor: Upper bound factor for bandwidth grid
 *     - max_iterations: Maximum fitting iterations
 *     - ridge_lambda: Ridge regularization parameter
 *     - tolerance: Convergence tolerance
 *
 * @throws R error if inputs are invalid or computation fails
 *
 * @note Protected SEXP objects are properly unprotected before return
 * @note Beta coefficients are only returned when using a fixed pilot bandwidth
 */
SEXP S_maelog(
    SEXP x_r,
    SEXP y_r,
    SEXP fit_quadratic_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP min_points_r,
    SEXP cv_folds_r,
    SEXP n_bws_r,
    SEXP min_bw_factor_r,
    SEXP max_bw_factor_r,
    SEXP max_iterations_r,
    SEXP ridge_lambda_r,
    SEXP tolerance_r,
    SEXP with_errors_r,
    SEXP with_bw_preditions_r
    //SEXP parallel_r,
    //SEXP verbose_r
    ) {

    int n_protected = 0;

    try {
        // Check for NULL inputs
        if (x_r == R_NilValue || y_r == R_NilValue || fit_quadratic_r == R_NilValue ||
            pilot_bandwidth_r == R_NilValue || kernel_type_r == R_NilValue ||
            min_points_r == R_NilValue || cv_folds_r == R_NilValue ||
            n_bws_r == R_NilValue || min_bw_factor_r == R_NilValue ||
            max_bw_factor_r == R_NilValue || max_iterations_r == R_NilValue ||
            ridge_lambda_r == R_NilValue ||
            tolerance_r == R_NilValue || with_errors_r == R_NilValue ||
            with_bw_preditions_r == R_NilValue) {
            Rf_error("Input arguments cannot be NULL");
        }

        int n_points = LENGTH(x_r);

        // Convert inputs
        std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
        std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);
        bool fit_quadratic = LOGICAL(fit_quadratic_r)[0];
        double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
        int kernel_type = INTEGER(kernel_type_r)[0];
        int min_points = INTEGER(min_points_r)[0];
        int cv_folds = INTEGER(cv_folds_r)[0];
        int n_bws = INTEGER(n_bws_r)[0];
        double min_bw_factor = REAL(min_bw_factor_r)[0];
        double max_bw_factor = REAL(max_bw_factor_r)[0];
        int max_iterations = INTEGER(max_iterations_r)[0];
        double ridge_lambda = REAL(ridge_lambda_r)[0];
        double tolerance = REAL(tolerance_r)[0];
        bool with_errors = LOGICAL(with_errors_r)[0];
        bool with_bw_preditions = LOGICAL(with_bw_preditions_r)[0];
        //bool parallel = LOGICAL(parallel_r)[0];
        //bool verbose = LOGICAL(verbose_r)[0];

        auto result = maelog(x, y, fit_quadratic, pilot_bandwidth, kernel_type,
                                min_points, cv_folds, n_bws, min_bw_factor, max_bw_factor,
                                max_iterations, ridge_lambda, tolerance,
                                with_errors, with_bw_preditions);

        // Call appropriate implementation
        #if 0
        auto result = parallel ?
            maelog_mp(x, y, fit_quadratic, pilot_bandwidth, kernel_type,
                         min_points, cv_folds, n_bws, min_bw_factor, max_bw_factor,
                         max_iterations, ridge_lambda, tolerance,
                         with_errors, with_bw_preditions, verbose) :
            maelog(x, y, fit_quadratic, pilot_bandwidth, kernel_type,
                      min_points, cv_folds, n_bws, min_bw_factor, max_bw_factor,
                      max_iterations, ridge_lambda, tolerance,
                      with_errors, with_bw_preditions);
        #endif

        // Create return list with updated size
        SEXP result_list = PROTECT(Rf_allocVector(VECSXP, 7)); n_protected++;

        // Updated names for list elements
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 7)); n_protected++;
        SET_STRING_ELT(names, 0, Rf_mkChar("bw_predictions"));
        SET_STRING_ELT(names, 1, Rf_mkChar("mean_brier_errors"));
        SET_STRING_ELT(names, 2, Rf_mkChar("opt_brier_bw_idx"));
        SET_STRING_ELT(names, 3, Rf_mkChar("beta1s"));
        SET_STRING_ELT(names, 4, Rf_mkChar("beta2s"));
        SET_STRING_ELT(names, 5, Rf_mkChar("bws")); // candidate_bandwidths
        SET_STRING_ELT(names, 6, Rf_mkChar("fit_info"));

        // Bandwidth predictions
        SEXP bw_pred_r;
        if (!result.bw_predictions.empty()) {
            bw_pred_r = PROTECT(Rf_allocMatrix(REALSXP, n_points, result.bw_predictions.size())); n_protected++;
            double* ptr = REAL(bw_pred_r);
            for (size_t i = 0; i < result.bw_predictions.size(); ++i) {
                if (!result.bw_predictions[i].empty()) {
                    std::copy(result.bw_predictions[i].begin(), result.bw_predictions[i].end(),
                             ptr + i * n_points);
                }
            }
        } else {
            bw_pred_r = PROTECT(Rf_allocVector(REALSXP, 0)); n_protected++;
        }
        SET_VECTOR_ELT(result_list, 0, bw_pred_r);

        // Mean errors
        SEXP mean_brier_errors_r = PROTECT(Rf_allocVector(REALSXP, result.mean_brier_errors.size())); n_protected++;
        std::copy(result.mean_brier_errors.begin(), result.mean_brier_errors.end(), REAL(mean_brier_errors_r));
        SET_VECTOR_ELT(result_list, 1, mean_brier_errors_r);

        // Optimal indices
        SEXP opt_brier_bw_idx_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
        INTEGER(opt_brier_bw_idx_r)[0] = result.opt_brier_bw_idx + 1; // 1-base
        SET_VECTOR_ELT(result_list, 2, opt_brier_bw_idx_r);

        // Beta coefficients (only when pilot_bandwidth > 0)
        SEXP beta1s_r;
        if (pilot_bandwidth > 0 && !result.beta1s.empty()) {
            beta1s_r = PROTECT(Rf_allocVector(REALSXP, result.beta1s.size())); n_protected++;
            std::copy(result.beta1s.begin(), result.beta1s.end(), REAL(beta1s_r));
        } else {
            beta1s_r = PROTECT(Rf_allocVector(REALSXP, 0)); n_protected++;
        }
        SET_VECTOR_ELT(result_list, 3, beta1s_r);

        SEXP beta2s_r;
        if (pilot_bandwidth > 0 && fit_quadratic && !result.beta2s.empty()) {
            beta2s_r = PROTECT(Rf_allocVector(REALSXP, result.beta2s.size())); n_protected++;
            std::copy(result.beta2s.begin(), result.beta2s.end(), REAL(beta2s_r));
        } else {
            beta2s_r = PROTECT(Rf_allocVector(REALSXP, 0)); n_protected++;
        }
        SET_VECTOR_ELT(result_list, 4, beta2s_r);

        // Candidate bandwidths
        SEXP candidate_bandwidths_r = PROTECT(Rf_allocVector(REALSXP, result.candidate_bandwidths.size())); n_protected++;
        std::copy(result.candidate_bandwidths.begin(), result.candidate_bandwidths.end(),
                 REAL(candidate_bandwidths_r));
        SET_VECTOR_ELT(result_list, 5, candidate_bandwidths_r);

        // Fit info as a named list
        SEXP fit_info = PROTECT(Rf_allocVector(VECSXP, 9)); n_protected++;

        SEXP fit_info_names = PROTECT(Rf_allocVector(STRSXP, 9)); n_protected++;
        SET_STRING_ELT(fit_info_names, 0, Rf_mkChar("fit_quadratic"));
        SET_STRING_ELT(fit_info_names, 1, Rf_mkChar("pilot_bandwidth"));
        SET_STRING_ELT(fit_info_names, 2, Rf_mkChar("kernel_type"));
        SET_STRING_ELT(fit_info_names, 3, Rf_mkChar("cv_folds"));
        SET_STRING_ELT(fit_info_names, 4, Rf_mkChar("min_bw_factor"));
        SET_STRING_ELT(fit_info_names, 5, Rf_mkChar("max_bw_factor"));
        SET_STRING_ELT(fit_info_names, 6, Rf_mkChar("max_iterations"));
        SET_STRING_ELT(fit_info_names, 7, Rf_mkChar("ridge_lambda"));
        SET_STRING_ELT(fit_info_names, 8, Rf_mkChar("tolerance"));

        SEXP fit_quad_r = PROTECT(Rf_allocVector(LGLSXP, 1));
        LOGICAL(fit_quad_r)[0] = result.fit_quadratic;
        SET_VECTOR_ELT(fit_info, 0, fit_quad_r); n_protected++;

        SEXP pilot_bw_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(pilot_bw_r)[0] = result.pilot_bandwidth;
        SET_VECTOR_ELT(fit_info, 1, pilot_bw_r); n_protected++;

        SEXP kernel_type_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(kernel_type_out_r)[0] = result.kernel_type;
        SET_VECTOR_ELT(fit_info, 2, kernel_type_out_r); n_protected++;

        SEXP cv_folds_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(cv_folds_out_r)[0] = result.cv_folds;
        SET_VECTOR_ELT(fit_info, 3, cv_folds_out_r); n_protected++;

        SEXP min_bw_factor_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(min_bw_factor_out_r)[0] = result.min_bw_factor;
        SET_VECTOR_ELT(fit_info, 4, min_bw_factor_out_r); n_protected++;

        SEXP max_bw_factor_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(max_bw_factor_out_r)[0] = result.max_bw_factor;
        SET_VECTOR_ELT(fit_info, 5, max_bw_factor_out_r); n_protected++;

        SEXP max_iter_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(max_iter_out_r)[0] = result.max_iterations;
        SET_VECTOR_ELT(fit_info, 6, max_iter_out_r); n_protected++;

        SEXP ridge_lambda_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(ridge_lambda_out_r)[0] = result.ridge_lambda;
        SET_VECTOR_ELT(fit_info, 7, ridge_lambda_out_r); n_protected++;

        SEXP tolerance_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(tolerance_out_r)[0] = result.tolerance;
        SET_VECTOR_ELT(fit_info, 8, tolerance_out_r); n_protected++;

        Rf_setAttrib(fit_info, R_NamesSymbol, fit_info_names);
        SET_VECTOR_ELT(result_list, 6, fit_info);

        // Set names for the main list
        Rf_setAttrib(result_list, R_NamesSymbol, names);

        UNPROTECT(n_protected);
        return result_list;
    }
    catch (const std::exception& e) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("C++ error in maelog: %s", e.what());
    }
    catch (...) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("Unknown error in maelog");
    }
}
