#include <R.h>
#include <Rinternals.h>

#include <Eigen/Dense>
#include <vector>
#include <numeric>    // for std::iota
#include <random>     // for std::mt19937
#include <algorithm>  // for std::shuffle

#include "ulogit.hpp"
#include "kernels.h"

double compute_std_dev(const std::vector<double>& x);
double compute_iqr(const std::vector<double>& x);

extern "C" {
    SEXP S_mabwlogit(
        SEXP x_r,
        SEXP y_r,
        SEXP fit_quadratic_r,
        SEXP pilot_bandwidth_r,
        SEXP kernel_type_r,
        SEXP cv_folds_r);
}

/**
 * @brief Selects bandwidth for local logistic regression
 *
 * Provides two methods for bandwidth selection:
 * 1. Modified rule-of-thumb based on logistic regression theory
 * 2. Cross-validation based selection (if cv_folds > 0)
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param kernel_type Integer specifying kernel function
 * @param cv_folds Number of cross-validation folds (0 for rule-of-thumb)
 * @return Selected bandwidth value
 */
double logistic_bandwidth_select(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int kernel_type,
    int cv_folds = 0,
    int n_bws = 100,
    double min_bw_factor = 0.5,
    double max_bw_factor = 2.0) {

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

    // Return rule-of-thumb if no cross-validation requested
    if (cv_folds <= 0) {
        return h_rot;
    }

    // Cross-validation based selection
    // Creating a uniform grid of bandwidths
    std::vector<double> candidate_bandwidths(n_bws);
    double start = h_rot * min_bw_factor;
    double end   = h_rot * max_bw_factor;
    double dx = (end - start) / (n_bws - 1);
    for(int i = 0; i < n_bws; i++)
        candidate_bandwidths[i] = start + i * dx;

    std::vector<double> cv_errors(n_bws);
    // double best_bandwidth = h_rot;
    // double min_error = std::numeric_limits<double>::infinity();

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Cross-validation loop
    for (double h : candidate_bandwidths) {
        double cv_error = 0.0;

        // Create fold indices
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);  // Modern replacement for random_shuffle

        for (int fold = 0; fold < cv_folds; ++fold) {
            std::vector<double> train_x, train_y, test_x, test_y, weights;

            // Split data into training and validation
            for (int i = 0; i < n; ++i) {
                if (i % cv_folds == fold) {
                    test_x.push_back(x[indices[i]]);
                    test_y.push_back(y[indices[i]]);
                } else {
                    train_x.push_back(x[indices[i]]);
                    train_y.push_back(y[indices[i]]);
                }
            }

            // Compute predictions for test points
            for (size_t i = 0; i < test_x.size(); ++i) {
                // Compute kernel weights
                weights.clear();
                for (double train_point : train_x) {
                    double dist = std::abs(test_x[i] - train_point) / h;
                    double weight;
                    kernel_fn(&dist, 1, &weight);
                    weights.push_back(weight);
                }

                // Fit local model and get prediction
                eigen_ulogit_t fit = eigen_ulogit_fit(
                    train_x.data(),
                    train_y.data(),
                    weights,
                    true
                );

                if (fit.converged) {
                    // Compute prediction error
                    double pred = fit.predictions[0];  // Prediction for test point
                    cv_error += -test_y[i] * std::log(pred + 1e-10) -
                               (1 - test_y[i]) * std::log(1 - pred + 1e-10);
                }
            }
        }

        // Update best bandwidth if current error is lower
        if (cv_error < min_error) {
            min_error = cv_error;
            best_bandwidth = h;
        }
    }

    return best_bandwidth;
}


/**
 * @brief Performs model-averaged bandwidth logistic regression
 *
 * This function implements a model-averaged approach to local logistic regression,
 * using multiple bandwidth scales to capture both local and global patterns in
 * the data. For each point, it:
 * 1. Fits local logistic models using eigen_ulogit_fit()
 * 2. Computes weighted predictions based on kernel weights
 * 3. At each point averages predictions using all models with the given point in their support
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param fit_quadratic If true, fits quadratic term in addition to linear
 * @param pilot_bandwidth Initial bandwidth for local fitting (if <= 0, automatically selected)
 * @param kernel_type Integer specifying the kernel function type
 * @return Vector of model-averaged predicted probabilities for each observation
 * @throw std::invalid_argument if input vectors are empty or of unequal length
 * @throw std::runtime_error if no successful model fits are achieved
 */
std::vector<double> mabwlogit(
    const std::vector<double>& x,
    const std::vector<double>& y,
    bool fit_quadratic,
    double pilot_bandwidth,
    int kernel_type,
    int cv_folds = 0) {

    if (x.empty() || y.empty()) {
        error("Invalid input parameters");
    }

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = logistic_bandwidth_select(x, y, kernel_type, cv_folds);
    }

    initialize_kernel(kernel_type, 1.0);
    int n_pts = x.size();

    // Pre-allocate vectors
    std::vector<double> local_x, local_y, local_w;
    local_x.reserve(n_pts);
    local_y.reserve(n_pts);
    local_w.reserve(n_pts);

    // Storage for predictions and weights of every model with the given point in its support
    std::vector<std::vector<std::pair<double,double>>> pt_pred_weights(n_pts);

    bool any_successful_fit = false;  // Track if we get any fits at all

    // Loop over each point as the center of a local neighborhood
    for (int pti = 0; pti < n_pts; pti++) {
        double center = x[pti];

        // Clear vectors for reuse
        local_x.clear();
        local_y.clear();
        local_w.clear();

        // Store indices along with local data for accurate mapping
        std::vector<int> local_indices;
        local_indices.reserve(n_pts);

        // Collect local points and compute weights
        double sum_weights = 0.0;
        for (int i = 0; i < n_pts; ++i) {
            double dist = std::abs(x[i] - center) / pilot_bandwidth;
            if (dist < 1.0) {
                local_x.push_back(x[i] - center);
                local_y.push_back(y[i]);
                local_indices.push_back(i);

                double weight;
                kernel_fn(&dist, 1, &weight);
                local_w.push_back(weight);
                sum_weights += weight;
            }
        }

        // Skip if too few points for stable estimation
        if (local_x.size() < 6) {
            continue;
        }

        // Normalize weights
        for (auto& w : local_w) {
            w /= sum_weights;
        }

        // Fit local logistic model using eigen_ulogit_fit
        eigen_ulogit_t fit_result = eigen_ulogit_fit(
            local_x.data(),
            local_y.data(),
            local_w,
            fit_quadratic
        );

        if (fit_result.converged) {
            any_successful_fit = true;
        }

        // Store predictions and weights for each point in the neighborhood
        for (size_t i = 0; i < local_x.size(); ++i) {
            int orig_idx = local_indices[i];
            pt_pred_weights[orig_idx].push_back(
                std::make_pair(fit_result.predictions[i], local_w[i])
            );
        }
    }

    // Warn if no fits converged
    if (!any_successful_fit) {
        warning("No fits achieved full convergence");
    }

    // Compute weighted averages
    std::vector<double> predictions(n_pts);
    for (int pti = 0; pti < n_pts; pti++) {
        const auto& v = pt_pred_weights[pti];
        if (v.empty()) {
            predictions[pti] = std::numeric_limits<double>::quiet_NaN();  // No valid estimates
            continue;
        }

        double total_weight = 0.0;
        double weighted_sum = 0.0;
        for (const auto& p : v) {
            weighted_sum += p.first * p.second;
            total_weight += p.second;
        }

        predictions[pti] = total_weight > 0.0 ? weighted_sum / total_weight :
                                              std::numeric_limits<double>::quiet_NaN();
    }

    return predictions;
}

/**
 * @brief R interface for model averaged bandwidth logistic regression
 *
 * This function provides an R interface to the C++ implementation of model averaged
 * bandwidth logistic regression (mabwlogit). It handles conversion between R and C++
 * data types, performs input validation, and manages R's garbage collection protection.
 *
 * @param x_r SEXP (NumericVector) Predictor variable values
 * @param y_r SEXP (NumericVector) Binary response variable values (0 or 1)
 * @param fit_quadratic_r SEXP (LogicalVector) Whether to include quadratic terms
 * @param pilot_bandwidth_r SEXP (NumericVector) Initial bandwidth value (if <= 0, automatically selected)
 * @param kernel_type_r SEXP (IntegerVector) Kernel function type (1-7)
 * @param cv_folds_r SEXP (IntegerVector) Number of cross-validation folds (0 for no CV)
 * @return SEXP (NumericVector) Vector of predicted probabilities
 * @throws R error if inputs are invalid or computation fails
 *
 * @note The function uses R's error handling system via Rf_error()
 * @note Protected SEXP objects are properly unprotected before return
 */
SEXP S_mabwlogit(
    SEXP x_r,
    SEXP y_r,
    SEXP fit_quadratic_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP cv_folds_r) {

    // Ensure we're running in R's protection stack
    int n_protected = 0;

    try {
        // Check for NULL inputs
        if (x_r == R_NilValue || y_r == R_NilValue || fit_quadratic_r == R_NilValue ||
            pilot_bandwidth_r == R_NilValue || kernel_type_r == R_NilValue || cv_folds_r == R_NilValue) {
            Rf_error("Input arguments cannot be NULL");
        }

        int n_points = LENGTH(x_r);

        // Convert inputs to C++ vectors
        std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
        std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);
        bool fit_quadratic = LOGICAL(fit_quadratic_r)[0];
        double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
        int kernel_type = INTEGER(kernel_type_r)[0];
        int cv_folds = INTEGER(cv_folds_r)[0];

        // Call the C++ implementation
        auto result = mabwlogit(x, y, fit_quadratic, pilot_bandwidth, kernel_type, cv_folds);

        // Allocate result vector
        SEXP result_r = PROTECT(Rf_allocVector(REALSXP, n_points));
        n_protected++;

        // Copy results back to R vector
        std::copy(result.begin(), result.end(), REAL(result_r));

        UNPROTECT(n_protected);
        return result_r;
    }
    catch (const std::exception& e) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("C++ error in mabwlogit: %s", e.what());
    }
    catch (...) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("Unknown error in mabwlogit");
    }
}
