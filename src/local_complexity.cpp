#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "density.hpp"
#include "kernels.h"
#include "ulogit.hpp"

extern "C" {
    SEXP S_estimate_local_complexity(
        SEXP x_r,
        SEXP y_r,
        SEXP center_idx_r,
        SEXP pilot_bandwidth_r,
        SEXP kernel_type_r);

    SEXP S_estimate_binary_local_complexity(
        SEXP x_r,
        SEXP y_r,
        SEXP center_idx_r,
        SEXP pilot_bandwidth_r,
        SEXP kernel_type_r,
        SEXP method_r);

    SEXP S_estimate_ma_binary_local_complexity_quadratic(SEXP x_r,
                                                         SEXP y_r,
                                                         SEXP pilot_bandwidth_r,
                                                         SEXP kernel_type_r);
}


/**
 * @brief Estimates the local complexity of data at a specified point using weighted quadratic regression
 *
 * This function fits a local quadratic model around a center point to estimate the
 * local complexity of the relationship between x and y. The complexity is measured
 * by the absolute value of the quadratic coefficient, scaled by the bandwidth.
 *
 * @param x Vector of predictor values
 * @param y Vector of response values
 * @param center_idx Index of the center point in x around which to estimate complexity
 * @param pilot_bandwidth Bandwidth parameter for the local window. If <= 0, a default
 *                       bandwidth will be computed using kernel_specific_bandwidth()
 *
 * @return The estimated local complexity (absolute value of scaled quadratic coefficient)
 *
 * @throws std::invalid_argument If input vectors are empty, of unequal length, or if
 *                              center_idx is out of bounds
 *
 * @note The function uses kernel weights defined by kernel_fn() and requires at least
 *       4 points within the bandwidth window for stable estimation
 *
 * @see kernel_specific_bandwidth(), kernel_fn()
 */
double estimate_local_complexity(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth,
    int kernel_type) {

    if (x.empty() || y.empty() || x.size() != y.size() ||
        center_idx < 0 || center_idx >= static_cast<int>(x.size())) {
        Rf_error("Invalid input parameters");
    }

    double center = x[center_idx];
    int n = x.size();

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = kernel_specific_bandwidth(x, kernel_type);
    }

    initialize_kernel(kernel_type, 1.0);

    // Pre-allocate vectors with estimated size
    std::vector<double> local_x, local_y, local_w;
    local_x.reserve(n);
    local_y.reserve(n);
    local_w.reserve(n);

    // Collect local points and compute weights in one pass
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {  // Within bandwidth window
            local_x.push_back((x[i] - center) / pilot_bandwidth);  // Scale x values
            local_y.push_back(y[i]);

            double weight;
            kernel_fn(&dist, 1, &weight);  // Compute kernel weight
            local_w.push_back(std::sqrt(weight));  // Pre-compute square root of weights
        }
    }

    // Check if we have enough points for meaningful estimation
    if (local_x.size() < 4) {  // Require more points for stable estimation
        return 0.0;
    }

    // Construct design matrix with scaled x values
    Eigen::MatrixXd X(local_x.size(), 3);
    for (size_t i = 0; i < local_x.size(); ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = local_x[i];
        X(i, 2) = local_x[i] * local_x[i];
    }

    // Map vectors to Eigen types
    Eigen::VectorXd Y = Eigen::Map<Eigen::VectorXd>(local_y.data(), local_y.size());
    Eigen::VectorXd W = Eigen::Map<Eigen::VectorXd>(local_w.data(), local_w.size());

    // Compute weighted least squares with QR decomposition for better stability
    // Note: We now use the pre-computed square root of weights
    Eigen::MatrixXd WX = W.asDiagonal() * X;
    Eigen::VectorXd WY = W.asDiagonal() * Y;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(WX);
    if (qr.rank() < 3) {
        return 0.0;  // Not enough independent columns for quadratic fit
    }

    Eigen::Vector3d beta = qr.solve(WY);

    // Return scaled quadratic coefficient
    return std::abs(beta(2)) / (pilot_bandwidth * pilot_bandwidth);
}


/**
 * @brief R interface for local complexity estimation
 *
 * This function provides an interface between R and the C++ estimate_local_complexity
 * function. It handles R/C++ data conversion and ensures proper memory management
 * using R's protection mechanism.
 *
 * @param x_r SEXP containing numeric vector of predictor values
 * @param y_r SEXP containing numeric vector of response values
 * @param center_idx_r SEXP containing integer scalar index (1-based) of center point
 * @param pilot_bandwidth_r SEXP containing numeric scalar bandwidth parameter
 *
 * @return SEXP containing a single numeric value representing the estimated complexity
 *
 * @note This function converts R's 1-based indexing to C++'s 0-based indexing
 *
 * @throws R Rf_error if:
 *         - Input types are incorrect (non-numeric x or y, non-integer center_idx)
 *         - Input vectors have unequal length
 *         - Any C++ exception occurs during computation
 */
SEXP S_estimate_local_complexity(
    SEXP x_r,
    SEXP y_r,
    SEXP center_idx_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r) {

    // Check input types
    if (!Rf_isReal(x_r) || !Rf_isReal(y_r) || !Rf_isInteger(center_idx_r) || !Rf_isReal(pilot_bandwidth_r)) {
        Rf_error("Invalid input types");
    }

    // Get input lengths and check consistency
    int n_points = LENGTH(x_r);
    if (n_points != LENGTH(y_r)) {
        Rf_error("Length of x and y must match");
    }

    // Create vectors from R data
    std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
    std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);

    // Get scalar inputs
    int center_idx = INTEGER(center_idx_r)[0] - 1;  // Convert from R 1-based to C++ 0-based
    double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
    int kernel_type = INTEGER(kernel_type_r)[0];

    // Allocate output
    SEXP result_r = PROTECT(Rf_allocVector(REALSXP, 1));

    // Compute result
    try {
        REAL(result_r)[0] = estimate_local_complexity(x, y, center_idx, pilot_bandwidth, kernel_type);
    } catch (const std::exception& e) {
        UNPROTECT(1);
        Rf_error("Error in estimate_local_complexity: %s", e.what());
    }

    UNPROTECT(1);
    return result_r;
}

/**
 * @brief Estimates local complexity by comparing linear and quadratic logistic fits
 *
 * @details Fits both linear and quadratic logistic models in a local neighborhood
 * and uses the difference between their predictions as a measure of local complexity.
 * Both models are fitted using univariate logistic regression with ridge regularization.
 *
 * @param x Vector of predictor values
 * @param y Vector of binary response values (should contain only 0s and 1s)
 * @param center_idx Index of the center point for local estimation
 * @param pilot_bandwidth Bandwidth for kernel weights (if ≤ 0, computed automatically)
 * @param kernel_type Integer specifying the kernel function type
 *
 * @return Estimated local complexity (weighted average of prediction differences,
 *         scaled by squared bandwidth)
 *
 * @throws std::invalid_argument If input vectors are empty or of unequal length,
 *                              or if center_idx is out of range
 *
 * @note Returns 0.0 if insufficient points for estimation
 * @note The complexity measure is normalized by the total weight and bandwidth
 */
double estimate_binary_local_complexity_comparative(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth,
    int kernel_type) {

    if (x.empty() || y.empty() || x.size() != y.size() ||
        center_idx < 0 || center_idx >= static_cast<int>(x.size())) {
        Rf_error("Invalid input parameters");
    }

    double center = x[center_idx];
    int n = x.size();

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = kernel_specific_bandwidth(x, kernel_type);
    }

    initialize_kernel(kernel_type, 1.0);

    // Collect local points and weights
    std::vector<double> local_x, local_y, local_w;
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {
            local_x.push_back((x[i] - center) / pilot_bandwidth);
            local_y.push_back(y[i]);

            double weight;
            kernel_fn(&dist, 1, &weight);
            local_w.push_back(weight);
        }
    }

    if (local_x.size() < 6) return 0.0;

    // Fit local linear model
    int max_iterations  = 100;
    double ridge_lambda = 0.0;
    double max_beta     = 100.0;
    double tolerance    = 1e-8;
    std::vector<double> pred1 = ulogit_predict(local_x.data(),
                                               local_y.data(),
                                               local_w,
                                               max_iterations,
                                               ridge_lambda,
                                               max_beta,
                                               tolerance);

    // Fit local quadratic model by adding squared terms
    std::vector<double> local_x2 = local_x;
    for (auto& val : local_x2) val = val * val;
    std::vector<double> pred2 = ulogit_predict(local_x2.data(),
                                               local_y.data(),
                                               local_w,
                                               max_iterations,
                                               ridge_lambda,
                                               max_beta,
                                               tolerance);

    // Compare models to estimate complexity
    double complexity = 0.0;
    double total_weight = 0.0;
    for (size_t i = 0; i < local_w.size(); ++i) {
        double diff = std::abs(pred2[i] - pred1[i]);
        complexity += local_w[i] * diff;
        total_weight += local_w[i];
    }

    return complexity / (total_weight * pilot_bandwidth * pilot_bandwidth);
}

/**
 * @brief Estimates local complexity using quadratic logistic regression
 *
 * @details Fits a quadratic logistic regression model in a local neighborhood and
 * uses the magnitude of the quadratic coefficient as a measure of local complexity.
 * The model is fitted using Newton-Raphson optimization with ridge regularization.
 *
 * The fitted model has the form: logit(p) = β₀ + β₁(x - x̄) + β₂(x - x̄)²
 *
 * @param x Vector of predictor values
 * @param y Vector of binary response values (should contain only 0s and 1s)
 * @param center_idx Index of the center point for local estimation
 * @param pilot_bandwidth Bandwidth for kernel weights (if ≤ 0, computed automatically)
 * @param kernel_type Integer specifying the kernel function type
 *
 * @return Estimated local complexity (absolute value of scaled quadratic coefficient)
 *
 * @throws std::invalid_argument If input vectors are empty or of unequal length,
 *                              or if center_idx is out of range
 *
 * @note Returns 0.0 if insufficient points for estimation
 * @note The quadratic coefficient is scaled by the square of the bandwidth
 */
#if 0
// init version
double estimate_binary_local_complexity_quadratic(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth,
    int kernel_type) {

    if (x.empty() || y.empty() || x.size() != y.size() ||
        center_idx < 0 || center_idx >= static_cast<int>(x.size())) {
        Rf_error("Invalid input parameters");
    }

    double center = x[center_idx];
    int n = x.size();

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = kernel_specific_bandwidth(x, kernel_type);
    }

    initialize_kernel(kernel_type, 1.0);

    // Pre-allocate vectors
    std::vector<double> local_x, local_y, local_w;
    local_x.reserve(n);
    local_y.reserve(n);
    local_w.reserve(n);

    // Collect local points and compute weights
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {
            local_x.push_back((x[i] - center) / pilot_bandwidth);  // Scale x values
            local_y.push_back(y[i]);

            double weight;
            kernel_fn(&dist, 1, &weight);
            local_w.push_back(weight);
        }
    }

    // Need more points for stable estimation of quadratic logistic regression
    if (local_x.size() < 6) {
        return 0.0;
    }

    // Initialize parameters for quadratic logistic regression
    Eigen::Vector3d beta = Eigen::Vector3d::Zero();
    const int max_iter = 100;
    const double tolerance = 1e-8;
    bool converged = false;

    // Newton-Raphson iteration for logistic regression
    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::Vector3d score = Eigen::Vector3d::Zero();
        Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();

        // Compute score and hessian
        for (size_t i = 0; i < local_x.size(); ++i) {
            double x_i = local_x[i];
            double x_i2 = x_i * x_i;

            // Design vector for current point
            Eigen::Vector3d z(1.0, x_i, x_i2);

            // Compute logistic probability
            double eta = beta.dot(z);
            double p = 1.0 / (1.0 + std::exp(-eta));

            // Weight for this observation
            double w = local_w[i];

            // Update score and hessian
            score += w * (local_y[i] - p) * z;
            hessian -= w * p * (1.0 - p) * (z * z.transpose());
        }

        // Solve for update
        Eigen::Vector3d delta = hessian.ldlt().solve(score);
        beta -= delta;

        // Check convergence
        if (delta.norm() < tolerance) {
            converged = true;
            break;
        }
    }

    return std::abs(beta(2)) / (pilot_bandwidth * pilot_bandwidth);
}
#endif


// I'll review the quadratic logistic regression implementation and check for correctness. Let me analyze it step by step:

// 1. The function's mathematical formulation and setup look correct:
// - Using logit(p) = β₀ + β₁(x - x̄) + β₂(x - x̄)²
// - Appropriate input validation
// - Proper kernel weight calculation
// - Good bandwidth handling

// 2. The data preparation is mostly correct:
// - Local points are collected within bandwidth
// - X values are properly centered and scaled by bandwidth
// - Weights are computed using the kernel function

// 3. In the quadratic logistic regression fitting, I notice a few potential issues:

// Key improvements and corrections made:

// 1. Added weight normalization:
//    - Weights now sum to 1.0 to ensure proper scaling
//    - Prevents numerical issues with very large or small weights

// 2. Added numerical stability improvements:
//    - Checks for extreme values in logistic calculation
//    - Prevents overflow/underflow in probability computation

// 3. Added ridge regularization:
//    - Helps with potential collinearity issues
//    - Makes the optimization more stable
//    - Added small constant to diagonal of Hessian

// 4. Improved linear system solver:
//    - Added fallback method if LDLT decomposition fails
//    - Uses complete orthogonal decomposition as backup
//    - More robust to ill-conditioned matrices

// 5. Maintained proper scaling:
//    - Correctly handles the bandwidth scaling for the quadratic coefficient
//    - Returns properly scaled complexity measure


double estimate_binary_local_complexity_quadratic(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth,
    int kernel_type) {

    if (x.empty() || y.empty() || x.size() != y.size() ||
        center_idx < 0 || center_idx >= static_cast<int>(x.size())) {
        Rf_error("Invalid input parameters");
    }

    double center = x[center_idx];
    int n = x.size();

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = kernel_specific_bandwidth(x, kernel_type);
    }

    initialize_kernel(kernel_type, 1.0);

    // Pre-allocate vectors
    std::vector<double> local_x, local_y, local_w;
    local_x.reserve(n);
    local_y.reserve(n);
    local_w.reserve(n);

    // Collect local points and compute weights
    double sum_weights = 0.0;
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {
            local_x.push_back((x[i] - center) / pilot_bandwidth);  // Scale x values
            local_y.push_back(y[i]);

            double weight;
            kernel_fn(&dist, 1, &weight);
            local_w.push_back(weight);
            sum_weights += weight;
        }
    }

    // Need more points for stable estimation of quadratic logistic regression
    if (local_x.size() < 6) {
        return 0.0;
    }

    // Initialize parameters for quadratic logistic regression
    Eigen::Vector3d beta = Eigen::Vector3d::Zero();
    const int max_iter = 100;
    const double tolerance = 1e-8;
    const double ridge_lambda = 1e-4;  // Ridge regularization parameter
    // Newton-Raphson iteration for logistic regression
    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::Vector3d score = Eigen::Vector3d::Zero();
        Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();

        // Add ridge regularization to diagonal of hessian
        hessian.diagonal() += ridge_lambda * Eigen::Vector3d::Ones();

        // Compute score and hessian
        for (size_t i = 0; i < local_x.size(); ++i) {
            double x_i = local_x[i];
            double x_i2 = x_i * x_i;

            // Design vector for current point
            Eigen::Vector3d z(1.0, x_i, x_i2);

            // Compute logistic probability with numerical stability
            double eta = beta.dot(z);
            double p;
            if (eta > 15.0) {  // Prevent overflow
                p = 1.0;
            } else if (eta < -15.0) {
                p = 0.0;
            } else {
                p = 1.0 / (1.0 + std::exp(-eta));
            }

            // Weight for this observation (normalized)
            double w = local_w[i] / sum_weights;

            // Update score and hessian
            score += w * (local_y[i] - p) * z;
            hessian -= w * p * (1.0 - p) * (z * z.transpose());
        }

        // Solve for update with regularized system
        Eigen::Vector3d delta;
        bool solved = false;

        // Try LDLT first
        Eigen::LDLT<Eigen::Matrix3d> ldlt(hessian);
        if (ldlt.info() == Eigen::Success) {
            delta = ldlt.solve(score);
            solved = true;
        }

        // Fall back to more stable but slower method if needed
        if (!solved) {
            delta = hessian.completeOrthogonalDecomposition().solve(score);
        }

        beta -= delta;

        // Check convergence
        if (delta.norm() < tolerance) {
            break;
        }
    }

    // Scale the quadratic coefficient back by bandwidth squared
    return std::abs(beta(2)) / (pilot_bandwidth * pilot_bandwidth);
}


//
// with model avaraging - over all points
//
/**
 * @brief Estimates local complexity using model-averaged quadratic logistic regression
 *
 * @details Computes local complexity estimates for every point in the input data
 * using model averaging. For each point, a quadratic logistic regression model is
 * fitted in its local neighborhood, and the magnitude of the quadratic coefficient
 * is used as a measure of local complexity. The final complexity estimate at each
 * point is a weighted average of the estimates from all models where that point
 * appears in the local neighborhood.
 *
 * The fitted models have the form: logit(p) = β₀ + β₁(x - x₀)/h + β₂((x - x₀)/h)²
 *
 * @param x Vector of predictor values
 * @param y Vector of binary response values (should contain only 0s and 1s)
 * @param pilot_bandwidth Bandwidth for kernel weights (if ≤ 0, computed automatically)
 * @param kernel_type Integer specifying the kernel function type
 *
 * @return Vector of estimated local complexity values, one for each input point
 *
 * @throws std::invalid_argument If input vectors are empty or of unequal length
 * @throws std::runtime_error If no successful model fits are achieved
 *
 * @note Returns 0.0 for points with insufficient neighbors for estimation
 * @note The quadratic coefficients are scaled by the square of the bandwidth
 * @note Model averaging helps reduce the variance of the complexity estimates
 * @note A minimum of 6 points in each local neighborhood is required for stable estimation
 */
std::vector<double> estimate_ma_binary_local_complexity_quadratic(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double pilot_bandwidth,
    int kernel_type) {

    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::invalid_argument("Invalid input parameters");
    }

    // Automatic bandwidth selection if needed
    if (pilot_bandwidth <= 0) {
        pilot_bandwidth = kernel_specific_bandwidth(x, kernel_type);
    }

    initialize_kernel(kernel_type, 1.0);

    int n_pts = x.size();
    double pilot_bandwidth2 = pilot_bandwidth * pilot_bandwidth;

    // Pre-allocate vectors with proper clearing mechanism
    std::vector<double> local_x, local_y, local_w;
    local_x.reserve(n_pts);
    local_y.reserve(n_pts);
    local_w.reserve(n_pts);

    // Storage for complexity estimates and their weights
    std::vector<std::vector<std::pair<double,double>>> pt_nbhds(n_pts);

    // Temporary storage for convergence checking
    bool any_successful_fit = false;

    for (int pti = 0; pti < n_pts; pti++) {
        double center = x[pti];

        // Clear vectors for reuse
        local_x.clear();
        local_y.clear();
        local_w.clear();

        // Collect local points and compute weights
        double sum_weights = 0.0;
        for (int i = 0; i < n_pts; ++i) {
            double dist = std::abs(x[i] - center) / pilot_bandwidth;
            if (dist < 1.0) {
                local_x.push_back(x[i] - center);
                local_y.push_back(y[i]);

                double weight;
                kernel_fn(&dist, 1, &weight);
                local_w.push_back(weight);
                sum_weights += weight;
            }
        }

        // Need more points for stable estimation
        if (local_x.size() < 6) {
            continue;
        }

        // Normalize weights once
        for (auto& w : local_w) {
            w /= sum_weights;
        }

        // Store the normalized weight for the center point
        double pt_weight = 0.0;
        for (size_t i = 0; i < local_x.size(); ++i) {
            if (std::abs(local_x[i]) < 1e-10) {  // Check if this is the center point
                pt_weight = local_w[i];
                break;
            }
        }

        // Initialize parameters
        Eigen::Vector3d beta = Eigen::Vector3d::Zero();
        const int max_iter = 100;
        const double tolerance = 1e-8;
        const double ridge_lambda = 1e-4;
        bool converged = false;

        // Newton-Raphson iteration
        for (int iter = 0; iter < max_iter && !converged; ++iter) {
            Eigen::Vector3d score = Eigen::Vector3d::Zero();
            Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();

            // Add ridge regularization
            hessian.diagonal() += ridge_lambda * Eigen::Vector3d::Ones();

            // Compute score and hessian
            for (size_t i = 0; i < local_x.size(); ++i) {
                double x_i = local_x[i];
                double x_i2 = x_i * x_i;
                Eigen::Vector3d z(1.0, x_i, x_i2);

                // Compute logistic probability
                double eta = beta.dot(z);
                double p = (eta > 15.0) ? 1.0 :
                          (eta < -15.0) ? 0.0 :
                          1.0 / (1.0 + std::exp(-eta));

                score += local_w[i] * (local_y[i] - p) * z;
                hessian -= local_w[i] * p * (1.0 - p) * (z * z.transpose());
            }

            // Solve system with fallback
            Eigen::Vector3d delta;
            Eigen::LDLT<Eigen::Matrix3d> ldlt(hessian);

            if (ldlt.info() == Eigen::Success) {
                delta = ldlt.solve(score);
            } else {
                delta = hessian.completeOrthogonalDecomposition().solve(score);
            }

            // Update parameters
            beta -= delta;

            // Check convergence
            if (delta.norm() < tolerance) {
                converged = true;
                any_successful_fit = true;

                // Store result only if converged
                double complexity = std::abs(beta(2)) / pilot_bandwidth2;
                pt_nbhds[pti].push_back(std::make_pair(pt_weight, complexity));
            }
        }
    }

    // Check if we had any successful fits
    if (!any_successful_fit) {
        Rf_error("No successful model fits achieved");
    }

    // Compute weighted averages
    std::vector<double> cx(n_pts, 0.0);
    for (int pti = 0; pti < n_pts; pti++) {
        const auto& v = pt_nbhds[pti];
        if (v.empty()) {
            cx[pti] = 0.0;  // No valid estimates for this point
            continue;
        }

        double total_weight = 0.0;
        for (const auto& p : v) {
            cx[pti] += p.first * p.second;
            total_weight += p.first;
        }

        if (total_weight > 0.0) {
            cx[pti] /= total_weight;
        }
    }

    return cx;
}

/**
 * @brief R interface to model-averaged binary local complexity estimation
 *
 * @details This function provides an R interface to the C++ implementation of
 * model-averaged binary local complexity estimation using quadratic logistic regression.
 * It handles conversion between R and C++ data types and provides R-appropriate
 * Rf_error handling.
 *
 * @param x_r SEXP (numeric vector) Predictor values
 * @param y_r SEXP (numeric vector) Binary response values (should contain only 0s and 1s)
 * @param pilot_bandwidth_r SEXP (numeric scalar) Bandwidth for kernel weights
 * @param kernel_type_r SEXP (integer scalar) Integer specifying the kernel function type
 *
 * @return SEXP (numeric vector) Estimated local complexity values for each point
 *
 * @throws R Rf_error if inputs are invalid or if computation fails
 *
 * @note This function is not intended to be called directly from C/C++
 * @note R_NilValue inputs will trigger an Rf_error
 * @note Non-matching vector lengths will trigger an Rf_error
 * @note Non-numeric or non-integer inputs will trigger an Rf_error
 */
SEXP S_estimate_ma_binary_local_complexity_quadratic(SEXP x_r,
                                                     SEXP y_r,
                                                     SEXP pilot_bandwidth_r,
                                                     SEXP kernel_type_r) {
    // Check for NULL inputs
    if (x_r == R_NilValue || y_r == R_NilValue ||
        pilot_bandwidth_r == R_NilValue || kernel_type_r == R_NilValue) {
        Rf_error("Input arguments cannot be NULL");
    }

    // Type checking with more specific Rf_error messages
    if (!Rf_isReal(x_r)) Rf_error("'x' must be a numeric vector");
    if (!Rf_isReal(y_r)) Rf_error("'y' must be a numeric vector");
    if (!Rf_isReal(pilot_bandwidth_r)) Rf_error("'pilot_bandwidth' must be numeric");
    if (!Rf_isInteger(kernel_type_r)) Rf_error("'kernel_type' must be integer");

    // Length checking
    int n_points = LENGTH(x_r);
    if (n_points < 1) Rf_error("Input vectors cannot be empty");
    if (n_points != LENGTH(y_r)) Rf_error("Length of x and y must match");
    if (LENGTH(pilot_bandwidth_r) != 1) Rf_error("pilot_bandwidth must be a scalar");
    if (LENGTH(kernel_type_r) != 1) Rf_error("kernel_type must be a scalar");

    // Convert inputs to C++ vectors
    std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
    std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);
    double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
    int kernel_type = INTEGER(kernel_type_r)[0];

    // Allocate result vector
    SEXP result_r = PROTECT(Rf_allocVector(REALSXP, n_points));

    try {
        auto result = estimate_ma_binary_local_complexity_quadratic(
            x, y, pilot_bandwidth, kernel_type);

        // Copy results back to R vector
        std::copy(result.begin(), result.end(), REAL(result_r));

    } catch (const std::exception& e) {
        UNPROTECT(1);
        Rf_error("Error in estimate_ma_binary_local_complexity: %s", e.what());
    }

    UNPROTECT(1);
    return result_r;
}


/**
 * @brief R interface for local complexity estimation of binary response
 *
 * @details Provides R interface to both quadratic and comparative methods of
 * estimating local complexity. Handles R-specific memory management and
 * type conversions.
 *
 * @param x_r SEXP containing numeric vector of predictor values
 * @param y_r SEXP containing numeric vector of binary responses
 * @param center_idx_r SEXP containing integer center point index (1-based)
 * @param pilot_bandwidth_r SEXP containing numeric bandwidth
 * @param kernel_type_r SEXP containing integer kernel type
 * @param method_r SEXP containing integer method selector (1=quadratic, 2=comparative)
 *
 * @return SEXP containing single numeric value of estimated complexity
 *
 * @note Converts R's 1-based indices to C++'s 0-based indices
 * @note Handles R's Rf_error reporting via Rf_error
 */
SEXP S_estimate_binary_local_complexity(
    SEXP x_r,
    SEXP y_r,
    SEXP center_idx_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP method_r) {  // 1 for quadratic, 2 for comparative

    if (!Rf_isReal(x_r) || !Rf_isReal(y_r) || !Rf_isInteger(center_idx_r) || !Rf_isReal(pilot_bandwidth_r)) {
        Rf_error("Invalid input types");
    }

    int n_points = LENGTH(x_r);
    if (n_points != LENGTH(y_r)) {
        Rf_error("Length of x and y must match");
    }

    std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
    std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);

    int center_idx = INTEGER(center_idx_r)[0] - 1;  // Convert from R 1-based to C++ 0-based
    double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
    int kernel_type = INTEGER(kernel_type_r)[0];
    int method = INTEGER(method_r)[0];

    SEXP result_r = PROTECT(Rf_allocVector(REALSXP, 1));

    try {
        double result;
        if (method == 1) {
            result = estimate_binary_local_complexity_quadratic(x, y, center_idx,
                                                             pilot_bandwidth, kernel_type);
        } else if (method == 2) {
            result = estimate_binary_local_complexity_comparative(x, y, center_idx,
                                                               pilot_bandwidth, kernel_type);
        } else {
            Rf_error("Invalid method: must be 1 (quadratic) or 2 (comparative)");
        }
        REAL(result_r)[0] = result;
    } catch (const std::exception& e) {
        UNPROTECT(1);
        Rf_error("Error in estimate_binary_local_complexity: %s", e.what());
    }

    UNPROTECT(1);
    return result_r;
}

#if 0
SEXP S_estimate_binary_local_complexity(
    SEXP x_r,
    SEXP y_r,
    SEXP center_idx_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r) {

    int n_protected = 0;

    if (!Rf_isReal(x_r) || !Rf_isReal(y_r) || !Rf_isInteger(center_idx_r) || !Rf_isReal(pilot_bandwidth_r)) {
        Rf_error("Invalid input types");
    }

    int n_points = LENGTH(x_r);
    if (n_points != LENGTH(y_r)) {
        Rf_error("Length of x and y must match");
    }

    std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
    std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);

    int center_idx = INTEGER(center_idx_r)[0] - 1;  // Convert from R 1-based to C++ 0-based
    double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
    int kernel_type = INTEGER(kernel_type_r)[0];

    SEXP result_r = PROTECT(Rf_allocVector(REALSXP, 1));

    try {
        REAL(result_r)[0] = estimate_binary_local_complexity(x, y, center_idx, pilot_bandwidth, kernel_type);
    } catch (const std::exception& e) {
        UNPROTECT(1);
        Rf_error("Error in estimate_binary_local_complexity: %s", e.what());
    }

    UNPROTECT(1);
    return result_r;
}
#endif
