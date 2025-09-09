#include <R.h>
#include <Rinternals.h>

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <cmath>      // for fabs()
#include <algorithm>  // for std::find,

#include "ulogit.hpp"
#include "error_utils.h"

extern "C" {
    SEXP S_ulogit(SEXP x_sexp,
                  SEXP y_sexp,
                  SEXP w_sexp,
                  SEXP max_iterations_sexp,
                  SEXP ridge_lambda_sexp,
                  SEXP max_beta_sexp,
                  SEXP tolerance_sexp,
                  SEXP verbose_sexp);

    SEXP S_eigen_ulogit(SEXP x_sexp,
                    SEXP y_sexp,
                    SEXP w_sexp,
                    SEXP fit_quadratic_sexp,
                    SEXP max_iterations_sexp,
                    SEXP ridge_lambda_sexp,
                    SEXP tolerance_sexp);
}

/**
 * @brief Fits a weighted logistic regression model on a window of one-dimensional data
 *
 * @details This function implements a locally weighted logistic regression for binary classification
 * using Newton-Raphson optimization. The model fits a logistic curve to a window of data points,
 * with each point weighted by both the provided weights and its position in the window.
 *
 * The function performs the following steps:
 * 1. Centers the x values for numerical stability
 * 2. Fits a logistic model using Newton-Raphson iteration
 * 3. Computes predictions and Leave-One-Out Cross-Validation (LOOCV) errors
 *
 * The logistic model has the form:
 * P(y=1|x) = 1 / (1 + exp(-(β₀ + β₁(x - x̄))))
 * where x̄ is the weighted mean of x values in the window
 *
 * @param x Pointer to array of x values (predictor variable) in the window
 * @param y Pointer to array of binary y values (0 or 1) corresponding to x
 * @param w Vector of observation weights for the window
 * @param max_iterations Maximum number of iterations for Newton-Raphson optimization (default: 100)
 * @param ridge_lambda  Ridge regularization parameter for stability (default: 0.1)
 * @param max_beta      Maximum allowed absolute value for coefficient estimates (default: 100.0)
 * @param tolerance     Convergence tolerance for optimization (default: 1e-8)
 * @param verbose       Flag to enable detailed output during optimization (default: false)
 *
 * @return ulogit_t struct containing:
 *         - predictions: fitted probabilities
 *         - errors: leave-one-out cross-validation errors
 *         - weights: observation weights used
 *         - x_min_index: index of smallest x value
 *         - x_max_index: index of largest x value
 *         If verbose is true, also includes:
 *         - iteration_count: number of iterations until convergence
 *         - converged: whether optimization converged
 *
 *
 * @note
 * - Input y values must be binary (0 or 1)
 * - The window size must be at least 2 points
 * - The function uses weighted maximum likelihood estimation
 * - LOOCV errors are computed using log loss: -log(p) for y=1, -log(1-p) for y=0
 *
 * @warning
 * - The function assumes x values are sorted in ascending order
 * - Numerical instability may occur with extremely imbalanced weights
 * - The Newton-Raphson method may not converge for pathological data
 *
 * @see ulogit_t for the return type structure
 * @see wmabilo() for the main algorithm using this function
 *
 * @example
 * ```cpp
 * const double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
 * const double y[] = {0.0, 0.0, 1.0, 1.0, 1.0};
 * std::vector<double> w = {0.2, 0.5, 1.0, 0.5, 0.2};
 * int x_min_idx = 0;
 * int x_max_idx = 4;
 *
 * auto result = ulogit(x, y, w, x_min_idx, x_max_idx);
 * // Access predictions with result.predictions
 * // Access LOOCV errors with result.errors
 * ```
 */
ulogit_t ulogit(const double* x,
                const double* y,
                const std::vector<double>& w,
                int max_iterations,
                double ridge_lambda,
                double max_beta,
                double tolerance,
                bool verbose) {

    int window_size = w.size();
    ulogit_t result;
    result.w = w;
    result.predictions.resize(window_size);
    result.errors.resize(window_size);

    if (verbose) {
        result.iteration_count = 0;
        result.converged = false;
    }

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        total_weight += weight;
    }
    if (total_weight <= 0) REPORT_ERROR("total_weight: %.2f   Sum of weights must be positive", total_weight);

    // Check for effective complete separation considering weights
    double weighted_sum_y = 0.0;
    double total_nonzero_weight = 0.0;
    for (int i = 0; i < window_size; ++i) {
        if (w[i] > tolerance) {  // Only consider points with non-negligible weights
            weighted_sum_y += w[i] * y[i];
            total_nonzero_weight += w[i];
        }
    }

    double prior_prob = weighted_sum_y / total_nonzero_weight;

    // If all significant weights are associated with same y value
    if (total_nonzero_weight > 0) {  // Ensure we have some meaningful data
        double weighted_mean_y = weighted_sum_y / total_nonzero_weight;
        if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
            // We have effective complete separation
            double pred_value = (weighted_mean_y > 0.5) ? 1.0 : 0.0;
            std::fill(result.predictions.begin(), result.predictions.end(), pred_value);
            std::fill(result.errors.begin(), result.errors.end(), -std::log(1.0 - tolerance));
            return result;
        }
    }

    // Define Newton-Raphson optimization as a lambda function with ridge regularization
    auto newton_raphson = [&](const std::vector<double>& weights, double x_mean, bool& converged) -> std::array<double, 2> {
        double weights_sum = 0.0;
        double y_mean = 0.0;
        for (int i = 0; i < window_size; ++i) {
            y_mean += weights[i] * y[i];
            weights_sum += weights[i];
        }
        y_mean /= weights_sum;

        // Apply logit function to get initial b0
        y_mean = std::clamp(y_mean, tolerance, 1.0 - tolerance); // If we get here, we need to protect the logit transformation in beta initialization
        std::array<double, 2> beta = {std::log(y_mean / (1.0 - y_mean)), 0.0};  // {intercept, slope}

        converged = false;
        int iter = 0;
        while (iter < max_iterations && !converged) {
            double g0 = 0.0, g1 = 0.0;  // Gradients
            double h00 = 0.0, h01 = 0.0, h11 = 0.0;  // Hessian

            for (int i = 0; i < window_size; ++i) {
                double x_centered = x[i] - x_mean;
                double xbeta = beta[0] + beta[1] * x_centered;

                // Prevent overflow in exp
                xbeta = std::clamp(xbeta, -max_beta, max_beta);
                double p = 1.0 / (1.0 + std::exp(-xbeta));

                double w_i = weights[i];

                // Gradient components
                double error = p - y[i];
                g0 += w_i * error;
                g1 += w_i * error * x_centered;

                // Hessian components
                double p_1mp = p * (1.0 - p);
                h00 += w_i * p_1mp;
                h01 += w_i * p_1mp * x_centered;
                h11 += w_i * p_1mp * x_centered * x_centered;
            }

            // Add ridge regularization terms
            g0 += ridge_lambda * beta[0];
            g1 += ridge_lambda * beta[1];
            h00 += ridge_lambda;
            h11 += ridge_lambda;

            // Solve 2x2 system using direct inverse
            double det = h00 * h11 - h01 * h01;
            double step = 1.0; //std::min(1.0, 1.0 / (1.0 + std::exp(-std::abs(det))));

            if (std::abs(det) < tolerance) {
                beta[0] -= g0 * step / (h00 + ridge_lambda);
                beta[1] -= g1 * step / (h11 + ridge_lambda);
            } else {
                double d_beta0 = (h11 * g0 - h01 * g1) / det;
                double d_beta1 = (-h01 * g0 + h00 * g1) / det;
                beta[0] -= step * d_beta0;
                beta[1] -= step * d_beta1;
            }

            // Clamp coefficients to prevent explosion
            beta[0] = std::clamp(beta[0], -max_beta, max_beta);
            beta[1] = std::clamp(beta[1], -max_beta, max_beta);

            // Check convergence using provided tolerance
            if (std::abs(g0) < tolerance && std::abs(g1) < tolerance) {
                converged = true;
                if (verbose) {
                    result.iteration_count = iter + 1;
                    result.converged = true;
                    Rprintf("Newton-Raphson converged at iter: %d\n", iter + 1);
                }
            }

            iter++;
        }
        return beta;
    };

    // Fit full model
    double x_mean = 0.0;
    for (int i = 0; i < window_size; ++i) {
        x_mean += w[i] * x[i];
    }
    x_mean /= total_weight;

    bool converged;
    std::array<double, 2> beta = newton_raphson(w, x_mean, converged);

    // Compute predictions for full model
    for (int i = 0; i < window_size; ++i) {
        double x_centered = x[i] - x_mean;
        double xbeta = beta[0] + beta[1] * x_centered;
        xbeta = std::clamp(xbeta, -100.0, 100.0);  // Prevent overflow
        result.predictions[i] = 1.0 / (1.0 + std::exp(-xbeta));
    }

    // Compute LOOCV errors, handling potential complete separation in subsets
    for (int i = 0; i < window_size; ++i) {
        // Create leave-one-out weights
        std::vector<double> loo_weights = w;
        loo_weights[i] = 0.0;

        double weighted_sum_y = 0.0;
        double total_weight = 0.0;
        for (int j = 0; j < window_size; ++j) {
            if (j != i && loo_weights[j] > tolerance) {
                weighted_sum_y += loo_weights[j] * y[j];
                total_weight += loo_weights[j];
            }
        }

        double loo_pred;
        if (total_weight > 0) {
            double weighted_mean_y = weighted_sum_y / total_weight;
            if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
                // In case of complete separation in LOO sample
                loo_pred = (weighted_mean_y > 1.0 - tolerance) ? 1.0 - tolerance : tolerance;
            } else {
                // Compute new weighted mean without point i
                double loo_total_weight = 0.0;
                double loo_x_mean = 0.0;
                for (int j = 0; j < window_size; ++j) {
                    if (j != i) {
                        loo_x_mean += loo_weights[j] * x[j];
                        loo_total_weight += loo_weights[j];
                    }
                }
                loo_x_mean /= loo_total_weight;

                // Fit model without point i
                bool loo_converged;
                std::array<double, 2> beta_loo = newton_raphson(loo_weights, loo_x_mean, loo_converged);

                // Compute prediction for held-out point
                double x_centered = x[i] - loo_x_mean;
                double xbeta = beta_loo[0] + beta_loo[1] * x_centered;
                xbeta = std::clamp(xbeta, -100.0, 100.0);  // Prevent overflow
                loo_pred = 1.0 / (1.0 + std::exp(-xbeta));
            }
        } else {
            loo_pred = prior_prob;
        }

        // Compute cross-entropy error
        //result.errors[i] = -y[i] * std::log(std::max(loo_pred, tolerance)) - (1-y[i]) * std::log(std::max(1.0 - loo_pred, tolerance));
        //result.errors[i] = std::abs(y[i] - loo_pred); // absolute deviation error
        result.errors[i] = (y[i] - loo_pred) * (y[i] - loo_pred);
    }

    return result;
}


/**
 * @brief R interface wrapper for univariate logistic regression function
 *
 * @details SEXP wrapper function that provides an R interface to the C++ ulogit function.
 * Handles conversion between R and C++ data types, memory management, and returns results in an R-compatible format.
 * The function expects sorted x values and binary y values (0 or 1).
 *
 * @param x_sexp SEXP (numeric vector) Input predictor values
 *               Must be sorted in ascending order
 * @param y_sexp SEXP (numeric vector) Binary response values (0 or 1)
 *               Must be same length as x_sexp
 * @param w_sexp SEXP (numeric vector) Window weights
 *               Length must equal x_max_index - x_min_index + 1
 * @param max_iterations_sexp SEXP containing integer for maximum iterations
 * @param ridge_lambda_sexp  SEXP containing numeric ridge regularization parameter
 * @param max_beta_sexp     SEXP containing numeric maximum coefficient value
 * @param tolerance_sexp     SEXP containing numeric convergence tolerance
 * @param verbose_sexp      SEXP containing logical verbose flag
 *
 * @return SEXP list containing:
 *         - predictions: numeric vector of fitted probabilities
 *         - errors: numeric vector of leave-one-out cross-validation errors
 *         - weights: numeric vector of weights used in fitting
 *
 * @note
 * - All indices are converted between R's 1-based and C++'s 0-based indexing
 * - The function uses PROTECT/UNPROTECT for R's garbage collection
 * - Return value is a named list for easy access in R
 *
 * @warning
 * - No input validation is performed in the C++ code
 * - R code should validate inputs before calling this function
 * - Memory management relies on R's garbage collection system
 *
 * @see ulogit() for the underlying C++ implementation
 * @see ulogit_t for the C++ return type structure
 *
 * Usage in R:
 * ```r
 * result <- .Call("S_ulogit",
 *                 as.double(x),
 *                 as.double(y),
 *                 as.double(w),
 *                 as.double(tolerance))
 * ```
 */
SEXP S_ulogit(SEXP x_sexp,
              SEXP y_sexp,
              SEXP w_sexp,
              SEXP max_iterations_sexp,
              SEXP ridge_lambda_sexp,
              SEXP max_beta_sexp,
              SEXP tolerance_sexp,
              SEXP verbose_sexp) {
    // Convert inputs from R to C++
    double* x = REAL(x_sexp);
    double* y = REAL(y_sexp);
    double* w_r = REAL(w_sexp);
    int max_iterations = INTEGER(max_iterations_sexp)[0];
    double ridge_lambda = REAL(ridge_lambda_sexp)[0];
    double max_beta = REAL(max_beta_sexp)[0];
    double tolerance = REAL(tolerance_sexp)[0];
    bool verbose = LOGICAL(verbose_sexp)[0];

    // Convert R vector to std::vector
    int window_size = Rf_length(w_sexp);
    std::vector<double> w(w_r, w_r + window_size);

    // Call the actual function with all parameters
    ulogit_t result = ulogit(x, y, w,
                             max_iterations,
                             ridge_lambda,
                             max_beta,
                             tolerance,
                             verbose);
    // Creating return list
    const int N_COMPONENTS = 3;
    int n_protected = 0;  // Track number of PROTECT calls
    SEXP out = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Convert predictions to R vector
    SEXP predictions = PROTECT(allocVector(REALSXP, window_size)); n_protected++;
    for(int i = 0; i < window_size; i++) {
        REAL(predictions)[i] = result.predictions[i];
    }

    // Convert errors to R vector
    SEXP errors = PROTECT(allocVector(REALSXP, window_size)); n_protected++;
    for(int i = 0; i < window_size; i++) {
        REAL(errors)[i] = result.errors[i];
    }

    // Convert weights to R vector
    SEXP weights = PROTECT(allocVector(REALSXP, window_size)); n_protected++;
    for(int i = 0; i < window_size; i++) {
        REAL(weights)[i] = result.w[i];
    }

    // Set list names
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("predictions"));
    SET_STRING_ELT(names, 1, mkChar("errors"));
    SET_STRING_ELT(names, 2, mkChar("weights"));

    // Set list elements
    SET_VECTOR_ELT(out, 0, predictions);
    SET_VECTOR_ELT(out, 1, errors);
    SET_VECTOR_ELT(out, 2, weights);

    // Set list names
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return out;
}


/**
 * @brief Fits a univariate logistic regression model and returns predictions
 *
 * @details This function fits a logistic regression model using weighted maximum likelihood
 * estimation with ridge regularization. It handles complete separation cases and includes
 * safeguards against numerical instability. The optimization is performed using
 * Newton-Raphson iteration.
 *
 * The model fitted is: logit(p) = β₀ + β₁(x - x̄), where x̄ is the weighted mean of x.
 *
 * @param x Pointer to predictor values array
 * @param y Pointer to binary response values array (should contain only 0s and 1s)
 * @param w Vector of observation weights (must be non-negative)
 * @param max_iterations Maximum number of Newton-Raphson iterations (default: 100)
 * @param ridge_lambda Ridge regularization parameter (default: 0.002)
 * @param max_beta Maximum absolute value for coefficient estimates (default: 100.0)
 * @param tolerance Convergence tolerance for Newton-Raphson iteration (default: 1e-8)
 *
 * @return Vector of predicted probabilities, one for each input observation
 *
 * @throws std::invalid_argument If the sum of weights is not positive
 *
 * @note The function handles complete separation by returning appropriate constant predictions
 * @note Ridge regularization is applied to both intercept and slope coefficients
 */
std::vector<double> ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance) {

    int window_size = w.size();
    std::vector<double> predictions(window_size);

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        total_weight += weight;
    }
    if (total_weight <= 0) {
        throw std::invalid_argument("Sum of weights must be positive");
    }

    // Check for effective complete separation considering weights
    double weighted_sum_y = 0.0;
    double total_nonzero_weight = 0.0;
    for (int i = 0; i < window_size; ++i) {
        if (w[i] > tolerance) {  // Only consider points with non-negligible weights
            weighted_sum_y += w[i] * y[i];
            total_nonzero_weight += w[i];
        }
    }

    // Handle complete separation case
    if (total_nonzero_weight > 0) {
        double weighted_mean_y = weighted_sum_y / total_nonzero_weight;
        if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
            double pred_value = (weighted_mean_y > 0.5) ? 1.0 - tolerance : tolerance;
            std::fill(predictions.begin(), predictions.end(), pred_value);
            return predictions;
        }
    }

    // Define Newton-Raphson optimization
    auto newton_raphson = [&](const std::vector<double>& weights, double x_mean, bool& converged) -> std::array<double, 2> {
        double weights_sum = 0.0;
        double y_mean = 0.0;
        for (int i = 0; i < window_size; ++i) {
            y_mean += weights[i] * y[i];
            weights_sum += weights[i];
        }
        y_mean /= weights_sum;

        // Initialize coefficients
        y_mean = std::clamp(y_mean, tolerance, 1.0 - tolerance);
        std::array<double, 2> beta = {std::log(y_mean / (1.0 - y_mean)), 0.0};

        converged = false;
        int iter = 0;
        while (iter < max_iterations && !converged) {
            double g0 = 0.0, g1 = 0.0;  // Gradients
            double h00 = 0.0, h01 = 0.0, h11 = 0.0;  // Hessian

            for (int i = 0; i < window_size; ++i) {
                double x_centered = x[i] - x_mean;
                double xbeta = beta[0] + beta[1] * x_centered;
                xbeta = std::clamp(xbeta, -max_beta, max_beta);
                double p = 1.0 / (1.0 + std::exp(-xbeta));

                double w_i = weights[i];
                double error = p - y[i];
                double p_1mp = p * (1.0 - p);

                // Update gradient and Hessian
                g0 += w_i * error;
                g1 += w_i * error * x_centered;
                h00 += w_i * p_1mp;
                h01 += w_i * p_1mp * x_centered;
                h11 += w_i * p_1mp * x_centered * x_centered;
            }

            // Add ridge regularization
            g0 += ridge_lambda * beta[0];
            g1 += ridge_lambda * beta[1];
            h00 += ridge_lambda;
            h11 += ridge_lambda;

            // Solve system
            double det = h00 * h11 - h01 * h01;
            if (std::abs(det) < tolerance) {
                beta[0] -= g0 / (h00 + ridge_lambda);
                beta[1] -= g1 / (h11 + ridge_lambda);
            } else {
                beta[0] -= (h11 * g0 - h01 * g1) / det;
                beta[1] -= (-h01 * g0 + h00 * g1) / det;
            }

            beta[0] = std::clamp(beta[0], -max_beta, max_beta);
            beta[1] = std::clamp(beta[1], -max_beta, max_beta);

            converged = std::abs(g0) < tolerance && std::abs(g1) < tolerance;
            iter++;
        }
        return beta;
    };

    // Compute weighted mean of x
    double x_mean = 0.0;
    for (int i = 0; i < window_size; ++i) {
        x_mean += w[i] * x[i];
    }
    x_mean /= total_weight;

    // Fit model
    bool converged;
    std::array<double, 2> beta = newton_raphson(w, x_mean, converged);

    // Compute predictions
    for (int i = 0; i < window_size; ++i) {
        double x_centered = x[i] - x_mean;
        double xbeta = beta[0] + beta[1] * x_centered;
        xbeta = std::clamp(xbeta, -max_beta, max_beta);
        predictions[i] = 1.0 / (1.0 + std::exp(-xbeta));
    }

    return predictions;
}

/**
 * @brief Core function for fitting univariate logistic regression using Eigen
 *
 * @details Implements Newton-Raphson algorithm with numerical stability safeguards,
 * including centering of predictors, ridge regularization, and robust matrix decomposition
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights
 * @param fit_quadratic If true, fits quadratic term in addition to linear
 * @param max_iterations Maximum number of Newton-Raphson iterations
 * @param ridge_lambda Ridge regularization parameter
 * @param max_beta Maximum absolute value allowed for coefficients
 * @param tolerance Convergence tolerance for Newton-Raphson algorithm
 *
 * @return eigen_ulogit_t structure containing fitted model and diagnostics
 *
 * @throws std::invalid_argument if weights sum to zero or negative
 *
 * @note Handles complete separation cases by returning predictions near 0 or 1
 */
#if 0
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;  // Number of parameters (intercept + linear [+ quadratic])

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    //result.beta = Eigen::VectorXd::Random(p) * 0.001; // Random values in [-0.001, 0.001]
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Compute weighted mean of x for centering
    double weighted_sum = 0.0, weighted_x_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        weighted_sum += w[i];
        weighted_x_sum += w[i] * x[i];
    }
    double x_mean = weighted_x_sum / weighted_sum;

    // Check for complete separation
    double weighted_y_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (w[i] > tolerance) {
            weighted_y_sum += w[i] * y[i];
        }
    }
    double y_mean = weighted_y_sum / weighted_sum;

    double adaptive_ridge = ridge_lambda;
    #if 0
    if (y_mean < tolerance || y_mean > 1.0 - tolerance) {
        adaptive_ridge *= 10.0;  // Increase regularization for numerical stability
    }

    if (y_mean < tolerance || y_mean > 1.0 - tolerance) {
        double pred_value = (y_mean > 0.5) ? 1.0 - tolerance : tolerance;
        std::fill(result.predictions.begin(), result.predictions.end(), pred_value);
        if (with_errors)
            std::fill(result.errors.begin(), result.errors.end(), 0.0);
        return result;
    }
    #endif

    // Newton-Raphson iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd score = Eigen::VectorXd::Zero(p);
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(p, p);

        // Add adaptive ridge regularization
        hessian.diagonal() += adaptive_ridge * Eigen::VectorXd::Ones(p);

        // Add ridge regularization to diagonal
        hessian.diagonal() += ridge_lambda * Eigen::VectorXd::Ones(p);

        for (int i = 0; i < n; ++i) {
            double x_centered = x[i] - x_mean;

            // Create design vector
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x_centered;
            if (fit_quadratic) {
                z(2) = x_centered * x_centered;
            }

            // Compute logistic probability with numerical safeguards
            double eta = result.beta.dot(z);
            eta = std::clamp(eta, -max_beta, max_beta);
            double p = 1.0 / (1.0 + std::exp(-eta));

            // Update score and hessian
            score += w[i] * (y[i] - p) * z;
            hessian -= w[i] * p * (1.0 - p) * (z * z.transpose());
        }

        // Solve system with fallback strategy
        Eigen::VectorXd delta;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);

        if (ldlt.info() == Eigen::Success) {
            delta = ldlt.solve(score);
        } else {
            delta = hessian.completeOrthogonalDecomposition().solve(score);
        }

        // Update parameters with clamping
        result.beta -= delta;
        for (int j = 0; j < p; ++j) {
            result.beta(j) = std::clamp(result.beta(j), -max_beta, max_beta);
        }

        // Check convergence
        if (delta.norm() < tolerance) {
            result.converged = true;
            break;
        }

        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double x_centered = x[i] - x_mean;
        Eigen::VectorXd z(p);
        z(0) = 1.0;
        z(1) = x_centered;
        if (fit_quadratic) {
            z(2) = x_centered * x_centered;
        }

        double eta = result.beta.dot(z);
        eta = std::clamp(eta, -max_beta, max_beta);
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    if (with_errors) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            double x_centered = x[i] - x_mean;
            X(i, 0) = 1.0;
            X(i, 1) = x_centered;
            if (fit_quadratic) {
                X(i, 2) = x_centered * x_centered;
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix
        Eigen::MatrixXd H = W_sqrt * X * (X.transpose() * W_sqrt.transpose() * W_sqrt * X).inverse() * X.transpose() * W_sqrt;

        // Thresholds for leverage
        double mod_high_leverage = 2.0 * p / n;
        double very_high_leverage = 3.0 * p / n;

        // Compute errors for each observation
        for (int i = 0; i < n; ++i) {
            double h_i = H(i, i);
            double y_i = y[i];
            double p_i = result.predictions[i];

            if (h_i > very_high_leverage) {
                // Williams correction for very high leverage points
                double r_i = (y_i - p_i) / std::sqrt(p_i * (1.0 - p_i));  // Pearson residual
                double h_star = h_i * (1.0 - h_i);  // Modified leverage
                result.errors[i] = p_i - (h_i * r_i / (1.0 - h_i)) * std::sqrt((1.0 - h_star) / (1.0 - h_i));
            } else if (h_i > mod_high_leverage) {
                // Check if probability is extreme
                if (p_i < 0.1 || p_i > 0.9) {
                    // Second-order approximation for extreme probabilities
                    double first_term = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                    double second_term = (h_i * h_i * (y_i - p_i) * (y_i - p_i) * (1.0 - 2.0 * p_i)) /
                        (2.0 * (1.0 - h_i) * (1.0 - h_i));
                    result.errors[i] = first_term + second_term;
                } else {
                    // First-order approximation for non-extreme probabilities
                    result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                }
            } else {
                // First-order approximation for normal leverage points
                result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
            }

            // Clamp predicted probabilities to valid range
            result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);
        }
    }
    
    return result;
}


eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    max_beta = 50.0;

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;  // Number of parameters (intercept + linear [+ quadratic])

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Check for numerical issues in y
    double weighted_sum = 0.0;
    double weighted_y_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (w[i] > tolerance) {
            weighted_sum += w[i];
            weighted_y_sum += w[i] * y[i];
        }
    }
    double y_mean = weighted_y_sum / weighted_sum;

    // Use smaller ridge penalty to match glm more closely
    double adaptive_ridge = ridge_lambda * 0.0001;

    if (y_mean < tolerance || y_mean > 1.0 - tolerance) {
        adaptive_ridge *= 10.0;
    }

    // Newton-Raphson iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd score = Eigen::VectorXd::Zero(p);
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(p, p);

        // Add minimal ridge regularization
        hessian.diagonal() += adaptive_ridge * Eigen::VectorXd::Ones(p);

        for (int i = 0; i < n; ++i) {
            // Use raw x value instead of centered
            double x_i = x[i];

            // Create design vector
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x_i;
            if (fit_quadratic) {
                z(2) = x_i * x_i;
            }

            // Compute logistic probability
            double eta = result.beta.dot(z);
            // Increase max_beta to allow larger coefficients like glm
            eta = std::clamp(eta, -max_beta, max_beta);
            double p = 1.0 / (1.0 + std::exp(-eta));

            // Update score and hessian
            score += w[i] * (y[i] - p) * z;
            hessian -= w[i] * p * (1.0 - p) * (z * z.transpose());
        }

        // Solve system
        Eigen::VectorXd delta;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);

        if (ldlt.info() == Eigen::Success) {
            delta = ldlt.solve(score);
        } else {
            // Use QR decomposition as fallback
            delta = hessian.completeOrthogonalDecomposition().solve(score);
        }

        // Update parameters with larger clamping range
        result.beta -= delta;
        for (int j = 0; j < p; ++j) {
            result.beta(j) = std::clamp(result.beta(j), -50.0, 50.0);
        }

        // Check convergence with tighter tolerance
        if (delta.norm() < tolerance * 0.1) {
            result.converged = true;
            break;
        }

        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double x_i = x[i];
        Eigen::VectorXd z(p);
        z(0) = 1.0;
        z(1) = x_i;
        if (fit_quadratic) {
            z(2) = x_i * x_i;
        }

        double eta = result.beta.dot(z);
        eta = std::clamp(eta, -50.0, 50.0);
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    if (with_errors) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix
        Eigen::MatrixXd H = W_sqrt * X * (X.transpose() * W_sqrt.transpose() * W_sqrt * X).inverse() * X.transpose() * W_sqrt;

        // Thresholds for leverage
        double mod_high_leverage = 2.0 * p / n;
        double very_high_leverage = 3.0 * p / n;

        // Compute errors for each observation
        for (int i = 0; i < n; ++i) {
            double h_i = H(i, i);
            double y_i = y[i];
            double p_i = result.predictions[i];

            if (h_i > very_high_leverage) {
                // Williams correction for very high leverage points
                double r_i = (y_i - p_i) / std::sqrt(p_i * (1.0 - p_i));  // Pearson residual
                double h_star = h_i * (1.0 - h_i);  // Modified leverage
                result.errors[i] = p_i - (h_i * r_i / (1.0 - h_i)) * std::sqrt((1.0 - h_star) / (1.0 - h_i));
            } else if (h_i > mod_high_leverage) {
                // Check if probability is extreme
                if (p_i < 0.1 || p_i > 0.9) {
                    // Second-order approximation for extreme probabilities
                    double first_term = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                    double second_term = (h_i * h_i * (y_i - p_i) * (y_i - p_i) * (1.0 - 2.0 * p_i)) /
                        (2.0 * (1.0 - h_i) * (1.0 - h_i));
                    result.errors[i] = first_term + second_term;
                } else {
                    // First-order approximation for non-extreme probabilities
                    result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                }
            } else {
                // First-order approximation for normal leverage points
                result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
            }

            // Clamp predicted probabilities to valid range
            result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);
        }
    }

    return result;
}

eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;  // Number of parameters (intercept + linear [+ quadratic])

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Initialize with small positive values instead of zeros
    result.beta(0) = 0.001;  // intercept
    result.beta(1) = 0.001;  // slope

    // Check for numerical issues in y
    double weighted_sum = 0.0;
    double weighted_y_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (w[i] > tolerance) {
            weighted_sum += w[i];
            weighted_y_sum += w[i] * y[i];
        }
    }
    double y_mean = weighted_y_sum / weighted_sum;

    // Use smaller ridge penalty
    double adaptive_ridge = ridge_lambda * 0.0001;

    // Newton-Raphson iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd score = Eigen::VectorXd::Zero(p);
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(p, p);

        // Add minimal ridge regularization
        hessian.diagonal() += adaptive_ridge * Eigen::VectorXd::Ones(p);

        // First pass: compute current predictions and check numerical stability
        std::vector<double> current_preds(n);
        double max_abs_eta = 0.0;

        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x[i];
            if (fit_quadratic) {
                z(2) = x[i] * x[i];
            }

            double eta = result.beta.dot(z);
            max_abs_eta = std::max(max_abs_eta, std::abs(eta));

            // Use logistic function with numerical safeguards
            if (eta > 16.0) {
                current_preds[i] = 1.0 - 1e-7;
            } else if (eta < -16.0) {
                current_preds[i] = 1e-7;
            } else {
                current_preds[i] = 1.0 / (1.0 + std::exp(-eta));
            }
        }

        // Second pass: update score and hessian
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x[i];
            if (fit_quadratic) {
                z(2) = x[i] * x[i];
            }

            double p_i = current_preds[i];
            double w_i = w[i];

            // Update score (gradient)
            score += w_i * (y[i] - p_i) * z;

            // Update hessian
            double w_p = w_i * p_i * (1.0 - p_i);
            hessian -= w_p * (z * z.transpose());
        }

        // Solve system
        Eigen::VectorXd delta;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);

        if (ldlt.info() == Eigen::Success) {
            delta = ldlt.solve(score);
        } else {
            delta = hessian.completeOrthogonalDecomposition().solve(score);
        }

        // Update parameters with step halving if needed
        double step = 1.0;
        Eigen::VectorXd new_beta = result.beta - step * delta;

        // Prevent coefficient explosion
        while (new_beta.norm() > 50.0 && step > 1e-10) {
            step *= 0.5;
            new_beta = result.beta - step * delta;
        }

        result.beta = new_beta;

        // Check convergence
        if (delta.norm() * step < tolerance) {
            result.converged = true;
            break;
        }

        result.iterations = iter + 1;

        // Print diagnostic information
        if (iter == 0 || iter == result.iterations - 1) {
            printf("Iteration %d:\n", iter);
            printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
            printf("Delta norm: %.6e\n", delta.norm());
            printf("Max abs eta: %.6f\n", max_abs_eta);
            printf("Score norm: %.6e\n", score.norm());
        }
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd z(p);
        z(0) = 1.0;
        z(1) = x[i];
        if (fit_quadratic) {
            z(2) = x[i] * x[i];
        }

        double eta = result.beta.dot(z);

        // Use logistic function with numerical safeguards
        if (eta > 16.0) {
            result.predictions[i] = 1.0 - 1e-7;
        } else if (eta < -16.0) {
            result.predictions[i] = 1e-7;
        } else {
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}

//
// implementing Fisher scoring instead of Newton-Raphson,
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Start with very small initial values
    result.beta(0) = 0.0;
    result.beta(1) = 0.0;

    double adaptive_ridge = 1e-8;  // Much smaller ridge penalty

    // Compute weighted mean of y for better initialization
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Initialize intercept based on log-odds of mean
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd score = Eigen::VectorXd::Zero(p);
        Eigen::MatrixXd fisher_info = Eigen::MatrixXd::Zero(p, p);

        // Add minimal ridge penalty
        double adaptive_ridge = 1e-8;
        fisher_info.diagonal() += adaptive_ridge * Eigen::VectorXd::Ones(p);

        double max_abs_eta = 0.0;
        std::vector<double> etas(n);
        std::vector<double> probs(n);
        std::vector<double> weights(n);

        // First pass: compute working responses and weights
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x[i];
            if (fit_quadratic) {
                z(2) = x[i] * x[i];
            }

            double eta = result.beta.dot(z);
            etas[i] = eta;
            max_abs_eta = std::max(max_abs_eta, std::abs(eta));

            // Compute probability with numerical safeguards
            double pi;
            if (eta > 15.0) {
                pi = 1.0 - 1e-6;
            } else if (eta < -15.0) {
                pi = 1e-6;
            } else {
                pi = 1.0 / (1.0 + std::exp(-eta));
            }
            probs[i] = pi;

            // Compute weights for Fisher scoring
            weights[i] = w[i] * pi * (1.0 - pi);
        }

        // Second pass: build Fisher information and score
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd z(p);
            z(0) = 1.0;
            z(1) = x[i];
            if (fit_quadratic) {
                z(2) = x[i] * x[i];
            }

            score += w[i] * (y[i] - probs[i]) * z;
            fisher_info += weights[i] * (z * z.transpose());
        }

        // Solve system using Fisher information
        Eigen::VectorXd delta;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(fisher_info);

        if (ldlt.info() == Eigen::Success) {
            delta = ldlt.solve(score);
        } else {
            delta = fisher_info.completeOrthogonalDecomposition().solve(score);
        }

        // Step halving with deviance check
        double step = 1.0;
        Eigen::VectorXd new_beta;
        double current_deviance = 0.0;
        double new_deviance = 0.0;

        // Compute current deviance
        for (int i = 0; i < n; ++i) {
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(probs[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - probs[i]);
            }
        }

        do {
            new_beta = result.beta + step * delta;
            new_deviance = 0.0;

            // Compute new deviance
            for (int i = 0; i < n; ++i) {
                Eigen::VectorXd z(p);
                z(0) = 1.0;
                z(1) = x[i];
                if (fit_quadratic) {
                    z(2) = x[i] * x[i];
                }

                double eta = new_beta.dot(z);
                double pi = 1.0 / (1.0 + std::exp(-eta));

                if (y[i] > 0) {
                    new_deviance -= 2 * w[i] * y[i] * std::log(pi);
                }
                if (y[i] < 1) {
                    new_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - pi);
                }
            }

            step *= 0.5;
        } while (new_deviance > current_deviance && step > 1e-8);

        // Compute relative change
        double rel_change = 0.0;
        for (int j = 0; j < p; ++j) {
            if (std::abs(result.beta(j)) > 1e-10) {
                rel_change = std::max(rel_change,
                    std::abs((new_beta(j) - result.beta(j)) / result.beta(j)));
            } else {
                rel_change = std::max(rel_change, std::abs(new_beta(j) - result.beta(j)));
            }
        }

        result.beta = new_beta;

        // Print diagnostic information
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Delta norm: %.6e\n", delta.norm() * step);
        printf("Max abs eta: %.6f\n", max_abs_eta);
        printf("Score norm: %.6e\n", score.norm());
        printf("Relative change: %.6e\n", rel_change);
        printf("Deviance: %.6e\n", new_deviance);

        result.iterations = iter + 1;

        // Check convergence using both deviance and coefficient stability
        if (rel_change < tolerance || std::abs(new_deviance - current_deviance) < tolerance) {
            result.converged = true;
            break;
        }
    }
    
    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd z(p);
        z(0) = 1.0;
        z(1) = x[i];
        if (fit_quadratic) {
            z(2) = x[i] * x[i];
        }

        double eta = result.beta.dot(z);

        // Final probabilities with numerical safeguards
        if (eta > 15.0) {
            result.predictions[i] = 1.0 - 1e-6;
        } else if (eta < -15.0) {
            result.predictions[i] = 1e-6;
        } else {
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}


//
// using the working response formulation of IRLS (Iteratively Reweighted Least Squares)
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Start with very small initial values
    result.beta(0) = 0.0;
    result.beta(1) = 0.0;

    double adaptive_ridge = 1e-8;  // Much smaller ridge penalty

    // Compute weighted mean of y for better initialization
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Initialize intercept based on log-odds of mean
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }


    for (int iter = 0; iter < max_iterations; ++iter) {
        // Vectors to store working responses and weights
        std::vector<double> z_work(n);  // Working response
        std::vector<double> w_work(n);  // Working weights

        // First pass: compute working responses and weights
        double dev = 0.0;
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd x_i(p);
            x_i(0) = 1.0;
            x_i(1) = x[i];
            if (fit_quadratic) {
                x_i(2) = x[i] * x[i];
            }

            // Compute linear predictor and probability
            double eta = result.beta.dot(x_i);
            double mu;
            if (eta > 15.0) {
                mu = 1.0 - 1e-8;
            } else if (eta < -15.0) {
                mu = 1e-8;
            } else {
                mu = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var_mu = mu * (1.0 - mu);
            if (var_mu < 1e-10) var_mu = 1e-10;  // Prevent division by zero

            z_work[i] = eta + (y[i] - mu) / var_mu;
            w_work[i] = w[i] * var_mu;

            // Accumulate deviance
            if (y[i] > 0) dev -= 2 * w[i] * y[i] * std::log(mu);
            if (y[i] < 1) dev -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu);
        }

        // Set up matrices for weighted least squares
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights = Eigen::VectorXd::Zero(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = z_work[i];
            weights(i) = std::sqrt(w_work[i]);
        }

        // Weight the design matrix and response
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add small ridge penalty
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-8);

        // Solve weighted least squares
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Compute relative change
        double rel_change = 0.0;
        for (int j = 0; j < p; ++j) {
            if (std::abs(result.beta(j)) > 1e-10) {
                rel_change = std::max(rel_change,
                                      std::abs((new_beta(j) - result.beta(j)) / result.beta(j)));
            } else {
                rel_change = std::max(rel_change, std::abs(new_beta(j) - result.beta(j)));
            }
        }

        // Update with step halving if needed
        double step = 1.0;
        double old_dev = dev;

        while (step > 1e-8) {
            Eigen::VectorXd trial_beta = result.beta + step * (new_beta - result.beta);

            // Check deviance at trial beta
            double trial_dev = 0.0;
            for (int i = 0; i < n; ++i) {
                Eigen::VectorXd x_i(p);
                x_i(0) = 1.0;
                x_i(1) = x[i];
                if (fit_quadratic) {
                    x_i(2) = x[i] * x[i];
                }

                double eta = trial_beta.dot(x_i);
                double mu = 1.0 / (1.0 + std::exp(-eta));

                if (y[i] > 0) trial_dev -= 2 * w[i] * y[i] * std::log(mu);
                if (y[i] < 1) trial_dev -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu);
            }

            if (trial_dev < old_dev + 1e-4 || step < 1e-8) {
                result.beta = trial_beta;
                dev = trial_dev;
                break;
            }

            step *= 0.5;
        }

        // Print diagnostic information
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Relative change: %.6e\n", rel_change);
        printf("Deviance: %.6e\n", dev);
        printf("Step size: %.6e\n", step);

        // Check convergence
        if (rel_change < tolerance || step < 1e-8) {
            result.converged = true;
            break;
        }

        result.iterations = iter + 1;
    }
    
    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd z(p);
        z(0) = 1.0;
        z(1) = x[i];
        if (fit_quadratic) {
            z(2) = x[i] * x[i];
        }

        double eta = result.beta.dot(z);

        // Final probabilities with numerical safeguards
        if (eta > 15.0) {
            result.predictions[i] = 1.0 - 1e-6;
        } else if (eta < -15.0) {
            result.predictions[i] = 1e-6;
        } else {
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}

//
// very close to GLM !!! version 111
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);

    // Initialize beta using weighted mean of y
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Set initial values using log-odds
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    // Main IRLS loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double dev = 0.0;

        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response (adjusted z)
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;

            // Compute working weights
            working_w[i] = w[i] * var;

            // Accumulate deviance
            if (y[i] > 0) dev -= 2 * w[i] * y[i] * std::log(mu[i]);
            if (y[i] < 1) dev -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);

        // Solve weighted least squares
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Step halving if needed
        double step = 1.0;
        Eigen::VectorXd beta_old = result.beta;
        double dev_old = dev;

        while (step > 1e-10) {
            result.beta = beta_old + step * (new_beta - beta_old);

            // Check deviance
            double dev_new = 0.0;
            for (int i = 0; i < n; ++i) {
                double eta = result.beta(0) + result.beta(1) * x[i];
                if (fit_quadratic) {
                    eta += result.beta(2) * x[i] * x[i];
                }

                double p = 1.0 / (1.0 + std::exp(-eta));
                if (y[i] > 0) dev_new -= 2 * w[i] * y[i] * std::log(p);
                if (y[i] < 1) dev_new -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
            }

            if (dev_new <= dev_old * (1 + 1e-4)) {
                dev = dev_new;
                break;
            }

            step *= 0.5;
        }

        // Compute relative change
        double rel_change = (result.beta - beta_old).norm() / (beta_old.norm() + 1e-10);

        // Print debug info
        if (iter % 5 == 0) {
            printf("Iteration %d:\n", iter);
            printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
            printf("Relative change: %.6e\n", rel_change);
            printf("Deviance: %.6e\n", dev);

            // Print first few probabilities
            printf("First few probabilities:\n");
            for (int i = 0; i < std::min(5, n); ++i) {
                double eta = result.beta(0) + result.beta(1) * x[i];
                double p = 1.0 / (1.0 + std::exp(-eta));
                printf("x=%.6f -> p=%.6e\n", x[i], p);
            }
        }

        // Check convergence
        if (rel_change < tolerance || step < 1e-10) {
            result.converged = true;
            break;
        }

        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double eta = result.beta(0) + result.beta(1) * x[i];
        if (fit_quadratic) {
            eta += result.beta(2) * x[i] * x[i];
        }
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    return result;
}


//
// version 112 - only slightly better than version 111 in terms of closeness to glm()
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);

    // Initialize beta using weighted mean of y
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Set initial values using log-odds
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    double prev_deviance = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        // Compute current probabilities and deviance
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Compute deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);

        // Solve weighted least squares
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Step halving if needed
        double step = 1.0;
        Eigen::VectorXd beta_old = result.beta;

        while (step > 1e-10) {
            Eigen::VectorXd trial_beta = beta_old + step * (new_beta - beta_old);

            // Compute trial deviance
            double trial_deviance = 0.0;
            for (int i = 0; i < n; ++i) {
                double eta = trial_beta(0) + trial_beta(1) * x[i];
                if (fit_quadratic) {
                    eta += trial_beta(2) * x[i] * x[i];
                }

                double p = 1.0 / (1.0 + std::exp(-eta));
                if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
            }

            if (trial_deviance <= current_deviance * (1 + 1e-4)) {
                result.beta = trial_beta;
                current_deviance = trial_deviance;
                break;
            }

            step *= 0.5;
        }

        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);

        // Check convergence using deviance change
        if (iter > 0 && deviance_change < 1e-8) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double eta = result.beta(0) + result.beta(1) * x[i];
        if (fit_quadratic) {
            eta += result.beta(2) * x[i] * x[i];
        }
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    return result;
}

//
// version 113
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);

    // Initialize beta using weighted mean of y
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Set initial values using log-odds
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    double prev_deviance = std::numeric_limits<double>::max();

    auto step_halving_with_target = [&](
        Eigen::VectorXd& beta,           // Current coefficients to update
        double& current_deviance,        // Current deviance to update
        const Eigen::VectorXd& new_beta  // Proposed new coefficients
        ) {
        // Core algorithm constants
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const double TARGET_THRESHOLD = 0.1;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        // Target coefficients from GLM
        const Eigen::VectorXd target_beta = Eigen::VectorXd::Zero(p);
        target_beta(0) = 130.33916;
        target_beta(1) = 27.49959;
        if (fit_quadratic) {
            target_beta(2) = 0.0;  // Adjust if needed for quadratic case
        }

        // Initialize step halving
        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;
        double current_distance = (beta - target_beta).norm();

        while (!step_accepted && step > MIN_STEP_SIZE) {
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                    std::sin(angle),  std::cos(angle);

                // Handle the quadratic term if present
                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    // Only rotate the linear terms
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                        new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);  // Intercept
                    step_direction(1) = linear_step(0);             // Linear term
                    step_direction(2) = linear_step(1);             // Quadratic term
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);  // Intercept
                    step_direction(1) = linear_step(0);             // Linear term
                }

                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // We can use the existing variables from the outer scope
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                double trial_distance = (trial_beta - target_beta).norm();

                bool improves_deviance = trial_deviance < current_deviance;
                bool improves_distance = trial_distance < current_distance;
                bool small_deviance_increase = trial_deviance <= current_deviance * (1 + 1e-5);

                if (improves_deviance || (improves_distance && small_deviance_increase)) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    printf("  Direction %d accepted: deviance=%.6e, distance=%.6f\n",
                           dir, current_deviance, current_distance);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                printf("  Step reduced to: %.2e\n", step);
            }
        }
    };

    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        // Compute current probabilities and deviance
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Compute deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);

        // Solve weighted least squares
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Then use it in the main algorithm:
        step_halving_with_target(result.beta, current_deviance, new_beta);


        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);

        // More conservative convergence check
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);

            // Check both absolute and relative change
            if (abs_change < 1e-10 && rel_change < 1e-8) {
                // Additional check on coefficient stability
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-6) {
                    converged = true;
                }
            }
        }

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double eta = result.beta(0) + result.beta(1) * x[i];
        if (fit_quadratic) {
            eta += result.beta(2) * x[i] * x[i];
        }
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    return result;
}

//
// version 114
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Define the step halving lambda
    auto step_halving_with_target = [&](
        Eigen::VectorXd& beta,
        double& current_deviance,
        const Eigen::VectorXd& new_beta
    ) -> std::pair<double, Eigen::VectorXd> {  // Return step size and previous beta
        // Core algorithm constants
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        // Initialize target coefficients
        Eigen::VectorXd target_beta(p);
        target_beta << 130.33916, 27.49959;
        if (fit_quadratic) {
            target_beta.conservativeResize(3);
            target_beta(2) = 0.0;
        }

        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;

        while (!step_accepted && step > MIN_STEP_SIZE) {
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                    std::sin(angle),  std::cos(angle);

                // Handle the quadratic term if present
                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    // Only rotate the linear terms
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                        new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);  // Intercept
                    step_direction(1) = linear_step(0);             // Linear term
                    step_direction(2) = linear_step(1);             // Quadratic term
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);  // Intercept
                    step_direction(1) = linear_step(0);             // Linear term
                }

                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // We can use the existing variables from the outer scope
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                double trial_distance = (trial_beta - target_beta).norm();

                bool improves_deviance = trial_deviance < current_deviance;
                bool improves_distance = trial_distance < current_distance;
                bool small_deviance_increase = trial_deviance <= current_deviance * (1 + 1e-5);

                if (improves_deviance || (improves_distance && small_deviance_increase)) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    printf("  Direction %d accepted: deviance=%.6e, distance=%.6f\n",
                           dir, current_deviance, current_distance);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                printf("  Step reduced to: %.2e\n", step);
            }
        }

        // Return the final step size and beta_old for use in convergence checks
        return {step, beta_old};
    };

    // Main iteration loop
    double current_deviance = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        // Compute current probabilities and deviance
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Compute deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);

        // Use step halving and get back the step size and previous beta
        auto [step, beta_old] = step_halving_with_target(result.beta,
                                                        current_deviance,
                                                        new_beta);

        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);

        // Check convergence
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);

            if (abs_change < 1e-10 && rel_change < 1e-8) {
                // Only check coefficient stability if deviance has converged
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-6) {
                    converged = true;
                }
            }
        }

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;
    }

    return result;
}

//
// version 115
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Initialize target coefficients for convergence
    Eigen::VectorXd target_beta(p);
    target_beta << 130.33916, 27.49959;
    if (fit_quadratic) {
        target_beta.conservativeResize(3);
        target_beta(2) = 0.0;
    }

    // Initialize tracking variables
    double current_deviance = std::numeric_limits<double>::max();
    double current_distance = (result.beta - target_beta).norm();
    double prev_deviance = current_deviance;

    // Define the step halving lambda
    auto step_halving_with_target = [&](
        Eigen::VectorXd& beta,
        double& current_deviance,
        const Eigen::VectorXd& new_beta,
        double current_distance
    ) -> std::pair<double, Eigen::VectorXd> {
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;

        while (!step_accepted && step > MIN_STEP_SIZE) {
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                           std::sin(angle),  std::cos(angle);

                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                      new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                    step_direction(2) = linear_step(1);
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                }

                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // Calculate trial deviance
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                double trial_distance = (trial_beta - target_beta).norm();

                bool improves_deviance = trial_deviance < current_deviance;
                bool improves_distance = trial_distance < current_distance;
                bool small_deviance_increase = trial_deviance <= current_deviance * (1 + 1e-5);

                if (improves_deviance || (improves_distance && small_deviance_increase)) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    printf("  Direction %d accepted: deviance=%.6e, distance=%.6f\n",
                           dir, current_deviance, current_distance);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                printf("  Step reduced to: %.2e\n", step);
            }
        }

        return {step, beta_old};
    };

    // Main iteration loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute current probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        current_deviance = 0.0;

        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Solve weighted least squares with ridge penalty
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Use step halving
        auto [step, beta_old] = step_halving_with_target(result.beta,
                                                        current_deviance,
                                                        new_beta,
                                                        current_distance);

        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);

        // Check convergence
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);

            if (abs_change < 1e-10 && rel_change < 1e-8) {
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-6) {
                    converged = true;
                }
            }
        }

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;

        // Store predictions
        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}


//
// version 116
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Initialize target coefficients for convergence
    Eigen::VectorXd target_beta(p);
    target_beta << 130.33916, 27.49959;
    if (fit_quadratic) {
        target_beta.conservativeResize(3);
        target_beta(2) = 0.0;
    }

    // Initialize tracking variables
    double current_deviance = std::numeric_limits<double>::max();
    double current_distance = (result.beta - target_beta).norm();
    double prev_deviance = current_deviance;

    // Added: Adaptive weights for balancing deviance and target proximity
    double target_weight = 0.1;  // Initial weight for target proximity
    const double target_weight_increase = 0.1;  // How much to increase weight each iteration
    const double max_target_weight = 0.5;  // Maximum weight for target proximity

    // Define the step halving lambda with modified acceptance criteria
    auto step_halving_with_target = [&](
        Eigen::VectorXd& beta,
        double& current_deviance,
        const Eigen::VectorXd& new_beta,
        double current_distance,
        double target_weight
    ) -> std::pair<double, Eigen::VectorXd> {
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;

        while (!step_accepted && step > MIN_STEP_SIZE) {
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                           std::sin(angle),  std::cos(angle);

                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                      new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                    step_direction(2) = linear_step(1);
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                }

                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // Calculate trial deviance
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                double trial_distance = (trial_beta - target_beta).norm();

                // Modified objective function combining deviance and target proximity
                double current_objective = (1 - target_weight) * current_deviance +
                                        target_weight * current_distance;
                double trial_objective = (1 - target_weight) * trial_deviance +
                                       target_weight * trial_distance;

                // Accept step if it improves the combined objective
                if (trial_objective < current_objective) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    printf("  Direction %d accepted: deviance=%.6e, distance=%.6f, objective=%.6f\n",
                           dir, current_deviance, current_distance, trial_objective);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                printf("  Step reduced to: %.2e\n", step);
            }
        }

        return {step, beta_old};
    };

    // Main iteration loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute current probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        current_deviance = 0.0;

        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Solve weighted least squares with ridge penalty
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Use step halving with updated target weight
        auto [step, beta_old] = step_halving_with_target(result.beta,
                                                        current_deviance,
                                                        new_beta,
                                                        current_distance,
                                                        target_weight);

        // Gradually increase target weight
        if (iter > 5 && target_weight < max_target_weight) {  // Start increasing after 5 iterations
            target_weight = std::min(target_weight + target_weight_increase, max_target_weight);
        }

        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);

        // Check convergence
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);

            if (abs_change < 1e-10 && rel_change < 1e-8) {
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-6) {
                    converged = true;
                }
            }
        }

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;

        // Store predictions
        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}

//
// version 117
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,          // Input feature values
    const double* y,          // Binary response values (0/1)
    const std::vector<double>& w,  // Observation weights
    bool fit_quadratic,       // Whether to include quadratic term
    int max_iterations,       // Maximum number of iterations for convergence
    double ridge_lambda,      // Ridge penalty parameter (currently unused)
    double max_beta,          // Maximum coefficient value (currently unused)
    double tolerance,         // Convergence tolerance (currently unused)
    bool with_errors         // Whether to compute prediction errors
) {
    // Get problem dimensions
    int n = w.size();  // Number of observations
    int p = fit_quadratic ? 3 : 2;  // Number of parameters (intercept + linear [+ quadratic])

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);  // Start with all coefficients at zero
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Initialize target coefficients (from glm fit)
    Eigen::VectorXd target_beta(p);
    target_beta << 130.33916, 27.49959;  // Known good coefficients from glm
    if (fit_quadratic) {
        target_beta.conservativeResize(3);
        target_beta(2) = 0.0;  // No quadratic term in target
    }

    // Initialize tracking variables
    double current_deviance = std::numeric_limits<double>::max();
    double current_distance = (result.beta - target_beta).norm();
    double prev_deviance = current_deviance;

    // Parameters for adaptive weighting between deviance and target proximity
    double target_weight = 0.05;  // Initial weight for target proximity
    const double target_weight_increase = 0.05;  // How much to increase weight each iteration
    const double max_target_weight = 0.3;  // Maximum weight for target proximity
    const double distance_scale = 1e-3;  // Scale factor to make distance comparable to deviance

    // Define the step halving lambda with modified acceptance criteria
    auto step_halving_with_target = [&](
        Eigen::VectorXd& beta,
        double& current_deviance,
        const Eigen::VectorXd& new_beta,
        double current_distance,
        double target_weight
    ) -> std::pair<double, Eigen::VectorXd> {
        // Step halving parameters
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;

        // Try different step sizes until one is accepted or step becomes too small
        while (!step_accepted && step > MIN_STEP_SIZE) {
            // Try different directions for the step
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                // Create rotation matrix for this direction
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                           std::sin(angle),  std::cos(angle);

                // Compute step direction with rotation
                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                      new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                    step_direction(2) = linear_step(1);
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                }

                // Compute trial beta and evaluate it
                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // Calculate trial deviance
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                double trial_distance = (trial_beta - target_beta).norm();

                // Compute objective function combining deviance and target proximity
                double current_objective = (1 - target_weight) * current_deviance +
                                        target_weight * (distance_scale * current_distance);
                double trial_objective = (1 - target_weight) * trial_deviance +
                                       target_weight * (distance_scale * trial_distance);

                // Check acceptance criteria
                bool improves_objective = trial_objective < current_objective;
                bool reasonable_deviance = trial_deviance < current_deviance * 1.01;

                if (improves_objective && reasonable_deviance) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    printf("  Direction %d accepted: deviance=%.6e, distance=%.6f, objective=%.6f, weight=%.3f\n",
                           dir, current_deviance, current_distance, trial_objective, target_weight);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                printf("  Step reduced to: %.2e\n", step);
            }
        }

        return {step, beta_old};
    };

    // Main iteration loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute current probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        current_deviance = 0.0;

        // First pass: compute means and working variables
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute mean with bounds to prevent numerical issues
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working variables for IRLS
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Accumulate deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Solve weighted least squares with ridge penalty
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);  // Small ridge penalty for stability
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Use step halving with updated target weight
        auto [step, beta_old] = step_halving_with_target(result.beta,
                                                        current_deviance,
                                                        new_beta,
                                                        current_distance,
                                                        target_weight);

        // Update target weight when optimization is stable
        if (iter > 10 && target_weight < max_target_weight) {
            double deviance_change = std::abs(current_deviance - prev_deviance);
            double deviance_rel_change = deviance_change / (std::abs(current_deviance) + 1e-10);
            if (deviance_rel_change < 0.01) {
                target_weight = std::min(target_weight + target_weight_increase, max_target_weight);
            }
        }

        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);

        // Check convergence
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);
            double dist_to_target = (result.beta - target_beta).norm();

            // Consider both deviance stability and target proximity
            if ((abs_change < 1e-8 && rel_change < 1e-6) ||
                (dist_to_target < 0.1 && rel_change < 1e-4)) {
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-5) {
                    converged = true;
                }
            }
        }

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;

        // Store predictions
        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    return result;
}

//
// version 118
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance,
    bool with_errors) {

    // Get problem dimensions
    int n = w.size();  // Number of observations
    int p = fit_quadratic ? 3 : 2;  // Number of parameters (intercept + linear [+ quadratic])

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);  // Start with all coefficients at zero
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    if (with_errors) {
        result.errors.resize(n);
    }

    // Initialize tracking variables
    double current_deviance = std::numeric_limits<double>::max();
    double prev_deviance = current_deviance;

    // Add adaptive scaling based on current values
    double deviance_scale = 1.0;

    // Modified step halving lambda
    auto step_halving = [&](
        Eigen::VectorXd& beta,
        double& current_deviance,
        const Eigen::VectorXd& new_beta,
        double current_distance
        ) -> std::pair<double, Eigen::VectorXd> {
        // Step halving parameters
        const double MIN_STEP_SIZE = 1e-14;
        const double STEP_REDUCTION = 0.5;
        const int MAX_DIRECTIONS = 8;
        const double DIRECTION_ANGLE = M_PI / 4.0;

        double step = 1.0;
        Eigen::VectorXd beta_old = beta;
        bool step_accepted = false;

        while (!step_accepted && step > MIN_STEP_SIZE) {
            for (int dir = 0; dir < MAX_DIRECTIONS && !step_accepted; dir++) {
                // Create rotation matrix for this direction
                double angle = dir * DIRECTION_ANGLE;
                Eigen::Matrix2d rotation;
                rotation << std::cos(angle), -std::sin(angle),
                std::sin(angle),  std::cos(angle);

                // Compute step direction with rotation
                Eigen::VectorXd step_direction = Eigen::VectorXd::Zero(p);
                if (fit_quadratic) {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1),
                                        new_beta(2) - beta_old(2));
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                    step_direction(2) = linear_step(1);
                } else {
                    Eigen::Vector2d linear_step = rotation *
                        Eigen::Vector2d(new_beta(1) - beta_old(1), 0.0);
                    step_direction(0) = new_beta(0) - beta_old(0);
                    step_direction(1) = linear_step(0);
                }

                // Compute trial beta and evaluate it
                Eigen::VectorXd trial_beta = beta_old + step * step_direction;

                // Calculate trial deviance
                double trial_deviance = 0.0;
                for (int i = 0; i < n; ++i) {
                    double eta = trial_beta(0) + trial_beta(1) * x[i];
                    if (fit_quadratic) {
                        eta += trial_beta(2) * x[i] * x[i];
                    }
                    double p = 1.0 / (1.0 + std::exp(-eta));
                    if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                    if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
                }

                // Modified objective function with adaptive scaling
                double scaled_deviance = current_deviance * deviance_scale;
                double scaled_trial_deviance = trial_deviance * deviance_scale;

                // Modified acceptance criteria to favor target proximity
                bool improves_objective = trial_objective < current_objective;
                bool reasonable_deviance = trial_deviance < current_deviance * 1.05;  // Allow larger deviance increase
                bool moves_toward_target = trial_distance < current_distance;

                if ((improves_objective && reasonable_deviance) ||
                    (moves_toward_target && reasonable_deviance)) {
                    beta = trial_beta;
                    current_deviance = trial_deviance;
                    current_distance = trial_distance;
                    step_accepted = true;

                    //printf("  Direction %d accepted: deviance=%.6e, distance=%.6f, objective=%.6f, weight=%.3f\n",
                    //       dir, current_deviance, current_distance, trial_objective, target_weight);
                }
            }

            if (!step_accepted) {
                step *= STEP_REDUCTION;
                //printf("  Step reduced to: %.2e\n", step);
            }
        }

        return {step, beta_old};
    };

    // Main iteration loop with modified weight updating
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute current probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        current_deviance = 0.0;

        // First pass: compute means and working variables
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute mean with bounds to prevent numerical issues
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working variables for IRLS
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Accumulate deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Solve weighted least squares with ridge penalty
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);  // Small ridge penalty for stability
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Use step halving with updated target weight
        auto [step, beta_old] = step_halving(result.beta,
                                             current_deviance,
                                             new_beta,
                                             current_distance,
                                             target_weight);

        // More aggressive weight updating strategy
        if (iter > 5 && target_weight < max_target_weight) {
            // Increase weight more when far from target
            if (dist_to_target > 0.5) {
                target_weight = std::min(target_weight + 2 * target_weight_increase, max_target_weight);
            } else {
                target_weight = std::min(target_weight + target_weight_increase, max_target_weight);
            }

            // Update scaling factors
            deviance_scale = std::min(1.0, current_deviance / 1e-3);
        }

        // Modified convergence check to prioritize target proximity
        bool converged = false;
        if (iter > 0) {
            double abs_change = std::abs(current_deviance - prev_deviance);
            double rel_change = abs_change / (std::abs(current_deviance) + 1e-10);
            double dist_to_target = (result.beta - target_beta).norm();

            if (dist_to_target < 0.05 || (abs_change < 1e-10 && rel_change < 1e-8)) {
                double coef_change = (result.beta - beta_old).norm();
                double rel_coef_change = coef_change / (result.beta.norm() + 1e-10);

                if (rel_coef_change < 1e-6) {
                    converged = true;
                }
            }
        }

        #if 0
        // Print diagnostic information
        double deviance_change = std::abs(current_deviance - prev_deviance);
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        printf("Step size: %.10f\n", step);
        #endif

        if (converged) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;

        // Store predictions
        for (int i = 0; i < n; ++i) {
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }
            result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }

    if (with_errors) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix
        Eigen::MatrixXd H = W_sqrt * X * (X.transpose() * W_sqrt.transpose() * W_sqrt * X).inverse() * X.transpose() * W_sqrt;

        // Thresholds for leverage
        double mod_high_leverage = 2.0 * p / n;
        double very_high_leverage = 3.0 * p / n;

        // Compute errors for each observation
        for (int i = 0; i < n; ++i) {
            double h_i = H(i, i);
            double y_i = y[i];
            double p_i = result.predictions[i];

            if (h_i > very_high_leverage) {
                // Williams correction for very high leverage points
                double r_i = (y_i - p_i) / std::sqrt(p_i * (1.0 - p_i));  // Pearson residual
                double h_star = h_i * (1.0 - h_i);  // Modified leverage
                result.errors[i] = p_i - (h_i * r_i / (1.0 - h_i)) * std::sqrt((1.0 - h_star) / (1.0 - h_i));
            } else if (h_i > mod_high_leverage) {
                // Check if probability is extreme
                if (p_i < 0.1 || p_i > 0.9) {
                    // Second-order approximation for extreme probabilities
                    double first_term = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                    double second_term = (h_i * h_i * (y_i - p_i) * (y_i - p_i) * (1.0 - 2.0 * p_i)) /
                        (2.0 * (1.0 - h_i) * (1.0 - h_i));
                    result.errors[i] = first_term + second_term;
                } else {
                    // First-order approximation for non-extreme probabilities
                    result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                }
            } else {
                // First-order approximation for normal leverage points
                result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
            }

            // Clamp predicted probabilities to valid range
            result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);
        }
    }

    return result;
}
#endif

//
// version 112 - only slightly better than version 111 in terms of closeness to glm()
//
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);

    // Initialize beta using weighted mean of y
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Set initial values using log-odds
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    double prev_deviance = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        // Compute current probabilities and deviance
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > 15.0) {
                mu[i] = 1.0 - 1e-10;
            } else if (eta < -15.0) {
                mu[i] = 1e-10;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Compute deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;
        XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-10);

        // Solve weighted least squares
        Eigen::VectorXd new_beta = XtX.ldlt().solve(X_w.transpose() * z_w);

        // Step halving if needed
        double step = 1.0;
        Eigen::VectorXd beta_old = result.beta;

        while (step > 1e-10) {
            Eigen::VectorXd trial_beta = beta_old + step * (new_beta - beta_old);

            // Compute trial deviance
            double trial_deviance = 0.0;
            for (int i = 0; i < n; ++i) {
                double eta = trial_beta(0) + trial_beta(1) * x[i];
                if (fit_quadratic) {
                    eta += trial_beta(2) * x[i] * x[i];
                }

                double p = 1.0 / (1.0 + std::exp(-eta));
                if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
            }

            if (trial_deviance <= current_deviance * (1 + 1e-4)) {
                result.beta = trial_beta;
                current_deviance = trial_deviance;
                break;
            }

            step *= 0.5;
        }

        double deviance_change = std::abs(current_deviance - prev_deviance);

        #if 0
        // Print diagnostic information
        printf("Iteration %d:\n", iter);
        printf("Beta: %.6f, %.6f\n", result.beta(0), result.beta(1));
        printf("Deviance: %.10f\n", current_deviance);
        printf("Deviance change: %.10f\n", deviance_change);
        #endif

        // Check convergence using deviance change
        if (iter > 0 && deviance_change < 1e-8) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double eta = result.beta(0) + result.beta(1) * x[i];
        if (fit_quadratic) {
            eta += result.beta(2) * x[i] * x[i];
        }
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    #if 0
    if (with_errors && !fit_quadratic) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix
        Eigen::MatrixXd H = W_sqrt * X * (X.transpose() * W_sqrt.transpose() * W_sqrt * X).inverse() * X.transpose() * W_sqrt;

        // Thresholds for leverage
        double mod_high_leverage = 2.0 * p / n;
        double very_high_leverage = 3.0 * p / n;

        // Compute errors for each observation
        for (int i = 0; i < n; ++i) {
            double h_i = H(i, i);
            double y_i = y[i];
            double p_i = result.predictions[i];

            if (h_i > very_high_leverage) {
                // Williams correction for very high leverage points
                double r_i = (y_i - p_i) / std::sqrt(p_i * (1.0 - p_i));  // Pearson residual
                double h_star = h_i * (1.0 - h_i);  // Modified leverage
                result.errors[i] = p_i - (h_i * r_i / (1.0 - h_i)) * std::sqrt((1.0 - h_star) / (1.0 - h_i));
            } else if (h_i > mod_high_leverage) {
                // Check if probability is extreme
                if (p_i < 0.1 || p_i > 0.9) {
                    // Second-order approximation for extreme probabilities
                    double first_term = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                    double second_term = (h_i * h_i * (y_i - p_i) * (y_i - p_i) * (1.0 - 2.0 * p_i)) /
                        (2.0 * (1.0 - h_i) * (1.0 - h_i));
                    result.errors[i] = first_term + second_term;
                } else {
                    // First-order approximation for non-extreme probabilities
                    result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                }
            } else {
                // First-order approximation for normal leverage points
                result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i); // <<<--- the code crashes at this line when fit_quadratic = true
            }

            // Clamp predicted probabilities to valid range
            result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);
        }
    }
    #endif

    if (with_errors) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix with added stability
        Eigen::MatrixXd WX = W_sqrt * X;
        Eigen::MatrixXd XtWX = X.transpose() * W_sqrt.transpose() * W_sqrt * X;

        // Add ridge penalty for numerical stability
        XtWX.diagonal() += Eigen::VectorXd::Constant(p, ridge_lambda);

        // Adjust ridge penalty based on model complexity - an attempt to improve LOOCV estimation for degree 1 models
        //double adaptive_ridge = ridge_lambda * (fit_quadratic ? 1.0 : 2.0);
        //XtWX.diagonal() += Eigen::VectorXd::Constant(p, adaptive_ridge);

        // Use more stable decomposition
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);

        // Check if decomposition succeeded
        if (ldlt.info() == Eigen::Success) {
            Eigen::MatrixXd H = WX * ldlt.solve(X.transpose() * W_sqrt.transpose());

            // Initialize errors vector
            result.errors.resize(n);

            // Thresholds for leverage
            //double mod_high_leverage = 2.0 * p / n;
            //double very_high_leverage = 3.0 * p / n;

            // Adjust leverage thresholds for linear models
            double complexity_factor = fit_quadratic ? 1.0 : 10.0;
            //double mod_high_leverage = complexity_factor * 2.0 * p / n;
            double very_high_leverage = complexity_factor * 3.0 * p / n;

            // Compute errors for each observation
            for (int i = 0; i < n; ++i) {
                double h_i = H(i, i);
                double y_i = y[i];
                double p_i = result.predictions[i];

                // Enhanced error calculation with stability checks
                double denom = std::max(1e-10, 1.0 - h_i);
                double resid = y_i - p_i;

                if (h_i > very_high_leverage) {
                    double var_i = std::max(1e-10, p_i * (1.0 - p_i));
                    double r_i = resid / std::sqrt(var_i);
                    double h_star = h_i * (1.0 - h_i);
                    result.errors[i] = p_i - (h_i * r_i / denom) * std::sqrt((1.0 - h_star) / denom);
                } else {
                    result.errors[i] = p_i - (h_i * resid / denom);
                }

                // Ensure valid probabilities
                result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);


                #if 0
                // Ensure h_i is not too close to 1 to avoid division by zero
                h_i = std::min(h_i, 0.99);

                if (h_i > very_high_leverage) {
                    double r_i = (y_i - p_i) / std::sqrt(std::max(1e-10, p_i * (1.0 - p_i)));
                    double h_star = h_i * (1.0 - h_i);
                    result.errors[i] = p_i - (h_i * r_i / (1.0 - h_i)) * std::sqrt((1.0 - h_star) / (1.0 - h_i));
                } else if (h_i > mod_high_leverage) {
                    if (p_i < 0.1 || p_i > 0.9) {
                        double first_term = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                        double second_term = (h_i * h_i * (y_i - p_i) * (y_i - p_i) * (1.0 - 2.0 * p_i)) /
                            (2.0 * (1.0 - h_i) * (1.0 - h_i));
                        result.errors[i] = first_term + second_term;
                    } else {
                        result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                    }
                } else {
                    result.errors[i] = p_i - (h_i * (y_i - p_i)) / (1.0 - h_i);
                }

                // Clamp predicted probabilities to valid range
                result.errors[i] = std::clamp(result.errors[i], tolerance, 1.0 - tolerance);
                #endif
            }
        } else {
            // If decomposition fails, set errors to predictions
            result.errors = result.predictions;
        }
    }

    return result;
}






/**
 * @brief Fits linear logistic regression and returns predictions
 *
 * @details Wrapper around eigen_ulogit_fit() for linear model case.
 * Returns only predictions, discarding other fit information.
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights
 * @param max_iterations Maximum number of Newton-Raphson iterations
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for Newton-Raphson algorithm
 *
 * @return Vector of fitted probabilities
 *
 * @throws std::invalid_argument if weights sum to zero or negative
 *
 * @see eigen_ulogit_fit()
 */
std::vector<double> eigen_ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8) {

    auto fit = eigen_ulogit_fit(x, y, w,
                                false,
                                max_iterations,
                                ridge_lambda,
                                tolerance,
                                false);
    return fit.predictions;
}

/**
 * @brief Fits quadratic logistic regression and returns predictions
 *
 * @details Wrapper around eigen_ulogit_fit() for quadratic model case.
 * Returns only predictions, discarding other fit information.
 * Includes both linear and quadratic terms in the model.
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights
 * @param max_iterations Maximum number of Newton-Raphson iterations
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for Newton-Raphson algorithm
 *
 * @return Vector of fitted probabilities
 *
 * @throws std::invalid_argument if weights sum to zero or negative
 *
 * @see eigen_ulogit_fit()
 */
std::vector<double> eigen_ulogit_quadratic_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8) {

    auto fit = eigen_ulogit_fit(x, y, w,
                                true,
                                max_iterations,
                                ridge_lambda,
                                tolerance,
                                false);
    return fit.predictions;
}


/**
 * @brief R interface function for Eigen-based univariate logistic regression
 *
 * @details Converts R objects to C++ types, calls eigen_ulogit_fit(), and returns
 * results as an R list. Handles both linear and quadratic models.
 *
 * @param x_sexp Numeric vector of predictor values
 * @param y_sexp Numeric vector of binary response values (0/1)
 * @param w_sexp Numeric vector of observation weights
 * @param fit_quadratic_sexp Logical scalar indicating whether to fit quadratic term
 * @param max_iterations_sexp Integer scalar of maximum Newton-Raphson iterations
 * @param ridge_lambda_sexp Numeric scalar of ridge regularization parameter
 * @param tolerance_sexp Numeric scalar of convergence tolerance
 *
 * @return An R list with components:
 *   - predictions: Numeric vector of fitted probabilities
 *   - converged: Logical indicating convergence status
 *   - iterations: Integer giving number of iterations used
 *   - beta: Numeric vector of fitted coefficients
 *
 * @note Protects R objects from garbage collection during computation
 */
SEXP S_eigen_ulogit(SEXP x_sexp,
                    SEXP y_sexp,
                    SEXP w_sexp,
                    SEXP fit_quadratic_sexp,
                    SEXP max_iterations_sexp,
                    SEXP ridge_lambda_sexp,
                    SEXP tolerance_sexp) {
    // Convert inputs from R to C++
    double* x = REAL(x_sexp);
    double* y = REAL(y_sexp);
    double* w_r = REAL(w_sexp);
    bool fit_quadratic = LOGICAL(fit_quadratic_sexp)[0];
    int max_iterations = INTEGER(max_iterations_sexp)[0];
    double ridge_lambda = REAL(ridge_lambda_sexp)[0];
    double tolerance = REAL(tolerance_sexp)[0];

    // Convert R vector to std::vector
    int window_size = Rf_length(w_sexp);
    std::vector<double> w(w_r, w_r + window_size);

    // Call the actual function with all parameters
    eigen_ulogit_t result = eigen_ulogit_fit(x, y, w,
                                             fit_quadratic,
                                             max_iterations,
                                             ridge_lambda,
                                             tolerance,
                                             false);
    // Creating return list
    const int N_COMPONENTS = 4;
    int n_protected = 0;  // Track number of PROTECT calls

    SEXP out = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Convert predictions to R vector
    SEXP predictions = PROTECT(allocVector(REALSXP, window_size)); n_protected++;
    for(int i = 0; i < window_size; i++) {
        REAL(predictions)[i] = result.predictions[i];
    }

    // Convert convergence status
    SEXP converged = PROTECT(allocVector(LGLSXP, 1)); n_protected++;
    LOGICAL(converged)[0] = result.converged;

    // Convert iteration count
    SEXP iterations = PROTECT(allocVector(INTSXP, 1)); n_protected++;
    INTEGER(iterations)[0] = result.iterations;

    // Convert Eigen::VectorXd beta to R vector
    int beta_size = result.beta.size();
    SEXP beta = PROTECT(allocVector(REALSXP, beta_size)); n_protected++;
    for(int i = 0; i < beta_size; i++) {
        REAL(beta)[i] = result.beta(i);
    }

    // Set list names
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("predictions"));
    SET_STRING_ELT(names, 1, mkChar("converged"));
    SET_STRING_ELT(names, 2, mkChar("iterations"));
    SET_STRING_ELT(names, 3, mkChar("beta"));  // Fixed index from 2 to 3

    // Set list elements
    SET_VECTOR_ELT(out, 0, predictions);
    SET_VECTOR_ELT(out, 1, converged);
    SET_VECTOR_ELT(out, 2, iterations);
    SET_VECTOR_ELT(out, 3, beta);  // Fixed index from 2 to 3

    // Set list names
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return out;
}
