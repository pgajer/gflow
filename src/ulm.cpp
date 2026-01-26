#include <Eigen/Dense>
#include <vector>
#include <cmath>      // For fabs()
#include <algorithm>  // For std::find, std::clamp
#include <numeric>    // For std::accumulate
#include <limits>

#include "ulm.hpp"
#include "error_utils.h"

#include <R.h>
#include <Rinternals.h>

/**
 * @file ulm.cpp
 * @brief Weighted univariate linear model (ULM) with analytic LOOCV errors.
 *
 * This implementation fits a weighted linear regression with intercept:
 *   y_i ≈ beta0 + beta1 * (x_i - x_wmean),
 * where x_wmean and y_wmean are the weighted means. It returns fitted values
 * at each input x_i and per-point LOOCV squared errors using the standard
 * leverage-based formula:
 *
 *   e_{-i} = (y_i - yhat_i) / (1 - h_ii),   LOOCV_i = e_{-i}^2,
 *
 * with special handling for degenerate neighborhoods where the weighted
 * variance of x is (near) zero. In that case, the model reduces to a
 * weighted mean (local-constant), and the hat diagonal is:
 *
 *   h_ii = w_i / sum(w).
 *
 * Robustness and safety features:
 *   - Handles empty neighborhoods (n=0) safely.
 *   - Validates x/y pointers when n>0.
 *   - Validates weights for finiteness, non-negativity, and positive total weight.
 *   - Sizes output vectors (predictions and errors) consistently (prevents OOB writes).
 *   - Guards against division-by-zero / near-singular leverage computations.
 *   - Clamps leverage to [0, 1) for numerical stability.
 *
 * Notes:
 *   - This routine treats the provided weights as fixed. If used inside an
 *     IRLS robustification loop (e.g., cleveland_ulm), the returned LOOCV
 *     errors correspond to the final weighted WLS fit conditional on the final
 *     weights (a common approximation).
 */
ulm_t ulm(const double* x,
          const double* y,
          const std::vector<double>& w,
          bool y_binary,
          double epsilon)
{
    const int n_points = static_cast<int>(w.size());

    ulm_t results;
    results.slope = 0.0;
    results.x_wmean = 0.0;
    results.y_wmean = 0.0;

    // Empty neighborhood: valid upstream outcome; return empty safely.
    if (n_points == 0) {
        results.predictions.clear();
        results.errors.clear();
        return results;
    }

    // Defensive pointer validation when n > 0.
    if (x == nullptr || y == nullptr) {
        Rf_error("ulm(): null x/y pointer with n_points=%d", n_points);
    }

    // Validate weights and compute total weight.
    double total_weight = 0.0;
    for (int i = 0; i < n_points; ++i) {
        const double wi = w[i];
        if (!std::isfinite(wi)) {
            Rf_error("ulm(): non-finite weight w[%d]=%g", i, wi);
        }
        if (wi < 0.0) {
            Rf_error("ulm(): negative weight w[%d]=%g", i, wi);
        }
        total_weight += wi;
    }

    if (!std::isfinite(total_weight) || total_weight <= epsilon) {
        Rf_error("ulm(): invalid total_weight=%g (n_points=%d)", total_weight, n_points);
    }

    // Working copy of x for centering and leverage computation.
    std::vector<double> x_working(x, x + n_points);

    // Weighted means.
    double x_wmean = 0.0;
    double y_wmean = 0.0;

    for (int i = 0; i < n_points; ++i) {
        x_wmean += w[i] * x_working[i];
        y_wmean += w[i] * y[i];
    }
    x_wmean /= total_weight;
    y_wmean /= total_weight;

    results.x_wmean = x_wmean;
    results.y_wmean = y_wmean;

    // Center x around weighted mean and compute weighted sum of squares.
    double sum_wx_squared = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_working[i] -= x_wmean;
        sum_wx_squared += w[i] * x_working[i] * x_working[i];
    }

    // Fit slope (local linear). If degenerate, slope stays 0.
    double slope = 0.0;
    if (sum_wx_squared > epsilon) {
        double wxy_sum = 0.0;
        for (int i = 0; i < n_points; ++i) {
            wxy_sum += w[i] * x_working[i] * y[i];
        }
        slope = wxy_sum / sum_wx_squared;
    }
    results.slope = slope;

    // Size outputs (critical to avoid OOB writes).
    results.predictions.resize(n_points);
    results.errors.resize(n_points);

    // Predictions at original x.
    for (int i = 0; i < n_points; ++i) {
        double yhat = y_wmean + slope * (x[i] - x_wmean);
        if (y_binary) {
            yhat = std::clamp(yhat, 0.0, 1.0);
        }
        results.predictions[i] = yhat;
    }

    // Analytic LOOCV squared errors using hat diagonals.
    for (int i = 0; i < n_points; ++i) {

        double h_i = 0.0;

        if (sum_wx_squared <= epsilon) {
            // Degenerate-x case: local-constant weighted mean.
            h_i = (w[i] > 0.0) ? (w[i] / total_weight) : 0.0;
        } else {
            // Linear smoother hat diagonal with centered x:
            // h_i = w_i * (1/sum(w) + x_i^2 / sum(w * x^2))
            h_i = w[i] * (1.0 / total_weight + (x_working[i] * x_working[i]) / sum_wx_squared);
        }

        // Numerical hygiene.
        if (!std::isfinite(h_i)) h_i = 0.0;
        if (h_i < 0.0) h_i = 0.0;
        if (h_i > 0.999999) h_i = 0.999999;

        const double denom = 1.0 - h_i;

        if (denom > epsilon) {
            const double residual = (y[i] - results.predictions[i]) / denom;

            if (!std::isfinite(residual)) {
                // During debugging, fail fast to preserve the first bad context.
                Rf_error("ulm(): non-finite LOOCV residual at i=%d (denom=%g, h=%g, y=%g, yhat=%g)",
                         i, denom, h_i, y[i], results.predictions[i]);
            }

            results.errors[i] = residual * residual;
        } else {
            // Leverage too close to 1: LOOCV undefined/infinite.
            results.errors[i] = std::numeric_limits<double>::infinity();
        }
    }

    return results;
}


/**
 * @brief Results structure for univariate linear and polynomial models
 *
 * This structure holds the results of fitting either a linear (degree 1) or
 * quadratic (degree 2) univariate model to data, including predictions,
 * leave-one-out cross-validation errors, and model coefficients.
 *
 * @field predictions Vector of fitted values for each observation
 * @field errors Vector of leave-one-out cross-validation (LOOCV) prediction errors
 * @field coefficients Vector of model coefficients. For degree 1: [intercept, slope].
 *                     For degree 2: [intercept, linear term, quadratic term]
 * @field w Vector of observation weights used in the model fitting
 */
struct ulm_results {
    std::vector<double> predictions;
    std::vector<double> errors;
    std::vector<double> coefficients;
    std::vector<double> w;
};

/**
 * @brief Fits a univariate linear model using direct computation
 *
 * Implements fast, numerically stable fitting of degree 1 linear models
 * using direct computation methods instead of matrix operations.
 * Includes weighted least squares and leave-one-out cross-validation.
 *
 * @param x Pointer to predictor variable array
 * @param y Pointer to response variable array
 * @param w Vector of observation weights
 * @param y_binary Boolean indicating if response is binary (0/1)
 * @param epsilon Small number for numerical stability checks (default: 1e-10)
 *
 * @return ulm_results structure containing:
 *   - predictions: Fitted values
 *   - errors: LOOCV squared prediction errors
 *   - coefficients: [intercept, slope]
 *   - w: Input weights
 *
 * @throws std::runtime_error if sum of weights is not positive
 *
 * @note This implementation is optimized for degree 1 models by:
 *   - Avoiding matrix operations
 *   - Minimizing memory allocations
 *   - Using centered data for numerical stability
 *   - Computing statistics in single passes where possible
 */
ulm_results fit_linear_direct(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool y_binary,
    double epsilon = 1e-10) {

    int n_points = w.size();
    ulm_results results;
    results.w = w;
    results.coefficients.resize(2);  // [intercept, slope]

    // Calculate weighted means in one pass
    double total_weight = 0.0;
    double x_wmean = 0.0;
    double y_wmean = 0.0;

    for (int i = 0; i < n_points; ++i) {
        total_weight += w[i];
        x_wmean += w[i] * x[i];
        y_wmean += w[i] * y[i];
    }

    if (total_weight <= 0) {
        Rf_error("Sum of weights must be positive");
    }

    x_wmean /= total_weight;
    y_wmean /= total_weight;

    // Calculate slope using centered data
    double sum_wx_squared = 0.0;
    double wxy_sum = 0.0;

    for (int i = 0; i < n_points; ++i) {
        double x_centered = x[i] - x_wmean;
        sum_wx_squared += w[i] * x_centered * x_centered;
        wxy_sum += w[i] * x_centered * (y[i] - y_wmean);
    }

    // Calculate coefficients
    double slope = (sum_wx_squared > epsilon) ? wxy_sum / sum_wx_squared : 0.0;
    double intercept = y_wmean - slope * x_wmean;

    results.coefficients[0] = intercept;
    results.coefficients[1] = slope;

    // Calculate predictions and errors
    results.predictions.resize(n_points);
    results.errors.resize(n_points);

    for (int i = 0; i < n_points; i++) {
        // Calculate prediction
        results.predictions[i] = intercept + slope * x[i];
        if (y_binary) {
            results.predictions[i] = std::clamp(results.predictions[i], 0.0, 1.0);
        }

        // Calculate leverage
        double x_centered = x[i] - x_wmean;
        double h_i = w[i] * (1.0/total_weight + (x_centered * x_centered) / sum_wx_squared);

        // Calculate LOOCV Rf_error
        if (1.0 - h_i > epsilon) {
            double residual = (y[i] - results.predictions[i]) / (1.0 - h_i);
            results.errors[i] = residual * residual;
        } else {
            results.errors[i] = std::numeric_limits<double>::infinity();
        }
    }

    return results;
}


/**
 * @brief Fits a polynomial model using weighted least squares
 *
 * Implements polynomial model fitting using Eigen-based weighted least squares
 * with ridge regularization for stability. Particularly suited for
 * degree 2 (quadratic) models.
 *
 * @param x Pointer to predictor variable array
 * @param y Pointer to response variable array
 * @param w Vector of observation weights
 * @param degree Polynomial degree (1 or 2)
 * @param y_binary Boolean indicating if response is binary (0/1)
 * @param ridge_lambda Ridge regularization parameter (default: 1e-10)
 *
 * @return ulm_results structure containing:
 *   - predictions: Fitted values
 *   - errors: LOOCV squared prediction errors
 *   - coefficients: Model coefficients [intercept, linear term, quadratic term]
 *   - w: Input weights
 * 
 * @note This implementation:
 *   - Uses Eigen for efficient matrix operations
 *   - Adds ridge penalty for numerical stability
 *   - Handles polynomial terms systematically
 *   - Computes LOOCV errors using hat matrix
 */
// WLS implementation for degree 2 (quadratic) models using Eigen
ulm_results fit_polynomial_wls(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int degree,
    bool y_binary,
    double ridge_lambda = 1e-10) {

    int n = w.size();
    int p = degree + 1;  // number of parameters (including intercept)

    ulm_results results;
    results.w = w;
    results.coefficients.resize(p);

    // Set up design matrix
    Eigen::MatrixXd X(n, p);
    Eigen::VectorXd Y(n);
    Eigen::VectorXd W = Eigen::Map<const Eigen::VectorXd>(w.data(), n);

    // Fill design matrix with powers of x
    for (int i = 0; i < n; ++i) {
        X(i, 0) = 1.0;  // intercept
        double x_power = 1.0;
        for (int j = 1; j < p; ++j) {
            x_power *= x[i];
            X(i, j) = x_power;
        }
        Y(i) = y[i];
    }

    // Weight the matrices
    Eigen::MatrixXd X_w = W.asDiagonal() * X;
    Eigen::VectorXd Y_w = W.cwiseProduct(Y);

    // Add ridge penalty for stability
    Eigen::MatrixXd XtX = X_w.transpose() * X_w;
    XtX.diagonal() += Eigen::VectorXd::Constant(p, ridge_lambda);

    // Solve weighted least squares with LDLT decomposition
    Eigen::LDLT<Eigen::MatrixXd> solver(XtX);
    Eigen::VectorXd beta = solver.solve(X_w.transpose() * Y_w);

    // Copy coefficients
    for (int i = 0; i < p; ++i) {
        results.coefficients[i] = beta(i);
    }

    // Compute predictions
    results.predictions.resize(n);
    for (int i = 0; i < n; ++i) {
        double pred = beta(0);
        double x_power = 1.0;
        for (int j = 1; j < p; ++j) {
            x_power *= x[i];
            pred += beta(j) * x_power;
        }
        results.predictions[i] = y_binary ? std::clamp(pred, 0.0, 1.0) : pred;
    }

    // Compute LOOCV errors using hat matrix
    Eigen::MatrixXd H = X_w * solver.solve(X.transpose() * W.asDiagonal());

    results.errors.resize(n);
    for (int i = 0; i < n; ++i) {
        double h_i = H(i, i);
        if (h_i < 1.0 - 1e-10) {
            double residual = (y[i] - results.predictions[i]) / (1.0 - h_i);
            results.errors[i] = residual * residual;
        } else {
            results.errors[i] = std::numeric_limits<double>::infinity();
        }
    }

    return results;
}

/**
 * @brief Unified interface for fitting univariate polynomial models
 *
 * Wrapper function that automatically selects the most appropriate
 * implementation based on the polynomial degree:
 *   - For degree 1: Uses fast direct computation
 *   - For degree 2: Uses WLS with matrix operations
 *
 * @param x Pointer to predictor variable array
 * @param y Pointer to response variable array
 * @param w Vector of observation weights
 * @param degree Polynomial degree (1 or 2)
 * @param y_binary Boolean indicating if response is binary (0/1)
 * @param epsilon Small number for numerical stability in linear case (default: 1e-10)
 * @param ridge_lambda Ridge parameter for polynomial case (default: 1e-10)
 *
 * @return ulm_results structure containing model results
 *
 * @throws std::runtime_error if:
 *   - degree is not 1 or 2
 *   - sum of weights is not positive
 *
 * @note Automatically chooses between:
 *   - Direct computation for linear models (faster, less memory)
 *   - WLS matrix approach for quadratic models (more stable)
 */
ulm_results fit_ulm(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int degree,
    bool y_binary = false,
    double epsilon = 1e-10,
    double ridge_lambda = 1e-10) {

    if (degree < 1 || degree > 2) {
        Rf_error("Only degrees 1 and 2 are supported");
    }

    if (degree == 1) {
        return fit_linear_direct(x, y, w, y_binary, epsilon);
    } else {
        return fit_polynomial_wls(x, y, w, degree, y_binary, ridge_lambda);
    }
}



struct ortho_polynomial_t {  // More descriptive name
    std::vector<double> predictions;
    std::vector<double> errors;
    std::vector<double> coefficients;
    std::vector<double> w;

    // Store orthogonalization information
    struct orthogonal_basis_t {
        std::vector<double> norms;  // Normalization factors
        std::vector<std::vector<double>> proj_coeffs;  // Projection coefficients
    } basis;

    std::vector<double> predict(const std::vector<double>& x) const;
};

std::vector<double> ortho_polynomial_t::predict(const std::vector<double>& x) const {
    int n = x.size();
    int p = coefficients.size();
    std::vector<double> predictions(n);
    std::vector<std::vector<double>> basis_values(p, std::vector<double>(n));

    // Compute first basis (constant)
    for (int i = 0; i < n; ++i) {
        basis_values[0][i] = 1.0 / basis.norms[0];
    }

    // Compute subsequent bases using stored orthogonalization coefficients
    for (int j = 1; j < p; ++j) {
        for (int i = 0; i < n; ++i) {
            // Start with x * previous basis
            basis_values[j][i] = x[i] * basis_values[j-1][i];

            // Subtract projections
            for (int k = 0; k < j; ++k) {
                basis_values[j][i] -= basis.proj_coeffs[j][k] * basis_values[k][i];
            }

            // Normalize
            basis_values[j][i] /= basis.norms[j];
        }
    }

    // Compute predictions
    for (int i = 0; i < n; ++i) {
        double pred = 0.0;
        for (int j = 0; j < p; ++j) {
            pred += coefficients[j] * basis_values[j][i];
        }
        predictions[i] = pred;
    }

    return predictions;
}

ortho_polynomial_t fit_ortho_polynomial_wls(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int degree,
    bool y_binary,
    double ridge_lambda = 1e-10) {

    int n = w.size();
    int p = degree + 1;
    ortho_polynomial_t results;
    results.w = w;
    results.coefficients.resize(p);
    results.basis.norms.resize(p);
    results.basis.proj_coeffs.resize(p);

    Eigen::MatrixXd X(n, p);
    Eigen::VectorXd Y = Eigen::Map<const Eigen::VectorXd>(y, n);
    Eigen::VectorXd W = Eigen::Map<const Eigen::VectorXd>(w.data(), n);

    // First basis (constant)
    X.col(0).setOnes();
    results.basis.norms[0] = std::sqrt(W.sum());
    X.col(0) /= results.basis.norms[0];

    // For each subsequent polynomial, orthogonalize against previous ones
    for (int j = 1; j < p; ++j) {
        results.basis.proj_coeffs[j].resize(j);

        // Compute next polynomial (x * previous basis)
        Eigen::VectorXd new_col = X.col(j-1).cwiseProduct(Eigen::Map<const Eigen::VectorXd>(x, n));

        // Gram-Schmidt with respect to w
        for (int k = 0; k < j; ++k) {
            double numerator = (new_col.cwiseProduct(X.col(k))).dot(W);
            double denominator = X.col(k).squaredNorm() * W.sum();  // Optimized denominator
            results.basis.proj_coeffs[j][k] = numerator / denominator;
            new_col -= results.basis.proj_coeffs[j][k] * X.col(k);
        }

        // Normalize
        results.basis.norms[j] = std::sqrt(new_col.cwiseProduct(new_col).dot(W));
        X.col(j) = new_col / results.basis.norms[j];
    }

    // Weight the matrices
    Eigen::MatrixXd X_w = W.asDiagonal() * X;
    Eigen::VectorXd Y_w = W.cwiseProduct(Y);

    // Since X is orthogonal with respect to W, XtWX is diagonal
    // We can solve the system more efficiently
    Eigen::VectorXd XtWX_diag = X_w.cwiseProduct(X).colwise().sum();
    XtWX_diag.array() += ridge_lambda;

    // Solve system directly using diagonal structure
    Eigen::VectorXd beta = (X_w.transpose() * Y_w).cwiseQuotient(XtWX_diag);

    // Store coefficients
    results.coefficients = std::vector<double>(beta.data(), beta.data() + p);

    // Compute predictions
    results.predictions.resize(n);
    Eigen::VectorXd preds = X * beta;
    for (int i = 0; i < n; ++i) {
        results.predictions[i] = y_binary ? std::clamp(preds(i), 0.0, 1.0) : preds(i);
    }

    // Compute LOOCV errors using simplified hat matrix computation
    // Since X is orthogonal, H is simpler to compute
    results.errors.resize(n);
    for (int i = 0; i < n; ++i) {
        double h_i = 0.0;
        for (int j = 0; j < p; ++j) {
            h_i += w[i] * X(i, j) * X(i, j) / XtWX_diag(j);
        }

        if (h_i < 1.0 - 1e-10) {
            double residual = (y[i] - results.predictions[i]) / (1.0 - h_i);
            results.errors[i] = residual * residual;
        } else {
            results.errors[i] = std::numeric_limits<double>::infinity();
        }
    }

    return results;
}


/**
 * @brief Robust local linear model fitting using Cleveland's iterative reweighting
 *
 * @param x Array of x-coordinates
 * @param y Array of y-coordinates
 * @param w Initial weights for each data point
 * @param y_binary Whether y is binary (0/1) - affects prediction constraints
 * @param tolerance Convergence tolerance for linear model fitting
 * @param n_iter Number of robustness iterations (typically 1-3)
 * @param robust_scale Scale factor for residuals (Cleveland recommends 6.0)
 * @return ulm_t Fitted model with robust weights
 */
ulm_t cleveland_ulm(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool y_binary,
    double tolerance,
    int n_iter,
    double robust_scale) {

    // Initial non-robust fit
    ulm_t fit = ulm(x, y, w, y_binary, tolerance);

    // Get data size
    size_t n = w.size();

    // Current working weights (initialize to original weights)
    std::vector<double> current_weights = w;

    // Iterative robust fitting
    for (int iter = 0; iter < n_iter; ++iter) {
        // Compute residuals
        std::vector<double> residuals(n);
        for (size_t i = 0; i < n; ++i) {
            residuals[i] = std::abs(y[i] - fit.predictions[i]);
        }

        // Find median absolute residual for scaling
        if (residuals.empty()) continue;

        std::vector<double> abs_residuals = residuals;
        std::nth_element(abs_residuals.begin(),
                        abs_residuals.begin() + abs_residuals.size()/2,
                        abs_residuals.end());
        double median_abs_residual = abs_residuals[abs_residuals.size()/2];

        // Convert MAD to sigma estimate under Normality (avoid zero).
        const double mad_sigma = (median_abs_residual / 0.6745);

        // Avoid division by zero
        if (median_abs_residual < 1e-10) {
            break; // No need for further iterations, fit is already good
        }

        // Scale residuals
        // double scale = robust_scale * median_abs_residual;
        const double scale = robust_scale * (mad_sigma + tolerance);

        for (auto& r : residuals) {
            r /= scale;
        }

        // Apply bisquare weights
        std::vector<double> robust_weights(n);
        for (size_t i = 0; i < n; ++i) {
            double u = residuals[i];
            if (u >= 1.0) {
                robust_weights[i] = 0.0;
            } else {
                double tmp = 1.0 - u*u;
                robust_weights[i] = tmp*tmp;
            }
        }

        // Combine with original weights
        for (size_t i = 0; i < n; ++i) {
            current_weights[i] = w[i] * robust_weights[i];
        }

        // Normalize weights if they sum to very small value
        double weight_sum = std::accumulate(current_weights.begin(), current_weights.end(), 0.0);
        if (weight_sum < 1e-10) {
            // Revert to original weights if robust weights become too small
            current_weights = w;
        } else if (std::abs(weight_sum - 1.0) > 1e-6) {
            // Normalize weights to sum to 1
            for (auto& weight : current_weights) {
                weight /= weight_sum;
            }
        }

        // Refit model with new weights
        fit = ulm(x, y, current_weights, y_binary, tolerance);
    }

    // Return the robustly fitted model
    return fit;
}
