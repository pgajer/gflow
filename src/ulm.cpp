
#include <Eigen/Dense>
#include <vector>
#include <cmath>      // For fabs()
#include <algorithm>  // For std::find, std::clamp
#include <numeric>    // For std::accumulate

#include "ulm.hpp"
#include "error_utils.h"

/**
 * @brief Fits a weighted linear model and computes predictions with LOOCV errors
 *
 * This function fits a weighted linear regression model to one-dimensional data and computes
 * both predictions and their Leave-One-Out Cross-Validation (LOOCV) errors. For each point i,
 * the LOOCV squared error is computed as:
 * \f[
 *     \text{error}_i = \left(\frac{y_i - \hat{y}_i}{1 - h_i}\right)^2
 * \f]
 * where \f$h_i\f$ is the leverage statistic for weighted linear regression.
 *
 * For the 1D case with predictor x and weights w, the weighted leverage \f$h_i\f$ at point i is:
 * \f[
 *     h_i = w_i\left(\frac{1}{\sum w_j} + \frac{(x_i - \bar{x}_w)^2}{\sum_{j=1}^n w_j(x_j - \bar{x}_w)^2}\right)
 * \f]
 * where \f$\bar{x}_w = \frac{\sum w_j x_j}{\sum w_j}\f$ is the weighted mean of x.
 *
 * The fitted linear model makes predictions using:
 * \f[
 *     \hat{y}(x) = \bar{y}_w + \beta(x - \bar{x}_w)
 * \f]
 * where \f$\beta\f$ is the slope coefficient and \f$\bar{y}_w\f$ is the weighted mean of y.
 *
 * @param y Vector of response variables
 * @param x Vector of predictor variables
 * @param w Vector of weights for weighted least squares regression
 * @param x_min_index Index of smallest x value in the larger dataset
 * @param x_max_index Index of largest x value in the larger dataset
 * @param y_binary Whether y values should be clamped to [0,1] range (default: false)
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return ulm_t structure containing:
 *         Data components:
 *         - x: Original predictor values
 *         - w: Original weights
 *         - x_min_index: Index of smallest x value
 *         - x_max_index: Index of largest x value
 *
 *         Model components:
 *         - slope: Fitted slope coefficient β
 *         - y_wmean: Weighted mean of response variable
 *         - x_wmean: Weighted mean of predictor variable
 *
 *         Model evaluation:
 *         - predictions: Vector of fitted values at each x point
 *         - errors: Vector of LOOCV squared errors at each point
 *
 *         Methods:
 *         - predict(pt_index): Returns prediction at given index
 *         - predict_with_error(pt_index): Returns pair of prediction and LOOCV error
 *
 * @throws REPORT_ERROR if:
 *         - Sum of weights is not positive
 *         - Prediction is requested for index outside [x_min_index, x_max_index]
 *
 * @note For constant x values (weighted variance less than epsilon), returns weighted mean of y
 *       with slope set to 0
 * @note When leverage h_i is within epsilon of 1, the LOOCV error is set to infinity
 */
ulm_t ulm(const double* x,
          const double* y,
          const std::vector<double>& w,
          bool y_binary,
          double epsilon) {

    int n_points = w.size();

    ulm_t results;  // We'll populate this with more information

    // Create working copy of x for computations
    std::vector<double> x_working(x, x + n_points);

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        total_weight += weight;
    }
    if (total_weight <= 0) REPORT_ERROR("total_weight: %.2f   Sum of weights must be positive", total_weight);

    // Calculate weighted means - MODIFY THIS SECTION
    // Store these in the results structure
    results.x_wmean = 0.0;
    results.y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        results.x_wmean += w[i] * x_working[i] / total_weight;
        results.y_wmean += w[i] * y[i] / total_weight;
    }

    // Center x_working around weighted mean for leverage calculation
    double sum_wx_squared = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_working[i] -= results.x_wmean;  // Now x is centered around its weighted mean
        sum_wx_squared += w[i] * x_working[i] * x_working[i];
    }

    // Calculate slope if x has sufficient variation - MODIFY THIS SECTION
    results.slope = 0.0;
    if (sum_wx_squared > epsilon) {
        double wxy_sum = 0.0;
        for (int i = 0; i < n_points; ++i) {
            wxy_sum += w[i] * x_working[i] * y[i];
        }
        results.slope = wxy_sum / sum_wx_squared;
    }

    // Calculate predicted values - NO CHANGE NEEDED
    results.predictions.resize(n_points);
    for (int i = 0; i < n_points; i++) {
        results.predictions[i] = results.y_wmean + results.slope * (x[i] - results.x_wmean);
        if (y_binary) {
            results.predictions[i] = std::clamp(results.predictions[i], 0.0, 1.0);
        }
    }

    // Calculate LOOCV components - NO CHANGE NEEDED
    results.errors.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
        // Calculate weighted leverage h_i
        double h_i = w[i] * (1.0/total_weight + (x_working[i] * x_working[i]) / sum_wx_squared);
        // Calculate LOOCV prediction error
        if (std::abs(1.0 - h_i) > epsilon) {
            double residual = (y[i] - results.predictions[i]) / (1.0 - h_i);
            results.errors[i] = residual * residual;
        } else {
            // Handle case where leverage is close to 1
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
        throw std::runtime_error("Sum of weights must be positive");
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

        // Calculate LOOCV error
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
        throw std::runtime_error("Only degrees 1 and 2 are supported");
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

        // Avoid division by zero
        if (median_abs_residual < 1e-10) {
            break; // No need for further iterations, fit is already good
        }

        // Scale residuals
        double scale = robust_scale * median_abs_residual;
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
