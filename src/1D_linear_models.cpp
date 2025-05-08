#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <algorithm>  // for std::find, std::clamp
#include <vector>       // for std::vector
#include <utility>      // for std::pair
#include <cmath>        // for fabs()

#include "kernels.h"
#include "1D_linear_models.h"
#include "error_utils.h"

#if 0
/**
 * @brief Computes kernel-weighted linear regression for data points positioned relative to a reference point.
 *
 * @details This function performs weighted linear regression on data points arranged relative
 * to a reference point (at x=0). The data points are provided as separate left and right sequences,
 * representing points on either side of the reference. The regression incorporates both
 * user-provided weights and kernel-based weights that depend on the distance from the reference point.
 *
 * The algorithm follows these steps:
 * 1. Combines left and right sequences into unified arrays
 * 2. Normalizes distances for kernel weight calculation
 * 3. Applies kernel weights and combines them with user weights
 * 4. Performs weighted least squares regression
 *
 * For numerical stability, when the weighted sum of squared x-deviations is very small
 * (less than epsilon), the function returns the weighted mean of y as the intercept
 * and 0 as the slope.
 *
 * @param Ld Vector of distances for points left of the reference point
 * @param Rd Vector of distances for points right of the reference point
 * @param Ly Vector of y-values corresponding to left points
 * @param Ry Vector of y-values corresponding to right points
 * @param Lw Vector of weights for left points
 * @param Rw Vector of weights for right points
 * @param kernel_index Integer specifying the kernel function to use:
 *                     0 = Constant
 *                     1 = Epanechnikov
 *                     2 = Triangular
 *                     3 = Truncated exponential
 *                     4 = Normal
 * @param normalization_factor Factor used to normalize distances (default = 1.01)
 * @param epsilon Threshold for detecting numerically unstable cases (default = 1e-10)
 *
 * @return std::pair<double,double> First element is the y-intercept (predicted value at x=0),
 *                                 second element is the slope of the fitted line
 *
 * @throws Rf_error if input vectors on either side have different lengths
 *
 * @note The function is particularly optimized for predicting values at the reference point (x=0),
 *       which corresponds to the y-intercept of the fitted line.
 *
 * @warning The kernel_index parameter must be initialized with initialize_kernel() before calling
 *          this function.
 *
 * @pre All input vectors for each side (left/right) must have the same length
 * @pre kernel_index must be valid (0-4)
 * @pre normalization_factor should be positive
 * @pre epsilon should be positive and small
 *
 * @see initialize_kernel()
 * @see epanechnikov_kernel()
 * @see triangular_kernel()
 * @see tr_exponential_kernel()
 * @see normal_kernel()
 */
std::pair<double,double> linear_model_1d(
    const std::vector<double>& Ld,
    const std::vector<double>& Rd,
    const std::vector<double>& Ly,
    const std::vector<double>& Ry,
    const std::vector<double>& Lw,
    const std::vector<double>& Rw,
    int kernel_index,
    double normalization_factor,
    double epsilon) {

    // Check that all Left vectors have the same length
    int Lsize = Ld.size();
    if (Lsize != Ly.size() || Lsize != Lw.size()) {
        Rprintf("Ld.size: %d\n", Lsize);
        Rprintf("Ly.size: %d\n", (int)Ld.size());
        Rprintf("Lw.size: %d\n", (int)Lw.size());
        Rf_error("All left arrays have to have the same length");
    }

    // Check that all Right vectors have the same length
    int Rsize = Rd.size();
    if (Rsize != Ry.size() || Rsize != Rw.size()) {
        Rprintf("Rd.size: %d\n", Rsize);
        Rprintf("Ry.size: %d\n", (int)Rd.size());
        Rprintf("Rw.size: %d\n", (int)Rw.size());
        Rf_error("All right arrays have to have the same length");
    }

    // Defining dist, y and w vectors
    int n_points = Lsize + Rsize;
    std::vector<double> x = std::vector<double>(n_points);
    std::vector<double> y = std::vector<double>(n_points);
    std::vector<double> w = std::vector<double>(n_points);
    for (int i = 0; i < Lsize; ++i) {
        x[i] = Ld[i];
        y[i] = Ly[i];
        w[i] = Lw[i];
    }
    for (int i = 0; i < Rsize; ++i) {
        int shifted_i = i + Lsize;
        x[shifted_i] = Rd[i];
        y[shifted_i] = Ry[i];
        w[shifted_i] = Rw[i];
    }

    // Normalizing x's
    double max_x = 0.0;
    for (int i = 0; i < n_points; ++i) {
        if (x[i] > max_x)
            max_x = x[i];
    }
    if (max_x == 0) max_x = 1;  // Avoid division by zero
    max_x *= normalization_factor;
    for (int i = 0; i < n_points; ++i) {
        x[i] /= max_x;
    }

    // Creating kernel weights
    double scale = 1.0;
    initialize_kernel(kernel_index, scale);
    std::vector<double> kernel_weight = std::vector<double>(n_points);
    kernel_fn(x.data(), n_points, kernel_weight.data());

    // Multiplying weights by kernel weights <<--- NOTE that w's are adjusted by kernel weights !!!
    double total_weight = 0;
    for (int i = 0; i < n_points; ++i) {
        w[i] *= kernel_weight[i];
        total_weight += w[i];
    }

    // Normalizing the new weights so that their sum is 1
    for (int i = 0; i < n_points; ++i) {
        w[i] /= total_weight;
    }

    // x, y weighted means
    double x_wmean = 0;
    double y_wmean = 0;
    for (int i = 0; i < n_points; ++i) {
        x_wmean += w[i] * x[i];
        y_wmean += w[i] * y[i];
    }

    // slope
    double slope = 0;
    double xx_wmean = 0;
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w[i] * xx * xx; // <<--- correct formula
        slope += w[i] * xx * (y[i] - y_wmean);
    }

    if (fabs(xx_wmean) < epsilon) {
        return std::make_pair(y_wmean, 0.0);
    } else {
        slope /= xx_wmean;
    }

    // y-intercept
    double y_intercept = y_wmean - slope * x_wmean;

    return std::make_pair(y_intercept, slope);
}
#endif

/**
* @brief Predicts the value at the reference point (x=0) using kernel-weighted linear regression.
*
* @details Performs weighted linear regression on data points arranged relative to a reference point,
* specifically optimized for predicting the value at x=0. The data points are provided as separate
* left and right sequences from the reference point. The regression incorporates both user-provided
* weights and kernel-based weights that depend on the distance from the reference point.
*
* The algorithm steps are:
* 1. Validates input data dimensions
* 2. Normalizes distances for kernel weight calculation
* 3. Computes and combines kernel weights with user weights
* 4. Calculates weighted means and regression coefficients
* 5. Returns the predicted value at x=0 (y-intercept)
*
* For numerical stability, when the weighted sum of squared x-deviations is very small
* (less than epsilon), the function returns the weighted mean of y values.
*
* @param Ld Vector of distances for points left of the reference point
* @param Rd Vector of distances for points right of the reference point
* @param Ly Vector of y-values corresponding to left points
* @param Ry Vector of y-values corresponding to right points
* @param Lw Vector of weights for left points
* @param Rw Vector of weights for right points
* @param kernel_index Integer specifying the kernel function:
*                     0 = Constant
*                     1 = Epanechnikov
*                     2 = Triangular
*                     3 = Truncated exponential
*                     4 = Normal
* @param normalization_factor Factor for distance normalization (default = 1.01)
* @param epsilon Numerical stability threshold (default = 1e-8)
*
* @return double Predicted value at the reference point (x=0)
*
* @throws Rf_error if input vectors on either side have inconsistent lengths
*
* @note This function is optimized for predicting at x=0 and is more efficient than
*       calculating both slope and intercept when only the prediction at the reference
*       point is needed.
*
* @note When xx_wmean < epsilon, the function returns the weighted mean of y values,
*       which is appropriate when all influential points are very close to x=0.
*
* @warning The kernel_index must be initialized with initialize_kernel() before calling
*          this function.
*
* @pre All input vectors for each side (left/right) must have the same length
* @pre kernel_index must be valid (0-4)
* @pre normalization_factor should be positive
* @pre epsilon should be positive and small
*
* Example:
* @code
* std::vector<double> Ld = {1.0, 2.0};     // Left distances
* std::vector<double> Rd = {1.5, 2.5};     // Right distances
* std::vector<double> Ly = {1.1, 1.2};     // Left y-values
* std::vector<double> Ry = {0.9, 0.8};     // Right y-values
* std::vector<double> Lw = {1.0, 1.0};     // Left weights
* std::vector<double> Rw = {1.0, 1.0};     // Right weights
* int kernel = 1;                          // Epanechnikov kernel
* double prediction = predict_1d_linear_model(Ld, Rd, Ly, Ry, Lw, Rw, kernel);
* @endcode
*
* @see initialize_kernel()
* @see epanechnikov_kernel()
* @see triangular_kernel()
* @see tr_exponential_kernel()
* @see normal_kernel()
*/
double predict_linear_model_1d(
    const std::vector<double>& Ld,
    const std::vector<double>& Rd,
    const std::vector<double>& Ly,
    const std::vector<double>& Ry,
    const std::vector<double>& Lw,
    const std::vector<double>& Rw,
    int kernel_index,
    double normalization_factor = 1.01,
    double epsilon = 1e-8) {

    // Check that all Left vectors have the same length
    int Lsize = Ld.size();
    if (Lsize != Ly.size() || Lsize != Lw.size()) {
        Rprintf("Ld.size: %d\n", Lsize);
        Rprintf("Ly.size: %d\n", (int)Ld.size());
        Rprintf("Lw.size: %d\n", (int)Lw.size());
        Rf_error("All left arrays have to have the same length");
    }

    // Check that all Right vectors have the same length
    int Rsize = Rd.size();
    if (Rsize != Ry.size() || Rsize != Rw.size()) {
        Rprintf("Rd.size: %d\n", Rsize);
        Rprintf("Ry.size: %d\n", (int)Rd.size());
        Rprintf("Rw.size: %d\n", (int)Rw.size());
        Rf_error("All right arrays have to have the same length");
    }

    int n_points = Lsize + Rsize;

    // We could potentially avoid creating the y vector and combine steps
    std::vector<double> x(n_points);
    std::vector<double> w(n_points);

    // Combine x and w arrays and compute maximum x
    double max_x = 0.0;
    for (int i = 0; i < Lsize; ++i) {
        x[i] = Ld[i];
        w[i] = Lw[i];
        max_x = std::max(max_x, std::abs(x[i]));
    }
    for (int i = 0; i < Rsize; ++i) {
        int shifted_i = i + Lsize;
        x[shifted_i] = Rd[i];
        w[shifted_i] = Rw[i];
        max_x = std::max(max_x, std::abs(x[shifted_i]));
    }

    // Normalize x values
    if (max_x == 0) max_x = 1;
    max_x *= normalization_factor;
    for (double& xi : x) xi /= max_x;

    // Apply kernel weights
    double scale = 1.0;
    initialize_kernel(kernel_index, scale);
    std::vector<double> kernel_weight(n_points);
    kernel_fn(x.data(), n_points, kernel_weight.data());

    // Combine all weights and compute weighted means in one pass
    double total_weight = 0.0;
    double x_wmean = 0.0;
    double y_wmean = 0.0;

    for (int i = 0; i < Lsize; ++i) {
        w[i] *= kernel_weight[i];
        total_weight += w[i];
        x_wmean += w[i] * x[i];
        y_wmean += w[i] * Ly[i];
    }
    for (int i = 0; i < Rsize; ++i) {
        int shifted_i = i + Lsize;
        w[shifted_i] *= kernel_weight[shifted_i];
        total_weight += w[shifted_i];
        x_wmean += w[shifted_i] * x[shifted_i];
        y_wmean += w[shifted_i] * Ry[i];
    }

    // Normalize means
    x_wmean /= total_weight;
    y_wmean /= total_weight;

    // Compute slope components
    double xx_wmean = 0.0;
    double xy_wmean = 0.0;

    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w[i] * xx * xx;
        xy_wmean += w[i] * xx * ((i < Lsize) ? Ly[i] : Ry[i-Lsize]);
    }

    if (fabs(xx_wmean) < epsilon) {
        return y_wmean;
    }

    double slope = xy_wmean / xx_wmean;
    return y_wmean - slope * x_wmean;
}


/**
 * @brief Predicts the value of a weighted linear regression model at a reference point
 *
 * This function fits a weighted linear regression model to the provided data points
 * and returns the predicted value at the reference point. The prediction is computed
 * using the following steps:
 * 1. Centers the x-coordinates around the reference point
 * 2. Normalizes the weights to sum to 1
 * 3. Computes weighted means of x and y
 * 4. Calculates the weighted regression slope
 * 5. Returns the predicted y-value at x = 0 (the centered reference point)
 *
 * @param y Vector of dependent variable values
 * @param x Vector of independent variable values (will be modified during computation)
 * @param w Vector of weights for each observation
 * @param ref_index Index of the reference point in the input vectors
 * @param epsilon Small value to handle numerical stability when xx_wmean is close to zero
 *
 * @pre All input vectors (y, x, w) must have the same length
 * @pre ref_index must be within bounds: 0 <= ref_index < x.size()
 * @pre All weights must be non-negative
 * @pre The sum of weights must be positive
 *
 * @throws Rf_error if input vectors have different lengths
 * @throws Rf_error if ref_index is out of bounds
 * @throws Rf_error if weights are invalid
 *
 * @note This function modifies the input vector x by centering it around the reference point
 *
 * @return Predicted y-value at the reference point
 *
 * Example usage:
 * @code
 * std::vector<double> y = {1.0, 2.0, 3.0};
 * std::vector<double> x = {0.0, 1.0, 2.0};
 * std::vector<double> w = {1.0, 1.0, 1.0};
 * int ref_idx = 0;
 * double prediction = predict_lm_1d(y, x, w, ref_idx);
 * @endcode
 */
double predict_lm_1d(const std::vector<double>& y,
                    std::vector<double>& x,
                    const std::vector<double>& w,
                    int ref_index,
                    double epsilon = 1e-8) {
    // Check that all vectors have the same length
    int n_points = x.size();
    if (n_points != y.size() || n_points != w.size()) {
        Rprintf("x.size: %d\n", n_points);
        Rprintf("y.size: %d\n", (int)y.size());
        Rprintf("w.size: %d\n", (int)w.size());
        Rf_error("All vectors have to have the same length");
    }

    // Validate ref_index
    if (ref_index < 0 || ref_index >= n_points) {
        Rf_error("Reference index out of bounds");
    }

    // Validate weights
    double total_weight = 0.0;
    for (int i = 0; i < n_points; ++i) {
        if (w[i] < 0) {
            Rf_error("Negative weights are not allowed");
        }
        total_weight += w[i];
    }
    if (total_weight <= 0) {
        Rf_error("Sum of weights must be positive");
    }

    // Create normalized weights
    std::vector<double> w_normalized(w);
    for (int i = 0; i < n_points; ++i) {
        w_normalized[i] /= total_weight;
    }

    // Center x values around reference point
    double x_ref = x[ref_index];
    for (int i = 0; i < n_points; ++i) {
        x[i] -= x_ref;
    }

    // Calculate weighted means
    double x_wmean = 0.0;
    double y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_wmean += w_normalized[i] * x[i];
        y_wmean += w_normalized[i] * y[i];
    }

    // Compute slope components
    double xx_wmean = 0.0;
    double xy_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w_normalized[i] * xx * xx;
        xy_wmean += w_normalized[i] * xx * y[i];
    }

    // Handle case where x values are effectively constant
    if (fabs(xx_wmean) <= epsilon) {
        return y_wmean;
    }

    // Calculate and return prediction at reference point
    double slope = xy_wmean / xx_wmean;

    return y_wmean - slope * x_wmean;
}




/**
 * @brief Predicts the value and calculates the derivative of a weighted linear regression model at a reference point
 *
 * This function fits a weighted linear regression model to the provided data points and returns both
 * the predicted value and the partial derivative ∂μᵢ/∂yᵢ at the reference point.
 *
 * @details
 * For local linear regression, the estimate μᵢ at point i is given by:
 * μᵢ = β₀ + β₁(xᵢ - x_ref)
 *
 * Where β₀, β₁ are obtained by minimizing:
 * ∑ⱼ wⱼ(yⱼ - β₀ - β₁(xⱼ - x_ref))²
 *
 * The solution can be written in matrix form as:
 * μᵢ = e₁ᵀ(XᵀWX)⁻¹XᵀWy
 *
 * where:
 * - X is the design matrix with rows (1, xⱼ - x_ref)
 * - W is the diagonal matrix of normalized weights
 * - e₁ᵀ = (1,0)
 *
 * The partial derivative ∂μᵢ/∂yᵢ is the i-th diagonal element of the hat matrix:
 * H = X(XᵀWX)⁻¹XᵀW
 *
 * For the simple case of local linear regression, this reduces to:
 * ∂μᵢ/∂yᵢ = wᵢ(1 - (xᵢ - x̄_w)²/∑ⱼwⱼ(xⱼ - x̄_w)²)
 *
 * where:
 * - x̄_w is the weighted mean of x values
 * - wᵢ is the normalized weight at point i
 *
 * When x values are effectively constant (∑ⱼwⱼ(xⱼ - x̄_w)² ≈ 0),
 * the model reduces to weighted averaging and:
 * ∂μᵢ/∂yᵢ = wᵢ
 *
 * @param y Vector of dependent variable values
 * @param x Vector of independent variable values (will be modified during computation)
 * @param w Vector of weights for each observation
 * @param ref_index Index of the reference point in the input vectors
 * @param epsilon Small value to handle numerical stability when variance in x is close to zero
 *
 * @return PredictionResult struct containing:
 *         - prediction: Predicted y-value at the reference point
 *         - derivative: Value of ∂μᵢ/∂yᵢ at the reference point
 *
 * @pre All input vectors (y, x, w) must have the same length
 * @pre ref_index must be within bounds: 0 <= ref_index < x.size()
 * @pre All weights must be non-negative
 * @pre The sum of weights must be positive
 *
 * @throws Rf_error if input vectors have different lengths
 * @throws Rf_error if ref_index is out of bounds
 * @throws Rf_error if weights are invalid
 *
 * @note This function modifies the input vector x by centering it around the reference point
 *
 * Example usage:
 * @code
 * std::vector<double> y = {1.0, 2.0, 3.0};
 * std::vector<double> x = {0.0, 1.0, 2.0};
 * std::vector<double> w = {1.0, 1.0, 1.0};
 * int ref_idx = 0;
 * auto result = predict_lm_1d(y, x, w, ref_idx);
 * double prediction = result.prediction;  // predicted y-value
 * double derivative = result.derivative;  // ∂μᵢ/∂yᵢ at reference point
 * @endcode
 */
lm_dmudy_t predict_lm_dmudy(const std::vector<double>& y,
                            std::vector<double>& x,
                            const std::vector<double>& w,
                            int ref_index,
                            double epsilon = 1e-8) {
    // Check that all vectors have the same length
    int n_points = x.size();
    if (n_points != y.size() || n_points != w.size()) {
        Rprintf("x.size: %d\n", n_points);
        Rprintf("y.size: %d\n", (int)y.size());
        Rprintf("w.size: %d\n", (int)w.size());
        Rf_error("All vectors have to have the same length");
    }

    // Validate ref_index
    if (ref_index < 0 || ref_index >= n_points) {
        Rf_error("Reference index out of bounds");
    }

    // Validate weights
    double total_weight = 0.0;
    for (int i = 0; i < n_points; ++i) {
        if (w[i] < 0) {
            Rf_error("Negative weights are not allowed");
        }
        total_weight += w[i];
    }
    if (total_weight <= 0) {
        Rf_error("Sum of weights must be positive");
    }

    // Create normalized weights
    std::vector<double> w_normalized(w);
    for (int i = 0; i < n_points; ++i) {
        w_normalized[i] /= total_weight;
    }

    // Center x values around reference point
    double x_ref = x[ref_index];
    for (int i = 0; i < n_points; ++i) {
        x[i] -= x_ref;
    }

    // Calculate weighted means
    double x_wmean = 0.0;
    double y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_wmean += w_normalized[i] * x[i];
        y_wmean += w_normalized[i] * y[i];
    }

    // Compute slope components
    double xx_wmean = 0.0;
    double xy_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w_normalized[i] * xx * xx;
    }

    // Handle case where x values are effectively constant
    if (fabs(xx_wmean) <= epsilon) {
        // In this case, the estimate is just a weighted average
        return {y_wmean, w_normalized[ref_index]};
    }

    // Calculate the derivative ∂μᵢ/∂yᵢ
    // For local linear regression, this is the i-th element of the hat matrix
    double x_centered = x[ref_index] - x_wmean;
    double derivative = w_normalized[ref_index] * (1.0 - (x_centered * x_centered) / xx_wmean);

    // Calculate prediction
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xy_wmean += w_normalized[i] * xx * y[i];
    }
    double slope = xy_wmean / xx_wmean;
    double prediction = y_wmean - slope * x_wmean;

    return {prediction, derivative};
}



/**
 * @brief Fits a weighted local linear model and computes various goodness-of-fit metrics
 *
 * This function fits a weighted linear regression model to one-dimensional data and
 * returns both the predicted value at a reference point and several measures of fit quality.
 * The function centers the x-values around the reference point before fitting.
 *
 * @param y Vector of response variables
 * @param x Vector of predictor variables (modified during execution - centered around reference point)
 * @param w Vector of weights for weighted least squares regression
 * @param ref_index Index of the reference point in x and y vectors
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return LMFitResults structure containing:
 *         - predicted_value: Fitted value at the reference point
 *         - r_squared: Coefficient of determination (R²), ranges from 0 to 1
 *         - mae: Mean Absolute Error of residuals, weighted by w
 *         - rmse: Root Mean Square Error of residuals, weighted by w
 *         - aic: Akaike Information Criterion, using 2 parameters (slope and intercept)
 *
 * @throws Rf_error if:
 *         - Input vectors have different lengths
 *         - ref_index is out of bounds
 *         - Any weight is negative
 *         - Sum of weights is not positive
 *
 * @note The x vector is modified during execution (centered around reference point)
 * @note Weights are normalized internally to sum to 1
 * @note For constant x values (variance less than epsilon), returns mean of y
 *
 * @see LMFitResults for the structure definition
 * @see predict_lm_1d for the basic version without goodness-of-fit metrics
 *
 * @example
 * std::vector<double> y = {1.0, 2.0, 3.0, 4.0};
 * std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
 * std::vector<double> w = {1.0, 1.0, 1.0, 1.0};
 * auto results = predict_lm_1d_with_metrics(y, x, w, 1);
 * // results contains predicted value and goodness-of-fit metrics
 */
lm_fit_t predict_lm_1d_with_metrics(
    const std::vector<double>& y,
    std::vector<double>& x,
    const std::vector<double>& w,
    int ref_index,
    double epsilon = 1e-8) {

    lm_fit_t results;
    int n_points = x.size();

    // Original prediction logic
    double total_weight = 0.0;
    for (const auto& weight : w) total_weight += weight;

    std::vector<double> w_normalized(w);
    for (auto& weight : w_normalized) weight /= total_weight;

    double x_ref = x[ref_index];
    for (auto& xi : x) xi -= x_ref;

    // Calculate weighted means
    double x_wmean = 0.0, y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_wmean += w_normalized[i] * x[i];
        y_wmean += w_normalized[i] * y[i];
    }

    // Compute slope components
    double xx_wmean = 0.0, xy_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w_normalized[i] * xx * xx;
        xy_wmean += w_normalized[i] * xx * y[i];
    }

    // Calculate prediction and slope
    double slope = (fabs(xx_wmean) <= epsilon) ? 0.0 : xy_wmean / xx_wmean;
    results.predicted_value = y_wmean - slope * x_wmean;

    // Calculate residuals and goodness-of-fit measures
    double ss_res = 0.0;  // Sum of squared residuals
    double ss_tot = 0.0;  // Total sum of squares
    double sum_abs_res = 0.0;  // Sum of absolute residuals

    for (int i = 0; i < n_points; ++i) {
        double predicted = y_wmean + slope * (x[i] - x_wmean);
        double residual = y[i] - predicted;

        ss_res += w_normalized[i] * residual * residual;
        ss_tot += w_normalized[i] * (y[i] - y_wmean) * (y[i] - y_wmean);
        sum_abs_res += w_normalized[i] * std::abs(residual);
    }

    // Calculate R-squared
    results.r_squared = (ss_tot > epsilon) ? 1.0 - (ss_res / ss_tot) : 0.0;

    // Calculate MAE (Mean Absolute Error)
    results.mae = sum_abs_res;

    // Calculate RMSE (Root Mean Square Error)
    results.rmse = std::sqrt(ss_res);

    // Calculate AIC (Akaike Information Criterion)
    // k = 2 parameters (slope and intercept)
    results.aic = n_points * std::log(ss_res / n_points) + 2 * 2;

    return results;
}

/**
 * @brief Fits a weighted local linear model and computes Mean Absolute Error
 *
 * This function fits a weighted linear regression model to one-dimensional data and
 * returns both the predicted value at a reference point and the Mean Absolute Error (MAE)
 * as a measure of fit quality. The function centers the x-values around the reference point
 * before fitting.
 *
 * The MAE is calculated as the weighted average of absolute residuals using the normalized weights:
 * MAE = sum(w_i * |y_i - ŷ_i|), where w_i are normalized weights and ŷ_i are predicted values.
 *
 * @param y Vector of response variables
 * @param x Vector of predictor variables (modified during execution - centered around reference point)
 * @param w Vector of weights for weighted least squares regression
 * @param ref_index Index of the reference point in x and y vectors
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return lm_mae_t structure containing:
 *         - predicted_value: Fitted value at the reference point
 *         - mae: Mean Absolute Error of residuals, weighted by normalized weights
 *
 * @throws Rf_error if:
 *         - Input vectors have different lengths
 *         - ref_index is out of bounds
 *         - Any weight is negative
 *         - Sum of weights is not positive
 *
 * @note The x vector is modified during execution (centered around reference point)
 * @note Weights are normalized internally to sum to 1
 * @note For constant x values (variance less than epsilon), returns mean of y
 */
lm_mae_t predict_lm_1d_with_mae(const std::vector<double>& y,
                                std::vector<double>& x,
                                const std::vector<double>& w,
                                int ref_index,
                                double epsilon = 1e-8) {

    lm_mae_t results;
    int n_points = x.size();

    // Input validation
    if (n_points != y.size() || n_points != w.size()) {
        Rf_error("All vectors must have the same length");
    }
    if (ref_index < 0 || ref_index >= n_points) {
        Rf_error("Reference index out of bounds");
    }

    // Original prediction logic
    double total_weight = 0.0;
    for (const auto& weight : w) {
        if (weight < 0) Rf_error("Negative weights are not allowed");
        total_weight += weight;
    }
    if (total_weight <= 0) Rf_error("Sum of weights must be positive");

    std::vector<double> w_normalized(w);
    for (auto& weight : w_normalized) weight /= total_weight;

    double x_ref = x[ref_index];
    for (auto& xi : x) xi -= x_ref;

    // Calculate weighted means
    double x_wmean = 0.0, y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_wmean += w_normalized[i] * x[i];
        y_wmean += w_normalized[i] * y[i];
    }

    // Compute slope components
    double xx_wmean = 0.0, xy_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double xx = x[i] - x_wmean;
        xx_wmean += w_normalized[i] * xx * xx;
        xy_wmean += w_normalized[i] * xx * y[i];
    }

    // Calculate prediction and slope
    double slope = (fabs(xx_wmean) <= epsilon) ? 0.0 : xy_wmean / xx_wmean;
    results.predicted_value = y_wmean - slope * x_wmean;

    // Calculate MAE using normalized weights
    double sum_abs_res = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double predicted = y_wmean + slope * (x[i] - x_wmean);
        double residual = y[i] - predicted;
        sum_abs_res += w_normalized[i] * std::abs(residual);
    }
    results.mae = sum_abs_res;

    return results;
}

/**
 * @brief Fits a weighted local linear model and computes model components and LOOCV errors
 *
 * This function fits a weighted linear regression model to one-dimensional data and returns
 * both the fitted model components and LOOCV errors. The model uses the efficient
 * formula for linear model LOOCV:
 * \f[
 *     CV_{(n)} = \frac{1}{n} \sum_{i = 1}^{n} \left(\frac{y_{i} - \hat{y}_{i}}{1 - h_{i}} \right)^{2}
 * \f]
 * where \f$h_i\f$ is the leverage statistic for weighted linear regression.
 *
 * For the 1D case with predictor x and weights w, the weighted leverage \f$h_i\f$ at point i is:
 * \f[
 *     h_i = w_i\left(\frac{1}{W} + \frac{(x_i - \bar{x}_w)^2}{\sum_{j=1}^n w_j(x_j - \bar{x}_w)^2}\right)
 * \f]
 * where:
 * - \f$W = \sum_{j=1}^n w_j\f$ is the sum of weights
 * - \f$\bar{x}_w = \frac{1}{W}\sum_{j=1}^n w_j x_j\f$ is the weighted mean of x
 *
 * The fitted linear model has the form:
 * \f[
 *     \hat{y}(x) = \bar{y}_w + \beta(x - x_{ref} - \bar{x}_w)
 * \f]
 * where \f$\beta\f$ is the slope coefficient and \f$x_{ref}\f$ is the reference point value and \f$\bar{x}_w\f$, \f$\bar{y}_w\f$ are weighted means of \f$x\f$ and \f$y\f$ values.
 *
 * @param y Vector of response variables
 * @param x Vector of predictor variables (modified during execution - centered around reference point)
 * @param w Vector of weights for weighted least squares regression
 * @param vertex_indices Vector of indices corresponding to vertices in the training set
 * @param ref_index Index of the reference point in x and y vectors. If equal to -1, then predicted_value and loocv_at_ref_vertex are set to std::numeric_limits<double>::quiet_NaN()
 * @param epsilon Small positive number for numerical stability (default: 1e-8)
 *
 * @return lm_loocv_t structure containing:
 *         Model components:
 *         - slope: Fitted slope coefficient β
 *         - y_wmean: Weighted mean of response variable
 *         - x_wmean: Weighted mean of centered predictor variable
 *         - x_ref: Reference point value
 *         - vertex_indices: Indices of vertices in the training set
 *         - x_values: Original x values for these vertices
 *         - w_values: Weight values for the x_values
 *
 *         Model evaluation:
 *         - predicted_value: Fitted value at the reference point
 *         - loocv: Leave-One-Out Cross-Validation Mean Squared Error across all points
 *         - loocv_at_ref_vertex: LOOCV squared error at the reference point
 *
 *         Methods:
 *         - predict(x): Computes model prediction at any input value x
 *
 * @throws Rf_error if:
 *         - Input vectors have different lengths
 *         - ref_index is out of bounds
 *         - Any weight is negative
 *         - Sum of weights is not positive
 *
 * @note The x vector is modified during execution (centered first around reference point,
 *       then around weighted mean)
 * @note For constant x values (weighted variance less than epsilon), returns weighted mean of y
 *       with slope set to 0
 */
lm_loocv_t predict_lm_1d_loocv(const std::vector<double>& y,
                               const std::vector<double>& x,
                               const std::vector<double>& w,
                               const std::vector<int>& vertex_indices,
                               int ref_index,
                               bool y_binary = false,
                               double epsilon = 1e-8) {

    int n_points = x.size();
    if (ref_index < 0 || ref_index >= n_points) {
        Rf_error("Reference index out of bounds");
    }

    // Store vertex indices, x values, and weights
    lm_loocv_t results;
    results.vertex_indices = vertex_indices;
    results.x_values = x;
    results.w_values = w;  // Store the weight values
    results.ref_index = ref_index;

    #define DEBUG__predict_lm_1d_loocv 0
    #if DEBUG__predict_lm_1d_loocv
    Rprintf("\nIn predict_lm_1d_loocv()\n");
    Rprintf("n_points: %d\n", n_points);
    Rprintf("ref_index: %d\n", ref_index);
    print_vect(vertex_indices, "vertex_indices");
    print_vect(results.vertex_indices, "results.vertex_indices");
    print_vect(results.x_values, "results.x_values");
    print_vect(results.w_values, "results.w_values");
    Rprintf("results.ref_index: %d\n", results.ref_index);
    #endif

    // Create working copy of x for computations
    std::vector<double> x_working = x;

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        //if (weight < 0) Rf_error("Negative weights are not allowed");
        total_weight += weight;
    }
    if (total_weight <= 0) Rf_error("Sum of weights must be positive");

    // Store reference point
    results.x_ref = x[ref_index];

    // Center x_working around reference point
    for (auto& xi : x_working) xi -= results.x_ref;

    // Calculate weighted means
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

    // Calculate slope if x has sufficient variation
    results.slope = 0.0;
    if (sum_wx_squared > epsilon) {
        double wxy_sum = 0.0;
        for (int i = 0; i < n_points; ++i) {
            wxy_sum += w[i] * x_working[i] * y[i];
        }
        results.slope = wxy_sum / sum_wx_squared;
    }

    #if DEBUG__predict_lm_1d_loocv
    Rprintf("Just before results.predicted_value = results.predict(vertex_indices[ref_index])\n");
    print_vect(vertex_indices, "vertex_indices");
    Rprintf("ref_index: %d\n", ref_index);
    Rprintf("vertex_indices[ref_index]: %d\n", vertex_indices[ref_index]);
    #endif

    // Calculate predicted value at reference point
    results.predicted_value = results.predict(vertex_indices[ref_index]);

    if (y_binary) {
        results.predicted_value = std::clamp(results.predicted_value, 0.0, 1.0);
    }


    // Calculate LOOCV components
    double loocv_sum = 0.0;
    for (int i = 0; i < n_points; ++i) {
        // Calculate weighted leverage h_i
        double h_i = w[i] * (1.0/total_weight + (x_working[i] * x_working[i]) / sum_wx_squared);

        // Calculate fitted value using the predict method

        #if DEBUG__predict_lm_1d_loocv
        Rprintf("Just before y_hat = results.predict(vertex_indices[i])\n");
        print_vect(vertex_indices, "vertex_indices");
        Rprintf("i: %d\n", i);
        Rprintf("vertex_indices[i]: %d\n", vertex_indices[i]);
        #endif

        double y_hat = results.predict(vertex_indices[i]);
        if (y_binary) {
            y_hat = std::clamp(y_hat, 0.0, 1.0);
        }

        // Calculate LOOCV prediction error
        double residual = (y[i] - y_hat) / (1.0 - h_i);
        double squared_error = residual * residual;

        // Add to overall LOOCV sum
        loocv_sum += squared_error;

        // Store specific error for reference vertex
        if (i == ref_index) {
            results.loocv_at_ref_vertex = squared_error;
        }
    }

    // LOOCV MSE is the average of squared prediction errors
    results.loocv = loocv_sum / n_points;

    return results;
}


/**
 * @brief Computes linear models for multiple reference points along a path simultaneously
 *
 * @details This function is an optimized version of predict_lm_1d_loocv that computes
 * linear models for all positions in a path at once. It avoids redundant calculations
 * by processing all models for a path in a single pass. Each model uses a different
 * reference point and corresponding weight vector from w_list.
 *
 * @param y Vector of response variables
 * @param x Vector of predictor variables (distances along path)
 * @param vertex_indices Vector of vertex indices corresponding to path positions
 * @param w_list List of weight vectors, one for each reference point position
 * @param epsilon Small number for numerical stability (default: 1e-8)
 *
 * @return Vector of lm_loocv_t objects, one for each reference point
 *
 * @pre All input vectors (y, x, vertex_indices) must have the same length
 * @pre w_list must have the same length as other input vectors
 * @pre Each vector in w_list must sum to a positive value
 * @pre All vectors in w_list must have the same length as other input vectors
 *
 * @note This function assumes parameter validation is done by the calling function
 */
std::vector<lm_loocv_t> predict_lms_1d_loocv(const std::vector<double>& y,
                                             const std::vector<double>& x,
                                             const std::vector<int>& vertex_indices,
                                             const std::vector<std::vector<double>>& w_list,
                                             bool y_binary,
                                             double epsilon = 1e-8) {

    int n_points = x.size();
    if (n_points != w_list.size()) { // It is assumed that w_list has as many elements as the length of x (= length(y) = length(vertex_indices)) and so ref_index is any number between 0 and (n_points - 1)
        REPORT_ERROR("ref_index_list.size(): %d\tw_list.size(): %d\tref_index_list and w_list sizes have to be the same.", n_points, (int)w_list.size());
    }

    // Store vertex indices, x values, and weights
    std::vector<lm_loocv_t> results(n_points);
    std::vector<double> x_working(n_points);
    std::vector<double> w(n_points);
    for (int ref_index = 0; ref_index < n_points; ++ref_index) {

        results[ref_index].vertex_indices = vertex_indices;
        results[ref_index].x_values = x;
        results[ref_index].ref_index = ref_index;

        // Create working copy of x for computations
        x_working = x;

        // Center x_working around reference point
        double x_ref = x[ref_index];
        for (auto& xi : x_working) xi -= x_ref;
        results[ref_index].x_ref = x_ref;

        // Weight validation is done in the parent function
        double total_weight = 0.0;
        w = std::move(w_list[ref_index]);
        results[ref_index].w_values = w;
        for (const auto& weight : w) {
            total_weight += weight;
        }

        // Calculate weighted means
        results[ref_index].x_wmean = 0.0;
        results[ref_index].y_wmean = 0.0;
        for (int i = 0; i < n_points; ++i) {
            results[ref_index].x_wmean += w[i] * x_working[i] / total_weight;
            results[ref_index].y_wmean += w[i] * y[i] / total_weight;
        }

        // Center x_working around weighted mean for leverage calculation
        double sum_wx_squared = 0.0;
        for (int i = 0; i < n_points; ++i) {
            x_working[i] -= results[ref_index].x_wmean;  // Now x is centered around its weighted mean
            sum_wx_squared += w[i] * x_working[i] * x_working[i];
        }

        // Calculate slope if x has sufficient variation
        results[ref_index].slope = 0.0;
        if (sum_wx_squared > epsilon) {
            double wxy_sum = 0.0;
            for (int i = 0; i < n_points; ++i) {
                wxy_sum += w[i] * x_working[i] * y[i];
            }
            results[ref_index].slope = wxy_sum / sum_wx_squared;
        }

        // Calculate predicted value at reference point
        results[ref_index].predicted_value = results[ref_index].predict(vertex_indices[ref_index]);

        if (y_binary) {
            results[ref_index].predicted_value = std::clamp(results[ref_index].predicted_value, 0.0, 1.0);
        }

        // Calculate LOOCV components
        double loocv_sum = 0.0;
        for (int i = 0; i < n_points; ++i) {
            // Calculate weighted leverage h_i
            double h_i = w[i] * (1.0/total_weight + (x_working[i] * x_working[i]) / sum_wx_squared);

            // Calculate fitted value using the predict method
            double y_hat = results[ref_index].predict(vertex_indices[i]);
            if (y_binary) {
                y_hat = std::clamp(y_hat, 0.0, 1.0);
            }

            // Calculate LOOCV prediction error
            double residual = (y[i] - y_hat) / (1.0 - h_i);
            double squared_error = residual * residual;

            // Add to overall LOOCV sum
            loocv_sum += squared_error;

            // Store specific error for reference vertex
            if (i == ref_index) {
                results[ref_index].loocv_at_ref_vertex = squared_error;
            }
        }

        // LOOCV MSE is the average of squared prediction errors
        results[ref_index].loocv = loocv_sum / n_points;
    }

    return results;
}


