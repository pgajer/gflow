#include <vector>
#include <cmath>

double calculate_weighted_correlation(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& weights) {

    double sum_weights = 0.0;
    double weighted_mean_x = 0.0;
    double weighted_mean_y = 0.0;

    // Calculate weighted means
    for (size_t i = 0; i < x.size(); i++) {
        sum_weights += weights[i];
        weighted_mean_x += weights[i] * x[i];
        weighted_mean_y += weights[i] * y[i];
    }

    weighted_mean_x /= sum_weights;
    weighted_mean_y /= sum_weights;

    // Calculate weighted covariance and variances
    double weighted_cov = 0.0;
    double weighted_var_x = 0.0;
    double weighted_var_y = 0.0;

    for (size_t i = 0; i < x.size(); i++) {
        double x_diff = x[i] - weighted_mean_x;
        double y_diff = y[i] - weighted_mean_y;

        weighted_cov += weights[i] * x_diff * y_diff;
        weighted_var_x += weights[i] * x_diff * x_diff;
        weighted_var_y += weights[i] * y_diff * y_diff;
    }

    weighted_cov /= sum_weights;
    weighted_var_x /= sum_weights;
    weighted_var_y /= sum_weights;

    // Calculate correlation
	// Check if variance is very small
	if (weighted_var_x <= 1e-10 || weighted_var_y <= 1e-10) {
		return 0.0; // No meaningful correlation with near-constant data
	}

    return weighted_cov / (sqrt(weighted_var_x) * sqrt(weighted_var_y));
}


/**
 * @brief Computes correlation coefficients for increasing subsequences of data points
 *
 * @details This function calculates Pearson correlation coefficients between two sequences
 *          for each prefix of length i (1 ≤ i ≤ n), using an iterative algorithm that
 *          reuses calculations from previous steps. The implementation uses running sums
 *          to achieve O(n) time complexity.
 *
 * @param[in] d_values A vector of distance values, sorted in ascending order
 * @param[in] y_values A vector of corresponding y values
 *
 * @pre d_values and y_values must have the same size
 * @pre d_values should be sorted in ascending order (d₁ < d₂ < ... < dₙ)
 *
 * @return A vector of correlation coefficients, where the i-th element represents
 *         the correlation between d_values[0:i] and y_values[0:i]
 *
 * @throws None, but behavior is undefined if input vectors have different sizes
 *
 * @note When there is no variance in either subsequence (all elements identical),
 *       the correlation for that subsequence is set to 0.0
 *
 * @complexity Time complexity: O(n), Space complexity: O(n)
 *
 * @code
 * // Example usage:
 * std::vector<double> distances = {0.0, 2.0, 3.0, 5.0, 8.0};
 * std::vector<double> values = {1.2, 2.4, 0.5, 3.1, 2.8};
 * std::vector<double> correlations = iterative_correlation(distances, values);
 * @endcode
 */
#include <vector>
#include <cmath>

std::vector<double> iterative_correlation(const std::vector<double>& d_values,
                                         const std::vector<double>& y_values) {
    size_t n = d_values.size();
    std::vector<double> correlations(n, 0.0);

    // Initialize running sums
    double sum_d = 0.0;
    double sum_y = 0.0;
    double sum_d_squared = 0.0;
    double sum_y_squared = 0.0;
    double sum_dy = 0.0;

    for (size_t i = 0; i < n; ++i) {
        // Update running sums
        sum_d += d_values[i];
        sum_y += y_values[i];
        sum_d_squared += d_values[i] * d_values[i];
        sum_y_squared += y_values[i] * y_values[i];
        sum_dy += d_values[i] * y_values[i];

        // Calculate correlation for points 1 through i+1
        double count = static_cast<double>(i + 1);
        double numerator = count * sum_dy - sum_d * sum_y;
        double denominator = std::sqrt((count * sum_d_squared - sum_d * sum_d) *
                                      (count * sum_y_squared - sum_y * sum_y));

        if (denominator != 0.0) {
            correlations[i] = numerator / denominator;
        } else {
            correlations[i] = 0.0;  // Handle the case where there's no variance
        }
    }

    return correlations;
}
