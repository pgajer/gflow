#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <optional>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "adaptive_local_logistic.h"

// Constructor implementation
adaptive_params_t::adaptive_params_t(
    double gb,
    double md,
    double mc,
    double min_bf,
    double max_bf
) : global_bandwidth(gb),
    median_density(md),
    median_complexity(mc),
    min_bandwidth_factor(min_bf),
    max_bandwidth_factor(max_bf) {

    if (global_bandwidth <= 0.0) {
        Rf_error("global_bandwidth must be positive");
    }
    if (min_bandwidth_factor <= 0.0 || min_bandwidth_factor > 1.0) {
        Rf_error("min_bandwidth_factor must be in (0,1]");
    }
    if (max_bandwidth_factor < 1.0) {
        Rf_error("max_bandwidth_factor must be >= 1.0");
    }
}

/**
 * @brief Estimates local data density using kernel density estimation
 *
 * @param x Sorted predictor values
 * @param center_idx Index of the point where density is estimated
 * @param pilot_bandwidth Initial bandwidth for density estimation
 * @return Estimated local density at the center point
 */
double estimate_local_density(
    const std::vector<double>& x,
    int center_idx,
    double pilot_bandwidth) {

    double center = x[center_idx];
    double density = 0.0;
    int n = x.size();

    // Use tricube kernel for density estimation
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {
            double w = std::pow(1.0 - std::pow(dist, 3.0), 3.0);
            density += w;
        }
    }

    return density / (n * pilot_bandwidth);
}

/**
 * @brief Estimates local complexity using local polynomial regression
 *
 * @details Fits a local quadratic model and uses the magnitude of
 * the second derivative as a measure of local complexity
 */
double estimate_local_complexity(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth) {

    double center = x[center_idx];
    int n = x.size();

    // Collect local points for quadratic fit
    std::vector<double> local_x, local_y, local_w;
    for (int i = 0; i < n; ++i) {
        double dist = std::abs(x[i] - center) / pilot_bandwidth;
        if (dist < 1.0) {
            local_x.push_back(x[i] - center);  // Center x values
            local_y.push_back(y[i]);
            local_w.push_back(std::pow(1.0 - std::pow(dist, 3.0), 3.0));
        }
    }

    // Fit weighted quadratic regression
    // y = b0 + b1*x + b2*x^2
    // The absolute value of b2 indicates local complexity

    if (local_x.size() < 3) return 0.0;  // Not enough points for quadratic fit

    // Construct design matrix
    Eigen::MatrixXd X(local_x.size(), 3);
    Eigen::VectorXd Y = Eigen::Map<Eigen::VectorXd>(local_y.data(), local_y.size());
    Eigen::VectorXd W = Eigen::Map<Eigen::VectorXd>(local_w.data(), local_w.size());

    for (size_t i = 0; i < local_x.size(); ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = local_x[i];
        X(i, 2) = local_x[i] * local_x[i];
    }

    // Solve weighted least squares
    Eigen::MatrixXd WX = W.asDiagonal() * X;
    Eigen::VectorXd WY = W.asDiagonal() * Y;
    Eigen::Vector3d beta = (WX.transpose() * X).ldlt().solve(WX.transpose() * Y);

    // Return absolute value of quadratic coefficient
    return std::abs(beta(2));
}

/**
 * @brief Computes adaptive bandwidth and weights for a given point
 *
 * @details Combines density and complexity estimates to determine
 * appropriate local bandwidth, then computes kernel weights using
 * this adaptive bandwidth.
 */
adaptive_bandwidth_t compute_adaptive_bandwidth(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    const adaptive_params_t& params) {

    // Start with pilot bandwidth estimation
    double pilot_bandwidth = params.global_bandwidth;

    // Get local density and complexity estimates
    double density = estimate_local_density(x, center_idx, pilot_bandwidth);
    double complexity = estimate_local_complexity(x, y, center_idx, pilot_bandwidth);

    // Adjust local bandwidth based on density and complexity
    double density_factor = std::pow(density / params.median_density, -0.2);
    double complexity_factor = std::pow(complexity / params.median_complexity, -0.2);

    // Combine adjustments with bounds
    double bandwidth_factor = std::clamp(
        density_factor * complexity_factor,
        params.min_bandwidth_factor,
        params.max_bandwidth_factor
    );

    double local_bandwidth = pilot_bandwidth * bandwidth_factor;

    // Compute weights using adapted bandwidth
    std::vector<double> weights(x.size());
    double weight_sum = 0.0;
    double center = x[center_idx];

    for (size_t i = 0; i < x.size(); ++i) {
        double dist = std::abs(x[i] - center) / local_bandwidth;
        if (dist < 1.0) {
            weights[i] = std::pow(1.0 - std::pow(dist, 3.0), 3.0);
            weight_sum += weights[i];
        } else {
            weights[i] = 0.0;
        }
    }

    // Normalize weights
    if (weight_sum > 0.0) {
        for (auto& w : weights) {
            w /= weight_sum;
        }
    }

    return {
        local_bandwidth / (x.back() - x.front()),  // Scale to [0,1] range
        std::move(weights),
        weight_sum
    };
}
