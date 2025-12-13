#include "set_wgraph.hpp"
#include "lslope.hpp"
#include "lcor.hpp"
#include "error_utils.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

// Sentinel value for invalid vertex
static const size_t LSLOPE_INVALID_VERTEX = std::numeric_limits<size_t>::max();

/**
 * @brief Compute adaptive epsilon from data for log-ratio transformations
 */
static double compute_adaptive_epsilon_lslope(const std::vector<double>& f) {
    double min_nonzero = std::numeric_limits<double>::max();

    for (double val : f) {
        if (val > 0.0 && val < min_nonzero) {
            min_nonzero = val;
        }
    }

    return (min_nonzero < std::numeric_limits<double>::max())
        ? 1e-6 * min_nonzero
        : 1e-6;
}


/**
 * @brief Compute edge difference based on type
 */
static inline double compute_edge_diff_lslope(
    double f_u,
    double f_v,
    edge_diff_type_t diff_type,
    double epsilon
) {
    switch (diff_type) {
        case edge_diff_type_t::DIFFERENCE:
            return f_u - f_v;

        case edge_diff_type_t::LOGRATIO:
            return std::log((f_u + epsilon) / (f_v + epsilon));

        default:
            return f_u - f_v;
    }
}


/**
 * @brief Find the gradient edge for a vertex
 *
 * The gradient edge at vertex v is the edge to the neighbor u that maximizes
 * the difference Δy = y(u) - y(v). If all neighbors have y <= y(v), then v
 * is a local maximum and has no ascending gradient edge.
 *
 * @param v Current vertex
 * @param y Function values at all vertices
 * @param y_diff_type Type of edge difference (DIFFERENCE or LOGRATIO)
 * @param epsilon_y Pseudocount for log-ratio (if applicable)
 * @param ascending If true, find ascending gradient (for basins of attraction to maxima)
 *                  If false, find descending gradient (for basins of attraction to minima)
 *
 * @return Pair of (gradient neighbor index, delta_y along gradient edge)
 *         If v is a local extremum, returns (INVALID_VERTEX, 0.0)
 */
std::pair<size_t, double> set_wgraph_t::find_gradient_edge(
    size_t v,
    const std::vector<double>& y,
    edge_diff_type_t y_diff_type,
    double epsilon_y,
    bool ascending
) const {
    const auto& neighbors = adjacency_list[v];

    if (neighbors.empty()) {
        return {LSLOPE_INVALID_VERTEX, 0.0};
    }

    size_t best_neighbor = LSLOPE_INVALID_VERTEX;
    double best_delta_y = 0.0;

    for (const auto& edge_info : neighbors) {
        size_t u = edge_info.vertex;
        double delta_y = compute_edge_diff_lslope(y[u], y[v], y_diff_type, epsilon_y);

        if (ascending) {
            // Looking for the steepest ascent: maximum positive delta_y
            if (delta_y > best_delta_y) {
                best_delta_y = delta_y;
                best_neighbor = u;
            }
        } else {
            // Looking for the steepest descent: minimum (most negative) delta_y
            if (delta_y < best_delta_y) {
                best_delta_y = delta_y;
                best_neighbor = u;
            }
        }
    }

    // If no improvement found, v is a local extremum
    if (best_neighbor == LSLOPE_INVALID_VERTEX) {
        return {LSLOPE_INVALID_VERTEX, 0.0};
    }

    return {best_neighbor, best_delta_y};
}


/**
 * @brief Compute gradient-restricted local slope measures
 *
 * For each vertex v, computes the local slope of z with respect to y along
 * the gradient direction of y. The gradient edge is the edge from v to the
 * neighbor u that maximizes Δy = y(u) - y(v).
 *
 * Three variants are supported:
 * - GRADIENT_SLOPE: Raw ratio Δz/Δy
 * - GRADIENT_SLOPE_NORMALIZED: Sigmoid-normalized ratio tanh(α·Δz/Δy)
 * - GRADIENT_SIGN: Sign of Δz along gradient direction
 *
 * @param y Directing function values at vertices (length = num_vertices)
 * @param z Response function values at vertices (length = num_vertices)
 * @param slope_type Type of slope measure to compute
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for z (DIFFERENCE or LOGRATIO)
 * @param epsilon Pseudocount for log-ratio (0 = adaptive)
 * @param sigmoid_alpha Scale for sigmoid normalization (0 = calibrate from data)
 * @param sigmoid_type Type of sigmoid function
 * @param ascending If true, use ascending gradient; if false, use descending
 *
 * @return lslope_result_t structure with coefficients and diagnostics
 */
lslope_result_t set_wgraph_t::lslope_gradient(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lslope_type_t slope_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double sigmoid_alpha,
    sigmoid_type_t sigmoid_type,
    bool ascending
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DELTA_Y = 1e-10;

    // Validate input lengths
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (z.size() != n_vertices) {
        REPORT_ERROR("Length of z (%zu) does not match number of vertices (%zu)\n",
                     z.size(), n_vertices);
    }

    // Compute adaptive epsilon if needed
    double epsilon_y = epsilon;
    double epsilon_z = epsilon;

    if (y_diff_type == edge_diff_type_t::LOGRATIO && epsilon_y <= 0.0) {
        epsilon_y = compute_adaptive_epsilon_lslope(y);
    }
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0) {
        epsilon_z = compute_adaptive_epsilon_lslope(z);
    }

    // Initialize result structure
    lslope_result_t result;
    result.vertex_coefficients.resize(n_vertices, 0.0);
    result.gradient_neighbors.resize(n_vertices, LSLOPE_INVALID_VERTEX);
    result.gradient_delta_y.resize(n_vertices, 0.0);
    result.gradient_delta_z.resize(n_vertices, 0.0);
    result.is_local_extremum.resize(n_vertices, false);

    // First pass: find gradient edges and compute raw slopes
    std::vector<double> raw_slopes;
    raw_slopes.reserve(n_vertices);

    for (size_t v = 0; v < n_vertices; ++v) {
        auto [grad_neighbor, delta_y] = find_gradient_edge(v, y, y_diff_type, epsilon_y, ascending);

        result.gradient_neighbors[v] = grad_neighbor;
        result.gradient_delta_y[v] = delta_y;

        if (grad_neighbor == LSLOPE_INVALID_VERTEX) {
            // v is a local extremum
            result.is_local_extremum[v] = true;
            if (ascending) {
                result.n_local_maxima++;
            } else {
                result.n_local_minima++;
            }
            result.vertex_coefficients[v] = 0.0;
        } else {
            // Compute delta_z along gradient edge
            double delta_z = compute_edge_diff_lslope(
                z[grad_neighbor], z[v], z_diff_type, epsilon_z
            );
            result.gradient_delta_z[v] = delta_z;

            // Compute raw slope
            if (std::abs(delta_y) > MIN_DELTA_Y) {
                double raw_slope = delta_z / delta_y;
                raw_slopes.push_back(raw_slope);

                // Store based on slope type
                if (slope_type == lslope_type_t::GRADIENT_SIGN) {
                    // Just the sign
                    if (delta_z > 0.0) {
                        result.vertex_coefficients[v] = 1.0;
                    } else if (delta_z < 0.0) {
                        result.vertex_coefficients[v] = -1.0;
                    } else {
                        result.vertex_coefficients[v] = 0.0;
                    }
                } else {
                    // Store raw slope for now (will normalize in second pass if needed)
                    result.vertex_coefficients[v] = raw_slope;
                }
            } else {
                // Near-zero delta_y: treat as local extremum for numerical stability
                result.is_local_extremum[v] = true;
                result.vertex_coefficients[v] = 0.0;
            }
        }
    }

    // Second pass: apply sigmoid normalization if requested
    if (slope_type == lslope_type_t::GRADIENT_SLOPE_NORMALIZED) {

        // Calibrate sigmoid_alpha based on data when user requests automatic calibration
        // (indicated by sigmoid_alpha <= 0)
        //
        // For each sigmoid type σ, we solve for alpha such that:
        //   σ(alpha * median_abs) = 0.5

        // The median represents a "typical" slope magnitude in the sense that
        // half of the absolute slopes are smaller and half are larger. Mapping
        // this typical value to the midpoint of the sigmoid's positive range
        // ensures balanced utilization of the output range. Slopes with
        // absolute value smaller than the median_abs produce outputs in (-0.5,
        // 0.5), while slopes larger than the median produce outputs in
        // (-1, -0.5) \cup (0.5, 1).
        // This prevents both saturation (where α is too large and nearly all
        // outputs cluster near ±1) and compression (where α is too small and
        // all outputs cluster near 0).

        if (sigmoid_alpha <= 0.0) {
            double median_abs = compute_median_abs(raw_slopes);

            // Guard against degenerate case of zero or near-zero slopes
            const double min_median = 1e-10;
            if (median_abs < min_median) {
                median_abs = min_median;
            }

            if (sigmoid_type == sigmoid_type_t::TANH) {
                // Hyperbolic tangent: sigma(x) = tanh(alpha * x)
                //
                // Note for those familiar with logistic regression: the rescaled logistic
                // sigma_L(x) = 2/(1 + exp(-x)) - 1, which has range (-1, 1), is equivalent
                // to tanh with a scale factor of 1/2. To see this:
                //
                //   2/(1 + exp(-x)) - 1 = (2 - 1 - exp(-x)) / (1 + exp(-x))
                //                       = (1 - exp(-x)) / (1 + exp(-x))
                //
                // Multiplying numerator and denominator by exp(x/2):
                //
                //                       = (exp(x/2) - exp(-x/2)) / (exp(x/2) + exp(-x/2))
                //                       = tanh(x/2)
                //
                // Thus tanh(alpha * x) = 2/(1 + exp(-2*alpha*x)) - 1, and there is no need
                // for a separate rescaled logistic sigmoid type.
                //
                // Solve: tanh(alpha * m) = 0.5
                // => alpha * m = arctanh(0.5)
                // => alpha = arctanh(0.5) / m
                sigmoid_alpha = std::atanh(0.5) / median_abs;

            } else if (sigmoid_type == sigmoid_type_t::ARCTAN) {
                // Scaled arctan: sigma(x) = (2/pi) * arctan(x)
                // Solve: (2/pi) * arctan(alpha * m) = 0.5
                // => arctan(alpha * m) = pi/4
                // => alpha * m = tan(pi/4) = 1
                // => alpha = 1 / m
                sigmoid_alpha = 1.0 / median_abs;

            } else if (sigmoid_type == sigmoid_type_t::ALGEBRAIC) {
                // Algebraic sigmoid: sigma(x) = x / sqrt(alpha^{-2} + x^2)
                // Solve: m / sqrt(alpha^{-2} + m^2) = 0.5
                // => m^2 / (alpha^{-2} + m^2) = 0.25
                // => 4 * m^2 = alpha^{-2} + m^2
                // => 3 * m^2 = alpha^{-2}
                // => alpha = 1 / (sqrt(3) * m)
                sigmoid_alpha = 1.0 / (std::sqrt(3.0) * median_abs);
            }
        }

        result.sigmoid_alpha = sigmoid_alpha;

        for (size_t v = 0; v < n_vertices; ++v) {
            if (!result.is_local_extremum[v]) {
                result.vertex_coefficients[v] = apply_sigmoid(
                    result.vertex_coefficients[v],
                    sigmoid_alpha,
                    sigmoid_type
                );
            }
        }
    }

    // Compute summary statistics
    size_t n_valid = 0;
    double sum = 0.0;
    std::vector<double> valid_coeffs;
    valid_coeffs.reserve(n_vertices);

    for (size_t v = 0; v < n_vertices; ++v) {
        if (!result.is_local_extremum[v]) {
            double coef = result.vertex_coefficients[v];
            sum += coef;
            n_valid++;
            valid_coeffs.push_back(coef);

            if (coef > 0.0) {
                result.n_positive++;
            } else if (coef < 0.0) {
                result.n_negative++;
            } else {
                result.n_zero++;
            }
        }
    }

    if (n_valid > 0) {
        result.mean_coefficient = sum / static_cast<double>(n_valid);
        result.median_coefficient = compute_median(valid_coeffs);
    }

    return result;
}


/**
 * @brief Compute neighborhood local regression coefficient
 *
 * For each vertex v, computes the local regression coefficient of z on y
 * using all edges in the neighborhood:
 *
 *   β_loc(z; y, w)(v) = Σ w_e Δ_e y · Δ_e z / Σ w_e (Δ_e y)²
 *
 * This equals lcor(y, z) × sd_loc(z) / sd_loc(y).
 *
 * @param y Directing function values at vertices
 * @param z Response function values at vertices
 * @param weight_type Weighting scheme (UNIT or DERIVATIVE)
 * @param y_diff_type Edge difference type for y
 * @param z_diff_type Edge difference type for z
 * @param epsilon Pseudocount for log-ratio (0 = adaptive)
 * @param winsorize_quantile Quantile for winsorization (0 = none)
 *
 * @return lslope_nbhd_result_t structure with coefficients and diagnostics
 */
lslope_nbhd_result_t set_wgraph_t::lslope_neighborhood(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;
    const double MIN_EDGE_LENGTH = 1e-10;

    // Validate input lengths
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (z.size() != n_vertices) {
        REPORT_ERROR("Length of z (%zu) does not match number of vertices (%zu)\n",
                     z.size(), n_vertices);
    }

    // Compute adaptive epsilon if needed
    double epsilon_y = epsilon;
    double epsilon_z = epsilon;

    if (y_diff_type == edge_diff_type_t::LOGRATIO && epsilon_y <= 0.0) {
        epsilon_y = compute_adaptive_epsilon_lslope(y);
    }
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0) {
        epsilon_z = compute_adaptive_epsilon_lslope(z);
    }

    // Initialize result structure
    lslope_nbhd_result_t result;
    result.vertex_coefficients.resize(n_vertices, 0.0);
    result.sd_y.resize(n_vertices, 0.0);
    result.sd_z.resize(n_vertices, 0.0);
    result.lcor.resize(n_vertices, 0.0);

    // Compute coefficients at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;      // Σ w_e Δ_e y · Δ_e z
        double sum_y_squared = 0.0;  // Σ w_e (Δ_e y)²
        double sum_z_squared = 0.0;  // Σ w_e (Δ_e z)²
        double sum_weights = 0.0;    // Σ w_e

        for (const auto& edge_info : adjacency_list[v]) {
            size_t u = edge_info.vertex;

            // Compute weight
            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;

                case lcor_type_t::DERIVATIVE:
                    if (edge_info.weight > MIN_EDGE_LENGTH) {
                        weight = 1.0 / (edge_info.weight * edge_info.weight);
                    } else {
                        continue;  // Skip degenerate edges
                    }
                    break;

                case lcor_type_t::SIGN:
                    weight = 1.0;
                    break;
            }

            // Compute edge differences
            double delta_y = compute_edge_diff_lslope(y[u], y[v], y_diff_type, epsilon_y);
            double delta_z = compute_edge_diff_lslope(z[u], z[v], z_diff_type, epsilon_z);

            // Accumulate
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
            sum_weights += weight;
        }

        // Compute local regression coefficient: β = Σ(w·Δy·Δz) / Σ(w·Δy²)
        if (sum_y_squared > MIN_DENOMINATOR) {
            result.vertex_coefficients[v] = numerator / sum_y_squared;
        } else {
            result.vertex_coefficients[v] = 0.0;
        }

        // Compute local standard deviations
        if (sum_weights > MIN_DENOMINATOR) {
            result.sd_y[v] = std::sqrt(sum_y_squared / sum_weights);
            result.sd_z[v] = std::sqrt(sum_z_squared / sum_weights);
        }

        // Compute local correlation for reference
        double denom_y = std::sqrt(sum_y_squared);
        double denom_z = std::sqrt(sum_z_squared);
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result.lcor[v] = numerator / (denom_y * denom_z);
            // Clamp to [-1, 1]
            if (result.lcor[v] > 1.0) result.lcor[v] = 1.0;
            if (result.lcor[v] < -1.0) result.lcor[v] = -1.0;
        }
    }

    // Compute summary statistics
    double sum = 0.0;
    std::vector<double> coeffs = result.vertex_coefficients;
    for (size_t v = 0; v < n_vertices; ++v) {
        sum += result.vertex_coefficients[v];
    }
    result.mean_coefficient = sum / static_cast<double>(n_vertices);
    result.median_coefficient = compute_median(coeffs);

    return result;
}


/**
 * @brief Compute local slope (production version returning only coefficients)
 *
 * Unified entry point that dispatches to gradient-restricted or neighborhood
 * implementations based on slope_type parameter.
 *
 * @param y Directing function values
 * @param z Response function values
 * @param slope_type Type of slope measure
 * @param weight_type Weighting scheme (for NEIGHBORHOOD_SLOPE)
 * @param y_diff_type Edge difference type for y
 * @param z_diff_type Edge difference type for z
 * @param epsilon Pseudocount for log-ratio
 * @param sigmoid_alpha Scale for sigmoid (0 = auto-calibrate)
 * @param ascending Use ascending gradient (true) or descending (false)
 *
 * @return Vector of local slope coefficients
 */
std::vector<double> set_wgraph_t::lslope(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lslope_type_t slope_type,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double sigmoid_alpha,
    bool ascending
) const {
    if (slope_type == lslope_type_t::NEIGHBORHOOD_SLOPE) {
        auto result = lslope_neighborhood(
            y, z, weight_type, y_diff_type, z_diff_type, epsilon, 0.0
        );
        return result.vertex_coefficients;
    } else {
        auto result = lslope_gradient(
            y, z, slope_type, y_diff_type, z_diff_type,
            epsilon, sigmoid_alpha, sigmoid_type_t::TANH, ascending
        );
        return result.vertex_coefficients;
    }
}
