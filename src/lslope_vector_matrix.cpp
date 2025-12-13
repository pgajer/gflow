// lslope_vector_matrix.cpp
//
// Efficient computation of local slope between a vector y and each column of matrix Z.
// This implementation pre-computes gradient edges from y once and reuses them for all
// columns, with optional OpenMP parallelization over columns.
//
// This file should be appended to lslope.cpp or compiled separately.

#include "set_wgraph.hpp"
#include "lslope.hpp"
#include "lcor.hpp"
#include "error_utils.h"

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Sentinel value for invalid vertex (must match lslope.cpp)
static const size_t LSLOPE_VM_INVALID_VERTEX = std::numeric_limits<size_t>::max();

/**
 * @brief Compute adaptive epsilon from data for log-ratio transformations
 */
static double compute_adaptive_epsilon_vm(const std::vector<double>& f) {
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
 * @brief Compute adaptive epsilon from matrix data
 */
static double compute_adaptive_epsilon_matrix(const Eigen::MatrixXd& Z) {
    double min_nonzero = std::numeric_limits<double>::max();
    const size_t n_rows = static_cast<size_t>(Z.rows());
    const size_t n_cols = static_cast<size_t>(Z.cols());

    for (size_t col = 0; col < n_cols; ++col) {
        for (size_t row = 0; row < n_rows; ++row) {
            double val = Z(row, col);
            if (val > 0.0 && val < min_nonzero) {
                min_nonzero = val;
            }
        }
    }

    return (min_nonzero < std::numeric_limits<double>::max())
        ? 1e-6 * min_nonzero
        : 1e-6;
}

/**
 * @brief Compute edge difference based on type
 */
static inline double compute_edge_diff_vm(
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
 * @brief Compute median absolute value (for sigmoid calibration)
 */
static double compute_median_abs_vm(const std::vector<double>& values) {
    if (values.empty()) return 1.0;

    std::vector<double> abs_values;
    abs_values.reserve(values.size());

    for (double v : values) {
        if (std::isfinite(v) && v != 0.0) {
            abs_values.push_back(std::abs(v));
        }
    }

    if (abs_values.empty()) return 1.0;

    size_t n = abs_values.size();
    size_t mid = n / 2;
    std::nth_element(abs_values.begin(), abs_values.begin() + mid, abs_values.end());

    if (n % 2 == 0) {
        double upper = abs_values[mid];
        std::nth_element(abs_values.begin(), abs_values.begin() + mid - 1, abs_values.end());
        return (abs_values[mid - 1] + upper) / 2.0;
    } else {
        return abs_values[mid];
    }
}

/**
 * @brief Compute local slope between vector y and each column of matrix Z
 *
 * This method efficiently computes gradient-restricted local slope coefficients
 * between a single directing function y and multiple response functions stored
 * as columns of matrix Z.
 *
 * ALGORITHM:
 * Phase 1: Pre-compute gradient edges from y (once)
 *   - For each vertex v, find the neighbor u that maximizes (ascending) or
 *     minimizes (descending) the y-difference
 *   - Store gradient neighbor, delta_y, and extremum status
 *   - Calibrate sigmoid alpha from raw slopes of first valid column
 *
 * Phase 2: Process each column of Z (parallelized with OpenMP)
 *   - For each vertex v, compute delta_z along the pre-computed gradient edge
 *   - Apply slope transformation (raw, normalized, or sign)
 *
 * OPTIMIZATION:
 * The gradient structure depends only on y, so it's computed once and reused
 * for all q columns of Z. This gives O(n*d + q*n) complexity instead of O(q*n*d)
 * where d is average vertex degree.
 *
 * @param y Directing function values at vertices (length = num_vertices)
 * @param Z Response matrix (rows = num_vertices, columns = number of features)
 * @param slope_type Type of slope measure (GRADIENT_SLOPE, GRADIENT_SLOPE_NORMALIZED, GRADIENT_SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for Z columns (DIFFERENCE or LOGRATIO)
 * @param epsilon Pseudocount for log-ratio transformations (0 = adaptive)
 * @param sigmoid_alpha Scale for sigmoid normalization (0 = auto-calibrate)
 * @param ascending If true, use ascending gradient; if false, use descending
 * @param n_threads Number of OpenMP threads (0 = use default)
 *
 * @return lslope_vector_matrix_result_t with coefficient matrix and metadata
 */
lslope_vector_matrix_result_t set_wgraph_t::lslope_vector_matrix(
    const std::vector<double>& y,
    const Eigen::MatrixXd& Z,
    lslope_type_t slope_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double sigmoid_alpha,
    bool ascending,
    int n_threads
) const {
    const size_t n_vertices = num_vertices();
    const size_t n_columns = static_cast<size_t>(Z.cols());
    const double MIN_DELTA_Y = 1e-10;

    // ---- Input validation ----
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (static_cast<size_t>(Z.rows()) != n_vertices) {
        REPORT_ERROR("Number of rows in Z (%zu) does not match number of vertices (%zu)\n",
                     static_cast<size_t>(Z.rows()), n_vertices);
    }

    // ---- Compute adaptive epsilon if needed ----
    double epsilon_y = epsilon;
    if (y_diff_type == edge_diff_type_t::LOGRATIO && epsilon_y <= 0.0) {
        epsilon_y = compute_adaptive_epsilon_vm(y);
    }

    double epsilon_z = epsilon;
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0 && n_columns > 0) {
        epsilon_z = compute_adaptive_epsilon_matrix(Z);
    }

    // ---- Initialize result ----
    lslope_vector_matrix_result_t result;
    result.coefficients.resize(n_vertices, n_columns);
    result.coefficients.setZero();
    result.gradient_neighbors.resize(n_vertices, LSLOPE_VM_INVALID_VERTEX);
    result.gradient_delta_y.resize(n_vertices, 0.0);
    result.is_local_extremum.resize(n_vertices, false);
    result.n_local_maxima = 0;
    result.n_local_minima = 0;

    // ---- Phase 1: Pre-compute gradient edges from y ----

    // Find gradient edge for each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& neighbors = adjacency_list[v];

        if (neighbors.empty()) {
            result.is_local_extremum[v] = true;
            if (ascending) {
                result.n_local_maxima++;
            } else {
                result.n_local_minima++;
            }
            continue;
        }

        size_t best_neighbor = LSLOPE_VM_INVALID_VERTEX;
        double best_delta_y = 0.0;

        for (const auto& edge_info : neighbors) {
            size_t u = edge_info.vertex;
            double delta_y = compute_edge_diff_vm(y[u], y[v], y_diff_type, epsilon_y);

            if (ascending) {
                // Looking for steepest ascent: maximum positive delta_y
                if (delta_y > best_delta_y) {
                    best_delta_y = delta_y;
                    best_neighbor = u;
                }
            } else {
                // Looking for steepest descent: minimum (most negative) delta_y
                if (delta_y < best_delta_y) {
                    best_delta_y = delta_y;
                    best_neighbor = u;
                }
            }
        }

        if (best_neighbor == LSLOPE_VM_INVALID_VERTEX) {
            // No improvement found - v is local extremum
            result.is_local_extremum[v] = true;
            if (ascending) {
                result.n_local_maxima++;
            } else {
                result.n_local_minima++;
            }
        } else {
            result.gradient_neighbors[v] = best_neighbor;
            result.gradient_delta_y[v] = best_delta_y;
        }
    }

    // ---- Calibrate sigmoid alpha if needed ----
    double calibrated_alpha = sigmoid_alpha;

    if (slope_type == lslope_type_t::GRADIENT_SLOPE_NORMALIZED && sigmoid_alpha <= 0.0) {
        // Collect raw slopes from first column (or a sample) for calibration
        std::vector<double> raw_slopes;
        raw_slopes.reserve(n_vertices);

        size_t calibration_col = 0;
        if (n_columns > 0) {
            for (size_t v = 0; v < n_vertices; ++v) {
                if (!result.is_local_extremum[v]) {
                    size_t u = result.gradient_neighbors[v];
                    double delta_y = result.gradient_delta_y[v];

                    if (std::abs(delta_y) > MIN_DELTA_Y) {
                        double delta_z = compute_edge_diff_vm(
                            Z(u, calibration_col), Z(v, calibration_col),
                            z_diff_type, epsilon_z
                        );
                        raw_slopes.push_back(delta_z / delta_y);
                    }
                }
            }
        }

        if (!raw_slopes.empty()) {
            double median_abs = compute_median_abs_vm(raw_slopes);
            if (median_abs > 1e-10) {
                // Calibrate so tanh(alpha * median) ≈ 0.5
                // arctanh(0.5) ≈ 0.549
                calibrated_alpha = 0.549 / median_abs;
            } else {
                calibrated_alpha = 1.0;
            }
        } else {
            calibrated_alpha = 1.0;
        }
    }
    result.sigmoid_alpha = calibrated_alpha;

    // ---- Phase 2: Process each column of Z (OpenMP parallelized) ----

#ifdef _OPENMP
    if (n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
#endif

    #pragma omp parallel for schedule(dynamic)
    for (size_t col = 0; col < n_columns; ++col) {
        for (size_t v = 0; v < n_vertices; ++v) {
            if (result.is_local_extremum[v]) {
                result.coefficients(v, col) = 0.0;
                continue;
            }

            size_t u = result.gradient_neighbors[v];
            double delta_y = result.gradient_delta_y[v];

            // Skip if delta_y is too small (near-extremum)
            if (std::abs(delta_y) <= MIN_DELTA_Y) {
                result.coefficients(v, col) = 0.0;
                continue;
            }

            // Compute delta_z along gradient edge
            double delta_z = compute_edge_diff_vm(
                Z(u, col), Z(v, col), z_diff_type, epsilon_z
            );

            // Apply transformation based on slope type
            double coeff = 0.0;

            switch (slope_type) {
                case lslope_type_t::GRADIENT_SLOPE:
                    // Raw slope
                    coeff = delta_z / delta_y;
                    break;

                case lslope_type_t::GRADIENT_SLOPE_NORMALIZED:
                    // Sigmoid-normalized slope
                    coeff = std::tanh(calibrated_alpha * delta_z / delta_y);
                    break;

                case lslope_type_t::GRADIENT_SIGN:
                    // Sign only
                    if (delta_z > 0.0) {
                        coeff = 1.0;
                    } else if (delta_z < 0.0) {
                        coeff = -1.0;
                    } else {
                        coeff = 0.0;
                    }
                    break;

                default:
                    coeff = delta_z / delta_y;
                    break;
            }

            result.coefficients(v, col) = coeff;
        }
    }

    return result;
}
