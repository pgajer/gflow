#include "set_wgraph.hpp"
#include "lcor.hpp"
#include "error_utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

/**
 * @brief Compute adaptive epsilon from data
 *
 * For log-ratio transformations, we need a pseudocount ε to handle zeros.
 * This function computes ε = 1e-6 * min(non-zero values), which is small
 * enough to minimally perturb the data but large enough to avoid numerical
 * issues with log(0).
 *
 * @param f Vector of function values
 * @return Adaptive epsilon (1e-6 if all values are zero)
 */
static double compute_adaptive_epsilon(const std::vector<double>& f) {
    double min_nonzero = std::numeric_limits<double>::max();
    
    for (double val : f) {
        if (val > 0.0 && val < min_nonzero) {
            min_nonzero = val;
        }
    }
    
    // If no positive values found, return a default small epsilon
    return (min_nonzero < std::numeric_limits<double>::max()) 
        ? 1e-6 * min_nonzero 
        : 1e-6;
}

/**
 * @brief Compute edge difference based on type
 *
 * Computes the directional difference Δ_e f along an edge e = [v,u]:
 * - DIFFERENCE: Δ_e f = f(u) - f(v)
 * - LOGRATIO: Δ_e f = log((f(u) + ε) / (f(v) + ε))
 *
 * @param f_u Function value at vertex u (edge target)
 * @param f_v Function value at vertex v (edge source)
 * @param diff_type Type of difference to compute
 * @param epsilon Pseudocount for log-ratio (ignored for DIFFERENCE)
 * @return Edge difference Δ_e f
 */
static inline double compute_edge_diff(
    double f_u, 
    double f_v,
    edge_diff_type_t diff_type,
    double epsilon
) {
    switch (diff_type) {
        case edge_diff_type_t::DIFFERENCE:
            return f_u - f_v;
            
        case edge_diff_type_t::LOGRATIO:
            // log((f_u + ε) / (f_v + ε))
            // Both f_u and f_v should be non-negative for compositional data
            // Adding epsilon ensures we never take log(0)
            return std::log((f_u + epsilon) / (f_v + epsilon));
            
        default:
            // Should never reach here, but default to difference
            return f_u - f_v;
    }
}

/**
 * @brief Compute winsorization bounds from data
 *
 * Winsorization clips extreme values to reduce outlier influence.
 * Given a quantile q, values below the q-th percentile are clipped to
 * that percentile, and values above the (1-q)-th percentile are clipped
 * to that percentile.
 *
 * For example, q = 0.05 clips to the 5th and 95th percentiles.
 *
 * @param values Vector of values to winsorize (copied for sorting)
 * @param quantile Quantile for clipping (0 = no clipping, 0.05 = 5%/95%)
 * @return Pair of (lower_bound, upper_bound)
 */
static std::pair<double, double> compute_winsorize_bounds(
    std::vector<double> values,  // Copy for sorting
    double quantile
) {
    // No winsorization requested or invalid quantile
    if (quantile <= 0.0 || quantile >= 0.5) {
        return {-std::numeric_limits<double>::max(), 
                std::numeric_limits<double>::max()};
    }
    
    if (values.empty()) {
        return {0.0, 0.0};
    }
    
    // Sort to find percentiles
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    
    // Compute indices for quantiles
    // Use floor for lower, ceil for upper to be conservative
    size_t lower_idx = static_cast<size_t>(std::floor(n * quantile));
    size_t upper_idx = static_cast<size_t>(std::ceil(n * (1.0 - quantile))) - 1;
    
    // Bounds checking
    if (lower_idx >= n) lower_idx = 0;
    if (upper_idx >= n) upper_idx = n - 1;
    
    return {values[lower_idx], values[upper_idx]};
}

/**
 * @brief Apply winsorization clipping to a value
 *
 * @param val Value to clip
 * @param lower Lower bound
 * @param upper Upper bound
 * @return Clipped value in [lower, upper]
 */
static inline double winsorize_clip(double val, double lower, double upper) {
    if (val < lower) return lower;
    if (val > upper) return upper;
    return val;
}

/**
 * @brief Compute local correlation coefficients with diagnostic output (one-pass, no winsorization)
 *
 * Instrumented single-pass implementation that captures edge-level diagnostic data
 * in addition to vertex coefficients. Use this version for algorithm validation,
 * parameter tuning, and detailed analysis of edge-wise contributions.
 *
 * This function stores all intermediate edge differences and weights, enabling
 * post-hoc analysis of:
 * - Distribution of edge differences for y and z
 * - Per-vertex breakdown of contributing edges
 * - Identification of influential edges or outliers
 *
 * COMPUTATIONAL OVERHEAD:
 * Compared to the production lcor_one_pass(), this version allocates O(|E|)
 * additional memory for storing edge differences and weights.
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme (UNIT, DERIVATIVE, or SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for z (DIFFERENCE or LOGRATIO)
 * @param epsilon_y Pseudocount for y when using LOGRATIO (pre-computed)
 * @param epsilon_z Pseudocount for z when using LOGRATIO (pre-computed)
 *
 * @return lcor_result_t structure containing:
 *   - vertex_coefficients: Vector of length n with lcor(y,z)(v) for each vertex v
 *   - vertex_delta_y: Edge differences for y at each vertex (vector of vectors)
 *   - vertex_delta_z: Edge differences for z at each vertex (vector of vectors)
 *   - vertex_weights: Edge weights at each vertex (vector of vectors)
 *   - all_delta_y: All edge differences for y across the graph
 *   - all_delta_z: All edge differences for z across the graph
 *   - y_lower, y_upper, z_lower, z_upper: Winsorization bounds (set to ±max for this version)
 *
 * @note This is an internal helper called by lcor_instrumented(). For production
 *       use without diagnostic overhead, use lcor() instead.
 *
 * @see lcor_instrumented() for the public entry point
 * @see lcor_one_pass() for the streamlined production version
 */
lcor_result_t set_wgraph_t::lcor_one_pass_instrumented(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon_y,
    double epsilon_z
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;
    const double MIN_EDGE_LENGTH = 1e-10;

    std::vector<double> result(n_vertices, 0.0);
    std::vector<std::vector<double>> vertex_delta_y(n_vertices);
    std::vector<std::vector<double>> vertex_delta_z(n_vertices);
    std::vector<std::vector<double>> vertex_weights(n_vertices);

    std::vector<double> all_delta_y;
    std::vector<double> all_delta_z;

    // Estimate size for reservation
    size_t total_edges = 0;
    for (size_t v = 0; v < n_vertices; ++v) {
        total_edges += adjacency_list[v].size();
    }
    all_delta_y.reserve(total_edges);
    all_delta_z.reserve(total_edges);

    // Compute coefficient at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;
        double sum_y_squared = 0.0;
        double sum_z_squared = 0.0;

        vertex_delta_y[v].reserve(adjacency_list[v].size());
        vertex_delta_z[v].reserve(adjacency_list[v].size());
        vertex_weights[v].reserve(adjacency_list[v].size());

        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;

            const double delta_y = compute_edge_diff(
                y[u], y[v], y_diff_type, epsilon_y
            );
            const double delta_z = compute_edge_diff(
                z[u], z[v], z_diff_type, epsilon_z
            );

            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;

                case lcor_type_t::DERIVATIVE:
                    if (edge_info.weight > MIN_EDGE_LENGTH) {
                        weight = 1.0 / (edge_info.weight * edge_info.weight);
                    } else {
                        continue;
                    }
                    break;

                case lcor_type_t::SIGN:
                    weight = 1.0;
                    break;
            }

            // Store per-vertex data
            vertex_delta_y[v].push_back(delta_y);
            vertex_delta_z[v].push_back(delta_z);
            vertex_weights[v].push_back(weight);

            // Store global data for diagnostics
            all_delta_y.push_back(delta_y);
            all_delta_z.push_back(delta_z);

            // Accumulate weighted sums
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
        }

        // [rest of coefficient computation unchanged]
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);

        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result[v] = numerator / (denom_y * denom_z);
            if (result[v] > 1.0) {
                result[v] = 1.0;
            } else if (result[v] < -1.0) {
                result[v] = -1.0;
            }
        } else {
            result[v] = 0.0;
        }
    }

    lcor_result_t lcor_result;
    lcor_result.vertex_coefficients = std::move(result);
    lcor_result.vertex_delta_y = std::move(vertex_delta_y);
    lcor_result.vertex_delta_z = std::move(vertex_delta_z);
    lcor_result.vertex_weights = std::move(vertex_weights);
    lcor_result.all_delta_y = std::move(all_delta_y);
    lcor_result.all_delta_z = std::move(all_delta_z);

    return lcor_result;
}

/**
 * @brief Compute local correlation coefficients with diagnostic output (two-pass, with winsorization)
 *
 * Instrumented two-pass implementation that captures edge-level diagnostic data
 * including raw (pre-winsorization) edge differences and computed winsorization
 * bounds. Use this version for validating winsorization behavior and analyzing
 * the impact of outlier clipping.
 *
 * Pass 1: Compute all edge differences, store for diagnostics, determine clipping bounds
 * Pass 2: Compute coefficients using winsorized differences
 *
 * DIAGNOSTIC OUTPUT:
 * - all_delta_y, all_delta_z: Raw edge differences BEFORE winsorization
 * - vertex_delta_y, vertex_delta_z: Edge differences (also before winsorization)
 * - y_lower, y_upper, z_lower, z_upper: Computed winsorization bounds
 *
 * This enables analysis of:
 * - How many edges were clipped by winsorization
 * - The distribution of raw vs. clipped differences
 * - Sensitivity of results to the winsorization quantile
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme (UNIT, DERIVATIVE, or SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for z (DIFFERENCE or LOGRATIO)
 * @param epsilon_y Pseudocount for y when using LOGRATIO (pre-computed)
 * @param epsilon_z Pseudocount for z when using LOGRATIO (pre-computed)
 * @param winsorize_quantile Quantile for clipping (e.g., 0.05 for 5%/95%)
 *
 * @return lcor_result_t structure containing:
 *   - vertex_coefficients: Vector of length n with lcor(y,z)(v) using winsorized differences
 *   - vertex_delta_y: Raw edge differences for y at each vertex (before winsorization)
 *   - vertex_delta_z: Raw edge differences for z at each vertex (before winsorization)
 *   - vertex_weights: Edge weights at each vertex
 *   - all_delta_y: All raw edge differences for y across the graph
 *   - all_delta_z: All raw edge differences for z across the graph
 *   - y_lower, y_upper: Winsorization bounds for y differences
 *   - z_lower, z_upper: Winsorization bounds for z differences
 *
 * @note This is an internal helper called by lcor_instrumented(). For production
 *       use without diagnostic overhead, use lcor() instead.
 *
 * @see lcor_instrumented() for the public entry point
 * @see lcor_two_pass() for the streamlined production version
 */
lcor_result_t set_wgraph_t::lcor_two_pass_instrumented(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon_y,
    double epsilon_z,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;
    
    // Pass 1: Collect all edge differences
    std::vector<double> all_delta_y;
    std::vector<double> all_delta_z;
    
    // Estimate total number of edges for efficient allocation
    size_t total_edges = 0;
    for (size_t v = 0; v < n_vertices; ++v) {
        total_edges += adjacency_list[v].size();
    }
    all_delta_y.reserve(total_edges);
    all_delta_z.reserve(total_edges);
    
    // Store edge differences and weights for each vertex
    // Using vectors-of-vectors to maintain vertex-neighborhood structure
    std::vector<std::vector<double>> vertex_delta_y(n_vertices);
    std::vector<std::vector<double>> vertex_delta_z(n_vertices);
    std::vector<std::vector<double>> vertex_weights(n_vertices);
    
    for (size_t v = 0; v < n_vertices; ++v) {
        const size_t n_neighbors = adjacency_list[v].size();
        vertex_delta_y[v].reserve(n_neighbors);
        vertex_delta_z[v].reserve(n_neighbors);
        vertex_weights[v].reserve(n_neighbors);
        
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;

            // Compute weight
            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;
                    
                case lcor_type_t::DERIVATIVE:
                    if (edge_info.weight > 1e-10) {
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
            const double delta_y = compute_edge_diff(
                y[u], y[v], y_diff_type, epsilon_y
            );
            const double delta_z = compute_edge_diff(
                z[u], z[v], z_diff_type, epsilon_z
            );
            
            // Store for this vertex
            vertex_delta_y[v].push_back(delta_y);
            vertex_delta_z[v].push_back(delta_z);
            vertex_weights[v].push_back(weight);
            
            // Collect for global winsorization bounds
            all_delta_y.push_back(delta_y);
            all_delta_z.push_back(delta_z);
        }
    }
    
    // Compute winsorization bounds across all edges
    auto [y_lower, y_upper] = compute_winsorize_bounds(
        all_delta_y, winsorize_quantile
    );
    auto [z_lower, z_upper] = compute_winsorize_bounds(
        all_delta_z, winsorize_quantile
    );
    
    // Pass 2: Compute coefficients with winsorized differences
    std::vector<double>  result(n_vertices, 0.0);
    
    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;
        double sum_y_squared = 0.0;
        double sum_z_squared = 0.0;
        
        const size_t n_edges = vertex_delta_y[v].size();
        for (size_t i = 0; i < n_edges; ++i) {
            // Apply winsorization clipping
            const double delta_y = winsorize_clip(
                vertex_delta_y[v][i], y_lower, y_upper
            );
            const double delta_z = winsorize_clip(
                vertex_delta_z[v][i], z_lower, z_upper
            );
            const double weight = vertex_weights[v][i];
            
            // Accumulate weighted sums
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
        }
        
        // Compute denominators
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);
        
        // Compute coefficient with numerical stability check
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result[v] = numerator / (denom_y * denom_z);
            
            // Clamp to [-1, 1]
            if (result[v] > 1.0) {
                result[v] = 1.0;
            } else if (result[v] < -1.0) {
                result[v] = -1.0;
            }
        } else {
            result[v] = 0.0;
        }
    }

    lcor_result_t lcor_result;
    lcor_result.vertex_coefficients = std::move(result);
    lcor_result.vertex_delta_y      = std::move(vertex_delta_y);
    lcor_result.vertex_delta_z      = std::move(vertex_delta_z);
    lcor_result.vertex_weights      = std::move(vertex_weights);
    lcor_result.all_delta_y         = std::move(all_delta_y);
    lcor_result.all_delta_z         = std::move(all_delta_z);
    lcor_result.y_lower             = y_lower;
    lcor_result.y_upper             = y_upper;
    lcor_result.z_lower             = z_lower;
    lcor_result.z_upper             = z_upper;

    return lcor_result;
}

/**
 * @brief Compute local correlation coefficients with full diagnostic output (instrumented version)
 *
 * Computes vertex-level correlation-type coefficients measuring the alignment
 * of directional changes between functions y and z:
 *
 *   lcor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * where the sum is over edges incident to vertex v, and Δ_e f represents the
 * directional difference along edge e = [v,u]. The key innovation is the flexible
 * choice of how to compute edge differences:
 *
 * STANDARD DIFFERENCES (edge_diff_type_t::DIFFERENCE):
 *   Δ_e f = f(u) - f(v)
 *   Appropriate for continuous data in Euclidean space where additive changes
 *   are meaningful (e.g., temperature, pH, expression levels).
 *
 * LOG-RATIOS (edge_diff_type_t::LOGRATIO):
 *   Δ_e f = log((f(u) + ε) / (f(v) + ε))
 *   Appropriate for compositional data (relative abundances, proportions) where
 *   multiplicative changes are more meaningful. The log transformation:
 *   - Maps ratios to a symmetric additive scale
 *   - Corresponds to the Aitchison distance on the simplex
 *   - Makes correlation interpretable as gradient alignment in the natural geometry
 *
 * This flexibility allows appropriate treatment of different data types:
 * - y and z both continuous: Use DIFFERENCE for both
 * - y continuous, z compositional: Use DIFFERENCE for y, LOGRATIO for z
 * - y and z both compositional: Use LOGRATIO for both
 *
 * GEOMETRIC INTERPRETATION:
 * For smooth functions on manifolds with derivative weighting (w_e = 1/ℓ_e²),
 * this coefficient converges to cos(θ), where θ is the angle between gradient
 * vectors ∇y and ∇z (or their appropriate generalizations for log-ratios).
 * - Values near +1: functions increase together (parallel gradients)
 * - Values near -1: functions vary oppositely (anti-parallel gradients)
 * - Values near 0: functions vary independently (orthogonal gradients)
 *
 * WINSORIZATION FOR ROBUSTNESS:
 * When winsorize_quantile > 0, extreme edge differences are clipped to reduce
 * outlier influence. For example, winsorize_quantile = 0.05 clips differences
 * to the 5th and 95th percentiles. This trades some sensitivity to extreme
 * changes for robustness against technical artifacts or measurement noise.
 * Note: Winsorization uses a two-pass algorithm and is more computationally
 * expensive than the default one-pass approach.
 *
 * NUMERICAL STABILITY:
 * Returns 0 when either function has insufficient variation at a vertex
 * (denominator < 1e-10), which occurs when the function is essentially constant
 * across all incident edges. For log-ratios, the pseudocount ε prevents log(0).
 * If ε = 0 (default), it is computed adaptively as 1e-6 times the minimum
 * non-zero value in the data.
 *
 * COMPARISON TO comono_cor():
 * This function generalizes comono_cor() by allowing flexible edge difference
 * types. When y_diff_type = z_diff_type = DIFFERENCE, it is equivalent to
 * comono_cor() with the same weight_type.
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme:
 *   - UNIT: w_e = 1 (combinatorial, ignores edge lengths)
 *   - DERIVATIVE: w_e = 1/ℓ_e² (geometric, compares directional derivatives)
 * @param y_diff_type How to compute edge differences for y:
 *   - DIFFERENCE: Δ_e y = y(u) - y(v)
 *   - LOGRATIO: Δ_e y = log((y(u) + ε) / (y(v) + ε))
 * @param z_diff_type How to compute edge differences for z (same options as y)
 * @param epsilon Pseudocount for log-ratio transformations (default 0 = adaptive)
 *   If 0, computed as 1e-6 * min(non-zero values) for each function.
 *   Only used when corresponding diff_type is LOGRATIO.
 * @param winsorize_quantile Quantile for winsorization clipping (default 0 = none)
 *   If > 0, clips edge differences to [q, 1-q] percentiles to reduce outlier
 *   influence. Common values: 0.05 (5%/95% clip), 0.10 (10%/90% clip).
 *   Uses more expensive two-pass algorithm.
 *
 * @return lcor_result_t structure containing:
 *   - vertex_coefficients: Vector of length n with lcor(y,z)(v) for each vertex v
 *   - vertex_delta_y: Edge differences for y at each vertex
 *   - vertex_delta_z: Edge differences for z at each vertex
 *   - vertex_weights: Edge weights at each vertex
 *   - y_lower, y_upper: Winsorization bounds for y (only meaningful with winsorization)
 *   - z_lower, z_upper: Winsorization bounds for z (only meaningful with winsorization)
 *
 * EXAMPLE USAGE:
 * @code
 * set_wgraph_t graph(adj_list, weight_list);
 * std::vector<double> response = ...;           // Continuous outcome
 * std::vector<double> abundance = ...;          // Relative abundance (compositional)
 *
 * // Compare continuous response to compositional feature using appropriate transformations
 * auto result = graph.lcor_instrumented(
 *     response,                              // y: continuous data
 *     abundance,                             // z: compositional data
 *     lcor_type_t::DERIVATIVE,            // Geometric weighting
 *     edge_diff_type_t::DIFFERENCE,         // Standard differences for y
 *     edge_diff_type_t::LOGRATIO,           // Log-ratios for z (compositional)
 *     0.0,                                   // Adaptive epsilon
 *     0.0                                    // No winsorization
 * );
 *
 * // Examine results
 * std::cout << "Mean correlation: " << result.mean_coefficient << std::endl;
 * std::cout << "Positive vertices: " << result.n_positive << std::endl;
 *
 * // With winsorization for robustness
 * auto robust_result = graph.lcor_intrumented(
 *     response, abundance,
 *     lcor_type_t::DERIVATIVE,
 *     edge_diff_type_t::DIFFERENCE,
 *     edge_diff_type_t::LOGRATIO,
 *     1e-6,        // Explicit epsilon
 *     0.05         // Clip to 5th/95th percentiles
 * );
 * @endcode
 *
 * WHEN TO USE DIFFERENT DIFFERENCE TYPES:
 * 
 * Both DIFFERENCE:
 * - Standard regression scenarios
 * - Gene expression vs. clinical outcomes
 * - Physical measurements (temperature, pressure, etc.)
 *
 * Mixed DIFFERENCE/LOGRATIO:
 * - Clinical outcome (DIFFERENCE) vs. microbiome abundance (LOGRATIO)
 * - Continuous predictor (DIFFERENCE) vs. compositional response (LOGRATIO)
 *
 * Both LOGRATIO:
 * - Co-occurrence analysis in compositional data
 * - Relative abundance vs. relative abundance
 * - Any two compositional features
 *
 * @see comono_cor() for the standard difference-only version
 * @see edge_diff_type_t for details on difference type selection
 */
lcor_result_t set_wgraph_t::lcor_instrumented(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();
    
    // Validate input lengths
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (z.size() != n_vertices) {
        REPORT_ERROR("Length of z (%zu) does not match number of vertices (%zu)\n",
                     z.size(), n_vertices);
    }
    
    // Compute adaptive epsilon if needed (epsilon = 0 means adaptive)
    double epsilon_y = epsilon;
    double epsilon_z = epsilon;
    
    if (y_diff_type == edge_diff_type_t::LOGRATIO && epsilon_y <= 0.0) {
        epsilon_y = compute_adaptive_epsilon(y);
    }
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0) {
        epsilon_z = compute_adaptive_epsilon(z);
    }
    
    // Choose implementation based on whether winsorization is requested
    const bool need_winsorization = (winsorize_quantile > 0.0 && 
                                     winsorize_quantile < 0.5);
    
    if (need_winsorization) {
        // Two-pass algorithm: collect differences, compute bounds, then correlate
        return lcor_two_pass_instrumented(
            y, z,
            weight_type, y_diff_type, z_diff_type,
            epsilon_y, epsilon_z, winsorize_quantile
        );
    } else {
        // One-pass algorithm: compute correlations directly (more efficient)
        return lcor_one_pass_instrumented(
            y, z,
            weight_type, y_diff_type, z_diff_type,
            epsilon_y, epsilon_z
        );
    }
}

// ----------------------------------------------------------------------------
// production lcor
// ----------------------------------------------------------------------------

/**
 * @brief Compute local correlation coefficients (one-pass, no winsorization) - production version
 *
 * Efficient single-pass implementation for computing local correlation coefficients
 * when winsorization is not required. This streamlined version returns only the
 * vertex coefficients without diagnostic data, minimizing memory allocation and
 * computational overhead.
 *
 * The local correlation at each vertex v is computed as:
 *
 *   lcor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * where the sum is over edges incident to vertex v.
 *
 * COMPUTATIONAL COMPLEXITY:
 * Time: O(|E|) where |E| is the total number of edge-vertex incidences
 * Space: O(n) for the output vector only
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme:
 *   - UNIT: w_e = 1 (combinatorial, ignores edge lengths)
 *   - DERIVATIVE: w_e = 1/ℓ_e² (geometric, compares directional derivatives)
 *   - SIGN: w_e = 1 (falls back to unit weighting)
 * @param y_diff_type How to compute edge differences for y:
 *   - DIFFERENCE: Δ_e y = y(u) - y(v)
 *   - LOGRATIO: Δ_e y = log((y(u) + ε) / (y(v) + ε))
 * @param z_diff_type How to compute edge differences for z (same options as y)
 * @param epsilon_y Pseudocount for y when using LOGRATIO (pre-computed, not adaptive)
 * @param epsilon_z Pseudocount for z when using LOGRATIO (pre-computed, not adaptive)
 *
 * @return Vector of length n containing lcor(y,z)(v) for each vertex v.
 *         Values are in [-1, 1], with 0 returned when either function has
 *         insufficient variation at a vertex.
 *
 * @note This is an internal helper called by lcor(). For diagnostic output
 *       including edge differences and weights, use lcor_instrumented() instead.
 *
 * @see lcor() for the public entry point with input validation and adaptive epsilon
 * @see lcor_one_pass_instrumented() for the version with diagnostic output
 */
std::vector<double> set_wgraph_t::lcor_one_pass(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon_y,
    double epsilon_z
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;
    const double MIN_EDGE_LENGTH = 1e-10;

    std::vector<double> result(n_vertices, 0.0);

    // Compute coefficient at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;       // Σ w_e Δ_e y · Δ_e z
        double sum_y_squared = 0.0;   // Σ w_e (Δ_e y)²
        double sum_z_squared = 0.0;   // Σ w_e (Δ_e z)²

        // Iterate over neighbors of vertex v
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;

            // Compute edge differences based on type
            const double delta_y = compute_edge_diff(
                y[u], y[v], y_diff_type, epsilon_y
            );
            const double delta_z = compute_edge_diff(
                z[u], z[v], z_diff_type, epsilon_z
            );

            // Compute weight based on type
            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;

                case lcor_type_t::DERIVATIVE:
                    // w_e = 1/ℓ_e² for derivative-like weighting
                    if (edge_info.weight > MIN_EDGE_LENGTH) {
                        weight = 1.0 / (edge_info.weight * edge_info.weight);
                    } else {
                        // Skip degenerate edges (would have infinite weight)
                        continue;
                    }
                    break;

                case lcor_type_t::SIGN:
                    // Not applicable for correlation-type coefficients
                    // Fall back to unit weighting
                    weight = 1.0;
                    break;
            }

            // Accumulate weighted sums
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
        }

        // Compute denominators (square roots of sum of squares)
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);

        // Check numerical stability
        // If either function is essentially constant in this neighborhood,
        // the coefficient is undefined (set to 0 by convention)
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result[v] = numerator / (denom_y * denom_z);

            // Clamp to [-1, 1] to handle numerical errors
            // (mathematically should be in [-1,1], but floating-point...)
            if (result[v] > 1.0) {
                result[v] = 1.0;
            } else if (result[v] < -1.0) {
                result[v] = -1.0;
            }
        } else {
            // One or both functions essentially constant
            result[v] = 0.0;
        }
    }

    return result;
}

/**
 * @brief Compute local correlation coefficients (two-pass, with winsorization) - production version
 *
 * Two-pass implementation for computing local correlation coefficients with
 * winsorization for outlier robustness. This streamlined version returns only
 * the vertex coefficients without exposing intermediate diagnostic data.
 *
 * Pass 1: Compute all edge differences and determine winsorization bounds
 * Pass 2: Compute coefficients using clipped (winsorized) differences
 *
 * WINSORIZATION:
 * Edge differences are clipped to [q, 1-q] percentiles where q is the
 * winsorize_quantile. For example, q = 0.05 clips to the 5th and 95th
 * percentiles. This reduces the influence of extreme outliers while
 * preserving the bulk of the distribution.
 *
 * COMPUTATIONAL COMPLEXITY:
 * Time: O(|E| + |E|log|E|) due to sorting for percentile computation
 * Space: O(|E|) for storing edge differences during winsorization
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme (UNIT, DERIVATIVE, or SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for z (DIFFERENCE or LOGRATIO)
 * @param epsilon_y Pseudocount for y when using LOGRATIO (pre-computed)
 * @param epsilon_z Pseudocount for z when using LOGRATIO (pre-computed)
 * @param winsorize_quantile Quantile for clipping (e.g., 0.05 for 5%/95%)
 *
 * @return Vector of length n containing lcor(y,z)(v) for each vertex v.
 *         Values are in [-1, 1], computed using winsorized edge differences.
 *
 * @note This is an internal helper called by lcor(). For diagnostic output
 *       including raw edge differences and winsorization bounds, use
 *       lcor_instrumented() instead.
 *
 * @see lcor() for the public entry point with input validation
 * @see lcor_two_pass_instrumented() for the version with diagnostic output
 */
std::vector<double> set_wgraph_t::lcor_two_pass(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon_y,
    double epsilon_z,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;

    // Pass 1: Collect all edge differences
    std::vector<double> all_delta_y;
    std::vector<double> all_delta_z;

    // Estimate total number of edges for efficient allocation
    size_t total_edges = 0;
    for (size_t v = 0; v < n_vertices; ++v) {
        total_edges += adjacency_list[v].size();
    }
    all_delta_y.reserve(total_edges);
    all_delta_z.reserve(total_edges);

    // Store edge differences and weights for each vertex
    // Using vectors-of-vectors to maintain vertex-neighborhood structure
    std::vector<std::vector<double>> vertex_delta_y(n_vertices);
    std::vector<std::vector<double>> vertex_delta_z(n_vertices);
    std::vector<std::vector<double>> vertex_weights(n_vertices);

    for (size_t v = 0; v < n_vertices; ++v) {
        const size_t n_neighbors = adjacency_list[v].size();
        vertex_delta_y[v].reserve(n_neighbors);
        vertex_delta_z[v].reserve(n_neighbors);
        vertex_weights[v].reserve(n_neighbors);

        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;

            // Compute weight
            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;

                case lcor_type_t::DERIVATIVE:
                    if (edge_info.weight > 1e-10) {
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
            const double delta_y = compute_edge_diff(
                y[u], y[v], y_diff_type, epsilon_y
            );
            const double delta_z = compute_edge_diff(
                z[u], z[v], z_diff_type, epsilon_z
            );

            // Store for this vertex
            vertex_delta_y[v].push_back(delta_y);
            vertex_delta_z[v].push_back(delta_z);
            vertex_weights[v].push_back(weight);

            // Collect for global winsorization bounds
            all_delta_y.push_back(delta_y);
            all_delta_z.push_back(delta_z);
        }
    }

    // Compute winsorization bounds across all edges
    auto [y_lower, y_upper] = compute_winsorize_bounds(
        all_delta_y, winsorize_quantile
    );
    auto [z_lower, z_upper] = compute_winsorize_bounds(
        all_delta_z, winsorize_quantile
    );

    // Pass 2: Compute coefficients with winsorized differences
    std::vector<double>  result(n_vertices, 0.0);

    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;
        double sum_y_squared = 0.0;
        double sum_z_squared = 0.0;

        const size_t n_edges = vertex_delta_y[v].size();
        for (size_t i = 0; i < n_edges; ++i) {
            // Apply winsorization clipping
            const double delta_y = winsorize_clip(
                vertex_delta_y[v][i], y_lower, y_upper
            );
            const double delta_z = winsorize_clip(
                vertex_delta_z[v][i], z_lower, z_upper
            );
            const double weight = vertex_weights[v][i];

            // Accumulate weighted sums
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
        }

        // Compute denominators
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);

        // Compute coefficient with numerical stability check
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result[v] = numerator / (denom_y * denom_z);

            // Clamp to [-1, 1]
            if (result[v] > 1.0) {
                result[v] = 1.0;
            } else if (result[v] < -1.0) {
                result[v] = -1.0;
            }
        } else {
            result[v] = 0.0;
        }
    }

    return result;
}

/**
 * @brief Compute local correlation coefficients between two functions on a graph
 *
 * Production entry point for computing vertex-level correlation coefficients
 * measuring the alignment of directional changes between functions y and z.
 * This streamlined version returns only the coefficient vector, making it
 * suitable for production pipelines where diagnostic data is not needed.
 *
 * The local correlation coefficient at each vertex v is:
 *
 *   lcor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * where the sum is over edges incident to vertex v, and Δ_e f represents the
 * directional difference along edge e = [v,u].
 *
 * EDGE DIFFERENCE TYPES:
 * - DIFFERENCE: Δ_e f = f(u) - f(v) for continuous data
 * - LOGRATIO: Δ_e f = log((f(u) + ε)/(f(v) + ε)) for compositional data
 *
 * WEIGHTING SCHEMES:
 * - UNIT: w_e = 1, treats all edges equally
 * - DERIVATIVE: w_e = 1/ℓ_e², normalizes by edge length for geometric interpretation
 *
 * GEOMETRIC INTERPRETATION:
 * With derivative weighting, the coefficient converges to cos(θ) where θ is
 * the angle between gradient vectors. Values near +1 indicate parallel gradients,
 * near -1 indicate anti-parallel gradients, and near 0 indicate orthogonality.
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param weight_type Weighting scheme (UNIT, DERIVATIVE, or SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for z (DIFFERENCE or LOGRATIO)
 * @param epsilon Pseudocount for log-ratio transformations (default 0 = adaptive).
 *   If 0, computed as 1e-6 * min(non-zero values) for each function.
 * @param winsorize_quantile Quantile for winsorization (default 0 = none).
 *   If > 0, uses two-pass algorithm with clipping to [q, 1-q] percentiles.
 *
 * @return Vector of length n containing lcor(y,z)(v) for each vertex v.
 *         Values are in [-1, 1].
 *
 * @throws REPORT_ERROR if y or z length does not match num_vertices
 *
 * @note For diagnostic output including edge differences, weights, and
 *       winsorization bounds, use lcor_instrumented() instead.
 *
 * @see lcor_instrumented() for the version with full diagnostic output
 * @see comono_cor() for the standard difference-only version
 */
std::vector<double> set_wgraph_t::lcor(
    const std::vector<double>& y,
    const std::vector<double>& z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();

    // Validate input lengths
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (z.size() != n_vertices) {
        REPORT_ERROR("Length of z (%zu) does not match number of vertices (%zu)\n",
                     z.size(), n_vertices);
    }

    // Compute adaptive epsilon if needed (epsilon = 0 means adaptive)
    double epsilon_y = epsilon;
    double epsilon_z = epsilon;

    if (y_diff_type == edge_diff_type_t::LOGRATIO && epsilon_y <= 0.0) {
        epsilon_y = compute_adaptive_epsilon(y);
    }
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0) {
        epsilon_z = compute_adaptive_epsilon(z);
    }

    // Choose implementation based on whether winsorization is requested
    const bool need_winsorization = (winsorize_quantile > 0.0 &&
                                     winsorize_quantile < 0.5);

    if (need_winsorization) {
        return lcor_two_pass(
            y, z,
            weight_type, y_diff_type, z_diff_type,
            epsilon_y, epsilon_z, winsorize_quantile
        );
    } else {
        return lcor_one_pass(
            y, z,
            weight_type, y_diff_type, z_diff_type,
            epsilon_y, epsilon_z
        );
    }
}

/**
 * @brief Compute local correlation between a vector y and each column of matrix Z
 *
 * This method efficiently computes vertex-level correlation coefficients between
 * a single response function y and multiple feature functions stored as columns
 * of matrix Z. The implementation pre-computes y-dependent quantities once and
 * reuses them across all columns.
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param Z Feature matrix (rows = num_vertices, columns = number of features)
 * @param weight_type Weighting scheme (UNIT, DERIVATIVE, or SIGN)
 * @param y_diff_type Edge difference type for y (DIFFERENCE or LOGRATIO)
 * @param z_diff_type Edge difference type for Z columns (DIFFERENCE or LOGRATIO)
 * @param epsilon Pseudocount for log-ratio transformations (0 = adaptive)
 * @param winsorize_quantile Quantile for winsorization (0 = none)
 * @return lcor_vector_matrix_result_t with coefficient matrix and bounds
 */
lcor_vector_matrix_result_t set_wgraph_t::lcor_vector_matrix(
    const std::vector<double>& y,
    const Eigen::MatrixXd& Z,
    lcor_type_t weight_type,
    edge_diff_type_t y_diff_type,
    edge_diff_type_t z_diff_type,
    double epsilon,
    double winsorize_quantile
) const {
    const size_t n_vertices = num_vertices();
    const size_t n_columns = static_cast<size_t>(Z.cols());
    const double MIN_DENOMINATOR = 1e-10;

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
        epsilon_y = compute_adaptive_epsilon(y);
    }

    double epsilon_z = epsilon;
    if (z_diff_type == edge_diff_type_t::LOGRATIO && epsilon_z <= 0.0 && n_columns > 0) {
        double min_nonzero = std::numeric_limits<double>::max();
        for (size_t col = 0; col < n_columns; ++col) {
            for (size_t i = 0; i < n_vertices; ++i) {
                double val = Z(i, col);
                if (val > 0.0 && val < min_nonzero) {
                    min_nonzero = val;
                }
            }
        }
        epsilon_z = (min_nonzero < std::numeric_limits<double>::max())
            ? 1e-6 * min_nonzero
            : 1e-6;
    }

    // ---- Initialize result ----
    lcor_vector_matrix_result_t result;
    result.coefficients.resize(n_vertices, n_columns);
    result.coefficients.setZero();
    result.z_lower.resize(n_columns, -std::numeric_limits<double>::max());
    result.z_upper.resize(n_columns, std::numeric_limits<double>::max());

    // ---- Phase 1: Pre-compute y-dependent quantities ----

    std::vector<std::vector<size_t>> neighbor_indices(n_vertices);
    std::vector<std::vector<double>> y_edge_diffs(n_vertices);
    std::vector<std::vector<double>> edge_weights(n_vertices);

    std::vector<double> all_delta_y;
    const bool need_winsorization = (winsorize_quantile > 0.0 &&
                                     winsorize_quantile < 0.5);

    // Estimate total edges for efficient allocation
    size_t total_edges = 0;
    for (size_t v = 0; v < n_vertices; ++v) {
        total_edges += adjacency_list[v].size();
    }
    if (need_winsorization) {
        all_delta_y.reserve(total_edges);
    }

    // Pre-compute y quantities at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        const size_t n_neighbors = adjacency_list[v].size();
        neighbor_indices[v].reserve(n_neighbors);
        y_edge_diffs[v].reserve(n_neighbors);
        edge_weights[v].reserve(n_neighbors);

        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;
            const double edge_length = edge_info.weight;

            neighbor_indices[v].push_back(u);

            const double delta_y = compute_edge_diff(
                y[u], y[v], y_diff_type, epsilon_y
            );
            y_edge_diffs[v].push_back(delta_y);

            if (need_winsorization) {
                all_delta_y.push_back(delta_y);
            }

            double weight = 1.0;
            switch (weight_type) {
                case lcor_type_t::UNIT:
                    weight = 1.0;
                    break;
                case lcor_type_t::DERIVATIVE:
                    weight = (edge_length > 1e-10)
                        ? 1.0 / (edge_length * edge_length)
                        : 0.0;
                    break;
                case lcor_type_t::SIGN:
                    weight = 1.0;
                    break;
            }
            edge_weights[v].push_back(weight);
        }
    }

    // Compute y winsorization bounds
    double y_lower = -std::numeric_limits<double>::max();
    double y_upper = std::numeric_limits<double>::max();
    if (need_winsorization) {
        auto bounds = compute_winsorize_bounds(all_delta_y, winsorize_quantile);
        y_lower = bounds.first;
        y_upper = bounds.second;
    }
    result.y_lower = y_lower;
    result.y_upper = y_upper;

    // Pre-compute y normalization factors at each vertex
    std::vector<double> denom_y(n_vertices, 0.0);
    for (size_t v = 0; v < n_vertices; ++v) {
        double sum_y_squared = 0.0;
        const size_t n_neighbors = neighbor_indices[v].size();

        for (size_t i = 0; i < n_neighbors; ++i) {
            double delta_y = y_edge_diffs[v][i];
            if (need_winsorization) {
                delta_y = winsorize_clip(delta_y, y_lower, y_upper);
            }
            const double weight = edge_weights[v][i];
            sum_y_squared += weight * delta_y * delta_y;
        }
        denom_y[v] = std::sqrt(sum_y_squared);
    }

    // ---- Phase 2: Process each column of Z ----
    for (size_t col = 0; col < n_columns; ++col) {

        // If winsorizing, collect z edge differences for this column
        std::vector<double> all_delta_z;
        double z_lower_col = -std::numeric_limits<double>::max();
        double z_upper_col = std::numeric_limits<double>::max();

        if (need_winsorization) {
            all_delta_z.reserve(total_edges);
            for (size_t v = 0; v < n_vertices; ++v) {
                for (size_t i = 0; i < neighbor_indices[v].size(); ++i) {
                    const size_t u = neighbor_indices[v][i];
                    const double delta_z = compute_edge_diff(
                        Z(u, col), Z(v, col), z_diff_type, epsilon_z
                    );
                    all_delta_z.push_back(delta_z);
                }
            }
            auto bounds = compute_winsorize_bounds(all_delta_z, winsorize_quantile);
            z_lower_col = bounds.first;
            z_upper_col = bounds.second;
            result.z_lower[col] = z_lower_col;
            result.z_upper[col] = z_upper_col;
        }

        // Compute coefficient at each vertex for this column
        for (size_t v = 0; v < n_vertices; ++v) {
            const size_t n_neighbors = neighbor_indices[v].size();

            if (n_neighbors == 0 || denom_y[v] < MIN_DENOMINATOR) {
                result.coefficients(v, col) = 0.0;
                continue;
            }

            double numerator = 0.0;
            double sum_z_squared = 0.0;

            for (size_t i = 0; i < n_neighbors; ++i) {
                const size_t u = neighbor_indices[v][i];
                const double weight = edge_weights[v][i];

                double delta_y = y_edge_diffs[v][i];
                if (need_winsorization) {
                    delta_y = winsorize_clip(delta_y, y_lower, y_upper);
                }

                double delta_z = compute_edge_diff(
                    Z(u, col), Z(v, col), z_diff_type, epsilon_z
                );
                if (need_winsorization) {
                    delta_z = winsorize_clip(delta_z, z_lower_col, z_upper_col);
                }

                numerator += weight * delta_y * delta_z;
                sum_z_squared += weight * delta_z * delta_z;
            }

            const double denom_z = std::sqrt(sum_z_squared);
            if (denom_z > MIN_DENOMINATOR) {
                double coeff = numerator / (denom_y[v] * denom_z);
                // Clamp to [-1, 1]
                if (coeff > 1.0) coeff = 1.0;
                else if (coeff < -1.0) coeff = -1.0;
                result.coefficients(v, col) = coeff;
            } else {
                result.coefficients(v, col) = 0.0;
            }
        }
    }

    return result;
}
