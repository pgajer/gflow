/**
 * @file comono_coefficient.cpp
 * @brief Implementation of co-monotonicity coefficient computation
 *
 * This file implements methods for computing co-monotonicity coefficients
 * between two functions defined on graph vertices. Co-monotonicity measures
 * quantify the extent to which two functions vary together across the edges
 * of a graph, providing a graph-based analogue of correlation that respects
 * the local geometric structure.
 *
 * MOTIVATION
 *
 * Consider two real-valued functions y and z defined on the vertices of a
 * weighted graph G = (V, E). In classical statistics, we might measure their
 * association using Pearson correlation, which treats all pairs of observations
 * as equally related. However, when the data possess inherent geometric structure
 * encoded by the graph, we can define more nuanced measures that respect the
 * local neighborhood relationships.
 *
 * The co-monotonicity coefficient addresses the question: do the functions y and z
 * tend to increase or decrease together as we move along edges of the graph? This
 * question is natural in spatial statistics, where we expect nearby observations
 * to be more strongly related, and in network analysis, where the graph structure
 * encodes meaningful relationships.
 *
 * MATHEMATICAL FRAMEWORK
 *
 * For an edge e = [v,u] ∈ E connecting vertices v and u, we define the edge
 * difference operator Δ_e acting on a function f: V → ℝ by
 *
 *   Δ_e f = f(u) - f(v).
 *
 * The product Δ_e y · Δ_e z captures the co-variation of y and z along edge e:
 * it is positive when both functions increase or both decrease, negative when
 * they change in opposite directions, and zero when at least one function
 * remains constant.
 *
 * VERTEX-LEVEL CO-MONOTONICITY
 *
 * At each vertex v, we aggregate the edge-wise co-variation over its neighborhood
 * N(v) to obtain a local co-monotonicity measure. The general form is
 *
 *   comono(y,z;w)(v) = Σ_{u ∈ N(v)} w_e Δ_e y Δ_e z / Σ_{u ∈ N(v)} w_e |Δ_e y Δ_e z|,
 *
 * where w_e ≥ 0 are edge weights and the denominator ensures the coefficient
 * lies in [-1, 1]. The value equals +1 when y and z always change in the same
 * direction within the neighborhood, -1 when they always change in opposite
 * directions, and values near 0 indicate little systematic relationship.
 *
 * THREE WEIGHTING SCHEMES
 *
 * We implement three standard weighting schemes, each suited to different contexts:
 *
 * 1. Unit weights (w_e = 1): Treats all edges equally, appropriate when edge
 *    lengths are roughly comparable or when we wish to count directional
 *    agreements without geometric normalization.
 *
 * 2. Derivative weights (w_e = 1/(Δ_e)²): Normalizes by edge length squared,
 *    making the measure analogous to comparing derivatives rather than absolute
 *    changes. This is natural when the functions represent continuous quantities
 *    sampled at irregular spatial positions.
 *
 * 3. Sign-based weights: Uses only the sign of Δ_e y · Δ_e z, counting the
 *    proportion of edges where the functions agree in direction. This robust
 *    measure is insensitive to outliers and magnitude of change.
 *
 * COMPUTATIONAL CONSIDERATIONS
 *
 * The algorithm iterates over each vertex and its incident edges, computing the
 * weighted sum of co-variations. The time complexity is O(|E|), linear in the
 * number of edges. For sparse graphs with bounded degree, this is also O(|V|).
 *
 * Special care is taken to handle numerical edge cases:
 * - When the denominator is zero at a vertex (all products Δ_e y · Δ_e z = 0),
 *   the coefficient is defined as 0.
 * - For the derivative weighting, we guard against division by zero for edges
 *   with negligible length.
 */

#include "set_wgraph.hpp"
#include "comono_coefficient.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

/**
 * @brief Compute co-monotonicity coefficients between two functions on graph vertices
 *
 * Given two real-valued functions y and z defined on the vertices of a weighted
 * graph, this method computes the vertex-wise co-monotonicity coefficient, which
 * measures the extent to which the two functions vary together across graph edges.
 *
 * The co-monotonicity coefficient at vertex v quantifies the agreement in
 * directional changes of y and z within the neighborhood of v. A coefficient
 * near +1 indicates that y and z tend to increase or decrease together, while
 * a coefficient near -1 indicates they change in opposite directions. Values
 * near 0 suggest no systematic relationship.
 *
 * WEIGHTING SCHEMES
 *
 * The method supports three weighting schemes controlled by the type parameter:
 *
 * - UNIT: Uses unit weights (w_e = 1) for all edges. This treats all edges
 *   equally regardless of their length.
 *
 * - DERIVATIVE: Uses weights w_e = 1/(Δ_e)² where Δ_e is the edge length.
 *   This normalizes by edge length squared, making the measure analogous to
 *   comparing derivatives of y and z.
 *
 * - SIGN: Uses only the sign of the product Δ_e y · Δ_e z, computing the
 *   proportion of edges where the functions agree in direction. This is a
 *   robust measure insensitive to outliers.
 *
 * INTERPRETATION
 *
 * The vertex-wise coefficients provide a spatially-resolved view of co-monotonicity,
 * revealing where in the graph the two functions are most strongly (or weakly)
 * associated. The summary statistics (mean, median, proportion positive/negative)
 * characterize the global pattern of association.
 *
 * USAGE EXAMPLE
 *
 * Suppose we have fitted a regression model ŷ to observed data y on a spatial
 * graph. We might compute comono(y, ŷ) to assess where the model captures the
 * local directional behavior of y. Vertices with low co-monotonicity indicate
 * regions where the model may poorly represent local trends.
 *
 * @param y First function values at vertices (length must equal num_vertices())
 * @param z Second function values at vertices (length must equal num_vertices())
 * @param type Type of co-monotonicity coefficient (UNIT, DERIVATIVE, or SIGN)
 *
 * @return comono_result_t structure containing:
 *   - vertex_coefficients: Co-monotonicity at each vertex
 *   - mean_coefficient: Mean over all vertices
 *   - median_coefficient: Median over all vertices
 *   - n_positive: Count of vertices with positive co-monotonicity
 *   - n_negative: Count of vertices with negative co-monotonicity
 *   - n_zero: Count of vertices with zero co-monotonicity
 *
 * @throws std::invalid_argument if y or z has incorrect length
 *
 * @note Time complexity: O(|E|) where |E| is the number of edges
 * @note Space complexity: O(|V|) for storing vertex coefficients
 *
 * @see comono_type_t for detailed description of weighting schemes
 */
comono_result_t set_wgraph_t::comono(
    const std::vector<double>& y,
    const std::vector<double>& z,
    comono_type_t type
) const {
    const size_t n_vertices = num_vertices();

    // Validate input lengths
    if (y.size() != n_vertices) {
        throw std::invalid_argument(
            "Length of y (" + std::to_string(y.size()) +
            ") does not match number of vertices (" + std::to_string(n_vertices) + ")"
        );
    }
    if (z.size() != n_vertices) {
        throw std::invalid_argument(
            "Length of z (" + std::to_string(z.size()) +
            ") does not match number of vertices (" + std::to_string(n_vertices) + ")"
        );
    }

    // Initialize result structure
    comono_result_t result;
    result.vertex_coefficients.resize(n_vertices, 0.0);

    // Compute co-monotonicity at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        double weighted_sum = 0.0;      // Numerator: Σ w_e Δ_e y Δ_e z
        double weighted_abs_sum = 0.0;  // Denominator: Σ w_e |Δ_e y Δ_e z|
        size_t n_neighbors = 0;         // Neighborhood size (for SIGN type)

        // Iterate over neighbors of vertex v
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;
            const double edge_length = edge_info.weight;

            // Compute edge differences
            const double delta_y = y[u] - y[v];
            const double delta_z = z[u] - z[v];
            const double product = delta_y * delta_z;

            // Compute weight based on type
            double weight = 1.0;
            switch (type) {
                case comono_type_t::UNIT:
                    weight = 1.0;
                    break;

                case comono_type_t::DERIVATIVE:
                    // Guard against division by zero
                    if (edge_length > 1e-10) {
                        weight = 1.0 / (edge_length * edge_length);
                    } else {
                        weight = 0.0;
                    }
                    break;

                case comono_type_t::SIGN:
                    // For sign-based measure, we count agreements
                    if (product > 0.0) {
                        weighted_sum += 1.0;
                    } else if (product < 0.0) {
                        weighted_sum -= 1.0;
                    }
                    ++n_neighbors;
                    continue;  // Skip the general weight accumulation below
            }

            // Accumulate weighted sums (for UNIT and DERIVATIVE types)
            weighted_sum += weight * product;
            weighted_abs_sum += weight * std::abs(product);
        }

        // Compute vertex-level coefficient
        if (type == comono_type_t::SIGN) {
            // For sign-based: proportion of agreeing edges
            result.vertex_coefficients[v] = (n_neighbors > 0) ?
                weighted_sum / static_cast<double>(n_neighbors) : 0.0;
        } else {
            // For UNIT and DERIVATIVE: ratio of weighted sums
            result.vertex_coefficients[v] = (weighted_abs_sum > 1e-10) ?
                weighted_sum / weighted_abs_sum : 0.0;
        }
    }

    // Compute summary statistics
    result.mean_coefficient = std::accumulate(
        result.vertex_coefficients.begin(),
        result.vertex_coefficients.end(),
        0.0
    ) / n_vertices;

    // Compute median
    std::vector<double> sorted_coeffs = result.vertex_coefficients;
    std::sort(sorted_coeffs.begin(), sorted_coeffs.end());
    if (n_vertices % 2 == 0) {
        result.median_coefficient = 0.5 * (
            sorted_coeffs[n_vertices / 2 - 1] + sorted_coeffs[n_vertices / 2]
        );
    } else {
        result.median_coefficient = sorted_coeffs[n_vertices / 2];
    }

    // Count positive, negative, and zero coefficients
    const double epsilon = 1e-10;
    result.n_positive = 0;
    result.n_negative = 0;
    result.n_zero = 0;

    for (double coeff : result.vertex_coefficients) {
        if (coeff > epsilon) {
            ++result.n_positive;
        } else if (coeff < -epsilon) {
            ++result.n_negative;
        } else {
            ++result.n_zero;
        }
    }

    return result;
}

/**
 * @brief Compute global co-monotonicity coefficient (graph-wide summary)
 *
 * This convenience method computes a single scalar co-monotonicity coefficient
 * summarizing the association between y and z across the entire graph. It is
 * defined as the mean of the vertex-wise coefficients.
 *
 * For a quick assessment of overall co-monotonicity, use this method. For
 * spatially-resolved analysis, use the full comono() method which returns
 * vertex-wise coefficients.
 *
 * @param y First function values at vertices
 * @param z Second function values at vertices
 * @param type Type of co-monotonicity coefficient
 *
 * @return Scalar co-monotonicity coefficient in [-1, 1]
 *
 * @throws std::invalid_argument if y or z has incorrect length
 */
double set_wgraph_t::comono_global(
    const std::vector<double>& y,
    const std::vector<double>& z,
    comono_type_t type
) const {
    comono_result_t result = comono(y, z, type);
    return result.mean_coefficient;
}

/**
 * @brief Compute co-monotonicity coefficients between a vector and matrix columns
 *
 * @details
 * This function efficiently computes vertex-wise co-monotonicity coefficients
 * between a vector y and each column of a matrix Z. The implementation optimizes
 * performance by pre-computing graph-dependent quantities that remain constant
 * across all columns, then reusing them for each column computation.
 *
 * MOTIVATION
 *
 * In multivariate statistical analysis on graphs, we often need to assess the
 * association between a scalar response y and multiple predictors z_1, ..., z_q.
 * Examples include:
 * - Feature selection: identifying which predictors are most concordant with y
 * - Model comparison: evaluating multiple candidate predictions simultaneously
 * - Multivariate regression diagnostics: checking co-monotonicity of residuals
 *   with multiple covariates
 * - Time series analysis: assessing temporal concordance across multiple series
 *
 * Computing these associations independently using repeated calls to comono()
 * would redundantly recompute graph structure information q times. This matrix
 * version eliminates redundancy by performing setup once and reusing it.
 *
 * ALGORITHM
 *
 * The algorithm proceeds in two phases:
 *
 * Phase 1 - Pre-computation (performed once):
 *   For each vertex v and each incident edge e = [v,u]:
 *     - Compute Δ_e y = y(u) - y(v)
 *     - Compute edge weight w_e based on type
 *     - Store neighbor index u
 *
 * Phase 2 - Column processing (repeated for each column):
 *   For each column j of Z:
 *     For each vertex v:
 *       For each incident edge e = [v,u]:
 *         - Compute Δ_e z_j = z_j(u) - z_j(v)
 *         - Form product with pre-computed Δ_e y
 *         - Accumulate using pre-computed weights
 *       Compute vertex coefficient from accumulated values
 *     Compute summary statistics for column j
 *
 * This structure ensures that expensive graph traversal and weight computation
 * occur only once, while column-specific computations reuse the cached information.
 *
 * WEIGHTING SCHEMES
 *
 * Three weighting schemes are supported, identical to the vector version:
 *
 * UNIT (w_e = 1): Treats all edges equally, counting directional agreements
 * without geometric normalization. Appropriate when edges are comparable in
 * length or when emphasizing combinatorial structure over geometry.
 *
 * DERIVATIVE (w_e = 1/Δ_e²): Normalizes by squared edge length, making the
 * measure analogous to comparing derivatives. Natural for continuous functions
 * sampled at irregular positions where absolute changes should be scaled by
 * spatial separation.
 *
 * SIGN: Uses only direction of change (sign of Δ_e y · Δ_e z_j), computing
 * proportion of edges where functions agree. Robust to outliers and insensitive
 * to magnitude, useful when only ordinal relationships matter.
 *
 * INTERPRETATION
 *
 * The results provide both fine-grained (vertex-wise) and aggregated (mean, median)
 * views of co-monotonicity for each column. This multi-scale perspective enables:
 * - Global assessment via mean coefficients across columns
 * - Local diagnostics via vertex coefficients for specific columns
 * - Comparative analysis across multiple predictors
 * - Identification of spatial patterns in association strength
 *
 * For feature selection, columns with high mean absolute co-monotonicity indicate
 * strong concordance with y. For model assessment, comparing co-monotonicity of
 * multiple models' predictions reveals which best captures y's directional behavior.
 *
 * NUMERICAL CONSIDERATIONS
 *
 * The implementation guards against numerical issues:
 * - Division by zero: denominators checked against threshold (1e-10)
 * - Edge length validation: derivative weights computed only for edges with
 *   length exceeding threshold
 * - Epsilon comparison: classifications (positive/negative/zero) use tolerance
 *   to handle floating-point imprecision
 *
 * For very large matrices, consider processing columns in batches to manage
 * memory usage, as results require O(q · |V|) storage.
 *
 * @param y Vector of function values at each vertex, length must equal num_vertices()
 * @param Z Matrix of function values, each column is a function on vertices,
 *          rows must equal num_vertices()
 * @param type Weighting scheme: UNIT, DERIVATIVE, or SIGN
 *
 * @return comono_matrix_result_t containing:
 *   - column_coefficients: vector of vectors, [j][v] gives coefficient at vertex v
 *                          for column j
 *   - mean_coefficients: vector of length q, mean coefficient for each column
 *   - median_coefficients: vector of length q, median coefficient for each column
 *   - n_positive: vector of length q, count of positive vertices per column
 *   - n_negative: vector of length q, count of negative vertices per column
 *   - n_zero: vector of length q, count of zero vertices per column
 *   - n_vertices: number of vertices in graph
 *   - n_columns: number of columns in Z (q)
 *
 * @throws std::invalid_argument if y.size() ≠ num_vertices()
 * @throws std::invalid_argument if Z.rows() ≠ num_vertices()
 *
 * @note Time complexity: O(q · |E| + q · |V| log |V|) where q is number of columns
 * @note Space complexity: O(|E| + q · |V|)
 *
 * @see comono for single vector-vector co-monotonicity
 * @see comono_type_t for weighting scheme details
 * @see comono_matrix_result_t for result structure
 *
 * @example
 * // Assess multiple model predictions
 * std::vector<double> y = observed_values;
 * Eigen::MatrixXd predictions(n, 3);  // 3 models
 * predictions.col(0) = model1_predictions;
 * predictions.col(1) = model2_predictions;
 * predictions.col(2) = model3_predictions;
 *
 * comono_matrix_result_t result = graph.comono_matrix(y, predictions,
 *                                                      comono_type_t::UNIT);
 *
 * // Compare models by mean co-monotonicity
 * for (size_t j = 0; j < 3; ++j) {
 *     std::cout << "Model " << j << " mean comono: "
 *               << result.mean_coefficients[j] << std::endl;
 * }
 *
 * // Identify vertices where model 0 performs poorly
 * for (size_t v = 0; v < n; ++v) {
 *     if (result.column_coefficients[0][v] < 0.3) {
 *         std::cout << "Poor fit at vertex " << v << std::endl;
 *     }
 * }
 */
comono_matrix_result_t set_wgraph_t::comono_matrix(
    const std::vector<double>& y,
    const Eigen::MatrixXd& Z,
    comono_type_t type
) const {
    const size_t n_vertices = num_vertices();
    const size_t n_columns = Z.cols();

    // Validate
    if (y.size() != n_vertices) {
        throw std::invalid_argument("y size mismatch");
    }
    if (Z.rows() != static_cast<Eigen::Index>(n_vertices)) {
        throw std::invalid_argument("Z rows mismatch");
    }

    // Initialize result
    comono_matrix_result_t result;
    result.column_coefficients.resize(n_columns,
                                      std::vector<double>(n_vertices, 0.0));
    result.n_vertices = n_vertices;
    result.n_columns = n_columns;

    // Pre-compute y edge differences and weights for each vertex
    // This avoids redundant computation across columns
    std::vector<std::vector<double>> y_edge_diffs(n_vertices);
    std::vector<std::vector<double>> edge_weights(n_vertices);
    std::vector<std::vector<size_t>> neighbor_indices(n_vertices);

    for (size_t v = 0; v < n_vertices; ++v) {
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;
            const double edge_length = edge_info.weight;

            // Store y difference
            y_edge_diffs[v].push_back(y[u] - y[v]);

            // Store neighbor index
            neighbor_indices[v].push_back(u);

            // Compute and store weight based on type
            double weight = 1.0;
            if (type == comono_type_t::DERIVATIVE) {
                weight = (edge_length > 1e-10) ? 1.0 / (edge_length * edge_length) : 0.0;
            }
            edge_weights[v].push_back(weight);
        }
    }

    // Process each column of Z
    for (size_t col = 0; col < n_columns; ++col) {
        std::vector<double> coeffs(n_vertices, 0.0);

        // Compute coefficient for each vertex
        for (size_t v = 0; v < n_vertices; ++v) {
            double weighted_sum = 0.0;
            double weighted_abs_sum = 0.0;
            size_t n_neighbors = neighbor_indices[v].size();

            for (size_t i = 0; i < n_neighbors; ++i) {
                const size_t u = neighbor_indices[v][i];
                const double delta_y = y_edge_diffs[v][i];
                const double delta_z = Z(u, col) - Z(v, col);
                const double product = delta_y * delta_z;

                if (type == comono_type_t::SIGN) {
                    if (product > 0.0) {
                        weighted_sum += 1.0;
                    } else if (product < 0.0) {
                        weighted_sum -= 1.0;
                    }
                } else {
                    const double weight = edge_weights[v][i];
                    weighted_sum += weight * product;
                    weighted_abs_sum += weight * std::abs(product);
                }
            }

            // Compute vertex coefficient
            if (type == comono_type_t::SIGN) {
                coeffs[v] = (n_neighbors > 0) ?
                    weighted_sum / static_cast<double>(n_neighbors) : 0.0;
            } else {
                coeffs[v] = (weighted_abs_sum > 1e-10) ?
                    weighted_sum / weighted_abs_sum : 0.0;
            }
        }

        result.column_coefficients[col] = coeffs;

        // Compute summary statistics for this column
        double mean = std::accumulate(coeffs.begin(), coeffs.end(), 0.0) / n_vertices;
        result.mean_coefficients.push_back(mean);

        // Median
        std::vector<double> sorted = coeffs;
        std::sort(sorted.begin(), sorted.end());
        double median = (n_vertices % 2 == 0) ?
            0.5 * (sorted[n_vertices/2 - 1] + sorted[n_vertices/2]) :
            sorted[n_vertices/2];
        result.median_coefficients.push_back(median);

        // Count positive/negative/zero
        size_t n_pos = 0, n_neg = 0, n_zero = 0;
        const double eps = 1e-10;
        for (double c : coeffs) {
            if (c > eps) ++n_pos;
            else if (c < -eps) ++n_neg;
            else ++n_zero;
        }
        result.n_positive.push_back(n_pos);
        result.n_negative.push_back(n_neg);
        result.n_zero.push_back(n_zero);
    }

    return result;
}
