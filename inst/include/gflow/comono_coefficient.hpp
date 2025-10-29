#ifndef COMONO_COEFFICIENT_HPP
#define COMONO_COEFFICIENT_HPP

/**
 * @file comono_coefficient.hpp
 * @brief Co-monotonicity coefficient types and structures
 *
 * This file defines the types and structures for computing co-monotonicity
 * coefficients between functions defined on graph vertices. These measures
 * quantify the agreement in directional changes across graph edges.
 */

/**
 * @enum comono_type_t
 * @brief Types of co-monotonicity coefficients
 *
 * Defines the different weighting schemes for computing co-monotonicity between
 * two functions y and z defined on graph vertices. Each type represents a
 * different way of aggregating the edge-wise co-variation.
 */
enum class comono_type_t {
    /**
     * Unit weights (w_e = 1 for all edges)
     *
     * The co-monotonicity at vertex v is:
     *   comono(y,z)(v) = Σ(Δ_e y · Δ_e z) / Σ|Δ_e y · Δ_e z|
     * where the sum is over all edges incident to v, and Δ_e y = y(u) - y(v)
     * for edge e = [v,u].
     */
    UNIT,

    /**
     * Derivative-like weights (w_e = 1/(Δ_e)²)
     *
     * The co-monotonicity at vertex v is:
     *   comono_∂(y,z)(v) = Σ(w_e · Δ_e y · Δ_e z) / Σ(w_e · |Δ_e y · Δ_e z|)
     * where w_e = 1/(Δ_e)² and Δ_e is the edge length. This weighting normalizes
     * by edge length squared, making it analogous to comparing derivatives rather
     * than absolute changes.
     */
    DERIVATIVE,

    /**
     * Sign-based co-monotonicity (using only direction of change)
     *
     * The co-monotonicity at vertex v is:
     *   comono_±(y,z)(v) = Σ sign(Δ_e y · Δ_e z) / |N(v)|
     * where sign(x) = +1 if x > 0, -1 if x < 0, and 0 if x = 0.
     * This measure only captures whether the functions change in the same
     * direction, ignoring the magnitude of change.
     */
    SIGN
};

/**
 * @struct comono_result_t
 * @brief Result structure for co-monotonicity computation
 *
 * Contains the vertex-wise co-monotonicity coefficients along with summary
 * statistics and diagnostic information.
 */
struct comono_result_t {
    std::vector<double> vertex_coefficients;  ///< Co-monotonicity coefficient at each vertex
    double mean_coefficient;                  ///< Mean co-monotonicity across all vertices
    double median_coefficient;                ///< Median co-monotonicity across all vertices
    size_t n_positive;                        ///< Number of vertices with positive co-monotonicity
    size_t n_negative;                        ///< Number of vertices with negative co-monotonicity
    size_t n_zero;                            ///< Number of vertices with zero co-monotonicity

    /**
     * @brief Default constructor
     */
    comono_result_t()
        : mean_coefficient(0.0),
          median_coefficient(0.0),
          n_positive(0),
          n_negative(0),
          n_zero(0)
    {}
};

struct comono_matrix_result_t {
    std::vector<std::vector<double>> column_coefficients;  // [col][vertex]
    std::vector<double> mean_coefficients;      // One per column
    std::vector<double> median_coefficients;    // One per column
    std::vector<size_t> n_positive;             // One per column
    std::vector<size_t> n_negative;             // One per column
    std::vector<size_t> n_zero;                 // One per column
    size_t n_vertices;
    size_t n_columns;
};

#endif // COMONO_COEFFICIENT_HPP
