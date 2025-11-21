#include "set_wgraph.hpp"
#include "comono_coefficient.hpp"
#include "error_utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>

#define DEBUG_COMONO_COR 1  // Set to 1 to enable correlation debugging

/**
 * @brief Compute correlation-type co-monotonicity coefficients
 *
 * Computes vertex-level co-monotonicity using Pearson-style normalization:
 *
 *   cm_cor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * where the sum is over edges incident to vertex v. This normalization separates
 * the variation in y and z, making the measure robust to sparse signal patterns
 * where one function varies on only a subset of edges.
 *
 * GEOMETRIC INTERPRETATION:
 * For smooth functions on manifolds with derivative weighting (w_e = 1/ℓ_e²),
 * this measure converges to cos(θ), where θ is the angle between gradient
 * vectors ∇y and ∇z in the tangent space. Values near +1 indicate parallel
 * gradients (functions increase together), -1 indicates anti-parallel gradients
 * (functions vary oppositely), and 0 indicates orthogonal gradients (independent
 * variation).
 *
 * ADVANTAGES:
 * - Naturally down-weights sparse signals (if z varies on 1/30 edges,
 *   coefficient scales by ~1/√30)
 * - Rigorous geometric interpretation via manifold theory
 * - Compatible with correlation-based downstream analyses
 * - Self-regularizing: no need for explicit thresholds
 *
 * NUMERICAL STABILITY:
 * Returns 0 when either function has insufficient variation (denominator < 1e-10).
 * This occurs when one or both functions are essentially constant across all
 * edges incident to a vertex.
 *
 * @param y Response function values at vertices (length = num_vertices)
 * @param z Feature function values at vertices (length = num_vertices)
 * @param type Weighting scheme:
 *   - UNIT: w_e = 1 (combinatorial, ignores edge lengths)
 *   - DERIVATIVE: w_e = 1/ℓ_e² (geometric, compares directional derivatives)
 *
 * @return comono_result_t structure containing:
 *   - vertex_coefficients: Vector of length n with cm_cor(y,z)(v) for each vertex v
 *   - mean_coefficient: Mean of vertex_coefficients
 *   - median_coefficient: Median of vertex_coefficients
 *   - n_positive: Count of vertices with coefficient > 1e-10
 *   - n_negative: Count of vertices with coefficient < -1e-10
 *   - n_zero: Count of vertices with |coefficient| ≤ 1e-10
 *
 *
 * WHEN TO USE:
 * - Default choice for most applications
 * - Continuous functions with real variation across neighborhoods
 * - When gradient correlation is the scientific question
 * - Microbiome data with continuous abundance values after transformation
 *
 * COMPARISON TO OTHER VARIANTS:
 * - vs. comono(): cm_cor has smaller magnitude when signal is sparse; more robust
 * - vs. comono_proportion(): cm_cor is magnitude-weighted; proportion is count-based
 *
 * EXAMPLE USAGE:
 * @code
 * set_wgraph_t graph(adj_list, weight_list);
 * std::vector<double> y = ...; // response values
 * std::vector<double> z = ...; // feature values
 *
 * // Compute with derivative weighting (recommended)
 * auto result = graph.comono_cor(y, z, comono_type_t::DERIVATIVE);
 *
 * // Examine vertex-level coefficients
 * for (size_t v = 0; v < graph.num_vertices(); ++v) {
 *     std::cout << "Vertex " << v << ": "
 *               << result.vertex_coefficients[v] << std::endl;
 * }
 * @endcode
 *
 * @see comono_proportion() for threshold-based variant
 * @see comono() for absolute-value normalization (legacy)
 */

#if 0
comono_result_t set_wgraph_t::comono_cor(
    const std::vector<double>& y,
    const std::vector<double>& z,
    comono_type_t type
) const {
    const size_t n_vertices = num_vertices();

    // Numerical stability threshold for denominators
    // If sqrt(sum of squares) < this, function is essentially constant
    const double MIN_DENOMINATOR = 1e-10;

    // Validate input lengths
    if (y.size() != n_vertices) {
        REPORT_ERROR("Length of y (%zu) does not match number of vertices (%zu)\n",
                     y.size(), n_vertices);
    }
    if (z.size() != n_vertices) {
        REPORT_ERROR("Length of z (%zu) does not match number of vertices (%zu)\n",
                     z.size(), n_vertices);
    }

    // Initialize result structure
    comono_result_t result;
    result.vertex_coefficients.resize(n_vertices, 0.0);

    // Compute co-monotonicity at each vertex
    for (size_t v = 0; v < n_vertices; ++v) {
        double numerator = 0.0;       // Σ w_e Δ_e y · Δ_e z
        double sum_y_squared = 0.0;   // Σ w_e (Δ_e y)²
        double sum_z_squared = 0.0;   // Σ w_e (Δ_e z)²

        // Iterate over neighbors of vertex v
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;
            const double edge_length = edge_info.weight;

            // Compute edge differences
            const double delta_y = y[u] - y[v];
            const double delta_z = z[u] - z[v];

            // Compute weight based on type
            double weight = 1.0;
            switch (type) {
                case comono_type_t::UNIT:
                    // Unit weighting: w_e = 1
                    weight = 1.0;
                    break;

                case comono_type_t::DERIVATIVE:
                    // Derivative weighting: w_e = 1/ℓ_e²
                    // Guard against division by zero for degenerate edges
                    if (edge_length > 1e-10) {
                        weight = 1.0 / (edge_length * edge_length);
                    } else {
                        // Skip degenerate edges (zero or near-zero length)
                        // These would have infinite weight and dominate the calculation
                        continue;
                    }
                    break;

                case comono_type_t::SIGN:
                    // Sign-based type not applicable for correlation normalization
                    // Fall back to unit weighting
                    weight = 1.0;
                    break;
            }

            // Accumulate weighted sums
            // Note: ALL edges contribute to all three sums
            // This is key difference from absolute-value normalization
            numerator += weight * delta_y * delta_z;
            sum_y_squared += weight * delta_y * delta_y;
            sum_z_squared += weight * delta_z * delta_z;
        }

        // Compute denominators (square roots of sum of squares)
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);

        // Check numerical stability
        // If either function is essentially constant in this neighborhood,
        // the coefficient is undefined (or 0 by convention)
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            // Standard case: both functions vary sufficiently
            result.vertex_coefficients[v] = numerator / (denom_y * denom_z);

            // Clamp to [-1, 1] to handle numerical errors
            // (should not be necessary mathematically, but floating-point arithmetic...)
            if (result.vertex_coefficients[v] > 1.0) {
                result.vertex_coefficients[v] = 1.0;
            } else if (result.vertex_coefficients[v] < -1.0) {
                result.vertex_coefficients[v] = -1.0;
            }
        } else {
            // One or both functions essentially constant
            // Return 0 by convention (undefined correlation)
            result.vertex_coefficients[v] = 0.0;
        }
    }

    // Compute summary statistics
    result.mean_coefficient = std::accumulate(
        result.vertex_coefficients.begin(),
        result.vertex_coefficients.end(),
        0.0
    ) / static_cast<double>(n_vertices);

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
#endif


#if DEBUG_COMONO_COR
#include <fstream>
#include <sstream>
#include <map>

// Structure to hold debug pair information for correlation comparison
struct debug_pair_cor_t {
    size_t vertex_idx;
    size_t phylotype_idx;
    std::string phylotype_name;
    double unit_value;
    double deriv_value;
    double difference;
    std::string case_type;
};

// Global storage for debug pairs
static std::vector<debug_pair_cor_t> g_debug_pairs_cor;
static bool g_debug_pairs_cor_loaded = false;

// Function to load debug pairs from CSV
void load_debug_pairs_cor(const std::string& filepath) {
    if (g_debug_pairs_cor_loaded) {
        return;  // Already loaded
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        Rprintf("Warning: Could not open debug file: %s\n", filepath.c_str());
        Rprintf("Make sure to run find.discrepant.pairs.cor() first!\n");
        return;
    }

    std::string line;
    std::getline(file, line);  // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        debug_pair_cor_t pair;

        // Parse CSV: vertex_idx,phylotype_idx,phylotype_name,unit_value,deriv_value,difference,case
        std::getline(ss, item, ',');
        pair.vertex_idx = std::stoull(item);

        std::getline(ss, item, ',');
        pair.phylotype_idx = std::stoull(item);

        std::getline(ss, item, ',');
        // Handle quoted strings
        if (item.front() == '"' && item.back() == '"') {
            item = item.substr(1, item.length() - 2);
        }
        pair.phylotype_name = item;

        std::getline(ss, item, ',');
        pair.unit_value = std::stod(item);

        std::getline(ss, item, ',');
        pair.deriv_value = std::stod(item);

        std::getline(ss, item, ',');
        pair.difference = std::stod(item);

        std::getline(ss, item, ',');
        // Handle quoted strings
        if (item.front() == '"' && item.back() == '"') {
            item = item.substr(1, item.length() - 2);
        }
        pair.case_type = item;

        g_debug_pairs_cor.push_back(pair);
    }

    file.close();
    g_debug_pairs_cor_loaded = true;

    Rprintf("========================================\n");
    Rprintf("Loaded %zu correlation debug pairs from %s\n",
            g_debug_pairs_cor.size(), filepath.c_str());
    Rprintf("========================================\n\n");
}

// Check if a vertex should be debugged
bool should_debug_vertex_cor(size_t vertex_idx) {
    for (const auto& pair : g_debug_pairs_cor) {
        if (pair.vertex_idx == vertex_idx) {
            return true;
        }
    }
    return false;
}

// Get debug pair info for printing context
const debug_pair_cor_t* get_debug_pair_info_cor(size_t vertex_idx) {
    for (const auto& pair : g_debug_pairs_cor) {
        if (pair.vertex_idx == vertex_idx) {
            return &pair;
        }
    }
    return nullptr;
}
#endif  // DEBUG_COMONO_COR


comono_result_t set_wgraph_t::comono_cor(
    const std::vector<double>& y,
    const std::vector<double>& z,
    comono_type_t type
) const {
    const size_t n_vertices = num_vertices();
    const double MIN_DENOMINATOR = 1e-10;

#if DEBUG_COMONO_COR
    // Load debug pairs on first call
    if (!g_debug_pairs_cor_loaded) {
        const char* debug_file_path = "/tmp/comono_debugging_dir/comono_cor_debug_pairs.csv";
        load_debug_pairs_cor(debug_file_path);
    }
#endif

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
        double numerator = 0.0;       // Σ w_e Δ_e y · Δ_e z
        double sum_y_squared = 0.0;   // Σ w_e (Δ_e y)²
        double sum_z_squared = 0.0;   // Σ w_e (Δ_e z)²

#if DEBUG_COMONO_COR
        bool debug_this_vertex = should_debug_vertex_cor(v);
        const debug_pair_cor_t* debug_info = nullptr;

        if (debug_this_vertex) {
            debug_info = get_debug_pair_info_cor(v);
            Rprintf("\n========================================\n");
            Rprintf("DEBUG COMONO_COR: Vertex %zu (phylotype: %s)\n",
                    v, debug_info->phylotype_name.c_str());
            Rprintf("Expected values - Unit: %.6f, Deriv: %.6f, Diff: %.6f\n",
                    debug_info->unit_value, debug_info->deriv_value, debug_info->difference);
            Rprintf("Case: %s\n", debug_info->case_type.c_str());
            Rprintf("========================================\n");
            Rprintf("Type: %s\n", (type == comono_type_t::UNIT) ? "UNIT" : "DERIVATIVE");
            Rprintf("Response value at vertex: y[%zu] = %.6f\n", v, y[v]);
            Rprintf("Feature value at vertex: z[%zu] = %.6f\n", v, z[v]);
            Rprintf("Number of neighbors: %zu\n", adjacency_list[v].size());
            Rprintf("\nFormula: cm_cor = Σ(w·Δy·Δz) / √[Σ(w·(Δy)²)] √[Σ(w·(Δz)²)]\n\n");
        }
#endif

        // Iterate over neighbors of vertex v
        size_t neighbor_idx = 0;
        for (const auto& edge_info : adjacency_list[v]) {
            const size_t u = edge_info.vertex;
            const double edge_length = edge_info.weight;

            // Compute edge differences
            const double delta_y = y[u] - y[v];
            const double delta_z = z[u] - z[v];

#if DEBUG_COMONO_COR
            if (debug_this_vertex) {
                Rprintf("  Neighbor %zu (vertex %zu):\n", neighbor_idx, u);
                Rprintf("    Edge length: %.6f\n", edge_length);
                Rprintf("    y[%zu] = %.6f, delta_y = %.6f\n", u, y[u], delta_y);
                Rprintf("    z[%zu] = %.6f, delta_z = %.6f\n", u, z[u], delta_z);
            }
#endif

            // Compute weight based on type
            double weight = 1.0;
            switch (type) {
                case comono_type_t::UNIT:
                    weight = 1.0;
#if DEBUG_COMONO_COR
                    if (debug_this_vertex) {
                        Rprintf("    weight (UNIT) = 1.0\n");
                    }
#endif
                    break;

                case comono_type_t::DERIVATIVE:
                    if (edge_length > 1e-10) {
                        weight = 1.0 / (edge_length * edge_length);
                    } else {
#if DEBUG_COMONO_COR
                        if (debug_this_vertex) {
                            Rprintf("    SKIPPING: edge length too small (%.10f)\n", edge_length);
                        }
#endif
                        ++neighbor_idx;
                        continue;
                    }
#if DEBUG_COMONO_COR
                    if (debug_this_vertex) {
                        Rprintf("    weight (DERIVATIVE) = 1/(%.6f)² = %.6f\n", edge_length, weight);
                    }
#endif
                    break;

                case comono_type_t::SIGN:
                    // Not applicable for cm_cor
                    weight = 1.0;
                    break;
            }

            // Accumulate weighted sums
            const double weighted_product = weight * delta_y * delta_z;
            const double weighted_y_sq = weight * delta_y * delta_y;
            const double weighted_z_sq = weight * delta_z * delta_z;

            numerator += weighted_product;
            sum_y_squared += weighted_y_sq;
            sum_z_squared += weighted_z_sq;

#if DEBUG_COMONO_COR
            if (debug_this_vertex) {
                Rprintf("    product (Δy · Δz) = %.6f\n", delta_y * delta_z);
                Rprintf("    weighted product = %.6f * %.6f = %.6f\n",
                        weight, delta_y * delta_z, weighted_product);
                Rprintf("    weighted (Δy)² = %.6f * %.6f = %.6f\n",
                        weight, delta_y * delta_y, weighted_y_sq);
                Rprintf("    weighted (Δz)² = %.6f * %.6f = %.6f\n",
                        weight, delta_z * delta_z, weighted_z_sq);
                Rprintf("    Running numerator: %.6f\n", numerator);
                Rprintf("    Running Σ(w·(Δy)²): %.6f\n", sum_y_squared);
                Rprintf("    Running Σ(w·(Δz)²): %.6f\n", sum_z_squared);
                Rprintf("\n");
            }
#endif
            ++neighbor_idx;
        }

        // Compute denominators (square roots of sum of squares)
        const double denom_y = std::sqrt(sum_y_squared);
        const double denom_z = std::sqrt(sum_z_squared);

#if DEBUG_COMONO_COR
        if (debug_this_vertex) {
            Rprintf("----------------------------------------\n");
            Rprintf("FINAL CALCULATION:\n");
            Rprintf("  Numerator: Σ(w·Δy·Δz) = %.6f\n", numerator);
            Rprintf("  Sum of y-squares: Σ(w·(Δy)²) = %.6f\n", sum_y_squared);
            Rprintf("  Sum of z-squares: Σ(w·(Δz)²) = %.6f\n", sum_z_squared);
            Rprintf("  Denominator y: √[Σ(w·(Δy)²)] = %.6f\n", denom_y);
            Rprintf("  Denominator z: √[Σ(w·(Δz)²)] = %.6f\n", denom_z);
            Rprintf("  Product of denominators: %.6f * %.6f = %.6f\n",
                    denom_y, denom_z, denom_y * denom_z);
        }
#endif

        // Check numerical stability
        if (denom_y > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
            result.vertex_coefficients[v] = numerator / (denom_y * denom_z);

            // Clamp to [-1, 1]
            if (result.vertex_coefficients[v] > 1.0) {
                result.vertex_coefficients[v] = 1.0;
            } else if (result.vertex_coefficients[v] < -1.0) {
                result.vertex_coefficients[v] = -1.0;
            }

#if DEBUG_COMONO_COR
            if (debug_this_vertex) {
                Rprintf("  Coefficient: %.6f / %.6f = %.6f\n",
                        numerator, denom_y * denom_z, result.vertex_coefficients[v]);
                if (type == comono_type_t::UNIT) {
                    Rprintf("  Expected (UNIT): %.6f\n", debug_info->unit_value);
                    Rprintf("  Error: %.6f\n",
                            result.vertex_coefficients[v] - debug_info->unit_value);
                } else {
                    Rprintf("  Expected (DERIVATIVE): %.6f\n", debug_info->deriv_value);
                    Rprintf("  Error: %.6f\n",
                            result.vertex_coefficients[v] - debug_info->deriv_value);
                }
                Rprintf("========================================\n\n");
            }
#endif
        } else {
            result.vertex_coefficients[v] = 0.0;

#if DEBUG_COMONO_COR
            if (debug_this_vertex) {
                Rprintf("  Coefficient: 0.0 (insufficient variation)\n");
                if (denom_y <= MIN_DENOMINATOR) {
                    Rprintf("  Reason: y denominator too small (%.10f)\n", denom_y);
                }
                if (denom_z <= MIN_DENOMINATOR) {
                    Rprintf("  Reason: z denominator too small (%.10f)\n", denom_z);
                }
                Rprintf("========================================\n\n");
            }
#endif
        }
    }

    // Compute summary statistics (same as before)
    result.mean_coefficient = std::accumulate(
        result.vertex_coefficients.begin(),
        result.vertex_coefficients.end(),
        0.0
    ) / n_vertices;

    std::vector<double> sorted_coeffs = result.vertex_coefficients;
    std::sort(sorted_coeffs.begin(), sorted_coeffs.end());
    if (n_vertices % 2 == 0) {
        result.median_coefficient = 0.5 * (
            sorted_coeffs[n_vertices / 2 - 1] + sorted_coeffs[n_vertices / 2]
        );
    } else {
        result.median_coefficient = sorted_coeffs[n_vertices / 2];
    }

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
