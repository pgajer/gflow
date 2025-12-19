#include "gfassoc_polarity.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>


/**
 * @brief Compute theta for a single cell
 *
 * The normalized height within a cell is:
 *   theta = (y(v) - y_min) / (y_max - y_min + epsilon)
 *
 * This is guaranteed to be in [0, 1] when y_min <= y(v) <= y_max,
 * which holds by the gradient flow property.
 */
double compute_cell_theta(
    double y_v,
    double y_max,
    double y_min,
    double epsilon
) {
    double range = y_max - y_min;

    // Check for flat cell
    if (range < epsilon) {
        return -1.0;  // Invalid marker
    }

    double theta = (y_v - y_min) / (range + epsilon);

    // Clamp to [0, 1] for numerical safety
    if (theta < 0.0) theta = 0.0;
    if (theta > 1.0) theta = 1.0;

    return theta;
}


/**
 * @brief Compute polarity coordinates for a function
 *
 * For each vertex, we compute:
 * 1. For each cell containing the vertex, compute cell-specific theta
 * 2. Average theta values weighted by cell membership
 * 3. Convert to polarity: p = 2*theta - 1
 *
 * Time complexity: O(n * average_cell_multiplicity)
 */
polarity_result_t compute_polarity(
    const std::vector<double>& y,
    const basin_membership_t& membership,
    const cell_membership_t& cells,
    double epsilon
) {
    polarity_result_t result;
    size_t n_vertices = y.size();

    result.theta.resize(n_vertices, 0.0);
    result.polarity.resize(n_vertices, 0.0);
    result.range.resize(n_vertices, 0.0);
    result.is_valid.resize(n_vertices, false);
    result.epsilon = epsilon;

    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& cell_idx = cells.cell_indices[v];
        const auto& cell_wt = cells.cell_membership[v];

        // Skip vertices not in any cell
        if (cell_idx.empty()) {
            result.is_valid[v] = false;
            continue;
        }

        double weighted_theta = 0.0;
        double weighted_range = 0.0;
        double total_valid_weight = 0.0;

        for (size_t k = 0; k < cell_idx.size(); ++k) {
            size_t max_idx = cell_idx[k].first;
            size_t min_idx = cell_idx[k].second;

            double y_max = membership.max_values[max_idx];
            double y_min = membership.min_values[min_idx];
            double cell_range = y_max - y_min;

            double theta_k = compute_cell_theta(y[v], y_max, y_min, epsilon);

            // Skip invalid (flat) cells
            if (theta_k < 0.0) {
                continue;
            }

            weighted_theta += cell_wt[k] * theta_k;
            weighted_range += cell_wt[k] * cell_range;
            total_valid_weight += cell_wt[k];
        }

        // Normalize if we have valid contributions
        if (total_valid_weight > epsilon) {
            result.theta[v] = weighted_theta / total_valid_weight;
            result.range[v] = weighted_range / total_valid_weight;
            result.polarity[v] = 2.0 * result.theta[v] - 1.0;
            result.is_valid[v] = true;
        } else {
            result.is_valid[v] = false;
        }
    }

    return result;
}


/**
 * @brief Compute polarity using rank transformation
 *
 * For each cell, we rank the vertices by their y values and use the
 * normalized rank as theta:
 *   theta_rank(v) = (rank(v) - 1) / (n_cell - 1)
 *
 * This provides invariance under monotone transformations.
 */
polarity_result_t compute_polarity_rank(
    const std::vector<double>& y,
    const basin_membership_t& membership,
    const cell_membership_t& cells,
    double epsilon
) {
    polarity_result_t result;
    size_t n_vertices = y.size();

    result.theta.resize(n_vertices, 0.0);
    result.polarity.resize(n_vertices, 0.0);
    result.range.resize(n_vertices, 0.0);
    result.is_valid.resize(n_vertices, false);
    result.epsilon = epsilon;

    // Build a map from cell (i,j) -> list of (vertex, y_value)
    std::unordered_map<size_t, std::vector<std::pair<size_t, double>>> cell_vertex_map;

    // Use a simple hash for (i, j) pairs
    auto cell_hash = [&](size_t i, size_t j) -> size_t {
        return i * membership.n_min_basins + j;
    };

    // First pass: collect vertices for each cell
    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& cell_idx = cells.cell_indices[v];
        for (const auto& [max_idx, min_idx] : cell_idx) {
            size_t hash = cell_hash(max_idx, min_idx);
            cell_vertex_map[hash].emplace_back(v, y[v]);
        }
    }

    // Second pass: compute ranks within each cell
    // Store rank-based theta for each (vertex, cell) pair
    std::vector<std::unordered_map<size_t, double>> vertex_cell_theta(n_vertices);

    for (auto& [hash, vertex_list] : cell_vertex_map) {
        size_t n_cell = vertex_list.size();

        // Skip cells that are too small
        if (n_cell < 2) {
            continue;
        }

        // Sort by y value
        std::sort(vertex_list.begin(), vertex_list.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });

        // Assign ranks (handling ties by average rank)
        std::vector<double> ranks(n_cell);
        size_t i = 0;
        while (i < n_cell) {
            size_t j = i;
            // Find group of ties
            while (j < n_cell && std::abs(vertex_list[j].second - vertex_list[i].second) < epsilon) {
                ++j;
            }
            // Average rank for ties
            double avg_rank = 0.5 * (i + j - 1);  // 0-based average
            for (size_t k = i; k < j; ++k) {
                ranks[k] = avg_rank;
            }
            i = j;
        }

        // Normalize ranks to [0, 1]
        double max_rank = static_cast<double>(n_cell - 1);
        for (size_t k = 0; k < n_cell; ++k) {
            size_t v = vertex_list[k].first;
            double theta_rank = ranks[k] / max_rank;
            vertex_cell_theta[v][hash] = theta_rank;
        }
    }

    // Third pass: average theta values weighted by cell membership
    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& cell_idx = cells.cell_indices[v];
        const auto& cell_wt = cells.cell_membership[v];

        if (cell_idx.empty()) {
            result.is_valid[v] = false;
            continue;
        }

        double weighted_theta = 0.0;
        double weighted_range = 0.0;
        double total_valid_weight = 0.0;

        for (size_t k = 0; k < cell_idx.size(); ++k) {
            size_t max_idx = cell_idx[k].first;
            size_t min_idx = cell_idx[k].second;
            size_t hash = cell_hash(max_idx, min_idx);

            // Check if we have a rank for this cell
            auto it = vertex_cell_theta[v].find(hash);
            if (it == vertex_cell_theta[v].end()) {
                continue;  // Cell was too small
            }

            double theta_k = it->second;
            double cell_range = membership.max_values[max_idx] - membership.min_values[min_idx];

            weighted_theta += cell_wt[k] * theta_k;
            weighted_range += cell_wt[k] * cell_range;
            total_valid_weight += cell_wt[k];
        }

        if (total_valid_weight > epsilon) {
            result.theta[v] = weighted_theta / total_valid_weight;
            result.range[v] = weighted_range / total_valid_weight;
            result.polarity[v] = 2.0 * result.theta[v] - 1.0;
            result.is_valid[v] = true;
        } else {
            result.is_valid[v] = false;
        }
    }

    return result;
}


/**
 * @brief Compute vertex-level association from polarities
 *
 * The association score is simply the product of polarities:
 *   a_pol(v) = p_y(v) * p_z(v)
 *
 * This is positive when both polarities have the same sign (both high or both low),
 * and negative when they have opposite signs.
 */
vertex_association_t compute_vertex_association(
    const polarity_result_t& pol_y,
    const polarity_result_t& pol_z
) {
    vertex_association_t result;
    size_t n_vertices = pol_y.polarity.size();

    if (pol_z.polarity.size() != n_vertices) {
        throw std::invalid_argument(
            "Polarity vectors must have the same length"
        );
    }

    result.a_pol.resize(n_vertices, 0.0);
    result.sign_pol.resize(n_vertices, 0.0);
    result.confidence.resize(n_vertices, 0.0);
    result.is_valid.resize(n_vertices, false);

    for (size_t v = 0; v < n_vertices; ++v) {
        // Both must be valid for association to be defined
        if (!pol_y.is_valid[v] || !pol_z.is_valid[v]) {
            result.is_valid[v] = false;
            continue;
        }

        double p_y = pol_y.polarity[v];
        double p_z = pol_z.polarity[v];

        result.a_pol[v] = p_y * p_z;
        result.confidence[v] = std::abs(result.a_pol[v]);

        // Compute sign with tolerance
        if (result.a_pol[v] > 1e-10) {
            result.sign_pol[v] = 1.0;
        } else if (result.a_pol[v] < -1e-10) {
            result.sign_pol[v] = -1.0;
        } else {
            result.sign_pol[v] = 0.0;
        }

        result.is_valid[v] = true;
    }

    return result;
}


/**
 * @brief Compute global association statistics
 *
 * Aggregates vertex-level association:
 *   A_pol = sum_v m_0(v) * a_pol(v) / sum_v m_0(v)
 *   kappa_pol = sum_v m_0(v) * sign(a_pol(v)) / sum_v m_0(v)
 *
 * Only valid vertices contribute to the sums.
 */
global_association_t compute_global_association(
    const vertex_association_t& va,
    const std::vector<double>& vertex_mass
) {
    global_association_t result;
    size_t n_vertices = va.a_pol.size();

    result.A_pol = 0.0;
    result.kappa_pol = 0.0;
    result.n_positive = 0;
    result.n_negative = 0;
    result.n_zero = 0;
    result.n_invalid = 0;
    result.total_mass = 0.0;

    double sum_weighted_a_pol = 0.0;
    double sum_weighted_sign = 0.0;

    for (size_t v = 0; v < n_vertices; ++v) {
        double mass = vertex_mass.empty() ? 1.0 : vertex_mass[v];

        if (!va.is_valid[v]) {
            result.n_invalid++;
            continue;
        }

        result.total_mass += mass;
        sum_weighted_a_pol += mass * va.a_pol[v];
        sum_weighted_sign += mass * va.sign_pol[v];

        // Count signs
        if (va.sign_pol[v] > 0) {
            result.n_positive++;
        } else if (va.sign_pol[v] < 0) {
            result.n_negative++;
        } else {
            result.n_zero++;
        }
    }

    // Normalize
    if (result.total_mass > 1e-15) {
        result.A_pol = sum_weighted_a_pol / result.total_mass;
        result.kappa_pol = sum_weighted_sign / result.total_mass;
    }

    return result;
}


/**
 * @brief Compute basin association character
 *
 * For each basin, compute the average polarity of the other function
 * within that basin.
 *
 * For y-max basin i:
 *   chi^{y,+}_i = sum_v m_0(v) * mu^{y,+}_i(v) * p_z(v) / sum_v m_0(v) * mu^{y,+}_i(v)
 */
basin_character_t compute_basin_character(
    const basin_membership_t& y_membership,
    const basin_membership_t& z_membership,
    const polarity_result_t& pol_y,
    const polarity_result_t& pol_z,
    const std::vector<double>& vertex_mass
) {
    basin_character_t result;
    size_t n_vertices = y_membership.n_vertices;

    // Initialize
    result.chi_y_max.resize(y_membership.n_max_basins, 0.0);
    result.chi_y_min.resize(y_membership.n_min_basins, 0.0);
    result.chi_z_max.resize(z_membership.n_max_basins, 0.0);
    result.chi_z_min.resize(z_membership.n_min_basins, 0.0);

    result.mass_y_max.resize(y_membership.n_max_basins, 0.0);
    result.mass_y_min.resize(y_membership.n_min_basins, 0.0);
    result.mass_z_max.resize(z_membership.n_max_basins, 0.0);
    result.mass_z_min.resize(z_membership.n_min_basins, 0.0);

    // Accumulators for weighted sums
    std::vector<double> sum_y_max(y_membership.n_max_basins, 0.0);
    std::vector<double> sum_y_min(y_membership.n_min_basins, 0.0);
    std::vector<double> sum_z_max(z_membership.n_max_basins, 0.0);
    std::vector<double> sum_z_min(z_membership.n_min_basins, 0.0);

    // Iterate over vertices
    for (size_t v = 0; v < n_vertices; ++v) {
        double mass = vertex_mass.empty() ? 1.0 : vertex_mass[v];

        // y-max basins: accumulate z-polarity
        if (pol_z.is_valid[v]) {
            const auto& y_max_idx = y_membership.max_basin_indices[v];
            const auto& y_max_wt = y_membership.max_membership[v];

            for (size_t k = 0; k < y_max_idx.size(); ++k) {
                size_t i = y_max_idx[k];
                double weighted_mass = mass * y_max_wt[k];
                sum_y_max[i] += weighted_mass * pol_z.polarity[v];
                result.mass_y_max[i] += weighted_mass;
            }
        }

        // y-min basins: accumulate z-polarity
        if (pol_z.is_valid[v]) {
            const auto& y_min_idx = y_membership.min_basin_indices[v];
            const auto& y_min_wt = y_membership.min_membership[v];

            for (size_t k = 0; k < y_min_idx.size(); ++k) {
                size_t j = y_min_idx[k];
                double weighted_mass = mass * y_min_wt[k];
                sum_y_min[j] += weighted_mass * pol_z.polarity[v];
                result.mass_y_min[j] += weighted_mass;
            }
        }

        // z-max basins: accumulate y-polarity
        if (pol_y.is_valid[v]) {
            const auto& z_max_idx = z_membership.max_basin_indices[v];
            const auto& z_max_wt = z_membership.max_membership[v];

            for (size_t k = 0; k < z_max_idx.size(); ++k) {
                size_t i = z_max_idx[k];
                double weighted_mass = mass * z_max_wt[k];
                sum_z_max[i] += weighted_mass * pol_y.polarity[v];
                result.mass_z_max[i] += weighted_mass;
            }
        }

        // z-min basins: accumulate y-polarity
        if (pol_y.is_valid[v]) {
            const auto& z_min_idx = z_membership.min_basin_indices[v];
            const auto& z_min_wt = z_membership.min_membership[v];

            for (size_t k = 0; k < z_min_idx.size(); ++k) {
                size_t j = z_min_idx[k];
                double weighted_mass = mass * z_min_wt[k];
                sum_z_min[j] += weighted_mass * pol_y.polarity[v];
                result.mass_z_min[j] += weighted_mass;
            }
        }
    }

    // Normalize to get character values
    for (size_t i = 0; i < y_membership.n_max_basins; ++i) {
        if (result.mass_y_max[i] > 1e-15) {
            result.chi_y_max[i] = sum_y_max[i] / result.mass_y_max[i];
        }
    }

    for (size_t j = 0; j < y_membership.n_min_basins; ++j) {
        if (result.mass_y_min[j] > 1e-15) {
            result.chi_y_min[j] = sum_y_min[j] / result.mass_y_min[j];
        }
    }

    for (size_t i = 0; i < z_membership.n_max_basins; ++i) {
        if (result.mass_z_max[i] > 1e-15) {
            result.chi_z_max[i] = sum_z_max[i] / result.mass_z_max[i];
        }
    }

    for (size_t j = 0; j < z_membership.n_min_basins; ++j) {
        if (result.mass_z_min[j] > 1e-15) {
            result.chi_z_min[j] = sum_z_min[j] / result.mass_z_min[j];
        }
    }

    return result;
}


/**
 * @brief Get vertices belonging to a specific cell
 */
std::vector<size_t> get_cell_vertices(
    const cell_membership_t& cells,
    size_t max_idx,
    size_t min_idx
) {
    std::vector<size_t> vertices;

    for (size_t v = 0; v < cells.n_vertices; ++v) {
        for (const auto& [i, j] : cells.cell_indices[v]) {
            if (i == max_idx && j == min_idx) {
                vertices.push_back(v);
                break;
            }
        }
    }

    return vertices;
}


/**
 * @brief Compute weighted average range for a vertex
 */
double compute_weighted_range(
    size_t v,
    const basin_membership_t& membership,
    const cell_membership_t& cells
) {
    const auto& cell_idx = cells.cell_indices[v];
    const auto& cell_wt = cells.cell_membership[v];

    if (cell_idx.empty()) {
        return 0.0;
    }

    double weighted_range = 0.0;

    for (size_t k = 0; k < cell_idx.size(); ++k) {
        size_t max_idx = cell_idx[k].first;
        size_t min_idx = cell_idx[k].second;

        double cell_range = membership.max_values[max_idx] - membership.min_values[min_idx];
        weighted_range += cell_wt[k] * cell_range;
    }

    return weighted_range;
}
