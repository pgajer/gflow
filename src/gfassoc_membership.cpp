#include "gfassoc_membership.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

/**
 * @brief Compute soft basin membership from gradient basin structures
 *
 * This function iterates through all gradient basins and builds membership
 * vectors by recording which basins contain each vertex. The membership
 * weights are then normalized so that weights sum to 1 for each vertex.
 *
 * The algorithm proceeds in two passes:
 * 1. First pass: collect all basin indices containing each vertex
 * 2. Second pass: normalize weights to sum to 1
 *
 * Time complexity: O(sum of basin sizes)
 */
basin_membership_t compute_basin_membership(
    const std::vector<gradient_basin_t>& max_basins,
    const std::vector<gradient_basin_t>& min_basins,
    size_t n_vertices
) {
    basin_membership_t result;
    result.n_vertices = n_vertices;
    result.n_max_basins = max_basins.size();
    result.n_min_basins = min_basins.size();

    // Initialize storage
    result.max_basin_indices.resize(n_vertices);
    result.min_basin_indices.resize(n_vertices);
    result.max_membership.resize(n_vertices);
    result.min_membership.resize(n_vertices);

    result.max_vertices.resize(max_basins.size());
    result.min_vertices.resize(min_basins.size());
    result.max_values.resize(max_basins.size());
    result.min_values.resize(min_basins.size());

    // Process max basins (descending flow from maxima)
    for (size_t i = 0; i < max_basins.size(); ++i) {
        const gradient_basin_t& basin = max_basins[i];

        // Store extremum information
        result.max_vertices[i] = basin.vertex;
        result.max_values[i] = basin.value;

        // Add this basin index to all vertices in the basin
        for (const auto& [vertex, hop_dist] : basin.hop_dist_map) {
            result.max_basin_indices[vertex].push_back(i);
        }
    }

    // Process min basins (ascending flow from minima)
    for (size_t j = 0; j < min_basins.size(); ++j) {
        const gradient_basin_t& basin = min_basins[j];

        // Store extremum information
        result.min_vertices[j] = basin.vertex;
        result.min_values[j] = basin.value;

        // Add this basin index to all vertices in the basin
        for (const auto& [vertex, hop_dist] : basin.hop_dist_map) {
            result.min_basin_indices[vertex].push_back(j);
        }
    }

    // Normalize membership weights (uniform weighting: each basin gets 1/k)
    for (size_t v = 0; v < n_vertices; ++v) {
        // Max basins
        size_t k_max = result.max_basin_indices[v].size();
        if (k_max > 0) {
            double weight = 1.0 / static_cast<double>(k_max);
            result.max_membership[v].assign(k_max, weight);
        }

        // Min basins
        size_t k_min = result.min_basin_indices[v].size();
        if (k_min > 0) {
            double weight = 1.0 / static_cast<double>(k_min);
            result.min_membership[v].assign(k_min, weight);
        }
    }

    return result;
}


/**
 * @brief Compute cell membership from basin membership
 *
 * A vertex v belongs to cell (i, j) if and only if it belongs to both max
 * basin i and min basin j. The cell membership is the product of basin
 * memberships, normalized to sum to 1.
 *
 * For a vertex belonging to max basins {i1, ..., ip} and min basins {j1, ..., jq},
 * it belongs to p*q cells with unnormalized weights:
 *   c_{ik,jl}(v) = 1  for all combinations
 *
 * After normalization:
 *   gamma_{ik,jl}(v) = 1/(p*q)
 */
cell_membership_t compute_cell_membership(
    const basin_membership_t& membership
) {
    cell_membership_t result;
    result.n_vertices = membership.n_vertices;
    result.n_max_basins = membership.n_max_basins;
    result.n_min_basins = membership.n_min_basins;

    result.cell_indices.resize(membership.n_vertices);
    result.cell_membership.resize(membership.n_vertices);

    for (size_t v = 0; v < membership.n_vertices; ++v) {
        const auto& max_indices = membership.max_basin_indices[v];
        const auto& min_indices = membership.min_basin_indices[v];

        // Skip vertices not in any cell
        if (max_indices.empty() || min_indices.empty()) {
            continue;
        }

        // Generate all (max, min) pairs
        size_t n_cells = max_indices.size() * min_indices.size();
        result.cell_indices[v].reserve(n_cells);
        result.cell_membership[v].reserve(n_cells);

        double weight = 1.0 / static_cast<double>(n_cells);

        for (size_t i : max_indices) {
            for (size_t j : min_indices) {
                result.cell_indices[v].emplace_back(i, j);
                result.cell_membership[v].push_back(weight);
            }
        }
    }

    return result;
}


/**
 * @brief Compute soft overlap matrices between two basin structures
 *
 * The overlap between y-basin i and z-basin j is computed as:
 *   O_{ij} = sum_v m_0(v) * mu^y_i(v) * mu^z_j(v)
 *
 * This iterates over all vertices, accumulating contributions based on
 * the membership weights in both structures.
 *
 * Time complexity: O(n * avg_multiplicity_y * avg_multiplicity_z)
 */
overlap_matrices_t compute_soft_overlap(
    const basin_membership_t& y_membership,
    const basin_membership_t& z_membership,
    const std::vector<double>& vertex_mass
) {
    overlap_matrices_t result;

    size_t n_vertices = y_membership.n_vertices;

    // Validate dimensions
    if (z_membership.n_vertices != n_vertices) {
        throw std::invalid_argument(
            "y_membership and z_membership must have same number of vertices"
        );
    }

    // Initialize overlap matrices
    result.O_pp = Eigen::MatrixXd::Zero(y_membership.n_max_basins, z_membership.n_max_basins);
    result.O_mm = Eigen::MatrixXd::Zero(y_membership.n_min_basins, z_membership.n_min_basins);
    result.O_pm = Eigen::MatrixXd::Zero(y_membership.n_max_basins, z_membership.n_min_basins);
    result.O_mp = Eigen::MatrixXd::Zero(y_membership.n_min_basins, z_membership.n_max_basins);

    result.total_mass = 0.0;

    // Iterate over all vertices
    for (size_t v = 0; v < n_vertices; ++v) {
        double mass = vertex_mass.empty() ? 1.0 : vertex_mass[v];
        result.total_mass += mass;

        const auto& y_max_idx = y_membership.max_basin_indices[v];
        const auto& y_min_idx = y_membership.min_basin_indices[v];
        const auto& y_max_wt = y_membership.max_membership[v];
        const auto& y_min_wt = y_membership.min_membership[v];

        const auto& z_max_idx = z_membership.max_basin_indices[v];
        const auto& z_min_idx = z_membership.min_basin_indices[v];
        const auto& z_max_wt = z_membership.max_membership[v];
        const auto& z_min_wt = z_membership.min_membership[v];

        // O_pp: y-max with z-max
        for (size_t a = 0; a < y_max_idx.size(); ++a) {
            for (size_t b = 0; b < z_max_idx.size(); ++b) {
                result.O_pp(y_max_idx[a], z_max_idx[b]) +=
                    mass * y_max_wt[a] * z_max_wt[b];
            }
        }

        // O_mm: y-min with z-min
        for (size_t a = 0; a < y_min_idx.size(); ++a) {
            for (size_t b = 0; b < z_min_idx.size(); ++b) {
                result.O_mm(y_min_idx[a], z_min_idx[b]) +=
                    mass * y_min_wt[a] * z_min_wt[b];
            }
        }

        // O_pm: y-max with z-min
        for (size_t a = 0; a < y_max_idx.size(); ++a) {
            for (size_t b = 0; b < z_min_idx.size(); ++b) {
                result.O_pm(y_max_idx[a], z_min_idx[b]) +=
                    mass * y_max_wt[a] * z_min_wt[b];
            }
        }

        // O_mp: y-min with z-max
        for (size_t a = 0; a < y_min_idx.size(); ++a) {
            for (size_t b = 0; b < z_max_idx.size(); ++b) {
                result.O_mp(y_min_idx[a], z_max_idx[b]) +=
                    mass * y_min_wt[a] * z_max_wt[b];
            }
        }
    }

    return result;
}


/**
 * @brief Compute deviation from independence for an overlap matrix
 *
 * Under independence, expected overlap is:
 *   E_ij = (sum_k O_ik) * (sum_l O_lj) / (sum_{k,l} O_kl)
 *
 * The standardized deviation follows Pearson residual form:
 *   zeta_ij = (O_ij - E_ij) / sqrt(E_ij * (1 - r_i) * (1 - c_j))
 *
 * where r_i = (row sum i) / total and c_j = (col sum j) / total.
 */
basin_deviation_t compute_basin_deviation(
    const Eigen::MatrixXd& O
) {
    basin_deviation_t result;

    size_t n_rows = O.rows();
    size_t n_cols = O.cols();

    // Handle edge cases
    if (n_rows == 0 || n_cols == 0) {
        result.delta = Eigen::MatrixXd::Zero(n_rows, n_cols);
        result.zeta = Eigen::MatrixXd::Zero(n_rows, n_cols);
        result.expected = Eigen::MatrixXd::Zero(n_rows, n_cols);
        return result;
    }

    // Compute marginals
    Eigen::VectorXd row_sums = O.rowwise().sum();
    Eigen::VectorXd col_sums = O.colwise().sum();
    double total = O.sum();

    // Handle zero total
    if (total < 1e-15) {
        result.delta = Eigen::MatrixXd::Zero(n_rows, n_cols);
        result.zeta = Eigen::MatrixXd::Zero(n_rows, n_cols);
        result.expected = Eigen::MatrixXd::Zero(n_rows, n_cols);
        return result;
    }

    // Compute expected values under independence
    result.expected = (row_sums * col_sums.transpose()) / total;

    // Compute raw deviation
    result.delta = O - result.expected;

    // Compute standardized deviation (Pearson residuals)
    result.zeta = Eigen::MatrixXd::Zero(n_rows, n_cols);

    Eigen::VectorXd row_props = row_sums / total;
    Eigen::VectorXd col_props = col_sums / total;

    for (size_t i = 0; i < n_rows; ++i) {
        for (size_t j = 0; j < n_cols; ++j) {
            double E_ij = result.expected(i, j);
            double r_i = row_props(i);
            double c_j = col_props(j);

            // Avoid division by zero
            double denom = E_ij * (1.0 - r_i) * (1.0 - c_j);
            if (denom > 1e-15) {
                result.zeta(i, j) = result.delta(i, j) / std::sqrt(denom);
            } else {
                result.zeta(i, j) = 0.0;
            }
        }
    }

    return result;
}


/**
 * @brief Get core vertices of a basin by trimming boundary
 *
 * A vertex is on the k-hop boundary if its minimum hop distance to any
 * vertex outside the basin is <= k. We identify core vertices as those
 * with hop distance > boundary_hops from all boundary vertices.
 *
 * Implementation uses the hop_dist_map from the basin structure:
 * vertices with hop_dist <= basin.hop_idx - boundary_hops are considered core.
 */
std::vector<size_t> get_basin_core_vertices(
    const gradient_basin_t& basin,
    size_t boundary_hops
) {
    std::vector<size_t> core_vertices;

    // If basin is too small to have a core, return empty
    if (basin.hop_idx <= boundary_hops) {
        return core_vertices;
    }

    size_t max_core_hop = basin.hop_idx - boundary_hops;

    for (const auto& [vertex, hop_dist] : basin.hop_dist_map) {
        if (hop_dist <= max_core_hop) {
            core_vertices.push_back(vertex);
        }
    }

    return core_vertices;
}


/**
 * @brief Identify boundary vertices of a cell
 *
 * A vertex v in cell C is on the boundary if there exists a neighbor u
 * such that u is not in cell C (either u is in a different cell or u
 * has zero membership in C).
 */
std::vector<size_t> get_cell_boundary_vertices(
    const std::pair<size_t, size_t>& cell_id,
    const cell_membership_t& cells,
    const std::vector<std::vector<size_t>>& adjacency_list
) {
    std::vector<size_t> boundary_vertices;

    // First, collect all vertices that belong to this cell
    std::unordered_set<size_t> cell_vertices;
    for (size_t v = 0; v < cells.n_vertices; ++v) {
        for (const auto& cell_idx : cells.cell_indices[v]) {
            if (cell_idx.first == cell_id.first &&
                cell_idx.second == cell_id.second) {
                cell_vertices.insert(v);
                break;
            }
        }
    }

    // Identify boundary vertices
    for (size_t v : cell_vertices) {
        bool is_boundary = false;
        for (size_t u : adjacency_list[v]) {
            if (cell_vertices.find(u) == cell_vertices.end()) {
                is_boundary = true;
                break;
            }
        }
        if (is_boundary) {
            boundary_vertices.push_back(v);
        }
    }

    return boundary_vertices;
}
