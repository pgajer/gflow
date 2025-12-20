#include "set_wgraph.hpp"
#include "gfc.hpp"
#include <queue>
#include <algorithm>
#include <cmath>
#include <R_ext/Print.h>

/**
 * @brief Compute edge length density weights for gradient flow modulation
 *
 * Uses kernel density estimation on edge lengths, normalized so max = 1.
 * Edges near the mode of the length distribution receive weight ≈ 1,
 * while very short or very long edges are down-weighted.
 *
 * @param bandwidth KDE bandwidth; if <= 0, uses Silverman's rule
 * @return Map from (min_vertex, max_vertex) pair to density weight in [0,1]
 */
std::unordered_map<size_t, std::unordered_map<size_t, double>>
set_wgraph_t::compute_edge_length_weights(double bandwidth) const {

    // Collect all edge lengths
    std::vector<double> all_lengths;
	size_t n_vertices = adjacency_list.size();

    for (size_t v = 0; v < n_vertices; ++v) {
        for (const auto& edge : adjacency_list[v]) {
            if (edge.vertex > v) {  // Avoid double-counting
                all_lengths.push_back(edge.weight);
            }
        }
    }

    if (all_lengths.empty()) {
        return {};
    }

    // Compute bandwidth using Silverman's rule if not specified
    if (bandwidth <= 0) {
        double n = static_cast<double>(all_lengths.size());
        double mean = 0.0, var = 0.0;
        for (double len : all_lengths) mean += len;
        mean /= n;
        for (double len : all_lengths) var += (len - mean) * (len - mean);
        var /= (n - 1);
        double sd = std::sqrt(var);

        // Silverman's rule of thumb
        bandwidth = 1.06 * sd * std::pow(n, -0.2);
    }

    // Compute KDE at each edge length
    std::unordered_map<size_t, std::unordered_map<size_t, double>> weights;
    double max_density = 0.0;

    // First pass: compute densities
    std::vector<double> densities;
    densities.reserve(all_lengths.size());

    for (double len : all_lengths) {
        double density = 0.0;
        for (double other_len : all_lengths) {
            double z = (len - other_len) / bandwidth;
            density += std::exp(-0.5 * z * z);  // Gaussian kernel
        }
        density /= (all_lengths.size() * bandwidth * std::sqrt(2.0 * M_PI));
        densities.push_back(density);
        if (density > max_density) max_density = density;
    }

    // Second pass: assign normalized weights to edges
    size_t edge_idx = 0;
    for (size_t v = 0; v < n_vertices; ++v) {
        for (const auto& edge : adjacency_list[v]) {
            if (edge.vertex > v) {
                double normalized = densities[edge_idx] / max_density;
                weights[v][edge.vertex] = normalized;
                weights[edge.vertex][v] = normalized;
                ++edge_idx;
            }
        }
    }

    return weights;
}

/**
 * @brief Compute vertex density from nearest neighbor distances
 *
 * Estimates density as the conditional expectation of d_1^{-1},
 * where d_1(v) is the distance to the nearest neighbor of v.
 *
 * @return Vector of density values, normalized to [0,1]
 */
std::vector<double> set_wgraph_t::compute_vertex_density() const {

	size_t n_vertices = adjacency_list.size();
    std::vector<double> density(n_vertices, 0.0);
    double max_density = 0.0;

    for (size_t v = 0; v < n_vertices; ++v) {
        if (adjacency_list[v].empty()) {
            density[v] = 0.0;
            continue;
        }

        // Find minimum edge weight (nearest neighbor distance)
        double min_dist = std::numeric_limits<double>::max();
        for (const auto& edge : adjacency_list[v]) {
            if (edge.weight < min_dist) {
                min_dist = edge.weight;
            }
        }

        // Density is inverse of nearest neighbor distance
        density[v] = (min_dist > 0) ? (1.0 / min_dist) : 0.0;
        if (density[v] > max_density) {
            max_density = density[v];
        }
    }

    // Normalize to [0,1]
    if (max_density > 0) {
        for (double& d : density) {
            d /= max_density;
        }
    }

    return density;
}

/**
 * @brief Find the modulated gradient direction from a vertex
 *
 * Returns the neighbor that maximizes the modulated gradient:
 * - NONE: max Δŷ
 * - DENSITY: max ρ(u) · Δŷ
 * - EDGELEN: max dl([v,u]) · Δŷ
 * - DENSITY_EDGELEN: max ρ(u) · dl([v,u]) · Δŷ
 *
 * @param v Current vertex
 * @param y Function values
 * @param ascending If true, find steepest ascent; if false, steepest descent
 * @param modulation Modulation type
 * @param density Vertex density values (for DENSITY modulation)
 * @param edge_weights Edge length weights (for EDGELEN modulation)
 * @return Pair of (next_vertex, is_extremum). If is_extremum, v is a local extremum.
 */
std::pair<size_t, bool> set_wgraph_t::find_modulated_gradient_neighbor(
    size_t v,
    const std::vector<double>& y,
    bool ascending,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    const std::unordered_map<size_t, std::unordered_map<size_t, double>>& edge_weights
) const {

    double best_score = 0.0;
    size_t best_neighbor = v;
    bool found_improvement = false;

    for (const auto& edge : adjacency_list[v]) {
        size_t u = edge.vertex;

        // Compute Δŷ
        double delta_y = ascending ? (y[u] - y[v]) : (y[v] - y[u]);

        if (delta_y <= 0) continue;  // Must be strictly improving

        // Compute modulation factor
        double mod_factor = 1.0;

        if (modulation == gflow_modulation_t::DENSITY ||
            modulation == gflow_modulation_t::DENSITY_EDGELEN) {
            if (!density.empty()) {
                mod_factor *= density[u];
            }
        }

        if (modulation == gflow_modulation_t::EDGELEN ||
            modulation == gflow_modulation_t::DENSITY_EDGELEN) {
            auto it_v = edge_weights.find(v);
            if (it_v != edge_weights.end()) {
                auto it_u = it_v->second.find(u);
                if (it_u != it_v->second.end()) {
                    mod_factor *= it_u->second;
                }
            }
        }

        double score = mod_factor * delta_y;

        if (score > best_score) {
            best_score = score;
            best_neighbor = u;
            found_improvement = true;
        }
    }

    return {best_neighbor, !found_improvement};
}

/**
 * @brief Compute basin of attraction with modulated gradient flow
 *
 * @param extremum_vertex The local extremum (attractor)
 * @param y Function values
 * @param is_maximum True for ascending flow to maximum, false for descending to minimum
 * @param modulation Gradient flow modulation type
 * @param density Vertex densities (for DENSITY modulation)
 * @param edge_weights Edge length weights (for EDGELEN modulation)
 * @return bbasin_t structure with vertices and boundary
 */
bbasin_t set_wgraph_t::compute_single_basin(
    size_t extremum_vertex,
    const std::vector<double>& y,
    bool is_maximum,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    const std::unordered_map<size_t, std::unordered_map<size_t, double>>& edge_weights
) const {

    bbasin_t basin;
    basin.extremum_vertex = extremum_vertex;
    basin.value = y[extremum_vertex];
    basin.is_maximum = is_maximum;

    // Track which vertices flow to this extremum
    std::unordered_set<size_t> basin_set;
    std::unordered_set<size_t> boundary_set;

    // For each vertex, trace gradient flow to see if it reaches this extremum
	size_t n_vertices = adjacency_list.size();
    for (size_t v = 0; v < n_vertices; ++v) {
        size_t current = v;
        std::vector<size_t> trajectory;
        trajectory.push_back(current);

        // Follow gradient flow
        while (true) {
            auto [next, is_extremum] = find_modulated_gradient_neighbor(
                current, y, is_maximum, modulation, density, edge_weights
            );

            if (is_extremum || next == current) {
                // Reached a local extremum
                if (current == extremum_vertex) {
                    // This trajectory belongs to our basin
                    for (size_t traj_v : trajectory) {
                        basin_set.insert(traj_v);
                    }
                }
                break;
            }

            // Check for cycles (shouldn't happen with strict inequality)
            if (std::find(trajectory.begin(), trajectory.end(), next) != trajectory.end()) {
                break;
            }

            trajectory.push_back(next);
            current = next;
        }
    }

    // Convert basin set to vector
    basin.vertices.assign(basin_set.begin(), basin_set.end());
    std::sort(basin.vertices.begin(), basin.vertices.end());

    // Compute boundary: neighbors of basin vertices that are not in basin
    for (size_t v : basin_set) {
        for (const auto& edge : adjacency_list[v]) {
            size_t u = edge.vertex;
            if (basin_set.find(u) == basin_set.end()) {
                boundary_set.insert(u);
            }
        }
    }

    basin.boundary.assign(boundary_set.begin(), boundary_set.end());
    std::sort(basin.boundary.begin(), basin.boundary.end());

    return basin;
}

/**
 * @brief Compute all basins of attraction with modulated gradient flow
 *
 * Identifies all local extrema and computes their basins of attraction
 * using the specified gradient flow modulation.
 *
 * @param y Function values at vertices
 * @param params Basin computation parameters (modulation type, density, etc.)
 * @param verbose Print progress information
 * @return Map from extremum vertex to its basin structure
 */
std::unordered_map<size_t, bbasin_t> set_wgraph_t::compute_gfc_basins(
    const std::vector<double>& y,
    const gfc_basin_params_t& params,
    bool verbose
) const {

    std::unordered_map<size_t, bbasin_t> result;

    // Prepare modulation data
    std::vector<double> density;
    std::unordered_map<size_t, std::unordered_map<size_t, double>> edge_weights;

    if (params.modulation == gflow_modulation_t::DENSITY ||
        params.modulation == gflow_modulation_t::DENSITY_EDGELEN) {
        if (params.density.empty()) {
            if (verbose) {
                Rprintf("Computing vertex density from nearest neighbor distances...\n");
            }
            density = compute_vertex_density();
        } else {
            density = params.density;
        }
    }

    if (params.modulation == gflow_modulation_t::EDGELEN ||
        params.modulation == gflow_modulation_t::DENSITY_EDGELEN) {
        if (verbose) {
            Rprintf("Computing edge length distribution weights...\n");
        }
        edge_weights = compute_edge_length_weights(params.edgelen_bandwidth);
    }

    // Find all local extrema using modulated gradient
    if (verbose) {
        Rprintf("Finding local extrema with %s modulation...\n",
                params.modulation == gflow_modulation_t::NONE ? "no" :
                params.modulation == gflow_modulation_t::DENSITY ? "density" :
                params.modulation == gflow_modulation_t::EDGELEN ? "edge-length" :
                "density+edge-length");
    }

    std::vector<size_t> local_minima;
    std::vector<size_t> local_maxima;

	size_t n_vertices = adjacency_list.size();
    for (size_t v = 0; v < n_vertices; ++v) {
        // Check if v is a local minimum (no descending neighbor)
        auto [desc_nbr, is_min] = find_modulated_gradient_neighbor(
            v, y, false, params.modulation, density, edge_weights
        );
        if (is_min) {
            local_minima.push_back(v);
        }

        // Check if v is a local maximum (no ascending neighbor)
        auto [asc_nbr, is_max] = find_modulated_gradient_neighbor(
            v, y, true, params.modulation, density, edge_weights
        );
        if (is_max) {
            local_maxima.push_back(v);
        }
    }

    if (verbose) {
        Rprintf("Found %zu local minima, %zu local maxima\n",
                local_minima.size(), local_maxima.size());
    }

    // Compute basin for each minimum
    if (verbose) {
        Rprintf("Computing basins of attraction for minima...\n");
    }
    for (size_t min_v : local_minima) {
        bbasin_t basin = compute_single_basin(
            min_v, y, false, params.modulation, density, edge_weights
        );
        result[min_v] = basin;
    }

    // Compute basin for each maximum
    if (verbose) {
        Rprintf("Computing basins of attraction for maxima...\n");
    }
    for (size_t max_v : local_maxima) {
        bbasin_t basin = compute_single_basin(
            max_v, y, true, params.modulation, density, edge_weights
        );
        result[max_v] = basin;
    }

    if (verbose) {
        Rprintf("Basin computation complete: %zu total basins\n", result.size());
    }

    return result;
}
