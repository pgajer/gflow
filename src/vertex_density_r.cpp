/**
 * @file vertex_density_r.cpp
 * @brief SEXP interface for vertex density computation
 */

#include <R.h>
#include <Rinternals.h>
#include "set_wgraph.hpp"

#include <vector>

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

/**
 * @brief Compute raw vertex density from nearest neighbor distances
 *
 * Returns d₁⁻¹ for each vertex, where d₁(v) is the minimum edge weight
 * (nearest neighbor distance) from vertex v. Values are NOT normalized.
 *
 * @param s_adj_list Adjacency list (0-based from R conversion)
 * @param s_weight_list Edge weight list (edge lengths)
 * @param s_normalize If TRUE, normalize to [0,1]; if FALSE, return raw d₁⁻¹
 *
 * @return Numeric vector of density values
 */
extern "C" SEXP S_compute_vertex_density(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_normalize
) {
    // Convert graph structure
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    const size_t n = adj_list.size();
    bool normalize = Rf_asLogical(s_normalize);

    // Compute raw d₁⁻¹
    std::vector<double> density(n, 0.0);
    double max_density = 0.0;

    for (size_t v = 0; v < n; ++v) {
        if (graph.adjacency_list[v].empty()) {
            density[v] = 0.0;
            continue;
        }

        // Find minimum edge weight (nearest neighbor distance)
        double min_dist = std::numeric_limits<double>::max();
        for (const auto& edge : graph.adjacency_list[v]) {
            double weight = graph.get_edge_weight(v, edge.vertex);
            if (weight < min_dist) {
                min_dist = weight;
            }
        }

        // Density is inverse of nearest neighbor distance
        density[v] = (min_dist > 0) ? (1.0 / min_dist) : 0.0;
        if (density[v] > max_density) {
            max_density = density[v];
        }
    }

    // Normalize to [0,1] if requested
    if (normalize && max_density > 0) {
        for (double& d : density) {
            d /= max_density;
        }
    }

    // Convert to R vector
    SEXP s_result = PROTECT(Rf_allocVector(REALSXP, static_cast<int>(n)));
    double* p_result = REAL(s_result);
    for (size_t i = 0; i < n; ++i) {
        p_result[i] = density[i];
    }

    UNPROTECT(1);
    return s_result;
}
