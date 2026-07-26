#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

namespace {

using adjacency_t = std::vector<std::vector<std::size_t>>;

std::size_t compute_hop_extremp_radius(
    std::size_t vertex,
    const adjacency_t& adjacency,
    const double* vertex_densities,
    const double* y,
    double p_threshold,
    bool detect_maxima,
    std::size_t max_hop
) {
    const std::size_t n = adjacency.size();

    double global_min = y[0];
    double global_max = y[0];
    for (std::size_t i = 1; i < n; ++i) {
        global_min = std::min(global_min, y[i]);
        global_max = std::max(global_max, y[i]);
    }

    const double epsilon =
        1e-14 * std::max(std::abs(global_min), std::abs(global_max));
    const bool is_global_min = std::abs(y[vertex] - global_min) <= epsilon;
    const bool is_global_max = std::abs(y[vertex] - global_max) <= epsilon;

    if ((detect_maxima && is_global_max) ||
        (!detect_maxima && is_global_min)) {
        return std::numeric_limits<std::size_t>::max();
    }

    for (const std::size_t neighbor : adjacency[vertex]) {
        if ((detect_maxima && y[neighbor] >= y[vertex]) ||
            (!detect_maxima && y[neighbor] <= y[vertex])) {
            return 0;
        }
    }

    std::vector<bool> visited(n, false);
    std::queue<std::size_t> current_level;
    std::queue<std::size_t> next_level;
    visited[vertex] = true;
    current_level.push(vertex);

    double total_density = 0.0;
    double satisfying_density = 0.0;
    std::size_t current_hop = 0;
    std::size_t last_valid_hop = 0;

    while (!current_level.empty() && current_hop < max_hop) {
        while (!current_level.empty()) {
            const std::size_t current = current_level.front();
            current_level.pop();

            for (const std::size_t neighbor : adjacency[current]) {
                if (visited[neighbor]) {
                    continue;
                }

                visited[neighbor] = true;
                next_level.push(neighbor);
                total_density += vertex_densities[neighbor];

                const bool satisfies =
                    detect_maxima
                        ? y[vertex] > y[neighbor]
                        : y[vertex] < y[neighbor];
                if (satisfies) {
                    satisfying_density += vertex_densities[neighbor];
                }
            }
        }

        ++current_hop;
        if (total_density < std::numeric_limits<double>::epsilon()) {
            return last_valid_hop;
        }

        if (satisfying_density / total_density >= p_threshold) {
            last_valid_hop = current_hop;
        } else {
            return last_valid_hop;
        }

        std::swap(current_level, next_level);
    }

    return last_valid_hop;
}

}  // namespace

extern "C" SEXP S_compute_hop_extremp_radii_batch(
    SEXP s_adj_list,
    SEXP s_edge_densities,
    SEXP s_vertex_densities,
    SEXP s_candidates,
    SEXP s_y,
    SEXP s_p_threshold,
    SEXP s_detect_maxima,
    SEXP s_max_hop
) {
    if (!Rf_isNewList(s_adj_list)) {
        Rf_error("S_compute_hop_extremp_radii_batch: adj_list must be a list");
    }
    if (!Rf_isReal(s_edge_densities)) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: edge_densities must be numeric"
        );
    }
    if (!Rf_isReal(s_vertex_densities)) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: vertex_densities must be numeric"
        );
    }
    if (!Rf_isInteger(s_candidates)) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: candidates must be an integer vector"
        );
    }
    if (!Rf_isReal(s_y)) {
        Rf_error("S_compute_hop_extremp_radii_batch: y must be a numeric vector");
    }
    if (!Rf_isReal(s_p_threshold) || Rf_length(s_p_threshold) != 1) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: p_threshold must be a numeric scalar"
        );
    }
    if (!Rf_isLogical(s_detect_maxima) ||
        Rf_length(s_detect_maxima) != 1) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: detect_maxima must be a logical scalar"
        );
    }
    if (!Rf_isInteger(s_max_hop) || Rf_length(s_max_hop) != 1) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: max_hop must be an integer scalar"
        );
    }

    const std::size_t n = static_cast<std::size_t>(Rf_length(s_adj_list));
    if (n == 0) {
        Rf_error("S_compute_hop_extremp_radii_batch: graph must not be empty");
    }
    if (static_cast<std::size_t>(Rf_length(s_vertex_densities)) != n ||
        static_cast<std::size_t>(Rf_length(s_y)) != n) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: y and vertex_densities must match graph size"
        );
    }

    const double p_threshold = REAL(s_p_threshold)[0];
    const int max_hop_value = INTEGER(s_max_hop)[0];
    if (p_threshold <= 0.0 || p_threshold > 1.0) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: p_threshold must be in (0,1]"
        );
    }
    if (max_hop_value <= 0) {
        Rf_error(
            "S_compute_hop_extremp_radii_batch: max_hop must be positive"
        );
    }

    adjacency_t adjacency(n);
    for (std::size_t i = 0; i < n; ++i) {
        SEXP s_neighbors = VECTOR_ELT(s_adj_list, i);
        if (!Rf_isInteger(s_neighbors)) {
            Rf_error(
                "S_compute_hop_extremp_radii_batch: adjacency entries must be integer vectors"
            );
        }

        const R_xlen_t degree = XLENGTH(s_neighbors);
        adjacency[i].reserve(static_cast<std::size_t>(degree));
        for (R_xlen_t j = 0; j < degree; ++j) {
            const int neighbor = INTEGER(s_neighbors)[j];
            if (neighbor < 1 || neighbor > static_cast<int>(n)) {
                Rf_error(
                    "S_compute_hop_extremp_radii_batch: neighbor index out of bounds"
                );
            }
            adjacency[i].push_back(static_cast<std::size_t>(neighbor - 1));
        }
    }

    const R_xlen_t n_candidates = XLENGTH(s_candidates);
    SEXP result = PROTECT(Rf_allocVector(INTSXP, n_candidates));
    const bool detect_maxima = LOGICAL(s_detect_maxima)[0] != 0;
    const double* vertex_densities = REAL(s_vertex_densities);
    const double* y = REAL(s_y);

    for (R_xlen_t i = 0; i < n_candidates; ++i) {
        const int candidate = INTEGER(s_candidates)[i];
        if (candidate < 1 || candidate > static_cast<int>(n)) {
            UNPROTECT(1);
            Rf_error(
                "S_compute_hop_extremp_radii_batch: candidate index out of bounds"
            );
        }

        const std::size_t radius = compute_hop_extremp_radius(
            static_cast<std::size_t>(candidate - 1),
            adjacency,
            vertex_densities,
            y,
            p_threshold,
            detect_maxima,
            static_cast<std::size_t>(max_hop_value)
        );

        INTEGER(result)[i] =
            radius == std::numeric_limits<std::size_t>::max()
                ? -1
                : static_cast<int>(radius);
    }

    UNPROTECT(1);
    return result;
}
