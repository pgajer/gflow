#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "error_utils.h"

#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>

#include <R.h>
#include <Rinternals.h>

namespace {

constexpr double kTubeLensTol = 1e-12;

std::pair<std::vector<size_t>, double> compute_shortest_path_trace(
    const set_wgraph_t& graph,
    size_t source,
    size_t target
) {
    const size_t n = graph.adjacency_list.size();
    if (source >= n || target >= n) {
        REPORT_ERROR("Source or target vertex out of range");
    }

    if (source == target) {
        return {{source}, 0.0};
    }

    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);

    using queue_item_t = std::pair<double, size_t>;
    std::priority_queue<queue_item_t,
                        std::vector<queue_item_t>,
                        std::greater<queue_item_t>> pq;

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        const auto [curr_dist, u] = pq.top();
        pq.pop();

        if (curr_dist > dist[u] + kTubeLensTol) {
            continue;
        }
        if (u == target) {
            break;
        }

        for (const auto& edge : graph.adjacency_list[u]) {
            const size_t v = edge.vertex;
            const double next_dist = curr_dist + edge.weight;
            if (next_dist + kTubeLensTol < dist[v]) {
                dist[v] = next_dist;
                parent[v] = static_cast<int>(u);
                pq.push({next_dist, v});
            }
        }
    }

    if (!std::isfinite(dist[target])) {
        return {{}, std::numeric_limits<double>::infinity()};
    }

    std::vector<size_t> path_vertices;
    for (int curr = static_cast<int>(target); curr != -1; curr = parent[static_cast<size_t>(curr)]) {
        path_vertices.push_back(static_cast<size_t>(curr));
    }
    std::reverse(path_vertices.begin(), path_vertices.end());
    return {path_vertices, dist[target]};
}

std::vector<double> compute_shortest_path_distances_within_radius(
    const set_wgraph_t& graph,
    size_t source,
    double radius
) {
    const size_t n = graph.adjacency_list.size();
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());

    using queue_item_t = std::pair<double, size_t>;
    std::priority_queue<queue_item_t,
                        std::vector<queue_item_t>,
                        std::greater<queue_item_t>> pq;

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        const auto [curr_dist, u] = pq.top();
        pq.pop();

        if (curr_dist > dist[u] + kTubeLensTol) {
            continue;
        }
        if (curr_dist > radius + kTubeLensTol) {
            break;
        }

        for (const auto& edge : graph.adjacency_list[u]) {
            const size_t v = edge.vertex;
            const double next_dist = curr_dist + edge.weight;
            if (next_dist > radius + kTubeLensTol) {
                continue;
            }
            if (next_dist + kTubeLensTol < dist[v]) {
                dist[v] = next_dist;
                pq.push({next_dist, v});
            }
        }
    }

    return dist;
}

SEXP size_t_vector_to_r_int(const std::vector<size_t>& values) {
    SEXP out = Rf_allocVector(INTSXP, static_cast<R_xlen_t>(values.size()));
    int* out_ptr = INTEGER(out);
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(values.size()); ++i) {
        out_ptr[i] = static_cast<int>(values[static_cast<size_t>(i)] + 1);
    }
    return out;
}

SEXP double_vector_to_r_num(const std::vector<double>& values) {
    SEXP out = Rf_allocVector(REALSXP, static_cast<R_xlen_t>(values.size()));
    double* out_ptr = REAL(out);
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(values.size()); ++i) {
        out_ptr[i] = values[static_cast<size_t>(i)];
    }
    return out;
}

}  // namespace

tube_lens_corridor_result_t set_wgraph_t::compute_tube_lens_corridor(
    size_t source,
    size_t target,
    double path_relative_radius,
    double excess_tolerance,
    bool compute_excess
) const {
    if (source >= adjacency_list.size() || target >= adjacency_list.size()) {
        REPORT_ERROR("Source or target vertex out of range");
    }
    if (!std::isfinite(path_relative_radius) || path_relative_radius < 0.0) {
        REPORT_ERROR("path_relative_radius must be a finite non-negative number");
    }
    if (compute_excess && std::isfinite(excess_tolerance) && excess_tolerance < 0.0) {
        REPORT_ERROR("excess_tolerance must be non-negative when supplied");
    }

    tube_lens_corridor_result_t result;

    const auto [path_vertices, path_length] = compute_shortest_path_trace(*this, source, target);
    if (path_vertices.empty()) {
        REPORT_ERROR("No weighted shortest path exists between the selected vertices");
    }

    result.path_vertices = path_vertices;
    auto [path_arc_length, total_length] = compute_arc_length_coords(result.path_vertices);
    result.path_arc_length = path_arc_length;
    result.path_length = total_length;
    result.tube_radius = path_relative_radius * result.path_length;

    if (compute_excess) {
        if (std::isfinite(excess_tolerance)) {
            result.excess_tolerance = excess_tolerance;
        } else {
            result.excess_tolerance = result.tube_radius;
        }
    }

    const tubular_neighborhood_t tube_nbhd = compute_tubular_neighborhood_geodesic(
        result.path_vertices,
        result.tube_radius
    );
    result.tube_vertices = tube_nbhd.vertices;
    result.tube_geodesic_distances = tube_nbhd.geodesic_distances;

    const std::vector<double> dist_from_source =
        compute_shortest_path_distances_within_radius(*this, source, result.path_length);
    const std::vector<double> dist_from_target =
        compute_shortest_path_distances_within_radius(*this, target, result.path_length);

    std::unordered_set<size_t> trajectory_set(result.path_vertices.begin(), result.path_vertices.end());
    std::unordered_map<size_t, double> trajectory_coord_map;
    for (size_t i = 0; i < result.path_vertices.size(); ++i) {
        trajectory_coord_map[result.path_vertices[i]] = result.path_arc_length[i];
    }

    std::vector<double> corridor_initial_coords;
    corridor_initial_coords.reserve(result.tube_vertices.size());

    result.corridor_vertices.reserve(result.tube_vertices.size());
    result.t_balance.reserve(result.tube_vertices.size());
    result.distance_to_path.reserve(result.tube_vertices.size());
    result.excess.reserve(result.tube_vertices.size());
    if (compute_excess) {
        result.excess_vertices.reserve(result.tube_vertices.size());
    }

    for (size_t idx = 0; idx < result.tube_vertices.size(); ++idx) {
        const size_t vertex = result.tube_vertices[idx];
        const double du = dist_from_source[vertex];
        const double dv = dist_from_target[vertex];

        if (!std::isfinite(du) || !std::isfinite(dv)) {
            continue;
        }
        if (du > result.path_length + kTubeLensTol || dv > result.path_length + kTubeLensTol) {
            continue;
        }

        result.corridor_vertices.push_back(vertex);
        result.distance_to_path.push_back(result.tube_geodesic_distances[idx]);
        const double excess = du + dv - result.path_length;
        result.excess.push_back(excess);
        if (result.path_length > 0) {
            result.t_balance.push_back((result.path_length + du - dv) / (2.0 * result.path_length));
        } else {
            result.t_balance.push_back(0.0);
        }
        corridor_initial_coords.push_back(result.path_arc_length[tube_nbhd.nearest_traj_idx[idx]]);

        if (compute_excess) {
            if (excess <= result.excess_tolerance + kTubeLensTol) {
                result.excess_vertices.push_back(vertex);
            }
        }
    }

    int n_iterations = 0;
    double final_max_change = 0.0;
    result.harmonic_t = solve_harmonic_extension(
        result.corridor_vertices,
        trajectory_set,
        trajectory_coord_map,
        corridor_initial_coords,
        true,
        1000,
        1e-8,
        n_iterations,
        final_max_change
    );

    return result;
}

extern "C" SEXP S_compute_tube_lens_corridor(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_start,
    SEXP s_end,
    SEXP s_path_relative_radius,
    SEXP s_excess_tolerance,
    SEXP s_compute_excess
) {
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    const size_t start = static_cast<size_t>(Rf_asInteger(s_start));
    const size_t end = static_cast<size_t>(Rf_asInteger(s_end));
    const double path_relative_radius = Rf_asReal(s_path_relative_radius);
    const double excess_tolerance = Rf_asReal(s_excess_tolerance);
    const bool compute_excess = Rf_asLogical(s_compute_excess) == TRUE;

    const tube_lens_corridor_result_t result = graph.compute_tube_lens_corridor(
        start,
        end,
        path_relative_radius,
        excess_tolerance,
        compute_excess
    );

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 13));
    SEXP r_path_vertices = PROTECT(size_t_vector_to_r_int(result.path_vertices));
    SEXP r_path_arc_length = PROTECT(double_vector_to_r_num(result.path_arc_length));
    SEXP r_path_length = PROTECT(Rf_ScalarReal(result.path_length));
    SEXP r_tube_radius = PROTECT(Rf_ScalarReal(result.tube_radius));
    SEXP r_excess_tolerance = PROTECT(Rf_ScalarReal(
        std::isfinite(result.excess_tolerance) ? result.excess_tolerance : NA_REAL
    ));
    SEXP r_tube_vertices = PROTECT(size_t_vector_to_r_int(result.tube_vertices));
    SEXP r_tube_geodesic_distances = PROTECT(double_vector_to_r_num(result.tube_geodesic_distances));
    SEXP r_corridor_vertices = PROTECT(size_t_vector_to_r_int(result.corridor_vertices));
    SEXP r_t_balance = PROTECT(double_vector_to_r_num(result.t_balance));
    SEXP r_harmonic_t = PROTECT(double_vector_to_r_num(result.harmonic_t));
    SEXP r_distance_to_path = PROTECT(double_vector_to_r_num(result.distance_to_path));
    SEXP r_excess = PROTECT(double_vector_to_r_num(result.excess));
    SET_VECTOR_ELT(out, 0, r_path_vertices);
    SET_VECTOR_ELT(out, 1, r_path_arc_length);
    SET_VECTOR_ELT(out, 2, r_path_length);
    SET_VECTOR_ELT(out, 3, r_tube_radius);
    SET_VECTOR_ELT(out, 4, r_excess_tolerance);
    SET_VECTOR_ELT(out, 5, r_tube_vertices);
    SET_VECTOR_ELT(out, 6, r_tube_geodesic_distances);
    SET_VECTOR_ELT(out, 7, r_corridor_vertices);
    SET_VECTOR_ELT(out, 8, r_t_balance);
    SET_VECTOR_ELT(out, 9, r_harmonic_t);
    SET_VECTOR_ELT(out, 10, r_distance_to_path);
    SET_VECTOR_ELT(out, 11, r_excess);
    if (compute_excess) {
        SEXP r_excess_vertices = PROTECT(size_t_vector_to_r_int(result.excess_vertices));
        SET_VECTOR_ELT(out, 12, r_excess_vertices);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(out, 12, R_NilValue);
    }

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 13));
    SET_STRING_ELT(names, 0, Rf_mkChar("path.vertices"));
    SET_STRING_ELT(names, 1, Rf_mkChar("path.arc.length"));
    SET_STRING_ELT(names, 2, Rf_mkChar("path.length"));
    SET_STRING_ELT(names, 3, Rf_mkChar("tube.radius"));
    SET_STRING_ELT(names, 4, Rf_mkChar("excess.tolerance"));
    SET_STRING_ELT(names, 5, Rf_mkChar("tube.vertices"));
    SET_STRING_ELT(names, 6, Rf_mkChar("tube.geodesic.distances"));
    SET_STRING_ELT(names, 7, Rf_mkChar("corridor.vertices"));
    SET_STRING_ELT(names, 8, Rf_mkChar("t.balance"));
    SET_STRING_ELT(names, 9, Rf_mkChar("harmonic.t"));
    SET_STRING_ELT(names, 10, Rf_mkChar("distance.to.path"));
    SET_STRING_ELT(names, 11, Rf_mkChar("excess"));
    SET_STRING_ELT(names, 12, Rf_mkChar("excess.vertices"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(14);
    return out;
}
