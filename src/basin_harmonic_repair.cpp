#include "set_wgraph.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Basin-construction support retained in gflow. This is intentionally separate
// from the retired public harmonic-smoothing family.
void set_wgraph_t::perform_harmonic_repair(
    std::vector<double>& harmonic_predictions,
    const basin_t& absorbing_basin,
    const basin_t& absorbed_basin
) const {
    std::unordered_set<size_t> basin_vertices;
    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        basin_vertices.insert(v_info.vertex);
    }

    std::unordered_map<size_t, double> boundary_values;
    const size_t absorbing_extremum = absorbing_basin.extremum_vertex;
    if (absorbing_extremum < harmonic_predictions.size()) {
        boundary_values[absorbing_extremum] =
            harmonic_predictions[absorbing_extremum];
    }

    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        const size_t v = v_info.vertex;
        if (v == absorbed_basin.extremum_vertex ||
            v >= adjacency_list.size()) {
            continue;
        }
        for (const auto& edge : adjacency_list[v]) {
            if (basin_vertices.count(edge.vertex) == 0) {
                if (v < harmonic_predictions.size()) {
                    boundary_values[v] = harmonic_predictions[v];
                }
                break;
            }
        }
    }

    for (const auto& v_info : absorbing_basin.reachability_map.sorted_vertices) {
        const size_t v = v_info.vertex;
        if (basin_vertices.count(v) > 0 &&
            v != absorbed_basin.extremum_vertex &&
            v < harmonic_predictions.size()) {
            boundary_values[v] = harmonic_predictions[v];
        }
    }

    constexpr int max_iterations = 100;
    constexpr double tolerance = 1e-6;
    std::vector<double> new_values = harmonic_predictions;

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool converged = true;
        for (const auto& v_info :
             absorbed_basin.reachability_map.sorted_vertices) {
            const size_t v = v_info.vertex;
            if (boundary_values.count(v) > 0 ||
                v >= adjacency_list.size() ||
                v >= harmonic_predictions.size()) {
                continue;
            }

            double sum = 0.0;
            double weight_sum = 0.0;
            for (const auto& edge : adjacency_list[v]) {
                if (edge.vertex >= harmonic_predictions.size()) {
                    continue;
                }
                const double weight = 1.0 / (edge.weight + 1e-10);
                sum += harmonic_predictions[edge.vertex] * weight;
                weight_sum += weight;
            }
            if (weight_sum > 0) {
                const double weighted_average = sum / weight_sum;
                if (std::abs(new_values[v] - weighted_average) > tolerance) {
                    converged = false;
                }
                new_values[v] = weighted_average;
            }
        }

        for (const auto& v_info :
             absorbed_basin.reachability_map.sorted_vertices) {
            const size_t v = v_info.vertex;
            if (boundary_values.count(v) == 0 &&
                v < harmonic_predictions.size()) {
                harmonic_predictions[v] = new_values[v];
            }
        }
        if (converged) {
            break;
        }
    }

    const size_t absorbed_extremum = absorbed_basin.extremum_vertex;
    if (absorbed_extremum >= adjacency_list.size() ||
        absorbed_extremum >= harmonic_predictions.size()) {
        return;
    }

    double neighbor_average = 0.0;
    double weight_sum = 0.0;
    for (const auto& edge : adjacency_list[absorbed_extremum]) {
        if (edge.vertex >= harmonic_predictions.size()) {
            continue;
        }
        const double weight = 1.0 / (edge.weight + 1e-10);
        neighbor_average += harmonic_predictions[edge.vertex] * weight;
        weight_sum += weight;
    }

    if (weight_sum > 0) {
        neighbor_average /= weight_sum;
        double extrema_distance = 1.0;
        const auto distance_it =
            absorbing_basin.reachability_map.distances.find(absorbed_extremum);
        if (distance_it !=
            absorbing_basin.reachability_map.distances.end()) {
            extrema_distance = distance_it->second;
        } else {
            extrema_distance = std::max(
                1.0,
                compute_shortest_path_distance(
                    absorbing_basin.extremum_vertex,
                    absorbed_extremum
                )
            );
        }
        const double blend_factor = 1.0 / (1.0 + extrema_distance);
        harmonic_predictions[absorbed_extremum] =
            blend_factor * harmonic_predictions[absorbing_extremum] +
            (1.0 - blend_factor) * neighbor_average;
    }
}
