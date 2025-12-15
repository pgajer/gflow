/**
 * @file gfc.cpp
 * @brief Implementation of Gradient Flow Complex computation with refinement pipeline
 *
 * This file implements the GFC computation pipeline, including basin computation,
 * filtering, clustering, and expansion stages. The implementation leverages
 * existing set_wgraph_t methods for core operations while adding the refinement
 * logic specific to the association analysis framework.
 */

#include "gfc.hpp"
#include "set_wgraph.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <queue>
#include <limits>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

// For R integration
#include <R.h>

// ============================================================================
// Helper Function Implementations
// ============================================================================

Eigen::MatrixXd compute_overlap_distance_matrix(
    const std::vector<std::vector<size_t>>& basin_vertices
) {
    const size_t n_basins = basin_vertices.size();

    if (n_basins == 0) {
        return Eigen::MatrixXd(0, 0);
    }

    Eigen::MatrixXd dist(n_basins, n_basins);
    dist.setZero();

    // Precompute sets for efficient intersection
    std::vector<std::unordered_set<size_t>> basin_sets(n_basins);
    for (size_t i = 0; i < n_basins; ++i) {
        basin_sets[i].insert(basin_vertices[i].begin(), basin_vertices[i].end());
    }

    // Compute pairwise overlap distances
    for (size_t i = 0; i < n_basins; ++i) {
        const auto& set_i = basin_sets[i];
        const size_t size_i = set_i.size();

        for (size_t j = i + 1; j < n_basins; ++j) {
            const auto& vertices_j = basin_vertices[j];
            const size_t size_j = vertices_j.size();

            // Count intersection by iterating over smaller set
            size_t intersection_size = 0;
            if (size_i <= size_j) {
                for (size_t v : basin_vertices[i]) {
                    if (basin_sets[j].count(v)) {
                        ++intersection_size;
                    }
                }
            } else {
                for (size_t v : vertices_j) {
                    if (set_i.count(v)) {
                        ++intersection_size;
                    }
                }
            }

            // Szymkiewicz-Simpson overlap coefficient
            const size_t min_size = std::min(size_i, size_j);
            double d = (min_size > 0)
                ? 1.0 - static_cast<double>(intersection_size) / static_cast<double>(min_size)
                : 1.0;

            dist(i, j) = d;
            dist(j, i) = d;
        }
    }

    return dist;
}

std::vector<std::vector<int>> create_threshold_graph(
    const Eigen::MatrixXd& dist_matrix,
    double threshold
) {
    const int n = static_cast<int>(dist_matrix.rows());
    std::vector<std::vector<int>> adj_list(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && dist_matrix(i, j) < threshold) {
                adj_list[i].push_back(j);
            }
        }
    }

    return adj_list;
}

std::vector<int> find_connected_components(
    const std::vector<std::vector<int>>& adj_list
) {
    const int n = static_cast<int>(adj_list.size());

    if (n == 0) {
        return std::vector<int>();
    }

    std::vector<int> component(n, -1);
    int current_component = 0;

    for (int start = 0; start < n; ++start) {
        if (component[start] >= 0) {
            continue;  // Already assigned
        }

        // BFS from this vertex
        std::queue<int> q;
        q.push(start);
        component[start] = current_component;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v : adj_list[u]) {
                if (component[v] < 0) {
                    component[v] = current_component;
                    q.push(v);
                }
            }
        }

        ++current_component;
    }

    return component;
}

double compute_edge_length_threshold(
    const set_wgraph_t& graph,
    double quantile
) {
    // Collect all edge lengths
    std::vector<double> edge_lengths;
    const size_t n = graph.num_vertices();

    for (size_t u = 0; u < n; ++u) {
        for (const auto& edge : graph.adjacency_list[u]) {
            if (edge.vertex > u) {  // Count each edge once
                edge_lengths.push_back(edge.weight);
            }
        }
    }

    if (edge_lengths.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    // Sort and find quantile
    std::sort(edge_lengths.begin(), edge_lengths.end());

    const size_t idx = static_cast<size_t>(
        std::min(quantile, 1.0) * static_cast<double>(edge_lengths.size() - 1)
    );

    return edge_lengths[idx];
}

std::vector<int> expand_basins_to_cover(
    const set_wgraph_t& graph,
    const std::vector<std::vector<size_t>>& basin_vertices,
    size_t n_vertices
) {
    const int n_basins = static_cast<int>(basin_vertices.size());

    if (n_basins == 0) {
        return std::vector<int>(n_vertices, -1);
    }

    // Initialize assignment: -1 means unassigned
    std::vector<int> assignment(n_vertices, -1);

    // Mark vertices already in basins
    std::vector<bool> covered(n_vertices, false);
    for (int b = 0; b < n_basins; ++b) {
        for (size_t v : basin_vertices[b]) {
            if (v < n_vertices) {
                covered[v] = true;
                // If vertex is in multiple basins, assign to first one found
                if (assignment[v] < 0) {
                    assignment[v] = b;
                }
            }
        }
    }

    // Find uncovered vertices
    std::vector<size_t> uncovered;
    for (size_t v = 0; v < n_vertices; ++v) {
        if (!covered[v]) {
            uncovered.push_back(v);
        }
    }

    if (uncovered.empty()) {
        return assignment;
    }

    // For each uncovered vertex, find nearest basin via Dijkstra
    // This is done by running multi-source Dijkstra from all basin vertices

    // Initialize distances
    std::vector<double> dist(n_vertices, std::numeric_limits<double>::infinity());
    std::vector<int> nearest_basin(n_vertices, -1);

    // Priority queue: (distance, vertex)
    using pq_entry = std::pair<double, size_t>;
    std::priority_queue<pq_entry, std::vector<pq_entry>, std::greater<pq_entry>> pq;

    // Initialize with all basin vertices as sources
    for (int b = 0; b < n_basins; ++b) {
        for (size_t v : basin_vertices[b]) {
            if (v < n_vertices && dist[v] > 0.0) {
                dist[v] = 0.0;
                nearest_basin[v] = b;
                pq.push({0.0, v});
            }
        }
    }

    // Run Dijkstra
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) {
            continue;  // Stale entry
        }

        for (const auto& edge : graph.adjacency_list[u]) {
            size_t v = edge.vertex;
            double new_dist = dist[u] + edge.weight;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                nearest_basin[v] = nearest_basin[u];
                pq.push({new_dist, v});
            }
        }
    }

    // Assign uncovered vertices to nearest basin
    for (size_t v : uncovered) {
        assignment[v] = nearest_basin[v];
    }

    return assignment;
}

// ============================================================================
// Internal Pipeline Functions
// ============================================================================

namespace gfc_internal {

void filter_by_relvalue(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double min_rel_value_max,
    double max_rel_value_min,
    bool is_maximum
) {
    if (summaries.size() != basins.size()) {
        return;  // Safety check
    }

    std::vector<extremum_summary_t> filtered_summaries;
    std::vector<basin_compact_t> filtered_basins;

    for (size_t i = 0; i < summaries.size(); ++i) {
        bool keep = false;

        if (is_maximum) {
            keep = (summaries[i].rel_value >= min_rel_value_max);
        } else {
            keep = (summaries[i].rel_value <= max_rel_value_min);
        }

        if (keep) {
            filtered_summaries.push_back(std::move(summaries[i]));
            filtered_basins.push_back(std::move(basins[i]));
        }
    }

    summaries = std::move(filtered_summaries);
    basins = std::move(filtered_basins);
}

void cluster_and_merge_basins(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double overlap_threshold,
    bool is_maximum
) {
    const size_t n = basins.size();

    if (n <= 1) {
        return;  // Nothing to cluster
    }

    // Extract vertex sets from basins
    std::vector<std::vector<size_t>> basin_vertices(n);
    for (size_t i = 0; i < n; ++i) {
        basin_vertices[i] = basins[i].vertices;
    }

    // Compute overlap distance matrix
    Eigen::MatrixXd dist_matrix = compute_overlap_distance_matrix(basin_vertices);

    // Create threshold graph and find connected components
    std::vector<std::vector<int>> threshold_graph = create_threshold_graph(
        dist_matrix, overlap_threshold
    );
    std::vector<int> cluster_assignment = find_connected_components(threshold_graph);

    // Find number of clusters
    int n_clusters = 0;
    for (int c : cluster_assignment) {
        n_clusters = std::max(n_clusters, c + 1);
    }

    // Merge basins within each cluster
    std::vector<extremum_summary_t> merged_summaries;
    std::vector<basin_compact_t> merged_basins;

    for (int c = 0; c < n_clusters; ++c) {
        // Find all basins in this cluster
        std::vector<size_t> cluster_members;
        for (size_t i = 0; i < n; ++i) {
            if (cluster_assignment[i] == c) {
                cluster_members.push_back(i);
            }
        }

        if (cluster_members.empty()) {
            continue;
        }

        // Find representative (most extreme value)
        size_t rep_idx = cluster_members[0];
        for (size_t idx : cluster_members) {
            if (is_maximum) {
                if (summaries[idx].value > summaries[rep_idx].value) {
                    rep_idx = idx;
                }
            } else {
                if (summaries[idx].value < summaries[rep_idx].value) {
                    rep_idx = idx;
                }
            }
        }

        // Create merged basin
        basin_compact_t merged_basin = basins[rep_idx];

        if (cluster_members.size() > 1) {
            // Union all vertices from cluster members
            std::unordered_set<size_t> vertex_union;
            int max_hop = 0;

            for (size_t idx : cluster_members) {
                for (size_t v : basins[idx].vertices) {
                    vertex_union.insert(v);
                }
                max_hop = std::max(max_hop, basins[idx].max_hop_distance);
            }

            // Update merged basin vertices
            merged_basin.vertices.assign(vertex_union.begin(), vertex_union.end());
            std::sort(merged_basin.vertices.begin(), merged_basin.vertices.end());

            // Recompute hop distances (simplified: assign max_hop + 1 to new vertices)
            merged_basin.hop_distances.resize(merged_basin.vertices.size());
            std::unordered_set<size_t> original_vertices(
                basins[rep_idx].vertices.begin(),
                basins[rep_idx].vertices.end()
            );

            for (size_t i = 0; i < merged_basin.vertices.size(); ++i) {
                size_t v = merged_basin.vertices[i];
                // Find original hop distance if vertex was in representative basin
                auto it = std::find(basins[rep_idx].vertices.begin(),
                                    basins[rep_idx].vertices.end(), v);
                if (it != basins[rep_idx].vertices.end()) {
                    size_t orig_idx = std::distance(basins[rep_idx].vertices.begin(), it);
                    merged_basin.hop_distances[i] = basins[rep_idx].hop_distances[orig_idx];
                } else {
                    merged_basin.hop_distances[i] = max_hop + 1;
                }
            }

            merged_basin.max_hop_distance = max_hop + 1;
        }

        // Update summary for merged basin
        extremum_summary_t merged_summary = summaries[rep_idx];
        merged_summary.basin_size = static_cast<int>(merged_basin.vertices.size());
        merged_summary.hop_index = merged_basin.max_hop_distance;

        merged_summaries.push_back(std::move(merged_summary));
        merged_basins.push_back(std::move(merged_basin));
    }

    summaries = std::move(merged_summaries);
    basins = std::move(merged_basins);
}

void filter_by_geometry(
    std::vector<extremum_summary_t>& summaries,
    std::vector<basin_compact_t>& basins,
    double p_mean_hopk_dist_threshold,
    double p_deg_threshold,
    int min_basin_size
) {
    if (summaries.size() != basins.size()) {
        return;
    }

    std::vector<extremum_summary_t> filtered_summaries;
    std::vector<basin_compact_t> filtered_basins;

    for (size_t i = 0; i < summaries.size(); ++i) {
        bool keep = true;

        // Check geometric criteria
        if (summaries[i].mean_hopk_dist >= p_mean_hopk_dist_threshold) {
            keep = false;
        }
        if (summaries[i].deg_percentile >= p_deg_threshold) {
            keep = false;
        }
        if (summaries[i].basin_size < min_basin_size) {
            keep = false;
        }

        if (keep) {
            filtered_summaries.push_back(std::move(summaries[i]));
            filtered_basins.push_back(std::move(basins[i]));
        }
    }

    summaries = std::move(filtered_summaries);
    basins = std::move(filtered_basins);
}

std::vector<extremum_summary_t> compute_basin_summaries(
    const set_wgraph_t& graph,
    const std::vector<basin_compact_t>& basins,
    const std::vector<double>& y,
    double y_median,
    int hop_k
) {
    const size_t n_basins = basins.size();
    const size_t n_vertices = graph.num_vertices();

    std::vector<extremum_summary_t> summaries(n_basins);

    // Compute vertex degrees for percentile calculation
    std::vector<int> degrees(n_vertices);
    for (size_t v = 0; v < n_vertices; ++v) {
        degrees[v] = static_cast<int>(graph.adjacency_list[v].size());
    }

    // Sort degrees for percentile computation
    std::vector<int> sorted_degrees = degrees;
    std::sort(sorted_degrees.begin(), sorted_degrees.end());

    for (size_t b = 0; b < n_basins; ++b) {
        const auto& basin = basins[b];
        extremum_summary_t& summary = summaries[b];

        summary.vertex = basin.extremum_vertex;
        summary.value = basin.extremum_value;
        summary.is_maximum = basin.is_maximum;
        summary.basin_size = static_cast<int>(basin.vertices.size());
        summary.hop_index = basin.max_hop_distance;

        // Relative value
        summary.rel_value = (std::abs(y_median) > 1e-10)
            ? basin.extremum_value / y_median
            : 1.0;

        // Degree percentile
        int deg = degrees[basin.extremum_vertex];
        auto it = std::lower_bound(sorted_degrees.begin(), sorted_degrees.end(), deg);
        summary.deg_percentile = static_cast<double>(
            std::distance(sorted_degrees.begin(), it)
        ) / static_cast<double>(n_vertices);

        // Mean distance to neighbors (simplified computation)
        double sum_nbr_dist = 0.0;
        int n_nbrs = 0;
        for (const auto& edge : graph.adjacency_list[basin.extremum_vertex]) {
            sum_nbr_dist += edge.weight;
            ++n_nbrs;
        }
        summary.mean_nbrs_dist = (n_nbrs > 0) ? sum_nbr_dist / n_nbrs : 0.0;

        // Mean hop-k distance percentile (simplified: use mean neighbor distance as proxy)
        // A full implementation would compute actual hop-k neighborhoods
        summary.mean_hopk_dist = summary.mean_nbrs_dist;
    }

    // Convert mean_hopk_dist to percentiles
    if (n_basins > 1) {
        std::vector<double> hopk_values(n_basins);
        for (size_t b = 0; b < n_basins; ++b) {
            hopk_values[b] = summaries[b].mean_hopk_dist;
        }
        std::vector<double> sorted_hopk = hopk_values;
        std::sort(sorted_hopk.begin(), sorted_hopk.end());

        for (size_t b = 0; b < n_basins; ++b) {
            auto it = std::lower_bound(sorted_hopk.begin(), sorted_hopk.end(),
                                       hopk_values[b]);
            summaries[b].mean_hopk_dist = static_cast<double>(
                std::distance(sorted_hopk.begin(), it)
            ) / static_cast<double>(n_basins);
        }
    }

    return summaries;
}

std::vector<std::vector<int>> build_membership_vectors(
    const std::vector<basin_compact_t>& basins,
    size_t n_vertices
) {
    std::vector<std::vector<int>> membership(n_vertices);

    for (int b = 0; b < static_cast<int>(basins.size()); ++b) {
        for (size_t v : basins[b].vertices) {
            if (v < n_vertices) {
                membership[v].push_back(b);
            }
        }
    }

    return membership;
}

} // namespace gfc_internal

// ============================================================================
// Main GFC Computation
// ============================================================================

gfc_result_t compute_gfc(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const gfc_params_t& params,
    bool verbose
) {
    gfc_result_t result;
    result.n_vertices = graph.num_vertices();
    result.params = params;

    const size_t n = result.n_vertices;

    if (n == 0 || y.size() != n) {
        return result;
    }

    // Compute median of y
    std::vector<double> y_sorted = y;
    std::sort(y_sorted.begin(), y_sorted.end());
    result.y_median = y_sorted[n / 2];

    if (verbose) {
        Rprintf("GFC computation: %zu vertices, median = %.4f\n", n, result.y_median);
    }

    // ========================================================================
    // Step 1: Compute initial basins
    // ========================================================================

    if (verbose) {
        Rprintf("Step 1: Computing initial basins of attraction...\n");
    }

    // Find local extrema
    auto [lmin_vertices, lmax_vertices] = graph.find_nbr_extrema(y);

    // Compute basins for each extremum
    std::vector<basin_compact_t> max_basins;
    std::vector<basin_compact_t> min_basins;

    // Maximum basins
    for (size_t vertex : lmax_vertices) {
        auto basin = graph.find_local_extremum_bfs_basin(vertex, y, true);

        basin_compact_t compact;
        compact.extremum_vertex = vertex;
        compact.extremum_value = y[vertex];
        compact.is_maximum = true;

        // Extract vertices and hop distances from basin_t
        for (const auto& vi : basin.reachability_map.sorted_vertices) {
            compact.vertices.push_back(vi.vertex);
            compact.hop_distances.push_back(static_cast<int>(vi.distance));
        }

        compact.max_hop_distance = compact.hop_distances.empty()
            ? 0
            : *std::max_element(compact.hop_distances.begin(), compact.hop_distances.end());

        max_basins.push_back(std::move(compact));
    }

    // Minimum basins
    for (size_t vertex : lmin_vertices) {
        auto basin = graph.find_local_extremum_bfs_basin(vertex, y, false);

        basin_compact_t compact;
        compact.extremum_vertex = vertex;
        compact.extremum_value = y[vertex];
        compact.is_maximum = false;

        for (const auto& vi : basin.reachability_map.sorted_vertices) {
            compact.vertices.push_back(vi.vertex);
            compact.hop_distances.push_back(static_cast<int>(vi.distance));
        }

        compact.max_hop_distance = compact.hop_distances.empty()
            ? 0
            : *std::max_element(compact.hop_distances.begin(), compact.hop_distances.end());

        min_basins.push_back(std::move(compact));
    }

    // Compute initial summaries
    auto max_summaries = gfc_internal::compute_basin_summaries(
        graph, max_basins, y, result.y_median, params.hop_k
    );
    auto min_summaries = gfc_internal::compute_basin_summaries(
        graph, min_basins, y, result.y_median, params.hop_k
    );

    // Record initial counts
    result.stage_history.emplace_back(
        "initial",
        static_cast<int>(max_basins.size()),
        static_cast<int>(max_basins.size()),
        static_cast<int>(min_basins.size()),
        static_cast<int>(min_basins.size())
    );

    if (verbose) {
        Rprintf("  Found %zu maxima and %zu minima\n",
                max_basins.size(), min_basins.size());
    }

    // ========================================================================
    // Step 2: Filter by relative values
    // ========================================================================

    if (params.apply_relvalue_filter) {
        if (verbose) {
            Rprintf("Step 2: Filtering by relative values (max >= %.2f, min <= %.2f)...\n",
                    params.min_rel_value_max, params.max_rel_value_min);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::filter_by_relvalue(
            max_summaries, max_basins,
            params.min_rel_value_max, params.max_rel_value_min, true
        );
        gfc_internal::filter_by_relvalue(
            min_summaries, min_basins,
            params.min_rel_value_max, params.max_rel_value_min, false
        );

        result.stage_history.emplace_back(
            "relvalue",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Retained %zu maxima and %zu minima\n",
                    max_basins.size(), min_basins.size());
        }
    }

    // ========================================================================
    // Step 3: Cluster and merge maxima
    // ========================================================================

    if (params.apply_maxima_clustering && max_basins.size() > 1) {
        if (verbose) {
            Rprintf("Step 3: Clustering maxima (overlap threshold = %.2f)...\n",
                    params.max_overlap_threshold);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::cluster_and_merge_basins(
            max_summaries, max_basins,
            params.max_overlap_threshold, true
        );

        result.stage_history.emplace_back(
            "merge.max",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Result: %zu maxima after merging\n", max_basins.size());
        }
    }

    // ========================================================================
    // Step 4: Cluster and merge minima
    // ========================================================================

    if (params.apply_minima_clustering && min_basins.size() > 1) {
        if (verbose) {
            Rprintf("Step 4: Clustering minima (overlap threshold = %.2f)...\n",
                    params.min_overlap_threshold);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::cluster_and_merge_basins(
            min_summaries, min_basins,
            params.min_overlap_threshold, false
        );

        result.stage_history.emplace_back(
            "merge.min",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Result: %zu minima after merging\n", min_basins.size());
        }
    }

    // ========================================================================
    // Step 5: Filter by geometric characteristics
    // ========================================================================

    if (params.apply_geometric_filter) {
        if (verbose) {
            Rprintf("Step 5: Geometric filtering (hopk < %.2f, deg < %.2f, size >= %d)...\n",
                    params.p_mean_hopk_dist_threshold,
                    params.p_deg_threshold,
                    params.min_basin_size);
        }

        int n_max_before = static_cast<int>(max_basins.size());
        int n_min_before = static_cast<int>(min_basins.size());

        gfc_internal::filter_by_geometry(
            max_summaries, max_basins,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size
        );

        gfc_internal::filter_by_geometry(
            min_summaries, min_basins,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size
        );

        result.stage_history.emplace_back(
            "geometric",
            n_max_before, static_cast<int>(max_basins.size()),
            n_min_before, static_cast<int>(min_basins.size())
        );

        if (verbose) {
            Rprintf("  Retained %zu maxima and %zu minima\n",
                    max_basins.size(), min_basins.size());
        }
    }

    // ========================================================================
    // Step 6: Build membership vectors
    // ========================================================================

    result.max_membership = gfc_internal::build_membership_vectors(max_basins, n);
    result.min_membership = gfc_internal::build_membership_vectors(min_basins, n);

    // ========================================================================
    // Step 7: Expand basins to cover all vertices (optional)
    // ========================================================================

    if (params.expand_basins) {
        if (verbose) {
            Rprintf("Step 6: Expanding basins to cover all vertices...\n");
        }

        // Extract vertex sets
        std::vector<std::vector<size_t>> max_vertex_sets(max_basins.size());
        std::vector<std::vector<size_t>> min_vertex_sets(min_basins.size());

        for (size_t i = 0; i < max_basins.size(); ++i) {
            max_vertex_sets[i] = max_basins[i].vertices;
        }
        for (size_t i = 0; i < min_basins.size(); ++i) {
            min_vertex_sets[i] = min_basins[i].vertices;
        }

        result.expanded_max_assignment = expand_basins_to_cover(
            graph, max_vertex_sets, n
        );
        result.expanded_min_assignment = expand_basins_to_cover(
            graph, min_vertex_sets, n
        );

        // Count coverage
        int n_covered_max = 0, n_covered_min = 0;
        for (size_t v = 0; v < n; ++v) {
            if (result.expanded_max_assignment[v] >= 0) ++n_covered_max;
            if (result.expanded_min_assignment[v] >= 0) ++n_covered_min;
        }

        result.stage_history.emplace_back(
            "expand",
            static_cast<int>(max_basins.size()), n_covered_max,
            static_cast<int>(min_basins.size()), n_covered_min
        );

        if (verbose) {
            Rprintf("  Max coverage: %d/%zu, Min coverage: %d/%zu\n",
                    n_covered_max, n, n_covered_min, n);
        }
    }

    // ========================================================================
    // Finalize result
    // ========================================================================

    result.max_basins = std::move(max_basins);
    result.min_basins = std::move(min_basins);
    result.max_summaries = std::move(max_summaries);
    result.min_summaries = std::move(min_summaries);

    if (verbose) {
        Rprintf("GFC computation complete: %zu maxima, %zu minima\n",
                result.max_basins.size(), result.min_basins.size());
    }

    return result;
}

std::vector<gfc_result_t> compute_gfc_matrix(
    const set_wgraph_t& graph,
    const Eigen::MatrixXd& Y,
    const gfc_params_t& params,
    int n_cores,
    bool verbose
) {
    const int n = static_cast<int>(Y.rows());
    const int p = static_cast<int>(Y.cols());

    if (verbose) {
        Rprintf("GFC matrix computation: %d functions, %d vertices\n", p, n);
#ifdef _OPENMP
        Rprintf("  OpenMP threads: %d\n", n_cores);
#else
        Rprintf("  OpenMP: not available (sequential execution)\n");
#endif
    }

    std::vector<gfc_result_t> results(p);

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    int progress_interval = std::max(1, p / 20);

    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
    for (int j = 0; j < p; ++j) {
        // Extract column j
        std::vector<double> y_j(n);
        for (int i = 0; i < n; ++i) {
            y_j[i] = Y(i, j);
        }

        // Compute GFC for this column (verbose = false in parallel)
        results[j] = compute_gfc(graph, y_j, params, false);

        // Progress reporting (only in sequential mode)
        #pragma omp critical
        {
            if (verbose && n_cores == 1 && (j + 1) % progress_interval == 0) {
                Rprintf("  Processed %d/%d functions (%.0f%%)\n",
                        j + 1, p, 100.0 * (j + 1) / p);
            }
        }
    }

    if (verbose) {
        Rprintf("  Complete.\n");
    }

    return results;
}
