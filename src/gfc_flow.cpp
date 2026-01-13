/**
 * @file gfc_flow.cpp
 * @brief Implementation of trajectory-based Gradient Flow Complex computation
 *
 * This implementation preserves ALL extrema and basins, marking spurious ones
 * rather than deleting them. This enables complete trajectory analysis and
 * harmonic repair of spurious regions.
 *
 * Key changes from original:
 * 1. All basins stored in *_basins_all vectors
 * 2. Filtering marks is_spurious=true instead of removing
 * 3. Labels assigned to all extrema (m/M for retained, sm/sM for spurious)
 * 4. Trajectories track endpoint spurious status
 * 5. Membership computed for both retained and spurious basins
 */

#include "gfc_flow.hpp"
#include "gfc.hpp"
#include "set_wgraph.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>

// ============================================================================
// Helper Functions for Filtering with Spurious Tracking
// ============================================================================

namespace gfc_flow_internal {

    /**
     * @brief Mark basins as spurious based on relative value filter
     *
     * Instead of removing basins, marks them as spurious with appropriate stage.
     */
    void mark_spurious_by_relvalue(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        double min_rel_value_max,
        double max_rel_value_min,
        bool is_maximum
        ) {
        for (size_t i = 0; i < summaries.size(); ++i) {
            if (summaries[i].is_spurious) continue;  // Already marked
        
            bool should_filter = false;
            if (is_maximum) {
                should_filter = (summaries[i].rel_value < min_rel_value_max);
            } else {
                should_filter = (summaries[i].rel_value > max_rel_value_min);
            }
        
            if (should_filter) {
                summaries[i].is_spurious = true;
                summaries[i].filter_stage = filter_stage_t::RELVALUE;
                basins[i].is_spurious = true;
                basins[i].filter_stage = filter_stage_t::RELVALUE;
            }
        }
    }

    /**
     * @brief Mark basins as spurious based on geometric filter
     */
    void mark_spurious_by_geometry(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        double p_mean_nbrs_dist_threshold,
        double p_mean_hopk_dist_threshold,
        double p_deg_threshold,
        int min_basin_size
        ) {
        for (size_t i = 0; i < summaries.size(); ++i) {
            if (summaries[i].is_spurious) continue;  // Already marked
        
            bool should_filter = false;
        
            // Check basin size
            if (summaries[i].basin_size < min_basin_size) {
                should_filter = true;
            }
        
            if (summaries[i].p_mean_nbrs_dist > p_mean_nbrs_dist_threshold) {
                should_filter = true;
            }

            if (summaries[i].p_mean_hopk_dist > p_mean_hopk_dist_threshold) {
                should_filter = true;
            }
        
            if (summaries[i].deg_percentile > p_deg_threshold) {
                should_filter = true;
            }
        
            if (should_filter) {
                summaries[i].is_spurious = true;
                summaries[i].filter_stage = filter_stage_t::GEOMETRIC;
                basins[i].is_spurious = true;
                basins[i].filter_stage = filter_stage_t::GEOMETRIC;
            }
        }
    }

    Eigen::MatrixXd get_overlap_distance_matrix(std::vector<basin_extended_t>& basins) {
        const size_t n = basins.size();

        // Extract vertex sets from basins
        std::vector<std::vector<size_t>> basin_vertices(n);
        for (size_t i = 0; i < n; ++i) {
            basin_vertices[i] = basins[i].vertices;
        }

        // Compute overlap distance matrix
        Eigen::MatrixXd dist_matrix = compute_overlap_distance_matrix(basin_vertices);

        return dist_matrix;
    }

    void mark_spurious_by_min_basin_size(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        int min_basin_size
        ) {
        if (min_basin_size <= 0) return;

        for (size_t i = 0; i < summaries.size(); ++i) {
            if (summaries[i].is_spurious) continue;

            if (summaries[i].basin_size < min_basin_size) {
                summaries[i].is_spurious = true;
                summaries[i].filter_stage = filter_stage_t::MIN_BASIN_SIZE;
                basins[i].is_spurious = true;
                basins[i].filter_stage = filter_stage_t::MIN_BASIN_SIZE;
            }
        }
    }

    void mark_spurious_by_min_n_trajectories(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        const std::unordered_map<size_t, int>& n_traj_by_extremum,
        int min_n_trajectories
        ) {
        if (min_n_trajectories <= 0) return;

        for (size_t i = 0; i < summaries.size(); ++i) {
            if (summaries[i].is_spurious) continue;

            size_t ext = basins[i].extremum_vertex;
            auto it = n_traj_by_extremum.find(ext);
            int n_traj = (it == n_traj_by_extremum.end()) ? 0 : it->second;

            if (n_traj < min_n_trajectories) {
                summaries[i].is_spurious = true;
                summaries[i].filter_stage = filter_stage_t::MIN_N_TRAJ;
                basins[i].is_spurious = true;
                basins[i].filter_stage = filter_stage_t::MIN_N_TRAJ;
            }
        }
    }

    /**
     * @brief Cluster and merge basins based on overlap, marking merged ones as spurious
     *
     * Uses the same clustering algorithm as gfc_internal::cluster_and_merge_basins():
     * - Computes overlap distance matrix using Szymkiewicz-Simpson coefficient
     * - Creates threshold graph and finds connected components
     * - For each cluster, selects representative (most extreme value)
     * - Merges all vertices from cluster members into representative basin
     * - Marks non-representative basins as spurious with tracking info
     *
     * @param summaries Vector of extremum summaries (modified in place)
     * @param basins Vector of basins (modified in place)
     * @param overlap_threshold Distance threshold for clustering (edges where d < threshold)
     * @param is_maximum True for maxima (pick highest), false for minima (pick lowest)
     * @return Overlap distance matrix for the active (non-spurious) basins at time of call
     */
    Eigen::MatrixXd cluster_and_merge_basins(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        double overlap_threshold,
        bool is_maximum
        ) {
        // Get indices of currently non-spurious basins
        std::vector<int> active_indices;
        for (size_t i = 0; i < basins.size(); ++i) {
            if (!basins[i].is_spurious) {
                active_indices.push_back(static_cast<int>(i));
            }
        }

        const size_t n = active_indices.size();

        if (n <= 1) {
            // Return empty or 1x1 matrix
            Eigen::MatrixXd dist(n, n);
            dist.setZero();
            return dist;
        }

        // Extract vertex vectors for active basins (for compute_overlap_distance_matrix)
        std::vector<std::vector<size_t>> basin_vertices(n);
        for (size_t i = 0; i < n; ++i) {
            basin_vertices[i] = basins[active_indices[i]].vertices;
        }

        // Compute overlap distance matrix using the standard function
        Eigen::MatrixXd dist_matrix = compute_overlap_distance_matrix(basin_vertices);

        // Create threshold graph and find connected components
        std::vector<std::vector<int>> threshold_graph = create_threshold_graph(
            dist_matrix, overlap_threshold);
        std::vector<int> cluster_assignment = find_connected_components(threshold_graph);

        // Find number of clusters
        int n_clusters = 0;
        for (int c : cluster_assignment) {
            n_clusters = std::max(n_clusters, c + 1);
        }

        // Group basins by cluster
        std::vector<std::vector<size_t>> clusters(n_clusters);
        for (size_t i = 0; i < n; ++i) {
            clusters[cluster_assignment[i]].push_back(i);
        }

        // Process each cluster
        for (int c = 0; c < n_clusters; ++c) {
            const auto& cluster_members = clusters[c];

            if (cluster_members.size() <= 1) {
                continue;  // Nothing to merge
            }

            // Find representative (most extreme value)
            size_t rep_local = cluster_members[0];
            int rep_global = active_indices[rep_local];
            double rep_value = summaries[rep_global].value;

            for (size_t local_idx : cluster_members) {
                int global_idx = active_indices[local_idx];
                double val = summaries[global_idx].value;
                if ((is_maximum && val > rep_value) || (!is_maximum && val < rep_value)) {
                    rep_value = val;
                    rep_local = local_idx;
                    rep_global = global_idx;
                }
            }

            // Collect hop distance info from all member basins BEFORE modifying anything
            std::unordered_map<size_t, int> vertex_min_hop;
            int max_hop = 0;

            for (size_t local_idx : cluster_members) {
                int global_idx = active_indices[local_idx];
                const auto& member_basin = basins[global_idx];
                for (size_t vi = 0; vi < member_basin.vertices.size(); ++vi) {
                    size_t v = member_basin.vertices[vi];
                    int hop = member_basin.hop_distances[vi];
                    auto it = vertex_min_hop.find(v);
                    if (it == vertex_min_hop.end() || hop < it->second) {
                        vertex_min_hop[v] = hop;
                    }
                    max_hop = std::max(max_hop, hop);
                }
            }

            // Merge vertices from all cluster members into representative
            std::unordered_set<size_t> vertex_union;
            for (size_t local_idx : cluster_members) {
                int global_idx = active_indices[local_idx];
                for (size_t v : basins[global_idx].vertices) {
                    vertex_union.insert(v);
                }
            }

            // Update representative basin with merged vertices
            basin_extended_t& rep_basin = basins[rep_global];

            // Update vertices
            rep_basin.vertices.assign(vertex_union.begin(), vertex_union.end());
            std::sort(rep_basin.vertices.begin(), rep_basin.vertices.end());

            // Assign hop distances using pre-collected info
            rep_basin.hop_distances.resize(rep_basin.vertices.size());
            int new_max_hop = 0;
            for (size_t i = 0; i < rep_basin.vertices.size(); ++i) {
                size_t v = rep_basin.vertices[i];
                auto it = vertex_min_hop.find(v);
                if (it != vertex_min_hop.end()) {
                    rep_basin.hop_distances[i] = it->second;
                } else {
                    rep_basin.hop_distances[i] = max_hop + 1;
                }
                new_max_hop = std::max(new_max_hop, rep_basin.hop_distances[i]);
            }
            rep_basin.max_hop_distance = new_max_hop;

            // Update representative summary
            summaries[rep_global].basin_size = static_cast<int>(rep_basin.vertices.size());
            summaries[rep_global].hop_index = rep_basin.max_hop_distance;

            // Mark non-representative basins as spurious
            for (size_t local_idx : cluster_members) {
                int global_idx = active_indices[local_idx];
                if (global_idx != rep_global) {
                    summaries[global_idx].is_spurious = true;
                    summaries[global_idx].filter_stage = filter_stage_t::CLUSTER_MERGE;
                    summaries[global_idx].merged_into = rep_global;
                    basins[global_idx].is_spurious = true;
                    basins[global_idx].filter_stage = filter_stage_t::CLUSTER_MERGE;
                    basins[global_idx].merged_into = rep_global;
                }
            }
        }

        return dist_matrix;
    }

    /**
     * @brief Cluster and merge basins, marking merged ones as spurious
     *
     * This function was used initially in place of cluster_and_merge_basins.
     */
    void cluster_and_mark_spurious(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        double overlap_threshold,
        bool is_maximum
        ) {
        // Get indices of currently non-spurious basins
        std::vector<int> active_indices;
        for (size_t i = 0; i < basins.size(); ++i) {
            if (!basins[i].is_spurious) {
                active_indices.push_back(static_cast<int>(i));
            }
        }
    
        if (active_indices.size() <= 1) return;
    
        // Build vertex sets for active basins
        std::vector<std::unordered_set<size_t>> basin_sets(active_indices.size());
        for (size_t i = 0; i < active_indices.size(); ++i) {
            int idx = active_indices[i];
            basin_sets[i].insert(basins[idx].vertices.begin(), basins[idx].vertices.end());
        }
    
        // Compute overlap distance matrix
        size_t n = active_indices.size();
        Eigen::MatrixXd dist(n, n);
        dist.setZero();
    
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                // Count intersection
                size_t intersection = 0;
                for (size_t v : basin_sets[i]) {
                    if (basin_sets[j].count(v)) ++intersection;
                }
            
                // Overlap distance = 1 - |A ∩ B| / min(|A|, |B|)
                size_t min_size = std::min(basin_sets[i].size(), basin_sets[j].size());
                double overlap = (min_size > 0) ?
                    static_cast<double>(intersection) / min_size : 0.0;
                double d = 1.0 - overlap;
            
                dist(i, j) = d;
                dist(j, i) = d;
            }
        }
    
        // Single-linkage clustering
        std::vector<int> cluster_id(n);
        std::iota(cluster_id.begin(), cluster_id.end(), 0);
    
        // Find pairs to merge (distance < threshold means overlap > 1-threshold)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (dist(i, j) < overlap_threshold) {
                    // Merge clusters
                    int old_cluster = cluster_id[j];
                    int new_cluster = cluster_id[i];
                    for (size_t k = 0; k < n; ++k) {
                        if (cluster_id[k] == old_cluster) {
                            cluster_id[k] = new_cluster;
                        }
                    }
                }
            }
        }
    
        // For each cluster, pick representative (highest/lowest value) and mark others spurious
        std::map<int, std::vector<int>> clusters;
        for (size_t i = 0; i < n; ++i) {
            clusters[cluster_id[i]].push_back(static_cast<int>(i));
        }
    
        for (auto& [cid, members] : clusters) {
            if (members.size() <= 1) continue;
        
            // Find representative: highest value for maxima, lowest for minima
            int rep_local = members[0];
            double rep_value = summaries[active_indices[rep_local]].value;
        
            for (int m : members) {
                double val = summaries[active_indices[m]].value;
                if ((is_maximum && val > rep_value) || (!is_maximum && val < rep_value)) {
                    rep_value = val;
                    rep_local = m;
                }
            }
        
            // Mark non-representatives as spurious
            int rep_global = active_indices[rep_local];
            for (int m : members) {
                if (m != rep_local) {
                    int global_idx = active_indices[m];
                    summaries[global_idx].is_spurious = true;
                    summaries[global_idx].filter_stage = filter_stage_t::CLUSTER_MERGE;
                    summaries[global_idx].merged_into = rep_global;
                    basins[global_idx].is_spurious = true;
                    basins[global_idx].filter_stage = filter_stage_t::CLUSTER_MERGE;
                    basins[global_idx].merged_into = rep_global;
                }
            }
        }
    }

    /**
     * @brief Assign labels to all extrema
     *
     * Retained: m1, m2, ... / M1, M2, ...
     * Spurious: sm1, sm2, ... / sM1, sM2, ...
     */
    void assign_labels(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        bool is_maximum
        ) {
        std::string retained_prefix = is_maximum ? "M" : "m";
        std::string spurious_prefix = is_maximum ? "sM" : "sm";

        int retained_count = 0;
        int spurious_count = 0;

        // Sort by basin size (largest first) for consistent labeling.
        // Tie-breaks:
        //   - maxima: higher value first
        //   - minima: lower value first
        //   - finally: smaller vertex index first (deterministic)
        std::vector<size_t> order(summaries.size());
        std::iota(order.begin(), order.end(), 0);

        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {

            int ba = summaries[a].basin_size;
            int bb = summaries[b].basin_size;

            if (ba != bb) return ba > bb;  // largest basin first

            double va = summaries[a].value;
            double vb = summaries[b].value;

            if (va != vb) {
                if (is_maximum) {
                    return va > vb;  // maxima: higher first
                } else {
                    return va < vb;  // minima: lower first
                }
            }

            // last tie-break: vertex id (stable, deterministic)
            return basins[a].extremum_vertex < basins[b].extremum_vertex;
        });

        for (size_t idx : order) {
            std::string label;
            if (summaries[idx].is_spurious) {
                ++spurious_count;
                label = spurious_prefix + std::to_string(spurious_count);
            } else {
                ++retained_count;
                label = retained_prefix + std::to_string(retained_count);
            }
            summaries[idx].label = label;
            basins[idx].label = label;
        }
    }

    #if 0
    void assign_labels(
        std::vector<extremum_summary_extended_t>& summaries,
        std::vector<basin_extended_t>& basins,
        bool is_maximum
        ) {
        std::string retained_prefix = is_maximum ? "M" : "m";
        std::string spurious_prefix = is_maximum ? "sM" : "sm";
    
        int retained_count = 0;
        int spurious_count = 0;
    
        // Sort by value for consistent labeling
        std::vector<size_t> order(summaries.size());
        std::iota(order.begin(), order.end(), 0);
    
        if (is_maximum) {
            // Maxima: highest value first
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return summaries[a].value > summaries[b].value;
            });
        } else {
            // Minima: lowest value first
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return summaries[a].value < summaries[b].value;
            });
        }
    
        for (size_t idx : order) {
            std::string label;
            if (summaries[idx].is_spurious) {
                ++spurious_count;
                label = spurious_prefix + std::to_string(spurious_count);
            } else {
                ++retained_count;
                label = retained_prefix + std::to_string(retained_count);
            }
            summaries[idx].label = label;
            basins[idx].label = label;
        }
    }
    #endif

    /**
     * @brief Build membership vectors for all basins
     */
    std::vector<std::vector<int>> build_membership_all(
        const std::vector<basin_extended_t>& basins,
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

    /**
     * @brief Build membership vectors for retained basins only
     *
     * Returns indices into the retained_indices vector, not into basins_all
     */
    std::vector<std::vector<int>> build_membership_retained(
        const std::vector<basin_extended_t>& basins,
        const std::vector<int>& retained_indices,
        size_t n_vertices
        ) {
        std::vector<std::vector<int>> membership(n_vertices);
    
        for (int r = 0; r < static_cast<int>(retained_indices.size()); ++r) {
            int b = retained_indices[r];
            for (size_t v : basins[b].vertices) {
                if (v < n_vertices) {
                    membership[v].push_back(r);  // Index into retained_indices
                }
            }
        }
    
        return membership;
    }

} // namespace gfc_flow_internal

// ============================================================================
// Member Function Implementations for set_wgraph_t
// ============================================================================

int set_wgraph_t::check_nbr_extremum_type(
    size_t vertex,
    const std::vector<double>& y
) const {
    if (vertex >= adjacency_list.size()) {
        return 0;
    }

    const auto& neighbors = adjacency_list[vertex];
    if (neighbors.empty()) {
        return 0;
    }

    bool is_max = true;
    bool is_min = true;
    double y_v = y[vertex];

    for (const auto& edge : neighbors) {
        double y_u = y[edge.vertex];
        if (y_u >= y_v) is_max = false;
        if (y_u <= y_v) is_min = false;
        if (!is_max && !is_min) return 0;
    }

    if (is_max) return 1;
    if (is_min) return -1;
    return 0;
}

static int select_next_vertex(
    const set_wgraph_t& graph,
    size_t current,
    const std::vector<double>& y,
    bool ascending,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    const edge_weight_map_t& edge_length_weights,
    double edge_length_thld
) {
    const size_t n = graph.adjacency_list.size();
    if (current >= n || y.size() != n) return -1;

    bool needs_density = (modulation == gflow_modulation_t::DENSITY ||
                          modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_density && density.size() != n) {
        modulation = (modulation == gflow_modulation_t::DENSITY) ?
            gflow_modulation_t::NONE : gflow_modulation_t::EDGELEN;
        needs_density = false;
    }

    bool needs_edge_weights = (modulation == gflow_modulation_t::EDGELEN ||
                               modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_edge_weights && edge_length_weights.empty()) {
        modulation = (modulation == gflow_modulation_t::EDGELEN) ?
            gflow_modulation_t::NONE : gflow_modulation_t::DENSITY;
        needs_edge_weights = (modulation == gflow_modulation_t::EDGELEN ||
                              modulation == gflow_modulation_t::DENSITY_EDGELEN);
        needs_density = (modulation == gflow_modulation_t::DENSITY ||
                         modulation == gflow_modulation_t::DENSITY_EDGELEN);
    }

    size_t best_next = INVALID_VERTEX;

    // CLOSEST: prefer shortest ascending neighbor with edge <= thld, else any ascending neighbor
    if (modulation == gflow_modulation_t::CLOSEST) {
        double min_short = std::numeric_limits<double>::infinity();
        double min_any   = std::numeric_limits<double>::infinity();
        size_t best_short = INVALID_VERTEX;
        size_t best_any = INVALID_VERTEX;

        for (const auto& edge : graph.adjacency_list[current]) {
            size_t u = edge.vertex;
            double edge_len = edge.weight;

            double delta_y = y[u] - y[current];
            if (ascending && delta_y <= 0) continue;
            if (!ascending && delta_y >= 0) continue;

            if (edge_len < min_any) {
                min_any = edge_len;
                best_any = u;
            }
            if (edge_len <= edge_length_thld && edge_len < min_short) {
                min_short = edge_len;
                best_short = u;
            }
        }

        best_next = (best_short != INVALID_VERTEX) ? best_short : best_any;
        return (best_next == INVALID_VERTEX) ? -1 : static_cast<int>(best_next);
    }

    // Score-based modulations
    double best_score = ascending ? -std::numeric_limits<double>::infinity()
                                  :  std::numeric_limits<double>::infinity();

    for (const auto& edge : graph.adjacency_list[current]) {
        size_t u = edge.vertex;
        double edge_len = edge.weight;

        // For non-CLOSEST modulations, long edges are not allowed
        if (edge_len > edge_length_thld) continue;

        double delta_y = y[u] - y[current];
        if (ascending && delta_y <= 0) continue;
        if (!ascending && delta_y >= 0) continue;

        double target_dens = needs_density ? density[u] : 1.0;

        double edge_weight = 1.0;
        if (needs_edge_weights) {
            auto it_v = edge_length_weights.find(current);
            if (it_v != edge_length_weights.end()) {
                auto it_u = it_v->second.find(u);
                if (it_u != it_v->second.end()) {
                    edge_weight = it_u->second;
                }
            }
        }

        double score = compute_modulated_score(delta_y, edge_weight, target_dens, modulation);

        bool is_better = ascending ? (score > best_score) : (score < best_score);
        if (is_better) {
            best_score = score;
            best_next = u;
        }
    }

    return (best_next == INVALID_VERTEX) ? -1 : static_cast<int>(best_next);
}

/**
 * @brief Follow gradient trajectory from a starting vertex to a local extremum
 *
 * Traces a path through the graph following the gradient direction
 * (ascending or descending) with optional modulation. The trajectory continues
 * until reaching a true local extremum where no further progress can be made.
 *
 * Modulation options:
 * - NONE: Steepest ascent/descent (max |Δy|)
 * - DENSITY: Density-weighted steepest (max ρ(u) · Δy)
 * - EDGELEN: Edge-length-weighted steepest (max w(d) · Δy)
 * - DENSITY_EDGELEN: Combined density and edge-length weighting
 * - CLOSEST: Lexicographic rule - among ascending neighbors, pick closest
 *
 * The CLOSEST modulation implements a lexicographic ordering:
 *   Level 1: Filter to A(v) = {u : y(u) > y(v)} (ascending) or
 *                    D(v) = {u : y(u) < y(v)} (descending)
 *   Level 2: Among the filtered set, select the neighbor with minimum edge weight
 *
 * This approach minimizes basin-jumping by preferring small steps over steep ones.
 */
std::vector<size_t> set_wgraph_t::follow_gradient_trajectory(
    size_t start_vertex,
    const std::vector<double>& y,
    bool ascending,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    const edge_weight_map_t& edge_length_weights,
    double edge_length_thld,
    size_t max_length
) const {
    std::vector<size_t> trajectory;
    const size_t n = adjacency_list.size();

    // Validate inputs
    if (start_vertex >= n || y.size() != n) {
        return trajectory;
    }

    // Check if density is required but not provided
    bool needs_density = (modulation == gflow_modulation_t::DENSITY ||
                          modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_density && density.size() != n) {
        // Fall back to NONE or EDGELEN
        if (modulation == gflow_modulation_t::DENSITY) {
            modulation = gflow_modulation_t::NONE;
        } else {
            modulation = gflow_modulation_t::EDGELEN;
        }
    }

    // Check if edge length weights are required but not provided
    // Note: CLOSEST does NOT need precomputed edge_length_weights - it uses raw edge.weight
    bool needs_edge_weights = (modulation == gflow_modulation_t::EDGELEN ||
                               modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_edge_weights && edge_length_weights.empty()) {
        // Fall back to NONE or DENSITY
        if (modulation == gflow_modulation_t::EDGELEN) {
            modulation = gflow_modulation_t::NONE;
        } else {
            modulation = gflow_modulation_t::DENSITY;
        }
    }

    // Track vertices in this trajectory to avoid self-loops
    std::unordered_set<size_t> trajectory_set;

    trajectory.push_back(start_vertex);
    trajectory_set.insert(start_vertex);

    size_t current = start_vertex;

    while (trajectory.size() < max_length) {
        size_t best_next = INVALID_VERTEX;

        // ====================================================================
        // CLOSEST modulation: Lexicographic selection
        // ====================================================================
        if (modulation == gflow_modulation_t::CLOSEST) {
            double min_distance_short = std::numeric_limits<double>::infinity();
            double min_distance_any = std::numeric_limits<double>::infinity();
            size_t best_short = INVALID_VERTEX;
            size_t best_any = INVALID_VERTEX;

            for (const auto& edge : adjacency_list[current]) {
                size_t u = edge.vertex;
                double edge_len = edge.weight;

                // Skip if already in this trajectory (avoid loops)
                if (trajectory_set.count(u) > 0) {
                    continue;
                }

                double delta_y = y[u] - y[current];

                // Filter to strictly ascending/descending neighbors
                if (ascending && delta_y <= 0) {
                    continue;
                }
                if (!ascending && delta_y >= 0) {
                    continue;
                }

                // Track best among ALL valid neighbors (fallback)
                if (edge_len < min_distance_any) {
                    min_distance_any = edge_len;
                    best_any = u;
                }

                // Track best among SHORT valid neighbors (preferred)
                if (edge_len <= edge_length_thld && edge_len < min_distance_short) {
                    min_distance_short = edge_len;
                    best_short = u;
                }
            }

            // Prefer short edge, but use any valid edge if no short ones exist
            best_next = (best_short != INVALID_VERTEX) ? best_short : best_any;
        }

        // ====================================================================
        // Other modulations: Score-based selection
        // ====================================================================
        else {
            double best_score = ascending ? -std::numeric_limits<double>::infinity()
                                          : std::numeric_limits<double>::infinity();

            for (const auto& edge : adjacency_list[current]) {
                size_t u = edge.vertex;
                double edge_len = edge.weight;

                // Skip if edge is too long
                if (edge_len > edge_length_thld) {
                    continue;
                }

                // Skip if already in this trajectory (avoid loops within the trajectory)
                if (trajectory_set.count(u) > 0) {
                    continue;
                }

                double delta_y = y[u] - y[current];

                // Skip if not moving in desired direction
                if (ascending && delta_y <= 0) {
                    continue;
                }
                if (!ascending && delta_y >= 0) {
                    continue;
                }

                // Get density for modulation (default 1.0 if not used)
                double target_dens = needs_density ? density[u] : 1.0;

                // Get edge length weight for modulation (default 1.0 if not used)
                double edge_weight = 1.0;
                if (needs_edge_weights) {
                    // Look up the precomputed edge length weight
                    auto it_v = edge_length_weights.find(current);
                    if (it_v != edge_length_weights.end()) {
                        auto it_u = it_v->second.find(u);
                        if (it_u != it_v->second.end()) {
                            edge_weight = it_u->second;
                        }
                    }
                }

                // Compute modulated score
                double score = compute_modulated_score(delta_y, edge_weight, target_dens, modulation);

                // Check if this is better
                bool is_better = ascending ? (score > best_score) : (score < best_score);
                if (is_better) {
                    best_score = score;
                    best_next = u;
                }
            }
        }

        // If no valid neighbor found, we've reached a local extremum
        if (best_next == INVALID_VERTEX) {
            break;
        }

        // Move to best neighbor
        trajectory.push_back(best_next);
        trajectory_set.insert(best_next);
        current = best_next;
    }

    return trajectory;
}

gflow_trajectory_t set_wgraph_t::join_trajectories_at_vertex(
    size_t vertex,
    const std::vector<double>& y,
    gflow_modulation_t modulation,
    const std::vector<double>& density,
    const edge_weight_map_t& edge_length_weights,
    double edge_length_thld,
    size_t max_length
) const {
    gflow_trajectory_t result;
    const size_t n = adjacency_list.size();

    if (vertex >= n) return result;

    std::vector<size_t> desc_traj = follow_gradient_trajectory(
        vertex, y, false, modulation, density, edge_length_weights, edge_length_thld, max_length
    );

    std::vector<size_t> asc_traj = follow_gradient_trajectory(
        vertex, y, true, modulation, density, edge_length_weights, edge_length_thld, max_length
    );

    std::reverse(desc_traj.begin(), desc_traj.end());

    result.vertices = std::move(desc_traj);
    for (size_t i = 1; i < asc_traj.size(); ++i) {
        result.vertices.push_back(asc_traj[i]);
    }

    if (result.vertices.empty()) return result;

    result.start_vertex = result.vertices.front();
    result.end_vertex = result.vertices.back();
    result.starts_at_lmin = (check_nbr_extremum_type(result.start_vertex, y) == -1);
    result.ends_at_lmax = (check_nbr_extremum_type(result.end_vertex, y) == 1);
    result.total_change = y[result.end_vertex] - y[result.start_vertex];

    return result;
}

// ============================================================================
// Main Computation Function
// ============================================================================

gfc_flow_result_t compute_gfc_flow(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density,
    bool verbose
) {
    gfc_flow_result_t result;
    result.n_vertices = graph.num_vertices();
    result.params = params;

    const size_t n = result.n_vertices;

    if (n == 0 || y.size() != n) {
        return result;
    }

    // Compute median
    std::vector<double> y_sorted = y;
    std::sort(y_sorted.begin(), y_sorted.end());
    result.y_median = y_sorted[n / 2];

    if (verbose) {
        Rprintf("GFC Flow computation: %zu vertices, median = %.4f\n", n, result.y_median);
    }

    // ========================================================================
    // Step 1: Find local extrema and prepare trajectory computation
    // ========================================================================

    if (verbose) {
        Rprintf("Step 1: Finding local extrema...\n");
    }

    auto [lmin_vertices, lmax_vertices] = graph.find_nbr_extrema(y);

    std::unordered_set<size_t> lmin_set(lmin_vertices.begin(), lmin_vertices.end());
    std::unordered_set<size_t> lmax_set(lmax_vertices.begin(), lmax_vertices.end());

    double edge_length_thld = graph.compute_quantile_edge_length(params.edge_length_quantile_thld);
    result.edge_length_thld = edge_length_thld;

    edge_weight_map_t edge_length_weights;
    bool needs_edge_weights = (params.modulation == gflow_modulation_t::EDGELEN ||
                               params.modulation == gflow_modulation_t::DENSITY_EDGELEN);
    if (needs_edge_weights) {
        edge_length_weights = graph.compute_edge_length_weights(-1.0);
    }

    if (verbose) {
        Rprintf("  Found %zu local minima, %zu local maxima\n",
                lmin_vertices.size(), lmax_vertices.size());
    }

    // ========================================================================
    // Precompute ascent map next_up (single-step pointers)
    // ========================================================================
    result.next_up.assign(n, -1);
    for (size_t v = 0; v < n; ++v) {
        result.next_up[v] = select_next_vertex(
            graph, v, y, true,
            params.modulation, density,
            edge_length_weights, edge_length_thld
            );
    }

    // ========================================================================
    // Step 2: Trace trajectories and build basin assignments
    // ========================================================================

    if (verbose) {
        Rprintf("Step 2: Tracing lmin trajectories...\n");
    }

    // Multi-valued membership tracking during trajectory tracing
    std::vector<std::set<size_t>> vertex_to_lmins(n);
    std::vector<std::set<size_t>> vertex_to_lmaxs(n);
    std::vector<bool> covered(n, false);
    std::vector<std::vector<int>> vertex_to_trajs(n);

    int trajectory_id = 0;
    result.n_lmin_trajectories = 0;

    // Sort minima by value for deterministic processing
    std::vector<size_t> lmin_sorted = lmin_vertices;
    std::sort(lmin_sorted.begin(), lmin_sorted.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });

    for (size_t lmin : lmin_sorted) {
        std::vector<size_t> traj = graph.follow_gradient_trajectory(
            lmin, y, true,
            params.modulation, density, edge_length_weights,
            edge_length_thld, params.max_trajectory_length
            );

        if (traj.empty()) continue;

        size_t lmax = traj.back();

        // Check if trajectory ends at a true local maximum
        bool ends_at_true_lmax = (lmax_set.count(lmax) > 0);

        // Only contribute to basin assignments if trajectory reaches true extremum
        if (ends_at_true_lmax) {
            for (size_t v : traj) {
                vertex_to_lmins[v].insert(lmin);
                vertex_to_lmaxs[v].insert(lmax);
                vertex_to_trajs[v].push_back(trajectory_id);
                covered[v] = true;
            }
        }

        // Always store trajectory if requested (for diagnostic purposes)
        if (params.store_trajectories) {
            gflow_trajectory_t gft;
            gft.vertices = std::move(traj);
            gft.start_vertex = lmin;
            gft.end_vertex = lmax;
            gft.starts_at_lmin = true;
            gft.ends_at_lmax = ends_at_true_lmax;
            gft.total_change = y[lmax] - y[lmin];
            gft.trajectory_id = trajectory_id;
            result.trajectories.push_back(std::move(gft));
        }

        ++trajectory_id;
        ++result.n_lmin_trajectories;
    }

    // ========================================================================
    // Step 3: Handle uncovered vertices
    // ========================================================================

    if (verbose) {
        Rprintf("Step 3: Processing uncovered vertices...\n");
    }

    result.n_join_trajectories = 0;

    std::vector<size_t> uncovered;
    for (size_t v = 0; v < n; ++v) {
        if (!covered[v]) uncovered.push_back(v);
    }

    std::sort(uncovered.begin(), uncovered.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });

    for (size_t v : uncovered) {
        if (covered[v]) continue;

        gflow_trajectory_t joined = graph.join_trajectories_at_vertex(
            v, y, params.modulation, density, edge_length_weights,
            edge_length_thld, params.max_trajectory_length
            );

        if (joined.vertices.empty()) {
            covered[v] = true;
            continue;
        }

        size_t lmin = joined.start_vertex;
        size_t lmax = joined.end_vertex;

        // Check if trajectory endpoints are valid
        bool starts_at_true_lmin = (lmin_set.count(lmin) > 0);
        bool ends_at_true_lmax = (lmax_set.count(lmax) > 0);

        // Only contribute to basin assignments if both endpoints are valid
        if (starts_at_true_lmin && ends_at_true_lmax) {
            for (size_t u : joined.vertices) {
                vertex_to_lmins[u].insert(lmin);
                vertex_to_lmaxs[u].insert(lmax);
                vertex_to_trajs[u].push_back(trajectory_id);
                covered[u] = true;
            }
        }

        // Always store trajectory if requested (for diagnostic purposes)
        if (params.store_trajectories) {
            joined.trajectory_id = trajectory_id;
            joined.starts_at_lmin = starts_at_true_lmin;
            // ends_at_lmax already set in join_trajectories_at_vertex
            result.trajectories.push_back(std::move(joined));
        }

        ++trajectory_id;
        ++result.n_join_trajectories;
    }

    if (verbose) {
        size_t n_covered = std::count(covered.begin(), covered.end(), true);
        size_t n_uncovered_final = n - n_covered;
        Rprintf("  Total trajectories: %zu, coverage: %zu/%zu (%.1f%%)\n",
                result.trajectories.size(), n_covered, n,
                100.0 * n_covered / n);
        if (n_uncovered_final > 0) {
            Rprintf("  Note: %zu vertices uncovered (isolated by edge length threshold)\n",
                    n_uncovered_final);
        }
    }

    // ========================================================================
    // Step 4: Build ALL basins from trajectory assignments
    // ========================================================================

    if (verbose) {
        Rprintf("Step 4: Building basin structures...\n");
    }

    // Accumulate vertices for each extremum
    std::map<size_t, std::set<size_t>> ascending_basins;
    std::map<size_t, std::set<size_t>> descending_basins;

    for (size_t v = 0; v < n; ++v) {
        for (size_t lmin : vertex_to_lmins[v]) {
            ascending_basins[lmin].insert(v);
        }
        for (size_t lmax : vertex_to_lmaxs[v]) {
            descending_basins[lmax].insert(v);
        }
    }

    // Build extended basin structures for minima
    for (const auto& [lmin, vertices] : ascending_basins) {
        basin_extended_t basin;
        basin.extremum_vertex = lmin;
        basin.extremum_value = y[lmin];
        basin.is_maximum = false;
        basin.vertices.assign(vertices.begin(), vertices.end());
        
        // Compute hop distances
        basin.hop_distances.resize(basin.vertices.size(), 0);
        std::unordered_map<size_t, int> hop_map;
        hop_map[lmin] = 0;
        
        std::queue<size_t> bfs_queue;
        bfs_queue.push(lmin);
        int max_hop = 0;
        
        while (!bfs_queue.empty()) {
            size_t u = bfs_queue.front();
            bfs_queue.pop();
            int u_hop = hop_map[u];
            
            for (const auto& edge : graph.adjacency_list[u]) {
                size_t nbr = edge.vertex;
                if (vertices.count(nbr) && hop_map.find(nbr) == hop_map.end()) {
                    hop_map[nbr] = u_hop + 1;
                    max_hop = std::max(max_hop, u_hop + 1);
                    bfs_queue.push(nbr);
                }
            }
        }
        
        for (size_t i = 0; i < basin.vertices.size(); ++i) {
            auto it = hop_map.find(basin.vertices[i]);
            basin.hop_distances[i] = (it != hop_map.end()) ? it->second : -1;
        }
        basin.max_hop_distance = max_hop;
        
        result.min_basins_all.push_back(std::move(basin));
    }

    // Build extended basin structures for maxima
    for (const auto& [lmax, vertices] : descending_basins) {
        basin_extended_t basin;
        basin.extremum_vertex = lmax;
        basin.extremum_value = y[lmax];
        basin.is_maximum = true;
        basin.vertices.assign(vertices.begin(), vertices.end());
        
        basin.hop_distances.resize(basin.vertices.size(), 0);
        std::unordered_map<size_t, int> hop_map;
        hop_map[lmax] = 0;
        
        std::queue<size_t> bfs_queue;
        bfs_queue.push(lmax);
        int max_hop = 0;
        
        while (!bfs_queue.empty()) {
            size_t u = bfs_queue.front();
            bfs_queue.pop();
            int u_hop = hop_map[u];
            
            for (const auto& edge : graph.adjacency_list[u]) {
                size_t nbr = edge.vertex;
                if (vertices.count(nbr) && hop_map.find(nbr) == hop_map.end()) {
                    hop_map[nbr] = u_hop + 1;
                    max_hop = std::max(max_hop, u_hop + 1);
                    bfs_queue.push(nbr);
                }
            }
        }
        
        for (size_t i = 0; i < basin.vertices.size(); ++i) {
            auto it = hop_map.find(basin.vertices[i]);
            basin.hop_distances[i] = (it != hop_map.end()) ? it->second : -1;
        }
        basin.max_hop_distance = max_hop;
        
        result.max_basins_all.push_back(std::move(basin));
    }

    // Compute summaries for all basins
    auto min_basins_compact = std::vector<basin_compact_t>(
        result.min_basins_all.begin(), result.min_basins_all.end());
    auto max_basins_compact = std::vector<basin_compact_t>(
        result.max_basins_all.begin(), result.max_basins_all.end());
    
    auto min_summaries_base = gfc_internal::compute_basin_summaries(
        graph, min_basins_compact, y, result.y_median, params.hop_k);
    auto max_summaries_base = gfc_internal::compute_basin_summaries(
        graph, max_basins_compact, y, result.y_median, params.hop_k);

    // Convert to extended summaries
    result.min_summaries_all.reserve(min_summaries_base.size());
    for (const auto& s : min_summaries_base) {
        result.min_summaries_all.emplace_back(s);
    }
    result.max_summaries_all.reserve(max_summaries_base.size());
    for (const auto& s : max_summaries_base) {
        result.max_summaries_all.emplace_back(s);
    }

    // Record initial counts
    result.stage_history.emplace_back(
        "initial",
        static_cast<int>(result.max_basins_all.size()),
        static_cast<int>(result.max_basins_all.size()),
        static_cast<int>(result.min_basins_all.size()),
        static_cast<int>(result.min_basins_all.size())
    );

    if (verbose) {
        Rprintf("  Built %zu min basins, %zu max basins\n",
                result.min_basins_all.size(), result.max_basins_all.size());
    }

    // ========================================================================
    // Step 5: Apply filtering (marking spurious, not removing)
    // ========================================================================

    // Count non-spurious before each stage
    auto count_retained = [](const auto& basins) {
        return std::count_if(basins.begin(), basins.end(),
                             [](const auto& b) { return !b.is_spurious; });
    };

    // 5a: Relative value filter
    if (params.apply_relvalue_filter) {
        if (verbose) {
            Rprintf("Step 5a: Relative value filter...\n");
        }
        
        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);
        
        gfc_flow_internal::mark_spurious_by_relvalue(
            result.max_summaries_all, result.max_basins_all,
            params.min_rel_value_max, params.max_rel_value_min, true);
        gfc_flow_internal::mark_spurious_by_relvalue(
            result.min_summaries_all, result.min_basins_all,
            params.min_rel_value_max, params.max_rel_value_min, false);
        
        int n_max_after = count_retained(result.max_basins_all);
        int n_min_after = count_retained(result.min_basins_all);
        
        result.stage_history.emplace_back("relvalue",
            n_max_before, n_max_after, n_min_before, n_min_after);
        
        if (verbose) {
            Rprintf("  Maxima: %d -> %d, Minima: %d -> %d\n",
                    n_max_before, n_max_after, n_min_before, n_min_after);
        }
    }

    // 5b: Cluster maxima
    if (params.apply_maxima_clustering) {
        if (verbose) {
            Rprintf("Step 5b: Clustering maxima...\n");
        }

        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);

        result.max_overlap_dist = gfc_flow_internal::cluster_and_merge_basins(
            result.max_summaries_all, result.max_basins_all,
            params.max_overlap_threshold, true);

        int n_max_after = count_retained(result.max_basins_all);

        result.stage_history.emplace_back("merge.max",
            n_max_before, n_max_after, n_min_before, n_min_before);

        if (verbose) {
            Rprintf("  Maxima: %d -> %d\n", n_max_before, n_max_after);
        }
    } else {

        result.max_overlap_dist = gfc_flow_internal::get_overlap_distance_matrix(result.max_basins_all);
    }

    // 5c: Cluster minima
    if (params.apply_minima_clustering) {
        if (verbose) {
            Rprintf("Step 5c: Clustering minima...\n");
        }

        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);

        result.min_overlap_dist = gfc_flow_internal::cluster_and_merge_basins(
            result.min_summaries_all, result.min_basins_all,
            params.min_overlap_threshold, false);

        int n_min_after = count_retained(result.min_basins_all);

        result.stage_history.emplace_back("merge.min",
            n_max_before, n_max_before, n_min_before, n_min_after);

        if (verbose) {
            Rprintf("  Minima: %d -> %d\n", n_min_before, n_min_after);
        }
    } else {

        result.min_overlap_dist = gfc_flow_internal::get_overlap_distance_matrix(result.min_basins_all);
    }

    // 5d: Minimum basin size filter (independent of geometric filter flag)
    {
        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);

        gfc_flow_internal::mark_spurious_by_min_basin_size(
            result.max_summaries_all, result.max_basins_all, params.min_basin_size);
        gfc_flow_internal::mark_spurious_by_min_basin_size(
            result.min_summaries_all, result.min_basins_all, params.min_basin_size);

        int n_max_after = count_retained(result.max_basins_all);
        int n_min_after = count_retained(result.min_basins_all);

        result.stage_history.emplace_back("min_basin_size",
                                          n_max_before, n_max_after, n_min_before, n_min_after);

        if (verbose) {
            Rprintf("Step 5d: Min basin size filter... Maxima: %d -> %d, Minima: %d -> %d\n",
                    n_max_before, n_max_after, n_min_before, n_min_after);
        }
    }

    // 5e: Minimum number of trajectories filter
    if (params.min_n_trajectories > 0) {

        if (!params.store_trajectories) {
            Rf_error("min_n_trajectories > 0 requires store_trajectories = TRUE.");
        }

        std::unordered_map<size_t, int> max_traj_counts;
        std::unordered_map<size_t, int> min_traj_counts;

        for (const auto& tr : result.trajectories) {
            if (tr.ends_at_lmax) {
                max_traj_counts[tr.end_vertex] += 1;
            }
            if (tr.starts_at_lmin) {
                min_traj_counts[tr.start_vertex] += 1;
            }
        }

        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);

        gfc_flow_internal::mark_spurious_by_min_n_trajectories(
            result.max_summaries_all, result.max_basins_all,
            max_traj_counts, params.min_n_trajectories);
        gfc_flow_internal::mark_spurious_by_min_n_trajectories(
            result.min_summaries_all, result.min_basins_all,
            min_traj_counts, params.min_n_trajectories);

        int n_max_after = count_retained(result.max_basins_all);
        int n_min_after = count_retained(result.min_basins_all);

        result.stage_history.emplace_back("min_n_trajectories",
                                          n_max_before, n_max_after, n_min_before, n_min_after);

        if (verbose) {
            Rprintf("Step 5e: Min n trajectories filter... Maxima: %d -> %d, Minima: %d -> %d\n",
                    n_max_before, n_max_after, n_min_before, n_min_after);
        }
    }

    // 5f: Geometric filter
    if (params.apply_geometric_filter) {
        if (verbose) {
            Rprintf("Step 5d: Geometric filter...\n");
        }
        
        int n_max_before = count_retained(result.max_basins_all);
        int n_min_before = count_retained(result.min_basins_all);
        
        gfc_flow_internal::mark_spurious_by_geometry(
            result.max_summaries_all, result.max_basins_all,
            params.p_mean_nbrs_dist_threshold,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size);
        gfc_flow_internal::mark_spurious_by_geometry(
            result.min_summaries_all, result.min_basins_all,
            params.p_mean_nbrs_dist_threshold,
            params.p_mean_hopk_dist_threshold,
            params.p_deg_threshold,
            params.min_basin_size);
        
        int n_max_after = count_retained(result.max_basins_all);
        int n_min_after = count_retained(result.min_basins_all);
        
        result.stage_history.emplace_back("geometric",
            n_max_before, n_max_after, n_min_before, n_min_after);
        
        if (verbose) {
            Rprintf("  Maxima: %d -> %d, Minima: %d -> %d\n",
                    n_max_before, n_max_after, n_min_before, n_min_after);
        }
    }

    // ========================================================================
    // Step 6: Assign labels and build index vectors
    // ========================================================================

    if (verbose) {
        Rprintf("Step 6: Assigning labels...\n");
    }

    gfc_flow_internal::assign_labels(
        result.max_summaries_all, result.max_basins_all, true);
    gfc_flow_internal::assign_labels(
        result.min_summaries_all, result.min_basins_all, false);

    // Build retained/spurious index vectors
    for (int i = 0; i < static_cast<int>(result.max_basins_all.size()); ++i) {
        if (result.max_basins_all[i].is_spurious) {
            result.spurious_max_indices.push_back(i);
        } else {
            result.retained_max_indices.push_back(i);
        }
    }
    for (int i = 0; i < static_cast<int>(result.min_basins_all.size()); ++i) {
        if (result.min_basins_all[i].is_spurious) {
            result.spurious_min_indices.push_back(i);
        } else {
            result.retained_min_indices.push_back(i);
        }
    }

    result.n_max_retained = static_cast<int>(result.retained_max_indices.size());
    result.n_min_retained = static_cast<int>(result.retained_min_indices.size());
    result.n_max_spurious = static_cast<int>(result.spurious_max_indices.size());
    result.n_min_spurious = static_cast<int>(result.spurious_min_indices.size());

    // ========================================================================
    // Step 7: Update trajectory endpoint classification
    // ========================================================================

    // Build vertex-to-basin-index maps
    std::unordered_map<size_t, int> lmin_to_basin_idx;
    std::unordered_map<size_t, int> lmax_to_basin_idx;
    
    for (int i = 0; i < static_cast<int>(result.min_basins_all.size()); ++i) {
        lmin_to_basin_idx[result.min_basins_all[i].extremum_vertex] = i;
    }
    for (int i = 0; i < static_cast<int>(result.max_basins_all.size()); ++i) {
        lmax_to_basin_idx[result.max_basins_all[i].extremum_vertex] = i;
    }

    for (auto& traj : result.trajectories) {
        // Start endpoint
        auto it_min = lmin_to_basin_idx.find(traj.start_vertex);
        if (it_min != lmin_to_basin_idx.end()) {
            traj.start_basin_idx = it_min->second;
            traj.start_is_spurious = result.min_basins_all[it_min->second].is_spurious;
        }
        
        // End endpoint
        auto it_max = lmax_to_basin_idx.find(traj.end_vertex);
        if (it_max != lmax_to_basin_idx.end()) {
            traj.end_basin_idx = it_max->second;
            traj.end_is_spurious = result.max_basins_all[it_max->second].is_spurious;
        }
    }

    // ========================================================================
    // Step 8: Build membership vectors
    // ========================================================================

    if (verbose) {
        Rprintf("Step 8: Building membership vectors...\n");
    }

    result.max_membership_all = gfc_flow_internal::build_membership_all(
        result.max_basins_all, n);
    result.min_membership_all = gfc_flow_internal::build_membership_all(
        result.min_basins_all, n);

    result.max_membership_retained = gfc_flow_internal::build_membership_retained(
        result.max_basins_all, result.retained_max_indices, n);
    result.min_membership_retained = gfc_flow_internal::build_membership_retained(
        result.min_basins_all, result.retained_min_indices, n);

    // Single-valued assignments (first retained basin)
    result.max_assignment.resize(n, -1);
    result.min_assignment.resize(n, -1);

    for (size_t v = 0; v < n; ++v) {
        if (!result.max_membership_retained[v].empty()) {
            result.max_assignment[v] = result.max_membership_retained[v][0];
        }
        if (!result.min_membership_retained[v].empty()) {
            result.min_assignment[v] = result.min_membership_retained[v][0];
        }
    }

    if (verbose) {
        Rprintf("GFC Flow complete: %d retained maxima, %d retained minima\n",
                result.n_max_retained, result.n_min_retained);
        Rprintf("                   %d spurious maxima, %d spurious minima\n",
                result.n_max_spurious, result.n_min_spurious);
    }

    return result;
}

std::vector<gfc_flow_result_t> compute_gfc_flow_matrix(
    const set_wgraph_t& graph,
    const Eigen::MatrixXd& Y,
    const gfc_flow_params_t& params,
    const std::vector<double>& density,
    int n_cores,
    bool verbose
) {
    const int n = static_cast<int>(Y.rows());
    const int p = static_cast<int>(Y.cols());

    if (verbose) {
        Rprintf("GFC Flow matrix: %d functions, %d vertices\n", p, n);
    }

    std::vector<gfc_flow_result_t> results(p);

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
    for (int j = 0; j < p; ++j) {
        std::vector<double> y_j(n);
        for (int i = 0; i < n; ++i) {
            y_j[i] = Y(i, j);
        }
        results[j] = compute_gfc_flow(graph, y_j, params, density, false);
    }

    return results;
}
