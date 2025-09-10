#include <vector>        // For std::vector
#include <algorithm>     // For std::sort
#include <utility>       // For std::pair
#include <cmath>         // For mathematical operations
#include <limits>        // For std::numeric_limits
#include <queue>         // For std::priority_queue
#include <functional>    // For std::greater
#include <mutex>         // For std::mutex
#include <execution>     // For parallel execution

#include "exec_policy.hpp"
#include "edge_pruning_stats.hpp" // For edge_pruning_stats_t
#include "set_wgraph.hpp"     // For set_wgraph_t definition
#include "error_utils.h"      // For REPORT_ERROR

/**
 * @brief Compute the length of the shortest path between two vertices while excluding the direct edge between them
 *
 * This function implements a bidirectional Dijkstra algorithm that finds the shortest path
 * between the vertices of an edge, while ignoring the edge connecting them.
 *
 * @param source The starting vertex index
 * @param target The destination vertex index
 * @return The length of the shortest path from source to target excluding the direct edge.
 *         Returns infinity if no alternative path exists.
 */
double set_wgraph_t::bidirectional_dijkstra_excluding_edge(
    size_t source,
    size_t target) const
{
    const size_t n = adjacency_list.size();

    // Initialize distance arrays for forward and backward search
    std::vector<double> dist_forward(n, std::numeric_limits<double>::infinity());
    std::vector<double> dist_backward(n, std::numeric_limits<double>::infinity());

    dist_forward[source] = 0.0;
    dist_backward[target] = 0.0;

    // Initialize priority queues for forward and backward search
    std::priority_queue<
        std::pair<double, size_t>,
        std::vector<std::pair<double, size_t>>,
        std::greater<std::pair<double, size_t>>
    > pq_forward, pq_backward;

    pq_forward.push({0.0, source});
    pq_backward.push({0.0, target});

    // Track visited vertices in both directions
    std::vector<bool> visited_forward(n, false);
    std::vector<bool> visited_backward(n, false);

    // Best path length found so far
    double best_path_length = std::numeric_limits<double>::infinity();

    // Main bidirectional search loop
    while (!pq_forward.empty() && !pq_backward.empty()) {
        // Check termination condition
        if (best_path_length <= pq_forward.top().first + pq_backward.top().first) {
            break;  // No better path possible
        }

        // Process forward direction
        {
            auto [dist, vertex] = pq_forward.top();
            pq_forward.pop();

            if (visited_forward[vertex]) continue;
            visited_forward[vertex] = true;

            // Check if this vertex has been reached from both directions
            if (dist_backward[vertex] != std::numeric_limits<double>::infinity()) {
                double path_length = dist_forward[vertex] + dist_backward[vertex];
                if (path_length < best_path_length) {
                    best_path_length = path_length;
                }
            }

            // Relax outgoing edges, excluding the direct edge to target
            for (const auto& edge : adjacency_list[vertex]) {
                size_t neighbor = edge.vertex;

                // Skip the direct edge between source and target
                if (vertex == source && neighbor == target) {
                    continue;
                }

                double weight = edge.weight;

                if (dist_forward[vertex] + weight < dist_forward[neighbor]) {
                    dist_forward[neighbor] = dist_forward[vertex] + weight;
                    pq_forward.push({dist_forward[neighbor], neighbor});
                }
            }
        }

        // Process backward direction
        {
            auto [dist, vertex] = pq_backward.top();
            pq_backward.pop();

            if (visited_backward[vertex]) continue;
            visited_backward[vertex] = true;

            // Check if this vertex has been reached from both directions
            if (dist_forward[vertex] != std::numeric_limits<double>::infinity()) {
                double path_length = dist_forward[vertex] + dist_backward[vertex];
                if (path_length < best_path_length) {
                    best_path_length = path_length;
                }
            }

            // Relax incoming edges (backward direction), excluding the direct edge from source
            for (const auto& edge : adjacency_list[vertex]) {
                size_t neighbor = edge.vertex;

                // Skip the direct edge between target and source
                if (vertex == target && neighbor == source) {
                    continue;
                }

                double weight = edge.weight;

                if (dist_backward[vertex] + weight < dist_backward[neighbor]) {
                    dist_backward[neighbor] = dist_backward[vertex] + weight;
                    pq_backward.push({dist_backward[neighbor], neighbor});
                }
            }
        }
    }

    return best_path_length;
}


/**
 * @brief Compute statistics for potential edge pruning based on geometric criteria
 *
 * This function analyzes edges longer than a specified threshold (default: median edge length)
 * to determine how pruning each edge would affect the graph's geometry. For each candidate edge,
 * it computes the ratio of the alternative geodesic path length to the original edge length.
 *
 * @param threshold_percentile Percentile threshold for edge Rf_length(0.0-1.0).
 *        Default 0.5 (median). Only edges with length greater than or equal to
 *        the value at this percentile are analyzed.
 * @return edge_pruning_stats_t structure containing statistics for each analyzed edge
 */
edge_pruning_stats_t set_wgraph_t::compute_edge_pruning_stats(double threshold_percentile) const {
    // Structure to represent an edge with source, target, and weight
    struct Edge {
        size_t source;
        size_t target;
        double weight;

        Edge(size_t s, size_t t, double w) : source(s), target(t), weight(w) {}

        bool operator<(const Edge& other) const {
            return weight < other.weight;  // Ascending order for finding median
        }
    };

    // 1. Identify all edges in the graph
    std::vector<Edge> all_edges;
    all_edges.reserve(adjacency_list.size() * 5);  // Rough estimate of edges

    // Using for_each with parallel execution policy
    std::vector<size_t> vertices(adjacency_list.size());
    std::iota(vertices.begin(), vertices.end(), 0);

    std::mutex edges_mutex;  // Mutex for thread-safe access to all_edges

    //std::for_each(std::execution::par_unseq, vertices.begin(), vertices.end(),
    gflow::for_each(GFLOW_EXEC_POLICY, vertices.begin(), vertices.end(),
        [this, &all_edges, &edges_mutex](size_t vertex) {
            std::vector<Edge> local_edges;

            for (const auto& edge_info : adjacency_list[vertex]) {
                size_t target = edge_info.vertex;

                // Only add each edge once (when source < target)
                if (vertex < target) {
                    local_edges.emplace_back(vertex, target, edge_info.weight);
                }
            }

            // Thread-safe merging of local edges into all_edges
            if (!local_edges.empty()) {
                std::lock_guard<std::mutex> lock(edges_mutex);
                all_edges.insert(all_edges.end(), local_edges.begin(), local_edges.end());
            }
        }
    );

    // Return early if no edges
    if (all_edges.empty()) {
        edge_pruning_stats_t empty_stats;
        empty_stats.median_edge_length = 0.0;
        return empty_stats;
    }

    // 2. Sort edges by weight in ascending order for finding median
    std::sort(all_edges.begin(), all_edges.end());

    // 3. Compute the median edge weight (or other threshold)
    double threshold_weight;
    if (threshold_percentile <= 0.0) {
        threshold_weight = all_edges.front().weight - 1.0;  // Analyze all edges
    } else if (threshold_percentile >= 1.0) {
        threshold_weight = all_edges.back().weight + 1.0;  // Analyze no edges
    } else {
        size_t threshold_index = static_cast<size_t>(all_edges.size() * threshold_percentile);
        threshold_index = std::min(threshold_index, all_edges.size() - 1);
        threshold_weight = all_edges[threshold_index].weight;
    }

    // Initialize result structure
    edge_pruning_stats_t result;
    result.median_edge_length = threshold_weight;

    // 4. Analyze edges longer than or equal to the threshold
    for (const auto& edge : all_edges) {
        if (edge.weight >= threshold_weight) {
            // Find the length of alternative geodesic path (excluding this edge)
            double alt_path_length = bidirectional_dijkstra_excluding_edge(
                edge.source, edge.target);
            
            // Only include edges that have an alternative path
            if (alt_path_length != std::numeric_limits<double>::infinity()) {
                result.stats.emplace_back(
                    edge.source, edge.target, edge.weight, alt_path_length);
            }
        }
    }

    return result;
}

/**
 * @brief Prune edges geometrically based on alternative path ratio
 *
 * This function prunes edges whose alternative geodesic path length is within a specified ratio
 * of the original edge length. This approach preserves the graph's geometric structure while
 * removing redundant edges.
 *
 * @param max_ratio_threshold Maximum acceptable ratio of alternative path length to edge length
 *        (default: 1.2). Edges with ratio <= this value will be pruned.
 * @param threshold_percentile Percentile threshold for edge Rf_length(0.0-1.0).
 *        Default 0.5 (median). Only edges with length greater than or equal to
 *        the value at this percentile are considered for pruning.
 * @return A new set_wgraph_t with redundant edges removed
 */
set_wgraph_t set_wgraph_t::prune_edges_geometrically(
    double max_ratio_threshold,
    double threshold_percentile) const
{
    // Compute edge pruning statistics
    edge_pruning_stats_t stats = compute_edge_pruning_stats(threshold_percentile);
    
    // Create a copy of the graph to prune
    set_wgraph_t pruned_graph(*this);
    
    // Get edges that can be pruned based on the ratio threshold
    std::vector<std::pair<size_t, size_t>> prunable_edges = 
        stats.get_prunable_edges(max_ratio_threshold);
    
    // Remove all prunable edges
    for (const auto& [source, target] : prunable_edges) {
        pruned_graph.remove_edge(source, target);
    }
    
    return pruned_graph;
}


/**
 * @brief Get the list of edges that can be pruned based on a length ratio threshold
 *
 * This method filters the edge statistics to find edges whose alternative path to edge length
 * ratio is below or equal to the specified threshold.
 *
 * @param max_ratio_threshold Maximum acceptable ratio of alternative path length to edge length
 * @return Vector of (source, target) pairs representing prunable edges
 */
std::vector<std::pair<size_t, size_t>> edge_pruning_stats_t::get_prunable_edges(
    double max_ratio_threshold
    ) const {

    std::vector<std::pair<size_t, size_t>> prunable_edges;
    for (const auto& stat : stats) {
        if (stat.length_ratio <= max_ratio_threshold) {
            prunable_edges.emplace_back(stat.source, stat.target);
        }
    }
    return prunable_edges;
}
