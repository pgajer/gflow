/**
 * @file madag.cpp
 * @brief Implementation of Monotonic Ascent DAG construction and analysis
 *
 * This file implements the algorithms for building MADAGs from local minima,
 * computing cell supports, counting and enumerating trajectories, and
 * analyzing trajectory structure.
 */

#include "madag.hpp"
#include "set_wgraph.hpp"

#include <algorithm>
#include <queue>
#include <stack>
#include <numeric>
#include <cmath>
#include <random>
#include <sstream>

#include <R.h>

// ============================================================================
// ms_cell_t Implementation
// ============================================================================

std::string ms_cell_t::id_string() const {
    std::ostringstream oss;
    oss << "(m" << min_vertex << ", M" << max_vertex << ")";
    return oss.str();
}

// ============================================================================
// MADAG Construction
// ============================================================================

madag_t construct_madag(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    size_t source_vertex,
    const madag_params_t& params,
    bool verbose
) {
    madag_t madag;
    const size_t n = graph.num_vertices();
    
    if (source_vertex >= n || y.size() != n) {
        if (verbose) {
            Rprintf("ERROR: Invalid source vertex or y vector size\n");
        }
        return madag;
    }
    
    madag.source_vertex = source_vertex;
    madag.source_value = y[source_vertex];
    
    if (verbose) {
        Rprintf("Constructing MADAG from source vertex %zu (y = %.4f)\n",
                source_vertex, madag.source_value);
    }
    
    // Compute edge length threshold if specified
    double edge_length_thld = std::numeric_limits<double>::infinity();
    if (params.edge_length_quantile_thld < 1.0) {
        edge_length_thld = graph.compute_quantile_edge_length(params.edge_length_quantile_thld);
        if (verbose) {
            Rprintf("  Edge length threshold (%.2f quantile): %.4f\n",
                    params.edge_length_quantile_thld, edge_length_thld);
        }
    }
    
    // Phase 1: BFS-like exploration to find all reachable vertices
    // Unlike standard BFS, we allow a vertex to be "discovered" multiple times
    // through different predecessors, since we need to track all predecessors.
    
    // Queue contains vertices whose successors need to be explored
    std::queue<size_t> to_explore;
    to_explore.push(source_vertex);
    
    // Track which vertices have been added to the queue at least once
    std::unordered_set<size_t> queued;
    queued.insert(source_vertex);
    
    // Initialize source vertex
    madag.reachable_set.insert(source_vertex);
    madag.predecessors[source_vertex] = {};  // Source has no predecessors
    madag.successors[source_vertex] = {};    // Will be populated below
    
    while (!to_explore.empty()) {
        size_t u = to_explore.front();
        to_explore.pop();
        
        double y_u = y[u];
        
        // Explore all neighbors of u
        for (const auto& edge : graph.adjacency_list[u]) {
            size_t v = edge.vertex;
            double edge_len = edge.weight;
            
            // Skip if edge is too long
            if (edge_len > edge_length_thld) {
                continue;
            }
            
            // Check monotonicity: only follow ascending edges
            if (y[v] > y_u) {
                // v is a monotonic successor of u
                madag.successors[u].push_back(v);
                
                // u is a monotonic predecessor of v
                madag.predecessors[v].push_back(u);
                
                // Mark v as reachable
                madag.reachable_set.insert(v);
                
                // Add v to exploration queue if not already queued
                if (queued.find(v) == queued.end()) {
                    to_explore.push(v);
                    queued.insert(v);
                    
                    // Initialize successor list for v
                    madag.successors[v] = {};
                }
            }
        }
    }
    
    // Convert reachable set to vector
    madag.reachable_vertices.assign(madag.reachable_set.begin(), madag.reachable_set.end());
    
    // Sort reachable vertices by y value for consistent ordering
    std::sort(madag.reachable_vertices.begin(), madag.reachable_vertices.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });
    
    // Phase 2: Identify sinks (local maxima reachable from source)
    for (size_t v : madag.reachable_vertices) {
        if (madag.is_sink(v)) {
            madag.reachable_maxima.push_back(v);
        }
    }
    
    // Sort maxima by y value
    std::sort(madag.reachable_maxima.begin(), madag.reachable_maxima.end(),
              [&y](size_t a, size_t b) { return y[a] < y[b]; });
    
    if (verbose) {
        Rprintf("  Reachable vertices: %zu\n", madag.reachable_vertices.size());
        Rprintf("  Reachable maxima: %zu\n", madag.reachable_maxima.size());
    }
    
    // Phase 3: Compute topological ordering
    madag.topological_order = compute_topological_order(madag);
    madag.reverse_topological_order = madag.topological_order;
    std::reverse(madag.reverse_topological_order.begin(), 
                 madag.reverse_topological_order.end());
    
    // Phase 4: Compute path counts if requested
    if (params.compute_path_counts) {
        compute_path_counts_from_source(madag);
        compute_path_counts_to_sinks(madag);
        
        if (verbose) {
            Rprintf("  Path counts computed\n");
        }
    }
    
    // Phase 5: Build cell structures
    for (size_t i = 0; i < madag.reachable_maxima.size(); ++i) {
        size_t max_v = madag.reachable_maxima[i];
        
        ms_cell_t cell;
        cell.min_vertex = source_vertex;
        cell.max_vertex = max_v;
        cell.min_value = y[source_vertex];
        cell.max_value = y[max_v];
        
        // Compute support
        cell.support = compute_cell_support(madag, max_v);
        
        // Count trajectories
        if (params.compute_path_counts) {
            cell.n_trajectories = count_cell_trajectories(madag, max_v);
        }
        
        // Check if enumeration is feasible
        if (params.enumerate_trajectories && 
            cell.n_trajectories <= params.max_trajectories_per_cell) {
            cell.trajectories = enumerate_cell_trajectories(madag, y, max_v, 0);
            cell.explicitly_enumerated = true;
        } else {
            cell.explicitly_enumerated = false;
        }
        
        // Store cell
        madag.max_to_cell_idx[max_v] = madag.cells.size();
        madag.cells.push_back(std::move(cell));
    }
    
    if (verbose) {
        Rprintf("  Cells constructed: %zu\n", madag.cells.size());
        for (const auto& cell : madag.cells) {
            Rprintf("    %s: support=%zu, trajectories=%zu%s\n",
                    cell.id_string().c_str(),
                    cell.support.size(),
                    cell.n_trajectories,
                    cell.explicitly_enumerated ? " (enumerated)" : " (not enumerated)");
        }
    }
    
    return madag;
}

// ============================================================================
// Topological Ordering
// ============================================================================

std::vector<size_t> compute_topological_order(const madag_t& madag) {
    std::vector<size_t> order;
    order.reserve(madag.reachable_vertices.size());
    
    // Compute in-degrees
    std::unordered_map<size_t, size_t> in_degree;
    for (size_t v : madag.reachable_vertices) {
        auto it = madag.predecessors.find(v);
        in_degree[v] = (it != madag.predecessors.end()) ? it->second.size() : 0;
    }
    
    // Initialize queue with vertices having in-degree 0 (should be just source)
    std::queue<size_t> ready;
    for (size_t v : madag.reachable_vertices) {
        if (in_degree[v] == 0) {
            ready.push(v);
        }
    }
    
    // Kahn's algorithm
    while (!ready.empty()) {
        size_t u = ready.front();
        ready.pop();
        order.push_back(u);
        
        auto it = madag.successors.find(u);
        if (it != madag.successors.end()) {
            for (size_t v : it->second) {
                --in_degree[v];
                if (in_degree[v] == 0) {
                    ready.push(v);
                }
            }
        }
    }
    
    return order;
}

// ============================================================================
// Path Counting
// ============================================================================

void compute_path_counts_from_source(madag_t& madag) {
    madag.path_count_from_source.clear();
    
    // Source has exactly 1 path to itself (the empty path)
    madag.path_count_from_source[madag.source_vertex] = 1;
    
    // Process vertices in topological order
    for (size_t u : madag.topological_order) {
        size_t count_u = madag.path_count_from_source[u];
        
        auto it = madag.successors.find(u);
        if (it != madag.successors.end()) {
            for (size_t v : it->second) {
                madag.path_count_from_source[v] += count_u;
            }
        }
    }
}

void compute_path_counts_to_sinks(madag_t& madag) {
    madag.path_count_to_sinks.clear();
    
    // Sinks have exactly 1 path to themselves
    for (size_t v : madag.reachable_maxima) {
        madag.path_count_to_sinks[v] = 1;
    }
    
    // Process vertices in reverse topological order
    for (size_t v : madag.reverse_topological_order) {
        // Skip if already set (sinks)
        if (madag.path_count_to_sinks.find(v) != madag.path_count_to_sinks.end()) {
            continue;
        }
        
        size_t count = 0;
        auto it = madag.successors.find(v);
        if (it != madag.successors.end()) {
            for (size_t w : it->second) {
                count += madag.path_count_to_sinks[w];
            }
        }
        madag.path_count_to_sinks[v] = count;
    }
}

// ============================================================================
// Cell Analysis
// ============================================================================

std::vector<size_t> compute_cell_support(
    const madag_t& madag,
    size_t max_vertex
) {
    // We need to find vertices that:
    // 1. Are reachable from source (already guaranteed by being in MADAG)
    // 2. Can reach max_vertex via monotonic paths
    
    // Backward BFS from max_vertex to find vertices that can reach it
    std::unordered_set<size_t> can_reach_max;
    std::queue<size_t> bfs_queue;
    bfs_queue.push(max_vertex);
    can_reach_max.insert(max_vertex);
    
    while (!bfs_queue.empty()) {
        size_t v = bfs_queue.front();
        bfs_queue.pop();
        
        auto it = madag.predecessors.find(v);
        if (it != madag.predecessors.end()) {
            for (size_t u : it->second) {
                if (can_reach_max.find(u) == can_reach_max.end()) {
                    can_reach_max.insert(u);
                    bfs_queue.push(u);
                }
            }
        }
    }
    
    // The support is all vertices that can reach max_vertex,
    // excluding the source and max_vertex themselves
    std::vector<size_t> support;
    for (size_t v : can_reach_max) {
        if (v != madag.source_vertex && v != max_vertex) {
            support.push_back(v);
        }
    }
    
    return support;
}

size_t count_cell_trajectories(
    const madag_t& madag,
    size_t max_vertex
) {
    // The number of trajectories from source to max_vertex is simply
    // the path count from source to max_vertex
    auto it = madag.path_count_from_source.find(max_vertex);
    if (it != madag.path_count_from_source.end()) {
        return it->second;
    }
    return 0;
}

std::vector<size_t> identify_bottlenecks(
    const madag_t& madag,
    size_t max_vertex,
    double min_fraction
) {
    std::vector<size_t> bottlenecks;
    
    // Total trajectories in the cell
    size_t total_trajectories = count_cell_trajectories(madag, max_vertex);
    if (total_trajectories == 0) {
        return bottlenecks;
    }
    
    // For each vertex in the cell support, compute fraction of trajectories
    // passing through it: (paths from source to v) * (paths from v to max) / total
    
    // First, compute paths from each vertex to max_vertex specifically
    // This requires a separate backward DP from max_vertex
    std::unordered_map<size_t, size_t> paths_to_max;
    paths_to_max[max_vertex] = 1;
    
    // Process in reverse topological order, but only for vertices that can reach max
    std::unordered_set<size_t> can_reach_max;
    can_reach_max.insert(max_vertex);
    
    // Backward BFS to find vertices that can reach max
    std::queue<size_t> bfs_queue;
    bfs_queue.push(max_vertex);
    while (!bfs_queue.empty()) {
        size_t v = bfs_queue.front();
        bfs_queue.pop();
        
        auto it = madag.predecessors.find(v);
        if (it != madag.predecessors.end()) {
            for (size_t u : it->second) {
                if (can_reach_max.find(u) == can_reach_max.end()) {
                    can_reach_max.insert(u);
                    bfs_queue.push(u);
                }
            }
        }
    }
    
    // Compute paths to max in reverse topological order
    for (auto it = madag.reverse_topological_order.rbegin(); 
         it != madag.reverse_topological_order.rend(); ++it) {
        size_t v = *it;
        
        if (can_reach_max.find(v) == can_reach_max.end()) {
            continue;
        }
        if (paths_to_max.find(v) != paths_to_max.end()) {
            continue;  // Already computed (max_vertex)
        }
        
        size_t count = 0;
        auto succ_it = madag.successors.find(v);
        if (succ_it != madag.successors.end()) {
            for (size_t w : succ_it->second) {
                if (can_reach_max.find(w) != can_reach_max.end()) {
                    count += paths_to_max[w];
                }
            }
        }
        paths_to_max[v] = count;
    }
    
    // Check each vertex for bottleneck status
    for (size_t v : can_reach_max) {
        if (v == madag.source_vertex || v == max_vertex) {
            continue;  // Skip endpoints
        }
        
        auto from_source_it = madag.path_count_from_source.find(v);
        auto to_max_it = paths_to_max.find(v);
        
        if (from_source_it != madag.path_count_from_source.end() &&
            to_max_it != paths_to_max.end()) {
            size_t paths_through = from_source_it->second * to_max_it->second;
            double fraction = static_cast<double>(paths_through) / total_trajectories;
            
            if (fraction >= min_fraction) {
                bottlenecks.push_back(v);
            }
        }
    }
    
    return bottlenecks;
}

// ============================================================================
// Trajectory Enumeration
// ============================================================================

std::vector<madag_trajectory_t> enumerate_cell_trajectories(
    const madag_t& madag,
    const std::vector<double>& y,
    size_t max_vertex,
    size_t max_trajectories
) {
    std::vector<madag_trajectory_t> trajectories;
    
    // First, identify which vertices can reach max_vertex
    std::unordered_set<size_t> can_reach_max;
    std::queue<size_t> bfs_queue;
    bfs_queue.push(max_vertex);
    can_reach_max.insert(max_vertex);
    
    while (!bfs_queue.empty()) {
        size_t v = bfs_queue.front();
        bfs_queue.pop();
        
        auto it = madag.predecessors.find(v);
        if (it != madag.predecessors.end()) {
            for (size_t u : it->second) {
                if (can_reach_max.find(u) == can_reach_max.end()) {
                    can_reach_max.insert(u);
                    bfs_queue.push(u);
                }
            }
        }
    }
    
    // DFS enumeration of paths
    // Stack contains: (current_vertex, current_path)
    std::stack<std::pair<size_t, std::vector<size_t>>> dfs_stack;
    dfs_stack.push({madag.source_vertex, {madag.source_vertex}});
    
    while (!dfs_stack.empty()) {
        auto [v, path] = dfs_stack.top();
        dfs_stack.pop();
        
        if (v == max_vertex) {
            // Complete trajectory found
            madag_trajectory_t traj;
            traj.vertices = std::move(path);
            traj.source_min = madag.source_vertex;
            traj.sink_max = max_vertex;
            traj.total_ascent = y[max_vertex] - y[madag.source_vertex];
            traj.cluster_id = -1;
            
            trajectories.push_back(std::move(traj));
            
            // Check limit
            if (max_trajectories > 0 && trajectories.size() >= max_trajectories) {
                return trajectories;
            }
            continue;
        }
        
        // Explore successors that can reach max_vertex
        auto it = madag.successors.find(v);
        if (it != madag.successors.end()) {
            for (size_t w : it->second) {
                if (can_reach_max.find(w) != can_reach_max.end()) {
                    std::vector<size_t> new_path = path;
                    new_path.push_back(w);
                    dfs_stack.push({w, std::move(new_path)});
                }
            }
        }
    }
    
    return trajectories;
}

std::vector<madag_trajectory_t> sample_cell_trajectories(
    const madag_t& madag,
    const std::vector<double>& y,
    size_t max_vertex,
    size_t n_samples,
    unsigned int seed
) {
    std::vector<madag_trajectory_t> trajectories;
    trajectories.reserve(n_samples);
    
    // Initialize random generator
    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    
    // First, identify which vertices can reach max_vertex
    std::unordered_set<size_t> can_reach_max;
    std::queue<size_t> bfs_queue;
    bfs_queue.push(max_vertex);
    can_reach_max.insert(max_vertex);
    
    while (!bfs_queue.empty()) {
        size_t v = bfs_queue.front();
        bfs_queue.pop();
        
        auto it = madag.predecessors.find(v);
        if (it != madag.predecessors.end()) {
            for (size_t u : it->second) {
                if (can_reach_max.find(u) == can_reach_max.end()) {
                    can_reach_max.insert(u);
                    bfs_queue.push(u);
                }
            }
        }
    }
    
    // Sample trajectories by random walks
    for (size_t i = 0; i < n_samples; ++i) {
        std::vector<size_t> path;
        path.push_back(madag.source_vertex);
        
        size_t current = madag.source_vertex;
        while (current != max_vertex) {
            // Get successors that can reach max_vertex
            std::vector<size_t> valid_successors;
            auto it = madag.successors.find(current);
            if (it != madag.successors.end()) {
                for (size_t w : it->second) {
                    if (can_reach_max.find(w) != can_reach_max.end()) {
                        valid_successors.push_back(w);
                    }
                }
            }
            
            if (valid_successors.empty()) {
                // Should not happen if can_reach_max is correct
                break;
            }
            
            // Uniform random selection
            std::uniform_int_distribution<size_t> dist(0, valid_successors.size() - 1);
            current = valid_successors[dist(rng)];
            path.push_back(current);
        }
        
        if (current == max_vertex) {
            madag_trajectory_t traj;
            traj.vertices = std::move(path);
            traj.source_min = madag.source_vertex;
            traj.sink_max = max_vertex;
            traj.total_ascent = y[max_vertex] - y[madag.source_vertex];
            traj.cluster_id = -1;
            
            trajectories.push_back(std::move(traj));
        }
    }
    
    return trajectories;
}

// ============================================================================
// Trajectory Similarity
// ============================================================================

double trajectory_jaccard_similarity(
    const madag_trajectory_t& traj1,
    const madag_trajectory_t& traj2
) {
    std::unordered_set<size_t> set1(traj1.vertices.begin(), traj1.vertices.end());
    std::unordered_set<size_t> set2(traj2.vertices.begin(), traj2.vertices.end());
    
    // Compute intersection size
    size_t intersection_size = 0;
    for (size_t v : set1) {
        if (set2.count(v) > 0) {
            ++intersection_size;
        }
    }
    
    // Compute union size
    size_t union_size = set1.size() + set2.size() - intersection_size;
    
    if (union_size == 0) {
        return 1.0;  // Both empty
    }
    
    return static_cast<double>(intersection_size) / union_size;
}

std::vector<std::vector<double>> compute_trajectory_similarity_matrix(
    const std::vector<madag_trajectory_t>& trajectories,
    const std::string& similarity_type
) {
    size_t n = trajectories.size();
    std::vector<std::vector<double>> sim_matrix(n, std::vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; ++i) {
        sim_matrix[i][i] = 1.0;
        for (size_t j = i + 1; j < n; ++j) {
            double sim = trajectory_jaccard_similarity(trajectories[i], trajectories[j]);
            sim_matrix[i][j] = sim;
            sim_matrix[j][i] = sim;
        }
    }
    
    return sim_matrix;
}

// ============================================================================
// Utility Functions
// ============================================================================

void print_madag_summary(const madag_t& madag, bool verbose) {
    Rprintf("MADAG Summary\n");
    Rprintf("=============\n");
    Rprintf("Source vertex: %zu (y = %.4f)\n", madag.source_vertex, madag.source_value);
    Rprintf("Reachable vertices: %zu\n", madag.reachable_vertices.size());
    Rprintf("Reachable maxima: %zu\n", madag.reachable_maxima.size());
    Rprintf("Cells: %zu\n", madag.cells.size());
    
    if (verbose && !madag.cells.empty()) {
        Rprintf("\nCell details:\n");
        for (const auto& cell : madag.cells) {
            print_cell_summary(cell, false);
        }
    }
}

void print_cell_summary(const ms_cell_t& cell, bool verbose) {
    Rprintf("  %s: support=%zu, trajectories=%zu",
            cell.id_string().c_str(),
            cell.support.size(),
            cell.n_trajectories);
    if (cell.explicitly_enumerated) {
        Rprintf(" (enumerated)");
    }
    if (cell.n_clusters > 0) {
        Rprintf(", clusters=%d", cell.n_clusters);
    }
    Rprintf("\n");
    
    if (verbose) {
        Rprintf("    y range: [%.4f, %.4f], change = %.4f\n",
                cell.min_value, cell.max_value, cell.total_change());
        if (!cell.bottlenecks.empty()) {
            Rprintf("    bottlenecks: ");
            for (size_t i = 0; i < cell.bottlenecks.size(); ++i) {
                if (i > 0) Rprintf(", ");
                Rprintf("%zu", cell.bottlenecks[i]);
            }
            Rprintf("\n");
        }
    }
}
