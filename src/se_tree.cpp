/**
 * @file se_tree.cpp
 * @brief Implementation of Spurious Extrema tree construction
 */

#include "se_tree.hpp"
#include "set_wgraph.hpp"
#include "gfc.hpp"
#include <queue>
#include <algorithm>
#include <R_ext/Print.h>

/**
 * @brief Build SE tree rooted at a single spurious extremum
 *
 * Starting from a spurious extremum, constructs a tree by:
 * 1. Finding opposite-polarity spurious extrema within the root's basin
 * 2. Recursively expanding from each discovered spurious extremum
 * 3. Terminating branches at non-spurious extrema
 *
 * @param root_vertex The spurious extremum to use as root
 * @param basins Pre-computed basins from compute_gfc_basins()
 * @param spurious_min Set of spurious minimum vertices
 * @param spurious_max Set of spurious maximum vertices
 * @param verbose Print progress information
 * @return SE tree structure
 */
se_tree_t set_wgraph_t::build_se_tree(
    size_t root_vertex,
    const std::unordered_map<size_t, bbasin_t>& basins,
    const std::unordered_set<size_t>& spurious_min,
    const std::unordered_set<size_t>& spurious_max,
    bool verbose
) const {

    se_tree_t tree;
    tree.root_vertex = root_vertex;

    // Determine root type
    auto root_it = basins.find(root_vertex);
    if (root_it == basins.end()) {
        if (verbose) {
            Rprintf("Warning: root vertex %zu not found in basins\n", root_vertex + 1);
        }
        tree.classification = se_tree_class_t::NO_TERMINALS;
        return tree;
    }

    tree.root_is_maximum = root_it->second.is_maximum;

    // Verify root is spurious
    bool root_is_spurious = tree.root_is_maximum ?
        (spurious_max.count(root_vertex) > 0) :
        (spurious_min.count(root_vertex) > 0);

    if (!root_is_spurious) {
        if (verbose) {
            Rprintf("Warning: root vertex %zu is not spurious\n", root_vertex + 1);
        }
        tree.classification = se_tree_class_t::NO_TERMINALS;
        return tree;
    }

    // Track discovered vertices to avoid cycles
    std::unordered_set<size_t> discovered;

    // BFS queue: (vertex, parent_vertex)
    std::queue<std::pair<size_t, size_t>> queue;
    queue.push({root_vertex, SIZE_MAX});
    discovered.insert(root_vertex);

    while (!queue.empty()) {
        auto [current_vertex, parent_vertex] = queue.front();
        queue.pop();

        // Get basin for current vertex
        auto basin_it = basins.find(current_vertex);
        if (basin_it == basins.end()) {
            continue;
        }

        const bbasin_t& current_basin = basin_it->second;

        // Create node for current vertex
        se_node_t node;
        node.vertex = current_vertex;
        node.is_maximum = current_basin.is_maximum;
        node.parent = parent_vertex;
        node.basin_vertices = current_basin.vertices;
        node.basin_boundary = current_basin.boundary;

        // Determine if current is spurious
        node.is_spurious = node.is_maximum ?
            (spurious_max.count(current_vertex) > 0) :
            (spurious_min.count(current_vertex) > 0);

        // If non-spurious, this is a terminal (leaf)
        if (!node.is_spurious) {
            if (node.is_maximum) {
                tree.ns_max_terminals.push_back(current_vertex);
            } else {
                tree.ns_min_terminals.push_back(current_vertex);
            }
            tree.nodes[current_vertex] = node;
            continue;
        }

        // Find opposite-polarity extrema within this basin
        const std::unordered_set<size_t>& target_set =
            node.is_maximum ? spurious_min : spurious_max;
        const std::unordered_set<size_t>& ns_target_set =
            node.is_maximum ?
                std::unordered_set<size_t>{} :  // Will check individually
                std::unordered_set<size_t>{};

        // Convert basin vertices to set for fast lookup
        std::unordered_set<size_t> basin_vertex_set(
            current_basin.vertices.begin(),
            current_basin.vertices.end()
        );

        // Find all opposite-polarity extrema in this basin
        for (size_t v : current_basin.vertices) {
            // Skip if already discovered
            if (discovered.count(v)) continue;

            // Check if v is an extremum of opposite polarity
            auto v_basin_it = basins.find(v);
            if (v_basin_it == basins.end()) continue;

            // Must be opposite polarity
            if (v_basin_it->second.is_maximum == node.is_maximum) continue;

            // Check if spurious or non-spurious
            bool v_is_spurious = v_basin_it->second.is_maximum ?
                (spurious_max.count(v) > 0) : (spurious_min.count(v) > 0);

            // Add as child and queue for expansion
            node.children.push_back(v);
            discovered.insert(v);
            queue.push({v, current_vertex});

            if (verbose) {
                Rprintf("  %s %zu -> %s %zu (%s)\n",
                        node.is_maximum ? "SMAX" : "SMIN", current_vertex + 1,
                        v_basin_it->second.is_maximum ? "MAX" : "MIN", v + 1,
                        v_is_spurious ? "spurious" : "non-spurious");
            }
        }

        tree.nodes[current_vertex] = node;
    }

    // Classify tree
    bool has_ns_min = !tree.ns_min_terminals.empty();
    bool has_ns_max = !tree.ns_max_terminals.empty();

    if (has_ns_min && has_ns_max) {
        tree.classification = se_tree_class_t::BOTH_TYPES;
    } else if (has_ns_min) {
        tree.classification = se_tree_class_t::ONLY_MIN;
    } else if (has_ns_max) {
        tree.classification = se_tree_class_t::ONLY_MAX;
    } else {
        tree.classification = se_tree_class_t::NO_TERMINALS;
    }

    if (verbose) {
        Rprintf("SE tree rooted at %s %zu: %zu nodes, %zu NSMIN terminals, %zu NSMAX terminals\n",
                tree.root_is_maximum ? "SMAX" : "SMIN", root_vertex + 1,
                tree.nodes.size(),
                tree.ns_min_terminals.size(),
                tree.ns_max_terminals.size());
        Rprintf("Classification: %s\n",
                tree.classification == se_tree_class_t::BOTH_TYPES ? "BOTH_TYPES" :
                tree.classification == se_tree_class_t::ONLY_MIN ? "ONLY_MIN" :
                tree.classification == se_tree_class_t::ONLY_MAX ? "ONLY_MAX" :
                "NO_TERMINALS");
    }

    return tree;
}

/**
 * @brief Extract path from root to a terminal in the SE tree
 *
 * @param tree The SE tree
 * @param terminal_vertex The terminal (non-spurious) vertex
 * @return Vector of vertices from root to terminal
 */
std::vector<size_t> extract_se_tree_path(
    const se_tree_t& tree,
    size_t terminal_vertex
) {
    std::vector<size_t> path;

    // Trace back from terminal to root
    size_t current = terminal_vertex;
    while (current != SIZE_MAX) {
        path.push_back(current);
        auto it = tree.nodes.find(current);
        if (it == tree.nodes.end()) break;
        current = it->second.parent;
    }

    // Reverse to get root-to-terminal order
    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Compute harmonic repair support region from SE tree
 *
 * For BOTH_TYPES trees: Union of basins along paths to one NSMIN and one NSMAX
 * For ONLY_MIN/ONLY_MAX trees: Union of basins along path to nearest NS terminal
 *
 * @param tree The SE tree (modified in place to set hr_support fields)
 * @param verbose Print progress information
 */
void compute_se_tree_hr_support(
    se_tree_t& tree,
    bool verbose
) {
    std::unordered_set<size_t> support_set;
    std::unordered_set<size_t> boundary_set;

    std::vector<size_t> paths_to_use;

    switch (tree.classification) {
        case se_tree_class_t::BOTH_TYPES: {
            // Use shortest path to NSMIN and shortest path to NSMAX
            // For now, just use the first of each
            if (!tree.ns_min_terminals.empty()) {
                auto path = extract_se_tree_path(tree, tree.ns_min_terminals[0]);
                paths_to_use.insert(paths_to_use.end(), path.begin(), path.end());
            }
            if (!tree.ns_max_terminals.empty()) {
                auto path = extract_se_tree_path(tree, tree.ns_max_terminals[0]);
                paths_to_use.insert(paths_to_use.end(), path.begin(), path.end());
            }
            break;
        }

        case se_tree_class_t::ONLY_MIN: {
            if (!tree.ns_min_terminals.empty()) {
                auto path = extract_se_tree_path(tree, tree.ns_min_terminals[0]);
                paths_to_use.insert(paths_to_use.end(), path.begin(), path.end());
            }
            break;
        }

        case se_tree_class_t::ONLY_MAX: {
            if (!tree.ns_max_terminals.empty()) {
                auto path = extract_se_tree_path(tree, tree.ns_max_terminals[0]);
                paths_to_use.insert(paths_to_use.end(), path.begin(), path.end());
            }
            break;
        }

        case se_tree_class_t::NO_TERMINALS: {
            // Use all nodes in tree
            for (const auto& [vertex, node] : tree.nodes) {
                paths_to_use.push_back(vertex);
            }
            break;
        }
    }

    // Remove duplicates from paths_to_use
    std::sort(paths_to_use.begin(), paths_to_use.end());
    paths_to_use.erase(
        std::unique(paths_to_use.begin(), paths_to_use.end()),
        paths_to_use.end()
    );

    // Collect union of basins along selected paths
    for (size_t vertex : paths_to_use) {
        auto it = tree.nodes.find(vertex);
        if (it == tree.nodes.end()) continue;

        const se_node_t& node = it->second;

        // Add basin vertices to support
        for (size_t v : node.basin_vertices) {
            support_set.insert(v);
        }
    }

    // Compute boundary: neighbors of support that are not in support
    // This requires access to the graph, so we store basin boundaries
    // and compute the union boundary from those
    for (size_t vertex : paths_to_use) {
        auto it = tree.nodes.find(vertex);
        if (it == tree.nodes.end()) continue;

        const se_node_t& node = it->second;

        // Add boundary vertices that are not in support
        for (size_t v : node.basin_boundary) {
            if (support_set.find(v) == support_set.end()) {
                boundary_set.insert(v);
            }
        }
    }

    // Convert to vectors
    tree.hr_support_vertices.assign(support_set.begin(), support_set.end());
    std::sort(tree.hr_support_vertices.begin(), tree.hr_support_vertices.end());

    tree.hr_support_boundary.assign(boundary_set.begin(), boundary_set.end());
    std::sort(tree.hr_support_boundary.begin(), tree.hr_support_boundary.end());

    if (verbose) {
        Rprintf("HR support: %zu interior vertices, %zu boundary vertices\n",
                tree.hr_support_vertices.size(),
                tree.hr_support_boundary.size());
    }
}
