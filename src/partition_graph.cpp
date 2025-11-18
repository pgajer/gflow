#include "partition_graph.hpp"

void compute_partition_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<int>& partition,
    const std::string& weight_type,
    std::vector<std::vector<int>>& adj_list_out,
    std::vector<std::vector<double>>& weight_list_out
) {
    const int n_vertices = adj_list.size();

    // Map partition labels to contiguous cell indices
    std::set<int> unique_labels(partition.begin(), partition.end());
    std::map<int, int> label_to_cell;
    int cell_idx = 0;
    for (int label : unique_labels) {
        label_to_cell[label] = cell_idx++;
    }
    const int n_cells = cell_idx;

    // Invert partition: which vertices belong to each cell
    std::vector<std::vector<int>> cell_vertices(n_cells);
    for (int v = 0; v < n_vertices; ++v) {
        int cell = label_to_cell[partition[v]];
        cell_vertices[cell].push_back(v);
    }

    // Build edge counts between cells using map for efficiency
    std::vector<std::map<int, double>> cell_edges(n_cells);

    // For each vertex, examine its neighbors
    for (int v = 0; v < n_vertices; ++v) {
        int cell_v = label_to_cell[partition[v]];
        const std::vector<int>& neighbors = adj_list[v];
        const std::vector<double>& weights = weight_list[v];

        for (size_t i = 0; i < neighbors.size(); ++i) {
            int u = neighbors[i];
            int cell_u = label_to_cell[partition[u]];

            // Only count edges between different cells
            if (cell_v != cell_u) {
                cell_edges[cell_v][cell_u] += weights[i];
            }
        }
    }

    // Compute normalization factors if needed
    std::vector<int> cell_sizes(n_cells);
    for (int c = 0; c < n_cells; ++c) {
        cell_sizes[c] = cell_vertices[c].size();
    }

    // Compute boundary sizes for Jaccard index
    std::vector<std::map<int, std::set<int>>> boundary_vertices;
    if (weight_type == "jaccard") {
        boundary_vertices.resize(n_cells);

        for (int v = 0; v < n_vertices; ++v) {
            int cell_v = label_to_cell[partition[v]];
            const std::vector<int>& neighbors = adj_list[v];

            for (int u : neighbors) {
                int cell_u = label_to_cell[partition[u]];
                if (cell_v != cell_u) {
                    boundary_vertices[cell_v][cell_u].insert(v);
                }
            }
        }
    }

    // Convert to adjacency list format with appropriate weights
    adj_list_out.resize(n_cells);
    weight_list_out.resize(n_cells);

    for (int c = 0; c < n_cells; ++c) {
        std::vector<std::pair<int, double>> edges;

        for (const auto& kv : cell_edges[c]) {
            int c_neighbor = kv.first;
            double edge_count = kv.second;
            double weight = edge_count;

            if (weight_type == "normalized") {
                // Normalize by geometric mean of cell sizes
                double norm = std::sqrt(
                    static_cast<double>(cell_sizes[c]) *
                    static_cast<double>(cell_sizes[c_neighbor])
                );
                if (norm > 0) {
                    weight = edge_count / norm;
                }
            } else if (weight_type == "jaccard") {
                // Jaccard index based on boundary vertices
                const std::set<int>& boundary_c = boundary_vertices[c][c_neighbor];
                const std::set<int>& boundary_c_neighbor = boundary_vertices[c_neighbor][c];

                // Union size: all vertices in either cell that touch the other cell
                std::set<int> union_set = boundary_c;
                union_set.insert(boundary_c_neighbor.begin(), boundary_c_neighbor.end());

                // Intersection: vertices that are adjacent across the boundary
                // This is approximated by the minimum of the two boundary sizes
                int intersection_size = std::min(boundary_c.size(), boundary_c_neighbor.size());
                int union_size = union_set.size();

                if (union_size > 0) {
                    weight = static_cast<double>(intersection_size) / static_cast<double>(union_size);
                } else {
                    weight = 0.0;
                }
            }
            // else weight_type == "count": use edge_count as is

            if (weight > 0) {
                edges.push_back(std::make_pair(c_neighbor, weight));
            }
        }

        // Sort edges by neighbor index for consistency
        std::sort(edges.begin(), edges.end());

        adj_list_out[c].reserve(edges.size());
        weight_list_out[c].reserve(edges.size());

        for (const auto& edge : edges) {
            adj_list_out[c].push_back(edge.first);
            weight_list_out[c].push_back(edge.second);
        }
    }
}
