#ifndef UNIFORM_GRID_GRAPH_H_
#define UNIFORM_GRID_GRAPH_H_

// Using a pair to store both vertex and weight information
using edge_info_t = std::pair<int, double>;  // (neighbor_vertex, edge_weight)

// Custom comparator for edge_info_t to ensure proper ordering in the set
struct edge_info_compare_t {
    bool operator()(const edge_info_t& a, const edge_info_t& b) const {
        // Order by vertex index only as each edge is uniquely defined by its end vertices
        return a.first < b.first;
    }
};

struct uniform_grid_graph_t {
    std::vector<std::set<edge_info_t, edge_info_compare_t>> ugg;
    std::unordered_set<int> grid_vertices;
    int n_original_vertices;
    //size_t size() const { return ugg.size(); }
};

uniform_grid_graph_t create_uniform_grid_graph(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    int grid_size,
    int start_vertex = 0);

#endif // UNIFORM_GRID_GRAPH_H_


#if 0
// example of uniform_grid_graph_t with dist and prev members
class uniform_grid_graph_t : public set_wgraph_t {
private:
    std::vector<double> dist;
    std::vector<size_t> prev;

public:
    struct radius_bounded_paths_t {
        std::vector<vertex_info_t> sorted_vertices;
        size_t ref_vertex;
    };

    uniform_grid_graph_t(/* ... */) {
        // Initialize vectors with appropriate size
        dist.resize(adjacency_list.size(), INFINITY);
        prev.resize(adjacency_list.size(), INVALID_VERTEX);
    }

    radius_bounded_paths_t compute_reachability_map(size_t ref_vertex, double radius) {
        // Reset the arrays instead of creating new ones
        std::fill(dist.begin(), dist.end(), INFINITY);
        std::fill(prev.begin(), prev.end(), INVALID_VERTEX);

        // Implement Dijkstra's algorithm using member arrays
        // ... implementation details ...

        // Return only the sorted vertices and reference vertex
        radius_bounded_paths_t result;
        result.ref_vertex = ref_vertex;

        // Populate sorted_vertices based on dist and prev arrays
        for (size_t v = 0; v < dist.size(); ++v) {
            if (v != ref_vertex &&
                is_original_vertex(v, n_original_vertices) &&
                dist[v] <= radius) {
                result.sorted_vertices.push_back({v, dist[v]});
            }
        }

        std::sort(result.sorted_vertices.begin(), result.sorted_vertices.end(),
                 [](const vertex_info_t& a, const vertex_info_t& b) {
                     return a.distance > b.distance;
                 });

        return result;
    }

    // Add methods to access distance and predecessor information if needed
    double get_distance(size_t vertex) const { return dist[vertex]; }
    size_t get_predecessor(size_t vertex) const { return prev[vertex]; }
};

#endif
