#ifndef PATH_GRAPHS_H_
#define PATH_GRAPHS_H_

#include <vector>   // for std::vector
#include <map>      // for std::map
#include <utility>  // for std::pair
#include <algorithm>  // for std::min, std::max

// Basic path graph structure
struct path_graph_t {
    std::vector<std::vector<int>> adj_list;
    std::vector<std::vector<double>> weight_list;
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths;
};

// Augmented path graph structure with hop information
struct path_graph_plus_t {
    std::vector<std::vector<int>> adj_list;
    std::vector<std::vector<double>> weight_list;
    std::vector<std::vector<int>> hop_list;
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths;
};

struct vertex_path_info_t {
    std::vector<std::pair<int,int>> containing_paths;  // List of (start,end) pairs for paths that contain the given vertex
    std::vector<int> position_in_path;                 // Position of vertex in each path

    // Helper method to get full paths of specific length
    std::vector<std::pair<std::vector<int>, int>>
    get_full_paths(int desired_length,
                  const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths) const {
        std::vector<std::pair<std::vector<int>, int>> result;

        for (size_t i = 0; i < containing_paths.size(); ++i) {
            const auto& [start, end] = containing_paths[i];
            auto path_it = shortest_paths.find({std::min(start, end), std::max(start, end)});
            if (path_it != shortest_paths.end()) {
                const auto& path = path_it->second;
                if (path.size() == desired_length) {
                    result.emplace_back(path, position_in_path[i]);
                }
            }
        }
        return result;
    }
};

// Modified path_graph_plus_t
struct path_graph_plm_t {
    int h;
    std::vector<std::vector<int>> adj_list;
    std::vector<std::vector<double>> weight_list;
    std::vector<std::vector<int>> hop_list;
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths;
    std::vector<vertex_path_info_t> vertex_paths;

    // Helper method to get paths containing a vertex
    std::vector<std::pair<std::vector<int>, int>>
    get_paths_containing(int vertex, int desired_length) const {
        return vertex_paths[vertex].get_full_paths(desired_length, shortest_paths);
    }
};

#endif // PATH_GRAPHS_H_
