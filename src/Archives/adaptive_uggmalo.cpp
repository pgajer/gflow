#include "adaptive_uggmalo.h"

        #if 0
        reachability_map_t reachability_map = grid_graph.compute_reachability_map(
            grid_vertex,
            max_bw
            );
        #endif

        std::vector<grid_vertex_path_t> grid_vertex_paths = grid_graph.reconstruct_paths(reachability_map);
        std::vector<compose_path_t> composite_paths = grid_graph.create_composite_paths(grid_vertex_paths);


        if (!grid_graph.has_min_size_path(composite_paths, min_path_size, min_num_grid_vertices)) {
            composite_paths = grid_graph.find_min_size_composite_paths(
                grid_vertex,
                min_path_size,
                max_bw,
                result.graph_diameter
                );
        }
