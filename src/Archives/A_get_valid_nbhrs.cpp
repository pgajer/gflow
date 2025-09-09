#include "A_get_valid_nbhrs.h"

    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<hop_neighbor_info_t> neighbors;
        int min_valid_hop = -1;

        // First pass: check for extrema in neighbors and find minimum hop distance
        for (size_t i = 0; i < adj_list[current_vertex].size(); i++) {
            int neighbor = adj_list[current_vertex][i];
            int hop_dist = hop_list[current_vertex][i];

            bool neighbor_is_max = ms_cx.local_maxima.count(neighbor);
            bool neighbor_is_min = ms_cx.local_minima.count(neighbor);

            // If we find a local maximum while ascending, return it as the only neighbor
            if (ascending && neighbor_is_max) {
                neighbors.clear();
                auto path_key = std::make_pair(
                    std::min(current_vertex, neighbor),
                    std::max(current_vertex, neighbor)
                    );

                auto path_it = shortest_paths.find(path_key);
                if (path_it == shortest_paths.end()) {
                    Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                            current_vertex, neighbor);
                    error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                }

                const auto& path = path_it->second;
                if (path.empty()) {
                    Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                            current_vertex, neighbor);
                    error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                }

                neighbors.emplace_back(hop_neighbor_info_t{
                        neighbor,
                        hop_dist,
                        Ey[neighbor] - Ey[current_vertex],
                        path_it->second,
                        true,   // is_lmax
                        false   // is_lmin
                    });

                return neighbors;
            }

            // If we find a local minimum while descending, return it as the only neighbor
            if (!ascending && neighbor_is_min) {
                neighbors.clear();
                auto path_key = std::make_pair(
                    std::min(current_vertex, neighbor),
                    std::max(current_vertex, neighbor)
                    );

                auto path_it = shortest_paths.find(path_key);
                if (path_it == shortest_paths.end()) {
                    Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                            current_vertex, neighbor);
                    error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                }

                const auto& path = path_it->second;
                if (path.empty()) {
                    Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                            current_vertex, neighbor);
                    error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                }

                neighbors.emplace_back(hop_neighbor_info_t{
                        neighbor,
                        hop_dist,
                        Ey[current_vertex] - Ey[neighbor],
                        path_it->second,
                        false,  // is_lmax
                        true   // is_lmin
                    });

                return neighbors;
            }

            // For non-extrema neighbors, track minimum hop distance
            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0) {
                if (min_valid_hop == -1 || hop_dist < min_valid_hop) {
                    min_valid_hop = hop_dist;
                }
            }
        }

        // Second pass: collect all neighbors at minimum hop distance
        if (min_valid_hop != -1) {
            for (size_t i = 0; i < adj_list[current_vertex].size(); i++) {
                int neighbor = adj_list[current_vertex][i];
                int hop_dist = hop_list[current_vertex][i];

                if (hop_dist == min_valid_hop) {
                    double diff = ascending ?
                        Ey[neighbor] - Ey[current_vertex] :
                        Ey[current_vertex] - Ey[neighbor];

                    if (diff > 0) {
                        auto path_key = std::make_pair(
                            std::min(current_vertex, neighbor),
                            std::max(current_vertex, neighbor)
                            );

                        auto path_it = shortest_paths.find(path_key);
                        if (path_it == shortest_paths.end()) {
                            Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                                    current_vertex, neighbor);
                            error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                        }

                        const auto& path = path_it->second;
                        if (path.empty()) {
                            Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                                    current_vertex, neighbor);
                            error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                        }

                        neighbors.emplace_back(hop_neighbor_info_t{
                                neighbor,
                                hop_dist,
                                diff,
                                path_it->second,
                                ms_cx.local_maxima.count(neighbor) > 0,  // is_lmax
                                ms_cx.local_minima.count(neighbor) > 0   // is_lmin
                            });
                    }
                }
            }

            // Only sort if we have multiple neighbors
            if (neighbors.size() > 1) {
                // Sort by difference value in descending order
                std::sort(neighbors.begin(), neighbors.end(),
                          [](const hop_neighbor_info_t& a, const hop_neighbor_info_t& b) {
                              return a.diff_value > b.diff_value;
                          });
            }
        }

        return neighbors;
    };
