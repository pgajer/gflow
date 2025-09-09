#include "A_path_graph.h"


#if 0
path_graph_plm_t create_path_graph_plm(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    Rprintf("\nStarting path graph creation\n");
    const int n_vertices = adj_list.size();
    Rprintf("Number of vertices: %d\n", n_vertices);

    // Print input adjacency list and weights
    Rprintf("\nInput adjacency lists:\n");
    for (int i = 0; i < n_vertices; ++i) {
        Rprintf("Vertex %d connected to:", i);
        for (int j : adj_list[i]) {
            Rprintf(" %d", j);
        }
        Rprintf("\n");
    }

    path_graph_plm_t result;
    result.h = h;
    result.adj_list.resize(n_vertices);
    result.weight_list.resize(n_vertices);
    result.hop_list.resize(n_vertices);

    std::vector<double> distances(n_vertices);
    std::vector<int> hops(n_vertices);
    std::vector<int> parent(n_vertices);
    std::vector<bool> in_queue(n_vertices);
    std::vector<neighbor_info_t> current_neighbors;
    current_neighbors.reserve(n_vertices / 2);

    for (int start = 0; start < n_vertices; ++start) {
        Rprintf("\nProcessing paths from vertex %d\n", start);

        std::fill(distances.begin(), distances.end(), std::numeric_limits<double>::infinity());
        std::fill(hops.begin(), hops.end(), std::numeric_limits<int>::max());
        std::fill(parent.begin(), parent.end(), -1);
        std::fill(in_queue.begin(), in_queue.end(), false);
        current_neighbors.clear();

        distances[start] = 0;
        hops[start] = 0;

        std::deque<int> q{start};
        in_queue[start] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop_front();
            in_queue[current] = false;

            if (hops[current] >= h) {
                Rprintf("  Reached hop limit at vertex %d\n", current);
                continue;
            }

            Rprintf("  Processing vertex %d (distance=%.2f, hops=%d)\n",
                   current, distances[current], hops[current]);

            const double current_dist = distances[current];
            const int current_hops = hops[current];
            const int next_hops = current_hops + 1;

            const auto& current_adj = adj_list[current];
            const auto& current_weights = weight_list[current];
            const size_t n_neighbors = current_adj.size();

            for (size_t i = 0; i < n_neighbors; ++i) {
                const int neighbor = current_adj[i];
                const double new_distance = current_dist + current_weights[i];

                Rprintf("    Checking neighbor %d (new_distance=%.2f, current_best=%.2f)\n",
                       neighbor, new_distance, distances[neighbor]);

                if (new_distance < distances[neighbor]) {
                    Rprintf("    Updating path to %d: distance %.2f->%.2f, hops %d->%d\n",
                           neighbor, distances[neighbor], new_distance,
                           hops[neighbor], next_hops);

                    distances[neighbor] = new_distance;
                    hops[neighbor] = next_hops;
                    parent[neighbor] = current;

                    if (!in_queue[neighbor] && next_hops < h) {
                        q.push_back(neighbor);
                        in_queue[neighbor] = true;
                    }
                }
            }
        }

        // Store paths
        for (int v = 0; v < n_vertices; ++v) {
            if (v != start && hops[v] <= h) {
                current_neighbors.emplace_back(v, distances[v], hops[v]);

                if (start < v) {
                    Rprintf("\nStoring path from %d to %d:\n", start, v);
                    auto& path = result.shortest_paths[{start, v}];
                    path.clear();
                    path.reserve(hops[v] + 1);

                    // Print path reconstruction
                    Rprintf("  Path: ");
                    for (int current = v; current != -1; current = parent[current]) {
                        path.push_back(current);
                        Rprintf("%d ", current);
                    }

                    std::reverse(path.begin(), path.end());

                    Rprintf("\n  Final path after reverse: ");
                    for (int vertex : path) {
                        Rprintf("%d ", vertex);
                    }
                    Rprintf("\n");
                }
            }
        }

        Rprintf("\nPaths found from vertex %d:\n", start);
        for (const auto& info : current_neighbors) {
            Rprintf("  To %d: distance=%.2f, hops=%d\n",
                   info.vertex, info.distance, info.hops);
        }
    }

    result.vertex_paths.resize(n_vertices);

    // Print final shortest paths
    Rprintf("\nFinal shortest paths:\n");
    for (const auto& path_entry : result.shortest_paths) {
        Rprintf("Path %d->%d: ", path_entry.first.first, path_entry.first.second);
        for (int v : path_entry.second) {
            Rprintf("%d ", v);
        }
        Rprintf("\n");
    }

    return result;
}
#endif
