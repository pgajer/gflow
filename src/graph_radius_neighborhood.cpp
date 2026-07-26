#include "set_wgraph.hpp"

#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

// Generic radius-neighborhood primitive retained for current basin and
// Laplacian diagnostics. Its eventual owner is dgraphs.
std::unordered_map<size_t, double>
set_wgraph_t::find_vertices_within_radius(size_t vertex, double radius) const {
    std::unordered_map<size_t, double> result;
    result[vertex] = 0.0;

    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<
        queue_entry,
        std::vector<queue_entry>,
        std::greater<queue_entry>
    > queue;

    std::vector<double> distances(
        adjacency_list.size(),
        std::numeric_limits<double>::infinity()
    );
    distances[vertex] = 0.0;
    queue.push({0.0, vertex});

    while (!queue.empty()) {
        const auto [current_distance, current_vertex] = queue.top();
        queue.pop();
        if (current_distance > radius) {
            break;
        }
        if (current_distance > distances[current_vertex]) {
            continue;
        }
        result[current_vertex] = current_distance;
        for (const auto& edge : adjacency_list[current_vertex]) {
            const double candidate = current_distance + edge.weight;
            if (candidate < distances[edge.vertex]) {
                distances[edge.vertex] = candidate;
                queue.push({candidate, edge.vertex});
            }
        }
    }

    return result;
}

std::unordered_map<size_t, double>
set_wgraph_t::compute_shortest_path_distances(
    size_t from,
    const std::unordered_set<size_t>& to_set
) const {
    std::vector<double> distances(
        adjacency_list.size(),
        std::numeric_limits<double>::infinity()
    );
    distances[from] = 0.0;

    std::unordered_map<size_t, double> result;
    result.reserve(to_set.size());
    for (const size_t vertex : to_set) {
        result[vertex] = std::numeric_limits<double>::infinity();
    }

    std::unordered_set<size_t> remaining = to_set;
    if (remaining.erase(from) > 0) {
        result[from] = 0.0;
    }

    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<
        queue_entry,
        std::vector<queue_entry>,
        std::greater<queue_entry>
    > queue;
    queue.push({0.0, from});

    while (!queue.empty() && !remaining.empty()) {
        const auto [distance, vertex] = queue.top();
        queue.pop();
        if (distance > distances[vertex]) {
            continue;
        }
        if (remaining.erase(vertex) > 0) {
            result[vertex] = distance;
        }
        for (const auto& edge : adjacency_list[vertex]) {
            const double candidate = distance + edge.weight;
            if (candidate < distances[edge.vertex]) {
                distances[edge.vertex] = candidate;
                queue.push({candidate, edge.vertex});
            }
        }
    }

    return result;
}
