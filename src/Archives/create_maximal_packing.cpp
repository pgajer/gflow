#include "create_maximal_packing.h"

void set_wgraph_t::trace_exploration_path(
    size_t from_vertex,
    size_t to_vertex,
    double radius) const {

    // Rerun Dijkstra's algorithm with detailed logging
    std::vector<double> dist(adjacency_list.size(), INFINITY);
    std::vector<size_t> prev(adjacency_list.size(), INVALID_VERTEX);
    dist[from_vertex] = 0;

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, from_vertex});

    Rprintf("TRACE: Starting path exploration from %zu to %zu\n", from_vertex + 1, to_vertex + 1);

    while (!pq.empty()) {
        auto top = pq.top();
        double d = -top.first;  // Note the negation here!
        size_t u = top.second;
        pq.pop();

        if (u == to_vertex) {
            Rprintf("TRACE: Reached target vertex %zu with distance %.4f\n", to_vertex + 1, d);
            break;
        }

        if (d > dist[u]) continue;

        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double w = edge.weight;
            double new_dist = dist[u] + w;

            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                pq.push({-new_dist, v});  // Note the negation!

                if (v == to_vertex || new_dist < radius) {
                    Rprintf("TRACE: Updated vertex %zu distance to %.4f (via %zu)\n",
                            v + 1, new_dist, u + 1);
                }

                // Key debugging for the issue
                if (new_dist < radius) {
                    Rprintf("TRACE: Vertex %zu should be marked as explored (dist=%.4f < radius=%.4f)\n",
                            v + 1, new_dist, radius);
                }
            }
        }
    }

    // Reconstruct and print the shortest path
    if (prev[to_vertex] != INVALID_VERTEX) {
        std::vector<size_t> path;
        for (size_t at = to_vertex; at != from_vertex; at = prev[at]) {
            path.push_back(at);
        }
        path.push_back(from_vertex);
        std::reverse(path.begin(), path.end());

        Rprintf("TRACE: Shortest path from %zu to %zu: ", from_vertex + 1, to_vertex + 1);
        for (size_t i = 0; i < path.size(); i++) {
            Rprintf("%zu", path[i] + 1);
            if (i < path.size() - 1) {
                Rprintf(" -> ");
            }
        }
        Rprintf(" (total distance: %.4f)\n", dist[to_vertex]);
    }
}

std::vector<size_t> set_wgraph_t::create_maximal_packing(
    double radius,
    size_t start_vertex) const {

    size_t n_vertices = adjacency_list.size();
    std::vector<size_t> packing;
    packing.reserve(n_vertices);
    packing.push_back(start_vertex);

    Rprintf("DEBUG: Started with packing vertex %zu\n", start_vertex + 1);

    explored_tracker_t explored_tracker(n_vertices, start_vertex);
    std::queue<size_t> boundary_queue;

    // Initial boundary vertices from the starting point
    std::vector<size_t> boundary_vertices = find_boundary_vertices_outside_radius(start_vertex, radius, explored_tracker);

    Rprintf("DEBUG: Found %zu initial boundary vertices\n", boundary_vertices.size());

    // Add all initial boundary vertices to the queue
    for (size_t vertex : boundary_vertices) {
        boundary_queue.push(vertex);
    }

    // Main packing algorithm loop
    int iteration = 0;
    while (!boundary_queue.empty() && explored_tracker.n_explored < n_vertices) {
        iteration++;

        // Get the next candidate vertex from the queue
        size_t current_vertex = boundary_queue.front();
        boundary_queue.pop();

        // Skip if this vertex has already been explored
        if (explored_tracker.is_explored(current_vertex)) {
            Rprintf("DEBUG [%d]: Vertex %zu already explored, skipping\n", iteration, current_vertex + 1);
            continue;
        }

        Rprintf("DEBUG [%d]: Examining vertex %zu for packing (explored so far: %zu/%zu)\n",
                iteration, current_vertex + 1, explored_tracker.n_explored, n_vertices);

        // Check if vertex maintains minimum distance from all packing vertices
        bool violation_found = false;
        for (size_t i = 0; i < packing.size(); i++) {
            size_t packing_vertex = packing[i];

            // Skip self-comparison
            if (packing_vertex == current_vertex) continue;

            // Compute shortest path distance
            double distance = compute_shortest_path_distance(current_vertex, packing_vertex);

            if (distance < radius) {
                // CRITICAL ERROR: We found a vertex that should have been marked as explored
                Rprintf("ERROR [%d]: Vertex %zu (dist=%.4f) is too close to packing vertex %zu\n",
                        iteration, current_vertex + 1, distance, packing_vertex + 1);
                Rprintf("ERROR [%d]: This indicates a failure in the explored_tracker mechanism\n", iteration);
                Rprintf("ERROR [%d]: When processing packing vertex %zu, vertex %zu should have been marked as explored\n",
                        iteration, packing_vertex + 1, current_vertex + 1);

                // Let's trace the exploration path from packing_vertex to current_vertex
                Rprintf("DEBUG [%d]: Tracing path from packing vertex %zu to current vertex %zu\n",
                        iteration, packing_vertex + 1, current_vertex + 1);

                // You could add a separate function here to trace the path
                trace_exploration_path(packing_vertex, current_vertex, radius);

                violation_found = true;

				REPORT_ERROR("DEBUGGING\n");

                break;
            }
        }

        // For debugging only: add the vertex even if it's too close
        // This will help us identify more issues in a single run
        if (violation_found) {
            Rprintf("DEBUG [%d]: Adding vertex %zu to packing DESPITE VIOLATION (for debugging)\n",
                    iteration, current_vertex + 1);
        } else {
            Rprintf("DEBUG [%d]: Adding vertex %zu to packing (valid addition)\n",
                    iteration, current_vertex + 1);
        }

        packing.push_back(current_vertex);
        explored_tracker.explored[current_vertex] = true;
        explored_tracker.n_explored++;

        // Find new boundary vertices from this new packing vertex
        boundary_vertices = find_boundary_vertices_outside_radius(current_vertex, radius, explored_tracker);

        Rprintf("DEBUG [%d]: Found %zu boundary vertices from vertex %zu\n",
                iteration, boundary_vertices.size(), current_vertex + 1);

        // Add newly discovered boundary vertices to the queue
        for (size_t vertex : boundary_vertices) {
            if (!explored_tracker.is_explored(vertex)) {
                boundary_queue.push(vertex);
            }
        }
    }

    return packing;
}
