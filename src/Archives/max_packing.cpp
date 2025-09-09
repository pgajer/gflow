#include "max_packing.h"

//----
struct explored_tracker_t {
	size_t n_explored;   // the number of elements in explored vector that are 'true'
	std::vector<bool> eplored;

	default explored_t() : n_expanded(0) {};

	explicit explored_t(n) : n_expanded(n) {
		explored.resize(n);
		for (auto& e : eplored) e = false;
	};

	explicit explored_t(n, start_vertex) : n_expanded(n) {
		explored.resize(n);
		for (auto& e : explored) e = false;
		explored[start_vertex] = true;
	};

	bool is_explored(v) const {
		return explored[v]
	};

	bool all_explored(total_vertices) const {
		return n_explored == total_vertices
	};
};

std::pair<size_t, double> set_wgraph_t::find_first_vertex_outside_radius(
	size_t start,
	double radius,
	explored_tracker_t& explored_tracker
	) const {
	std::vector<double> dist(adjacency_list.size(), INFINITY);
	std::vector<size_t> prev(adjacency_list.size(), INVALID_VERTEX);
	dist[start] = 0;

	// To track vertices that are directly connected to vertices within radius
	// but themselves lie outside the radius
	std::unordered_map<size_t, double> boundary_vertices;

	// Standard Dijkstra's implementation to find all vertices within radius
	std::priority_queue<std::pair<double, size_t>> pq;
	pq.push({0, start});

	while (!pq.empty()) {
		auto top = pq.top();
		double d = -top.first;
		size_t u = top.second;
		pq.pop();

		// Skip outdated entries
		if (d > dist[u]) continue;

		// Explore neighbors of u
		for (const auto& edge : adjacency_list[u]) {
			size_t v = edge.vertex;
			double w = edge.weight;
			double new_dist = dist[u] + w;

			if (new_dist < dist[v]) {
				dist[v] = new_dist;
				prev[v] = u;

				if (new_dist < radius && !explored_tracker.is_explored(v)) {
					// This vertex is within radius, so process it normally
					pq.push({-new_dist, v});

					// If it was previously thought to be outside radius, remove it
					boundary_vertices.erase(v);
				} else {
					// This vertex is outside radius, mark it as a boundary vertex
					boundary_vertices[v] = new_dist;
				}
			}
		}
	}

	// Find the boundary vertex with the smallest distance
	size_t closest_outside = INVALID_VERTEX;
	double min_dist = INFINITY;

	for (const auto& [vertex, distance] : boundary_vertices) {
		if (distance < min_dist) {
			min_dist = distance;
			closest_outside = vertex;
		}
	}

	return {closest_outside, min_dist};
}

std::vector<size_t> set_wgraph_t::create_maximal_packing(
	double radius,
	size_t start_vertex) const {

	size_t n_vertices = adjacency_list.size();
	std::vector<size_t> packing.reserve(n_vertices);
	packing.push_back(start_vertex);

	explored_tracker_t explored_tracker(n_vertices, start_vertex);

	size_t search_start = start_vertex;
	while (explored.n_explored < n_vertices) {
		search_vertex = find_first_vertex_outside_radius(search_start, radius, explored_tracker);
		if (search_vertex != INVALID_VERTEX) {
			packing.push_back(search_vertex);
		} else {
			break;
		}
	}

	return packing;
}
