/**
 * @brief Creates a uniform grid graph based on a maximal ε-packing
 *
 * @param input_adj_list Input graph adjacency list
 * @param input_weight_list Input graph edge weights
 * @param desired_grid_size Approximate number of grid points desired
 * @param start_vertex Starting vertex for the first grid point
 *
 * @return uniform_grid_graph_t A graph with grid_vertices marking the ε-packing
 */

double estimate_diameter(
	const std::vector<std::vector<int>>& adj_list,
	const std::vector<std::vector<double>>& weight_list,
	size_t start_vertex);

uniform_grid_graph_t create_maximal_packing_grid(
	const std::vector<std::vector<int>>& input_adj_list,
	const std::vector<std::vector<double>>& input_weight_list,
	size_t desired_grid_size,
	size_t start_vertex) {

// Initialize the graph directly using the inherited constructor
	uniform_grid_graph_t grid_graph(input_adj_list, input_weight_list);

// Estimate the graph diameter using Dijkstra from the start vertex
	double graph_diameter = estimate_diameter(input_adj_list, input_weight_list, start_vertex);
	Rprintf("graph_diameter: %f\n", graph_diameter);

	graph_diameter = get_vertex_eccentricity(input_adj_list, input_weight_list, start_vertex);
	Rprintf("graph_diameter: %f\n", graph_diameter);

// Estimate epsilon based on desired grid size
	double epsilon = graph_diameter / desired_grid_size;
	Rprintf("epsilon: %f\n", epsilon);


// Create the maximal ε-packing
	std::unordered_set<size_t> packing;
	packing.insert(start_vertex);

	std::vector<bool> explored(input_adj_list.size(), false); // vertices that are within epsilon radius of all already found packing vertices
	explored[start_vertex] = true;

	size_t search_start = start_vertex;
	while (true) {
		size_t new_vertex = find_distant_vertex(
			input_adj_list, input_weight_list, packing, epsilon, search_start);

		if (new_vertex == SIZE_MAX) {
// No more vertices can be added to the packing at this epsilon
			break;
		}

		packing.insert(new_vertex);
		search_start = new_vertex; // Start next search from the newly added vertex
	}

// Set the grid vertices in the output graph
	grid_graph.grid_vertices = packing;

	return grid_graph;
}

/**
 * @brief Estimates the diameter of a graph using Dijkstra's algorithm
 *
 * @param adj_list Adjacency list representation of the graph
 * @param weight_list Edge weights corresponding to the adjacency list
 * @param start_vertex Vertex to start the search from
 *
 * @return double Estimated diameter of the graph
 */
double estimate_diameter(
	const std::vector<std::vector<int>>& adj_list,
	const std::vector<std::vector<double>>& weight_list,
	size_t start_vertex) {

	size_t n = adj_list.size();

// Run Dijkstra from start_vertex
	std::vector<double> dist(n, std::numeric_limits<double>::infinity());
	std::vector<bool> visited(n, false);
	dist[start_vertex] = 0.0;

	std::priority_queue<std::pair<double, size_t>,
						std::vector<std::pair<double, size_t>>,
						std::greater<std::pair<double, size_t>>> pq;
	pq.push({0.0, start_vertex});

	while (!pq.empty()) {
		double curr_dist = pq.top().first;
		size_t curr = pq.top().second;
		pq.pop();

		if (visited[curr]) {
			continue;
		}

		visited[curr] = true;

		for (size_t i = 0; i < adj_list[curr].size(); ++i) {
			int next = adj_list[curr][i];
			double weight = weight_list[curr][i];

			if (!visited[next] && curr_dist + weight < dist[next]) {
				dist[next] = curr_dist + weight;
				pq.push({dist[next], next});
			}
		}
	}

// Find the furthest vertex from start_vertex
	size_t furthest_vertex = start_vertex;
	double max_dist = 0.0;

	for (size_t i = 0; i < n; ++i) {
		if (dist[i] != std::numeric_limits<double>::infinity() && dist[i] > max_dist) {
			max_dist = dist[i];
			furthest_vertex = i;
		}
	}

// Run Dijkstra from the furthest vertex
	std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::infinity());
	std::fill(visited.begin(), visited.end(), false);
	dist[furthest_vertex] = 0.0;

	pq.push({0.0, furthest_vertex});

	while (!pq.empty()) {
		double curr_dist = pq.top().first;
		size_t curr = pq.top().second;
		pq.pop();

		if (visited[curr]) {
			continue;
		}

		visited[curr] = true;

		for (size_t i = 0; i < adj_list[curr].size(); ++i) {
			int next = adj_list[curr][i];
			double weight = weight_list[curr][i];

			if (!visited[next] && curr_dist + weight < dist[next]) {
				dist[next] = curr_dist + weight;
				pq.push({dist[next], next});
			}
		}
	}

// Find the maximum distance, which is an estimate of the diameter
	double diameter = 0.0;
	for (size_t i = 0; i < n; ++i) {
		if (dist[i] != std::numeric_limits<double>::infinity() && dist[i] > diameter) {
			diameter = dist[i];
		}
	}

return diameter;
}
