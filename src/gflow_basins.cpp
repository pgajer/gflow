#include <unordered_map>
#include <queue>
#include <algorithm>      // For std::find()

#include <filesystem>     // For DEBUGGING only !!!
#include "cpp_utils.hpp"

#include "reachability_map.hpp"
#include "set_wgraph.hpp"
#include "basin.hpp"
#include "amagelo.hpp"
#include "error_utils.h"   // For REPORT_ERROR()

/**
 * @brief Writes neighborhood extrema to files for debugging
 *
 * @param out_dir Output directory path
 * @param nbr_lmin Vector of local minima vertex indices
 * @param nbr_lmax Vector of local maxima vertex indices
 * @param offset Value to add to vertex indices (0 for 0-based, 1 for 1-based R indexing)
 * @param prefix Optional prefix for output filenames
 */
void write_nbr_extrema(
    const std::string& out_dir,
    const std::vector<size_t>& nbr_lmin,
    const std::vector<size_t>& nbr_lmax,
    size_t offset,
    const std::string& prefix
	) {
    // Ensure output directory exists
    if (!std::filesystem::exists(out_dir)) {
        if (!std::filesystem::create_directories(out_dir)) {
            REPORT_ERROR("ERROR: Failed to create output directory: %s\n", out_dir.c_str());
            return;
        }
    }

    // Generate filenames
    std::string lmin_file = out_dir + "/" + prefix + "nbr_lmin.txt";
    std::string lmax_file = out_dir + "/" + prefix + "nbr_lmax.txt";

    // Write local minima vertices
    {
        std::ofstream lmin_out(lmin_file);
        if (!lmin_out) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmin_file.c_str());
            return;
        }

        lmin_out << "# Neighborhood Local Minima\n";
        lmin_out << "# Total count: " << nbr_lmin.size() << "\n";

        for (const auto& vertex : nbr_lmin) {
            lmin_out << (vertex + offset) << "\n";
        }

        lmin_out.close();
        Rprintf("Written neighborhood local minima to: %s\n", lmin_file.c_str());
    }

    // Write local maxima vertices
    {
        std::ofstream lmax_out(lmax_file);
        if (!lmax_out) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmax_file.c_str());
            return;
        }

        lmax_out << "# Neighborhood Local Maxima\n";
        lmax_out << "# Total count: " << nbr_lmax.size() << "\n";

        for (const auto& vertex : nbr_lmax) {
            lmax_out << (vertex + offset) << "\n";
        }

        lmax_out.close();
        Rprintf("Written neighborhood local maxima to: %s\n", lmax_file.c_str());
    }
}


/**
 * @param y                Vector of function values at each vertex
 * @param min_basin_size  The size of the neighborhood where the maximum (or minium) condition has to hold.
 * @param min_path_size   The minimal number of vertices of a shortest path for which we are going to estimate the monotonicity of y over that path
 * @param detect_maxima   If true, detect local maxima; if false, detect local minima
 *
 * @return basin_t structures. If 'vertex' is not a local extremum, 'basin' component of basin_t will be an empty vector
 */
struct geodesic_t {
	std::vector<size_t> vertices;  ///< vertices of the path
	std::vector<double> distances; ///< distance of each vertex to the reference vertex
	std::vector<double> y;         ///< y values at vertices
};

basin_t set_wgraph_t::find_gflow_basin(
    size_t vertex,
    const std::vector<double>& y,
    size_t min_basin_size,
    size_t min_path_size, // the minimum size of the shortest path for which non-linear regression condition for basin expansion is going to be applied
	double q_edge_thld,
	bool detect_maxima
    ) const {

	// amagelo() model parameter values
	// maybe in the future these may be passed via function arguments or a small struct, to avoid hard‐coding them deep inside this search routine.
	size_t grid_size = 100;
	double min_bw_factor = 0.05;
	double max_bw_factor = 0.5;
	size_t n_bws = 20;
	bool use_global_bw_grid = true;
	bool with_bw_predictions = false;
	bool log_grid = true;
	size_t domain_min_size = 4;
	size_t kernel_type = 7;
	double dist_normalization_factor = 1.1;
	size_t n_cleveland_iterations = 1;
	double blending_coef = 0.0;
	bool use_linear_blending = false;
	double precision = 1e-4;
	double small_depth_threshold = 0.05;
    double depth_similarity_tol  = 0.0001;
	bool verbose = false;

	// Initializing output basin
	basin_t basin;
	// Record the seed (reference) vertex
	basin.reachability_map.ref_vertex = vertex;
	basin.reachability_map.distances[vertex]    = 0.0;
	basin.reachability_map.predecessors[vertex] = INVALID_VERTEX;
	// Include the seed in the basin’s vertex list
	basin.reachability_map.sorted_vertices.push_back({
			basin.reachability_map.ref_vertex,       // the seed
			0.0                                      // zero distance to itself
		});
	// Store its value and extremum‐type
	basin.value      = y[vertex];
	basin.is_maximum = detect_maxima;


	// Initialize priority queue for Dijkstra's algorithm
	// Using min-heap with pairs of (distance, grid_vertex)
    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<>> pq;

	// Implement a bounded Dijkstra's algorithm
    size_t n_vertices = adjacency_list.size();
    std::vector<double> dist(n_vertices, INFINITY);
    std::vector<size_t> prev(n_vertices, INVALID_VERTEX);  // Using INVALID_VERTEX as sentinel
    dist[vertex] = 0;

    pq.push({0, vertex});

    while (!pq.empty()) {
		auto [d, u] = pq.top();
		pq.pop();

        if (d > dist[u]) continue;  // Skip if we've found a better path already

		// Explore neighbors - this allows alternative paths through u to be explored even when monotonicity condition at u is not met
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double w = edge.weight;
			double new_dist = dist[u] + w;
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                pq.push({dist[v], v});
            }
		}

		if (u != vertex) {  // Don't include the ref vertex itself

			// Calculate delta_y with appropriate sign
            double delta_y = y[u] - y[vertex];

			// Check extremum condition based on detect_maxima flag
            bool condition_met = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (condition_met) {

				// Determine the number of elements shortst_path_len in the
				// shortest path - should we enhance pq to be a pq of a struct
				// that holds this an possibly the whole path info of each
				// vertex or we simply reconstruct the path on the fly?

				// Reconstruct the shortest path from u to 'vertex'
				geodesic_t geodesic;
				{
					std::vector<size_t> path;
					for (size_t curr = u; curr != INVALID_VERTEX; curr = prev[curr]) {
						path.push_back(curr);
						if (curr == vertex) break;
					}
					std::reverse(path.begin(), path.end());

					size_t n_path_vertices = path.size();
					geodesic.y.resize(n_path_vertices);
					geodesic.distances.resize(n_path_vertices);
					for (size_t i = 0; i < n_path_vertices; ++i) {
						geodesic.y[i] = y[path[i]];
						geodesic.distances[i] = dist[path[i]];
					}
					geodesic.vertices = std::move(path);

					double max_edge_weight = 0;
					for (size_t i = 1; i < n_path_vertices; ++i) {
						double edge_weight = geodesic.distances[i] - geodesic.distances[i-1];
						if (edge_weight > max_edge_weight) {
							max_edge_weight = edge_weight;
						}
					}

					if (max_edge_weight > q_edge_thld) {
						continue;
					}
				}

				if (geodesic.vertices.size() >= min_path_size) {
					// Fit a non-linear regression model to (distances, y) components of 'geodesic'

					amagelo_t amagelo_fit = amagelo(
						geodesic.distances,
						geodesic.y,
						grid_size,
						min_bw_factor,
						max_bw_factor,
						n_bws,
						use_global_bw_grid,
						with_bw_predictions,
						log_grid,
						domain_min_size,
						kernel_type,
						dist_normalization_factor,
						n_cleveland_iterations,
						blending_coef,
						use_linear_blending,
						precision,
						small_depth_threshold,
						depth_similarity_tol,
						verbose
						);

					// Check monotonicity of amagelo_fit.predictions from that model
					// Note that amagelo sorts the input data in the increasing order of x-value, so (distances, y) will be sorted so that distances_sorted is in the increasing order

					bool is_monotonic = true;
					if (detect_maxima) {
						// the orientation of the (distances_sorted, y_sorted)
						// data is in the increasing distance from the local
						// maximum vertex, so we want to test if the differences
						// of all consecutive predictions of the model are
						// negative
						for (size_t i = 1; i < amagelo_fit.predictions.size(); ++i) {
							if (amagelo_fit.predictions[i] > amagelo_fit.predictions[i-1]) {
								is_monotonic = false;
								break;
							}
						}
					} else {
						for (size_t i = 1; i < amagelo_fit.predictions.size(); ++i) {
							if (amagelo_fit.predictions[i] < amagelo_fit.predictions[i-1]) {
								is_monotonic = false;
								break;
							}
						}
					}

					// If distance_predictions is stricly monotonically decreasing. when detect_maxima = true, we record u and explore its neighbors.
					// Otherwise, we not add u to the basin and we do not explore its neighbors

                    if (is_monotonic) {
						// Store vertex with its distance
						basin.reachability_map.sorted_vertices.push_back({u, d}); // the shortest path between u and 'vertex' has at least min_path_size vertices and y is monotonic along that path so we add it to the reachability map
						basin.reachability_map.distances[u]    = d;
						basin.reachability_map.predecessors[u] = prev[u];

                    } else {
                        continue; // we are not going to explore the neighbors of u as we already know that y cannot be monotonic along a longer geodesic containing u
                    }
                } else {
					// Store vertex with its distance
					basin.reachability_map.sorted_vertices.push_back({u, d}); // the shortest path between u and vertex has less than min_path_size vertices, so u is small hop distance from 'vertex', although the shortest path to u may have long edges !!! Explore it when testing
					basin.reachability_map.distances[u]    = d;
					basin.reachability_map.predecessors[u] = prev[u];
                }

            } else { // condition_met is not met
                continue; // we want to continuous as the basin of 'vertex' may be of irregular shape and so we need to allow Dijkstra to explore other paths
            }
        }
    }

	if (basin.reachability_map.sorted_vertices.size() < min_basin_size) {
		basin.reachability_map.sorted_vertices.clear();
		return basin;
	}

	// Sort vertices by distance in descending order
    std::sort(basin.reachability_map.sorted_vertices.begin(), basin.reachability_map.sorted_vertices.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    return basin;
}

/**
 * @brief Detects local minima and maxima of a response variable together with their basins
 *
 * @param y               Vector of function values at each vertex
 * @param min_basin_size  The size of the neighborhood where the maximum (or minium) condition has to hold.
 * @param min_path_size   The minimal number of vertices of a shortest path for which we are going to estimate the monotonicity of y over that path
 *
 * @return A pair <lmin basins, lmax basins>
 */
std::pair<std::vector<basin_t>, std::vector<basin_t>> set_wgraph_t::find_gflow_basins(
    const std::vector<double>& y,
	size_t min_basin_size,
    size_t min_path_size,
	double q_edge_thld
    ) const {

	// Phase 1:  Process each vertex as a potential extremum; checking local extremum conditions only within the set of each vertex's neighbors
	std::vector<size_t> lmin_candidates; // indices of vertices that satisfy local min condition over the set of its neighbors
    std::vector<size_t> lmax_candidates; // indices of vertices that satisfy local max condition over the set of its neighbors

    for (size_t vertex = 0; vertex < adjacency_list.size(); ++vertex) {

		// Checking if vertex is a local minimum within the set of its neighbors
        bool lmin_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
            if (y[edge.vertex] <= y[vertex]) { // vertex cannot be a lmin of y
				lmin_condition_met = false;
                break;
            }
        }

		if(lmin_condition_met) {
			lmin_candidates.push_back(vertex);
			continue; // there is no point of checking if it is a local maximum
		}

		// Checking if vertex is a local minimum within the set of its neighbors
		bool lmax_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
			if (y[edge.vertex] >= y[vertex]) { // vertex cannot be a lmax of y
				lmax_condition_met = false;
                break;
            }
        }

		if(lmax_condition_met) {
			lmax_candidates.push_back(vertex);
		}
    }

	std::vector<basin_t> lmin_basins;
	std::vector<basin_t> lmax_basins;

	// Phase 2: For each vertex of lmin_candidates check if it is a local
	// minimum within the set of core_basin_size neighbors
	bool detect_maxima = false;
	for (const auto& vertex : lmin_candidates) {
		auto basin =  find_gflow_basin(
			vertex,
			y,
			min_basin_size,
			min_path_size,
			q_edge_thld,
			detect_maxima);

		if (basin.reachability_map.sorted_vertices.size()) {
			lmin_basins.push_back(basin);
		}
	}


	// Phase 3: For each vertex of lmax_candidates check if it is a local
	// maximum within the set of core_basin_size neighbors
	detect_maxima = true;
	for (const auto& vertex : lmax_candidates) {
		auto basin =  find_gflow_basin(
			vertex,
			y,
			min_basin_size,
			min_path_size,
			q_edge_thld,
			detect_maxima);

		if (basin.reachability_map.sorted_vertices.size()) {
			lmax_basins.push_back(basin);
		}
	}

	return std::make_pair(lmin_basins, lmax_basins);
}


/**
 * @brief Finds a local extremum and its basin of size min_basin_size by bounded Dijkstra.
 *
 * @param vertex           Reference vertex to test as local extremum
 * @param y                Vector of function values at each vertex
 * @param min_basin_size   Number of vertices (including seed) to include in basin
 * @param detect_maxima    If true, detect local maximum; if false, local minimum
 * @return basin_t         If vertex is not an extremum or basin too small, sorted_vertices is empty
 */
basin_t set_wgraph_t::find_local_extremum(
    size_t vertex,
    const std::vector<double>& y,
    size_t min_basin_size,
    bool detect_maxima
	) const {
    basin_t basin;
    basin.value      = y[vertex];
    basin.is_maximum = detect_maxima;
    basin.reachability_map.ref_vertex = vertex;

#if 0
	// 1. Check immediate extremum condition
    for (const auto& edge : adjacency_list[vertex]) {
        double neigh_val = y[edge.vertex];
        if (detect_maxima ? (neigh_val > y[vertex]) : (neigh_val < y[vertex])) {
			// Not a local extremum
            return basin;  // sorted_vertices empty
        }
    }
#endif

	// 2. Initialize Dijkstra structures
    size_t n = adjacency_list.size();
    std::vector<double> dist(n, INFINITY);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<queue_entry>> pq;

    dist[vertex] = 0.0;
    pq.push({0.0, vertex});

	// Seed the basin with the reference vertex
    basin.reachability_map.distances[vertex]    = 0.0;
    basin.reachability_map.predecessors[vertex] = INVALID_VERTEX;
    basin.reachability_map.sorted_vertices.push_back({vertex, 0.0});

	// 3. Perform bounded Dijkstra until we have min_basin_size vertices
    size_t count = 1;  // already have the seed
    while (!pq.empty() && count < min_basin_size) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;  // stale entry

		// Explore neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double nd = d + edge.weight;
            if (nd < dist[v]) {
                dist[v] = nd;
                prev[v] = u;
                pq.push({nd, v});
            }
        }

        if (u != vertex) {
			// Add to basin
            basin.reachability_map.sorted_vertices.push_back({u, d});
            basin.reachability_map.distances[u]    = d;
            basin.reachability_map.predecessors[u] = prev[u];
            ++count;
        }
    }

	// 4. If we didn't reach the required size, clear and return empty
    if (count < min_basin_size) {
        basin.reachability_map.sorted_vertices.clear();
        return basin;
    }

	// 5. Sort vertices by descending distance
    std::sort(
        basin.reachability_map.sorted_vertices.begin(),
        basin.reachability_map.sorted_vertices.end(),
        [](const vertex_info_t& a, const vertex_info_t& b) {
            return a.distance > b.distance;
        }
		);

    return basin;
}

/**
 * @brief  Find all local extrema basins of a graph‐valued function.
 *
 * For each vertex, test whether it is a strict local minimum or maximum
 * among its immediate neighbors.  Then for each candidate, grow a basin
 * of size `min_basin_size` using a bounded Dijkstra search
 * (via `find_local_extremum`).  Basins that fail to reach the required
 * size return empty and are discarded.
 *
 * @param[in] y               Vector of function values at each graph vertex.
 * @param[in] min_basin_size  Minimum number of vertices (including the seed)
 *                            required for a basin to be accepted.
 * @return A pair of vectors:
 *   - first:  basins around local minima
 *   - second: basins around local maxima
 * Each `basin_t` in the output has its `reachability_map.sorted_vertices`
 * populated only if the basin grew to at least `min_basin_size` vertices.
 */
std::pair<std::vector<basin_t>, std::vector<basin_t>> set_wgraph_t::find_local_extrema(
    const std::vector<double>& y,
	size_t min_basin_size
    ) const {

	// Phase 1:  Process each vertex as a potential extremum; checking local extremum conditions only within the set of each vertex's neighbors
	std::vector<size_t> lmin_candidates; // indices of vertices that satisfy local min condition over the set of its neighbors
    std::vector<size_t> lmax_candidates; // indices of vertices that satisfy local max condition over the set of its neighbors

    for (size_t vertex = 0; vertex < adjacency_list.size(); ++vertex) {

		// Checking if vertex is a local minimum within the set of its neighbors
        bool lmin_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
            if (y[edge.vertex] <= y[vertex]) { // vertex cannot be a lmin of y
				lmin_condition_met = false;
                break;
            }
        }

		if(lmin_condition_met) {
			lmin_candidates.push_back(vertex);
			continue; // there is no point of checking if it is a local maximum
		}

		// Checking if vertex is a local minimum within the set of its neighbors
		bool lmax_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
			if (y[edge.vertex] >= y[vertex]) { // vertex cannot be a lmax of y
				lmax_condition_met = false;
                break;
            }
        }

		if(lmax_condition_met) {
			lmax_candidates.push_back(vertex);
		}
    }

	std::vector<basin_t> lmin_basins;
	std::vector<basin_t> lmax_basins;

	// Phase 2: For each vertex of lmin_candidates check if it is a local
	// minimum within the set of core_basin_size neighbors
	bool detect_maxima = false;
	for (const auto& vertex : lmin_candidates) {
		auto basin = find_local_extremum(
			vertex,
			y,
			min_basin_size,
			detect_maxima);

		if (basin.reachability_map.sorted_vertices.size()) {
			lmin_basins.push_back(basin);
		}
	}


	// Phase 3: For each vertex of lmax_candidates check if it is a local
	// maximum within the set of core_basin_size neighbors
	detect_maxima = true;
	for (const auto& vertex : lmax_candidates) {
		auto basin = find_local_extremum(
			vertex,
			y,
			min_basin_size,
			detect_maxima);

		if (basin.reachability_map.sorted_vertices.size()) {
			lmax_basins.push_back(basin);
		}
	}

	return std::make_pair(lmin_basins, lmax_basins);

}


/**
 * @brief Finds vertices that are local extrema with respect to their neighbors
 *
 * For each vertex in the graph, this function determines if its value in the given
 * array is either a local minimum or a local maximum compared to all its adjacent
 * vertices (neighbors).
 *
 * A vertex is considered a local minimum if its value is strictly less than the values
 * of all its neighbors. Similarly, a vertex is a local maximum if its value is strictly
 * greater than the values of all its neighbors.
 *
 * Note: A vertex cannot be both a local minimum and a local maximum simultaneously.
 *
 * @param[in] y A vector of double values, where y[i] represents the value at vertex i
 *
 * @return A pair of vectors containing:
 *         - first: Indices of vertices that are local minima
 *         - second: Indices of vertices that are local maxima
 *
 * @pre The size of vector y must be at least equal to the number of vertices in the graph
 *
 * @complexity O(E), where E is the total number of edges in the graph
 */
std::pair<std::vector<size_t>, std::vector<size_t>> set_wgraph_t::find_nbr_extrema(
    const std::vector<double>& y
    ) const {

	//  Check for each vertex if y is a local extremum over the completion of the set of its neighbors
	std::vector<size_t> nbr_lmin; // indices of vertices that satisfy local min condition over the set of its neighbors
    std::vector<size_t> nbr_lmax; // indices of vertices that satisfy local max condition over the set of its neighbors

    for (size_t vertex = 0; vertex < adjacency_list.size(); ++vertex) {

		// Checking if vertex is a local minimum within the set of its neighbors
        bool lmin_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
            if (y[edge.vertex] <= y[vertex]) { // vertex cannot be a lmin of y
				lmin_condition_met = false;
                break;
            }
        }

		if(lmin_condition_met) {
			nbr_lmin.push_back(vertex);
			continue; // there is no point of checking if it is a local maximum
		}

		// Checking if vertex is a local maximum within the set of its neighbors
		bool lmax_condition_met = true;
		for (const auto& edge : adjacency_list[vertex]) {
			if (y[edge.vertex] >= y[vertex]) { // vertex cannot be a lmax of y
				lmax_condition_met = false;
                break;
            }
        }

		if(lmax_condition_met) {
			nbr_lmax.push_back(vertex);
		}
    }

	return std::make_pair(nbr_lmin, nbr_lmax);
}


/**
 * @brief Segment the graph via an edge-weighted watershed using local extrema from \c find_local_extrema().
 *
 * First identifies significant local minima of size >= min_basin_size by calling
 * \c find_local_extrema.  These minima become the markers for a priority-flood watershed.
 * The barrier on each edge (u,v) is
 *   B(u,v) = max(|yhat[u] - yhat[v]|, q_uv),
 * where yhat is a min-max normalization of y to [0,1] and q_uv is the empirical
 * quantile of the edge weight w_{uv} among all edges.
 *
 * @param y               Input values at vertices (size = n_vertices).
 * @param min_basin_size  Minimum basin size to accept a local minimum marker.
 * @return vector<int>    1-based basin ID for each vertex, or 0 if unassigned.
 */
std::vector<int> set_wgraph_t::watershed_edge_weighted(
    const std::vector<double>& y,
    size_t min_basin_size
	) const {
    size_t n = adjacency_list.size();
	// 1. Normalize y -> yhat in [0,1]
    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    double range = (y_max > y_min ? y_max - y_min : 1.0);
    std::vector<double> yhat(n);
    for (size_t i = 0; i < n; ++i) {
        yhat[i] = (y[i] - y_min) / range;
    }

	// 2. Gather & sort all edge weights for quantile computation
    std::vector<double> all_w;
    all_w.reserve(n * 4); // rough estimate
    for (size_t u = 0; u < n; ++u) {
        for (auto const &e : adjacency_list[u]) {
            all_w.push_back(e.weight);
        }
    }
    std::sort(all_w.begin(), all_w.end());
    size_t m = all_w.size();

	// 3. Build per-vertex barrier map: barrier_map[u][v] = B(u,v)
    std::vector<std::unordered_map<size_t,double>> barrier_map(n);
    for (size_t u = 0; u < n; ++u) {
        for (auto const &e : adjacency_list[u]) {
            size_t v = e.vertex;
            double diff = std::abs(yhat[u] - yhat[v]);
            auto it = std::lower_bound(all_w.begin(), all_w.end(), e.weight);
            size_t idx = std::distance(all_w.begin(), it);
            double q = (m > 1 ? double(idx) / (m - 1) : 0.0);
            barrier_map[u][v] = std::max(diff, q);
        }
    }

	// 4. Find true local minima basins (size >= min_basin_size)
    auto bas_pair = find_local_extrema(y, min_basin_size);
    const auto &lmin_basins = bas_pair.first;

	// 5. Initialize labels (0 = unassigned)
    std::vector<int> label(n, 0);
    int next_label = 1;
    for (auto const &b : lmin_basins) {
        size_t u = b.reachability_map.ref_vertex;
        label[u] = next_label++;
    }

	// 6. Priority-flood from each minima marker
    struct Entry { double pr; size_t v; int reg; };
    struct Cmp { bool operator()(Entry const &a, Entry const &b) const {
        return a.pr > b.pr;
    }};
    std::priority_queue<Entry, std::vector<Entry>, Cmp> pq;

	// Seed neighbors of each minima
    for (auto const &b : lmin_basins) {
        size_t u = b.reachability_map.ref_vertex;
        int reg = label[u];
        for (auto const &e : adjacency_list[u]) {
            size_t v = e.vertex;
            if (label[v] == 0) {
                double pr = barrier_map[u][v];
                pq.push({pr, v, reg});
            }
        }
    }

	// Flood assign
    while (!pq.empty()) {
        auto cur = pq.top(); pq.pop();
        size_t v = cur.v; int reg = cur.reg; double pr = cur.pr;
        if (label[v] == 0) {
            label[v] = reg;
            for (auto const &e : adjacency_list[v]) {
                size_t w = e.vertex;
                if (label[w] == 0) {
                    double pr2 = std::max(pr, barrier_map[v][w]);
                    pq.push({pr2, w, reg});
                }
            }
        }
    }

    return label;
}

/**
 * @brief Find the basin of attraction for a local extremum (maximum or minimum)
 *
 * @details This function constructs the basin of attraction around a local extremum by
 *          exploring the graph while maintaining monotonicity constraints. It identifies
 *          the prominence/depth of the extremum and all boundary vertices of the basin.
 *
 * The basin construction follows these key principles:
 * 1. For a maximum: function values must strictly decrease as we move away from the extremum
 * 2. For a minimum: function values must strictly increase as we move away from the extremum
 * 3. The basin includes all vertices reachable from the extremum while maintaining monotonicity
 * 4. The prominence (depth) is the smallest vertical distance from the extremum to the point
 *    of monotonicity violation
 *
 * The function identifies two types of boundary vertices:
 * 1. Monotonicity violation boundaries: vertices where the monotonicity condition fails
 * 2. Graph extent boundaries: vertices within the basin that have neighbors outside the basin
 *
 * @param vertex The vertex to treat as a local extremum
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, detect a maximum basin; if false, detect a minimum basin
 * @param edge_length_thld Edge length threshold for basin construction
 *
 * @return A basin_t structure containing:
 *         - value: Function value at the extremum
 *         - is_maximum: Whether this is a maximum (true) or minimum (false)
 *         - reachability_map: Information to reconstruct paths within the basin
 *         - depth: Prominence of the extremum (minimal vertical distance from extremum to last valid point)
 *         - depth_vertex: Vertex at which the depth/prominence is realized
 *         - boundary_vertices_map: Map of distances to boundary vertices, sorted by increasing distance
 *
 * @details The edge_length_thld parameter imposes a geometric constraint on
 *          basin growth by excluding edges whose Riemannian length exceeds
 *          the specified threshold. This prevents the "basin jumping" pathology
 *          where long edges cause trajectories to skip over intermediate critical
 *          points or violate local monotonicity structure.
 *
 *          During basin exploration, an edge [u,v] is only traversed if both:
 *          1. The monotonicity condition is satisfied (function values increase
 *             for minima or decrease for maxima), AND
 *          2. The edge length satisfies: length([u,v]) <= edge_length_thld
 *
 *          When an edge fails the length constraint, vertex u is treated as a
 *          boundary vertex of the basin, exactly as if the monotonicity condition
 *          had been violated.
 *
 * @section setting Setting the Threshold
 *
 *          The threshold is typically set as a quantile of the empirical edge
 *          length distribution. Common choices:
 *
 *          - **50th percentile (median)**: Conservative, excludes longest half of edges
 *          - **75th percentile**: Moderate, balances completeness with pathology prevention
 *          - **90th percentile**: Liberal, only excludes the longest 10% of edges
 *          - **std::numeric_limits<double>::infinity()**: No constraint (original behavior)
 *
 *          Example calculation in R:
 *          @code{.r}
 *          ## Compute edge length quantile from graph
 *          edge_lengths <- get_edge_lengths(graph)
 *          edge_length_thld <- quantile(edge_lengths, probs = 0.75)
 *          @endcode
 *
 * @section behavior Behavioral Notes
 *
 *          - **No constraint**: Pass std::numeric_limits<double>::infinity() to
 *            disable edge length filtering entirely, recovering original basin
 *            construction behavior
 *
 *          - **Non-existent edges**: The get_edge_weight() function returns
 *            INFINITY for non-existent edges, which automatically fails the
 *            threshold test (correct behavior)
 *
 *          - **Zero threshold**: Passing 0.0 would restrict basins to single
 *            vertices (typically not useful)
 *
 *          - **Negative threshold**: Not recommended; interpreted as "no edges
 *            satisfy constraint"
 *
 * @section motivation Geometric Motivation
 *
 *          In discrete geometric data analysis, long edges often indicate
 *          regions where the graph structure poorly represents the underlying
 *          geometry. A long edge [v_i, v_j] with y(v_j) > y(v_i) may appear
 *          to satisfy gradient ascent locally, yet span a region where the
 *          response function attains intermediate extrema. By excluding such
 *          edges, we enforce that gradient flows respect the local structure
 *          of the response surface at the resolution captured by shorter edges.
 *
 * @warning This constraint is applied uniformly across the entire graph. In
 *          regions with naturally varying edge density, a global quantile may
 *          be too restrictive in sparse regions or too permissive in dense regions.
 *          Consider the spatial distribution of edge lengths when selecting the
 *          threshold.
 *
 * @see compute_geodesic_basin() Wrapper that passes threshold to basin computation
 *
 * @note This function assumes the caller has already verified that 'vertex' is a local extremum.
 *       It does not perform the initial extremum check.
 */
basin_t set_wgraph_t::find_local_extremum_geodesic_basin(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima,
	double edge_length_thld
    ) const {

	// Initialize basin structure with extremum properties
    basin_t basin;
    basin.value      = y[vertex];                 // Function value at extremum
    basin.is_maximum = detect_maxima;             // Type of extremum (max or min)
    basin.reachability_map.ref_vertex = vertex;   // Reference vertex for paths
    basin.min_monotonicity_span = INFINITY;       // needed for making sure min_monotonicity_span is set when the reason for basin termination is that we got to an end of the graph
    basin.min_span_vertex       = INVALID_VERTEX; // if no violations were found during exploration we never find basin.min_span_vertex
	// Seed the basin with the reference vertex (extremum)
    basin.reachability_map.distances[vertex]    = 0.0;
    basin.reachability_map.predecessors[vertex] = INVALID_VERTEX;
    basin.reachability_map.sorted_vertices.push_back({vertex, 0.0});

	// Initialize Dijkstra's algorithm data structures
    size_t n = adjacency_list.size();             // Total number of vertices in graph
    std::vector<double> dist(n, INFINITY);        // Shortest path distances
    std::vector<size_t> prev(n, INVALID_VERTEX);  // Predecessor vertices in shortest paths
    std::vector<bool> is_in_pq(n, false);         // Track which vertices are in priority queue

    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<queue_entry>> pq;

	// Initialize starting vertex
    dist[vertex] = 0.0;
    pq.push({0.0, vertex});
    is_in_pq[vertex] = true;  // Mark extremum as in queue

	// Perform monotonicity-constrained Dijkstra's algorithm
    while (!pq.empty()) {
		// Get closest vertex from priority queue
        auto [d, u] = pq.top(); pq.pop();
        is_in_pq[u] = false;  // No longer in queue

		// Skip if we've already found a shorter path to u
        if (d > dist[u]) continue;

		// For all vertices except the extremum itself
        if (u != vertex) {
			// Check if monotonicity condition is maintained
            double delta_y = y[u] - y[prev[u]];

			// For maxima: values must decrease (delta_y < 0)
			// For minima: values must increase (delta_y > 0)
			bool monotonicity_ok = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);
			double edge_length = get_edge_weight(prev[u], u);
			bool edge_length_ok = (edge_length <= edge_length_thld);
			bool condition_met = monotonicity_ok && edge_length_ok;

            if (!condition_met) {
				// *** MONOTONICITY OR EDGE LENGTH VIOLATION BOUNDARY ***
				// Boundary reached when either:
				// 1. Monotonicity is violated, OR
				// 2. Edge length exceeds threshold

				// Could add diagnostic info about which constraint failed
				// if (!monotonicity_ok) { /* monotonicity violation */ }
				// if (!edge_length_ok) { /* edge length violation */ }

				// When monotonicity or edge length is violated, calculate monotonicity span

				// For maxima: span = y[vertex] - y[prev[u]]
				// For minima: span = y[prev[u]] - y[vertex]
                double curr_span = detect_maxima ? (y[vertex] - y[prev[u]]) : (y[prev[u]] - y[vertex]);

				// Store the monotonicity span for this boundary vertex (prev[u])
                basin.boundary_monotonicity_spans_map[prev[u]] = curr_span;

				// Update minimum span if this violation gives a smaller valid span
                constexpr double eps = 1e-12;  // Small epsilon to handle floating-point precision
                if (curr_span > eps && curr_span < basin.min_monotonicity_span) {
                    basin.min_monotonicity_span = curr_span;
                    basin.min_span_vertex = prev[u];  // Last valid vertex before violation
                }

				// Add the last valid vertex to boundary map
				// (the vertex before monotonicity violation)
                basin.boundary_vertices_map[dist[prev[u]]] = prev[u];

				// Skip exploring neighbors of violating vertex
                continue;
            }

			// If monotonicity is maintained, add vertex to basin
            basin.reachability_map.sorted_vertices.push_back({u, d});
            basin.reachability_map.distances[u] = d;
            basin.reachability_map.predecessors[u] = prev[u];
        }

		// Explore neighbors of acceptable vertices (including extremum)
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double nd = d + edge.weight;  // New distance to v

			// If we found a shorter path to v
            if (nd < dist[v]) {
                dist[v] = nd;         // Update distance
                prev[v] = u;          // Update predecessor
                pq.push({nd, v});     // Add to priority queue
                is_in_pq[v] = true;   // Mark as in queue
            }
        }

		// *** GRAPH EXTENT BOUNDARY ***
        // Check if this vertex is a boundary due to having neighbors outside the basin or its degree is 1
		// We need this to correctly compute vertex's monotonicity span
		// Even though at an end of a graph we are not violating the monotonicity condition, usually it will be a local extremum of y and so if 'vertex' is a spurious local extremum the end vertex can be used to cut or fill it the values of y within this spurios basin
        bool is_boundary = false;
		if (adjacency_list[u].size() == 1) {
			is_boundary = true;
		}  else {
			for (const auto& edge : adjacency_list[u]) {
				size_t v = edge.vertex;

				// A vertex is a boundary if it has a neighbor that:
				// 1. Is not currently in the basin AND
				// 2. Is not waiting in the priority queue
				//
				// This identifies vertices that are truly at the edge of the basin,
				// where the basin cannot extend further due to graph structure
				if (basin.reachability_map.distances.find(v) == basin.reachability_map.distances.end() &&
					!is_in_pq[v]) {
					is_boundary = true;
					break;
				}
			}
		}

        // If this is a boundary vertex, add it to the boundary map
        if (is_boundary) {
            basin.boundary_vertices_map[d] = u;
			basin.boundary_monotonicity_spans_map[u] = std::abs(y[vertex] - y[u]);
		}
    }

	// Sort vertices by descending distance for easier downstream analysis
    std::sort(
        basin.reachability_map.sorted_vertices.begin(),
        basin.reachability_map.sorted_vertices.end(),
        [](const vertex_info_t& a, const vertex_info_t& b) {
            return a.distance > b.distance;
        }
        );

    return basin;
}

/**
 * @brief Find the basin of attraction for a local extremum using breadth-first search
 *
 * @details This function constructs the basin of attraction around a local extremum by
 *          exploring the graph using breadth-first search (BFS) rather than Dijkstra's algorithm.
 *          Unlike the geodesic basin approach which only considers shortest paths, this method
 *          explores all possible paths from the extremum while maintaining monotonicity constraints.
 *
 *          The key differences from find_local_extremum_geodesic_basin():
 *          1. Uses BFS instead of Dijkstra's algorithm (treating all edges as unit weight).
 *          2. Considers all paths that maintain monotonicity, not just geodesic paths.
 *          3. Typically results in larger, more complete basins that better reflect the shape of the response y.
 *
 * The basin construction follows these key principles:
 * 1. For a maximum: function values must strictly decrease as we move away from the extremum
 * 2. For a minimum: function values must strictly increase as we move away from the extremum
 * 3. The basin includes all vertices reachable from the extremum while maintaining monotonicity
 * 4. The basin boundary consists of vertices where monotonicity would be violated or at the graph boundary
 *
 * @param vertex The vertex to treat as a local extremum
 * @param y Vector of function values at each vertex
 * @param detect_maxima If true, detect a maximum basin; if false, detect a minimum basin
 *
 * @return A basin_t structure containing:
 *         - value: Function value at the extremum
 *         - is_maximum: Whether this is a maximum (true) or minimum (false)
 *         - extremum_vertex: The vertex index of the extremum
 *         - reachability_map: Information to reconstruct paths within the basin
 *         - min_monotonicity_span: Minimum vertical distance from extremum to any boundary
 *         - min_span_vertex: Vertex at which the minimum monotonicity span is realized
 *         - boundary_vertices_map: Map of hop distances to boundary vertices
 *         - boundary_monotonicity_spans_map: Monotonicity spans at each boundary vertex
 *
 * @note This function produces more comprehensive basins that include all vertices
 *       reachable via any monotonic path, not just those on shortest paths. This is particularly
 *       useful for complex graphs with multiple potential paths between vertices.
 *
 * @see find_local_extremum_geodesic_basin - The geodesic (Dijkstra) version of this function
 * @see create_basin_cx - The basin complex construction that uses extremum basins
 */
basin_t set_wgraph_t::find_local_extremum_bfs_basin(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
    ) const {
    // Initialize basin structure with extremum properties
    basin_t basin;
    basin.value                 = y[vertex];      // Function value at extremum
    basin.is_maximum            = detect_maxima;  // Type of extremum (max or min)
    basin.reachability_map.ref_vertex = vertex;   // Reference vertex for paths
    basin.min_monotonicity_span = INFINITY;       // Initialize to infinity
    basin.min_span_vertex       = INVALID_VERTEX; // Invalid until a boundary is found
    basin.extremum_vertex       = vertex;         // Store the extremum vertex explicitly

    // Seed the basin with the reference vertex (extremum)
    basin.reachability_map.distances[vertex]    = 0.0;
    basin.reachability_map.predecessors[vertex] = INVALID_VERTEX;
    basin.reachability_map.sorted_vertices.push_back({vertex, 0.0});

    // Initialize BFS data structures
    size_t n = adjacency_list.size();              // Total number of vertices in graph
    std::vector<int> hop_distance(n, -1);          // BFS hop distances (-1 means not visited)
    std::vector<size_t> prev(n, INVALID_VERTEX);   // Predecessor vertices
    std::queue<size_t> q;                          // BFS queue
    std::vector<bool> in_queue(n, false);          // Track which vertices are in queue
    std::vector<bool> in_basin(n, false);          // Track which vertices are in basin

    // Initialize starting vertex
    hop_distance[vertex] = 0;
    q.push(vertex);
    in_queue[vertex] = true;
    in_basin[vertex] = true;

    // Perform breadth-first search with monotonicity constraints
    while (!q.empty()) {
        // Get next vertex from queue
        size_t u = q.front();
        q.pop();
        in_queue[u] = false;

        // For vertices other than the extremum
        if (u != vertex) {
            // Check if monotonicity condition is maintained along the path
            double delta_y = y[u] - y[prev[u]];

            // For maxima: values must decrease (delta_y < 0)
            // For minima: values must increase (delta_y > 0)
            bool condition_met = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!condition_met) {
                // *** MONOTONICITY VIOLATION BOUNDARY ***

                // Calculate monotonicity span from extremum to the last valid vertex
                // For maxima: span = y[vertex] - y[prev[u]]
                // For minima: span = y[prev[u]] - y[vertex]
                double curr_span = detect_maxima ? (y[vertex] - y[prev[u]]) : (y[prev[u]] - y[vertex]);

                // Store the monotonicity span for this boundary vertex
                basin.boundary_monotonicity_spans_map[prev[u]] = curr_span;

                // Update minimum span if this violation gives a smaller valid span
                constexpr double eps = 1e-12;  // Small epsilon for floating-point precision
                if (curr_span > eps && curr_span < basin.min_monotonicity_span) {
                    basin.min_monotonicity_span = curr_span;
                    basin.min_span_vertex = prev[u];
                }

                // Add to boundary vertices map (using hop distance instead of geodesic distance)
                basin.boundary_vertices_map[hop_distance[prev[u]]] = prev[u];

                // Skip exploring this violating vertex
                continue;
            }

            // If monotonicity is maintained, add vertex to basin
            double dist = static_cast<double>(hop_distance[u]); // Use hop count as distance
            basin.reachability_map.sorted_vertices.push_back({u, dist});
            basin.reachability_map.distances[u] = dist;
            basin.reachability_map.predecessors[u] = prev[u];
            in_basin[u] = true;
        }

        // Explore neighbors (all edge weights treated equally in BFS)
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;

            // If not visited or found a better path
            if (hop_distance[v] == -1) {
                hop_distance[v] = hop_distance[u] + 1;
                prev[v] = u;

                // Add to queue if not already there
                if (!in_queue[v]) {
                    q.push(v);
                    in_queue[v] = true;
                }
            }
        }
    } // END OF while (!q.empty())

	// Identify the boundary of the basin
	for (const auto& v_info : basin.reachability_map.sorted_vertices) {
		size_t u = v_info.vertex;

		// Check if vertex is a boundary due to having neighbors outside the basin
		bool is_boundary = false;

		// Leaf vertices (degree 1) are boundaries
		if (adjacency_list[u].size() == 1) {
			is_boundary = true;
		} else {
			// Check if any neighbors are outside the basin
			for (const auto& edge : adjacency_list[u]) {
				size_t v = edge.vertex;

				// A vertex is a boundary if it has a neighbor that is not in the basin
				if (basin.reachability_map.distances.find(v) == basin.reachability_map.distances.end()) {
					is_boundary = true;
					break;
				}
			}
		}

		// If this is a boundary vertex, add it to the boundary maps
		if (is_boundary) {
			double dist = basin.reachability_map.distances[u];
			basin.boundary_vertices_map[dist] = u;

			// Calculate monotonicity span for this boundary
			double span = detect_maxima ? (y[vertex] - y[u]) : (y[u] - y[vertex]);
			basin.boundary_monotonicity_spans_map[u] = span;
		}
	}

    // Sort vertices by descending distance for easier downstream analysis
    std::sort(
        basin.reachability_map.sorted_vertices.begin(),
        basin.reachability_map.sorted_vertices.end(),
        [](const vertex_info_t& a, const vertex_info_t& b) {
            return a.distance > b.distance;
        }
		);

    return basin;
}


/**
 * @brief Constructs a graph gradient flow basin complex from a scalar function on vertices
 *
 * @details This function identifies local extrema (minima and maxima) on the graph and
 *          constructs monotonic basins around each extremum. It then processes basin pairs
 *          that satisfy cancellation criteria, absorbing less significant basins into more
 *          significant ones based on relative monotonicity span. The function performs harmonic
 *          interpolation to repair function values across absorbed basins, eliminating
 *          topologically insignificant features.
 *
 * The algorithm proceeds through these key steps:
 * 1. Identify local extrema based on neighborhood comparison
 * 2. Grow monotonic basins around each extremum using modified Dijkstra
 * 3. Sort basins by their relative monotonicity span
 * 4. Process cancellation pairs (A, B) where:
 *    - A_basin.min_span_vertex = B_basin.extremum_vertex
 *    - B_basin.min_span_vertex = A_basin.extremum_vertex
 * 5. Find the closest "cousin" basin for each member of a cancellation pair
 * 6. Absorb the less significant basin into its cousin basin
 * 7. Update function values using harmonic repair
 * 8. Continue until all basins with rel_min_monotonicity_span < threshold are processed
 *
 * @param y The scalar function values defined at each vertex of the graph
 * @param rel_min_monotonicity_span_thld Threshold for relative monotonicity span below which
 *                             basins are considered for cancellation
 *
 * @return A basin_cx_t structure containing:
 *         - harmonic_predictions: Updated function values after repair
 *         - lmin_basins_map: Map of local minima to their basins
 *         - lmax_basins_map: Map of local maxima to their basins
 *         - sorted_basins: Basins sorted by increasing relative monotonicity span
 *         - cancellation_pairs: Pairs of basins that were identified for cancellation
 *         - repair_lextr: Map tracking which basin absorbed which during cancellation
 *
 * @see basin_t
 * @see find_local_extremum_geodesic_basin
 * @see find_nbr_extrema
 */
basin_cx_t set_wgraph_t::create_basin_cx(
    const std::vector<double>& y
    ) const {

	// 1. find nbr-local extrema
    auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

    basin_cx_t basin_cx;
    basin_cx.harmonic_predictions = y;  // Initialize with original values

    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    double y_range = y_max - y_min;

	// 2. for each nbr lmin of nbr_lmin find its basin
    bool detect_maxima = false;
    for (const auto& vertex : nbr_lmin) {
        auto basin = find_local_extremum_bfs_basin(
            vertex,
            y,
            detect_maxima
            );

		basin.extremum_hop_index = compute_extremum_hop_index(
			vertex,
			y,
			detect_maxima);


        if (basin.reachability_map.sorted_vertices.size()) {
			// Check if basin.min_monotonicity_span was computed correctly. Note
			// that basin.min_monotonicity_span was computed using only
			// monotonicity violation vertices, not the end of the graph
			// vertices that could be also local extrema, in which case we
			// should include them in the span calculation. Now, that we have
			// access to the local extrema information we can include also end
			// of graph vertices in the calculation of
			// basin.min_monotonicity_span if they are local maxima (for the
			// current vertex local min case)
			double min_span = INFINITY;
			size_t min_span_vertex = INVALID_VERTEX;
			double max_span = 0;
			size_t max_span_vertex = INVALID_VERTEX;
			for (const auto& [u, span] : basin.boundary_monotonicity_spans_map) {
				if (span < min_span) {
					min_span = span;
					min_span_vertex = u;
				}
				if (span > max_span) {
					max_span = span;
					max_span_vertex = u;
				}
			}

			basin.min_span_vertex           = min_span_vertex;
			basin.min_monotonicity_span     = min_span;
			basin.rel_min_monotonicity_span = (y_range > 0.0) ? min_span / y_range : 0.0;

			basin.max_span_vertex           = max_span_vertex;
			basin.max_monotonicity_span     = max_span;
			basin.rel_max_monotonicity_span = (y_range > 0.0) ? max_span / y_range : 0.0;

			basin.delta_rel_span = basin.rel_max_monotonicity_span - basin.rel_min_monotonicity_span;

			basin.rel_size = static_cast<double>(basin.reachability_map.sorted_vertices.size()) / y.size();

            basin.extremum_vertex = vertex;
            basin_cx.lmin_basins_map[vertex] = std::move(basin);
        }
    }

	// 3. for each nbr lmax of nbr_lmax find its basin
    detect_maxima = true;
    for (const auto& vertex : nbr_lmax) {
        auto basin = find_local_extremum_bfs_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
			// See comments in the local min block for this part
			double min_span = INFINITY;
			size_t min_span_vertex = INVALID_VERTEX;
			double max_span = 0;
			size_t max_span_vertex = INVALID_VERTEX;
			for (const auto& [u, span] : basin.boundary_monotonicity_spans_map) {
				if (span < min_span) {
					min_span = span;
					min_span_vertex = u;
				}
				if (span > max_span) {
					max_span = span;
					max_span_vertex = u;
				}
			}

			basin.min_span_vertex           = min_span_vertex;
			basin.min_monotonicity_span     = min_span;
			basin.rel_min_monotonicity_span = (y_range > 0.0) ? min_span / y_range : 0.0;

			basin.max_span_vertex           = max_span_vertex;
			basin.max_monotonicity_span     = max_span;
			basin.rel_max_monotonicity_span = (y_range > 0.0) ? max_span / y_range : 0.0;

			basin.delta_rel_span = basin.rel_max_monotonicity_span - basin.rel_min_monotonicity_span;

			basin.rel_size = static_cast<double>(basin.reachability_map.sorted_vertices.size()) / y.size();

            basin.extremum_vertex = vertex;
            basin_cx.lmax_basins_map[vertex] = std::move(basin);
        }
    }

	// Populate basin_cx.init_rel_min_monotonicity_spans
	basin_cx.init_rel_min_monotonicity_spans.clear();

	// Add spans from minima
	for (const auto& [vertex, basin] : basin_cx.lmin_basins_map) {
		basin_cx.init_rel_min_monotonicity_spans.push_back(basin.rel_min_monotonicity_span);
	}

	// Add spans from maxima
	for (const auto& [vertex, basin] : basin_cx.lmax_basins_map) {
		basin_cx.init_rel_min_monotonicity_spans.push_back(basin.rel_min_monotonicity_span);
	}

	// Sort the spans for easier analysis
	std::sort(basin_cx.init_rel_min_monotonicity_spans.begin(),
			  basin_cx.init_rel_min_monotonicity_spans.end());

    return basin_cx;
}

