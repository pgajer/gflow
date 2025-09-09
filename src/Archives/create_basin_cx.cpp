// Tue A29 7:18PM
basin_cx_t set_wgraph_t::create_basin_cx(
    const std::vector<double>& y,
    double rel_min_monotonicity_span_thld
    ) const {

#define DEBUG__create_basin_cx 1

	// For DEBUGGING only !!!
	std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data";
	if (!std::filesystem::exists(debug_dir)) {
		if (!std::filesystem::create_directories(debug_dir)) {
			REPORT_ERROR("ERROR: Failed to create debug directory: %s\n", debug_dir.c_str());
		}
	}

	// 1. find nbr-local extrema
    auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

    basin_cx_t basin_cx;
    basin_cx.rel_min_monotonicity_span_thld = rel_min_monotonicity_span_thld;
    basin_cx.harmonic_predictions = y;  // Initialize with original values

    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    double y_range = y_max - y_min;

	// 2. for each nbr lmin of nbr_lmin find its basin
    bool detect_maxima = false;
    for (const auto& vertex : nbr_lmin) {
        auto basin = find_local_extremum_basin(
            vertex,
            y,
            detect_maxima
            );

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

			// Rprintf("\nnbr_lmin vertex %zu\n", vertex + 1);

			for (const auto& [u, span] : basin.boundary_monotonicity_spans_map) {

				// Rprintf("[%zu, %.3f]\n", u + 1, span);
				// bool is_in_nbr_lmax = std::find(nbr_lmax.begin(), nbr_lmax.end(), u) != nbr_lmax.end();
				// Rprintf("is_in_nbr_lmax: %s\n", is_in_nbr_lmax ? "TRUE" : "FALSE");

				if (std::find(nbr_lmax.begin(), nbr_lmax.end(), u) != nbr_lmax.end() && span < min_span) {
					min_span = span;
					min_span_vertex = u;
				}
			}

			basin.min_span_vertex           = min_span_vertex;
			basin.min_monotonicity_span     = min_span;
			basin.rel_min_monotonicity_span = (y_range > 0.0) ? min_span / y_range : 0.0;

			// Rprintf("basin.min_span_vertex: %zu  basin.min_monotonicity_span: %.3f\n",
			//  		basin.min_span_vertex + 1, basin.min_monotonicity_span);

            basin.extremum_vertex = vertex;
            basin_cx.lmin_basins_map[vertex] = std::move(basin);
        }
    }

	// debugging - testing if min_span_vertex values are still OK
	Rprintf("\nScannig basin_cx.lmin_basins_map\n");
	for (const auto& [vertex, basin] : basin_cx.lmin_basins_map) {
		Rprintf("vertex: %zu  basin.min_span_vertex: %zu  basin.min_monotonicity_span: %.3f\n",
				vertex + 1, basin.min_span_vertex + 1, basin.min_monotonicity_span);
    }
	Rprintf("\n");


	// 3. for each nbr lmax of nbr_lmax find its basin
    detect_maxima = true;
    for (const auto& vertex : nbr_lmax) {
        auto basin = find_local_extremum_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
			// See comments in the local min block for this part
			double min_span = INFINITY;
			size_t min_span_vertex = INVALID_VERTEX;

			// Rprintf("\nnbr_lmax vertex %zu\n", vertex + 1);

			for (const auto& [u, span] : basin.boundary_monotonicity_spans_map) {

				// Rprintf("[%zu, %.3f]\n", u + 1, span);
				// bool is_in_nbr_lmin = std::find(nbr_lmin.begin(), nbr_lmin.end(), u) != nbr_lmin.end();
				// Rprintf("is_in_nbr_lmin: %s\n", is_in_nbr_lmin ? "TRUE" : "FALSE");

				if (std::find(nbr_lmin.begin(), nbr_lmin.end(), u) != nbr_lmin.end() && span < min_span) {
					min_span = span;
					min_span_vertex = u;
				}
			}

			basin.min_span_vertex           = min_span_vertex;
			basin.min_monotonicity_span     = min_span;
			basin.rel_min_monotonicity_span = (y_range > 0.0) ? min_span / y_range : 0.0;

			// Rprintf("basin.min_span_vertex: %zu  basin.min_monotonicity_span: %.3f\n",
			// 		basin.min_span_vertex + 1, basin.min_monotonicity_span);

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


	// debugging
	basin_cx.write_basins_map(debug_dir);
	//	error("DEBUGGING\n");

#if 0
	// 4. Create a priority queue of basins sorted by increasing rel_min_monotonicity_span
    struct basin_entry_t {
        size_t extremum_vertex;
        double rel_min_monotonicity_span;
		size_t min_span_vertex;
		bool is_maximum;  ///< True if the extremum is a maximum; false if a minimum

        bool operator>(const basin_entry_t& other) const {
            return rel_min_monotonicity_span > other.rel_min_monotonicity_span;
        }

		explicit basin_entry_t(size_t e, double r, size_t m, bool ismax)
			: extremum_vertex(e),
			  rel_min_monotonicity_span(r),
			  min_span_vertex(m),
			  is_maximum(ismax)
			{}
    };

    std::priority_queue<basin_entry_t, std::vector<basin_entry_t>, std::greater<basin_entry_t>> basin_queue;

	// Add all local minima to the queue
    for (const auto& [vertex, basin] : basin_cx.lmin_basins_map) {
        basin_queue.emplace(vertex, basin.rel_min_monotonicity_span, basin.min_span_vertex, false);
    }

	// Add all local maxima to the queue
    for (const auto& [vertex, basin] : basin_cx.lmax_basins_map) {
        basin_queue.emplace(vertex, basin.rel_min_monotonicity_span, basin.min_span_vertex, true);
    }

	// Track processed extrema to avoid re-processing
    std::unordered_set<size_t> processed_extrema;

	// 5. Process basins in order of increasing rel_min_monotonicity_span
    while (!basin_queue.empty()) {
        auto basin_entry = basin_queue.top();
        basin_queue.pop();

        size_t A_basin_extremum = basin_entry.extremum_vertex;

		// Skip if this extremum has already been processed
        if (processed_extrema.count(A_basin_extremum) > 0) {
            continue;
        }

		// Break if rel_min_monotonicity_span exceeds threshold
        if (basin_entry.rel_min_monotonicity_span > rel_min_monotonicity_span_thld) {
            break;
        }

        basin_t A_basin;
		// Determine if A_basin is a minimum or maximum basin
        if (basin_entry.is_maximum) {
            A_basin = basin_cx.lmax_basins_map.at(A_basin_extremum);
        } else {
            A_basin = basin_cx.lmin_basins_map.at(A_basin_extremum);
        }

		// Skip if the entry is stale meaning
		// 1) it is not even found on basin_cx.lmin/lmax_basins_map due to being absorbed and so removed from basin_cx.lmin/lmax_basins_map
		// 2) basin_entrpy min_span_vertex is out of sync with min_span_vertex on the current basin_cx.lmin/lmax_basins_map entry of that extremum
		if (
			( // cannot be found in any basins_map
				basin_cx.lmin_basins_map.find(A_basin.extremum_vertex) == basin_cx.lmin_basins_map.end() &&
				basin_cx.lmax_basins_map.find(A_basin.extremum_vertex) == basin_cx.lmax_basins_map.end()
				) || (
					(basin_entry.is_maximum && basin_entry.min_span_vertex != basin_cx.lmax_basins_map.at(basin_entry.extremum_vertex).min_span_vertex) ||
					(!basin_entry.is_maximum && basin_entry.min_span_vertex != basin_cx.lmin_basins_map.at(basin_entry.extremum_vertex).min_span_vertex)
					)
			) {
			continue;
		}

		Rprintf("\n----------------\n"
				"Processing A_basin_extremum: %zu which is %s with A_basin.min_span_vertex: %zu\n",
				A_basin_extremum + 1, basin_entry.is_maximum ? "lmax" : "lmin", A_basin.min_span_vertex + 1);

		// Check if the next in queue basin is forms a cancellation pair with A_basin
		// Search for cancellation basin of the opposite sign
        auto B_basin_entry        = basin_queue.top();
		size_t B_basin_extremum = B_basin_entry.extremum_vertex;
		basin_t B_basin;
		// Determine if B_basin is a minimum or maximum basin
        if (B_basin_entry.is_maximum) {
            B_basin = basin_cx.lmax_basins_map.at(B_basin_extremum);
        } else {
            B_basin = basin_cx.lmin_basins_map.at(B_basin_extremum);
        }

		basin_t* absorbing_basin_ptr = nullptr;
		if (A_basin.min_span_vertex == B_basin.extremum_vertex && B_basin.min_span_vertex == A_basin.extremum_vertex) {
			basin_queue.pop();
			absorbing_basin_ptr = basin_cx.process_cancellation_pair(A_basin, B_basin, y, *this);
		} else {
			absorbing_basin_ptr = basin_cx.process_single_basin(A_basin, y, *this);
		}

		basin_queue.emplace(absorbing_basin_ptr->extremum_vertex,
							absorbing_basin_ptr->rel_min_monotonicity_span,
							absorbing_basin_ptr->min_span_vertex,
							absorbing_basin_ptr->is_maximum);
	} // END OF while (!basin_queue.empty())
#endif

    return basin_cx;
}


// Mon A28 12:32PM
basin_cx_t set_wgraph_t::create_basin_cx(
    const std::vector<double>& y,
    double range_rel_depth_thld
    ) const {

    // 1. find nbr-local extrema
    auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

    basin_cx_t basin_cx;
    basin_cx.range_rel_depth_thld = range_rel_depth_thld;

    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    double y_range = y_max - y_min;

    // 2. for each nbr lmin of nbr_lmin find its basin
    bool detect_maxima = false;
    for (const auto& vertex : nbr_lmin) {
        auto basin = find_local_extremum_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
            basin.rel_min_monotonicity_span = (y_range > 0.0) ? basin.min_monotonicity_span / y_range : 0.0;
            basin.extremum_vertex = vertex;
            basin_cx.lmin_basins_map[vertex] = basin;
        }
    }

    // 3. for each nbr lmax of nbr_lmax find its basin
    detect_maxima = true;
    for (const auto& vertex : nbr_lmax) {
        auto basin = find_local_extremum_basin(
            vertex,
            y,
            detect_maxima
            );

        if (basin.reachability_map.sorted_vertices.size()) {
            basin.rel_min_monotonicity_span = (y_range > 0.0) ? basin.min_monotonicity_span / y_range : 0.0;
            basin.extremum_vertex = vertex;
            basin_cx.lmax_basins_map[vertex] = basin;
        }
    }

    // 4. Populated the vector of nbr local extrema sorted in the ascending order w/r rel_min_monotonicity_span
    basin_cx.sorted_basins.clear();

    // Add all local minima and maxima vertices to the sorted_basins vector
    for (const auto& [vertex, basin] : basin_cx.lmin_basins_map) {
        basin_cx.sorted_basins.push_back(vertex);
    }
    for (const auto& [vertex, basin] : basin_cx.lmax_basins_map) {
        basin_cx.sorted_basins.push_back(vertex);
    }

    // Sort the basins by increasing rel_min_monotonicity_span
    std::sort(basin_cx.sorted_basins.begin(), basin_cx.sorted_basins.end(),
              [&](size_t v1, size_t v2) {
                  double span1 = basin_cx.lmin_basins_map.count(v1) ?
                      basin_cx.lmin_basins_map.at(v1).rel_min_monotonicity_span :
                      basin_cx.lmax_basins_map.at(v1).rel_min_monotonicity_span;

                  double span2 = basin_cx.lmin_basins_map.count(v2) ?
                      basin_cx.lmin_basins_map.at(v2).rel_min_monotonicity_span :
                      basin_cx.lmax_basins_map.at(v2).rel_min_monotonicity_span;

                  return span1 < span2;
              });

    // Initialize harmonic predictions with original values
    basin_cx.harmonic_predictions = y;

    // Track processed extrema to avoid re-processing
    std::unordered_set<size_t> processed_extrema;

    // 5. Identify cancelling pairs and perform repairs
    for (const auto& A_basin_extremum : basin_cx.sorted_basins) {
        // Skip if this extremum has already been processed
        if (processed_extrema.count(A_basin_extremum) > 0) {
            continue;
        }

        basin_t A_basin;
        // Determine if A_basin is a minimum or maximum basin
        if (basin_cx.lmin_basins_map.count(A_basin_extremum) > 0) {
            A_basin = basin_cx.lmin_basins_map.at(A_basin_extremum);
        } else {
            A_basin = basin_cx.lmax_basins_map.at(A_basin_extremum);
        }

        if (A_basin_extremum != A_basin.extremum_vertex) { // sanity check
            REPORT_ERROR("A_basin_extremum: %zu\t"
                         "A_basin.extremum_vertex: %zu should be the same!!!\n",
                         A_basin_extremum,
                         A_basin.extremum_vertex);
        }

        if (A_basin.rel_min_monotonicity_span > range_rel_depth_thld) {
            // We only process basins with rel_min_monotonicity_span <= range_rel_depth_thld
            break;
        }

        // Skip if min_span_vertex is invalid
        if (A_basin.min_span_vertex == INVALID_VERTEX) {
            continue;
        }

        // Search for cancellation basin of the opposite sign
        basin_t B_basin;
        bool found_B_basin = false;

        if (A_basin.is_maximum) {
            // If A is a maximum, look for minimum basin at min_span_vertex
            auto it = basin_cx.lmin_basins_map.find(A_basin.min_span_vertex);
            if (it != basin_cx.lmin_basins_map.end()) {
                B_basin = it->second;
                found_B_basin = true;
            }
        } else {
            // If A is a minimum, look for maximum basin at min_span_vertex
            auto it = basin_cx.lmax_basins_map.find(A_basin.min_span_vertex);
            if (it != basin_cx.lmax_basins_map.end()) {
                B_basin = it->second;
                found_B_basin = true;
            }
        }

        if (!found_B_basin) {
            continue; // No cancellation partner found
        }

        // Check if B's min_span_vertex points back to A
        if (B_basin.min_span_vertex != A_basin.extremum_vertex) {
            continue; // Not a mutual cancellation pair
        }

        // Found a valid cancellation pair - record it
        basin_cx.cancellation_pairs.push_back({
            A_basin.is_maximum ? B_basin.extremum_vertex : A_basin.extremum_vertex,
            A_basin.is_maximum ? A_basin.extremum_vertex : B_basin.extremum_vertex
        });

        // 1) Find cousin basin of opposite type for A_basin
        basin_t A_cousin_basin;
        double A_cousin_distance = INFINITY;
        bool found_A_cousin_basin = false;
        for (const auto& [d, u] : A_basin.boundary_vertices_map) {
            if (u == B_basin.extremum_vertex) {
                continue; // Skip the cancellation partner
            }

            // Check if u is an extremum of opposite type to A
            if (A_basin.is_maximum) {
                auto it = basin_cx.lmin_basins_map.find(u);
                if (it != basin_cx.lmin_basins_map.end()) {
                    A_cousin_basin = it->second;
                    A_cousin_distance = d;
                    found_A_cousin_basin = true;
                    break;
                }
            } else {
                auto it = basin_cx.lmax_basins_map.find(u);
                if (it != basin_cx.lmax_basins_map.end()) {
                    A_cousin_basin = it->second;
                    A_cousin_distance = d;
                    found_A_cousin_basin = true;
                    break;
                }
            }
        }

        // 2) Find cousin basin of opposite type for B_basin
        basin_t B_cousin_basin;
        double B_cousin_distance = INFINITY;
        bool found_B_cousin_basin = false;
        for (const auto& [d, u] : B_basin.boundary_vertices_map) {
            if (u == A_basin.extremum_vertex) {
                continue; // Skip the cancellation partner
            }

            // Check if u is an extremum of opposite type to B
            if (B_basin.is_maximum) {
                auto it = basin_cx.lmin_basins_map.find(u);
                if (it != basin_cx.lmin_basins_map.end()) {
                    B_cousin_basin = it->second;
                    B_cousin_distance = d;
                    found_B_cousin_basin = true;
                    break;
                }
            } else {
                auto it = basin_cx.lmax_basins_map.find(u);
                if (it != basin_cx.lmax_basins_map.end()) {
                    B_cousin_basin = it->second;
                    B_cousin_distance = d;
                    found_B_cousin_basin = true;
                    break;
                }
            }
        }

        // 3) Determine which cousin basin to use
        size_t absorbing_extremum = INVALID_VERTEX;
        size_t absorbed_extremum = INVALID_VERTEX;

        if (found_A_cousin_basin && found_B_cousin_basin) {
            if (A_cousin_distance < B_cousin_distance) {
                absorbing_extremum = A_cousin_basin.extremum_vertex;
                absorbed_extremum = A_basin.extremum_vertex;
            } else {
                absorbing_extremum = B_cousin_basin.extremum_vertex;
                absorbed_extremum = B_basin.extremum_vertex;
            }
        } else if (found_A_cousin_basin) {
            absorbing_extremum = A_cousin_basin.extremum_vertex;
            absorbed_extremum = A_basin.extremum_vertex;
        } else if (found_B_cousin_basin) {
            absorbing_extremum = B_cousin_basin.extremum_vertex;
            absorbed_extremum = B_basin.extremum_vertex;
        } else {
            // No cousin basin found - we would need a default strategy here
            // For now, just skip this pair
            continue;
        }

        // Record the extrema for repair
        basin_cx.repair_lextr[{A_basin.extremum_vertex, B_basin.extremum_vertex}] = absorbing_extremum;

        // Mark both extrema as processed
        processed_extrema.insert(A_basin.extremum_vertex);
        processed_extrema.insert(B_basin.extremum_vertex);
    }

    return basin_cx;
}


// Mon A28 12:28PM
basin_t set_wgraph_t::find_local_extremum_basin(
    size_t vertex,
    const std::vector<double>& y,
    bool detect_maxima
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
            bool condition_met = detect_maxima ? (delta_y < 0.0) : (delta_y > 0.0);

            if (!condition_met) {
                // *** MONOTONICITY VIOLATION BOUNDARY ***
                // When monotonicity is violated, calculate depth/prominence

                // For maxima: depth = y[vertex] - y[prev[u]]
                // For minima: depth = y[prev[u]] - y[vertex]
                double curr_depth = detect_maxima ? (y[vertex] - y[prev[u]]) : (y[prev[u]] - y[vertex]);

                // Update minimum depth if this violation gives a smaller valid depth
                constexpr double eps = 1e-12;  // Small epsilon to handle floating-point precision
                if (curr_depth > eps && curr_depth < basin.min_monotonicity_span) {
                    basin.min_monotonicity_span = curr_depth;
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
        // Check if this vertex is a boundary due to having neighbors outside the basin
        bool is_boundary = false;
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

        // If this is a boundary vertex, add it to the boundary map
        if (is_boundary) {
            basin.boundary_vertices_map[d] = u;
			boundary_vertices_span_map[u] = std::abs(y[vertex] - y[u]);
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

#if 0
	// If we do this, we will not be able to determine if the basin should have a cancellation partner
	// If no violations were found during exploration, calculate basin.min_monotonicity_span from boundary vertices
	if (basin.min_monotonicity_span == INFINITY) {
		// Calculate min_monotonicity_span using all boundary vertices
		for (const auto& [distance, boundary_vertex] : basin.boundary_vertices_map) {
			double span = detect_maxima ?
				(y[vertex] - y[boundary_vertex]) :
				(y[boundary_vertex] - y[vertex]);

			if (span < basin.min_monotonicity_span) {
				basin.min_monotonicity_span = span;
				basin.min_span_vertex = boundary_vertex;
			}
		}
	}
#endif

    return basin;
}


// Mon A28 12:06PM
basin_cx_t set_wgraph_t::create_basin_cx(
    const std::vector<double>& y,
	double range_rel_depth_thld
    ) const {

	// 1. find nbr-local extrema
	//std::pair<std::vector<size_t>, std::vector<size_t>>
	auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

	basin_cx_t basin_cx;
	basin_cx.range_rel_depth_thld = range_rel_depth_thld;

	double y_min = *std::min_element(y.begin(), y.end());
	double y_max = *std::max_element(y.begin(), y.end());
	double y_range = y_max - y_min;

	// 2. for each nbr lmin of nbr_lmin find its basin
	bool detect_maxima = false;
	for (const auto& vertex : nbr_lmin) {
		auto basin = find_local_extremum_basin(
			vertex,
			y,
			detect_maxima
			);

		if (basin.reachability_map.sorted_vertices.size()) {
			basin.rel_min_monotonicity_span  = (y_range > 0.0) ? basin.min_monotonicity_span / y_range : 0.0;
			basin.extremum_vertex = vertex;
			basin_cx.lmin_basins_map[vertex] = basin;
		}
	}

	// 3. for each nbr lmax of nbr_lmax find its basin
	detect_maxima = true;
	for (const auto& vertex : nbr_lmax) {
		auto basin = find_local_extremum_basin(
			vertex,
			y,
			detect_maxima
			);

		if (basin.reachability_map.sorted_vertices.size()) {
			basin.rel_min_monotonicity_span  = (y_range > 0.0) ? basin.min_monotonicity_span / y_range : 0.0;
			basin.extremum_vertex = vertex;
			basin_cx.lmax_basins_map[vertex] = basin;
		}
	}

	// 4. Populated the vector of nbr local extrema sorted in the ascending order w/r rel_min_monotonicity_span
	basin_cx.sorted_basins.clear();

	// Add all local minima and maxima vertices to the sorted_basins vector
	for (const auto& [vertex, basin] : basin_cx.lmin_basins_map) {
		basin_cx.sorted_basins.push_back(vertex);
	}
	for (const auto& [vertex, basin] : basin_cx.lmax_basins_map) {
		basin_cx.sorted_basins.push_back(vertex);
	}

	// Sort the basins by increasing rel_min_monotonicity_span
	std::sort(basin_cx.sorted_basins.begin(), basin_cx.sorted_basins.end(),
			  [&](size_t v1, size_t v2) {
				  double span1 = basin_cx.lmin_basins_map.count(v1) ?
					  basin_cx.lmin_basins_map.at(v1).rel_min_monotonicity_span :
					  basin_cx.lmax_basins_map.at(v1).rel_min_monotonicity_span;

				  double span2 = basin_cx.lmin_basins_map.count(v2) ?
					  basin_cx.lmin_basins_map.at(v2).rel_min_monotonicity_span :
					  basin_cx.lmax_basins_map.at(v2).rel_min_monotonicity_span;

				  return span1 < span2;
			  });


	// 5. Identify cancelling pairs
	basin_t A_basin;
	basin_t B_basin;
	for (const auto& A_basin_extremum : basin_cx.sorted_basins) {

		A_basin = basin_cx.lmin_basins_map.count(A_basin_extremum) ?
			basin_cx.lmin_basins_map.at(A_basin_extremum) :
			basin_cx.lmax_basins_map.at(A_basin_extremum);

		if (A_basin_extremum != A_basin.extremum_vertex) { // sanity check
			REPORT_ERROR("A_basin_extremum: %zu\t"
						 "A_basin.extremum_vertex: %zu should be the same!!!\n",
						 A_basin_extremum,
						 A_basin.extremum_vertex);
		}

		if (A_basin.rel_min_monotonicity_span > range_rel_depth_thld) { // we only process basins with  rel_min_monotonicity_span <= range_rel_depth_thld
			break;
		}

		// if (A_basin.min_span_vertex != INVALID_VERTEX) { // since A_basin.rel_min_monotonicity_span is finite this should always be true
		// search for cancellation basin of the opposite sign
		// We are searching for a cancellation_cousin such that
		// A_basin.min_span_vertex = B_basin.extremum_vertex AND B_basin.min_span_vertex = A_basin.extremum_vertex
		if (A_basin.is_maximum) {
			auto it = basin_cx.lmin_basins_map.find(A_basin.min_span_vertex);
			if (it != basin_cx.lmin_basins_map.end()) {
				B_basin = it->second;
			} else {
				REPORT_ERROR("lmax %zu basin with rel_min_monotonicity_span %.3f "
							 "does not have B_basin\n",
							 A_basin_extremum, A_basin.rel_min_monotonicity_span);
			}

			if (B_basin.min_span_vertex != A_basin.extremum_vertex) {
				Rprintf("A_basin.rel_min_monotonicity_span: %.3f\n"
						"range_rel_depth_thld: %.3f\n"
						"B_basin.rel_min_monotonicity_span: %.3f\n",
						A_basin.rel_min_monotonicity_span,
						range_rel_depth_thld,
						B_basin.rel_min_monotonicity_span);
				REPORT_ERROR("lmin B_basin %zu basin's min_span_vertex: %zu is not equalt to A_basin.extremum_vertex: %zu\n",
							 B_basin.extremum_vertex,
							 B_basin.min_span_vertex,
							 A_basin.extremum_vertex);
			}

			// if we got here, it means that
			// A_basin.min_span_vertex = B_basin.extremum_vertex AND B_basin.min_span_vertex = A_basin.extremum_vertex
			// 1) Find lmin vertex (that is NOT B_basin) at the boundary of 'A_basin' that is closest to A_basin_extremum
			basin_t A_cousin_basin;
			double A_cousin_basin_distance = 0;
			bool found_A_cousin_basin = false;
			for (const auto& [d,u] : A_basin.boundary_vertices_map ) {
				if (u == B_basin.extremum_vertex) {
					continue;
				}

				// let's check if u is a lmin vertex
				auto it = basin_cx.lmin_basins_map.find(u);
				if (it != basin_cx.lmin_basins_map.end()) {
					A_cousin_basin = it->second;
					A_cousin_basin_distance = d;
					break;
				}
			}

			// 2) Find lmax vertex (that is NOT basin_etremum) at the boundary of B_basin that is closest to B_basin.extremum_vertex
			basin_t B_cousin_basin;
			double B_cousin_basin_distance = 0;
			bool found_B_cousin_basin = false;
			for (const auto& [d,u] : B_basin.boundary_vertices_map ) {
				if (u == A_basin.extremum_vertex) {
					continue;
				}

				// let's check if u is a lmax vertex
				auto it = basin_cx.lmax_basins_map.find(u);
				if (it != basin_cx.lmax_basins_map.end()) {
					B_cousin_basin = it->second;
					B_cousin_basin_distance = d;
					break;
				}
			}

			// 3) pick the one with the shortest distance
			basin_t absorbing_basin;
			basin_t absorbed_basin;
			if (found_A_cousin_basin && found_B_cousin_basin) {
				if (A_cousin_basin_distance < B_cousin_basin_distance) {
					absorbing_basin = std::move(A_cousin_basin);
					absorbed_basin  = std::move(A_basin);
				} else {
					absorbing_basin = std::move(B_cousin_basin);
					absorbed_basin  = std::move(B_basin);
				}
			} else if (!found_A_cousin_basin && found_B_cousin_basin) {
				absorbing_basin = std::move(B_cousin_basin);
			} else if (found_A_cousin_basin && !found_B_cousin_basin) {
				absorbing_basin = std::move(A_cousin_basin);
			} else if(!found_A_cousin_basin && !found_B_cousin_basin) {
				REPORT_ERROR("We didn't find found_A_cousin_basin for A_basin_extremum: %zu AND "
							 "We didn't find found_B_cousin_basin for B_basin.extremum_vertex: %zu\n",
							 A_basin_extremum,
							 B_basin.extremum_vertex);
			}

			// 4) Add the basin of the winner of the cancellation pair (the one with the shortest distance) to the closest found basin


			// 5) Perform harmonic repair/smoothing of the values of y within the winner basin
		} else {

		}
	}

	return basin_cx;
}


// Mon A28 11:31

// Helper function for harmonic repair of function values in absorbed basin
void set_wgraph_t::perform_harmonic_repair(
    std::vector<double>& harmonic_predictions,
    const basin_t& absorbing_basin,
    const basin_t& absorbed_basin
	) const {
	// Collect all vertices in the absorbed basin
    std::unordered_set<size_t> basin_vertices;
    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        basin_vertices.insert(v_info.vertex);
    }

	// Identify boundary vertices (vertices in the basin with neighbors outside)
    std::unordered_map<size_t, double> boundary_values;

	// Add the absorbing basin's extremum as a boundary vertex with fixed value
    boundary_values[absorbing_basin.extremum_vertex] = harmonic_predictions[absorbing_basin.extremum_vertex];

	// Add other boundary vertices from the absorbed basin
    for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
        size_t v = v_info.vertex;

		// Skip the extremum itself - it will be recalculated
        if (v == absorbed_basin.extremum_vertex) {
            continue;
        }

		// Check if this vertex has neighbors outside the basin
        for (const auto& edge : adjacency_list[v]) {
            if (basin_vertices.count(edge.vertex) == 0) {
				// This is a boundary vertex - keep its original value
                boundary_values[v] = harmonic_predictions[v];
                break;
            }
        }
    }

	// Include common vertices between absorbing and absorbed basins as boundary values
    for (const auto& v_info : absorbing_basin.reachability_map.sorted_vertices) {
        size_t v = v_info.vertex;
        if (basin_vertices.count(v) > 0 && v != absorbed_basin.extremum_vertex) {
            boundary_values[v] = harmonic_predictions[v];
        }
    }

	// Simple Jacobi iteration for harmonic interpolation
    const int MAX_ITERATIONS = 100;
    const double TOLERANCE = 1e-6;

    std::vector<double> new_values = harmonic_predictions;

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        bool converged = true;

        for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
            size_t v = v_info.vertex;

			// Skip boundary vertices - their values are fixed
            if (boundary_values.count(v) > 0) {
                continue;
            }

			// Average the values of neighbors
            double sum = 0.0;
            double weight_sum = 0.0;

            for (const auto& edge : adjacency_list[v]) { // line (gflow_basins.cpp:1324)
				// Use inverse of edge weight as weighting factor
                double weight = 1.0 / (edge.weight + 1e-10); // Avoid division by zero
                sum += harmonic_predictions[edge.vertex] * weight;
                weight_sum += weight;
            }

            if (weight_sum > 0) {
                double weighted_avg = sum / weight_sum;
                if (std::abs(new_values[v] - weighted_avg) > TOLERANCE) {
                    converged = false;
                }
                new_values[v] = weighted_avg;
            }
        }

		// Update the values for next iteration
        for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
            size_t v = v_info.vertex;
            if (boundary_values.count(v) == 0) {
                harmonic_predictions[v] = new_values[v];
            }
        }

        if (converged) {
            break;
        }
    }

	// Make a special adjustment for the absorbed extremum - interpolate between
	// absorbing extremum value and the average of the absorbed extremum's neighbors
    size_t absorbed_extremum = absorbed_basin.extremum_vertex;
    double neighbor_avg = 0.0;
    double weight_sum = 0.0;

    for (const auto& edge : adjacency_list[absorbed_extremum]) {
        double weight = 1.0 / (edge.weight + 1e-10);
        neighbor_avg += harmonic_predictions[edge.vertex] * weight;
        weight_sum += weight;
    }

    if (weight_sum > 0) {
        neighbor_avg /= weight_sum;

		// Calculate distance between extrema
        double extrema_distance = 1.0;
        if (absorbing_basin.reachability_map.distances.count(absorbed_extremum) > 0) {
            extrema_distance = absorbing_basin.reachability_map.distances.at(absorbed_extremum);
        } else {
			// Compute approximate distance if not directly available
            double d = compute_shortest_path_distance(
                absorbing_basin.extremum_vertex,
                absorbed_extremum
				);
            extrema_distance = std::max(1.0, d);
        }

		// Blend based on distance
        double blend_factor = 1.0 / (1.0 + extrema_distance);
        harmonic_predictions[absorbed_extremum] =
            blend_factor * harmonic_predictions[absorbing_basin.extremum_vertex] +
            (1.0 - blend_factor) * neighbor_avg;
    }
}

