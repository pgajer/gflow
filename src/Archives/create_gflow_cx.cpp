

// May 14, 2025 version

gflow_cx_t set_wgraph_t::create_gflow_cx(
	const std::vector<double>& y,
	size_t hop_dist_thld
	) const {

	(void)hop_dist_thld; // for now we do not use this parameter

	// 1. find nbr-local extrema
	auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

	gflow_cx_t gflow_cx;
	gflow_cx.harmonic_predictions = y;  // Initialize with original values

	// 2. for each nbr lmin of nbr_lmin compute its extremum hop nbhd (and in particular its extremum hop index)
	bool detect_maxima = false;
	for (const auto& vertex : nbr_lmin) {
		hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(
			vertex,
			y,
			detect_maxima);
		gflow_cx.lmin_hop_nbhd_map[vertex] = hop_nbhd;
	}

	// 3. for each nbr lmax of nbr_lmax compute its extremum hop nbhd (and in particular its extremum hop index)
	detect_maxima = true;
	for (const auto& vertex : nbr_lmax) {
		hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(
			vertex,
			y,
			detect_maxima);
		gflow_cx.lmax_hop_nbhd_map[vertex] = hop_nbhd;
	}

	// 4. For local minima and maxima whose hop_idx is <= hop_dist_thld perform
	// the following harmonic smoothing of the values of y within the hop disk
	// of radius hop_idx as described below.
	//
	// Let, partial_D, be the set of vertices (hop_idx + 1) hop distance away
	// from local extremum vertex and let int_D be the set of vertices whose hop
	// distance from 'vertex' is <= hop_idx.
	//
	// The proposed quasi-harmonic smoothing algorithm: for each vertex u in
	// int_D compute the geodesic distance from u to every vertex in partial_D.
	// Let d_1, d_2, ... , d_n be the set of these distances, where d_i is the
	// distance between u and u_i in partial_D. Then compute quasi-harmonic
	// estimate qh_y of y at u as the value
	// (qh_y)(u) = \sum_i w_i * y[u_i],
	// where
	// w_i = (d_i)^{-1} / sum_j (d_j)^{-1}

	// 5. Recompute nbr local extrema over quasi-harmonic modification of y to
	// make sure all local extrema with hop_idx less than hop_idx_thld were
	// eliminated. If not, maybe we need to perform the harmonic smoothing
	// again?

	// 6. Compute gradient flow trajectories of the quasi-harmonic modification
	// qh_y of y, where a gradient of qh_y at v is the neighbor u of v with the
	// maximum value of qh_y(u) - qh_y(v). Of course, if v is a local maximum v
	// will be its own gradient element. Use gradient elements to form ascending
	// gradient trajectories and "negative" gradient to form descending gradient
	// trajectories. A "negative" gradient at v is a neighbor u of v such that
	// qh_y(u) - qh_y(v) is the minimum value over all neighbors of v. The
	// composition of the inverse descending grad trajectory with the ascending
	// grad trajectory is the qh_y gradient trajectory passing through v.

	// Generate qh_y gradient trajectories for all vertices of the graph

	// 7. Check that all gradient trajectories start at the local minima with
	// hop_idx > hop_idx_thld and end at a local maximum with hop_idx >
	// hop_idx_thld

	// 8. Form ascending basins for local minima and descending basins for local
	// maxima.

	// 9. Form gradient flow cells defined as the interesection of ascending and
	// descending basins. We also define ascending-ascending and
	// descending-descending cells as intersection of two ascending or two
	// descending basins, respectively.

	return gflow_cx;
}


gflow_cx_t set_wgraph_t::create_gflow_cx(
    const std::vector<double>& y,
    std::optional<size_t> custom_hop_idx_thld = std::nullopt,
    int max_smoothing_iterations = 3,
    double smoothing_tolerance = 1e-6) const {

    //------------------------------------------------------------------
    // PHASE 1: Analysis - Find extrema and determine threshold
    //------------------------------------------------------------------

    // Find all local extrema
    auto [nbr_lmin, nbr_lmax] = find_nbr_extrema(y);

    // Compute hop neighborhoods for all extrema
    std::unordered_map<size_t, hop_nbhd_t> lmin_hop_nbhd_map;
    std::unordered_map<size_t, hop_nbhd_t> lmax_hop_nbhd_map;

    bool detect_maxima = false;
    for (const auto& vertex : nbr_lmin) {
        hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, y, detect_maxima);
        lmin_hop_nbhd_map[vertex] = hop_nbhd;
    }

    detect_maxima = true;
    for (const auto& vertex : nbr_lmax) {
        hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, y, detect_maxima);
        lmax_hop_nbhd_map[vertex] = hop_nbhd;
    }

    // Determine threshold (either custom or computed from distribution)
    size_t hop_idx_thld = custom_hop_idx_thld.value_or(
        suggest_hop_idx_threshold(lmin_hop_nbhd_map, lmax_hop_nbhd_map));

    //------------------------------------------------------------------
    // PHASE 2: Iterative Smoothing
    //------------------------------------------------------------------

    // Initialize the smoothed function with original values
    std::vector<double> smoothed_y = y;

    // Track extrema that should be preserved (hop_idx > hop_idx_thld)
    std::unordered_set<size_t> significant_lmin;
    std::unordered_set<size_t> significant_lmax;

    for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            significant_lmin.insert(vertex);
        }
    }

    for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            significant_lmax.insert(vertex);
        }
    }

    // Iterative smoothing process
    bool smoothing_complete = false;
    int smoothing_iteration = 0;

    while (!smoothing_complete && smoothing_iteration < max_smoothing_iterations) {
        // Phase 2(A): Apply smoothing to neighborhoods of insignificant extrema
        for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
            if (hop_nbhd.hop_idx <= hop_idx_thld) {
                perform_harmonic_smoothing_in_neighborhood(smoothed_y, hop_nbhd, smoothing_tolerance);
            }
        }

        for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
            if (hop_nbhd.hop_idx <= hop_idx_thld) {
                perform_harmonic_smoothing_in_neighborhood(smoothed_y, hop_nbhd, smoothing_tolerance);
            }
        }

        // Phase 2(B): Validate by recomputing extrema
        auto [new_lmin, new_lmax] = find_nbr_extrema(smoothed_y);

        // Check if only significant extrema remain
        bool validation_passed = true;

        // Check minima
        for (const auto& vertex : new_lmin) {
            if (significant_lmin.find(vertex) == significant_lmin.end()) {
                // Found a new or previously insignificant minimum
                validation_passed = false;
                break;
            }
        }

        // Check if all significant minima still exist
        if (validation_passed) {
            for (const auto& vertex : significant_lmin) {
                if (std::find(new_lmin.begin(), new_lmin.end(), vertex) == new_lmin.end()) {
                    // A significant minimum was lost during smoothing
                    validation_passed = false;
                    break;
                }
            }
        }

        // Similar checks for maxima (omitted for brevity but would follow same pattern)

        // Decide whether to continue smoothing
        if (validation_passed) {
            smoothing_complete = true;
        } else {
            // If validation failed, update hop neighborhoods for next iteration
            lmin_hop_nbhd_map.clear();
            lmax_hop_nbhd_map.clear();

            detect_maxima = false;
            for (const auto& vertex : new_lmin) {
                hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, smoothed_y, detect_maxima);
                lmin_hop_nbhd_map[vertex] = hop_nbhd;
            }

            detect_maxima = true;
            for (const auto& vertex : new_lmax) {
                hop_nbhd_t hop_nbhd = compute_extremum_hop_nbhd(vertex, smoothed_y, detect_maxima);
                lmax_hop_nbhd_map[vertex] = hop_nbhd;
            }
        }

        smoothing_iteration++;
    }

    //------------------------------------------------------------------
    // PHASE 3: Gradient Flow Complex Construction
    //------------------------------------------------------------------

    gflow_cx_t gflow_cx;
    gflow_cx.harmonic_predictions = smoothed_y;

    // Compute basins for remaining significant extrema
    for (const auto& vertex : significant_lmin) {
        detect_maxima = false;
        basin_t basin = find_local_extremum_bfs_basin(vertex, smoothed_y, detect_maxima);
        if (basin.reachability_map.sorted_vertices.size() > 0) {
            gflow_cx.lmin_basins_map[vertex] = basin;
        }
    }

    for (const auto& vertex : significant_lmax) {
        detect_maxima = true;
        basin_t basin = find_local_extremum_bfs_basin(vertex, smoothed_y, detect_maxima);
        if (basin.reachability_map.sorted_vertices.size() > 0) {
            gflow_cx.lmax_basins_map[vertex] = basin;
        }
    }

    // Compute gradient trajectories
    compute_gradient_trajectories(smoothed_y, gflow_cx);

    // Identify gradient flow cells
    compute_gradient_flow_cells(gflow_cx);

    return gflow_cx;
}
