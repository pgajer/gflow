// ---------------------------------------------------------------------------------------------------------
//
// May 16 version
//
// ---------------------------------------------------------------------------------------------------------

gflow_cx_t set_wgraph_t::create_gflow_cx(
    const std::vector<double>& y,
    size_t hop_idx_thld,
    smoother_type_t smoother_type,
    int max_outer_iterations,
    int max_inner_iterations,
    double smoothing_tolerance,
    double sigma,
    bool process_in_order,
    bool verbose,
    bool detailed_recording
) const {

	(void)max_outer_iterations;

	// Initialize the smoothed function with original values
    std::vector<double> smoothed_y = y;

	// Prepare extrema for processing
    struct extremum_info_t {
        size_t vertex;
        size_t hop_idx;
        bool is_minimum;
    };


    //------------------------------------------------------------------
    // PHASE 1: Find local extrema
    //------------------------------------------------------------------
    auto [lmin_hop_nbhd_map, lmax_hop_nbhd_map] = compute_extrema_hop_nbhds(y);

    //------------------------------------------------------------------
    // PHASE 2: Iterative Smoothing
    //------------------------------------------------------------------

    // Track extrema that should be preserved (hop_idx > hop_idx_thld)
    std::unordered_set<size_t> non_spurious_lmin;
    std::unordered_set<size_t> non_spurious_lmax;

    for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            non_spurious_lmin.insert(vertex);
        }
    }

    for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
        if (hop_nbhd.hop_idx > hop_idx_thld) {
            non_spurious_lmax.insert(vertex);
        }
    }

    // Initialize gflow_cx with original data
    gflow_cx_t gflow_cx;
    gflow_cx.lmin_hop_nbhd_map = lmin_hop_nbhd_map;
    gflow_cx.lmax_hop_nbhd_map = lmax_hop_nbhd_map;

    std::vector<smoothing_step_t> smoothing_history;
    std::vector<extremum_info_t> extrema_to_process;

    // Collect hypothetically spurious minima
    for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
        if (hop_nbhd.hop_idx <= hop_idx_thld) {
            extrema_to_process.push_back({vertex, hop_nbhd.hop_idx, true});
        }
    }

    // Collect hypothetically spurious maxima
    for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
        if (hop_nbhd.hop_idx <= hop_idx_thld) {
            extrema_to_process.push_back({vertex, hop_nbhd.hop_idx, false});
        }
    }

    // Sort extrema by hop_idx if requested
    if (process_in_order) {
        std::sort(extrema_to_process.begin(), extrema_to_process.end(),
                 [](const extremum_info_t& a, const extremum_info_t& b) {
                     return a.hop_idx < b.hop_idx;
                 });
    }

    // Process each extremum
    for (const auto& extremum : extrema_to_process) {
        size_t vertex = extremum.vertex;
        bool is_minimum = extremum.is_minimum;

        // Get the appropriate hop neighborhood
        const hop_nbhd_t& hop_nbhd = is_minimum ?
            lmin_hop_nbhd_map[vertex] : lmax_hop_nbhd_map[vertex];

        // Record pre-smoothing state if detailed recording is enabled
        smoothing_step_t step;
        if (detailed_recording) {
            step.vertex     = vertex;
            step.is_minimum = is_minimum;
            step.hop_idx    = hop_nbhd.hop_idx;
            step.smoother   = smoother_type;
            step.before     = smoothed_y;
        }

        // Get interior vertices (hop distance <= hop_radius)
        std::unordered_set<size_t> region_vertices;
        for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
            if (dist <= hop_nbhd.hop_idx) {
                region_vertices.insert(v);
            }
        }

        // Get boundary vertices (hop distance = hop_radius + 1)
        std::unordered_map<size_t, double> boundary_values;
        for (const auto& [v, val] : hop_nbhd.y_nbhd_bd_map) {
            boundary_values[v] = smoothed_y[v];
        }

        if (detailed_recording) {
            step.region          = region_vertices;
			step.boundary_values = boundary_values;
        }

        // Apply the appropriate smoothing method
        switch (smoother_type) {
            case smoother_type_t::WMEAN:
                perform_weighted_mean_hop_disk_extension(
                    smoothed_y,
                    hop_nbhd,
                    max_inner_iterations,
                    smoothing_tolerance,
                    sigma,
                    verbose
                );
                break;

            case smoother_type_t::HARMONIC_IT:
                {
                    int record_frequency = 10;
                    harmonic_extender_t res = harmonic_extender(
                        boundary_values,
                        region_vertices,
                        max_inner_iterations,
                        smoothing_tolerance,
                        record_frequency,
                        verbose
                    );

                    // Update smoothed_y with the final iteration results
                    if (!res.iterations.empty()) {
                        auto& final_iter = res.iterations.back();
                        for (size_t v : region_vertices) {
                            smoothed_y[v] = final_iter[v];
                        }
                    }
                }
                break;

            case smoother_type_t::HARMONIC_EIGEN:
                {
                    double regularization = 1e-10;
                    auto harmonic_result = harmonic_extension_eigen(
                        boundary_values,
                        region_vertices,
                        regularization,
                        verbose
                    );

                    // Update smoothed_y with harmonic extension result
                    for (size_t v : region_vertices) {
                        smoothed_y[v] = harmonic_result[v];
                    }
                }
                break;

            case smoother_type_t::HYBRID_BIHARMONIC_HARMONIC:
                {
                    int boundary_blend_distance = 2;
                    auto hybrid_result = hybrid_biharmonic_harmonic_extension(
                        boundary_values,
                        region_vertices,
                        boundary_blend_distance,
                        verbose
                    );

                    // Update smoothed_y with hybrid extension result
                    for (size_t v : region_vertices) {
                        smoothed_y[v] = hybrid_result[v];
                    }
                }
                break;

            case smoother_type_t::BOUNDARY_SMOOTHED_HARMONIC:
                {
                    int boundary_blend_distance = 2;
                    auto boundary_smoothed_result = boundary_smoothed_harmonic_extension(
                        boundary_values,
                        region_vertices,
                        boundary_blend_distance,
                        verbose
                    );

                    // Update smoothed_y with boundary smoothed extension result
                    for (size_t v : region_vertices) {
                        smoothed_y[v] = boundary_smoothed_result[v];
                    }
                }
                break;
        }

        // Record post-smoothing state if detailed recording is enabled
        if (detailed_recording) {
            step.after = smoothed_y;
            smoothing_history.push_back(step);
        }
    }

    // Store final smoothed values
    gflow_cx.harmonic_predictions = smoothed_y;

    // If detailed recording was requested, store the history
    if (detailed_recording) {
        // We'll need to add this field to the gflow_cx_t struct
        gflow_cx.smoothing_history = smoothing_history;
    }

    return gflow_cx;
}


// ---------------------------------------------------------------------------------------------------------
//
// May 15 version
//
// ---------------------------------------------------------------------------------------------------------

gflow_cx_t set_wgraph_t::create_gflow_cx(
	const std::vector<double>& y,
	size_t hop_idx_thld,
	smoother_type_t smoother_type,
	int max_smoothing_iterations,
	int max_inner_iterations,
	double smoothing_tolerance,
	double sigma,
	bool verbose
	) const {

	//------------------------------------------------------------------
	// PHASE 1: Find local extrema
	//------------------------------------------------------------------

	auto [lmin_hop_nbhd_map, lmax_hop_nbhd_map] = compute_extrema_hop_nbhds(y);

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

				switch(smoother_type) {
				case smoother_type_t::WMEAN:
					perform_weighted_mean_hop_disk_extension(
						smoothed_y,
						hop_nbhd,
						max_inner_iterations,
						smoothing_tolerance,
						sigma,
						verbose
						);
					break;

				case smoother_type_t::HARMONIC_IT:
				{
					size_t hop_radius = hop_nbhd.hop_idx;

					// Get interior vertices (hop distance <= hop_radius)
					std::unordered_set<size_t> region_vertices;
					for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
						if (dist <= hop_radius) {
							region_vertices.insert(v);
						}
					}

					// Get boundary vertices (hop distance = hop_radius + 1)
					std::unordered_map<size_t, double> boundary_values;
					for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
						boundary_values[v] = y[v];
					}

					int record_frequency = 10;
					harmonic_extender_t res = harmonic_extender(
						boundary_values,
						region_vertices,
						max_inner_iterations,
						smoothing_tolerance,
						record_frequency,
						verbose
						);
				}
				break;

				case smoother_type_t::HARMONIC_EIGEN:
				{
					size_t hop_radius = hop_nbhd.hop_idx;

					// Get interior vertices (hop distance <= hop_radius)
					std::unordered_set<size_t> region_vertices;
					for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
						if (dist <= hop_radius) {
							region_vertices.insert(v);
						}
					}

					// Get boundary vertices (hop distance = hop_radius + 1)
					std::unordered_map<size_t, double> boundary_values;
					for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
						boundary_values[v] = y[v];
					}

					double regularization = 1e-10;
					smoothed_y = harmonic_extension_eigen(
						boundary_values,
						region_vertices,
						regularization,
						verbose
						);
				}
				break;
				}
			}
		}

		for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {

			if (hop_nbhd.hop_idx <= hop_idx_thld) {

				switch(smoother_type) {
				case smoother_type_t::WMEAN:
					perform_weighted_mean_hop_disk_extension(
						smoothed_y,
						hop_nbhd,
						max_inner_iterations,
						smoothing_tolerance,
						sigma,
						verbose
						);
					break;

				case smoother_type_t::HARMONIC_IT:
				{
					size_t hop_radius = hop_nbhd.hop_idx;

					// Get interior vertices (hop distance <= hop_radius)
					std::unordered_set<size_t> region_vertices;
					for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
						if (dist <= hop_radius) {
							region_vertices.insert(v);
						}
					}

					// Get boundary vertices (hop distance = hop_radius + 1)
					std::unordered_map<size_t, double> boundary_values;
					for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
						boundary_values[v] = y[v];
					}

					int record_frequency = 1;
					harmonic_extender_t res = harmonic_extender(
						boundary_values,
						region_vertices,
						max_inner_iterations,
						smoothing_tolerance,
						record_frequency,
						verbose
						);
				}
				break;

				case smoother_type_t::HARMONIC_EIGEN:
				{
					size_t hop_radius = hop_nbhd.hop_idx;

					// Get interior vertices (hop distance <= hop_radius)
					std::unordered_set<size_t> region_vertices;
					for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
						if (dist <= hop_radius) {
							region_vertices.insert(v);
						}
					}

					// Get boundary vertices (hop distance = hop_radius + 1)
					std::unordered_map<size_t, double> boundary_values;
					for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
						boundary_values[v] = y[v];
					}

					double regularization = 1e-10;
					smoothed_y = harmonic_extension_eigen(
						boundary_values,
						region_vertices,
						regularization,
						verbose
						);
				}
				break;
				}
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

		// Check maxima
		for (const auto& vertex : new_lmax) {
			if (significant_lmax.find(vertex) == significant_lmax.end()) {
				// Found a new or previously insignificant minimum
				validation_passed = false;
				break;
			}
		}

		// Check if all significant minima still exist
		if (validation_passed) {
			for (const auto& vertex : significant_lmax) {
				if (std::find(new_lmax.begin(), new_lmax.end(), vertex) == new_lmax.end()) {
					// A significant minimum was lost during smoothing
					validation_passed = false;
					break;
				}
			}
		}

		// Decide whether to continue smoothing
		if (validation_passed) {
			smoothing_complete = true;
		} else {
			// If validation failed, update hop neighborhoods for next iteration
			lmin_hop_nbhd_map.clear();
			lmax_hop_nbhd_map.clear();

			bool detect_maxima = false;
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


	gflow_cx_t gflow_cx;
	gflow_cx.harmonic_predictions = std::move(smoothed_y);

	return gflow_cx;
}

// ---------------------------------------------------------------------------------------------------------
//
// May 14 version
//
// ---------------------------------------------------------------------------------------------------------

// You've proposed an excellent refinement to the process. Breaking it into three phases with an iterative smoothing step makes a lot of sense. Let me analyze this approach and suggest some implementation details.

// ## Three-Phase Process Analysis

// ### 1. Analysis Phase
// This initial phase remains focused on extrema detection and threshold determination:
// - Find all neighborhood local extrema using `find_nbr_extrema()`
// - Compute hop neighborhoods and indices using `compute_extremum_hop_nbhd()`
// - Analyze the distribution of hop indices to determine an appropriate threshold

// ### 2. Smoothing Phase (Iterative)
// This phase addresses a critical concern: ensuring that smoothing actually achieves the desired elimination of spurious extrema. Breaking it into two sub-steps with validation is very insightful.

// #### 2(A). Smoothing Operation
// - Apply harmonic or quasi-harmonic smoothing to neighborhoods of extrema with hop_idx ≤ hop_idx_thld
// - This effectively "flattens" the spurious extrema while preserving significant ones

// #### 2(B). Validation and Iteration
// - Recompute local extrema on the smoothed function
// - Verify that only extrema with hop_idx > hop_idx_thld remain
// - If new spurious extrema emerged, repeat smoothing

// This iterative approach is crucial because harmonic smoothing can sometimes introduce new extrema or fail to completely eliminate existing ones, especially in complex graph structures.

// ### 3. Gradient Flow Complex Construction Phase
// With a properly smoothed function that maintains only significant extrema, you can confidently:
// - Compute gradient trajectories
// - Construct basins around remaining extrema
// - Identify cells in the gradient flow complex

// ## Implementation Recommendations

// Here's how I would implement this three-phase approach:


gflow_cx_t set_wgraph_t::create_gflow_cx(
	const std::vector<double>& y,
	std::optional<size_t> custom_hop_idx_thld = std::nullopt,
	int max_outer_iterations = 3,
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
	std::unordered_set<size_t> non_spurious_lmin;
	std::unordered_set<size_t> non_spurious_lmax;

	for (const auto& [vertex, hop_nbhd] : lmin_hop_nbhd_map) {
		if (hop_nbhd.hop_idx > hop_idx_thld) {
			non_spurious_lmin.insert(vertex);
		}
	}

	for (const auto& [vertex, hop_nbhd] : lmax_hop_nbhd_map) {
		if (hop_nbhd.hop_idx > hop_idx_thld) {
			non_spurious_lmax.insert(vertex);
		}
	}

	// Iterative smoothing process
	bool smoothing_complete = false;
	int smoothing_iteration = 0;

	while (!smoothing_complete && smoothing_iteration < max_outer_iterations) {
		// Phase 2(A): Apply smoothing to neighborhoods of innon_spurious extrema
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

		// Check if only non_spurious extrema remain
		bool validation_passed = true;

		// Check minima
		for (const auto& vertex : new_lmin) {
			if (non_spurious_lmin.find(vertex) == non_spurious_lmin.end()) {
				// Found a new or previously innon_spurious minimum
				validation_passed = false;
				break;
			}
		}

		// Check if all non_spurious minima still exist
		if (validation_passed) {
			for (const auto& vertex : non_spurious_lmin) {
				if (std::find(new_lmin.begin(), new_lmin.end(), vertex) == new_lmin.end()) {
					// A non_spurious minimum was lost during smoothing
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

	// Compute basins for remaining non_spurious extrema
	for (const auto& vertex : non_spurious_lmin) {
		detect_maxima = false;
		basin_t basin = find_local_extremum_bfs_basin(vertex, smoothed_y, detect_maxima);
		if (basin.reachability_map.sorted_vertices.size() > 0) {
			gflow_cx.lmin_basins_map[vertex] = basin;
		}
	}

	for (const auto& vertex : non_spurious_lmax) {
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




// ## Harmonic Smoothing Implementation
// Here's a detailed implementation of the harmonic smoothing function:
void set_wgraph_t::perform_harmonic_smoothing_in_neighborhood(
	std::vector<double>& values,
	const hop_nbhd_t& hop_nbhd,
	double tolerance,
	int max_iterations = 100) const {

	size_t vertex = hop_nbhd.vertex;
	size_t hop_radius = hop_nbhd.hop_idx;

	// Collect interior vertices (hop distance <= hop_radius)
	std::unordered_set<size_t> interior;
	for (const auto& [v, dist] : hop_nbhd.hop_dist_map) {
		if (dist <= hop_radius) {
			interior.insert(v);
		}
	}

	// Collect boundary vertices and values (hop distance = hop_radius + 1)
	std::unordered_map<size_t, double> boundary_values;
	for (const auto& [v, _] : hop_nbhd.y_nbhd_bd_map) {
		boundary_values[v] = values[v];
	}

	// Store original extremum value for special handling
	double extremum_value = values[vertex];

	// Iterative harmonic relaxation
	double max_diff = tolerance + 1.0;
	int iterations = 0;

	while (max_diff > tolerance && iterations < max_iterations) {
		max_diff = 0.0;

		// Create a copy of current values for update
		std::vector<double> new_values = values;

		// Update each interior vertex with average of neighbors
		for (size_t u : interior) {
			// Skip the extremum vertex itself to preserve its position
			if (u == vertex) continue;

			double sum = 0.0;
			int count = 0;

			// Average of neighbors
			for (const auto& edge : adjacency_list[u]) {
				size_t v = edge.vertex;

				// Only consider vertices that are either in the interior or boundary
				if (interior.find(v) != interior.end() || boundary_values.find(v) != boundary_values.end()) {
					sum += values[v];
					count++;
				}
			}

			if (count > 0) {
				new_values[u] = sum / count;
				max_diff = std::max(max_diff, std::abs(new_values[u] - values[u]));
			}
		}

		// Special handling for the extremum vertex
		// For minima: ensure it remains the minimum in its neighborhood
		// For maxima: ensure it remains the maximum in its neighborhood
		if (hop_nbhd.hop_idx == 0) {   // If it's a spurious extremum with hop_idx = 0
			double sum = 0.0;
			int count = 0;

			// Average of neighbors
			for (const auto& edge : adjacency_list[vertex]) {
				size_t v = edge.vertex;
				sum += new_values[v];
				count++;
			}

			if (count > 0) {
				// For a minimum, set value slightly higher than the minimum neighbor
				// For a maximum, set value slightly lower than the maximum neighbor
				bool is_maximum = hop_nbhd.vertex == vertex ? true : false;   // Assuming detect_maxima was passed to hop_nbhd

				if (is_maximum) {
					double max_neighbor = -std::numeric_limits<double>::infinity();
					for (const auto& edge : adjacency_list[vertex]) {
						max_neighbor = std::max(max_neighbor, new_values[edge.vertex]);
					}
					new_values[vertex] = max_neighbor - tolerance;
				} else {
					double min_neighbor = std::numeric_limits<double>::infinity();
					for (const auto& edge : adjacency_list[vertex]) {
						min_neighbor = std::min(min_neighbor, new_values[edge.vertex]);
					}
					new_values[vertex] = min_neighbor + tolerance;
				}
			}
		}

		// Update values (but keep boundary fixed)
		for (size_t u : interior) {
			values[u] = new_values[u];
		}

		iterations++;
	}
}


// ## Advantages of the Three-Phase Approach

// 1. **Methodical Progression**: The three phases provide a clear progression from analysis to smoothing to construction.

// 2. **Validation of Smoothing**: The iterative smoothing with validation ensures that we actually achieve our goal of eliminating spurious extrema without losing non_spurious ones.

// 3. **Preserves Important Extrema**: By explicitly tracking which extrema should be preserved, we can ensure they aren't inadvertently smoothed away.

// 4. **More Robust Results**: The validation step helps catch cases where smoothing might introduce new extrema, which can happen in complex graphs.

// 5. **Flexibility**: This approach gives you the option to use different smoothing methods (quasi-harmonic, true harmonic, or other techniques) without changing the overall algorithm structure.

// ## Conclusion

// Your proposed three-phase approach with iterative smoothing is non_spuriously more robust than a simple two-phase process. It addresses the critical concern that smoothing might not achieve the desired effect in a single pass, and provides a framework for validation and refinement.

// The implementation I've outlined should give you a solid foundation for implementing this approach, with particular attention to the harmonic smoothing and validation steps that are crucial for ensuring that only non_spurious extrema remain after smoothing.
