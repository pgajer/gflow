#include <string>       // For std::string
#include <fstream>      // For file I/O (std::ofstream)
#include <filesystem>   // For std::filesystem operations
#include <algorithm>    // For std::sort
#include <map>          // For std::map used in basin_t
#include <unordered_map>// For std::unordered_map in basin_cx_t
#include <utility>      // For std::pair
#include <vector>       // For std::vector
#include <cstddef>      // For size_t

#include "basin.hpp"    // For basin_t and basin_cx_t definitions
#include "invalid_vertex.hpp" // For INVALID_VERTEX constant
#include "error_utils.h"      // For REPORT_ERROR macro
#include "set_wgraph.hpp"

/**
 * @brief Absorbs one basin into another during topological simplification
 *
 * @details This function performs the critical operation of basin absorption during
 *          the construction of a gradient flow basin complex. When a cancellation pair
 *          of basins is identified below the significance threshold, one basin is
 *          absorbed by its "cousin" basin. This function:
 *
 *          1. Merges vertices from the absorbed basin into the absorbing basin
 *          2. Updates distance calculations in the reachability map
 *          3. Removes the absorbed basin from the basin maps
 *          4. Updates boundary relationships
 *          5. Updates min monotonicity span and vertex using new basin boundary
 *          6. Performs harmonic repair of function values in the absorbed region
 *
 * The absorption process is key to eliminating topologically insignificant
 * features while preserving the geometric structure of the function.
 *
 * @param absorbing_basin    The basin that will absorb another basin
 * @param absorbed_basin     The basin being absorbed
 * @param y_range            Range of function values (max - min) for relative calculations
 * @param graph              Reference to the graph for geometric operations
 *
 * @post The absorbed basin is removed from the basin maps, its vertices are incorporated
 *       into the absorbing basin, and function values are harmonically repaired.
 *
 * @see find_local_extremum_basin
 * @see create_basin_cx
 * @see performHarmonicRepair
 */
void basin_cx_t::absorb_basin(
	basin_t& absorbing_basin,
	const basin_t& absorbed_basin,
	const std::vector<double>& y,
	const set_wgraph_t& graph
    ) {

	double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    double y_range = y_max - y_min;

	// 1. Get the extremum vertices
	size_t absorbing_extremum = absorbing_basin.extremum_vertex;
	size_t absorbed_extremum = absorbed_basin.extremum_vertex;

	// 2. Determine which maps to update
	auto& absorbing_basin_map = absorbing_basin.is_maximum ?
		lmax_basins_map : lmin_basins_map;
	auto& absorbed_basin_map = absorbed_basin.is_maximum ?
		lmax_basins_map : lmin_basins_map;

	// 3. Get mutable reference to the absorbing basin
	basin_t& absorbing_basin_ref = absorbing_basin_map[absorbing_extremum];

	// 4. Update the absorbing basin with vertices from the absorbed basin
	for (const auto& v_info : absorbed_basin.reachability_map.sorted_vertices) {
		// Add vertices if not already in the absorbing basin
		if (absorbing_basin_ref.reachability_map.distances.find(v_info.vertex) ==
			absorbing_basin_ref.reachability_map.distances.end()) {

			// Calculate distance from absorbing extremum
			double dist = absorbing_basin_ref.reachability_map.distances[absorbed_extremum] +
				v_info.distance;

			// Add to absorbing basin's reachability map
			absorbing_basin_ref.reachability_map.distances[v_info.vertex] = dist;
			absorbing_basin_ref.reachability_map.sorted_vertices.push_back({v_info.vertex, dist});

			// Set predecessor (simplified - could be improved with actual path reconstruction)
			absorbing_basin_ref.reachability_map.predecessors[v_info.vertex] = absorbed_extremum;
		}
	}

	// 5. Re-sort vertices by distance
	std::sort(
		absorbing_basin_ref.reachability_map.sorted_vertices.begin(),
		absorbing_basin_ref.reachability_map.sorted_vertices.end(),
		[](const vertex_info_t& a, const vertex_info_t& b) {
			return a.distance > b.distance;
		}
        );

	// 6. Remove the absorbed basin from its map
	absorbed_basin_map.erase(absorbed_extremum);

	// 7. Update the absorbing basin's boundary information We need to find the
	// new boundary of the absorbing basin after the absorbing process; this
	// needs to be done de novo as there is no other way to figure out what
	// vertices are at the boundary after the merge
	// The vertex u is at the boundary of absorbing_basin either if it is of
	// degree 1 or it has a neighbor that is not in the basin
	for (auto& v_info : absorbing_basin.reachability_map.sorted_vertices) {

		auto& u = v_info.vertex;
		bool is_boundary = false;
		if (graph.adjacency_list[u].size() == 1) {
			is_boundary = true;
		}  else {
			for (const auto& edge : graph.adjacency_list[u]) {
				size_t v = edge.vertex;
				if (absorbing_basin.reachability_map.distances.find(v) == absorbing_basin.reachability_map.distances.end()) {
					is_boundary = true;
					break;
				}
			}
		}

        // If this is a boundary vertex, add it to the boundary map
        if (is_boundary) {
			// computing correct distance for boundary vertices
			double dist = graph.bidirectional_dijkstra(absorbing_extremum, u);
            absorbing_basin.boundary_vertices_map[dist] = u;
			absorbing_basin.boundary_monotonicity_spans_map[u] = std::abs(y[absorbing_extremum] - y[u]);
		}
	}

	// 8. Update min_monotonicity_span, rel_min_monotonicity_span and min_span_vertex
	double min_span = INFINITY;
	size_t min_span_vertex = INVALID_VERTEX;
	for (const auto& [u, span] : absorbing_basin.boundary_monotonicity_spans_map) {
		if (
			( // only use for the min_span calculation boundary vertices that come from local extrema of the opposite orientation than absorbing_basin.is_maximum
				(absorbing_basin.is_maximum && lmin_basins_map.find(u) != lmin_basins_map.end())
				|| (!absorbing_basin.is_maximum && lmax_basins_map.find(u) != lmax_basins_map.end())
				)
			&& span < min_span
			) {
			min_span = span;
			min_span_vertex = u;
		}
	}

	absorbing_basin.min_span_vertex           = min_span_vertex;
	absorbing_basin.min_monotonicity_span     = min_span;
	absorbing_basin.rel_min_monotonicity_span = (y_range > 0.0) ? min_span / y_range : 0.0;

	// 9. Perform harmonic repair
	graph.perform_harmonic_repair(harmonic_predictions, absorbing_basin, absorbed_basin);
}

/**
 * @brief Writes comprehensive basin data to files with simplified format
 *
 * @param out_dir Output directory path
 * @param prefix Optional prefix for output filenames
 */
void basin_cx_t::write_basins_map(
    const std::string& out_dir,
    const std::string& prefix
) const {
    // Fixed parameters
    constexpr size_t offset = 1;  // Always use 1-based indexing (for R)
    const std::string delimiter = ",";  // Always use comma as delimiter

    // Ensure output directory exists
    if (!std::filesystem::exists(out_dir)) {
        if (!std::filesystem::create_directories(out_dir)) {
            REPORT_ERROR("ERROR: Failed to create output directory: %s\n", out_dir.c_str());
            return;
        }
    }

    // Generate filenames
    std::string lmin_file = out_dir + "/" + prefix + "lmin_basins.txt";
    std::string lmax_file = out_dir + "/" + prefix + "lmax_basins.txt";
    std::string lmin_boundaries_file = out_dir + "/" + prefix + "lmin_boundaries.txt";
    std::string lmax_boundaries_file = out_dir + "/" + prefix + "lmax_boundaries.txt";
    std::string summary_file = out_dir + "/" + prefix + "basin_summary.txt";

    // Write local minima basins
    {
        std::ofstream lmin_out(lmin_file);
        if (!lmin_out) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmin_file.c_str());
            return;
        }

        // Extract and sort minima indices
        std::vector<size_t> lmins;
        for (const auto& [vertex, _] : lmin_basins_map) {
            lmins.push_back(vertex);
        }
        std::sort(lmins.begin(), lmins.end());

        // Write header and sorted list of minima
        lmin_out << "## Local Minima\n";
        lmin_out << "lmins <- c(";
        for (size_t i = 0; i < lmins.size(); ++i) {
            lmin_out << (lmins[i] + offset);
            if (i < lmins.size() - 1) {
                lmin_out << delimiter;
            }
        }
        lmin_out << ") # vertex indices of local minima sorted in ascending order\n";

        // Write basin list header
        lmin_out << "## Local Minima Basins list\n";
        lmin_out << "lmin.basins <- list()\n";

        // Write each basin
        for (const auto& min_idx : lmins) {
            const basin_t& basin = lmin_basins_map.at(min_idx);

            // Collect and sort basin vertices
            std::vector<size_t> basin_vertices;
            for (const auto& v_info : basin.reachability_map.sorted_vertices) {
                basin_vertices.push_back(v_info.vertex);
            }
            std::sort(basin_vertices.begin(), basin_vertices.end());

            // Write basin entry with metadata comment
            lmin_out << "lmin.basins[[\"" << (min_idx + offset) << "\"]] <- c(";
            for (size_t i = 0; i < basin_vertices.size(); ++i) {
                lmin_out << (basin_vertices[i] + offset);
                if (i < basin_vertices.size() - 1) {
                    lmin_out << delimiter;
                }
            }
            lmin_out << ") # value: " << basin.value
                     << " min_span: " << basin.min_monotonicity_span
                     << " rel_min_span: " << basin.rel_min_monotonicity_span
                     << " min_span_vertex: " << (basin.min_span_vertex == INVALID_VERTEX ?
                                                -1 : basin.min_span_vertex + offset)
                     << "\n";
        }

        lmin_out.close();
        Rprintf("Written local minima basins to: %s\n", lmin_file.c_str());
    }

    // Write local maxima basins
    {
        std::ofstream lmax_out(lmax_file);
        if (!lmax_out) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmax_file.c_str());
            return;
        }

        // Extract and sort maxima indices
        std::vector<size_t> lmaxs;
        for (const auto& [vertex, _] : lmax_basins_map) {
            lmaxs.push_back(vertex);
        }
        std::sort(lmaxs.begin(), lmaxs.end());

        // Write header and sorted list of maxima
        lmax_out << "## Local Maxima\n";
        lmax_out << "lmaxs <- c(";
        for (size_t i = 0; i < lmaxs.size(); ++i) {
            lmax_out << (lmaxs[i] + offset);
            if (i < lmaxs.size() - 1) {
                lmax_out << delimiter;
            }
        }
        lmax_out << ") # vertex indices of local maxima sorted in ascending order\n";

        // Write basin list header
        lmax_out << "## Local Maxima Basins list\n";
        lmax_out << "lmax.basins <- list()\n";

        // Write each basin
        for (const auto& max_idx : lmaxs) {
            const basin_t& basin = lmax_basins_map.at(max_idx);

            // Collect and sort basin vertices
            std::vector<size_t> basin_vertices;
            for (const auto& v_info : basin.reachability_map.sorted_vertices) {
                basin_vertices.push_back(v_info.vertex);
            }
            std::sort(basin_vertices.begin(), basin_vertices.end());

            // Write basin entry with metadata comment
            lmax_out << "lmax.basins[[\"" << (max_idx + offset) << "\"]] <- c(";
            for (size_t i = 0; i < basin_vertices.size(); ++i) {
                lmax_out << (basin_vertices[i] + offset);
                if (i < basin_vertices.size() - 1) {
                    lmax_out << delimiter;
                }
            }
            lmax_out << ") # value: " << basin.value
                     << " min_span: " << basin.min_monotonicity_span
                     << " rel_min_span: " << basin.rel_min_monotonicity_span
                     << " min_span_vertex: " << (basin.min_span_vertex == INVALID_VERTEX ?
                                                -1 : basin.min_span_vertex + offset)
                     << "\n";
        }

        lmax_out.close();
        Rprintf("Written local maxima basins to: %s\n", lmax_file.c_str());
    }

    // Write local minima boundary information (keep the original format)
    {
        std::ofstream lmin_boundaries(lmin_boundaries_file);
        if (!lmin_boundaries) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmin_boundaries_file.c_str());
            return;
        }

        lmin_boundaries << "# Local Minima Basin Boundaries\n";
        lmin_boundaries << "# Format: extremum_vertex" << delimiter << "num_boundary_vertices" << delimiter << "[boundary_vertex span_value]...\n";

        for (const auto& [vertex, basin] : lmin_basins_map) {
            lmin_boundaries << (vertex + offset) << delimiter
                          << basin.boundary_vertices_map.size();

            // Output boundary vertices sorted by distance
            for (const auto& [distance, bv] : basin.boundary_vertices_map) {
                lmin_boundaries << delimiter << (bv + offset);
            }

            lmin_boundaries << "\n# Boundary monotonicity spans for extremum " << (vertex + offset) << ":\n";
            lmin_boundaries << "# vertex" << delimiter << "monotonicity_span" << delimiter << "is_extremum" << delimiter << "extremum_type\n";

            // Output detailed boundary monotonicity spans
            for (const auto& [bv, span] : basin.boundary_monotonicity_spans_map) {
                bool is_extremum = false;
                bool is_maximum = false;

                // Check if this boundary vertex is an extremum
                if (lmax_basins_map.find(bv) != lmax_basins_map.end()) {
                    is_extremum = true;
                    is_maximum = true;
                } else if (lmin_basins_map.find(bv) != lmin_basins_map.end()) {
                    is_extremum = true;
                    is_maximum = false;
                }

                lmin_boundaries << (bv + offset) << delimiter << span << delimiter
                              << (is_extremum ? 1 : 0) << delimiter
                              << (is_maximum ? "max" : "min") << "\n";
            }

            lmin_boundaries << "\n";
        }

        lmin_boundaries.close();
        Rprintf("Written local minima boundary information to: %s\n", lmin_boundaries_file.c_str());
    }

    // Write local maxima boundary information (keep the original format)
    {
        std::ofstream lmax_boundaries(lmax_boundaries_file);
        if (!lmax_boundaries) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", lmax_boundaries_file.c_str());
            return;
        }

        lmax_boundaries << "# Local Maxima Basin Boundaries\n";
        lmax_boundaries << "# Format: extremum_vertex" << delimiter << "num_boundary_vertices" << delimiter << "[boundary_vertex span_value]...\n";

        for (const auto& [vertex, basin] : lmax_basins_map) {
            lmax_boundaries << (vertex + offset) << delimiter
                          << basin.boundary_vertices_map.size();

            // Output boundary vertices sorted by distance
            for (const auto& [distance, bv] : basin.boundary_vertices_map) {
                lmax_boundaries << delimiter << (bv + offset);
            }

            lmax_boundaries << "\n# Boundary monotonicity spans for extremum " << (vertex + offset) << ":\n";
            lmax_boundaries << "# vertex" << delimiter << "monotonicity_span" << delimiter << "is_extremum" << delimiter << "extremum_type\n";

            // Output detailed boundary monotonicity spans
            for (const auto& [bv, span] : basin.boundary_monotonicity_spans_map) {
                bool is_extremum = false;
                bool is_maximum = false;

                // Check if this boundary vertex is an extremum
                if (lmax_basins_map.find(bv) != lmax_basins_map.end()) {
                    is_extremum = true;
                    is_maximum = true;
                } else if (lmin_basins_map.find(bv) != lmin_basins_map.end()) {
                    is_extremum = true;
                    is_maximum = false;
                }

                lmax_boundaries << (bv + offset) << delimiter << span << delimiter
                              << (is_extremum ? 1 : 0) << delimiter
                              << (is_maximum ? "max" : "min") << "\n";
            }

            lmax_boundaries << "\n";
        }

        lmax_boundaries.close();
        Rprintf("Written local maxima boundary information to: %s\n", lmax_boundaries_file.c_str());
    }

    // Write a summary file with basin statistics and cancellation pairs
    {
        std::ofstream summary_out(summary_file);
        if (!summary_out) {
            REPORT_ERROR("ERROR: Could not open file for writing: %s\n", summary_file.c_str());
            return;
        }

        summary_out << "# Basin Complex Summary\n";
        summary_out << "Number of local minima basins: " << lmin_basins_map.size() << "\n";
        summary_out << "Number of local maxima basins: " << lmax_basins_map.size() << "\n";
        summary_out << "Number of cancellation pairs: " << cancellation_pairs.size() << "\n\n";

        // Write all cancellation pairs with more details
        summary_out << "# Cancellation Pairs (lmin_vertex, lmax_vertex)\n";
        summary_out << "# Format: lmin_vertex" << delimiter << "lmax_vertex" << delimiter << "lmin_span" << delimiter << "lmax_span\n";
        for (const auto& [lmin, lmax] : cancellation_pairs) {
            double lmin_span = lmin_basins_map.count(lmin) > 0 ?
                              lmin_basins_map.at(lmin).rel_min_monotonicity_span : -1.0;
            double lmax_span = lmax_basins_map.count(lmax) > 0 ?
                              lmax_basins_map.at(lmax).rel_min_monotonicity_span : -1.0;

            summary_out << (lmin + offset) << delimiter << (lmax + offset) << delimiter
                       << lmin_span << delimiter << lmax_span << "\n";
        }

        summary_out.close();
        Rprintf("Written basin summary to: %s\n", summary_file.c_str());
    }
}

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

basin_t* basin_cx_t::process_cancellation_pair(
	basin_t& A_basin,
	basin_t& B_basin,
	const std::vector<double>& y,
	const set_wgraph_t& graph
	) {

#define DEBUG__process_cancellation_pair 0

	size_t A_basin_extremum = A_basin.extremum_vertex;
	(void)A_basin_extremum;

	// Track processed extrema to avoid re-processing
    std::unordered_set<size_t> processed_extrema;

	bool found_B_basin = false;
	if (A_basin.is_maximum) {
		// If A is a maximum, look for minimum basin at min_span_vertex
		auto it = lmin_basins_map.find(A_basin.min_span_vertex);
		if (it != lmin_basins_map.end()) {
			B_basin = it->second;
			found_B_basin = true;
		}
	} else {
		// If A is a minimum, look for maximum basin at min_span_vertex
		auto it = lmax_basins_map.find(A_basin.min_span_vertex);
		if (it != lmax_basins_map.end()) {
			B_basin = it->second;
			found_B_basin = true;
		}
	}

	if (!found_B_basin) {
#if DEBUG__process_cancellation_pair
		// I believe this should never happen, so I want to stop the process if it happens to understand that case
		REPORT_ERROR("A_basin %zu basin with rel_min_monotonicity_span %.3f "
					 "does not have B_basin cancellation partner\n",
					 A_basin_extremum,
					 A_basin.rel_min_monotonicity_span);
#endif
	}

	// Check if B's min_span_vertex points back to A
	if (B_basin.min_span_vertex != A_basin.extremum_vertex) {
#if DEBUG__process_cancellation_pair
		// I believe this should never happen, so I want to stop the process if it happens to understand that case
		Rprintf("\nA_basin_extremum: %zu\n"
				"A_basin.extremum_vertex: %zu\n"
				"A_basin.min_span_vertex: %zu\n"
				"B_basin.extremum_vertex: %zu\n"
				"B_basin.min_span_vertex: %zu\n\n",
				A_basin_extremum + 1,
				A_basin.extremum_vertex + 1,
				A_basin.min_span_vertex + 1,
				B_basin.extremum_vertex + 1,
				B_basin.min_span_vertex + 1);

		Rprintf("A_basin.rel_min_monotonicity_span: %.3f\n"
				"rel_min_monotonicity_span_thld: %.3f\n"
				"B_basin.rel_min_monotonicity_span: %.3f\n",
				A_basin.rel_min_monotonicity_span,
				rel_min_monotonicity_span_thld,
				B_basin.rel_min_monotonicity_span);
		REPORT_ERROR("lmin B_basin %zu basin's min_span_vertex: %zu is not equal to A_basin.extremum_vertex: %zu\n",
					 B_basin.extremum_vertex + 1,
					 B_basin.min_span_vertex + 1,
					 A_basin.extremum_vertex + 1);
#endif
	}

	// Found a valid cancellation pair - record it
	cancellation_pairs.emplace_back(
		A_basin.is_maximum ? B_basin.extremum_vertex : A_basin.extremum_vertex,
		A_basin.is_maximum ? A_basin.extremum_vertex : B_basin.extremum_vertex
		);

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
			auto it = lmin_basins_map.find(u);
			if (it != lmin_basins_map.end() &&
				processed_extrema.count(u) == 0) {
				A_cousin_basin = it->second;
				A_cousin_distance = d;
				found_A_cousin_basin = true;
				break;
			}
		} else {
			auto it = lmax_basins_map.find(u);
			if (it != lmax_basins_map.end() &&
				processed_extrema.count(u) == 0) {
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
			auto it = lmin_basins_map.find(u);
			if (it != lmin_basins_map.end() &&
				processed_extrema.count(u) == 0) {
				B_cousin_basin = it->second;
				B_cousin_distance = d;
				found_B_cousin_basin = true;
				break;
			}
		} else {
			auto it = lmax_basins_map.find(u);
			if (it != lmax_basins_map.end() &&
				processed_extrema.count(u) == 0) {
				B_cousin_basin = it->second;
				B_cousin_distance = d;
				found_B_cousin_basin = true;
				break;
			}
		}
	}

	// 3) Determine which basin to absorb
	size_t absorbing_extremum = INVALID_VERTEX;
	size_t absorbed_extremum = INVALID_VERTEX;
	basin_t* absorbing_basin_ptr = nullptr;
	basin_t* absorbed_basin_ptr = nullptr;

	if (found_A_cousin_basin && found_B_cousin_basin) {
		if (A_cousin_distance < B_cousin_distance) {
			absorbing_extremum = A_cousin_basin.extremum_vertex;
			absorbed_extremum = A_basin.extremum_vertex;
			absorbing_basin_ptr = A_cousin_basin.is_maximum ?
				&lmax_basins_map[absorbing_extremum] :
				&lmin_basins_map[absorbing_extremum];
			absorbed_basin_ptr = A_basin.is_maximum ?
				&lmax_basins_map[absorbed_extremum] :
				&lmin_basins_map[absorbed_extremum];
		} else {
			absorbing_extremum = B_cousin_basin.extremum_vertex;
			absorbed_extremum = B_basin.extremum_vertex;
			absorbing_basin_ptr = B_cousin_basin.is_maximum ?
				&lmax_basins_map[absorbing_extremum] :
				&lmin_basins_map[absorbing_extremum];
			absorbed_basin_ptr = B_basin.is_maximum ?
				&lmax_basins_map[absorbed_extremum] :
				&lmin_basins_map[absorbed_extremum];
		}
	} else if (found_A_cousin_basin) {
		absorbing_extremum = A_cousin_basin.extremum_vertex;
		absorbed_extremum = A_basin.extremum_vertex;
		absorbing_basin_ptr = A_cousin_basin.is_maximum ?
			&lmax_basins_map[absorbing_extremum] :
			&lmin_basins_map[absorbing_extremum];
		absorbed_basin_ptr = A_basin.is_maximum ?
			&lmax_basins_map[absorbed_extremum] :
			&lmin_basins_map[absorbed_extremum];
	} else if (found_B_cousin_basin) {
		absorbing_extremum = B_cousin_basin.extremum_vertex;
		absorbed_extremum = B_basin.extremum_vertex;
		absorbing_basin_ptr = B_cousin_basin.is_maximum ?
			&lmax_basins_map[absorbing_extremum] :
			&lmin_basins_map[absorbing_extremum];
		absorbed_basin_ptr = B_basin.is_maximum ?
			&lmax_basins_map[absorbed_extremum] :
			&lmin_basins_map[absorbed_extremum];
	} else {
		// No cousin basin found - we would need a default strategy here
		// For now, just skip this pair
#if DEBUG__process_cancellation_pair
		// I believe this should never happen, so I want to stop the process if it happens to understand that case
		REPORT_ERROR("We didn't find A_cousin_basin for A_basin_extremum: %zu AND "
					 "We didn't find B_cousin_basin for B_basin.extremum_vertex: %zu\n",
					 A_basin_extremum + 1,
					 B_basin.extremum_vertex + 1);
#endif
	}

	// Record the extrema for repair
	repair_lextr[{A_basin.extremum_vertex, B_basin.extremum_vertex}] = absorbing_extremum;

	// Mark the absorbed basin's extremum as processed
	processed_extrema.insert(absorbed_extremum);

	// 4) Update the non-absorbed partner of the cancellation pair
	// If A was absorbed, update B; if B was absorbed, update A
	size_t partner_extremum = (absorbed_extremum == A_basin.extremum_vertex) ?
		B_basin.extremum_vertex : A_basin.extremum_vertex;
	(void)partner_extremum;

#if 0
	basin_t* partner_basin_ptr = nullptr;
	if (lmin_basins_map.count(partner_extremum) > 0) {
		partner_basin_ptr = &lmin_basins_map[partner_extremum];
	} else {
		partner_basin_ptr = &lmax_basins_map[partner_extremum];
	}
#endif

	// Remove absorbed extremum from partner's boundary mapping
	// this has to be totally rewritten
	// partner_basin_ptr->boundary_monotonicity_spans_map.erase(absorbed_extremum);
	// Update partner's min_span and recalculate
	// partner_basin_ptr->update_min_span(y_range);
	// basin_queue.emplace(partner_extremum,
	// 					partner_basin_ptr->rel_min_monotonicity_span,
	// 					partner_basin_ptr->min_span_vertex,
	// 					partner_basin_ptr->is_maximum);

	// 5) Perform basin absorption and harmonic repair
	absorb_basin(*absorbing_basin_ptr, *absorbed_basin_ptr, y, graph);
	// basin_queue.emplace(absorbing_extremum,
	// 					absorbing_basin_ptr->rel_min_monotonicity_span,
	// 					absorbing_basin_ptr->min_span_vertex,
	// 					absorbing_basin_ptr->is_maximum);

	return absorbing_basin_ptr;
}



// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

basin_t* basin_cx_t::process_single_basin(
	basin_t& A_basin,
	const std::vector<double>& y,
	const set_wgraph_t& graph
	) {

	(void)A_basin;
	(void)y;
	(void)graph;

	basin_t* absorbing_basin_ptr = nullptr;

	return absorbing_basin_ptr;
}
