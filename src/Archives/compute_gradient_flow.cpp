#include "compute_gradient_flow.h"

gradient_flow_t set_wgraph_t::compute_gradient_flow(
    std::vector<double>& y,
    std::vector<double>& scale) const {

    #define DEBUG__compute_gradient_flow 1

    gradient_flow_t result;
    result.scale = std::move(scale);

    #if 1
    double noise_magnitude = 0.01;
    if (break_duplicate_values(y, noise_magnitude)) {
        // print_vect(y,"y - after adding noise");
        result.messages = "Small noise added to break ties in duplicate values";
        // REPORT_WARNING("Warning: small noise added to break ties in duplicate values\n");
    }
    #endif

    // Clear cache and initialize unprocessed vertices
    paths_cache.clear();
    unprocessed_vertices.clear();
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        unprocessed_vertices.insert(i);
    }

    // Find local extrema
    for (size_t v : unprocessed_vertices) {
        auto shortest_paths = find_graph_paths_within_radius(v, result.scale[v]);
        auto extremum = check_local_extremum(v, shortest_paths, y);
        if (extremum) {
            result.local_extrema[extremum->first] = extremum->second;
        }
    }

    #if DEBUG__compute_gradient_flow
    {
        Rprintf("\nIn gradient_flow_t set_wgraph_t::compute_gradient_flow()\n");
        Rprintf("vertex\tis_lmax\n");
        for (const auto& [v, is_lmax] : result.local_extrema) {
            Rprintf("%zu\t%d\n", v, (int)is_lmax);
        }
        Rprintf("\n");
        //print_vect_pair(result.local_extrema, "local_extrema", true);
        // for (size_t i = 0; i < result.local_extrema.size(); ++i) {
        //     size_t v = result.local_extrema[i].first;
        //     auto shortest_paths = find_graph_paths_within_radius(v, result.scale[v]);
        //     Rprintf("vertex: %zu\n", v);
        //     print_uset(shortest_paths.reachable_vertices,"shortest_paths.reachable_vertices");
        // }
    }
    #endif

    // Construct trajectories
    #if DEBUG__compute_gradient_flow
    size_t while_counter = 0;
    #endif
    while (!unprocessed_vertices.empty()) {
        size_t v = *unprocessed_vertices.begin();

        gradient_trajectory_t ascending_traj  = construct_trajectory(v, true, result.scale, y, result.local_extrema);
        gradient_trajectory_t descending_traj = construct_trajectory(v, false, result.scale, y, result.local_extrema);

        // In extreme cases we we discover trajectories connecting lmin with lmin and lmax with lmax

#if DEBUG__compute_gradient_flow
        Rprintf("After construct_trajectory() calls\n");
        Rprintf("v: %zu\n",v);
        print_vect(ascending_traj.path,"ascending_traj.path");
        Rprintf("ascending_traj.ends_at_critical: %d\tascending_traj.ends_at_lmax: %d\n",
                (int)ascending_traj.ends_at_critical, (int)ascending_traj.ends_at_lmax);
        print_vect(descending_traj.path,"descending_traj.path");
        Rprintf("descending_traj.ends_at_critical: %d\tdescending_traj.ends_at_lmax: %d\n",
                (int)descending_traj.ends_at_critical, (int)descending_traj.ends_at_lmax);
#endif

        // this needs to be done with care
        // what if ascending_traj turns out to be a descending one - finds itself at a lmin?
        // what if descending_traj turns out to be a ascending one - finds itself at a lmax?
        gradient_flow_t::trajectory_t trajectory;
        trajectory.vertices.insert(
            trajectory.vertices.end(),
            descending_traj.path.rbegin(),
            descending_traj.path.rend() - 1
        );
        trajectory.vertices.insert(
            trajectory.vertices.end(),
            ascending_traj.path.begin(),
            ascending_traj.path.end()
        );

        // Determine trajectory type
        if ((descending_traj.ends_at_critical && !descending_traj.ends_at_lmax) &&
            (ascending_traj.ends_at_critical  && ascending_traj.ends_at_lmax)) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_LMAX;
        } else if (
            (descending_traj.ends_at_critical && !descending_traj.ends_at_lmax && !ascending_traj.ends_at_critical) ||
            (ascending_traj.ends_at_critical  && !ascending_traj.ends_at_lmax  && !descending_traj.ends_at_critical)
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_ONLY;
        } else if (
            (descending_traj.ends_at_critical && descending_traj.ends_at_lmax && !ascending_traj.ends_at_critical) ||
            (ascending_traj.ends_at_critical  && ascending_traj.ends_at_lmax  && !descending_traj.ends_at_critical)
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMAX_ONLY;
        } else if (
            descending_traj.ends_at_critical && !descending_traj.ends_at_lmax &&
            ascending_traj.ends_at_critical  && !ascending_traj.ends_at_lmax
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMIN_LMIN;
        } else if (
            descending_traj.ends_at_critical && descending_traj.ends_at_lmax &&
            ascending_traj.ends_at_critical  && ascending_traj.ends_at_lmax
            ) {
            trajectory.trajectory_type = gradient_flow_t::LMAX_LMAX;
        } else {
            // This shouldn't happen if we're properly finding critical points
            REPORT_WARNING("Found trajectory with no critical endpoints");
            trajectory.trajectory_type = gradient_flow_t::UNKNOWN;
        }

#if DEBUG__compute_gradient_flow
        print_vect(trajectory.vertices,"trajectory.vertices");

        switch (trajectory.trajectory_type) {
        case gradient_flow_t::LMIN_LMAX:
            Rprintf("trajectory.trajectory_type: LMIN_LMAX\n");
            break;

        case gradient_flow_t::LMIN_ONLY:
            Rprintf("trajectory.trajectory_type: LMIN_ONLY\n");
            break;

        case gradient_flow_t::LMAX_ONLY:
            Rprintf("trajectory.trajectory_type: LMAX_ONLY\n");
            break;

        case gradient_flow_t::UNKNOWN:
            Rprintf("trajectory.trajectory_type: UNKNOWN\n");
            break;

        default:
            Rprintf("trajectory.trajectory_type: ???\n");
            break;
        }

        if (while_counter == 2)
            REPORT_ERROR("DEBUGGING compute_gradient_flow\n");
        while_counter++;
#endif

        result.trajectories.push_back(trajectory);

        // Remove processed vertices
        for (size_t vertex : trajectory.vertices) {
            unprocessed_vertices.erase(vertex);
        }
    }


    // Computing ascending and descending basins
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];

        // Handle each trajectory based on its type
        switch (traj.trajectory_type) {
        case gradient_flow_t::LMIN_LMAX:
            // Add all vertices to both basins using range insertion
            result.ascending_basin_map[traj.vertices.front()].insert(
                traj.vertices.begin(), traj.vertices.end());
            result.descending_basin_map[traj.vertices.back()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::LMIN_ONLY:
            // Only min is a critical point
            result.ascending_basin_map[traj.vertices.front()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::LMAX_ONLY:
            // Only max is a critical point
            result.descending_basin_map[traj.vertices.back()].insert(
                traj.vertices.begin(), traj.vertices.end());
            break;

        case gradient_flow_t::UNKNOWN:
            REPORT_WARNING("Found trajectory of UNKNOWN type");
            break;

        default:
            REPORT_WARNING("Found trajectory with no type assignment");
            break;
        }
    }

    // Create map of (min, max) pairs to the set = the union of the correspoing trajectory vertices
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_LMAX) {
            result.procell_map[{traj.vertices.front(), traj.vertices.back()}].insert(
                traj.vertices.begin(), traj.vertices.end());
        }
    }
    // After all proper pro-cells are populated we add to them corresponding dangling trajectories
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_ONLY) {
            // How to insert traj.vertices into all pro-cells that have traj.vertices.front() as the local minimum?
            // loop over all pro-cell with traj.vertices.front() as the lmin
            // procell_map[{traj.vertices.front(), ???}].insert(
            //     traj.vertices.begin(), traj.vertices.end());
        }
    }

    // Create map of (min, max) pairs to the set = the union of the correspoing trajectory vertices
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];
        if (traj.trajectory_type == gradient_flow_t::LMIN_LMAX) {
            result.procell_map[{traj.vertices.front(), traj.vertices.back()}].insert(
                traj.vertices.begin(), traj.vertices.end());
        }
    }

    // After all proper pro-cells are populated we add to them corresponding dangling trajectories
    for (size_t i = 0; i < result.trajectories.size(); ++i) {
        const auto& traj = result.trajectories[i];

        if (traj.trajectory_type == gradient_flow_t::LMIN_ONLY) {
            // Add vertices to all pro-cells that have the same local minimum
            size_t lmin = traj.vertices.front();

            for (auto& [key, vertices] : result.procell_map) {
                if (key.first == lmin) {
                    // This pro-cell has the same local minimum
                    vertices.insert(traj.vertices.begin(), traj.vertices.end());
                }
            }
        }
        else if (traj.trajectory_type == gradient_flow_t::LMAX_ONLY) {
            // Add vertices to all pro-cells that have the same local maximum
            size_t lmax = traj.vertices.back();

            for (auto& [key, vertices] : result.procell_map) {
                if (key.second == lmax) {
                    // This pro-cell has the same local maximum
                    vertices.insert(traj.vertices.begin(), traj.vertices.end());
                }
            }
        }
    }

    return result;
}
