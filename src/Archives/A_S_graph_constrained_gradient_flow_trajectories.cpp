extern "C" {
    SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    // SEXP S_graph_gradient_flow_trajectories(SEXP s_graph, SEXP s_Ey);
    // SEXP S_graph_hop_gradient_flow_trajectories(SEXP s_graph, SEXP s_Ey, SEXP s_hop_size);
    // SEXP S_graph_hop_gradient_flow_trajectories_optimized(SEXP s_graph, SEXP s_Ey, SEXP s_hop_size);
    SEXP S_graph_MS_cx(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey, SEXP s_hop_size, SEXP s_method);
    //SEXP S_graph_MS_cx_min_cell(SEXP s_graph, SEXP s_Ey, SEXP s_hop_radius, SEXP s_method, SEXP s_min_cell_size);
}

/**
 * @brief Computes constrained gradient flow trajectories on a graph while preventing basin-jumping
 *
 * @details This function computes gradient flow trajectories of a function Ey defined on vertices
 * of a graph G = pIG_k^h(X), which is the h-th power of a pruned intersection kNN graph.
 * For each vertex v, it constructs:
 * 1. An ascending trajectory that follows the steepest ascent
 * 2. A descending trajectory that follows the steepest descent
 *
 * To prevent basin-jumping between different local maxima regions, the function:
 * - Uses a core graph (pIG_k(X)) to find valid paths between vertices
 * - Enforces monotonicity along paths in the core graph
 * - Prioritizes shorter paths when selecting next vertices
 *
 * The trajectory construction process:
 * 1. For each vertex v in the graph:
 *    a. Finds all neighbors w in G connected to v
 *    b. For each neighbor w:
 *       - Computes shortest path γ_vw in core graph
 *       - Verifies Ey is strictly increasing/decreasing along γ_vw
 *       - Sorts valid neighbors by path length (descending) and gradient (descending)
 *    c. Selects neighbor with valid path and maximum gradient
 *    d. Repeats until no valid neighbors exist
 *
 * @param s_graph SEXP containing h-th power graph pIG_k^h(X) as adjacency lists
 * @param s_core_graph SEXP containing core graph pIG_k(X) as adjacency lists
 * @param s_Ey SEXP containing function values at vertices
 *
 * @return SEXP List containing:
 *         - lext: Matrix (n_vertices × 2) with local minima and maxima for each vertex
 *         - trajectories: List of integer vectors, each containing complete trajectory
 *           (descending path reversed + ascending path) for corresponding vertex
 *         - local_maxima: Integer vector containing indices of all local maxima
 *         - local_minima: Integer vector containing indices of all local minima
 *
 * @throws R error if input lengths don't match
 */

// Helper struct to store neighbor information
struct neighbor_info_t {
    int vertex;
    std::vector<int> path;
    int path_length;
    double diff_value;
};

// Helper function to check if path is monotonic in Ey
bool is_path_monotonic_increasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] <= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_path_monotonic_decreasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] >= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_local_maximum(int vertex, const std::vector<std::vector<int>>& graph, const double* Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] < Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

bool is_local_minimum(int vertex, const std::vector<std::vector<int>>& graph, const double* Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] > Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {
    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);

    int nprot = 0;
    SEXP local_extrema = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    int* local_extrema_ptr = INTEGER(local_extrema);

    // Sets to store unique local maxima and minima
    std::set<int> local_maxima_set;
    std::set<int> local_minima_set;

    // Helper function to get valid neighbors with their paths
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : graph[current_vertex]) {
            std::vector<int> path = shortest_path_by_hops(core_graph, current_vertex, neighbor);
            if (path.empty()) continue;

            double diff = ascending ?
                         Ey[neighbor] - Ey[current_vertex] :
                         Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                is_path_monotonic_increasing(path, Ey) :
                is_path_monotonic_decreasing(path, Ey))) {
                neighbors.emplace_back(neighbor_info_t{neighbor, path, static_cast<int>(path.size() - 1), diff});
            }
        }

        // Sort by difference value in descending order
        std::sort(neighbors.begin(), neighbors.end(),
                 [](const neighbor_info_t& a, const neighbor_info_t& b) {
                     return a.diff_value > b.diff_value;
                 });

        return neighbors;
    };

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        bool is_max = is_local_maximum(vertex, graph, Ey);

        // Ascending trajectory
        std::vector<int> ascending_trajectory = {vertex};
        int current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, true);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                ascending_trajectory.insert(ascending_trajectory.end(),
                                         path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }

            // Check if endpoint is a local maximum
            if (is_local_maximum(current_vertex, graph, Ey)) {
                local_maxima_set.insert(current_vertex);
            }
        } else {
            local_maxima_set.insert(vertex);
        }

        // Descending trajectory
        std::vector<int> descending_trajectory = {vertex};
        current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, false);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                descending_trajectory.insert(descending_trajectory.end(),
                                          path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }

            // Check if endpoint is a local minimum
            if (is_local_minimum(current_vertex, graph, Ey)) {
                local_minima_set.insert(current_vertex);
            }
        }

        if (is_max) {
            // Store local extrema
            local_extrema_ptr[vertex] = -1;
            local_extrema_ptr[vertex + n_vertices] = vertex;

            // Single vertex trajectory for local maximum
            SEXP trajectory = PROTECT(allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
            // Store local extrema
            local_extrema_ptr[vertex] = descending_trajectory.back();
            local_extrema_ptr[vertex + n_vertices] = ascending_trajectory.back();

            // Construct full trajectory
            descending_trajectory.erase(descending_trajectory.begin());
            int tr_size = ascending_trajectory.size() + descending_trajectory.size();
            SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
            int* trajectory_ptr = INTEGER(trajectory);

            std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
            std::copy(ascending_trajectory.begin(), ascending_trajectory.end(),
                     trajectory_ptr + descending_trajectory.size());

            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        }
    }

    // Convert sets to R vectors
    SEXP local_maxima = PROTECT(allocVector(INTSXP, local_maxima_set.size())); nprot++;
    SEXP local_minima = PROTECT(allocVector(INTSXP, local_minima_set.size())); nprot++;

    int i = 0;
    for (int max_idx : local_maxima_set) {
        INTEGER(local_maxima)[i++] = max_idx;
    }

    i = 0;
    for (int min_idx : local_minima_set) {
        INTEGER(local_minima)[i++] = min_idx;
    }

    // Construct result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(result_list, 0, local_extrema);
    SET_VECTOR_ELT(result_list, 1, trajectories);
    SET_VECTOR_ELT(result_list, 2, local_maxima);
    SET_VECTOR_ELT(result_list, 3, local_minima);

    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    SET_STRING_ELT(names, 2, mkChar("local_maxima"));
    SET_STRING_ELT(names, 3, mkChar("local_minima"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return result_list;
}


/// previous version

/**
 * @brief Computes constrained gradient flow trajectories on a graph while preventing basin-jumping
 *
 * @details This function computes gradient flow trajectories of a function Ey defined on vertices
 * of a graph G = pIG_k^h(X), which is the h-th power of a pruned intersection kNN graph.
 * For each vertex v, it constructs:
 * 1. An ascending trajectory that follows the steepest ascent
 * 2. A descending trajectory that follows the steepest descent
 *
 * To prevent basin-jumping between different local maxima regions, the function:
 * - Uses a core graph (pIG_k(X)) to find valid paths between vertices
 * - Enforces monotonicity along paths in the core graph
 * - Prioritizes shorter paths when selecting next vertices
 *
 * The trajectory construction process:
 * 1. For each vertex v in the graph:
 *    a. Finds all neighbors w in G connected to v
 *    b. For each neighbor w:
 *       - Computes shortest path γ_vw in core graph
 *       - Verifies Ey is strictly increasing/decreasing along γ_vw
 *       - Sorts valid neighbors by path length (descending) and gradient (descending)
 *    c. Selects neighbor with valid path and maximum gradient
 *    d. Repeats until no valid neighbors exist
 *
 * @param s_graph SEXP containing h-th power graph pIG_k^h(X) as adjacency lists
 * @param s_core_graph SEXP containing core graph pIG_k(X) as adjacency lists
 * @param s_Ey SEXP containing function values at vertices
 *
 * @return SEXP List containing:
 *         - lext: Matrix (n_vertices × 2) with local minima and maxima for each vertex
 *         - trajectories: List of integer vectors, each containing complete trajectory
 *           (descending path reversed + ascending path) for corresponding vertex
 *
 * @throws R error if input lengths don't match
 *
 * @note The function ensures that trajectories cannot jump between basins of attraction
 *       of different local extrema by enforcing strict monotonicity along paths in
 *       the core graph.
 */

// Helper struct to store neighbor information
struct neighbor_info_t {
    int vertex;
    std::vector<int> path;
    int path_length;
    double diff_value;
};

// Helper function to check if path is monotonic in Ey
bool is_path_monotonic_increasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] <= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_path_monotonic_decreasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] >= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_local_maximum(int vertex, const std::vector<std::vector<int>>& graph, const double* Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] < Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {

    #define DEBUG__graph_constrained_gradient_flow_trajectories 0

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

#if DEBUG__graph_constrained_gradient_flow_trajectories
    Rprintf("In S_graph_constrained_gradient_flow_trajectories\n");
    Rprintf("n_vertices: %d\n", n_vertices);
#endif

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);

#if DEBUG__graph_constrained_gradient_flow_trajectories
    print_vect_vect(graph, "graph");
    print_vect_vect(core_graph, "core_graph");
#endif

    int nprot = 0;
    SEXP local_extrema = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    int* local_extrema_ptr = INTEGER(local_extrema);

    // Helper function to get valid neighbors with their paths
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : graph[current_vertex]) {
            std::vector<int> path = shortest_path_by_hops(core_graph, current_vertex, neighbor);
            if (path.empty()) continue;

            double diff = ascending ?
                         Ey[neighbor] - Ey[current_vertex] :
                         Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                is_path_monotonic_increasing(path, Ey) :
                is_path_monotonic_decreasing(path, Ey))) {
                neighbors.emplace_back(neighbor_info_t{neighbor, path, static_cast<int>(path.size() - 1), diff});
            }
        }

        // Sort by difference value in the descending order
        std::sort(neighbors.begin(), neighbors.end(),
         [](const neighbor_info_t& a, const neighbor_info_t& b) {
             return a.diff_value > b.diff_value;
         });

        return neighbors;
    };

    for (int vertex = 0; vertex < n_vertices; vertex++) {

        bool is_max = is_local_maximum(vertex, graph, Ey);

#if DEBUG__graph_constrained_gradient_flow_trajectories
        Rprintf("\n-------\nvertex: %d\t is_max: %s\n", vertex, is_max ? "true" : "false");
#endif

        // Ascending trajectory
        std::vector<int> ascending_trajectory = {vertex};
        int current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, true);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                // Skip first vertex and append rest of path
                // for (size_t i = 1; i < path.size(); ++i)
                //     ascending_trajectory.push_back(path[i]);
                ascending_trajectory.insert(ascending_trajectory.end(),
                                            path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }
        }

        // Descending trajectory
        std::vector<int> descending_trajectory = {vertex};
        current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, false);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                // Skip first vertex and append rest of path
                descending_trajectory.insert(descending_trajectory.end(),
                                            path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }
        }

#if DEBUG__graph_constrained_gradient_flow_trajectories
        print_vect(ascending_trajectory, "ascending_trajectory");
        print_vect(descending_trajectory, "descending_trajectory");
#endif

        if (is_max) {
            // Store local extrema
            local_extrema_ptr[vertex] = -1; // local maximum has not local minimum associated with it
            local_extrema_ptr[vertex + n_vertices] = vertex;

            // Single vertex trajectory for local maximum
            SEXP trajectory = PROTECT(allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
            // Store local extrema
            local_extrema_ptr[vertex] = descending_trajectory.back();
            local_extrema_ptr[vertex + n_vertices] = ascending_trajectory.back();

            // Construct full trajectory
            descending_trajectory.erase(descending_trajectory.begin());
            int tr_size = ascending_trajectory.size() + descending_trajectory.size();
            SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
            int* trajectory_ptr = INTEGER(trajectory);

            std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
            std::copy(ascending_trajectory.begin(), ascending_trajectory.end(),
                      trajectory_ptr + descending_trajectory.size());

            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        }
    }

    // Construct result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(result_list, 0, local_extrema);
    SET_VECTOR_ELT(result_list, 1, trajectories);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return result_list;
}
