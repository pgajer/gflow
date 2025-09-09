// ##
// ##
// ##

// Nov 4, 2024

#if 0
//
// The uniquness of solutions of ODE's implies that gradient field trajectories on
// a smooth manifold are either identical or disjoint. This does not have to be the
// case in a discreate case of gradient trajectories on graphs. For example,
// consider a three arm star graph with a function with a local miniumum at the end
// vertex of one arm and two local maxima at the ends of the other two arms, then
// the gradient trajectories of this function of points on the two local maxima
// arms will merge over the local minium arm. Thus, they are not disjoint.

// In the following graph_MS_cx function gradient trajectories of Ey are computed
// and stored in the trajectory_t struct. What do you think about streamlining

// 1) the return struct of this funcition and
// 2) the trajectory algorithm


// by

// 1) Merging ascending and descending trajectory of the given vertext into one trajectory

//         // Remove the first element i from descending_trajectory as it is
//         // already in ascending_trajectory
//         descending_trajectory.erase(descending_trajectory.begin());

//         // Create an R integer vector for the trajectory
//         int tr_size = ascending_trajectory.size() + descending_trajectory.size();
//         std::vector<int> trajectory(tr_size);

//         // Copy elements using std::copy instead of manually copying them using for loop
//         std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory);
//         std::copy(ascending_trajectory.begin(), ascending_trajectory.end(), trajectory + descending_trajectory.size());


// 2) Storing unique trajectories in a std::set<std::vector<int>> or any other structure that autmatically removes duplicate trajectories

//    std::set<std::vector<int>> trajectories_set;

// This will substantially decrease memory usage as we store on unique trajectories

// 3) Keeping track for each vertex  which unique trajectories contain that vertex

// One possible data structure for this could be std::map<int, std::set<int>>
// vertex_trajectories, where the set stores the indices of the trajectories. The
// issue here is that as unique trajectories are store (at least) initially in a
// std::set, we don't have indices of the trajectories in the set. How would you
// resolve this issue in order to simultenously keep track of unique trajectories
// and map each vertex to a set of trajectories which contain that vertex.

// Note that once a trajectory of the given vertex is found. For all vertices of
// that trajectory we can update vertex_trajectories map adding the given trajector
// to each vertex.

// 4) I am not sure if these structures can lead to any time gains in the
// trajectory detection algorithm. I am not sure as the issue is that a given
// vertex may beloong to multiple gradient trajectories.

// Please propose a modified version of graph_MS_cx that incorporates the above suggestions and also puts all the objects that it returns into one MS_complex_t struct

/**
 * @brief A structure representing a Morse-Smale complex
 *
 * The Morse-Smale complex captures the topological relationships between critical points
 * (local maxima and minima) of a scalar function defined on a graph, along with the
 * connecting trajectories and cells between these critical points.
 */
struct MS_complex_t {
    std::map<int, std::set<int>> lmax_to_lmin;    ///< Maps local maxima to their connected local minima
    std::map<int, std::set<int>> lmin_to_lmax;    ///< Maps local minima to their connected local maxima
    std::set<int> local_maxima;                   ///< Set of all local maxima vertices
    std::set<int> local_minima;                   ///< Set of all local minima vertices
    std::map<std::pair<int,int>, std::set<int>> procells;    ///< Maps (max,min) pairs to their proto-cells
    std::map<std::pair<int,int>, std::vector<std::set<int>>> cells;  ///< Maps (max,min) pairs to their decomposed cells
};

/**
 * @brief A structure representing ascending and descending trajectories from a vertex
 *
 * Stores the complete paths from a vertex to its associated local maximum (ascending)
 * and local minimum (descending).
 */
struct trajectory_t {
    std::vector<int> ascending;     ///< Path from vertex to local maximum
    std::vector<int> descending;    ///< Path from vertex to local minimum
};


std::pair<MS_complex_t, std::vector<trajectory_t>>
graph_MS_cx(const std::vector<std::vector<int>>& adj_list,
            const std::vector<std::vector<int>>& core_adj_list,
            const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
            const std::vector<double>& Ey) {

    MS_complex_t ms_cx;
    std::vector<trajectory_t> trajectories(adj_list.size());

    // Helper function to get valid neighbors
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : adj_list[current_vertex]) {
            auto path_key = std::make_pair(current_vertex, neighbor);
            auto path_it = shortest_paths.find(path_key);
            if (path_it == shortest_paths.end()) continue;

            const auto& path = path_it->second;
            if (path.empty()) continue;

            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                is_path_monotonic_increasing(path, Ey.data()) :
                is_path_monotonic_decreasing(path, Ey.data()))) {
                neighbors.push_back({neighbor, path, static_cast<int>(path.size() - 1), diff});
            }
        }

        // Sort by difference value in descending order
        std::sort(neighbors.begin(), neighbors.end(),
                 [](const neighbor_info_t& a, const neighbor_info_t& b) {
                     return a.diff_value > b.diff_value;
                 });

        return neighbors;
    };

    // Process each vertex
    for (size_t vertex = 0; vertex < adj_list.size(); vertex++) {
        bool is_max = is_local_maximum(vertex, adj_list, Ey.data());
        bool is_min = is_local_minimum(vertex, adj_list, Ey.data());

        if (is_max) {
            ms_cx.local_maxima.insert(vertex);
        }
        if (is_min) {
            ms_cx.local_minima.insert(vertex);
        }

        // Compute ascending and descending trajectories
        std::vector<int> ascending = {static_cast<int>(vertex)};
        std::vector<int> descending = {static_cast<int>(vertex)};
        int asc_vertex = vertex;
        int desc_vertex = vertex;

        // Compute ascending trajectory if not at local maximum
        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(asc_vertex, true);
                if (valid_neighbors.empty()) break;
                auto& path = valid_neighbors[0].path;
                ascending.insert(ascending.end(), path.begin() + 1, path.end());
                asc_vertex = valid_neighbors[0].vertex;
            }

            if (is_local_maximum(asc_vertex, adj_list, Ey.data())) {
                ms_cx.local_maxima.insert(asc_vertex);
            }
        }

        // Compute descending trajectory if not at local minimum
        if (!is_min) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(desc_vertex, false);
                if (valid_neighbors.empty()) break;
                auto& path = valid_neighbors[0].path;
                descending.insert(descending.end(), path.begin() + 1, path.end());
                desc_vertex = valid_neighbors[0].vertex;
            }

            if (is_local_minimum(desc_vertex, adj_list, Ey.data())) {
                ms_cx.local_minima.insert(desc_vertex);
            }
        }

        // Store trajectories
        trajectories[vertex].ascending = ascending;
        trajectories[vertex].descending = descending;

        // Update connectivity maps and pro-cells if trajectory connects local extrema
        if (is_local_minimum(desc_vertex, adj_list, Ey.data()) &&
            is_local_maximum(asc_vertex, adj_list, Ey.data())) {

            // Update connectivity maps
            ms_cx.lmax_to_lmin[asc_vertex].insert(desc_vertex);
            ms_cx.lmin_to_lmax[desc_vertex].insert(asc_vertex);

            // Update pro-cell
            std::pair<int,int> cell_key(asc_vertex, desc_vertex);

            // Add all vertices from the trajectory to the pro-cell
            ms_cx.procells[cell_key].insert(descending.begin(), descending.end());
            ms_cx.procells[cell_key].insert(ascending.begin(), ascending.end());
        }
    }

    // Compute MS cells from pro-cells
    for (const auto& [key, procell] : ms_cx.procells) {
        if (procell.size() > 2) {
            // Remove extrema from the set for component calculation
            std::set<int> cell_vertices = procell;
            cell_vertices.erase(key.first);   // Remove local maximum
            cell_vertices.erase(key.second);  // Remove local minimum

            // Find connected components
            auto components = count_subgraph_set_components(core_adj_list, cell_vertices);

            // Store components as MS cells
            std::vector<std::set<int>> cell_components;
            for (const auto& [comp_id, comp_vertices] : components) {
                std::set<int> new_comp = comp_vertices;
                new_comp.insert(key.first);
                new_comp.insert(key.second);
                cell_components.push_back(new_comp);
            }
            ms_cx.cells[key] = cell_components;
        } else {
            std::vector<std::set<int>> cell_components;
            cell_components.push_back(procell);
            ms_cx.cells[key] = cell_components;
        }
    }

    return {ms_cx, trajectories};
}

#endif


// Oct 31, 2024

// - This struct provides a hash function for `std::pair<T1, T2>` objects.
// - It combines the hash values of the first and second elements of the pair using the XOR (^) operator.
// - This hash function is used in the `std::unordered_map` to efficiently store and retrieve key-value pairs.
struct pair_hash {
    template<typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

/*

I've reviewed the function implementation and documentation. Overall, the implementation seems to be correct and follows the described functionality. However, there are a few areas where potential improvements could be made:

1. Correctness:
   The implementation appears to correctly calculate the gradient flow trajectories and local extrema as described in the documentation.

2. Areas for potential improvement:

   a. Index conversion:
   The documentation states that vertex indices in the input should be 1-based, but the implementation treats them as 0-based. This discrepancy should be addressed either by adjusting the implementation or updating the documentation.

   b. Error handling:
   The function checks if the lengths of s_graph and s_Ey are the same, but it doesn't verify if the graph is empty or if s_Ey contains valid numeric values. Adding these checks could improve robustness.

   c. Memory management:
   The function uses PROTECT/UNPROTECT for memory management, which is correct. However, consider using PROTECT_WITH_INDEX and REPROTECT for more flexible protection, especially when dealing with loops that allocate memory.

   d. Optimization:
   The function performs two separate traversals (ascending and descending) for each vertex. It might be possible to combine these into a single traversal to improve efficiency.

   e. Code duplication:
   The logic for finding ascending and descending trajectories is very similar. Consider refactoring this into a separate function to reduce code duplication.

   f. Use of C++11 features:
   Since the code uses C++11 features (like std::move), consider using more modern C++ features such as auto, range-based for loops, and std::vector::emplace_back() for potential performance improvements.

   g. Naming convention:
   The function name and some variable names use different naming conventions (e.g., S_graph_gradient_flow_trajectories vs s_Ey). Consider standardizing the naming convention for better readability.

   h. Comments:
   Adding more inline comments explaining the logic, especially for the trajectory calculations, could improve code maintainability.

   i. Return value adjustment:
   The function returns 1-based indices for the local extrema (as they are stored directly in an R matrix), but 0-based indices for the trajectories. This inconsistency might be confusing for users. Consider adjusting one or the other for consistency.

Overall, the implementation is functional and correct, but these improvements could enhance its robustness, efficiency, and maintainability.

*/

/**
 * Estimates graph gradient flow trajectories of a response variable for each vertex of a graph.
 *
 * This function estimates the graph gradient flow trajectories of a response
 * variable for each vertex of the graph. It takes an R list representing the
 * graph adjacency list (`s_graph`) and an R numeric vector containing the
 * estimates of a conditional expectation of a random variable over the graph
 * vertices (`s_Ey`). It returns an R list with two components: "lext" (an
 * integer matrix containing local minimum and maximum indices) and
 * "trajectories" (a list of integer vectors representing the gradient flow
 * trajectories). The function performs the main computation by finding the
 * ascending and descending trajectories for each vertex based on the gradient
 * of the response variable.
 *
 * @param s_graph An R list representing the graph adjacency list. Each element of
 *               the list should be an integer vector specifying the neighboring
 *               vertex indices for a given vertex. The vertex indices should be 1-based.
 *
 * @param s_Ey    An R numeric vector containing estimates of a conditional expectation
 *               of a random variable over the graph vertices.
 *
 * @return An R list with two components:
 *         1) "lext": An integer matrix with dimensions (n_vertices, 2), where the first
 *                    column contains the local minimum indices and the second column
 *                    contains the local maximum indices for each vertex.
 *         2) "trajectories": A list of integer vectors representing the gradient flow
 *                           trajectories for each vertex of the graph. Each trajectory
 *                           includes both the descending and ascending paths from the
 *                           vertex to its local minimum and maximum, respectively.
 */
SEXP S_graph_gradient_flow_trajectories(SEXP s_graph, SEXP s_Ey) {

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices) {
        error("The lengths of s_graph and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    double* Ey = REAL(s_Ey);

    int nprot = 0;
    // Create list elements
    SEXP local_extrema = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    int* local_extrema_ptr = INTEGER(local_extrema);

    // Main computation
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        // Find ascending trajectory
        std::vector<int> ascending_trajectory = {vertex};
        int current_vertex = vertex;
        double max_diff = 0.0;

        while (true) {
            int next_vertex = -1;
            max_diff = 0.0;
            for (int neighbor : graph[current_vertex]) {
                double diff = Ey[neighbor] - Ey[current_vertex];
                if (diff > max_diff) {
                    max_diff = diff;
                    next_vertex = neighbor;
                }
            }

            if (next_vertex == -1) break;
            ascending_trajectory.push_back(next_vertex);
            current_vertex = next_vertex;
        }

        // Find descending trajectory (similar logic)
        std::vector<int> descending_trajectory = {vertex};
        current_vertex = vertex;


        while (true) {
            int next_vertex = -1;

            max_diff = 0.0;
            for (int neighbor : graph[current_vertex]) {
                double diff = Ey[current_vertex] - Ey[neighbor];

                if (diff > max_diff) {
                    max_diff = diff;
                    next_vertex = neighbor;
                }
            }

            if (next_vertex == -1) break;
            descending_trajectory.push_back(next_vertex);
            current_vertex = next_vertex;
        }

        // Store local extrema in the matrix
        local_extrema_ptr[vertex] = descending_trajectory.back();
        local_extrema_ptr[vertex + n_vertices] = ascending_trajectory.back();

        // Trajectory
        // Remove the first element i from descending_trajectory as it is
        // already in ascending_trajectory
        descending_trajectory.erase(descending_trajectory.begin());

        // Create an R integer vector for the trajectory
        int tr_size = ascending_trajectory.size() + descending_trajectory.size();
        SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
        int* trajectory_ptr = INTEGER(trajectory);

        // Copy elements using std::copy instead of manually copying them using for loop
        std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
        std::copy(ascending_trajectory.begin(), ascending_trajectory.end(), trajectory_ptr + descending_trajectory.size());

        SET_VECTOR_ELT(trajectories, vertex, trajectory);
        UNPROTECT(1);
    }

    // Allocate memory for the result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 2)); nprot++;
    // Set elements of the result list
    SET_VECTOR_ELT(result_list, 0, local_extrema);
    SET_VECTOR_ELT(result_list, 1, trajectories);

    // Add names
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result_list;
}

/**
 * @brief Estimates graph gradient flow trajectories using a hop-based neighborhood approach.
 *
 * This function estimates the graph gradient flow trajectories of a response
 * variable for each vertex of the graph using an extended neighborhood defined
 * by a hop radius. It provides more robust trajectory estimation by considering
 * vertices within a specified number of hops, potentially avoiding small-scale
 * local extrema.
 *
 * @param s_graph An R list representing the graph adjacency list. Each element of
 *               the list should be an integer vector specifying the neighboring
 *               vertex indices for a given vertex. The vertex indices should be 1-based.
 * @param s_Ey    An R numeric vector containing estimates of a conditional expectation
 *               of a random variable over the graph vertices.
 * @param s_hop_radius An R integer specifying the number of hops to consider for the
 *                  neighborhood of each vertex. Default is 1, which is equivalent
 *                  to the original algorithm considering only immediate neighbors.
 *
 * @return An R list with two components:
 *         1) "lext": An integer matrix with dimensions (n_vertices, 2), where the first
 *                    column contains the local minimum indices and the second column
 *                    contains the local maximum indices for each vertex.
 *         2) "trajectories": A list of integer vectors representing the gradient flow
 *                            trajectories for each vertex of the graph. Each trajectory
 *                            includes both the descending and ascending paths from the
 *                            vertex to its local minimum and maximum, respectively.
 *
 * @details The function uses a hop-based approach to define the neighborhood of each
 *          vertex. For each vertex, it considers all vertices that can be reached
 *          within 'hop_radius' steps. This extended neighborhood is used to determine
 *          the direction of steepest ascent/descent for constructing the gradient
 *          flow trajectories.
 *
 *          The hop-based approach allows for more robust trajectory estimation,
 *          particularly in graphs with noise or small-scale fluctuations. By
 *          considering a larger neighborhood, the algorithm can potentially bypass
 *          small local extrema and capture more significant trends in the data.
 *
 * @note The time complexity of this algorithm increases with larger hop sizes,
 *       as it needs to explore a larger neighborhood for each vertex.
 *
 * @see S_graph_gradient_flow_trajectories for the original version of this function
 *      that considers only immediate neighbors.
 *
 * @warning For large graphs or large hop sizes, this function may be computationally
 *          expensive. Consider the trade-off between robustness and computation time
 *          when choosing the hop size.
 */
SEXP S_graph_hop_gradient_flow_trajectories(SEXP s_graph, SEXP s_Ey, SEXP s_hop_radius) {

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices) {
        error("The lengths of s_graph and s_Ey must be the same");
    }

    int hop_radius = INTEGER(s_hop_radius)[0];
    if (hop_radius < 0) {
        error("hop_radius must be non-negative");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    double* Ey = REAL(s_Ey);

    int nprot = 0;
    // Create list elements
    SEXP local_extrema = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    int* local_extrema_ptr = INTEGER(local_extrema);

    /**
     * @brief Lambda function to find all vertices within a specified hop distance from a start vertex.
     *
     * This lambda function performs a breadth-first search (BFS) to identify all vertices
     * that are within a specified number of hops from the start vertex in the given graph.
     *
     * @param[in] graph      A const reference to a vector of vector of integers representing
     *                       the adjacency list of the graph. graph[i] contains the indices
     *                       of vertices adjacent to vertex i.
     * @param[in] start_vertex The index of the starting vertex for the neighborhood search.
     * @param[in] hop_radius   The maximum number of hops to consider in the neighborhood.
     *
     * @return A unique pointer to a vector of integers containing the indices of all
     *         vertices within the specified hop distance from the start vertex.
     *         The start vertex itself is not included in this vector.
     *
     * @note The function uses a breadth-first search approach, which ensures that
     *       vertices are discovered in order of their distance from the start vertex.
     *
     * @warning The function assumes that the graph is well-formed and that all vertex
     *          indices in the graph are valid. It does not perform bounds checking on
     *          the vertex indices.
     *
     * @see S_graph_hop_gradient_flow_trajectories, which uses this function to compute
     *      extended neighborhoods for robust gradient flow trajectory estimation.
     *
     * Example usage:
     * @code
     * auto get_neighborhood = [](const std::vector<std::vector<int>>& graph, int start_vertex, int hop_radius) -> std::unique_ptr<std::vector<int>> {
     *     // Lambda function implementation...
     * };
     *
     * // Using the lambda function
     * auto neighborhood_ptr = get_neighborhood(graph, current_vertex, hop_radius);
     * for (int neighbor : *neighborhood_ptr) {
     *     // Process each neighbor
     * }
     * @endcode
     */
    auto get_neighborhood = [](const std::vector<std::vector<int>>& graph, int start_vertex, int hop_radius) -> std::unique_ptr<std::vector<int>> {

        auto neighborhood = std::make_unique<std::vector<int>>();
        std::vector<bool> visited(graph.size(), false);
        std::queue<std::pair<int, int>> q;  // vertex, distance

        q.push({start_vertex, 0});
        visited[start_vertex] = true;

        while (!q.empty()) {
            int vertex = q.front().first;
            int distance = q.front().second;
            q.pop();

            if (distance > 0) {  // Don't include the start vertex
                neighborhood->push_back(vertex);
            }

            if (distance < hop_radius) {
                for (int neighbor : graph[vertex]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        q.push({neighbor, distance + 1});
                    }
                }
            }
        }

        return neighborhood;
    };

    std::vector<std::vector<int>> neighborhoods(n_vertices);
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        auto neighborhood_ptr = get_neighborhood(graph, vertex, hop_radius);
        neighborhoods[vertex] = *neighborhood_ptr;
    }

    // Main computation
    for (int vertex = 0; vertex < n_vertices; vertex++) {

        // Finding the ascending trajectory of 'vertex'
        std::vector<int> ascending_trajectory = {vertex};
        int current_vertex = vertex;
        double max_diff = 0.0;

        while (true) {
            int next_vertex = -1;
            max_diff = 0.0;
            for (int neighbor : neighborhoods[current_vertex]) {
                double diff = Ey[neighbor] - Ey[current_vertex];
                if (diff > max_diff) {
                    max_diff = diff;
                    next_vertex = neighbor;
                }
            }

            if (next_vertex == -1) break;
            ascending_trajectory.push_back(next_vertex);
            current_vertex = next_vertex;
        }

        // Finding the descending trajectory of 'vertex' (similar logic)
        std::vector<int> descending_trajectory = {vertex};
        current_vertex = vertex;

        while (true) {
            int next_vertex = -1;
            max_diff = 0.0;
            for (int neighbor : neighborhoods[current_vertex]) {
                double diff = Ey[current_vertex] - Ey[neighbor];
                if (diff > max_diff) {
                    max_diff = diff;
                    next_vertex = neighbor;
                }
            }

            if (next_vertex == -1) break;
            descending_trajectory.push_back(next_vertex);
            current_vertex = next_vertex;
        }

        // Store local extrema in the matrix
        local_extrema_ptr[vertex] = descending_trajectory.back();
        local_extrema_ptr[vertex + n_vertices] = ascending_trajectory.back();

        // Trajectory
        // Remove the first element i from descending_trajectory as it is
        // already in ascending_trajectory
        descending_trajectory.erase(descending_trajectory.begin());

        // Create an R integer vector for the trajectory
        int tr_size = ascending_trajectory.size() + descending_trajectory.size();
        SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
        int* trajectory_ptr = INTEGER(trajectory);

        // Copy elements using std::copy instead of manually copying them using for loop
        std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
        std::copy(ascending_trajectory.begin(), ascending_trajectory.end(), trajectory_ptr + descending_trajectory.size());

        SET_VECTOR_ELT(trajectories, vertex, trajectory);
        UNPROTECT(1);
    }

    // Allocate memory for the result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 2)); nprot++;
    // Set elements of the result list
    SET_VECTOR_ELT(result_list, 0, local_extrema);
    SET_VECTOR_ELT(result_list, 1, trajectories);
    // Add names
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result_list;
}



//
// version 2 utilizing a lambda fuction for creating trajectories
//
SEXP S_graph_hop_gradient_flow_trajectories_optimized(SEXP s_graph, SEXP s_Ey, SEXP s_hop_radius) {

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices) {
        error("The lengths of s_graph and s_Ey must be the same");
    }

    if (TYPEOF(s_graph) != VECSXP || TYPEOF(s_Ey) != REALSXP || TYPEOF(s_hop_radius) != INTSXP) {
        error("Invalid input types");
    }

    if (LENGTH(s_hop_radius) != 1) {
        error("hop_radius must be a single integer");
    }

    int hop_radius = INTEGER(s_hop_radius)[0];
    if (hop_radius < 0) {
        error("hop_radius must be non-negative");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    double* Ey = REAL(s_Ey);

    /**
     * @brief Lambda function to find all vertices within a specified hop distance from a start vertex.
     *
     * This lambda function performs a breadth-first search (BFS) to identify all vertices
     * that are within a specified number of hops from the start vertex in the given graph.
     *
     * @param[in] graph      A const reference to a vector of vector of integers representing
     *                       the adjacency list of the graph. graph[i] contains the indices
     *                       of vertices adjacent to vertex i.
     * @param[in] start_vertex The index of the starting vertex for the neighborhood search.
     * @param[in] hop_radius   The maximum number of hops to consider in the neighborhood.
     *
     * @return A unique pointer to a vector of integers containing the indices of all
     *         vertices within the specified hop distance from the start vertex.
     *         The start vertex itself is not included in this vector.
     *
     * @note The function uses a breadth-first search approach, which ensures that
     *       vertices are discovered in order of their distance from the start vertex.
     *
     * @warning The function assumes that the graph is well-formed and that all vertex
     *          indices in the graph are valid. It does not perform bounds checking on
     *          the vertex indices.
     *
     * @see S_graph_hop_gradient_flow_trajectories, which uses this function to compute
     *      extended neighborhoods for robust gradient flow trajectory estimation.
     *
     * Example usage:
     * @code
     * auto get_neighborhood = [](const std::vector<std::vector<int>>& graph, int start_vertex, int hop_radius) -> std::unique_ptr<std::vector<int>> {
     *     // Lambda function implementation...
     * };
     *
     * // Using the lambda function
     * auto neighborhood_ptr = get_neighborhood(graph, current_vertex, hop_radius);
     * for (int neighbor : *neighborhood_ptr) {
     *     // Process each neighbor
     * }
     * @endcode
     */
    auto get_neighborhood = [](const std::vector<std::vector<int>>& graph, int start_vertex, int hop_radius) -> std::unique_ptr<std::vector<int>> {

        auto neighborhood = std::make_unique<std::vector<int>>();
        std::vector<bool> visited(graph.size(), false);
        std::queue<std::pair<int, int>> q;  // vertex, distance

        q.push({start_vertex, 0});
        visited[start_vertex] = true;

        while (!q.empty()) {
            int vertex = q.front().first;
            int distance = q.front().second;
            q.pop();

            if (distance > 0) {  // Don't include the start vertex
                neighborhood->push_back(vertex);
            }

            if (distance < hop_radius) {
                for (int neighbor : graph[vertex]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        q.push({neighbor, distance + 1});
                    }
                }
            }
        }

        return neighborhood;
    };

    #if 0
    // this is a different implementation of the previous lambda to test correctness of the previous lambda
    // using std::unordered_set<int> as neighbors container
    auto uset_get_neighbors = [](const std::vector<std::vector<int>>& graph, int start_vertex, int hop_radius) -> std::unique_ptr<std::unordered_set<int>> {

        if (start_vertex < 0 || start_vertex >= graph.size()) {
            Rprintf("start_vertex: %d\n", start_vertex);
            error("start_vertex out of bounds");
        }

        std::unordered_set<int> neighbors(graph[start_vertex].begin(), graph[start_vertex].end());
        std::unordered_set<int> prev_neighbors;
        prev_neighbors.insert(start_vertex);

        for (int hop_counter = 1; hop_counter < hop_radius; hop_counter++) {
            //Rprintf("hop_counter: %d\n", hop_counter);
            std::unordered_set<int> neighbors_boundary = std::unordered_set<int>(neighbors);

            // neighbors_boundary.erase(prev_neighbors.begin(), prev_neighbors.end());
            // Safer way to remove elements
            for (const auto& prev : prev_neighbors) {
                neighbors_boundary.erase(prev);
            }

            prev_neighbors = neighbors;

            // adding kNN's of each vertex of the boundary of neighbors to neighbors
            for (const auto& nbr_vertex : neighbors_boundary) {
                for (const auto& neighbor : graph[nbr_vertex]) {
                    neighbors.insert(neighbor);
                }
            }
        }

        neighbors.erase(start_vertex);

        return std::make_unique<std::unordered_set<int>>(std::move(neighbors));
    };

    // A block testing if get_neighborhood and uset_get_neighbors functions generate the same sets of vertex indices for each vertex
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        auto neighborhood_ptr = get_neighborhood(graph, vertex, hop_radius);
        auto neighborhs_ptr = uset_get_neighbors(graph, vertex, hop_radius);

        // Convert *neighborhood_ptr to std::unordered_set<int>
        std::unordered_set<int> neighborhood_set(neighborhood_ptr->begin(), neighborhood_ptr->end());

        // Check if the sets are the same
        if (neighborhood_set != *neighborhs_ptr) {
            // If they're not the same, prepare an error message
            std::string error_msg = "Mismatch for vertex " + std::to_string(vertex) + ":\n";

            // Show content of neighborhood_set
            error_msg += "get_neighborhood: { ";
            for (int n : neighborhood_set) {
                error_msg += std::to_string(n) + " ";
            }
            error_msg += "}\n";

            // Show content of *neighborhs_ptr
            error_msg += "uset_get_neighbors: { ";
            for (int n : *neighborhs_ptr) {
                error_msg += std::to_string(n) + " ";
            }
            error_msg += "}\n";

            // Calculate and show the difference
            error_msg += "Difference:\n";
            error_msg += "In get_neighborhood but not in uset_get_neighbors: { ";
            for (int n : neighborhood_set) {
                if (neighborhs_ptr->find(n) == neighborhs_ptr->end()) {
                    error_msg += std::to_string(n) + " ";
                }
            }
            error_msg += "}\n";
            error_msg += "In uset_get_neighbors but not in get_neighborhood: { ";
            for (int n : *neighborhs_ptr) {
                if (neighborhood_set.find(n) == neighborhood_set.end()) {
                    error_msg += std::to_string(n) + " ";
                }
            }
            error_msg += "}\n";

            // Throw an error with the prepared message
            Rprintf("ERROR: %s\n", error_msg.c_str());
        }
    }
    #endif

    std::vector<std::vector<int>> neighborhoods(n_vertices);
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        auto neighborhood_ptr = get_neighborhood(graph, vertex, hop_radius);
        neighborhoods[vertex] = *neighborhood_ptr;
    }


    /**
     * @brief Creates a function for computing gradient flow trajectories.
     *
     * This higher-order lambda function generates another lambda function that computes
     * either ascending or descending trajectories in the graph based on the provided comparator.
     * It utilizes pre-computed neighborhoods for efficient trajectory calculation.
     *
     * @param[in] graph A const reference to a vector of vectors representing the graph's adjacency list.
     * @param[in] Ey A pointer to an array of doubles representing the values at each vertex.
     * @param[in] hop_radius The maximum number of hops to consider in the neighborhood (not used directly in this function, but relevant for context).
     * @param[in] comparator A function object defining the comparison criteria for determining the trajectory direction.
     *
     * @return A lambda function that takes a start vertex and returns a vector of integers representing the trajectory.
     *
     * @details
     * The outer lambda, create_trajectory_function, is a factory function that creates and returns an inner lambda.
     * This inner lambda is the actual trajectory computation function.
     *
     * The inner lambda:
     * - Takes a single parameter: the start_vertex (int).
     * - Returns a std::vector<int> representing the computed trajectory.
     * - Captures by reference:
     *   - neighborhoods: The pre-computed neighborhoods for each vertex.
     *   - Ey: The array of values at each vertex.
     *   - comparator: The comparison function for determining trajectory direction.
     *
     * The inner lambda computes the trajectory by iteratively finding the next vertex
     * in the neighborhood that satisfies the comparator condition, until a local extremum is reached.
     *
     * Usage example:
     * @code
     * auto ascending_trajectory = create_trajectory_function(graph, Ey, hop_radius,
     *                                                        [](double a, double b) { return a < b; });
     * std::vector<int> trajectory = ascending_trajectory(start_vertex);
     * @endcode
     *
     * @note This function assumes that the 'neighborhoods' vector is already computed and available in the enclosing scope.
     * @warning The inner lambda captures 'neighborhoods' by reference. Ensure 'neighborhoods' remains valid for the lifetime of the returned function.
     */
    auto create_trajectory_function = [&neighborhoods](const std::vector<std::vector<int>>& graph,
                                                       const double* Ey,
                                                       int hop_radius,
                                                       const std::function<bool(double, double)>& comparator) {
        return [&neighborhoods, Ey, comparator](int start_vertex) -> std::vector<int> {
            std::vector<int> trajectory = {start_vertex};
            int current_vertex = start_vertex;

            while (true) {
                auto& neighborhood = neighborhoods[current_vertex];
                if (neighborhood.empty()) {
                    break;  // No neighbors, end the trajectory
                }
                auto next_it = std::max_element(neighborhood.begin(), neighborhood.end(),
                                                [Ey, comparator](int a, int b) { return comparator(Ey[a], Ey[b]); });

                if (next_it == neighborhood.end() || !comparator(Ey[current_vertex], Ey[*next_it])) {
                    break; // Local extremum found
                }

                int next_vertex = *next_it;
                trajectory.push_back(next_vertex);
                current_vertex = next_vertex;
            }

            return trajectory;
        };
    };

    auto ascending_trajectory = create_trajectory_function(graph, Ey, hop_radius,
                                                           [](double a, double b) { return a < b; });

    auto descending_trajectory = create_trajectory_function(graph, Ey, hop_radius,
                                                            [](double a, double b) { return a > b; });

    int nprot = 0;
    // Create list elements
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    SEXP local_extrema = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    int* local_extrema_ptr = INTEGER(local_extrema);

    // Main computation
    for (int vertex = 0; vertex < n_vertices; vertex++) {

        // Finding the ascending and descending trajectoris of the 'vertex'
        std::vector<int> asc_traj = ascending_trajectory(vertex);
        std::vector<int> desc_traj = descending_trajectory(vertex);

        // Store local extrema in the matrix
        local_extrema_ptr[vertex] = desc_traj.back();
        local_extrema_ptr[vertex + n_vertices] = asc_traj.back();

        // Trajectory
        // Remove the first element i from desc_traj as it is
        // already in asc_traj
        desc_traj.erase(desc_traj.begin());

        // Create an R integer vector for the trajectory
        int tr_size = asc_traj.size() + desc_traj.size();
        SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
        int* trajectory_ptr = INTEGER(trajectory);

        // Copy elements using std::copy instead of manually copying them using for loop
        std::copy(desc_traj.rbegin(), desc_traj.rend(), trajectory_ptr);
        std::copy(asc_traj.begin(), asc_traj.end(), trajectory_ptr + desc_traj.size());

        SET_VECTOR_ELT(trajectories, vertex, trajectory);
        UNPROTECT(1);
    }

    // Allocate memory for the result list
    SEXP result_list = PROTECT(allocVector(VECSXP, 2)); nprot++;
    // Set elements of the result list
    SET_VECTOR_ELT(result_list, 0, local_extrema);
    SET_VECTOR_ELT(result_list, 1, trajectories);
    // Add names
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result_list;
}

/**
 * @brief Computes constrained gradient flow trajectories on a graph while preventing basin-jumping
 *
 * @param s_graph SEXP containing h-th power graph pIG_k^h(X) as adjacency lists
 * @param s_core_graph SEXP containing core graph pIG_k(X) as adjacency lists
 * @param s_Ey SEXP containing function values at vertices
 *
 * @return SEXP List of integer vectors, each containing complete trajectory
 *         (descending path reversed + ascending path) for corresponding vertex
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

SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {
    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);

    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices));

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
        }

        if (is_max) {
            // Single vertex trajectory for local maximum
            SEXP trajectory = PROTECT(allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
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

    UNPROTECT(1);
    return trajectories;
}


/**
 * Constructs the Morse-Smale complex of an estimate of a conditional expectation of a random variable over a graph.
 *
 * This function identifies the basins of attraction for local minima and local
 * maxima in the gradient flow of a random variable over a graph, and computes
 * connected components of the sub-graph derived from vertiecs of all
 * trajectories with the same local minimum and local maximum. These connected
 * components define Morse-Smale cells of the Morse-Smale complex.
 *
 * The function first calls `S_graph_gradient_flow_trajectories` to obtain the gradient flow trajectories for each
 * vertex in the graph. It then constructs the Morse-Smale complex cell assignments by mapping unique pairs of
 * local minimum and local maximum indices to cell IDs using an `std::unordered_map`.
 *
 * The local extrema information and cell assignments are combined into a single matrix ("MS_complex_info") with
 * columns for local minimum indices, local maximum indices, cell IDs, and cell connected component IDs. Column names
 * are assigned to this matrix for clarity.
 *
 * The function constructs the adjacency lists of the cells of the Morse-Smale complex by iterating over each vertex
 * and its neighbors, adding the neighbors belonging to the same cell to the corresponding cell subgraph.
 *
 * Additionally, the function constructs the adjacency list of the Morse-Smale directed graph, where the vertices are
 * local minima and local maxima, and the edges represent the connections between them based on the Morse-Smale cells.
 * The number of edges connecting a local minimum and a local maximum is equal to the number of connected components
 * of the set difference of the set of all vertices of the cell minus the local minimum and local maximum of the cell.
 *
 * @param s_graph An R list representing the graph adjacency list. Each element of the list should be an integer vector
 *               specifying the neighboring vertex indices for a given vertex. The vertex indices should be 1-based.
 * @param s_Ey    An R numeric vector containing estimates of a conditional expectation of a random variable over the
 *               graph vertices.
 * @param s_hop_radius An R integer specifying the number of hops to consider for the
 *                  neighborhood of each vertex. Default is 1, which is equivalent
 *                  to the original algorithm considering only immediate neighbors.
 *
 * @param s_method An R integer specifying the method to use for gradient flow trajectory calculation:
 *                 - 0: Constrained method (S_graph_constrained_gradient_flow_trajectories)
 *                 - 1: Original method (S_graph_gradient_flow_trajectories)
 *                 - 2: Hop-based method (S_graph_hop_gradient_flow_trajectories)
 *                 - 3: Optimized hop-based method (S_graph_hop_gradient_flow_trajectories_optimized)
 *
 * @return An R list with the following components:
 *         - "MS_cx": An integer matrix with dimensions (n_vertices, 4), where the first column contains the local
 *                    minimum indices, the second column contains the local maximum indices, the third column contains
 *                    the Morse-Smale complex cell assignments for each vertex, and the fourth column contains the
 *                    cell connected component IDs.
 *         - "trajectories": A list of integer vectors representing the gradient flow trajectories for each vertex
 *                           of the graph. Each trajectory includes both the descending and ascending paths from the
 *                           vertex to its local minimum and maximum, respectively.
 *         - "MS_cell_graphs": A list of the graphs (adjacency lists) of the cells of the constructed Morse-Smale complex.
 *         - "MS_graph": An R list representing the adjacency list of the Morse-Smale directed graph, where the vertices
 *                       are local minima and local maxima, and the edges represent the connections between them based on
 *                       the Morse-Smale cells.
 */
SEXP S_graph_MS_cx(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey, SEXP s_hop_radius, SEXP s_method) {

#define DEBUG__S_graph_MS_cx 0

    // Check if the input is a list
    if (!isNewList(s_graph))
        error("s_graph must be a list");

    int method = INTEGER(s_method)[0];
    if (method < 0) {
        error("method must be non-negative");
    } else if (method > 3) {
        error("method must be less than 3");
    }

    // Gradient flow method: 0 = original, 1 = hop-based, 2 = optimized hop-based
    enum gradient_flow_method_t {
        CONSTRAINED = 0,
        ORIGINAL = 1,
        HOP_BASED = 2,
        OPTIMIZED = 3
    };

    gradient_flow_method_t gradient_method = static_cast<gradient_flow_method_t>(method);

    if ((gradient_method == HOP_BASED || gradient_method == OPTIMIZED) && s_hop_radius == R_NilValue) {
        error("s_hop_radius must be provided for hop-based methods");
    }

    int nprot = 0;
    SEXP grad_flow_trajectories_Rlist;
    switch(gradient_method) {
    case CONSTRAINED:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_constrained_gradient_flow_trajectories(s_graph, s_core_graph, s_Ey));
        UNPROTECT(1);
        break;
    case ORIGINAL:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_gradient_flow_trajectories(s_graph, s_Ey));
        UNPROTECT(1);
        break;
    case HOP_BASED:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_hop_gradient_flow_trajectories(s_graph, s_Ey, s_hop_radius));
        UNPROTECT(1);
        break;
    case OPTIMIZED:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_hop_gradient_flow_trajectories_optimized(s_graph, s_Ey, s_hop_radius));
        UNPROTECT(1);
        break;
    default:
        error("Unexpected method value: method = %d. method can be only 0, 1, 2 or 3.", method);
    }

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);

    SEXP lext = VECTOR_ELT(grad_flow_trajectories_Rlist, 0);
    int* lext_ptr = INTEGER(lext);
    SEXP trajectory_Rlist = VECTOR_ELT(grad_flow_trajectories_Rlist, 1);

    int n_vertices = INTEGER(getAttrib(lext, R_DimSymbol))[0];

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the MS_cells integer vector of the Morse-Smale complex cell
    // assignments.
    //
    // -------------------------------------------------------------------------------------------------------------------
    std::unordered_map<std::pair<int, int>, int, pair_hash> precell_lminlmax_to_precell_index_map;
    std::unordered_map<int, std::pair<int,int>> precell_index_to_precell_lminlmax_map;

    SEXP MS_precells = PROTECT(allocVector(INTSXP, n_vertices)); nprot++;
    int* MS_precells_ptr = INTEGER(MS_precells);
    int precell_id = 0;

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        int lmin = lext_ptr[vertex];
        int lmax = lext_ptr[vertex + n_vertices];
        std::pair<int, int> key(lmin, lmax);

        auto it = precell_lminlmax_to_precell_index_map.find(key);
        if (it == precell_lminlmax_to_precell_index_map.end()) { // if precell_lminlmax_to_precell_index_map has no value for 'key', we add the value (precell_lminlmax_to_precell_index_map[key] = precell_id;) and increment 'precell_id'
            precell_lminlmax_to_precell_index_map[key] = precell_id;
            MS_precells_ptr[vertex] = precell_id;
            precell_index_to_precell_lminlmax_map[precell_id] = key;
            precell_id++;
        } else {
            MS_precells_ptr[vertex] = it->second;
        }
    }

    int n_precells = precell_id;

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Creating an unordered_map that maps a __trajectory precell__ to the set of its vertices
    //
    // -------------------------------------------------------------------------------------------------------------------

    std::unordered_map<int, std::set<int>> trajectory_precell_vertices_map;

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        int precell = MS_precells_ptr[vertex];
        SEXP trajectory = VECTOR_ELT(trajectory_Rlist, vertex);
        int* trajectory_ptr = INTEGER(trajectory);
        int trajectory_length = LENGTH(trajectory);  // Get the length of the trajectory vector

        // Insert the vertices of trajectory_ptr into trajectory_precell_vertices_map[precell] set
        for (int i = 0; i < trajectory_length; ++i)
            trajectory_precell_vertices_map[precell].insert(trajectory_ptr[i]);
    }

    //SEXP trajectory_precell_vertices_Rlist = Cpp_map_int_set_int_to_Rlist(trajectory_precell_vertices_map);

    // Creating a precell vertex index maps. Suppose {32, 42, 57} is the sorted set
    // of the vertices of one of the MS precells. This precell's index map maps
    // 32 to 0,
    // 42 to 1 and
    // 57 to 2
    // This map will be utilized in the construction of MS precell graph
    std::unordered_map<int, std::map<int,int>> trajectory_precell_vertices_index_map;
    for (const auto& [precell, vertex_set] : trajectory_precell_vertices_map) {
        int index = 0;
        for (const auto& vertex : vertex_set) {
            trajectory_precell_vertices_index_map[precell][vertex] = index++;
        }
    }

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing a Morse-Smale complex graph. The vertices of the graph are
    // local minima and local maxima of Ey. A local minimum, lmin_i, is
    // connected with a local maximum, lmax_j, if there is a Morse-Smale precell
    // (lmin_i, lmax_j). If such precell exists, the number of directed edges
    // connecting lmin_i with lmax_j is equal to the number of connected
    // components of the set difference of the set of all vertices of the precell
    // minus the local miminum and local maximum of the precell.
    //
    // -------------------------------------------------------------------------------------------------------------------

    // Creating indices of local minima and local maxima
    std::set<int> lmin_set;
    std::set<int> lmax_set;
    for (int precell = 0; precell < n_precells; precell++) {
        auto lminlmax_pair = precell_index_to_precell_lminlmax_map[precell];
        lmin_set.insert(lminlmax_pair.first);
        lmax_set.insert(lminlmax_pair.second);
    }

    // Creating indices of local minima and local maxima
    std::map<int,int> lmin_map;
    std::map<int,int> lmax_map;
    int index = 0;
    for (const auto& e : lmin_set)
        lmin_map[e] = index++;
    for (const auto& e : lmax_set)
        lmax_map[e] = index++;

    // Creating a vector of MS_graph vertex labels
    std::vector<int> MS_graph_vertex_label = std::vector<int>(lmin_set.size() + lmax_set.size());
    for (const auto& [key,val] : lmin_map)
        MS_graph_vertex_label[val] = key;
    for (const auto& [key,val] : lmax_map)
        MS_graph_vertex_label[val] = key;

    SEXP MS_graph_vertex_label_Rlist = PROTECT(Cpp_vector_int_to_Rvector(MS_graph_vertex_label)); nprot++;


    // MS graph is a bipartite graph; assigning lmin vertices type index 1 and lmax vertices type index 2
    std::vector<int> MS_graph_vertex_type = std::vector<int>(lmin_set.size() + lmax_set.size());
    for (const auto& [k,v] : lmin_map)
        MS_graph_vertex_type[v] = 1;
    for (const auto& [k,v] : lmax_map)
        MS_graph_vertex_type[v] = 2;

    SEXP MS_graph_vertex_type_Rlist = PROTECT(Cpp_vector_int_to_Rvector(MS_graph_vertex_type)); nprot++;

    // Finding connected components of the set difference of the set of all
    // vertices of the precell minus the local miminum and local maximum of the
    // precell.

    // Creating an adjacency list of the Morse-Smale directed graph in the form
    // of a unordered map first
    std::unordered_map<int, std::vector<int>> MS_graph_map;
    std::unordered_map<int, std::vector<int>> MS_graph_edge_label; // edge labels are indices of precell connected components
    // Creating a map of each components vertices
    std::unordered_map<int, std::set<int>> MS_cell_vertices_map;
    int precell_cc_index = 0;
    for (int precell = 0; precell < n_precells; precell++) {

        std::set<int> precell_vertices = trajectory_precell_vertices_map[precell];

        // Subtracting the corresponding local mimimum and local maximum vertices from this set
        auto lminlmax_pair = precell_index_to_precell_lminlmax_map[precell];

#if DEBUG__S_graph_MS_cx
        Rprintf("precell: %d\tlmin: %d\tlmax: %d\n", precell, lminlmax_pair.first, lminlmax_pair.second);
        print_set(precell_vertices, "precell_vertices");
#endif

        MS_graph_map[lmin_map[lminlmax_pair.first]].push_back(lmax_map[lminlmax_pair.second]); // adding an edge from lmin to lmax
        MS_graph_edge_label[lmin_map[lminlmax_pair.first]].push_back(precell_cc_index + 1);

        precell_vertices.erase(lminlmax_pair.first);  // removing local minimum vertex
        precell_vertices.erase(lminlmax_pair.second); // removing local maximum vertex

#if DEBUG__S_graph_MS_cx
        print_set(precell_vertices, "precell_vertices AFTER ERASE");
#endif

        // Finding connected components of the resulting set of vertices within the graph
        std::vector<int> precell_vertices_vect(precell_vertices.begin(), precell_vertices.end()); // in the future I can create a version of count_subgraph_components that accepts std::set<int> as the seond parameter
        std::unique_ptr<std::unordered_map<int, int>> precell_comps_uptr = count_subgraph_components(graph, precell_vertices_vect);
#if DEBUG__S_graph_MS_cx
        print_umap(*precell_comps_uptr, "*precell_comps_uptr");
#endif
        // Creating a map of each components vertices
        std::unordered_map<int, std::set<int>> local_MS_cell_vertices_map;
        for (const auto& [vertex, comp] : *precell_comps_uptr) {
            local_MS_cell_vertices_map[comp].insert(vertex);
            //MS_precell_cc[vertex].insert(precell_cc_index++);
        }

        std::unordered_map<int, int> local_to_global_cc_map;
        for (const auto& [comp, set] : local_MS_cell_vertices_map) {
            MS_cell_vertices_map[precell_cc_index].insert(set.begin(), set.end());
            MS_cell_vertices_map[precell_cc_index].insert(lminlmax_pair.first);
            MS_cell_vertices_map[precell_cc_index].insert(lminlmax_pair.second);
            local_to_global_cc_map[comp] = precell_cc_index++;
        }

        // If the number of connected components of precell_vertices_vect is more
        // than 1, insert additional edges to the graph so that the number of
        // edges connectint the given (lmin, lmax) pair is equal to the number
        // of the connected components of precell_vertices_vect
        int cc_counter = 0;
        for (const auto& [local_cc_index, global_cc_index] : local_to_global_cc_map) {
            if ( cc_counter > 0 ) {
                MS_graph_map[lmin_map[lminlmax_pair.first]].push_back(lmax_map[lminlmax_pair.second]);
                MS_graph_edge_label[lmin_map[lminlmax_pair.first]].push_back(global_cc_index + 1);
            }
            cc_counter++;
        }
    }

#if DEBUG__S_graph_MS_cx
    // print_umap_to_vect(MS_graph_map, "MS_graph_map");
    print_umap_to_set(MS_cell_vertices_map, "MS_cell_vertices_map");
    // print_vect(MS_graph_vertex_label, "MS_graph_vertex_label");
    // print_vect(MS_graph_vertex_type, "MS_graph_vertex_type");
    //error("Stopping function execution for debugging purposes in S_graph_MS_cx()");
#endif


    // Creating an R adjacency list corresponding to MS_graph_map
    SEXP MS_graph_Rlist = PROTECT(Cpp_map_int_vector_int_to_Rlist(MS_graph_map)); nprot++;

    // Creating an R list of edge labels
    SEXP MS_graph_edge_label_Rlist = PROTECT(Cpp_map_int_vector_int_to_Rlist(MS_graph_edge_label)); nprot++;

    // Creating an R list of trajectory precell's connected components
    SEXP MS_cell_vertices_Rlist = PROTECT(Cpp_map_int_set_int_to_Rlist(MS_cell_vertices_map)); nprot++;


    // ------------------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the MS_cell_graphs list of the graphs using MS cells of the constructed Morse-Smale complex
    //
    // ------------------------------------------------------------------------------------------------------------------------------

    int n_cells = precell_cc_index;

    // Creating a precell vertex index maps. Suppose {32, 42, 57} is the sorted set
    // of the vertices of one of the MS precells. This precell's index map maps
    // 32 to 0,
    // 42 to 1 and
    // 57 to 2
    // This map will be utilized in the construction of MS precell graph
    std::unordered_map<int, std::map<int,int>> MS_cell_vertices_index_map;
    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        int index = 0;
        for (const auto& vertex : vertex_set) {
            MS_cell_vertices_index_map[precell][vertex] = index++;
        }
    }

    // - For each vertex `precell` and the corresponding `vertex_set`, iterate over the vertices of the `vertex_set`
    // - Iterate over each neighbor of vertex `vertex`
    // - If the neighbor belongs to the same precell as vertex `vertex` (checked by `vertex_set.find(neighbor) != vertex_set.end()`), add the neighbor to the adjacency list of vertex `vertex` within the corresponding precell subgraph.
    std::vector<std::vector<std::vector<int>>> MS_cell_graphs(n_cells, std::vector<std::vector<int>>(n_vertices)); // `MS_cell_graphs` is initialized as a vector of size `n_cells`, where each element is a vector of size `n_vertices`. This ensures that `MS_cell_graphs[precell]` has a size equal to the number of vertices in the graph.

    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        for (const auto& vertex : vertex_set) {
            for (int neighbor : graph[vertex]) {
                if (vertex_set.find(neighbor) != vertex_set.end()) {
                    int vertex_i = MS_cell_vertices_index_map[precell][vertex];
                    int neighbor_i = MS_cell_vertices_index_map[precell][neighbor];
                    MS_cell_graphs[precell][vertex_i].push_back(neighbor_i);
                }
            }
        }
    }

    SEXP MS_cell_graphs_Rlist = PROTECT(allocVector(VECSXP, n_cells)); nprot++;

    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        int vertex_set_size = vertex_set.size();
        SEXP cell_graph = PROTECT(allocVector(VECSXP, vertex_set_size));

        for (int j = 0; j < vertex_set_size; j++) {
            int n_neighbors = MS_cell_graphs[precell][j].size();
            SEXP neighbors = PROTECT(allocVector(INTSXP, n_neighbors));
            int* neighbors_ptr = INTEGER(neighbors);

            for (int k = 0; k < n_neighbors; k++) {
                neighbors_ptr[k] = MS_cell_graphs[precell][j][k] + 1; // Convert to 1-based index
            }

            SET_VECTOR_ELT(cell_graph, j, neighbors);
            UNPROTECT(1);
        }

        if (precell < 0 || precell >= n_cells) {
            Rprintf("precell: %d", precell);
            Rprintf("n_cells: %d", n_cells);
            error("precell >= n_cells");
        }
        SET_VECTOR_ELT(MS_cell_graphs_Rlist, precell, cell_graph);
        UNPROTECT(1);
    }

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Combining lext, MS_precells and precell connected components into a single matrix
    //
    // -------------------------------------------------------------------------------------------------------------------

    SEXP MS_complex_info = PROTECT(allocMatrix(INTSXP, n_vertices, 4)); nprot++;
    int* MS_complex_info_ptr = INTEGER(MS_complex_info);

    for (int i = 0; i < n_vertices; i++) {
        MS_complex_info_ptr[i]                  = lext_ptr[i] + 1;               // lmin as a 1-based integer
        MS_complex_info_ptr[i + n_vertices]     = lext_ptr[i + n_vertices] + 1;  // lmax as a 1-based integer
        MS_complex_info_ptr[i + 2 * n_vertices] = MS_precells_ptr[i] + 1;        // cell ID as a 1-based integer
        MS_complex_info_ptr[i + 3 * n_vertices] = 0;                             // placeholder for connected component ID
    }

    // Update column names
    SEXP colnames = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(colnames, 0, mkChar("local_min"));
    SET_STRING_ELT(colnames, 1, mkChar("local_max"));
    SET_STRING_ELT(colnames, 2, mkChar("cell_id"));
    SET_STRING_ELT(colnames, 3, mkChar("component_id"));

    #if 0
    SEXP MS_complex_info = PROTECT(allocMatrix(INTSXP, n_vertices, 2)); nprot++;
    int* MS_complex_info_ptr = INTEGER(MS_complex_info);

    for (int i = 0; i < n_vertices; i++) {
        MS_complex_info_ptr[i]              = lext_ptr[i] + 1;               // lmin as a 1-based integer
        MS_complex_info_ptr[i + n_vertices] = lext_ptr[i + n_vertices] + 1;  // lmax as a 1-based integer
    }

    // Assigning column names to the MS_complex_info matrix
    SEXP colnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(colnames, 0, mkChar("local_min"));
    SET_STRING_ELT(colnames, 1, mkChar("local_max"));
    #endif

    // Wrapping colnames in a list
    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue); // No row names
    SET_VECTOR_ELT(dimnames, 1, colnames);

    // Assigning the dimnames to the matrix
    setAttrib(MS_complex_info, R_DimNamesSymbol, dimnames);

    // Unprotecting colnames and dimnames as they are no longer needed
    UNPROTECT(1); // unprotect dimnames
    UNPROTECT(1); // unprotect colnames

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the result list
    //
    // -------------------------------------------------------------------------------------------------------------------

    int result_len = 8;
    SEXP result = PROTECT(allocVector(VECSXP, result_len)); nprot++;

    SET_VECTOR_ELT(result, 0, MS_complex_info);
    SET_VECTOR_ELT(result, 1, MS_graph_Rlist);
    SET_VECTOR_ELT(result, 2, MS_graph_vertex_label_Rlist);
    SET_VECTOR_ELT(result, 3, MS_graph_vertex_type_Rlist);
    SET_VECTOR_ELT(result, 4, MS_graph_edge_label_Rlist);
    SET_VECTOR_ELT(result, 5, MS_cell_graphs_Rlist);
    //SET_VECTOR_ELT(result, 6, trajectory_precell_vertices_Rlist);
    SET_VECTOR_ELT(result, 6, MS_cell_vertices_Rlist);
    SET_VECTOR_ELT(result, 7, trajectory_Rlist);

    // Add names to list elements
    SEXP names = PROTECT(allocVector(STRSXP, result_len)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("MS_cx"));
    SET_STRING_ELT(names, 1, mkChar("MS_graph"));
    SET_STRING_ELT(names, 2, mkChar("MS_graph_labels"));
    SET_STRING_ELT(names, 3, mkChar("MS_graph_types"));
    SET_STRING_ELT(names, 4, mkChar("MS_graph_edge_label"));
    SET_STRING_ELT(names, 5, mkChar("MS_cell_graphs"));
    //SET_STRING_ELT(names, 6, mkChar("MS_precell_vertices"));
    SET_STRING_ELT(names, 6, mkChar("MS_cell_vertices"));
    SET_STRING_ELT(names, 7, mkChar("trajectories"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result;
}

/**
 * @brief Constructs the Morse-Smale complex of an estimate of a conditional expectation of a random variable over a graph, with a minimum cell size constraint.
 *
 * This function identifies the basins of attraction for local minima and local maxima in the gradient flow of a random
 * variable over a graph, and computes connected components of the sub-graph derived from vertices of all trajectories
 * with the same local minimum and local maximum. These connected components define Morse-Smale cells of the Morse-Smale
 * complex. The function incorporates a minimum cell size constraint, merging or eliminating cells smaller than the
 * specified threshold.
 *
 * The function first calls a gradient flow trajectory calculation method to obtain the gradient flow trajectories for
 * each vertex in the graph. It then constructs the Morse-Smale complex cell assignments, filtering out cells smaller
 * than the minimum size and merging their vertices with neighboring valid cells.
 *
 * The local extrema information and cell assignments are combined into a single matrix ("MS_complex_info") with columns
 * for local minimum indices, local maximum indices, cell IDs, and cell connected component IDs.
 *
 * The function constructs the adjacency lists of the cells of the Morse-Smale complex and the adjacency list of the
 * Morse-Smale directed graph, considering only valid cells that meet the minimum size requirement.
 *
 * @param s_graph An R list representing the graph adjacency list. Each element of the list should be an integer vector
 *                specifying the neighboring vertex indices for a given vertex. The vertex indices should be 1-based.
 * @param s_Ey An R numeric vector containing estimates of a conditional expectation of a random variable over the
 *             graph vertices.
 * @param s_hop_radius An R integer specifying the number of hops to consider for the neighborhood of each vertex.
 *                     Default is 1, which is equivalent to the original algorithm considering only immediate neighbors.
 * @param s_method An R integer specifying the method to use for gradient flow trajectory calculation:
 *                 - 0: Original method (S_graph_gradient_flow_trajectories)
 *                 - 1: Hop-based method (S_graph_hop_gradient_flow_trajectories)
 *                 - 2: Optimized hop-based method (S_graph_hop_gradient_flow_trajectories_optimized)
 * @param s_min_cell_size An R integer specifying the minimum number of vertices a cell must contain to be considered valid.
 *                        Cells smaller than this threshold will be merged with neighboring cells or eliminated.
 *
 * @return An R list with the following components:
 *         - "MS_cx": An integer matrix with dimensions (n_vertices, 4), where the first column contains the local
 *                    minimum indices, the second column contains the local maximum indices, the third column contains
 *                    the Morse-Smale complex cell assignments for each vertex, and the fourth column contains the
 *                    cell connected component IDs. Only includes information for valid cells meeting the minimum size requirement.
 *         - "MS_graph": An R list representing the adjacency list of the Morse-Smale directed graph, where the vertices
 *                       are local minima and local maxima of valid cells, and the edges represent the connections between
 *                       them based on the Morse-Smale cells.
 *         - "MS_graph_labels": An R integer vector containing the labels (original vertex indices) of the vertices in the MS_graph.
 *         - "MS_graph_types": An R integer vector indicating the type (1 for local minimum, 2 for local maximum) of each vertex in the MS_graph.
 *         - "MS_graph_edge_label": An R list of integer vectors, where each vector contains the cell IDs associated with the edges in MS_graph.
 *         - "MS_cell_graphs": A list of the graphs (adjacency lists) of the valid cells of the constructed Morse-Smale complex.
 *         - "MS_cell_vertices": An R list of integer vectors, where each vector contains the vertices belonging to a valid cell.
 *         - "trajectories": A list of integer vectors representing the gradient flow trajectories for each vertex
 *                           of the graph. Each trajectory includes both the descending and ascending paths from the
 *                           vertex to its local minimum and maximum, respectively.
 *         - "cell_sizes": An R integer vector containing the sizes (number of vertices) of all valid cells.
 *
 * @note This function modifies the original Morse-Smale complex construction by incorporating a minimum cell size constraint.
 *       Cells smaller than the specified threshold are merged with neighboring cells or eliminated, potentially altering
 *       the topology of the resulting complex.
 */
SEXP S_graph_MS_cx_min_cell(SEXP s_graph, SEXP s_Ey, SEXP s_hop_radius, SEXP s_method, SEXP s_min_cell_size) {

    // Check if the input is a list
    if (!isNewList(s_graph))
        error("s_graph must be a list");

    int method = INTEGER(s_method)[0];
    if (method < 0) {
        error("method must be non-negative");
    } else if (method > 2) {
        error("method must be less than 3");
    }

    int min_cell_size = INTEGER(s_min_cell_size)[0];

        // Gradient flow method: 0 = original, 1 = hop-based, 2 = optimized hop-based
    enum gradient_flow_method_t {
        ORIGINAL = 0,
        HOP_BASED = 1,
        OPTIMIZED = 2
    };

    gradient_flow_method_t gradient_method = static_cast<gradient_flow_method_t>(method);

    if ((gradient_method == HOP_BASED || gradient_method == OPTIMIZED) && s_hop_radius == R_NilValue) {
        error("s_hop_radius must be provided for hop-based methods");
    }

    int nprot = 0;
    SEXP grad_flow_trajectories_Rlist;
    switch(gradient_method) {
    case ORIGINAL:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_gradient_flow_trajectories(s_graph, s_Ey)); nprot++;
        break;
    case HOP_BASED:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_hop_gradient_flow_trajectories(s_graph, s_Ey, s_hop_radius)); nprot++;
        break;
    case OPTIMIZED:
        grad_flow_trajectories_Rlist = PROTECT(S_graph_hop_gradient_flow_trajectories_optimized(s_graph, s_Ey, s_hop_radius)); nprot++;
        break;
    default:
        error("Unexpected method value: method = %d. method can be only 0, 1 or 2.", method);
    }

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);

    SEXP lext = VECTOR_ELT(grad_flow_trajectories_Rlist, 0);
    int* lext_ptr = INTEGER(lext);
    SEXP trajectory_Rlist = VECTOR_ELT(grad_flow_trajectories_Rlist, 1);

    int n_vertices = INTEGER(getAttrib(lext, R_DimSymbol))[0];

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the MS_cells integer vector of the Morse-Smale complex cell
    // assignments.
    //
    // -------------------------------------------------------------------------------------------------------------------
    // Modify the MS complex construction
    std::unordered_map<std::pair<int, int>, int, pair_hash> precell_lminlmax_to_precell_index_map;
    std::unordered_map<int, std::pair<int,int>> precell_index_to_precell_lminlmax_map;
    std::unordered_map<int, int> cell_sizes;

    SEXP MS_precells = PROTECT(allocVector(INTSXP, n_vertices)); nprot++;
    int* MS_precells_ptr = INTEGER(MS_precells);
    int precell_id = 0;

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        int lmin = lext_ptr[vertex];
        int lmax = lext_ptr[vertex + n_vertices];
        std::pair<int, int> key(lmin, lmax);

        auto it = precell_lminlmax_to_precell_index_map.find(key);
        if (it == precell_lminlmax_to_precell_index_map.end()) {
            precell_lminlmax_to_precell_index_map[key] = precell_id;
            MS_precells_ptr[vertex] = precell_id;
            precell_index_to_precell_lminlmax_map[precell_id] = key;
            cell_sizes[precell_id] = 1;
            precell_id++;
        } else {
            MS_precells_ptr[vertex] = it->second;
            cell_sizes[it->second]++;
        }
    }

    // Filter cells based on min_cell_size
    std::unordered_set<int> valid_cells;
    for (const auto& [cell, size] : cell_sizes) {
        if (size >= min_cell_size) {
            valid_cells.insert(cell);
        }
    }

    // Reassign vertices to valid cells or merge with neighboring cells
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        int current_cell = MS_precells_ptr[vertex];
        if (valid_cells.find(current_cell) == valid_cells.end()) {
            // Find a valid neighboring cell
            int new_cell = -1;
            for (int neighbor : graph[vertex]) {
                int neighbor_cell = MS_precells_ptr[neighbor];
                if (valid_cells.find(neighbor_cell) != valid_cells.end()) {
                    new_cell = neighbor_cell;
                    break;
                }
            }
            if (new_cell != -1) {
                MS_precells_ptr[vertex] = new_cell;
                cell_sizes[new_cell]++;
            }
        }
    }

    int n_precells = precell_id;

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Creating an unordered_map that maps a __trajectory precell__ to the set of its vertices
    //
    // -------------------------------------------------------------------------------------------------------------------

    std::unordered_map<int, std::set<int>> trajectory_precell_vertices_map;

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        int precell = MS_precells_ptr[vertex];
        SEXP trajectory = VECTOR_ELT(trajectory_Rlist, vertex);
        int* trajectory_ptr = INTEGER(trajectory);
        int trajectory_length = LENGTH(trajectory);  // Get the length of the trajectory vector

        // Insert the vertices of trajectory_ptr into trajectory_precell_vertices_map[precell] set
        for (int i = 0; i < trajectory_length; ++i)
            trajectory_precell_vertices_map[precell].insert(trajectory_ptr[i]);
    }

    //SEXP trajectory_precell_vertices_Rlist = Cpp_map_int_set_int_to_Rlist(trajectory_precell_vertices_map);

    // Creating a precell vertex index maps. Suppose {32, 42, 57} is the sorted set
    // of the vertices of one of the MS precells. This precell's index map maps
    // 32 to 0,
    // 42 to 1 and
    // 57 to 2
    // This map will be utilized in the construction of MS precell graph
    std::unordered_map<int, std::map<int,int>> trajectory_precell_vertices_index_map;
    for (const auto& [precell, vertex_set] : trajectory_precell_vertices_map) {
        int index = 0;
        for (const auto& vertex : vertex_set) {
            trajectory_precell_vertices_index_map[precell][vertex] = index++;
        }
    }

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing a Morse-Smale complex graph. The vertices of the graph are
    // local minima and local maxima of Ey. A local minimum, lmin_i, is
    // connected with a local maximum, lmax_j, if there is a Morse-Smale precell
    // (lmin_i, lmax_j). If such precell exists, the number of directed edges
    // connecting lmin_i with lmax_j is equal to the number of connected
    // components of the set difference of the set of all vertices of the precell
    // minus the local miminum and local maximum of the precell.
    //
    // -------------------------------------------------------------------------------------------------------------------

    // Creating indices of local minima and local maxima
    std::set<int> lmin_set;
    std::set<int> lmax_set;
    for (int precell = 0; precell < n_precells; precell++) {
        auto lminlmax_pair = precell_index_to_precell_lminlmax_map[precell];
        lmin_set.insert(lminlmax_pair.first);
        lmax_set.insert(lminlmax_pair.second);
    }

    // Creating indices of local minima and local maxima
    std::map<int,int> lmin_map;
    std::map<int,int> lmax_map;
    int index = 0;
    for (const auto& e : lmin_set)
        lmin_map[e] = index++;
    for (const auto& e : lmax_set)
        lmax_map[e] = index++;

    // Creating a vector of MS_graph vertex labels
    std::vector<int> MS_graph_vertex_label = std::vector<int>(lmin_set.size() + lmax_set.size());
    for (const auto& [key,val] : lmin_map)
        MS_graph_vertex_label[val] = key;
    for (const auto& [key,val] : lmax_map)
        MS_graph_vertex_label[val] = key;

    SEXP MS_graph_vertex_label_Rlist = PROTECT(Cpp_vector_int_to_Rvector(MS_graph_vertex_label)); nprot++;


    // MS graph is a bipartite graph; assigning lmin vertices type index 1 and lmax vertices type index 2
    std::vector<int> MS_graph_vertex_type = std::vector<int>(lmin_set.size() + lmax_set.size());
    for (const auto& [k,v] : lmin_map)
        MS_graph_vertex_type[v] = 1;
    for (const auto& [k,v] : lmax_map)
        MS_graph_vertex_type[v] = 2;

    SEXP MS_graph_vertex_type_Rlist = PROTECT(Cpp_vector_int_to_Rvector(MS_graph_vertex_type)); nprot++;

    // Finding connected components of the set difference of the set of all
    // vertices of the precell minus the local miminum and local maximum of the
    // precell.

    // Creating an adjacency list of the Morse-Smale directed graph in the form
    // of a unordered map first
    std::unordered_map<int, std::vector<int>> MS_graph_map;
    std::unordered_map<int, std::vector<int>> MS_graph_edge_label; // edge labels are indices of precell connected components
    // Creating a map of each components vertices
    std::unordered_map<int, std::set<int>> MS_cell_vertices_map;
    int precell_cc_index = 0;
    for (int precell = 0; precell < n_precells; precell++) {

        std::set<int> precell_vertices = trajectory_precell_vertices_map[precell];

        // Subtracting the corresponding local mimimum and local maximum vertices from this set
        auto lminlmax_pair = precell_index_to_precell_lminlmax_map[precell];

#if DEBUG__S_graph_MS_cx
        Rprintf("precell: %d\tlmin: %d\tlmax: %d\n", precell, lminlmax_pair.first, lminlmax_pair.second);
        print_set(precell_vertices, "precell_vertices");
#endif

        MS_graph_map[lmin_map[lminlmax_pair.first]].push_back(lmax_map[lminlmax_pair.second]); // adding an edge from lmin to lmax
        MS_graph_edge_label[lmin_map[lminlmax_pair.first]].push_back(precell_cc_index + 1);

        precell_vertices.erase(lminlmax_pair.first);  // removing local minimum vertex
        precell_vertices.erase(lminlmax_pair.second); // removing local maximum vertex

#if DEBUG__S_graph_MS_cx
        print_set(precell_vertices, "precell_vertices AFTER ERASE");
#endif

        // Finding connected components of the resulting set of vertices within the graph
        std::vector<int> precell_vertices_vect(precell_vertices.begin(), precell_vertices.end()); // in the future I can create a version of count_subgraph_components that accepts std::set<int> as the seond parameter
        std::unique_ptr<std::unordered_map<int, int>> precell_comps_uptr = count_subgraph_components(graph, precell_vertices_vect);
#if DEBUG__S_graph_MS_cx
        print_umap(*precell_comps_uptr, "*precell_comps_uptr");
#endif
        // Creating a map of each components vertices
        std::unordered_map<int, std::set<int>> local_MS_cell_vertices_map;
        for (const auto& [vertex, comp] : *precell_comps_uptr) {
            local_MS_cell_vertices_map[comp].insert(vertex);
            //MS_precell_cc[vertex].insert(precell_cc_index++);
        }

        std::unordered_map<int, int> local_to_global_cc_map;
        for (const auto& [comp, set] : local_MS_cell_vertices_map) {
            MS_cell_vertices_map[precell_cc_index].insert(set.begin(), set.end());
            MS_cell_vertices_map[precell_cc_index].insert(lminlmax_pair.first);
            MS_cell_vertices_map[precell_cc_index].insert(lminlmax_pair.second);
            local_to_global_cc_map[comp] = precell_cc_index++;
        }

        // If the number of connected components of precell_vertices_vect is more
        // than 1, insert additional edges to the graph so that the number of
        // edges connectint the given (lmin, lmax) pair is equal to the number
        // of the connected components of precell_vertices_vect
        int cc_counter = 0;
        for (const auto& [local_cc_index, global_cc_index] : local_to_global_cc_map) {
            if ( cc_counter > 0 ) {
                MS_graph_map[lmin_map[lminlmax_pair.first]].push_back(lmax_map[lminlmax_pair.second]);
                MS_graph_edge_label[lmin_map[lminlmax_pair.first]].push_back(global_cc_index + 1);
            }
            cc_counter++;
        }
    }

#if DEBUG__S_graph_MS_cx
    // print_umap_to_vect(MS_graph_map, "MS_graph_map");
    print_umap_to_set(MS_cell_vertices_map, "MS_cell_vertices_map");
    // print_vect(MS_graph_vertex_label, "MS_graph_vertex_label");
    // print_vect(MS_graph_vertex_type, "MS_graph_vertex_type");
    //error("Stopping function execution for debugging purposes in S_graph_MS_cx()");
#endif


    // Creating an R adjacency list corresponding to MS_graph_map
    SEXP MS_graph_Rlist = PROTECT(Cpp_map_int_vector_int_to_Rlist(MS_graph_map)); nprot++;

    // Creating an R list of edge labels
    SEXP MS_graph_edge_label_Rlist = PROTECT(Cpp_map_int_vector_int_to_Rlist(MS_graph_edge_label)); nprot++;

    // Creating an R list of trajectory precell's connected components
    SEXP MS_cell_vertices_Rlist = PROTECT(Cpp_map_int_set_int_to_Rlist(MS_cell_vertices_map)); nprot++;


    // ------------------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the MS_cell_graphs list of the graphs using MS cells of the constructed Morse-Smale complex
    //
    // ------------------------------------------------------------------------------------------------------------------------------

    int n_cells = precell_cc_index;

    // Creating a precell vertex index maps. Suppose {32, 42, 57} is the sorted set
    // of the vertices of one of the MS precells. This precell's index map maps
    // 32 to 0,
    // 42 to 1 and
    // 57 to 2
    // This map will be utilized in the construction of MS precell graph
    std::unordered_map<int, std::map<int,int>> MS_cell_vertices_index_map;
    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        int index = 0;
        for (const auto& vertex : vertex_set) {
            MS_cell_vertices_index_map[precell][vertex] = index++;
        }
    }

    // - For each vertex `precell` and the corresponding `vertex_set`, iterate over the vertices of the `vertex_set`
    // - Iterate over each neighbor of vertex `vertex`
    // - If the neighbor belongs to the same precell as vertex `vertex` (checked by `vertex_set.find(neighbor) != vertex_set.end()`), add the neighbor to the adjacency list of vertex `vertex` within the corresponding precell subgraph.
    std::vector<std::vector<std::vector<int>>> MS_cell_graphs(n_cells, std::vector<std::vector<int>>(n_vertices)); // `MS_cell_graphs` is initialized as a vector of size `n_cells`, where each element is a vector of size `n_vertices`. This ensures that `MS_cell_graphs[precell]` has a size equal to the number of vertices in the graph.

    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        for (const auto& vertex : vertex_set) {
            for (int neighbor : graph[vertex]) {
                if (vertex_set.find(neighbor) != vertex_set.end()) {
                    int vertex_i = MS_cell_vertices_index_map[precell][vertex];
                    int neighbor_i = MS_cell_vertices_index_map[precell][neighbor];
                    MS_cell_graphs[precell][vertex_i].push_back(neighbor_i);
                }
            }
        }
    }

    SEXP MS_cell_graphs_Rlist = PROTECT(allocVector(VECSXP, n_cells)); nprot++;

    for (const auto& [precell, vertex_set] : MS_cell_vertices_map) {
        int vertex_set_size = vertex_set.size();
        SEXP cell_graph = PROTECT(allocVector(VECSXP, vertex_set_size));

        for (int j = 0; j < vertex_set_size; j++) {
            int n_neighbors = MS_cell_graphs[precell][j].size();
            SEXP neighbors = PROTECT(allocVector(INTSXP, n_neighbors));
            int* neighbors_ptr = INTEGER(neighbors);

            for (int k = 0; k < n_neighbors; k++) {
                neighbors_ptr[k] = MS_cell_graphs[precell][j][k] + 1; // Convert to 1-based index
            }

            SET_VECTOR_ELT(cell_graph, j, neighbors);
            UNPROTECT(1);
        }

        if (precell < 0 || precell >= n_cells) {
            Rprintf("precell: %d", precell);
            Rprintf("n_cells: %d", n_cells);
            error("precell >= n_cells");
        }
        SET_VECTOR_ELT(MS_cell_graphs_Rlist, precell, cell_graph);
        UNPROTECT(1);
    }

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Combining lext, MS_precells and precell connected components into a single matrix
    //
    // -------------------------------------------------------------------------------------------------------------------

    double* Ey = REAL(s_Ey);

    // Creating a mapping from old cell IDs to new cell IDs
    std::unordered_map<int, int> cell_id_map;
    int new_cell_id = 1;
    for (int cell : valid_cells) {
        cell_id_map[cell] = new_cell_id++;
    }

    // Creating the updated MS_complex_info matrix
    SEXP MS_complex_info = PROTECT(allocMatrix(INTSXP, n_vertices, 4)); nprot++;
    int* MS_complex_info_ptr = INTEGER(MS_complex_info);

    for (int i = 0; i < n_vertices; i++) {
        int current_cell = MS_precells_ptr[i];
        if (valid_cells.find(current_cell) != valid_cells.end()) {
            // Valid cell
            MS_complex_info_ptr[i]                  = lext_ptr[i] + 1;               // lmin as a 1-based integer
            MS_complex_info_ptr[i + n_vertices]     = lext_ptr[i + n_vertices] + 1;  // lmax as a 1-based integer
            MS_complex_info_ptr[i + 2 * n_vertices] = cell_id_map[current_cell];     // new cell ID as a 1-based integer
            MS_complex_info_ptr[i + 3 * n_vertices] = 1;                             // component ID (all 1 for now, update if needed)
        } else {
            // Invalid cell, find the nearest valid cell
            int nearest_valid_cell = -1;
            double min_distance = std::numeric_limits<double>::max();
            for (int valid_cell : valid_cells) {
                double distance = std::abs(Ey[i] - Ey[precell_index_to_precell_lminlmax_map[valid_cell].first]);
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_valid_cell = valid_cell;
                }
            }

            auto [lmin, lmax] = precell_index_to_precell_lminlmax_map[nearest_valid_cell];
            MS_complex_info_ptr[i]                  = lmin + 1;                          // lmin of nearest valid cell
            MS_complex_info_ptr[i + n_vertices]     = lmax + 1;                          // lmax of nearest valid cell
            MS_complex_info_ptr[i + 2 * n_vertices] = cell_id_map[nearest_valid_cell];   // cell ID of nearest valid cell
            MS_complex_info_ptr[i + 3 * n_vertices] = 1;                                 // component ID (all 1 for now, update if needed)
        }
    }

    // Update column names
    SEXP colnames = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(colnames, 0, mkChar("local_min"));
    SET_STRING_ELT(colnames, 1, mkChar("local_max"));
    SET_STRING_ELT(colnames, 2, mkChar("cell_id"));
    SET_STRING_ELT(colnames, 3, mkChar("component_id"));

    // Wrap colnames in a list
    SEXP dimnames = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(dimnames, 0, R_NilValue); // No row names
    SET_VECTOR_ELT(dimnames, 1, colnames);

    // Assign the dimnames to the matrix
    setAttrib(MS_complex_info, R_DimNamesSymbol, dimnames);

    // -------------------------------------------------------------------------------------------------------------------
    //
    // Constructing the result list
    //
    // -------------------------------------------------------------------------------------------------------------------
    // Modify the return list to include cell sizes
    SEXP cell_sizes_R = PROTECT(allocVector(INTSXP, valid_cells.size())); nprot++;
    int* cell_sizes_ptr = INTEGER(cell_sizes_R);
    int i = 0;
    for (int cell : valid_cells) {
        cell_sizes_ptr[i++] = cell_sizes[cell];
    }

    int result_len = 9;
    SEXP result = PROTECT(allocVector(VECSXP, result_len)); nprot++;

    SET_VECTOR_ELT(result, 0, MS_complex_info);
    SET_VECTOR_ELT(result, 1, MS_graph_Rlist);
    SET_VECTOR_ELT(result, 2, MS_graph_vertex_label_Rlist);
    SET_VECTOR_ELT(result, 3, MS_graph_vertex_type_Rlist);
    SET_VECTOR_ELT(result, 4, MS_graph_edge_label_Rlist);
    SET_VECTOR_ELT(result, 5, MS_cell_graphs_Rlist);
    //SET_VECTOR_ELT(result, 6, trajectory_precell_vertices_Rlist);
    SET_VECTOR_ELT(result, 6, MS_cell_vertices_Rlist);
    SET_VECTOR_ELT(result, 7, trajectory_Rlist);
    SET_VECTOR_ELT(result, result_len - 1, cell_sizes_R);

    // Add names to list elements
    SEXP names = PROTECT(allocVector(STRSXP, result_len)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("MS_cx"));
    SET_STRING_ELT(names, 1, mkChar("MS_graph"));
    SET_STRING_ELT(names, 2, mkChar("MS_graph_labels"));
    SET_STRING_ELT(names, 3, mkChar("MS_graph_types"));
    SET_STRING_ELT(names, 4, mkChar("MS_graph_edge_label"));
    SET_STRING_ELT(names, 5, mkChar("MS_cell_graphs"));
    //SET_STRING_ELT(names, 6, mkChar("MS_precell_vertices"));
    SET_STRING_ELT(names, 6, mkChar("MS_cell_vertices"));
    SET_STRING_ELT(names, 7, mkChar("trajectories"));
    SET_STRING_ELT(names, result_len - 1, mkChar("cell_sizes"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result;
}
