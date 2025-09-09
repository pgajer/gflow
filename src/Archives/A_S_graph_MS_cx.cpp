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
