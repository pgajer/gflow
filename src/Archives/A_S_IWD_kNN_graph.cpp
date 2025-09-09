// Below you will find a definition of the S_IWD_kNN_graph function that computes and prunes an intersection-weighted k-nearest neighbors graph and returns an SEXP object (a named list) containing:
//    - adj_list: List of integer vectors containing adjacency lists (1-based indices)
//    - isize_list: List of integer vectors with intersection sizes for each edge
//    - dist_list: List of numeric vectors with distances for each edge
//    - conn_comps: Integer vector identifying connected components
//    - pruned_adj_list: Adjacency lists after edge pruning (1-based indices)
//    - pruned_isize_list: Intersection sizes for edges in pruned graph
//    - pruned_dist_list: Distances for edges in pruned graph
//    - long_edge_isize: Intersection sizes of removed edges
//    - alt_path_lengths: Lengths of alternative paths for removed edges
//    - alt_path_total_isize: Total intersection sizes along alternative paths

// The key parameter of S_IWD_kNN_graph whose optimal value we are seeking is an
// optimal number of nearest neighbors (k). To find this optimal k, I run
// S_IWD_kNN_graph for a range of k values and then use graph edit distance between
// consecutive pruned graphs, as in the following example, to find the global
// minimum of the resulting graph edit distances


// ## Finding IkNN graphs for different values of k
// k.values <- 5 * (16:24)
// for (k in k.values) {
//     cat("k:",k,"\n")
//     S.ibs.pca99.graph <- IWD.kNN.graph(S.ibs.pca99, k = k, method = "sp", n.cores = 1, max.alt.path.length = 10)
//     file <- paste0("~/current_projects/msr2_paper/data/S_ag_graphs/S_ibs_pca99_IWD_k",k,"_NN_graph.rda")
//     save(S.ibs.pca99.graph, file = file)
// }

// ## computing edit distances between consecutive pruned graphs
// ag.pruned.graphs.offset.1.GED <- c()
// offset <- 1
// offset.1.indices <- 1:(length(k.values) - offset)
// offset.1.k.values <- k.values[offset.1.indices]
// for (i in offset.1.indices) {
//     cat("\ri:",i)
//     k1 <- k.values[i]
//     k2 <- k.values[i + 1]
//     ag.pruned.graphs.offset.1.GED[i] <- graph.edit.distance(ag.pruned.graph.adj.list[[k1]],
//                                                             ag.pruned.graph.dist.list[[k1]],
//                                                             ag.pruned.graph.adj.list[[k2]],
//                                                             ag.pruned.graph.dist.list[[k2]])
// }

// If X is a large data matrix, the computation of all the components is expensive
// and unnecessary as I am only going to use adj_list and dist_list for the purpose
// of finding the optimal k. Also saving all these cmponents of the IkNN graph is
// time consuming and takes unnecessary space as for large values of k the size of
// the binary R file staring the output of graph.edit.distance may be more than
// 600MB.

// What do you think about adding a flag, only_prune (propose a better name if you
// think this one is not the best), to S_IWD_kNN_graph that would indicate that in
// the output list only adj_list and dist_list are non-empty and the other
// components are NULL.

// Also I want to add to the output list the following summary statistics components

// - The number of edges of the IkNN graph - n(E)
// - The number of edges of the pruned IkNN graph - n(pE)
// - The difference n(E) - n(pE)
// - The ratio (n(E) - n(pE)) / n(E)

// If only_prune is TRUE, then in the calculations of all components except is
// suppressed adj_list and dist_list, but the aboe summary statistics are always
// computed.

// Please generate a new version of S_IWD_kNN_graph implementing all these changes.


SEXP S_IWD_kNN_graph(SEXP s_X, SEXP s_k, SEXP s_max_alt_path_length) {
    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];
    int max_alt_path_length = INTEGER(s_max_alt_path_length)[0];

    // Creating a kNN graph
    std::unique_ptr<std::vector<std::vector<IkNN_vertex_t>>> IWD_kNN_graph_res = IWD_kNN_graph(s_X, s_k);
    if (!IWD_kNN_graph_res) Rf_error("IWD_kNN_graph function returned an invalid result (null pointer).");

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(n_vertices);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(n_vertices);
    auto dist_vect = std::make_unique<std::vector<std::vector<double>>>(n_vertices);
    int IWD_kNN_graph_n_vertices = IWD_kNN_graph_res->size();
    for (int i = 0; i < IWD_kNN_graph_n_vertices; i++) {
        for (auto nn_vertex : (*IWD_kNN_graph_res)[i]) {
            (*adj_vect)[i].push_back(nn_vertex.nn_index);
            (*intersection_size_vect)[i].push_back(nn_vertex.I_index);
            (*dist_vect)[i].push_back(nn_vertex.dist);
        }
    }

    // Preparing adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, (*adj_vect)[i].size()));
        int* A = INTEGER(RA);
        for (auto neighbor : (*adj_vect)[i])
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Preparing intersection size list
    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, (*intersection_size_vect)[i].size()));
        int* W = INTEGER(RW);
        for (auto size : (*intersection_size_vect)[i])
            *W++ = size;
        SET_VECTOR_ELT(intersection_size_list, i, RW);
        UNPROTECT(1);
    }

    // Preparing distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, (*dist_vect)[i].size()));
        double* D = REAL(RD);
        for (auto dist : (*dist_vect)[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
        UNPROTECT(1);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the IWD graph by deriving a IW kNN graph and then pruning that graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph = IWD_to_IW_kNN_graph(*IWD_kNN_graph_res);

    std::vector<int> long_edge_isize;
    std::vector<int> alt_path_lengths;
    std::vector<int> alt_path_total_isize;

    std::vector<std::vector<std::pair<int, int>>> pruned_graph = prune_edges_with_alt_paths(*IW_kNN_graph,
                                                                                            long_edge_isize,
                                                                                            alt_path_lengths,
                                                                                            alt_path_total_isize,
                                                                                            max_alt_path_length);

    SEXP long_edge_isize_Rvector = Cpp_vector_int_to_Rvector(long_edge_isize);
    SEXP alt_path_lengths_Rvector = Cpp_vector_int_to_Rvector(alt_path_lengths);
    SEXP alt_path_total_isize_Rvector = Cpp_vector_int_to_Rvector(alt_path_total_isize);

    // Creating pruned graph adjacency and isizes vector
    int pruned_graph_n_vertices = pruned_graph.size();
    if (pruned_graph_n_vertices != n_vertices) {
        Rf_error("pruned_graph_n_vertices: %d  n_vertices: %d\n - they should be the same!", pruned_graph_n_vertices, n_vertices);
    }
    auto pruned_adj_vect = std::vector<std::vector<int>>(n_vertices);
    auto pruned_isize_vect = std::vector<std::vector<int>>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        for (auto neighbor : pruned_graph[vertex]) {
            pruned_adj_vect[vertex].push_back(neighbor.first);
            pruned_isize_vect[vertex].push_back(neighbor.second);
        }
    }

    // Preparing pruned adjacency list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_adj_vect.at(i).size()));
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
        UNPROTECT(1);
    }

    // Preparing pruned isize list
    SEXP pruned_isize_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_isize_vect.at(vertex).size()));
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_isize_vect.at(vertex))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_isize_list, vertex, RA);
        UNPROTECT(1);
    }

    // Creating pruned_dist_list
    auto pruned_dist_vect = std::vector<std::vector<double>>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        for (auto nn_vertex : pruned_graph[vertex]) {
            int neighbor = nn_vertex.first;
            // Finding the corresponding distance in dist_vect
            auto it = std::find((*adj_vect)[vertex].begin(), (*adj_vect)[vertex].end(), neighbor);
            if (it != (*adj_vect)[vertex].end()) {
                int index = std::distance((*adj_vect)[vertex].begin(), it);
                pruned_dist_vect[vertex].push_back((*dist_vect)[vertex][index]);
            }
        }
    }

    // Preparing pruned distance list
    SEXP pruned_dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, pruned_dist_vect[i].size()));
        double* D = REAL(RD);
        for (auto dist : pruned_dist_vect[i])
            *D++ = dist;
        SET_VECTOR_ELT(pruned_dist_list, i, RD);
        UNPROTECT(1);
    }

    // Preparing the result list
    SEXP res = PROTECT(allocVector(VECSXP, 10)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, dist_list);
    SET_VECTOR_ELT(res, 3, conn_comps);
    SET_VECTOR_ELT(res, 4, pruned_adj_list);
    SET_VECTOR_ELT(res, 5, pruned_isize_list);
    SET_VECTOR_ELT(res, 6, pruned_dist_list);
    SET_VECTOR_ELT(res, 7, long_edge_isize_Rvector);
    SET_VECTOR_ELT(res, 8, alt_path_lengths_Rvector);
    SET_VECTOR_ELT(res, 9, alt_path_total_isize_Rvector);

    SEXP names = PROTECT(allocVector(STRSXP, 10)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("dist_list"));
    SET_STRING_ELT(names, 3, mkChar("conn_comps"));
    SET_STRING_ELT(names, 4, mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 5, mkChar("pruned_isize_list"));
    SET_STRING_ELT(names, 6, mkChar("pruned_dist_list"));
    SET_STRING_ELT(names, 7, mkChar("long_edge_isize"));
    SET_STRING_ELT(names, 8, mkChar("alt_path_lengths"));
    SET_STRING_ELT(names, 9, mkChar("alt_path_total_isize"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
