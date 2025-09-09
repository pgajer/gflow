//
// Testing 'optimized' code
//

// Here are results of the tests. Old method is about 36 to 39% faster

// k <- 20
// microbenchmark(old.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 0), new.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 1), times = 20L)

// ## n = 2000, k = 20
// ## old.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 0)
// ## new.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 1)
// ##      min       lq     mean   median       uq       max neval
// ## 6.946578 7.013080 7.139296 7.051944 7.105131  8.052643    20
// ## 9.543653 9.751969 9.866410 9.834663 9.932958 10.328668    20
// ## 9.834663 / 7.051944 = 1.394603

// ## n = 500, k = 50
// ##                                                                    expr
// ## old.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 0)
// ## new.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 1)
// ##      min       lq     mean   median       uq      max neval
// ## 1.735550 1.755074 1.770798 1.763900 1.792986 1.817628    20
// ## 2.316841 2.357167 2.373138 2.376306 2.390353 2.417653    20
// ## 2.376306 / 1.763900 = 1.347189


// ## n = 500, k = 20                                                                expr
// ## old.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 0)
// ## new.method <- IWD.kNN.graph(X, k = k, compute.full = TRUE, version = 1)
// ##      min       lq     mean   median       uq      max neval
// ## 503.3314 523.1087 535.6012 534.4764 544.8658 592.8008   100
// ## 685.7874 715.2077 726.4161 727.5064 736.3039 776.3735   100
// ## 727 / 534 = 1.361423


// performance improvemented code
std::unique_ptr<std::vector<std::vector<iknn_vertex_t>>> optimized_create_iknn_graph(SEXP RX, SEXP Rk) {
    PROTECT(RX = coerceVector(RX, REALSXP));
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    UNPROTECT(1);
    int n_points = dimX[0];
    PROTECT(Rk = coerceVector(Rk, INTSXP));
    int k = INTEGER(Rk)[0];
    UNPROTECT(1);

    SEXP knn_res = PROTECT(S_kNN(RX, Rk));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));
    UNPROTECT(1);

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);
    std::vector<int> intersection;
    intersection.reserve(k);

    auto res = std::make_unique<std::vector<std::vector<iknn_vertex_t>>>(n_points);
    for (auto& vec : *res) {
        vec.reserve(k);
    }

    int n_points_minus_one = n_points - 1;
    for (int pt_i = 0; pt_i < n_points_minus_one; pt_i++) {
        // Copy indices for pt_i
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + n_points * j];
            sorted_nn_i[j] = nn_i[j];
        }
        std::sort(sorted_nn_i.begin(), sorted_nn_i.end());

        for (int pt_j = pt_i + 1; pt_j < n_points; pt_j++) {
            // Copy indices for pt_j
            for (int j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + n_points * j];
                sorted_nn_j[j] = nn_j[j];
            }
            std::sort(sorted_nn_j.begin(), sorted_nn_j.end());

            intersection.clear();
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(),
                                sorted_nn_j.begin(), sorted_nn_j.end(),
                                std::back_inserter(intersection));

            int common_count = intersection.size();
            if (common_count > 0) {
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    // Use linear search for small k
                    auto it_i = std::find(nn_i.begin(), nn_i.end(), x_k);
                    auto it_j = std::find(nn_j.begin(), nn_j.end(), x_k);

                    int idx_i = it_i - nn_i.begin();
                    int idx_j = it_j - nn_j.begin();

                    double dist_i_k = distances[pt_i + n_points * idx_i];
                    double dist_j_k = distances[pt_j + n_points * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                (*res)[pt_i].emplace_back(iknn_vertex_t{pt_j, common_count, min_dist});
                (*res)[pt_j].emplace_back(iknn_vertex_t{pt_i, common_count, min_dist});
            }
        }
    }
    return res;
}

SEXP S_create_iknn_graph_with_version(SEXP s_X, SEXP s_k, SEXP s_max_alt_path_length, SEXP s_compute_full, SEXP s_version) {
    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];
    int max_alt_path_length = INTEGER(s_max_alt_path_length)[0];
    int compute_full = LOGICAL(s_compute_full)[0];

    int version = INTEGER(s_version)[0];

    // Creating a kNN graph
    //std::unique_ptr<std::vector<std::vector<iknn_vertex_t>>> iknn_graph = create_iknn_graph(s_X, s_k, static_cast<bool>(compute_full));
    std::unique_ptr<std::vector<std::vector<iknn_vertex_t>>> iknn_graph;
    if (version == 0) {
        Rprintf("Running OLD IkNN\n");
        iknn_graph = create_iknn_graph(s_X, s_k);
    } else {
        Rprintf("Running OPTIMIZED IkNN\n");
        iknn_graph = optimized_create_iknn_graph(s_X, s_k);
    }


    if (!iknn_graph) Rf_error("create_iknn_graph function returned an invalid result (null pointer).");

    // Count total edges in original graph
    int total_edges = 0;
    for (const auto& vertex_edges : *iknn_graph) {
        total_edges += vertex_edges.size();
    }

    // Create basic graph vectors needed for pruning
    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(n_vertices);
    auto dist_vect = std::make_unique<std::vector<std::vector<double>>>(n_vertices);
    int iknn_graph_n_vertices = iknn_graph->size();
    for (int i = 0; i < iknn_graph_n_vertices; i++) {
        for (auto nn_vertex : (*iknn_graph)[i]) {
            (*adj_vect)[i].push_back(nn_vertex.index);
            (*dist_vect)[i].push_back(nn_vertex.dist);
        }
    }

    // Initialize all components as NULL
    SEXP adj_list                     = R_NilValue;
    SEXP intersection_size_list       = R_NilValue;
    SEXP dist_list                    = R_NilValue;
    SEXP conn_comps                   = R_NilValue;
    SEXP pruned_adj_list             = R_NilValue;
    SEXP pruned_dist_list            = R_NilValue;
    SEXP pruned_isize_list           = R_NilValue;
    SEXP long_edge_isize_Rvector     = R_NilValue;
    SEXP alt_path_lengths_Rvector    = R_NilValue;
    SEXP alt_path_total_isize_Rvector = R_NilValue;

    // Compute pruned graph (needed for both cases)
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph = IWD_to_IW_kNN_graph(*iknn_graph);
    std::vector<int> long_edge_isize;
    std::vector<int> alt_path_lengths;
    std::vector<int> alt_path_total_isize;
    std::vector<std::vector<std::pair<int, int>>> pruned_graph = prune_edges_with_alt_paths(*IW_kNN_graph,
                                                                                            long_edge_isize,
                                                                                            alt_path_lengths,
                                                                                            alt_path_total_isize,
                                                                                            max_alt_path_length);

    // Count edges in pruned graph
    int pruned_edges = 0;
    for (const auto& vertex_edges : pruned_graph) {
        pruned_edges += vertex_edges.size();
    }

    // Compute graph statistics
    SEXP stats = PROTECT(allocVector(REALSXP, 4)); nprot++;
    REAL(stats)[0] = total_edges;
    REAL(stats)[1] = pruned_edges;
    REAL(stats)[2] = total_edges - pruned_edges;
    REAL(stats)[3] = (double)(total_edges - pruned_edges) / total_edges;

    // Always create pruned graph adjacency and distance lists
    pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    pruned_dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    // Create and fill pruned graph vectors
    auto pruned_adj_vect = std::vector<std::vector<int>>(n_vertices);
    auto pruned_dist_vect = std::vector<std::vector<double>>(n_vertices);

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        for (auto neighbor_pair : pruned_graph[vertex]) {
            int neighbor = neighbor_pair.first;
            pruned_adj_vect[vertex].push_back(neighbor);

            // Find corresponding distance
            auto it = std::find((*adj_vect)[vertex].begin(), (*adj_vect)[vertex].end(), neighbor);
            if (it != (*adj_vect)[vertex].end()) {
                int index = std::distance((*adj_vect)[vertex].begin(), it);
                pruned_dist_vect[vertex].push_back((*dist_vect)[vertex][index]);
            }
        }
    }

    // Fill pruned lists
    for (int i = 0; i < n_vertices; i++) {
        {
            SEXP RA = PROTECT(allocVector(INTSXP, pruned_adj_vect[i].size()));
            int* A = INTEGER(RA);
            for (auto neighbor : pruned_adj_vect[i])
                *A++ = neighbor + 1;
            SET_VECTOR_ELT(pruned_adj_list, i, RA);
            UNPROTECT(1);
        }

        {
            SEXP RD = PROTECT(allocVector(REALSXP, pruned_dist_vect[i].size()));
            double* D = REAL(RD);
            for (auto dist : pruned_dist_vect[i])
                *D++ = dist;
            SET_VECTOR_ELT(pruned_dist_list, i, RD);
            UNPROTECT(1);
        }
    }

    if (compute_full) {
        // Create full components
        adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        intersection_size_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        pruned_isize_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

        // Fill original graph components
        for (int i = 0; i < n_vertices; i++) {
            // Original graph adjacency list
            {
                SEXP RA = PROTECT(allocVector(INTSXP, (*adj_vect)[i].size()));
                int* A = INTEGER(RA);
                for (auto neighbor : (*adj_vect)[i])
                    *A++ = neighbor + 1;
                SET_VECTOR_ELT(adj_list, i, RA);
                UNPROTECT(1);
            }

            // Original graph intersection sizes
            {
                SEXP RW = PROTECT(allocVector(INTSXP, (*adj_vect)[i].size()));
                int* W = INTEGER(RW);
                for (auto nn_vertex : (*iknn_graph)[i])
                    *W++ = nn_vertex.isize;
                SET_VECTOR_ELT(intersection_size_list, i, RW);
                UNPROTECT(1);
            }

            // Original graph distances
            {
                SEXP RD = PROTECT(allocVector(REALSXP, (*dist_vect)[i].size()));
                double* D = REAL(RD);
                for (auto dist : (*dist_vect)[i])
                    *D++ = dist;
                SET_VECTOR_ELT(dist_list, i, RD);
                UNPROTECT(1);
            }

            // Pruned graph intersection sizes
            {
                SEXP RI = PROTECT(allocVector(INTSXP, pruned_adj_vect[i].size()));
                int* I = INTEGER(RI);
                for (auto neighbor_pair : pruned_graph[i])
                    *I++ = neighbor_pair.second;
                SET_VECTOR_ELT(pruned_isize_list, i, RI);
                UNPROTECT(1);
            }
        }

        // Compute connected components
        std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
        conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
        std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

        // Convert vectors to R vectors with proper protection
        long_edge_isize_Rvector = PROTECT(convert_vector_int_to_R(long_edge_isize)); nprot++;
        alt_path_lengths_Rvector = PROTECT(convert_vector_int_to_R(alt_path_lengths)); nprot++;
        alt_path_total_isize_Rvector = PROTECT(convert_vector_int_to_R(alt_path_total_isize)); nprot++;
    }

    // Prepare result list
    SEXP res = PROTECT(allocVector(VECSXP, 14)); nprot++;
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
    SET_VECTOR_ELT(res, 10, ScalarReal(REAL(stats)[0]));
    SET_VECTOR_ELT(res, 11, ScalarReal(REAL(stats)[1]));
    SET_VECTOR_ELT(res, 12, ScalarReal(REAL(stats)[2]));
    SET_VECTOR_ELT(res, 13, ScalarReal(REAL(stats)[3]));

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, 14)); nprot++;
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
    SET_STRING_ELT(names, 10, mkChar("n_edges"));
    SET_STRING_ELT(names, 11, mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(names, 12, mkChar("n_removed_edges"));
    SET_STRING_ELT(names, 13, mkChar("edge_reduction_ratio"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}

#if 0
SEXP S_iknn_graphs(
    SEXP s_X,
    SEXP s_kmin,
    SEXP s_kmax,
    SEXP s_max_alt_path_length,
    SEXP s_compute_full) {

    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;
    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];

    int kmin = INTEGER(s_kmin)[0];
    int kmax = INTEGER(s_kmax)[0];
    int max_alt_path_length = INTEGER(s_max_alt_path_length)[0];
    int compute_full = LOGICAL(s_compute_full)[0];

    std::unordered_map<edge_t, birth_death_time_t, edge_hash> birth_death_map;
    std::vector<std::vector<double>> k_statistics;
    std::vector<SEXP> r_pruned_graphs;

    // Process each k value
    for (int k = kmin; k <= kmax; k++) {
        // Get the IWD-kNN graph directly
        auto graph_result = IWD_kNN_graph(s_X, ScalarInteger(k));
        if (!graph_result) Rf_error("IWD_kNN_graph returned null pointer for k=%d", k);

        // Count edges in original graph
        int total_edges = 0;
        for (const auto& vertex_edges : *graph_result) {
            total_edges += vertex_edges.size();
        }
        total_edges /= 2;

        auto adj_list = std::vector<std::vector<int>>>(n_vertices);
        auto dist_list = std::vector<std::vector<double>>(n_vertices);
        for (int i = 0; i < n_vertices; i++) {
            for (auto neighbor : (*graph_result)[i]) {
                adj_list[i].push_back(neighbor.index);
                dist_list[i].push_back(neighbor.dist);
            }
        }

        // Process and prune the graph
        std::vector<int> long_edge_isize;
        std::vector<int> alt_path_lengths;
        std::vector<int> alt_path_total_isize;

        auto pruned_graph = prune_graph(*graph_result,
                                        max_alt_path_length,
                                        long_edge_isize,
                                        alt_path_lengths,
                                        alt_path_total_isize);

        // Count pruned edges and update statistics
        int pruned_edges = 0;
        for (const auto& vertex_edges : pruned_graph) {
            pruned_edges += vertex_edges.size();
        }
        pruned_edges /= 2;

        // Create and fill pruned graph vectors
        auto pruned_adj_list = std::vector<std::vector<int>>(n_vertices);
        auto pruned_dist_list = std::vector<std::vector<double>>(n_vertices);
        for (int vertex = 0; vertex < n_vertices; vertex++) {
            for (auto neighbor_pair : pruned_graph[vertex]) {
                int neighbor = neighbor_pair.first;
                pruned_adj_list[vertex].push_back(neighbor);

                // Find corresponding distance
                auto it = std::find(adj_list[vertex].begin(), adj_list[vertex].end(), neighbor);
                if (it != adj_list[vertex].end()) {
                    int index = std::distance(adj_list[vertex].begin(), it);
                    pruned_dist_list[vertex].push_back(dist_list[vertex][index]);
                }
            }
        }

        // Fill pruned lists
        r_pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices));
        r_pruned_dist_list = PROTECT(allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            {
                SEXP RA = PROTECT(allocVector(INTSXP, pruned_adj_vect[i].size()));
                int* A = INTEGER(RA);
                for (auto neighbor : pruned_adj_vect[i])
                    *A++ = neighbor + 1;
                SET_VECTOR_ELT(r_pruned_adj_list, i, RA);
                UNPROTECT(1);
            }

            {
                SEXP RD = PROTECT(allocVector(REALSXP, pruned_dist_vect[i].size()));
                double* D = REAL(RD);
                for (auto dist : pruned_dist_vect[i])
                    *D++ = dist;
                SET_VECTOR_ELT(r_pruned_dist_list, i, RD);
                UNPROTECT(1);
            }
        }
        UNPROTECT(2); // for r_pruned_adj_list and r_pruned_dist_list


        k_statistics.push_back({
            (double)total_edges,
            (double)pruned_edges,
            (double)(total_edges - pruned_edges),
            (double)(total_edges - pruned_edges) / total_edges
        });

        // Update birth/death times
        for (size_t i = 0; i < pruned_graph.size(); i++) {
            for (const auto& edge_data : pruned_graph[i]) {
                int neighbor = edge_data.first;
                if (i < (size_t)neighbor) {
                    edge_t edge = {i, (size_t)neighbor};

                    auto it = birth_death_map.find(edge);
                    if (it == birth_death_map.end()) {
                        birth_death_map[edge] = {(size_t)k, (size_t)(k + 1)};
                    } else if (it->second.death_time == SIZE_MAX) {
                        it->second.death_time = (size_t)(k + 1);
                    }
                }
            }
        }

        // Create R representation of pruned graph if needed
        if (compute_full) {
            SEXP pruned_graph_R = PROTECT(create_R_graph_representation(pruned_graph, *graph_result));
            r_pruned_graphs.push_back(pruned_graph_R);
            nprot++;
        }
    }

    // Create statistics matrix
    SEXP k_stats_matrix = PROTECT(allocMatrix(REALSXP, k_statistics.size(), 5)); nprot++;
    double* data = REAL(k_stats_matrix);

    for (size_t i = 0; i < k_statistics.size(); i++) {
        data[i] = kmin + i;
        for (size_t j = 0; j < 4; j++) {
            data[i + (j + 1) * k_statistics.size()] = k_statistics[i][j];
        }
    }

    // Set column names for statistics matrix
    SEXP stat_colnames = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(stat_colnames, 0, mkChar("k"));
    SET_STRING_ELT(stat_colnames, 1, mkChar("n_edges"));
    SET_STRING_ELT(stat_colnames, 2, mkChar("n_pruned_edges"));
    SET_STRING_ELT(stat_colnames, 3, mkChar("n_removed_edges"));
    SET_STRING_ELT(stat_colnames, 4, mkChar("edge_reduction_ratio"));
    setAttrib(k_stats_matrix, R_NamesSymbol, stat_colnames);

    // Create birth/death matrix
    SEXP birth_death_matrix = PROTECT(convert_birth_death_map_to_matrix(birth_death_map)); nprot++;

    // Create pruned graphs list if compute_full is true
    SEXP pruned_graphs_list = R_NilValue;
    if (compute_full) {
        PROTECT(pruned_graphs_list = allocVector(VECSXP, r_pruned_graphs.size())); nprot++;
        for (size_t i = 0; i < r_pruned_graphs.size(); i++) {
            SET_VECTOR_ELT(pruned_graphs_list, i, r_pruned_graphs[i]);
        }
    }

    // Create return list
    SEXP result = PROTECT(allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(result, 0, birth_death_matrix);
    SET_VECTOR_ELT(result, 1, k_stats_matrix);
    SET_VECTOR_ELT(result, 2, pruned_graphs_list);

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("birth_death_matrix"));
    SET_STRING_ELT(names, 1, mkChar("k_statistics"));
    SET_STRING_ELT(names, 2, mkChar("pruned_graphs"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return result;
}
#endif




SEXP S_create_iknn_graph(SEXP s_X, SEXP s_k, SEXP s_max_alt_path_length, SEXP s_compute_full) {
    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];
    int max_alt_path_length = INTEGER(s_max_alt_path_length)[0];
    int compute_full = LOGICAL(s_compute_full)[0];

    // Creating a kNN graph
    auto iknn_graph = create_iknn_graph(s_X, s_k);

    // Count total edges in original graph
    int total_edges = 0;
    for (const auto& vertex_edges : iknn_graph.graph) {
        total_edges += vertex_edges.size();
    }
    total_edges /= 2;

    //Rprintf("total_edges: %d\n",total_edges);

    // Create basic graph vectors needed for pruning
    auto adj_vect = std::vector<std::vector<int>>(n_vertices);
    auto dist_vect = std::vector<std::vector<double>>(n_vertices);
    for (size_t i = 0; i < iknn_graph.graph.size(); i++) {
        for (auto nn_vertex : iknn_graph.graph[i]) {
            adj_vect[i].push_back(nn_vertex.index);
            dist_vect[i].push_back(nn_vertex.dist);
        }
    }

    // Initialize all components as NULL
    SEXP adj_list                     = R_NilValue;
    SEXP intersection_size_list       = R_NilValue;
    SEXP dist_list                    = R_NilValue;
    SEXP s_conn_comps                 = R_NilValue;
    SEXP pruned_adj_list              = R_NilValue;
    SEXP pruned_dist_list             = R_NilValue;
    SEXP pruned_isize_list            = R_NilValue;
    SEXP long_edge_isize_Rvector      = R_NilValue;
    SEXP alt_path_lengths_Rvector     = R_NilValue;
    SEXP alt_path_total_isize_Rvector = R_NilValue;

    // Compute pruned graph (needed for both cases)
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph = IWD_to_IW_kNN_graph(iknn_graph.graph);
    std::vector<int> long_edge_isize;
    std::vector<int> alt_path_lengths;
    std::vector<int> alt_path_total_isize;
    std::vector<std::vector<std::pair<int, int>>> pruned_graph = prune_edges_with_alt_paths(*IW_kNN_graph,
                                                                                            long_edge_isize,
                                                                                            alt_path_lengths,
                                                                                            alt_path_total_isize,
                                                                                            max_alt_path_length);

    // Count edges in pruned graph
    int pruned_edges = 0;
    for (const auto& vertex_edges : pruned_graph) {
        pruned_edges += vertex_edges.size();
    }
    pruned_edges /= 2;

    //Rprintf("pruned_edges: %d\n", pruned_edges);


    // Compute graph statistics
    SEXP stats = PROTECT(allocVector(REALSXP, 4)); nprot++;
    REAL(stats)[0] = total_edges;
    REAL(stats)[1] = pruned_edges;
    REAL(stats)[2] = total_edges - pruned_edges;
    REAL(stats)[3] = (double)(total_edges - pruned_edges) / total_edges;

    // Always create pruned graph adjacency and distance lists
    pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    pruned_dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    // Create and fill pruned graph vectors
    auto pruned_adj_vect = std::vector<std::vector<int>>(n_vertices);
    auto pruned_dist_vect = std::vector<std::vector<double>>(n_vertices);

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        for (auto neighbor_pair : pruned_graph[vertex]) {
            int neighbor = neighbor_pair.first;
            pruned_adj_vect[vertex].push_back(neighbor);

            // Find corresponding distance
            auto it = std::find(adj_vect[vertex].begin(), adj_vect[vertex].end(), neighbor);
            if (it != adj_vect[vertex].end()) {
                int index = std::distance(adj_vect[vertex].begin(), it);
                pruned_dist_vect[vertex].push_back(dist_vect[vertex][index]);
            }
        }
    }

    // Fill pruned lists
    for (int i = 0; i < n_vertices; i++) {
        {
            SEXP RA = PROTECT(allocVector(INTSXP, pruned_adj_vect[i].size()));
            int* A = INTEGER(RA);
            for (auto neighbor : pruned_adj_vect[i])
                *A++ = neighbor + 1;
            SET_VECTOR_ELT(pruned_adj_list, i, RA);
            UNPROTECT(1);
        }

        {
            SEXP RD = PROTECT(allocVector(REALSXP, pruned_dist_vect[i].size()));
            double* D = REAL(RD);
            for (auto dist : pruned_dist_vect[i])
                *D++ = dist;
            SET_VECTOR_ELT(pruned_dist_list, i, RD);
            UNPROTECT(1);
        }
    }

    if (compute_full) {
        // Create full components
        adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        intersection_size_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
        pruned_isize_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

        // Fill original graph components
        for (int i = 0; i < n_vertices; i++) {
            // Original graph adjacency list
            {
                SEXP RA = PROTECT(allocVector(INTSXP, adj_vect[i].size()));
                int* A = INTEGER(RA);
                for (auto neighbor : adj_vect[i])
                    *A++ = neighbor + 1;
                SET_VECTOR_ELT(adj_list, i, RA);
                UNPROTECT(1);
            }

            // Original graph intersection sizes
            {
                SEXP RW = PROTECT(allocVector(INTSXP, adj_vect[i].size()));
                int* W = INTEGER(RW);
                for (auto nn_vertex : iknn_graph.graph[i])
                    *W++ = nn_vertex.isize;
                SET_VECTOR_ELT(intersection_size_list, i, RW);
                UNPROTECT(1);
            }

            // Original graph distances
            {
                SEXP RD = PROTECT(allocVector(REALSXP, dist_vect[i].size()));
                double* D = REAL(RD);
                for (auto dist : dist_vect[i])
                    *D++ = dist;
                SET_VECTOR_ELT(dist_list, i, RD);
                UNPROTECT(1);
            }

            // Pruned graph intersection sizes
            {
                SEXP RI = PROTECT(allocVector(INTSXP, pruned_adj_vect[i].size()));
                int* I = INTEGER(RI);
                for (auto neighbor_pair : pruned_graph[i])
                    *I++ = neighbor_pair.second;
                SET_VECTOR_ELT(pruned_isize_list, i, RI);
                UNPROTECT(1);
            }
        }

        // Create a unique_ptr wrapping your vector
        // auto adj_vect_ptr = std::make_unique<std::vector<std::vector<int>>>(std::move(adj_vect)); // the original adj_vect is empty now as std::move transfered its ownership

        // Compute connected components
        std::vector<int> conn_comps = union_find(adj_vect);

        s_conn_comps = PROTECT(allocVector(INTSXP, conn_comps.size())); nprot++;
        std::copy(conn_comps.begin(), conn_comps.end(), INTEGER(s_conn_comps));

        // Convert vectors to R vectors with proper protection
        long_edge_isize_Rvector = PROTECT(convert_vector_int_to_R(long_edge_isize)); nprot++;
        alt_path_lengths_Rvector = PROTECT(convert_vector_int_to_R(alt_path_lengths)); nprot++;
        alt_path_total_isize_Rvector = PROTECT(convert_vector_int_to_R(alt_path_total_isize)); nprot++;
    }

    // Prepare result list
    SEXP res = PROTECT(allocVector(VECSXP, 14)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, dist_list);
    SET_VECTOR_ELT(res, 3, s_conn_comps);
    SET_VECTOR_ELT(res, 4, pruned_adj_list);
    SET_VECTOR_ELT(res, 5, pruned_isize_list);
    SET_VECTOR_ELT(res, 6, pruned_dist_list);
    SET_VECTOR_ELT(res, 7, long_edge_isize_Rvector);
    SET_VECTOR_ELT(res, 8, alt_path_lengths_Rvector);
    SET_VECTOR_ELT(res, 9, alt_path_total_isize_Rvector);
    SET_VECTOR_ELT(res, 10, ScalarReal(REAL(stats)[0]));
    SET_VECTOR_ELT(res, 11, ScalarReal(REAL(stats)[1]));
    SET_VECTOR_ELT(res, 12, ScalarReal(REAL(stats)[2]));
    SET_VECTOR_ELT(res, 13, ScalarReal(REAL(stats)[3]));

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, 14)); nprot++;
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
    SET_STRING_ELT(names, 10, mkChar("n_edges"));
    SET_STRING_ELT(names, 11, mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(names, 12, mkChar("n_removed_edges"));
    SET_STRING_ELT(names, 13, mkChar("edge_reduction_ratio"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}
