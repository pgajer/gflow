SEXP original_S_create_iknn_graphs(
    SEXP s_X,
    SEXP s_kmin,
    SEXP s_kmax,
    SEXP s_pruning_thld,
    SEXP s_long_edge_thld,
    SEXP s_compute_full,
    SEXP s_verbose) {

    auto total_start_time = std::chrono::steady_clock::now();

    double pruning_thld = REAL(s_pruning_thld)[0];
    int verbose = LOGICAL(s_verbose)[0];

    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;

    int kmin = INTEGER(s_kmin)[0];
    int kmax = INTEGER(s_kmax)[0];
    int compute_full = LOGICAL(s_compute_full)[0];
    double long_edge_thld = REAL(s_long_edge_thld)[0];

    int* dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];

    if (verbose) {
        Rprintf("Processing k values from %d to %d for %d vertices\n", kmin, kmax, n_vertices);
    }

    // Prepare vectors to store results
    std::vector<std::vector<double>> k_statistics(kmax - kmin + 1);

    // Create vector of k values
    std::vector<int> k_values(kmax - kmin + 1);
    std::iota(k_values.begin(), k_values.end(), kmin);

    // Parallel processing of graph creation and pruning
    if (verbose) Rprintf("Starting parallel graph processing\n");
    auto parallel_start_time = std::chrono::steady_clock::now();

    // Progress tracking
    std::atomic<int> progress_counter{0};
    const size_t n_k_values = kmax - kmin + 1;

    // Compute kNN once for maximum k
    auto knn_results = compute_knn(s_X, kmax);

    #define USE_GEOMETRIC_PRUNING 1

    #if USE_GEOMETRIC_PRUNING
    std::vector<set_wgraph_t> pruned_graphs(kmax - kmin + 1);
    // Add vector for double pruned graphs
    std::vector<set_wgraph_t> double_pruned_graphs(kmax - kmin + 1);
    #else
    int max_alt_path_length = 2; //INTEGER(s_max_alt_path_length)[0];
    std::vector<vect_wgraph_t> pruned_graphs(kmax - kmin + 1);
    #endif

    // Parallel processing using pre-computed kNN results
    std::for_each(std::execution::par_unseq,
                  k_values.begin(),
                  k_values.end(),
                  [&](int k) {

                      if (verbose) {
                          int current_count = ++progress_counter;
                          REprintf("\rProcessing %d %d%%",
                                   current_count,
                                   static_cast<int>((100.0 * current_count) / n_k_values));
                      }

                      auto iknn_graph = create_iknn_graph(knn_results, k);

                      // Count original edges
                      size_t n_edges = 0;
                      for (const auto& vertex_edges : iknn_graph.graph) {
                          n_edges += vertex_edges.size();
                      }
                      n_edges /= 2;


                      #if USE_GEOMETRIC_PRUNING
                      // geometric graph pruning
                      // transfering iknn_graph_t to set_wgraph_t
                      auto pruned_graph = set_wgraph_t(iknn_graph);

                      // Compute the deviations using the optimized method
                      auto rel_deviations = pruned_graph.compute_edge_weight_rel_deviations();

                      // Process edges
                      for (size_t i = 0; i < rel_deviations.size(); i++) {
                          if (rel_deviations[i].rel_deviation < pruning_thld) {
                              size_t source = rel_deviations[i].source;
                              size_t target = rel_deviations[i].target;

                              // Check if removing this edge would isolate any vertex
                              if (pruned_graph.adjacency_list[source].size() <= 1) {
                                  continue;
                              }

                              if (pruned_graph.adjacency_list[target].size() <= 1) {
                                  continue;
                              }

                              // Safe to remove the edge
                              pruned_graph.remove_edge(source, target);
                          }
                      }

                      // Apply the set_wgraph_t::prune_long_edges method
                      // Create a double-pruned graph by applying prune_long_edges to the pruned graph
                      set_wgraph_t double_pruned_graph = pruned_graph.prune_long_edges(long_edge_thld);

                      // Count edges in the pruned graph
                      size_t n_edges_in_pruned_graph = 0;
                      for (const auto& neighbors : pruned_graph.adjacency_list) {
                          n_edges_in_pruned_graph += neighbors.size();
                      }
                      n_edges_in_pruned_graph /= 2;

                      // Count edges in the double-pruned graph
                      size_t n_edges_in_double_pruned_graph = 0;
                      for (const auto& neighbors : double_pruned_graph.adjacency_list) {
                          n_edges_in_double_pruned_graph += neighbors.size();
                      }
                      n_edges_in_double_pruned_graph /= 2;

                      k_statistics[k - kmin] = {
                          (double)n_edges,
                          (double)n_edges_in_pruned_graph,
                          (double)(n_edges - n_edges_in_pruned_graph),
                          (double)(n_edges - n_edges_in_pruned_graph) / n_edges,
                          (double)n_edges_in_double_pruned_graph,
                          (double)(n_edges_in_pruned_graph - n_edges_in_double_pruned_graph),
                          (double)(n_edges_in_pruned_graph - n_edges_in_double_pruned_graph) / n_edges_in_pruned_graph
                      };

                      // Store the double-pruned graph
                      double_pruned_graphs[k - kmin] = std::move(double_pruned_graph);

                      // Store the first-level pruned graph
                      pruned_graphs[k - kmin] = std::move(pruned_graph);

                      #else
                      // isize graph pruning
                      vect_wgraph_t pruned_graph = iknn_graph.prune_graph(max_alt_path_length);

                      // Count edges in the pruned graph
                      size_t n_edges_in_pruned_graph = 0;
                      for (const auto& neighbors : pruned_graph.adjacency_list) {
                          n_edges_in_pruned_graph += neighbors.size();
                      }
                      n_edges_in_pruned_graph /= 2;

                      // Store results
                      pruned_graphs[k - kmin] = std::move(pruned_graph);
                      k_statistics[k - kmin] = {
                          (double)n_edges,
                          (double)n_edges_in_pruned_graph,
                          (double)(n_edges - n_edges_in_pruned_graph),
                          (double)(n_edges - n_edges_in_pruned_graph) / n_edges
                      };
                      #endif
                  });

    if (verbose) {
        elapsed_time(parallel_start_time, "\rParallel processing completed", true);
        Rprintf("Processing edge birth/death time ... ");
    }

    // Serial processing of results
    auto serial_start_time = std::chrono::steady_clock::now();

    // Process birth/death times and edge frequencies
    std::unordered_map<edge_t, birth_death_time_t, edge_hash> birth_death_map;
    std::unordered_map<edge_t, size_t, edge_hash> edge_freq_map;

    // Add maps for double-pruned graphs
    std::unordered_map<edge_t, birth_death_time_t, edge_hash> double_birth_death_map;
    std::unordered_map<edge_t, size_t, edge_hash> double_edge_freq_map;

    for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
        int k = kmin + k_idx;
        const auto& pruned_graph = pruned_graphs[k_idx];

        // Process edges for first-level pruned graph
        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};

                    // Update birth/death times
                    auto it = birth_death_map.find(edge_key);
                    if (it == birth_death_map.end()) {
                        birth_death_map[edge_key] = {(size_t)k, (size_t)(k + 1)};
                    } else {
                        it->second.death_time = (size_t)(k + 1);
                    }

                    // Update edge frequency
                    edge_freq_map[edge_key]++;
                }
            }
        }

        // Process edges for double-pruned graph
        #if USE_GEOMETRIC_PRUNING
        const auto& double_pruned_graph = double_pruned_graphs[k_idx];

        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : double_pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};

                    // Update birth/death times for double-pruned graph
                    auto it = double_birth_death_map.find(edge_key);
                    if (it == double_birth_death_map.end()) {
                        double_birth_death_map[edge_key] = {(size_t)k, (size_t)(k + 1)};
                    } else {
                        it->second.death_time = (size_t)(k + 1);
                    }

                    // Update edge frequency for double-pruned graph
                    double_edge_freq_map[edge_key]++;
                }
            }
        }
        #endif
    }

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        serial_start_time = std::chrono::steady_clock::now();
        Rprintf("Processing edge frequencies ... ");
    }

    // Compute relative frequencies and k-specific statistics
    std::unordered_map<edge_t, double, edge_hash> edge_rel_freq_map;
    std::vector<double> k_freq_ratios(n_k_values);

    // Add maps for double-pruned graphs
    std::unordered_map<edge_t, double, edge_hash> double_edge_rel_freq_map;
    std::vector<double> double_k_freq_ratios(n_k_values);

    for (const auto& [edge, freq] : edge_freq_map) {
        edge_rel_freq_map[edge] = static_cast<double>(freq) / n_k_values;
    }

    // Process for double-pruned graphs
    for (const auto& [edge, freq] : double_edge_freq_map) {
        double_edge_rel_freq_map[edge] = static_cast<double>(freq) / n_k_values;
    }

    // Calculate frequency ratios for each k
    for (int k_idx = 0; k_idx < n_k_values; k_idx++) {
        // Process for first-level pruned graph
        const auto& pruned_graph = pruned_graphs[k_idx];
        double total_rel_freq = 0.0;
        size_t edge_count = 0;

        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};
                    total_rel_freq += edge_rel_freq_map[edge_key];
                    edge_count++;
                }
            }
        }

        k_freq_ratios[k_idx] = edge_count > 0 ? total_rel_freq / edge_count : 0.0;

        // Process for double-pruned graph
        #if USE_GEOMETRIC_PRUNING
        const auto& double_pruned_graph = double_pruned_graphs[k_idx];
        double double_total_rel_freq = 0.0;
        size_t double_edge_count = 0;

        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : double_pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};
                    double_total_rel_freq += double_edge_rel_freq_map[edge_key];
                    double_edge_count++;
                }
            }
        }

        double_k_freq_ratios[k_idx] = double_edge_count > 0 ? double_total_rel_freq / double_edge_count : 0.0;
        #endif
    }

    // Add frequency statistics to k_statistics
    for (size_t i = 0; i < k_statistics.size(); i++) {
        k_statistics[i].push_back(k_freq_ratios[i]);
        #if USE_GEOMETRIC_PRUNING
        k_statistics[i].push_back(double_k_freq_ratios[i]);
        #endif
    }

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        serial_start_time = std::chrono::steady_clock::now();
        Rprintf("Creating return list objects ... ");
    }

    // Create R objects
    SEXP pruned_graphs_list = R_NilValue;
    SEXP double_pruned_graphs_list = R_NilValue;

    if (compute_full) {
        PROTECT(pruned_graphs_list = allocVector(VECSXP, kmax - kmin + 1)); nprot++;
        #if USE_GEOMETRIC_PRUNING
        PROTECT(double_pruned_graphs_list = allocVector(VECSXP, kmax - kmin + 1)); nprot++;
        #endif

        for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
            // Process first-level pruned graph
            const auto& pruned_graph = pruned_graphs[k_idx];

            SEXP r_pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices));
            SEXP r_pruned_weight_list = PROTECT(allocVector(VECSXP, n_vertices));

            for (int i = 0; i < n_vertices; i++) {
                const auto& edges = pruned_graph.adjacency_list[i];

                {
                    SEXP RA = PROTECT(allocVector(INTSXP, edges.size()));
                    int* A = INTEGER(RA);
                    for (const auto& edge : edges) {
                        *A++ = edge.vertex + 1;
                    }
                    SET_VECTOR_ELT(r_pruned_adj_list, i, RA);
                    UNPROTECT(1);
                }

                {
                    SEXP RD = PROTECT(allocVector(REALSXP, edges.size()));
                    double* D = REAL(RD);
                    for (const auto& edge : edges) {
                        *D++ = edge.weight;
                    }
                    SET_VECTOR_ELT(r_pruned_weight_list, i, RD);
                    UNPROTECT(1);
                }
            }

            SEXP r_pruned_graph = PROTECT(allocVector(VECSXP, 2));
            SET_VECTOR_ELT(r_pruned_graph, 0, r_pruned_adj_list);
            SET_VECTOR_ELT(r_pruned_graph, 1, r_pruned_weight_list);

            SEXP r_pruned_graph_names = PROTECT(allocVector(STRSXP, 2));
            SET_STRING_ELT(r_pruned_graph_names, 0, mkChar("pruned_adj_list"));
            SET_STRING_ELT(r_pruned_graph_names, 1, mkChar("pruned_weight_list"));
            setAttrib(r_pruned_graph, R_NamesSymbol, r_pruned_graph_names);

            SET_VECTOR_ELT(pruned_graphs_list, k_idx, r_pruned_graph);
            UNPROTECT(4);

            // Process double-pruned graph
            #if USE_GEOMETRIC_PRUNING
            const auto& double_pruned_graph = double_pruned_graphs[k_idx];

            SEXP r_double_pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices));
            SEXP r_double_pruned_weight_list = PROTECT(allocVector(VECSXP, n_vertices));

            for (int i = 0; i < n_vertices; i++) {
                const auto& edges = double_pruned_graph.adjacency_list[i];

                {
                    SEXP RA = PROTECT(allocVector(INTSXP, edges.size()));
                    int* A = INTEGER(RA);
                    for (const auto& edge : edges) {
                        *A++ = edge.vertex + 1;
                    }
                    SET_VECTOR_ELT(r_double_pruned_adj_list, i, RA);
                    UNPROTECT(1);
                }

                {
                    SEXP RD = PROTECT(allocVector(REALSXP, edges.size()));
                    double* D = REAL(RD);
                    for (const auto& edge : edges) {
                        *D++ = edge.weight;
                    }
                    SET_VECTOR_ELT(r_double_pruned_weight_list, i, RD);
                    UNPROTECT(1);
                }
            }

            SEXP r_double_pruned_graph = PROTECT(allocVector(VECSXP, 2));
            SET_VECTOR_ELT(r_double_pruned_graph, 0, r_double_pruned_adj_list);
            SET_VECTOR_ELT(r_double_pruned_graph, 1, r_double_pruned_weight_list);

            SEXP r_double_pruned_graph_names = PROTECT(allocVector(STRSXP, 2));
            SET_STRING_ELT(r_double_pruned_graph_names, 0, mkChar("pruned_adj_list"));
            SET_STRING_ELT(r_double_pruned_graph_names, 1, mkChar("pruned_weight_list"));
            setAttrib(r_double_pruned_graph, R_NamesSymbol, r_double_pruned_graph_names);

            SET_VECTOR_ELT(double_pruned_graphs_list, k_idx, r_double_pruned_graph);
            UNPROTECT(4);
            #endif
        }
    }

    // Create statistics matrix
    #if USE_GEOMETRIC_PRUNING
    // Now 9 columns instead of 6
    SEXP k_stats_matrix = PROTECT(allocMatrix(REALSXP, k_statistics.size(), 9)); nprot++;
    #else
    SEXP k_stats_matrix = PROTECT(allocMatrix(REALSXP, k_statistics.size(), 6)); nprot++;
    #endif

    double* data = REAL(k_stats_matrix);

    for (size_t i = 0; i < k_statistics.size(); i++) {
        data[i] = kmin + i;
        for (size_t j = 0; j < k_statistics[i].size(); j++) {
            data[i + (j + 1) * k_statistics.size()] = k_statistics[i][j];
        }
    }

    // Set column names
    SEXP k_stats_dimnames = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(k_stats_dimnames, 0, R_NilValue);

    #if USE_GEOMETRIC_PRUNING
    SEXP k_stats_colnames = PROTECT(allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(k_stats_colnames, 0, mkChar("k"));
    SET_STRING_ELT(k_stats_colnames, 1, mkChar("n_edges"));
    SET_STRING_ELT(k_stats_colnames, 2, mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 3, mkChar("n_removed_edges"));
    SET_STRING_ELT(k_stats_colnames, 4, mkChar("edge_reduction_ratio"));
    SET_STRING_ELT(k_stats_colnames, 5, mkChar("n_edges_in_double_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 6, mkChar("n_removed_edges_in_double_pruning"));
    SET_STRING_ELT(k_stats_colnames, 7, mkChar("double_edge_reduction_ratio"));
    SET_STRING_ELT(k_stats_colnames, 8, mkChar("frequency_ratio"));
    #else
    SEXP k_stats_colnames = PROTECT(allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(k_stats_colnames, 0, mkChar("k"));
    SET_STRING_ELT(k_stats_colnames, 1, mkChar("n_edges"));
    SET_STRING_ELT(k_stats_colnames, 2, mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 3, mkChar("n_removed_edges"));
    SET_STRING_ELT(k_stats_colnames, 4, mkChar("edge_reduction_ratio"));
    SET_STRING_ELT(k_stats_colnames, 5, mkChar("frequency_ratio"));
    #endif

    SET_VECTOR_ELT(k_stats_dimnames, 1, k_stats_colnames);
    setAttrib(k_stats_matrix, R_DimNamesSymbol, k_stats_dimnames);

    // Create birth/death matrix
    SEXP birth_death_matrix = PROTECT(convert_birth_death_map_to_matrix(birth_death_map)); nprot++;

    // Create double-pruned birth/death matrix
    #if USE_GEOMETRIC_PRUNING
    SEXP double_birth_death_matrix = PROTECT(convert_birth_death_map_to_matrix(double_birth_death_map)); nprot++;
    #endif

    // Create return list
    #if USE_GEOMETRIC_PRUNING
    SEXP result = PROTECT(allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(result, 0, birth_death_matrix);
    SET_VECTOR_ELT(result, 1, k_stats_matrix);
    SET_VECTOR_ELT(result, 2, pruned_graphs_list);
    SET_VECTOR_ELT(result, 3, double_birth_death_matrix);
    SET_VECTOR_ELT(result, 4, double_pruned_graphs_list);

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("birth_death_matrix"));
    SET_STRING_ELT(names, 1, mkChar("k_statistics"));
    SET_STRING_ELT(names, 2, mkChar("pruned_graphs"));
    SET_STRING_ELT(names, 3, mkChar("double_birth_death_matrix"));
    SET_STRING_ELT(names, 4, mkChar("double_pruned_graphs"));
    #else
    SEXP result = PROTECT(allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(result, 0, birth_death_matrix);
    SET_VECTOR_ELT(result, 1, k_stats_matrix);
    SET_VECTOR_ELT(result, 2, pruned_graphs_list);

    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("birth_death_matrix"));
    SET_STRING_ELT(names, 1, mkChar("k_statistics"));
    SET_STRING_ELT(names, 2, mkChar("pruned_graphs"));
    #endif

    setAttrib(result, R_NamesSymbol, names);

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        elapsed_time(total_start_time, "Total elapsed time", true);
    }

    UNPROTECT(nprot);
    return result;
}

/**
 * @brief Computes and analyzes a sequence of pruned intersection-weighted k-nearest neighbors graphs
 *
 * @details
 * For a given range of k values [kmin, kmax], this function:
 * 1. Constructs and prunes intersection-weighted k-nearest neighbors (IWD-kNN) graphs
 *    using the iknn_graph.prune_graph() method
 * 2. Tracks the birth and death times of edges across different k values
 * 3. Computes statistics about edge persistence, pruning efficiency, and graph stability
 * 4. Calculates frequency metrics to identify the most stable pruned graph configurations
 *
 * For each edge e, birth time b(e) and death time d(e) are defined as:
 * - b(e) = min{k : e ∈ pG(X,k)}
 * - d(e) = min{k : e ∈ pG(X,k-1) and e ∉ pG(X,s) for all s ≥ k}
 * where pG(X,k) is the pruned graph for parameter k
 *
 * Edge frequency metrics are calculated as follows:
 * - For each edge, its relative frequency is the proportion of k values where it appears
 * - For each k value, the frequency_ratio is the average relative frequency of all edges
 *   in that specific pruned graph
 * - Higher frequency_ratio values indicate more stable graph configurations where edges
 *   tend to persist across different k values
 *
 * @param s_X SEXP object representing the input data matrix. Must be a numeric matrix
 *           (not a data frame) where:
 *           - Rows represent data points (vertices)
 *           - Columns represent features
 *           - Values must be numeric (coercible to REALSXP)
 *
 * @param s_kmin SEXP object (integer) representing the minimum k value to consider
 *              Must be positive
 *
 * @param s_kmax SEXP object (integer) representing the maximum k value to consider
 *              Must be greater than kmin
 *
 * @param s_compute_full SEXP object (logical) controlling computation of optional components:
 *                      - TRUE: Store all pruned graphs in the output
 *                      - FALSE: Only return edge statistics and birth/death times
 *
 * @param s_pruning_thld SEXP object (double) controlling intensity of the geometric edge pruning.
 *        Edge weight relative deviation is computed as
 *
 *        rel_dev = (w_ij + w_jk) / w_ik - 1.0;
 *
 *        Geometric pruning is performed on all edges with rel_dev < pruning_thld
 *
 * @param s_verbose SEXP object (logical) controlling progress reporting during computation
 *
 * @return SEXP object (a named list) containing:
 * - birth_death_matrix: Matrix with columns:
 *   - start: Starting vertex of edge (1-based)
 *   - end: Ending vertex of edge (1-based)
 *   - birth_time: k value when edge first appears
 *   - death_time: k value when edge permanently disappears (SIZE_MAX if never dies)
 *
 * - k_statistics: Matrix with columns:
 *   - k: k value
 *   - n_edges: Total edges in original graph
 *   - n_edges_in_pruned_graph: Edges remaining after pruning
 *   - n_removed_edges: Edges removed during pruning
 *   - edge_reduction_ratio: Proportion of edges removed
 *   - frequency_ratio: Average persistence of edges across k values (stability metric)
 *
 * - pruned_graphs: If compute_full is TRUE, list of pruned graphs for each k value
 *                 Each graph contains:
 *   - pruned_adj_list: Adjacency lists after pruning
 *   - pruned_weight_list: Distances for pruned edges
 *                 If compute_full is FALSE, R_NilValue
 *
 * @note
 * - The function uses the create_iknn_graph function directly for efficiency
 * - Graph pruning is performed using iknn_graph.prune_graph()
 * - Edges are stored with start < end to avoid duplication
 * - Birth/death times use 0-based k values internally
 * - All indices in returned R objects are 1-based
 * - The frequency_ratio can be used to select the most stable pruned graph configuration
 * - Higher frequency_ratio values indicate graphs whose edges tend to persist across multiple k values
 *
 * @warning
 * - Input matrix must be numeric and cannot be a data frame
 * - kmin must be positive and less than kmax
 * - Large k values or dense datasets may require significant memory
 *
 * @see
 * - create_iknn_graph(): Computes single IWD-kNN graph
 * - iknn_graph_t::prune_graph(): Graph pruning method
 * - create_R_graph_representation(): Converts internal graph format to R
 *
 */
SEXP orig_S_create_iknn_graphs(
    SEXP s_X,
    SEXP s_kmin,
    SEXP s_kmax,
    SEXP s_pruning_thld,
    SEXP s_compute_full,
    SEXP s_verbose) {

    auto total_start_time = std::chrono::steady_clock::now();

    double pruning_thld = REAL(s_pruning_thld)[0];
    int verbose = LOGICAL(s_verbose)[0];

    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;

    int kmin = INTEGER(s_kmin)[0];
    int kmax = INTEGER(s_kmax)[0];
    int compute_full = LOGICAL(s_compute_full)[0];
    int* dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];

    if (verbose) {
        Rprintf("Processing k values from %d to %d for %d vertices\n", kmin, kmax, n_vertices);
    }

    // Prepare vectors to store results
    std::vector<std::vector<double>> k_statistics(kmax - kmin + 1);

    // Create vector of k values
    std::vector<int> k_values(kmax - kmin + 1);
    std::iota(k_values.begin(), k_values.end(), kmin);

    // Parallel processing of graph creation and pruning
    if (verbose) Rprintf("Starting parallel graph processing\n");
    auto parallel_start_time = std::chrono::steady_clock::now();

    // Progress tracking
    std::atomic<int> progress_counter{0};
    //size_t n_k_values = k_values.size();
    //const size_t progress_chunk = std::max<size_t>(1, n_k_values / 5);  // Report every 5% progress
    const size_t n_k_values = kmax - kmin + 1;

    // Compute kNN once for maximum k
    auto knn_results = compute_knn(s_X, kmax);

    #define USE_GEOMETRIC_PRUNING 1

    #if USE_GEOMETRIC_PRUNING
    std::vector<set_wgraph_t> pruned_graphs(kmax - kmin + 1);
    #else
    int max_alt_path_length = 2; //INTEGER(s_max_alt_path_length)[0];
    std::vector<vect_wgraph_t> pruned_graphs(kmax - kmin + 1);
    #endif

    // Parallel processing using pre-computed kNN results
    std::for_each(std::execution::par_unseq,
                  k_values.begin(),
                  k_values.end(),
                  [&](int k) {

                      if (verbose) {
                          int current_count = ++progress_counter;
                          //if (current_count % progress_chunk == 0) {
                          REprintf("\rProcessing %d %d%%",
                                   current_count,
                                   static_cast<int>((100.0 * current_count) / n_k_values));
                          //}
                      }

                      auto iknn_graph = create_iknn_graph(knn_results, k);

                      // Count original edges
                      size_t n_edges = 0;
                      for (const auto& vertex_edges : iknn_graph.graph) {
                          n_edges += vertex_edges.size();
                      }
                      n_edges /= 2;


                      #if USE_GEOMETRIC_PRUNING
                      // geometric graph pruning
                      // transfering iknn_graph_t to set_wgraph_t
                      auto pruned_graph = set_wgraph_t(iknn_graph);

                      // Compute the deviations using the optimized method
                      auto rel_deviations = pruned_graph.compute_edge_weight_rel_deviations();
                      //const double EPSILON = 1e-16; // Threshold for considering a value as zero

                      // Track removed edges
                      //int removed_count = 0;

                      // Process edges
                      for (size_t i = 0; i < rel_deviations.size(); i++) {

                          //if (rel_deviations[i].rel_deviation < EPSILON) {
                          if (rel_deviations[i].rel_deviation < pruning_thld) {
                              size_t source = rel_deviations[i].source;
                              size_t target = rel_deviations[i].target;
                              //size_t intermediate = rel_deviations[i].best_intermediate;

                              // Check if removing this edge would isolate any vertex
                              if (pruned_graph.adjacency_list[source].size() <= 1) {
                                  continue;
                                  // REPORT_ERROR("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                                  //              source + 1, target + 1, intermediate + 1, source + 1);
                              }

                              if (pruned_graph.adjacency_list[target].size() <= 1) {
                                  continue;
                                  // REPORT_ERROR("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                                  //              source + 1, target + 1, intermediate + 1, target + 1);
                              }

                              // Safe to remove the edge
                              pruned_graph.remove_edge(source, target);
                              //removed_count++;
                          }
                      }

                      /// <<<--- here insert the code for long edge pruning - create a graph derived from pruned_graph by running prune_long_edges() on it


                      pruned_graphs[k - kmin] = std::move(pruned_graph);

                      // Count edges in the pruned graph
                      size_t n_edges_in_pruned_graph = 0;
                      for (const auto& neighbors : pruned_graph.adjacency_list) {
                          n_edges_in_pruned_graph += neighbors.size();
                      }
                      n_edges_in_pruned_graph /= 2;

                      k_statistics[k - kmin] = {
                          (double)n_edges,
                          (double)n_edges_in_pruned_graph,
                          (double)(n_edges - n_edges_in_pruned_graph),
                          (double)(n_edges - n_edges_in_pruned_graph) / n_edges
                      };




                      #else
                      // isize graph pruning
                      vect_wgraph_t pruned_graph = iknn_graph.prune_graph(max_alt_path_length);

                      // Count edges in the pruned graph
                      size_t n_edges_in_pruned_graph = 0;
                      for (const auto& neighbors : pruned_graph.adjacency_list) {
                          n_edges_in_pruned_graph += neighbors.size();
                      }
                      n_edges_in_pruned_graph /= 2;

                      // Store results
                      pruned_graphs[k - kmin] = std::move(pruned_graph);
                      k_statistics[k - kmin] = {
                          (double)n_edges,
                          (double)n_edges_in_pruned_graph,
                          (double)(n_edges - n_edges_in_pruned_graph),
                          (double)(n_edges - n_edges_in_pruned_graph) / n_edges
                      };
                      #endif
                  });

    if (verbose) {
        elapsed_time(parallel_start_time, "\rParallel processing completed", true);
        Rprintf("Processing edge brith/death time ... ");
    }

    // Serial processing of results
    auto serial_start_time = std::chrono::steady_clock::now();

    // Process birth/death times and edge frequencies
    std::unordered_map<edge_t, birth_death_time_t, edge_hash> birth_death_map;
    std::unordered_map<edge_t, size_t, edge_hash> edge_freq_map;

    for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
        int k = kmin + k_idx;
        const auto& pruned_graph = pruned_graphs[k_idx];

        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};

                    // Update birth/death times
                    auto it = birth_death_map.find(edge_key);
                    if (it == birth_death_map.end()) {
                        birth_death_map[edge_key] = {(size_t)k, (size_t)(k + 1)};
                    } else {
                        it->second.death_time = (size_t)(k + 1);
                    }

                    // Update edge frequency
                    edge_freq_map[edge_key]++;
                }
            }
        }
    }

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        serial_start_time = std::chrono::steady_clock::now();
        Rprintf("Processing edge frequencies ... ");
    }

    // Compute relative frequencies and k-specific statistics
    std::unordered_map<edge_t, double, edge_hash> edge_rel_freq_map;
    std::vector<double> k_freq_ratios(n_k_values);

    for (const auto& [edge, freq] : edge_freq_map) {
        edge_rel_freq_map[edge] = static_cast<double>(freq) / n_k_values;
    }

    // Calculate frequency ratios for each k
    for (int k_idx = 0; k_idx < n_k_values; k_idx++) {
        const auto& pruned_graph = pruned_graphs[k_idx];
        double total_rel_freq = 0.0;
        size_t edge_count = 0;

        for (int i = 0; i < n_vertices; i++) {
            for (const auto& edge : pruned_graph.adjacency_list[i]) {
                int neighbor = edge.vertex;
                if (i < neighbor) {
                    edge_t edge_key = {(size_t)i, (size_t)neighbor};
                    total_rel_freq += edge_rel_freq_map[edge_key];
                    edge_count++;
                }
            }
        }

        k_freq_ratios[k_idx] = edge_count > 0 ? total_rel_freq / edge_count : 0.0;
    }

    // Add frequency statistics to k_statistics
    for (size_t i = 0; i < k_statistics.size(); i++) {
        k_statistics[i].push_back(k_freq_ratios[i]);
    }

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        serial_start_time = std::chrono::steady_clock::now();
        Rprintf("Creating return list objects ... ");
    }

    // Create R objects
    SEXP pruned_graphs_list = R_NilValue;
    if (compute_full) {
        PROTECT(pruned_graphs_list = allocVector(VECSXP, kmax - kmin + 1)); nprot++;

        for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
            const auto& pruned_graph = pruned_graphs[k_idx];

            SEXP r_pruned_adj_list = PROTECT(allocVector(VECSXP, n_vertices));
            SEXP r_pruned_weight_list = PROTECT(allocVector(VECSXP, n_vertices));

            for (int i = 0; i < n_vertices; i++) {
                const auto& edges = pruned_graph.adjacency_list[i];

                {
                    SEXP RA = PROTECT(allocVector(INTSXP, edges.size()));
                    int* A = INTEGER(RA);
                    for (const auto& edge : edges) {
                        *A++ = edge.vertex + 1;
                    }
                    SET_VECTOR_ELT(r_pruned_adj_list, i, RA);
                    UNPROTECT(1);
                }

                {
                    SEXP RD = PROTECT(allocVector(REALSXP, edges.size()));
                    double* D = REAL(RD);
                    for (const auto& edge : edges) {
                        *D++ = edge.weight;
                    }
                    SET_VECTOR_ELT(r_pruned_weight_list, i, RD);
                    UNPROTECT(1);
                }
            }

            SEXP r_pruned_graph = PROTECT(allocVector(VECSXP, 2));
            SET_VECTOR_ELT(r_pruned_graph, 0, r_pruned_adj_list);
            SET_VECTOR_ELT(r_pruned_graph, 1, r_pruned_weight_list);

            SEXP r_pruned_graph_names = PROTECT(allocVector(STRSXP, 2));
            SET_STRING_ELT(r_pruned_graph_names, 0, mkChar("pruned_adj_list"));
            SET_STRING_ELT(r_pruned_graph_names, 1, mkChar("pruned_weight_list"));
            setAttrib(r_pruned_graph, R_NamesSymbol, r_pruned_graph_names);

            SET_VECTOR_ELT(pruned_graphs_list, k_idx, r_pruned_graph);
            UNPROTECT(4);
        }
    }

    // Create statistics matrix
    // if (verbose) Rprintf("\nCreating final output...");

    SEXP k_stats_matrix = PROTECT(allocMatrix(REALSXP, k_statistics.size(), 6)); nprot++;
    double* data = REAL(k_stats_matrix);

    for (size_t i = 0; i < k_statistics.size(); i++) {
        data[i] = kmin + i;
        for (size_t j = 0; j < 5; j++) {  // Now 5 statistics instead of 4
            data[i + (j + 1) * k_statistics.size()] = k_statistics[i][j];
        }
    }

    // Set column names
    SEXP k_stats_dimnames = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(k_stats_dimnames, 0, R_NilValue);

    SEXP k_stats_colnames = PROTECT(allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(k_stats_colnames, 0, mkChar("k"));
    SET_STRING_ELT(k_stats_colnames, 1, mkChar("n_edges"));
    SET_STRING_ELT(k_stats_colnames, 2, mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 3, mkChar("n_removed_edges"));
    SET_STRING_ELT(k_stats_colnames, 4, mkChar("edge_reduction_ratio"));
    SET_STRING_ELT(k_stats_colnames, 5, mkChar("frequency_ratio"));
    SET_VECTOR_ELT(k_stats_dimnames, 1, k_stats_colnames);

    setAttrib(k_stats_matrix, R_DimNamesSymbol, k_stats_dimnames);

    // Create birth/death matrix
    SEXP birth_death_matrix = PROTECT(convert_birth_death_map_to_matrix(birth_death_map)); nprot++;

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

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        elapsed_time(total_start_time, "Total elapsed time", true);
    }

    UNPROTECT(nprot);
    return result;
}
