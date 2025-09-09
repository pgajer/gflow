std::vector<std::vector<double>> gplm_bb(const path_graph_plm_t& path_graph,
                                         const std::vector<double>& y,
                                         int n_bb,
                                         int n_cores = 1,
                                         int ikernel = 1,
                                         double dist_normalization_factor = 1.01,
                                         double epsilon = 1e-15) {
    // Input validation
    int n_vertices = y.size();

    // Define parameters missing from function signature but needed for wwlm
    int max_distance_deviation = 2;  // These should probably be function parameters

    // Initialize result vector before the main
    std::vector<std::vector<double>> bb_Ey(n_bb);
    for (int i = 0; i < n_bb; ++i) {
        bb_Ey[i].resize(n_vertices);
    }

    // Set number of threads if parallel processing is requested
    if (n_cores > 1) {

        Rprintf("In the if (n_cores > 1) block of gplm_bb\n");

        omp_set_num_threads(n_cores);

        // Pre-allocate vectors before the parallel region
        std::vector<std::vector<double>> thread_weights(n_cores, std::vector<double>(n_vertices));

        #pragma omp parallel for schedule(dynamic)
        for (int iboot = 0; iboot < n_bb; ++iboot) {
            // Get thread-specific pre-allocated vector
            int thread_id = omp_get_thread_num();
            auto& weights = thread_weights[thread_id];
#if 1
            C_runif_simplex(&n_vertices, weights.data());

            auto res = gplm_with_precomputed_path_models(path_graph,
                                                         y,
                                                         weights,
                                                         ikernel,
                                                         max_distance_deviation,
                                                         dist_normalization_factor,
                                                         epsilon);

            #if 0
            auto res = gplm_with_global_bb_weights(path_graph,
                                                   y,
                                                   weights,
                                                   ikernel,
                                                   max_distance_deviation,
                                                   dist_normalization_factor,
                                                   epsilon);
            #endif

            bb_Ey[iboot] = std::move(res.first);
#else
            int min_required_paths = 3;
            auto res = gplm_with_bb_weights(path_graph,
                                            y,
                                            ikernel,
                                            max_distance_deviation,
                                            min_required_paths,
                                            dist_normalization_factor,
                                            epsilon);
            bb_Ey[iboot] = std::move(res.first);
#endif
        }

    } else {

        Rprintf("In the if (n_cores == 1) block of gplm_bb\n");

        for (int iboot = 0; iboot < n_bb; ++iboot) {
#if 1
            std::vector<double> weights(n_vertices);
            C_runif_simplex(&n_vertices, weights.data());

            auto res = gplm_with_global_bb_weights(path_graph,
                                                   y,
                                                   weights,
                                                   ikernel,
                                                   max_distance_deviation,
                                                   dist_normalization_factor,
                                                   epsilon);
            bb_Ey[iboot] = std::move(res.first);
#else
            int min_required_paths = 3;
            auto res = gplm_with_bb_weights(path_graph,
                                            y,
                                            ikernel,
                                            max_distance_deviation,
                                            min_required_paths,
                                            dist_normalization_factor,
                                            epsilon);
            bb_Ey[iboot] = std::move(res.first);
#endif
        }
    }

    #if 0
    // Parallel for loop with private weights vector
    #pragma omp parallel for if(n_cores > 1) schedule(dynamic)
    for (int iboot = 0; iboot < n_bb; ++iboot) {
        try {
            #if 1
            std::vector<double> weights(n_vertices);
            C_runif_simplex(&n_vertices, weights.data());

            auto res = gplm_with_global_bb_weights(path_graph,
                                                   y,
                                                   weights,
                                                   ikernel,
                                                   max_distance_deviation,
                                                   min_required_paths,
                                                   dist_normalization_factor,
                                                   epsilon);
            bb_Ey[iboot] = std::move(res.first);
            #else
            auto res = gplm_with_bb_weights(path_graph,
                                            y,
                                            ikernel,
                                            max_distance_deviation,
                                            min_required_paths,
                                            dist_normalization_factor,
                                            epsilon);
            bb_Ey[iboot] = std::move(res.first);
           #endif
        } catch (const std::exception& e) {
            Rf_error("Error in bootstrap iteration %d: %s", iboot, e.what());
        }
    }
    #endif

    return bb_Ey;
}

std::vector<std::vector<double>> gplm_bb(const path_graph_plm_t& path_graph,
                                         const std::vector<double>& y,
                                         int n_bb,
                                         int n_cores = 1,
                                         int ikernel = 1,
                                         double dist_normalization_factor = 1.01,
                                         double epsilon = 1e-15) {
    // Input validation
    int n_vertices = y.size();

    // Initialize result vector before the main
    std::vector<std::vector<double>> bb_Ey(n_bb);
    for (int i = 0; i < n_bb; ++i) {
        bb_Ey[i].resize(n_vertices);
    }

    std::vector<double> weights(n_vertices);
    if (n_cores > 1) {
        Rprintf("In the if (n_cores > 1) block of gplm_bb\n");
        omp_set_num_threads(n_cores);

        #pragma omp parallel for schedule(dynamic)
        for (int iboot = 0; iboot < n_bb; ++iboot) {
            C_runif_simplex(&n_vertices, weights.data());
            auto res = gplm_with_precomputed_path_models(path_graph,
                                                         y,
                                                         weights,
                                                         ikernel,
                                                         max_distance_deviation,
                                                         dist_normalization_factor,
                                                         epsilon);
            bb_Ey[iboot] = std::move(res.first);
        }
    } else {
        Rprintf("In the if (n_cores == 1) block of gplm_bb\n");
        for (int iboot = 0; iboot < n_bb; ++iboot) {
            C_runif_simplex(&n_vertices, weights.data());
            auto res = gplm_with_precomputed_path_models(path_graph,
                                                         y,
                                                         weights,
                                                         ikernel,
                                                         max_distance_deviation,
                                                         dist_normalization_factor,
                                                         epsilon);
            bb_Ey[iboot] = std::move(res.first);
        }
    }

    return bb_Ey;
}

adaptive_nbhd_size_plm_t gplm(const std::vector<std::vector<int>>& neighbors,
                              const std::vector<std::vector<double>>& edge_lengths,
                              const std::vector<double>& y,
                              const std::vector<double>& y_true,
                              //int max_distance_deviation,
                              bool use_median = false,
                              int h_min = 3,
                              int h_max = 31,
                              double p = 0.95,
                              int n_bb = 500,
                              int ikernel = 1,
                              int n_cores = 1,
                              double dist_normalization_factor = 1.01,
                              double epsilon = 1e-15,
                              unsigned int seed = 0,
                              bool verbose = true) {

    double total_ptm = (double)clock() / CLOCKS_PER_SEC; // for total elapsed time

    int n_vertices = static_cast<int>(y.size());

    // Should add at start of function:
    if (h_min > h_max) {
        Rf_error("h_min must be less than or equal to h_max");
    }
    if (h_min % 2 == 0 || h_max % 2 == 0) {
        Rf_error("h_min and h_max must be odd numbers");
    }
    if (neighbors.size() != edge_lengths.size() || neighbors.size() != y.size()) {
        Rf_error("Inconsistent input dimensions");
    }
    if (p <= 0 || p >= 1) {
        Rf_error("p must be between 0 and 1");
    }

    #define DEBUG__gplm 0
    #if DEBUG__gplm
    Rprintf("\nIn gplm()\n");
    Rprintf("Number of vertices: %d\n", n_vertices);
    print_vect_vect(neighbors, "neighbors");
    print_vect_vect(edge_lengths, "edge_lengths");
    #endif

    adaptive_nbhd_size_plm_t results;
    int n_h_values = (h_max - h_min) / 2 + 1; // h_max - h_min + 1 is not correct as it would be OK only if all numbers between h_min and h_max were used, but we are using only odd numbers within that range
    results.graphs.resize(n_h_values);
    results.h_cv_errors.resize(n_h_values);
    results.h_values.resize(n_h_values);

    std::unordered_map<int, std::vector<double>> errors_map;
    std::unordered_map<int, std::vector<double>> predictions_map;

    double ptm; // for elapsed time
    std::vector<double> weights(n_vertices, 1.0 / n_vertices);
    for (int i = 0, h = h_min; h <= h_max; h += 2, i++) { // note that h is incremented by 2 so that all h values are odd
        results.h_values[i] = h;

        // Creating an path graph for the given value of h
        if (verbose) {
            Rprintf("\nProcessing h = %d\n", h);
            Rprintf("\tCreating path graph ... ");
            ptm = (double)clock() / CLOCKS_PER_SEC;
        }
        auto path_graph = create_path_graph_plm(neighbors, edge_lengths, h);
        if (verbose) elapsed_time(ptm,"DONE");

        #if DEBUG__gplm
        Rprintf("\ni: %d\th: %d\n", i, h);
        #endif

        // Computing predictions and errors
        if (verbose) {
            Rprintf("\tComputing predictions and errors ... ");
            ptm = (double)clock() / CLOCKS_PER_SEC;
        }

        int local_max_distance_deviation = (h - 1) / 2;
        auto [predictions, errors] = gplm_with_precomputed_path_models(path_graph,
                                                                       y,
                                                                       weights,
                                                                       ikernel,
                                                                       local_max_distance_deviation,
                                                                       dist_normalization_factor,
                                                                       epsilon);
        if (verbose) elapsed_time(ptm,"DONE");

        #if DEBUG__gplm
        print_vect(predictions,"predictions");
        print_vect(errors,"errors");
        #endif

        // Calculating mean error across vertices
        if (verbose) {
            Rprintf("\tCalculating mean error across vertices ... ");
            ptm = (double)clock() / CLOCKS_PER_SEC;
        }
        double total_error = std::accumulate(errors.begin(), errors.end(), 0.0);
        results.h_cv_errors[i] = total_error / n_vertices;
        results.graphs[i] = std::move(path_graph);
        errors_map[h] = errors;
        predictions_map[h] = predictions;
        if (verbose) elapsed_time(ptm,"DONE");
    }

    #if DEBUG__gplm
    Rprintf("After main loop\n");
    #endif

    // Find the optimal h (minimum CV error)
    auto min_it = std::min_element(results.h_cv_errors.begin(), results.h_cv_errors.end());
    int opt_h_idx = std::distance(results.h_cv_errors.begin(), min_it);
    results.opt_h = h_min + 2 * opt_h_idx; // we are multiplying by 2 as we advance h by 2 in the above for loop

    #if DEBUG__gplm
    Rprintf("opt_h_idx: %d\tresults.opt_h: %d\n", opt_h_idx, results.opt_h);
    //error("Debug");
    #endif

    // Store optimal graph
    results.opt_h_graph = results.graphs[opt_h_idx];

    // Store the conditional expectations using optimal graph
    results.opt_predictions = predictions_map[results.opt_h];

    #if DEBUG__gplm
    Rprintf("results.opt_h: %d\n", results.opt_h);
    print_vect(results.opt_predictions,"results.opt_predictions");
    #endif

    // Computing local conditional expectations
    if (verbose) {
        Rprintf("Computing local conditional expectations ... ");
        ptm = (double)clock() / CLOCKS_PER_SEC;
    }
    results.opt_local_predictions.resize(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        // Find h with minimum error for vertex i
        double min_error = std::numeric_limits<double>::infinity();
        int local_opt_h = h_min;
        for (int h_idx = 0; h_idx < n_h_values; h_idx++) {
            int h = results.h_values[h_idx];
            if (errors_map[h][i] < min_error) {
                min_error = errors_map[h][i];
                local_opt_h = h;
            }
        }
        results.opt_local_predictions[i] = predictions_map[local_opt_h][i];
    }
    if (verbose) elapsed_time(ptm,"DONE");

    #if DEBUG__gplm
    print_vect(results.opt_local_predictions,"results.opt_local_predictions");
    #endif

    // Cleaning the memory for the maps
    if (verbose) {
        Rprintf("Cleaning up temporary data structures ... ");
        ptm = (double)clock() / CLOCKS_PER_SEC;
    }
    errors_map = std::unordered_map<int, std::vector<double>>();      // Force deallocation by replacement
    predictions_map = std::unordered_map<int, std::vector<double>>();
    if (verbose) elapsed_time(ptm, "DONE", false);

    // Computing true errors
    if (!y_true.empty() && y_true.size() == n_vertices) {
        if (verbose) {
            Rprintf("Computing true errors ... ");
            ptm = (double)clock() / CLOCKS_PER_SEC;
        }
        results.opt_true_errors.resize(n_vertices);
        for (size_t i = 0; i < n_vertices; i++) {
            results.opt_true_errors[i] = std::abs(y_true[i] - results.opt_predictions[i]);
        }
        if (verbose) elapsed_time(ptm,"DONE");
    } else {
        results.opt_true_errors.clear();  // Ensure empty if no true values
    }

    // Optional: Computing bootstrap credible intervals
    if (n_bb > 0) {
        if (verbose) {
            Rprintf("Computing bootstrap credible intervals ... ");
            ptm = (double)clock() / CLOCKS_PER_SEC;
        }
        bb_cri_t bb_cri_results = gplm_bb_cri(results.opt_h_graph,
                                              y,
                                              p,
                                              n_bb,
                                              use_median,
                                              n_cores,
                                              ikernel,
                                              dist_normalization_factor,
                                              epsilon);

        results.opt_bb_predictions = std::move(bb_cri_results.bb_Ey);
        results.opt_ci_lower = std::move(bb_cri_results.cri_L);
        results.opt_ci_upper = std::move(bb_cri_results.cri_U);
        if (verbose) elapsed_time(ptm,"DONE");
    }

    if (verbose) {
        Rprintf("Optimal h: %d\n", results.opt_h);
        elapsed_time(total_ptm,"Total Elapsed Time");
    }

    return results;
}

SEXP S_univariate_gplm(SEXP s_x,
                       SEXP s_y,
                       SEXP s_y_true,
                       SEXP s_max_distance_deviation,
                       SEXP s_use_median,
                       SEXP s_h_min,
                       SEXP s_h_max,
                       SEXP s_p,
                       SEXP s_n_bb,
                       SEXP s_ikernel,
                       SEXP s_n_cores,
                       SEXP s_dist_normalization_factor,
                       SEXP s_epsilon,
                       SEXP s_seed,
                       SEXP s_verbose) {

    int n_protected = 0;  // Track number of PROTECT calls

    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(s_y_true) == n_points) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    int max_distance_deviation = INTEGER(s_max_distance_deviation)[0];
    bool use_median = LOGICAL(s_use_median)[0];
    int h_min = INTEGER(s_h_min)[0];
    int h_max = INTEGER(s_h_max)[0];
    double p = REAL(s_p)[0];
    int n_bb = INTEGER(s_n_bb)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    auto cpp_results = univariate_gplm(x,
                                       y,
                                       y_true,
                                       max_distance_deviation,
                                       use_median,
                                       h_min,
                                       h_max,
                                       p,
                                       n_bb,
                                       ikernel,
                                       n_cores,
                                       dist_normalization_factor,
                                       epsilon,
                                       seed,
                                       verbose);
    // Creating return list
    const int N_COMPONENTS = 11;
    SEXP result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    SEXP s_h_values = convert_vector_int_to_R(cpp_results.h_values);
    SET_VECTOR_ELT(result, 0, s_h_values);

    if (!cpp_results.h_cv_errors.empty() && n_points > 0) {
        SEXP s_cv_errors = convert_vector_double_to_R(cpp_results.h_cv_errors);
        SET_VECTOR_ELT(result, 1, s_cv_errors);
    } else {
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    SEXP s_opt_h = PROTECT(allocVector(REALSXP, 1)); n_protected++;
    REAL(s_opt_h)[0] = cpp_results.opt_h;
    SET_VECTOR_ELT(result, 2, s_opt_h);

    SEXP s_opt_h_graph_adj_list = convert_vector_vector_int_to_R(cpp_results.opt_h_graph.adj_list);
    SET_VECTOR_ELT(result, 3, s_opt_h_graph_adj_list);

    SEXP s_opt_h_graph_edge_lengths = convert_vector_vector_double_to_R(cpp_results.opt_h_graph.weight_list);
    SET_VECTOR_ELT(result, 4, s_opt_h_graph_edge_lengths);

    SEXP s_predictions = convert_vector_double_to_R(cpp_results.opt_predictions);
    SET_VECTOR_ELT(result, 5, s_predictions);

    SEXP s_local_predictions = convert_vector_double_to_R(cpp_results.opt_local_predictions);
    SET_VECTOR_ELT(result, 6, s_local_predictions);

    if (cpp_results.opt_bb_predictions.size() > 0) {
        SEXP s_bb_predictions = convert_vector_double_to_R(cpp_results.opt_bb_predictions);
        SET_VECTOR_ELT(result, 7, s_bb_predictions);

        SEXP s_opt_ci_lower = convert_vector_double_to_R(cpp_results.opt_ci_lower);
        SET_VECTOR_ELT(result, 8, s_opt_ci_lower);

        SEXP s_opt_ci_upper = convert_vector_double_to_R(cpp_results.opt_ci_upper);
        SET_VECTOR_ELT(result, 9, s_opt_ci_upper);
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }

    if (cpp_results.opt_true_errors.size() > 0) {
        SEXP s_true_error = PROTECT(allocVector(REALSXP, 1)); n_protected++;
        double mean_true_error = std::accumulate(cpp_results.opt_true_errors.begin(),
                                                 cpp_results.opt_true_errors.end(), 0.0) /  cpp_results.opt_true_errors.size();
        REAL(s_true_error)[0] = mean_true_error;
        SET_VECTOR_ELT(result, 10, s_true_error);
    } else {
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("h_values"));
    SET_STRING_ELT(names, 1, mkChar("h_cv_errors"));
    SET_STRING_ELT(names, 2, mkChar("opt_h"));
    SET_STRING_ELT(names, 3, mkChar("opt_graph_adj_list"));
    SET_STRING_ELT(names, 4, mkChar("opt_graph_edge_lengths"));
    SET_STRING_ELT(names, 5, mkChar("opt_predictions"));
    SET_STRING_ELT(names, 6, mkChar("opt_local_predictions"));
    SET_STRING_ELT(names, 7, mkChar("opt_bb_predictions"));
    SET_STRING_ELT(names, 8, mkChar("opt_ci_lower"));
    SET_STRING_ELT(names, 9, mkChar("opt_ci_upper"));
    SET_STRING_ELT(names, 10, mkChar("true_error"));

    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protected);

    return result;
}
