
/**
 * @brief Computes kernel-weighted means and Bayesian bootstrap errors for a graph function.
 *
 * @details For each vertex in the graph, this function:
 * 1. Calculates the kernel-weighted mean using neighboring vertices:
 *    kmean[i] = Σ(w[j] * K(d[i,j]) * y[j]) / Σ(w[j] * K(d[i,j]))
 *    where:
 *    - w[j] are Bayesian bootstrap weights
 *    - K(d[i,j]) are kernel weights based on graph distances
 *    - y[j] are function values
 *
 * 2. Computes the Bayesian bootstrap error:
 *    bb_error[i] = K(w[i]/max(w)) * |y[i] - kmean[i]|
 *
 * The error kernel K must satisfy K(0) = 1 and K(1) = 0.
 *
 * @param graph Adjacency lists representing the graph structure
 * @param edge_lengths Corresponding edge lengths
 * @param weights Bayesian bootstrap weights (positive, summing to 1)
 * @param y Function values at vertices
 * @param mean_kernel Kernel type for mean calculation
 * @param error_kernel Kernel type for error calculation (must satisfy K(0)=1, K(1)=0)
 * @param n_cores Number of cores to use for parallel computation:
 *                - n_cores = 1: serial execution
 *                - n_cores > 1: parallel execution with specified number of cores
 * @param dist_normalization_factor Factor for scaling distances (default: 1.01)
 * @param epsilon Numerical stability threshold (default: 1e-15)
 *
 * @return pair<vector<double>, double> containing:
 *         - first: kernel-weighted means at each vertex
 *         - second: total Bayesian bootstrap error (sum of individual errors)
 *
 * @throws std::invalid_argument if input sizes are inconsistent or if
 *         error kernel doesn't satisfy K(0)=1 and K(1)=0
 */
std::pair<std::vector<double>, double>
graph_kmean_with_bb_weights_and_error(
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<double>& weights,
    const std::vector<double>& y,
    int kmm_kernel,
    int error_kernel,
    int n_cores = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15) {

    const size_t n_vertices = graph.size();
    const size_t max_neighbors = std::max_element(
        graph.begin(), graph.end(),
        [](const auto& a, const auto& b) { return a.size() < b.size(); }
        )->size();

    // Validate inputs
    if (n_vertices != y.size() || n_vertices != weights.size()) {
        throw std::invalid_argument("Inconsistent input sizes");
    }

    // Store kernel function pointers separately
    auto error_kernel_fn = kernel_fn;
    auto kmm_kernel_fn = kernel_fn;

    // Initialize and validate error kernel
    initialize_kernel(error_kernel, 1.0);

    auto validate_error_kernel = [&](int kernel_type, double epsilon) {
        std::vector<double> x(1);
        std::vector<double> k_val(1);

        // Check K(0) = 1
        x[0] = 0.0;
        kernel_fn(x.data(), 1, k_val.data());
        if (std::fabs(1.0 - k_val[0]) > epsilon) {
            throw std::invalid_argument("Error kernel must satisfy K(0) = 1");
        }

        // Check K(1) = 0
        x[0] = 1.0;
        kernel_fn(x.data(), 1, k_val.data());
        if (std::fabs(k_val[0]) > epsilon) {
            throw std::invalid_argument("Error kernel must satisfy K(1) = 0");
        }
    };

    validate_error_kernel(error_kernel, epsilon);
    error_kernel_fn = kernel_fn;

    // Compute error weights using error kernel
    const double max_weight = *std::max_element(weights.begin(), weights.end());
    std::vector<double> normalized_weights(n_vertices);
    std::vector<double> error_weights(n_vertices);

    for (size_t i = 0; i < n_vertices; ++i) {
        normalized_weights[i] = weights[i] / max_weight;
    }
    error_kernel_fn(normalized_weights.data(), n_vertices, error_weights.data());

    // Initialize kmm kernel
    initialize_kernel(kmm_kernel, 1.0);
    kmm_kernel_fn = kernel_fn;

    // Preallocate vectors for results and temporary calculations
    std::vector<double> kmean(n_vertices);
    std::vector<double> bb_errors(n_vertices);
    std::vector<double> kmm_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

    // Define compute lambda
    auto compute_vertex_mean_and_error = [&, kmm_kernel_fn](size_t i) {
        distances[0] = 0;
        double max_dist = 0.0;
        for (size_t j = 0; j < graph[i].size(); ++j) {
            distances[j + 1] = edge_lengths[i][j];
            if (distances[j + 1] > max_dist)
                max_dist = distances[j + 1];
        }
        if (max_dist == 0) max_dist = 1;
        max_dist *= dist_normalization_factor;

        for (size_t j = 0; j < graph[i].size(); ++j)
            distances[j + 1] /= max_dist;

        int n = graph[i].size() + 1;
        kmm_kernel_fn(distances.data(), n, kmm_weights.data());

        double weight_sum = weights[i] * kmm_weights[0];
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            weight_sum += weights[neighbor] * kmm_weights[j + 1];
        }

        if (weight_sum < epsilon) {
            kmean[i] = y[i];
            bb_errors[i] = 0.0;
            return;
        }

        double weighted_sum = y[i] * weights[i] * kmm_weights[0];
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            weighted_sum += weights[neighbor] * kmm_weights[j + 1] * y[neighbor];
        }
        kmean[i] = weighted_sum / weight_sum;

        bb_errors[i] = error_weights[i] * std::abs(y[i] - kmean[i]);
    };

    // Main computation loop - switches to serial if n_cores = 1
    if (n_cores == 1) {
        // Serial execution
        for (size_t i = 0; i < n_vertices; ++i) {
            compute_vertex_mean_and_error(i);
        }
    } else {
        // Parallel execution
        omp_set_num_threads(n_cores); // Set number of threads for OpenMP
        #pragma omp parallel for private(kmm_weights, distances) if(n_vertices > 1000)
        for (size_t i = 0; i < n_vertices; ++i) {
            compute_vertex_mean_and_error(i);
        }
    }

    // Compute total error
    double total_bb_error = std::accumulate(bb_errors.begin(), bb_errors.end(), 0.0);

    return {kmean, total_bb_error};
}

/**
 * @brief Computes kernel-weighted means and Bayesian bootstrap errors for multiple bootstrap weights.
 *
 * @details This function implements a kernel-weighted mean calculation on a graph with
 * multiple sets of Bayesian bootstrap weights. For each bootstrap replicate and vertex,
 * it:
 * 1. Calculates the kernel-weighted mean using neighboring vertices:
 *    kmean[i] = Σ(w[j] * K(d[i,j]) * y[j]) / Σ(w[j] * K(d[i,j]))
 *    where:
 *    - w[j] are Bayesian bootstrap weights
 *    - K(d[i,j]) are kernel weights based on graph distances
 *    - y[j] are function values
 *
 * 2. Computes the Bayesian bootstrap error for each replicate:
 *    bb_error[i] = K(w[i]/max(w)) * |y[i] - kmean[i]|
 *
 * The error kernel K must satisfy K(0) = 1 and K(1) = 0.
 *
 * @param graph Adjacency lists representing the graph structure
 * @param edge_lengths Corresponding edge lengths for each adjacency in the graph
 * @param weights_vect Vector of bootstrap weight vectors, each summing to 1
 * @param y Function values at vertices
 * @param kmm_kernel Kernel type for mean calculation (see kernels.h for options)
 * @param error_kernel Kernel type for error calculation (must satisfy K(0)=1, K(1)=0)
 * @param n_cores Number of cores to use for parallel computation:
 *                - n_cores = 1: serial execution
 *                - n_cores > 1: parallel execution with specified number of cores
 * @param dist_normalization_factor Factor for scaling distances (default: 1.01)
 * @param epsilon Numerical stability threshold (default: 1e-15)
 *
 * @return std::pair containing:
 *         - first: Vector of vectors with kernel-weighted means for each bootstrap replicate
 *         - second: Vector of total Bayesian bootstrap errors for each replicate
 *
 * @throws R error if:
 *         - Input sizes are inconsistent
 *         - Bootstrap weights vectors have inconsistent sizes
 *         - Error kernel doesn't satisfy K(0)=1 and K(1)=0
 *
 * @note The function supports both serial and parallel execution based on n_cores parameter.
 *       Parallel execution uses OpenMP and is automatically enabled when n_cores > 1 and
 *       the problem size is sufficiently large.
 */
std::pair<std::vector<std::vector<double>>, std::vector<double>>
graph_kmean_with_bb_weights_and_errors(
    const std::vector<std::vector<int>>& graph,
    const std::vector<std::vector<double>>& edge_lengths,
    const std::vector<std::vector<double>>& weights_vect,
    const std::vector<double>& y,
    int kmm_kernel,
    int error_kernel,
    int n_cores = 1,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-15) {

    const size_t n_vertices = graph.size();
    const size_t n_bootstrap = weights_vect.size();
    const size_t max_neighbors = std::max_element(
        graph.begin(), graph.end(),
        [](const auto& a, const auto& b) { return a.size() < b.size(); }
        )->size();

    // Validate inputs
    if (n_vertices != y.size() || weights_vect[0].size() != n_vertices) {
        error("Inconsistent input sizes");
    }
    for (const auto& weights : weights_vect) {
        if (weights.size() != n_vertices) {
            error("Inconsistent bootstrap weights sizes");
        }
    }

    // Initialize kernels
    initialize_kernel(error_kernel, 1.0);
    auto error_kernel_fn = kernel_fn;
    // Validate error kernel properties
    {
        std::vector<double> x(1), k_val(1);

        x[0] = 0.0;
        error_kernel_fn(x.data(), 1, k_val.data());
        if (std::fabs(1.0 - k_val[0]) > epsilon) {
            error("Error kernel must satisfy K(0) = 1");
        }

        x[0] = 1.0;
        error_kernel_fn(x.data(), 1, k_val.data());
        if (std::fabs(k_val[0]) > epsilon) {
            error("Error kernel must satisfy K(1) = 0");
        }
    }

    initialize_kernel(kmm_kernel, 1.0);
    auto kmm_kernel_fn = kernel_fn;

    // Preallocate results
    std::vector<std::vector<double>> all_kmeans(n_bootstrap, std::vector<double>(n_vertices));
    std::vector<double> all_errors(n_bootstrap, 0.0);

    // Precompute max weights for each bootstrap replicate
    std::vector<double> max_weights(n_bootstrap);
    std::vector<std::vector<double>> all_error_weights(n_bootstrap, std::vector<double>(n_vertices));

    if (n_cores > 1) {
        #pragma omp parallel for num_threads(n_cores)
        for(size_t b = 0; b < n_bootstrap; ++b) {
            max_weights[b] = *std::max_element(weights_vect[b].begin(), weights_vect[b].end());
        }

        #pragma omp parallel for collapse(2) if(n_bootstrap * n_vertices > 1000) num_threads(n_cores)
        for(size_t b = 0; b < n_bootstrap; ++b) {
            for(size_t i = 0; i < n_vertices; ++i) {
                std::vector<double> x(1, weights_vect[b][i] / max_weights[b]);
                std::vector<double> k_val(1);
                error_kernel_fn(x.data(), 1, k_val.data());
                all_error_weights[b][i] = k_val[0];
            }
        }

        // Main parallel computation loop
        #pragma omp parallel if(n_bootstrap * n_vertices > 1000) num_threads(n_cores)
        {
            // Thread-local storage for temporary calculations
            std::vector<double> thread_kmm_weights(max_neighbors + 1);
            std::vector<double> thread_distances(max_neighbors + 1);

            // Thread-local computation lambda
            auto compute_vertex_mean_and_error = [&](size_t b, size_t i) {
                thread_distances[0] = 0;
                double max_dist = 0.0;

                // Compute and normalize distances
                for (size_t j = 0; j < graph[i].size(); ++j) {
                    thread_distances[j + 1] = edge_lengths[i][j];
                    if (thread_distances[j + 1] > max_dist) {
                        max_dist = thread_distances[j + 1];
                    }
                }
                if (max_dist == 0) max_dist = 1;
                max_dist *= dist_normalization_factor;

                for (size_t j = 0; j < graph[i].size(); ++j) {
                    thread_distances[j + 1] /= max_dist;
                }

                // Compute kernel weights
                int n = graph[i].size() + 1;
                kmm_kernel_fn(thread_distances.data(), n, thread_kmm_weights.data());

                // Calculate weighted sum and normalization
                double weight_sum = weights_vect[b][i] * thread_kmm_weights[0];
                double weighted_sum = y[i] * weights_vect[b][i] * thread_kmm_weights[0];

                for (size_t j = 0; j < graph[i].size(); ++j) {
                    int neighbor = graph[i][j];
                    double w = weights_vect[b][neighbor] * thread_kmm_weights[j + 1];
                    weight_sum += w;
                    weighted_sum += w * y[neighbor];
                }

                // Compute mean and error
                if (weight_sum < epsilon) {
                    all_kmeans[b][i] = y[i];
                    return 0.0;
                }

                all_kmeans[b][i] = weighted_sum / weight_sum;
                return all_error_weights[b][i] * std::abs(y[i] - all_kmeans[b][i]);
            };

            // Parallel computation of means and errors
            #pragma omp for collapse(2) schedule(dynamic)
            for(size_t b = 0; b < n_bootstrap; ++b) {
                for(size_t i = 0; i < n_vertices; ++i) {
                    double vertex_error = compute_vertex_mean_and_error(b, i);
                    #pragma omp atomic
                    all_errors[b] += vertex_error;
                }
            }
        }
    } else {
        // Serial execution
        std::vector<double> kmm_weights(max_neighbors + 1);
        std::vector<double> distances(max_neighbors + 1);

        // Compute max weights
        for(size_t b = 0; b < n_bootstrap; ++b) {
            max_weights[b] = *std::max_element(weights_vect[b].begin(), weights_vect[b].end());
        }

        // Compute error weights
        for(size_t b = 0; b < n_bootstrap; ++b) {
            for(size_t i = 0; i < n_vertices; ++i) {
                std::vector<double> x(1, weights_vect[b][i] / max_weights[b]);
                std::vector<double> k_val(1);
                error_kernel_fn(x.data(), 1, k_val.data());
                all_error_weights[b][i] = k_val[0];
            }
        }

        // Serial computation lambda
        auto compute_vertex_mean_and_error = [&](size_t b, size_t i) {
            distances[0] = 0;
            double max_dist = 0.0;

            // Compute and normalize distances
            for (size_t j = 0; j < graph[i].size(); ++j) {
                distances[j + 1] = edge_lengths[i][j];
                if (distances[j + 1] > max_dist) {
                    max_dist = distances[j + 1];
                }
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            for (size_t j = 0; j < graph[i].size(); ++j) {
                distances[j + 1] /= max_dist;
            }

            // Compute kernel weights
            int n = graph[i].size() + 1;
            kmm_kernel_fn(distances.data(), n, kmm_weights.data());

            // Calculate weighted sum and normalization
            double weight_sum = weights_vect[b][i] * kmm_weights[0];
            double weighted_sum = y[i] * weights_vect[b][i] * kmm_weights[0];

            for (size_t j = 0; j < graph[i].size(); ++j) {
                int neighbor = graph[i][j];
                double w = weights_vect[b][neighbor] * kmm_weights[j + 1];
                weight_sum += w;
                weighted_sum += w * y[neighbor];
            }

            // Compute mean and error
            if (weight_sum < epsilon) {
                all_kmeans[b][i] = y[i];
                return 0.0;
            }

            all_kmeans[b][i] = weighted_sum / weight_sum;
            return all_error_weights[b][i] * std::abs(y[i] - all_kmeans[b][i]);
        };

        // Main computation loop
        for(size_t b = 0; b < n_bootstrap; ++b) {
            for(size_t i = 0; i < n_vertices; ++i) {
                all_errors[b] += compute_vertex_mean_and_error(b, i);
            }
        }
    }

    return {all_kmeans, all_errors};
}

#if 0
std::vector<std::vector<double>> mrunif_simplex(int K, int n_bootstraps, int seed);

//  * @param graph The adjacency matrix representation of the graph as a vector of vectors
//  * @param edge_lengths Matrix containing the lengths/weights of edges between vertices
//  * @param y Vector of observed values associated with graph vertices
//    @param s_y_true SEXP containing ground truth values for performance evaluation; can be null, if it is not null, we compute true errors
//  * @param dist_normalization_factor Factor used to normalize distances (default: 1.01)
//  * @param epsilon Small value to prevent numerical instability (default: 1e-15)
//  * @param p Confidence level for the bootstrap confidence intervals (default: 0.95)

// * @param s_graph SEXP containing adjacency list representation of graph
// * @param s_edge_lengths SEXP containing edge lengths for each vertex
// * @param s_y SEXP containing initial values at vertices

// sampling_method:
// 0 - samples from the uniform distribution over the simplex
// 1 - homogeneous Dirichlet with parameter alpha
// 2 - beta radial

// if s_with_cv_errors = TRUE generate also cross-validation errors

// add y_true possibly null
// if not null, we produce true errors

// @param use_median If true, uses median for location estimation; if false, uses mean (default: false)

SEXP S_graph_kmean_with_bb_weights_and_errors(SEXP s_graph,
                                              SEXP s_edge_lengths,
                                              SEXP s_y,
                                              SEXP s_y_true,
                                              SEXP s_with_cv_errors,
                                              SEXP s_sampling_method,
                                              SEXP s_alpha,
                                              SEXP s_beta,
                                              SEXP s_error_method,
                                              SEXP s_use_median,
                                              SEXP s_p,
                                              SEXP s_kmm_kernel,
                                              SEXP s_error_kernel,
                                              SEXP s_n_bootstraps,
                                              SEXP s_n_cores,
                                              SEXP s_dist_normalization_factor,
                                              SEXP s_epsilon,
                                              SEXP s_seed) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);

    int n_vertices = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_vertices);

    int nprot = 0;
    PROTECT(s_y_true = coerceVector(s_y_true, REALSXP)); nprot++;
    std::vector<double> y_true(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));

    bool with_cv_errors = LOGICAL(s_with_cv_errors)[0];
    int sampling_method = INTEGER(s_sampling_method)[0];
    double alpha = REAL(s_alpha)[0];
    double beta = REAL(s_beta)[0];
    int error_method = INTEGER(s_error_method)[0];
    bool use_median = LOGICAL(s_use_median)[0];
    double p = REAL(s_p)[0];
    int kmm_kernel = INTEGER(s_kmm_kernel)[0];
    int error_kernel = INTEGER(s_error_kernel)[0];
    int n_bootstraps = INTEGER(s_n_bootstraps)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];

    // Generating Bayesian bootstrap weights
    std::vector<std::vector<double>> weights_vect;
    switch(sampling_method) {
    case 0:
        weights_vect = mrunif_simplex(n_vertices n_bootstraps, seed);
        // for (int i = 0; i < n_bootstraps; i++) {
        //     C_runif_simplex(&n_vertices, weights_vect[i].data());
        break;
    case 1:
        for (int i = 0; i < n_bootstraps; i++)
            weights_vect[i] = sample_symmetric_dirichlet(n_vertices, alpha);
        break;
    case 2:
        for (int i = 0; i < n_bootstraps; i++)
            weights_vect[i] = beta_radial(n_vertices, alpha, beta, seed);
        break;
    default:
        Rf_error("Invalid method specified");
    }

    // Computing kernel-weighted means and Bayesian bootstrap errors for weights_vect Bayesian bootstrap weights
    auto res = graph_kmean_with_bb_weights_and_errors(graph,
                                                      edge_lengths,
                                                      weights_vect,
                                                      y,
                                                      kmm_kernel,
                                                      error_kernel,
                                                      n_cores,
                                                      dist_normalization_factor,
                                                      epsilon);

    // Converting y bootstraps to a list of vectors
    auto y_bootstraps = res.first;
    SEXP s_y_bootstraps = convert_vector_vector_double_to_R(y_bootstraps);

    // Converting Bayesian bootstrap errors to an R vector of errors
    SEXP s_bb_errors = convert_vector_double_to_R(res.second);

    // Compute location of y_bootstraps
    // Helper function to compute location (mean or median) of a sequence
    auto compute_location = [use_median](const std::vector<double>& values) -> double {
        if (!use_median) {
            return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        } else {
            std::vector<double> sorted_values = values;
            std::sort(sorted_values.begin(), sorted_values.end());
            size_t n = sorted_values.size();
            if (n % 2 == 0) {
                return (sorted_values[n/2 - 1] + sorted_values[n/2]) / 2.0;
            } else {
                return sorted_values[n/2];
            }
        }
    };

    std::vector<double> y_locations(n_vertices);
    std::vector<double> vertex_estimates(n_bootstraps);
    for (int i = 0; i < n_vertices; ++i) {
        for (int b = 0; b < n_bootstraps; ++b) {
            vertex_estimates[b] = y_bootstraps[b][i];
        }
        y_locations[i] = compute_location(vertex_estimates);
    }

    // Turning y_locations into an R vector
    SEXP s_y_locations = convert_vector_double_to_R(y_locations);

    // Computing credible intervals
    auto cri = bb_cri(y_bootstraps, y_locations, p);

    // Turning cri into two R vectors
    SEXP s_cri_L = convert_vector_double_to_R(cri.first);
    SEXP s_cri_U = convert_vector_double_to_R(cri.second);

    SEXP s_cv_errors;
    if (with_cv_errors) { // Computing CV errors
        std::vector<double> cv_errors =  graph_kmean_cv(graph,
                                                        edge_lengths,
                                                        y,
                                                        ikernel,
                                                        dist_normalization_factor,
                                                        n_CVs,
                                                        n_CV_folds,
                                                        epsilon,
                                                        seed);

        s_cv_errors = convert_vector_double_to_R(cv_errors);
    }

    // Computing true error in y_true is given
    SEXP s_true_error;
    if (y_true.size() == n_vertices) {
        double true_error = 0;
        for (int i = 0; i < n_vertices; ++i) {
            true_error += std::abs(y[i] - y_true[i]);
        }
        true_error /= n_vertices;

        s_true_error = PROTECT(allocVector(REALSXP, 1));
        double* ptr = REAL(s_true_error);
        ptr[0] = true_error;
        UNPROTECT(1);
    }

    // Do all these SEXP's created so far need to be protected??? <<----

    // Creating a return list and its names
    const int N_RESULTS = 7;
    SEXP results = PROTECT(allocVector(VECSXP, N_RESULTS)); nprot++;
    SET_VECTOR_ELT(results, 0, s_y_bootstraps);
    SET_VECTOR_ELT(results, 1, s_bb_errors);
    SET_VECTOR_ELT(results, 2, s_y_locations);
    SET_VECTOR_ELT(results, 3, s_cri_L);
    SET_VECTOR_ELT(results, 4, s_cri_U);
    SET_VECTOR_ELT(results, 5, s_cv_errors);
    SET_VECTOR_ELT(results, 6, s_true_error);
    // Set names
    SEXP names = PROTECT(allocVector(STRSXP, N_RESULTS)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("y_bootstraps"));
    SET_STRING_ELT(names, 1, mkChar("bb_errors"));
    SET_STRING_ELT(names, 2, mkChar("y_locations"));
    SET_STRING_ELT(names, 3, mkChar("cri_L"));
    SET_STRING_ELT(names, 4, mkChar("cri_U"));
    SET_STRING_ELT(names, 5, mkChar("cv_errors"));
    SET_STRING_ELT(names, 6, mkChar("true_error"));

    setAttrib(results, R_NamesSymbol, names);
    UNPROTECT(nprot);

    return results;
}
#endif



/**
 * @brief Performs Bayesian bootstrap estimation of kernel-weighted means and MAD errors.
 *
 * @details This function performs multiple Bayesian bootstrap iterations to estimate the
 * distribution of both kernel-weighted means and Mean Absolute Deviation (MAD) errors
 * for each vertex in the graph. For each bootstrap iteration, it:
 * 1. Generates weights using the ordered differences method
 * 2. Computes kernel-weighted means and MAD errors
 * 3. Stores results in separate vectors for means and errors
 *
 * @param graph Graph structure as adjacency lists
 * @param edge_lengths Edge lengths corresponding to the graph structure
 * @param y Function values at vertices
 * @param n_bb Number of Bayesian bootstrap iterations
 * @param ikernel Kernel function identifier
 * @param dist_normalization_factor Factor for distance normalization (default: 1.01)
 * @param epsilon Threshold for zero weight sums (default: 1e-15)
 *
 * @return std::pair containing:
 *         - first: vector of n_bb vectors of kernel-weighted means
 *         - second: vector of n_bb vectors of MAD errors
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
graph_kmean_bb_and_mad(const std::vector<std::vector<int>>& graph,
                       const std::vector<std::vector<double>>& edge_lengths,
                       const std::vector<double>& y,
                       int n_bb,
                       int ikernel,
                       double dist_normalization_factor = 1.01,
               double epsilon = 1e-15) {
    int n_points = y.size();
    std::vector<double> weights(n_points);
    std::vector<std::vector<double>> bb_Ey(n_bb);
    std::vector<std::vector<double>> bb_MAD(n_bb);

    for (int iboot = 0; iboot < n_bb; iboot++) {
        C_runif_simplex(&n_points, weights.data());
        auto [means, mads] = graph_kmean_with_bb_weigths_and_mad(graph,
                                                                 edge_lengths,
                                                                 weights,
                                                                 y,
                                                                 ikernel,
                                                                 dist_normalization_factor,
                                                                 epsilon);
        bb_Ey[iboot] = std::move(means);
        bb_MAD[iboot] = std::move(mads);
    }

    return {bb_Ey, bb_MAD};
}

/**
 * @brief Computes the total Bayesian bootstrap mean absolute deviation error from Bayesian bootstrap results
 *
 * @param bb_results The pair of vectors returned by graph_kmean_bb_and_mad, where
 *                   the second component contains the MAD errors for each bootstrap iteration
 * @return The sum of all MAD errors across all bootstrap iterations and vertices
 */
double compute_total_bb_mad(
    const std::pair<std::vector<std::vector<double>>,
    std::vector<std::vector<double>>>& bb_results) {

    double total_error = 0.0;

    // Get the vector of MAD errors (second component)
    const auto& mad_errors = bb_results.second;

    // Sum over all bootstrap iterations
    for (const auto& iteration_errors : mad_errors) {
        // Sum over all vertices in this iteration
        total_error += std::accumulate(iteration_errors.begin(),
                                       iteration_errors.end(),
                                       0.0);
    }

    return total_error;
}

#if 0
error_t graph_kmean_errors(const std::vector<std::vector<int>>& graph,
                           const std::vector<std::vector<double>>& edge_lengths,
                           const std::vector<double>& y,
                           const std::vector<double>& y_true,
                           int n_bb,
                           int ikernel,
                           double dist_normalization_factor = 1.01,
                           double epsilon = 1e-15,
                           bool use_median = true,
                           int n_CVs = 0,
                           int n_CV_folds = 10,
                           unsigned int seed = 0,
                           const error_methods_t method = error_methods_t::ALL,
                           double p = 0.95) {

    int n_points = y.size();
    error_t results;

    std::vector<std::vector<double>> bb_Ey;
    std::pair<std::vector<double>, std::vector<double>> Ey_cri;
    if (method == error_methods_t::ALL ||
        method == error_methods_t::BB_COV_ERRORS ||
        method == error_methods_t::BB_WASSERSTEIN_ERROR) {

        // Computing Bayesian bootstraps
        bb_Ey = graph_kmean_bb(graph,
                               edge_lengths,
                               y,
                               n_bb,
                               ikernel,
                               dist_normalization_factor,
                               epsilon);

        // Helper function to compute location (mean or median) of a sequence
        auto compute_location = [use_median](const std::vector<double>& values) -> double {
            if (!use_median) {
                return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
            } else {
                std::vector<double> sorted_values = values;
                std::sort(sorted_values.begin(), sorted_values.end());
                size_t n = sorted_values.size();
                if (n % 2 == 0) {
                    return (sorted_values[n/2 - 1] + sorted_values[n/2]) / 2.0;
                } else {
                    return sorted_values[n/2];
                }
            }
        };

        // Computing location of bb_Ey
        std::vector<double> mu_location(n_points, 0.0);
        for (int i = 0; i < n_points; ++i) {
            std::vector<double> vertex_estimates(B);
            for (int b = 0; b < B; ++b) {
                vertex_estimates[b] = bb_Ey[b][i];
            }
            mu_location[i] = compute_location(vertex_estimates);
        }
        results.Ey = std::move(mu_location);

        // Computing credible intervals
        results.Ey_cri = bb_cri(bb_Ey, y, p);
    }

    switch(method) {
    case error_methods_t::ALL:
        std::pair<double, double> res = compute_bbcov_error(bb_Ey, results.Ey, y, use_median);
        results.bb_cov_errorA = res.first;
        results.bb_cov_errorB = res.second;
        results.bb_wasserstein_error = compute_bbwasserstein_error(bb_Ey, results.Ey, y);
        results.cv_error = graph_kmean_cv(graph,
                                          edge_lengths,
                                          y,
                                          ikernel,
                                          dist_normalization_factor,
                                          n_CVs,
                                          n_CV_folds,
                                          epsilon,
                                          seed);

        if (y_true.size() == y.size()) {


        }

        break;
    case error_methods_t::BB_COV_ERRORS:
        std::pair<double, double> res = compute_bbcov_error(bb_Ey, results.Ey, y, use_median);
        results.bb_cov_errorA = res.first;
        results.bb_cov_errorB = res.second;
        break;
    case error_methods_t::BB_WASSERSTEIN_ERROR:
        results.bb_wasserstein_error = compute_bbwasserstein_error(bb_Ey, results.Ey, y);
        break;
    case error_methods_t::CV_ERROR:
        results.cv_error = graph_kmean_cv(graph,
                                          edge_lengths,
                                          y,
                                          ikernel,
                                          dist_normalization_factor,
                                          n_CVs,
                                          n_CV_folds,
                                          epsilon,
                                          seed);
        break;
    default:
        Rf_error("Invalid method specified");
    }

    return results;
}
#endif


/**
 * @brief Computes kernel-weighted means and MAD errors for a graph with Bayesian Bootstrap weights.
 *
 * @details This function calculates both the kernel-weighted means and Mean Absolute Deviation (MAD)
 * errors for a function over graph vertices, incorporating Bayesian Bootstrap weights. For each vertex,
 * it computes:
 * 1. kmean[i] = sum(weights[j] * kernel_weights[j] * y[j]) / sum(weights[j] * kernel_weights[j])
 * 2. mad[i] = sum(weights[j] * kernel_weights[j] * |y[j] - kmean[i]|) / sum(weights[j] * kernel_weights[j])
 *
 * @param graph Graph structure as adjacency lists
 * @param edge_lengths Edge lengths corresponding to the graph structure
 * @param weights Bayesian bootstrap weights (positive, sum to 1)
 * @param y Function values at vertices
 * @param ikernel Kernel function identifier
 * @param dist_normalization_factor Factor for distance normalization (default: 1.01)
 * @param epsilon Threshold for zero weight sums (default: 1e-15)
 *
 * @return std::pair containing:
 *         - first: vector of kernel-weighted means
 *         - second: vector of MAD errors
 */
std::pair<std::vector<double>, std::vector<double>>
graph_kmean_with_bb_weigths_and_mad(const std::vector<std::vector<int>>& graph,
                                    const std::vector<std::vector<double>>& edge_lengths,
                                    const std::vector<double>& weights,
                                    const std::vector<double>& y,
                                    int ikernel,
                                    double dist_normalization_factor = 1.01,
                                    double epsilon = 1e-15) {

    auto kmean = std::vector<double>(y.size(), 0.0);
    auto mad = std::vector<double>(y.size(), 0.0);
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);

    for (size_t i = 0; i < graph.size(); ++i) {
        distances[0] = 0;
        double max_dist = 0.0;
        for (size_t j = 0; j < graph[i].size(); ++j) {
            distances[j + 1] = edge_lengths[i][j];
            if (distances[j + 1] > max_dist)
                max_dist = distances[j + 1];
        }
        if (max_dist == 0) max_dist = 1;
        max_dist *= dist_normalization_factor;
        for (size_t j = 0; j < graph[i].size(); ++j)
            distances[j + 1] /= max_dist;

        int n = graph[i].size() + 1;
        kernel_fn(distances.data(), n, kernel_weights.data());

        double weight_sum = weights[i] * kernel_weights[0];
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            weight_sum += weights[neighbor] * kernel_weights[j + 1];
        }

        if (weight_sum < epsilon) {
            kmean[i] = y[i];
            mad[i] = 0.0;
            continue;
        }

        double weighted_sum = y[i] * weights[i] * kernel_weights[0];
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            weighted_sum += weights[neighbor] * kernel_weights[j + 1] * y[neighbor];
        }
        kmean[i] = weighted_sum / weight_sum;

        // Calculate MAD
        double mad_sum = weights[i] * kernel_weights[0] * std::abs(y[i] - kmean[i]);
        for (size_t j = 0; j < graph[i].size(); ++j) {
            int neighbor = graph[i][j];
            mad_sum += weights[neighbor] * kernel_weights[j + 1] *
                      std::abs(y[neighbor] - kmean[i]);
        }
        mad[i] = mad_sum / weight_sum;
    }

    return {kmean, mad};
}




/**
 * @brief Computes graph kernel weighted mean using Bayesian bootstrap and BBMWD.
 *
 * This function performs the following steps:
 * 1. Computes the graph kernel weighted mean using Bayesian bootstrap and credible intervals (CrI).
 * 2. Calculates the Bayesian Bootstrap Mean Wasserstein Distance (BBMWD) between the bootstrap estimates and observed values.
 *
 * @param s_graph SEXP representing the graph structure as a list of integer vectors.
 *                Each vector contains the indices of neighboring vertices for a given vertex.
 * @param s_edge_lengths SEXP representing edge lengths as a list of double vectors.
 *                       Each vector contains the lengths of edges corresponding to the neighbors in the graph structure.
 * @param s_y SEXP representing a vector of observed values at each vertex of the graph.
 * @param s_n_bb SEXP representing the number of Bayesian bootstrap iterations to perform.
 * @param s_ikernel SEXP representing the integer identifier for the kernel function to use.
 * @param s_dist_normalization_factor SEXP representing the factor used to normalize distances in the kernel computation.
 *
 * @return SEXP A list with two named elements:
 *         - "bb_Ey": A matrix where each column represents a bootstrap sample and each row represents a vertex.
 *         - "bbmwd": A single numeric value representing the computed BBMWD.
 *
 * @note This function uses R's C API and assumes that appropriate R headers are included.
 * @note The function uses helper functions like Rgraph_to_vector and
 *       R_list_of_dvectors_to_cpp_vector_of_dvectors which should be defined elsewhere.
 *
 * @see graph_kmean_bb for the underlying Bayesian bootstrap computation.
 * @see compute_bbwasserstein_error for the BBMWD calculation.
 *
 * @warning This function may be computationally intensive for large graphs or high numbers of bootstrap iterations.
 * @warning Proper memory management with PROTECT/UNPROTECT is crucial when working with R's C API.
 */
SEXP S_graph_kmean_gof(SEXP s_graph,
                       SEXP s_edge_lengths,
                       SEXP s_y,
                       SEXP s_n_bb,
                       SEXP s_ikernel,
                       SEXP s_dist_normalization_factor) {
    // Converting R objects to C++ objects
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<double>> edge_lengths = Rweights_to_vector(s_edge_lengths);
    int y_length = LENGTH(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);
    int n_bb = INTEGER(s_n_bb)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];

    // Computing Bayesian bootstrap estimates
    std::vector<std::vector<double>> bb_Ey = graph_kmean_bb(graph,
                                                            edge_lengths,
                                                            y,
                                                            n_bb,
                                                            ikernel,
                                                            dist_normalization_factor);

    // Computing BBMWD
    double bbmwd = compute_bbwasserstein_error(bb_Ey, y);

    // Computing BBMWD vect
    // std::vector<double> bbmwd_vect = compute_bbwasserstein_error(bb_Ey, y);

    // Converting bb_Ey to SEXP (R matrix)
    int n_vertices = bb_Ey[0].size();
    SEXP s_bb_Ey = PROTECT(allocMatrix(REALSXP, n_vertices, n_bb)); // Protect the matrix allocation and population
    double *ptr = REAL(s_bb_Ey);
    for (int j = 0; j < n_bb; ++j) {
        for (int i = 0; i < n_vertices; ++i) {
            ptr[i + j * n_vertices] = bb_Ey[j][i];
        }
    }
    UNPROTECT(1); // s_bb_Ey is now referenced by the result list, safe to unprotect

    // Converting bbmwd to SEXP
    SEXP s_bbmwd = PROTECT(allocVector(REALSXP, 1));
    REAL(s_bbmwd)[0] = bbmwd;
    UNPROTECT(1);

    // Converting bbmwd_vect to SEXP
    #if 0
    SEXP s_bbmwd_vect = PROTECT(allocVector(REALSXP, n_vertices));
    for (int i = 0; i < n_vertices; ++i) {
        REAL(s_bbmwd_vect)[i] = bbmwd_vect[i];
    }
    UNPROTECT(1);
    #endif

    // Creating return list
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, s_bb_Ey);
    SET_VECTOR_ELT(result, 1, s_bbmwd);
    //SET_VECTOR_ELT(result, 2, s_bbmwd_vect);

    // Setting names for return list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("bb_Ey"));
    SET_STRING_ELT(names, 1, mkChar("bbmwd"));
    //SET_STRING_ELT(names, 2, mkChar("bbmwd_vect"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2);

    return result;
}
