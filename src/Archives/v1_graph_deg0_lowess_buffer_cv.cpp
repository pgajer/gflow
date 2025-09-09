////// old
graph_deg0_lowess_buffer_cv_t set_wgraph_t::graph_deg0_lowess_buffer_cv(
    const std::vector<double>& y,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool log_grid,
    size_t kernel_type,
    double dist_normalization_factor,
	bool use_uniform_weights,
    size_t buffer_hops,
    bool auto_buffer_hops,
    size_t n_folds,
    bool with_bw_predictions,
    double precision,
    bool verbose) {

    auto start_time = std::chrono::steady_clock::now();

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    size_t n_vertices = adjacency_list.size();

    if (verbose) {
        Rprintf("Starting graph_deg0_lowess_buffer_cv() with Buffer Zone Method\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Number of CV folds: %zu\n", n_folds);
        Rprintf("Buffer zone hops: %zu\n", buffer_hops);
    }

    // Initialize result structure
    graph_deg0_lowess_buffer_cv_t result;

    // If auto_buffer_hops is true, determine the optimal buffer size based on autocorrelation
    result.buffer_hops_used = buffer_hops;
    if (auto_buffer_hops) {
        if (verbose) {
            Rprintf("Auto-determining buffer hop distance based on spatial autocorrelation...\n");
        }
        buffer_hops = determine_optimal_buffer_hops(y, verbose);
        result.buffer_hops_used = buffer_hops;
        if (verbose) {
            Rprintf("Auto-determined buffer hop distance: %zu\n", buffer_hops);
        }
    }

    // 1. Create spatially stratified CV folds using maximal packing
    if (verbose) {
        Rprintf("Creating spatially stratified CV folds...\n");
    }

    // Calculate target number of seeds based on graph size
    double cv_seed_factor = n_vertices < 500 ? 0.2 :
                            n_vertices < 5000 ? 0.1 : 0.05;

    size_t target_seeds = std::max(size_t(n_vertices * cv_seed_factor), n_folds * 3);
    target_seeds = std::min(target_seeds, n_vertices / 10);

    if (verbose) {
        Rprintf("Creating maximal packing with target of %zu seeds for %zu folds\n",
                target_seeds, n_folds);
    }

    // Select seed points using maximal packing
    size_t max_iterations = 100;
    std::vector<size_t> fold_seeds = create_maximal_packing(
        target_seeds,
        max_iterations,
        precision);

    // Ensure we have enough seeds
    size_t actual_seeds = fold_seeds.size();
    if (verbose) {
        Rprintf("Created maximal packing with %zu seed points\n", actual_seeds);
    }

    // Determine actual number of folds based on available seeds
    size_t actual_folds = std::min(n_folds, actual_seeds);
    if (actual_folds < n_folds && verbose) {
        Rprintf("Warning: Could only create %zu folds (requested %zu)\n",
                actual_folds, n_folds);
    }

    // Create CV folds by grouping seeds if necessary
    std::vector<std::vector<size_t>> seed_groups;
    if (actual_seeds <= n_folds) {
        seed_groups.resize(actual_seeds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i].push_back(fold_seeds[i]);
        }
    } else {
        seed_groups.resize(n_folds);
        for (size_t i = 0; i < actual_seeds; i++) {
            seed_groups[i % n_folds].push_back(fold_seeds[i]);
        }
    }

    // Assign vertices to folds based on closest seed group
    std::vector<std::vector<size_t>> folds(actual_folds);
    std::vector<size_t> fold_assignment(n_vertices, SIZE_MAX);

    for (size_t vertex = 0; vertex < n_vertices; ++vertex) {
        double min_distance = std::numeric_limits<double>::infinity();
        size_t closest_fold = 0;

        for (size_t fold_idx = 0; fold_idx < seed_groups.size(); ++fold_idx) {
            for (size_t seed : seed_groups[fold_idx]) {
                double distance = compute_shortest_path_distance(vertex, seed);
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_fold = fold_idx;
                }
            }
        }

        folds[closest_fold].push_back(vertex);
        fold_assignment[vertex] = closest_fold;
    }

    if (verbose) {
        Rprintf("Created %zu folds with sizes: ", actual_folds);
        for (size_t i = 0; i < actual_folds; ++i) {
            Rprintf("%zu ", folds[i].size());
        }
        Rprintf("\n");
    }

    // 2. Define minimum and maximum bandwidth
    if (graph_diameter <= 0) {
        REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
    }

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    if (verbose) {
        Rprintf("Graph diameter: %f\n", graph_diameter);
        Rprintf("Bandwidth range: [%.4f, %.4f]\n", min_bw, max_bw);
    }

    // Generate bandwidth grid
    std::vector<double> bw_grid = get_candidate_bws(
        min_bw, max_bw, n_bws, log_grid, precision);

    if (verbose) {
        Rprintf("Testing %zu candidate bandwidths\n", bw_grid.size());
    }

    result.bws = bw_grid;
    result.bw_errors.resize(bw_grid.size(), 0.0);
    if (with_bw_predictions) {
        result.bw_predictions.resize(bw_grid.size(),
                                     std::vector<double>(n_vertices, 0.0));
    }

    // 3. Perform cross-validation with buffer zones for each fold
    if (verbose) {
        Rprintf("Starting cross-validation across %zu bandwidths with buffer zone method...\n",
                bw_grid.size());
    }

    // Setup for multi-threaded processing
    std::mutex mutex;
    progress_tracker_t progress(actual_folds, "Processing folds", 1);
    std::atomic<size_t> completed_tasks(0);

    // Process each fold
    std::vector<size_t> fold_indices(actual_folds);
    std::iota(fold_indices.begin(), fold_indices.end(), 0);

    std::vector<size_t> fold_test_counts(actual_folds, 0);
    std::vector<double> fold_total_errors(actual_folds, 0.0);

    std::for_each(std::execution::par_unseq, fold_indices.begin(), fold_indices.end(),
                  [&](size_t fold) {
        try {
            if (verbose) {
                std::lock_guard<std::mutex> lock(mutex);
                Rprintf("Processing fold %zu/%zu\n", fold+1, actual_folds);
            }

            // Get test vertices for this fold
            const std::vector<size_t>& test_vertices = folds[fold];

            // NEW: Create buffer zone for this fold
            std::unordered_set<size_t> buffer_zone = create_buffer_zone(
                test_vertices, buffer_hops);

            // Generate weights:
            // - 0.0 for test vertices and buffer zone (excluded from training)
            // - 1.0 for training vertices (outside buffer zone)
            std::vector<double> fold_weights(n_vertices, 1.0);
            for (size_t vertex : test_vertices) {
                fold_weights[vertex] = 0.0;  // Test vertex
            }
            for (size_t vertex : buffer_zone) {
                fold_weights[vertex] = 0.0;  // Buffer zone vertex
            }

            // Process each bandwidth
            for (size_t bw_idx = 0; bw_idx < bw_grid.size(); ++bw_idx) {
                double bandwidth = bw_grid[bw_idx];

                // Modified prediction function that properly handles buffer zones
                auto predictions = predict_with_buffer_zone(
                    bandwidth,
                    fold_weights,
                    y,
                    test_vertices,
                    buffer_zone,
                    dist_normalization_factor,
                    use_uniform_weights);

                // Calculate error on test vertices
                double fold_error = 0.0;
                size_t test_count = 0;

                for (size_t i = 0; i < test_vertices.size(); ++i) {
                    size_t vertex = test_vertices[i];
                    if (!predictions[i].is_excluded) {  // Ensure prediction was possible
                        fold_error += std::pow(predictions[i].value - y[vertex], 2);
                        test_count++;
                    }
                }

                // Store fold results
                if (test_count > 0) {
                    std::lock_guard<std::mutex> lock(mutex);
                    result.bw_errors[bw_idx] += fold_error;
                    fold_test_counts[fold] = test_count;
                    fold_total_errors[fold] += fold_error;
                }

                // Store predictions if requested
                if (with_bw_predictions) {
                    std::lock_guard<std::mutex> lock(mutex);
                    for (size_t i = 0; i < test_vertices.size(); ++i) {
                        size_t vertex = test_vertices[i];
                        if (!predictions[i].is_excluded) {
                            result.bw_predictions[bw_idx][vertex] = predictions[i].value;
                        }
                    }
                }
            }

            // Update progress
            if (verbose) {
                size_t tasks_completed = completed_tasks.fetch_add(1, std::memory_order_relaxed) + 1;
                std::lock_guard<std::mutex> lock(mutex);
                progress.update(tasks_completed);
            }
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(mutex);
            REPORT_ERROR("Error processing fold %zu: %s", fold, e.what());
        }
    });

    // Normalize errors by total test points
    size_t total_test_points = 0;
    for (auto count : fold_test_counts) {
        total_test_points += count;
    }

    if (total_test_points > 0) {
        for (auto& error : result.bw_errors) {
            error /= total_test_points;
        }
    }

    // Find the optimal bandwidth
    auto min_it = std::min_element(result.bw_errors.begin(), result.bw_errors.end());
    result.opt_bw_idx = min_it - result.bw_errors.begin();
    result.opt_bw = bw_grid[result.opt_bw_idx];

    if (verbose) {
        Rprintf("Optimal bandwidth index: %zu\n", result.opt_bw_idx);
        Rprintf("Optimal bandwidth: %.4f\n", result.opt_bw);
        Rprintf("CV error: %.6f\n", result.bw_errors[result.opt_bw_idx]);
    }

    // Final prediction with optimal bandwidth
    std::vector<double> uniform_weights(n_vertices, 1.0);  // Use all vertices
    result.predictions = predict_all_with_optimal_bandwidth(
        result.opt_bw, y, dist_normalization_factor, use_uniform_weights);

    if (verbose) {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        Rprintf("Completed graph_deg0_lowess_cv with Buffer Zone Method in %.2f seconds\n", elapsed);
    }

    return result;
}
