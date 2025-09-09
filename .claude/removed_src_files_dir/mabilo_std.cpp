
#if 0
/**
 * @brief Implements a parallel version of Model-Averaged LOWESS (MABILO) using C++17 Standard Library parallel algorithms
 *
 * @details This implementation parallelizes the MABILO algorithm using std::execution policies and parallel algorithms.
 * The algorithm consists of three main phases, all utilizing parallel processing:
 * 1. Parallel model computation for different window sizes (k values)
 * 2. Optimal k selection using parallel reduction
 * 3. Parallel model averaging with error-based filtering and weighting
 *
 * Key parallel optimizations include:
 * - Use of std::execution::par_unseq for computationally intensive operations
 * - Parallel processing of reference points within each k iteration
 * - Parallel reductions for statistical calculations
 * - Parallel transformations for weight computations
 * - Thread-safe updates using atomic operations and mutexes where necessary
 *
 * @param x Vector of ordered x values (predictor variable)
 * @param y Observed y values corresponding to x (response variable)
 * @param y_true Optional true y values for error calculation (empty vector if not available)
 * @param w Bayesian bootstrap weights for observations
 * @param k_min Minimal number of nearest neighbors for window construction (window size = 2k + 1)
 * @param k_max Maximal number of nearest neighbors for window construction
 * @param model_averaging_strategy Strategy for combining model predictions:
 *        - KERNEL_WEIGHTS_ONLY: Uses only distance-based kernel weights
 *        - ERROR_WEIGHTS_ONLY: Uses only error-based weights
 *        - KERNEL_AND_ERROR_WEIGHTS: Combines both weight types
 * @param error_filtering_strategy Method for filtering models based on their errors:
 *        - GLOBAL_PERCENTILE: Uses global error distribution
 *        - LOCAL_PERCENTILE: Uses error distribution at each point
 *        - COMBINED_PERCENTILE: Applies both global and local filtering
 *        - BEST_K_MODELS: Selects k best models by error
 * @param distance_kernel Kernel function for distance-based weights:
 *        - 0: Tricube
 *        - 1: Epanechnikov
 *        - 2: Exponential
 * @param model_kernel Kernel function for error-based weights (same options as distance_kernel)
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Numerical stability parameter for weight calculations (default: 1e-15)
 * @param verbose Enable progress messages and timing information
 *
 * @return mabilo_t structure containing:
 *         - opt_k: Optimal k value selected
 *         - opt_k_idx: Index of optimal k value
 *         - k_mean_sm_errors: Mean LOOCV errors for each k
 *         - k_true_errors: True errors for each k (if y_true provided)
 *         - k_mean_true_errors: Mean true errors for each k
 *         - sm_predictions: Predictions before model averaging
 *         - predictions: Final smoothed values after model averaging
 *         - errors: LOOCV errors for final predictions
 *         - k_predictions: Predictions for all k values
 *
 * @throws Rf_error if input vectors are empty or of different sizes
 *
 * @note
 * - The implementation requires C++17 or later
 * - Input x values must be ordered. Sort x and reorder y accordingly if necessary
 * - Performance may vary based on the hardware and compiler implementation of parallel algorithms
 * - Progress tracking in verbose mode is thread-safe but may affect performance
 *
 * @example
 * ```cpp
 * std::vector<double> x = {1, 2, 3, 4, 5};
 * std::vector<double> y = {1.2, 1.9, 3.2, 4.1, 4.8};
 * std::vector<double> w(x.size(), 1.0);  // Equal weights
 * std::vector<double> y_true;  // Empty for no true values
 *
 * mabilo_t result = mabilo_std(
 *     x, y, y_true, w,
 *     2, 4,  // k_min, k_max
 *     KERNEL_AND_ERROR_WEIGHTS,
 *     COMBINED_PERCENTILE,
 *     0, 0,  // Tricube kernel for both
 *     1.01, 1e-15,  // Default factors
 *     true  // Enable verbose output
 * );
 * ```
 *
 * @see mabilo_t for the structure definition
 * @see initialize_kernel() for kernel function details
 * @see ulm_plus_t for the structure of individual linear models
 */
mabilo_t mabilo_std(const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& y_true,
                    const std::vector<double>& w,
                    int k_min,
                    int k_max,
                    model_averaging_strategy_t model_averaging_strategy,
                    error_filtering_strategy_t error_filtering_strategy,
                    int distance_kernel,
                    int model_kernel,
                    double dist_normalization_factor,
                    double epsilon,
                    bool verbose) {

    if (x.empty() || y.empty()) {
        REPORT_ERROR("Input vectors x and y cannot be empty");
    }
    if (x.size() != y.size()) {
        REPORT_ERROR("Input vectors x and y must have the same size");
    }

    int n_points = x.size();
    auto total_ptm = std::chrono::steady_clock::now();
    memory_tracker_t mem_tracker("MABILO");

    if (verbose) {
        Rprintf("Starting parallel MABILO computation (STL)\n");
        Rprintf("Input size: %d points\n", n_points);
        Rprintf("k range: %d to %d\n", k_min, k_max);
    }

    initialize_kernel(distance_kernel, 1.0);

    mabilo_t results;
    auto models_ptm = std::chrono::steady_clock::now();

    if (verbose) {
        Rprintf("\nPhase 1: Computing models for different k values\n");
    }

    progress_tracker_t k_progress(k_max - k_min + 1, "Model computation");

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    bool y_true_exists = !y_true.empty();

    int n_k_values = k_max - k_min + 1;
    results.k_mean_sm_errors.resize(n_k_values);
    std::vector<std::vector<double>> k_predictions(n_k_values, std::vector<double>(n_points));
    std::vector<std::vector<double>> k_errors(n_k_values, std::vector<double>(n_points));

    std::vector<std::vector<std::vector<ulm_plus_t>>> pt_k_models(n_k_values,
        std::vector<std::vector<ulm_plus_t>>(n_points));
    results.k_true_errors.resize(n_k_values, std::vector<double>(n_points));
    results.k_mean_true_errors.resize(n_k_values);

    // Generate k indices
    std::vector<int> k_indices(n_k_values);
    std::iota(k_indices.begin(), k_indices.end(), 0);

    // Phase 1: Parallel processing for each k value
    //std::for_each(std::execution::par_unseq, k_indices.begin(), k_indices.end(),
    std::for_each(GFLOW_EXEC_POLICY, k_indices.begin(), k_indices.end(),
        [&](int k_index) {
            int k = k_min + k_index;
            auto k_ptm = std::chrono::steady_clock::now();

            if (verbose) {
                Rprintf("\nProcessing k=%d (%d/%d) ... ",
                       k, k_index + 1, k_max - k_min + 1);
            }

            int n_points_minus_k = n_points - k;
            int k_minus_one = k - 1;
            int two_k = 2 * k;
            int n_points_minus_one_minus_two_k = n_points - 1 - two_k;

            // Generate reference point indices
            std::vector<int> ref_indices(n_points);
            std::iota(ref_indices.begin(), ref_indices.end(), 0);

            // Parallel processing of reference points
            std::for_each(GFLOW_EXEC_POLICY, ref_indices.begin(), ref_indices.end(),
                [&](int ref_pt_index) {
                    int x_min_index, x_max_index;

                    // Window bounds calculation
                    if (ref_pt_index > k_minus_one && ref_pt_index < n_points_minus_k) {
                        x_min_index = ref_pt_index - k;
                        x_max_index = ref_pt_index + k;
                    } else if (ref_pt_index < k) {
                        x_min_index = 0;
                        x_max_index = two_k;
                    } else {
                        x_min_index = n_points_minus_one_minus_two_k;
                        x_max_index = n_points - 1;
                    }

                    // Computing window weights
                    int window_size = x_max_index - x_min_index + 1;
                    std::vector<double> dists(window_size);
                    std::vector<double> weights(window_size);

                    // Calculate distances
                    std::transform(std::execution::par_unseq,
                        x.begin() + x_min_index, x.begin() + x_max_index + 1,
                        dists.begin(),
                        [&](double xi) { return std::abs(xi - x[ref_pt_index]); });

                    double max_dist = *std::max_element(std::execution::par_unseq,
                        dists.begin(), dists.end());

                    if (max_dist) {
                        max_dist *= dist_normalization_factor;
                        std::transform(std::execution::par_unseq,
                            dists.begin(), dists.end(), dists.begin(),
                            [max_dist](double d) { return d / max_dist; });
                    }

                    kernel_fn(dists.data(), window_size, weights.data());

                    double total_weights = std::reduce(std::execution::par_unseq,
                        weights.begin(), weights.end());

                    std::transform(std::execution::par_unseq,
                        weights.begin(), weights.end(),
                        w.begin() + x_min_index,
                        weights.begin(),
                        [total_weights](double weight, double wi) {
                            return (weight / total_weights) * wi;
                        });

                    // Fitting weighted linear model
                    ulm_t fit = ulm(x.data() + x_min_index,
                                    y.data() + x_min_index,
                                    weights,
                                    y_binary,
                                    epsilon);

                    ulm_plus_t wlm_res;
                    wlm_res.predictions = std::move(fit.predictions);
                    wlm_res.errors      = std::move(fit.errors);
                    wlm_res.x_min_index = x_min_index;
                    wlm_res.x_max_index = x_max_index;

                    k_errors[k_index][ref_pt_index] = wlm_res.errors[ref_pt_index - x_min_index];
                    double y_prediction_at_ref_pt = wlm_res.predictions[ref_pt_index - x_min_index];
                    k_predictions[k_index][ref_pt_index] = y_prediction_at_ref_pt;

                    if (y_true_exists) {
                        results.k_true_errors[k_index][ref_pt_index] =
                            std::abs(y_true[ref_pt_index] - y_prediction_at_ref_pt);
                    }

                    wlm_res.w = weights;

                    // Use mutex for thread-safe updates to pt_k_models
                    static std::mutex models_mutex;
                    for (int i = x_min_index; i <= x_max_index; i++) {
                        std::lock_guard<std::mutex> lock(models_mutex);
                        pt_k_models[k_index][i].push_back(wlm_res);
                    }
                });

            // Parallel reduction for mean error calculation
            results.k_mean_sm_errors[k_index] = std::reduce(std::execution::par_unseq,
                k_errors[k_index].begin(), k_errors[k_index].end()) / n_points;

            if (y_true_exists) {
                results.k_mean_true_errors[k_index] = std::reduce(std::execution::par_unseq,
                    results.k_true_errors[k_index].begin(),
                    results.k_true_errors[k_index].end()) / n_points;
            }

            if (verbose) {
                elapsed_time(k_ptm, "Done");
                mem_tracker.report();
                k_progress.update(k_index + 1);
            }
        });

    if (verbose) {
        elapsed_time(models_ptm, "\nTotal model computation time: ");
    }

    // Phase 2: Finding optimal k
    auto opt_k_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 2: Finding optimal k ... ");
    }

    auto min_it = std::min_element(std::execution::par_unseq,
        results.k_mean_sm_errors.begin(), results.k_mean_sm_errors.end());
    results.opt_k_idx = std::distance(results.k_mean_sm_errors.begin(), min_it);
    results.opt_k = k_min + results.opt_k_idx;
    results.sm_predictions = k_predictions[results.opt_k_idx];
    results.k_predictions = std::move(k_predictions);
    auto pt_models = std::move(pt_k_models[results.opt_k_idx]);
    std::vector<std::vector<std::vector<ulm_plus_t>>>().swap(pt_k_models);

    if (verbose) {
        elapsed_time(opt_k_ptm, "Done");
        mem_tracker.report();
    }

    // Phase 3: Parallel model averaging
    auto avg_ptm = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("\nPhase 3: Model averaging\n");
    }

    progress_tracker_t avg_progress(n_points, "Averaging progress");

    results.predictions.resize(n_points);
    results.errors.resize(n_points);

    initialize_kernel(model_kernel, 1.0);

    // Parallel computation of global error statistics if needed
    std::vector<double> all_errors;
    if (error_filtering_strategy == GLOBAL_PERCENTILE ||
        error_filtering_strategy == COMBINED_PERCENTILE) {

        // Collect errors in parallel
        std::vector<double> temp_errors;
        for (size_t i = 0; i < n_points; ++i) {
            for (const auto& model : pt_models[i]) {
                temp_errors.push_back(model.errors[i - model.x_min_index]);
            }
        }

        // Sort errors in parallel
        std::sort(std::execution::par_unseq, temp_errors.begin(), temp_errors.end());
        all_errors = std::move(temp_errors);
    }

    // Generate indices for parallel processing
    std::vector<size_t> indices(n_points);
    std::iota(indices.begin(), indices.end(), 0);

    // Parallel model averaging
    std::for_each(GFLOW_EXEC_POLICY, indices.begin(), indices.end(),
        [&](size_t i) {
            double weighted_sum = 0.0;
            double weight_sum = 0.0;
            std::vector<std::pair<double, const ulm_plus_t*>> filtered_models;

            switch (model_averaging_strategy) {
            case KERNEL_WEIGHTS_ONLY: {
                for (const auto& model : pt_models[i]) {
                    int local_index = i - model.x_min_index;
                    double weight = model.w[local_index];
                    weighted_sum += weight * model.predictions[local_index];
                    weight_sum += weight;
                }
                if (weight_sum > epsilon) {
                    results.predictions[i] = weighted_sum / weight_sum;
                } else {
                    auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                        [i](const auto& a, const auto& b) {
                            return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                        });
                    results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                }
                break;
            }

            case ERROR_WEIGHTS_ONLY: {
                // Filter models based on chosen strategy
                switch (error_filtering_strategy) {
                case GLOBAL_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.25 * all_errors.size())];
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case LOCAL_PERCENTILE: {
                    // Collect local errors for the current point
                    std::vector<double> local_errors;
                    local_errors.reserve(pt_models[i].size());
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        local_errors.push_back(model.errors[local_index]);
                    }

                    if (!local_errors.empty()) {
                        // Sort local errors in parallel
                        std::sort(std::execution::par_unseq, local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        // Filter models
                        filtered_models.reserve(pt_models[i].size());
                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case COMBINED_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.5 * all_errors.size())];
                    std::vector<double> local_errors;
                    local_errors.reserve(pt_models[i].size());

                    // First apply global filter and collect local errors
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            local_errors.push_back(model.errors[local_index]);
                        }
                    }

                    if (!local_errors.empty()) {
                        // Sort local errors in parallel
                        std::sort(std::execution::par_unseq, local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        // Apply local filter
                        filtered_models.reserve(pt_models[i].size());
                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case BEST_K_MODELS: {
                    const int k = 3;  // Number of best models to use
                    filtered_models.reserve(pt_models[i].size());

                    // Collect all models with their errors
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        filtered_models.push_back({model.errors[local_index], &model});
                    }

                    if (filtered_models.size() > k) {
                        // Partial sort to get k best models
                        std::partial_sort(filtered_models.begin(),
                                          filtered_models.begin() + k,
                                          filtered_models.end(),
                                          [](const auto& a, const auto& b) {
                                              return a.first < b.first;
                                          });
                        filtered_models.resize(k);
                    }
                    break;
                }
                }

                // Compute error-based weights for filtered models
                if (!filtered_models.empty()) {
                    // Find min and max errors using parallel reduction
                    auto [min_error_it, max_error_it] = std::minmax_element(
                        std::execution::par_unseq,
                        filtered_models.begin(), filtered_models.end(),
                        [](const auto& a, const auto& b) { return a.first < b.first; }
                        );
                    double min_error = min_error_it->first;
                    double max_error = max_error_it->first;
                    double error_range = max_error - min_error;

                    if (error_range > epsilon) {
                        // Process all models in parallel using transforms and reductions
                        std::vector<double> weights(filtered_models.size());
                        std::vector<double> weighted_preds(filtered_models.size());

                        // Compute normalized errors and weights
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       weights.begin(),
                                       [min_error, error_range](const auto& model_pair) {
                                           double normalized_error = (model_pair.first - min_error) / error_range;
                                           std::vector<double> error_weight(1);
                                           kernel_fn(&normalized_error, 1, error_weight.data());
                                           return error_weight[0];
                                       }
                            );

                        // Compute weighted predictions
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       weights.begin(),
                                       weighted_preds.begin(),
                                       [i](const auto& model_pair, double weight) {
                                           int local_index = i - model_pair.second->x_min_index;
                                           return weight * model_pair.second->predictions[local_index];
                                       }
                            );

                        // Sum up weights and weighted predictions
                        weighted_sum = std::reduce(std::execution::par_unseq, weighted_preds.begin(), weighted_preds.end());
                        weight_sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end());

                        if (weight_sum > epsilon) {
                            results.predictions[i] = weighted_sum / weight_sum;
                            results.errors[i] = std::abs(results.predictions[i] - y[i]);
                        } else {
                            // Fallback if weights sum to zero
                            auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                               [i](const auto& a, const auto& b) {
                                                                   return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                               });
                            results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                            results.errors[i] = best_model->errors[i - best_model->x_min_index];
                        }
                    } else {
                        // All errors are essentially the same, use simple average
                        weighted_sum = std::transform_reduce(std::execution::par_unseq,
                                                             filtered_models.begin(), filtered_models.end(),
                                                             0.0,
                                                             std::plus<>(),
                                                             [i](const auto& model_pair) {
                                                                 int local_index = i - model_pair.second->x_min_index;
                                                                 return model_pair.second->predictions[local_index];
                                                             }
                            );

                        weight_sum = filtered_models.size();
                        results.predictions[i] = weighted_sum / weight_sum;
                        results.errors[i] = std::abs(results.predictions[i] - y[i]);
                    }
                } else {
                    // If no models passed filtering, use the best available model
                    auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                       [i](const auto& a, const auto& b) {
                                                           return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                       });
                    results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                    results.errors[i] = best_model->errors[i - best_model->x_min_index];
                }
                break;
            }

            case KERNEL_AND_ERROR_WEIGHTS: {
                // Filter models based on chosen strategy
                switch (error_filtering_strategy) {
                case GLOBAL_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.25 * all_errors.size())];
                    filtered_models.reserve(pt_models[i].size());
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            filtered_models.push_back({model.errors[local_index], &model});
                        }
                    }
                    break;
                }

                case LOCAL_PERCENTILE: {
                    // Collect and sort local errors
                    std::vector<double> local_errors;
                    local_errors.reserve(pt_models[i].size());
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        local_errors.push_back(model.errors[local_index]);
                    }

                    if (!local_errors.empty()) {
                        std::sort(std::execution::par_unseq, local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        filtered_models.reserve(pt_models[i].size());
                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case COMBINED_PERCENTILE: {
                    double global_threshold = all_errors[static_cast<int>(0.5 * all_errors.size())];
                    std::vector<double> local_errors;
                    local_errors.reserve(pt_models[i].size());

                    // Apply global filter
                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        if (model.errors[local_index] <= global_threshold) {
                            local_errors.push_back(model.errors[local_index]);
                        }
                    }

                    if (!local_errors.empty()) {
                        std::sort(std::execution::par_unseq, local_errors.begin(), local_errors.end());
                        double local_threshold = local_errors[static_cast<int>(0.25 * local_errors.size())];

                        filtered_models.reserve(pt_models[i].size());
                        for (const auto& model : pt_models[i]) {
                            int local_index = i - model.x_min_index;
                            if (model.errors[local_index] <= local_threshold) {
                                filtered_models.push_back({model.errors[local_index], &model});
                            }
                        }
                    }
                    break;
                }

                case BEST_K_MODELS: {
                    const int k = 3;
                    filtered_models.reserve(pt_models[i].size());

                    for (const auto& model : pt_models[i]) {
                        int local_index = i - model.x_min_index;
                        filtered_models.push_back({model.errors[local_index], &model});
                    }

                    if (filtered_models.size() > k) {
                        std::partial_sort(filtered_models.begin(),
                                          filtered_models.begin() + k,
                                          filtered_models.end(),
                                          [](const auto& a, const auto& b) {
                                              return a.first < b.first;
                                          });
                        filtered_models.resize(k);
                    }
                    break;
                }
                }

                // Compute combined kernel and error-based weights for filtered models
                if (!filtered_models.empty()) {
                    // Find min and max errors using parallel algorithm
                    auto [min_error_it, max_error_it] = std::minmax_element(
                        std::execution::par_unseq,
                        filtered_models.begin(), filtered_models.end(),
                        [](const auto& a, const auto& b) { return a.first < b.first; }
                        );
                    double min_error = min_error_it->first;
                    double max_error = max_error_it->first;
                    double error_range = max_error - min_error;

                    if (error_range > epsilon) {
                        // Create vectors for parallel processing
                        std::vector<double> error_weights(filtered_models.size());
                        std::vector<double> kernel_weights(filtered_models.size());
                        std::vector<double> combined_weights(filtered_models.size());
                        std::vector<double> predictions(filtered_models.size());

                        // Compute error weights
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       error_weights.begin(),
                                       [min_error, error_range](const auto& model_pair) {
                                           double normalized_error = (model_pair.first - min_error) / error_range;
                                           std::vector<double> weight(1);
                                           kernel_fn(&normalized_error, 1, weight.data());
                                           return weight[0];
                                       }
                            );

                        // Get kernel weights and predictions
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       kernel_weights.begin(),
                                       [i](const auto& model_pair) {
                                           int local_index = i - model_pair.second->x_min_index;
                                           return model_pair.second->w[local_index];
                                       }
                            );

                        // Combine weights
                        std::transform(std::execution::par_unseq,
                                       error_weights.begin(), error_weights.end(),
                                       kernel_weights.begin(),
                                       combined_weights.begin(),
                                       std::multiplies<>()
                            );

                        // Get predictions
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       predictions.begin(),
                                       [i](const auto& model_pair) {
                                           int local_index = i - model_pair.second->x_min_index;
                                           return model_pair.second->predictions[local_index];
                                       }
                            );

                        // Compute weighted sum using transform_reduce
                        weighted_sum = std::transform_reduce(
                            std::execution::par_unseq,
                            predictions.begin(), predictions.end(),
                            combined_weights.begin(),
                            0.0,
                            std::plus<>(),
                            std::multiplies<>()
                            );

                        weight_sum = std::reduce(
                            std::execution::par_unseq,
                            combined_weights.begin(), combined_weights.end()
                            );

                        if (weight_sum > epsilon) {
                            results.predictions[i] = weighted_sum / weight_sum;
                            results.errors[i] = std::abs(results.predictions[i] - y[i]);
                        } else {
                            // Fallback to best model if combined weights sum to zero
                            auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                               [i](const auto& a, const auto& b) {
                                                                   return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                               });
                            results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                            results.errors[i] = best_model->errors[i - best_model->x_min_index];
                        }
                    } else {
                        // When error range is too small, use only kernel weights
                        std::vector<double> kernel_weights(filtered_models.size());
                        std::vector<double> predictions(filtered_models.size());

                        // Get kernel weights and predictions in parallel
                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       kernel_weights.begin(),
                                       [i](const auto& model_pair) {
                                           int local_index = i - model_pair.second->x_min_index;
                                           return model_pair.second->w[local_index];
                                       }
                            );

                        std::transform(std::execution::par_unseq,
                                       filtered_models.begin(), filtered_models.end(),
                                       predictions.begin(),
                                       [i](const auto& model_pair) {
                                           int local_index = i - model_pair.second->x_min_index;
                                           return model_pair.second->predictions[local_index];
                                       }
                            );

                        weighted_sum = std::transform_reduce(
                            std::execution::par_unseq,
                            predictions.begin(), predictions.end(),
                            kernel_weights.begin(),
                            0.0,
                            std::plus<>(),
                            std::multiplies<>()
                            );

                        weight_sum = std::reduce(
                            std::execution::par_unseq,
                            kernel_weights.begin(), kernel_weights.end()
                            );

                        if (weight_sum > epsilon) {
                            results.predictions[i] = weighted_sum / weight_sum;
                            results.errors[i] = std::abs(results.predictions[i] - y[i]);
                        } else {
                            auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                               [i](const auto& a, const auto& b) {
                                                                   return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                               });
                            results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                            results.errors[i] = best_model->errors[i - best_model->x_min_index];
                        }
                    }
                } else {
                    // If no models passed filtering, use the best available model
                    auto best_model = std::min_element(pt_models[i].begin(), pt_models[i].end(),
                                                       [i](const auto& a, const auto& b) {
                                                           return a.errors[i - a.x_min_index] < b.errors[i - b.x_min_index];
                                                       });
                    results.predictions[i] = best_model->predictions[i - best_model->x_min_index];
                    results.errors[i] = best_model->errors[i - best_model->x_min_index];
                }
                break;
            }
            }
            // Thread-safe progress update
            if (verbose) {
                static std::atomic<size_t> progress{0};
                if (++progress % (n_points/100 + 1) == 0) {
                    avg_progress.update(progress);
                }
            }
        });

    if (verbose) {
        elapsed_time(avg_ptm, "\nModel averaging time: ");
        mem_tracker.report();
        elapsed_time(total_ptm, "\nTotal parallel MABILO computation time: ");
        Rprintf("Final ");
        mem_tracker.report();
    }

    return results;
}
#endif
