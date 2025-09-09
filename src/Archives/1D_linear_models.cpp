
lm_loocv_t predict_lm_1d_loocv(const std::vector<double>& y,
                               const std::vector<double>& x,
                               const std::vector<double>& w,
                               const std::vector<int>& vertex_indices,
                               int ref_index,
                               double epsilon = 1e-8) {
    int n_points = x.size();

    // Modified input validation
    // the commented out condition is checked in the parent function
    // if (n_points != y.size() || n_points != w.size() || n_points != vertex_indices.size()) {
    //     Rf_error("All vectors must have the same length");
    // }
    if (ref_index < -1 || ref_index >= n_points) {
        REPORT_ERROR("Reference index out of bounds");
    }

    // Store vertex indices, x values, and weights
    lm_loocv_t results;
    results.vertex_indices = vertex_indices;
    results.x_values = x;
    results.w_values = w;  // Store the weight values
    results.ref_index = ref_index;

    #define DEBUG__predict_lm_1d_loocv 0
    #if DEBUG__predict_lm_1d_loocv
    Rprintf("\nIn predict_lm_1d_loocv()\n");
    Rprintf("n_points: %d\n", n_points);
    Rprintf("ref_index: %d\n", ref_index);
    print_vect(vertex_indices, "vertex_indices");
    print_vect(results.vertex_indices, "results.vertex_indices");
    print_vect(results.x_values, "results.x_values");
    print_vect(results.w_values, "results.w_values");
    Rprintf("results.ref_index: %d\n", results.ref_index);
    #endif

    // Create working copy of x for computations
    std::vector<double> x_working = x;

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        //if (weight < 0) Rf_error("Negative weights are not allowed");
        total_weight += weight;
    }
    if (total_weight <= 0) Rf_error("Sum of weights must be positive");

    // Store reference point
    results.x_ref = x[ref_index];

    // Center x_working around reference point
    for (auto& xi : x_working) xi -= results.x_ref;

    // Calculate weighted means
    results.x_wmean = 0.0;
    results.y_wmean = 0.0;
    for (int i = 0; i < n_points; ++i) {
        results.x_wmean += w[i] * x_working[i] / total_weight;
        results.y_wmean += w[i] * y[i] / total_weight;
    }

    // Center x_working around weighted mean for leverage calculation
    double sum_wx_squared = 0.0;
    for (int i = 0; i < n_points; ++i) {
        x_working[i] -= results.x_wmean;  // Now x is centered around its weighted mean
        sum_wx_squared += w[i] * x_working[i] * x_working[i];
    }

    // Calculate slope if x has sufficient variation
    results.slope = 0.0;
    if (sum_wx_squared > epsilon) {
        double wxy_sum = 0.0;
        for (int i = 0; i < n_points; ++i) {
            wxy_sum += w[i] * x_working[i] * y[i];
        }
        results.slope = wxy_sum / sum_wx_squared;
    }

#if DEBUG__predict_lm_1d_loocv
    Rprintf("Just before results.predicted_value = results.predict(vertex_indices[ref_index])\n");
    print_vect(vertex_indices, "vertex_indices");
    Rprintf("ref_index: %d\n", ref_index);
    Rprintf("vertex_indices[ref_index]: %d\n", vertex_indices[ref_index]);
#endif

    // Calculate LOOCV components
    double loocv_sum = 0.0;
    for (int i = 0; i < n_points; ++i) {
        // Calculate weighted leverage h_i
        double h_i = w[i] * (1.0/total_weight + (x_working[i] * x_working[i]) / sum_wx_squared);

        // Calculate fitted value using the predict method

#if DEBUG__predict_lm_1d_loocv
        Rprintf("Just before y_hat = results.predict(vertex_indices[i])\n");
        print_vect(vertex_indices, "vertex_indices");
        Rprintf("i: %d\n", i);
        Rprintf("vertex_indices[i]: %d\n", vertex_indices[i]);
#endif

        double y_hat = results.predict(vertex_indices[i]);

        // Calculate LOOCV prediction error
        double residual = (y[i] - y_hat) / (1.0 - h_i);
        double squared_error = residual * residual;

        // Add to overall LOOCV sum
        loocv_sum += squared_error;

        // Store specific error for reference vertex
        if (i == ref_index) {
            results.loocv_at_ref_vertex = squared_error;
        }
    }

    // LOOCV MSE is the average of squared prediction errors

    // Handle ref_index = -1 case
    if (ref_index == -1) {
        results.predicted_value = std::numeric_limits<double>::quiet_NaN();
        results.loocv_at_ref_vertex = std::numeric_limits<double>::quiet_NaN();
    } else {
        // Calculate predicted value at reference point
        results.predicted_value = results.predict(vertex_indices[ref_index]);
    }

    return results;
}
