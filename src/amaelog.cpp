#include <R.h>
#include <Rinternals.h>

#undef length  // to resolve naming conflict between the R macro length defined in Rinternals.h and a member function in the C++ standard library's codecvt class
#undef eval

#include <execution>
#include <mutex>
#include <ANN/ANN.h>
#include <vector>
#include <numeric>    // for std::iota
#include <random>     // for std::mt19937
#include <algorithm>  // for std::shuffle
//#include <fstream>       // for debugging and elapsed.time

#include <Eigen/Dense>

#include "ulogit.hpp"
#include "kernels.h"
#include "maelog.hpp"
// #include "cpp_utils.hpp" // for debugging and elapsed.time

extern "C" {
    SEXP S_adaptive_maelog(
        SEXP x_r,
        SEXP y_r,
        SEXP max_adapt_iterations_r,
        SEXP convergence_threshold_r,
        SEXP c_min_r,
        SEXP c_max_r,
        SEXP power_r,
        SEXP kernel_type_r,
        SEXP min_points_r,
        SEXP max_iterations_r,
        SEXP ridge_lambda_r,
        SEXP tolerance_r);
}

struct adaptive_bandwidth_t {

    std::vector<double> predictions;
    std::vector<std::vector<double>> apredictions; // adaptive predictions

    double baseline_bw;                   // Local bandwidth parameter
    std::vector<double> var_bw;

    std::vector<double> beta1s;       ///< beta[1] coefficient of the local linear or quadratic model
    std::vector<double> beta2s;       ///< beta[2] coefficient of the local quadratic model

    // Diagnostic information
    std::vector<int> n_points_used;  // Number of points used in each local fit
    std::vector<bool> used_knn;      // Whether k-NN fallback was used

    // Input parameters
    bool fit_quadratic;      ///< Whether quadratic term was included in local models
    int kernel_type;         ///< Type of kernel function used for local weighting
    double ridge_lambda;     ///< Ridge regularization parameter
    double tolerance;
};

adaptive_bandwidth_t adaptive_maelog(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int max_adapt_iterations = 10,
    double convergence_threshold = 1e-4,
    double c_min = 0.2,
    double c_max = 5.0,
    double power = -0.2,
    int kernel_type = 7,
    int min_points = 6,
    int max_iterations = 100,
    double ridge_lambda = 1e-6,
    double tolerance = 1e-8);

struct vbw_local_logit_t {
    std::vector<double> preds;

    std::vector<double> beta1s;
    std::vector<double> beta2s;

    std::vector<double> brier_errors;

    std::vector<int> n_points_used;  // Number of points used in each local fit
    std::vector<bool> used_knn;      // Whether k-NN fallback was used

};

struct vbw_maelog_t {
    std::vector<double> predictions;
    std::vector<double> beta1s;                           ///< beta[1] coefficient of the local linear or quadratic model
    std::vector<double> beta2s;                           ///< beta[2] coefficient of the local quadratic model

    // Diagnostic information
    std::vector<int> n_points_used;  // Number of points used in each local fit
    std::vector<bool> used_knn;      // Whether k-NN fallback was used

    // Input parameters
    bool fit_quadratic;      ///< Whether quadratic term was included in local models
    std::vector<double> variable_bandwidth;  ///< Fixed bandwidth if > 0, otherwise bandwidth is selected by CV
    int kernel_type;         ///< Type of kernel function used for local weighting
    double ridge_lambda;     ///< Ridge regularization parameter
    double tolerance;
};

/**
 * @brief Performs local logistic regression with variable bandwidths
 *
 * @details This function implements a local logistic regression algorithm using variable
 * bandwidths for each point in the dataset. The algorithm fits local linear or quadratic
 * logistic models in the neighborhood of each data point, with the neighborhood size
 * determined by the provided variable bandwidths.
 *
 * The function uses a k-NN fallback mechanism when the bandwidth-based neighborhood
 * contains fewer than the minimum required points. This ensures numerical stability
 * and reliable estimates even in sparse regions of the data.
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param fit_quadratic If true, includes quadratic terms in local models
 * @param variable_bandwidth Vector of bandwidths, one for each data point
 * @param kernel_type Integer specifying kernel function (1=Gaussian, 2=Epanechnikov)
 * @param min_points Minimum number of points required in local neighborhood
 * @param max_iterations Maximum iterations for logistic regression fitting
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for model fitting
 *
 * @return vbw_maelog_t structure containing:
 *   - predictions: Fitted probabilities for each point
 *   - beta1s: Linear coefficients from local fits
 *   - beta2s: Quadratic coefficients (if fit_quadratic is true)
 *   - n_points_used: Number of points used in each local fit
 *   - used_knn: Boolean indicators for whether k-NN fallback was used
 *   - Input parameters are also stored in the return structure
 *
 * @throws R error if x and variable_bandwidth have different sizes
 *
 * @note This function is part of an R package and uses the error() function
 * for error handling instead of C++ exceptions
 */
vbw_maelog_t vbw_maelog(
    const std::vector<double>& x,
    const std::vector<double>& y,
    bool fit_quadratic,
    const std::vector<double>& variable_bandwidth,
    int kernel_type,
    int min_points,
    int max_iterations,
    double ridge_lambda,
    double tolerance) {

    // Consider adding checks like:
    if (x.size() != variable_bandwidth.size()) {
        error("x and variable_bandwidth must have same size");
    }

    vbw_maelog_t result;
    result.fit_quadratic = fit_quadratic;
    result.variable_bandwidth = variable_bandwidth;
    result.kernel_type = kernel_type;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    // Lambda function for fitting local logistic regression with a specific bandwidth
    auto fit_local_logistic = [](
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& var_bw,
        bool fit_quadratic,
        int kernel_type,
        int min_points,
        int max_iterations,
        double ridge_lambda,
        double tolerance) {

        // Create data points array for ANN
        int n_pts = x.size();
        ANNpointArray data_points = annAllocPts(n_pts, 1);  // 1 dimension
        for(int i = 0; i < n_pts; i++) {
            data_points[i][0] = x[i];
        }

        // Build kd-tree
        ANNkd_tree* kdtree = new ANNkd_tree(data_points, n_pts, 1);

        std::vector<std::vector<point_t>> pt_pred_weights(n_pts);
        std::vector<double> predictions(n_pts);
        std::vector<int> n_points_used(n_pts);  // Number of points used in each local fit
        std::vector<bool> used_knn(n_pts);      // Whether k-NN fallback was used

        std::vector<double> brier_errors;
        std::vector<double> beta1s(n_pts);
        std::vector<double> beta2s;
        if (fit_quadratic) beta2s.resize(n_pts);

        initialize_kernel(kernel_type, 1.0);

        // Pre-allocate arrays for ANN search
        ANNpoint query_point = annAllocPt(1);
        ANNidxArray nn_idx = new ANNidx[min_points];
        ANNdistArray nn_dists = new ANNdist[min_points];

        // Pre-allocate vectors
        std::vector<double> local_x, local_y, local_w, local_d;
        local_x.reserve(n_pts);
        local_y.reserve(n_pts);
        local_w.reserve(n_pts);
        local_d.reserve(n_pts);
        std::vector<int> local_indices; // Store indices along with local data for accurate mapping
        local_indices.reserve(n_pts);

        bool with_errors = false;

        // fitting individual models in a neighborhood of bandwidth 'bandwidth' of each point
        for (int pti = 0; pti < n_pts; pti++) {
            double center = x[pti];

            // Clear vectors for reuse
            local_x.clear();
            local_y.clear();
            local_w.clear();
            local_d.clear();
            local_indices.clear();

            // First try bandwidth-based neighborhood
            for (int i = 0; i < n_pts; ++i) {
                double dist = std::abs(x[i] - center) / var_bw[i];
                if (dist < 1.0) {
                    local_indices.push_back(i);
                }
            }

            int n_window_pts = local_indices.size();
            if (n_window_pts < min_points) { // If not enough points, switch to k-NN

                used_knn[pti] = true;
                n_points_used[pti] = min_points;

                query_point[0] = center;
                kdtree->annkSearch(query_point, min_points, nn_idx, nn_dists);

                // Calculate shifted x values and distances
                for (int i = 0; i < min_points; ++i) {
                    double shifted_x = x[nn_idx[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[nn_idx[i]]);
                    local_d.push_back(std::abs(shifted_x) / var_bw[nn_idx[i]]);
                }

                // Find maximum distance and normalize
                double max_dist = *std::max_element(local_d.begin(), local_d.end());
                if (max_dist < std::numeric_limits<double>::epsilon()) {
                    max_dist = std::numeric_limits<double>::epsilon();
                }

                for (auto& d : local_d) {
                    d /= max_dist;
                }

                n_window_pts = min_points;

            } else {
                n_points_used[pti] = n_window_pts;

                for (int i = 0; i < n_window_pts; ++i) {
                    double shifted_x = x[local_indices[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[local_indices[i]]);
                    local_d.push_back(std::abs(shifted_x) / var_bw[i]);
                }
            }

            local_w.resize(n_window_pts); // Resize local_w before using it
            kernel_fn(local_d.data(), n_window_pts, local_w.data());

            // Normalize weights with protection against zero sum
            double sum_weights = std::accumulate(local_w.begin(), local_w.end(), 0.0);
            if (sum_weights < std::numeric_limits<double>::epsilon()) {
                sum_weights = std::numeric_limits<double>::epsilon();
            }
            for (auto& w : local_w) {
                w /= sum_weights;
            }

            eigen_ulogit_t fit_result = eigen_ulogit_fit(
                local_x.data(),
                local_y.data(),
                local_w,
                fit_quadratic,
                max_iterations,
                ridge_lambda,
                tolerance,
                with_errors
            );

            // Store predictions and errors
            point_t pt;
            for (size_t i = 0; i < local_x.size(); ++i) {
                int orig_idx = local_indices[i];
                pt.w = local_w[i];
                pt.p = fit_result.predictions[i];
                pt.beta1 = fit_result.beta[1];
                if (fit_quadratic) pt.beta2 = fit_result.beta[2];
                pt_pred_weights[orig_idx].push_back(pt);
            }
        }

        // Compute weighted averages
        for (int pti = 0; pti < n_pts; pti++) {
            const auto& v = pt_pred_weights[pti];
            if (v.empty()) {
                predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double total_weight = 0.0;
            double weighted_prediction = 0.0;
            double weighted_beta1 = 0.0;
            double weighted_beta2 = 0.0;

            for (const auto& pt : v) {
                weighted_prediction += pt.w * pt.p;
                weighted_beta1 += pt.w * pt.beta1;
                if (fit_quadratic) weighted_beta2 += pt.w * pt.beta2;
                total_weight += pt.w;
            }

            if (total_weight > 0.0) {
                predictions[pti] = weighted_prediction / total_weight;
                beta1s[pti] = weighted_beta1 / total_weight;
                if (fit_quadratic) beta2s[pti] = weighted_beta2 / total_weight;
            } else {
                predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                beta1s[pti] = std::numeric_limits<double>::quiet_NaN();
                if (fit_quadratic) beta2s[pti] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // Clean up ANN allocations
        annDeallocPt(query_point);
        delete[] nn_idx;
        delete[] nn_dists;
        // Clean up ANN allocations
        annDeallocPts(data_points);
        delete kdtree;
        annClose();

        vbw_local_logit_t ll;
        ll.preds = std::move(predictions);
        ll.beta1s = std::move(beta1s);
        if (fit_quadratic) ll.beta2s = std::move(beta2s);
        ll.n_points_used = std::move(n_points_used);
        ll.used_knn = std::move(used_knn);

        return ll;
    };


    auto ll = fit_local_logistic(x, y, variable_bandwidth, fit_quadratic, kernel_type,
                                 min_points, max_iterations, ridge_lambda, tolerance);

    result.beta1s = std::move(ll.beta1s);
    if (fit_quadratic) {
        result.beta2s = std::move(ll.beta2s);
    }

    result.predictions = std::move(ll.preds);
    result.n_points_used = std::move(ll.n_points_used);
    result.used_knn = std::move(ll.used_knn);

    return result;
}

/**
 * @brief Performs adaptive bandwidth local logistic regression
 *
 * @details This function implements an adaptive bandwidth selection strategy for local
 * logistic regression. It first estimates a global pilot bandwidth using cross-validation,
 * then adapts the bandwidth locally based on the complexity of the fitted surface
 * (measured by the quadratic coefficient beta2).
 *
 * The adaptive bandwidth formula is:
 * h(x) = h_0 * min{max{(|beta_2(x)|/tilde{beta_2})^(-0.2), c_min}, c_max}
 *
 * where:
 * - h_0 is the pilot bandwidth selected by cross-validation
 * - beta_2(x) is the local quadratic coefficient
 * - tilde{beta_2} is the median of |beta_2(x)|
 * - c_min, c_max are boundary factors to prevent extreme adaptation
 *
 * @param x Vector of predictor variables
 * @param y Vector of binary response variables (0 or 1)
 * @param c_min Minimum allowed adaptation factor
 * @param c_max Maximum allowed adaptation factor
 * @param power Power for adaptation formula
 * @param kernel_type Integer specifying kernel function (1=Gaussian, 2=Epanechnikov)
 * @param min_points Minimum number of points required in local neighborhood
 * @param max_iterations Maximum iterations for logistic regression fitting
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for model fitting
 *
 * @return adaptive_bandwidth_t structure containing:
 *   - baseline_bw: Initial pilot bandwidth
 *   - predictions: Initial predictions using pilot bandwidth
 *   - apredictions: Predictions using adaptive bandwidths
 *   - beta1s, beta2s: Coefficients from local fits
 *   - n_points_used: Number of points used in each local fit
 *   - used_knn: Boolean indicators for whether k-NN fallback was used
 *   - Input parameters are also stored in the return structure
 */
adaptive_bandwidth_t adaptive_maelog(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int max_adapt_iterations,
    double convergence_threshold,
    double c_min,
    double c_max,
    double power,
    int kernel_type,
    int min_points,
    int max_iterations,
    double ridge_lambda,
    double tolerance) {

    adaptive_bandwidth_t result;
    bool fit_quadratic = true;
    double pilot_bandwidth = -1.0;
    int cv_folds = 5;
    int n_bws = 50;
    double min_bw_factor = 0.05;
    double max_bw_factor = 0.9;
    bool with_errors = false;
    bool with_bw_predictions = false;

    // Get initial fit with CV-selected bandwidth
    maelog_t res = maelog(
        x,
        y,
        fit_quadratic,
        pilot_bandwidth,
        kernel_type,
        min_points,
        cv_folds,
        n_bws,
        min_bw_factor,
        max_bw_factor,
        max_iterations,
        ridge_lambda,
        tolerance,
        with_errors,
        with_bw_predictions);

    double baseline_bw = res.candidate_bandwidths[res.opt_brier_bw_idx];
    result.baseline_bw = baseline_bw;
    result.predictions = std::move(res.bw_predictions[res.opt_brier_bw_idx]);  // Store initial predictions

    int n_pts = x.size();
    std::vector<double> var_bw(n_pts);
    vbw_maelog_t prev_vbw_res;
    bool first_iteration = true;

    for(int iter = 0; iter < max_adapt_iterations; iter++) {
        // Compute median of absolute beta2 values from previous fit
        std::vector<double> abs_beta2;
        abs_beta2.reserve(n_pts);

        const std::vector<double>& beta2s_to_use = first_iteration ?
            res.beta2s : prev_vbw_res.beta2s;

        for (double b2 : beta2s_to_use) {
            abs_beta2.push_back(std::abs(b2));
        }

        size_t n = abs_beta2.size() / 2;
        std::nth_element(abs_beta2.begin(), abs_beta2.begin() + n, abs_beta2.end());
        double median_abs_beta2 = abs_beta2[n];

        // Store previous bandwidths for convergence check
        std::vector<double> prev_var_bw;
        if (!first_iteration) {
            prev_var_bw = var_bw;
        }

        // Compute adaptive bandwidths
        for (int i = 0; i < n_pts; i++) {
            double rel_complexity = std::abs(beta2s_to_use[i]) / median_abs_beta2;
            if (rel_complexity < std::numeric_limits<double>::epsilon()) {
                rel_complexity = std::numeric_limits<double>::epsilon();
            }
            double adaptation_factor = std::pow(rel_complexity, power);
            adaptation_factor = std::min(std::max(adaptation_factor, c_min), c_max);
            var_bw[i] = baseline_bw * adaptation_factor;
        }

        // Fit model with current adaptive bandwidths
        vbw_maelog_t vbw_res = vbw_maelog(
            x,
            y,
            fit_quadratic,
            var_bw,
            kernel_type,
            min_points,
            max_iterations,
            ridge_lambda,
            tolerance);

        // Store predictions for this iteration
        result.apredictions.push_back(vbw_res.predictions);

        // Check for convergence if not first iteration
        if (!first_iteration) {
            double max_bw_change = 0.0;
            for (size_t i = 0; i < n_pts; i++) {
                double rel_change = std::abs(var_bw[i] - prev_var_bw[i]) / prev_var_bw[i];
                max_bw_change = std::max(max_bw_change, rel_change);
            }

            if (max_bw_change < convergence_threshold) {
                break;  // Convergence achieved
            }
        }

        prev_vbw_res = vbw_res;
        first_iteration = false;
    }

    // Store final results
    result.beta1s = prev_vbw_res.beta1s;
    result.beta2s = prev_vbw_res.beta2s;
    result.n_points_used = prev_vbw_res.n_points_used;
    result.used_knn = prev_vbw_res.used_knn;
    result.fit_quadratic = fit_quadratic;
    result.var_bw = var_bw;
    result.kernel_type = kernel_type;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    return result;
}



/**
 * @brief R interface for adaptive bandwidth local logistic regression
 *
 * @details This function provides an R interface to the C++ adaptive_maelog function.
 * It converts R objects to C++ types, calls the core implementation, and converts
 * the results back to R objects.
 *
 * @param x_r SEXP containing numeric vector of predictor variables
 * @param y_r SEXP containing numeric vector of binary response variables (0 or 1)
 * @param c_min_r SEXP containing numeric minimum allowed adaptation factor
 * @param c_max_r SEXP containing numeric maximum allowed adaptation factor
 * @param power_r SEXP containing numeric power for adaptation formula
 * @param kernel_type_r SEXP containing integer kernel function type
 * @param min_points_r SEXP containing integer minimum points for local fits
 * @param max_iterations_r SEXP containing integer maximum iterations for fitting
 * @param ridge_lambda_r SEXP containing numeric ridge regularization parameter
 * @param tolerance_r SEXP containing numeric convergence tolerance
 *
 * @return SEXP (list) containing:
 *   - predictions: Initial predictions using pilot bandwidth
 *   - apredictions: Matrix of predictions using adaptive bandwidths
 *   - baseline_bw: Initial pilot bandwidth
 *   - var_bw: Vector of variable bandwidths
 *   - beta1s: Linear coefficients
 *   - beta2s: Quadratic coefficients
 *   - n_points_used: Number of points used in each local fit
 *   - used_knn: Indicators for k-NN fallback usage
 *   - fit_info: List of fitting parameters and information
 */
SEXP S_adaptive_maelog(
    SEXP x_r,
    SEXP y_r,
    SEXP max_adapt_iterations_r,
    SEXP convergence_threshold_r,
    SEXP c_min_r,
    SEXP c_max_r,
    SEXP power_r,
    SEXP kernel_type_r,
    SEXP min_points_r,
    SEXP max_iterations_r,
    SEXP ridge_lambda_r,
    SEXP tolerance_r
    ) {

    int n_protected = 0;

    try {
        // Check for NULL inputs
        if (x_r == R_NilValue ||
            y_r == R_NilValue ||
            kernel_type_r == R_NilValue ||
            min_points_r == R_NilValue ||
            max_iterations_r == R_NilValue ||
            ridge_lambda_r == R_NilValue ||
            tolerance_r == R_NilValue) {
            Rf_error("Input arguments cannot be NULL");
        }

        int n_points = LENGTH(x_r);

        // Convert inputs
        std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
        std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);
        double max_adapt_iterations = REAL(max_adapt_iterations_r)[0];
        double convergence_threshold = REAL(convergence_threshold_r)[0];
        double c_min = REAL(c_min_r)[0];
        double c_max = REAL(c_max_r)[0];
        double power = REAL(power_r)[0];
        int kernel_type = INTEGER(kernel_type_r)[0];
        int min_points = INTEGER(min_points_r)[0];
        int max_iterations = INTEGER(max_iterations_r)[0];
        double ridge_lambda = REAL(ridge_lambda_r)[0];
        double tolerance = REAL(tolerance_r)[0];

        auto result = adaptive_maelog(x,
                                         y,
                                         max_adapt_iterations,
                                         convergence_threshold,
                                         c_min,
                                         c_max,
                                         power,
                                         kernel_type,
                                         min_points,
                                         max_iterations,
                                         ridge_lambda,
                                         tolerance);

        // Create return list with updated size
        const int N_RETURN_COMPS = 9;
        SEXP result_list = PROTECT(Rf_allocVector(VECSXP, N_RETURN_COMPS)); n_protected++;

        // Updated names for list elements
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_RETURN_COMPS)); n_protected++;
        SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 1, Rf_mkChar("apredictions"));
        SET_STRING_ELT(names, 2, Rf_mkChar("baseline_bw"));
        SET_STRING_ELT(names, 3, Rf_mkChar("var_bw"));
        SET_STRING_ELT(names, 4, Rf_mkChar("beta1s"));
        SET_STRING_ELT(names, 5, Rf_mkChar("beta2s"));
        SET_STRING_ELT(names, 6, Rf_mkChar("n_points_used"));
        SET_STRING_ELT(names, 7, Rf_mkChar("used_knn"));
        SET_STRING_ELT(names, 8, Rf_mkChar("fit_info"));

        // predictions
        SEXP pred_r = PROTECT(Rf_allocVector(REALSXP, result.predictions.size())); n_protected++;
        std::copy(result.predictions.begin(), result.predictions.end(), REAL(pred_r));
        SET_VECTOR_ELT(result_list, 0, pred_r);

        // Convert matrix of adaptive predictions
        SEXP apred_r = PROTECT(Rf_allocMatrix(REALSXP,
                                              result.apredictions[0].size(),
                                              result.apredictions.size())); n_protected++;
        for(size_t i = 0; i < result.apredictions.size(); i++) {
            std::copy(result.apredictions[i].begin(),
                      result.apredictions[i].end(),
                      REAL(apred_r) + i * result.apredictions[0].size());
        }
        SET_VECTOR_ELT(result_list, 1, apred_r);

        // baseline bandwidth
        SEXP baseline_bw_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(baseline_bw_r)[0] = result.baseline_bw;
        SET_VECTOR_ELT(result_list, 2, baseline_bw_r);

        // var_bw
        SEXP var_bw_r = PROTECT(Rf_allocVector(REALSXP, result.var_bw.size())); n_protected++;
        std::copy(result.var_bw.begin(), result.var_bw.end(),
                 REAL(var_bw_r));
        SET_VECTOR_ELT(result_list, 3, var_bw_r);

        // Beta coefficients
        SEXP beta1s_r;
        beta1s_r = PROTECT(Rf_allocVector(REALSXP, result.beta1s.size())); n_protected++;
        std::copy(result.beta1s.begin(), result.beta1s.end(), REAL(beta1s_r));
        SET_VECTOR_ELT(result_list, 4, beta1s_r);

        SEXP beta2s_r;
        beta2s_r = PROTECT(Rf_allocVector(REALSXP, result.beta2s.size())); n_protected++;
        SET_VECTOR_ELT(result_list, 5, beta2s_r);

        // n_points_used and used_knn
        SEXP n_points_used_r = PROTECT(Rf_allocVector(INTSXP, result.n_points_used.size())); n_protected++;
        std::copy(result.n_points_used.begin(), result.n_points_used.end(), INTEGER(n_points_used_r));
        SET_VECTOR_ELT(result_list, 6, n_points_used_r);

        SEXP used_knn_r = PROTECT(Rf_allocVector(LGLSXP, result.used_knn.size())); n_protected++;
        std::copy(result.used_knn.begin(), result.used_knn.end(), LOGICAL(used_knn_r));
        SET_VECTOR_ELT(result_list, 7, used_knn_r);


        SEXP fit_info_names = PROTECT(Rf_allocVector(STRSXP, 8)); n_protected++;
        SET_STRING_ELT(fit_info_names, 0, Rf_mkChar("c_min"));
        SET_STRING_ELT(fit_info_names, 1, Rf_mkChar("c_max"));
        SET_STRING_ELT(fit_info_names, 2, Rf_mkChar("power"));
        SET_STRING_ELT(fit_info_names, 3, Rf_mkChar("kernel_type"));
        SET_STRING_ELT(fit_info_names, 4, Rf_mkChar("min_points"));
        SET_STRING_ELT(fit_info_names, 5, Rf_mkChar("max_iterations"));
        SET_STRING_ELT(fit_info_names, 6, Rf_mkChar("ridge_lambda"));
        SET_STRING_ELT(fit_info_names, 7, Rf_mkChar("tolerance"));

        // Fit info
        SEXP fit_info = PROTECT(Rf_allocVector(VECSXP, 8)); n_protected++;

        // Create new SEXPs for each parameter
        SEXP c_min_out_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(c_min_out_r)[0] = c_min;
        SET_VECTOR_ELT(fit_info, 0, c_min_out_r);

        SEXP c_max_out_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(c_max_out_r)[0] = c_max;
        SET_VECTOR_ELT(fit_info, 1, c_max_out_r);

        SEXP power_out_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(power_out_r)[0] = power;
        SET_VECTOR_ELT(fit_info, 2, power_out_r);

        SEXP kernel_type_out_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
        INTEGER(kernel_type_out_r)[0] = kernel_type;
        SET_VECTOR_ELT(fit_info, 3, kernel_type_out_r);

        SEXP min_points_out_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
        INTEGER(min_points_out_r)[0] = min_points;
        SET_VECTOR_ELT(fit_info, 4, min_points_out_r);

        SEXP max_iterations_out_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
        INTEGER(max_iterations_out_r)[0] = max_iterations;
        SET_VECTOR_ELT(fit_info, 5, max_iterations_out_r);

        SEXP ridge_lambda_out_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(ridge_lambda_out_r)[0] = ridge_lambda;
        SET_VECTOR_ELT(fit_info, 6, ridge_lambda_out_r);

        SEXP tolerance_out_r = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
        REAL(tolerance_out_r)[0] = tolerance;
        SET_VECTOR_ELT(fit_info, 7, tolerance_out_r);

        // Set names for fit_info
        Rf_setAttrib(fit_info, R_NamesSymbol, fit_info_names);

        // Set names for the main list
        Rf_setAttrib(result_list, R_NamesSymbol, names);

        UNPROTECT(n_protected);
        return result_list;
    }
    catch (const std::exception& e) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("C++ error in maelog: %s", e.what());
    }
    catch (...) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("Unknown error in maelog");
    }
}
