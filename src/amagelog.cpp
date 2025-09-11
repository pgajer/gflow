//
// Model average local logistic regression with models positioned over a uniform grid
//
#include <R.h>
#include <Rinternals.h>

#undef length  // to resolve naming conflict between the R macro length defined in Rinternals.h and a member function in the C++ standard library's codecvt class
#undef Rf_eval

#include <execution>
#include <mutex>
#include <vector>
#include <numeric>    // for std::iota
#include <random>     // for std::mt19937
#include <algorithm>  // for std::shuffle
#include <fstream>    // for debugging and elapsed.time
#include <ANN/ANN.h>

#include <Eigen/Dense>

#include "ulogit.hpp"
#include "kernels.h"
#include "maelog.hpp"
#include "cpp_utils.hpp" // for debugging and elapsed.time

extern "C" {
    SEXP S_amagelogit(
    SEXP x_r,
    SEXP y_r,
    SEXP grid_size_r,  // Add this parameter
    SEXP fit_quadratic_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP min_points_r,
    SEXP cv_folds_r,
    SEXP n_bws_r,
    SEXP min_bw_factor_r,
    SEXP max_bw_factor_r,
    SEXP max_iterations_r,
    SEXP ridge_lambda_r,
    SEXP tolerance_r,
    SEXP with_bw_predictions_r);
}


/**
 * @brief Performs local logistic regression with automatic bandwidth selection
 *
 * @details Fits local logistic regression models centered at points of a uniform grid spanning
 * the range of input x values. The models can be either linear or quadratic, with
 * bandwidth selection performed via K-fold cross-validation when pilot_bandwidth <= 0.
 *
 * The function uses
 *  Brier score: $$(\hat{p}_{(-i)} - y_i)^2$$
 * for bandwidth selection, where $$\hat{p}_{(-i)}$$ represents the leave-one-out prediction for observation i.
 *
 * The algorithm proceeds as follows:
 * 1. If pilot_bandwidth > 0, uses this fixed bandwidth for all local models
 * 2. Otherwise:
 *    a. Creates a grid of candidate bandwidths
 *    b. For each bandwidth:
 *       - Fits local models using either LOOCV or k-fold CV
 *       - Computes all three Rf_error measures
 *    c. Selects optimal bandwidths minimizing each Rf_error measure
 * 3. Returns predictions and model parameters for optimal or specified bandwidth(s)
 *
 * @param x Vector of predictor values
 * @param y Vector of binary response values (0 or 1)
 * @param grid_size Number of points in the uniform grid where models are centered
 * @param fit_quadratic Whether to include quadratic term in local models
 * @param pilot_bandwidth Fixed bandwidth if > 0, otherwise bandwidth is selected by CV
 * @param kernel_type Type of kernel function for local weighting (0:Gaussian, 1:Epanechnikov)
 * @param min_points Minimum number of points required for local fitting
 * @param cv_folds Number of cross-validation folds (0 for approximate LOOCV)
 * @param n_bws Number of bandwidths to test in CV
 * @param min_bw_factor Lower bound factor for bandwidth grid relative to range(x)
 * @param max_bw_factor Upper bound factor for bandwidth grid relative to range(x)
 * @param max_iterations Maximum iterations for logistic regression fitting
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for logistic regression
 * @param with_bw_grid_predictions Whether to return predictions for all tested bandwidths
 *
 * @return amagelogit_t structure containing:
 *         - grid_size: number of grid points
 *         - x_grid: uniform grid points where models are centered
 *         - predictions: interpolated predictions at original x points
 *         - bw_grid_predictions: predictions at grid points for each bandwidth
 *         - candidate_bandwidths: tested bandwidth values
 *         - mean_brier_errors: cross-validation errors for each bandwidth
 *         - opt_brier_bw_idx: index of optimal bandwidth
 *         - other fitting parameters
 *
 * Additional members store input parameters and configuration:
 * - fit_quadratic: Whether quadratic terms were included
 * - pilot_bandwidth: Fixed bandwidth if specified
 * - kernel_type: Type of kernel function used
 * - cv_folds: Number of CV folds used
 * - min_bw_factor: Lower bound factor for bandwidth grid
 * - max_bw_factor: Upper bound factor for bandwidth grid
 * - max_iterations: Maximum fitting iterations
 * - ridge_lambda: Ridge regularization parameter
 * - tolerance: Convergence tolerance
 *
 * @throws std::runtime_error if input x values are effectively constant
 *
 * @note If with_bw_predictions=false, only predictions for unique optimal
 * bandwidths are retained in bw_predictions to conserve memory.
 *
 * @Rf_warning Input vectors x and y must be the same length and y must contain
 * only binary values (0 or 1). The function uses the ANN library for efficient
 * nearest neighbor searches, which must be properly initialized before calling.
 *
 * @see eigen_ulogit_fit For details on the local logistic regression fitting
 * @see initialize_kernel For kernel function initialization
 *
 */

#if 0
struct logit_model_t {
   /** @brief Model coefficients vector
    *  @details Contains:
    *   - beta[0]: Intercept term
    *   - beta[1]: Linear term coefficient
    *   - beta[2]: Quadratic term coefficient (if fit_quadratic=true)
    */
    Eigen::VectorXd beta;
    //bool fit_quadratic; // needed for predict fn; do we actually need that function?

    int grid_start_idx; // grid index of the start grid point of the window of the model
    int grid_end_idx;   // grid index of the end grid point of the window of the defined
    //int ref_index; ///< Reference point within the window
    //double bw;     // bandwidth of the model
    std::vector<double> grid_weights; // weights of the grid points within

    #if 0
public:
    std::vector<double> predict(const std::vector<double>& x) const {
        // Check for empty input
        if (x.empty()) {
            return std::vector<double>();
        }

        // Verify beta has correct dimensions based on model type
        int expected_size = fit_quadratic ? 3 : 2;
        if (beta.size() != expected_size) {
            Rf_error("Model coefficients vector has incorrect size");
        }

        int n = x.size();
        std::vector<double> predictions(n);

        for (int i = 0; i < n; ++i) {
            double eta = beta(0) + beta(1) * x[i];
            if (fit_quadratic) {
                eta += beta(2) * x[i] * x[i];
            }
            predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }

        return predictions;
    }
    #endif
};

struct local_logit_t {
    std::vector<double> grid_predictions;
    std::vector<logit_model_t> grid_models;
};

struct amagelogit_t {
    // bandwidth grid
    std::vector<double> candidate_bandwidths;                   ///< Grid of bandwidths tested during optimization

    // Mean errors and optimal indices
    std::vector<double> mean_brier_errors;                      ///< Mean Brier Rf_error for each candidate bandwidth
    int opt_brier_bw_idx;                                       ///< Index of bandwidth with minimal mean Brier Rf_error

    // grid-based members
    std::vector<double> x_grid;                                ///< Uniform grid over the range of x values, models are estimated at these locations
    std::vector<std::vector<double>> bw_grid_predictions;      ///< Predictions for each bandwidth in LOOCV or CV estimation
    std::vector<std::vector<double>> bw_grid_errors;

    std::vector<double> predictions;                           ///< Predictions at x points

    // Input parameters
    bool fit_quadratic;      ///< Whether quadratic term was included in local models
    double pilot_bandwidth;  ///< Fixed bandwidth if > 0, otherwise bandwidth is selected by CV
    int kernel_type;         ///< Type of kernel function used for local weighting
    int cv_folds;            ///< Number of CV folds (0 for LOOCV approximation)
    double min_bw_factor;    ///< Lower bound factor for bandwidth grid relative to h_rot
    double max_bw_factor;    ///< Upper bound factor for bandwidth grid relative to h_rot
    int max_iterations;      ///< Maximum iterations for logistic regression fitting
    double ridge_lambda;     ///< Ridge regularization parameter
    double tolerance;        ///< Number of points in the evaluation grid
};


amagelogit_t amagelogit(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int grid_size,
    bool fit_quadratic,
    double pilot_bandwidth,
    int kernel_type,
    int min_points,
    int cv_folds,
    int n_bws,
    double min_bw_factor,
    double max_bw_factor,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_bw_grid_predictions) {

    if (grid_size < 2) {
        Rf_error("grid_size must be at least 2");
    }
    if (x.size() != y.size()) {
        Rf_error("x and y must have the same size");
    }
    if (x.empty()) {
        Rf_error("Input vectors cannot be empty");
    }
    if (min_points > static_cast<int>(x.size())) {
        Rf_error("min_points cannot be larger than the number of data points");
    }

    amagelogit_t result;
    result.fit_quadratic = fit_quadratic;
    result.pilot_bandwidth = pilot_bandwidth;
    result.kernel_type = kernel_type;
    result.cv_folds = cv_folds;
    result.min_bw_factor = min_bw_factor;
    result.max_bw_factor = max_bw_factor;
    result.max_iterations = max_iterations;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    int n_pts = x.size();
    result.predictions.resize(n_pts);

    // Create bandwidth grid
    result.candidate_bandwidths.resize(n_bws);
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    double x_range = x_max - x_min;

    if (x_range < std::numeric_limits<double>::epsilon()) {
        Rf_error("Input x values are effectively constant");
    }

    double min_bw = min_bw_factor * x_range;
    double max_bw = max_bw_factor * x_range;
    double dx = (max_bw - min_bw) / (n_bws - 1);
    for(int i = 0; i < n_bws; i++) {
        result.candidate_bandwidths[i] = min_bw + i * dx;
    }

    // Create x_grid
    result.x_grid.resize(grid_size);
    double grid_dx = x_range / (grid_size - 1);
    for(int i = 0; i < grid_size; i++) {
        result.x_grid[i] = x_min + i * grid_dx;
    }

    std::vector<std::vector<logit_model_t>> bw_grid_models; //  bw_grid_models[bw_idx] stores model parameters for every model centered the grid points

    // Lambda function for fitting local logistic regression with a specific bandwidth - models are estimated only at the x_grid points
    auto fit_local_logit = [&cv_folds,&result,&grid_size](
        const std::vector<double>& x,
        const std::vector<double>& y,
        double bandwidth,
        bool fit_quadratic,
        int kernel_type,
        int min_points,
        int max_iterations,
        double ridge_lambda,
        double tolerance,
        bool get_grid_models) {

        int n_pts = x.size();
        // Create data points array for ANN
        ANNpointArray data_points = annAllocPts(n_pts, 1);  // 1 dimension
        for(int i = 0; i < n_pts; i++) {
            data_points[i][0] = x[i];
        }
        // Maybe instead of using kdtree we should sort x and the identify (bi-)kNN's ??? <<<---
        // Build kd-tree
        ANNkd_tree* kdtree = new ANNkd_tree(data_points, n_pts, 1);

        // Instead of computing at x points, compute at grid points
        std::vector<double> grid_predictions(grid_size);
        std::vector<std::vector<point_t>> grid_pt_pred_weights(grid_size);
        for(auto& v : grid_pt_pred_weights) {
            v.reserve(grid_size); // Pre-allocate reasonable size
        }
        std::vector<logit_model_t> grid_models;
        if (get_grid_models) grid_models.resize(grid_size);

        initialize_kernel(kernel_type, 1.0);

        // Pre-allocate arrays for ANN search
        ANNpoint query_point = annAllocPt(1);
        ANNidxArray nn_idx = new ANNidx[min_points];
        ANNdistArray nn_dists = new ANNdist[min_points];

        // Pre-allocate local data vectors for model input
        std::vector<double> local_x, local_y, local_w, local_d;
        local_x.reserve(n_pts);
        local_y.reserve(n_pts);
        local_w.reserve(n_pts);
        local_d.reserve(n_pts);
        std::vector<int> local_indices; // Store indices along with local data for accurate mapping
        local_indices.reserve(n_pts);

        // Pre-allocate local x_grid vectors for model predictions at the grid locations
        std::vector<double> local_x_grid, local_w_grid, local_d_grid;
        local_x_grid.reserve(grid_size);
        local_w_grid.reserve(grid_size);
        local_d_grid.reserve(grid_size);
        std::vector<int> local_grid_indices; // Store indices along with local data for accurate mapping
        local_grid_indices.reserve(grid_size);

        // fitting individual models in a neighborhood of bandwidth 'bandwidth' of each point
        for (int pti = 0; pti < grid_size; pti++) {
            double center = result.x_grid[pti];

            // Clear vectors for reuse
            local_x.clear();
            local_y.clear();
            local_w.clear();
            local_d.clear();
            local_indices.clear();
            //
            local_x_grid.clear();
            local_w_grid.clear();
            local_d_grid.clear();
            local_grid_indices.clear();

            // First try bandwidth-based neighborhood
            for (int i = 0; i < n_pts; ++i) {
                double dist = std::abs(x[i] - center) / bandwidth;
                if (dist < 1.0) {
                    local_indices.push_back(i);
                }
            }
            for (int i = 0; i < grid_size; ++i) {
                double dist = std::abs(result.x_grid[i] - center) / bandwidth;
                if (dist < 1.0) {
                    local_grid_indices.push_back(i);
                }
            }

            int n_window_pts = local_indices.size();
            if (n_window_pts < min_points) { // If not enough points, switch to k-NN
                local_indices.clear();  // Clear existing indices
                query_point[0] = center;
                kdtree->annkSearch(query_point, min_points, nn_idx, nn_dists);

                // Calculate shifted x values and distances
                for (int i = 0; i < min_points; ++i) {
                    local_indices.push_back(nn_idx[i]);  // Store the new indices
                    double shifted_x = x[nn_idx[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[nn_idx[i]]);
                    local_d.push_back(std::abs(shifted_x));
                }

                // Find maximum distance and normalize
                double max_dist = *std::max_element(local_d.begin(), local_d.end());
                if (max_dist < std::numeric_limits<double>::epsilon()) {
                    max_dist = std::numeric_limits<double>::epsilon();
                }

                for (auto& d : local_d) d /= max_dist;

                n_window_pts = min_points;

                // updating local grid members
                local_grid_indices.clear();
                for (int i = 0; i < grid_size; ++i) {
                    double dist = std::abs(result.x_grid[i] - center);
                    if (dist < max_dist) {
                        local_grid_indices.push_back(i);
                    }
                }

                for (size_t i = 0; i < local_grid_indices.size(); ++i) {
                    double shifted_x_grid = result.x_grid[local_grid_indices[i]] - center;
                    local_x_grid.push_back(shifted_x_grid);
                    local_d_grid.push_back(std::abs(shifted_x_grid) / bandwidth);
                }

                // Normalize distances
                for (auto& d : local_d_grid) d /= max_dist;

            } else {
                for (int i = 0; i < n_window_pts; ++i) {
                    double shifted_x = x[local_indices[i]] - center;
                    local_x.push_back(shifted_x);
                    local_y.push_back(y[local_indices[i]]);
                    local_d.push_back(std::abs(shifted_x) / bandwidth);
                }

                for (size_t i = 0; i < local_grid_indices.size(); ++i) {
                    double shifted_x_grid = result.x_grid[local_grid_indices[i]] - center;
                    local_x_grid.push_back(shifted_x_grid);
                    local_d_grid.push_back(std::abs(shifted_x_grid) / bandwidth);
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

            bool with_errors = false;
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

            std::vector<double> local_p_grid = fit_result.predict(local_x_grid); // predicting conditional expectation values at the grid point within the support of the model

            // computing local_w_grid - needed for model avaraging
            local_w_grid.resize(local_x_grid.size()); // Resize
            kernel_fn(local_d_grid.data(), (int)local_x_grid.size(), local_w_grid.data());
            double sum_grid_weights = std::accumulate(local_w_grid.begin(), local_w_grid.end(), 0.0);
            if (sum_grid_weights < std::numeric_limits<double>::epsilon()) {
                sum_grid_weights = std::numeric_limits<double>::epsilon();
            }
            for (auto& w : local_w_grid) w /= sum_grid_weights;

            if (get_grid_models) {
                grid_models[pti].beta = fit_result.beta;
                grid_models[pti].grid_start_idx = local_grid_indices[0];
                grid_models[pti].grid_end_idx   = local_grid_indices[local_grid_indices.size() - 1];
                grid_models[pti].grid_weights   = local_w_grid;
            }

            // Store predictions and errors at grid points
            point_t pt;
            for (size_t i = 0; i < local_grid_indices.size(); ++i) {
                int orig_idx = local_grid_indices[i];
                pt.w = local_w_grid[i];
                pt.p = local_p_grid[i];
                //if (with_errors) pt.brier_error = fit_result.loocv_brier_errors[i];
                grid_pt_pred_weights[orig_idx].push_back(pt);
            }
        }

        // Compute weighted averages
        for (int pti = 0; pti < grid_size; pti++) {
            const auto& v = grid_pt_pred_weights[pti];
            if (v.empty()) {
                grid_predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double total_weight = 0.0;
            double weighted_prediction = 0.0;
            //double weighted_brier_error = 0.0;

            for (const auto& pt : v) {
                weighted_prediction += pt.w * pt.p;
                //if (with_errors) weighted_brier_error += pt.w * pt.brier_error;
                total_weight += pt.w;
            }

            if (total_weight > 0.0) {
                grid_predictions[pti] = weighted_prediction / total_weight;
                //if (with_errors) brier_errors[pti] = weighted_brier_error / total_weight;
            } else {
                grid_predictions[pti] = std::numeric_limits<double>::quiet_NaN();
                //if (with_errors) brier_errors[pti] = std::numeric_limits<double>::quiet_NaN();
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

        local_logit_t ll;
        ll.grid_predictions = grid_predictions;
        if (get_grid_models) ll.grid_models = grid_models;

        return ll;
    };

    auto interpolate_grid = [](double x, const std::vector<double>& x_grid,
                               const std::vector<double>& y_grid) -> double {
        if (x <= x_grid.front()) return y_grid.front();
        if (x >= x_grid.back()) return y_grid.back();

        auto upper = std::upper_bound(x_grid.begin(), x_grid.end(), x);
        size_t upper_idx = std::distance(x_grid.begin(), upper);
        size_t lower_idx = upper_idx - 1;

        double lambda = (x - x_grid[lower_idx]) / (x_grid[upper_idx] - x_grid[lower_idx]);
        return (1.0 - lambda) * y_grid[lower_idx] + lambda * y_grid[upper_idx];
    };

    // If pilot_bandwidth > 0, use it directly
    if (pilot_bandwidth > 0) {
        auto ll = fit_local_logit(x, y, pilot_bandwidth, fit_quadratic, kernel_type,
                                  min_points, max_iterations, ridge_lambda, tolerance, false);


        for(int i = 0; i < n_pts; i++)
            result.predictions[i] = interpolate_grid(x[i], result.x_grid, ll.grid_predictions);

        result.bw_grid_predictions.resize(1);
        result.bw_grid_predictions[0] = std::move(ll.grid_predictions);
        result.opt_brier_bw_idx = 0;

        return result;
    }


    // Initialize Rf_error vectors
    result.mean_brier_errors.resize(n_bws);
    result.bw_grid_predictions.resize(n_bws); // this vector is always initialized, even if with_bw_grid_predictions = false, in which case only optimal bw predictions will be non-empty

    result.bw_grid_errors.resize(n_bws);
    for(auto& v : result.bw_grid_errors) {
        v.resize(grid_size);
    }

    // Perform bandwidth selection using K-fold cross-validation
    int n = x.size();
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.end(), gen);

    int fold_size = n / cv_folds;
    for (int i = 0; i < n_bws; i++) {
        double total_brier_error    = 0.0;
        int valid_count = 0;

        for (int fold = 0; fold < cv_folds; fold++) {
            // Create training and validation sets
            std::vector<double> train_x, train_y, valid_x, valid_y, valid_idx;
            int start_idx = fold * fold_size;
            int end_idx = (fold == cv_folds - 1) ? n : (fold + 1) * fold_size;

            for (int j = 0; j < n; j++) {
                if (j >= start_idx && j < end_idx) {
                    valid_x.push_back(x[indices[j]]);
                    valid_y.push_back(y[indices[j]]);
                    valid_idx.push_back(indices[j]);
                } else {
                    train_x.push_back(x[indices[j]]);
                    train_y.push_back(y[indices[j]]);
                }
            }

            // Sort training data for efficient interpopolation
            std::vector<size_t> train_sort_idx(train_x.size());
            std::iota(train_sort_idx.begin(), train_sort_idx.end(), 0);
            std::sort(train_sort_idx.begin(), train_sort_idx.end(),
                      [&train_x](size_t i1, size_t i2) { return train_x[i1] < train_x[i2]; });

            std::vector<double> sorted_train_x, sorted_train_y;
            for (size_t idx : train_sort_idx) {
                sorted_train_x.push_back(train_x[idx]);
                sorted_train_y.push_back(train_y[idx]);
            }

            // Fit on training set
            auto ll = fit_local_logit(sorted_train_x,
                                      sorted_train_y,
                                      result.candidate_bandwidths[i],
                                      fit_quadratic,
                                      kernel_type,
                                      min_points,
                                      max_iterations,
                                      ridge_lambda,
                                      tolerance,
                                      false);

            auto train_preds = std::move(ll.grid_predictions); // these predictions take values at all x_grid points !!!

            // Compute validation Rf_error using linear interpolation
            for (size_t j = 0; j < valid_x.size(); j++) {
                double pred_j = interpolate_grid(valid_x[j], result.x_grid, train_preds);
                double brier_error = std::pow(pred_j - valid_y[j], 2);
                total_brier_error  += brier_error;
                result.bw_grid_errors[i][valid_idx[j]] = brier_error;
                valid_count++;
            }
        }

        result.mean_brier_errors[i] = valid_count > 0 ? total_brier_error / valid_count :
            std::numeric_limits<double>::infinity();
    }

    // Find optimal bandwidths for each Rf_error measure
    auto find_opt_idx = [](const std::vector<double>& errors) {
        return std::distance(errors.begin(),
                           std::min_element(errors.begin(), errors.end()));
    };


    result.opt_brier_bw_idx = find_opt_idx(result.mean_brier_errors);

    // debug
    // Rprintf("\nresult.opt_brier_bw_idx: %d\n", result.opt_brier_bw_idx);
    // print_vect(result.mean_brier_errors, "result.mean_brier_errors");

    std::set<int> unique_opt_indices = {
        result.opt_brier_bw_idx
    };


    if (!with_bw_grid_predictions) {
        // Store coefficients for each optimal bandwidth
        for (int idx : unique_opt_indices) {

            auto preds = fit_local_logit(x, y, result.candidate_bandwidths[idx],
                                         fit_quadratic, kernel_type, min_points,
                                         max_iterations, ridge_lambda, tolerance);
            result.bw_grid_predictions[idx] = std::move(preds);
        }
    } else {
        for (size_t i = 0; i < result.candidate_bandwidths.size(); i++) {
            auto preds = fit_local_logit(x, y, result.candidate_bandwidths[i],
                                         fit_quadratic, kernel_type, min_points,
                                         max_iterations, ridge_lambda, tolerance);
            result.bw_grid_predictions[i] = std::move(preds);
        }
    }

    // Identifying for each grid point the model of bw different from opt_i with
    // smaller CV Rf_error smaller than the CV of the opt_i model at that grid
    // point (if such model exist)
    int opt_i = result.opt_brier_bw_idx;
    double opt_bw = result.candidate_bandwidths[opt_i];

    std::vector<double> min_error_bw(grid_size, opt_bw);
    std::vector<int> min_error_bw_idx(grid_size, opt_i);
    for (int pti = 0; pti < grid_size; pti++) {

        double smallest_error = result.bw_grid_errors[opt_i][pti];
        int smallest_error_idx = opt_i;
        for (int i = 0; i < n_bws; i++) {
            if (result.bw_grid_errors[i][pti] < smallest_error) { // also check if the predicted values differ enough to warant the change of the model
                smallest_error = result.bw_grid_errors[i][pti];
                smallest_error_idx = i;
            }
        }

        min_error_bw[pti]     = smallest_error;
        min_error_bw_idx[pti] = smallest_error_idx;
    }

    // Smooth min_error_bw
    auto sm_res = mabilo(result.candidate_bandwidths, min_error_bw);
    auto sm_min_error_bw = std::move(sm_res.predictions);

    // clamp the smoothing to the range [min(min_error_bw), max(min_error_bw)]
    double min_error = *std::min_element(min_error_bw.begin(), min_error_bw.end());
    double max_error = *std::max_element(min_error_bw.begin(), min_error_bw.end());
    sm_min_error_bw = std::clamp(sm_min_error_bw, min_error, max_error);

    // Find bw's corresponding to the values of the smoothed min_error_bw
    std::vector<int> sm_min_error_bw_idx; // sm_min_error_bw_idx[pti] is the
                                          // index within
                                          // result.candidate_bandwidths of that
                                          // bandwith model averaged model with
                                          // the smallest CV Rf_error smaller than
                                          // opt_be Rf_error
    for (int i = 0; i < n_bws; i++) {
        // How to find the element of result.candidate_bandwidths closest to sm_min_error_bw[i] ??? <<---
    }

    // Compute the models needed for the model avareging using variable bw's identified above

    // Define a lambda that fits the model of specified bw at a specified grid location
    //auto get_grid_model <<---

    // grid_pt_models
    for (int pti = 0; pti < grid_size; pti++) {
        auto grid_model = get_grid_model();

        // for each point in the support of grid_model push that model (or the min info needed) to grid_pt_models <<---
    }

    // Generate final grid predicitons using model averaging
    // Use grid_pt_models to generate final grid_predictions

    // Find predictions at x locations
    for(int i = 0; i < n_pts; i++)
        result.predictions[i] = interpolate_grid(x[i], result.x_grid, result.bw_grid_predictions[result.opt_brier_bw_idx]);

    return result;
}


/**
 * @brief R interface for amagelogit local logistic regression
 *
 * @details Converts R inputs to C++ types, calls amagelogit(), and converts results
 * back to R objects. Returns a list containing grid predictions, Rf_error measures,
 * optimal bandwidths, and fitting information.
 *
 * @param x_r Numeric vector of predictor values
 * @param y_r Numeric vector of binary response values (0 or 1)
 * @param fit_quadratic_r Logical scalar for quadratic term inclusion
 * @param pilot_bandwidth_r Numeric scalar for fixed bandwidth (0 for CV)
 * @param kernel_type_r Integer scalar for kernel type
 * @param min_points_r Integer scalar for minimum local points
 * @param cv_folds_r Integer scalar for number of CV folds
 * @param n_bws_r Integer scalar for number of bandwidths to test
 * @param min_bw_factor_r Numeric scalar for minimum bandwidth factor
 * @param max_bw_factor_r Numeric scalar for maximum bandwidth factor
 * @param max_iterations_r Integer scalar for maximum iterations
 * @param ridge_lambda_r Numeric scalar for ridge parameter
 * @param tolerance_r Numeric scalar for convergence tolerance
 * @param with_bw_predictions_r Logical scalar for returning all bandwidths' predictions
 *
 * @return R list containing:
 *         - bw_grid_predictions: Matrix of predictions for each bandwidth
 *         - bw_grid_errors: Matrix of errors for each bandwidth
 *         - mean_brier_errors: Vector of CV errors
 *         - opt_brier_bw_idx: Optimal bandwidth index (1-based)
 *         - bws: Vector of tested bandwidths
 *         - fit_info: List of fitting parameters
 *
 * @throws Rf_error on NULL inputs or C++ exceptions
 */
SEXP S_amagelogit(
    SEXP x_r,
    SEXP y_r,
    SEXP grid_size_r,  // Add this parameter
    SEXP fit_quadratic_r,
    SEXP pilot_bandwidth_r,
    SEXP kernel_type_r,
    SEXP min_points_r,
    SEXP cv_folds_r,
    SEXP n_bws_r,
    SEXP min_bw_factor_r,
    SEXP max_bw_factor_r,
    SEXP max_iterations_r,
    SEXP ridge_lambda_r,
    SEXP tolerance_r,
    SEXP with_bw_predictions_r) {

    int n_protected = 0;

    try {
        if (x_r == R_NilValue || y_r == R_NilValue || grid_size_r == R_NilValue ||
            fit_quadratic_r == R_NilValue || pilot_bandwidth_r == R_NilValue ||
            kernel_type_r == R_NilValue || min_points_r == R_NilValue ||
            cv_folds_r == R_NilValue || n_bws_r == R_NilValue ||
            min_bw_factor_r == R_NilValue || max_bw_factor_r == R_NilValue ||
            max_iterations_r == R_NilValue || ridge_lambda_r == R_NilValue ||
            tolerance_r == R_NilValue || with_bw_predictions_r == R_NilValue) {
            Rf_error("Input arguments cannot be NULL");
        }

        int n_points = LENGTH(x_r);

        std::vector<double> x(REAL(x_r), REAL(x_r) + n_points);
        std::vector<double> y(REAL(y_r), REAL(y_r) + n_points);
        int grid_size = INTEGER(grid_size_r)[0];
        bool fit_quadratic = (LOGICAL(fit_quadratic_r)[0] == 1);
        double pilot_bandwidth = REAL(pilot_bandwidth_r)[0];
        int kernel_type = INTEGER(kernel_type_r)[0];
        int min_points = INTEGER(min_points_r)[0];
        int cv_folds = INTEGER(cv_folds_r)[0];
        int n_bws = INTEGER(n_bws_r)[0];
        double min_bw_factor = REAL(min_bw_factor_r)[0];
        double max_bw_factor = REAL(max_bw_factor_r)[0];
        int max_iterations = INTEGER(max_iterations_r)[0];
        double ridge_lambda = REAL(ridge_lambda_r)[0];
        double tolerance = REAL(tolerance_r)[0];
        bool with_bw_predictions = (LOGICAL(with_bw_predictions_r)[0] == 1);

        auto result = amagelogit(x, y, grid_size, fit_quadratic, pilot_bandwidth, kernel_type,
                                min_points, cv_folds, n_bws, min_bw_factor, max_bw_factor,
                                max_iterations, ridge_lambda, tolerance,
                                with_bw_predictions);

        // Create return list with updated size to include x_grid
        SEXP result_list = PROTECT(Rf_allocVector(VECSXP, 8)); n_protected++;

        // Updated names including x_grid
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 8)); n_protected++;
        SET_STRING_ELT(names, 0, Rf_mkChar("x_grid"));
        SET_STRING_ELT(names, 1, Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 2, Rf_mkChar("bw_grid_predictions"));
        SET_STRING_ELT(names, 3, Rf_mkChar("bw_grid_errors"));
        SET_STRING_ELT(names, 4, Rf_mkChar("mean_brier_errors"));
        SET_STRING_ELT(names, 5, Rf_mkChar("opt_brier_bw_idx"));
        SET_STRING_ELT(names, 6, Rf_mkChar("bws")); // candidate_bandwidths
        SET_STRING_ELT(names, 7, Rf_mkChar("fit_info"));

        // Add x_grid to result
        SEXP x_grid_r = PROTECT(Rf_allocVector(REALSXP, grid_size)); n_protected++;
        std::copy(result.x_grid.begin(), result.x_grid.end(), REAL(x_grid_r));
        SET_VECTOR_ELT(result_list, 0, x_grid_r);

        // Add predictions at original x points
        SEXP predictions_r = PROTECT(Rf_allocVector(REALSXP, n_points)); n_protected++;
        std::copy(result.predictions.begin(), result.predictions.end(), REAL(predictions_r));
        SET_VECTOR_ELT(result_list, 1, predictions_r);

        // Bandwidth grid predictions
        SEXP bw_pred_r;
        bw_pred_r = PROTECT(Rf_allocMatrix(REALSXP, grid_size,
                                           result.bw_grid_predictions.size())); n_protected++;
        for (size_t i = 0; i < result.bw_grid_predictions.size(); ++i) {
            if (!result.bw_grid_predictions[i].empty()) {
                // Copy column by column for column-major order
                std::copy(result.bw_grid_predictions[i].begin(),
                          result.bw_grid_predictions[i].end(),
                          REAL(bw_pred_r) + i * grid_size);
            }
        }
        SET_VECTOR_ELT(result_list, 2, bw_pred_r);

        // Bandwidth grid errors
        SEXP bw_err_r;
        bw_err_r = PROTECT(Rf_allocMatrix(REALSXP, grid_size,
                                           result.bw_grid_errors.size())); n_protected++;
        for (size_t i = 0; i < result.bw_grid_errors.size(); ++i) {
            if (!result.bw_grid_errors[i].empty()) {
                // Copy column by column for column-major order
                std::copy(result.bw_grid_errors[i].begin(),
                          result.bw_grid_errors[i].end(),
                          REAL(bw_err_r) + i * grid_size);
            }
        }
        SET_VECTOR_ELT(result_list, 3, bw_err_r);

        // Mean errors
        SEXP mean_brier_errors_r = PROTECT(Rf_allocVector(REALSXP, result.mean_brier_errors.size())); n_protected++;
        std::copy(result.mean_brier_errors.begin(), result.mean_brier_errors.end(), REAL(mean_brier_errors_r));
        SET_VECTOR_ELT(result_list, 4, mean_brier_errors_r);

        // Optimal indices
        SEXP opt_brier_bw_idx_r = PROTECT(Rf_allocVector(INTSXP, 1)); n_protected++;
        INTEGER(opt_brier_bw_idx_r)[0] = result.opt_brier_bw_idx + 1; // 1-base
        SET_VECTOR_ELT(result_list, 5, opt_brier_bw_idx_r);

        // Candidate bandwidths
        SEXP candidate_bandwidths_r = PROTECT(Rf_allocVector(REALSXP, result.candidate_bandwidths.size())); n_protected++;
        std::copy(result.candidate_bandwidths.begin(), result.candidate_bandwidths.end(),
                  REAL(candidate_bandwidths_r));
        SET_VECTOR_ELT(result_list, 6, candidate_bandwidths_r);

        // Fit info as a named list
        SEXP fit_info = PROTECT(Rf_allocVector(VECSXP, 9)); n_protected++;

        SEXP fit_info_names = PROTECT(Rf_allocVector(STRSXP, 9)); n_protected++;
        SET_STRING_ELT(fit_info_names, 0, Rf_mkChar("fit_quadratic"));
        SET_STRING_ELT(fit_info_names, 1, Rf_mkChar("pilot_bandwidth"));
        SET_STRING_ELT(fit_info_names, 2, Rf_mkChar("kernel_type"));
        SET_STRING_ELT(fit_info_names, 3, Rf_mkChar("cv_folds"));
        SET_STRING_ELT(fit_info_names, 4, Rf_mkChar("min_bw_factor"));
        SET_STRING_ELT(fit_info_names, 5, Rf_mkChar("max_bw_factor"));
        SET_STRING_ELT(fit_info_names, 6, Rf_mkChar("max_iterations"));
        SET_STRING_ELT(fit_info_names, 7, Rf_mkChar("ridge_lambda"));
        SET_STRING_ELT(fit_info_names, 8, Rf_mkChar("tolerance"));

        SEXP fit_quad_r = PROTECT(Rf_allocVector(LGLSXP, 1));
        (LOGICAL(fit_quad_r)[0] == 1) = result.fit_quadratic;
        SET_VECTOR_ELT(fit_info, 0, fit_quad_r); n_protected++;

        SEXP pilot_bw_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(pilot_bw_r)[0] = result.pilot_bandwidth;
        SET_VECTOR_ELT(fit_info, 1, pilot_bw_r); n_protected++;

        SEXP kernel_type_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(kernel_type_out_r)[0] = result.kernel_type;
        SET_VECTOR_ELT(fit_info, 2, kernel_type_out_r); n_protected++;

        SEXP cv_folds_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(cv_folds_out_r)[0] = result.cv_folds;
        SET_VECTOR_ELT(fit_info, 3, cv_folds_out_r); n_protected++;

        SEXP min_bw_factor_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(min_bw_factor_out_r)[0] = result.min_bw_factor;
        SET_VECTOR_ELT(fit_info, 4, min_bw_factor_out_r); n_protected++;

        SEXP max_bw_factor_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(max_bw_factor_out_r)[0] = result.max_bw_factor;
        SET_VECTOR_ELT(fit_info, 5, max_bw_factor_out_r); n_protected++;

        SEXP max_iter_out_r = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(max_iter_out_r)[0] = result.max_iterations;
        SET_VECTOR_ELT(fit_info, 6, max_iter_out_r); n_protected++;

        SEXP ridge_lambda_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(ridge_lambda_out_r)[0] = result.ridge_lambda;
        SET_VECTOR_ELT(fit_info, 7, ridge_lambda_out_r); n_protected++;

        SEXP tolerance_out_r = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(tolerance_out_r)[0] = result.tolerance;
        SET_VECTOR_ELT(fit_info, 8, tolerance_out_r); n_protected++;

        Rf_setAttrib(fit_info, R_NamesSymbol, fit_info_names);
        SET_VECTOR_ELT(result_list, 7, fit_info);

        // Set names for the main list
        Rf_setAttrib(result_list, R_NamesSymbol, names);

        UNPROTECT(n_protected);

        return result_list;
    }
    catch (const std::exception& e) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("C++ Rf_error in amagelogit: %s", e.what());
    }
    catch (...) {
        if (n_protected > 0) UNPROTECT(n_protected);
        Rf_error("Unknown Rf_error in amagelogit");
    }
}
#endif
