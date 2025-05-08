#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <unordered_set>
#include <queue>
#include <utility>
#include <set>
#include <unordered_set>
#include <numeric>
#include <tuple>

#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "graph_utils.hpp"
#include "kernels.h"
#include "ulogit.hpp"

struct gmagelog_t {
    // Additional members specific to graph-based implementation
    uniform_grid_graph_t grid_graph;          // Stores the expanded uniform grid graph
    std::map<int, double> grid_predictions;  // Maps grid vertices to their computed values
    double graph_diameter;

    // bandwidth grid
    std::vector<double> candidate_bandwidths;                   ///< Grid of bandwidths tested during optimization

    // grid-based members
    int opt_bw_idx;                                             ///< Index of bandwidth with minimal mean Brier error
    std::vector<std::map<int, double>> bw_grid_predictions;      ///< Predictions for each bandwidth in LOOCV or CV estimation
    std::vector<double> bw_mean_errors;                          ///< Mean Brier error for each candidate bandwidth

    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> predictions;                           ///< Predictions at input graph vertices

    std::vector<std::string> warnings;  // Vector to store warning messages

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

struct ugg_local_nbhd_t {
    // non-grid elements of a local neighborhood - used for model building
    std::vector<int> local_indices;      // indices of non-grid vertieces of result.grid_graph within the given neighhborhood;
    std::vector<double> local_d;         // distance (accounting for edge lengths) from grid_vertex to all non-grid vertices such that the distance is not more than the bandwidth value
    std::vector<double> local_y;         // y values at local_idices
    std::vector<double> local_w;         // weights computed using a kernel function applied to local_x
    std::vector<double> local_errors;    // LOOCV prediction errors estimated by a local model
    std::vector<double> local_predictions;

    // grid elements of a local neighborhood - used for model averaging
    std::vector<int> local_grid_indices; // indices of grid vertieces of result.grid_graph within the given neighhborhood; local_grid_idices[i] is the index of the vertex whose distance to grid_vertex is local_grid_x[i]
    std::vector<double> local_grid_d;    // distance (accounting for edge lengths) from grid_vertex to grid vertices such that the distance is not more than the bandwidth value
    std::vector<double> local_grid_w;    // weights computed using a kernel function applied to local_grid_x
    std::vector<double> local_grid_predictions; // predictions from a local model at the grid points of the local neibhborhood
};

gmagelog_t gmagelog(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    const std::vector<double>& y,
    int grid_size,
    bool fit_quadratic = false,
    double pilot_bandwidth = -1.0,
    int kernel_type = 0,
    int cv_folds = 0,
    int n_bws = 15,
    double min_bw_factor = 0.1,
    double max_bw_factor = 2.0,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8,
    int start_vertex = 0) {

    bool with_errors = true; // cv_folds > 0; for now we use LOOCV prediction error estimates from eigen_ulogit_fit()

    // Initialize result structure
    gmagelog_t result;
    result.fit_quadratic = fit_quadratic;
    result.pilot_bandwidth = pilot_bandwidth;
    result.kernel_type = kernel_type;
    result.cv_folds = cv_folds;
    result.min_bw_factor = min_bw_factor;
    result.max_bw_factor = max_bw_factor;
    result.max_iterations = max_iterations;
    result.ridge_lambda = ridge_lambda;
    result.tolerance = tolerance;

    // Create uniform grid graph
    double snap_tolerance = 0.1;
    result.grid_graph = create_uniform_grid_graph(
        input_adj_list,
        input_weight_list,
        grid_size,
        start_vertex,
        snap_tolerance
        );

    initialize_kernel(kernel_type, 1.0);

    int n_data_vertices = y.size();
    int n_grid_vertices = result.grid_graph.grid_vertices.size(); // note that the actual number of grid points in grid_graph does not have to be grid_size - but it should be close to this number

    auto fit_local_model = [&y, &n_data_vertices, &n_grid_vertices, &fit_quadratic, &max_iterations, &ridge_lambda, &tolerance, &with_errors](
        const uniform_grid_graph_t& grid_graph,
        int grid_vertex,
        double bandwidth) -> ugg_local_nbhd_t {

        ugg_local_nbhd_t result;
        // Pre-reserve space in vectors to avoid reallocations
        result.local_indices.reserve(n_data_vertices);
        result.local_d.reserve(n_data_vertices);
        result.local_y.reserve(n_data_vertices);
        result.local_w.reserve(n_data_vertices);

        result.local_grid_indices.reserve(n_grid_vertices);
        result.local_grid_d.reserve(n_grid_vertices);
        result.local_grid_w.reserve(n_grid_vertices);

        std::unordered_map<int, double> vertex_distances;
        vertex_distances[grid_vertex] = 0.0;

        std::queue<std::pair<int, double>> bfs_queue;
        bfs_queue.push({grid_vertex, 0.0});

        while (!bfs_queue.empty()) {
            auto [current_vertex, current_distance] = bfs_queue.front();
            bfs_queue.pop();

            // If vertex has an observation
            if (current_vertex < n_data_vertices) {
                result.local_indices.push_back(current_vertex);
                result.local_d.push_back(current_distance);
                result.local_y.push_back(y[current_vertex]);
            }

            // Check if it's a grid vertex (using unordered_set)
            if (grid_graph.grid_vertices.contains(current_vertex) &&
                current_vertex != grid_vertex) {
                result.local_grid_indices.push_back(current_vertex);
                result.local_grid_d.push_back(current_distance);
            }

            // Explore neighbors within bandwidth
            for (const auto& edge_info : grid_graph.adjacency_list[current_vertex]) {
                int neighbor = edge_info.vertex;
                double edge_length = edge_info.weight;
                double neighbor_distance = current_distance + edge_length;

                if (neighbor_distance > bandwidth ||
                    vertex_distances.find(neighbor) != vertex_distances.end()) {
                    continue;
                }

                vertex_distances[neighbor] = neighbor_distance;
                bfs_queue.push({neighbor, neighbor_distance});
            }
        }

        // computing non-grid points weights
        int n_nbhd_nongrid_vertices = result.local_d.size();
        result.local_w.resize(n_nbhd_nongrid_vertices); // Resize local_w before using it
        kernel_fn(result.local_d.data(), n_nbhd_nongrid_vertices, result.local_w.data());

        // Normalize weights with protection against zero sum
        double sum_weights = std::accumulate(result.local_w.begin(), result.local_w.end(), 0.0);
        if (sum_weights < std::numeric_limits<double>::epsilon()) {
            sum_weights = std::numeric_limits<double>::epsilon();
        }
        for (auto& w : result.local_w) {
            w /= sum_weights;
        }

        // Fitting local weighted logistic regression model
        eigen_ulogit_t fit_result = eigen_ulogit_fit(
            result.local_d.data(), // should be local_x.data() but
            result.local_y.data(),
            result.local_w,
            fit_quadratic,
            max_iterations,
            ridge_lambda,
            tolerance,
            with_errors
            );

        result.local_predictions = fit_result.predictions;
        result.local_errors = fit_result.loocv_brier_errors; // those are computed at local non-grid points. What if estimated error[grid_v] as fit_result.loocv_brier_errors at the non-grid vertex closest to grid_v?
        result.local_grid_predictions = fit_result.predict(result.local_grid_d); // ideally we want here to use local_grid_x not local_grid_d; predicting conditional expectation values at the grid point within the support of the model

        int n_nbhd_grid_vertices = result.local_grid_d.size();
        result.local_grid_w.resize(n_nbhd_grid_vertices); // Resize local_w before using it
        kernel_fn(result.local_grid_d.data(), n_nbhd_grid_vertices, result.local_grid_w.data());

        // Normalize weights with protection against zero sum
        sum_weights = std::accumulate(result.local_grid_w.begin(), result.local_grid_w.end(), 0.0);
        if (sum_weights < std::numeric_limits<double>::epsilon()) {
            sum_weights = std::numeric_limits<double>::epsilon();
        }
        for (auto& w : result.local_grid_w) {
            w /= sum_weights;
        }

        return result;
    };

    auto estimate_model_avgd_predictions_and_errors = [&](
        double bandwidth,
        const uniform_grid_graph_t& grid_graph,
        const auto& fit_local_model_fn,
        std::vector<std::string>& warnings
        ) {

        // Track problematic vertices
        std::vector<int> no_prediction_grid_vertices;
        std::vector<int> no_prediction_data_vertices;

        // Collecting prediction and error data from different models for model averaging and estimation of the mean error for the given bandwidth
        std::map<int, std::vector<std::pair<double,double>>> grid_wpreds_map; // grid_wpreds_map[grid_idx][i].first/second = weight/prediction of the i-th local model whose support contain 'grid_idx' grid vertex
        std::vector<std::vector<std::tuple<double,double,double>>> w_pred_errors(n_data_vertices); // w_pred_errors[i] a vector of {weight,prediction,error} pairs for non-grid vertices within the support of different models
        for (int grid_vertex : grid_graph.grid_vertices) {
            ugg_local_nbhd_t nbhd_model = fit_local_model_fn(
                result.grid_graph,
                grid_vertex,
                bandwidth
            );

            for (size_t i = 0; i < nbhd_model.local_grid_indices.size(); i++) {
                grid_wpreds_map[nbhd_model.local_grid_indices[i]].emplace_back(
                    nbhd_model.local_grid_w[i],
                    nbhd_model.local_grid_predictions[i]
                    );
            }

            for (size_t i = 0; i < nbhd_model.local_indices.size(); i++) {
                w_pred_errors[nbhd_model.local_indices[i]].emplace_back(
                    nbhd_model.local_w[i],
                    nbhd_model.local_predictions[i],
                    nbhd_model.local_errors[i]
                    );
            }
        }

        // Compute grid vertices weighted averaged predictions
        std::map<int, double> grid_predictions;
        for (const auto& grid_idx : grid_graph.grid_vertices) {
            const auto& v = grid_wpreds_map[grid_idx];
            if (v.empty()) {
                grid_predictions[grid_idx] = std::numeric_limits<double>::quiet_NaN();
                no_prediction_grid_vertices.push_back(grid_idx);
                continue;
            }

            double total_weight = 0.0;
            double weighted_prediction = 0.0;
            for (const auto& p : v) {
                weighted_prediction += p.first * p.second;
                total_weight += p.first;
            }

            if (total_weight > 0.0) {
                grid_predictions[grid_idx] = weighted_prediction / total_weight;
            } else {
                grid_predictions[grid_idx] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // Compute non-grid vertices weighted averaged LOOCV errors
        std::vector<double> errors(n_data_vertices);
        std::vector<double> predictions(n_data_vertices);
        for (int i = 0; i < n_data_vertices; i++) {
            const auto& v = w_pred_errors[i];

            if (v.empty()) {
                errors[i] = std::numeric_limits<double>::quiet_NaN();
                predictions[i] = std::numeric_limits<double>::quiet_NaN();
                no_prediction_data_vertices.push_back(i);
                continue;
            }

            double total_weight = 0.0;
            double weighted_error = 0.0;
            double weighted_prediction = 0.0;
            for (const auto& [w, p, e] : v) {
                weighted_prediction += w * p;
                weighted_error += w * e;
                total_weight += w;
            }

            if (total_weight > 0.0) {
                errors[i] = weighted_error / total_weight;
                predictions[i] = weighted_prediction / total_weight;
            } else {
                errors[i] = std::numeric_limits<double>::quiet_NaN();
                predictions[i] = errors[i];
            }
        }

        double sum = 0.0;
        int count = 0;
        for (double error : errors) {
            if (!std::isnan(error)) {
                sum += error;
                count++;
            }
        }
        double mean_error = count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();

        // Generate warnings based on coverage
        if (!no_prediction_grid_vertices.empty()) {
            std::string msg = "No predictions available for " +
                std::to_string(no_prediction_grid_vertices.size()) +
                " grid vertices (" +
                std::to_string(100.0 * no_prediction_grid_vertices.size() / grid_graph.grid_vertices.size()) +
                "% of grid). This may indicate insufficient bandwidth or disconnected graph components.";
            warnings.push_back(msg);

            if (no_prediction_grid_vertices.size() > grid_graph.grid_vertices.size() / 2) {
                warnings.push_back("More than 50% of grid vertices have no predictions. Consider increasing bandwidth or checking graph connectivity.");
            }
        }

        if (!no_prediction_data_vertices.empty()) {
            std::string msg = "No predictions available for " +
                std::to_string(no_prediction_data_vertices.size()) +
                " data vertices (" +
                std::to_string(100.0 * no_prediction_data_vertices.size() / n_data_vertices) +
                "% of data). This may indicate isolated vertices or insufficient bandwidth coverage.";
            warnings.push_back(msg);

            if (no_prediction_data_vertices.size() > n_data_vertices / 2) {
                warnings.push_back("More than 50% of data vertices have no predictions. Consider increasing bandwidth or checking data vertex distribution.");
            }
        }


        return std::tuple(grid_predictions, predictions, mean_error);
    };

    // If pilot bandwidth provided, use it directly
    if (pilot_bandwidth > 0) {

        double mean_error;  // Declare first
        std::tie(result.grid_predictions,
                 result.predictions,
                 mean_error) = estimate_model_avgd_predictions_and_errors(
                     pilot_bandwidth,
                     result.grid_graph,
                     fit_local_model,
                     result.warnings
                     );
    } else {
        // Creating candidate_bandwidths
        // Compute the graph diameter from the perspective of start_vertex - that is the longest shortest path starting at start_vertex
        result.graph_diameter = get_vertex_eccentricity(
            input_adj_list,
            input_weight_list,
            start_vertex
            );

        if (result.graph_diameter <= 0) {
            error("result.graph_diameter: %f - Invalid graph diameter", result.graph_diameter);
        }

        result.candidate_bandwidths.resize(n_bws);
        double min_bw = min_bw_factor * result.graph_diameter;
        double max_bw = max_bw_factor * result.graph_diameter;
        double bw_step = (max_bw - min_bw) / (n_bws - 1);

        for(int i = 0; i < n_bws; i++) {
            result.candidate_bandwidths[i] = min_bw + i * bw_step;
        }

        result.bw_grid_predictions.resize(n_bws);
        result.bw_predictions.resize(n_bws);
        result.bw_mean_errors.resize(n_bws);

        for(int i = 0; i < n_bws; i++) {
            double bw = result.candidate_bandwidths[i];
            std::tie(result.bw_grid_predictions[i],
                     result.bw_predictions[i],
                     result.bw_mean_errors[i]) = estimate_model_avgd_predictions_and_errors(
                         bw,
                         result.grid_graph,
                         fit_local_model,
                         result.warnings
                         );
        }
    }

    // Cross-validation code would reuse the same collect_local_neighborhood lambda
    // without recreating it for each bandwidth evaluation

    // At the end of gmagelog() function, before return:
    if (!result.warnings.empty()) {
        Rprintf("\nWarnings:\n");
        for (const auto& warning : result.warnings) {
            Rprintf("  - %s\n", warning.c_str());
        }
    }

    return result;
}

