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

#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.h"
#include "cpp_utils.h"



    // Define the kernel weight computation lambda before collect_local_neighborhood
    auto compute_kernel_weight = [](double normalized_distance, int kernel_type) -> double {
        // Ensure normalized distance is non-negative
        normalized_distance = std::abs(normalized_distance);

        // Early return for points outside the kernel support
        if (normalized_distance >= 1.0) {
            return 0.0;
        }

        // Compute kernel weight based on the selected type
        switch (kernel_type) {
        case 0: {  // Epanechnikov kernel: K(u) = 3/4 * (1 - u^2) for |u| < 1
            double u_squared = normalized_distance * normalized_distance;
            return 0.75 * (1.0 - u_squared);
        }

        case 1: {  // Triangular kernel: K(u) = (1 - |u|) for |u| < 1
            return 1.0 - normalized_distance;
        }

        case 2: {  // Gaussian kernel: K(u) = exp(-u^2/2) / sqrt(2π)
            // Note: We don't restrict to |u| < 1 for Gaussian
            static const double SQRT_2PI = std::sqrt(2.0 * M_PI);
            double u_squared = normalized_distance * normalized_distance;
            return std::exp(-0.5 * u_squared) / SQRT_2PI;
        }

        case 3: {  // Quartic/Biweight kernel: K(u) = 15/16 * (1 - u^2)^2 for |u| < 1
            double u_squared = normalized_distance * normalized_distance;
            double one_minus_u_squared = 1.0 - u_squared;
            return (15.0/16.0) * one_minus_u_squared * one_minus_u_squared;
        }

        case 4: {  // Tricube kernel: K(u) = 70/81 * (1 - |u|^3)^3 for |u| < 1
            double u_cubed = normalized_distance * normalized_distance * normalized_distance;
            double one_minus_u_cubed = 1.0 - u_cubed;
            return (70.0/81.0) * one_minus_u_cubed * one_minus_u_cubed * one_minus_u_cubed;
        }

        case 5: {  // Cosine kernel: K(u) = π/4 * cos(πu/2) for |u| < 1
            static const double PI_OVER_4 = M_PI / 4.0;
            return PI_OVER_4 * std::cos(M_PI_2 * normalized_distance);
        }

        default: {  // Default to Epanechnikov kernel
            double u_squared = normalized_distance * normalized_distance;
            return 0.75 * (1.0 - u_squared);
        }
        }
    };


struct gmagelog_t {
    // Additional members specific to graph-based implementation
    uniform_grid_graph_t grid_graph;          // Stores the expanded uniform grid graph
    std::map<int, double> grid_predictions;  // Maps grid vertices to their computed values

    // bandwidth grid
    std::vector<double> candidate_bandwidths;                   ///< Grid of bandwidths tested during optimization

    // Mean errors and optimal indices
    std::vector<double> mean_brier_errors;                      ///< Mean Brier error for each candidate bandwidth
    int opt_brier_bw_idx;                                       ///< Index of bandwidth with minimal mean Brier error

    // grid-based members
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

struct ugg_local_nbhd_t {
    std::vector<int> local_indices;      // indices of non-grid vertieces of result.grid_graph within the given neighhborhood;
    std::vector<double> local_x;         // distance (accounting for edge lengths) from grid_vertex to all non-grid vertices such that the distance is not more than the bandwidth value
    std::vector<double> local_y;         // y values at local_idices
    std::vector<double> local_w;         // weights computed using a kernel function applied to local_x
    std::vector<int> local_grid_indices; // indices of grid vertieces of result.grid_graph within the given neighhborhood; local_grid_idices[i] is the index of the vertex whose distance to grid_vertex is local_grid_x[i]
    std::vector<double> local_grid_x;    // distance (accounting for edge lengths) from grid_vertex to grid vertices such that the distance is not more than the bandwidth value
    std::vector<double> local_grid_w;    // weights computed using a kernel function applied to local_grid_x
};

gmagelog_t gmagelog(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    const std::vector<int>& data_vertices,
    const std::vector<double>& y,
    int grid_size,
    bool fit_quadratic = false,
    double pilot_bandwidth = -1.0,
    int kernel_type = 0,
    int min_points = 5,
    int cv_folds = 0,
    int n_bws = 15,
    double min_bw_factor = 0.1,
    double max_bw_factor = 2.0,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8,
    bool with_bw_grid_predictions = false,
    int start_vertex = 0) {

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

    // Define the kernel weight computation lambda before collect_local_neighborhood
    auto compute_kernel_weight = [](double normalized_distance, int kernel_type) -> double {
        // Ensure normalized distance is non-negative
        normalized_distance = std::abs(normalized_distance);

        // Early return for points outside the kernel support
        if (normalized_distance >= 1.0) {
            return 0.0;
        }

        // Compute kernel weight based on the selected type
        switch (kernel_type) {
        case 0: {  // Epanechnikov kernel: K(u) = 3/4 * (1 - u^2) for |u| < 1
            double u_squared = normalized_distance * normalized_distance;
            return 0.75 * (1.0 - u_squared);
        }

        case 1: {  // Triangular kernel: K(u) = (1 - |u|) for |u| < 1
            return 1.0 - normalized_distance;
        }

        case 2: {  // Gaussian kernel: K(u) = exp(-u^2/2) / sqrt(2π)
            // Note: We don't restrict to |u| < 1 for Gaussian
            static const double SQRT_2PI = std::sqrt(2.0 * M_PI);
            double u_squared = normalized_distance * normalized_distance;
            return std::exp(-0.5 * u_squared) / SQRT_2PI;
        }

        case 3: {  // Quartic/Biweight kernel: K(u) = 15/16 * (1 - u^2)^2 for |u| < 1
            double u_squared = normalized_distance * normalized_distance;
            double one_minus_u_squared = 1.0 - u_squared;
            return (15.0/16.0) * one_minus_u_squared * one_minus_u_squared;
        }

        case 4: {  // Tricube kernel: K(u) = 70/81 * (1 - |u|^3)^3 for |u| < 1
            double u_cubed = normalized_distance * normalized_distance * normalized_distance;
            double one_minus_u_cubed = 1.0 - u_cubed;
            return (70.0/81.0) * one_minus_u_cubed * one_minus_u_cubed * one_minus_u_cubed;
        }

        case 5: {  // Cosine kernel: K(u) = π/4 * cos(πu/2) for |u| < 1
            static const double PI_OVER_4 = M_PI / 4.0;
            return PI_OVER_4 * std::cos(M_PI_2 * normalized_distance);
        }

        default: {  // Default to Epanechnikov kernel
            double u_squared = normalized_distance * normalized_distance;
            return 0.75 * (1.0 - u_squared);
        }
        }
    };

    // Define neighborhood collection lambda at the outer scope
    // This lambda is created only once and can be reused throughout the function
    auto collect_local_neighborhood = [&kernel_type](
        const uniform_grid_graph_t& grid_graph,
        int grid_vertex,
        double bandwidth,
        const std::vector<int>& data_vertices,
        const std::vector<double>& responses) -> ugg_local_nbhd_t {

        ugg_local_nbhd_t result;

        // Track visited vertices and their distances from grid_vertex
        std::unordered_map<int, double> vertex_distances;
        vertex_distances[grid_vertex] = 0.0;

        std::queue<std::pair<int, double>> bfs_queue;
        bfs_queue.push({grid_vertex, 0.0});

        // Create lookup maps/sets only once per neighborhood collection
        std::unordered_map<int, size_t> data_vertex_indices;
        data_vertex_indices.reserve(data_vertices.size());  // Prevent rehashing
        for (size_t i = 0; i < data_vertices.size(); ++i) {
            data_vertex_indices[data_vertices[i]] = i;
        }

        std::unordered_set<int> grid_vertex_set(
            grid_graph.grid_vertices.begin(),
            grid_graph.grid_vertices.end()
        );

        // Pre-reserve space in vectors to avoid reallocations
        result.local_indices.reserve(data_vertices.size());
        result.local_x.reserve(data_vertices.size());
        result.local_y.reserve(data_vertices.size());
        result.local_w.reserve(data_vertices.size());
        result.local_grid_indices.reserve(grid_graph.grid_vertices.size());
        result.local_grid_x.reserve(grid_graph.grid_vertices.size());
        result.local_grid_w.reserve(grid_graph.grid_vertices.size());

        while (!bfs_queue.empty()) {
            auto [current_vertex, current_distance] = bfs_queue.front();
            bfs_queue.pop();

            // Process current vertex if it's a data point
            auto data_it = data_vertex_indices.find(current_vertex);
            if (data_it != data_vertex_indices.end()) {
                result.local_indices.push_back(current_vertex);
                result.local_x.push_back(current_distance);
                result.local_y.push_back(responses[data_it->second]);
                result.local_w.push_back(compute_kernel_weight(
                    current_distance / bandwidth,
                    kernel_type
                ));
            }

            // Process current vertex if it's a grid point (except the center)
            if (grid_vertex_set.count(current_vertex) && current_vertex != grid_vertex) {
                result.local_grid_indices.push_back(current_vertex);
                result.local_grid_x.push_back(current_distance);
                result.local_grid_w.push_back(compute_kernel_weight(
                    current_distance / bandwidth,
                    kernel_type
                ));
            }

            // Explore neighbors within bandwidth
            for (const auto& edge_info : grid_graph.ugg[current_vertex]) {
                int neighbor = edge_info.target;
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

        return result;
    };

    // Define the local logistic regression fitting function
    // Now it takes collect_local_neighborhood as a parameter
    auto fit_local_logistic_graph = [&](
        const std::vector<int>& vertices,
        const std::vector<double>& responses,
        double bandwidth,
        const uniform_grid_graph_t& grid_graph,
        const auto& neighborhood_collector) {

        std::map<int, double> grid_predictions;

        for (int grid_vertex : grid_graph.grid_vertices) {
            ugg_local_nbhd_t nbhd = neighborhood_collector(
                grid_graph,
                grid_vertex,
                bandwidth,
                vertices,
                responses
            );

            if (nbhd.local_x.size() >= min_points) {
                eigen_ulogit_t fit = eigen_ulogit_fit(
                    nbhd.local_x.data(),
                    nbhd.local_y.data(),
                    nbhd.local_w,
                    fit_quadratic,
                    max_iterations,
                    ridge_lambda,
                    tolerance
                );

                grid_predictions[grid_vertex] = fit.predictions[0];
            }
        }

        return grid_predictions;
    };

    // If pilot bandwidth provided, use it directly
    if (pilot_bandwidth > 0) {
        result.grid_predictions = fit_local_logistic_graph(
            data_vertices,
            y,
            pilot_bandwidth,
            result.grid_graph,
            collect_local_neighborhood  // Pass the collector function
        );

        // ... rest of the pilot bandwidth code ...
    }

    // Cross-validation code would reuse the same collect_local_neighborhood lambda
    // without recreating it for each bandwidth evaluation

    return result;
}

#if 0
gmagelog_t v0gmagelog(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    const std::vector<int>& data_vertices,    // Vertices where we have observations
    const std::vector<double>& y,             // Binary response values
    int grid_size,                            // Number of grid points to create
    bool fit_quadratic = false,
    double pilot_bandwidth = -1.0,
    int kernel_type = 7,
    int min_points = 5,
    int cv_folds = 0,
    int n_bws = 25,
    double min_bw_factor = 0.1,
    double max_bw_factor = 2.0,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8,
    bool with_bw_grid_predictions = false,
    int start_vertex = 0) {

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
    result.grid_graph = create_uniform_grid_graph(
        input_adj_list,
        input_weight_list,
        grid_size,
        start_vertex
    );

    // Compute shortest path distances from start_vertex to all vertices
    result.vertex_distances = compute_shortest_paths(
        input_adj_list,
        input_weight_list,
        start_vertex
    );

    // Create bandwidth grid based on graph diameter
    double graph_diameter = *std::max_element(
        result.vertex_distances.begin(),
        result.vertex_distances.end()
    );

    result.candidate_bandwidths.resize(n_bws);
    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;
    double bw_step = (max_bw - min_bw) / (n_bws - 1);

    for(int i = 0; i < n_bws; i++) {
        result.candidate_bandwidths[i] = min_bw + i * bw_step;
    }


    // Function to collect local neighborhood data within bandwidth using BFS
    ugg_local_nbhd_t collect_local_neighborhood(
        const uniform_grid_graph_t& grid_graph,
        int grid_vertex,
        double bandwidth,
        const std::vector<int>& data_vertices,
        const std::vector<double>& responses,
        int kernel_type) {

        ugg_local_nbhd_t result;

        // Track visited vertices and their distances from grid_vertex
        std::unordered_map<int, double> vertex_distances;
        vertex_distances[grid_vertex] = 0.0;

        // BFS queue stores pairs of (vertex, cumulative_distance)
        std::queue<std::pair<int, double>> bfs_queue;
        bfs_queue.push({grid_vertex, 0.0});

        // Create index mapping for data vertices for efficient lookup
        std::unordered_map<int, size_t> data_vertex_indices;
        for (size_t i = 0; i < data_vertices.size(); ++i) {
            data_vertex_indices[data_vertices[i]] = i;
        }

        // Create set of grid vertices for efficient lookup
        std::unordered_set<int> grid_vertex_set(
            grid_graph.grid_vertices.begin(),
            grid_graph.grid_vertices.end()
            );

        while (!bfs_queue.empty()) {
            auto [current_vertex, current_distance] = bfs_queue.front();
            bfs_queue.pop();

            // Process current vertex if it's a data point
            auto data_it = data_vertex_indices.find(current_vertex);
            if (data_it != data_vertex_indices.end()) {
                result.local_indices.push_back(current_vertex);
                result.local_x.push_back(current_distance);
                result.local_y.push_back(responses[data_it->second]);
                result.local_w.push_back(compute_kernel_weight(
                                             current_distance / bandwidth,
                                             kernel_type
                                             ));
            }

            // Process current vertex if it's a grid point (except the center)
            if (grid_vertex_set.count(current_vertex) && current_vertex != grid_vertex) {
                result.local_grid_indices.push_back(current_vertex);
                result.local_grid_x.push_back(current_distance);
                result.local_grid_w.push_back(compute_kernel_weight(
                                                  current_distance / bandwidth,
                                                  kernel_type
                                                  ));
            }

            // Explore neighbors within bandwidth
            for (const auto& edge_info : grid_graph.ugg[current_vertex]) {
                int neighbor = edge_info.target;
                double edge_length = edge_info.weight;
                double neighbor_distance = current_distance + edge_length;

                // Skip if we've seen this vertex or if it's too far
                if (neighbor_distance > bandwidth ||
                    vertex_distances.find(neighbor) != vertex_distances.end()) {
                    continue;
                }

                // Record distance and add to queue
                vertex_distances[neighbor] = neighbor_distance;
                bfs_queue.push({neighbor, neighbor_distance});
            }
        }

        return result;
    }

// Updated fit_local_logistic_graph function using the new neighborhood collection
    auto fit_local_logistic_graph = [&](
        const std::vector<int>& vertices,
        const std::vector<double>& responses,
        double bandwidth,
        const uniform_grid_graph_t& grid_graph) {

        std::map<int, double> grid_predictions;

        // For each grid vertex
        for (int grid_vertex : grid_graph.grid_vertices) {
            // Collect local neighborhood using BFS
            ugg_local_nbhd_t nbhd = collect_local_neighborhood(
                grid_graph,
                grid_vertex,
                bandwidth,
                vertices,
                responses,
                kernel_type
                );

            // Fit local model if enough points
            if (nbhd.local_x.size() >= min_points) {
                eigen_ulogit_t fit = eigen_ulogit_fit(
                    nbhd.local_x.data(),
                    nbhd.local_y.data(),
                    nbhd.local_w,
                    fit_quadratic,
                    max_iterations,
                    ridge_lambda,
                    tolerance
                    );

                // Store prediction for this grid vertex
                grid_predictions[grid_vertex] = fit.predictions[0];
            }
        }

        return grid_predictions;
    };


    // Lambda function for fitting local logistic regression on graph
    auto fit_local_logistic_graph = [&](
        const std::vector<int>& vertices,
        const std::vector<double>& responses,
        double bandwidth,
        const uniform_grid_graph_t& grid_graph) {

        std::map<int, double> grid_predictions;

        // For each grid vertex
        for(int grid_vertex : grid_graph.grid_vertices) {
            // Collect local data within bandwidth neighborhood
            std::vector<double> local_x, local_y, local_w;

            for(size_t i = 0; i < vertices.size(); i++) {
                double dist = graph_distance(grid_vertex, vertices[i]);
                if(dist <= bandwidth) {
                    local_x.push_back(dist);  // Use distance as predictor
                    local_y.push_back(responses[i]);

                    // Compute kernel weight
                    double kernel_weight = compute_kernel_weight(
                        dist / bandwidth,
                        kernel_type
                    );
                    local_w.push_back(kernel_weight);
                }
            }

            // Fit local model if enough points
            if(local_x.size() >= min_points) {
                eigen_ulogit_t fit = eigen_ulogit_fit(
                    local_x.data(),
                    local_y.data(),
                    local_w,
                    fit_quadratic,
                    max_iterations,
                    ridge_lambda,
                    tolerance
                );

                // Store prediction for this grid vertex
                grid_predictions[grid_vertex] = fit.predictions[0];  // Predict at center
            }
        }

        return grid_predictions;
    };

    // If pilot bandwidth provided, use it directly
    if(pilot_bandwidth > 0) {
        result.grid_predictions = fit_local_logistic_graph(
            data_vertices,
            y,
            pilot_bandwidth,
            result.grid_graph
        );

        // Interpolate predictions for all vertices
        result.predictions = interpolate_graph_values(
            result.grid_predictions,
            result.grid_graph,
            data_vertices
        );

        return result;
    }

    // Otherwise perform bandwidth selection via cross-validation
    // ... implement CV logic similar to magelog() but using graph distances

    return result;
}

// Helper function to compute graph distances
double graph_distance(
    const uniform_grid_graph_t& graph,
    int vertex1,
    int vertex2) {
    // Implement Dijkstra's algorithm to find shortest path distance
    // between vertex1 and vertex2 on the expanded grid graph
}

// Helper function to interpolate values on graph
std::vector<double> interpolate_graph_values(
    const std::map<int, double>& grid_values,
    const uniform_grid_graph_t& graph,
    const std::vector<int>& query_vertices) {
    // Implement interpolation of values from grid vertices
    // to arbitrary vertices using graph structure
}
#endif
