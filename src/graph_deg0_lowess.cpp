#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions
// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <algorithm>                // For std::for_each
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions

#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++

extern "C" {
    SEXP S_graph_deg0_lowess(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_bandwidth,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_verbose );
}

/**
 * @brief Performs degree 0 LOWESS (locally weighted average) on graph data
 *
 * @details This function implements a simplified version of LOWESS for graph data
 * where only degree 0 local models (weighted averages) are fit. The algorithm:
 * 1. For each vertex, identifies neighboring vertices within a specified bandwidth
 * 2. Computes weighted averages of response values using kernel weights
 * 3. Returns smoothed predictions
 *
 * @param y Response values at each vertex in the graph
 * @param bandwidth The fixed bandwidth (radius) to use for all local neighborhoods
 * @param kernel_type Type of kernel function for weighting (e.g., Gaussian, triangular)
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param verbose Whether to print progress information
 *
 * @return Vector of smoothed predictions for each vertex
 */
std::vector<double> set_wgraph_t::graph_deg0_lowess(
    const std::vector<double>& y,
    double bandwidth,
    size_t kernel_type,
    double dist_normalization_factor,
    bool verbose) const {

    auto start_time = std::chrono::steady_clock::now();

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    size_t n_vertices = adjacency_list.size();

    if (verbose) {
        Rprintf("Starting graph_deg0_lowess() algorithm\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Fixed bandwidth: %.4f\n", bandwidth);
    }

    // Initialize result vector
    std::vector<double> predictions(n_vertices);

    // Progress counter for parallel processing
    std::atomic<size_t> progress_counter{0};
    const size_t progress_step = std::max(size_t(1), n_vertices / 20);

    // Process each vertex
    std::vector<size_t> vertices(n_vertices);
    std::iota(vertices.begin(), vertices.end(), 0);

    // Use std::execution::par_unseq for parallel processing if available
    std::for_each(std::execution::seq, vertices.begin(), vertices.end(),
        [&](size_t vertex) {
            try {
                // Find vertices within bandwidth radius
                auto vertex_map = find_vertices_within_radius(vertex, bandwidth);

                // Distance normalization strategy:
                // 1. Find the maximum distance (max_dist) in the vertex neighborhood
                // 2. Scale max_dist by dist_normalization_factor (typically 1.1)
                // 3. Divide all distances by this scaled maximum
                // This approach ensures all normalized distances fall within [0, 1/dist_normalization_factor],
                // which is approximately [0, 0.91] when dist_normalization_factor = 1.1
                // This keeps all distances within the effective support of the kernel functions,
                // as most kernels in this implementation have support on [-1, 1]

                size_t nbhd_size = vertex_map.size();
                std::vector<double> normalized_dists(nbhd_size);
                std::vector<size_t> neighbors(nbhd_size);
                size_t counter = 0;
                for (const auto& [neighbor, distance] : vertex_map) {
                    // Normalize distance and compute kernel weight
                    normalized_dists[counter] = distance;
                    neighbors[counter++]      = neighbor;
                }

                double max_dist = 0.0;
                for (size_t k = 0; k < nbhd_size; ++k) {
                    max_dist = std::max(max_dist, normalized_dists[k]);
                }

                if (max_dist == 0) max_dist = 1;
                max_dist *= dist_normalization_factor;

                for (size_t k = 0; k < nbhd_size; ++k) {
                    normalized_dists[k] /= max_dist;
                }

                std::vector<double> weights(nbhd_size);
                kernel_fn(normalized_dists.data(), nbhd_size, weights.data());

                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                for (size_t k = 0; k < nbhd_size; ++k) {
                    weighted_sum += weights[k] * y[neighbors[k]];
                    weight_sum += weights[k];
                }

                // Compute weighted average
                if (weight_sum > 0) {
                    predictions[vertex] = weighted_sum / weight_sum;
                } else {
                    predictions[vertex] = std::numeric_limits<double>::quiet_NaN();
                }

                // Update progress counter
                if (verbose) {
                    size_t current = ++progress_counter;
                    if (current % progress_step == 0 || current == n_vertices) {
                        double percentage = 100.0 * current / n_vertices;
                        Rprintf("\rProcessing vertices: %.1f%% complete (%zu/%zu)",
                                percentage, current, n_vertices);
                        R_FlushConsole();
                    }
                }
            } catch (const std::exception& e) {
                REPORT_ERROR("Error processing vertex %zu: %s", vertex, e.what());
            }
        }
    );

    if (verbose) {
        Rprintf("\nCompleted graph_deg0_lowess in %.2f seconds\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count());
    }

    return predictions;
}

/**
 * @brief R interface for graph degree 0 LOWESS
 *
 * @param s_adj_list Graph adjacency list
 * @param s_weight_list Graph weight list
 * @param s_y Response values
 * @param s_bandwidth Fixed bandwidth for all vertices
 * @param s_kernel_type Kernel function type
 * @param s_dist_normalization_factor Distance normalization factor
 * @param s_verbose Whether to display progress
 *
 * @return A vector of smoothed predictions
 */
SEXP S_graph_deg0_lowess(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_bandwidth,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_verbose) {

    // Convert input parameters
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Get scalar parameters
    double bandwidth = REAL(s_bandwidth)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool verbose = LOGICAL(s_verbose)[0];

    // Create graph and run algorithm
    set_wgraph_t graph = set_wgraph_t(adj_list, weight_list);

    std::vector<double> predictions = graph.graph_deg0_lowess(
        y,
        bandwidth,
        kernel_type,
        dist_normalization_factor,
        verbose
    );

    // Convert result to R vector
    SEXP result = PROTECT(Rf_allocVector(REALSXP, predictions.size()));
    double* result_ptr = REAL(result);
    std::copy(predictions.begin(), predictions.end(), result_ptr);

    UNPROTECT(1);
    return result;
}
