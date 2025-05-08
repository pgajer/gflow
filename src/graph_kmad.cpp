#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>     // for std::vector
#include <algorithm>  // for std::sort, std::max
#include <stdexcept> // for std::invalid_argument, std::runtime_error
#include <cmath>     // for std::abs
#include <utility>   // for std::pair (used in weighted version only)

#include "SEXP_cpp_conversion_utils.hpp"
#include "kernels.h"

extern "C" {
    SEXP S_graph_mad(SEXP s_graph, SEXP s_y);
}

/**
* @brief Computes Kernel-Weighted Median Absolute Deviation (MAD) on a Graph
*
* @details
* This function calculates the kernel-weighted median absolute deviation of values associated
* with vertices in a graph. The MAD is a robust measure of statistical dispersion that is more
* resilient to outliers than standard deviation.
*
* For each vertex, the function:
* 1. Collects the vertex's value and its neighbors' values with their kernel weights
* 2. Computes the weighted median of these values
* 3. Calculates absolute deviations from this median
* 4. Computes the weighted median of these deviations
*
* The kernel weights are determined by:
* - Normalizing edge lengths by the maximum distance in the local neighborhood
* - Applying a specified kernel function to these normalized distances
*
* For vertices with no neighbors, their MAD is set to 0 as there is no variation to measure.
*
* @param graph A list of integer vectors where each element represents a vertex and contains
*        indices of its neighboring vertices
* @param edge_lengths A list of numeric vectors containing the lengths of edges to neighbors.
*        Must match the structure of the `graph` parameter
* @param y A numeric vector of values associated with each vertex in the graph
* @param ikernel An integer specifying the kernel function to use for weighting
* @param dist_normalization_factor A numeric value used to normalize distances in the graph.
*        Default is 1.01
*
* @return A numeric vector containing the kernel-weighted MAD for each vertex in the graph
*
* @throws std::invalid_argument If input vectors have inconsistent sizes
* @throws std::runtime_error If kernel initialization fails
*
* @note The function uses a weighted median calculation which preserves the influence of
*       kernel weights on the final MAD value
*
* @see initialize_kernel(), kernel_fn()
*/
std::vector<double> graph_kernel_weighted_mad(const std::vector<std::vector<int>>& graph,
                                            const std::vector<std::vector<double>>& edge_lengths,
                                            const std::vector<double>& y,
                                            int ikernel,
                                            double dist_normalization_factor = 1.01) {
    // Validate input
    if (graph.size() != y.size() || graph.size() != edge_lengths.size()) {
        Rf_error("Input vectors must have consistent sizes");
    }

    auto kmad = std::vector<double>(y.size(), 0.0);
    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Determine the maximum number of neighbors
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Initialize working vectors
    std::vector<double> kernel_weights(max_neighbors + 1);
    std::vector<double> distances(max_neighbors + 1);
    std::vector<double> weighted_values(max_neighbors + 1);
    std::vector<double> weights(max_neighbors + 1);
    std::vector<double> deviations(max_neighbors + 1);

    // Helper function to compute weighted median
    auto weighted_median = [](const std::vector<double>& values,
                            const std::vector<double>& weights,
                            size_t n) -> double {
        // Create pairs of values and weights
        std::vector<std::pair<double, double>> paired_data(n);
        for (size_t i = 0; i < n; ++i) {
            paired_data[i] = {values[i], weights[i]};
        }

        // Sort by values
        std::sort(paired_data.begin(), paired_data.begin() + n,
                 [](const auto& a, const auto& b) { return a.first < b.first; });

        // Find weighted median
        double total_weight = 0.0;
        for (size_t i = 0; i < n; ++i) {
            total_weight += paired_data[i].second;
        }

        double weight_sum = 0.0;
        double half_weight = total_weight / 2.0;

        for (size_t i = 0; i < n; ++i) {
            weight_sum += paired_data[i].second;
            if (weight_sum >= half_weight) {
                return paired_data[i].first;
            }
        }
        return paired_data[0].first; // Fallback for edge cases
    };

    // Process each vertex
    for (size_t i = 0; i < graph.size(); ++i) {
        if (!graph[i].empty()) {
            // Compute kernel weights
            distances[0] = 0;  // Distance to self
            double max_dist = 0.0;
            for (size_t j = 0; j < graph[i].size(); ++j) {
                distances[j + 1] = edge_lengths[i][j];
                max_dist = std::max(max_dist, distances[j + 1]);
            }

            max_dist = (max_dist == 0) ? 1 : max_dist * dist_normalization_factor;

            for (size_t j = 0; j < graph[i].size() + 1; ++j) {
                distances[j] /= max_dist;
            }

            size_t n = graph[i].size() + 1;
            kernel_fn(distances.data(), n, kernel_weights.data());

            // Collect values and weights
            weighted_values[0] = y[i];
            weights[0] = kernel_weights[0];
            for (size_t j = 0; j < graph[i].size(); ++j) {
                weighted_values[j + 1] = y[graph[i][j]];
                weights[j + 1] = kernel_weights[j + 1];
            }

            // Compute weighted median
            double median = weighted_median(weighted_values, weights, n);

            // Compute absolute deviations
            for (size_t j = 0; j < n; ++j) {
                deviations[j] = std::abs(weighted_values[j] - median);
            }

            // Compute MAD as weighted median of absolute deviations
            kmad[i] = weighted_median(deviations, weights, n);
        }
        // For vertices with no neighbors, MAD remains 0
    }

    return kmad;
}


/**
* @brief Computes Unweighted Median Absolute Deviation (MAD) on a Graph
*
* @details
* This function calculates the unweighted median absolute deviation of values associated
* with vertices in a graph. For each vertex, it considers its direct neighbors in the graph
* topology without accounting for edge weights or distances.
*
* The MAD computation for each vertex follows these steps:
* 1. Collects values from the vertex and its immediate neighbors
* 2. Computes the median of these values
* 3. Calculates absolute deviations from this median
* 4. Computes the median of these absolute deviations
*
* All vertices in each neighborhood (including the central vertex) are treated with equal
* importance, making this a purely topological measure of local variation.
*
* For vertices with no neighbors:
* - Their MAD value is set to 0 as there is no variation to measure
* - This reflects the inability to compute dispersion for a single point
*
* @param graph A list of integer vectors where each element represents a vertex and contains
*        indices of its neighboring vertices. The graph structure is assumed to be valid
*        (i.e., all vertex indices are within bounds)
* @param y A numeric vector of values associated with each vertex in the graph
*
* @return A numeric vector containing the unweighted MAD for each vertex in the graph
*
* @throws std::invalid_argument If the size of y doesn't match the number of vertices in the graph
*                              or if the graph contains invalid vertex indices
*
* @note This implementation is more efficient than the weighted version as it:
*       - Doesn't require kernel computations
*       - Uses simpler median calculations
*       - Requires less memory for intermediate calculations
*
* @see graph_kernel_weighted_mad() for the weighted version of this computation
*/
std::vector<double> graph_mad(const std::vector<std::vector<int>>& graph,
                              const std::vector<double>& y) {
    // Input validation
    if (graph.size() != y.size()) {
        Rf_error("Number of vertices in graph must match size of y");
    }

    // For each vertex, validate neighbor indices
    for (size_t i = 0; i < graph.size(); ++i) {
        for (int neighbor : graph[i]) {
            if (neighbor < 0 || static_cast<size_t>(neighbor) >= y.size()) {
                Rf_error("Graph contains invalid vertex indices");
            }
        }
    }

    auto mad = std::vector<double>(y.size(), 0.0);

    // Determine maximum neighborhood size for pre-allocation
    size_t max_neighbors = 0;
    for (const auto& neighbors : graph) {
        max_neighbors = std::max(max_neighbors, neighbors.size());
    }

    // Pre-allocate vectors for calculations
    std::vector<double> neighborhood_values(max_neighbors + 1);  // +1 for the vertex itself
    std::vector<double> deviations(max_neighbors + 1);

    // Helper function to compute median
    auto compute_median = [](std::vector<double>& values, size_t n) -> double {
        std::sort(values.begin(), values.begin() + n);
        if (n % 2 == 0) {
            return (values[n/2 - 1] + values[n/2]) / 2.0;
        }
        return values[n/2];
    };

    // Process each vertex
    for (size_t i = 0; i < graph.size(); ++i) {
        const auto& neighbors = graph[i];

        if (!neighbors.empty()) {
            // Collect values from neighborhood
            size_t n = neighbors.size() + 1;  // Include the vertex itself
            neighborhood_values[0] = y[i];    // Start with the vertex's value

            // Add neighbor values
            for (size_t j = 0; j < neighbors.size(); ++j) {
                neighborhood_values[j + 1] = y[neighbors[j]];
            }

            // Compute median of neighborhood
            double median = compute_median(neighborhood_values, n);

            // Compute absolute deviations
            for (size_t j = 0; j < n; ++j) {
                deviations[j] = std::abs(neighborhood_values[j] - median);
            }

            // Compute MAD as median of absolute deviations
            mad[i] = compute_median(deviations, n);
        }
        // For vertices with no neighbors, MAD remains 0
    }

    return mad;
}


/**
* @brief R Interface for Graph Kernel MAD Calculation
*
* @details
* This function serves as an interface between R and the C++ implementation of the graph
* kernel MAD calculation. It handles the conversion of R objects to C++ data structures,
* calls the core computation function, and converts the results back to R format.
*
* The function performs the following steps:
* 1. Validates input parameters
* 2. Converts R graph structure to C++ vectors
* 3. Converts R numeric vector to C++ vector
* 4. Calls the core MAD calculation function
* 5. Converts results back to R numeric vector
*
* @param s_graph An R object (SEXP) representing the graph structure. Must be convertible
*        by Rgraph_to_vector to a vector of integer vectors.
* @param s_y An R object (SEXP) representing the numeric vector of values associated
*        with vertices. Must be a numeric vector (REALSXP).
*
* @return An R numeric vector (SEXP) containing the MAD values for each vertex
*
* @throws R error if:
*         - s_y is not a numeric vector
*         - Lengths of graph and y are inconsistent
*         - Memory allocation fails
*
* @note This function handles R's garbage collection by using PROTECT/UNPROTECT
*       macros appropriately.
*
* @see graph_mad()
*/
SEXP S_graph_mad(SEXP s_graph, SEXP s_y) {
    // Input validation
    if (!isReal(s_y)) {
        error("'y' must be a numeric vector");
    }

    // Convert graph structure
    std::vector<std::vector<int>> graph = convert_adj_list_from_R(s_graph);

    // Get and validate dimensions
    int y_length = LENGTH(s_y);
    if (static_cast<size_t>(y_length) != graph.size()) {
        error("Length of 'y' must match number of vertices in graph");
    }

    // Convert y vector
    std::vector<double> y(REAL(s_y), REAL(s_y) + y_length);

    // Compute MAD
    std::vector<double> result;
    try {
        result = graph_mad(graph, y);
    } catch (const std::exception& e) {
        error("Error in MAD calculation: %s", e.what());
    }

    // Convert result to R vector
    SEXP s_result = PROTECT(allocVector(REALSXP, result.size()));
    if (s_result == R_NilValue) {
        UNPROTECT(1);
        error("Memory allocation failed");
    }

    // Copy results
    std::copy(result.begin(), result.end(), REAL(s_result));

    UNPROTECT(1);
    return s_result;
}
