#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length

#include <Eigen/Dense>                 // For Eigen matrix operations
#include <queue>                       // For std::priority_queue
#include <limits>                      // For std::numeric_limits, INFINITY
#include <algorithm>                   // For std::max_element
#include <numeric>                     // For std::accumulate
#include <map>                         // For std::map

#include "set_wgraph.hpp"
#include "error_utils.h"               // For REPORT_ERROR

/**
 * @brief Finds all vertices within a specified radius of a reference vertex
 *
 * @param vertex The reference vertex from which distances are measured
 * @param radius Maximum distance threshold for inclusion
 * @return std::unordered_map<size_t, double> Map of vertices to their distances from the reference vertex
 *
 * @details This function computes the shortest path distance from the reference vertex
 * to all other vertices in the graph, and returns those vertices whose distance
 * is less than or equal to the specified radius along with their distances.
 * The reference vertex is always included in the result with distance 0.
 */
std::unordered_map<size_t, double> set_wgraph_t::find_vertices_within_radius(
    size_t vertex,
    double radius) const {

    // Result will map vertices to their distances from the reference vertex
    std::unordered_map<size_t, double> result;

    // Always include the reference vertex with distance 0
    result[vertex] = 0.0;

    // Initialize priority queue for Dijkstra's algorithm
    // Using min-heap with pairs of (distance, vertex)
    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<>> pq;

    // Track distances
    std::vector<double> distances(adjacency_list.size(), INFINITY);
    distances[vertex] = 0.0;

    // Initialize with reference vertex
    pq.push({0.0, vertex});

    // Run Dijkstra's algorithm
    while (!pq.empty()) {
        auto [current_dist, current_vertex] = pq.top();
        pq.pop();

        // If we've exceeded the radius, we can terminate early
        if (current_dist > radius) {
            break;
        }

        // Skip if we've already found a shorter path
        if (current_dist > distances[current_vertex]) {
            continue;
        }

        // Add this vertex to the result with its distance
        result[current_vertex] = current_dist;

        // Explore neighbors
        for (const auto& edge : adjacency_list[current_vertex]) {
            size_t neighbor = edge.vertex;
            double new_dist = current_dist + edge.weight;

            // Update distance if we found a shorter path
            if (new_dist < distances[neighbor]) {
                distances[neighbor] = new_dist;
                pq.push({new_dist, neighbor});
            }
        }
    }

    return result;
}

/**
 * @brief Creates a spectral embedding of vertices using Laplacian eigenvectors
 *
 * @param vertex_map Map of vertices to their distances from a reference vertex
 * @param eigenvectors Matrix of Laplacian eigenvectors (each column is an eigenvector)
 * @param n_evectors Number of eigenvectors to use for the embedding
 *
 * @return Eigen::MatrixXd Matrix where each row corresponds to a vertex's coordinates in the embedding space
 *
 * @details This function takes a set of vertices and creates a lower-dimensional embedding using
 * the first n_evectors eigenvectors of the graph Laplacian. The resulting embedding places
 * vertices in a space where Euclidean distances approximate geodesic distances in the
 * graph structure, enabling effective local linear modeling.
 *
 * The eigenvectors should be sorted by increasing eigenvalue (the first eigenvector typically
 * corresponds to the constant function and may be excluded from the embedding).
 *
 * @note The first column of the returned matrix contains a constant (1.0) to account for
 * the intercept term in subsequent linear regression models.
 */
Eigen::MatrixXd create_spectral_embedding(
    const std::map<size_t, double>& vertex_map,
    const Eigen::MatrixXd& eigenvectors,
    size_t n_evectors) {

    // Ensure n_evectors is valid
    size_t n_available = eigenvectors.cols();
    if (n_evectors > n_available) {
        n_evectors = n_available;
        Rprintf("Warning: Requested %zu eigenvectors but only %zu are available. Using %zu.\n",
                n_evectors, n_available, n_evectors);
    }

    // Initialize the embedding matrix with one extra column for the intercept
    //size_t n_vertices = vertex_map.size();
    Eigen::MatrixXd embedding(vertex_map.size(), n_evectors + 1);

    // Set first column to 1.0 for the intercept term
    embedding.col(0).setOnes();

    // Fill in the embedding coordinates from eigenvectors
    // We'll skip the first eigenvector (constant function) and use the next n_evectors
    size_t row_idx = 0;
    for (const auto& [vertex_idx, distance] : vertex_map) {
        // Vertex coordinates in the embedding space are the corresponding entries in the eigenvectors
        for (size_t j = 0; j < n_evectors; ++j) {
            // Add 1 to j to skip first eigenvector
            embedding(row_idx, j + 1) = eigenvectors(vertex_idx, j + 1); 
        }
        row_idx++;
    }

    return embedding;
}
