#include <queue>                       // For std::priority_queue
#include <limits>                      // For std::numeric_limits, INFINITY
#include <algorithm>                   // For std::max_element
#include <numeric>                     // For std::accumulate
#include <Eigen/Dense>                 // For Eigen matrix operations


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
 * @brief Finds the minimum radius required to include at least domain_min_size vertices
 *
 * @param vertex The center vertex from which to measure distances
 * @param lower_bound Initial lower bound for binary search
 * @param upper_bound Initial upper bound for binary search
 * @param domain_min_size Minimum number of vertices required in the neighborhood
 * @param precision Precision threshold for terminating binary search
 *
 * @return double The minimum radius that ensures at least domain_min_size vertices
 *                are within the neighborhood of the vertex
 *
 * @details This function performs a binary search to find the smallest radius value
 * for which the vertex neighborhood (set of all vertices within the radius) contains
 * at least domain_min_size elements. This is required to ensure sufficient data
 * points for linear model fitting.
 */
double set_wgraph_t::find_minimum_radius_for_domain_min_size(
    size_t vertex,
    double lower_bound,
    double upper_bound,
    size_t domain_min_size,
    double precision) const {

    // Check if the lower bound already satisfies the condition
    std::unordered_map<size_t, double> neighborhood = find_vertices_within_radius(vertex, lower_bound);
    if (neighborhood.size() >= domain_min_size) {
        return lower_bound;
    }

    // Check if the upper bound fails to satisfy the condition
    neighborhood = find_vertices_within_radius(vertex, upper_bound);
    if (neighborhood.size() < domain_min_size) {
        Rprintf("\n---------------------\nERROR\nvertex: %zu\n"
                "Not enough vertices (found: %zu, needed: %zu) within upper_bound: %.4f\n",
                vertex + 1, neighborhood.size(), domain_min_size, upper_bound);

        // Additional debugging similar to the original function could be added here

        REPORT_ERROR("ERROR: Insufficient vertices in maximum radius\n");
    }

    // Binary search for the minimum radius
    while ((upper_bound - lower_bound) > precision * graph_diameter) {
        double mid = (lower_bound + upper_bound) / 2.0;
        neighborhood = find_vertices_within_radius(vertex, mid);

        if (neighborhood.size() >= domain_min_size) {
            // If condition is satisfied, try a smaller radius
            upper_bound = mid;
        } else {
            // If condition is not satisfied, try a larger radius
            lower_bound = mid;
        }
    }

    // Return the upper bound as the minimum radius that satisfies the condition
    return upper_bound;
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
    size_t n_vertices = vertex_map.size();
    Eigen::MatrixXd embedding(n_vertices, n_evectors + 1);

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

