/**
 * @brief Computes graph path linear models
 *
 * @details This function fits local linear models at each vertex of the graph using paths
 * of length h that contain the vertex. For each vertex, it:
 * 1. Identifies paths of length h containing the vertex
 * 2. Selects paths where the vertex is closest to the midpoint
 * 3. Fits weighted linear models along these paths
 * 4. Averages predictions from multiple paths if they exist
 *
 * The weighting scheme combines:
 * - Optional input weights for each vertex
 * - Kernel weights based on distance along paths
 * - Normalization to handle different path scales
 *
 * @param path_graph Graph structure containing path and connectivity information
 * @param y Response values at each vertex
 * @param weights Optional vertex weights (must be empty or match y length)
 * @param h Path length for local linear models (must be odd: h = 2s + 1)
 * @param ikernel Integer specifying the kernel type for weight computation
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 *
 * @return std::vector<double> Predicted values at each vertex based on local linear models
 *
 * @throws Rf_error if:
 * - h is not positive
 * - h is not odd
 * - ikernel is negative
 * - y.size() doesn't match number of vertices
 * - weights is non-empty and weights.size() doesn't match number of vertices
 * - dist_normalization_factor is not positive
 * - no valid paths are found for any vertex
 * - internal linear model fitting fails
 *
 * @pre h must be positive odd integer (h = 2s + 1 for some s ≥ 0)
 * @pre y.size() must equal path_graph.vertex_paths.size()
 * @pre weights must be empty or weights.size() must equal y.size()
 * @pre ikernel must be non-negative
 * @pre dist_normalization_factor must be positive
 *
 * @note Time Complexity: O(V * P * h) where:
 *       V is number of vertices
 *       P is average number of h-length paths per vertex
 *       h is path length
 *
 * @note The function implements a weighted local linear regression where:
 * - Each vertex prediction is based on paths of length h containing it
 * - Only paths where the vertex is closest to midpoint are used
 * - Weights combine user-provided weights and kernel-based distance weights
 * - Distances along paths are normalized by maximum path distance
 *
 * Example usage:
 * @code
 *     path_graph_llm_t graph = create_path_graph_llm(adj_list, weight_list, h);
 *     std::vector<double> y = {1.0, 2.0, 3.0}; // response values
 *     std::vector<double> weights = {}; // empty for uniform weights
 *     int h = 3; // path length
 *     int ikernel = 1; // kernel type
 *     auto predictions = graph_kpath_lm(graph, y, weights, h, ikernel);
 * @endcode
 *
 * @see path_graph_llm_t
 * @see predict_lm_1d
 * @see initialize_kernel
 * @see kernel_fn
 */
std::vector<double> slow_graph_kpath_lm(const path_graph_llm_t& path_graph,
                                   const std::vector<double>& y,
                                   const std::vector<double>& weights, // this vector can either be empty or has to have the same number of elements as y
                                   int ikernel,
                                   double dist_normalization_factor = 1.01) {

    int h = path_graph.h;
    if (h <= 0) {
        Rf_error("Path length h must be positive");
    }
    if (h % 2 == 0) { // Precondition: h must be odd (h = 2s + 1 for some integer s)
        Rf_error("Path length h must be odd");
    }
    if (ikernel < 0) {
        Rf_error("Invalid kernel type");
    }

    int n_vertices = path_graph.vertex_paths.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of y must match number of vertices");
    }
    if (!weights.empty() && weights.size() != n_vertices) {
        Rf_error("Length of weights must be 0 or match number of vertices");
    }
    if (dist_normalization_factor <= 0) {
        Rf_error("Distance normalization factor must be positive");
    }

    std::vector<double> Ey(n_vertices); // output vector of the mean local linear models

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

        // Extracting full paths of length h, with the position of the given vertex within the path
        std::vector<std::pair<std::vector<int>, int>> vertex_paths = vertex_path_graph_info.get_full_paths(h, path_graph.shortest_paths);

        int min_dist_to_mid_pt = h;
        std::unordered_set<int> indices_of_min_dist_to_mid_pt;
        int mid_pt = (h - 1) / 2;
        int dist_to_mid_pt;
        int n_vertex_paths = vertex_paths.size();

        for (int path_i = 0; path_i < n_vertex_paths; ++path_i) {
            if ((dist_to_mid_pt = std::abs(vertex_paths[path_i].second - mid_pt)) < min_dist_to_mid_pt) {
                min_dist_to_mid_pt = dist_to_mid_pt;
                indices_of_min_dist_to_mid_pt.clear();
                indices_of_min_dist_to_mid_pt.insert(path_i);
            } else if (dist_to_mid_pt == min_dist_to_mid_pt) {
                indices_of_min_dist_to_mid_pt.insert(path_i);
            }
        }

        if (indices_of_min_dist_to_mid_pt.empty()) {
            Rf_error("No valid paths were found");
            // Ey[vertex_i] = std::numeric_limits<double>::quiet_NaN();
            // continue;
        }

        double vertex_Ey = 0;
        for (const auto& path_i : indices_of_min_dist_to_mid_pt) {

            auto ends = vertex_path_graph_info.containing_paths[path_i];            // end points of the path contating the given vertex
            std::vector<int> path = path_graph.shortest_paths.at(ends);
            int path_n_vertices = path.size();
            int position_in_path = vertex_path_graph_info.position_in_path[path_i]; // position of the given vertex, vertex_i, within the given path with index path_i

            // computing x, y and w for a weighted linear 1d model
            std::vector<double> y_path(path_n_vertices); // y restricted to the path
            std::vector<double> x_path(path_n_vertices); // distance from the initial vertex of the path to each vertex within the path
            std::vector<double> d_path(path_n_vertices); // distance from the reference vertex to any other vertex - used in computing kernel weights
            std::vector<double> w_path(path_n_vertices, 1.0); // weights restricted to the path
            x_path[0] = 0;
            for (int i = 0; i < path_n_vertices; ++i) {
                y_path[i] = y[path[i]];
                if (i > 0) {
                    auto neighbors = path_graph.adj_list[path[i - 1]]; // neighbors of path[i - 1]
                    auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                    // finding path[i] in neighbors
                    auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                    int neighbor_i = it - neighbors.begin(); // index of path[i] in neighbors
                    x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i]; // neighbor_weights[neighbor_i] is the weight/distance from path[i-1] to path[i]
                }
            }

            if (weights.size() == n_vertices) {
                for (int i = 0; i < path_n_vertices; ++i)
                    w_path[i] = weights[path[i]];
            }

            std::vector<double> kernel_w_path(path_n_vertices); // kernel weights; these weights will be adjust by input weights if the input weights vector is non-empty
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[position_in_path]);
                if (d_path[i] > max_dist)
                    max_dist = d_path[i];
            }
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;
            for (int i = 0; i < path_n_vertices; ++i)
                d_path[i] /= max_dist;


            // creating kernel weights
            initialize_kernel(ikernel);
            kernel_fn(d_path.data(), path_n_vertices, kernel_w_path.data());

            for (int i = 0; i < path_n_vertices; ++i)
                w_path[i] *= kernel_w_path[i];

            vertex_Ey += predict_lm_1d(y_path, x_path, w_path, position_in_path);
        }

        Ey[vertex_i] = vertex_Ey / indices_of_min_dist_to_mid_pt.size();

    }

    return Ey;
}
