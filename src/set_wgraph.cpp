#include <R.h>
#include <Rinternals.h>

#include <queue>
#include <random>     // for std::mt19937
#include <sstream>
#include <unordered_set>
#include <functional>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>

#include "set_wgraph.hpp"
#include "error_utils.h" // for REPORT_ERROR()
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
    SEXP S_compute_edge_weight_deviations(SEXP s_adj_list, SEXP s_weight_list);
    SEXP S_compute_edge_weight_rel_deviations(SEXP s_adj_list, SEXP s_weight_list);
    SEXP S_remove_redundant_edges(SEXP s_adj_list, SEXP s_weight_list);
    SEXP S_find_graph_paths_within_radius(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_start,
        SEXP s_radius);
}

void set_wgraph_t::compute_graph_diameter() {
    // Estimate the diameter of the graph
    auto [end1, diam] = get_vertex_eccentricity(0);  // Start from vertex 0
    auto [end2, diameter] = get_vertex_eccentricity(end1);
    graph_diameter = diameter;
}

/**
 * @brief Creates a subgraph containing only the specified vertices
 *
 * This method creates a new graph that includes only the vertices specified in the input vector.
 * The vertex indices in the new graph are renumbered sequentially starting from 0.
 * Edge weights between the remaining vertices are preserved.
 *
 * @param vertices Vector of vertex indices to include in the subgraph
 * @return set_wgraph_t A new graph containing only the specified vertices
 */
set_wgraph_t set_wgraph_t::create_subgraph(const std::vector<size_t>& vertices) const {
    // Validate input
    for (size_t v : vertices) {
        if (v >= adjacency_list.size()) {
            throw std::out_of_range("Vertex index out of range");
        }
    }

    // Create a mapping from old indices to new indices
    std::unordered_map<size_t, size_t> old_to_new_index;
    for (size_t i = 0; i < vertices.size(); ++i) {
        old_to_new_index[vertices[i]] = i;
    }

    // Create a new graph with the appropriate number of vertices
    set_wgraph_t subgraph;
    subgraph.adjacency_list.resize(vertices.size());

    // Copy edges between vertices that are in the subgraph
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t old_vertex = vertices[i];

        // Iterate through all edges of the original vertex
        for (const auto& edge : adjacency_list[old_vertex]) {
            // Check if the target vertex is also in the subgraph
            auto it = old_to_new_index.find(edge.vertex);
            if (it != old_to_new_index.end()) {
                // Add the edge with the new vertex index
                subgraph.adjacency_list[i].insert({it->second, edge.weight});
            }
        }
    }

    // Copy other relevant properties of the graph
    // (assuming there are properties like graph_diameter in the actual implementation)
    if (adjacency_list.size() > 0 && subgraph.adjacency_list.size() > 0) {

        subgraph.max_packing_radius = -1.0;

        // Find diameter endpoints
        auto [end1, diam] = subgraph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = subgraph.get_vertex_eccentricity(end1);
        subgraph.graph_diameter = diameter;
    }

    return subgraph;
}


std::vector<int> union_find(const std::vector<std::vector<int>>& adj_vect);

/**
 * @brief Counts the number of connected components in the graph
 *
 * @details Uses the union-find algorithm to identify connected components.
 * Two vertices are in the same connected component if there exists a path between them.
 *
 * Time complexity: O(α(V)×E) where α is the inverse Ackermann function, V is number of vertices,
 * and E is number of edges.
 * Space complexity: O(V)
 *
 * @return The number of connected components in the graph
 */
size_t set_wgraph_t::count_connected_components() const {
    // Get the number of vertices in the graph
    const size_t n_vertices = adjacency_list.size();

    // Convert the set_wgraph_t adjacency list to the format expected by union_find
    std::vector<std::vector<int>> adj_vect(n_vertices);

    // Reserve space for each adjacency list to avoid reallocations
    for (size_t i = 0; i < n_vertices; ++i) {
        adj_vect[i].reserve(adjacency_list[i].size());
    }

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            adj_vect[i].push_back(static_cast<int>(edge.vertex));
        }
    }

    // Apply union_find to get component identifiers for each vertex
    std::vector<int> components = union_find(adj_vect);

    // Count unique component identifiers using a set
    std::unordered_set<int> unique_components;
    for (int component_id : components) {
        unique_components.insert(component_id);
    }

    return unique_components.size();
}

/**
 * @brief Get all connected components in the graph
 *
 * This method identifies all connected components in the graph and returns
 * a vector of vectors, where each inner vector contains the vertex indices
 * belonging to one connected component.
 *
 * @return std::vector<std::vector<size_t>> Vector of connected components where each component
 *                                         is a vector of vertex indices
 */
std::vector<std::vector<size_t>> set_wgraph_t::get_connected_components() const {
    // Get the number of vertices in the graph
    const size_t n_vertices = adjacency_list.size();

    // Convert the set_wgraph_t adjacency list to the format expected by union_find
    std::vector<std::vector<int>> adj_vect(n_vertices);

    // Reserve space for each adjacency list to avoid reallocations
    for (size_t i = 0; i < n_vertices; ++i) {
        adj_vect[i].reserve(adjacency_list[i].size());
    }

    // Fill the adjacency vectors
    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            adj_vect[i].push_back(static_cast<int>(edge.vertex));
        }
    }

    // Apply union_find to get component identifiers for each vertex
    std::vector<int> component_ids = union_find(adj_vect);

    // Create a map from component_id to list of vertices in that component
    std::unordered_map<int, std::vector<size_t>> component_map;

    // Group vertices by their component identifiers
    for (size_t vertex = 0; vertex < component_ids.size(); ++vertex) {
        int component_id = component_ids[vertex];
        component_map[component_id].push_back(vertex);
    }

    // Convert the map to a vector of vectors
    std::vector<std::vector<size_t>> components;
    components.reserve(component_map.size());

    for (const auto& [component_id, vertices] : component_map) {
        components.push_back(vertices);
    }

    return components;
}


/**
 * @brief Prints the weighted graph structure to the R console
 *
 * @details This function prints a detailed representation of the weighted graph including:
 * - Total number of vertices
 * - Total number of edges
 * - Average vertex degree
 * - Complete edge list with weights
 *
 * For each edge, it prints a tuple (neighbor_vertex, weight) where:
 * - neighbor_vertex is the index of the neighboring vertex (adjusted by vertex_index_shift)
 * - weight is the weight (length) of the edge
 *
 * @param vertex_index_shift Integer value added to all vertex indices in the output
 *        (useful for converting between 0-based and 1-based indexing systems)
 * @param name Optional name identifier for the graph
 *
 * @note The function uses Rprintf for output, making it suitable for R package integration
 * @note Edge counts are divided by 2 as each edge is stored twice in the undirected graph
 */
void vect_wgraph_t::print(const std::string& name,
                          size_t vertex_index_shift) const {

    if (adjacency_list.empty()) {
        Rprintf("\nEmpty graph\n");
        return;
    }

    // Calculate total edges
    size_t total_edges = 0;
    for (const auto& vertex_edges : adjacency_list) {
        total_edges += vertex_edges.size();
    }
    total_edges /= 2;  // Each edge is counted twice in undirected graph

    // Print graph summary
    //Rprintf("\nWeighted graph");
    if (!name.empty()) {
        Rprintf("\n\nvect_wgraph_t: '%s'\n", name.c_str());
    } else {
        Rprintf("\n");
    }

    Rprintf("Number of vertices: %zu\n", adjacency_list.size());
    Rprintf("Number of edges: %zu\n", total_edges);
    Rprintf("Average degree: %.2f\n",
            adjacency_list.empty() ? 0.0 : (2.0 * total_edges / adjacency_list.size()));

    // Print edge list
    Rprintf("Edge List:\n");
    for (size_t vertex_id = 0; vertex_id < adjacency_list.size(); ++vertex_id) {
        Rprintf("Vertex %zu:\n", vertex_id + vertex_index_shift);
        for (const auto& edge : adjacency_list[vertex_id]) {
            Rprintf("\t(%zu, %.3f)\n",
                   edge.vertex + vertex_index_shift,
                   edge.weight);
        }
    }
}


/**
 * @brief Constructs a weighted graph from adjacency and weight lists
 *
 * Creates a graph representation where each vertex's edges are stored as a set of edge_info_t
 * structures in the adjacency_list. The sets ensure no duplicate edges are stored for each vertex.
 *
 * @param adj_list Vector of adjacency lists where adj_list[i] contains the indices of vertices
 *                 connected to vertex i
 * @param weight_list Vector of weight lists where weight_list[i][j] contains the weight of the edge
 *                    from vertex i to vertex adj_list[i][j]
 *
 * @throws Rf_error if:
 *         - adj_list and weight_list have different sizes
 *         - For any vertex i, adj_list[i] and weight_list[i] have different sizes
 *         - Any target vertex index is negative or >= number of vertices
 *
 * @note The constructor assumes 0-based indexing for vertices
 *
 * @example
 * std::vector<std::vector<int>> adj = {{1, 2}, {0}, {0}};     // Vertex 0 connects to 1,2; 1 to 0; 2 to 0
 * std::vector<std::vector<double>> w = {{1.0, 2.0}, {1.0}, {2.0}}; // Corresponding edge weights
 * set_wgraph_t graph(adj, w);
 */
set_wgraph_t::set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                           const std::vector<std::vector<double>>& weight_list) {

    graph_diameter = -1.0;
    max_packing_radius = -1.0;

    // Validate input sizes
    if (adj_list.size() != weight_list.size()) {
        REPORT_ERROR("Adjacency and weight lists must have the same size");
    }

    // Initialize adjacency list with the correct number of vertices
    adjacency_list.resize(adj_list.size());

    // Populate the graph
    for (size_t i = 0; i < adj_list.size(); ++i) {
        if (adj_list[i].size() != weight_list[i].size()) {
            REPORT_ERROR("Mismatch between adjacency and weight list sizes at vertex %d", i);
        }

        // Add each edge
        for (size_t j = 0; j < adj_list[i].size(); ++j) {
            size_t target = static_cast<size_t>(adj_list[i][j]);
            double weight = weight_list[i][j];

            // Validate vertex index
            if (target >= adj_list.size()) {
                REPORT_ERROR("ERROR: Invalid index of neighbor vertex %zu at vertex %zu.\nMake sure adj_list is modified to be 0-based!", target, i);
            }

            // Create and insert edge
            edge_info_t edge{target, weight};
            adjacency_list[i].insert(edge);
        }
    }
}

/**
 * Constructor that builds a weighted graph from an iknn_graph_t.
 * Each vertex in the iknn_graph corresponds to a vertex in this graph.
 * Each nearest neighbor in the iknn_graph is converted to an edge in this graph,
 * with the target vertex being the neighbor's index and the weight being the distance.
 *
 * @param iknn_graph The input k-nearest neighbors graph
 */
set_wgraph_t::set_wgraph_t(const iknn_graph_t& iknn_graph) {

    graph_diameter = -1.0;
    max_packing_radius = -1.0;

    // Initialize the adjacency list with the same number of vertices as the input graph
    adjacency_list.resize(iknn_graph.graph.size());

    // Process each vertex in the input graph
    for (size_t source_vertex = 0; source_vertex < iknn_graph.graph.size(); ++source_vertex) {
        // Process each neighbor of the current vertex
        for (const auto& neighbor : iknn_graph.graph[source_vertex]) {
            // Create an edge from the source vertex to the neighbor
            edge_info_t edge;
            edge.vertex = neighbor.index;  // Target vertex is the neighbor's index
            edge.weight = neighbor.dist;   // Weight is the distance to the neighbor

            // Add the edge to the adjacency list of the source vertex
            adjacency_list[source_vertex].insert(edge);
        }
    }
}

/**
 * @brief Prints the graph structure
 *
 * @param split If true, prints adjacency and weight lists separately
 * @param shift Optional vertex ID offset for printing
 *
 * The shift parameter is useful when vertex IDs need to be adjusted for display
 * (e.g., for 1-based indexing in R).
 */
void set_wgraph_t::print(const std::string& name = "",
                         bool split = false,
                         size_t shift = 0) const {

    int n_edges = 0;
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        n_edges += adjacency_list[i].size();
    }

    if (!name.empty()) {
        Rprintf("\n\nset_wgraph_t: '%s'\n", name.c_str());
    } else {
        Rprintf("\n");
    }

    Rprintf("\nGraph Structure:\n");
    Rprintf("Number of vertices: %zu\n", adjacency_list.size());
    Rprintf("Number of edges: %d\n", n_edges);

    if (split) {
        std::vector<std::vector<int>> adj_list(adjacency_list.size());
        std::vector<std::vector<double>> weight_list(adjacency_list.size());

        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                adj_list[i].push_back(edge.vertex + shift);
                weight_list[i].push_back(edge.weight);
            }
        }

        Rprintf("Adjacency Lists:\n");
        for (size_t i = 0; i < adj_list.size(); ++i) {
            Rprintf("%zu: ", i + shift);
            for (size_t j = 0; j < adj_list[i].size(); ++j) {
                Rprintf("(%d, %.6f) ", adj_list[i][j], weight_list[i][j]);
            }
            Rprintf("\n");
        }
    } else {
        Rprintf("Edge List:\n");
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                if (i < edge.vertex) { // Print each edge only once
                    Rprintf("(%zu, %zu): %.6f\n", i + shift, edge.vertex + shift, edge.weight);
                }
            }
        }
    }
}


/**
 * @brief R interface for finding shortest paths within a radius from a start vertex
 *
 * This function serves as an interface between R and the C++ implementation of
 * find_graph_paths_within_radius
 * executes the algorithm, and returns the results as R objects.
 *
 * @param s_adj_list SEXP containing the adjacency list representation of the graph
 * @param s_weight_list SEXP containing weights for each edge in the graph
 * @param s_start SEXP containing an integer specifying the starting vertex (0-based)
 * @param s_radius SEXP containing a numeric value for the maximum path distance
 *
 * @return SEXP A named list containing:
 *         - paths: A list of integer vectors, each representing a path from start vertex
 *         - reachable_vertices: An integer vector of all vertices reachable within the radius
 *         - vertex_to_path_map: A numeric matrix with three columns:
 *           * First column: vertex indices
 *           * Second column: corresponding path indices in the paths list
 *           * Third column: total path length
 *
 * @details The returned vertex_to_path_map enables O(1) lookups to find the path
 *          to any specific vertex. This is particularly useful for applications
 *          that need to quickly retrieve paths to specific destinations.
 *
 * @note This function uses 0-based vertex indexing in the C++ implementation.
 *       If working with 1-based indexing in R, users should adjust indices accordingly.
 *
 * @see find_graph_paths_within_radius (the underlying C++ implementation)
 */
// SEXP S_find_shortest_paths_within_radius(
SEXP S_find_graph_paths_within_radius(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_start,
    SEXP s_radius) {

    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    set_wgraph_t graph(adj_list, weight_list);

    size_t start_vertex = (size_t)INTEGER(s_start)[0];
    double radius = REAL(s_radius)[0];

    shortest_paths_t shortest_paths = graph.find_graph_paths_within_radius(start_vertex, radius);

    // Creating return list
    // return list component names
    const char* names[] = {
        "paths",
        "reachable_vertices",
        "vertex_to_path_map",
        NULL
    };

    size_t n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    SEXP result = PROTECT(allocVector(VECSXP, n_elements));

    SEXP result_names = PROTECT(allocVector(STRSXP, n_elements));
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, mkChar(names[i]));
    }
    setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(1); // for result_names

    // Creating paths list
    // Helper function to convert vector to SEXP
    auto convert_int_vector_to_R = [](const std::vector<size_t>& vec) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(INTSXP, vec.size()));
        int* ptr = INTEGER(r_vec);

        // std::copy(vec.begin(), vec.end(), ptr); // Direct copy without index adjustment

        for (size_t i = 0; i < vec.size(); i++) {
            ptr[i] = static_cast<int>(vec[i]) + 1;  // Convert to int and add 1
        }

        return r_vec;
    };

    size_t n_paths = shortest_paths.paths.size();
    SEXP r_paths_list = PROTECT(allocVector(VECSXP, n_paths));
    for (size_t i = 0; i < n_paths; i++) {
        SEXP r_path = convert_int_vector_to_R(shortest_paths.paths[i].vertices);
        SET_VECTOR_ELT(r_paths_list, i, r_path);
        UNPROTECT(1); // for r_path
    }
    SET_VECTOR_ELT(result, 0, r_paths_list);
    UNPROTECT(1); // for r_paths_list

    // Creating reachable_vertices vector
    auto convert_int_unordered_set_to_R = [](const std::unordered_set<int>& set) -> SEXP {
        SEXP r_vec = PROTECT(allocVector(INTSXP, set.size()));
        int* ptr = INTEGER(r_vec);  // Corrected: int* instead of double*
        //std::copy(set.begin(), set.end(), ptr);
        size_t i = 0;
        for (const auto& val : set) {
            ptr[i++] = val + 1;  // Convert to 1-based indexing
        }

        return r_vec;
    };

    SEXP r_reachable_vertices = convert_int_unordered_set_to_R(shortest_paths.reachable_vertices);
    SET_VECTOR_ELT(result, 1, r_reachable_vertices);
    UNPROTECT(1); // for r_reachable_vertices

    // Helper function to convert subpath_t map to a matrix
    auto convert_vertex_subpath_map_to_matrix = [](
        const std::unordered_map<size_t, subpath_t>& map) -> SEXP {

        if (map.empty()) return R_NilValue;
        size_t nrow = map.size();
        size_t ncol = 3;  // vertex, path_idx, vertex_idx

        SEXP r_mat = PROTECT(allocMatrix(REALSXP, nrow, ncol));
        double* ptr = REAL(r_mat);

        size_t i = 0;
        size_t nrow_times_1 = nrow;
        size_t nrow_times_2 = 2 * nrow;

        for (const auto& [vertex, subpath] : map) {
            ptr[i]                 = vertex + 1;
            ptr[i + nrow_times_1]  = subpath.path_idx + 1;
            ptr[i + nrow_times_2]  = subpath.vertex_idx + 1;
            i++;
        }

        // Add column names
        SEXP colnames = PROTECT(allocVector(STRSXP, ncol));
        SET_STRING_ELT(colnames, 0, mkChar("vertex"));
        SET_STRING_ELT(colnames, 1, mkChar("path_idx"));
        SET_STRING_ELT(colnames, 2, mkChar("vertex_idx"));

        SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);  // No row names
        SET_VECTOR_ELT(dimnames, 1, colnames);

        setAttrib(r_mat, R_DimNamesSymbol, dimnames);

        UNPROTECT(2);  // for colnames and dimnames
        return r_mat;
    };

    SEXP r_vertex_to_path_mat = convert_vertex_subpath_map_to_matrix(shortest_paths.vertex_to_path_map);
    SET_VECTOR_ELT(result, 2, r_vertex_to_path_mat);
    UNPROTECT(1); // for r_vertex_to_path_mat

    UNPROTECT(1); // for result
    return result;
}


/**
 * @brief Identifies connected components within a pro-cell of the graph
 *
 * This method performs a depth-first search to find connected components among
 * the vertices belonging to a Morse-Smale pro-cell. Two vertices are considered
 * connected if they are adjacent in the original graph and both belong to the
 * pro-cell.
 *
 * @param procell_vertices Vector of vertex indices belonging to the pro-cell
 * @param graph Reference to the original weighted graph
 * @return std::vector<std::vector<size_t>> List of connected components, where each
 *         component is represented as a vector of vertex indices
 */
std::vector<std::vector<size_t>> set_wgraph_t::find_connected_components(
    const std::vector<size_t>& procell_vertices) const {

    std::unordered_map<size_t, size_t> vertex_to_component;
    size_t current_component = 0;

    // Helper function for DFS
    std::function<void(size_t, size_t)> dfs = [&](size_t vertex, size_t component) {
        vertex_to_component[vertex] = component;

        // Get neighbors from the graph's adjacency list
        for (const auto& edge : adjacency_list[vertex]) {
            size_t neighbor = edge.vertex;
            // Only consider vertices within the pro-cell
            if (std::find(procell_vertices.begin(), procell_vertices.end(), neighbor) != procell_vertices.end()
                && vertex_to_component.find(neighbor) == vertex_to_component.end()) {
                dfs(neighbor, component);
            }
        }
    };

    // Find components using DFS
    for (size_t vertex : procell_vertices) {
        if (vertex_to_component.find(vertex) == vertex_to_component.end()) {
            dfs(vertex, current_component++);
        }
    }

    // Organize vertices by component
    std::vector<std::vector<size_t>> components(current_component);
    for (const auto& [vertex, component] : vertex_to_component) {
        components[component].push_back(vertex);
    }

    return components;
}

/**
 * @brief Calculates the absolute and relative deviations between direct edges
 *        and their corresponding two-step paths through intermediate nodes.
 *
 * @details For each edge (i,k), this function finds common neighbors j
 *          and computes the minimum deviation between w(i,k) and w(i,j)+w(j,k).
 *          These deviations can be used to identify geometrically redundant edges
 *          using either percentile or statistical outlier detection methods.
 *
 * @param adj_list   Adjacency list where adj_list[i] contains neighbors of node i
 * @param weight_list Weight list where weight_list[i][j] is the weight of edge
 *                   from i to adj_list[i][j]
 *
 * @return edge_weight_deviations_t Structure containing vectors of absolute and
 *                                 relative deviations for all valid edges
 */
edge_weight_deviations_t set_wgraph_t::compute_edge_weight_deviations() const {

    size_t n_vertices = adjacency_list.size();
    edge_weight_deviations_t result;

    // For each node i
    for (size_t i = 0; i < n_vertices; ++i) {
        const auto& i_neighbors = adjacency_list[i];

        // Process each edge (i,k)
        for (const auto& edge_ik : i_neighbors) {
            size_t k = edge_ik.vertex;
            double w_ik = edge_ik.weight;

            // Only process each edge once
            if (i >= k) continue;

            double min_abs_deviation = std::numeric_limits<double>::infinity();
            double min_rel_deviation = std::numeric_limits<double>::infinity();

            const auto& k_neighbors = adjacency_list[k];

            // Check all potential intermediate nodes (common neighbors)
            for (const auto& edge_ij : i_neighbors) {
                size_t j = edge_ij.vertex;
                double w_ij = edge_ij.weight;

                if (j == k) continue;  // Skip direct edge

                // Create a search key for the set lookup
                edge_info_t search_key{j, 0.0};  // Weight doesn't matter for search
                auto it = k_neighbors.find(search_key);

                if (it != k_neighbors.end()) {  // j is a common neighbor
                    double w_jk = it->weight;

                    // Compute absolute deviation
                    double abs_dev = std::abs(w_ij + w_jk - w_ik);

                    // Compute relative deviation (with safeguard against division by zero)
                    double sum = w_ij + w_jk;
                    double rel_dev = (sum != 0.0) ? abs_dev / sum : std::numeric_limits<double>::infinity(); // the fraction of the indirect path that's "wasted" when compared to the direct path

                    // Update minimum deviations
                    if (abs_dev < min_abs_deviation) {
                        min_abs_deviation = abs_dev;
                        min_rel_deviation = rel_dev;
                    }
                }
            }

            // Store results if we found at least one alternative path
            if (min_abs_deviation != std::numeric_limits<double>::infinity()) {
                result.absolute_deviations.push_back(min_abs_deviation);
                result.relative_deviations.push_back(min_rel_deviation);
            }
        }
    }

    return result;
}

/**
* Compute Relative Deviations for Graph Edge Weights
*
* @description
* Calculates the relative deviation for each edge in a weighted graph to identify
* potentially redundant edges. For each edge (i,k), the function finds the best
* alternative path (i,j,k) through a common neighbor j and calculates how well
* the triangle inequality holds.
*
* @details
* The relative deviation is defined as:
* \deqn{\frac{|w(i,j) + w(j,k) - w(i,k)|}{w(i,j) + w(j,k)}}
*
* This measures the fraction of the indirect path length that differs from
* the direct path. Values close to 0 indicate potential geometric redundancy,
* meaning the edge (i,k) follows almost perfectly from the edges (i,j) and (j,k).
* The relative deviation always falls within the range [0,1] (or infinity if no path exists).
*
* For each edge, the function examines all possible intermediate nodes (common
* neighbors of both endpoints) and finds the one that minimizes the relative
* deviation.
*
* Each edge is processed only once (when source < target) to avoid duplicate work.
*
* @return A vector of \code{edge_weight_rel_deviation_t} structures containing:
*   \describe{
*     \item{rel_deviation}{The minimum relative deviation found for the edge}
*     \item{source}{The source node i of the edge}
*     \item{target}{The target node k of the edge}
*     \item{best_intermediate}{The intermediate node j that gives the minimum deviation}
*   }
*
* @examples
* # Create a weighted graph
* graph <- set_wgraph_t$new(adjacency_matrix)
*
* # Compute relative deviations
* rel_deviations <- graph$compute_edge_weight_rel_deviations()
*
* # Identify potentially redundant edges using a threshold
* epsilon <- 0.1
* redundant_edges <- rel_deviations[rel_deviations$rel_deviation < epsilon, ]
*
* # Examine the distribution of relative deviations
* hist(rel_deviations$rel_deviation, breaks = 30,
*      main = "Edge Weight Relative Deviations")
*
* @seealso
* The relative deviation concept relates to metric space embedding and spanner graphs.
* This analysis can be useful for graph simplification and identifying structural redundancies.
*
* @export
*/
std::vector<edge_weight_rel_deviation_t> set_wgraph_t::compute_edge_weight_rel_deviations() const {

    size_t n_vertices = adjacency_list.size();
    std::vector<edge_weight_rel_deviation_t> result;

    // For each node i
    for (size_t i = 0; i < n_vertices; ++i) {
        const auto& i_neighbors = adjacency_list[i];
        // Process each edge (i,k)
        for (const auto& edge_ik : i_neighbors) {
            size_t k = edge_ik.vertex;
            double w_ik = edge_ik.weight;
            // Only process each edge once
            if (i >= k) continue;
            double min_rel_deviation = std::numeric_limits<double>::infinity();
            size_t best_intermediate = 0; // Initialize best intermediate
            const auto& k_neighbors = adjacency_list[k];
            // Check all potential intermediate nodes (common neighbors)
            for (const auto& edge_ij : i_neighbors) {
                size_t j = edge_ij.vertex;
                double w_ij = edge_ij.weight;
                if (j == k) continue;  // Skip direct edge
                // Create a search key for the set lookup
                edge_info_t search_key{j, 0.0};  // Weight doesn't matter for search
                auto it = k_neighbors.find(search_key);
                if (it != k_neighbors.end()) {  // j is a common neighbor
                    double w_jk = it->weight;

                    double rel_dev = (w_ij + w_jk) / w_ik - 1.0;

                    #if 0
                    // Compute absolute deviation
                    double abs_dev = std::abs(w_ij + w_jk - w_ik);
                    // Compute relative deviation
                    double sum = w_ij + w_jk;
                    double rel_dev = (sum != 0.0) ? abs_dev / sum : std::numeric_limits<double>::infinity();
                    #endif

                    // Update minimum deviations
                    if (rel_dev < min_rel_deviation) {
                        min_rel_deviation = rel_dev;
                        best_intermediate = j;
                    }
                }
            }
            // Store results if we found at least one alternative path
            if (min_rel_deviation != std::numeric_limits<double>::infinity()) {
                result.emplace_back(min_rel_deviation, i, k, best_intermediate);
            }
        }
    }
    return result;
}


#if 0
edge_weight_deviations_t set_wgraph_t::compute_edge_weight_deviations() const {
    size_t n_vertices = adjacency_list.size();
    edge_weight_deviations_t result;

    // Define a threshold for "large" deviations
    const double LARGE_DEVIATION_THRESHOLD = 0.01;

    // For each node i
    for (size_t i = 0; i < n_vertices; ++i) {
        const auto& i_neighbors = adjacency_list[i];

        // Process each edge (i,k)
        for (const auto& edge_ik : i_neighbors) {
            size_t k = edge_ik.vertex;
            double w_ik = edge_ik.weight;

            if (i >= k) continue;  // Skip if i >= k to process each edge once

            double min_abs_deviation = std::numeric_limits<double>::infinity();
            double min_rel_deviation = std::numeric_limits<double>::infinity();
            size_t best_intermediate = n_vertices; // Invalid vertex as default

            const auto& k_neighbors = adjacency_list[k];

            // Check all potential intermediate nodes (common neighbors)
            for (const auto& edge_ij : i_neighbors) {
                size_t j = edge_ij.vertex;
                double w_ij = edge_ij.weight;

                if (j == k) continue;  // Skip direct edge

                // Create a temporary edge_info_t to use as key for lookup
                edge_info_t search_key{j, 0.0};  // weight doesn't matter for search
                auto it = k_neighbors.find(search_key);

                if (it != k_neighbors.end()) {  // j is a common neighbor
                    double w_jk = it->weight;

                    // Compute absolute deviation
                    double abs_dev = std::abs(w_ij + w_jk - w_ik);

                    // Compute relative deviation
                    double rel_dev = (w_ik > 1e-10) ? abs_dev / w_ik : std::numeric_limits<double>::infinity();

                    // Update minimum deviations
                    if (abs_dev < min_abs_deviation) {
                        min_abs_deviation = abs_dev;
                        min_rel_deviation = rel_dev;
                        best_intermediate = j;
                    }
                }
            }

            // Store results if we found at least one alternative path
            if (min_abs_deviation != std::numeric_limits<double>::infinity()) {
                result.absolute_deviations.push_back(min_abs_deviation);
                result.relative_deviations.push_back(min_rel_deviation);

                // Print triangles with large deviations for debugging
                if (min_abs_deviation > LARGE_DEVIATION_THRESHOLD) {
                    Rprintf("Large deviation found: Triangle (%zu, %zu, %zu)\n",
                           i+1, best_intermediate+1, k+1); // +1 for R indexing
                    Rprintf("  Weights: w(%zu,%zu) = %.10f, w(%zu,%zu) = %.10f, w(%zu,%zu) = %.10f\n",
                           i+1, best_intermediate+1,
                           adjacency_list[i].find(edge_info_t{best_intermediate, 0.0})->weight,
                           best_intermediate+1, k+1,
                           adjacency_list[best_intermediate].find(edge_info_t{k, 0.0})->weight,
                           i+1, k+1, w_ik);
                    Rprintf("  Absolute deviation: %.10f\n", min_abs_deviation);
                    Rprintf("  Relative deviation: %.10f\n\n", min_rel_deviation);
                }
            }
        }
    }

    return result;
}
#endif


/**
 * @brief Computes geometric edge weight deviations for graph structure analysis.
 *
 * @details This function serves as an interface between R and C++ to compute
 *          absolute and relative deviations between direct edges and their
 *          corresponding two-step paths through intermediate nodes. For each
 *          edge (i,k), it finds the minimum deviation between w(i,k) and
 *          w(i,j)+w(j,k) for all common neighbors j.
 *
 *          These deviations can be used to identify geometrically redundant edges
 *          for graph pruning, using either percentile or statistical outlier
 *          detection methods. This can simplify the graph structure while
 *          preserving essential geometric properties needed for gradient trajectory
 *          analysis.
 *
 * @param s_adj_list An R list representing the adjacency list structure where
 *                   s_adj_list[[i]] contains indices of vertices adjacent to vertex i
 *                   (VECSXP of INTSXP vectors)
 * @param s_weight_list An R list representing the weight list where s_weight_list[[i]][j]
 *                      is the weight of edge from vertex i to s_adj_list[[i]][j]
 *                      (VECSXP of REALSXP vectors)
 *
 * @return An R list (VECSXP) with two named elements:
 *         - "absolute_deviations": A numeric vector (REALSXP) of absolute deviations
 *           |w(i,j) + w(j,k) - w(i,k)|
 *         - "relative_deviations": A numeric vector (REALSXP) of relative deviations
 *           |w(i,j) + w(j,k) - w(i,k)| / w(i,k)
 *
 * @note This function properly handles R's garbage collection by using PROTECT/UNPROTECT
 *       mechanism for all allocated R objects.
 */
SEXP S_compute_edge_weight_deviations(SEXP s_adj_list, SEXP s_weight_list) {
    // Convert R objects to C++ vectors
    std::vector<std::vector<int>> adj_vect       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

    // Create the graph once (fixed variable names)
    set_wgraph_t graph(adj_vect, weight_vect);

    // Compute the deviations using the optimized method
    edge_weight_deviations_t deviations = graph.compute_edge_weight_deviations();

    // Create the return list using R's C API
    const char* names[] = {
        "absolute_deviations",
        "relative_deviations",
        NULL
    };
    int n_elements = 2;

    // Track protection count explicitly
    int protect_count = 0;
    SEXP result = PROTECT(allocVector(VECSXP, n_elements));
    protect_count++;

    SEXP result_names = PROTECT(allocVector(STRSXP, n_elements));
    protect_count++;

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, mkChar(names[i]));
    }
    setAttrib(result, R_NamesSymbol, result_names);

    // Create and set the absolute deviations vector
    SEXP abs_dev = PROTECT(allocVector(REALSXP, deviations.absolute_deviations.size()));
    protect_count++;
    double* abs_ptr = REAL(abs_dev);
    std::copy(deviations.absolute_deviations.begin(), deviations.absolute_deviations.end(), abs_ptr);
    SET_VECTOR_ELT(result, 0, abs_dev);

    // Create and set the relative deviations vector
    SEXP rel_dev = PROTECT(allocVector(REALSXP, deviations.relative_deviations.size()));
    protect_count++;
    double* rel_ptr = REAL(rel_dev);
    std::copy(deviations.relative_deviations.begin(), deviations.relative_deviations.end(), rel_ptr);
    SET_VECTOR_ELT(result, 1, rel_dev);

    UNPROTECT(protect_count);
    return result;
}


/**
 * @brief Compute relative deviations for graph edge weights
 *
 * @details
 * This function calculates the relative deviation for each edge in a weighted graph
 * to identify potentially redundant edges. For each edge (i,k), it finds the best
 * alternative path (i,j,k) through a common neighbor j and calculates how well
 * the triangle inequality holds using the formula:
 *
 * \f[ \frac{|w(i,j) + w(j,k) - w(i,k)|}{w(i,j) + w(j,k)} \f]
 *
 * This measure ranges from 0 to 1 (or infinity if no alternative path exists), where:
 * - Values close to 0 indicate edges that are nearly perfectly redundant (geometric redundancy)
 * - Values close to 1 indicate edges that cannot be well-approximated by alternative paths
 *
 * The function serves as an interface between R and the C++ implementation of
 * compute_edge_weight_rel_deviations().
 *
 * @param graph_sexp An R external pointer to a C++ set_wgraph_t object
 *
 * @return An R matrix with four columns:
 *   - source: The source node ID of the edge
 *   - target: The target node ID of the edge
 *   - best_intermediate: The intermediate node that gives the minimum deviation
 *   - rel_deviation: The minimum relative deviation found for the edge
 *
 * @see compute_edge_weight_rel_deviations
 */
SEXP S_compute_edge_weight_rel_deviations(SEXP s_adj_list, SEXP s_weight_list) {

    // Convert R objects to C++ vectors
    std::vector<std::vector<int>> adj_vect       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

    // Create the graph once (fixed variable names)
    set_wgraph_t graph(adj_vect, weight_vect);

    // Compute the deviations using the optimized method
    auto results = graph.compute_edge_weight_rel_deviations();

    // Create the result matrix
    size_t n_edges = results.size();
    size_t protect_count = 0;
    SEXP result_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_edges, 4));
    protect_count++;
    double* result_data = REAL(result_matrix);

    // Fill the matrix with data
    size_t n_edges_times_two   = 2 * n_edges;
    size_t n_edges_times_three = 3 * n_edges;
    for (size_t i = 0; i < n_edges; i++) {
        // In column-major order (as R uses)
        result_data[i]                       = static_cast<double>(results[i].source + 1); // Column 1: source
        result_data[i + n_edges]             = static_cast<double>(results[i].target + 1); // Column 2: target
        result_data[i + n_edges_times_two]   = static_cast<double>(results[i].best_intermediate + 1); // Column 3: best_intermediate
        result_data[i + n_edges_times_three] = results[i].rel_deviation;                // Column 4: rel_deviation
    }

    // Set column names
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 4));
    protect_count++;
    SET_STRING_ELT(colnames, 0, Rf_mkChar("source"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("target"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("best_intermediate"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("rel_deviation"));

    // Set the colnames attribute
    Rf_setAttrib(result_matrix, R_DimNamesSymbol,
                 Rf_list2(R_NilValue, colnames));

    UNPROTECT(protect_count);
    return result_matrix;
}

/**
 * @brief Removes an edge between two vertices in the graph
 *
 * @param graph The graph to modify
 * @param v1 First vertex of the edge
 * @param v2 Second vertex of the edge
 *
 * Removes the edge in both directions (undirected graph).
 * Prints debug information about the removed edge.
 */
void set_wgraph_t::remove_edge(size_t v1, size_t v2) {

    // Create a dummy edge_info_t with any weight (0.0) since our comparator only checks vertices
    edge_info_t dummy{v2, 0.0};
    // Use set's efficient find operation
    auto it1 = adjacency_list[v1].find(dummy);
    if (it1 != adjacency_list[v1].end()) {
        adjacency_list[v1].erase(it1);
        // Do the same for the reverse edge
        edge_info_t reverse_dummy{v1, 0.0};
        adjacency_list[v2].erase(reverse_dummy);
        //Rprintf("DEBUG: Removed edge (%d, %d) with weight %.6f\n", v1, v2, weight);
    }
}

#if 0
SEXP S_remove_redundant_edges(SEXP s_adj_list, SEXP s_weight_list) {

    // Convert R objects to C++ vectors
    std::vector<std::vector<int>> adj_vect       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

    // Create the graph once (fixed variable names)
    set_wgraph_t graph(adj_vect, weight_vect);

    // Compute the deviations using the optimized method
    auto rel_deviations = graph.compute_edge_weight_rel_deviations();

    const double EPSILON = 1e-16; // when a value is considered 0 in R ???

    // Create the result matrix
    size_t n_edges = rel_deviations.size();
    for (size_t i = 0; i < n_edges; i++) {
        if (rel_deviations[i].rel_deviation < EPSILON) {
            graph.remove_edge(rel_deviations[i].source, rel_deviations[i].target);
        }
    }

    // Creating return list
    int n_protected = 0;
    const int N_COMPONENTS = 2;
    SEXP r_result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;
    int n_vertices = graph.adjacency_list.size();
    SEXP r_adj_list = PROTECT(allocVector(VECSXP, n_vertices)); n_protected++;
    SEXP r_weight_list = PROTECT(allocVector(VECSXP, n_vertices)); n_protected++;

    // Convert the set-based representation back to R lists
    for (int i = 0; i < n_vertices; ++i) {
        const auto& neighbors = graph.adjacency_list[i];

        // Create vectors for this vertex's adjacency list and weights
        SEXP r_adj = PROTECT(allocVector(INTSXP, neighbors.size()));
        SEXP r_weights = PROTECT(allocVector(REALSXP, neighbors.size()));

        // Fill the vectors
        int idx = 0;
        for (const auto& [neighbor, weight] : neighbors) {
            // Convert to 1-based indices for R
            INTEGER(r_adj)[idx] = neighbor + 1;
            REAL(r_weights)[idx] = weight;
            ++idx;
        }

        SET_VECTOR_ELT(r_adj_list, i, r_adj);
        SET_VECTOR_ELT(r_weight_list, i, r_weights);
        UNPROTECT(2); // for r_adj and r_weights
    }

    // Set components in the result list
    SET_VECTOR_ELT(r_result, 0, r_adj_list);
    SET_VECTOR_ELT(r_result, 1, r_weight_list);

    // Set names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("weight_list"));
    setAttrib(r_result, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return r_result;
}
#endif

/**
 * @brief Remove redundant edges from a weighted graph while preserving connectivity
 *
 * @details
 * This function identifies and removes edges with relative deviation values extremely
 * close to zero, indicating that they are geometrically redundant. An edge (i,k) is
 * considered redundant if there exists an alternative path (i,j,k) such that
 * w(i,k) ≈ w(i,j) + w(j,k).
 *
 * The function ensures that no vertex becomes disconnected (isolated) from the graph
 * as a result of edge removal. If removing an edge would leave any vertex with no
 * neighbors, the function halts with an error message that identifies:
 * - The source and target vertices of the edge being removed
 * - The intermediate vertex that provides the alternative path
 * - Which vertex would become isolated
 *
 * Edges are processed in their original order, and each redundant edge is removed
 * if it passes the connectivity check. The algorithm only addresses the extreme case
 * of vertex isolation (where a vertex would have zero neighbors) as specified in
 * the requirements.
 *
 * @param s_adj_list R list of integer vectors representing adjacency lists for each vertex
 * @param s_weight_list R list of numeric vectors representing edge weights corresponding to adjacencies
 *
 * @return An R list with two components:
 *         - adj_list: The updated adjacency list after redundant edge removal
 *         - weight_list: The updated weight list corresponding to the new graph structure
 *
 * @throws Rf_error If removing an edge would disconnect a vertex from the graph
 * @throws Rf_error If any other error occurs during processing
 *
 * @note The function uses a threshold of 1e-16 to determine when a relative deviation
 *       is considered zero.
 *
 * @see compute_edge_weight_rel_deviations
 * @see set_wgraph_t::remove_edge
 */
SEXP S_remove_redundant_edges(SEXP s_adj_list, SEXP s_weight_list) {

    size_t n_protected = 0;

    try {
        // Convert R objects to C++ vectors
        std::vector<std::vector<int>> adj_vect       = convert_adj_list_from_R(s_adj_list);
        std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

        // Create the graph once
        set_wgraph_t graph(adj_vect, weight_vect);

        // Compute the deviations using the optimized method
        auto rel_deviations = graph.compute_edge_weight_rel_deviations();
        const double EPSILON = 1e-16; // Threshold for considering a value as zero

        // Track removed edges
        //int removed_count = 0;

        // Process edges
        for (size_t i = 0; i < rel_deviations.size(); i++) {
            if (rel_deviations[i].rel_deviation < EPSILON) {
                size_t source = rel_deviations[i].source;
                size_t target = rel_deviations[i].target;
                size_t intermediate = rel_deviations[i].best_intermediate;

                // Check if removing this edge would isolate any vertex
                if (graph.adjacency_list[source].size() <= 1) {
                    Rf_error("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                             source + 1, target + 1, intermediate + 1, source + 1);
                }

                if (graph.adjacency_list[target].size() <= 1) {
                    Rf_error("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                             source + 1, target + 1, intermediate + 1, target + 1);
                }

                // Safe to remove the edge
                graph.remove_edge(source, target);
                //removed_count++;
            }
        }

        const int N_COMPONENTS = 2;
        SEXP r_result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;
        int n_vertices = graph.adjacency_list.size();
        SEXP r_adj_list = PROTECT(allocVector(VECSXP, n_vertices)); n_protected++;
        SEXP r_weight_list = PROTECT(allocVector(VECSXP, n_vertices)); n_protected++;

        // Convert the set-based representation back to R lists
        for (int i = 0; i < n_vertices; ++i) {
            const auto& neighbors = graph.adjacency_list[i];

            // Create vectors for this vertex's adjacency list and weights
            SEXP r_adj = PROTECT(allocVector(INTSXP, neighbors.size()));
            SEXP r_weights = PROTECT(allocVector(REALSXP, neighbors.size()));

            // Fill the vectors
            int idx = 0;
            for (const auto& [neighbor, weight] : neighbors) {
                // Convert to 1-based indices for R
                INTEGER(r_adj)[idx] = neighbor + 1;
                REAL(r_weights)[idx] = weight;
                ++idx;
            }

            SET_VECTOR_ELT(r_adj_list, i, r_adj);
            SET_VECTOR_ELT(r_weight_list, i, r_weights);
            UNPROTECT(2); // for r_adj and r_weights
        }

        // Set components in the result list
        SET_VECTOR_ELT(r_result, 0, r_adj_list);
        SET_VECTOR_ELT(r_result, 1, r_weight_list);

        // Set names for return list
        SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
        SET_STRING_ELT(names, 0, mkChar("adj_list"));
        SET_STRING_ELT(names, 1, mkChar("weight_list"));
        setAttrib(r_result, R_NamesSymbol, names);

        // Rprintf("Removed %zu redundant edges\n", removed_count);

        UNPROTECT(n_protected);
        return r_result;

    } catch (const std::exception& e) {
        if (n_protected > 0) UNPROTECT(n_protected);
        error("Error removing redundant edges: %s", e.what());
    } catch (...) {
        if (n_protected > 0) UNPROTECT(n_protected);
        error("Unknown error occurred while removing redundant edges");
    }
}

/**
 * @brief Precomputes all edge weights in the graph for efficient lookup
 *
 * Populates the internal edge_weights map for O(1) weight lookup by vertex pair.
 * For undirected graphs, only one entry is stored per edge with the
 * convention that vertex indices are ordered (i < j).
 *
 * @note This method is const as it only modifies mutable cache members
 */

/**
 * @brief Precomputes all edge weights in the graph for efficient lookup
 *
 * This function creates a hash map that stores weights for all edges in the graph,
 * enabling O(1) weight lookup by vertex pair. For undirected graphs, only one entry
 * is stored per edge with the convention that the vertex indices are ordered (i < j).
 *
 * The resulting map uses a custom hash function for std::pair<size_t, size_t> to
 * ensure efficient storage and retrieval. This precomputation is particularly useful
 * for algorithms that require frequent edge weight lookups.
 *
 * @return edge_weights_t A hash map where keys are vertex pairs (i,j) with i < j and
 *         values are the corresponding edge weights
 *
 * @note The function automatically reserves appropriate space in the hash map to
 *       avoid rehashing operations during population
 * @note Time complexity: O(E) where E is the number of edges in the graph
 * @note Space complexity: O(E) for storing the edge weights
 *
 * @see edge_weights_t For the hash map type definition
 * @see size_t_pair_hash_t For the hash function used by the map
 */
void set_wgraph_t::precompute_edge_weights() const {
        edge_weights.clear();

        // Reserve space to avoid rehashing
        size_t estimated_edges = 0;
        for (const auto& adj : adjacency_list) {
            estimated_edges += adj.size();
        }
        estimated_edges /= 2;  // Divide by 2 since we store unique edges (i < j)
        edge_weights.reserve(estimated_edges);

        // Go through adjacency list to compute weights
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                size_t j = edge.vertex;
                if (i < j) {
                    edge_weights[{i, j}] = edge.weight;
                }
            }
        }

        edge_weights_computed = true;
}

/**
 * @brief Computes multiple percentiles of edge weights in the graph
 *
 * This function calculates percentiles of the edge weight distribution by directly
 * accessing the adjacency list. For each requested probability value (0.0 to 1.0),
 * it computes the corresponding percentile of edge weights.
 *
 * @param probs Vector of probability values (between 0.0 and 1.0) for which to compute percentiles
 * @return std::vector<double> Vector of computed percentiles, in the same order as the input probabilities
 *
 * @throws std::invalid_argument If any probability value is outside the [0,1] range
 * @throws std::runtime_error If the graph has no edges
 *
 * @note Time complexity: O(E log E) where E is the number of edges (dominated by sorting)
 * @note Space complexity: O(E) for storing all edge weights
 * @note For undirected graphs, each edge is counted only once
 *
 * @example
 * // Get 25th, 50th, 75th and 90th percentiles
 * std::vector<double> percentiles = graph.compute_weight_percentiles({0.25, 0.5, 0.75, 0.9});
 */
std::vector<double> set_wgraph_t::compute_weight_percentiles(const std::vector<double>& probs) const {
    // Validate input probabilities
    for (const double p : probs) {
        if (p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Probability values must be between 0.0 and 1.0");
        }
    }

    // Collect all edge weights
    std::vector<double> weights;
    weights.reserve(num_vertices() * 2); // Rough estimate to avoid frequent reallocations

    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        for (const auto& edge : adjacency_list[i]) {
            // For undirected graphs, only collect each edge once
            if (i < edge.vertex) {
                weights.push_back(edge.weight);
            }
        }
    }

    // Check if we have any edges
    if (weights.empty()) {
        throw std::runtime_error("Cannot compute percentiles: graph has no edges");
    }

    // Sort the weights
    std::sort(weights.begin(), weights.end());

    // Compute requested percentiles
    std::vector<double> result;
    result.reserve(probs.size());

    const size_t n = weights.size();
    for (double p : probs) {
        // Calculate the position
        double pos = p * (n - 1);

        // For exact positions, just take the value
        if (std::floor(pos) == pos) {
            result.push_back(weights[static_cast<size_t>(pos)]);
        }
        // For positions between two values, interpolate
        else {
            size_t lower_idx = static_cast<size_t>(std::floor(pos));
            size_t upper_idx = static_cast<size_t>(std::ceil(pos));
            double fraction = pos - std::floor(pos);

            double lower_val = weights[lower_idx];
            double upper_val = weights[upper_idx];

            // Linear interpolation
            result.push_back(lower_val + fraction * (upper_val - lower_val));
        }
    }

    return result;
}

/**
* @brief Extracts edge lengths from a path in the graph
*
* This function converts a path (sequence of vertices) into a vector of edge lengths
* by looking up the weight of each edge in the graph's adjacency list.
*
* @param path Path structure containing ordered sequence of vertices
* @return std::vector<double> Vector of edge lengths corresponding to each segment of the path
*
* @note Edge lengths are determined by the weight attribute in the adjacency list
* @note If an edge between consecutive vertices doesn't exist in the adjacency list,
*       the length will be set to 0.0
*
* @pre The path object must contain at least one vertex
* @pre Consecutive vertices in the path must be directly connected in the graph
*/
std::vector<double> set_wgraph_t::extract_edge_lengths(const path_t& path) const {
    std::vector<double> edge_lengths;
    edge_lengths.reserve(path.vertices.size() - 1);

    for (size_t i = 0; i < path.vertices.size() - 1; ++i) {
        size_t v1 = path.vertices[i];
        size_t v2 = path.vertices[i + 1];

        // Find the edge weight between v1 and v2
        double edge_length = 0.0;

        // Check which adjacency list is smaller and iterate over that one
        if (adjacency_list[v1].size() <= adjacency_list[v2].size()) {
            // Iterate over v1's adjacency list if it's smaller
            for (const auto& edge : adjacency_list[v1]) {
                if (edge.vertex == v2) {
                    edge_length = edge.weight;
                    break;
                }
            }
        } else {
            // Iterate over v2's adjacency list if it's smaller
            for (const auto& edge : adjacency_list[v2]) {
                if (edge.vertex == v1) {
                    edge_length = edge.weight;
                    break;
                }
            }
        }

        edge_lengths.push_back(edge_length);
    }

    return edge_lengths;
}
