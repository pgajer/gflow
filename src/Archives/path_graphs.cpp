#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <limits>
#include <map>
#include <utility>

#include "path_graphs.h"
#include "SEXP_cpp_conversion_utils.h"

extern "C" {
    SEXP S_create_path_graph_series(SEXP s_adj_list,
                                    SEXP s_weight_list,
                                    SEXP s_h_values);

    SEXP S_create_path_graph_plus(SEXP s_adj_list,
                                  SEXP s_edge_length_list,
                                  SEXP s_h);

    SEXP S_create_path_graph_plm(SEXP s_adj_list,
                                 SEXP s_edge_length_list,
                                 SEXP s_h);
}

/**
 * @brief Creates a path graph P_{\bullet \leq h}(G) with hop counts from a given weighted undirected graph G
 *
 * @details The path graph augments G^{\leq h} by storing, for each edge (v,w),
 *          the shortest path connecting v to w in the original graph G.
 *          This implementation also stores hop counts for efficient subgraph generation.
 *
 *          Important Notes:
 *          - This implementation assumes an undirected graph as input
 *          - For memory efficiency, shortest paths are stored only for vertex pairs (v,w) where v < w
 *          - To retrieve the path between vertices a and b, always query shortest_paths[{min(a,b), max(a,b)}]
 *          - Each stored path is a vector of vertices representing the sequence from source to target
 *
 * @param adj_list The adjacency list of the input undirected graph G
 * @param weight_list The weight list of the input undirected graph G. weight_list[v][i] represents
 *                    the weight of the edge between vertex v and its i-th neighbor in adj_list[v]
 * @param h The maximum number of hops allowed
 * @return path_graph_plus_t The h-hop path graph with additional hop count information and shortest paths
 *
 * @note Time Complexity: O(V * (V + E)) where V is the number of vertices and E is the number of edges
 * @note Space Complexity: O(V^2) in worst case for storing all paths
 */
path_graph_plus_t create_path_graph_plus(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    const int n_vertices = adj_list.size();
    path_graph_plus_t result;
    result.adj_list.resize(n_vertices);
    result.weight_list.resize(n_vertices);
    result.hop_list.resize(n_vertices);

    // Pre-allocate vectors to avoid reallocations
    std::vector<double> distances(n_vertices);
    std::vector<int> hops(n_vertices);
    std::vector<int> parent(n_vertices);
    std::vector<bool> in_queue(n_vertices);

    // Temporary storage for neighbors with their properties
    struct neighbor_info_t {
        int vertex;
        double distance;
        int hops;
        neighbor_info_t(int v, double d, int h) : vertex(v), distance(d), hops(h) {}
    };
    std::vector<neighbor_info_t> current_neighbors;
    current_neighbors.reserve(n_vertices / 2);  // Conservative estimate

    for (int start = 0; start < n_vertices; ++start) {
        // Reset arrays (faster than reinitializing)
        std::fill(distances.begin(), distances.end(), std::numeric_limits<double>::infinity());
        std::fill(hops.begin(), hops.end(), std::numeric_limits<int>::max());
        std::fill(parent.begin(), parent.end(), -1);
        std::fill(in_queue.begin(), in_queue.end(), false);
        current_neighbors.clear();

        // Initialize start vertex
        distances[start] = 0;
        hops[start] = 0;

        // Use deque for better cache locality
        std::deque<int> q{start};
        in_queue[start] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop_front();
            in_queue[current] = false;

            // Early termination if we've reached hop limit
            if (hops[current] >= h) continue;

            // Cache current values
            const double current_dist = distances[current];
            const int current_hops = hops[current];
            const int next_hops = current_hops + 1;

            // Process neighbors
            const auto& current_adj = adj_list[current];
            const auto& current_weights = weight_list[current];
            const size_t n_neighbors = current_adj.size();

            for (size_t i = 0; i < n_neighbors; ++i) {
                const int neighbor = current_adj[i];

                if (neighbor == current) continue;  // Skip self-loops

                const double new_distance = current_dist + current_weights[i];

                // Update if we found a shorter path or same distance with fewer hops
                if (new_distance < distances[neighbor]) {
                    distances[neighbor] = new_distance;
                    hops[neighbor] = next_hops;
                    parent[neighbor] = current;

                    // Add to queue if not already in it
                    if (!in_queue[neighbor] && next_hops <= h) {
                        q.push_back(neighbor);
                        in_queue[neighbor] = true;
                    }
                }
            }
        }

        // Collect neighbors and their properties
        for (int v = 0; v < n_vertices; ++v) {
            if (v != start && hops[v] <= h) {
                current_neighbors.emplace_back(v, distances[v], hops[v]);

                // Store path for pairs where start < v
                if (start < v) {
                    auto& path = result.shortest_paths[{start, v}];
                    path.clear();
                    path.reserve(hops[v] + 1);  // Exact size known

                    // Reconstruct path
                    for (int current = v; current != -1; current = parent[current]) {
                        path.push_back(current);
                    }
                    std::reverse(path.begin(), path.end());
                }
            }
        }

        // Bulk insert neighbors and their properties
        const size_t n_new_neighbors = current_neighbors.size();
        result.adj_list[start].resize(n_new_neighbors);
        result.weight_list[start].resize(n_new_neighbors);
        result.hop_list[start].resize(n_new_neighbors);

        for (size_t i = 0; i < n_new_neighbors; ++i) {
            const auto& info = current_neighbors[i];
            result.adj_list[start][i] = info.vertex;
            result.weight_list[start][i] = info.distance;
            result.hop_list[start][i] = info.hops;
        }
    }

    return result;
}

/**
 * @brief Creates a k-hop neighborhood (hHN) graph from the input graph.
 *
 * This function generates a new graph where each vertex is connected to all other vertices
 * that are within k hops in the original graph. The edge weights in the new graph represent
 * the shortest path distances in the original graph.
 *
 * @param adj_list A vector of vectors representing the adjacency list of the input graph.
 *                 adj_list[i] contains the indices of vertices adjacent to vertex i.
 * @param weight_list A vector of vectors representing the edge weights of the input graph.
 *                    weight_list[i][j] is the weight of the edge from vertex i to its j-th neighbor.
 * @param h The maximum number of hops to consider for neighborhood connectivity.
 *
 * @return A pair containing:
 *         - First: A vector of vectors representing the adjacency list of the hHN graph.
 *         - Second: A vector of vectors representing the edge weights of the hHN graph.
 *
 * @note The function assumes that the input adjacency list and weight list are valid and consistent.
 * @note The time complexity is O(n * (m + n log n)), where n is the number of vertices and m is the number of edges.
 * @note The space complexity is O(n^2) in the worst case, where the hHN graph becomes a complete graph.
 *
 * @warning If k is greater than or equal to the diameter of the graph, the resulting hHN graph may be complete.
 *
 * @see For background on k-hop neighborhood graphs, refer to:
 *      - Asano, T., & Uno, T. (2018). Efficient construction of k-hop neighborhood graphs.
 *        In 2018 Proceedings of the Twentieth Workshop on Algorithm Engineering and Experiments (ALENEX) (pp. 156-167).
 *
 * @example
 * std::vector<std::vector<int>> adj_list = {{1, 2}, {0, 2}, {0, 1, 3}, {2}};
 * std::vector<std::vector<double>> weight_list = {{1.0, 2.0}, {1.0, 1.5}, {2.0, 1.5, 3.0}, {3.0}};
 * int h = 2;
 * auto [hhn_adj_list, hhn_weight_list] = create_hHN_graph(adj_list, weight_list, h);
 */
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> create_hHN_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    int n_vertices = adj_list.size();

    std::vector<std::vector<int>> hhn_adj_list(n_vertices);
    std::vector<std::vector<double>> hhn_weight_list(n_vertices);

    // Pre-allocate vectors to avoid reallocations
    std::vector<double> distances(n_vertices);
    std::vector<int> hops(n_vertices);
    std::vector<bool> in_queue(n_vertices);

    for (int start = 0; start < n_vertices; ++start) {
        // Reset arrays
        std::fill(distances.begin(), distances.end(), std::numeric_limits<double>::infinity());
        std::fill(hops.begin(), hops.end(), std::numeric_limits<int>::max());
        std::fill(in_queue.begin(), in_queue.end(), false);

        // Initialize start vertex
        distances[start] = 0;
        hops[start] = 0;

        std::queue<int> q;
        q.push(start);
        in_queue[start] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop();
            in_queue[current] = false;

            // Early termination if we've reached hop limit
            if (hops[current] >= h) continue;

            // Process neighbors
            const auto& current_adj = adj_list[current];
            const auto& current_weights = weight_list[current];

            for (size_t i = 0; i < current_adj.size(); ++i) {
                const int neighbor = current_adj[i];

                if (neighbor == current) continue;  // Skip self-loops

                const double new_distance = distances[current] + current_weights[i];
                const int next_hops = hops[current] + 1;

                if (new_distance < distances[neighbor]) {
                    distances[neighbor] = new_distance;
                    hops[neighbor] = next_hops;

                    // Only add to queue if not already in queue and within hop limit
                    if (!in_queue[neighbor] && next_hops <= h) {
                        q.push(neighbor);
                        in_queue[neighbor] = true;
                    }
                }
            }
        }

        // Collect neighbors within h hops
        for (int v = 0; v < n_vertices; ++v) {
            if (v != start && hops[v] <= h) {
                hhn_adj_list[start].push_back(v);
                hhn_weight_list[start].push_back(distances[v]);
            }
        }
    }

    return {hhn_adj_list, hhn_weight_list};
}


/**
 * @brief Creates a path graph for a smaller h value from an existing path graph with larger h
 *
 * @details Given P_{\bullet \leq h_max}(G) and h ≤ h_max, constructs P_{\bullet \leq h}(G)
 *          using the stored hop counts in the source graph
 *
 * @param source_graph The source path graph computed with a larger h value
 * @param h The desired (smaller) number of hops
 * @return path_graph_t The path graph for the smaller h value
 *
 * @pre h must be less than or equal to the h value used to generate source_graph
 */
path_graph_t create_sub_path_graph(const path_graph_plus_t& source_graph, int h) {
    path_graph_t result;
    const int n_vertices = source_graph.adj_list.size();
    result.adj_list.resize(n_vertices);
    result.weight_list.resize(n_vertices);

    for (int v = 0; v < n_vertices; ++v) {
        const auto& src_adj = source_graph.adj_list[v];
        const auto& src_weights = source_graph.weight_list[v];
        const auto& src_hops = source_graph.hop_list[v];

        // Pre-allocate with conservative estimate
        result.adj_list[v].reserve(src_adj.size());
        result.weight_list[v].reserve(src_adj.size());

        // Filter edges based on hop count
        for (size_t i = 0; i < src_adj.size(); ++i) {
            if (src_hops[i] <= h) {
                result.adj_list[v].push_back(src_adj[i]);
                result.weight_list[v].push_back(src_weights[i]);

                // Only store path for smaller vertex index to larger
                if (v < src_adj[i]) {
                    result.shortest_paths[{v, src_adj[i]}] =
                        source_graph.shortest_paths.at({v, src_adj[i]});
                }
            }
        }
    }

    return result;
}

/**
 * @brief Converts a path_graph_plus_t to a path_graph_t by dropping the hop_list information
 *
 * @param plus_graph The augmented path graph with hop information
 * @return path_graph_t The basic path graph structure
 */
path_graph_t convert_to_path_graph(const path_graph_plus_t& plus_graph) {
    path_graph_t result;
    result.adj_list = plus_graph.adj_list;
    result.weight_list = plus_graph.weight_list;
    result.shortest_paths = plus_graph.shortest_paths;
    return result;
}

/**
 * @brief Generates a series of path graphs for different h values efficiently
 *
 * @param adj_list The adjacency list of the input graph
 * @param weight_list The weight list of the input graph
 * @param h_values Vector of h values to generate graphs for
 * @return std::vector<path_graph_t> Vector of path graphs for each requested h value
 *
 * @pre h_values must not be empty
 * @pre h_values must be sorted in ascending order for optimal performance
 */
std::vector<path_graph_t> create_path_graph_series(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<int>& h_values) {

    if (h_values.empty()) {
        Rf_error("h_values must not be empty");
    }

    // Use the maximum value in h_values as h_max
    const int h_max = *std::max_element(h_values.begin(), h_values.end());

    // First, create the graph for maximum h
    path_graph_plus_t max_graph = create_path_graph_plus(adj_list, weight_list, h_max);

    std::vector<path_graph_t> result;
    result.reserve(h_values.size());

    // Generate graphs for each requested h value
    for (int h : h_values) {
        if (h == h_max) {
            result.push_back(convert_to_path_graph(max_graph));
        } else {
            result.push_back(create_sub_path_graph(max_graph, h));
        }
    }

    return result;
}

/**
 * @brief Converts a path_graph_t struct to an R list representation
 *
 * @details Creates an R list containing the components of a path_graph_t struct,
 * converting all vertex indices from 0-based (C++) to 1-based (R) indexing.
 * The function handles memory protection for all created R objects.
 *
 * @param path_graph [in] Reference to a path_graph_t struct containing:
 *                       - adj_list: Vector of vectors containing adjacent vertices
 *                       - weight_list: Vector of vectors containing edge weights
 *                       - shortest_paths: Map of vertex pairs to path vectors
 *
 * @return SEXP A named R list with three components:
 *              - adj_list: List of integer vectors containing adjacent vertices
 *              - edge_length_list: List of numeric vectors with edge lengths
 *              - shortest_paths: List containing:
 *                * i: Integer vector of source vertices
 *                * j: Integer vector of target vertices
 *                * paths: List of integer vectors containing paths
 *
 * @note All vertex indices in the returned structure are converted to 1-based indexing
 *
 * @throw R_exception If memory allocation fails
 */
SEXP path_graph_from_path_graph_t(path_graph_t& path_graph) {

    std::vector<std::vector<int>>    path_graph_adj_vect                     = path_graph.adj_list;
    std::vector<std::vector<double>> path_graph_weight_vect                  = path_graph.weight_list;
    std::map<std::pair<int,int>, std::vector<int>> path_graph_shortest_paths = path_graph.shortest_paths;

    int n_vertices = static_cast<int>(path_graph_adj_vect.size());
    int nprot = 0; // Protection counter

    // Creating an R list for the path graph adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, path_graph_adj_vect[i].size()));
        int* A = INTEGER(RA);
        for (const auto& neighbor : path_graph_adj_vect[i])
            *A++ = neighbor + 1; // Converting to 1-based indexing
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Creating an R list for the path graph edge length list
    SEXP edge_length_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, path_graph_weight_vect[i].size()));
        double* D = REAL(RD);
        for (const auto& dist : path_graph_weight_vect[i])
            *D++ = dist;
        SET_VECTOR_ELT(edge_length_list, i, RD);
        UNPROTECT(1);
    }

    // Creating an R list for the path graph shortest paths map
    std::vector<std::pair<std::pair<int,int>, std::vector<int>>> path_pairs;
    path_pairs.reserve(path_graph_shortest_paths.size());

    // Convert map to vector for easier processing
    for (const auto& pair : path_graph_shortest_paths) {
        path_pairs.push_back(pair);
    }

    // Create R list for shortest paths
    SEXP shortest_paths = PROTECT(allocVector(VECSXP, 3)); nprot++;

    // First component: i coordinates
    SEXP i_coords = PROTECT(allocVector(INTSXP, path_pairs.size())); nprot++;
    int* i_ptr = INTEGER(i_coords);

    // Second component: j coordinates
    SEXP j_coords = PROTECT(allocVector(INTSXP, path_pairs.size())); nprot++;
    int* j_ptr = INTEGER(j_coords);

    // Third component: paths list
    SEXP paths = PROTECT(allocVector(VECSXP, path_pairs.size())); nprot++;

    for (size_t idx = 0; idx < path_pairs.size(); ++idx) {
        // Store i,j coordinates (add 1 for R indexing)
        i_ptr[idx] = path_pairs[idx].first.first + 1;
        j_ptr[idx] = path_pairs[idx].first.second + 1;

        // Create vector for this path
        const auto& path = path_pairs[idx].second;
        SEXP path_vec = PROTECT(allocVector(INTSXP, path.size()));
        int* path_ptr = INTEGER(path_vec);

        // Copy path values (add 1 for R indexing)
        for (size_t k = 0; k < path.size(); ++k) {
            path_ptr[k] = path[k] + 1;
        }

        SET_VECTOR_ELT(paths, idx, path_vec);
        UNPROTECT(1); // path_vec
    }

    // Set components of shortest_paths list
    SET_VECTOR_ELT(shortest_paths, 0, i_coords);
    SET_VECTOR_ELT(shortest_paths, 1, j_coords);
    SET_VECTOR_ELT(shortest_paths, 2, paths);

    // Names for shortest_paths components
    SEXP sp_names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(sp_names, 0, mkChar("i"));
    SET_STRING_ELT(sp_names, 1, mkChar("j"));
    SET_STRING_ELT(sp_names, 2, mkChar("paths"));
    setAttrib(shortest_paths, R_NamesSymbol, sp_names);


    // Creating the final result list
    SEXP res = PROTECT(allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, edge_length_list);
    SET_VECTOR_ELT(res, 2, shortest_paths);

    SEXP names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("edge_length_list"));
    SET_STRING_ELT(names, 2, mkChar("shortest_paths"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}


/**
 * @brief Creates a series of path graphs for different hop limits
 *
 * @details This function computes multiple path graphs from the same input graph,
 * each with a different maximum hop limit specified in h_values. It efficiently
 * computes these graphs by leveraging the maximum hop limit (h_max) to avoid
 * redundant calculations.
 *
 * @param s_adj_list [in] SEXP representing an R list where each element i contains
 *                        an integer vector of vertices adjacent to vertex i.
 *                        Vertices are 1-based in R format.
 *
 * @param s_weight_list [in] SEXP representing an R list where each element i contains
 *                          a numeric vector of edge weights corresponding to the
 *                          adjacencies in s_adj_list.
 *
 * @param s_h_values [in] SEXP representing an integer vector of hop limits for which to compute path graphs.
 *
 * @return SEXP A list where each element corresponds to a path graph computed for
 *              the corresponding hop limit in h_values. Each path graph is itself
 *              a named list with three components:
 *              - adj_list: List of integer vectors containing reachable vertices
 *              - edge_length_list: List of numeric vectors with path lengths
 *              - shortest_paths: List containing path information
 *
 * @note All vertex indices in the returned structure are 1-based (R convention)
 *
 * @throw R_exception If memory allocation fails or input validation fails
 *
 * @see create_path_graph_series
 * @see path_graph_from_path_graph_t
 */
SEXP S_create_path_graph_series(SEXP s_adj_list,
                                SEXP s_weight_list,
                                SEXP s_h_values) {
    // Converting R inputs to C++ types
    std::vector<std::vector<int>> adj_vect       = Rgraph_to_vector(s_adj_list);
    std::vector<std::vector<double>> weight_vect = Rweights_to_vector(s_weight_list);

    int nprot = 0; // Protection counter

    PROTECT(s_h_values = coerceVector(s_h_values, INTSXP)); ++nprot;
    int* h_values_array = INTEGER(s_h_values);
    int h_values_len = LENGTH(s_h_values);
    std::vector<int> h_values(h_values_array, h_values_array + h_values_len);

    std::vector<path_graph_t> graph_series = create_path_graph_series(adj_vect, weight_vect, h_values);

    SEXP res = PROTECT(allocVector(VECSXP, h_values_len)); nprot++;
    for (int i = 0; i < h_values_len; ++i) {
        SEXP path_graph = path_graph_from_path_graph_t(graph_series[i]);
        SET_VECTOR_ELT(res, i, path_graph);
    }

    UNPROTECT(nprot);
    return res;
}


/**
 * @brief Creates a path graph structure from adjacency and weight lists
 *
 * @details This function converts a graph representation from R to C++, processes it using create_path_graph_plus(),
 * and returns a complex R list structure containing the path graph information. The function handles
 * the conversion between 0-based (C++) and 1-based (R) indexing.
 *
 * @param s_adj_list [in] SEXP representing an R list where each element i contains an integer vector
 *                        of vertices adjacent to vertex i. Vertices are 1-based in R format.
 *
 * @param s_edge_length_list [in] SEXP representing an R list where each element i contains a numeric vector
 *                                of edge lengths corresponding to the adjacencies in s_adj_list.
 *                                Must have the same structure as s_adj_list.
 *
 * @param s_h [in] SEXP representing an integer scalar specifying the maximum number of hops allowed
 *                 in the path graph.
 *
 * @return SEXP A named R list with four components:
 *         - adj_list: List of integer vectors, where adj_list[[i]] contains vertices reachable from i
 *         - edge_length_list: List of numeric vectors with corresponding path lengths
 *         - hop_list: List of integer vectors with number of hops for each path
 *         - shortest_paths: List with three components:
 *           * i: Integer vector of source vertices
 *           * j: Integer vector of destination vertices
 *           * paths: List of integer vectors containing the vertex sequences for each path
 *         All vertex indices in the returned structure are 1-based (R convention)
 *
 * @note Memory management is handled using PROTECT/UNPROTECT with a counter (nprot)
 *       to ensure proper cleanup in case of errors.
 *
 * @throw R_exception If memory allocation fails or input validation fails
 *
 * @see create_path_graph_plus
 * @see Rgraph_to_vector
 * @see Rweights_to_vector
 */
SEXP S_create_path_graph_plus(SEXP s_adj_list,
                              SEXP s_edge_length_list,
                              SEXP s_h) {

    // Converting R inputs to C++ types
    std::vector<std::vector<int>> adj_vect = Rgraph_to_vector(s_adj_list);
    std::vector<std::vector<double>> weight_vect = Rweights_to_vector(s_edge_length_list);
    int h = INTEGER(s_h)[0];

    path_graph_plus_t path_graph = create_path_graph_plus(adj_vect, weight_vect, h);

    // Extracting results
    std::vector<std::vector<int>>    path_graph_adj_vect                     = path_graph.adj_list;
    std::vector<std::vector<double>> path_graph_weight_vect                  = path_graph.weight_list;
    std::vector<std::vector<int>>    path_graph_hop_vect                     = path_graph.hop_list;
    std::map<std::pair<int,int>, std::vector<int>> path_graph_shortest_paths = path_graph.shortest_paths;

    int n_vertices = static_cast<int>(path_graph_adj_vect.size());
    int nprot = 0; // Protection counter

    // Creating an R list for the path graph adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, path_graph_adj_vect[i].size()));
        int* A = INTEGER(RA);
        for (const auto& neighbor : path_graph_adj_vect[i])
            *A++ = neighbor + 1; // Converting to 1-based indexing
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Creating an R list for the path graph edge length list
    SEXP edge_length_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, path_graph_weight_vect[i].size()));
        double* D = REAL(RD);
        for (const auto& dist : path_graph_weight_vect[i])
            *D++ = dist;
        SET_VECTOR_ELT(edge_length_list, i, RD);
        UNPROTECT(1);
    }

    // Creating an R list for the path graph hop list
    SEXP hop_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, path_graph_hop_vect[i].size()));
        int* A = INTEGER(RA);
        for (const auto& hops : path_graph_hop_vect[i])
            *A++ = hops;  // No conversion needed - hop counts are already 1-based
        SET_VECTOR_ELT(hop_list, i, RA);
        UNPROTECT(1);
    }

    // Creating an R list for the path graph shortest paths map
    std::vector<std::pair<std::pair<int,int>, std::vector<int>>> path_pairs;
    path_pairs.reserve(path_graph_shortest_paths.size());

    // Convert map to vector for easier processing
    for (const auto& pair : path_graph_shortest_paths) {
        path_pairs.push_back(pair);
    }

    // Create R list for shortest paths
    SEXP shortest_paths = PROTECT(allocVector(VECSXP, 3)); nprot++;

    // First component: i coordinates
    SEXP i_coords = PROTECT(allocVector(INTSXP, path_pairs.size())); nprot++;
    int* i_ptr = INTEGER(i_coords);

    // Second component: j coordinates
    SEXP j_coords = PROTECT(allocVector(INTSXP, path_pairs.size())); nprot++;
    int* j_ptr = INTEGER(j_coords);

    // Third component: paths list
    SEXP paths = PROTECT(allocVector(VECSXP, path_pairs.size())); nprot++;

    for (size_t idx = 0; idx < path_pairs.size(); ++idx) {
        // Store i,j coordinates (add 1 for R indexing)
        i_ptr[idx] = path_pairs[idx].first.first + 1;
        j_ptr[idx] = path_pairs[idx].first.second + 1;

        // Create vector for this path
        const auto& path = path_pairs[idx].second;
        SEXP path_vec = PROTECT(allocVector(INTSXP, path.size()));
        int* path_ptr = INTEGER(path_vec);

        // Copy path values (add 1 for R indexing)
        for (size_t k = 0; k < path.size(); ++k) {
            path_ptr[k] = path[k] + 1;
        }

        SET_VECTOR_ELT(paths, idx, path_vec);
        UNPROTECT(1); // path_vec
    }

    // Set components of shortest_paths list
    SET_VECTOR_ELT(shortest_paths, 0, i_coords);
    SET_VECTOR_ELT(shortest_paths, 1, j_coords);
    SET_VECTOR_ELT(shortest_paths, 2, paths);

    // Names for shortest_paths components
    SEXP sp_names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(sp_names, 0, mkChar("i"));
    SET_STRING_ELT(sp_names, 1, mkChar("j"));
    SET_STRING_ELT(sp_names, 2, mkChar("paths"));
    setAttrib(shortest_paths, R_NamesSymbol, sp_names);


    // Creating the final result list
    SEXP res = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, edge_length_list);
    SET_VECTOR_ELT(res, 2, hop_list);
    SET_VECTOR_ELT(res, 3, shortest_paths);

    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("edge_length_list"));
    SET_STRING_ELT(names, 2, mkChar("hop_list"));
    SET_STRING_ELT(names, 3, mkChar("shortest_paths"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

/**
 * @brief Creates a path graph structure optimized for path linear model computation
 *
 * @details This function constructs a path_graph_plm_t structure from an input graph
 * represented by adjacency and weight lists. The resulting structure contains:
 * - h-hop neighborhood information for each vertex
 * - All shortest paths between connected vertices
 * - For each vertex, information about all paths containing it
 *
 * The algorithm uses a modified BFS (Breadth-First Search) approach to:
 * 1. Find all vertices within h hops of each start vertex
 * 2. Compute shortest paths between connected vertices
 * 3. Track distances and hop counts
 * 4. Build vertex-to-path mapping for efficient local linear model computation
 *
 * @param adj_list Vector of adjacency lists where adj_list[i] contains indices
 *                 of vertices adjacent to vertex i
 * @param weight_list Vector of weight lists where weight_list[i] contains weights
 *                    of edges incident to vertex i. Must match adj_list structure.
 * @param h Maximum number of hops to consider when building the path graph
 *
 * @return path_graph_plm_t A structure containing:
 *         - adj_list: h-hop neighborhood adjacency lists
 *         - weight_list: accumulated weights to h-hop neighbors
 *         - hop_list: hop counts to neighbors
 *         - shortest_paths: map from vertex pairs to shortest paths between them
 *         - vertex_paths: for each vertex, lists of paths containing it and positions within them
 *
 * @throws std::invalid_argument If adj_list and weight_list sizes don't match
 * @throws std::invalid_argument If h < 0
 *
 * @pre adj_list and weight_list must have the same size
 * @pre For each vertex i, adj_list[i] and weight_list[i] must have the same size
 * @pre h must be non-negative
 *
 * @note Time Complexity: O(V * (V + E)) where V is number of vertices and E is number of edges
 * @note Space Complexity: O(V^2) in worst case for storing all paths
 *
 * Example usage:
 * @code
 *     std::vector<std::vector<int>> adj_list = {{1,2}, {0,2}, {0,1}};
 *     std::vector<std::vector<double>> weight_list = {{1.0,1.0}, {1.0,1.0}, {1.0,1.0}};
 *     int h = 2;
 *     auto path_graph = create_path_graph_plm(adj_list, weight_list, h);
 * @endcode
 *
 * @see path_graph_plm_t
 * @see vertex_path_info_t
 */
path_graph_plm_t create_path_graph_plm(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    #define DEBUG_PLM_GRAPH 0

    const int n_vertices = adj_list.size();

    #if DEBUG_PLM_GRAPH
    Rprintf("\nStarting path graph creation\n");
    Rprintf("Number of vertices: %d\n", n_vertices);

    // Print input adjacency list and weights
    Rprintf("\nInput adjacency lists:\n");
    for (int i = 0; i < n_vertices; ++i) {
        Rprintf("Vertex %d connected to:", i);
        for (int j : adj_list[i]) {
            Rprintf(" %d", j);
        }
        Rprintf("\n");
    }
    #endif

    path_graph_plm_t result;
    result.h = h;
    result.adj_list.resize(n_vertices);
    result.weight_list.resize(n_vertices);
    result.hop_list.resize(n_vertices);

    // Pre-allocate vectors to avoid reallocations
    std::vector<double> distances(n_vertices);
    std::vector<int> hops(n_vertices);
    std::vector<int> parent(n_vertices);
    std::vector<bool> in_queue(n_vertices);

    // Temporary storage for neighbors with their properties
    struct neighbor_info_t {
        int vertex;
        double distance;
        int hops;
        neighbor_info_t(int v, double d, int h) : vertex(v), distance(d), hops(h) {}
    };
    std::vector<neighbor_info_t> current_neighbors;
    current_neighbors.reserve(n_vertices / 2);  // Conservative estimate

    for (int start = 0; start < n_vertices; ++start) {
        #if DEBUG_PLM_GRAPH
        Rprintf("\nProcessing paths from vertex %d\n", start);
        #endif

        // Reset arrays (faster than reinitializing)
        std::fill(distances.begin(), distances.end(), std::numeric_limits<double>::infinity());
        std::fill(hops.begin(), hops.end(), std::numeric_limits<int>::max());
        std::fill(parent.begin(), parent.end(), -1);
        std::fill(in_queue.begin(), in_queue.end(), false);
        current_neighbors.clear();

        // Initialize start vertex
        distances[start] = 0;
        hops[start] = 0;

        // Use deque for better cache locality
        std::deque<int> q{start};
        in_queue[start] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop_front();
            in_queue[current] = false;

            // Early termination if we've reached hop limit
            if (hops[current] >= h) {
                #if DEBUG_PLM_GRAPH
                Rprintf("  Reached hop limit at vertex %d\n", current);
                #endif
                continue;
            }

            // Cache current values
            const double current_dist = distances[current];
            const int current_hops = hops[current];
            const int next_hops = current_hops + 1;

            // Process neighbors
            const auto& current_adj = adj_list[current];
            const auto& current_weights = weight_list[current];
            const size_t n_neighbors = current_adj.size();

            for (size_t i = 0; i < n_neighbors; ++i) {
                const int neighbor = current_adj[i];

                if (neighbor == current) continue;  // Skip self-loops

                const double new_distance = current_dist + current_weights[i];

                #if DEBUG_PLM_GRAPH
                Rprintf("    Checking neighbor %d (new_distance=%.2f, current_best=%.2f)\n",
                        neighbor, new_distance, distances[neighbor]);
                #endif

                // Update if we found a shorter path or same distance with fewer hops
                if (new_distance < distances[neighbor]) {
                    #if DEBUG_PLM_GRAPH
                    Rprintf("    Updating path to %d: distance %.2f->%.2f, hops %d->%d\n",
                            neighbor, distances[neighbor], new_distance,
                            hops[neighbor], next_hops);
                    #endif

                    distances[neighbor] = new_distance;
                    hops[neighbor] = next_hops;
                    parent[neighbor] = current;

                    // Add to queue if not already in it
                    if (!in_queue[neighbor] && next_hops <= h) {
                        q.push_back(neighbor);
                        in_queue[neighbor] = true;
                    }
                }
            }
        }

        // Collect neighbors and their properties
        for (int v = 0; v < n_vertices; ++v) {
            if (v != start && hops[v] <= h) {
                current_neighbors.emplace_back(v, distances[v], hops[v]);

                // Store path for pairs where start < v
                if (start < v) {
                    #if DEBUG_PLM_GRAPH
                    Rprintf("\nStoring path from %d to %d:\n", start, v);
                    #endif

                    auto& path = result.shortest_paths[{start, v}];
                    path.clear();
                    path.reserve(hops[v] + 1);  // Exact size known

                    // Reconstruct path
                    // Print path reconstruction
                    #if DEBUG_PLM_GRAPH
                    Rprintf("  Path: ");
                    #endif
                    for (int current = v; current != -1; current = parent[current]) {
                        path.push_back(current);
                        #if DEBUG_PLM_GRAPH
                        Rprintf("%d ", current);
                        #endif
                    }

                    std::reverse(path.begin(), path.end());

                    #if DEBUG_PLM_GRAPH
                    Rprintf("\n  Final path after reverse: ");
                    for (int vertex : path) {
                        Rprintf("%d ", vertex);
                    }
                    Rprintf("\n");
                    #endif
                }
            }

            #if DEBUG_PLM_GRAPH
            Rprintf("\nPaths found from vertex %d:\n", start);
            for (const auto& info : current_neighbors) {
                Rprintf("  To %d: distance=%.2f, hops=%d\n",
                        info.vertex, info.distance, info.hops);
            }
            #endif
        }

        // Bulk insert neighbors and their properties
        const size_t n_new_neighbors = current_neighbors.size();
        result.adj_list[start].resize(n_new_neighbors);
        result.weight_list[start].resize(n_new_neighbors);
        result.hop_list[start].resize(n_new_neighbors);

        for (size_t i = 0; i < n_new_neighbors; ++i) {
            const auto& info = current_neighbors[i];
            result.adj_list[start][i] = info.vertex;
            result.weight_list[start][i] = info.distance;
            result.hop_list[start][i] = info.hops;
        }
    }

    result.vertex_paths.resize(n_vertices);

    // Update vertex_paths with endpoint pairs instead of full paths
    for (const auto& path_entry : result.shortest_paths) {
        const auto& endpoints = path_entry.first;
        const auto& path = path_entry.second;
        for (size_t pos = 0; pos < path.size(); ++pos) {
            int vertex = path[pos];
            result.vertex_paths[vertex].containing_paths.push_back(endpoints);
            result.vertex_paths[vertex].position_in_path.push_back(pos);
        }
    }

    return result;
}


/**
 * @brief Creates an R list representation of a path_graph_plm_t object
 *
 * This function creates a path graph PLM (Path Length Matrix) representation that can be used in R.
 * It extends the functionality of S_create_path_graph_plus by adding vertex path information.
 *
 * @param s_adj_list SEXP containing the adjacency list representation of the graph.
 *        Each element i contains the vertices adjacent to vertex i (1-based indexing in R).
 * @param s_edge_length_list SEXP containing the edge weights corresponding to the adjacency list.
 *        Each element i contains the weights of edges connected to vertex i.
 * @param s_h SEXP containing an integer specifying the maximum path length to consider.
 *
 * @return SEXP containing a list with the following components:
 *   - adj_list: List of integer vectors representing the path graph's adjacency list
 *   - edge_length_list: List of numeric vectors containing edge weights
 *   - hop_list: List of integer vectors containing hop counts
 *   - shortest_paths: List with components:
 *     - i: Integer vector of start vertices
 *     - j: Integer vector of end vertices
 *     - paths: List of integer vectors representing paths
 *   - vertex_paths: List of matrices, one per vertex, where each matrix has columns:
 *     - start: Integer vector of path start vertices
 *     - end: Integer vector of path end vertices
 *     - position: Integer vector of vertex positions in paths
 *
 * @note All vertex indices in the returned object use 1-based indexing (R convention)
 * @note The function handles memory protection internally using PROTECT/UNPROTECT
 * @note This function requires the graph to be undirected and connected
 *
 * @see S_create_path_graph_plus
 * @see path_graph_plm_t
 */
SEXP S_create_path_graph_plm(SEXP s_adj_list,
                             SEXP s_edge_length_list,
                             SEXP s_h) {
    // Converting R inputs to C++ types
    std::vector<std::vector<int>> adj_vect = Rgraph_to_vector(s_adj_list);
    std::vector<std::vector<double>> weight_vect = Rweights_to_vector(s_edge_length_list);
    int h = INTEGER(s_h)[0];

    // Create path graph PLM
    path_graph_plm_t path_graph = create_path_graph_plm(adj_vect, weight_vect, h);

    // First, create the base path_graph_plus components using the existing function
    SEXP base_components = S_create_path_graph_plus(s_adj_list, s_edge_length_list, s_h);
    int nprot = 0;
    PROTECT(base_components); nprot++;

    // Now create the vertex_paths component
    int n_vertices = path_graph.vertex_paths.size();
    SEXP vertex_paths_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    // For each vertex
    for (int v = 0; v < n_vertices; v++) {
        const auto& vertex_info = path_graph.vertex_paths[v];
        int n_paths = vertex_info.containing_paths.size();

        // Create matrix with 3 columns: start, end, position_in_path
        SEXP path_matrix = PROTECT(allocMatrix(INTSXP, n_paths, 3));
        int* matrix_ptr = INTEGER(path_matrix);

        // Fill the matrix
        for (int i = 0; i < n_paths; i++) {
            // start vertex (adding 1 for R's 1-based indexing)
            matrix_ptr[i] = vertex_info.containing_paths[i].first + 1;
            // end vertex
            matrix_ptr[i + n_paths] = vertex_info.containing_paths[i].second + 1;
            // position in path
            matrix_ptr[i + 2 * n_paths] = vertex_info.position_in_path[i] + 1;
        }

        // Set column names
        SEXP colnames = PROTECT(allocVector(STRSXP, 3));
        SET_STRING_ELT(colnames, 0, mkChar("start"));
        SET_STRING_ELT(colnames, 1, mkChar("end"));
        SET_STRING_ELT(colnames, 2, mkChar("position"));
        setAttrib(path_matrix, R_DimNamesSymbol,
                 PROTECT(list2(R_NilValue, colnames)));

        SET_VECTOR_ELT(vertex_paths_list, v, path_matrix);
        UNPROTECT(3); // path_matrix, colnames, and dimnames list
    }

    // Create the final result list by combining base_components and vertex_paths
    SEXP res = PROTECT(allocVector(VECSXP, 5)); nprot++;

    // Copy all elements from base_components
    for (int i = 0; i < 4; i++) {
        SET_VECTOR_ELT(res, i, VECTOR_ELT(base_components, i));
    }

    // Add vertex_paths as the fifth element
    SET_VECTOR_ELT(res, 4, vertex_paths_list);

    // Set names for the result list
    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;

    // Copy names from base_components
    SEXP base_names = getAttrib(base_components, R_NamesSymbol);
    for (int i = 0; i < 4; i++) {
        SET_STRING_ELT(names, i, STRING_ELT(base_names, i));
    }

    // Add name for vertex_paths
    SET_STRING_ELT(names, 4, mkChar("vertex_paths"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}

/**
 * @brief Converts an R list representation of a path graph PLM to a C++ path_graph_plm_t object
 *
 * This function takes the R list output from S_create_path_graph_plm and reconstructs
 * the corresponding C++ path_graph_plm_t object. It handles the conversion of all
 * components including adjacency lists, weights, hop counts, shortest paths, and
 * vertex path information.
 *
 * @param s_path_graph SEXP containing the R list representation of the path graph PLM.
 *        Must contain the following named components:
 *        - adj_list: List of integer vectors (adjacency lists)
 *        - edge_length_list: List of numeric vectors (edge weights)
 *        - hop_list: List of integer vectors (hop counts)
 *        - shortest_paths: List containing path information
 *        - vertex_paths: List of matrices, where each matrix has columns:
 *          - start: Integer vector of path start vertices
 *          - end: Integer vector of path end vertices
 *          - position: Integer vector of vertex positions in paths
 *
 * @return path_graph_plm_t containing the converted C++ representation
 *
 * @note The function expects the input to use R's 1-based indexing and converts
 *       it to C++'s 0-based indexing
 * @note The function assumes the input SEXP is properly protected by the caller
 * @note If the input format is invalid, the function will throw a std::runtime_error
 *
 * @see S_create_path_graph_plm
 * @see path_graph_plm_t
 */
path_graph_plm_t sexp_to_path_graph_plm(SEXP s_path_graph) {
    path_graph_plm_t result;

    // Verify that input is a list with required components
    if (!isNewList(s_path_graph)) {
        Rf_error("Input must be an R list");
    }

    // Get list components
    SEXP s_adj_list = VECTOR_ELT(s_path_graph, 0);
    SEXP s_weight_list = VECTOR_ELT(s_path_graph, 1);
    SEXP s_hop_list = VECTOR_ELT(s_path_graph, 2);
    SEXP s_shortest_paths = VECTOR_ELT(s_path_graph, 3);
    SEXP s_vertex_paths = VECTOR_ELT(s_path_graph, 4);

    // First initialize h (not available in R structure, will be derived from paths)
    result.h = 0;  // Will be updated when processing paths

    // Convert adjacency list
    int n_vertices = Rf_length(s_adj_list);
    result.adj_list.resize(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        SEXP adj_vec = VECTOR_ELT(s_adj_list, i);
        int* adj_ptr = INTEGER(adj_vec);
        int adj_len = Rf_length(adj_vec);

        result.adj_list[i].resize(adj_len);
        for (int j = 0; j < adj_len; j++) {
            result.adj_list[i][j] = adj_ptr[j] - 1;  // Convert to 0-based indexing
        }
    }

    // Convert weight list
    result.weight_list.resize(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        SEXP weight_vec = VECTOR_ELT(s_weight_list, i);
        double* weight_ptr = REAL(weight_vec);
        int weight_len = Rf_length(weight_vec);

        result.weight_list[i].resize(weight_len);
        for (int j = 0; j < weight_len; j++) {
            result.weight_list[i][j] = weight_ptr[j];
        }
    }

    // Convert hop list
    result.hop_list.resize(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        SEXP hop_vec = VECTOR_ELT(s_hop_list, i);
        int* hop_ptr = INTEGER(hop_vec);
        int hop_len = Rf_length(hop_vec);

        result.hop_list[i].resize(hop_len);
        for (int j = 0; j < hop_len; j++) {
            result.hop_list[i][j] = hop_ptr[j];  // Hop counts don't need index conversion
        }
    }

    // Convert shortest paths
    SEXP s_i_coords = VECTOR_ELT(s_shortest_paths, 0);
    SEXP s_j_coords = VECTOR_ELT(s_shortest_paths, 1);
    SEXP s_paths = VECTOR_ELT(s_shortest_paths, 2);

    int n_paths = Rf_length(s_i_coords);
    int* i_ptr = INTEGER(s_i_coords);
    int* j_ptr = INTEGER(s_j_coords);

    for (int idx = 0; idx < n_paths; idx++) {
        int i = i_ptr[idx] - 1;  // Convert to 0-based indexing
        int j = j_ptr[idx] - 1;  // Convert to 0-based indexing

        SEXP path_vec = VECTOR_ELT(s_paths, idx);
        int* path_ptr = INTEGER(path_vec);
        int path_len = Rf_length(path_vec);

        // Update h if we find a longer path
        result.h = std::max(result.h, path_len - 1);

        std::vector<int> path(path_len);
        for (int k = 0; k < path_len; k++) {
            path[k] = path_ptr[k] - 1;  // Convert to 0-based indexing
        }

        result.shortest_paths[{std::min(i, j), std::max(i, j)}] = path;
    }

    // Convert vertex paths
    result.vertex_paths.resize(n_vertices);
    for (int v = 0; v < n_vertices; v++) {
        SEXP path_matrix = VECTOR_ELT(s_vertex_paths, v);
        int* matrix_ptr = INTEGER(path_matrix);
        int n_paths = nrows(path_matrix);

        // Get pointers to each column
        int* start_ptr = matrix_ptr;
        int* end_ptr = matrix_ptr + n_paths;
        int* pos_ptr = matrix_ptr + 2 * n_paths;

        // Resize the containing_paths and position_in_path vectors
        result.vertex_paths[v].containing_paths.resize(n_paths);
        result.vertex_paths[v].position_in_path.resize(n_paths);

        // Fill the vectors
        for (int i = 0; i < n_paths; i++) {
            // Convert to 0-based indexing for vertices
            result.vertex_paths[v].containing_paths[i] = {
                start_ptr[i] - 1,
                end_ptr[i] - 1
            };
            result.vertex_paths[v].position_in_path[i] = pos_ptr[i] - 1;
        }
    }

    return result;
}
