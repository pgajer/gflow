#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"

#include <vector>
#include <queue>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <utility>

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h);
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
 * @Rf_warning If k is greater than or equal to the diameter of the graph, the resulting hHN graph may be complete.
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


struct hvertex_t { // neighbor of a given vertex struct
    int idx;       // index of the neighbor
    double dist;   // distance from the given vertex
    double hops;   // number of hops from the given vertex
};

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> create_hHN_graph_hvertex_vector(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    const double EPSILON = 1e-9;  // Define epsilon for floating-point comparisons
    int n_vertices = adj_list.size();
    std::vector<std::vector<int>> hhn_adj_list(n_vertices);
    std::vector<std::vector<double>> hhn_weight_list(n_vertices);
    hvertex_t hvertex;

    for (int start = 0; start < n_vertices; ++start) {
        std::vector<hvertex_t> neighbors;
        std::queue<hvertex_t> q;
        hvertex.idx = start;
        hvertex.dist = 0;
        hvertex.hops = 0;
        q.push(hvertex);

        while (!q.empty()) {
            hvertex = q.front();
            int current = hvertex.idx;
            q.pop();

            if (hvertex.hops == h) continue;

            bool all_neighbors_processed = true;
            for (size_t i = 0; i < adj_list[current].size(); ++i) {
                int neighbor = adj_list[current][i];

                if (neighbor == current) continue;  // Skip self-loops

                double new_distance = hvertex.dist + weight_list[current][i];

                auto it = std::find_if(neighbors.begin(), neighbors.end(),
                    [neighbor](const hvertex_t& v) { return v.idx == neighbor; });

                if (it == neighbors.end()) {
                    // Neighbor not found, add it
                    hvertex_t new_neighbor;
                    new_neighbor.idx = neighbor;
                    new_neighbor.dist = new_distance;
                    new_neighbor.hops = hvertex.hops + 1;
                    q.push(new_neighbor);
                    neighbors.push_back(new_neighbor);
                    all_neighbors_processed = false;
                } else if (new_distance < it->dist - EPSILON) {
                    // Update existing neighbor if new distance is shorter
                    it->dist = new_distance;
                    it->hops = hvertex.hops + 1;
                    q.push(*it);
                    all_neighbors_processed = false;
                }
            }

            // Early termination if all neighbors have been processed
            if (all_neighbors_processed && hvertex.hops == h - 1) {
                break;
            }
        }

        for (const auto& neighbor : neighbors) {
            hhn_adj_list[start].push_back(neighbor.idx);
            hhn_weight_list[start].push_back(neighbor.dist);
        }
    }

    return {hhn_adj_list, hhn_weight_list};
}


std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> create_hHN_graph_hvertex_hashmap(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int h) {

    const double EPSILON = 1e-9;
    int n_vertices = adj_list.size();
    std::vector<std::vector<int>> hhn_adj_list(n_vertices);
    std::vector<std::vector<double>> hhn_weight_list(n_vertices);
    hvertex_t hvertex;

    for (int start = 0; start < n_vertices; ++start) {
        std::unordered_map<int, hvertex_t> neighbors;
        std::queue<hvertex_t> q;
        hvertex.idx = start;
        hvertex.dist = 0;
        hvertex.hops = 0;
        q.push(hvertex);
        neighbors[start] = hvertex;

        while (!q.empty()) {
            hvertex = q.front();
            int current = hvertex.idx;
            q.pop();

            if (hvertex.hops == h) continue;

            bool all_neighbors_processed = true;
            for (size_t i = 0; i < adj_list[current].size(); ++i) {
                int neighbor = adj_list[current][i];

                if (neighbor == current) continue;  // Skip self-loops

                double new_distance = hvertex.dist + weight_list[current][i];

                auto it = neighbors.find(neighbor);

                if (it == neighbors.end()) {
                    // Neighbor not found, add it
                    hvertex_t new_neighbor;
                    new_neighbor.idx = neighbor;
                    new_neighbor.dist = new_distance;
                    new_neighbor.hops = hvertex.hops + 1;
                    q.push(new_neighbor);
                    neighbors[neighbor] = new_neighbor;
                    all_neighbors_processed = false;
                } else if (new_distance < it->second.dist - EPSILON) {
                    // Update existing neighbor if new distance is shorter
                    it->second.dist = new_distance;
                    it->second.hops = hvertex.hops + 1;
                    q.push(it->second);
                    all_neighbors_processed = false;
                }
            }

            // Early termination if all neighbors have been processed
            if (all_neighbors_processed && hvertex.hops == h - 1) {
                break;
            }
        }

        for (const auto& [idx, neighbor] : neighbors) {
            if (idx != start) {
                hhn_adj_list[start].push_back(idx);
                hhn_weight_list[start].push_back(neighbor.dist);
            }
        }
    }

    return {hhn_adj_list, hhn_weight_list};
}



/**
 * @brief R interface for creating a k-hop neighborhood (hHN) graph from an input graph.
 *
 * This function serves as a wrapper for the C++ implementation of create_hHN_graph,
 * allowing it to be called from R. It handles the conversion of R data structures to
 * C++ and back, creating a k-hop neighborhood graph based on the input adjacency list
 * and edge weights.
 *
 * @param s_adj_list An R list of integer vectors representing the adjacency list of the input graph.
 *                   Each element s_adj_list[[i]] contains the indices of vertices adjacent to vertex i.
 *                   Note: R uses 1-based indexing, which is converted to 0-based indexing for C++.
 * @param s_weight_list An R list of double vectors representing the edge weights of the input graph.
 *                      s_weight_list[[i]][j] is the weight of the edge from vertex i to its j-th neighbor.
 * @param s_h An R integer scalar representing the maximum number of hops to consider for neighborhood connectivity.
 *
 * @return An R list containing two elements:
 *         - adj_list: A list of integer vectors representing the adjacency list of the hHN graph.
 *                     Indices in this list are 1-based to conform with R's indexing.
 *         - dist_list: A list of double vectors representing the edge weights (distances) of the hHN graph.
 *
 * @note This function uses the R API and handles memory management with PROTECT/UNPROTECT.
 * @note The time complexity is O(n * (m + n log n)), where n is the number of vertices and m is the number of edges.
 * @note The space complexity is O(n^2) in the worst case, where the hHN graph becomes a complete graph.
 *
 * @Rf_warning Ensure that s_k is a scalar integer and non-negative. The function currently does not perform
 *          extensive Rf_error checking on inputs.
 *
 * @see create_hHN_graph for the underlying C++ implementation.
 *
 * @example
 * In R:
 * result <- .Call("S_create_hHN_graph", adj_list, weight_list, h)
 * hhn_adj_list <- result$adj_list
 * hhn_dist_list <- result$dist_list
 */
SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h) {
    // Correct: Converts R inputs to C++ types
    std::vector<std::vector<int>> adj_vv       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vv = convert_weight_list_from_R(s_weight_list);

    // Correct: Coerces s_k to integer and extracts its value
    PROTECT(s_h = Rf_coerceVector(s_h, INTSXP));
    int h = INTEGER(s_h)[0];
    UNPROTECT(1);

    #if 0
    PROTECT(s_type = Rf_coerceVector(s_type, INTSXP));
    int type = INTEGER(s_type)[0];
    UNPROTECT(1);

    // benchmark_hHN_graph results:
    // type0   2.355964   2.406112   2.454287   2.420903   2.433518   2.654935     5
    // type1 215.185024 220.344316 222.421501 220.622162 227.920013 228.035989     5
    // type2   8.643084   8.769095   8.790769   8.794763   8.836087   8.910816     5

    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> hhn_graph;
    switch(type) {
    case 0:
        hhn_graph = create_hHN_graph(adj_vv, weight_vv, h);
        break;
    case 1:
        hhn_graph = create_hHN_graph_hvertex_vector(adj_vv, weight_vv, h);
        break;
    case 2:
        hhn_graph = create_hHN_graph_hvertex_hashmap(adj_vv, weight_vv, h);
        break;
    default:
        Rf_error("Invalid type specified.");
    }
    #endif

    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> hhn_graph = create_hHN_graph(adj_vv, weight_vv, h);

    // Correct: Extracts results
    std::vector<std::vector<int>> hhn_adj_vv       = hhn_graph.first;
    std::vector<std::vector<double>> hhn_weight_vv = hhn_graph.second;

    // Correct: Prepares R list for adjacency list
    int n_vertices = static_cast<int>(hhn_adj_vv.size());
    SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(Rf_allocVector(INTSXP, hhn_adj_vv[i].size()));
        int* A = INTEGER(RA);
        for (const auto& neighbor : hhn_adj_vv[i])
            *A++ = neighbor + 1;  // Correct: Adjusts for R's 1-based indexing
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Correct: Prepares R list for distance list
    SEXP dist_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(Rf_allocVector(REALSXP, hhn_weight_vv[i].size()));
        double* D = REAL(RD);
        for (const auto& dist : hhn_weight_vv[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
        UNPROTECT(1);
    }

    // Correct: Prepares the final result list
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, dist_list);

    UNPROTECT(2); // Unprotect adj_list and dist_list

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("dist_list"));
    Rf_setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(2); // Unprotect res and names

    return res;
}
