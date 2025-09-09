#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <utility>

#include "SEXP_cpp_conversion_utils.h"
#include "cpp_utils.h"

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

    #define DEBUG__create_hHN_graph 0
    #if DEBUG__create_hHN_graph
    Rprintf("Starting create_hHN_graph with h = %d\n", h);
    Rprintf("Number of vertices: %d\n", n_vertices);
    #endif

    std::vector<std::vector<int>> hhn_adj_list(n_vertices);
    std::vector<std::vector<double>> hhn_weight_list(n_vertices);

    // Pre-allocate vectors to avoid reallocations
    std::vector<double> distances(n_vertices);
    std::vector<int> hops(n_vertices);
    std::vector<bool> in_queue(n_vertices);

    for (int start = 0; start < n_vertices; ++start) {
        #if DEBUG__create_hHN_graph
        if (start == 0) {  // Print debug info for first vertex only to avoid spam
            Rprintf("\nProcessing start vertex %d:\n", start);

            // Print its adjacency list
            Rprintf("Adjacency list for vertex 0: ");
            for (size_t i = 0; i < adj_list[0].size(); ++i) {
                Rprintf("%d ", adj_list[0][i]);
            }
            Rprintf("\n");

            // Print its weights
            Rprintf("Weights for vertex 0: ");
            for (size_t i = 0; i < weight_list[0].size(); ++i) {
                Rprintf("%.4f ", weight_list[0][i]);
            }
            Rprintf("\n");
        }
        #endif

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
            #if DEBUG__create_hHN_graph
            if (start == 0) {
                Rprintf("Processing vertex %d with %d hops\n",
                        current, hops[current]);
            }
            #endif

            // Early termination if we've reached hop limit
            if (hops[current] > h) continue;
            #if DEBUG__create_hHN_graph
            if (hops[current] > h) {
                if (start == 0) Rprintf("Skipping vertex %d - exceeded hop limit\n", current);
                continue;
            }
            #endif

            // Process neighbors
            const auto& current_adj = adj_list[current];
            const auto& current_weights = weight_list[current];

            for (size_t i = 0; i < current_adj.size(); ++i) {
                const int neighbor = current_adj[i];

                if (neighbor == current) continue;  // Skip self-loops

                const double new_distance = distances[current] + current_weights[i];
                const int next_hops = hops[current] + 1;

                #if DEBUG__create_hHN_graph
                if (start == 0) {  // Debug for first vertex only
                    Rprintf("Checking neighbor %d: current distance %.4f, new distance %.4f, hops %d\n",
                            neighbor, distances[neighbor], new_distance, next_hops);
                }
                #endif

                if (new_distance < distances[neighbor]) {
                    #if DEBUG__create_hHN_graph
                    if (start == 0) {  // Debug for first vertex only
                        Rprintf("Updating neighbor %d: distance %.4f, hops %d\n",
                               neighbor, new_distance, next_hops);
                    }
                    #endif

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

        #if DEBUG__create_hHN_graph
        if (start == 0) {  // Debug for first vertex only
            Rprintf("\nFinal distances and hops for vertex 0:\n");
            for (int v = 0; v < std::min(5, n_vertices); v++) {  // Print first 5 vertices only
                Rprintf("Vertex %d: distance = %.4f, hops = %d\n",
                        v, distances[v], hops[v]);
            }
        }
        #endif

        // Collect neighbors within h hops
        for (int v = 0; v < n_vertices; ++v) {
            if (v != start && hops[v] <= h) {
                hhn_adj_list[start].push_back(v);
                hhn_weight_list[start].push_back(distances[v]);

                #if DEBUG__create_hHN_graph
                if (start == 0) {  // Debug for first vertex only
                    Rprintf("Added vertex %d to hHN graph with distance %.4f\n",
                            v, distances[v]);
                }
                #endif
            }
        }
    }

    #if DEBUG__create_hHN_graph
    // Add verification before return
    Rprintf("\nVerification before return:\n");
    Rprintf("hhn_adj_list[0].size(): %d\n", (int)hhn_adj_list[0].size());
    Rprintf("hhn_weight_list[0].size(): %d\n", (int)hhn_weight_list[0].size());

    auto result = std::make_pair(hhn_adj_list, hhn_weight_list);

    // Verify the result pair
    Rprintf("result.first[0].size(): %d\n", (int)result.first[0].size());
    Rprintf("result.second[0].size(): %d\n", (int)result.second[0].size());

    return result;
    #endif

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
 * @warning Ensure that s_k is a scalar integer and non-negative. The function currently does not perform
 *          extensive error checking on inputs.
 *
 * @see create_hHN_graph for the underlying C++ implementation.
 * @see Rgraph_to_vector and R_list_of_dvectors_to_cpp_vector_of_dvectors
 *      for the helper functions that convert R lists to C++ vectors.
 *
 * @example
 * In R:
 * result <- .Call("S_create_hHN_graph", adj_list, weight_list, h)
 * hhn_adj_list <- result$adj_list
 * hhn_dist_list <- result$dist_list
 */
SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h) {
    // Correct: Converts R inputs to C++ types
    std::vector<std::vector<int>> adj_vv = Rgraph_to_vector(s_adj_list);
    std::vector<std::vector<double>> weight_vv = Rweights_to_vector(s_weight_list);

    // Correct: Coerces s_k to integer and extracts its value
    PROTECT(s_h = coerceVector(s_h, INTSXP));
    int h = INTEGER(s_h)[0];
    UNPROTECT(1);

    #if 0
    PROTECT(s_type = coerceVector(s_type, INTSXP));
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
        error("Invalid type specified.");
    }
    #endif

    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> hhn_graph = create_hHN_graph(adj_vv, weight_vv, h);

    // Correct: Extracts results
    std::vector<std::vector<int>> hhn_adj_vv       = hhn_graph.first;
    std::vector<std::vector<double>> hhn_weight_vv = hhn_graph.second;

    // Correct: Prepares R list for adjacency list
    int n_vertices = static_cast<int>(hhn_adj_vv.size());
    SEXP adj_list = PROTECT(allocVector(VECSXP, n_vertices));
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, hhn_adj_vv[i].size()));
        int* A = INTEGER(RA);
        for (const auto& neighbor : hhn_adj_vv[i])
            *A++ = neighbor + 1;  // Correct: Adjusts for R's 1-based indexing
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Correct: Prepares R list for distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, n_vertices));
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, hhn_weight_vv[i].size()));
        double* D = REAL(RD);
        for (const auto& dist : hhn_weight_vv[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
        UNPROTECT(1);
    }

    // Correct: Prepares the final result list
    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, dist_list);

    UNPROTECT(2); // Unprotect adj_list and dist_list

    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("dist_list"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(2); // Unprotect res and names

    return res;
}
