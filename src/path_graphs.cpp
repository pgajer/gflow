#include "path_graphs.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include <vector>
#include <queue>
#include <limits>
#include <map>
#include <utility>
#include <mutex>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

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
 * @brief Convert a path_graph_t into an R list representation.
 *
 * @details Builds a named R list with adjacency, edge lengths, and shortest-path
 * information from @p path_graph. All vertex indices are converted from 0-based
 * (C++) to 1-based (R).
 *
 * @param path_graph [in] Reference to a path_graph_t containing:
 *   - adj_list:       std::vector<std::vector<int>>
 *   - weight_list:    std::vector<std::vector<double>>
 *   - shortest_paths: std::map<std::pair<int,int>, std::vector<int>>
 *
 * @return SEXP (PROTECTED) A named R list with components:
 *   - "adj_list":         list of integer vectors (1-based)
 *   - "edge_length_list": list of numeric vectors
 *   - "shortest_paths":   list with:
 *       * "i":    integer vector of source vertices (1-based)
 *       * "j":    integer vector of target vertices (1-based)
 *       * "paths": list of integer vectors (each path, 1-based)
 *
 * @note The returned SEXP is in a PROTECTED state. Callers must UNPROTECT(1)
 *       after anchoring it (e.g., via SET_VECTOR_ELT) or otherwise making it safe.
 * @error On allocation failure or inconsistent inputs, calls Rf_error() and does not return.
 */
SEXP path_graph_from_path_graph_t(path_graph_t& path_graph) {

    std::vector<std::vector<int>>    path_graph_adj_vect                     = path_graph.adj_list;
    std::vector<std::vector<double>> path_graph_weight_vect                  = path_graph.weight_list;
    std::map<std::pair<int,int>, std::vector<int>> path_graph_shortest_paths = path_graph.shortest_paths;

    // Creating the final result list
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));

    // res slot names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
        SET_STRING_ELT(names, 1, Rf_mkChar("edge_length_list"));
        SET_STRING_ELT(names, 2, Rf_mkChar("shortest_paths"));
        Rf_setAttrib(res, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    int n_vertices = static_cast<int>(path_graph_adj_vect.size());

    // Creating an R list for the path graph adjacency list
    {
        SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RA = PROTECT(Rf_allocVector(INTSXP, path_graph_adj_vect[i].size()));
            int* A = INTEGER(RA);
            for (const auto& neighbor : path_graph_adj_vect[i])
                *A++ = neighbor + 1; // Converting to 1-based indexing
            SET_VECTOR_ELT(adj_list, i, RA);
            UNPROTECT(1);
        }

        // Creating an R list for the path graph edge length list
        SEXP edge_length_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RD = PROTECT(Rf_allocVector(REALSXP, path_graph_weight_vect[i].size()));
            double* D = REAL(RD);
            for (const auto& dist : path_graph_weight_vect[i])
                *D++ = dist;
            SET_VECTOR_ELT(edge_length_list, i, RD);
            UNPROTECT(1);
        }

        SET_VECTOR_ELT(res, 0, adj_list);
        SET_VECTOR_ELT(res, 1, edge_length_list);
        UNPROTECT(2); // adj_list, edge_length_list
    }

    // Creating an R list for the path graph shortest paths map
    {
        const size_t n_pairs = path_graph_shortest_paths.size();
        SEXP shortest_paths = PROTECT(Rf_allocVector(VECSXP, 3));

        // First component: i coordinates
        SEXP i_coords = PROTECT(Rf_allocVector(INTSXP, n_pairs));
        int* i_ptr = INTEGER(i_coords);

        // Second component: j coordinates
        SEXP j_coords = PROTECT(Rf_allocVector(INTSXP, n_pairs));
        int* j_ptr = INTEGER(j_coords);

        // Third component: paths list
        SEXP paths = PROTECT(Rf_allocVector(VECSXP, n_pairs));

        size_t idx = 0;
        for (const auto& kv : path_graph_shortest_paths) {
            const int i = kv.first.first;
            const int j = kv.first.second;
            const std::vector<int>& path = kv.second;

            i_ptr[idx] = i + 1; // 1-based for R
            j_ptr[idx] = j + 1;

            SEXP path_vec = PROTECT(Rf_allocVector(INTSXP, path.size()));
            int* path_ptr = INTEGER(path_vec);
            for (size_t k = 0; k < path.size(); ++k)
                path_ptr[k] = path[k] + 1;

            SET_VECTOR_ELT(paths, idx, path_vec);
            UNPROTECT(1); // path_vec
            ++idx;
        }

        // Set components of shortest_paths list
        SET_VECTOR_ELT(shortest_paths, 0, i_coords);
        SET_VECTOR_ELT(shortest_paths, 1, j_coords);
        SET_VECTOR_ELT(shortest_paths, 2, paths);
        UNPROTECT(3); // i_coords, j_coords, paths

        // Names for shortest_paths components
        {
            SEXP sp_names = PROTECT(Rf_allocVector(STRSXP, 3));
            SET_STRING_ELT(sp_names, 0, Rf_mkChar("i"));
            SET_STRING_ELT(sp_names, 1, Rf_mkChar("j"));
            SET_STRING_ELT(sp_names, 2, Rf_mkChar("paths"));
            Rf_setAttrib(shortest_paths, R_NamesSymbol, sp_names);
            UNPROTECT(1); // sp_names
        }

        SET_VECTOR_ELT(res, 2, shortest_paths);
        UNPROTECT(1); // shortest_paths
    }

    return res;
}


/**
 * @brief Creates a series of path graphs for different hop limits.
 *
 * @details Computes multiple path graphs from the same input graph for each hop
 * limit in @p s_h_values. Implementations may reuse work up to the maximum hop
 * limit to avoid redundant computation.
 *
 * @param s_adj_list    [in] R list; element i is an integer vector of neighbors of vertex i (1-based).
 * @param s_weight_list [in] R list; element i is a numeric vector of edge weights aligned with s_adj_list[[i]].
 * @param s_h_values    [in] Integer vector of hop limits for which to compute path graphs (will be coerced to INTSXP).
 *
 * @return SEXP An R list of length length(@p s_h_values). Each element is a named list with:
 *   - "adj_list":         list of integer vectors (reachable vertices, 1-based)
 *   - "edge_length_list": list of numeric vectors (path lengths)
 *   - "shortest_paths":   list with components "i", "j", and "paths"
 *
 * @note All vertex indices in the result use R’s 1-based convention.
 * @error On invalid inputs, calls Rf_error() and does not return.
 *
 * @see create_path_graph_series
 * @see path_graph_from_path_graph_t
 */
SEXP S_create_path_graph_series(SEXP s_adj_list,
                                SEXP s_weight_list,
                                SEXP s_h_values) {
    // Convert lists (no R allocations here)
    std::vector<std::vector<int>>    adj_vect    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);

    // Require integer input (the R caller already uses as.integer())
    if (!Rf_isInteger(s_h_values)) {
        Rf_error("h_values must be an integer vector (use as.integer()).");
    }

    const int n = LENGTH(s_h_values);           // LENGTH-first
    if ((size_t)n > (size_t)INT_MAX)            // overflow guard per your policy
        Rf_error("too large");

    const int* hv = INTEGER(s_h_values);        // safe after type check
    std::vector<int> h_values(hv, hv + n);      // no PROTECTs needed

    std::vector<path_graph_t> graph_series = create_path_graph_series(adj_vect, weight_vect, h_values);

    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n));
    for (int i = 0; i < n; ++i) {
        SEXP path_graph = PROTECT(path_graph_from_path_graph_t(graph_series[i]));
        SET_VECTOR_ELT(r_result, i, path_graph);
        UNPROTECT(1); // path_graph
    }

    UNPROTECT(1); // r_result
    return r_result;
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
 * @note Returns PROTECTed object !!! user has to UNPROTECT it!!!
 *
 */
SEXP S_create_path_graph_plus(SEXP s_adj_list,
                              SEXP s_edge_length_list,
                              SEXP s_h) {

    // ---- Convert inputs (no R allocations here) ----
    std::vector<std::vector<int>>    adj_vect    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_edge_length_list);

    // LENGTH-first scalar read for h
    if (!Rf_isInteger(s_h) || LENGTH(s_h) < 1)
        Rf_error("h must be an integer scalar.");
    const int h = INTEGER(s_h)[0];

    // ---- Build structure ----
    path_graph_plus_t path_graph = create_path_graph_plus(adj_vect, weight_vect, h);

    // Extract components (C++ side)
    std::vector<std::vector<int>>           path_graph_adj_vect    = path_graph.adj_list;
    std::vector<std::vector<double>>        path_graph_weight_vect = path_graph.weight_list;
    std::vector<std::vector<int>>           path_graph_hop_vect    = path_graph.hop_list;
    std::map<std::pair<int,int>, std::vector<int>> shortest_paths  = path_graph.shortest_paths;

    // LENGTH-first + overflow guard for vertex count
    const size_t nv_sz = path_graph_adj_vect.size();
    if (nv_sz > (size_t)INT_MAX) Rf_error("too large");
    const int n_vertices = (int) nv_sz;

    // ---- Parent result list [adj_list, edge_length_list, hop_list, shortest_paths] ----
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 4));
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
        SET_STRING_ELT(names, 1, Rf_mkChar("edge_length_list"));
        SET_STRING_ELT(names, 2, Rf_mkChar("hop_list"));
        SET_STRING_ELT(names, 3, Rf_mkChar("shortest_paths"));
        Rf_setAttrib(r_result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    // ---- 0) adj_list: list of integer vectors (1-based neighbors) ----
    {
        SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; ++i) {
            const size_t deg_sz = path_graph_adj_vect[(size_t)i].size();
            if (deg_sz > (size_t)INT_MAX) { UNPROTECT(1); Rf_error("too large"); }
            SEXP RA = PROTECT(Rf_allocVector(INTSXP, (int)deg_sz));
            int* A = INTEGER(RA);
            for (int k = 0, dk = (int)deg_sz; k < dk; ++k) {
                A[k] = (int) path_graph_adj_vect[(size_t)i][(size_t)k] + 1; // 1-based
            }
            SET_VECTOR_ELT(adj_list, i, RA);
            UNPROTECT(1); // RA
        }
        SET_VECTOR_ELT(r_result, 0, adj_list);
        UNPROTECT(1); // adj_list
    }

    // ---- 1) edge_length_list: list of numeric vectors (weights) ----
    {
        SEXP edge_length_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; ++i) {
            const size_t deg_sz = path_graph_weight_vect[(size_t)i].size();
            if (deg_sz > (size_t)INT_MAX) { UNPROTECT(1); Rf_error("too large"); }
            SEXP RD = PROTECT(Rf_allocVector(REALSXP, (int)deg_sz));
            double* D = REAL(RD);
            if (deg_sz > 0) {
                const auto& src = path_graph_weight_vect[(size_t)i];
                std::copy(src.begin(), src.end(), D);
            }
            SET_VECTOR_ELT(edge_length_list, i, RD);
            UNPROTECT(1); // RD
        }
        SET_VECTOR_ELT(r_result, 1, edge_length_list);
        UNPROTECT(1); // edge_length_list
    }

    // ---- 2) hop_list: list of integer vectors (hop counts) ----
    {
        SEXP hop_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; ++i) {
            const size_t hop_sz = path_graph_hop_vect[(size_t)i].size();
            if (hop_sz > (size_t)INT_MAX) { UNPROTECT(1); Rf_error("too large"); }
            SEXP RA = PROTECT(Rf_allocVector(INTSXP, (int)hop_sz));
            int* A = INTEGER(RA);
            for (int k = 0, dk = (int)hop_sz; k < dk; ++k) {
                A[k] = (int) path_graph_hop_vect[(size_t)i][(size_t)k]; // already 1-based per spec
            }
            SET_VECTOR_ELT(hop_list, i, RA);
            UNPROTECT(1); // RA
        }
        SET_VECTOR_ELT(r_result, 2, hop_list);
        UNPROTECT(1); // hop_list
    }

    // ---- 3) shortest_paths: list with {i, j, paths} ----
    {
        // Convert map -> vector to have a contiguous container and count
        std::vector<std::pair<std::pair<int,int>, std::vector<int>>> path_pairs;
        path_pairs.reserve(shortest_paths.size());
        for (const auto& p : shortest_paths) path_pairs.push_back(p);

        if (path_pairs.size() > (size_t)INT_MAX) { UNPROTECT(1); Rf_error("too large"); }
        const int npairs = (int) path_pairs.size();

        SEXP shortest_paths_r = PROTECT(Rf_allocVector(VECSXP, 3));
        {
            SEXP sp_names = PROTECT(Rf_allocVector(STRSXP, 3));
            SET_STRING_ELT(sp_names, 0, Rf_mkChar("i"));
            SET_STRING_ELT(sp_names, 1, Rf_mkChar("j"));
            SET_STRING_ELT(sp_names, 2, Rf_mkChar("paths"));
            Rf_setAttrib(shortest_paths_r, R_NamesSymbol, sp_names);
            UNPROTECT(1); // sp_names
        }

        // Components
        SEXP i_coords = PROTECT(Rf_allocVector(INTSXP, npairs));
        SEXP j_coords = PROTECT(Rf_allocVector(INTSXP, npairs));
        SEXP paths    = PROTECT(Rf_allocVector(VECSXP,  npairs));

        int* ip = INTEGER(i_coords);
        int* jp = INTEGER(j_coords);

        for (int idx = 0; idx < npairs; ++idx) {
            const auto& ij  = path_pairs[(size_t)idx].first;
            const auto& path_vec_src = path_pairs[(size_t)idx].second;

            ip[idx] = (int) ij.first  + 1; // 1-based
            jp[idx] = (int) ij.second + 1; // 1-based

            if (path_vec_src.size() > (size_t)INT_MAX) { UNPROTECT(3); UNPROTECT(1); Rf_error("too large"); }
            const int npath = (int) path_vec_src.size();

            SEXP path_vec = PROTECT(Rf_allocVector(INTSXP, npath));
            int* pp = INTEGER(path_vec);
            for (int k = 0; k < npath; ++k) pp[k] = (int) path_vec_src[(size_t)k] + 1; // 1-based

            SET_VECTOR_ELT(paths, idx, path_vec);
            UNPROTECT(1); // path_vec
        }

        SET_VECTOR_ELT(shortest_paths_r, 0, i_coords);
        SET_VECTOR_ELT(shortest_paths_r, 1, j_coords);
        SET_VECTOR_ELT(shortest_paths_r, 2, paths);
        UNPROTECT(3); // i_coords, j_coords, paths

        SET_VECTOR_ELT(r_result, 3, shortest_paths_r);
        UNPROTECT(1); // shortest_paths_r
    }

    UNPROTECT(1); // r_result
    return r_result;
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

    const int n_vertices = adj_list.size();

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

    std::mutex shortest_paths_mutex;
    
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
            if (hops[current] >= h) {
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
                    std::lock_guard<std::mutex> lock(shortest_paths_mutex);
                    auto& path = result.shortest_paths[{start, v}];
                    path.clear();
                    path.reserve(hops[v] + 1);  // Exact size known

                    // Reconstruct path
                    for (int current = v; current != -1; current = parent[current]) {
                        path.push_back(current);
                    }

                    std::reverse(path.begin(), path.end());

                    if (hops[v] == h)
                        result.longest_paths.push_back({start, v});
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
    // Convert inputs (no R allocations here)
    std::vector<std::vector<int>>    adj_vect    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_edge_length_list);

    // LENGTH-first for scalar h
    if (!Rf_isInteger(s_h) || LENGTH(s_h) < 1)
        Rf_error("h must be an integer scalar.");
    const int h = INTEGER(s_h)[0];

    // Build structure
    path_graph_plm_t path_graph = create_path_graph_plm(adj_vect, weight_vect, h);

    // Parent list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 5));                 // [P1]

    // ---- Base components (must PROTECT the callee's return) ----
    {
        SEXP base_components = PROTECT(S_create_path_graph_plus(        // [P2]
            s_adj_list, s_edge_length_list, s_h));

        // Copy first 4 elements into parent (no extra PROTECT needed for children)
        for (int i = 0; i < 4; ++i) {
            SET_VECTOR_ELT(r_result, i, VECTOR_ELT(base_components, i));
        }

        // Names: copy base names + add "vertex_paths"
        {
            SEXP names       = PROTECT(Rf_allocVector(STRSXP, 5));      // [P3]
            SEXP base_names  = PROTECT(Rf_getAttrib(base_components,    // [P4]
                                                    R_NamesSymbol));
            for (int i = 0; i < 4; ++i)
                SET_STRING_ELT(names, i, STRING_ELT(base_names, i));
            SET_STRING_ELT(names, 4, Rf_mkChar("vertex_paths"));

            Rf_setAttrib(r_result, R_NamesSymbol, names);
            UNPROTECT(2); // base_names [P4], names [P3]
        }

        UNPROTECT(1); // base_components [P2]
    }

    // ---- vertex_paths (slot 4) ----
    {
        // LENGTH-first + overflow guard
        const size_t nv_sz = path_graph.vertex_paths.size();
        if (nv_sz > (size_t)INT_MAX) { UNPROTECT(1); Rf_error("too large"); }
        const int n_vertices = (int) nv_sz;

        SEXP vertex_paths_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); // [P5]

        for (int v = 0; v < n_vertices; ++v) {
            const auto& vp = path_graph.vertex_paths[(size_t)v];
            const size_t np_sz = vp.containing_paths.size();
            if (np_sz > (size_t)INT_MAX) { UNPROTECT(2); Rf_error("too large"); }
            const int n_paths = (int) np_sz;

            SEXP path_matrix = PROTECT(Rf_allocMatrix(INTSXP, n_paths, 3));    // [P6]

            // Column names
            {
                SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 3));            // [P7]
                SET_STRING_ELT(colnames, 0, Rf_mkChar("start"));
                SET_STRING_ELT(colnames, 1, Rf_mkChar("end"));
                SET_STRING_ELT(colnames, 2, Rf_mkChar("position"));

                SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));            // [P8]
                SET_VECTOR_ELT(dimnames, 0, R_NilValue);
                SET_VECTOR_ELT(dimnames, 1, colnames);
                Rf_setAttrib(path_matrix, R_DimNamesSymbol, dimnames);
                UNPROTECT(2); // dimnames [P8], colnames [P7]
            }

            // Fill columns (1-based indices for R)
            int* m = INTEGER(path_matrix);
            for (int i = 0; i < n_paths; ++i) {
                m[i]               = (int) vp.containing_paths[(size_t)i].first  + 1;
                m[i + n_paths]     = (int) vp.containing_paths[(size_t)i].second + 1;
                m[i + 2*n_paths]   = (int) vp.position_in_path[(size_t)i]        + 1;
            }

            SET_VECTOR_ELT(vertex_paths_list, v, path_matrix);
            UNPROTECT(1); // path_matrix [P6]
        }

        SET_VECTOR_ELT(r_result, 4, vertex_paths_list);
        UNPROTECT(1); // vertex_paths_list [P5]
    }

    UNPROTECT(1); // r_result [P1]
    return r_result;
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
    if (!Rf_isNewList(s_path_graph)) {
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
        int n_paths = Rf_nrows(path_matrix);

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


/**
 * @brief Validates that each vertex in the path graph has at least one path of length h containing it
 *
 * @param path_graph The path graph structure to validate
 *
 * @return bool True if each vertex has at least one path of length h containing it, false otherwise
 *
 * @throws None
 *
 * @pre path_graph must be a valid path_graph_plm_t structure with h >= 0
 *
 * @note Time Complexity: O(V * P) where V is number of vertices and P is average number of paths per vertex
 *
 * Example usage:
 * @code
 *     path_graph_plm_t pg = create_path_graph_plm(adj_list, weight_list, h);
 *     if (!validate_vertex_paths(pg)) {
 *         Rf_error("Path validation failed");
 *     }
 * @endcode
 *
 * @see path_graph_plm_t
 * @see create_path_graph_plm
 */
bool validate_vertex_paths(const path_graph_plm_t& path_graph) {
    const int path_n_vertices = path_graph.h + 1;
    const int n_vertices = path_graph.adj_list.size();

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        auto vertex_paths = path_graph.get_paths_containing(vertex, path_n_vertices);
        if (vertex_paths.empty()) {
            return false;
        }
    }
    return true;
}
