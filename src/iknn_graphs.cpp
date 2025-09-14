/*!
  Creating Intersetion Weighted kNN graph (IWkNN graph)

  Two points are connected by an edge in this graph if and only if they share at
  least one common point within their respective k-nearest neighbor sets. In
  other words, there's an intersection between their k-nearest neighborhoods.

  \textbf{Characteristics:}

  * \textbf{Reciprocity:} This connection is inherently reciprocal. If point A is within the k-nearest neighbors of point B, then point B will also be within the k-nearest neighbors of point A. This generally results in an undirected graph.
  * \textbf{Sensitivity to k:} The choice of the parameter k significantly influences the structure of this graph. Larger values of k tend to create denser graphs with more connections.
  * \textbf{Focus on Shared Neighbors:} This variant emphasizes points that are "locally similar," meaning they share commonalities within their local neighborhoods.

  \textbf{Potential Use Cases:}

  * \textbf{Clustering and Community Detection:} This type of graph could be valuable in identifying densely connected clusters of points that share neighbors, indicating groups exhibiting similar characteristics.
  * \textbf{Outlier Detection:}  Points with abnormally few connections in this graph might indicate potential outliers, as they don't consistently appear in the k-nearest neighborhoods of other points.
  * \textbf{Collaborative Filtering:} This construct might be applicable in recommendation systems where you want to emphasize items that have been "liked" by users with overlapping tastes.

  \textbf{Caveats:}

  * \textbf{Computational Cost:} Constructing this graph can be computationally more demanding than a standard k-NN graph because you need to compare k-nearest neighbor sets.
  * \textbf{Data Density:}  In very sparse datasets, this construct might lead to highly disconnected graphs.

  \textbf{Relationship with Other Graphs:}

  It's interesting to note that this graph would likely be denser than a Mutual k-Nearest Neighbor Graph (Mk-NN).  In an Mk-NN, two points need to be *each other's* k-nearest neighbors, while in this variant they just need to share at least one common neighbor in their sets.

*/

#include "omp_compat.h"
#include "iknn_graphs.hpp"
#include "set_wgraph.hpp"
#include "kNN_r.h"            // for S_kNN()
#include "kNN.h"              // for struct iknn_vertex_tt
#include "cpp_utils.hpp"      // for debugging
#include "progress_utils.hpp" // for elapsed.time
#include "SEXP_cpp_conversion_utils.hpp"
#include "edge_pruning_stats.hpp"

#include <vector>
#include <unordered_set>
#include <set>
#include <memory>
#include <utility>
#include <limits>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <queue>
#include <numeric>    // for std::iota

#include <ANN/ANN.h>  // ANN library header

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {

    SEXP S_verify_pruning(SEXP s_X,
                          SEXP s_k,
                          SEXP s_max_alt_path_length);

    SEXP S_create_single_iknn_graph(
        SEXP RX,
        SEXP Rk,
        SEXP s_pruning_thld,
        SEXP s_compute_full
        );

    SEXP S_create_iknn_graphs(
        SEXP s_X,
        SEXP s_kmin,
        SEXP s_kmax,
        // pruning parameters
        SEXP s_max_path_edge_ratio_thld,
        SEXP s_path_edge_ratio_percentile,
        // other
        SEXP s_compute_full,
        SEXP s_verbose
        );
}

iknn_graph_t create_iknn_graph(SEXP RX, SEXP Rk);

std::vector<int> union_find(const std::vector<std::vector<int>>& adj_vect);
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> prune_long_edges(const std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& wgraph, int version);
std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph);
std::vector<std::vector<std::pair<int, int>>> prune_edges_with_alt_paths(const std::vector<std::vector<std::pair<int, int>>>& graph,
                                                                         std::vector<int>& long_edge_isize,
                                                                         std::vector<int>& alt_path_lengths,
                                                                         std::vector<int>& alt_path_total_isize,
                                                                         int max_alt_path_length);


/**
 * Transforms an Intersection-Weighted-Distance k-Nearest Neighbor graph to an Intersection-Weighted k-Nearest Neighbor graph.
 *
 * This function converts a graph representation where each edge contains a distance measurement alongside the intersection size
 * to a representation that only considers the intersection size, discarding the distance information.
 *
 * @param graph The original graph represented as a vector of vectors of iknn_vertex_t, where each iknn_vertex_t contains
 *              the nearest neighbor index, intersection size, and distance.
 * @return A unique pointer to a vector of vectors of pairs, where each pair consists of a nearest neighbor index and an intersection size.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IWD_to_IW_kNN_graph(const std::vector<std::vector<iknn_vertex_t>>& graph) {
    auto iw_graph = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(graph.size());

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto& vertex : graph[i]) {
            (*iw_graph)[i].emplace_back(vertex.index, vertex.isize);
        }
    }

    return iw_graph;
}

// ------------------------------------------------------------------------------------------
//
// iknn_graph_t members
//
// ------------------------------------------------------------------------------------------

vect_wgraph_t iknn_graph_t::prune_graph(int max_alt_path_length) const {

    // Create result graph with same size as input
    vect_wgraph_t result;
    result.adjacency_list.resize(graph.size());

    // Collect all edges with their distances
    std::set<weighted_edge_t> edges;
    for (size_t i = 0; i < graph.size(); i++) {
        for (const auto& vertex : graph[i]) {
            if (i < vertex.index) {
                edges.insert(
                    {i,
                     vertex.index,
                     vertex.dist,
                     vertex.isize}
                    );
            }
        }
    }

    // Process each edge
    for (const auto& edge : edges) {
        auto alt_path_exists = find_alternative_path(edge.start,
                                                     edge.end,
                                                     edge.isize,
                                                     max_alt_path_length);

        // If no alternative path exists or it's too long, keep the edge
        if (!alt_path_exists) {
            result.adjacency_list[edge.start].push_back({edge.end, edge.dist});
            result.adjacency_list[edge.end].push_back({edge.start, edge.dist});
        }
    }

    return result;
}

bool iknn_graph_t::find_alternative_path(
    int start,
    int end,
    int edge_isize,
    int max_path_length) const {

    size_t n_vertices = graph.size();

    total_isize_cache.assign(n_vertices, -1); // Track total isize
    num_edges_cache.assign(n_vertices, 0);    // Track number of edges in path

    std::queue<int> q;

    // Initialize BFS from source
    q.push(start);  // Should be 'start'
    total_isize_cache[start] = 0;
    num_edges_cache[start] = 0;

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        // Add after popping from queue:
        if (max_path_length > 0 && num_edges_cache[current] >= max_path_length) {
            continue;
        }

        // Check all neighbors of current vertex
        for (const auto& neighbor_vertex : graph[current]) {
            int neighbor = neighbor_vertex.index;
            int neighbor_isize = neighbor_vertex.isize;

            // Skip if intersection size requirement is not met
            if (neighbor_isize <= edge_isize) {
                continue;
            }

            // If we haven't visited this neighbor yet
            if (total_isize_cache[neighbor] == -1) {
                total_isize_cache[neighbor] = total_isize_cache[current] + neighbor_isize;
                num_edges_cache[neighbor] = num_edges_cache[current] + 1;
                q.push(neighbor);

                // If we've reached the target, check if path length is valid
                if (neighbor == end && num_edges_cache[neighbor] >= 2) {
                    return true;
                }
            }
        }
    }
    return false;
}


/**
 * @brief Prints the graph structure to the R console
 *
 * This function prints a detailed representation of the graph including:
 * - Total number of vertices
 * - Total number of edges
 * - Average vertex degree
 * - Complete edge list with intersection sizes and distances
 *
 * For each edge, it prints a tuple (neighbor_index, intersection_size, distance)
 * where:
 * - neighbor_index is the shifted index of the neighboring vertex
 * - intersection_size is the size of the intersection of kNN sets
 * - distance is the minimum indirect distance between vertices
 *
 * @param vertex_index_shift An integer value added to all vertex indices in the output.
 *        Useful for converting between 0-based and 1-based indexing systems.
 * @param name Optional name to label the output (defaults to empty string)
 *
 * @note The function uses Rprintf for output, making it suitable for R package integration
 * @note Edge counts are divided by 2 as each edge is stored twice in the undirected graph
 *
 * @see iknn_vertex_t
 */
void iknn_graph_t::print(size_t vertex_index_shift,
                         const std::string& name) const {
    if (graph.empty()) {
        Rprintf("\nEmpty graph\n");
        return;
    }

    // Calculate total edges (existing code is correct but could be more readable)
    size_t total_edges = 0;
    for (const auto& vertex_neighbors : graph) {
        total_edges += vertex_neighbors.size();
    }
    total_edges /= 2;  // Each edge is counted twice

    if (!name.empty()) {
        Rprintf("\n\niknn_graph_t: '%s'\n", name.c_str());
    } else {
        Rprintf("\n");
    }

    // Print graph summary
    Rprintf("iknn graph structure:\n");
    Rprintf("Number of vertices: %zu\n", graph.size());
    Rprintf("Number of edges: %zu\n", total_edges);
    Rprintf("Average degree: %.2f\n", graph.empty() ? 0.0 : (2.0 * total_edges / graph.size()));

    // Print detailed edge list
    Rprintf("Edge List:\n");
    for (size_t vertex_id = 0; vertex_id < graph.size(); ++vertex_id) {
        Rprintf("Vertex %zu:\n", vertex_id + vertex_index_shift);
        for (const auto& neighbor : graph[vertex_id]) {
            Rprintf("\t(%zu, %zu, %.3f)\n",
                   neighbor.index + vertex_index_shift,
                   neighbor.isize,
                   neighbor.dist);
        }
    }
}

/**
 * @brief Prints the structure of an integer-indexed k-nearest neighbors graph
 *
 * @details This function prints a detailed representation of an undirected graph where each edge
 * is represented by a vertex index and an intersection size. The output includes:
 * - Graph name (if provided)
 * - Total number of vertices
 * - Total number of edges
 * - Average vertex degree
 * - Complete edge list with intersection sizes
 *
 * For each edge, it prints a tuple (neighbor_index, intersection_size) where:
 * - neighbor_index is the index of the neighboring vertex (adjusted by vertex_index_shift)
 * - intersection_size is the size of the intersection set between the vertex and its neighbor
 *
 * @param graph The graph structure represented as a vector of vectors of pairs,
 *        where each pair contains (vertex_index, intersection_size)
 * @param name Optional name identifier for the graph (defaults to empty string)
 * @param vertex_index_shift Integer value added to all vertex indices in the output
 *        (useful for converting between 0-based and 1-based indexing systems)
 *
 * @note The function uses Rprintf for output, making it suitable for R package integration
 * @note Edge counts are divided by 2 as each edge is stored twice in the undirected graph
 */
void print_iiknn_graph(
    const std::vector<std::vector<std::pair<int, int>>>& graph,
    const std::string& name = "",
    int vertex_index_shift = 0) {

    if (graph.empty()) {
        Rprintf("\nEmpty graph\n");
        return;
    }

    if (!name.empty()) {
        Rprintf("\n\n%s:\n", name.c_str());
    }

    // Calculate total edges
    size_t total_edges = 0;
    for (const auto& vertex_neighbors : graph) {
        total_edges += vertex_neighbors.size();
    }
    total_edges /= 2;  // Each edge is counted twice in undirected graph

    // Print graph summary
    Rprintf("iknn graph structure:\n");
    Rprintf("Number of vertices: %zu\n", graph.size());
    Rprintf("Number of edges: %zu\n", total_edges);
    Rprintf("Average degree: %.2f\n", graph.empty() ? 0.0 : (2.0 * total_edges / graph.size()));

    // Print edge list
    Rprintf("Edge List:\n");
    for (size_t vertex_id = 0; vertex_id < graph.size(); ++vertex_id) {
        Rprintf("Vertex %zu:\n", vertex_id + vertex_index_shift);
        for (const auto& [neighbor_index, isize] : graph[vertex_id]) {
            Rprintf("\t(%d, %d)\n", neighbor_index + vertex_index_shift, isize);
        }
    }
}


/**
 * @brief Verifies the equivalence between old and new graph pruning implementations
 *
 * @details This function compares the results of the old edge pruning implementation
 * (prune_edges_with_alt_paths) with the new implementation (iknn_graph.prune_graph).
 * It constructs kNN graphs using both methods and performs a detailed comparison of
 * the resulting graph structures, identifying any discrepancies in edges or weights.
 *
 * @param s_X SEXP containing a numeric matrix of input data points
 * @param s_k SEXP containing the number of nearest neighbors (k)
 * @param s_max_alt_path_length SEXP containing the maximum alternative path length for pruning
 *
 * @return SEXP A list containing three elements:
 *   - identical: logical, TRUE if no discrepancies found
 *   - total_discrepancies: integer, number of vertices with differences
 *   - discrepancies: list of length n_vertices, where each non-NULL element contains:
 *     - vertex index
 *     - missing edges (present in old but not in new implementation)
 *     - extra edges (present in new but not in old implementation)
 *     Each edge is represented as a pair (vertex_index, weight)
 *
 * @throws Rf_error if X cannot be coerced to a numeric matrix
 *
 * @note The function performs a comprehensive comparison including both topology
 * (edge existence) and edge weights
 *
 * @see create_iknn_graph
 * @see prune_edges_with_alt_paths
 * @see iknn_graph_t::prune_graph
 */
SEXP S_verify_pruning(SEXP s_X,
                      SEXP s_k,
                      SEXP s_max_alt_path_length) {
    int nprot = 0;
    PROTECT(s_X = Rf_coerceVector(s_X, REALSXP)); nprot++;
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }
    int *dimX = INTEGER(Rf_getAttrib(s_X, R_DimSymbol));
    size_t n_vertices = dimX[0];
    int max_alt_path_length = INTEGER(s_max_alt_path_length)[0];

    // Creating a kNN graph
    iknn_graph_t iknn_graph = create_iknn_graph(s_X, s_k);

    Rprintf("\n\nIn S_verify_pruning()\n");
    iknn_graph.print(0,"original iknn_graph");

    // Phase 1: Getting results from old implementation
    auto IW_graph = IWD_to_IW_kNN_graph(iknn_graph.graph);
    std::vector<int> long_edge_isize, alt_path_lengths, alt_path_total_isize;
    std::vector<std::vector<std::pair<int, int>>> old_pruned_graph = prune_edges_with_alt_paths(*IW_graph,
                                                                                                long_edge_isize,
                                                                                                alt_path_lengths,
                                                                                                alt_path_total_isize,
                                                                                                max_alt_path_length);
    print_iiknn_graph(old_pruned_graph, "Old pruned graph", 0);

    // Create vect_wgraph_t from old implementation results
    vect_wgraph_t old_pruned_vect_wgraph;
    old_pruned_vect_wgraph.adjacency_list.resize(n_vertices);
    for (size_t vertex = 0; vertex < n_vertices; vertex++) {
        for (auto neighbor_pair : old_pruned_graph[vertex]) {
            size_t neighbor = neighbor_pair.first;
            for (const auto& iknn_neighbor : iknn_graph.graph[vertex]) {
                if (iknn_neighbor.index == neighbor) {
                    old_pruned_vect_wgraph.adjacency_list[vertex].emplace_back(neighbor, iknn_neighbor.dist);
                }
            }
        }
    }

    old_pruned_vect_wgraph.print("old_pruned_vect_wgraph");

    // Phase 2: Get results from new implementation
    vect_wgraph_t new_pruned_vect_wgraph = iknn_graph.prune_graph(max_alt_path_length);

    new_pruned_vect_wgraph.print("new_pruned_vect_wgraph");


    // Comparison logic
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SEXP discrepancies = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
    SEXP vertex_names = PROTECT(Rf_allocVector(STRSXP, n_vertices)); nprot++;
    int total_discrepancies = 0;

    for (size_t vertex = 0; vertex < n_vertices; vertex++) {
        std::vector<edge_info_t>& old_edges = old_pruned_vect_wgraph.adjacency_list[vertex];
        std::vector<edge_info_t>& new_edges = new_pruned_vect_wgraph.adjacency_list[vertex];

        // Create sets of edges for easy comparison
        std::set<std::pair<int, double>> old_edge_set, new_edge_set;
        for (const auto& edge : old_edges) {
            old_edge_set.insert({edge.vertex, edge.weight});
        }
        for (const auto& edge : new_edges) {
            new_edge_set.insert({edge.vertex, edge.weight});
        }

        // Find differences
        std::vector<std::pair<int, double>> missing_in_new, extra_in_new;

        for (const auto& edge : old_edge_set) {
            if (new_edge_set.find(edge) == new_edge_set.end()) {
                missing_in_new.push_back(edge);
            }
        }

        for (const auto& edge : new_edge_set) {
            if (old_edge_set.find(edge) == old_edge_set.end()) {
                extra_in_new.push_back(edge);
            }
        }

        // If discrepancies found, create a report for this vertex
        if (!missing_in_new.empty() || !extra_in_new.empty()) {
            total_discrepancies++;

            SEXP vertex_report = PROTECT(Rf_allocVector(VECSXP, 3));

            // Missing edges
            SEXP missing = PROTECT(Rf_allocMatrix(REALSXP, missing_in_new.size(), 2));
            double* missing_ptr = REAL(missing);
            for (size_t i = 0; i < missing_in_new.size(); i++) {
                missing_ptr[i] = missing_in_new[i].first;
                missing_ptr[i + missing_in_new.size()] = missing_in_new[i].second;
            }

            // Extra edges
            SEXP extra = PROTECT(Rf_allocMatrix(REALSXP, extra_in_new.size(), 2));
            double* extra_ptr = REAL(extra);
            for (size_t i = 0; i < extra_in_new.size(); i++) {
                extra_ptr[i] = extra_in_new[i].first;
                extra_ptr[i + extra_in_new.size()] = extra_in_new[i].second;
            }

            SET_VECTOR_ELT(vertex_report, 0, Rf_ScalarInteger(vertex));
            SET_VECTOR_ELT(vertex_report, 1, missing);
            SET_VECTOR_ELT(vertex_report, 2, extra);

            SET_VECTOR_ELT(discrepancies, vertex, vertex_report);
            SET_STRING_ELT(vertex_names, vertex, Rf_mkChar(std::to_string(vertex).c_str()));

            UNPROTECT(3);  // vertex_report, missing, extra
        } else {
            SET_VECTOR_ELT(discrepancies, vertex, R_NilValue);
            SET_STRING_ELT(vertex_names, vertex, Rf_mkChar(std::to_string(vertex).c_str()));
        }
    }

    // Set names for the discrepancies list
    Rf_setAttrib(discrepancies, R_NamesSymbol, vertex_names);

    // Create final result
    SET_VECTOR_ELT(result, 0, Rf_ScalarLogical(total_discrepancies == 0));  // TRUE if no discrepancies
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(total_discrepancies));
    SET_VECTOR_ELT(result, 2, discrepancies);

    // Set names for the result list
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(result_names, 0, Rf_mkChar("identical"));
    SET_STRING_ELT(result_names, 1, Rf_mkChar("total_discrepancies"));
    SET_STRING_ELT(result_names, 2, Rf_mkChar("discrepancies"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    UNPROTECT(nprot);
    return result;
}



// ------------------------------------------------------------------------------------------
//
// print iknn graph
//
// ------------------------------------------------------------------------------------------

/**
 * @brief Prints the structure of an IkNN (Intersection k-Nearest Neighbors) graph with configurable detail levels
 *
 * @details This function displays the graph structure where each vertex maintains a list of its nearest neighbors
 * along with intersection information and distances. The output detail can be controlled via the level parameter.
 *
 * The function supports three levels of detail:
 * - Level 1: Prints only the neighbor indices (basic adjacency list)
 * - Level 2: Prints neighbor indices and distances
 * - Level 3: Prints complete information including neighbor indices, intersection indices, and distances
 *
 * @param graph A vector of vectors where each inner vector contains iknn_vertex_t structures representing
 *              the nearest neighbors of a vertex. Each iknn_vertex_t contains:
 *              - index: Index of the nearest neighbor
 *              - isize: Size of the intersection between N(x) and N(x_j)
 *              - dist: Minimum distance computed as min_k{d(x,x_k) + d(x_k,x_j)}
 * @param level Detail level for printing (default = 3):
 *              - 1: Print only index
 *              - 2: Print index and dist
 *              - 3: Print all components (index, isize, dist)
 *
 * @note Output is formatted using Rprintf for R package compatibility
 *
 * @example
 * auto g = std::make_unique<std::vector<std::vector<iknn_vertex_t>>>;
 * // Print full details
 * print_iknn_graph(*g);
 * // Print only adjacency list
 * print_iknn_graph(*g, 1);
 */
void print_iknn_graph(
    const std::vector<std::vector<iknn_vertex_t>>& graph,
    int level = 3
    ) {
    if (level < 1 || level > 3) {
        Rprintf("Error: level must be 1, 2, or 3\n");
        return;
    }

    if (graph.empty()) {
        Rprintf("Empty graph\n");
        return;
    }

    // Print header
    Rprintf("\nIkNN Graph Structure:\n");
    Rprintf("%s\n", std::string(80, '-').c_str());

    // For each vertex in the graph
    for (size_t vertex_idx = 0; vertex_idx < graph.size(); ++vertex_idx) {
        Rprintf("Vertex %zu neighbors:\n", vertex_idx);

        // If this vertex has no neighbors
        if (graph[vertex_idx].empty()) {
            Rprintf("  No neighbors\n");
            continue;
        }

        // Print column headers based on level
        if (level == 1) {
            Rprintf("%15s\n", "NN Index");
            Rprintf("%s\n", std::string(15, '-').c_str());
        } else if (level == 2) {
            Rprintf("%15s %15s\n", "NN Index", "Distance");
            Rprintf("%s\n", std::string(30, '-').c_str());
        } else { // level == 3
            Rprintf("%15s %15s %15s\n", "NN Index", "I Index", "Distance");
            Rprintf("%s\n", std::string(45, '-').c_str());
        }

        // Print each neighbor's information
        for (const auto& neighbor : graph[vertex_idx]) {
            if (level == 1) {
                Rprintf("%15zu\n", neighbor.index);
            } else if (level == 2) {
                std::ostringstream dist_stream;
                dist_stream << std::fixed << std::setprecision(6) << neighbor.dist;
                Rprintf("%15zu %15s\n",
                       neighbor.index,
                       dist_stream.str().c_str());
            } else { // level == 3
                std::ostringstream dist_stream;
                dist_stream << std::fixed << std::setprecision(6) << neighbor.dist;
                Rprintf("%15zu %15zu %15s\n",
                       neighbor.index,
                       neighbor.isize,
                       dist_stream.str().c_str());
            }
        }

        // Add spacing between vertices
        Rprintf("\n");
    }
}


// ------------------------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------------------------


/**
 * @brief Constructs an Intersection-Weighted-Distance k-Nearest Neighbors (kNN) graph.
 *
 * This function builds an Intersection-Weighted-Distance kNN graph from a given data matrix.
 * For each point, it finds its k nearest neighbors, then constructs edges between points
 * based on the number of common nearest neighbors and computes the distance as the
 * minimum sum of distances through common neighbors.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @return A unique pointer to a 2D vector of iknn_vertex_t structures representing the kNN graph.
 *
 * Each iknn_vertex_t structure contains:
 * - `index`: Index of the nearest neighbor.
 * - `isize`: Number of common elements in the kNN sets.
 * - `dist`: Distance metric computed as the minimum sum of distances through common neighbors.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 */
iknn_graph_t create_iknn_graph(SEXP RX, SEXP Rk) {

    #define DEBUG__create_iknn_graph 0

    PROTECT(RX = Rf_coerceVector(RX, REALSXP));
    int *dimX = INTEGER(Rf_getAttrib(RX, R_DimSymbol));
    UNPROTECT(1);
    size_t n_points = dimX[0];

    PROTECT(Rk = Rf_coerceVector(Rk, INTSXP));
    size_t k = INTEGER(Rk)[0];
    UNPROTECT(1);

    #if DEBUG__create_iknn_graph
    Rprintf("\nIn create_iknn_graph\tk: %d\tn_points: %d\n", k, n_points);
    #endif

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));
    UNPROTECT(1);

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);

    // Creating iknn_graph_t object directly with the right size
    iknn_graph_t res(n_points);

    // Perform k-NN search for each point
    size_t n_points_minus_one = n_points - 1;
    std::vector<int> intersection;
    for (size_t pt_i = 0; pt_i < n_points_minus_one; pt_i++) {
        // Copying indices of kNN of the pt_i point to nn_i
        for (size_t j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + n_points * j];
            sorted_nn_i[j] = nn_i[j];
        }
        std::sort(sorted_nn_i.begin(), sorted_nn_i.end()); // Ensure sorted for set intersection

        #if DEBUG__create_iknn_graph
        {
            //Rprintf("pt_i %zu: nn_i: ", pt_i);
            //for (int j = 0; j < k; j++) Rprintf("%d ", nn_i[j]);
            //Rprintf("\n");
            Rprintf("\n-----------------\n");
            Rprintf("pt_i %zu: sorted_nn_i: ", pt_i);
            for (int j = 0; j < k; j++) Rprintf("%d ", sorted_nn_i[j]);
            Rprintf("\n");
        }
        #endif

        for (size_t pt_j = pt_i + 1; pt_j < n_points; pt_j++) {
            // Copying indices of kNN of the pt_j point to nn_j
            for (size_t j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + n_points * j];
                sorted_nn_j[j] = nn_j[j];
            }
            std::sort(sorted_nn_j.begin(), sorted_nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(), sorted_nn_j.begin(), sorted_nn_j.end(), std::back_inserter(intersection));

            size_t common_count = intersection.size();

            #if DEBUG__create_iknn_graph
            {
                // Rprintf("\npt_j %zu: nn_j: ", pt_j);
                // for (int j = 0; j < k; j++) Rprintf("%d ", nn_j[j]);
                // Rprintf("\n");
                Rprintf("\npt_j %zu: sorted_nn_j: ", pt_j);
                for (int j = 0; j < k; j++) Rprintf("%d ", sorted_nn_j[j]);
                Rprintf("\n");
                Rprintf("common_count: %zu\n", common_count);
                Rprintf("intersection: ");
                for (int j = 0; j < common_count; j++) Rprintf("%d ", intersection[j]);
                Rprintf("\n");
            }
            #endif

            if (common_count > 0) {
                // Computing the minimum of d(x,x_k) + d(x_k,x_j)
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    size_t idx_i = std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin(); // The std::find function returns an iterator, not an index. We need to subtract the beginning iterator to get the index.
                    size_t idx_j = std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = distances[pt_i + n_points * idx_i];
                    double dist_j_k = distances[pt_j + n_points * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);

                    #if 0 // DEBUG__create_iknn_graph
                    Rprintf("pt_i %zu, pt_j %zu, x_k %d, idx_i %zu, idx_j %zu, dist_i_k %.6f, dist_j_k %.6f, min_dist %.6f\n",
                            pt_i, pt_j, x_k, idx_i, idx_j, dist_i_k, dist_j_k, min_dist);
                    #endif
                }

                #if 0 // DEBUG__create_iknn_graph
                Rprintf("pt_i %zu, pt_j %zu: final min_dist %.6f\n", pt_i, pt_j, min_dist);
                #endif

                // Add edge from pt_i to pt_j and from pt_j to pt_i
                res.graph[pt_i].emplace_back(iknn_vertex_t{pt_j, common_count, min_dist});
                res.graph[pt_j].emplace_back(iknn_vertex_t{pt_i, common_count, min_dist});
            }
        }
    }

    #if DEBUG__create_iknn_graph
    print_iknn_graph(res.graph, 3);
    #endif

    return res;
}

/**
 * @brief Computes and prunes an intersection-weighted k-nearest neighbors graph,
 *        with optional computation of full graph components.
 *
 * @details
 * This function performs the following operations:
 * 1. Constructs an intersection-weighted k-nearest neighbors (IWD-kNN) graph
 * 2. Computes intersection sizes and distances for each edge
 * 3. Prunes the graph using iknn_graph.prune_graph()
 * 4. If compute_full is TRUE:
 *    - Computes all graph components including original graph metrics
 * 5. Always computes:
 *    - Pruned graph adjacency and distance lists
 *    - Basic graph statistics including edge counts and reduction ratios
 *
 * The pruning process removes edges where:
 * - An alternative path exists with >= 2 edges
 * - All edges in the alternative path have larger intersection sizes
 *
 * @param s_X SEXP object representing the input data matrix. Must be a numeric matrix
 *           (not a data frame) where:
 *           - Rows represent data points (vertices)
 *           - Columns represent features
 *           - Values must be numeric (coercible to REALSXP)
 *
 * @param s_k SEXP object representing k, the number of nearest neighbors to consider
 *           for each point. Must be a positive integer.
 *
 * @param s_compute_full SEXP object (logical) controlling computation of optional components:
 *        - TRUE: Compute all graph components and metrics
 *        - FALSE: Compute only pruned graph components and statistics
 *
 * @return SEXP object (a named list) containing:
 *   Always computed components:
 *   - pruned_adj_list: Adjacency lists after edge pruning (1-based indices)
 *   - pruned_weight_list: Distances for edges in pruned graph
 *   - n_edges: Total number of edges in the original graph
 *   - n_edges_in_pruned_graph: Number of edges after pruning
 *   - n_removed_edges: Number of edges removed during pruning
 *   - edge_reduction_ratio: Proportion of edges removed (n_removed_edges / n_edges)
 *
 *   Components computed only when compute_full is TRUE:
 *   - adj_list: List of integer vectors containing adjacency lists (1-based indices)
 *   - isize_list: List of integer vectors with intersection sizes for each edge
 *   - weight_list: List of numeric vectors with distances for each edge
 *   - conn_comps: Integer vector identifying connected components
 *   - connected_components: Integer vector identifying connected components
 *
 * @throws Rf_error If:
 *   - Input matrix cannot be coerced to numeric
 *   - create_iknn_graph returns null
 *   - Pruned graph size doesn't match input size
 *
 * @pre
 *   - s_X must be a valid R matrix object
 *   - s_k must be a valid R integer
 *   - s_compute_full must be a valid R logical value
 *   - Input data should not contain NA/NaN values
 *
 * @post
 *   - All returned indices are 1-based (R convention)
 *   - Graph connectivity is preserved after pruning
 *   - Memory is properly freed (PROTECT/UNPROTECT balanced)
 *   - When compute_full is FALSE, optional components are set to NULL
 *
 * @note
 *   - Time complexity:
 *     * O(V * k * log(V)) for initial graph construction and pruning
 *     * Additional O(V * k) for full component computation when compute_full is TRUE
 *   - Space complexity:
 *     * O(V * k) in all cases (storage for pruned graph)
 *   - All graph operations maintain symmetry of edges
 *
 * @Rf_warning
 *   - Large k values relative to dataset size may impact performance
 *   - The function assumes dense numeric data; sparse matrices require preprocessing
 *   - When compute_full is FALSE, accessing optional components will return NULL
 *
 * @see
 *   - create_iknn_graph for initial graph construction
 *   - iknn_graph_t::prune_graph for graph pruning implementation
 */
extern "C" SEXP S_create_single_iknn_graph(SEXP s_X,
                                SEXP s_k,
                                SEXP s_pruning_thld,
                                SEXP s_compute_full) {
    int nprot = 0;
    PROTECT(s_X = Rf_coerceVector(s_X, REALSXP)); nprot++;
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(Rf_getAttrib(s_X, R_DimSymbol));
    double pruning_thld = REAL(s_pruning_thld)[0];
    int n_vertices = dimX[0];
    int compute_full = (LOGICAL(s_compute_full)[0] == 1);

    // Creating a kNN graph
    auto iknn_graph = create_iknn_graph(s_X, s_k);

    // Count total edges in original graph
    int n_edges = 0;
    for (const auto& vertex_edges : iknn_graph.graph) {
        n_edges += vertex_edges.size();
    }
    n_edges /= 2;

    // Create basic graph vectors needed for original graph
    auto adj_vect = std::vector<std::vector<int>>(n_vertices);
    auto dist_vect = std::vector<std::vector<double>>(n_vertices);
    for (size_t i = 0; i < iknn_graph.graph.size(); i++) {
        for (auto nn_vertex : iknn_graph.graph[i]) {
            adj_vect[i].push_back(nn_vertex.index);
            dist_vect[i].push_back(nn_vertex.dist);
        }
    }

    // Initialize all components as NULL
    SEXP adj_list = R_NilValue;
    SEXP intersection_size_list = R_NilValue;
    SEXP weight_list = R_NilValue;
    SEXP s_conn_comps = R_NilValue;
    SEXP pruned_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
    SEXP pruned_weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;

    #define USE_GEOMETRIC_PRUNING_IN_SINGLE_IKNN_GRAPH 1

    size_t n_edges_in_pruned_graph = 0;

    #if USE_GEOMETRIC_PRUNING_IN_SINGLE_IKNN_GRAPH
    // transfering iknn_graph_t to set_wgraph_t
    auto pruned_graph = set_wgraph_t(iknn_graph);

    // Compute the deviations using the optimized method
    auto rel_deviations = pruned_graph.compute_edge_weight_rel_deviations();
    //const double EPSILON = 1e-16; // Threshold for considering a value as zero

    // Track removed edges
    //int removed_count = 0;

    // Process edges
    for (size_t i = 0; i < rel_deviations.size(); i++) {

        if (rel_deviations[i].rel_deviation < pruning_thld) {

            size_t source = rel_deviations[i].source;
            size_t target = rel_deviations[i].target;
            // size_t intermediate = rel_deviations[i].best_intermediate;

            // Check if removing this edge would isolate any vertex
            if (pruned_graph.adjacency_list[source].size() <= 1) {
                continue;
                // REPORT_ERROR("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                //              source + 1, target + 1, intermediate + 1, source + 1);
            }

            if (pruned_graph.adjacency_list[target].size() <= 1) {
                continue;
                // REPORT_ERROR("Cannot remove edge (%zu,%zu) through intermediate %zu: Vertex %zu would become isolated with 0 neighbors",
                //              source + 1, target + 1, intermediate + 1, target + 1);
            }

            // Safe to remove the edge
            pruned_graph.remove_edge(source, target);
            //removed_count++;
        }
    }

    // Count edges in the pruned graph
    for (const auto& neighbors : pruned_graph.adjacency_list) {
        n_edges_in_pruned_graph += neighbors.size();
    }
    n_edges_in_pruned_graph /= 2;

    #else
    int max_alt_path_length = 2; //INTEGER(s_max_alt_path_length)[0];
    vect_wgraph_t pruned_graph = iknn_graph.prune_graph(max_alt_path_length);

    // Count edges in pruned graph
    int n_edges_in_pruned_graph = 0;
    for (const auto& vertex_edges : pruned_graph.adjacency_list) {
        n_edges_in_pruned_graph += vertex_edges.size();
    }
    n_edges_in_pruned_graph /= 2;

    #endif


    // Compute graph statistics
    SEXP stats = PROTECT(Rf_allocVector(REALSXP, 4)); nprot++;
    REAL(stats)[0] = n_edges;
    REAL(stats)[1] = n_edges_in_pruned_graph;
    REAL(stats)[2] = n_edges - n_edges_in_pruned_graph;
    REAL(stats)[3] = (double)(n_edges - n_edges_in_pruned_graph) / n_edges;

    // Fill pruned lists from pruned_graph
    for (int i = 0; i < n_vertices; i++) {
        const auto& edges = pruned_graph.adjacency_list[i];

        // Adjacency list
        SEXP RA = PROTECT(Rf_allocVector(INTSXP, edges.size()));
        int* A = INTEGER(RA);
        for (const auto& edge : edges) {
            *A++ = edge.vertex + 1;  // Convert to 1-based indexing for R
        }
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
        UNPROTECT(1);

        // Distance list
        SEXP RD = PROTECT(Rf_allocVector(REALSXP, edges.size()));
        double* D = REAL(RD);
        for (const auto& edge : edges) {
            *D++ = edge.weight;
        }
        SET_VECTOR_ELT(pruned_weight_list, i, RD);
        UNPROTECT(1);
    }

    if (compute_full) {
        // Create full components
        adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
        intersection_size_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
        weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;

        // Fill original graph components
        for (int i = 0; i < n_vertices; i++) {
            // Original graph adjacency list
            {
                SEXP RA = PROTECT(Rf_allocVector(INTSXP, adj_vect[i].size()));
                int* A = INTEGER(RA);
                for (auto neighbor : adj_vect[i])
                    *A++ = neighbor + 1;
                SET_VECTOR_ELT(adj_list, i, RA);
                UNPROTECT(1);
            }

            // Original graph intersection sizes
            {
                SEXP RW = PROTECT(Rf_allocVector(INTSXP, adj_vect[i].size()));
                int* W = INTEGER(RW);
                for (auto nn_vertex : iknn_graph.graph[i])
                    *W++ = nn_vertex.isize;
                SET_VECTOR_ELT(intersection_size_list, i, RW);
                UNPROTECT(1);
            }

            // Original graph distances
            {
                SEXP RD = PROTECT(Rf_allocVector(REALSXP, dist_vect[i].size()));
                double* D = REAL(RD);
                for (auto dist : dist_vect[i])
                    *D++ = dist;
                SET_VECTOR_ELT(weight_list, i, RD);
                UNPROTECT(1);
            }
        }

        // Compute connected components
        std::vector<int> conn_comps = union_find(adj_vect);
        s_conn_comps = PROTECT(Rf_allocVector(INTSXP, conn_comps.size())); nprot++;
        std::copy(conn_comps.begin(), conn_comps.end(), INTEGER(s_conn_comps));
    }

    // Prepare result list - note we removed some elements that were specific to the old pruning method
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 11)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, weight_list);
    SET_VECTOR_ELT(res, 3, s_conn_comps);
    SET_VECTOR_ELT(res, 4, pruned_adj_list);
    SET_VECTOR_ELT(res, 5, pruned_weight_list);
    SET_VECTOR_ELT(res, 6, Rf_ScalarReal(REAL(stats)[0]));
    SET_VECTOR_ELT(res, 7, Rf_ScalarReal(REAL(stats)[1]));
    SET_VECTOR_ELT(res, 8, Rf_ScalarReal(REAL(stats)[2]));
    SET_VECTOR_ELT(res, 9, Rf_ScalarReal(REAL(stats)[3]));
    SET_VECTOR_ELT(res, 10, s_conn_comps);  // Adding connected components to the end

    // Set names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 11)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("isize_list"));
    SET_STRING_ELT(names, 2, Rf_mkChar("weight_list"));
    SET_STRING_ELT(names, 3, Rf_mkChar("conn_comps"));
    SET_STRING_ELT(names, 4, Rf_mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 5, Rf_mkChar("pruned_weight_list"));
    SET_STRING_ELT(names, 6, Rf_mkChar("n_edges"));
    SET_STRING_ELT(names, 7, Rf_mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(names, 8, Rf_mkChar("n_removed_edges"));
    SET_STRING_ELT(names, 9, Rf_mkChar("edge_reduction_ratio"));
    SET_STRING_ELT(names, 10, Rf_mkChar("connected_components"));

    Rf_setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}


// -----------------------------------------------------------------------------------------------
//
//  Implementation of S_create_iknn_graphs() with edge birth/death time statistics
//
// -----------------------------------------------------------------------------------------------

struct edge_t {
    size_t start; ///< start index of an edge
    size_t end;   ///< end index of an edge; start < end
};

struct birth_death_time_t {
    size_t birth_time;
    size_t death_time;
};

// Hash function for edge_t
struct edge_hash {
    std::size_t operator()(const edge_t& e) const {
        return std::hash<size_t>()(e.start) ^ (std::hash<size_t>()(e.end) << 1);
    }
};

// Equality comparison for edge_t
bool operator==(const edge_t& lhs, const edge_t& rhs) {
    return lhs.start == rhs.start && lhs.end == rhs.end;
}

// Helper function to create R matrix from birth_death_map
SEXP convert_birth_death_map_to_matrix(
    const std::unordered_map<edge_t, birth_death_time_t, edge_hash>& birth_death_map) {

    // Create matrix
    SEXP result;
    PROTECT(result = Rf_allocMatrix(REALSXP, birth_death_map.size(), 4));
    double* data = REAL(result);

    size_t row = 0;
    for (const auto& entry : birth_death_map) {
        // Convert to 1-based indexing for R
        data[row] = entry.first.start + 1;
        data[row + birth_death_map.size()] = entry.first.end + 1;
        data[row + 2 * birth_death_map.size()] = entry.second.birth_time - 1; // note that kmin, kmax values set in R are incremented in C++ by 1, to account for the fact that ANN library includes the ref vertex in its set of kNN's
        data[row + 3 * birth_death_map.size()] = entry.second.death_time - 1;
        row++;
    }

    // Set column names properly
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);  // NULL for row names

    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("start"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("end"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("birth_time"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("death_time"));
    SET_VECTOR_ELT(dimnames, 1, colnames);

    Rf_setAttrib(result, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return result;
}

// Helper function to convert pruned graph to R representation
SEXP create_R_graph_representation(
    const std::vector<std::vector<std::pair<int, int>>>& pruned_graph,
    const std::vector<std::vector<iknn_vertex_t>>& original_graph) {

    int n_vertices = pruned_graph.size();
    int nprot = 0;

    // Create lists for adjacency, distances, and intersection sizes
    SEXP pruned_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
    SEXP pruned_weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;
    SEXP pruned_isize_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); nprot++;

    for (int i = 0; i < n_vertices; i++) {
        size_t n_neighbors = pruned_graph[i].size();

        // Create vectors for this vertex
        SEXP adj = PROTECT(Rf_allocVector(INTSXP, n_neighbors));
        SEXP dist = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        SEXP isize = PROTECT(Rf_allocVector(INTSXP, n_neighbors));

        // Fill vectors
        for (size_t j = 0; j < n_neighbors; j++) {
            INTEGER(adj)[j] = pruned_graph[i][j].first + 1;  // 1-based indexing
            INTEGER(isize)[j] = pruned_graph[i][j].second;

            // Find corresponding distance in original graph
            double distance = 0.0;
            for (const auto& vertex : original_graph[i]) {
                if (vertex.index == (size_t)pruned_graph[i][j].first) {
                    distance = vertex.dist;
                    break;
                }
            }
            REAL(dist)[j] = distance;
        }

        SET_VECTOR_ELT(pruned_adj_list, i, adj);
        SET_VECTOR_ELT(pruned_weight_list, i, dist);
        SET_VECTOR_ELT(pruned_isize_list, i, isize);

        UNPROTECT(3);
    }

    // Create return list
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(result, 0, pruned_adj_list);
    SET_VECTOR_ELT(result, 1, pruned_weight_list);
    SET_VECTOR_ELT(result, 2, pruned_isize_list);

    // Set names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("pruned_weight_list"));
    SET_STRING_ELT(names, 2, Rf_mkChar("pruned_isize_list"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return result;
}

// Helper function to convert graph to pruned form
std::vector<std::vector<std::pair<int, int>>> prune_graph(
    const std::vector<std::vector<iknn_vertex_t>>& graph,
    int max_alt_path_length,
    std::vector<int>& long_edge_isize,
    std::vector<int>& alt_path_lengths,
    std::vector<int>& alt_path_total_isize) {

    // Convert to IW_kNN_graph format
    auto IW_graph = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(graph.size());
    for (size_t i = 0; i < graph.size(); i++) {
        for (const auto& vertex : graph[i]) {
            (*IW_graph)[i].push_back({vertex.index, vertex.isize});
        }
    }

    // Prune the graph
    return prune_edges_with_alt_paths(*IW_graph,
                                      long_edge_isize,
                                      alt_path_lengths,
                                      alt_path_total_isize,
                                      max_alt_path_length);
}


// New struct to hold kNN results
struct knn_search_result_t {
    std::vector<std::vector<int>> indices;     // [n_points][k]
    std::vector<std::vector<double>> distances; // [n_points][k]
    size_t n_points;
    size_t k;

    knn_search_result_t(size_t n, size_t k_val) : n_points(n), k(k_val) {
        indices.resize(n_points, std::vector<int>(k));
        distances.resize(n_points, std::vector<double>(k));
    }
};

// Function to compute kNN once for max k
knn_search_result_t compute_knn(SEXP RX, int k) {
    PROTECT(RX = Rf_coerceVector(RX, REALSXP));
    int *dimX = INTEGER(Rf_getAttrib(RX, R_DimSymbol));
    size_t n_points = dimX[0];

    SEXP Rk = PROTECT(Rf_ScalarInteger(k));
    SEXP knn_res = PROTECT(S_kNN(RX, Rk));

    int *indices_raw = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances_raw = REAL(VECTOR_ELT(knn_res, 1));

    // Create and fill result structure
    knn_search_result_t result(n_points, k);

    // Reorganize data into more convenient format
    for (size_t i = 0; i < n_points; i++) {
        for (int j = 0; j < k; j++) {
            result.indices[i][j] = indices_raw[i + n_points * j];
            result.distances[i][j] = distances_raw[i + n_points * j];
        }
    }

    UNPROTECT(3);
    return result;
}

// Modified create_iknn_graph that uses pre-computed kNN results
iknn_graph_t create_iknn_graph(const knn_search_result_t& knn_results, int k) {
    size_t n_points = knn_results.n_points;
    iknn_graph_t res(n_points);

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);
    std::vector<int> intersection;

    for (size_t pt_i = 0; pt_i < n_points - 1; pt_i++) {
        // Get k nearest neighbors for point i
        for (int j = 0; j < k; j++) {
            nn_i[j] = knn_results.indices[pt_i][j];
            sorted_nn_i[j] = nn_i[j];
        }
        std::sort(sorted_nn_i.begin(), sorted_nn_i.end());

        for (size_t pt_j = pt_i + 1; pt_j < n_points; pt_j++) {
            // Get k nearest neighbors for point j
            for (int j = 0; j < k; j++) {
                nn_j[j] = knn_results.indices[pt_j][j];
                sorted_nn_j[j] = nn_j[j];
            }
            std::sort(sorted_nn_j.begin(), sorted_nn_j.end());

            intersection.clear();
            std::set_intersection(
                sorted_nn_i.begin(), sorted_nn_i.end(),
                sorted_nn_j.begin(), sorted_nn_j.end(),
                std::back_inserter(intersection)
            );

            size_t common_count = intersection.size();
            if (common_count > 0) {
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    size_t idx_i = std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin();
                    size_t idx_j = std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = knn_results.distances[pt_i][idx_i];
                    double dist_j_k = knn_results.distances[pt_j][idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                res.graph[pt_i].emplace_back(iknn_vertex_t{pt_j, common_count, min_dist});
                res.graph[pt_j].emplace_back(iknn_vertex_t{pt_i, common_count, min_dist});
            }
        }
    }

    return res;
}

/**
 * @brief Computes and analyzes a series of intersection-weighted k-nearest neighbors graphs with multiple pruning strategies
 *
 * @details
 * For a given range of k values [kmin, kmax], this function:
 * 1. Constructs and analyzes intersection-weighted k-nearest neighbors (IWD-kNN) graphs
 * 2. Applies two different pruning strategies to each graph:
 *    - Geometric pruning: Removes edges where alternative paths exist with better geometric properties
 *    - Intersection-size pruning: Removes edges where alternative paths with larger intersection sizes exist
 * 3. Computes comprehensive statistics for each graph and pruning method
 * 4. Optionally returns the full pruned graph structures for further analysis
 *
 * The geometric pruning uses the path-to-edge ratio to determine which edges to remove:
 * - Edges are pruned when an alternative path exists with a better path-to-edge ratio
 * - The pruning threshold and percentile can be configured via parameters
 *
 * The intersection-size pruning evaluates alternative paths with at least 2 edges and prunes
 * edges when all edges in the alternative path have larger intersection sizes.
 *
 * @param s_X SEXP object representing the input data matrix. Must be a numeric matrix
 *           where rows represent data points (vertices) and columns represent features
 *
 * @param s_kmin SEXP object (integer) representing the minimum k value to consider
 *              Must be positive
 *
 * @param s_kmax SEXP object (integer) representing the maximum k value to consider
 *              Must be greater than or equal to kmin
 *
 * @param s_max_path_edge_ratio_thld SEXP object (double) Maximum acceptable ratio of
 *        alternative path length to edge length for geometric pruning.
 *        Edges with ratio <= this value will be pruned. If <= 0, this pruning stage is skipped.
 *
 * @param s_path_edge_ratio_percentile SEXP object (double) Percentile threshold (0.0-1.0)
 *        for edge lengths to consider in geometric pruning. Only edges with length
 *        greater than this percentile are considered for pruning.
 *
 * @param s_compute_full SEXP object (logical) controlling computation of optional components:
 *                      - TRUE: Store complete pruned graph structures for both pruning methods
 *                      - FALSE: Return only statistics without full graph structures
 *
 * @param s_verbose SEXP object (logical) controlling progress reporting during computation
 *
 * @return SEXP object (a named list) containing:
 * - k_statistics: Matrix with columns:
 *   - k: k value
 *   - n_edges: Total edges in original graph
 *   - n_edges_in_pruned_graph: Edges remaining after geometric pruning
 *   - n_removed_edges: Edges removed during geometric pruning
 *   - edge_reduction_ratio: Proportion of edges removed by geometric pruning
 *   - n_edges_in_isize_pruned_graph: Edges remaining after intersection-size pruning
 *   - n_removed_edges_in_double_pruning: Additional edges removed by intersection-size pruning
 *   - double_edge_reduction_ratio: Proportion of edges removed by intersection-size pruning
 *
 * - geom_pruned_graphs: If compute_full is TRUE, list of geometrically pruned graphs for each k value,
 *                       each containing adj_list and weight_list components. NULL if compute_full is FALSE.
 *
 * - isize_pruned_graphs: If compute_full is TRUE, list of intersection-size pruned graphs for each k value,
 *                        each containing adj_list and weight_list components. NULL if compute_full is FALSE.
 *
 * - edge_pruning_stats: List of matrices, one for each k value, containing edge pruning statistics:
 *   - edge_length: Length of each analyzed edge
 *   - length_ratio: Ratio of alternative path length to edge length
 *
 * @note
 * - The function computes k-nearest neighbors only once with the maximum k value for efficiency
 * - Parallelization is implemented using OpenMP for processing multiple k values simultaneously
 * - The function uses intersection size and geometric distance for edge pruning decisions
 * - Edge pruning preserves graph connectivity while reducing redundant connections
 * - All indices in returned R objects are 1-based (R convention)
 *
 * @see
 * - create_iknn_graph(): Creates a single intersection-weighted k-nearest neighbors graph
 * - set_wgraph_t::prune_edges_geometrically(): Geometric pruning implementation
 * - iknn_graph_t::prune_graph(): Intersection-size pruning implementation
 */
SEXP S_create_iknn_graphs(
    SEXP s_X,
    SEXP s_kmin,
    SEXP s_kmax,
    // pruning parameters
    SEXP s_max_path_edge_ratio_thld,
    SEXP s_path_edge_ratio_percentile,
    // other
    SEXP s_compute_full,
    SEXP s_verbose) {

    auto total_start_time = std::chrono::steady_clock::now();

    int kmin = INTEGER(s_kmin)[0];
    int kmax = INTEGER(s_kmax)[0];

    double max_path_edge_ratio_thld   = REAL(s_max_path_edge_ratio_thld)[0];
    double path_edge_ratio_percentile = REAL(s_path_edge_ratio_percentile)[0];

    int compute_full = (LOGICAL(s_compute_full)[0] == 1);
    int verbose = (LOGICAL(s_verbose)[0] == 1);

    int nprot = 0;
    PROTECT(s_X = Rf_coerceVector(s_X, REALSXP)); nprot++;
    int* dimX = INTEGER(Rf_getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];

    if (verbose) {
        Rprintf("Processing k values from %d to %d for %d vertices\n", kmin, kmax, n_vertices);
    }

    // Create vector of k values
    std::vector<int> k_values(kmax - kmin + 1);
    std::iota(k_values.begin(), k_values.end(), kmin);

    // Parallel processing of graph creation and pruning
    if (verbose) Rprintf("Starting parallel graph processing\n");
    auto parallel_start_time = std::chrono::steady_clock::now();

    // Progress tracking
    std::atomic<int> progress_counter{0};
    // const size_t n_k_values = kmax - kmin + 1;

    // Compute kNN once for maximum k
    auto knn_results = compute_knn(s_X, kmax);

    // Create a vector to store edge pruning stats for each k value
    int max_alt_path_length = 2; //INTEGER(s_max_alt_path_length)[0];
    double threshold_percentile = 0.5;

    // Pre-allocate all data structures
    std::vector<set_wgraph_t> geom_pruned_graphs(kmax - kmin + 1);
    geom_pruned_graphs.resize(kmax - kmin + 1);

    std::vector<vect_wgraph_t> isize_pruned_graphs(kmax - kmin + 1);
    isize_pruned_graphs.resize(kmax - kmin + 1);

    // Prepare vectors to store results
    std::vector<std::vector<double>> k_statistics(kmax - kmin + 1);
    k_statistics.resize(kmax - kmin + 1);

    std::vector<edge_pruning_stats_t> all_edge_pruning_stats(kmax - kmin + 1);
    all_edge_pruning_stats.resize(kmax - kmin + 1);

// Set number of threads (optional)
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
    if (verbose) Rprintf("Using %d OpenMP threads\n", num_threads);
#endif

// Parallel region
    #ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
    #endif
    for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
        int k = kmin + k_idx;

        // Report progress - use critical section to avoid output interleaving
        if (verbose) {
            #ifdef _OPENMP
#pragma omp critical
            #endif
            {
                REprintf("\rProcessing k=%d (%d of %d) - %d%%",
                         k, k_idx+1, kmax-kmin+1,
                         static_cast<int>((100.0 * (k_idx+1)) / (kmax-kmin+1)));
            }
        }

        // Create local graph and process
        auto iknn_graph = create_iknn_graph(knn_results, k);

        // Count original edges
        size_t n_edges = 0;
        for (const auto& vertex_edges : iknn_graph.graph) {
            n_edges += vertex_edges.size();
        }
        n_edges /= 2;

        // Transfer iknn_graph to set_wgraph_t
        auto pruned_graph = set_wgraph_t(iknn_graph);

        if (max_path_edge_ratio_thld > 0) {
            pruned_graph = pruned_graph.prune_edges_geometrically(
                max_path_edge_ratio_thld,
                path_edge_ratio_percentile
                );
        }

        // Count edges in pruned graph
        size_t n_edges_in_pruned_graph = 0;
        for (const auto& neighbors : pruned_graph.adjacency_list) {
            n_edges_in_pruned_graph += neighbors.size();
        }
        n_edges_in_pruned_graph /= 2;


        // isize graph pruning
        vect_wgraph_t isize_pruned_graph = iknn_graph.prune_graph(max_alt_path_length);

        // Count edges in the pruned graph
        size_t n_edges_in_isize_pruned_graph = 0;
        for (const auto& neighbors : isize_pruned_graph.adjacency_list) {
            n_edges_in_isize_pruned_graph += neighbors.size();
        }
        n_edges_in_isize_pruned_graph /= 2;

        // Store results
        k_statistics[k_idx] = {
            (double)n_edges,
            (double)n_edges_in_pruned_graph,
            (double)(n_edges - n_edges_in_pruned_graph),
            (double)(n_edges - n_edges_in_pruned_graph) / n_edges,
            (double)n_edges_in_isize_pruned_graph,
            (double)(n_edges_in_pruned_graph - n_edges_in_isize_pruned_graph),
            (double)(n_edges_in_pruned_graph - n_edges_in_isize_pruned_graph) / n_edges_in_pruned_graph
        };

        // Store the computed graphs and stats
        isize_pruned_graphs[k_idx]    = std::move(isize_pruned_graph);
        all_edge_pruning_stats[k_idx] = pruned_graph.compute_edge_pruning_stats(threshold_percentile);
        geom_pruned_graphs[k_idx]          = std::move(pruned_graph);
    }

    if (verbose) {
        elapsed_time(parallel_start_time, "\rParallel processing completed", true);
    }

    // Set up a start time for creating return list objects
    auto serial_start_time = std::chrono::steady_clock::now();
    if (verbose) {
        Rprintf("Creating return list objects ... ");
    }

    //
    // Create R objects
    //

    // Create a list to hold all the edge pruning stats matrices
    SEXP edge_stats_list = PROTECT(Rf_allocVector(VECSXP, kmax - kmin + 1)); nprot++;

    // Fill the list with edge stats matrices for each k value
    for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
        const auto& edge_pruning_stats = all_edge_pruning_stats[k_idx];

        // Create matrix for this k value's stats
        SEXP edge_stats_matrix = PROTECT(Rf_allocMatrix(REALSXP, edge_pruning_stats.stats.size(), 2));
        double* stats_data = REAL(edge_stats_matrix);

        // Fill the matrix with data - edge_length and length_ratio
        for (size_t i = 0; i < edge_pruning_stats.stats.size(); i++) {
            const auto& stat = edge_pruning_stats.stats[i];
            stats_data[i] = stat.edge_length;
            stats_data[i + edge_pruning_stats.stats.size()] = stat.length_ratio;
        }

        // Set column names for the matrix
        SEXP stats_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(stats_dimnames, 0, R_NilValue);  // No row names

        SEXP stats_colnames = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(stats_colnames, 0, Rf_mkChar("edge_length"));
        SET_STRING_ELT(stats_colnames, 1, Rf_mkChar("length_ratio"));
        SET_VECTOR_ELT(stats_dimnames, 1, stats_colnames);

        Rf_setAttrib(edge_stats_matrix, R_DimNamesSymbol, stats_dimnames);

        // Add this matrix to the list
        SET_VECTOR_ELT(edge_stats_list, k_idx, edge_stats_matrix);

        UNPROTECT(3); // Unprotect the matrix, dimnames, and colnames
    }

    // Set names for the edge_stats_list (k values)
    SEXP edge_stats_list_names = PROTECT(Rf_allocVector(STRSXP, kmax - kmin + 1)); nprot++;
    for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
        SET_STRING_ELT(edge_stats_list_names, k_idx, Rf_mkChar(std::to_string(kmin + k_idx - 1).c_str()));
    }
    Rf_setAttrib(edge_stats_list, R_NamesSymbol, edge_stats_list_names);

    
    SEXP geom_pruned_graphs_list = R_NilValue;
    SEXP isize_pruned_graphs_list = R_NilValue;

    if (compute_full) {
        PROTECT(geom_pruned_graphs_list = Rf_allocVector(VECSXP, kmax - kmin + 1)); nprot++;
        PROTECT(isize_pruned_graphs_list = Rf_allocVector(VECSXP, kmax - kmin + 1)); nprot++;

        for (int k_idx = 0; k_idx < kmax - kmin + 1; k_idx++) {
            // Process first-level pruned graph
            const auto& pruned_graph = geom_pruned_graphs[k_idx];

            SEXP r_pruned_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
            SEXP r_pruned_weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));

            for (int i = 0; i < n_vertices; i++) {
                const auto& edges = pruned_graph.adjacency_list[i];

                {
                    SEXP RA = PROTECT(Rf_allocVector(INTSXP, edges.size()));
                    int* A = INTEGER(RA);
                    for (const auto& edge : edges) {
                        *A++ = edge.vertex + 1;
                    }
                    SET_VECTOR_ELT(r_pruned_adj_list, i, RA);
                    UNPROTECT(1);
                }

                {
                    SEXP RD = PROTECT(Rf_allocVector(REALSXP, edges.size()));
                    double* D = REAL(RD);
                    for (const auto& edge : edges) {
                        *D++ = edge.weight;
                    }
                    SET_VECTOR_ELT(r_pruned_weight_list, i, RD);
                    UNPROTECT(1);
                }
            }

            SEXP r_pruned_graph = PROTECT(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(r_pruned_graph, 0, r_pruned_adj_list);
            SET_VECTOR_ELT(r_pruned_graph, 1, r_pruned_weight_list);

            SEXP r_pruned_graph_names = PROTECT(Rf_allocVector(STRSXP, 2));
            SET_STRING_ELT(r_pruned_graph_names, 0, Rf_mkChar("adj_list"));
            SET_STRING_ELT(r_pruned_graph_names, 1, Rf_mkChar("weight_list"));
            Rf_setAttrib(r_pruned_graph, R_NamesSymbol, r_pruned_graph_names);

            SET_VECTOR_ELT(geom_pruned_graphs_list, k_idx, r_pruned_graph);
            UNPROTECT(4);

            // Process isize-pruned graph
            const auto& isize_pruned_graph = isize_pruned_graphs[k_idx];
            SEXP r_isize_pruned_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
            SEXP r_isize_pruned_weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));

            for (int i = 0; i < n_vertices; i++) {
                const auto& edges = isize_pruned_graph.adjacency_list[i];

                {
                    SEXP RA = PROTECT(Rf_allocVector(INTSXP, edges.size()));
                    int* A = INTEGER(RA);
                    for (const auto& edge : edges) {
                        *A++ = edge.vertex + 1;
                    }
                    SET_VECTOR_ELT(r_isize_pruned_adj_list, i, RA);
                    UNPROTECT(1);
                }

                {
                    SEXP RD = PROTECT(Rf_allocVector(REALSXP, edges.size()));
                    double* D = REAL(RD);
                    for (const auto& edge : edges) {
                        *D++ = edge.weight;
                    }
                    SET_VECTOR_ELT(r_isize_pruned_weight_list, i, RD);
                    UNPROTECT(1);
                }
            }

            SEXP r_isize_pruned_graph = PROTECT(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(r_isize_pruned_graph, 0, r_isize_pruned_adj_list);
            SET_VECTOR_ELT(r_isize_pruned_graph, 1, r_isize_pruned_weight_list);

            SEXP r_isize_pruned_graph_names = PROTECT(Rf_allocVector(STRSXP, 2));
            SET_STRING_ELT(r_isize_pruned_graph_names, 0, Rf_mkChar("adj_list"));
            SET_STRING_ELT(r_isize_pruned_graph_names, 1, Rf_mkChar("weight_list"));
            Rf_setAttrib(r_isize_pruned_graph, R_NamesSymbol, r_isize_pruned_graph_names);

            SET_VECTOR_ELT(isize_pruned_graphs_list, k_idx, r_isize_pruned_graph);
            UNPROTECT(4);
        }
    }

    // Create statistics matrix without frequency ratios that were computed in the disabled section
    SEXP k_stats_matrix = PROTECT(Rf_allocMatrix(REALSXP, k_statistics.size(), 8)); nprot++;
    double* data = REAL(k_stats_matrix);
    for (size_t i = 0; i < k_statistics.size(); i++) {
        data[i] = kmin + i;
        for (size_t j = 0; j < 7; j++) {
            data[i + (j + 1) * k_statistics.size()] = k_statistics[i][j];
        }
    }

    // Set column names
    SEXP k_stats_dimnames = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(k_stats_dimnames, 0, R_NilValue);

    SEXP k_stats_colnames = PROTECT(Rf_allocVector(STRSXP, 8)); nprot++;
    SET_STRING_ELT(k_stats_colnames, 0, Rf_mkChar("k"));
    SET_STRING_ELT(k_stats_colnames, 1, Rf_mkChar("n_edges"));
    SET_STRING_ELT(k_stats_colnames, 2, Rf_mkChar("n_edges_in_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 3, Rf_mkChar("n_removed_edges"));
    SET_STRING_ELT(k_stats_colnames, 4, Rf_mkChar("edge_reduction_ratio"));
    SET_STRING_ELT(k_stats_colnames, 5, Rf_mkChar("n_edges_in_isize_pruned_graph"));
    SET_STRING_ELT(k_stats_colnames, 6, Rf_mkChar("n_removed_edges_in_double_pruning"));
    SET_STRING_ELT(k_stats_colnames, 7, Rf_mkChar("double_edge_reduction_ratio"));

    SET_VECTOR_ELT(k_stats_dimnames, 1, k_stats_colnames);
    Rf_setAttrib(k_stats_matrix, R_DimNamesSymbol, k_stats_dimnames);

    // Create return list
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(result, 0, k_stats_matrix);
    SET_VECTOR_ELT(result, 1, geom_pruned_graphs_list);
    SET_VECTOR_ELT(result, 2, isize_pruned_graphs_list);
    SET_VECTOR_ELT(result, 3, edge_stats_list);

    // Set names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("k_statistics"));
    SET_STRING_ELT(names, 1, Rf_mkChar("geom_pruned_graphs"));
    SET_STRING_ELT(names, 2, Rf_mkChar("isize_pruned_graphs"));
    SET_STRING_ELT(names, 3, Rf_mkChar("edge_pruning_stats"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    if (verbose) {
        elapsed_time(serial_start_time, "DONE", true);
        elapsed_time(total_start_time, "Total elapsed time", true);
    }

    UNPROTECT(nprot);
    return result;
}
