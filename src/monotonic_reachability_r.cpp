#include "set_wgraph.hpp"                    // set_wgraph_t
#include "reachability_map.hpp"              // monotonic_reachability_map_t
#include "SEXP_cpp_conversion_utils.hpp"     // convert_adj_list_from_R, convert_weight_list_from_R
#include "error_utils.h"                     // if any errors are checked via REPORT_ERROR()

#include <vector>       // std::vector
#include <utility>      // std::pair
#include <algorithm>    // std::sort

#include <R.h>
#include <Rinternals.h>

// R interface function to test monotonic reachability map

extern "C" {
    SEXP S_test_monotonic_reachability_map(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_ref_vertex,
        SEXP s_radius,
        SEXP s_ascending
    );
}

/**
 * @brief R interface to test the monotonic reachability map functionality
 *
 * This function creates a graph from the provided adjacency and weight lists,
 * computes the monotonic reachability map, and returns the results as an R list.
 *
 * @param s_adj_list SEXP containing adjacency lists for the graph
 * @param s_weight_list SEXP containing weight lists for the graph edges
 * @param s_y SEXP containing function values at vertices
 * @param s_ref_vertex SEXP containing reference vertex index
 * @param s_radius SEXP containing search radius
 * @param s_ascending SEXP containing boolean for direction (TRUE for ascending)
 *
 * @return SEXP R list containing monotonic reachability information
 */
SEXP S_test_monotonic_reachability_map(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_ref_vertex,
    SEXP s_radius,
    SEXP s_ascending
) {
    // Convert input from R format
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    std::vector<double> y(REAL(s_y), REAL(s_y) + Rf_length(s_y));
    size_t ref_vertex = static_cast<size_t>(INTEGER(s_ref_vertex)[0] - 1); // Convert from 1-based to 0-based
    double radius = REAL(s_radius)[0];
    bool ascending = (LOGICAL(s_ascending)[0] == 1);
    
    // Create graph
    set_wgraph_t graph(adj_list, weight_list);
    
    // Compute monotonic reachability map
    monotonic_reachability_map_t map = graph.compute_monotonic_reachability_map(
        ref_vertex,
        y,
        radius,
        ascending
        );
    
    // Find best gradient vertex
    double min_distance = 0.0;
    size_t min_path_length = 2;
    auto best_vertex_info = graph.find_best_gradient_vertex(
        map,
        min_distance,
        min_path_length,
        ascending
        );
    
    // Prepare results for R
    // Return list structure:
    // - vertices: Vector of all reachable vertices (0-based in C++, converted to 1-based for R)
    // - monotonicity: Vector of monotonicity indices for each vertex
    // - total_change: Vector of total y changes for each vertex
    // - distances: Vector of distances from reference vertex
    // - paths: List of paths to each vertex (as vertex indices)
    // - best_vertex: Index of vertex with best monotonicity index (1-based)
    // - best_monotonicity: Monotonicity index of best vertex
    
    const char* names[] = {
        "vertices",
        "monotonicity",
        "total_change",
        "distances",
        "paths",
        "best_vertex",
        "best_monotonicity",
        NULL
    };
    
    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;
    
    // Create result list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    
    // Set names
    {
        SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
        for (int i = 0; i < n_elements; i++) {
            SET_STRING_ELT(r_result_names, i, Rf_mkChar(names[i]));
        }
        Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
        UNPROTECT(1);
    }

    // Prepare data vectors
    std::vector<size_t> vertices;
    std::vector<double> monotonicity;
    std::vector<double> total_change;
    std::vector<double> distances;
    
    // Sort vertices by monotonicity index for better visualization
    std::vector<std::pair<size_t, double>> vertex_mono_pairs;
    for (const auto& [vertex, info] : map.info) {
        if (vertex != ref_vertex) { // Skip reference vertex
            vertex_mono_pairs.push_back({vertex, info.monotonicity_index});
        }
    }
    
    // Sort by monotonicity index (descending)
    std::sort(vertex_mono_pairs.begin(), vertex_mono_pairs.end(),
              [](const auto& a, const auto& b) {
                  return a.second > b.second;
              });
    
    // Fill data vectors in sorted order
    for (const auto& [vertex, _] : vertex_mono_pairs) {
        vertices.push_back(vertex);
        monotonicity.push_back(map.info.at(vertex).monotonicity_index);
        total_change.push_back(map.info.at(vertex).total_change);
        distances.push_back(map.info.at(vertex).distance);
    }
    
    // Convert to R vectors
    size_t n_vertices = vertices.size();
    
    // Vertices vector (convert to 1-based)
    {
        SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_vertices));
        for (size_t i = 0; i < n_vertices; i++) {
            INTEGER(r_vertices)[i] = vertices[i] + 1; // Convert to 1-based
        }
        SET_VECTOR_ELT(r_result, 0, r_vertices);
        UNPROTECT(1);
    }

    // Monotonicity vector
    {
        SEXP r_monotonicity = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        for (size_t i = 0; i < n_vertices; i++) {
            REAL(r_monotonicity)[i] = monotonicity[i];
        }
        SET_VECTOR_ELT(r_result, 1, r_monotonicity);
        UNPROTECT(1);
    }

    // Total change vector
    {
        SEXP r_total_change = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        for (size_t i = 0; i < n_vertices; i++) {
            REAL(r_total_change)[i] = total_change[i];
        }
        SET_VECTOR_ELT(r_result, 2, r_total_change);
        UNPROTECT(1);
    }

    // Distances vector
    {
        SEXP r_distances = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        for (size_t i = 0; i < n_vertices; i++) {
            REAL(r_distances)[i] = distances[i];
        }
        SET_VECTOR_ELT(r_result, 3, r_distances);
        UNPROTECT(1);
    }

    // Paths list
    {
        SEXP r_paths = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (size_t i = 0; i < n_vertices; i++) {
            // Reconstruct path to this vertex
            path_t path = graph.reconstruct_monotonic_path(map, vertices[i]);

            // Convert path to R vector (1-based)
            SEXP r_path = PROTECT(Rf_allocVector(INTSXP, path.vertices.size()));
            for (size_t j = 0; j < path.vertices.size(); j++) {
                INTEGER(r_path)[j] = path.vertices[j] + 1; // Convert to 1-based
            }

            SET_VECTOR_ELT(r_paths, i, r_path);
            UNPROTECT(1); // r_path
        }
        SET_VECTOR_ELT(r_result, 4, r_paths);
        UNPROTECT(1);
    }

    // Best vertex
    {
        SEXP r_best_vertex = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_best_vertex)[0] = best_vertex_info.first + 1; // Convert to 1-based
        SET_VECTOR_ELT(r_result, 5, r_best_vertex);
        UNPROTECT(1);
    }

    // Best monotonicity
    {
        SEXP r_best_monotonicity = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(r_best_monotonicity)[0] = best_vertex_info.second;
        SET_VECTOR_ELT(r_result, 6, r_best_monotonicity);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return r_result;
}
