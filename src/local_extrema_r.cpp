#include <R.h>
#include <Rinternals.h>
#undef length

#include <vector>

#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "error_utils.h"

extern "C" {
    SEXP S_detect_local_extrema(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_max_radius,
        SEXP s_min_neighborhood_size,
        SEXP s_detect_maxima
    );
}

/**
 * @brief R interface function for detecting local extrema
 *
 * This function creates a graph, runs the local extrema detection algorithm,
 * and returns the results as an R list.
 *
 * @param s_adj_list SEXP containing adjacency lists for the graph
 * @param s_weight_list SEXP containing weight lists for the graph edges
 * @param s_y SEXP containing function values at vertices
 * @param s_max_radius SEXP containing maximum radius for neighborhood search
 * @param s_min_neighborhood_size SEXP containing minimum required neighborhood size
 * @param s_detect_maxima SEXP boolean indicating whether to detect maxima (TRUE) or minima (FALSE)
 *
 * @return SEXP containing a list with extrema information
 */
SEXP S_detect_local_extrema(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_max_radius,
    SEXP s_min_neighborhood_size,
    SEXP s_detect_maxima
) {
    // Convert input from R format
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    std::vector<double> y(REAL(s_y), REAL(s_y) + Rf_length(s_y));
    double max_radius = REAL(s_max_radius)[0];
    size_t min_neighborhood_size = static_cast<size_t>(INTEGER(s_min_neighborhood_size)[0]);
    bool detect_maxima = LOGICAL(s_detect_maxima)[0];

    // Create graph
    set_wgraph_t graph(adj_list, weight_list);

    // Ensure graph diameter is available
    double graph_diameter = graph.graph_diameter;
    if (graph_diameter <= 0) {
        // Find diameter if not already set
        auto [end1, diam]     = graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = graph.get_vertex_eccentricity(end1);
        graph_diameter = diameter;
    }

    // Detect local extrema
    std::vector<local_extremum_t> extrema = graph.detect_local_extrema(
        y, max_radius, min_neighborhood_size, detect_maxima);

    // Prepare R output
    const char* names[] = {
        "vertices",
        "values",
        "radii",
        "neighborhood_sizes",
        "is_maxima",
        "neighborhood_vertices",
        "graph_diameter",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    int n_extrema = extrema.size();

    // Create result list
    int protect_count = 0;
    SEXP result = PROTECT(allocVector(VECSXP, n_elements)); protect_count++;

    // Set names
    SEXP result_names = PROTECT(allocVector(STRSXP, n_elements)); protect_count++;
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, mkChar(names[i]));
    }
    setAttrib(result, R_NamesSymbol, result_names);

    // Create data vectors
    SEXP r_vertices = PROTECT(allocVector(INTSXP, n_extrema)); protect_count++;
    SEXP r_values = PROTECT(allocVector(REALSXP, n_extrema)); protect_count++;
    SEXP r_radii = PROTECT(allocVector(REALSXP, n_extrema)); protect_count++;
    SEXP r_neighborhood_sizes = PROTECT(allocVector(INTSXP, n_extrema)); protect_count++;
    SEXP r_is_maxima = PROTECT(allocVector(LGLSXP, n_extrema)); protect_count++;

    // Create list to hold neighborhood vertices for each extremum
    SEXP r_neighborhood_vertices = PROTECT(allocVector(VECSXP, n_extrema)); protect_count++;

    // Fill data vectors
    for (int i = 0; i < n_extrema; i++) {
        INTEGER(r_vertices)[i] = extrema[i].vertex + 1;  // Convert to 1-based
        REAL(r_values)[i] = extrema[i].value;
        REAL(r_radii)[i] = extrema[i].radius;
        INTEGER(r_neighborhood_sizes)[i] = extrema[i].neighborhood_size;
        LOGICAL(r_is_maxima)[i] = extrema[i].is_maximum;

        // Add neighborhood vertices for this extremum
        size_t n_neighborhood = extrema[i].vertices.size();
        SEXP r_neighborhood = PROTECT(allocVector(INTSXP, n_neighborhood));

        for (size_t j = 0; j < n_neighborhood; j++) {
            INTEGER(r_neighborhood)[j] = extrema[i].vertices[j] + 1;  // Convert to 1-based
        }

        SET_VECTOR_ELT(r_neighborhood_vertices, i, r_neighborhood);
        UNPROTECT(1);  // Unprotect r_neighborhood after adding it to the list
    }

    // Create graph diameter and max packing radius values
    SEXP r_graph_diameter = PROTECT(allocVector(REALSXP, 1)); protect_count++;
    REAL(r_graph_diameter)[0] = graph_diameter;

    // Add vectors to result list
    SET_VECTOR_ELT(result, 0, r_vertices);
    SET_VECTOR_ELT(result, 1, r_values);
    SET_VECTOR_ELT(result, 2, r_radii);
    SET_VECTOR_ELT(result, 3, r_neighborhood_sizes);
    SET_VECTOR_ELT(result, 4, r_is_maxima);
    SET_VECTOR_ELT(result, 5, r_neighborhood_vertices);
    SET_VECTOR_ELT(result, 6, r_graph_diameter);

    UNPROTECT(protect_count);
    return result;
}
