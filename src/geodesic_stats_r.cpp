#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <numeric> // for std::iota()

#include "geodesic_stats.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

/**
 * R interface for compute_geodesic_stats C++ function
 */
extern "C" SEXP S_compute_geodesic_stats(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP min_radius_sexp,
    SEXP max_radius_sexp,
    SEXP n_steps_sexp,
    SEXP n_packing_vertices_sexp,
    SEXP max_packing_iterations_sexp,
    SEXP packing_precision_sexp,
    SEXP verbose_sexp
) {
    // Ensure R won't interrupt our C++ code
    R_CheckUserInterrupt();

    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    // Extract scalar parameters
    double min_radius = Rf_asReal(min_radius_sexp);
    double max_radius = Rf_asReal(max_radius_sexp);
    size_t n_steps = Rf_asInteger(n_steps_sexp);
    size_t n_packing_vertices = Rf_asInteger(n_packing_vertices_sexp);
    size_t max_packing_iterations = Rf_asInteger(max_packing_iterations_sexp);
    double packing_precision = Rf_asReal(packing_precision_sexp);
    bool verbose = Rf_asLogical(verbose_sexp);

    uniform_grid_graph_t grid_graph;
    size_t n_vertices = adj_list.size();
    if (n_packing_vertices < n_vertices) {
        // The returned uniform_grid_graph_t object inherits both the graph_diameter
        // and max_packing_radius values from the intermediate set_wgraph_t object,
        // making these calculated values available for further analysis.
        grid_graph = create_maximal_packing(adj_list,
                                            weight_list,
                                            n_packing_vertices,
                                            max_packing_iterations,
                                            packing_precision);
    } else {
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), 0);

        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);

        // Find diameter endpoints
        auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
        grid_graph.graph_diameter = diameter;
    }

    // Compute geodesic statistics
    geodesic_stats_t stats = compute_geodesic_stats(
        grid_graph,
        min_radius,
        max_radius,
        n_steps,
        verbose
    );

    // Create R list to return results
    size_t n_protected = 0;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 6)); n_protected++;

        // Get the size of the data
    size_t n_radii = stats.radii.size();

    // Convert grid_vertices to vector for more convenient handling
    std::vector<size_t> grid_vertices_vec(
        grid_graph.grid_vertices.begin(),
        grid_graph.grid_vertices.end()
        );

    // 1. Radii
    SEXP radii_sexp = PROTECT(Rf_allocVector(REALSXP, n_radii));
    double* radii_ptr = REAL(radii_sexp);
    for (size_t i = 0; i < n_radii; i++) {
        radii_ptr[i] = stats.radii[i];
    }
    SET_VECTOR_ELT(result, 0, radii_sexp);
    UNPROTECT(1);
    
    // 2. Grid vertices (returning 1-based indices for R)
    size_t n_grid_vertices = grid_vertices_vec.size();
    SEXP grid_vertices_sexp = PROTECT(Rf_allocVector(INTSXP, n_grid_vertices));
    int* grid_vertices_ptr = INTEGER(grid_vertices_sexp);
    for (size_t i = 0; i < n_grid_vertices; i++) {
        grid_vertices_ptr[i] = grid_vertices_vec[i] + 1; // 1-based for R
    }
    SET_VECTOR_ELT(result, 1, grid_vertices_sexp);
    UNPROTECT(1);

    // 3. Create matrices for geodesic rays, composite geodesics, and path overlap ratio
    // First create the matrices with all zeros/NA

    // Geodesic rays matrix (vertices x radii)
    SEXP geodesic_rays_matrix = PROTECT(Rf_allocMatrix(INTSXP, n_grid_vertices, n_radii)); n_protected++;
    int* geodesic_rays_ptr = INTEGER(geodesic_rays_matrix);

    // Composite geodesics matrix
    SEXP composite_geodesics_matrix = PROTECT(Rf_allocMatrix(INTSXP, n_grid_vertices, n_radii)); n_protected++;
    int* composite_geodesics_ptr = INTEGER(composite_geodesics_matrix);

    // Overlap statistics matrices - 7 statistics per vertex/radius
    SEXP overlap_min_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_p05_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_p25_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_median_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_p75_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_p95_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;
    SEXP overlap_max_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_grid_vertices, n_radii)); n_protected++;

    double* overlap_min_ptr = REAL(overlap_min_matrix);
    double* overlap_p05_ptr = REAL(overlap_p05_matrix);
    double* overlap_p25_ptr = REAL(overlap_p25_matrix);
    double* overlap_median_ptr = REAL(overlap_median_matrix);
    double* overlap_p75_ptr = REAL(overlap_p75_matrix);
    double* overlap_p95_ptr = REAL(overlap_p95_matrix);
    double* overlap_max_ptr = REAL(overlap_max_matrix);

    // Initialize matrices (set all elements to 0)
    for (size_t i = 0; i < n_grid_vertices * n_radii; i++) {
        geodesic_rays_ptr[i] = -1; // Default to 1 ray
        composite_geodesics_ptr[i] = 0;

        // Initialize overlap statistics to 0
        overlap_min_ptr[i] = 0.0;
        overlap_p05_ptr[i] = 0.0;
        overlap_p25_ptr[i] = 0.0;
        overlap_median_ptr[i] = 0.0;
        overlap_p75_ptr[i] = 0.0;
        overlap_p95_ptr[i] = 0.0;
        overlap_max_ptr[i] = 0.0;
    }

    // Fill matrix values
    for (size_t r = 0; r < n_radii; r++) {
        // For each radius, get the maps
        const auto& rays_map = stats.geodesic_rays[r];
        const auto& comp_map = stats.composite_geodesics[r];
        const auto& overlap_map = stats.paths_overlap[r];

        // For each grid vertex in our vectors
        for (size_t i = 0; i < n_grid_vertices; i++) {
            size_t vertex = grid_vertices_vec[i];
            // Matrix position - column-major layout in R
            size_t pos = i + r * n_grid_vertices;

            // Set values if they exist in the maps, otherwise leave default values
            if (rays_map.count(vertex) > 0) {
                geodesic_rays_ptr[pos] = rays_map.at(vertex);
            }

            if (comp_map.count(vertex) > 0) {
                composite_geodesics_ptr[pos] = comp_map.at(vertex);
            }

            if (overlap_map.count(vertex) > 0) {
                const auto& stats = overlap_map.at(vertex);
                overlap_min_ptr[pos] = stats.min;
                overlap_p05_ptr[pos] = stats.p05;
                overlap_p25_ptr[pos] = stats.p25;
                overlap_median_ptr[pos] = stats.median;
                overlap_p75_ptr[pos] = stats.p75;
                overlap_p95_ptr[pos] = stats.p95;
                overlap_max_ptr[pos] = stats.max;
            }
        }
    }
    
    // Set row names to grid vertices
    SEXP rownames = PROTECT(Rf_allocVector(STRSXP, n_grid_vertices)); n_protected++;
    for (size_t i = 0; i < n_grid_vertices; i++) {
        SET_STRING_ELT(rownames, i, Rf_mkChar(std::to_string(grid_vertices_vec[i] + 1).c_str()));
    }

    // Set column names to radii
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, n_radii)); n_protected++;
    for (size_t r = 0; r < n_radii; r++) {
        char radius_str[32];
        //snprintf(radius_str, sizeof(radius_str), "%.4f", stats.radii[r]);
        snprintf(radius_str, sizeof(radius_str), "%zu", r + 1);
        SET_STRING_ELT(colnames, r, Rf_mkChar(radius_str));
    }

    // Create dimnames list (row names, column names)
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2)); n_protected++;
    SET_VECTOR_ELT(dimnames, 0, rownames);
    SET_VECTOR_ELT(dimnames, 1, colnames);

    // Set dimnames for all matrices
    Rf_setAttrib(geodesic_rays_matrix, R_DimNamesSymbol, dimnames);
    Rf_setAttrib(composite_geodesics_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_min_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_p05_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_p25_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_median_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_p75_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_p95_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));
    Rf_setAttrib(overlap_max_matrix, R_DimNamesSymbol, Rf_duplicate(dimnames));

    // Combine overlap statistics matrices into a list
    SEXP overlap_list = PROTECT(Rf_allocVector(VECSXP, 7)); n_protected++;
    SET_VECTOR_ELT(overlap_list, 0, overlap_min_matrix);
    SET_VECTOR_ELT(overlap_list, 1, overlap_p05_matrix);
    SET_VECTOR_ELT(overlap_list, 2, overlap_p25_matrix);
    SET_VECTOR_ELT(overlap_list, 3, overlap_median_matrix);
    SET_VECTOR_ELT(overlap_list, 4, overlap_p75_matrix);
    SET_VECTOR_ELT(overlap_list, 5, overlap_p95_matrix);
    SET_VECTOR_ELT(overlap_list, 6, overlap_max_matrix);

    // Set names for overlap list
    SEXP overlap_names = PROTECT(Rf_allocVector(STRSXP, 7)); n_protected++;
    SET_STRING_ELT(overlap_names, 0, Rf_mkChar("min"));
    SET_STRING_ELT(overlap_names, 1, Rf_mkChar("p05"));
    SET_STRING_ELT(overlap_names, 2, Rf_mkChar("p25"));
    SET_STRING_ELT(overlap_names, 3, Rf_mkChar("median"));
    SET_STRING_ELT(overlap_names, 4, Rf_mkChar("p75"));
    SET_STRING_ELT(overlap_names, 5, Rf_mkChar("p95"));
    SET_STRING_ELT(overlap_names, 6, Rf_mkChar("max"));
    Rf_setAttrib(overlap_list, R_NamesSymbol, overlap_names);

    // Add matrices to result
    SET_VECTOR_ELT(result, 2, geodesic_rays_matrix);
    SET_VECTOR_ELT(result, 3, composite_geodesics_matrix);
    SET_VECTOR_ELT(result, 4, overlap_list);

    SEXP radius_overlaps_list = PROTECT(Rf_allocVector(VECSXP, n_radii)); n_protected++;

    for (size_t r = 0; r < n_radii; r++) {
        size_t n_values = stats.radius_overlaps[r].size();
        SEXP overlaps_vector = PROTECT(Rf_allocVector(REALSXP, n_values));

        double* overlaps_ptr = REAL(overlaps_vector);
        for (size_t i = 0; i < n_values; i++) {
            overlaps_ptr[i] = stats.radius_overlaps[r][i];
        }

        SET_VECTOR_ELT(radius_overlaps_list, r, overlaps_vector);
        UNPROTECT(1);
    }

    SET_VECTOR_ELT(result, 5, radius_overlaps_list);


    // Set names for result list
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, 6)); n_protected++;
    SET_STRING_ELT(result_names, 0, Rf_mkChar("radii"));
    SET_STRING_ELT(result_names, 1, Rf_mkChar("grid_vertices"));
    SET_STRING_ELT(result_names, 2, Rf_mkChar("geodesic_rays"));
    SET_STRING_ELT(result_names, 3, Rf_mkChar("composite_geodesics"));
    SET_STRING_ELT(result_names, 4, Rf_mkChar("path_overlap"));
    SET_STRING_ELT(result_names, 5, Rf_mkChar("radius_overlaps"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // Clean up
    UNPROTECT(n_protected); // All PROTECT calls

    return result;
}

/**
 * R interface for compute_vertex_geodesic_stats C++ function
 */
extern "C" SEXP S_compute_vertex_geodesic_stats(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP grid_vertex_sexp,
    SEXP min_radius_sexp,
    SEXP max_radius_sexp,
    SEXP n_steps_sexp,
    SEXP n_packing_vertices_sexp,
    SEXP packing_precision_sexp
) {
    // Ensure R won't interrupt our C++ code
    R_CheckUserInterrupt();

    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    // Extract scalar parameters
    size_t grid_vertex = Rf_asInteger(grid_vertex_sexp);
    double min_radius = Rf_asReal(min_radius_sexp);
    double max_radius = Rf_asReal(max_radius_sexp);
    int n_steps = Rf_asInteger(n_steps_sexp);
    size_t n_packing_vertices = Rf_asInteger(n_packing_vertices_sexp);
    double packing_precision = Rf_asReal(packing_precision_sexp);

    uniform_grid_graph_t grid_graph;
    size_t n_vertices = adj_list.size();
    if (n_packing_vertices < n_vertices) {
        size_t max_packing_iterations = 20;
        // The returned uniform_grid_graph_t object inherits both the graph_diameter
        // and max_packing_radius values from the intermediate set_wgraph_t object,
        // making these calculated values available for further analysis.
        grid_graph = create_maximal_packing(adj_list,
                                            weight_list,
                                            n_packing_vertices,
                                            max_packing_iterations,
                                            packing_precision);
    } else {
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), 0);

        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);

        // Find diameter endpoints
        auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
        grid_graph.graph_diameter = diameter;
    }

    // Compute vertex-specific geodesic statistics
    auto stats = compute_vertex_geodesic_stats(
        grid_graph,
        grid_vertex,
        min_radius,
        max_radius,
        n_steps
    );

    // Get the size of the data
    size_t n_radii = stats.size();

    // Create result matrix with 9 columns:
    // [radius, geodesic_rays, composite_geodesics, overlap_min, overlap_p05, overlap_p25, overlap_median, overlap_p75, overlap_p95, overlap_max]
    SEXP result_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_radii, 10));
    double* result_ptr = REAL(result_matrix);

    // Fill the matrix
    for (size_t i = 0; i < n_radii; i++) {
        auto [radius, geodesic_rays, composite_geodesics, overlap_stats] = stats[i];

        // Column-major layout in R
        result_ptr[i] = radius;  // Column 1: radius
        result_ptr[i + n_radii] = static_cast<double>(geodesic_rays);  // Column 2: geodesic_rays
        result_ptr[i + 2 * n_radii] = static_cast<double>(composite_geodesics);  // Column 3: composite_geodesics

        // Overlap statistics
        result_ptr[i + 3 * n_radii] = overlap_stats.min;  // Column 4: overlap_min
        result_ptr[i + 4 * n_radii] = overlap_stats.p05;  // Column 5: overlap_p05
        result_ptr[i + 5 * n_radii] = overlap_stats.p25;  // Column 6: overlap_p25
        result_ptr[i + 6 * n_radii] = overlap_stats.median;  // Column 7: overlap_median
        result_ptr[i + 7 * n_radii] = overlap_stats.p75;  // Column 8: overlap_p75
        result_ptr[i + 8 * n_radii] = overlap_stats.p95;  // Column 9: overlap_p95
        result_ptr[i + 9 * n_radii] = overlap_stats.max;  // Column 10: overlap_max
    }

    // Set column names
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 10));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("radius"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("rays"));
    SET_STRING_ELT(colnames, 2, Rf_mkChar("composite_geodesics"));
    SET_STRING_ELT(colnames, 3, Rf_mkChar("overlap_min"));
    SET_STRING_ELT(colnames, 4, Rf_mkChar("overlap_p05"));
    SET_STRING_ELT(colnames, 5, Rf_mkChar("overlap_p25"));
    SET_STRING_ELT(colnames, 6, Rf_mkChar("overlap_median"));
    SET_STRING_ELT(colnames, 7, Rf_mkChar("overlap_p75"));
    SET_STRING_ELT(colnames, 8, Rf_mkChar("overlap_p95"));
    SET_STRING_ELT(colnames, 9, Rf_mkChar("overlap_max"));

    // Create dimnames list (NULL rownames, column names)
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(dimnames, 1, colnames);

    // Set dimnames for matrix
    Rf_setAttrib(result_matrix, R_DimNamesSymbol, dimnames);

    // Create list to return the matrix and additional info
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, result_matrix);
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(grid_vertex + 1)); // Return 1-based vertex

    // Set names for result list
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, Rf_mkChar("data"));
    SET_STRING_ELT(result_names, 1, Rf_mkChar("vertex"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // Clean up
    UNPROTECT(5); // result_matrix, colnames, dimnames, result, result_names

    return result;
}
