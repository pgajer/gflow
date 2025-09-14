/**
 * @brief Fixed version of centered_paths.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_ugg_get_path_data
 * 
 * Issues fixed:
 * 1. S_ugg_get_path_data (line 2700): UNPROTECT(variable) - unsupported
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constants
 * 2. Used container-first pattern consistently
 * 3. Fixed all PROTECT/UNPROTECT imbalances
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_matrix_double_to_R(const std::vector<std::vector<double>>&);

// Core structures and declarations
struct path_data_t {
    std::vector<std::vector<double>> path_coordinates;
    std::vector<int> path_indices;
    std::vector<double> path_distances;
    std::vector<double> path_curvatures;
    double total_length;
    int n_points;
    bool is_closed;
};

// Core computation function (assumed available)
path_data_t ugg_get_path_data(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<std::vector<double>>& coordinates,
    int start_vertex,
    int end_vertex,
    double smoothing_bandwidth,
    bool compute_curvature,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_ugg_get_path_data
 * Fixes: UNPROTECT(variable) at line 2700
 * Solution: Use container-first pattern, literal UNPROTECT
 */
SEXP S_ugg_get_path_data(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP coordinates_sexp,
    SEXP start_vertex_sexp,
    SEXP end_vertex_sexp,
    SEXP smoothing_bandwidth_sexp,
    SEXP compute_curvature_sexp,
    SEXP verbose_sexp
) {
    // Input coercion block with PROTECT_WITH_INDEX
    PROTECT_INDEX ipx;
    SEXP coords_real = coordinates_sexp;
    PROTECT_WITH_INDEX(coords_real, &ipx);
    
    if (TYPEOF(coords_real) != REALSXP) {
        REPROTECT(coords_real = Rf_coerceVector(coords_real, REALSXP), ipx);
    }
    
    // Get dimensions
    SEXP dims = Rf_getAttrib(coords_real, R_DimSymbol);
    if (Rf_isNull(dims) || LENGTH(dims) != 2) {
        UNPROTECT(1); // coords_real
        Rf_error("coordinates must be a matrix");
    }
    
    int* dimptr = INTEGER(dims);
    int n_points = dimptr[0];
    int n_dims = dimptr[1];
    
    // Convert coordinates to C++ format
    std::vector<std::vector<double>> coordinates(n_points);
    double* coord_ptr = REAL(coords_real);
    for (int i = 0; i < n_points; i++) {
        coordinates[i].resize(n_dims);
        for (int j = 0; j < n_dims; j++) {
            coordinates[i][j] = coord_ptr[i + n_points * j]; // Column-major to row-major
        }
    }
    
    // Convert other inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    // Extract scalar parameters
    int start_vertex = Rf_asInteger(start_vertex_sexp) - 1; // Convert to 0-based
    int end_vertex = Rf_asInteger(end_vertex_sexp) - 1; // Convert to 0-based
    double smoothing_bandwidth = Rf_asReal(smoothing_bandwidth_sexp);
    bool compute_curvature = (Rf_asLogical(compute_curvature_sexp) == TRUE);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Validate vertex indices
    if (start_vertex < 0 || start_vertex >= n_points) {
        UNPROTECT(1); // coords_real
        Rf_error("start_vertex out of bounds");
    }
    if (end_vertex < 0 || end_vertex >= n_points) {
        UNPROTECT(1); // coords_real
        Rf_error("end_vertex out of bounds");
    }
    
    // Core computation
    path_data_t path_data = ugg_get_path_data(
        adj_list,
        weight_list,
        coordinates,
        start_vertex,
        end_vertex,
        smoothing_bandwidth,
        compute_curvature,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: path_coordinates (matrix)
    {
        SEXP path_coords = PROTECT(convert_matrix_double_to_R(path_data.path_coordinates));
        SET_VECTOR_ELT(result, 0, path_coords);
        UNPROTECT(1);
    }
    
    // 1: path_indices (convert to 1-based for R)
    {
        std::vector<int> r_indices = path_data.path_indices;
        for (size_t i = 0; i < r_indices.size(); i++) {
            r_indices[i]++; // Convert to 1-based
        }
        SEXP indices = PROTECT(convert_vector_int_to_R(r_indices));
        SET_VECTOR_ELT(result, 1, indices);
        UNPROTECT(1);
    }
    
    // 2: path_distances
    {
        SEXP distances = PROTECT(convert_vector_double_to_R(path_data.path_distances));
        SET_VECTOR_ELT(result, 2, distances);
        UNPROTECT(1);
    }
    
    // 3: path_curvatures
    {
        SEXP curvatures = PROTECT(convert_vector_double_to_R(path_data.path_curvatures));
        SET_VECTOR_ELT(result, 3, curvatures);
        UNPROTECT(1);
    }
    
    // 4: total_length
    {
        SEXP length = PROTECT(Rf_ScalarReal(path_data.total_length));
        SET_VECTOR_ELT(result, 4, length);
        UNPROTECT(1);
    }
    
    // 5: n_points
    {
        SEXP n_pts = PROTECT(Rf_ScalarInteger(path_data.n_points));
        SET_VECTOR_ELT(result, 5, n_pts);
        UNPROTECT(1);
    }
    
    // 6: is_closed
    {
        SEXP closed = PROTECT(Rf_ScalarLogical(path_data.is_closed ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 6, closed);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("path_coordinates"));
        SET_STRING_ELT(names, 1, Rf_mkChar("path_indices"));
        SET_STRING_ELT(names, 2, Rf_mkChar("path_distances"));
        SET_STRING_ELT(names, 3, Rf_mkChar("path_curvatures"));
        SET_STRING_ELT(names, 4, Rf_mkChar("total_length"));
        SET_STRING_ELT(names, 5, Rf_mkChar("n_points"));
        SET_STRING_ELT(names, 6, Rf_mkChar("is_closed"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(2); // coords_real, result
    return result;
}

} // extern "C"