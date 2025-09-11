#include "uniform_grid_graph.hpp"
#include "graph_spectral_lowess_mat.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_graph_spectral_lowess_mat(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_Y,
        SEXP s_n_evectors,
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_precision,
        SEXP s_n_cleveland_iterations,
        SEXP s_with_errors,
        SEXP s_with_scale,
        SEXP s_verbose
    );
}

/**
 * @brief R interface for the matrix version of spectral-based locally weighted regression
 * 
 * This function serves as the bridge between R and the C++ implementation of the
 * graph_spectral_lowess_mat algorithm, handling data conversion and memory management.
 * 
 * @param s_adj_list R list representing the adjacency list of the graph
 * @param s_weight_list R list representing the edge weights
 * @param s_Y R list of numeric vectors or matrix, where each element is a response variable
 * @param s_n_evectors R integer for number of eigenvectors to use
 * @param s_n_bws R integer for number of bandwidths to evaluate
 * @param s_log_grid R logical for logarithmic bandwidth spacing
 * @param s_min_bw_factor R numeric for minimum bandwidth factor
 * @param s_max_bw_factor R numeric for maximum bandwidth factor
 * @param s_dist_normalization_factor R numeric for distance normalization
 * @param s_kernel_type R integer for kernel type
 * @param s_precision R numeric for precision in calculations
 * @param s_n_cleveland_iterations R integer for number of robust fitting iterations
 * @param s_with_errors R logical for including Rf_error estimates
 * @param s_with_scale R logical for including scale estimates
 * @param s_verbose R logical for verbose output
 * 
 * @return R list containing:
 *         - predictions: List of numeric vectors, one for each response variable
 *         - errors: List of numeric vectors of errors (if with_errors is TRUE)
 *         - scale: List of numeric vectors of scales (if with_scale is TRUE)
 */
SEXP S_graph_spectral_lowess_mat(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_n_evectors,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    SEXP s_precision,
    SEXP s_n_cleveland_iterations,
    SEXP s_with_errors,
    SEXP s_with_scale,
    SEXP s_verbose
) {
    // Convert R data structures to C++ types
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    
    // Convert Y matrix/list from R
    std::vector<std::vector<double>> Y;
    
    // Check if s_Y is a matrix or a list
    if (Rf_isMatrix(s_Y)) {
        // Handle matrix format (convert column-major R matrix to row-major C++ vector of vectors)
        SEXP s_Y_dim = Rf_getAttrib(s_Y, R_DimSymbol);
        int nrow = INTEGER(s_Y_dim)[0];
        int ncol = INTEGER(s_Y_dim)[1];
        
        // Each column of the matrix becomes a response variable vector
        Y.resize(ncol);
        double* Y_ptr = REAL(s_Y);
        
        for (int j = 0; j < ncol; j++) {
            Y[j].resize(nrow);
            for (int i = 0; i < nrow; i++) {
                // Access column-major matrix: Y_ptr[i + j*nrow]
                Y[j][i] = Y_ptr[i + j*nrow];
            }
        }
    } else if (Rf_isNewList(s_Y)) {
        // Handle list format
        int n_cols = LENGTH(s_Y);
        Y.resize(n_cols);
        
        for (int j = 0; j < n_cols; j++) {
            SEXP s_y_j = VECTOR_ELT(s_Y, j);
            if (!Rf_isReal(s_y_j)) {
                Rf_error("All elements of Y must be numeric vectors");
            }
            
            int n_rows = LENGTH(s_y_j);
            Y[j].resize(n_rows);
            double* y_ptr = REAL(s_y_j);
            
            for (int i = 0; i < n_rows; i++) {
                Y[j][i] = y_ptr[i];
            }
        }
    } else {
        Rf_error("Y must be a numeric matrix or a list of numeric vectors");
    }
    
    // Extract scalar parameters
    size_t n_evectors = static_cast<size_t>(Rf_asInteger(s_n_evectors));
    size_t n_bws = static_cast<size_t>(Rf_asInteger(s_n_bws));
    bool log_grid = Rf_asLogical(s_log_grid);
    double min_bw_factor = Rf_asReal(s_min_bw_factor);
    double max_bw_factor = Rf_asReal(s_max_bw_factor);
    double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    size_t kernel_type = static_cast<size_t>(Rf_asInteger(s_kernel_type));
    double precision = Rf_asReal(s_precision);
    size_t n_cleveland_iterations = static_cast<size_t>(Rf_asInteger(s_n_cleveland_iterations));
    bool with_errors = Rf_asLogical(s_with_errors);
    bool with_scale = Rf_asLogical(s_with_scale);
    bool verbose = Rf_asLogical(s_verbose);
    
    // Create the graph object
    uniform_grid_graph_t grid_graph = uniform_grid_graph_t(adj_list, weight_list);

    // Find diameter endpoints
    auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
    auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
    grid_graph.graph_diameter = diameter;

    // Call the C++ function for matrix-based spectral lowess
    graph_spectral_lowess_mat_t result = grid_graph.graph_spectral_lowess_mat(
        Y,
        n_evectors,
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        dist_normalization_factor,
        kernel_type,
        precision,
        n_cleveland_iterations,
        with_errors,
        with_scale,
        verbose
    );
    
    // Initialize return list with proper components
    // Basic components + optional Rf_error and scale components
    int n_components = 1 + (with_errors ? 1 : 0) + (with_scale ? 1 : 0);
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));
    
    // Set component names
    SET_STRING_ELT(s_names, 0, Rf_mkChar("predictions"));
    int name_idx = 1;
    if (with_errors) {
        SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("errors"));
    }
    if (with_scale) {
        SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("scale"));
    }
    Rf_setAttrib(s_result, R_NamesSymbol, s_names);
    
    // Convert predictions to R
    size_t n_response_vars = result.predictions.size();
    size_t n_vertices = result.predictions[0].size();
    
    SEXP s_predictions = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
    double* pred_ptr = REAL(s_predictions);
    
    for (size_t j = 0; j < n_response_vars; j++) {
        for (size_t i = 0; i < n_vertices; i++) {
            // Column-major format for R matrix
            pred_ptr[i + j*n_vertices] = result.predictions[j][i];
        }
    }
    SET_VECTOR_ELT(s_result, 0, s_predictions);
    
    // Convert errors to R if requested
    int result_idx = 1;
    if (with_errors) {
        SEXP s_errors = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
        double* err_ptr = REAL(s_errors);
        
        for (size_t j = 0; j < n_response_vars; j++) {
            for (size_t i = 0; i < n_vertices; i++) {
                // Column-major format for R matrix
                err_ptr[i + j*n_vertices] = result.errors[j][i];
            }
        }
        SET_VECTOR_ELT(s_result, result_idx++, s_errors);
    }
    
    // Convert scale to R if requested
    if (with_scale) {
        SEXP s_scale = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
        double* scale_ptr = REAL(s_scale);
        
        for (size_t j = 0; j < n_response_vars; j++) {
            for (size_t i = 0; i < n_vertices; i++) {
                // Column-major format for R matrix
                scale_ptr[i + j*n_vertices] = result.scale[j][i];
            }
        }
        SET_VECTOR_ELT(s_result, result_idx++, s_scale);
    }
    
    // Calculate total number of protected objects for unprotecting
    int n_protected = 2; // s_result and s_names
    n_protected += 1;    // s_predictions
    if (with_errors) n_protected += 1; // s_errors
    if (with_scale) n_protected += 1;  // s_scale
    
    UNPROTECT(n_protected);
    return s_result;
}
