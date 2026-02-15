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
    // Convert R data structures to C++ types (no PROTECT activity)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert Y matrix/list from R
    std::vector<std::vector<double>> Y;

    if (Rf_isMatrix(s_Y)) {
        SEXP s_dim = PROTECT(Rf_getAttrib(s_Y, R_DimSymbol));
        if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) < 2) {
            UNPROTECT(1);
            Rf_error("Y must be a matrix with a valid integer 'dim' attribute.");
        }
        const int nrow = INTEGER(s_dim)[0];
        const int ncol = INTEGER(s_dim)[1];
        UNPROTECT(1); // s_dim

        Y.resize(ncol);
        const double* Y_ptr = REAL(s_Y);
        for (int j = 0; j < ncol; ++j) {
            Y[j].resize(nrow);
            for (int i = 0; i < nrow; ++i) {
                Y[j][i] = Y_ptr[i + j * nrow]; // column-major
            }
        }
    } else if (Rf_isNewList(s_Y)) {
        const int n_cols = LENGTH(s_Y);
        Y.resize(n_cols);
        for (int j = 0; j < n_cols; ++j) {
            SEXP s_yj = VECTOR_ELT(s_Y, j);
            if (!Rf_isReal(s_yj)) Rf_error("All elements of Y must be numeric vectors");
            const int n_rows = LENGTH(s_yj);
            Y[j].resize(n_rows);
            const double* y_ptr = REAL(s_yj);
            for (int i = 0; i < n_rows; ++i) Y[j][i] = y_ptr[i];
        }
    } else {
        Rf_error("Y must be a numeric matrix or a list of numeric vectors");
    }

    // Scalars
    const int n_evectors = Rf_asInteger(s_n_evectors);
    const int n_bws      = Rf_asInteger(s_n_bws);
    const int log_grid   = Rf_asLogical(s_log_grid);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const int kernel_type = Rf_asInteger(s_kernel_type);
    const double precision = Rf_asReal(s_precision);
    const int n_cleveland_iterations = Rf_asInteger(s_n_cleveland_iterations);
    const int with_errors = Rf_asLogical(s_with_errors);
    const int with_scale  = Rf_asLogical(s_with_scale);
    const int verbose     = Rf_asLogical(s_verbose);

    // Build graph and compute
    uniform_grid_graph_t grid_graph(adj_list, weight_list);
    auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);
    auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
    grid_graph.graph_diameter = diameter;

    graph_spectral_lowess_mat_t result = grid_graph.graph_spectral_lowess_mat(
        Y,
        (size_t)n_evectors,
        (size_t)n_bws,
        (bool)log_grid,
        min_bw_factor,
        max_bw_factor,
        dist_normalization_factor,
        (size_t)kernel_type,
        precision,
        (size_t)n_cleveland_iterations,
        (bool)with_errors,
        (bool)with_scale,
        (bool)verbose
    );

    // Parent list
    int n_components = 1 + (with_errors ? 1 : 0) + (with_scale ? 1 : 0);
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));

    // Names
    {
        SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));
        int name_idx = 0;
        SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("predictions"));
        if (with_errors) SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("errors"));
        if (with_scale)  SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("scale"));
        Rf_setAttrib(s_result, R_NamesSymbol, s_names);
        UNPROTECT(1); // s_names
    }

    // Dimensions for matrices
    const int n_response_vars = (int) result.predictions.size();
    const int n_vertices = n_response_vars > 0 ? (int) result.predictions[0].size() : 0;

    // predictions
    {
        SEXP s_predictions = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
        double* pred_ptr = REAL(s_predictions);
        for (int j = 0; j < n_response_vars; ++j) {
            for (int i = 0; i < n_vertices; ++i) {
                pred_ptr[i + (size_t)j * n_vertices] = result.predictions[(size_t)j][(size_t)i];
            }
        }
        SET_VECTOR_ELT(s_result, 0, s_predictions);
        UNPROTECT(1); // s_predictions
    }

    int slot = 1;

    // errors (optional)
    if (with_errors) {
        SEXP s_errors = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
        double* err_ptr = REAL(s_errors);
        for (int j = 0; j < n_response_vars; ++j) {
            for (int i = 0; i < n_vertices; ++i) {
                err_ptr[i + (size_t)j * n_vertices] = result.errors[(size_t)j][(size_t)i];
            }
        }
        SET_VECTOR_ELT(s_result, slot++, s_errors);
        UNPROTECT(1); // s_errors
    }

    // scale (optional)
    if (with_scale) {
        SEXP s_scale = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_response_vars));
        double* scale_ptr = REAL(s_scale);
        for (int j = 0; j < n_response_vars; ++j) {
            for (int i = 0; i < n_vertices; ++i) {
                scale_ptr[i + (size_t)j * n_vertices] = result.scale[(size_t)j][(size_t)i];
            }
        }
        SET_VECTOR_ELT(s_result, slot++, s_scale);
        UNPROTECT(1); // s_scale
    }

    UNPROTECT(1); // s_result
    return s_result;
}
