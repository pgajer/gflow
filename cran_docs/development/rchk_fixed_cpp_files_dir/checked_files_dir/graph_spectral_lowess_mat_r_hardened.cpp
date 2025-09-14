/**
 * @brief Hardened rchk-safe drop-in of S_graph_spectral_lowess_mat
 *
 * Changes vs. provided version:
 * - Keep `s_result` and `s_names` protected until the tail; finish with UNPROTECT(2).
 * - Maintain container-first local PROTECT/UNPROTECT(1) for each element.
 * - Defensive scalars via Rf_as*; use explicit == TRUE for logicals.
 * - Use R_xlen_t casts when allocating matrices from size_t sizes.
 */
#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers
#include "graph_spectral_lowess.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations
struct graph_spectral_lowess_mat_t {
    std::vector<std::vector<double>> predictions;
    std::vector<std::vector<double>> errors;
    std::vector<std::vector<double>> scale;
};

class uniform_grid_graph_t {
public:
    double graph_diameter;
    
    uniform_grid_graph_t(const std::vector<std::vector<int>>& adj_list,
                         const std::vector<std::vector<double>>& weight_list);
    
    std::pair<size_t, double> get_vertex_eccentricity(size_t vertex);
    
    graph_spectral_lowess_mat_t graph_spectral_lowess_mat(
        const std::vector<std::vector<double>>& Y,
        size_t n_evectors,
        size_t n_bws,
        bool log_grid,
        double min_bw_factor,
        double max_bw_factor,
        double dist_normalization_factor,
        size_t kernel_type,
        double precision,
        size_t n_cleveland_iterations,
        bool with_errors,
        bool with_scale,
        bool verbose);
};

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
) {
    // Convert R data structures to C++ types
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    
    // Convert Y matrix/list from R
    std::vector<std::vector<double>> Y;
    
    if (Rf_isMatrix(s_Y)) {
        SEXP s_Y_dim = Rf_getAttrib(s_Y, R_DimSymbol);
        if (Rf_isNull(s_Y_dim) || TYPEOF(s_Y_dim) != INTSXP || LENGTH(s_Y_dim) != 2) {
            Rf_error("S_graph_spectral_lowess_mat: 'Y' must be a numeric matrix with proper dimensions");
        }
        const int nrow = INTEGER(s_Y_dim)[0];
        const int ncol = INTEGER(s_Y_dim)[1];
        if (nrow <= 0 || ncol <= 0) {
            Rf_error("S_graph_spectral_lowess_mat: invalid dimensions for 'Y'");
        }
        Y.resize((size_t)ncol);
        const double* Y_ptr = REAL(s_Y);
        for (int j = 0; j < ncol; ++j) {
            Y[(size_t)j].resize((size_t)nrow);
            for (int i = 0; i < nrow; ++i) {
                Y[(size_t)j][(size_t)i] = Y_ptr[i + (size_t)j * (size_t)nrow];
            }
        }
    } else if (Rf_isNewList(s_Y)) {
        const int n_cols = LENGTH(s_Y);
        if (n_cols <= 0) {
            Rf_error("S_graph_spectral_lowess_mat: 'Y' list must be non-empty");
        }
        Y.resize((size_t)n_cols);
        for (int j = 0; j < n_cols; ++j) {
            SEXP s_y_j = VECTOR_ELT(s_Y, j);
            if (!Rf_isReal(s_y_j)) {
                Rf_error("S_graph_spectral_lowess_mat: all elements of 'Y' list must be numeric vectors");
            }
            const R_xlen_t n_rows = XLENGTH(s_y_j);
            Y[(size_t)j].resize((size_t)n_rows);
            const double* y_ptr = REAL(s_y_j);
            for (R_xlen_t i = 0; i < n_rows; ++i) {
                Y[(size_t)j][(size_t)i] = y_ptr[i];
            }
        }
    } else {
        Rf_error("S_graph_spectral_lowess_mat: 'Y' must be a numeric matrix or a list of numeric vectors");
    }
    
    // Extract scalar parameters (defensive)
    const size_t n_evectors = (size_t) Rf_asInteger(s_n_evectors);
    const size_t n_bws = (size_t) Rf_asInteger(s_n_bws);
    const bool log_grid = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type = (size_t) Rf_asInteger(s_kernel_type);
    const double precision = Rf_asReal(s_precision);
    const size_t n_cleveland_iterations = (size_t) Rf_asInteger(s_n_cleveland_iterations);
    const bool with_errors = (Rf_asLogical(s_with_errors) == TRUE);
    const bool with_scale = (Rf_asLogical(s_with_scale) == TRUE);
    const bool verbose = (Rf_asLogical(s_verbose) == TRUE);
    
    // Create the graph object
    uniform_grid_graph_t grid_graph = uniform_grid_graph_t(adj_list, weight_list);

    // Find diameter endpoints
    auto [end1, diam]     = grid_graph.get_vertex_eccentricity(0);
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
    
    // Prepare result container (container-first). Keep names protected until tail.
    const int n_components = 1 + (with_errors ? 1 : 0) + (with_scale ? 1 : 0);
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));

    // Names
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));
    SET_STRING_ELT(s_names, 0, Rf_mkChar("predictions"));
    int name_idx = 1;
    if (with_errors) SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("errors"));
    if (with_scale)  SET_STRING_ELT(s_names, name_idx++, Rf_mkChar("scale"));
    Rf_setAttrib(s_result, R_NamesSymbol, s_names);
    
    // Dimensions
    const size_t n_response_vars = result.predictions.size();
    const size_t n_vertices      = n_response_vars ? result.predictions[0].size() : 0;

    // 0: predictions
    {
        SEXP s_predictions = PROTECT(Rf_allocMatrix(REALSXP, (R_xlen_t)n_vertices, (R_xlen_t)n_response_vars));
        double* pred_ptr = REAL(s_predictions);
        for (size_t j = 0; j < n_response_vars; ++j) {
            if (result.predictions[j].size() != n_vertices) {
                UNPROTECT(3); // s_predictions, s_result, s_names
                Rf_error("S_graph_spectral_lowess_mat: inconsistent prediction sizes");
            }
            for (size_t i = 0; i < n_vertices; ++i) {
                pred_ptr[i + j * n_vertices] = result.predictions[j][i];
            }
        }
        SET_VECTOR_ELT(s_result, 0, s_predictions);
        UNPROTECT(1); // s_predictions
    }
    
    int result_idx = 1;

    // errors (optional)
    if (with_errors) {
        SEXP s_errors = PROTECT(Rf_allocMatrix(REALSXP, (R_xlen_t)n_vertices, (R_xlen_t)n_response_vars));
        double* err_ptr = REAL(s_errors);
        for (size_t j = 0; j < n_response_vars; ++j) {
            if (result.errors[j].size() != n_vertices) {
                UNPROTECT(3); // s_errors, s_result, s_names
                Rf_error("S_graph_spectral_lowess_mat: inconsistent error sizes");
            }
            for (size_t i = 0; i < n_vertices; ++i) {
                err_ptr[i + j * n_vertices] = result.errors[j][i];
            }
        }
        SET_VECTOR_ELT(s_result, result_idx++, s_errors);
        UNPROTECT(1); // s_errors
    }
    
    // scale (optional)
    if (with_scale) {
        SEXP s_scale = PROTECT(Rf_allocMatrix(REALSXP, (R_xlen_t)n_vertices, (R_xlen_t)n_response_vars));
        double* scale_ptr = REAL(s_scale);
        for (size_t j = 0; j < n_response_vars; ++j) {
            if (result.scale[j].size() != n_vertices) {
                UNPROTECT(3); // s_scale, s_result, s_names
                Rf_error("S_graph_spectral_lowess_mat: inconsistent scale sizes");
            }
            for (size_t i = 0; i < n_vertices; ++i) {
                scale_ptr[i + j * n_vertices] = result.scale[j][i];
            }
        }
        SET_VECTOR_ELT(s_result, result_idx++, s_scale);
        UNPROTECT(1); // s_scale
    }
    
    UNPROTECT(2); // s_result, s_names
    return s_result;
}

} // extern "C"
