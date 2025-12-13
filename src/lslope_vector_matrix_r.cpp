// lslope_vector_matrix_r.cpp
//
// R interface for lslope_vector_matrix functionality.
// This should be appended to lslope_r.cpp or compiled separately.

#include "set_wgraph.hpp"
#include "lslope.hpp"
#include "lcor.hpp"

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include <R.h>
#include <Rinternals.h>

// Forward declare conversion utilities
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// Sentinel for invalid vertex
static const size_t LSLOPE_R_VM_INVALID_VERTEX = std::numeric_limits<size_t>::max();

/**
 * @brief SEXP interface for computing local slope between vector y and matrix Z
 *
 * @details
 * This function provides an R interface to the optimized vector-matrix local
 * slope computation. It pre-computes gradient edges from y once and reuses
 * them for all columns of Z, with OpenMP parallelization over columns.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector of directing function values
 * @param s_Z R numeric matrix of response function values
 * @param s_type R character: "slope", "normalized", or "sign"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
 * @param s_sigmoid_alpha R numeric: sigmoid scale (0 = auto-calibrate)
 * @param s_ascending R logical: use ascending (TRUE) or descending (FALSE) gradient
 * @param s_n_threads R integer: number of OpenMP threads (0 = default)
 *
 * @return R list with components:
 *   - coefficients: Matrix (n_vertices x n_columns) of local slopes
 *   - gradient.neighbors: Integer vector of gradient edge neighbors (NA if extremum)
 *   - gradient.delta.y: Numeric vector of delta y along gradient edges
 *   - is.local.extremum: Logical vector indicating local extrema
 *   - n.local.maxima: Count of local maxima
 *   - n.local.minima: Count of local minima
 *   - sigmoid.alpha: Sigmoid scale used (for normalized type)
 */
extern "C" SEXP S_lslope_vector_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_Z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_sigmoid_alpha,
    SEXP s_ascending,
    SEXP s_n_threads
) {
    // ---- Input validation and conversion ----

    // Convert adjacency and weight lists
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y vector
    if (!Rf_isReal(s_y)) {
        Rf_error("s_y must be a numeric vector");
    }
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    // Convert Z matrix
    if (!Rf_isMatrix(s_Z) || !Rf_isReal(s_Z)) {
        Rf_error("s_Z must be a numeric matrix");
    }
    int n_rows = Rf_nrows(s_Z);
    int n_cols = Rf_ncols(s_Z);
    double* Z_ptr = REAL(s_Z);

    // Create Eigen matrix (column-major like R)
    Eigen::Map<Eigen::MatrixXd> Z(Z_ptr, n_rows, n_cols);

    // Validate dimensions
    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices) {
        Rf_error("Length of y (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_y), static_cast<int>(n_vertices));
    }
    if (static_cast<size_t>(n_rows) != n_vertices) {
        Rf_error("Number of rows in Z (%d) does not match number of vertices (%d)",
                 n_rows, static_cast<int>(n_vertices));
    }

    // Convert type string
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);

    lslope_type_t slope_type;
    if (type_str == "slope") {
        slope_type = lslope_type_t::GRADIENT_SLOPE;
    } else if (type_str == "normalized") {
        slope_type = lslope_type_t::GRADIENT_SLOPE_NORMALIZED;
    } else if (type_str == "sign") {
        slope_type = lslope_type_t::GRADIENT_SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'slope', 'normalized', or 'sign'",
                 type_cstr);
    }

    // Convert y_diff_type
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    const char* y_diff_cstr = CHAR(STRING_ELT(s_y_diff_type, 0));
    std::string y_diff_str(y_diff_cstr);
    edge_diff_type_t y_diff_type = (y_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    // Convert z_diff_type
    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    const char* z_diff_cstr = CHAR(STRING_ELT(s_z_diff_type, 0));
    std::string z_diff_str(z_diff_cstr);
    edge_diff_type_t z_diff_type = (z_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    // Convert epsilon
    double epsilon = 0.0;
    if (Rf_isReal(s_epsilon) && LENGTH(s_epsilon) == 1) {
        epsilon = REAL(s_epsilon)[0];
    }

    // Convert sigmoid_alpha
    double sigmoid_alpha = 0.0;
    if (Rf_isReal(s_sigmoid_alpha) && LENGTH(s_sigmoid_alpha) == 1) {
        sigmoid_alpha = REAL(s_sigmoid_alpha)[0];
    }

    // Convert ascending
    bool ascending = true;
    if (Rf_isLogical(s_ascending) && LENGTH(s_ascending) == 1) {
        ascending = (LOGICAL(s_ascending)[0] == TRUE);
    }

    // Convert n_threads
    int n_threads = 0;
    if (Rf_isInteger(s_n_threads) && LENGTH(s_n_threads) == 1) {
        n_threads = INTEGER(s_n_threads)[0];
    } else if (Rf_isReal(s_n_threads) && LENGTH(s_n_threads) == 1) {
        n_threads = static_cast<int>(REAL(s_n_threads)[0]);
    }

    // ---- Build graph ----
    set_wgraph_t graph(adj_list, weight_list);

    // ---- Compute local slopes ----
    lslope_vector_matrix_result_t result = graph.lslope_vector_matrix(
        y,
        Z,
        slope_type,
        y_diff_type,
        z_diff_type,
        epsilon,
        sigmoid_alpha,
        ascending,
        n_threads
    );

    // ---- Build R return object ----

    // Create result list with 7 components
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 7));

    // 1. coefficients matrix (n_vertices x n_cols)
    SEXP r_coeffs = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_cols));
    double* coeffs_ptr = REAL(r_coeffs);
    // Copy from Eigen matrix (column-major, same as R)
    for (int col = 0; col < n_cols; ++col) {
        for (size_t row = 0; row < n_vertices; ++row) {
            coeffs_ptr[col * n_vertices + row] = result.coefficients(row, col);
        }
    }
    SET_VECTOR_ELT(r_result, 0, r_coeffs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("coefficients"));
    UNPROTECT(1);

    // 2. gradient.neighbors (1-based, NA for extrema)
    SEXP r_neighbors = PROTECT(Rf_allocVector(INTSXP, n_vertices));
    int* neighbors_ptr = INTEGER(r_neighbors);
    for (size_t i = 0; i < n_vertices; ++i) {
        if (result.gradient_neighbors[i] == LSLOPE_R_VM_INVALID_VERTEX) {
            neighbors_ptr[i] = NA_INTEGER;
        } else {
            neighbors_ptr[i] = static_cast<int>(result.gradient_neighbors[i] + 1);
        }
    }
    SET_VECTOR_ELT(r_result, 1, r_neighbors);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("gradient.neighbors"));
    UNPROTECT(1);

    // 3. gradient.delta.y
    SEXP r_delta_y = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* delta_y_ptr = REAL(r_delta_y);
    for (size_t i = 0; i < n_vertices; ++i) {
        delta_y_ptr[i] = result.gradient_delta_y[i];
    }
    SET_VECTOR_ELT(r_result, 2, r_delta_y);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("gradient.delta.y"));
    UNPROTECT(1);

    // 4. is.local.extremum
    SEXP r_extremum = PROTECT(Rf_allocVector(LGLSXP, n_vertices));
    int* extremum_ptr = LOGICAL(r_extremum);
    for (size_t i = 0; i < n_vertices; ++i) {
        extremum_ptr[i] = result.is_local_extremum[i] ? TRUE : FALSE;
    }
    SET_VECTOR_ELT(r_result, 3, r_extremum);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("is.local.extremum"));
    UNPROTECT(1);

    // 5. n.local.maxima
    SET_VECTOR_ELT(r_result, 4, Rf_ScalarInteger(static_cast<int>(result.n_local_maxima)));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("n.local.maxima"));

    // 6. n.local.minima
    SET_VECTOR_ELT(r_result, 5, Rf_ScalarInteger(static_cast<int>(result.n_local_minima)));
    SET_STRING_ELT(r_names, 5, Rf_mkChar("n.local.minima"));

    // 7. sigmoid.alpha
    SET_VECTOR_ELT(r_result, 6, Rf_ScalarReal(result.sigmoid_alpha));
    SET_STRING_ELT(r_names, 6, Rf_mkChar("sigmoid.alpha"));

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(2);

    return r_result;
}
