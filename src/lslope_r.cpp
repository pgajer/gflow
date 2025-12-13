#include "set_wgraph.hpp"
#include "lslope.hpp"
#include "lcor.hpp"

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
static const size_t LSLOPE_R_INVALID_VERTEX = std::numeric_limits<size_t>::max();

/**
 * @brief SEXP interface for computing gradient-restricted local slope (instrumented)
 *
 * @details
 * This function provides an R interface to the gradient-restricted local slope
 * computation with full diagnostic output. It computes asymmetric association
 * measures between a directing function y and response function z along the
 * gradient direction of y.
 *
 * Three measures are supported:
 * - "slope": Raw ratio Δz/Δy along gradient edge
 * - "normalized": Sigmoid-normalized ratio tanh(α·Δz/Δy)
 * - "sign": Sign of Δz along gradient direction
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector of directing function values
 * @param s_z R numeric vector of response function values
 * @param s_type R character: "slope", "normalized", or "sign"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
 * @param s_sigmoid_alpha R numeric: sigmoid scale (0 = auto-calibrate)
 * @param s_sigmoid_type R character: "tanh", "arctan", or "algebraic"
 * @param s_ascending R logical: use ascending (TRUE) or descending (FALSE) gradient
 *
 * @return R list with components:
 *   - vertex.coefficients: Local slope at each vertex
 *   - gradient.neighbors: Gradient edge neighbor (NA if extremum)
 *   - gradient.delta.y: Δy along gradient edge
 *   - gradient.delta.z: Δz along gradient edge
 *   - is.local.extremum: Logical vector
 *   - n.local.maxima: Count of local maxima
 *   - n.local.minima: Count of local minima
 *   - sigmoid.alpha: Sigmoid scale used (for normalized type)
 *   - mean.coefficient: Mean of non-extremum coefficients
 *   - median.coefficient: Median of non-extremum coefficients
 */
extern "C" SEXP S_lslope_gradient_instrumented(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_sigmoid_alpha,
    SEXP s_sigmoid_type,
    SEXP s_ascending
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

    // Convert z vector
    if (!Rf_isReal(s_z)) {
        Rf_error("s_z must be a numeric vector");
    }
    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    // Validate dimensions
    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices) {
        Rf_error("Length of y (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_y), static_cast<int>(n_vertices));
    }
    if (n_z != n_vertices) {
        Rf_error("Length of z (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_z), static_cast<int>(n_vertices));
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

    edge_diff_type_t y_diff_type;
    if (y_diff_str == "difference") {
        y_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (y_diff_str == "logratio") {
        y_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid y_diff_type '%s'. Must be 'difference' or 'logratio'",
                 y_diff_cstr);
    }

    // Convert z_diff_type
    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    const char* z_diff_cstr = CHAR(STRING_ELT(s_z_diff_type, 0));
    std::string z_diff_str(z_diff_cstr);

    edge_diff_type_t z_diff_type;
    if (z_diff_str == "difference") {
        z_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (z_diff_str == "logratio") {
        z_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid z_diff_type '%s'. Must be 'difference' or 'logratio'",
                 z_diff_cstr);
    }

    // Convert epsilon
    if (!Rf_isReal(s_epsilon) || LENGTH(s_epsilon) != 1) {
        Rf_error("s_epsilon must be a single numeric value");
    }
    double epsilon = REAL(s_epsilon)[0];

    // Convert sigmoid_alpha
    if (!Rf_isReal(s_sigmoid_alpha) || LENGTH(s_sigmoid_alpha) != 1) {
        Rf_error("s_sigmoid_alpha must be a single numeric value");
    }
    double sigmoid_alpha = REAL(s_sigmoid_alpha)[0];

    // Convert sigmoid_type
    if (!Rf_isString(s_sigmoid_type) || LENGTH(s_sigmoid_type) != 1) {
        Rf_error("s_sigmoid_type must be a single character string");
    }
    const char* sigmoid_type_cstr = CHAR(STRING_ELT(s_sigmoid_type, 0));
    std::string sigmoid_type_str(sigmoid_type_cstr);

    sigmoid_type_t sigmoid_type;
    if (sigmoid_type_str == "tanh") {
        sigmoid_type = sigmoid_type_t::TANH;
    } else if (sigmoid_type_str == "arctan") {
        sigmoid_type = sigmoid_type_t::ARCTAN;
    } else if (sigmoid_type_str == "algebraic") {
        sigmoid_type = sigmoid_type_t::ALGEBRAIC;
    } else {
        Rf_error("Invalid sigmoid_type '%s'. Must be 'tanh', 'arctan', or 'algebraic'",
                 sigmoid_type_cstr);
    }

    // Convert ascending
    if (!Rf_isLogical(s_ascending) || LENGTH(s_ascending) != 1) {
        Rf_error("s_ascending must be a single logical value");
    }
    bool ascending = (LOGICAL(s_ascending)[0] == TRUE);

    // ---- Build graph ----
    set_wgraph_t graph(adj_list, weight_list);

    // ---- Compute local slope ----
    lslope_result_t result = graph.lslope_gradient(
        y,
        z,
        slope_type,
        y_diff_type,
        z_diff_type,
        epsilon,
        sigmoid_alpha,
        sigmoid_type,
        ascending
    );

    // ---- Build R return object ----

    // Create result list with 11 components
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 11));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 11));

    // 1. vertex.coefficients
    SEXP r_coeffs = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* coeffs_ptr = REAL(r_coeffs);
    for (size_t i = 0; i < n_vertices; ++i) {
        coeffs_ptr[i] = result.vertex_coefficients[i];
    }
    SET_VECTOR_ELT(r_result, 0, r_coeffs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("vertex.coefficients"));

    // 2. gradient.neighbors (1-based, NA for extrema)
    SEXP r_neighbors = PROTECT(Rf_allocVector(INTSXP, n_vertices));
    int* neighbors_ptr = INTEGER(r_neighbors);
    for (size_t i = 0; i < n_vertices; ++i) {
        if (result.gradient_neighbors[i] == LSLOPE_R_INVALID_VERTEX) {
            neighbors_ptr[i] = NA_INTEGER;
        } else {
            neighbors_ptr[i] = static_cast<int>(result.gradient_neighbors[i] + 1);
        }
    }
    SET_VECTOR_ELT(r_result, 1, r_neighbors);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("gradient.neighbors"));

    // 3. gradient.delta.y
    SEXP r_delta_y = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* delta_y_ptr = REAL(r_delta_y);
    for (size_t i = 0; i < n_vertices; ++i) {
        delta_y_ptr[i] = result.gradient_delta_y[i];
    }
    SET_VECTOR_ELT(r_result, 2, r_delta_y);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("gradient.delta.y"));

    // 4. gradient.delta.z
    SEXP r_delta_z = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* delta_z_ptr = REAL(r_delta_z);
    for (size_t i = 0; i < n_vertices; ++i) {
        delta_z_ptr[i] = result.gradient_delta_z[i];
    }
    SET_VECTOR_ELT(r_result, 3, r_delta_z);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("gradient.delta.z"));

    // 5. is.local.extremum
    SEXP r_extremum = PROTECT(Rf_allocVector(LGLSXP, n_vertices));
    int* extremum_ptr = LOGICAL(r_extremum);
    for (size_t i = 0; i < n_vertices; ++i) {
        extremum_ptr[i] = result.is_local_extremum[i] ? TRUE : FALSE;
    }
    SET_VECTOR_ELT(r_result, 4, r_extremum);
    SET_STRING_ELT(r_names, 4, Rf_mkChar("is.local.extremum"));

    // 6. n.local.maxima
    SEXP r_n_max = PROTECT(Rf_ScalarInteger(static_cast<int>(result.n_local_maxima)));
    SET_VECTOR_ELT(r_result, 5, r_n_max);
    SET_STRING_ELT(r_names, 5, Rf_mkChar("n.local.maxima"));

    // 7. n.local.minima
    SEXP r_n_min = PROTECT(Rf_ScalarInteger(static_cast<int>(result.n_local_minima)));
    SET_VECTOR_ELT(r_result, 6, r_n_min);
    SET_STRING_ELT(r_names, 6, Rf_mkChar("n.local.minima"));

    // 8. sigmoid.alpha
    SEXP r_alpha = PROTECT(Rf_ScalarReal(result.sigmoid_alpha));
    SET_VECTOR_ELT(r_result, 7, r_alpha);
    SET_STRING_ELT(r_names, 7, Rf_mkChar("sigmoid.alpha"));

    // 9. mean.coefficient
    SEXP r_mean = PROTECT(Rf_ScalarReal(result.mean_coefficient));
    SET_VECTOR_ELT(r_result, 8, r_mean);
    SET_STRING_ELT(r_names, 8, Rf_mkChar("mean.coefficient"));

    // 10. median.coefficient
    SEXP r_median = PROTECT(Rf_ScalarReal(result.median_coefficient));
    SET_VECTOR_ELT(r_result, 9, r_median);
    SET_STRING_ELT(r_names, 9, Rf_mkChar("median.coefficient"));

    // 11. n.vertices (for convenience)
    SEXP r_n_vert = PROTECT(Rf_ScalarInteger(static_cast<int>(n_vertices)));
    SET_VECTOR_ELT(r_result, 10, r_n_vert);
    SET_STRING_ELT(r_names, 10, Rf_mkChar("n.vertices"));

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(13);  // r_result, r_names, + 11 components

    return r_result;
}


/**
 * @brief SEXP interface for computing gradient-restricted local slope (production)
 *
 * Streamlined version that returns only the coefficient vector.
 */
extern "C" SEXP S_lslope_gradient(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_sigmoid_alpha,
    SEXP s_sigmoid_type,
    SEXP s_ascending
) {
    // ---- Input validation and conversion ----

    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    if (!Rf_isReal(s_y)) Rf_error("s_y must be a numeric vector");
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    if (!Rf_isReal(s_z)) Rf_error("s_z must be a numeric vector");
    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices || n_z != n_vertices) {
        Rf_error("Length mismatch between y, z, and number of vertices");
    }

    // Parse type
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    std::string type_str(CHAR(STRING_ELT(s_type, 0)));

    lslope_type_t slope_type;
    if (type_str == "slope") {
        slope_type = lslope_type_t::GRADIENT_SLOPE;
    } else if (type_str == "normalized") {
        slope_type = lslope_type_t::GRADIENT_SLOPE_NORMALIZED;
    } else if (type_str == "sign") {
        slope_type = lslope_type_t::GRADIENT_SIGN;
    } else {
        Rf_error("Invalid type. Must be 'slope', 'normalized', or 'sign'");
    }

    // Parse diff types
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    std::string y_diff_str(CHAR(STRING_ELT(s_y_diff_type, 0)));
    edge_diff_type_t y_diff_type = (y_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    std::string z_diff_str(CHAR(STRING_ELT(s_z_diff_type, 0)));
    edge_diff_type_t z_diff_type = (z_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    // Parse scalars
    double epsilon = Rf_isReal(s_epsilon) ? REAL(s_epsilon)[0] : 0.0;
    double sigmoid_alpha = Rf_isReal(s_sigmoid_alpha) ? REAL(s_sigmoid_alpha)[0] : 0.0;

    // Convert sigmoid_type
    if (!Rf_isString(s_sigmoid_type) || LENGTH(s_sigmoid_type) != 1) {
        Rf_error("s_sigmoid_type must be a single character string");
    }
    const char* sigmoid_type_cstr = CHAR(STRING_ELT(s_sigmoid_type, 0));
    std::string sigmoid_type_str(sigmoid_type_cstr);

    sigmoid_type_t sigmoid_type;
    if (sigmoid_type_str == "tanh") {
        sigmoid_type = sigmoid_type_t::TANH;
    } else if (sigmoid_type_str == "arctan") {
        sigmoid_type = sigmoid_type_t::ARCTAN;
    } else if (sigmoid_type_str == "algebraic") {
        sigmoid_type = sigmoid_type_t::ALGEBRAIC;
    } else {
        Rf_error("Invalid sigmoid_type '%s'. Must be 'tanh', 'arctan', or 'algebraic'",
                 sigmoid_type_cstr);
    }

    // Convert ascending
    if (!Rf_isLogical(s_ascending) || LENGTH(s_ascending) != 1) {
        Rf_error("s_ascending must be a single logical value");
    }

    bool ascending = Rf_isLogical(s_ascending) ? (LOGICAL(s_ascending)[0] == TRUE) : true;

    // Build graph and compute
    set_wgraph_t graph(adj_list, weight_list);

    lslope_result_t result = graph.lslope_gradient(
        y,
        z,
        slope_type,
        y_diff_type,
        z_diff_type,
        epsilon,
        sigmoid_alpha,
        sigmoid_type,
        ascending
        );

    // Return vector
    SEXP r_coeffs = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* coeffs_ptr = REAL(r_coeffs);
    for (size_t i = 0; i < n_vertices; ++i) {
        coeffs_ptr[i] = result.vertex_coefficients[i];
    }
    UNPROTECT(1);

    return r_coeffs;
}


/**
 * @brief SEXP interface for computing neighborhood local regression coefficient
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector of directing function values
 * @param s_z R numeric vector of response function values
 * @param s_weight_type R character: "unit" or "derivative"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
 *
 * @return R list with components:
 *   - vertex.coefficients: Local regression coefficient at each vertex
 *   - sd.y: Local standard deviation of y
 *   - sd.z: Local standard deviation of z
 *   - lcor: Local correlation (for reference)
 *   - mean.coefficient: Mean coefficient
 *   - median.coefficient: Median coefficient
 */
extern "C" SEXP S_lslope_neighborhood(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_weight_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon
) {
    // ---- Input validation and conversion ----

    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    if (!Rf_isReal(s_y)) Rf_error("s_y must be a numeric vector");
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    if (!Rf_isReal(s_z)) Rf_error("s_z must be a numeric vector");
    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices || n_z != n_vertices) {
        Rf_error("Length mismatch between y, z, and number of vertices");
    }

    // Parse weight_type
    if (!Rf_isString(s_weight_type) || LENGTH(s_weight_type) != 1) {
        Rf_error("s_weight_type must be a single character string");
    }
    std::string wt_str(CHAR(STRING_ELT(s_weight_type, 0)));
    lcor_type_t weight_type = (wt_str == "derivative")
        ? lcor_type_t::DERIVATIVE : lcor_type_t::UNIT;

    // Parse diff types
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    std::string y_diff_str(CHAR(STRING_ELT(s_y_diff_type, 0)));
    edge_diff_type_t y_diff_type = (y_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    std::string z_diff_str(CHAR(STRING_ELT(s_z_diff_type, 0)));
    edge_diff_type_t z_diff_type = (z_diff_str == "logratio")
        ? edge_diff_type_t::LOGRATIO : edge_diff_type_t::DIFFERENCE;

    double epsilon = Rf_isReal(s_epsilon) ? REAL(s_epsilon)[0] : 0.0;

    // Build graph and compute
    set_wgraph_t graph(adj_list, weight_list);
    lslope_nbhd_result_t result = graph.lslope_neighborhood(
        y, z, weight_type, y_diff_type, z_diff_type, epsilon, 0.0
    );

    // ---- Build R return object ----

    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 6));

    // 1. vertex.coefficients
    SEXP r_coeffs = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    for (size_t i = 0; i < n_vertices; ++i) {
        REAL(r_coeffs)[i] = result.vertex_coefficients[i];
    }
    SET_VECTOR_ELT(r_result, 0, r_coeffs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("vertex.coefficients"));
    UNPROTECT(1);

    // 2. sd.y
    SEXP r_sd_y = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    for (size_t i = 0; i < n_vertices; ++i) {
        REAL(r_sd_y)[i] = result.sd_y[i];
    }
    SET_VECTOR_ELT(r_result, 1, r_sd_y);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("sd.y"));
    UNPROTECT(1);

    // 3. sd.z
    SEXP r_sd_z = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    for (size_t i = 0; i < n_vertices; ++i) {
        REAL(r_sd_z)[i] = result.sd_z[i];
    }
    SET_VECTOR_ELT(r_result, 2, r_sd_z);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("sd.z"));
    UNPROTECT(1);

    // 4. lcor
    SEXP r_lcor = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    for (size_t i = 0; i < n_vertices; ++i) {
        REAL(r_lcor)[i] = result.lcor[i];
    }
    SET_VECTOR_ELT(r_result, 3, r_lcor);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("lcor"));
    UNPROTECT(1);

    // 5. mean.coefficient
    SET_VECTOR_ELT(r_result, 4, Rf_ScalarReal(result.mean_coefficient));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("mean.coefficient"));

    // 6. median.coefficient
    SET_VECTOR_ELT(r_result, 5, Rf_ScalarReal(result.median_coefficient));
    SET_STRING_ELT(r_names, 5, Rf_mkChar("median.coefficient"));

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(2);

    return r_result;
}
