#include "gfassoc_r.h"
#include "gfassoc_types.hpp"
#include "gfassoc_membership.hpp"
#include "gfassoc_polarity.hpp"
#include "gradient_basin.hpp"

#include <vector>
#include <string>
#include <unordered_map>
#include <cstring>

/**
 * @file gfassoc_r.cpp
 * @brief SEXP wrapper implementations for gradient flow association functions
 *
 * These functions convert between R SEXP types and C++ types, call the
 * appropriate C++ implementations, and convert results back to SEXP.
 */

// ============================================================================
// Helper functions for SEXP conversion
// ============================================================================

/**
 * @brief Helper to get named element from R list
 */
static SEXP get_list_element(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (Rf_isNull(names)) {
        return R_NilValue;
    }

    for (int i = 0; i < Rf_length(list); ++i) {
        if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}


/**
 * @brief Convert R basin list structure to gradient_basin_t
 *
 * Expected R structure for each basin (named list):
 *   - vertex: integer (1-based, converted to 0-based)
 *   - value: numeric
 *   - hop_idx: integer
 *   - basin_df: matrix with columns (vertex, hop_distance)
 */
static gradient_basin_t convert_basin_from_R(SEXP s_basin) {
    gradient_basin_t basin;

    // Get vertex (convert from 1-based to 0-based)
    SEXP s_vertex = get_list_element(s_basin, "vertex");
    if (!Rf_isNull(s_vertex)) {
        basin.vertex = static_cast<size_t>(Rf_asInteger(s_vertex) - 1);
    } else {
        // Fallback to positional access
        basin.vertex = static_cast<size_t>(Rf_asInteger(VECTOR_ELT(s_basin, 0)) - 1);
    }

    // Get value
    SEXP s_value = get_list_element(s_basin, "value");
    if (!Rf_isNull(s_value)) {
        basin.value = Rf_asReal(s_value);
    } else {
        basin.value = Rf_asReal(VECTOR_ELT(s_basin, 1));
    }

    // Get hop_idx
    SEXP s_hop_idx = get_list_element(s_basin, "hop_idx");
    if (!Rf_isNull(s_hop_idx)) {
        basin.hop_idx = static_cast<size_t>(Rf_asInteger(s_hop_idx));
    } else {
        basin.hop_idx = static_cast<size_t>(Rf_asInteger(VECTOR_ELT(s_basin, 2)));
    }

    // Get basin_df (matrix with vertex and hop_distance columns)
    SEXP s_basin_df = get_list_element(s_basin, "basin_df");
    if (Rf_isNull(s_basin_df)) {
        s_basin_df = VECTOR_ELT(s_basin, 3);
    }

    if (!Rf_isNull(s_basin_df) && Rf_isMatrix(s_basin_df)) {
        int n_rows = Rf_nrows(s_basin_df);

        // Handle both integer and numeric matrices
        if (TYPEOF(s_basin_df) == INTSXP) {
            int* vertices = INTEGER(s_basin_df);
            int* hop_dists = vertices + n_rows;  // Second column

            for (int i = 0; i < n_rows; ++i) {
                size_t v = static_cast<size_t>(vertices[i] - 1);  // Convert to 0-based
                size_t h = static_cast<size_t>(hop_dists[i]);
                basin.hop_dist_map[v] = h;
            }
        } else if (TYPEOF(s_basin_df) == REALSXP) {
            double* vertices = REAL(s_basin_df);
            double* hop_dists = vertices + n_rows;  // Second column

            for (int i = 0; i < n_rows; ++i) {
                size_t v = static_cast<size_t>(static_cast<int>(vertices[i]) - 1);  // Convert to 0-based
                size_t h = static_cast<size_t>(static_cast<int>(hop_dists[i]));
                basin.hop_dist_map[v] = h;
            }
        }
    }

    basin.is_maximum = true;  // Will be set by caller if needed

    return basin;
}


/**
 * @brief Convert R list of basins to vector of gradient_basin_t
 */
static std::vector<gradient_basin_t> convert_basin_list_from_R(SEXP s_basins, bool is_maximum) {
    std::vector<gradient_basin_t> basins;

    if (Rf_isNull(s_basins) || Rf_length(s_basins) == 0) {
        return basins;
    }

    int n_basins = Rf_length(s_basins);
    basins.reserve(n_basins);

    for (int i = 0; i < n_basins; ++i) {
        SEXP s_basin = PROTECT(VECTOR_ELT(s_basins, i));
        gradient_basin_t basin = convert_basin_from_R(s_basin);
        basin.is_maximum = is_maximum;
        basins.push_back(basin);
        UNPROTECT(1);
    }

    return basins;
}


/**
 * @brief Convert basin_membership_t to R list
 */
static SEXP convert_membership_to_R(const basin_membership_t& membership) {
    int n_protect = 0;

    // Create output list
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 10));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
    n_protect++;

    // max_basin_indices (list of integer vectors, 0-based)
    SEXP s_max_idx = PROTECT(Rf_allocVector(VECSXP, membership.n_vertices));
    n_protect++;
    for (size_t v = 0; v < membership.n_vertices; ++v) {
        const auto& indices = membership.max_basin_indices[v];
        SEXP s_v_idx = PROTECT(Rf_allocVector(INTSXP, indices.size()));
        for (size_t k = 0; k < indices.size(); ++k) {
            INTEGER(s_v_idx)[k] = static_cast<int>(indices[k]);  // Keep 0-based
        }
        SET_VECTOR_ELT(s_max_idx, v, s_v_idx);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 0, s_max_idx);
    SET_STRING_ELT(names, 0, Rf_mkChar("max_basin_indices"));

    // min_basin_indices
    SEXP s_min_idx = PROTECT(Rf_allocVector(VECSXP, membership.n_vertices));
    n_protect++;
    for (size_t v = 0; v < membership.n_vertices; ++v) {
        const auto& indices = membership.min_basin_indices[v];
        SEXP s_v_idx = PROTECT(Rf_allocVector(INTSXP, indices.size()));
        for (size_t k = 0; k < indices.size(); ++k) {
            INTEGER(s_v_idx)[k] = static_cast<int>(indices[k]);
        }
        SET_VECTOR_ELT(s_min_idx, v, s_v_idx);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 1, s_min_idx);
    SET_STRING_ELT(names, 1, Rf_mkChar("min_basin_indices"));

    // max_membership (list of numeric vectors)
    SEXP s_max_wt = PROTECT(Rf_allocVector(VECSXP, membership.n_vertices));
    n_protect++;
    for (size_t v = 0; v < membership.n_vertices; ++v) {
        const auto& weights = membership.max_membership[v];
        SEXP s_v_wt = PROTECT(Rf_allocVector(REALSXP, weights.size()));
        std::copy(weights.begin(), weights.end(), REAL(s_v_wt));
        SET_VECTOR_ELT(s_max_wt, v, s_v_wt);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 2, s_max_wt);
    SET_STRING_ELT(names, 2, Rf_mkChar("max_membership"));

    // min_membership
    SEXP s_min_wt = PROTECT(Rf_allocVector(VECSXP, membership.n_vertices));
    n_protect++;
    for (size_t v = 0; v < membership.n_vertices; ++v) {
        const auto& weights = membership.min_membership[v];
        SEXP s_v_wt = PROTECT(Rf_allocVector(REALSXP, weights.size()));
        std::copy(weights.begin(), weights.end(), REAL(s_v_wt));
        SET_VECTOR_ELT(s_min_wt, v, s_v_wt);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 3, s_min_wt);
    SET_STRING_ELT(names, 3, Rf_mkChar("min_membership"));

    // max_vertices (1-based for R)
    SEXP s_max_v = PROTECT(Rf_allocVector(INTSXP, membership.n_max_basins));
    n_protect++;
    for (size_t i = 0; i < membership.n_max_basins; ++i) {
        INTEGER(s_max_v)[i] = static_cast<int>(membership.max_vertices[i] + 1);
    }
    SET_VECTOR_ELT(result, 4, s_max_v);
    SET_STRING_ELT(names, 4, Rf_mkChar("max_vertices"));

    // min_vertices (1-based)
    SEXP s_min_v = PROTECT(Rf_allocVector(INTSXP, membership.n_min_basins));
    n_protect++;
    for (size_t i = 0; i < membership.n_min_basins; ++i) {
        INTEGER(s_min_v)[i] = static_cast<int>(membership.min_vertices[i] + 1);
    }
    SET_VECTOR_ELT(result, 5, s_min_v);
    SET_STRING_ELT(names, 5, Rf_mkChar("min_vertices"));

    // max_values
    SEXP s_max_val = PROTECT(Rf_allocVector(REALSXP, membership.n_max_basins));
    n_protect++;
    std::copy(membership.max_values.begin(), membership.max_values.end(), REAL(s_max_val));
    SET_VECTOR_ELT(result, 6, s_max_val);
    SET_STRING_ELT(names, 6, Rf_mkChar("max_values"));

    // min_values
    SEXP s_min_val = PROTECT(Rf_allocVector(REALSXP, membership.n_min_basins));
    n_protect++;
    std::copy(membership.min_values.begin(), membership.min_values.end(), REAL(s_min_val));
    SET_VECTOR_ELT(result, 7, s_min_val);
    SET_STRING_ELT(names, 7, Rf_mkChar("min_values"));

    // n_max_basins
    SEXP s_n_max = PROTECT(Rf_ScalarInteger(static_cast<int>(membership.n_max_basins)));
    n_protect++;
    SET_VECTOR_ELT(result, 8, s_n_max);
    SET_STRING_ELT(names, 8, Rf_mkChar("n_max_basins"));

    // n_min_basins
    SEXP s_n_min = PROTECT(Rf_ScalarInteger(static_cast<int>(membership.n_min_basins)));
    n_protect++;
    SET_VECTOR_ELT(result, 9, s_n_min);
    SET_STRING_ELT(names, 9, Rf_mkChar("n_min_basins"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


/**
 * @brief Convert R membership list back to basin_membership_t
 */
static basin_membership_t convert_membership_from_R(SEXP s_membership) {
    basin_membership_t result;

    SEXP s_n_max = PROTECT(VECTOR_ELT(s_membership, 8));
    SEXP s_n_min = PROTECT(VECTOR_ELT(s_membership, 9));
    result.n_max_basins = static_cast<size_t>(INTEGER(s_n_max)[0]);
    result.n_min_basins = static_cast<size_t>(INTEGER(s_n_min)[0]);
    UNPROTECT(2);

    SEXP s_max_idx = PROTECT(VECTOR_ELT(s_membership, 0));
    result.n_vertices = static_cast<size_t>(Rf_length(s_max_idx));
    UNPROTECT(1);

    result.max_basin_indices.resize(result.n_vertices);
    result.min_basin_indices.resize(result.n_vertices);
    result.max_membership.resize(result.n_vertices);
    result.min_membership.resize(result.n_vertices);

    // max_basin_indices
    s_max_idx = PROTECT(VECTOR_ELT(s_membership, 0));
    for (size_t v = 0; v < result.n_vertices; ++v) {
        SEXP s_v = PROTECT(VECTOR_ELT(s_max_idx, v));
        int n = Rf_length(s_v);
        result.max_basin_indices[v].resize(n);
        for (int k = 0; k < n; ++k) {
            result.max_basin_indices[v][k] = static_cast<size_t>(INTEGER(s_v)[k]);
        }
        UNPROTECT(1);
    }
    UNPROTECT(1);

    // min_basin_indices
    SEXP s_min_idx = PROTECT(VECTOR_ELT(s_membership, 1));
    for (size_t v = 0; v < result.n_vertices; ++v) {
        SEXP s_v = PROTECT(VECTOR_ELT(s_min_idx, v));
        int n = Rf_length(s_v);
        result.min_basin_indices[v].resize(n);
        for (int k = 0; k < n; ++k) {
            result.min_basin_indices[v][k] = static_cast<size_t>(INTEGER(s_v)[k]);
        }
        UNPROTECT(1);
    }
    UNPROTECT(1);

    // max_membership
    SEXP s_max_wt = PROTECT(VECTOR_ELT(s_membership, 2));
    for (size_t v = 0; v < result.n_vertices; ++v) {
        SEXP s_v = PROTECT(VECTOR_ELT(s_max_wt, v));
        int n = Rf_length(s_v);
        result.max_membership[v].resize(n);
        std::copy(REAL(s_v), REAL(s_v) + n, result.max_membership[v].begin());
        UNPROTECT(1);
    }
    UNPROTECT(1);

    // min_membership
    SEXP s_min_wt = PROTECT(VECTOR_ELT(s_membership, 3));
    for (size_t v = 0; v < result.n_vertices; ++v) {
        SEXP s_v = PROTECT(VECTOR_ELT(s_min_wt, v));
        int n = Rf_length(s_v);
        result.min_membership[v].resize(n);
        std::copy(REAL(s_v), REAL(s_v) + n, result.min_membership[v].begin());
        UNPROTECT(1);
    }
    UNPROTECT(1);

    // max_vertices
    SEXP s_max_v = PROTECT(VECTOR_ELT(s_membership, 4));
    result.max_vertices.resize(result.n_max_basins);
    for (size_t i = 0; i < result.n_max_basins; ++i) {
        result.max_vertices[i] = static_cast<size_t>(INTEGER(s_max_v)[i] - 1);  // Convert to 0-based
    }
    UNPROTECT(1);

    // min_vertices
    SEXP s_min_v = PROTECT(VECTOR_ELT(s_membership, 5));
    result.min_vertices.resize(result.n_min_basins);
    for (size_t i = 0; i < result.n_min_basins; ++i) {
        result.min_vertices[i] = static_cast<size_t>(INTEGER(s_min_v)[i] - 1);
    }
    UNPROTECT(1);

    // max_values
    SEXP s_max_val = PROTECT(VECTOR_ELT(s_membership, 6));
    result.max_values.resize(result.n_max_basins);
    std::copy(REAL(s_max_val), REAL(s_max_val) + result.n_max_basins, result.max_values.begin());
    UNPROTECT(1);

    // min_values
    SEXP s_min_val = PROTECT(VECTOR_ELT(s_membership, 7));
    result.min_values.resize(result.n_min_basins);
    std::copy(REAL(s_min_val), REAL(s_min_val) + result.n_min_basins, result.min_values.begin());
    UNPROTECT(1);

    return result;
}


/**
 * @brief Convert polarity_result_t to R list
 */
static SEXP convert_polarity_to_R(const polarity_result_t& pol) {
    int n_protect = 0;
    size_t n = pol.theta.size();

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    n_protect++;

    // theta
    SEXP s_theta = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(pol.theta.begin(), pol.theta.end(), REAL(s_theta));
    SET_VECTOR_ELT(result, 0, s_theta);
    SET_STRING_ELT(names, 0, Rf_mkChar("theta"));

    // polarity
    SEXP s_polarity = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(pol.polarity.begin(), pol.polarity.end(), REAL(s_polarity));
    SET_VECTOR_ELT(result, 1, s_polarity);
    SET_STRING_ELT(names, 1, Rf_mkChar("polarity"));

    // range
    SEXP s_range = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(pol.range.begin(), pol.range.end(), REAL(s_range));
    SET_VECTOR_ELT(result, 2, s_range);
    SET_STRING_ELT(names, 2, Rf_mkChar("range"));

    // is_valid
    SEXP s_valid = PROTECT(Rf_allocVector(LGLSXP, n));
    n_protect++;
    for (size_t i = 0; i < n; ++i) {
        LOGICAL(s_valid)[i] = pol.is_valid[i] ? TRUE : FALSE;
    }
    SET_VECTOR_ELT(result, 3, s_valid);
    SET_STRING_ELT(names, 3, Rf_mkChar("is_valid"));

    // epsilon
    SEXP s_eps = PROTECT(Rf_ScalarReal(pol.epsilon));
    n_protect++;
    SET_VECTOR_ELT(result, 4, s_eps);
    SET_STRING_ELT(names, 4, Rf_mkChar("epsilon"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


/**
 * @brief Convert R polarity list back to polarity_result_t
 */
static polarity_result_t convert_polarity_from_R(SEXP s_pol) {
    polarity_result_t result;

    SEXP s_theta = PROTECT(VECTOR_ELT(s_pol, 0));
    size_t n = static_cast<size_t>(Rf_length(s_theta));
    result.theta.resize(n);
    std::copy(REAL(s_theta), REAL(s_theta) + n, result.theta.begin());
    UNPROTECT(1);

    SEXP s_polarity = PROTECT(VECTOR_ELT(s_pol, 1));
    result.polarity.resize(n);
    std::copy(REAL(s_polarity), REAL(s_polarity) + n, result.polarity.begin());
    UNPROTECT(1);

    SEXP s_range = PROTECT(VECTOR_ELT(s_pol, 2));
    result.range.resize(n);
    std::copy(REAL(s_range), REAL(s_range) + n, result.range.begin());
    UNPROTECT(1);

    SEXP s_valid = PROTECT(VECTOR_ELT(s_pol, 3));
    result.is_valid.resize(n);
    for (size_t i = 0; i < n; ++i) {
        result.is_valid[i] = (LOGICAL(s_valid)[i] == TRUE);
    }
    UNPROTECT(1);

    SEXP s_eps = PROTECT(VECTOR_ELT(s_pol, 4));
    result.epsilon = REAL(s_eps)[0];
    UNPROTECT(1);

    return result;
}


// ============================================================================
// SEXP Interface Implementations
// ============================================================================

extern "C" {

SEXP S_gfassoc_membership(
    SEXP s_lmax_basins,
    SEXP s_lmin_basins,
    SEXP s_n_vertices
) {
    // Convert inputs
    std::vector<gradient_basin_t> max_basins = convert_basin_list_from_R(s_lmax_basins, true);
    std::vector<gradient_basin_t> min_basins = convert_basin_list_from_R(s_lmin_basins, false);
    size_t n_vertices = static_cast<size_t>(INTEGER(s_n_vertices)[0]);

    // Compute membership
    basin_membership_t membership = compute_basin_membership(max_basins, min_basins, n_vertices);

    // Compute cell membership
    cell_membership_t cells = compute_cell_membership(membership);

    // Convert membership to R
    SEXP s_membership = PROTECT(convert_membership_to_R(membership));

    // Add cell membership to the result
    // For now, we'll add cell_indices as an additional element
    int n_protect = 1;

    // Extend the list to include cell information
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 12));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 12));
    n_protect++;

    // Copy existing elements
    for (int i = 0; i < 10; ++i) {
        SET_VECTOR_ELT(result, i, VECTOR_ELT(s_membership, i));
    }

    // Copy names
    SEXP old_names = PROTECT(Rf_getAttrib(s_membership, R_NamesSymbol));
    n_protect++;
    for (int i = 0; i < 10; ++i) {
        SET_STRING_ELT(names, i, STRING_ELT(old_names, i));
    }

    // Add cell_indices (list of 2-column matrices)
    SEXP s_cell_idx = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    n_protect++;
    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& cell_idx = cells.cell_indices[v];
        int n_cells = static_cast<int>(cell_idx.size());
        SEXP s_v_cells = PROTECT(Rf_allocMatrix(INTSXP, n_cells, 2));
        for (int k = 0; k < n_cells; ++k) {
            INTEGER(s_v_cells)[k] = static_cast<int>(cell_idx[k].first);          // 0-based max idx
            INTEGER(s_v_cells)[k + n_cells] = static_cast<int>(cell_idx[k].second); // 0-based min idx
        }
        SET_VECTOR_ELT(s_cell_idx, v, s_v_cells);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 10, s_cell_idx);
    SET_STRING_ELT(names, 10, Rf_mkChar("cell_indices"));

    // Add cell_membership (list of numeric vectors)
    SEXP s_cell_wt = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    n_protect++;
    for (size_t v = 0; v < n_vertices; ++v) {
        const auto& weights = cells.cell_membership[v];
        SEXP s_v_wt = PROTECT(Rf_allocVector(REALSXP, weights.size()));
        std::copy(weights.begin(), weights.end(), REAL(s_v_wt));
        SET_VECTOR_ELT(s_cell_wt, v, s_v_wt);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(result, 11, s_cell_wt);
    SET_STRING_ELT(names, 11, Rf_mkChar("cell_membership"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    // Set class
    SEXP class_name = PROTECT(Rf_mkString("gfassoc_membership"));
    n_protect++;
    Rf_setAttrib(result, R_ClassSymbol, class_name);

    UNPROTECT(n_protect);
    return result;
}


SEXP S_gfassoc_polarity(
    SEXP s_y,
    SEXP s_membership,
    SEXP s_polarity_scale,
    SEXP s_epsilon
) {
    // Convert inputs
    size_t n = static_cast<size_t>(Rf_length(s_y));
    std::vector<double> y(REAL(s_y), REAL(s_y) + n);

    basin_membership_t membership = convert_membership_from_R(s_membership);

    const char* scale_str = CHAR(STRING_ELT(s_polarity_scale, 0));
    bool use_rank = (std::strcmp(scale_str, "rank") == 0);

    double epsilon = REAL(s_epsilon)[0];

    // Compute cell membership (needed for polarity)
    cell_membership_t cells = compute_cell_membership(membership);

    // Compute polarity
    polarity_result_t pol;
    if (use_rank) {
        pol = compute_polarity_rank(y, membership, cells, epsilon);
    } else {
        pol = compute_polarity(y, membership, cells, epsilon);
    }

    // Convert to R
    SEXP result = PROTECT(convert_polarity_to_R(pol));

    // Set class
    SEXP class_name = PROTECT(Rf_mkString("gfassoc_polarity"));
    Rf_setAttrib(result, R_ClassSymbol, class_name);

    UNPROTECT(2);
    return result;
}


SEXP S_gfassoc_association(
    SEXP s_pol_y,
    SEXP s_pol_z,
    SEXP s_vertex_mass
) {
    // Convert inputs
    polarity_result_t pol_y = convert_polarity_from_R(s_pol_y);
    polarity_result_t pol_z = convert_polarity_from_R(s_pol_z);

    std::vector<double> vertex_mass;
    if (!Rf_isNull(s_vertex_mass)) {
        size_t n = static_cast<size_t>(Rf_length(s_vertex_mass));
        vertex_mass.assign(REAL(s_vertex_mass), REAL(s_vertex_mass) + n);
    }

    // Compute vertex association
    vertex_association_t va = compute_vertex_association(pol_y, pol_z);

    // Compute global association
    global_association_t ga = compute_global_association(va, vertex_mass);

    // Build result
    int n_protect = 0;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    n_protect++;

    // Vertex-level results
    size_t n = va.a_pol.size();
    SEXP s_vertex = PROTECT(Rf_allocVector(VECSXP, 4));
    n_protect++;

    SEXP s_vertex_names = PROTECT(Rf_allocVector(STRSXP, 4));
    n_protect++;

    SEXP s_a_pol = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(va.a_pol.begin(), va.a_pol.end(), REAL(s_a_pol));
    SET_VECTOR_ELT(s_vertex, 0, s_a_pol);
    SET_STRING_ELT(s_vertex_names, 0, Rf_mkChar("a_pol"));

    SEXP s_sign_pol = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(va.sign_pol.begin(), va.sign_pol.end(), REAL(s_sign_pol));
    SET_VECTOR_ELT(s_vertex, 1, s_sign_pol);
    SET_STRING_ELT(s_vertex_names, 1, Rf_mkChar("sign_pol"));

    SEXP s_conf = PROTECT(Rf_allocVector(REALSXP, n));
    n_protect++;
    std::copy(va.confidence.begin(), va.confidence.end(), REAL(s_conf));
    SET_VECTOR_ELT(s_vertex, 2, s_conf);
    SET_STRING_ELT(s_vertex_names, 2, Rf_mkChar("confidence"));

    SEXP s_valid = PROTECT(Rf_allocVector(LGLSXP, n));
    n_protect++;
    for (size_t i = 0; i < n; ++i) {
        LOGICAL(s_valid)[i] = va.is_valid[i] ? TRUE : FALSE;
    }
    SET_VECTOR_ELT(s_vertex, 3, s_valid);
    SET_STRING_ELT(s_vertex_names, 3, Rf_mkChar("is_valid"));

    Rf_setAttrib(s_vertex, R_NamesSymbol, s_vertex_names);
    SET_VECTOR_ELT(result, 0, s_vertex);
    SET_STRING_ELT(names, 0, Rf_mkChar("vertex"));

    // Global results
    SEXP s_global = PROTECT(Rf_allocVector(VECSXP, 6));
    n_protect++;

    SEXP s_global_names = PROTECT(Rf_allocVector(STRSXP, 6));
    n_protect++;

    SET_VECTOR_ELT(s_global, 0, Rf_ScalarReal(ga.A_pol));
    SET_STRING_ELT(s_global_names, 0, Rf_mkChar("A_pol"));

    SET_VECTOR_ELT(s_global, 1, Rf_ScalarReal(ga.kappa_pol));
    SET_STRING_ELT(s_global_names, 1, Rf_mkChar("kappa_pol"));

    SET_VECTOR_ELT(s_global, 2, Rf_ScalarInteger(static_cast<int>(ga.n_positive)));
    SET_STRING_ELT(s_global_names, 2, Rf_mkChar("n_positive"));

    SET_VECTOR_ELT(s_global, 3, Rf_ScalarInteger(static_cast<int>(ga.n_negative)));
    SET_STRING_ELT(s_global_names, 3, Rf_mkChar("n_negative"));

    SET_VECTOR_ELT(s_global, 4, Rf_ScalarInteger(static_cast<int>(ga.n_zero)));
    SET_STRING_ELT(s_global_names, 4, Rf_mkChar("n_zero"));

    SET_VECTOR_ELT(s_global, 5, Rf_ScalarInteger(static_cast<int>(ga.n_invalid)));
    SET_STRING_ELT(s_global_names, 5, Rf_mkChar("n_invalid"));

    Rf_setAttrib(s_global, R_NamesSymbol, s_global_names);
    SET_VECTOR_ELT(result, 1, s_global);
    SET_STRING_ELT(names, 1, Rf_mkChar("global"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


SEXP S_gfassoc_basin_character(
    SEXP s_y_membership,
    SEXP s_z_membership,
    SEXP s_pol_y,
    SEXP s_pol_z,
    SEXP s_vertex_mass
) {
    // Convert inputs
    basin_membership_t y_membership = convert_membership_from_R(s_y_membership);
    basin_membership_t z_membership = convert_membership_from_R(s_z_membership);
    polarity_result_t pol_y = convert_polarity_from_R(s_pol_y);
    polarity_result_t pol_z = convert_polarity_from_R(s_pol_z);

    std::vector<double> vertex_mass;
    if (!Rf_isNull(s_vertex_mass)) {
        size_t n = static_cast<size_t>(Rf_length(s_vertex_mass));
        vertex_mass.assign(REAL(s_vertex_mass), REAL(s_vertex_mass) + n);
    }

    // Compute basin character
    basin_character_t bc = compute_basin_character(
        y_membership, z_membership, pol_y, pol_z, vertex_mass
    );

    // Build result
    int n_protect = 0;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 8));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 8));
    n_protect++;

    // chi_y_max
    SEXP s_chi_y_max = PROTECT(Rf_allocVector(REALSXP, bc.chi_y_max.size()));
    n_protect++;
    std::copy(bc.chi_y_max.begin(), bc.chi_y_max.end(), REAL(s_chi_y_max));
    SET_VECTOR_ELT(result, 0, s_chi_y_max);
    SET_STRING_ELT(names, 0, Rf_mkChar("chi_y_max"));

    // chi_y_min
    SEXP s_chi_y_min = PROTECT(Rf_allocVector(REALSXP, bc.chi_y_min.size()));
    n_protect++;
    std::copy(bc.chi_y_min.begin(), bc.chi_y_min.end(), REAL(s_chi_y_min));
    SET_VECTOR_ELT(result, 1, s_chi_y_min);
    SET_STRING_ELT(names, 1, Rf_mkChar("chi_y_min"));

    // chi_z_max
    SEXP s_chi_z_max = PROTECT(Rf_allocVector(REALSXP, bc.chi_z_max.size()));
    n_protect++;
    std::copy(bc.chi_z_max.begin(), bc.chi_z_max.end(), REAL(s_chi_z_max));
    SET_VECTOR_ELT(result, 2, s_chi_z_max);
    SET_STRING_ELT(names, 2, Rf_mkChar("chi_z_max"));

    // chi_z_min
    SEXP s_chi_z_min = PROTECT(Rf_allocVector(REALSXP, bc.chi_z_min.size()));
    n_protect++;
    std::copy(bc.chi_z_min.begin(), bc.chi_z_min.end(), REAL(s_chi_z_min));
    SET_VECTOR_ELT(result, 3, s_chi_z_min);
    SET_STRING_ELT(names, 3, Rf_mkChar("chi_z_min"));

    // mass_y_max
    SEXP s_mass_y_max = PROTECT(Rf_allocVector(REALSXP, bc.mass_y_max.size()));
    n_protect++;
    std::copy(bc.mass_y_max.begin(), bc.mass_y_max.end(), REAL(s_mass_y_max));
    SET_VECTOR_ELT(result, 4, s_mass_y_max);
    SET_STRING_ELT(names, 4, Rf_mkChar("mass_y_max"));

    // mass_y_min
    SEXP s_mass_y_min = PROTECT(Rf_allocVector(REALSXP, bc.mass_y_min.size()));
    n_protect++;
    std::copy(bc.mass_y_min.begin(), bc.mass_y_min.end(), REAL(s_mass_y_min));
    SET_VECTOR_ELT(result, 5, s_mass_y_min);
    SET_STRING_ELT(names, 5, Rf_mkChar("mass_y_min"));

    // mass_z_max
    SEXP s_mass_z_max = PROTECT(Rf_allocVector(REALSXP, bc.mass_z_max.size()));
    n_protect++;
    std::copy(bc.mass_z_max.begin(), bc.mass_z_max.end(), REAL(s_mass_z_max));
    SET_VECTOR_ELT(result, 6, s_mass_z_max);
    SET_STRING_ELT(names, 6, Rf_mkChar("mass_z_max"));

    // mass_z_min
    SEXP s_mass_z_min = PROTECT(Rf_allocVector(REALSXP, bc.mass_z_min.size()));
    n_protect++;
    std::copy(bc.mass_z_min.begin(), bc.mass_z_min.end(), REAL(s_mass_z_min));
    SET_VECTOR_ELT(result, 7, s_mass_z_min);
    SET_STRING_ELT(names, 7, Rf_mkChar("mass_z_min"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


SEXP S_gfassoc_overlap(
    SEXP s_y_membership,
    SEXP s_z_membership,
    SEXP s_vertex_mass
) {
    // Convert inputs
    basin_membership_t y_membership = convert_membership_from_R(s_y_membership);
    basin_membership_t z_membership = convert_membership_from_R(s_z_membership);

    std::vector<double> vertex_mass;
    if (!Rf_isNull(s_vertex_mass)) {
        size_t n = static_cast<size_t>(Rf_length(s_vertex_mass));
        vertex_mass.assign(REAL(s_vertex_mass), REAL(s_vertex_mass) + n);
    }

    // Compute overlap
    overlap_matrices_t overlap = compute_soft_overlap(y_membership, z_membership, vertex_mass);

    // Build result
    int n_protect = 0;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    n_protect++;

    // Helper lambda to convert Eigen matrix to R matrix
    auto eigen_to_R = [&n_protect](const Eigen::MatrixXd& M) -> SEXP {
        int nrow = static_cast<int>(M.rows());
        int ncol = static_cast<int>(M.cols());
        SEXP s_M = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        n_protect++;
        // Eigen is column-major like R
        std::copy(M.data(), M.data() + nrow * ncol, REAL(s_M));
        return s_M;
    };

    SET_VECTOR_ELT(result, 0, eigen_to_R(overlap.O_pp));
    SET_STRING_ELT(names, 0, Rf_mkChar("O_pp"));

    SET_VECTOR_ELT(result, 1, eigen_to_R(overlap.O_mm));
    SET_STRING_ELT(names, 1, Rf_mkChar("O_mm"));

    SET_VECTOR_ELT(result, 2, eigen_to_R(overlap.O_pm));
    SET_STRING_ELT(names, 2, Rf_mkChar("O_pm"));

    SET_VECTOR_ELT(result, 3, eigen_to_R(overlap.O_mp));
    SET_STRING_ELT(names, 3, Rf_mkChar("O_mp"));

    SET_VECTOR_ELT(result, 4, Rf_ScalarReal(overlap.total_mass));
    SET_STRING_ELT(names, 4, Rf_mkChar("total_mass"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


SEXP S_gfassoc_deviation(
    SEXP s_overlap_matrix
) {
    // Convert input
    int nrow = Rf_nrows(s_overlap_matrix);
    int ncol = Rf_ncols(s_overlap_matrix);

    Eigen::MatrixXd O(nrow, ncol);
    std::copy(REAL(s_overlap_matrix), REAL(s_overlap_matrix) + nrow * ncol, O.data());

    // Compute deviation
    basin_deviation_t dev = compute_basin_deviation(O);

    // Build result
    int n_protect = 0;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
    n_protect++;

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    n_protect++;

    // Helper lambda
    auto eigen_to_R = [&n_protect](const Eigen::MatrixXd& M) -> SEXP {
        int nr = static_cast<int>(M.rows());
        int nc = static_cast<int>(M.cols());
        SEXP s_M = PROTECT(Rf_allocMatrix(REALSXP, nr, nc));
        n_protect++;
        std::copy(M.data(), M.data() + nr * nc, REAL(s_M));
        return s_M;
    };

    SET_VECTOR_ELT(result, 0, eigen_to_R(dev.delta));
    SET_STRING_ELT(names, 0, Rf_mkChar("delta"));

    SET_VECTOR_ELT(result, 1, eigen_to_R(dev.zeta));
    SET_STRING_ELT(names, 1, Rf_mkChar("zeta"));

    SET_VECTOR_ELT(result, 2, eigen_to_R(dev.expected));
    SET_STRING_ELT(names, 2, Rf_mkChar("expected"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(n_protect);
    return result;
}


SEXP S_gfcor(
    SEXP s_y_hat,
    SEXP s_z_hat,
    SEXP s_y_lmax_basins,
    SEXP s_y_lmin_basins,
    SEXP s_z_lmax_basins,
    SEXP s_z_lmin_basins,
    SEXP s_vertex_mass,
    SEXP s_options
) {
    // Extract options
    SEXP s_polarity_scale = PROTECT(VECTOR_ELT(s_options, 0));
    SEXP s_epsilon = PROTECT(VECTOR_ELT(s_options, 1));

    size_t n_vertices = static_cast<size_t>(Rf_length(s_y_hat));

    // Compute memberships
    SEXP s_n_vertices = PROTECT(Rf_ScalarInteger(static_cast<int>(n_vertices)));
    SEXP s_y_membership = PROTECT(S_gfassoc_membership(s_y_lmax_basins, s_y_lmin_basins, s_n_vertices));
    SEXP s_z_membership = PROTECT(S_gfassoc_membership(s_z_lmax_basins, s_z_lmin_basins, s_n_vertices));

    // Compute polarities
    SEXP s_pol_y = PROTECT(S_gfassoc_polarity(s_y_hat, s_y_membership, s_polarity_scale, s_epsilon));
    SEXP s_pol_z = PROTECT(S_gfassoc_polarity(s_z_hat, s_z_membership, s_polarity_scale, s_epsilon));

    // Compute association
    SEXP s_association = PROTECT(S_gfassoc_association(s_pol_y, s_pol_z, s_vertex_mass));

    // Compute basin character
    SEXP s_basin_char = PROTECT(S_gfassoc_basin_character(
        s_y_membership, s_z_membership, s_pol_y, s_pol_z, s_vertex_mass
    ));

    // Compute overlap
    SEXP s_overlap = PROTECT(S_gfassoc_overlap(s_y_membership, s_z_membership, s_vertex_mass));

    // Build comprehensive result
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));

    SET_VECTOR_ELT(result, 0, VECTOR_ELT(s_association, 1));  // global
    SET_STRING_ELT(names, 0, Rf_mkChar("global"));

    SET_VECTOR_ELT(result, 1, VECTOR_ELT(s_association, 0));  // vertex
    SET_STRING_ELT(names, 1, Rf_mkChar("vertex"));

    SET_VECTOR_ELT(result, 2, s_pol_y);
    SET_STRING_ELT(names, 2, Rf_mkChar("polarity_y"));

    SET_VECTOR_ELT(result, 3, s_pol_z);
    SET_STRING_ELT(names, 3, Rf_mkChar("polarity_z"));

    SET_VECTOR_ELT(result, 4, s_basin_char);
    SET_STRING_ELT(names, 4, Rf_mkChar("basin_character"));

    SET_VECTOR_ELT(result, 5, s_overlap);
    SET_STRING_ELT(names, 5, Rf_mkChar("overlap"));

    // Store memberships for downstream use
    SEXP s_memberships = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP s_mem_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_VECTOR_ELT(s_memberships, 0, s_y_membership);
    SET_STRING_ELT(s_mem_names, 0, Rf_mkChar("y"));
    SET_VECTOR_ELT(s_memberships, 1, s_z_membership);
    SET_STRING_ELT(s_mem_names, 1, Rf_mkChar("z"));
    Rf_setAttrib(s_memberships, R_NamesSymbol, s_mem_names);
    SET_VECTOR_ELT(result, 6, s_memberships);
    SET_STRING_ELT(names, 6, Rf_mkChar("membership"));

    Rf_setAttrib(result, R_NamesSymbol, names);

    // Set class
    SEXP class_name = PROTECT(Rf_mkString("gfcor"));
    Rf_setAttrib(result, R_ClassSymbol, class_name);

    UNPROTECT(15);
    return result;
}

} // extern "C"
