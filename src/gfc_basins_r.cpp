/**
 * @file gfc_basins_r.cpp
 * @brief SEXP interface for modulated gradient flow basin computation
 */

#include <R.h>
#include <Rinternals.h>
#include "set_wgraph.hpp"
#include "gfc.hpp"

#include <vector>
#include <unordered_map>

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

/**
 * @brief Convert bbasin_t to R list
 */
static SEXP bbasin_to_R(const bbasin_t& basin) {
    const int n_components = 5;
    SEXP s_basin = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // extremum.vertex (1-based)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("extremum.vertex"));
    SET_VECTOR_ELT(s_basin, idx++,
                   Rf_ScalarInteger(static_cast<int>(basin.extremum_vertex) + 1));

    // value
    SET_STRING_ELT(s_names, idx, Rf_mkChar("value"));
    SET_VECTOR_ELT(s_basin, idx++, Rf_ScalarReal(basin.value));

    // is.maximum
    SET_STRING_ELT(s_names, idx, Rf_mkChar("is.maximum"));
    SET_VECTOR_ELT(s_basin, idx++, Rf_ScalarLogical(basin.is_maximum ? TRUE : FALSE));

    // vertices (1-based)
    const int n_verts = static_cast<int>(basin.vertices.size());
    SEXP s_vertices = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* p_vertices = INTEGER(s_vertices);
    for (int i = 0; i < n_verts; ++i) {
        p_vertices[i] = static_cast<int>(basin.vertices[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("vertices"));
    SET_VECTOR_ELT(s_basin, idx++, s_vertices);
    UNPROTECT(1);

    // boundary (1-based)
    const int n_boundary = static_cast<int>(basin.boundary.size());
    SEXP s_boundary = PROTECT(Rf_allocVector(INTSXP, n_boundary));
    int* p_boundary = INTEGER(s_boundary);
    for (int i = 0; i < n_boundary; ++i) {
        p_boundary[i] = static_cast<int>(basin.boundary[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("boundary"));
    SET_VECTOR_ELT(s_basin, idx++, s_boundary);
    UNPROTECT(1);

    Rf_setAttrib(s_basin, R_NamesSymbol, s_names);
    UNPROTECT(2);  // s_basin, s_names

    return s_basin;
}

extern "C" SEXP S_compute_gfc_basins(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_modulation_type,
    SEXP s_density,
    SEXP s_edgelen_bandwidth,
    SEXP s_verbose
) {
    // Convert graph structure
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    // Convert y vector
    const int n = LENGTH(s_y);
    std::vector<double> y(n);
    const double* p_y = REAL(s_y);
    for (int i = 0; i < n; ++i) {
        y[i] = p_y[i];
    }

    // Parse modulation type
    const char* mod_str = CHAR(STRING_ELT(s_modulation_type, 0));
    gfc_basin_params_t params;

    if (strcmp(mod_str, "none") == 0) {
        params.modulation = gflow_modulation_t::NONE;
    } else if (strcmp(mod_str, "density") == 0) {
        params.modulation = gflow_modulation_t::DENSITY;
    } else if (strcmp(mod_str, "edgelen") == 0) {
        params.modulation = gflow_modulation_t::EDGELEN;
    } else if (strcmp(mod_str, "density_edgelen") == 0) {
        params.modulation = gflow_modulation_t::DENSITY_EDGELEN;
    } else {
        Rf_error("Unknown modulation type: %s", mod_str);
    }

    // Convert optional density vector
    if (!Rf_isNull(s_density) && LENGTH(s_density) == n) {
        const double* p_density = REAL(s_density);
        params.density.assign(p_density, p_density + n);
    }

    // Edge length bandwidth
    params.edgelen_bandwidth = Rf_asReal(s_edgelen_bandwidth);

    bool verbose = Rf_asLogical(s_verbose);

    // Compute basins
    auto basins = graph.compute_gfc_basins(y, params, verbose);

    // Convert to R list
    const int n_basins = static_cast<int>(basins.size());
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_basins));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_basins));

    int idx = 0;
    for (const auto& [vertex, basin] : basins) {
        // Create name: "min_123" or "max_123" (1-based)
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%s_%d",
                 basin.is_maximum ? "max" : "min",
                 static_cast<int>(vertex) + 1);
        SET_STRING_ELT(s_names, idx, Rf_mkChar(name_buf));

        SEXP s_basin = PROTECT(bbasin_to_R(basin));
        SET_VECTOR_ELT(s_result, idx, s_basin);
        UNPROTECT(1);

        ++idx;
    }

    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    // Add class
    SEXP s_class = PROTECT(Rf_mkString("gfc_basins"));
    Rf_setAttrib(s_result, R_ClassSymbol, s_class);

    UNPROTECT(3);  // s_result, s_names, s_class
    return s_result;
}
