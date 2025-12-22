/**
 * @file gfc_flow_r.cpp
 * @brief R interface implementation for trajectory-based GFC computation
 *
 * Implements SEXP wrapper functions that convert between R and C++ types,
 * call the C++ implementation, and convert results back to R structures.
 */

#include "gfc_flow_r.h"
#include "gfc_flow.hpp"
#include "set_wgraph.hpp"

#include <vector>
#include <string>
#include <cstring>

#include <R.h>
#include <Rinternals.h>

// Forward declarations for utility functions (defined in other _r.cpp files)
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// Helper: Get named element from R list, returns R_NilValue if not found
static SEXP get_list_element(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (Rf_isNull(names)) return R_NilValue;

    for (R_xlen_t i = 0; i < Rf_xlength(list); ++i) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

// Helper: Get double from named list element with default
static double get_list_double(SEXP list, const char* name, double default_val) {
    SEXP elem = get_list_element(list, name);
    if (Rf_isNull(elem) || Rf_length(elem) == 0) return default_val;
    return Rf_asReal(elem);
}

// Helper: Get int from named list element with default
static int get_list_int(SEXP list, const char* name, int default_val) {
    SEXP elem = get_list_element(list, name);
    if (Rf_isNull(elem) || Rf_length(elem) == 0) return default_val;
    return Rf_asInteger(elem);
}

// Helper: Get bool from named list element with default
static bool get_list_bool(SEXP list, const char* name, bool default_val) {
    SEXP elem = get_list_element(list, name);
    if (Rf_isNull(elem) || Rf_length(elem) == 0) return default_val;
    return Rf_asLogical(elem) != 0;
}

// Helper: Get string from named list element with default
static std::string get_list_string(SEXP list, const char* name, const char* default_val) {
    SEXP elem = get_list_element(list, name);
    if (Rf_isNull(elem) || Rf_length(elem) == 0) return default_val;
    return CHAR(STRING_ELT(elem, 0));
}

// Helper: Parse parameters from R list
static gfc_flow_params_t parse_gfc_flow_params(SEXP s_params) {
    gfc_flow_params_t params;

    if (Rf_isNull(s_params) || !Rf_isNewList(s_params)) {
        return params;  // Return defaults
    }

    // Base gfc_params_t fields
    params.edge_length_quantile_thld = get_list_double(s_params, "edge_length_quantile_thld", 0.9);
    params.apply_relvalue_filter = get_list_bool(s_params, "apply_relvalue_filter", true);
    params.min_rel_value_max = get_list_double(s_params, "min_rel_value_max", 1.1);
    params.max_rel_value_min = get_list_double(s_params, "max_rel_value_min", 0.9);
    params.apply_maxima_clustering = get_list_bool(s_params, "apply_maxima_clustering", true);
    params.apply_minima_clustering = get_list_bool(s_params, "apply_minima_clustering", true);
    params.max_overlap_threshold = get_list_double(s_params, "max_overlap_threshold", 0.15);
    params.min_overlap_threshold = get_list_double(s_params, "min_overlap_threshold", 0.15);
    params.apply_geometric_filter = get_list_bool(s_params, "apply_geometric_filter", true);
    params.p_mean_nbrs_dist_threshold = get_list_double(s_params, "p_mean_nbrs_dist_threshold", 0.9);
    params.p_mean_hopk_dist_threshold = get_list_double(s_params, "p_mean_hopk_dist_threshold", 0.9);
    params.p_deg_threshold = get_list_double(s_params, "p_deg_threshold", 0.9);
    params.min_basin_size = get_list_int(s_params, "min_basin_size", 10);
    params.hop_k = get_list_int(s_params, "hop_k", 2);

    // gfc_flow_params_t specific fields
    std::string mod_str = get_list_string(s_params, "modulation", "NONE");
    params.modulation = string_to_gflow_modulation(mod_str);
    params.store_trajectories = get_list_bool(s_params, "store_trajectories", true);
    params.max_trajectory_length = static_cast<size_t>(
        get_list_int(s_params, "max_trajectory_length", 10000));

    return params;
}

// Helper: Convert single basin to R list
static SEXP basin_to_R(const basin_compact_t& basin) {
    // Basin list: extremum_vertex, extremum_value, is_maximum, vertices, hop_distances, max_hop_distance
    SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 6));

    // extremum_vertex (1-based)
    SEXP r_extremum = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(r_extremum)[0] = static_cast<int>(basin.extremum_vertex) + 1;
    SET_VECTOR_ELT(r_basin, 0, r_extremum);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("extremum_vertex"));
    UNPROTECT(1);

    // extremum_value
    SEXP r_value = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_value)[0] = basin.extremum_value;
    SET_VECTOR_ELT(r_basin, 1, r_value);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("extremum_value"));
    UNPROTECT(1);

    // is_maximum
    SEXP r_is_max = PROTECT(Rf_allocVector(LGLSXP, 1));
    LOGICAL(r_is_max)[0] = basin.is_maximum ? TRUE : FALSE;
    SET_VECTOR_ELT(r_basin, 2, r_is_max);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("is_maximum"));
    UNPROTECT(1);

    // vertices (1-based)
    R_xlen_t n_verts = static_cast<R_xlen_t>(basin.vertices.size());
    SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* vptr = INTEGER(r_vertices);
    for (R_xlen_t i = 0; i < n_verts; ++i) {
        vptr[i] = static_cast<int>(basin.vertices[i]) + 1;
    }
    SET_VECTOR_ELT(r_basin, 3, r_vertices);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("vertices"));
    UNPROTECT(1);

    // hop_distances
    SEXP r_hop = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* hptr = INTEGER(r_hop);
    for (R_xlen_t i = 0; i < n_verts; ++i) {
        hptr[i] = basin.hop_distances[i];
    }
    SET_VECTOR_ELT(r_basin, 4, r_hop);
    SET_STRING_ELT(r_names, 4, Rf_mkChar("hop_distances"));
    UNPROTECT(1);

    // max_hop_distance
    SEXP r_max_hop = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(r_max_hop)[0] = basin.max_hop_distance;
    SET_VECTOR_ELT(r_basin, 5, r_max_hop);
    SET_STRING_ELT(r_names, 5, Rf_mkChar("max_hop_distance"));
    UNPROTECT(1);

    Rf_setAttrib(r_basin, R_NamesSymbol, r_names);
    UNPROTECT(2);  // r_basin, r_names

    return r_basin;
}

// Helper: Convert trajectory to R list
static SEXP trajectory_to_R(const gflow_trajectory_t& traj) {
    // Trajectory list: vertices, start_vertex, end_vertex, starts_at_lmin, ends_at_lmax, total_change, trajectory_id
    SEXP r_traj = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 7));

    // vertices (1-based)
    R_xlen_t n_verts = static_cast<R_xlen_t>(traj.vertices.size());
    SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* vptr = INTEGER(r_vertices);
    for (R_xlen_t i = 0; i < n_verts; ++i) {
        vptr[i] = static_cast<int>(traj.vertices[i]) + 1;
    }
    SET_VECTOR_ELT(r_traj, 0, r_vertices);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("vertices"));
    UNPROTECT(1);

    // start_vertex (1-based)
    SEXP r_start = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(r_start)[0] = static_cast<int>(traj.start_vertex) + 1;
    SET_VECTOR_ELT(r_traj, 1, r_start);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("start_vertex"));
    UNPROTECT(1);

    // end_vertex (1-based)
    SEXP r_end = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(r_end)[0] = static_cast<int>(traj.end_vertex) + 1;
    SET_VECTOR_ELT(r_traj, 2, r_end);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("end_vertex"));
    UNPROTECT(1);

    // starts_at_lmin
    SEXP r_starts_lmin = PROTECT(Rf_allocVector(LGLSXP, 1));
    LOGICAL(r_starts_lmin)[0] = traj.starts_at_lmin ? TRUE : FALSE;
    SET_VECTOR_ELT(r_traj, 3, r_starts_lmin);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("starts_at_lmin"));
    UNPROTECT(1);

    // ends_at_lmax
    SEXP r_ends_lmax = PROTECT(Rf_allocVector(LGLSXP, 1));
    LOGICAL(r_ends_lmax)[0] = traj.ends_at_lmax ? TRUE : FALSE;
    SET_VECTOR_ELT(r_traj, 4, r_ends_lmax);
    SET_STRING_ELT(r_names, 4, Rf_mkChar("ends_at_lmax"));
    UNPROTECT(1);

    // total_change
    SEXP r_change = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_change)[0] = traj.total_change;
    SET_VECTOR_ELT(r_traj, 5, r_change);
    SET_STRING_ELT(r_names, 5, Rf_mkChar("total_change"));
    UNPROTECT(1);

    // trajectory_id
    SEXP r_id = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(r_id)[0] = traj.trajectory_id + 1;  // 1-based
    SET_VECTOR_ELT(r_traj, 6, r_id);
    SET_STRING_ELT(r_names, 6, Rf_mkChar("trajectory_id"));
    UNPROTECT(1);

    Rf_setAttrib(r_traj, R_NamesSymbol, r_names);
    UNPROTECT(2);  // r_traj, r_names

    return r_traj;
}

// Helper: Convert gfc_flow_result_t to R list
static SEXP gfc_flow_result_to_R(const gfc_flow_result_t& result) {
    // Main result list with 12 elements
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 12));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 12));

    int idx = 0;

    // 0: max_basins
    {
        R_xlen_t n_max = static_cast<R_xlen_t>(result.max_basins.size());
        SEXP r_max_basins = PROTECT(Rf_allocVector(VECSXP, n_max));
        for (R_xlen_t i = 0; i < n_max; ++i) {
            SEXP r_basin = PROTECT(basin_to_R(result.max_basins[i]));
            SET_VECTOR_ELT(r_max_basins, i, r_basin);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_max_basins);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max_basins"));
        UNPROTECT(1);
        ++idx;
    }

    // 1: min_basins
    {
        R_xlen_t n_min = static_cast<R_xlen_t>(result.min_basins.size());
        SEXP r_min_basins = PROTECT(Rf_allocVector(VECSXP, n_min));
        for (R_xlen_t i = 0; i < n_min; ++i) {
            SEXP r_basin = PROTECT(basin_to_R(result.min_basins[i]));
            SET_VECTOR_ELT(r_min_basins, i, r_basin);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_min_basins);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min_basins"));
        UNPROTECT(1);
        ++idx;
    }

    // 2: max_assignment (1-based, NA for -1)
    {
        R_xlen_t n = static_cast<R_xlen_t>(result.max_assignment.size());
        SEXP r_max_assign = PROTECT(Rf_allocVector(INTSXP, n));
        int* ptr = INTEGER(r_max_assign);
        for (R_xlen_t i = 0; i < n; ++i) {
            ptr[i] = (result.max_assignment[i] >= 0) ? result.max_assignment[i] + 1 : NA_INTEGER;
        }
        SET_VECTOR_ELT(r_result, idx, r_max_assign);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max_assignment"));
        UNPROTECT(1);
        ++idx;
    }

    // 3: min_assignment (1-based, NA for -1)
    {
        R_xlen_t n = static_cast<R_xlen_t>(result.min_assignment.size());
        SEXP r_min_assign = PROTECT(Rf_allocVector(INTSXP, n));
        int* ptr = INTEGER(r_min_assign);
        for (R_xlen_t i = 0; i < n; ++i) {
            ptr[i] = (result.min_assignment[i] >= 0) ? result.min_assignment[i] + 1 : NA_INTEGER;
        }
        SET_VECTOR_ELT(r_result, idx, r_min_assign);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min_assignment"));
        UNPROTECT(1);
        ++idx;
    }

    // 4: trajectories
    {
        R_xlen_t n_traj = static_cast<R_xlen_t>(result.trajectories.size());
        SEXP r_trajectories = PROTECT(Rf_allocVector(VECSXP, n_traj));
        for (R_xlen_t i = 0; i < n_traj; ++i) {
            SEXP r_traj = PROTECT(trajectory_to_R(result.trajectories[i]));
            SET_VECTOR_ELT(r_trajectories, i, r_traj);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_trajectories);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("trajectories"));
        UNPROTECT(1);
        ++idx;
    }

    // 5: vertex_trajectory (1-based, NA for -1)
    {
        R_xlen_t n = static_cast<R_xlen_t>(result.vertex_trajectory.size());
        SEXP r_vtraj = PROTECT(Rf_allocVector(INTSXP, n));
        int* ptr = INTEGER(r_vtraj);
        for (R_xlen_t i = 0; i < n; ++i) {
            ptr[i] = (result.vertex_trajectory[i] >= 0) ? result.vertex_trajectory[i] + 1 : NA_INTEGER;
        }
        SET_VECTOR_ELT(r_result, idx, r_vtraj);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("vertex_trajectory"));
        UNPROTECT(1);
        ++idx;
    }

    // 6: n_lmin_trajectories
    {
        SEXP r_nlmin = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_nlmin)[0] = result.n_lmin_trajectories;
        SET_VECTOR_ELT(r_result, idx, r_nlmin);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("n_lmin_trajectories"));
        UNPROTECT(1);
        ++idx;
    }

    // 7: n_join_trajectories
    {
        SEXP r_njoin = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_njoin)[0] = result.n_join_trajectories;
        SET_VECTOR_ELT(r_result, idx, r_njoin);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("n_join_trajectories"));
        UNPROTECT(1);
        ++idx;
    }

    // 8: n_vertices
    {
        SEXP r_nv = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_nv)[0] = static_cast<int>(result.n_vertices);
        SET_VECTOR_ELT(r_result, idx, r_nv);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("n_vertices"));
        UNPROTECT(1);
        ++idx;
    }

    // 9: y_median
    {
        SEXP r_med = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(r_med)[0] = result.y_median;
        SET_VECTOR_ELT(r_result, idx, r_med);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("y_median"));
        UNPROTECT(1);
        ++idx;
    }

    // 10: stage_history (data.frame-like structure)
    {
        R_xlen_t n_stages = static_cast<R_xlen_t>(result.stage_history.size());
        SEXP r_stages = PROTECT(Rf_allocVector(VECSXP, 5));
        SEXP r_stage_names = PROTECT(Rf_allocVector(STRSXP, 5));

        SEXP r_sname = PROTECT(Rf_allocVector(STRSXP, n_stages));
        SEXP r_max_before = PROTECT(Rf_allocVector(INTSXP, n_stages));
        SEXP r_max_after = PROTECT(Rf_allocVector(INTSXP, n_stages));
        SEXP r_min_before = PROTECT(Rf_allocVector(INTSXP, n_stages));
        SEXP r_min_after = PROTECT(Rf_allocVector(INTSXP, n_stages));

        for (R_xlen_t i = 0; i < n_stages; ++i) {
            SET_STRING_ELT(r_sname, i, Rf_mkChar(result.stage_history[i].stage_name.c_str()));
            INTEGER(r_max_before)[i] = result.stage_history[i].n_max_before;
            INTEGER(r_max_after)[i] = result.stage_history[i].n_max_after;
            INTEGER(r_min_before)[i] = result.stage_history[i].n_min_before;
            INTEGER(r_min_after)[i] = result.stage_history[i].n_min_after;
        }

        SET_VECTOR_ELT(r_stages, 0, r_sname);
        SET_VECTOR_ELT(r_stages, 1, r_max_before);
        SET_VECTOR_ELT(r_stages, 2, r_max_after);
        SET_VECTOR_ELT(r_stages, 3, r_min_before);
        SET_VECTOR_ELT(r_stages, 4, r_min_after);

        SET_STRING_ELT(r_stage_names, 0, Rf_mkChar("stage"));
        SET_STRING_ELT(r_stage_names, 1, Rf_mkChar("n_max_before"));
        SET_STRING_ELT(r_stage_names, 2, Rf_mkChar("n_max_after"));
        SET_STRING_ELT(r_stage_names, 3, Rf_mkChar("n_min_before"));
        SET_STRING_ELT(r_stage_names, 4, Rf_mkChar("n_min_after"));

        Rf_setAttrib(r_stages, R_NamesSymbol, r_stage_names);

        // Make it a data.frame
        SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(r_rownames)[0] = NA_INTEGER;
        INTEGER(r_rownames)[1] = -static_cast<int>(n_stages);
        Rf_setAttrib(r_stages, R_RowNamesSymbol, r_rownames);
        Rf_setAttrib(r_stages, R_ClassSymbol, Rf_mkString("data.frame"));

        SET_VECTOR_ELT(r_result, idx, r_stages);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("stage_history"));
        UNPROTECT(8);  // r_stages, r_stage_names, r_sname, 4 int vectors, r_rownames
        ++idx;
    }

    // 11: modulation (string)
    {
        SEXP r_mod = PROTECT(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(r_mod, 0, Rf_mkChar(gflow_modulation_to_string(result.params.modulation).c_str()));
        SET_VECTOR_ELT(r_result, idx, r_mod);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("modulation"));
        UNPROTECT(1);
        ++idx;
    }

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(2);  // r_result, r_names
    return r_result;
}

// ============================================================================
// Main SEXP Interface Functions
// ============================================================================

extern "C" {

SEXP S_compute_gfc_flow(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_density,
    SEXP s_params,
    SEXP s_verbose
) {
    // Convert adjacency list and weights
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Convert y
    const R_xlen_t n = Rf_xlength(s_y);
    std::vector<double> y(n);
    double* y_ptr = REAL(s_y);
    for (R_xlen_t i = 0; i < n; ++i) {
        y[i] = y_ptr[i];
    }

    // Convert density (may be empty)
    std::vector<double> density;
    if (!Rf_isNull(s_density) && Rf_xlength(s_density) > 0) {
        const R_xlen_t nd = Rf_xlength(s_density);
        density.resize(nd);
        double* d_ptr = REAL(s_density);
        for (R_xlen_t i = 0; i < nd; ++i) {
            density[i] = d_ptr[i];
        }
    }

    // Parse parameters
    gfc_flow_params_t params = parse_gfc_flow_params(s_params);

    // Get verbose flag
    bool verbose = Rf_asLogical(s_verbose) != 0;

    // Call C++ implementation
    gfc_flow_result_t result = compute_gfc_flow(graph, y, params, density, verbose);

    // Convert result to R
    return gfc_flow_result_to_R(result);
}

SEXP S_compute_gfc_flow_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_density,
    SEXP s_params,
    SEXP s_n_cores,
    SEXP s_verbose
) {
    // Convert adjacency list and weights
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Convert Y matrix
    SEXP dims = Rf_getAttrib(s_Y, R_DimSymbol);
    int nrow = INTEGER(dims)[0];
    int ncol = INTEGER(dims)[1];

    Eigen::MatrixXd Y(nrow, ncol);
    double* Y_ptr = REAL(s_Y);
    for (int j = 0; j < ncol; ++j) {
        for (int i = 0; i < nrow; ++i) {
            Y(i, j) = Y_ptr[i + j * nrow];
        }
    }

    // Convert density
    std::vector<double> density;
    if (!Rf_isNull(s_density) && Rf_xlength(s_density) > 0) {
        const R_xlen_t nd = Rf_xlength(s_density);
        density.resize(nd);
        double* d_ptr = REAL(s_density);
        for (R_xlen_t i = 0; i < nd; ++i) {
            density[i] = d_ptr[i];
        }
    }

    // Parse parameters
    gfc_flow_params_t params = parse_gfc_flow_params(s_params);

    // Get other arguments
    int n_cores = Rf_asInteger(s_n_cores);
    bool verbose = Rf_asLogical(s_verbose) != 0;

    // Call C++ implementation
    std::vector<gfc_flow_result_t> results = compute_gfc_flow_matrix(
        graph, Y, params, density, n_cores, verbose
    );

    // Convert results to R list
    R_xlen_t n_results = static_cast<R_xlen_t>(results.size());
    SEXP r_results = PROTECT(Rf_allocVector(VECSXP, n_results));

    for (R_xlen_t i = 0; i < n_results; ++i) {
        SEXP r_result = PROTECT(gfc_flow_result_to_R(results[i]));
        SET_VECTOR_ELT(r_results, i, r_result);
        UNPROTECT(1);
    }

    UNPROTECT(1);  // r_results
    return r_results;
}

}  // extern "C"
