/**
 * @file gfc_flow_r.cpp
 * @brief R interface for trajectory-based GFC with spurious extrema tracking
 *
 * Exports ALL basins (retained + spurious) with full classification.
 */

#include "gfc_flow_r.h"
#include "gfc_flow.hpp"
#include "set_wgraph.hpp"

#include <vector>
#include <string>
#include <cstring>

#include <R.h>
#include <Rinternals.h>

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// ============================================================================
// Helper Functions
// ============================================================================

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

static double get_list_double(SEXP list, const char* name, double def) {
    SEXP elem = get_list_element(list, name);
    return (Rf_isNull(elem) || Rf_length(elem) == 0) ? def : Rf_asReal(elem);
}

static int get_list_int(SEXP list, const char* name, int def) {
    SEXP elem = get_list_element(list, name);
    return (Rf_isNull(elem) || Rf_length(elem) == 0) ? def : Rf_asInteger(elem);
}

static bool get_list_bool(SEXP list, const char* name, bool def) {
    SEXP elem = get_list_element(list, name);
    return (Rf_isNull(elem) || Rf_length(elem) == 0) ? def : (Rf_asLogical(elem) != 0);
}

static std::string get_list_string(SEXP list, const char* name, const char* def) {
    SEXP elem = get_list_element(list, name);
    return (Rf_isNull(elem) || Rf_length(elem) == 0) ? def : CHAR(STRING_ELT(elem, 0));
}

static gfc_flow_params_t parse_gfc_flow_params(SEXP s_params) {
    gfc_flow_params_t params;
    if (Rf_isNull(s_params) || !Rf_isNewList(s_params)) return params;

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

    std::string mod_str = get_list_string(s_params, "modulation", "NONE");
    params.modulation = string_to_gflow_modulation(mod_str);
    params.store_trajectories = get_list_bool(s_params, "store_trajectories", true);
    params.max_trajectory_length = static_cast<size_t>(
        get_list_int(s_params, "max_trajectory_length", 10000));

    return params;
}

// Convert extended basin to R list
static SEXP basin_extended_to_R(const basin_extended_t& basin) {
    // 10 elements: base 6 + is_spurious, filter_stage, merged_into, label
    SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 10));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 10));

    int idx = 0;

    // extremum_vertex (1-based)
    SEXP r_ext = PROTECT(Rf_ScalarInteger(static_cast<int>(basin.extremum_vertex) + 1));
    SET_VECTOR_ELT(r_basin, idx, r_ext);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("extremum.vertex"));
    UNPROTECT(1); ++idx;

    // extremum_value
    SEXP r_val = PROTECT(Rf_ScalarReal(basin.extremum_value));
    SET_VECTOR_ELT(r_basin, idx, r_val);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("extremum.value"));
    UNPROTECT(1); ++idx;

    // is_maximum
    SEXP r_max = PROTECT(Rf_ScalarLogical(basin.is_maximum ? TRUE : FALSE));
    SET_VECTOR_ELT(r_basin, idx, r_max);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("is.maximum"));
    UNPROTECT(1); ++idx;

    // vertices (1-based)
    R_xlen_t nv = static_cast<R_xlen_t>(basin.vertices.size());
    SEXP r_verts = PROTECT(Rf_allocVector(INTSXP, nv));
    int* vptr = INTEGER(r_verts);
    for (R_xlen_t i = 0; i < nv; ++i) {
        vptr[i] = static_cast<int>(basin.vertices[i]) + 1;
    }
    SET_VECTOR_ELT(r_basin, idx, r_verts);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertices"));
    UNPROTECT(1); ++idx;

    // hop_distances
    SEXP r_hop = PROTECT(Rf_allocVector(INTSXP, nv));
    int* hptr = INTEGER(r_hop);
    for (R_xlen_t i = 0; i < nv; ++i) {
        hptr[i] = basin.hop_distances[i];
    }
    SET_VECTOR_ELT(r_basin, idx, r_hop);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("hop.distances"));
    UNPROTECT(1); ++idx;

    // max_hop_distance
    SEXP r_mhop = PROTECT(Rf_ScalarInteger(basin.max_hop_distance));
    SET_VECTOR_ELT(r_basin, idx, r_mhop);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("max.hop.distance"));
    UNPROTECT(1); ++idx;

    // is_spurious
    SEXP r_spur = PROTECT(Rf_ScalarLogical(basin.is_spurious ? TRUE : FALSE));
    SET_VECTOR_ELT(r_basin, idx, r_spur);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("is.spurious"));
    UNPROTECT(1); ++idx;

    // filter_stage
    SEXP r_stage = PROTECT(Rf_mkString(filter_stage_to_string(basin.filter_stage).c_str()));
    SET_VECTOR_ELT(r_basin, idx, r_stage);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("filter.stage"));
    UNPROTECT(1); ++idx;

    // merged_into (1-based or NA)
    SEXP r_merged = PROTECT(Rf_ScalarInteger(
        basin.merged_into >= 0 ? basin.merged_into + 1 : NA_INTEGER));
    SET_VECTOR_ELT(r_basin, idx, r_merged);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("merged.into"));
    UNPROTECT(1); ++idx;

    // label
    SEXP r_label = PROTECT(Rf_mkString(basin.label.c_str()));
    SET_VECTOR_ELT(r_basin, idx, r_label);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("label"));
    UNPROTECT(1); ++idx;

    Rf_setAttrib(r_basin, R_NamesSymbol, r_names);
    UNPROTECT(2);
    return r_basin;
}

// Convert trajectory with endpoint classification to R
static SEXP trajectory_to_R(const gflow_trajectory_t& traj) {
    // 11 elements: base 7 + endpoint classification
    SEXP r_traj = PROTECT(Rf_allocVector(VECSXP, 11));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 11));

    int idx = 0;

    // vertices (1-based)
    R_xlen_t nv = static_cast<R_xlen_t>(traj.vertices.size());
    SEXP r_verts = PROTECT(Rf_allocVector(INTSXP, nv));
    int* vptr = INTEGER(r_verts);
    for (R_xlen_t i = 0; i < nv; ++i) {
        vptr[i] = static_cast<int>(traj.vertices[i]) + 1;
    }
    SET_VECTOR_ELT(r_traj, idx, r_verts);
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertices"));
    UNPROTECT(1); ++idx;

    // start_vertex (1-based)
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarInteger(static_cast<int>(traj.start_vertex) + 1));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("start.vertex"));
    ++idx;

    // end_vertex (1-based)
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarInteger(static_cast<int>(traj.end_vertex) + 1));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("end.vertex"));
    ++idx;

    // starts_at_lmin
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarLogical(traj.starts_at_lmin ? TRUE : FALSE));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("starts.at.lmin"));
    ++idx;

    // ends_at_lmax
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarLogical(traj.ends_at_lmax ? TRUE : FALSE));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("ends.at.lmax"));
    ++idx;

    // total_change
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarReal(traj.total_change));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("total.change"));
    ++idx;

    // trajectory_id (1-based)
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarInteger(traj.trajectory_id + 1));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("trajectory.id"));
    ++idx;

    // start_is_spurious
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarLogical(traj.start_is_spurious ? TRUE : FALSE));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("start.is.spurious"));
    ++idx;

    // end_is_spurious
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarLogical(traj.end_is_spurious ? TRUE : FALSE));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("end.is.spurious"));
    ++idx;

    // start_basin_idx (1-based or NA)
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarInteger(
        traj.start_basin_idx >= 0 ? traj.start_basin_idx + 1 : NA_INTEGER));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("start.basin.idx"));
    ++idx;

    // end_basin_idx (1-based or NA)
    SET_VECTOR_ELT(r_traj, idx, Rf_ScalarInteger(
        traj.end_basin_idx >= 0 ? traj.end_basin_idx + 1 : NA_INTEGER));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("end.basin.idx"));
    ++idx;

    Rf_setAttrib(r_traj, R_NamesSymbol, r_names);
    UNPROTECT(2);
    return r_traj;
}

// Convert extended summary to data frame
static SEXP summaries_extended_to_R(const std::vector<extremum_summary_extended_t>& summaries) {
    R_xlen_t n = static_cast<R_xlen_t>(summaries.size());

    // 14 columns: base 10 + is_spurious, filter_stage, merged_into, label
    SEXP r_vertex = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_value = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_rel_value = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_is_maximum = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP r_basin_size = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_hop_index = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_p_mean_nbrs = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_p_mean_hopk = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_degree = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_deg_pct = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_is_spurious = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP r_filter_stage = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP r_merged_into = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_label = PROTECT(Rf_allocVector(STRSXP, n));

    for (R_xlen_t i = 0; i < n; ++i) {
        const auto& s = summaries[i];
        INTEGER(r_vertex)[i] = static_cast<int>(s.vertex) + 1;
        REAL(r_value)[i] = s.value;
        REAL(r_rel_value)[i] = s.rel_value;
        LOGICAL(r_is_maximum)[i] = s.is_maximum ? TRUE : FALSE;
        INTEGER(r_basin_size)[i] = s.basin_size;
        INTEGER(r_hop_index)[i] = s.hop_index;
        REAL(r_p_mean_nbrs)[i] = s.p_mean_nbrs_dist;
        REAL(r_p_mean_hopk)[i] = s.p_mean_hopk_dist;
        INTEGER(r_degree)[i] = s.degree;
        REAL(r_deg_pct)[i] = s.deg_percentile;
        LOGICAL(r_is_spurious)[i] = s.is_spurious ? TRUE : FALSE;
        SET_STRING_ELT(r_filter_stage, i, 
            Rf_mkChar(filter_stage_to_string(s.filter_stage).c_str()));
        INTEGER(r_merged_into)[i] = s.merged_into >= 0 ? s.merged_into + 1 : NA_INTEGER;
        SET_STRING_ELT(r_label, i, Rf_mkChar(s.label.c_str()));
    }

    SEXP r_df = PROTECT(Rf_allocVector(VECSXP, 14));
    SET_VECTOR_ELT(r_df, 0, r_vertex);
    SET_VECTOR_ELT(r_df, 1, r_value);
    SET_VECTOR_ELT(r_df, 2, r_rel_value);
    SET_VECTOR_ELT(r_df, 3, r_is_maximum);
    SET_VECTOR_ELT(r_df, 4, r_basin_size);
    SET_VECTOR_ELT(r_df, 5, r_hop_index);
    SET_VECTOR_ELT(r_df, 6, r_p_mean_nbrs);
    SET_VECTOR_ELT(r_df, 7, r_p_mean_hopk);
    SET_VECTOR_ELT(r_df, 8, r_degree);
    SET_VECTOR_ELT(r_df, 9, r_deg_pct);
    SET_VECTOR_ELT(r_df, 10, r_is_spurious);
    SET_VECTOR_ELT(r_df, 11, r_filter_stage);
    SET_VECTOR_ELT(r_df, 12, r_merged_into);
    SET_VECTOR_ELT(r_df, 13, r_label);

    SEXP r_colnames = PROTECT(Rf_allocVector(STRSXP, 14));
    SET_STRING_ELT(r_colnames, 0, Rf_mkChar("vertex"));
    SET_STRING_ELT(r_colnames, 1, Rf_mkChar("value"));
    SET_STRING_ELT(r_colnames, 2, Rf_mkChar("rel.value"));
    SET_STRING_ELT(r_colnames, 3, Rf_mkChar("is.maximum"));
    SET_STRING_ELT(r_colnames, 4, Rf_mkChar("basin.size"));
    SET_STRING_ELT(r_colnames, 5, Rf_mkChar("hop.index"));
    SET_STRING_ELT(r_colnames, 6, Rf_mkChar("p.mean.nbrs.dist"));
    SET_STRING_ELT(r_colnames, 7, Rf_mkChar("p.mean.hopk.dist"));
    SET_STRING_ELT(r_colnames, 8, Rf_mkChar("degree"));
    SET_STRING_ELT(r_colnames, 9, Rf_mkChar("deg.percentile"));
    SET_STRING_ELT(r_colnames, 10, Rf_mkChar("is.spurious"));
    SET_STRING_ELT(r_colnames, 11, Rf_mkChar("filter.stage"));
    SET_STRING_ELT(r_colnames, 12, Rf_mkChar("merged.into"));
    SET_STRING_ELT(r_colnames, 13, Rf_mkChar("label"));
    Rf_setAttrib(r_df, R_NamesSymbol, r_colnames);

    SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(r_rownames)[0] = NA_INTEGER;
    INTEGER(r_rownames)[1] = -static_cast<int>(n);
    Rf_setAttrib(r_df, R_RowNamesSymbol, r_rownames);
    Rf_setAttrib(r_df, R_ClassSymbol, Rf_mkString("data.frame"));

    UNPROTECT(17);
    return r_df;
}

// Convert membership list to R
static SEXP membership_to_R(const std::vector<std::vector<int>>& membership, bool one_based = true) {
    R_xlen_t n = static_cast<R_xlen_t>(membership.size());
    SEXP r_mem = PROTECT(Rf_allocVector(VECSXP, n));

    for (R_xlen_t i = 0; i < n; ++i) {
        R_xlen_t nb = static_cast<R_xlen_t>(membership[i].size());
        SEXP r_basins = PROTECT(Rf_allocVector(INTSXP, nb));
        int* ptr = INTEGER(r_basins);
        for (R_xlen_t j = 0; j < nb; ++j) {
            ptr[j] = one_based ? membership[i][j] + 1 : membership[i][j];
        }
        SET_VECTOR_ELT(r_mem, i, r_basins);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return r_mem;
}

// Convert int vector to R with optional 1-based conversion
static SEXP int_vec_to_R(const std::vector<int>& vec, bool one_based = true) {
    R_xlen_t n = static_cast<R_xlen_t>(vec.size());
    SEXP r_vec = PROTECT(Rf_allocVector(INTSXP, n));
    int* ptr = INTEGER(r_vec);
    for (R_xlen_t i = 0; i < n; ++i) {
        if (vec[i] >= 0) {
            ptr[i] = one_based ? vec[i] + 1 : vec[i];
        } else {
            ptr[i] = NA_INTEGER;
        }
    }
    UNPROTECT(1);
    return r_vec;
}

// Main result conversion
static SEXP gfc_flow_result_to_R(const gfc_flow_result_t& result) {
    // 30 elements total
    const int n_elements = 32;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_elements));

    int idx = 0;

    // 0: max.basins.all
    {
        R_xlen_t nb = static_cast<R_xlen_t>(result.max_basins_all.size());
        SEXP r_basins = PROTECT(Rf_allocVector(VECSXP, nb));
        for (R_xlen_t i = 0; i < nb; ++i) {
            SEXP r_b = PROTECT(basin_extended_to_R(result.max_basins_all[i]));
            SET_VECTOR_ELT(r_basins, i, r_b);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_basins);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.basins.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 1: min.basins.all
    {
        R_xlen_t nb = static_cast<R_xlen_t>(result.min_basins_all.size());
        SEXP r_basins = PROTECT(Rf_allocVector(VECSXP, nb));
        for (R_xlen_t i = 0; i < nb; ++i) {
            SEXP r_b = PROTECT(basin_extended_to_R(result.min_basins_all[i]));
            SET_VECTOR_ELT(r_basins, i, r_b);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_basins);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.basins.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 2: max.summaries.all
    {
        SEXP r_sum = PROTECT(summaries_extended_to_R(result.max_summaries_all));
        SET_VECTOR_ELT(r_result, idx, r_sum);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.summaries.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 3: min.summaries.all
    {
        SEXP r_sum = PROTECT(summaries_extended_to_R(result.min_summaries_all));
        SET_VECTOR_ELT(r_result, idx, r_sum);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.summaries.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 4: retained.max.indices (1-based)
    {
        SEXP r_idx = PROTECT(int_vec_to_R(result.retained_max_indices, true));
        SET_VECTOR_ELT(r_result, idx, r_idx);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("retained.max.indices"));
        UNPROTECT(1);
        ++idx;
    }

    // 5: retained.min.indices (1-based)
    {
        SEXP r_idx = PROTECT(int_vec_to_R(result.retained_min_indices, true));
        SET_VECTOR_ELT(r_result, idx, r_idx);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("retained.min.indices"));
        UNPROTECT(1);
        ++idx;
    }

    // 6: spurious.max.indices (1-based)
    {
        SEXP r_idx = PROTECT(int_vec_to_R(result.spurious_max_indices, true));
        SET_VECTOR_ELT(r_result, idx, r_idx);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("spurious.max.indices"));
        UNPROTECT(1);
        ++idx;
    }

    // 7: spurious.min.indices (1-based)
    {
        SEXP r_idx = PROTECT(int_vec_to_R(result.spurious_min_indices, true));
        SET_VECTOR_ELT(r_result, idx, r_idx);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("spurious.min.indices"));
        UNPROTECT(1);
        ++idx;
    }

    // 8: max.membership.all
    {
        SEXP r_mem = PROTECT(membership_to_R(result.max_membership_all, true));
        SET_VECTOR_ELT(r_result, idx, r_mem);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.membership.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 9: min.membership.all
    {
        SEXP r_mem = PROTECT(membership_to_R(result.min_membership_all, true));
        SET_VECTOR_ELT(r_result, idx, r_mem);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.membership.all"));
        UNPROTECT(1);
        ++idx;
    }

    // 10: max.membership.retained
    {
        SEXP r_mem = PROTECT(membership_to_R(result.max_membership_retained, true));
        SET_VECTOR_ELT(r_result, idx, r_mem);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.membership.retained"));
        UNPROTECT(1);
        ++idx;
    }

    // 11: min.membership.retained
    {
        SEXP r_mem = PROTECT(membership_to_R(result.min_membership_retained, true));
        SET_VECTOR_ELT(r_result, idx, r_mem);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.membership.retained"));
        UNPROTECT(1);
        ++idx;
    }

    // 12: max.assignment
    {
        SEXP r_assign = PROTECT(int_vec_to_R(result.max_assignment, true));
        SET_VECTOR_ELT(r_result, idx, r_assign);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.assignment"));
        UNPROTECT(1);
        ++idx;
    }

    // 13: min.assignment
    {
        SEXP r_assign = PROTECT(int_vec_to_R(result.min_assignment, true));
        SET_VECTOR_ELT(r_result, idx, r_assign);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.assignment"));
        UNPROTECT(1);
        ++idx;
    }

    // 14: trajectories
    {
        R_xlen_t nt = static_cast<R_xlen_t>(result.trajectories.size());
        SEXP r_trajs = PROTECT(Rf_allocVector(VECSXP, nt));
        for (R_xlen_t i = 0; i < nt; ++i) {
            SEXP r_t = PROTECT(trajectory_to_R(result.trajectories[i]));
            SET_VECTOR_ELT(r_trajs, i, r_t);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(r_result, idx, r_trajs);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("trajectories"));
        UNPROTECT(1);
        ++idx;
    }

    // 15: n.lmin.trajectories
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_lmin_trajectories));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.lmin.trajectories"));
    ++idx;

    // 16: n.join.trajectories
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_join_trajectories));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.join.trajectories"));
    ++idx;

    // 17: n.vertices
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(static_cast<int>(result.n_vertices)));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.vertices"));
    ++idx;

    // 18: y.median
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarReal(result.y_median));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("y.median"));
    ++idx;

    // 19: modulation
    SET_VECTOR_ELT(r_result, idx, 
        Rf_mkString(gflow_modulation_to_string(result.params.modulation).c_str()));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("modulation"));
    ++idx;

    // 20: stage_history
    {
        R_xlen_t ns = static_cast<R_xlen_t>(result.stage_history.size());
        SEXP r_stages = PROTECT(Rf_allocVector(VECSXP, 5));
        SEXP r_sname = PROTECT(Rf_allocVector(STRSXP, ns));
        SEXP r_max_before = PROTECT(Rf_allocVector(INTSXP, ns));
        SEXP r_max_after = PROTECT(Rf_allocVector(INTSXP, ns));
        SEXP r_min_before = PROTECT(Rf_allocVector(INTSXP, ns));
        SEXP r_min_after = PROTECT(Rf_allocVector(INTSXP, ns));

        for (R_xlen_t i = 0; i < ns; ++i) {
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

        SEXP r_stage_names = PROTECT(Rf_allocVector(STRSXP, 5));
        SET_STRING_ELT(r_stage_names, 0, Rf_mkChar("stage"));
        SET_STRING_ELT(r_stage_names, 1, Rf_mkChar("n.max.before"));
        SET_STRING_ELT(r_stage_names, 2, Rf_mkChar("n.max.after"));
        SET_STRING_ELT(r_stage_names, 3, Rf_mkChar("n.min.before"));
        SET_STRING_ELT(r_stage_names, 4, Rf_mkChar("n.min.after"));
        Rf_setAttrib(r_stages, R_NamesSymbol, r_stage_names);

        SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(r_rownames)[0] = NA_INTEGER;
        INTEGER(r_rownames)[1] = -static_cast<int>(ns);
        Rf_setAttrib(r_stages, R_RowNamesSymbol, r_rownames);
        Rf_setAttrib(r_stages, R_ClassSymbol, Rf_mkString("data.frame"));

        SET_VECTOR_ELT(r_result, idx, r_stages);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("stage.history"));
        UNPROTECT(8);
        ++idx;
    }

    // 21: n_max_retained
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_max_retained));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.max.retained"));
    ++idx;

    // 22: n.min.retained
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_min_retained));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.min.retained"));
    ++idx;

    // 23: n.max.spurious
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_max_spurious));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.max.spurious"));
    ++idx;

    // 24: n.min.spurious
    SET_VECTOR_ELT(r_result, idx, Rf_ScalarInteger(result.n_min_spurious));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.min.spurious"));
    ++idx;

    // ========================================================================
    // 25: n.max.all (total max basins = retained + spurious)
    // ========================================================================
    SET_VECTOR_ELT(r_result, idx,
                   Rf_ScalarInteger(static_cast<int>(result.max_basins_all.size())));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.max.all"));
    ++idx;

    // ========================================================================
    // 26: n.min.all (total min basins = retained + spurious)
    // ========================================================================
    SET_VECTOR_ELT(r_result, idx,
                   Rf_ScalarInteger(static_cast<int>(result.min_basins_all.size())));
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.min.all"));
    ++idx;

    // ========================================================================
    // 27: summary.all (combined data frame of ALL extrema: min + max)
    // ========================================================================
    {
        R_xlen_t n_min = static_cast<R_xlen_t>(result.min_summaries_all.size());
        R_xlen_t n_max = static_cast<R_xlen_t>(result.max_summaries_all.size());
        R_xlen_t n_total = n_min + n_max;

        // 14 columns
        const int n_cols = 14;
        SEXP r_df = PROTECT(Rf_allocVector(VECSXP, n_cols));
        SEXP r_col_names = PROTECT(Rf_allocVector(STRSXP, n_cols));

        SEXP r_label = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_vertex = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_value = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_rel_value = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_type = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_is_spurious = PROTECT(Rf_allocVector(LGLSXP, n_total));
        SEXP r_filter_stage = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_merged_into = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_basin_size = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_hop_index = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_p_mean_nbrs_dist = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_p_mean_hopk_dist = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_degree = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_deg_percentile = PROTECT(Rf_allocVector(REALSXP, n_total));

        R_xlen_t row = 0;

        // Fill minima first
        for (R_xlen_t i = 0; i < n_min; ++i) {
            const auto& s = result.min_summaries_all[i];
            SET_STRING_ELT(r_label, row, Rf_mkChar(s.label.c_str()));
            INTEGER(r_vertex)[row] = static_cast<int>(s.vertex) + 1;
            REAL(r_value)[row] = s.value;
            REAL(r_rel_value)[row] = s.rel_value;
            SET_STRING_ELT(r_type, row, Rf_mkChar("min"));
            LOGICAL(r_is_spurious)[row] = s.is_spurious ? TRUE : FALSE;
            SET_STRING_ELT(r_filter_stage, row,
                           Rf_mkChar(filter_stage_to_string(s.filter_stage).c_str()));
            INTEGER(r_merged_into)[row] = (s.merged_into >= 0) ?
                (static_cast<int>(s.merged_into) + 1) : NA_INTEGER;
            INTEGER(r_basin_size)[row] = s.basin_size;
            INTEGER(r_hop_index)[row] = s.hop_index;
            REAL(r_p_mean_nbrs_dist)[row] = s.p_mean_nbrs_dist;
            REAL(r_p_mean_hopk_dist)[row] = s.p_mean_hopk_dist;
            INTEGER(r_degree)[row] = s.degree;
            REAL(r_deg_percentile)[row] = s.deg_percentile;
            ++row;
        }

        // Fill maxima
        for (R_xlen_t i = 0; i < n_max; ++i) {
            const auto& s = result.max_summaries_all[i];
            SET_STRING_ELT(r_label, row, Rf_mkChar(s.label.c_str()));
            INTEGER(r_vertex)[row] = static_cast<int>(s.vertex) + 1;
            REAL(r_value)[row] = s.value;
            REAL(r_rel_value)[row] = s.rel_value;
            SET_STRING_ELT(r_type, row, Rf_mkChar("max"));
            LOGICAL(r_is_spurious)[row] = s.is_spurious ? TRUE : FALSE;
            SET_STRING_ELT(r_filter_stage, row,
                           Rf_mkChar(filter_stage_to_string(s.filter_stage).c_str()));
            INTEGER(r_merged_into)[row] = (s.merged_into >= 0) ?
                (static_cast<int>(s.merged_into) + 1) : NA_INTEGER;
            INTEGER(r_basin_size)[row] = s.basin_size;
            INTEGER(r_hop_index)[row] = s.hop_index;
            REAL(r_p_mean_nbrs_dist)[row] = s.p_mean_nbrs_dist;
            REAL(r_p_mean_hopk_dist)[row] = s.p_mean_hopk_dist;
            INTEGER(r_degree)[row] = s.degree;
            REAL(r_deg_percentile)[row] = s.deg_percentile;
            ++row;
        }

        // Set columns in data frame
        SET_VECTOR_ELT(r_df, 0, r_label);
        SET_VECTOR_ELT(r_df, 1, r_vertex);
        SET_VECTOR_ELT(r_df, 2, r_value);
        SET_VECTOR_ELT(r_df, 3, r_rel_value);
        SET_VECTOR_ELT(r_df, 4, r_type);
        SET_VECTOR_ELT(r_df, 5, r_is_spurious);
        SET_VECTOR_ELT(r_df, 6, r_filter_stage);
        SET_VECTOR_ELT(r_df, 7, r_merged_into);
        SET_VECTOR_ELT(r_df, 8, r_basin_size);
        SET_VECTOR_ELT(r_df, 9, r_hop_index);
        SET_VECTOR_ELT(r_df, 10, r_p_mean_nbrs_dist);
        SET_VECTOR_ELT(r_df, 11, r_p_mean_hopk_dist);
        SET_VECTOR_ELT(r_df, 12, r_degree);
        SET_VECTOR_ELT(r_df, 13, r_deg_percentile);

        // Column names (dot.snake format)
        SET_STRING_ELT(r_col_names, 0, Rf_mkChar("label"));
        SET_STRING_ELT(r_col_names, 1, Rf_mkChar("vertex"));
        SET_STRING_ELT(r_col_names, 2, Rf_mkChar("value"));
        SET_STRING_ELT(r_col_names, 3, Rf_mkChar("rel.value"));
        SET_STRING_ELT(r_col_names, 4, Rf_mkChar("type"));
        SET_STRING_ELT(r_col_names, 5, Rf_mkChar("is.spurious"));
        SET_STRING_ELT(r_col_names, 6, Rf_mkChar("filter.stage"));
        SET_STRING_ELT(r_col_names, 7, Rf_mkChar("merged.into"));
        SET_STRING_ELT(r_col_names, 8, Rf_mkChar("basin.size"));
        SET_STRING_ELT(r_col_names, 9, Rf_mkChar("hop.index"));
        SET_STRING_ELT(r_col_names, 10, Rf_mkChar("p.mean.nbrs.dist"));
        SET_STRING_ELT(r_col_names, 11, Rf_mkChar("p.mean.hopk.dist"));
        SET_STRING_ELT(r_col_names, 12, Rf_mkChar("degree"));
        SET_STRING_ELT(r_col_names, 13, Rf_mkChar("deg.percentile"));
        Rf_setAttrib(r_df, R_NamesSymbol, r_col_names);

        // Row names (compact form)
        SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(r_rownames)[0] = NA_INTEGER;
        INTEGER(r_rownames)[1] = -static_cast<int>(n_total);
        Rf_setAttrib(r_df, R_RowNamesSymbol, r_rownames);

        // Set class
        Rf_setAttrib(r_df, R_ClassSymbol, Rf_mkString("data.frame"));

        SET_VECTOR_ELT(r_result, idx, r_df);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("summary.all"));
        UNPROTECT(17);  // r_df, r_col_names, 14 columns, r_rownames
        ++idx;
    }

    // ========================================================================
    // 28: summary (combined data frame of RETAINED extrema only)
    // ========================================================================
    {
        R_xlen_t n_min_ret = static_cast<R_xlen_t>(result.retained_min_indices.size());
        R_xlen_t n_max_ret = static_cast<R_xlen_t>(result.retained_max_indices.size());
        R_xlen_t n_total = n_min_ret + n_max_ret;

        const int n_cols = 14;
        SEXP r_df = PROTECT(Rf_allocVector(VECSXP, n_cols));
        SEXP r_col_names = PROTECT(Rf_allocVector(STRSXP, n_cols));

        SEXP r_label = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_vertex = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_value = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_rel_value = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_type = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_is_spurious = PROTECT(Rf_allocVector(LGLSXP, n_total));
        SEXP r_filter_stage = PROTECT(Rf_allocVector(STRSXP, n_total));
        SEXP r_merged_into = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_basin_size = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_hop_index = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_p_mean_nbrs_dist = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_p_mean_hopk_dist = PROTECT(Rf_allocVector(REALSXP, n_total));
        SEXP r_degree = PROTECT(Rf_allocVector(INTSXP, n_total));
        SEXP r_deg_percentile = PROTECT(Rf_allocVector(REALSXP, n_total));

        R_xlen_t row = 0;

        // Fill retained minima
        for (R_xlen_t i = 0; i < n_min_ret; ++i) {
            int summary_idx = result.retained_min_indices[i];
            const auto& s = result.min_summaries_all[summary_idx];
            SET_STRING_ELT(r_label, row, Rf_mkChar(s.label.c_str()));
            INTEGER(r_vertex)[row] = static_cast<int>(s.vertex) + 1;
            REAL(r_value)[row] = s.value;
            REAL(r_rel_value)[row] = s.rel_value;
            SET_STRING_ELT(r_type, row, Rf_mkChar("min"));
            LOGICAL(r_is_spurious)[row] = FALSE;
            SET_STRING_ELT(r_filter_stage, row, Rf_mkChar("none"));
            INTEGER(r_merged_into)[row] = NA_INTEGER;
            INTEGER(r_basin_size)[row] = s.basin_size;
            INTEGER(r_hop_index)[row] = s.hop_index;
            REAL(r_p_mean_nbrs_dist)[row] = s.p_mean_nbrs_dist;
            REAL(r_p_mean_hopk_dist)[row] = s.p_mean_hopk_dist;
            INTEGER(r_degree)[row] = s.degree;
            REAL(r_deg_percentile)[row] = s.deg_percentile;
            ++row;
        }

        // Fill retained maxima
        for (R_xlen_t i = 0; i < n_max_ret; ++i) {
            int summary_idx = result.retained_max_indices[i];
            const auto& s = result.max_summaries_all[summary_idx];
            SET_STRING_ELT(r_label, row, Rf_mkChar(s.label.c_str()));
            INTEGER(r_vertex)[row] = static_cast<int>(s.vertex) + 1;
            REAL(r_value)[row] = s.value;
            REAL(r_rel_value)[row] = s.rel_value;
            SET_STRING_ELT(r_type, row, Rf_mkChar("max"));
            LOGICAL(r_is_spurious)[row] = FALSE;
            SET_STRING_ELT(r_filter_stage, row, Rf_mkChar("none"));
            INTEGER(r_merged_into)[row] = NA_INTEGER;
            INTEGER(r_basin_size)[row] = s.basin_size;
            INTEGER(r_hop_index)[row] = s.hop_index;
            REAL(r_p_mean_nbrs_dist)[row] = s.p_mean_nbrs_dist;
            REAL(r_p_mean_hopk_dist)[row] = s.p_mean_hopk_dist;
            INTEGER(r_degree)[row] = s.degree;
            REAL(r_deg_percentile)[row] = s.deg_percentile;
            ++row;
        }

        // Set columns
        SET_VECTOR_ELT(r_df, 0, r_label);
        SET_VECTOR_ELT(r_df, 1, r_vertex);
        SET_VECTOR_ELT(r_df, 2, r_value);
        SET_VECTOR_ELT(r_df, 3, r_rel_value);
        SET_VECTOR_ELT(r_df, 4, r_type);
        SET_VECTOR_ELT(r_df, 5, r_is_spurious);
        SET_VECTOR_ELT(r_df, 6, r_filter_stage);
        SET_VECTOR_ELT(r_df, 7, r_merged_into);
        SET_VECTOR_ELT(r_df, 8, r_basin_size);
        SET_VECTOR_ELT(r_df, 9, r_hop_index);
        SET_VECTOR_ELT(r_df, 10, r_p_mean_nbrs_dist);
        SET_VECTOR_ELT(r_df, 11, r_p_mean_hopk_dist);
        SET_VECTOR_ELT(r_df, 12, r_degree);
        SET_VECTOR_ELT(r_df, 13, r_deg_percentile);

        // Column names (dot.snake format)
        SET_STRING_ELT(r_col_names, 0, Rf_mkChar("label"));
        SET_STRING_ELT(r_col_names, 1, Rf_mkChar("vertex"));
        SET_STRING_ELT(r_col_names, 2, Rf_mkChar("value"));
        SET_STRING_ELT(r_col_names, 3, Rf_mkChar("rel.value"));
        SET_STRING_ELT(r_col_names, 4, Rf_mkChar("type"));
        SET_STRING_ELT(r_col_names, 5, Rf_mkChar("is.spurious"));
        SET_STRING_ELT(r_col_names, 6, Rf_mkChar("filter.stage"));
        SET_STRING_ELT(r_col_names, 7, Rf_mkChar("merged.into"));
        SET_STRING_ELT(r_col_names, 8, Rf_mkChar("basin.size"));
        SET_STRING_ELT(r_col_names, 9, Rf_mkChar("hop.index"));
        SET_STRING_ELT(r_col_names, 10, Rf_mkChar("p.mean.nbrs.dist"));
        SET_STRING_ELT(r_col_names, 11, Rf_mkChar("p.mean.hopk.dist"));
        SET_STRING_ELT(r_col_names, 12, Rf_mkChar("degree"));
        SET_STRING_ELT(r_col_names, 13, Rf_mkChar("deg.percentile"));
        Rf_setAttrib(r_df, R_NamesSymbol, r_col_names);

        // Row names
        SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(r_rownames)[0] = NA_INTEGER;
        INTEGER(r_rownames)[1] = -static_cast<int>(n_total);
        Rf_setAttrib(r_df, R_RowNamesSymbol, r_rownames);

        // Class
        Rf_setAttrib(r_df, R_ClassSymbol, Rf_mkString("data.frame"));

        SET_VECTOR_ELT(r_result, idx, r_df);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("summary"));
        UNPROTECT(17);
        ++idx;
    }

    // ========================================================================
    // 29: edge.length.thld
    // ========================================================================
    {
        SET_VECTOR_ELT(r_result, idx, Rf_ScalarReal(result.edge_length_thld));
        SET_STRING_ELT(r_names, idx, Rf_mkChar("edge.length.thld"));
        ++idx;
    }

    // ========================================================================
    // 30: max.overlap.dist (Eigen matrix -> R matrix)
    // ========================================================================
    {
        R_xlen_t nrow = result.max_overlap_dist.rows();
        R_xlen_t ncol = result.max_overlap_dist.cols();
        SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* mat_ptr = REAL(r_mat);

        // Copy data (Eigen is column-major like R)
        for (R_xlen_t j = 0; j < ncol; ++j) {
            for (R_xlen_t i = 0; i < nrow; ++i) {
                mat_ptr[i + j * nrow] = result.max_overlap_dist(i, j);
            }
        }

        // Add row/column names using basin labels
        if (nrow > 0) {
            SEXP r_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
            SEXP r_rownames = PROTECT(Rf_allocVector(STRSXP, nrow));
            SEXP r_colnames = PROTECT(Rf_allocVector(STRSXP, ncol));

            for (R_xlen_t i = 0; i < nrow; ++i) {
                SET_STRING_ELT(r_rownames, i,
                               Rf_mkChar(result.max_summaries_all[i].label.c_str()));
            }
            for (R_xlen_t j = 0; j < ncol; ++j) {
                SET_STRING_ELT(r_colnames, j,
                               Rf_mkChar(result.max_summaries_all[j].label.c_str()));
            }

            SET_VECTOR_ELT(r_dimnames, 0, r_rownames);
            SET_VECTOR_ELT(r_dimnames, 1, r_colnames);
            Rf_setAttrib(r_mat, R_DimNamesSymbol, r_dimnames);
            UNPROTECT(3);  // r_dimnames, r_rownames, r_colnames
        }

        SET_VECTOR_ELT(r_result, idx, r_mat);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("max.overlap.dist"));
        UNPROTECT(1);  // r_mat
        ++idx;
    }

    // ========================================================================
    // 31: min.overlap.dist (Eigen matrix -> R matrix)
    // ========================================================================
    {
        R_xlen_t nrow = result.min_overlap_dist.rows();
        R_xlen_t ncol = result.min_overlap_dist.cols();
        SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* mat_ptr = REAL(r_mat);

        // Copy data (Eigen is column-major like R)
        for (R_xlen_t j = 0; j < ncol; ++j) {
            for (R_xlen_t i = 0; i < nrow; ++i) {
                mat_ptr[i + j * nrow] = result.min_overlap_dist(i, j);
            }
        }

        // Add row/column names using basin labels
        if (nrow > 0) {
            SEXP r_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
            SEXP r_rownames = PROTECT(Rf_allocVector(STRSXP, nrow));
            SEXP r_colnames = PROTECT(Rf_allocVector(STRSXP, ncol));

            for (R_xlen_t i = 0; i < nrow; ++i) {
                SET_STRING_ELT(r_rownames, i,
                               Rf_mkChar(result.min_summaries_all[i].label.c_str()));
            }
            for (R_xlen_t j = 0; j < ncol; ++j) {
                SET_STRING_ELT(r_colnames, j,
                               Rf_mkChar(result.min_summaries_all[j].label.c_str()));
            }

            SET_VECTOR_ELT(r_dimnames, 0, r_rownames);
            SET_VECTOR_ELT(r_dimnames, 1, r_colnames);
            Rf_setAttrib(r_mat, R_DimNamesSymbol, r_dimnames);
            UNPROTECT(3);  // r_dimnames, r_rownames, r_colnames
        }

        SET_VECTOR_ELT(r_result, idx, r_mat);
        SET_STRING_ELT(r_names, idx, Rf_mkChar("min.overlap.dist"));
        UNPROTECT(1);  // r_mat
        ++idx;
    }

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    UNPROTECT(2);
    return r_result;
}

// ============================================================================
// Main SEXP Interface
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
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    set_wgraph_t graph(adj_list, weight_list);

    const R_xlen_t n = Rf_xlength(s_y);
    std::vector<double> y(n);
    double* y_ptr = REAL(s_y);
    for (R_xlen_t i = 0; i < n; ++i) {
        y[i] = y_ptr[i];
    }

    std::vector<double> density;
    if (!Rf_isNull(s_density) && Rf_xlength(s_density) > 0) {
        const R_xlen_t nd = Rf_xlength(s_density);
        density.resize(nd);
        double* d_ptr = REAL(s_density);
        for (R_xlen_t i = 0; i < nd; ++i) {
            density[i] = d_ptr[i];
        }
    }

    gfc_flow_params_t params = parse_gfc_flow_params(s_params);
    bool verbose = Rf_asLogical(s_verbose) != 0;

    gfc_flow_result_t result = compute_gfc_flow(graph, y, params, density, verbose);

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
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    set_wgraph_t graph(adj_list, weight_list);

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

    std::vector<double> density;
    if (!Rf_isNull(s_density) && Rf_xlength(s_density) > 0) {
        const R_xlen_t nd = Rf_xlength(s_density);
        density.resize(nd);
        double* d_ptr = REAL(s_density);
        for (R_xlen_t i = 0; i < nd; ++i) {
            density[i] = d_ptr[i];
        }
    }

    gfc_flow_params_t params = parse_gfc_flow_params(s_params);
    int n_cores = Rf_asInteger(s_n_cores);
    bool verbose = Rf_asLogical(s_verbose) != 0;

    std::vector<gfc_flow_result_t> results = compute_gfc_flow_matrix(
        graph, Y, params, density, n_cores, verbose);

    R_xlen_t n_results = static_cast<R_xlen_t>(results.size());
    SEXP r_results = PROTECT(Rf_allocVector(VECSXP, n_results));

    for (R_xlen_t i = 0; i < n_results; ++i) {
        SEXP r_result = PROTECT(gfc_flow_result_to_R(results[i]));
        SET_VECTOR_ELT(r_results, i, r_result);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return r_results;
}

}  // extern "C"
