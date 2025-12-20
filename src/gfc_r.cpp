/**
 * @file gfc_r.cpp
 * @brief SEXP interface implementation for Gradient Flow Complex computation
 *
 * This file implements the R-callable functions that convert between R's
 * SEXP types and C++ types, call the core GFC implementation, and convert
 * results back to R structures.
 */

#include "gfc_r.h"
#include "gfc.hpp"
#include "set_wgraph.hpp"

#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <string>

// Forward declarations for conversion utilities
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// ============================================================================
// Helper: Convert basin_compact_t to R list
// ============================================================================

static SEXP basin_compact_to_R(const basin_compact_t& basin, bool with_trajectories) {
    // Determine number of components
    const int n_components = with_trajectories ? 7 : 5;
    SEXP s_basin = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // vertex (1-based for R)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("vertex"));
    SET_VECTOR_ELT(s_basin, idx++, Rf_ScalarInteger(static_cast<int>(basin.extremum_vertex) + 1));

    // value
    SET_STRING_ELT(s_names, idx, Rf_mkChar("value"));
    SET_VECTOR_ELT(s_basin, idx++, Rf_ScalarReal(basin.extremum_value));

    // vertices (1-based for R)
    const int n_verts = static_cast<int>(basin.vertices.size());
    SEXP s_vertices = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* p_vertices = INTEGER(s_vertices);
    for (int i = 0; i < n_verts; ++i) {
        p_vertices[i] = static_cast<int>(basin.vertices[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("vertices"));
    SET_VECTOR_ELT(s_basin, idx++, s_vertices);
    UNPROTECT(1);  // s_vertices

    // hop.distances
    SEXP s_hop = PROTECT(Rf_allocVector(INTSXP, n_verts));
    int* p_hop = INTEGER(s_hop);
    for (int i = 0; i < n_verts; ++i) {
        p_hop[i] = (i < static_cast<int>(basin.hop_distances.size()))
            ? basin.hop_distances[i] : 0;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("hop.distances"));
    SET_VECTOR_ELT(s_basin, idx++, s_hop);
    UNPROTECT(1);  // s_hop

    // max.hop.distance
    SET_STRING_ELT(s_names, idx, Rf_mkChar("max.hop.distance"));
    SET_VECTOR_ELT(s_basin, idx++, Rf_ScalarInteger(basin.max_hop_distance));

    // Trajectory data (only if with_trajectories)
    if (with_trajectories) {
        // terminal.extrema (1-based)
        const int n_terminals = static_cast<int>(basin.terminal_extrema.size());
        SEXP s_terminals = PROTECT(Rf_allocVector(INTSXP, n_terminals));
        int* p_terminals = INTEGER(s_terminals);
        for (int i = 0; i < n_terminals; ++i) {
            p_terminals[i] = static_cast<int>(basin.terminal_extrema[i]) + 1;
        }
        SET_STRING_ELT(s_names, idx, Rf_mkChar("terminal.extrema"));
        SET_VECTOR_ELT(s_basin, idx++, s_terminals);
        UNPROTECT(1);

        // trajectory.sets
        const int n_traj_sets = static_cast<int>(basin.trajectory_sets.size());
        SEXP s_traj_sets = PROTECT(Rf_allocVector(VECSXP, n_traj_sets));

        for (int i = 0; i < n_traj_sets; ++i) {
            const trajectory_set_t& traj_set = basin.trajectory_sets[i];

            SEXP s_tset = PROTECT(Rf_allocVector(VECSXP, 2));
            SEXP s_tset_names = PROTECT(Rf_allocVector(STRSXP, 2));
            SET_STRING_ELT(s_tset_names, 0, Rf_mkChar("terminal.vertex"));
            SET_STRING_ELT(s_tset_names, 1, Rf_mkChar("trajectories"));
            Rf_setAttrib(s_tset, R_NamesSymbol, s_tset_names);

            // terminal.vertex (1-based)
            SEXP s_tvert = PROTECT(Rf_ScalarInteger(static_cast<int>(traj_set.terminal_vertex) + 1));
            SET_VECTOR_ELT(s_tset, 0, s_tvert);

            // trajectories (list of integer vectors, 1-based)
            const int n_trajs = static_cast<int>(traj_set.trajectories.size());
            SEXP s_trajs = PROTECT(Rf_allocVector(VECSXP, n_trajs));
            for (int j = 0; j < n_trajs; ++j) {
                const std::vector<size_t>& traj = traj_set.trajectories[j];
                const int traj_len = static_cast<int>(traj.size());
                SEXP s_traj = PROTECT(Rf_allocVector(INTSXP, traj_len));
                int* p_traj = INTEGER(s_traj);
                for (int k = 0; k < traj_len; ++k) {
                    p_traj[k] = static_cast<int>(traj[k]) + 1;
                }
                SET_VECTOR_ELT(s_trajs, j, s_traj);
                UNPROTECT(1);
            }
            SET_VECTOR_ELT(s_tset, 1, s_trajs);

            SET_VECTOR_ELT(s_traj_sets, i, s_tset);
            UNPROTECT(4);  // s_tset, s_tset_names, s_tvert, s_trajs
        }

        SET_STRING_ELT(s_names, idx, Rf_mkChar("trajectory.sets"));
        SET_VECTOR_ELT(s_basin, idx++, s_traj_sets);
        UNPROTECT(1);  // s_traj_sets
    }

    Rf_setAttrib(s_basin, R_NamesSymbol, s_names);
    UNPROTECT(2);  // s_basin, s_names

    return s_basin;
}

// ============================================================================
// Helper: Convert single joined trajectory to R list
// ============================================================================

static SEXP joined_trajectory_to_R(const joined_trajectory_t& jt) {
    const int n_components = 6;
    SEXP s_traj = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // min.vertex (1-based)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("min.vertex"));
    SET_VECTOR_ELT(s_traj, idx++, Rf_ScalarInteger(static_cast<int>(jt.min_vertex) + 1));

    // max.vertex (1-based)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("max.vertex"));
    SET_VECTOR_ELT(s_traj, idx++, Rf_ScalarInteger(static_cast<int>(jt.max_vertex) + 1));

    // path (integer vector, 1-based)
    const int path_len = static_cast<int>(jt.path.size());
    SEXP s_path = PROTECT(Rf_allocVector(INTSXP, path_len));
    int* p_path = INTEGER(s_path);
    for (int i = 0; i < path_len; ++i) {
        p_path[i] = static_cast<int>(jt.path[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("path"));
    SET_VECTOR_ELT(s_traj, idx++, s_path);
    UNPROTECT(1);

    // intermediate.extrema (integer vector, 1-based)
    const int n_inter = static_cast<int>(jt.intermediate_extrema.size());
    SEXP s_inter = PROTECT(Rf_allocVector(INTSXP, n_inter));
    int* p_inter = INTEGER(s_inter);
    for (int i = 0; i < n_inter; ++i) {
        p_inter[i] = static_cast<int>(jt.intermediate_extrema[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("intermediate.extrema"));
    SET_VECTOR_ELT(s_traj, idx++, s_inter);
    UNPROTECT(1);

    // total.change
    SET_STRING_ELT(s_names, idx, Rf_mkChar("total.change"));
    SET_VECTOR_ELT(s_traj, idx++, Rf_ScalarReal(jt.total_change));

    // path.length
    SET_STRING_ELT(s_names, idx, Rf_mkChar("path.length"));
    SET_VECTOR_ELT(s_traj, idx++, Rf_ScalarReal(jt.path_length));

    Rf_setAttrib(s_traj, R_NamesSymbol, s_names);
    UNPROTECT(2);  // s_traj, s_names

    return s_traj;
}

// ============================================================================
// Helper: Convert all joined trajectories to R list
// ============================================================================

static SEXP joined_trajectories_to_R(const std::vector<joined_trajectory_t>& trajectories) {
    const int n = static_cast<int>(trajectories.size());
    SEXP s_list = PROTECT(Rf_allocVector(VECSXP, n));

    for (int i = 0; i < n; ++i) {
        SET_VECTOR_ELT(s_list, i, joined_trajectory_to_R(trajectories[i]));
    }

    UNPROTECT(1);
    return s_list;
}

// ============================================================================
// Helper: Convert cell map to R named list
// ============================================================================

/**
 * Converts the cell map to an R named list where names are "min_vertex-max_vertex"
 * and values are integer vectors of trajectory indices (1-based).
 */
static SEXP cell_map_to_R(
    const std::map<std::pair<size_t, size_t>, std::vector<size_t>>& cell_map
) {
    const int n_cells = static_cast<int>(cell_map.size());
    SEXP s_list = PROTECT(Rf_allocVector(VECSXP, n_cells));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_cells));

    int idx = 0;
    for (const auto& [key, indices] : cell_map) {
        // Create name: "minVertex-maxVertex" (1-based)
        char name_buf[64];
        snprintf(name_buf, sizeof(name_buf), "%d-%d",
                 static_cast<int>(key.first) + 1,
                 static_cast<int>(key.second) + 1);
        SET_STRING_ELT(s_names, idx, Rf_mkChar(name_buf));

        // Create integer vector of indices (1-based)
        const int n_idx = static_cast<int>(indices.size());
        SEXP s_indices = PROTECT(Rf_allocVector(INTSXP, n_idx));
        int* p = INTEGER(s_indices);
        for (int i = 0; i < n_idx; ++i) {
            p[i] = static_cast<int>(indices[i]) + 1;
        }
        SET_VECTOR_ELT(s_list, idx, s_indices);
        UNPROTECT(1);

        ++idx;
    }

    Rf_setAttrib(s_list, R_NamesSymbol, s_names);
    UNPROTECT(2);

    return s_list;
}

// ============================================================================
// Helper: Convert extremum_summary_t to R data frame row
// ============================================================================

static SEXP summaries_to_dataframe(
    const std::vector<extremum_summary_t>& max_summaries,
    const std::vector<extremum_summary_t>& min_summaries
) {
    const int n_max = static_cast<int>(max_summaries.size());
    const int n_min = static_cast<int>(min_summaries.size());
    const int n_total = n_max + n_min;

    if (n_total == 0) {
        // Return empty data frame
        SEXP s_df = PROTECT(Rf_allocVector(VECSXP, 0));
        SEXP s_class = PROTECT(Rf_mkString("data.frame"));
        Rf_setAttrib(s_df, R_ClassSymbol, s_class);
        UNPROTECT(2);
        return s_df;
    }

    // Create columns (11 columns matching R output)
    const int n_cols = 11;
    SEXP s_df = PROTECT(Rf_allocVector(VECSXP, n_cols));
    SEXP s_colnames = PROTECT(Rf_allocVector(STRSXP, n_cols));

    // Allocate column vectors
    SEXP s_label = PROTECT(Rf_allocVector(STRSXP, n_total));
    SEXP s_vertex = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP s_value = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP s_rel_value = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP s_type = PROTECT(Rf_allocVector(STRSXP, n_total));
    SEXP s_hop_idx = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP s_basin_size = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP s_p_mean_nbrs = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP s_p_mean_hopk = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP s_deg = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP s_p_deg = PROTECT(Rf_allocVector(REALSXP, n_total));

    int* p_vertex = INTEGER(s_vertex);
    double* p_value = REAL(s_value);
    double* p_rel_value = REAL(s_rel_value);
    int* p_hop_idx = INTEGER(s_hop_idx);
    int* p_basin_size = INTEGER(s_basin_size);
    double* p_mean_nbrs = REAL(s_p_mean_nbrs);
    double* p_mean_hopk = REAL(s_p_mean_hopk);
    int* p_deg = INTEGER(s_deg);
    double* p_p_deg = REAL(s_p_deg);

    int row = 0;

    // Fill minima (m1, m2, ... in order of increasing value)
    std::vector<size_t> min_order(n_min);
    for (int i = 0; i < n_min; ++i) min_order[i] = i;
    std::sort(min_order.begin(), min_order.end(), [&](size_t a, size_t b) {
        return min_summaries[a].value < min_summaries[b].value;
    });

    for (int i = 0; i < n_min; ++i) {
        const auto& s = min_summaries[min_order[i]];
        char label[16];
        snprintf(label, sizeof(label), "m%d", i + 1);

        SET_STRING_ELT(s_label, row, Rf_mkChar(label));
        p_vertex[row] = static_cast<int>(s.vertex) + 1;  // 1-based
        p_value[row] = s.value;
        p_rel_value[row] = s.rel_value;
        SET_STRING_ELT(s_type, row, Rf_mkChar("min"));
        p_hop_idx[row] = s.hop_index;
        p_basin_size[row] = s.basin_size;
        p_mean_nbrs[row] = s.p_mean_nbrs_dist;
        p_mean_hopk[row] = s.p_mean_hopk_dist;
        p_deg[row] = s.degree;
        p_p_deg[row] = s.deg_percentile;
        ++row;
    }

    // Fill maxima (M1, M2, ... in order of decreasing value)
    std::vector<size_t> max_order(n_max);
    for (int i = 0; i < n_max; ++i) max_order[i] = i;
    std::sort(max_order.begin(), max_order.end(), [&](size_t a, size_t b) {
        return max_summaries[a].value > max_summaries[b].value;
    });

    for (int i = 0; i < n_max; ++i) {
        const auto& s = max_summaries[max_order[i]];
        char label[16];
        snprintf(label, sizeof(label), "M%d", i + 1);

        SET_STRING_ELT(s_label, row, Rf_mkChar(label));
        p_vertex[row] = static_cast<int>(s.vertex) + 1;
        p_value[row] = s.value;
        p_rel_value[row] = s.rel_value;
        SET_STRING_ELT(s_type, row, Rf_mkChar("max"));
        p_hop_idx[row] = s.hop_index;
        p_basin_size[row] = s.basin_size;
        p_mean_nbrs[row] = s.p_mean_nbrs_dist;
        p_mean_hopk[row] = s.p_mean_hopk_dist;
        p_deg[row] = s.degree;
        p_p_deg[row] = s.deg_percentile;
        ++row;
    }

    // Set columns
    int col = 0;
    SET_STRING_ELT(s_colnames, col, Rf_mkChar("label"));
    SET_VECTOR_ELT(s_df, col++, s_label);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("vertex"));
    SET_VECTOR_ELT(s_df, col++, s_vertex);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("value"));
    SET_VECTOR_ELT(s_df, col++, s_value);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("rel.value"));
    SET_VECTOR_ELT(s_df, col++, s_rel_value);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("type"));
    SET_VECTOR_ELT(s_df, col++, s_type);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("hop.idx"));
    SET_VECTOR_ELT(s_df, col++, s_hop_idx);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("basin.size"));
    SET_VECTOR_ELT(s_df, col++, s_basin_size);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("p.mean.nbrs.dist"));
    SET_VECTOR_ELT(s_df, col++, s_p_mean_nbrs);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("p.mean.hopk.dist"));
    SET_VECTOR_ELT(s_df, col++, s_p_mean_hopk);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("deg"));
    SET_VECTOR_ELT(s_df, col++, s_deg);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("p.deg"));
    SET_VECTOR_ELT(s_df, col++, s_p_deg);

    // Set names
    Rf_setAttrib(s_df, R_NamesSymbol, s_colnames);

    // Set row names
    SEXP s_rownames = PROTECT(Rf_allocVector(INTSXP, n_total));
    int* p_rownames = INTEGER(s_rownames);
    for (int i = 0; i < n_total; ++i) {
        p_rownames[i] = i + 1;
    }
    Rf_setAttrib(s_df, R_RowNamesSymbol, s_rownames);

    // Set class
    SEXP s_class = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(s_df, R_ClassSymbol, s_class);

    UNPROTECT(15);  // All allocations

    return s_df;
}

// ============================================================================
// Helper: Convert stage history to R data frame
// ============================================================================

static SEXP stage_history_to_dataframe(const std::vector<stage_counts_t>& history) {
    const int n = static_cast<int>(history.size());

    if (n == 0) {
        SEXP s_df = PROTECT(Rf_allocVector(VECSXP, 0));
        SEXP s_class = PROTECT(Rf_mkString("data.frame"));
        Rf_setAttrib(s_df, R_ClassSymbol, s_class);
        UNPROTECT(2);
        return s_df;
    }

    const int n_cols = 5;
    SEXP s_df = PROTECT(Rf_allocVector(VECSXP, n_cols));
    SEXP s_colnames = PROTECT(Rf_allocVector(STRSXP, n_cols));

    SEXP s_stage = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP s_n_max_before = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP s_n_max_after = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP s_n_min_before = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP s_n_min_after = PROTECT(Rf_allocVector(INTSXP, n));

    int* p_max_before = INTEGER(s_n_max_before);
    int* p_max_after = INTEGER(s_n_max_after);
    int* p_min_before = INTEGER(s_n_min_before);
    int* p_min_after = INTEGER(s_n_min_after);

    for (int i = 0; i < n; ++i) {
        SET_STRING_ELT(s_stage, i, Rf_mkChar(history[i].stage_name.c_str()));
        p_max_before[i] = history[i].n_max_before;
        p_max_after[i] = history[i].n_max_after;
        p_min_before[i] = history[i].n_min_before;
        p_min_after[i] = history[i].n_min_after;
    }

    int col = 0;
    SET_STRING_ELT(s_colnames, col, Rf_mkChar("stage"));
    SET_VECTOR_ELT(s_df, col++, s_stage);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("n.max.before"));
    SET_VECTOR_ELT(s_df, col++, s_n_max_before);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("n.max.after"));
    SET_VECTOR_ELT(s_df, col++, s_n_max_after);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("n.min.before"));
    SET_VECTOR_ELT(s_df, col++, s_n_min_before);

    SET_STRING_ELT(s_colnames, col, Rf_mkChar("n.min.after"));
    SET_VECTOR_ELT(s_df, col++, s_n_min_after);

    Rf_setAttrib(s_df, R_NamesSymbol, s_colnames);

    SEXP s_rownames = PROTECT(Rf_allocVector(INTSXP, n));
    for (int i = 0; i < n; ++i) {
        INTEGER(s_rownames)[i] = i + 1;
    }
    Rf_setAttrib(s_df, R_RowNamesSymbol, s_rownames);

    SEXP s_class = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(s_df, R_ClassSymbol, s_class);

    UNPROTECT(9);

    return s_df;
}

// ============================================================================
// Helper: Convert gfc_result_t to R list
// ============================================================================

static SEXP gfc_result_to_R(const gfc_result_t& result) {
    // Determine number of components (12 if trajectories, 10 otherwise)
    const bool has_trajectories = !result.joined_trajectories.empty();
    const int n_components = has_trajectories ? 12 : 10;

    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;
    const bool with_traj = result.params.with_trajectories;

    // max.basins - list of basins
    const int n_max = static_cast<int>(result.max_basins.size());
    SEXP s_max_basins = PROTECT(Rf_allocVector(VECSXP, n_max));
    for (int i = 0; i < n_max; ++i) {
        SET_VECTOR_ELT(s_max_basins, i, basin_compact_to_R(result.max_basins[i], with_traj));
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("max.basins"));
    SET_VECTOR_ELT(s_result, idx++, s_max_basins);
    UNPROTECT(1);

    // min.basins - list of basins
    const int n_min = static_cast<int>(result.min_basins.size());
    SEXP s_min_basins = PROTECT(Rf_allocVector(VECSXP, n_min));
    for (int i = 0; i < n_min; ++i) {
        SET_VECTOR_ELT(s_min_basins, i, basin_compact_to_R(result.min_basins[i], with_traj));
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("min.basins"));
    SET_VECTOR_ELT(s_result, idx++, s_min_basins);
    UNPROTECT(1);

    // summary - data frame
    SEXP s_summary = PROTECT(summaries_to_dataframe(
        result.max_summaries, result.min_summaries
    ));
    SET_STRING_ELT(s_names, idx, Rf_mkChar("summary"));
    SET_VECTOR_ELT(s_result, idx++, s_summary);
    UNPROTECT(1);

    // max.membership - list of integer vectors
    const int n = static_cast<int>(result.n_vertices);
    SEXP s_max_memb = PROTECT(Rf_allocVector(VECSXP, n));
    for (int v = 0; v < n; ++v) {
        const auto& memb = result.max_membership[v];
        const int n_memb = static_cast<int>(memb.size());
        SEXP s_memb_v = PROTECT(Rf_allocVector(INTSXP, n_memb));
        int* p = INTEGER(s_memb_v);
        for (int j = 0; j < n_memb; ++j) {
            p[j] = memb[j] + 1;  // 1-based basin index
        }
        SET_VECTOR_ELT(s_max_memb, v, s_memb_v);
        UNPROTECT(1);
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("max.membership"));
    SET_VECTOR_ELT(s_result, idx++, s_max_memb);
    UNPROTECT(1);

    // min.membership - list of integer vectors
    SEXP s_min_memb = PROTECT(Rf_allocVector(VECSXP, n));
    for (int v = 0; v < n; ++v) {
        const auto& memb = result.min_membership[v];
        const int n_memb = static_cast<int>(memb.size());
        SEXP s_memb_v = PROTECT(Rf_allocVector(INTSXP, n_memb));
        int* p = INTEGER(s_memb_v);
        for (int j = 0; j < n_memb; ++j) {
            p[j] = memb[j] + 1;
        }
        SET_VECTOR_ELT(s_min_memb, v, s_memb_v);
        UNPROTECT(1);
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("min.membership"));
    SET_VECTOR_ELT(s_result, idx++, s_min_memb);
    UNPROTECT(1);

    // expanded.max.assignment - integer vector (1-based, NA if unassigned)
    SEXP s_exp_max = PROTECT(Rf_allocVector(INTSXP, n));
    int* p_exp_max = INTEGER(s_exp_max);
    for (int v = 0; v < n; ++v) {
        int val = (v < static_cast<int>(result.expanded_max_assignment.size()))
            ? result.expanded_max_assignment[v] : -1;
        p_exp_max[v] = (val >= 0) ? val + 1 : NA_INTEGER;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("expanded.max.assignment"));
    SET_VECTOR_ELT(s_result, idx++, s_exp_max);
    UNPROTECT(1);

    // expanded.min.assignment - integer vector
    SEXP s_exp_min = PROTECT(Rf_allocVector(INTSXP, n));
    int* p_exp_min = INTEGER(s_exp_min);
    for (int v = 0; v < n; ++v) {
        int val = (v < static_cast<int>(result.expanded_min_assignment.size()))
            ? result.expanded_min_assignment[v] : -1;
        p_exp_min[v] = (val >= 0) ? val + 1 : NA_INTEGER;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("expanded.min.assignment"));
    SET_VECTOR_ELT(s_result, idx++, s_exp_min);
    UNPROTECT(1);

    // stage.history - data frame
    SEXP s_history = PROTECT(stage_history_to_dataframe(result.stage_history));
    SET_STRING_ELT(s_names, idx, Rf_mkChar("stage.history"));
    SET_VECTOR_ELT(s_result, idx++, s_history);
    UNPROTECT(1);

    // n.vertices
    SET_STRING_ELT(s_names, idx, Rf_mkChar("n.vertices"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarInteger(static_cast<int>(result.n_vertices)));

    // y.median
    SET_STRING_ELT(s_names, idx, Rf_mkChar("y.median"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarReal(result.y_median));

    // NEW: joined.trajectories (if available)
    if (has_trajectories) {
        SEXP s_joined = PROTECT(joined_trajectories_to_R(result.joined_trajectories));
        SET_STRING_ELT(s_names, idx, Rf_mkChar("joined.trajectories"));
        SET_VECTOR_ELT(s_result, idx++, s_joined);
        UNPROTECT(1);

        // NEW: cell.map
        SEXP s_cell_map = PROTECT(cell_map_to_R(result.cell_map));
        SET_STRING_ELT(s_names, idx, Rf_mkChar("cell.map"));
        SET_VECTOR_ELT(s_result, idx++, s_cell_map);
        UNPROTECT(1);
    }

    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    UNPROTECT(2);  // s_result, s_names

    return s_result;
}

// ============================================================================
// Main SEXP Interface Functions
// ============================================================================

extern "C" SEXP S_compute_gfc(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_min_rel_value_max,
    SEXP s_max_rel_value_min,
    SEXP s_max_overlap_threshold,
    SEXP s_min_overlap_threshold,
    SEXP s_p_mean_nbrs_dist_threshold,
    SEXP s_p_mean_hopk_dist_threshold,
    SEXP s_p_deg_threshold,
    SEXP s_min_basin_size,
    SEXP s_expand_basins,
    SEXP s_apply_relvalue_filter,
    SEXP s_apply_maxima_clustering,
    SEXP s_apply_minima_clustering,
    SEXP s_apply_geometric_filter,
    SEXP s_hop_k,
    SEXP s_with_trajectories,
    SEXP s_max_chain_depth,
    SEXP s_verbose
) {
    // Extract parameters
    gfc_params_t params;
    params.edge_length_quantile_thld = Rf_asReal(s_edge_length_quantile_thld);
    params.min_rel_value_max = Rf_asReal(s_min_rel_value_max);
    params.max_rel_value_min = Rf_asReal(s_max_rel_value_min);
    params.max_overlap_threshold = Rf_asReal(s_max_overlap_threshold);
    params.min_overlap_threshold = Rf_asReal(s_min_overlap_threshold);
    params.p_mean_nbrs_dist_threshold = Rf_asReal(s_p_mean_nbrs_dist_threshold);
    params.p_mean_hopk_dist_threshold = Rf_asReal(s_p_mean_hopk_dist_threshold);
    params.p_deg_threshold = Rf_asReal(s_p_deg_threshold);
    params.min_basin_size = Rf_asInteger(s_min_basin_size);
    params.expand_basins = Rf_asLogical(s_expand_basins);
    params.apply_relvalue_filter = Rf_asLogical(s_apply_relvalue_filter);
    params.apply_maxima_clustering = Rf_asLogical(s_apply_maxima_clustering);
    params.apply_minima_clustering = Rf_asLogical(s_apply_minima_clustering);
    params.apply_geometric_filter = Rf_asLogical(s_apply_geometric_filter);
    params.hop_k = Rf_asInteger(s_hop_k);
    params.with_trajectories = Rf_asLogical(s_with_trajectories);
    params.max_chain_depth = Rf_asInteger(s_max_chain_depth);

    bool verbose = Rf_asLogical(s_verbose);

    // Convert adjacency list and weights
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Extract y values
    const int n = LENGTH(s_y);
    std::vector<double> y(n);
    const double* p_y = REAL(s_y);
    for (int i = 0; i < n; ++i) {
        y[i] = p_y[i];
    }

    // Compute GFC
    gfc_result_t result = compute_gfc(graph, y, params, verbose);

    // Convert to R
    return gfc_result_to_R(result);
}

extern "C" SEXP S_compute_gfc_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_min_rel_value_max,
    SEXP s_max_rel_value_min,
    SEXP s_max_overlap_threshold,
    SEXP s_min_overlap_threshold,
    SEXP s_p_mean_nbrs_dist_threshold,
    SEXP s_p_mean_hopk_dist_threshold,
    SEXP s_p_deg_threshold,
    SEXP s_min_basin_size,
    SEXP s_expand_basins,
    SEXP s_apply_relvalue_filter,
    SEXP s_apply_maxima_clustering,
    SEXP s_apply_minima_clustering,
    SEXP s_apply_geometric_filter,
    SEXP s_hop_k,
    SEXP s_n_cores,
    SEXP s_verbose
) {
    // Extract parameters
    gfc_params_t params;
    params.edge_length_quantile_thld = Rf_asReal(s_edge_length_quantile_thld);
    params.min_rel_value_max = Rf_asReal(s_min_rel_value_max);
    params.max_rel_value_min = Rf_asReal(s_max_rel_value_min);
    params.max_overlap_threshold = Rf_asReal(s_max_overlap_threshold);
    params.min_overlap_threshold = Rf_asReal(s_min_overlap_threshold);
    params.p_mean_nbrs_dist_threshold = Rf_asReal(s_p_mean_nbrs_dist_threshold);
    params.p_mean_hopk_dist_threshold = Rf_asReal(s_p_mean_hopk_dist_threshold);
    params.p_deg_threshold = Rf_asReal(s_p_deg_threshold);
    params.min_basin_size = Rf_asInteger(s_min_basin_size);
    params.expand_basins = Rf_asLogical(s_expand_basins);
    params.apply_relvalue_filter = Rf_asLogical(s_apply_relvalue_filter);
    params.apply_maxima_clustering = Rf_asLogical(s_apply_maxima_clustering);
    params.apply_minima_clustering = Rf_asLogical(s_apply_minima_clustering);
    params.apply_geometric_filter = Rf_asLogical(s_apply_geometric_filter);
    params.hop_k = Rf_asInteger(s_hop_k);

    int n_cores = Rf_asInteger(s_n_cores);
    bool verbose = Rf_asLogical(s_verbose);

    // Convert adjacency list and weights
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Extract Y matrix dimensions
    SEXP s_dim = Rf_getAttrib(s_Y, R_DimSymbol);
    const int n = INTEGER(s_dim)[0];
    const int p = INTEGER(s_dim)[1];

    // Convert Y to Eigen matrix
    Eigen::MatrixXd Y(n, p);
    const double* p_Y = REAL(s_Y);
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i < n; ++i) {
            Y(i, j) = p_Y[i + j * n];  // Column-major
        }
    }

    // Compute GFC for all columns
    std::vector<gfc_result_t> results = compute_gfc_matrix(
        graph, Y, params, n_cores, verbose
    );

    // Convert to R list of results
    SEXP s_results = PROTECT(Rf_allocVector(VECSXP, p));
    for (int j = 0; j < p; ++j) {
        SET_VECTOR_ELT(s_results, j, gfc_result_to_R(results[j]));
    }

    // Add names if Y has column names
    SEXP s_dimnames = Rf_getAttrib(s_Y, R_DimNamesSymbol);
    if (!Rf_isNull(s_dimnames) && LENGTH(s_dimnames) >= 2) {
        SEXP s_colnames = VECTOR_ELT(s_dimnames, 1);
        if (!Rf_isNull(s_colnames)) {
            Rf_setAttrib(s_results, R_NamesSymbol, s_colnames);
        }
    }

    UNPROTECT(1);

    return s_results;
}

void debug_extension_search(
    size_t query_vertex,
    bool query_is_min,
    const std::unordered_map<size_t, gradient_basin_t>& all_max_basins,
    const std::unordered_map<size_t, gradient_basin_t>& all_min_basins,
    const std::unordered_set<size_t>& non_spurious_max,
    const std::unordered_set<size_t>& non_spurious_min,
    int max_depth
    );

extern "C" SEXP S_debug_extension_search(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_query_vertex,      // 1-based
    SEXP s_query_is_min,
    SEXP s_edge_length_quantile_thld,
    SEXP s_max_depth,
    SEXP s_non_spurious_max,  // 1-based integer vector
    SEXP s_non_spurious_min   // 1-based integer vector
) {
    // Convert inputs
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    const int n = LENGTH(s_y);
    std::vector<double> y(n);
    const double* p_y = REAL(s_y);
    for (int i = 0; i < n; ++i) y[i] = p_y[i];

    size_t query_vertex = static_cast<size_t>(Rf_asInteger(s_query_vertex) - 1);
    bool query_is_min = Rf_asLogical(s_query_is_min);
    double edge_thld = graph.compute_quantile_edge_length(Rf_asReal(s_edge_length_quantile_thld));
    int max_depth = Rf_asInteger(s_max_depth);

    // Build non-spurious sets
    std::unordered_set<size_t> non_spurious_max, non_spurious_min;
    int* p_max = INTEGER(s_non_spurious_max);
    for (int i = 0; i < LENGTH(s_non_spurious_max); ++i) {
        non_spurious_max.insert(static_cast<size_t>(p_max[i] - 1));
    }
    int* p_min = INTEGER(s_non_spurious_min);
    for (int i = 0; i < LENGTH(s_non_spurious_min); ++i) {
        non_spurious_min.insert(static_cast<size_t>(p_min[i] - 1));
    }

    // Compute ALL basins with trajectories
    auto [lmin_vertices, lmax_vertices] = graph.find_nbr_extrema(y);

    std::unordered_map<size_t, gradient_basin_t> all_max_basins, all_min_basins;

    for (size_t v : lmax_vertices) {
        auto gb = graph.compute_geodesic_basin(v, y, true, edge_thld, true);
        if (gb.hop_idx != std::numeric_limits<size_t>::max()) {
            all_max_basins[v] = std::move(gb);
        }
    }
    for (size_t v : lmin_vertices) {
        auto gb = graph.compute_geodesic_basin(v, y, false, edge_thld, true);
        if (gb.hop_idx != std::numeric_limits<size_t>::max()) {
            all_min_basins[v] = std::move(gb);
        }
    }

    Rprintf("Computed %zu max basins, %zu min basins\n",
            all_max_basins.size(), all_min_basins.size());

    // Run debug
    debug_extension_search(
        query_vertex, query_is_min,
        all_max_basins, all_min_basins,
        non_spurious_max, non_spurious_min,
        max_depth
    );

    return R_NilValue;
}
