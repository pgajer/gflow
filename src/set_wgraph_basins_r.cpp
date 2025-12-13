#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include <vector>

#include <R.h>
#include <Rinternals.h>

/**
 * @brief SEXP interface for computing gradient basins of attraction for local extrema.
 *
 * @details This function provides an R interface to the gradient basin computation.
 *          It identifies all local extrema in the graph and computes their gradient basins,
 *          returning separate lists for minima and maxima basins with complete trajectory
 *          information including predecessors and terminal extrema.
 *
 *          For each local extremum, the basin contains all vertices reachable via
 *          monotone paths (ascending for minima, descending for maxima), along with
 *          predecessor information for trajectory reconstruction and identification
 *          of terminal extrema where trajectories terminate.
 *
 * @param s_adj_list R list of integer vectors representing the adjacency list (0-based in C++)
 * @param s_weight_list R list of numeric vectors representing edge weights
 * @param s_y R numeric vector of function values at each vertex
 * @param s_edge_length_quantile_thld Edge length threshold for basin construction
 * @param s_with_trajectories Set to true to return gradient trajectories
 *
 * @return R list with two components:
 *         - lmin_basins: list of basin structures for local minima
 *         - lmax_basins: list of basin structures for local maxima
 *         Each basin structure contains:
 *         - vertex: 1-based vertex index of the extremum
 *         - value: function value at the extremum
 *         - hop_idx: maximum hop distance within the basin
 *         - basin_df: matrix with columns (vertex, hop_distance) for all basin members
 *         - basin_bd_df: matrix with columns (vertex, y_value) for boundary vertices
 *         - predecessors: integer vector where predecessors[i] is predecessor of vertex i (0 = no predecessor)
 *         - terminal_extrema: integer vector of terminal extrema vertices (1-based)
 */
extern "C" SEXP S_compute_basins_of_attraction(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_with_trajectories
    ) {

    // Input validation and conversion
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // ------ s_edge_length_quantile_thld validation/conversion
    if (!Rf_isReal(s_edge_length_quantile_thld) || LENGTH(s_edge_length_quantile_thld) != 1) {
            Rf_error("edge_length_quantile_thld must be a single numeric value");
    }

    double edge_length_quantile_thld = REAL(s_edge_length_quantile_thld)[0];

    if (std::isnan(edge_length_quantile_thld)) {
        Rf_error("edge_length_quantile_thld cannot be NA");
    }

    if (edge_length_quantile_thld < 0.0) {
        Rf_error("edge_length_quantile_thld must be non-negative");
    }

    // ------ s_with_trajectories conversion
    bool with_trajectories = Rf_asLogical(s_with_trajectories);

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    double edge_length_thld = graph.compute_quantile_edge_length(edge_length_quantile_thld);

    // Identify local extrema (same as before)
    size_t n = y.size();
    std::vector<size_t> local_minima;
    std::vector<size_t> local_maxima;

    for (size_t v = 0; v < n; ++v) {
        bool is_local_min = true;
        bool is_local_max = true;

        for (const auto& edge : graph.adjacency_list[v]) {
            size_t u = edge.vertex;
            if (y[u] >= y[v]) is_local_max = false;
            if (y[u] <= y[v]) is_local_min = false;
        }

        if (is_local_min) local_minima.push_back(v);
        if (is_local_max) local_maxima.push_back(v);
    }

    // Compute basins using GEODESIC method
    std::vector<gradient_basin_t> min_basins;
    for (size_t m : local_minima) {
        gradient_basin_t basin = graph.compute_geodesic_basin(m, y, false, edge_length_thld, with_trajectories);
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            min_basins.push_back(basin);
        }
    }

    std::vector<gradient_basin_t> max_basins;
    for (size_t M : local_maxima) {
        gradient_basin_t basin = graph.compute_geodesic_basin(M, y, true, edge_length_thld, with_trajectories);
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            max_basins.push_back(basin);
        }
    }

    // Convert to R format
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_result_names, 0, Rf_mkChar("lmin_basins"));
    SET_STRING_ELT(r_result_names, 1, Rf_mkChar("lmax_basins"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);

    // Helper lambda to convert basin to R
    auto basin_to_r = [&](const gradient_basin_t& basin) -> SEXP {
        bool is_empty_basin = (basin.hop_idx == std::numeric_limits<size_t>::max());

        // Determine number of fields
        int n_fields = with_trajectories ? 8 : 7;

        SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, n_fields));
        SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, n_fields));

        int field_idx = 0;

        // Standard fields (0-6)
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("vertex"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("value"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("hop_idx"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("basin_df"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("basin_bd_df"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("terminal_extrema"));

        if (with_trajectories) {
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("all_predecessors"));
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("trajectory_sets"));
        } else {
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("all_predecessors"));
        }

        Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);

        field_idx = 0;

        // vertex (1-based)
        SEXP r_vertex = PROTECT(Rf_ScalarInteger(basin.vertex + 1));
        SET_VECTOR_ELT(r_basin, field_idx++, r_vertex);

        // value
        SEXP r_value = PROTECT(Rf_ScalarReal(basin.value));
        SET_VECTOR_ELT(r_basin, field_idx++, r_value);

        // hop_idx
        SEXP r_hop_idx;
        if (is_empty_basin) {
            r_hop_idx = PROTECT(Rf_ScalarInteger(NA_INTEGER));
        } else {
            r_hop_idx = PROTECT(Rf_ScalarInteger(basin.hop_idx));
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_hop_idx);

        // basin_df matrix
        SEXP r_basin_df;
        if (is_empty_basin) {
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.hop_dist_map.size();
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
            double *pr = REAL(r_basin_df);
            size_t i = 0;
            for (const auto& [v, d] : basin.hop_dist_map) {
                pr[i]     = v + 1;
                pr[i + m] = d;
                i++;
            }
        }

        SEXP basin_colnames = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(basin_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_colnames, 1, Rf_mkChar("hop_distance"));
        SEXP basin_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(basin_dimnames, 0, R_NilValue);
        SET_VECTOR_ELT(basin_dimnames, 1, basin_colnames);
        Rf_setAttrib(r_basin_df, R_DimNamesSymbol, basin_dimnames);
        SET_VECTOR_ELT(r_basin, field_idx++, r_basin_df);
        UNPROTECT(2);

        // basin_bd_df matrix
        SEXP r_basin_bd_df;
        if (is_empty_basin) {
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.y_nbhd_bd_map.size();
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
            double *pr_bd = REAL(r_basin_bd_df);
            size_t i = 0;
            for (const auto& [v, y_val] : basin.y_nbhd_bd_map) {
                pr_bd[i]     = v + 1;
                pr_bd[i + m] = y_val;
                i++;
            }
        }

        SEXP basin_bd_colnames = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(basin_bd_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_bd_colnames, 1, Rf_mkChar("y_value"));
        SEXP basin_bd_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(basin_bd_dimnames, 0, R_NilValue);
        SET_VECTOR_ELT(basin_bd_dimnames, 1, basin_bd_colnames);
        Rf_setAttrib(r_basin_bd_df, R_DimNamesSymbol, basin_bd_dimnames);
        SET_VECTOR_ELT(r_basin, field_idx++, r_basin_bd_df);
        UNPROTECT(2);

        // terminal_extrema vector
        SEXP r_term_extrema = PROTECT(Rf_allocVector(INTSXP, basin.terminal_extrema.size()));
        int *p_term = INTEGER(r_term_extrema);
        for (size_t i = 0; i < basin.terminal_extrema.size(); ++i) {
            p_term[i] = basin.terminal_extrema[i] + 1;
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_term_extrema);

        // all_predecessors (empty if not with_trajectories)
        SEXP r_all_pred;
        if (with_trajectories) {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, n));
            for (size_t v = 0; v < n; ++v) {
                auto it = basin.all_predecessors.find(v);
                if (it != basin.all_predecessors.end() && !it->second.empty()) {
                    const std::vector<size_t>& preds = it->second;
                    SEXP r_preds_v = PROTECT(Rf_allocVector(INTSXP, preds.size()));
                    int *p_preds_v = INTEGER(r_preds_v);
                    for (size_t j = 0; j < preds.size(); ++j) {
                        p_preds_v[j] = preds[j] + 1;
                    }
                    SET_VECTOR_ELT(r_all_pred, v, r_preds_v);
                    UNPROTECT(1);
                } else {
                    SET_VECTOR_ELT(r_all_pred, v, Rf_allocVector(INTSXP, 0));
                }
            }
        } else {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, 0));
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_all_pred);

        // trajectory_sets (only if with_trajectories)
        if (with_trajectories) {
            SEXP r_traj_sets = PROTECT(Rf_allocVector(VECSXP, basin.trajectory_sets.size()));

            for (size_t i = 0; i < basin.trajectory_sets.size(); ++i) {
                const trajectory_set_t& traj_set = basin.trajectory_sets[i];

                SEXP r_tset = PROTECT(Rf_allocVector(VECSXP, 2));
                SEXP r_tset_names = PROTECT(Rf_allocVector(STRSXP, 2));
                SET_STRING_ELT(r_tset_names, 0, Rf_mkChar("terminal_vertex"));
                SET_STRING_ELT(r_tset_names, 1, Rf_mkChar("trajectories"));
                Rf_setAttrib(r_tset, R_NamesSymbol, r_tset_names);

                // terminal_vertex
                SEXP r_tvert = PROTECT(Rf_ScalarInteger(traj_set.terminal_vertex + 1));
                SET_VECTOR_ELT(r_tset, 0, r_tvert);

                // trajectories (list of integer vectors)
                SEXP r_trajs = PROTECT(Rf_allocVector(VECSXP, traj_set.trajectories.size()));
                for (size_t j = 0; j < traj_set.trajectories.size(); ++j) {
                    const std::vector<size_t>& traj = traj_set.trajectories[j];
                    SEXP r_traj = PROTECT(Rf_allocVector(INTSXP, traj.size()));
                    int *p_traj = INTEGER(r_traj);
                    for (size_t k = 0; k < traj.size(); ++k) {
                        p_traj[k] = traj[k] + 1;
                    }
                    SET_VECTOR_ELT(r_trajs, j, r_traj);
                    UNPROTECT(1);
                }
                SET_VECTOR_ELT(r_tset, 1, r_trajs);

                SET_VECTOR_ELT(r_traj_sets, i, r_tset);
                UNPROTECT(4); // r_tset, r_tset_names, r_tvert, r_trajs
            }

            SET_VECTOR_ELT(r_basin, field_idx++, r_traj_sets);
            UNPROTECT(1); // r_traj_sets
        }

        UNPROTECT(9); // r_basin, r_basin_names, r_vertex, r_value, r_hop_idx, r_basin_df, r_basin_bd_df, r_term_extrema, r_all_pred
        return r_basin;
    };

    // Convert basins
    SEXP r_min_basins = PROTECT(Rf_allocVector(VECSXP, min_basins.size()));
    for (size_t i = 0; i < min_basins.size(); ++i) {
        SEXP r_basin = basin_to_r(min_basins[i]);
        SET_VECTOR_ELT(r_min_basins, i, r_basin);
    }
    SET_VECTOR_ELT(r_result, 0, r_min_basins);

    SEXP r_max_basins = PROTECT(Rf_allocVector(VECSXP, max_basins.size()));
    for (size_t i = 0; i < max_basins.size(); ++i) {
        SEXP r_basin = basin_to_r(max_basins[i]);
        SET_VECTOR_ELT(r_max_basins, i, r_basin);
    }
    SET_VECTOR_ELT(r_result, 1, r_max_basins);

    UNPROTECT(4);
    return r_result;
}


#if 0
extern "C" SEXP S_compute_basins_of_attraction(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_with_trajectories,
    SEXP s_k_paths
    ) {

    // Input validation and conversion
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    bool with_trajectories = Rf_asLogical(s_with_trajectories);
    size_t k_paths = static_cast<size_t>(Rf_asInteger(s_k_paths));

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Identify local extrema
    size_t n = y.size();
    std::vector<size_t> local_minima;
    std::vector<size_t> local_maxima;

    for (size_t v = 0; v < n; ++v) {
        bool is_local_min = true;
        bool is_local_max = true;

        for (const auto& edge : graph.adjacency_list[v]) {
            size_t u = edge.vertex;
            if (y[u] >= y[v]) is_local_max = false;
            if (y[u] <= y[v]) is_local_min = false;
        }

        if (is_local_min) local_minima.push_back(v);
        if (is_local_max) local_maxima.push_back(v);
    }

    // Compute basins
    std::vector<gradient_basin_t> min_basins;
    for (size_t m : local_minima) {
        gradient_basin_t basin = graph.compute_basin_of_attraction(m, y, false, with_trajectories, k_paths);
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            min_basins.push_back(basin);
        }
    }

    std::vector<gradient_basin_t> max_basins;
    for (size_t M : local_maxima) {
        gradient_basin_t basin = graph.compute_basin_of_attraction(M, y, true, with_trajectories, k_paths);
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            max_basins.push_back(basin);
        }
    }

    // Convert to R format
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_result_names, 0, Rf_mkChar("lmin_basins"));
    SET_STRING_ELT(r_result_names, 1, Rf_mkChar("lmax_basins"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);

    // Helper lambda to convert basin to R
    auto basin_to_r = [&](const gradient_basin_t& basin) -> SEXP {
        bool is_empty_basin = (basin.hop_idx == std::numeric_limits<size_t>::max());

        // Determine number of fields
        int n_fields = with_trajectories ? 8 : 7;

        SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, n_fields));
        SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, n_fields));

        int field_idx = 0;

        // Standard fields (0-6)
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("vertex"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("value"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("hop_idx"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("basin_df"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("basin_bd_df"));
        SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("terminal_extrema"));

        if (with_trajectories) {
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("all_predecessors"));
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("trajectory_sets"));
        } else {
            SET_STRING_ELT(r_basin_names, field_idx++, Rf_mkChar("all_predecessors"));
        }

        Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);

        field_idx = 0;

        // vertex (1-based)
        SEXP r_vertex = PROTECT(Rf_ScalarInteger(basin.vertex + 1));
        SET_VECTOR_ELT(r_basin, field_idx++, r_vertex);

        // value
        SEXP r_value = PROTECT(Rf_ScalarReal(basin.value));
        SET_VECTOR_ELT(r_basin, field_idx++, r_value);

        // hop_idx
        SEXP r_hop_idx;
        if (is_empty_basin) {
            r_hop_idx = PROTECT(Rf_ScalarInteger(NA_INTEGER));
        } else {
            r_hop_idx = PROTECT(Rf_ScalarInteger(basin.hop_idx));
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_hop_idx);

        // basin_df matrix
        SEXP r_basin_df;
        if (is_empty_basin) {
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.hop_dist_map.size();
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
            double *pr = REAL(r_basin_df);
            size_t i = 0;
            for (const auto& [v, d] : basin.hop_dist_map) {
                pr[i]     = v + 1;
                pr[i + m] = d;
                i++;
            }
        }

        SEXP basin_colnames = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(basin_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_colnames, 1, Rf_mkChar("hop_distance"));
        SEXP basin_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(basin_dimnames, 0, R_NilValue);
        SET_VECTOR_ELT(basin_dimnames, 1, basin_colnames);
        Rf_setAttrib(r_basin_df, R_DimNamesSymbol, basin_dimnames);
        SET_VECTOR_ELT(r_basin, field_idx++, r_basin_df);
        UNPROTECT(2);

        // basin_bd_df matrix
        SEXP r_basin_bd_df;
        if (is_empty_basin) {
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.y_nbhd_bd_map.size();
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
            double *pr_bd = REAL(r_basin_bd_df);
            size_t i = 0;
            for (const auto& [v, y_val] : basin.y_nbhd_bd_map) {
                pr_bd[i]     = v + 1;
                pr_bd[i + m] = y_val;
                i++;
            }
        }

        SEXP basin_bd_colnames = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(basin_bd_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_bd_colnames, 1, Rf_mkChar("y_value"));
        SEXP basin_bd_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(basin_bd_dimnames, 0, R_NilValue);
        SET_VECTOR_ELT(basin_bd_dimnames, 1, basin_bd_colnames);
        Rf_setAttrib(r_basin_bd_df, R_DimNamesSymbol, basin_bd_dimnames);
        SET_VECTOR_ELT(r_basin, field_idx++, r_basin_bd_df);
        UNPROTECT(2);

        // terminal_extrema vector
        SEXP r_term_extrema = PROTECT(Rf_allocVector(INTSXP, basin.terminal_extrema.size()));
        int *p_term = INTEGER(r_term_extrema);
        for (size_t i = 0; i < basin.terminal_extrema.size(); ++i) {
            p_term[i] = basin.terminal_extrema[i] + 1;
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_term_extrema);

        // all_predecessors (empty if not with_trajectories)
        SEXP r_all_pred;
        if (with_trajectories) {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, n));
            for (size_t v = 0; v < n; ++v) {
                auto it = basin.all_predecessors.find(v);
                if (it != basin.all_predecessors.end() && !it->second.empty()) {
                    const std::vector<size_t>& preds = it->second;
                    SEXP r_preds_v = PROTECT(Rf_allocVector(INTSXP, preds.size()));
                    int *p_preds_v = INTEGER(r_preds_v);
                    for (size_t j = 0; j < preds.size(); ++j) {
                        p_preds_v[j] = preds[j] + 1;
                    }
                    SET_VECTOR_ELT(r_all_pred, v, r_preds_v);
                    UNPROTECT(1);
                } else {
                    SET_VECTOR_ELT(r_all_pred, v, Rf_allocVector(INTSXP, 0));
                }
            }
        } else {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, 0));
        }
        SET_VECTOR_ELT(r_basin, field_idx++, r_all_pred);

        // trajectory_sets (only if with_trajectories)
        if (with_trajectories) {
            SEXP r_traj_sets = PROTECT(Rf_allocVector(VECSXP, basin.trajectory_sets.size()));

            for (size_t i = 0; i < basin.trajectory_sets.size(); ++i) {
                const trajectory_set_t& traj_set = basin.trajectory_sets[i];

                SEXP r_tset = PROTECT(Rf_allocVector(VECSXP, 2));
                SEXP r_tset_names = PROTECT(Rf_allocVector(STRSXP, 2));
                SET_STRING_ELT(r_tset_names, 0, Rf_mkChar("terminal_vertex"));
                SET_STRING_ELT(r_tset_names, 1, Rf_mkChar("trajectories"));
                Rf_setAttrib(r_tset, R_NamesSymbol, r_tset_names);

                // terminal_vertex
                SEXP r_tvert = PROTECT(Rf_ScalarInteger(traj_set.terminal_vertex + 1));
                SET_VECTOR_ELT(r_tset, 0, r_tvert);

                // trajectories (list of integer vectors)
                SEXP r_trajs = PROTECT(Rf_allocVector(VECSXP, traj_set.trajectories.size()));
                for (size_t j = 0; j < traj_set.trajectories.size(); ++j) {
                    const std::vector<size_t>& traj = traj_set.trajectories[j];
                    SEXP r_traj = PROTECT(Rf_allocVector(INTSXP, traj.size()));
                    int *p_traj = INTEGER(r_traj);
                    for (size_t k = 0; k < traj.size(); ++k) {
                        p_traj[k] = traj[k] + 1;
                    }
                    SET_VECTOR_ELT(r_trajs, j, r_traj);
                    UNPROTECT(1);
                }
                SET_VECTOR_ELT(r_tset, 1, r_trajs);

                SET_VECTOR_ELT(r_traj_sets, i, r_tset);
                UNPROTECT(4); // r_tset, r_tset_names, r_tvert, r_trajs
            }

            SET_VECTOR_ELT(r_basin, field_idx++, r_traj_sets);
            UNPROTECT(1); // r_traj_sets
        }

        UNPROTECT(9); // r_basin, r_basin_names, r_vertex, r_value, r_hop_idx, r_basin_df, r_basin_bd_df, r_term_extrema, r_all_pred
        return r_basin;
    };

    // Convert basins
    SEXP r_min_basins = PROTECT(Rf_allocVector(VECSXP, min_basins.size()));
    for (size_t i = 0; i < min_basins.size(); ++i) {
        SEXP r_basin = basin_to_r(min_basins[i]);
        SET_VECTOR_ELT(r_min_basins, i, r_basin);
    }
    SET_VECTOR_ELT(r_result, 0, r_min_basins);

    SEXP r_max_basins = PROTECT(Rf_allocVector(VECSXP, max_basins.size()));
    for (size_t i = 0; i < max_basins.size(); ++i) {
        SEXP r_basin = basin_to_r(max_basins[i]);
        SET_VECTOR_ELT(r_max_basins, i, r_basin);
    }
    SET_VECTOR_ELT(r_result, 1, r_max_basins);

    UNPROTECT(4);
    return r_result;
}
#endif
