// graph_gradient_flow_hardened.cpp
// Hardened drop-in for S_construct_graph_gradient_flow
// - Defensive type/length checks
// - Container-first PROTECT pattern everywhere
// - Size-safe casts (R_xlen_t -> int) with overflow guards
// - No UNPROTECT after errors; no variable UNPROTECT counts
// - Avoids PROTECT gaps when allocating nested structures

#include <vector>
#include <unordered_map>
#include <set>
#include <string>
#include <limits>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available in the package)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Headers from the package
#include "set_wgraph.hpp"
#include "gradient_flow.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Structures (mirroring originals)
struct gradient_flow_t {
    enum TrajectoryType { LMIN_LMAX, LMIN_ONLY, LMAX_ONLY };
    struct trajectory_t {
        std::vector<size_t> vertices;
        TrajectoryType trajectory_type;
        double quality_metric;
        double total_change;
    };
    std::unordered_map<size_t, bool> local_extrema;
    std::vector<trajectory_t> trajectories;
    std::unordered_map<size_t, std::set<size_t>> ascending_basin_map;
    std::unordered_map<size_t, std::set<size_t>> descending_basin_map;
    std::unordered_map<std::pair<size_t, size_t>, std::set<size_t>,
                       std::hash<std::pair<size_t, size_t>>> cell_map;
    std::string messages;
};

class set_wgraph_t {
public:
    set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                 const std::vector<std::vector<double>>& weight_list);
    gradient_flow_t compute_gradient_flow(const std::vector<double>& y,
                                          const std::vector<double>& scale,
                                          double quantile_scale_thld);
};

// ---- helpers ----
static inline int safe_len_cast(R_xlen_t x, const char* what) {
    if (x < 0 || x > std::numeric_limits<int>::max()) {
        Rf_error("S_construct_graph_gradient_flow(): %s length (%lld) exceeds INT_MAX or is negative",
                 what, (long long)x);
    }
    return static_cast<int>(x);
}

extern "C" SEXP S_construct_graph_gradient_flow(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_scale,
    SEXP s_quantile_scale_thld,
    SEXP s_with_trajectories
) {
    // ---- defensive checks (no allocations yet) ----
    if (TYPEOF(s_y) != REALSXP)  Rf_error("S_construct_graph_gradient_flow(): 'y' must be a double vector");
    if (TYPEOF(s_scale) != REALSXP) Rf_error("S_construct_graph_gradient_flow(): 'scale' must be a double vector");
    if (TYPEOF(s_quantile_scale_thld) != REALSXP || Rf_length(s_quantile_scale_thld) != 1)
        Rf_error("S_construct_graph_gradient_flow(): 'quantile_scale_thld' must be a length-1 double");
    int with_traj = Rf_asLogical(s_with_trajectories);
    if (with_traj == NA_LOGICAL)
        Rf_error("S_construct_graph_gradient_flow(): 'with_trajectories' must be TRUE/FALSE");

    R_xlen_t n_vertices_x = Rf_xlength(s_y);
    if (n_vertices_x <= 0) Rf_error("S_construct_graph_gradient_flow(): 'y' must be non-empty");
    if (Rf_xlength(s_scale) != n_vertices_x)
        Rf_error("S_construct_graph_gradient_flow(): 'scale' length (%lld) must equal 'y' length (%lld)",
                 (long long)Rf_xlength(s_scale), (long long)n_vertices_x);

    // Convert lists (these helpers do their own validation)
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // ---- copy inputs into std::vector BEFORE any allocations ----
    const R_xlen_t n_vertices_x2 = n_vertices_x;
    std::vector<double> y;
    y.reserve(static_cast<size_t>(n_vertices_x2));
    {
        const double* y_ptr = REAL(s_y);
        y.assign(y_ptr, y_ptr + n_vertices_x2);
    }
    std::vector<double> scale;
    scale.reserve(static_cast<size_t>(n_vertices_x2));
    {
        const double* sc_ptr = REAL(s_scale);
        scale.assign(sc_ptr, sc_ptr + n_vertices_x2);
    }
    const double quantile_scale_thld = REAL(s_quantile_scale_thld)[0];
    const bool with_trajectories = (with_traj == TRUE);

    // ---- core compute (no R allocs) ----
    set_wgraph_t graph(adj_list, weight_list);
    gradient_flow_t result = graph.compute_gradient_flow(y, scale, quantile_scale_thld);

    // ---- build return object (container-first) ----
    const int n_elements = 5;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));

    // names
    {
        SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
        SET_STRING_ELT(r_names, 0, Rf_mkChar("local_extrema"));
        SET_STRING_ELT(r_names, 1, Rf_mkChar("trajectories"));
        SET_STRING_ELT(r_names, 2, Rf_mkChar("basins"));
        SET_STRING_ELT(r_names, 3, Rf_mkChar("cells"));
        SET_STRING_ELT(r_names, 4, Rf_mkChar("messages"));
        Rf_setAttrib(r_result, R_NamesSymbol, r_names);
        UNPROTECT(1); // r_names
    }

    // 0: local_extrema (matrix n x 2)
    {
        const R_xlen_t n_extrema_x = static_cast<R_xlen_t>(result.local_extrema.size());
        const int n_extrema = safe_len_cast(n_extrema_x, "local_extrema");
        SEXP r_extrema = PROTECT(Rf_allocMatrix(INTSXP, n_extrema, 2));
        int* extrema_ptr = INTEGER(r_extrema);
        int i = 0;
        for (const auto& kv : result.local_extrema) {
            const size_t vertex = kv.first;
            const bool is_maximum = kv.second;
            extrema_ptr[i] = static_cast<int>(vertex) + 1;
            extrema_ptr[i + n_extrema] = is_maximum ? 1 : 0;
            ++i;
        }
        SET_VECTOR_ELT(r_result, 0, r_extrema);
        UNPROTECT(1); // r_extrema
    }

    // 1: trajectories (list or NULL)
    if (with_trajectories) {
        const R_xlen_t n_traj_x = static_cast<R_xlen_t>(result.trajectories.size());
        const int n_traj = safe_len_cast(n_traj_x, "trajectories");
        SEXP r_trajectories = PROTECT(Rf_allocVector(VECSXP, n_traj));
        for (int i = 0; i < n_traj; ++i) {
            const auto& traj = result.trajectories[static_cast<size_t>(i)];
            SEXP r_traj = PROTECT(Rf_allocVector(VECSXP, 4)); // vertices, type, quality_metric, total_change)

            // vertices
            {
                const R_xlen_t n_v_x = static_cast<R_xlen_t>(traj.vertices.size());
                const int n_v = safe_len_cast(n_v_x, "trajectory vertices");
                SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_v));
                int* vptr = INTEGER(r_vertices);
                for (int j = 0; j < n_v; ++j) {
                    vptr[j] = static_cast<int>(traj.vertices[static_cast<size_t>(j)]) + 1;
                }
                SET_VECTOR_ELT(r_traj, 0, r_vertices);
                UNPROTECT(1); // r_vertices
            }

            // type
            {
                const char* type_str = "UNKNOWN";
                switch (traj.trajectory_type) {
                    case gradient_flow_t::LMIN_LMAX: type_str = "LMIN_LMAX"; break;
                    case gradient_flow_t::LMIN_ONLY: type_str = "LMIN_ONLY"; break;
                    case gradient_flow_t::LMAX_ONLY: type_str = "LMAX_ONLY"; break;
                    default: break;
                }
                SEXP r_type = PROTECT(Rf_allocVector(STRSXP, 1));
                SET_STRING_ELT(r_type, 0, Rf_mkChar(type_str));
                SET_VECTOR_ELT(r_traj, 1, r_type);
                UNPROTECT(1); // r_type
            }

            // quality_metric
            {
                SEXP r_quality = PROTECT(Rf_allocVector(REALSXP, 1));
                REAL(r_quality)[0] = traj.quality_metric;
                SET_VECTOR_ELT(r_traj, 2, r_quality);
                UNPROTECT(1); // r_quality
            }

            // total_change
            {
                SEXP r_change = PROTECT(Rf_allocVector(REALSXP, 1));
                REAL(r_change)[0] = traj.total_change;
                SET_VECTOR_ELT(r_traj, 3, r_change);
                UNPROTECT(1); // r_change
            }

            // names of trajectory
            {
                SEXP r_traj_names = PROTECT(Rf_allocVector(STRSXP, 4));
                SET_STRING_ELT(r_traj_names, 0, Rf_mkChar("vertices"));
                SET_STRING_ELT(r_traj_names, 1, Rf_mkChar("type"));
                SET_STRING_ELT(r_traj_names, 2, Rf_mkChar("quality_metric"));
                SET_STRING_ELT(r_traj_names, 3, Rf_mkChar("total_change"));
                Rf_setAttrib(r_traj, R_NamesSymbol, r_traj_names);
                UNPROTECT(1); // r_traj_names
            }

            SET_VECTOR_ELT(r_trajectories, i, r_traj);
            UNPROTECT(1); // r_traj
        }
        SET_VECTOR_ELT(r_result, 1, r_trajectories);
        UNPROTECT(1); // r_trajectories
    } else {
        SET_VECTOR_ELT(r_result, 1, R_NilValue);
    }

    // 2: basins (list: ascending, descending)
    {
        SEXP r_basins = PROTECT(Rf_allocVector(VECSXP, 2));

        // ascending
        {
            const R_xlen_t n_asc_x = static_cast<R_xlen_t>(result.ascending_basin_map.size());
            const int n_asc = safe_len_cast(n_asc_x, "ascending basins");
            SEXP r_ascending = PROTECT(Rf_allocVector(VECSXP, n_asc));
            int i = 0;
            for (const auto& kv : result.ascending_basin_map) {
                const size_t min_vertex = kv.first;
                const std::set<size_t>& basin = kv.second;

                SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 2)); // local_min, vertices

                SEXP r_min = PROTECT(Rf_allocVector(INTSXP, 1));
                INTEGER(r_min)[0] = static_cast<int>(min_vertex) + 1;

                const R_xlen_t n_b_x = static_cast<R_xlen_t>(basin.size());
                const int n_b = safe_len_cast(n_b_x, "ascending basin vertices");
                SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_b));
                {
                    int j = 0; int* vptr = INTEGER(r_vertices);
                    for (size_t v : basin) vptr[j++] = static_cast<int>(v) + 1;
                }

                SET_VECTOR_ELT(r_basin, 0, r_min);
                SET_VECTOR_ELT(r_basin, 1, r_vertices);
                UNPROTECT(2); // r_min, r_vertices

                SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, 2));
                SET_STRING_ELT(r_basin_names, 0, Rf_mkChar("local_min"));
                SET_STRING_ELT(r_basin_names, 1, Rf_mkChar("vertices"));
                Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);
                UNPROTECT(1); // r_basin_names

                SET_VECTOR_ELT(r_ascending, i++, r_basin);
                UNPROTECT(1); // r_basin
            }
            SET_VECTOR_ELT(r_basins, 0, r_ascending);
            UNPROTECT(1); // r_ascending
        }

        // descending
        {
            const R_xlen_t n_desc_x = static_cast<R_xlen_t>(result.descending_basin_map.size());
            const int n_desc = safe_len_cast(n_desc_x, "descending basins");
            SEXP r_descending = PROTECT(Rf_allocVector(VECSXP, n_desc));
            int i = 0;
            for (const auto& kv : result.descending_basin_map) {
                const size_t max_vertex = kv.first;
                const std::set<size_t>& basin = kv.second;

                SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, 2)); // local_max, vertices

                SEXP r_max = PROTECT(Rf_allocVector(INTSXP, 1));
                INTEGER(r_max)[0] = static_cast<int>(max_vertex) + 1;

                const R_xlen_t n_b_x = static_cast<R_xlen_t>(basin.size());
                const int n_b = safe_len_cast(n_b_x, "descending basin vertices");
                SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_b));
                {
                    int j = 0; int* vptr = INTEGER(r_vertices);
                    for (size_t v : basin) vptr[j++] = static_cast<int>(v) + 1;
                }

                SET_VECTOR_ELT(r_basin, 0, r_max);
                SET_VECTOR_ELT(r_basin, 1, r_vertices);
                UNPROTECT(2); // r_max, r_vertices

                SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, 2));
                SET_STRING_ELT(r_basin_names, 0, Rf_mkChar("local_max"));
                SET_STRING_ELT(r_basin_names, 1, Rf_mkChar("vertices"));
                Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);
                UNPROTECT(1); // r_basin_names

                SET_VECTOR_ELT(r_descending, i++, r_basin);
                UNPROTECT(1); // r_basin
            }
            SET_VECTOR_ELT(r_basins, 1, r_descending);
            UNPROTECT(1); // r_descending
        }

        SEXP r_basins_names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(r_basins_names, 0, Rf_mkChar("ascending"));
        SET_STRING_ELT(r_basins_names, 1, Rf_mkChar("descending"));
        Rf_setAttrib(r_basins, R_NamesSymbol, r_basins_names);
        UNPROTECT(1); // r_basins_names

        SET_VECTOR_ELT(r_result, 2, r_basins);
        UNPROTECT(1); // r_basins
    }

    // 3: cells
    {
        const R_xlen_t n_cells_x = static_cast<R_xlen_t>(result.cell_map.size());
        const int n_cells = safe_len_cast(n_cells_x, "cells");
        SEXP r_cells = PROTECT(Rf_allocVector(VECSXP, n_cells));
        int i = 0;
        for (const auto& kv : result.cell_map) {
            const std::pair<size_t,size_t>& key_pair = kv.first;
            const std::set<size_t>& vertices = kv.second;

            SEXP r_cell = PROTECT(Rf_allocVector(VECSXP, 3)); // local_min, local_max, vertices

            SEXP r_min = PROTECT(Rf_allocVector(INTSXP, 1));
            INTEGER(r_min)[0] = static_cast<int>(key_pair.first) + 1;

            SEXP r_max = PROTECT(Rf_allocVector(INTSXP, 1));
            INTEGER(r_max)[0] = static_cast<int>(key_pair.second) + 1;

            const R_xlen_t n_v_x = static_cast<R_xlen_t>(vertices.size());
            const int n_v = safe_len_cast(n_v_x, "cell vertices");
            SEXP r_vertices = PROTECT(Rf_allocVector(INTSXP, n_v));
            {
                int j = 0; int* vptr = INTEGER(r_vertices);
                for (size_t v : vertices) vptr[j++] = static_cast<int>(v) + 1;
            }

            SET_VECTOR_ELT(r_cell, 0, r_min);
            SET_VECTOR_ELT(r_cell, 1, r_max);
            SET_VECTOR_ELT(r_cell, 2, r_vertices);
            UNPROTECT(3); // r_min, r_max, r_vertices

            SEXP r_cell_names = PROTECT(Rf_allocVector(STRSXP, 3));
            SET_STRING_ELT(r_cell_names, 0, Rf_mkChar("local_min"));
            SET_STRING_ELT(r_cell_names, 1, Rf_mkChar("local_max"));
            SET_STRING_ELT(r_cell_names, 2, Rf_mkChar("vertices"));
            Rf_setAttrib(r_cell, R_NamesSymbol, r_cell_names);
            UNPROTECT(1); // r_cell_names

            SET_VECTOR_ELT(r_cells, i++, r_cell);
            UNPROTECT(1); // r_cell
        }
        SET_VECTOR_ELT(r_result, 3, r_cells);
        UNPROTECT(1); // r_cells
    }

    // 4: messages
    {
        SEXP r_messages = PROTECT(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(r_messages, 0, Rf_mkChar(result.messages.c_str()));
        SET_VECTOR_ELT(r_result, 4, r_messages);
        UNPROTECT(1); // r_messages
    }

    UNPROTECT(1); // r_result
    return r_result;
}