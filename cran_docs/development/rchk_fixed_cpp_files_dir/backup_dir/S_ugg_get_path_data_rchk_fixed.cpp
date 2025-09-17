
#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <string>
#include <utility>
#include <optional>
#include <numeric>

// ---- Project-provided types & functions (assumed available) ----
struct path_data_t {
    std::vector<size_t> vertices;
    size_t ref_vertex;
    double rel_center_offset;
    double total_weight;
    std::vector<double> x_path;
    std::vector<double> w_path;
    std::vector<double> y_path;
};

struct uniform_grid_graph_t {
    std::vector<std::vector<std::pair<size_t,double>>> adjacency_list;
    std::vector<size_t> grid_vertices;
};

using edge_weights_t = std::vector<std::vector<double>>;

extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

extern uniform_grid_graph_t create_uniform_grid_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    size_t grid_size,
    size_t ref_vertex,
    double snap_tolerance);

extern edge_weights_t precompute_edge_weights(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list);

extern std::vector<path_data_t> ugg_get_path_data(
    const uniform_grid_graph_t& ugg,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    const edge_weights_t& edge_weights,
    bool verbose);

// =========================== rchk-safe drop-in ===========================
extern "C" SEXP S_ugg_get_path_data(
    SEXP adj_list_s,
    SEXP weight_list_s,
    SEXP grid_size_s,
    SEXP y_s,
    SEXP ref_vertex_s,
    SEXP bandwidth_s,
    SEXP dist_normalization_factor_s,
    SEXP min_path_size_s,
    SEXP diff_threshold_s,
    SEXP kernel_type_s,
    SEXP verbose_s
) {
    // ---- Convert adjacency / weights with protection to guard allocating helpers ----
    std::vector<std::vector<int>>    adj_list;
    std::vector<std::vector<double>> weight_list;
    {
        SEXP a = adj_list_s, w = weight_list_s;
        PROTECT(a);
        PROTECT(w);
        adj_list    = convert_adj_list_from_R(a);
        weight_list = convert_weight_list_from_R(w);
        UNPROTECT(2);
    }

    // ---- y: coerce + copy under indexed PROTECT (long-vector safe) ----
    std::vector<double> y;
    {
        SEXP sy = y_s;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        const R_xlen_t ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + (size_t)ny);
        UNPROTECT(1);
    }

    // ---- Scalars / flags (defensive coercion) ----
    size_t grid_size                = (size_t) Rf_asInteger(grid_size_s);
    size_t ref_vertex               = (size_t) Rf_asInteger(ref_vertex_s);
    double bandwidth                = Rf_asReal(bandwidth_s);
    double dist_normalization_factor= Rf_asReal(dist_normalization_factor_s);
    size_t min_path_size            = (size_t) Rf_asInteger(min_path_size_s);
    size_t diff_threshold           = (size_t) Rf_asInteger(diff_threshold_s);
    size_t kernel_type              = (size_t) Rf_asInteger(kernel_type_s);
    bool   verbose                  = (Rf_asLogical(verbose_s) == TRUE);

    // ---- Core compute (no R allocations) ----
    const double snap_tolerance = 0.1;
    uniform_grid_graph_t uniform_grid_graph = create_uniform_grid_graph(
        adj_list, weight_list, grid_size, ref_vertex, snap_tolerance);

    edge_weights_t edge_weights = precompute_edge_weights(adj_list, weight_list);

    std::vector<path_data_t> paths = ugg_get_path_data(
        uniform_grid_graph, y, ref_vertex, bandwidth, dist_normalization_factor,
        min_path_size, diff_threshold, kernel_type, edge_weights, verbose);

    const R_xlen_t n_paths = (R_xlen_t)paths.size();

    // ---- Result container + names (kept protected to the end) ----
    const R_xlen_t RESULT_LIST_SIZE = n_paths + 3;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, RESULT_LIST_SIZE));
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, RESULT_LIST_SIZE));

    for (R_xlen_t i = 0; i < n_paths; ++i) {
        std::string path_id = std::string("path") + std::to_string((size_t)i);
        SET_STRING_ELT(result_names, i, Rf_mkChar(path_id.c_str()));
    }
    SET_STRING_ELT(result_names, n_paths + 0, Rf_mkChar("ugg_adj_list"));
    SET_STRING_ELT(result_names, n_paths + 1, Rf_mkChar("ugg_weight_list"));
    SET_STRING_ELT(result_names, n_paths + 2, Rf_mkChar("ugg_grid_vertices"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // ---- Path component names (protected across the loop) ----
    const char* path_comps_names[] = {
        "vertices", "ref_vertex", "rel_center_offset", "total_weight",
        "x_path", "w_path", "y_path"
    };
    const int N_PATH_COMPS = 7;
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_PATH_COMPS));
    for (int i = 0; i < N_PATH_COMPS; ++i) {
        SET_STRING_ELT(names, i, Rf_mkChar(path_comps_names[i]));
    }

    // ---- Populate path entries ----
    for (R_xlen_t i = 0; i < n_paths; ++i) {
        const path_data_t& pd = paths[(size_t)i];

        SEXP path = PROTECT(Rf_allocVector(VECSXP, N_PATH_COMPS));

        // 0: vertices (1-based)
        {
            const R_xlen_t nv = (R_xlen_t)pd.vertices.size();
            SEXP s = PROTECT(Rf_allocVector(INTSXP, nv));
            int* ip = INTEGER(s);
            for (R_xlen_t j = 0; j < nv; ++j) ip[j] = (int)pd.vertices[(size_t)j] + 1;
            SET_VECTOR_ELT(path, 0, s);
            UNPROTECT(1);
        }

        // 1: ref_vertex (1-based)
        {
            SEXP s = PROTECT(Rf_ScalarInteger((int)pd.ref_vertex + 1));
            SET_VECTOR_ELT(path, 1, s);
            UNPROTECT(1);
        }

        // 2: rel_center_offset
        {
            SEXP s = PROTECT(Rf_ScalarReal(pd.rel_center_offset));
            SET_VECTOR_ELT(path, 2, s);
            UNPROTECT(1);
        }

        // 3: total_weight
        {
            SEXP s = PROTECT(Rf_ScalarReal(pd.total_weight));
            SET_VECTOR_ELT(path, 3, s);
            UNPROTECT(1);
        }

        // helper to set a numeric vector
        auto set_num_vec = [&](int idx, const std::vector<double>& v){
            SEXP s = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)v.size()));
            double* rp = REAL(s);
            for (size_t k = 0; k < v.size(); ++k) rp[k] = v[k];
            SET_VECTOR_ELT(path, idx, s);
            UNPROTECT(1);
        };

        // 4..6: x_path, w_path, y_path
        set_num_vec(4, pd.x_path);
        set_num_vec(5, pd.w_path);
        set_num_vec(6, pd.y_path);

        // names for the path sublist
        Rf_setAttrib(path, R_NamesSymbol, names);

        // place in result and release path
        SET_VECTOR_ELT(result, i, path);
        UNPROTECT(1); // path
    }

    // ---- Append uniform grid graph exports ----
    // a) ugg_adj_list
    {
        const R_xlen_t nV = (R_xlen_t)uniform_grid_graph.adjacency_list.size();
        SEXP r_adj_list = PROTECT(Rf_allocVector(VECSXP, nV));
        for (R_xlen_t i = 0; i < nV; ++i) {
            const auto& neigh = uniform_grid_graph.adjacency_list[(size_t)i];
            const R_xlen_t deg = (R_xlen_t)neigh.size();
            SEXP r_adj = PROTECT(Rf_allocVector(INTSXP, deg));
            int* ip = INTEGER(r_adj);
            for (R_xlen_t j = 0; j < deg; ++j) ip[j] = (int)neigh[(size_t)j].first + 1;
            SET_VECTOR_ELT(r_adj_list, i, r_adj);
            UNPROTECT(1); // r_adj
        }
        SET_VECTOR_ELT(result, n_paths + 0, r_adj_list);
        UNPROTECT(1); // r_adj_list
    }

    // b) ugg_weight_list
    {
        const R_xlen_t nV = (R_xlen_t)uniform_grid_graph.adjacency_list.size();
        SEXP r_weight_list = PROTECT(Rf_allocVector(VECSXP, nV));
        for (R_xlen_t i = 0; i < nV; ++i) {
            const auto& neigh = uniform_grid_graph.adjacency_list[(size_t)i];
            const R_xlen_t deg = (R_xlen_t)neigh.size();
            SEXP r_w = PROTECT(Rf_allocVector(REALSXP, deg));
            double* wp = REAL(r_w);
            for (R_xlen_t j = 0; j < deg; ++j) wp[j] = neigh[(size_t)j].second;
            SET_VECTOR_ELT(r_weight_list, i, r_w);
            UNPROTECT(1); // r_w
        }
        SET_VECTOR_ELT(result, n_paths + 1, r_weight_list);
        UNPROTECT(1); // r_weight_list
    }

    // c) ugg_grid_vertices (1-based)
    {
        const R_xlen_t n = (R_xlen_t)uniform_grid_graph.grid_vertices.size();
        SEXP r_grid_vertices = PROTECT(Rf_allocVector(INTSXP, n));
        int* ip = INTEGER(r_grid_vertices);
        for (R_xlen_t i = 0; i < n; ++i) ip[i] = (int)uniform_grid_graph.grid_vertices[(size_t)i] + 1;
        SET_VECTOR_ELT(result, n_paths + 2, r_grid_vertices);
        UNPROTECT(1); // r_grid_vertices
    }

    // Release: result + result_names + names (path sublist names)
    UNPROTECT(3);
    return result;
}
