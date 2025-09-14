/**
 * @file set_wgraph_rchk_hardened.cpp
 * @brief rchk-hardened S_... wrappers for set_wgraph, with defensive checks.
 *
 * This file provides hardened versions of the three exported S_... functions:
 *  - S_find_graph_paths_within_radius
 *  - S_compute_edge_weight_rel_deviations
 *  - S_remove_redundant_edges
 *
 * Design:
 *  - Container-first PROTECT pattern.
 *  - Only literal UNPROTECT counts.
 *  - No PROTECTs inside lambdas.
 *  - Defensive type/NA/finiteness checks for scalar inputs.
 *  - Bounds checks for indices derived from R inputs.
 *  - R_xlen_t used for alloc sizes (long vectors).
 */

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <cmath>

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

static inline void ensure_integerish(SEXP x, const char* name) {
    if (!(TYPEOF(x) == INTSXP || TYPEOF(x) == REALSXP))
        Rf_error("%s must be integer-like", name);
    if (TYPEOF(x) == INTSXP) {
        if (Rf_length(x) != 1) Rf_error("%s must be a scalar", name);
        if (INTEGER(x)[0] == NA_INTEGER) Rf_error("%s is NA", name);
    } else { // REALSXP
        if (Rf_length(x) != 1) Rf_error("%s must be a scalar", name);
        double v = REAL(x)[0];
        if (!R_finite(v)) Rf_error("%s must be finite", name);
        // Require it to be close to an integer to avoid surprises
        if (std::floor(v) != v) Rf_error("%s must be an integer value", name);
    }
}

static inline void ensure_numeric_scalar(SEXP x, const char* name) {
    if (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP)
        Rf_error("%s must be numeric", name);
    if (Rf_length(x) != 1) Rf_error("%s must be a scalar", name);
    double v = (TYPEOF(x) == REALSXP) ? REAL(x)[0] : (double)INTEGER(x)[0];
    if (!R_finite(v)) Rf_error("%s must be finite", name);
}

extern "C" {

SEXP S_find_graph_paths_within_radius(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_start,
    SEXP s_radius)
{
    // Convert adjacency/weights (assume converters validate structural invariants)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    // Validate and extract start (1-based) and radius
    ensure_integerish(s_start, "start");
    ensure_numeric_scalar(s_radius, "radius");
    const int start_1based_i = Rf_asInteger(s_start);
    if (start_1based_i <= 0) Rf_error("start must be >= 1");
    const size_t nV = graph.adjacency_list.size();
    if ((size_t)start_1based_i > nV)
        Rf_error("start out of range [1, %zu]", nV);

    const size_t start_vertex = (size_t)start_1based_i - 1;
    const double radius = Rf_asReal(s_radius);
    if (radius < 0.0) Rf_error("radius must be non-negative");

    // Compute
    shortest_paths_t shortest_paths = graph.find_graph_paths_within_radius(start_vertex, radius);

    // Build return list: c(paths, reachable_vertices, vertex_to_path_map)
    const R_xlen_t n_elements = 3;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));          // [1]

    // 0) paths : list<int>
    {
        const R_xlen_t n_paths = (R_xlen_t)shortest_paths.paths.size();
        SEXP r_paths_list = PROTECT(Rf_allocVector(VECSXP, n_paths));   // [2]
        for (R_xlen_t i = 0; i < n_paths; ++i) {
            const auto& verts = shortest_paths.paths[(size_t)i].vertices;
            const R_xlen_t m = (R_xlen_t)verts.size();
            SEXP r_path = PROTECT(Rf_allocVector(INTSXP, m));           // [3]
            int* p = INTEGER(r_path);
            for (R_xlen_t j = 0; j < m; ++j) p[j] = (int)verts[(size_t)j] + 1;
            SET_VECTOR_ELT(r_paths_list, i, r_path);
            UNPROTECT(1); // r_path                                    // [-1] -> [2]
        }
        SET_VECTOR_ELT(result, 0, r_paths_list);
        UNPROTECT(1); // r_paths_list                                   // [-1] -> [1]
    }

    // 1) reachable_vertices : int vector
    {
        const auto& reachable = shortest_paths.reachable_vertices;
        const R_xlen_t m = (R_xlen_t)reachable.size();
        SEXP r_reach = PROTECT(Rf_allocVector(INTSXP, m));              // [2]
        int* p = INTEGER(r_reach);
        R_xlen_t idx = 0;
        for (const auto& v : reachable) p[idx++] = (int)v + 1;
        SET_VECTOR_ELT(result, 1, r_reach);
        UNPROTECT(1); // r_reach                                        // [-1] -> [1]
    }

    // 2) vertex_to_path_map : numeric matrix (nrow x 3) or NULL
    {
        const auto& map = shortest_paths.vertex_to_path_map;
        if (!map.empty()) {
            const R_xlen_t nrow = (R_xlen_t)map.size();
            const R_xlen_t ncol = 3;
            SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));  // [2]
            double* M = REAL(r_mat);
            R_xlen_t i = 0;
            for (const auto& kv : map) {
                const size_t vertex = kv.first;
                const auto& sub    = kv.second;
                M[i]          = (double)vertex + 1.0;
                M[i + nrow]   = (double)sub.path_idx + 1.0;
                M[i + 2*nrow] = (double)sub.vertex_idx + 1.0;
                ++i;
            }
            SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol));      // [3]
            SET_STRING_ELT(colnames, 0, Rf_mkChar("vertex"));
            SET_STRING_ELT(colnames, 1, Rf_mkChar("path_idx"));
            SET_STRING_ELT(colnames, 2, Rf_mkChar("vertex_idx"));
            SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));         // [4]
            SET_VECTOR_ELT(dimnames, 0, R_NilValue);
            SET_VECTOR_ELT(dimnames, 1, colnames);
            Rf_setAttrib(r_mat, R_DimNamesSymbol, dimnames);
            SET_VECTOR_ELT(result, 2, r_mat);
            UNPROTECT(3); // r_mat, colnames, dimnames                   // [-3] -> [1]
        } else {
            SET_VECTOR_ELT(result, 2, R_NilValue);
        }
    }

    // names(result)
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, n_elements));              // [2]
    SET_STRING_ELT(nm, 0, Rf_mkChar("paths"));
    SET_STRING_ELT(nm, 1, Rf_mkChar("reachable_vertices"));
    SET_STRING_ELT(nm, 2, Rf_mkChar("vertex_to_path_map"));
    Rf_setAttrib(result, R_NamesSymbol, nm);

    UNPROTECT(2); // result, nm
    return result;
}

SEXP S_compute_edge_weight_rel_deviations(SEXP s_adj_list, SEXP s_weight_list)
{
    std::vector<std::vector<int>>    adj_vect    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_vect, weight_vect);

    // Compute
    auto results = graph.compute_edge_weight_rel_deviations();

    const R_xlen_t n_edges = (R_xlen_t)results.size();
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n_edges, 4));            // [1]
    double* X = REAL(out);

    for (R_xlen_t i = 0; i < n_edges; ++i) {
        const auto& r = results[(size_t)i];
        X[i]              = (double)r.source + 1.0;
        X[i + n_edges]    = (double)r.target + 1.0;
        X[i + 2*n_edges]  = (double)r.best_intermediate + 1.0;
        X[i + 3*n_edges]  = r.rel_deviation;
    }

    SEXP col = PROTECT(Rf_allocVector(STRSXP, 4));                      // [2]
    SET_STRING_ELT(col, 0, Rf_mkChar("source"));
    SET_STRING_ELT(col, 1, Rf_mkChar("target"));
    SET_STRING_ELT(col, 2, Rf_mkChar("best_intermediate"));
    SET_STRING_ELT(col, 3, Rf_mkChar("rel_deviation"));

    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));                       // [3]
    SET_VECTOR_ELT(dn, 0, R_NilValue);
    SET_VECTOR_ELT(dn, 1, col);
    Rf_setAttrib(out, R_DimNamesSymbol, dn);

    UNPROTECT(3); // out, col, dn
    return out;
}

SEXP S_remove_redundant_edges(SEXP s_adj_list, SEXP s_weight_list)
{
    std::vector<std::vector<int>>    adj_vect    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vect = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_vect, weight_vect);

    const auto rel_devs = graph.compute_edge_weight_rel_deviations();
    constexpr double EPS = 1e-16;

    for (size_t i = 0; i < rel_devs.size(); ++i) {
        const auto& e = rel_devs[i];
        if (e.rel_deviation < EPS) {
            const size_t u = e.source;
            const size_t v = e.target;
            const size_t w = e.best_intermediate;
            // Note: Rf_error below occurs before any PROTECT in this function,
            // so no PROTECT imbalance is possible on these paths.
            if (graph.adjacency_list[u].size() <= 1) {
                Rf_error("Cannot remove edge (%zu,%zu) via %zu: vertex %zu would become isolated",
                         u + 1, v + 1, w + 1, u + 1);
            }
            if (graph.adjacency_list[v].size() <= 1) {
                Rf_error("Cannot remove edge (%zu,%zu) via %zu: vertex %zu would become isolated",
                         u + 1, v + 1, w + 1, v + 1);
            }
            graph.remove_edge(u, v);
        }
    }

    // Build return list: list(adj_list, weight_list)
    const R_xlen_t N = 2;
    SEXP res = PROTECT(Rf_allocVector(VECSXP, N));                      // [1]

    const R_xlen_t nV = (R_xlen_t)graph.adjacency_list.size();

    // adj_list
    {
        SEXP adj_out = PROTECT(Rf_allocVector(VECSXP, nV));             // [2]
        for (R_xlen_t i = 0; i < nV; ++i) {
            const auto& nbrs = graph.adjacency_list[(size_t)i];
            const R_xlen_t m = (R_xlen_t)nbrs.size();
            SEXP vi = PROTECT(Rf_allocVector(INTSXP, m));               // [3]
            int* p = INTEGER(vi);
            R_xlen_t k = 0;
            for (const auto& pr : nbrs) p[k++] = (int)pr.first + 1;
            SET_VECTOR_ELT(adj_out, i, vi);
            UNPROTECT(1); // vi                                          // [-1] -> [2]
        }
        SET_VECTOR_ELT(res, 0, adj_out);
        UNPROTECT(1); // adj_out                                         // [-1] -> [1]
    }

    // weight_list
    {
        SEXP w_out = PROTECT(Rf_allocVector(VECSXP, nV));               // [2]
        for (R_xlen_t i = 0; i < nV; ++i) {
            const auto& nbrs = graph.adjacency_list[(size_t)i];
            const R_xlen_t m = (R_xlen_t)nbrs.size();
            SEXP wi = PROTECT(Rf_allocVector(REALSXP, m));              // [3]
            double* p = REAL(wi);
            R_xlen_t k = 0;
            for (const auto& pr : nbrs) p[k++] = pr.second;
            SET_VECTOR_ELT(w_out, i, wi);
            UNPROTECT(1); // wi                                          // [-1] -> [2]
        }
        SET_VECTOR_ELT(res, 1, w_out);
        UNPROTECT(1); // w_out                                           // [-1] -> [1]
    }

    SEXP nm = PROTECT(Rf_allocVector(STRSXP, N));                       // [2]
    SET_STRING_ELT(nm, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(nm, 1, Rf_mkChar("weight_list"));
    Rf_setAttrib(res, R_NamesSymbol, nm);

    UNPROTECT(2); // res, nm
    return res;
}

} // extern "C"
