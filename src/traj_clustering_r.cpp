/**
 * @file traj_clustering_r.cpp
 * @brief SEXP interface for trajectory-graph construction within a single cell.
 *
 * This file provides an R-callable entry point to build a kNN-sparsified
 * trajectory similarity graph using (IDF-weighted) Jaccard overlap.
 *
 * Leiden clustering is intentionally not implemented in C++ at this time;
 * the R wrapper should run Leiden via igraph.
 */

#include <R.h>
#include <Rinternals.h>

#include "traj_clustering.hpp"

#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static std::string as_std_string(SEXP s) {
    if (TYPEOF(s) != STRSXP || Rf_length(s) < 1) {
        throw std::runtime_error("Expected a character scalar");
    }
    SEXP elt = STRING_ELT(s, 0);
    if (elt == NA_STRING) {
        throw std::runtime_error("Character scalar is NA");
    }
    return std::string(CHAR(elt));
}

static int as_int_scalar(SEXP s, const char* name) {
    if (!Rf_isInteger(s) || Rf_length(s) != 1) {
        throw std::runtime_error(std::string("Expected integer scalar for ") + name);
    }
    return INTEGER(s)[0];
}

static double as_real_scalar(SEXP s, const char* name) {
    if (!Rf_isReal(s) || Rf_length(s) != 1) {
        throw std::runtime_error(std::string("Expected numeric scalar for ") + name);
    }
    return REAL(s)[0];
}

static bool as_bool_scalar(SEXP s, const char* name) {
    if (!Rf_isLogical(s) || Rf_length(s) != 1) {
        throw std::runtime_error(std::string("Expected logical scalar for ") + name);
    }
    return LOGICAL(s)[0] == TRUE;
}

static std::vector<std::vector<int>> convert_traj_list_from_R(SEXP s_traj_list) {
    if (TYPEOF(s_traj_list) != VECSXP) {
        throw std::runtime_error("s_traj_list must be a list of integer vectors");
    }

    const R_len_t n = Rf_length(s_traj_list);
    std::vector<std::vector<int>> out;
    out.resize(static_cast<size_t>(n));

    for (R_len_t i = 0; i < n; ++i) {
        SEXP s_traj = VECTOR_ELT(s_traj_list, i);
        if (!Rf_isInteger(s_traj)) {
            throw std::runtime_error("Each trajectory must be an integer vector");
        }
        const R_len_t len = Rf_length(s_traj);
        const int* p = INTEGER(s_traj);
        out[static_cast<size_t>(i)].assign(p, p + len);
    }
    return out;
}

// -----------------------------------------------------------------------------
// R entry point
// -----------------------------------------------------------------------------

/**
 * @brief Build a kNN-sparsified trajectory similarity graph for clustering.
 *
 * @param s_traj_list List of integer vectors (0-based vertex indices).
 * @param s_similarity_type "idf.jaccard" or "jaccard".
 * @param s_overlap_mode "any", "min.shared", "min.frac", "min.subpath".
 * @param s_min_shared Integer scalar.
 * @param s_min_frac Numeric scalar.
 * @param s_min_subpath Integer scalar.
 * @param s_exclude_endpoints Logical scalar.
 * @param s_idf_smooth Numeric scalar.
 * @param s_k Integer scalar (kNN).
 * @param s_symmetrize "mutual", "union", "none".
 * @param s_knn_select "weighted" or "raw".
 * @param s_n_threads Integer scalar.
 *
 * @return R list with edge_from, edge_to (1-based trajectory indices) and edge_weight.
 */
extern "C" SEXP S_cluster_cell_trajectories(SEXP s_traj_list,
                                            SEXP s_similarity_type,
                                            SEXP s_overlap_mode,
                                            SEXP s_min_shared,
                                            SEXP s_min_frac,
                                            SEXP s_min_subpath,
                                            SEXP s_exclude_endpoints,
                                            SEXP s_idf_smooth,
                                            SEXP s_k,
                                            SEXP s_symmetrize,
                                            SEXP s_knn_select,
                                            SEXP s_n_threads) {

    try {
        const std::vector<std::vector<int>> traj_list = convert_traj_list_from_R(s_traj_list);

        const std::string similarity_type = as_std_string(s_similarity_type);
        const std::string overlap_mode = as_std_string(s_overlap_mode);
        const std::string symmetrize = as_std_string(s_symmetrize);
        const std::string knn_select = as_std_string(s_knn_select);

        traj_similarity_spec_t spec;
        if (similarity_type == "idf.jaccard") {
            spec.similarity_type = traj_similarity_spec_t::similarity_t::idf_jaccard;
        } else if (similarity_type == "jaccard") {
            spec.similarity_type = traj_similarity_spec_t::similarity_t::jaccard;
        } else {
            throw std::runtime_error("similarity_type must be 'idf.jaccard' or 'jaccard'");
        }

        if (overlap_mode == "any") {
            spec.overlap_mode = traj_similarity_spec_t::overlap_mode_t::any;
        } else if (overlap_mode == "min.shared") {
            spec.overlap_mode = traj_similarity_spec_t::overlap_mode_t::min_shared;
        } else if (overlap_mode == "min.frac") {
            spec.overlap_mode = traj_similarity_spec_t::overlap_mode_t::min_frac;
        } else if (overlap_mode == "min.subpath") {
            spec.overlap_mode = traj_similarity_spec_t::overlap_mode_t::min_subpath;
        } else {
            throw std::runtime_error("overlap_mode must be one of: any, min.shared, min.frac, min.subpath");
        }

        spec.min_shared = as_int_scalar(s_min_shared, "min_shared");
        spec.min_frac = as_real_scalar(s_min_frac, "min_frac");
        spec.min_subpath = as_int_scalar(s_min_subpath, "min_subpath");
        spec.exclude_endpoints = as_bool_scalar(s_exclude_endpoints, "exclude_endpoints");
        spec.idf_smooth = as_real_scalar(s_idf_smooth, "idf_smooth");

        const int k = as_int_scalar(s_k, "k");
        const int n_threads = as_int_scalar(s_n_threads, "n_threads");

        traj_graph_t g = build_traj_knn_graph(traj_list, spec, k, symmetrize, knn_select, n_threads);

        const R_len_t n_edges = static_cast<R_len_t>(g.from.size());

        SEXP s_out = PROTECT(Rf_allocVector(VECSXP, 4));
        SEXP s_names = PROTECT(Rf_allocVector(STRSXP, 4));

        SET_STRING_ELT(s_names, 0, Rf_mkChar("n.traj"));
        SET_STRING_ELT(s_names, 1, Rf_mkChar("edge.from"));
        SET_STRING_ELT(s_names, 2, Rf_mkChar("edge.to"));
        SET_STRING_ELT(s_names, 3, Rf_mkChar("edge.weight"));

        SET_VECTOR_ELT(s_out, 0, Rf_ScalarInteger(g.n_traj));

        SEXP s_from = PROTECT(Rf_allocVector(INTSXP, n_edges));
        SEXP s_to = PROTECT(Rf_allocVector(INTSXP, n_edges));
        SEXP s_w = PROTECT(Rf_allocVector(REALSXP, n_edges));

        int* p_from = INTEGER(s_from);
        int* p_to = INTEGER(s_to);
        double* p_w = REAL(s_w);

        for (R_len_t e = 0; e < n_edges; ++e) {
            p_from[e] = g.from[static_cast<size_t>(e)] + 1; // to 1-based
            p_to[e] = g.to[static_cast<size_t>(e)] + 1;
            p_w[e] = g.weight[static_cast<size_t>(e)];
        }

        SET_VECTOR_ELT(s_out, 1, s_from);
        SET_VECTOR_ELT(s_out, 2, s_to);
        SET_VECTOR_ELT(s_out, 3, s_w);

        Rf_setAttrib(s_out, R_NamesSymbol, s_names);

        UNPROTECT(5);
        return s_out;

    } catch (const std::exception& e) {
        Rf_error("S_cluster_cell_trajectories: %s", e.what());
    } catch (...) {
        Rf_error("S_cluster_cell_trajectories: unknown error");
    }

    return R_NilValue;
}
