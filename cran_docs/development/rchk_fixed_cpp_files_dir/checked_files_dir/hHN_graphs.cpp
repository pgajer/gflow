/**
 * @brief Fixed version of hHN_graphs.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_create_hHN_graph
 * 
 * Issues fixed:
 * 1. S_create_hHN_graph (lines 387-389): unprotected 'res' during allocations
 * 2. Incorrect UNPROTECT order at lines 385 and 392
 * 
 * Changes made:
 * 1. Used container-first pattern - protect 'res' before other allocations
 * 2. Fixed PROTECT/UNPROTECT balance and order
 * 3. Ensured all allocations are properly protected
 */

#include <vector>
#include <utility>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Core function declaration (assumed available)
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> 
create_hHN_graph(const std::vector<std::vector<int>>& adj_vv,
                 const std::vector<std::vector<double>>& weight_vv,
                 int h);

extern "C" {

/**
 * Fixed version of S_create_hHN_graph
 * Fixes: unprotected 'res' during allocations at lines 387-389
 * Solution: Use container-first pattern, protect res before setting names
 */
SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_vv       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_vv = convert_weight_list_from_R(s_weight_list);

    // Coerce s_h to integer and extract its value
    PROTECT_INDEX ipx;
    SEXP sh = s_h;
    PROTECT_WITH_INDEX(sh, &ipx);
    if (TYPEOF(sh) != INTSXP) {
        REPROTECT(sh = Rf_coerceVector(sh, INTSXP), ipx);
    }
    int h = INTEGER(sh)[0];
    UNPROTECT(1); // sh

    // Call core function
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> hhn_graph = 
        create_hHN_graph(adj_vv, weight_vv, h);

    // Extract results
    std::vector<std::vector<int>> hhn_adj_vv       = hhn_graph.first;
    std::vector<std::vector<double>> hhn_weight_vv = hhn_graph.second;

    int n_vertices = static_cast<int>(hhn_adj_vv.size());

    // Create result list (container-first pattern)
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));

    // Build adjacency list
    {
        SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RA = PROTECT(Rf_allocVector(INTSXP, hhn_adj_vv[i].size()));
            int* A = INTEGER(RA);
            for (const auto& neighbor : hhn_adj_vv[i])
                *A++ = neighbor + 1;  // Adjust for R's 1-based indexing
            SET_VECTOR_ELT(adj_list, i, RA);
            UNPROTECT(1); // RA
        }
        SET_VECTOR_ELT(res, 0, adj_list);
        UNPROTECT(1); // adj_list
    }

    // Build distance/weight list
    {
        SEXP dist_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RD = PROTECT(Rf_allocVector(REALSXP, hhn_weight_vv[i].size()));
            double* D = REAL(RD);
            for (const auto& dist : hhn_weight_vv[i])
                *D++ = dist;
            SET_VECTOR_ELT(dist_list, i, RD);
            UNPROTECT(1); // RD
        }
        SET_VECTOR_ELT(res, 1, dist_list);
        UNPROTECT(1); // dist_list
    }

    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
        SET_STRING_ELT(names, 1, Rf_mkChar("dist_list"));
        Rf_setAttrib(res, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    UNPROTECT(1); // res
    return res;
}

} // extern "C"