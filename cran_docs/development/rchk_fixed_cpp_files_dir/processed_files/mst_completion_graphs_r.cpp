/**
 * @brief Fixed version of mst_completion_graphs_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_create_mst_completion_graph
 * 
 * Issues fixed:
 * 1. S_create_mst_completion_graph (line 216): UNPROTECT(variable) - nprot
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constant
 * 2. Used container-first pattern consistently
 * 3. Fixed PROTECT/UNPROTECT balance
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Macro for error reporting (assumed available)
#ifndef REPORT_ERROR
#define REPORT_ERROR(msg) Rf_error(msg)
#endif

// Forward declarations for helper functions (assumed available)
extern std::unique_ptr<std::vector<std::vector<double>>> Rmatrix_to_cpp(SEXP);

// Core function structures and declarations
struct mst_completion_graph_t {
    struct graph_t {
        // Graph structure details
    } mstree, completed_mstree;
    std::vector<double> mstree_edge_weights;
};

// Core computation function (assumed available)
mst_completion_graph_t create_mst_completion_graph(
    const std::vector<std::vector<double>>& X,
    double q_thld,
    bool verbose);

// Helper function (assumed available)
std::pair<SEXP, SEXP> create_r_graph_from_set_wgraph(const mst_completion_graph_t::graph_t& graph);

extern "C" {

/**
 * Fixed version of S_create_mst_completion_graph
 * Fixes: UNPROTECT(variable) at line 216
 * Solution: Use container-first pattern, use literal UNPROTECT
 */
SEXP S_create_mst_completion_graph(
    SEXP s_X,
    SEXP s_q_thld,
    SEXP s_verbose
) {
    // Input checking
    if (!Rf_isReal(s_X) || !Rf_isMatrix(s_X)) {
        REPORT_ERROR("X must be a numeric matrix");
    }
    if (!Rf_isReal(s_q_thld) || Rf_length(s_q_thld) != 1) {
        REPORT_ERROR("q_thld must be a numeric scalar");
    }
    if (!Rf_isLogical(s_verbose) || Rf_length(s_verbose) != 1) {
        REPORT_ERROR("verbose must be a logical scalar");
    }

    double q_thld = REAL(s_q_thld)[0];
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));

    // Call core computation
    mst_completion_graph_t res = create_mst_completion_graph(X, q_thld, verbose);

    // Create R return list (container-first pattern)
    SEXP r_list = PROTECT(Rf_allocVector(VECSXP, 5));

    // Get MST graph components
    auto [r_mst_adj_list, r_mst_weights_list] = create_r_graph_from_set_wgraph(res.mstree);
    PROTECT(r_mst_adj_list);
    PROTECT(r_mst_weights_list);
    SET_VECTOR_ELT(r_list, 0, r_mst_adj_list);
    SET_VECTOR_ELT(r_list, 1, r_mst_weights_list);
    UNPROTECT(2); // r_mst_adj_list, r_mst_weights_list

    // Get completed MST graph components
    auto [r_cmst_adj_list, r_cmst_weights_list] = create_r_graph_from_set_wgraph(res.completed_mstree);
    PROTECT(r_cmst_adj_list);
    PROTECT(r_cmst_weights_list);
    SET_VECTOR_ELT(r_list, 2, r_cmst_adj_list);
    SET_VECTOR_ELT(r_list, 3, r_cmst_weights_list);
    UNPROTECT(2); // r_cmst_adj_list, r_cmst_weights_list

    // MST edge weights
    {
        size_t m = res.mstree_edge_weights.size();
        SEXP r_mst_weights = PROTECT(Rf_allocVector(REALSXP, m));
        std::copy(res.mstree_edge_weights.begin(), res.mstree_edge_weights.end(), REAL(r_mst_weights));
        SET_VECTOR_ELT(r_list, 4, r_mst_weights);
        UNPROTECT(1); // r_mst_weights
    }

    // Set names
    {
        SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 5));
        SET_STRING_ELT(r_names, 0, Rf_mkChar("mst_adj_list"));
        SET_STRING_ELT(r_names, 1, Rf_mkChar("mst_weight_list"));
        SET_STRING_ELT(r_names, 2, Rf_mkChar("cmst_adj_list"));
        SET_STRING_ELT(r_names, 3, Rf_mkChar("cmst_weight_list"));
        SET_STRING_ELT(r_names, 4, Rf_mkChar("mst_edge_weights"));
        Rf_setAttrib(r_list, R_NamesSymbol, r_names);
        UNPROTECT(1); // r_names
    }

    UNPROTECT(1); // r_list
    return r_list;
}

} // extern "C"