/**
 * @brief Hardened rchk-safe drop-in of S_wgraph_prune_long_edges
 *
 * Changes vs. provided version:
 * - Keep `result` and `names` protected until the end; finish with UNPROTECT(2).
 * - Maintain container-first local PROTECT/UNPROTECT(1) for each element.
 * - Defensive scalars via Rf_as*; validate before allocations.
 * - Use R_xlen_t where allocating from size_t.
 */
#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);

// Core structures and declarations
struct pruned_graph_result_t {
    std::vector<std::vector<int>> pruned_adj_list;
    std::vector<std::vector<double>> pruned_edge_lengths;
    std::vector<double> removed_edge_lengths;
    std::vector<double> alt_path_lengths;
    std::vector<int> removed_edges_source;
    std::vector<int> removed_edges_target;
    int n_edges_removed;
    int n_edges_remaining;
    double total_length_removed;
    double pruning_threshold;
};

// Core computation function (assumed available)
pruned_graph_result_t wgraph_prune_long_edges(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    double length_threshold,
    double ratio_threshold,
    int max_path_length,
    bool preserve_connectivity,
    bool verbose);

extern "C" {

SEXP S_wgraph_prune_long_edges(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP length_threshold_sexp,
    SEXP ratio_threshold_sexp,
    SEXP max_path_length_sexp,
    SEXP preserve_connectivity_sexp,
    SEXP verbose_sexp
) {
    // Convert inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    // Extract scalar parameters (defensive)
    const double length_threshold = Rf_asReal(length_threshold_sexp);
    const double ratio_threshold  = Rf_asReal(ratio_threshold_sexp);
    const int    max_path_length  = Rf_asInteger(max_path_length_sexp);
    const bool   preserve_connectivity = (Rf_asLogical(preserve_connectivity_sexp) == TRUE);
    const bool   verbose               = (Rf_asLogical(verbose_sexp) == TRUE);

    // Validate parameters before any PROTECTs
    if (!(length_threshold > 0.0)) {
        Rf_error("S_wgraph_prune_long_edges: 'length_threshold' must be positive");
    }
    if (!(ratio_threshold > 1.0)) {
        Rf_error("S_wgraph_prune_long_edges: 'ratio_threshold' must be greater than 1.0");
    }

    // Core computation
    pruned_graph_result_t result_data = wgraph_prune_long_edges(
        adj_list,
        weight_list,
        length_threshold,
        ratio_threshold,
        max_path_length,
        preserve_connectivity,
        verbose
    );

    // Create result container
    const int N = 10;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));

    // 0: pruned_adj_list
    {
        SEXP s0 = PROTECT(convert_vector_vector_int_to_R(result_data.pruned_adj_list));
        SET_VECTOR_ELT(result, 0, s0);
        UNPROTECT(1);
    }
    // 1: pruned_edge_lengths
    {
        SEXP s1 = PROTECT(convert_vector_vector_double_to_R(result_data.pruned_edge_lengths));
        SET_VECTOR_ELT(result, 1, s1);
        UNPROTECT(1);
    }
    // 2: removed_edge_lengths
    {
        SEXP s2 = PROTECT(convert_vector_double_to_R(result_data.removed_edge_lengths));
        SET_VECTOR_ELT(result, 2, s2);
        UNPROTECT(1);
    }
    // 3: alt_path_lengths
    {
        SEXP s3 = PROTECT(convert_vector_double_to_R(result_data.alt_path_lengths));
        SET_VECTOR_ELT(result, 3, s3);
        UNPROTECT(1);
    }
    // 4: removed_edges_source (1-based)
    {
        const R_xlen_t m = (R_xlen_t) result_data.removed_edges_source.size();
        SEXP s4 = PROTECT(Rf_allocVector(INTSXP, m));
        int* p = INTEGER(s4);
        for (R_xlen_t i = 0; i < m; ++i) p[i] = result_data.removed_edges_source[(size_t)i] + 1;
        SET_VECTOR_ELT(result, 4, s4);
        UNPROTECT(1);
    }
    // 5: removed_edges_target (1-based)
    {
        const R_xlen_t m = (R_xlen_t) result_data.removed_edges_target.size();
        SEXP s5 = PROTECT(Rf_allocVector(INTSXP, m));
        int* p = INTEGER(s5);
        for (R_xlen_t i = 0; i < m; ++i) p[i] = result_data.removed_edges_target[(size_t)i] + 1;
        SET_VECTOR_ELT(result, 5, s5);
        UNPROTECT(1);
    }
    // 6: n_edges_removed
    {
        SEXP s6 = PROTECT(Rf_ScalarInteger(result_data.n_edges_removed));
        SET_VECTOR_ELT(result, 6, s6);
        UNPROTECT(1);
    }
    // 7: n_edges_remaining
    {
        SEXP s7 = PROTECT(Rf_ScalarInteger(result_data.n_edges_remaining));
        SET_VECTOR_ELT(result, 7, s7);
        UNPROTECT(1);
    }
    // 8: total_length_removed
    {
        SEXP s8 = PROTECT(Rf_ScalarReal(result_data.total_length_removed));
        SET_VECTOR_ELT(result, 8, s8);
        UNPROTECT(1);
    }
    // 9: pruning_threshold
    {
        SEXP s9 = PROTECT(Rf_ScalarReal(result_data.pruning_threshold));
        SET_VECTOR_ELT(result, 9, s9);
        UNPROTECT(1);
    }

    // names — keep protected until tail
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("pruned_edge_lengths"));
    SET_STRING_ELT(names, 2, Rf_mkChar("removed_edge_lengths"));
    SET_STRING_ELT(names, 3, Rf_mkChar("alt_path_lengths"));
    SET_STRING_ELT(names, 4, Rf_mkChar("removed_edges_source"));
    SET_STRING_ELT(names, 5, Rf_mkChar("removed_edges_target"));
    SET_STRING_ELT(names, 6, Rf_mkChar("n_edges_removed"));
    SET_STRING_ELT(names, 7, Rf_mkChar("n_edges_remaining"));
    SET_STRING_ELT(names, 8, Rf_mkChar("total_length_removed"));
    SET_STRING_ELT(names, 9, Rf_mkChar("pruning_threshold"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
