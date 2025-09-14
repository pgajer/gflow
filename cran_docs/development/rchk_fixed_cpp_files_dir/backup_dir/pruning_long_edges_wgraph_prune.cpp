/**
 * @brief Fixed version of pruning_long_edges.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_wgraph_prune_long_edges
 * 
 * Issues fixed:
 * 1. S_wgraph_prune_long_edges (lines 1307, 1309, 1320, 1330, 1337): 
 *    unprotected pruned_adj_list, pruned_edge_lengths_list, R_path_lengths, R_edge_lengths during allocations
 * 
 * Changes made:
 * 1. Protected all intermediate allocations
 * 2. Used container-first pattern consistently
 * 3. Fixed all PROTECT/UNPROTECT imbalances
 * 4. Used literal UNPROTECT values
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

/**
 * Fixed version of S_wgraph_prune_long_edges
 * Fixes: unprotected allocations at lines 1307, 1309, 1320, 1330, 1337
 * Solution: Protect all allocations, use container-first pattern
 */
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
    
    // Extract scalar parameters
    double length_threshold = Rf_asReal(length_threshold_sexp);
    double ratio_threshold = Rf_asReal(ratio_threshold_sexp);
    int max_path_length = Rf_asInteger(max_path_length_sexp);
    bool preserve_connectivity = (Rf_asLogical(preserve_connectivity_sexp) == TRUE);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Validate parameters
    if (length_threshold <= 0) {
        Rf_error("length_threshold must be positive");
    }
    if (ratio_threshold <= 1.0) {
        Rf_error("ratio_threshold must be greater than 1.0");
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
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 10;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: pruned_adj_list
    {
        SEXP pruned_adj_list = PROTECT(convert_vector_vector_int_to_R(result_data.pruned_adj_list));
        SET_VECTOR_ELT(result, 0, pruned_adj_list);
        UNPROTECT(1);
    }
    
    // 1: pruned_edge_lengths_list
    {
        SEXP pruned_edge_lengths_list = PROTECT(convert_vector_vector_double_to_R(result_data.pruned_edge_lengths));
        SET_VECTOR_ELT(result, 1, pruned_edge_lengths_list);
        UNPROTECT(1);
    }
    
    // 2: removed_edge_lengths
    {
        SEXP R_edge_lengths = PROTECT(convert_vector_double_to_R(result_data.removed_edge_lengths));
        SET_VECTOR_ELT(result, 2, R_edge_lengths);
        UNPROTECT(1);
    }
    
    // 3: alt_path_lengths
    {
        SEXP R_path_lengths = PROTECT(convert_vector_double_to_R(result_data.alt_path_lengths));
        SET_VECTOR_ELT(result, 3, R_path_lengths);
        UNPROTECT(1);
    }
    
    // 4: removed_edges_source (convert to 1-based for R)
    {
        SEXP sources = PROTECT(Rf_allocVector(INTSXP, result_data.removed_edges_source.size()));
        int* src_ptr = INTEGER(sources);
        for (size_t i = 0; i < result_data.removed_edges_source.size(); i++) {
            src_ptr[i] = result_data.removed_edges_source[i] + 1; // Convert to 1-based
        }
        SET_VECTOR_ELT(result, 4, sources);
        UNPROTECT(1);
    }
    
    // 5: removed_edges_target (convert to 1-based for R)
    {
        SEXP targets = PROTECT(Rf_allocVector(INTSXP, result_data.removed_edges_target.size()));
        int* tgt_ptr = INTEGER(targets);
        for (size_t i = 0; i < result_data.removed_edges_target.size(); i++) {
            tgt_ptr[i] = result_data.removed_edges_target[i] + 1; // Convert to 1-based
        }
        SET_VECTOR_ELT(result, 5, targets);
        UNPROTECT(1);
    }
    
    // 6: n_edges_removed
    {
        SEXP n_removed = PROTECT(Rf_ScalarInteger(result_data.n_edges_removed));
        SET_VECTOR_ELT(result, 6, n_removed);
        UNPROTECT(1);
    }
    
    // 7: n_edges_remaining
    {
        SEXP n_remaining = PROTECT(Rf_ScalarInteger(result_data.n_edges_remaining));
        SET_VECTOR_ELT(result, 7, n_remaining);
        UNPROTECT(1);
    }
    
    // 8: total_length_removed
    {
        SEXP total_removed = PROTECT(Rf_ScalarReal(result_data.total_length_removed));
        SET_VECTOR_ELT(result, 8, total_removed);
        UNPROTECT(1);
    }
    
    // 9: pruning_threshold
    {
        SEXP threshold = PROTECT(Rf_ScalarReal(result_data.pruning_threshold));
        SET_VECTOR_ELT(result, 9, threshold);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
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
        UNPROTECT(1); // names
    }
    
    UNPROTECT(1); // result
    return result;
}

} // extern "C"