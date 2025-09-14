/**
 * @brief Fixed versions of pruning_long_edges.cpp functions to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected versions of:
 * - S_find_shortest_alt_path
 * - S_wgraph_prune_long_edges
 * 
 * Issues fixed:
 * 1. S_find_shortest_alt_path (lines 794, 797): Negative depth and over-unprotect
 * 2. S_wgraph_prune_long_edges (lines 1307-1337): Unprotected allocations
 * 
 * Changes made:
 * 1. Fixed unbalanced PROTECT/UNPROTECT pairs
 * 2. Protected all allocations immediately
 * 3. Used container-first pattern consistently
 * 4. Ensured all UNPROTECTs use literal constants
 */

#include <vector>
#include <queue>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);

// Include necessary headers
#include "pruning_long_edges.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Helper function assumed to be available
std::vector<std::vector<std::pair<int, int>>> convert_to_int_weighted_adj_list(
    const std::vector<std::vector<int>>& adj_vect,
    const std::vector<std::vector<int>>& weight_vect);

// Core function declarations
std::vector<int> find_shortest_alt_path(
    const std::vector<std::vector<std::pair<int, int>>>& graph,
    int source, int target, int edge_isize);

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>
graph_prune_long_edges(const std::vector<std::vector<int>>& adj_vect,
                       const std::vector<std::vector<double>>& edge_length_vect,
                       std::vector<double>& path_lengths,
                       std::vector<double>& edge_lengths,
                       double alt_path_len_ratio_thld,
                       bool use_total_length_constraint,
                       bool verbose);

extern "C" {

/**
 * Fixed version of S_find_shortest_alt_path
 * Fixes: Negative depth at lines 794, 797 - UNPROTECT without matching PROTECT
 */
SEXP S_find_shortest_alt_path(SEXP s_adj_list,
                              SEXP s_isize_list,
                              SEXP s_source,
                              SEXP s_target,
                              SEXP s_edge_isize) {
    
    // Convert inputs (no PROTECT needed for conversion helpers)
    std::vector<std::vector<int>> adj_vect   = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<int>> isize_vect = convert_adj_list_from_R(s_isize_list);
    
    // Extract scalar parameters using defensive coercion
    int source = Rf_asInteger(s_source);
    int target = Rf_asInteger(s_target);
    int edge_isize = Rf_asInteger(s_edge_isize);
    
    // Core computation (no R allocations inside)
    std::vector<std::vector<std::pair<int, int>>> iigraph = 
        convert_to_int_weighted_adj_list(adj_vect, isize_vect);
    
    std::vector<int> path = find_shortest_alt_path(iigraph, source, target, edge_isize);
    
    // Convert result to R (convert_vector_int_to_R returns unprotected SEXP)
    // We need to check if this function PROTECTs internally
    SEXP s_path = PROTECT(convert_vector_int_to_R(path));
    UNPROTECT(1);
    
    return s_path;
}

/**
 * Fixed version of S_wgraph_prune_long_edges
 * Fixes: Unprotected allocations at lines 1307-1337
 * The original had unprotected pruned_adj_list and pruned_edge_lengths_list
 * while calling allocating functions
 */
SEXP S_wgraph_prune_long_edges(SEXP s_adj_list,
                               SEXP s_edge_length_list,
                               SEXP s_alt_path_len_ratio_thld,
                               SEXP s_use_total_length_constraint,
                               SEXP s_verbose) {

    // Convert inputs
    std::vector<std::vector<int>> adj_vect            = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> edge_length_vect = convert_weight_list_from_R(s_edge_length_list);

    // Extract scalar parameters
    double alt_path_len_ratio_thld = Rf_asReal(s_alt_path_len_ratio_thld);
    bool use_total_length_constraint = (Rf_asLogical(s_use_total_length_constraint) == TRUE);
    bool verbose = (Rf_asLogical(s_verbose) == TRUE);

    // Core computation
    std::vector<double> path_lengths;
    std::vector<double> edge_lengths;

    auto pruning_res = graph_prune_long_edges(adj_vect,
                                               edge_length_vect,
                                               path_lengths,
                                               edge_lengths,
                                               alt_path_len_ratio_thld,
                                               use_total_length_constraint,
                                               verbose);

    auto pruned_adj_vect          = pruning_res.first;
    auto pruned_edge_lengths_vect = pruning_res.second;

    // Build result (container-first pattern)
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 4));
    
    int n_vertices = static_cast<int>(pruned_adj_vect.size());
    
    // 0: pruned_adj_list
    {
        SEXP pruned_adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RA = PROTECT(Rf_allocVector(INTSXP, pruned_adj_vect[i].size()));
            int* A = INTEGER(RA);
            for (size_t j = 0; j < pruned_adj_vect[i].size(); j++) {
                A[j] = pruned_adj_vect[i][j] + 1; // Convert to 1-based
            }
            SET_VECTOR_ELT(pruned_adj_list, i, RA);
            UNPROTECT(1); // RA
        }
        SET_VECTOR_ELT(res, 0, pruned_adj_list);
        UNPROTECT(1); // pruned_adj_list
    }
    
    // 1: pruned_edge_lengths_list
    {
        SEXP pruned_edge_lengths_list = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SEXP RD = PROTECT(Rf_allocVector(REALSXP, pruned_edge_lengths_vect[i].size()));
            double* D = REAL(RD);
            for (size_t j = 0; j < pruned_edge_lengths_vect[i].size(); j++) {
                D[j] = pruned_edge_lengths_vect[i][j];
            }
            SET_VECTOR_ELT(pruned_edge_lengths_list, i, RD);
            UNPROTECT(1); // RD
        }
        SET_VECTOR_ELT(res, 1, pruned_edge_lengths_list);
        UNPROTECT(1); // pruned_edge_lengths_list
    }
    
    // 2: path_lengths
    {
        int n_path_lengths = static_cast<int>(path_lengths.size());
        SEXP R_path_lengths = PROTECT(Rf_allocVector(REALSXP, n_path_lengths));
        double* path_lengths_array = REAL(R_path_lengths);
        for (int i = 0; i < n_path_lengths; i++) {
            path_lengths_array[i] = path_lengths[i];
        }
        SET_VECTOR_ELT(res, 2, R_path_lengths);
        UNPROTECT(1); // R_path_lengths
    }
    
    // 3: edge_lengths
    {
        int n_edge_lengths = static_cast<int>(edge_lengths.size());
        int n_path_lengths = static_cast<int>(path_lengths.size());
        if (n_path_lengths != n_edge_lengths) {
            UNPROTECT(1); // res
            Rf_error("n_edge_lengths is not equal to n_path_lengths.");
        }
        SEXP R_edge_lengths = PROTECT(Rf_allocVector(REALSXP, n_edge_lengths));
        double* edge_lengths_array = REAL(R_edge_lengths);
        for (int i = 0; i < n_edge_lengths; i++) {
            edge_lengths_array[i] = edge_lengths[i];
        }
        SET_VECTOR_ELT(res, 3, R_edge_lengths);
        UNPROTECT(1); // R_edge_lengths
    }
    
    // Set names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("edge_lengths_list"));
    SET_STRING_ELT(names, 2, Rf_mkChar("path_lengths"));
    SET_STRING_ELT(names, 3, Rf_mkChar("edge_lengths"));
    Rf_setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(2); // res, names
    return res;
}

} // extern "C"