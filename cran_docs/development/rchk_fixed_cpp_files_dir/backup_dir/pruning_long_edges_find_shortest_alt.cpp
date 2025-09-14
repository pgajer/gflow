/**
 * @brief Fixed version of pruning_long_edges.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_find_shortest_alt_path
 * 
 * Issues fixed:
 * 1. S_find_shortest_alt_path (lines 794, 797): negative depth, over-unprotect, stack imbalance
 * 
 * Changes made:
 * 1. Used container-first pattern consistently
 * 2. Fixed all PROTECT/UNPROTECT imbalances
 * 3. Used literal UNPROTECT values
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);

// Core structures and declarations
struct alt_path_result_t {
    std::vector<int> path;
    double path_length;
    int n_edges;
    bool path_found;
    double length_ratio;  // ratio to direct edge length
};

// Core computation function (assumed available)
alt_path_result_t find_shortest_alt_path(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int source,
    int target,
    double direct_edge_length,
    int max_path_length,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_find_shortest_alt_path
 * Fixes: negative depth, over-unprotect at lines 794, 797
 * Solution: Use container-first pattern, literal UNPROTECT
 */
SEXP S_find_shortest_alt_path(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP source_sexp,
    SEXP target_sexp,
    SEXP direct_edge_length_sexp,
    SEXP max_path_length_sexp,
    SEXP verbose_sexp
) {
    // Convert inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    // Extract scalar parameters
    int source = Rf_asInteger(source_sexp) - 1; // Convert to 0-based
    int target = Rf_asInteger(target_sexp) - 1; // Convert to 0-based
    double direct_edge_length = Rf_asReal(direct_edge_length_sexp);
    int max_path_length = Rf_asInteger(max_path_length_sexp);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Validate indices
    int n_vertices = static_cast<int>(adj_list.size());
    if (source < 0 || source >= n_vertices) {
        Rf_error("source vertex out of bounds");
    }
    if (target < 0 || target >= n_vertices) {
        Rf_error("target vertex out of bounds");
    }
    
    // Core computation
    alt_path_result_t result_data = find_shortest_alt_path(
        adj_list,
        weight_list,
        source,
        target,
        direct_edge_length,
        max_path_length,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 5;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: path (convert to 1-based for R)
    {
        std::vector<int> r_path = result_data.path;
        for (size_t i = 0; i < r_path.size(); i++) {
            r_path[i]++; // Convert to 1-based
        }
        SEXP path = PROTECT(convert_vector_int_to_R(r_path));
        SET_VECTOR_ELT(result, 0, path);
        UNPROTECT(1);
    }
    
    // 1: path_length
    {
        SEXP length = PROTECT(Rf_ScalarReal(result_data.path_length));
        SET_VECTOR_ELT(result, 1, length);
        UNPROTECT(1);
    }
    
    // 2: n_edges
    {
        SEXP n_edges = PROTECT(Rf_ScalarInteger(result_data.n_edges));
        SET_VECTOR_ELT(result, 2, n_edges);
        UNPROTECT(1);
    }
    
    // 3: path_found
    {
        SEXP found = PROTECT(Rf_ScalarLogical(result_data.path_found ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 3, found);
        UNPROTECT(1);
    }
    
    // 4: length_ratio
    {
        SEXP ratio = PROTECT(Rf_ScalarReal(result_data.length_ratio));
        SET_VECTOR_ELT(result, 4, ratio);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("path"));
        SET_STRING_ELT(names, 1, Rf_mkChar("path_length"));
        SET_STRING_ELT(names, 2, Rf_mkChar("n_edges"));
        SET_STRING_ELT(names, 3, Rf_mkChar("path_found"));
        SET_STRING_ELT(names, 4, Rf_mkChar("length_ratio"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(1); // result
    return result;
}

} // extern "C"