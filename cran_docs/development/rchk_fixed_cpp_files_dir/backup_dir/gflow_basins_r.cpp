/**
 * @brief Fixed version of gflow_basins_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_create_basin_cx
 * 
 * Issues fixed:
 * 1. S_create_basin_cx (line 692): UNPROTECT(variable) - unsupported
 * 2. Lambda at line 384: possible stack imbalance
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constants
 * 2. Removed lambda function that captured protect_count
 * 3. Used container-first pattern consistently
 * 4. Fixed all PROTECT/UNPROTECT imbalances
 */

#include <vector>
#include <unordered_map>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);

// Core structures and declarations
struct basin_complex_t {
    std::vector<int> basin_labels;
    std::vector<double> basin_scores;
    std::vector<int> basin_sizes;
    std::vector<std::vector<int>> basin_members;
    std::vector<std::vector<double>> basin_boundaries;
    int n_basins;
    double flow_threshold;
};

// Core computation function (assumed available)
basin_complex_t create_basin_cx(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& vertex_values,
    double flow_threshold,
    int min_basin_size,
    bool merge_small_basins,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_create_basin_cx
 * Fixes: UNPROTECT(variable) at line 692, lambda at line 384
 * Solution: Remove lambda, use container-first pattern, literal UNPROTECT
 */
SEXP S_create_basin_cx(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP vertex_values_sexp,
    SEXP flow_threshold_sexp,
    SEXP min_basin_size_sexp,
    SEXP merge_small_basins_sexp,
    SEXP verbose_sexp
) {
    // Input coercion block with PROTECT_WITH_INDEX
    PROTECT_INDEX ipx;
    SEXP values_real = vertex_values_sexp;
    PROTECT_WITH_INDEX(values_real, &ipx);
    
    if (TYPEOF(values_real) != REALSXP) {
        REPROTECT(values_real = Rf_coerceVector(values_real, REALSXP), ipx);
    }
    
    // Convert inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    R_xlen_t n_values = XLENGTH(values_real);
    const double* values_ptr = REAL(values_real);
    std::vector<double> vertex_values;
    vertex_values.assign(values_ptr, values_ptr + static_cast<size_t>(n_values));
    
    // Extract scalar parameters
    double flow_threshold = Rf_asReal(flow_threshold_sexp);
    int min_basin_size = Rf_asInteger(min_basin_size_sexp);
    bool merge_small_basins = (Rf_asLogical(merge_small_basins_sexp) == TRUE);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Core computation
    basin_complex_t basin_cx = create_basin_cx(
        adj_list,
        weight_list,
        vertex_values,
        flow_threshold,
        min_basin_size,
        merge_small_basins,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: basin_labels
    {
        SEXP labels = PROTECT(convert_vector_int_to_R(basin_cx.basin_labels));
        SET_VECTOR_ELT(result, 0, labels);
        UNPROTECT(1);
    }
    
    // 1: basin_scores
    {
        SEXP scores = PROTECT(convert_vector_double_to_R(basin_cx.basin_scores));
        SET_VECTOR_ELT(result, 1, scores);
        UNPROTECT(1);
    }
    
    // 2: basin_sizes
    {
        SEXP sizes = PROTECT(convert_vector_int_to_R(basin_cx.basin_sizes));
        SET_VECTOR_ELT(result, 2, sizes);
        UNPROTECT(1);
    }
    
    // 3: basin_members (list of integer vectors)
    {
        SEXP members_list = PROTECT(Rf_allocVector(VECSXP, basin_cx.n_basins));
        for (int i = 0; i < basin_cx.n_basins; i++) {
            SEXP members = PROTECT(convert_vector_int_to_R(basin_cx.basin_members[i]));
            SET_VECTOR_ELT(members_list, i, members);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(result, 3, members_list);
        UNPROTECT(1);
    }
    
    // 4: basin_boundaries (list of numeric vectors)
    {
        SEXP boundaries_list = PROTECT(Rf_allocVector(VECSXP, basin_cx.n_basins));
        for (int i = 0; i < basin_cx.n_basins; i++) {
            SEXP boundary = PROTECT(convert_vector_double_to_R(basin_cx.basin_boundaries[i]));
            SET_VECTOR_ELT(boundaries_list, i, boundary);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(result, 4, boundaries_list);
        UNPROTECT(1);
    }
    
    // 5: n_basins
    {
        SEXP n_basins = PROTECT(Rf_ScalarInteger(basin_cx.n_basins));
        SET_VECTOR_ELT(result, 5, n_basins);
        UNPROTECT(1);
    }
    
    // 6: flow_threshold
    {
        SEXP threshold = PROTECT(Rf_ScalarReal(basin_cx.flow_threshold));
        SET_VECTOR_ELT(result, 6, threshold);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("basin_labels"));
        SET_STRING_ELT(names, 1, Rf_mkChar("basin_scores"));
        SET_STRING_ELT(names, 2, Rf_mkChar("basin_sizes"));
        SET_STRING_ELT(names, 3, Rf_mkChar("basin_members"));
        SET_STRING_ELT(names, 4, Rf_mkChar("basin_boundaries"));
        SET_STRING_ELT(names, 5, Rf_mkChar("n_basins"));
        SET_STRING_ELT(names, 6, Rf_mkChar("flow_threshold"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(2); // values_real, result
    return result;
}

} // extern "C"