/**
 * @brief Fixed version of graph_kernel_smoother_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_graph_kernel_smoother
 * 
 * Issues fixed:
 * 1. S_graph_kernel_smoother (line 285): multiple pointer protection counters, UNPROTECT(variable)
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constants
 * 2. Used PROTECT_WITH_INDEX for conditional coercion
 * 3. Used container-first pattern consistently
 * 4. Fixed all PROTECT/UNPROTECT imbalances
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);

// Core structures and declarations
struct kernel_smooth_result_t {
    std::vector<double> smoothed_values;
    std::vector<double> kernel_weights;
    double bandwidth;
    double effective_df;
    std::vector<double> influence_diagonal;
    int kernel_type;
    bool converged;
};

// Core computation function (assumed available)
kernel_smooth_result_t graph_kernel_smoother(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    double bandwidth,
    int kernel_type,
    bool adaptive_bandwidth,
    int max_iterations,
    double tolerance,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_graph_kernel_smoother
 * Fixes: multiple pointer protection counters, UNPROTECT(variable) at line 285
 * Solution: Use PROTECT_WITH_INDEX, container-first pattern, literal UNPROTECT
 */
SEXP S_graph_kernel_smoother(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP y_sexp,
    SEXP weights_sexp,
    SEXP bandwidth_sexp,
    SEXP kernel_type_sexp,
    SEXP adaptive_bandwidth_sexp,
    SEXP max_iterations_sexp,
    SEXP tolerance_sexp,
    SEXP verbose_sexp
) {
    // Input coercion block with PROTECT_WITH_INDEX
    PROTECT_INDEX ipy, ipw;
    SEXP y_real = y_sexp;
    SEXP weights_real = weights_sexp;
    
    PROTECT_WITH_INDEX(y_real, &ipy);
    PROTECT_WITH_INDEX(weights_real, &ipw);
    
    if (TYPEOF(y_real) != REALSXP) {
        REPROTECT(y_real = Rf_coerceVector(y_real, REALSXP), ipy);
    }
    
    if (TYPEOF(weights_real) != REALSXP) {
        REPROTECT(weights_real = Rf_coerceVector(weights_real, REALSXP), ipw);
    }
    
    // Get lengths and validate
    R_xlen_t ny = XLENGTH(y_real);
    R_xlen_t nw = XLENGTH(weights_real);
    
    if (ny != nw) {
        UNPROTECT(2); // y_real, weights_real
        Rf_error("y and weights must have the same length");
    }
    
    // Convert inputs to C++ vectors
    const double* y_ptr = REAL(y_real);
    const double* weights_ptr = REAL(weights_real);
    
    std::vector<double> y;
    std::vector<double> weights;
    y.assign(y_ptr, y_ptr + static_cast<size_t>(ny));
    weights.assign(weights_ptr, weights_ptr + static_cast<size_t>(nw));
    
    // Convert other inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    // Validate graph size matches data
    if (static_cast<size_t>(ny) != adj_list.size()) {
        UNPROTECT(2); // y_real, weights_real
        Rf_error("y length must match graph size");
    }
    
    // Extract scalar parameters
    double bandwidth = Rf_asReal(bandwidth_sexp);
    int kernel_type = Rf_asInteger(kernel_type_sexp);
    bool adaptive_bandwidth = (Rf_asLogical(adaptive_bandwidth_sexp) == TRUE);
    int max_iterations = Rf_asInteger(max_iterations_sexp);
    double tolerance = Rf_asReal(tolerance_sexp);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Core computation
    kernel_smooth_result_t result_data = graph_kernel_smoother(
        adj_list,
        weight_list,
        y,
        weights,
        bandwidth,
        kernel_type,
        adaptive_bandwidth,
        max_iterations,
        tolerance,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: smoothed_values
    {
        SEXP smoothed = PROTECT(convert_vector_double_to_R(result_data.smoothed_values));
        SET_VECTOR_ELT(result, 0, smoothed);
        UNPROTECT(1);
    }
    
    // 1: kernel_weights
    {
        SEXP kernel_wts = PROTECT(convert_vector_double_to_R(result_data.kernel_weights));
        SET_VECTOR_ELT(result, 1, kernel_wts);
        UNPROTECT(1);
    }
    
    // 2: bandwidth
    {
        SEXP bw = PROTECT(Rf_ScalarReal(result_data.bandwidth));
        SET_VECTOR_ELT(result, 2, bw);
        UNPROTECT(1);
    }
    
    // 3: effective_df
    {
        SEXP edf = PROTECT(Rf_ScalarReal(result_data.effective_df));
        SET_VECTOR_ELT(result, 3, edf);
        UNPROTECT(1);
    }
    
    // 4: influence_diagonal
    {
        SEXP influence = PROTECT(convert_vector_double_to_R(result_data.influence_diagonal));
        SET_VECTOR_ELT(result, 4, influence);
        UNPROTECT(1);
    }
    
    // 5: kernel_type
    {
        SEXP ktype = PROTECT(Rf_ScalarInteger(result_data.kernel_type));
        SET_VECTOR_ELT(result, 5, ktype);
        UNPROTECT(1);
    }
    
    // 6: converged
    {
        SEXP conv = PROTECT(Rf_ScalarLogical(result_data.converged ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 6, conv);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("smoothed_values"));
        SET_STRING_ELT(names, 1, Rf_mkChar("kernel_weights"));
        SET_STRING_ELT(names, 2, Rf_mkChar("bandwidth"));
        SET_STRING_ELT(names, 3, Rf_mkChar("effective_df"));
        SET_STRING_ELT(names, 4, Rf_mkChar("influence_diagonal"));
        SET_STRING_ELT(names, 5, Rf_mkChar("kernel_type"));
        SET_STRING_ELT(names, 6, Rf_mkChar("converged"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(3); // y_real, weights_real, result
    return result;
}

} // extern "C"