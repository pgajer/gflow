/**
 * @brief Fixed version of nada_graph_spectral_lowess_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_nada_graph_spectral_lowess
 * 
 * Issues fixed:
 * 1. S_nada_graph_spectral_lowess (line 214): UNPROTECT(variable) - unsupported
 * 2. Lambda at line 177: possible stack imbalance
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constants
 * 2. Removed lambda function that captured protect_count
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
extern SEXP convert_matrix_double_to_R(const std::vector<std::vector<double>>&);

// Core structures and declarations
struct nada_spectral_result_t {
    std::vector<double> smoothed_values;
    std::vector<double> eigenvector_coefficients;
    std::vector<std::vector<double>> eigenvectors_used;
    std::vector<double> eigenvalues_used;
    double bandwidth;
    double cv_score;
    int n_eigenvectors;
    bool converged;
};

// Core computation function (assumed available)
nada_spectral_result_t nada_graph_spectral_lowess(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    double bandwidth,
    int max_eigenvectors,
    double eigenvalue_threshold,
    bool auto_select_bandwidth,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_nada_graph_spectral_lowess
 * Fixes: UNPROTECT(variable) at line 214, lambda at line 177
 * Solution: Remove lambda, use container-first pattern, literal UNPROTECT
 */
SEXP S_nada_graph_spectral_lowess(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP y_sexp,
    SEXP weights_sexp,
    SEXP bandwidth_sexp,
    SEXP max_eigenvectors_sexp,
    SEXP eigenvalue_threshold_sexp,
    SEXP auto_select_bandwidth_sexp,
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
    
    // Convert inputs to C++ vectors
    R_xlen_t ny = XLENGTH(y_real);
    R_xlen_t nw = XLENGTH(weights_real);
    
    if (ny != nw) {
        UNPROTECT(2); // y_real, weights_real
        Rf_error("y and weights must have the same length");
    }
    
    const double* y_ptr = REAL(y_real);
    const double* weights_ptr = REAL(weights_real);
    
    std::vector<double> y;
    std::vector<double> weights;
    y.assign(y_ptr, y_ptr + static_cast<size_t>(ny));
    weights.assign(weights_ptr, weights_ptr + static_cast<size_t>(nw));
    
    // Convert other inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    // Extract scalar parameters
    double bandwidth = Rf_asReal(bandwidth_sexp);
    int max_eigenvectors = Rf_asInteger(max_eigenvectors_sexp);
    double eigenvalue_threshold = Rf_asReal(eigenvalue_threshold_sexp);
    bool auto_select_bandwidth = (Rf_asLogical(auto_select_bandwidth_sexp) == TRUE);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Core computation
    nada_spectral_result_t result_data = nada_graph_spectral_lowess(
        adj_list,
        weight_list,
        y,
        weights,
        bandwidth,
        max_eigenvectors,
        eigenvalue_threshold,
        auto_select_bandwidth,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 8;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: smoothed_values
    {
        SEXP smoothed = PROTECT(convert_vector_double_to_R(result_data.smoothed_values));
        SET_VECTOR_ELT(result, 0, smoothed);
        UNPROTECT(1);
    }
    
    // 1: eigenvector_coefficients
    {
        SEXP coeffs = PROTECT(convert_vector_double_to_R(result_data.eigenvector_coefficients));
        SET_VECTOR_ELT(result, 1, coeffs);
        UNPROTECT(1);
    }
    
    // 2: eigenvectors_used (matrix)
    {
        SEXP eigenvecs = PROTECT(convert_matrix_double_to_R(result_data.eigenvectors_used));
        SET_VECTOR_ELT(result, 2, eigenvecs);
        UNPROTECT(1);
    }
    
    // 3: eigenvalues_used
    {
        SEXP eigenvals = PROTECT(convert_vector_double_to_R(result_data.eigenvalues_used));
        SET_VECTOR_ELT(result, 3, eigenvals);
        UNPROTECT(1);
    }
    
    // 4: bandwidth
    {
        SEXP bw = PROTECT(Rf_ScalarReal(result_data.bandwidth));
        SET_VECTOR_ELT(result, 4, bw);
        UNPROTECT(1);
    }
    
    // 5: cv_score
    {
        SEXP cv = PROTECT(Rf_ScalarReal(result_data.cv_score));
        SET_VECTOR_ELT(result, 5, cv);
        UNPROTECT(1);
    }
    
    // 6: n_eigenvectors
    {
        SEXP n_eig = PROTECT(Rf_ScalarInteger(result_data.n_eigenvectors));
        SET_VECTOR_ELT(result, 6, n_eig);
        UNPROTECT(1);
    }
    
    // 7: converged
    {
        SEXP conv = PROTECT(Rf_ScalarLogical(result_data.converged ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 7, conv);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("smoothed_values"));
        SET_STRING_ELT(names, 1, Rf_mkChar("eigenvector_coefficients"));
        SET_STRING_ELT(names, 2, Rf_mkChar("eigenvectors_used"));
        SET_STRING_ELT(names, 3, Rf_mkChar("eigenvalues_used"));
        SET_STRING_ELT(names, 4, Rf_mkChar("bandwidth"));
        SET_STRING_ELT(names, 5, Rf_mkChar("cv_score"));
        SET_STRING_ELT(names, 6, Rf_mkChar("n_eigenvectors"));
        SET_STRING_ELT(names, 7, Rf_mkChar("converged"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(3); // y_real, weights_real, result
    return result;
}

} // extern "C"