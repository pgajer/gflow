/**
 * @brief Fixed version of parameterize_circular_graph_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_parameterize_circular_graph
 * 
 * Issues fixed:
 * 1. S_parameterize_circular_graph (line 89): UNPROTECT(variable) - unsupported
 * 2. Lambda at line 73: possible stack imbalance
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

// Core structures and declarations
struct circular_parameterization_t {
    std::vector<double> angles;
    std::vector<double> arc_lengths;
    double total_length;
    std::vector<int> ordering;
    double residual_error;
    int n_iterations;
    bool converged;
};

// Core computation function (assumed available)
circular_parameterization_t parameterize_circular_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start_vertex,
    double tolerance,
    int max_iterations,
    bool normalize_angles,
    bool verbose);

extern "C" {

/**
 * Fixed version of S_parameterize_circular_graph
 * Fixes: UNPROTECT(variable) at line 89, lambda at line 73
 * Solution: Remove lambda, use container-first pattern, literal UNPROTECT
 */
SEXP S_parameterize_circular_graph(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP start_vertex_sexp,
    SEXP tolerance_sexp,
    SEXP max_iterations_sexp,
    SEXP normalize_angles_sexp,
    SEXP verbose_sexp
) {
    // Convert inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    
    // Extract scalar parameters
    int start_vertex = Rf_asInteger(start_vertex_sexp) - 1; // Convert to 0-based
    double tolerance = Rf_asReal(tolerance_sexp);
    int max_iterations = Rf_asInteger(max_iterations_sexp);
    bool normalize_angles = (Rf_asLogical(normalize_angles_sexp) == TRUE);
    bool verbose = (Rf_asLogical(verbose_sexp) == TRUE);
    
    // Validate start_vertex
    int n_vertices = static_cast<int>(adj_list.size());
    if (start_vertex < 0 || start_vertex >= n_vertices) {
        Rf_error("start_vertex out of bounds");
    }
    
    // Core computation
    circular_parameterization_t param = parameterize_circular_graph(
        adj_list,
        weight_list,
        start_vertex,
        tolerance,
        max_iterations,
        normalize_angles,
        verbose
    );
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: angles
    {
        SEXP angles = PROTECT(convert_vector_double_to_R(param.angles));
        SET_VECTOR_ELT(result, 0, angles);
        UNPROTECT(1);
    }
    
    // 1: arc_lengths
    {
        SEXP arc_lengths = PROTECT(convert_vector_double_to_R(param.arc_lengths));
        SET_VECTOR_ELT(result, 1, arc_lengths);
        UNPROTECT(1);
    }
    
    // 2: total_length
    {
        SEXP total_length = PROTECT(Rf_ScalarReal(param.total_length));
        SET_VECTOR_ELT(result, 2, total_length);
        UNPROTECT(1);
    }
    
    // 3: ordering (convert to 1-based for R)
    {
        SEXP ordering = PROTECT(Rf_allocVector(INTSXP, param.ordering.size()));
        int* ord_ptr = INTEGER(ordering);
        for (size_t i = 0; i < param.ordering.size(); i++) {
            ord_ptr[i] = param.ordering[i] + 1; // Convert to 1-based
        }
        SET_VECTOR_ELT(result, 3, ordering);
        UNPROTECT(1);
    }
    
    // 4: residual_error
    {
        SEXP residual = PROTECT(Rf_ScalarReal(param.residual_error));
        SET_VECTOR_ELT(result, 4, residual);
        UNPROTECT(1);
    }
    
    // 5: n_iterations
    {
        SEXP iterations = PROTECT(Rf_ScalarInteger(param.n_iterations));
        SET_VECTOR_ELT(result, 5, iterations);
        UNPROTECT(1);
    }
    
    // 6: converged
    {
        SEXP converged = PROTECT(Rf_ScalarLogical(param.converged ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 6, converged);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("angles"));
        SET_STRING_ELT(names, 1, Rf_mkChar("arc_lengths"));
        SET_STRING_ELT(names, 2, Rf_mkChar("total_length"));
        SET_STRING_ELT(names, 3, Rf_mkChar("ordering"));
        SET_STRING_ELT(names, 4, Rf_mkChar("residual_error"));
        SET_STRING_ELT(names, 5, Rf_mkChar("n_iterations"));
        SET_STRING_ELT(names, 6, Rf_mkChar("converged"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    UNPROTECT(1); // result
    return result;
}

} // extern "C"