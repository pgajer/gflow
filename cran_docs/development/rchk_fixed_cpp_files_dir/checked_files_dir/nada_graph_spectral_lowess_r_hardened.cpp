/**
 * @brief Hardened rchk-safe drop-in of S_nada_graph_spectral_lowess
 *
 * Changes vs. provided version:
 * - Keep `result` and `names` protected until the tail; finish with UNPROTECT(4) (y, weights, result, names).
 * - Maintain container-first local PROTECT/UNPROTECT(1) for each element.
 * - Defensive scalars via Rf_as*; length check with exact UNPROTECT before error.
 * - Use XLENGTH/R_xlen_t for long-vector safety on inputs.
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
    // --- Input coercion blocks (indexed protects) ---
    PROTECT_INDEX ipy, ipw;
    SEXP y_real = y_sexp;
    SEXP weights_real = weights_sexp;
    PROTECT_WITH_INDEX(y_real, &ipy);
    PROTECT_WITH_INDEX(weights_real, &ipw);
    if (TYPEOF(y_real) != REALSXP)       REPROTECT(y_real = Rf_coerceVector(y_real, REALSXP), ipy);
    if (TYPEOF(weights_real) != REALSXP) REPROTECT(weights_real = Rf_coerceVector(weights_real, REALSXP), ipw);

    const R_xlen_t ny = XLENGTH(y_real);
    const R_xlen_t nw = XLENGTH(weights_real);
    if (ny != nw) {
        UNPROTECT(2); // y_real, weights_real
        Rf_error("S_nada_graph_spectral_lowess: 'y' and 'weights' must have the same length");
    }

    // Copy to C++ containers while protected
    std::vector<double> y, weights;
    y.assign(REAL(y_real), REAL(y_real) + (size_t)ny);
    weights.assign(REAL(weights_real), REAL(weights_real) + (size_t)nw);

    // Other inputs via helpers
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    // Defensive scalars
    const double bandwidth            = Rf_asReal(bandwidth_sexp);
    const int    max_eigenvectors     = Rf_asInteger(max_eigenvectors_sexp);
    const double eigenvalue_threshold = Rf_asReal(eigenvalue_threshold_sexp);
    const bool   auto_select_bandwidth = (Rf_asLogical(auto_select_bandwidth_sexp) == TRUE);
    const bool   verbose               = (Rf_asLogical(verbose_sexp) == TRUE);

    // Core computation (no R allocations here)
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

    // --- Assemble result (container-first). Keep result & names protected until tail. ---
    const int N = 8;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));

    // 0: smoothed_values
    {
        SEXP s0 = PROTECT(convert_vector_double_to_R(result_data.smoothed_values));
        SET_VECTOR_ELT(result, 0, s0);
        UNPROTECT(1);
    }
    // 1: eigenvector_coefficients
    {
        SEXP s1 = PROTECT(convert_vector_double_to_R(result_data.eigenvector_coefficients));
        SET_VECTOR_ELT(result, 1, s1);
        UNPROTECT(1);
    }
    // 2: eigenvectors_used (matrix)
    {
        SEXP s2 = PROTECT(convert_matrix_double_to_R(result_data.eigenvectors_used));
        SET_VECTOR_ELT(result, 2, s2);
        UNPROTECT(1);
    }
    // 3: eigenvalues_used
    {
        SEXP s3 = PROTECT(convert_vector_double_to_R(result_data.eigenvalues_used));
        SET_VECTOR_ELT(result, 3, s3);
        UNPROTECT(1);
    }
    // 4: bandwidth
    {
        SEXP s4 = PROTECT(Rf_ScalarReal(result_data.bandwidth));
        SET_VECTOR_ELT(result, 4, s4);
        UNPROTECT(1);
    }
    // 5: cv_score
    {
        SEXP s5 = PROTECT(Rf_ScalarReal(result_data.cv_score));
        SET_VECTOR_ELT(result, 5, s5);
        UNPROTECT(1);
    }
    // 6: n_eigenvectors
    {
        SEXP s6 = PROTECT(Rf_ScalarInteger(result_data.n_eigenvectors));
        SET_VECTOR_ELT(result, 6, s6);
        UNPROTECT(1);
    }
    // 7: converged
    {
        SEXP s7 = PROTECT(Rf_ScalarLogical(result_data.converged ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 7, s7);
        UNPROTECT(1);
    }

    // names — keep protected until tail
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("smoothed_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("eigenvector_coefficients"));
    SET_STRING_ELT(names, 2, Rf_mkChar("eigenvectors_used"));
    SET_STRING_ELT(names, 3, Rf_mkChar("eigenvalues_used"));
    SET_STRING_ELT(names, 4, Rf_mkChar("bandwidth"));
    SET_STRING_ELT(names, 5, Rf_mkChar("cv_score"));
    SET_STRING_ELT(names, 6, Rf_mkChar("n_eigenvectors"));
    SET_STRING_ELT(names, 7, Rf_mkChar("converged"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(4); // y_real, weights_real, result, names
    return result;
}

} // extern "C"
