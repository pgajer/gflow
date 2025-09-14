// graph_smoothing_functions_hardened.cpp — rchk-hardened drop-in for S_graph_kernel_smoother
//
// Hardenings applied:
// - Coercion block for y with PROTECT_WITH_INDEX/REPROTECT; copy to STL then UNPROTECT immediately.
// - Defensive scalar extraction via Rf_as* + NA/range checks (bandwidth>0, tolerance>0, max_iterations>0, kernel_type>=0, logicals not NA).
// - Container-first conversions before any result allocations.
// - Result assembly keeps only `result` + `names` protected to tail → final UNPROTECT(2).
// - Long-vector safety: XLENGTH for input sizes; explicit casts when copying.
//
// Assumptions: helper converters and the core algorithm exist and match the declarations below.

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
    // If your original returned extra fields (e.g., bw_predictions), they can be added here
    // and attached below in the result assembly in the same local-protect style.
};

// Core computation function (assumed available)
kernel_smooth_result_t graph_kernel_smoother(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    const std::vector<double>& weights,   // may be empty to indicate implicit weights
    double bandwidth,
    int kernel_type,
    bool adaptive_bandwidth,
    int max_iterations,
    double tolerance,
    bool verbose);

extern "C" {

SEXP S_graph_kernel_smoother(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP y_sexp,
    SEXP weights_sexp,              // may be NULL/Nil to indicate equal weights
    SEXP bandwidth_sexp,
    SEXP kernel_type_sexp,
    SEXP adaptive_bandwidth_sexp,
    SEXP max_iterations_sexp,
    SEXP tolerance_sexp,
    SEXP verbose_sexp
) {
    // --- y coercion block ---
    std::vector<double> y;
    R_xlen_t ny = 0;
    {
        SEXP y_real = y_sexp;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(y_real, &py);
        if (TYPEOF(y_real) != REALSXP) {
            REPROTECT(y_real = Rf_coerceVector(y_real, REALSXP), py);
        }
        ny = XLENGTH(y_real);
        if (ny <= 0) { UNPROTECT(1); Rf_error("S_graph_kernel_smoother(): 'y' must be a non-empty numeric vector."); }
        y.assign(REAL(y_real), REAL(y_real) + (size_t)ny);
        UNPROTECT(1); // y_real
    }

    // --- optional weights coercion block ---
    std::vector<double> weights;
    if (weights_sexp != R_NilValue) {
        SEXP w_real = weights_sexp;
        PROTECT_INDEX pw;
        PROTECT_WITH_INDEX(w_real, &pw);
        if (TYPEOF(w_real) != REALSXP) {
            REPROTECT(w_real = Rf_coerceVector(w_real, REALSXP), pw);
        }
        const R_xlen_t nw = XLENGTH(w_real);
        if (nw != ny) { UNPROTECT(1); Rf_error("S_graph_kernel_smoother(): 'weights' must have the same length as 'y'."); }
        weights.assign(REAL(w_real), REAL(w_real) + (size_t)nw);
        UNPROTECT(1); // w_real
    }

    // --- Container-first graph inputs ---
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    if ((size_t)ny != adj_list.size()) {
        Rf_error("S_graph_kernel_smoother(): length(y) must match number of graph vertices (%llu).",
                 (unsigned long long)adj_list.size());
    }

    // --- Defensive scalars ---
    const double bandwidth        = Rf_asReal(bandwidth_sexp);
    const int    kernel_type_i    = Rf_asInteger(kernel_type_sexp);
    const int    adaptive_bw_i    = Rf_asLogical(adaptive_bandwidth_sexp);
    const int    max_iterations_i = Rf_asInteger(max_iterations_sexp);
    const double tolerance        = Rf_asReal(tolerance_sexp);
    const int    verbose_i        = Rf_asLogical(verbose_sexp);

    // --- NA / range checks ---
    if (ISNAN(bandwidth) || bandwidth <= 0.0) {
        Rf_error("S_graph_kernel_smoother(): 'bandwidth' must be > 0.");
    }
    if (kernel_type_i == NA_INTEGER || kernel_type_i < 0) {
        Rf_error("S_graph_kernel_smoother(): 'kernel_type' must be a non-negative integer.");
    }
    if (adaptive_bw_i == NA_LOGICAL) {
        Rf_error("S_graph_kernel_smoother(): 'adaptive_bandwidth' must be TRUE/FALSE.");
    }
    if (max_iterations_i == NA_INTEGER || max_iterations_i <= 0) {
        Rf_error("S_graph_kernel_smoother(): 'max_iterations' must be a positive integer.");
    }
    if (ISNAN(tolerance) || tolerance <= 0.0) {
        Rf_error("S_graph_kernel_smoother(): 'tolerance' must be > 0.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_graph_kernel_smoother(): 'verbose' must be TRUE/FALSE.");
    }

    const int  kernel_type        = kernel_type_i;
    const bool adaptive_bandwidth = (adaptive_bw_i == TRUE);
    const bool verbose            = (verbose_i == TRUE);

    // --- Core computation ---
    kernel_smooth_result_t res = graph_kernel_smoother(
        adj_list, weight_list, y, weights,
        bandwidth, kernel_type, adaptive_bandwidth, max_iterations_i, tolerance, verbose
    );

    // --- Result assembly (list + names) ---
    const int N = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("smoothed_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("kernel_weights"));
    SET_STRING_ELT(names, 2, Rf_mkChar("bandwidth"));
    SET_STRING_ELT(names, 3, Rf_mkChar("effective_df"));
    SET_STRING_ELT(names, 4, Rf_mkChar("influence_diagonal"));
    SET_STRING_ELT(names, 5, Rf_mkChar("kernel_type"));
    SET_STRING_ELT(names, 6, Rf_mkChar("converged"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0: smoothed_values
    {
        SEXP s = PROTECT(convert_vector_double_to_R(res.smoothed_values));
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1);
    }
    // 1: kernel_weights
    {
        SEXP s = PROTECT(convert_vector_double_to_R(res.kernel_weights));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }
    // 2: bandwidth
    {
        SEXP s = PROTECT(Rf_ScalarReal(res.bandwidth));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }
    // 3: effective_df
    {
        SEXP s = PROTECT(Rf_ScalarReal(res.effective_df));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }
    // 4: influence_diagonal
    {
        SEXP s = PROTECT(convert_vector_double_to_R(res.influence_diagonal));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }
    // 5: kernel_type
    {
        SEXP s = PROTECT(Rf_ScalarInteger(res.kernel_type));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }
    // 6: converged
    {
        SEXP s = PROTECT(Rf_ScalarLogical(res.converged ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
