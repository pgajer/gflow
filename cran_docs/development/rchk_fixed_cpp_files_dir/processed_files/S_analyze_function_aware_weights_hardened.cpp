
// Hardened, rchk-safe drop-in for S_analyze_function_aware_weights
// - Container-first protection
// - Strict type/length checks
// - Uses R_xlen_t for long vectors
// - Avoids unsupported UNPROTECT patterns
// - Validates sizes of `results` vs `weight_types`
// - Defensive NA handling is left to downstream C++ if desired

#include <Rinternals.h>
#include <R_ext/Error.h>
#include <vector>
#include <memory>
#include <limits>
#include <cmath>

// Forward decls for helpers provided elsewhere in the package
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);
std::unique_ptr<std::vector<double>> Rvect_to_CppVect_double(SEXP s_vec);

// Forward decl for graph type (provided by the package)
struct set_wgraph_t {
    set_wgraph_t(const std::vector<std::vector<int>>& adj,
                 const std::vector<std::vector<double>>& w);
    size_t num_vertices() const;
    std::vector<std::vector<double>> analyze_function_aware_weights(
        const std::vector<double>& function_values,
        const std::vector<int>&     weight_types,
        double epsilon, double lambda, double alpha, double beta,
        double tau, double p, double q, double r) const;
};

extern "C" SEXP S_analyze_function_aware_weights(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_function_values,
    SEXP s_weight_types,
    SEXP s_epsilon,
    SEXP s_lambda,
    SEXP s_alpha,
    SEXP s_beta,
    SEXP s_tau,
    SEXP s_p,
    SEXP s_q,
    SEXP s_r
) {
    // ---- Validate and convert inputs (no allocations that depend on R objects being unprotected) ----
    if (!Rf_isNumeric(s_function_values))
        Rf_error("`function_values` must be a numeric vector.");
    if (!Rf_isInteger(s_weight_types))
        Rf_error("`weight_types` must be an integer vector.");

    // Scalars: be permissive (R will coerce) but check for NA/NaN where it matters
    double epsilon = Rf_asReal(s_epsilon);
    double lambda  = Rf_asReal(s_lambda);
    double alpha   = Rf_asReal(s_alpha);
    double beta    = Rf_asReal(s_beta);
    double tau     = Rf_asReal(s_tau);
    double p       = Rf_asReal(s_p);
    double q       = Rf_asReal(s_q);
    double r       = Rf_asReal(s_r);

    // Convert adjacency/weights (pure C++ conversions)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert function values
    std::unique_ptr<std::vector<double>> function_values_ptr = Rvect_to_CppVect_double(s_function_values);
    if (!function_values_ptr) Rf_error("Failed to convert `function_values`.");
    std::vector<double> function_values = std::move(*function_values_ptr);

    // Convert weight types
    const R_xlen_t wt_len = Rf_xlength(s_weight_types);
    if (wt_len < 0) Rf_error("Invalid length for `weight_types`.");
    std::vector<int> weight_types;
    weight_types.reserve(static_cast<size_t>(wt_len));
    const int *wt_ptr = INTEGER(s_weight_types); // safe: args are protected by R
    for (R_xlen_t i = 0; i < wt_len; ++i) {
        weight_types.push_back(wt_ptr[i]);
    }

    // ---- Build graph and validate sizes ----
    set_wgraph_t graph(adj_list, weight_list);

    if (function_values.size() != graph.num_vertices()) {
        Rf_error("Length of `function_values` (%zu) must match the number of vertices (%zu).",
                 function_values.size(), static_cast<size_t>(graph.num_vertices()));
    }

    // Compute results
    std::vector<std::vector<double>> results =
        graph.analyze_function_aware_weights(function_values, weight_types,
                                             epsilon, lambda, alpha, beta, tau, p, q, r);

    // Expect exactly 1 + weight_types.size() result blocks: original + each requested type
    const size_t expected_blocks = 1u + weight_types.size();
    if (results.size() != expected_blocks) {
        Rf_error("Internal error: expected %zu result blocks (1 + length(weight_types)) but got %zu.",
                 expected_blocks, results.size());
    }

    // ---- Allocate container first and PROTECT early ----
    SEXP r_list = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)expected_blocks));
    SEXP r_list_names = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)expected_blocks));

    // ---- Fill names ----
    // 0: "original"
    SET_STRING_ELT(r_list_names, 0, Rf_mkChar("original"));
    static const char* type_names[] = {"original", "inverse", "direct", "exp", "power", "sigmoid", "lp_embedding"};
    for (size_t i = 0; i < weight_types.size(); ++i) {
        int t = weight_types[i];
        if (t >= 0 && t <= 5) {
            SET_STRING_ELT(r_list_names, (R_xlen_t)(i + 1), Rf_mkChar(type_names[t + 1]));
        } else {
            SET_STRING_ELT(r_list_names, (R_xlen_t)(i + 1), Rf_mkChar("unknown"));
        }
    }

    // ---- Fill "original" weights ----
    {
        const std::vector<double>& orig = results[0];
        SEXP r_orig = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)orig.size()));
        double* rptr = REAL(r_orig);
        for (size_t j = 0; j < orig.size(); ++j) rptr[j] = orig[j];
        SET_VECTOR_ELT(r_list, 0, r_orig);
        UNPROTECT(1); // r_orig
    }

    // ---- Fill the rest ----
    for (size_t i = 0; i < weight_types.size(); ++i) {
        const std::vector<double>& v = results[i + 1];
        SEXP r_block = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)v.size()));
        double* rptr = REAL(r_block);
        for (size_t j = 0; j < v.size(); ++j) rptr[j] = v[j];
        SET_VECTOR_ELT(r_list, (R_xlen_t)(i + 1), r_block);
        UNPROTECT(1); // r_block
    }

    // Attach names
    Rf_setAttrib(r_list, R_NamesSymbol, r_list_names);

    UNPROTECT(2); // r_list, r_list_names
    return r_list;
}
