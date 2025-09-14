// pruning_long_edges_find_shortest_alt_hardened.cpp — rchk-hardened drop-in
//
// Hardenings vs provided version:
// - Defensive scalar extraction via Rf_as* with NA/range checks (source/target not NA, direct_edge_length>=0, max_path_length>0, verbose not NA).
// - Use R_xlen_t for vertex counts (long-vector safety).
// - Keep only `result` + `names` protected until function tail → final UNPROTECT(2).
// - Container-first preserved; path indices returned 1-based for R.
// - Local PROTECT/UNPROTECT(1) around each temporary allocation.
//
// Assumptions: helper converters and core algorithm exist with the same signatures.

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Helpers (provided elsewhere)
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);

// Core structures
struct alt_path_result_t {
    std::vector<int> path;
    double path_length;
    int    n_edges;
    bool   path_found;
    double length_ratio;  // ratio to direct edge length
};

// Core algorithm (provided elsewhere)
alt_path_result_t find_shortest_alt_path(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int source,
    int target,
    double direct_edge_length,
    int max_path_length,
    bool verbose);

extern "C" {

SEXP S_find_shortest_alt_path(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP source_sexp,
    SEXP target_sexp,
    SEXP direct_edge_length_sexp,
    SEXP max_path_length_sexp,
    SEXP verbose_sexp
) {
    // --- Container-first conversions ---
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    const R_xlen_t n_vertices_x = (R_xlen_t)adj_list.size();

    // --- Defensive scalars ---
    const int    source_i  = Rf_asInteger(source_sexp);
    const int    target_i  = Rf_asInteger(target_sexp);
    const double direct_len= Rf_asReal(direct_edge_length_sexp);
    const int    maxlen_i  = Rf_asInteger(max_path_length_sexp);
    const int    verbose_i = Rf_asLogical(verbose_sexp);

    // --- NA / range checks ---
    if (source_i == NA_INTEGER || target_i == NA_INTEGER) {
        Rf_error("S_find_shortest_alt_path(): 'source' and 'target' cannot be NA.");
    }
    if (ISNAN(direct_len) || direct_len < 0.0) {
        Rf_error("S_find_shortest_alt_path(): 'direct_edge_length' must be >= 0.");
    }
    if (maxlen_i == NA_INTEGER || maxlen_i <= 0) {
        Rf_error("S_find_shortest_alt_path(): 'max_path_length' must be a positive integer.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_find_shortest_alt_path(): 'verbose' must be TRUE/FALSE.");
    }

    // Convert to 0-based
    const int source = source_i - 1;
    const int target = target_i - 1;

    if (source < 0 || (R_xlen_t)source >= n_vertices_x) {
        Rf_error("S_find_shortest_alt_path(): 'source' out of bounds (1..%lld).", (long long)n_vertices_x);
    }
    if (target < 0 || (R_xlen_t)target >= n_vertices_x) {
        Rf_error("S_find_shortest_alt_path(): 'target' out of bounds (1..%lld).", (long long)n_vertices_x);
    }

    const bool verbose = (verbose_i == TRUE);

    // --- Core computation ---
    alt_path_result_t result_data = find_shortest_alt_path(
        adj_list,
        weight_list,
        source,
        target,
        direct_len,
        maxlen_i,
        verbose
    );

    // --- Result assembly (list + names) ---
    const int N = 5;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("path"));
    SET_STRING_ELT(names, 1, Rf_mkChar("path_length"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_edges"));
    SET_STRING_ELT(names, 3, Rf_mkChar("path_found"));
    SET_STRING_ELT(names, 4, Rf_mkChar("length_ratio"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0: path (1-based indices for R)
    {
        std::vector<int> r_path = result_data.path;
        for (size_t i = 0; i < r_path.size(); ++i) r_path[i] += 1;
        SEXP v = PROTECT(convert_vector_int_to_R(r_path));
        SET_VECTOR_ELT(result, 0, v);
        UNPROTECT(1);
    }
    // 1: path_length
    {
        SEXP s = PROTECT(Rf_ScalarReal(result_data.path_length));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }
    // 2: n_edges
    {
        SEXP s = PROTECT(Rf_ScalarInteger(result_data.n_edges));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }
    // 3: path_found
    {
        SEXP s = PROTECT(Rf_ScalarLogical(result_data.path_found ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }
    // 4: length_ratio
    {
        SEXP s = PROTECT(Rf_ScalarReal(result_data.length_ratio));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
