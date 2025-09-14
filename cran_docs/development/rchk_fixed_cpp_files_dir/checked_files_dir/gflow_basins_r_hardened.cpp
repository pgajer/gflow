// gflow_basins_r_hardened.cpp — rchk-hardened drop-in for S_create_basin_cx
//
// Hardenings vs provided file:
// - Coercion block for vertex_values with PROTECT_WITH_INDEX/REPROTECT; copy to STL then UNPROTECT immediately.
// - Defensive scalar extraction via Rf_as* + NA/range checks (flow_threshold>=0, min_basin_size>0, logicals not NA).
// - Result assembly keeps only `result` + `names` protected to tail → final UNPROTECT(2).
// - Long-vector safety: XLENGTH for input sizes; explicit casts when copying; loops use size_t/R_xlen_t appropriately.
// - Container-first preserved.
//
// Assumptions: helper converters and core algorithm exist with the same signatures as in the original.

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

SEXP S_create_basin_cx(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP vertex_values_sexp,
    SEXP flow_threshold_sexp,
    SEXP min_basin_size_sexp,
    SEXP merge_small_basins_sexp,
    SEXP verbose_sexp
) {
    // --- Coercion block for vertex_values ---
    std::vector<double> vertex_values;
    {
        SEXP values_real = vertex_values_sexp;
        PROTECT_INDEX pv;
        PROTECT_WITH_INDEX(values_real, &pv);
        if (TYPEOF(values_real) != REALSXP) {
            REPROTECT(values_real = Rf_coerceVector(values_real, REALSXP), pv);
        }
        const R_xlen_t n_values = XLENGTH(values_real);
        if (n_values <= 0) {
            UNPROTECT(1);
            Rf_error("S_create_basin_cx(): 'vertex_values' must be a non-empty numeric vector.");
        }
        vertex_values.assign(REAL(values_real), REAL(values_real) + (size_t)n_values);
        UNPROTECT(1); // values_real
    }

    // --- Container-first for graph inputs ---
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    if (vertex_values.size() != adj_list.size()) {
        Rf_error("S_create_basin_cx(): length(vertex_values) must match number of graph vertices (%llu).",
                 (unsigned long long)adj_list.size());
    }

    // --- Defensive scalars ---
    const double flow_threshold     = Rf_asReal(flow_threshold_sexp);
    const int    min_basin_size_i   = Rf_asInteger(min_basin_size_sexp);
    const int    merge_small_i      = Rf_asLogical(merge_small_basins_sexp);
    const int    verbose_i          = Rf_asLogical(verbose_sexp);

    // --- NA / range checks ---
    if (ISNAN(flow_threshold) || flow_threshold < 0.0) {
        Rf_error("S_create_basin_cx(): 'flow_threshold' must be >= 0.");
    }
    if (min_basin_size_i == NA_INTEGER || min_basin_size_i <= 0) {
        Rf_error("S_create_basin_cx(): 'min_basin_size' must be a positive integer.");
    }
    if (merge_small_i == NA_LOGICAL) {
        Rf_error("S_create_basin_cx(): 'merge_small_basins' must be TRUE/FALSE.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_create_basin_cx(): 'verbose' must be TRUE/FALSE.");
    }

    const int  min_basin_size    = min_basin_size_i;
    const bool merge_small_basins= (merge_small_i == TRUE);
    const bool verbose           = (verbose_i == TRUE);

    // --- Core computation ---
    basin_complex_t basin_cx = create_basin_cx(
        adj_list,
        weight_list,
        vertex_values,
        flow_threshold,
        min_basin_size,
        merge_small_basins,
        verbose
    );

    // --- Result assembly (list + names) ---
    const int N = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("basin_labels"));
    SET_STRING_ELT(names, 1, Rf_mkChar("basin_scores"));
    SET_STRING_ELT(names, 2, Rf_mkChar("basin_sizes"));
    SET_STRING_ELT(names, 3, Rf_mkChar("basin_members"));
    SET_STRING_ELT(names, 4, Rf_mkChar("basin_boundaries"));
    SET_STRING_ELT(names, 5, Rf_mkChar("n_basins"));
    SET_STRING_ELT(names, 6, Rf_mkChar("flow_threshold"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0: basin_labels
    {
        SEXP s = PROTECT(convert_vector_int_to_R(basin_cx.basin_labels));
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1);
    }
    // 1: basin_scores
    {
        SEXP s = PROTECT(convert_vector_double_to_R(basin_cx.basin_scores));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }
    // 2: basin_sizes
    {
        SEXP s = PROTECT(convert_vector_int_to_R(basin_cx.basin_sizes));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }
    // 3: basin_members (list of integer vectors)
    {
        const int nb = basin_cx.n_basins;
        SEXP lst = PROTECT(Rf_allocVector(VECSXP, nb));
        for (int i = 0; i < nb; ++i) {
            SEXP v = PROTECT(convert_vector_int_to_R(basin_cx.basin_members[(size_t)i]));
            SET_VECTOR_ELT(lst, i, v);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(result, 3, lst);
        UNPROTECT(1);
    }
    // 4: basin_boundaries (list of numeric vectors)
    {
        const int nb = basin_cx.n_basins;
        SEXP lst = PROTECT(Rf_allocVector(VECSXP, nb));
        for (int i = 0; i < nb; ++i) {
            SEXP v = PROTECT(convert_vector_double_to_R(basin_cx.basin_boundaries[(size_t)i]));
            SET_VECTOR_ELT(lst, i, v);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(result, 4, lst);
        UNPROTECT(1);
    }
    // 5: n_basins
    {
        SEXP s = PROTECT(Rf_ScalarInteger(basin_cx.n_basins));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }
    // 6: flow_threshold
    {
        SEXP s = PROTECT(Rf_ScalarReal(basin_cx.flow_threshold));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
