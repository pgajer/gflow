// density_hardened.cpp — rchk-hardened drop-in for S_estimate_local_density_over_grid
//
// Hardenings vs provided file:
// - Coercion block for x with PROTECT_WITH_INDEX/REPROTECT; copy to STL then UNPROTECT immediately.
// - Defensive scalar extraction via Rf_as* + NA/range checks (nx>0, grid_size>0, pilot_bandwidth>0, kernel_type>=0, verbose not NA, poffset not NA).
/* - Result assembly keeps only `result` + `names` protected to tail → final UNPROTECT(2).
 * - Long-vector safety: XLENGTH for input sizes; explicit casts when copying.
 * - Container-first preserved.
 */
//
// Assumptions: helper converters and core algorithm exist with the same signatures as in the original.

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declaration for helper function (assumed available)
extern SEXP convert_vector_double_to_R(const std::vector<double>&);

// Core function structure and declaration
struct gdensity_t {
    std::vector<double> density;
    double bandwidth;
    bool auto_selected;
    double offset;
    double start;
    double end;
};

// Core computation function (assumed available)
gdensity_t estimate_local_density_over_grid(
    const std::vector<double>& x,
    int grid_size,
    double poffset,
    double pilot_bandwidth,
    int kernel_type,
    bool verbose);

extern "C" {

SEXP S_estimate_local_density_over_grid(SEXP s_x,
                                        SEXP s_grid_size,
                                        SEXP s_poffset,
                                        SEXP s_pilot_bandwidth,
                                        SEXP s_kernel_type,
                                        SEXP s_verbose) {
    // --- Coerce x to REAL and copy into std::vector<double> ---
    std::vector<double> x;
    {
        SEXP sx = s_x;
        PROTECT_INDEX ipx;
        PROTECT_WITH_INDEX(sx, &ipx);
        if (TYPEOF(sx) != REALSXP) {
            REPROTECT(sx = Rf_coerceVector(sx, REALSXP), ipx);
        }
        const R_xlen_t nx = XLENGTH(sx);
        if (nx <= 0) {
            UNPROTECT(1);
            Rf_error("S_estimate_local_density_over_grid(): 'x' must be a non-empty numeric vector.");
        }
        x.assign(REAL(sx), REAL(sx) + (size_t)nx);
        UNPROTECT(1); // sx
    }

    // --- Scalars (defensive) ---
    const int    grid_size_i      = Rf_asInteger(s_grid_size);
    const double poffset          = Rf_asReal(s_poffset);
    const double pilot_bandwidth  = Rf_asReal(s_pilot_bandwidth);
    const int    kernel_type_i    = Rf_asInteger(s_kernel_type);
    const int    verbose_i        = Rf_asLogical(s_verbose);

    // --- NA / range checks ---
    if (grid_size_i == NA_INTEGER || grid_size_i <= 0) {
        Rf_error("S_estimate_local_density_over_grid(): 'grid_size' must be a positive integer.");
    }
    if (ISNAN(poffset)) {
        Rf_error("S_estimate_local_density_over_grid(): 'poffset' cannot be NA.");
    }
    if (ISNAN(pilot_bandwidth) || pilot_bandwidth <= 0.0) {
        Rf_error("S_estimate_local_density_over_grid(): 'pilot_bandwidth' must be > 0.");
    }
    if (kernel_type_i == NA_INTEGER || kernel_type_i < 0) {
        Rf_error("S_estimate_local_density_over_grid(): 'kernel_type' must be a non-negative integer.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_estimate_local_density_over_grid(): 'verbose' must be TRUE/FALSE.");
    }

    const int  grid_size     = grid_size_i;
    const int  kernel_type   = kernel_type_i;
    const bool verbose       = (verbose_i == TRUE);

    // --- Core computation (no R allocations inside) ---
    gdensity_t gdens_res = estimate_local_density_over_grid(
        x, grid_size, poffset, pilot_bandwidth, kernel_type, verbose);

    // --- Build result list (container-first) ---
    const int N_COMPONENTS = 6;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
    SET_STRING_ELT(names, 0, Rf_mkChar("y"));
    SET_STRING_ELT(names, 1, Rf_mkChar("bw"));
    SET_STRING_ELT(names, 2, Rf_mkChar("bw_auto_selected"));
    SET_STRING_ELT(names, 3, Rf_mkChar("offset"));
    SET_STRING_ELT(names, 4, Rf_mkChar("start"));
    SET_STRING_ELT(names, 5, Rf_mkChar("end"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0: y (density)
    {
        SEXP el0 = PROTECT(convert_vector_double_to_R(gdens_res.density));
        SET_VECTOR_ELT(result, 0, el0);
        UNPROTECT(1);
    }
    // 1: bw
    {
        SEXP s_bw = PROTECT(Rf_ScalarReal(gdens_res.bandwidth));
        SET_VECTOR_ELT(result, 1, s_bw);
        UNPROTECT(1);
    }
    // 2: bw_auto_selected
    {
        SEXP s_auto = PROTECT(Rf_ScalarLogical(gdens_res.auto_selected ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 2, s_auto);
        UNPROTECT(1);
    }
    // 3: offset
    {
        SEXP s_off = PROTECT(Rf_ScalarReal(gdens_res.offset));
        SET_VECTOR_ELT(result, 3, s_off);
        UNPROTECT(1);
    }
    // 4: start
    {
        SEXP s_start = PROTECT(Rf_ScalarReal(gdens_res.start));
        SET_VECTOR_ELT(result, 4, s_start);
        UNPROTECT(1);
    }
    // 5: end
    {
        SEXP s_end = PROTECT(Rf_ScalarReal(gdens_res.end));
        SET_VECTOR_ELT(result, 5, s_end);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
