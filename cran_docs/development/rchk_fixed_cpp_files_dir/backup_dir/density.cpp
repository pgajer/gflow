/**
 * @brief Fixed version of density.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_estimate_local_density_over_grid
 * 
 * Issues fixed:
 * 1. S_estimate_local_density_over_grid (lines 463, 528): negative depth, over-/under-protect, stack imbalance
 * 2. Variable UNPROTECT with n_protected tracking
 * 3. Complex conditional PROTECT/UNPROTECT pattern
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

/**
 * Fixed version of S_estimate_local_density_over_grid
 * Fixes: negative depth, stack imbalance at lines 463, 528
 * Solution: Use PROTECT_WITH_INDEX, container-first pattern, literal UNPROTECT
 */
SEXP S_estimate_local_density_over_grid(SEXP s_x,
                                        SEXP s_grid_size,
                                        SEXP s_poffset,
                                        SEXP s_pilot_bandwidth,
                                        SEXP s_kernel_type,
                                        SEXP s_verbose) {
    
    // --- Coerce x to REAL and copy into std::vector<double> ---
    PROTECT_INDEX ipx;
    SEXP sx = s_x;
    PROTECT_WITH_INDEX(sx, &ipx);
    
    if (TYPEOF(sx) != REALSXP) {
        REPROTECT(sx = Rf_coerceVector(sx, REALSXP), ipx);
    }
    
    const R_xlen_t nx = XLENGTH(sx);
    const double* px = REAL(sx);
    std::vector<double> x;
    x.assign(px, px + static_cast<size_t>(nx));
    
    // --- Scalars ---
    const int    grid_size       = Rf_asInteger(s_grid_size);
    const double poffset         = Rf_asReal(s_poffset);
    const double pilot_bandwidth = Rf_asReal(s_pilot_bandwidth);
    const int    kernel_type     = Rf_asInteger(s_kernel_type);
    const bool   verbose         = (Rf_asLogical(s_verbose) == TRUE);

    // --- Core computation (no R allocations inside) ---
    gdensity_t gdens_res = estimate_local_density_over_grid(
        x, grid_size, poffset, pilot_bandwidth, kernel_type, verbose);

    // --- Build result list (container-first) ---
    const int N_COMPONENTS = 6;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));

    // 0: y (density)
    {
        SEXP el0 = PROTECT(convert_vector_double_to_R(gdens_res.density));
        SET_VECTOR_ELT(result, 0, el0);
        UNPROTECT(1); // el0
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

    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("y"));
        SET_STRING_ELT(names, 1, Rf_mkChar("bw"));
        SET_STRING_ELT(names, 2, Rf_mkChar("bw_auto_selected"));
        SET_STRING_ELT(names, 3, Rf_mkChar("offset"));
        SET_STRING_ELT(names, 4, Rf_mkChar("start"));
        SET_STRING_ELT(names, 5, Rf_mkChar("end"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    UNPROTECT(2); // sx, result
    return result;
}

} // extern "C"