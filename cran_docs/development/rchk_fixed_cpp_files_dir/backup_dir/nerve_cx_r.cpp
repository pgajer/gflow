/**
 * @brief Fixed version of nerve_cx_r.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_create_nerve_complex
 * 
 * Issues fixed:
 * 1. S_create_nerve_complex (lines 58, 65): unprotected 'complex_ptr' during allocations
 * 2. complex_ptr not protected when allocating result list
 * 
 * Changes made:
 * 1. Used container-first pattern - protect result list first
 * 2. Protected complex_ptr before setting it in result
 * 3. Fixed PROTECT/UNPROTECT balance
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// External pointer finalizer (assumed available)
extern void nerve_complex_finalizer(SEXP);

// Core class declaration (assumed available)
class nerve_complex_t {
public:
    nerve_complex_t(const std::vector<std::vector<double>>& coords, int k, int max_dim);
    int num_simplices(int dim);
    size_t num_vertices() { return n_vertices; }
private:
    size_t n_vertices;
};

extern "C" {

/**
 * Fixed version of S_create_nerve_complex
 * Fixes: unprotected 'complex_ptr' during allocations at lines 58, 65
 * Solution: Use container-first pattern, protect all allocations properly
 */
SEXP S_create_nerve_complex(SEXP s_coords, SEXP s_k, SEXP s_max_dim) {
    // Convert inputs
    PROTECT_INDEX ipx;
    SEXP coords_real = s_coords;
    PROTECT_WITH_INDEX(coords_real, &ipx);
    if (TYPEOF(coords_real) != REALSXP) {
        REPROTECT(coords_real = Rf_coerceVector(coords_real, REALSXP), ipx);
    }
    
    int* dimX = INTEGER(Rf_getAttrib(coords_real, R_DimSymbol));
    int n_points = dimX[0];
    int n_dims = dimX[1];
    double* X = REAL(coords_real);

    int k = INTEGER(s_k)[0];
    int max_dim = INTEGER(s_max_dim)[0];

    // Convert coordinates to C++ format
    std::vector<std::vector<double>> coords(n_points);
    for (int i = 0; i < n_points; i++) {
        coords[i].resize(n_dims);
        for (int j = 0; j < n_dims; j++) {
            coords[i][j] = X[i + n_points * j]; // Column-major to row-major
        }
    }

    // Create the nerve complex
    nerve_complex_t* complex = new nerve_complex_t(coords, k, max_dim);

    // Create external pointer
    SEXP complex_ptr = PROTECT(R_MakeExternalPtr(complex, R_NilValue, R_NilValue));

    // Register finalizer
    R_RegisterCFinalizerEx(complex_ptr, (R_CFinalizer_t) nerve_complex_finalizer, TRUE);

    // Create return list (container-first pattern)
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));

    // Set complex_ptr as first element
    SET_VECTOR_ELT(result, 0, complex_ptr);

    // Set n_points
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(n_points));

    // Get simplex counts
    {
        SEXP simplex_counts = PROTECT(Rf_allocVector(INTSXP, max_dim + 1));
        int* counts = INTEGER(simplex_counts);
        for (int d = 0; d <= max_dim; d++) {
            counts[d] = complex->num_simplices(d);
        }
        SET_VECTOR_ELT(result, 2, simplex_counts);
        UNPROTECT(1); // simplex_counts
    }

    // Set max_dim
    SET_VECTOR_ELT(result, 3, Rf_ScalarInteger(max_dim));

    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("complex_ptr"));
        SET_STRING_ELT(names, 1, Rf_mkChar("n_vertices"));
        SET_STRING_ELT(names, 2, Rf_mkChar("simplex_counts"));
        SET_STRING_ELT(names, 3, Rf_mkChar("max_dimension"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    UNPROTECT(3); // coords_real, complex_ptr, result
    return result;
}

} // extern "C"