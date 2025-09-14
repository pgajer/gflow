/**
 * @brief Hardened rchk-safe drop-in of S_create_nerve_complex
 *
 * Changes vs. provided version:
 * - Keep `result` and `names` protected until the tail; finish with UNPROTECT(4) (coords_real, complex_ptr, result, names).
 * - Defensive scalar parsing via Rf_asInteger; basic dimension validation.
 * - Maintain container-first assembly with local PROTECT/UNPROTECT(1) blocks.
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

SEXP S_create_nerve_complex(SEXP s_coords, SEXP s_k, SEXP s_max_dim) {
    // --- Coerce coords to REAL matrix (indexed protect) ---
    PROTECT_INDEX px;
    SEXP coords_real = s_coords;
    PROTECT_WITH_INDEX(coords_real, &px);
    if (TYPEOF(coords_real) != REALSXP) {
        REPROTECT(coords_real = Rf_coerceVector(coords_real, REALSXP), px);
    }

    // Validate and read dimensions defensively
    SEXP dims = Rf_getAttrib(coords_real, R_DimSymbol);
    if (Rf_isNull(dims) || LENGTH(dims) != 2 || TYPEOF(dims) != INTSXP) {
        UNPROTECT(1); // coords_real
        Rf_error("S_create_nerve_complex: 'coords' must be a numeric matrix");
    }
    const int n_points = INTEGER(dims)[0];
    const int n_dims   = INTEGER(dims)[1];
    if (n_points <= 0 || n_dims <= 0) {
        UNPROTECT(1); // coords_real
        Rf_error("S_create_nerve_complex: invalid matrix dimensions");
    }
    const double* X = REAL(coords_real);

    // Defensive scalars
    const int k       = Rf_asInteger(s_k);
    const int max_dim = Rf_asInteger(s_max_dim);
    if (k < 1) {
        UNPROTECT(1); // coords_real
        Rf_error("S_create_nerve_complex: 'k' must be >= 1");
    }
    if (max_dim < 0) {
        UNPROTECT(1); // coords_real
        Rf_error("S_create_nerve_complex: 'max_dim' must be >= 0");
    }

    // Convert coordinates to C++ format (col-major -> row-major)
    std::vector<std::vector<double>> coords((size_t)n_points, std::vector<double>((size_t)n_dims));
    for (int i = 0; i < n_points; ++i) {
        for (int j = 0; j < n_dims; ++j) {
            coords[(size_t)i][(size_t)j] = X[i + n_points * j];
        }
    }

    // Create the nerve complex (no R allocations here)
    nerve_complex_t* complex = new nerve_complex_t(coords, k, max_dim);

    // External pointer (protected)
    SEXP complex_ptr = PROTECT(R_MakeExternalPtr((void*)complex, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(complex_ptr, (R_CFinalizer_t) nerve_complex_finalizer, TRUE);

    // Result container (keep protected until tail)
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));

    // 0: complex_ptr
    SET_VECTOR_ELT(result, 0, complex_ptr);

    // 1: n_vertices (expose number of vertices; fall back to n_points if method not available at this stage)
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(n_points));

    // 2: simplex_counts (0..max_dim)
    {
        SEXP simplex_counts = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)(max_dim + 1)));
        int* counts = INTEGER(simplex_counts);
        for (int d = 0; d <= max_dim; ++d) {
            counts[d] = complex->num_simplices(d);
        }
        SET_VECTOR_ELT(result, 2, simplex_counts);
        UNPROTECT(1); // simplex_counts
    }

    // 3: max_dimension
    SET_VECTOR_ELT(result, 3, Rf_ScalarInteger(max_dim));

    // names — keep protected until tail
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("complex_ptr"));
    SET_STRING_ELT(names, 1, Rf_mkChar("n_vertices"));
    SET_STRING_ELT(names, 2, Rf_mkChar("simplex_counts"));
    SET_STRING_ELT(names, 3, Rf_mkChar("max_dimension"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(4); // coords_real, complex_ptr, result, names
    return result;
}

} // extern "C"
