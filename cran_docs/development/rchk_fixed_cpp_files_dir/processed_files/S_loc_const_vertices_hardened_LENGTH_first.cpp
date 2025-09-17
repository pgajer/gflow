
// Hardened, rchk-safe drop-in for S_loc_const_vertices
// Policy: Pragmatic LENGTH-first (no R_xlen_t / XLENGTH).
// - Type checks and simple coercion blocks
// - Fixed-count UNPROTECTs (no variable UNPROTECT)
// - Uses LENGTH() + int (R_len_t) for sizes; size_t only for STL
// - Validates y length == number of vertices
// - try/catch around compute; clear error messages

#include <Rinternals.h>
#include <R_ext/Error.h>
#include <vector>
#include <memory>
#include <stdexcept>

// Forward decls expected from the package
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::unique_ptr<std::vector<int>> loc_const_vertices(const std::vector<std::vector<int>>& adj_list,
                                                     const std::vector<double>& y,
                                                     double prec);

extern "C" SEXP S_loc_const_vertices(SEXP r_adj_list, SEXP Ry, SEXP Rprec) {
    // ---- Basic validation (cheap checks before any conversion) ----
    if (!Rf_isVectorList(r_adj_list)) {
        Rf_error("`adj.list` must be a list (of integer vectors).");
    }
    if (!Rf_isNumeric(Ry)) {
        Rf_error("`y` must be a numeric vector.");
    }
    if (!Rf_isNumeric(Rprec) || Rf_length(Rprec) < 1) {
        Rf_error("`prec` must be a numeric scalar.");
    }

    // ---- Scalars ----
    double prec = Rf_asReal(Rprec);
    if (!(prec > 0)) {
        Rf_error("`prec` must be positive.");
    }

    // ---- Convert adjacency ----
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(r_adj_list);

    // ---- Coercion block for y (REAL) with container-first protection ----
    SEXP sy = Ry;
    PROTECT_INDEX py;
    PROTECT_WITH_INDEX(sy, &py);
    if (TYPEOF(sy) != REALSXP) {
        REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
    }
    const int ny = LENGTH(sy); // LENGTH-first policy
    if (ny < 0) {
        UNPROTECT(1);
        Rf_error("Negative length for `y`.");
    }
    // Copy while protected
    std::vector<double> y(REAL(sy), REAL(sy) + (size_t)ny);
    UNPROTECT(1); // sy

    // ---- Cross-check sizes ----
    if ((size_t)ny != adj_list.size()) {
        Rf_error("Length of `y` (%d) must equal number of vertices (%zu).", ny, adj_list.size());
    }

    // ---- Compute ----
    std::unique_ptr<std::vector<int>> locs;
    try {
        locs = loc_const_vertices(adj_list, y, prec);
    } catch (const std::exception& e) {
        Rf_error("loc_const_vertices() failed: %s", e.what());
    } catch (...) {
        Rf_error("loc_const_vertices() failed with an unknown error.");
    }
    if (!locs) {
        Rf_error("Internal error: computation returned null result.");
    }

    // ---- Result assembly (fixed UNPROTECT pattern) ----
    const int nres = (int)locs->size(); // safe: we do not support long vectors by policy
    SEXP result = PROTECT(Rf_allocVector(INTSXP, nres));
    int* rptr = INTEGER(result);
    for (int i = 0; i < nres; ++i) rptr[i] = (*locs)[(size_t)i];

    UNPROTECT(1);
    return result;
}
