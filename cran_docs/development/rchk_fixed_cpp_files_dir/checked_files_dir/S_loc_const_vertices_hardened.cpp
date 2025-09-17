
// Hardened, rchk-safe drop-in for S_loc_const_vertices
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
    // ---- Validate inputs to avoid UB when accessing REAL()/INTEGER() ----
    if (!Rf_isVectorList(r_adj_list)) {
        Rf_error("`adj.list` must be a list (of integer vectors).");
    }
    if (!Rf_isNumeric(Ry)) {
        Rf_error("`y` must be a numeric vector.");
    }
    if (!Rf_isNumeric(Rprec) || Rf_length(Rprec) < 1) {
        Rf_error("`prec` must be a numeric scalar.");
    }

    // Coerce scalar and check positivity as the R wrapper enforces
    double prec = Rf_asReal(Rprec);
    if (!(prec > 0)) {
        Rf_error("`prec` must be positive.");
    }

    // Convert structures
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(r_adj_list);

    const R_xlen_t ny = Rf_xlength(Ry);
    std::vector<double> y;
    y.reserve((size_t)ny);
    const double* yptr = REAL(Ry);
    for (R_xlen_t i = 0; i < ny; ++i) y.push_back(yptr[i]);

    if ((size_t)ny != adj_list.size()) {
        Rf_error("Length of `y` (%ld) must equal number of vertices in the graph (%zu).",
                 (long)ny, adj_list.size());
    }

    // Compute
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

    // Allocate result vector and fill; container protected before assigning
    const R_xlen_t nres = (R_xlen_t)locs->size();
    SEXP result = PROTECT(Rf_allocVector(INTSXP, nres));
    int* rptr = INTEGER(result);
    for (R_xlen_t i = 0; i < nres; ++i) rptr[i] = (*locs)[(size_t)i];

    UNPROTECT(1);
    return result;
}
