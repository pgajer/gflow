/**
 * @brief Hardened rchk-safe version of graph_spectral_lowess_r.cpp
 *
 * Changes vs. previous corrected version:
 *  - Defensive coercion of s_y using PROTECT_WITH_INDEX/REPROTECT.
 *  - Long-vector safety with XLENGTH() and R_xlen_t.
 *  - Scalars extracted via Rf_asInteger / Rf_asReal / Rf_asLogical.
 *  - Fixed-count UNPROTECT only; container-first pattern preserved.
 */

#include <vector>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers
#include "graph_spectral_lowess.hpp"
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations (as in the existing codebase)
struct graph_spectral_lowess_t {
    std::vector<double> predictions;
    std::vector<double> errors;
    std::vector<double> scale;
};

class uniform_grid_graph_t {
public:
    double graph_diameter;

    uniform_grid_graph_t(const std::vector<std::vector<int>>& adj_list,
                         const std::vector<std::vector<double>>& weight_list);

    std::pair<size_t, double> get_vertex_eccentricity(size_t vertex);

    graph_spectral_lowess_t graph_spectral_lowess(
        const std::vector<double>& y,
        size_t n_evectors,
        size_t n_bws,
        bool log_grid,
        double min_bw_factor,
        double max_bw_factor,
        double dist_normalization_factor,
        size_t kernel_type,
        double precision,
        size_t n_cleveland_iterations,
        bool verbose);
};

extern "C" {

SEXP S_graph_spectral_lowess(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_n_evectors,
    // bw parameters
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    // kernel parameters
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    // other
    SEXP s_precision,
    SEXP s_n_cleveland_iterations,
    SEXP s_verbose
    ) {

    // Convert adjacency and weights (helpers assumed to be pure C++ without R allocations)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Defensive coercion block for y (REAL), long-vector safe
    std::vector<double> y;
    {
        SEXP sy = s_y;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        const R_xlen_t ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + (size_t)ny);
        UNPROTECT(1); // sy
    }

    // Scalars / flags (defensive extraction)
    const size_t n_evectors              = (size_t) Rf_asInteger(s_n_evectors);
    const size_t n_bws                   = (size_t) Rf_asInteger(s_n_bws);
    const bool   log_grid                = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor           = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor           = Rf_asReal(s_max_bw_factor);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type             = (size_t) Rf_asInteger(s_kernel_type);
    const double precision               = Rf_asReal(s_precision);
    const size_t n_cleveland_iterations  = (size_t) Rf_asInteger(s_n_cleveland_iterations);
    const bool   verbose                 = (Rf_asLogical(s_verbose) == TRUE);

    // Build grid graph and compute diameter endpoints
    uniform_grid_graph_t grid_graph(adj_list, weight_list);
    auto ecc1 = grid_graph.get_vertex_eccentricity(0);
    auto ecc2 = grid_graph.get_vertex_eccentricity(ecc1.first);
    grid_graph.graph_diameter = ecc2.second;

    // Run core computation (no R allocations inside)
    graph_spectral_lowess_t res = grid_graph.graph_spectral_lowess(
        y,
        n_evectors,
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        dist_normalization_factor,
        kernel_type,
        precision,
        n_cleveland_iterations,
        verbose
    );

    // Assemble result: list(predictions, errors, scale, graph_diameter)
    const int n_elements = 4;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements)); // [1]

    // names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, n_elements)); // [2]
        SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 1, Rf_mkChar("errors"));
        SET_STRING_ELT(names, 2, Rf_mkChar("scale"));
        SET_STRING_ELT(names, 3, Rf_mkChar("graph_diameter"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names -> [1]
    }

    // 0: predictions
    {
        const R_xlen_t n = (R_xlen_t) res.predictions.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.predictions.begin(), res.predictions.end(), REAL(v));
        SET_VECTOR_ELT(result, 0, v);
        UNPROTECT(1); // v -> [1]
    }

    // 1: errors
    {
        const R_xlen_t n = (R_xlen_t) res.errors.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.errors.begin(), res.errors.end(), REAL(v));
        SET_VECTOR_ELT(result, 1, v);
        UNPROTECT(1); // v -> [1]
    }

    // 2: scale
    {
        const R_xlen_t n = (R_xlen_t) res.scale.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.scale.begin(), res.scale.end(), REAL(v));
        SET_VECTOR_ELT(result, 2, v);
        UNPROTECT(1); // v -> [1]
    }

    // 3: graph_diameter (scalar)
    {
        SEXP v = PROTECT(Rf_allocVector(REALSXP, 1)); // [2]
        REAL(v)[0] = grid_graph.graph_diameter;
        SET_VECTOR_ELT(result, 3, v);
        UNPROTECT(1); // v -> [1]
    }

    UNPROTECT(1); // result
    return result;
}

} // extern "C"
