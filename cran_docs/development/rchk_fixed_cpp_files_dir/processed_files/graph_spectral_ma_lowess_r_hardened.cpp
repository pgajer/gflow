#include <vector>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers
#include "set_wgraph.hpp"
#include "graph_spectral_lowess.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations
struct graph_spectral_lowess_t {
    std::vector<double> predictions;
    std::vector<double> errors;
    std::vector<double> scale;
};

class set_wgraph_t {
public:
    double graph_diameter;

    set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                 const std::vector<std::vector<double>>& weight_list);

    std::pair<size_t, double> get_vertex_eccentricity(size_t vertex);

    graph_spectral_lowess_t graph_spectral_ma_lowess(
        const std::vector<double>& y,
        size_t n_evectors,
        size_t n_bws,
        bool log_grid,
        double min_bw_factor,
        double max_bw_factor,
        double dist_normalization_factor,
        size_t kernel_type,
        double blending_coef,
        double precision,
        bool verbose);
};

extern "C" {

/**
 * Hardened rchk-clean drop-in for S_graph_spectral_ma_lowess
 * - Uses PROTECT_WITH_INDEX + REPROTECT for y coercion
 * - Uses XLENGTH/R_xlen_t for long-vector safety
 * - Uses Rf_as* for scalar extraction
 * - Container-first result assembly with fixed-count UNPROTECTs
 */
SEXP S_graph_spectral_ma_lowess(
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
    // model parameters
    SEXP s_blending_coef,
    // other
    SEXP s_precision,
    SEXP s_verbose
) {
    // -------- Convert adjacency/weights (helpers assumed to avoid allocations on their own) --------
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // -------- Coerce y to REAL and copy (long-vector safe) --------
    std::vector<double> y;
    {
        SEXP sy = s_y;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        const R_xlen_t ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));
        UNPROTECT(1); // sy
    }

    // -------- Scalars / flags (defensive extraction) --------
    const size_t n_evectors = static_cast<size_t>(Rf_asInteger(s_n_evectors));

    // bw parameters
    const size_t n_bws      = static_cast<size_t>(Rf_asInteger(s_n_bws));
    const bool   log_grid   = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);

    // kernel parameters
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type = static_cast<size_t>(Rf_asInteger(s_kernel_type));

    // model parameters
    const double blending_coef = Rf_asReal(s_blending_coef);

    // other
    const double precision = Rf_asReal(s_precision);
    const bool   verbose   = (Rf_asLogical(s_verbose) == TRUE);

    // -------- Create the graph --------
    set_wgraph_t graph(adj_list, weight_list);

    // -------- Ensure diameter is computed --------
    if (graph.graph_diameter <= 0) {
        auto [end1, diam] = graph.get_vertex_eccentricity(0);  // Start from vertex 0
        (void)diam;
        auto [end2, diameter] = graph.get_vertex_eccentricity(end1);
        (void)end2;
        graph.graph_diameter = diameter;
    }

    // -------- Core computation (no R allocations inside) --------
    graph_spectral_lowess_t res = graph.graph_spectral_ma_lowess(
        y,
        n_evectors,
        // bw parameters
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        // kernel parameters
        dist_normalization_factor,
        kernel_type,
        // model parameters
        blending_coef,
        // other
        precision,
        verbose
    );

    // -------- Build result (container-first; fixed UNPROTECT counts) --------
    const int n_elements = 4;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));

    // names
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    SET_STRING_ELT(result_names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(result_names, 1, Rf_mkChar("errors"));
    SET_STRING_ELT(result_names, 2, Rf_mkChar("scale"));
    SET_STRING_ELT(result_names, 3, Rf_mkChar("graph_diameter"));
    Rf_setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(1); // result_names

    // 0: predictions
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)res.predictions.size()));
        double* ptr = REAL(r_vec);
        std::copy(res.predictions.begin(), res.predictions.end(), ptr);
        SET_VECTOR_ELT(result, 0, r_vec);
        UNPROTECT(1);
    }

    // 1: errors
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)res.errors.size()));
        double* ptr = REAL(r_vec);
        std::copy(res.errors.begin(), res.errors.end(), ptr);
        SET_VECTOR_ELT(result, 1, r_vec);
        UNPROTECT(1);
    }

    // 2: scale
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)res.scale.size()));
        double* ptr = REAL(r_vec);
        std::copy(res.scale.begin(), res.scale.end(), ptr);
        SET_VECTOR_ELT(result, 2, r_vec);
        UNPROTECT(1);
    }

    // 3: graph_diameter
    {
        SEXP diam = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(diam)[0] = graph.graph_diameter;
        SET_VECTOR_ELT(result, 3, diam);
        UNPROTECT(1);
    }

    UNPROTECT(1); // result
    return result;
}

} // extern "C"
