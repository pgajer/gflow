/**
 * @brief Hardened rchk-safe drop-in of S_graph_bw_adaptive_spectral_smoother
 *
 * Changes vs. provided version:
 * - Use PROTECT_WITH_INDEX/REPROTECT for s_y coercion (container-first).
 * - Use XLENGTH/R_xlen_t for long-vector safety.
 * - Use Rf_asInteger/Rf_asReal/Rf_asLogical for scalar extraction.
 * - Result assembly keeps `result` and `names` protected until the tail UNPROTECT(2).
 * - All UNPROTECT calls are constant literals.
 */
#include <vector>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers (assumed present in the build)
#include "set_wgraph.hpp"
#include "graph_bw_adaptive_spectral_smoother.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations (assumed available at link time)
struct graph_bw_adaptive_spectral_smoother_t {
    std::vector<double> predictions;
    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> bw_mean_abs_errors;
    std::vector<double> vertex_min_bws;
    size_t opt_bw_idx;
};

class set_wgraph_t {
public:
    set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                 const std::vector<std::vector<double>>& weight_list);
    
    graph_bw_adaptive_spectral_smoother_t graph_bw_adaptive_spectral_smoother(
        const std::vector<double>& y,
        size_t n_evectors,
        double min_bw_factor,
        double max_bw_factor,
        size_t n_bws,
        bool log_grid,
        size_t kernel_type,
        double dist_normalization_factor,
        double precision,
        bool use_global_bw_grid,
        bool with_bw_predictions,
        bool with_vertex_bw_errors,
        bool verbose);
};

extern "C" {

SEXP S_graph_bw_adaptive_spectral_smoother(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_n_evectors,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_precision,
    SEXP s_use_global_bw_grid,
    SEXP s_with_bw_predictions,
    SEXP s_with_vertex_bw_errors,
    SEXP s_verbose
) {
    // Convert adjacency & weight lists via trusted helpers
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // --- Coercion block for y (container-first + indexed protect) ---
    std::vector<double> y;
    {
        SEXP sy = s_y;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        const R_xlen_t ny = XLENGTH(sy);
        if (ny <= 0) {
            UNPROTECT(1); // sy
            Rf_error("S_graph_bw_adaptive_spectral_smoother: 'y' must be a non-empty numeric vector");
        }
        y.assign(REAL(sy), REAL(sy) + (size_t)ny);
        UNPROTECT(1); // sy
    }

    // Defensive scalar extraction
    const size_t n_evectors = (size_t) Rf_asInteger(s_n_evectors);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);
    const size_t n_bws = (size_t) Rf_asInteger(s_n_bws);
    const bool log_grid = Rf_asLogical(s_log_grid) == TRUE;
    const size_t kernel_type = (size_t) Rf_asInteger(s_kernel_type);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const double precision = Rf_asReal(s_precision);
    const bool use_global_bw_grid = Rf_asLogical(s_use_global_bw_grid) == TRUE;
    const bool with_bw_predictions = Rf_asLogical(s_with_bw_predictions) == TRUE;
    const bool with_vertex_bw_errors = Rf_asLogical(s_with_vertex_bw_errors) == TRUE;
    const bool verbose = Rf_asLogical(s_verbose) == TRUE;

    // Optionally add simple validations (no PROTECTs held here)
    if (n_bws == 0) Rf_error("S_graph_bw_adaptive_spectral_smoother: 'n_bws' must be >= 1");
    if (!(min_bw_factor > 0.0) || !(max_bw_factor > 0.0) || !(max_bw_factor >= min_bw_factor)) {
        Rf_error("S_graph_bw_adaptive_spectral_smoother: invalid bw factors (ensure 0 < min <= max)");
    }

    // Compute
    set_wgraph_t graph(adj_list, weight_list);
    graph_bw_adaptive_spectral_smoother_t result =
        graph.graph_bw_adaptive_spectral_smoother(
            y,
            n_evectors,
            min_bw_factor,
            max_bw_factor,
            n_bws,
            log_grid,
            kernel_type,
            dist_normalization_factor,
            precision,
            use_global_bw_grid,
            with_bw_predictions,
            with_vertex_bw_errors,
            verbose
        );

    // --- Result assembly (fixed-count UNPROTECT) ---
    const int N = 5;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N));

    // 0: predictions
    {
        SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) result.predictions.size()));
        std::copy(result.predictions.begin(), result.predictions.end(), REAL(x));
        SET_VECTOR_ELT(r_result, 0, x);
        UNPROTECT(1);
    }

    // 1: bw_predictions (matrix, column-major)
    {
        size_t ncol = result.bw_predictions.size();
        size_t nrow = ncol ? result.bw_predictions[0].size() : 0;
        SEXP x = PROTECT(Rf_allocMatrix(REALSXP, (R_xlen_t) nrow, (R_xlen_t) ncol));
        double* ptr = REAL(x);
        for (size_t j = 0; j < ncol; ++j) {
            if (result.bw_predictions[j].size() != nrow) {
                UNPROTECT(1); // x
                Rf_error("S_graph_bw_adaptive_spectral_smoother: inconsistent 'bw_predictions' dimensions");
            }
            for (size_t i = 0; i < nrow; ++i) {
                ptr[i + j * nrow] = result.bw_predictions[j][i];
            }
        }
        SET_VECTOR_ELT(r_result, 1, x);
        UNPROTECT(1);
    }

    // 2: bw_mean_abs_errors
    {
        SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) result.bw_mean_abs_errors.size()));
        std::copy(result.bw_mean_abs_errors.begin(), result.bw_mean_abs_errors.end(), REAL(x));
        SET_VECTOR_ELT(r_result, 2, x);
        UNPROTECT(1);
    }

    // 3: vertex_min_bws
    {
        SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) result.vertex_min_bws.size()));
        std::copy(result.vertex_min_bws.begin(), result.vertex_min_bws.end(), REAL(x));
        SET_VECTOR_ELT(r_result, 3, x);
        UNPROTECT(1);
    }

    // 4: opt_bw_idx (1-based for R)
    {
        SEXP x = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(x)[0] = (int) result.opt_bw_idx + 1;
        SET_VECTOR_ELT(r_result, 4, x);
        UNPROTECT(1);
    }

    // names (keep protected until tail)
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(r_names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("bw_predictions"));
    SET_STRING_ELT(r_names, 2, Rf_mkChar("bw_mean_abs_errors"));
    SET_STRING_ELT(r_names, 3, Rf_mkChar("vertex_min_bws"));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("opt_bw_idx"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(2); // r_result, r_names
    return r_result;
}

} // extern "C"
