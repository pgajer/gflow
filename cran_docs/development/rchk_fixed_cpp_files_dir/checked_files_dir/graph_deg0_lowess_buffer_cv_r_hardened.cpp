// graph_deg0_lowess_buffer_cv_r_hardened.cpp — rchk-hardened drop-in
//
// Hardenings applied:
// - Defensive scalar coercion via Rf_as* helpers (no direct REAL/INTEGER/LOGICAL indexing)
// - NA and basic range checks for all scalars
// - Long-vector safety for y (XLENGTH(), explicit casts)
// - Coercion block with PROTECT_WITH_INDEX for s_y (copy to std::vector before UNPROTECT)
// - Result assembly uses names protected until function tail (UNPROTECT(2) final)
// - Container-first: convert adj/weights and y before any result allocations
//
// Assumptions:
// - convert_adj_list_from_R(s_adj_list) and convert_weight_list_from_R(s_weight_list) exist and are safe.
// - set_wgraph_t and graph_deg0_lowess_buffer_cv_t interfaces match the declarations below.
// - kernel_type validity domain is runtime-defined; we only check non-negativity.

#include <vector>
#include <algorithm>
#include <cmath>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Include necessary headers (as in original file)
#include "set_wgraph.hpp"
#include "graph_deg0_lowess_buffer_cv.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core result structure — mirrored from original
struct graph_deg0_lowess_buffer_cv_t {
    std::vector<double> predictions;
    std::vector<std::vector<double>> bw_predictions;
    std::vector<double> bw_errors;
    std::vector<double> bws;
    double opt_bw;
    size_t opt_bw_idx;
    size_t buffer_hops_used;
};

class set_wgraph_t {
public:
    set_wgraph_t(const std::vector<std::vector<int>>& adj_list,
                 const std::vector<std::vector<double>>& weight_list);
    
    graph_deg0_lowess_buffer_cv_t graph_deg0_lowess_buffer_cv(
        const std::vector<double>& y,
        double min_bw_factor,
        double max_bw_factor,
        size_t n_bws,
        bool log_grid,
        size_t kernel_type,
        double dist_normalization_factor,
        bool use_uniform_weights,
        size_t buffer_hops,
        bool auto_buffer_hops,
        size_t n_folds,
        bool with_bw_predictions,
        double precision,
        bool verbose);
};

extern "C" {

SEXP S_graph_deg0_lowess_buffer_cv(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_use_uniform_weights,
    SEXP s_buffer_hops,
    SEXP s_auto_buffer_hops,
    SEXP s_n_folds,
    SEXP s_with_bw_predictions,
    SEXP s_precision,
    SEXP s_verbose) 
{
    // --- Container-first conversions ---
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // s_y: allow type-flexible input, coerce if needed, copy out, then unprotect
    std::vector<double> y;
    {
        SEXP sy = s_y;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) {
            REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        }
        const R_xlen_t ny = XLENGTH(sy);
        if (ny <= 0) {
            UNPROTECT(1);
            Rf_error("S_graph_deg0_lowess_buffer_cv(): 'y' must be a non-empty numeric vector.");
        }
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));
        UNPROTECT(1); // sy
    }

    // --- Defensive scalar coercions ---
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);
    const int n_bws_i = Rf_asInteger(s_n_bws);
    const int log_grid_i = Rf_asLogical(s_log_grid);
    const int kernel_type_i = Rf_asInteger(s_kernel_type);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const int use_uniform_weights_i = Rf_asLogical(s_use_uniform_weights);
    const int buffer_hops_i = Rf_asInteger(s_buffer_hops);
    const int auto_buffer_hops_i = Rf_asLogical(s_auto_buffer_hops);
    const int n_folds_i = Rf_asInteger(s_n_folds);
    const int with_bw_predictions_i = Rf_asLogical(s_with_bw_predictions);
    const double precision = Rf_asReal(s_precision);
    const int verbose_i = Rf_asLogical(s_verbose);

    // --- NA / range checks ---
    if (ISNAN(min_bw_factor) || ISNAN(max_bw_factor)) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): min_bw_factor/max_bw_factor cannot be NA.");
    }
    if (min_bw_factor <= 0.0 || max_bw_factor <= 0.0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): min_bw_factor and max_bw_factor must be > 0.");
    }
    if (max_bw_factor < min_bw_factor) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): max_bw_factor must be >= min_bw_factor.");
    }
    if (n_bws_i == NA_INTEGER || n_bws_i <= 0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): n_bws must be a positive integer.");
    }
    if (log_grid_i == NA_LOGICAL) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): log_grid must be TRUE/FALSE.");
    }
    if (kernel_type_i == NA_INTEGER || kernel_type_i < 0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): kernel_type must be a non-negative integer.");
    }
    if (ISNAN(dist_normalization_factor) || dist_normalization_factor < 0.0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): dist_normalization_factor must be >= 0.");
    }
    if (use_uniform_weights_i == NA_LOGICAL) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): use_uniform_weights must be TRUE/FALSE.");
    }
    if (buffer_hops_i == NA_INTEGER || buffer_hops_i < 0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): buffer_hops must be a non-negative integer.");
    }
    if (auto_buffer_hops_i == NA_LOGICAL) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): auto_buffer_hops must be TRUE/FALSE.");
    }
    if (n_folds_i == NA_INTEGER || n_folds_i <= 0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): n_folds must be a positive integer.");
    }
    if (with_bw_predictions_i == NA_LOGICAL) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): with_bw_predictions must be TRUE/FALSE.");
    }
    if (ISNAN(precision) || precision <= 0.0) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): precision must be > 0.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_graph_deg0_lowess_buffer_cv(): verbose must be TRUE/FALSE.");
    }

    // Cast to runtime types
    const size_t n_bws = static_cast<size_t>(n_bws_i);
    const bool log_grid = (log_grid_i == TRUE);
    const size_t kernel_type = static_cast<size_t>(kernel_type_i);
    const bool use_uniform_weights = (use_uniform_weights_i == TRUE);
    const size_t buffer_hops = static_cast<size_t>(buffer_hops_i);
    const bool auto_buffer_hops = (auto_buffer_hops_i == TRUE);
    const size_t n_folds = static_cast<size_t>(n_folds_i);
    const bool with_bw_predictions = (with_bw_predictions_i == TRUE);
    const bool verbose = (verbose_i == TRUE);

    // --- Core call ---
    set_wgraph_t graph(adj_list, weight_list);
    graph_deg0_lowess_buffer_cv_t result =
        graph.graph_deg0_lowess_buffer_cv(
            y,
            min_bw_factor,
            max_bw_factor,
            n_bws,
            log_grid,
            kernel_type,
            dist_normalization_factor,
            use_uniform_weights,
            buffer_hops,
            auto_buffer_hops,
            n_folds,
            with_bw_predictions,
            precision,
            verbose
        );

    // --- Result assembly (list + names) ---
    const int n_elements = 7;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    SEXP r_names  = PROTECT(Rf_allocVector(STRSXP, n_elements));

    // names
    SET_STRING_ELT(r_names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("bw_predictions"));
    SET_STRING_ELT(r_names, 2, Rf_mkChar("bw_errors"));
    SET_STRING_ELT(r_names, 3, Rf_mkChar("bws"));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("opt_bw"));
    SET_STRING_ELT(r_names, 5, Rf_mkChar("opt_bw_idx"));
    SET_STRING_ELT(r_names, 6, Rf_mkChar("buffer_hops_used"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    // 0: predictions
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)result.predictions.size()));
        std::copy(result.predictions.begin(), result.predictions.end(), REAL(r_vec));
        SET_VECTOR_ELT(r_result, 0, r_vec);
        UNPROTECT(1);
    }

    // 1: bw_predictions (may be empty)
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        SEXP r_bw_predictions = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t)result.bw_predictions.size()));
        for (size_t i = 0; i < result.bw_predictions.size(); ++i) {
            const std::vector<double>& row = result.bw_predictions[i];
            SEXP r_row = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)row.size()));
            std::copy(row.begin(), row.end(), REAL(r_row));
            SET_VECTOR_ELT(r_bw_predictions, (R_xlen_t)i, r_row);
            UNPROTECT(1); // r_row
        }
        SET_VECTOR_ELT(r_result, 1, r_bw_predictions);
        UNPROTECT(1); // r_bw_predictions
    } else {
        SET_VECTOR_ELT(r_result, 1, Rf_allocVector(VECSXP, 0)); // empty list
    }

    // 2: bw_errors
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)result.bw_errors.size()));
        std::copy(result.bw_errors.begin(), result.bw_errors.end(), REAL(r_vec));
        SET_VECTOR_ELT(r_result, 2, r_vec);
        UNPROTECT(1);
    }

    // 3: bws
    {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)result.bws.size()));
        std::copy(result.bws.begin(), result.bws.end(), REAL(r_vec));
        SET_VECTOR_ELT(r_result, 3, r_vec);
        UNPROTECT(1);
    }

    // 4: opt_bw
    {
        SEXP r_opt_bw = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(r_opt_bw)[0] = result.opt_bw;
        SET_VECTOR_ELT(r_result, 4, r_opt_bw);
        UNPROTECT(1);
    }

    // 5: opt_bw_idx (1-based for R)
    {
        SEXP r_opt_bw_idx = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_opt_bw_idx)[0] = (int)(result.opt_bw_idx + 1);
        SET_VECTOR_ELT(r_result, 5, r_opt_bw_idx);
        UNPROTECT(1);
    }

    // 6: buffer_hops_used
    {
        SEXP r_bhu = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(r_bhu)[0] = (int)result.buffer_hops_used;
        SET_VECTOR_ELT(r_result, 6, r_bhu);
        UNPROTECT(1);
    }

    UNPROTECT(2); // r_result, r_names
    return r_result;
}

} // extern "C"
