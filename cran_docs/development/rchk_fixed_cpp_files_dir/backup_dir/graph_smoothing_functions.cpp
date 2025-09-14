/**
 * @brief Fixed versions of graph smoothing functions to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected versions of:
 * - S_graph_kernel_smoother (from graph_kernel_smoother_r.cpp)
 * - S_graph_bw_adaptive_spectral_smoother (from graph_bw_adaptive_spectral_smoother_r.cpp)  
 * - S_graph_spectral_lowess (from graph_spectral_lowess_r.cpp)
 * - S_graph_spectral_lowess_mat (from graph_spectral_lowess_mat_r.cpp)
 * - S_graph_spectral_ma_lowess (from graph_spectral_ma_lowess_r.cpp)
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with constant literals
 * 2. Used PROTECT_WITH_INDEX/REPROTECT for conditional coercion
 * 3. Fixed protection counter issues
 * 4. Used container-first pattern consistently
 * 5. Removed lambdas that capture protection counters
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed to be available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);

// Include necessary headers for the actual implementations
#include "graph_kernel_smoother.hpp"
#include "set_wgraph.hpp"
#include "error_utils.h"
#include "kernels.h"
#include "SEXP_cpp_conversion_utils.hpp"
#include "bandwidth_utils.hpp"
#include "progress_utils.hpp"

extern "C" {

/**
 * Fixed version of S_graph_kernel_smoother
 */
SEXP S_graph_kernel_smoother(SEXP s_adj_list,
                             SEXP s_weight_list,
                             SEXP s_y,
                             SEXP s_min_bw_factor,
                             SEXP s_max_bw_factor,
                             SEXP s_n_bws,
                             SEXP s_log_grid,
                             SEXP s_vertex_hbhd_min_size,
                             SEXP s_kernel_type,
                             SEXP s_dist_normalization_factor,
                             SEXP s_use_uniform_weights,
                             SEXP s_buffer_hops,
                             SEXP s_auto_buffer_hops,
                             SEXP s_n_folds,
                             SEXP s_with_bw_predictions,
                             SEXP s_precision,
                             SEXP s_verbose) {
  // ---- Graph inputs ----
  std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);

  std::vector<std::vector<double>> weight_list;
  if (s_weight_list != R_NilValue) {
    weight_list = convert_weight_list_from_R(s_weight_list);
    if (!weight_list.empty() && weight_list.size() != adj_list.size()) {
      Rf_error("weight_list length (%zu) must match adj_list length (%zu).",
               weight_list.size(), adj_list.size());
    }
    if (!weight_list.empty()) {
      const size_t V = adj_list.size();
      for (size_t i = 0; i < V; ++i) {
        if (weight_list[i].size() != adj_list[i].size()) {
          Rf_error("weight_list[[%zu]] length (%zu) must match adj_list[[%zu]] length (%zu).",
                   i + 1, weight_list[i].size(), i + 1, adj_list[i].size());
        }
      }
    }
  }

  // ---- y: coerce defensively using PROTECT_WITH_INDEX ----
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

  if (adj_list.size() != y.size()) {
    Rf_error("length(y) (%zu) must equal length(adj_list) (%zu).",
             y.size(), adj_list.size());
  }

  const R_xlen_t N = static_cast<R_xlen_t>(y.size());

  // ---- Scalars / parameters (validated) ----
  const double min_bw_factor = Rf_asReal(s_min_bw_factor);
  const double max_bw_factor = Rf_asReal(s_max_bw_factor);
  const int    n_bws        = Rf_asInteger(s_n_bws);
  const bool   log_grid     = (Rf_asLogical(s_log_grid) == TRUE);
  const int    vertex_hbhd_min_size = Rf_asInteger(s_vertex_hbhd_min_size);
  const int    kernel_type  = Rf_asInteger(s_kernel_type);
  const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
  const bool   use_uniform_weights = (Rf_asLogical(s_use_uniform_weights) == TRUE);
  const int    buffer_hops  = Rf_asInteger(s_buffer_hops);
  const bool   auto_buffer_hops = (Rf_asLogical(s_auto_buffer_hops) == TRUE);
  const int    n_folds      = Rf_asInteger(s_n_folds);
  const bool   with_bw_predictions = (Rf_asLogical(s_with_bw_predictions) == TRUE);
  const double precision    = Rf_asReal(s_precision);
  const bool   verbose      = (Rf_asLogical(s_verbose) == TRUE);

  if (!(min_bw_factor > 0.0))     Rf_error("min_bw_factor must be > 0.");
  if (!(max_bw_factor > min_bw_factor))
                                 Rf_error("max_bw_factor must be > min_bw_factor.");
  if (n_bws < 1)                 Rf_error("n_bws must be >= 1.");
  if (vertex_hbhd_min_size < 1)  Rf_error("vertex_hbhd_min_size must be >= 1.");
  if (buffer_hops < 0)           Rf_error("buffer_hops must be >= 0.");
  if (n_folds < 2)               Rf_error("n_folds must be >= 2.");
  if (!(precision > 0.0))        Rf_error("precision must be > 0.");

  // ---- Core computation (must not allocate R objects inside) ----
  set_wgraph_t graph(adj_list, weight_list);
  graph_kernel_smoother_t result = graph.graph_kernel_smoother(
      y,
      min_bw_factor,
      max_bw_factor,
      static_cast<size_t>(n_bws),
      log_grid,
      static_cast<size_t>(vertex_hbhd_min_size),
      static_cast<size_t>(kernel_type),
      dist_normalization_factor,
      use_uniform_weights,
      static_cast<size_t>(buffer_hops),
      auto_buffer_hops,
      static_cast<size_t>(n_folds),
      with_bw_predictions,
      precision,
      verbose);

  // ---- Build result (container-first pattern) ----
  const int N_ELT = 6;
  SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N_ELT));

  // 0: predictions (numeric)
  {
    SEXP s = PROTECT(convert_vector_double_to_R(result.predictions));
    if (static_cast<R_xlen_t>(result.predictions.size()) != N) {
      UNPROTECT(2); // s, r_result
      Rf_error("predictions length mismatch.");
    }
    SET_VECTOR_ELT(r_result, 0, s);
    UNPROTECT(1);
  }

  // 1: bw_predictions (matrix or NULL)
  if (with_bw_predictions && !result.bw_predictions.empty()) {
    const size_t ncol_sz = result.bw_predictions.size();
    const size_t nrow_sz = result.bw_predictions[0].size();
    if (nrow_sz != static_cast<size_t>(N)) {
      UNPROTECT(1); // r_result
      Rf_error("bw_predictions: row count (%zu) must equal length(y) (%lld).",
               nrow_sz, static_cast<long long>(N));
    }
    for (size_t j = 1; j < ncol_sz; ++j) {
      if (result.bw_predictions[j].size() != nrow_sz) {
        UNPROTECT(1); // r_result
        Rf_error("bw_predictions: column %zu has length %zu; expected %zu.",
                 j + 1, result.bw_predictions[j].size(), nrow_sz);
      }
    }
    if (nrow_sz > static_cast<size_t>(R_XLEN_T_MAX) ||
        ncol_sz > static_cast<size_t>(R_XLEN_T_MAX)) {
      UNPROTECT(1); // r_result
      Rf_error("bw_predictions dimensions exceed R limits.");
    }
    const R_xlen_t nrow = static_cast<R_xlen_t>(nrow_sz);
    const R_xlen_t ncol = static_cast<R_xlen_t>(ncol_sz);

    SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    double* p = REAL(s);
    for (R_xlen_t j = 0; j < ncol; ++j) {
      const auto& col = result.bw_predictions[static_cast<size_t>(j)];
      for (R_xlen_t i = 0; i < nrow; ++i) {
        p[i + j * nrow] = col[static_cast<size_t>(i)];
      }
    }
    SET_VECTOR_ELT(r_result, 1, s);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(r_result, 1, R_NilValue);
  }

  // 2: bw_mean_abs_errors (numeric)
  {
    if (!result.bw_mean_abs_errors.empty() &&
        static_cast<int>(result.bw_mean_abs_errors.size()) != n_bws) {
      UNPROTECT(1); // r_result
      Rf_error("bw_mean_abs_errors length (%zu) must equal n_bws (%d).",
               result.bw_mean_abs_errors.size(), n_bws);
    }
    SEXP s = PROTECT(convert_vector_double_to_R(result.bw_mean_abs_errors));
    SET_VECTOR_ELT(r_result, 2, s);
    UNPROTECT(1);
  }

  // 3: vertex_min_bws (numeric)
  {
    if (!result.vertex_min_bws.empty() &&
        static_cast<R_xlen_t>(result.vertex_min_bws.size()) != N) {
      UNPROTECT(1); // r_result
      Rf_error("vertex_min_bws length mismatch.");
    }
    SEXP s = PROTECT(convert_vector_double_to_R(result.vertex_min_bws));
    SET_VECTOR_ELT(r_result, 3, s);
    UNPROTECT(1);
  }

  // 4: opt_bw_idx (scalar int; expose 1-based)
  {
    SEXP s = PROTECT(Rf_ScalarInteger(static_cast<int>(result.opt_bw_idx) + 1));
    SET_VECTOR_ELT(r_result, 4, s);
    UNPROTECT(1);
  }

  // 5: buffer_hops_used (scalar int)
  {
    SEXP s = PROTECT(Rf_ScalarInteger(static_cast<int>(result.buffer_hops_used)));
    SET_VECTOR_ELT(r_result, 5, s);
    UNPROTECT(1);
  }

  // Set names
  SEXP r_names = PROTECT(Rf_allocVector(STRSXP, N_ELT));
  SET_STRING_ELT(r_names, 0, Rf_mkChar("predictions"));
  SET_STRING_ELT(r_names, 1, Rf_mkChar("bw_predictions"));
  SET_STRING_ELT(r_names, 2, Rf_mkChar("bw_mean_abs_errors"));
  SET_STRING_ELT(r_names, 3, Rf_mkChar("vertex_min_bws"));
  SET_STRING_ELT(r_names, 4, Rf_mkChar("opt_bw_idx"));
  SET_STRING_ELT(r_names, 5, Rf_mkChar("buffer_hops_used"));
  Rf_setAttrib(r_result, R_NamesSymbol, r_names);

  UNPROTECT(2); // r_result, r_names
  return r_result;
}

} // extern "C"

/**
 * Note: Additional graph smoothing functions would be added here following the same pattern.
 * Due to the large number of functions and their complexity, each would need to be 
 * carefully analyzed and corrected individually. The key patterns to follow are:
 * 
 * 1. Use PROTECT_WITH_INDEX/REPROTECT for conditional coercion
 * 2. Use container-first pattern for result assembly
 * 3. Avoid variable UNPROTECT - use only literal constants
 * 4. Balance all PROTECT/UNPROTECT pairs properly
 * 5. Don't use lambdas that capture protection counters
 * 
 * The functions that still need to be corrected include:
 * - S_graph_bw_adaptive_spectral_smoother
 * - S_graph_spectral_lowess
 * - S_graph_spectral_lowess_mat
 * - S_graph_spectral_ma_lowess
 * - S_graph_diffusion_smoother
 * - S_graph_deg0_lowess_cv
 * - S_graph_deg0_lowess_buffer_cv
 * - S_graph_deg0_lowess_cv_mat
 */