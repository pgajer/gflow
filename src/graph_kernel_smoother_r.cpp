#include "graph_kernel_smoother.hpp" // For graph_kernel_smoother_t
#include "set_wgraph.hpp"           // For the set_wgraph_t class
#include "error_utils.h"            // For REPORT_ERROR
#include "kernels.h"                // For kernel functions and initialization
#include "SEXP_cpp_conversion_utils.hpp" // For converting R objects to C++
#include "bandwidth_utils.hpp"      // For get_candidate_bws()
#include "progress_utils.hpp"       // For progress_tracker_t

#include <vector>                   // For std::vector
#include <numeric>                  // For std::iota
#include <execution>                // For std::execution::seq/par
#include <atomic>                   // For std::atomic
#include <chrono>                   // For timing
#include <cmath>                    // For math functions
#include <mutex>                    // For std::mutex
#include <execution>                // For std::execution::par_unseq
#include <atomic>                   // For std::atomic

#include <R.h>                      // For R_FlushConsole, Rprintf
#include <Rinternals.h>             // For R C API functions

extern "C" {
    SEXP S_graph_kernel_smoother(
    SEXP s_adj_list,
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
    SEXP s_verbose);
}

/**
 * @brief  R interface for graph kernel smoother with buffer‐zone cross‐validation
 *
 * @details Converts R inputs into C++ types, invokes set_wgraph_t::graph_kernel_smoother,
 * and returns only the fields in graph_kernel_smoother_t:
 *   - predictions:          Optimal‐bandwidth predictions (length |V|)
 *   - bw_predictions:       Per‐bandwidth predictions (|V|×n_bws matrix)
 *   - bw_mean_abs_errors:   Mean absolute Rf_error for each bandwidth (length n_bws)
 *   - vertex_min_bws:       Per‐vertex minimum bandwidths (length |V|)
 *   - opt_bw_idx:           1‐based index of the chosen bandwidth in R
 *   - buffer_hops_used:     Number of buffer hops actually applied
 *
 * @param s_adj_list                R list of integer vectors for adjacency
 * @param s_weight_list             R list of numeric vectors for edge weights
 * @param s_y                       Numeric vector of responses
 * @param s_min_bw_factor           Numeric scalar: min bandwidth factor
 * @param s_max_bw_factor           Numeric scalar: max bandwidth factor
 * @param s_n_bws                   Integer scalar: number of bandwidths
 * @param s_log_grid                Logical scalar: log‐spaced grid?
 * @param s_vertex_hbhd_min_size    Integer: min neighborhood size per vertex
 * @param s_kernel_type             Integer: kernel type code
 * @param s_dist_normalization_factor Numeric: distance normalization factor
 * @param s_use_uniform_weights     Logical: uniform vs. kernel weights
 * @param s_buffer_hops             Integer: buffer‐zone hop distance
 * @param s_auto_buffer_hops        Logical: auto‐determine buffer hops?
 * @param s_n_folds                 Integer: number of CV folds
 * @param s_with_bw_predictions     Logical: store per‐bandwidth predictions?
 * @param s_precision               Numeric: precision for bandwidth grid
 * @param s_verbose                 Logical: verbose output?
 *
 * @return SEXP R list with components:
 *         "predictions", "bw_predictions", "bw_mean_abs_errors",
 *         "vertex_min_bws", "opt_bw_idx", "buffer_hops_used".
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

  // ---- y: coerce defensively, copy out; long-vector safe ----
  std::vector<double> y;
  {
    int tprot = 0;
    SEXP sy = s_y;
    if (TYPEOF(sy) != REALSXP) { sy = PROTECT(Rf_coerceVector(sy, REALSXP)); ++tprot; }
    const R_xlen_t ny = XLENGTH(sy);
    y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));
    if (tprot) UNPROTECT(tprot);
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

  // ---- Build result (container-first; per-element PROTECT/UNPROTECT) ----
  int nprot = 0;
  const int N_ELT = 6;
  SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N_ELT)); ++nprot;

  // names while r_result is protected
  {
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, N_ELT)); ++nprot;
    SET_STRING_ELT(r_names, 0, Rf_mkChar("predictions"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("bw_predictions"));
    SET_STRING_ELT(r_names, 2, Rf_mkChar("bw_mean_abs_errors"));
    SET_STRING_ELT(r_names, 3, Rf_mkChar("vertex_min_bws"));
    SET_STRING_ELT(r_names, 4, Rf_mkChar("opt_bw_idx"));
    SET_STRING_ELT(r_names, 5, Rf_mkChar("buffer_hops_used"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    // keep names protected only until attribute is set
    UNPROTECT(1); --nprot;
  }

  // 0: predictions (numeric)
  {
    SEXP s = PROTECT(convert_vector_double_to_R(result.predictions));
    if (static_cast<R_xlen_t>(result.predictions.size()) != N) {
      UNPROTECT(1); // s
      UNPROTECT(nprot);
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
      UNPROTECT(nprot);
      Rf_error("bw_predictions: row count (%zu) must equal length(y) (%lld).",
               nrow_sz, static_cast<long long>(N));
    }
    for (size_t j = 1; j < ncol_sz; ++j) {
      if (result.bw_predictions[j].size() != nrow_sz) {
        UNPROTECT(nprot);
        Rf_error("bw_predictions: column %zu has length %zu; expected %zu.",
                 j + 1, result.bw_predictions[j].size(), nrow_sz);
      }
    }
    if (nrow_sz > static_cast<size_t>(R_XLEN_T_MAX) ||
        ncol_sz > static_cast<size_t>(R_XLEN_T_MAX)) {
      UNPROTECT(nprot);
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
      UNPROTECT(nprot);
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
      UNPROTECT(nprot);
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

  UNPROTECT(nprot); // r_result
  return r_result;
}
