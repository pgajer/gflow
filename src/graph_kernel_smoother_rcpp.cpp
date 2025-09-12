// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;

/**
 * @brief Rcpp wrapper for graph_kernel_smoother core routine.
 *
 * Converts R-level adjacency/weights/signal into C++ containers, calls the
 * core C++ implementation, and returns a named R list.
 *
 * @param adj List of integer vectors. adj[i] are 1-based neighbor indices of vertex i+1.
 * @param w   List of numeric vectors matching adj (or NULL). If NULL, unit weights are assumed.
 * @param y   Numeric vector of length length(adj); the signal values on vertices.
 * @param bandwidth Integer bandwidth parameter (>= 1).
 * @param with_details If true, include additional diagnostic output in the result.
 *
 * @return R list with components:
 *   - fitted     : numeric vector of smoothed values
 *   - bandwidth  : integer bandwidth actually used
 *   - details    : (optional) diagnostic list when with_details = true
 *
 * @throws Rcpp::exception if inputs are malformed.
 *
 * @note Indexing: This wrapper passes neighbor indices through as-is.
 *       If your core expects 0-based indices, uncomment the indicated block to convert.
 * @note Thread-safety: This wrapper performs no shared-state mutation;
 *       concurrency depends solely on the core routine.
 */
// Core C++ API you already have (example signature shown for context):
// std::vector<double> graph_kernel_smoother(
//     const std::vector<std::vector<int>>& adj,
//     const std::vector<std::vector<double>>& w,
//     const std::vector<double>& y,
//     int bandwidth,
//     bool with_details);

// [[Rcpp::export]]
List Rcpp_graph_kernel_smoother(List adj, List w, NumericVector y,
                                int bandwidth, bool with_details = false) {
  const R_xlen_t n = adj.size();
  if (n == 0) {
    stop("`adj` must be a non-empty list.");
  }
  if (y.size() != n) {
    stop("length(y) must match length(adj).");
  }
  if (bandwidth < 1) {
    stop("`bandwidth` must be >= 1.");
  }

  // Convert y
  std::vector<double> y_cpp(y.begin(), y.end());

  // Convert adjacency
  std::vector<std::vector<int>> adj_cpp(static_cast<size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    IntegerVector ai = adj[i];
    adj_cpp[static_cast<size_t>(i)] = std::vector<int>(ai.begin(), ai.end());
  }

  // Optional: convert 1-based -> 0-based (UNCOMMENT if your core expects 0-based)
  // for (auto& nb : adj_cpp) {
  //   for (int& idx : nb) { --idx; } // assumes idx >= 1
  // }

  // Convert weights: if NULL, use unit weights per neighbor
  std::vector<std::vector<double>> w_cpp(static_cast<size_t>(n));
  const bool have_w = !w.isNULL();
  if (have_w) {
    if (w.size() != n) {
      stop("`w` must be NULL or a list of the same length as `adj`.");
    }
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    const size_t deg = static_cast<size_t>(adj_cpp[static_cast<size_t>(i)].size());
    if (have_w) {
      NumericVector wi = w[i];
      if (static_cast<size_t>(wi.size()) != deg) {
        stop("`w[[%d]]` must have the same length as `adj[[%d]]`.", (int)(i+1), (int)(i+1));
      }
      w_cpp[static_cast<size_t>(i)] = std::vector<double>(wi.begin(), wi.end());
    } else {
      w_cpp[static_cast<size_t>(i)].assign(deg, 1.0);
    }
  }

  // ---- Call your core implementation ----
  // Replace this placeholder with your actual call and, if applicable, populate details.
  // auto fitted = graph_kernel_smoother(adj_cpp, w_cpp, y_cpp, bandwidth, with_details);
  std::vector<double> fitted = y_cpp; // placeholder: identity

  // ---- Build result ----
  List out = List::create(
    _["fitted"]    = NumericVector(fitted.begin(), fitted.end()),
    _["bandwidth"] = bandwidth
  );

  if (with_details) {
    // Populate real diagnostics from your core; placeholder is empty list.
    out["details"] = List::create();
  }

  return out;
}
