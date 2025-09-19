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
Rcpp::List Rcpp_graph_kernel_smoother(Rcpp::List adj, Rcpp::List w, Rcpp::NumericVector y,
                                      int bandwidth, bool with_details = false) {
  const int n = (int) adj.size();
  if (n == 0) Rcpp::stop("`adj` must be a non-empty list.");
  if ((int)y.size() != n) Rcpp::stop("length(y) must match length(adj).");
  if (bandwidth < 1) Rcpp::stop("`bandwidth` must be >= 1.");

  // Convert y
  std::vector<double> y_cpp(y.begin(), y.end());

  // Convert adjacency
  std::vector<std::vector<int>> adj_cpp((size_t)n);
  for (int i = 0; i < n; ++i) {
    Rcpp::IntegerVector ai = adj[i];
    adj_cpp[(size_t)i] = std::vector<int>(ai.begin(), ai.end());
  }

  // Convert weights (NULL => unit weights)
  std::vector<std::vector<double>> w_cpp((size_t)n);
  const bool have_w = !w.isNULL();
  if (have_w && (int)w.size() != n) {
    Rcpp::stop("`w` must be NULL or a list of the same length as `adj`.");
  }
  for (int i = 0; i < n; ++i) {
    const size_t deg = adj_cpp[(size_t)i].size();
    if (have_w) {
      Rcpp::NumericVector wi = w[i];
      if ((size_t)wi.size() != deg) {
        Rcpp::stop("`w[[%d]]` must have the same length as `adj[[%d]]`.", i+1, i+1);
      }
      w_cpp[(size_t)i] = std::vector<double>(wi.begin(), wi.end());
    } else {
      w_cpp[(size_t)i].assign(deg, 1.0);
    }
  }

  // ---- Call core implementation ----
  // std::vector<double> fitted = graph_kernel_smoother(adj_cpp, w_cpp, y_cpp, bandwidth, with_details);
  std::vector<double> fitted = y_cpp; // placeholder

  // ---- Build result without Named push-back ----
  const int out_len = with_details ? 3 : 2;
  Rcpp::List out(out_len);
  Rcpp::CharacterVector nm(out_len);

  out[0] = Rcpp::NumericVector(fitted.begin(), fitted.end());
  nm[0]  = "fitted";

  out[1] = bandwidth;
  nm[1]  = "bandwidth";

  if (with_details) {
    out[2] = Rcpp::List();  // fill with real diagnostics if available
    nm[2]  = "details";
  }

  out.attr("names") = nm;
  return out;
}
