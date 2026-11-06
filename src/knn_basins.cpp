#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// Exact neighbors in the supplied coordinates. Keep only one distance row
// while searching; a dense n-by-n distance matrix is never constructed.
// [[Rcpp::export(rng = false)]]
Rcpp::List rcpp_knn_basins_neighbors(const Rcpp::NumericMatrix& X,
                                    const int k,
                                    const bool linfty) {
    const int n = X.nrow(), p = X.ncol();
    if (n < 2 || p < 1 || k < 1 || k >= n) {
        Rcpp::stop("Require at least two rows, one column, and 1 <= k < nrow(X).");
    }
    // Row-major storage makes the exhaustive distance scan contiguous.
    std::vector<double> x(static_cast<size_t>(n) * p);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < p; ++c) {
            if (!std::isfinite(X(i, c))) Rcpp::stop("X must be finite.");
            x[static_cast<size_t>(i) * p + c] = X(i, c);
        }
    }
    using Neighbor = std::pair<double, int>;
    std::vector<Neighbor> candidates;
    candidates.reserve(n - 1);
    std::vector<std::vector<std::pair<int, double>>> edges(n);
    Rcpp::IntegerMatrix index(n, k);
    Rcpp::NumericMatrix distance(n, k);
    for (int i = 0; i < n; ++i) {
        Rcpp::checkUserInterrupt();
        candidates.clear();
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (j % 1024 == 0) Rcpp::checkUserInterrupt();
            double scale = 0.0, sumsq = 1.0;
            for (int c = 0; c < p; ++c) {
                const double delta = std::abs(x[static_cast<size_t>(i) * p + c] -
                                              x[static_cast<size_t>(j) * p + c]);
                if (!std::isfinite(delta)) {
                    Rcpp::stop("Pairwise distances overflow; rescale X explicitly.");
                }
                if (linfty) {
                    scale = std::max(scale, delta);
                } else if (delta > 0.0) {
                    // Scaled sum of squares avoids overflow and underflow
                    // when the Euclidean distance itself is representable.
                    if (delta > scale) {
                        const double ratio = scale / delta;
                        sumsq = 1.0 + sumsq * ratio * ratio;
                        scale = delta;
                    } else {
                        const double ratio = delta / scale;
                        sumsq += ratio * ratio;
                    }
                }
            }
            const double d = linfty ? scale : scale * std::sqrt(sumsq);
            if (!std::isfinite(d)) {
                Rcpp::stop("Pairwise distances overflow; rescale X explicitly.");
            }
            candidates.emplace_back(d, j);
        }
        // Distance ties use row index, including ties at the kth neighbor.
        std::partial_sort(candidates.begin(), candidates.begin() + k,
                          candidates.end());
        for (int rank = 0; rank < k; ++rank) {
            const int j = candidates[rank].second;
            const double d = candidates[rank].first;
            index(i, rank) = j + 1;
            distance(i, rank) = d;
            edges[i].emplace_back(j, d);
            edges[j].emplace_back(i, d);
        }
    }
    Rcpp::List adjacency(n), lengths(n);
    for (int i = 0; i < n; ++i) {
        auto& row = edges[i];
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end(),
                              [](const auto& a, const auto& b) {
                                  return a.first == b.first;
                              }), row.end());
        Rcpp::IntegerVector neighbors(row.size());
        Rcpp::NumericVector weights(row.size());
        for (size_t j = 0; j < row.size(); ++j) {
            neighbors[j] = row[j].first + 1;
            weights[j] = row[j].second;
        }
        adjacency[i] = neighbors;
        lengths[i] = weights;
    }
    Rcpp::IntegerVector component(n);
    int count = 0;
    std::vector<int> queue;
    for (int i = 0; i < n; ++i) {
        if (component[i] != 0) continue;
        Rcpp::checkUserInterrupt();
        ++count;
        queue.clear();
        queue.push_back(i);
        component[i] = count;
        for (size_t head = 0; head < queue.size(); ++head) {
            for (const auto& edge : edges[queue[head]]) {
                const int j = edge.first;
                if (component[j] == 0) {
                    component[j] = count;
                    queue.push_back(j);
                }
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::_["knn.index"] = index,
        Rcpp::_["knn.distance"] = distance,
        Rcpp::_["adj.list"] = adjacency,
        Rcpp::_["edge.length.list"] = lengths,
        Rcpp::_["component"] = component,
        Rcpp::_["n.components"] = count
    );
}
