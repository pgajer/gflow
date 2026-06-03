#include <Rcpp.h>

#include <R_ext/Lapack.h>

#include <ANN/ANN.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

enum class klp_kernel_t {
    gaussian,
    tricube,
    epanechnikov,
    triangular
};

klp_kernel_t parse_kernel(const std::string& kernel) {
    if (kernel == "gaussian") return klp_kernel_t::gaussian;
    if (kernel == "tricube") return klp_kernel_t::tricube;
    if (kernel == "epanechnikov") return klp_kernel_t::epanechnikov;
    if (kernel == "triangular") return klp_kernel_t::triangular;
    Rcpp::stop("Unsupported kernel: %s", kernel);
}

int design_ncol(const int degree, const int chart_dim) {
    if (degree == 0) return 1;
    if (degree == 1) return 1 + chart_dim;
    if (degree == 2) return 1 + chart_dim + chart_dim * (chart_dim + 1) / 2;
    Rcpp::stop("Unsupported local polynomial degree: %d", degree);
}

double kernel_weight(const double distance, const double bandwidth,
                     const klp_kernel_t kernel) {
    const double denom = bandwidth + std::sqrt(std::numeric_limits<double>::epsilon());
    const double u = distance / denom;
    if (!std::isfinite(u)) return 0.0;
    switch (kernel) {
    case klp_kernel_t::gaussian:
        return std::exp(-0.5 * u * u);
    case klp_kernel_t::tricube:
        if (u >= 1.0) return 0.0;
        return std::pow(1.0 - u * u * u, 3.0);
    case klp_kernel_t::epanechnikov:
        return std::max(0.0, 1.0 - u * u);
    case klp_kernel_t::triangular:
        return std::max(0.0, 1.0 - u);
    }
    return 0.0;
}

double weighted_mean(const std::vector<double>& y,
                     const std::vector<double>& weights,
                     const int n) {
    double sw = 0.0;
    double swy = 0.0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (std::isfinite(wi) && wi > 0.0 && std::isfinite(yi)) {
            sw += wi;
            swy += wi * yi;
        }
    }
    if (sw <= 0.0 || !std::isfinite(sw)) {
        double sy = 0.0;
        int ny = 0;
        for (int i = 0; i < n; ++i) {
            const double yi = y[static_cast<size_t>(i)];
            if (std::isfinite(yi)) {
                sy += yi;
                ++ny;
            }
        }
        return ny > 0 ? sy / static_cast<double>(ny) : NA_REAL;
    }
    return swy / sw;
}

void fill_features(const NumericMatrix& X,
                   const std::vector<int>& original_index,
                   const std::vector<double>& center,
                   const int local_row,
                   const int degree,
                   std::vector<double>& features) {
    const int p = X.ncol();
    features[0] = 1.0;
    if (degree == 0) return;

        const int obs = original_index[static_cast<size_t>(local_row)];
        for (int j = 0; j < p; ++j) {
        features[static_cast<size_t>(1 + j)] = X(obs, j) - center[static_cast<size_t>(j)];
    }
    if (degree == 1) return;

    int col = 1 + p;
    for (int a = 0; a < p; ++a) {
        const double za = features[static_cast<size_t>(1 + a)];
        for (int b = a; b < p; ++b) {
            const double zb = features[static_cast<size_t>(1 + b)];
            features[static_cast<size_t>(col++)] = za * zb;
        }
    }
}

bool solve_spd(std::vector<double>& a, std::vector<double>& b, const int n) {
    int info = 0;
    int nrhs = 1;
    F77_CALL(dpotrf)("L", &n, a.data(), &n, &info FCONE);
    if (info != 0) return false;
    F77_CALL(dpotrs)("L", &n, &nrhs, a.data(), &n, b.data(), &n, &info FCONE);
    return info == 0 && std::isfinite(b[0]);
}

double fit_intercept(const NumericMatrix& X,
                     const std::vector<int>& original_index,
                     const std::vector<double>& center,
                     const std::vector<double>& y,
                     const std::vector<double>& weights,
                     const int support_size,
                     const int degree) {
    const int p = X.ncol();
    const int q = design_ncol(degree, p);
    int n_ok = 0;
    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (std::isfinite(wi) && wi > 0.0 && std::isfinite(yi)) {
            ++n_ok;
        }
    }
    if (n_ok < q) {
        return weighted_mean(y, weights, support_size);
    }
    if (degree == 0) {
        return weighted_mean(y, weights, support_size);
    }

    std::vector<double> xtwx(static_cast<size_t>(q) * static_cast<size_t>(q), 0.0);
    std::vector<double> xtwy(static_cast<size_t>(q), 0.0);
    std::vector<double> f(static_cast<size_t>(q), 0.0);

    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        fill_features(X, original_index, center, i, degree, f);
        for (int a = 0; a < q; ++a) {
            xtwy[static_cast<size_t>(a)] += wi * f[static_cast<size_t>(a)] * yi;
            for (int b = 0; b <= a; ++b) {
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)] +=
                    wi * f[static_cast<size_t>(a)] * f[static_cast<size_t>(b)];
            }
        }
    }
    for (int a = 0; a < q; ++a) {
        for (int b = 0; b < a; ++b) {
            xtwx[static_cast<size_t>(b) + static_cast<size_t>(q) * static_cast<size_t>(a)] =
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)];
        }
    }
    if (!solve_spd(xtwx, xtwy, q)) {
        return weighted_mean(y, weights, support_size);
    }
    return xtwy[0];
}

class AnnTree {
public:
    AnnTree(const NumericMatrix& X, const std::vector<int>& rows)
        : n_(static_cast<int>(rows.size())), p_(X.ncol()),
          data_(nullptr), tree_(nullptr) {
        if (n_ <= 0 || p_ <= 0) {
            throw std::runtime_error("ANN tree needs positive dimensions");
        }
        data_ = annAllocPts(n_, p_);
        try {
            for (int i = 0; i < n_; ++i) {
                const int row = rows[static_cast<size_t>(i)];
                for (int j = 0; j < p_; ++j) {
                    data_[i][j] = X(row, j);
                }
            }
            tree_ = new ANNkd_tree(data_, n_, p_);
        } catch (...) {
            cleanup();
            throw;
        }
    }

    ~AnnTree() {
        cleanup();
    }

    void search(const std::vector<double>& center, const int k,
                std::vector<ANNidx>& nn_idx,
                std::vector<ANNdist>& nn_dist) const {
        ANNpoint query = annAllocPt(p_);
        try {
            for (int j = 0; j < p_; ++j) query[j] = center[static_cast<size_t>(j)];
            tree_->annkSearch(query, k, nn_idx.data(), nn_dist.data(), 0.0);
        } catch (...) {
            annDeallocPt(query);
            throw;
        }
        annDeallocPt(query);
    }

private:
    int n_;
    int p_;
    ANNpointArray data_;
    ANNkd_tree* tree_;

    void cleanup() {
        if (tree_ != nullptr) {
            delete tree_;
            tree_ = nullptr;
        }
        if (data_ != nullptr) {
            annDeallocPts(data_);
            data_ = nullptr;
        }
        annClose();
    }
};

NumericVector predict_coordinates_cpp(const NumericMatrix& X_train,
                                      const NumericVector& y_train,
                                      const NumericMatrix& X_eval,
                                      const IntegerVector& train_rows,
                                      const int support_size,
                                      const int degree,
                                      const klp_kernel_t kernel) {
    const int n_train = train_rows.size();
    const int n_eval = X_eval.nrow();
    const int k = std::min(support_size, n_train);
    if (k <= 0) Rcpp::stop("No training rows available.");

    std::vector<int> rows(static_cast<size_t>(n_train));
    for (int i = 0; i < n_train; ++i) {
        rows[static_cast<size_t>(i)] = train_rows[i] - 1;
    }
    AnnTree tree(X_train, rows);
    NumericVector out(n_eval);
    std::vector<ANNidx> nn_idx(static_cast<size_t>(k));
    std::vector<ANNdist> nn_dist(static_cast<size_t>(k));
    std::vector<double> local_y(static_cast<size_t>(k));
    std::vector<double> local_w(static_cast<size_t>(k));
    std::vector<int> original_index(static_cast<size_t>(k));
    std::vector<double> center(static_cast<size_t>(X_train.ncol()));

    for (int i = 0; i < n_eval; ++i) {
        for (int j = 0; j < X_train.ncol(); ++j) {
            center[static_cast<size_t>(j)] = X_eval(i, j);
        }
        tree.search(center, k, nn_idx, nn_dist);
        double bandwidth = 0.0;
        for (int j = 0; j < k; ++j) {
            bandwidth = std::max(bandwidth, std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)])));
        }
        if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
        for (int j = 0; j < k; ++j) {
            const int local = nn_idx[static_cast<size_t>(j)];
            const int original = rows[static_cast<size_t>(local)];
            original_index[static_cast<size_t>(j)] = original;
            local_y[static_cast<size_t>(j)] = y_train[original];
            const double d = std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            local_w[static_cast<size_t>(j)] = kernel_weight(d, bandwidth, kernel);
        }
        out[i] = fit_intercept(X_train, original_index, center, local_y,
                               local_w, k, degree);
    }
    return out;
}

} // namespace

//' Kernel local polynomial CV RMSE for ambient coordinates
//'
//' Internal C++ backend for `kernel.local.polynomial.cv()`.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_cv_coordinates(
        const NumericMatrix& X,
        const NumericVector& y,
        const IntegerVector& foldid,
        const IntegerVector& support_size,
        const IntegerVector& degree,
        const CharacterVector& kernel) {
    const int n = X.nrow();
    const int n_cand = support_size.size();
    if (y.size() != n || foldid.size() != n ||
        degree.size() != n_cand || kernel.size() != n_cand) {
        Rcpp::stop("Inconsistent input lengths.");
    }
    NumericVector sse(n_cand, 0.0);
    IntegerVector count(n_cand);

    std::vector<int> folds(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) folds[static_cast<size_t>(i)] = foldid[i];
    std::sort(folds.begin(), folds.end());
    folds.erase(std::unique(folds.begin(), folds.end()), folds.end());

    std::vector<klp_kernel_t> kernels(static_cast<size_t>(n_cand));
    int max_support = 0;
    for (int rr = 0; rr < n_cand; ++rr) {
        kernels[static_cast<size_t>(rr)] =
            parse_kernel(Rcpp::as<std::string>(kernel[rr]));
        max_support = std::max(max_support, support_size[rr]);
    }

    for (int ff = 0; ff < static_cast<int>(folds.size()); ++ff) {
        const int fold = folds[static_cast<size_t>(ff)];
        std::vector<int> train_zero;
        std::vector<int> test_zero;
        for (int i = 0; i < n; ++i) {
            if (foldid[i] == fold) {
                test_zero.push_back(i);
            } else {
                train_zero.push_back(i);
            }
        }
        if (train_zero.empty()) continue;
        const int fold_k = std::min(max_support, static_cast<int>(train_zero.size()));
        AnnTree tree(X, train_zero);
        std::vector<ANNidx> nn_idx(static_cast<size_t>(fold_k));
        std::vector<ANNdist> nn_dist(static_cast<size_t>(fold_k));
        std::vector<double> local_y(static_cast<size_t>(fold_k));
        std::vector<double> local_w(static_cast<size_t>(fold_k));
        std::vector<int> original_index(static_cast<size_t>(fold_k));
        std::vector<double> center(static_cast<size_t>(X.ncol()));

        for (const int target : test_zero) {
            for (int j = 0; j < X.ncol(); ++j) {
                center[static_cast<size_t>(j)] = X(target, j);
            }
            tree.search(center, fold_k, nn_idx, nn_dist);
            std::vector<double> distances(static_cast<size_t>(fold_k));
            for (int j = 0; j < fold_k; ++j) {
                const int local = nn_idx[static_cast<size_t>(j)];
                const int original = train_zero[static_cast<size_t>(local)];
                original_index[static_cast<size_t>(j)] = original;
                local_y[static_cast<size_t>(j)] = y[original];
                distances[static_cast<size_t>(j)] =
                    std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            }
            for (int rr = 0; rr < n_cand; ++rr) {
                const int k = std::min(support_size[rr], fold_k);
                double bandwidth = 0.0;
                for (int j = 0; j < k; ++j) {
                    bandwidth = std::max(bandwidth, distances[static_cast<size_t>(j)]);
                }
                if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
                for (int j = 0; j < k; ++j) {
                    local_w[static_cast<size_t>(j)] = kernel_weight(
                        distances[static_cast<size_t>(j)],
                        bandwidth,
                        kernels[static_cast<size_t>(rr)]
                    );
                }
                const double pred = fit_intercept(X, original_index, center,
                                                  local_y, local_w, k,
                                                  degree[rr]);
                const double err = pred - y[target];
                if (std::isfinite(err)) {
                    sse[rr] += err * err;
                    count[rr] += 1;
                }
            }
        }
    }

    NumericVector rmse(n_cand);
    for (int rr = 0; rr < n_cand; ++rr) {
        rmse[rr] = count[rr] > 0 ? std::sqrt(sse[rr] / count[rr]) : NA_REAL;
    }
    return rmse;
}

//' Kernel local polynomial predictions for ambient coordinates
//'
//' Internal C++ backend for `kernel.local.polynomial.cv()`.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_predict_coordinates(
        const NumericMatrix& X_train,
        const NumericVector& y_train,
        const NumericMatrix& X_eval,
        const int support_size,
        const int degree,
        const std::string& kernel) {
    if (y_train.size() != X_train.nrow()) {
        Rcpp::stop("'y_train' must have length nrow(X_train).");
    }
    if (X_eval.ncol() != X_train.ncol()) {
        Rcpp::stop("'X_eval' must have ncol(X_train) columns.");
    }
    IntegerVector train_rows(X_train.nrow());
    for (int i = 0; i < X_train.nrow(); ++i) train_rows[i] = i + 1;
    return predict_coordinates_cpp(X_train, y_train, X_eval, train_rows,
                                   support_size, degree,
                                   parse_kernel(kernel));
}
