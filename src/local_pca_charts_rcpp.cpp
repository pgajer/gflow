#include "local_pca_charts.hpp"

#include <Rcpp.h>

#include <limits>

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::Nullable;

namespace {

Eigen::MatrixXd numeric_matrix_to_eigen(const NumericMatrix& X) {
    Eigen::MatrixXd out(X.nrow(), X.ncol());
    for (int j = 0; j < X.ncol(); ++j) {
        for (int i = 0; i < X.nrow(); ++i) {
            out(i, j) = X(i, j);
        }
    }
    return out;
}

NumericMatrix eigen_matrix_to_numeric(const Eigen::MatrixXd& X) {
    NumericMatrix out(X.rows(), X.cols());
    for (int j = 0; j < X.cols(); ++j) {
        for (int i = 0; i < X.rows(); ++i) {
            out(i, j) = X(i, j);
        }
    }
    return out;
}

NumericVector eigen_vector_to_numeric(const Eigen::VectorXd& x) {
    NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i) out[i] = x(i);
    return out;
}

} // namespace

//' Local PCA chart
//'
//' Internal C++ backend for shared local-PCA chart construction.
//'
//' @keywords internal
// [[Rcpp::export]]
List rcpp_local_pca_chart(
        const NumericMatrix& X_support,
        const NumericVector& center,
        const int chart_dim,
        const std::string& center_mode = "anchor",
        const std::string& dim_rule = "fixed",
        const double eigen_tolerance = 0.9,
        Nullable<NumericVector> weights = R_NilValue,
        const bool rebase_to_anchor = true,
        const bool orient_basis = false) {
    if (center.size() != X_support.ncol()) {
        Rcpp::stop("'center' must have ncol(X_support) entries.");
    }
    Eigen::MatrixXd local = numeric_matrix_to_eigen(X_support);
    Eigen::RowVectorXd anchor(center.size());
    for (int j = 0; j < center.size(); ++j) anchor(j) = center[j];

    Eigen::VectorXd weight_vec;
    Eigen::VectorXd* weight_ptr = nullptr;
    if (weights.isNotNull()) {
        NumericVector w(weights);
        if (w.size() != X_support.nrow()) {
            Rcpp::stop("'weights' must have nrow(X_support) entries.");
        }
        weight_vec.resize(w.size());
        for (int i = 0; i < w.size(); ++i) weight_vec(i) = w[i];
        weight_ptr = &weight_vec;
    }

    gflow::local_pca_chart_result_t chart =
        gflow::compute_local_pca_chart(
            local,
            anchor,
            chart_dim,
            dim_rule,
            eigen_tolerance,
            center_mode,
            weight_ptr,
            rebase_to_anchor,
            orient_basis
        );

    return List::create(
        Rcpp::Named("coordinates") =
            eigen_matrix_to_numeric(chart.coordinates),
        Rcpp::Named("basis") = eigen_matrix_to_numeric(chart.basis),
        Rcpp::Named("singular.values") =
            eigen_vector_to_numeric(chart.singular_values),
        Rcpp::Named("chart.dim") = chart.selected_dim,
        Rcpp::Named("total.variance") = chart.total_variance,
        Rcpp::Named("selected.variance.ratio") =
            chart.selected_variance_ratio
    );
}
