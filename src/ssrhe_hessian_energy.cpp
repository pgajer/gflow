#include "ssrhe_hessian_energy_r.h"

#include <Rcpp.h>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct triplet_accumulator_t {
    std::vector<int> i;
    std::vector<int> j;
    std::vector<double> x;

    void add(int row, int col, double val) {
        if (std::isfinite(val) && std::abs(val) > 0.0) {
            i.push_back(row + 1);
            j.push_back(col + 1);
            x.push_back(val);
        }
    }
};

struct cubic_component_t {
    int a;
    int b;
    int c;
};

Eigen::MatrixXd r_matrix_to_eigen(const Rcpp::NumericMatrix& X) {
    Eigen::MatrixXd out(X.nrow(), X.ncol());
    for (int j = 0; j < X.ncol(); ++j) {
        for (int i = 0; i < X.nrow(); ++i) {
            out(i, j) = X(i, j);
        }
    }
    return out;
}

Eigen::MatrixXd pseudo_inverse(
    const Eigen::MatrixXd& M,
    double tol,
    int* rank_out,
    double* condition_out
) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        M, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    const Eigen::VectorXd sing = svd.singularValues();
    const double smax = sing.size() > 0 ? sing.maxCoeff() : 0.0;
    const double cutoff = std::max(tol, 0.0) *
        static_cast<double>(std::max(M.rows(), M.cols())) *
        std::max(smax, 1.0);

    Eigen::VectorXd inv = Eigen::VectorXd::Zero(sing.size());
    int rank = 0;
    double smin_kept = std::numeric_limits<double>::infinity();
    for (int i = 0; i < sing.size(); ++i) {
        if (sing(i) > cutoff) {
            inv(i) = 1.0 / sing(i);
            ++rank;
            smin_kept = std::min(smin_kept, sing(i));
        }
    }
    if (rank_out != nullptr) {
        *rank_out = rank;
    }
    if (condition_out != nullptr) {
        *condition_out = (rank > 0 && std::isfinite(smin_kept) && smin_kept > 0.0) ?
            smax / smin_kept : std::numeric_limits<double>::infinity();
    }

    return svd.matrixV() * inv.asDiagonal() * svd.matrixU().transpose();
}

Rcpp::IntegerMatrix compute_knn_including_self(const Eigen::MatrixXd& X, int k) {
    const int n = X.rows();
    const int d = X.cols();
    Rcpp::IntegerMatrix nn(n, k);

    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double, int> > dist;
        dist.reserve(n);
        for (int j = 0; j < n; ++j) {
            double ss = 0.0;
            for (int a = 0; a < d; ++a) {
                const double diff = X(i, a) - X(j, a);
                ss += diff * diff;
            }
            dist.emplace_back(ss, j);
        }
        std::stable_sort(dist.begin(), dist.end(),
            [](const std::pair<double, int>& lhs,
               const std::pair<double, int>& rhs) {
                if (lhs.first == rhs.first) return lhs.second < rhs.second;
                return lhs.first < rhs.first;
            });
        for (int h = 0; h < k; ++h) {
            nn(i, h) = dist[h].second + 1;
        }
    }

    return nn;
}

Rcpp::List supports_to_list(const std::vector<std::vector<int> >& supports) {
    Rcpp::List out(supports.size());
    for (std::size_t i = 0; i < supports.size(); ++i) {
        Rcpp::IntegerVector ids(supports[i].size());
        for (std::size_t h = 0; h < supports[i].size(); ++h) {
            ids[h] = supports[i][h] + 1;
        }
        out[i] = ids;
    }
    return out;
}

std::vector<std::pair<int, int> > hessian_components(int m) {
    std::vector<std::pair<int, int> > comps;
    comps.reserve(m * (m + 1) / 2);
    for (int a = 0; a < m; ++a) {
        for (int b = a; b < m; ++b) {
            comps.emplace_back(a, b);
        }
    }
    return comps;
}

std::vector<cubic_component_t> third_derivative_components(int m) {
    std::vector<cubic_component_t> comps;
    comps.reserve(m * (m + 1) * (m + 2) / 6);
    for (int a = 0; a < m; ++a) {
        for (int b = a; b < m; ++b) {
            for (int c = b; c < m; ++c) {
                comps.push_back({a, b, c});
            }
        }
    }
    return comps;
}

int third_derivative_component_count(int m) {
    return m * (m + 1) * (m + 2) / 6;
}

int quadratic_component_count(int m) {
    return m * (m + 1) / 2;
}

double third_derivative_scale(const cubic_component_t& comp) {
    if (comp.a == comp.b && comp.b == comp.c) {
        return 6.0;
    }
    if (comp.a == comp.b || comp.a == comp.c || comp.b == comp.c) {
        return 2.0;
    }
    return 1.0;
}

int third_derivative_tensor_multiplicity(const cubic_component_t& comp) {
    if (comp.a == comp.b && comp.b == comp.c) {
        return 1;
    }
    if (comp.a == comp.b || comp.a == comp.c || comp.b == comp.c) {
        return 3;
    }
    return 6;
}

Eigen::MatrixXd build_reduced_design(const Eigen::MatrixXd& Z, int derivative_order) {
    const int k = Z.rows();
    const int m = Z.cols();
    const int q2 = quadratic_component_count(m);
    const int q3 = derivative_order == 3 ? third_derivative_component_count(m) : 0;
    Eigen::MatrixXd design(k, q3 + q2 + m);

    int col = 0;
    if (derivative_order == 3) {
        for (int a = 0; a < m; ++a) {
            for (int b = a; b < m; ++b) {
                for (int c = b; c < m; ++c) {
                    design.col(col) =
                        Z.col(a).array() * Z.col(b).array() * Z.col(c).array();
                    ++col;
                }
            }
        }
    }
    for (int a = 0; a < m; ++a) {
        for (int b = a; b < m; ++b) {
            design.col(col) = Z.col(a).array() * Z.col(b).array();
            ++col;
        }
    }
    for (int a = 0; a < m; ++a) {
        design.col(col) = Z.col(a);
        ++col;
    }
    return design;
}

double compact_kernel_weight(double radius_sq, double d2) {
    if (!(radius_sq > 0.0) || d2 >= radius_sq) {
        return 0.0;
    }
    return std::exp(radius_sq / (d2 - radius_sq));
}

Rcpp::List triplets_to_list(
    const triplet_accumulator_t& trip,
    int nrow,
    int ncol
) {
    return Rcpp::List::create(
        Rcpp::Named("i") = trip.i,
        Rcpp::Named("j") = trip.j,
        Rcpp::Named("x") = trip.x,
        Rcpp::Named("dim") = Rcpp::IntegerVector::create(nrow, ncol)
    );
}

Rcpp::List bs_map_to_triplets(
    const std::map<std::pair<int, int>, double>& values,
    int n
) {
    std::vector<int> ii;
    std::vector<int> jj;
    std::vector<double> xx;
    ii.reserve(values.size());
    jj.reserve(values.size());
    xx.reserve(values.size());

    for (const auto& item : values) {
        if (std::isfinite(item.second) && std::abs(item.second) > 0.0) {
            ii.push_back(item.first.first + 1);
            jj.push_back(item.first.second + 1);
            xx.push_back(item.second);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("i") = ii,
        Rcpp::Named("j") = jj,
        Rcpp::Named("x") = xx,
        Rcpp::Named("dim") = Rcpp::IntegerVector::create(n, n)
    );
}

} // namespace

extern "C" SEXP S_ssrhe_hessian_operator(
    SEXP s_X,
    SEXP s_k,
    SEXP s_tangent_dim,
    SEXP s_nn_index,
    SEXP s_support_index,
    SEXP s_tangent_dim_rule,
    SEXP s_eigen_tolerance,
    SEXP s_derivative_order,
    SEXP s_stabilizer,
    SEXP s_pinv_tol,
    SEXP s_verbose
) {
BEGIN_RCPP
    Rcpp::NumericMatrix Xr(s_X);
    if (Xr.nrow() < 2 || Xr.ncol() < 1) {
        Rcpp::stop("X must be a numeric matrix with at least two rows and one column.");
    }
    const Eigen::MatrixXd X = r_matrix_to_eigen(Xr);
    const int n = X.rows();
    const int ambient_dim = X.cols();
    const int k = Rcpp::as<int>(s_k);
    const int requested_tangent_dim = Rcpp::as<int>(s_tangent_dim);
    const std::string tangent_dim_rule = Rcpp::as<std::string>(s_tangent_dim_rule);
    const double eigen_tolerance = Rcpp::as<double>(s_eigen_tolerance);
    const int derivative_order = Rcpp::as<int>(s_derivative_order);
    const bool stabilizer = Rcpp::as<bool>(s_stabilizer);
    const double pinv_tol = Rcpp::as<double>(s_pinv_tol);
    const bool verbose = Rcpp::as<bool>(s_verbose);

    if (derivative_order != 2 && derivative_order != 3) {
        Rcpp::stop("derivative.order must be 2 or 3.");
    }
    if (derivative_order == 3 && stabilizer) {
        Rcpp::stop("stabilizer is currently only supported for derivative.order = 2.");
    }
    if (!(pinv_tol >= 0.0) || !std::isfinite(pinv_tol)) {
        Rcpp::stop("pinv.tol must be a finite nonnegative scalar.");
    }

    Rcpp::IntegerMatrix nn;
    std::vector<std::vector<int> > supports;
    supports.resize(n);
    const bool has_support_index = !Rf_isNull(s_support_index);

    if (has_support_index) {
        Rcpp::List support_list(s_support_index);
        if (support_list.size() != n) {
            Rcpp::stop("support.index must have length nrow(X).");
        }
        for (int center = 0; center < n; ++center) {
            Rcpp::IntegerVector ids_r = support_list[center];
            if (ids_r.size() < 2) {
                Rcpp::stop("Each support.index element must contain at least two vertices.");
            }
            bool has_center = false;
            supports[center].reserve(ids_r.size());
            std::vector<int> seen(n, 0);
            for (int h = 0; h < ids_r.size(); ++h) {
                const int idx = ids_r[h] - 1;
                if (idx < 0 || idx >= n) {
                    Rcpp::stop("support.index contains a vertex outside 1:nrow(X).");
                }
                if (seen[idx]) {
                    Rcpp::stop("support.index elements must not contain duplicate vertices.");
                }
                seen[idx] = 1;
                supports[center].push_back(idx);
                if (idx == center) {
                    has_center = true;
                }
            }
            if (!has_center) {
                Rcpp::stop("Each support.index element must contain its center vertex.");
            }
        }
    } else {
        if (k < 2 || k > n) {
            Rcpp::stop("k must be between 2 and nrow(X).");
        }
        if (Rf_isNull(s_nn_index)) {
            nn = compute_knn_including_self(X, k);
        } else {
            nn = Rcpp::IntegerMatrix(s_nn_index);
            if (nn.nrow() != n || nn.ncol() != k) {
                Rcpp::stop("nn.index must be an nrow(X) by k integer matrix.");
            }
        }
        for (int center = 0; center < n; ++center) {
            supports[center].reserve(k);
            for (int h = 0; h < k; ++h) {
                supports[center].push_back(nn(center, h) - 1);
            }
        }
    }

    triplet_accumulator_t Atrip;
    std::map<std::pair<int, int>, double> BStrip;

    std::vector<int> row_center;
    std::vector<int> row_component;
    std::vector<int> row_a;
    std::vector<int> row_b;
    std::vector<int> row_c;
    std::vector<int> row_diagonal;
    std::vector<int> row_derivative_order;
    std::vector<int> row_tensor_multiplicity;
    std::vector<double> row_derivative_scale;
    std::vector<double> row_scale;

    std::vector<int> diag_center;
    std::vector<int> diag_k;
    std::vector<int> diag_tangent_dim;
    std::vector<int> diag_design_ncol;
    std::vector<int> diag_design_rank;
    std::vector<double> diag_condition;
    std::vector<double> diag_local_variance_sum;
    std::vector<double> diag_local_variance_ratio;

    int global_row = 0;

    for (int center = 0; center < n; ++center) {
        std::vector<int> ids = supports[center];
        const int ki = static_cast<int>(ids.size());
        int base_local = -1;
        for (int h = 0; h < ki; ++h) {
            const int idx = ids[h];
            if (idx < 0 || idx >= n) {
                Rcpp::stop("Local support contains a vertex outside 1:nrow(X).");
            }
            if (idx == center && base_local < 0) {
                base_local = h;
            }
        }
        if (base_local < 0) {
            Rcpp::stop("Each local support must contain the center vertex.");
        }

        Eigen::MatrixXd local(ki, ambient_dim);
        for (int h = 0; h < ki; ++h) {
            local.row(h) = X.row(ids[h]);
        }
        Eigen::RowVectorXd local_mean = local.colwise().mean();
        Eigen::MatrixXd centered = local.rowwise() - local_mean;

        Eigen::JacobiSVD<Eigen::MatrixXd> pca(
            centered, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        const Eigen::VectorXd sing = pca.singularValues();
        const int max_dim = std::min(static_cast<int>(sing.size()), ambient_dim);

        double total_var = 0.0;
        for (int a = 0; a < sing.size(); ++a) {
            total_var += sing(a) * sing(a);
        }

        int m = requested_tangent_dim;
        if (tangent_dim_rule == "eigen.cumulative") {
            if (requested_tangent_dim > 0) {
                m = std::min(requested_tangent_dim, max_dim);
            } else {
                m = 1;
                double cum = 0.0;
                for (int a = 0; a < max_dim; ++a) {
                    cum += sing(a) * sing(a);
                    if (total_var <= 0.0 || cum / total_var >= eigen_tolerance) {
                        m = a + 1;
                        break;
                    }
                }
            }
        }
        if (m < 1 || m > max_dim) {
            Rcpp::stop("The selected tangent dimension must be between 1 and min(k, ncol(X)).");
        }

        Eigen::MatrixXd coords = centered * pca.matrixV().leftCols(m);
        const Eigen::RowVectorXd base_coord = coords.row(base_local);
        coords.rowwise() -= base_coord;
        Eigen::MatrixXd design = build_reduced_design(coords, derivative_order);
        const int q = derivative_order == 3 ?
            third_derivative_component_count(m) : quadratic_component_count(m);
        const int p = design.cols();
        if (ki < p) {
            Rcpp::stop("Local support is too small for the selected tangent dimension and derivative order.");
        }

        int rank = 0;
        double condition = 0.0;
        Eigen::MatrixXd pinv = pseudo_inverse(design, pinv_tol, &rank, &condition);
        Eigen::MatrixXd reg = pinv;
        reg.col(base_local) -= pinv.rowwise().sum();

        if (derivative_order == 2) {
            std::vector<std::pair<int, int> > comps = hessian_components(m);
            for (int comp = 0; comp < q; ++comp) {
                const int a = comps[comp].first;
                const int b = comps[comp].second;
                const bool diagonal = (a == b);
                const double scale = diagonal ? std::sqrt(2.0) : 1.0;

                for (int h = 0; h < ki; ++h) {
                    Atrip.add(global_row, ids[h], scale * reg(comp, h));
                }

                row_center.push_back(center + 1);
                row_component.push_back(comp + 1);
                row_a.push_back(a + 1);
                row_b.push_back(b + 1);
                row_c.push_back(NA_INTEGER);
                row_diagonal.push_back(diagonal ? 1 : 0);
                row_derivative_order.push_back(2);
                row_tensor_multiplicity.push_back(diagonal ? 1 : 2);
                row_derivative_scale.push_back(NA_REAL);
                row_scale.push_back(scale);
                ++global_row;
            }
        } else {
            std::vector<cubic_component_t> comps = third_derivative_components(m);
            for (int comp = 0; comp < q; ++comp) {
                const cubic_component_t cc = comps[comp];
                const bool diagonal = (cc.a == cc.b && cc.b == cc.c);
                const double derivative_scale = third_derivative_scale(cc);
                const int tensor_multiplicity = third_derivative_tensor_multiplicity(cc);
                const double scale = derivative_scale *
                    std::sqrt(static_cast<double>(tensor_multiplicity));

                for (int h = 0; h < ki; ++h) {
                    Atrip.add(global_row, ids[h], scale * reg(comp, h));
                }

                row_center.push_back(center + 1);
                row_component.push_back(comp + 1);
                row_a.push_back(cc.a + 1);
                row_b.push_back(cc.b + 1);
                row_c.push_back(cc.c + 1);
                row_diagonal.push_back(diagonal ? 1 : 0);
                row_derivative_order.push_back(3);
                row_tensor_multiplicity.push_back(tensor_multiplicity);
                row_derivative_scale.push_back(derivative_scale);
                row_scale.push_back(scale);
                ++global_row;
            }
        }

        if (stabilizer) {
            Eigen::MatrixXd T = design * pinv;
            T.diagonal().array() -= 1.0;
            T.col(base_local) -= T.rowwise().sum();

            double radius_sq = 0.0;
            for (int h = 0; h < ki; ++h) {
                radius_sq = std::max(radius_sq, coords.row(h).squaredNorm());
            }
            radius_sq += 0.1e-9;

            std::vector<double> weights(ki, 0.0);
            for (int h = 0; h < ki; ++h) {
                weights[h] = compact_kernel_weight(radius_sq, coords.row(h).squaredNorm());
            }

            Eigen::MatrixXd local_bs = Eigen::MatrixXd::Zero(ki, ki);
            for (int h = 0; h < ki; ++h) {
                if (weights[h] != 0.0) {
                    local_bs.noalias() += weights[h] * T.row(h).transpose() * T.row(h);
                }
            }
            for (int a = 0; a < ki; ++a) {
                for (int b = 0; b < ki; ++b) {
                    const double val = local_bs(a, b);
                    if (std::isfinite(val) && std::abs(val) > 0.0) {
                        BStrip[std::make_pair(ids[a], ids[b])] += val;
                    }
                }
            }
        }

        double selected_var = 0.0;
        for (int a = 0; a < m && a < sing.size(); ++a) {
            selected_var += sing(a) * sing(a);
        }
        diag_center.push_back(center + 1);
        diag_k.push_back(ki);
        diag_tangent_dim.push_back(m);
        diag_design_ncol.push_back(p);
        diag_design_rank.push_back(rank);
        diag_condition.push_back(condition);
        diag_local_variance_sum.push_back(total_var);
        diag_local_variance_ratio.push_back(total_var > 0.0 ? selected_var / total_var : 1.0);
    }

    if (verbose) {
        Rcpp::Rcout << "Constructed SSRHE Hessian operator with "
                    << global_row << " rows and " << n << " columns." << std::endl;
    }

    Rcpp::List row_table = Rcpp::List::create(
        Rcpp::Named("row") = Rcpp::seq(1, global_row),
        Rcpp::Named("center") = row_center,
        Rcpp::Named("component") = row_component,
        Rcpp::Named("a") = row_a,
        Rcpp::Named("b") = row_b,
        Rcpp::Named("c") = row_c,
        Rcpp::Named("diagonal") = row_diagonal,
        Rcpp::Named("derivative.order") = row_derivative_order,
        Rcpp::Named("derivative.scale") = row_derivative_scale,
        Rcpp::Named("tensor.multiplicity") = row_tensor_multiplicity,
        Rcpp::Named("scale") = row_scale
    );

    Rcpp::List diagnostics = Rcpp::List::create(
        Rcpp::Named("center") = diag_center,
        Rcpp::Named("k") = diag_k,
        Rcpp::Named("tangent.dim") = diag_tangent_dim,
        Rcpp::Named("design.ncol") = diag_design_ncol,
        Rcpp::Named("design.rank") = diag_design_rank,
        Rcpp::Named("design.condition") = diag_condition,
        Rcpp::Named("local.variance.sum") = diag_local_variance_sum,
        Rcpp::Named("local.variance.ratio") = diag_local_variance_ratio
    );

    Rcpp::List parity = Rcpp::List::create(
        Rcpp::Named("knn.includes.self") = true,
        Rcpp::Named("local.pca.on.neighborhood") = true,
        Rcpp::Named("coordinates.centered.at.base.point") = true,
        Rcpp::Named("design.columns") = derivative_order == 3 ?
            "cubic-first-then-quadratic-then-linear-constant-dropped" :
            "quadratic-first-then-linear-constant-dropped",
        Rcpp::Named("fixed.intercept") = "RegMat = XInv - XInv * IndMat",
        Rcpp::Named("derivative.order") = derivative_order,
        Rcpp::Named("diagonal.hessian.scale") = std::sqrt(2.0),
        Rcpp::Named("offdiagonal.hessian.scale") = 1.0,
        Rcpp::Named("third.derivative.scale") =
            "factorial monomial derivative scale times sqrt(symmetric tensor multiplicity)"
    );

    Rcpp::RObject bs_triplet = stabilizer ?
        Rcpp::RObject(bs_map_to_triplets(BStrip, n)) :
        Rcpp::RObject(R_NilValue);

    return Rcpp::List::create(
        Rcpp::Named("A.triplet") = triplets_to_list(Atrip, global_row, n),
        Rcpp::Named("BS.triplet") = bs_triplet,
        Rcpp::Named("nn.index") = has_support_index ?
            Rcpp::RObject(R_NilValue) : Rcpp::RObject(nn),
        Rcpp::Named("support.index") = supports_to_list(supports),
        Rcpp::Named("row.table") = row_table,
        Rcpp::Named("diagnostics") = diagnostics,
        Rcpp::Named("parity") = parity,
        Rcpp::Named("parameters") = Rcpp::List::create(
            Rcpp::Named("k") = k,
            Rcpp::Named("support.variable") = has_support_index,
            Rcpp::Named("tangent.dim") = requested_tangent_dim,
            Rcpp::Named("tangent.dim.rule") = tangent_dim_rule,
            Rcpp::Named("eigen.tolerance") = eigen_tolerance,
            Rcpp::Named("derivative.order") = derivative_order,
            Rcpp::Named("stabilizer") = stabilizer,
            Rcpp::Named("pinv.tol") = pinv_tol
        )
    );
END_RCPP
}
