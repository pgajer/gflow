#include "metric_graph_lowpass_r.h"
#include "omp_compat.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

using spmat_t = Eigen::SparseMatrix<double>;

struct metric_edge_t {
    int u;
    int v;
    double length;
    double conductance;
};

struct metric_operator_t {
    int n = 0;
    std::vector<metric_edge_t> edges;
    std::vector<double> degree;
    std::vector<double> local_scale;
    spmat_t laplacian;
    int n_nonfinite_conductance = 0;
    int n_clipped_conductance = 0;
    double global_sigma = NA_REAL;
};

std::string scalar_string(SEXP x, const char* name) {
    if (TYPEOF(x) != STRSXP || Rf_length(x) != 1) {
        Rf_error("%s must be a character scalar", name);
    }
    return std::string(CHAR(STRING_ELT(x, 0)));
}

double scalar_real(SEXP x, const char* name) {
    if (TYPEOF(x) != REALSXP || Rf_length(x) != 1) {
        Rf_error("%s must be a numeric scalar", name);
    }
    return REAL(x)[0];
}

int scalar_int(SEXP x, const char* name) {
    if (TYPEOF(x) != INTSXP || Rf_length(x) != 1) {
        Rf_error("%s must be an integer scalar", name);
    }
    int out = INTEGER(x)[0];
    if (out == NA_INTEGER) {
        Rf_error("%s cannot be NA", name);
    }
    return out;
}

bool scalar_logical(SEXP x, const char* name) {
    if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
        Rf_error("%s must be a logical scalar", name);
    }
    int out = LOGICAL(x)[0];
    if (out == NA_LOGICAL) {
        Rf_error("%s cannot be NA", name);
    }
    return out != 0;
}

std::vector<std::vector<int>> parse_adj_list(SEXP s_adj_list) {
    if (TYPEOF(s_adj_list) != VECSXP) {
        Rf_error("adj.list must be a list");
    }
    const int n = Rf_length(s_adj_list);
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
        SEXP elt = VECTOR_ELT(s_adj_list, i);
        if (TYPEOF(elt) != INTSXP) {
            Rf_error("adj.list entries must be integer vectors");
        }
        const int m = Rf_length(elt);
        adj[i].resize(m);
        for (int k = 0; k < m; ++k) {
            int v = INTEGER(elt)[k];
            if (v == NA_INTEGER || v < 0 || v >= n) {
                Rf_error("adj.list contains an out-of-range 0-based index");
            }
            adj[i][k] = v;
        }
    }
    return adj;
}

std::vector<std::vector<double>> parse_weight_list(SEXP s_weight_list, int n) {
    if (TYPEOF(s_weight_list) != VECSXP || Rf_length(s_weight_list) != n) {
        Rf_error("weight.list must be a list with the same length as adj.list");
    }
    std::vector<std::vector<double>> weights(n);
    for (int i = 0; i < n; ++i) {
        SEXP elt = VECTOR_ELT(s_weight_list, i);
        if (TYPEOF(elt) != REALSXP) {
            Rf_error("weight.list entries must be numeric vectors");
        }
        const int m = Rf_length(elt);
        weights[i].resize(m);
        for (int k = 0; k < m; ++k) {
            double w = REAL(elt)[k];
            if (!std::isfinite(w) || w <= 0.0) {
                Rf_error("weight.list entries must be finite positive edge lengths");
            }
            weights[i][k] = w;
        }
    }
    return weights;
}

std::vector<double> sorted_copy(std::vector<double> x) {
    std::sort(x.begin(), x.end());
    return x;
}

double empirical_quantile(std::vector<double> x, double q) {
    if (x.empty()) {
        Rf_error("Cannot compute sigma from an empty edge-length set");
    }
    if (!std::isfinite(q) || q < 0.0 || q > 1.0) {
        Rf_error("conductance.sigma.quantile must be in [0, 1]");
    }
    x = sorted_copy(std::move(x));
    if (x.size() == 1) return x[0];
    const double pos = q * static_cast<double>(x.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);
    return x[lo] * (1.0 - frac) + x[hi] * frac;
}

double vector_median(std::vector<double> x) {
    return empirical_quantile(std::move(x), 0.5);
}

void sort_eigenpairs(Eigen::VectorXd& values, Eigen::MatrixXd& vectors) {
    std::vector<int> idx(values.size());
    for (int i = 0; i < values.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return values[a] < values[b];
    });
    Eigen::VectorXd v2(values.size());
    Eigen::MatrixXd U2(vectors.rows(), vectors.cols());
    for (int j = 0; j < values.size(); ++j) {
        v2[j] = values[idx[j]];
        U2.col(j) = vectors.col(idx[j]);
    }
    values = v2;
    vectors = U2;
}

double sparse_laplacian_shift(const spmat_t& L) {
    double max_diag = 0.0;
    for (int k = 0; k < L.outerSize(); ++k) {
        for (spmat_t::InnerIterator it(L, k); it; ++it) {
            if (it.row() == it.col()) {
                max_diag = std::max(max_diag, std::abs(it.value()));
            }
        }
    }
    return -std::max(max_diag, 1.0) * 1e-8;
}

metric_operator_t build_metric_operator(
    SEXP s_adj_list,
    SEXP s_weight_list,
    const std::string& conductance_rule,
    double epsilon,
    double alpha,
    double sigma,
    const std::string& sigma_rule,
    double sigma_quantile,
    int local_k,
    const std::string& laplacian_type
) {
    if (laplacian_type != "unnormalized" && laplacian_type != "symmetric.normalized") {
        Rf_error("laplacian.type must be 'unnormalized' or 'symmetric.normalized'");
    }
    if (!std::isfinite(epsilon) || epsilon <= 0.0) {
        Rf_error("conductance.epsilon must be finite and positive");
    }
    if (!std::isfinite(alpha) || alpha <= 0.0) {
        Rf_error("conductance.alpha must be finite and positive");
    }
    if (local_k < 1) {
        Rf_error("conductance.local.k must be a positive integer");
    }

    std::vector<std::vector<int>> adj = parse_adj_list(s_adj_list);
    const int n = static_cast<int>(adj.size());
    std::vector<std::vector<double>> weights = parse_weight_list(s_weight_list, n);

    metric_operator_t op;
    op.n = n;
    op.degree.assign(n, 0.0);
    op.local_scale.assign(n, NA_REAL);

    std::vector<double> all_lengths;
    std::vector<std::vector<double>> incident_lengths(n);

    for (int i = 0; i < n; ++i) {
        if (adj[i].size() != weights[i].size()) {
            Rf_error("adj.list and weight.list entry lengths differ");
        }
        for (size_t k = 0; k < adj[i].size(); ++k) {
            const int j = adj[i][k];
            const double ell = weights[i][k];
            if (i < j) {
                op.edges.push_back({i, j, ell, NA_REAL});
                all_lengths.push_back(ell);
            }
            incident_lengths[i].push_back(ell);
        }
    }

    if (op.edges.empty()) {
        Rf_error("Graph must contain at least one undirected edge");
    }

    if (conductance_rule == "exp.length" || conductance_rule == "exp.length.squared") {
        if (std::isnan(sigma)) {
            if (sigma_rule == "edge.quantile") {
                sigma = empirical_quantile(all_lengths, sigma_quantile);
            } else if (sigma_rule == "median") {
                sigma = vector_median(all_lengths);
            } else {
                Rf_error("conductance.sigma.rule = 'local.k' is only valid for self.tuned.gaussian");
            }
        }
        if (!std::isfinite(sigma) || sigma <= 0.0) {
            Rf_error("conductance.sigma must be finite and positive");
        }
        op.global_sigma = sigma;
    }

    const double fallback_scale = vector_median(all_lengths);
    if (conductance_rule == "self.tuned.gaussian") {
        for (int i = 0; i < n; ++i) {
            if (incident_lengths[i].empty()) {
                op.local_scale[i] = fallback_scale;
            } else {
                std::vector<double> vals = sorted_copy(incident_lengths[i]);
                const size_t idx = std::min(static_cast<size_t>(local_k - 1), vals.size() - 1);
                op.local_scale[i] = std::max(vals[idx], epsilon);
            }
        }
    }

    constexpr double min_conductance = 1e-300;
    constexpr double max_conductance = 1e300;

    for (auto& e : op.edges) {
        double c = NA_REAL;
        const double ell = e.length;
        if (conductance_rule == "inverse.length.power") {
            c = std::pow(ell + epsilon, -alpha);
        } else if (conductance_rule == "exp.length") {
            c = std::exp(-ell / sigma);
        } else if (conductance_rule == "exp.length.squared") {
            c = std::exp(-(ell * ell) / (sigma * sigma));
        } else if (conductance_rule == "self.tuned.gaussian") {
            const double denom = std::max(op.local_scale[e.u] * op.local_scale[e.v], epsilon);
            c = std::exp(-(ell * ell) / denom);
        } else {
            Rf_error("Unsupported conductance.rule");
        }

        if (!std::isfinite(c)) {
            op.n_nonfinite_conductance++;
            c = (c > 0.0) ? max_conductance : min_conductance;
        }
        if (c < min_conductance) {
            c = min_conductance;
            op.n_clipped_conductance++;
        } else if (c > max_conductance) {
            c = max_conductance;
            op.n_clipped_conductance++;
        }
        e.conductance = c;
        op.degree[e.u] += c;
        op.degree[e.v] += c;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(4 * op.edges.size() + n);
    if (laplacian_type == "unnormalized") {
        for (const auto& e : op.edges) {
            const double c = e.conductance;
            triplets.emplace_back(e.u, e.v, -c);
            triplets.emplace_back(e.v, e.u, -c);
        }
        for (int i = 0; i < n; ++i) {
            triplets.emplace_back(i, i, op.degree[i]);
        }
    } else {
        for (const auto& e : op.edges) {
            if (op.degree[e.u] > 0.0 && op.degree[e.v] > 0.0) {
                const double c_norm = e.conductance / std::sqrt(op.degree[e.u] * op.degree[e.v]);
                triplets.emplace_back(e.u, e.v, -c_norm);
                triplets.emplace_back(e.v, e.u, -c_norm);
            }
        }
        for (int i = 0; i < n; ++i) {
            if (op.degree[i] > 0.0) {
                triplets.emplace_back(i, i, 1.0);
            }
        }
    }

    op.laplacian.resize(n, n);
    op.laplacian.setFromTriplets(triplets.begin(), triplets.end());
    op.laplacian.makeCompressed();
    return op;
}

SEXP make_real_vector(const std::vector<double>& x) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, x.size()));
    for (size_t i = 0; i < x.size(); ++i) REAL(out)[i] = x[i];
    UNPROTECT(1);
    return out;
}

SEXP make_int_vector(const std::vector<int>& x) {
    SEXP out = PROTECT(Rf_allocVector(INTSXP, x.size()));
    for (size_t i = 0; i < x.size(); ++i) INTEGER(out)[i] = x[i];
    UNPROTECT(1);
    return out;
}

SEXP make_summary(const std::vector<double>& x) {
    const int n_fields = 8;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    const char* nms[] = {"min", "q10", "median", "q90", "max", "mean", "n.nonfinite", "n.clipped"};
    for (int i = 0; i < n_fields; ++i) SET_STRING_ELT(names, i, Rf_mkChar(nms[i]));

    std::vector<double> vals = x;
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(*std::min_element(vals.begin(), vals.end())));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(empirical_quantile(vals, 0.10)));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(empirical_quantile(vals, 0.50)));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(empirical_quantile(vals, 0.90)));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(*std::max_element(vals.begin(), vals.end())));
    double sum = 0.0;
    for (double v : vals) sum += v;
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(sum / static_cast<double>(vals.size())));
    SET_VECTOR_ELT(out, 6, Rf_ScalarInteger(0));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger(0));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

SEXP operator_to_sexp(const metric_operator_t& op,
                      const std::string& conductance_rule,
                      double epsilon,
                      double alpha,
                      const std::string& sigma_rule,
                      double sigma_quantile,
                      int local_k,
                      const std::string& laplacian_type) {
    const int n_fields = 8;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    const char* nms[] = {"graph", "edge.table", "conductance", "degree", "laplacian", "laplacian.type", "parameters", "diagnostics"};
    for (int i = 0; i < n_fields; ++i) SET_STRING_ELT(names, i, Rf_mkChar(nms[i]));

    const int nedges = static_cast<int>(op.edges.size());
    std::vector<int> u(nedges), v(nedges);
    std::vector<double> len(nedges), cond(nedges);
    for (int i = 0; i < nedges; ++i) {
        u[i] = op.edges[i].u + 1;
        v[i] = op.edges[i].v + 1;
        len[i] = op.edges[i].length;
        cond[i] = op.edges[i].conductance;
    }

    SEXP graph = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP graph_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(graph_names, 0, Rf_mkChar("n.vertices"));
    SET_STRING_ELT(graph_names, 1, Rf_mkChar("n.edges"));
    SET_VECTOR_ELT(graph, 0, Rf_ScalarInteger(op.n));
    SET_VECTOR_ELT(graph, 1, Rf_ScalarInteger(nedges));
    Rf_setAttrib(graph, R_NamesSymbol, graph_names);
    SET_VECTOR_ELT(out, 0, graph);
    UNPROTECT(2);

    SEXP edge_table = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP edge_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(edge_names, 0, Rf_mkChar("u"));
    SET_STRING_ELT(edge_names, 1, Rf_mkChar("v"));
    SET_STRING_ELT(edge_names, 2, Rf_mkChar("length"));
    SET_STRING_ELT(edge_names, 3, Rf_mkChar("conductance"));
    SET_VECTOR_ELT(edge_table, 0, make_int_vector(u));
    SET_VECTOR_ELT(edge_table, 1, make_int_vector(v));
    SET_VECTOR_ELT(edge_table, 2, make_real_vector(len));
    SET_VECTOR_ELT(edge_table, 3, make_real_vector(cond));
    Rf_setAttrib(edge_table, R_NamesSymbol, edge_names);
    Rf_setAttrib(edge_table, R_ClassSymbol, Rf_mkString("data.frame"));
    SEXP rn = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -nedges;
    Rf_setAttrib(edge_table, R_RowNamesSymbol, rn);
    SET_VECTOR_ELT(out, 1, edge_table);
    UNPROTECT(3);

    SEXP conductance = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP conductance_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(conductance_names, 0, Rf_mkChar("rule"));
    SET_STRING_ELT(conductance_names, 1, Rf_mkChar("values"));
    SET_STRING_ELT(conductance_names, 2, Rf_mkChar("local.scales"));
    SET_STRING_ELT(conductance_names, 3, Rf_mkChar("summary"));
    SET_VECTOR_ELT(conductance, 0, Rf_mkString(conductance_rule.c_str()));
    SET_VECTOR_ELT(conductance, 1, make_real_vector(cond));
    SET_VECTOR_ELT(conductance, 2, make_real_vector(op.local_scale));
    SEXP cond_summary = PROTECT(make_summary(cond));
    SET_VECTOR_ELT(cond_summary, 6, Rf_ScalarInteger(op.n_nonfinite_conductance));
    SET_VECTOR_ELT(cond_summary, 7, Rf_ScalarInteger(op.n_clipped_conductance));
    SET_VECTOR_ELT(conductance, 3, cond_summary);
    Rf_setAttrib(conductance, R_NamesSymbol, conductance_names);
    SET_VECTOR_ELT(out, 2, conductance);
    UNPROTECT(3);

    SET_VECTOR_ELT(out, 3, make_real_vector(op.degree));

    std::vector<int> li, lj;
    std::vector<double> lx;
    li.reserve(op.laplacian.nonZeros());
    lj.reserve(op.laplacian.nonZeros());
    lx.reserve(op.laplacian.nonZeros());
    for (int col = 0; col < op.laplacian.outerSize(); ++col) {
        for (spmat_t::InnerIterator it(op.laplacian, col); it; ++it) {
            li.push_back(it.row());
            lj.push_back(it.col());
            lx.push_back(it.value());
        }
    }
    SEXP lap = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP lap_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(lap_names, 0, Rf_mkChar("dim"));
    SET_STRING_ELT(lap_names, 1, Rf_mkChar("i"));
    SET_STRING_ELT(lap_names, 2, Rf_mkChar("j"));
    SET_STRING_ELT(lap_names, 3, Rf_mkChar("x"));
    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dims)[0] = op.n;
    INTEGER(dims)[1] = op.n;
    SET_VECTOR_ELT(lap, 0, dims);
    SET_VECTOR_ELT(lap, 1, make_int_vector(li));
    SET_VECTOR_ELT(lap, 2, make_int_vector(lj));
    SET_VECTOR_ELT(lap, 3, make_real_vector(lx));
    Rf_setAttrib(lap, R_NamesSymbol, lap_names);
    SET_VECTOR_ELT(out, 4, lap);
    UNPROTECT(3);

    SET_VECTOR_ELT(out, 5, Rf_mkString(laplacian_type.c_str()));

    SEXP params = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP param_names = PROTECT(Rf_allocVector(STRSXP, 7));
    const char* pn[] = {"conductance.epsilon", "conductance.alpha", "conductance.sigma",
                        "conductance.sigma.rule", "conductance.sigma.quantile",
                        "conductance.local.k", "laplacian.type"};
    for (int i = 0; i < 7; ++i) SET_STRING_ELT(param_names, i, Rf_mkChar(pn[i]));
    SET_VECTOR_ELT(params, 0, Rf_ScalarReal(epsilon));
    SET_VECTOR_ELT(params, 1, Rf_ScalarReal(alpha));
    SET_VECTOR_ELT(params, 2, Rf_ScalarReal(op.global_sigma));
    SET_VECTOR_ELT(params, 3, Rf_mkString(sigma_rule.c_str()));
    SET_VECTOR_ELT(params, 4, Rf_ScalarReal(sigma_quantile));
    SET_VECTOR_ELT(params, 5, Rf_ScalarInteger(local_k));
    SET_VECTOR_ELT(params, 6, Rf_mkString(laplacian_type.c_str()));
    Rf_setAttrib(params, R_NamesSymbol, param_names);
    SET_VECTOR_ELT(out, 6, params);
    UNPROTECT(2);

    SEXP diagnostics = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP diag_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(diag_names, 0, Rf_mkChar("degree.summary"));
    SET_STRING_ELT(diag_names, 1, Rf_mkChar("length.summary"));
    SET_VECTOR_ELT(diagnostics, 0, make_summary(op.degree));
    SET_VECTOR_ELT(diagnostics, 1, make_summary(len));
    Rf_setAttrib(diagnostics, R_NamesSymbol, diag_names);
    SET_VECTOR_ELT(out, 7, diagnostics);
    UNPROTECT(2);

    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

void compute_spectrum(const spmat_t& L,
                      int n_eigenpairs,
                      const std::string& eigen_solver,
                      int dense_eigen_threshold,
                      int dense_fallback_threshold,
                      const std::string& dense_fallback,
                      bool verbose,
                      Eigen::VectorXd& eigenvalues,
                      Eigen::MatrixXd& eigenvectors,
                      std::string& backend) {
    const int n = L.rows();
    if (n < 2) Rf_error("Need at least two vertices for eigendecomposition");
    if (n_eigenpairs < 1) Rf_error("n.eigenpairs must be positive");
    n_eigenpairs = std::min(n_eigenpairs, n);

    const bool force_dense = (eigen_solver == "dense");
    const bool force_sparse = (eigen_solver == "sparse");
    const bool auto_dense = (eigen_solver == "auto" && n <= dense_eigen_threshold);

    auto dense_decomp = [&]() -> bool {
        Eigen::MatrixXd Ldense = Eigen::MatrixXd(L);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ldense);
        if (es.info() != Eigen::Success) return false;
        const int m = std::min(n_eigenpairs, n);
        eigenvalues = es.eigenvalues().head(m);
        eigenvectors = es.eigenvectors().leftCols(m);
        backend = "dense";
        return true;
    };

    if (force_dense || auto_dense) {
        if (!dense_decomp()) Rf_error("Dense eigendecomposition failed");
        return;
    }

    const int nev = std::min(n_eigenpairs, std::max(1, n - 2));
    const int ncv = std::min(std::max(2 * nev + 10, 20), n);
    bool sparse_ok = false;
    if (nev < ncv) {
        try {
            Spectra::SparseSymShiftSolve<double> op(L);
            const double sigma = sparse_laplacian_shift(L);
            Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, nev, ncv, sigma);
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                eigenvalues = eigs.eigenvalues();
                eigenvectors = eigs.eigenvectors();
                sort_eigenpairs(eigenvalues, eigenvectors);
                backend = "sparse.shift";
                sparse_ok = true;
            }
        } catch (...) {
            sparse_ok = false;
        }

        if (!sparse_ok) {
            Spectra::SparseSymMatProd<double> op(L);
            Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
            eigs.init();
            eigs.compute(Spectra::SortRule::SmallestAlge, 1000, 1e-8);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                eigenvalues = eigs.eigenvalues();
                eigenvectors = eigs.eigenvectors();
                sort_eigenpairs(eigenvalues, eigenvectors);
                backend = "sparse";
                sparse_ok = true;
            }
        }
    }

    if (sparse_ok) return;

    const bool fallback_allowed =
        dense_fallback == "always" ||
        (dense_fallback == "auto" && n <= dense_fallback_threshold && !force_sparse);
    if (fallback_allowed) {
        if (verbose) Rprintf("Sparse eigendecomposition failed; trying dense fallback\n");
        if (dense_decomp()) {
            backend = "dense.fallback";
            return;
        }
    }
    Rf_error("Sparse eigendecomposition failed and dense fallback was not used");
}

SEXP spectrum_to_sexp(const Eigen::VectorXd& values,
                      const Eigen::MatrixXd& vectors,
                      const std::string& backend) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("eigenvalues"));
    SET_STRING_ELT(names, 1, Rf_mkChar("eigenvectors"));
    SET_STRING_ELT(names, 2, Rf_mkChar("backend"));

    SEXP vals = PROTECT(Rf_allocVector(REALSXP, values.size()));
    for (int i = 0; i < values.size(); ++i) REAL(vals)[i] = values[i];
    SET_VECTOR_ELT(out, 0, vals);
    UNPROTECT(1);

    SEXP vecs = PROTECT(Rf_allocMatrix(REALSXP, vectors.rows(), vectors.cols()));
    for (int j = 0; j < vectors.cols(); ++j) {
        for (int i = 0; i < vectors.rows(); ++i) {
            REAL(vecs)[i + j * vectors.rows()] = vectors(i, j);
        }
    }
    SET_VECTOR_ELT(out, 1, vecs);
    UNPROTECT(1);

    SET_VECTOR_ELT(out, 2, Rf_mkString(backend.c_str()));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

} // namespace

extern "C" SEXP S_metric_graph_lowpass_operator(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_conductance_rule,
    SEXP s_conductance_epsilon,
    SEXP s_conductance_alpha,
    SEXP s_conductance_sigma,
    SEXP s_sigma_rule,
    SEXP s_sigma_quantile,
    SEXP s_local_k,
    SEXP s_laplacian_type
) {
    const std::string conductance_rule = scalar_string(s_conductance_rule, "conductance.rule");
    const double epsilon = scalar_real(s_conductance_epsilon, "conductance.epsilon");
    const double alpha = scalar_real(s_conductance_alpha, "conductance.alpha");
    const double sigma = scalar_real(s_conductance_sigma, "conductance.sigma");
    const std::string sigma_rule = scalar_string(s_sigma_rule, "conductance.sigma.rule");
    const double sigma_quantile = scalar_real(s_sigma_quantile, "conductance.sigma.quantile");
    const int local_k = scalar_int(s_local_k, "conductance.local.k");
    const std::string laplacian_type = scalar_string(s_laplacian_type, "laplacian.type");

    metric_operator_t op = build_metric_operator(
        s_adj_list, s_weight_list, conductance_rule, epsilon, alpha, sigma,
        sigma_rule, sigma_quantile, local_k, laplacian_type
    );
    return operator_to_sexp(op, conductance_rule, epsilon, alpha, sigma_rule,
                            sigma_quantile, local_k, laplacian_type);
}

extern "C" SEXP S_metric_graph_lowpass_spectrum(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_conductance_rule,
    SEXP s_conductance_epsilon,
    SEXP s_conductance_alpha,
    SEXP s_conductance_sigma,
    SEXP s_sigma_rule,
    SEXP s_sigma_quantile,
    SEXP s_local_k,
    SEXP s_laplacian_type,
    SEXP s_n_eigenpairs,
    SEXP s_eigen_solver,
    SEXP s_dense_eigen_threshold,
    SEXP s_dense_fallback_threshold,
    SEXP s_dense_fallback,
    SEXP s_verbose
) {
    const std::string conductance_rule = scalar_string(s_conductance_rule, "conductance.rule");
    const double epsilon = scalar_real(s_conductance_epsilon, "conductance.epsilon");
    const double alpha = scalar_real(s_conductance_alpha, "conductance.alpha");
    const double sigma = scalar_real(s_conductance_sigma, "conductance.sigma");
    const std::string sigma_rule = scalar_string(s_sigma_rule, "conductance.sigma.rule");
    const double sigma_quantile = scalar_real(s_sigma_quantile, "conductance.sigma.quantile");
    const int local_k = scalar_int(s_local_k, "conductance.local.k");
    const std::string laplacian_type = scalar_string(s_laplacian_type, "laplacian.type");
    const int n_eigenpairs = scalar_int(s_n_eigenpairs, "n.eigenpairs");
    const std::string eigen_solver = scalar_string(s_eigen_solver, "eigen.solver");
    const int dense_eigen_threshold = scalar_int(s_dense_eigen_threshold, "dense.eigen.threshold");
    const int dense_fallback_threshold = scalar_int(s_dense_fallback_threshold, "dense.fallback.threshold");
    const std::string dense_fallback = scalar_string(s_dense_fallback, "dense.fallback");
    const bool verbose = scalar_logical(s_verbose, "verbose");

    unsigned int available_threads = gflow_get_max_threads();
    if (available_threads == 0) available_threads = 1;
    Eigen::setNbThreads(available_threads);

    metric_operator_t op = build_metric_operator(
        s_adj_list, s_weight_list, conductance_rule, epsilon, alpha, sigma,
        sigma_rule, sigma_quantile, local_k, laplacian_type
    );

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    std::string backend;
    compute_spectrum(op.laplacian, n_eigenpairs, eigen_solver,
                     dense_eigen_threshold, dense_fallback_threshold,
                     dense_fallback, verbose, eigenvalues, eigenvectors, backend);

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("operator"));
    SET_STRING_ELT(names, 1, Rf_mkChar("spectral"));
    SET_VECTOR_ELT(out, 0, operator_to_sexp(op, conductance_rule, epsilon, alpha,
                                            sigma_rule, sigma_quantile, local_k,
                                            laplacian_type));
    SET_VECTOR_ELT(out, 1, spectrum_to_sexp(eigenvalues, eigenvectors, backend));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}
