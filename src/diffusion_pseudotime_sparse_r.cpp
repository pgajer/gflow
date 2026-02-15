#include <R.h>
#include <Rinternals.h>

#include "SEXP_cpp_conversion_utils.hpp"
#include "Eigen_utils.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

enum class weight_mode_t {
    identity = 0,
    inverse = 1,
    exp_neg = 2
};

struct transition_build_t {
    Eigen::SparseMatrix<double> P;
    std::vector<double> strength;
};

int as_scalar_int_checked(SEXP x, const char* name, int min_value = 0) {
    if (TYPEOF(x) != INTSXP || Rf_length(x) != 1) {
        Rf_error("%s must be an integer scalar.", name);
    }
    const int v = INTEGER(x)[0];
    if (v == NA_INTEGER || v < min_value) {
        Rf_error("%s must be >= %d.", name, min_value);
    }
    return v;
}

double as_scalar_real_checked(SEXP x, const char* name, double lower_open = -INFINITY) {
    if (TYPEOF(x) != REALSXP || Rf_length(x) != 1) {
        Rf_error("%s must be a numeric scalar.", name);
    }
    const double v = REAL(x)[0];
    if (!std::isfinite(v) || v <= lower_open) {
        Rf_error("%s must be finite and > %.6g.", name, lower_open);
    }
    return v;
}

std::vector<int> parse_index_vector_0based(SEXP s_idx, int n_vertices, const char* name) {
    if (TYPEOF(s_idx) != INTSXP) {
        Rf_error("%s must be an integer vector (0-based indices).", name);
    }
    const R_xlen_t n = XLENGTH(s_idx);
    if (n < 1) {
        Rf_error("%s must have at least one index.", name);
    }
    std::vector<int> out;
    out.reserve(static_cast<size_t>(n));
    const int* p = INTEGER(s_idx);
    for (R_xlen_t i = 0; i < n; ++i) {
        const int v = p[i];
        if (v == NA_INTEGER || v < 0 || v >= n_vertices) {
            Rf_error("%s contains invalid index %d at position %lld (range: 0..%d).",
                     name, v, static_cast<long long>(i + 1), n_vertices - 1);
        }
        out.push_back(v);
    }
    return out;
}

std::vector<double> parse_root_weights(SEXP s_root_weights, size_t n_roots) {
    std::vector<double> out(n_roots, 1.0);
    if (s_root_weights == R_NilValue) {
        const double inv = 1.0 / static_cast<double>(n_roots);
        std::fill(out.begin(), out.end(), inv);
        return out;
    }
    if (TYPEOF(s_root_weights) != REALSXP) {
        Rf_error("root.weights must be NULL or numeric vector.");
    }
    if (XLENGTH(s_root_weights) != static_cast<R_xlen_t>(n_roots)) {
        Rf_error("root.weights length (%lld) must match number of roots (%lld).",
                 static_cast<long long>(XLENGTH(s_root_weights)),
                 static_cast<long long>(n_roots));
    }
    const double* p = REAL(s_root_weights);
    double s = 0.0;
    for (size_t i = 0; i < n_roots; ++i) {
        if (!std::isfinite(p[i]) || p[i] < 0.0) {
            Rf_error("root.weights must be finite and non-negative.");
        }
        out[i] = p[i];
        s += out[i];
    }
    if (!(s > 0.0)) {
        Rf_error("root.weights must have positive sum.");
    }
    for (double& w : out) {
        w /= s;
    }
    return out;
}

double transform_weight(double raw_w, weight_mode_t mode, double param) {
    if (!std::isfinite(raw_w) || raw_w < 0.0) {
        return 0.0;
    }
    if (mode == weight_mode_t::identity) {
        return raw_w;
    }
    if (mode == weight_mode_t::inverse) {
        const double eps = std::max(param, 1e-12);
        return 1.0 / std::max(raw_w, eps);
    }
    const double tau = std::max(param, 1e-12);
    return std::exp(-raw_w / tau);
}

transition_build_t build_sparse_transition_matrix(
    const std::vector<std::vector<int>>& adj,
    const std::vector<std::vector<double>>& weights,
    weight_mode_t mode,
    double weight_param,
    double lazy
) {
    const int n = static_cast<int>(adj.size());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(n) * 12U);

    std::vector<double> strength(static_cast<size_t>(n), 0.0);

    for (int i = 0; i < n; ++i) {
        const auto& nbr = adj[static_cast<size_t>(i)];
        const bool has_weights = !weights.empty();
        if (has_weights && weights[static_cast<size_t>(i)].size() != nbr.size()) {
            Rf_error("weight.list[[%d]] length must match adj.list[[%d]] length.", i + 1, i + 1);
        }

        std::vector<double> w_row(nbr.size(), 0.0);
        double row_sum = 0.0;
        for (size_t k = 0; k < nbr.size(); ++k) {
            const int j = nbr[k];
            if (j < 0 || j >= n) {
                Rf_error("adj.list contains out-of-range index %d at vertex %d.", j, i + 1);
            }
            const double raw_w = has_weights ? weights[static_cast<size_t>(i)][k] : 1.0;
            const double w = transform_weight(raw_w, mode, weight_param);
            if (std::isfinite(w) && w > 0.0) {
                w_row[k] = w;
                row_sum += w;
            }
        }

        strength[static_cast<size_t>(i)] = row_sum;

        if (row_sum <= 0.0) {
            triplets.emplace_back(i, i, 1.0);
            continue;
        }

        const double inv = 1.0 / row_sum;
        for (size_t k = 0; k < nbr.size(); ++k) {
            if (w_row[k] > 0.0) {
                triplets.emplace_back(i, nbr[k], w_row[k] * inv);
            }
        }
    }

    Eigen::SparseMatrix<double> P(n, n);
    P.setFromTriplets(
        triplets.begin(),
        triplets.end(),
        [](const double a, const double b) { return a + b; }
    );
    P.makeCompressed();

    if (lazy < 1.0) {
        const double stay = 1.0 - lazy;
        P *= lazy;
        for (int i = 0; i < n; ++i) {
            P.coeffRef(i, i) += stay;
        }
        P.prune(0.0);
        P.makeCompressed();
    }

    return transition_build_t{std::move(P), std::move(strength)};
}

SEXP sparse_to_csc_slots(const Eigen::SparseMatrix<double>& M, const std::vector<double>& strength) {
    Eigen::SparseMatrix<double> A = M;
    A.makeCompressed();

    const int nnz = A.nonZeros();
    const int nrow = A.rows();
    const int ncol = A.cols();

    SEXP s_i = PROTECT(Rf_allocVector(INTSXP, nnz));
    SEXP s_p = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(ncol) + 1));
    SEXP s_x = PROTECT(Rf_allocVector(REALSXP, nnz));
    SEXP s_dim = PROTECT(Rf_allocVector(INTSXP, 2));
    SEXP s_strength = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(strength.size())));

    std::copy(A.innerIndexPtr(), A.innerIndexPtr() + nnz, INTEGER(s_i));
    std::copy(A.outerIndexPtr(), A.outerIndexPtr() + ncol + 1, INTEGER(s_p));
    std::copy(A.valuePtr(), A.valuePtr() + nnz, REAL(s_x));

    INTEGER(s_dim)[0] = nrow;
    INTEGER(s_dim)[1] = ncol;
    std::copy(strength.begin(), strength.end(), REAL(s_strength));

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(nms, 0, Rf_mkChar("i"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("p"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("x"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("Dim"));
    SET_STRING_ELT(nms, 4, Rf_mkChar("strength"));
    SET_VECTOR_ELT(out, 0, s_i);
    SET_VECTOR_ELT(out, 1, s_p);
    SET_VECTOR_ELT(out, 2, s_x);
    SET_VECTOR_ELT(out, 3, s_dim);
    SET_VECTOR_ELT(out, 4, s_strength);
    Rf_setAttrib(out, R_NamesSymbol, nms);

    UNPROTECT(7);
    return out;
}

Eigen::VectorXd minmax_normalize(const Eigen::VectorXd& x) {
    if (x.size() < 1) {
        return x;
    }
    const double xmin = x.minCoeff();
    const double xmax = x.maxCoeff();
    if (!std::isfinite(xmin) || !std::isfinite(xmax) || xmax <= xmin) {
        return Eigen::VectorXd::Zero(x.size());
    }
    return (x.array() - xmin) / (xmax - xmin);
}

SEXP make_result_with_transition(
    SEXP s_core,
    const Eigen::SparseMatrix<double>& P,
    const std::vector<double>& strength,
    bool with_transition
) {
    if (!with_transition) {
        return s_core;
    }

    SEXP s_transition = PROTECT(sparse_to_csc_slots(P, strength));
    const int n_core = Rf_length(s_core);

    SEXP out = PROTECT(Rf_allocVector(VECSXP, n_core + 1));
    SEXP nms_old = PROTECT(Rf_getAttrib(s_core, R_NamesSymbol));
    SEXP nms_new = PROTECT(Rf_allocVector(STRSXP, n_core + 1));

    for (int i = 0; i < n_core; ++i) {
        SET_VECTOR_ELT(out, i, VECTOR_ELT(s_core, i));
        SET_STRING_ELT(nms_new, i, STRING_ELT(nms_old, i));
    }
    SET_VECTOR_ELT(out, n_core, s_transition);
    SET_STRING_ELT(nms_new, n_core, Rf_mkChar("transition"));
    Rf_setAttrib(out, R_NamesSymbol, nms_new);

    UNPROTECT(4);
    return out;
}

} // namespace

extern "C" SEXP S_build_sparse_transition(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy
) {
    const auto adj = convert_adj_list_from_R(s_adj_list);
    const auto w = convert_weight_list_from_R(s_weight_list);

    const int mode_i = as_scalar_int_checked(s_weight_mode, "weight.mode", 0);
    if (mode_i < 0 || mode_i > 2) {
        Rf_error("weight.mode must be 0 (identity), 1 (inverse), or 2 (exp_neg).");
    }
    const double weight_param = as_scalar_real_checked(s_weight_param, "weight.param", 0.0);
    const double lazy = as_scalar_real_checked(s_lazy, "lazy", -1e-12);
    if (lazy < 0.0 || lazy > 1.0) {
        Rf_error("lazy must be in [0, 1].");
    }

    const auto tr = build_sparse_transition_matrix(
        adj, w, static_cast<weight_mode_t>(mode_i), weight_param, lazy
    );
    return sparse_to_csc_slots(tr.P, tr.strength);
}

extern "C" SEXP S_compute_diffusion_pseudotime_sparse(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_root_vertices,
    SEXP s_root_weights,
    SEXP s_t_steps,
    SEXP s_n_probes,
    SEXP s_seed,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy,
    SEXP s_normalize,
    SEXP s_return_transition
) {
    const auto adj = convert_adj_list_from_R(s_adj_list);
    const auto w = convert_weight_list_from_R(s_weight_list);
    const int n = static_cast<int>(adj.size());
    if (n < 2) {
        Rf_error("Graph must have at least 2 vertices.");
    }

    const int mode_i = as_scalar_int_checked(s_weight_mode, "weight.mode", 0);
    if (mode_i < 0 || mode_i > 2) {
        Rf_error("weight.mode must be 0 (identity), 1 (inverse), or 2 (exp_neg).");
    }
    const double weight_param = as_scalar_real_checked(s_weight_param, "weight.param", 0.0);
    const double lazy = as_scalar_real_checked(s_lazy, "lazy", -1e-12);
    if (lazy < 0.0 || lazy > 1.0) {
        Rf_error("lazy must be in [0, 1].");
    }

    const int t_steps = as_scalar_int_checked(s_t_steps, "t.steps", 1);
    // Kept for forward-compatible API stability; currently not used.
    (void) as_scalar_int_checked(s_n_probes, "n.probes", 1);
    (void) as_scalar_int_checked(s_seed, "seed", 1);
    const bool normalize = Rf_asLogical(s_normalize) == TRUE;
    const bool return_transition = Rf_asLogical(s_return_transition) == TRUE;

    const auto tr = build_sparse_transition_matrix(
        adj, w, static_cast<weight_mode_t>(mode_i), weight_param, lazy
    );

    const auto roots = parse_index_vector_0based(s_root_vertices, n, "root.vertices");
    const auto root_w = parse_root_weights(s_root_weights, roots.size());

    // Rooted diffusion score: q_t = (P^T)^t q_0
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);
    for (size_t r = 0; r < roots.size(); ++r) {
        q[roots[r]] += root_w[r];
    }
    for (int t = 0; t < t_steps; ++t) {
        q = tr.P.transpose() * q;
    }

    // Distance proxy in diffusion-score space: farther vertices have lower q.
    const double q_max = q.maxCoeff();
    Eigen::VectorXd d = (q_max - q.array()).matrix();

    Eigen::VectorXd pt = normalize ? minmax_normalize(d) : d;

    SEXP s_pt = PROTECT(EigenVectorXd_to_SEXP(pt));
    SEXP s_d = PROTECT(EigenVectorXd_to_SEXP(d));
    SEXP s_q = PROTECT(EigenVectorXd_to_SEXP(q));

    SEXP core = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, Rf_mkChar("pseudotime"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("diffusion.distance"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("root.score"));
    SET_VECTOR_ELT(core, 0, s_pt);
    SET_VECTOR_ELT(core, 1, s_d);
    SET_VECTOR_ELT(core, 2, s_q);
    Rf_setAttrib(core, R_NamesSymbol, nms);

    SEXP out = PROTECT(make_result_with_transition(
        core, tr.P, tr.strength, return_transition
    ));

    UNPROTECT(6);
    return out;
}

extern "C" SEXP S_compute_potential_pseudotime_sparse(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_root_vertices,
    SEXP s_root_weights,
    SEXP s_t_steps,
    SEXP s_potential_eps,
    SEXP s_landmark_vertices,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy,
    SEXP s_normalize,
    SEXP s_return_transition
) {
    const auto adj = convert_adj_list_from_R(s_adj_list);
    const auto w = convert_weight_list_from_R(s_weight_list);
    const int n = static_cast<int>(adj.size());
    if (n < 2) {
        Rf_error("Graph must have at least 2 vertices.");
    }

    const int mode_i = as_scalar_int_checked(s_weight_mode, "weight.mode", 0);
    if (mode_i < 0 || mode_i > 2) {
        Rf_error("weight.mode must be 0 (identity), 1 (inverse), or 2 (exp_neg).");
    }
    const double weight_param = as_scalar_real_checked(s_weight_param, "weight.param", 0.0);
    const double lazy = as_scalar_real_checked(s_lazy, "lazy", -1e-12);
    if (lazy < 0.0 || lazy > 1.0) {
        Rf_error("lazy must be in [0, 1].");
    }

    const int t_steps = as_scalar_int_checked(s_t_steps, "t.steps", 1);
    const double potential_eps = as_scalar_real_checked(s_potential_eps, "potential.eps", 0.0);
    const bool normalize = Rf_asLogical(s_normalize) == TRUE;
    const bool return_transition = Rf_asLogical(s_return_transition) == TRUE;

    const auto tr = build_sparse_transition_matrix(
        adj, w, static_cast<weight_mode_t>(mode_i), weight_param, lazy
    );

    const auto roots = parse_index_vector_0based(s_root_vertices, n, "root.vertices");
    const auto root_w = parse_root_weights(s_root_weights, roots.size());

    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    for (size_t r = 0; r < roots.size(); ++r) {
        v[roots[r]] += root_w[r];
    }
    for (int t = 0; t < t_steps; ++t) {
        v = tr.P.transpose() * v;
    }

    Eigen::VectorXd potential = (-v.array().max(potential_eps).log()).matrix();
    Eigen::VectorXd pt = normalize ? minmax_normalize(potential) : potential;

    SEXP s_pt = PROTECT(EigenVectorXd_to_SEXP(pt));
    SEXP s_score = PROTECT(EigenVectorXd_to_SEXP(v));
    SEXP s_potential = PROTECT(EigenVectorXd_to_SEXP(potential));

    SEXP core;
    SEXP nms;

    if (s_landmark_vertices == R_NilValue || XLENGTH(s_landmark_vertices) == 0) {
        core = PROTECT(Rf_allocVector(VECSXP, 3));
        nms = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(nms, 0, Rf_mkChar("pseudotime"));
        SET_STRING_ELT(nms, 1, Rf_mkChar("root.score"));
        SET_STRING_ELT(nms, 2, Rf_mkChar("root.potential"));
        SET_VECTOR_ELT(core, 0, s_pt);
        SET_VECTOR_ELT(core, 1, s_score);
        SET_VECTOR_ELT(core, 2, s_potential);
        Rf_setAttrib(core, R_NamesSymbol, nms);
    } else {
        const auto landmarks = parse_index_vector_0based(s_landmark_vertices, n, "landmark.vertices");
        const int n_landmarks = static_cast<int>(landmarks.size());
        Eigen::MatrixXd U_landmark(n, n_landmarks);

        for (int col = 0; col < n_landmarks; ++col) {
            Eigen::VectorXd u = Eigen::VectorXd::Zero(n);
            u[landmarks[static_cast<size_t>(col)]] = 1.0;
            for (int t = 0; t < t_steps; ++t) {
                u = tr.P.transpose() * u;
            }
            U_landmark.col(col) = (-u.array().max(potential_eps).log()).matrix();
        }

        SEXP s_landmark_potential = PROTECT(EigenMatrixXd_to_SEXP(U_landmark));

        core = PROTECT(Rf_allocVector(VECSXP, 4));
        nms = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(nms, 0, Rf_mkChar("pseudotime"));
        SET_STRING_ELT(nms, 1, Rf_mkChar("root.score"));
        SET_STRING_ELT(nms, 2, Rf_mkChar("root.potential"));
        SET_STRING_ELT(nms, 3, Rf_mkChar("landmark.potential"));
        SET_VECTOR_ELT(core, 0, s_pt);
        SET_VECTOR_ELT(core, 1, s_score);
        SET_VECTOR_ELT(core, 2, s_potential);
        SET_VECTOR_ELT(core, 3, s_landmark_potential);
        Rf_setAttrib(core, R_NamesSymbol, nms);

        UNPROTECT(1); // s_landmark_potential
    }

    SEXP out = PROTECT(make_result_with_transition(
        core, tr.P, tr.strength, return_transition
    ));

    UNPROTECT(6);
    return out;
}
