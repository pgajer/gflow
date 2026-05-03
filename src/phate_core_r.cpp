#include "phate_core_r.h"

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

int scalar_int(SEXP x, const char* name, int min_value) {
    if (TYPEOF(x) != INTSXP || Rf_length(x) != 1) {
        Rf_error("%s must be an integer scalar.", name);
    }
    const int v = INTEGER(x)[0];
    if (v == NA_INTEGER || v < min_value) {
        Rf_error("%s must be >= %d.", name, min_value);
    }
    return v;
}

double scalar_real(SEXP x, const char* name, double lower_open) {
    if (TYPEOF(x) != REALSXP || Rf_length(x) != 1) {
        Rf_error("%s must be a numeric scalar.", name);
    }
    const double v = REAL(x)[0];
    if (!std::isfinite(v) || v <= lower_open) {
        Rf_error("%s must be finite and > %.6g.", name, lower_open);
    }
    return v;
}

bool scalar_bool(SEXP x, const char* name) {
    if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
        Rf_error("%s must be TRUE/FALSE.", name);
    }
    const int v = LOGICAL(x)[0];
    if (v == NA_LOGICAL) {
        Rf_error("%s must be TRUE/FALSE.", name);
    }
    return v != 0;
}

int nrow_matrix(SEXP x, const char* name) {
    if (!Rf_isMatrix(x) || TYPEOF(x) != REALSXP) {
        Rf_error("%s must be a numeric matrix.", name);
    }
    SEXP dim = PROTECT(Rf_getAttrib(x, R_DimSymbol));
    if (dim == R_NilValue || TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
        UNPROTECT(1);
        Rf_error("%s must have valid matrix dimensions.", name);
    }
    const int nr = INTEGER(dim)[0];
    const int nc = INTEGER(dim)[1];
    UNPROTECT(1);
    if (nr != nc) {
        Rf_error("%s must be square.", name);
    }
    if (nr < 2) {
        Rf_error("%s must have at least 2 rows.", name);
    }
    return nr;
}

double mat_get(const double* x, int n, int i, int j) {
    return x[i + n * j];
}

void mat_set(double* x, int n, int i, int j, double v) {
    x[i + n * j] = v;
}

} // namespace

extern "C" SEXP S_phate_build_kernel(SEXP s_D,
                                     SEXP s_k,
                                     SEXP s_decay,
                                     SEXP s_thresh,
                                     SEXP s_bandwidth_scale,
                                     SEXP s_symm_code,
                                     SEXP s_diag_one) {
    const int n = nrow_matrix(s_D, "D");
    const int k_in = scalar_int(s_k, "k", 1);
    const double decay = scalar_real(s_decay, "decay", 0.0);
    const double thresh = scalar_real(s_thresh, "thresh", -std::numeric_limits<double>::infinity());
    const double bandwidth_scale = scalar_real(s_bandwidth_scale, "bandwidth.scale", 0.0);
    const int symm_code = scalar_int(s_symm_code, "symm.code", 0);
    const bool diag_one = scalar_bool(s_diag_one, "diag.one");

    if (symm_code > 2) {
        Rf_error("symm.code must be 0 (mean), 1 (max), or 2 (none).");
    }
    if (thresh < 0.0 || thresh >= 1.0) {
        Rf_error("thresh must be in [0, 1).");
    }

    const double eps = std::numeric_limits<double>::epsilon();
    const int k = std::min(k_in, n - 1);
    const double* D = REAL(s_D);

    SEXP s_K_dir = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP s_K = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP s_P = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP s_bandwidth = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP s_row_sum = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP s_n_retained = PROTECT(Rf_allocVector(INTSXP, n));

    double* K_dir = REAL(s_K_dir);
    double* K = REAL(s_K);
    double* P = REAL(s_P);
    double* bandwidth = REAL(s_bandwidth);
    double* row_sum = REAL(s_row_sum);
    int* n_retained = INTEGER(s_n_retained);

    std::fill(K_dir, K_dir + static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    std::fill(K, K + static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    std::fill(P, P + static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);

    for (int i = 0; i < n; ++i) {
        std::vector<double> d_row;
        d_row.reserve(static_cast<size_t>(n - 1));
        for (int j = 0; j < n; ++j) {
            const double dij_raw = mat_get(D, n, i, j);
            if (!std::isfinite(dij_raw) || dij_raw < 0.0) {
                Rf_error("D contains non-finite or negative values.");
            }
            if (i != j) {
                d_row.push_back(dij_raw);
            }
        }
        std::sort(d_row.begin(), d_row.end());
        const double bw = std::max(d_row[static_cast<size_t>(k - 1)] * bandwidth_scale, eps);
        bandwidth[i] = bw;

        int retained = 0;
        for (int j = 0; j < n; ++j) {
            double val = 0.0;
            if (i == j) {
                val = diag_one ? 1.0 : 0.0;
            } else {
                const double dij = mat_get(D, n, i, j);
                const double scaled = dij / bw;
                val = std::exp(-std::pow(scaled, decay));
                if (!std::isfinite(val) || val < thresh) {
                    val = 0.0;
                }
            }
            if (val > 0.0) {
                ++retained;
            }
            mat_set(K_dir, n, i, j, val);
        }
        n_retained[i] = retained;
    }

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const double a = mat_get(K_dir, n, i, j);
            const double b = mat_get(K_dir, n, j, i);
            double v = a;
            if (symm_code == 0) {
                v = 0.5 * (a + b);
            } else if (symm_code == 1) {
                v = std::max(a, b);
            }
            mat_set(K, n, i, j, v);
        }
    }

    for (int i = 0; i < n; ++i) {
        double rs = 0.0;
        for (int j = 0; j < n; ++j) {
            const double v = mat_get(K, n, i, j);
            if (!std::isfinite(v) || v < 0.0) {
                Rf_error("Internal kernel construction produced invalid affinity.");
            }
            rs += v;
        }
        row_sum[i] = rs;
        if (rs > 0.0) {
            for (int j = 0; j < n; ++j) {
                mat_set(P, n, i, j, mat_get(K, n, i, j) / rs);
            }
        } else {
            mat_set(P, n, i, i, 1.0);
        }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, Rf_mkChar("K.directed"));
    SET_STRING_ELT(names, 1, Rf_mkChar("K"));
    SET_STRING_ELT(names, 2, Rf_mkChar("P"));
    SET_STRING_ELT(names, 3, Rf_mkChar("bandwidth"));
    SET_STRING_ELT(names, 4, Rf_mkChar("row_sum"));
    SET_STRING_ELT(names, 5, Rf_mkChar("n_retained_directed"));
    SET_VECTOR_ELT(out, 0, s_K_dir);
    SET_VECTOR_ELT(out, 1, s_K);
    SET_VECTOR_ELT(out, 2, s_P);
    SET_VECTOR_ELT(out, 3, s_bandwidth);
    SET_VECTOR_ELT(out, 4, s_row_sum);
    SET_VECTOR_ELT(out, 5, s_n_retained);
    Rf_setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(8);
    return out;
}
