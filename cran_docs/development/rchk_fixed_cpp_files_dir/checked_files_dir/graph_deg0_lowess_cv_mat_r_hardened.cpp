/**
 * @file graph_deg0_lowess_cv_mat_r_hardened.cpp
 * @brief Hardened rchk-safe R entrypoint for degree-0 LOWESS CV over multiple columns.
 *
 * This is a drop-in replacement for S_graph_deg0_lowess_cv_mat, written to follow
 * rchk-safe patterns:
 *  - Container-first PROTECT
 *  - Only fixed-count UNPROTECTs
 *  - No lambdas that manipulate protection state
 *  - Defensive coercion of inputs via PROTECT_WITH_INDEX/REPROTECT
 *  - Long-vector safety using XLENGTH/R_xlen_t
 *  - Robust scalar extraction via Rf_asInteger/Rf_asReal/Rf_asLogical
 *
 * Assumptions about core compute layer (kept as externs):
 *   - convert_adj_list_from_R / convert_weight_list_from_R return plain C++ containers.
 *   - graph_deg0_lowess_cv_mat(...) performs no R allocations.
 *   - The result struct matches the fields used below.
 */

#include <vector>
#include <algorithm>

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

// Helpers (provided elsewhere in the project)
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Core compute result & function (provided elsewhere)
struct graph_deg0_lowess_cv_mat_result_t {
    // Column-major CV error matrix of size (n_bws x n_columns)
    std::vector<double> cv_errors;     // length = n_bws * n_cols
    int n_bws;
    int n_cols;

    // Per-column mean CV errors across folds (length n_bws or n_cols depending on convention).
    // Here we assume mean across folds for each (bw, column) has already been taken, so
    // mean_cv_errors is a vector of length n_bws for each column, flattened column-major:
    // length = n_bws * n_cols.
    std::vector<double> mean_cv_errors;

    // For each column, the 0-based index of the optimal bandwidth on the grid
    std::vector<int> opt_bw_idxs;      // length = n_cols

    // The actual bandwidth factors used on the grid (length n_bws)
    std::vector<double> bw_factors;
};

graph_deg0_lowess_cv_mat_result_t graph_deg0_lowess_cv_mat(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<std::vector<double>>& Y_cols, // Y as list of column vectors
    int n_bws,
    int log_grid, // 0/1
    double min_bw_factor,
    double max_bw_factor,
    double dist_normalization_factor,
    int ikernel,
    double precision,
    int n_folds,
    int buffer_hops,
    int auto_buffer_hops,
    int use_uniform_weights,
    int with_bw_predictions, // not returned here; compute layer may ignore
    int verbose);

// Hardened R entrypoint
extern "C" SEXP S_graph_deg0_lowess_cv_mat(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,                         // matrix (n x p) or list of length p with numeric vectors length n
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    SEXP s_precision,
    SEXP s_n_folds,
    SEXP s_buffer_hops,
    SEXP s_auto_buffer_hops,
    SEXP s_use_uniform_weights,
    SEXP s_with_bw_predictions,       // accepted for signature parity; not returned here
    SEXP s_verbose)
{
    // Convert graph inputs
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Extract Y (either matrix or list) into vector< vector<double> > of columns.
    std::vector<std::vector<double>> Y_cols;

    if (Rf_isMatrix(s_Y)) {
        // Coerce to REAL and copy column-major
        SEXP s_Y_mat = s_Y;
        PROTECT_INDEX pY;
        PROTECT_WITH_INDEX(s_Y_mat, &pY);
        if (TYPEOF(s_Y_mat) != REALSXP) REPROTECT(s_Y_mat = Rf_coerceVector(s_Y_mat, REALSXP), pY);

        SEXP s_dim = Rf_getAttrib(s_Y, R_DimSymbol);
        if (TYPEOF(s_dim) != INTSXP || LENGTH(s_dim) != 2) {
            UNPROTECT(1); // s_Y_mat
            Rf_error("Y must be a numeric matrix with dim attribute");
        }
        const int n_rows = INTEGER(s_dim)[0];
        const int n_cols = INTEGER(s_dim)[1];

        Y_cols.resize(static_cast<size_t>(n_cols));
        const double* Yp = REAL(s_Y_mat);
        for (int j = 0; j < n_cols; ++j) {
            Y_cols[j].resize(static_cast<size_t>(n_rows));
            const double* col = Yp + static_cast<R_xlen_t>(j) * static_cast<R_xlen_t>(n_rows);
            for (int i = 0; i < n_rows; ++i) {
                Y_cols[j][static_cast<size_t>(i)] = col[i];
            }
        }
        UNPROTECT(1); // s_Y_mat
    } else if (Rf_isNewList(s_Y)) {
        // List of numeric vectors; coerce each element to REAL and copy
        const R_xlen_t p = XLENGTH(s_Y);
        Y_cols.resize(static_cast<size_t>(p));
        for (R_xlen_t j = 0; j < p; ++j) {
            SEXP s_y_j = VECTOR_ELT(s_Y, j);
            PROTECT_INDEX py;
            PROTECT_WITH_INDEX(s_y_j, &py);
            if (TYPEOF(s_y_j) != REALSXP) REPROTECT(s_y_j = Rf_coerceVector(s_y_j, REALSXP), py);
            const R_xlen_t n = XLENGTH(s_y_j);
            Y_cols[static_cast<size_t>(j)].resize(static_cast<size_t>(n));
            const double* yp = REAL(s_y_j);
            for (R_xlen_t i = 0; i < n; ++i) {
                Y_cols[static_cast<size_t>(j)][static_cast<size_t>(i)] = yp[i];
            }
            UNPROTECT(1); // s_y_j
        }
    } else {
        Rf_error("Y must be either a numeric matrix (n x p) or a list of numeric vectors");
    }

    // Scalars (defensive extraction)
    const int   n_bws                     = Rf_asInteger(s_n_bws);
    const int   log_grid                  = (Rf_asLogical(s_log_grid) == TRUE) ? 1 : 0;
    const double min_bw_factor            = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor            = Rf_asReal(s_max_bw_factor);
    const double dist_norm_factor         = Rf_asReal(s_dist_normalization_factor);
    const int   ikernel                   = Rf_asInteger(s_kernel_type);
    const double precision                = Rf_asReal(s_precision);
    const int   n_folds                   = Rf_asInteger(s_n_folds);
    const int   buffer_hops               = Rf_asInteger(s_buffer_hops);
    const int   auto_buffer_hops          = (Rf_asLogical(s_auto_buffer_hops) == TRUE) ? 1 : 0;
    const int   use_uniform_weights       = (Rf_asLogical(s_use_uniform_weights) == TRUE) ? 1 : 0;
    const int   with_bw_predictions       = (Rf_asLogical(s_with_bw_predictions) == TRUE) ? 1 : 0;
    const int   verbose                   = (Rf_asLogical(s_verbose) == TRUE) ? 1 : 0;

    // Call core computation
    graph_deg0_lowess_cv_mat_result_t res = graph_deg0_lowess_cv_mat(
        adj_list,
        weight_list,
        Y_cols,
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        dist_norm_factor,
        ikernel,
        precision,
        n_folds,
        buffer_hops,
        auto_buffer_hops,
        use_uniform_weights,
        with_bw_predictions,
        verbose
    );

    // Assemble R result: list(cv_errors, mean_cv_errors, opt_bw_idxs, bw_factors)
    const int N = 4;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N)); // [1]

    // names
    {
        SEXP nm = PROTECT(Rf_allocVector(STRSXP, N)); // [2]
        SET_STRING_ELT(nm, 0, Rf_mkChar("cv_errors"));
        SET_STRING_ELT(nm, 1, Rf_mkChar("mean_cv_errors"));
        SET_STRING_ELT(nm, 2, Rf_mkChar("opt_bw_idxs"));
        SET_STRING_ELT(nm, 3, Rf_mkChar("bw_factors"));
        Rf_setAttrib(r_result, R_NamesSymbol, nm);
        UNPROTECT(1); // nm -> [1]
    }

    // 0) cv_errors: matrix (n_bws x n_cols), column-major
    {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.n_bws);
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.n_cols);
        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); // [2]
        double* Mp = REAL(mat);
        const size_t expected = static_cast<size_t>(nrow) * static_cast<size_t>(ncol);
        const size_t have = res.cv_errors.size();
        const size_t copyN = std::min(expected, have);
        for (size_t k = 0; k < copyN; ++k) Mp[k] = res.cv_errors[k];
        // If have < expected, remaining cells stay zero-initialized.
        SET_VECTOR_ELT(r_result, 0, mat);
        UNPROTECT(1); // mat -> [1]
    }

    // 1) mean_cv_errors: matrix (n_bws x n_cols) to mirror cv_errors shape
    {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.n_bws);
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.n_cols);
        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); // [2]
        double* Mp = REAL(mat);
        const size_t expected = static_cast<size_t>(nrow) * static_cast<size_t>(ncol);
        const size_t have = res.mean_cv_errors.size();
        const size_t copyN = std::min(expected, have);
        for (size_t k = 0; k < copyN; ++k) Mp[k] = res.mean_cv_errors[k];
        SET_VECTOR_ELT(r_result, 1, mat);
        UNPROTECT(1); // mat -> [1]
    }

    // 2) opt_bw_idxs: integer vector (length = n_cols), 1-based for R
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.opt_bw_idxs.size());
        SEXP v = PROTECT(Rf_allocVector(INTSXP, n)); // [2]
        int* ip = INTEGER(v);
        for (R_xlen_t j = 0; j < n; ++j) {
            int idx0 = res.opt_bw_idxs[static_cast<size_t>(j)];
            ip[j] = (idx0 >= 0) ? (idx0 + 1) : NA_INTEGER;
        }
        SET_VECTOR_ELT(r_result, 2, v);
        UNPROTECT(1); // v -> [1]
    }

    // 3) bw_factors: numeric vector (length = n_bws)
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.bw_factors.size());
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(res.bw_factors.begin(), res.bw_factors.end(), REAL(v));
        SET_VECTOR_ELT(r_result, 3, v);
        UNPROTECT(1); // v -> [1]
    }

    UNPROTECT(1); // r_result
    return r_result;
}
