/**
 * @file graph_diffusion_smoother_hardened.cpp
 * @brief Hardened rchk-safe S_graph_diffusion_smoother (drop-in replacement).
 *
 * Changes vs. corrected version:
 *  - Defensive coercion of s_y via PROTECT_WITH_INDEX / REPROTECT;
 *    copy to std::vector<double>, then UNPROTECT(1).
 *  - Long-vector safety using XLENGTH() and R_xlen_t.
 *  - Scalar extraction via Rf_asInteger / Rf_asReal / Rf_asLogical.
 *  - Container-first assembly, only fixed-count UNPROTECTs.
 *
 * Assumptions:
 *  - convert_adj_list_from_R / convert_weight_list_from_R do not allocate R objects
 *    or, if they do, they manage PROTECT/UNPROTECT internally and return plain C++ containers.
 *  - graph_diffusion_smoother(...) performs no R allocations.
 */

#include <vector>
#include <algorithm>

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

// Helpers provided elsewhere
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Core result & function (declared elsewhere)
struct graph_diffusion_smoother_result_t {
    std::vector<std::vector<double>> y_traj;
    std::vector<double> cv_errors;
    std::vector<double> mean_cv_errors;
    std::vector<double> y_optimal;
    int n_time_steps;
    int n_CVs;
    int optimal_time_step;
    double min_cv_error;
};

graph_diffusion_smoother_result_t graph_diffusion_smoother(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<double>& y,
    int n_time_steps,
    double step_factor,
    double binary_threshold,
    int ikernel,
    double dist_normalization_factor,
    int n_CVs,
    int n_CV_folds,
    double epsilon,
    bool verbose,
    unsigned int seed);

extern "C" SEXP S_graph_diffusion_smoother(SEXP s_adj_list,
                                           SEXP s_weight_list,
                                           SEXP s_y,
                                           SEXP s_n_time_steps,
                                           SEXP s_step_factor,
                                           SEXP s_binary_threshold,
                                           SEXP s_ikernel,
                                           SEXP s_dist_normalization_factor,
                                           SEXP s_n_CVs,
                                           SEXP s_n_CV_folds,
                                           SEXP s_epsilon,
                                           SEXP s_verbose,
                                           SEXP s_seed)
{
    // Convert graph inputs (assumed to be pure C++)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Defensive coercion for y (REAL), long-vector safe
    std::vector<double> y;
    {
        SEXP sy = s_y;
        PROTECT_INDEX py;
        PROTECT_WITH_INDEX(sy, &py);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);
        const R_xlen_t ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));
        UNPROTECT(1); // sy
    }

    // Scalars (defensive extraction)
    const int   n_time_steps             = Rf_asInteger(s_n_time_steps);
    const double step_factor             = Rf_asReal(s_step_factor);
    const double binary_threshold        = Rf_asReal(s_binary_threshold);
    const int   ikernel                  = Rf_asInteger(s_ikernel);
    const double dist_normalization_fact = Rf_asReal(s_dist_normalization_factor);
    const int   n_CVs                    = Rf_asInteger(s_n_CVs);
    const int   n_CV_folds               = Rf_asInteger(s_n_CV_folds);
    const double epsilon                 = Rf_asReal(s_epsilon);
    const bool  verbose                  = (Rf_asLogical(s_verbose) == TRUE);
    const unsigned int seed              = static_cast<unsigned int>(Rf_asInteger(s_seed));

    // Call core computation
    graph_diffusion_smoother_result_t result = graph_diffusion_smoother(
        adj_list,
        weight_list,
        y,
        n_time_steps,
        step_factor,
        binary_threshold,
        ikernel,
        dist_normalization_fact,
        n_CVs,
        n_CV_folds,
        epsilon,
        verbose,
        seed
    );

    // Assemble return value
    const int N = 8;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N)); // [1]

    // names
    {
        SEXP nm = PROTECT(Rf_allocVector(STRSXP, N)); // [2]
        SET_STRING_ELT(nm, 0, Rf_mkChar("y_traj"));
        SET_STRING_ELT(nm, 1, Rf_mkChar("cv_errors"));
        SET_STRING_ELT(nm, 2, Rf_mkChar("mean_cv_errors"));
        SET_STRING_ELT(nm, 3, Rf_mkChar("y_optimal"));
        SET_STRING_ELT(nm, 4, Rf_mkChar("n_time_steps"));
        SET_STRING_ELT(nm, 5, Rf_mkChar("n_CVs"));
        SET_STRING_ELT(nm, 6, Rf_mkChar("optimal_time_step"));
        SET_STRING_ELT(nm, 7, Rf_mkChar("min_cv_error"));
        Rf_setAttrib(r_result, R_NamesSymbol, nm);
        UNPROTECT(1); // nm -> [1]
    }

    // 0: y_traj — list of numeric vectors
    {
        const R_xlen_t T = static_cast<R_xlen_t>(result.y_traj.size());
        SEXP r_y_traj = PROTECT(Rf_allocVector(VECSXP, T)); // [2]
        for (R_xlen_t t = 0; t < T; ++t) {
            const R_xlen_t n = static_cast<R_xlen_t>(result.y_traj[static_cast<size_t>(t)].size());
            SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, n)); // [3]
            std::copy(result.y_traj[static_cast<size_t>(t)].begin(),
                      result.y_traj[static_cast<size_t>(t)].end(),
                      REAL(r_vec));
            SET_VECTOR_ELT(r_y_traj, t, r_vec);
            UNPROTECT(1); // r_vec -> [2]
        }
        SET_VECTOR_ELT(r_result, 0, r_y_traj);
        UNPROTECT(1); // r_y_traj -> [1]
    }

    // 1: cv_errors — matrix (n_time_steps x n_CVs) or NULL
    if (result.n_CVs > 0 && result.n_time_steps > 0) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(result.n_time_steps);
        const R_xlen_t ncol = static_cast<R_xlen_t>(result.n_CVs);
        SEXP r_cv = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); // [2]
        double* p = REAL(r_cv);
        for (R_xlen_t i = 0; i < nrow; ++i) {
            for (R_xlen_t j = 0; j < ncol; ++j) {
                p[i + j * nrow] = result.cv_errors[static_cast<size_t>(i + j * nrow)];
            }
        }
        SET_VECTOR_ELT(r_result, 1, r_cv);
        UNPROTECT(1); // r_cv -> [1]
    } else {
        SET_VECTOR_ELT(r_result, 1, R_NilValue);
    }

    // 2: mean_cv_errors — numeric vector
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.mean_cv_errors.size());
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(result.mean_cv_errors.begin(), result.mean_cv_errors.end(), REAL(v));
        SET_VECTOR_ELT(r_result, 2, v);
        UNPROTECT(1); // v -> [1]
    }

    // 3: y_optimal — numeric vector
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.y_optimal.size());
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n)); // [2]
        std::copy(result.y_optimal.begin(), result.y_optimal.end(), REAL(v));
        SET_VECTOR_ELT(r_result, 3, v);
        UNPROTECT(1); // v -> [1]
    }

    // 4: n_time_steps — scalar int
    {
        SEXP v = PROTECT(Rf_ScalarInteger(result.n_time_steps)); // [2]
        SET_VECTOR_ELT(r_result, 4, v);
        UNPROTECT(1); // v -> [1]
    }

    // 5: n_CVs — scalar int
    {
        SEXP v = PROTECT(Rf_ScalarInteger(result.n_CVs)); // [2]
        SET_VECTOR_ELT(r_result, 5, v);
        UNPROTECT(1); // v -> [1]
    }

    // 6: optimal_time_step — scalar int
    {
        SEXP v = PROTECT(Rf_ScalarInteger(result.optimal_time_step)); // [2]
        SET_VECTOR_ELT(r_result, 6, v);
        UNPROTECT(1); // v -> [1]
    }

    // 7: min_cv_error — scalar double
    {
        SEXP v = PROTECT(Rf_ScalarReal(result.min_cv_error)); // [2]
        SET_VECTOR_ELT(r_result, 7, v);
        UNPROTECT(1); // v -> [1]
    }

    UNPROTECT(1); // r_result
    return r_result;
}
