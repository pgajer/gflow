// pgmalo_hardened.cpp — rchk-hardened drop-ins for S_pgmalo and S_upgmalo
//
// Hardenings vs current file:
// - Keep only `result` + `names` protected until the tail (final UNPROTECT(2)).
// - Use PROTECT_WITH_INDEX/REPROTECT coercion blocks for numeric inputs; copy to STL then UNPROTECT immediately.
// - Defensive scalar extraction via Rf_as* + NA/range checks across all scalars.
// - Container-first conversions before any result allocations.
// - Local PROTECT/UNPROTECT(1) around each element inserted into the result list.
// - 1-based indices returned to R where applicable.
//
// Assumptions: core algorithms and converters exist with the same signatures used here.

#include <vector>
#include <numeric>
#include <R.h>
#include <Rinternals.h>

// Converters (provided elsewhere)
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);

// Core API (provided elsewhere)
struct pgmalo_results {
    std::vector<int>    h_values;
    int                 opt_h;
    int                 opt_h_idx;            // 0-based
    std::vector<double> h_cv_errors;          // aka h_errors
    std::vector<double> true_errors;
    std::vector<double> predictions;
    std::vector<double> local_predictions;
    std::vector<std::vector<double>> h_predictions;
    std::vector<double> bb_predictions;
    std::vector<double> ci_lower;
    std::vector<double> ci_upper;
};

// Graph-based version
pgmalo_results pgmalo(const std::vector<std::vector<int>>& neighbors,
                      const std::vector<std::vector<double>>& edge_lengths,
                      const std::vector<double>& y,
                      const std::vector<double>* y_true_opt, // nullptr if none
                      int use_median,
                      int h_min,
                      int h_max,
                      int p,
                      int n_bb,
                      double bb_max_distance_deviation,
                      int n_CVs,
                      int n_CV_folds,
                      int seed,
                      int kernel_type,
                      double dist_normalization_factor,
                      double epsilon,
                      int verbose);

// Univariate (x,y) version
pgmalo_results upgmalo(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<double>* y_true_opt, // nullptr if none
                       int use_median,
                       int h_min,
                       int h_max,
                       int p,
                       int n_bb,
                       double bb_max_distance_deviation,
                       int n_CVs,
                       int n_CV_folds,
                       int seed,
                       int kernel_type,
                       double dist_normalization_factor,
                       double epsilon,
                       int verbose);

extern "C" {

// ------------------------------- S_pgmalo ----------------------------------
SEXP S_pgmalo(SEXP neighbors_r,
              SEXP edge_lengths_r,
              SEXP y_r,
              SEXP y_true_r,
              SEXP use_median_r,
              SEXP h_min_r,
              SEXP h_max_r,
              SEXP p_r,
              SEXP n_bb_r,
              SEXP bb_max_distance_deviation_r,
              SEXP n_CVs_r,
              SEXP n_CV_folds_r,
              SEXP seed_r,
              SEXP kernel_type_r,
              SEXP dist_normalization_factor_r,
              SEXP epsilon_r,
              SEXP verbose_r) {

    // Container-first conversions
    std::vector<std::vector<int>>    neighbors    = convert_adj_list_from_R(neighbors_r);
    std::vector<std::vector<double>> edge_lengths = convert_weight_list_from_R(edge_lengths_r);

    // y (coercion block)
    std::vector<double> y;
    {
        SEXP sy = y_r;
        PROTECT_INDEX ipy;
        PROTECT_WITH_INDEX(sy, &ipy);
        if (TYPEOF(sy) != REALSXP) {
            REPROTECT(sy = Rf_coerceVector(sy, REALSXP), ipy);
        }
        const R_xlen_t Ny = XLENGTH(sy);
        if (Ny <= 0) { UNPROTECT(1); Rf_error("S_pgmalo(): 'y' must be a non-empty numeric vector."); }
        y.assign(REAL(sy), REAL(sy) + (size_t)Ny);
        UNPROTECT(1); // sy
    }

    // optional y_true
    std::vector<double> y_true;
    const std::vector<double>* y_true_ptr = nullptr;
    if (y_true_r != R_NilValue) {
        SEXP syt = y_true_r;
        PROTECT_INDEX ipyt;
        PROTECT_WITH_INDEX(syt, &ipyt);
        if (TYPEOF(syt) != REALSXP) {
            REPROTECT(syt = Rf_coerceVector(syt, REALSXP), ipyt);
        }
        const R_xlen_t Nyt = XLENGTH(syt);
        y_true.assign(REAL(syt), REAL(syt) + (size_t)Nyt);
        UNPROTECT(1); // syt
        if ((R_xlen_t)y_true.size() != (R_xlen_t)y.size()) {
            Rf_error("S_pgmalo(): 'y_true' must have the same length as 'y'.");
        }
        y_true_ptr = &y_true;
    }

    if (neighbors.size() != y.size()) {
        Rf_error("S_pgmalo(): length(y) must equal number of graph vertices (%llu).",
                 (unsigned long long)neighbors.size());
    }
    if (edge_lengths.size() != neighbors.size()) {
        Rf_error("S_pgmalo(): 'edge_lengths' must match 'neighbors' in length.");
    }

    // Scalars (defensive)
    const int    use_median  = Rf_asLogical(use_median_r);
    const int    h_min_i     = Rf_asInteger(h_min_r);
    const int    h_max_i     = Rf_asInteger(h_max_r);
    const int    p_i         = Rf_asInteger(p_r);
    const int    n_bb_i      = Rf_asInteger(n_bb_r);
    const double bb_dev      = Rf_asReal(bb_max_distance_deviation_r);
    const int    n_CVs_i     = Rf_asInteger(n_CVs_r);
    const int    n_folds_i   = Rf_asInteger(n_CV_folds_r);
    const int    seed_i      = Rf_asInteger(seed_r);
    const int    ikernel_i   = Rf_asInteger(kernel_type_r);
    const double dnorm_fac   = Rf_asReal(dist_normalization_factor_r);
    const double eps         = Rf_asReal(epsilon_r);
    const int    verbose_i   = Rf_asLogical(verbose_r);

    // NA / range checks
    if (use_median  == NA_LOGICAL) Rf_error("S_pgmalo(): 'use_median' must be TRUE/FALSE.");
    if (h_min_i     == NA_INTEGER || h_min_i <= 0) Rf_error("S_pgmalo(): 'h_min' must be a positive integer.");
    if (h_max_i     == NA_INTEGER || h_max_i <  h_min_i) Rf_error("S_pgmalo(): 'h_max' must be >= h_min.");
    if (p_i         == NA_INTEGER || p_i <= 0) Rf_error("S_pgmalo(): 'p' must be a positive integer.");
    if (n_bb_i      == NA_INTEGER || n_bb_i < 0) Rf_error("S_pgmalo(): 'n_bb' must be >= 0.");
    if (ISNAN(bb_dev) || bb_dev < 0.0) Rf_error("S_pgmalo(): 'bb_max_distance_deviation' must be >= 0.");
    if (n_CVs_i     == NA_INTEGER || n_CVs_i <= 0) Rf_error("S_pgmalo(): 'n_CVs' must be a positive integer.");
    if (n_folds_i   == NA_INTEGER || n_folds_i <= 0) Rf_error("S_pgmalo(): 'n_CV_folds' must be a positive integer.");
    if (seed_i      == NA_INTEGER) Rf_error("S_pgmalo(): 'seed' cannot be NA.");
    if (ikernel_i   == NA_INTEGER || ikernel_i < 0) Rf_error("S_pgmalo(): 'kernel_type' must be a non-negative integer.");
    if (ISNAN(dnorm_fac) || dnorm_fac < 0.0) Rf_error("S_pgmalo(): 'dist_normalization_factor' must be >= 0.");
    if (ISNAN(eps) || eps <= 0.0) Rf_error("S_pgmalo(): 'epsilon' must be > 0.");
    if (verbose_i   == NA_LOGICAL) Rf_error("S_pgmalo(): 'verbose' must be TRUE/FALSE.");

    // Core computation
    pgmalo_results cpp_results = pgmalo(
        neighbors, edge_lengths, y, y_true_ptr,
        use_median, h_min_i, h_max_i, p_i, n_bb_i, bb_dev,
        n_CVs_i, n_folds_i, seed_i, ikernel_i, dnorm_fac, eps, verbose_i
    );

    // Build result
    const int N_COMPONENTS = 11;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
    SET_STRING_ELT(names, 0,  Rf_mkChar("h_values"));
    SET_STRING_ELT(names, 1,  Rf_mkChar("opt_h"));
    SET_STRING_ELT(names, 2,  Rf_mkChar("opt_h_idx"));
    SET_STRING_ELT(names, 3,  Rf_mkChar("h_errors"));
    SET_STRING_ELT(names, 4,  Rf_mkChar("true_errors"));
    SET_STRING_ELT(names, 5,  Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 6,  Rf_mkChar("local_predictions"));
    SET_STRING_ELT(names, 7,  Rf_mkChar("h_predictions"));
    SET_STRING_ELT(names, 8,  Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 9,  Rf_mkChar("ci_lower"));
    SET_STRING_ELT(names, 10, Rf_mkChar("ci_upper"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0 h_values
    {
        SEXP s = PROTECT(convert_vector_int_to_R(cpp_results.h_values));
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1);
    }
    // 1 opt_h
    {
        SEXP s = PROTECT(Rf_ScalarInteger(cpp_results.opt_h));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }
    // 2 opt_h_idx (1-based)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(cpp_results.opt_h_idx + 1));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }
    // 3 h_errors
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.h_cv_errors));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }
    // 4 true_errors
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.true_errors));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }
    // 5 predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.predictions));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }
    // 6 local_predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.local_predictions));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }
    // 7 h_predictions (list of numeric)
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(cpp_results.h_predictions));
        SET_VECTOR_ELT(result, 7, s);
        UNPROTECT(1);
    }
    // 8..10: optional bb/ci
    if (!cpp_results.bb_predictions.empty()) {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.bb_predictions));
        SET_VECTOR_ELT(result, 8, s);
        UNPROTECT(1);
        s = PROTECT(convert_vector_double_to_R(cpp_results.ci_lower));
        SET_VECTOR_ELT(result, 9, s);
        UNPROTECT(1);
        s = PROTECT(convert_vector_double_to_R(cpp_results.ci_upper));
        SET_VECTOR_ELT(result, 10, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 8,  R_NilValue);
        SET_VECTOR_ELT(result, 9,  R_NilValue);
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    UNPROTECT(2); // result, names
    return result;
}

// ------------------------------ S_upgmalo ----------------------------------
SEXP S_upgmalo(SEXP s_x,
               SEXP s_y,
               SEXP s_y_true,
               SEXP s_use_median,
               SEXP s_h_min,
               SEXP s_h_max,
               SEXP s_p,
               SEXP s_n_bb,
               SEXP s_bb_max_distance_deviation,
               SEXP s_n_CVs,
               SEXP s_n_CV_folds,
               SEXP s_seed,
               SEXP s_ikernel,
               SEXP s_dist_normalization_factor,
               SEXP s_epsilon,
               SEXP s_verbose) {

    // x, y coercion blocks
    std::vector<double> x, y;
    {
        SEXP sx = s_x;
        PROTECT_INDEX ipx;
        PROTECT_WITH_INDEX(sx, &ipx);
        if (TYPEOF(sx) != REALSXP) {
            REPROTECT(sx = Rf_coerceVector(sx, REALSXP), ipx);
        }
        const R_xlen_t nx = XLENGTH(sx);
        if (nx <= 0) { UNPROTECT(1); Rf_error("S_upgmalo(): 'x' must be a non-empty numeric vector."); }
        x.assign(REAL(sx), REAL(sx) + (size_t)nx);
        UNPROTECT(1); // sx

        SEXP sy = s_y;
        PROTECT_INDEX ipy;
        PROTECT_WITH_INDEX(sy, &ipy);
        if (TYPEOF(sy) != REALSXP) {
            REPROTECT(sy = Rf_coerceVector(sy, REALSXP), ipy);
        }
        const R_xlen_t ny = XLENGTH(sy);
        if (ny <= 0) { UNPROTECT(1); Rf_error("S_upgmalo(): 'y' must be a non-empty numeric vector."); }
        y.assign(REAL(sy), REAL(sy) + (size_t)ny);
        UNPROTECT(1); // sy

        if (nx != ny) {
            Rf_error("S_upgmalo(): 'x' and 'y' must have the same length.");
        }
    }

    // optional y_true
    std::vector<double> y_true;
    const std::vector<double>* y_true_ptr = nullptr;
    if (s_y_true != R_NilValue) {
        SEXP syt = s_y_true;
        PROTECT_INDEX ipyt;
        PROTECT_WITH_INDEX(syt, &ipyt);
        if (TYPEOF(syt) != REALSXP) {
            REPROTECT(syt = Rf_coerceVector(syt, REALSXP), ipyt);
        }
        const R_xlen_t Nyt = XLENGTH(syt);
        y_true.assign(REAL(syt), REAL(syt) + (size_t)Nyt);
        UNPROTECT(1); // syt
        if ((R_xlen_t)y_true.size() != (R_xlen_t)y.size()) {
            Rf_error("S_upgmalo(): 'y_true' must have the same length as 'y'.");
        }
        y_true_ptr = &y_true;
    }

    // Scalars (defensive)
    const int    use_median  = Rf_asLogical(s_use_median);
    const int    h_min_i     = Rf_asInteger(s_h_min);
    const int    h_max_i     = Rf_asInteger(s_h_max);
    const int    p_i         = Rf_asInteger(s_p);
    const int    n_bb_i      = Rf_asInteger(s_n_bb);
    const double bb_dev      = Rf_asReal(s_bb_max_distance_deviation);
    const int    n_CVs_i     = Rf_asInteger(s_n_CVs);
    const int    n_folds_i   = Rf_asInteger(s_n_CV_folds);
    const int    seed_i      = Rf_asInteger(s_seed);
    const int    ikernel_i   = Rf_asInteger(s_ikernel);
    const double dnorm_fac   = Rf_asReal(s_dist_normalization_factor);
    const double eps         = Rf_asReal(s_epsilon);
    const int    verbose_i   = Rf_asLogical(s_verbose);

    // NA / range checks
    if (use_median  == NA_LOGICAL) Rf_error("S_upgmalo(): 'use_median' must be TRUE/FALSE.");
    if (h_min_i     == NA_INTEGER || h_min_i <= 0) Rf_error("S_upgmalo(): 'h_min' must be a positive integer.");
    if (h_max_i     == NA_INTEGER || h_max_i <  h_min_i) Rf_error("S_upgmalo(): 'h_max' must be >= h_min.");
    if (p_i         == NA_INTEGER || p_i <= 0) Rf_error("S_upgmalo(): 'p' must be a positive integer.");
    if (n_bb_i      == NA_INTEGER || n_bb_i < 0) Rf_error("S_upgmalo(): 'n_bb' must be >= 0.");
    if (ISNAN(bb_dev) || bb_dev < 0.0) Rf_error("S_upgmalo(): 'bb_max_distance_deviation' must be >= 0.");
    if (n_CVs_i     == NA_INTEGER || n_CVs_i <= 0) Rf_error("S_upgmalo(): 'n_CVs' must be a positive integer.");
    if (n_folds_i   == NA_INTEGER || n_folds_i <= 0) Rf_error("S_upgmalo(): 'n_CV_folds' must be a positive integer.");
    if (seed_i      == NA_INTEGER) Rf_error("S_upgmalo(): 'seed' cannot be NA.");
    if (ikernel_i   == NA_INTEGER || ikernel_i < 0) Rf_error("S_upgmalo(): 'kernel_type' must be a non-negative integer.");
    if (ISNAN(dnorm_fac) || dnorm_fac < 0.0) Rf_error("S_upgmalo(): 'dist_normalization_factor' must be >= 0.");
    if (ISNAN(eps) || eps <= 0.0) Rf_error("S_upgmalo(): 'epsilon' must be > 0.");
    if (verbose_i   == NA_LOGICAL) Rf_error("S_upgmalo(): 'verbose' must be TRUE/FALSE.");

    // Core computation
    pgmalo_results cpp_results = upgmalo(
        x, y, y_true_ptr,
        use_median, h_min_i, h_max_i, p_i, n_bb_i, bb_dev,
        n_CVs_i, n_folds_i, seed_i, ikernel_i, dnorm_fac, eps, verbose_i
    );

    // Build result (includes scalar mean true_error)
    const int N = 13;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0,  Rf_mkChar("h_values"));
    SET_STRING_ELT(names, 1,  Rf_mkChar("opt_h"));
    SET_STRING_ELT(names, 2,  Rf_mkChar("opt_h_idx"));
    SET_STRING_ELT(names, 3,  Rf_mkChar("h_errors"));
    SET_STRING_ELT(names, 4,  Rf_mkChar("true_errors"));
    SET_STRING_ELT(names, 5,  Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 6,  Rf_mkChar("local_predictions"));
    SET_STRING_ELT(names, 7,  Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 8,  Rf_mkChar("ci_lower"));
    SET_STRING_ELT(names, 9,  Rf_mkChar("ci_upper"));
    SET_STRING_ELT(names, 10, Rf_mkChar("true_error"));
    SET_STRING_ELT(names, 11, Rf_mkChar("h_predictions"));
    SET_STRING_ELT(names, 12, Rf_mkChar("x")); // echo x grid back if desired; optional
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0..3 basics
    {
        SEXP s = PROTECT(convert_vector_int_to_R(cpp_results.h_values));
        SET_VECTOR_ELT(result, 0, s); UNPROTECT(1);
    }
    { SEXP s = PROTECT(Rf_ScalarInteger(cpp_results.opt_h)); SET_VECTOR_ELT(result, 1, s); UNPROTECT(1); }
    { SEXP s = PROTECT(Rf_ScalarInteger(cpp_results.opt_h_idx + 1)); SET_VECTOR_ELT(result, 2, s); UNPROTECT(1); }
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.h_cv_errors));
        SET_VECTOR_ELT(result, 3, s); UNPROTECT(1);
    }

    // 4 true_errors (vector)
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.true_errors));
        SET_VECTOR_ELT(result, 4, s); UNPROTECT(1);
    }

    // 5 predictions
    { SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.predictions)); SET_VECTOR_ELT(result, 5, s); UNPROTECT(1); }

    // 6 local_predictions
    { SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.local_predictions)); SET_VECTOR_ELT(result, 6, s); UNPROTECT(1); }

    // 7..9: bb + ci (optional)
    if (!cpp_results.bb_predictions.empty()) {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.bb_predictions)); SET_VECTOR_ELT(result, 7, s); UNPROTECT(1);
        s = PROTECT(convert_vector_double_to_R(cpp_results.ci_lower));          SET_VECTOR_ELT(result, 8, s); UNPROTECT(1);
        s = PROTECT(convert_vector_double_to_R(cpp_results.ci_upper));          SET_VECTOR_ELT(result, 9, s); UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 7, R_NilValue);
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }

    // 10: true_error (mean of vector)
    {
        double mean_err = NA_REAL;
        if (!cpp_results.true_errors.empty()) {
            mean_err = std::accumulate(cpp_results.true_errors.begin(), cpp_results.true_errors.end(), 0.0)
                     / (double)cpp_results.true_errors.size();
        }
        SEXP s = PROTECT(Rf_ScalarReal(mean_err));
        SET_VECTOR_ELT(result, 10, s);
        UNPROTECT(1);
    }

    // 11: h_predictions (list of numeric)
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(cpp_results.h_predictions));
        SET_VECTOR_ELT(result, 11, s);
        UNPROTECT(1);
    }

    // 12: x (echo back; optional)
    {
        SEXP s = PROTECT(convert_vector_double_to_R(x));
        SET_VECTOR_ELT(result, 12, s);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
