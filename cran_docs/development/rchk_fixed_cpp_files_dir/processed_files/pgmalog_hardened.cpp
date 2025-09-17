/**
 * @brief Hardened rchk-safe drop-in of S_upgmalog
 *
 * Changes vs. provided version:
 * - Replaced direct INTEGER/REAL/LOGICAL indexing with Rf_asInteger/Rf_asReal/Rf_asLogical.
 * - Kept `result` and `names` protected until the tail; now finishes with UNPROTECT(2).
 * - Preserved container-first assembly with local PROTECT/UNPROTECT(1) for each element.
 * - Used XLENGTH/R_xlen_t for input vector lengths.
 */
#include <vector>
#include <numeric>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);

// Include necessary headers
#include "pgmalog.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations (assumed unchanged)
struct upgmalog_results {
    std::vector<int> h_values;
    int opt_h_idx;
    double opt_h;
    std::vector<double> h_cv_errors;
    std::vector<double> true_errors;
    std::vector<double> predictions;
    std::vector<double> local_predictions;
    std::vector<std::vector<double>> h_predictions;
    std::vector<double> bb_predictions;
    std::vector<double> ci_lower;
    std::vector<double> ci_upper;
    
    struct {
        std::vector<std::vector<int>> adj_list;
        std::vector<std::vector<double>> weight_list;
    } graph;
};

// Core computation function (assumed available)
upgmalog_results upgmalog(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const std::vector<double>& y_true,
                          bool use_median, int h_min, int h_max, double p,
                          int n_bb, int bb_max_distance_deviation,
                          int n_CVs, int n_CV_folds, unsigned int seed,
                          int ikernel, int n_cores,
                          double dist_normalization_factor,
                          double epsilon, bool verbose);

extern "C" {

SEXP S_upgmalog(SEXP s_x,
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
                SEXP s_n_cores,
                SEXP s_dist_normalization_factor,
                SEXP s_epsilon,
                SEXP s_verbose) {

    // Copy inputs into C++ containers (no coercion; assume numeric/logical inputs are correct types)
    const R_xlen_t nx = XLENGTH(s_x);
    const R_xlen_t ny = XLENGTH(s_y);
    if (nx != ny) {
        Rf_error("S_upgmalog: 'x' and 'y' must have the same length");
    }
    std::vector<double> x(REAL(s_x), REAL(s_x) + (size_t)nx);
    std::vector<double> y(REAL(s_y), REAL(s_y) + (size_t)ny);

    // Optional y_true
    std::vector<double> y_true;
    const R_xlen_t ny_true = XLENGTH(s_y_true);
    if (ny_true == nx) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + (size_t)ny_true);
    }

    // Defensive scalar parsing
    const bool  use_median = (Rf_asLogical(s_use_median) == TRUE);
    const int   h_min = Rf_asInteger(s_h_min);
    const int   h_max = Rf_asInteger(s_h_max);
    const double p = Rf_asReal(s_p);
    const int   n_bb = Rf_asInteger(s_n_bb);
    const int   bb_max_distance_deviation = Rf_asInteger(s_bb_max_distance_deviation);
    const int   n_CVs = Rf_asInteger(s_n_CVs);
    const int   n_CV_folds = Rf_asInteger(s_n_CV_folds);
    const unsigned int seed = (unsigned int) Rf_asInteger(s_seed);
    const int   ikernel = Rf_asInteger(s_ikernel);
    const int   n_cores = Rf_asInteger(s_n_cores);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const double epsilon = Rf_asReal(s_epsilon);
    const bool  verbose = (Rf_asLogical(s_verbose) == TRUE);

    // Core computation
    upgmalog_results cpp_results = upgmalog(x,
                                            y,
                                            y_true,
                                            use_median,
                                            h_min,
                                            h_max,
                                            p,
                                            n_bb,
                                            bb_max_distance_deviation,
                                            n_CVs,
                                            n_CV_folds,
                                            seed,
                                            ikernel,
                                            n_cores,
                                            dist_normalization_factor,
                                            epsilon,
                                            verbose);
    
    // Creating return list (container-first pattern)
    const int N_COMPONENTS = 13;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));

    // 0: h_values
    {
        SEXP s = PROTECT(convert_vector_int_to_R(cpp_results.h_values));
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1);
    }

    // 1: h_errors (or NULL)
    if (!cpp_results.h_cv_errors.empty() && nx > 0) {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.h_cv_errors));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    // 2: opt_h_idx (1-based, keep as REALSXP to match original behavior)
    {
        SEXP s_opt_h_idx = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s_opt_h_idx)[0] = (double)(cpp_results.opt_h_idx + 1);
        SET_VECTOR_ELT(result, 2, s_opt_h_idx);
        UNPROTECT(1);
    }

    // 3: opt_h
    {
        SEXP s_opt_h = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s_opt_h)[0] = cpp_results.opt_h;
        SET_VECTOR_ELT(result, 3, s_opt_h);
        UNPROTECT(1);
    }

    // 4: graph_adj_list
    {
        SEXP s = PROTECT(convert_vector_vector_int_to_R(cpp_results.graph.adj_list));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }

    // 5: graph_edge_lengths
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(cpp_results.graph.weight_list));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }

    // 6: predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.predictions));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    // 7: local_predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.local_predictions));
        SET_VECTOR_ELT(result, 7, s);
        UNPROTECT(1);
    }

    // 8-10: bootstrap results (or NULLs)
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
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    // 11: true_error (mean of vector, or NULL)
    if (!cpp_results.true_errors.empty()) {
        SEXP s_true_error = PROTECT(Rf_allocVector(REALSXP, 1));
        double mean_true_error = std::accumulate(cpp_results.true_errors.begin(),
                                                 cpp_results.true_errors.end(), 0.0) / 
                                 (double) cpp_results.true_errors.size();
        REAL(s_true_error)[0] = mean_true_error;
        SET_VECTOR_ELT(result, 11, s_true_error);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 11, R_NilValue);
    }

    // 12: h_predictions
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(cpp_results.h_predictions));
        SET_VECTOR_ELT(result, 12, s);
        UNPROTECT(1);
    }

    // Setting names for return list — keep protected until tail
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
    SET_STRING_ELT(names, 0, Rf_mkChar("h_values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("h_errors"));
    SET_STRING_ELT(names, 2, Rf_mkChar("opt_h_idx"));
    SET_STRING_ELT(names, 3, Rf_mkChar("opt_h"));
    SET_STRING_ELT(names, 4, Rf_mkChar("graph_adj_list"));
    SET_STRING_ELT(names, 5, Rf_mkChar("graph_edge_lengths"));
    SET_STRING_ELT(names, 6, Rf_mkChar("predictions"));
    SET_STRING_ELT(names, 7, Rf_mkChar("local_predictions"));
    SET_STRING_ELT(names, 8, Rf_mkChar("bb_predictions"));
    SET_STRING_ELT(names, 9, Rf_mkChar("ci_lower"));
    SET_STRING_ELT(names, 10, Rf_mkChar("ci_upper"));
    SET_STRING_ELT(names, 11, Rf_mkChar("true_error"));
    SET_STRING_ELT(names, 12, Rf_mkChar("h_predictions"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
