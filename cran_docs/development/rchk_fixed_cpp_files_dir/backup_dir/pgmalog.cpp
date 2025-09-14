/**
 * @brief Fixed version of pgmalog.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - S_upgmalog
 * 
 * Issues fixed:
 * 1. S_upgmalog (line 1269): UNPROTECT(variable) - n_protected
 * 2. S_upgmalog: Multiple unprotected allocations (lines 1206-1248)
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with literal constants
 * 2. Protected all allocations immediately
 * 3. Used container-first pattern consistently
 * 4. Fixed the pattern of incrementing n_protected and calling convert functions
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

// Core function structures and declarations
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

/**
 * Fixed version of S_upgmalog
 * Fixes: UNPROTECT(variable) at line 1269, unprotected allocations
 * Solution: Use container-first pattern with immediate protection
 */
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

    // Extract data vectors (no coercion needed if we trust input types)
    int n_points = LENGTH(s_x);
    std::vector<double> x(REAL(s_x), REAL(s_x) + n_points);
    std::vector<double> y(REAL(s_y), REAL(s_y) + n_points);

    // Handle empty y_true vector
    std::vector<double> y_true;
    if (LENGTH(s_y_true) == n_points) {
        y_true.assign(REAL(s_y_true), REAL(s_y_true) + LENGTH(s_y_true));
    }

    // Extract scalar parameters
    bool use_median = (LOGICAL(s_use_median)[0] == 1);
    int h_min = INTEGER(s_h_min)[0];
    int h_max = INTEGER(s_h_max)[0];
    double p = REAL(s_p)[0];
    int n_bb = INTEGER(s_n_bb)[0];
    int bb_max_distance_deviation = INTEGER(s_bb_max_distance_deviation)[0];
    int n_CVs = INTEGER(s_n_CVs)[0];
    int n_CV_folds = INTEGER(s_n_CV_folds)[0];
    unsigned int seed = (unsigned int)INTEGER(s_seed)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    int n_cores = INTEGER(s_n_cores)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    double epsilon = REAL(s_epsilon)[0];
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    // Core computation
    auto cpp_results = upgmalog(x,
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
    if (!cpp_results.h_cv_errors.empty() && n_points > 0) {
        SEXP s = PROTECT(convert_vector_double_to_R(cpp_results.h_cv_errors));
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    // 2: opt_h_idx (1-based)
    {
        SEXP s_opt_h_idx = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s_opt_h_idx)[0] = cpp_results.opt_h_idx + 1;
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
    if (cpp_results.bb_predictions.size() > 0) {
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
    if (cpp_results.true_errors.size() > 0) {
        SEXP s_true_error = PROTECT(Rf_allocVector(REALSXP, 1));
        double mean_true_error = std::accumulate(cpp_results.true_errors.begin(),
                                                 cpp_results.true_errors.end(), 0.0) / 
                                 cpp_results.true_errors.size();
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

    // Setting names for return list
    {
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
        UNPROTECT(1); // names
    }

    UNPROTECT(1); // result
    return result;
}

} // extern "C"