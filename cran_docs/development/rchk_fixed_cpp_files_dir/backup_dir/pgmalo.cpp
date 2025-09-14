/**
 * @brief Fixed versions of pgmalo.cpp functions to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected versions of:
 * - S_pgmalo
 * - S_upgmalo
 * 
 * Issues fixed:
 * 1. S_pgmalo (lines 1238, 1251): Variable in UNPROTECT (tprot)
 * 2. S_upgmalo (lines 615, 628): Variable in UNPROTECT (tprot)
 * 
 * Changes made:
 * 1. Replaced variable UNPROTECT with PROTECT_WITH_INDEX/REPROTECT pattern
 * 2. Used literal constants for all UNPROTECT calls
 * 3. Maintained container-first pattern (already present)
 */

#include <vector>
#include <numeric>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);

// Include necessary headers
#include "pgmalo.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

// Core function structures and declarations
struct pgmalo_results {
    std::vector<int> h_values;
    int opt_h;
    int opt_h_idx;
    std::vector<double> h_cv_errors;
    std::vector<double> true_errors;
    std::vector<double> predictions;
    std::vector<double> local_predictions;
    std::vector<std::vector<double>> h_predictions;
    std::vector<double> bb_predictions;
    std::vector<double> ci_lower;
    std::vector<double> ci_upper;
    
    bool has_true_errors() const { return !true_errors.empty(); }
    bool has_bootstrap_results() const { return !bb_predictions.empty(); }
};

struct upgmalo_results {
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

// Core computation functions (assumed available)
pgmalo_results pgmalo(const std::vector<std::vector<int>>& neighbors,
                      const std::vector<std::vector<double>>& edge_lengths,
                      const std::vector<double>& y,
                      const std::vector<double>& y_true,
                      bool use_median, int h_min, int h_max, double p,
                      int n_bb, int bb_max_distance_deviation,
                      int n_CVs, int n_CV_folds, unsigned int seed,
                      int kernel_type, double dist_normalization_factor,
                      double epsilon, bool verbose);

upgmalo_results upgmalo(const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& y_true,
                        bool use_median, int h_min, int h_max, double p,
                        int n_bb, int bb_max_distance_deviation,
                        int n_CVs, int n_CV_folds, unsigned int seed,
                        int ikernel, double dist_normalization_factor,
                        double epsilon, bool verbose);

extern "C" {

/**
 * Fixed version of S_pgmalo
 * Fixes: Variable UNPROTECT at lines 1238, 1251
 * Solution: Use PROTECT_WITH_INDEX/REPROTECT pattern for conditional coercion
 */
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
    
    // ---- Graph inputs (pure reads; no PROTECT needed) ----
    std::vector<std::vector<int>>    neighbors    = convert_adj_list_from_R(neighbors_r);
    std::vector<std::vector<double>> edge_lengths;
    if (edge_lengths_r != R_NilValue) {
        edge_lengths = convert_weight_list_from_R(edge_lengths_r);
        if (!edge_lengths.empty() && edge_lengths.size() != neighbors.size()) {
            Rf_error("edge_lengths length must be 0 or equal to neighbors length.");
        }
        // Per-vertex length check when weights are supplied
        if (!edge_lengths.empty()) {
            const size_t V = neighbors.size();
            for (size_t i = 0; i < V; ++i) {
                if (edge_lengths[i].size() != neighbors[i].size()) {
                    Rf_error("edge_lengths[[%zu]] length (%zu) must equal neighbors[[%zu]] length (%zu).",
                             i + 1, edge_lengths[i].size(), i + 1, neighbors[i].size());
                }
            }
        }
    }

    // ---- y / y_true (coerce defensively using PROTECT_WITH_INDEX) ----
    std::vector<double> y, y_true;
    {
        PROTECT_INDEX ipx;
        SEXP sy = y_r;
        PROTECT_WITH_INDEX(sy, &ipx);
        
        if (TYPEOF(sy) != REALSXP) {
            REPROTECT(sy = Rf_coerceVector(sy, REALSXP), ipx);
        }
        const R_xlen_t Ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(Ny));

        if (neighbors.size() != static_cast<size_t>(Ny)) {
            UNPROTECT(1); // sy
            Rf_error("length(y) must equal length(neighbors).");
        }

        if (y_true_r != R_NilValue) {
            PROTECT_INDEX ipx2;
            SEXP syt = y_true_r;
            PROTECT_WITH_INDEX(syt, &ipx2);
            
            if (TYPEOF(syt) != REALSXP) {
                REPROTECT(syt = Rf_coerceVector(syt, REALSXP), ipx2);
            }
            const R_xlen_t Nyt = XLENGTH(syt);
            if (Nyt == Ny) {
                y_true.assign(REAL(syt), REAL(syt) + static_cast<size_t>(Nyt));
            } // else: leave empty (treated as unavailable)
            
            UNPROTECT(1); // syt
        }

        UNPROTECT(1); // sy
    }

    const R_xlen_t N = static_cast<R_xlen_t>(y.size());

    // ---- Scalars / parameters (validated) ----
    const bool   use_median = (Rf_asLogical(use_median_r) == TRUE);
    const int    h_min      = Rf_asInteger(h_min_r);
    const int    h_max      = Rf_asInteger(h_max_r);
    const double p          = Rf_asReal(p_r);
    const int    n_bb       = Rf_asInteger(n_bb_r);
    const int    bb_max_distance_deviation = Rf_asInteger(bb_max_distance_deviation_r);
    const int    n_CVs      = Rf_asInteger(n_CVs_r);
    const int    n_CV_folds = Rf_asInteger(n_CV_folds_r);
    const unsigned int seed  = static_cast<unsigned int>(Rf_asInteger(seed_r));
    const int    kernel_type = Rf_asInteger(kernel_type_r);
    const double dist_normalization_factor = Rf_asReal(dist_normalization_factor_r);
    const double epsilon     = Rf_asReal(epsilon_r);
    const bool   verbose     = (Rf_asLogical(verbose_r) == TRUE);

    if (h_min < 1)                 Rf_error("h_min must be >= 1.");
    if (h_max < h_min)             Rf_error("h_max must be >= h_min.");
    if (static_cast<R_xlen_t>(h_max) > N)
                                   Rf_error("h_max (%d) cannot exceed N (%lld).",
                                            h_max, static_cast<long long>(N));
    if (!(p > 0.0 && p <= 1.0))    Rf_error("p must be in (0, 1].");
    if (n_bb < 0)                  Rf_error("n_bb must be >= 0.");
    if (bb_max_distance_deviation < 0)
                                   Rf_error("bb_max_distance_deviation must be >= 0.");
    if (n_CVs < 1)                 Rf_error("n_CVs must be >= 1.");
    if (n_CV_folds < 2)            Rf_error("n_CV_folds must be >= 2.");
    if (!(epsilon > 0.0))          Rf_error("epsilon must be > 0.");

    // ---- Core computation (no R allocations inside) ----
    auto results = pgmalo(neighbors,
                          edge_lengths,
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
                          kernel_type,
                          dist_normalization_factor,
                          epsilon,
                          verbose);

    // ---- Build result (container-first; per-element PROTECT/UNPROTECT) ----
    const int N_COMPONENTS = 11;
    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));

    // 0: h_values
    {
        SEXP s = PROTECT(convert_vector_int_to_R(results.h_values));
        SET_VECTOR_ELT(result_r, 0, s);
        UNPROTECT(1);
    }

    // 1: opt_h (scalar int)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(results.opt_h));
        SET_VECTOR_ELT(result_r, 1, s);
        UNPROTECT(1);
    }

    // 2: opt_h_idx (if your convention is 1-based in R, add +1 here)
    {
        SEXP s = PROTECT(Rf_ScalarInteger(results.opt_h_idx));
        SET_VECTOR_ELT(result_r, 2, s);
        UNPROTECT(1);
    }

    // 3: h_errors
    {
        SEXP s = PROTECT(convert_vector_double_to_R(results.h_cv_errors));
        SET_VECTOR_ELT(result_r, 3, s);
        UNPROTECT(1);
    }

    // 4: true_errors (or NULL)
    if (results.has_true_errors()) {
        SEXP s = PROTECT(convert_vector_double_to_R(results.true_errors));
        SET_VECTOR_ELT(result_r, 4, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result_r, 4, R_NilValue);
    }

    // 5: predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(results.predictions));
        SET_VECTOR_ELT(result_r, 5, s);
        UNPROTECT(1);
    }

    // 6: local_predictions
    {
        SEXP s = PROTECT(convert_vector_double_to_R(results.local_predictions));
        SET_VECTOR_ELT(result_r, 6, s);
        UNPROTECT(1);
    }

    // 7: h_predictions (list<numeric>)
    {
        const size_t H = results.h_predictions.size();
        SEXP hp = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(H)));
        for (size_t i = 0; i < H; ++i) {
            SEXP vec = PROTECT(convert_vector_double_to_R(results.h_predictions[i]));
            SET_VECTOR_ELT(hp, static_cast<R_xlen_t>(i), vec);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(result_r, 7, hp);
        UNPROTECT(1); // hp
    }

    // 8–10: bootstrap outputs (or NULLs)
    if (results.has_bootstrap_results()) {
        SEXP s = PROTECT(convert_vector_double_to_R(results.bb_predictions));
        SET_VECTOR_ELT(result_r, 8, s); 
        UNPROTECT(1);

        s = PROTECT(convert_vector_double_to_R(results.ci_lower));
        SET_VECTOR_ELT(result_r, 9, s); 
        UNPROTECT(1);

        s = PROTECT(convert_vector_double_to_R(results.ci_upper));
        SET_VECTOR_ELT(result_r, 10, s); 
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result_r, 8,  R_NilValue);
        SET_VECTOR_ELT(result_r, 9,  R_NilValue);
        SET_VECTOR_ELT(result_r, 10, R_NilValue);
    }

    // names while result_r is protected
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
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
        Rf_setAttrib(result_r, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    UNPROTECT(1); // result_r
    return result_r;
}

/**
 * Fixed version of S_upgmalo
 * Fixes: Variable UNPROTECT at lines 615, 628
 * Solution: Use PROTECT_WITH_INDEX/REPROTECT pattern for conditional coercion
 */
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

    // ---- Coerce x/y (and optional y_true) to REAL and copy (using PROTECT_WITH_INDEX) ----
    std::vector<double> x, y, y_true;
    {
        PROTECT_INDEX ipx_x, ipx_y;
        SEXP sx = s_x;
        PROTECT_WITH_INDEX(sx, &ipx_x);
        
        if (TYPEOF(sx) != REALSXP) {
            REPROTECT(sx = Rf_coerceVector(sx, REALSXP), ipx_x);
        }
        const R_xlen_t nx = XLENGTH(sx);
        x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(nx));

        SEXP sy = s_y;
        PROTECT_WITH_INDEX(sy, &ipx_y);
        
        if (TYPEOF(sy) != REALSXP) {
            REPROTECT(sy = Rf_coerceVector(sy, REALSXP), ipx_y);
        }
        const R_xlen_t ny = XLENGTH(sy);
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(ny));

        if (nx != ny) {
            UNPROTECT(2); // sx, sy
            Rf_error("x and y must have the same length.");
        }

        if (s_y_true != R_NilValue) {
            PROTECT_INDEX ipx_yt;
            SEXP syt = s_y_true;
            PROTECT_WITH_INDEX(syt, &ipx_yt);
            
            if (TYPEOF(syt) != REALSXP) {
                REPROTECT(syt = Rf_coerceVector(syt, REALSXP), ipx_yt);
            }
            const R_xlen_t nyt = XLENGTH(syt);
            if (nyt == nx) {
                y_true.assign(REAL(syt), REAL(syt) + static_cast<size_t>(nyt));
            } // else: treat as unavailable (leave empty)
            
            UNPROTECT(1); // syt
        }

        UNPROTECT(2); // sx, sy
    }

    const R_xlen_t n_points = static_cast<R_xlen_t>(x.size());

    // ---- Scalars / parameters (validated) ----
    const bool use_median = (Rf_asLogical(s_use_median) == TRUE);
    const int  h_min = Rf_asInteger(s_h_min);
    const int  h_max = Rf_asInteger(s_h_max);
    const double p   = Rf_asReal(s_p);
    const int  n_bb  = Rf_asInteger(s_n_bb);
    const int  bb_max_distance_deviation = Rf_asInteger(s_bb_max_distance_deviation);
    const int  n_CVs      = Rf_asInteger(s_n_CVs);
    const int  n_CV_folds = Rf_asInteger(s_n_CV_folds);
    const unsigned int seed = static_cast<unsigned int>(Rf_asInteger(s_seed));
    const int    ikernel = Rf_asInteger(s_ikernel);
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const double epsilon = Rf_asReal(s_epsilon);
    const bool   verbose = (Rf_asLogical(s_verbose) == TRUE);

    if (h_min < 1)                Rf_error("h_min must be >= 1.");
    if (h_max < h_min)            Rf_error("h_max must be >= h_min.");
    if (static_cast<R_xlen_t>(h_max) > n_points)
                                  Rf_error("h_max (%d) cannot exceed N (%lld).",
                                           h_max, static_cast<long long>(n_points));
    if (!(p > 0.0 && p <= 1.0))   Rf_error("p must be in (0, 1].");
    if (n_bb < 0)                 Rf_error("n_bb must be >= 0.");
    if (bb_max_distance_deviation < 0)
                                  Rf_error("bb_max_distance_deviation must be >= 0.");
    if (n_CVs < 1)                Rf_error("n_CVs must be >= 1.");
    if (n_CV_folds < 2)           Rf_error("n_CV_folds must be >= 2.");
    if (!(epsilon > 0.0))         Rf_error("epsilon must be > 0.");

    // ---- Core computation (no R allocations inside) ----
    auto cpp_results = upgmalo(x, y, y_true,
                               use_median,
                               h_min, h_max,
                               p,
                               n_bb, bb_max_distance_deviation,
                               n_CVs, n_CV_folds,
                               seed,
                               ikernel,
                               dist_normalization_factor,
                               epsilon,
                               verbose);

    // ---- Build result (container-first; per-element PROTECT/UNPROTECT) ----
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
        SEXP s = PROTECT(Rf_ScalarInteger(cpp_results.opt_h_idx + 1));
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }

    // 3: opt_h
    {
        SEXP s = PROTECT(Rf_ScalarReal(cpp_results.opt_h));
        SET_VECTOR_ELT(result, 3, s);
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

    // 8–10: bootstrap CI pieces (or NULLs)
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

    // 11: true_error (mean of vector, or NULL)
    if (!cpp_results.true_errors.empty()) {
        const double mean_true_error =
            std::accumulate(cpp_results.true_errors.begin(),
                            cpp_results.true_errors.end(), 0.0) /
            static_cast<double>(cpp_results.true_errors.size());
        SEXP s = PROTECT(Rf_ScalarReal(mean_true_error));
        SET_VECTOR_ELT(result, 11, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 11, R_NilValue);
    }

    // 12: h_predictions (list<numeric>)
    {
        SEXP s = PROTECT(convert_vector_vector_double_to_R(cpp_results.h_predictions));
        SET_VECTOR_ELT(result, 12, s);
        UNPROTECT(1);
    }

    // names while result is protected
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0,  Rf_mkChar("h_values"));
        SET_STRING_ELT(names, 1,  Rf_mkChar("h_errors"));
        SET_STRING_ELT(names, 2,  Rf_mkChar("opt_h_idx"));
        SET_STRING_ELT(names, 3,  Rf_mkChar("opt_h"));
        SET_STRING_ELT(names, 4,  Rf_mkChar("graph_adj_list"));
        SET_STRING_ELT(names, 5,  Rf_mkChar("graph_edge_lengths"));
        SET_STRING_ELT(names, 6,  Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 7,  Rf_mkChar("local_predictions"));
        SET_STRING_ELT(names, 8,  Rf_mkChar("bb_predictions"));
        SET_STRING_ELT(names, 9,  Rf_mkChar("ci_lower"));
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