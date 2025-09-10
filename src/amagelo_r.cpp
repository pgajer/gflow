#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros from R headers
#undef length
#undef Rf_eval

#include <vector>
#include <unordered_set>
#include <numeric>                // For std::iota()

#include "amagelo.hpp"
#include "uniform_grid_graph.hpp" // For create_chain_graph, etc.

/**
 * @brief R interface to the AMAGELO algorithm
 *
 * This function serves as the interface between R and the C++ implementation
 * of AMAGELO. It handles conversion between R and C++ data structures, input
 * validation, and packaging results into an R list.
 *
 * @param s_x                     R numeric vector of predictor values
 * @param s_y                     R numeric vector of response values
 * @param s_grid_size             R integer for grid size
 * @param s_min_bw_factor         R numeric for minimum bandwidth factor
 * @param s_max_bw_factor         R numeric for maximum bandwidth factor
 * @param s_n_bws                 R integer for number of bandwidths to try
 * @param s_use_global_bw_grid    R logical for global bandwidth grid usage
 * @param s_with_bw_predictions   R logical for storing per-bandwidth predictions
 * @param s_log_grid              R logical for logarithmic bandwidth grid
 * @param s_domain_min_size       R integer for minimum domain size
 * @param s_kernel_type           R integer code for kernel function
 * @param s_dist_normalization_factor R numeric scale factor for distances
 * @param s_n_cleveland_iterations R integer for robustness iterations
 * @param s_blending_coef         R numeric for model averaging weight
 * @param s_use_linear_blending   R logical for linear vs. power blending
 * @param s_precision             R numeric for optimization precision
 * @param s_small_depth_threshold R numeric for wiggle detection threshold
 * @param s_depth_similarity_tol  R numeric for depth similarity tolerance
 * @param s_verbose               R logical for verbose output
 *
 * @return An R list containing the following components:
 *   - x_sorted, y_sorted: Sorted data vectors
 *   - order: 1-based indices of original data order
 *   - grid_coords: Grid point x-coordinates
 *   - predictions: Smoothed values at optimal bandwidth
 *   - bw_predictions: Matrix of predictions at different bandwidths (if requested)
 *   - grid_predictions: Predictions at grid points
 *   - harmonic_predictions: Predictions after triplet harmonic smoothing
 *   - local_extrema: Matrix of local extrema information
 *   - monotonic_interval_proportions: Proportions of monotonic intervals
 *   - change_scaled_monotonicity_index: Weighted signed average of directional changes, quantifying monotonicity strength and directionality; values close to \(+1\) or \(-1\) indicate strong global monotonic trends.
 *   - bw_errors: Cross-validation errors per bandwidth
 *   - opt_bw_idx: 1-based index of optimal bandwidth
 *   - min_bw, max_bw: Bandwidth range values
 *   - bws: Vector of evaluated bandwidths
 */
extern "C" SEXP S_amagelo(
    SEXP s_x,
    SEXP s_y,
    SEXP s_grid_size,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_n_bws,
    SEXP s_use_global_bw_grid,
    SEXP s_with_bw_predictions,
    SEXP s_log_grid,
    SEXP s_domain_min_size,
    SEXP s_kernel_type,
    SEXP s_dist_normalization_factor,
    SEXP s_n_cleveland_iterations,
    SEXP s_blending_coef,
    SEXP s_use_linear_blending,
    SEXP s_precision,
    SEXP s_small_depth_threshold,
    SEXP s_depth_similarity_tol,
    SEXP s_verbose
) {
    // --- 1) Unmarshal inputs ---
    if (TYPEOF(s_x) != REALSXP || TYPEOF(s_y) != REALSXP)
        Rf_error("x and y must be real vectors");

    R_xlen_t n = XLENGTH(s_x);
    if (XLENGTH(s_y) != n)
        Rf_error("Rf_length(x) must equal Rf_length(y)");

    double* x_ptr = REAL(s_x);
    double* y_ptr = REAL(s_y);
    std::vector<double> x(x_ptr, x_ptr + n);
    std::vector<double> y(y_ptr, y_ptr + n);

    size_t grid_size               = (size_t) INTEGER(s_grid_size)[0];
    double min_bw_factor           = REAL(s_min_bw_factor)[0];
    double max_bw_factor           = REAL(s_max_bw_factor)[0];
    size_t n_bws                   = (size_t) INTEGER(s_n_bws)[0];
    bool   use_global_bw_grid      = LOGICAL(s_use_global_bw_grid)[0];
    bool   with_bw_predictions     = LOGICAL(s_with_bw_predictions)[0];
    bool   log_grid                = LOGICAL(s_log_grid)[0];
    size_t domain_min_size         = (size_t) INTEGER(s_domain_min_size)[0];
    size_t kernel_type             = (size_t) INTEGER(s_kernel_type)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    size_t n_cleveland_iterations  = (size_t) INTEGER(s_n_cleveland_iterations)[0];
    double blending_coef           = REAL(s_blending_coef)[0];
    bool   use_linear_blending     = LOGICAL(s_use_linear_blending)[0];
    double precision               = REAL(s_precision)[0];
    double small_depth_threshold   = REAL(s_small_depth_threshold)[0];
    double depth_similarity_tol    = REAL(s_depth_similarity_tol)[0];
    bool   verbose                 = LOGICAL(s_verbose)[0];

    // --- 2) Call C++ backend ---
    amagelo_t result = amagelo(
        x, y,
        grid_size,
        min_bw_factor,
        max_bw_factor,
        n_bws,
        use_global_bw_grid,
        with_bw_predictions,
        log_grid,
        domain_min_size,
        kernel_type,
        dist_normalization_factor,
        n_cleveland_iterations,
        blending_coef,
        use_linear_blending,
        precision,
        small_depth_threshold,
        depth_similarity_tol,
        verbose
    );

    // --- 3) Prepare R return list ---
    const char* names[] = {
        "x_sorted",
        "y_sorted",
        "order",
        "grid_coords",
        // predictions
        "predictions",
        "bw_predictions",
        "grid_predictions",
        "harmonic_predictions",
        // local extrema
        "local_extrema",
        "harmonic_predictions_local_extrema",
        // monotonicity measures
        // "tvmi",
        "monotonic_interval_proportions",
        "change_scaled_monotonicity_index",
        // "simpson_index",
        // errors etc
        "bw_errors",
        "opt_bw_idx",
        "min_bw",
        "max_bw",
        "bws",
        NULL
    };
    int n_el = 0;
    while (names[n_el]) ++n_el;

    int protect_count = 0;
    SEXP r_list  = PROTECT(Rf_allocVector(VECSXP, n_el)); ++protect_count;
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_el)); ++protect_count;
    for (int i = 0; i < n_el; ++i)
        SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
    Rf_setAttrib(r_list, R_NamesSymbol, r_names);

    // Helpers
    auto vec_real = [&](const std::vector<double>& v) {
        SEXP ans = PROTECT(Rf_allocVector(REALSXP, v.size())); ++protect_count;
        std::copy(v.begin(), v.end(), REAL(ans));
        return ans;
    };
    auto vec_int = [&](const std::vector<size_t>& v) {
        SEXP ans = PROTECT(Rf_allocVector(INTSXP, v.size())); ++protect_count;
        for (size_t i = 0; i < v.size(); ++i)
            INTEGER(ans)[i] = (int)(v[i] + 1);  // 1‐based
        return ans;
    };
    auto mat_real = [&](const std::vector<std::vector<double>>& M, R_xlen_t nrow) {
        R_xlen_t ncol = M.size();
        SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); ++protect_count;
        double* p = REAL(ans);
        for (R_xlen_t j = 0; j < ncol; ++j) {
            if ((R_xlen_t)M[j].size() != nrow)
                Rf_error("Inconsistent matrix column length");
            for (R_xlen_t i = 0; i < nrow; ++i)
                p[i + j * nrow] = M[j][i];
        }
        return ans;
    };

    // Fill in
    int i = 0;
    SET_VECTOR_ELT(r_list, i++, vec_real(result.x_sorted));
    SET_VECTOR_ELT(r_list, i++, vec_real(result.y_sorted));
    SET_VECTOR_ELT(r_list, i++, vec_int (result.order));
    SET_VECTOR_ELT(r_list, i++, vec_real(result.grid_coords));

    SET_VECTOR_ELT(r_list, i++, vec_real(result.predictions));
    if (with_bw_predictions) {
        SET_VECTOR_ELT(r_list, i++, mat_real(result.bw_predictions,   (R_xlen_t)result.predictions.size()));
    } else {
        SET_VECTOR_ELT(r_list, i++, R_NilValue);
    }
    SET_VECTOR_ELT(r_list, i++, vec_real(result.grid_predictions));
    SET_VECTOR_ELT(r_list, i++, vec_real(result.harmonic_predictions));

    // local_extrema matrix: columns = [idx, x, y, is_max, depth, depth_idx, rel_depth, range_rel_depth]
    {
        const auto& extrema = result.local_extrema;
        R_xlen_t nrow = extrema.size();
        R_xlen_t ncol = 8;

        // Allocate matrix
        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); ++protect_count;
        double* p = REAL(mat);

        // Fill matrix by sorted order
        for (R_xlen_t row = 0; row < nrow; ++row) {
            const auto& e = extrema[row];
            p[row + 0 * nrow] = static_cast<double>(e.idx + 1);
            p[row + 1 * nrow] = e.x;
            p[row + 2 * nrow] = e.y;
            p[row + 3 * nrow] = e.is_max ? 1.0 : 0.0;
            p[row + 4 * nrow] = e.depth;
            p[row + 5 * nrow] = static_cast<double>(e.depth_idx + 1);
            p[row + 6 * nrow] = e.rel_depth;
            p[row + 7 * nrow] = e.range_rel_depth;
        }

        // Add column names
        SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol)); ++protect_count;
        const char* colnames_cstr[] = {
            "idx", "x", "y", "is_max", "depth", "depth_idx",
            "rel_depth", "range_rel_depth"
        };
        for (R_xlen_t j = 0; j < ncol; ++j)
            SET_STRING_ELT(colnames, j, Rf_mkChar(colnames_cstr[j]));

        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2)); ++protect_count;
        SET_VECTOR_ELT(dimnames, 1, colnames);  // colnames
        Rf_setAttrib(mat, R_DimNamesSymbol, dimnames);

        SET_VECTOR_ELT(r_list, i++, mat);  // local_extrema = 8th element in names[]
    }

    // harmonic_predictions_local_extrema matrix: columns = [idx, x, y, is_max, depth, depth_idx, rel_depth, range_rel_depth]
    {
        const auto& extrema = result.harmonic_predictions_local_extrema;
        R_xlen_t nrow = extrema.size();
        R_xlen_t ncol = 8;

        // Allocate matrix
        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); ++protect_count;
        double* p = REAL(mat);

        // Fill matrix by sorted order
        for (R_xlen_t row = 0; row < nrow; ++row) {
            const auto& e = extrema[row];
            p[row + 0 * nrow] = static_cast<double>(e.idx + 1);
            p[row + 1 * nrow] = e.x;
            p[row + 2 * nrow] = e.y;
            p[row + 3 * nrow] = e.is_max ? 1.0 : 0.0;
            p[row + 4 * nrow] = e.depth;
            p[row + 5 * nrow] = static_cast<double>(e.depth_idx + 1);
            p[row + 6 * nrow] = e.rel_depth;
            p[row + 7 * nrow] = e.range_rel_depth;
        }

        // Add column names
        SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol)); ++protect_count;
        const char* colnames_cstr[] = {
            "idx", "x", "y", "is_max", "depth", "depth_idx",
            "rel_depth", "range_rel_depth"
        };
        for (R_xlen_t j = 0; j < ncol; ++j)
            SET_STRING_ELT(colnames, j, Rf_mkChar(colnames_cstr[j]));

        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2)); ++protect_count;
        SET_VECTOR_ELT(dimnames, 1, colnames);  // colnames
        Rf_setAttrib(mat, R_DimNamesSymbol, dimnames);

        SET_VECTOR_ELT(r_list, i++, mat);  // local_extrema = 8th element in names[]
    }

    // SET_VECTOR_ELT(r_list, i++, vec_real(std::vector<double>{ result.tvmi }));
    SET_VECTOR_ELT(r_list, i++, vec_real(result.monotonic_interval_proportions));
    // SET_VECTOR_ELT(r_list, i++, vec_real(std::vector<double>{ result.simpson_index }));
    SET_VECTOR_ELT(r_list, i++, vec_real(std::vector<double>{ result.change_scaled_monotonicity_index }));

    SET_VECTOR_ELT(r_list, i++, vec_real(result.bw_errors));

    // opt_bw_idx
    SEXP sx = PROTECT(Rf_allocVector(INTSXP, 1)); ++protect_count;
    INTEGER(sx)[0] = (int)(result.opt_bw_idx + 1); // 1-based
    SET_VECTOR_ELT(r_list, i++, sx);

    // min_bw, max_bw
    SET_VECTOR_ELT(r_list, i++, vec_real(std::vector<double>{ result.min_bw }));
    SET_VECTOR_ELT(r_list, i++, vec_real(std::vector<double>{ result.max_bw }));

    // bws
    SET_VECTOR_ELT(r_list, i++, vec_real(result.bws));

    UNPROTECT(protect_count);
    return r_list;
}
