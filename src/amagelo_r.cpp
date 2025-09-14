#include "amagelo.hpp"
#include "uniform_grid_graph.hpp" // For create_chain_graph, etc.

#include <vector>
#include <unordered_set>
#include <numeric>                // For std::iota()

#include <R.h>
#include <Rinternals.h>

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
    // --- 1) Unmarshal inputs (coercion block) ---
    std::vector<double> x, y;
    {
        SEXP sx = s_x, sy = s_y;
        PROTECT_INDEX px, py;
        PROTECT_WITH_INDEX(sx, &px);
        PROTECT_WITH_INDEX(sy, &py);

        if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), px);
        if (TYPEOF(sy) != REALSXP) REPROTECT(sy = Rf_coerceVector(sy, REALSXP), py);

        R_xlen_t n = XLENGTH(sx);
        if (XLENGTH(sy) != n) {
            UNPROTECT(2);
            Rf_error("length(x) must equal length(y)");
        }

        x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(n));
        y.assign(REAL(sy), REAL(sy) + static_cast<size_t>(n));
        UNPROTECT(2); // sx, sy
    }

    // Extract scalar parameters using defensive coercion
    size_t grid_size               = static_cast<size_t>(Rf_asInteger(s_grid_size));
    double min_bw_factor           = Rf_asReal(s_min_bw_factor);
    double max_bw_factor           = Rf_asReal(s_max_bw_factor);
    size_t n_bws                   = static_cast<size_t>(Rf_asInteger(s_n_bws));
    bool   use_global_bw_grid      = (Rf_asLogical(s_use_global_bw_grid) == TRUE);
    bool   with_bw_predictions     = (Rf_asLogical(s_with_bw_predictions) == TRUE);
    bool   log_grid                = (Rf_asLogical(s_log_grid) == TRUE);
    size_t domain_min_size         = static_cast<size_t>(Rf_asInteger(s_domain_min_size));
    size_t kernel_type             = static_cast<size_t>(Rf_asInteger(s_kernel_type));
    double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    size_t n_cleveland_iterations  = static_cast<size_t>(Rf_asInteger(s_n_cleveland_iterations));
    double blending_coef           = Rf_asReal(s_blending_coef);
    bool   use_linear_blending     = (Rf_asLogical(s_use_linear_blending) == TRUE);
    double precision               = Rf_asReal(s_precision);
    double small_depth_threshold   = Rf_asReal(s_small_depth_threshold);
    double depth_similarity_tol    = Rf_asReal(s_depth_similarity_tol);
    bool   verbose                 = (Rf_asLogical(s_verbose) == TRUE);

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

    // --- 3) Prepare R return list (container-first pattern) ---
    const int n_el = 17; // Fixed number of elements
    SEXP r_list = PROTECT(Rf_allocVector(VECSXP, n_el));

    int idx = 0;

    // x_sorted
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.x_sorted.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.x_sorted.begin(), result.x_sorted.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // y_sorted
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.y_sorted.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.y_sorted.begin(), result.y_sorted.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // order (1-based)
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.order.size());
        SEXP s = PROTECT(Rf_allocVector(INTSXP, n));
        int* p = INTEGER(s);
        for (R_xlen_t i = 0; i < n; ++i)
            p[i] = static_cast<int>(result.order[static_cast<size_t>(i)] + 1);
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // grid_coords
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.grid_coords.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.grid_coords.begin(), result.grid_coords.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // predictions
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.predictions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.predictions.begin(), result.predictions.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // bw_predictions (matrix or NULL)
    if (with_bw_predictions && !result.bw_predictions.empty()) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(result.predictions.size());
        const R_xlen_t ncol = static_cast<R_xlen_t>(result.bw_predictions.size());
        for (R_xlen_t j = 0; j < ncol; ++j) {
            if (static_cast<R_xlen_t>(result.bw_predictions[static_cast<size_t>(j)].size()) != nrow) {
                Rf_error("Inconsistent matrix column length");
            }
        }

        // Allocate and fill after validation
        SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(s);
        for (R_xlen_t j = 0; j < ncol; ++j)
            for (R_xlen_t i = 0; i < nrow; ++i)
                p[i + j * nrow] = result.bw_predictions[static_cast<size_t>(j)][static_cast<size_t>(i)];
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(r_list, idx++, R_NilValue);
    }

    // grid_predictions
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.grid_predictions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.grid_predictions.begin(), result.grid_predictions.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // harmonic_predictions
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.harmonic_predictions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.harmonic_predictions.begin(), result.harmonic_predictions.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // local_extrema matrix
    {
        const auto& extrema = result.local_extrema;
        const R_xlen_t nrow = static_cast<R_xlen_t>(extrema.size());
        const R_xlen_t ncol = 8;

        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(mat);

        for (R_xlen_t row = 0; row < nrow; ++row) {
            const auto& e = extrema[static_cast<size_t>(row)];
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
        SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol));
        const char* colnames_cstr[] = {
            "idx", "x", "y", "is_max", "depth", "depth_idx",
            "rel_depth", "range_rel_depth"
        };
        for (R_xlen_t j = 0; j < ncol; ++j)
            SET_STRING_ELT(colnames, j, Rf_mkChar(colnames_cstr[j]));

        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);  // rownames
        SET_VECTOR_ELT(dimnames, 1, colnames);    // colnames
        Rf_setAttrib(mat, R_DimNamesSymbol, dimnames);

        SET_VECTOR_ELT(r_list, idx++, mat);
        UNPROTECT(3); // mat, colnames, dimnames
    }

    // harmonic_predictions_local_extrema matrix
    {
        const auto& extrema = result.harmonic_predictions_local_extrema;
        const R_xlen_t nrow = static_cast<R_xlen_t>(extrema.size());
        const R_xlen_t ncol = 8;

        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(mat);

        for (R_xlen_t row = 0; row < nrow; ++row) {
            const auto& e = extrema[static_cast<size_t>(row)];
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
        SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol));
        const char* colnames_cstr[] = {
            "idx", "x", "y", "is_max", "depth", "depth_idx",
            "rel_depth", "range_rel_depth"
        };
        for (R_xlen_t j = 0; j < ncol; ++j)
            SET_STRING_ELT(colnames, j, Rf_mkChar(colnames_cstr[j]));

        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);  // rownames
        SET_VECTOR_ELT(dimnames, 1, colnames);    // colnames
        Rf_setAttrib(mat, R_DimNamesSymbol, dimnames);

        SET_VECTOR_ELT(r_list, idx++, mat);
        UNPROTECT(3); // mat, colnames, dimnames
    }

    // monotonic_interval_proportions
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.monotonic_interval_proportions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.monotonic_interval_proportions.begin(),
                  result.monotonic_interval_proportions.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // change_scaled_monotonicity_index
    {
        SEXP s = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s)[0] = result.change_scaled_monotonicity_index;
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // bw_errors
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.bw_errors.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.bw_errors.begin(), result.bw_errors.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // opt_bw_idx (1-based)
    {
        SEXP s = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(s)[0] = static_cast<int>(result.opt_bw_idx + 1);
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // min_bw
    {
        SEXP s = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s)[0] = result.min_bw;
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // max_bw
    {
        SEXP s = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s)[0] = result.max_bw;
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // bws
    {
        const R_xlen_t n = static_cast<R_xlen_t>(result.bws.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(result.bws.begin(), result.bws.end(), REAL(s));
        SET_VECTOR_ELT(r_list, idx++, s);
        UNPROTECT(1);
    }

    // Set names
    const char* names[] = {
        "x_sorted",
        "y_sorted",
        "order",
        "grid_coords",
        "predictions",
        "bw_predictions",
        "grid_predictions",
        "harmonic_predictions",
        "local_extrema",
        "harmonic_predictions_local_extrema",
        "monotonic_interval_proportions",
        "change_scaled_monotonicity_index",
        "bw_errors",
        "opt_bw_idx",
        "min_bw",
        "max_bw",
        "bws"
    };

    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_el));
    for (int i = 0; i < n_el; ++i)
        SET_STRING_ELT(r_names, i, Rf_mkChar(names[i]));
    Rf_setAttrib(r_list, R_NamesSymbol, r_names);

    UNPROTECT(2); // r_list, r_names
    return r_list;
}
