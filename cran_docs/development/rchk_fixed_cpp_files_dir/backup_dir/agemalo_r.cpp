/**
 * @brief Implementation of Adaptive GEMALO algorithm
 */

#include "agemalo.hpp"
#include "graph_maximal_packing.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"

#include <vector>
#include <queue>
#include <chrono>
#include <numeric> // for std::iota()

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_agemalo(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_min_path_size,
        SEXP s_n_packing_vertices,
        SEXP s_max_packing_iterations,
        SEXP s_packing_precision,
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_model_tolerance,
        SEXP s_n_bb,
        SEXP s_cri_probability,
        SEXP s_n_perms,
        SEXP s_blending_coef,
        SEXP s_verbose
        );
}

/**
 * @brief R interface function for the Adaptive Graph Geodesic Model-Averaged LOcal Linear regression algorithm
 *
 * @details This function serves as the bridge between R and the C++ implementation of the Adaptive-GEMALO
 * algorithm. It handles the conversion of R data structures to C++ types, manages memory protection for R
 * objects, and transforms the C++ results back into an R-compatible format. The function is designed to be
 * called from R using .Call() interface.
 *
 * The function performs several key operations:
 * 1. Data Structure Conversion
 *    - Converts R lists representing the graph structure to C++ vectors
 *    - Transforms R numeric vectors to C++ vectors while preserving data integrity
 *    - Handles R's 1-based indexing to C++'s 0-based indexing conversion
 *
 * 2. Memory Management
 *    - Implements proper PROTECT/UNPROTECT mechanisms for R objects
 *    - Ensures memory safety during R-to-C++ and C++-to-R conversions
 *    - Manages temporary storage for intermediate calculations
 *
 * 3. Result Processing
 *    - Creates an R list with named components for return values
 *    - Handles optional bootstrap and permutation test results
 *    - Manages conversion of C++ data structures back to R format
 *
 * @param s_adj_list [SEXP] R list of integer vectors representing graph adjacency lists
 * @param s_weight_list [SEXP] R list of numeric vectors containing edge weights
 * @param s_y [SEXP] R numeric vector of response values at each vertex
 *
 * @param s_min_path_size [SEXP] R integer for minimum required path size
 *
 * @param s_n_packing_vertices [SEXP] R integer specifying number of grid vertices
 * @param s_max_packing_iterations [SEXP] R integer specifying maxima number of maximal packing construction iterations
 * @param s_packing_precision [SEXP] R numeric scaling factor for the precision parameter in the maximal packing construction
 *
 * @param s_n_bws [SEXP] R integer for number of bandwidth values to evaluate
 * @param s_log_grid [SEXP] R logical If true, use logarithmic spacing; if false, use linear spacing
 * @param s_min_bw_factor [SEXP] R numeric scaling factor for minimum bandwidth
 * @param s_max_bw_factor [SEXP] R numeric scaling factor for maximum bandwidth
 *
 * @param s_dist_normalization_factor [SEXP] R numeric factor for distance normalization
 * @param s_kernel_type [SEXP] R integer specifying kernel function type
 *
 * @param s_model_tolerance [SEXP] R numeric convergence tolerance for model fitting
 * @param s_use_mean_error [SEXP] R logical for using mean Rf_error in model averaging
 *
 * @param s_n_bb [SEXP] R integer for number of bootstrap iterations
 * @param s_cri_probability [SEXP] R numeric confidence level for intervals
 *
 * @param s_n_perms [SEXP] R integer for number of permutation test iterations
 *
 * @param s_verbose [SEXP] R logical controlling progress output
 *
 * @return SEXP R list containing:
 *         - graph_diameter: Numeric value of computed graph diameter
 *         - grid_opt_bw: Numeric vector of optimal bandwidths for grid vertices
 *         - predictions: Numeric vector of predictions for original vertices
 *         - grid_predictions: Numeric vector of predictions for grid vertices
 *         - bb_predictions: Matrix of bootstrap predictions (if n_bb > 0)
 *         - cri_lower: Numeric vector of lower confidence bounds (if n_bb > 0)
 *         - cri_upper: Numeric vector of upper confidence bounds (if n_bb > 0)
 *         - null_predictions: Matrix of permutation test predictions (if n_perms > 0)
 *         - p_values: Numeric vector of vertex-wise p-values (if n_perms > 0)
 *         - effect_sizes: Numeric vector of effect sizes (if n_perms > 0)
 *         - significant_vertices: Logical vector indicating significance (if n_perms > 0)
 *
 * @note Implementation Details:
 * - Uses helper functions for converting between R and C++ data structures
 * - Implements careful Rf_error checking for R object types
 * - Maintains proper memory protection counting
 * - Handles matrix transposition for R's column-major format
 *
 * @Rf_warning
 * - All input parameters must be of the correct R type (numeric, integer, logical)
 * - Vertex indices in s_adj_list must be 1-based (R convention)
 * - The function assumes input validation has been performed in R
 * - Large graphs may require significant memory for bootstrap/permutation results
 *
 * @note Memory Management:
 * - Each created R object is protected using PROTECT()
 * - Protection stack is balanced using UNPROTECT() before return
 * - Temporary objects are protected during matrix/vector conversions
 *
 * @see adaptive_uggmalo For the underlying C++ implementation
 * @see convert_R_graph_to_vector For adjacency list conversion details
 * @see convert_R_weights_to_vector For weight list conversion details
 */
SEXP S_agemalo(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    // geodesic parameter
    SEXP s_min_path_size,
    // packing parameters
    SEXP s_n_packing_vertices,
    SEXP s_max_packing_iterations,
    SEXP s_packing_precision,
    // bw parameters
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    // kernel parameters
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    // model parameters
    SEXP s_model_tolerance,
    SEXP s_model_blending_coef,
    // Bayesian bootstrap parameters
    SEXP s_n_bb,
    SEXP s_cri_probability,
    // permutation parameters
    SEXP s_n_perms,
    // verbose
    SEXP s_verbose
) {
    // --- Convert inputs (no PROTECT needed for pure reads) ---
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // y: coerce to REAL safely using indexed protect, then copy out
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

    const size_t n_vertices = adj_list.size();
    if (y.size() != n_vertices) {
        Rf_error("Length of 'y' (%zu) must match number of vertices (%zu).",
                 y.size(), n_vertices);
    }

    // --- Scalars (use Rf_as* for defensive coercion) ---
    // geodesic
    const size_t min_path_size = static_cast<size_t>(Rf_asInteger(s_min_path_size));

    // packing
    const size_t n_packing_vertices   = static_cast<size_t>(Rf_asInteger(s_n_packing_vertices));
    const size_t max_packing_iterations = static_cast<size_t>(Rf_asInteger(s_max_packing_iterations));
    const double packing_precision     = Rf_asReal(s_packing_precision);

    // bandwidth grid
    const size_t n_bws        = static_cast<size_t>(Rf_asInteger(s_n_bws));
    const bool   log_grid     = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);

    // kernel
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type = static_cast<size_t>(Rf_asInteger(s_kernel_type));

    // model
    const double model_tolerance     = Rf_asReal(s_model_tolerance);
    const double model_blending_coef = Rf_asReal(s_model_blending_coef);

    // Bayesian bootstrap
    const size_t n_bb           = static_cast<size_t>(Rf_asInteger(s_n_bb));
    const double cri_probability = Rf_asReal(s_cri_probability);

    // permutations
    const size_t n_perms = static_cast<size_t>(Rf_asInteger(s_n_perms));

    // verbosity
    const bool verbose = (Rf_asLogical(s_verbose) == TRUE);

    // --- Build (or reuse) a uniform grid graph ---
    uniform_grid_graph_t grid_graph;
    if (n_packing_vertices < n_vertices) {
        grid_graph = create_maximal_packing(adj_list, weight_list,
                                            n_packing_vertices,
                                            max_packing_iterations,
                                            packing_precision);
    } else {
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), 0);
        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);
        // approximate diameter from two BFS sweeps
        auto e1 = grid_graph.get_vertex_eccentricity(0);
        auto e2 = grid_graph.get_vertex_eccentricity(e1.first);
        grid_graph.graph_diameter = e2.second;
    }

    // --- Core computation (no R allocation inside) ---
    agemalo_result_t res = agemalo(
        grid_graph,
        y,
        // geodesic
        min_path_size,
        // bw
        n_bws, log_grid, min_bw_factor, max_bw_factor,
        // kernel
        dist_normalization_factor, kernel_type,
        // model
        model_tolerance, model_blending_coef,
        // bayesian bootstrap
        n_bb, cri_probability,
        // permutations
        n_perms,
        // verbose
        verbose
    );

    // --- Assemble return list ---
    // Names (fixed, constant count)
    const char* names[] = {
        "graph_diameter",            // 0
        "packing_radius",            // 1
        "packing_vertices",          // 2
        "grid_opt_bw",               // 3
        "predictions",               // 4
        "errors",                    // 5
        "scale",                     // 6
        "grid_predictions",          // 7
        "bb_predictions",            // 8
        "cri_lower",                 // 9
        "cri_upper",                 // 10
        "null_predictions",          // 11
        "null_predictions_cri_lower",// 12
        "null_predictions_cri_upper",// 13
        "p_values",                  // 14
        "effect_sizes",              // 15
        "significant_vertices",      // 16
        "null_distances"             // 17
    };
    const int N = 18;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));

    // 0: graph_diameter
    {
        SEXP s = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s)[0] = grid_graph.graph_diameter;
        SET_VECTOR_ELT(result, 0, s);
        UNPROTECT(1);
    }
    // 1: packing_radius
    {
        SEXP s = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s)[0] = grid_graph.max_packing_radius;
        SET_VECTOR_ELT(result, 1, s);
        UNPROTECT(1);
    }
    // 2: packing_vertices (1-based)
    {
        const R_xlen_t m = static_cast<R_xlen_t>(grid_graph.grid_vertices.size());
        SEXP s = PROTECT(Rf_allocVector(INTSXP, m));
        int* p = INTEGER(s);
        for (R_xlen_t i = 0; i < m; ++i) p[i] = grid_graph.grid_vertices[static_cast<size_t>(i)] + 1;
        SET_VECTOR_ELT(result, 2, s);
        UNPROTECT(1);
    }
    // 3: grid_opt_bw
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.grid_opt_bw.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.grid_opt_bw.begin(), res.grid_opt_bw.end(), REAL(s));
        SET_VECTOR_ELT(result, 3, s);
        UNPROTECT(1);
    }
    // 4: predictions
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.predictions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.predictions.begin(), res.predictions.end(), REAL(s));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }
    // 5: errors
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.errors.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.errors.begin(), res.errors.end(), REAL(s));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }
    // 6: scale
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.scale.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.scale.begin(), res.scale.end(), REAL(s));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }
    // 7: grid_predictions (assuming map_to_vector already applied in core; if not, adapt)
    {
        const R_xlen_t n = static_cast<R_xlen_t>(res.grid_predictions.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.grid_predictions.begin(), res.grid_predictions.end(), REAL(s));
        SET_VECTOR_ELT(result, 7, s);
        UNPROTECT(1);
    }
    // 8: bb_predictions (matrix or NULL)
    if (!res.bb_predictions.empty()) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.bb_predictions.size());
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.bb_predictions.front().size());
        SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(s);
        for (R_xlen_t i = 0; i < nrow; ++i)
            for (R_xlen_t j = 0; j < ncol; ++j)
                p[i + j * nrow] = res.bb_predictions[static_cast<size_t>(i)][static_cast<size_t>(j)];
        SET_VECTOR_ELT(result, 8, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 8, R_NilValue);
    }
    // 9: cri_lower
    if (!res.cri_lower.empty()) {
        const R_xlen_t n = static_cast<R_xlen_t>(res.cri_lower.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.cri_lower.begin(), res.cri_lower.end(), REAL(s));
        SET_VECTOR_ELT(result, 9, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 9, R_NilValue);
    }
    // 10: cri_upper
    if (!res.cri_upper.empty()) {
        const R_xlen_t n = static_cast<R_xlen_t>(res.cri_upper.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.cri_upper.begin(), res.cri_upper.end(), REAL(s));
        SET_VECTOR_ELT(result, 10, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }
    // 11: null_predictions
    if (!res.null_predictions.empty()) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.null_predictions.size());
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.null_predictions.front().size());
        SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(s);
        for (R_xlen_t i = 0; i < nrow; ++i)
            for (R_xlen_t j = 0; j < ncol; ++j)
                p[i + j * nrow] = res.null_predictions[static_cast<size_t>(i)][static_cast<size_t>(j)];
        SET_VECTOR_ELT(result, 11, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 11, R_NilValue);
    }
    // 12: null_predictions_cri_lower
    if (!res.null_predictions_cri_lower.empty()) {
        const R_xlen_t n = static_cast<R_xlen_t>(res.null_predictions_cri_lower.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.null_predictions_cri_lower.begin(), res.null_predictions_cri_lower.end(), REAL(s));
        SET_VECTOR_ELT(result, 12, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 12, R_NilValue);
    }
    // 13: null_predictions_cri_upper
    if (!res.null_predictions_cri_upper.empty()) {
        const R_xlen_t n = static_cast<R_xlen_t>(res.null_predictions_cri_upper.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(res.null_predictions_cri_upper.begin(), res.null_predictions_cri_upper.end(), REAL(s));
        SET_VECTOR_ELT(result, 13, s);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 13, R_NilValue);
    }
    // 14..17: permutation test outputs (or NULL)
    if (res.permutation_tests) {
        // 14: p_values
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->p_values.size());
            SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
            std::copy(res.permutation_tests->p_values.begin(), res.permutation_tests->p_values.end(), REAL(s));
            SET_VECTOR_ELT(result, 14, s);
            UNPROTECT(1);
        }
        // 15: effect_sizes
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->effect_sizes.size());
            SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
            std::copy(res.permutation_tests->effect_sizes.begin(), res.permutation_tests->effect_sizes.end(), REAL(s));
            SET_VECTOR_ELT(result, 15, s);
            UNPROTECT(1);
        }
        // 16: significant_vertices (logical)
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->significant_vertices.size());
            SEXP s = PROTECT(Rf_allocVector(LGLSXP, n));
            int* p = LOGICAL(s);
            for (R_xlen_t i = 0; i < n; ++i)
                p[i] = res.permutation_tests->significant_vertices[static_cast<size_t>(i)] ? 1 : 0;
            SET_VECTOR_ELT(result, 16, s);
            UNPROTECT(1);
        }
        // 17: null_distances
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->null_distances.size());
            SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
            std::copy(res.permutation_tests->null_distances.begin(),
                      res.permutation_tests->null_distances.end(), REAL(s));
            SET_VECTOR_ELT(result, 17, s);
            UNPROTECT(1);
        }
    } else {
        for (int i = 14; i <= 17; ++i) SET_VECTOR_ELT(result, i, R_NilValue);
    }

    // names
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, N));
    for (int i = 0; i < N; ++i) SET_STRING_ELT(nms, i, Rf_mkChar(names[i]));
    Rf_setAttrib(result, R_NamesSymbol, nms);

    UNPROTECT(2); // result, nms
    return result;
}
