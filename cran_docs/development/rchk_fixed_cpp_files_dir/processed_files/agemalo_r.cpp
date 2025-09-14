/**
 * @brief Fixed version of S_agemalo to address rchk PROTECT/UNPROTECT issues
 * 
 * Changes made:
 * 1. Removed lambdas that capture protect_count by reference
 * 2. Fixed PROTECT/UNPROTECT balance using container-first pattern
 * 3. Used PROTECT_WITH_INDEX/REPROTECT for conditional coercion
 * 4. Ensured constant literal UNPROTECT at function end
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
        SEXP s_model_blending_coef,
        SEXP s_n_bb,
        SEXP s_cri_probability,
        SEXP s_n_perms,
        SEXP s_verbose
        );
}

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

    // --- Assemble return list (container-first pattern) ---
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
    // 7: grid_predictions
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
    
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, N));
    for (int i = 0; i < N; ++i) SET_STRING_ELT(nms, i, Rf_mkChar(names[i]));
    Rf_setAttrib(result, R_NamesSymbol, nms);

    UNPROTECT(2); // result, nms
    return result;
}


/**
 * @brief Fixed version of S_amagelo to address rchk PROTECT/UNPROTECT issues
 * 
 * Changes made:
 * 1. Removed lambdas that capture protect_count by reference
 * 2. Fixed UNPROTECT(variable) issue at line 303 - now uses constant 2
 * 3. Properly balanced PROTECT/UNPROTECT for all allocations
 * 4. Used container-first pattern consistently
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
