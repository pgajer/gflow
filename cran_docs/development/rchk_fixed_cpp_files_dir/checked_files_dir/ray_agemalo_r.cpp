/**
 * @brief Fixed version of S_ray_agemalo to address rchk PROTECT/UNPROTECT issues
 * 
 * Changes made:
 * 1. Removed lambdas that capture and modify protection counters
 * 2. Fixed incorrect protection counting (lines 384, 386)
 * 3. Used coercion block with PROTECT_WITH_INDEX/REPROTECT for y
 * 4. Properly balanced all PROTECT/UNPROTECT pairs
 * 5. Used container-first pattern consistently
 * 6. Ensured constant literal UNPROTECT at function end
 */

#include "ray_agemalo.hpp"
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
    SEXP S_ray_agemalo(
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

SEXP S_ray_agemalo(
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
    // --- Convert inputs (no PROTECT needed for pure reads from helper functions) ---
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // --- Coerce y to REAL safely using indexed protect, then copy out ---
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

    // --- Extract scalar parameters using defensive coercion ---
    // geodesic parameter
    const size_t min_path_size = static_cast<size_t>(Rf_asInteger(s_min_path_size));
    // packing parameters
    const size_t n_packing_vertices = static_cast<size_t>(Rf_asInteger(s_n_packing_vertices));
    const size_t max_packing_iterations = static_cast<size_t>(Rf_asInteger(s_max_packing_iterations));
    const double packing_precision = Rf_asReal(s_packing_precision);

    // bw parameters
    const size_t n_bws = static_cast<size_t>(Rf_asInteger(s_n_bws));
    const bool log_grid = (Rf_asLogical(s_log_grid) == TRUE);
    const double min_bw_factor = Rf_asReal(s_min_bw_factor);
    const double max_bw_factor = Rf_asReal(s_max_bw_factor);

    // kernel parameters
    const double dist_normalization_factor = Rf_asReal(s_dist_normalization_factor);
    const size_t kernel_type = static_cast<size_t>(Rf_asInteger(s_kernel_type));
    // model parameters
    const double model_tolerance = Rf_asReal(s_model_tolerance);
    const double model_blending_coef = Rf_asReal(s_model_blending_coef);
    // Bayesian bootstrap parameters
    const size_t n_bb = static_cast<size_t>(Rf_asInteger(s_n_bb));
    const double cri_probability = Rf_asReal(s_cri_probability);
    // permutation parameters
    const size_t n_perms = static_cast<size_t>(Rf_asInteger(s_n_perms));
    // verbose
    const bool verbose = (Rf_asLogical(s_verbose) == TRUE);

    // --- Build grid graph ---
    uniform_grid_graph_t grid_graph;
    const size_t n_vertices = adj_list.size();
    if (n_packing_vertices < n_vertices) {
        grid_graph = create_maximal_packing(adj_list,
                                            weight_list,
                                            n_packing_vertices,
                                            max_packing_iterations,
                                            packing_precision);
    } else {
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), 0);
        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);
        
        // Find diameter endpoints
        auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);
        auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
        grid_graph.graph_diameter = diameter;
    }

    // --- Call the C++ function (no R allocations inside) ---
    agemalo_result_t res = ray_agemalo(
        grid_graph,
        y,
        // geodesic parameter
        min_path_size,
        // bw parameters
        n_bws,
        log_grid,
        min_bw_factor,
        max_bw_factor,
        // kernel parameters
        dist_normalization_factor,
        kernel_type,
        // model parameters
        model_tolerance,
        model_blending_coef,
        // Bayesian bootstrap parameters
        n_bb,
        cri_probability,
        // permutation parameters
        n_perms,
        // verbose
        verbose
    );

    // --- Create the return list (container-first pattern) ---
    const int n_elements = 18;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));

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

    // 2: packing_vertices (1-based indices)
    {
        const R_xlen_t n = static_cast<R_xlen_t>(grid_graph.grid_vertices.size());
        SEXP s = PROTECT(Rf_allocVector(INTSXP, n));
        int* p = INTEGER(s);
        for (R_xlen_t i = 0; i < n; ++i) {
            p[i] = static_cast<int>(grid_graph.grid_vertices[static_cast<size_t>(i)] + 1);
        }
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

    // 7: grid_predictions (convert map to vector if needed)
    {
        std::vector<double> grid_pred_vec = map_to_vector(res.grid_predictions_map);
        const R_xlen_t n = static_cast<R_xlen_t>(grid_pred_vec.size());
        SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
        std::copy(grid_pred_vec.begin(), grid_pred_vec.end(), REAL(s));
        SET_VECTOR_ELT(result, 7, s);
        UNPROTECT(1);
    }

    // 8: bb_predictions (matrix or NULL)
    if (!res.bb_predictions.empty()) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.bb_predictions.size());
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.bb_predictions[0].size());
        SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(s);
        for (R_xlen_t i = 0; i < nrow; ++i) {
            for (R_xlen_t j = 0; j < ncol; ++j) {
                p[i + j * nrow] = res.bb_predictions[static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
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

    // 11: null_predictions (matrix or NULL)
    if (!res.null_predictions.empty()) {
        const R_xlen_t nrow = static_cast<R_xlen_t>(res.null_predictions.size());
        const R_xlen_t ncol = static_cast<R_xlen_t>(res.null_predictions[0].size());
        SEXP s = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* p = REAL(s);
        for (R_xlen_t i = 0; i < nrow; ++i) {
            for (R_xlen_t j = 0; j < ncol; ++j) {
                p[i + j * nrow] = res.null_predictions[static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
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

    // 14-17: permutation test results if available
    if (res.permutation_tests) {
        // 14: p_values
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->p_values.size());
            SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
            std::copy(res.permutation_tests->p_values.begin(), 
                      res.permutation_tests->p_values.end(), REAL(s));
            SET_VECTOR_ELT(result, 14, s);
            UNPROTECT(1);
        }
        
        // 15: effect_sizes
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->effect_sizes.size());
            SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
            std::copy(res.permutation_tests->effect_sizes.begin(), 
                      res.permutation_tests->effect_sizes.end(), REAL(s));
            SET_VECTOR_ELT(result, 15, s);
            UNPROTECT(1);
        }
        
        // 16: significant_vertices (logical)
        {
            const R_xlen_t n = static_cast<R_xlen_t>(res.permutation_tests->significant_vertices.size());
            SEXP s = PROTECT(Rf_allocVector(LGLSXP, n));
            int* p = LOGICAL(s);
            for (R_xlen_t i = 0; i < n; ++i) {
                p[i] = res.permutation_tests->significant_vertices[static_cast<size_t>(i)] ? TRUE : FALSE;
            }
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
        for (int i = 14; i <= 17; ++i) {
            SET_VECTOR_ELT(result, i, R_NilValue);
        }
    }

    // Set names
    const char* names[] = {
        "graph_diameter",
        "packing_radius",
        "packing_vertices",
        "grid_opt_bw",
        "predictions",
        "errors",
        "scale",
        "grid_predictions",
        "bb_predictions",
        "cri_lower",
        "cri_upper",
        "null_predictions",
        "null_predictions_cri_lower",
        "null_predictions_cri_upper",
        "p_values",
        "effect_sizes",
        "significant_vertices",
        "null_distances"
    };

    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    for (int i = 0; i < n_elements; ++i) {
        SET_STRING_ELT(result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    UNPROTECT(2); // result, result_names
    return result;
}