/**
 * @brief Implementation of Adaptive GEMALO algorithm
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

    // Convert input parameters using R's C API
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Extract scalar parameters using R's C API

    // geodesic parameter
    size_t min_path_size = (size_t)INTEGER(s_min_path_size)[0];
    // packing parameters
    size_t n_packing_vertices = (size_t)INTEGER(s_n_packing_vertices)[0];
    size_t max_packing_iterations = (size_t)INTEGER(s_max_packing_iterations)[0];
    double packing_precision = REAL(s_packing_precision)[0];

    // bw parameters
    size_t n_bws = (size_t)INTEGER(s_n_bws)[0];
    bool log_grid = (LOGICAL(s_log_grid)[0] == 1);
    double min_bw_factor = REAL(s_min_bw_factor)[0];
    double max_bw_factor = REAL(s_max_bw_factor)[0];

    // kernel parameters
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    size_t kernel_type = (size_t)INTEGER(s_kernel_type)[0];
    // model parameters
    double model_tolerance = REAL(s_model_tolerance)[0];
    double model_blending_coef = REAL(s_model_blending_coef)[0];
    // Bayesian bootstrap parameters
    size_t n_bb = (size_t)INTEGER(s_n_bb)[0];
    double cri_probability = REAL(s_cri_probability)[0];
    // permutation parameters
    size_t n_perms = (size_t)INTEGER(s_n_perms)[0];
    // verbose
    bool verbose = (LOGICAL(s_verbose)[0] == 1);

    uniform_grid_graph_t grid_graph;
    size_t n_vertices = adj_list.size();
    if (n_packing_vertices < n_vertices) {
        // The returned uniform_grid_graph_t object inherits both the graph_diameter
        // and max_packing_radius values from the intermediate set_wgraph_t object,
        // making these calculated values available for further analysis.
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
        auto [end1, diam] = grid_graph.get_vertex_eccentricity(0);  // Start from vertex 0
        auto [end2, diameter] = grid_graph.get_vertex_eccentricity(end1);
        grid_graph.graph_diameter = diameter;
    }

    // Call the C++ function
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

    // Create the return list using R's C API
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
        "null_distances",
        NULL
    };

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    // Create list and protect it
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_elements));
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));

    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, result_names);

    // Helper function to convert vector to SEXP
    auto create_numeric_vector = [](const std::vector<double>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size()));
        double* ptr = REAL(r_vec);
        std::copy(vec.begin(), vec.end(), ptr);
        return r_vec;
    };

    // Helper function to convert matrix to SEXP
    auto create_numeric_matrix = [](const std::vector<std::vector<double>>& mat) -> SEXP {
        if (mat.empty()) return R_NilValue;
        size_t nrow = mat.size();
        size_t ncol = mat[0].size();

        SEXP r_mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
        double* ptr = REAL(r_mat);

        for (size_t i = 0; i < nrow; ++i) {
            for (size_t j = 0; j < ncol; ++j) {
                ptr[i + j * nrow] = mat[i][j];  // Column-major order for R
            }
        }
        return r_mat;
    };

    // Helper function to convert boolean vector to SEXP
    auto create_logical_vector = [](const std::vector<bool>& vec) -> SEXP {
        SEXP r_vec = PROTECT(Rf_allocVector(LGLSXP, vec.size()));
        int* ptr = LOGICAL(r_vec);
        for (size_t i = 0; i < vec.size(); ++i) {
            ptr[i] = vec[i];
        }
        return r_vec;
    };

    // Set graph diameter
    SEXP diam = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(diam)[0] = grid_graph.graph_diameter;
    SET_VECTOR_ELT(result, 0, diam);

    // Set packing radius
    SEXP packing_radius = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(packing_radius)[0] = grid_graph.max_packing_radius;
    SET_VECTOR_ELT(result, 1, packing_radius);

    // Create packing vertices vector (1-based indices)
    int actual_n_packing_vertices = grid_graph.grid_vertices.size();
    SEXP packing_vertices = PROTECT(Rf_allocVector(INTSXP, actual_n_packing_vertices));

    int counter = 0;
    for (const auto& i : grid_graph.grid_vertices) {
        // Convert to 1-based indices for R
        INTEGER(packing_vertices)[counter++] = i + 1;
    }
    SET_VECTOR_ELT(result, 2, packing_vertices);

    // Set other numeric vectors
    SET_VECTOR_ELT(result, 3, create_numeric_vector(res.grid_opt_bw));
    SET_VECTOR_ELT(result, 4, create_numeric_vector(res.predictions));
    SET_VECTOR_ELT(result, 5, create_numeric_vector(res.errors));
    SET_VECTOR_ELT(result, 6, create_numeric_vector(res.scale));
    SET_VECTOR_ELT(result, 7, create_numeric_vector(map_to_vector(res.grid_predictions_map)));

    // Set bootstrap fields if available
    if (!res.bb_predictions.empty()) {
        SET_VECTOR_ELT(result, 8, create_numeric_matrix(res.bb_predictions));
        SET_VECTOR_ELT(result, 9, create_numeric_vector(res.cri_lower));
        SET_VECTOR_ELT(result, 10, create_numeric_vector(res.cri_upper));
    } else {
        // Set to NULL if not available
        SET_VECTOR_ELT(result, 8, R_NilValue);
        SET_VECTOR_ELT(result, 9, R_NilValue);
        SET_VECTOR_ELT(result, 10, R_NilValue);
    }

    // Set permutation fields if available
    if (!res.null_predictions.empty()) {
        SET_VECTOR_ELT(result, 11, create_numeric_matrix(res.null_predictions));
        SET_VECTOR_ELT(result, 12, create_numeric_vector(res.null_predictions_cri_lower));
        SET_VECTOR_ELT(result, 13, create_numeric_vector(res.null_predictions_cri_upper));
    } else {
        SET_VECTOR_ELT(result, 11, R_NilValue);
        SET_VECTOR_ELT(result, 12, R_NilValue);
        SET_VECTOR_ELT(result, 13, R_NilValue);
    }

    // Set permutation test results if available
    if (res.permutation_tests) {
        SET_VECTOR_ELT(result, 14, create_numeric_vector(res.permutation_tests->p_values));
        SET_VECTOR_ELT(result, 15, create_numeric_vector(res.permutation_tests->effect_sizes));
        SET_VECTOR_ELT(result, 16, create_logical_vector(res.permutation_tests->significant_vertices));
        SET_VECTOR_ELT(result, 17, create_numeric_vector(res.permutation_tests->null_distances));
    } else {
        for (int i = 14; i <= 17; i++) {
            SET_VECTOR_ELT(result, i, R_NilValue);
        }
    }

    // Count protected objects: result, result_names, diam, and all non-NULL vectors/matrices
    int total_protected = 2;  // Start with result, result_names
    for (int i = 0; i < n_elements; i++) {
        if (VECTOR_ELT(result, i) != R_NilValue) total_protected++;
    }

    UNPROTECT(total_protected);
    return result;
}
