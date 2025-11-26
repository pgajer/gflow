#include "set_wgraph.hpp"
#include "lcor.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#include <R.h>
#include <Rinternals.h>

// Forward declare conversion utilities available from SEXP_cpp_conversion_utils.cpp
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

/**
 * @brief SEXP interface for computing local correlation coefficients with flexible edge difference types
 *
 * @details
 * This function provides an R interface to the local correlation (lcor) coefficient computation,
 * which generalizes the co-monotonicity framework to support flexible edge difference types.
 * It measures the extent to which two functions y and z vary together across the edges of a
 * weighted graph, with appropriate handling for both continuous and compositional data.
 *
 * The key innovation is the ability to specify how edge differences are computed:
 * - DIFFERENCE: Standard differences Δ_e f = f(u) - f(v) for continuous data
 * - LOGRATIO: Log-ratios Δ_e f = log((f(u) + ε) / (f(v) + ε)) for compositional data
 *
 * This flexibility enables proper geometric treatment of mixed data types, such as:
 * - Continuous response vs. compositional features (e.g., disease severity vs. bacterial abundance)
 * - Compositional co-occurrence analysis (e.g., two taxa in microbiome data)
 * - Any combination of continuous and compositional variables
 *
 * The local correlation coefficient at each vertex quantifies the alignment of directional
 * changes within the vertex's neighborhood:
 *
 *   lcor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * where the sum is over edges incident to vertex v. This correlation-style normalization
 * makes the coefficient scale-invariant and interpretable as the cosine of the angle between
 * gradient vectors (or their log-ratio analogs for compositional data).
 *
 * WEIGHTING SCHEMES:
 * Three weighting schemes are supported via the type parameter:
 *
 * - "unit": Uses w_e = 1 for all edges, treating all edges equally regardless of length.
 *   Appropriate for combinatorial graphs or when edge lengths are comparable.
 *
 * - "derivative": Uses w_e = 1/(ℓ_e)² where ℓ_e is the edge length. This normalizes by
 *   edge length squared, making the measure analogous to comparing directional derivatives.
 *   Natural for functions on geometric graphs where positions are meaningful.
 *
 * - "sign": Uses only the sign of the product, computing proportion of edges with directional
 *   agreement. Robust but loses magnitude information. (Primarily for compatibility; typically
 *   not used with correlation-style normalization.)
 *
 * EDGE DIFFERENCE TYPES:
 * The y.diff.type and z.diff.type parameters control how edge differences are computed:
 *
 * - "difference": Standard differences Δ_e f = f(u) - f(v)
 *   Use for continuous data in Euclidean space (temperatures, expression levels, clinical scores)
 *
 * - "logratio": Log-ratios Δ_e f = log((f(u) + ε) / (f(v) + ε))
 *   Use for compositional data (relative abundances, proportions)
 *   The log transformation maps multiplicative changes to additive scale, corresponding to
 *   the Aitchison distance on the simplex
 *
 * PSEUDOCOUNT (EPSILON):
 * For log-ratio transformations, a pseudocount ε is added to prevent log(0):
 * - If epsilon = 0 (default): Computed adaptively as 1e-6 × min(non-zero values) for each function
 * - If epsilon > 0: Uses the specified value explicitly
 * Only applied when the corresponding diff.type is "logratio"
 *
 * WINSORIZATION:
 * For robustness against outliers, extreme edge differences can be clipped:
 * - If winsorize.quantile = 0 (default): No winsorization, uses efficient one-pass algorithm
 * - If winsorize.quantile > 0: Clips differences to [q, 1-q] percentiles (e.g., 0.05 for 5th/95th)
 *   Requires two-pass algorithm with extra memory overhead
 *
 * RETURN VALUE:
 * The function returns detailed results including vertex-level coefficients, edge-level data
 * for diagnostics, and winsorization bounds (when applicable). This enables both spatial
 * analysis and validation/tuning of the method.
 *
 * COMPUTATIONAL COMPLEXITY:
 * - One-pass (no winsorization): O(|E|) time, O(n) space
 * - Two-pass (with winsorization): O(|E| + |E|log|E|) time, O(|E|) space
 * For sparse graphs with bounded degree, both are effectively O(n) in the number of vertices.
 *
 * @param s_adj_list R list of integer vectors representing the graph's adjacency structure.
 *   Element i contains 0-based indices of vertices adjacent to vertex i. The graph is
 *   undirected, so edge [i,j] should appear in both adjacency lists. Must have length
 *   equal to the number of vertices.
 *
 * @param s_weight_list R list of numeric vectors containing edge weights (lengths).
 *   Element i contains weights corresponding to the edges in s_adj_list[[i]]. Must have
 *   the same structure as s_adj_list. Weights should be positive.
 *
 * @param s_y R numeric vector of response function values at each vertex. Length must
 *   equal the number of vertices (i.e., length(s_adj_list)).
 *
 * @param s_z R numeric vector of feature function values at each vertex. Length must
 *   equal the number of vertices.
 *
 * @param s_type R character scalar specifying the weighting scheme. Must be one of:
 *   "unit", "derivative", or "sign". Case-sensitive.
 *
 * @param s_y_diff_type R character scalar specifying how to compute edge differences for y.
 *   Must be one of: "difference" or "logratio". Case-sensitive.
 *
 * @param s_z_diff_type R character scalar specifying how to compute edge differences for z.
 *   Must be one of: "difference" or "logratio". Case-sensitive.
 *
 * @param s_epsilon R numeric scalar specifying the pseudocount for log-ratio transformations.
 *   If 0 (default), computed adaptively as 1e-6 × min(non-zero values). Only used when
 *   the corresponding diff.type is "logratio".
 *
 * @param s_winsorize_quantile R numeric scalar specifying the winsorization quantile.
 *   If 0 (default), no winsorization. If > 0 (e.g., 0.05), clips edge differences to
 *   [q, 1-q] percentiles for robustness against outliers.
 *
 * @return R list with components:
 *   - vertex.coefficients: Numeric vector of local correlation coefficients at each vertex
 *   - vertex.delta.y: R list of numeric vectors, element i contains edge differences for y
 *     at vertex i (one value per neighbor)
 *   - vertex.delta.z: R list of numeric vectors, element i contains edge differences for z
 *     at vertex i (one value per neighbor)
 *   - vertex.weights: R list of numeric vectors, element i contains edge weights at vertex i
 *     (one value per neighbor)
 *   - all.delta.y: Numeric vector containing all edge differences for y across the entire graph
 *     (before winsorization if applicable). Useful for diagnostics and distribution analysis.
 *   - all.delta.z: Numeric vector containing all edge differences for z across the entire graph
 *     (before winsorization if applicable). Useful for diagnostics and distribution analysis.
 *   - y.lower: Scalar giving the lower winsorization bound for y differences (or -Inf if no winsorization)
 *   - y.upper: Scalar giving the upper winsorization bound for y differences (or +Inf if no winsorization)
 *   - z.lower: Scalar giving the lower winsorization bound for z differences (or -Inf if no winsorization)
 *   - z.upper: Scalar giving the upper winsorization bound for z differences (or +Inf if no winsorization)
 *
 * @note All vertex indices in input are 0-based (C++ convention). The R wrapper function
 *   handles conversion to/from R's 1-based indexing.
 *
 * @note Memory is managed using R's PROTECT/UNPROTECT mechanism. Each allocated SEXP
 *   must be protected before use and unprotected when no longer needed.
 *
 * @note The vertex.delta.y, vertex.delta.z, and vertex.weights lists maintain vertex-neighborhood
 *   structure. For a vertex with k neighbors, these lists contain k values corresponding to the
 *   k edges incident to that vertex. After winsorization, vertex.delta.y/z contain clipped values.
 *
 * @note The all.delta.y and all.delta.z vectors contain edge differences BEFORE winsorization,
 *   enabling analysis of the raw distribution. For an undirected graph with |E| edges, these
 *   vectors have length 2|E| because each edge is traversed from both endpoints.
 *
 * @throws Rf_error if s_adj_list is not a list
 * @throws Rf_error if s_weight_list is not a list
 * @throws Rf_error if s_y is not a numeric vector
 * @throws Rf_error if s_z is not a numeric vector
 * @throws Rf_error if s_type is not a single character string
 * @throws Rf_error if s_y_diff_type is not a single character string
 * @throws Rf_error if s_z_diff_type is not a single character string
 * @throws Rf_error if s_epsilon is not a numeric scalar
 * @throws Rf_error if s_winsorize_quantile is not a numeric scalar
 * @throws Rf_error if s_type is not "unit", "derivative", or "sign"
 * @throws Rf_error if s_y_diff_type is not "difference" or "logratio"
 * @throws Rf_error if s_z_diff_type is not "difference" or "logratio"
 * @throws Rf_error if length(s_y) does not match number of vertices
 * @throws Rf_error if length(s_z) does not match number of vertices
 *
 * @see lcor_type_t for detailed description of weighting schemes
 * @see edge_diff_type_t for detailed description of edge difference types
 * @see set_wgraph_t::lcor for the underlying C++ implementation
 * @see lcor_result_t for the C++ result structure
 *
 * @example
 * // In R (with 0-based indexing for adjacency):
 * // adj.list <- list(c(1,2), c(0,2), c(0,1))  # Triangle graph
 * // weight.list <- list(c(1.0,1.0), c(1.0,1.0), c(1.0,1.0))
 * // 
 * // # Continuous response vs compositional feature
 * // severity <- c(0.1, 0.5, 0.9)  # Disease severity
 * // abundance <- c(0.01, 0.10, 0.50)  # Bacterial abundance
 * // 
 * // result <- .Call(S_lcor, adj.list, weight.list, 
 * //                 severity, abundance,
 * //                 "derivative", "difference", "logratio",
 * //                 0.0, 0.0)
 * //
 * // # Access results
 * // result$vertex.coefficients  # Correlation at each vertex
 * // result$all.delta.y         # Distribution of y differences
 * // result$all.delta.z         # Distribution of log-ratios for z
 */
extern "C" SEXP S_lcor_instrumented(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_winsorize_quantile
) {
    // ---- Input validation and conversion ----
    
    // Convert adjacency and weight lists from R to C++
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    
    // Convert y vector
    if (!Rf_isReal(s_y)) {
        Rf_error("s_y must be a numeric vector");
    }
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);
    
    // Convert z vector
    if (!Rf_isReal(s_z)) {
        Rf_error("s_z must be a numeric vector");
    }
    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);
    
    // Validate dimensions
    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices) {
        Rf_error("Length of y (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_y), static_cast<int>(n_vertices));
    }
    if (n_z != n_vertices) {
        Rf_error("Length of z (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_z), static_cast<int>(n_vertices));
    }
    
    // Convert type string (weighting scheme)
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);
    
    lcor_type_t lcor_type;
    if (type_str == "unit") {
        lcor_type = lcor_type_t::UNIT;
    } else if (type_str == "derivative") {
        lcor_type = lcor_type_t::DERIVATIVE;
    } else if (type_str == "sign") {
        lcor_type = lcor_type_t::SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'unit', 'derivative', or 'sign'",
                 type_cstr);
    }
    
    // Convert y_diff_type string
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    const char* y_diff_type_cstr = CHAR(STRING_ELT(s_y_diff_type, 0));
    std::string y_diff_type_str(y_diff_type_cstr);
    
    edge_diff_type_t y_diff_type;
    if (y_diff_type_str == "difference") {
        y_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (y_diff_type_str == "logratio") {
        y_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid y_diff_type '%s'. Must be 'difference' or 'logratio'",
                 y_diff_type_cstr);
    }
    
    // Convert z_diff_type string
    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    const char* z_diff_type_cstr = CHAR(STRING_ELT(s_z_diff_type, 0));
    std::string z_diff_type_str(z_diff_type_cstr);
    
    edge_diff_type_t z_diff_type;
    if (z_diff_type_str == "difference") {
        z_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (z_diff_type_str == "logratio") {
        z_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid z_diff_type '%s'. Must be 'difference' or 'logratio'",
                 z_diff_type_cstr);
    }
    
    // Convert epsilon (scalar)
    if (!Rf_isReal(s_epsilon) || LENGTH(s_epsilon) != 1) {
        Rf_error("s_epsilon must be a single numeric value");
    }
    double epsilon = REAL(s_epsilon)[0];
    
    // Convert winsorize_quantile (scalar)
    if (!Rf_isReal(s_winsorize_quantile) || LENGTH(s_winsorize_quantile) != 1) {
        Rf_error("s_winsorize_quantile must be a single numeric value");
    }
    double winsorize_quantile = REAL(s_winsorize_quantile)[0];
    
    // ---- Build graph ----
    set_wgraph_t graph(adj_list, weight_list);
    
    // ---- Compute local correlation ----
    lcor_result_t result = graph.lcor_instrumented(
        y, z,
        lcor_type,
        y_diff_type,
        z_diff_type,
        epsilon,
        winsorize_quantile
    );
    
    // ---- Build R return object ----
    
    // Create list with 10 components
    const int n_components = 10;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_components));
    
    int idx = 0;
    
    // Component 1: vertex.coefficients (numeric vector)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertex.coefficients"));
    SEXP r_vertex_coeffs = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* vertex_coeffs_ptr = REAL(r_vertex_coeffs);
    for (size_t i = 0; i < n_vertices; ++i) {
        vertex_coeffs_ptr[i] = result.vertex_coefficients[i];
    }
    SET_VECTOR_ELT(r_result, idx++, r_vertex_coeffs);
    UNPROTECT(1); // r_vertex_coeffs
    
    // Component 2: vertex.delta.y (list of numeric vectors)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertex.delta.y"));
    SEXP r_vertex_delta_y = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    for (size_t v = 0; v < n_vertices; ++v) {
        size_t n_neighbors = result.vertex_delta_y[v].size();
        SEXP r_delta_y_v = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        double* delta_y_ptr = REAL(r_delta_y_v);
        for (size_t i = 0; i < n_neighbors; ++i) {
            delta_y_ptr[i] = result.vertex_delta_y[v][i];
        }
        SET_VECTOR_ELT(r_vertex_delta_y, v, r_delta_y_v);
        UNPROTECT(1); // r_delta_y_v
    }
    SET_VECTOR_ELT(r_result, idx++, r_vertex_delta_y);
    UNPROTECT(1); // r_vertex_delta_y
    
    // Component 3: vertex.delta.z (list of numeric vectors)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertex.delta.z"));
    SEXP r_vertex_delta_z = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    for (size_t v = 0; v < n_vertices; ++v) {
        size_t n_neighbors = result.vertex_delta_z[v].size();
        SEXP r_delta_z_v = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        double* delta_z_ptr = REAL(r_delta_z_v);
        for (size_t i = 0; i < n_neighbors; ++i) {
            delta_z_ptr[i] = result.vertex_delta_z[v][i];
        }
        SET_VECTOR_ELT(r_vertex_delta_z, v, r_delta_z_v);
        UNPROTECT(1); // r_delta_z_v
    }
    SET_VECTOR_ELT(r_result, idx++, r_vertex_delta_z);
    UNPROTECT(1); // r_vertex_delta_z
    
    // Component 4: vertex.weights (list of numeric vectors)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("vertex.weights"));
    SEXP r_vertex_weights = PROTECT(Rf_allocVector(VECSXP, n_vertices));
    for (size_t v = 0; v < n_vertices; ++v) {
        size_t n_neighbors = result.vertex_weights[v].size();
        SEXP r_weights_v = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        double* weights_ptr = REAL(r_weights_v);
        for (size_t i = 0; i < n_neighbors; ++i) {
            weights_ptr[i] = result.vertex_weights[v][i];
        }
        SET_VECTOR_ELT(r_vertex_weights, v, r_weights_v);
        UNPROTECT(1); // r_weights_v
    }
    SET_VECTOR_ELT(r_result, idx++, r_vertex_weights);
    UNPROTECT(1); // r_vertex_weights
    
    // Component 5: all.delta.y (numeric vector - all edge differences for y)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("all.delta.y"));
    size_t n_all_delta_y = result.all_delta_y.size();
    SEXP r_all_delta_y = PROTECT(Rf_allocVector(REALSXP, n_all_delta_y));
    double* all_delta_y_ptr = REAL(r_all_delta_y);
    for (size_t i = 0; i < n_all_delta_y; ++i) {
        all_delta_y_ptr[i] = result.all_delta_y[i];
    }
    SET_VECTOR_ELT(r_result, idx++, r_all_delta_y);
    UNPROTECT(1); // r_all_delta_y
    
    // Component 6: all.delta.z (numeric vector - all edge differences for z)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("all.delta.z"));
    size_t n_all_delta_z = result.all_delta_z.size();
    SEXP r_all_delta_z = PROTECT(Rf_allocVector(REALSXP, n_all_delta_z));
    double* all_delta_z_ptr = REAL(r_all_delta_z);
    for (size_t i = 0; i < n_all_delta_z; ++i) {
        all_delta_z_ptr[i] = result.all_delta_z[i];
    }
    SET_VECTOR_ELT(r_result, idx++, r_all_delta_z);
    UNPROTECT(1); // r_all_delta_z
    
    // Component 7: y.lower (scalar - winsorization lower bound for y)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("y.lower"));
    SEXP r_y_lower = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_y_lower)[0] = result.y_lower;
    SET_VECTOR_ELT(r_result, idx++, r_y_lower);
    UNPROTECT(1); // r_y_lower
    
    // Component 8: y.upper (scalar - winsorization upper bound for y)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("y.upper"));
    SEXP r_y_upper = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_y_upper)[0] = result.y_upper;
    SET_VECTOR_ELT(r_result, idx++, r_y_upper);
    UNPROTECT(1); // r_y_upper
    
    // Component 9: z.lower (scalar - winsorization lower bound for z)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("z.lower"));
    SEXP r_z_lower = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_z_lower)[0] = result.z_lower;
    SET_VECTOR_ELT(r_result, idx++, r_z_lower);
    UNPROTECT(1); // r_z_lower
    
    // Component 10: z.upper (scalar - winsorization upper bound for z)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("z.upper"));
    SEXP r_z_upper = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(r_z_upper)[0] = result.z_upper;
    SET_VECTOR_ELT(r_result, idx++, r_z_upper);
    UNPROTECT(1); // r_z_upper
    
    // Set names attribute on result list
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    
    // Unprotect the result list and names vector
    UNPROTECT(2); // r_result, r_names
    
    return r_result;
}

/**
 * @brief SEXP interface for computing local correlation coefficients (production version)
 *
 * @details
 * Streamlined R interface that returns only the vertex-level local correlation
 * coefficients without diagnostic data. This is the preferred interface for
 * production use where edge-level details are not needed.
 *
 * The local correlation coefficient at each vertex quantifies the alignment of
 * directional changes between functions y and z within the vertex's neighborhood:
 *
 *   lcor(y,z)(v) = Σ w_e Δ_e y · Δ_e z / √(Σ w_e (Δ_e y)²) √(Σ w_e (Δ_e z)²)
 *
 * PERFORMANCE:
 * This version avoids the O(|E|) memory overhead of storing edge differences,
 * making it more efficient for large graphs when diagnostic data is not required.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency structure)
 * @param s_weight_list R list of numeric vectors (edge weights/lengths)
 * @param s_y R numeric vector of response function values
 * @param s_z R numeric vector of feature function values
 * @param s_type R character: "unit", "derivative", or "sign"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric scalar: pseudocount (0 = adaptive)
 * @param s_winsorize_quantile R numeric scalar: winsorization quantile (0 = none)
 *
 * @return R numeric vector of local correlation coefficients (one per vertex)
 *
 * @throws Rf_error for invalid input types or dimension mismatches
 *
 * @note For diagnostic output including edge differences and winsorization bounds,
 *       use S_lcor_instrumented() instead.
 *
 * @see S_lcor_instrumented() for the version with full diagnostic output
 * @see set_wgraph_t::lcor() for the underlying C++ implementation
 */
extern "C" SEXP S_lcor(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_winsorize_quantile
) {
    // ---- Input validation and conversion ----

    // Convert adjacency and weight lists from R to C++
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y vector
    if (!Rf_isReal(s_y)) {
        Rf_error("s_y must be a numeric vector");
    }
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    // Convert z vector
    if (!Rf_isReal(s_z)) {
        Rf_error("s_z must be a numeric vector");
    }
    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    // Validate dimensions
    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices) {
        Rf_error("Length of y (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_y), static_cast<int>(n_vertices));
    }
    if (n_z != n_vertices) {
        Rf_error("Length of z (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_z), static_cast<int>(n_vertices));
    }

    // Convert type string (weighting scheme)
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);

    lcor_type_t lcor_type;
    if (type_str == "unit") {
        lcor_type = lcor_type_t::UNIT;
    } else if (type_str == "derivative") {
        lcor_type = lcor_type_t::DERIVATIVE;
    } else if (type_str == "sign") {
        lcor_type = lcor_type_t::SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'unit', 'derivative', or 'sign'",
                 type_cstr);
    }

    // Convert y_diff_type string
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    const char* y_diff_type_cstr = CHAR(STRING_ELT(s_y_diff_type, 0));
    std::string y_diff_type_str(y_diff_type_cstr);

    edge_diff_type_t y_diff_type;
    if (y_diff_type_str == "difference") {
        y_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (y_diff_type_str == "logratio") {
        y_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid y_diff_type '%s'. Must be 'difference' or 'logratio'",
                 y_diff_type_cstr);
    }

    // Convert z_diff_type string
    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    const char* z_diff_type_cstr = CHAR(STRING_ELT(s_z_diff_type, 0));
    std::string z_diff_type_str(z_diff_type_cstr);

    edge_diff_type_t z_diff_type;
    if (z_diff_type_str == "difference") {
        z_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (z_diff_type_str == "logratio") {
        z_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid z_diff_type '%s'. Must be 'difference' or 'logratio'",
                 z_diff_type_cstr);
    }

    // Convert epsilon (scalar)
    if (!Rf_isReal(s_epsilon) || LENGTH(s_epsilon) != 1) {
        Rf_error("s_epsilon must be a single numeric value");
    }
    double epsilon = REAL(s_epsilon)[0];

    // Convert winsorize_quantile (scalar)
    if (!Rf_isReal(s_winsorize_quantile) || LENGTH(s_winsorize_quantile) != 1) {
        Rf_error("s_winsorize_quantile must be a single numeric value");
    }
    double winsorize_quantile = REAL(s_winsorize_quantile)[0];

    // ---- Build graph ----
    set_wgraph_t graph(adj_list, weight_list);

    // ---- Compute local correlation ----
    auto vertex_coefficients = graph.lcor(
        y, z,
        lcor_type,
        y_diff_type,
        z_diff_type,
        epsilon,
        winsorize_quantile
    );

    // ---- Build R return object: vertex.coefficients (numeric vector) ----

    SEXP r_vertex_coeffs = PROTECT(Rf_allocVector(REALSXP, n_vertices));
    double* vertex_coeffs_ptr = REAL(r_vertex_coeffs);
    for (size_t i = 0; i < n_vertices; ++i) {
        vertex_coeffs_ptr[i] = vertex_coefficients[i];
    }
    UNPROTECT(1); // r_vertex_coeffs

    return r_vertex_coeffs;
}

/**
 * @brief SEXP interface for computing local correlation between a vector and matrix columns
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency structure)
 * @param s_weight_list R list of numeric vectors (edge weights/lengths)
 * @param s_y R numeric vector of response function values
 * @param s_Z R numeric matrix of feature function values
 * @param s_type R character: "unit", "derivative", or "sign"
 * @param s_y_diff_type R character: "difference" or "logratio"
 * @param s_z_diff_type R character: "difference" or "logratio"
 * @param s_epsilon R numeric scalar: pseudocount (0 = adaptive)
 * @param s_winsorize_quantile R numeric scalar: winsorization quantile (0 = none)
 * @return R list with column.coefficients, mean.coefficients, etc.
 */
extern "C" SEXP S_lcor_vector_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_Z,
    SEXP s_type,
    SEXP s_y_diff_type,
    SEXP s_z_diff_type,
    SEXP s_epsilon,
    SEXP s_winsorize_quantile
) {
    // ---- Input validation and conversion ----

    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y vector
    if (!Rf_isReal(s_y)) {
        Rf_error("s_y must be a numeric vector");
    }
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    // Convert Z matrix
    if (!Rf_isMatrix(s_Z) || !Rf_isReal(s_Z)) {
        Rf_error("s_Z must be a numeric matrix");
    }

    int* dims = INTEGER(Rf_getAttrib(s_Z, R_DimSymbol));
    int n_rows = dims[0];
    int n_cols = dims[1];

    double* Z_ptr = REAL(s_Z);
    Eigen::Map<Eigen::MatrixXd> Z(Z_ptr, n_rows, n_cols);

    // Validate dimensions
    size_t n_vertices = adj_list.size();
    if (n_y != n_vertices) {
        Rf_error("Length of y (%d) does not match number of vertices (%d)",
                 static_cast<int>(n_y), static_cast<int>(n_vertices));
    }
    if (static_cast<size_t>(n_rows) != n_vertices) {
        Rf_error("Number of rows in Z (%d) does not match number of vertices (%d)",
                 n_rows, static_cast<int>(n_vertices));
    }

    // Convert type string
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);

    lcor_type_t lcor_type;
    if (type_str == "unit") {
        lcor_type = lcor_type_t::UNIT;
    } else if (type_str == "derivative") {
        lcor_type = lcor_type_t::DERIVATIVE;
    } else if (type_str == "sign") {
        lcor_type = lcor_type_t::SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'unit', 'derivative', or 'sign'",
                 type_cstr);
    }

    // Convert y_diff_type string
    if (!Rf_isString(s_y_diff_type) || LENGTH(s_y_diff_type) != 1) {
        Rf_error("s_y_diff_type must be a single character string");
    }
    const char* y_diff_type_cstr = CHAR(STRING_ELT(s_y_diff_type, 0));
    std::string y_diff_type_str(y_diff_type_cstr);

    edge_diff_type_t y_diff_type;
    if (y_diff_type_str == "difference") {
        y_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (y_diff_type_str == "logratio") {
        y_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid y_diff_type '%s'. Must be 'difference' or 'logratio'",
                 y_diff_type_cstr);
    }

    // Convert z_diff_type string
    if (!Rf_isString(s_z_diff_type) || LENGTH(s_z_diff_type) != 1) {
        Rf_error("s_z_diff_type must be a single character string");
    }
    const char* z_diff_type_cstr = CHAR(STRING_ELT(s_z_diff_type, 0));
    std::string z_diff_type_str(z_diff_type_cstr);

    edge_diff_type_t z_diff_type;
    if (z_diff_type_str == "difference") {
        z_diff_type = edge_diff_type_t::DIFFERENCE;
    } else if (z_diff_type_str == "logratio") {
        z_diff_type = edge_diff_type_t::LOGRATIO;
    } else {
        Rf_error("Invalid z_diff_type '%s'. Must be 'difference' or 'logratio'",
                 z_diff_type_cstr);
    }

    // Convert epsilon
    if (!Rf_isReal(s_epsilon) || LENGTH(s_epsilon) != 1) {
        Rf_error("s_epsilon must be a single numeric value");
    }
    double epsilon = REAL(s_epsilon)[0];

    // Convert winsorize_quantile
    if (!Rf_isReal(s_winsorize_quantile) || LENGTH(s_winsorize_quantile) != 1) {
        Rf_error("s_winsorize_quantile must be a single numeric value");
    }
    double winsorize_quantile = REAL(s_winsorize_quantile)[0];

    // ---- Build graph and compute ----
    set_wgraph_t graph(adj_list, weight_list);

    lcor_vector_matrix_result_t result = graph.lcor_vector_matrix(
        y, Z, lcor_type, y_diff_type, z_diff_type, epsilon, winsorize_quantile
    );

    // ---- Build R return object ----

    const int n_components = 8;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // column.coefficients (list of vectors)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("column.coefficients"));
    SEXP r_col_coeffs = PROTECT(Rf_allocVector(VECSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        SEXP r_col = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        double* col_ptr = REAL(r_col);
        for (size_t i = 0; i < n_vertices; ++i) {
            col_ptr[i] = result.column_coefficients[col][i];
        }
        SET_VECTOR_ELT(r_col_coeffs, col, r_col);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(r_result, idx++, r_col_coeffs);
    UNPROTECT(1);

    // mean.coefficients
    SET_STRING_ELT(r_names, idx, Rf_mkChar("mean.coefficients"));
    SEXP r_means = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_means)[col] = result.mean_coefficients[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_means);
    UNPROTECT(1);

    // median.coefficients
    SET_STRING_ELT(r_names, idx, Rf_mkChar("median.coefficients"));
    SEXP r_medians = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_medians)[col] = result.median_coefficients[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_medians);
    UNPROTECT(1);

    // n.positive
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.positive"));
    SEXP r_n_positive = PROTECT(Rf_allocVector(INTSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        INTEGER(r_n_positive)[col] = static_cast<int>(result.n_positive[col]);
    }
    SET_VECTOR_ELT(r_result, idx++, r_n_positive);
    UNPROTECT(1);

    // n.negative
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.negative"));
    SEXP r_n_negative = PROTECT(Rf_allocVector(INTSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        INTEGER(r_n_negative)[col] = static_cast<int>(result.n_negative[col]);
    }
    SET_VECTOR_ELT(r_result, idx++, r_n_negative);
    UNPROTECT(1);

    // n.zero
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.zero"));
    SEXP r_n_zero = PROTECT(Rf_allocVector(INTSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        INTEGER(r_n_zero)[col] = static_cast<int>(result.n_zero[col]);
    }
    SET_VECTOR_ELT(r_result, idx++, r_n_zero);
    UNPROTECT(1);

    // y.lower
    SET_STRING_ELT(r_names, idx, Rf_mkChar("y.lower"));
    SEXP r_y_lower = PROTECT(Rf_ScalarReal(result.y_lower));
    SET_VECTOR_ELT(r_result, idx++, r_y_lower);
    UNPROTECT(1);

    // y.upper
    SET_STRING_ELT(r_names, idx, Rf_mkChar("y.upper"));
    SEXP r_y_upper = PROTECT(Rf_ScalarReal(result.y_upper));
    SET_VECTOR_ELT(r_result, idx++, r_y_upper);
    UNPROTECT(1);

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    UNPROTECT(2);  // r_result, r_names

    return r_result;
}
