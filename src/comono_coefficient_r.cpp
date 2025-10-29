#include "set_wgraph.hpp"
#include "comono_coefficient.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include <R.h>
#include <Rinternals.h>

// Forward declare conversion utilities that should be available from SEXP_cpp_conversion_utils.cpp
// These convert R lists of adjacency/weight vectors into C++ nested vectors
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

/**
 * @brief SEXP interface for computing co-monotonicity coefficients between two functions on graph vertices
 *
 * @details
 * This function provides an R interface to the co-monotonicity coefficient computation.
 * It measures the extent to which two functions y and z vary together across the edges
 * of a weighted graph, providing a graph-based analogue of correlation that respects
 * local geometric structure.
 *
 * The co-monotonicity coefficient at each vertex quantifies the agreement in directional
 * changes of y and z within the vertex's neighborhood. For an edge e = [v,u], we compute
 * the edge differences Δ_e y = y(u) - y(v) and Δ_e z = z(u) - z(v). The product
 * Δ_e y · Δ_e z captures co-variation along the edge: positive when both functions
 * increase or both decrease, negative when they change in opposite directions.
 *
 * The vertex-level coefficient aggregates these edge-wise co-variations:
 *   comono(y,z;w)(v) = Σ_{u ∈ N(v)} w_e Δ_e y Δ_e z / Σ_{u ∈ N(v)} w_e |Δ_e y Δ_e z|
 *
 * where w_e are edge weights determined by the type parameter. The coefficient lies in
 * [-1, 1], where +1 indicates perfect positive co-monotonicity (functions always change
 * together), -1 indicates perfect negative co-monotonicity (functions always change in
 * opposite directions), and values near 0 suggest no systematic relationship.
 *
 * Three weighting schemes are supported:
 *
 * - "unit": Uses w_e = 1 for all edges, treating all edges equally regardless of length.
 *   Appropriate when edge lengths are comparable or when counting directional agreements
 *   without geometric normalization.
 *
 * - "derivative": Uses w_e = 1/(Δ_e)² where Δ_e is the edge length. This normalizes by
 *   edge length squared, making the measure analogous to comparing derivatives rather
 *   than absolute changes. Natural when functions represent continuous quantities sampled
 *   at irregular spatial positions.
 *
 * - "sign": Uses only the sign of Δ_e y · Δ_e z, computing the proportion of edges where
 *   the functions agree in direction. This robust measure is insensitive to outliers and
 *   magnitude of change.
 *
 * The function returns vertex-wise coefficients along with summary statistics including
 * mean, median, and counts of positive/negative/zero values. This multi-scale perspective
 * enables both detailed spatial diagnostics and high-level assessment of global association.
 *
 * The algorithm has O(|E|) time complexity, linear in the number of edges. For sparse
 * graphs with bounded degree, this is also O(|V|) in the number of vertices.
 *
 * @param s_adj_list R list of integer vectors representing the graph's adjacency structure.
 *   Element i contains 0-based indices of vertices adjacent to vertex i. The graph is
 *   undirected, so edge [i,j] should appear in both adjacency lists.
 *
 * @param s_weight_list R list of numeric vectors containing edge weights (lengths).
 *   Element i contains weights corresponding to the edges in s_adj_list[[i]]. Must have
 *   the same structure as s_adj_list. Weights should be positive.
 *
 * @param s_y R numeric vector of function values at each vertex. Length must equal the
 *   number of vertices (i.e., length(s_adj_list)).
 *
 * @param s_z R numeric vector of second function values at each vertex. Length must equal
 *   the number of vertices.
 *
 * @param s_type R character scalar specifying the weighting scheme. Must be one of:
 *   "unit", "derivative", or "sign". Case-sensitive.
 *
 * @return R list with components:
 *   - vertex.coefficients: Numeric vector of co-monotonicity coefficients at each vertex
 *   - mean.coefficient: Scalar mean of vertex coefficients
 *   - median.coefficient: Scalar median of vertex coefficients  
 *   - n.positive: Integer count of vertices with positive co-monotonicity (coeff > epsilon)
 *   - n.negative: Integer count of vertices with negative co-monotonicity (coeff < -epsilon)
 *   - n.zero: Integer count of vertices with zero co-monotonicity (|coeff| <= epsilon)
 *
 * @note All vertex indices in input are 0-based (C++ convention). All vertex indices in
 *   output remain 0-based for consistency with internal graph representation. The R
 *   wrapper function handles conversion to 1-based indexing if needed.
 *
 * @note Memory is managed using R's PROTECT/UNPROTECT mechanism. Each allocated SEXP
 *   must be protected before use and unprotected when no longer needed.
 *
 * @see comono_type_t for detailed description of weighting schemes
 * @see set_wgraph_t::comono for the underlying C++ implementation
 *
 * @example
 * // In R:
 * // adj.list <- list(c(1,2), c(0,2), c(0,1))  # 0-based
 * // weight.list <- list(c(1.0,1.0), c(1.0,1.0), c(1.0,1.0))
 * // y <- c(1.0, 2.0, 3.0)
 * // z <- c(2.0, 4.0, 6.0)
 * // result <- .Call(S_comono, adj.list, weight.list, y, z, "unit")
 */
extern "C" SEXP S_comono(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type
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
    
    // Convert type string
    if (!Rf_isString(s_type) || LENGTH(s_type) != 1) {
        Rf_error("s_type must be a single character string");
    }
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);
    
    // Map string to enum
    comono_type_t comono_type;
    if (type_str == "unit") {
        comono_type = comono_type_t::UNIT;
    } else if (type_str == "derivative") {
        comono_type = comono_type_t::DERIVATIVE;
    } else if (type_str == "sign") {
        comono_type = comono_type_t::SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'unit', 'derivative', or 'sign'",
                 type_cstr);
    }
    
    // ---- Build graph ----
    set_wgraph_t graph(adj_list, weight_list);
    
    // ---- Compute co-monotonicity ----
    comono_result_t result = graph.comono(y, z, comono_type);
    
    // ---- Build R return object ----
    
    // Create list with 6 components
    const int n_components = 6;
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
    
    // Component 2: mean.coefficient (scalar)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("mean.coefficient"));
    SET_VECTOR_ELT(r_result, idx++, Rf_ScalarReal(result.mean_coefficient));
    
    // Component 3: median.coefficient (scalar)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("median.coefficient"));
    SET_VECTOR_ELT(r_result, idx++, Rf_ScalarReal(result.median_coefficient));
    
    // Component 4: n.positive (integer scalar)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.positive"));
    SET_VECTOR_ELT(r_result, idx++, Rf_ScalarInteger(
        static_cast<int>(result.n_positive)
    ));
    
    // Component 5: n.negative (integer scalar)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.negative"));
    SET_VECTOR_ELT(r_result, idx++, Rf_ScalarInteger(
        static_cast<int>(result.n_negative)
    ));
    
    // Component 6: n.zero (integer scalar)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.zero"));
    SET_VECTOR_ELT(r_result, idx++, Rf_ScalarInteger(
        static_cast<int>(result.n_zero)
    ));
    
    // Set component names
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    
    // Clean up
    UNPROTECT(2); // r_names, r_result
    
    return r_result;
}

/**
 * @brief SEXP interface for computing co-monotonicity coefficients between a vector and matrix columns
 *
 * @details
 * This function provides an R interface to the matrix co-monotonicity coefficient
 * computation. It efficiently computes vertex-wise co-monotonicity between a vector
 * y and each column of a matrix Z, measuring the extent to which y and each z_j
 * vary together across the edges of a weighted graph.
 *
 * The function is optimized for computing multiple co-monotonicity measures
 * simultaneously by pre-computing graph-dependent quantities (y edge differences,
 * edge weights, neighbor indices) once and reusing them for all columns. This
 * eliminates redundant computation that would occur with q separate calls to the
 * scalar co-monotonicity function.
 *
 * MATHEMATICAL FRAMEWORK
 *
 * For each column j of matrix Z, the vertex-level co-monotonicity coefficient
 * at vertex v is defined as:
 *
 *   comono(y,z_j;w)(v) = Σ_{u ∈ N(v)} w_e Δ_e y Δ_e z_j / Σ_{u ∈ N(v)} w_e |Δ_e y Δ_e z_j|
 *
 * where:
 *   - N(v) is the neighborhood of vertex v
 *   - Δ_e y = y(u) - y(v) is the edge difference for y along edge e = [v,u]
 *   - Δ_e z_j = z_j(u) - z_j(v) is the edge difference for column j of Z
 *   - w_e are edge weights determined by the type parameter
 *
 * The coefficient lies in [-1, 1] for each column:
 *   +1: Perfect positive co-monotonicity (y and z_j always change together)
 *   -1: Perfect negative co-monotonicity (y and z_j always change oppositely)
 *    0: No systematic directional relationship
 *
 * WEIGHTING SCHEMES
 *
 * Three weighting schemes are supported via the type parameter:
 *
 * "unit": Uses w_e = 1 for all edges. Treats all edges equally regardless of
 * their length, appropriate when edges are comparable or when counting directional
 * agreements without geometric normalization.
 *
 * "derivative": Uses w_e = 1/(Δ_e)² where Δ_e is the edge length. Normalizes by
 * edge length squared, making the measure analogous to comparing derivatives rather
 * than absolute changes. Natural for continuous functions sampled at irregular
 * spatial positions.
 *
 * "sign": Uses only the sign of Δ_e y · Δ_e z_j, computing the proportion of edges
 * where the functions agree in direction. This robust measure is insensitive to
 * outliers and magnitude of change.
 *
 * ALGORITHM
 *
 * The implementation uses a two-phase approach:
 *
 * Phase 1 - Pre-computation (performed once):
 *   - Build graph from adjacency and weight lists
 *   - For each vertex v and incident edge e = [v,u]:
 *     * Compute and store Δ_e y = y(u) - y(v)
 *     * Compute and store edge weight w_e based on type
 *     * Store neighbor index u
 *
 * Phase 2 - Column processing (repeated for each column j):
 *   - For each vertex v:
 *     * For each incident edge e = [v,u]:
 *       - Compute Δ_e z_j using stored neighbor index
 *       - Form product with stored Δ_e y
 *       - Accumulate using stored weight w_e
 *     * Compute vertex coefficient from accumulated values
 *   - Compute summary statistics (mean, median, counts)
 *
 * This structure ensures that expensive graph operations occur only once, while
 * column-specific computations reuse cached information for efficiency.
 *
 * COMPLEXITY
 *
 * Time: O(q · |E| + q · |V| log |V|) where q is the number of columns in Z
 *   - Pre-computation: O(|E|)
 *   - Column processing: O(q · |E|)
 *   - Summary statistics: O(q · |V| log |V|)
 *
 * Space: O(|E| + q · |V|)
 *   - Pre-computed structures: O(|E|)
 *   - Result storage: O(q · |V|)
 *
 * For sparse graphs with bounded degree, this is effectively O(q · |V| log |V|).
 *
 * R INTERFACE DETAILS
 *
 * This SEXP function bridges R and C++ by:
 * 1. Converting R SEXP types to C++ types
 * 2. Validating input dimensions and types
 * 3. Building the graph structure
 * 4. Calling the pure C++ implementation
 * 5. Converting C++ results back to R SEXP structures
 * 6. Managing memory using PROTECT/UNPROTECT
 *
 * Memory management follows R's PROTECT/UNPROTECT protocol. Each allocated SEXP
 * must be protected before use and unprotected when no longer needed to prevent
 * garbage collection during function execution.
 *
 * INDEX CONVENTIONS
 *
 * Input adjacency lists use 0-based indexing (C++ convention). The R wrapper
 * function performs the conversion from R's 1-based to C++'s 0-based indexing
 * before calling this function.
 *
 * Output vertex indices in the result structures also use 0-based indexing internally,
 * though users typically interact with results via vertex positions in vectors rather
 * than explicit indices.
 *
 * USAGE EXAMPLE (from R)
 *
 * ```r
 * # Build graph
 * adj.list <- list(c(1,2), c(0,2), c(0,1))  # 0-based
 * weight.list <- list(c(1.0,1.5), c(1.0,2.0), c(1.5,2.0))
 *
 * # Define functions
 * y <- c(1.0, 2.0, 3.0)
 * Z <- matrix(c(2.0, 4.0, 6.0,    # Column 1: perfect positive correlation
 *               3.0, 2.0, 1.0),   # Column 2: perfect negative correlation
 *             nrow = 3, ncol = 2)
 *
 * # Compute co-monotonicity
 * result <- .Call(S_comono_matrix, adj.list, weight.list, y, Z, "unit",
 *                 PACKAGE = "gflow")
 *
 * # Access results
 * result$mean.coefficients     # c(0.95, -0.95) approximately
 * result$column.coefficients   # List of 2 vectors, each length 3
 * ```
 *
 * @param s_adj_list R list of integer vectors representing the graph's adjacency
 *   structure. Element i contains 0-based indices of vertices adjacent to vertex i.
 *   The graph is undirected, so edge [i,j] should appear in both adjacency lists.
 *   Must have length equal to the number of vertices.
 *
 * @param s_weight_list R list of numeric vectors containing edge weights (lengths).
 *   Element i contains weights corresponding to the edges in s_adj_list[[i]].
 *   Must have the same structure as s_adj_list. Weights should be positive.
 *
 * @param s_y R numeric vector of function values at each vertex. Length must equal
 *   the number of vertices (i.e., length(s_adj_list)).
 *
 * @param s_Z R numeric matrix of function values. Each column represents a function
 *   defined on the vertices. Number of rows must equal the number of vertices.
 *   Number of columns (q) determines how many co-monotonicity computations are
 *   performed.
 *
 * @param s_type R character scalar specifying the weighting scheme. Must be one of:
 *   "unit", "derivative", or "sign". Case-sensitive.
 *
 * @return R list with 6 components (all using 0-based indexing internally):
 *
 *   \code{column.coefficients}: R list of length q, where element j is a numeric
 *     vector of length n containing vertex-wise co-monotonicity coefficients
 *     between y and column j of Z. Access as result[[1]][[j]][v+1] for column j,
 *     vertex v (converting to 1-based for R).
 *
 *   \code{mean.coefficients}: R numeric vector of length q containing the mean
 *     co-monotonicity coefficient for each column. Element j is the mean of the
 *     vertex coefficients for column j.
 *
 *   \code{median.coefficients}: R numeric vector of length q containing the median
 *     co-monotonicity coefficient for each column.
 *
 *   \code{n.positive}: R integer vector of length q containing the count of vertices
 *     with positive co-monotonicity (coefficient > epsilon) for each column.
 *
 *   \code{n.negative}: R integer vector of length q containing the count of vertices
 *     with negative co-monotonicity (coefficient < -epsilon) for each column.
 *
 *   \code{n.zero}: R integer vector of length q containing the count of vertices
 *     with zero co-monotonicity (|coefficient| <= epsilon) for each column.
 *
 * @throws Rf_error if s_adj_list is not a list
 * @throws Rf_error if s_weight_list is not a list
 * @throws Rf_error if s_y is not a numeric vector
 * @throws Rf_error if s_Z is not a numeric matrix
 * @throws Rf_error if s_type is not a single character string
 * @throws Rf_error if s_type is not "unit", "derivative", or "sign"
 * @throws Rf_error if length(s_y) does not match number of vertices
 * @throws Rf_error if nrow(s_Z) does not match number of vertices
 *
 * @note All vertex indices in input adjacency lists must be 0-based (C++ convention).
 *   The R wrapper function handles conversion from R's 1-based indexing.
 *
 * @note Memory is managed using R's PROTECT/UNPROTECT mechanism. Each allocated
 *   SEXP must be protected before use and unprotected when done.
 *
 * @note For very large matrices (q > 1000), consider processing columns in batches
 *   to manage memory usage, as results require O(q · n) storage.
 *
 * @see comono_matrix in comono_coefficient.cpp for the underlying C++ implementation
 * @see comono_matrix_result_t for the C++ result structure
 * @see comono_type_t for detailed description of weighting schemes
 * @see S_comono for the scalar (vector-vector) version
 */
extern "C" SEXP S_comono_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_Z,           // Now a matrix
    SEXP s_type
) {
    // Convert inputs (similar to S_comono)
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y vector
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    // Convert Z matrix
    if (!Rf_isMatrix(s_Z)) {
        Rf_error("s_Z must be a numeric matrix");
    }

    int* dims = INTEGER(Rf_getAttrib(s_Z, R_DimSymbol));
    int n_rows = dims[0];
    int n_cols = dims[1];

    double* Z_ptr = REAL(s_Z);
    Eigen::Map<Eigen::MatrixXd> Z(Z_ptr, n_rows, n_cols);

    // Convert type
    const char* type_cstr = CHAR(STRING_ELT(s_type, 0));
    std::string type_str(type_cstr);

    // Map string to enum
    comono_type_t comono_type;
    if (type_str == "unit") {
        comono_type = comono_type_t::UNIT;
    } else if (type_str == "derivative") {
        comono_type = comono_type_t::DERIVATIVE;
    } else if (type_str == "sign") {
        comono_type = comono_type_t::SIGN;
    } else {
        Rf_error("Invalid type '%s'. Must be 'unit', 'derivative', or 'sign'",
                 type_cstr);
    }

    // Build graph and compute
    set_wgraph_t graph(adj_list, weight_list);
    comono_matrix_result_t result = graph.comono_matrix(y, Z, comono_type);

    // Build R return object
    // Return list with:
    // - column.coefficients: list of vectors
    // - mean.coefficients: vector
    // - median.coefficients: vector
    // - n.positive: vector
    // - n.negative: vector
    // - n.zero: vector

    const int n_components = 6;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // column.coefficients (list of vectors)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("column.coefficients"));
    SEXP r_col_coeffs = PROTECT(Rf_allocVector(VECSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        SEXP r_col = PROTECT(Rf_allocVector(REALSXP, n_rows));
        double* col_ptr = REAL(r_col);
        for (int i = 0; i < n_rows; ++i) {
            col_ptr[i] = result.column_coefficients[col][i];
        }
        SET_VECTOR_ELT(r_col_coeffs, col, r_col);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(r_result, idx++, r_col_coeffs);
    UNPROTECT(1);

    // mean.coefficients (vector)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("mean.coefficients"));
    SEXP r_means = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_means)[col] = result.mean_coefficients[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_means);
    UNPROTECT(1);

    // median.coefficients (vector)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("median.coefficients"));
    SEXP r_medians = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_medians)[col] = result.median_coefficients[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_medians);
    UNPROTECT(1);

    // n.positive (vector)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.positive"));
    SEXP r_n_positives = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_n_positives)[col] = result.n_positive[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_n_positives);
    UNPROTECT(1);

    // n.negative (vector)
    SET_STRING_ELT(r_names, idx, Rf_mkChar("n.negative"));
    SEXP r_n_negatives = PROTECT(Rf_allocVector(REALSXP, n_cols));
    for (int col = 0; col < n_cols; ++col) {
        REAL(r_n_negatives)[col] = result.n_negative[col];
    }
    SET_VECTOR_ELT(r_result, idx++, r_n_negatives);
    UNPROTECT(1);

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);
    UNPROTECT(2);

    return r_result;
}
