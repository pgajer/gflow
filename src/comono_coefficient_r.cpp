#include "set_wgraph.hpp"
#include "lcor.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>    // For std::accumulate

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


#if 0
/**
 * @brief Efficiently compute co-monotonicity for multiple z vectors
 *
 * @details
 * This function is optimized for permutation testing. It pre-computes
 * all y-dependent quantities once and reuses them for multiple z vectors.
 * This is much faster than calling comono() repeatedly.
 *
 * @param y Vector of response values
 * @param Z_batch Matrix where each column is a z vector to test
 * @param type Co-monotonicity type
 *
 * @return Matrix where column j contains vertex coefficients for Z_batch[,j]
 */
Eigen::MatrixXd set_wgraph_t::comono_batch(
    const std::vector<double>& y,
    const Eigen::MatrixXd& Z_batch,
    comono_type_t type = comono_type_t::UNIT
) const {
    const size_t n_vertices = num_vertices();
    const size_t n_vectors = Z_batch.cols();

    // Pre-compute y edge differences and weights (reused for all z)
    std::vector<std::vector<double>> y_edge_diffs(n_vertices);
    std::vector<std::vector<double>> edge_weights(n_vertices);
    std::vector<std::vector<size_t>> neighbor_indices(n_vertices);

    // ... [same pre-computation as in comono_matrix] ...

    // Compute for each z vector
    Eigen::MatrixXd result(n_vertices, n_vectors);

    #pragma omp parallel for if(n_vectors > 10)
    for (size_t col = 0; col < n_vectors; ++col) {
        for (size_t v = 0; v < n_vertices; ++v) {
            // ... [same per-vertex computation] ...
            result(v, col) = vertex_coeff;
        }
    }

    return result;
}
#endif


/**
 * R interface for correlation-type co-monotonicity
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency list)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector (response values)
 * @param s_z R numeric vector (feature values)
 * @param s_type R string ("unit" or "derivative")
 * @return R list with vertex_coefficients, mean, median, counts
 */
extern "C" SEXP S_comono_cor(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_type
) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y and z
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    // Parse type string
    const char* type_str = CHAR(STRING_ELT(s_type, 0));
    comono_type_t type;
    if (std::string(type_str) == "unit") {
        type = comono_type_t::UNIT;
    } else if (std::string(type_str) == "derivative") {
        type = comono_type_t::DERIVATIVE;
    } else {
        Rf_error("Unknown type '%s'. Must be 'unit' or 'derivative'", type_str);
    }

    // Build graph and compute co-monotonicity
    set_wgraph_t graph(adj_list, weight_list);
    comono_result_t result = graph.comono_cor(y, z, type);

    // Build R return list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 5));

    // vertex.coefficients
    SEXP r_coeffs = PROTECT(Rf_allocVector(REALSXP, result.vertex_coefficients.size()));
    std::copy(result.vertex_coefficients.begin(), result.vertex_coefficients.end(), REAL(r_coeffs));
    SET_VECTOR_ELT(r_result, 0, r_coeffs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("vertex.coefficients"));

    // mean.coefficient
    SEXP r_mean = PROTECT(Rf_ScalarReal(result.mean_coefficient));
    SET_VECTOR_ELT(r_result, 1, r_mean);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("mean.coefficient"));

    // median.coefficient
    SEXP r_median = PROTECT(Rf_ScalarReal(result.median_coefficient));
    SET_VECTOR_ELT(r_result, 2, r_median);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("median.coefficient"));

    // n.positive, n.negative, n.zero
    SEXP r_counts = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(r_counts)[0] = static_cast<int>(result.n_positive);
    INTEGER(r_counts)[1] = static_cast<int>(result.n_negative);
    INTEGER(r_counts)[2] = static_cast<int>(result.n_zero);
    SET_VECTOR_ELT(r_result, 3, r_counts);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("counts"));

    // Add count names
    SEXP r_count_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(r_count_names, 0, Rf_mkChar("n.positive"));
    SET_STRING_ELT(r_count_names, 1, Rf_mkChar("n.negative"));
    SET_STRING_ELT(r_count_names, 2, Rf_mkChar("n.zero"));
    Rf_setAttrib(r_counts, R_NamesSymbol, r_count_names);

    // Set list names
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(7);
    return r_result;
}

/**
 * R interface for proportion-based co-monotonicity
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency list)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector (response values)
 * @param s_z R numeric vector (feature values)
 * @param s_tau_y R numeric scalar (threshold for |Delta_y|)
 * @param s_tau_z R numeric scalar (threshold for |Delta_z|)
 * @return R list with vertex_coefficients, mean, median, counts
 */
extern "C" SEXP S_comono_proportion(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_z,
    SEXP s_tau_y,
    SEXP s_tau_z
) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y and z
    double* y_ptr = REAL(s_y);
    size_t n_y = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n_y);

    double* z_ptr = REAL(s_z);
    size_t n_z = LENGTH(s_z);
    std::vector<double> z(z_ptr, z_ptr + n_z);

    // Extract thresholds
    double tau_y = Rf_asReal(s_tau_y);
    double tau_z = Rf_asReal(s_tau_z);

    // Build graph and compute co-monotonicity
    set_wgraph_t graph(adj_list, weight_list);
    comono_result_t result = graph.comono_proportion(y, z, tau_y, tau_z);

    // Build R return list (same structure as above)
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 5));

    SEXP r_coeffs = PROTECT(Rf_allocVector(REALSXP, result.vertex_coefficients.size()));
    std::copy(result.vertex_coefficients.begin(), result.vertex_coefficients.end(), REAL(r_coeffs));
    SET_VECTOR_ELT(r_result, 0, r_coeffs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("vertex.coefficients"));

    SEXP r_mean = PROTECT(Rf_ScalarReal(result.mean_coefficient));
    SET_VECTOR_ELT(r_result, 1, r_mean);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("mean.coefficient"));

    SEXP r_median = PROTECT(Rf_ScalarReal(result.median_coefficient));
    SET_VECTOR_ELT(r_result, 2, r_median);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("median.coefficient"));

    SEXP r_counts = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(r_counts)[0] = static_cast<int>(result.n_positive);
    INTEGER(r_counts)[1] = static_cast<int>(result.n_negative);
    INTEGER(r_counts)[2] = static_cast<int>(result.n_zero);
    SET_VECTOR_ELT(r_result, 3, r_counts);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("counts"));

    SEXP r_count_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(r_count_names, 0, Rf_mkChar("n.positive"));
    SET_STRING_ELT(r_count_names, 1, Rf_mkChar("n.negative"));
    SET_STRING_ELT(r_count_names, 2, Rf_mkChar("n.zero"));
    Rf_setAttrib(r_counts, R_NamesSymbol, r_count_names);

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(7);
    return r_result;
}

/**
 * R interface for correlation-type co-monotonicity with multiple features
 *
 * Computes co-monotonicity between response y and each column of feature matrix Z.
 * More efficient than calling S_comono_cor repeatedly because edge differences
 * for y are computed only once and reused across all features.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency list)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector (response values, length n)
 * @param s_Z R numeric matrix (feature values, n x m)
 * @param s_type R string ("unit" or "derivative")
 * @return R list with:
 *   - column.coefficients: n x m matrix of vertex-level coefficients
 *   - column.means: vector of length m with mean coefficient per feature
 *   - column.medians: vector of length m with median coefficient per feature
 *   - column.counts: m x 3 matrix with (n.positive, n.negative, n.zero) per feature
 */
extern "C" SEXP S_comono_cor_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_Z,
    SEXP s_type
) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y
    double* y_ptr = REAL(s_y);
    size_t n = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n);

    // Convert Z matrix
    // R stores matrices in column-major order
    double* Z_ptr = REAL(s_Z);
    int* Z_dims = INTEGER(Rf_getAttrib(s_Z, R_DimSymbol));
    size_t n_rows = Z_dims[0];
    size_t n_cols = Z_dims[1];

    if (n_rows != n) {
        Rf_error("Number of rows in Z (%zu) must match length of y (%zu)", n_rows, n);
    }

    // Parse type string
    const char* type_str = CHAR(STRING_ELT(s_type, 0));
    comono_type_t type;
    if (std::string(type_str) == "unit") {
        type = comono_type_t::UNIT;
    } else if (std::string(type_str) == "derivative") {
        type = comono_type_t::DERIVATIVE;
    } else {
        Rf_error("Unknown type '%s'. Must be 'unit' or 'derivative'", type_str);
    }

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Precompute y-dependent quantities for efficiency
    // For each vertex, store: (weight, delta_y, delta_y^2) for each edge
    const double MIN_DENOMINATOR = 1e-10;

    struct EdgeYInfo {
        size_t neighbor;
        double weight;
        double delta_y;
        double delta_y_sq;
    };

    std::vector<std::vector<EdgeYInfo>> y_info(n);
    std::vector<double> y_denom(n);  // sqrt(sum w (delta_y)^2) per vertex

    // Precompute for all vertices
    for (size_t v = 0; v < n; ++v) {
        double sum_y_sq = 0.0;

        for (const auto& edge_info : adj_list[v]) {
            size_t u = edge_info;
            double edge_length = weight_list[v][&edge_info - &adj_list[v][0]];

            double delta_y = y[u] - y[v];

            // Compute weight
            double weight = 1.0;
            if (type == comono_type_t::DERIVATIVE) {
                if (edge_length > 1e-10) {
                    weight = 1.0 / (edge_length * edge_length);
                } else {
                    continue;  // Skip degenerate edges
                }
            }

            EdgeYInfo info;
            info.neighbor = u;
            info.weight = weight;
            info.delta_y = delta_y;
            info.delta_y_sq = weight * delta_y * delta_y;

            y_info[v].push_back(info);
            sum_y_sq += info.delta_y_sq;
        }

        y_denom[v] = std::sqrt(sum_y_sq);
    }

    // Allocate output matrix (n x m) for vertex-level coefficients
    SEXP r_coeffs_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
    double* coeffs_ptr = REAL(r_coeffs_matrix);

    // Allocate output vectors for summary statistics
    SEXP r_means = PROTECT(Rf_allocVector(REALSXP, n_cols));
    SEXP r_medians = PROTECT(Rf_allocVector(REALSXP, n_cols));
    SEXP r_counts_matrix = PROTECT(Rf_allocMatrix(INTSXP, n_cols, 3));

    double* means_ptr = REAL(r_means);
    double* medians_ptr = REAL(r_medians);
    int* counts_ptr = INTEGER(r_counts_matrix);

    // Compute co-monotonicity for each feature (column of Z)
    for (size_t j = 0; j < n_cols; ++j) {
        // Extract column j from Z
        std::vector<double> z(n);
        for (size_t i = 0; i < n; ++i) {
            z[i] = Z_ptr[i + j * n_rows];  // Column-major indexing
        }

        // Compute coefficient for each vertex
        std::vector<double> vertex_coeffs(n);
        size_t n_pos = 0, n_neg = 0, n_zero = 0;
        const double epsilon = 1e-10;

        for (size_t v = 0; v < n; ++v) {
            double numerator = 0.0;
            double sum_z_sq = 0.0;

            // Use precomputed y info
            for (const auto& info : y_info[v]) {
                double delta_z = z[info.neighbor] - z[v];

                numerator += info.weight * info.delta_y * delta_z;
                sum_z_sq += info.weight * delta_z * delta_z;
            }

            double denom_z = std::sqrt(sum_z_sq);

            // Check numerical stability
            if (y_denom[v] > MIN_DENOMINATOR && denom_z > MIN_DENOMINATOR) {
                vertex_coeffs[v] = numerator / (y_denom[v] * denom_z);

                // Clamp to [-1, 1]
                if (vertex_coeffs[v] > 1.0) vertex_coeffs[v] = 1.0;
                else if (vertex_coeffs[v] < -1.0) vertex_coeffs[v] = -1.0;
            } else {
                vertex_coeffs[v] = 0.0;
            }

            // Store in output matrix (column-major)
            coeffs_ptr[v + j * n_rows] = vertex_coeffs[v];

            // Count signs
            if (vertex_coeffs[v] > epsilon) ++n_pos;
            else if (vertex_coeffs[v] < -epsilon) ++n_neg;
            else ++n_zero;
        }

        // Compute mean
        means_ptr[j] = std::accumulate(vertex_coeffs.begin(), vertex_coeffs.end(), 0.0) / n;

        // Compute median
        std::sort(vertex_coeffs.begin(), vertex_coeffs.end());
        if (n % 2 == 0) {
            medians_ptr[j] = 0.5 * (vertex_coeffs[n/2 - 1] + vertex_coeffs[n/2]);
        } else {
            medians_ptr[j] = vertex_coeffs[n/2];
        }

        // Store counts (row-major in counts matrix)
        counts_ptr[j + 0 * n_cols] = static_cast<int>(n_pos);
        counts_ptr[j + 1 * n_cols] = static_cast<int>(n_neg);
        counts_ptr[j + 2 * n_cols] = static_cast<int>(n_zero);
    }

    // Build return list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 4));

    SET_VECTOR_ELT(r_result, 0, r_coeffs_matrix);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("column.coefficients"));

    SET_VECTOR_ELT(r_result, 1, r_means);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("column.means"));

    SET_VECTOR_ELT(r_result, 2, r_medians);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("column.medians"));

    SET_VECTOR_ELT(r_result, 3, r_counts_matrix);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("column.counts"));

    // Add row/column names to counts matrix
    SEXP r_count_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_count_rownames = PROTECT(Rf_allocVector(STRSXP, n_cols));
    SEXP r_count_colnames = PROTECT(Rf_allocVector(STRSXP, 3));

    // Feature names as row names (if Z has colnames)
    SEXP Z_dimnames = Rf_getAttrib(s_Z, R_DimNamesSymbol);
    if (!Rf_isNull(Z_dimnames) && LENGTH(Z_dimnames) >= 2) {
        SEXP Z_colnames = VECTOR_ELT(Z_dimnames, 1);
        if (!Rf_isNull(Z_colnames)) {
            for (size_t j = 0; j < n_cols; ++j) {
                SET_STRING_ELT(r_count_rownames, j, STRING_ELT(Z_colnames, j));
            }
        }
    }

    SET_STRING_ELT(r_count_colnames, 0, Rf_mkChar("n.positive"));
    SET_STRING_ELT(r_count_colnames, 1, Rf_mkChar("n.negative"));
    SET_STRING_ELT(r_count_colnames, 2, Rf_mkChar("n.zero"));

    SET_VECTOR_ELT(r_count_dimnames, 0, r_count_rownames);
    SET_VECTOR_ELT(r_count_dimnames, 1, r_count_colnames);
    Rf_setAttrib(r_counts_matrix, R_DimNamesSymbol, r_count_dimnames);

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(9);
    return r_result;
}

/**
 * R interface for proportion-based co-monotonicity with multiple features
 *
 * Computes threshold-filtered co-monotonicity between response y and each
 * column of feature matrix Z. Accepts per-feature thresholds for adaptive
 * filtering.
 *
 * @param s_adj_list R list of integer vectors (0-based adjacency list)
 * @param s_weight_list R list of numeric vectors (edge weights)
 * @param s_y R numeric vector (response values, length n)
 * @param s_Z R numeric matrix (feature values, n x m)
 * @param s_tau_y R numeric scalar (threshold for |Delta_y|)
 * @param s_tau_z R numeric vector of length m (thresholds for |Delta_z| per feature)
 *               Can also be scalar, which is replicated for all features
 * @return R list with:
 *   - column.coefficients: n x m matrix of vertex-level coefficients
 *   - column.means: vector of length m with mean coefficient per feature
 *   - column.medians: vector of length m with median coefficient per feature
 *   - column.counts: m x 3 matrix with (n.positive, n.negative, n.zero) per feature
 */
extern "C" SEXP S_comono_proportion_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_Z,
    SEXP s_tau_y,
    SEXP s_tau_z
) {
    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert y
    double* y_ptr = REAL(s_y);
    size_t n = LENGTH(s_y);
    std::vector<double> y(y_ptr, y_ptr + n);

    // Convert Z matrix
    double* Z_ptr = REAL(s_Z);
    int* Z_dims = INTEGER(Rf_getAttrib(s_Z, R_DimSymbol));
    size_t n_rows = Z_dims[0];
    size_t n_cols = Z_dims[1];

    if (n_rows != n) {
        Rf_error("Number of rows in Z (%zu) must match length of y (%zu)", n_rows, n);
    }

    // Extract tau_y (scalar)
    double tau_y = Rf_asReal(s_tau_y);
    if (tau_y < 0.0) {
        Rf_error("tau_y must be non-negative");
    }

    // Extract tau_z (vector or scalar)
    std::vector<double> tau_z_vec(n_cols);
    size_t tau_z_len = LENGTH(s_tau_z);

    if (tau_z_len == 1) {
        // Scalar: replicate for all features
        double tau_z_scalar = Rf_asReal(s_tau_z);
        if (tau_z_scalar < 0.0) {
            Rf_error("tau_z must be non-negative");
        }
        std::fill(tau_z_vec.begin(), tau_z_vec.end(), tau_z_scalar);
    } else if (tau_z_len == n_cols) {
        // Vector: one threshold per feature
        double* tau_z_ptr = REAL(s_tau_z);
        for (size_t j = 0; j < n_cols; ++j) {
            if (tau_z_ptr[j] < 0.0) {
                Rf_error("All elements of tau_z must be non-negative");
            }
            tau_z_vec[j] = tau_z_ptr[j];
        }
    } else {
        Rf_error("tau_z must be scalar or vector of length %zu (number of features)", n_cols);
    }

    // Build graph
    set_wgraph_t graph(adj_list, weight_list);

    // Precompute y-dependent quantities
    struct EdgeInfo {
        size_t neighbor;
        double abs_delta_y;
        double delta_y;
    };

    std::vector<std::vector<EdgeInfo>> edge_info(n);
    std::vector<size_t> n_edges(n);  // Total edges per vertex

    for (size_t v = 0; v < n; ++v) {
        n_edges[v] = adj_list[v].size();

        for (size_t idx = 0; idx < adj_list[v].size(); ++idx) {
            size_t u = adj_list[v][idx];
            double delta_y = y[u] - y[v];

            EdgeInfo info;
            info.neighbor = u;
            info.delta_y = delta_y;
            info.abs_delta_y = std::abs(delta_y);

            edge_info[v].push_back(info);
        }
    }

    // Allocate output structures
    SEXP r_coeffs_matrix = PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
    SEXP r_means = PROTECT(Rf_allocVector(REALSXP, n_cols));
    SEXP r_medians = PROTECT(Rf_allocVector(REALSXP, n_cols));
    SEXP r_counts_matrix = PROTECT(Rf_allocMatrix(INTSXP, n_cols, 3));

    double* coeffs_ptr = REAL(r_coeffs_matrix);
    double* means_ptr = REAL(r_means);
    double* medians_ptr = REAL(r_medians);
    int* counts_ptr = INTEGER(r_counts_matrix);

    // Compute co-monotonicity for each feature
    for (size_t j = 0; j < n_cols; ++j) {
        // Extract column j from Z
        std::vector<double> z(n);
        for (size_t i = 0; i < n; ++i) {
            z[i] = Z_ptr[i + j * n_rows];
        }

        double tau_z = tau_z_vec[j];

        // Compute coefficient for each vertex
        std::vector<double> vertex_coeffs(n);
        size_t n_pos = 0, n_neg = 0, n_zero = 0;
        const double epsilon = 1e-10;

        for (size_t v = 0; v < n; ++v) {
            if (n_edges[v] == 0) {
                vertex_coeffs[v] = 0.0;
                coeffs_ptr[v + j * n_rows] = 0.0;
                ++n_zero;
                continue;
            }

            int agreement_count = 0;
            int disagreement_count = 0;

            for (const auto& info : edge_info[v]) {
                double delta_z = z[info.neighbor] - z[v];
                double abs_delta_z = std::abs(delta_z);

                // Check thresholds
                bool y_meaningful = info.abs_delta_y > tau_y;
                bool z_meaningful = abs_delta_z > tau_z;

                if (y_meaningful && z_meaningful) {
                    double product = info.delta_y * delta_z;

                    if (product > 0.0) {
                        ++agreement_count;
                    } else if (product < 0.0) {
                        ++disagreement_count;
                    }
                }
            }

            // Compute coefficient
            double net = static_cast<double>(agreement_count - disagreement_count);
            double total = static_cast<double>(n_edges[v]);
            vertex_coeffs[v] = net / total;

            // Store in matrix
            coeffs_ptr[v + j * n_rows] = vertex_coeffs[v];

            // Count signs
            if (vertex_coeffs[v] > epsilon) ++n_pos;
            else if (vertex_coeffs[v] < -epsilon) ++n_neg;
            else ++n_zero;
        }

        // Compute summary statistics
        means_ptr[j] = std::accumulate(vertex_coeffs.begin(), vertex_coeffs.end(), 0.0) / n;

        std::sort(vertex_coeffs.begin(), vertex_coeffs.end());
        if (n % 2 == 0) {
            medians_ptr[j] = 0.5 * (vertex_coeffs[n/2 - 1] + vertex_coeffs[n/2]);
        } else {
            medians_ptr[j] = vertex_coeffs[n/2];
        }

        counts_ptr[j + 0 * n_cols] = static_cast<int>(n_pos);
        counts_ptr[j + 1 * n_cols] = static_cast<int>(n_neg);
        counts_ptr[j + 2 * n_cols] = static_cast<int>(n_zero);
    }

    // Build return list
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 4));

    SET_VECTOR_ELT(r_result, 0, r_coeffs_matrix);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("column.coefficients"));

    SET_VECTOR_ELT(r_result, 1, r_means);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("column.means"));

    SET_VECTOR_ELT(r_result, 2, r_medians);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("column.medians"));

    SET_VECTOR_ELT(r_result, 3, r_counts_matrix);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("column.counts"));

    // Add dimension names to counts matrix
    SEXP r_count_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_count_rownames = PROTECT(Rf_allocVector(STRSXP, n_cols));
    SEXP r_count_colnames = PROTECT(Rf_allocVector(STRSXP, 3));

    SEXP Z_dimnames = Rf_getAttrib(s_Z, R_DimNamesSymbol);
    if (!Rf_isNull(Z_dimnames) && LENGTH(Z_dimnames) >= 2) {
        SEXP Z_colnames = VECTOR_ELT(Z_dimnames, 1);
        if (!Rf_isNull(Z_colnames)) {
            for (size_t j = 0; j < n_cols; ++j) {
                SET_STRING_ELT(r_count_rownames, j, STRING_ELT(Z_colnames, j));
            }
        }
    }

    SET_STRING_ELT(r_count_colnames, 0, Rf_mkChar("n.positive"));
    SET_STRING_ELT(r_count_colnames, 1, Rf_mkChar("n.negative"));
    SET_STRING_ELT(r_count_colnames, 2, Rf_mkChar("n.zero"));

    SET_VECTOR_ELT(r_count_dimnames, 0, r_count_rownames);
    SET_VECTOR_ELT(r_count_dimnames, 1, r_count_colnames);
    Rf_setAttrib(r_counts_matrix, R_DimNamesSymbol, r_count_dimnames);

    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    UNPROTECT(9);
    return r_result;
}
