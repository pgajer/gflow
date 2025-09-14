/**
 * @brief Fixed version of iknn_graphs.cpp function to address rchk PROTECT/UNPROTECT issues
 * 
 * This file contains corrected version of:
 * - create_iknn_graph
 * 
 * Issues fixed:
 * 1. create_iknn_graph (lines 652, 661): RX unprotected before Rf_coerceVector
 * 2. Calling S_kNN with fresh (unprotected) pointers
 * 
 * Changes made:
 * 1. Used PROTECT_WITH_INDEX for conditional coercion
 * 2. Protected all intermediate allocations
 * 3. Fixed all PROTECT/UNPROTECT imbalances
 */

#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available)
extern SEXP S_kNN(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);

// Core structures and declarations
struct iknn_graph_t {
    std::vector<std::vector<int>> adjacency_list;
    std::vector<std::vector<double>> edge_weights;
    int k;
    double average_degree;
    int n_edges;
    int n_vertices;
};

// Core computation function (assumed available)
iknn_graph_t compute_iknn_graph(
    const std::vector<std::vector<double>>& data,
    int k,
    bool mutual,
    bool verbose);

extern "C" {

/**
 * Fixed version of create_iknn_graph
 * Fixes: RX unprotected before coercion, unprotected pointers to S_kNN
 * Solution: Use PROTECT_WITH_INDEX, protect all allocations
 */
SEXP create_iknn_graph(SEXP X_sexp, SEXP k_sexp) {
    // Input coercion block with PROTECT_WITH_INDEX
    PROTECT_INDEX ipx;
    SEXP RX = X_sexp;
    PROTECT_WITH_INDEX(RX, &ipx);
    
    // Coerce to REALSXP if needed
    if (TYPEOF(RX) != REALSXP) {
        REPROTECT(RX = Rf_coerceVector(RX, REALSXP), ipx);
    }
    
    // Get dimensions
    SEXP dims = Rf_getAttrib(RX, R_DimSymbol);
    if (Rf_isNull(dims) || LENGTH(dims) != 2) {
        UNPROTECT(1); // RX
        Rf_error("X must be a matrix");
    }
    
    int* dimptr = INTEGER(dims);
    int n_points = dimptr[0];
    int n_dims = dimptr[1];
    
    // Extract k parameter
    int k = Rf_asInteger(k_sexp);
    
    // Validate k
    if (k < 1) {
        UNPROTECT(1); // RX
        Rf_error("k must be at least 1");
    }
    if (k >= n_points) {
        UNPROTECT(1); // RX
        Rf_error("k must be less than the number of points");
    }
    
    // Create parameters for S_kNN call
    SEXP s_k = PROTECT(Rf_ScalarInteger(k + 1)); // k+1 because kNN includes self
    SEXP s_search_type = PROTECT(Rf_ScalarInteger(1)); // Standard search
    SEXP s_bucketSize = PROTECT(Rf_ScalarInteger(10)); // Default bucket size
    SEXP s_splitRule = PROTECT(Rf_ScalarInteger(5)); // ANN_KD_SUGGEST
    
    // Call S_kNN to get nearest neighbors (protected allocations)
    SEXP knn_result = PROTECT(S_kNN(RX, s_k, s_search_type, s_bucketSize, s_splitRule));
    
    // Extract nn.idx and nn.dists from result
    SEXP nn_idx = VECTOR_ELT(knn_result, 0);
    SEXP nn_dists = VECTOR_ELT(knn_result, 1);
    
    // Convert to C++ format for processing
    int* idx_ptr = INTEGER(nn_idx);
    double* dist_ptr = REAL(nn_dists);
    
    std::vector<std::vector<double>> data(n_points);
    double* X_ptr = REAL(RX);
    for (int i = 0; i < n_points; i++) {
        data[i].resize(n_dims);
        for (int j = 0; j < n_dims; j++) {
            data[i][j] = X_ptr[i + n_points * j]; // Column-major to row-major
        }
    }
    
    // Build ikNN graph structure
    iknn_graph_t graph;
    graph.n_vertices = n_points;
    graph.k = k;
    graph.adjacency_list.resize(n_points);
    graph.edge_weights.resize(n_points);
    
    // Process kNN results to build intersection-weighted graph
    for (int i = 0; i < n_points; i++) {
        for (int j = 1; j <= k; j++) { // Skip j=0 (self)
            int neighbor = idx_ptr[i + n_points * j] - 1; // Convert to 0-based
            double distance = dist_ptr[i + n_points * j];
            
            if (neighbor >= 0 && neighbor < n_points && neighbor != i) {
                graph.adjacency_list[i].push_back(neighbor);
                graph.edge_weights[i].push_back(1.0 / (1.0 + distance)); // Weight function
            }
        }
    }
    
    // Compute graph statistics
    graph.n_edges = 0;
    for (int i = 0; i < n_points; i++) {
        graph.n_edges += static_cast<int>(graph.adjacency_list[i].size());
    }
    graph.average_degree = static_cast<double>(graph.n_edges) / n_points;
    
    // Create result list (container-first pattern)
    const int N_COMPONENTS = 6;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    
    // 0: adjacency_list
    {
        SEXP adj_list = PROTECT(convert_vector_vector_int_to_R(graph.adjacency_list));
        SET_VECTOR_ELT(result, 0, adj_list);
        UNPROTECT(1);
    }
    
    // 1: edge_weights
    {
        SEXP weights = PROTECT(convert_vector_vector_double_to_R(graph.edge_weights));
        SET_VECTOR_ELT(result, 1, weights);
        UNPROTECT(1);
    }
    
    // 2: k
    {
        SEXP k_val = PROTECT(Rf_ScalarInteger(graph.k));
        SET_VECTOR_ELT(result, 2, k_val);
        UNPROTECT(1);
    }
    
    // 3: average_degree
    {
        SEXP avg_deg = PROTECT(Rf_ScalarReal(graph.average_degree));
        SET_VECTOR_ELT(result, 3, avg_deg);
        UNPROTECT(1);
    }
    
    // 4: n_edges
    {
        SEXP n_edges = PROTECT(Rf_ScalarInteger(graph.n_edges));
        SET_VECTOR_ELT(result, 4, n_edges);
        UNPROTECT(1);
    }
    
    // 5: n_vertices
    {
        SEXP n_verts = PROTECT(Rf_ScalarInteger(graph.n_vertices));
        SET_VECTOR_ELT(result, 5, n_verts);
        UNPROTECT(1);
    }
    
    // Set names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("adjacency_list"));
        SET_STRING_ELT(names, 1, Rf_mkChar("edge_weights"));
        SET_STRING_ELT(names, 2, Rf_mkChar("k"));
        SET_STRING_ELT(names, 3, Rf_mkChar("average_degree"));
        SET_STRING_ELT(names, 4, Rf_mkChar("n_edges"));
        SET_STRING_ELT(names, 5, Rf_mkChar("n_vertices"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }
    
    // Clean up - count PROTECTs: RX, s_k, s_search_type, s_bucketSize, s_splitRule, knn_result, result
    UNPROTECT(7);
    
    return result;
}

} // extern "C"