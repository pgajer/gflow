/**
 * @brief Hardened rchk-safe drop-in of create_iknn_graph
 *
 * Key patterns:
 * - Input matrix coercion via PROTECT_WITH_INDEX/REPROTECT; copy/reads while protected.
 * - Defensive scalar parsing (Rf_asInteger), early validation with balanced UNPROTECT(1) on RX.
 * - Local PROTECT/UNPROTECT(1) blocks for all temporaries; no running protection counters.
 * - Keep `result` and `names` protected until the tail; end with UNPROTECT(2).
 */
#include <vector>
#include <R.h>
#include <Rinternals.h>

// Forward declarations (assumed available at link time)
extern SEXP S_kNN(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);
extern SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>&);

struct iknn_graph_t {
    std::vector<std::vector<int>> adjacency_list;
    std::vector<std::vector<double>> edge_weights;
    int k;
    double average_degree;
    int n_edges;
    int n_vertices;
};

extern "C" {

SEXP create_iknn_graph(SEXP X_sexp, SEXP k_sexp) {
    // --- Coerce X to REAL matrix (indexed protect) ---
    PROTECT_INDEX px;
    SEXP RX = X_sexp;
    PROTECT_WITH_INDEX(RX, &px);
    if (TYPEOF(RX) != REALSXP) {
        REPROTECT(RX = Rf_coerceVector(RX, REALSXP), px);
    }

    // Validate dimensions
    SEXP dims = Rf_getAttrib(RX, R_DimSymbol);
    if (Rf_isNull(dims) || LENGTH(dims) != 2 || TYPEOF(dims) != INTSXP) {
        UNPROTECT(1); // RX
        Rf_error("create_iknn_graph: 'X' must be a numeric matrix");
    }
    const int n_points = INTEGER(dims)[0];
    const int n_dims   = INTEGER(dims)[1];
    if (n_points <= 1 || n_dims <= 0) {
        UNPROTECT(1); // RX
        Rf_error("create_iknn_graph: invalid matrix dimensions");
    }

    // Defensive scalar parse
    const int k = Rf_asInteger(k_sexp);
    if (k < 1) {
        UNPROTECT(1); // RX
        Rf_error("create_iknn_graph: 'k' must be >= 1");
    }
    if (k >= n_points) {
        UNPROTECT(1); // RX
        Rf_error("create_iknn_graph: 'k' must be < number of points");
    }

    // --- Call S_kNN (protect inputs; unprotect locals with a fixed count) ---
    SEXP nn_idx = R_NilValue;
    SEXP nn_dists = R_NilValue;
    {
        SEXP s_k          = PROTECT(Rf_ScalarInteger(k + 1)); // include self
        SEXP s_search     = PROTECT(Rf_ScalarInteger(1));     // standard search
        SEXP s_bucketSize = PROTECT(Rf_ScalarInteger(10));    // default bucket size
        SEXP s_splitRule  = PROTECT(Rf_ScalarInteger(5));     // ANN_KD_SUGGEST (as in prior code)

        SEXP knn_result = PROTECT(S_kNN(RX, s_k, s_search, s_bucketSize, s_splitRule));
        // Extract needed components (assumed positions 0 and 1)
        nn_idx   = VECTOR_ELT(knn_result, 0);
        nn_dists = VECTOR_ELT(knn_result, 1);
        // Release local temporaries (s_k, s_search, s_bucketSize, s_splitRule, knn_result)
        UNPROTECT(5);
    }

    // Convert RX (column-major) to row-major if needed for further processing
    // (We use RX directly to read neighbor indices/distances; data copy is not strictly necessary here.
    //  Keep as minimal as possible.)
    const double* X_ptr = REAL(RX);

    // Build ikNN graph from kNN outputs
    iknn_graph_t graph;
    graph.n_vertices = n_points;
    graph.k = k;
    graph.adjacency_list.assign(n_points, std::vector<int>());
    graph.edge_weights.assign(n_points, std::vector<double>());

    // Accessors from kNN
    const int*    idx_ptr  = INTEGER(nn_idx);
    const double* dist_ptr = REAL(nn_dists);

    for (int i = 0; i < n_points; ++i) {
        for (int j = 1; j <= k; ++j) { // skip j=0 (self)
            const int neighbor = idx_ptr[i + n_points * j] - 1; // 1-based -> 0-based
            const double distance = dist_ptr[i + n_points * j];
            if (neighbor >= 0 && neighbor < n_points && neighbor != i) {
                graph.adjacency_list[i].push_back(neighbor);
                graph.edge_weights[i].push_back(1.0 / (1.0 + distance)); // simple weight
            }
        }
    }

    // Compute stats
    graph.n_edges = 0;
    for (int i = 0; i < n_points; ++i) {
        graph.n_edges += (int) graph.adjacency_list[i].size();
    }
    graph.average_degree = (double) graph.n_edges / (double) n_points;

    // --- Assemble result (container-first). Keep result & names protected until tail. ---
    const int N = 6;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));

    // 0: adjacency_list
    {
        SEXP s_adj = PROTECT(convert_vector_vector_int_to_R(graph.adjacency_list));
        SET_VECTOR_ELT(result, 0, s_adj);
        UNPROTECT(1);
    }
    // 1: edge_weights
    {
        SEXP s_w = PROTECT(convert_vector_vector_double_to_R(graph.edge_weights));
        SET_VECTOR_ELT(result, 1, s_w);
        UNPROTECT(1);
    }
    // 2: k
    {
        SEXP s_kval = PROTECT(Rf_ScalarInteger(graph.k));
        SET_VECTOR_ELT(result, 2, s_kval);
        UNPROTECT(1);
    }
    // 3: average_degree
    {
        SEXP s_avg = PROTECT(Rf_ScalarReal(graph.average_degree));
        SET_VECTOR_ELT(result, 3, s_avg);
        UNPROTECT(1);
    }
    // 4: n_edges
    {
        SEXP s_ne = PROTECT(Rf_ScalarInteger(graph.n_edges));
        SET_VECTOR_ELT(result, 4, s_ne);
        UNPROTECT(1);
    }
    // 5: n_vertices
    {
        SEXP s_nv = PROTECT(Rf_ScalarInteger(graph.n_vertices));
        SET_VECTOR_ELT(result, 5, s_nv);
        UNPROTECT(1);
    }

    // names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("adjacency_list"));
    SET_STRING_ELT(names, 1, Rf_mkChar("edge_weights"));
    SET_STRING_ELT(names, 2, Rf_mkChar("k"));
    SET_STRING_ELT(names, 3, Rf_mkChar("average_degree"));
    SET_STRING_ELT(names, 4, Rf_mkChar("n_edges"));
    SET_STRING_ELT(names, 5, Rf_mkChar("n_vertices"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // Release RX (input) and then result+names at tail
    UNPROTECT(1); // RX
    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
