#include "partition_graph_r.h"
#include "partition_graph.hpp"
#include <string>

#include <R.h>
#include <Rinternals.h>

// Forward declarations of helper functions (assumed to exist in package)
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

extern "C" SEXP S_partition_graph(SEXP s_adj_list, SEXP s_weight_list,
                                   SEXP s_partition, SEXP s_weight_type) {
    // Convert adjacency list
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);

    // Convert weight list
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert partition vector
    int n_vertices = LENGTH(s_partition);
    int* partition_ptr = INTEGER(s_partition);
    std::vector<int> partition(partition_ptr, partition_ptr + n_vertices);

    // Convert weight type string
    const char* weight_type_cstr = CHAR(STRING_ELT(s_weight_type, 0));
    std::string weight_type(weight_type_cstr);

    // Compute partition graph
    std::vector<std::vector<int>> adj_list_out;
    std::vector<std::vector<double>> weight_list_out;

    compute_partition_graph(adj_list, weight_list, partition, weight_type,
                           adj_list_out, weight_list_out);

    // Convert results to SEXP
    int n_cells = adj_list_out.size();

    SEXP s_adj_list_out = PROTECT(Rf_allocVector(VECSXP, n_cells));
    SEXP s_weight_list_out = PROTECT(Rf_allocVector(VECSXP, n_cells));

    for (int c = 0; c < n_cells; ++c) {
        int n_neighbors = adj_list_out[c].size();

        // Adjacency list for cell c (convert to 1-based for R)
        SEXP s_neighbors = PROTECT(Rf_allocVector(INTSXP, n_neighbors));
        int* neighbors_ptr = INTEGER(s_neighbors);
        for (int i = 0; i < n_neighbors; ++i) {
            neighbors_ptr[i] = adj_list_out[c][i] + 1; // Convert to 1-based
        }
        SET_VECTOR_ELT(s_adj_list_out, c, s_neighbors);
        UNPROTECT(1);

        // Weight list for cell c
        SEXP s_weights = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        double* weights_ptr = REAL(s_weights);
        for (int i = 0; i < n_neighbors; ++i) {
            weights_ptr[i] = weight_list_out[c][i];
        }
        SET_VECTOR_ELT(s_weight_list_out, c, s_weights);
        UNPROTECT(1);
    }

    // Create named list for return
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, 2));

    SET_STRING_ELT(s_names, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(s_names, 1, Rf_mkChar("weight_list"));
    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    SET_VECTOR_ELT(s_result, 0, s_adj_list_out);
    SET_VECTOR_ELT(s_result, 1, s_weight_list_out);

    UNPROTECT(4); // s_result, s_names, s_adj_list_out, s_weight_list_out

    return s_result;
}
