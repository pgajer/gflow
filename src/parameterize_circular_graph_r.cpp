
#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include <R.h>
#include <Rinternals.h>

extern "C" {
	SEXP S_parameterize_circular_graph(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_use_edge_lengths
		);
}

/**
 * @brief Parameterizes a circular graph structure using spectral methods and returns the result to R
 *
 * This function takes adjacency list and weight list from R, constructs a weighted graph,
 * applies spectral parameterization, and returns the result to R as a named list.
 *
 * The parameterization is done using eigenvectors of the graph Laplacian to embed the
 * graph in a circle. The resulting angles represent the positions of vertices on the circle.
 *
 * @param s_adj_list SEXP containing an R list of integer vectors representing adjacency lists
 * @param s_weight_list SEXP containing an R list of numeric vectors representing edge weights
 * @param s_use_edge_lengths SEXP containing a logical value indicating whether to use edge weights
 *
 * @return SEXP containing a list with three elements:
 *         - angles: numeric vector of angles (in radians) for each vertex
 *         - eig_vec2: numeric vector containing the second eigenvector
 *         - eig_vec3: numeric vector containing the third eigenvector
 *
 * @note The function assumes that the input adjacency list and weight list have the same structure
 *       and that the graph is connected.
 */
SEXP S_parameterize_circular_graph(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_use_edge_lengths
) {
    // Convert input parameters using R's C API
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    bool use_edge_lengths = (LOGICAL(s_use_edge_lengths)[0] == 1);

    set_wgraph_t graph(adj_list, weight_list);

    circular_param_result_t res = graph.parameterize_circular_graph(use_edge_lengths);

    // We will always return a 6-element list:
    // c("angles","eig_vec2","eig_vec3","eig_vec4","eig_vec5","eig_vec6")
    static const char* kNames[] = {
        "angles","eig_vec2","eig_vec3","eig_vec4","eig_vec5","eig_vec6"
    };
    constexpr int K = 6;

    int protect_count = 0;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, K));                              // (1)
    // names
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, K));                        // (2)
    for (int i = 0; i < K; ++i) SET_STRING_ELT(result_names, i, Rf_mkChar(kNames[i]));
    Rf_setAttrib(result, R_NamesSymbol, result_names);
    UNPROTECT(1); // result_names

    // Helpers that PROTECT so we can UNPROTECT at the end.
    auto create_numeric_vector = [&](const std::vector<double>& v) -> SEXP {
        SEXP r = PROTECT(Rf_allocVector(REALSXP, v.size()));
        ++protect_count;
        double* p = REAL(r);
        std::copy(v.begin(), v.end(), p);
        return r;
    };
    auto maybe_numeric_or_null = [&](const std::vector<double>& v) -> SEXP {
        if (v.empty()) return R_NilValue;
        return create_numeric_vector(v);
    };

    // Set elements
    SET_VECTOR_ELT(result, 0, create_numeric_vector(res.angles));  // angles
    SET_VECTOR_ELT(result, 1, maybe_numeric_or_null(res.eig_vec2));
    SET_VECTOR_ELT(result, 2, maybe_numeric_or_null(res.eig_vec3));
    SET_VECTOR_ELT(result, 3, maybe_numeric_or_null(res.eig_vec4));
    SET_VECTOR_ELT(result, 4, maybe_numeric_or_null(res.eig_vec5));
    SET_VECTOR_ELT(result, 5, maybe_numeric_or_null(res.eig_vec6));

    // Unprotect all vectors we created via helpers (result itself stays protected for return)
    UNPROTECT(protect_count + 1);  // +1 accounts for 'result' at the top (1)

    return result;
}
