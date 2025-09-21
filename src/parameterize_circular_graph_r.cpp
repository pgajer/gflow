
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
    // Convert inputs (no allocations via R API)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    const int use_edge_lengths = LOGICAL(s_use_edge_lengths)[0] ? 1 : 0;

    set_wgraph_t graph(adj_list, weight_list);
    circular_param_result_t res = graph.parameterize_circular_graph(use_edge_lengths != 0);

    // Prepare result list and names
    static const char* kNames[] = {
        "angles","eig_vec2","eig_vec3","eig_vec4","eig_vec5","eig_vec6"
    };
    enum { K = 6 };

    SEXP result = PROTECT(Rf_allocVector(VECSXP, K));
    {
        SEXP result_names = PROTECT(Rf_allocVector(STRSXP, K));
        for (int i = 0; i < K; ++i) {
            SET_STRING_ELT(result_names, i, Rf_mkChar(kNames[i]));
        }
        Rf_setAttrib(result, R_NamesSymbol, result_names);
        UNPROTECT(1); // result_names
    }

    // angles
    {
        const int n = (int) res.angles.size();
        SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
        if (n > 0) {
            double* p = REAL(v);
            std::copy(res.angles.begin(), res.angles.end(), p);
        }
        SET_VECTOR_ELT(result, 0, v);
        UNPROTECT(1); // v
    }

    // helper macro for remaining eigenvectors (handles empty => NULL)
    #define SET_EIG_VEC(slot_idx, vec_ref)                                 \
        do {                                                               \
            if ((vec_ref).empty()) {                                       \
                SET_VECTOR_ELT(result, (slot_idx), R_NilValue);            \
            } else {                                                       \
                const int n_ = (int) (vec_ref).size();                     \
                SEXP v_ = PROTECT(Rf_allocVector(REALSXP, n_));            \
                double* p_ = REAL(v_);                                     \
                std::copy((vec_ref).begin(), (vec_ref).end(), p_);         \
                SET_VECTOR_ELT(result, (slot_idx), v_);                    \
                UNPROTECT(1); /* v_ */                                     \
            }                                                              \
        } while (0)

    SET_EIG_VEC(1, res.eig_vec2);
    SET_EIG_VEC(2, res.eig_vec3);
    SET_EIG_VEC(3, res.eig_vec4);
    SET_EIG_VEC(4, res.eig_vec5);
    SET_EIG_VEC(5, res.eig_vec6);

    #undef SET_EIG_VEC

    UNPROTECT(1); // result
    return result;
}
