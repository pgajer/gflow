#undef match
#include <omp.h>
#define match Rf_match

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm> // Required for std::set_intersection and std::sort
#include <memory>
#include <utility>
#include <limits>
#include <omp.h>

#include "msr2.h"
#include "kNN.h"         // for struct iknn_vertex_tt and kNN_search
#include "iknn_graphs.hpp"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
    SEXP S_create_iknn_graph_mp(SEXP RX, SEXP Rk, SEXP Rn_cores);
}

std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<int>>>& adj_vect);
std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& adj_vect);
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> prune_long_edges(const std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& wgraph, int version);
std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph);
std::vector<std::vector<std::pair<int, int>>> prune_edges_with_alt_paths(const std::vector<std::vector<std::pair<int, int>>>& graph,
                                                                         std::vector<int>& long_edge_isize,
                                                                         std::vector<int>& alt_path_lengths,
                                                                         std::vector<int>& alt_path_total_isize,
                                                                         int max_alt_path_length);


/**
 * @brief Constructs a k-Nearest Neighbors graph using the Intersection Weighted Distance method with multi-processing.
 *
 * This function creates a k-NN graph based on the Intersection Weighted Distance (IWD) metric.
 * It uses OpenMP for parallel processing to improve performance on multi-core systems.
 *
 * @param s_X An R matrix (SEXP) containing the input data. Each row represents a data point,
 *           and each column represents a feature.
 * @param s_k An R integer (SEXP) specifying the number of nearest neighbors to consider for each point.
 * @param Rn_cores An R integer (SEXP) specifying the number of cores to use for parallel processing.
 *                 If set to 0 or a value greater than the available cores, it will use all available cores.
 *
 * @return SEXP A list of lists, where each sublist corresponds to a data point and contains:
 *         - indices: Integer vector of indices of the neighboring points.
 *         - common_counts: Integer vector of the number of common neighbors for each connection.
 *         - min_dists: Double vector of the minimum IWD distances for each connection.
 *
 * @note This function assumes that the input matrix s_X contains valid numeric data and that s_k is a positive integer.
 * @warning The function may consume significant memory for large datasets. Ensure sufficient system resources are available.
 *
 * @see S_create_iknn_graph_mp for a higher-level interface to this function.
 */
std::unique_ptr<std::vector<std::vector<iknn_vertex_tt>>> create_iknn_graph_mp(SEXP s_X, SEXP s_k, SEXP Rn_cores) {

#define DEBUG__create_iknn_graph_mp 0

    PROTECT(s_X = coerceVector(s_X, REALSXP));
    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    UNPROTECT(1);
    int n_vertices = dimX[0];

    PROTECT(s_k = coerceVector(s_k, INTSXP));
    int k = INTEGER(s_k)[0];
    UNPROTECT(1);

    PROTECT(Rn_cores = coerceVector(Rn_cores, INTSXP));
    int n_cores = INTEGER(Rn_cores)[0];
    UNPROTECT(1);

    // Perform kNN search
    SEXP knn_res = PROTECT(S_kNN(s_X, s_k));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));
    UNPROTECT(1);

    // Prepare result container
    auto res = std::make_unique<std::vector<std::vector<iknn_vertex_tt>>>(n_vertices);

    // Set up OpenMP
    int max_threads = omp_get_max_threads();
    int num_threads = (n_cores > 0 && n_cores <= max_threads) ? n_cores : max_threads;
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        std::vector<int> nn_i(k), nn_j(k), sorted_nn_i(k), sorted_nn_j(k), intersection;

        #pragma omp for schedule(dynamic)
        for (int pt_i = 0; pt_i < n_vertices - 1; pt_i++) {
            for (int j = 0; j < k; j++) {
                nn_i[j] = indices[pt_i + n_vertices * j];
                sorted_nn_i[j] = nn_i[j];
            }
            std::sort(sorted_nn_i.begin(), sorted_nn_i.end());

            #if DEBUG__create_iknn_graph_mp
            #pragma omp critical
            {
                Rprintf("Thread %d, pt_i %d\nnn_i: ", omp_get_thread_num(), pt_i);
                for (int j = 0; j < k; j++) Rprintf("%d ", nn_i[j]);
                Rprintf("\n");
                Rprintf("sorted_nn_i: ");
                for (int j = 0; j < k; j++) Rprintf("%d ", sorted_nn_i[j]);
                Rprintf("\n");
            }
            #endif

            for (int pt_j = pt_i + 1; pt_j < n_vertices; pt_j++) {
                for (int j = 0; j < k; j++) {
                    nn_j[j] = indices[pt_j + n_vertices * j];
                    sorted_nn_j[j] = nn_j[j];
                }
                std::sort(sorted_nn_j.begin(), sorted_nn_j.end());

                intersection.clear();
                std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(), sorted_nn_j.begin(), sorted_nn_j.end(), std::back_inserter(intersection));

                int common_count = intersection.size();

                #if DEBUG__create_iknn_graph_mp
                #pragma omp critical
                {
                    Rprintf("\nThread %d, pt_j %d\nnn_j: ", omp_get_thread_num(), pt_j);
                    for (int j = 0; j < k; j++) Rprintf("%d ", nn_j[j]);
                    Rprintf("\n");
                    Rprintf("sorted_nn_j: ");
                    for (int j = 0; j < k; j++) Rprintf("%d ", sorted_nn_j[j]);
                    Rprintf("\n");
                    Rprintf("common_count: %d\n", common_count);
                    Rprintf("intersection: ");
                    for (int j = 0; j < common_count; j++) Rprintf("%d ", intersection[j]);
                    Rprintf("\n");
                }
                #endif

                if (common_count > 0) {
                    double min_dist = std::numeric_limits<double>::max();
                    for (int x_k : intersection) {
                        int idx_i = static_cast<int>(std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin()); // The std::find function returns an iterator, not an index. We need to subtract the beginning iterator to get the index.
                        int idx_j = static_cast<int>(std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin());
                        double dist_i_k = distances[pt_i + n_vertices * idx_i];
                        double dist_j_k = distances[pt_j + n_vertices * idx_j];
                        double new_dist = dist_i_k + dist_j_k;
                        min_dist = std::min(min_dist, new_dist);

                        #if DEBUG__create_iknn_graph_mp
                        #pragma omp critical
                        {
                            Rprintf("Thread %d, pt_i %d, pt_j %d, x_k %d, idx_i %d, idx_j %d, dist_i_k %.6f, dist_j_k %.6f, min_dist %.6f\n",
                                    omp_get_thread_num(), pt_i, pt_j, x_k, idx_i, idx_j, dist_i_k, dist_j_k, min_dist);
                        }
                        #endif
                    }

                    #if DEBUG__create_iknn_graph_mp
                    #pragma omp critical
                    {
                        Rprintf("Thread %d, pt_i %d, pt_j %d: final min_dist %.6f\n",
                                omp_get_thread_num(), pt_i, pt_j, min_dist);
                    }
                    #endif

                    #pragma omp critical
                    {
                        (*res)[pt_i].emplace_back(iknn_vertex_tt{pt_j, common_count, min_dist});
                        (*res)[pt_j].emplace_back(iknn_vertex_tt{pt_i, common_count, min_dist});
                    }
                }
            }
        }
    }

    return res;
}


/**
 * @brief Creates a comprehensive k-Nearest Neighbors graph structure using the IWD method with multi-processing.
 *
 * This function serves as a higher-level interface to create_iknn_graph_mp. It processes the raw output
 * from create_iknn_graph_mp and creates a more structured representation of the k-NN graph, including
 * adjacency lists, intersection sizes, distances, and connected components.
 *
 * @param s_X An R matrix (SEXP) containing the input data. Each row represents a data point,
 *           and each column represents a feature.
 * @param s_k An R integer (SEXP) specifying the number of nearest neighbors to consider for each point.
 * @param Rn_cores An R integer (SEXP) specifying the number of cores to use for parallel processing.
 *                 If set to 0 or a value greater than the available cores, it will use all available cores.
 *
 * @return SEXP A named list containing:
 *         - adj_list: A list of integer vectors, each representing the adjacency list for a point.
 *         - isize_list: A list of integer vectors, each containing the intersection sizes for connections.
 *         - dist_list: A list of double vectors, each containing the IWD distances for connections.
 *         - conn_comps: An integer vector representing the connected components of the graph.
 *
 * @note The function converts 0-based indices from C++ to 1-based indices for R in the adjacency list.
 * @warning This function may require significant memory for large datasets. Ensure your system has sufficient resources.
 *
 * @see create_iknn_graph_mp for the underlying graph construction algorithm.
 * @see union_find for the algorithm used to identify connected components.
 */
SEXP S_create_iknn_graph_mp(SEXP s_X, SEXP s_k, SEXP Rn_cores) {

    int nprot = 0;
    PROTECT(s_X = coerceVector(s_X, REALSXP));
    if (TYPEOF(s_X) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_vertices = dimX[0];

    // Call create_iknn_graph_mp (assuming it returns SEXP)
    std::unique_ptr<std::vector<std::vector<iknn_vertex_tt>>> create_iknn_graph_res = create_iknn_graph_mp(s_X, s_k, Rn_cores);
    if (!create_iknn_graph_res) Rf_error("create_iknn_graph_mp function returned an invalid result (null pointer).");
    UNPROTECT(1); // Unprotecting s_X

    // Prepare vectors to hold the data
    std::vector<std::vector<int>> adj_vect(n_vertices);
    std::vector<std::vector<int>> intersection_size_vect(n_vertices);
    std::vector<std::vector<double>> dist_vect(n_vertices);

    // Extract data from create_iknn_graph_res
    for (int i = 0; i < n_vertices; i++) {
        for (auto nn_vertex : (*create_iknn_graph_res)[i]) {
            adj_vect[i].push_back(nn_vertex.index);
            intersection_size_vect[i].push_back(nn_vertex.common_count);
            dist_vect[i].push_back(nn_vertex.min_dist);
        }
    }

    // Preparing adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect[i].size()));
        int* A = INTEGER(RA);
        for (auto neighbor : adj_vect[i])
            *A++ = neighbor + 1;  // Adding 1 to convert to 1-based indexing for R
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Preparing intersection size list
    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect[i].size()));
        int* W = INTEGER(RW);
        for (auto size : intersection_size_vect[i])
            *W++ = size;
        SET_VECTOR_ELT(intersection_size_list, i, RW);
        UNPROTECT(1);
    }

    // Preparing distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;
    for (int i = 0; i < n_vertices; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, dist_vect[i].size()));
        double* D = REAL(RD);
        for (auto dist : dist_vect[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
        UNPROTECT(1);
    }

    // Modify this part
    auto adj_vect_ptr = std::make_unique<std::vector<std::vector<int>>>(std::move(adj_vect));

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect_ptr);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    // Preparing the result list
    SEXP res = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, dist_list);
    SET_VECTOR_ELT(res, 3, conn_comps);

    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("dist_list"));
    SET_STRING_ELT(names, 3, mkChar("conn_comps"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
