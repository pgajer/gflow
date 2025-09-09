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
#include "kNN.h"         // for struct IkNN_vertex_tt and kNN_search
#include "IWD_kNN_graphs.h"
#include "cpp_utils.h"
#include "SEXP_cpp_conversion_utils.h"

extern "C" {
    SEXP S_I_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_IW_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion, SEXP Rprune_version);
    SEXP S_IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho);
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

// ------------------------------------------------------------------------------------------
//
// I_kNN_graph
//
// ------------------------------------------------------------------------------------------

/*!
  Creates an intersection kNN graph without weights

  This function creates an undirected unweighted graph that connects two points
  of X with an edge if the sets of kNN of these points intersect. Thus, the
  graph is the 1-skeleton of the Cech simplicial complex of the kNN covering of
  X.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a pointer to a vector of vectors, adj_kNN, such that
  adj_kNN[i] is the set of neighbors of the i-th vertex.
*/
std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *nn_indices = INTEGER(VECTOR_ELT(knn_res, 0));

    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX);

    int n_points_minus_one = nrX - 1;
    std::set<int> knn_set_i;
    std::set<int> knn_set_j;
    std::vector<int> intersection; // For storing the intersection result
    for (int point_i = 0; point_i < n_points_minus_one; point_i++) {
        knn_set_i.clear();
        for (int j = 0; j < k; j++)
            knn_set_i.insert(nn_indices[point_i + nrX * j]);

        for (int point_j = point_i + 1; point_j < nrX; point_j++) {
            knn_set_j.clear();
            for (int j = 0; j < k; j++)  // j = 0 is the point whose NN we are computing, so we skip it
                knn_set_j.insert(nn_indices[point_j + nrX * j]);

            intersection.clear();
            std::set_intersection(knn_set_i.begin(), knn_set_i.end(),
                                  knn_set_j.begin(), knn_set_j.end(),
                                  std::back_inserter(intersection));

            if (!intersection.empty()) {
                (*res)[point_i].push_back(point_j);
                (*res)[point_j].push_back(point_i);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates an Intersection kNN graph (with no weights)

  This function calls an internal C++ function I_kNN_graph() to create an interesection kNN graph.

  \param RX          A matrix of points.
  \param Rk          The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns an adjacency list of point indices adjacent to each vertex.
*/
SEXP S_I_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) { // This is here only in case the R function that calls this one is changed in a way that removes the conversion of X to a matrix
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Creating an intersection kNN graph
    std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph_res = I_kNN_graph(RX, Rk);

    // Preparing adj_list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, I_kNN_graph_res->at(i).size())); nprot++; // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : I_kNN_graph_res->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> I_kNN_graph_conn_comp = union_find(I_kNN_graph_res);

    SEXP conn_comps = PROTECT(allocVector(INTSXP, I_kNN_graph_conn_comp->size())); nprot++; // Graph connected components
    std::copy(I_kNN_graph_conn_comp->begin(), I_kNN_graph_conn_comp->end(), INTEGER(conn_comps));

    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, conn_comps);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("conn_comps"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}


// ------------------------------------------------------------------------------------------
//
// IW_kNN_graph
//
// ------------------------------------------------------------------------------------------

/*!
  Creates an intersection weighted kNN graph

  This function creates a weighted graph that connects two points of X with an
  edge if the sets of kNN of these points intersect. Thus, the graph is the
  1-skeleton of the nerve simplicial complex of the kNN covering of X. The weight
  of an edge is the number of points in common between two sets of kNN's.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a vector of vector of pairs, adj_kNN, such that adj_kNN[i] is
  a vector of pairs: (i1, w1), (i2, w2), ... , (iN, wN), with each pair being
  the index of the neighbor and the corresponding edge weight.
*/
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    int common_count = 0;
    std::vector<int> intersection;
    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        for (int j = 0; j < k; j++)  // Copying indices of kNN of the pt_i point to nn_i
            nn_i[j] = indices[pt_i + nrX * j];

        std::sort(nn_i.begin(), nn_i.end()); // Ensure sorted for set intersection

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {

            for (int j = 0; j < k; j++) // Copying indices of kNN of the pt_j point to nn_j
                nn_j[j] = indices[pt_j +  nrX * j];

            std::sort(nn_j.begin(), nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(nn_i.begin(), nn_i.end(), nn_j.begin(), nn_j.end(), std::back_inserter(intersection));

            common_count = intersection.size();
            if (common_count > 0) {
                // Add edge from pt_i to pt_j and from pt_j to pt_i with the weight
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

//
// Version 2 utilizing std::unordered_set<int> set_i/set_j
//
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_v2(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    int common_count = 0;

    std::unordered_set<int> set_i;
    std::unordered_set<int> set_j;

    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        set_i.clear();
        // Populate set_i for pt_i
        for (int j = 0; j < k; j++) {
            set_i.insert(indices[pt_i + nrX * j]);
        }

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            set_j.clear();
            for (int j = 0; j < k; j++) {
                set_j.insert(indices[pt_j + nrX * j]);
            }

#if 0
            Rprintf("i: %d\n", (int)pt_i);
            Rprintf("j: %d\n", (int)pt_j);
            print_uset(set_i, "set_i");
            print_uset(set_j, "set_j");
            error("Stopping for debugging purposes");
#endif

            common_count = 0;
            // If set_i is smaller, iterate over it for efficiency
            if (set_i.size() < set_j.size()) {
                for (const int& elem : set_i) {
                    if (set_j.count(elem) > 0) {
                        common_count++;
                    }
                }
            } else {
                for (const int& elem : set_j) {
                    if (set_i.count(elem) > 0) {
                        common_count++;
                    }
                }
            }

            if (common_count > 0) {
                // Add edge from pt_i to pt_j and from pt_j to pt_i with the weight
                // Assuming res is a properly initialized structure to hold these values
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates an intersection weighted kNN graph

  This function calls an internal C++ function kNN_graph() to create a kNN graph.

  \param RX               A matrix of points.
  \param Rk               The number of nearest neighbors to use for the construction of the kNN graph.
  \param Rversion         A version of the IW_kNN_graph() function.
  \param Rprune_version   A version of the prune_long_edges() function to use.

  \return A list with four components:
  1. adj_list: An adjacency list of point indices adjacent to each vertex.
  2. isize_list: A list of the sizes of the intersections of kNN sets.
  3. conn_comps: A vector of graph connected component indicators.
  4. pruned_adj_list: The pruned graph with 'long' edges removed.

*/
SEXP S_IW_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion, SEXP Rprune_version) {

    int version = INTEGER(Rversion)[0];
    if (version < 1 || version > 2) {
        Rf_error("Invalid 'Rversion' value. Must be 1 or 2");
    }

    int prune_version = INTEGER(Rprune_version)[0];
    if (prune_version < 1 || prune_version > 2) {
        Rf_error("Invalid 'Rprune_version' value. Must be 1 or 2");
    }

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    //
    // Creating a kNN graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_res;
    if (version == 1) {
        IW_kNN_graph_res = IW_kNN_graph(RX, Rk);
        if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph function returned an invalid result (null pointer).");
    } else if (version == 2) {
        IW_kNN_graph_res = IW_kNN_graph_v2(RX, Rk);
        if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph_v2 function returned an invalid result (null pointer).");
    }

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    int IW_kNN_graph_n_vertices = IW_kNN_graph_res->size();
    for (int i = 0; i < IW_kNN_graph_n_vertices; i++) {
        for (auto nn_pair : IW_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_pair.first);
            (*intersection_size_vect)[i].push_back(nn_pair.second);
        }
    }

    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);

        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);

        for (auto dist : intersection_size_vect->at(i))
            *W++ = dist;

        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    //
    // Identifying connected components of the graph
    //
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph_res, prune_version);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_pair : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_pair.first);
        }
    }

    // Preparing pruned_I_kNN_graph adj list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    //
    // Finding the sizes of cycles within the pruned graph
    //
    #if 0
    std::unique_ptr<std::vector<int>> pruned_graph_cycle_lengths_ptr = cycle_sizes(pruned_graph_adj_vect);
    SEXP pruned_graph_cycle_lengths_Rvector = PROTECT(allocVector(INTSXP, pruned_graph_cycle_lengths_ptr->size())); nprot++;
    std::copy(pruned_graph_cycle_lengths_ptr->begin(), pruned_graph_cycle_lengths_ptr->end(), INTEGER(pruned_graph_cycle_lengths_Rvector));
    #endif

    SEXP res = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, conn_comps);
    SET_VECTOR_ELT(res, 3, pruned_adj_list);
    //SET_VECTOR_ELT(res, 4, pruned_graph_cycle_lengths_Rvector);

    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("conn_comps"));
    SET_STRING_ELT(names, 3, mkChar("pruned_adj_list"));
    //SET_STRING_ELT(names, 4, mkChar("pruned_cycles_sizes"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

// ------------------------------------------------------------------------------------------
//
// IW_rho_kNN_graph functions
//
// ------------------------------------------------------------------------------------------

/**
 * @brief Constructs a rho-corrected intersection-weighted k-Nearest Neighbors (kNN) graph.
 *
 * This function builds a rho-corrected intersection-weighted kNN graph from a given data matrix.
 * For each point, it finds its k nearest neighbors with a correction factor rho, then constructs
 * edges between points based on the number of common nearest neighbors.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @param Rrho A numeric vector (SEXP) representing the correction factors for each point.
 * @return A unique pointer to a 2D vector of pairs representing the kNN graph. Each pair contains:
 *         - The index of the nearest neighbor.
 *         - The number of common elements in the kNN sets.
 *
 * Each point in the graph is connected to its nearest neighbors, with weights determined by the
 * intersection of their kNN sets. The distance between points is adjusted by a correction factor
 * rho, which normalizes the distances.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 *
 * @note The kNN search is performed with 2*k neighbors initially to allow for rho correction.
 * This function efficiently handles the case where the number of points is less than 2*k.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;
    PROTECT(Rrho = coerceVector(Rrho, REALSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];
    int k = INTEGER(Rk)[0];
    double *rho = REAL(Rrho);

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X with K = 2 * k
    SEXP Rtwo_k = PROTECT(allocVector(INTSXP, 1)); nprot++;
    int two_k = 2 * k;
    if (two_k >= nrX) {
        two_k = nrX - 1;
    }
    INTEGER(Rtwo_k)[0] = two_k;

    // Finding kNN for K = 2*k
    SEXP knn_res = PROTECT(S_kNN(RX, Rtwo_k)); nprot++;
    int *two_k_indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *two_k_distances = REAL(VECTOR_ELT(knn_res, 1));

    // Indices of kNN's computed with d1 distance corrections
    std::vector<int> k_indices(k * nrX);
    std::vector<std::pair<int,double>> nn_idist(two_k); // ANN includes the query point as the first of the kNN's, we are sorting kNN's excluding the query point
    int nn_j;

    for (int i = 0; i < nrX; i++) {
        nn_idist.clear();  // Reset for the next point
        for (int j = 1; j < two_k; j++) {
            nn_j = two_k_indices[i + nrX * j];
            nn_idist.emplace_back(nn_j, two_k_distances[i + nrX * j] / rho[nn_j]);
        }
        // Sorting nn_idist in the ascending order with respect to the second component of each pair
        std::sort(nn_idist.begin(), nn_idist.end(),
                  [&](std::pair<int,double> a, std::pair<int,double> b) {
            return a.second < b.second;
        });

        // Populating k_indices with d1-corrected kNN's indices
        k_indices[i] = i;
        for (int j = 1; j < k; j++)
            k_indices[i + nrX * j] = nn_idist[j - 1].first;
    }

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    int common_count = 0;

    std::unordered_set<int> set_i;
    std::unordered_set<int> set_j;

    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        set_i.clear();
        // Populate set_i for pt_i
        for (int j = 0; j < k; j++)
            set_i.insert(k_indices[pt_i + nrX * j]);

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            set_j.clear();
            for (int j = 0; j < k; j++)
                set_j.insert(k_indices[pt_j + nrX * j]);

            common_count = 0;
            // If set_i is smaller, iterate over it for efficiency
            if (set_i.size() < set_j.size()) {
                for (const int& elem : set_i) {
                    if (set_j.count(elem) > 0) {
                        common_count++;
                    }
                }
            } else {
                for (const int& elem : set_j) {
                    if (set_i.count(elem) > 0) {
                        common_count++;
                    }
                }
            }

            if (common_count > 0) {
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/**
 * @brief Constructs a rho-corrected intersection-weighted k-Nearest Neighbors (kNN) graph and processes it.
 *
 * This function builds a rho-corrected intersection-weighted kNN graph from a given data matrix and processes
 * the graph to produce several outputs: adjacency list, intersection sizes, connected components, pruned
 * adjacency list, and cycle sizes in the pruned graph.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 *           It is coerced to a numeric matrix if not already in the correct format.
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @param Rrho A numeric vector (SEXP) representing the correction factors for each point.
 * @return A list (SEXP) with the following components:
 *         - `adj_list`: A list of integer vectors, each representing the adjacency list for a point.
 *         - `isize_list`: A list of integer vectors, each representing the intersection size with neighbors.
 *         - `conn_comps`: An integer vector representing the connected components of the graph.
 *         - `pruned_adj_list`: A list of integer vectors, each representing the adjacency list of the pruned graph.
 *         - `pruned_cycles_sizes`: An integer vector representing the sizes of cycles within the pruned graph.
 *
 * The function performs the following steps:
 * 1. Constructs a rho-corrected intersection-weighted kNN graph using `IW_rho_kNN_graph`.
 * 2. Extracts the adjacency lists and intersection sizes from the graph.
 * 3. Identifies connected components of the graph using a union-find algorithm.
 * 4. Prunes long edges from the graph using `prune_long_edges`.
 * 5. Finds the sizes of cycles within the pruned graph using `cycle_sizes`.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 *
 * @note The kNN search is performed with 2*k neighbors initially to allow for rho correction.
 * This function efficiently handles the case where the number of points is less than 2*k.
 *
 * @throws Rf_error if RX cannot be coerced to a numeric matrix or if `IW_rho_kNN_graph` returns a null pointer.
 */
SEXP S_IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    //
    // Creating a kNN graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_res = IW_rho_kNN_graph(RX, Rk, Rrho);
    if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph function returned an invalid result (null pointer).");

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    int IW_kNN_graph_n_vertices = IW_kNN_graph_res->size();
    for (int i = 0; i < IW_kNN_graph_n_vertices; i++) {
        for (auto nn_pair : IW_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_pair.first);
            (*intersection_size_vect)[i].push_back(nn_pair.second);
        }
    }

    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);

        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);

        for (auto dist : intersection_size_vect->at(i))
            *W++ = dist;

        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    //
    // Identifying connected components of the graph
    //
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph_res, 1);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_pair : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_pair.first);
        }
    }

    // Preparing pruned_I_kNN_graph adj list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    //
    // Finding the sizes of cycles within the pruned graph
    //
    std::unique_ptr<std::vector<int>> pruned_graph_cycle_lengths_ptr = cycle_sizes(pruned_graph_adj_vect);
    SEXP pruned_graph_cycle_lengths_Rvector = PROTECT(allocVector(INTSXP, pruned_graph_cycle_lengths_ptr->size())); nprot++;
    std::copy(pruned_graph_cycle_lengths_ptr->begin(), pruned_graph_cycle_lengths_ptr->end(), INTEGER(pruned_graph_cycle_lengths_Rvector));


    SEXP res = PROTECT(allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, conn_comps);
    SET_VECTOR_ELT(res, 3, pruned_adj_list);
    SET_VECTOR_ELT(res, 4, pruned_graph_cycle_lengths_Rvector);

    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("conn_comps"));
    SET_STRING_ELT(names, 3, mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 4, mkChar("pruned_cycles_sizes"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

