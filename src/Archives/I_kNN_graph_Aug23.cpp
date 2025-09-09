/*!
  Creating Intersetion Weighted kNN graph (IWkNN graph)

  Two points are connected by an edge in this graph if and only if they share at
  least one common point within their respective k-nearest neighbor sets. In
  other words, there's an intersection between their k-nearest neighborhoods.

  \textbf{Characteristics:}

  * \textbf{Reciprocity:} This connection is inherently reciprocal. If point A is within the k-nearest neighbors of point B, then point B will also be within the k-nearest neighbors of point A. This generally results in an undirected graph.
  * \textbf{Sensitivity to k:} The choice of the parameter k significantly influences the structure of this graph. Larger values of k tend to create denser graphs with more connections.
  * \textbf{Focus on Shared Neighbors:} This variant emphasizes points that are "locally similar," meaning they share commonalities within their local neighborhoods.

  \textbf{Potential Use Cases:}

  * \textbf{Clustering and Community Detection:} This type of graph could be valuable in identifying densely connected clusters of points that share neighbors, indicating groups exhibiting similar characteristics.
  * \textbf{Outlier Detection:}  Points with abnormally few connections in this graph might indicate potential outliers, as they don't consistently appear in the k-nearest neighborhoods of other points.
  * \textbf{Collaborative Filtering:} This construct might be applicable in recommendation systems where you want to emphasize items that have been "liked" by users with overlapping tastes.

  \textbf{Caveats:}

  * \textbf{Computational Cost:} Constructing this graph can be computationally more demanding than a standard k-NN graph because you need to compare k-nearest neighbor sets.
  * \textbf{Data Density:}  In very sparse datasets, this construct might lead to highly disconnected graphs.

  \textbf{Relationship with Other Graphs:}

  It's interesting to note that this graph would likely be denser than a Mutual k-Nearest Neighbor Graph (Mk-NN).  In an Mk-NN, two points need to be *each other's* k-nearest neighbors, while in this variant they just need to share at least one common neighbor in their sets.

*/

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
#include "I_kNN_graph.h"
#include "cpp_utils.h"

extern "C" {
    SEXP S_I_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_IW_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion, SEXP Rprune_version);
    SEXP S_IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho);
    SEXP S_IWD_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_IWD_kNN_graph_mp(SEXP RX, SEXP Rk);
}

std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<int>>>& adj_vect);
std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& adj_vect);
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> prune_long_edges(const std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& wgraph, int version);
std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph);


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

// ------------------------------------------------------------------------------------------
//
// IWD_kNN_graph
//
// ------------------------------------------------------------------------------------------

/**
 * Transforms an Intersection-Weighted-Distance k-Nearest Neighbor graph to an Intersection-Weighted k-Nearest Neighbor graph.
 *
 * This function converts a graph representation where each edge contains a distance measurement alongside the intersection size
 * to a representation that only considers the intersection size, discarding the distance information.
 *
 * @param graph The original graph represented as a vector of vectors of IkNN_vertex_t, where each IkNN_vertex_t contains
 *              the nearest neighbor index, intersection size, and distance.
 * @return A unique pointer to a vector of vectors of pairs, where each pair consists of a nearest neighbor index and an intersection size.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IWD_to_IW_kNN_graph(const std::vector<std::vector<IkNN_vertex_t>>& graph) {
    auto iw_graph = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(graph.size());

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto& vertex : graph[i]) {
            (*iw_graph)[i].emplace_back(vertex.nn_index, vertex.I_index);
        }
    }

    return iw_graph;
}


/**
 * @brief Constructs an Intersection-Weighted-Distance k-Nearest Neighbors (kNN) graph.
 *
 * This function builds an Intersection-Weighted-Distance kNN graph from a given data matrix.
 * For each point, it finds its k nearest neighbors, then constructs edges between points
 * based on the number of common nearest neighbors and computes the distance as the
 * minimum sum of distances through common neighbors.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @return A unique pointer to a 2D vector of IkNN_vertex_t structures representing the kNN graph.
 *
 * Each IkNN_vertex_t structure contains:
 * - `nn_index`: Index of the nearest neighbor.
 * - `I_index`: Number of common elements in the kNN sets.
 * - `dist`: Distance metric computed as the minimum sum of distances through common neighbors.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 */
std::unique_ptr<std::vector<std::vector<IkNN_vertex_t>>> IWD_kNN_graph(SEXP RX, SEXP Rk) {
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
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);

    auto res = std::make_unique<std::vector<std::vector<IkNN_vertex_t>>>(nrX);

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    std::vector<int> intersection;
    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        // Copying indices of kNN of the pt_i point to nn_i
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + nrX * j];
            sorted_nn_i[j] = nn_i[j];
        }

        std::sort(sorted_nn_i.begin(), sorted_nn_i.end()); // Ensure sorted for set intersection

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            // Copying indices of kNN of the pt_j point to nn_j
            for (int j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + nrX * j];
                sorted_nn_j[j] = nn_j[j];
            }

            std::sort(sorted_nn_j.begin(), sorted_nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(), sorted_nn_j.begin(), sorted_nn_j.end(), std::back_inserter(intersection));

            int common_count = intersection.size();
            if (common_count > 0) {
                // Computing the minimum of d(x,x_k) + d(x_k,x_j)
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    auto idx_i = std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin(); // The std::find function returns an iterator, not an index. We need to subtract the beginning iterator to get the index.
                    auto idx_j = std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = distances[pt_i + nrX * idx_i];
                    double dist_j_k = distances[pt_j + nrX * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                // Add edge from pt_i to pt_j and from pt_j to pt_i
                (*res)[pt_i].emplace_back(IkNN_vertex_t{pt_j, common_count, min_dist});
                (*res)[pt_j].emplace_back(IkNN_vertex_t{pt_i, common_count, min_dist});
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
 * S_IWD_kNN_graph
 *
 * This function computes a weighted k-nearest neighbors graph based on the input data matrix.
 * It returns a list containing the adjacency list of the graph, intersection sizes, distances,
 * connected components, pruned adjacency list, and sizes of cycles within the pruned graph.
 *
 * @param RX A numeric matrix (in column-major order) representing the data points.
 * @param Rk An integer vector of length 1 representing the number of nearest neighbors (k).
 *
 * @return A list with the following elements:
 *   - adj_list: A list of integer vectors where each vector represents the indices of the k-nearest neighbors for each point.
 *   - isize_list: A list of integer vectors where each vector represents the intersection sizes for the k-nearest neighbors.
 *   - dist_list: A list of numeric vectors where each vector represents the distances to the k-nearest neighbors.
 *   - conn_comps: An integer vector representing the connected components of the graph.
 *   - pruned_adj_list: A list of integer vectors where each vector represents the indices of the k-nearest neighbors for each point in the pruned graph.
 *
 * @note The function uses the base R C++ API and does not rely on Rcpp.
 *
 * @error Throws an error if the input data matrix cannot be coerced to a numeric matrix.
 * @error Throws an error if the k-nearest neighbors graph function returns a null pointer.
 *
 * Example:
 * \code
 * SEXP result = S_IWD_kNN_graph(RX, Rk);
 * \endcode
 */
#if 0
SEXP S_IWD_kNN_graph(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Creating a kNN graph
    std::unique_ptr<std::vector<std::vector<IkNN_vertex_t>>> IWD_kNN_graph_res = IWD_kNN_graph(RX, Rk);
    if (!IWD_kNN_graph_res) Rf_error("IWD_kNN_graph function returned an invalid result (null pointer).");

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto dist_vect = std::make_unique<std::vector<std::vector<double>>>(nrX);
    int IWD_kNN_graph_n_vertices = IWD_kNN_graph_res->size();
    for (int i = 0; i < IWD_kNN_graph_n_vertices; i++) {
        for (auto nn_vertex : IWD_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_vertex.nn_index);
            (*intersection_size_vect)[i].push_back(nn_vertex.I_index);
            (*dist_vect)[i].push_back(nn_vertex.dist);
        }
    }

    // Preparing adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Preparing intersection size list
    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);
        for (auto size : intersection_size_vect->at(i))
            *W++ = size;
        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    // Preparing distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, dist_vect->at(i).size())); nprot++;
        double* D = REAL(RD);
        for (auto dist : dist_vect->at(i))
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the IWD graph by deriving a IW kNN graph and then pruning that graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph = IWD_to_IW_kNN_graph(*IWD_kNN_graph_res);

    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph, 1);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_vertex : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_vertex.first);
        }
    }

    // Preparing pruned adjacency list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    // Preparing the result list
    SEXP res = PROTECT(allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, dist_list);
    SET_VECTOR_ELT(res, 3, conn_comps);
    SET_VECTOR_ELT(res, 4, pruned_adj_list);

    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("dist_list"));
    SET_STRING_ELT(names, 3, mkChar("conn_comps"));
    SET_STRING_ELT(names, 4, mkChar("pruned_adj_list"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
#endif



struct IkNN_vertex_tt {
    int index;
    int common_count;
    double min_dist;
};

// Helper function to perform kNN search
SEXP kNN_search(ANNpointArray data_pts, int nrX, int ncX, int k) {
    ANNkd_tree* kdTree = new ANNkd_tree(data_pts, nrX, ncX);
    std::vector<int> nnIdx(k);
    std::vector<double> dists(k);

    SEXP indices = PROTECT(allocMatrix(INTSXP, nrX, k));
    SEXP distances = PROTECT(allocMatrix(REALSXP, nrX, k));
    int *indicesPtr = INTEGER(indices);
    double *distancesPtr = REAL(distances);

    for (int i = 0; i < nrX; i++) {
        kdTree->annkSearch(data_pts[i], k, nnIdx.data(), dists.data());
        for (int j = 0; j < k; j++) {
            indicesPtr[i + nrX * j] = nnIdx[j];
            distancesPtr[i + nrX * j] = dists[j];
        }
    }

    delete kdTree;

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, indices);
    SET_VECTOR_ELT(result, 1, distances);

    UNPROTECT(3);
    return result;
}

SEXP IWD_kNN_graph_mp(SEXP RX, SEXP Rk) {
    SEXP dims = PROTECT(getAttrib(RX, R_DimSymbol));
    int nrX = INTEGER(dims)[0];
    int ncX = INTEGER(dims)[1];
    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    double *X = REAL(RX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j];
        }
    }

    // Perform kNN search
    SEXP knn_res = PROTECT(kNN_search(data_pts, nrX, ncX, k));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    // Prepare result container
    std::vector<std::vector<IkNN_vertex_tt>> res(nrX);

    // Temporary vectors for set intersection
    std::vector<int> nn_i(k), nn_j(k), intersection;

    // Main loop for graph construction
    omp_set_num_threads(10);
    #pragma omp parallel for schedule(dynamic) private(nn_i, nn_j, intersection)
    for (int pt_i = 0; pt_i < nrX - 1; pt_i++) {
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + nrX * j];
        }
        std::sort(nn_i.begin(), nn_i.end());

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            for (int j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + nrX * j];
            }
            std::sort(nn_j.begin(), nn_j.end());

            intersection.clear();
            std::set_intersection(nn_i.begin(), nn_i.end(), nn_j.begin(), nn_j.end(),
                                  std::back_inserter(intersection));

            int common_count = intersection.size();
            if (common_count > 0) {
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    auto idx_i = std::lower_bound(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin();
                    auto idx_j = std::lower_bound(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = distances[pt_i + nrX * idx_i];
                    double dist_j_k = distances[pt_j + nrX * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                #pragma omp critical
                {
                    res[pt_i].emplace_back(IkNN_vertex_tt{pt_j, common_count, min_dist});
                    res[pt_j].emplace_back(IkNN_vertex_tt{pt_i, common_count, min_dist});
                }
            }
        }
    }

    // Convert result to R list
    SEXP R_res = PROTECT(allocVector(VECSXP, nrX));
    for (int i = 0; i < nrX; i++) {
        int n = res[i].size();
        SEXP indices = PROTECT(allocVector(INTSXP, n));
        SEXP common_counts = PROTECT(allocVector(INTSXP, n));
        SEXP min_dists = PROTECT(allocVector(REALSXP, n));

        for (int j = 0; j < n; j++) {
            INTEGER(indices)[j] = res[i][j].index;
            INTEGER(common_counts)[j] = res[i][j].common_count;
            REAL(min_dists)[j] = res[i][j].min_dist;
        }

        SEXP sublist = PROTECT(allocVector(VECSXP, 3));
        SET_VECTOR_ELT(sublist, 0, indices);
        SET_VECTOR_ELT(sublist, 1, common_counts);
        SET_VECTOR_ELT(sublist, 2, min_dists);

        SET_VECTOR_ELT(R_res, i, sublist);
        UNPROTECT(4);
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose();

    UNPROTECT(3);
    return R_res;
}


SEXP S_IWD_kNN_graph_mp(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Call IWD_kNN_graph_mp (assuming it returns SEXP)
    SEXP IWD_kNN_graph_res = PROTECT(IWD_kNN_graph_mp(RX, Rk)); nprot++;

    // Prepare vectors to hold the data
    std::vector<std::vector<int>> adj_vect(nrX);
    std::vector<std::vector<int>> intersection_size_vect(nrX);
    std::vector<std::vector<double>> dist_vect(nrX);

    // Extract data from IWD_kNN_graph_res
    for (int i = 0; i < nrX; i++) {
        SEXP sublist = VECTOR_ELT(IWD_kNN_graph_res, i);
        SEXP indices = VECTOR_ELT(sublist, 0);
        SEXP common_counts = VECTOR_ELT(sublist, 1);
        SEXP min_dists = VECTOR_ELT(sublist, 2);

        int n = Rf_length(indices);
        for (int j = 0; j < n; j++) {
            adj_vect[i].push_back(INTEGER(indices)[j]);
            intersection_size_vect[i].push_back(INTEGER(common_counts)[j]);
            dist_vect[i].push_back(REAL(min_dists)[j]);
        }
    }

    // Preparing adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect[i].size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : adj_vect[i])
            *A++ = neighbor + 1;  // Adding 1 to convert to 1-based indexing for R
        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Preparing intersection size list
    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect[i].size())); nprot++;
        int* W = INTEGER(RW);
        for (auto size : intersection_size_vect[i])
            *W++ = size;
        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    // Preparing distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, dist_vect[i].size())); nprot++;
        double* D = REAL(RD);
        for (auto dist : dist_vect[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
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
