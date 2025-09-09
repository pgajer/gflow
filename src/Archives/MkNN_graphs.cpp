/*!
  Creating Mutual kNN graph (MkNN graph)

  Two points are connected if and only if they are within each other's k-nearest neighbor lists.

*/

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <unordered_set>
#include <utility>
#include <memory>

#include "msr2.h"

extern "C" {

    SEXP S_M_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion);
    SEXP S_MW_kNN_graph(SEXP RX, SEXP Rk);
}

/*!
  Creates a mutual kNN graph

  This function creates a mutual kNN graph that connects two points of X with an
  edge if and only if they are within each other's k-nearest neighbor lists.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a vector, A, of int vectors, such that A[j] is the
  vector of indices of X that are connected with X[j] by an edge.
*/
std::unique_ptr<std::vector<std::vector<int>>> M_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_M_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    // Error Checking:
    if (indices == NULL) {
        Rf_error("M_kNN_graph: Error in S_kNN: kNN index extraction failed.");
    }

    std::vector<std::unordered_set<int>> nn_sets(nrX);

    // Populate nn_sets with the k-nearest neighbors of each point
    for (int pt_i = 0; pt_i < nrX; ++pt_i) {
        for (int j = 0; j < k; ++j) {
            nn_sets[pt_i].insert(indices[pt_i + nrX * j]);
        }
    }

    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX); // Initialize a pointer to an empty adjacency list

#if DEBUG_M_kNN_graph
    Rprintf("In M_kNN_graph before the main loop\n");
#endif

    // Find mutual k-nearest neighbors and add edges
    for (int pt_i = 0; pt_i < nrX; ++pt_i) {

#if DEBUG_M_kNN_graph
        Rprintf("pt_i: %d\n", pt_i);
        Rprintf("nn_sets[pt_i]: ");
        for (int pt_j : nn_sets[pt_i])
            Rprintf(" %d,", pt_j);
        Rprintf("\n");
#endif

        for (int pt_j : nn_sets[pt_i]) {
            if (pt_i < pt_j && nn_sets[pt_j].count(pt_i) > 0) {
                (*res)[pt_i].push_back(pt_j);
                (*res)[pt_j].push_back(pt_i);
            }
        }
    }

    // Clean up
    annDeallocPts(dataPts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

//
// Version 2
//
std::unique_ptr<std::vector<std::vector<int>>> M_kNN_graph_v2(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    // Error Checking:
    if (indices == NULL) {
        Rf_error("M_kNN_graph_v2: Error in S_kNN: kNN index extraction failed.");
    }

    std::vector<std::unordered_set<int>> nn_sets(nrX);

    // Populate nn_sets with the k-nearest neighbors of each point
    for (int i = 0; i < nrX; ++i) {
        for (int j = 0; j < k; ++j) {
            nn_sets[i].insert(indices[i + nrX * j]);
        }
    }

    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX); // Initialize a pointer to an empty adjacency list

    // Create the mutual kNN graph
    for (int i = 0; i < nrX; ++i) {
        for (int neighbor : nn_sets[i]) {
            // Check if the relationship is mutual
            if (i < neighbor && nn_sets[neighbor].find(i) != nn_sets[neighbor].end()) {
                (*res)[i].push_back(neighbor);
                (*res)[neighbor].push_back(i);
            }
        }
    }

    // Clean up
    annDeallocPts(dataPts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}


/*!
  Creates a Mutual kNN graph

  This function calls an internal C++ function M_kNN_graph() to create a kNN graph.

  \param RX          A matrix of points.
  \param Rk          The number of nearest neighbors to use for the construction of the kNN graph.
  \param Rversion    A version of the IW_kNN_graph() function.

  \return Returns an adjacency list of point indices adjacent to each vertex.
*/
SEXP S_M_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion) {

    int version = INTEGER(Rversion)[0];
    if (version < 1 || version > 2) {
        Rf_error("Invalid 'Rversion' value. Must be 1 or 2");
    }

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) { // This is here only in case the R function that calls this one is changed in a way that removes the conversion of X to a matrix
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    std::unique_ptr<std::vector<std::vector<int>>> M_kNN_graph_res;
    if (version == 1) {
        M_kNN_graph_res = M_kNN_graph(RX, Rk);
        if (!M_kNN_graph_res) Rf_error("M_kNN_graph function returned an invalid result (null pointer).");
    } else if (version == 2) {
        M_kNN_graph_res = M_kNN_graph_v2(RX, Rk);
        if (!M_kNN_graph_res) Rf_error("M_kNN_graph_v2 function returned an invalid result (null pointer).");
    }

#if 0
    // Count the total number of edges in the graph
    int num_edges = 0;
    for (int i = 0; i < nrX; i++) {
        num_edges += M_kNN_graph_res->at(i).size();
    }

    // Returns an m-by-2 matrix edge index pairs, where m is the number of edges.
    // Creates from M_kNN_graph_res an m-by-2 matrix with indices and weights of edges of the graph
    SEXP edge_matrix = PROTECT(allocMatrix(INTSXP, num_edges, 2)); nprot++;
    int *edge_matrix_ptr = INTEGER(edge_matrix);
    int edge_index = 0;
    for (int i = 0; i < nrX; i++) {
        for (int neighbor : M_kNN_graph_res->at(i)) {
            edge_matrix_ptr[edge_index] = i + 1;
            edge_matrix_ptr[edge_index + num_edges] = neighbor + 1;
            edge_index++;
        }
    }
#endif

    // Preparing return list
    SEXP result = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, M_kNN_graph_res->at(i).size())); // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : M_kNN_graph_res->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(result, i, RA);
        UNPROTECT(1);
    }

    UNPROTECT(nprot);

    return result;
}

/*!
  Creates a mutual kNN graph with weights

  This function creates a mutual kNN graph that connects two points of X with an
  edge if and only if they are within each other's k-nearest neighbor lists. The
  weight of the edge is the distance between the vertices in X.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a list with two components: 1) adjacency list, A, such that
  A[[i]] is a vector of vertices connected with the i-th vertex by an edge; 2)
  weights list, W, such that W[[i]] is a vector of distances between the i-th
  vertex and the neighbor vertices of that vertex.
*/
SEXP S_MW_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_MW_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    if (indices == NULL) {
        Rf_error("MW_kNN_graph: Error in S_kNN: kNN index extraction failed.");
    }

    if (distances == NULL) {
        Rf_error("MW_kNN_graph: Error in S_kNN: kNN distances extraction failed.");
    }

    std::vector<std::unordered_set<int>> nn_sets(nrX);
    // Populate nn_sets with the k-nearest neighbors of each point
    for (int pt_i = 0; pt_i < nrX; ++pt_i) {
        for (int j = 0; j < k; ++j) {
            nn_sets[pt_i].insert(indices[pt_i + nrX * j]);
        }
    }

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);       // Adjacency vector
    auto dists_vect = std::make_unique<std::vector<std::vector<double>>>(nrX);  // Distances vector

    // Find mutual k-nearest neighbors and add edges
    double d;
    for (int pt_i = 0; pt_i < nrX; ++pt_i) {
        for (int pt_j : nn_sets[pt_i]) {
            if (pt_i < pt_j && nn_sets[pt_j].count(pt_i) > 0) {
                (*adj_vect)[pt_i].push_back(pt_j);
                (*adj_vect)[pt_j].push_back(pt_i);
                d = distances[pt_i + nrX * pt_j];
                (*dists_vect)[pt_i].push_back(d);
                (*dists_vect)[pt_j].push_back(d);
            }
        }
    }

    // Clean up
    annDeallocPts(dataPts);
    annClose(); // Close ANN

    // Creating adj_list associated with adj_vect
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }

    // Creating dists_list associated with dists_vect
    SEXP dists_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(REALSXP, dists_vect->at(i).size())); // Dists list of the i-th vertex
        double* W = REAL(RW);

        for (auto dist : dists_vect->at(i))
            *W++ = dist;

        SET_VECTOR_ELT(dists_list, i, RW);
        UNPROTECT(1);
    }

    // Prepare return list
    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++; // List with 2 elements
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, dists_list);

    // Add names to list elements
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("weight_list"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
