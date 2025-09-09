/*!
  Creating kNN graph
*/

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <set>
#include <memory>
#include <vector>

#include "msr2.h"

/*!
  Creates a kNN graph

  There are different constructs of graphs associated with the set of kNN of the
  given set of points. This function connects two vertices with an edge if the
  corresponding sets of kNN of these points intersect. Thus, the graph is the
  1-skeleton of the Cech simplicial complex of the kNN covering of X.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a vector of sets, adj_kNN, such that adj_kNN[j] is the
  set of indices, i1, i2, ... , iN, of points of X such that j belongs to
  N(i_j,k), for all these indices. Thus, adj_kNN[j] is the set of vertices
  connected with j by an edge.
*/
std::unique_ptr<std::vector<std::set<int>>> kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int k = INTEGER(Rk)[0];
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    // Convert RX matrix to ANNpointArray
    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    ANNkd_tree *kdTree = new ANNkd_tree(dataPts, nrX, ncX);

    // Allocate memory for nearest neighbor indices and distances
    ANNidxArray nnIdx = new ANNidx[k];
    ANNdistArray nnDist = new ANNdist[k];

    auto adj_kNN = std::make_unique<std::vector<std::set<int>>>(nrX); // Initialize adj_kNN with nrX empty sets

    // Perform k-NN search for each point
    for (int i = 0; i < nrX; i++) {
        kdTree->annkSearch(dataPts[i], k, nnIdx, nnDist, 0); // 0 for exact search

        for (int j = 0; j < k; j++) {
            int nnIdx_j = nnIdx[j];
            if (nnIdx_j != i) { // Exclude self-loops
                (*adj_kNN)[i].insert(nnIdx_j);
                (*adj_kNN)[nnIdx_j].insert(i); // Add the symmetric edge
            }
        }
    }

    // Clean up
    delete[] nnIdx;
    delete[] nnDist;
    delete kdTree;
    annDeallocPts(dataPts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return adj_kNN;
}

/*!
Creates a kNN graph

This function calls an internal C++ function kNN_graph() to create a kNN graph.

\param X A matrix of points.
\param k The number of nearest neighbors to use for the construction of the kNN graph.

\return Returns list with two components: 1) the vector of vertex indices 2)
an m-by-2 matrix of edge indices.
*/
SEXP S_kNN_graph(SEXP RX, SEXP Rk) {

    std::unique_ptr<std::vector<std::set<int>>> res = kNN_graph(RX, Rk);

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Creates a vector seq(nrX) of the graph indices
    SEXP vertex_indices = PROTECT(allocVector(INTSXP, nrX)); nprot++;
    int *vertex_indices_ptr = INTEGER(vertex_indices);
    for (int i = 0; i < nrX; i++) {
        vertex_indices_ptr[i] = i + 1;
    }

    // Count the total number of edges in the graph
    int num_edges = 0;
    for (int i = 0; i < nrX; i++) {
        num_edges += res->at(i).size();
    }

    // Creates from res an m-by-2 matrix with indices of edges of the graph
    SEXP edge_indices = PROTECT(allocMatrix(INTSXP, num_edges, 2)); nprot++;
    int *edge_indices_ptr = INTEGER(edge_indices);
    int edge_index = 0;
    for (int i = 0; i < nrX; i++) {
        for (int neighbor : res->at(i)) {
            edge_indices_ptr[edge_index] = i + 1;
            edge_indices_ptr[edge_index + num_edges] = neighbor + 1;
            edge_index++;
        }
    }

    // Creates a list with two components: 1) the vector of vertex indices 2) an m-by-2 matrix of edge indices.
    SEXP result = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(result, 0, vertex_indices);
    SET_VECTOR_ELT(result, 1, edge_indices);

    // Add names to list elements
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("vertices"));
    SET_STRING_ELT(names, 1, mkChar("edges"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result;
}
