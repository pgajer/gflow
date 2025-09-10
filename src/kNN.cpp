#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include <vector>
#include <cmath>
#include <ANN/ANN.h>  // ANN library header

#include "kNN.h"

extern "C" {
    SEXP S_kNN(SEXP RX, SEXP Rk);
    SEXP S_kNN_v2(SEXP RX, SEXP Rk);
}

/**
 * @brief Performs k-Nearest Neighbors search on a set of points.
 *
 * This function uses the ANN (Approximate Nearest Neighbor) library to perform
 * an efficient k-NN search on the input data points.
 *
 * @param X A vector of vectors representing the input data points.
 *          Each inner vector represents a single data point.
 * @param k The number of nearest neighbors to find for each point.
 *
 * @return A struct containing two arrays:
 *         - indices: An array of size nrX * k containing the indices of the k nearest neighbors
 *                    for each point. The neighbors for point i are stored in positions
 *                    [i*k] to [i*k + k-1].
 *         - distances: An array of size nrX * k containing the distances to the k nearest neighbors
 *                      for each point, in the same order as the indices.
 *
 * @note The caller is responsible for freeing the memory allocated for
 *       indices and distances in the returned struct.
 *
 * @Rf_warning This function modifies the global state of the ANN library.
 *          Make sure to call annClose() when you're done using ANN functions.
 */
knn_result_t kNN(const std::vector<std::vector<double>>& X, int k) {
    int nrX = X.size();
    int ncX = X[0].size();

    // Convert X to ANNpointArray
    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i][j];
        }
    }

    ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrX, ncX);

    // Allocate memory for nearest neighbor indices and distances
    knn_result_t result;
    result.indices = new int[nrX * k];
    result.distances = new double[nrX * k];

    ANNidxArray nnIdx = new ANNidx[k];
    ANNdistArray nnDist = new ANNdist[k];

    // Perform k-NN search for each point
    for (int i = 0; i < nrX; i++) {
        kdTree->annkSearch(dataPts[i], k, nnIdx, nnDist, 0); // 0 for Rf_error bound
        for (int j = 0; j < k; j++) {
            result.indices[i * k + j] = nnIdx[j];
            result.distances[i * k + j] = std::sqrt(nnDist[j]);
        }
    }

    // Clean up
    delete[] nnIdx;
    delete[] nnDist;
    delete kdTree;
    annDeallocPts(dataPts);
    annClose(); // Close ANN

    return result;
}

/**
 * @brief R interface for k-Nearest Neighbors search.
 *
 * This function provides an interface between R and the C++ kNN function.
 * It takes R objects as input, performs the k-NN search using the kNN function,
 * and returns the results in a format that can be easily used in R.
 *
 * @param RX An R matrix (REALSXP) containing the input data points.
 *           Each row represents a single data point.
 * @param Rk An R integer (INTSXP) specifying the number of nearest neighbors to find.
 *
 * @return An R list containing two elements:
 *         - indices: An integer matrix of size nrX x k containing the indices of the k nearest neighbors
 *                    for each point.
 *         - distances: A numeric matrix of size nrX x k containing the distances to the k nearest neighbors
 *                      for each point.
 *
 * @note This function uses R's protection mechanism to prevent garbage collection
 *       of temporary R objects. The number of protections is balanced with unprotections
 *       at the end of the function.
 *
 * @Rf_warning This function assumes that the input matrix RX has dimension attributes.
 *          Unexpected behavior may occur if RX is not a proper R matrix.
 */
SEXP S_kNN_v2(SEXP RX, SEXP Rk) {
    // Protect R objects from garbage collection
    int nprot = 0;
    PROTECT(RX = Rf_coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = Rf_coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int k = INTEGER(Rk)[0];
    int *dimX = INTEGER(Rf_getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    // Convert RX matrix to std::vector<std::vector<double>>
    std::vector<std::vector<double>> X_vec(nrX, std::vector<double>(ncX));
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            X_vec[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Call kNN function
    knn_result_t result = kNN(X_vec, k);

    // Create R objects to return
    SEXP nn_i = PROTECT(Rf_allocMatrix(INTSXP, nrX, k)); nprot++;
    SEXP nn_d = PROTECT(Rf_allocMatrix(REALSXP, nrX, k)); nprot++;

    // Copy results to R objects
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < k; j++) {
            INTEGER(nn_i)[i + nrX * j] = result.indices[i * k + j];
            REAL(nn_d)[i + nrX * j] = result.distances[i * k + j];
        }
    }

    // Clean up
    delete[] result.indices;
    delete[] result.distances;

    // Prepare return list
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++; // List with 2 elements
    SET_VECTOR_ELT(res, 0, nn_i);
    SET_VECTOR_ELT(res, 1, nn_d);

    // Add names to list elements
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("indices"));
    SET_STRING_ELT(names, 1, Rf_mkChar("distances"));
    Rf_setAttrib(res, R_NamesSymbol, names);

    // Unprotect R objects
    UNPROTECT(nprot);
    return res;
}



/**
 * @brief Computes k-Nearest Neighbors using the ANN library
 *
 * This function calculates the k-Nearest Neighbors for a given set of points using
 * the Approximate Nearest Neighbor (ANN) library. It constructs a kd-tree for
 * efficient searching and returns both the indices and distances of the k nearest
 * neighbors for each input point.
 *
 * @param RX An R matrix (REALSXP) containing the input points. Each column
 *           represents a point, and each row represents a dimension.
 * @param Rk An R integer (INTSXP) specifying the number of nearest neighbors to find.
 *
 * @return An R list containing two elements:
 *         - indices: An integer matrix (nrX x k) where each column contains the
 *                    indices of the k nearest neighbors for each input point.
 *         - distances: A numeric matrix (nrX x k) where each column contains the
 *                      Euclidean distances to the k nearest neighbors for each input point.
 *
 * @note This function uses the ANN library for efficient nearest neighbor searching.
 *       It constructs a kd-tree from the input points and performs exact (non-approximate)
 *       searches.
 *
 * @Rf_warning The function assumes that the input matrix RX is a valid numeric matrix
 *          and that Rk is a positive integer. No extensive Rf_error checking is performed.
 *
 * @see For more information on the ANN library, visit:
 *      http://www.cs.umd.edu/~mount/ANN/
 */
SEXP S_kNN(SEXP RX, SEXP Rk) {

        // Protect R objects from garbage collection
        int nprot = 0;
        PROTECT(RX = Rf_coerceVector(RX, REALSXP)); nprot++;
        PROTECT(Rk = Rf_coerceVector(Rk, INTSXP)); nprot++;

        double *X = REAL(RX);
        int k = INTEGER(Rk)[0];
        int *dimX = INTEGER(Rf_getAttrib(RX, R_DimSymbol));
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

        // Create R objects to return
        SEXP nn_i = PROTECT(Rf_allocMatrix(INTSXP, nrX, k)); nprot++;
        SEXP nn_d = PROTECT(Rf_allocMatrix(REALSXP, nrX, k)); nprot++;

        // Perform k-NN search for each point
        for (int i = 0; i < nrX; i++) {
                kdTree->annkSearch(dataPts[i], k, nnIdx, nnDist, 0); // 0 for Rf_error bound

                // Copy results to R objects
                for (int j = 0; j < k; j++) {
                        INTEGER(nn_i)[i + nrX * j] = nnIdx[j];
                        REAL(nn_d)[i + nrX * j] = sqrt(nnDist[j]);
                }
        }

        // Clean up
        delete[] nnIdx;
        delete[] nnDist;
        delete kdTree;
        annDeallocPts(dataPts);
        annClose(); // Close ANN

        // Prepare return list
        SEXP res = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++; // List with 2 elements
        SET_VECTOR_ELT(res, 0, nn_i);
        SET_VECTOR_ELT(res, 1, nn_d);

        // Add names to list elements
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
        SET_STRING_ELT(names, 0, Rf_mkChar("indices"));
        SET_STRING_ELT(names, 1, Rf_mkChar("distances"));
        Rf_setAttrib(res, R_NamesSymbol, names);

        // Unprotect R objects
        UNPROTECT(nprot); // RX, Rk, nn_i, nn_d, res, names

        return res;
}
