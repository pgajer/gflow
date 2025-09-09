/*! Fits linear model over the projection of the data within kNN of the given
  point to the first m PCA principal components (PCs), where m is the smallest
  integer such that the m PCs explain at least, say 90\%, of the total variance.
*/

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

// #include <RcppEigen.h>
#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <algorithm>

#include "msr2.h"

extern "C" {
    SEXP S_pclm(SEXP RX,
                SEXP Ry,
                SEXP Rk,
                SEXP Rp,
                SEXP Rsubset);
}

/*!
  Estimates gradient flow trajectories for all points of the given state space; A version with Rsubset parameter !!!

  \param RX       A data matrix.
  \param Ry       A response variable.
  \param Rk       The number of nearest neighbors to use for the estimate of the trajectory.
  \param Rp       The scaling factor for the step in the mean shift algorithm.
  \param Rsubset  A vector of indices of X for which the trajectories are to be found.

  \return Returns a vector of Principla Component Regression estimate of the
  response variable at the center of the set of kNN of a given point (that is
  the center).
*/
SEXP S_pclm(SEXP RX,
            SEXP Ry,
            SEXP Rk,
            SEXP Rp,
            SEXP Rsubset) {

#define DEBUG_S_pclm 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); ++nprot;
    PROTECT(Ry = coerceVector(Ry, REALSXP)); ++nprot;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); ++nprot;
    PROTECT(Rp = coerceVector(Rp, REALSXP)); ++nprot;

    // Check if Rsubset is provided and coerce it to integer vector
    int *subset = NULL;
    int subsetLength = 0;
    if (Rsubset != R_NilValue) {
        PROTECT(Rsubset = coerceVector(Rsubset, INTSXP)); ++nprot;
        subset = INTEGER(Rsubset);
        subsetLength = LENGTH(Rsubset);
    }

    int k = INTEGER(Rk)[0];
    double *y = REAL(Ry);
    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    ANNpointArray dataPts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            dataPts[i][j] = X[i + nrX * j];
        }
    }

#if DEBUG_S_pclm
    Rprintf("S_pclm(): nrX: %d  ncX: %d\n", nrX, ncX);
#endif

    // Creating a kd-tree for X
    ANNkd_tree *kdTree = new ANNkd_tree(dataPts, nrX, ncX);

    // Determine the size of the trajectory list based on Rsubset
    int nPredictions = (subset != NULL) ? subsetLength : nrX;
    //SEXP predictions = PROTECT(allocVector(REALSXP, nPredictions)); ++nprot;

#if DEBUG_S_pclm
    Rprintf("nPredictions: %d\n", nPredictions);
#endif

    // Allocate memory for nearest neighbor indices and distances
    ANNidxArray nn_i = new ANNidx[k];
    ANNdistArray nn_d = new ANNdist[k];

    std::vector<double> nny(k);          // values of y over N(x,k) = the set of kNN inculding x
    std::vector<double> nnX(k * ncX);    // X restricted to N(x,k) - it represents a matrix with k rows and ncX columns;

    for (int i = 0; i < nPredictions; i++) {
        int index = (subset != NULL) ? subset[i] - 1 : i;

        // Check if the index is within valid range
        if (index < 0 || index >= nrX) {
            error("Invalid index in Rsubset");
        }

#if DEBUG_S_pclm
        Rprintf("i: %d  index: %d\n", i, index);
#endif

        // Finding k-NN's of dataPts[i]
        kdTree->annkSearch(dataPts[i], k, nn_i, nn_d, 0); // Identifying kNN's of dataPts[i], including the point dataPts[i].

        // Perform PCA on the subspace span by the vectors x_i - x, where x_i is the i-th nearest neighbor of x = dataPts[i]
        // Identify the smallest number, m, of PC's that explain p proportion of the total variance of nnX

        // Perform linear regression of y versus the projection of nnX onto the first m PCs of the nnX data

        // Construct the nnX matrix (X restricted to kNN)
        for (int j = 0; j < k; j++) {
            int nnIndex = nn_i[j];
            for (int col = 0; col < ncX; col++) {
                nnX[j * ncX + col] = X[nnIndex + nrX * col];
            }
            // Include response values over kNN's
            nny[j] = y[nnIndex];
        }

        // PCA using RcppEigen
#if 0
        Eigen::MatrixXd nnX_eigen(k, ncX);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < ncX; j++) {
                nnX_eigen(i, j) = nnX[i * ncX + j];
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(nnX_eigen, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Variance calculation (using eigenvalues from SVD)
        double totalVariance = svd.singularValues().squaredNorm();
        // ... rest of variance calculation logic ...

        // Projection
        Eigen::MatrixXd projections = nnX_eigen * svd.matrixU().leftCols(m);
        // ... convert projections back to std::vector if needed ...


        // // **PCA**
        // PCA_library::PCA pca(nnX, ncX); // Assuming your PCA library has a 'PCA' class
        // pca.compute();

        // // **Variance Calculation**
        // double totalVariance = pca.getTotalVariance();
        // double cumulativeVariance = 0.0;
        // int m = 0;
        // while (cumulativeVariance < pct * totalVariance && m < k) {
        //     cumulativeVariance += pca.getEigenvalue(m);
        //     m++;
        // }

        // **Projection onto m PCs**
        // std::vector<double> projectedNNX(k * m);
        // pca.project(nnX, projectedNNX, m);

        // **Linear Regression**
        LinReg_library::LinearRegression linReg;
        linReg.fit(projectedNNX, nny, m);

        // **Prediction (Replace with your prediction logic)**
        double xCenter = 0.0; // Calculate the appropriate center representation
        std::vector<double> xCenterProjection = pca.project(xCenter, m); // Project center
        double prediction = linReg.predict(xCenterProjection);

        REAL(predictions)[i] = prediction;
#endif
    }

    // Free memory allocated by ANN
    delete kdTree;
    annDeallocPts(dataPts);

    SEXP res = PROTECT(allocVector(REALSXP, LENGTH(Ry))); nprot++;

    UNPROTECT(nprot);

    return res;
}
