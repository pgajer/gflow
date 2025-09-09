//
// Functions that generate a Morse-Smale complex given an estimate of a
// conditional expectation Ey of a random variable y over a graph or weighted
// graph associated with some data X.
//

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

std::unique_ptr<std::vector<std::vector<int>>> M_kNN_graph(SEXP RX, SEXP Rk);
std::unique_ptr<std::vector<std::vector<int>>> M_kNN_graph_v2(SEXP RX, SEXP Rk);


extern "C" {
    SEXP S_M_kNN_graph_MS_cx(SEXP REy, SEXP RX, SEXP Rk, SEXP Rversion);
}

/*!
  Creates a Morse-Smale Complex Using a Mutual kNN Graph

  This function analyzes a dataset using a kNN graph and a conditional
  expectation estimate to construct the Morse-Smale complex of the
  conditional expectation.

  \param REy An estimate of the conditional expectation of a random variable y over X.
  \param RX A matrix where each row represents a data point.
  \param Rk The number of nearest neighbors for the kNN graph construction.
  \param Rversion (Unused at present)

  \return A list with two components:
  * **lext:**  An nrX-by-2 matrix of local extrema. Each row contains the
  indices of the local minimum and local maximum associated with the
  corresponding data point.
  * **trajectories:** A list of gradient flow trajectories. Each element is a
  vector of data point indices representing the path from a local
  minimum to a local maximum.

  \details For each point, the function finds ascending and descending gradient
  flow trajectories based on the conditional expectation estimate. These
  trajectories are used to identify local minima and maxima, which form the
  basis of the Morse-Smale complex representation.
*/
SEXP S_M_kNN_graph_MS_cx(SEXP REy, SEXP RX, SEXP Rk, SEXP Rversion) {

    int version = INTEGER(Rversion)[0];
    if (version < 1 || version > 2) {
        Rf_error("Invalid 'Rversion' value. Must be 1 or 2");
    }

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    PROTECT(REy = coerceVector(REy, REALSXP)); nprot++;
    double *Ey = REAL(REy);

    std::unique_ptr<std::vector<std::vector<int>>> graph_res;
    if (version == 1) {
        graph_res = M_kNN_graph(RX, Rk);
        if (!graph_res) Rf_error("M_kNN_graph function returned an invalid result (null pointer).");
    } else if (version == 2) {
        graph_res = M_kNN_graph_v2(RX, Rk);
        if (!graph_res) Rf_error("M_kNN_graph_v2 function returned an invalid result (null pointer).");
    }

    // Creating a list of gradient flow trajectories
    SEXP trajectories = PROTECT(allocVector(VECSXP, nrX)); nprot++;

    // Creating a list of the indices of local minimum and local maximum of each point
    SEXP lmin_lmax_matrix = PROTECT(allocMatrix(INTSXP, nrX, 2)); nprot++;
    int *lmin_lmax_matrix_ptr = INTEGER(lmin_lmax_matrix);

    int max_diff_i = 0;  // index of the maximum difference w
    double max_diff = 0; // value of the maximum difference
    double d;
    int lmax_i;          // local maximum index of the current point
    int lmin_i;          // local minimum index of the current point
    for (int i = 0; i < nrX; i++) {

        // Finding ascending Ey gradient trajectory of the i-th point of X
        std::vector<int> ascending_trajectory;
        max_diff_i = i;
        while ( max_diff_i > -1 ) {
            ascending_trajectory.push_back(max_diff_i);
            // Finding the edge with the larges positive difference Ey[neighbor] - Ey[i]
            max_diff_i = -1;
            max_diff = 0;
            for (int neighbor : graph_res->at(i)) {
                d = Ey[neighbor] - Ey[i];
                if ( d > max_diff ) {
                    max_diff = d;
                    max_diff_i = neighbor;
                }
            }
        }
        lmax_i = max_diff_i;

        // Finding descending Ey gradient trajectory of the i-th point of X
        std::vector<int> descending_trajectory; // At the end I need to
        // remove the first
        // element i from
        // descending_trajectory
        // as it is already in
        // ascending_trajectory
        max_diff_i = i;
        while ( max_diff_i > -1 ) {
            ascending_trajectory.push_back(max_diff_i);
            // Finding the edge with the larges positive difference Ey[i] - Ey[neighbor]
            max_diff_i = -1;
            max_diff = 0;
            for (int neighbor : graph_res->at(i)) {
                d = Ey[i] - Ey[neighbor];
                if ( d > max_diff ) {
                    max_diff = d;
                    max_diff_i = neighbor;
                }
            }
        }
        lmin_i = max_diff_i;

        // Adding lmin_i and lmax_i to lmin_lmax_matrix_ptr
        lmin_lmax_matrix_ptr[i] = lmin_i;
        lmin_lmax_matrix_ptr[i + nrX] = lmax_i;

        // Trajectory
        int tr_size = ascending_trajectory.size() + descending_trajectory.size() - 1; // removing i from the descending_trajectory trajectory
        SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
        int n = descending_trajectory.size() - 1;
        for (int j = n; j > 0; j--)
            REAL(trajectory)[n - j] = descending_trajectory[j];

        int ascending_trajectory_size = ascending_trajectory.size();
        for (int j = 0; j < ascending_trajectory_size; j++)
            REAL(trajectory)[j + n] = ascending_trajectory[j];

        SET_VECTOR_ELT(trajectories, i, trajectory);
        UNPROTECT(1);

    } // END OF for (int i = 0; i < nrX; i++)

    // Creating a return list
    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++; // List with 2 elements
    SET_VECTOR_ELT(res, 0, lmin_lmax_matrix);
    SET_VECTOR_ELT(res, 1, trajectories);

    // Add names to list elements
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("lext"));
    SET_STRING_ELT(names, 1, mkChar("trajectories"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
