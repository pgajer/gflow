/*!
  Methods for kNN local maxima search and associated MS complex construct
*/

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <algorithm>
#include <memory>

#include "msr2.h"

extern "C" {

    SEXP S_knn_MS_cx(SEXP RX,
                     SEXP Ry,
                     SEXP Rk);

    SEXP S_knn_lmax_search(SEXP RX, SEXP Ry, SEXP Rk);
}

typedef struct
{
  std::vector<int> trajectory;
  int lmax;
} idx_trajectory_t;

/*!
  Estimates PL gradient flow trajectory of a response variable for each point of the give state space

  \param RX             The data matrix.
  \param Ry             The response variable.
  \param Rk             The number of nearest neighbors to use for the estimate of the trajectory.

  \return Returns a vector idx_trajectory_t one for each point of X.
*/
std::unique_ptr<std::vector<idx_trajectory_t>> knn_lmax_search(SEXP RX,
                                                               SEXP Ry,
                                                               SEXP Rk) {
    #define DEBUG_knn_lmax_search 0

    #if DEBUG_knn_lmax_search
    Rprintf("In knn_lmax_search()\n");
    #endif

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *y = REAL(Ry);

    #if DEBUG_knn_lmax_search
    Rprintf("first 10 y:\n");
    for (int i = 0; i < 10; i++)
        Rprintf("%.8f\n",y[i]);
    Rprintf("\n");
    #endif

    int k = INTEGER(Rk)[0];
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // This holds index trajectory structures for all points
    auto res = std::make_unique<std::vector<idx_trajectory_t>>(); // return object
    res->resize(nrX); // Ensure lmax has the correct size

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    // ANNidxArray nn_i = new ANNidx[k];
    std::vector<int> nn_i(k);
    std::vector<double> nn_d(k);

    // The main loop - Finding the PL gradient flow trajectories for all points of X
    for (int pt_i = 0; pt_i < nrX; pt_i++) {

        #if DEBUG_knn_lmax_search
        Rprintf("pt_i: %d\n", pt_i);
        #endif

        // Access the indices of the kNNs for the current point
        //int *nn_i = indices + pt_i * k;
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + nrX * j];
            nn_d[j] = distances[pt_i + nrX * j];
        }

        #if DEBUG_knn_lmax_search
        Rprintf("nn_i:\n");
        print_int_array(nn_i, k);

        Rprintf("nny:\n");
        for (int j = 0; j < k; j++)
            Rprintf("%.8f\n",y[nn_i[j]]);
        Rprintf("\n");
        #endif

        // Finding a neighbor with the maximum positive difference y[nn_i[j]] - y[nn_i[0]]
        #if 0
        double max_delta = 0;
        int max_index = -1;
        double d = 0;
        for (int j = 1; j < k; j++) {
            d = y[nn_i[j]] - y[nn_i[0]];
            if (d > max_delta) {
                max_delta = d;
                max_index = nn_i[j];
            }
        }

        // Finding a __closest__ nearest neighbor with the positive difference y[nn_i[j]] - y[nn_i[0]]
        int first_positive_delta_index = -1;
        double d = 0;
        for (int j = 1; j < k; j++) {
            d = y[nn_i[j]] - y[nn_i[0]];
            if ( d > 0 ) {
                first_positive_delta_index = nn_i[j];
                break;
            }
        }
        #endif

        // Finding the nearest neighbor with the largest positive difference ratio (y[nn_i[j]] - y[nn_i[0]])  / nn_d[j]
        int max_rate_of_change_index = -1;
        double d = 0;
        double max_d = 0;
        for (int j = 1; j < k; j++) {
            d = (y[nn_i[j]] - y[nn_i[0]]) / nn_d[j];
            if ( d > max_d ) {
                max_rate_of_change_index = nn_i[j];
                max_d = d;
                break;
            }
        }

        #if DEBUG_knn_lmax_search
        Rprintf("max_delta: %.8f\n", max_delta);
        Rprintf("max_index: %d\n", max_index);
        #endif

        // Finding remaining indices of the gradient flow of y trajectory for the pt_i-th pt
        #if DEBUG_knn_lmax_search
        int counter = 0;
        #endif
        //while ( max_index > -1 ) {
        //while ( first_positive_delta_index > -1 ) {
        while ( max_rate_of_change_index > -1  ) {

            #if DEBUG_knn_lmax_search
            Rprintf("counter: %d\n", counter);
            counter++;
            #endif

            // (*res)[pt_i].trajectory.push_back(max_index);
            // (*res)[pt_i].lmax = max_index;

            #if 0
            (*res)[pt_i].trajectory.push_back(first_positive_delta_index);
            (*res)[pt_i].lmax = first_positive_delta_index;
            #endif

            (*res)[pt_i].trajectory.push_back(max_rate_of_change_index);
            (*res)[pt_i].lmax = max_rate_of_change_index;

            // Looking for the next point (index) of the trajectory
            //nn_i = indices + max_index * k;
            for (int j = 0; j < k; j++)
                nn_i[j] = indices[max_rate_of_change_index + nrX * j];
                //nn_i[j] = indices[first_positive_delta_index + nrX * j];
                //nn_i[j] = indices[max_index + nrX * j];

            #if 0
            max_delta = 0;
            max_index = -1;
            for (int j = 1; j < k; j++) {
                d = y[nn_i[j]] - y[nn_i[0]];
                if (d > max_delta) {
                    max_delta = d;
                    max_index = nn_i[j];
                }
            }

            first_positive_delta_index = -1;
            for (int j = 1; j < k; j++) {
                d = y[nn_i[j]] - y[nn_i[0]];
                if ( d > 0 ) {
                    first_positive_delta_index = nn_i[j];
                    break;
                }
            }
            #endif

            max_rate_of_change_index = -1;
            d = 0;
            max_d = 0;
            for (int j = 1; j < k; j++) {
                d = (y[nn_i[j]] - y[nn_i[0]]) / nn_d[j];
                if ( d > max_d ) {
                    max_rate_of_change_index = nn_i[j];
                    max_d = d;
                    break;
                }
            }

            #if DEBUG_knn_lmax_search
            Rprintf("max_delta: %.8f\n", max_delta);
            Rprintf("max_index: %d\n", max_index);
            #endif

        } // END Of while (max_index > -1)
    } // END OF for (int pt_i = 0; pt_i < nrX; pt_i++)

    UNPROTECT(nprot);

    return res;
}




/*!
  Estimates PL gradient flow trajectory of a response variable for each point of the give state space

  \param RX             The data matrix.
  \param Ry             The response variable.
  \param Rk             The number of nearest neighbors to use for the estimate of the trajectory.

  \return Returns a list with two components: 1) a list of integer vectors representing the gradient flow trajectories,
  2) an integer vector of local maximum indices for each point of X.
*/
SEXP S_knn_lmax_search(SEXP RX, SEXP Ry, SEXP Rk) {

#define DEBUG_S_knn_lmax_search 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *y = REAL(Ry);

#if DEBUG_S_knn_lmax_search
    Rprintf("first 10 y:\n");
    for (int i = 0; i < 10; i++)
        Rprintf("%.f\n",y[i]);
    Rprintf("\n");
#endif


    int k = INTEGER(Rk)[0];
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Create a list to store the gradient flow trajectories
    SEXP trajectories = PROTECT(allocVector(VECSXP, nrX)); nprot++;

    // Create an integer vector to store the local maximum indices
    SEXP lmax_indices = PROTECT(allocVector(INTSXP, nrX)); nprot++;

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    SEXP indices = VECTOR_ELT(knn_res, 0);

    // -------------------------------------------------------------------------------------
    //
    // The main loop - Finding the PL gradient flow trajectories for all points of X
    //
    // -------------------------------------------------------------------------------------
    for (int pt_i = 0; pt_i < nrX; pt_i++) {
        // Access the indices of the kNNs for the current point
        int* nn_i = INTEGER(VECTOR_ELT(indices, pt_i));

        // Create an integer vector to store the trajectory indices
        SEXP trajectory = PROTECT(allocVector(INTSXP, 0));

        // Finding a neighbor with the maximum positive difference y[nn_i[j]] - y[nn_i[0]]
        double max_delta = 0;
        int max_index = -1;
        double d = 0;
        for (int j = 1; j < k; j++) {
            d = y[nn_i[j]] - y[nn_i[0]];
            if (d > max_delta) {
                max_delta = d;
                max_index = nn_i[j];
            }
        }

        // Finding remaining indices of the gradient flow of y trajectory for the pt_i-th pt
        while (max_index > -1) {
            // Extend trajectory
            int current_length = LENGTH(trajectory);
            trajectory = lengthgets(trajectory, current_length + 1);
            INTEGER(trajectory)[current_length] = max_index + 1; // Convert to 1-based index
            INTEGER(lmax_indices)[pt_i] = max_index + 1; // Store the local maximum index (1-based)

            // Looking for the next point (index) of the trajectory
            nn_i = INTEGER(VECTOR_ELT(indices, max_index));
            max_delta = 0;
            max_index = -1;
            d = 0;
            for (int j = 1; j < k; j++) {
                d = y[nn_i[j]] - y[nn_i[0]];
                if (d > max_delta) {
                    max_delta = d;
                    max_index = nn_i[j];
                }
            }
        } // END of while (max_index > -1)

        SET_VECTOR_ELT(trajectories, pt_i, trajectory);
        UNPROTECT(1);
    } // END of for (int pt_i = 0; pt_i < nrX; pt_i++)

    // Create a list to store the result
    SEXP result = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(result, 0, trajectories);
    SET_VECTOR_ELT(result, 1, lmax_indices);

    // Set names for the list components
    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("trajectories"));
    SET_STRING_ELT(names, 1, mkChar("lmax_indices"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result;
}

/*!
  Finds the Morse-Smale complex of a response variable over the give state space

  \param RX             The data matrix representing points of some state space.
  \param Ry             The response variable.
  \param Rk             The number of nearest neighbors to use for the estimate of the trajectory.

  \return Returns a list with two components 1) a list of trajectoris for all points of X 2) A matrix local maximum associated with all point of X
*/
SEXP S_knn_MS_cx(SEXP RX,
                 SEXP Ry,
                 SEXP Rk) {

#define DEBUG_knn_MS_cx 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double* y = REAL(Ry);
    int nrX = INTEGER(getAttrib(RX, R_DimSymbol))[0];

#if DEBUG_knn_MS_cx
    Rprintf("In S_knn_MS_cx()\n");

    Rprintf("first 10 y:\n");
    for (int i = 0; i < 10; i++)
        Rprintf("%.8f\n",y[i]);
    Rprintf("\n");
#endif

#if 1
#if DEBUG_knn_MS_cx
    Rprintf("Before lmax_res = knn_lmax_search(RX, Ry, Rk)\n");
#endif

    // Gradient flow trajectories for y
    std::unique_ptr<std::vector<idx_trajectory_t>> lmax_res = knn_lmax_search(RX, Ry, Rk);

#if DEBUG_knn_MS_cx
    Rprintf("After lmax_res = knn_lmax_search(RX, Ry, Rk)\n");
    Rprintf("Creating ascending_trajectories ... ");
#endif

    // Create a list to store the ascending trajectories
    SEXP ascending_trajectories = PROTECT(allocVector(VECSXP, nrX)); nprot++;

    // Create an integer vector to store the local maximum indices
    SEXP lmax_indices = PROTECT(allocVector(INTSXP, nrX)); nprot++;

    // Iterate over each idx_trajectory_t structure
    for (int i = 0; i < nrX; i++) {
        const std::vector<int>& trajectory = (*lmax_res)[i].trajectory;
        int lmax = (*lmax_res)[i].lmax;

        // Create an integer vector for the current trajectory
        SEXP trajectory_vec = PROTECT(allocVector(INTSXP, trajectory.size())); nprot++;

        // Copy the trajectory indices to the integer vector
        for (size_t j = 0; j < trajectory.size(); j++) {
            INTEGER(trajectory_vec)[j] = trajectory[j] + 1; // Convert to 1-based index
        }

        // Set the corresponding element of ascending_trajectories
        SET_VECTOR_ELT(ascending_trajectories, i, trajectory_vec);

        // Set the corresponding element of lmax_indices
        INTEGER(lmax_indices)[i] = lmax + 1; // Convert to 1-based index
    }

#if DEBUG_knn_MS_cx
    Rprintf("DONE\n");
#endif

#endif


    //
    // Gradient flow trajectories for -y
    //

#if DEBUG_knn_MS_cx
    Rprintf("Before lmin_res = knn_lmax_search(RX, Rnegy, Rk)\n");
#endif

    // Create a new SEXP for the negative of y
    SEXP Rnegy = PROTECT(allocVector(REALSXP, LENGTH(Ry))); nprot++;
    double* negy = REAL(Rnegy);


    // Copy the negated values of y to negy
    for (int i = 0; i < LENGTH(Ry); i++) {
        negy[i] = -y[i];
    }

#if DEBUG_knn_MS_cx
    Rprintf("first 10 y:\n");
    for (int i = 0; i < 10; i++)
        Rprintf("%.8f\n",y[i]);
    Rprintf("\n");

    Rprintf("first 10 negy:\n");
    for (int i = 0; i < 10; i++)
        Rprintf("%.8f\n",negy[i]);
    Rprintf("\n");

    //error("Stopping function execution for debugging purposes.");
#endif

    // Call S_knn_lmax_search() with the negated y
    std::unique_ptr<std::vector<idx_trajectory_t>> lmin_res = knn_lmax_search(RX, Rnegy, Rk);

#if DEBUG_knn_MS_cx
    Rprintf("After lmin_res = knn_lmax_search(RX, Rnegy, Rk)\n");
    Rprintf("Creating descending_trajectories ... ");
#endif


    // Create a list to store the descending trajectories
    SEXP descending_trajectories = PROTECT(allocVector(VECSXP, nrX)); nprot++;

    // Create an integer vector to store the local minimum indices
    SEXP lmin_indices = PROTECT(allocVector(INTSXP, nrX)); nprot++;

    // Iterate over each idx_trajectory_t structure
    for (int i = 0; i < nrX; i++) {
        const std::vector<int>& trajectory = (*lmin_res)[i].trajectory;
        int lmin = (*lmin_res)[i].lmax;

        // Create an integer vector for the current trajectory
        SEXP trajectory_vec = PROTECT(allocVector(INTSXP, trajectory.size())); nprot++;

        // Copy the trajectory indices to the integer vector
        for (size_t j = 0; j < trajectory.size(); j++) {
            INTEGER(trajectory_vec)[j] = trajectory[j] + 1; // Convert to 1-based index
        }

        // Set the corresponding element of descending_trajectories
        SET_VECTOR_ELT(descending_trajectories, i, trajectory_vec);

        // Set the corresponding element of lmin_indices
        INTEGER(lmin_indices)[i] = lmin + 1; // Convert to 1-based index
    }

#if DEBUG_knn_MS_cx
    Rprintf("DONE\n");
    Rprintf("Creating output list ... ");
#endif


#if 0
    SEXP result = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(result, 0, descending_trajectories);
    SET_VECTOR_ELT(result, 1, lmin_indices);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("descending_trajectories"));
    SET_STRING_ELT(names, 1, mkChar("lmin_basins"));
#endif


#if 1
    // Create a list to store the result
    SEXP result = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(result, 0, ascending_trajectories);
    SET_VECTOR_ELT(result, 1, lmax_indices);
    SET_VECTOR_ELT(result, 2, descending_trajectories);
    SET_VECTOR_ELT(result, 3, lmin_indices);

    // Set names for the list components
    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("ascending_trajectories"));
    SET_STRING_ELT(names, 1, mkChar("lmax_basins"));
    SET_STRING_ELT(names, 2, mkChar("descending_trajectories"));
    SET_STRING_ELT(names, 3, mkChar("lmin_basins"));
    setAttrib(result, R_NamesSymbol, names);
#endif

    UNPROTECT(nprot);

#if DEBUG_knn_MS_cx
    Rprintf("DONE\n");
#endif

    return result;
}
