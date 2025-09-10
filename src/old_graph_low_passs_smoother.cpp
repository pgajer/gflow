#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

//#include <omp.h>
#include "omp_compat.h"

#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <cmath>
#include <limits>
// #include <iostream>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <map>
#include <numeric>
#include <random>
#include <chrono>

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "Eigen_utils.h"
#include "msr2.h"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_diffusion_smoother.hpp"
#include "kernels.h"
#include "error_utils.h" // for REPORT_ERROR()

void scale_to_range(std::vector<double>& x, double ymin, double ymax);
double calculate_smallest_difference(const std::vector<double>& y);
std::unique_ptr<std::pair<std::vector<int>, std::vector<int>>> find_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y);
std::unique_ptr<std::vector<double>> prop_nbhrs_with_smaller_y(const std::vector<std::vector<int>>& graph,
                                                               const std::vector<double>& y );
std::unique_ptr<std::vector<double>> cv_imputation(const std::set<int>& test_set,
                                                   const std::vector<std::vector<int>>& graph,
                                                   const std::vector<std::vector<double>>& d,
                                                   const std::vector<double>& y,
                                                   bool y_binary,
                                                   imputation_method_t imputation_method,
                                                   iterative_imputation_params_t iterative_params,
                                                   bool apply_binary_threshold,
                                                   double binary_threshold,
                                                   int ikernel,
                                                   double dist_normalization_factor);

extern "C" {
    SEXP S_graph_spectral_smoother(SEXP Rgraph,
                                   SEXP Rd,
                                   SEXP Rweights,
                                   SEXP Ry,
                                   SEXP Rimputation_method,
                                   SEXP Rmax_iterations,
                                   SEXP Rconvergence_threshold,
                                   SEXP Rapply_binary_threshold,
                                   SEXP Rbinary_threshold,
                                   SEXP Rikernel,
                                   SEXP Rdist_normalization_factor,
                                   SEXP Rn_CVs,
                                   SEXP Rn_CV_folds,
                                   SEXP Repsilon,
                                   SEXP Rmin_plambda,
                                   SEXP Rmax_plambda,
                                   SEXP Rseed);
}



double mean_absolute_deviation(const Eigen::VectorXd& y_t,
                               const Eigen::VectorXd& y) {
    return (y_t - y).cwiseAbs().mean();
}

Eigen::VectorXd compute_y_t(const Eigen::VectorXd& evalues,
                            const Eigen::MatrixXd& evectors,
                            const Eigen::VectorXd& gft,
                            double t,
                            int n_eigenvectors) {
    int n_vertices = gft.size();
    Eigen::VectorXd y_t = Eigen::VectorXd::Zero(n_vertices);
    for (int i = 0; i < n_eigenvectors; ++i) {
        y_t += std::exp(-evalues[i] * t) * gft[i] * evectors.col(i);
    }
    return y_t;
}

double optimize_t(const Eigen::VectorXd& evalues,
                  const Eigen::MatrixXd& evectors,
                  const Eigen::VectorXd& gft,
                  const Eigen::VectorXd& y,
                  int n_eigenvectors) {
    auto objective_function = [&](double t) -> double {
        Eigen::VectorXd y_t = compute_y_t(evalues, evectors, gft, t, n_eigenvectors);
        return mean_absolute_deviation(y_t, y);
    };

    // Nelder-Mead simplex algorithm
    double t0 = 1.0;  // Initial guess for t
    double best_t = t0;
    double best_value = objective_function(t0);

    double step = 0.1;
    double tol = 1e-6;
    int max_iter = 100;

    for (int iter = 0; iter < max_iter; ++iter) {
        double t1 = best_t + step;
        double value1 = objective_function(t1);

        if (value1 < best_value) {
            best_t = t1;
            best_value = value1;
        } else {
            step /= 2.0;
        }

        if (step < tol) {
            break;
        }
    }

    return best_t;
}


/**
 * @brief Applies a graph spectral smoothing algorithm using a low-pass filter
 *        approach with cross-validation to determine the optimal number of eigenvectors.
 *
 * This function smooths a function `y` defined over the vertices of a graph by
 * projecting it onto the subspace spanned by the eigenvectors of the graph Laplacian.
 * It uses sparse Laplacian matrix representation and supports various imputation methods
 * for cross-validation.
 *
 * The algorithm can be summarized as follows:
 *
 * 1. **Input Preparation**:
 *    - Determine if hop index or kernel distance strategy should be used based on `d`.
 *    - Check if `y` is binary and set the binary threshold if necessary.
 *
 * 2. **Graph Construction**:
 *    - Construct the adjacency matrix `A` from the input graph.
 *    - Compute the degree matrix `D` and the Laplacian matrix `L = D - A`.
 *
 * 3. **Eigenvalue Decomposition**:
 *    - Use the Spectra library to compute the smallest eigenvalues and their
 *      corresponding eigenvectors of the Laplacian matrix.
 *
 * 4. **Cross-Validation**:
 *    - Split the vertices into training and test sets for `n_CVs` cross-validation rounds.
 *    - Create a test graph for imputation based on the chosen strategy (hop index or kernel distance).
 *    - Impute values for test vertices using the specified imputation method.
 *    - Apply binary thresholding if necessary.
 *    - Compute the Graph Fourier Transform (GFT) of the imputed function values.
 *    - For each number of eigenvectors, compute the low-pass filtered version and calculate errors.
 *
 * 5. **Optimal Filter Selection**:
 *    - Compute the mean CV Rf_error for each number of eigenvectors.
 *    - Identify the number of eigenvectors that minimizes the mean CV Rf_error.
 *
 * 6. **Final Smoothing**:
 *    - Compute the low-pass filtered version of `y` using the optimal number of eigenvectors.
 *
 * 7. **Result Storage**:
 *    - Store the optimal number of eigenvectors, the smoothed function values,
 *      CV errors, and other relevant information in the result structure.
 *
 * @param graph The adjacency list of the graph.
 * @param d Distances between vertices. If empty, hop index strategy is used.
 * @param weights Weights associated with the vertices (currently unused).
 * @param y Function values defined on the vertices of the graph.
 * @param imputation_method Method for imputing values during cross-validation.
 * @param iterative_params Parameters for iterative imputation methods (only used with ITERATIVE_NEIGHBORHOOD_MATCHING)
 * @param apply_binary_threshold Whether to apply binary thresholding.
 * @param binary_threshold Threshold for binary classification.
 * @param ikernel Integer specifying the kernel function to use.
 * @param dist_normalization_factor A scaling factor applied to the maximum
 * distance between a vertex and its neighbors. This ensures non-zero weights
 * even when all distances are equal, by slightly increasing the normalization
 * range. Default value is 1.01.
 * @param n_CVs Number of cross-validation rounds.
 * @param n_CV_folds Number of folds in each cross-validation round.
 * @param epsilon Small constant for numerical stability in binary case.
 * @param use_low_pass_filter Whether to use low-pass filtering (currently unused).
 * @param preserve_local_extrema Whether to preserve local extrema (currently unused).
 * @param min_plambda Lower bound on the proportion of eigenvectors to use.
 * @param max_plambda Upper bound on the proportion of eigenvectors to use.
 * @param seed Seed for random number generation.
 *
 * @return A unique pointer to a `graph_spectral_smoother_result_t` structure
 *         containing the results of the smoothing.
 *
 * @note The function supports both binary and continuous data, with special
 *       handling for binary cases including different Rf_error calculations.
 */
std::unique_ptr<graph_spectral_smoother_result_t>
graph_spectral_smoother(const std::vector<std::vector<int>>& graph,
                        const std::vector<std::vector<double>>& d,
                        const std::vector<double>& weights,
                        const std::vector<double>& y,
                        imputation_method_t imputation_method = imputation_method_t::LOCAL_MEAN_THRESHOLD,
                        iterative_imputation_params_t iterative_params = {},
                        bool apply_binary_threshold = true,
                        double binary_threshold = 0.5,
                        int ikernel = 1,
                        double dist_normalization_factor = 1.01,
                        int n_CVs = 0,
                        int n_CV_folds = 10,
                        double epsilon = 1e-10,
                        double min_plambda = 0.01, // lower bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.01 that corresponds to 1% of eigenvectors with the smallest eigenvalues
                        double max_plambda = 0.30, // upper bound on the number of eigenvectors to use for smoothing expressed as proportion (of eigenvectors); default is 0.20 that corresponds to 20% of eigenvectors with the smallest eigenvalues
                        unsigned int seed = 0) {

    (void)weights;
    int n_vertices = y.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary && imputation_method == imputation_method_t::GLOBAL_MEAN_THRESHOLD) {
        binary_threshold = mean(y.data(), (int)y.size());
    }

    // Check if a seed was provided
    if (seed == 0) {
        // If no seed was provided, use the current time
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    // Create and seed the random number generator
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> uni(0, n_vertices - 1);

    // Adjacency matrix as a sparse matrix
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> tripletList;
    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        for (int neighbor : graph[vertex]) {
            if (vertex < neighbor) {  // Ensure each edge is added only once
                tripletList.push_back(Eigen::Triplet<double>(vertex, neighbor, 1.0));
                tripletList.push_back(Eigen::Triplet<double>(neighbor, vertex, 1.0));
            }
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // Compute Laplacian matrix as a sparse matrix
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            sum += it.value();
        }
        D.insert(k, k) = sum;
    }
    Eigen::SparseMatrix<double> L = D - A;

    // Eigenvalue decomposition
    Spectra::SparseSymMatProd<double> op(L);

    int min_num_eigenvectors = static_cast<int>(min_plambda * n_vertices);
    int max_num_eigenvectors = static_cast<int>(max_plambda * n_vertices);

    if (min_num_eigenvectors <= 0)
        min_num_eigenvectors = 1;

    if (max_num_eigenvectors < min_num_eigenvectors)
        max_num_eigenvectors = min_num_eigenvectors;

    // Ensure ncv is within bounds: nev < ncv <= n
    int nev = std::min(max_num_eigenvectors + 1, n_vertices);
    int ncv = std::min(2 * nev, n_vertices); // Adjust ncv to be within bounds
    // Ensure nev < ncv
    if (nev >= ncv) {
        nev = ncv - 1;
    }

    if (nev < 1) {
        nev = 1;
    }

    // Construct eigen solver object to find eigenvalues closest to 0
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    eigs.compute(Spectra::SortRule::SmallestAlge);

    if (eigs.info() != Spectra::CompInfo::Successful)
        REPORT_ERROR("Eigenvalue estimation with Spectra failed.");

    Eigen::VectorXd evalues = eigs.eigenvalues();
    Eigen::MatrixXd evectors = eigs.eigenvectors();

    int n_filters = max_num_eigenvectors - min_num_eigenvectors + 1;
    Eigen::MatrixXd cv_errors(n_filters, n_CVs);
    Eigen::MatrixXd low_pass_ys(n_vertices, n_filters);

    if (n_CVs > 0) {

        // Creating a set version of the adjacency matrix of the graph
        std::vector<std::set<int>> set_graph(n_vertices);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            set_graph[vertex].insert(graph[vertex].begin(), graph[vertex].end());
        }

        int fold_size = n_vertices / n_CV_folds;

        //
        // The main cross-validation loop
        //
        for (int cv = 0; cv < n_CVs; ++cv) {
            // Creating a test set
            std::set<int> test_set;
            if (fold_size == 1 && n_vertices == n_CVs) {
                test_set.insert(cv);
            } else {
                while ((int)test_set.size() < fold_size) {
                    int vertex = uni(rng);
                    test_set.insert(vertex);
                }
            }

            std::vector<double> cv_y = std::move(*cv_imputation(test_set,
                                                                graph,
                                                                d,
                                                                y,
                                                                y_binary,
                                                                imputation_method,
                                                                iterative_params,
                                                                apply_binary_threshold,
                                                                binary_threshold,
                                                                ikernel,
                                                                dist_normalization_factor));

            // Turn cv_y into an Eigen vector
            Eigen::VectorXd cv_y_evect = Eigen::Map<const Eigen::VectorXd>(cv_y.data(), cv_y.size());

            // Compute the Graph Fourier Transform (GFT) of cv_y
            Eigen::VectorXd gft_cv_y = evectors.transpose() * cv_y_evect;

            // print_Eigen_VectorXd(gft_cv_y, "gft_cv_y");

            // for each n_eigenvectors construct a low-pass filtered version of cv_y using the first 'n_eigenvectors' eigenvectors
            for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors;
                 n_eigenvectors <= max_num_eigenvectors;
                 ++n_eigenvectors, ++filter_index) {

                Eigen::VectorXd low_pass_cv_y = Eigen::VectorXd::Zero(n_vertices);
                for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
                    low_pass_cv_y += gft_cv_y[i] * evectors.col(i);

                // Computing MAD Rf_error at the test vertices
                double cv_error = 0.0;

                if (y_binary) {
                    for (const auto& vertex : test_set) {
                        double clipped_low_pass_cv_y = std::max(epsilon, std::min(1.0 - epsilon, low_pass_cv_y[vertex]));
                        cv_error += y[vertex] * log(clipped_low_pass_cv_y) + (1 - y[vertex]) * log(1 - clipped_low_pass_cv_y);
                    }
                    cv_error *= -1;
                } else {
                    for (const auto& vertex : test_set) {
                        cv_error += std::abs(low_pass_cv_y[vertex] - y[vertex]);
                    }
                    cv_error /= test_set.size();
                }
                cv_errors(filter_index, cv) = cv_error;
            }

        } // END OF for (int cv = 0; cv < n_CVs; ++cv)
    } // END OF if (n_CVs > 0)

    // computing row-wise means of the matrix cv_errors
    std::vector<double> mean_cv_errors(n_filters, 0.0);
    for (int filter_index = 0; filter_index < n_filters; ++filter_index) {
        for (int cv = 0; cv < n_CVs; ++cv) {
            mean_cv_errors[filter_index] += cv_errors(filter_index, cv); //cv_errors[filter_index + cv * n_filters];
        }
        mean_cv_errors[filter_index] /= n_CVs;
    }

    // finding the low-pass filter index with the smallest mean CV Rf_error
    double min_cv_error = mean_cv_errors[0];
    int opt_filter_index = 0;
    for (int filter_index = 1; filter_index < n_filters; ++filter_index) {
        if (mean_cv_errors[filter_index] < min_cv_error) {
            min_cv_error = mean_cv_errors[filter_index];
            opt_filter_index = filter_index;
        }
    }

    // computing the low-pass filtered version of y for the given opt_filter_index
    // Turning y into an Eigen vector
    Eigen::VectorXd y_evect = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());

    // Computing the Graph Fourier Transform (GFT) of y
    Eigen::VectorXd gft_y = evectors.transpose() * y_evect;

    Eigen::VectorXd opt_low_pass_y = Eigen::VectorXd::Zero(n_vertices);
    int n_eigenvectors = min_num_eigenvectors + opt_filter_index;
    for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
        opt_low_pass_y += gft_y[i] * evectors.col(i);

    for (int filter_index = 0, n_eigenvectors = min_num_eigenvectors; n_eigenvectors <= max_num_eigenvectors; ++n_eigenvectors, ++filter_index) {
        Eigen::VectorXd low_pass_y = Eigen::VectorXd::Zero(n_vertices);
        for (int ev_counter = 0, i = nev - 1; ev_counter < n_eigenvectors; ++ev_counter, --i)
            low_pass_y += gft_y[i] * evectors.col(i);

        for (int vertex = 0; vertex < n_vertices; ++vertex)
            low_pass_ys(vertex, filter_index) = low_pass_y[vertex];
    }

    auto result = std::make_unique<graph_spectral_smoother_result_t>();

    result->optimal_num_eigenvectors = n_eigenvectors;
    result->y_smoothed.assign(opt_low_pass_y.data(), opt_low_pass_y.data() + opt_low_pass_y.size());
    result->cv_errors = cv_errors;
    result->mean_cv_errors = mean_cv_errors;
    result->low_pass_ys = low_pass_ys;
    result->n_filters = n_filters;
    result->evalues = evalues;
    result->evectors = evectors;
    result->min_num_eigenvectors = min_num_eigenvectors;
    result->max_num_eigenvectors = max_num_eigenvectors;

    return result;
}

/**
 * @brief R interface for the graph spectral smoothing algorithm using a low-pass filter approach with cross-validation.
 *
 * This function serves as an interface between R and the C++ graph_spectral_smoother function.
 * It applies a spectral smoothing algorithm to a function defined on the vertices of a graph,
 * using cross-validation to determine the optimal number of eigenvectors for the low-pass filter.
 *
 * The algorithm includes the following main steps:
 * 1. Graph construction and Laplacian computation
 * 2. Eigenvalue decomposition of the Laplacian
 * 3. Cross-validation for determining optimal smoothing parameters
 * 4. Final smoothing of the input function
 *
 * @param Rgraph An R list of integer vectors representing the graph structure (1-based indexing).
 * @param Redge_length An R list of numeric vectors representing distances between vertices.
 *           If empty, a hop index strategy will be used for imputation.
 * @param Rweights An R numeric vector of weights associated with the vertices (currently unused).
 * @param Ry An R numeric vector of function values defined on the vertices of the graph.
 * @param Rimputation_method An R integer specifying the imputation method to use.
 * @param Rmax_iterations The number of iterations in the iterative matching method.
 * @param Rconvergence_threshold The convergence threshold in the iterative matching method.
 * @param Rapply_binary_threshold An R logical indicating whether to apply binary thresholding.
 * @param Rbinary_threshold An R numeric value specifying the threshold for binary classification.
 * @param Rikernel An R integer specifying the kernel function to use.
 * @param Rn_CVs An R integer specifying the number of cross-validation rounds.
 * @param Rn_CV_folds An R integer specifying the number of folds in each cross-validation round.
 * @param Repsilon An R numeric value specifying a small constant for numerical stability.
 * @param Rmin_plambda An R numeric value representing the lower bound on the proportion of eigenvectors to use.
 * @param Rmax_plambda An R numeric value representing the upper bound on the proportion of eigenvectors to use.
 * @param Rseed An R integer used as a seed for the random number generator.
 *
 * @return An R list with the following components:
 *         - evalues: Eigenvalues of the graph Laplacian
 *         - evectors: Eigenvectors of the graph Laplacian
 *         - optimal_num_eigenvectors: The optimal number of eigenvectors determined by cross-validation
 *         - y_smoothed: The smoothed function values
 *         - cv_errors: Matrix of cross-validation errors for each number of eigenvectors and CV round
 *         - mean_cv_errors: Vector of mean cross-validation errors for each number of eigenvectors
 *         - low_pass_ys: Matrix of low-pass filtered versions of y for different numbers of eigenvectors
 *         - min_num_eigenvectors: Minimum number of eigenvectors used
 *         - max_num_eigenvectors: Maximum number of eigenvectors used
 *
 * @note This function converts R objects to C++ types, calls the graph_spectral_smoother function,
 *       and then converts the results back to R objects. It uses PROTECT/UNPROTECT for proper
 *       memory management in R.
 */
SEXP S_graph_spectral_smoother(SEXP Rgraph,
                               SEXP Redge_length,
                               SEXP Rweights,
                               SEXP Ry,
                               SEXP Rimputation_method,
                               SEXP Rmax_iterations,
                               SEXP Rconvergence_threshold,
                               SEXP Rapply_binary_threshold,
                               SEXP Rbinary_threshold,
                               SEXP Rikernel,
                               SEXP Rdist_normalization_factor,
                               SEXP Rn_CVs,
                               SEXP Rn_CV_folds,
                               SEXP Repsilon,
                               SEXP Rmin_plambda,
                               SEXP Rmax_plambda,
                               SEXP Rseed) {

    std::vector<std::vector<int>> graph          = convert_adj_list_from_R(Rgraph);
    std::vector<std::vector<double>> edge_length = convert_weight_list_from_R(Redge_length);
    std::vector<double> weights = std::move(*Rvect_to_CppVect_double(Rweights));
    std::vector<double> y = std::move(*Rvect_to_CppVect_double(Ry));

    imputation_method_t imputation_method = static_cast<imputation_method_t>(INTEGER(Rimputation_method)[0]);
    int max_iterations = INTEGER(Rmax_iterations)[0];
    double convergence_threshold = REAL(Rconvergence_threshold)[0];
    bool apply_binary_threshold = LOGICAL(Rapply_binary_threshold)[0];
    double binary_threshold = REAL(Rbinary_threshold)[0];
    int ikernel = INTEGER(Rikernel)[0];
    double dist_normalization_factor = REAL(Rdist_normalization_factor)[0];
    int n_CVs = INTEGER(Rn_CVs)[0];
    int n_CV_folds = INTEGER(Rn_CV_folds)[0];
    double epsilon = REAL(Repsilon)[0];
    double min_plambda = REAL(Rmin_plambda)[0];
    double max_plambda = REAL(Rmax_plambda)[0];
    unsigned int seed = static_cast<unsigned int>(INTEGER(Rseed)[0]);

    iterative_imputation_params_t iterative_params;
    iterative_params.max_iterations = max_iterations;
    iterative_params.convergence_threshold = convergence_threshold;

    // Call the graph_spectral_smoother function
    std::unique_ptr<graph_spectral_smoother_result_t> result = graph_spectral_smoother(graph,
                                                                                       edge_length,
                                                                                       weights,
                                                                                       y,
                                                                                       imputation_method,
                                                                                       iterative_params,
                                                                                       apply_binary_threshold,
                                                                                       binary_threshold,
                                                                                       ikernel,
                                                                                       dist_normalization_factor,
                                                                                       n_CVs,
                                                                                       n_CV_folds,
                                                                                       epsilon,
                                                                                       min_plambda,
                                                                                       max_plambda,
                                                                                       seed);

    // Convert the results back to R types (this part remains largely the same)
    int nprot = 0;
    SEXP Rresult = PROTECT(Rf_allocVector(VECSXP, 9));  nprot++;

        // evalues
    int evalues_length = result->evalues.size();
    SEXP Revalues = PROTECT(Rf_allocVector(REALSXP, evalues_length)); nprot++;
    for (int i = 0; i < evalues_length; ++i) {
        REAL(Revalues)[i] = result->evalues[i];
    }
    SET_VECTOR_ELT(Rresult, 0, Revalues);

    // evectors
    int evectors_rows = result->evectors.rows();
    int evectors_cols = result->evectors.cols();
    SEXP Revectors = PROTECT(Rf_allocMatrix(REALSXP, evectors_rows, evectors_cols)); nprot++;
    for (int i = 0; i < evectors_rows; ++i) {
        for (int j = 0; j < evectors_cols; ++j) {
            REAL(Revectors)[i + evectors_rows * j] = result->evectors(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 1, Revectors);

    // optimal_num_eigenvectors
    SEXP Roptimal_num_eigenvectors = PROTECT(Rf_allocVector(INTSXP, 1)); nprot++;
    INTEGER(Roptimal_num_eigenvectors)[0] = result->optimal_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 2, Roptimal_num_eigenvectors);

    // y_smoothed
    int y_smoothed_length = result->y_smoothed.size();
    SEXP Ry_smoothed = PROTECT(Rf_allocVector(REALSXP, y_smoothed_length)); nprot++;
    for (int i = 0; i < y_smoothed_length; ++i) {
        REAL(Ry_smoothed)[i] = result->y_smoothed[i];
    }
    SET_VECTOR_ELT(Rresult, 3, Ry_smoothed);

    // cv_errors matrix
    int cv_errors_rows = result->cv_errors.rows();
    int cv_errors_cols = result->cv_errors.cols();
    SEXP Rcv_errors = PROTECT(Rf_allocMatrix(REALSXP, cv_errors_rows, cv_errors_cols)); nprot++;
    for (int i = 0; i < cv_errors_rows; ++i) {
        for (int j = 0; j < cv_errors_cols; ++j) {
            REAL(Rcv_errors)[i + cv_errors_rows * j] = result->cv_errors(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 4, Rcv_errors);

    // mean_cv_errors
    int mean_cv_errors_length = result->mean_cv_errors.size();
    SEXP Rmean_cv_errors = PROTECT(Rf_allocVector(REALSXP, mean_cv_errors_length)); nprot++;
    for (int i = 0; i < mean_cv_errors_length; ++i) {
        REAL(Rmean_cv_errors)[i] = result->mean_cv_errors[i];
    }
    SET_VECTOR_ELT(Rresult, 5, Rmean_cv_errors);

    // low_pass_ys matrix
    int low_pass_ys_rows = result->low_pass_ys.rows();
    int low_pass_ys_cols = result->low_pass_ys.cols();
    SEXP Rlow_pass_ys = PROTECT(Rf_allocMatrix(REALSXP, low_pass_ys_rows, low_pass_ys_cols)); nprot++;
    for (int i = 0; i < low_pass_ys_rows; ++i) {
        for (int j = 0; j < low_pass_ys_cols; ++j) {
            REAL(Rlow_pass_ys)[i + low_pass_ys_rows * j] = result->low_pass_ys(i, j);
        }
    }
    SET_VECTOR_ELT(Rresult, 6, Rlow_pass_ys);

    SEXP Rmin_num_eigenvectors = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
    REAL(Rmin_num_eigenvectors)[0] = result->min_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 7, Rmin_num_eigenvectors);

    SEXP Rmax_num_eigenvectors = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
    REAL(Rmax_num_eigenvectors)[0] = result->max_num_eigenvectors;
    SET_VECTOR_ELT(Rresult, 8, Rmax_num_eigenvectors);

    // Set the names of the list components
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("evalues"));
    SET_STRING_ELT(names, 1, Rf_mkChar("evectors"));
    SET_STRING_ELT(names, 2, Rf_mkChar("optimal_num_eigenvectors"));
    SET_STRING_ELT(names, 3, Rf_mkChar("y_smoothed"));
    SET_STRING_ELT(names, 4, Rf_mkChar("cv_errors"));
    SET_STRING_ELT(names, 5, Rf_mkChar("mean_cv_errors"));
    SET_STRING_ELT(names, 6, Rf_mkChar("low_pass_ys"));
    SET_STRING_ELT(names, 7, Rf_mkChar("min_num_eigenvectors"));
    SET_STRING_ELT(names, 8, Rf_mkChar("max_num_eigenvectors"));
    Rf_setAttrib(Rresult, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return Rresult;
}
