// Standard C++ library headers
#include <vector>
#include <unordered_map>
#include <cmath>        // For std::exp
#include <algorithm>    // For std::min
#include <cstdio>       // For Rprintf

#include <Eigen/Core>
#include <Eigen/Dense>  // For Eigen::MatrixXd
#include <Eigen/Sparse> // For Eigen::SparseMatrix, Triplet
#include <Spectra/SymEigsSolver.h>          // For SymEigsSolver
#include <Spectra/MatOp/SparseSymMatProd.h> // For SparseSymMatProd

// Project-specific headers
#include "set_wgraph.hpp"    // For set_wgraph_t struct
#include "edge_info.hpp"     // For edge_info_t if needed
#include "error_utils.h"     // For REPORT_ERROR macro

/**
 * @brief Computes the eigenvectors of the graph Laplacian.
 *
 * @details Constructs the Laplacian matrix from the adjacency list and edge weights of the graph.
 * Performs eigen decomposition to obtain the first `n_evectors` eigenvectors corresponding to the
 * smallest eigenvalues, which are useful for spectral embedding or smoothing tasks.
 *
 * @param n_evectors Number of eigenvectors to compute (excluding the trivial constant eigenvector).
 * @param verbose Whether to print progress and diagnostic information.
 * @return Matrix of eigenvectors. Each column corresponds to one eigenvector.
 */
Eigen::MatrixXd set_wgraph_t::compute_graph_laplacian_eigenvectors(
    size_t n_evectors,
    bool verbose
    ) const {

    size_t n_vertices = adjacency_list.size();

    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;
    triples.reserve(n_vertices * 2);

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            if (i < j) {
                triples.emplace_back(i, j, edge.weight);
                triples.emplace_back(j, i, edge.weight);
            }
        }
    }
    A.setFromTriplets(triples.begin(), triples.end());

    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            sum += it.value();
        }
        D.insert(k, k) = sum;
    }

    Eigen::SparseMatrix<double> L = D - A;
    constexpr double lap_reg = 1e-8;
    for (size_t i = 0; i < n_vertices; i++) {
        L.coeffRef(i, i) += lap_reg;  // Small regularization
    }

    Spectra::SparseSymMatProd<double> op(L);

    int nev = std::min<int>(2 * n_evectors + 5, n_vertices);
    int ncv = std::min<int>(2 * nev, n_vertices);  // Control parameter for the algorithm

    // Ensure nev < ncv
    if (nev >= ncv) {
        ncv = nev + 1;
    }

    // ncv rule-of-thumb: 3*nev is often safer for hard problems
    int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
    ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv_default);
    eigs.init();
    int maxit = 1000;
    double tol = 1e-6;
    eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

    Eigen::MatrixXd eigenvectors;

    if (eigs.info() != Spectra::CompInfo::Successful) {
        // Define fallback parameters to try
        std::vector<std::pair<int, double>> attempts = {
            {2000, 1e-8},    // More iterations, slightly relaxed tolerance
            {3000, 1e-6},    // Even more iterations, more relaxed tolerance
            {5000, 1e-4}     // Final attempt with very relaxed parameters
        };
        // Try with original ncv first
        bool success = false;
        for (const auto& params : attempts) {
            int adjusted_maxit = params.first;
            double adjusted_tol = params.second;
            eigs.init();
            eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                if (verbose) {
                    Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                            ncv, adjusted_maxit, adjusted_tol);
                }
                eigenvectors = eigs.eigenvectors();  // Add this line to extract eigenvectors
                success = true;
                break;
            }
        }
        // If still not successful, try with increased ncv values
        if (!success) {
            int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
            std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

            for (const auto& multiplier : ncv_multipliers) {
                int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
                // Create a new solver with adjusted ncv
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
                for (const auto& params : attempts) {
                    int adjusted_maxit = params.first;
                    double adjusted_tol = params.second;
                    adjusted_eigs.init();
                    adjusted_eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
                    if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
                        if (verbose) {
                            Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                                    adjusted_ncv, adjusted_maxit, adjusted_tol);
                        }
                        eigenvectors = adjusted_eigs.eigenvectors();
                        success = true;
                        break;
                    }
                }
                if (success) {
                    break;
                }
            }
        }
        // If all attempts failed, report an Rf_error
        if (!success) {
            REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
        }
    } else {
        // Get eigenvectors
        eigenvectors = eigs.eigenvectors();
    }

    return eigenvectors;
}

/**
 * @brief Computes the eigenvectors of the odd-powered shifted kernel graph Laplacian.
 *
 * @details Constructs the kernel Laplacian matrix using a Gaussian-like kernel based on
 * graph geodesic distances, then computes (I - L)^k, and performs eigen decomposition.
 * It preserves the fallback strategy for eigenvalue computations when convergence issues arise.
 *
 * @param n_evectors Number of eigenvectors to compute (excluding the trivial constant eigenvector).
 * @param tau Kernel bandwidth controlling neighborhood size (small positive value).
 * @param k Odd positive integer representing the power to which (I - L) is raised.
 * @param verbose Whether to print progress and diagnostic information.
 * @return Matrix of eigenvectors. Each column corresponds to one eigenvector.
 */
Eigen::MatrixXd set_wgraph_t::compute_graph_shifted_kernel_laplacian_eigenvectors(
    size_t n_evectors,
    double tau,
    size_t k,
    bool verbose
    ) const {
    if (k % 2 == 0 || k == 0) {
        REPORT_ERROR("k must be a positive odd integer.");
    }

    size_t n_vertices = adjacency_list.size();

    // Step 1: Build Kernel Weight Matrix K
    Eigen::SparseMatrix<double> K(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;

    for (size_t i = 0; i < n_vertices; ++i) {
        reachability_map_t rmap = compute_graph_reachability_map(i, tau);
        for (const auto& [j, dist] : rmap.distances) {
            if (i < j) {  // To avoid duplicating symmetric entries
                double weight = std::exp(-(dist * dist) / (tau * tau));
                triples.emplace_back(i, j, weight);
                triples.emplace_back(j, i, weight);
            }
        }
    }
    K.setFromTriplets(triples.begin(), triples.end());

    // Step 2: Build Degree Matrix D
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (size_t i = 0; i < (size_t)K.outerSize(); ++i) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it) {
            sum += it.value();
        }
        if (sum > 0) {
            D.insert(i, i) = sum;
        }
    }

    // Step 3: Compute Kernel Laplacian L = D - K
    Eigen::SparseMatrix<double> L = D - K;

    // Step 4: Regularize L
    constexpr double lap_reg = 1e-8;
    for (size_t i = 0; i < n_vertices; ++i) {
        L.coeffRef(i, i) += lap_reg;
    }

    // Step 5: Form (I - tau * L)
    Eigen::SparseMatrix<double> I_minus_L(n_vertices, n_vertices);
    I_minus_L.setIdentity();
    I_minus_L -= L;

    // Step 6: Compute (I - tau * L)^k
    Eigen::SparseMatrix<double> shifted_kernel_matrix = I_minus_L;
    for (size_t power = 1; power < k; ++power) {
        shifted_kernel_matrix = shifted_kernel_matrix * I_minus_L;
    }

    // Step 7: Eigen decomposition
    Spectra::SparseSymMatProd<double> op(shifted_kernel_matrix);

    int nev = std::min<int>(2 * n_evectors + 5, n_vertices);
    int ncv = std::min<int>(2 * nev, n_vertices);

    if (nev >= ncv) {
        ncv = nev + 1;
    }

    // ncv rule-of-thumb: 3*nev is often safer for hard problems
    int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
    ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv_default);
    eigs.init();
    int maxit = 5000;
    double tol = 1e-6;
    eigs.compute(Spectra::SortRule::LargestAlge, maxit, tol);

    Eigen::MatrixXd eigenvectors;

    if (eigs.info() != Spectra::CompInfo::Successful) {
        // Define fallback parameters to try
        std::vector<std::pair<int, double>> attempts = {
            {2000, 1e-8},
            {3000, 1e-6},
            {5000, 1e-4}
        };

        bool success = false;
        for (const auto& params : attempts) {
            int adjusted_maxit = params.first;
            double adjusted_tol = params.second;
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                if (verbose) {
                    Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                            ncv, adjusted_maxit, adjusted_tol);
                }
                eigenvectors = eigs.eigenvectors();
                success = true;
                break;
            }
        }

        if (!success) {
            int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
            std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

            for (const auto& multiplier : ncv_multipliers) {
                int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
                for (const auto& params : attempts) {
                    int adjusted_maxit = params.first;
                    double adjusted_tol = params.second;
                    adjusted_eigs.init();
                    adjusted_eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
                    if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
                        if (verbose) {
                            Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                                    adjusted_ncv, adjusted_maxit, adjusted_tol);
                        }
                        eigenvectors = adjusted_eigs.eigenvectors();
                        success = true;
                        break;
                    }
                }
                if (success) {
                    break;
                }
            }
        }

        if (!success) {
            REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
        }
    }
    else {
        eigenvectors = eigs.eigenvectors();
    }

    return eigenvectors;
}


/**
 * @brief Computes the spectrum (eigenvalues and eigenvectors) of the shifted kernel graph Laplacian.
 *
 * @details Constructs the kernel Laplacian matrix using a Gaussian-like kernel based on
 * graph geodesic distances, then computes (I - L)^k, and performs eigen decomposition.
 * It preserves the fallback strategy for eigenvalue computations when convergence issues arise.
 *
 * @param n_evectors Number of eigenvectors to compute (excluding the trivial constant eigenvector).
 * @param tau Kernel bandwidth controlling neighborhood size (small positive value).
 * @param laplacian_power Odd positive integer representing the power to which (I - L) is raised.
 * @param verbose Whether to print progress and diagnostic information.
 *
 * @return A std::pair consisting of:
 *   - first:  Eigen::VectorXd of eigenvalues (sorted ascending),
 *   - second: Eigen::MatrixXd whose columns are the corresponding eigenvectors.
 */
std::pair<Eigen::VectorXd, Eigen::MatrixXd>
set_wgraph_t::compute_graph_shifted_kernel_laplacian_spectrum(
    size_t n_evectors,
    double tau,
    double radius_factor,
    size_t laplacian_power,
    bool verbose
    ) const {
    // if (laplacian_power % 2 == 0 || laplacian_power == 0) {
    //     REPORT_ERROR("laplacian_power must be a positive odd integer.");
    // }

    (void)radius_factor;

    size_t n_vertices = adjacency_list.size();

    // Step 1: Build Kernel Weight Matrix K
    Eigen::SparseMatrix<double> K(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;

    double tau_squared = tau * tau;
    (void)tau_squared;
    for (size_t i = 0; i < n_vertices; ++i) {

        for (const auto& edge : adjacency_list[i]) {
            if(i < edge.vertex) {
                triples.emplace_back(i, edge.vertex, edge.weight);
                triples.emplace_back(edge.vertex, i, edge.weight);
            }
        }

#if 0
        reachability_map_t rmap = compute_graph_reachability_map(i, radius_factor * tau);
        for (const auto& [j, dist] : rmap.distances) {
            if (i < j) {  // To avoid duplicating symmetric entries
                double weight = std::exp(-(dist * dist) / tau_squared);
                triples.emplace_back(i, j, weight);
                triples.emplace_back(j, i, weight);
            }
        }
#endif
    }
    K.setFromTriplets(triples.begin(), triples.end());

    // Step 2: Build Degree Matrix D
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (size_t i = 0; i < (size_t)K.outerSize(); ++i) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, i); it; ++it) {
            sum += it.value();
        }
        if (sum > 0) {
            D.insert(i, i) = sum;
        }
    }

    // Step 3: Compute Kernel Laplacian L = D - (tau * K)
    // Eigen::SparseMatrix<double> L = D - (tau * K);
    Eigen::SparseMatrix<double> L = D - K;

    // Step 4: Regularize L
    constexpr double lap_reg = 1e-8;
    for (size_t i = 0; i < n_vertices; ++i) {
        L.coeffRef(i, i) += lap_reg;
    }

    // Step 5: Form (I - L)
    Eigen::SparseMatrix<double> I_minus_L(n_vertices, n_vertices);
    I_minus_L.setIdentity();
    I_minus_L -= L;

    // Step 6: Compute (I - L)^laplacian_power
    Eigen::SparseMatrix<double> shifted_kernel_matrix = I_minus_L;
    for (size_t power_index = 1; power_index < laplacian_power; ++power_index) {
        shifted_kernel_matrix = shifted_kernel_matrix * I_minus_L;
    }

    // Step 7: Eigen decomposition
    Spectra::SparseSymMatProd<double> op(shifted_kernel_matrix);

    int nev = std::min<int>(2 * n_evectors + 5, n_vertices);
    int ncv = std::min<int>(2 * nev, n_vertices);

    if (nev >= ncv) {
        ncv = nev + 1;
    }

    // ncv rule-of-thumb: 3*nev is often safer for hard problems
    int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
    ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv_default);
    eigs.init();
    int maxit = 5000;
    double tol = 1e-6;
    eigs.compute(Spectra::SortRule::LargestAlge, maxit, tol);

    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;

    if (eigs.info() != Spectra::CompInfo::Successful) {
        // Define fallback parameters to try
        std::vector<std::pair<int, double>> attempts = {
            {2000, 1e-8},
            {3000, 1e-6},
            {5000, 1e-4}
        };

        bool success = false;
        for (const auto& params : attempts) {
            int adjusted_maxit = params.first;
            double adjusted_tol = params.second;
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                if (verbose) {
                    Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                            ncv, adjusted_maxit, adjusted_tol);
                }
                eigenvectors = eigs.eigenvectors();
                eigenvalues = eigs.eigenvalues();
                success = true;
                break;
            }
        }

        if (!success) {
            int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
            std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

            for (const auto& multiplier : ncv_multipliers) {
                int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
                for (const auto& params : attempts) {
                    int adjusted_maxit = params.first;
                    double adjusted_tol = params.second;
                    adjusted_eigs.init();
                    adjusted_eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
                    if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
                        if (verbose) {
                            Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                                    adjusted_ncv, adjusted_maxit, adjusted_tol);
                        }
                        eigenvectors = adjusted_eigs.eigenvectors();
                        eigenvalues = adjusted_eigs.eigenvalues();
                        success = true;
                        break;
                    }
                }
                if (success) {
                    break;
                }
            }
        }

        if (!success) {
            REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
        }
    }
    else {
        eigenvectors = eigs.eigenvectors();
        eigenvalues = eigs.eigenvalues();
    }

    // Return eigenvalues and eigenvectors as a pair
    return {eigenvalues, eigenvectors};
}


std::pair<Eigen::VectorXd, Eigen::MatrixXd>
set_wgraph_t::compute_graph_laplacian_spectrum_generic(
    const std::vector<double>& y,
    size_t n_evectors,
    laplacian_type_t laplacian_type,
    const kernel_params_t& kernel_params,
    size_t power,
    bool verbose
    ) const {
    size_t n_vertices = adjacency_list.size();

    // Step 1: Construct the appropriate Laplacian matrix
    Eigen::SparseMatrix<double> L(n_vertices, n_vertices);
    bool is_shifted = false;

    if (verbose) {
        Rprintf("Computing %s Laplacian spectrum...\n",
                laplacian_type_to_string(laplacian_type).c_str());
    }

    switch (laplacian_type) {
    case laplacian_type_t::STANDARD:
        L = construct_standard_laplacian(kernel_params);
        break;

    case laplacian_type_t::NORMALIZED:
        L = construct_normalized_laplacian(kernel_params);
        break;

    case laplacian_type_t::RANDOM_WALK:
        L = construct_random_walk_laplacian(kernel_params);
        break;

    case laplacian_type_t::KERNEL:
        L = construct_kernel_laplacian(kernel_params);
        break;

    case laplacian_type_t::NORMALIZED_KERNEL:
        L = construct_normalized_kernel_laplacian(kernel_params);
        break;

    case laplacian_type_t::ADAPTIVE_KERNEL:
        L = construct_adaptive_kernel_laplacian(kernel_params);
        break;

    case laplacian_type_t::SHIFTED:
        L = construct_standard_laplacian(kernel_params);
        is_shifted = true;
        break;

    case laplacian_type_t::SHIFTED_KERNEL:
        L = construct_kernel_laplacian(kernel_params);
        is_shifted = true;
        break;

    case laplacian_type_t::REGULARIZED:
        L = construct_regularized_laplacian(kernel_params, 1e-8);
        break;

    case laplacian_type_t::REGULARIZED_KERNEL:
        L = construct_regularized_kernel_laplacian(kernel_params, 1e-8);
        break;

    case laplacian_type_t::MULTI_SCALE:
        L = construct_multi_scale_laplacian(kernel_params);
        break;

    case laplacian_type_t::PATH:
        L = construct_path_laplacian(y, kernel_params, verbose);
        break;

    default:
        REPORT_ERROR("Unsupported Laplacian type");
    }

    // Step 2: Apply shift if needed
    Eigen::SparseMatrix<double> op_matrix;
    if (is_shifted) {
        // Create I - L
        Eigen::SparseMatrix<double> I(n_vertices, n_vertices);
        I.setIdentity();
        op_matrix = I - L;
    } else {
        op_matrix = L;
    }

    // Step 3: Power the matrix if needed (we'll handle this implicitly through eigenvalues)

    // Step 4: Compute eigendecomposition
    Spectra::SparseSymMatProd<double> op(op_matrix);

    // Setup Spectra parameters with additional validation
    // First, ensure n_evectors doesn't exceed matrix size
    size_t max_evectors = std::min(n_evectors, n_vertices - 1);

    // Calculate nev and ncv with bounds checking
    int nev = std::min<int>(max_evectors, n_vertices - 1);
    // Add a small buffer for convergence but don't exceed n_vertices
    int ncv = std::min<int>(2 * nev + 5, n_vertices);

    // Ensure nev < ncv and both are valid
    if (nev >= ncv) {
        ncv = nev + 1;
    }

    // Final safety check
    if (ncv > (int)n_vertices) {
        ncv = n_vertices;
        nev = std::max(1, static_cast<int>(ncv - 1));
    }

    if (verbose) {
        Rprintf("Eigendecomposition parameters: nev=%d, ncv=%d, n_vertices=%zu\n",
                nev, ncv, n_vertices);
    }

    // Create solver and compute
    // ncv rule-of-thumb: 3*nev is often safer for hard problems
    int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
    ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv_default);

    eigs.init();
    int maxit = 1000;
    double tol = 1e-10;

    // Sort based on whether we're using shifted or standard Laplacian
    Spectra::SortRule sort_rule = is_shifted ?
        Spectra::SortRule::LargestAlge :
        Spectra::SortRule::SmallestAlge;

    eigs.compute(sort_rule, maxit, tol);

    // Get eigendecomposition
    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;

    // Handle computation failures with fallback strategies
    if (eigs.info() != Spectra::CompInfo::Successful) {
        // Define fallback parameters to try
        std::vector<std::pair<int, double>> attempts = {
            {2000, 1e-8},
            {3000, 1e-6},
            {5000, 1e-4}
        };

        bool success = false;
        for (const auto& params : attempts) {
            int adjusted_maxit = params.first;
            double adjusted_tol = params.second;
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                if (verbose) {
                    Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                            ncv, adjusted_maxit, adjusted_tol);
                }
                eigenvectors = eigs.eigenvectors();
                eigenvalues = eigs.eigenvalues();
                success = true;
                break;
            }
        }

        if (!success) {
            int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
            std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

            for (const auto& multiplier : ncv_multipliers) {
                int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
                for (const auto& params : attempts) {
                    int adjusted_maxit = params.first;
                    double adjusted_tol = params.second;
                    adjusted_eigs.init();
                    adjusted_eigs.compute(Spectra::SortRule::LargestAlge, adjusted_maxit, adjusted_tol);
                    if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
                        if (verbose) {
                            Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                                    adjusted_ncv, adjusted_maxit, adjusted_tol);
                        }
                        eigenvectors = adjusted_eigs.eigenvectors();
                        eigenvalues = adjusted_eigs.eigenvalues();
                        success = true;
                        break;
                    }
                }
                if (success) {
                    break;
                }
            }
        }

        if (!success) {
            REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
        }

    } else {
        eigenvectors = eigs.eigenvectors();
        eigenvalues = eigs.eigenvalues();
    }

    // If power > 1, apply power to eigenvalues
    if (power > 1) {
        for (size_t i = 0; i < (size_t)eigenvalues.size(); ++i) {
            eigenvalues(i) = std::pow(eigenvalues(i), power);
        }
    }

    if (verbose) {
        Rprintf("Computed %d eigenpairs (requested %zu)\n",
                (int)eigenvalues.size(), n_evectors);
    }

    return {eigenvalues, eigenvectors};
}


/**
 * @brief Constructs the standard combinatorial Laplacian matrix L = D - A
 *
 * @return Eigen::SparseMatrix<double> The standard Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_standard_laplacian(
    const kernel_params_t& params
    ) const {

    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 2); // Rough estimate for size

    // Construct adjacency matrix A and degree matrix D simultaneously
    std::vector<double> degrees(n_vertices, 0.0);

    // Validate graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    double tau = params.tau_factor * graph_diameter;

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            // Compute kernel weight based on distance
            double weight = compute_kernel_weight(edge.weight, tau, params.kernel_type);

            // Add to adjacency matrix (only once per edge)
            if (i < j) {
                triplets.emplace_back(i, j, -weight);
                triplets.emplace_back(j, i, -weight);
            }

            // Accumulate degree
            degrees[i] += weight;
        }
    }

    // Add diagonal degree entries
    for (size_t i = 0; i < n_vertices; ++i) {
        triplets.emplace_back(i, i, degrees[i]);
    }

    // Construct the sparse Laplacian matrix
    Eigen::SparseMatrix<double> L(n_vertices, n_vertices);
    L.setFromTriplets(triplets.begin(), triplets.end());
    L.makeCompressed();

    return L;
}

/**
 * @brief Constructs the normalized Laplacian matrix L_norm = D^(-1/2) L D^(-1/2)
 *
 * @param params Parameters controlling kernel weights
 * @return Eigen::SparseMatrix<double> The normalized Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_normalized_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 2);

    // Validate graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    double tau = params.tau_factor * graph_diameter;

    // First compute degrees and store connection weights
    std::vector<double> degrees(n_vertices, 0.0);
    std::vector<std::tuple<size_t, size_t, double>> connections;

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            // Compute kernel weight based on distance
            double weight = compute_kernel_weight(edge.weight, tau, params.kernel_type);

            // Store the connection if it's the first time we see it
            if (i < j) {
                connections.emplace_back(i, j, weight);
            }

            // Accumulate degree
            degrees[i] += weight;
        }
    }

    // Compute D^(-1/2)
    std::vector<double> d_neg_half(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        d_neg_half[i] = (degrees[i] > 0) ? 1.0 / std::sqrt(degrees[i]) : 0.0;
    }

    // Construct normalized Laplacian
    for (size_t i = 0; i < n_vertices; ++i) {
        // Diagonal: 1 for any vertex with non-zero degree
        if (degrees[i] > 0) {
            triplets.emplace_back(i, i, 1.0);
        }
    }

    // Off-diagonal: -w_ij / sqrt(d_i * d_j)
    for (const auto& [i, j, weight] : connections) {
        if (degrees[i] > 0 && degrees[j] > 0) {
            double norm_weight = -weight * d_neg_half[i] * d_neg_half[j];
            triplets.emplace_back(i, j, norm_weight);
            triplets.emplace_back(j, i, norm_weight);
        }
    }

    // Construct the sparse normalized Laplacian matrix
    Eigen::SparseMatrix<double> L_norm(n_vertices, n_vertices);
    L_norm.setFromTriplets(triplets.begin(), triplets.end());
    L_norm.makeCompressed();

    return L_norm;
}

/**
 * @brief Constructs the random walk Laplacian matrix L_rw = D^(-1) L
 *
 * @param params Parameters controlling kernel weights
 * @return Eigen::SparseMatrix<double> The random walk Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_random_walk_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 2);

    // Validate graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    double tau = params.tau_factor * graph_diameter;

    // First compute degrees
    std::vector<double> degrees(n_vertices, 0.0);

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            // Compute kernel weight based on distance
            double weight = compute_kernel_weight(edge.weight, tau, params.kernel_type);

            // Accumulate degree
            degrees[i] += weight;
        }
    }

    // Construct random walk Laplacian
    for (size_t i = 0; i < n_vertices; ++i) {
        // Diagonal: 1 for any vertex with non-zero degree
        if (degrees[i] > 0) {
            triplets.emplace_back(i, i, 1.0);
        }

        // Off-diagonal: -w_ij / d_i
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            double dist = edge.weight;
            double weight = compute_kernel_weight(dist, tau, params.kernel_type);

            if (degrees[i] > 0) {
                double norm_weight = -weight / degrees[i];
                triplets.emplace_back(i, j, norm_weight);
            }
        }
    }

    // Construct the sparse random walk Laplacian matrix
    Eigen::SparseMatrix<double> L_rw(n_vertices, n_vertices);
    L_rw.setFromTriplets(triplets.begin(), triplets.end());
    L_rw.makeCompressed();

    return L_rw;
}


/**
 * @brief Constructs a kernel-based Laplacian matrix using distance-based kernel weights
 *
 * @param params Parameters controlling kernel construction
 * @return Eigen::SparseMatrix<double> The kernel Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_kernel_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 4); // Estimate - may need more for dense graphs

    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    double min_radius = params.min_radius_factor * graph_diameter;
    double max_radius = params.max_radius_factor * graph_diameter;
    double tau = params.tau_factor * graph_diameter;
    double radius = params.radius_factor * tau;

    // Compute kernel weights
    std::vector<double> degrees(n_vertices, 0.0);

    for (size_t i = 0; i < n_vertices; ++i) {

        // Making sure we have at least params.domain_min_size vertices in the disk of given radius
        double vertex_min_radius = find_minimum_radius_for_domain_min_size(
            i,
            min_radius,
            max_radius,
            params.domain_min_size,
            params.precision
            );

        if (vertex_min_radius > radius) {
            radius = vertex_min_radius;
        }

        // Find vertices within radius
        std::unordered_map<size_t, double> nearby = find_vertices_within_radius(i, radius);

        for (const auto& [j, dist] : nearby) {
            // Only process each edge once
            if (i < j) {
                // Compute kernel weight based on distance
                double weight = compute_kernel_weight(dist, tau, params.kernel_type);

                if (weight > 0) {
                    // Add to adjacency matrix
                    triplets.emplace_back(i, j, -weight);
                    triplets.emplace_back(j, i, -weight);

                    // Accumulate degree
                    degrees[i] += weight;
                    degrees[j] += weight;
                }
            }
        }
    }

    // Add diagonal degree entries
    for (size_t i = 0; i < n_vertices; ++i) {
        triplets.emplace_back(i, i, degrees[i]);
    }

    // Construct the sparse kernel Laplacian matrix
    Eigen::SparseMatrix<double> L_kernel(n_vertices, n_vertices);
    L_kernel.setFromTriplets(triplets.begin(), triplets.end());
    L_kernel.makeCompressed();

    return L_kernel;
}

/**
 * @brief Constructs a normalized kernel Laplacian matrix
 *
 * @param params Parameters controlling kernel construction
 * @return Eigen::SparseMatrix<double> The normalized kernel Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_normalized_kernel_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    // First compute the kernel weights and degrees
    double min_radius = params.min_radius_factor * graph_diameter;
    double max_radius = params.max_radius_factor * graph_diameter;
    double tau        = params.tau_factor * graph_diameter;
    double radius     = params.radius_factor * tau;

    // Store kernel weights and degrees
    std::vector<std::tuple<size_t, size_t, double>> kernel_weights;
    std::vector<double> degrees(n_vertices, 0.0);

    for (size_t i = 0; i < n_vertices; ++i) {

        // Making sure we have at least params.domain_min_size vertices in the disk of given radius
        double vertex_min_radius = find_minimum_radius_for_domain_min_size(
            i,
            min_radius,
            max_radius,
            params.domain_min_size,
            params.precision
            );

        if (vertex_min_radius > radius) {
            radius = vertex_min_radius;
        }

        // Find vertices within radius
        std::unordered_map<size_t, double> nearby = find_vertices_within_radius(i, radius);

        for (const auto& [j, dist] : nearby) {
            // Only process each edge once
            if (i < j) {
                // Compute kernel weight based on distance
                double weight = compute_kernel_weight(dist, tau, params.kernel_type);

                if (weight > 0) {
                    // Store the weight
                    kernel_weights.emplace_back(i, j, weight);

                    // Accumulate degree
                    degrees[i] += weight;
                    degrees[j] += weight;
                }
            }
        }
    }

    // Compute D^(-1/2)
    std::vector<double> d_neg_half(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        d_neg_half[i] = (degrees[i] > 0) ? 1.0 / std::sqrt(degrees[i]) : 0.0;
    }

    // Create triplets for normalized Laplacian
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices + 2 * kernel_weights.size());

    // Diagonal entries (all 1)
    for (size_t i = 0; i < n_vertices; ++i) {
        if (degrees[i] > 0) {
            triplets.emplace_back(i, i, 1.0);
        }
    }

    // Off-diagonal entries
    for (const auto& [i, j, weight] : kernel_weights) {
        double norm_weight = -weight * d_neg_half[i] * d_neg_half[j];
        triplets.emplace_back(i, j, norm_weight);
        triplets.emplace_back(j, i, norm_weight);
    }

    // Construct the sparse normalized kernel Laplacian matrix
    Eigen::SparseMatrix<double> L_norm_kernel(n_vertices, n_vertices);
    L_norm_kernel.setFromTriplets(triplets.begin(), triplets.end());
    L_norm_kernel.makeCompressed();

    return L_norm_kernel;
}

/**
 * @brief Constructs an adaptive kernel Laplacian with locally varying bandwidths
 *
 * @param params Base parameters for kernel construction
 * @return Eigen::SparseMatrix<double> The adaptive kernel Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_adaptive_kernel_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 4);

    if (graph_diameter <= 0) {
        REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    }

    // First compute the kernel weights and degrees
    double min_radius = params.min_radius_factor * graph_diameter;
    double max_radius = params.max_radius_factor * graph_diameter;
    double base_tau = params.tau_factor * graph_diameter;
    double radius = params.radius_factor * base_tau;

    // Compute local density to adapt bandwidth
    std::vector<double> local_density(n_vertices, 0.0);
    std::vector<size_t> neighbor_counts(n_vertices, 0);

    // First pass: count neighbors and compute average distance
    std::vector<double> vertex_radius(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        // Making sure we have at least params.domain_min_size vertices in the disk of given radius
        double vertex_min_radius = find_minimum_radius_for_domain_min_size(
            i,
            min_radius,
            max_radius,
            params.domain_min_size,
            params.precision
            );

        if (vertex_min_radius > radius) {
            radius = vertex_min_radius;
        }
        vertex_radius[i] = radius;

        std::unordered_map<size_t, double> nearby = find_vertices_within_radius(i, radius);

        neighbor_counts[i] = nearby.size();
        if (!nearby.empty()) {
            double avg_dist = 0.0;
            for (const auto& [j, dist] : nearby) {
                avg_dist += dist;
            }
            local_density[i] = avg_dist / nearby.size();
        }
    }

    // Compute adaptive tau for each vertex
    std::vector<double> adaptive_tau(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        // Higher density (lower average distance) -> smaller tau
        // Lower density (higher average distance) -> larger tau
        if (neighbor_counts[i] > 0 && local_density[i] > 0) {
            adaptive_tau[i] = base_tau * std::sqrt(local_density[i]);
        } else {
            adaptive_tau[i] = base_tau;
        }
    }

    // Compute kernel weights with adaptive bandwidth
    std::vector<double> degrees(n_vertices, 0.0);

    for (size_t i = 0; i < n_vertices; ++i) {
        std::unordered_map<size_t, double> nearby = find_vertices_within_radius(i, vertex_radius[i]);

        for (const auto& [j, dist] : nearby) {
            if (i < j) {
                // Use average of the two adaptive bandwidths
                double tau_ij = (adaptive_tau[i] + adaptive_tau[j]) / 2.0;

                // Compute kernel weight based on distance and adaptive bandwidth
                double weight = compute_kernel_weight(dist, tau_ij, params.kernel_type);

                if (weight > 0) {
                    // Add to adjacency matrix
                    triplets.emplace_back(i, j, -weight);
                    triplets.emplace_back(j, i, -weight);

                    // Accumulate degree
                    degrees[i] += weight;
                    degrees[j] += weight;
                }
            }
        }
    }

    // Add diagonal degree entries
    for (size_t i = 0; i < n_vertices; ++i) {
        triplets.emplace_back(i, i, degrees[i]);
    }

    // Construct the sparse adaptive kernel Laplacian matrix
    Eigen::SparseMatrix<double> L_adaptive(n_vertices, n_vertices);
    L_adaptive.setFromTriplets(triplets.begin(), triplets.end());
    L_adaptive.makeCompressed();

    return L_adaptive;
}

/**
 * @brief Constructs a regularized Laplacian L + ε*I to ensure positive definiteness
 *
 * @param epsilon Regularization parameter (small positive value)
 * @return Eigen::SparseMatrix<double> The regularized Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_regularized_laplacian(
    const kernel_params_t& params,
    double epsilon
    ) const {

    // Start with the standard Laplacian
    Eigen::SparseMatrix<double> L = construct_standard_laplacian(params);

    // Add regularization to diagonal
    for (size_t k = 0; k < (size_t)L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            if (it.row() == it.col()) {
                it.valueRef() += epsilon;
            }
        }
    }

    return L;
}

/**
 * @brief Constructs a regularized kernel Laplacian L_kernel + ε*I
 *
 * @param params Parameters for kernel construction
 * @param epsilon Regularization parameter (small positive value)
 * @return Eigen::SparseMatrix<double> The regularized kernel Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_regularized_kernel_laplacian(
    const kernel_params_t& params,
    double epsilon) const {

    // Start with the kernel Laplacian
    Eigen::SparseMatrix<double> L_kernel = construct_kernel_laplacian(params);

    // Add regularization to diagonal
    for (size_t k = 0; k < (size_t)L_kernel.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L_kernel, k); it; ++it) {
            if (it.row() == it.col()) {
                it.valueRef() += epsilon;
            }
        }
    }

    return L_kernel;
}

/**
 * @brief Constructs a multi-scale Laplacian that combines kernel Laplacians at different scales
 *
 * @param params Base parameters for kernel construction
 * @return Eigen::SparseMatrix<double> The multi-scale Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_multi_scale_laplacian(const kernel_params_t& params) const {
    size_t n_vertices = adjacency_list.size();

    // Define a series of scales (bandwidths)
    // For example, use 3 scales: small, medium, large

    // if (graph_diameter <= 0) {
    // REPORT_ERROR("ERROR: graph_diameter has to be greater than 0!");
    // }

    // double tau = params.tau_factor * graph_diameter;

    std::vector<double> scales = {params.tau_factor * 0.5, params.tau_factor, params.tau_factor * 2.0};
    std::vector<double> weights = {0.25, 0.5, 0.25}; // Weights for each scale

    // Initialize the combined Laplacian
    Eigen::SparseMatrix<double> L_multi(n_vertices, n_vertices);
    L_multi.setZero();

    // Compute and combine Laplacians at different scales
    for (size_t s = 0; s < scales.size(); ++s) {
        kernel_params_t scale_params = params;
        scale_params.tau_factor = scales[s];

        // Compute Laplacian at this scale
        Eigen::SparseMatrix<double> L_scale = construct_kernel_laplacian(scale_params);

        // Add to combined Laplacian with appropriate weight
        L_multi += weights[s] * L_scale;
    }

    L_multi.makeCompressed();
    return L_multi;
}

/**
 * @brief Creates a composite path for a graph with only a single geodesic ray
 *
 * This function handles the special case where only one geodesic ray is found from a reference vertex.
 * It creates a composite path by mirroring the single ray through the reference vertex, allowing
 * for meaningful spectral analysis on simple path graphs and chains.
 *
 * The function performs these steps:
 * 1. Creates a mirror path by reversing the original path (excluding the reference vertex)
 * 2. Appends the reference vertex and the original path to form a complete mirror path
 * 3. Calculates appropriate distance values, with negative distances for the mirrored segment
 *    and positive distances for the original segment
 * 4. Adds the resulting composite path to the paths collection and registers it with a special marker
 *
 * @param ref_vertex The reference vertex that serves as the mirroring point
 * @param shortest_paths The original shortest_paths_t structure containing a single path
 * @param composite_paths The composite_shortest_paths_t structure to which the mirrored path will be added
 * @param stats Statistics structure to track composite geodesic creation
 *
 * @note This function assumes that shortest_paths contains exactly one path
 * @note The mirrored path uses INVALID_VERTEX as a special marker in the composite_paths structure
 * @note Distances along the mirrored path are negative before the reference vertex and positive after it
 *
 * @see composite_shortest_paths_t Structure for tracking composite geodesic paths
 * @see get_edge_weight() For calculating distances between adjacent vertices
 */
void set_wgraph_t::handle_single_ray_geodesic(
    shortest_paths_t& shortest_paths,
    composite_shortest_paths_t& composite_paths
    ) const {

    // Get the single path
    const auto& path = shortest_paths.paths[0];

    // Create a mirror path
    std::vector<size_t> mirror_path;

    // First add reversed path (excluding reference vertex)
    for (int i = path.vertices.size() - 1; i > 0; --i) {
        mirror_path.push_back(path.vertices[i]);
    }

    // Then add reference vertex and forward path
    mirror_path.insert(mirror_path.end(), path.vertices.begin(), path.vertices.end());

    // Create an artificial composite path
    path_t composite_path;
    composite_path.vertices = mirror_path;

    // Compute distances along the mirror path
    composite_path.distances.resize(mirror_path.size());
    composite_path.ref_vertex_index = path.vertices.size() - 1; // Position of the reference vertex

    // Compute distances (negative for mirror part, positive for forward part)
    double cumulative_dist = 0.0;

    // Forward distances from ref_vertex
    for (size_t i = composite_path.ref_vertex_index + 1; i < mirror_path.size(); ++i) {
        size_t u = mirror_path[i-1];
        size_t v = mirror_path[i];
        double edge_dist = get_edge_weight(u, v);
        cumulative_dist += edge_dist;
        composite_path.distances[i] = cumulative_dist;
    }

    // Backward distances from ref_vertex
    cumulative_dist = 0.0;
    for (int i = composite_path.ref_vertex_index - 1; i >= 0; --i) {
        size_t u = mirror_path[i];
        size_t v = mirror_path[i+1];
        double edge_dist = get_edge_weight(u, v);
        cumulative_dist -= edge_dist;
        composite_path.distances[i] = cumulative_dist;
    }

    // Add to composite paths using a special marker
    shortest_paths.paths.push_back(composite_path);
    composite_paths.add_composite_shortest_path(0, INVALID_VERTEX); // Special marker for mirror pathb
}

/**
 * @brief Computes weights for second derivative estimation using a mirror-extended path approach
 *
 * @details
 * This function estimates weights for approximating the second derivative at a reference vertex
 * using a specially constructed mirror-extended path. This approach is designed for boundary vertices
 * where only a single ray geodesic exists, creating a mirror reflection of the path to handle
 * the boundary condition appropriately.
 *
 * The function implements a constrained quadratic model with the following properties:
 * - Centers the function at the reference vertex by shifting all y values
 * - Uses a model of the form y(t) = βt + γt²/2 where t=0 at reference vertex
 * - Enforces odd symmetry for the linear term (βt) across the reference point
 * - Enforces even symmetry for the quadratic term (γt²/2) across the reference point
 *
 * This symmetry handling ensures that the second derivative estimate remains stable and
 * mathematically coherent at boundary vertices, which would otherwise suffer from
 * directional bias and instability.
 *
 * @param[in] path Vector of vertex indices representing the path (reference vertex should be in this path)
 * @param[in] center_vertex The reference vertex at which to estimate the second derivative
 * @param[in] y Vector of function values for all vertices in the graph
 * @param[in] params Kernel parameters controlling bandwidth and other kernel properties
 *
 * @return Vector of weights for each vertex in the path. When these weights are applied to the
 *         function values, they approximate the second derivative at the reference vertex.
 *
 * @note This function is specifically designed for boundary vertices where traditional
 *       quadratic fitting methods would fail due to asymmetric neighborhood information.
 *
 * @see compute_quadratic_weights For the standard approach used with interior vertices
 * @see handle_single_ray_geodesic For the function that constructs the mirror-extended path
 *
 * @pre The center_vertex must be present in the path
 * @pre y must have values for all vertices in the path
 * @pre path must be a mirror-extended path with center_vertex at the mirror point
 */
std::vector<double> set_wgraph_t::compute_mirror_quadratic_weights(
    const std::vector<size_t>& path,
    size_t center_vertex,
    const std::vector<double>& y,
    const kernel_params_t& params
    ) const {

    size_t n = path.size();
    std::vector<double> weights(n, 0.0);

    // Find center vertex position
    size_t center_idx = std::find(path.begin(), path.end(), center_vertex) - path.begin();

    // Compute distances along path from center
    std::vector<double> distances(n, 0.0);

    // Forward distances
    double cumulative_dist = 0.0;
    for (size_t i = center_idx + 1; i < n; ++i) {
        double edge_weight = get_edge_weight(path[i-1], path[i]);
        cumulative_dist += edge_weight;
        distances[i] = cumulative_dist;
    }

    // Backward distances
    cumulative_dist = 0.0;
    for (int i = center_idx - 1; i >= 0; --i) {
        double edge_weight = get_edge_weight(path[i], path[i+1]);
        cumulative_dist += edge_weight;
        distances[i] = -cumulative_dist;  // Note the negative sign for mirror part
    }

    // Extract y values and apply the shift
    std::vector<double> shifted_y(n);
    double center_y = y[center_vertex];

    for (size_t i = 0; i < n; ++i) {
        shifted_y[i] = y[path[i]] - center_y;
    }

    // Compute kernel weights
    std::vector<double> kernel_weights(n);
    double max_dist = std::max(
        std::abs(*std::min_element(distances.begin(), distances.end())),
        std::abs(*std::max_element(distances.begin(), distances.end()))
        );
    double bandwidth = params.tau_factor * max_dist;

    for (size_t i = 0; i < n; ++i) {
        kernel_weights[i] = compute_kernel_weight(
            std::abs(distances[i]),
            bandwidth,
            params.kernel_type
            );
    }

    // Now create design matrix for the constrained model
    // Here's the key insight: For a function f(t) with f(0) = 0,
    // if f(t) = βt + γt²/2 for t > 0,
    // and we want f(-t) = -f(t) for t < 0 (odd function symmetry for β term)
    // and f(-t) = f(t) for t < 0 (even function symmetry for γ term)
    // Then our design matrix needs to account for this

    Eigen::MatrixXd X(n, 2);
    for (size_t i = 0; i < n; ++i) {
        // For t term: use sign(t) to account for odd symmetry
        X(i, 0) = distances[i];  // Will be negative for mirror part automatically

        // For t² term: always positive (even symmetry)
        X(i, 1) = distances[i] * distances[i] / 2.0;
    }

    // Weighted least squares with shifted_y
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(n);
    W.diagonal() = Eigen::Map<const Eigen::VectorXd>(kernel_weights.data(), n);

    Eigen::VectorXd y_eigen = Eigen::Map<const Eigen::VectorXd>(shifted_y.data(), n);

    Eigen::MatrixXd XtWX = X.transpose() * W * X;
    XtWX.diagonal() += Eigen::VectorXd::Constant(2, 1e-10);  // Regularization

    Eigen::VectorXd coeffs = XtWX.inverse() * (X.transpose() * W * y_eigen);

    // Extract weights for second derivative (the coefficient of t²/2)
    // γ = coeffs[1], so y'' = 2γ
    // The weights are: 2 * W * X * (X^T W X)^{-1} * e2
    Eigen::VectorXd e2 = Eigen::VectorXd::Zero(2);
    e2(1) = 1.0;

    Eigen::VectorXd weight_vector = 2.0 * W * X * XtWX.inverse() * e2;

    // Copy weights to output vector
    for (size_t i = 0; i < n; ++i) {
        weights[i] = weight_vector(i);
    }

    return weights;
}

/**
 * @brief Constructs a path Laplacian matrix
 *
 * @param params Parameters controlling kernel weights construction
 * @return Eigen::SparseMatrix<double> The path Laplacian matrix
 */
Eigen::SparseMatrix<double>
set_wgraph_t::construct_path_laplacian(
    const std::vector<double>& y,
    const kernel_params_t& params,
    bool verbose
    ) const {

    size_t n_vertices = adjacency_list.size();

    // Create triplets for sparse matrix construction
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 10); // Rough estimate

    // Validate graph diameter
    if (graph_diameter <= 0) {
        const_cast<set_wgraph_t*>(this)->compute_graph_diameter();
        if (graph_diameter <= 0) {
            REPORT_ERROR("ERROR: graph_diameter must be greater than 0!");
        }
    }

    // Set up radius parameters
    double min_radius = params.min_radius_factor * graph_diameter;
    double max_radius = params.max_radius_factor * graph_diameter;
    double tau = params.tau_factor * graph_diameter;
    double radius = params.radius_factor * tau;

    // Structure to track geodesic statistics if needed
    struct geodesic_stats_t {
        size_t total_paths = 0;
        size_t composite_geodesics = 0;
        size_t total_vertices_processed = 0;
    };
    geodesic_stats_t stats;

    // Process each vertex
    for (size_t i = 0; i < n_vertices; ++i) {
        if (verbose && i % 100 == 0) {
            Rprintf("Processing vertex %zu/%zu\n", i, n_vertices);
        }

        // Ensure sufficient vertices in the disk
        double vertex_radius = find_minimum_radius_for_domain_min_size(
            i,
            min_radius,
            max_radius,
            params.domain_min_size,
            params.precision
            );

        if (vertex_radius > radius) {
            radius = vertex_radius;
        }

        // Find shortest paths within radius
        shortest_paths_t shortest_paths = find_graph_paths_within_radius_and_path_min_size(i, radius, params.domain_min_size);

        composite_shortest_paths_t composite_shortest_paths(shortest_paths);

        if (shortest_paths.paths.size() > 1) {
            // Calculate composite geodesics (symmetric construction)
            for (size_t j = 0; j < shortest_paths.paths.size(); j++) {
                for (size_t k = j + 1; k < shortest_paths.paths.size(); k++) {
                    // Check both directions for symmetry
                    if (is_composite_path_geodesic(j, k, shortest_paths)) {
                        composite_shortest_paths.add_composite_shortest_path(j, k);
                        stats.composite_geodesics++;
                    }
                    if (is_composite_path_geodesic(k, j, shortest_paths)) {
                        composite_shortest_paths.add_composite_shortest_path(k, j);
                        stats.composite_geodesics++;
                    }
                }
            }
        } else if (shortest_paths.paths.size() == 1) {
            // Special handling for single ray geodesic
            // Use mirror path approach
            handle_single_ray_geodesic(shortest_paths, composite_shortest_paths);
        } else {
            REPORT_ERROR("ERROR: No single ray geodesic found for vertex: %zu\n", i);
        }

        stats.total_paths += shortest_paths.paths.size();

        // Weight accumulator for vertex i
        std::unordered_map<size_t, double> weight_accumulator;
        size_t valid_composite_paths = 0;

        // Process each composite geodesic path
        for (const auto& [j, k] : composite_shortest_paths.composite_paths) {
            if (k == INVALID_VERTEX) {
                // Handle mirror path case
                // Use the single-ray mirror path that was added to the paths collection. More specifically, here is what has happened in this case:

                // 1) Initially, when shortest_paths.paths.size() == 1, there's only one path (the original single ray path)
                // 2) The handle_single_ray_geodesic() function adds a new mirror path to shortest_paths.paths
                //    This makes shortest_paths.paths.size() == 2
                // 3) Therefore, the mirror path is at index 1, which is also shortest_paths.paths.size() - 1

                // Get the mirror path
                size_t mirror_path_idx = shortest_paths.paths.size() - 1;
                const auto& mirror_path = shortest_paths.paths[mirror_path_idx];

                // Compute weights using mirror quadratic weights approach
                std::vector<double> path_weights = compute_mirror_quadratic_weights(
                    mirror_path.vertices, i, y, params
                    );

                // Accumulate weights
                for (size_t v = 0; v < mirror_path.vertices.size(); ++v) {
                    size_t vertex = mirror_path.vertices[v];
                    weight_accumulator[vertex] += path_weights[v];
                }

                valid_composite_paths++;
                continue;
            }

            // Regular composite path handling (unchanged)
            // Construct composite path (forward direction)
            std::vector<size_t> composite_path;
            // double composite_path_length = 0;

            // First part: path j
            const auto& path1 = composite_shortest_paths.paths[j];
            composite_path.insert(composite_path.end(),
                                  path1.vertices.begin(),
                                  path1.vertices.end());
            // composite_path_length += path1.total_weight;

            // Second part: path k in reverse (excluding Rf_duplicate vertex i)
            const auto& path2 = composite_shortest_paths.paths[k];
            for (int idx = path2.vertices.size() - 1; idx >= 0; idx--) {
                if (path2.vertices[idx] != i) {
                    composite_path.push_back(path2.vertices[idx]);
                }
            }
            // composite_path_length += path2.total_weight;

            // Compute weights for second derivative along composite path
            std::vector<double> path_weights = compute_quadratic_weights(
                composite_path, i, params
                );

            // Accumulate weights
            for (size_t v = 0; v < composite_path.size(); ++v) {
                size_t vertex = composite_path[v];
                weight_accumulator[vertex] += path_weights[v];
            }

            valid_composite_paths++;
        }

        // Normalize weights and add to triplets
        if (valid_composite_paths > 0) {
            double normalization = 1.0 / valid_composite_paths;

            for (const auto& [j, weight] : weight_accumulator) {
                double normalized_weight = weight * normalization;
                if (std::abs(normalized_weight) > 1e-12) {
                    triplets.emplace_back(i, j, normalized_weight);
                }
            }
        }

        stats.total_vertices_processed++;
    }

    if (verbose) {
        Rprintf("Path Laplacian construction completed:\n");
        Rprintf("  Total paths explored: %zu\n", stats.total_paths);
        Rprintf("  Composite geodesics found: %zu\n", stats.composite_geodesics);
        Rprintf("  Vertices processed: %zu\n", stats.total_vertices_processed);
    }

    // Construct sparse matrix
    Eigen::SparseMatrix<double> L_path(n_vertices, n_vertices);
    L_path.setFromTriplets(triplets.begin(), triplets.end());
    L_path.makeCompressed();

    if (verbose) {
        Eigen::MatrixXd L_dense = Eigen::MatrixXd(L_path);
        double asymmetry = (L_dense - L_dense.transpose()).norm();
        if (asymmetry > 1e-10) {
            REPORT_WARNING("Warning: Path Laplacian asymmetry: %e\n", asymmetry);
        }
    }

    // Always symmetrize (using proper sparse operation)
    {
        Eigen::SparseMatrix<double> L_trans = L_path.transpose();
        std::vector<Eigen::Triplet<double>> sym_triplets;
        sym_triplets.reserve(L_path.nonZeros() * 2);

        // Iterate through both matrices and create a half-sum
        for (int k = 0; k < L_path.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L_path, k); it; ++it) {
                sym_triplets.emplace_back(it.row(), it.col(), it.value() * 0.5);
            }

            for (Eigen::SparseMatrix<double>::InnerIterator it(L_trans, k); it; ++it) {
                sym_triplets.emplace_back(it.row(), it.col(), it.value() * 0.5);
            }
        }

        // Create symmetrized matrix
        Eigen::SparseMatrix<double> L_sym(n_vertices, n_vertices);
        L_sym.setFromTriplets(sym_triplets.begin(), sym_triplets.end());
        L_sym.makeCompressed();

        // Replace L_path with symmetrized version
        L_path = L_sym;
    }

    return L_path;
}

// Helper function to compute quadratic weights
std::vector<double>
set_wgraph_t::compute_quadratic_weights(
    const std::vector<size_t>& path,
    size_t center_vertex,
    const kernel_params_t& params
    ) const {

    size_t n = path.size();
    std::vector<double> weights(n, 0.0);

    // Find center vertex position in path
    size_t center_idx = std::find(path.begin(), path.end(), center_vertex) - path.begin();
    if (center_idx >= n) {
        REPORT_ERROR("Center vertex not found in path");
    }

    // Compute distances along path from center
    std::vector<double> distances(n, 0.0);

    // Forward distances
    double cumulative_dist = 0.0;
    for (size_t i = center_idx + 1; i < n; ++i) {
        double edge_weight = get_edge_weight(path[i-1], path[i]);
        cumulative_dist += edge_weight;
        distances[i] = cumulative_dist;
    }

    // Backward distances
    cumulative_dist = 0.0;
    for (int i = center_idx - 1; i >= 0; --i) {
        double edge_weight = get_edge_weight(path[i], path[i+1]);
        cumulative_dist += edge_weight;
        distances[i] = -cumulative_dist;
    }

    // Compute kernel bandwidth
    double max_dist = std::max(
        std::abs(*std::min_element(distances.begin(), distances.end())),
        std::abs(*std::max_element(distances.begin(), distances.end()))
        );
    double bandwidth = params.tau_factor * max_dist;

    // Compute kernel weights
    std::vector<double> kernel_weights(n);
    for (size_t i = 0; i < n; ++i) {
        kernel_weights[i] = compute_kernel_weight(
            std::abs(distances[i]),
            bandwidth,
            params.kernel_type
            );
    }

    // Design matrix for quadratic fit
    Eigen::MatrixXd X(n, 3);
    for (size_t i = 0; i < n; ++i) {
        X(i, 0) = 1.0;
        X(i, 1) = distances[i];
        X(i, 2) = distances[i] * distances[i] / 2.0;
    }

    // Weight matrix
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(n);
    W.diagonal() = Eigen::Map<const Eigen::VectorXd>(kernel_weights.data(), n);

    // Weighted least squares solution
    Eigen::MatrixXd XtWX = X.transpose() * W * X;
    Eigen::MatrixXd XtW = X.transpose() * W;

    // Add regularization for numerical stability
    XtWX.diagonal() += Eigen::VectorXd::Constant(3, 1e-10);

    // Solve for coefficients
    Eigen::MatrixXd inv_XtWX = XtWX.inverse();

    // Extract weights for second derivative (third coefficient)
    Eigen::VectorXd e3 = Eigen::VectorXd::Zero(3);
    e3(2) = 1.0;

    Eigen::VectorXd weight_vector = W * X * inv_XtWX * e3;

    // Copy weights to output vector
    for (size_t i = 0; i < n; ++i) {
        weights[i] = weight_vector(i);
    }

    return weights;
}

// Helper function to get edge weight between two vertices
double
set_wgraph_t::get_edge_weight(size_t u, size_t v) const {
    ensure_edge_weights_computed();

    // Create ordered pair for symmetric lookup
    std::pair<size_t, size_t> edge_key = (u < v) ?
        std::make_pair(u, v) : std::make_pair(v, u);

    auto it = edge_weights.find(edge_key);
    if (it != edge_weights.end()) {
        return it->second;
    }

    // If edge doesn't exist, return infinity
    return std::numeric_limits<double>::infinity();
}

/**
 * @brief Computes kernel weight based on distance and kernel type
 *
 * @param distance Distance between vertices
 * @param tau Kernel bandwidth parameter
 * @param kernel_type Type of kernel function to use
 * @return double The computed kernel weight
 */
double
set_wgraph_t::compute_kernel_weight(
    double distance,
    double tau,
    kernel_type_t kernel_type) const {

    switch (kernel_type) {
    case kernel_type_t::S_INVERSE:
        if (distance < 1e-10) return 1.0 / (1e-10 * tau);
        return 1.0 / (distance * tau);

    case kernel_type_t::S_GAUSSIAN:
        return std::exp(-(distance * distance) / (tau * tau));

    case kernel_type_t::S_EXPONENTIAL:
        return std::exp(-distance / tau);

    case kernel_type_t::S_HEAT:
        return std::exp(-(distance * distance) / (4 * tau));

    case kernel_type_t::S_TRICUBE: {
        double ratio = distance / tau;
        if (ratio < 1.0) {
            double tmp = 1.0 - ratio * ratio * ratio;
            return tmp * tmp * tmp;
        }
        return 0.0;
    }

    case kernel_type_t::S_EPANECHNIKOV: {
        double ratio = distance / tau;
        if (ratio < 1.0) {
            return 1.0 - ratio * ratio;
        }
        return 0.0;
    }

    case kernel_type_t::S_UNIFORM:
        return (distance < tau) ? 1.0 : 0.0;

    case kernel_type_t::S_TRIANGULAR: {
        double ratio = distance / tau;
        if (ratio < 1.0) {
            return 1.0 - ratio;
        }
        return 0.0;
    }

    case kernel_type_t::S_QUARTIC: {
        double ratio = distance / tau;
        if (ratio < 1.0) {
            double tmp = 1.0 - ratio * ratio;
            return tmp * tmp;
        }
        return 0.0;
    }

    case kernel_type_t::S_TRIWEIGHT: {
        double ratio = distance / tau;
        if (ratio < 1.0) {
            double tmp = 1.0 - ratio * ratio;
            return tmp * tmp * tmp;
        }
        return 0.0;
    }

    default:
        return std::exp(-(distance * distance) / (tau * tau)); // Default to Gaussian
    }
}

/**
 * @brief Converts a Laplacian type enum to a string description
 *
 * @param type The Laplacian type enum value
 * @return std::string A string description of the Laplacian type
 */
std::string
set_wgraph_t::laplacian_type_to_string(laplacian_type_t type) const {
    switch (type) {
    case laplacian_type_t::STANDARD:
        return "Standard";
    case laplacian_type_t::NORMALIZED:
        return "Normalized";
    case laplacian_type_t::RANDOM_WALK:
        return "Random Walk";
    case laplacian_type_t::KERNEL:
        return "Kernel";
    case laplacian_type_t::NORMALIZED_KERNEL:
        return "Normalized Kernel";
    case laplacian_type_t::ADAPTIVE_KERNEL:
        return "Adaptive Kernel";
    case laplacian_type_t::SHIFTED:
        return "Shifted";
    case laplacian_type_t::SHIFTED_KERNEL:
        return "Shifted Kernel";
    case laplacian_type_t::REGULARIZED:
        return "Regularized";
    case laplacian_type_t::REGULARIZED_KERNEL:
        return "Regularized Kernel";
    case laplacian_type_t::MULTI_SCALE:
        return "Multi-Scale";
    default:
        return "Unknown";
    }
}
