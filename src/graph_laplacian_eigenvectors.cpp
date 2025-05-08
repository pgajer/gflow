// Standard C++ library headers
#include <vector>
#include <unordered_map>
#include <cmath>        // For std::exp
#include <algorithm>    // For std::min
#include <cstdio>       // For Rprintf

// Eigen headers
#include <Eigen/Core>
#include <Eigen/Dense>  // For Eigen::MatrixXd
#include <Eigen/Sparse> // For Eigen::SparseMatrix, Triplet

// Spectra headers
#include <SymEigsSolver.h>          // For SymEigsSolver
#include <MatOp/SparseSymMatProd.h> // For SparseSymMatProd

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
	for (int i = 0; i < n_vertices; i++) {
		L.coeffRef(i, i) += lap_reg;  // Small regularization
	}

	Spectra::SparseSymMatProd<double> op(L);

	int nev = std::min<int>(2 * n_evectors + 5, n_vertices);
	int ncv = std::min<int>(2 * nev, n_vertices);  // Control parameter for the algorithm

	// Ensure nev < ncv
	if (nev >= ncv) {
		ncv = nev + 1;
	}

	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
	eigs.init();
	int maxit = 1000;
	double tol = 1e-10;
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
		// If all attempts failed, report an error
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
    for (int i = 0; i < K.outerSize(); ++i) {
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
    for (int i = 0; i < n_vertices; ++i) {
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

    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    int maxit = 1000;
    double tol = 1e-10;
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
    for (int i = 0; i < K.outerSize(); ++i) {
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
    for (int i = 0; i < n_vertices; ++i) {
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

    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    int maxit = 1000;
    double tol = 1e-10;
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
