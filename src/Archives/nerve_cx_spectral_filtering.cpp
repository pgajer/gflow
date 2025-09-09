#include "nerve_cx_spectral_filtering.h"

/**
 * @brief Helper function to compute eigendecomposition of a sparse matrix
 *
 * @param L Sparse matrix to decompose
 * @param n_evectors Number of eigenvectors to compute
 * @param smallest_first Whether to compute smallest eigenvalues first (true) or largest (false)
 * @param verbose Whether to print progress information
 * @return std::pair<Eigen::VectorXd, Eigen::MatrixXd> Eigenvalues and eigenvectors
 */


//////
std::pair<Eigen::VectorXd, Eigen::MatrixXd>
compute_matrix_spectrum(
    const Eigen::SparseMatrix<double>& L,
    size_t n_evectors,
    bool smallest_first,
    bool verbose
) {
    size_t n_vertices = L.rows();

    // Compute eigendecomposition
    Spectra::SparseSymMatProd<double> op(L);

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
    if (ncv > n_vertices) {
        ncv = n_vertices;
        nev = std::max(1, static_cast<int>(ncv - 1));
    }

    if (verbose) {
        Rprintf("Eigendecomposition parameters: nev=%d, ncv=%d, n_vertices=%zu\n",
                nev, ncv, n_vertices);
    }

    // Create solver and compute
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>>
        eigs(op, nev, ncv);

    eigs.init();
    int maxit = 1000;
    double tol = 1e-10;

    // Sort based on whether we want smallest or largest eigenvalues
    Spectra::SortRule sort_rule = smallest_first ?
        Spectra::SortRule::SmallestAlge :
        Spectra::SortRule::LargestAlge;

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

    // // If power > 1, apply power to eigenvalues
    // if (power > 1) {
    //     for (int i = 0; i < eigenvalues.size(); ++i) {
    //         eigenvalues(i) = std::pow(eigenvalues(i), power);
    //     }
    // }

    if (verbose) {
        Rprintf("Computed %d eigenpairs (requested %zu)\n",
                (int)eigenvalues.size(), n_evectors);
    }

    return {eigenvalues, eigenvectors};
}
