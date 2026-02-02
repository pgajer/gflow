#include "set_wgraph.hpp"              // set_wgraph_t, adjacency_list
#include "klaps_low_pass_smoother.hpp" // klaps_low_pass_smoother_t declaration
#include "bandwidth_utils.hpp"         // get_candidate_ks()
#include "error_utils.h"               // REPORT_ERROR()

#include <vector>
#include <utility>                 // std::pair
#include <algorithm>               // std::min, std::accumulate
#include <numeric>                 // std::accumulate
#include <cmath>                   // std::log2, std::exp
#include <limits>                  // for numeric_limits

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <R.h>                     // Rprintf

std::vector<size_t> get_candidate_ks(
	size_t k_min,
	size_t k_max,
	size_t n_candidates,
	bool   log_grid
	);
/**
 * @brief Smooths a vertex‐valued signal by projecting onto low‐frequency Laplacian modes.
 *
 * Implements a family of low‐pass filters on the graph by retaining only the first
 * \(k\) eigencomponents of the Laplacian spectrum, for a range of candidate \(k\)'s.
 * Three criteria are computed to choose the optimal \(k\):
 *   - **Eigengap**: maximize \(\lambda_{k+1} - \lambda_k\),
 *   - **Generalized Cross‐Validation (GCV)**:
 *       \(\mathrm{GCV}(k) = \frac{\|y - \hat y^{(k)}\|^2}{(n - k)^2}\),
 *   - **Spectral‐Energy Threshold**: minimal \(k\) with cumulative energy ≥ threshold.
 *
 * Optionally returns the full family of reconstructed signals for each \(k\)
 * if `with_k_predictions` is set to true.
 *
 * @param y
 *   A length‑n vector of observed values at each graph vertex.
 *
 * @param n_evectors_to_compute
 *   Number of Laplacian eigenvectors to precompute (should exceed `max_num_eigenvectors`).
 *
 * @param min_num_eigenvectors
 *   Smallest \(k\) to evaluate (must be ≥ 1).
 *
 * @param max_num_eigenvectors
 *   Largest \(k\) to evaluate (≤ `n_evectors_to_compute`).
 *
 * @param tau_factor
 *     Positive real scaling factor that determines the kernel bandwidth \(\tau\) as a fraction of the graph diameter.
 *     The kernel bandwidth is computed as:
 *         \[
 *         \tau = \text{tau\_factor} \times \text{graph diameter}
 *         \]
 *     Smaller `tau_factor` values result in more localized smoothing; larger values make smoothing more global.
 *
 * @param radius_factor A real scaling factor of tau radius that is not less than 1.
 *
 * @param laplacian_power
 *     Positive odd integer specifying the power to which (I - L) is raised.
 *     Higher values apply stronger smoothing by repeatedly reinforcing low-pass filtering.
 *     Typically, larger `laplacian_power` requires smaller `tau_factor` to maintain the smoothing scale.
 *
 * @param n_candidates
 *   Number of \(k\)-values to test between `min_num_eigenvectors` and `max_num_eigenvectors`.
 *   If `log_grid` is false, they are equally spaced; otherwise log‑spaced.
 *
 * @param log_grid
 *   Whether to place candidate \(k\)'s on a logarithmic scale.
 *
 * @param energy_threshold
 *   Fraction of total spectral energy (sum of eigenvalues) at which to stop for the
 *   energy‐threshold criterion (e.g. 0.90 for 90%).
 *
 * @param with_k_predictions
 *   If true, returns `k_predictions`, a length‑`n_candidates` array of full
 *   reconstructed signals for each tested \(k\).  Otherwise, only the final
 *   chosen `predictions` vector is stored.
 *
 * @param verbose
 *   If true, prints progress messages (via Rprintf) during:
 *     - Eigen‐decomposition,
 *     - Looping over candidate \(k\)'s, and
 *     - Final criterion selection.
 *
 * @return
 *   A fully populated `klaps_low_pass_smoother_t` containing:
 *     - `evalues`, `evectors`,
 *     - `candidate_ks`,
 *     - criterion vectors: `eigengaps`, `gcv_scores`, `spectral_energy`,
 *     - optimal indices: `opt_k_eigengap`, `opt_k_gcv`, `opt_k_spectral_energy`,
 *     - `used_method` indicator,
 *     - `predictions`: the chosen low‑pass smooth,
 *     - optionally `k_predictions`: all intermediate smooths.
 *
 * @see klaps_low_pass_smoother_t for detailed field descriptions.
 */
klaps_low_pass_smoother_t
set_wgraph_t::klaps_low_pass_smoother(
	const std::vector<double>& y,
	size_t n_evectors_to_compute,
	size_t min_num_eigenvectors,
	size_t max_num_eigenvectors,
	double tau_factor,
	double radius_factor,
	size_t laplacian_power,
	size_t n_candidates,
	bool   log_grid,
	double energy_threshold,
	bool   with_k_predictions,
	bool   verbose
	) const {

	size_t n_vertices = adjacency_list.size();
	klaps_low_pass_smoother_t result;

	// 1) Compute Laplacian eigen‐decomposition
	double tau = tau_factor * graph_diameter; // (std::sqrt(k) * 0.95 * n_evectors_to_compute);
	// spectrum.first  is an Eigen::VectorXd of length m
	// spectrum.second is an Eigen::MatrixXd of size (n_vertices × m)
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> spectrum = compute_graph_shifted_kernel_laplacian_spectrum(
		n_evectors_to_compute,
		tau,
		radius_factor,
		laplacian_power,
		verbose
		);

	// Copy eigenvalues into std::vector<double>
	const Eigen::VectorXd& ev = spectrum.first;
	result.evalues.assign(ev.data(), ev.data() + ev.size());

	// Copy eigenvectors into the Eigen::MatrixXd
	result.evectors = std::move(spectrum.second);

	// 2) Build candidate_ks between [min, max]
	result.candidate_ks = get_candidate_ks(
		min_num_eigenvectors,
		max_num_eigenvectors,
		n_candidates,
		log_grid
		);

	// 3) Precompute Graph Fourier Transform of y
	Eigen::VectorXd y_ev = Eigen::Map<const Eigen::VectorXd>(y.data(), n_vertices);
	Eigen::VectorXd gft   = result.evectors.transpose() * y_ev;

	// 4) For each k in candidate_ks, compute:
	//- low-pass y_k = ∑_{i<k} gft[i] * result.evectors.col(i)
	//- residual r_k = y_ev − y_k
	//- squared Rf_error ||r_k||²
	std::vector<Eigen::VectorXd> low_passes;
	low_passes.reserve(result.candidate_ks.size());
	std::vector<double> sq_errors, trace_S;

	//int nev = result.evectors.size();
	for (size_t idx = 0; idx < result.candidate_ks.size(); ++idx) {
		size_t k = result.candidate_ks[idx];
		// size_t m = result.evectors.cols();
		Eigen::VectorXd yk = Eigen::VectorXd::Zero(n_vertices);
		for (size_t j = 0; j < k; ++j) {
			//size_t idx = m - 1 - j;                     // starts at nev‑1, then nev‑2, …
			size_t idx = j;
			yk += gft[idx] * result.evectors.col(idx);
		}

		low_passes.push_back(yk);

		Eigen::VectorXd resid = y_ev - yk;
		sq_errors.push_back(resid.squaredNorm());
		trace_S.push_back(static_cast<double>(k));   // Smoothing matrix trace ≈ k
	}

	// 5) Eigengap: λ_{i+1} − λ_i
	{
		auto& gaps = result.eigengaps;
		gaps.resize(result.evalues.size() - 1);
		for (size_t i = 0; i + 1 < result.evalues.size(); ++i)
			gaps[i] = result.evalues[i+1] - result.evalues[i];

		// Find the index in candidate_ks whose gap is largest
		double max_gap = -std::numeric_limits<double>::infinity();
		size_t best_idx = 0;
		for (size_t idx = 0; idx < result.candidate_ks.size(); ++idx) {
			size_t k = result.candidate_ks[idx];
			// gap between lambda_{i+1} and lambda_{i} is stored at gaps[k-1]
			double g = gaps[k - 1];
			if (g > max_gap) {
				max_gap = g;
				best_idx = idx;
			}
		}

		result.opt_k_eigengap = best_idx;
	}

	// 6) GCV: Rf_error(k) = sq_errors[k] / (n_vertices − k)²
	{
		auto& scores = result.gcv_scores;
		scores.resize(low_passes.size());
		for (size_t i = 0; i < scores.size(); ++i) {
			double denom = n_vertices - trace_S[i];
			scores[i] = sq_errors[i] / (denom*denom);
		}

		// Pick the index of the minimum GCV score
		auto it = std::min_element(result.gcv_scores.begin(), result.gcv_scores.end());
		result.opt_k_gcv = std::distance(result.gcv_scores.begin(), it);
	}

	// 7) Spectral‐energy: cumulative ∑_{i<k} λ_i / ∑ λ_i
	{
		double total = std::accumulate(result.evalues.begin(), result.evalues.end(), 0.0);
		auto& energy = result.spectral_energy;
		energy.resize(result.candidate_ks.size());
		for (size_t i = 0; i < energy.size(); ++i) {
			size_t k = result.candidate_ks[i];
			double cum = std::accumulate(result.evalues.begin(),
										 result.evalues.begin()+k, 0.0);
			energy[i] = cum / total;
		}

		// First index where energy ≥ energy_threshold
		size_t best_idx = result.spectral_energy.size() - 1;
		for (size_t i = 0; i < result.spectral_energy.size(); ++i) {
			if (result.spectral_energy[i] >= energy_threshold) {
				best_idx = i;
				break;
			}
		}
		result.opt_k_spectral_energy = best_idx;// choose first k with energy ≥ energy_threshold
	}

	// 8) Build final predictions (default: GCV)
	{
		// default to the GCV‐chosen candidate
		size_t choice_idx = result.opt_k_gcv;
		result.predictions.assign(
			low_passes[choice_idx].data(),
			low_passes[choice_idx].data() + n_vertices
			);
		result.used_method = klaps_low_pass_smoother_t::method_t::GCV;
	}

	// Optionally, if the user requested all intermediate smooths:
	if (with_k_predictions) {
		result.k_predictions.resize(low_passes.size());
		for (size_t i = 0; i < low_passes.size(); ++i) {
			result.k_predictions[i].assign(
				low_passes[i].data(),
				low_passes[i].data() + n_vertices
				);
		}
	}

	auto is_binary01 = [](const std::vector<double>& yy, double tol = 1e-12) -> bool {
        for (double v : yy) {
            if (!(std::fabs(v) <= tol || std::fabs(v - 1.0) <= tol)) {
                return false;
            }
        }
        return true;
    };

    const bool y_binary = is_binary01(y);

	if (y_binary) {
		for (size_t i = 0; i < result.predictions.size(); ++i) {
			result.predictions[i] = std::clamp(result.predictions[i], 0.0, 1.0);
		}

		if (with_k_predictions) {
			for (size_t j = 0; j < low_passes.size(); ++j) {
				for (size_t i = 0; i < result.predictions.size(); ++i) {
					result.k_predictions[j][i] = std::clamp(result.k_predictions[j][i], 0.0, 1.0);
				}
			}
		}
	}

	return result;
}

/**
 * @brief Computes the partial eigenspectrum (eigenvalues & eigenvectors) of the graph Laplacian.
 *
 * This routine constructs the combinatorial Laplacian matrix \(L = D - A\) from the
 * adjacency_list of the current graph, applies a small diagonal regularization
 * to ensure positive definiteness, and then invokes the Spectra library to compute
 * the smallest eigenvalues and corresponding eigenvectors.
 *
 * The user can request a limited number of eigenpairs (`n_evectors`), and
 * Spectra’s internal parameters (`nev` and `ncv`) are chosen adaptively to ensure
 * convergence.  If the initial solve fails, several fallback strategies are attempted:
 * increasing the maximum number of iterations, relaxing the tolerance, and finally
 * increasing the Krylov subspace dimension (`ncv`).
 *
 * @note All sparse matrices (A, D, and L) are compressed before decomposition
 *   to satisfy Spectra’s CSR input requirements.
 *
 * @param n_evectors
 *   The number of nontrivial eigenpairs to compute.  Internally this will
 *   request up to `2*n_evectors + 5` eigenvalues to improve robustness.
 *
 * @param verbose
 *   If true, prints diagnostic messages via `Rprintf()` when:
 * - The fallback strategies succeed with adjusted parameters, and
 * - The final number of eigenpairs actually computed.
 *
 * @return
 *   A `std::pair` consisting of:
 *   - `first`  : an Eigen::VectorXd of eigenvalues (length = actual # computed), sorted ascending,
 *   - `second` : an Eigen::MatrixXd whose columns are the matching eigenvectors.
 *
 * @throws std::runtime_error (via REPORT_ERROR) if all fallback strategies fail.
 */
std::pair<Eigen::VectorXd, Eigen::MatrixXd>
set_wgraph_t::compute_graph_laplacian_spectrum(
	size_t n_evectors,
	bool   verbose
	) const {

	size_t n_vertices = adjacency_list.size();

	Rprintf("\n\nEntering compute_graph_laplacian_spectrum(): n_vertices: %zu\n", n_vertices);

	// 1) Build sparse A, degree D, Laplacian L = D − A
	Eigen::SparseMatrix<double> A(n_vertices,n_vertices), D(n_vertices,n_vertices);
	std::vector<Eigen::Triplet<double>> T;
	T.reserve(n_vertices * 4);
	for (size_t i = 0; i < n_vertices; ++i) {
		for (auto & edge : adjacency_list[i]) {
			size_t j = edge.vertex;
			if (i < j) {
				T.emplace_back(i,j, edge.weight);
				T.emplace_back(j,i, edge.weight);
			}
		}
	}
	A.setFromTriplets(T.begin(), T.end());

	for (int k = 0; k < A.outerSize(); ++k) {
		double sum = 0;
		for (auto it = Eigen::SparseMatrix<double>::InnerIterator(A,k); it; ++it)
			sum += it.value();
		D.insert(k,k) = sum;
	}

	Eigen::SparseMatrix<double> L = D - A;
	for (size_t i = 0; i < n_vertices; ++i)
		L.coeffRef(i,i) += 1e-8;   // regularization

	A.makeCompressed();
	D.makeCompressed();
	L.makeCompressed();

	// 2) Spectra setup
	Spectra::SparseSymMatProd<double> op(L);

	// Spectra requires  1 <= nev < ncv <= n_vertices

	// (A) force nev ≤ n_vertices-2
	int max_nev = (n_vertices > 2 ? n_vertices - 2 : 1);
	// request enough for convergence, but never above max_nev
	int nev     = std::min<int>(static_cast<int>(2 * n_evectors + 5), max_nev);
	nev = std::max<int>(1, nev);   // always at least one

	// (B)) now build ncv in (nev, n_vertices-1]
	int max_ncv = (n_vertices > 1 ? n_vertices - 1 : 1);
	int ncv     = std::min<int>(2 * nev, max_ncv);
	// if that didn't give us ncv>nev, bump it up by exactly one
	if (ncv <= nev) ncv = std::min<int>(nev + 1, max_ncv);

	if (verbose) {
		Rprintf("Spectra: nev=%d, ncv=%d, n_vertices=%d\n",
				nev, ncv, (int)n_vertices);
	}

	if (nev >= ncv) {
		REPORT_ERROR("nev: %d  ncv: %d - Spectra requires nev < ncv\n", nev, ncv);
	}

	if (ncv > (int)n_vertices) {
		REPORT_ERROR("ncv: %d n_vertices: %d - Spectra requires nev < ncv\n", ncv, (int)n_vertices);
	}

	// finally safe to ask Spectra for (nev,ncv)

	// ncv rule-of-thumb: 3*nev is often safer for hard problems
	int ncv_default = std::max(2 * nev + 10, 150);  // for nev=50 => at least 150
	ncv_default = std::min(ncv_default, (int)L.rows()); // cannot exceed n

	Spectra::SymEigsSolver<decltype(op)> eigs(op, nev, ncv_default);
	eigs.init();
	int maxit = 1000;
	double tol = 1e-10;
	eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

	Eigen::VectorXd evalues;
	Eigen::MatrixXd evectors;

	if (eigs.info() != Spectra::CompInfo::Successful) {
		// Define fallback parameters to try
		std::vector<std::pair<int, double>> attempts = {
			{2000, 1e-8},// More iterations, slightly relaxed tolerance
			{3000, 1e-6},// Even more iterations, more relaxed tolerance
			{5000, 1e-4} // Final attempt with very relaxed parameters
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

				evalues  = eigs.eigenvalues();
				evectors = eigs.eigenvectors();

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

						evalues  = adjusted_eigs.eigenvalues();
						evectors = adjusted_eigs.eigenvectors();

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
		evalues  = eigs.eigenvalues();
		evectors = eigs.eigenvectors();
	}

	if (verbose) {
		Rprintf("Computed %d eigenpairs (first %zu requested)\n",
				(int)evalues.size(), n_evectors);
	}

	return { evalues, evectors };
}

/**
 * @brief Generate a list of candidate k values between k_min and k_max.
 *
 * If log_grid is false, values are evenly spaced in integer steps;
 * if true, they are spaced uniformly on the log scale of k.
 * Duplicate values (due to rounding) are removed, and the list is guaranteed
 * to include k_min and k_max.
 *
 * @param k_min          Minimum k (must be >= 1)
 * @param k_max          Maximum k (must be >= k_min)
 * @param n_candidates  Number of k values to generate (>= 1)
 * @param log_grid       Whether to space on a logarithmic scale
 * @return               Sorted, unique vector of candidate k's
 * @throws std::invalid_argument if n_candidates == 0 or k_min < 1
 */
inline std::vector<size_t> get_candidate_ks(
	size_t k_min,
	size_t k_max,
	size_t n_candidates,
	bool   log_grid
	) {

	// Trivial case: single candidate
	if (n_candidates == 1 || k_min == k_max) {
		return { k_min };
	}

	std::vector<size_t> ks;
	ks.reserve(n_candidates);

	if (!log_grid) {
		// Linear spacing
		double span = static_cast<double>(k_max - k_min);
		for (size_t i = 0; i < n_candidates; ++i) {
			double t = static_cast<double>(i) / static_cast<double>(n_candidates - 1);
			size_t k = k_min + static_cast<size_t>(std::round(t * span));
			ks.push_back(k);
		}
	} else {
		// Logarithmic spacing
		double log_min = std::log(static_cast<double>(k_min));
		double log_max = std::log(static_cast<double>(k_max));
		double span= log_max - log_min;
		for (size_t i = 0; i < n_candidates; ++i) {
			double t   = static_cast<double>(i) / static_cast<double>(n_candidates - 1);
			double val = std::exp(log_min + t * span);
			size_t k   = static_cast<size_t>(std::round(val));
			ks.push_back(k);
		}
	}

	// Ensure boundaries
	ks.front() = k_min;
	ks.back()  = k_max;

	// Sort and dedupe
	std::sort(ks.begin(), ks.end());
	ks.erase(std::unique(ks.begin(), ks.end()), ks.end());

	return ks;
}
