#ifndef KLAPS_SPECTRAL_FILTER_H_
#define KLAPS_SPECTRAL_FILTER_H_

#include <vector>
#include <Eigen/Core>

// Define the enum for filter types
enum class filter_type_t {
    HEAT,           // Standard heat kernel filter: exp(-t*λ)
    GAUSSIAN,       // Gaussian spectral filter: exp(-t*λ²)
    NON_NEGATIVE    // Non-negative truncated filter: exp(-t*max(λ,0))
};

/**
 * @brief Results of a graph heat‐kernel (spectral diffusion) smoothing.
 */
struct klaps_spectral_filter_t {
	std::vector<double>   evalues;        ///< Laplacian eigenvalues (ascending)
	Eigen::MatrixXd       evectors;       ///< Laplacian eigenvectors (columns)
	std::vector<double>   candidate_ts;   ///< Diffusion times tested
	std::vector<double>   gcv_scores;     ///< GCV score for each t
	size_t                opt_t_idx;      ///< Index of optimal t in candidate_ts
	std::vector<double>   predictions;    ///< y_{t*} (optimal smooth)
	std::vector<std::vector<double>> t_predictions;  ///< Optional: full family y_{t_j}, size = candidate_ts.size()
};

#endif // KLAPS_SPECTRAL_FILTER_H_
