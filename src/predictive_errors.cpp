#include "sampling.h"  // For C_runif_simplex()
#include "predictive_errors.hpp" // For bb_cri_t
#include "wasserstein_dist.h" // for C_wasserstein_distance_1D()

#include <vector>      // For std::vector usage throughout
#include <algorithm>   // For std::sort, std::fill
#include <numeric>     // For std::accumulate
#include <limits>      // For numeric limits

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Computes Bayesian credible intervals with configurable central location measure.
 *
 * @details For each vertex, computes:
 *          1. Central location (mean or median) of bootstrap estimates
 *          2. Two-sided Bayesian credible intervals from the bootstrap distribution
 *
 *          The central location can be either:
 *          - Mean (arithmetic average)
 *          - Median (middle value or average of two middle values)
 *
 * @param bb_Ey Bootstrap estimates matrix where:
 *              - Each row (outer vector) represents one bootstrap iteration
 *              - Each column (inner vector element) corresponds to a vertex
 * @param y Original observed values at each vertex (used for size validation)
 * @param use_median If true, use median as central location; if false, use mean (default: false)
 * @param p Probability level for credible intervals [0 < p < 1] (default: 0.95)
 *
 * @return bb_cri_t struct containing:
 *         - bb_Ey: Vector of central locations (mean/median) for each vertex
 *         - cri_L: Vector of lower credible interval bounds
 *         - cri_U: Vector of upper credible interval bounds
 *
 * @throws std::invalid_argument if:
 *         - bb_Ey is empty
 *         - Any bootstrap sample has different size than y
 *         - p is not in (0,1)
 */

/**
 * @brief Computes credible intervals from Bayesian bootstrap estimates
 *
 * @param bb_Ey Bootstrap estimates [n_bootstrap × n_vertices]
 * @param y Original response values [n_vertices]
 * @param use_median If true, uses median instead of mean for central location
 * @param p Probability level for credible intervals
 *
 * @return bb_cri_t struct containing central estimates and interval bounds
 */
bb_cri_t bb_cri(const std::vector<std::vector<double>>& bb_Ey,
                bool use_median = false,
                double p = 0.95) {
    // Input validation
    if (bb_Ey.empty()) {
        Rf_error("Bootstrap estimates cannot be empty");
    }
    if (p <= 0.0 || p >= 1.0) {
        Rf_error("Probability level p must be in (0,1)");
    }

    const size_t n_bb = bb_Ey.size();
    const size_t n_points = bb_Ey[0].size();

    // Validate dimensions
    for (const auto& sample : bb_Ey) {
        if (sample.size() != n_points) {
            Rf_error("Inconsistent dimensions in bootstrap samples");
        }
    }

    // Calculate indices for quantiles
    size_t lower_idx = static_cast<size_t>(std::floor((1.0 - p) / 2.0 * n_bb));
    size_t upper_idx = static_cast<size_t>(std::ceil((1.0 + p) / 2.0 * n_bb));

    // Initialize result structure
    bb_cri_t results;
    results.bb_Ey.resize(n_points);
    results.cri_L.resize(n_points);
    results.cri_U.resize(n_points);

    // Process each vertex
    std::vector<double> vertex_estimates(n_bb);
    for (size_t i = 0; i < n_points; ++i) {
        // Collect bootstrap estimates for current vertex
        for (size_t j = 0; j < n_bb; ++j) {
            vertex_estimates[j] = bb_Ey[j][i];
        }

        // Compute central location
        if (!use_median) {
            results.bb_Ey[i] = std::accumulate(vertex_estimates.begin(),
                                             vertex_estimates.end(), 0.0) / n_bb;
        } else {
            auto mid = vertex_estimates.begin() + n_bb/2;
            std::nth_element(vertex_estimates.begin(), mid, vertex_estimates.end());
            results.bb_Ey[i] = n_bb % 2 == 0 ?
                (*mid + *std::min_element(mid, vertex_estimates.end())) / 2.0 : *mid;
        }

        // Compute interval bounds
        std::vector<double> sorted_estimates = vertex_estimates;  // Make copy for sorting
        std::sort(sorted_estimates.begin(), sorted_estimates.end());
        results.cri_L[i] = sorted_estimates[lower_idx];
        results.cri_U[i] = sorted_estimates[upper_idx];
    }

    return results;
}

#if 0
bb_cri_t bb_cri(const std::vector<std::vector<double>>& bb_Ey,
                const std::vector<double>& y,
                bool use_median = false,
                double p = 0.95) {
    // Input validation
    if (bb_Ey.empty()) {
        Rf_error("bb_Ey cannot be empty");
    }
    if (p <= 0.0 || p >= 1.0) {
        Rprintf("In bb_cri()  p: %.3f\n", p);
        Rf_error("p must be in (0,1)");
    }

    const size_t n_points = y.size();
    const size_t n_bb = bb_Ey.size();

    // Validate dimensions
    for (const auto& sample : bb_Ey) {
        if (sample.size() != n_points) {
            Rf_error("All bootstrap samples must have same size as y");
        }
    }

    // Calculate indices for quantiles using safe arithmetic
    size_t lower_idx = static_cast<size_t>(std::floor((1.0 - p) / 2.0 * (n_bb - 1)));
    size_t upper_idx = static_cast<size_t>(std::ceil((1.0 + p) / 2.0 * (n_bb - 1)));

    // Initialize result vectors
    bb_cri_t results;
    results.cri_L.resize(n_points);
    results.cri_U.resize(n_points);
    results.bb_Ey.resize(n_points);

    // Helper function to compute location (mean or median) of a sequence
    auto compute_location = [use_median](const std::vector<double>& values) -> double {
        if (!use_median) {
            return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        } else {
            std::vector<double> sorted_values = values;
            std::sort(sorted_values.begin(), sorted_values.end());
            size_t n = sorted_values.size();
            if (n % 2 == 0) {
                return (sorted_values[n/2 - 1] + sorted_values[n/2]) / 2.0;
            } else {
                return sorted_values[n/2];
            }
        }
    };

    std::vector<double> vertex_estimates(n_bb);
    for (size_t i = 0; i < n_points; ++i) {
        // Collecting all bootstrap estimates for this vertex
        for (size_t j = 0; j < n_bb; ++j) {
            vertex_estimates[j] = bb_Ey[j][i];
        }

        // Computing the central location for this vertex
        results.bb_Ey[i] = compute_location(vertex_estimates);

        // Sort for quantiles (reuse the vector since we're done with location calculation)
        std::sort(vertex_estimates.begin(), vertex_estimates.end());

        // Store bounds
        results.cri_L[i] = vertex_estimates[lower_idx];
        results.cri_U[i] = vertex_estimates[upper_idx];
    }

    return results;
}
#endif

/**
 * @brief Computes Wasserstein distances between bootstrap estimates and observed values.
 *
 * For each vertex i, computes the Wasserstein distance between the empirical distribution
 * of its bootstrap estimates and the Dirac delta distribution at its observed value y[i].
 *
 * @param bb_Ey Bootstrap estimates where bb_Ey[j][i] is the estimate for vertex i
 *              in bootstrap sample j.
 * @param y Observed values where y[i] is the value at vertex i.
 *
 * @return Vector of Wasserstein distances, one per vertex.
 *
 * @throws std::invalid_argument If inputs are empty or dimensions mismatch.
 *
 * @pre bb_Ey must not be empty and all its inner vectors must have size equal to y.size().
 * @pre y must not be empty.
 *
 * @details
 * The function performs the following steps:
 * 1. Validates input dimensions and non-emptiness.
 * 2. For each vertex or point:
 *    a. Extracts bootstrap estimates for that vertex from all samples.
 *    b. Creates a vector of the observed value repeated for each bootstrap sample.
 *    c. Computes the Wasserstein distance between the bootstrap estimates
 *       and the repeated observed value using C_wasserstein_distance_1D.
 * 3. Returns the vector of computed Wasserstein distances.
 *
 * @note This function assumes the existence of a C_wasserstein_distance_1D function
 *       for computing 1D Wasserstein distances.
 *
 * @Rf_warning The function modifies the input vectors vertex_estimates and y_dirac
 *          in each iteration. Ensure these are not used externally if persistence
 *          is required.
 *
 * @pre bb_Ey must not be empty and all its inner vectors must have the same size as y.
 * @pre y must not be empty.
 *
 * @see C_wasserstein_distance_1D
 *
 * @example
 * std::vector<std::vector<double>> bb_Ey = {{1.0, 2.0, 3.0}, {1.5, 2.5, 3.5}};
 * std::vector<double> y = {1.2, 2.2, 3.2};
 * try {
 *     std::vector<double> distances = compute_bbwasserstein_errors(bb_Ey, y);
 *     // Use the computed distances...
 * } catch (const std::invalid_argument& e) {
 *     std::cerr << "Error: " << e.what() << std::endl;
 * }
 */
std::vector<double> compute_bbwasserstein_errors(const std::vector<std::vector<double>>& bb_Ey,
                                                 const std::vector<double>& y) {
    if (bb_Ey.empty() || y.empty()) {
        Rf_error("Input vectors must not be empty");
    }

    int n_points = static_cast<int>(y.size());
    int B = static_cast<int>(bb_Ey.size());

    if (bb_Ey[0].size() != y.size()) {
        Rf_error("Mismatch in dimensions between bb_Ey and y");
    }

    // Computing Bayesian bootstrap version of y
    std::vector<std::vector<double>> bb_y(B, std::vector<double>(n_points));
    std::vector<double> weights(n_points);
    for (int b = 0; b < B; ++b) {
        C_runif_simplex(&n_points, weights.data());
        for (int i = 0; i < n_points; ++i)
            bb_y[b][i] = n_points * weights[i] * y[i];
    }

    // Computing Wasserstein distances
    std::vector<double> dist(n_points);
    std::vector<double> vertex_estimates(B);
    std::vector<double> vertex_bb_y(B);

    for (int i = 0; i < n_points; ++i) {
        // Extract estimates and bb_y for vertex i
        for (int b = 0; b < B; ++b) {
            vertex_estimates[b] = bb_Ey[b][i];
            vertex_bb_y[b] = bb_y[b][i];
        }

        double wasserstein_dist;
        C_wasserstein_distance_1D(vertex_estimates.data(), vertex_bb_y.data(), &B, &wasserstein_dist);
        dist[i] = wasserstein_dist;
    }

    return dist;
}


/**
* @brief Computes Wasserstein distances between bootstrap estimates and bootstrap observations.
*
* For each vertex i, computes the Wasserstein distance between two empirical distributions:
* 1. The distribution of bootstrap estimates \f$\hat{\mu}_i^{(b)}\f$
* 2. The distribution of bootstrap observations \f$y_i^{*bb,b} = n w_i^b y_i\f$
* where \f$w_i^b\f$ are Bayesian bootstrap weights.
*
* @param bb_Ey Matrix of bootstrap estimates where bb_Ey[b][i] is the estimate
*              \f$\hat{\mu}_i^{(b)}\f$ for vertex i in bootstrap sample b
* @param y Vector of observed values where y[i] is the value at vertex i
* @param bb_y Matrix of bootstrap observations where bb_y[b][i] is
*              \f$y_i^{*bb,b} = n w_i^b y_i\f$ for vertex i in bootstrap sample b
*
* @return Vector of Wasserstein distances, one per vertex
*
* @throws std::invalid_argument If inputs are empty or dimensions mismatch
*
* @pre bb_Ey must not be empty and all its inner vectors must have size equal to y.size()
* @pre bb_y must have the same dimensions as bb_Ey
* @pre y must not be empty
*
* @note This function assumes bb_y has been computed using the same number of bootstrap
*       samples as bb_Ey and with properly generated Bayesian bootstrap weights.
*
* @see C_wasserstein_distance_1D for the underlying distance computation
* @see graph_kmean_predictive_error where bb_y is typically generated
*/
std::vector<double> compute_bbwasserstein_errors(const std::vector<std::vector<double>>& bb_Ey,
                                                 const std::vector<double>& y,
                                                 const std::vector<std::vector<double>>& bb_y) {
    if (bb_Ey.empty() || y.empty() || bb_y.empty()) {
        Rf_error("Input vectors must not be empty");
    }
    size_t n_points = static_cast<int>(y.size());
    int B = static_cast<int>(bb_Ey.size());

    if (bb_Ey[0].size() != n_points || bb_y.size() != B || bb_y[0].size() != n_points) {
        Rf_error("Mismatch in dimensions between inputs");
    }

    // Computing Wasserstein distances
    std::vector<double> dist(n_points);
    std::vector<double> vertex_estimates(B);
    std::vector<double> vertex_bb_y(B);

    for (size_t i = 0; i < n_points; ++i) {
        // Extract estimates and bb_y for vertex i
        for (int b = 0; b < B; ++b) {
            vertex_estimates[b] = bb_Ey[b][i];
            vertex_bb_y[b] = bb_y[b][i];
        }
        double wasserstein_dist;
        C_wasserstein_distance_1D(vertex_estimates.data(), vertex_bb_y.data(), &B, &wasserstein_dist);
        dist[i] = wasserstein_dist;
    }
    return dist;
}

/**
 * @brief Computes the Bayesian Bootstrap Mean Wasserstein Distance (BBMWD) Rf_error using Kahan summation.
 *
 * This function calculates the BBMWD between the observed values and their
 * Bayesian bootstrap estimates for each vertex in a graph. The BBMWD is defined as:
 *
 * \f[
 * \text{BBMWD}(E(y|G,d)|y) = \frac{1}{n}\sum_{i=1}^n d_W(\delta_{y_i}, \{p_i^w\}_{w \in W})
 * \f]
 *
 * where:
 * - \f$n\f$ is the number of vertices in the graph
 * - \f$y_i\f$ is the observed value at vertex \f$i\f$
 * - \f$\delta_{y_i}\f$ is the Dirac delta distribution at \f$y_i\f$
 * - \f$\{p_i^w\}_{w \in W}\f$ is the distribution of estimates at vertex \f$i\f$ from the Bayesian bootstrap
 * - \f$d_W\f$ is the Wasserstein distance
 *
 * The computation process is as follows:
 * 1. For each vertex \f$i\f$ in the graph:
 *    a. Extract the bootstrap estimates for vertex \f$i\f$ from all bootstrap samples.
 *    b. Create a Dirac delta distribution at the observed value \f$y_i\f$.
 *    c. Compute the Wasserstein distance between the bootstrap estimates and the Dirac delta distribution.
 * 2. Calculate the mean of these Wasserstein distances across all vertices.
 *
 * The function uses Kahan summation for improved numerical stability in the total distance calculation,
 * which is particularly important for large graphs or when dealing with widely varying Wasserstein distances.
 *
 * @param bb_Ey A vector of vectors containing Bayesian bootstrap estimates.
 *                     Each inner vector represents estimates for all vertices for one bootstrap sample.
 *                     bb_Ey[j][i] is the estimate for vertex i in bootstrap sample j.
 * @param y A vector of observed values for each vertex. y[i] is the observed value at vertex i.
 * @return The computed BBMWD Rf_error.
 *
 * @pre bb_Ey and y must not be empty.
 * @pre All inner vectors in bb_Ey must have the same size as y.
 *
 * @throws std::invalid_argument if preconditions are not met.
 *
 * @note The Wasserstein distance is computed using the C_wasserstein_distance_1D function,
 *       which should be implemented to calculate the 1D Wasserstein distance between two empirical distributions.
 *
 * @Rf_warning This function may be computationally intensive for large graphs or high numbers of bootstrap samples.
 *
 * @see C_wasserstein_distance_1D for the underlying Wasserstein distance computation.
 */
double compute_bbwasserstein_error(const std::vector<std::vector<double>>& bb_Ey,
                                   const std::vector<double>& y) {

    std::vector<double> distances = compute_bbwasserstein_errors(bb_Ey, y);

    double total_distance = 0.0;
    double c = 0.0;  // Kahan summation compensation
    double eps = 1e-10;
    int n_non_zero_dists = 0;

    // Compute mean using Kahan summation, counting non-zero distances
    for (double dist : distances) {
        if (dist > eps) {
            double y = dist - c;
            double t = total_distance + y;
            c = (t - total_distance) - y;
            total_distance = t;
            n_non_zero_dists++;
        }
    }

    return total_distance / n_non_zero_dists;
}


/**
 * @brief Computes predictive errors for graph kernel mean estimation using mean or median for location.
 *
 * @param bb_Ey Bootstrap estimates of conditional expectations
 * @param y Observed values
 * @param use_median If true, uses median for location estimation; if false, uses mean (default: false)
 * @return Pair of predictive errors (errorA, errorB)
 */
std::pair<double, double> compute_bbcov_error(const std::vector<std::vector<double>>& bb_Ey,
                                              const std::vector<double>& y,
                                              bool use_median = false) {

    if (bb_Ey.empty() || y.empty() || bb_Ey[0].size() != y.size()) {
        Rf_error("Invalid input dimensions");
    }
    int n_points = static_cast<int>(y.size());
    int B = static_cast<int>(bb_Ey.size());
    std::vector<double> weights(n_points);

    // Helper function to compute location (mean or median) of a sequence
    auto compute_location = [use_median](const std::vector<double>& values) -> double {
        if (!use_median) {
            return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        } else {
            std::vector<double> sorted_values = values;
            std::sort(sorted_values.begin(), sorted_values.end());
            size_t n = sorted_values.size();
            if (n % 2 == 0) {
                return (sorted_values[n/2 - 1] + sorted_values[n/2]) / 2.0;
            } else {
                return sorted_values[n/2];
            }
        }
    };

    // Compute location of bb_Ey
    std::vector<double> mu_location(n_points, 0.0);
    for (int i = 0; i < n_points; ++i) {
        std::vector<double> vertex_estimates(B);
        for (int b = 0; b < B; ++b) {
            vertex_estimates[b] = bb_Ey[b][i];
        }
        mu_location[i] = compute_location(vertex_estimates);
    }

    // Computing the apparent Rf_error
    double apparent_error = 0.0;
    for (int i = 0; i < n_points; ++i)
        apparent_error += std::pow(y[i] - mu_location[i], 2);

    // Computing Bayesian bootstrap version of y
    std::vector<std::vector<double>> bb_y(B, std::vector<double>(n_points));
    for (int b = 0; b < B; ++b) {
        C_runif_simplex(&n_points, weights.data());
        for (int i = 0; i < n_points; ++i)
            bb_y[b][i] = n_points * weights[i] * y[i];
    }

    // Computing location of bb_y
    std::vector<double> bb_y_location(n_points, 0.0);
    for (int i = 0; i < n_points; ++i) {
        std::vector<double> vertex_bb_y(B);
        for (int b = 0; b < B; ++b) {
            vertex_bb_y[b] = bb_y[b][i];
        }
        bb_y_location[i] = compute_location(vertex_bb_y);
    }

    // Computing both covariance penalties
    double covariance_penaltyA = 0.0;
    double covariance_penaltyB = 0.0;

    for (int i = 0; i < n_points; ++i) {
        double cov_i_A = 0.0;
        double cov_i_B = 0.0;
        for (int b = 0; b < B; ++b) {
            double y_centered = bb_y[b][i] - bb_y_location[i];
            cov_i_A += bb_Ey[b][i] * y_centered;
            cov_i_B += (bb_Ey[b][i] - mu_location[i]) * y_centered;
        }
        covariance_penaltyA += cov_i_A / (B - 1);
        covariance_penaltyB += cov_i_B / (B - 1);
    }

    return std::make_pair(
        apparent_error + 2 * covariance_penaltyA,
        apparent_error + 2 * covariance_penaltyB
    );
}
