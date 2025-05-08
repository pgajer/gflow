#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

#include "kernels.h"      // for kernel_fn
#include "error_utils.h"  // for REPORT_ERROR()

/**
 * @brief Compute normalized kernel weights from a distance vector
 *
 * Given a vector of distances from a reference point, this function normalizes
 * the distances and evaluates a pre-initialized kernel function to compute weights.
 * The weights are then normalized to sum to 1. This function is typically used in
 * kernel smoothing or local regression settings.
 *
 * @details
 * The kernel function must be initialized in advance using:
 * \code
 *     initialize_kernel(kernel_type, 1.0);
 * \endcode
 * The distances are scaled by the largest observed distance multiplied by
 * a normalization factor. After computing raw kernel weights, the function
 * normalizes them to sum to 1. If all weights are zero, they are renormalized
 * using a small epsilon threshold to prevent division by zero.
 *
 * @param dists Vector of distances from a reference point (modified in-place)
 * @param dist_normalization_factor Multiplicative scaling factor for normalization
 *
 * @return A vector of normalized kernel weights corresponding to each input distance
 *
 * @throws REPORT_ERROR if all distances are zero (i.e., max_dist == 0)
 *
 * @note The input distance vector is modified in-place during normalization.
 * @note This function assumes that the global `kernel_fn` pointer has been initialized.
 * @see initialize_kernel
 */
std::vector<double> get_weights(
    std::vector<double>& dists,
    double dist_normalization_factor
    ) {

    std::vector<double> weights;
    size_t n_dists = dists.size();
    if (n_dists == 0) {
        return weights;
    } else if (n_dists == 1) {
        weights.push_back(1.0);
        return weights;
    }

    weights.resize(n_dists);
    double max_dist = *std::max_element(dists.begin(), dists.end());

    if (max_dist == 0) {
        REPORT_ERROR("ERROR: max_dist = 0. At least one distance has to be greater than 0");
    }

    max_dist *= dist_normalization_factor;

    // Normalize by the scaled maximum
    for (size_t i = 0; i < n_dists; ++i) {
        dists[i] /= max_dist;
    }

    // Use the kernel function to calculate weights
    kernel_fn(dists.data(), static_cast<int>(n_dists), weights.data());

    // Normalize weights with protection against zero sum
    double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum_weights < std::numeric_limits<double>::epsilon()) {
        sum_weights = std::numeric_limits<double>::epsilon();
    }
    for (auto& w : weights) {
        w /= sum_weights;
    }

    return weights;
}
