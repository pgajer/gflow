#include "bandwidth_utils.hpp"
#include <algorithm>
#include <cmath>

#include "error_utils.h" // For REPORT_WARNING

/**
 * @brief Generates a set of candidate bandwidth values over the specified range
 *
 * @param min_bw Minimum bandwidth value (validated elsewhere)
 * @param max_bw Maximum bandwidth value (validated elsewhere)
 * @param n_bws Number of bandwidth values to generate
 * @param log_grid Whether to use logarithmic spacing (true) or linear spacing (false)
 * @param min_spacing Minimum meaningful difference between consecutive bandwidths
 *
 * @return Vector of candidate bandwidth values
 */
std::vector<double> get_candidate_bws(
    double min_bw,
    double max_bw,
    size_t n_bws,
    bool log_grid,
    double min_spacing) {

    // Handle the case where min_bw and max_bw are very close
    double bw_range = max_bw - min_bw;

    // If range is too small for requested number of bandwidths
    if (bw_range < (n_bws - 1) * min_spacing) {
        // Reduce the number of bandwidths based on the available range
        size_t effective_n_bws = std::max(size_t(2), size_t(bw_range / min_spacing) + 1);
        if (effective_n_bws < n_bws) {
            REPORT_WARNING("Warning: Bandwidth range (%.6f) is too small for %zu distinct values. "
                           "Reducing to %zu bandwidths.\n",
                           bw_range, n_bws, effective_n_bws);
            n_bws = effective_n_bws;
        }
    }

    // Generate candidate bandwidths using specified spacing strategy
    std::vector<double> bws;
    bws.reserve(n_bws);  // Reserve space but don't resize

    // Always include min_bw
    bws.push_back(min_bw);

    if (n_bws > 2) {
        // Generate intermediate bandwidths
        if (log_grid && min_bw > 0) {
            // Logarithmic spacing
            double log_min = std::log(min_bw);
            double log_max = std::log(max_bw);
            double log_step = (log_max - log_min) / (n_bws - 1);

            for (size_t i = 1; i < n_bws - 1; ++i) {
                double log_bw = log_min + i * log_step;
                bws.push_back(std::exp(log_bw));
            }
        } else {
            // Linear spacing
            double step = bw_range / (n_bws - 1);

            for (size_t i = 1; i < n_bws - 1; ++i) {
                bws.push_back(min_bw + i * step);
            }
        }
    }

    // Always include max_bw
    if (n_bws > 1) {
        bws.push_back(max_bw);
    }

    // Ensure we have the right number of bandwidths
    if (bws.size() != n_bws) {
        REPORT_WARNING("Warning: Generated %zu bandwidths instead of requested %zu\n",
                       bws.size(), n_bws);
    }

    return bws;
}
