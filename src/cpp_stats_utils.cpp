#include "error_utils.h" // For REPORT_ERROR()

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <memory>
#include <set>
#include <cmath>         // For std::exp

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_ecdf(SEXP x);
}


/*!
   Calculates the empirical cumulative distribution function (ECDF)

   This function calculates the ECDF of a given vector of data points. The ECDF
   at a value x is defined as the proportion of data points in the vector that
   are less than or equal to x.

   \param x A vector of numeric data points.

   \return A unique pointer to a vector containing the ECDF values
           corresponding to the elements in `x`.
*/
std::unique_ptr<std::vector<double>> ecdf(const std::vector<double> &x) {

#define DEBUG_ecdf 0

    double n = x.size();
    std::map<double, int> mult_map;
    for (int i = 0; i < n; i++)
        mult_map[x[i]] = 0;

    for (int i = 0; i < n; i++)
        mult_map[x[i]] += 1;

#if DEBUG_ecdf
    Rprintf("mult_map:\n");
    for (const auto& pair : mult_map)
       Rprintf("(%f, %d)\n", pair.first, pair.second);
    Rprintf("\n");
#endif

    // Creating a map to store the ECDF values for each unique value in x
    std::unordered_map<double, double> ecdf_map;

    // Calculating the empirical cumulative distribution function (ECDF)
    double cumsum = 0;
    for (const auto& pair : mult_map) { // this effectively loops over the unique values of x
        double value = pair.first;
        cumsum += pair.second;
        ecdf_map[value] = cumsum / n;
    }

    // Creating the result vector and populating it with ECDF values in the order of the original x
    auto result = std::make_unique<std::vector<double>>(x.size());
    for (int i = 0; i < n; i++) {
        result->at(i) = ecdf_map[x[i]];
    }

    return result;
}

/*!
  Calculates the empirical cumulative distribution function (ECDF) for a given vector of data.

  The ECDF is a step function that represents the cumulative distribution of a sample of data.
  It is defined as the proportion of observations in the sample that are less than or equal to
  a given value. This function calculates the ECDF for each unique value in the input vector.

  The function takes an SEXP (S expression) representing a numeric vector as input and returns
  an SEXP representing the ECDF values. The ECDF values are computed using the C++ `ecdf` function.

  \param x An SEXP representing a numeric vector of data.

  \return An SEXP representing the ECDF values for each unique value in the input vector.
  The return value is a numeric vector with the same length as the number of unique values
  in the input vector. The names attribute of the return value is set to "ecdf".

  \note The input vector is sorted internally to compute the ECDF. If the input vector contains
  missing values (NA), they will be removed before computing the ECDF.

  Example usage in R:
  \code{
  # Assuming the shared library containing S_ecdf is loaded
  x <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  result <- S_ecdf(x)
  print(result)
  }

  Output:
  \code{
  ecdf
  0.2 0.4 0.6 0.8 1.0
  }
*/
SEXP S_ecdf(SEXP x) {

    // Create a vector from the data
    PROTECT(x = Rf_coerceVector(x, REALSXP));
    std::vector<double> vec_x(REAL(x), REAL(x) + (size_t)XLENGTH(x));
    UNPROTECT(1); // x

    // Call the C++ ecdf function
    std::unique_ptr<std::vector<double>> result = ecdf(vec_x);

    // Convert the result to an SEXP
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, result->size()));
    std::copy(result->begin(), result->end(), REAL(ans));

    UNPROTECT(1);
    return ans;
}


/**
 * @brief Linear transformation of 'x' values so that the min(x) is sent to ymin and max(x) to ymax
 *
 * This function scales the input vector `x` such that its minimum value is mapped to `ymin`
 * and its maximum value is mapped to `ymax`. The scaling is done in place, modifying the input vector.
 *
 * @param x Input vector of double values to be scaled.
 * @param ymin The value to which the minimum of `x` is to be mapped.
 * @param ymax The value to which the maximum of `x` is to be mapped.
 *
 * @Rf_error if `x` is empty or if `ymin` is equal to `ymax`.
 */
void scale_to_range(std::vector<double>& x, double ymin, double ymax) {
    if (x.empty()) {
        Rf_error("Input vector x is empty.");
    }

    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());

    if (xmin == xmax) {
        std::fill(x.begin(), x.end(), (ymin + ymax) / 2);
        return;
    }

    double slope = (ymax - ymin) / (xmax - xmin);
    for (auto& xi : x) {
        xi = slope * (xi - xmin) + ymin;
    }
}


/**
 * @brief Performs a running window average smoothing on a sequence using sliding window
 *
 * @details This function computes the average of values within a window centered
 *          at each position, using an efficient sliding window approach that reuses
 *          previous calculations. For each new position, it removes the leftmost value
 *          and adds the rightmost value of the new window.
 *
 * @param[in] values The input sequence to be smoothed
 * @param[in] window_size The number of neighbors to include on each side
 *
 * @return A vector of smoothed values corresponding to the input sequence
 *
 * @note A window_size of 0 returns the original sequence unchanged
 * @note Time complexity: O(n), where n is the length of the input sequence
 *
 * @example
 * // Example usage:
 * std::vector<double> data = {1.0, 3.0, 7.0, 2.0, 9.0};
 * std::vector<double> smoothed = running_window_average(data, 1);
 * // smoothed = {2.0, 3.67, 4.0, 6.0, 5.5}
 */
std::vector<double> running_window_average(const std::vector<double>& values, int window_size) {
    size_t n = values.size();
    std::vector<double> smoothed(n, 0.0);

    if (n == 0) return smoothed;

    // Calculate the initial window sum
    double window_sum = 0.0;
    size_t left = 0;
    size_t right = std::min(window_size, static_cast<int>(n) - 1);

    for (size_t j = left; j <= right; ++j) {
        window_sum += values[j];
    }

    // First window average
    size_t window_count = right - left + 1;
    smoothed[0] = window_sum / window_count;

    // Slide the window for each position
    for (size_t i = 1; i < n; ++i) {
        // Update window boundaries
        size_t prev_left = (i - 1 >= window_size) ? i - 1 - window_size : 0;
        size_t new_right = std::min(i + window_size, n - 1);

        // Update window sum by removing leftmost value (if it's no longer in window)
        // and adding rightmost value (if it's newly in window)
        if (i > window_size) {
            window_sum -= values[prev_left];
        }

        if (new_right > right) {
            window_sum += values[new_right];
            right = new_right;
        }

        // Compute window size (handles edge cases)
        left = (i >= window_size) ? i - window_size : 0;
        window_count = right - left + 1;

        // Compute and store the average
        smoothed[i] = window_sum / window_count;
    }

    return smoothed;
}

/**
 * @brief Calculates the Jaccard similarity index between two sets
 *
 * The Jaccard index measures the similarity between finite sample sets by comparing
 * the size of their intersection to the size of their union. The result is a value
 * between 0 and 1, where 0 indicates no overlap and 1 indicates identical sets.
 *
 * @param A First set of integers represented as a vector
 * @param B Second set of integers represented as a vector
 * @return The Jaccard similarity coefficient as a double in the range [0,1]
 *
 * @note This implementation handles input vectors that may contain duplicates by
 *       converting them to sets.
 *
 * @see https://en.wikipedia.org/wiki/Jaccard_index
 */
double jaccard_index(const std::vector<int>& A,
                     const std::vector<int>& B) {
    std::set<int> set_A(A.begin(), A.end());
    std::set<int> set_B(B.begin(), B.end());

    std::set<int> union_set = set_A;
    union_set.insert(set_B.begin(), set_B.end());

    std::set<int> intersection;
    std::set_intersection(set_A.begin(), set_A.end(),
                          set_B.begin(), set_B.end(),
                          std::inserter(intersection, intersection.begin()));

    return static_cast<double>(intersection.size()) / union_set.size();
}


/**
 * @brief Sigmoidal function with lower and upper thresholds
 *
 * Creates a sigmoidal function that smoothly transitions from approximately 0 to
 * approximately 1 between two specified threshold values. This is useful for
 * creating soft thresholds in data processing or for probability-based filtering.
 *
 * @param x Input value to transform
 * @param lower_threshold The lower threshold value where the function begins its
 *        transition from 0 to 1. For values well below this threshold, the function
 *        returns values close to 0.
 * @param upper_threshold The upper threshold value where the function completes its
 *        transition from 0 to 1. For values well above this threshold, the function
 *        returns values close to 1.
 * @param steepness Controls the steepness of the transition between thresholds.
 *        Higher values create a sharper transition. Default is 10.0.
 *
 * @return A value ranging from approximately 0 to approximately 1
 */
double thresholded_sigmoid(
    double x,
    double lower_threshold,
    double upper_threshold,
    double steepness = 10.0
    ) {
    // Calculate midpoint between thresholds
    double midpoint = (lower_threshold + upper_threshold) / 2.0;

    // Calculate scale factor to control steepness of transition
    // Higher steepness values create sharper transitions
    double scale_factor = steepness / (upper_threshold - lower_threshold);

    // Apply logistic function
    return 1.0 / (1.0 + std::exp(-scale_factor * (x - midpoint)));
}

/**
 * @brief Vector version of the thresholded sigmoid function
 *
 * Applies the thresholded sigmoid function to each element of an input vector.
 *
 * @param x Vector of input values to transform
 * @param lower_threshold The lower threshold value where the function begins its
 *        transition from 0 to 1
 * @param upper_threshold The upper threshold value where the function completes its
 *        transition from 0 to 1
 * @param steepness Controls the steepness of the transition between thresholds.
 *        Default is 10.0.
 *
 * @return A vector of values ranging from approximately 0 to approximately 1
 */
std::vector<double> vthresholded_sigmoid(
    const std::vector<double>& x,
    double lower_threshold,
    double upper_threshold,
    double steepness = 10.0
    ) {
    std::vector<double> result(x.size());

    // Calculate midpoint between thresholds
    double midpoint = (lower_threshold + upper_threshold) / 2.0;

    // Calculate scale factor to control steepness of transition
    double scale_factor = steepness / (upper_threshold - lower_threshold);

    // Apply logistic function to each element
    for (size_t i = 0; i < x.size(); ++i) {
        result[i] = 1.0 / (1.0 + std::exp(-scale_factor * (x[i] - midpoint)));
    }

    return result;
}

/**
* @brief Piece-wise linear thresholding function
*
* Creates a thresholding function that transitions linearly from 0 to 1
* between two specified threshold values. This function returns:
*   - 0 for values less than or equal to the lower threshold
*   - 1 for values greater than or equal to the upper threshold
*   - A linear interpolation for values between the two thresholds
*
* This is useful for creating simple threshold transitions in data processing
* or for filtering based on threshold values.
*
* @param x Input value to transform
* @param lower_threshold The lower threshold value below which the function returns 0
* @param upper_threshold The upper threshold value above which the function returns 1
*
* @return A value ranging from 0 to 1
*/
double thresholded_linear(
    double x,
    double lower_threshold,
    double upper_threshold
    ) {
   // For values below or equal to lower threshold, return 0
   if (x <= lower_threshold) {
       return 0.0;
   }

   // For values above or equal to upper threshold, return 1
   if (x >= upper_threshold) {
       return 1.0;
   }

   // For values between thresholds, return linear interpolation
   return (x - lower_threshold) / (upper_threshold - lower_threshold);
}

/**
 * @brief Calculates evenness of path segments using Shannon entropy
 *
 * This function computes the evenness of a path by measuring how uniformly the total
 * path length is distributed across its constituent segments. It uses normalized
 * Shannon entropy as the evenness measure, where 1.0 represents perfectly even
 * distribution and lower values indicate increasing unevenness.
 *
 * @param edge_lengths Vector of individual segment lengths along the path
 * @return double Evenness value in range [0,1] where:
 *         - 1.0: Perfect evenness (all segments have equal length)
 *         - Near 0: High unevenness (one segment dominates the total length)
 *
 * @note Returns 1.0 for single-edge paths or zero-length paths by definition
 * @note Uses log base 2 for entropy calculation
 */
double calculate_path_evenness(
    const std::vector<double>& edge_lengths
    ) {

    // Early return for single-edge paths
    if (edge_lengths.size() <= 1) {
        return 1.0;  // By definition, a single edge path has perfect evenness
    }

    // Get total length for normalization
    double total_length = 0.0;
    for (const auto& length : edge_lengths) {
        total_length += length;
    }

    // Safety check to avoid division by zero
    if (total_length <= 0.0) {
        return 1.0;  // Default to perfect evenness for zero-length paths
    }

    // Calculate normalized lengths (creating a probability distribution)
    std::vector<double> normalized_lengths;
    normalized_lengths.reserve(edge_lengths.size());
    for (const auto& length : edge_lengths) {
        normalized_lengths.push_back(length / total_length);
    }

    // Calculate entropy
    double entropy = 0.0;
    for (const auto& p : normalized_lengths) {
        if (p > 0.0) {  // Avoid log(0)
            entropy -= p * std::log2(p);
        }
    }

    // Calculate maximum possible entropy for this number of edges
    // For n equally likely outcomes, max entropy is log2(n)
    double max_entropy = std::log2(edge_lengths.size());

    // Safety check to avoid division by zero
    if (max_entropy <= 0.0) {
        return 1.0;
    }

    // Return normalized entropy (evenness)
    return entropy / max_entropy;
}

/**
 * @brief Calculates maximum thresholded edge length in a path
 *
 * This function determines the maximum value among all edge lengths after applying
 * a threshold function to each. It represents an L^∞ (L-infinity) norm analogue to
 * path evenness, capturing the worst-case scenario in path segment distribution.
 *
 * @param edge_lengths Vector of individual segment lengths along the path
 * @param lower_threshold Minimum threshold value below which lengths are mapped to 0
 * @param upper_threshold Maximum threshold value above which lengths are mapped to 1
 * @return double Maximum thresholded value across all edges, in range [0,1]
 *
 * @note Values between lower_threshold and upper_threshold are linearly interpolated
 * @note If all edges are below lower_threshold, returns 0
 * @see thresholded_linear() For the specific thresholding applied to each length
 */
double calculate_path_max_threshold(
    const std::vector<double>& edge_lengths,
    double lower_threshold,
    double upper_threshold
) {
    double max_thldd_q = 0;
    for (const auto& length : edge_lengths) {
        double thldd_q = thresholded_linear(length, lower_threshold, upper_threshold);
        if (thldd_q > max_thldd_q) {
            max_thldd_q = thldd_q;
        }
    }

    return max_thldd_q;
}

double calculate_monotonicity_index(double total_change, double cumulative_absolute_changes) {
    // Safety check to avoid division by zero
    if (cumulative_absolute_changes <= 0.0) {
        return 1.0;  // Default to perfect monotonicity for zero changes
    }

    // Calculate monotonicity index
    return std::abs(total_change) / cumulative_absolute_changes;
}
