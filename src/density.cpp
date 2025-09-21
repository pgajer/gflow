#include <R.h>
#include <R_ext/Rdynload.h>

// #include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/LU>

#include <stdexcept>
#include <limits>
#include <map>

#include "density.hpp"
#include "kernels.h"
#include "SEXP_cpp_conversion_utils.hpp"

static const char* kernel_name(int t) {
    switch (t) {
    case 1: return "Gaussian";
    case 2: return "Epanechnikov";
    case 3: return "Uniform";
    case 4: return "Triangular";
    case 5: return "Biweight";
    case 6: return "Triweight";
    case 7: return "Cosine";
    case 8: return "Logistic";   // adjust names to your package’s convention
    default: return "Unknown";
    }
}

extern "C" {
    SEXP S_estimate_local_density_over_grid(SEXP s_x,
                                            SEXP s_grid_size,
                                            SEXP s_poffset,
                                            SEXP s_pilot_bandwidth,
                                            SEXP s_kernel_type,
                                            SEXP s_verbose);
}

double integrate_kernel_mass(int kernel_type) {
    const int n_points = 10000;  // High precision
    const double dx = 1.0 / n_points;
    initialize_kernel(kernel_type, 1.0);

    std::vector<double> x(n_points + 1);
    std::vector<double> w(n_points + 1);
    for(int i = 0; i <= n_points; i++)
        x[i] = i * dx;

    kernel_fn(x.data(), n_points + 1, w.data());

    double sum = w[0] + w[n_points];
    for(int i = 1; i < n_points; i++)
        sum += 2.0 * w[i];

    return sum * dx;
}


/**
 * @brief Computes the standard deviation of a vector of doubles
 *
 * @param x Input vector of data points
 * @return Standard deviation of the data
 * @pre x must not be empty
 */
double compute_std_dev(const std::vector<double>& x) {
    int n = x.size();
    if (n < 1) return 0.0;

    // Compute mean
    double mean = 0.0;
    for (double val : x) {
        mean += val;
    }
    mean /= n;

    // Compute variance
    double variance = 0.0;
    for (double val : x) {
        double diff = val - mean;
        variance += diff * diff;
    }
    variance /= (n - 1);  // Use n-1 for sample standard deviation

    return std::sqrt(variance);
}

/**
 * @brief Computes the interquartile range (IQR) of a vector of doubles
 *
 * @param x Input vector of data points
 * @return Interquartile range (Q3 - Q1)
 * @pre x must not be empty
 */
double compute_iqr(const std::vector<double>& x) {
    std::vector<double> sorted_x = x;  // Create a copy for sorting
    std::sort(sorted_x.begin(), sorted_x.end());

    int n = sorted_x.size();
    if (n < 4) return 0.0;  // Need at least 4 points for meaningful quartiles

    // Find Q1 (25th percentile)
    double q1_pos = 0.25 * (n - 1);
    int q1_index = static_cast<int>(q1_pos);
    double q1 = sorted_x[q1_index] + (q1_pos - q1_index) *
                (sorted_x[q1_index + 1] - sorted_x[q1_index]);

    // Find Q3 (75th percentile)
    double q3_pos = 0.75 * (n - 1);
    int q3_index = static_cast<int>(q3_pos);
    double q3 = sorted_x[q3_index] + (q3_pos - q3_index) *
                (sorted_x[q3_index + 1] - sorted_x[q3_index]);

    return q3 - q1;
}

/**
 * @brief Computes optimal bandwidth using Silverman's rule of thumb
 *
 * @details Implements Silverman's rule of thumb for bandwidth selection:
 * h = 0.9 * min(σ, IQR/1.34) * n^(-1/5)
 * where:
 * - σ is the standard deviation
 * - IQR is the interquartile range
 * - n is the sample size
 *
 * @param x Input vector of data points
 * @return Optimal bandwidth according to Silverman's rule
 * @pre x must not be empty
 * @throw std::invalid_argument if x is empty
 *
 * @note If IQR is zero (e.g., due to identical quartiles),
 * falls back to using only the standard deviation
 */
double silverman_bandwidth(const std::vector<double>& x) {
    if (x.empty()) {
        Rf_error("Input vector must not be empty");
    }

    int n = x.size();
    double sd = compute_std_dev(x);
    double iqr = compute_iqr(x);

    // If IQR is zero, use only sd-based estimate
    if (iqr < std::numeric_limits<double>::epsilon()) {
        return 0.9 * sd * std::pow(n, -0.2);
    }

    double iqr_based = iqr / 1.34;  // Scale IQR as per Silverman's rule
    double min_spread = std::min(sd, iqr_based);

    return 0.9 * min_spread * std::pow(n, -0.2);
}


/**
 * @brief Computes optimal bandwidth using Silverman's rule adjusted for different kernels
 *
 * @details Implements a modified Silverman's rule for different kernel types:
 * h = C(K) * min(σ, IQR/1.34) * n^(-1/5)
 * where:
 * - C(K) is the kernel-specific constant
 * - σ is the standard deviation
 * - IQR is the interquartile range
 * - n is the sample size
 *
 * Kernel adjustment factors from theoretical calculations:
 * - Normal (Gaussian): 0.9
 * - Epanechnikov: 2.34
 * - Biweight: 2.78
 * - Triangular: 2.58
 * - Rectangular: 2.00
 *
 * @param x Input vector of data points
 * @param kernel_type Type of kernel being used
 * @return Optimal bandwidth adjusted for the specified kernel
 */
double kernel_specific_bandwidth(const std::vector<double>& x, int kernel_type) {
    // Constants for different kernels
    const double NORMAL_CONST = 0.9;
    const double EPANECHNIKOV_CONST = 2.34;
    const double BIWEIGHT_CONST = 2.78;
    const double TRIANGULAR_CONST = 2.58;
    const double TRICUBE_CONST = 3.12;    // Approximate
    const double LAPLACE_CONST = 1.30;    // Approximate

    // Select kernel-specific constant
    double kernel_const;
    switch(kernel_type) {
        case 1: kernel_const = EPANECHNIKOV_CONST; break;
        case 2: kernel_const = TRIANGULAR_CONST; break;
        case 4: kernel_const = LAPLACE_CONST; break;
        case 5: kernel_const = NORMAL_CONST; break;
        case 6: kernel_const = BIWEIGHT_CONST; break;
        case 7: kernel_const = TRICUBE_CONST; break;
        default: kernel_const = NORMAL_CONST;  // Default to normal kernel constant
    }

    int n = x.size();
    double sd = compute_std_dev(x);
    double iqr = compute_iqr(x);

    // If IQR is zero, use only sd-based estimate
    if (iqr < std::numeric_limits<double>::epsilon()) {
        return kernel_const * sd * std::pow(n, -0.2);
    }

    double iqr_based = iqr / 1.34;
    double min_spread = std::min(sd, iqr_based);

    return kernel_const * min_spread * std::pow(n, -0.2);
}

/**
* @brief Estimates the local density at a specified point using kernel density estimation
* with optional automatic bandwidth selection
*
* @details This function implements kernel density estimation with either manual or
* automatic bandwidth selection based on Silverman's rule adjusted for different kernels.
* If pilot_bandwidth <= 0, automatic bandwidth selection is used.
*
* @param x Vector of data points for density estimation
* @param center_idx Index of the point at which to estimate the density
* @param pilot_bandwidth Bandwidth parameter (if <= 0, automatic selection is used)
* @param kernel_type Integer specifying the kernel function (1-8)
* @param verbose If true, prints debugging information
*
* @return density_t structure containing density estimate, bandwidth used, and selection flag
* @throw std::invalid_argument if input parameters are invalid
*/
density_t estimate_local_density(
    const std::vector<double>& x,
    int center_idx,
    double pilot_bandwidth,
    int kernel_type,
    bool verbose) {

    // Input validation
    if (x.empty())                          Rf_error("Input vector must not be empty");
    if (center_idx < 0 || center_idx >= (int)x.size()) Rf_error("center_idx out of range");
    if (kernel_type < 1 || kernel_type > 8) Rf_error("Invalid kernel_type (must be 1-8)");

    density_t result;
    result.auto_selected = (pilot_bandwidth <= 0);
    try {
        result.bandwidth = result.auto_selected
            ? kernel_specific_bandwidth(x, kernel_type)
            : pilot_bandwidth;
    } catch (const std::exception& e) {
        Rf_error("Bandwidth selection failed: %s", e.what());  // no newline
    }

    if (verbose) {
        Rprintf(
            "Density estimation parameters:\n"
            "  - Kernel type: %s\n"
            "  - Bandwidth: %.6g%s\n"
            "  - Sample size: %lld\n"
            "  - Estimation point index: %lld\n",
            kernel_name(kernel_type),
            result.bandwidth,
            (result.auto_selected ? " (automatically selected)" : " (user-provided)"),
            (long long)x.size(),
            (long long)center_idx
            );
    }

    // Initialize the chosen kernel
    initialize_kernel(kernel_type, 1.0);
    double kernel_total_mass = get_kernel_mass(kernel_type);
    double center = x[center_idx];
    double density = 0.0;
    int n = x.size();

    // Compute density using the selected kernel
    std::vector<double> dists(n);
    std::vector<double> weights(n);

    for (int i = 0; i < n; ++i) {
        dists[i] = std::abs(x[i] - center) / result.bandwidth;
    }

    // Apply the selected kernel
    kernel_fn(dists.data(), n, weights.data());

    // Sum the weights
    for (int i = 0; i < n; ++i) {
        density += weights[i];
    }

    result.density = density / (n * result.bandwidth * kernel_total_mass);

    if (verbose) {
        Rprintf(
            "Density estimation results:\n"
            "  - Estimated density: %.6g\n"
            "  - Total kernel mass: %.6g\n",
            result.density,
            kernel_total_mass
            );
    }

    return result;
}

/**
 * @brief Estimates local density of 1D data over a uniform grid using kernel density estimation
 *
 * @param x Input vector of data points
 * @param grid_size Number of points in the output grid (must be positive)
 * @param poffset Proportion of data width to add as offset on each end of the grid
 * @param pilot_bandwidth Kernel bandwidth (if <= 0, bandwidth is automatically selected)
 * @param kernel_type Type of kernel to use (1-8)
 * @param verbose If true, prints detailed information about the estimation process
 *
 * @return gdensity_t Structure containing:
 *         - grid_size: size of the output grid
 *         - offset: actual offset added to data range
 *         - start: starting point of the grid
 *         - end: ending point of the grid
 *         - density: vector of density estimates at grid points
 *         - bandwidth: bandwidth used in the estimation
 *         - auto_selected: whether bandwidth was automatically selected
 *
 * @throw std::invalid_argument if:
 *        - input vector is empty
 *        - grid_size is not positive
 *        - kernel_type is not in range [1,8]
 *
 * @note The function uses external helper functions:
 *       - kernel_specific_bandwidth() for automatic bandwidth selection
 *       - initialize_kernel() to set up the kernel
 *       - get_kernel_mass() to get kernel normalization factor
 *       - kernel_fn() to evaluate the kernel
 */
gdensity_t estimate_local_density_over_grid(
   const std::vector<double>& x,
   int grid_size,
   double poffset,
   double pilot_bandwidth,
   int kernel_type,
   bool verbose) {

   gdensity_t result(grid_size);
   result.auto_selected = (pilot_bandwidth <= 0);

   // Automatic bandwidth selection if needed
   result.bandwidth = result.auto_selected ?
           kernel_specific_bandwidth(x, kernel_type) : pilot_bandwidth;

   double min_x = *std::min_element(x.begin(), x.end());
   double max_x = *std::max_element(x.begin(), x.end());
   double x_width = max_x - min_x;
   result.offset = poffset * x_width;
   result.start = min_x - result.offset;
   result.end = max_x + result.offset;

   if (verbose) {
       Rprintf(
           "Density estimation parameters:\n"
           "  - Sample size: %lld\n"
           "  - Grid Size: %lld\n"
           "  - p-offset: %.6g\n"
           "  - Offset: %.6g\n"
           "  - Grid Start: %.6g\n"
           "  - Grid End: %.6g\n"
           "  - Kernel type: %s\n"
           "  - Bandwidth: %.6g%s\n",
           (long long)x.size(),
           (long long)grid_size,
           poffset,
           result.offset,
           result.start,
           result.end,
           kernel_name(kernel_type),
           result.bandwidth,
           (result.auto_selected ? " (automatically selected)" : " (user-provided)")
           );
   }

   // Initialize the chosen kernel
   initialize_kernel(kernel_type, 1.0);
   double kernel_total_mass = get_kernel_mass(kernel_type);

   // Create uniform grid on [start, end]
   double dx = (result.end - result.start) / (grid_size - 1);
   int n = x.size();
   std::vector<double> dists(n);
   std::vector<double> weights(n);
   double norm = n * result.bandwidth * kernel_total_mass;

   for(int gi = 0; gi < grid_size; gi++) {
       double center = result.start + gi * dx;
       result.density[gi] = 0.0;

       // Compute density using the selected kernel
       for (int i = 0; i < n; ++i)
           dists[i] = std::abs(x[i] - center) / result.bandwidth;

       // Apply the selected kernel
       kernel_fn(dists.data(), n, weights.data());

       // Sum the weights
       for (int i = 0; i < n; ++i) {
           result.density[gi] += weights[i];
       }

       result.density[gi] /= norm;
   }

   if (verbose) {
       Rprintf("Finished density estimation\n");
   }

   return result;
}

/**
 * @brief R interface for kernel density estimation over a uniform grid
 *
 * @details This function wraps the C++ density estimation function for use in R,
 * handling conversion between R and C++ data structures and memory protection.
 *
 * @param s_x SEXP containing numeric vector of data points
 * @param s_grid_size SEXP containing integer scalar for grid size
 * @param s_poffset SEXP containing numeric scalar for proportion of data width to add as offset
 * @param s_pilot_bandwidth SEXP containing numeric scalar for kernel bandwidth (<=0 for automatic selection)
 * @param s_kernel_type SEXP containing integer scalar for kernel type (1-8)
 * @param s_verbose SEXP containing logical scalar for verbose output
 *
 * @return SEXP (list) containing:
 *   - density: numeric vector of density estimates at grid points
 *   - bw: numeric scalar of bandwidth used
 *   - bw_auto_selected: logical indicating if bandwidth was automatically selected
 *   - offset: numeric scalar of actual offset used
 *   - start: numeric scalar of grid start point
 *   - end: numeric scalar of grid end point
 *
 * @note Memory protection is handled internally for all SEXP objects
 */
SEXP S_estimate_local_density_over_grid(SEXP s_x,
                                        SEXP s_grid_size,
                                        SEXP s_poffset,
                                        SEXP s_pilot_bandwidth,
                                        SEXP s_kernel_type,
                                        SEXP s_verbose) {
    // --- Coerce x to REAL and copy into std::vector<double> ---
    std::vector<double> x;
    {
        SEXP sx = s_x;
        PROTECT_INDEX ipx;
        PROTECT_WITH_INDEX(sx, &ipx);
        if (TYPEOF(sx) != REALSXP) {
            REPROTECT(sx = Rf_coerceVector(sx, REALSXP), ipx);
        }
        const R_xlen_t nx = XLENGTH(sx);
        if (nx <= 0) {
            UNPROTECT(1);
            Rf_error("S_estimate_local_density_over_grid(): 'x' must be a non-empty numeric vector.");
        }
        x.assign(REAL(sx), REAL(sx) + (size_t)nx);
        UNPROTECT(1); // sx
    }

    // --- Scalars (defensive) ---
    const int    grid_size_i      = Rf_asInteger(s_grid_size);
    const double poffset          = Rf_asReal(s_poffset);
    const double pilot_bandwidth  = Rf_asReal(s_pilot_bandwidth);
    const int    kernel_type_i    = Rf_asInteger(s_kernel_type);
    const int    verbose_i        = Rf_asLogical(s_verbose);

    // --- NA / range checks ---
    if (grid_size_i == NA_INTEGER || grid_size_i <= 0) {
        Rf_error("S_estimate_local_density_over_grid(): 'grid_size' must be a positive integer.");
    }
    if (ISNAN(poffset)) {
        Rf_error("S_estimate_local_density_over_grid(): 'poffset' cannot be NA.");
    }
    if (ISNAN(pilot_bandwidth) || pilot_bandwidth < 0.0) {
        Rf_error("S_estimate_local_density_over_grid(): 'pilot_bandwidth' must be >= 0.");
    }
    if (kernel_type_i == NA_INTEGER || kernel_type_i < 0) {
        Rf_error("S_estimate_local_density_over_grid(): 'kernel_type' must be a non-negative integer.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_estimate_local_density_over_grid(): 'verbose' must be TRUE/FALSE.");
    }

    const int  grid_size     = grid_size_i;
    const int  kernel_type   = kernel_type_i;
    const bool verbose       = (verbose_i == TRUE);

    // --- Core computation (no R allocations inside) ---
    gdensity_t gdens_res = estimate_local_density_over_grid(
        x, grid_size, poffset, pilot_bandwidth, kernel_type, verbose);

    // --- Build result list (container-first) ---
    const int N_COMPONENTS = 6;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    {
        SEXP r_names  = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(r_names, 0, Rf_mkChar("y"));
        SET_STRING_ELT(r_names, 1, Rf_mkChar("bw"));
        SET_STRING_ELT(r_names, 2, Rf_mkChar("bw_auto_selected"));
        SET_STRING_ELT(r_names, 3, Rf_mkChar("offset"));
        SET_STRING_ELT(r_names, 4, Rf_mkChar("start"));
        SET_STRING_ELT(r_names, 5, Rf_mkChar("end"));
        Rf_setAttrib(r_result, R_NamesSymbol, r_names);
        UNPROTECT(1); // r_names
    }

    // 0: y (density)
    {
        SEXP el0 = PROTECT(convert_vector_double_to_R(gdens_res.density));
        SET_VECTOR_ELT(r_result, 0, el0);
        UNPROTECT(1);
    }
    // 1: bw
    {
        SEXP s_bw = PROTECT(Rf_ScalarReal(gdens_res.bandwidth));
        SET_VECTOR_ELT(r_result, 1, s_bw);
        UNPROTECT(1);
    }
    // 2: bw_auto_selected
    {
        SEXP s_auto = PROTECT(Rf_ScalarLogical(gdens_res.auto_selected ? TRUE : FALSE));
        SET_VECTOR_ELT(r_result, 2, s_auto);
        UNPROTECT(1);
    }
    // 3: offset
    {
        SEXP s_off = PROTECT(Rf_ScalarReal(gdens_res.offset));
        SET_VECTOR_ELT(r_result, 3, s_off);
        UNPROTECT(1);
    }
    // 4: start
    {
        SEXP s_start = PROTECT(Rf_ScalarReal(gdens_res.start));
        SET_VECTOR_ELT(r_result, 4, s_start);
        UNPROTECT(1);
    }
    // 5: end
    {
        SEXP s_end = PROTECT(Rf_ScalarReal(gdens_res.end));
        SET_VECTOR_ELT(r_result, 5, s_end);
        UNPROTECT(1);
    }

    UNPROTECT(1); // r_result
    return r_result;
}
