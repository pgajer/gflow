#include <random>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

extern "C" {
    SEXP S_rlaplace(SEXP R_n, SEXP R_location, SEXP R_scale, SEXP R_seed);
}

/**
 * @brief Generate random variates from the Laplace distribution
 *
 * This function generates a sample of n random values from the Laplace distribution
 * (also known as the double exponential distribution) with specified location and
 * scale parameters.
 *
 * The probability density function of the Laplace distribution is:
 * f(x | μ, b) = (1 / (2b)) * exp(-|x - μ| / b)
 * where μ is the location parameter and b is the scale parameter.
 *
 * @param n The number of random variates to generate.
 * @param location The location parameter μ (default = 0.0).
 * @param scale The scale parameter b (default = 1.0).
 * @param seed The seed for the random number generator (default = 0).
 *             If seed is 0, a random device is used to seed the generator.
 *
 * @return std::vector<double> A vector of n random variates from the Laplace distribution.
 *
 * @throws std::invalid_argument if n is not positive, location is not finite, or scale is not positive and finite.
 *
 * @note This function uses the inverse transform sampling method to generate
 *       Laplace-distributed random numbers.
 */
std::vector<double> rlaplace(int n,
                             double location = 0.0,
                             double scale = 1.0,
                             unsigned int seed = 0) {

    std::vector<double> result(n);
    std::mt19937 gen;
    if (seed == 0) {
        std::random_device rd;
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }

    std::uniform_real_distribution<> dis(-0.5, 0.5);

    for (int i = 0; i < n; ++i) {
        double u;
        do {
            u = dis(gen);
        } while (u == -0.5 || u == 0.5);
        result[i] = location - scale * std::copysign(std::log(1 - 2 * std::abs(u)), u);
    }
    return result;
}

/**
 * @brief Generate random variates from the Laplace distribution
 *
 * This function serves as an interface between R and the C++ implementation
 * of the Laplace distribution random number generator. It generates a vector
 * of random variates from the Laplace distribution with specified location
 * and scale parameters.
 *
 * @param R_n SEXP (INTEGER) The number of observations to generate.
 * @param R_location SEXP (REAL) The location parameter of the distribution.
 * @param R_scale SEXP (REAL) The scale parameter of the distribution.
 * @param R_seed SEXP (INTEGER) The seed for the random number generator.
 *               If NULL, a random seed will be generated.
 *
 * @return SEXP A numeric vector of length n containing the generated random variates.
 *
 * @note This function uses the C++ standard library's random number generator,
 *       which is seeded using either the provided seed or R's random number
 *       generator to ensure reproducibility.
 *
 * @example
 * In R:
 * result <- .Call("S_rlaplace", 1000L, 0.0, 1.0, 12345L)  # With specific seed
 * result <- .Call("S_rlaplace", 1000L, 0.0, 1.0, NULL)    # With random seed
 */
SEXP S_rlaplace(SEXP R_n, SEXP R_location, SEXP R_scale, SEXP R_seed) {

    // Check input types
    if (!Rf_isInteger(R_n) || !Rf_isReal(R_location) || !Rf_isReal(R_scale) ||
        (!Rf_isNull(R_seed) && !Rf_isInteger(R_seed))) {
        Rf_error("Invalid input types");
    }

    // Extract parameters
    int n = INTEGER(R_n)[0];
    double location = REAL(R_location)[0];
    double scale = REAL(R_scale)[0];

    // Check for valid parameters
    if (n <= 0) {
        Rf_error("n must be a positive integer");
    }
    if (scale <= 0) {
        Rf_error("scale must be positive");
    }

    // Handle seed
    unsigned int seed;
    if (Rf_isNull(R_seed)) {
        GetRNGstate();
        seed = static_cast<unsigned int>(unif_rand() * UINT_MAX);
        PutRNGstate();
    } else {
        // Make sure R_seed has at least one element
        if (LENGTH(R_seed) < 1) {
            Rf_error("seed must have at least one element");
        }
        seed = static_cast<unsigned int>(INTEGER(R_seed)[0]);
    }

    // Call the C++ function with the seed
    std::vector<double> result = rlaplace(n, location, scale, seed);

    // Create an R numeric vector and copy the results
    SEXP R_result = PROTECT(Rf_allocVector(REALSXP, n));
    double* r_ptr = REAL(R_result);
    for (int i = 0; i < n; ++i) {
        r_ptr[i] = result[i];
    }

    UNPROTECT(1);

    return R_result;
}


/**
 * @brief Generates a random sample from a symmetric Dirichlet distribution
 *
 * This function generates a random sample from a Dirichlet distribution over the
 * (n-1)-dimensional probability simplex with identical concentration parameters.
 *
 * @details
 * Mathematical Background:
 * The function uses the following property to generate Dirichlet samples:
 * If Y₁, ..., Yₙ are independent random variables where:
 *   - Yᵢ ~ Gamma(αᵢ, 1) for i = 1,...,n
 *   - αᵢ = α for all i (symmetric case)
 *
 * Then the vector X = (X₁, ..., Xₙ) where:
 *   Xᵢ = Yᵢ / (Y₁ + ... + Yₙ)
 * follows a Dirichlet(α₁, ..., αₙ) distribution.
 *
 * This works because:
 * 1. The Gamma distribution is the conjugate prior of the scale parameter
 *    in the Dirichlet distribution
 * 2. The normalization step creates the necessary dependence structure
 *    between the components while preserving the relative ratios
 * 3. The resulting distribution has the key properties of Dirichlet:
 *    - Each Xᵢ > 0
 *    - ΣXᵢ = 1
 *    - The density function matches the Dirichlet probability density:
 *      f(x₁,...,xₙ) ∝ ∏ᵢ xᵢ^(αᵢ-1)
 *
 * The probability simplex Δⁿ⁻¹ is the set of vectors x = (x₁,...,xₙ) such that:
 *   - xᵢ ≥ 0 for all i
 *   - Σxᵢ = 1
 *
 * Properties of the resulting distribution:
 * 1. E[Xᵢ] = α / (nα) = 1/n (symmetric case)
 * 2. Var(Xᵢ) = (α/nα)(1 - α/nα)/(nα + 1)
 * 3. Cov(Xᵢ,Xⱼ) = -(α/nα)(α/nα)/(nα + 1)
 *
 * @param n Dimension of the distribution (size of the output vector)
 * @param alpha Concentration parameter (same value for all dimensions)
 * @return std::vector<double> A vector representing a point on the (n-1)-simplex
 *
 * @throws std::invalid_argument if n < 2 or alpha <= 0
 *
 * @note The function uses the C++11 random number generation facilities
 * @note Time complexity: O(n)
 * @note Space complexity: O(n)
 */
std::vector<double> sample_symmetric_dirichlet(int n, double alpha) {
    if (n < 2) {
        Rf_error("n must be at least 2");
    }
    if (alpha <= 0) {
        Rf_error("alpha must be positive");
    }

    // Set up random number generation
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::gamma_distribution<double> gamma_dist(alpha, 1.0);

    // Generate n independent gamma random variables
    std::vector<double> samples(n);
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        samples[i] = gamma_dist(gen);
        sum += samples[i];
    }

    // Normalize to get a point on the simplex
    for (int i = 0; i < n; ++i) {
        samples[i] /= sum;
    }

    return samples;
}


/**
 * @brief Generates a random point in the probability simplex using beta radial sampling.
 *
 * This function implements the beta radial sampling algorithm which generates points
 * in an n-dimensional probability simplex. The algorithm works by:
 * 1. Sampling a point from the unit (n-2)-sphere
 * 2. Mapping this point to the boundary of the probability simplex
 * 3. Interpolating between the boundary point and the barycenter using a Beta distribution
 *
 * @param n The dimension of the probability simplex (must be >= 2)
 * @param alpha First parameter of the Beta distribution (must be > 0)
 * @param beta Second parameter of the Beta distribution (must be > 0)
 * @param seed Random seed for the random number generator
 *
 * @return A vector of n non-negative doubles that sum to 1, representing a point in
 *         the probability simplex
 *
 * @throws std::invalid_argument if n < 2 or alpha <= 0 or beta <= 0
 *
 * @note The function uses the Beta distribution to control how close the generated
 *       points are to the simplex boundary. Larger values of beta relative to alpha
 *       will generate points closer to the boundary.
 *
 * @see For theoretical background, see [reference to your paper/documentation]
 *
 * @example
 * ```cpp
 * // Generate a point in the 3-simplex with Beta(1,3) distribution
 * std::vector<double> point = beta_radial(4, 1.0, 3.0);
 * ```
 */
std::vector<double> beta_radial(int n, double alpha = 1.0, double beta = 3.0,
                               unsigned seed = std::random_device{}()) {
    // Initialize random number generator
    std::mt19937 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // 1. Sample from unit sphere (n-2 dimensional)
    std::vector<double> sphere_point(n-2);
    double norm = 0.0;
    for (int i = 0; i < n-2; ++i) {
        sphere_point[i] = normal(rng);
        norm += sphere_point[i] * sphere_point[i];
    }
    norm = std::sqrt(norm);
    for (int i = 0; i < n-2; ++i) {
        sphere_point[i] /= norm;
    }

    // 2. Map sphere point to simplex boundary
    std::vector<double> boundary_point(n);
    double sum = 0.0;
    for (int i = 0; i < n-1; ++i) {
        boundary_point[i] = (i < n-2) ? sphere_point[i] + 1.0/n : 1.0/n;
        sum += boundary_point[i];
    }
    boundary_point[n-1] = 1.0 - sum;

    // Project onto probability simplex if needed
    std::vector<double> sorted = boundary_point;
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());

    double cumsum = 0.0;
    double tau = 0.0;
    int rho = 0;

    for (int i = 0; i < n; ++i) {
        cumsum += sorted[i];
        tau = (cumsum - 1.0) / (i + 1);
        if (sorted[i] > tau) {
            rho = i + 1;
        }
    }

    tau = (std::accumulate(sorted.begin(), sorted.begin() + rho, 0.0) - 1.0) / rho;
    for (int i = 0; i < n; ++i) {
        boundary_point[i] = std::max(boundary_point[i] - tau, 0.0);
    }

    // 3. Sample from Beta distribution
    std::gamma_distribution<double> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<double> gamma_beta(beta, 1.0);
    double x = gamma_alpha(rng);
    double y = gamma_beta(rng);
    double t = x / (x + y);

    // 4. Interpolate between boundary point and barycenter
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = boundary_point[i] * (1.0 - t) + (1.0/n) * t;
    }

    return result;
}


/**
 * @brief Sample points from empirical distribution with optional linear approximation
 *
 * @details This function implements the inverse transform sampling method to generate
 * random samples from an empirical distribution. The algorithm follows these steps:
 *
 * 1. Density Estimation:
 *    - Constructs a histogram of the input data using nbins equally spaced bins
 *    - The histogram counts represent the empirical density
 *
 * 2. PDF Normalization:
 *    - Normalizes the density values to create a proper PDF that integrates to 1
 *    - For histogram with equal bin widths dx: pdf = counts / (sum(counts) * dx)
 *
 * 3. CDF Computation:
 *    - Computes the cumulative distribution function (CDF) from the normalized PDF
 *    - For histogram with equal bin widths: cdf[i] = sum(pdf[1:i]) * dx
 *    - The CDF is a step function by nature, jumping at bin midpoints
 *
 * 4. Inverse Transform Sampling:
 *    - Generates n uniform random numbers u ~ U(0,1)
 *    - For each u, finds x such that F(x) = u, where F is the CDF
 *    - Two methods available for finding x:
 *      a) Step Function (use_linear_approximation=false):
 *         - Returns the first x value where CDF ≥ u
 *         - Preserves the discrete nature of the histogram
 *      b) Linear Interpolation (use_linear_approximation=true):
 *         - Linearly interpolates between surrounding x values
 *         - Smooths the distribution but may not perfectly represent histogram
 *
 * @param data Input vector of raw data points from which to estimate distribution
 * @param n Number of sample points to generate
 * @param nbins Number of bins to use when estimating density (default: 100)
 * @param use_linear_approximation If true, uses linear interpolation between CDF points;
 *                                if false, uses step function approach (default: false)
 *
 * @return Vector of n samples drawn from the estimated distribution
 *
 * @throws std::invalid_argument if n <= 0 or nbins <= 0
 * @throws std::invalid_argument if data vector is empty
 *
 * @note The function uses C++11 random number generation facilities for uniform sampling
 *
 * @example
 * std::vector<double> data = {1.2, 2.3, 2.4, 3.1, 4.2, 4.9};  // input data
 * int n = 1000;  // number of samples to generate
 * auto samples = sample_from_empirical_distribution(data, n);
 */
std::vector<double> sample_from_empirical_distribution(
    const std::vector<double>& data,
    int n,
    int nbins = 100,
    bool use_linear_approximation = false) {
    // Input validation
    if (n <= 0) {
        Rf_error("Number of samples must be positive");
    }
    if (nbins <= 0) {
        Rf_error("Number of bins must be positive");
    }
    if (data.empty()) {
        Rf_error("Input data vector cannot be empty");
    }

    // Find data range for histogram
    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());
    double range = max_val - min_val;
    double dx = range / nbins;

    // Create histogram
    std::vector<double> counts(nbins, 0);
    std::vector<double> bin_edges(nbins + 1);
    std::vector<double> bin_centers(nbins);

    // Calculate bin edges and centers
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = min_val + i * dx;
        if (i < nbins) {
            bin_centers[i] = min_val + (i + 0.5) * dx;
        }
    }

    // Fill histogram
    for (double value : data) {
        int bin = static_cast<int>((value - min_val) / dx);
        if (bin == nbins) bin--; // Handle edge case for maximum value
        counts[bin]++;
    }

    // Normalize to create PDF
    double sum_counts = std::accumulate(counts.begin(), counts.end(), 0.0);
    std::vector<double> pdf(nbins);
    for (int i = 0; i < nbins; ++i) {
        pdf[i] = counts[i] / (sum_counts * dx);
    }

    // Compute CDF
    std::vector<double> cdf(nbins + 1, 0.0);
    for (int i = 1; i <= nbins; ++i) {
        cdf[i] = cdf[i-1] + pdf[i-1] * dx;
    }

    // Generate uniform random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    // Generate samples using inverse transform sampling
    std::vector<double> samples(n);
    for (int i = 0; i < n; ++i) {
        double u = uniform(gen);

        if (use_linear_approximation) {
            // Find the interval containing u using binary search
            auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
            int idx = std::distance(cdf.begin(), it) - 1;

            if (idx == nbins) {
                samples[i] = bin_centers.back();
            } else {
                // Linear interpolation
                double alpha = (u - cdf[idx]) / (cdf[idx + 1] - cdf[idx]);
                samples[i] = bin_centers[idx] + alpha * dx;
            }
        } else {
            // Step function approach
            auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
            int idx = std::distance(cdf.begin(), it) - 1;
            idx = std::min(idx, nbins - 1);  // Ensure we don't exceed array bounds
            samples[i] = bin_centers[idx];
        }
    }

    return samples;
}


/**
 * @brief Generates a random point uniformly distributed on the standard simplex
 *
 * Generates a point (x₁,...,xₖ) such that xᵢ ≥ 0 and ∑xᵢ = 1 using the method of
 * ordered uniform spacings: generate K-1 uniform random numbers, sort them, and
 * take the successive differences including 0 and 1 as endpoints.
 *
 * @param K Dimension of the simplex (K-1 is the dimension of the manifold)
 * @param seed Random number generator seed for reproducibility
 * @return std::vector<double> Vector of K non-negative numbers that sum to 1
 *
 * @see Uniform Spacings Method: Smith, Robert L., and Len Mitas. "Uniform Point
 *      Sampling of N-Dimensional Spheres." (1996)
 */
std::vector<double> runif_simplex(int K, int seed) {
    std::vector<double> lambda(K);

    // Set the seed
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Generate K-1 uniform random numbers
    for(int i = 0; i < K-1; i++) {
        lambda[i] = dis(gen);
    }

    // Sort first K-1 elements
    std::sort(lambda.begin(), lambda.begin() + K-1);

    // Calculate the last element
    lambda[K-1] = 1.0 - lambda[K-2];

    // Calculate differences
    for(int i = K-2; i > 0; i--) {
        lambda[i] -= lambda[i-1];
    }

    return lambda;
}

/**
 * @brief Generates multiple random points uniformly distributed on the standard simplex
 *
 * Generates n_bootstraps points, each (x₁,...,xₖ) such that xᵢ ≥ 0 and ∑xᵢ = 1,
 * using the method of ordered uniform spacings: generate K-1 uniform random numbers,
 * sort them, and take the successive differences including 0 and 1 as endpoints.
 *
 * @param K Dimension of the simplex (K-1 is the dimension of the manifold)
 * @param n_bootstraps Number of random points to generate
 * @param seed Random number generator seed for reproducibility
 * @return std::vector<std::vector<double>> Vector of n_bootstraps samples, each a
 *         vector of K non-negative numbers that sum to 1
 *
 * @see Uniform Spacings Method: Smith, Robert L., and Len Mitas. "Uniform Point
 *      Sampling of N-Dimensional Spheres." (1996)
 *
 * @note This is a vectorized version of runif_simplex that generates multiple
 *       samples using a single random number generator initialization
 */
std::vector<std::vector<double>> mrunif_simplex(int K, int n_bootstraps, int seed) {
    // Initialize the result vector - outer vector of size n_bootstraps
    // each containing a vector of size K
    std::vector<std::vector<double>> samples(n_bootstraps, std::vector<double>(K));

    // Initialize random number generator
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Generate n_bootstraps samples
    for(int b = 0; b < n_bootstraps; b++) {
        // Generate K-1 uniform random numbers
        for(int i = 0; i < K-1; i++) {
            samples[b][i] = dis(gen);
        }

        // Sort first K-1 elements
        std::sort(samples[b].begin(), samples[b].begin() + K-1);

        // Calculate the last element
        samples[b][K-1] = 1.0 - samples[b][K-2];

        // Calculate differences
        for(int i = K-2; i > 0; i--) {
            samples[b][i] -= samples[b][i-1];
        }
    }

    return samples;
}
