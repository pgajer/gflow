/**
 * @file lcor_posterior_r.cpp
 * @brief Memory-efficient local correlation computation with posterior uncertainty propagation
 *
 * This implementation computes lcor statistics with full posterior uncertainty
 * quantification without materializing all posterior samples in memory. For each
 * feature, posterior samples are generated, lcor is computed for each sample,
 * summary statistics are computed, and samples are discarded before moving to
 * the next feature.
 *
 * Memory complexity: O(n × B) per feature, O(n × p) total output
 * vs. O(n × B × p) if all samples were stored
 *
 * OpenMP parallelization is supported for the outer loop over features.
 */

#include "set_wgraph.hpp"
#include "lcor.hpp"

#include <R.h>
#include <Rinternals.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// Type aliases
using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

/**
 * @brief Parse lcor type string to enum
 */
static lcor_type_t parse_lcor_type(const char* type_str) {
    std::string type(type_str);
    if (type == "unit") return lcor_type_t::UNIT;
    if (type == "derivative") return lcor_type_t::DERIVATIVE;
    if (type == "sign") return lcor_type_t::SIGN;
    Rf_error("Unknown lcor type: %s", type_str);
    return lcor_type_t::DERIVATIVE; // Never reached
}

/**
 * @brief Compute filter weight for a single eigenvalue
 */
static double compute_filter_weight(double lambda, double eta, const char* filter_type) {
    std::string ftype(filter_type);

    if (ftype == "heat_kernel") {
        return std::exp(-eta * lambda);
    } else if (ftype == "tikhonov") {
        return 1.0 / (1.0 + eta * lambda);
    } else if (ftype == "cubic_spline") {
        return 1.0 / (1.0 + eta * lambda * lambda);
    } else if (ftype == "gaussian") {
        return std::exp(-eta * lambda * lambda);
    } else if (ftype == "exponential") {
        return std::exp(-eta * std::sqrt(std::max(lambda, 0.0)));
    } else if (ftype == "butterworth") {
        double x = lambda / eta;
        return 1.0 / (1.0 + x * x * x * x);
    }

    // Default to heat kernel
    return std::exp(-eta * lambda);
}

/**
 * @brief Compute filter weights for all eigenvalues
 */
static vec_t compute_filter_weights(const vec_t& eigenvalues, double eta,
                                    const char* filter_type) {
    const int m = eigenvalues.size();
    vec_t weights(m);
    for (int j = 0; j < m; ++j) {
        weights[j] = compute_filter_weight(eigenvalues[j], eta, filter_type);
    }
    return weights;
}

/**
 * @brief Select optimal eta via GCV for a single response
 *
 * @param V Eigenvector matrix (n x m)
 * @param eigenvalues Raw eigenvalues (length m)
 * @param z Response vector (length n)
 * @param filter_type Filter type string
 * @param n_candidates Number of eta candidates
 *
 * @return Optimal eta value
 */
static double select_eta_gcv(const mat_t& V, const vec_t& eigenvalues,
                             const vec_t& z, const char* filter_type,
                             int n_candidates) {
    const int n = V.rows();
    // const int m = V.cols();

    // Project onto spectral basis
    vec_t Vt_z = V.transpose() * z;

    // Generate log-spaced eta grid
    const double eps = 1e-10;
    double lambda_max = eigenvalues.maxCoeff();
    double eta_min = eps;
    double eta_max = (lambda_max > 0) ? -std::log(eps) / lambda_max : 1.0;

    std::vector<double> eta_grid(n_candidates);
    for (int k = 0; k < n_candidates; ++k) {
        double t = static_cast<double>(k + 1) / n_candidates;
        eta_grid[k] = std::exp(std::log(eta_min) + t * (std::log(eta_max) - std::log(eta_min)));
    }

    // Find optimal eta
    double best_gcv = std::numeric_limits<double>::max();
    double best_eta = eta_grid[0];

    for (int k = 0; k < n_candidates; ++k) {
        double eta = eta_grid[k];
        vec_t f_lambda = compute_filter_weights(eigenvalues, eta, filter_type);

        // Filtered coefficients
        vec_t alpha_filtered = f_lambda.array() * Vt_z.array();

        // Fitted values
        vec_t z_hat = V * alpha_filtered;

        // RSS
        double rss = (z - z_hat).squaredNorm();

        // Effective df = sum of filter weights
        double trace_S = f_lambda.sum();

        // GCV score
        double denom = std::max(n - trace_S, 1.0);
        double gcv = rss / (denom * denom);

        if (gcv < best_gcv) {
            best_gcv = gcv;
            best_eta = eta;
        }
    }

    return best_eta;
}

/**
 * @brief Result structure for single-feature lcor posterior
 */
struct lcor_posterior_single_t {
    vec_t mean;      // Posterior mean at each vertex
    vec_t sd;        // Posterior SD at each vertex
    vec_t lower;     // Lower credible bound (2.5%)
    vec_t upper;     // Upper credible bound (97.5%)
    double eta;      // Eta used for this feature
    double eff_df;   // Effective degrees of freedom
};

/**
 * @brief Compute lcor posterior summary for a single feature
 *
 * This function generates posterior samples of the smoothed feature z_hat,
 * computes lcor for each sample, and returns summary statistics.
 * Samples are not retained after computation.
 *
 * @param graph Graph structure for lcor computation
 * @param V Eigenvector matrix (n x m)
 * @param eigenvalues Raw eigenvalues (length m)
 * @param y_hat Smoothed response (length n)
 * @param z Original feature values (length n)
 * @param eta Smoothing parameter for this feature
 * @param filter_type Filter type string
 * @param lcor_type Type of lcor weighting
 * @param n_samples Number of posterior samples
 * @param credible_level Credible interval level (e.g., 0.95)
 * @param seed Random seed
 *
 * @return lcor_posterior_single_t with summary statistics
 */
static lcor_posterior_single_t compute_lcor_posterior_single(
    const set_wgraph_t& graph,
    const mat_t& V,
    const vec_t& eigenvalues,
    const vec_t& y_hat,
    const vec_t& z,
    double eta,
    const char* filter_type,
    lcor_type_t lcor_type,
    int n_samples,
    double credible_level,
    unsigned int seed
) {
    const int n = V.rows();
    const int m = V.cols();

    // Compute filtered eigenvalues for this eta
    vec_t f_lambda = compute_filter_weights(eigenvalues, eta, filter_type);

    // Project z onto spectral basis
    vec_t Vt_z = V.transpose() * z;

    // Compute point estimate z_hat
    vec_t alpha_mean = f_lambda.array() * Vt_z.array();
    vec_t z_hat = V * alpha_mean;

    // Estimate residual variance
    vec_t residuals = z - z_hat;
    double eff_df = f_lambda.sum();

    double sigma_hat = 0.0;
    if (n - eff_df > 1.0) {
        sigma_hat = std::sqrt(residuals.squaredNorm() / (n - eff_df));
    } else {
        // Fallback when eff_df is too large
        sigma_hat = std::sqrt(residuals.squaredNorm() / std::max(1.0, n * 0.1));
    }

    // Posterior SD for spectral coefficients
    vec_t alpha_sd(m);
    for (int j = 0; j < m; ++j) {
        double posterior_precision = 1.0 + eta * eigenvalues[j];
        alpha_sd[j] = sigma_hat / std::sqrt(posterior_precision);
    }

    // Storage for lcor samples at each vertex
    mat_t lcor_samples(n, n_samples);

    // Random number generator
    std::mt19937 rng(seed);
    std::normal_distribution<double> std_normal(0.0, 1.0);

    // Convert y_hat to std::vector for lcor interface
    std::vector<double> y_hat_vec(n);
    for (int i = 0; i < n; ++i) {
        y_hat_vec[i] = y_hat[i];
    }

    // Generate samples and compute lcor for each
    for (int b = 0; b < n_samples; ++b) {
        // Draw spectral coefficients from posterior
        vec_t alpha_sample(m);
        for (int j = 0; j < m; ++j) {
            alpha_sample[j] = alpha_mean[j] + alpha_sd[j] * std_normal(rng);
        }

        // Transform to vertex domain
        vec_t z_hat_sample = V * alpha_sample;

        // Convert to std::vector
        std::vector<double> z_hat_vec(n);
        for (int i = 0; i < n; ++i) {
            z_hat_vec[i] = z_hat_sample[i];
        }

        // Compute lcor for this sample using existing set_wgraph_t::lcor()
        std::vector<double> lcor_result = graph.lcor(
            y_hat_vec, z_hat_vec,
            lcor_type,
            edge_diff_type_t::DIFFERENCE,
            edge_diff_type_t::DIFFERENCE,
            0.0,    // epsilon (auto)
            0.0     // winsorize_quantile (none)
        );

        // Store vertex coefficients
        for (int v = 0; v < n; ++v) {
            lcor_samples(v, b) = lcor_result[v];
        }
    }

    // Compute summary statistics at each vertex
    lcor_posterior_single_t result;
    result.mean.resize(n);
    result.sd.resize(n);
    result.lower.resize(n);
    result.upper.resize(n);
    result.eta = eta;
    result.eff_df = eff_df;

    double alpha_lower = (1.0 - credible_level) / 2.0;
    double alpha_upper = (1.0 + credible_level) / 2.0;

    for (int v = 0; v < n; ++v) {
        // Extract samples for this vertex
        std::vector<double> vertex_samples(n_samples);
        for (int b = 0; b < n_samples; ++b) {
            vertex_samples[b] = lcor_samples(v, b);
        }

        // Mean
        double sum = 0.0;
        for (double val : vertex_samples) sum += val;
        result.mean[v] = sum / n_samples;

        // SD
        double sq_sum = 0.0;
        for (double val : vertex_samples) {
            sq_sum += (val - result.mean[v]) * (val - result.mean[v]);
        }
        result.sd[v] = std::sqrt(sq_sum / (n_samples - 1));

        // Quantiles (sort for empirical quantiles)
        std::sort(vertex_samples.begin(), vertex_samples.end());

        int idx_lower = static_cast<int>(std::floor(alpha_lower * (n_samples - 1)));
        int idx_upper = static_cast<int>(std::floor(alpha_upper * (n_samples - 1)));
        idx_lower = std::max(0, std::min(idx_lower, n_samples - 1));
        idx_upper = std::max(0, std::min(idx_upper, n_samples - 1));

        result.lower[v] = vertex_samples[idx_lower];
        result.upper[v] = vertex_samples[idx_upper];
    }

    return result;
}

/**
 * @brief SEXP interface for memory-efficient lcor with posterior propagation
 *
 * @param s_adj_list Adjacency list (0-based from R after conversion)
 * @param s_weight_list Edge weight list
 * @param s_y_hat Smoothed response vector (length n)
 * @param s_Z Original (unsmoothed) feature matrix (n x p)
 * @param s_V Eigenvector matrix from fitted model (n x m)
 * @param s_eigenvalues Raw eigenvalues (length m)
 * @param s_eta_fixed Fixed eta value (used when per_column_gcv = FALSE)
 * @param s_lcor_type "derivative", "unit", or "sign"
 * @param s_filter_type Filter type string
 * @param s_per_column_gcv Logical: select eta per column via GCV
 * @param s_n_gcv_candidates Number of GCV candidates
 * @param s_n_posterior_samples Number of posterior samples
 * @param s_credible_level Credible level (e.g., 0.95)
 * @param s_seed Base random seed
 * @param s_n_cores Number of OpenMP threads
 * @param s_verbose Logical: print progress
 *
 * @return R list with posterior summary matrices
 */
extern "C" SEXP S_lcor_with_posterior_internal(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y_hat,
    SEXP s_Z,
    SEXP s_V,
    SEXP s_eigenvalues,
    SEXP s_eta_fixed,
    SEXP s_lcor_type,
    SEXP s_filter_type,
    SEXP s_per_column_gcv,
    SEXP s_n_gcv_candidates,
    SEXP s_n_posterior_samples,
    SEXP s_credible_level,
    SEXP s_seed,
    SEXP s_n_cores,
    SEXP s_verbose
) {
    // ========================================================================
    // EXTRACT PARAMETERS
    // ========================================================================

    const char* lcor_type_str = CHAR(STRING_ELT(s_lcor_type, 0));
    const char* filter_type_str = CHAR(STRING_ELT(s_filter_type, 0));
    lcor_type_t lcor_type = parse_lcor_type(lcor_type_str);

    bool per_column_gcv = Rf_asLogical(s_per_column_gcv);
    int n_gcv_candidates = Rf_asInteger(s_n_gcv_candidates);
    int n_posterior_samples = Rf_asInteger(s_n_posterior_samples);
    double credible_level = Rf_asReal(s_credible_level);
    unsigned int base_seed = static_cast<unsigned int>(Rf_asInteger(s_seed));
    int n_cores = Rf_asInteger(s_n_cores);
    bool verbose = Rf_asLogical(s_verbose);

    double eta_fixed = Rf_asReal(s_eta_fixed);

    // ========================================================================
    // EXTRACT DIMENSIONS AND BUILD STRUCTURES
    // ========================================================================

    // Get dimensions of Z
    SEXP s_Z_dim = Rf_getAttrib(s_Z, R_DimSymbol);
    const int n = INTEGER(s_Z_dim)[0];
    const int p = INTEGER(s_Z_dim)[1];

    // Get dimensions of V
    SEXP s_V_dim = Rf_getAttrib(s_V, R_DimSymbol);
    const int m = INTEGER(s_V_dim)[1];

    if (verbose) {
        Rprintf("lcor with posterior propagation (C++ implementation)\n");
        Rprintf("  Features: %d, Vertices: %d, Eigenpairs: %d\n", p, n, m);
        Rprintf("  Posterior samples: %d, Per-column GCV: %s\n",
                n_posterior_samples, per_column_gcv ? "yes" : "no");
#ifdef _OPENMP
        Rprintf("  OpenMP threads: %d\n", n_cores);
#else
        Rprintf("  OpenMP: not available (sequential execution)\n");
#endif
    }

    // Build graph
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    // Extract V matrix
    mat_t V(n, m);
    const double* p_V = REAL(s_V);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            V(i, j) = p_V[i + j * n];
        }
    }

    // Extract eigenvalues
    vec_t eigenvalues(m);
    const double* p_eig = REAL(s_eigenvalues);
    for (int j = 0; j < m; ++j) {
        eigenvalues[j] = p_eig[j];
    }

    // Extract y_hat
    vec_t y_hat(n);
    const double* p_y_hat = REAL(s_y_hat);
    for (int i = 0; i < n; ++i) {
        y_hat[i] = p_y_hat[i];
    }

    // Pointer to Z data
    const double* p_Z = REAL(s_Z);

    // ========================================================================
    // ALLOCATE OUTPUT MATRICES
    // ========================================================================

    SEXP s_mean = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    SEXP s_sd = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    SEXP s_lower = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    SEXP s_upper = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    SEXP s_eta_used = PROTECT(Rf_allocVector(REALSXP, p));
    SEXP s_eff_df = PROTECT(Rf_allocVector(REALSXP, p));

    double* out_mean = REAL(s_mean);
    double* out_sd = REAL(s_sd);
    double* out_lower = REAL(s_lower);
    double* out_upper = REAL(s_upper);
    double* out_eta = REAL(s_eta_used);
    double* out_eff_df = REAL(s_eff_df);

    // ========================================================================
    // MAIN LOOP OVER FEATURES (PARALLELIZABLE)
    // ========================================================================

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

    // Note: Progress reporting is tricky with OpenMP, so we skip it in parallel mode
    int progress_interval = std::max(1, p / 20);

    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
    for (int j = 0; j < p; ++j) {

        // Extract column j of Z
        vec_t z_j(n);
        for (int i = 0; i < n; ++i) {
            z_j[i] = p_Z[i + j * n];  // Column-major
        }

        // Determine eta for this feature
        double eta_j;
        if (per_column_gcv) {
            eta_j = select_eta_gcv(V, eigenvalues, z_j, filter_type_str, n_gcv_candidates);
        } else {
            eta_j = eta_fixed;
        }

        // Compute posterior summary for this feature
        // Use different seed for each feature to ensure independence
        unsigned int seed_j = base_seed + static_cast<unsigned int>(j);

        lcor_posterior_single_t result = compute_lcor_posterior_single(
            graph, V, eigenvalues, y_hat, z_j,
            eta_j, filter_type_str, lcor_type,
            n_posterior_samples, credible_level, seed_j
        );

        // Store results (row j in output matrices)
        for (int v = 0; v < n; ++v) {
            out_mean[j + v * p] = result.mean[v];      // Column-major: (j, v)
            out_sd[j + v * p] = result.sd[v];
            out_lower[j + v * p] = result.lower[v];
            out_upper[j + v * p] = result.upper[v];
        }
        out_eta[j] = result.eta;
        out_eff_df[j] = result.eff_df;

        // Progress (only in sequential mode)
        #pragma omp critical
        {
            if (verbose && n_cores == 1 && (j + 1) % progress_interval == 0) {
                Rprintf("  Processed %d/%d features (%.0f%%)\n",
                        j + 1, p, 100.0 * (j + 1) / p);
            }
        }
    }

    if (verbose) {
        Rprintf("  Complete.\n");
    }

    // ========================================================================
    // BUILD RESULT LIST
    // ========================================================================

    const int n_components = 8;
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    SET_STRING_ELT(s_names, idx, Rf_mkChar("mean"));
    SET_VECTOR_ELT(s_result, idx++, s_mean);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("sd"));
    SET_VECTOR_ELT(s_result, idx++, s_sd);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("lower"));
    SET_VECTOR_ELT(s_result, idx++, s_lower);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("upper"));
    SET_VECTOR_ELT(s_result, idx++, s_upper);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("eta.used"));
    SET_VECTOR_ELT(s_result, idx++, s_eta_used);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("effective.df"));
    SET_VECTOR_ELT(s_result, idx++, s_eff_df);

    SET_STRING_ELT(s_names, idx, Rf_mkChar("n.samples"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarInteger(n_posterior_samples));

    SET_STRING_ELT(s_names, idx, Rf_mkChar("credible.level"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarReal(credible_level));

    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    UNPROTECT(8);  // s_mean, s_sd, s_lower, s_upper, s_eta_used, s_eff_df, s_result, s_names

    return s_result;
}
