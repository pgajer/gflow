#ifdef _OPENMP
    #include <omp.h>
#endif

// Undefine R's match macro if it exists
#ifdef match
    #undef match
#endif

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <execution>
#include <atomic>
#include <mutex>
#include <vector>
#include <algorithm> // for std::max
#include <random>
#include <cmath>         // for fabs()
#include <unordered_set>
#include <unordered_map>

#include "pgmalog.hpp"
#include "sampling.h" // for C_runif_simplex()
#include "error_utils.h"
#include "pglm.h"
#include "msr2.h"
#include "path_graphs.hpp"
#include "cpp_utils.hpp"                 // for elapsed_time
#include "SEXP_cpp_conversion_utils.hpp"
#include "kernels.h"
#include "1D_linear_models.h"
#include "predictive_errors.hpp"

/**
 * @brief Performs parallel bootstrap calculations for Graph Path Linear Model using STL parallel algorithms
 *
 * @param path_graph The path graph structure containing adjacency lists
 * @param y Vector of observed values
 * @param n_bb Number of bootstrap iterations
 * @param max_distance_deviation Maximum allowed deviation in distance calculations
 * @param ikernel Kernel function identifier
 * @param dist_normalization_factor Factor for normalizing distances
 * @param epsilon Numerical tolerance parameter
 *
 * @return Vector of vectors containing bootstrap predictions for each iteration
 */
std::vector<std::vector<double>> pgmalog_bb(const path_graph_plm_t& path_graph,
                                         const std::vector<double>& y,
                                         int n_bb,
                                         int max_distance_deviation,
                                         int ikernel,
                                         double dist_normalization_factor = 1.01,
                                         double epsilon = 1e-8) {

    int n_vertices = static_cast<int>(y.size());

    // Initialize results vector
    std::vector<std::vector<double>> bb_Ey(n_bb);
    for (auto& Ey : bb_Ey) {
        Ey.resize(n_vertices);
    }

    // Create indices for parallel iteration
    std::vector<int> bb_indices(n_bb);
    std::iota(bb_indices.begin(), bb_indices.end(), 0);

    // Mutex for thread-safe random number generation
    std::mutex rng_mutex;

    // Parallel execution of bootstrap iterations
    std::for_each(std::execution::par_unseq,
                  bb_indices.begin(),
                  bb_indices.end(),
                  [&](int iboot) {
        // Thread-local weight vector
        std::vector<double> weights(n_vertices);

        // Generate weights in a thread-safe manner
        {
            std::lock_guard<std::mutex> lock(rng_mutex);
            C_runif_simplex(&n_vertices, weights.data());
        }

        // Compute predictions for this bootstrap iteration
        auto [predictions, errors] = pgmalog(path_graph,
                                            y,
                                            weights,
                                            ikernel,
                                            max_distance_deviation,
                                            dist_normalization_factor,
                                            epsilon);

        // Store results - no need for mutex as each thread writes to its own index
        bb_Ey[iboot] = std::move(predictions);
    });

    return bb_Ey;
}

/**
 * @brief Performs Bayesian bootstrap estimation with credible intervals for graph path linear models
 *
 * @details This function implements a complete Bayesian bootstrap analysis in three steps:
 * 1. Performs multiple iterations of weighted graph path linear model estimation
 * 2. Computes central location estimates (mean or median) for each vertex
 * 3. Calculates credible intervals at the specified probability level
 *
 * The implementation supports parallel processing for large graphs and provides
 * robust uncertainty quantification through bootstrap resampling.
 *
 * @param path_graph PLM graph structure containing adjacency lists and path information
 * @param y Response variable values at each vertex
 * @param p Probability level for credible intervals (default: 0.95)
 * @param n_bb Number of bootstrap iterations (default: 100)
 * @param max_distance_deviation Maximum allowed deviation from optimal center position
 * @param use_median If true, uses median instead of mean for central location (default: false)
 * @param ikernel Kernel function selector (default: 1)
 * @param dist_normalization_factor Distance normalization factor (default: 1.01)
 * @param epsilon Numerical stability parameter (default: 1e-15)
 *
 * @return bb_cri_t struct containing:
 *         - bb_Ey: Central location estimates (mean/median) for each vertex
 *         - cri_L: Lower bounds of credible intervals
 *         - cri_U: Upper bounds of credible intervals
 *
 * @throws Rf_error if:
 *         - Input dimensions are inconsistent
 *         - Parameters are invalid (p ∉ (0,1), n_bb ≤ 0, n_cores ≤ 0)
 *         - Numerical instability occurs during computation
 *         - Memory allocation fails
 */
bb_cri_t pgmalog_bb_cri(const path_graph_plm_t& path_graph,
                     const std::vector<double>& y,
                     double p,
                     int n_bb,
                     int max_distance_deviation,
                     bool use_median,
                     int ikernel,
                     double dist_normalization_factor,
                     double epsilon) {

    // Perform bootstrap iterations
    std::vector<std::vector<double>> bb_Eys = pgmalog_bb(path_graph,
                                                      y,
                                                      n_bb,
                                                      max_distance_deviation,
                                                      ikernel,
                                                      dist_normalization_factor,
                                                      epsilon);

    // Calculate credible intervals
    return bb_cri(bb_Eys, use_median, p);
}
