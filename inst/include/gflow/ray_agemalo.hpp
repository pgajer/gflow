#ifndef RAY_AGEMALO_H_
#define RAY_AGEMALO_H_

#include "uniform_grid_graph.hpp"
#include "wasserstein_perm_test.hpp"
#include "agemalo.hpp"

agemalo_result_t ray_agemalo(
    const uniform_grid_graph_t& grid_graph,
    const std::vector<double>& y,
    // geodesic parameter
    size_t min_path_size,
    // bw parameters
    size_t n_bws,
    bool log_grid,
    double min_bw_factor,
    double max_bw_factor,
    // kernel parameters
    double dist_normalization_factor,
    size_t kernel_type,
    // model parameters
    double model_tolerance,
    double model_blending_coef,
    // Bayesian bootstrap parameters
    size_t n_bb,
    double cri_probability,
    // permutation parameters
    size_t n_perms,
    // verbose
    bool verbose
    );

#endif // RAY_AGEMALO_H_
