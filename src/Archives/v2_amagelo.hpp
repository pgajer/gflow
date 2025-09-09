#ifndef AMAGELO_HPP
#define AMAGELO_HPP

#include <vector>

struct extremum_t {
    size_t idx;                 ///< Index of extremum in prediction vector
    double x;                   ///< x-coordinate
    double y;                   ///< predicted ŷ value at extremum
    bool is_max;                ///< true = max, false = min
    double depth;               ///< vertical prominence
    size_t depth_idx;           ///< index of point where min/max descent terminates

    double rel_depth;           ///< depth / total depth
    double range_rel_depth;     ///< depth / (max(ŷ) - min(ŷ))
    // double q_rel_depth;         ///< quantile of rel_depth
    // double q_range_rel_depth;   ///< quantile of range_rel_depth

    extremum_t(size_t idx_,
               double x_,
               double y_,
               bool is_max_,
               double depth_     = 0.0,
               size_t depth_idx_ = 0,
               double rel_depth_ = 0.0,
               double range_rel_depth_ = 0.0)
        : idx(idx_), x(x_), y(y_), is_max(is_max_),
          depth(depth_), depth_idx(depth_idx_),
          rel_depth(rel_depth_),
          range_rel_depth(range_rel_depth_) {}
};


/**
 * @brief Return structure for amagelo smoother
 */
struct amagelo_t {
    std::vector<double> x_sorted;
    std::vector<double> y_sorted;
    std::vector<std::size_t> order;
    std::vector<double> grid_coords;

    std::vector<double> predictions;                 ///< predictions[i] is an estimate of E(Y|G) at the i-th vertex at the optimal bw, so predictions = bw_predictions[opt_bw_idx]
    std::vector<std::vector<double>> bw_predictions; ///< bw_predictions[bw_idx] are predictions for the bandwidth index 'bw_idx'
    std::vector<double> grid_predictions;            ///< predictions of y values over the grid vertices
    std::vector<double> harmonic_predictions;        ///< predictions after triplet harmonic smoothing

    std::vector<extremum_t> local_extrema;

    double tvmi;                                         ///< Total‐Variation Monotonicity Index (net change / total abs‐variation)
    std::vector<double> monotonic_interval_proportions;  ///< the p_i' - nonotnicity interval relative lengths - this verctor allows us to see where on the domain the smoother changes direction, and by how much.
    double simpson_index;                                ///< ∑ p_i^2 Simpson index - (sometimes called the “concentration” or “dominance” index) is highest (=1) when one interval dominates (very monotonic) and lowest when all intervals are equal (maximally wiggly).

    std::vector<double> bw_errors;                   ///< bw_errors[bw_idx] is the mean LOOCV error estimate over the models of the bw_idx-th bandwidth
    size_t opt_bw_idx;                               ///< the index of opt_bw in bws - the index is the same for all vertices
    double min_bw;
    double max_bw;
    std::vector<double> bws;
};

amagelo_t amagelo(
    const std::vector<double>& x,
    const std::vector<double>& y,
    size_t grid_size,
    double min_bw_factor,
    double max_bw_factor,
    size_t n_bws,
    bool use_global_bw_grid,
    bool with_bw_predictions,
    bool log_grid,
    size_t domain_min_size,
    size_t kernel_type,
    double dist_normalization_factor,
    size_t n_cleveland_iterations,
    double blending_coef,
    bool use_linear_blending,
    double precision,
    bool verbose
    );

#endif // AMAGELO_HPP
