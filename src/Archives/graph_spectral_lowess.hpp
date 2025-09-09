#ifndef GRAPH_SPECTRAL_LOWESS_HPP
#define GRAPH_SPECTRAL_LOWESS_HPP

#include <vector>
#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

/**
 * @brief Result structure for graph-based spectral LOWESS smoothing
 */
struct graph_spectral_lowess_t {
    std::vector<double> predictions;  ///< Smoothed predictions for each vertex
    std::vector<double> errors;       ///< Estimated prediction errors for each vertex
    std::vector<double> scale;        ///< Local scale/bandwidth for each vertex
};

/**
 * @brief Extension of set_wgraph_t with spectral LOWESS methods
 *
 * @note These methods would typically be included in set_wgraph.hpp
 */
class set_wgraph_t {
public:
    /**
     * @brief Original spectral LOWESS implementation (no model averaging)
     */
    graph_spectral_lowess_t graph_spectral_lowess(
        const std::vector<double>& y,
        size_t n_evectors,
        // bw parameters
        size_t n_bws,
        bool log_grid,
        double min_bw_factor,
        double max_bw_factor,
        // kernel parameters
        double dist_normalization_factor,
        size_t kernel_type,
        // other
        double precision,
        bool verbose
    ) const;

    /**
     * @brief Model-averaged spectral LOWESS implementation
     */
    graph_spectral_lowess_t graph_spectral_ma_lowess(
        const std::vector<double>& y,
        size_t n_evectors,
        // bw parameters
        size_t n_bws,
        bool log_grid,
        double min_bw_factor,
        double max_bw_factor,
        // kernel parameters
        double dist_normalization_factor,
        size_t kernel_type,
        // model parameters
        double model_blending_coef,
        // other
        double precision,
        bool verbose
    ) const;

    /**
     * @brief Find vertices within a specified radius of a reference vertex
     */
    std::unordered_map<size_t, double> find_vertices_within_radius(
        size_t vertex,
        double radius
    ) const;

    /**
     * @brief Find minimum radius that includes a specified number of vertices
     */
    double find_minimum_radius_for_domain_min_size(
        size_t vertex,
        double lower_bound,
        double upper_bound,
        size_t domain_min_size,
        double precision
    ) const;
};

/**
 * @brief Creates a spectral embedding for a set of vertices using Laplacian eigenvectors
 */
Eigen::MatrixXd create_spectral_embedding(
    const std::map<size_t, double>& vertex_map,
    const Eigen::MatrixXd& eigenvectors,
    size_t n_evectors
);

/**
 * @brief Fits a weighted linear model in the spectral embedding space
 */
lm_t fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    size_t kernel_type,
    double dist_normalization_factor
);

#endif // GRAPH_SPECTRAL_LOWESS_HPP
