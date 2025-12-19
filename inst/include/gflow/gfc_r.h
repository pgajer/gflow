/**
 * @file gfc_r.h
 * @brief SEXP interface declarations for Gradient Flow Complex computation
 *
 * This header declares the R-callable C functions that bridge the gap
 * between R and the C++ GFC implementation. These functions are registered
 * in init.c for use with .Call().
 */

#ifndef GFC_R_H
#define GFC_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute refined gradient flow complex for a single function
 *
 * @param s_adj_list Adjacency list (0-based after R conversion)
 * @param s_weight_list Edge weight list
 * @param s_y Function values (length n)
 * @param s_edge_length_quantile_thld Edge length quantile threshold
 * @param s_min_rel_value_max Minimum relative value for maxima
 * @param s_max_rel_value_min Maximum relative value for minima
 * @param s_max_overlap_threshold Overlap threshold for maxima clustering
 * @param s_min_overlap_threshold Overlap threshold for minima clustering
 * @param s_p_mean_nbrs_dist_threshold Mean neighbor distance percentile threshold (maxima only)
 * @param s_p_mean_hopk_dist_threshold Hop-k distance percentile threshold
 * @param s_p_deg_threshold Degree percentile threshold
 * @param s_min_basin_size Minimum basin size
 * @param s_expand_basins Whether to expand basins
 * @param s_apply_relvalue_filter Whether to apply relative value filtering
 * @param s_apply_maxima_clustering Whether to cluster maxima
 * @param s_apply_minima_clustering Whether to cluster minima
 * @param s_apply_geometric_filter Whether to apply geometric filtering
 * @param s_hop_k Hop distance for summary statistics
 * @param s_verbose Print progress messages
 *
 * @return R list containing GFC results
 */
SEXP S_compute_gfc(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_min_rel_value_max,
    SEXP s_max_rel_value_min,
    SEXP s_max_overlap_threshold,
    SEXP s_min_overlap_threshold,
    SEXP s_p_mean_nbrs_dist_threshold,
    SEXP s_p_mean_hopk_dist_threshold,
    SEXP s_p_deg_threshold,
    SEXP s_min_basin_size,
    SEXP s_expand_basins,
    SEXP s_apply_relvalue_filter,
    SEXP s_apply_maxima_clustering,
    SEXP s_apply_minima_clustering,
    SEXP s_apply_geometric_filter,
    SEXP s_hop_k,
	SEXP s_with_trajectories,
    SEXP s_verbose
);

/**
 * @brief Compute GFC for multiple functions over the same graph
 *
 * @param s_adj_list Adjacency list (0-based after R conversion)
 * @param s_weight_list Edge weight list
 * @param s_Y Matrix of function values (n x p)
 * @param s_edge_length_quantile_thld Edge length quantile threshold
 * @param s_min_rel_value_max Minimum relative value for maxima
 * @param s_max_rel_value_min Maximum relative value for minima
 * @param s_max_overlap_threshold Overlap threshold for maxima clustering
 * @param s_min_overlap_threshold Overlap threshold for minima clustering
 * @param s_p_mean_nbrs_dist_threshold Mean neighbor distance percentile threshold
 * @param s_p_mean_hopk_dist_threshold Hop-k distance percentile threshold
 * @param s_p_deg_threshold Degree percentile threshold
 * @param s_min_basin_size Minimum basin size
 * @param s_expand_basins Whether to expand basins
 * @param s_apply_relvalue_filter Whether to apply relative value filtering
 * @param s_apply_maxima_clustering Whether to cluster maxima
 * @param s_apply_minima_clustering Whether to cluster minima
 * @param s_apply_geometric_filter Whether to apply geometric filtering
 * @param s_hop_k Hop distance for summary statistics
 * @param s_n_cores Number of OpenMP threads
 * @param s_verbose Print progress messages
 *
 * @return R list containing vector of GFC results
 */
SEXP S_compute_gfc_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_min_rel_value_max,
    SEXP s_max_rel_value_min,
    SEXP s_max_overlap_threshold,
    SEXP s_min_overlap_threshold,
    SEXP s_p_mean_nbrs_dist_threshold,
    SEXP s_p_mean_hopk_dist_threshold,
    SEXP s_p_deg_threshold,
    SEXP s_min_basin_size,
    SEXP s_expand_basins,
    SEXP s_apply_relvalue_filter,
    SEXP s_apply_maxima_clustering,
    SEXP s_apply_minima_clustering,
    SEXP s_apply_geometric_filter,
    SEXP s_hop_k,
    SEXP s_n_cores,
    SEXP s_verbose
);

#ifdef __cplusplus
}
#endif

#endif // GFC_R_H
