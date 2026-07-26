#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rdefines.h>
#include <stdlib.h>

#include "cpp_mstrees_r.h"
#include "cpp_stats_utils_r.h"
#include "fns_over_graphs_r.h"
#include "graph_conn_components_r.h"
#include "graph_cycles_r.h"
#include "graph_edit_distance_r.h"
#include "graph_ms_cx_r.h"
#include "graph_shortest_path_r.h"
#include "graph_spectrum_r.h"
#include "graph_utils_r.h"
#include "grids.h"
#include "lm.h"
#include "mstree.h"
#include "mstree_total_length_r.h"
#include "pruning_long_edges_r.h"
#include "random_sampling_r.h"
#include "sampling.h"  // For C_runif_simplex()
#include "centered_paths_r.h"
#include "iknn_graphs_r.h"
#include "graph_gradient_flow_r.h"
#include "kNN_r.h"       // for S_kNN()
#include "set_wgraph_r.h"
#include "parameterize_circular_graph_r.h"
#include "geodesic_stats_r.h"
#include "monotonic_reachability_r.h"
#include "local_extrema_r.h"
#include "hop_extremp_radii_r.h"
#include "gflow_basins_r.h"
#include "gflow_cx_r.h"
#include "fn_graphs_r.h"
#include "nerve_cx_r.h"
#include "stats_utils.h"
#include "kernels.h"
#include "riem_dcx_r.h"
#include "lcor_r.h"
#include "partition_graph_r.h"
#include "lslope_r.h"
#include "gfc_r.h"
#include "basin_summary_r.h"
#include "gfassoc_r.h"
#include "harmonic_extension_r.h"
#include "gfc_flow_r.h"
#include "madag_r.h"
#include "traj_clustering_r.h"
#include "graph_core_endpoints_r.h"
#include "linf_simplex_knn_r.h"
#include "local_geodesic_pruning_r.h"

static R_NativePrimitiveArgType create_ED_grid_2D_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType create_ED_grid_3D_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType create_ED_grid_xD_type[] = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType create_ENPs_grid_2D_type[] = {INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType create_ENPs_grid_3D_type[] = {INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType mstree_type[] = {INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType runif_simplex_type[] = {INTSXP, REALSXP};

static R_NativePrimitiveArgType matrix_wmeans_type[]            = {REALSXP, INTSXP, INTSXP,  INTSXP,  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_type[]         = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_type[]      = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_qCrI_type[] = {INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_CrI_type[]  = {REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType normalize_dist_type[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType pearson_wcor_type[] = {REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType pearson_wcor_BB_qCrI_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType kernel_eval_type[] = { INTSXP, REALSXP, INTSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType kernel_type[] = { REALSXP, INTSXP, REALSXP, INTSXP, REALSXP };

static const R_CMethodDef cMethods[] = {
  {"C_epanechnikov_kernel_with_stop", (DL_FUNC) &C_epanechnikov_kernel_with_stop, 5, kernel_type},
  {"C_triangular_kernel_with_stop", (DL_FUNC) &C_triangular_kernel_with_stop, 5, kernel_type},
  {"C_tr_exponential_kernel_with_stop", (DL_FUNC) &C_tr_exponential_kernel_with_stop, 5, kernel_type},
  {"C_kernel_eval", (DL_FUNC) &C_kernel_eval, 5, kernel_eval_type},
  {"C_runif_simplex", (DL_FUNC) &C_runif_simplex, 2, runif_simplex_type},
  {"C_create_ED_grid_2D", (DL_FUNC) &C_create_ED_grid_2D, 6, create_ED_grid_2D_type},
  {"C_create_ED_grid_3D", (DL_FUNC) &C_create_ED_grid_3D, 8, create_ED_grid_3D_type},
  {"C_create_ED_grid_xD", (DL_FUNC) &C_create_ED_grid_xD, 6, create_ED_grid_xD_type},
  {"C_create_ENPs_grid_2D", (DL_FUNC) &C_create_ENPs_grid_2D, 7, create_ENPs_grid_2D_type},
  {"C_create_ENPs_grid_3D", (DL_FUNC) &C_create_ENPs_grid_3D, 9, create_ENPs_grid_3D_type},
  {"C_mstree", (DL_FUNC) &C_mstree, 7, mstree_type},
  {"C_matrix_wmeans", (DL_FUNC) &C_matrix_wmeans, 9, matrix_wmeans_type},
  {"C_columnwise_wmean", (DL_FUNC) &C_columnwise_wmean, 6, columnwise_wmean_type},
  {"C_columnwise_wmean_BB", (DL_FUNC) &C_columnwise_wmean_BB, 7, columnwise_wmean_BB_type},
  {"C_columnwise_wmean_BB_qCrI", (DL_FUNC) &C_columnwise_wmean_BB_qCrI, 9, columnwise_wmean_BB_qCrI_type},
  {"C_columnwise_wmean_BB_CrI_1", (DL_FUNC) &C_columnwise_wmean_BB_CrI_1, 8, columnwise_wmean_BB_CrI_type},
  {"C_columnwise_wmean_BB_CrI_2", (DL_FUNC) &C_columnwise_wmean_BB_CrI_2, 8, columnwise_wmean_BB_CrI_type},
  {"C_normalize_dist", (DL_FUNC) &C_normalize_dist, 7, normalize_dist_type},
  {"C_pearson_wcor", (DL_FUNC) &C_pearson_wcor, 5, pearson_wcor_type},
  {"C_pearson_wcor_BB_qCrI", (DL_FUNC) &C_pearson_wcor_BB_qCrI, 10, pearson_wcor_BB_qCrI_type},
  {NULL, NULL, 0, NULL}
};

#ifdef __cplusplus
extern "C" {
#endif
SEXP _gflow_rcpp_local_pca_chart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
#ifdef __cplusplus
}
#endif

SEXP S_gflow_openmp_diag();

static const R_CallMethodDef CallMethods[] = {
  // diagnostic
  {"S_gflow_openmp_diag", (DL_FUNC) &S_gflow_openmp_diag, 0},

  // clustering cell trajectories
  {"S_cluster_cell_trajectories", (DL_FUNC) &S_cluster_cell_trajectories, 12},

  // =========================================================================
  //  Harmonic extension
  // =========================================================================
  {"S_compute_harmonic_extension", (DL_FUNC) &S_compute_harmonic_extension, 10},
  {"S_select_max_density_trajectory", (DL_FUNC) &S_select_max_density_trajectory, 2},

  // =========================================================================
  //  Gradient flow association functions
  // =========================================================================
  {"S_gfassoc_membership",      (DL_FUNC) &S_gfassoc_membership,      3},
  {"S_gfassoc_polarity",        (DL_FUNC) &S_gfassoc_polarity,        4},
  {"S_gfassoc_basin_character", (DL_FUNC) &S_gfassoc_basin_character, 5},
  {"S_gfassoc_overlap",         (DL_FUNC) &S_gfassoc_overlap,         3},
  {"S_gfassoc_deviation",       (DL_FUNC) &S_gfassoc_deviation,       1},
  {"S_gfcor",                   (DL_FUNC) &S_gfcor,                   8},

  // =========================================================================
  // gradient flow
  // =========================================================================
  {"S_construct_madag", (DL_FUNC) &S_construct_madag, 6},
  {"S_enumerate_cell_trajectories", (DL_FUNC) &S_enumerate_cell_trajectories, 4},
  {"S_sample_cell_trajectories", (DL_FUNC) &S_sample_cell_trajectories, 5},
  {"S_trajectory_similarity_matrix", (DL_FUNC) &S_trajectory_similarity_matrix, 2},
  {"S_identify_bottlenecks", (DL_FUNC) &S_identify_bottlenecks, 3},

  {"S_compute_gfc_flow",        (DL_FUNC) &S_compute_gfc_flow,        6},
  {"S_compute_gfc_flow_matrix", (DL_FUNC) &S_compute_gfc_flow_matrix, 7},
  {"S_compute_vertex_density",  (DL_FUNC) &S_compute_vertex_density, 3},
  {"S_compute_gfc_basins",      (DL_FUNC) &S_compute_gfc_basins, 7},
  {"S_compute_basins_of_attraction", (DL_FUNC) &S_compute_basins_of_attraction, 5},
  {"S_compute_basins_of_attraction_rtcb", (DL_FUNC) &S_compute_basins_of_attraction_rtcb, 6},
  {"S_precompute_basin_vertex_metrics", (DL_FUNC) &S_precompute_basin_vertex_metrics, 3},
  {"S_summary_basins_of_attraction_cpp", (DL_FUNC) &S_summary_basins_of_attraction_cpp, 5},
  {"S_compute_hop_extremp_radii_batch", (DL_FUNC) &S_compute_hop_extremp_radii_batch, 8},
  {"S_partition_graph",        (DL_FUNC) &S_partition_graph, 4},

  {"S_create_basin_cx", (DL_FUNC) &S_create_basin_cx, 3},
  {"S_find_gflow_basins", (DL_FUNC) &S_find_gflow_basins, 6},
  {"S_construct_graph_gradient_flow", (DL_FUNC) &S_construct_graph_gradient_flow, 6},

  {"S_analyze_function_aware_weights", (DL_FUNC) &S_analyze_function_aware_weights, 12},
  {"S_create_gflow_cx", (DL_FUNC) &S_create_gflow_cx, 12},
  {"S_compute_extrema_hop_nbhds", (DL_FUNC) &S_compute_extrema_hop_nbhds, 3},
  {"S_find_local_extrema", (DL_FUNC) &S_find_local_extrema, 4},
  {"S_detect_local_extrema", (DL_FUNC) &S_detect_local_extrema, 6},
  {"S_graph_MS_cx_with_path_search", (DL_FUNC) &S_graph_MS_cx_with_path_search, 3},
  {"S_graph_MS_cx_using_short_h_hops", (DL_FUNC) &S_graph_MS_cx_using_short_h_hops, 4},
  {"S_graph_constrained_gradient_flow_trajectories", (DL_FUNC) &S_graph_constrained_gradient_flow_trajectories, 3},

  // =========================================================================
  // lslope
  // =========================================================================
  {"S_lslope_gradient_instrumented", (DL_FUNC) &S_lslope_gradient_instrumented, 11},
  {"S_lslope_gradient", (DL_FUNC) &S_lslope_gradient, 11},
  {"S_lslope_neighborhood", (DL_FUNC) &S_lslope_neighborhood, 8},
  {"S_lslope_vector_matrix", (DL_FUNC) &S_lslope_vector_matrix, 11},

  // =========================================================================
  // lcor
  // =========================================================================
  {"S_lcor", (DL_FUNC) &S_lcor, 9},
  {"S_lcor_instrumented", (DL_FUNC) &S_lcor_instrumented, 9},
  {"S_lcor_vector_matrix", (DL_FUNC) &S_lcor_vector_matrix, 10},

  // =========================================================================
  // data graphs
  // =========================================================================
  {"S_prune_graph_local_geodesic", (DL_FUNC) &S_prune_graph_local_geodesic, 6},
  {"S_prune_graph_global_geodesic_ratio", (DL_FUNC) &S_prune_graph_global_geodesic_ratio, 5},
  {"S_construct_function_aware_graph", (DL_FUNC) &S_construct_function_aware_graph, 14},
  {"S_extract_skeleton_graph", (DL_FUNC) &S_extract_skeleton_graph, 1},
  {"S_get_simplex_counts", (DL_FUNC) &S_get_simplex_counts, 1},

  {"S_kNN", (DL_FUNC) &S_kNN, 2},

  {"S_create_nerve_complex", (DL_FUNC) &S_create_nerve_complex, 3},

  // =========================================================================
  // graph utilities
  // =========================================================================
  {"S_detect_major_arms", (DL_FUNC) &S_detect_major_arms, 11},
  {"S_solve_full_laplacian", (DL_FUNC) &S_solve_full_laplacian, 3},
  {"S_parameterize_circular_graph", (DL_FUNC) &S_parameterize_circular_graph, 3},
  {"S_compute_tube_lens_corridor", (DL_FUNC) &S_compute_tube_lens_corridor, 7},
  {"S_get_path_data", (DL_FUNC) &S_get_path_data, 10},
  {"S_ugg_get_path_data", (DL_FUNC) &S_ugg_get_path_data, 11},

  {"_gflow_rcpp_local_pca_chart", (DL_FUNC) &_gflow_rcpp_local_pca_chart, 9},

  // =========================================================================
  // stat utilities
  // =========================================================================
  {"S_ecdf", (DL_FUNC) &S_ecdf, 1},
  {"S_rlaplace", (DL_FUNC) &S_rlaplace, 4},

  // =========================================================================
  // other
  // =========================================================================
  {"S_set_function_values", (DL_FUNC) &S_set_function_values, 2},
  {"S_set_weight_scheme", (DL_FUNC) &S_set_weight_scheme, 3},

  // =========================================================================
  // to remove
  // =========================================================================
  //{"S_kNN_v2", (DL_FUNC) &S_kNN_v2, 2},
  //
  {NULL, NULL, 0}
};

void R_init_gflow(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, CallMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);     // it will allow Rcpp’s generated .Call("...") work
}
