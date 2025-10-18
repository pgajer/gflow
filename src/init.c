#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rdefines.h>
#include <stdlib.h>

#include "msr2.h"
#include "sampling.h"  // For C_runif_simplex()
#include "density_r.h"
#include "local_complexity_r.h"
#include "ulogit_r.h"
#include "maelog_r.h"
#include "adaptive_maelog_r.h"
#include "magelog_r.h"
#include "pgmalo_r.h"
#include "pgmalog_r.h"
#include "uniform_grid_graph_r.h"
#include "centered_paths_r.h"
#include "uggmalog_r.h"
#include "uggmalo_r.h"
#include "iknn_graphs_r.h"
#include "mknn_graphs_r.h"
#include "graph_gradient_flow_r.h"
#include "kNN_r.h"       // for S_kNN()
#include "wasserstein_dist.h" // for C_wasserstein_distance_1D()
#include "adaptive_uggmalo_r.h"
#include "set_wgraph_r.h"
#include "parameterize_circular_graph_r.h"
#include "graph_diffusion_smoother_r.h"
#include "graph_maximal_packing_r.h"
#include "agemalo_r.h"
#include "ray_agemalo_r.h"
#include "geodesic_stats_r.h"
#include "graph_spectral_lowess_r.h"
#include "graph_spectral_ma_lowess_r.h"
#include "graph_spectral_lowess_mat_r.h"
#include "spectral_lowess_graph_smoothing_r.h"
#include "monotonic_reachability_r.h"
#include "local_extrema_r.h"
#include "nada_graph_spectral_lowess_r.h"
#include "graph_deg0_lowess_r.h"
#include "graph_deg0_lowess_cv_r.h"
#include "graph_deg0_lowess_cv_mat_r.h"
#include "deg0_lowess_graph_smoothing_r.h"
#include "graph_deg0_lowess_buffer_cv_r.h"
#include "graph_kernel_smoother_r.h"
#include "graph_bw_adaptive_spectral_smoother_r.h"
#include "klaps_low_pass_smoother_r.h"
#include "graph_spectral_filter_r.h"
#include "mst_completion_graphs_r.h"
#include "amagelo_r.h"
#include "gflow_basins_r.h"
#include "harmonic_smoother_r.h"
#include "gflow_cx_r.h"
#include "fn_graphs_r.h"
#include "nerve_cx_r.h"
#include "stats_utils.h"
#include "kernels.h"
#include "mean_shift_smoother_r.h"
#include "cv_deg0.h"
#include "riem_dcx_r.h"

static R_NativePrimitiveArgType create_ED_grid_2D_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType create_ED_grid_3D_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType create_ED_grid_xD_type[] = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType create_ENPs_grid_2D_type[] = {INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType create_ENPs_grid_3D_type[] = {INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType mstree_type[] = {INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType wasserstein_distance_1D_type[] = {REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType runif_simplex_type[] = {INTSXP, REALSXP};

static R_NativePrimitiveArgType llm_1D_beta_type[] = {REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType llm_1D_beta_perms_type[] = {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType predict_1D_type[] = {REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType wpredict_1D_type[] = {REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType llm_1D_fit_and_predict_type[] = {INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType llm_1D_fit_and_predict_BB_CrI_type[]  = {INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType llm_1D_fit_and_predict_BB_qCrI_type[] = {INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType llm_1D_fit_and_predict_BB_type[]     = {INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mllm_1D_fit_and_predict_type[] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType loo_llm_1D_type[] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType deg0_loo_llm_1D_type[] = {INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType cv_mae_1D_type[] = {INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType get_BB_Eyg_type[] = {INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType get_Eyg_CrI_type[] = {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType get_Eygs_type[] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType get_bws_type[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_weighting_type[]    = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_eval_type[]    = {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_TS_norm_type[] = {REALSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType matrix_wmeans_type[]            = {REALSXP, INTSXP, INTSXP,  INTSXP,  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_type[]         = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_type[]      = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_qCrI_type[] = {INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType columnwise_wmean_BB_CrI_type[]  = {REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType quantiles_type[] = {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType modified_columnwise_wmean_BB_type[] = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mat_columnwise_divide_type[] = {REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType normalize_dist_type[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType normalize_dist_with_minK_a_type[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType samplewr_type[] = {INTSXP, INTSXP};
static R_NativePrimitiveArgType vpermute_type[] = {INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType v_get_folds_type[] = {INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType winsorize_type[] = {REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType pdistr_type[] = {REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType pearson_cor_type[] = {REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType wcov_type[] = {REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType pearson_wcor_type[] = {REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType pearson_wcor_BB_qCrI_type[] = {REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType density_distance_type[] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType rmatrix_type[] = {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType kernel_eval_type[] = { INTSXP, REALSXP, INTSXP, REALSXP, REALSXP };
static R_NativePrimitiveArgType kernel_type[] = { REALSXP, INTSXP, REALSXP, INTSXP, REALSXP };

static R_NativePrimitiveArgType cv_type[] = { INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};

static const R_CMethodDef cMethods[] = {
  {"C_cv_deg0_binloss", (DL_FUNC) &C_cv_deg0_binloss, 14, cv_type},
  {"C_cv_deg0_mae", (DL_FUNC) &C_cv_deg0_mae, 14, cv_type},
  {"C_epanechnikov_kernel_with_stop", (DL_FUNC) &C_epanechnikov_kernel_with_stop, 5, kernel_type},
  {"C_triangular_kernel_with_stop", (DL_FUNC) &C_triangular_kernel_with_stop, 5, kernel_type},
  {"C_tr_exponential_kernel_with_stop", (DL_FUNC) &C_tr_exponential_kernel_with_stop, 5, kernel_type},
  {"C_kernel_eval", (DL_FUNC) &C_kernel_eval, 5, kernel_eval_type},
  {"C_llm_1D_beta", (DL_FUNC) &C_llm_1D_beta, 8, llm_1D_beta_type},
  {"C_llm_1D_beta_perms", (DL_FUNC) &C_llm_1D_beta_perms, 11, llm_1D_beta_perms_type},
  {"C_predict_1D", (DL_FUNC) &C_predict_1D, 11, predict_1D_type},
  {"C_wpredict_1D", (DL_FUNC) &C_wpredict_1D, 11, wpredict_1D_type},
  {"C_llm_1D_fit_and_predict", (DL_FUNC) &C_llm_1D_fit_and_predict, 12, llm_1D_fit_and_predict_type},
  {"C_mllm_1D_fit_and_predict", (DL_FUNC) &C_mllm_1D_fit_and_predict, 12, mllm_1D_fit_and_predict_type},
  {"C_llm_1D_fit_and_predict_BB_CrI", (DL_FUNC) &C_llm_1D_fit_and_predict_BB_CrI, 13, llm_1D_fit_and_predict_BB_CrI_type},
  {"C_llm_1D_fit_and_predict_global_BB_CrI", (DL_FUNC) &C_llm_1D_fit_and_predict_global_BB_CrI, 13, llm_1D_fit_and_predict_BB_CrI_type},
  {"C_llm_1D_fit_and_predict_global_BB", (DL_FUNC) &C_llm_1D_fit_and_predict_global_BB, 11, llm_1D_fit_and_predict_BB_type},
  {"C_llm_1D_fit_and_predict_global_BB_qCrI", (DL_FUNC) &C_llm_1D_fit_and_predict_global_BB_qCrI, 13, llm_1D_fit_and_predict_BB_qCrI_type},
  {"C_loo_llm_1D", (DL_FUNC) &C_loo_llm_1D, 11, loo_llm_1D_type},
  {"C_deg0_loo_llm_1D", (DL_FUNC) &C_deg0_loo_llm_1D, 7, deg0_loo_llm_1D_type},
  {"C_cv_mae_1D", (DL_FUNC) &C_cv_mae_1D, 16, cv_mae_1D_type},
  {"C_get_BB_Eyg", (DL_FUNC) &C_get_BB_Eyg, 19, get_BB_Eyg_type},
  {"C_get_Eyg_CrI", (DL_FUNC) &C_get_Eyg_CrI, 19, get_Eyg_CrI_type},
  {"C_get_Eygs", (DL_FUNC) &C_get_Eygs, 17, get_Eygs_type},
  {"C_runif_simplex", (DL_FUNC) &C_runif_simplex, 2, runif_simplex_type},
  {"C_create_ED_grid_2D", (DL_FUNC) &C_create_ED_grid_2D, 6, create_ED_grid_2D_type},
  {"C_create_ED_grid_3D", (DL_FUNC) &C_create_ED_grid_3D, 8, create_ED_grid_3D_type},
  {"C_create_ED_grid_xD", (DL_FUNC) &C_create_ED_grid_xD, 6, create_ED_grid_xD_type},
  {"C_create_ENPs_grid_2D", (DL_FUNC) &C_create_ENPs_grid_2D, 7, create_ENPs_grid_2D_type},
  {"C_create_ENPs_grid_3D", (DL_FUNC) &C_create_ENPs_grid_3D, 9, create_ENPs_grid_3D_type},
  {"C_mstree", (DL_FUNC) &C_mstree, 7, mstree_type},
  {"C_wasserstein_distance_1D", (DL_FUNC) &C_wasserstein_distance_1D, 4, wasserstein_distance_1D_type},
  {"C_get_bws",             (DL_FUNC) &C_get_bws, 6, get_bws_type},
  {"C_get_bws_with_minK_a", (DL_FUNC) &C_get_bws_with_minK_a, 6, get_bws_type},
  {"C_columnwise_weighting", (DL_FUNC) &C_columnwise_weighting, 7, columnwise_weighting_type},
  {"C_columnwise_eval", (DL_FUNC) &C_columnwise_eval, 5, columnwise_eval_type},
  {"C_columnwise_TS_norm", (DL_FUNC) &C_columnwise_TS_norm, 4, columnwise_TS_norm_type},
  {"C_matrix_wmeans", (DL_FUNC) &C_matrix_wmeans, 9, matrix_wmeans_type},
  {"C_columnwise_wmean", (DL_FUNC) &C_columnwise_wmean, 6, columnwise_wmean_type},
  {"C_columnwise_wmean_BB", (DL_FUNC) &C_columnwise_wmean_BB, 7, columnwise_wmean_BB_type},
  {"C_columnwise_wmean_BB_qCrI", (DL_FUNC) &C_columnwise_wmean_BB_qCrI, 9, columnwise_wmean_BB_qCrI_type},
  {"C_columnwise_wmean_BB_CrI_1", (DL_FUNC) &C_columnwise_wmean_BB_CrI_1, 8, columnwise_wmean_BB_CrI_type},
  {"C_columnwise_wmean_BB_CrI_2", (DL_FUNC) &C_columnwise_wmean_BB_CrI_2, 8, columnwise_wmean_BB_CrI_type},
  {"C_quantiles", (DL_FUNC) &C_quantiles, 5, quantiles_type},
  {"C_modified_columnwise_wmean_BB", (DL_FUNC) &C_modified_columnwise_wmean_BB, 7, modified_columnwise_wmean_BB_type},
  {"C_mat_columnwise_divide", (DL_FUNC) &C_mat_columnwise_divide, 4, mat_columnwise_divide_type},
  {"C_normalize_dist", (DL_FUNC) &C_normalize_dist, 7, normalize_dist_type},
  {"C_normalize_dist_with_minK_a", (DL_FUNC) &C_normalize_dist_with_minK_a, 6, normalize_dist_with_minK_a_type},
  {"C_samplewr", (DL_FUNC) &C_samplewr, 2, samplewr_type},
  {"C_vpermute", (DL_FUNC) &C_vpermute, 3, vpermute_type},
  {"C_v_get_folds", (DL_FUNC) &C_v_get_folds, 3, v_get_folds_type},
  {"C_winsorize", (DL_FUNC) &C_winsorize, 4, winsorize_type},
  {"C_pdistr", (DL_FUNC) &C_pdistr, 4, pdistr_type},
  {"C_pearson_cor", (DL_FUNC) &C_pearson_cor, 4, pearson_cor_type},
  {"C_wcov", (DL_FUNC) &C_wcov, 5, wcov_type},
  {"C_pearson_wcor", (DL_FUNC) &C_pearson_wcor, 5, pearson_wcor_type},
  {"C_pearson_wcor_BB_qCrI", (DL_FUNC) &C_pearson_wcor_BB_qCrI, 10, pearson_wcor_BB_qCrI_type},
  {"C_density_distance", (DL_FUNC) &C_density_distance, 5, density_distance_type},
  {"C_rmatrix", (DL_FUNC) &C_rmatrix, 5, rmatrix_type},
  {NULL, NULL, 0, NULL}
};

#ifdef __cplusplus
extern "C" {
#endif
SEXP _gflow_Rcpp_graph_kernel_smoother(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _gflow_rcpp_adaptive_mean_shift_gfa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _gflow_rcpp_knn_adaptive_mean_shift_gfa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
#ifdef __cplusplus
}
#endif

static const R_CallMethodDef CallMethods[] = {
  {"S_compute_hop_extremp_radii_batch", (DL_FUNC) &S_compute_hop_extremp_radii_batch, 8},
  {"S_fit_knn_riem_graph_regression", (DL_FUNC) &S_fit_knn_riem_graph_regression, 21},
  {"S_nerve_cx_spectral_filter", (DL_FUNC) &S_nerve_cx_spectral_filter, 12},
  {"S_create_nerve_complex", (DL_FUNC) &S_create_nerve_complex, 3},
  {"S_set_function_values", (DL_FUNC) &S_set_function_values, 2},
  {"S_set_weight_scheme", (DL_FUNC) &S_set_weight_scheme, 3},
  {"S_solve_full_laplacian", (DL_FUNC) &S_solve_full_laplacian, 3},
  {"S_get_simplex_counts", (DL_FUNC) &S_get_simplex_counts, 1},
  {"S_extract_skeleton_graph", (DL_FUNC) &S_extract_skeleton_graph, 1},
  {"S_construct_function_aware_graph", (DL_FUNC) &S_construct_function_aware_graph, 14},
  {"S_analyze_function_aware_weights", (DL_FUNC) &S_analyze_function_aware_weights, 12},
  {"S_apply_harmonic_extension", (DL_FUNC) &S_apply_harmonic_extension, 11},
  {"S_create_gflow_cx", (DL_FUNC) &S_create_gflow_cx, 12},
  {"S_compute_extrema_hop_nbhds", (DL_FUNC) &S_compute_extrema_hop_nbhds, 3},
  {"S_harmonic_smoother", (DL_FUNC) &S_harmonic_smoother, 9},
  {"S_perform_harmonic_smoothing", (DL_FUNC) &S_perform_harmonic_smoothing, 6},
  {"S_create_basin_cx", (DL_FUNC) &S_create_basin_cx, 3},
  {"S_find_local_extrema", (DL_FUNC) &S_find_local_extrema, 4},
  {"S_find_gflow_basins", (DL_FUNC) &S_find_gflow_basins, 6},
  {"S_amagelo", (DL_FUNC) &S_amagelo, 19},
  {"S_create_mst_completion_graph", (DL_FUNC) &S_create_mst_completion_graph, 3},
  {"S_graph_spectral_filter", (DL_FUNC) &S_graph_spectral_filter, 19},
  {"S_klaps_low_pass_smoother", (DL_FUNC) &S_klaps_low_pass_smoother, 14},
  {"S_graph_bw_adaptive_spectral_smoother", (DL_FUNC) &S_graph_bw_adaptive_spectral_smoother, 15},
  {"S_graph_kernel_smoother", (DL_FUNC) &S_graph_kernel_smoother, 17},
  {"S_graph_deg0_lowess_buffer_cv", (DL_FUNC) &S_graph_deg0_lowess_buffer_cv, 16},
  {"S_deg0_lowess_graph_smoothing", (DL_FUNC) &S_deg0_lowess_graph_smoothing, 20},
  {"S_graph_deg0_lowess_cv_mat", (DL_FUNC) &S_graph_deg0_lowess_cv_mat, 14},
  {"S_graph_deg0_lowess_cv", (DL_FUNC) &S_graph_deg0_lowess_cv, 14},
  {"S_graph_deg0_lowess", (DL_FUNC) &S_graph_deg0_lowess, 7},
  {"S_detect_local_extrema", (DL_FUNC) &S_detect_local_extrema, 6}, 
  {"S_test_monotonic_reachability_map", (DL_FUNC) &S_test_monotonic_reachability_map, 6},
  {"S_spectral_lowess_graph_smoothing", (DL_FUNC) &S_spectral_lowess_graph_smoothing, 20},
  {"S_graph_spectral_lowess_mat", (DL_FUNC) &S_graph_spectral_lowess_mat, 15},
  {"S_graph_spectral_ma_lowess", (DL_FUNC) &S_graph_spectral_ma_lowess, 13},
  {"S_graph_spectral_lowess", (DL_FUNC) &S_graph_spectral_lowess, 13},
  {"S_nada_graph_spectral_lowess", (DL_FUNC) &S_nada_graph_spectral_lowess, 13},
  {"S_compute_geodesic_stats", (DL_FUNC) &S_compute_geodesic_stats, 9},
  {"S_compute_vertex_geodesic_stats", (DL_FUNC) &S_compute_vertex_geodesic_stats, 8},
  {"S_agemalo", (DL_FUNC) &S_agemalo, 19},
  {"S_ray_agemalo", (DL_FUNC) &S_ray_agemalo, 19},
  {"S_parameterize_circular_graph", (DL_FUNC) &S_parameterize_circular_graph, 3},
  {"S_create_maximal_packing", (DL_FUNC) &S_create_maximal_packing, 5},
  {"S_find_graph_paths_within_radius", (DL_FUNC) &S_find_graph_paths_within_radius, 4},
  {"S_remove_redundant_edges", (DL_FUNC) &S_remove_redundant_edges, 2},
  {"S_compute_edge_weight_rel_deviations", (DL_FUNC) &S_compute_edge_weight_rel_deviations, 2},
  {"S_compute_edge_weight_deviations", (DL_FUNC) &S_compute_edge_weight_deviations, 2},
  {"S_adaptive_uggmalo", (DL_FUNC) &S_adaptive_uggmalo, 18},
  {"S_verify_pruning", (DL_FUNC) &S_verify_pruning, 3},
  {"S_construct_graph_gradient_flow", (DL_FUNC) &S_construct_graph_gradient_flow, 6},
  {"S_create_single_iknn_graph", (DL_FUNC) &S_create_single_iknn_graph, 5},
  {"S_create_iknn_graphs", (DL_FUNC) &S_create_iknn_graphs, 8},
  {"S_create_mknn_graph", (DL_FUNC) &S_create_mknn_graph, 2},
  {"S_create_mknn_graphs", (DL_FUNC) &S_create_mknn_graphs, 7},
  {"S_uggmalog", (DL_FUNC) &S_uggmalog, 19},
  {"S_uggmalo", (DL_FUNC) &S_uggmalo, 20},
  {"S_get_path_data", (DL_FUNC) &S_get_path_data, 10},
  {"S_ugg_get_path_data", (DL_FUNC) &S_ugg_get_path_data, 11},
  {"S_create_uniform_grid_graph", (DL_FUNC) &S_create_uniform_grid_graph, 5},
  {"S_compute_graph_analysis_sequence", (DL_FUNC) &S_compute_graph_analysis_sequence, 6},
  {"S_find_shortest_alt_path", (DL_FUNC) &S_find_shortest_alt_path, 5},
  {"S_shortest_alt_path_length", (DL_FUNC) &S_shortest_alt_path_length, 5},
  {"S_wgraph_prune_long_edges", (DL_FUNC) &S_wgraph_prune_long_edges, 5},
  {"S_graph_edit_distance", (DL_FUNC) &S_graph_edit_distance, 6},
  {"S_join_graphs", (DL_FUNC) &S_join_graphs, 4},
  {"S_convert_adjacency_to_edge_matrix", (DL_FUNC) &S_convert_adjacency_to_edge_matrix, 2},
  {"S_convert_adjacency_to_edge_matrix_set", (DL_FUNC) &S_convert_adjacency_to_edge_matrix_set, 1},
  {"S_convert_adjacency_to_edge_matrix_unordered_set", (DL_FUNC) &S_convert_adjacency_to_edge_matrix_unordered_set, 1},
  {"S_angular_wasserstein_index", (DL_FUNC) &S_angular_wasserstein_index, 3},
  {"S_compute_mstree_total_length", (DL_FUNC) &S_compute_mstree_total_length, 1},
  {"S_graph_mad", (DL_FUNC) &S_graph_mad, 2},
  {"S_graph_kmean", (DL_FUNC) &S_graph_kmean, 5},
  {"S_univariate_gkmm", (DL_FUNC) &S_univariate_gkmm, 15},
  {"S_upgmalo", (DL_FUNC) &S_upgmalo, 16},
  {"S_pgmalo", (DL_FUNC) &S_pgmalo, 17},
  {"S_upgmalog", (DL_FUNC) &S_upgmalog, 17},
  {"S_ulogit", (DL_FUNC) &S_ulogit, 8},
  {"S_eigen_ulogit", (DL_FUNC) &S_eigen_ulogit, 8},
  {"S_graph_kmean_wmad_cv", (DL_FUNC) &S_graph_kmean_wmad_cv, 9},
  {"S_graph_kmean_cv", (DL_FUNC) &S_graph_kmean_cv, 8},
  {"S_wmabilog", (DL_FUNC) &S_wmabilog, 13},
  {"S_maelog", (DL_FUNC) &S_maelog, 15},
  {"S_magelog", (DL_FUNC) &S_magelog, 15},
  {"S_adaptive_maelog", (DL_FUNC) &S_adaptive_maelog, 12},
  {"S_mabilog", (DL_FUNC) &S_mabilog, 14},
  {"S_mabilog_with_smoothed_errors", (DL_FUNC) &S_mabilog_with_smoothed_errors, 13},
  {"S_mabilo_plus", (DL_FUNC) &S_mabilo_plus, 13},
  {"S_wmabilo", (DL_FUNC) &S_wmabilo, 10},
  {"S_mabilo", (DL_FUNC) &S_mabilo, 11},
  {"S_mabilo_with_smoothed_errors", (DL_FUNC) &S_mabilo_with_smoothed_errors, 10},
  {"S_cv_imputation", (DL_FUNC) &S_cv_imputation, 12},
  {"S_prop_nbhrs_with_smaller_y", (DL_FUNC) &S_prop_nbhrs_with_smaller_y, 2},
  {"S_graph_spectrum", (DL_FUNC) &S_graph_spectrum, 2},
  {"S_graph_spectrum_plus", (DL_FUNC) &S_graph_spectrum_plus, 3},
  {"S_loc_const_vertices", (DL_FUNC) &S_loc_const_vertices, 3},
  {"S_make_response_locally_non_const", (DL_FUNC) &S_make_response_locally_non_const, 7},
  {"S_graph_diffusion_smoother", (DL_FUNC) &S_graph_diffusion_smoother, 13},
  {"S_ext_graph_diffusion_smoother", (DL_FUNC) &S_ext_graph_diffusion_smoother, 23},
  {"S_instrumented_gds", (DL_FUNC) &S_instrumented_gds, 14},
  {"S_mean_shift_data_smoother", (DL_FUNC) &S_mean_shift_data_smoother, 11},
  {"S_mean_shift_data_smoother_with_grad_field_averaging", (DL_FUNC) &S_mean_shift_data_smoother_with_grad_field_averaging, 8},
  {"S_mean_shift_data_smoother_adaptive", (DL_FUNC) &S_mean_shift_data_smoother_adaptive, 8},
  {"S_graph_spectral_smoother", (DL_FUNC) &S_graph_spectral_smoother, 17},
  {"S_graph_constrained_gradient_flow_trajectories", (DL_FUNC) &S_graph_constrained_gradient_flow_trajectories, 3},
  {"S_graph_MS_cx_with_path_search", (DL_FUNC) &S_graph_MS_cx_with_path_search, 3},
  {"S_graph_MS_cx_using_short_h_hops", (DL_FUNC) &S_graph_MS_cx_using_short_h_hops, 4},
  {"S_shortest_path", (DL_FUNC) &S_shortest_path, 3},
  {"S_mstree", (DL_FUNC) &S_mstree, 1},
  {"S_create_hHN_graph", (DL_FUNC) &S_create_hHN_graph, 3},
  {"S_create_path_graph_plus", (DL_FUNC) &S_create_path_graph_plus, 3},
  {"S_create_path_graph_plm", (DL_FUNC) &S_create_path_graph_plm, 3},
  {"S_create_path_graph_series", (DL_FUNC) &S_create_path_graph_series, 3},
  {"S_kNN", (DL_FUNC) &S_kNN, 2},
  //{"S_kNN_v2", (DL_FUNC) &S_kNN_v2, 2},
  {"S_cycle_sizes", (DL_FUNC) &S_cycle_sizes, 1},
  {"S_ecdf", (DL_FUNC) &S_ecdf, 1},
  {"S_rlaplace", (DL_FUNC) &S_rlaplace, 4},
  {"S_graph_connected_components", (DL_FUNC) &S_graph_connected_components, 1},
  {"S_estimate_local_density_over_grid", (DL_FUNC) &S_estimate_local_density_over_grid, 6},
  {"S_estimate_local_complexity", (DL_FUNC) &S_estimate_local_complexity, 5},
  {"S_estimate_binary_local_complexity", (DL_FUNC) &S_estimate_binary_local_complexity, 6},
  {"S_estimate_ma_binary_local_complexity_quadratic", (DL_FUNC) &S_estimate_ma_binary_local_complexity_quadratic, 4},
  {"S_pdistr", (DL_FUNC) &S_pdistr, 2},
  {"S_lwcor", (DL_FUNC) &S_lwcor, 3},
  {"S_lwcor_yY", (DL_FUNC) &S_lwcor_yY, 4},
  {"_gflow_Rcpp_graph_kernel_smoother", (DL_FUNC) &_gflow_Rcpp_graph_kernel_smoother, 5},
  {"_gflow_rcpp_adaptive_mean_shift_gfa", (DL_FUNC) &_gflow_rcpp_adaptive_mean_shift_gfa, 11},
  {"_gflow_rcpp_knn_adaptive_mean_shift_gfa", (DL_FUNC) &_gflow_rcpp_knn_adaptive_mean_shift_gfa, 8},
  {NULL, NULL, 0}
};

void R_init_gflow(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, CallMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);     // it will allow Rcpp’s generated .Call("...") work
}
