#ifndef MSC2_H_
#define MSC2_H_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

/* Functions from stats_utils.c */
void C_permute(int *x, int n);

/*!
    \fn CHECK_PTR(p)

    \brief A macro testing a pointer for NULL value

    Exits if the pointer is NULL. Replacing
    fprintf(stderr,"ERROR: Out of memory in file %s at line %d\n",__FILE__,__LINE__);
    exit (EXIT_FAILURE);

    \param p   - a pointer
*/
#define CHECK_PTR(p) \
    do { \
        if ((p) == NULL) { \
            Rf_error("Memory allocation failed in file %s at line %d.\n", __FILE__, __LINE__); \
        } \
    } while(0)

#ifdef __cplusplus
extern "C" {
#endif

    // Beta functions
    void C_llm_1D_beta(const double *nn_x,
                       const double *nn_y,
                       double *nn_w,
                       const int    *maxK,
                       const int    *rK,
                       const int    *rng,
                       const int    *rdeg,
                       double *beta);

    void C_llm_1D_beta_perms(const double *Tnn_x,
                             const int    *Tnn_i,
                             const double *y,
                             const int    *ny,
                             double       *Tnn_w,
                             const int    *maxK,
                             const int    *rnrTnn,
                             const int    *rncTnn,
                             const int    *rdeg,
                             const int    *rn_perms,
                             double *beta_perms);

// Prediction functions
    void C_predict_1D(const  double *beta,
                      const int    *nn_i,
                      const double *nn_d,
                      const double *nn_x,
                      const int    *maxK,
                      const int    *rK,
                      const int    *rng,
                      const int    *rdeg,
                      const int    *rikernel,
                      const int    *rnx,
                      double *Ey);

    void C_wpredict_1D(const double *beta,
                       const int    *nn_i,
                       const double *nn_w,
                       const double *nn_x,
                       const int    *maxK,
                       const int    *rK,
                       const int    *rng,
                       const int    *rdeg,
                       const int    *rnx,
                       const int    *rybinary,
                       double *Ey);

// Fit and predict functions
    void C_llm_1D_fit_and_predict(const int    *nn_i,
                                  double *nn_w,
                                  const double *nn_x,
                                  const double *nn_y,
                                  const int    *rybinary,
                                  const int    *maxK,
                                  const int    *rK,
                                  const int    *rng,
                                  const int    *rnx,
                                  const int    *rdeg,
                                  double *Ey,
                                  double *beta);

    void C_mllm_1D_fit_and_predict(const double *TY,
                                   const int    *nrTY,
                                   const int    *ncTY,
                                   const int    *rybinary,
                                   const int    *Tnn_i,
                                   double *Tnn_w,
                                   const double *Tnn_x,
                                   const int    *nrTnn,
                                   const int    *ncTnn,
                                   const int    *maxK,
                                   const int    *deg,
                                   double *ETY);

// LOO functions
    void C_loo_llm_1D(const double *x,
                      const int    *rnx,
                      const int    *nn_i,
                      double *nn_w,
                      const double *nn_x,
                      const double *nn_y,
                      const int    *rng,
                      const int    *rK,
                      const int    *rdeg,
                      const int    *rmax_nmodels,
                      double *Ey);

    void C_deg0_loo_llm_1D(const int    *rnx,
                           const int    *nn_i,
                           double *nn_w,
                           const double *nn_y,
                           const int    *rK,
                           const int    *rng,
                           double *Ey);

// CV functions
    void C_cv_mae_1D(const int    *rnfolds,
                     const int    *rnreps,
                     const double *y,
                     const int    *nn_i,
                     const double *nn_d,
                     const double *nn_x,
                     const double *nn_y,
                     const int    *rybinary,
                     const int    *rK,
                     const int    *rng,
                     const int    *rnx,
                     const int    *rdeg,
                     const double *rbw,
                     const int    *rminK,
                     const int    *rikernel,
                     double *rMAE);

// Bootstrap confidence interval functions
    void C_llm_1D_fit_and_predict_BB_CrI(const int    *nn_i,
                                         const double *nn_w,
                                         const double *nn_x,
                                         const double *nn_y,
                                         const int    *rybinary,
                                         const int    *maxK,
                                         const int    *rK,
                                         const int    *rng,
                                         const int    *rnx,
                                         const int    *rdeg,
                                         const int    *rnBB,
                                         const double *Ey,
                                         double *Ey_CI);

    void C_llm_1D_fit_and_predict_global_BB_CrI(const int    *nn_i,
                                                const double *nn_w,
                                                const double *nn_x,
                                                const double *nn_y,
                                                const int    *rybinary,
                                                const int    *maxK,
                                                const int    *rK,
                                                const int    *rng,
                                                const int    *rnx,
                                                const int    *rdeg,
                                                const int    *rnBB,
                                                const double *Ey,
                                                double *Ey_CI);

    void C_llm_1D_fit_and_predict_global_BB(const int    *nn_i,
                                            const double *nn_w,
                                            const double *nn_x,
                                            const double *nn_y,
                                            const int    *rybinary,
                                            const int    *maxK,
                                            const int    *rK,
                                            const int    *rng,
                                            const int    *rnx,
                                            const int    *rdeg,
                                            const int    *rnBB,
                                            double *gbbEy);

    void C_llm_1D_fit_and_predict_global_BB_qCrI(const int    *rybinary,
                                                 const int    *nn_i,
                                                 const double *nn_w,
                                                 const double *nn_x,
                                                 const double *nn_y,
                                                 const int    *maxK,
                                                 const int    *rK,
                                                 const int    *rng,
                                                 const int    *rnx,
                                                 const int    *rdeg,
                                                 const int    *rnBB,
                                                 const double *ralpha,
                                                 double *Ey_CI);

// Get expected values functions
    void C_get_BB_Eyg(const int    *rn_BB,
                      const int    *Tnn_i,
                      const double *Tnn_x,
                      const double *Tnn_y,
                      const int    *rybinary,
                      const double *Tnn_w,
                      const int    *rnx,
                      const int    *rnrTnn,
                      const int    *rncTnn,
                      const int    *max_K,
                      const int    *rdegree,
                      const int    *Tgrid_nn_i,
                      const double *Tgrid_nn_x,
                      const double *Tgrid_nn_w,
                      const int    *rnrTgrid_nn,
                      const int    *rncTgrid_nn,
                      const int    *grid_max_K,
                      double *bb_beta,
                      double *bb_Eyg);

    void C_get_Eyg_CrI(const int    *rybinary,
                       const int    *rn_BB,
                       const int    *Tnn_i,
                       const double *Tnn_x,
                       const double *Tnn_y,
                       const double *Tnn_w,
                       const int    *rnx,
                       const int    *rnrTnn,
                       const int    *rncTnn,
                       const int    *max_K,
                       const int    *rdegree,
                       const int    *Tgrid_nn_i,
                       const double *Tgrid_nn_x,
                       const double *Tgrid_nn_w,
                       const int    *rnrTgrid_nn,
                       const int    *rncTgrid_nn,
                       const int    *grid_max_K,
                       const double *ralpha,
                       double *Eyg_CrI);

    void C_get_Eygs(const double *bws,
                    const int    *rn_bws,
                    const int    *Tnn_i,
                    const double *Tnn_d,
                    const double *Tnn_x,
                    const double *Tnn_y,
                    const int    *rybinary,
                    const int    *rnrTnn,
                    const int    *rncTnn,
                    const int    *rdegree,
                    const int    *rminK,
                    const int    *Tgrid_nn_i,
                    const double *Tgrid_nn_d,
                    const double *Tgrid_nn_x,
                    const int    *rnrTgrid_nn,
                    const int    *rncTgrid_nn,
                    double *Eygs);

    void C_mstree(const int *riinit,
                  const int *nn_i,
                  const double *nn_d,
                  const double *rldist,
                  const int *rn,
                  int *edges,
                  double *edge_lens );

    // ------------------------------------------------------------
    // C_create_ED_grid_*D
    // ------------------------------------------------------------

    void C_create_ED_grid_2D(const double *rdx,
                             const double *rx1L,
                             const int    *rn1,
                             const double *rx2L,
                             const int    *rn2,
                             double *grid);

    void C_create_ED_grid_3D(const double *rdx,
                             const double *rx1L,
                             const int    *rn1,
                             const double *rx2L,
                             const int    *rn2,
                             const double *rx3L,
                             const int    *rn3,
                             double *grid);

    void C_create_ED_grid_xD(const double *rw,
                             const double *L,
                             const int    *rdim,
                             const int    *size_d,
                             const int    *rTotElts,
                             double *grid);

    // ------------------------------------------------------------
    // C_create_ENPs_grid_*D
    // ------------------------------------------------------------

    void C_create_ENPs_grid_2D(const int *rn,
                               const double *rx1L,
                               const double *rx1R,
                               const double *rx2L,
                               const double *rx2R,
                               const double *rf,
                               double *grid);


    void C_create_ENPs_grid_3D(const int *rn,
                               const double *rx1L,
                               const double *rx1R,
                               const double *rx2L,
                               const double *rx2R,
                               const double *rx3L,
                               const double *rx3R,
                               const double *rf,
                               double *grid);

    void print_2d_double_array_first_n(const double *x, int nr, int nc, int N, int precision);
    void print_int_array(const int *x, int n);
    void print_double_array(const double *x, int n);
    void print_double_array_with_precision(const double *x, int n, int precision);

    SEXP S_angular_wasserstein_index(SEXP s_X, SEXP s_Y, SEXP s_k);

    SEXP S_compute_mstree_total_length(SEXP s_X);

    SEXP S_mstree(SEXP X);

    SEXP S_mst_kNN_graph(SEXP RX, SEXP Rk);

    SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h);

    SEXP S_create_path_graph_series(SEXP s_adj_list,
                                    SEXP s_weight_list,
                                    SEXP s_h_values);

    SEXP S_create_path_graph_plus(SEXP s_adj_list,
                                  SEXP s_edge_length_list,
                                  SEXP s_h);

    SEXP S_create_path_graph_plm(SEXP s_adj_list,
                                 SEXP s_edge_length_list,
                                 SEXP s_h);

    SEXP S_shortest_path(SEXP s_graph, SEXP s_edge_lengths, SEXP s_vertices);
    SEXP S_join_graphs(SEXP Rgraph1, SEXP Rgraph2, SEXP Ri1, SEXP Ri2);
    SEXP S_convert_adjacency_to_edge_matrix(SEXP s_graph, SEXP s_weights);
    SEXP S_convert_adjacency_to_edge_matrix_set(SEXP s_graph);
    SEXP S_convert_adjacency_to_edge_matrix_unordered_set(SEXP s_graph);

    SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    SEXP S_graph_MS_cx_with_path_search(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    SEXP S_graph_MS_cx_using_short_h_hops(SEXP s_graph,
                                          SEXP s_hop_list,
                                          SEXP s_core_graph,
                                          SEXP s_Ey);
    SEXP S_make_response_locally_non_const(SEXP Rgraph,
                                           SEXP Ry,
                                           SEXP Rweights,
                                           SEXP Rstep_factor,
                                           SEXP Rprec,
                                           SEXP Rn_itrs,
                                           SEXP Rmean_adjust);

    SEXP S_graph_edit_distance(SEXP s_graph1_adj_list,
                               SEXP s_graph1_weights_list,
                               SEXP s_graph2_adj_list,
                               SEXP s_graph2_weights_list,
                               SEXP s_edge_cost,
                               SEXP s_weight_cost_factor);

    SEXP S_graph_spectral_smoother(SEXP Rgraph,
                                   SEXP Rd,
                                   SEXP Rweights,
                                   SEXP Ry,
                                   SEXP Rimputation_method,
                                   SEXP Rmax_iterations,
                                   SEXP Rconvergence_threshold,
                                   SEXP Rapply_binary_threshold,
                                   SEXP Rbinary_threshold,
                                   SEXP Rikernel,
                                   SEXP Rdist_normalization_factor,
                                   SEXP Rn_CVs,
                                   SEXP Rn_CV_folds,
                                   SEXP Repsilon,
                                   SEXP Rmin_plambda,
                                   SEXP Rmax_plambda,
                                   SEXP Rseed);

    SEXP S_graph_spectrum(SEXP Rgraph, SEXP Rnev);
    SEXP S_graph_spectrum_plus(SEXP Rgraph, SEXP Rnev, SEXP Rreturn_dense);

    SEXP S_mabilog(SEXP s_x,
                   SEXP s_y,
                   SEXP s_y_true,
                   SEXP s_max_iterations,
                   SEXP s_ridge_lambda,
                   SEXP s_max_beta,
                   SEXP s_tolerance,
                   SEXP s_k_min,
                   SEXP s_k_max,
                   SEXP s_n_bb,
                   SEXP s_p,
                   SEXP s_distance_kernel,
                   SEXP s_dist_normalization_factor,
                   SEXP s_verbose);

    SEXP S_mabilo(SEXP s_x,
                  SEXP s_y,
                  SEXP s_y_true,
                  SEXP s_k_min,
                  SEXP s_k_max,
                  SEXP s_n_bb,
                  SEXP s_p,
                  SEXP s_distance_kernel,
                  SEXP s_dist_normalization_factor,
                  SEXP s_epsilon,
                  SEXP s_verbose);

    SEXP S_mabilo_plus(SEXP s_x,
                       SEXP s_y,
                       SEXP s_y_true,
                       SEXP s_w,
                       SEXP s_k_min,
                       SEXP s_k_max,
                       SEXP s_model_averaging_strategy,
                       SEXP s_error_filtering_strategy,
                       SEXP s_distance_kernel,
                       SEXP s_model_kernel,
                       SEXP s_dist_normalization_factor,
                       SEXP s_epsilon,
                       SEXP s_verbose);

    SEXP S_prop_nbhrs_with_smaller_y(SEXP Rgraph, SEXP Ry);

    SEXP S_find_shortest_alt_path(SEXP s_adj_list,
                                  SEXP s_isize_list,
                                  SEXP s_source,
                                  SEXP s_target,
                                  SEXP s_edge_isize);

    SEXP S_shortest_alt_path_length(SEXP s_adj_list,
                                    SEXP s_isize_list,
                                    SEXP s_source,
                                    SEXP s_target,
                                    SEXP s_edge_isize);

    SEXP S_wgraph_prune_long_edges(SEXP s_adj_list,
                                   SEXP s_edge_length_list,
                                   SEXP s_alt_path_len_ratio_thld,
                                   SEXP s_use_total_length_constraint,
                                   SEXP s_verbose);

    SEXP S_cycle_sizes(SEXP RA);

    void C_flm(const double *x,
               const double *y,
               const int    *rnr,
               const int    *rnc,
               double *beta);

    // Stats utility functions (in stats_utils.c)
    double mean(const double *x, int n );
    double wmean(const double *x, const double *w, int n);
    double median(double *a, int n );

    void C_quantiles(const double *x, const int *rn, const double *probs, const int *rnprobs, double *quants);
    void C_winsorize(const double *y, const int *rn, const double *rp, double *wy);

    SEXP S_ecdf(SEXP x);

    SEXP S_rlaplace(SEXP R_n, SEXP R_location, SEXP R_scale, SEXP R_seed);
    SEXP S_graph_connected_components(SEXP R_graph);

    void C_get_bws(const double *d,
                   const int *rnr,
                   const int *rnc,
                   const int *rminK,
                   const double *rbw,
                   double *bws);

    void C_get_bws_with_minK_a(const double *d,
                               const int *rnr,
                               const int *rnc,
                               const int *minK,
                               const double *rbw,
                               double *bws);

    void C_columnwise_weighting(const double *x,
                            const int    *rnr,
                            const int    *rnc,
                            const int    *rikernel,
                            const double *bws,
                                  int    *maxK,
                                  double *w);

    void C_columnwise_TS_norm(const double *x,
                              const int *rnr,
                              const int *rnc,
                              double *nx);

#ifdef __cplusplus
}
#endif


#endif // MSC2_H_
