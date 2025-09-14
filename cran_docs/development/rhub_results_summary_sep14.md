# Summary of rhub rchk results

Here is a list of gflow package functions that rchk flaged, grouped by file and annotated with the key issue type(s) rchk reported (using your file\:line hints where present).

# Functions to fix (from rchk output)

## agemalo / ray_agemalo / amagelo family

 S_agemalo (`src/agemalo_r.cpp:384, 386`) — negative depth; UNPROTECT(<counter>) over-unprotects; stack imbalance.
  Lambdas: `:285`, `:303`, `:312` — possible stack imbalance.
 S_ray_agemalo (`src/ray_agemalo_r.cpp:384, 386`) — same pattern as above (multiple over-unprotects; negative depth).
  Lambdas: `:285`, `:303`, `:312` — possible stack imbalance.
 S_amagelo (`src/amagelo_r.cpp:303`) — unsupported UNPROTECT(variable).
  Lambdas: `:169`, `:175`, `:187` — possible stack imbalance.

## Graph smoothing / spectral / diffusion

 S_graph_kernel_smoother (`src/graph_kernel_smoother_r.cpp:285`) — multiple pointer protection counters; unsupported UNPROTECT(variable).
 S_graph_bw_adaptive_spectral_smoother (`src/graph_bw_adaptive_spectral_smoother_r.cpp:162`) — unsupported UNPROTECT(variable).
  Lambdas: `:135`, `:148`, `:153` — possible stack imbalance.
 S_graph_spectral_lowess (`src/graph_spectral_lowess_r.cpp:187`) — unsupported UNPROTECT(variable).
  Lambda: `:174` — possible stack imbalance.
 S_graph_spectral_lowess_mat (`src/graph_spectral_lowess_mat_r.cpp:226, 228`) — negative depth; over-unprotect; stack imbalance.
 S_graph_spectral_ma_lowess (`src/graph_spectral_ma_lowess_r.cpp:158, 160`) — negative depth; stack imbalance.
  Lambda: `:147` — possible stack imbalance.
 S_graph_diffusion_smoother (`src/graph_diffusion_smoother.cpp:3291`) — unsupported UNPROTECT(variable).
  Lambdas: `:3235`, `:3251` — possible stack imbalance.
 S_graph_deg0_lowess_cv (`src/graph_deg0_lowess_cv.cpp:559`) — unsupported UNPROTECT(variable).
  Lambda: `:521` — possible stack imbalance.
 S_graph_deg0_lowess_buffer_cv (`src/graph_deg0_lowess_buffer_cv_r.cpp:209`) — unsupported UNPROTECT(variable).
  Lambda: `:165` — possible stack imbalance.
 S_graph_deg0_lowess_cv_mat (`src/graph_deg0_lowess_cv_mat_r.cpp:198`) — protect stack too deep; unsupported UNPROTECT(variable).
  Lambda: `:142` — possible stack imbalance.

## Graph construction / utilities / paths

 S_construct_graph_gradient_flow (`src/graph_gradient_flow.cpp:1518–1534, 1683`) — unprotected `r_traj` during allocations; unsupported UNPROTECT(not const/variable).
 S_join_graphs (`src/graph_utils.cpp:275, 278`) — negative depth; over-unprotect; stack imbalance.
 S_find_graph_paths_within_radius (`src/set_wgraph.cpp:490, 492, 495, 510, 512, 552, 554, 558`) — unprotected intermediates (`r_paths_list`, `result`); negative depth; over-unprotect; stack imbalance.
  Lambdas: `:484`, `:507`, `:550` — possible stack imbalance.
 S_find_shortest_alt_path (`src/pruning_long_edges.cpp:794, 797`) — negative depth; over-unprotect; stack imbalance.
 S_remove_redundant_edges (`src/set_wgraph.cpp:1248`) — unsupported UNPROTECT(not const/variable).
 S_compute_edge_weight_rel_deviations (`src/set_wgraph.cpp:1031`) — unsupported UNPROTECT(not const/variable).

## Weighted graph pruning / H/HN graphs / MST / nerve

 S_wgraph_prune_long_edges (`src/pruning_long_edges.cpp:1307, 1309, 1320, 1330, 1337`) — unprotected `pruned_adj_list`, `pruned_edge_lengths_list`, `R_path_lengths`, `R_edge_lengths` during allocations.
 S_create_hHN_graph (`src/hHN_graphs.cpp:387–389`) — unprotected `res` during allocations (`Rf_allocVector`, `Rf_mkChar`).
 S_create_mst_completion_graph (`src/mst_completion_graphs_r.cpp:216`) — unsupported UNPROTECT(not const/variable).
 S_create_nerve_complex (`src/nerve_cx_r.cpp:58, 65`) — unprotected `complex_ptr` during allocations.

## Density / stats

 S_estimate_local_density_over_grid (`src/density.cpp:463, 528`) — negative depth; over-/under-protect; stack imbalance.
 S_compute_geodesic_stats (`src/geodesic_stats_r.cpp:268`) — unsupported UNPROTECT(not const/variable).

## Basins / centered paths

 S_create_basin_cx (`src/gflow_basins_r.cpp:692`) — unsupported UNPROTECT(not const/variable).
  Lambda: `:384` — possible stack imbalance.
 S_ugg_get_path_data (`src/centered_paths.cpp:2700`) — unsupported UNPROTECT(not const/variable).

## Parameterization

 S_parameterize_circular_graph (`src/parameterize_circular_graph_r.cpp:89`) — unsupported UNPROTECT(variable).
  Lambda: `:73` — possible stack imbalance.


## PGMALO / UPGMALO / UPGMALOG

 S_pgmalo (`src/pgmalo.cpp:1406`) — multiple pointer protection counters; unsupported UNPROTECT(variable).
 S_upgmalo (`src/pgmalo.cpp:791`) — multiple pointer protection counters; unsupported UNPROTECT(variable).
 S_upgmalog (`src/pgmalog.cpp:1269, 1272`) — negative depth; over-unprotect; stack imbalance.

## “nada” spectral lowess

 S_nada_graph_spectral_lowess (`src/nada_graph_spectral_lowess_r.cpp:214`) — unsupported UNPROTECT(variable).
  Lambda: `:177` — possible stack imbalance.

## kNN graphs / Rcpp wrappers (C API)

 create_iknn_graph(SEXPREC, SEXPREC)\\ (`src/iknn_graphs.cpp:652, 661`) — `RX` unprotected before `Rf_coerceVector`; calling `S_kNN` with fresh (unprotected) pointers.
 _gflow_Rcpp_graph_kernel_smoother, _gflow_rcpp_adaptive_mean_shift_gfa, _gflow_rcpp_knn_adaptive_mean_shift_gfa — rchk “address taken → ignored” (likely benign in Rcpp glue, but worth a quick glance).
 create_r_graph_from_set_wgraph(set_wgraph_t const&) — “address taken → ignored” (likely benign).

---

## Notes about non-package noise (safe to ignore for fixes)

 `objdump: Warning: Unrecognized form: 0x23` and the DIE/abbrev warnings: toolchain/DWARF noise.
 “too many states (abstraction error?)” in `strptime_internal`, `bcEval_loop`, `RunGenCollect`: these are in R internals; ignore.
 Rcpp `Armor`, `Shield`, `Vector::push_back_name__impl`: these are header-level generic patterns; treat as noise unless you’re directly using `Rcpp::Armor`/`Shield` incorrectly.

---

## Triage (where to start)

1. Hard errors (negative depth / over-UNPROTECT / unsupported UNPROTECT(variable))

    `S_agemalo`, `S_ray_agemalo`, `S_amagelo`
    `S_graph_` functions listed above (esp. `_lowess_mat`, `_diffusion_`)
    `S_construct_graph_gradient_flow`, `S_find_graph_paths_within_radius`, `S_wgraph_prune_long_edges`
    `S_join_graphs`, `S_estimate_local_density_over_grid`
2. Unprotected variables during allocation

    `S_construct_graph_gradient_flow`, `S_find_graph_paths_within_radius`, `S_wgraph_prune_long_edges`, `S_create_hHN_graph`, `S_create_nerve_complex`, `create_iknn_graph`
3. “multiple pointer protection counters” (refactor to single `nprot` counter; container-first protection)

    `S_graph_kernel_smoother`, `S_mabilo`, `S_mabilog`, `S_pgmalo`, `S_upgmalo`

If you want, I can generate patch skeletons for a couple of the worst offenders (e.g., `S_agemalo`, `S_graph_kernel_smoother`, `S_construct_graph_gradient_flow`) showing container-first PROTECT patterns and constant-count `UNPROTECT(nprot)` cleanups.


