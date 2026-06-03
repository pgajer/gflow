# gflow Split Audit SA1: Inventory

Generated: 2026-06-03

This inventory is a first-pass file-level classification of the current `gflow`
repository into the proposed future packages:

- `gflow`: graph, geodesic, geometry, gradient-flow, and Morse-Smale core;
- `geosmooth`: conditional-expectation, smoothing, regression, trend-filtering,
  and fitted-response methods;
- `gflowx`: experimental, application-specific, retired, or GitHub-only methods;
- `shared`: metadata, generated registration, tests/helpers, third-party
  headers, or utilities whose final home requires dependency analysis;
- `undecided`: assets that need manual review in SA2/SA3.

Detailed inventory files:

- `split_audit/gflow_split_inventory.csv`
- `split_audit/gflow_split_inventory_summary.csv`
- `split_audit/scripts/build_sa1_inventory.R`

## Inventory Scope

The inventory covers:

- `DESCRIPTION` and `NAMESPACE`;
- all top-level `R/*.R` files;
- all top-level `src/*` C/C++/header/build files;
- all top-level `tests/testthat/*.R` tests and helpers;
- `vignettes/*`;
- selected package-level `inst/*` files, with vendored third-party headers
  represented as aggregate rows;
- generated `man/*.Rd` pages represented as one aggregate row.

Generated Rd pages are not classified one-by-one because they should follow the
roxygen source files after physical extraction.

## Current Package Surface

The current package has:

- 493 `export(...)` entries in `NAMESPACE`;
- 206 `S3method(...)` registrations in `NAMESPACE`;
- 168 top-level R source files;
- 186 top-level `src` files;
- 56 top-level test files/helpers.

Representative exported smoothing/regression functions include:

- `kernel.local.polynomial.cv`;
- `fit.lpl.tf`, `refit.lpl.tf`, `lpl.tf.operator`;
- `fit.slpl.tf`, `refit.slpl.tf`, `slpl.tf.operator`;
- `fit.malps`, `refit.malps`, `malps.gcv`, `malps.smoother.matrix`;
- `fit.graph.trend.filtering`, `graph.trend.filtering.operator`;
- `fit.metric.graph.lowpass`, `metric.graph.lowpass.operator`;
- `fit.ssrhe.hessian.regression`, `fit.ssrhe.hessian.regression.cv`,
  `fit.ssrhe.hessian.regression.gcv`, `fit.ssrhe.hessian.l1.regression`;
- `harmonic.smoother`, `compute.harmonic.extension`,
  `apply.harmonic.extension`;
- `data.smoother`, `meanshift.data.smoother`.

## First-Pass Counts

| Proposed package | Asset count |
|---|---:|
| `geosmooth` | 65 |
| `gflow` | 174 |
| `gflowx` | 63 |
| `shared` | 65 |
| `undecided` | 63 |

Counts are intentionally conservative.  The undecided bucket is not a failure;
it identifies the files most in need of SA2 dependency analysis.

## R Source Files

### Proposed `geosmooth`

- `R/2d_smooth_morse_smale.R`
- `R/bayes_bootstrap_rdgraph.R`
- `R/data_smoother.R`
- `R/graph_trend_filtering.R`
- `R/harmonic_extension.R`
- `R/harmonic_smoother.R`
- `R/kernel_local_polynomial_cv.R`
- `R/lpl_tf.R`
- `R/malps.R`
- `R/mean_shift_smoother.R`
- `R/metric_graph_lowpass.R`
- `R/pttf_fit.R`
- `R/pttf_geometry.R`
- `R/pttf_operator.R`
- `R/refit_rdgraph_regression.R`
- `R/riem_dcx_regression.R`
- `R/slpl_tf.R`
- `R/ssrhe_hessian_energy.R`
- `R/ulogit.R`

### Proposed `gflow`

Includes graph/geodesic/gradient-flow files such as:

- `R/compute_gfc.R`, `R/gfc_flow.R`, `R/gflow_basins.R`,
  `R/gflow_graph.R`, and related basin/cell/trajectory files;
- `R/geodesic_iknn_graphs.R`, `R/geodesic_stats.R`, `R/geodesics.R`;
- `R/graph_generators.R`, `R/graph_geodesic_distances.R`,
  `R/graph_gradient_flow.R`, `R/graph_shortest_path.R`,
  `R/graph_spectrum.R`, `R/graph_utils.R`;
- `R/iknn_graph.R`, `R/iknn_graphs.R`, `R/mknn_graphs.R`,
  `R/radius_graphs.R`, `R/sknn_graphs.R`;
- `R/isometry_deviation.R`, `R/local_geodesic_pruning.R`,
  `R/mst_completion_graphs.R`, `R/nerve_cx.R`, `R/phate_core.R`,
  `R/quadform_geodesics.R`, `R/uniform_grid_graph.R`.

### Proposed `gflowx`

Includes application or experimental layers such as:

- association/correlation testing: `R/fassoc_test.R`, `R/fassoc0_test_paired.R`,
  `R/fassoc1_test_paired.R`, `R/gfassoc_utils.R`, `R/lcor.R`,
  `R/lslope.R`;
- clustering and microbiome/application workflow files:
  `R/clustering.R`, `R/consensus_clustering.R`, `R/cst_colors.R`,
  `R/cst_graph_mixing_stats.R`, `R/select_phylotypes_for_assoc.R`;
- trajectory and visualization workflow files:
  `R/trajectory_clustering.R`, `R/trajectory_graph_clustering.R`,
  `R/plot_trajectory_hurdle_association.R`,
  `R/visualize_mgcp_cluster_single.R`, `R/visualize_mgcp_clusters.R`;
- older/exploratory runtime helpers: `R/malo_runtime_helpers.R`.

### Shared R Infrastructure

- `R/RcppExports.R`
- `R/compat_helpers.R`
- `R/globals.R`
- `R/kernels.R`
- `R/knn_cache_helpers.R`
- `R/knn_metric_helpers.R`
- `R/local_pca_chart_dim.R`
- `R/plot_utils.R`
- `R/random_sampling.R`
- `R/row_wise_operations.R`
- `R/synthetic_data_utils.R`
- `R/utils.R`

### Undecided R Files

These files need SA2 dependency analysis and/or an explicit product decision:

- `R/annotation.R`
- `R/centered_paths.R`
- `R/density.R`
- `R/detect_binary_partition.R`
- `R/distance_quantile_bin_analysis.R`
- `R/divergences.R`
- `R/find_discrepant_pairs.R`
- `R/fit_rho_randomwalk.R`
- `R/generate_knn_example_data.R`
- `R/gflow-package.R`
- `R/grids.R`
- `R/madag.R`
- `R/major_arms.R`
- `R/partition_heatmap.R`
- `R/plot_feature_distance_diagnostics.R`
- `R/preprocess_matrix.R`
- `R/scan_resolution.R`
- `R/scan_stability.R`
- `R/se_tree.R`
- `R/select_pts.R`
- `R/subject_neighborhood_stats.R`
- `R/tube_lens_corridor.R`
- `R/vertex_density.R`
- `R/zzz.R`

## Compiled-Code Observations

The `src` inventory already shows why SA3 is necessary.  There are clear
compiled-code clusters:

- graph/geodesic infrastructure likely to remain in `gflow`;
- smoother/regression kernels likely to move with `geosmooth`;
- association/experimental/application kernels likely to move to `gflowx`;
- shared numerical helpers such as kNN, local PCA charts, kernels, weighted
  linear algebra, Rcpp registration, and OpenMP diagnostics.

No compiled files should be moved before SA2/SA3 map which R functions call
which native symbols.

## Immediate SA1 Conclusions

1. `geosmooth` is a coherent package candidate.  It has a recognizable method
   family centered on geometric smoothers and conditional expectation
   estimators.
2. `gflow` remains coherent if it keeps graph construction, graph geodesics,
   gradient-flow cells, basins, trajectories, Morse-Smale objects, and geometry
   diagnostics.
3. `gflowx` is useful.  There are enough experimental/application-specific
   assets to justify moving them out of the CRAN-facing core instead of deleting
   them.
4. The first extraction should not be LPL-TF or SLPLiFT.  They are too entangled
   with recently changing local chart, support, selector, and compiled
   infrastructure.
5. The recommended pilot remains `kernel.local.polynomial.cv`, after SA2/SA3
   confirm its dependency closure.

## SA1 Exit Criteria

SA1 is complete for audit purposes when:

- the inventory CSV exists and can be regenerated;
- the top-level package assignment counts are recorded;
- the main R source families are classified;
- undecided files are explicitly listed;
- no physical source movement has occurred.
