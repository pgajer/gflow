# Quadform Utilities And Graph Constructors

This note is the shared reference for quadform-derived graph construction and
validation work in `gflow`, `grip`, `geodesicMDS`, and the SIMODS manuscript
support projects. The canonical copy lives in:

```text
/Users/pgajer/current_projects/gflow/dev/quadforms/quadform_gflow_graph_constructor_notes.md
```

The `grip` repository should point to this file by symlink:

```text
/Users/pgajer/current_projects/grip/dev/quadforms/quadform_gflow_graph_constructor_notes.md
```

The guiding rule for grip validation is simple: do not reimplement quadform
geometry, exact edge lengths, reference geodesics, graph constructors, or MST
repair in grip. Use the exported `gflow` API whenever possible.

## Current Dataset API

`quadform.sample.dataset()` now handles both 2D and 3D sample/reference
datasets through its `dim` parameter.

Source:

```text
/Users/pgajer/current_projects/gflow/R/quadform_geodesics.R
```

Signature:

```r
quadform.sample.dataset(
  n,
  index.k,
  coefficients = NULL,
  domain.radius = 1,
  sample.method = NULL,
  grid.size = 51,
  sample.connection.k = 8,
  seed = NULL,
  dim = 2,
  domain.shape = NULL,
  n.ref = 5000,
  candidate.multiplier = 6,
  boundary.fraction = 0.2,
  epsilon = NULL,
  edge.length.factor = 4,
  delaunay.backend = c("cpp", "geometry"),
  qhull.options = "Qt Qbb Qc"
)
```

For `dim = 2`, the function preserves the legacy grid-reference workflow:

```r
ds2 <- gflow::quadform.sample.dataset(
  n = 120,
  dim = 2,
  index.k = 1,
  coefficients = c(1, 2),
  domain.radius = 1,
  domain.shape = "disk",
  sample.method = "uniform.parameter.disk",
  grid.size = 101,
  sample.connection.k = 8,
  seed = 1
)
```

For `dim = 3`, the function samples parameter points in a ball or cube and
uses `quadform.delaunay.geodesic.distances()` for the reference:

```r
ds3 <- gflow::quadform.sample.dataset(
  n = 120,
  dim = 3,
  index.k = 1,
  coefficients = c(1, 2, 4),
  domain.radius = 1,
  domain.shape = "ball",
  sample.method = "uniform.parameter.ball",
  n.ref = 5000,
  candidate.multiplier = 6,
  boundary.fraction = 0.25,
  edge.length.factor = 4,
  delaunay.backend = "cpp",
  seed = 1
)
```

Return structure is intentionally similar for both dimensions:

- `X_param`: parameter coordinates.
- `X_embed`: embedded quadratic graph coordinates.
- `q`: quadratic-form values.
- `D_geodesic` and `distances`: sample-by-sample reference distances.
- `reference`: reference graph/grid payload.
- `metadata`: dataset parameters, including `dim`, `index_k`,
  `coefficients`, `domain_shape`, `sample_method`, and reference settings.

## Public Quadform Functions

All public quadform functions are defined in:

```text
/Users/pgajer/current_projects/gflow/R/quadform_geodesics.R
```

Core geometry:

- `quadform.embed(X, index.k, coefficients = NULL)`
- `quadform.gradient(X, index.k, coefficients = NULL)`
- `quadform.metric(X, index.k, coefficients = NULL)`

Exact segment lengths:

- `quadform.edge.length(u, v, index.k, coefficients = NULL, tol = ...)`
- `quadform.edge.lengths(U, V, index.k, coefficients = NULL)`

2D reference geodesics:

- `quadform.reference.geodesics()`
- `quadform.grid.geodesic.distances()`
- `quadform.grid.geodesic.calibration()`

3D reference geodesics:

- `quadform.delaunay.geodesic.distances()`

Dataset facade:

- `quadform.sample.dataset(..., dim = 2)`
- `quadform.sample.dataset(..., dim = 3)`

## Internal Quadform Helpers

The exported dataset and reference functions are the preferred API. Internal
helpers are listed here for audit and provenance, not as grip dependencies.

Validation and shared helpers:

- `.validate.quadform.index()`
- `.validate.quadform.data.matrix()`
- `.quadform.signs()`
- `.validate.quadform.coefficients()`
- `.validate.quadform.domain.shape()`
- `.quadform.default.domain.radius()`
- `.validate.quadform.reference.domain.shape()`
- `.quadform.points.inside.domain()`
- `.quadform.points.inside.domain.nd()`
- `.quadform.domain.label()`
- `.quadform.value()`
- `.with.quadform.seed()`

2D sampling and grid-reference helpers:

- `.quadform.reference.vertices.2d()`
- `.quadform.reference.grid.edges.2d()`
- `.quadform.reference.sample.edges()`
- `.quadform.sample.uniform.parameter.disk()`
- `.quadform.sample.radial.parameter.disk()`
- `.quadform.sample.uniform.parameter.square()`
- `.quadform.sample.radial.parameter.square()`
- `.quadform.sample.parameter()`
- `.quadform.default.oracle.tube.radius()`
- `.quadform.reference.grid.2d()`
- `.sample.quadform.reference.pairs()`
- `.summarize.quadform.grid.error()`

3D sampling, epsilon-net, and Delaunay-reference helpers:

- `.quadform.sample.parameter.3d()`
- `.quadform.sample.boundary.3d()`
- `.quadform.sample.parameter.3d.method()`
- `.quadform.domain.volume.3d()`
- `.quadform.cell.key()`
- `.quadform.greedy.epsilon.net()`
- `.quadform.epsilon.net.3d()`
- `.quadform.delaunay.edges.3d.geometry()`
- `.quadform.delaunay.edges.3d.cpp()`
- `.quadform.delaunay.edges.3d()`
- `.quadform.edge.list.to.adj()`
- `.quadform.adj.components()`
- `.quadform.filter.delaunay.edges()`

## C++ And Rcpp Assets

Rcpp wrappers:

```text
/Users/pgajer/current_projects/gflow/R/RcppExports.R
```

- `rcpp_quadform_delaunay_edges_3d()`
- `rcpp_quadform_edge_lengths()`
- `rcpp_quadform_grid_pair_distances()`
- `rcpp_quadform_grid_geodesic_distances()`

C++ sources:

```text
/Users/pgajer/current_projects/gflow/src/quadform_grid_geodesics.cpp
/Users/pgajer/current_projects/gflow/src/quadform_delaunay_edges.cpp
```

These provide the exact vectorized quadratic segment-length kernel, 2D
grid-reference distance computation, and 3D Delaunay edge extraction.

## Tests

Primary package tests:

```text
/Users/pgajer/current_projects/gflow/tests/testthat/test-quadform-geodesics.R
```

Coverage includes:

- embedding, gradient, and metric identities;
- scalar and vectorized exact edge lengths in 2D and 3D;
- 2D grid-reference geodesics;
- 3D Delaunay-reference geodesics;
- C++ versus `geometry` Delaunay parity;
- edge-filtering diagnostics;
- `quadform.grid.geodesic.calibration()`;
- `quadform.sample.dataset()` for 2D and 3D;
- reproducible local seed behavior;
- input validation.

## Graph Constructors For Quadform Samples

Use exported graph constructors from:

```text
/Users/pgajer/current_projects/gflow/R/radius_graphs.R
```

Preferred constructors for the first grip MISF edge-KK validation pass:

- `create.adaptive.radius.graph()`
- `create.cknn.graph()`

Both support component repair:

```r
connect.components = TRUE
connect.method = "component.mst"
```

Alternative repair methods:

```r
connect.method = "component.mst.ann"
connect.method = "global.mst"
```

The default exact validation option should be `"component.mst"` unless graph
size forces an ANN bridge mode.

### Continuous kNN

```r
g <- gflow::create.cknn.graph(
  X,
  k.scale = 8L,
  delta = 1,
  prune.method = "none",
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

`create.cknn.graph()` is a wrapper around `create.adaptive.radius.graph()` with
`radius.rule = "geomean"` and `radius.factor = delta`.

### Adaptive Radius

```r
g <- gflow::create.adaptive.radius.graph(
  X,
  k.scale = 8L,
  radius.factor = 1,
  radius.rule = "max",
  prune.method = "none",
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

The `radius.rule` setting can be swept later (`"max"`, `"min"`, `"geomean"`).
For first-pass grip validation, use `"max"` for adaptive-radius and
`create.cknn.graph()` for continuous kNN.

### Output Fields To Capture

Final graph:

```r
g$edge_matrix
g$edge_weight
g$adj_list
g$weight_list
```

MST repair diagnostics:

```r
g$n_components_before
g$n_components_after
g$n_mst_edges_added
g$mst_edge_matrix
g$mst_edge_weight
g$connect_components
g$connect_method
g$bridge_method
g$bridge_exact_fallback_used
```

For validation reports, record:

- components before repair;
- components after repair;
- number of MST edges added;
- median and maximum MST repair edge weight;
- repair-edge weight relative to median original edge weight.

## Edge Weights For grip

Graph construction can be performed on embedded coordinates or parameter
coordinates depending on the scientific target. For grip layouts, edge weights
should normally use exact quadratic-hypersurface segment lengths from parameter
endpoints:

```r
edge_weights <- gflow::quadform.edge.lengths(
  U = X_param[g$edge_matrix[, 1L], , drop = FALSE],
  V = X_param[g$edge_matrix[, 2L], , drop = FALSE],
  index.k = index_k,
  coefficients = coefficients
)
```

This keeps connectivity and MST repair delegated to `gflow`, while grip
receives only the repaired edge matrix and scientifically meaningful edge
weights.

## gflow Quadform Benchmark Scripts

These scripts live in gflow and are the first sources to reuse or mirror.

2D and mixed 2D benchmark scripts:

- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R`
  - Uses `quadform.sample.dataset()`.
  - Builds graph families with MST repair.
  - Uses `quadform.grid.geodesic.distances()`.
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-curvature-tier1-benchmark/run_quadform_curvature_tier1_benchmark.R`
  - 2D curvature/index sweep.
  - Uses `quadform.embed()` and `quadform.grid.geodesic.distances()`.
  - Uses graph constructors with `connect.components = TRUE` and
    `connect.method = "component.mst"`.
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-tier2-domain-sampling-benchmark/run_quadform_tier2_domain_sampling_benchmark.R`
  - 2D domain-shape and sampling-profile sweep.
  - Contains conventions now mirrored by `quadform.sample.dataset()`.
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-surface-unpruned-benchmark/generate_surface_unpruned_report.R`
  - Uses cKNN graph construction with component-MST repair.
- `/Users/pgajer/current_projects/gflow/dev/quadform-grid-geodesic-calibration/run_quadform_grid_calibration_report.R`
  - Calibrates 2D grid geodesic reference accuracy.

3D benchmark and oracle scripts:

- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-3d-smoke-benchmark/run_quadform_3d_smoke_benchmark.R`
  - 3D ball/cube quadform benchmark.
  - Uses local sampling conventions, `quadform.embed()`,
    `quadform.delaunay.geodesic.distances()`, and MST-repaired graph
    constructors.
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/run_quadform_3d_delaunay_reference_report.R`
  - Focused 3D Delaunay-reference oracle validation.
  - Useful for selecting `n.ref`, `candidate.multiplier`,
    `boundary.fraction`, and `edge.length.factor`.
- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R`
  - Stress tests for the 3D Delaunay oracle.
  - Useful for boundary/interior/mixed placements and boundary-stratified
    errors.

gflow design notes and prompts:

- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/gflowui_quadform_benchmark_agent_prompt.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/gflowui_quadform_benchmark_explorer_brief.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_cpp_backend_decision.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_cpp_regression_fixtures.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_factor4_recommendation.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_agent_prompt.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_brief.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_followup_prompt.md`

Other gflow tools using quadform assets:

- `/Users/pgajer/current_projects/gflow/tools/phate_phase6_generate_baselines.py`
  - Defines `quadform_graph(seed, dim, index, n)`.
  - Provides PHATE phase-6 fixtures for 2D and 3D quadform graphs:
    `quadform_graph_n2_k0`, `quadform_graph_n2_k1`,
    `quadform_graph_n2_k2`, `quadform_graph_n3_k0`,
    `quadform_graph_n3_k1`, `quadform_graph_n3_k2`, and
    `quadform_graph_n3_k3`.
- `/Users/pgajer/current_projects/gflow/tools/phate_phase6_generate_report.R`
  - Reports on the quadform PHATE phase-6 fixtures.

## geodesic_data_geometry Mirror Scripts And Notes

The geodesic-data-geometry project mirrors many gflow benchmark assets. Treat
the gflow copy as canonical when both exist.

Benchmark scripts:

- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-first-benchmark/run_quadform_first_benchmark.R`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-curvature-tier1-benchmark/run_quadform_curvature_tier1_benchmark.R`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-tier2-domain-sampling-benchmark/run_quadform_tier2_domain_sampling_benchmark.R`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-3d-smoke-benchmark/run_quadform_3d_smoke_benchmark.R`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-3d-delaunay-reference/run_quadform_3d_delaunay_reference_report.R`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/quadform-benchmarks/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R`

Design notes:

- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/gflowui_quadform_benchmark_agent_prompt.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/gflowui_quadform_benchmark_explorer_brief.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_cpp_backend_decision.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_cpp_regression_fixtures.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_factor4_recommendation.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_agent_prompt.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_brief.md`
- `/Users/pgajer/current_projects/geodesic_data_geometry/experiments/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_stress_test_followup_prompt.md`

## geodesicMDS And SIMODS Consumers

geodesicMDS experiment scripts:

- `/Users/pgajer/current_projects/geodesicMDS/experiments/quadform-curvature-initialization-sensitivity/run_quadform_curvature_initialization_sensitivity.R`
- `/Users/pgajer/current_projects/geodesicMDS/experiments/quadform-curvature-tier1-gmdsui-bundle/build_quadform_curvature_tier1_gmdsui_bundle.R`
- `/Users/pgajer/current_projects/geodesicMDS/experiments/quadform-initialization-timing/run_quadform_initialization_timing.R`
- `/Users/pgajer/current_projects/geodesicMDS/experiments/edge-isometric-repulsive-unfolding-initial-suite/R/fixtures.R`

geodesicMDS reporting helpers:

- `/Users/pgajer/current_projects/geodesicMDS/notes/scripts/build_quadform_initialization_sensitivity_figures.R`
- `/Users/pgajer/current_projects/geodesicMDS/notes/scripts/build_quadform_initialization_timing_figures.R`
- `/Users/pgajer/current_projects/geodesicMDS/notes/cross_package_function_inventory.md`
- `/Users/pgajer/current_projects/geodesicMDS/notes/geodesic_mds_cross_package_function_inventory.md`

SIMODS manuscript planning and handoff notes:

- `/Users/pgajer/current_projects/geodesicMDS/simods_manuscript/notes/non_oracle_graph_selection_benchmark_plan.md`
- `/Users/pgajer/current_projects/geodesicMDS/simods_manuscript/codex_project/project_status/handoffs/H003_experiment-engineer_non-oracle-selection-benchmark-plan.md`

The SIMODS notes recommend 3D ball-domain smoke cases with:

- `index.k` classes 0 and 1;
- coefficients `(1,1,1)` and `(1,2,4)`;
- `n = 80, 120`;
- `edge.length.factor = 4`;
- `delaunay.backend = "cpp"` when parity remains acceptable.

## gflowui Benchmark Asset Helpers

If existing benchmark assets should be consumed rather than regenerated, use:

```text
/Users/pgajer/current_projects/gflowui/R/quadform_benchmark_helpers.R
```

Relevant helpers:

- `quadform_discover_benchmark_artifacts()`
- `quadform_parse_graph_asset()`
- `quadform_parse_layout_asset()`
- `quadform_generated_layout_cache_path()`
- `quadform_generate_weighted_layout()`
- `quadform_weighted_layout_fun()`
- `quadform_parse_layout_asset_structure()`

UI consumers:

- `/Users/pgajer/current_projects/gflowui/R/app_server.R`
  - Renders the quadform graph stage and quadform benchmark views.

These helpers are useful for report/UI asset reuse. For first-pass grip MISF
edge-KK validation, regenerating small deterministic graphs from exported
gflow functions is cleaner and easier to audit.

## Recommended grip Validation Flow

1. Create the dataset with `gflow::quadform.sample.dataset(dim = 2, ...)` or
   `gflow::quadform.sample.dataset(dim = 3, ...)`.
2. Build local graphs using `create.adaptive.radius.graph()` and
   `create.cknn.graph()` with `connect.components = TRUE` and
   `connect.method = "component.mst"`.
3. Replace constructor Euclidean weights with exact quadform segment lengths
   from `quadform.edge.lengths()` when the validation target is surface
   geometry.
4. Record MST repair diagnostics in every report.
5. Feed only the repaired `edge_matrix` and exact `edge_weight` to grip.

This keeps quadform geometry, graph construction, and MST repair owned by
gflow, while grip remains focused on layout behavior.
