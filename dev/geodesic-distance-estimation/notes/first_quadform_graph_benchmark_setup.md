# First Quadratic-Surface Graph Benchmark Setup

Date: 2026-05-06

This note fixes the first experimental setup for testing data-to-graph
geodesic reconstruction in `gflow`. The purpose is to turn the current design
discussion into an executable benchmark specification: the next step should be
implementation of a runner, not another round of implicit choices.

The benchmark focuses on the **data geodesic geometric reconstruction
problem**. Given a finite sample

$$
X=\{x_1,\ldots,x_n\}\subset \mathbb{R}^3
$$

from a known quadratic graph surface, and a weighted graph construction
\(G(X)\), we ask how close the graph shortest-path metric \(d_G\) is to the
geodesic geometry induced by the surface and to the best finite-sample path
suggested by the true surface geodesic.

## Scientific Question

The primary question is:

> Which data-to-graph construction gives graph shortest-path distances closest
> to the true geodesic geometry of 2D quadratic graph surfaces?

The graph construction includes both topology and edge lengths. All graph
families in this first benchmark use ambient Euclidean edge lengths on the
embedded sample points. The benchmark is not testing graph embeddings yet; the
weighted GRIP layouts in the report are visual diagnostics for the graph
geometry, not the target metric.

The benchmark should distinguish two notions of target distance.

1. **Continuum surface target.** This is the numerical approximation to the
   Riemannian geodesic distance on the underlying quadratic graph surface:

   $$
   d_M(x_i,x_j).
   $$

2. **Finite-sample oracle target.** This is the best sample-realizable path
   length suggested by the true surface geodesic. If \(\gamma_{ij}\) is the
   estimated surface geodesic from \(x_i\) to \(x_j\), sample points near
   \(\gamma_{ij}\) are ordered by their projection arclength along the path and
   ambient embedded distances are summed between consecutive selected sample
   points:

   $$
   d_X^{\mathrm{oracle}}(x_i,x_j)
   =
   \sum_{\ell=0}^{m-1}
   \|x_{a_{\ell+1}}-x_{a_\ell}\|_{\mathbb{R}^3}.
   $$

The continuum target asks whether the graph recovers the surface geometry. The
sample oracle target asks whether the graph is close to the best route that can
be built from the finite observed sample once the true geodesic route is known.

## Datasets

The first benchmark uses two quadratic graph surfaces over the unit parameter
disk

$$
\Omega=\{(u,v)\in\mathbb{R}^2:u^2+v^2\le 1\}.
$$

The graph embedding is

$$
F(u,v)=(u,v,q(u,v)).
$$

The two surfaces are:

$$
q_{\mathrm{paraboloid}}(u,v)=u^2+v^2,
$$

and

$$
q_{\mathrm{saddle}}(u,v)=u^2-v^2.
$$

In `gflow` terms:

| surface | `index.k` | function |
|---|---:|---|
| paraboloid | 2 | \(u^2+v^2\) |
| saddle | 1 | \(u^2-v^2\) |

Initial benchmark scope:

```r
surfaces <- c("paraboloid", "saddle")
n.values <- c(50, 100, 200)
seeds <- 1:5
domain.radius <- 1
sample.method <- "uniform.parameter.disk"
```

The first full runner may include `n = 300` as a heavier confirmation tier, but
the initial implementation should be stable at `n = 50, 100, 200` before the
larger case is added.

Datasets should be generated with:

```r
ds <- quadform.sample.dataset(
  n = n,
  index.k = index.k,
  domain.radius = 1,
  sample.method = "uniform.parameter.disk",
  seed = seed
)
```

The benchmark should use `ds$X_param` for reference-geodesic calculations and
`ds$X_embed` for data graph construction and visualization.

## Reference Distances

Reference distances should be computed with the C++ grid engine:

```r
ref <- quadform.grid.geodesic.distances(
  X = ds$X_param,
  index.k = index.k,
  domain.radius = 1,
  grid.size = 501,
  sample.connection.k = 8,
  oracle = "sample.path",
  oracle.tube.k = 8,
  oracle.tube.radius = NULL,
  return.oracle.paths = FALSE
)
```

This returns:

```r
D.surface <- ref$distances
D.oracle  <- ref$oracle_distances
```

The `501` grid size is based on the calibration report:

```text
dev/quadform-grid-geodesic-calibration/report/quadform_grid_geodesic_calibration_report.html
```

That report compared candidate grids against a 1001-grid numerical reference.
For both the paraboloid and saddle, `grid.size = 501` satisfied conservative
final-report thresholds:

$$
E_F < 0.3\%, \qquad Q_{0.95}(|D_{501}-D_{1001}|/D_{1001}) < 1\%,
\qquad Q_{0.99}<2.5\%.
$$

The surface target and sample oracle target should both be stored in the result
object. The benchmark report should always state which target each metric uses.

## Graph Families

All graph constructors should use embedded sample coordinates:

```r
X.graph <- ds$X_embed
```

The first benchmark includes five graph families.

### Symmetric kNN

```r
create.sknn.graph(
  X.graph,
  k = k,
  prune.method = prune.method,
  prune.tau = 1.05,
  prune.local.k = k,
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

Parameter grid:

```r
k <- 3:10
prune.method <- c("none", "local.geodesic")
```

### Mutual kNN

```r
create.mknn.graph(
  X.graph,
  k = k,
  prune.method = prune.method,
  prune.tau = 1.05,
  prune.local.k = k,
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

Parameter grid:

```r
k <- 3:10
prune.method <- c("none", "local.geodesic")
```

### Iterated kNN

```r
create.single.iknn.graph(
  X.graph,
  k = k,
  prune.method = prune.method,
  prune.tau = 1.05,
  prune.local.k = k,
  threshold.percentile = 0,
  connect.components = TRUE,
  connect.method = "component.mst",
  with.lifecycle.branches = TRUE,
  pca.dim = NULL,
  verbose = FALSE
)
```

Parameter grid:

```r
k <- 3:10
prune.method <- c("none", "local.geodesic")
```

The legacy global iKNN pruning path should not be used in this first benchmark.
The goal here is a parallel local-geodesic pruning comparison across all graph
families.

### Fixed-Radius Graph

For each dataset, compute nonzero pairwise ambient distances on `X.graph`.
Define radius by an expected-degree-like rank:

```r
radius.rank <- 3:10
radius.prob <- radius.rank / (n - 1)
radius <- as.numeric(stats::quantile(
  pairwise.nonzero.distances,
  probs = radius.prob,
  names = FALSE,
  type = 7
))
```

Then build:

```r
create.radius.graph(
  X.graph,
  radius = radius,
  prune.method = prune.method,
  prune.tau = 1.05,
  prune.local.k = radius.rank,
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

Parameter grid:

```r
radius.rank <- 3:10
prune.method <- c("none", "local.geodesic")
```

### Adaptive-Radius Graph

```r
create.adaptive.radius.graph(
  X.graph,
  k.scale = k.scale,
  radius.rule = radius.rule,
  radius.factor = radius.factor,
  prune.method = prune.method,
  prune.tau = 1.05,
  prune.local.k = k.scale,
  connect.components = TRUE,
  connect.method = "component.mst"
)
```

Parameter grid:

```r
k.scale <- 3:10
radius.rule <- c("min", "max")
radius.factor <- c(1, 1.25, 1.5)
prune.method <- c("none", "local.geodesic")
```

## Shared Evaluation Stages

The graph constructors expose lifecycle stages. This benchmark should evaluate
graphs after MST repair, and it should use the same operational order across
all graph families:

```text
raw graph -> MST component repair -> optional local geometric pruning
```

Therefore the evaluation stage is determined by pruning state:

| pruning state | graph stage passed to `graph.geodesic.distances()` |
|---|---|
| `none` | `"raw.repaired"` |
| `local.geodesic` | `"repaired.pruned"` |

Example:

```r
stage <- if (identical(prune.method, "none")) {
  "raw.repaired"
} else {
  "repaired.pruned"
}

D.graph <- graph.geodesic.distances(g, stage = stage)
```

The runner should not use `stage = "final"` for the primary benchmark, because
the primary comparison is explicitly repair-first.

## Metrics

For every dataset and graph setting, compute deviation from isometry against
both reference targets:

```r
surface.summary <- summarize.isometry.deviation(D.graph, D.surface)
oracle.summary  <- summarize.isometry.deviation(D.graph, D.oracle)
```

The existing summary includes:

| metric | meaning |
|---|---|
| `scale` | optimal scalar calibration from graph distances to reference distances |
| `rel_rms_error` | relative RMS distance-matrix error after optional scalar calibration |
| `rel_abs_error_median` | median pairwise relative absolute error |
| `rel_abs_error_q95` | 95% pairwise relative absolute error |
| `distortion_q05` | 5% quantile of calibrated multiplicative distortion |
| `distortion_median` | median calibrated multiplicative distortion |
| `distortion_q95` | 95% quantile of calibrated multiplicative distortion |
| `pearson_cor` | Pearson correlation between graph and reference distances |
| `spearman_cor` | Spearman correlation between graph and reference distances |

The runner should add a `target` column with values:

```r
target <- c("surface", "sample_oracle")
```

This makes downstream tables and plots target-aware.

## Graph Diagnostics

For each graph object, store the following diagnostics when available:

```r
n_vertices
n_edges_in_raw_graph
n_edges_in_raw_repaired_graph
n_edges_in_pruned_repaired_graph
n_edges_in_repaired_pruned_graph
n_components_raw
n_components_raw_repaired
n_components_pruned
n_components_pruned_repaired
n_components_repaired_pruned
n_components_before
n_components_after
n_mst_edges_added
bridge_method
bridge_k
bridge_k_max
bridge_k_used
bridge_exact_fallback_used
n_edges_before_pruning
n_edges_after_pruning
n_pruned_edges
prune_tau
prune_local_k
```

Diagnostics should be stored even when the graph build fails. Failed graph
settings should produce a row with status metadata, error message, and `NA`
metrics rather than stopping the whole benchmark.

## Benchmark Runner Output

The first runner should write a compact but complete artifact set:

```text
results.rds
metrics.csv
graph_diagnostics.csv
dataset_manifest.csv
run_config.json
report.html
```

The recommended directory is:

```text
dev/data-geodesic-reconstruction/quadform-first-benchmark/
```

The `dev/` directory is ignored by git, so report outputs are local artifacts.
Only reusable package code and source scripts should be committed.

## HTML Report Design

The report should be human-facing and self-contained. It should explain:

1. the two surfaces;
2. the two target distances \(D_{\mathrm{surface}}\) and
   \(D_{\mathrm{oracle}}\);
3. why `grid.size = 501` is used;
4. the graph families and parameter grid;
5. the repair-first evaluation stage convention;
6. how scalar calibration and deviation metrics are computed;
7. how to interpret the figures and tables.

Recommended static report sections:

- executive summary with best graph family/parameter by target and surface;
- full metric tables;
- rank plots by `rel_rms_error`;
- parameter heatmaps for each graph family;
- pruning impact plots comparing `"raw.repaired"` and `"repaired.pruned"`;
- connectivity and bridge diagnostics;
- comparison of surface-target and sample-oracle-target rankings.

Recommended interactive 3D section:

- original/input surface sample \(X\);
- weighted GRIP layout of the selected graph;
- optional overlay or side-by-side comparison for selected best/median/worst
  settings.

The 3D report should follow the same general interaction pattern as:

```text
dev/phate-graph-comparison/report/phate_iknn_graph_comparison_report.html
```

but it should not embed every possible widget directly into the page. With the
full graph grid, that would be unnecessarily heavy. Instead, precompute the
payloads and use dropdown selectors such as:

```text
surface
n
seed
graph_family
graph_parameter
pruning_state
target
metric_sort
```

The report should also include quick selectors for:

```text
best
median
worst
```

according to a chosen metric, usually `rel_rms_error`.

## First-Run Scope and Expansion Plan

The recommended first implementation run is:

```r
n.values <- c(50, 100, 200)
seeds <- 1:5
grid.size <- 501
```

Total graph settings per dataset:

| family | configurations |
|---|---:|
| sKNN | \(8 \times 2 = 16\) |
| mKNN | \(8 \times 2 = 16\) |
| iKNN | \(8 \times 2 = 16\) |
| fixed radius | \(8 \times 2 = 16\) |
| adaptive radius | \(8 \times 2 \times 3 \times 2 = 96\) |

Total:

$$
16+16+16+16+96=160
$$

graph configurations per dataset/seed/surface.

For two surfaces, three sample sizes, and five seeds, the first full run has:

$$
2 \times 3 \times 5 \times 160 = 4800
$$

graph evaluations, each with two target summaries. This is a meaningful run and
should be implemented with caching/checkpointing.

Before the full run, the runner should support a smoke configuration:

```r
n.values <- c(50)
seeds <- 1
k.values <- 3:4
radius.rank <- 3:4
k.scale <- 3:4
radius.rule <- c("min", "max")
radius.factor <- c(1, 1.25)
```

The smoke run should produce the same artifact types and a valid HTML report.

## Implementation Checklist

The benchmark runner should proceed in this order:

1. build dataset manifest;
2. generate or load cached `quadform.sample.dataset()` objects;
3. compute or load cached C++ reference distances with
   `quadform.grid.geodesic.distances(..., oracle = "sample.path")`;
4. expand graph parameter grid;
5. for each graph setting, build graph with error capture;
6. extract graph diagnostics;
7. compute graph geodesics for the chosen lifecycle stage;
8. summarize isometry deviation against `D.surface`;
9. summarize isometry deviation against `D.oracle`;
10. save incremental results;
11. compute weighted GRIP layouts for selected report settings;
12. generate the HTML report.

The runner should be resumable. If a dataset/reference/graph result already
exists and the stored configuration hash matches, it should be reused.

