# gflowui Quadratic-Surface Benchmark Explorer Brief

This brief covers background, assets, and proposed tasks for extending the
`gflowui` Shiny app so it can interactively inspect graph constructions from
the quadratic-surface data geodesic reconstruction benchmark.

## Background

The current `gflow` work is focused on the data geodesic geometric
reconstruction problem:

Given a finite sample \(X \subset \mathbb{R}^p\) drawn from a structured
geometric object, construct a weighted graph \(G(X)\) whose graph shortest-path
distances approximate reference geodesic distances on the sampled object. In
the first benchmark, \(X\) is sampled from two quadratic graph surfaces:

\[
  z = u^2 + v^2
\]

and

\[
  z = u^2 - v^2.
\]

The benchmark compares data-to-graph constructions including `sknn`, `mknn`,
`iknn`, fixed-radius graphs, and adaptive-radius graphs, with optional
connectivity repair and geometric pruning stages.

The static HTML benchmark report is useful for narrative results, tables, and
figures, but it is not the right tool for full interactive exploration. A
`file://` HTML report cannot safely run `Rscript` to generate missing layouts on
demand. `gflowui`, as a Shiny app, can run server-side R code, cache results,
and display interactive 3D views from in-memory or file-backed assets.

## Current Static Report

The current report is:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/runs/full/report/quadform_first_benchmark_report.html
```

The top-level redirect/index is:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/quadform_first_benchmark_full_report.html
```

The current report specification is:

```text
/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md
```

The report runner is:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R
```

The current static report contains a single interactive 3D row with two panels:

1. original sampled data;
2. weighted GRIP layout of the selected graph stage.

This static version precomputes only a small subset of widget HTML files. The
complete interactive exploration problem is better suited to `gflowui`.

## Relevant gflow Benchmark Assets

Benchmark run directory:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/runs/full
```

Important files:

```text
metrics.csv
graph_diagnostics.csv
dataset_manifest.csv
results.rds
run_config.json
progress.log
report/widget_index.json
```

Dataset/reference caches currently exist under:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/runs/full/cache
```

The benchmark runner has been updated to optionally persist graph and layout
assets during the benchmark pass. New runs can write dataset assets,
stage-specific graph assets containing `adj_list` / `weight_list`, weighted
GRIP 3D layout assets, and consolidated CSV/RDS/JSON manifests for downstream
tools such as `gflowui`.

Relevant command options:

```text
--save.graph.assets=true
--save.layout.assets=true
--asset.stages=benchmark|available|all
```

Use `asset.stages=benchmark` for the primary benchmark stage only, or
`asset.stages=all` when `gflowui` should be able to load all lifecycle stages
without rebuilding graph objects.

## Relevant gflow Package Functions

Graph constructors and distance helpers used in this project include:

```r
gflow::create.sknn.graph()
gflow::create.mknn.graph()
gflow::create.iknn.graph()
gflow::create.radius.graph()
gflow::create.adaptive.radius.graph()
gflow::graph.geodesic.distances()
gflow::quadform.sample.dataset()
gflow::quadform.reference.geodesics()
gflow::summarize.isometry.deviation()
```

The graph constructors now expose explicit lifecycle-style outputs so callers
can inspect multiple stages without rebuilding the graph repeatedly. The stages
used in the benchmark report are:

```text
raw
raw.repaired
pruned
pruned.repaired
repaired.pruned
final
```

The benchmark convention is:

| Pruning option | Primary benchmark stage |
|---|---|
| `none` | `raw.repaired` |
| `local.geodesic` | `repaired.pruned` |

## Relevant gflowui Assets

The `gflowui` repository is:

```text
/Users/pgajer/current_projects/gflowui
```

Important files:

```text
/Users/pgajer/current_projects/gflowui/README.md
/Users/pgajer/current_projects/gflowui/DESCRIPTION
/Users/pgajer/current_projects/gflowui/app.R
/Users/pgajer/current_projects/gflowui/R/app_main.R
/Users/pgajer/current_projects/gflowui/R/app_ui.R
/Users/pgajer/current_projects/gflowui/R/app_server.R
/Users/pgajer/current_projects/gflowui/R/app_server_renderer_helpers.R
/Users/pgajer/current_projects/gflowui/R/app_server_graph_structure_helpers.R
/Users/pgajer/current_projects/gflowui/R/app_server_graph_helpers.R
/Users/pgajer/current_projects/gflowui/R/project_registry_api.R
/Users/pgajer/current_projects/gflowui/docs/gflowui_project_asset_contract.md
/Users/pgajer/current_projects/gflowui/tests/testthat/test-renderer-helpers.R
/Users/pgajer/current_projects/gflowui/tests/testthat/test-graph-structure-helpers.R
/Users/pgajer/current_projects/gflowui/tests/testthat/test-project-registry.R
```

`gflowui` already supports project manifests, graph set selection, and multiple
3D renderers. The README describes three renderer modes:

1. `RGL (live)`: on-the-fly WebGL rendering from in-memory layout data;
2. `HTML`: prebuilt HTML artifact rendering;
3. `Plotly`: reactive Plotly-based 3D rendering.

The existing project asset contract is described here:

```text
/Users/pgajer/current_projects/gflowui/docs/gflowui_project_asset_contract.md
```

The implementation should extend the existing manifest-driven design where
possible rather than creating a disconnected one-off Shiny app.

## Existing gflowui Project-Creation Utilities

`gflowui` already has utilities for building and registering custom projects
from externally generated graph, layout, metadata, and related assets. The
quadratic-surface benchmark explorer should reuse this machinery instead of
building a parallel project system.

Most relevant exported functions:

```r
gflowui::build_project_spec_iknn_3x3()
gflowui::register_project()
gflowui::discover_project_artifacts()
gflowui::list_projects()
gflowui::unregister_project()
```

The `register_project()` function supports explicit custom assets via:

```r
gflowui::register_project(
  project_root = ...,
  project_id = ...,
  project_name = ...,
  profile = "custom",
  scan_results = FALSE,
  graph_sets = ...,
  condexp_sets = ...,
  endpoint_runs = ...,
  defaults = ...,
  overwrite = TRUE
)
```

It also supports a richer project-spec route:

```r
spec <- gflowui::build_project_spec_iknn_3x3(
  project_root = ...,
  graph_sets = ...,
  X = ...,
  defaults = ...,
  metadata = ...
)

gflowui::register_project(
  project_root = ...,
  project_id = ...,
  project_name = ...,
  profile = "iknn_3x3",
  project_spec = spec,
  overwrite = TRUE
)
```

The agent should audit which of these two existing paths is the smallest clean
extension point for the quadratic-surface benchmark. The likely answer is to
add a benchmark-specific project-spec builder or profile only if the existing
custom/project-spec route cannot represent the benchmark controls cleanly.

## Desired gflowui Capability

The desired end-user workflow is:

1. Open `gflowui`.
2. Select the quadratic-surface benchmark project.
3. Choose graph parameters:
   - surface;
   - sample size \(n\);
   - seed;
   - graph family;
   - family-specific parameters;
   - pruning option;
   - lifecycle stage.
4. View the original data and the selected graph layout side by side.
5. If the weighted GRIP layout is already cached, load it immediately.
6. If it is missing, compute it server-side, save it, and then display it.
7. Use benchmark metrics to jump to useful presets such as best-by-surface,
   best-by-sample-oracle, median, or worst cases.

## Recommended UI Design

Add a benchmark-explorer mode or project-specific panel that exposes the
quadratic-surface graph settings directly.

The central visualization should have two 3D panels:

| Left panel | Right panel |
|---|---|
| Original sampled surface/data | Weighted GRIP layout of selected graph/stage |

This two-panel mode is preferable to a single toggle because the scientific
task is visual comparison between the sampled geometry and the graph-induced
layout. A one-panel fallback may be kept for narrow screens or for compatibility
with the existing app shell.

The controls should include:

```text
preset
surface
n
seed
graph_family
prune_method
stage
```

Family-specific controls:

For `sknn`, `mknn`, and `iknn`:

```text
k
```

For `fixed_radius`:

```text
radius_rank
```

For `adaptive_radius`:

```text
k_scale
radius_rule
radius_factor
```

The `preset` control should jump to selected rows from the benchmark metrics,
but it must not hide manual access to any available graph setting.

Recommended presets:

```text
manual
best_by_surface
best_by_sample_oracle
median_by_surface
median_by_sample_oracle
worst_by_surface
worst_by_sample_oracle
```

## Asset Contract Recommendation

The benchmark should produce a manifest that `gflowui` can consume. The exact
shape may be adjusted to match existing `gflowui` patterns, but it should
contain enough information to resolve every graph setting and stage without
guesswork.

Recommended project-level files:

```text
quadform_benchmark_manifest.json
quadform_benchmark_manifest.rds
metrics.csv
graph_diagnostics.csv
dataset_manifest.csv
```

Recommended dataset assets:

```text
datasets/<dataset_id>.rds
```

Each dataset asset should contain at least:

```r
list(
  dataset_id = ...,
  surface = ...,
  n = ...,
  seed = ...,
  X = ...,          # embedded coordinates in R^3
  params = ...      # parameter-domain coordinates if available
)
```

Recommended graph assets:

```text
graphs/<dataset_id>/<setting_id>/<stage>.rds
```

Each graph-stage asset should contain at least:

```r
list(
  dataset_id = ...,
  setting_id = ...,
  stage = ...,
  graph_family = ...,
  parameters = list(...),
  adj_list = ...,
  weight_list = ...,
  diagnostics = list(...)
)
```

Recommended layout assets:

```text
layouts/<dataset_id>/<setting_id>/<stage>/weighted_grip_3d.rds
```

Each layout asset should contain at least:

```r
list(
  dataset_id = ...,
  setting_id = ...,
  stage = ...,
  method = "weighted_grip",
  coords = ...,
  created_at = ...,
  params = list(...)
)
```

## On-Demand Layout Generation

The Shiny app should not rely on the static HTML report to generate missing
layouts. Instead, when a user selects a graph-stage combination with no cached
weighted GRIP layout, the server should:

1. load the graph-stage asset;
2. run weighted GRIP on the selected `adj_list` / `weight_list`;
3. save the layout asset;
4. update the right-hand 3D panel.

The UI should show a clear running state and should avoid blocking the entire
session where practical. If a future async/job mechanism already exists in
`gflowui`, use it. Otherwise, implement the simplest reliable synchronous
version first and leave async execution as a separate follow-up.

## Testing Expectations

The `gflowui` changes should be tested with unit tests and, if practical, a
small app-construction smoke test.

Suggested tests:

1. Manifest parsing:
   - a small synthetic benchmark manifest resolves valid datasets, settings,
     stages, and layout paths.

2. Selector logic:
   - manual selections resolve the exact requested graph row;
   - preset selections map to the expected metric row;
   - family-specific parameter controls are represented in the resolved key.

3. Stage resolution:
   - available stages return the expected graph asset;
   - unavailable stages produce an explicit missing-stage state, not a silent
     fallback.

4. Layout loading:
   - existing layout assets are loaded without recomputation.

5. On-demand layout generation:
   - missing layout assets are generated from a small test graph and saved to
     the expected path.

6. Two-panel data preparation:
   - original-data coordinates and weighted-GRIP coordinates are returned as
     separate renderable objects with matching vertex counts.

7. App construction:
   - the app can be constructed after registering a minimal benchmark project.

## Suggested Implementation Phases

1. Inspect `gflowui` project registry, manifest loading, and existing renderer
   helpers.
2. Define a minimal benchmark manifest extension compatible with the existing
   `gflowui` project contract.
3. Add benchmark-selector resolution helpers with unit tests.
4. Add layout load/generate/cache helpers with unit tests.
5. Add the two-panel 3D viewer UI/server path.
6. Add a small toy benchmark project under tests or dev fixtures.
7. Run focused `gflowui` tests and a local Shiny smoke check.

## Non-Goals For The First Pass

The first pass does not need to:

1. generate all weighted GRIP layouts up front;
2. rewrite the existing `gflowui` workflow around conditional expectation;
3. support every possible future benchmark family;
4. make the static HTML report run background R commands;
5. optimize local geometric pruning.

The important first milestone is a reliable manifest-driven `gflowui` path that
can display original data and one selected graph-stage weighted-GRIP layout
side by side.
