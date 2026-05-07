# Quadratic-Surface Benchmark HTML Report Specification

Date: 2026-05-06

This document specifies the human-facing HTML report for the first
quadratic-surface data-to-graph geodesic reconstruction benchmark. It refines
the report requirements in
`dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
and fixes the interactive 3D visualization design.

The report has two distinct layers:

1. **Evaluation layer.** This layer compares graph geodesic distances with
   reference distances and ranks graph settings.
2. **Graph visualization layer.** This layer shows a selected graph construction
   and graph lifecycle stage in 3D.

These layers must remain separate. In particular, the reference target
(`surface` or `sample_oracle`) is used only for evaluation and preset
selection. It is not a graph-construction parameter.

## Report Purpose

The report should answer:

> For samples from the paraboloid and saddle quadratic graph surfaces, which
> data-to-graph constructions produce graph shortest-path distances closest to
> the reference geodesic geometry?

The report should also support visual inspection of any graph setting generated
by the benchmark, not only the best setting.

## Required Artifacts

The report run directory should contain:

```text
results.rds
metrics.csv
graph_diagnostics.csv
dataset_manifest.csv
run_config.json
report/quadform_first_benchmark_report.html
report/figures/*.png
report/widgets/*.html
report/widget_cache/*.rds        # optional, for cached weighted GRIP layouts
```

The top-level HTML file

```text
quadform_first_benchmark_full_report.html
```

may be a redirect/index into `runs/full/report/quadform_first_benchmark_report.html`.
The report itself should live inside the `report/` directory so that relative
paths to figures and widgets are stable.

## Report Structure

The report should read as a research report, not a table dump. Figures should
precede tables whenever they summarize the same result.

Required sections:

1. **Context and Question**
   - Define the data geodesic geometric reconstruction problem.
   - State that graph edges use ambient Euclidean lengths in the embedded
     sample.
   - Identify the two surfaces:
     \[
     q(u,v)=u^2+v^2
     \]
     and
     \[
     q(u,v)=u^2-v^2.
     \]

2. **Evaluation Targets**
   - Define the continuum surface target:
     \[
     D_{\mathrm{surface}} = \left(d_M(x_i,x_j)\right)_{ij}.
     \]
   - Define the finite-sample oracle target:
     \[
     D_{\mathrm{oracle}} =
     \left(d_X^{\mathrm{oracle}}(x_i,x_j)\right)_{ij}.
     \]
   - Explicitly state that targets affect ranking and presets only; they do not
     change graph construction.

3. **Metric Definitions**
   - Define scalar calibration.
   - Define relative RMS error.
   - Define the relative absolute error quantiles and distortion quantiles used
     in the tables.

4. **Execution Summary**
   - Configured datasets.
   - Datasets represented in current results.
   - Number of graph settings represented.
   - Metric rows.
   - Successful metric rows.
   - Error metric rows.
   - Whether an elapsed-time watchdog was used.

5. **Sample-Oracle Target Results**
   - Figure: method performance by graph family across `n`.
   - Figure: distribution of relative RMS error by graph family.
   - Table: best settings by surface and `n`.
   - Table: family-level median performance.

6. **Surface Target Results**
   - Same figure/table structure as sample-oracle target.

7. **Pruning and Connectivity Diagnostics**
   - Figure: pruning effect by target and graph family.
   - Figure: MST bridge and pruning diagnostics.
   - Table: diagnostic summaries.

8. **Interactive 3D Graph Diagnostics**
   - One interactive 3D row with controls.
   - The row should show one selected graph setting at a time.

9. **Artifacts and Reproducibility**
   - List CSV/RDS/JSON artifacts.
   - Show command used to generate the report.

## Definitions: Successful and Error Rows

A **successful metric row** is one row in `metrics.csv` with `status = "ok"`.
It means:

1. graph construction completed;
2. the selected graph lifecycle stage was available;
3. graph shortest-path distances were computed;
4. the graph distance matrix was compared with one evaluation target.

Each graph setting normally contributes two successful metric rows, one for
`surface` and one for `sample_oracle`.

An **error metric row** has `status = "error"`. It is retained for bookkeeping
and should include the graph setting, target, and error message. If an
elapsed-time watchdog is used, timeouts must be reported as computational
timeouts, not as mathematical graph failures.

## Graph Construction Controls

The 3D visualization controls should expose the graph construction layer
directly. They should not hide all but the best graph.

Required controls:

```text
preset
surface
n
seed
graph_family
pruning
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

The UI should show only controls relevant to the selected graph family, or keep
irrelevant controls disabled.

## Preset Control

The `preset` control is a convenience selector. It should not remove manual
access to any tested graph.

Required preset values:

```text
manual
best_by_surface
best_by_sample_oracle
median_by_surface
worst_by_surface
```

Optional preset values:

```text
median_by_sample_oracle
worst_by_sample_oracle
```

When `preset = manual`, the user can select every graph setting represented in
the benchmark results.

When a non-manual preset is selected, the controls should jump to the
corresponding graph setting. After the jump, the user should still be able to
change any construction parameter manually.

## Stage Control

`stage` selects the graph lifecycle version to visualize or evaluate.

Allowed values:

```text
raw
raw.repaired
pruned
pruned.repaired
repaired.pruned
final
```

The benchmark convention is:

| pruning | primary benchmark stage |
|---|---|
| `none` | `raw.repaired` |
| `local.geodesic` | `repaired.pruned` |

The report should default to this benchmark convention when a graph setting is
selected. The user may override `stage` manually to inspect repair/pruning
order effects.

If a selected stage is unavailable for a graph object, the visualization panel
should display a clear message rather than silently falling back.

## 3D Visualization Row

The report should contain exactly one interactive 3D row. The row should be
conceptually similar to the interactive rows in
`dev/phate-graph-comparison/report/phate_iknn_graph_comparison_report.html`.

Recommended layout:

```text
controls
metrics/metadata box
3D panel row
```

The 3D panel row should have either two or three panels.

Minimum two-panel design:

1. **Original data**
   - Points shown at the embedded surface coordinates
     \[
     (u_i,v_i,q(u_i,v_i)).
     \]
   - Points rendered as spheres.

2. **Weighted GRIP layout**
   - Weighted GRIP embedding of the selected graph and selected stage.
   - Points rendered as spheres.
   - Edges drawn when the edge count is small enough, or subsampled for dense
     graphs.

Preferred three-panel design:

1. **Original data**
2. **Selected graph on original coordinates**
3. **Weighted GRIP layout**

The three-panel design is better for interpretation because it separates graph
support geometry from graph-layout geometry.

## 3D Widget Computation Model

The report is static HTML, so it cannot run R/GRIP computation in the browser.
Therefore, weighted GRIP layouts must be precomputed or generated by an
external report-refresh step.

Recommended implementation:

1. Precompute widgets for:
   - all small settings, such as `n = 50`;
   - all preset-selected settings;
   - any setting already cached from a previous view.

2. For non-precomputed manual selections:
   - show original data and graph-on-original-coordinate panels if available;
   - show a clear message in the weighted-GRIP panel:
     `Weighted GRIP layout not precomputed for this setting.`

3. Provide an R helper/report mode that accepts a selected graph key, computes
   the missing weighted GRIP layout, caches it, and regenerates the report.

The layout cache key should include all fields that change the graph stage:

```text
surface
n
seed
graph_family
k
radius_rank
k_scale
radius_rule
radius_factor
pruning
stage
layout_algorithm
layout_parameters
gflow_version_or_git_sha
```

## Data Model For Interactive 3D Selection

The report should embed a compact JSON index with one row per graph setting
represented in the benchmark.

Required fields:

```text
graph_key
dataset_id
surface
n
seed
graph_family
k
radius_rank
k_scale
radius_rule
radius_factor
pruning
default_stage
available_stages
surface_rel_rms_error
sample_oracle_rel_rms_error
n_vertices
n_edges_stage
n_components_stage
n_mst_edges_added
n_pruned_edges
original_widget
graph_original_widget
grip_widget
grip_cached
```

`graph_key` should be stable and deterministic. It should be safe as a filename
stem after replacing punctuation with underscores.

## Metrics/Metadata Box

The 3D section should include a small metrics box below the controls and above
the panels.

Required fields:

```text
graph family
selected parameters
pruning
stage
surface rel_rms_error
sample-oracle rel_rms_error
n_edges
n_components
n_mst_edges_added
n_pruned_edges
```

If the selected graph has no metric for one target, show `NA`.

## On-Demand Layout Refresh

The static report cannot compute missing layouts by itself. The recommended
workflow is:

1. User selects an uncached graph in the report.
2. Report displays the graph key and a message that the weighted-GRIP layout is
   missing.
3. A companion R script can be run with that graph key:

```bash
Rscript dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R \
  --mode=full \
  --report.only=true \
  --layouts=true \
  --layout.graph.key=<graph_key>
```

4. The script computes the missing weighted GRIP layout, updates the widget
   cache, and regenerates the report.

This keeps the report static and reproducible while still supporting manual
exploration of every tested graph setting.

## Default 3D View

The initial 3D view should use:

```text
preset = best_by_surface
surface = paraboloid
n = smallest available n
seed = first available seed
graph_family = adaptive_radius if available, otherwise first available family
stage = benchmark default for that setting
```

This is only the initial view. Manual controls must expose all represented graph
settings.

## Implementation Notes

- Use relative paths from the report HTML to figures/widgets.
- Do not URL-encode absolute local paths into `src` attributes.
- Prefer a single visible widget row and JavaScript controls that switch iframe
  sources.
- Keep generated widget artifacts under `report/widgets/`.
- Keep reusable source scripts and specifications under version control; do not
  commit large generated widget/report output unless explicitly requested.

