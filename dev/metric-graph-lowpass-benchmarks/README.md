# Metric Graph Low-Pass Benchmarks

Development benchmarks and smoke reports for `fit.metric.graph.lowpass()`.

Benchmark result tables use a shared failure-status convention:

- `status = "ok"`: the fit returned and fitted/predicted values were finite.
- `status = "score_error"`: the fit returned, but fitted/predicted values were
  missing or non-finite, so accuracy metrics are not computed.
- `status = "fit_error"`: the fit call threw an error.

Runtime and accuracy summaries should be computed on scored `ok` rows only,
with fit and score-error counts reported beside the summaries.

Render the single-scenario report from any working directory with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/render_single_scenario_smoke.R")
```

The report compares a parameter sweep of `fit.metric.graph.lowpass()` against a
`fit.rdgraph.regression()` baseline on a fixed 1D two-Gaussian mixture example
with `n = 250`.

Compute the full curated multi-scenario fit table with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/compute_multi_scenario_1d_mixture_results.R")
```

or from the shell:

```sh
Rscript /Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/compute_multi_scenario_1d_mixture_results.R
```

Then render the curated multi-scenario benchmark report from the precomputed
fit table with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/render_multi_scenario_1d_mixture.R")
```

The multi-scenario report uses 16 curated two- or three-Gaussian mixture shapes,
3 sampling/noise variants, and 5 random replicates per scenario family.

Run the path-graph multi-method smoke benchmark with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/compute_multi_method_path_1d_mixture_smoke.R")
```

Compute the full path-graph multi-method benchmark with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/compute_multi_method_path_1d_mixture_results.R")
```

Then render the path-graph multi-method benchmark report with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/render_multi_method_path_1d_mixture.R")
```

This benchmark uses the sorted 1D path graph as the primary graph and compares
`fit.metric.graph.lowpass()`, `fit.rdgraph.regression()`, classical 1D
smoothers, `GSD::gsmoothing()`, and `gasper` SGWT/SURE denoising.

Compute the graph trend-filtering curated subset benchmark with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/compute_graph_trend_filtering_subset.R")
```

Then render its report with:

```r
source("/Users/pgajer/current_projects/gflow/dev/metric-graph-lowpass-benchmarks/render_graph_trend_filtering_subset.R")
```

This subset report compares phase-2 `fit.graph.trend.filtering()` on the sorted
path graph against the earlier benchmark winner
`genlasso::trendfilter(..., ord = 2)`. It includes \(k=0,1,2\) graph
trend-filtering fits with unweighted path-graph penalties plus selected
metric-conductance-weighted variants using `weight.rule = "conductance"` and
`"sqrt.conductance"`.
