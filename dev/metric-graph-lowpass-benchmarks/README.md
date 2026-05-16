# Metric Graph Low-Pass Benchmarks

Development benchmarks and smoke reports for `fit.metric.graph.lowpass()`.

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
