# Graph Trend Filtering Development Project

This directory is the development home for adding graph trend filtering to
`gflow`.

The package implementation should live in the normal R package locations:

- `R/graph_trend_filtering.R` for the R API and validation.
- `src/graph_trend_filtering.cpp` for computational kernels.
- `src/graph_trend_filtering.h` or `inst/include/` only if shared C++ helpers
  become useful.
- `tests/testthat/test-graph-trend-filtering.R` for package-level regression
  tests.

This `dev/graph-trend-filtering/` folder is for research and implementation
scaffolding that should remain separate from the package API:

- `papers/`: local literature notes and optional ignored PDF cache.
- `scripts/`: exploratory scripts, smoke benchmarks, and report builders.
- `reports/`: generated HTML/PDF benchmark reports.
- `examples/`: small runnable examples for path graphs and weighted graphs.
- `fixtures/`: tiny deterministic graphs used by exploratory scripts.

## Initial Target

The first graph trend filtering comparator should probably start with weighted
graph fused lasso:

```math
\hat\beta_\lambda =
\arg\min_{\beta \in \mathbb R^n}
\frac{1}{2}\|y-\beta\|_2^2 + \lambda \|D_w \beta\|_1,
```

where \(D_w\) is a weighted oriented incidence matrix. This is the graph
\(k=0\) case and is the cleanest place to settle graph validation, edge
weighting conventions, solver behavior, tuning, and return objects.

Higher-order graph trend filtering can then follow the Wang--Sharpnack--Smola--
Tibshirani construction:

```math
\hat\beta_\lambda =
\arg\min_{\beta \in \mathbb R^n}
\frac{1}{2}\|y-\beta\|_2^2 + \lambda \|\Delta_w^{(k+1)}\beta\|_1.
```

The main open design question is the weighted operator convention, especially
whether metric edge lengths should enter as conductances, inverse lengths, or
another scale in \(D_w\).

## Candidate R API Names

Possible entry points:

- `fit.graph.trend.filtering()`
- `refit.graph.trend.filtering()`
- `graph.trend.filtering.operator()`

The API should be a comparator to `fit.metric.graph.lowpass()` and
`fit.rdgraph.regression()`, not a replacement for either.

## Literature Starting Points

See `papers/README.md` for links to the relevant graph trend filtering papers
and `literature_notes.md` for implementation-facing notes.

## Report Build Metadata

LaTeX reports in this project should use generated build metadata rather than
hard-coded timestamp commands. The shared Codex convention is documented in
`/Users/pgajer/.codex/notes/report_build_conventions.md`.
