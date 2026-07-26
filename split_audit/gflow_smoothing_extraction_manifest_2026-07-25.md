# Smoothing Extraction Manifest

Date: 2026-07-25

This manifest records the concrete Phase 1 disposition. It is narrower and
more operational than the cleanup action plan.

## Moved to gflowx

### Public estimator families

- rdgraph fit/refit, fitted/residual/summary/print methods, and diagnostics;
- rdgraph permutation test and Bayesian bootstrap;
- data smoother and mean-shift smoother variants;
- ulogit and eigen ulogit;
- Lp response-graph selection and optimal-graph extraction;
- smoothed vertex density based on rdgraph refitting;
- the applied ikNN graph-response pipeline.

### Native implementation

`gflowx` now owns and registers:

- `S_fit_rdgraph_regression`;
- `S_compute_posterior_summary`;
- `S_compute_hop_extremp_radii_batch`;
- `S_ulogit` and `S_eigen_ulogit`;
- the three mean-shift `.Call` entry points;
- the two adaptive mean-shift Rcpp entry points.

Graph construction inside the archived rdgraph wrapper is performed by
`dgraphs::create.single.iknn.graph()`. The archive does not use private `gflow`
functions or symbols.

## Removed from gflow as obsolete duplicate implementation

The following native-only method families had no current R API in `gflow` and
were already owned elsewhere:

- MALO/AGEMALO/PGMALO/PGMALOG/UGGMALO implementation remnants: canonical
  maintained owner is `malo`;
- graph LOWESS, graph kernel smoothing, adaptive spectral smoothing, spectral
  filter, KLAPS low-pass, and standalone harmonic smoothing remnants:
  smoothing-package territory (`geosmooth` or archive), not `gflow`;
- local-complexity and orthogonal-logit remnants coupled to the retired
  regression code.

These were removed rather than copied into `gflowx` because copying inactive
duplicates would create a second archive of code already owned by `malo` or
the smoothing-package split. Historical versions remain recoverable from git.

## Retained in gflow

- harmonic coordinate extension used by flow-complex and trajectory
  exploration;
- basin harmonic repair, extracted into `src/basin_harmonic_repair.cpp` so
  removing the standalone smoother does not alter basin construction;
- graph radius and multi-target shortest-path primitives, extracted into
  `src/graph_radius_neighborhood.cpp`; their future owner is `dgraphs`;
- trajectory spline fitting and its current private support;
- local/flow-aware association consumers of archived fitted objects.
- `compute.pextrema.nbhds()` and its hop-extremp native kernel; this remains
  extrema computation even though it accepts the archived `riem.dcx` class.

The retained items are support for the current basin/flow-complex product
boundary, not a reclassification of general smoothing as core.

`gflowx` currently contains a private copy of the hop-extremp implementation
because it is compiled in the same historical rdgraph translation unit. Its R
entry point is de-exported. Removing that dead compatibility copy is deferred
to the basin/extrema consolidation phase so this split does not refactor
extrema algorithms opportunistically.

## Tests and historical workflows

- focused rdgraph testthat coverage moved to `gflowx/tests/testthat`;
- the rdgraph correctness-report suite moved to
  `gflowx/tests/correctness`;
- rdgraph denoising/manual workflows moved to `gflowx/tests/manual`;
- local-association manual tests remain in `gflow` and explicitly call
  `gflowx` when they need an archived fitted model.

## Verification target

The split should leave:

- no rdgraph/ulogit/mean-shift estimator export or registration in `gflow`;
- no mandatory `gflow` dependency in `gflowx`;
- no direct graph-construction implementation in the archived rdgraph R
  wrapper;
- passing focused estimator tests in `gflowx`;
- unchanged basin/complex tests in `gflow`.
