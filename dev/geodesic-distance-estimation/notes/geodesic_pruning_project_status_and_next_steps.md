# Geodesic Pruning Project Status and Next Steps

Date: 2026-05-06

This note summarizes the local and global geodesic-pruning work completed so
far in the `gflow` R package, the validation that was run, and the recommended
next steps. It is intended to be a handoff document for continuing the pruning
project without having to reconstruct the implementation history from commits
and chat context.

## Executive Summary

The pruning project has completed three implementation tasks and one follow-up
benchmark/parity-report task.

1. The original R local-pruning helper has been moved into C++ for all graph
   families that use the shared pruning path.
2. The new C++ local-pruning helper has been optimized enough for medium
   synthetic graph work.
3. A separately named whole-graph pruning method, `global.geodesic.ratio`, has
   been added without changing the semantics of `local.geodesic`.
4. A parity and runtime benchmark report has been generated for
   `global.geodesic.ratio`.

The current evidence supports the following interpretation:

- The C++ local-pruning helper preserves the original R behavior and now runs
  fast enough for the intended medium benchmark pass.
- The C++ global-ratio helper matches the older R whole-graph reference on the
  tested parity scenarios.
- On medium synthetic runtime benchmarks, `global.geodesic.ratio` was usually
  only modestly slower than `local.geodesic`, but also generally pruned fewer
  edges under the tested settings. Runtime and edge counts alone are not enough
  to decide whether it is geometrically preferable.

The recommended next task is to add downstream graph-geodesic fidelity metrics
so that pruning methods can be compared by geometric quality, not only runtime
and sparsity.

## Relevant Commits

- `076f455 Move local geodesic pruning to C++`
- `71e5401 Optimize local geodesic pruning C++`
- `3543c62 Add global geodesic ratio pruning method`
- `84c6030 Add global ratio pruning parity benchmark report`

## Task 1: Move Local Pruning Fully Into C++

### Goal

The first implementation task was to move the existing R helper
`.prune.graph.local.geodesic()` fully into C++ while preserving the R-level API
and result contract.

### What Changed

The shared C++ helper `S_prune_graph_local_geodesic` was added in
`src/local_geodesic_pruning.cpp` and registered in `src/init.c`. The private R
helper `.prune.graph.local.geodesic()` now validates/coerces inputs and delegates
the actual pruning work to C++.

The implementation preserves the important behavioral details of the old R
algorithm:

- Candidate edges are processed longest first.
- Ties are resolved by smaller `from`, then smaller `to`.
- Local neighborhoods are built using exact Euclidean kNN.
- kNN tie handling is deterministic by vertex index.
- The candidate edge itself is excluded from the alternative local path search.
- The pruning rule is
  `d_L(u, v) <= prune.tau * edge_length + 1e-12`.
- Pruning is sequential: each accepted deletion affects later candidate checks.
- Diagnostics and optional `pruned_edge_stats` keep the existing result shape.

### Graph-Family Routing

The C++ local helper is now used through the existing R dispatch paths for graph
families that support local geodesic pruning, including lifecycle branches such
as `repaired.pruned`.

Relevant call paths include:

- `create.sknn.graph()`
- `create.mknn.graph()`
- `create.radius.graph()`
- `create.adaptive.radius.graph()`
- `create.single.iknn.graph()` for `prune.method = "local.geodesic"`
- `.add.graph.lifecycle.branches()` through `.prune.graph.by.method()`

### Tests Added

Focused tests were added in:

- `tests/testthat/test-local-geodesic-pruning-cpp.R`

The tests cover:

- hand-built examples with known expected edge sets
- no-alternative-path cases
- locality-dependent pruning cases
- sequential pruning behavior
- deterministic tie behavior
- stats-on and stats-off behavior
- malformed graph validation
- graph-family call paths
- lifecycle branches such as `repaired.pruned`

## Task 2: Optimize the C++ Local Pruning Helper

### Goal

The second task was a conservative optimization pass on the C++ local helper,
before using it in broader synthetic graph sweeps. The aim was to improve
runtime without changing the local-pruning algorithm or result contract.

### Profiling Findings Before Optimization

The profiling pass showed that the main costs were not R/C++ call overhead.
They were inside the C++ helper:

- repeated Dijkstra workspace allocation and initialization
- adjacency rebuilding/removal work during sequential pruning
- full sorting of all candidate kNN distances when only the first `k` neighbors
  are needed
- repeated scans of small structures inside the Dijkstra loop

### Optimizations Implemented

The C++ local helper was updated to use:

- a reusable Dijkstra workspace for distance vectors, heaps, and touched nodes
- stable edge IDs with an `alive` flag rather than physically rebuilding the
  graph on every deletion
- adjacency references that point to edge IDs, allowing Dijkstra to skip inactive
  edges
- partial kNN selection with `std::nth_element`, followed by sorting only the
  first `k` entries to preserve deterministic order
- tighter reuse of temporary containers in the pruning loop

These optimizations kept the algorithm conservative: the pruning decision,
candidate order, local neighborhoods, and returned diagnostics remained the
same.

### Benchmark Script

The Task 2 benchmark script is:

- `dev/geodesic-distance-estimation/bench_local_pruning_cpp.R`

Representative optimized timings from the Task 2 wrap-up were:

| kind | n | graph_k | local_k | stats | edges_before | pruned_edges | elapsed_sec |
|---|---:|---:|---:|---|---:|---:|---:|
| line | 1000 | 20 | 20 | FALSE | 10693 | 9694 | 0.094 |
| circle | 1000 | 20 | 20 | FALSE | 10344 | 7062 | 0.108 |
| paraboloid | 1000 | 20 | 20 | FALSE | 11195 | 6165 | 0.134 |
| circle | 2000 | 20 | 20 | FALSE | 21352 | 13401 | 0.301 |
| paraboloid | 2000 | 20 | 20 | FALSE | 22284 | 12227 | 0.337 |
| circle | 2000 | 40 | 40 | FALSE | 41158 | 33065 | 0.864 |
| paraboloid | 2000 | 40 | 40 | FALSE | 44100 | 32681 | 1.194 |
| paraboloid | 2000 | 40 | 80 | TRUE | 44103 | 32767 | 1.474 |

### Validation

After the optimization pass, the focused pruning tests, requested graph tests,
and full test suite passed. The full suite result at that point was:

- `2765` passing tests
- `13` skipped tests

The optimization pass was committed and pushed as:

- `71e5401 Optimize local geodesic pruning C++`

## Task 3: Add a Separately Named Global-Ratio Method

### Goal

The third task was to add a new fast whole-graph pruning variant without
changing `local.geodesic`. The chosen method name was:

- `global.geodesic.ratio`

The intent was to make the behavior explicit and avoid overloading
`local.geodesic` with a global alternative-path interpretation.

### What Changed

A new C++ entry point was added:

- `S_prune_graph_global_geodesic_ratio`

The R helper `.prune.graph.global.geodesic.ratio()` wraps that C++ entry point.
The generic pruning dispatcher `.prune.graph.by.method()` now recognizes the new
method name.

The new method uses the whole graph as the alternative-path search domain. It
uses the same global-ratio predicate shape as the legacy R global helper:

- build candidate set by edge-length percentile
- candidate precheck uses an alternative whole-graph path with the candidate edge
  excluded
- candidate pruning order is longest first, then smaller `from`, then smaller
  `to`
- pruning is sequential and candidates are rechecked against the current graph
- optional `pruned_edge_stats` use the existing result format

### Graph-Family Routing

The new method was routed through:

- `create.sknn.graph()`
- `create.mknn.graph()`
- `create.radius.graph()`
- `create.adaptive.radius.graph()`
- `create.single.iknn.graph()`
- lifecycle branches through `.add.graph.lifecycle.branches()`

For `create.single.iknn.graph()`, `global.geodesic.ratio` is currently a
separately named alias for the existing native iKNN whole-graph global-ratio
pruning path. This was a conservative choice: it gives users the clearer method
name without disturbing the older `global.geodesic` behavior.

### Documentation

Exported function documentation was updated for constructors whose public
argument lists changed. `make document` was run, generating the relevant `.Rd`
updates.

### Tests Added

Focused tests were added in:

- `tests/testthat/test-global-geodesic-ratio-pruning.R`

The initial Task 3 tests covered:

- direct whole-graph shortcut pruning
- no-alternative-path cases
- candidate percentile filtering
- sequential pruning behavior
- deterministic tie order
- stats-on and stats-off result shape
- malformed graph validation
- routing through all graph-family constructors
- lifecycle branches such as `repaired.pruned`

### Validation

Task 3 validation included:

- focused global-ratio pruning tests
- existing graph-family tests
- full `tests/testthat` run
- `make check-fast`

The full suite passed with:

- `2819` passing tests
- `13` skipped tests

`make check-fast` completed with:

- `0` errors
- `0` warnings
- `3` existing-style NOTES

Task 3 was committed and pushed as:

- `3543c62 Add global geodesic ratio pruning method`

## Follow-Up: Parity Tests and Runtime Benchmark Report

### Goal

After Task 3, an additional pass was done to check that the new C++
`global.geodesic.ratio` helper matches the older R whole-graph implementation,
and to benchmark pruning runtime and edge reduction on medium synthetic graphs.

### Parity Tests

The test file `tests/testthat/test-global-geodesic-ratio-pruning.R` was extended
with R-vs-C++ parity checks. These tests compare:

- final edge tables
- edge counts before/after pruning
- number of pruned edges
- optional `pruned_edge_stats`

The parity scenarios vary:

- random seed
- path-edge ratio percentile
- max ratio threshold
- stats on/off

The focused test file passed with:

- `414` passing checks

### Benchmark and Report Script

The benchmark/report script is:

- `dev/geodesic-distance-estimation/bench_global_ratio_pruning_report.R`

The generated report is:

- `dev/geodesic-distance-estimation/reports/global_ratio_pruning_benchmark/index.html`

The generated report artifacts include:

- `raw_graphs.csv`
- `benchmark_replicates.csv`
- `benchmark_summary.csv`
- `parity_audit.csv`
- PNG figures for synthetic data, raw graph density, elapsed time, pruning
  fraction, throughput, and parity results

### Benchmark Design

The benchmark builds raw graphs once per case/family and then times only the
pruning dispatcher on the same raw adjacency and weight lists. This separates
graph construction cost from pruning cost.

Synthetic cases:

- noisy circle, `n = 500`
- noisy circle, `n = 900`
- Swiss-roll-style data, `n = 600`
- three-branch data, `n = 700`

Graph families:

- sKNN
- mKNN
- fixed radius
- adaptive radius
- iKNN

Pruning methods compared:

- `none`
- `local.geodesic`
- `global.geodesic.ratio`

Settings:

- `max.path.edge.ratio.deviation.thld = 0.1`
- `path.edge.ratio.percentile = 0.5`
- `prune.tau = 1.05`
- two timing replicates per graph/method combination

### Benchmark Findings

The independent report parity audit passed:

- `120 / 120` scenarios matched between R and C++

On the medium synthetic runtime benchmark:

- median global/local runtime ratio was about `1.08x`
- `global.geodesic.ratio` was faster than `local.geodesic` in `4 / 20`
  case/family pairs
- median global-minus-local edge-reduction fraction was about `-0.098`
- `global.geodesic.ratio` pruned a larger edge fraction in `4 / 20`
  case/family pairs

The practical interpretation is that `global.geodesic.ratio` is correct and
usable, but under these settings it is often a slightly slower and more
conservative pruning rule than `local.geodesic`. The Swiss-roll-style cases were
the main setting where global-ratio sometimes pruned more or ran slightly
faster.

The benchmark/report pass was committed and pushed as:

- `84c6030 Add global ratio pruning parity benchmark report`

## Current State of the Code

The pruning system currently has three meaningful method names:

- `none`: no geometric pruning
- `local.geodesic`: local-neighborhood weighted-geodesic pruning
- `global.geodesic.ratio`: whole-graph alternative-path ratio pruning

`create.single.iknn.graph()` also retains the historical method name:

- `global.geodesic`

For iKNN, `global.geodesic` and `global.geodesic.ratio` currently reach the same
native global-ratio pruning behavior. This is intentionally conservative.

## Suggested Next Step: Downstream Geodesic-Quality Experiments

### Why This Should Be Next

The runtime and edge-count benchmark answers whether pruning is fast and how much
it sparsifies the graph. It does not answer the more important geometric
question:

Does pruning improve, preserve, or damage graph-geodesic fidelity?

This matters because the downstream package use case is not sparse graph
construction for its own sake. The graph is a geometric estimator. A pruning
method that removes many edges quickly may still be harmful if it distorts
intrinsic distances, disconnects important regions, or changes branch geometry.
Conversely, a method that prunes fewer edges may be better if it preserves
manifold distances while removing shortcuts.

### Proposed Experiment

Add a benchmark/report pass that compares graph-geodesic distances against known
or approximate oracle distances on synthetic geometries.

Recommended synthetic cases:

- noisy circle, using circular arc distance as the oracle
- Swiss roll, using unrolled latent coordinates as the oracle
- branches/roots, using path distance along branch coordinates as the oracle
- optional noisy/outlier variants to test robustness

Recommended graph families:

- sKNN
- mKNN
- fixed radius
- adaptive radius
- iKNN

Recommended pruning methods:

- `none`
- `local.geodesic`
- `global.geodesic.ratio`

Recommended fidelity metrics:

- Spearman correlation between graph-geodesic and oracle distances
- Pearson correlation after optional scale alignment
- RMSE after fitting a single scale factor
- median and upper-quantile relative distance error
- local-neighborhood distance distortion
- disconnected-pair rate
- edge reduction and runtime as supporting metrics

Recommended figures:

- fidelity vs runtime
- fidelity vs edge-reduction fraction
- method ranking by graph family
- error distributions by synthetic geometry
- scatterplots of oracle distance vs graph-geodesic distance for representative
  cases

Recommended report output:

- an HTML report under
  `dev/geodesic-distance-estimation/reports/`
- CSV files for replicate-level metrics and summarized metrics
- a short recommendation section identifying where each pruning method is useful

### Expected Decision Value

This experiment should tell us whether:

- `global.geodesic.ratio` has a geometry-quality advantage that justifies its
  extra cost in some families or data geometries
- `local.geodesic` is the better default candidate because it is faster and at
  least as faithful
- either pruning method is harmful on branching or noisy geometries
- threshold tuning is worth doing before any API/default recommendations

## Additional Suggested Follow-Ups

### 1. Unify iKNN Naming Semantics Later

Right now iKNN keeps the legacy method name `global.geodesic` and adds
`global.geodesic.ratio` as a separately named alias. This was the right
conservative choice for Task 3 because it avoids breaking existing users and
keeps legacy behavior stable.

Eventually, we should decide what name is preferred.

Possible policy:

- keep `global.geodesic` as a legacy-compatible alias
- document `global.geodesic.ratio` as the preferred explicit name
- add tests proving the two names are equivalent for iKNN under the same inputs
- optionally emit no warning for now, to avoid annoying existing users

This does not need to happen before the downstream fidelity benchmark. It is a
semantic cleanup and documentation clarity task.

### 2. Run Downstream Geodesic-Quality Experiments

This is the recommended immediate next step, described above. The central
question is whether local vs global pruning improves graph-geodesic fidelity on
noisy circle, Swiss roll, branches/roots, and related synthetic geometries.

Runtime is only half the story. Pruning can change graph geometry in subtle
ways, especially by removing shortcuts, changing branch-to-branch distances, or
altering local density effects. The next benchmark should measure those effects
directly.

### 3. Clean Up Method-Name Helpers

The helper `.normalize.local.prune.method()` now supports a global method, so
its name is stale. It still works, but it no longer reads honestly.

A later internal cleanup could rename it to:

- `.normalize.graph.prune.method()`

or a similar package-private helper name.

This should be done carefully because it is a shared internal utility. A small
cleanup pass should:

- add the new helper name
- update internal call sites
- optionally leave a compatibility wrapper if useful
- avoid changing public APIs
- rerun graph constructor tests

This cleanup is not urgent and should be kept separate from the downstream
fidelity benchmark so that any behavior changes remain easy to review.

## Recommended Priority Order

1. Run the downstream graph-geodesic fidelity benchmark and report.
2. Use the fidelity results to decide whether threshold tuning is worthwhile.
3. Clarify iKNN naming semantics in docs/tests.
4. Rename stale internal pruning-method helper(s).

The key principle for the next phase is to avoid optimizing or renaming based on
runtime alone. The pruning methods should be judged by their effect on graph
geometry first, with runtime and sparsity as secondary constraints.
