# Local Geodesic Pruning C++ Implementation Brief

This brief covers background information, assets, implementation context, and
tasks related to moving and optimizing `gflow` local geometric pruning.

## Background

The current data-to-graph geodesic reconstruction work compares graph
shortest-path distances against reference geodesic distances on quadratic graph
surfaces. The benchmark sweeps graph families and parameters:

- symmetric kNN (`create.sknn.graph()`);
- mutual kNN (`create.mknn.graph()`);
- intersection kNN (`create.single.iknn.graph()`);
- fixed-radius graphs (`create.radius.graph()`);
- adaptive-radius graphs (`create.adaptive.radius.graph()`).

For these graph families, local geometric pruning is an experimental opt-in
stage intended to remove locally redundant edges. It is central to the current
benchmark design because it lets us ask whether sparse graph support can retain
or improve geodesic-distance approximation.

The current large benchmark runs showed that local pruning is much slower than
expected, even on moderate sample sizes. This is blocking large parameter sweeps
and interactive report refreshes. The slow behavior appears to be an
implementation issue rather than an inherent algorithmic impossibility: the
iKNN global pruning path is much faster, even though it is conceptually more
global.

## Current Local Pruning Definition

For a candidate edge \(\{u,v\}\) with current edge length \(\ell_{uv}\), the
local pruning rule forms a local vertex set

$$
L(u,v) =
\{u,v\}\cup N_{k_{\mathrm{local}}}(u)\cup N_{k_{\mathrm{local}}}(v).
$$

The direct edge is temporarily excluded, and a shortest alternative path length
is computed inside the induced local graph:

$$
d_{L(u,v)\setminus\{uv\}}(u,v).
$$

For pruning tolerance \(\tau>1\), the edge is removed when

$$
d_{L(u,v)\setminus\{uv\}}(u,v)
\le
\tau\,\ell_{uv}.
$$

Candidate edges are processed longest-first, with vertex-index tie breaking.
This ordering matters because pruning is sequential: once an edge is removed, it
is unavailable as an alternative path for later candidates.

The expected diagnostics are:

- `n_edges_before_pruning`;
- `n_edges_after_pruning`;
- `n_pruned_edges`;
- `prune_tau`;
- `prune_local_k`;
- `with_pruned_edge_stats`;
- optional `pruned_edge_stats` with columns:
  - `u`;
  - `v`;
  - `edge_length`;
  - `alt_path_length`;
  - `path_edge_ratio`.

## Current Implementation State

### R Shared Helper

The shared R implementation is in:

- `/Users/pgajer/current_projects/gflow/R/local_geodesic_pruning.R`

Important helpers:

- `.normalize.local.prune.method()`
- `.normalize.local.prune.controls()`
- `.exact.knn.index()`
- `.local.dijkstra.distance()`
- `.remove.undirected.edge()`
- `.prune.graph.local.geodesic()`
- `.prune.graph.global.geodesic()`
- `.prune.graph.by.method()`

The R helper is algorithmically clear, but slow. In particular,
`.prune.graph.local.geodesic()` repeatedly:

- rebuilds the edge table with `.graph.edge.table()` inside the candidate loop;
- copies and modifies R adjacency lists when removing edges;
- allocates fresh local masks, distance vectors, and visited vectors per
  Dijkstra call;
- computes exact local kNN neighborhoods in R;
- scans R vectors with `which()` and other allocation-heavy operations.

This is the main implementation bottleneck to remove.

### Existing sKNN C++ Pruning

sKNN already has a C++ pruning implementation for the native prune-first path:

- `/Users/pgajer/current_projects/gflow/src/sknn_graphs.cpp`
- `/Users/pgajer/current_projects/gflow/R/sknn_graphs.R`

Important C++ pieces:

- `pruned_edge_t`
- `prune_edges_locally()`
- `local_alternative_path_length()`
- `make_pruned_edge_stats()`
- `S_create_sknn_graph()`

However, the package is not yet fully using C++ local pruning everywhere.
`create.sknn.graph()` still calls `.add.graph.lifecycle.branches()`, and that
helper uses `.prune.graph.by.method()` to generate the `repaired.pruned`
lifecycle branch. Therefore even sKNN still reaches the R helper for part of
the lifecycle output.

### mKNN

mKNN graph construction is in:

- `/Users/pgajer/current_projects/gflow/R/mknn_graphs.R`
- `/Users/pgajer/current_projects/gflow/src/mknn_graphs.cpp`

`create.mknn.graph()` builds the raw graph through `S_create_mknn_graph()`, but
when `prune.method = "local.geodesic"` it calls the R helper
`.prune.graph.local.geodesic()`. Its lifecycle branches also call the R helper
through `.add.graph.lifecycle.branches()`.

### iKNN

iKNN graph construction is in:

- `/Users/pgajer/current_projects/gflow/R/iknn_graph.R`
- `/Users/pgajer/current_projects/gflow/src/iknn_graphs.cpp`

iKNN has historical whole-graph pruning behavior exposed as
`prune.method = "global.geodesic"`. It also supports the newer
`prune.method = "local.geodesic"`, but the local path currently calls the R
helper `.prune.graph.local.geodesic()`.

The iKNN global predicate is a useful model for a future faster pruning variant,
but it must not be conflated with the current local weighted-geodesic pruning
rule.

### Fixed-Radius And Adaptive-Radius Graphs

Radius graph construction is in:

- `/Users/pgajer/current_projects/gflow/R/radius_graphs.R`

Relevant functions:

- `create.radius.graph()`
- `create.adaptive.radius.graph()`
- `.finalize.radius.graph()`

Both graph families currently route pruning through `.prune.graph.by.method()`,
which calls the R local pruning helper when `prune.method = "local.geodesic"`.

### Lifecycle Branches

Graph constructors now expose lifecycle fields:

- `raw_adj_list`, `raw_weight_list`;
- `raw_repaired_adj_list`, `raw_repaired_weight_list`;
- `pruned_adj_list`, `pruned_weight_list`;
- `pruned_repaired_adj_list`, `pruned_repaired_weight_list`;
- `repaired_pruned_adj_list`, `repaired_pruned_weight_list`;
- final `adj_list`, `weight_list`.

The lifecycle logic lives in:

- `/Users/pgajer/current_projects/gflow/R/component_mst_connectivity.R`

Important helper:

- `.add.graph.lifecycle.branches()`

This helper is a major reason the shared pruning implementation needs to be in
C++: every graph family can request a `repaired.pruned` branch, and that branch
currently calls `.prune.graph.by.method()`.

## Main Tasks

### Task 1: Move The R Local-Pruning Helper Fully Into C++

Implement a shared C/C++ callable helper for local weighted-geodesic pruning and
route all graph families through it.

Recommended shape:

```r
.prune.graph.local.geodesic <- function(X,
                                        adj.list,
                                        weight.list,
                                        k,
                                        prune.tau = 1.05,
                                        prune.local.k = NULL,
                                        with.pruned.edge.stats = FALSE) {
    # validate controls in R, then call C++
    .Call("S_prune_graph_local_geodesic", ...)
}
```

Recommended C entry point:

```c
SEXP S_prune_graph_local_geodesic(SEXP s_X,
                                  SEXP s_adj_list,
                                  SEXP s_weight_list,
                                  SEXP s_prune_tau,
                                  SEXP s_prune_local_k,
                                  SEXP s_with_pruned_edge_stats);
```

This entry point should:

1. validate and convert the R adjacency/weight lists;
2. validate that the graph is undirected and weights match both directions, or
   fail with a clear error;
3. compute exact local kNN neighborhoods using Euclidean distance, with the same
   deterministic tie handling as `.exact.knn.index()`;
4. build the initial undirected edge table once;
5. process candidate edges longest-first, using the same ordering as the R
   reference:

   ```r
   edges[order(-edges$weight, edges$from, edges$to), ]
   ```

6. remove an edge only if the currently retained graph contains an alternative
   local path satisfying the \(\tau\)-rule;
7. return the same result contract as the current R helper.

The first C++ implementation should prioritize exact algorithmic parity and
testability over clever optimization. Once parity is proven, optimize.

Important: this is not only for mKNN/radius/adaptive/iKNN. It should also
remove R-helper use from lifecycle branches such as `repaired.pruned`, including
for sKNN.

### Task 2: Optimize The C++ Local Pruning Before Big Sweeps

After the shared C++ implementation is correct, optimize it before using local
pruning in large benchmarks.

Recommended optimizations:

1. No repeated edge-table rebuilds.
   - Build a vector of edge records once.
   - Keep stable edge IDs.
   - Sort candidate edge IDs once.

2. Edge-alive flags.
   - Instead of erasing edges from adjacency vectors, keep `alive[edge_id]`.
   - Dijkstra skips inactive edges.
   - Candidate processing skips inactive edge IDs.

3. Reusable work arrays.
   - Reuse `dist`, heap buffers, predecessor/visited arrays, and local masks.
   - Avoid allocating these per edge.

4. Timestamped local masks.
   - Use integer mark arrays and an incrementing timestamp instead of clearing
     Boolean vectors of length \(n\) for every edge.

5. Candidate filtering.
   - Skip candidates that cannot possibly have an alternative path, such as
     endpoints with no other alive local incident edges.
   - Consider degree/local-degree filters before running Dijkstra.

6. Cutoff-aware Dijkstra.
   - Stop as soon as the smallest queue distance exceeds
     \(\tau \ell_{uv}\).
   - Do not push paths already above cutoff.

7. Efficient local neighborhood construction.
   - Precompute local kNN lists once.
   - Avoid sorting/unique allocations for every edge if possible. Timestamped
     masks can form the local set directly.

8. Optional stats without overhead.
   - Only allocate and store pruned-edge stats when requested.

The implementation should be designed so the unoptimized correctness tests
continue to pass after optimization.

### Task 3: Add A Separately Named Fast Pruning Variant

Add an opt-in faster pruning variant modeled on the iKNN global predicate. This
variant must be clearly named separately from the current local weighted
geodesic pruning.

Do not call this `"local.geodesic"`.

Candidate names to consider:

- `"global.geodesic.ratio"`;
- `"global.shortcut"`;
- `"shortcut.ratio"`;
- `"iknn.global.geodesic"` only if it is intentionally limited to iKNN.

The preferred naming direction is probably `"global.geodesic.ratio"` because it
communicates both scope and criterion.

This future method should be documented as a faster, coarser pruning heuristic,
not an implementation of the current local weighted-geodesic rule.

## Correctness Requirements For Task 1

The C++ implementation must be extensively tested for both algorithmic and
implementation correctness. The goal is not merely that existing tests pass; new
tests should specifically prove that the C++ helper matches the intended local
pruning algorithm across graph families and edge cases.

### Reference Behavior To Preserve

The new C++ helper should match the current R helper on:

- candidate ordering: longest edge first, then smaller `from`, then smaller
  `to`;
- exact local neighborhood construction;
- inclusion of endpoints in the local vertex set;
- exclusion of the direct candidate edge from the alternative path;
- cutoff rule `alt <= prune.tau * edge.length + 1e-12`;
- sequential pruning behavior;
- optional `pruned_edge_stats`;
- no-pruning behavior when no eligible alternative path exists;
- empty-graph behavior;
- component preservation for prune-before-repair stages.

### Unit Test Assets

Existing tests that currently touch local pruning:

- `/Users/pgajer/current_projects/gflow/tests/testthat/test-sknn-graphs.R`
- `/Users/pgajer/current_projects/gflow/tests/testthat/test-mknn-graphs.R`
- `/Users/pgajer/current_projects/gflow/tests/testthat/test-radius-graphs.R`
- `/Users/pgajer/current_projects/gflow/tests/testthat/test-component-mst-connectivity.R`
- `/Users/pgajer/current_projects/gflow/tests/testthat/test-graph-geodesic-distances.R`

Add a new focused test file, recommended:

- `/Users/pgajer/current_projects/gflow/tests/testthat/test-local-geodesic-pruning-cpp.R`

Suggested test cases:

1. Hand-built line with redundant chords.
   - Verify exact remaining edge set.
   - Verify `n_edges_before_pruning`, `n_edges_after_pruning`,
     `n_pruned_edges`.
   - Verify `pruned_edge_stats` values.

2. Square or diamond with diagonal.
   - Verify diagonal is pruned when a two-edge path has equal or near-equal
     length.
   - Verify pruning changes when `prune.tau` is too small.

3. No alternative path.
   - Verify no edge is pruned when the candidate edge is a bridge in the local
     graph.

4. Locality matters.
   - Construct a graph where a valid global alternative path exists but is
     outside the local vertex set.
   - Verify local pruning does not remove the edge.

5. Sequential dependence.
   - Construct a graph where pruning an earlier long edge prevents pruning a
     later edge.
   - Verify C++ reproduces the sequential result.

6. Deterministic tie handling.
   - Use symmetric points with equal distances.
   - Verify local kNN and pruned edge order are deterministic.

7. All graph families.
   - Run the same simple geometry through `create.sknn.graph()`,
     `create.mknn.graph()`, `create.single.iknn.graph()`,
     `create.radius.graph()`, and `create.adaptive.radius.graph()`.
   - Verify local pruning diagnostics and lifecycle fields are populated.

8. Lifecycle branch correctness.
   - Verify `raw`, `raw.repaired`, `pruned`, `pruned.repaired`,
     `repaired.pruned`, and `final` remain selectable through
     `graph.geodesic.distances()`.
   - Verify `repaired.pruned` is genuinely pruning the repaired raw graph.

9. Stats disabled.
   - Verify `with.pruned.edge.stats = FALSE` returns an empty data frame and
     does not alter pruning decisions.

10. Validation failures.
    - Invalid adjacency/weight lengths.
    - Non-positive weights.
    - Out-of-range vertex IDs.
    - Asymmetric adjacency or mismatched reciprocal weights.

11. R-reference parity during transition.
    - Preserve a test-only R reference implementation, or call the old helper
      under a temporary internal name during the port.
    - Compare C++ and R reference on several small random graphs with fixed
      seeds.
    - Remove or keep the reference only as a test helper after parity is proven.

### Test Commands

Targeted tests:

```bash
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-local-geodesic-pruning-cpp.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-sknn-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-mknn-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-radius-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-component-mst-connectivity.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-graph-geodesic-distances.R")'
```

Broader package test:

```bash
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
```

If exported signatures or roxygen blocks change, regenerate docs:

```bash
make document
```

## Implementation Notes

### C++ Placement

Avoid duplicating the pruning implementation separately in every graph family.
A shared C++ source file is preferable, for example:

- `src/local_geodesic_pruning.cpp`
- `src/local_geodesic_pruning.h`

Then use the shared kernel from graph constructors or through a common `.Call`
entry point.

`src/sknn_graphs.cpp` already has useful code, but it is local to that file.
Consider moving reusable pieces into shared files rather than copying them.

### R Wrapper Strategy

Keep the R function `.prune.graph.local.geodesic()` as the package-side entry
point, but make it call C++. That minimizes changes to callers:

- `create.mknn.graph()`;
- `create.single.iknn.graph()`;
- `.finalize.radius.graph()`;
- `.add.graph.lifecycle.branches()`;
- any other path through `.prune.graph.by.method()`.

This approach should automatically move all graph families to the C++
implementation while preserving the existing R-level API.

### Rcpp Or C API

The package already uses both `.Call` C interfaces and Rcpp-generated exports.
For minimum integration risk, a plain `.Call` entry point similar to
`S_create_sknn_graph` is acceptable. If Rcpp is used, update generated exports
in the package's normal workflow.

### Graph Validation

The pruning helper should assume undirected graphs but should validate enough to
catch bad inputs. At minimum:

- `length(adj.list) == length(weight.list) == nrow(X)`;
- every neighbor ID is a valid 1-based integer;
- weights are finite and non-negative, preferably positive for edges;
- each undirected edge appears in both directions;
- reciprocal weights match within a small tolerance.

Validation should fail clearly rather than silently pruning malformed graphs.

### Edge Lengths

Local pruning uses the supplied `weight.list` as graph edge lengths. It should
not recompute graph edge weights from `X`. `X` is only needed for local kNN
neighborhoods.

### Local kNN

The current R helper computes local neighborhoods from Euclidean distances in
`X`, not from graph distances and not from the input adjacency. The C++ helper
must preserve this unless an explicit new method is added later.

### Performance Measurement

After correctness is established, add a small dev benchmark or testthat
performance smoke test outside CRAN-sensitive checks. The target is not a hard
timing assertion, but a reproducible way to compare:

- old R reference;
- first C++ implementation;
- optimized C++ implementation.

Useful benchmark examples:

- line/chord synthetic;
- noisy circle;
- quadratic paraboloid sample with \(n=100\);
- saddle sample with \(n=100\).

## Related Benchmark And Report Assets

Current benchmark specification:

- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/first_quadform_graph_benchmark_setup.md`
- `/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/quadform_benchmark_html_report_spec.md`

Current benchmark runner:

- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/run_quadform_first_benchmark.R`

Current HTML report:

- `/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/quadform-first-benchmark/runs/full/report/quadform_first_benchmark_report.html`

The report currently supports all graph settings through a JSON index, but
weighted-GRIP layout generation for selected graph settings can be slow because
report-only widget generation rebuilds graph objects and can trigger local
pruning. This is another reason to complete and optimize the C++ pruning port.

## Recommended Work Order

1. Add a shared C++ local pruning entry point with exact parity to the R helper.
2. Route `.prune.graph.local.geodesic()` through that C++ entry point.
3. Add focused unit tests for the C++ helper and all graph-family call paths.
4. Run existing graph tests and fix any lifecycle/diagnostic regressions.
5. Only after correctness is established, optimize the C++ kernel.
6. Add the separately named fast global-ratio pruning variant.
7. Update docs and report notes after behavior and naming stabilize.

