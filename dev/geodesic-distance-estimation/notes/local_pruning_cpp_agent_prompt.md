# Agent Prompt: Move Local Geodesic Pruning Into C++

You are working in the `gflow` R package repository:

`/Users/pgajer/current_projects/gflow`

Before making code changes, read this brief carefully:

`/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/local_pruning_cpp_optimization_brief.md`

Your first task is:

> Move the R local-pruning helper fully into C++ for all graph families.

Focus only on the first task from the brief for this implementation pass. The
performance optimization pass and the new fast global-ratio pruning variant are
separate follow-up tasks.

Implementation goal:

- Add a shared C++ implementation of the current local weighted-geodesic pruning
  algorithm.
- Route the existing R helper `.prune.graph.local.geodesic()` through the C++
  implementation.
- Preserve the current R-level API and result contract.
- Ensure all graph families that currently call `.prune.graph.local.geodesic()`
  or `.prune.graph.by.method()` use the C++ implementation, including lifecycle
  branches such as `repaired.pruned`.

Correctness requirements:

- Candidate edge order must match the current R algorithm: longest edge first,
  then smaller `from`, then smaller `to`.
- Local neighborhoods must match the current exact Euclidean kNN helper,
  including deterministic tie handling by vertex index.
- The candidate edge itself must be excluded from the alternative local path.
- The pruning rule must match:

  \[
  d_{L(u,v)\setminus\{uv\}}(u,v) \le \tau\,\ell_{uv} + 10^{-12}.
  \]

- Sequential pruning behavior must match the current implementation.
- Diagnostics and optional `pruned_edge_stats` must match the existing result
  format.

Testing requirements:

- Add focused unit tests for the C++ local pruning helper, preferably in:

  `tests/testthat/test-local-geodesic-pruning-cpp.R`

- Tests must cover simple hand-built examples with known expected edge sets,
  no-alternative-path cases, locality-dependent cases, sequential pruning,
  deterministic ties, stats-on/stats-off behavior, malformed graph validation,
  and all graph-family call paths.
- Existing graph tests must continue to pass:

```bash
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-sknn-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-mknn-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-radius-graphs.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-component-mst-connectivity.R")'
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-graph-geodesic-distances.R")'
```

After targeted tests pass, run:

```bash
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
```

If exported function documentation or roxygen blocks change, run:

```bash
make document
```

Use conservative engineering judgment: preserve existing behavior first, then
optimize in a later phase.

