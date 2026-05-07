# Prompt For The gflowui Agent

Please work in the `gflowui` repository:

```text
/Users/pgajer/current_projects/gflowui
```

Before making changes, read this brief:

```text
/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/gflowui_quadform_benchmark_explorer_brief.md
```

Your task is to extend `gflowui` so it can serve as an interactive explorer for
the `gflow` quadratic-surface data geodesic reconstruction benchmark.

First, inspect the current `gflowui` project registry, manifest loading,
renderer helpers, and tests. Then generate a phased action plan that covers:

1. how to reuse or minimally extend the existing `gflowui` project creation
   utilities, especially `build_project_spec_iknn_3x3()` and
   `register_project(profile = "custom" / "iknn_3x3")`;
2. the benchmark manifest/project asset contract needed by `gflowui`;
3. selector logic for surface, `n`, seed, graph family, graph parameters,
   pruning option, and lifecycle stage;
4. a two-panel 3D viewer showing original data and weighted GRIP graph layout
   side by side;
5. loading cached weighted GRIP layouts;
6. generating and caching missing weighted GRIP layouts server-side;
7. unit and smoke tests needed to prove correctness.

Do not implement the plan until the plan has been reviewed and approved. After
approval, execute the plan with high attention to correctness and keep changes
well scoped to the existing `gflowui` architecture.

Testing expectations after implementation:

```r
pkgload::load_all(".", quiet = TRUE)
testthat::test_dir("tests/testthat")
```

If exported APIs or roxygen blocks are changed, regenerate documentation using
the package's existing documentation workflow.
