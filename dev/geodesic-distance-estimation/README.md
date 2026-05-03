# Geodesic Distance Estimation Development Workspace

This directory contains source assets for studying geodesic distance estimation
from sampled data. The current work is graph-based: CMST/MST completion,
iKNN, mKNN, geometric pruning, component repair, oracle truth construction,
and edge/path attribution are treated as candidate tools for recovering useful
intrinsic distance structure.

The project began as a Minimal Spanning Tree Completion Graph development
workspace for `create.cmst.graph()`. That history remains important, but CMST
is now one method family inside a broader geodesic-distance-estimation study.

The files here are intentionally outside the R package build path. They are
meant to live in the GitHub repository as development research material, while
generated reports, figures, caches, and other run outputs stay untracked.

## Layout

- `benchmark_plan.md`: project objective, benchmark ladder, and evaluation
  metrics.
- `cases/`: source scripts for synthetic and semi-synthetic benchmark cases.
- `R/`: reusable benchmark helpers for graph construction, geodesic metrics,
  oracle truth construction, preprocessing, and plotting.
- `reports/`: generated HTML or Markdown reports; ignored by git.
- `figures/`: generated plots and rendered image artifacts; ignored by git.
- `cache/`: temporary simulation and benchmark outputs; ignored by git.

Package-level invariants for `create.cmst.graph()` belong in
`tests/testthat/`. Larger experimental runs and visual correctness studies
belong here.
