# gflow 0.2.0 downstream breaking-release audit

Audit date: 2026-07-26.

This audit distinguishes migrations that can be made against a verified
successor from experimental APIs for which no maintained successor currently
exists. The latter must not be changed by mechanical namespace substitution.

## Migrated in this release window

| Repository | Change | Status |
|---|---|---|
| `AGP` | `gflow::detect.graph.endpoints()` changed to the verified `dgraphs::detect.graph.endpoints()` owner | Updated locally; this directory is not a Git repository |
| `CT_clearance` | Replaced `gflow::cluster.graph.louvain()` and `gflow::congruence.with.labels()` with script-local `igraph` Louvain membership and a tested Hubert-Arabie ARI calculation | Committed and pushed as `6faf8c3` |

## Calls that cannot yet be migrated safely

| Repository | Qualified calls | Family | Reason |
|---|---:|---|---|
| `ZB` | 14 | PHATE and sparse diffusion/potential pseudotime | No verified maintained successor; worktree also contains extensive unrelated changes |
| `geodesicMDS` | 2 | quadratic-form embedding/sample fixtures | No verified maintained successor; experiment worktree is dirty |
| `geodesic_data_geometry` | 2 | quadratic-form Delaunay reference geodesics | No verified maintained successor; experiment worktree is dirty |
| `geosmooth` | 2 | quadratic-form migration fixtures | Fixtures explicitly compare against the removed implementation; no replacement fixture owner has been published |
| `graph_modularity` | 4 | `phate.core()` | No verified maintained successor; worktree is dirty |
| `gflow-w2` | 2 | interactive `select3D` manual drivers | Active secondary `gflow` worktree explicitly excluded from edits by the downstream registry |
| `gflowx` | 3 documentation references plus private runtime lookup | retired `plot3D.*.widget` UI | `gflowui` does not yet provide the named successor API; changing the lookup would invent a nonexistent contract |

## Acceptance consequence

The `gflow` package-local acceptance gates pass, including the namespace,
S3, dependency, native, protected-surface, source-test, example, vignette,
manual, and full package checks. The cleanup plan's final requirement that
all downstream checks pass is **not yet satisfied** for the unresolved rows
above.

Those repositories must either:

1. adopt a subsequently verified successor;
2. vendor an explicitly owned experimental implementation; or
3. pin the analysis environment to `gflow` 0.1.x.

Until one of those actions is selected per family, replacing the calls would
silently change scientific behavior and is therefore not an acceptable
migration.
