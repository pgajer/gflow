# gflow 0.2.0 breaking-release migration guide

Version 0.2.0 is a coordinated pre-1.0 breaking release. The package no longer
keeps compatibility wrappers for APIs owned by another package or for
experimental families that lack a maintained owner.

## Generic graph APIs

| Removed from `gflow` | Migration |
|---|---|
| Generic graph construction, adjacency/edge conversion, paths, distances, components, spectra, embeddings, diagnostics, graph selection, and graph plotting wrappers retired in Phases 2–6 | Call the same canonical API from `dgraphs`; see `dgraphs/GRAPH_MIGRATION_ACTION_PLAN.md` and its ownership ledger |
| `detect.graph.endpoints()` | `dgraphs::detect.graph.endpoints()` |
| `cluster.graph.louvain()`, `congruence.with.labels()` | No `gflow` compatibility wrapper. Use `igraph` clustering and an explicitly chosen label-agreement implementation until a verified `dgraphs` successor is published |
| `find.shortest.paths.within.radius()`, `compute.local.distance.fidelity()`, `create.pruned.isize.list()`, `create.vertex.labels()`, `create.vertex.label.indicators()` | No supported `gflow` replacement; migrate generic graph logic to `dgraphs` or application code |
| `triangle.plot()`, `itriangle.plot()`, `plot3D.graph()` and graph-degree/grid demonstration helpers | Use `dgraphs` plotting/diagnostic APIs or application plotting code |

## Smoothing and experimental geometry

| Removed from `gflow` | Migration |
|---|---|
| `fit.rdgraph.regression()`, `refit.rdgraph.regression()`, graph-response smoothers, conditional-expectation estimators | Use the archive package `gflowx` when reproducing legacy analyses |
| `weighted.p.value()`, `weighted.p.value.summary()` | Choose and document an application-level inferential procedure; these generic helpers have no core adapter |
| `compute.diffusion.pseudotime.sparse()`, `compute.potential.pseudotime.sparse()` | Removed experimental family; no maintained successor is claimed |
| `phate()`, `phate.core()`, `phate.embed()` | Use a dedicated maintained PHATE implementation; no `gflow` successor is claimed |
| `quadform.*` geodesic/embedding functions | Removed experimental family; no maintained successor is claimed |
| `fit.rho.randomwalk()`, `distance.quantile.bin.analysis()`, `E.geodesic.X()`, `estimate_density()`, `o.inv.fn()`, `vars.approx.monotonically.assoc.with.geodesic()` | Removed generic smoothing/geometry experiments; no core adapter |

## Applied and UI families

| Removed/de-exported from `gflow` | Owner or next action |
|---|---|
| Interactive `select3D.*` and `plot3D.*.widget` workflows | Relocate to `gflowui`; no public `gflow` adapter |
| Compositional/microbiome helpers | Relocate to `microbiome.utils` |
| Two-factor/CST analyses | Relocate to `gcstflow` or the owning analysis repository |
| Metabolon preprocessing and diagnostics | Relocate to the owning analysis repository |
| MGCP and partition heatmap visualization | Relocate to `gflowui` |
| Clustering summaries and DBSCAN representative extraction | Keep in the consuming analysis repository; they are no longer public `gflow` APIs |

## Association consolidation

| Historical entry point/family | Canonical API |
|---|---|
| Shape-specific local-correlation helpers | `lcor()` |
| Local directed-response aliases | `lslope()` or `lslope.neighborhood()` |
| Fitted-model posterior correlation adapters | Supply field draws to `lcor.with.posterior()`; fit legacy models in `gflowx` |
| Generic `fassoc*` conditional-mean tests | No public core adapter; use `gfcor()` only for the documented flow-aware estimand |
| Flow-membership/polarity/overlap/deviation summaries | `gfassoc.membership()`, `gfassoc.polarity()`, `gfassoc.overlap()`, `gfassoc.deviation()` |

## Release policy

The 0.2.0 release removes public names in one step because `gflow` remains
pre-1.0 and the known downstream repositories are inventoried in
`downstream-repositories.txt`. Downstream call sites must be updated to a
verified successor or pinned to 0.1.x when they depend on an experimental
family for which this guide states that no maintained successor exists.
