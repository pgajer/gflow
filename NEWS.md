# gflow 0.2.0 (unreleased)

## Bug fixes

* `lcor.with.posterior()` now validates every feature before computation and
  rejects inconsistent draw counts instead of silently truncating later
  features. Inputs must be finite numeric fields with a common positive draw
  count, and the credible level must be strictly between zero and one.
  Supplied feature and vertex names are validated and preserved in results.
  A single draw remains supported, with an undefined (`NA`) standard deviation.

## Further craftsmanship improvements

- Canonical trajectory complexes with both directions can now be passed to
  `gfcor()` and `gfassoc.membership()` with an explicit raw/retained support choice.
  Graph, vertex identity, field and refinement compatibility are checked before
  calling the existing association kernel. Archived input behavior is preserved.
- Local-correlation help now defines its finite-graph edge cosine, anisotropy
  conditions and logratio limitations. The legacy `sign` mode is documented as
  unit weighting. `pairs` computes requested matrix-column pairs in their
  supplied order, allowing bounded blocks without a full association array.
- Seeded permutation tests restore RNG state, support explicit within-group
  permutations, forward statistic settings and print the null/statistic contract.
- R >= 4.1.0 and dgraphs >= 0.2.0 are required. Graph boundary adapters support
  published list and development object APIs; optional viewer diagnostics name
  the required ivue version and installation route. Push/PR checks cover operating
  systems, the minimum R series, and serial/OpenMP builds.
- The guide links its complete catalog to help and opens with four task routes.
  Installed `basin-parameters` help derives defaults from constructor definitions.
  A developer map records implementation ownership, memory costs and interrupt
  limits without reorganizing stable numerical kernels.


* Removed the installed metric-graph low-pass Shiny demo: it called a retired
  estimator that is no longer exported by gflow. Its historical source remains
  in Git history; smoothing belongs to the separate smoothing projects.

## User experience

* Added a complete-help source installation route, a package-help quick start,
  runnable local-correlation examples, and a visual introduction shared by the
  README and installed task guide. Guides now adapt to narrow screens.
* Basin summaries print bounded rankings, diagnostics, and distinct raw,
  retained, and assigned vertex coverage. Explicit plot views show merge trees,
  primary assignments, or overlap counts. Failed results are labeled in field
  views and rejected by analytical views; numerical construction is unchanged.
* Added a task-organized documentation website built from the maintained help
  and vignettes, with the package marked as unreleased.

## User guides

* Added installed task-oriented function and example-graph vignettes, with
  complete export/S3 coverage checks and public migration help.
* Corrected the README constructor call, documented archived-object limits in
  flow-aware association, and made overlap assignment and permutation
  assumptions explicit in the existing workflows.
* Corrected RTCB terminology to relaxed trajectory-constrained basins, matching
  the implemented constrained-path search.
* Earlier extraction milestones below are history, not evidence that removed
  forwarders or experimental APIs remain available.

## Earlier development changes

## Retire duplicate generic 3D renderers

* Removed the internal `plot3D.plain()`, `plot3D.cont()`, `plot3D.cltrs()`,
  `plot3D.tree()`, and `plot3D.path()` implementations. Use `ivue`'s
  `plot3D.plain()`, `plot3D.cont()`, `plot3D.groups()`, `layer3D.edges()`,
  and `layer3D.path()`, respectively.
* Cluster highlights, selected-sample displays, and disk-embedding displays
  now return `ivue` browser widgets. Their additional scene controls use the
  `ivue` API; native additions belong in callback layers before widget capture.
  No compatibility aliases were added. The function-height `plot3D.graph()` was also removed in the graph extraction;
  use documented dgraphs or application plotting code. Domain-specific 2D
  functions remain available as listed in the function guide.

## Separate adaptive extrema from dgraphs

* Renamed `detect.local.extrema()` to `detect.adaptive.extrema()`, returning
  class `gflow_local_extrema` with separate summary, print, plot, and vertex
  methods. The old function name is no longer exported by `gflow`.
* Re-exported `dgraphs::vertices()` as the shared generic so vertex extraction
  works for both packages' objects regardless of attachment order.
* Fixed center inclusion in `vertices(..., include.center = TRUE)` and retained
  the requested maxima/minima setting for summaries of empty results.

## Migrate graph work to dgraphs

* Removed the remaining generic graph construction, conversion, clustering,
  endpoint, path, diagnostic, and plotting exports. Use `dgraphs` where the
  migration guide identifies a verified successor.

## Move retired estimators out of the core package

* Removed PHATE, sparse diffusion/potential pseudotime, random-walk smoothing,
  distance-quantile analysis, quadratic-form geodesics, and their native
  support from `gflow`.
* The archived graph-regression and conditional-expectation estimators remain
  available from `gflowx`; `gflowx` is no longer a `gflow` dependency.

## Use the consolidated analysis surface

* Centered the public package story on canonical basin/flow objects,
  post-construction exploration, local association, and flow-aware
  association.
* De-exported generic weighted-p-value helpers, clustering summaries,
  interactive selection/widgets, and remaining support utilities.
* Added breaking-release migration guidance, now accessible in installed
  `gflow-migration` help and the function guide.

* Moved fitted-model local-slope testing and subject-neighborhood diagnostics
  to `gflowx`; restricted posterior local correlation to supplied field draws;
  and made the package spline wrapper internal.
* Consolidated association entry points around `lcor()`, newly exported
  `lslope()`, `gfcor()`, and their explicit inference layers. Shape-specific
  `lcor` helpers and generic conditional-mean `fassoc*` support are now private.
* De-exported generic statistics, preprocessing, grid, sampling, histogram,
  divergence, Wasserstein, synthetic-data, and miscellaneous plotting helpers.
  Subsequent renderer cleanup removed the generic widget exports; native
  domain-specific extrema-labeling functions remain public. Duplicate matrix
  preprocessing helpers were removed.

* Retired `compute.gfc()`, `compute.basins.of.attraction()`,
  `compute.gfc.trajectory()`, `compute.gfc.flow()`, and `create.basin.cx()`
  with classed migration errors. Added archived `gfc.flow`,
  `basins_of_attraction`, and `basin_cx` conversion methods.
* Migrated canonical basin construction to private geodesic, trajectory, and
  overlap-cell backends so active package and downstream code no longer
  depends on legacy exported constructors.
* Added canonical post-construction basin refinement for relative-value
  filtering, extrema clustering, geometric filtering, support filtering, and
  basin expansion, with explicit per-stage provenance and preserved raw
  membership.
* Removed the temporary generic graph compatibility exports. Use the
  corresponding `dgraphs::` functions directly for graph construction,
  conversion, connected components, paths, distances, endpoint diagnostics,
  graph selection, and stability summaries.
* Added canonical RTCB and overlap-cell-complex adapters, including complete
  RTCB parameter provenance and preservation of merged/unmerged basin, cell,
  overlap, cluster-mapping, graph, and simplified-field structures.
* Added deterministic exact-plateau superlevel and sublevel merge trees to
  `create.basin.complex()`, including component-wise roots, elder ties,
  hierarchical support, branch assignment, merge events, and nonnegative
  persistence.
* Added canonical `trajectory_flow` and `geodesic_reachability` adapters to
  `create.basin.complex()`. The adapters preserve raw overlap separately from
  primary assignment, retain legacy backend objects for provenance, normalize
  set memberships, and isolate seeded tie breaking from the global RNG state.
* Added the Phase B canonical basin-complex API:
  - `create.basin.complex()` now validates graph geometry, fields, mass,
    density, method parameters, and refinement parameters into a stable
    `basin_complex` schema.
  - Added non-colliding `get.basin.*()` table, membership, assignment,
    merge-tree, trajectory, and cell accessors plus print, summary, plot,
    conversion, and data-frame methods.
  - All five method adapters, canonical refinement, downstream migration,
    archived-object conversion, and strict legacy-constructor retirement are
    now implemented. See the canonical basin-complex workflow vignette for
    runnable examples and lifecycle guidance.
* Phase 1 legacy-1D decoupling pass completed:
  - Internal consumers that previously depended on legacy `magelo`/`mabilo`-style
    1D smoothers now use spline-based replacements.
  - `magelo.with.external.BB()` now uses a spline backend while preserving its
    legacy output structure for compatibility.
  - Added consistent deprecation warnings across legacy exported 1D regression
    entry points (`magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`,
    `magelog`, `fit.pwlm*`, and `get.magelo.MAB`).
* At that intermediate stage, legacy 1D model-averaging APIs were marked
  experimental; their extraction and hard removal were subsequently completed.
* Phase 3 extraction started with new `malo` package:
  - `fit.pwlm*` and the 1D `get.*MAB*` benchmarking/model-comparison families now
    have native implementations in `malo`.
  - gflow legacy entry points for these families now delegate to `malo` when
    available (with in-package fallback retained for compatibility).
* Hard removal of legacy 1D `get.*MAB*` and related helpers from `gflow`:
  - `bias_utils` compatibility wrappers were removed from `gflow`.
  - These APIs are now available only from `malo`.
* Hard migration of remaining legacy 1D exported APIs from `gflow` to `malo`:
  - `magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`, `magelog`,
    `magelo.with.external.BB`, `generate.dirichlet.weights`, and `fit.pwlm*`
    are now thin forwarders to `malo` with no in-package fallback.
* Native/header cleanup after 1D migration:
  - Removed legacy 1D `src` implementations and symbol registrations from
    `gflow`.
  - Removed obsolete legacy 1D headers from `inst/include/gflow` and pruned
    stale 1D declarations from `msr2.h`.
* Migration operational cleanup:
  - Updated collaborator installation docs to require installing `malo`
    before using legacy 1D forwarders in `gflow`.
  - CI workflow now installs `malo` explicitly in rdgraph test jobs.
  - Forwarder-focused tests now skip cleanly when `malo` is unavailable.
  - Removed unused legacy 1D header artifacts from
    `inst/include/gflow` (`lm.h`, `ray_agemalo.hpp`, `uggmalog.hpp`).
  - Removed legacy `maelog` implementation and native symbols from `gflow`
    (R API, native sources, headers, and Rd docs).
  - Hard-cut removal of all remaining 1D forwarders from `gflow`:
    `magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`, `magelog`,
    `magelo.with.external.BB`, `generate.dirichlet.weights`, `fit.pwlm`,
    and `fit.pwlm.optimal` are no longer exported by `gflow`.

# gflow 0.1.0 — 2025-09-21

* Historical initial-release milestone; a CRAN publication for this version
  was not found in the current or archived CRAN listings checked on 2026-09-15.
* Implements geometric tools for intrinsic-structure modeling, Morse–Smale regression, and related utilities.
