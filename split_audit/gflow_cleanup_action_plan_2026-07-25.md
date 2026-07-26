# gflow Cleanup Action Plan

Date: 2026-07-25

## Target package identity

`gflow` should be organized around one scientific object and the analyses that
are native to it:

- basin/gradient-flow complex construction;
- complex, basin, trajectory, and flow exploration;
- local and flow-aware association methods whose estimand is defined by the
  complex or its flow;
- only the graph operations that are specific to those objects.

This is narrower than the package's current historical identity. In
particular, “uses a graph” is not sufficient reason for an implementation to
remain in `gflow`.

## Ownership rules

| Family | Long-term owner | Rule |
|---|---|---|
| Basin/gradient-flow complexes, trajectories, flow-aware summaries | `gflow` | Core |
| Local and flow-aware association defined on those objects | `gflow` | Core, subject to an API review |
| Generic graph construction, traversal, spectra, distances, and repair | `dgraphs` | Move or replace with imports |
| Retired response smoothing and conditional-expectation estimators | `gflowx` | Archive, installable but not CRAN-facing |
| Current maintained geometric smoothers | a dedicated smoothing package, not `gflow` | Already largely outside the current tree; do not reintroduce |
| Application-specific, microbiome-specific, reporting, and one-off workflows | `gflowx` or a domain package | De-export first, then relocate |
| Generic statistical and plotting helpers | lowest-level appropriate package or private helpers | De-export; do not make `gflow` a utilities package |

The existing `/Users/pgajer/current_projects/gflowx` repository is the archive
package for this split. It should contain its own implementations rather than
calling unexported `gflow` internals.

## Phase 1: extract deprecated estimators now

Move the complete implementation and documentation closure for these families
from `gflow` to `gflowx`:

### Rdgraph conditional-expectation family

- `fit.rdgraph.regression()`
- `fit.rdgraph.regression.semiy()`
- `refit.rdgraph.regression()`
- fitted, residual, summary, and print methods for rdgraph fit/refit objects
- `permutation.test.rdgraph()` and `perm.test.audit()`
- `bayes.bootstrap.rdgraph()`
- rdgraph filter/GCV helpers and native implementation

### Older response-smoothing estimators

- `data.smoother()`
- `meanshift.data.smoother()` and the adaptive/GFA variants
- `ulogit()` and `eigen.ulogit()` with their print/predict methods
- required native kernels and private cache/PCA helpers

### Conditional-expectation graph-selection helpers

- `compute.Lp.graphs()`
- `load.Lp.results()`
- `summary.Lp.results()`
- `extract.optimal.Lp.graph()`
- `compute.smoothed.density()`

`gflowx` should use public `dgraphs` graph constructors. The archive must not
retain a private copy of generic graph construction solely to support the old
estimators.

### Applied wrapper coupled to the archived estimator

- `iknn.graph.response.pipeline()` moves with the estimator because it directly
  orchestrates `fit.rdgraph.regression()`.

If its optional visualizations still require `gflow`, that dependency should
be a `Suggests` dependency and isolated behind an explicit runtime check. The
estimator itself must remain usable without loading `gflow`.

## Explicitly not part of Phase 1

### Basin and extrema computation

Do not merge, rename, or otherwise clean up the basin/extrema constructors in
this extraction. This includes the old two-dimensional Morse--Smale utilities,
basin assembly, extrema discovery, basin merging, and their compiled kernels.
They require a separate conceptual consolidation, not opportunistic edits
during package splitting.

Retain `compute.pextrema.nbhds()` in `gflow`. Although it consumes the archived
`riem.dcx` fit class, its estimand and output are extrema neighborhoods. Phase
1 therefore restores it as a standalone extrema API with its own native
hop-radius kernel. The implementation copy bundled with the archived native
rdgraph source is private compatibility code in `gflowx`; remove that duplicate
when Phase 5 establishes the canonical extrema object model.

### Complex and trajectory exploration support

Retain for now:

- harmonic coordinate extension used to extend or compare coordinates on
  graph-flow substructures;
- trajectory spline fitting used by complex/trajectory visualization;
- `gflow.smooth.spline()` only as private support until its callers are
  rationalized.

These are not endorsements of their current public API. The generic spline
wrapper should be de-exported in the support cleanup phase. Harmonic extension
should be revisited after basin/trajectory APIs are consolidated; if it is
primarily a response imputer rather than a complex operation, it should then
move to a smoothing package.

### Local and flow-aware association

Retain local/flow-aware association functions when their neighborhoods,
weights, or estimands are intrinsically defined by a `gflow` complex or
trajectory. Do not move them merely because an earlier inventory grouped all
association tests with application code.

The later review should distinguish:

- core: local association on basins, trajectories, or flow neighborhoods;
- generic: graph association that belongs beside `dgraphs`;
- applied: phenotype-, microbiome-, or report-specific orchestration that
  belongs in `gflowx` or a domain package.

## Phase 2: generic graph infrastructure to dgraphs

This phase should happen after the archive extraction is passing, and before
basin/extrema consolidation.

1. Produce an exact export-by-export overlap report for `gflow` and `dgraphs`.
2. For each overlap, designate `dgraphs` as the implementation owner.
3. Change `gflow` internals to call public `dgraphs` APIs.
4. Keep temporary compatibility wrappers only where downstream usage warrants
   them; mark them with a lifecycle state and a removal release.
5. Remove duplicate native graph builders only after call-site and ABI audits.
6. Add boundary tests asserting that basin/complex results are unchanged when
   graph objects originate in `dgraphs`.

Priority families include kNN/radius graph construction, generic shortest
paths and distances, connected components, spectra/embeddings, graph repair,
and adjacency/weight conversion.

## Phase 3: support-family API contraction

Review every exported helper using three questions:

1. Is the function meaningful to a user without a basin/flow-complex object?
2. Is it used by more than one public core workflow?
3. Is its behavior stable enough to support as an API?

If the answers do not justify a public API:

- first remove `@export` and document it as internal;
- replace cross-file use with a leading-dot private helper where practical;
- relocate generic implementations to `dgraphs` or another owning package;
- delete obsolete helpers only after repository-wide and downstream searches.

Initial candidates include generic spline wrappers, PCA/cache helpers, generic
density/divergence utilities, graph-only plotting utilities, random/synthetic
data helpers, grid helpers, and general statistical transforms.

## Phase 4: applied-family relocation

Clearly label application families before moving them. Candidate labels are:

- `application-microbiome`;
- `application-trajectory-reporting`;
- `application-clustering`;
- `application-visualization`;
- `experimental-association`;
- `legacy-workflow`.

For each family:

1. add a source-level ownership/lifecycle note;
2. de-export functions that are not supported as a stable API;
3. move reusable but non-core workflows to `gflowx`;
4. move domain-specific workflows to a domain package when one exists;
5. keep only thin examples or vignettes in `gflow`, not full analysis
   pipelines.

High-priority candidates are microbiome preprocessing/selection, consensus and
application-specific clustering, one-off trajectory reports, MGCP
visualization/export, concordance reporting, and two-factor application
analysis.

## Phase 5: basin/extrema consolidation

This is intentionally a separate future project. Start only after Phases 1--4
have reduced package noise.

Required deliverables:

- a vocabulary and object-model decision for extrema, basins, cells,
  trajectories, and complexes;
- a call graph of all competing constructors;
- equivalence tests on representative point clouds;
- one canonical constructor pipeline with compatibility adapters;
- removal plan for redundant R and native implementations.

## Phase 1 completion gates

The archive extraction is complete only when:

- `gflowx` installs and loads without importing `gflow`;
- archived estimators call local native symbols under package `gflowx`;
- graph construction is delegated to public `dgraphs` APIs;
- representative fit/refit, smoothing, prediction, and diagnostics tests pass
  in `gflowx`;
- the moved functions and S3 registrations disappear from `gflow`;
- `gflow` documentation no longer advertises conditional-expectation
  estimation as a core responsibility;
- `make document`, focused tests, and `make check-fast` pass in `gflow`;
- equivalent documentation/tests/checks pass in `gflowx`;
- remaining references to moved functions are either updated, explicitly
  optional examples, or listed as blockers.

## Main tradeoff

Making `gflowx` self-contained duplicates a bounded amount of legacy numerical
code and support headers. That is preferable to keeping private cross-package
calls or making the core package depend on its archive. The condition is that
the duplicate code is frozen as archived code: bug fixes may be accepted, but
new graph infrastructure or new estimator development should occur in the
appropriate owning package.
