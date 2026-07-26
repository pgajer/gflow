# gflow Cleanup Execution Plan

## Mandate

Shape `gflow` into a focused package whose public API supports:

- exploration and analysis of basin, trajectory, and gradient-flow objects;
- local association methods;
- flow-aware association methods; and
- narrowly scoped adapters needed to apply those analyses to `gflow` objects.

Generic graph infrastructure must be consumed from `dgraphs`. Retired
smoothing and conditional-expectation estimation must remain outside `gflow`.
Application-specific workflows and generic support utilities must not remain
public merely because they are present in the package.

The phases below contain only forward actions, deliverables, and exit gates.

## Non-negotiable protected boundary

Basin construction and its related computational machinery are outside this
cleanup. The cleanup must not merge, rename, de-export, relocate, redesign, or
delete:

- extrema discovery or classification;
- basin, cell, trajectory, or gradient-flow construction;
- basin assembly, filtering, merging, or persistence calculations;
- basin-complex schemas, constructors, validators, adapters, or merge trees;
- native kernels used by any of the above; or
- S3 methods whose primary responsibility is to construct or maintain those
  objects.

Analysis functions that consume an already constructed object may be reviewed,
but their review must not change the construction contract. Any function on the
transitive call path of protected construction code must remain unchanged until
it is explicitly shown that the proposed edit cannot alter construction
results.

Create a machine-readable protected-surface list before changing package code.
Every cleanup commit must pass both a protected-file diff check and the existing
basin/complex boundary tests.

## Ownership labels

Assign exactly one of these labels to every exported function, S3 registration,
dataset, vignette, native symbol, and user-facing workflow:

| Label | Meaning | Required action |
|---|---|---|
| `CORE-ANALYSIS` | Exploration or analysis native to a `gflow` object | Keep public and test as supported API |
| `CORE-ASSOCIATION` | Local or flow-aware association with a clear maintained estimand | Keep public after API consolidation |
| `CORE-PRIVATE` | Required implementation support with no independent user contract | De-export and document as internal |
| `DGRAPHS` | Generic graph construction, conversion, traversal, distance, spectrum, embedding, diagnostics, or plotting | Replace with `dgraphs::` use and remove from the `gflow` API |
| `GFLOWX` | Retired smoothing, conditional-expectation estimation, or an adapter inseparable from those estimators | Relocate to `gflowx` and remove from `gflow` |
| `DOMAIN` | Application- or data-domain-specific workflow | Relocate to a named domain package |
| `EXAMPLE` | Demonstration code rather than maintained package functionality | Move to a vignette or external analysis repository |
| `DELETE` | Dead, duplicated, or unsupported code | Remove after call-site and downstream audits |
| `PROTECTED` | Basin construction or related machinery excluded from this cleanup | Do not change |

No item may be labeled `CORE-ANALYSIS` merely because it uses a graph. No item
may be moved to another package merely to avoid deciding whether it deserves a
maintained API.

## Phase 0: establish the cleanup ledger and guardrails

### Actions

1. Generate an ownership ledger from `NAMESPACE`, roxygen source, registered S3
   methods, `RcppExports`, vignettes, and datasets.
2. Give every ledger row:
   - its ownership label;
   - intended package;
   - public or private status;
   - direct and transitive callers;
   - tests and documentation that cover it;
   - known downstream callers;
   - lifecycle action and target release; and
   - whether it touches the protected construction call graph.
3. Generate the protected-surface list from construction entry points and their
   transitive R and native-code dependencies.
4. Add a check that fails when:
   - a new export or S3 method has no ledger row;
   - a protected file changes in an in-scope cleanup commit;
   - a `DGRAPHS`, `GFLOWX`, or `DOMAIN` item gains new internal callers; or
   - a removed export reappears during documentation generation.
5. Record downstream repositories and scripts that must be searched before an
   API is de-exported or removed.

### Deliverables

- `cleanup/api-ownership.csv`
- `cleanup/protected-basin-surface.txt`
- `cleanup/downstream-repositories.txt`
- an automated namespace/ledger consistency test
- a protected-surface diff check

### Exit gate

Every current export and S3 registration has an owner and disposition, every
protected construction dependency is listed, and no cleanup item remains
classified only as “miscellaneous” or “utility.”

## Phase 1: remove the generic graph compatibility surface

### Actions

1. Search package code, tests, vignettes, documentation, known downstream
   repositories, and analysis scripts for calls to the temporary graph
   wrappers.
2. Replace those calls with their explicit `dgraphs::` equivalents.
3. Remove the following wrappers from the `gflow` namespace and source:
   - `as_igraph()`
   - `graph.connected.components()`
   - `geodesic.core.endpoints()`
   - `graph.geodesic.distances()`
   - `shortest.path()`
   - `summarize.isometry.deviation()`
   - `build.iknn.graphs.and.selectk()`
   - `create.adaptive.radius.graph()`
   - `create.cknn.graph()`
   - `create.cmst.graph()`
   - `create.geodesic.iknn.graph()`
   - `create.iknn.graphs()`
   - `create.iterated.iknn.graphs()`
   - `create.mknn.graph()`
   - `create.mknn.graphs()`
   - `create.radius.graph()`
   - `create.single.iknn.graph()`
   - `create.sknn.graph()`
   - `compute.stability.metrics()`
4. Remove wrapper-specific documentation and tests. Retain boundary tests that
   assert `gflow` calls supported `dgraphs` APIs with compatible inputs.
5. Audit remaining in-scope source for local implementations of generic graph
   conversion, traversal, component, path, distance, spectrum, embedding,
   diagnostic, selection, and plotting operations.
6. Add an overlap check that rejects a new `gflow` export matching a `dgraphs`
   export unless the name is explicitly allowlisted as `PROTECTED`.
7. Leave any graph-like function on the protected construction call graph
   untouched; place it on the protected allowlist for later review outside this
   plan.

### Exit gate

No in-scope `gflow` export is generic graph infrastructure, no in-scope file
contains a second implementation of a public `dgraphs` operation, and any
remaining namespace overlap is explicitly protected.

## Phase 2: enforce the smoothing and archive boundary

### Actions

1. Audit every in-scope reference to `gflowx`, `riem.dcx`,
   `fit.rdgraph.regression()`, and `refit.rdgraph.regression()`.
2. For association code coupled to an archived fitted-model class:
   - move the adapter to `gflowx` when the archived model is essential; or
   - separate a model-independent core calculation from a thin `gflowx`
     adapter.
3. Apply that decision to posterior local-correlation helpers, paired
   local-slope helpers, subject-neighborhood summaries, permutation helpers,
   and their examples and S3 methods.
4. De-export generic spline and smoothing helpers. Keep a helper private only
   when an in-scope core analysis still requires it.
5. Remove smoothing and conditional-expectation language from the in-scope
   public reference index, examples, and vignettes.
6. Do not modify archive hooks that are reachable only through protected basin
   construction or extrema code. Record those hooks as `PROTECTED` rather than
   using this phase to redesign them.
7. Recompute dependency ownership. Retain `gflowx` in `Suggests` only while an
   explicit public adapter or protected method needs it.

### Exit gate

No in-scope exported `gflow` API requires an archived smoother class, and no
in-scope documentation presents smoothing or conditional-expectation
estimation as a responsibility of `gflow`.

## Phase 3: define and consolidate the association API

### Actions

1. Build an estimand matrix for the `fassoc`, `gfassoc`, `gfcor`, `lcor`,
   `lslope`, paired-test, posterior, permutation, and trajectory-association
   families.
2. For each family, document:
   - the scientific estimand;
   - whether locality is defined by a basin, trajectory, flow, graph
     neighborhood, or generic distance;
   - required inputs and accepted object classes;
   - null hypothesis and uncertainty calculation;
   - output class and stable fields; and
   - overlap with another family.
3. Keep an association method public only when its locality or weighting is
   meaningfully local or flow-aware and its contract can be supported.
4. Assign generic graph association to a named owning package or de-export it;
   do not place statistical association infrastructure in `dgraphs` merely
   because it consumes a graph.
5. Select one canonical entry point per estimand. Convert alternate spellings
   and duplicated front ends into lifecycle-marked aliases, then remove them at
   the scheduled breaking release.
6. Standardize common arguments for data, locality, weights, alternatives,
   adjustment, permutation control, missing values, and reproducibility without
   changing the protected construction API.
7. Standardize maintained result objects and their `print()`, `summary()`,
   `coef()`, derivative-extraction, and plotting methods.
8. Make numerical and resampling helpers private.
9. Add tests for input validation, deterministic seeded behavior, edge cases,
   estimand equivalence, S3 structure, and use with representative constructed
   objects.

### Deliverables

- `cleanup/association-estimand-matrix.md`
- a short public association API reference
- lifecycle mappings from aliases to canonical entry points
- focused association tests independent of archived smoother classes

### Exit gate

Each maintained association estimand has one clear public entry point and one
documented result contract. Duplicate tests, aliases, and S3 methods have a
scheduled disposition.

## Phase 4: contract the generic support API

### Actions

Review support-heavy files first, including `stats_utils.R`, `plot_utils.R`,
`synthetic_data_utils.R`, `grids.R`, `random_sampling.R`,
`preprocess_matrix.R`, `hist_utils.R`, `boxcox_utils.R`, `divergences.R`, and
`wasserstein_dist.R`.

For every function in those families, apply this order:

1. If it has no live caller, mark it `DELETE`.
2. If it supports one core public function, make it a private helper beside
   that caller.
3. If it is generic graph infrastructure, assign it to `DGRAPHS`.
4. If it exists only for an archived estimator, assign it to `GFLOWX`.
5. If it is a generic statistical utility with multiple genuine consumers,
   name a suitable owning package before relocation.
6. Otherwise, de-export it; do not maintain a public utilities layer inside
   `gflow`.

Then:

- split large catch-all files by maintained responsibility;
- replace public helper-to-helper coupling with private implementation modules;
- remove duplicated transforms, plotting primitives, samplers, and synthetic
  generators;
- move test fixtures and synthetic examples under `tests/testthat/helper-*` or
  vignette-local code when they are not user APIs; and
- preserve any helper on the protected construction call graph without editing
  it.

### Exit gate

Every remaining public support function has at least one named core workflow,
a stable user contract, and dedicated tests. Catch-all utility files no longer
act as undocumented public subsystems.

## Phase 5: de-export or relocate applied families

### Actions

1. Create an applied-family ownership table covering at least:
   - microbiome and compositional preprocessing;
   - phylotype selection;
   - metabolon two-plane processing;
   - CST-specific colors and summaries;
   - ISA, consensus, and application-specific clustering;
   - PHATE, diffusion-pseudotime, and ordination workflows;
   - MGCP visualization and export;
   - partition heatmaps and report-oriented plots;
   - two-factor and categorical-proportion analyses;
   - concordance and distance-diagnostic reports; and
   - dataset-specific trajectory analyses.
2. Give every family a concrete destination:
   - a named domain package;
   - `gflowx` for frozen legacy material;
   - a vignette or external analysis repository for examples; or
   - deletion.
3. De-export unstable applied functions before relocation unless a coordinated
   downstream move can be completed in the same release.
4. Move each family as a complete closure: source, documentation, S3 methods,
   tests, fixtures, datasets, native symbols, and dependencies.
5. Keep a `gflow` adapter only when interoperability with a maintained core
   object is itself a supported use case. The adapter must not contain the
   application workflow.
6. Remove optional dependencies after their last in-scope caller moves.
7. Add reverse-dependency tests for any domain package that imports a retained
   `gflow` adapter.

### Exit gate

No public API in `gflow` is justified only by one data domain, one paper, one
report, or one analysis pipeline. Every retained adapter is small, documented,
and tied to a core object contract.

## Phase 6: clean the in-scope S3 and namespace surface

### Actions

1. Reconcile registered S3 methods with actual constructors and documented
   classes.
2. Remove orphan methods for relocated or deleted classes.
3. Collapse duplicate `print()` and `summary()` implementations around shared
   private formatters where the output contracts are genuinely the same.
4. Ensure every retained class has:
   - one constructor or clearly documented creation path;
   - a validator;
   - stable fields;
   - concise `print()` output;
   - a useful `summary()` contract when warranted; and
   - conversion methods only where semantics are unambiguous.
5. Remove accidental exports of internal generics and helper methods.
6. Regenerate `NAMESPACE` and Rd files from roxygen; never hand-edit generated
   namespace or documentation files.
7. Exclude protected basin-construction classes and methods from consolidation.

### Exit gate

All in-scope S3 registrations are reachable, documented, tested, and owned by a
maintained class. No relocated class leaves methods or aliases behind.

## Phase 7: prune dependencies and native support

### Actions

1. Map every `Imports`, `Suggests`, `LinkingTo`, and system requirement entry to
   its remaining callers, examples, vignettes, tests, or native symbols.
2. Remove a dependency when its last in-scope caller has moved or disappeared.
3. Place optional features behind `requireNamespace()` and test both dependency
   presence and absence.
4. Audit compiled registrations and headers outside the protected construction
   closure.
5. Move native code with the family that owns it; remove unreferenced symbols
   only after R-call, native-call, and downstream ABI searches.
6. Keep protected construction kernels and their build requirements unchanged.
7. Re-run portability checks for optional OpenMP behavior and all compiled
   in-scope code.

### Exit gate

Every declared dependency and non-protected native symbol has a named owner and
a live tested caller. Package installation does not pull in application or
archive dependencies for core analysis unless they are explicitly requested.

## Phase 8: rebuild the public package story and execute the breaking release

### Actions

1. Rewrite the package-level help, README API map, reference index, and
   vignettes around:
   - consuming constructed basin/flow objects;
   - exploring and summarizing those objects;
   - local association; and
   - flow-aware association.
2. Remove pages, examples, and cross-references for deleted or relocated APIs.
3. Add migration tables for:
   - `gflow` graph wrappers to `dgraphs`;
   - archived estimator adapters to `gflowx`;
   - relocated applied families to their named owners; and
   - consolidated association aliases to canonical functions.
4. Choose the release policy before removing public names:
   - use one announced deprecation cycle when known downstream users need time;
     or
   - make one coordinated pre-1.0 breaking release when all known downstream
     repositories can be updated together.
5. Update all known downstream repositories in the same release window.
6. Publish a concise changelog organized by user migration action, not by
   internal source-file movement.

### Exit gate

A new user can understand the supported package from the README and reference
index without encountering generic graph infrastructure, retired smoothing,
or domain-specific analysis pipelines.

## Validation required after every implementation family

Run these gates after each independently committed family:

1. confirm the protected-surface diff is empty;
2. run the focused tests for the changed family;
3. run the basin/complex boundary tests without modifying their expected
   scientific results;
4. regenerate documentation with `make document`;
5. run the full test suite:

   ```sh
   Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'
   ```

6. run `make check-fast`;
7. run downstream checks for packages or repositories named in the ledger; and
8. commit and push only when the family and boundary gates are green.

Run `make check` at the end of each phase and again from a clean checkout before
the breaking release.

## Final acceptance criteria

The cleanup is complete only when:

- every public function is `CORE-ANALYSIS` or `CORE-ASSOCIATION`;
- every non-public implementation helper is `CORE-PRIVATE` or `PROTECTED`;
- no in-scope generic graph API or implementation remains in `gflow`;
- no in-scope public function depends on retired smoothing or
  conditional-expectation estimators;
- no public function is application-specific without an explicitly approved
  core adapter rationale;
- all in-scope S3 methods, dependencies, native symbols, documentation, and
  tests have named owners;
- protected basin construction and related methods are unchanged; and
- full package, boundary, and downstream checks pass from clean checkouts.

## Decision rule for borderline functions

Keep a borderline function public only if a user can state its stable
`gflow`-specific contract without referring to its current file, historical
caller, or implementation convenience. If that contract cannot be stated, the
default action is de-export first, then relocate or delete after call-site and
downstream review.
