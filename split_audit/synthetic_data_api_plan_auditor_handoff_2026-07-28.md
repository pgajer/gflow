# Synthetic Data API Architecture and Implementation Handoff

Status: ready for independent implementation audit and scoped continuation  
Role: architecture-plan author and implementation worker  
Date: 2026-07-28  

Repositories and implemented revisions:

```text
gflow
  worktree: /Users/pgajer/current_projects/gflow
  branch: main
  report-source commit: 478946d7fdb98b89dd426f1077b63c43d4f86dae

geosmooth
  worktree: /Users/pgajer/current_projects/geosmooth
  branch: main
  final implementation commit: 6396af1579f63e48afdd0f26a6e04a9bcdf0d797
  final git status: clean

trend_filtering
  worktree: /Users/pgajer/current_projects/trend_filtering
  branch: main
  consumer-migration commit: 987822c8f27f4fec75b74e16e4c41822cc8a1d85
  final git status: contains unrelated pre-existing README, PDF, LaTeX,
                    bibliography, build-script, and citation-report changes
```

The `gflow` worktree also contains an unrelated modification to `AGENTS.md`.
It was neither staged nor modified in this work.

## Goal

Address the two nonblocking comments in:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_third_reaudit_2026-07-28.md`

and implement the accepted synthetic-data architecture far enough to establish
the canonical contract, absorb G1--G7 and the shared SSRHE regression suite,
migrate the primary SSRHE consumer, and create explicit ownership and
continuation records for the remaining statistical families. Basin and extrema
construction remained outside scope.

## Third Re-audit Disposition

The third re-audit accepted the architecture and left two nonblocking comments.
The implementation now validates every path-specific policy-derived seed before
RNG state changes or any draw is consumed. The report was rebuilt after its
canonical source was committed, and the rendered revision block records
`gflow` `478946d7` and `geosmooth` `6396af1`.

The finding response is:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_third_reaudit_response_2026-07-28.md`

## Work Completed

### Canonical geosmooth framework

`geosmooth` commit `232e412` added:

- serializable geometry, sampling, truth, and response components;
- `synthetic.spec` and `synthetic_dataset` contracts;
- canonical XDR/SHA-256 specification and content identity;
- locale-independent canonical ordering and normalization of incidental
  in-memory reference sharing;
- `named.stream.v1` and pinned legacy RNG policies;
- pre-draw validation of base and policy-derived seeds;
- `materialize.synthetic()` and frozen-instance materialization;
- analytic quadform embedding, gradient, metric, and segment length;
- sphere-cap, helix, torus-patch, simplex, and stratified geometry;
- thin deprecated G1--G7 argument translators; and
- a 16-case frozen G1--G7 parity fixture.

Arbitrary truth closures and G7 `zero.parts = integer(0)` are documented
breaking exclusions rather than silently changed compatibility behavior.

### Shared SSRHE suite

`geosmooth` commits `2ffa515` and `b9beaf4` added:

- all 16 S-shape rows and three V mechanisms;
- 48 S/V recipes;
- three maintained flat recipes;
- 45 maintained quadform recipes over intrinsic dimensions 2--4;
- 106 registry recipes in total after including G1--G7;
- an 846-case parity fixture covering small and ordinary sizes, three
  replicates, and the flat/quadform lanes;
- normalized recipe, geometry, sampling, truth, response, instance, and
  checksum tables;
- registry regeneration and foreign-key/hash tests; and
- byte-for-byte idempotent registry regeneration.

`trend_filtering` commit `987822c` converted
`development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R` into
thin adapters over the `geosmooth` registry. The helper no longer loads
`gflow` or contains copied dataset-generation algorithms.

### Remaining maintained statistical families

`geosmooth` commits `12e981a` and `6396af1` added:

- a finite-design occupation-mixture truth with explicit design-maximum
  normalization;
- draw-free circle and trefoil geometries;
- deterministic interval-grid sampling that preserves the historical circle
  and trefoil endpoint conventions;
- a common `plot.synthetic_dataset()` method;
- a user vignette; and
- a 16-row family disposition ledger.

The ledger is:

`/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/legacy_family_disposition.csv`

Its status counts at the implemented revision are:

```text
implemented                 5
implemented-core            3
implemented-g5-subset       1
implemented-noise-free      1
partial-private             1
approved-pending            3
pending-decision            1
out-of-scope                1
```

The partial and pending rows are admissions, not retirement authorizations.

## Files Changed or Created

Primary `geosmooth` source:

```text
DESCRIPTION
NAMESPACE
R/dgp_library.R
R/synthetic_dataset.R
R/synthetic_geometry.R
R/synthetic_materialize.R
R/synthetic_registry.R
R/synthetic_response.R
R/synthetic_sampling.R
R/synthetic_spec.R
R/synthetic_truth.R
```

Registry and documentation:

```text
inst/synthetic_registry/checksums.csv
inst/synthetic_registry/gaussian_mixture_1d_shapes.csv
inst/synthetic_registry/gaussian_mixture_1d_variants.csv
inst/synthetic_registry/geometries.csv
inst/synthetic_registry/instances.csv
inst/synthetic_registry/legacy_family_disposition.csv
inst/synthetic_registry/recipes.csv
inst/synthetic_registry/registry_manifest.md
inst/synthetic_registry/responses.csv
inst/synthetic_registry/samplings.csv
inst/synthetic_registry/truths.csv
scripts/update_synthetic_registry.R
vignettes/synthetic-datasets.Rmd
```

Parity fixtures and tests:

```text
tests/fixtures/g1_g7_legacy_parity_v1.rds
tests/fixtures/ssrhe_legacy_parity_v1.rds
tests/testthat/test-synthetic-api.R
tests/testthat/test-synthetic-curves.R
tests/testthat/test-synthetic-geometry-quadform.R
tests/testthat/test-synthetic-guardrails.R
tests/testthat/test-synthetic-legacy-parity.R
tests/testthat/test-synthetic-legacy-ssrhe.R
tests/testthat/test-synthetic-occupation.R
tests/testthat/test-synthetic-registry.R
```

Consumer source:

```text
/Users/pgajer/current_projects/trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R
```

Canonical report and audit records:

```text
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_audit_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_audit_response_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_reaudit_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_reaudit_response_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_second_reaudit_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_second_reaudit_response_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_third_reaudit_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_architecture_third_reaudit_response_2026-07-28.md
```

## Generated Artifacts

The canonical report source is:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`

The generated self-contained report is:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`

Final report facts:

```text
Rmd
  lines: 2,376
  bytes: 99,220
  SHA-256: 576a89d8c896139e49bc4784b7eb5371111f9789589dc698becf81d9002c1f8b

HTML
  lines: 4,375
  bytes: 1,363,357
  SHA-256: 2da67ae41ad09c2415a8dfe5bb3b8e4f71fd9a6aa9c6e801a34caf696ed461d8

Third re-audit response
  lines: 59
  bytes: 2,802
  SHA-256: 3d8505f1262f9ce02335eafc2ae2e8192bf20f36d6cb8f2c2cb6a2845ad9e67e
```

The HTML displays build time `2026-07-28 22:28:29 EDT`.

## Commands Run

Commands were run from `/Users/pgajer/current_projects/geosmooth` unless a
different directory is stated.

Documentation and focused tests:

```sh
make document

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  testthat::test_file("tests/testthat/test-synthetic-occupation.R")'

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  testthat::test_file("tests/testthat/test-synthetic-curves.R")'

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  files <- Sys.glob("tests/testthat/test-synthetic-*.R");
  for (f in files) testthat::test_file(f, stop_on_failure = TRUE)'

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  testthat::test_file(
    "tests/testthat/test-synthetic-guardrails.R",
    stop_on_failure = TRUE);
  rmarkdown::render(
    "vignettes/synthetic-datasets.Rmd",
    output_file = tempfile(fileext = ".html"),
    quiet = TRUE)'
```

Registry regeneration was run from the same directory:

```sh
Rscript scripts/update_synthetic_registry.R
```

The resulting registry files were compared byte-for-byte with their pre-run
copies before commit.

Package QA:

```sh
make check-fast
_R_CHECK_FORCE_SUGGESTS_=false make check-fast
```

Report render from `/Users/pgajer/current_projects/gflow`:

```sh
Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file =
    "synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html",
  quiet = TRUE
)'
```

A temporary non-self-contained render was also completed for asset validation.
Artifact checks used `rg`, `wc`, `shasum -a 256`, `git diff --check`, and
repository status/log commands.

## Validation

- All synthetic-focused test files passed after the final synthetic source
  changes, including G1--G7, SSRHE, registry, guardrail, occupation, curve, and
  quadform tests.
- The G1--G7 fixture contains 16 frozen cases.
- The shared SSRHE fixture contains 846 frozen cases.
- Registry regeneration is byte-for-byte idempotent at the committed revision.
- The normalized registry contains 106 recipes and no unresolved component
  foreign keys.
- The synthetic vignette renders successfully.
- The primary `trend_filtering` adapter was parsed, loaded, and exercised
  against representative S16.V3 and quadform recipes with exact round trips.
- The default `make check-fast` stopped at its dependency gate because
  `genlasso` is unavailable.
- `_R_CHECK_FORCE_SUGGESTS_=false make check-fast` completed installation,
  namespace, load/unload, S3 registration, source syntax, dependency, compiled
  code, vignette, and vignette-rebuild checks.
- The forced-Suggests check ended with two warnings: the non-mainstream local
  dependencies and an older `normalize.density` usage/documentation mismatch.
  Its five notes concern bundled ANN object files, `.github`, nonstandard
  `audit_artifacts` and `tmp`, an unimported `tail`, and an older lost-brace Rd
  diagnostic. Unavailable `genlasso` is additionally reported as INFO when
  forced Suggests are disabled.
- The final self-contained and non-self-contained report renders completed.
- Structural HTML checks confirmed the implementation section, accepted seed
  rule, locale-independent canonical ordering, reference-sharing
  normalization, current repository revisions, and 26 acceptance criteria.
- The report source was not modified after the successful final render. Only
  this handoff was updated afterward.

## Canonical and Generated File Notes

- The Rmd report is canonical; the HTML report is generated by
  `rmarkdown::render()` and was not edited by hand.
- `NAMESPACE` is generated by roxygen from R source. The repository ignores
  generated `man/*.Rd`; those files were regenerated locally but are not part
  of the commits.
- `recipes.csv`, component tables, and `checksums.csv` are generated by
  `scripts/update_synthetic_registry.R`.
- The S01--S16 and V1--V3 catalog CSVs, `instances.csv`,
  `legacy_family_disposition.csv`, and `registry_manifest.md` are maintained
  source assets.
- The two RDS fixtures are frozen parity evidence and should not be regenerated
  without an explicit compatibility-contract change and provenance update.
- No package source was changed after its final focused tests or package QA.

## Limitations and Unverified Claims

- Phases 6--9 are not complete. Noisy predictor curves, the mixed
  Gaussian/line sampler, star-response assembly, spline and bump truths,
  secondary consumers with distinct seed policies, and the two ambiguous
  torus-knot formulas remain partial, pending, or undecided.
- The deterministic `synthetic.point.line.junction()` implementation does not
  establish parity for the stochastic legacy `generate.mixed.points()`.
- The current clustered sampler establishes the maintained G5 subset, not the
  full unequal-size, supplied-mean, general-covariance legacy helper.
- Only the primary SSRHE consumer helper was migrated. Other research scripts
  use different historical seed formulas and were not silently redirected.
- Old de-exported support remains in `gflow`; this implementation does not
  claim those files are safe to delete.
- Generic graph topology for a future star-response dataset has not been
  implemented in `geosmooth`; its intended owner remains `dgraphs`.
- Exact same-environment fixtures were generated and tested on the current
  macOS/R environment. Cross-platform scientific tolerance has not been
  calibrated on independent CI systems.
- A clean default package check was not obtained because `genlasso` is absent.
  The forced-Suggests run still has two warnings and five notes, including
  pre-existing source-package and documentation issues.
- The full package test suite has two known failures outside the synthetic
  tranche: a CSD4 PCA-support call-count expectation and a PTTF path requiring
  unavailable `genlasso`.
- Direct interactive visual QA of the final local HTML was not performed:
  the in-app browser security policy rejected the `file://` page. The prior
  report revision had responsive QA at 1,440, 768, and 390 pixels; this revision
  changes prose and hashes but not the report CSS or table structure.
- Temporary non-self-contained report assets are under the operating system
  temporary directory and may be removed automatically.
- The unrelated dirty files in `trend_filtering` and the unrelated `gflow`
  `AGENTS.md` modification were not validated or committed as part of this
  work.

## Reusable Workflow Capture

Classification: script/template already present.

Rationale: registry normalization and checksum regeneration are order-sensitive
and reusable. `scripts/update_synthetic_registry.R` is the committed
regeneration script, while package tests enforce foreign keys, hashes,
disposition statuses, wrapper thinness, and the absence of cross-repository
source calls.

No additional skill or workflow note was created.

## Current State

The accepted architecture, implementation commits, generated report, parity
fixtures, registry, family disposition ledger, validation evidence, and
limitations are ready for independent implementation audit. The remaining
pending families require scoped continuation rather than bulk transfer or
deletion.
