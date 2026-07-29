# Synthetic Data Registry Schema Correction Handoff

Status: ready for independent implementation re-audit

Role: implementation worker responding to an independent re-audit

Date: 2026-07-29

Repositories and revisions:

```text
geosmooth
  worktree: /Users/pgajer/current_projects/geosmooth
  branch: main
  base commit:
    819996d43bba7b25980d073dce046e64cb9d3bed
  correction commit:
    a950b27b5ae73f34ea937a718a250f1aac9c8523
  remote state:
    pushed to origin/main
  final git status:
    clean

trend_filtering
  worktree: /Users/pgajer/current_projects/trend_filtering
  branch: main
  unchanged consumer revision:
    2062282986bd58d36c82c4154159f7df95059b83
  final git status:
    unrelated pre-existing README, PDF, LaTeX, bibliography,
    build-script, and citation-report changes remain

gflow
  worktree: /Users/pgajer/current_projects/gflow
  branch: main
  audit-response/report-source commit:
    d561ab1e1e51eebdfc50270313cc6f518ff8dd66
  rendered-report commit:
    c0576731844da477f0f14557f168e804d61ffad5
  final status before this handoff was added:
    unrelated pre-existing modification to AGENTS.md
```

The commit containing this document is the delivery commit for the handoff
itself. The code, report-source, and rendered-report revisions are fixed above
so their provenance does not depend on a self-referential handoff hash.

## Goal

Address ART-R1, QA-R1, and QA-R2 in:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_2026-07-29.md`

while preserving the accepted architecture in:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`

The work did not include basin/extrema construction or later synthetic
families recorded as partial, pending, or undecided.

## Work Completed

### JSON component and recipe records

The component ledgers no longer store R-XDR payloads. Each component row now
contains typed JSON `parameters`. The private JSON codec preserves R atomic
types, names, lists, `NULL`, typed missing values, full-precision doubles, and
matrices with declared dimensions and row-major nested values.

The normalized resolver reconstructs each component from JSON, derives its
versioned class from component kind and family, and verifies its canonical
component SHA-256. Recipe compatibility and metadata are also JSON. Dataset
content checksums continue to use the separate accepted
`synthetic-content-v1` / `R-xdr-v3` contract.

### Environment-fingerprint normalization and completeness

`checksums.csv` now references
`environment_fingerprints.csv` through an exact
`environment.fingerprint.id` foreign key. The normative required fields are
stored separately in `environment_fingerprint_schema.csv`.

The fingerprint records:

- full R build, platform, architecture, compiler, and endianness;
- OS name, release, version, and machine;
- RNG policy, all three RNG kinds, and the explicitly established RNG version;
- absolute BLAS and LAPACK paths and file SHA-256 values;
- LAPACK version;
- linked math-runtime identity, including `libR` path and digest;
- `geosmooth`, registry-schema, and evaluator-registry versions; and
- dependency versions.

Each environment row has a full SHA-256 and an ID derived from that digest.
Runtime matching verifies the row identity before comparing every normative
field.

### Deterministic verification-policy coverage

Production frozen-instance verification was factored into a private helper.
The registry suite passes the live fingerprint to exercise the exact branch,
then changes the copied BLAS digest and requires the same verifier to take the
scientific-parity branch. No public option or process-wide verification bypass
was added.

### Grouped test-runner exit status

`scripts/run_test_group.R` now calls each `testthat::test_file()` with
`stop_on_failure = TRUE`. The real `make test` lane now stops and exits
nonzero at the known coupled-KD expectation rather than returning false
success. A source contract test protects this behavior.

### Report consistency

The canonical report's implementation summary now describes JSON components,
normalized environment foreign keys, and the complete scope fields. Its
normative registry schema and implementation summary no longer describe
different designs.

The finding response is:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_response_2026-07-29.md`

## Files Changed Or Created

### geosmooth canonical source

```text
/Users/pgajer/current_projects/geosmooth/DESCRIPTION
/Users/pgajer/current_projects/geosmooth/R/synthetic_materialize.R
/Users/pgajer/current_projects/geosmooth/R/synthetic_registry.R
/Users/pgajer/current_projects/geosmooth/scripts/run_test_group.R
/Users/pgajer/current_projects/geosmooth/scripts/update_synthetic_registry.R
/Users/pgajer/current_projects/geosmooth/tests/testthat/test-synthetic-registry.R
/Users/pgajer/current_projects/geosmooth/tests/testthat/test-test-runner-contract.R
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/environment_fingerprint_schema.csv
```

### geosmooth generated registry artifacts

```text
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/checksums.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/environment_fingerprints.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/geometries.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/recipes.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/responses.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/samplings.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/truths.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/registry_manifest.md
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/fixtures/
```

The fixture directory contains the ten regenerated G-family frozen-instance
RDS files.

### gflow audit and report artifacts

```text
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_2026-07-29.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_response_2026-07-29.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_second_reaudit_handoff_2026-07-29.md
```

No `trend_filtering` source was changed in this tranche.

## Generated Artifacts

The registry was regenerated twice. The binary Git-diff digest for all changed
`inst/synthetic_registry` artifacts was identical before and after the second
run:

```text
5eb86f627ea1e369751750b4ce47b2fb084d4cb020424a133d446c08ae03ed6e
```

Report artifact measurements after the final render:

```text
canonical R Markdown
  lines:    2,431
  bytes:  102,652
  SHA-256:
    0ca8747a3ffc2f7eaa19df82cd73e13cc0dc8dc129ad3722b95471a20e8cb70d

self-contained HTML
  lines:      4,435
  bytes:  1,367,013
  SHA-256:
    808fa81058fc618faaeb98ff5add263d33a92b4f4aad7fa8219eb2ffc045ff60

implementation re-audit
  lines:      491
  bytes:   19,559
  SHA-256:
    20aac371a9c2cc2f3e97d7577a2664531f78e73f4ed6dde08c6e3610fcf87ab0

implementation re-audit response
  lines:      200
  bytes:    7,954
  SHA-256:
    6aa8b6eeb4610099bfded55a337f02b7718c95c8129cd2caa439850cec829a62
```

## Commands Run

From `/Users/pgajer/current_projects/geosmooth`:

```sh
Rscript -e '
  pkgload::load_all(".", quiet=TRUE)
  # Round-trip every component parameter, compatibility, and metadata value
  # for all 106 source-owned recipes through the typed JSON codec.
  # Construct and inspect the complete environment fingerprint.
'

Rscript scripts/update_synthetic_registry.R

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  testthat::test_file(
    "tests/testthat/test-synthetic-registry.R",
    reporter="summary", stop_on_failure=TRUE)'

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  files <- c(
    Sys.glob("tests/testthat/test-synthetic-*.R"),
    "tests/testthat/test-test-runner-contract.R"
  );
  for (f in files)
    testthat::test_file(
      f, reporter="summary", stop_on_failure=TRUE)'

git diff --binary -- inst/synthetic_registry | shasum -a 256
Rscript scripts/update_synthetic_registry.R
git diff --binary -- inst/synthetic_registry | shasum -a 256

make test
make document
_R_CHECK_FORCE_SUGGESTS_=false make check-fast
make clean
git diff --check
```

From `/Users/pgajer/current_projects/trend_filtering`:

```sh
Rscript development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
```

From `/Users/pgajer/current_projects/gflow`:

```sh
Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file =
    "synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html",
  quiet = TRUE
)'

wc -l -c \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html \
  split_audit/synthetic_data_api_plan_implementation_reaudit_2026-07-29.md \
  split_audit/synthetic_data_api_plan_implementation_reaudit_response_2026-07-29.md

shasum -a 256 \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html \
  split_audit/synthetic_data_api_plan_implementation_reaudit_2026-07-29.md \
  split_audit/synthetic_data_api_plan_implementation_reaudit_response_2026-07-29.md
```

## Validation

All nine `test-synthetic-*.R` files passed serially after the schema
correction. This includes the G1--G7 and full SSRHE frozen scientific-parity
suites. The registry suite reconstructed all 106 recipes and validated every
component, recipe, instance, checksum, fixture, environment-fingerprint, and
normative-schema relationship. The controlled exact and scientific-fallback
branches both passed. The grouped-runner contract test passed.

The SSRHE consumer integration returned:

```text
SSRHE synthetic dataset adapter identity integration: PASS
```

The real grouped smoke lane returned nonzero as intended:

```text
make test
  failing source:
    tests/testthat/test-coupled-kd-selection-csd4.R:183
  counted$calls: 76
  expected:       3
  Rscript failure propagated to make exit 2
```

This is evidence that the runner correction works. It is not a passing smoke
suite.

`make document` completed successfully. With forced Suggests disabled,
`make check-fast` completed with two warnings and five notes. The reported
items were the existing non-mainstream dependency warning, existing
documentation warning, missing optional `genlasso` information, source and
repository artifact notes, an existing unimported `tail()` note, and an
existing Rd formatting note. Package installation, dependency declarations
for `jsonlite`, namespace loading, R syntax, S3 checks, code checks, compiled
code, installed files, and vignettes completed.

## Canonical/Generated File Notes

- `R/synthetic_registry.R`, `R/synthetic_materialize.R`, the registry update
  script, test runner, tests, `DESCRIPTION`, and
  `environment_fingerprint_schema.csv` are canonical inputs.
- The component, recipe, checksum, environment-fingerprint, and instance
  registry CSV files and frozen RDS fixtures are generated by
  `scripts/update_synthetic_registry.R`.
- XDR remains canonical only for dataset content checksums and RDS fixture
  storage; component and recipe registry values are JSON.
- Roxygen outputs were regenerated with `make document`. No new public R
  function was exported.
- The R Markdown report is canonical. Its sibling HTML was rendered from the
  committed source revision and was not edited by hand.
- No Markdown note under `/Users/pgajer/.codex/notes` was changed.

No `geosmooth` source was modified after the passing complete synthetic suite.
The second registry generation was byte-identical to the first. No
`trend_filtering` source was modified. The canonical R Markdown was committed
before the final HTML render and was not modified afterward.

## Limitations And Unverified Claims

- This worker authored the corrections and has not independently audited them.
  Findings are represented as addressed by implementation and tests, not as
  independently closed.
- No second physical operating system, architecture, BLAS, or LAPACK runtime
  was available. Scientific fallback was exercised by a controlled mutation
  of the recorded current fingerprint.
- An ordinary forced-Suggests package check was not completed because the
  optional `genlasso` package is unavailable. The forced-Suggests-disabled
  check is not equivalent to that missing lane.
- A full `make check` with examples, tests, and manual was not run.
- The default smoke suite is not green because of the unrelated coupled-KD
  expectation. The runner now reports that failure truthfully; the expectation
  itself remains unresolved.
- The environment schema records the `RNGversion()` explicitly established by
  the materializer. R does not expose a getter for a caller's historical
  `RNGversion`; reproducibility therefore depends on the new explicit
  establishment step rather than inference from caller state.
- The HTML was rebuilt and its schema statements, revision markers, size, and
  digest were checked. A new interactive browser/layout inspection was not
  performed.
- The regenerated frozen instance fixtures are implementation-owned artifacts,
  not independent scientific oracles. G1--G7 and SSRHE parity remain supported
  by the separate legacy parity fixtures.
- No claim is made about the unrelated dirty files in `trend_filtering` or the
  unrelated `gflow/AGENTS.md` modification.

## Reusable Workflow Capture

Classification: no new external workflow artifact needed.

The normative environment-schema CSV, JSON round-trip tests, deterministic
verification seam, registry regeneration script, and grouped-runner contract
test now encode the reusable safety rules inside the package.

## Current State

The implementation correction, generated registry, response, canonical report
source, and rendered report are committed. The `geosmooth` correction is
pushed. The work is ready for independent re-audit. This handoff makes no claim
about the auditor's eventual verdict.
