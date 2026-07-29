# Synthetic Data API Implementation Correction Handoff

Status: ready for independent implementation re-audit
Role: implementation worker responding to an independent audit
Date: 2026-07-28

Repositories and revisions:

```text
gflow
  worktree: /Users/pgajer/current_projects/gflow
  branch: main
  audit-response/report-source commit:
    38a4f59cf917566db807d979c52a94a244ab327e
  rendered-report commit:
    8baf5cf52e36c1f2113d4daa602c48b586ea5a3b
  final status before this handoff was added:
    unrelated pre-existing modification to AGENTS.md

geosmooth
  worktree: /Users/pgajer/current_projects/geosmooth
  branch: main
  correction commit:
    819996d43bba7b25980d073dce046e64cb9d3bed
  remote state:
    pushed to origin/main
  final status:
    clean

trend_filtering
  worktree: /Users/pgajer/current_projects/trend_filtering
  branch: main
  consumer correction commit:
    2062282986bd58d36c82c4154159f7df95059b83
  remote state:
    pushed to origin/main
  final status:
    scoped adapter changes are committed; unrelated pre-existing README,
    PDF, LaTeX, bibliography, build-script, and citation-report changes remain
```

The commit containing this handoff is the delivery commit for the handoff
itself. The exact code, report-source, and rendered-report revisions are fixed
above so the evidence does not depend on a self-referential commit hash.

## Goal

Address the five findings in:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md`

without weakening or revising the accepted architecture in:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`

The findings were IMP-1, IMP-2, ART-1, CON-1, and ART-2. Basin and extrema
construction and the later partial/pending family migrations remained outside
scope.

## Work Completed

### Dataset contract

The `synthetic_dataset` validator now reconstructs ordinary dataset identity
exactly, binds frozen identity to one registry instance, verifies each RNG
policy and its recomputed streams, checks dimension and region metadata,
checks aliases and masks, and validates stored frames, predictor
reconstruction, and realized support.

The scientific comparator now validates both operands before comparison,
compares contractual identity and metadata fields exactly, and restricts
orientation-invariant predictor-distance comparison to valid random frames
and compatible truth scope.

### Registry and canonical serialization

The registry now stores serialized XDR component payloads and reconstructs
specifications through exact component foreign keys. Instance checksum IDs
resolve exactly once, frozen fixtures are committed and nonempty, and checksum
enforcement uses a full recorded environment fingerprint. Materialization
requires an exact checksum in the recorded environment and otherwise performs
a scientific comparison with the committed fixture under recorded
tolerances.

The process-wide registry verification option was removed. The regeneration
script uses private source builders instead.

Named atomic vectors now use locale-independent radix ordering unless their
field is one of the two explicit position-contractual regional dimension
vectors. A literal serializer fixture and a reordered-polynomial regression
test were added, followed by deliberate regeneration of registry artifacts and
frozen checksums.

### Consumer identity

The SSRHE adapters retain the canonical `dataset.id`. The historical
experiment label is stored as `legacy.dataset.id`, while
`canonical.dataset.id` remains available explicitly in consumer result rows.
An integration test covers one-dimensional, flat, and quadform adapters and
the consumer round-trip.

### Audit and report artifacts

The finding-by-finding response is:

`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_audit_response_2026-07-28.md`

The canonical report was updated to describe the correction state without
claiming that later migration phases are complete. Its derived self-contained
HTML was rebuilt from the committed R Markdown source.

## Files Changed Or Created

### geosmooth source and tests

```text
/Users/pgajer/current_projects/geosmooth/.gitignore
/Users/pgajer/current_projects/geosmooth/R/synthetic_dataset.R
/Users/pgajer/current_projects/geosmooth/R/synthetic_registry.R
/Users/pgajer/current_projects/geosmooth/R/synthetic_spec.R
/Users/pgajer/current_projects/geosmooth/scripts/update_synthetic_registry.R
/Users/pgajer/current_projects/geosmooth/tests/fixtures/canonical_serializer_v1.csv
/Users/pgajer/current_projects/geosmooth/tests/testthat/test-synthetic-api.R
/Users/pgajer/current_projects/geosmooth/tests/testthat/test-synthetic-contract-adversarial.R
/Users/pgajer/current_projects/geosmooth/tests/testthat/test-synthetic-registry.R
```

### geosmooth generated registry artifacts

```text
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/checksums.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/geometries.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/instances.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/recipes.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/responses.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/samplings.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/truths.csv
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/registry_manifest.md
/Users/pgajer/current_projects/geosmooth/inst/synthetic_registry/fixtures/
```

The fixture directory contains ten committed RDS files, one for each frozen
G-family instance.

### trend_filtering source and test

```text
/Users/pgajer/current_projects/trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R
/Users/pgajer/current_projects/trend_filtering/development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
```

### gflow audit and report artifacts

```text
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_audit_response_2026-07-28.md
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html
/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_handoff_2026-07-28.md
```

## Generated Artifacts

The following sizes and SHA-256 digests were measured after the final report
render:

```text
canonical R Markdown
  lines:    2,413
  bytes:  101,427
  SHA-256:
    683b0dee1c8c24b09667a76c473e5d4159058c2abba4a57b1b37554defb3c443

self-contained HTML
  lines:      4,415
  bytes:  1,365,714
  SHA-256:
    6e6c0233d7ffb96932f0314b00be641c7cee21c6cf6d1576cb977f2bf541c6e4

implementation audit
  lines:      517
  bytes:   21,078
  SHA-256:
    b062f12b93a9a9267b6f31d6a21a7b18c16d22711efd6ca2685f51d17597c2f1

implementation audit response
  lines:      219
  bytes:   10,801
  SHA-256:
    8a465acc1179721cf3381b38a19649bc6fb2888660ff92985ab9db7867385921
```

The audit digest differs from the auditor's original copy only because four
Markdown hard-break trailing spaces were removed before committing it; its
substantive text was not changed.

## Commands Run

From `/Users/pgajer/current_projects/geosmooth`:

```sh
make document

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  files <- Sys.glob("tests/testthat/test-synthetic-*.R");
  for (f in files)
    testthat::test_file(f, reporter="summary", stop_on_failure=TRUE)'

Rscript scripts/update_synthetic_registry.R

Rscript scripts/update_synthetic_registry.R
# generated outputs compared byte-for-byte with the preceding run

make test

make check-fast

_R_CHECK_FORCE_SUGGESTS_=false make check-fast

make clean

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  files <- c(
    "tests/testthat/test-synthetic-contract-adversarial.R",
    "tests/testthat/test-synthetic-registry.R"
  );
  for (f in files)
    testthat::test_file(f, reporter="summary", stop_on_failure=TRUE)'

make clean
git diff --check
```

From `/Users/pgajer/current_projects/trend_filtering`:

```sh
Rscript development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
git diff --check -- \
  development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R \
  development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
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
  split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md \
  split_audit/synthetic_data_api_plan_implementation_audit_response_2026-07-28.md

shasum -a 256 \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html \
  split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md \
  split_audit/synthetic_data_api_plan_implementation_audit_response_2026-07-28.md
```

## Validation

The final serial audit-specific run reported:

```text
synthetic-contract-adversarial
  expectations: 34
  failures:      0

synthetic-registry
  expectations: 1,219
  failures:      0
```

The complete focused synthetic group, including the frozen G1--G7 and SSRHE
scientific-parity files, passed before the final serial run. Consecutive
registry regenerations were byte-identical.

The consumer integration test reported:

```text
SSRHE synthetic dataset adapter identity integration: PASS
```

`make document` completed successfully.

The ordinary `make check-fast` stopped at its dependency gate because the
optional Suggested package `genlasso` was unavailable. With forced Suggests
disabled, the package check completed with two warnings and five notes. The
reported items concerned local/non-mainstream package dependencies, existing
documentation drift, object/repository artifacts, an unimported call, and an
existing Rd formatting issue; no new audit-correction failure was reported.

The broader `make test` smoke lane returned exit status zero because of the
current runner behavior, but its output contained one unrelated failing
expectation:

```text
tests/testthat/test-coupled-kd-selection-csd4.R:183
  counted$calls: 76
  expected:       3
```

That result is recorded as a broader baseline failure, not as a passing test.

One attempted final focused rerun was started concurrently with the
`trend_filtering` integration test. Both processes tried to compile/load the
same local `geosmooth` tree, causing a transient missing-object make error.
After `make clean`, the two audit-specific files were rerun serially and passed
with the counts above. The consumer integration was also rerun serially and
passed.

## Canonical/Generated File Notes

- `geosmooth/R/*.R`, the registry update script, and test sources are
  canonical. The registry CSV files and frozen RDS fixtures are generated by
  `scripts/update_synthetic_registry.R` and are committed package artifacts.
- Roxygen documentation was regenerated with `make document`. No private
  registry or canonicalization helper was exported.
- The R Markdown report is canonical. The sibling HTML is generated by the
  render command recorded above and was rendered after the source commit.
- The implementation audit is auditor-authored evidence. Only four trailing
  Markdown hard-break spaces were removed when it was committed.
- No Markdown note under `/Users/pgajer/.codex/notes` was modified, so the
  shared notes HTML builder was not applicable.

No `geosmooth` source was modified after the final serial adversarial and
registry validation. No scoped `trend_filtering` source was modified after the
passing integration test. The R Markdown source was committed before the final
HTML render and was not modified afterward.

## Limitations And Unverified Claims

- This worker is the author of the corrections and has not independently
  audited them. The audit findings are represented as addressed by the
  implementation and tests, not as independently closed.
- The complete focused synthetic group passed before the final small
  compiler-fingerprint fallback adjustment. After that adjustment, the
  adversarial and registry files were rerun serially; the entire focused group
  was not rerun a second time.
- An ordinary forced-Suggests package check was not completed because
  `genlasso` was unavailable. The forced-Suggests-disabled check is not
  equivalent to that missing lane.
- The broader package suite is not fully green because of the unrelated
  coupled-KD call-count expectation described above.
- No cross-operating-system or cross-BLAS materialization was available during
  this work. The scientific fixture fallback is tested by controlled
  environment-scope mismatch, not by a second physical platform.
- The HTML was rebuilt and its content markers and digest were checked. A new
  interactive browser/layout inspection was not performed.
- Registry fixtures were regenerated by the implementation worker and are not
  an independent oracle. The separate frozen G1--G7 and SSRHE parity fixtures
  remain the scientific parity evidence.
- Later synthetic families marked partial, pending, or undecided in the family
  disposition ledger were not implemented or retired.
- No claim is made that the unrelated dirty files in `trend_filtering` or
  `gflow/AGENTS.md` are correct, complete, or part of this work.

## Reusable Workflow Capture

Classification: no new reusable artifact needed.

The registry regeneration script, literal serializer fixture, adversarial
contract tests, and package QA targets already capture the repeatable parts of
this correction. The existing handoff and package-QA workflows cover the
remaining process requirements.

## Current State

The correction commits, response, canonical report source, and rendered report
are committed. The `geosmooth` and scoped `trend_filtering` commits are pushed.
The implementation correction is ready for independent re-audit. No claim of
audit acceptance is made in this handoff.
