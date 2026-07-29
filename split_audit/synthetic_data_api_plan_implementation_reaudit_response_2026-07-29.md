# Response to the Synthetic Data API Implementation Re-audit

Date: 2026-07-29

Re-audit:
`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_reaudit_2026-07-29.md`

## Outcome

The remaining blocker, ART-R1, and both nonblocking QA comments were addressed
in:

```text
geosmooth
  a950b27b5ae73f34ea937a718a250f1aac9c8523
```

The accepted architecture was not revised. The implementation was changed to
match its JSON component schema and normalized environment-fingerprint
relationship.

| Finding | Disposition | Principal evidence |
|---|---|---|
| ART-R1 — registry schema and exact-environment scope diverge from the accepted contract | Addressed | Component and recipe values are typed JSON; checksums reference normalized environment rows; every accepted fingerprint field is recorded and compared; the canonical report now states the implemented JSON schema consistently. |
| QA-R1 — cross-environment fallback lacks deterministic package coverage | Addressed | A private verification seam is used by production materialization and tests; the registry suite forces both matching and deliberately mismatched scope records on every run. |
| QA-R2 — grouped Makefile lanes can return success after failed expectations | Addressed | `scripts/run_test_group.R` calls `testthat::test_file(..., stop_on_failure = TRUE)`; `make test` now exits nonzero at the known unrelated coupled-KD failure. |

## ART-R1: accepted JSON component schema

The component ledgers now have:

```text
component.id, kind, family, version, parameters, component.sha256
```

`parameters` is typed JSON rather than serialized R data. The codec preserves:

- integer, double, logical, and character types;
- named and unnamed vectors;
- named and unnamed lists;
- typed missing values and `NULL`;
- full-precision doubles needed to reproduce canonical specification hashes;
  and
- matrices as row-major nested arrays with declared `nrow` and `ncol`.

The public registry resolver reads the JSON parameters, reconstructs the
versioned component class from its kind and family, and verifies the canonical
component SHA-256. Recipe compatibility and metadata are also stored as JSON
in `compatibility.recipe` and `metadata`. Public recipe resolution no longer
depends on an R-XDR component or recipe payload.

XDR remains in one separate and intended role: the canonical dataset content
checksum still uses the accepted `synthetic-content-v1` / `R-xdr-v3` payload
contract.

## ART-R1: normalized environment fingerprints

`checksums.csv` now contains an `environment.fingerprint.id` foreign key rather
than inline environment fields. It resolves exactly once in:

```text
inst/synthetic_registry/environment_fingerprints.csv
```

The normative required-field contract is:

```text
inst/synthetic_registry/environment_fingerprint_schema.csv
```

That schema requires:

- full R build string, platform, architecture, and C++17 compiler identity;
- OS name, release, version, and machine;
- endianness;
- RNG policy, uniform kind, normal kind, sample kind, and the explicitly
  established `RNGversion`;
- absolute BLAS and LAPACK paths and SHA-256 digests;
- LAPACK version;
- linked `libm` or equivalent math-runtime identity;
- `geosmooth`, registry-schema, and evaluator-registry versions; and
- dependency versions.

The generated fingerprint on the registry host contains concrete values for
every required field. Its BLAS and LAPACK paths are absolute and each has a
file digest. On macOS, the math-runtime identity records the absolute `libR`
path and digest plus the linked `libSystem` runtime line.

Each environment row has a full content SHA-256 and an ID derived from that
hash. Runtime exact-scope matching resolves the checksum foreign key, verifies
the fingerprint row's hash and ID, constructs the current fingerprint, and
compares every field in the normative schema. A missing compiler or library on
another host therefore produces a scope mismatch and scientific fallback
rather than an incorrectly classified exact environment.

The materializer now explicitly establishes the recorded current
`RNGversion()` and all three RNG kinds before both legacy and named-stream
draws. Caller RNG kind and state preservation remains in force.

## QA-R1: deterministic exact and fallback branches

Frozen-instance verification was factored into the private
`.synthetic.verify.frozen.instance()` helper used by
`materialize.synthetic.instance()`. It accepts an optional current fingerprint
only as a private test seam; there is no process-wide or public verification
bypass.

The registry test now:

1. verifies the ordinary exact-environment branch with the live complete
   fingerprint;
2. makes a controlled copy whose BLAS SHA-256 is deliberately different;
3. requires scope matching to return false; and
4. runs the same production verifier and requires
   `verification.scope = "scientific-parity"`.

Thus both verification policies execute deterministically on the registry host
and do not depend on whether the live development machine happens to differ
from the recorded environment.

## QA-R2: truthful Makefile exit status

Every grouped test target using `scripts/run_test_group.R` now supplies
`stop_on_failure = TRUE` to `testthat::test_file()`. A contract test protects
that setting.

The correction was demonstrated through the real fast lane:

```text
make test
  stopped at tests/testthat/test-coupled-kd-selection-csd4.R:183
  counted$calls: 76
  expected:       3
  process exit:   nonzero
  make exit:      2
```

The runner no longer converts this failure into a successful package gate. The
unrelated coupled-KD expectation itself was not changed.

## Registry regeneration and validation

The registry was regenerated twice. The complete binary diff digest for all
changed `inst/synthetic_registry` artifacts was identical before and after the
second run:

```text
5eb86f627ea1e369751750b4ce47b2fb084d4cb020424a133d446c08ae03ed6e
```

All nine `test-synthetic-*.R` files passed serially with
`stop_on_failure = TRUE`, including:

- G1--G7 frozen scientific parity;
- the full SSRHE parity suite;
- adversarial validation and comparison;
- typed JSON and row-major matrix round trips;
- all 106 recipe reconstructions;
- all instance, checksum, fixture, and environment foreign keys; and
- deterministic exact and scientific-fallback verification.

The new grouped-runner contract test also passed.

From `/Users/pgajer/current_projects/trend_filtering`:

```sh
Rscript development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
```

returned:

```text
SSRHE synthetic dataset adapter identity integration: PASS
```

`make document` completed successfully.

With forced Suggests disabled, `make check-fast` completed with the same two
warnings and five notes previously associated with non-mainstream
dependencies, existing documentation drift, source/repository artifacts, an
unimported `tail()` call, and an existing Rd formatting issue. No new
registry-schema, JSON dependency, installation, namespace, code, compiled-code,
or vignette failure was reported. The ordinary forced-Suggests lane remains
blocked by the unavailable optional `genlasso` package.

## Canonical report correction

The implementation-state section in
`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd` no longer
describes R-XDR component payloads. It now records typed JSON parameters,
normalized checksum-to-environment foreign keys, the normative fingerprint
schema, and the complete scope fields. The normative registry section and
implementation summary therefore describe the same design.

## Scope boundary

These changes address ART-R1 and QA-R1/R2. They do not implement or retire
later synthetic families, change the accepted package-ownership boundary, or
include basin and extrema construction. The correction is ready for independent
re-audit; this response does not claim auditor acceptance.
