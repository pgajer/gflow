# Response to the Synthetic Data API Implementation Audit

Date: 2026-07-28
Audit:
`/Users/pgajer/current_projects/gflow/split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md`

## Outcome

All five implementation findings were addressed without revising the accepted
architecture. The corrections are implemented in:

```text
geosmooth
  819996d43bba7b25980d073dce046e64cb9d3bed

trend_filtering
  2062282986bd58d36c82c4154159f7df95059b83
```

The changes deliberately preserve the audit's distinction between scientific
parity and canonical-contract correctness. The pre-existing G1--G7 and SSRHE
parity results remain intact; the identity, validation, comparator, registry,
serialization, and consumer layers now enforce the contracts on which that
parity depends.

## Finding dispositions

| Finding | Disposition | Principal evidence |
|---|---|---|
| IMP-1 — validator accepts inconsistent canonical objects | Addressed | Exact ordinary-ID reconstruction; exact frozen-instance binding; policy-specific RNG replay; dimension, alias, mask, frame, reconstruction, and realized-support checks; adversarial checksum-refusal tests. |
| IMP-2 — comparator omits exact fields and frame invariants | Addressed | Both inputs are validated first; contractual identity and metadata are exact fields; random-frame equivalence is restricted to valid frames and compatible truth scope; negative comparison tests added. |
| ART-1 — registry is not normalized, referentially complete, or environment-scoped | Addressed | Serialized XDR component payloads; component and instance foreign-key resolution; real checksum IDs; committed frozen fixtures; complete environment fingerprints; exact-scope checksum and cross-scope scientific verification; public bypass removed. |
| CON-1 — SSRHE adapters overwrite canonical IDs | Addressed | Canonical `dataset.id` is retained, legacy display identity moved to `legacy.dataset.id`, result-row behavior made explicit, and adapter round-trip validation/checksum integration coverage added. |
| ART-2 — named atomic vectors are not canonicalized | Addressed | Named atomic vectors sort by locale-independent radix order except for two explicit order-contractual regional fields; a literal serializer fixture and reordered-polynomial hash regression were added; registries and frozen checksums were regenerated. |

## IMP-1: canonical dataset validation

`geosmooth/R/synthetic_dataset.R` now derives the expected ordinary
`dataset.id` from the dataset's recipe label, full specification SHA-256,
sample size, seed, and RNG policy and requires exact equality. A frozen ID is
accepted only when the object bears the private frozen-instance attribute and
its ID resolves to exactly one complete registry instance. The earlier
Boolean-only bypass is rejected.

The validator also now:

- recomputes the seed plan and requires exact `resolved.seeds`;
- checks the exact named-stream names, seven-integer stream states, and states
  recomputed from each resolved seed without changing caller RNG state;
- checks the legacy stream layout and values;
- validates scalar intrinsic, ambient, and codimension metadata and the named
  regional dimension/codimension vectors;
- checks declared and observed region identities, truth scope, estimand,
  latent masks, and transitional aliases;
- validates stored frame dimensions, finiteness, orthonormality, canonical or
  supplied-frame semantics, and reconstruction of predictors from latent
  coordinates; and
- checks realized box, disk, interval, grid, truncated, rectangular, gapped,
  stratified, and structural-zero support contracts.

`tests/testthat/test-synthetic-contract-adversarial.R` corrupts each relevant
field and requires both validation and canonical checksum issuance to fail.

## IMP-2: scientific comparison

`compare.synthetic.dataset()` now validates each operand before making any
scientific comparison. It compares the following contract categories exactly:

- dataset and specification identity;
- recipe and compatibility identity;
- sample size, dimensions, masks, truth scope, estimand, regions, seed, and
  RNG policy;
- RNG streams and all geometry, sampling, truth, and response specifications;
  and
- component versions and resolved metadata embedded in those specifications.

Latent coordinates, truth, response, and other numeric realized fields remain
tolerance comparisons. Predictor-distance equivalence is used only when both
objects have already passed their frame and reconstruction invariants and the
truth is not evaluated in predictor coordinates. Otherwise predictors are
compared directly. The adversarial suite covers exact-field corruption, an
invalid zero frame, a valid orientation-equivalent random frame, and the
corresponding checksum behavior.

## ART-1: normalized registry and verification scope

The component ledgers now contain a serialized XDR payload for each component,
in addition to its ID, expected kind, family, version, and content hash.
`synthetic.registry.spec()` reconstructs a recipe by resolving its geometry,
sampling, truth, and response foreign keys and deserializing those payloads.
It verifies component kind, family, version, and hash before constructing the
final specification and checking its expected digest. String-fragment
constructor dispatch is no longer the public resolution path; a private
source-builder remains only for deliberate registry regeneration and
parameterized transitional G-family overrides.

Every registered instance now:

- resolves to exactly one recipe and exactly one checksum row;
- names a real checksum ID rather than the old placeholder;
- points to a nonempty committed RDS fixture; and
- materializes with a canonical frozen ID and passes the complete dataset
  validator.

The checksum ledger records R version, platform, architecture, compiler
command and version, operating system, endian mode, RNG kinds, BLAS,
LAPACK/version, math runtime, registry/evaluator version, dependency versions,
and absolute and relative scientific tolerances. In the exact recorded scope,
materialization requires the exact content checksum. Outside that scope, it
compares the regenerated object scientifically with the committed fixture
under the ledger tolerances and records which verification path was used.

The process-wide `options(geosmooth.registry.verify = FALSE)` behavior was
removed. Registry regeneration uses private construction functions and cannot
disable verification for ordinary API calls.

`tests/testthat/test-synthetic-registry.R` now checks exact-one component,
recipe, instance, and checksum relationships; global ID uniqueness; expected
component kind and version; serialized-payload and hash agreement; every
recipe reconstruction; every nonempty fixture; environment-scope completeness;
exact-scope or scientific-fallback verification; and the absence and
ineffectiveness of the former public bypass.

## CON-1: SSRHE consumer identity

The three SSRHE adapters in
`trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R`
now leave `object$dataset.id` unchanged. They copy that canonical identity to
`canonical.dataset.id` for explicit reporting and place the old display label
in `legacy.dataset.id`. Dataset, score, and prediction rows continue to use the
legacy label where the existing experimental tables expect it, while their
separate canonical-ID column carries the replay identity.

`development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R` exercises
one one-dimensional, one flat, and one quadform adapter plus the consumer
round-trip. Every returned and replayed object must pass
`validate.synthetic.dataset()` and `synthetic.dataset.checksum()`, and the test
checks that canonical and legacy identities remain distinct.

## ART-2: canonical named-vector ordering

`.synthetic.canonicalize()` now sorts named atomic vectors by UTF-8 name with
radix ordering. The field name is propagated through recursive
canonicalization so that the two arrays whose order is explicitly contractual,
`intrinsic.dim.by.region` and `codimension.by.region`, retain positional
semantics. Matrix and array dimension order also remains positional and
dimnames remain noncanonical.

The committed serializer fixture records the literal canonical hash for a
reordered named coefficient vector. The adversarial test also constructs two
polynomial specifications with the same named coefficients in different input
orders and requires the same specification hash. All affected registry rows,
fixtures, and checksum relationships were then regenerated deliberately.

## Validation evidence

Commands were run from `/Users/pgajer/current_projects/geosmooth` unless
otherwise noted.

```sh
make document

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  files <- c(
    "tests/testthat/test-synthetic-contract-adversarial.R",
    "tests/testthat/test-synthetic-registry.R"
  );
  for (f in files)
    testthat::test_file(f, reporter="summary", stop_on_failure=TRUE)'

Rscript scripts/update_synthetic_registry.R
# repeated and compared byte-for-byte

make test

make check-fast

_R_CHECK_FORCE_SUGGESTS_=false make check-fast
```

The final serial audit-specific run completed with:

```text
synthetic-contract-adversarial: 34 expectations, 0 failures
synthetic-registry:           1219 expectations, 0 failures
registry regeneration:       byte-identical on consecutive runs
```

The complete focused synthetic test group, including the frozen G1--G7 and
SSRHE parity suites, also passed before the final serial confirmation.

From `/Users/pgajer/current_projects/trend_filtering`:

```sh
Rscript development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R
```

returned:

```text
SSRHE synthetic dataset adapter identity integration: PASS
```

The ordinary `make check-fast` dependency gate remains blocked because the
optional Suggested package `genlasso` is not installed. With forced Suggests
disabled, the check completed with two warnings and five notes already
associated with local/non-mainstream dependencies, existing documentation
drift, and repository artifacts; it did not identify a new failure in this
audit correction. The broader `make test` command exited successfully because
of the current runner behavior but printed one unrelated failing expectation
in `test-coupled-kd-selection-csd4.R` (`counted$calls` was 76 rather than 3).
That baseline test issue is not represented as passing or fixed here.

## Scope boundary

These corrections close the implementation-contract findings only. They do not
retire the partial or pending families in the disposition ledger, expand the
accepted ownership boundary, or claim completion of the later synthetic-data
migration phases. Basin and extrema construction remain outside this work.
