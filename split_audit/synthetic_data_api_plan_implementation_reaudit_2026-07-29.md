# Synthetic Data API Consolidation Implementation Re-audit

Re-audit date: 2026-07-29

Auditor role: independent package implementation auditor

Correction handoff:
`split_audit/synthetic_data_api_plan_implementation_reaudit_handoff_2026-07-28.md`

Prior implementation audit:
`split_audit/synthetic_data_api_plan_implementation_audit_2026-07-28.md`

Implementation response:
`split_audit/synthetic_data_api_plan_implementation_audit_response_2026-07-28.md`

Accepted architecture:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`

Revisions inspected:

```text
gflow
  report source: 38a4f59cf917566db807d979c52a94a244ab327e
  rendered report: 8baf5cf52e36c1f2113d4daa602c48b586ea5a3b
  current handoff revision: a6a34bfd5abf8dd5edde121bf7f3c9426891295f

geosmooth
  819996d43bba7b25980d073dce046e64cb9d3bed

trend_filtering
  2062282986bd58d36c82c4154159f7df95059b83
```

## Independence and scope

I did not author the correction commits, response, handoff, registry fixtures,
tests, or revised report. I used the handoff as an evidence index, not as the
definition of the audit.

The re-audit independently retested the original five findings and attempted
to break the corrected validator, comparator, registry, serialization,
cross-environment verification, and consumer identity paths. Measurement,
selection fairness, and statistical inference remain not applicable because
this tranche contains no estimator comparison or inferential claim.

## Verdict

**Revise.**

Four of the five original findings are resolved:

- IMP-1: dataset identity, RNG, dimension, support, alias, and frame validation;
- IMP-2: scientific-comparator exact fields and random-frame invariants;
- CON-1: canonical identity in the SSRHE adapters; and
- ART-2: named atomic-vector canonicalization.

ART-1 is substantially improved but remains a blocker. Recipe and instance
foreign keys now resolve, component payloads reconstruct, fixtures exist, and
exact and fallback verification execute. However:

1. the component registry uses serialized R XDR payloads while the accepted
   architecture still requires JSON parameter records;
2. the accepted checksum schema calls for an environment-fingerprint ID, but
   the implemented checksum rows contain selected fields inline; and
3. the implemented fingerprint omits several fields explicitly required for
   complete exact-environment scope, including BLAS/LAPACK file hashes, an
   effective RNG version, math-runtime identity, and the `geosmooth` version.

The correction report simultaneously describes serialized component payloads
as implemented and retains the incompatible JSON schema as the accepted
contract. The worker's stated goal was to fix the implementation without
revising that architecture. Therefore the remaining divergence has neither
been implemented away nor formally adjudicated as an architecture change.

## Prior finding dispositions

| Prior finding | Re-audit disposition | Evidence |
|---|---|---|
| IMP-1 — validator accepts inconsistent objects | Resolved | All original corruptions are rejected; ordinary and frozen identities, RNG metadata, dimensions, aliases, support, frames, and reconstruction are validated. |
| IMP-2 — comparator omits required fields/invariants | Resolved | Both operands validate first; identity/specification fields and parameters are compared; orientation equivalence requires valid reconstructing frames and non-predictor truth coordinates. |
| ART-1 — registry and environment scope | Partially resolved; blocker remains | Foreign keys, fixtures, payload reconstruction, and fallback work, but schema and complete environment-scope requirements still diverge from the accepted architecture. |
| CON-1 — consumer overwrites canonical ID | Resolved | `dataset.id` remains canonical; the historical label is separate; one-dimensional, flat, quadform, and replay integration pass. |
| ART-2 — named vectors are order-sensitive | Resolved | Named atomic vectors canonicalize by radix-sorted name, with explicit regional-dimension exceptions; reordered polynomial hashes now agree. |

## Findings by Audit Charter layer

### 1. Data-generating process

#### Frozen scientific parity remains intact

I independently replayed the scientific fields directly from both frozen
legacy fixtures:

```text
G1--G7
  cases:                 16
  compared fields:      144
  maximum finite numeric absolute difference:
                        4.440892e-16

SSRHE
  one-dimensional cases: 720
  surface cases:          126
  total cases:            846
  compared fields:       3384
  maximum finite numeric absolute difference:
                         8.881784e-16
```

These are the same machine-precision results obtained in the first
implementation audit. Canonicalization and registry regeneration did not
alter the maintained scientific parity surface.

All nine current `test-synthetic-*.R` files passed serially during this
re-audit, including the full G1--G7 and SSRHE suites after the final compiler-
fingerprint adjustment. This independently closes the handoff's limitation
that the entire focused group had not been rerun after that last adjustment.

### 2. Measurement

Not applicable. No performance or diagnostic measurement is reported.

### 3. Estimation and selection fairness

Not applicable. No estimator or tuning comparison is reported.

### 4. Statistical inference

Not applicable. No interval, test, posterior, or inferential ranking is
reported.

### 5. Artifacts and provenance

#### [Blocker ART-R1] Registry schema and exact-environment scope still diverge from the accepted contract

The accepted architecture states:

- component tables contain a JSON `parameters` column, with matrix parameters
  represented as row-major nested arrays
  (`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1975-1977`);
- `checksums.csv` records an environment-fingerprint ID (`:1978-1980`); and
- a complete fingerprint includes:
  - OS name, release, version, and machine;
  - RNG, normal, and sample kinds plus effective `RNGversion()`;
  - absolute BLAS and LAPACK paths plus file SHA-256 digests;
  - `libm` or equivalent math-runtime identity; and
  - `geosmooth`, registry, evaluator, and dependency versions
    (`:2028-2042`).

The implementation instead writes:

```text
component.id, kind, family, version, component.sha256,
payload.serialization, component.payload.hex
```

where `component.payload.hex` is an R version-3 XDR serialization
(`geosmooth/scripts/update_synthetic_registry.R:17-30`). Recipe compatibility
and metadata are also stored as serialized R payloads (`:33-51`). There is no
JSON parameter column.

The checksum table has no environment-fingerprint ID. Its inline fingerprint
is generated by
`geosmooth/R/synthetic_registry.R:106-165` and omits:

- `Sys.info()["machine"]`;
- effective `RNGversion()`;
- BLAS and LAPACK file digests;
- a concrete math-library/runtime identity;
- the `geosmooth` package version; and
- a separately identified environment-fingerprint record.

The `math.runtime` value records only long-double capability and size. The
BLAS field is a path, while LAPACK is a path plus version; neither has the
required file SHA-256. `dependency.versions` covers `dgraphs`, `digest`,
`MASS`, and `Matrix`, but not `geosmooth`.

`.synthetic.checksum.scope.matches()` compares only the fields returned by
that incomplete function (`geosmooth/R/synthetic_registry.R:167-179`).
Consequently all 10 current checksum rows are classified as exact-scope even
though the accepted scope has not been fully recorded:

```text
rows classified exact-environment: 10 / 10
```

The registry test defines `required.scope` from the same incomplete field list
(`geosmooth/tests/testthat/test-synthetic-registry.R:65-72`), so it cannot
detect the omitted architecture fields.

The revised report is internally contradictory. Its implementation summary
says the registry uses “committed serialized component payloads”
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:395-399`),
while its normative schema still requires JSON parameters at `:1975-1977`.
The response says the architecture was not revised, so the implementation
summary cannot silently replace the accepted schema.

Why this remains blocking:

- exact versus scientific-fallback verification is a release gate, not merely
  descriptive metadata;
- an unrecorded library/runtime/version change is misclassified as the exact
  environment, causing the wrong verification policy to run; and
- adopting R-XDR component payloads instead of the accepted inspectable JSON
  schema is a substantive registry-design decision that requires either
  implementation correction or explicit architecture adjudication.

Required correction:

1. Implement the accepted JSON component schema and environment-fingerprint
   relationship, or explicitly revise the architecture and obtain a fresh
   architecture acceptance for the R-XDR/inline design.
2. Record and compare every promised environment field, including file
   digests, effective RNG version, machine, math-runtime identity, and
   `geosmooth` version.
3. Make the package test's required scope derive from the accepted schema
   rather than mirroring the implementation's current field list.
4. Remove the contradictory schema or implementation statement from the
   canonical report after the design choice is adjudicated.

#### Registry integrity improvements independently confirmed

The other parts of ART-1 now work:

```text
registered recipes:                 106
instance checksum foreign keys:     10 / 10 resolve
committed fixture paths:            10 / 10 exist
current exact content checksums:     10 / 10 replay
former process-wide bypass:         absent from package source
```

Component payloads are deserialized and checked against expected kind, family,
version, and canonical component digest before use. Public default recipe
resolution follows the normalized foreign keys. String parsing remains only
in the private source builder used for registry generation and transitional
G-family overrides.

I also forced a checksum-scope mismatch in process without altering repository
files. `materialize.synthetic.instance("geosmooth.g1.default.v1")` selected
`verification.scope = "scientific-parity"` and returned a valid dataset. Thus
the fallback implementation itself works; the blocker is the incomplete
definition of when it should be selected.

#### Artifact and revision facts reconcile

The supplied artifact facts match independently:

```text
correction handoff
  lines:       378
  bytes:    14,650
  SHA-256: bfb974e54557b01f23bdfd2401cf0c08a246e5cca14ded4534e170e3af381d45

implementation response
  lines:       219
  bytes:    10,801
  SHA-256: 8a465acc1179721cf3381b38a19649bc6fb2888660ff92985ab9db7867385921

committed implementation audit
  lines:       517
  bytes:    21,078
  SHA-256: b062f12b93a9a9267b6f31d6a21a7b18c16d22711efd6ca2685f51d17597c2f1

canonical R Markdown
  lines:      2,413
  bytes:    101,427
  SHA-256: 683b0dee1c8c24b09667a76c473e5d4159058c2abba4a57b1b37554defb3c443

self-contained HTML
  lines:        4,415
  bytes:    1,365,714
  SHA-256: 6e6c0233d7ffb96932f0314b00be641c7cee21c6cf6d1576cb977f2bf541c6e4
```

`geosmooth` remained clean after the re-audit tests. Before this audit report
was written, the current `trend_filtering` and `gflow` dirt was limited to the
unrelated files disclosed by the handoff.

### 6. Estimator and implementation correctness

#### IMP-1 — resolved

The corrected validator now:

- reconstructs the complete ordinary ID exactly;
- binds a frozen ID to one registry instance and its recipe;
- recomputes resolved legacy seeds and named stream states;
- checks scalar and regional dimensions;
- checks truth scope/estimand, masks, regions, and transitional aliases;
- validates support before and after materialization; and
- validates stored frame shape, orthonormality, canonical/supplied semantics,
  and predictor reconstruction.

Independent repetitions of the original attacks returned:

```text
corrupt named RNG stream:
  REJECTED: rng.streams violates the named.stream.v1 state contract.

seed inconsistent with ordinary ID:
  REJECTED: Ordinary dataset.id does not exactly encode its canonical fields.

intrinsic dimension changed from 2 to 99:
  REJECTED: Dataset dimension and region metadata are inconsistent.

stored random frame replaced by zeros:
  REJECTED: Stored frame.matrix has invalid dimensions or is not orthonormal.
```

`synthetic.dataset.checksum()` also refuses these objects because validation
precedes checksum issuance.

#### IMP-2 — resolved

`compare.synthetic.dataset()` validates both inputs, compares contractual
identity/specification/RNG fields exactly, tolerance-compares realized numeric
content, and recursively compares recorded parameters. Random-frame
orientation equivalence is allowed only after both frame/reconstruction
validators pass and neither truth evaluator uses predictor coordinates
(`geosmooth/R/synthetic_dataset.R:517-591`).

The original zero-frame and false-dimension examples now fail validation
before comparison. A separately reconstructed, valid sign-flipped G2 frame
passes and compares scientifically equal through predictor distances.

#### CON-1 — resolved

The three SSRHE constructors preserve `object$dataset.id`, copy it to
`canonical.dataset.id`, and store the historical experiment label in
`legacy.dataset.id`
(`trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R:565-613`).
Replay checks the manifest's canonical ID before returning the dataset
(`:638-656`).

The independent integration run returned:

```text
SSRHE synthetic dataset adapter identity integration: PASS
```

It covered one-dimensional, flat, and quadform constructors plus manifest
round trips, validation, and exact checksums.

#### ART-2 — resolved

`.synthetic.canonicalize()` now propagates field names recursively, sorts
named atomic vectors using UTF-8 radix ordering, and preserves order for
`intrinsic.dim.by.region` and `codimension.by.region`
(`geosmooth/R/synthetic_spec.R:50-84`).

Independent reproduction returned:

```text
c(b0 = 1, b1 = 2, b11 = 3)
c(b11 = 3, b0 = 1, b1 = 2)

specification hashes identical: TRUE
```

The literal serializer fixture also passes.

#### [Nonblocking QA-R1] The cross-environment branch lacks a controlled package regression

The handoff says scientific fallback is tested by a controlled scope mismatch.
The committed test does not create one. It branches on the live result of
`.synthetic.checksum.scope.matches()` (`geosmooth/tests/testthat/test-synthetic-registry.R:109-124`).
On the current machine every row matches, so only the exact-environment branch
runs.

I independently forced a mismatch and the fallback passed, so this is not an
implementation blocker. Add a deterministic test seam or controlled row
mutation so both exact and scientific-fallback branches are exercised on every
supported development environment.

#### [Nonblocking QA-R2] Makefile test lanes can return success with failed expectations

The handoff correctly reports that `make test` exited zero while
`test-coupled-kd-selection-csd4.R` contained a failed expectation. The cause is
that `scripts/run_test_group.R:74-77` calls each `testthat::test_file()` without
`stop_on_failure = TRUE` or a final failure aggregation. This affects
`make test`, `make test-all`, and the focused targets that use the same runner.

The synthetic re-audit used explicit `stop_on_failure = TRUE`, so this runner
defect did not mask the results above. It should nevertheless be fixed because
the repository instructions designate these Make targets as package QA gates.

### 7. Rendering fidelity

The canonical R Markdown rerendered successfully to a temporary directory.
Stored and fresh HTML had identical:

```text
lines:       4,415
bytes:   1,365,714
h1:             20
h2:             20
h3:              9
```

The only textual differences were:

- the two build timestamps; and
- the dynamic `gflow` revision (`38a4f59c` in the stored report versus
  `a6a34bfd` at re-audit render time).

The `geosmooth`, `dgraphs`, and `trend_filtering` revisions were unchanged.
No structural rendering drift or missing report section was found. The
unresolved registry-schema contradiction is a content defect, not a rendering
defect.

## Reproduction and falsification evidence

The re-audit:

- replayed all 862 legacy scientific-parity cases directly from RDS fields;
- ran all nine synthetic test files after the final correction;
- repeated every original validator corruption;
- exercised valid and invalid random-frame comparisons;
- reproduced canonical equality for reordered named coefficients;
- resolved every recipe, instance, checksum, and fixture relationship;
- ran the SSRHE consumer integration test;
- forced the scientific-fallback verification branch;
- compared the implemented checksum scope field by field with the accepted
  architecture rather than with the worker-defined test list; and
- independently rerendered and diffed the HTML.

The corrected identity and parity claims survived. The claim of a complete,
architecture-conforming registry scope did not.

## Validation commands

Principal commands included:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  files <- Sys.glob("tests/testthat/test-synthetic-*.R");
  for (f in files)
    testthat::test_file(f, reporter="summary", stop_on_failure=TRUE)'

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  # Direct replay of every G1--G7 and SSRHE frozen scientific field.'

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  # Original identity/RNG/dimension/frame corruptions and canonicalization
  # probes.'

Rscript development/ssrhe_hessian_energy/test_synthetic_dataset_adapter.R

Rscript -e 'pkgload::load_all(".", quiet=TRUE);
  # Process-local registry-row wrapper forcing a compiler-scope mismatch,
  # followed by materialize.synthetic.instance().'

Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file =
    "synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html",
  output_dir = <temporary-directory>,
  quiet = TRUE
)'

git diff --check
wc -l -c <artifacts>
shasum -a 256 <artifacts>
git status --short
```

I did not rerun `make check-fast`. The default lane is known to stop at the
missing `genlasso` dependency gate, and the forced-Suggests-disabled result is
already recorded. The focused corrections were independently exercised at
source level, including the full synthetic suite and consumer integration.

No package source, registry, fixture, test, consumer source, or generated
report was changed by this re-audit.

## Re-audit gate

Another implementation re-audit is warranted after:

1. ART-R1 is resolved either by implementing the accepted JSON/fingerprint
   schema or by obtaining an explicit architecture decision accepting the
   XDR/inline replacement;
2. every promised environment-scope field is recorded and tested;
3. the canonical report contains one internally consistent registry schema;
4. the registry suite deterministically exercises both exact and fallback
   verification paths; and
5. the focused synthetic and consumer integration suites pass again.

The Makefile runner correction may be carried as a separate package-QA task,
but its status should not be represented as a passing default test gate until
failures produce a nonzero exit status.
