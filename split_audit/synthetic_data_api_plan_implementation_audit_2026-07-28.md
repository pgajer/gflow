# Synthetic Data API Consolidation Implementation Audit

Audit date: 2026-07-28
Auditor role: independent package implementation auditor
Handoff under audit:
`split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md`
Accepted architecture:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`
Prior architecture verdict:
`split_audit/synthetic_data_api_plan_architecture_third_reaudit_2026-07-28.md`

Implemented revisions inspected:

```text
gflow
  b3f33df445d412fa6d07d9a20bfba2530e97e513

geosmooth
  6396af1579f63e48afdd0f26a6e04a9bcdf0d797

trend_filtering
  987822c8f27f4fec75b74e16e4c41822cc8a1d85
```

## Independence and scope

I did not author the implementation, handoff, parity fixtures, registry
ledgers, canonical report, or migrated consumer. I treated the handoff as an
evidence index rather than as the audit agenda.

This audit covers the implemented synthetic-data tranche in `geosmooth`, the
primary SSRHE consumer migration in `trend_filtering`, the committed parity
fixtures and registries, and the current report/handoff artifacts. It does not
authorize retirement of the partial or pending families listed in the family
disposition ledger.

There is no estimator ranking, tuning comparison, interval estimate, or
inferential conclusion in this tranche. The measurement, selection-fairness,
and statistical-inference layers are therefore not applicable. The audit
concentrates on the DGP, implementation, identity, replay, registry, consumer,
and artifact contracts.

## Verdict

**Revise.**

The implemented generators reproduce the frozen G1--G7 and SSRHE scientific
fields, but four blocker-level contract failures remain:

1. `validate.synthetic.dataset()` accepts internally inconsistent identity,
   RNG, dimension, and random-frame metadata;
2. `compare.synthetic.dataset()` can report equality after scientifically
   required metadata and frame invariants have been corrupted;
3. the registry is not the normalized, referentially complete, environment-
   scoped registry specified by the accepted architecture; and
4. the migrated SSRHE adapters overwrite the canonical `dataset.id` and return
   objects that no longer satisfy the `synthetic_dataset` validator.

A major canonicalization defect also gives scientifically identical polynomial
specifications different specification hashes when their named coefficients
are supplied in a different order.

These are not objections to the scientific formulas or to the decision to
centralize reusable generators in `geosmooth`. They are failures in the
identity and replay layer that the architecture made acceptance gates. The
implementation should not be released as the canonical API, and the old
generator layer should not be retired, until the blockers are closed and
independently re-audited.

## Finding summary

| Severity | ID | Finding |
|---|---|---|
| Blocker | IMP-1 | Dataset validation does not enforce the stated identity, RNG-stream, dimension, frame, or post-materialization support contracts. |
| Blocker | IMP-2 | The scientific comparator omits required exact fields and random-frame invariants. |
| Blocker | ART-1 | Registry construction, foreign keys, frozen-fixture links, and environment scope do not match the accepted normalized-registry contract. |
| Blocker | CON-1 | The migrated SSRHE consumer overwrites canonical dataset IDs and returns invalid `synthetic_dataset` objects. |
| Major | ART-2 | Canonicalization does not sort named atomic vectors, so semantically identical polynomial specifications receive different hashes. |

## Findings by Audit Charter layer

### 1. Data-generating process

#### Frozen scientific parity is reproduced

The implemented formula and seed paths reproduce the committed reference
fields. I independently replayed every fixture case without using the
`testthat` expectations as the calculation:

```text
G1--G7
  cases:                 16
  compared fields:      144
  maximum absolute difference over finite numeric fields:
                        4.440892e-16

SSRHE
  one-dimensional cases: 720
  surface cases:          126
  total cases:            846
  compared fields:       3384
  maximum absolute difference over finite numeric fields:
                         8.881784e-16

Frozen instance checksums:
  exact matches:          10 / 10
```

The SSRHE fixture's recorded source SHA-256,
`d15979f0d031baa67c54ada49bf9808f7f90356c990e31dd035ae3278beef4fb`,
matches the parent revision of the migrated
`ssrhe_order3_l1_validation_helpers.R`. The G1--G7 fixture's recorded source
revision, `5ecf721eae2765d8cde01d7fb82b17ff8bde8599`, exists in the `geosmooth`
history.

All eight `test-synthetic-*.R` files also completed with zero failures,
warnings, or skips during this audit. Those results support the implemented
formulas and frozen algorithms. They do not resolve the contract failures
below because the current tests do not attempt the relevant corruptions.

### 2. Measurement

Not applicable. No measured performance or quality score is asserted.

### 3. Estimation and selection fairness

Not applicable. No estimator or tuning procedure is compared.

### 4. Statistical inference

Not applicable. No uncertainty or inferential claim is made.

### 5. Artifacts and provenance

#### [Blocker ART-1] The committed registry is not the accepted normalized registry

The accepted architecture requires:

- component tables with serializable parameters;
- `synthetic.registry.spec()` reconstruction by resolving component foreign
  keys, without dispatching on string fragments;
- a real `instances.checksum.id -> checksums.checksum.id` foreign key;
- fixture paths for frozen instances; and
- a complete environment fingerprint governing exact-checksum enforcement
  (`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1924-1948`
  and `:1991-2005`).

The implementation differs materially:

1. Every row in `instances.csv` has
   `checksum.id = "synthetic-content-v1"`. No such key exists in
   `checksums.csv`; all 10 checksum foreign keys are unresolved.
2. Every `fixture.path` in `instances.csv` is empty.
3. The component tables contain only component ID, kind, family, version, and
   digest. They do not contain the parameters needed to reconstruct a
   component from the ledger.
4. `synthetic.registry.spec()` selects constructors through regular-expression
   and string-fragment dispatch, then verifies only the resulting final hash
   (`geosmooth/R/synthetic_registry.R:344-380`). It does not resolve the
   `geometry.id`, `sampling.id`, `truth.id`, and `response.id` foreign keys in
   `recipes.csv`.
5. `checksums.csv` records R version, platform, architecture, and endianness,
   but omits the architecture's compiler, operating-system, RNG, BLAS/LAPACK,
   math-runtime, registry/evaluator, dependency-version, and scientific-
   tolerance scope.
6. Public verification can be disabled process-wide with
   `options(geosmooth.registry.verify = FALSE)`
   (`geosmooth/R/synthetic_registry.R:30-44`). That option is needed by the
   current regeneration script, but it is not confined to regeneration.

The test named “instance and checksum ledgers are normalized and scoped”
matches rows by `instance.id` and duplicated `content.sha256`; it never tests
the `checksum.id` foreign key
(`geosmooth/tests/testthat/test-synthetic-registry.R:41-55`). The component
foreign-key test checks that recipe IDs appear in a combined component table,
but does not reconstruct components from those rows (`:7-29`).

Independent foreign-key reproduction returned:

```text
unique instances.csv checksum.id:
  synthetic-content-v1

all instance checksum IDs resolve in checksums.csv:
  FALSE
```

Required correction:

- either implement the accepted normalized resolver and schema, including
  complete checksum scope and real foreign keys, or explicitly revise the
  architecture and obtain a new architecture acceptance before treating the
  code-owned constructor registry as equivalent;
- make each instance's checksum ID resolve exactly once;
- populate and validate fixture paths or remove that field through an approved
  schema change;
- test every registry foreign key, version, duplicate, expected kind, fixture
  path, and environment-scope rule; and
- confine any verification bypass to a private regeneration context that
  cannot alter ordinary API behavior.

#### Artifact facts are internally consistent

The handoff's line, byte, and digest claims independently match:

```text
handoff
  lines:       394
  bytes:    15,357
  SHA-256: a4dcaea5a50b3c363ace8a54aba65f255e4ae795f0914b661e53308e40bebba8

canonical R Markdown
  lines:     2,376
  bytes:    99,220
  SHA-256: 576a89d8c896139e49bc4784b7eb5371111f9789589dc698becf81d9002c1f8b

self-contained HTML
  lines:       4,375
  bytes:   1,363,357
  SHA-256: 2da67ae41ad09c2415a8dfe5bb3b8e4f71fd9a6aa9c6e801a34caf696ed461d8

third re-audit response
  lines:          59
  bytes:       2,802
  SHA-256: 3d8505f1262f9ce02335eafc2ae2e8192bf20f36d6cb8f2c2cb6a2845ad9e67e
```

`geosmooth` was clean at the stated revision before and after the audit.
`trend_filtering` contains the unrelated dirty files disclosed by the handoff,
and `gflow` contains only the disclosed pre-existing `AGENTS.md` modification
in addition to this audit report.

### 6. Estimator and implementation correctness

#### [Blocker IMP-1] Dataset validation accepts inconsistent canonical objects

The accepted contract requires `rng.policy` to agree with `dataset.id` or the
frozen instance row and to have a policy-specific `rng.streams` structure
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:779-795`).
Ordinary IDs must also bind the full specification hash, `n`, seed, and RNG
policy (`:1399-1413`). Validation is required to enforce geometry-specific
invariants and support compatibility before and after materialization
(`:791-795` and `:2299-2305`).

The implemented validator:

- verifies that the RNG policy is one of two strings, but does not inspect
  `rng.streams`;
- looks only for the specification-hash substring in an ordinary ID;
- does not compare the ID's `n`, seed, or RNG policy with the object;
- treats the mere presence of `attr(x, "frozen.instance") == TRUE` as adequate
  ID binding;
- does not validate intrinsic dimension, codimension, resolved seeds, stored
  frame dimensions/orthonormality/reconstruction, or policy-specific stream
  states; and
- does not perform post-materialization support or geometry checks
  (`geosmooth/R/synthetic_dataset.R:123-202`).

The support checker is called only before draws
(`geosmooth/R/synthetic_materialize.R:176-177`), despite acceptance criterion
22 requiring checks both before and after materialization.

Independent tampering produced:

```text
rng.streams = "corrupt"                         validates = TRUE
seed = 999 with dataset ID encoding seed 17     validates = TRUE
dataset ID encoding seed 999 with seed = 17     validates = TRUE
rng.policy = "legacy" with named-stream ID/data validates = TRUE
intrinsic.dim = 99 for a two-dimensional G2     validates = TRUE
stored G2 random frame replaced by all zeros    validates = TRUE
```

Because `synthetic.dataset.checksum()` calls this validator and then hashes the
object (`geosmooth/R/synthetic_dataset.R:239-245`), it can issue a canonical-
looking checksum for these inconsistent objects. That defeats the intended
role of validation as the gate before identity and replay.

Required correction:

- parse and exactly reconstruct every ordinary ID from the object's recipe
  label, specification digest, `n`, seed, and RNG policy;
- bind frozen objects to one complete instance row rather than to a Boolean
  attribute;
- validate the exact stream names, count, state types, and policy-specific
  semantics for both RNG policies;
- validate all scalar dimensions, codimensions, region dimensions, resolved
  seeds, aliases, stored frames, and geometry reconstruction rules;
- run support and geometry invariants on realized data; and
- add adversarial tests for each corruption above, including checksum refusal.

#### [Blocker IMP-2] The scientific comparator omits contractual fields and invariants

The accepted comparator requires exact comparison of IDs, integer/logical
fields, masks, regions, component versions, and seed policies. For random
frames, it must separately validate each frame's dimensions, orthonormality,
and reconstruction of its own predictors
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:2007-2030`).

The implementation compares exactly only:

```text
n, ambient.dim, truth.scope, region, declared.regions,
observed.regions, rng.policy
```

It then compares latent coordinates, truth, response, and either predictors or
predictor distances (`geosmooth/R/synthetic_dataset.R:253-292`). It omits the
ID, specification digest, recipe/tag, seed, RNG streams, latent mask,
intrinsic/codimension metadata, component versions, estimand, resolved design
statistics, offsets, and frame invariants.

Independent falsification returned:

```text
G2 stored random frame replaced by all zeros:
  validates                 TRUE
  comparator equal          TRUE
  reported mismatches       none

G2 intrinsic.dim changed from 2 to 99:
  validates                 TRUE
  comparator equal          TRUE
```

Required correction:

- implement the complete exact-field and tolerance-field contract;
- validate each input's random frame and predictor reconstruction before using
  orientation-invariant distances;
- enforce the truth evaluator's declared ambient-orthogonal invariance before
  allowing random-frame equivalence; and
- add negative comparator fixtures for every exact field and invariant.

#### [Blocker CON-1] The migrated consumer invalidates canonical dataset identity

The accepted architecture states that an ordinary `dataset.id` is derived and
cannot be overridden; only `materialize.synthetic.instance()` may return a
registry-controlled short ID
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1399-1413`).

The migrated SSRHE helpers save the canonical ID in a new field and then
overwrite `object$dataset.id` with a legacy display label:

```text
make.1d.dataset()       lines 570-575
make.flat.dataset()     lines 593-598
make.quadform.dataset() lines 606-612
materialize.dataset()   lines 644-648
```

in
`trend_filtering/development/ssrhe_hessian_energy/ssrhe_order3_l1_validation_helpers.R`.
The object retains class `synthetic_dataset`, but its new ID no longer contains
the specification digest, sample size, seed, or RNG policy. A subsequent
`validate.synthetic.dataset()` or `synthetic.dataset.checksum()` must reject
it. The added `canonical.dataset.id` field does not repair the required
`dataset.id`.

Reproducing the adapter's ID reassignment on `S16.V3` returned:

```text
ERROR: dataset.id and specification.sha256 are inconsistent.
```

Required correction:

- preserve the canonical `dataset.id`;
- store the historical consumer label in a separate field such as
  `legacy.dataset.id` or in the consumer's result row;
- update `dataset.row()` and round-trip logic to keep both concepts explicit;
  and
- add an integration test requiring every adapter-produced object to pass
  `validate.synthetic.dataset()` and `synthetic.dataset.checksum()`.

#### [Major ART-2] Named atomic vectors are not canonically ordered

The checksum contract says named vectors are sorted by name using
locale-independent radix ordering unless their order is scientifically
contractual
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1965-1977`).

`.synthetic.canonicalize()` sorts named lists but returns named atomic vectors
unchanged (`geosmooth/R/synthetic_spec.R:50-72`). Polynomial evaluation looks
coefficients up by name, so coefficient order is not scientifically
contractual.

Two specifications differing only in coefficient order produced identical
scientific datasets but different specification hashes:

```text
coefficients 1: c(b0 = 1, b1 = 2, b11 = 3)
coefficients 2: c(b11 = 3, b0 = 1, b1 = 2)

hash 1:
  96cf4dd0b547933257157dd2cce0a7b057d3171782b00ccb168dc982f06c49e0

hash 2:
  9f724dff7e60421642bda033050dc2410834f24d37e28f0c1850edcd4423f425

hashes identical:                    FALSE
compare.synthetic.dataset() equal:   TRUE
```

Required correction:

- canonicalize named atomic vectors where name, rather than position, defines
  meaning; and
- make order-contractual vectors explicit types or field-specific exceptions
  instead of leaving the distinction implicit.

Because changing canonicalization can change registered specification and
content digests, regenerate and deliberately re-freeze every affected registry
row and frozen checksum after adding a literal serializer regression fixture.

#### The earlier seed-overflow comment is resolved

The implementation computes legacy offsets in double precision and validates
all path-specific derived seeds before RNG state changes
(`geosmooth/R/synthetic_materialize.R:187-194`). Independent maximum-integer
tests pass for draw-free branches and fail before changing caller RNG state
when a frame or response offset is required. This closes the nonblocking seed
comment from the third architecture re-audit.

#### The `simulate()` S3 naming collision is resolved

The exported operation is `materialize.synthetic()`, not
`simulate.synthetic()`. No `simulate.synthetic` definition or export exists.
This avoids accidental S3-method interpretation while retaining the
repository's dot-delimited naming convention.

### 7. Rendering fidelity

No new rendering defect was found in the supplied report artifacts. The
self-contained HTML exists, its digest matches the handoff, and the canonical
R Markdown was not modified after the reported render. The handoff explicitly
discloses that final interactive browser QA was not repeated. That limitation
does not affect this implementation verdict because the blocker findings are
source- and data-contract failures, not presentation failures.

## Reproduction and falsification evidence

The independent audit attempted to disprove the handoff's readiness by:

- replaying all 16 G1--G7 and all 846 SSRHE cases directly from the frozen RDS
  scientific fields;
- replaying all 10 frozen instance content checksums;
- verifying the fixture source revision and pre-migration source digest;
- running all eight synthetic-focused test files;
- corrupting RNG stream metadata, seed/ID agreement, policy/ID agreement,
  intrinsic dimension, and stored random frames;
- comparing corrupted datasets through the public comparator;
- resolving every instance checksum foreign key against `checksums.csv`;
- tracing registry construction from the regeneration script through the
  exported resolver;
- tracing canonical IDs through the migrated consumer; and
- constructing two semantically identical polynomial specifications with
  different named-vector order.

The parity headline survived; the canonical identity and registry headline did
not.

## Validation commands

Principal commands included:

```sh
Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  files <- Sys.glob("tests/testthat/test-synthetic-*.R");
  for (f in files)
    testthat::test_file(f, stop_on_failure = TRUE)'

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  # Direct replay of every G1--G7 and SSRHE fixture field and all frozen
  # instance checksums.'

Rscript -e 'pkgload::load_all(".", quiet = TRUE);
  # Adversarial dataset validation, comparator, and named-vector hash probes.'

Rscript -e '
  instances <- read.csv("inst/synthetic_registry/instances.csv");
  checksums <- read.csv("inst/synthetic_registry/checksums.csv");
  all(instances$checksum.id %in% checksums$checksum.id)'

git show \
  987822c^:development/ssrhe_hessian_energy/\
ssrhe_order3_l1_validation_helpers.R |
  shasum -a 256

wc -l -c <artifacts>
shasum -a 256 <artifacts>
git status --short
git diff --check
```

No package source, registry, fixture, generated documentation, or consumer
source was changed by this audit. I did not rerun `make check-fast`: the handoff
already records its dependency-gate failure and forced-Suggests result, while
the blocker findings reproduce without a package build.

## Re-audit gate

A focused implementation re-audit should be requested only after:

1. closing IMP-1 with adversarial validator and checksum-refusal tests;
2. closing IMP-2 with complete exact-field and random-frame comparator tests;
3. resolving the registry architecture/schema mismatch and every checksum
   foreign key under ART-1;
4. preserving canonical IDs through the SSRHE consumer under CON-1;
5. fixing named-vector canonicalization and deliberately regenerating affected
   hashes under ART-2;
6. rerunning all frozen parity cases and instance checksums; and
7. rerunning the applicable `geosmooth` package QA lane and consumer
   integration test.

Until that gate passes, the correct status is implemented scientific parity
with unresolved canonical-contract blockers, not an accepted canonical API.
