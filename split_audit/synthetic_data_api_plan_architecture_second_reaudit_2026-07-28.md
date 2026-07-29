# Synthetic Data API Consolidation Plan Architecture Second Re-audit

Re-audit date: 2026-07-28  
Auditor role: independent package-architecture auditor  
Response under audit:
`split_audit/synthetic_data_api_plan_architecture_reaudit_response_2026-07-28.md`  
Prior re-audit:
`split_audit/synthetic_data_api_plan_architecture_reaudit_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`  
Audited `gflow` revision: `92a61c086f2fa1fa77223edfb02b74a1be3f1a28`

## Independence and scope

I did not author the response, revised canonical report, generated HTML, or
primary generator implementations. I treated the response as a disposition
ledger, not as an audit agenda. I independently traced the revised contracts
from data generation through identifiers, checksums, implementation
requirements, and rendering.

This remains an architecture-only phase. There are no estimator fits, measured
scores, selection comparisons, or inferential claims. The measurement,
estimation/selection, and inference layers are therefore not applicable rather
than passed.

## Verdict

**Revise before acceptance.**

The revision fully resolves six of the eight prior findings and resolves the
central part of the other two. In particular:

- V1--V3 sampling and response algorithms are now exact;
- finite-design truth is distinguished from population truth;
- declared and observed regions are separated;
- arbitrary materializable closures are removed;
- random-frame scientific comparison and environment scope are defined;
- the SSRHE seed-policy identifiers are corrected;
- geometry/sampling support ownership is explicit; and
- the uniform-box ordering contract is consistent.

Two blockers remain:

1. exact G1--G7 parity is still under-specified for the G2 random frame and G4
   stratified draw sequence; and
2. caller-supplied `recipe.id` values can give scientifically different
   specifications the same derived `dataset.id`.

Two major implementation/compatibility comments also require correction or
explicit human adjudication: `rng.policy` is consumed by the checksum contract
but absent from the required dataset fields, and the plan rejects a currently
valid G7 `zero.parts = integer(0)` call without declaring that compatibility
break.

## Disposition of the eight response claims

| Finding | Second re-audit disposition | Evidence |
|---|---|---|
| DGP-R1: exact random algorithms | Partially resolved; blocker remains | V1--V3, disk, clustered, torus, and Dirichlet paths are fixed. G2 frame construction and G4 stratum/noise traversal remain unspecified. |
| DGP-R2: population versus finite-design truth | Resolved, with a wording comment | The two scopes and their estimands are now explicit, and realized normalizers are recorded. |
| DGP-R3: declared and observed regions | Resolved | G4 requires both strata and fails before drawing; G7 permits absent endpoint regions. |
| ART-R1: evaluator identity and checksums | Resolved | Materializable closures are rejected and the specification-hash payload is defined. |
| ART-R2: frame/environment comparison | Resolved | Exact and scientific-comparison scopes, frame invariants, and environment fingerprint fields are explicit. |
| IMP-R1: SSRHE seed-policy identifiers | Resolved | The examples use `ssrhe.flat.v1` and `ssrhe.quadform.v1`. |
| IMP-R2: geometry/sampling support | Resolved | Static containment/equality and realized-row validation are specified. |
| IMP-R3: uniform-box ordering | Resolved | `"ascending"` is limited to interval and truncated-normal samplers. |

## Findings by Audit Charter layer

### 1. Data-generating process

#### [Blocker DGP-RR1] Exact legacy parity still does not define the G2 and G4 random algorithms

The revised report correctly specifies exact algorithms for the SSRHE V1--V3
paths and for several G-family samplers
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1614-1646`).
It also requires fixtures for every G1--G7 default and nondefault branch and
requires the frame draw order to be recorded (`:1925-1934`).

The G2 frame contract remains distributional rather than algorithmic:

- `frame = "random.orthonormal"` records only a policy and materializer-owned
  stream (`:586-601`);
- no geometry signature carries a `frame.algorithm` or equivalent version
  (`:581-589`, `:1357-1383`); and
- no prose fixes the primary source's exact operation:
  `matrix(rnorm(D * d), nrow = D, ncol = d)` followed by
  `qr.Q(qr(M))[, seq_len(d), drop = FALSE]`
  (`geosmooth/R/dgp_library.R:177-182`).

This matters even though the cross-environment scientific comparator correctly
treats ambient orientation as a nuisance. The same-environment exact checksum
contains the realized frame (`:1879-1888`), and the acceptance criteria require
zero unexplained legacy scientific-field mismatches (`:2137-2140`). A different
QR sign convention or another valid orthonormal-frame algorithm changes that
exact payload.

Independent falsification using the current G2 source produced:

```text
source Q[1,1]                         = -0.7863460467199768
sign-canonical alternative Q[1,1]    =  0.7863460467199768
maximum predictor difference         =  1.2839703276107632
source orthonormal error              =  2.22e-16
alternative orthonormal error         =  2.22e-16
maximum squared-distance difference  =  6.66e-16
```

Both frames satisfy the current geometry and scientific-equivalence contracts,
but only one reproduces the legacy exact payload.

G4 has the same unresolved issue. Its revised component signatures have no
versioned stratified-draw algorithm (`:1438-1453`). The report defines
allocation and support but not the order in which it traverses strata, draws
latent coordinates, and draws A's two ambient-noise columns (`:1565-1581`).
The primary source draws A latent coordinates, both A noise columns, then B's
two latent coordinates (`geosmooth/R/dgp_library.R:523-530`).

A plausible generic traversal that draws all strata's latent coordinates before
ambient noise satisfies the current plan but differs under the same seed:

```text
source G4 X[1,2]             =  0.0265959852584500
latent-first alternative     = -0.0307990008380742
identical predictor matrices = FALSE
```

Add versioned, exact compatibility algorithms for:

- G2's column-major Gaussian matrix, base `qr()`, `qr.Q()`, column selection,
  and any sign/pivot convention; and
- G4's stratum traversal, within-stratum coordinate order, ambient-noise column
  order, row assembly, and response reset.

If the plan intends invariants rather than exact scientific-field parity for
either path, narrow Phase 1 and acceptance criterion 12 explicitly instead.

#### [Major DGP-RR2] The G7 contract silently rejects a currently valid nondefault input

The revised plan requires `zero.parts` to be a nonempty proper subset
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1583-1590`).
The current exported `dgp.g7()` requires only a proper subset. In R,
`all(integer(0) >= 1L)` and `all(integer(0) <= D)` are true, so
`zero.parts = integer(0)` is accepted
(`geosmooth/R/dgp_library.R:687-708`).

Independent materialization confirmed:

```text
empty zero.parts accepted = TRUE
rows labeled "zero"       = 3
exact zero entries        = 0
maximum row-sum error     = 1.11e-16
```

The case is semantically unusual—a `"zero"` region has no structurally zero
part—but it is part of the current public input surface. The report expressly
excludes arbitrary historical truth closures from exact compatibility
(`:1717-1723`) but does not declare this G7 exclusion. That conflicts with
Phase 1's coverage of G1--G7 nondefault branches (`:1925-1932`) and the broad
legacy-parity acceptance criterion (`:2137-2138`).

Choose one policy:

- preserve the empty subset and define its face/region semantics;
- reject it in both the compatibility wrapper and new constructor with a
  lifecycle transition; or
- explicitly exclude it from the maintained parity surface and record the
  breaking change.

#### [Nonblocking DGP-RR3] Introductory wording still implies that every truth is population-scoped

The core scientific contract is corrected, but the DGP callout, ownership card,
current-library description, and five-part abstraction still say “population
truth” or “population-truth specification”
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:342-376`,
`:395-399`, `:493-500`, and `:518-526`).

Because the plan now intentionally includes finite-design truth, use “truth
estimand” or “truth specification” in these general descriptions. This is
terminology drift, not a failure of the detailed estimand contract.

### 2. Measurement

Not applicable. No measured score or diagnostic is reported.

### 3. Estimation and selection fairness

Not applicable. No estimator comparison or tuning selection is reported.

### 4. Statistical inference

Not applicable. No inferential comparison is reported.

### 5. Artifacts and provenance

#### [Blocker ART-RR1] Ordinary dataset identifiers can collide across different specifications

`synthetic.spec()` publicly accepts `recipe.id`, and the worked examples assign
those IDs directly rather than resolving them from the registry
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:915-930` and
`:956-974`). For ordinary materialization, the derived ID is:

```text
<recipe-or-spec-id>__n<n>__seed<seed>__rng<rng-policy>
```

When `recipe.id` is present, it replaces the specification hash
(`:1329-1338`). Therefore two specifications with:

- the same caller-supplied `recipe.id`,
- different geometry, sampling, truth, response, compatibility, or scientific
  metadata,
- the same `n`, seed, and RNG policy

receive the same `dataset.id`. Their content hashes differ, but their primary
dataset identifier collides.

Registry resolution validates foreign keys (`:1813-1816`), but the constructor
contract does not require a non-null `recipe.id` to resolve to a matching
registry record. The current statement that derived IDs prevent callers from
attaching frozen identity to different content is consequently too strong.

Use one of these contracts:

- include a specification-hash suffix in every ordinary dataset ID, even when a
  recipe label is present;
- permit `recipe.id` only on specifications returned by
  `synthetic.registry.spec()` and verify their canonical component payload; or
- separate a caller-facing recipe label from a content-bound recipe identity,
  with ad hoc specifications always using their hash.

This is a blocker because stable identifiers and provenance are central
architecture contracts, and an implementation that follows the current text
cannot reliably use `dataset.id` as an identity key.

### 6. Estimator and implementation correctness

#### [Major IMP-RR1] `rng.policy` is absent from the required dataset object fields

The required `synthetic_dataset` fields list includes `seed` and `rng.streams`
but not `rng.policy`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:664-676`).
Later contracts treat `rng.policy` as a top-level materialized field:

- the dataset checksum payload includes `seed, rng.policy, rng.streams`
  (`:1820-1831`);
- the scientific comparator compares seed policies exactly (`:1871-1873`);
- the derived dataset ID includes the RNG policy (`:1329-1333`); and
- acceptance requires the same specification, size, seed, and RNG policy to
  reproduce canonical content (`:2139-2140`).

Add `rng.policy` to the required object contract and to validation. Storing it
only indirectly in `parameters` or `provenance` would not satisfy the explicitly
top-level checksum payload.

### 7. Rendering fidelity

No new rendering blocker was found.

The canonical R Markdown rerendered successfully. The original and independent
HTML files both had 4,203 lines and 1,353,228 bytes. After replacing only the
two dynamic Eastern build timestamps, both had SHA-256:

```text
7f94f9343975d5e90208caea5d027a46c996cf4f98742575a276d174dc523b43
```

The generated body contains 20 `h1`, 20 `h2`, 10 `h3`, three report tables, and
25 preformatted blocks. The retained 390-pixel ledger screenshot was inspected
and shows all three columns with readable wrapping. The implementer's recorded
completed-DOM measurements report no document-level overflow at 1,440, 768, or
390 pixels. A fresh local-file browser measurement was unavailable in this
audit surface, so that last responsive measurement is corroborated by the
screenshot and deterministic render rather than independently repeated.

## Reproduction and falsification evidence

The principal falsification checks were:

- replacing the legacy G2 frame by an equally orthonormal sign convention and
  comparing exact predictors and invariant distances;
- changing the otherwise valid G4 traversal to draw all latent coordinates
  before ambient noise;
- materializing the current G7 implementation with
  `zero.parts = integer(0)`;
- tracing two different component payloads through the caller-supplied
  `recipe.id` branch of the dataset-ID rule;
- tracing `rng.policy` from the public materializer into the required object,
  checksum, comparator, and acceptance contracts; and
- independently rerendering and timestamp-normalizing the HTML.

There is no numerical performance headline, quality flag, pooled effect,
selection protocol, or dependence-sensitive interval to recompute in this
design-only phase.

## Artifact and provenance validation

Current artifacts:

```text
re-audit response
  lines: 143
  bytes: 6,149
  SHA-256: 5e81a5b0065c920aebecdcfb57cc2d38c1fad8c2316f851d6df10e15c271144f

R Markdown
  lines: 2,226
  bytes: 90,254
  SHA-256: dc7ef89a7d63f9d9263e6d48245957d9480292d6ce74f2e57f248cf7ad7f3493

HTML
  lines: 4,203
  bytes: 1,353,228
  SHA-256: 99cea8ac1c7e6570bcaccc920988a33fed4e51f8234ecf8e17088316dc8814b4

implementer handoff
  lines: 229
  bytes: 8,974
  SHA-256: 147f7131ebd0697fb08d6f36e0fc67feaac95b33479827f3bbdc2ca8427c9468
```

All six repository revisions recorded in the report and handoff still match
their current `HEAD` values:

```text
gflow       92a61c086f2fa1fa77223edfb02b74a1be3f1a28
geosmooth   5ecf721eae2765d8cde01d7fb82b17ff8bde8599
dgraphs     fcd3d23d9fe80efddfd6fa8d404f8a14b980c819
trend_filtering
            e213170a6752d04302eb585d2e7c207efaa538a1
grip        c5c85420f5bbebc82c53dfa9b27cacaae8fd284a
vaginal_community_trajectory_types
            258b0113d625cbb57b3300721f7d583d178338dc
```

The response itself follows the workflow: it maps each prior finding to a
disposition, acknowledges unimplemented work, and does not assert or pressure
an independent verdict.

## Validation commands

Principal commands included:

```sh
Rscript -e 'source(
  "/Users/pgajer/current_projects/geosmooth/R/dgp_library.R"
); ...'

Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file = "second-reaudit-render.html",
  output_dir = <temporary-directory>,
  quiet = TRUE
)'

shasum -a 256 \
  split_audit/synthetic_data_api_plan_architecture_reaudit_response_2026-07-28.md \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html \
  split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md

wc -l -c \
  split_audit/synthetic_data_api_plan_architecture_reaudit_response_2026-07-28.md \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html \
  split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md
```

I also used line-numbered source inspection, targeted searches, normalized HTML
hash comparison, structural HTML inspection, and inspection of the retained
mobile ledger screenshot.

No `make check-fast` or `make check` was run. No package source, export,
documentation, or test was changed, and package checks cannot validate this
unimplemented architecture.

## Required revision and re-audit

Before acceptance:

1. version and specify the exact G2 frame and G4 stratified random algorithms,
   or explicitly narrow the claimed exact-parity surface;
2. make ordinary `dataset.id` content-bound even when a recipe label is
   supplied;
3. add `rng.policy` to the required dataset object and validator;
4. preserve, transition, or explicitly exclude the G7 empty-`zero.parts`
   behavior;
5. replace the residual general “population truth” wording with scope-neutral
   terminology;
6. rerender and update the factual handoff/response; and
7. return the revised artifacts for another independent disposition.
