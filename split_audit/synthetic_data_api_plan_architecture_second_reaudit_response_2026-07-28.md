# Synthetic Data API Plan Architecture Second Re-audit Response

Response date: 2026-07-28  
Role: architecture-plan author  
Second re-audit:
`split_audit/synthetic_data_api_plan_architecture_second_reaudit_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`

## Disposition

All five residual findings were accepted. None was deferred or rejected. The
canonical plan and generated HTML were revised. This response records the
author's changes and does not assert a new independent verdict.

## Finding resolutions

### DGP-RR1 — Exact G2 and G4 random algorithms

Resolved by specifying both legacy paths at operation level.

G2 now records `frame.algorithm = "legacy.base.qr.gaussian.v1"`. The algorithm
resolves the explicit or default frame seed, constructs a column-major Gaussian
matrix, calls base `qr(..., LAPACK = FALSE)`, calls
`qr.Q(..., complete = FALSE)`, selects the first required columns, and performs
no sign, pivot, or other post-processing. The frame seed, algorithm ID, and
realized matrix are recorded. Exact parity covers default and explicit
frame-seed paths.

G4 now records `algorithm = "legacy.g4.stratum.sequential.v1"`. It draws A's
latent coordinate, A's second and third ambient-noise columns, B's first and
second latent coordinates, and then assembles A rows before B rows. It evaluates
truth without a draw and applies the exact conditional response reset: positive
response standard deviation resets to `seed + 900000L`; zero standard deviation
returns truth without a reset or draw.

Phase 1, Phase 2, Phase 4, and the acceptance criteria now require exact fixture
comparison for these algorithms in addition to scientific invariants.

### DGP-RR2 — Empty G7 `zero.parts`

Resolved as an explicit breaking exclusion.

The current acceptance of `zero.parts = integer(0)` results from vacuous
`all()` checks and labels rows `"zero"` without creating structural zeros.
That behavior conflicts with the constructor's scientific meaning and is not
promoted into the canonical API.

The new constructor and breaking-release wrapper reject the call with a
dedicated condition and remediation message. It is explicitly excluded from
the maintained parity surface. A pre-release consumer audit must locate every
such call and replace it with either ordinary Dirichlet sampling or actual
structural-zero indices. No temporary canonical recipe preserves the
misleading region labels.

### DGP-RR3 — Scope-neutral introductory terminology

Resolved by changing the general descriptions from “population truth” to
“truth estimand” or “truth specification.” The precise sections still use
“population truth” where they specifically define the population-scoped
estimand and contrast it with finite-design truth.

### ART-RR1 — Content-bound ordinary dataset IDs

Resolved by including the full canonical specification SHA-256 in every
ordinary dataset ID:

```text
<recipe-label-or-adhoc>__spec<full-specification-SHA256>__n<n>__seed<seed>__rng<rng-policy>
```

The readable `recipe.id` no longer substitutes for content identity. Two
different specifications with the same label therefore receive different IDs.
The full digest is also a required top-level `specification.sha256` field,
belongs to the checksum payload, and is recomputed by validation. Registry
recipes carry the expected digest and registry resolution verifies it.

Frozen instance IDs remain registry-controlled exceptions because
`materialize.synthetic.instance()` binds and verifies both specification and
content checksums.

### IMP-RR1 — Required `rng.policy`

Resolved by adding `rng.policy` to the required `synthetic_dataset` fields.
Validation requires a known scalar policy, agreement with the ordinary dataset
ID or frozen instance row, and the corresponding `rng.streams` structure.
Indirect storage in `parameters` or `provenance` is explicitly insufficient.

## Rendering

The second re-audit found no rendering defect. The canonical source was
rerendered after these architectural revisions and the completed MathJax DOM
was checked at desktop, tablet, and mobile widths. Final artifact hashes and
responsive evidence are recorded in the auditor handoff.

## Remaining limitations

- This remains an architecture plan; no proposed package API has been
  implemented.
- The exact G2/G4 operation sequences were checked against current source but
  have not yet been reproduced through the proposed materializer.
- The full specification digest increases ordinary dataset-ID length; the plan
  favors collision resistance and auditable identity over compact IDs.
- The G7 empty-subset exclusion is a deliberate compatibility break. The
  required consumer audit has not yet been run.
- Cross-environment numerical tolerances still require empirical calibration.
- Another independent re-audit is required before the plan can be called
  accepted.
