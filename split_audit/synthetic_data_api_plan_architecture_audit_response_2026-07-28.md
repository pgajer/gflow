# Synthetic Data API Plan Architecture Audit Response

Response date: 2026-07-28  
Role: architecture-plan author  
Audit:
`split_audit/synthetic_data_api_plan_architecture_audit_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`

## Disposition

All audit findings were accepted. No finding was deferred or rejected. The
canonical plan and generated HTML were revised. This response records the
changes; it does not assert that the independent audit verdict has changed.

## Finding resolutions

### DGP-1 — Legacy one-dimensional sampling and noise semantics

Resolved in the sampling/response contracts, absorption ledger, legacy RNG
policy, Phase 1 parity grid, Phase 5, and acceptance criteria.

The plan now specifies:

- ascending row order for V1, V2, and V3;
- deterministic V3 allocation with
  `floor(0.42 * n)` left rows and the remainder right;
- explicit multinomial, fixed-proportion, and fixed-count allocation modes;
- explicit allocation rounding policies;
- `minimum.outliers = 1L` for V3, yielding
  `max(1L, round(0.04 * n))`;
- contractual V3 RNG draw order: Laplace errors, contamination-row sampling,
  then Gaussian contamination draws;
- a parity grid including `n = c(1L, 2L, 10L, 11L, 160L)`, three replicates,
  and all S/V combinations.

### DGP-2 — G4 and G7 dataset representation

Resolved by making latent and dimension fields explicitly nullable and adding
heterogeneous-region metadata.

The common object now contains:

- `intrinsic.dim.by.region`;
- `codimension.by.region`;
- `latent.mask`;
- geometry-specific nullable-field invariants.

G4 uses a two-column latent matrix with masked `NA_real_` entries for the
absent A-stratum coordinate and records dimensions A=1 and B=2. G7 uses
`latent = NULL`, records simplex part count, active parts, and the interior and
structural-zero face dimensions. Transitional `U` and `d` aliases preserve the
current outputs.

### ART-1 — Recipe versus instance identity

Resolved by separating reusable recipes from frozen instances.

- `recipes.csv` contains component combinations and no `n`, seed, or checksum.
- `instances.csv` contains instance ID, recipe ID, `n`, seed, RNG policy,
  checksum ID, and fixture path.
- `materialize.synthetic()` creates a derived non-frozen ID from recipe/spec
  identity, `n`, seed, and RNG policy.
- `materialize.synthetic.instance()` consumes all frozen values from the
  instance row and permits no override.
- Clustered sampling fixes `n = cluster.count * observations.per.cluster`.
- Graph-vertex sampling without replacement fixes `n` to vertex count.

### ART-2 — Canonical checksum scope

Resolved by adding payload contract `synthetic-content-v1`.

The plan now fixes:

- payload field order;
- type and missing-value normalization;
- component and graph canonicalization;
- R XDR serialization version 3;
- SHA-256 computation;
- same-platform scope fields;
- elementwise cross-platform tolerances;
- orthogonal-invariance comparison for random-frame predictors;
- CRAN restrictions against platform-sensitive QR/BLAS fixture hashes.

### ART-3 — Omitted repository provenance

Resolved by recording source paths and render-time revisions for:

- `grip` at `c5c85420f5bbebc82c53dfa9b27cacaae8fd284a`;
- `vaginal_community_trajectory_types` at
  `258b0113d625cbb57b3300721f7d583d178338dc`.

The source list now includes `grip` family/graph files and the
occupation-density application generator.

### IMP-1 — Collision with `stats::simulate()`

Resolved by replacing the ordinary `simulate.synthetic()` proposal with:

```r
materialize.synthetic()
materialize.synthetic.instance()
```

The report states that `stats::simulate()` is an existing S3 generic with
formals `object, nsim, seed, ...` and that the new materializers are ordinary
exported functions, not S3 methods.

### IMP-2 — Frame RNG ownership

Resolved by making every component constructor draw-free.

`synthetic.quadform()` stores one of:

- a canonical-frame policy;
- a random-orthonormal-frame policy;
- a supplied frame matrix.

Only the materializer draws a random frame, using the `geometry.frame` stream.
Legacy G2 uses `rng.policy = "legacy"` and seed policy
`geosmooth.g1.g7.v1`, so the materializer derives the historical
`seed + 500000L` frame state.

### IMP-3 — Missing constructor contracts

Resolved by adding signatures and validation/compatibility semantics for:

- sphere cap, helix, torus patch;
- point, segment, rectangle, stratified, and point-line geometry;
- simplex and graph geometry;
- stratified and Dirichlet-with-zeros sampling;
- occupation-mixture and graph-signal truth;
- Bernoulli acceptance behavior.

The contracts state intrinsic/ambient dimensions, frame behavior, row
allocation, region/face metadata, graph alignment, tie-breaking, and legacy G4,
G7, occupation, and graph-signal behavior.

### REND-1 — Hidden absorption-ledger column

Resolved in CSS and verified in Chromium.

- The overriding `overflow: hidden` was removed.
- Table overflow is `auto`.
- Table and inline-code tokens may wrap.
- Code blocks and MathJax displays have bounded local overflow.
- The ownership-card HTML indentation was corrected; it now renders as four
  cards rather than escaped HTML code blocks.

At viewport widths 1,440, 768, and 390 pixels:

- page `scrollWidth == clientWidth`;
- all three tables report `overflow-x: auto`;
- all three header columns are present and reachable;
- the absorption ledger's “Result after migration” column is visible;
- four ownership cards and zero ownership code blocks are present.

Full-page and ledger screenshots were generated under:

```text
/tmp/gflow-final-responsive-qa.4XLbwv/
```

These temporary screenshots are validation artifacts and are not package
sources.

## Final report artifacts

```text
Rmd
  bytes: 80,608
  SHA-256: 19e04c598523ead9e3a9e61f1cfd2169e904d48394b0a960d8479f1c177ada24

HTML
  bytes: 1,342,153
  SHA-256: a7b0332e0b71f527d098e6ccb8ec2f6ee90b6dca1b2193f0af40e0ac7de7c17a
```

The R Markdown was rendered after the final source edit. Structural HTML
checks and desktop/tablet/mobile Chromium checks were then run against that
generated HTML.

## Re-audit status

An independent re-audit subsequently returned **Revise before acceptance**.
Its residual findings are addressed in:

`split_audit/synthetic_data_api_plan_architecture_reaudit_response_2026-07-28.md`

This original response remains the disposition record for the first audit; the
re-audit response supersedes its artifact hashes and remaining-limitations
snapshot.

## Remaining limitations at the time of this response

- This remains a design plan; none of the proposed package API is implemented.
- No legacy parity fixture has yet been generated by the proposed API.
- The checksum serializer and tolerance comparator are specified but
  unimplemented.
- The responsive rendering checks establish accessibility of the full table
  content, not semantic correctness of the future package code.
- Independent re-audit had not yet occurred when this response was written.
