# Synthetic Data API Plan Architecture Re-audit Response

Response date: 2026-07-28  
Role: architecture-plan author  
Re-audit:
`split_audit/synthetic_data_api_plan_architecture_reaudit_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`

## Disposition

All eight residual findings were accepted. None was deferred or rejected. The
canonical plan and generated HTML were revised. This response records the
changes and does not assert a new independent verdict.

## Finding resolutions

### DGP-R1 — Exact random algorithms

Resolved by making legacy algorithm identity and draw order contractual.

- V2 is fixed to `pnorm` bounds, one probability-scale `runif()` draw,
  `qnorm()`, and sorting.
- V3 draws the fixed left count first, the right count second, concatenates,
  and sorts.
- V3 Laplace noise uses the exact uniform inverse transform with `log1p()`.
- Disk sampling uses radial and angular uniform vectors in that order.
- Clustered sampling uses column-major center draws followed by column-major
  Gaussian offsets; the response stream draws cluster effects before
  residuals.
- Torus sampling draws all first-angle values before all second-angle values.
- Dirichlet sampling uses a column-major gamma matrix, row normalization,
  structural-zero application, and a second row normalization.

The relevant component signatures now carry versioned algorithm IDs. Phase 1,
Phase 3, Phase 5, and the acceptance criteria require fixtures for these exact
sequences, not only distributional equivalence.

### DGP-R2 — Population versus finite-design truth

Resolved by defining two truth scopes:

- `"population"` for a fixed \(m(x)=E[Y\mid X=x]\); and
- `"finite.design"` for
  \(m_i(X_{1:n})=E[Y_i\mid X_1,\ldots,X_n]\).

The dataset contract now includes `truth.scope` and `truth.estimand`.
`sample.max` Gaussian-mixture truth and design-maximum occupation truth are
explicitly finite-design estimands, and their realized normalizers must be
recorded. A population version requires a fixed recipe-level normalizer and a
distinct component version.

### DGP-R3 — Declared and observed regions

Resolved by separating `declared.regions` from `observed.regions`.
Heterogeneous dimension vectors are keyed by declared regions, while row labels
must be a subset of them.

The maintained G4 recipe requires at least one A and one B row. It computes
legacy counts before any draw and rejects a zero-row stratum with an explicit
validation error. G7 permits endpoint zero fractions; both face labels remain
declared while an absent label is omitted from `observed.regions`.

### ART-R1 — Truth evaluator identity and checksums

Resolved by removing arbitrary closures from all materializable
specifications.

`synthetic.truth.named(id, parameters)` resolves a versioned evaluator registry
record. Functions, closures, environments, external pointers, and unevaluated
language objects are rejected. Maintained compatibility defaults map to frozen
evaluator IDs; arbitrary historical `truth.fn` customizations are outside the
exact-parity surface and receive a lifecycle warning followed by an explanatory
error.

The plan also defines `synthetic-specification-v1`, the canonical payload used
for ad hoc specification IDs. It includes evaluator ID/version and serializable
parameters, so the ID and dataset-checksum contracts no longer disagree.

### ART-R2 — Random-frame and environment scope

Resolved by defining both sides of the comparison.

- Same-environment exact checksums include the realized frame matrix.
- Cross-environment scientific comparison does not compare random-frame
  entries elementwise when ambient orientation is a nuisance.
- Each frame must be orthonormal and must reconstruct its own predictors.
- Predictor pairwise squared distances are compared, while latent coordinates,
  truth, response, offsets, and design statistics use numeric tolerance.
- Orthogonal equivalence is permitted only for latent-coordinate or explicitly
  orthogonally invariant truth.

The exact-check environment fingerprint now includes the full R build and
platform, architecture, endianness, operating system, compiler, RNG details,
BLAS/LAPACK paths and hashes, math runtime, package versions, registry schema,
and evaluator-registry version.

### IMP-R1 — SSRHE seed-policy identifiers

Resolved in both worked examples:

- flat uses `ssrhe.flat.v1`;
- quadform uses `ssrhe.quadform.v1`.

Both values now match the implemented-policy names in the RNG contract.

### IMP-R2 — Geometry and sampling support

Resolved by assigning distinct responsibilities:

- geometry owns the hard admissible latent domain;
- sampling owns its actual distribution and support.

`synthetic.spec()` validates sampler-support containment. Frozen G3b, G3c, and
G3d require support equality; new recipes may deliberately use a narrower
support. Quadform, simplex, stratified, and graph compatibility rules are also
specified, and the materializer revalidates realized rows.

### IMP-R3 — Uniform-box order policy

Resolved by limiting `"ascending"` to the one-dimensional interval and
truncated-normal samplers. `uniform.box` retains `"draw"` and
`"ascending.first.coordinate"` exactly as shown in its signature.

## Rendering

The re-audit found no residual rendering defect. The canonical source was
rerendered after the architectural revisions and the completed MathJax DOM was
checked at desktop, tablet, and mobile widths. Final artifact hashes and
responsive evidence are recorded in the auditor handoff.

## Remaining limitations

- This is still an architecture plan; the proposed package API is not yet
  implemented.
- Legacy fixtures, the serializer, evaluator registry, and comparators remain
  implementation deliverables.
- Cross-environment tolerance values remain subject to empirical calibration,
  although the comparison policy and environment fingerprint are now explicit.
- A new independent re-audit is still required before the plan can be called
  accepted.

## Second re-audit status

A second independent re-audit subsequently returned **Revise before
acceptance**. Its residual findings are addressed in:

`split_audit/synthetic_data_api_plan_architecture_second_reaudit_response_2026-07-28.md`

This response remains the disposition record for the first re-audit; the
second response supersedes its artifact and remaining-limitations snapshot.
