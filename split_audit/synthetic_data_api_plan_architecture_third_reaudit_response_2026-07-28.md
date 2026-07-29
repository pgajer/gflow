# Synthetic Data API Plan Architecture Third Re-audit Response

Response date: 2026-07-28  
Role: architecture-plan author and implementation worker  
Third re-audit:
`split_audit/synthetic_data_api_plan_architecture_third_reaudit_2026-07-28.md`

## Disposition

The third re-audit accepted the architecture with two nonblocking comments.
Both are carried forward:

- the canonical contract now requires path-specific validation of every
  policy-derived seed before any random draw; and
- the report is regenerated as a rolling report against the current `gflow`
  revision before implementation handoff.

No accepted architecture decision was reversed.

## Derived-seed validation

Legacy seed offsets are computed in double precision and validated as finite,
nonmissing whole numbers in R's representable integer range before RNG state is
changed. Validation is path-specific: unused frame or response branches do not
require derived seeds. Implementation tests cover maximum-integer overflow,
explicit frame seeds, and zero-noise branches.

## Implementation transition

Architecture acceptance closed the design-audit phase. Implementation then
proceeded in separate, pushed commits:

- `geosmooth` `232e412` — canonical specification, dataset, geometry,
  materializer, G1--G7 wrappers, fixtures, and parity tests;
- `geosmooth` `2ffa515` — the shared SSRHE recipe suite and 846-case parity
  fixture;
- `geosmooth` `b9beaf4` — normalized registry tables, manifests, regeneration,
  and foreign-key/hash tests;
- `geosmooth` `12e981a` — finite-design occupation-mixture truth;
- `geosmooth` `6396af1` — circle and trefoil geometries, deterministic grid
  sampling, common plotting, the family disposition ledger, guardrails, and
  vignette; and
- `trend_filtering` `987822c` — the primary SSRHE helper converted to thin
  `geosmooth` registry adapters.

All `test-synthetic-*.R` files pass, including the frozen G1--G7 and SSRHE
fixtures. The synthetic vignette renders successfully. `make check-fast`
cannot pass its default dependency gate because optional `genlasso` is absent.
With `_R_CHECK_FORCE_SUGGESTS_=false`, the check completes with two warnings
and five notes that predate or lie outside the synthetic tranche: unavailable
or non-CRAN local dependencies, bundled ANN object files, `.github` and
nonstandard audit directories, one unimported `tail`, an older lost-brace Rd
diagnostic, and an older `normalize.density` usage/documentation mismatch.

The implementation is not a claim that every Phase 6--9 family has migrated.
The committed disposition ledger explicitly marks noisy predictor curves,
mixed Gaussian/line data, star-response assembly, spline/bump truths, secondary
seed-policy consumers, and ambiguous torus-knot formulas as partial, pending,
or undecided.
