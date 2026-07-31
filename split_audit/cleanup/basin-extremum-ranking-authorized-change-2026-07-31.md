# Basin Extremum Ranking Protected-Surface Authorization

Date: 2026-07-31

Scope: explicit package-owner request to remove automatic ranking from the
`gflowui` Basin Inspector, expose peak-value ranking, and make mass the
occupation-density default.

Baseline commit:
`5698a1fb`

## Adjudicated Protected Changes

### `R/basin_complex.R`

The public `summary.basin_complex()` ranking vocabulary is extended with
`rank.by = "extremum.value"`. This is an additive API change. Existing ranking
values, positional arguments, defaults, construction methods, and basin
assignments are unchanged.

### `.basin.summary.measure.availability`

Finite extremum values are accepted as a usable explicit ranking measure.
Unlike mass and support measures, signed and all-zero extremum values remain
valid because a generic conditional-expectation field may be signed and
deterministic basin-ID tie breaking remains defined.

The existing `"auto"` resolution hierarchy remains limited to canonical mass
and support measures. It cannot silently fall back to extremum value, which
keeps peak ranking an explicit scientific choice.

### `.basin.summary.direction`

Extremum-value ranking is direction-aware: maximum basins are ordered from
highest to lowest representative value, while minimum basins are ordered from
lowest to highest. Canonical basin ID remains the deterministic tie breaker.
All mass and support ranking order is unchanged.

## Ownership-Guardrail Action

The protected surface is regenerated with
`tools/build_cleanup_ledger.R` only after the focused canonical-summary test
passes. The resulting hashes record this reviewed additive summary change;
they do not relax protected-file coverage.

## Exclusions

No basin reconstruction, trajectory-flow, plateau, assignment, refinement,
mass calculation, ownership transfer, retirement, or native-symbol change is
authorized by this record.
