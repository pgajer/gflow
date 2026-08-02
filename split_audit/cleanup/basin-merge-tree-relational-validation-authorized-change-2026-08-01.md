# Basin Merge-Tree Relational Validation Protected-Surface Authorization

Date: 2026-08-01

Scope: correction of blocker `P1-IMP-1` from the independent audit of the
Phase 1 filtered merge-tree layout implementation.

Phase 1 implementation base:
`eaaa65ea0cba34bd804db2e1c92a4a21648be8d5`

Normative Revision 9 specification SHA-256:
`d017677a5a236ea2e772e76124f106f36774359dd2ef6b64a4767e65f2102884`

Historical audit:
`/Users/pgajer/current_projects/gflow/split_audit/basin_analysis_phase1_filtered_layout_implementation_audit_2026-08-01.md`

## Adjudicated Protected Changes

### `R/basin_merge_tree_public.R`

The canonical `basin.merge.tree` validation boundary now checks exact
relational agreement between every nonroot branch and its unique matching
merge event. The validation requires:

- branch parent equals event survivor;
- event birth equals losing-branch birth;
- event death and merge level equal losing-branch death;
- event persistence equals losing-branch persistence;
- maximum persistence equals `birth.level - death.level`;
- minimum persistence equals `death.level - birth.level`; and
- each component root has the finite component-field boundary as its death
  level.

Parent and survivor branches must belong to the same requested direction and
graph component. These checks reject contradictory mutable R objects before
complete or restricted layout construction and before plotting.

The change does not alter valid canonical objects or their branch tables,
event tables, parentage, survivor identity, restricted leaf order, or
coordinates.

### Tests

Focused adversarial tests independently corrupt branch parent, event survivor,
event merge, event birth, event death, event persistence, direction-specific
branch death/persistence relations, and the component-root death convention.
Every case is exercised for maximum and minimum trees through direct
validation, complete extraction, restricted extraction, and plotting. Valid
complete and restricted outputs remain covered for both directions.

## Ownership-Guardrail Action

After the focused public merge-tree suite passes, regenerate:

- `split_audit/cleanup/api-ownership.csv`; and
- `split_audit/cleanup/protected-basin-surface.txt`

with `tools/build_cleanup_ledger.R`. The regenerated protected-file hash
records this explicit validation correction without weakening protected
coverage.

## Exclusions

No canonical construction, plateau handling, density-value elder rule,
parentage, persistence calculation, component membership, assignment,
selection policy, or rendering algorithm is changed.
