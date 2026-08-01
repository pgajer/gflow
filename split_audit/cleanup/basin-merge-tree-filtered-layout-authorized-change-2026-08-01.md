# Basin Merge-Tree Filtered-Layout Protected-Surface Authorization

Date: 2026-08-01

Scope: explicit package-owner request to implement Phase 1 of the accepted
`gflowui` Basin Analysis Panel redesign.

Baseline commit:
`24a671c4927df6ab6e5ac10361aecfd87cfaa0cb`

Normative architecture:
`gflowui/dev/basin_merge_tree_adaptive_initial_filtering_spec_2026-07-31.md`
at SHA-256
`d017677a5a236ea2e772e76124f106f36774359dd2ef6b64a4767e65f2102884`.

## Adjudicated Protected Changes

### `R/basin_merge_tree_public.R`

The package now exports the pure
`get.basin.merge.tree.layout()` computational accessor. It validates one
canonical direction/component, accepts either the complete component or an
explicit canonical-basin-ID subset, optionally adds canonical ancestors, and
returns the selected canonical branch/event tables, component root,
crossing-free restricted leaf order, and branch/event coordinates without
drawing.

The accessor rejects empty, duplicate, unknown, mixed-direction,
mixed-component, component-mismatched, and non-ancestor-closed requests. It
also validates finite canonical vertical values and nonnegative persistence.
Ancestor closure uses canonical parentage only; no graph-induced tree or basin
reconstruction is introduced.

`plot.basin.merge.tree()` delegates its computational inputs to the public
accessor and adds `basin.ids` and `close.ancestors` after the complete legacy
positional signature. Existing positional arguments, canonical construction,
elder survival, parentage, events, and vertical values are unchanged.

### Tests and documentation

Focused tests cover complete, restricted, closure-added, empty, one-branch,
unknown, mixed-direction, mixed-component, component-mismatch, invalid-value,
restricted-order, plotting-parity, and positional-compatibility cases.
Roxygen source documents the new public API and the two additive plotting
arguments.

## Ownership-Guardrail Action

After focused tests and documentation generation pass, regenerate
`split_audit/cleanup/api-ownership.csv` and
`split_audit/cleanup/protected-basin-surface.txt` with
`tools/build_cleanup_ledger.R`. The updated rows and hashes record this
explicitly authorized additive protected-surface change; they do not relax
protected-file coverage.

## Exclusions

No canonical basin reconstruction, plateau handling, density-value elder
rule, basin identity, assignment, membership, persistence calculation,
component membership, ownership transfer, retirement, or native-symbol
behavior is changed.
