# Canonical Basin Provenance Protected-Surface Authorization

Date: 2026-07-28

Scope: explicit canonical basin-complex implementation and implementation-audit
remediation requested by the package owner.

Baseline commit:
`92a61c086f2fa1fa77223edfb02b74a1be3f1a28`

Initial implementation commit:
`5567e11f4904c50fb5829ae04f322a408ce571f3`

## Adjudicated Protected Changes

### `R/basin_complex.R`

Authorized changes append provenance and external-vertex-identity arguments
after the complete legacy positional signature, retain internal integer
indices, add external-ID companion fields, attach immutable build/runtime
identity, and strengthen `summary.basin_complex()` with independently ranked
maximum and minimum basin tables.

These changes are within the protected constructor's intended ownership. They
do not remove or reorder legacy arguments, change the fixed trajectory-flow
backend selected by existing callers, or weaken schema validation.

### `.apply.basin.refinements`

Authorized changes mark membership-allocated mass stale when a refinement
completes without recomputing that allocation and refresh external-ID companion
fields after support membership changes.

This is necessary to prevent ranking by a stale allocation and to keep the
canonical tables internally and externally aligned after refinement.

## Ownership-Ledger Action

The public `get.gflow.build.identity` export is package-owned canonical basin
infrastructure and must receive an explicit ownership row. The protected
surface and ownership ledger are regenerated with
`tools/build_cleanup_ledger.R` only after the changes above have passed focused
tests. `make audit-cleanup-boundary` remains the acceptance check; replacing
hashes without this written adjudication is not sufficient.

The regenerated ledger also reveals the already-exported and independently
tested `compute.tube.lens.corridor()` as `CORE-ANALYSIS`. It is therefore added
to the maintained non-protected release list; this records the export requested
before the canonical-basin work rather than introducing another API.

## Exclusions

No retirement, ownership transfer, breaking API migration, native-symbol
change, or downstream-call removal is authorized by this record.
