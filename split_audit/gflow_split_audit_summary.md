# gflow Split Audit SA0/SA1 Summary

Generated: 2026-06-03

Branch: `codex/gflow-split-audit-sa0-sa1`

## What Was Done

- Created the SA0 boundary contract:
  `split_audit/gflow_split_boundary_contract.md`
- Created a reproducible SA1 inventory script:
  `split_audit/scripts/build_sa1_inventory.R`
- Generated detailed inventory CSV files:
  `split_audit/gflow_split_inventory.csv`
  `split_audit/gflow_split_inventory_summary.csv`
- Created the SA1 inventory narrative:
  `split_audit/gflow_split_inventory.md`

No source files were moved and no package APIs were changed.

## Proposed Split

- `gflow`: graph construction, graph/geodesic machinery, gradient-flow cells,
  Morse-Smale graph objects, basins, extrema, trajectories, and geometry
  diagnostics.
- `geosmooth`: conditional-expectation and geometric smoothing/regression
  methods, including kernel local polynomial smoothing, MALPS, LPL-TF,
  SLPLiFT/S-LPL-TF, SSRHE, graph trend filtering, metric graph low-pass, and
  harmonic response smoothers.
- `gflowx`: experimental, retired, application-specific, microbiome-specific,
  association-test, clustering, visualization, or workflow-heavy methods that
  should be preserved but not carried by the CRAN-facing core.

## Inventory Counts

| Proposed package | Asset count |
|---|---:|
| `geosmooth` | 65 |
| `gflow` | 174 |
| `gflowx` | 63 |
| `shared` | 65 |
| `undecided` | 63 |

The current package surface has 493 exports and 206 S3 registrations.

## Main Decision Points for SA2

SA2 should focus on dependency closure rather than further naming discussion.
The most important questions are:

1. Which graph/kNN/geodesic APIs must `geosmooth` import from `gflow`?
2. Which local-chart, support-size, kernel, weighted least-squares, and CV
   helpers are truly shared across `gflow`, `geosmooth`, and `gflowx`?
3. Which native symbols are called by the proposed `geosmooth` R files?
4. Can `kernel.local.polynomial.cv` be extracted as the first pilot without
   copying graph infrastructure?
5. Which currently undecided files are required by the pilot extraction?

## Recommended Next Step

Proceed to **SA2: Dependency Graph** on a follow-up branch after reviewing the
SA0/SA1 documents.  SA2 should generate a call/dependency map for the proposed
`geosmooth` assets, starting with the pilot candidate
`kernel.local.polynomial.cv`.

The physical package split should wait until SA2 and SA3 are complete.
