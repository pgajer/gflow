# gflow implementation and verification map

Start with the public contract, then follow the relevant row. Existing modules
and numerical kernels remain in place; generated bindings are not editing targets.

| Responsibility | Source | Meaningful regression coverage |
|---|---|---|
| Canonical constructor, input validation, parameter resolution and schema | `R/basin_complex.R` | `test-basin-complex-phase-b.R`: graph/field checks, method applicability, table schemas |
| Backend adapters and retained historical reconstruction | `R/basin_complex_adapters.R`, `R/basin_complex_converters.R` | `test-basin-complex-phase-c.R`, `phase-g.R`, `phase-h.R`: overlapping supports, distinct modulation modes, lifecycle boundary |
| Connected exact plateaus and trajectory forest | `R/basin_complex_plateau_flow.R` | `test-basin-complex-plateau-flow.R`: plateau representatives, tie behavior and terminal flow |
| Exact graph merge trees | `R/basin_complex_merge_tree.R` | `test-basin-complex-phase-d.R`, `test-basin-merge-tree-public.R`: plateau topology, extrema orientation, cuts |
| Ordered refinement stages | `R/basin_complex_refinement.R` | `test-basin-complex-phase-f.R`: filtering, retained/raw distinction, unsupported settings |
| Support summaries, identities, display | `R/basin_summary.R`, `R/basin_identity.R`, `R/basin_complex.R` | `test-basin-complex-summary-identity.R`, `test-basin-displays.R`: IDs, mass, coverage, failed-result labels |
| Canonical flow association | `R/gfassoc_canonical.R`, `R/gfcor.R`, `R/gfassoc_utils.R`; `src/gfassoc_*` | `test-canonical-association.R`: independent archived-support fixture, polarity formula, reversed fields, uncovered vertices |
| Local edge-difference association | `R/lcor.R`; `src/lcor.cpp`, `src/lcor_r.cpp` | `test-local-association-science.R`: independent cosine, anisotropic star, affine changes, logratio, selected pairs; `test-lcor-hop-radius.R`: neighborhoods |
| Supplied-draw summaries and permutation inference | `R/lcor_with_posterior.R`, `R/permutation_test_lcor.R` | `test-lcor-posterior.R`, `test-permutation-contract.R`: all draws/identities, independent summaries, RNG restoration, grouped null |
| External graph representation | `R/dgraphs_boundary.R` | `test-dependency-contracts.R`, `test-basin-complex-phase-g.R`: components, hop lengths, overlap cells; run with published and development dgraphs |

The phase-named tests retain their filenames; their descriptive `test_that()`
names identify the scientific contracts. Search by subject before adding tests.

## Generated files and ownership

Edit roxygen in `R/`, then run `make document` to regenerate `NAMESPACE`, Rd and
catalog links. `tools/update_guide_links.R` resolves help aliases, including aliases
that share a topic. `make build` also refreshes the source manifest and includes all
four rendered guides. `make audit-api-guide audit-final-acceptance` checks catalog,
namespace and ownership. Run `make check` for the full package check.

- [API ownership](../split_audit/cleanup/api-ownership.csv)
- [Native ownership](../split_audit/cleanup/native-symbol-ownership.csv)
- [Dependency ownership](../split_audit/cleanup/dependency-ownership.csv)
- [Protected basin surface](../split_audit/cleanup/protected-basin-surface.txt)
- [Document locations](../DOCUMENT_LOCATIONS.md)

When an authorized basin change alters protected fingerprints, regenerate with
`cleanup.build.protected.surface()` from `tools/cleanup_ledger_lib.R`, inspect the
exact changed entries, and commit the record with the source change. Do not
weaken the guardrail to make a check pass. Generic graph algorithms belong to
dgraphs; optional browser rendering belongs to ivue; smoothing belongs to gflowx.

## Resource contracts

For n vertices, p response fields, q other fields and B supplied draws, doubles
occupy eight bytes. These are output-size calculations, not measured peak memory.

| Calculation | Numeric storage before inputs and temporary copies |
|---|---|
| lcor vector/vector | 8n bytes |
| lcor vector/matrix | 8np bytes |
| identical matrices, unique unordered pairs | 8n p(p−1)/2 bytes |
| distinct matrices, all pairs | 8npq bytes |
| explicit K pairs | 8nK bytes plus one temporary vertex result |
| posterior mean, SD, lower, upper | 32np bytes plus a working 8nB matrix |
| returned posterior samples | an additional 8npB bytes |
| permutation test without saved statistics | feature counts, observed values and a working local result; no B-by-p statistic matrix |
| saved permutation statistics | an additional 8Bp bytes |

At n=10,000 and p=q=100 the asymmetric lcor result alone is 800 MB (decimal).
Use `pairs` in small blocks, persist each result with its pair indices, and discard
it before computing the next block. The API preserves supplied order, duplicate
pairs and names; it never silently samples pairs. Explicit common logratio epsilon
is needed for comparisons across dispatch paths. Native vector-matrix adaptive
epsilon uses the whole matrix minimum, whereas pairwise paths adapt per field.

`construct.madag()` already exposes `max.trajectories.per.cell` and
`enumerate.trajectories = FALSE`. Keep those limits; counting reachable paths does
not require enumerating all of them. Expanded hop graphs and parallel worker copies
can dominate memory even when the returned result is small.

## Interrupt audit and limits

A source audit of reachable local association, flow association, basin trajectory,
RTCB and MADAG kernels found no in-loop R interrupt polling. Defining
`CHECK_INTERRUPT` in `gflow_macros.h` does not make those kernels interruptible:
it has no call sites there. The active `R_CheckUserInterrupt()` occurrence in
`src/linf_simplex_knn.cpp` is outside these core routes; a geodesic wrapper has
only a commented occurrence. R can process interrupts between serial pair or
permutation calls, but a single large native call may finish before cancellation.
This is a documented current limit, not a claimed responsiveness benchmark.

Do not insert R API calls into OpenMP worker loops. Adding native interruption
requires an exception-safe boundary and main-thread polling design; this review
preserves the existing kernels and provides bounded pair execution instead.
A stopped block sequence is incomplete: consumers must check the saved completed
pair indices against the requested list before reporting a complete analysis.
