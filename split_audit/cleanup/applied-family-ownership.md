# Phase 5 applied-family ownership

This table is the forward ownership contract for applied and generic-analysis
families that do not belong in the public `gflow` API. A function may remain
defined privately during a destination-repository conflict, but it must not be
exported or registered as an S3 method.

| Family | `gflow` action | Destination | Relocation status |
|---|---|---|---|
| Compositional transforms | De-export `arcsin.sqrt.transform()`, `clr.transform()`, and `clr.transform.matrix()` | `microbiome.utils` | Destination worktree has unrelated changes; relocate from a clean worktree |
| Phylotype selection | Remove the public `print.phylotype_selection()` registration; keep the selection implementation private until moved | `microbiome.utils` | Destination worktree has unrelated changes; relocate from a clean worktree |
| Metabolon two-plane processing | De-export preprocessing and remove its plotting S3 registrations | external metabolomics analysis repository | No maintained package owner; retain privately only until the source analysis is identified |
| CST summaries and colors | De-export CST summaries and weighted analysis | `gcstflow` | Ready for coordinated relocation |
| Generic two-factor/categorical analysis | De-export | external analysis repository | No core-object interoperability contract |
| MGCP visualization/export | De-export and remove plotting registrations | `gflowui` | Ready for coordinated relocation |
| Partition heatmaps | De-export | `gflowui` | Ready for coordinated relocation |
| Distance-diagnostic reports | De-export report generator and plotting registrations | external analysis repository | No core-object interoperability contract |
| Concordance reports | Remove public plotting/printing registrations | external analysis repository | Underlying analysis remains private pending source-analysis identification |
| Dataset-specific trajectory reports | Remove public hurdle-association plotting registration | external analysis repository | Data preparation and plotting remain private pending source-analysis identification |
| ISA and consensus clustering | Keep private; no public API | external analysis repository | Delete after final unqualified downstream search |
| Application clustering | De-export application-specific clustering entry points | external analysis repository | `cluster.graph.louvain()` and `congruence.with.labels()` remain temporary public migration exceptions because `CT_clearance` has qualified calls |
| PHATE | Keep `phate()`, `phate.core()`, and `phate.embed()` temporarily public; admit no new callers | `dgraphs` | Qualified calls remain in `ZB` and `graph_modularity`; destination worktree has unrelated changes, so relocate as a parity-tested closure |
| Diffusion and potential pseudotime | Keep the two pseudotime entry points temporarily public; de-export the generic transition builder | `dgraphs` | Qualified calls remain in `ZB` and `cell_cycle`; destination worktree has unrelated changes |
| Potential-metric graph adapters | Retain only adapters whose input or output is a protected `gflow` object; de-export generic constructors | `dgraphs` | Destination worktree has unrelated changes |
| Generic trajectory distances | De-export generic metrics; retain no report-specific API | `dgraphs` | Destination worktree has unrelated changes; relocate after graph-metric parity tests |

## Temporary public migration exceptions

The following names remain public solely to avoid breaking known qualified
downstream calls while the `dgraphs` destination worktree is unavailable:

- `phate()`, `phate.core()`, and `phate.embed()`;
- `compute.diffusion.pseudotime.sparse()` and
  `compute.potential.pseudotime.sparse()`; and
- `cluster.graph.louvain()` and `congruence.with.labels()`.

These are not approved additions to the long-term `gflow` surface. New internal
or downstream callers are prohibited by the ownership-ledger guard. They move
as parity-tested closures in a clean `dgraphs` migration, after which their
`gflow` exports are retired.

## Enforcement

- `split_audit/cleanup/api-ownership.csv` is regenerated after every family.
- Removed names are recorded in `retired-exports.txt`.
- `make audit-cleanup-boundary` fails if a removed name returns, a current
  namespace entry lacks an owner, or the protected construction surface
  changes without an explicit protected-surface refresh.
- Qualified downstream calls are captured in the ownership ledger. Unqualified
  calls from scripts that attach `gflow` must still be searched before physical
  deletion.
