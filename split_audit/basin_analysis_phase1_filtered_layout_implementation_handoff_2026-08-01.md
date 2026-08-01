# Basin Analysis Redesign Phase 1: Implementation Handoff

Date: 2026-08-01

Repository: `/Users/pgajer/current_projects/gflow`

Baseline:
`24a671c4927df6ab6e5ac10361aecfd87cfaa0cb`

Normative Revision 9 specification SHA-256:
`d017677a5a236ea2e772e76124f106f36774359dd2ef6b64a4767e65f2102884`

The implementation commit is the commit containing this handoff. This
document records implementer evidence only; it does not claim independent
audit acceptance.

## Implemented scope

- Added and exported the pure `get.basin.merge.tree.layout()` API.
- Supports a complete selected direction/component or an explicit nonempty
  canonical-basin-ID subset.
- Rejects duplicate, unknown, mixed-direction, mixed-component,
  component-mismatched, and non-ancestor-closed requests.
- With `close.ancestors = TRUE`, adds canonical ancestors and reports only the
  closure-added IDs.
- Validates finite branch birth/death/persistence values, finite event
  birth/death/merge/persistence values, and nonnegative persistence.
- Returns exact selected canonical branch/event tables, component root,
  restricted canonical leaf order, dendrogram inputs, branch/event
  coordinates, and `validation.status = "ok"` without drawing.
- Preserves canonical vertical values and deterministically verifies that a
  restricted leaf order is the restriction of the complete canonical order.
- Refactored `plot.basin.merge.tree()` to consume the public accessor.
- Added `basin.ids` and `close.ancestors` after the complete legacy positional
  plotting signature.
- Kept canonical construction, elder survival, parentage, persistence,
  component membership, and assignments unchanged.

## Files owned by this phase

- `R/basin_merge_tree_public.R`
- `NAMESPACE`
- `tests/testthat/test-basin-merge-tree-public.R`
- `inst/extdata/gflow-code-manifest.tsv`
- `split_audit/cleanup/api-ownership.csv`
- `split_audit/cleanup/protected-basin-surface.txt`
- `split_audit/cleanup/basin-merge-tree-filtered-layout-authorized-change-2026-08-01.md`
- this handoff

The pre-existing modified `AGENTS.md` is excluded from the phase commit.
Generated `man/` files are repository-ignored and are rebuilt from roxygen
source by `make document` and package QA.

## Test evidence

Focused public merge-tree suite:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file(
  "tests/testthat/test-basin-merge-tree-public.R",
  reporter="summary", stop_on_failure=TRUE)'
```

Result: 99 expectations, 0 failures, 0 errors, 0 warnings, 0 skips.

Full source-loaded suite:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_dir(
  "tests/testthat", reporter="summary", stop_on_failure=TRUE)'
```

Result: all expectations passed; 10 existing CRAN-gated/empty-test skips.

Package and ownership QA:

```sh
make document
make audit-cleanup-boundary
make check-fast
make check
git diff --check
```

Both CRAN-style checks completed with one expected incoming-feasibility
warning: local packages `dgraphs` and `grip` are not present in the CRAN or
Bioconductor repositories. Examples, tests, vignettes, PDF/HTML manuals,
namespace checks, and code/documentation consistency passed. Ownership audits
pass for 112 exports and 99 S3 methods.

## Required independent audit

Before Phase 2 begins, independently verify:

1. complete, restricted, closure-added, empty, and one-branch behavior;
2. unknown, mixed-direction, mixed-component, component-mismatch, missing-root,
   and invalid-vertical-value rejection;
3. exact preservation of canonical birth, death, merge, persistence, parent,
   event, survivor, and component values;
4. restricted leaf-order preservation and coordinate determinism;
5. equality of the accessor inputs returned through
   `plot.basin.merge.tree()`;
6. legacy plotting positional compatibility;
7. the protected-surface authorization and regenerated ownership ledgers; and
8. the focused/full/package QA evidence above.

Phase 2 must not begin until this implementation gate receives independent
acceptance.
