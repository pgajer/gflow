# Phase 4 support-family ownership

The audit rule is intentionally strict: a utility is not a `gflow` public API
merely because it is useful. Public support must have a named core workflow,
a stable contract, and focused tests.

| Source family | Phase 4 disposition | Named owner / destination |
|---|---|---|
| `stats_utils.R` | de-export; retain only private callers while deletion is staged | graph/flow-specific pieces remain private in `gflow`; generic transforms and Bayesian summaries are candidates for a future `gstats` package |
| `plot_utils.R` | keep only `plot3D.plain.widget()`, `plot3D.cont.widget()`, and `plot3D.cltrs.widget()` public; de-export aliases and miscellaneous plot helpers | `gflow` point-selection and exploration workflow |
| `synthetic_data_utils.R` | de-export all generators and S3 helpers | test-only cases move to `tests/testthat/helper-*`; reusable distribution generators are candidates for `gstats`; graph fixtures belong in `dgraphs` |
| `grids.R` | de-export | generic spatial/grid construction belongs in `dgraphs` when used for graph construction; otherwise delete or move to `gstats` |
| `random_sampling.R` | de-export | generic samplers are candidates for `gstats`; core callers use them privately meanwhile |
| `preprocess_matrix.R` | de-export and make the file the single private implementation of winsorization/robust-z helpers | future `gstats` preprocessing module; `gflow` retains private callers |
| `hist_utils.R` | de-export | delete after downstream audit; base graphics already owns the primitive |
| `boxcox_utils.R` | de-export | private legacy `fassoc*` support in `gflow`; the archive has its own private implementation in `gflowx` |
| `divergences.R` | de-export | future `gstats` distance/divergence module; graph-specific summaries should call a dependency explicitly |
| `wasserstein_dist.R` | de-export | use `transport` directly or move inference wrappers to future `gstats` |

## Evidence

- Caller counts were computed over all `R/*.R` files, excluding each
  definition file, and over `tests/`.
- Most former exports had no package caller and no focused test.
- Functions with one core caller are now private implementation details.
- The three retained widget APIs are used by `select.pts()` and covered by
  `test-plot3d-widget-aliases.R`.
- `winsorize.vec()` and `robust.zscore.vec()` previously had duplicate
  top-level definitions in `preprocess_metabolon_twoplane.R`; Phase 4 removes
  those copies and uses the private implementations in `preprocess_matrix.R`.

## Protected call graph

No basin construction, basin conversion, extrema computation, or protected
`gflowx` hook is changed by this phase. Support functions called by protected
code remain defined privately even when they are no longer exported.

## Follow-up deletion gate

De-exporting is complete in this phase. Actual deletion or cross-package
relocation requires a downstream repository search and, for a new `gstats`
package, an explicit package-creation decision. Until then, private definitions
preserve current internal behavior without presenting them as maintained
`gflow` APIs.
