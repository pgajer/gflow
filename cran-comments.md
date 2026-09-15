# gflow 0.2.0 candidate — local validation, 2026-09-15

This is a release-readiness record, not a claim of submission. The version is
unchanged. Earlier multi-platform check claims have not been carried forward.

## Current validation

* macOS Apple Silicon, R-devel 4.7.0 (2026-06-24 r90190), Apple clang C++17.
* `make check` (build plus `R CMD check --as-cran`) with published dgraphs
  0.2.0: 0 errors, 0 warnings, 1 NOTE. All four installed vignettes, examples,
  tests, and PDF/HTML manuals passed.
* Tests: 1,003 passed, 0 failed, 0 warnings, 22 skipped. Skips comprise nine
  CRAN-disabled nerve tests, twelve source-only audits excluded from the source
  package, and one empty test. Source ownership, cleanup, S3, and final-acceptance
  audits were run separately through the Makefile.
* Separate serial build with OpenMP compile/link flags explicitly empty:
  installation and both introductory vignette workflows passed; the loaded
  native diagnostic reported `openmp_compiled = FALSE`. This is a macOS serial
  check, not Linux or Windows validation.
* Compilation flags and compiled-code checks passed. Compiler output still
  contains unused-variable warnings in existing native sources; these did not
  produce an R CMD check WARNING.

## NOTE and remaining coordination

CRAN incoming feasibility reports “New submission” and optional dependency
`ivue` not in mainstream repositories. The package source on public GitHub and
local ivue 0.1.0 expose the required scene/layer functions, but no current CRAN
entry was found. Resolve its acceptable release/distribution path before upload.
An initial check also flagged a noncanonical Hmisc documentation URL; the
maintained vignette now links to the canonical package page.

No gflow entry was found in the current CRAN package index or source archive,
and no gflow candidate appeared in the publicly visible incoming queues checked
on 2026-09-15. These observations do not reveal private review correspondence.
Confirm pending submission status with the maintainer's records before upload.

## Follow-up before release

* Reconcile the declared R >= 3.5.0 with required dgraphs' R >= 4.1.0 floor;
  the oldest supported R version has not been tested for this candidate.
* Record an explicit dgraphs minimum: the needed public entry-point definitions
  are identical in archived 0.1.0 and published 0.2.0; runtime checks used 0.2.0.
  The local 0.3.0.9000 development API changes connected-component inputs and
  currently breaks the overlap-cell backend. Coordinate that future change.
* `gfcor()` and `gfassoc.membership()` still require archived
  `basins_of_attraction` objects. The new guides document this boundary; a
  canonical adapter would be a separate API change.
* Obtain current Windows/Linux and oldest-supported-R checks. Historical
  Win-builder, R-hub, rchk, valgrind, and sanitizer reports do not validate this
  candidate. No upload, website publication, or release messages were sent.
