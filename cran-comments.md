## gflow 0.2.0

Resubmission following the earlier 0.1.0 attempts. The Linux loading issue is resolved: gflow has no TBB dependency and uses published dgraphs (>= 0.4.0).

Repeated vignette fits are precomputed with full reproduction code. Scientific results, examples and figures are retained. The vignette also supports execution of its extracted R code outside knitr on older R versions.

The current archive passes macOS R-devel with 0 errors, 0 warnings and the expected "New submission" NOTE. All 2,537 test expectations pass, with 22 existing skips. The extracted vignette code also passes its numerical assertions.

Windows checks are pending for this replacement archive. The preceding R-devel attempt stopped at dependency screening because the available dgraphs version was unsuitable; the preceding R-release attempt is still awaiting a result. These attempts are not checks of the current archive.

The archive is 8,730,268 bytes; documentation totals 4,956,641 bytes and data 1,763,585 bytes.
