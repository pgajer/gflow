## gflow 0.2.0

Resubmission following the earlier 0.1.0 attempts. The Linux loading issue is resolved: gflow has no TBB dependency and uses published dgraphs (>= 0.4.0).

Repeated vignette fits are precomputed with full reproduction code. Scientific results, examples and figures are retained. The vignette supports execution of its extracted R code outside knitr on older R versions.

## Checks

The exact source archive passes macOS R-devel with 0 errors, 0 warnings and the expected "New submission" NOTE. All 2,537 test expectations pass, with 22 existing skips. All five vignettes and both manuals pass; the extracted vignette script also passes its numerical assertions.

Checks built from the same prepared package sources pass on Linux R 4.1.3 and on Linux, macOS and Windows R 4.6.2. All three R-release checks report Status: OK. R 4.1.3 has no errors or warnings and three notes: unavailable optional suggested packages, installed size with debug symbols, and declared GNU make use. Installed help and the dependency-to-basin workflow pass on every platform.

Windows checks are pending for this exact replacement archive. The preceding Win-builder R-devel attempt stopped at dependency screening because the available dgraphs version was unsuitable; the preceding R-release upload is still awaiting a result. These attempts are not checks of the current archive.

The archive is 8,730,268 bytes; documentation totals 4,956,641 bytes and data 1,763,585 bytes.
