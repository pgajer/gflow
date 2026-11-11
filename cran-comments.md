## gflow 0.2.0

Resubmission following the earlier 0.1.0 attempts. The Linux loading issue is resolved: gflow has no TBB dependency and uses published dgraphs (>= 0.4.0).

Repeated vignette fits are precomputed with full reproduction code. Scientific results, examples and figures are retained.

The rebuilt archive passes macOS R-devel with 0 errors, 0 warnings and the expected "New submission" NOTE. Win-builder R-devel stopped at dependency screening because the available dgraphs version is unsuitable. CRAN currently lists dgraphs 0.2.0 for Windows R-devel, below gflow's required 0.4.0. The fresh R-release check is still in progress.

The archive is 8,729,136 bytes; documentation totals 4,956,159 bytes and data 1,762,765 bytes. Earlier Windows passes apply to the preceding archive.
