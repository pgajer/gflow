## Test environments
* local macOS (Apple Silicon), R 4.x and R-devel
* win-builder: devel, release, oldrelease — 0 ERRORs / 0 WARNINGs / 2 NOTEs
* Linux (local): R CMD check --as-cran — clean
* R-hub:
  - ubuntu-release — clean
  - rchk — clean
  - valgrind — no definite leaks reported
  - windows / macos-arm64 / clang-asan / gcc-asan — initiated; some runs encountered CI runner issues unrelated to package code

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE: “New submission.” — This is the first CRAN submission of ‘gflow’.
* NOTE: “Possibly misspelled words in DESCRIPTION: Smale, gflow.”
  - “Smale” is a proper noun (mathematician Stephen Smale).
  - “gflow” is the package name.

* NOTE: “pragmas in C/C++ headers … DisableStupidWarnings.h.”
  - This header is from the upstream Eigen distribution and only disables noisy compiler diagnostics; it does not introduce non-portable behavior.

## Additional details
* C++17 is required (declared in DESCRIPTION and Makevars).
* OpenMP is optional; the package builds and runs without OpenMP where unsupported.
* No non-portable compiler flags (e.g., -march=native) are used.
* The package writes only to temp directories during examples/tests; no network access is required.
* Examples and tests complete in ~20 seconds locally.
