## Test environments
- macOS 15.6.1 (local), R Under development (2025-08-22 r88678), aarch64-apple-darwin20, GCC, C++17

## R CMD check results
0 errors | 0 warnings | 2 notes

* New submission.

* NOTE on pragmas in C/C++ headers:
  Reported file is 'inst/include/Eigen/src/Core/util/DisableStupidWarnings.h',
  which is part of the upstream Eigen library. We do not add any additional
  diagnostic-suppressing pragmas in our code. The header is included
  unmodified under Eigen’s MPL-2.0 license solely to enable our linear
  algebra paths.

* Third-party code:
  - Eigen and Spectra (MPL-2.0)
  - ANN (LGPL-2.1-or-later)
  License texts are in inst/licenses; attributions in inst/COPYRIGHTS.
  
* Parallelization:
  OpenMP usage is optional; the package builds and runs without OpenMP
  when it is unavailable.
