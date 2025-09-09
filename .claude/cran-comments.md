## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new submission.

* checking installed package size ... NOTE
  installed size is 11.6Mb
  sub-directories of 1Mb or more:
    include   7.0Mb
    libs      4.0Mb
    
  The package includes bundled C++ libraries (Eigen, ANN, Spectra) necessary for
  computational geometry and linear algebra operations. These are essential
  dependencies that are not available as separate R packages.

* checking pragmas in C/C++ headers and code ... NOTE
  File which contains pragma(s) suppressing diagnostics:
    'inst/include/eigen/src/Core/util/DisableStupidWarnings.h'
    
  This is from the bundled Eigen library (version 3.4.0), a widely-used 
  third-party linear algebra library. The pragma is part of the original
  Eigen source code and is necessary for cross-platform compilation.

* Eigen library compilation warnings:
  The package may generate warnings related to memcpy operations in Eigen's
  NEON architecture optimizations. These are from the bundled third-party
  Eigen library and do not affect package functionality.

## Test environments

* local macOS install, R 4.4.0
* ubuntu 22.04 (on GitHub Actions), R 4.4.0
* win-builder (devel and release)

## Downstream dependencies

There are currently no downstream dependencies for this package.