# gflow Split Audit SA4: Self-Contained geosmooth Native/Support Contract

Generated: 2026-06-03

Branch: `codex/gflow-split-audit-sa4`

## Decision Frozen in SA4

`geosmooth` should be a self-contained smoothing package for the first physical
split.  In particular:

> `geosmooth` will vendor ANN and own its ordinary local-neighborhood support
> computations, instead of depending on unexported `gflow` kNN or ANN-native
> internals.

This changes the extraction strategy from "import all local-neighborhood
services from `gflow`" to "duplicate the small amount of native infrastructure
needed for the smoother package to stand on its own."

No source files are moved in SA4.  This is a contract for the next physical
split phases.

## Target First geosmooth Payload

The first `geosmooth` payload should still move the coherent smoother cluster
together:

- LPS, currently `kernel.local.polynomial.cv`;
- MALPS;
- LPL-TF;
- SLPLiFT / S-LPL-TF;
- SSRHE.

Public function names should remain unchanged during the split.  The cleaner
post-split API name `kernel.local.polynomial` or `fit.lps` should be added only
after the package is installable and passing tests.

## What geosmooth Will Own

### Coordinate kNN and Local Supports

`geosmooth` will own local support construction for the smoothing methods:

- coordinate kNN supports;
- fixed-radius coordinate supports, if retained;
- adaptive-radius supports needed by LPS, MALPS, LPL-TF, SLPLiFT, and SSRHE;
- support metadata needed by CV, refitting, diagnostics, and reproducibility.

This requires vendoring ANN.

### ANN Vendor Payload

Copy from `gflow`:

- `src/ANN/*`
- `inst/licenses/ANN-Copyright-Notice.txt`

Current ANN source count in `gflow`: 21 top-level files under `src/ANN`.

The `geosmooth` source package must preserve the ANN license notice.  Because
ANN is LGPL-2.1-or-later and `gflow` is GPL (>= 3), this is compatible with a
GPL `geosmooth` package.

### Local PCA Charts

`geosmooth` should own local PCA chart construction for the first physical
split.  This avoids cross-package calls into unexported native symbols and
keeps LPS/LPL-TF/SLPLiFT/SSRHE behavior identical after extraction.

Copy from `gflow`:

- `src/local_pca_charts.cpp`
- `src/local_pca_charts.hpp`
- `src/local_pca_charts_rcpp.cpp`
- `R/local_pca_chart_dim.R`
- `tests/testthat/test-local-pca-charts.R`

Rationale: local PCA is geometric, but it is also a core smoother design
primitive.  Keeping it in `geosmooth` initially reduces transition risk.  A
later refactor can move the shared implementation back to a lower-level
geometry package if duplication becomes painful.

### Kernel/LPS Native Backend

Copy from `gflow`:

- `src/kernel_local_polynomial_cv_rcpp.cpp`
- `R/kernel_local_polynomial_cv.R`
- `tests/testthat/test-kernel-local-polynomial-cv.R`

Regenerate, do not copy by hand:

- `R/RcppExports.R`
- `src/RcppExports.cpp`
- `src/init.c` entries for Rcpp wrappers.

The C++ backend uses:

- Rcpp;
- ANN;
- LAPACK through `R_ext/Lapack.h`.

### SSRHE Native Backend

If SSRHE moves in the first cluster, copy:

- `R/ssrhe_hessian_energy.R`
- `src/ssrhe_hessian_energy.cpp`
- `src/ssrhe_hessian_energy_r.h`
- `tests/testthat/test-ssrhe-hessian-energy.R`

SSRHE uses Eigen and the shared local PCA helper.  If local PCA is copied into
`geosmooth`, SSRHE can remain native and self-contained.

### MALPS, LPL-TF, and SLPLiFT R Files

Copy:

- `R/malps.R`
- `R/lpl_tf.R`
- `R/slpl_tf.R`
- `tests/testthat/test-malps.R`
- `tests/testthat/test-lpl-tf.R`
- `tests/testthat/test-slpl-tf.R`

Because MALPS moves with LPS, `.malps.design.matrix` does not need to be
factored before the initial split.  It may be renamed later to a neutral helper
such as `.local.polynomial.design.matrix`.

## What geosmooth Should Not Initially Own

Do not copy the full graph/gradient-flow stack into `geosmooth`:

- gradient-flow cells;
- Morse-Smale complexes;
- basin machinery;
- graph visualization workflows;
- PHATE/diffusion graph infrastructure;
- microbiome/application-specific workflows.

Do not copy Qhull for the first `geosmooth` payload unless a moved method
proves it needs Qhull.  The LPS/MALPS/LPL/SLPL/SSRHE cluster does not require
the full qhull source tree for ordinary coordinate-support smoothing.

## Graph-Geodesic Support Policy

Some existing smoother code has graph-geodesic support modes.  These are not
the same as ordinary coordinate/adaptive-radius local supports.

For the first physical split:

1. Coordinate kNN/adaptive-radius supports should be self-contained in
   `geosmooth`.
2. Graph-geodesic support modes should be handled explicitly by one of these
   policies:
   - temporarily require/import `gflow` for graph-geodesic support;
   - mark graph-geodesic support experimental/deferred in `geosmooth`;
   - port a minimal graph-distance service only after the core package builds.
3. The package must not silently change support semantics during extraction.

Recommended first choice: keep graph-geodesic modes available only when `gflow`
is installed, with a clear error otherwise.  This keeps `geosmooth` practical
without forcing the entire graph engine into the first split.

## Build Contract

### DESCRIPTION

Initial `geosmooth` should use:

- `Package: geosmooth`
- `Title: Geometric Smoothing and Conditional Expectation Methods`
- `License: GPL (>= 3)`
- `SystemRequirements: C++17, GNU make`
- `LinkingTo: Rcpp`
- `Imports: Rcpp, stats, methods, graphics, grDevices, Matrix`
- `Suggests: testthat, genlasso, igraph, gflow`

Notes:

- `genlasso` remains a solver dependency for LPL-TF/SLPLiFT paths.
- `gflow` can be in `Suggests` initially for graph-geodesic modes and migration
  tests, not for ordinary coordinate support.
- If we choose to use `RcppEigen` instead of vendored Eigen later,
  `LinkingTo` should include `RcppEigen`.  For the lowest-friction first split,
  copy the existing vendored Eigen headers and `eigen_config.hpp`.

### Native Headers

Copy from `gflow` for the first split:

- `inst/include/Eigen/*`
- `inst/include/gflow/eigen_config.hpp`

Rename `inst/include/gflow/eigen_config.hpp` to a package-neutral or
package-specific path later, for example:

- `inst/include/geosmooth/eigen_config.hpp`

During initial extraction, keeping the file name temporarily is acceptable if
it minimizes compile churn, but the package-facing include path should not
permanently say `gflow`.

### src/Makevars

Initial Unix Makevars should include the same core ingredients as `gflow`, but
without qhull unless needed:

```make
PKG_CPPFLAGS = -I. -IANN
PKG_CPPFLAGS += -I../inst/include -DR_NO_REMAP
PKG_CPPFLAGS += -D_Alignof=__alignof__
PKG_CPPFLAGS += -DDLL_API=
PKG_CPPFLAGS += -include ../inst/include/geosmooth/eigen_config.hpp

CXX_STD = CXX17
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

SOURCES_CPP = $(wildcard *.cpp)
SOURCES_C = $(wildcard *.c)
ANN_SOURCES = $(wildcard ANN/*.cpp)
OBJECTS = $(SOURCES_CPP:.cpp=.o) $(SOURCES_C:.c=.o) $(ANN_SOURCES:.cpp=.o)

$(SHLIB): $(OBJECTS)
```

OpenMP should be added only if a moved native file actually needs it.  Avoid
carrying the full `gflow` build profile system into `geosmooth` unless required.

### src/Makevars.win

Initial Windows Makevars should mirror the ANN and Eigen settings:

```make
PKG_CPPFLAGS = -I. -IANN
PKG_CPPFLAGS += -I../inst/include -DR_NO_REMAP
PKG_CPPFLAGS += -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_PARALLELIZE
PKG_CPPFLAGS += -D_Alignof=__alignof__
PKG_CPPFLAGS += -DDLL_EXPORTS
PKG_CPPFLAGS += -DGFLOW_NO_TBB
PKG_CPPFLAGS += -include ../inst/include/geosmooth/eigen_config.hpp

CXX_STD = CXX17
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

SOURCES_CPP = $(wildcard *.cpp)
SOURCES_C = $(wildcard *.c)
ANN_SOURCES = $(wildcard ANN/*.cpp)
OBJECTS = $(SOURCES_CPP:.cpp=.o) $(SOURCES_C:.c=.o) $(ANN_SOURCES:.cpp=.o)

$(SHLIB): $(OBJECTS)
```

The `GFLOW_NO_TBB` macro name should eventually be renamed, but it may remain
temporarily if included Eigen/Spectra config expects it.

## Registration Contract

Regenerate Rcpp exports in `geosmooth`; do not copy generated files unchanged.

Expected Rcpp-generated wrappers:

- `rcpp_kernel_local_polynomial_cv_coordinates`
- `rcpp_kernel_local_polynomial_predict_coordinates`
- `rcpp_local_pca_chart`

Expected manual `.Call` registrations if SSRHE moves:

- `S_ssrhe_hessian_operator`

The generated native symbols should use the `geosmooth` namespace prefix, for
example:

- `_geosmooth_rcpp_kernel_local_polynomial_cv_coordinates`
- `_geosmooth_rcpp_kernel_local_polynomial_predict_coordinates`
- `_geosmooth_rcpp_local_pca_chart`

Do not keep `_gflow_*` symbols in `geosmooth`.

## Test Contract

Initial tests to move:

- `tests/testthat/test-kernel-local-polynomial-cv.R`
- `tests/testthat/test-local-pca-charts.R`
- `tests/testthat/test-malps.R`
- `tests/testthat/test-lpl-tf.R`
- `tests/testthat/test-slpl-tf.R`
- `tests/testthat/test-ssrhe-hessian-energy.R`

Initial test expectations:

- coordinate backend and R backend agree on small examples;
- local PCA chart output remains unchanged relative to current `gflow`;
- LPL-TF and SLPLiFT local-PCA paths produce the same smoke-test results;
- SSRHE operator tests pass after native registration is regenerated;
- graph-geodesic support either passes with `gflow` installed or gives a clear
  skip/error if intentionally deferred.

## First Physical Split Phases

### GE0: Create geosmooth Skeleton

Create package structure, DESCRIPTION, NAMESPACE scaffold, `R/`, `src/`,
`tests/testthat/`, and license files.  Add vendored ANN and Eigen/config
headers.  Do not move public methods yet.

### GE1: Move R-Level Smoother Cluster

Move/copy the five target R files and minimal helper R files.  Keep public
names unchanged.  Temporarily disable native backends only if needed to get the
namespace loading.

### GE2: Add Native LPS and Local PCA

Add the LPS C++ backend and local-PCA C++ backend.  Regenerate Rcpp exports and
run focused LPS/local-PCA tests.

### GE3: Add LPL-TF, SLPLiFT, and MALPS Tests

Bring up local-PCA dependent LPL-TF/SLPLiFT paths and MALPS tests.  Fix only
split-induced namespace or helper issues.

### GE4: Add SSRHE Native Backend

Move/register the SSRHE Hessian operator and run SSRHE tests.  If SSRHE creates
too much compile or registration friction, defer it behind a clear GE4 gate
rather than blocking GE1-G3.

### GE5: Compatibility Layer in gflow

After `geosmooth` passes focused checks, add `gflow` compatibility wrappers or
documentation pointing users to `geosmooth`.  Do not deprecate aggressively
until downstream project scripts have migrated.

## SA4 Verdict

Vendoring ANN into `geosmooth` is the right practical choice for the transition.
It makes the first smoother package self-contained and avoids designing a
cross-package native kNN API before we know how the split feels in practice.

The cost is duplication of ANN and local-PCA infrastructure.  That duplication
is acceptable for the first extraction because it reduces risk and keeps
ordinary local-neighborhood smoothing independent of `gflow`.

The next step is **GE0: create the `geosmooth` package skeleton** on a new
implementation branch, using this SA4 contract as the checklist.
