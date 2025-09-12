# BUILD.md — gflow build & install guide (R + C++)

This guide captures the exact steps to add/modify C++ (including Rcpp) and build locally.

## Command recipes

When you **add/rename/remove** any `// [[Rcpp::export]]` or change its signature:

# 1) Nuke generated files (safe: they'll be regenerated)
unlink(c("NAMESPACE", "R/RcppExports.R", "src/RcppExports.cpp"), force = TRUE)

# 2) Recreate Rcpp glue
Rcpp::compileAttributes()

# 3) Regenerate NAMESPACE + Rd from roxygen
devtools::document()

# 4) Load
devtools::load_all()


unlink(c("NAMESPACE", "R/RcppExports.R", "src/RcppExports.cpp"), force = TRUE); Rcpp::compileAttributes(); devtools::document(); devtools::load_all()



When you **only edit C++ implementation** (no export changes):
- `devtools::load_all()`
- (Windows DLL locked) `devtools::clean_dll(); devtools::load_all()`

When you **only edit R/roxygen**:
- `devtools::document()` then `devtools::load_all()`

Fresh rebuild if stale:
- Restart R → `devtools::clean_dll(); Rcpp::compileAttributes(); devtools::document(); devtools::load_all()`

## Prereqs
- R ≥ 4.3, C++17 toolchain
- DESCRIPTION:
  - `Imports: Rcpp`
  - `LinkingTo: Rcpp`
- NAMESPACE is roxygen-generated; ensure somewhere in R/:
  ```r
  #' @useDynLib gflow, .registration = TRUE
  #' @importFrom Rcpp evalCpp
  NULL
  ```

## Typical workflows
### Add a new Rcpp wrapper
```
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
```
### Edit only C++ internals
```
devtools::load_all()
```
### Move wrappers to src/rcpp/
Update `src/Makevars`:
```
CPP_SOURCES := $(wildcard *.cpp) $(wildcard rcpp/*.cpp)
OBJECTS = $(CPP_SOURCES:.cpp=.o)
CXX_STD = CXX17
```

## devtools helpers
```
devtools::load_all()
devtools::clean_dll(); devtools::load_all()
devtools::document()
devtools::install()
devtools::test()
devtools::check()
```
Build & check manually:
```
R CMD build .
R CMD check gflow_*.tar.gz --as-cran
```

## Platform notes
- Windows: DLL locks → use `clean_dll()` or restart R.
- macOS: verify SDK/BLAS flags in Makevars.
- Linux: run sanitizers outside CRAN checks.

## Pre-commit checklist
- `Rcpp::compileAttributes()` ran and `R/RcppExports.R` + `src/RcppExports.cpp` updated (if exports changed).
- `devtools::document()` updated NAMESPACE & Rd.
- `devtools::test()` green; slow tests gated on NOT_CRAN.
- `devtools::check()` clean locally.
