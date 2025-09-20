# Development Commands (gflow)

A quick-access cheat sheet of common commands used during package development,
debugging, and CRAN submission prep.

---

## Quick Access (most used)
```r
devtools::load_all()       # Load package without rebuilding
devtools::document()       # Update documentation (roxygen2)
make build                 # Build tarball with Makefile
R CMD build .              # Build tarball (manual)
R CMD check gflow_*.tar.gz --as-cran   # Full CRAN check
```

## Build/package cycle:
```r
Rcpp::compileAttributes()  # Regenerate Rcpp glue
devtools::document()       # Rebuild docs/NAMESPACE
devtools::load_all()       # Reload package in current session
devtools::install()        # Install locally
devtools::clean_dll()      # Clean compiled shared objects
```

### Makefile helpers:
```bash
make build          # Rebuild tarball
make build-verbose  # Rebuild with logs
make check          # Full check
make check-fast     # Faster dev check
```

### Checks / Tests
```bash
R CMD check --as-cran gflow_*.tar.gz       # Full CRAN check
R CMD check gflow_*.tar.gz --use-valgrind  # With valgrind
```

```r
rhub::check_for_cran()      # Run CRAN checks on multiple platforms
devtools::check()           # Local R CMD check wrapper
```

## Debugging with memory tools
### In R:
```r
gctorture(TRUE)             # Force GC on every alloc (protection bugs)
gctorture2(step = 10)       # GC every 10 allocations
```

### Shell:
```bash
R -d valgrind -f script.R   # Run script under valgrind
R -d "valgrind --tool=memcheck --leak-check=full" -f script.R
```

## Rcpp / C++ Development
```r
Rcpp::compileAttributes()   # Update Rcpp exports
devtools::clean_dll()
devtools::load_all()
```

## Compilation flags
```bash
R CMD SHLIB src/file.cpp   # Compile single file
R CMD INSTALL .            # Install full package
```

## Git / GitHub
```bash
git status
git add -p
git commit -m "Message"
git push
git tag -a v0.1.0 -m "Release v0.1.0"
git push --tags
```

## CRAN Submission Prep
```r
devtools::check_win_devel()   # Test on win-builder
devtools::release()           # Semi-automated CRAN release workflow
```

Checklist:
- Update cran_docs/submission/cran_comments.md
- Ensure NEWS, README, DESCRIPTION are current
- Confirm reverse dependencies (if any)

## Package housekeeping
```bash
rm -rf gflow.Rcheck/         # Clean check dir
find src -name "*.o" -delete # Remove object files
find src -name "*.so" -delete
```

## Testing examples and vignettes
```r
devtools::run_examples()
devtools::build_vignettes()
```

## Profiling
```r
Rprof("out.prof"); run_function(); Rprof(NULL)
summaryRprof("out.prof")
```


