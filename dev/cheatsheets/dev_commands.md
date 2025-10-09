# Development Commands (gflow)

## see /Users/pgajer/Library/Logs/DiagnosticReports/
## [~/current_projects/gflow]% ls -lt /Users/pgajer/Library/Logs/DiagnosticReports/*ips | head


A quick-access cheat sheet of common commands used during package development,
debugging, and CRAN submission prep.

---

## Quick Access (most used)
```r
devtools::load_all()       # Load package without rebuilding
devtools::document()       # Update documentation (roxygen2)
roxygen2::roxygenize()
make build                 # Build tarball with Makefile
R CMD build .              # Build tarball (manual)
R CMD check gflow_*.tar.gz --as-cran   # Full CRAN check

# Clean everything
devtools::clean_dll()

# Regenerate documentation and NAMESPACE
devtools::document()

# Try loading again
devtools::load_all()
```

## Testing for non-ASCII characters within a file
```
M-x occur RET [←τŷΓ²·≈ΛℓΣβλΔ→ρμσ₀₁₂ε≥±] RET
M-x occur RET [N̂←τŷΓ²≈ΛℓΣβλΔ→ρμσ₀₁₂ε≥±] RET
```

## Interactive Debugging Setup
```r
devtools::load_all()

data <- list(
  X = matrix(rnorm(100), 50, 2),
  y = rnorm(50)
)

dcx1 <- build.nerve.from.knn(data$X, data$y, k = 10, max.p = 1)

## see /Users/pgajer/Library/Logs/DiagnosticReports/
## [~/current_projects/gflow]% ls -lt /Users/pgajer/Library/Logs/DiagnosticReports/*ips | head


source("tests/manual/setup-debug.R")
data <- make_test_data(n = 50, type = "clustered")

# Test with max.p = 1 (should work)
dcx1 <- build.nerve.from.knn(data$X, data$y, k = 10, max.p = 1)
riem.dcx.summary(dcx1)

# Test with max.p = 2 (triangles)
dcx2 <- build.nerve.from.knn(data$X, data$y, k = 10, max.p = 2)
riem.dcx.summary(dcx2)

check_all_duplicates(dcx2)

# Only if max.p = 2 works, try max.p = 3
dcx3 <- build.nerve.from.knn(data$X, data$y, k = 10, max.p = 3)
riem.dcx.summary(dcx3)
check_all_duplicates(dcx3)

# Try k=20 (40% of points)
dcx_k20 <- build.nerve.from.knn(data$X, data$y, k = 20, max.p = 3)
riem.dcx.summary(dcx_k20)

# Directed k-NN is less restrictive
dcx.dir <- build.nerve.from.knn(data$X, data$y, k = 20, max.p = 3, directed.knn = TRUE)
riem.dcx.summary(dcx.dir)




# Start debugging session
source("tests/manual/setup-debug.R")

# Quick test with defaults
dcx <- quick_check()

# Test with clustered data
data <- make_test_data(n = 100, type = "clustered")
dcx <- quick_check(data$X, data$y, k = 15, max.p = 3)

# Detailed inspection
print_simplices(dcx, 2)  # Show triangles
check_all_duplicates(dcx, verbose = TRUE)

# Visualize (if you have ggplot2)
plot_nerve_2d(dcx, data$X, show_triangles = TRUE)

# After modifying C++ code
reload_cpp()
dcx <- quick_check()

# Test different scenarios
for (type in c("random", "clustered", "circle")) {
  cat("\n\n### Testing", type, "data ###\n")
  data <- make_test_data(type = type)
  dcx <- quick_check(data$X, data$y, k = 12, max.p = 3)
}
```

## Running Manual Tests
```r
# From package root
source("tests/manual/debug-nerve-duplicates.R")

# Or load the package first
devtools::load_all()
source("tests/manual/debug-nerve-duplicates.R")

# Or interactively step through
file.edit("tests/manual/debug-nerve-duplicates.R")
# Then run chunks with Ctrl+Enter or Cmd+Enter
```

## Run all tests
```r
devtools::test()
```

## Run specific test file
```r
devtools::test(filter = "nerve-construction")
```

## Run with detailed output
```r
devtools::test(reporter = "progress")
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


