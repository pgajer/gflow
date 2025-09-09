# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`gflow` is an R package implementing gradient flow-based data analysis methods. It provides tools for adaptive k-nearest neighbor graph construction, graph-based smoothing, Morse-Smale statistical gradient-flow based data domain decomposition, and geometric data analysis on high-dimensional data.

## CRAN Compliance Requirements
- All examples must complete in < 5 seconds (use \donttest{} for longer ones)
- No side effects outside of tempdir()
- Proper cleanup with on.exit() where needed
- All exported functions must be documented
- Package must pass with --as-cran flag

## Package Location
- Primary location: `~/current_projects/gflow`
- Parent package: `~/current_projects/msr2`

### Parent Package
This is an R package built based on the msr2 R package (location: ~/current_projects/msr2)
gflow contains all C/C++ functions of msr2, but only a subset of R functions in order to fix all CRAN issues on a smaller function base. Any missing R functions/files issues can be fixed by copying the corresponding files from ~/current_projects/msr2/R 

## Build and Development Commands

### Building the Package
```bash
cd ~/current_projects/gflow
# Clean previous builds
find src -name "*.o" -delete
find src -name "*.so" -delete

# Build the package
cd ~/current_projects
R CMD build gflow
```

### Checking the Package
```bash
# Full CRAN check
R CMD check gflow_0.1.0.tar.gz --as-cran

# Quick check without examples
R CMD check gflow_0.1.0.tar.gz --no-examples
```

### Installing the Package
```bash
R CMD INSTALL gflow_0.1.0.tar.gz
```

## Architecture and Core Components

### Graph Construction Framework
The package centers around adaptive k-nearest neighbor (ikNN) graph construction:
- **Main Implementation**: `src/iknn_graphs.cpp` - Core C++ implementation using ANN library for efficient nearest neighbor searches
- **R Interface**: `R/iknn_graphs.R` - R wrapper functions (`create.iknn.graphs`, `create.single.iknn.graph`)
- **Graph Types**: Implements intersection-weighted kNN graphs where edges exist between points sharing neighbors

### Key Algorithmic Components

1. **Graph-Based Methods** (`src/graph_*.cpp`):
   - Spectral analysis and filtering
   - Laplacian eigenvectors computation
   - Diffusion-based smoothing
   - Shortest path algorithms
   - Connected components analysis

2. **Smoothing Algorithms**:
   - LOWESS-based graph smoothing (`deg0_lowess_graph_smoothing.cpp`)
   - Spectral LOWESS (`spectral_lowess_graph_smoothing.cpp`)
   - Harmonic smoothing (`harmonic_smoother.cpp`)
   - Mean shift smoothing (`mean_shift_smoother.cpp`)

3. **Piecewise Linear Modeling** (`R/pwlm.R`):
   - Automatic breakpoint detection
   - Optimal segmentation algorithms
   - Model fitting and prediction

4. **Gradient Flow Computation** (`src/gflow_*.cpp`):
   - Basin detection and analysis
   - Complex-based computations
   - Monotonic reachability

### External Dependencies

The package bundles three major C++ libraries in `inst/include/`:
- **ANN**: Approximate Nearest Neighbor library for efficient kNN searches
- **Eigen**: Linear algebra library (requires version >= 3.4.0)
- **Spectra**: Sparse eigenvalue computation library

### R/C++ Interface Pattern

The package follows a consistent pattern for R/C++ integration:
1. C++ computational core in `src/*.cpp` files
2. R wrapper functions in `R/*.R` files
3. Registration via `src/init.c` for `.Call()` interfaces
4. Headers in `src/*.h` or `src/*.hpp` for shared definitions

## Development Notes

### Header File Placement
When developing R packages with C/C++ code:
- **Keep internal header files in `src/`** - Header files (`.h`, `.hpp`) used by your package's C/C++ source code should remain in the `src/` directory alongside your `.c` and `.cpp` files. This is the standard location.
- **Use `inst/include/` only for exported headers** - Only place header files in `inst/include/` if you're creating a package that provides a C/C++ API for other packages to use. These headers will be installed and made available to other packages.
- **Ignore cosmetic warnings** - R CMD check may produce warnings about header files in `src/`, but these are cosmetic and won't prevent CRAN submission. Do not move headers from `src/` to `inst/` to address these warnings.

### OpenMP Configuration
The package uses OpenMP for parallelization. Configuration is in `src/Makevars`:
```makefile
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

### R CMD Check Compliance
The codebase has been modified to comply with CRAN requirements:
- Replaced `fprintf(stderr, ...)` with R's `error()` function
- Replaced `rand()` with R's `unif_rand()` with proper RNG state management
- Removed direct `exit()` calls
- Fixed `.Call()` interfaces to use quoted strings

### Function Naming Convention
- **Standard functions**: Use dot-separated names (e.g., `create.iknn.graphs`)
- **S3 methods**: Use format `method.class_name` where class names may use underscores
  - Example: `plot.gaussian_mixture` for class "gaussian_mixture"
- **Internal functions**: Prefix with dot (e.g., `.internal.function`)

## Testing Commands
### Quick iteration
make check-fast

### Full CRAN check
make check
# or manually:
./check_cran.sh

### Platform-specific checks
R -e "rhub::check_for_cran()"
R -e "devtools::check_win_devel()"


## Handling Mathematical Formulas and Non-ASCII Symbols

### Overview
When documenting R functions that contain mathematical formulas or non-ASCII symbols, special care must be taken to ensure compatibility with R's documentation system (Rd format) and CRAN requirements.

### Key Rules

1. **Never use LaTeX dollar signs (`$...$` or `$$...$$`) in roxygen2 comments**
   - These will cause "unknown macro" errors during package checking
   - Use `\eqn{}` for inline math
   - Use `\deqn{}` for display (centered) math

2. **Avoid non-ASCII characters in comments and documentation**
   - Replace `≈` with `\approx` inside math environments
   - Replace `≤` with `\leq` or `<=`
   - Replace `≥` with `\geq` or `>=`
   - Replace `∞` with `\infty` inside math environments

### Conversion Guide

#### Inline Math
```r
# WRONG
#' The mean parameter $\mu$ controls the center

# CORRECT
#' The mean parameter \eqn{\mu} controls the center
```

#### Display Math
```r
# WRONG
#' $$f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$

# CORRECT
#' \deqn{f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}}
```

#### Complex Formulas
```r
# WRONG
#' $f'(t_i) ≈ \frac{-f(t_{i+2}) + 8f(t_{i+1}) - 8f(t_{i-1}) + f(t_{i-2})}{12\Delta t}$

# CORRECT
#' \deqn{f'(t_i) \approx \frac{-f(t_{i+2}) + 8f(t_{i+1}) - 8f(t_{i-1}) + f(t_{i-2})}{12\Delta t}}
```

### Common Mathematical Symbols

| Symbol | Wrong | Correct in Math | Correct in Text |
|--------|-------|-----------------|-----------------|
| Approximately | ≈ | `\approx` | `approximately` |
| Less than or equal | ≤ | `\leq` | `<=` |
| Greater than or equal | ≥ | `\geq` | `>=` |
| Not equal | ≠ | `\neq` | `!=` |
| Infinity | ∞ | `\infty` | `Inf` |
| Plus/minus | ± | `\pm` | `+/-` |
| Multiplication | × | `\times` | `*` |
| Summation | Σ | `\sum` | `sum` |
| Product | ∏ | `\prod` | `product` |

### Finding and Fixing Issues

#### Find all LaTeX math expressions
```bash
# Find inline math with $...$
grep -n -E '\$[^$]+\$' R/*.R

# Find specific non-ASCII characters
grep -n -E '[≈≤≥±∞×÷≠]' R/*.R

# Find all non-ASCII characters
grep -n -P '[^\x00-\x7F]' R/*.R
```

#### Automated checking before submission
```r
# Add to your pre-CRAN checklist
check_math_in_docs <- function() {
  files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  issues <- list()
  
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    
    # Check for $ math delimiters
    math_lines <- grep("\\$[^$]+\\$", lines)
    if (length(math_lines) > 0) {
      issues[[f]] <- c(issues[[f]], paste("Line", math_lines, ": Contains $...$ math"))
    }
    
    # Check for non-ASCII
    non_ascii <- grep("[^\x01-\x7F]", lines, perl = TRUE)
    if (length(non_ascii) > 0) {
      issues[[f]] <- c(issues[[f]], paste("Line", non_ascii, ": Contains non-ASCII characters"))
    }
  }
  
  if (length(issues) > 0) {
    cat("Documentation issues found:\n")
    for (f in names(issues)) {
      cat("\n", f, ":\n", sep = "")
      cat(paste("  ", issues[[f]], collapse = "\n"), "\n")
    }
    return(FALSE)
  }
  
  cat("No documentation issues found.\n")
  return(TRUE)
}
```

### Best Practices

1. **Always test documentation rendering**
   ```r
   devtools::document()
   devtools::check()
   ```

2. **Use proper roxygen2 structure for math-heavy functions**
   ```r
   #' @title Statistical Test Function
   #' @description 
   #' Performs a test with statistic:
   #' \deqn{T = \frac{\bar{X} - \mu_0}{s/\sqrt{n}}}
   #' 
   #' @param x Numeric vector of observations
   #' @param mu0 Null hypothesis mean \eqn{\mu_0}
   #' 
   #' @details
   #' Under the null hypothesis, the test statistic follows a 
   #' t-distribution with \eqn{n-1} degrees of freedom.
   #'
   #' @return 
   #' A list with components:
   #' \itemize{
   #'   \item \code{statistic}: The t-statistic
   #'   \item \code{p.value}: The p-value
   #'   \item \code{df}: Degrees of freedom
   #' }
   ```

3. **For ASCII art or formatted output in examples**
   ```r
   #' @examples
   #' \dontrun{
   #' # Result will display as:
   #' #   Estimate   Std.Error   t-value   p-value
   #' #   1.234      0.567       2.177     0.031
   #' }
   ```

4. **When math is in code comments (not documentation)**
   ```r
   # Inside function: Unicode is OK in comments
   # f'(t_i) ≈ (-f(t_{i+2}) + 8*f(t_{i+1}) - 8*f(t_{i-1}) + f(t_{i-2})) / (12*dt)
   
   # But for CRAN, consider using ASCII
   # f'(t_i) ~= (-f(t_{i+2}) + 8*f(t_{i+1}) - 8*f(t_{i-1}) + f(t_{i-2})) / (12*dt)
   ```

### Error Messages and Solutions

| Error | Solution |
|-------|----------|
| `unknown macro '\frac'` | Wrap math in `\eqn{}` or `\deqn{}` |
| `Unicode character ≈ (U+2248)` | Replace with `\approx` in math or `approximately` in text |
| `unexpected INCOMPLETE_STRING` | Check for unmatched quotes in math expressions |
| `unrecognized escape in character string` | Double-check backslashes in LaTeX commands |

### Quick Reference Card

```r
# Documentation header for math-heavy functions
#' @title Function with Mathematical Formulas
#' @description Brief description with inline math \eqn{\alpha = 0.05}
#' 
#' @param x Input vector
#' @param mu Mean parameter \eqn{\mu}
#' @param sigma Standard deviation \eqn{\sigma > 0}
#' 
#' @details
#' This function implements the formula:
#' \deqn{f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{x} is the input value
#'   \item \eqn{\mu} is the mean
#'   \item \eqn{\sigma} is the standard deviation
#' }
#' 
#' @return Numeric vector of densities
#' 
#' @examples
#' x <- seq(-3, 3, 0.1)
#' y <- my_function(x, mu = 0, sigma = 1)
#' plot(x, y, type = "l")
#' 
#' @export
```

### Important: Correct \itemize Syntax

Never use `\item{name}{description}` - this will cause "Lost braces" warnings.

**Incorrect:**
```r
\item{variable}{Description of variable}
```

**Correct:**
```r
\item \code{variable}: Description of variable
```

Here's a comprehensive CLAUDE.md instruction file for handling R documentation bracket issues:


# R Documentation Bracket Issues

## Overview
This section provides instructions for identifying and fixing bracket notation issues in R package documentation that cause "Missing link(s)" warnings during R CMD check.

## The Problem
When writing roxygen2 documentation, bracket notation like `[i,j]`, `[0,2pi]`, or any `[,]` patterns are interpreted as cross-references by Rd, causing warnings like:
```
Missing link(s) in Rd file 'function.Rd':
  'i,j'
```

## Quick Fix Guide

### For any bracket notation in roxygen2 comments:
1. **Use `\code{}`**: `\code{matrix[i,j]}` (preferred)
2. **Escape brackets**: `matrix\[i,j\]`

### Example Fix:
```r
# WRONG:
#' @return Matrix where element [i,j] represents...

# CORRECT:
#' @return Matrix where element \code{[i,j]} represents...
```

## Detection Instructions

### 1. Find All Potential Issues
Run this R script to scan all R files for potential bracket issues:

```r
# scan_brackets.R
files <- list.files(pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
issues_found <- FALSE

for (file in files) {
  lines <- readLines(file, warn = FALSE)
  line_numbers <- seq_along(lines)
  
  # Find roxygen lines with potential bracket issues
  roxygen_lines <- grep("^#'", lines)
  
  for (i in roxygen_lines) {
    line <- lines[i]
    # Check for unescaped brackets with comma
    if (grepl("[^\\\\]\\[.*,.*\\]", line) || grepl("^\\[.*,.*\\]", line)) {
      if (!issues_found) {
        cat("=== R Documentation Bracket Issues Found ===\n\n")
        issues_found <- TRUE
      }
      cat(sprintf("File: %s (line %d)\n", file, i))
      cat(sprintf("  %s\n", line))
      cat(sprintf("  Fix: Use \\code{} around bracket notation\n\n"))
    }
  }
}

if (!issues_found) {
  cat("No bracket issues found in R documentation.\n")
}
```

### 2. Common Patterns to Check
Look for these patterns in roxygen2 comments (#'):
- `[i,j]` - Matrix/array indices
- `[0,1]` - Numeric ranges
- `[a,b]` - Variable ranges
- `[0,2pi]` - Mathematical ranges
- `[-inf,inf]` - Infinite ranges
- `[start,end]` - Named ranges

### 3. Command Line Search
```bash
# Find all potential bracket issues in R files
grep -n "#'.*\[.*,.*\]" *.R

# Find specific patterns
grep -n "#'.*\[[0-9],.*\]" *.R  # Numeric indices
grep -n "#'.*\[[a-zA-Z],.*\]" *.R  # Variable indices
```

## Automated Fixes

### Replace All Instances (with caution):
```r
# fix_brackets.R
fix_roxygen_brackets <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  modified <- FALSE
  
  for (i in seq_along(lines)) {
    if (grepl("^#'", lines[i])) {
      # Replace [x,y] with \code{[x,y]} in roxygen lines
      new_line <- gsub("([^\\\\]|^)(\\[[^]]*,[^]]*\\])", "\\1\\\\code{\\2}", lines[i])
      if (new_line != lines[i]) {
        lines[i] <- new_line
        modified <- TRUE
      }
    }
  }
  
  if (modified) {
    writeLines(lines, file_path)
    cat("Fixed:", file_path, "\n")
  }
  
  return(modified)
}

# Apply to all R files
files <- list.files(pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
fixed_count <- sum(sapply(files, fix_roxygen_brackets))
cat("\nFixed", fixed_count, "files\n")
```

## Prevention Guidelines

### 1. Documentation Templates
Use these templates for common cases:

```r
#' @param matrix A matrix where \code{matrix[i,j]} represents the element at row i, column j
#' @param range A numeric vector in the interval \code{[0,1]}
#' @return A list where \code{result[[i]]} contains the i-th element
#' @details The algorithm works on the domain \code{[a,b]} where a < b
```

### 2. Pre-commit Check
Add this to your development workflow:

```r
# check_docs.R
check_bracket_issues <- function() {
  cmd <- "R CMD check --no-manual --no-tests --no-vignettes --no-build-vignettes ."
  result <- system(cmd, intern = TRUE)
  
  # Look for missing link warnings
  warnings <- grep("Missing link.*in Rd file", result, value = TRUE)
  
  if (length(warnings) > 0) {
    cat("Documentation warnings found:\n")
    cat(warnings, sep = "\n")
    
    # Check if any are bracket-related
    bracket_warnings <- grep("'[^']*,'", warnings, value = TRUE)
    if (length(bracket_warnings) > 0) {
      cat("\nLikely bracket notation issues detected!\n")
      cat("Run scan_brackets.R to find and fix them.\n")
    }
  }
}
```

### 3. Style Guide Rules
- Always wrap array/matrix notation in `\code{}`
- Use `\code{[lower, upper]}` for intervals
- Use `\code{object[[i]]}` for list extraction
- Use `\code{array[i,j,k]}` for multi-dimensional arrays

## Quick Reference Card

| Pattern | Wrong | Correct |
|---------|-------|---------|
| Matrix index | `[i,j]` | `\code{[i,j]}` |
| Range | `[0,1]` | `\code{[0,1]}` |
| List element | `[[i]]` | `\code{[[i]]}` |
| Array | `[i,j,k]` | `\code{[i,j,k]}` |
| Math interval | `[a,b)` | `\code{[a,b)}` |

## Testing Your Fixes

After making changes:
```bash
# Quick check for documentation issues only
R CMD build . && R CMD check --as-cran --no-tests *.tar.gz

# Or in R:
devtools::check_man()  # Check documentation only
devtools::document()   # Regenerate documentation
```

## Notes
- This issue only affects roxygen2 comments (lines starting with #')
- Regular code comments (#) are not affected
- The issue occurs because Rd format interprets [,] as cross-references
- Using `\code{}` is preferred over escaping as it also provides formatting


# Pre-Submission Checklist

- [ ] Run `grep -n -E '\$[^$]+\$' R/*.R` - should return nothing
- [ ] Run `grep -n -P '[^\x00-\x7F]' R/*.R` - check all matches
- [ ] Run `devtools::document()` - should complete without warnings
- [ ] Run `devtools::check()` - should have no Rd warnings
- [ ] Build and view PDF manual to verify math rendering
- [ ] Test on Windows if using any special characters

# Additional Resources

- [Writing R Extensions - Mathematics](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Mathematics)
- [Roxygen2 Documentation](https://roxygen2.r-lib.org/)
- [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html)


