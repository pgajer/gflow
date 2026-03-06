# Repository Instructions

## Scope
- This repository is the source of truth for the `gflow` R package API and CRAN-facing behavior.

## R Style
- Prefer dot-delimited function and variable names for new R code (`fit.rdgraph.regression` style).
- Use leading-dot names only for private helpers.
- Follow existing argument names in touched files; do not rename public APIs unless explicitly requested.

## Package Hygiene
- When modifying exported functions or roxygen blocks, regenerate docs:
  - `Rscript -e 'roxygen2::roxygenise(".")'`
- Validate changes with tests before broader checks:
  - `Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'`
- Run package build/check for release-readiness work:
  - `R CMD build .`
  - `R CMD check --as-cran gflow_*.tar.gz`

## Build Artifacts
- Do not edit files under `gflow.Rcheck/` as source code.
- Keep generated artifacts (`*.tar.gz`, local check outputs) out of functional commits unless the task is release-specific.
