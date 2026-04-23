# Repository Instructions

## Scope
- This repository is the source of truth for the `gflow` R package API and CRAN-facing behavior.

## Preferred Skills
- Prefer `$r-package-qa` for package QA, CRAN-style checks, documentation drift, and release-readiness work.

## R Style
- Prefer dot-delimited function and variable names for new R code (`fit.rdgraph.regression` style).
- Use leading-dot names only for private helpers.
- Follow existing argument names in touched files; do not rename public APIs unless explicitly requested.

## Package Hygiene
- When modifying exported functions or roxygen blocks, regenerate docs:
  - `make document`
- Validate changes with tests before broader checks:
  - `Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'`
- Run package QA via Makefile targets (preferred; ensures roxygen/doc generation runs first):
  - Fast QA: `make check-fast`
  - Full CRAN-style QA: `make check`
- Do not run `R CMD build` / `R CMD check` directly unless explicitly requested.

## Build Artifacts
- Do not edit files under `gflow.Rcheck/` as source code.
- Keep generated artifacts (`*.tar.gz`, local check outputs) out of functional commits unless the task is release-specific.
