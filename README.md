# gflow

`gflow` is an R package for gradient-flow-based analysis of high-dimensional data.

## Collaborator Installation (`grip` + `gflow`)

The recommended path is to clone both repositories and install from source.

### 1. Clone repositories

```bash
git clone https://github.com/pgajer/grip.git
git clone https://github.com/pgajer/gflow.git
```

### 2. Install base R helpers

```bash
R -q -e 'install.packages(c("remotes","Rcpp"))'
```

### 3. Install `grip`

Run from the parent directory that contains both cloned folders:

```bash
R -q -e 'remotes::install_local("grip", dependencies=TRUE, upgrade="never")'
```

### 4. Install `gflow` (OpenMP required in `dev` profile)

```bash
R -q -e 'Sys.setenv(GFLOW_BUILD_PROFILE="dev"); remotes::install_local("gflow", dependencies=TRUE, upgrade="never")'
```

## OpenMP Requirement

`gflow` default install profile (`dev`) requires OpenMP. This is intentional for
performance-critical workflows.

Detailed OS-specific setup instructions are in `INSTALL.md`.

- Linux: GCC-based toolchains usually work out of the box.
- macOS: users must configure an OpenMP-capable toolchain (for example Homebrew GCC
  or LLVM + `libomp`).
- Windows: users need OpenMP-enabled Rtools toolchain setup.

If OpenMP is not configured, installation in `dev` profile will fail with a clear
error message.

## Verify OpenMP After Install

Run:

```r
.Call("S_gflow_openmp_diag", PACKAGE = "gflow")
```

Expected: `openmp_compiled` is `TRUE`.
