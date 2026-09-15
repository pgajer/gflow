# Installation with complete help and vignettes

From a cloned source checkout, run `make install-user`. This uses your current
R libraries, installs missing core and documentation-building dependencies,
regenerates Rd help, builds vignettes, installs the built archive, and checks
installed help and the introductory workflow. R, GNU make, a C++17 compiler,
and Pandoc must be available. RStudio includes Pandoc; otherwise install it
separately and put it on `PATH`.

For explicit control, the equivalent commands are:

```sh
Rscript --vanilla tools/install_build_dependencies.R
make install R_ENV="env -u R_HOME"
Rscript --vanilla tools/check_installed_guides.R
```

These commands share the active R library settings. On a machine without a
writable default library, create a user library directory, set `R_LIBS_USER`
to that directory, and run the commands in that environment. To validate a
separate installation, set `R_LIBS` to an isolated library before running them.
Requires R >= 4.1.0 and dgraphs >= 0.2.0. Published dgraphs 0.2.0 and the
0.3.0.9000 object API have separate compatibility checks for components,
hop graphs and canonical overlap cells. Future dependency changes still need
the same checks. Optional 3D viewing requires ivue >= 0.1.0; use its
[documented development installation](https://pgajer.github.io/ivue/).

`make build` produces `gflow_0.2.0.tar.gz` with generated help and all four
guides. Once the required dependencies are installed, this archive can also be
installed using `R CMD INSTALL gflow_0.2.0.tar.gz`. Direct Git/local installers
that do not regenerate Rd and build vignettes are not the complete-help route.

After installation:

```r
help("gflow-package", package = "gflow")
vignette("function-guide", package = "gflow")
browseVignettes("gflow")
```

## Optional OpenMP toolchains (macOS / Linux / Windows)

This document focuses on OpenMP toolchain setup so `gflow` can be installed in
the `dev` profile (parallel-enabled and OpenMP-required).

## Migration note (`malo`)

Legacy 1D model-averaging APIs were migrated out of `gflow` and are now
provided only by `malo` (`magelo*`, `mabilo*`, `magelog`, `fit.pwlm*`).

## Why this matters

`gflow` default build profile is `cran-safe` (portable build, OpenMP optional).
`dev` profile requires OpenMP for practical performance on key workflows.

If OpenMP is not configured, installation fails with:

`gflow dev profile requires OpenMP ...`

## macOS

`clang` from Xcode does not provide OpenMP by default.
Use one of the two options below.

### Option A (recommended): Homebrew GCC toolchain

1. Install GCC:

```bash
brew install gcc
```

2. Create/update `~/.R/Makevars` (adjust compiler version suffix if needed):

```make
CC    = /opt/homebrew/opt/gcc/bin/gcc-15
CXX   = /opt/homebrew/opt/gcc/bin/g++-15
CXX17 = /opt/homebrew/opt/gcc/bin/g++-15
CXX20 = /opt/homebrew/opt/gcc/bin/g++-15

CXXFLAGS   += -fopenmp
CXX17FLAGS += -fopenmp
LDFLAGS    += -fopenmp
```

### Option B: LLVM clang + `libomp`

1. Install LLVM and OpenMP runtime:

```bash
brew install llvm libomp
```

2. Create/update `~/.R/Makevars`:

```make
CC    = /opt/homebrew/opt/llvm/bin/clang
CXX   = /opt/homebrew/opt/llvm/bin/clang++
CXX17 = /opt/homebrew/opt/llvm/bin/clang++

CPPFLAGS   += -I/opt/homebrew/opt/libomp/include
CXXFLAGS   += -Xpreprocessor -fopenmp
CXX17FLAGS += -Xpreprocessor -fopenmp
LDFLAGS    += -L/opt/homebrew/opt/libomp/lib -lomp
```

## Linux

Most Linux GCC toolchains support OpenMP out of the box.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y build-essential gfortran
```

### Fedora / RHEL / Rocky

```bash
sudo dnf install -y gcc gcc-c++ gcc-gfortran make
```

### If OpenMP flags are missing in R

Create/update `~/.R/Makevars`:

```make
SHLIB_OPENMP_CFLAGS   = -fopenmp
SHLIB_OPENMP_CXXFLAGS = -fopenmp
SHLIB_OPENMP_FFLAGS   = -fopenmp
SHLIB_OPENMP_FCFLAGS  = -fopenmp
```

## Windows

Install Rtools (matching your R major/minor line) and ensure the Rtools UCRT
toolchain is on PATH in R sessions.

If needed, create/update `%USERPROFILE%/.R/Makevars.ucrt`:

```make
CXXFLAGS   += -fopenmp
CXX17FLAGS += -fopenmp
LDFLAGS    += -fopenmp
```

## Verify OpenMP after install

In R:

```r
.Call("S_gflow_openmp_diag", PACKAGE = "gflow")
```

Expected on an OpenMP-enabled toolchain: `openmp_compiled` is `TRUE`.

## Troubleshooting

1. If install fails with OpenMP requirement error, your compile flags do not
enable OpenMP in the active R toolchain.
2. After changing `~/.R/Makevars` or `Makevars.ucrt`, restart R before
reinstalling.
3. On macOS, double-check that R is using your configured compiler:

```bash
R CMD config CXX17
R CMD config CXX17FLAGS
```

4. To force the portable profile explicitly:

```bash
R -q -e 'Sys.setenv(GFLOW_BUILD_PROFILE="cran-safe"); system("make install-user")'
```

## Continuous platform checks

Pushes and pull requests build generated documentation and run tests on Linux,
macOS and Windows. The minimum R series is tested on Linux. Serial builds set
`GFLOW_DISABLE_OPENMP=1` with the default `cran-safe` profile; Linux and Windows
parallel builds use `GFLOW_BUILD_PROFILE=dev` and verify that OpenMP was compiled.
`GFLOW_DISABLE_OPENMP=1` and `GFLOW_BUILD_PROFILE=dev` are deliberately incompatible.
The manual R-hub workflow remains available for additional release checks.
