# Installation Notes for OpenMP (macOS / Linux / Windows)

Basic package installation steps are in `README.md`.
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
R -q -e 'Sys.setenv(GFLOW_BUILD_PROFILE="cran-safe"); remotes::install_local("gflow", dependencies=TRUE, upgrade="never")'
```
