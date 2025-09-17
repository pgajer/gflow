
# rchk function correction instructions — Pragmatic (LENGTH‑first) Policy

You are an expert C/C++ engineer familiar with R’s C API, rchk, and writing
robust SEXP wrappers. This guideline is a **practical** variant for typical
statistical workloads where objects never(!) exceed 2^31−1 elements (especially
in 100+ dimensional matrices). It aims to keep code simple and rchk‑clean
without over‑engineering for long vectors that are infeasible in memory for most
use cases.

## Goals

- Produce **rchk‑clean** wrappers with simple, readable protection patterns.
- Prefer **LENGTH() + int / R_len_t** for sizes; use `size_t` for STL counts.
- Avoid long‑vector plumbing unless explicitly needed.

> Note: If a project later needs long‑vector support, enable it with a single macro switch
> (see `GFLOW_ENABLE_LONG_VECTORS` below) without rewriting call sites.

---

## 0) Optional long‑vector switch

Add this in a shared header (e.g., `rchk_helpers.hpp`):

```c
#ifndef GFLOW_ENABLE_LONG_VECTORS
#  define GFLOW_ENABLE_LONG_VECTORS 0
#endif

#if GFLOW_ENABLE_LONG_VECTORS
#  define RLEN_T  R_xlen_t
#  define RLEN(x) XLENGTH(x)
#else
#  define RLEN_T  int /* or R_len_t */
#  define RLEN(x) LENGTH(x)
#endif
```

- Default: `GFLOW_ENABLE_LONG_VECTORS=0` (pragmatic, LENGTH‑first).
- Future‑proof: flip to `1` to opt into `XLENGTH()`/`R_xlen_t` across the codebase.

---

## 1) Coercion blocks (container‑first + indexed protects)

**Pattern (unchanged, but shorter):**

```c
{
  SEXP sx = s_x;
  PROTECT_INDEX px;
  PROTECT_WITH_INDEX(sx, &px);
  if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), px);

  const RLEN_T n = RLEN(sx);                 // LENGTH() by default
  std::vector<double> x(REAL(sx), REAL(sx) + (size_t)n);

  UNPROTECT(1); // fixed count
}
```

Why: rchk‑safe (no variable UNPROTECT), minimal ceremony, STL copy happens while protected.

---

## 2) Scalars / flags

- Use `Rf_asInteger`, `Rf_asReal`, `Rf_asLogical`.
- Validate early; if you `Rf_error()` **inside** a coercion block, first `UNPROTECT(1)` (the block’s protect).

---

## 3) Lengths and size math (LENGTH‑first)

- Prefer `int`/`R_len_t` from `LENGTH()` for element counts.
- Use `size_t` in STL calls and for buffer sizes **after** a safe cast.
- Guard multi‑dimensional products with a **checked multiply** to avoid overflow.

Helpers:

```c
static inline size_t checked_mul_size(size_t a, size_t b) {
  if (a != 0 && b > (SIZE_MAX / a)) Rf_error("size overflow");
  return a * b;
}

static inline size_t as_size_t_len(SEXP s) {
  int n = LENGTH(s);           // short vectors only
  if (n < 0) Rf_error("negative length");
#if GFLOW_ENABLE_LONG_VECTORS
  if ((R_xlen_t)n != XLENGTH(s))
    Rf_error("long vectors currently unsupported in this build");
#endif
  return (size_t)n;
}
```

Usage:

```c
const int nrow = LENGTH(s_mat_nrow);
const int ncol = LENGTH(s_mat_ncol);
const size_t N = checked_mul_size((size_t)nrow, (size_t)ncol);
SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
```

---

## 4) Result assembly (fixed UNPROTECT pattern)

```c
SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
{
  SEXP a = PROTECT(Rf_allocVector(REALSXP, n));
  // fill a...
  SET_VECTOR_ELT(res, 0, a);
  UNPROTECT(1);
}
{
  SEXP b = PROTECT(Rf_allocVector(INTSXP, m));
  // fill b...
  SET_VECTOR_ELT(res, 1, b);
  UNPROTECT(1);
}
SEXP nm = PROTECT(Rf_allocVector(STRSXP, 3));
SET_STRING_ELT(nm, 0, Rf_mkChar("a"));
SET_STRING_ELT(nm, 1, Rf_mkChar("b"));
SET_STRING_ELT(nm, 2, Rf_mkChar("c"));
Rf_setAttrib(res, R_NamesSymbol, nm);
UNPROTECT(2);
return res;
```

- Never `UNPROTECT(variable)`.
- Each temporary is paired locally with `UNPROTECT(1)`.

---

## 5) NULL / optional fields

```c
if (have_weights) {
  SEXP w = PROTECT(Rf_allocVector(REALSXP, n));
  // fill...
  SET_VECTOR_ELT(res, 2, w);
  UNPROTECT(1);
} else {
  SET_VECTOR_ELT(res, 2, R_NilValue);
}
```

---

## 6) Matrices (column‑major)

```c
SEXP M = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
double* p = REAL(M);
for (int j = 0; j < ncol; ++j)
  for (int i = 0; i < nrow; ++i)
    p[i + (size_t)j * (size_t)nrow] = ...;
SET_VECTOR_ELT(res, idx, M);
UNPROTECT(1);
```

---

## 7) Strings and names

- Build a protected `STRSXP` and assign with `SET_STRING_ELT`.
- No need to PROTECT individual `Rf_mkChar` calls that are immediately stored.

---

## 8) What we deliberately **don’t** do here

- No project‑wide requirement to use `R_xlen_t`/`XLENGTH()`.
- No attempt to support objects > 2^31−1 elements by default.
- No variable‑count `UNPROTECT`.

---

## 9) When to enable long vectors

If a specific entry point might reasonably handle a single vector exceeding 2^31−1 (e.g., edge lists in huge sparse graphs on HPC):

1. Compile with `-DGFLOW_ENABLE_LONG_VECTORS=1`.
2. Rebuild. All `RLEN_T`/`RLEN()` call sites automatically switch to `R_xlen_t`/`XLENGTH()`.

---

## Example — Before → After

**Before (complex long‑vector plumbing everywhere):**
```c
const R_xlen_t n = XLENGTH(s_x);
std::vector<double> x(REAL(s_x), REAL(s_x) + (size_t)n);
```

**After (pragmatic):**
```c
const int n = LENGTH(s_x);
std::vector<double> x(REAL(s_x), REAL(s_x) + (size_t)n);
```

Optionally, make it uniform with the macro:

```c
const RLEN_T n = RLEN(s_x); // expands to int/LENGTH by default
std::vector<double> x(REAL(s_x), REAL(s_x) + (size_t)n);
```

---

## rchk checklist

- [x] All allocating calls happen with protected objects.
- [x] Only fixed‑count UNPROTECTs.
- [x] No lingering PROTECTs across early errors.
- [x] Scalars read via Rf_as*.
- [x] Buffer arithmetic guarded with `checked_mul_size`.

This guideline keeps code lean and readable, while preserving a low‑friction path to long vectors if you ever need them.
