# What CRAN actually runs

* **R CMD check** on multiple platforms (Linux, macOS, Windows; devel/release/oldrel).
* **rchk** static analysis (Debian) → PROTECT/UNPROTECT, stale pointers, long-vector & coercion issues.
* **ASAN/UBSAN** (address/undefined sanitizers) on gcc/clang in some lanes.
* **valgrind** runs in some lanes (heap leaks/invalid reads/writes).
* **Windows UCRT** builds (catching MSVC/MinGW idiosyncrasies).
* **macOS clang** (stricter UB diagnostics).

# C/C++ entry-point checklist (rchk-safe + CRAN-friendly)

## Interface & argument handling

* Use **`.Call`** with **registered routines** (NAMESPACE `useDynLib`, and `R_registerRoutines` in `init.c`).
* Treat `.Call` args as **already PROTECTed**; still **defensively coerce scalars** with `Rf_asInteger/Real/Logical`.
* If you rely on an R wrapper for checks (e.g., matrix/double), that’s fine—just be consistent. Consider minimal C-side checks if users might call `.Call` directly.

## PROTECT discipline

* **Container-first**: `PROTECT` the outer list/vector first; create children → `SET_VECTOR_ELT` → `UNPROTECT(1)` for the child.
* **Fixed counts only**: avoid `UNPROTECT(nvar)` where `nvar` is computed at runtime.
* **No stale pointers**: never keep `REAL(x)`/`INTEGER(x)` pointers across any call that can allocate R objects.
* For conditional coercion, use **`PROTECT_WITH_INDEX`/`REPROTECT`**.
* If a helper returns a SEXP: pick a convention and stick to it:

  * **Recommended**: helper returns **PROTECTed** SEXP; caller anchors then `UNPROTECT(1)`.

## Numerics & ALTREP

* For scalars/flags use **`Rf_as*`** (not `REAL(x)[0]`/`INTEGER(x)[0]`).
* For long vectors, prefer `XLENGTH` (if you might exceed 2^31−1; otherwise OK to stay with `LENGTH`).
* Don’t assume `is.double(x)` implies contiguous memory unless you coerced it.

## Errors & exceptions

* No C++ exceptions may escape to R. Wrap and convert to `Rf_error(...)`.
* Remember `Rf_error` longjmps: ensure no resource leaks (e.g., RAII is fine; don’t rely on code after `Rf_error`).
* Use `Rprintf`/`Rf_warning` (not `printf`/`std::cerr`) for user messages.

## RNG & parallelism

* If you use R’s RNG, bracket with **`GetRNGstate()` / `PutRNGstate()`**.
* OpenMP/threads: respect CRAN’s limits; read thread counts from env vars (e.g., `OMP_NUM_THREADS`), provide a way to set threads from R, and avoid spawning unbounded threads.

## Memory & resources

* Don’t mix `R_alloc` and `free`. Use `R_alloc` for temp memory you don’t need to free, or `Calloc/Free` consistently.
* External pointers: register **finalizers**; ensure they’re robust and handle double-finalization.
* No writing outside temp dirs; no network access during checks.

## Portability gotchas

* Avoid non-portable intrinsics; guard with `#ifdef`.
* Watch for `int`/`size_t` conversions in `Rf_allocVector/Matrix` dims (cast as needed).
* Endianness rarely an issue unless you serialize.

## Output objects

* Ensure **1-based indices** when returning graph structures to R.
* Name your list components (reviewers look for clean, documented structure).
* Don’t return a still-PROTECTed top-level object—**unprotect it once** before `return`.

# Quick self-test before submit

* `R CMD check --as-cran` locally (all examples, vignettes).
* Run **`rchk`** locally if you can (or at least scan with `clang` static analyzer).
* On Linux: `ASAN/UBSAN` builds (`-fsanitize=address,undefined`), run tests.
* Optionally `valgrind R -d valgrind -f <testscript.R>`.


