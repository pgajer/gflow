# mabilo issues

# What’s crashing (root cause)

Your Bayesian-bootstrap path runs **many wmabilo() fits in parallel**. Each wmabilo() (and also the smoothed variant) calls `initialize_kernel(distance_kernel, …)` to set up a **global kernel function** that is later used via `kernel_fn(...)`. Because multiple threads are calling `initialize_kernel()` concurrently, they race on that global state. That kind of data race is classic “heap corruption now, crash later” — which is exactly what you see (malloc\_zone\_error during an unrelated graphics allocation). You can see:

* parallel loop for bootstrap in `mabilo_bb()` (gflow::for\_each with `GFLOW_EXEC_POLICY`) and per-iteration call into `wmabilo(...)`&#x20;
* each wmabilo() / smoothed path re-`initialize_kernel(...)` before using `kernel_fn(...)` (global)&#x20;

**Fix options (pick one now, then refactor later):**

1. **Immediate safety (quickest)**: make the bootstrap sequential (use a sequential exec policy) until the kernel state is thread-safe.

2. **Make kernel state thread-local**: refactor `kernels.h` so `initialize_kernel(...)` returns a callable (lambda/functor) and pass that down; do **not** set any global function pointer. In `mabilo_bb()`, construct the kernel functor once per thread or (better) once outside and capture it by value in the worker.

3. **Coarse lock (stop-gap)**: wrap every call that mutates or reads the global kernel state with the same mutex. This preserves parallelism poorly (heavy contention), but avoids corruption.

Until you adopt #2, I’d strongly recommend flipping the bootstrap to sequential to confirm the crash disappears.

# Other correctness issues that can also bite

1. **Length check missing in S\_mabilo()**
   `S_wmabilo()` validates `nx == ny == nw`. `S_mabilo()` does **not** check `nx == ny`, which means you can flow mismatched vectors into the core and access out-of-bounds memory. Add the same length guard you already use in the weighted entry point.&#x20;

2. **Kernel code mismatch between R and C**
   Your R interface advertises/accepts `kernel.type` in **1..10** (Epanechnikov, Triangular, … Tricube default=7, etc.). The C/C++ side documents/uses a `distance_kernel` with **0..2** (Tricube, Epanechnikov, Exponential). Passing the default 7 from R into C will select an **undefined kernel id** unless you remap it; depending on how `initialize_kernel()` handles unknown codes, that can produce nonsense weights or hit undefined behavior. Either:

* map R’s 1..10 into C’s 0..2 in the R wrapper before `.Call(...)`, or
* extend `initialize_kernel()` to support all 10 and keep one consistent enumeration end-to-end.&#x20;

3. **Indexing inconsistency for `opt_k_idx`**
   `S_wmabilo()` returns `opt_k_idx` “as is”, while `S_mabilo()` returns `opt_k_idx + 1` (R-1-based). This inconsistency won’t crash, but it’s a nasty footgun for R callers. Make both R-level entry points return 1-based indexes.&#x20;

4. **Type/length handling of `y.true` from R**
   Your R wrapper sets `y.true <- numeric(0)` when missing; the C side treats a length-0 real vector as “not present” (it only populates `y_true` when `nyt == nx`). That’s fine (and intentional), but if you want the R list field to be `NULL` rather than a length-0 numeric, pass `NULL` from R when `y.true` is not supplied (minor polish).&#x20;

5. **Thread safety of RNG**
   You protected `C_runif_simplex()` with a mutex (good), but the **only** protected piece is RNG. Because kernel initialization is still global (unprotected), the code remains unsafe (root cause above).&#x20;

# R interface review (the .R function you posted)

The high-level R function looks mostly solid (input checks, sorting, parameter defaults), but two practical issues:

* **`kernel.type` default 7** conflicts with the C side (see mismatch above). If you’re not extending C to support 10 kernels immediately, change the R default to “Tricube” in the C code’s domain (e.g., 0) or remap in R before the `.Call`.
* Consider passing `NULL` (not `numeric(0)`) for `y.true` when absent; it keeps R/C semantics a touch clearer (and will make `k_mean_true_errors` reliably `NULL` in C without depending on a length check).

# Why your simple example crashed

```
x <- seq(0, 10, length.out = 100)
y <- sin(x) + rnorm(100, 0, 0.1)
fit <- mabilo(x, y, n.bb = 100)
```

* `n.bb=100` → `mabilo_bb()` launches many parallel wmabilo() fits.
* Each fit races on `initialize_kernel(...)`/`kernel_fn` global state → **heap corruption** → later **malloc abort** while R tries to start a new plot, which aligns with your backtrace (abort in `Rf_allocVector3` called from `graphics`). Fixing the kernel state thread-safety (or making the bootstrap sequential) should make this exact crash disappear.&#x20;

# Concrete patch suggestions

**Short term (quick, safe):**

* In `mabilo_bb(...)`, set a sequential policy while you land the kernel refactor:

  ```cpp
  gflow::for_each(std::execution::seq, bb_indices.begin(), bb_indices.end(), ...);
  ```
* OR call `initialize_kernel(...)` **once** outside the parallel region and make `kernel_fn` a **captured callable** that does not rely on global mutable state. (Even better: have `initialize_kernel` return a stateless functor and pass it through to `wmabilo`/`ulm`.)

**Also fix these now:**

* Add `if (nx != ny) Rf_error("x and y must have the same length.");` to `S_mabilo()` (you already did this pattern in `S_wmabilo()`).&#x20;
* Make `opt_k_idx` consistently 1-based in **both** `S_mabilo()` and `S_wmabilo()` at the R boundary.&#x20;
* Align `kernel.type` enumerations between R and C (temporary remap in R is acceptable if you don’t want to expand the C kernels immediately).&#x20;

# Side note on the MST completion file

Your MST/CMST R interface (`S_create_mst_completion_graph`) and helpers look PROTECT-balanced and correct. `create_r_graph_from_set_wgraph()` now unprotects its two top-level lists internally; the caller immediately PROTECTs them again before assembling the return list. That’s correct and won’t cause this crash.&#x20;

---

If you want, I can sketch the kernel refactor (function-object based, no globals) so the bootstrap stays parallel and safe.


