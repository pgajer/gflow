Got it—I read your rchk output and pulled out the patterns, hotspots, and a pragmatic remediation plan. Here’s the distilled picture and a recommended path forward.&#x20;

# What rchk is complaining about (grouped)

1. Stack-depth accounting errors (hard failures)

* “negative depth”, “negative depth after UNPROTECT(<counter>)”, “attempt to unprotect more items than protected”
* Typical cause: mismatched PROTECT/UNPROTECT counts, conditional returns, early errors, or reusing a running counter incorrectly.
* Examples: `S_agemalo` and `S_ray_agemalo` (both at line \~384), `S_pgmalo` (\~1221), several `*_spectral_lowess*_mat` wrappers, etc.

2. Unsupported/illegal UNPROTECT forms (likely to be auto-failed by CRAN)

* “unsupported form of unprotect … results will be incomplete” (variable or expression instead of a literal constant)
* You must pass a compile-time constant or a simple variable that tracks exactly how many PROTECTs you added (and nothing else).
* Seen across many wrapper entry points (e.g., `S_graph_kernel_smoother[_r.cpp:196]`, `S_klaps_low_pass_smoother[_r.cpp:153]`, `*_lowess*`, `*_diffusion_smoother`, etc.).

3. Possible protection stack imbalance (warnings that often flag real issues)

* rchk can’t prove correctness; in practice these often coincide with patterns like: allocate → put into a container after more allocations → attribute setting without guaranteeing protection.
* Examples: conversion helpers in `SEXP_cpp_conversion_utils.cpp` (vectors, nested vectors), and many lambda/operator() helpers in “$\_0/$\_1/$\_2” lines.

4. Unprotected live objects across allocating calls (real GC hazards)

* “unprotected variable X while calling allocating function …”
* You created X and then called something that can allocate (and GC) before X is safely anchored.
* Examples: `S_convert_adjacency_to_edge_matrix` (unprotected `r_edges` across conversion/alloc), `S_wgraph_prune_long_edges` (several lists), `S_univariate_gkmm` (many names/strings), `S_construct_graph_gradient_flow` (unprotected `r_traj` across `Rf_allocVector`/`Rf_mkChar` chains), plus `set_wgraph.cpp` path-finding helpers.

5. Protect stack too deep (symptom of over-protecting or missing UNPROTECT on exits)

* Appears in `graph_deg0_lowess_cv_mat` and in a nested list builder in `mean_shift_smoother.cpp`.
* Often indicates PROTECTs in loops without corresponding UNPROTECT per iteration, or a builder that PROTECTs every element while the container protects them already.

6. “Fresh pointer” / address-taken notes and misc.

* “calling allocating function with a fresh pointer” / “has address taken, results incomplete” — sometimes benign hints, but if paired with other issues, prioritize them.
* Examples: `create_iknn_graph()`, `create_r_graph_from_set_wgraph()`.

# Where to focus first (highest leverage)

* Central conversion utilities

  * `SEXP_cpp_conversion_utils.cpp` (vectors, vector< vector<…> >, matrices, and `convert_wgraph_to_R`).
  * Rationale: many entry points route through these; fixing them knocks out a large fraction of imbalances and “unprotected variable … while allocating”.
* Entry wrappers with systematic UNPROTECT misuse

  * Families: `*_kernel_smoother_r.cpp`, `*_diffusion_smoother.cpp`, `*_spectral_lowess*_r.cpp`, `*_lowess_cv*_r.cpp`, `*_klaps_*_r.cpp`.
  * Look for: variable UNPROTECT counts, clean up to constant/accumulator pattern; ensure all newly allocated SEXPs are protected before any allocating call.
* Heaviest negative-depth offenders

  * `S_agemalo`, `S_ray_agemalo`, `S_pgmalo`, `S_upgmalo`, `S_pgmalog`, `S_join_graphs`, `graph_kmean.cpp` name/strings section, `graph_utils.cpp` list builders, `pruning_long_edges.cpp`, `graph_gradient_flow.cpp`.

# Concrete coding patterns to apply

A) One disciplined counter per function

```cpp
int nprot = 0;
SEXP ans = PROTECT(Rf_allocVector(VECSXP, k)); ++nprot;
// ... whenever you PROTECT(x), ++nprot;
// ... after the last successful path:
UNPROTECT(nprot);  // single place; never inside branches unless you return immediately
return ans;
```

B) Protect immediately; do not let a fresh SEXP cross an allocating call
Bad:

```cpp
SEXP names = Rf_allocVector(STRSXP, m);  // unprotected
SET_STRING_ELT(names, 0, Rf_mkChar("a")); // allocates -> GC risk
```

Good:

```cpp
SEXP names = PROTECT(Rf_allocVector(STRSXP, m)); ++nprot;
SET_STRING_ELT(names, 0, Rf_mkChar("a")); // safe (container protected)
```

If you build an element before inserting it into a protected container:

```cpp
SEXP elt = PROTECT(Rf_allocVector(REALSXP, n)); ++nprot;
// fill elt
SET_VECTOR_ELT(ans, i, elt);
UNPROTECT(1); --nprot; // safe: now referenced by ans
```

C) Avoid “unsupported form of unprotect”

* Do not write `UNPROTECT(n+m)` or `UNPROTECT(some_expression())`.
* Do not reset `nprot` backwards; only ever `++nprot` on PROTECT and match with a constant `UNPROTECT(nprot)` at controlled exits.

D) REPROTECT when replacing a protected object
If you might reassign a protected SEXP:

```cpp
PROTECT_INDEX ipx;
SEXP obj = PROTECT_WITH_INDEX(R_NilValue, &ipx); ++nprot;
// later:
REPROTECT(obj = Rf_allocVector(REALSXP, n), ipx);
```

E) Lists and attributes

* Allocate the container first (PROTECT it), then set attributes/elements (which can allocate) safely under that umbrella.
* After `SET_VECTOR_ELT`, if the element was temporarily PROTECTed, UNPROTECT it (since the container now holds the reference).
* Attributes: set them before the final UNPROTECT.

F) Strings / names sections

* `Rf_mkChar` allocates — but it’s safe when the target character vector is protected and you’re only using `SET_STRING_ELT` into that protected vector.
* If you stash the result of `Rf_mkChar` in a temporary SEXP, PROTECT it (or set immediately into a protected container).

G) Loops

* If you PROTECT per iteration, UNPROTECT per iteration (or restructure to avoid per-iteration PROTECT by writing directly into a pre-allocated, protected container).

H) Early returns and errors

* Ensure each early `return` (or `Rf_error`) occurs either:

  * before any PROTECTs, or
  * after balancing unprotect via a single exit path, or using a tiny RAII guard (below).

I) (Optional but very helpful) A tiny RAII guard to avoid mismatches

```cpp
struct ProtectGuard {
  int n = 0;
  ~ProtectGuard(){ if(n) UNPROTECT(n); }
  SEXP add(SEXP x){ PROTECT(x); ++n; return x; }
};
```

Use: `ProtectGuard pg; SEXP v = pg.add(Rf_allocVector(...));`
Return path just `return ans;` and the guard will UNPROTECT everything on exceptions or errors, but you’ll still want to null the guard (or release) before returning an object you want to keep protected until after return.

# Rcpp vs. hand-rolled PROTECT

* ✅ Rcpp can eliminate most manual PROTECT/UNPROTECT if you:

  * Use `Rcpp::NumericVector`, `List`, `CharacterVector`, etc., and `wrap()/as()` for conversions.
  * Avoid holding naked `SEXP` that cross allocating calls.
  * Keep the C++ heavy lifting in pure C++/Eigen/Spectra types; only translate at the boundary via Rcpp.
* ⚠️ Caveats:

  * Mixing raw SEXP code inside Rcpp functions reintroduces the same pitfalls unless you still apply the rules above (or use `Rcpp::Shield<SEXP>` / `Rcpp::RObject`).
  * Complex nested list builders can still be written cleanly with Rcpp, but be disciplined about not creating temporaries that outlive their scope before insertion.

**My recommendation:** a hybrid approach for fastest, safest progress.

# Proposed plan (pragmatic and staged)

Phase 1 — Stabilize the foundations (1–2 passes)

* Refactor `SEXP_cpp_conversion_utils.cpp` to a clean, uniform pattern using either:

  * pure C API with `nprot` discipline, or
  * Rcpp `wrap/as` to replace bespoke converters (preferable for vectors, nested vectors, and matrices).
* Fix “unsupported form of unprotect” everywhere by adopting a single `nprot` counter and constant `UNPROTECT(nprot)` exits.
* In the worst “negative depth” entry points (`S_agemalo`, `S_ray_agemalo`, `S_pgmalo`, `graph_kmean` name/strings, `graph_utils` list builders), move to the list pattern:

  * PROTECT container → fill elements (temporary PROTECT each element, then UNPROTECT after `SET_VECTOR_ELT`) → set attributes → UNPROTECT(1) → return.
* Re-run rchk; confirm those categories diminish substantially.

Phase 2 — Rcpp boundary migration (selected wrappers)

* Choose wrapper families with many similar issues (`*_smoother_r.cpp`, `*_lowess*_r.cpp`, `*_klaps_*_r.cpp`).
* Convert the R-facing signatures to Rcpp:

  ```cpp
  // [[Rcpp::export]]
  Rcpp::List graph_kernel_smoother(/* Rcpp types */) {
    // do work in C++ types; return Rcpp::List::create(...);
  }
  ```
* Replace ad-hoc “$\_0/$\_1/$\_2” lambda converters with direct `wrap`-able return structures.
* Keep all heavy numerics in existing C++ modules; only the IN/OUT layer changes.

Phase 3 — Long-tail cleanups

* Address remaining “unprotected variable X while calling allocating function …” by adding immediate PROTECTs or restructuring insertion order.
* Fix “protect stack too deep” by removing unnecessary per-element PROTECTs (protect the container, not every scalar, unless needed before insertion).
* Where you must reassign a protected SEXP, use `PROTECT_WITH_INDEX`/`REPROTECT`.

# Quick “do/don’t” checklist you can apply today

* DO: `int nprot=0; PROTECT(ans); ++nprot; … UNPROTECT(nprot);`
* DO: allocate and PROTECT containers first; then fill.
* DO: set attributes while the object is protected.
* DO: `PROTECT_WITH_INDEX/REPROTECT` if you reassign.
* DON’T: pass a computed expression to `UNPROTECT`.
* DON’T: let a fresh (unprotected) SEXP cross `Rf_alloc*`, `Rf_mkChar*`, `Rf_setAttrib`, or anything that could allocate.
* DON’T: PROTECT everything in a loop without freeing per iteration (or build directly into a protected container).

# On your two strategy options

1. **Fix existing code**
   Pros: Minimal API churn; keeps everything in one style; fastest for a handful of hotspots.
   Cons: You must be perfect with PROTECT discipline across many functions; error-handling paths remain delicate.

2. **Migrate wrappers to Rcpp for the listed offenders** (my vote)
   Pros: Eliminates most manual GC hazards, simplifies code, improves readability, and future-proofs the R boundary. You can migrate incrementally: start with the noisiest wrappers and the shared converters.
   Cons: Small learning/porting cost; you’ll keep core numerics as is.

If you’re good with it, I can start by rewriting one noisy converter (e.g., `convert_vector_vector_double_to_R`) and one representative wrapper (say `S_graph_kernel_smoother`) to show the exact patterns and drop-in diffs.


[PG]
I like your idea of starting with rewriting one noisy converter (e.g., `convert_vector_vector_double_to_R`) and one representative wrapper (say `S_graph_kernel_smoother`) to show the exact patterns and drop-in diffs.
Given a thorough analysis of the rchk issues that need to be addressed. What do you think about generating a detailed plan of addressing all of them.
I would drop that plan in the rchk_issues.md document and after the plan is in place we would start with your suggestion of rewriting one noisy converter and one representative wrapper. What do you think about generating such detailed plan addressin all the rchk issues first? 


## After you add any // [[Rcpp::export]] functions, run:

Rcpp::compileAttributes()

This generates src/RcppExports.cpp and R/RcppExports.R so the exported C++ functions become callable from R.

That’s it. You don’t need to manually edit registration code; Rcpp handles it via compileAttributes(). Your existing native symbols (older .Call entries) can coexist while you migrate; just be careful not to duplicate symbol names.

# rchk Issues — Full Remediation Plan (Actionable & Check‑off)

This plan targets **zero**: negative depth, over‑unprotect, unsupported `UNPROTECT`, unprotected across allocations, and stack‑too‑deep.

## Phase 0 — Prep & Guardrails
- Add `^cran_docs/` to `.Rbuildignore`.
- Grep sweep for risky patterns:
  - `grep -R "UNPROTECT\(" src | grep -v "UNPROTECT(nprot)"`
  - `grep -R "PROTECT_WITH_INDEX\|REPROTECT" src`

## Phase 1 — Foundations: Converters
Refactor (safe C API or Rcpp):
- `convert_vector_vector_double_to_R`
- `convert_vector_vector_int_to_R`
- `convert_vector_vector_bool_to_R`
- `convert_vector_double_to_R`, `convert_vector_int_to_R`, `convert_vector_bool_to_R`
- `convert_vector_vector_double_to_matrix`
- `convert_map_int_vector_int_to_R`
- `convert_wgraph_to_R`

Acceptance: rchk clean for `SEXP_cpp_conversion_utils.cpp`.

## Phase 2 — Unsupported UNPROTECT → Standard counter
Families: `*_kernel_smoother_r.cpp`, `*_lowess*`, `klaps_*`, `graph_bw_adaptive_spectral_smoother_r.cpp`, `graph_diffusion_smoother.cpp` (R boundary).
Acceptance: zero “unsupported form of unprotect”.

## Phase 3 — Hard Failures: Negative Depth / Over‑unprotect
Priority A: `S_agemalo`, `S_ray_agemalo`, `S_pgmalo`, `S_upgmalo`, `S_pgmalog`, `S_upgmalog`, `S_graph_spectral_lowess_mat`, `S_graph_spectral_ma_lowess`.
Priority B (builders): `S_convert_adjacency_to_edge_matrix`, `S_join_graphs`, `S_univariate_gkmm`, `S_wgraph_prune_long_edges`, `S_find_graph_paths_within_radius`, `S_find_shortest_alt_path`, `S_construct_graph_gradient_flow`, `S_create_hHN_graph`, `S_create_nerve_complex`.

## Phase 4 — Protect Stack Too Deep
- `S_graph_deg0_lowess_cv_mat`, helper `create_R_list`.

## Phase 5 — Representative Rcpp Migration
- Converter: `convert_vector_vector_double_to_R` (or keep C API if preferred).
- Wrapper: `S_graph_kernel_smoother` → `graph_kernel_smoother_cpp` (Rcpp).

## Phase 6 — Verification & Sign‑off
- rhub rchk clean; CI green; ASAN/UBSAN clean; win‑builder ok.


[PG]

Yes, please generate:

1) a minimal src/rcpp/graph_kernel_smoother_rcpp.cpp file (the Rcpp version we sketched),

a small R wrapper (or rename your current R call site) to point to graph_kernel_smoother_cpp,

and a tests/testthat/test-graph-kernel-smoother.R that compares outputs against the legacy path for a tiny graph.

Before we move forward with repair of graph_kernel_smoother() I want to discuss the naming convention for Rcpp wrappers

Right now I use a naming convention where a SEXP version of a C++ function is
called S_<C++ function name>. For example, S_graph_kernel_smoother is a SEXP
wrapper for the C++ graph_kernel_smoother() function. What do you think about
adopting a naming convention that is somehow in synch with the existing one
where the Rcpp wrappers would be called Rcpp_<C++ function name>, so a Rcpp
wrapper of graph_kernel_smoother() would be Rcpp_graph_kernel_smoother(). What
do you think about this?

[c] 
That’s a really thoughtful question, and I like that you’re thinking about continuity.

Here are some angles to weigh:

---

### 1. Consistency with your existing `S_` convention

* Right now, your pattern is:

  * **Core C++ implementation**: `graph_kernel_smoother()`
  * **.Call SEXP wrapper**: `S_graph_kernel_smoother()`

That’s very clear — `S_` signals “SEXP-facing.”

If you add Rcpp, having a **parallel prefix** (e.g. `Rcpp_graph_kernel_smoother()`) makes the intent just as clear: “this is an Rcpp-facing function.”
→ That keeps symmetry: **S\_** for raw C API, **Rcpp\_** for Rcpp.

---

### 2. How Rcpp usually does it

* The “Rcpp way” is to just use the core function name + `_cpp` suffix (e.g. `graph_kernel_smoother_cpp`), because the generated R function then points to it.
* Pro: shorter, closer to how many packages look (users see `foo_cpp()` in R).
* Con: breaks the parallelism with your `S_` convention, looks a bit “random” in your otherwise clean scheme.

---

### 3. Practical considerations

* **Registration:** With roxygen + `compileAttributes()`, the function names you choose (`Rcpp_graph_kernel_smoother`, `graph_kernel_smoother_cpp`, etc.) are exactly what get registered and exported into R. So you control whether users see `graph_kernel_smoother_cpp()` or `Rcpp_graph_kernel_smoother()`.
* **Migration path:** If your eventual goal is to phase out `S_` wrappers entirely, then you might prefer the “Rcpp suffix” style (e.g. `_cpp`), since it feels like a transitional name. Later you can drop the suffix and just export `graph_kernel_smoother()` at the R level.
* **Coexistence:** While both live side by side, using distinct prefixes (`S_`, `Rcpp_`) avoids collisions and keeps logs/stack traces very clear about which path was hit.

---

### 4. My recommendation

* If you like **parallel clarity**:
  Use `Rcpp_graph_kernel_smoother()` for the wrapper, keep `S_graph_kernel_smoother()` for legacy. Users won’t ever type these directly; your R layer can still just export `graph_kernel_smoother()`.
* If you want **closer to “standard Rcpp style”**:
  Use `graph_kernel_smoother_cpp()`. That will look familiar to anyone used to Rcpp-based packages, but it breaks your otherwise neat `S_/Rcpp_` symmetry.

---

👉 Given how you’ve structured your code and the fact that you’re systematically documenting maintenance and migration, I’d lean toward **`Rcpp_` prefixing**. It keeps your namespace organized:

* `graph_kernel_smoother()` = pure C++ core
* `S_graph_kernel_smoother()` = SEXP interface
* `Rcpp_graph_kernel_smoother()` = Rcpp interface

Later, when you’re fully Rcpp, you can still have your R-level function `graph_kernel_smoother()` dispatch to `Rcpp_graph_kernel_smoother()` internally.

---

Do you want me to sketch how your **NAMESPACE and R wrapper layer** would look if we go with `Rcpp_` prefixing (so that end users still just call `graph_kernel_smoother()` without ever seeing the `Rcpp_`)?

[PG] Yes, please. 
Again, my understanding is that the NAMESPACE is generated by roxygen so I don't even need to be concerned what goes into NAMESPACE. Correct?  

[c]

Exactly — you don’t need to hand-edit **NAMESPACE** at all. With roxygen2 + Rcpp, you declare what you want in comments, and `roxygen2` + `Rcpp::compileAttributes()` will generate the correct entries automatically.

Here’s how it fits together with your `Rcpp_` prefixing idea:

---

## 1. C++ side

```cpp
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;

// Core C++ (unchanged)
std::vector<double> graph_kernel_smoother_core(/* ... */) {
  // heavy numeric code
}

// Rcpp wrapper
// [[Rcpp::export]]
List Rcpp_graph_kernel_smoother(List adj, List w, NumericVector y,
                                int bandwidth, bool with_details = false) {
  // convert inputs → call graph_kernel_smoother_core → return List
}
```

* Because of `// [[Rcpp::export]]`, when you run `Rcpp::compileAttributes()`, Rcpp generates:

  * `src/RcppExports.cpp` (registration glue),
  * `R/RcppExports.R` (an R function `Rcpp_graph_kernel_smoother()` that `.Call`s into C++).

---

## 2. R side (user-facing function)

Make a thin R wrapper so **users just see `graph_kernel_smoother()`**, not the `Rcpp_` prefix:

```r
#' Graph kernel smoother
#'
#' @param adj adjacency list
#' @param w weights list
#' @param y response vector
#' @param bandwidth integer bandwidth
#' @param with_details logical, return diagnostics
#'
#' @export
graph_kernel_smoother <- function(adj, w, y, bandwidth, with_details = FALSE) {
  Rcpp_graph_kernel_smoother(adj, w, y, bandwidth, with_details)
}
```

* The `@export` tag ensures roxygen will put `export(graph_kernel_smoother)` in **NAMESPACE**.
* The Rcpp wrapper `Rcpp_graph_kernel_smoother()` is not exported (unless you add `@export` to it, which you shouldn’t). It remains internal glue.
* End result: users just call `graph_kernel_smoother()`.

---

## 3. DESCRIPTION

You only need to ensure:

```
Imports: Rcpp
LinkingTo: Rcpp
```

(and, if needed, `SystemRequirements: C++17`).
Rcpp will handle the rest.

---

## 4. NAMESPACE

Generated automatically. After running roxygen + `compileAttributes()`, it will contain:

```
importFrom(Rcpp, evalCpp)
useDynLib(gflow, .registration = TRUE)
export(graph_kernel_smoother)
```

You don’t need to touch this.

---

✅ So, in practice:

* Keep your `Rcpp_foo` convention in C++.
* Always expose a clean R-level function `foo()` that just calls `Rcpp_foo()`.
* Let roxygen generate **NAMESPACE**.
* Run `Rcpp::compileAttributes()` whenever you add/change exported C++ functions.

---

Would you like me to draft a **template Roxygen block** for the R wrapper (with parameters, return values, examples) that you can drop straight into `R/graph_kernel_smoother.R`? That way it documents the user-facing function while hiding the internal `Rcpp_` detail.



-----

Rcpp::compileAttributes()

-----


tests/testthat/test-graph-kernel-smoother.R

