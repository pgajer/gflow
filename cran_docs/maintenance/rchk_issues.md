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

[~/current_projects/gflow/src]% ls -1 *_lowess*_r.cpp
deg0_lowess_graph_smoothing_r.cpp
graph_deg0_lowess_buffer_cv_r.cpp
graph_deg0_lowess_cv_mat_r.cpp
graph_spectral_lowess_mat_r.cpp
graph_spectral_lowess_r.cpp
graph_spectral_ma_lowess_r.cpp
nada_graph_spectral_lowess_r.cpp
spectral_lowess_graph_smoothing_r.cpp

klaps_low_pass_smoother_r.cpp
graph_bw_adaptive_spectral_smoother_r.cpp
graph_diffusion_smoother.cpp





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



tests/testthat/test-graph-kernel-smoother.R



------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
# ToDo's: 
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

1) Find all instances of 

convert_vector_double_to_R
convert_vector_int_to_R
convert_vector_bool_to_R

convert_vector_vector_int_to_R
convert_vector_vector_double_to_R
convert_vector_vector_bool_to_R

convert_vector_vector_double_to_matrix
convert_map_int_vector_int_to_R
convert_wgraph_to_R

calls and make sure the output is UNPROTECTED 

Great question. Short answer: the converter returns an **unprotected** SEXP by design. That’s correct and CRAN-safe—**as long as you immediately anchor it** in one of the standard ways.

Here’s exactly how to use `convert_vector_double_to_R()` safely in all common situations:

---

### 0) What the converter does

```cpp
SEXP convert_vector_double_to_R(const std::vector<double>& vec) {
  int nprot = 0;
  SEXP Rvec = PROTECT(Rf_allocVector(REALSXP, vec.size())); ++nprot;
  // ... fill ...
  UNPROTECT(nprot);              // <- unprotected on return (intended)
  return Rvec;
}
```

* It **PROTECTs during construction** so GC can’t touch it while being filled.
* It **UNPROTECTs before return** so the caller decides how to anchor it.
* When a `.Call` function returns its result to R, the **R evaluator preserves the return value**—so returning it directly is safe.

---

### 1) Case A — Return it directly from a `.Call` wrapper (no further allocations)

Safe and simplest: no extra PROTECT needed.

```cpp
extern "C" SEXP S_make_vec(SEXP /*args*/) {
  // do not allocate anything after this call if you rely on returning it directly
  return convert_vector_double_to_R(my_cpp_vec);
}
```

---

### 2) Case B — You’ll do more allocations after getting it (e.g., building a list)

**Immediately PROTECT** the returned object (or insert it into a protected container and then UNPROTECT).

Conservative, rchk-friendly pattern:

```cpp
extern "C" SEXP S_build_list(SEXP /*args*/) {
  int nprot = 0;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 1)); ++nprot;

  SEXP elt = PROTECT(convert_vector_double_to_R(my_cpp_vec)); ++nprot; // PROTECT immediately
  SET_VECTOR_ELT(out, 0, elt);
  UNPROTECT(1); --nprot;                    // element now referenced by 'out'

  // ... you can safely make more allocations here ...

  UNPROTECT(nprot);
  return out;
}
```

Why PROTECT the element even though you set it immediately?

* It’s the **safest** pattern: if anything allocates between creation and `SET_VECTOR_ELT` (directly or via inline helpers), GC can’t collect it.

---

### 3) Case C — Build container first, then set elements (per-element protect/unprotect)

This is the “container-first” idiom we’ve been standardizing:

```cpp
SEXP out = PROTECT(Rf_allocVector(VECSXP, K)); ++nprot;

SEXP v = PROTECT(convert_vector_double_to_R(vec)); ++nprot;
// ... optional transforms that might allocate ...
SET_VECTOR_ELT(out, i, v);
UNPROTECT(1); --nprot; // v now owned by 'out'
```

---

### 4) Case D — Rcpp boundary (use a shield)

If you’re inside an Rcpp function but calling a raw-C converter:

```cpp
Rcpp::Shield<SEXP> v(convert_vector_double_to_R(vec));
Rcpp::List out(1);
out[0] = v;  // safe; shield releases when v goes out of scope
```

---

### 5) Long-lived objects across calls (rare)

If you need to hold onto the SEXP **beyond the current `.Call`** (e.g., store globally, use later):

* Use `R_PreserveObject(sexp);` when you obtain it, and **later** `R_ReleaseObject(sexp);`.
* Do **not** rely on PROTECT/UNPROTECT for multi-call lifetimes; those are only for the current C stack.

---

### 6) What **not** to do

* Don’t let the **unprotected** return value cross any allocating call (e.g., `Rf_alloc*`, `Rf_mkChar`, `Rf_setAttrib`) before you PROTECT or insert into a protected container.
* Don’t return an object that’s still PROTECTed **without** balancing the protect stack (it will trip rchk and can cause “protect stack too deep”).

---

### TL;DR

* Returning it straight to R? **Just return** (safe).
* Doing anything else before returning? **PROTECT it immediately**, then either:

  * insert into a **protected** container and `UNPROTECT(1)`, or
  * keep it PROTECTed until your final `UNPROTECT(nprot)` just before `return`.

If you show me a specific call site where you use the converter, I’ll mark exactly where to PROTECT/UNPROTECT to keep rchk perfectly quiet.


    


-------- checking examples

# 1) Build the same "all examples" script that R CMD check runs
rds <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
con <- file("gflow-Ex.R", "w", encoding = "UTF-8")
for (f in rds) {
  rd <- tools::parse_Rd(f, encoding = "UTF-8")
  tools::Rd2ex(rd, out = con)
}
close(con)

# 2) Now run it to reproduce the crash
source("gflow-Ex.R")


rds <- list.files("man", pattern="\\.Rd$", full.names=TRUE)
env <- new.env(parent = baseenv())
cat("Running examples in ONE session...\n")
for (i in seq_along(rds)) {
  f <- rds[i]
  rd <- tools::parse_Rd(f, encoding="UTF-8")
  tf <- tempfile(fileext = ".R")
  tools::Rd2ex(rd, out = tf)
  cat(sprintf("\n==== [%d/%d] %s ====\n", i, length(rds), basename(f)))
  flush.console()
  # Optional: make GC more aggressive to catch PROTECT bugs earlier
  # gctorture2(enable = TRUE, inhibit_release = TRUE)
  source(tf, local = env)  # if it aborts, the last header printed is your offender
  # gctorture2(enable = FALSE)
}
cat("\nAll examples completed without abort in single session.\n")


for f in man/*.Rd; do
  echo ">>> Running $(basename "$f")"
  Rscript -e "rd<-tools::parse_Rd('$f', encoding='UTF-8');
              tf<-tempfile(fileext='.R'); tools::Rd2ex(rd, out=tf);
              source(tf)" >/tmp/ex.log 2>&1 || { echo 'FAILED'; exit 1; }
done



# In package root
Sys.setenv(R_ABORT_ON_PROTECT_ERROR = "1")  # make PROTECT bugs fatal where they occur

# Turn only "stack imbalance" warnings into errors and capture a dump
withCallingHandlers(
  {
    # run the same examples script
    source("gflow-Ex.R")
  },
  warning = function(w) {
    if (grepl("stack imbalance", conditionMessage(w), fixed = TRUE)) {
      dump.frames("stack_imbalance_dump", to.file = TRUE)
      stop(w)  # convert to error so we halt right at the offender
    }
  }
)


# --- Strict settings so we fail at the origin ---
Sys.setenv(R_ABORT_ON_PROTECT_ERROR = "1")
options(warn = 1)  # print warnings in order (not buffered)

# We'll write the dump to a guaranteed path and print it.
dump_path <- file.path(getwd(), "last.dump.rda")

# Evaluate gflow-Ex.R line-by-line to pinpoint the first imbalance.
lines <- readLines("gflow-Ex.R", warn = FALSE)
env <- new.env(parent = baseenv())

cat(sprintf("Running %d lines from gflow-Ex.R ...\n", length(lines)))
for (i in seq_along(lines)) {
  withCallingHandlers(
    {
      expr <- parse(text = lines[i])
      if (length(expr)) eval(expr, envir = env)
    },
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("stack imbalance", msg, fixed = TRUE)) {
        cat("\n### STACK IMBALANCE at line", i, "###\n")
        cat("Message: ", msg, "\n", sep = "")
        cat("Context (prev..next):\n",
            paste0("  ", lines[pmax(1, i-5):pmin(length(lines), i+5)]),
            sep = "\n")
        # Create a dump you can inspect with debugger()
        dump.frames("last.dump", to.file = TRUE)
        cat("\nWrote dump to: ", dump_path, "\n", sep = "")
        stop(w)  # convert to error so we halt right here
      } else {
        invokeRestart("muffleWarning")
      }
    }
  )
}
cat("\nAll lines completed with no stack imbalance.\n")


# Strict settings
Sys.setenv(R_ABORT_ON_PROTECT_ERROR = "1")
options(warn = 1)  # warnings are immediate

# Attach default packages like an interactive session
for (pkg in c("stats","graphics","grDevices","utils","datasets","methods")) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Ensure your package is attached (exports visible like in R CMD check)
suppressPackageStartupMessages(library(gflow))

# Prepare runner
dump_path <- file.path(getwd(), "last.dump.rda")
lines <- readLines("gflow-Ex.R", warn = FALSE)
env <- new.env(parent = .GlobalEnv)  # <-- key: inherits normal search path

cat(sprintf("Running %d lines from gflow-Ex.R ...\n", length(lines)))
for (i in seq_along(lines)) {
  withCallingHandlers(
    {
      expr <- parse(text = lines[i])
      if (length(expr)) eval(expr, envir = env)
    },
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("stack imbalance", msg, fixed = TRUE)) {
        cat("\n### STACK IMBALANCE at line", i, "###\n")
        cat("Message: ", msg, "\n", sep = "")
        cat("Context (prev..next):\n",
            paste0("  ", lines[pmax(1, i-5):pmin(length(lines), i+5)]),
            sep = "\n")
        dump.frames("last.dump", to.file = TRUE)  # -> last.dump.rda
        cat("\nWrote dump to: ", dump_path, "\n", sep = "")
        stop(w)
      } else {
        invokeRestart("muffleWarning")
      }
    }
  )
}
cat("\nAll lines completed with no stack imbalance.\n")

# --- strict settings ---
Sys.setenv(R_ABORT_ON_PROTECT_ERROR = "1")
options(warn = 1)

# Attach default base packages & gflow so examples see them
for (pkg in c("stats","graphics","grDevices","utils","datasets","methods"))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
suppressPackageStartupMessages(library(gflow))

# Parse the entire script with source refs
ex_file <- "gflow-Ex.R"
exprs <- parse(file = ex_file, keep.source = TRUE)
lines <- readLines(ex_file, warn = FALSE)
env <- new.env(parent = .GlobalEnv)
dump_path <- file.path(getwd(), "last.dump.rda")

cat(sprintf("Running %d expressions from %s ...\n", length(exprs), ex_file))
for (i in seq_along(exprs)) {
  withCallingHandlers(
    {
      eval(exprs[[i]], envir = env)
    },
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("stack imbalance", msg, fixed = TRUE)) {
        sr <- attr(exprs[[i]], "srcref")
        if (!is.null(sr)) {
          start <- sr[[1]]; end <- sr[[3]]
          cat("\n### STACK IMBALANCE at expression", i,
              sprintf("(lines %d-%d)", start, end), "###\n", sep = " ")
          ctx <- lines[pmax(1, start-5):pmin(length(lines), end+5)]
          cat("Message: ", msg, "\n", sep = "")
          cat("Context:\n", paste0("  ", ctx), sep = "\n")
        } else {
          cat("\n### STACK IMBALANCE at expression", i, "(no srcref) ###\n")
          cat("Message: ", msg, "\n", sep = "")
        }
        dump.frames("last.dump", to.file = TRUE)  # -> last.dump.rda
        cat("\nWrote dump to: ", dump_path, "\n", sep = "")
        stop(w)
      } else {
        invokeRestart("muffleWarning")
      }
    }
  )
}
cat("\nAll expressions completed with no stack imbalance.\n")



# Strict settings
Sys.setenv(R_ABORT_ON_PROTECT_ERROR = "1")
options(warn = 2)  # turn warnings into errors so we can trap them

# Attach default packages + gflow (like an interactive session)
for (pkg in c("stats","graphics","grDevices","utils","datasets","methods")) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
suppressPackageStartupMessages(library(gflow))

# Parse all expressions with source references
ex_file <- "gflow-Ex.R"
exprs <- parse(file = ex_file, keep.source = TRUE)
lines <- readLines(ex_file, warn = FALSE)
env <- new.env(parent = .GlobalEnv)

cat(sprintf("Running %d expressions from %s ...\n", length(exprs), ex_file))

handle_fail <- function(i, cond) {
  msg <- conditionMessage(cond)
  sr  <- attr(exprs[[i]], "srcref")
  if (!is.null(sr)) {
    start <- sr[[1]]; end <- sr[[3]]
    cat("\n### FAILURE at expression", i, sprintf("(lines %d-%d)", start, end), "###\n")
    cat("Message: ", msg, "\n", sep = "")
    ctx <- lines[pmax(1, start-5):pmin(length(lines), end+5)]
    cat("Context:\n", paste0("  ", ctx), sep = "\n")
  } else {
    cat("\n### FAILURE at expression", i, "(no srcref) ###\nMessage: ", msg, "\n", sep = "")
  }
  dump.frames("last.dump", to.file = TRUE)  # -> last.dump.rda in cwd
  cat("\nWrote dump to: ", file.path(getwd(), "last.dump.rda"), "\n", sep = "")
}

for (i in seq_along(exprs)) {
  # Catch both warnings-as-errors and regular errors
  tryCatch(
    {
      withCallingHandlers(
        {
          eval(exprs[[i]], envir = env)
        },
        warning = function(w) {
          # If warn=2, this usually won't fire; but keep it as belt-and-suspenders
          if (grepl("stack imbalance", conditionMessage(w), fixed = TRUE)) {
            handle_fail(i, w); stop(w)
          } else {
            invokeRestart("muffleWarning")
          }
        }
      )
    },
    error = function(e) {
      # We stop on ANY error now; report and dump for inspection
      handle_fail(i, e)
      stop(e)
    }
  )
}

cat("\nAll expressions completed with no stack imbalance.\n")


-----

# Strict like R CMD check
Sys.setenv(
  R_ABORT_ON_PROTECT_ERROR = "1",  # fatal on PROTECT bugs
  `_R_CHECK_INTERNALS2_`   = "1"   # harden internal checks → turn latent issues into aborts
)
options(warn = 2)                  # warnings → errors

# Attach default packages + your package (like a normal session)
for (pkg in c("stats","graphics","grDevices","utils","datasets","methods"))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
suppressPackageStartupMessages(library(gflow))

# Parse whole file to keep multi-line expressions intact
ex_file <- "gflow-Ex.R"
exprs   <- parse(file = ex_file, keep.source = TRUE)
lines   <- readLines(ex_file, warn = FALSE)
env     <- new.env(parent = .GlobalEnv)

handle_fail <- function(i, cond) {
  sr  <- attr(exprs[[i]], "srcref")
  msg <- conditionMessage(cond)
  if (!is.null(sr)) {
    start <- sr[[1]]; end <- sr[[3]]
    cat("\n### FAILURE at expression", i, sprintf("(lines %d-%d)", start, end), "###\n")
    cat("Message: ", msg, "\n", sep = "")
    ctx <- lines[pmax(1, start-5):pmin(length(lines), end+5)]
    cat("Context:\n", paste0("  ", ctx), sep = "\n")
  } else {
    cat("\n### FAILURE at expression", i, "(no srcref) ###\nMessage: ", msg, "\n", sep = "")
  }
  dump.frames("last.dump", to.file = TRUE)
  cat("\nWrote dump to: ", file.path(getwd(), "last.dump.rda"), "\n", sep = "")
}

cat(sprintf("Running %d expressions under check-like settings...\n", length(exprs)))
for (i in seq_along(exprs)) {
  tryCatch(
    {
      withCallingHandlers(
        {
            # Make GC more aggressive every so often to surface relocations
            if (i %% 10L == 0L) gctorture2(on = TRUE, inhibit_release = TRUE)
            on.exit(gctorture2(on = FALSE), add = TRUE)

            eval(exprs[[i]], envir = env)
        },
        warning = function(w) {
          if (grepl("stack imbalance", conditionMessage(w), fixed = TRUE)) {
            handle_fail(i, w); stop(w)
          } else invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) { handle_fail(i, e); stop(e) }
  )
}
cat("\nAll expressions completed without stack imbalance.\n")


## strict like R CMD check
Sys.setenv(`_R_CHECK_INTERNALS2_`="1", R_ABORT_ON_PROTECT_ERROR="1")
options(warn = 2)

## attach defaults + your pkg
for (p in c("stats","graphics","grDevices","utils","datasets","methods"))
  suppressPackageStartupMessages(library(p, character.only=TRUE))
suppressPackageStartupMessages(library(gflow))

## parse whole examples file (handles multi-line calls)
ex_file <- "gflow-Ex.R"
exprs   <- parse(file = ex_file, keep.source = TRUE)
lines   <- readLines(ex_file, warn = FALSE)
env     <- new.env(parent = .GlobalEnv)

handle_fail <- function(i, cond){
  msg <- conditionMessage(cond); sr <- attr(exprs[[i]], "srcref")
  if (!is.null(sr)) {
    s <- sr[[1]]; e <- sr[[3]]
    cat("\n### FAILURE at expr", i, sprintf("(lines %d-%d)", s, e), "###\n")
    cat("Message: ", msg, "\n", sep = "")
    cat("Context:\n", paste0("  ", lines[pmax(1, s-5):pmin(length(lines), e+5)]), sep="\n")
  } else cat("\n### FAILURE at expr", i, "(no srcref) ###\nMessage: ", msg, "\n", sep = "")
  dump.frames("last.dump", to.file = TRUE)
  cat("\nWrote dump to: ", file.path(getwd(), "last.dump.rda"), "\n", sep = "")
}

cat(sprintf("Running %d expressions under check-like settings...\n", length(exprs)))
for (i in seq_along(exprs)) {
  tryCatch(
    {
      ## make GC nasty every N expressions
      if (i %% 10L == 0L) gctorture2(step = 1, wait = 0, inhibit_release = TRUE)
      on.exit(gctorture2(step = 0), add = TRUE)

      withCallingHandlers(
        { eval(exprs[[i]], envir = env) },
        warning = function(w) {
          if (grepl("stack imbalance", conditionMessage(w), fixed = TRUE)) {
            handle_fail(i, w); stop(w)
          } else invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) { handle_fail(i, e); stop(e) }
  )
}
cat("\nAll expressions completed with no stack imbalance.\n")



# Read the generated examples file
examples_code <- readLines("gflow-Ex.R")

# Find all function names being tested
# Look for patterns like "### Name: functionname"
name_lines <- grep("^### Name:", examples_code)
function_names <- gsub("### Name: ", "", examples_code[name_lines])

# Run each example block separately
for (i in seq_along(function_names)) {
  cat("\n\nTesting function:", function_names[i], "\n")
  
  # Extract code for this specific example
  start_line <- name_lines[i]
  end_line <- if (i < length(name_lines)) name_lines[i+1] - 1 else length(examples_code)
  
  # Create temporary file with just this example
  example_code <- examples_code[start_line:end_line]
  writeLines(example_code, "temp_example.R")
  
  # Try to run it
  tryCatch({
    source("temp_example.R", echo = TRUE)
    cat("SUCCESS: No crash for", function_names[i], "\n")
  }, error = function(e) {
    cat("ERROR in", function_names[i], ":", e$message, "\n")
  })
}
