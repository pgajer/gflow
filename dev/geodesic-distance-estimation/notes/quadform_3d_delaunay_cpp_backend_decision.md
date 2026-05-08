# 3D Quadform Delaunay C++ Backend Decision Note

This note records the backend spike for porting the validated R/Qhull-backed
3D quadratic-hypersurface reference oracle toward a package-owned native
implementation.

## Current Behavioral Reference

The current reference implementation is:

```r
quadform.delaunay.geodesic.distances()
```

It constructs an approximate epsilon-net in a 3D parameter domain, forces the
sample points into the reference vertex set, uses `geometry::delaunayn()` to
build a Delaunay tessellation, extracts the one-skeleton, filters long edges,
weights retained edges by exact quadratic-hypersurface segment lengths, and
computes shortest-path distances between sample vertices.

The stress-test lane supports the following reference policy:

- provisional reference default: `edge.length.factor = 4`;
- no-filter control: `edge.length.factor = Inf`;
- hardest current cases: high-anisotropy mixed ball/cube cases, especially
  `cube_mixed_k0_c144`;
- target for native parity: R/Qhull behavior, not closed-form exact geodesic
  truth.

## Local Backend Inspection

Installed `geometry` package:

- version: `0.5.2`;
- license: `GPL (>= 3)`;
- imports/linking: `Rcpp`, `RcppProgress`, `magic`, `lpSolve`, `linprog`;
- Delaunay path: `geometry::delaunayn()` calls `C_delaunayn()`;
- `geometry` source vendors Qhull reentrant-style source files such as
  `libqhull_r.c`, `geom_r.c`, `poly_r.c`, `merge_r.c`, `qset_r.c`, and related
  headers.

The installed local system did not expose a convenient `pkg-config` Qhull
development target such as `qhull_r`, and no Homebrew Qhull headers/libraries
were found in the local search. Depending on a system Qhull installation would
therefore add configure complexity and would likely be less reproducible than
vendoring a known Qhull source snapshot.

The `geometry` package documents that its Qhull-derived files are covered by
the Qhull license and distributes the Qhull `COPYING.txt`. The `gflow` package
is GPL (>= 3), so vendoring compatible Qhull sources with copyright/license
notices preserved appears viable, subject to final license review.

## Backend Options

### Option A: Keep Calling `geometry::delaunayn()`

This is the current state.

Pros:

- no native-code port risk;
- already validated by stress tests;
- CRAN package dependency is straightforward as `Suggests`.

Cons:

- not a gflow-owned C++ implementation;
- the main oracle still depends on an optional R package at runtime;
- harder to add low-level diagnostics, reuse allocations, or later optimize
  reference graph construction;
- not suitable if the goal is a self-contained C++ baseline.

Recommendation: keep this as the behavioral reference and fallback while the
native backend is developed.

### Option B: Call `geometry`'s Native Symbol From gflow

`geometry` registers `C_delaunayn()` internally. In principle gflow could call
into the installed package's native routine.

Pros:

- closer to the current validated implementation;
- avoids vendoring Qhull immediately.

Cons:

- still not a gflow-owned implementation;
- fragile cross-package native-symbol dependency;
- awkward lifecycle/error handling;
- poor long-term fit for CRAN-facing package internals.

Recommendation: do not use this as the production direction.

### Option C: Depend On A System Qhull Library

Use system headers/libs for reentrant Qhull.

Pros:

- avoids vendoring Qhull source;
- can use the canonical C API directly.

Cons:

- local inspection found no convenient system Qhull dev target;
- would require configure logic and platform-specific handling;
- CRAN/macOS/Windows portability risk;
- system Qhull versions may differ from the R/Qhull behavioral reference.

Recommendation: not the first choice for gflow.

### Option D: Vendor A Minimal Reentrant Qhull Source Subset

Vendor the same style of reentrant Qhull source used by `geometry`, preserve
the Qhull license notices, and implement a gflow-owned wrapper that extracts
3D Delaunay tetrahedra/edges.

Pros:

- closest path to R/Qhull parity;
- reproducible across platforms;
- no system dependency;
- diagnostics and filtering can be implemented fully inside gflow;
- consistent with gflow already carrying native third-party code such as ANN.

Cons:

- adds bundled C source and license-maintenance responsibility;
- requires careful symbol isolation to avoid conflicts;
- Qhull error and output handling must be adapted safely;
- exact edge ordering/tie behavior may still differ unless options and input
  handling match the R/Qhull reference closely.

Recommendation: preferred implementation path.

## Recommended Native Implementation Path

1. Keep the R/Qhull oracle as the default/fallback while native work begins.

2. Vendor a minimal Qhull reentrant source subset under a clearly isolated path,
   for example:

   ```text
   src/qhull/
   ```

   Preserve Qhull copyright/license files and note the source snapshot.

3. Add a low-level internal native routine:

   ```r
   rcpp_quadform_delaunay_edges_3d(X, qhull.options = "Qt Qbb Qc")
   ```

   It should return the 1-based undirected Delaunay one-skeleton edge matrix,
   plus optional tetrahedron count and Qhull diagnostics.

4. First validate only edge extraction against `geometry::delaunayn()` on small
   fixtures:

   - easy ball interior;
   - high-anisotropy ball mixed;
   - high-anisotropy cube mixed;
   - degenerate/small cases where Qhull behavior is known.

5. Move filtering and diagnostics into C++ only after edge extraction is stable:

   - unfiltered edge count;
   - retained edge count/fraction;
   - attempted filter factors;
   - component count per attempt;
   - relaxation flag;
   - relaxed-to-`Inf` flag.

6. Reuse the existing C++ exact segment-length kernel for quadratic edge
   weights.

7. Add an explicit backend argument:

   ```r
   backend = c("qhull.r", "cpp")
   ```

   Keep `backend = "qhull.r"` as the default until native parity passes the
   fixture suite and stress-test-derived tolerances.

## Regression Fixture Inputs

Use the fixture design in:

```text
dev/geodesic-distance-estimation/notes/quadform_3d_delaunay_oracle_cpp_regression_fixtures.md
```

The native backend should be compared to R/Qhull using:

- exact or near-exact edge-set parity where Qhull tie behavior permits;
- `summarize.isometry.deviation(D_cpp, D_r_qhull, scale = TRUE)`;
- max absolute distance error;
- max relative distance error for nonzero distances;
- filter/retention/component diagnostics.

## Decision

Proceed with **Option D: vendor a minimal reentrant Qhull subset** as the
primary C++ implementation route, while retaining `geometry::delaunayn()` as the
validated R/Qhull behavioral reference and fallback.

The next coding phase should start with the low-level C++ Delaunay edge
extractor only. Do not move the whole oracle at once.
