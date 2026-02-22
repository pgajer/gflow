# gflow development

* Phase 1 legacy-1D decoupling pass completed:
  - Internal consumers that previously depended on legacy `magelo`/`mabilo`-style
    1D smoothers now use spline-based replacements.
  - `magelo.with.external.BB()` now uses a spline backend while preserving its
    legacy output structure for compatibility.
  - Added consistent deprecation warnings across legacy exported 1D regression
    entry points (`magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`,
    `magelog`, `fit.pwlm*`, and `get.magelo.MAB`).
* Package-level docs now mark legacy 1D model-averaging APIs as experimental and
  planned for extraction to a separate package.
* Phase 3 extraction started with new `malo` package:
  - `fit.pwlm*` and the 1D `get.*MAB*` benchmarking/model-comparison families now
    have native implementations in `malo`.
  - gflow legacy entry points for these families now delegate to `malo` when
    available (with in-package fallback retained for compatibility).
* Hard removal of legacy 1D `get.*MAB*` and related helpers from `gflow`:
  - `bias_utils` compatibility wrappers were removed from `gflow`.
  - These APIs are now available only from `malo`.
* Hard migration of remaining legacy 1D exported APIs from `gflow` to `malo`:
  - `magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`, `magelog`,
    `magelo.with.external.BB`, `generate.dirichlet.weights`, and `fit.pwlm*`
    are now thin forwarders to `malo` with no in-package fallback.
* Native/header cleanup after 1D migration:
  - Removed legacy 1D `src` implementations and symbol registrations from
    `gflow`.
  - Removed obsolete legacy 1D headers from `inst/include/gflow` and pruned
    stale 1D declarations from `msr2.h`.
* Migration operational cleanup:
  - Updated collaborator installation docs to require installing `malo`
    before using legacy 1D forwarders in `gflow`.
  - CI workflow now installs `malo` explicitly in rdgraph test jobs.
  - Forwarder-focused tests now skip cleanly when `malo` is unavailable.
  - Removed unused legacy 1D header artifacts from
    `inst/include/gflow` (`lm.h`, `ray_agemalo.hpp`, `uggmalog.hpp`).
  - Removed legacy `maelog` implementation and native symbols from `gflow`
    (R API, native sources, headers, and Rd docs).
  - Hard-cut removal of all remaining 1D forwarders from `gflow`:
    `magelo`, `amagelo`, `mabilo`, `mabilo.plus`, `mabilog`, `magelog`,
    `magelo.with.external.BB`, `generate.dirichlet.weights`, `fit.pwlm`,
    and `fit.pwlm.optimal` are no longer exported by `gflow`.

# gflow 0.1.0 — 2025-09-21

* Initial CRAN release.
* Implements geometric tools for intrinsic-structure modeling, Morse–Smale regression, and related utilities.
