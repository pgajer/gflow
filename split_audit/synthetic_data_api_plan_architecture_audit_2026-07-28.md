# Synthetic Data API Consolidation Plan Architecture Audit

Audit date: 2026-07-28  
Auditor role: independent package-architecture auditor  
Handoff:
`split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`  
Audited repository revision: `92a61c086f2fa1fa77223edfb02b74a1be3f1a28`

## Audit independence and scope

I did not author the handoff, report source, generated HTML, or source
implementations under audit. I treated the handoff as factual evidence rather
than as the audit agenda. The audit followed the data-outward charter in
`worker_auditor_workflow.md`.

This is a design-only phase. There are no benchmark results, fitted estimators,
selection comparisons, or uncertainty statements. The measurement, estimation,
selection, and statistical-inference layers are therefore not applicable; they
are not marked as passed. The operative layers are the data-generating-process
contract, artifacts and provenance, implementation design, and rendering.

## Verdict

**Revise before acceptance.**

The ownership conclusion and the quadform generalization are directionally
sound, and the generated HTML is synchronized with the canonical R Markdown.
However, the plan cannot yet deliver its own exact-parity and homogeneous-
contract acceptance criteria. The blocking issues are:

1. legacy one-dimensional sampling and noise semantics are not represented by
   the proposed component signatures;
2. the common dataset contract does not define valid representations for G4
   and G7;
3. `simulate.synthetic()` collides with R's existing `simulate()` S3 generic
   without being designed as a conforming method;
4. random-frame construction is assigned simultaneously to a supposedly
   draw-free constructor and to the materializer's RNG stream;
5. registry identity is ambiguous when registry rows store `n` and `seed` but
   the materializer requires callers to supply them again; and
6. the central absorption ledger is clipped in the rendered report.

## Findings by Audit Charter layer

### 1. Data-generating process

#### [Blocker DGP-1] The proposed samplers cannot reproduce the legacy 1-D row order and fixed-count mechanisms

The report maps `sample.x()` to generic uniform-box, truncated-normal, and
gapped-uniform sampling constructors
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:775-780`)
and gives signatures without an ordering or fixed-count policy
(`:1186-1199`). It then requires exact scientific-field parity
(`:1322-1335`).

The primary source has additional contractual behavior:

- all three variants sort the sampled coordinates
  (`ssrhe_order3_l1_validation_helpers.R:106-114`);
- V3 deterministically assigns `floor(0.42 * n)` observations to the left
  interval and all remaining observations to the right interval, rather than
  drawing interval membership from probabilities (`:113-114`);
- V3 always contaminates at least one row with
  `max(1L, round(0.04 * n))` (`:128-130`).

The report instead specifies only nonnegative interval probabilities and says
that V3 contaminates `round(outlier.fraction * n)` rows
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1231-1241`).
At `n = 10`, the source contaminates one row while the stated rule contaminates
zero. A generic probabilistic gapped sampler also cannot guarantee the source's
fixed left/right counts or sorted row order.

This is not merely a future-test omission: the proposed signatures lack the
information needed to implement parity. Add explicit ordering and allocation
semantics, including a versioned legacy V3 fixed-count rule and the
`max(1L, ...)` contamination rule. Exercise small `n`, multiple sample sizes,
and multiple replicates rather than freezing only one sample size and
replicate.

#### [Blocker DGP-2] The homogeneous dataset contract is undefined for G4 and G7

The required `synthetic_dataset` fields include scalar `intrinsic.dim`,
`codimension`, and `latent`; validation is said to check dimensions and finite
numeric content
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:609-630`).
The report nevertheless promises exact G4 and G7 absorption (`:599-602`) without
defining the exceptions those families require.

In the primary G4 source, intrinsic dimension varies by stratum, `d` is `NA`,
and the latent matrix deliberately contains `NA` in the absent second
coordinate for stratum A (`geosmooth/R/dgp_library.R:493-548`). In G7, `d` is
`NA` and `U` is `NULL` (`:658-725`). Consequently, literal parity conflicts
with the required fields and a blanket finite-content validator.

The revised plan must define one of the following:

- a first-class heterogeneous-latent representation with
  `intrinsic.dim.by.region`, documented missing-coordinate semantics, and
  geometry-specific validation; or
- a narrower common contract in which `latent`, `intrinsic.dim`, and
  `codimension` are explicitly nullable, with invariants stated for each
  geometry family.

It must also say where G7's simplex dimension and face-specific dimensions live.
Without that, Phase 4 cannot implement the promised common object faithfully.

### 2. Measurement

Not applicable. The report contains no measured performance, telemetry, score,
or diagnostic claim. No measurement layer was accepted by this audit.

### 3. Estimation and selection fairness

Not applicable. No estimators are compared or selected in this phase.

### 4. Statistical inference

Not applicable. No uncertainty statement or inferential comparison is made.

### 5. Artifacts and provenance

#### [Blocker ART-1] Registry identity, sample size, and seed have two competing sources of truth

`datasets.csv` stores `n` and `seed`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1287-1315`),
but `simulate.synthetic()` requires both values from the caller (`:1169-1181`).
The helper ledger likewise presents
`simulate.synthetic(synthetic.registry.spec(id), n, seed)` (`:789-790`).
The plan does not say whether caller values must equal the registry row, may
override it, or create a new dataset identity.

This ambiguity can attach a canonical `dataset.id`, registry tag, and checksum
provenance to a noncanonical payload. It also conflicts with G5: its sampling
specification contains `cluster.count` and `observations.per.cluster`
(`:1193-1195`), while the current source fixes `n = K * m` and has no independent
sample-size argument (`geosmooth/R/dgp_library.R:575-580`).

Separate a reusable recipe identity from a frozen dataset instance, or make
registry materialization consume `n` and `seed` only from the frozen row.
Define override behavior, ID derivation, checksum lookup, and the required
relationship between `n`, cluster counts, and graph-vertex counts.

#### [Major ART-2] The checksum contract is not canonical enough for package tests

The plan requires SHA-256 fixture checksums and package-test enforcement
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:629-630`,
`:1314-1315`, `:1481-1482`, and `:1529-1530`) but does not define the payload's
field order, serialization format/version, numeric encoding, platform scope, or
treatment of floating-point differences. This matters particularly for random
orthonormal frames produced by QR decomposition.

Define a versioned canonical checksum payload and serialization algorithm.
Distinguish bitwise same-platform reproducibility from tolerance-based
cross-platform scientific parity, and do not make a platform-sensitive binary
hash an unconditional CRAN test unless the supported scope is explicit.

#### [Nonblocking ART-3] Provenance is accurate for four repositories but incomplete for other source-backed claims

The recorded revisions for `gflow`, `geosmooth`, `dgraphs`, and
`trend_filtering` match their current `HEAD`s. The report also makes source-
backed ownership claims about `grip` and the occupation-density application
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1013-1029`
and `:914-949`) but does not record those repository revisions in its render-
time revision block (`:1571-1578`). The handoff's command log names principally
the geosmooth and SSRHE inspections, not the commands used to establish all
inventory claims (`synthetic_data_api_plan_auditor_handoff_2026-07-28.md:133-153`).

Record every repository revision that materially supports the ownership
inventory, or label the corresponding statements as note-derived rather than
source-verified. Include the omitted source-inspection commands in the revised
handoff.

### 6. Estimator and implementation correctness

#### [Blocker IMP-1] `simulate.synthetic()` is an unresolved S3 method collision

The report correctly notices that `synthetic.quadform()` resembles a method for
a hypothetical `synthetic()` generic
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1082-1090`)
but misses the non-hypothetical collision in its central verb. `stats::simulate`
already is an S3 generic with formals `object, nsim, seed, ...`. The proposed
ordinary function has formals `spec, n, seed, rng.policy, validate`
(`:1169-1175`) and is called directly throughout the report.

The preferred replacement is:

```r
materialize.synthetic(spec, n, seed, ...)
```

“Materialize” is the appropriate verb: the operation turns a reusable
specification into one realized synthetic dataset. There is no planned
`materialize()` S3 generic, so this name preserves the repository's
dot-delimited style without colliding with an existing generic. The function
should be documented explicitly as an ordinary exported function.

If the package instead intends to participate in R's simulation protocol, it
must define a registered method such as `simulate.synthetic_spec()` and call it
through `simulate()`, with semantics and formals compatible with the generic,
including a decision about what `nsim` means. That route is less suitable for
the current one-specification-to-one-dataset contract.

Implementing `simulate.synthetic()` as an ordinary exported function would make
namespace/roxygen behavior and user expectations fragile.

#### [Blocker IMP-2] Frame construction violates the plan's own RNG ownership rule

`synthetic.quadform()` says `frame = "random.orthonormal"` constructs and stores
a deterministic frame from `frame.seed`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:532-558`).
The exact component contract later says constructors make no random draws and
all randomness occurs only in the materializer (`:1150-1155`). The named-stream
policy also reserves a `geometry.frame` stream (`:1256-1269`).

These three statements cannot all hold. The current G2 implementation makes the
distinction concrete: it draws `Q` from `frame.seed`, then resets to the dataset
seed for latent sampling (`geosmooth/R/dgp_library.R:291-307`).

Choose one owner for frame randomness. A coherent design would let the geometry
constructor store only a frame policy or a supplied matrix, let the materializer
draw a frame on the named geometry stream, and let the legacy G2 policy derive
the historical `seed + 500000L` frame state. If frames are instead frozen inside
specifications, remove the geometry-frame stream and the claim that constructors
are draw-free.

#### [Major IMP-3] Several promised public constructors do not have implementable contracts

The proposed public surface includes sphere-cap, helix, torus-patch, stratified,
point-line junction, simplex, graph geometry, occupation-mixture truth, and
graph-signal truth constructors
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1093-1139`).
The section titled “Exact constructor contracts” supplies signatures for none
of those geometry constructors and omits the two latter truth constructors
(`:1150-1229`). Phase 2 says only to implement them using the same protocol
(`:1366-1369`).

This leaves the G1--G7 migration under-specified despite the handoff's claim of
an exact absorption map. Add signatures, validation rules, latent/ambient
dimension semantics, deterministic row-allocation rules, and compatibility
behavior for every constructor needed before its implementation phase begins.

### 7. Rendering fidelity

#### [Blocker REND-1] The absorption ledger's third column is hidden

The table CSS first sets `overflow-x: auto` and then overwrites it with
`overflow: hidden`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:165-175`).
The central helper absorption ledger is the third table (`:764-794`).

An independent Chromium render measured, at a 1,440-pixel viewport:

```text
table client width: 771 px
table scroll width: 1,055 px
computed overflow-x: hidden
```

The screenshot and DOM both show only “Current helper or object” and “Absorbed
by”; the entire “Result after migration” column is outside the clipped region.
At 768- and 390-pixel viewports all three tables have clipped content, and the
page itself develops horizontal overflow.

Remove the overriding `overflow: hidden`, preserve horizontal scrolling, and
consider allowing code tokens to wrap in this table. Repeat desktop, tablet,
and mobile full-page checks after rerendering.

## End-to-end reproduction and falsification checks

I evaluated the relevant raw helper definitions directly. The full SSRHE helper
cannot currently be sourced in this environment because it stops when
`genlasso` is unavailable, so I parsed and evaluated only the generator
definitions. For quadform materialization, I supplied geosmooth's current
vendored `.dgp.quadform.embed`, matching the stale cross-repository dependency
described by the report.

Raw outputs reproduced:

```text
flat d=3, n=180, replicate=1
seed   = 287801
U[1,1] = 0.1601577638648450
truth[1] = 0.9610513819826418
y[1] = 0.9762859385318423

quadform d=3, index.k=2, curvature=1.25, n=180, replicate=1
seed   = 296126
U[1,1] = -0.2430855613201857
X[1,4] = -0.7412328687512033
truth[1] = -1.6809230437965765
y[1] = -1.5469526006255947
```

Falsification checks attempted:

- Replaying V1's RNG stream without `sort()` produced the same values in a
  different row order; the source output was sorted and not identical.
- At `n = 10`, the source V3 rule produced one outlier while the report's
  `round(0.04 * n)` rule produced zero.
- Materializing G4 confirmed `d = NA` and `NA` latent entries with
  `dim.by.region = c(1, 2)`.
- Materializing G7 confirmed `d = NA` and `U = NULL`.
- Inspecting `stats::simulate` confirmed that it is an existing S3 generic with
  formals `object, nsim, seed, ...`.
- Full-page browser checks at 1,440, 768, and 390 pixels falsified the handoff's
  implicit assumption that top-viewport inspection was adequate for the core
  tables.

There was no pooled numerical headline, quality flag, dependence-sensitive
interval, or method ranking to stratify or recompute in this design-only phase.

## Artifact validation

The handoff's file sizes, hashes, repository status, and four recorded
repository revisions were reproduced exactly:

```text
Rmd SHA-256:
c7be675504a49dc0310ceb2f45c800f3d9d97569dfa5cded2e12d7531663cf42

HTML SHA-256:
9dcfcce54c409e5da714dcdfc91124101d73d21b03114b9b1d63af5b83026ec4
```

An independent render to a temporary directory produced the same title,
heading counts, three tables, 25 preformatted blocks, file size, and normalized
visible text. After replacing the two render-time timestamps with a fixed
token, the original and independent HTML files had the same SHA-256:

```text
6af96a0a5e0e7b78194607394b3453c5a5af0472b280f4501faeb4c32ab8de8d
```

Thus the HTML is current with its canonical R Markdown source. The differing
raw rerender hash is explained solely by the intentionally dynamic build
timestamp.

The handoff itself follows the facts-and-admissions structure: it identifies
canonical and generated files, records validation not performed, includes a
limitations section, and does not propose an audit verdict or checklist. Its
use of “exact” for the absorption map and ledger is not supported until the
blocking parity findings above are resolved.

## Validation commands run

Principal commands were:

```sh
git status --short
git rev-parse HEAD
shasum -a 256 \
  split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html

Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file = "audit-render.html",
  output_dir = <temporary-directory>,
  quiet = TRUE
)'
```

I also used:

- line-numbered inspection of the report, handoff, current geosmooth DGP
  library, SSRHE helper, and historical quadform source at `b8df4053`;
- targeted R materialization of flat, quadform, G4, G7, and one-dimensional
  datasets;
- Beautiful Soup comparison of the original and independent HTML;
- Playwright/Chromium DOM and full-page checks at desktop, tablet, and mobile
  widths; and
- repository revision checks for all named current projects.

No `make check-fast` or `make check` was run. This phase changes only untracked
report artifacts and contains no package implementation. Package checks would
not validate the proposed, unimplemented contract.

## Required revision and re-audit

Before acceptance:

1. resolve findings DGP-1, DGP-2, ART-1, IMP-1, IMP-2, IMP-3, and REND-1 in the
   canonical R Markdown;
2. define the checksum scope requested in ART-2;
3. rerender the self-contained HTML;
4. rerun source-to-HTML synchronization and responsive full-page table checks;
5. revise the factual handoff with new hashes, validation evidence, and
   limitations; and
6. return an audit response mapping every finding to a fix, deferral, or
   reasoned disagreement.
