# Synthetic Data API Consolidation Plan Architecture Re-audit

Re-audit date: 2026-07-28  
Auditor role: independent package-architecture auditor  
Revised handoff:
`split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md`  
Prior audit:
`split_audit/synthetic_data_api_plan_architecture_audit_2026-07-28.md`  
Audit response:
`split_audit/synthetic_data_api_plan_architecture_audit_response_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`  
Audited repository revision: `92a61c086f2fa1fa77223edfb02b74a1be3f1a28`

## Independence and scope

I did not author the revised handoff, audit response, report source, generated
HTML, or source implementations. I used the response as a finding-disposition
ledger, not as the scope of the re-audit. I independently re-read the revised
contracts, re-materialized raw legacy outputs, falsified underspecified
alternatives, rerendered the report, and repeated responsive browser checks.

This remains a design-only phase. There are no fitted estimators, metrics,
selection comparisons, or uncertainty statements, so the measurement,
estimation/selection, and inference layers are not applicable rather than
passed.

## Verdict

**Revise before acceptance.**

The revision resolves most structural findings from the first audit:
recipe/instance identity is now separated; `materialize.synthetic()` replaces
the conflicting `simulate.synthetic()` name; frame randomness belongs to the
materializer; G4/G7 have explicit nullable/heterogeneous representations; the
missing constructor signatures are substantially filled in; provenance is
expanded; and the rendered tables are readable.

Three blockers remain:

1. exact legacy parity is still under-specified at the random-algorithm level;
2. the worked SSRHE examples request a seed policy that the plan does not
   define; and
3. arbitrary truth closures cannot satisfy the plan's ad hoc ID and checksum
   contracts.

Additional major comments concern the meaning of design-normalized “population
truth,” cross-component domain validation, and the cross-platform checksum
scope.

## Disposition of prior findings

| Prior finding | Re-audit disposition | Reason |
|---|---|---|
| DGP-1: 1-D ordering/allocation/outlier semantics | Partially resolved; blocker remains | Counts, sorting, and minimum-outlier behavior are fixed, but the exact truncated-normal, gapped-uniform, and Laplace algorithms/draw order are not fully specified. |
| DGP-2: G4/G7 object representation | Resolved for ordinary maintained cases; edge comment remains | Nullable latent fields, masks, and region dimensions are defined. Degenerate/absent-stratum behavior still needs a validation rule. |
| ART-1: recipe versus instance identity | Resolved | Recipes no longer carry frozen `n`, seed, IDs, or checksums; frozen instances have a no-override materializer. |
| ART-2: checksum contract | Partially resolved; new blocker and major comment | Payload and serialization are defined, but arbitrary closures are incompatible with specification hashing, and random-frame/platform comparison needs tighter scope. |
| ART-3: missing provenance | Resolved | All six material source repositories now have paths and matching revisions. |
| IMP-1: `stats::simulate()` collision | Resolved | The ordinary operation is now `materialize.synthetic()`. |
| IMP-2: frame RNG ownership | Resolved | Constructors are draw-free and the materializer owns the geometry-frame stream. |
| IMP-3: missing constructor contracts | Substantially resolved; major compatibility comment remains | Signatures are present, but some cross-component support/range invariants are absent. |
| REND-1: hidden ledger column | Resolved | All three columns render and remain reachable at the tested widths. |

## Findings by Audit Charter layer

### 1. Data-generating process

#### [Blocker DGP-R1] “Exact” legacy parity still permits multiple RNG algorithms with different outputs

The revised plan correctly adds sorting, fixed interval counts, the
minimum-one-outlier rule, and top-level V3 response draw order
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1522-1560`).
It does not, however, fix the algorithms that turn those draws into values.

The raw SSRHE source is more specific:

- V2 computes `p0`/`p1`, draws uniforms on that probability interval, applies
  `qnorm()`, and sorts
  (`ssrhe_order3_l1_validation_helpers.R:108-111`);
- V3 calls `runif()` first for the left interval and then for the right
  interval before sorting (`:113-114`);
- the Laplace draw is the exact inverse transform
  `-scale * sign(u) * log1p(-2 * abs(u))` with
  `u = runif(n) - 0.5` (`:117-120`).

The proposed signatures
`synthetic.sampling.truncated.normal()` and
`synthetic.response.laplace.outlier()` do not carry an algorithm/version field
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1381-1383`
and `:1457-1463`), and the prose says only to draw truncated-normal or Laplace
values. Rejection sampling is a valid truncated-normal algorithm, and a
difference of exponential variables is a valid Laplace algorithm, but both
consume different RNG sequences and fail the required field-by-field parity.

Independent falsification at the V2 legacy seed produced:

```text
source inverse-CDF V2 x[1]       = 0.2909025801394634
plausible rejection V2 x[1]     = 0.0824320274374329
identical vectors                = FALSE
```

For V3 Laplace errors:

```text
source inverse-transform error[1] = -0.1197540009954813
difference-of-exponentials error[1] = 0.1540832518526744
identical vectors                  = FALSE
```

Both alternatives satisfy the current public distribution-level signatures.
Therefore parity tests would detect a failure, but the plan still does not tell
the implementer which canonical private algorithm must pass.

Add versioned legacy algorithm fields or explicit compatibility rules for:

- V2 inverse-CDF sampling and its `pnorm`/`runif`/`qnorm` order;
- V3 left-interval draw followed by right-interval draw; and
- V3's exact uniform inverse-transform Laplace draw.

The same principle should be applied to G-family paths where bitwise parity is
required, notably disk sampling, clustered draws, and Dirichlet gamma
normalization. A seed policy alone does not define a random algorithm.

#### [Major DGP-R2] Design-normalized truth is not the stated population conditional expectation

The report defines its target as \(m(x)=E[Y\mid X=x]\) and calls the occupation
probability a population truth
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:356-365`
and `:1022-1024`). The proposed occupation recipe then uses
`normalization = "design.maximum"`:

\[
\pi_i = \pi_{\max}
\frac{\operatorname{density}(x_i)^\gamma}
{\max_j\operatorname{density}(x_j)^\gamma},
\]

as specified at `:1562-1571`. The primary Python source confirms that the
denominator is recomputed from each sampled design
(`118_eod_flat2d_gaussian_mixture_examples.py:444-451`).

Thus the probability assigned to the same \(x_i\) can change when the other
sampled design points, sample size, or replicate changes. It is a finite-design
conditional probability \(E[Y_i\mid X_1,\ldots,X_n]\), not a fixed population
function \(E[Y\mid X=x]\). The `"sample.max"` Gaussian-mixture normalization has
the same conceptual qualification.

Either:

- label these explicitly as finite-design normalized benchmark truths and
  record the design-level estimand in `truth.spec`; or
- use a fixed, recipe-level normalizing constant (for example, a validated
  analytic/numerical maximum) when claiming a population conditional
  expectation.

The ownership conclusion need not change, but the scientific claim and
estimand metadata must.

#### [Major DGP-R3] G4 absent-stratum behavior is inconsistent with the heterogeneous-region invariant

For heterogeneous geometry, the plan requires `*.by.region` fields to be keyed
exactly by observed region labels
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:674-679`).
The G4 allocation permits `minimum.counts = c(A = 1L, B = 0L)` and states
`nB = n - max(1L, round(fracA * n))` (`:1507-1513`). Consequently some allowed
inputs have no observed B row while the specified dimension vector still
contains A and B.

The current source is not a clean fallback: `dgp.g4(n = 2, fracA = 0.9)`
errors during matrix assembly because `nB = 0`. The new contract should make a
deliberate decision rather than inherit an accidental error:

- require at least one row per declared G4 stratum and reject incompatible
  `n`/`fracA` early; or
- allow declared-but-unobserved strata and change the “keyed exactly by
  observed labels” invariant accordingly.

The same declared-versus-observed distinction should be stated for G7 when
`zero.fraction` is 0 or 1.

### 2. Measurement

Not applicable. No measured score or diagnostic is reported.

### 3. Estimation and selection fairness

Not applicable. No estimator comparison or tuning selection is reported.

### 4. Statistical inference

Not applicable. No inferential comparison is reported.

### 5. Artifacts and provenance

#### [Blocker ART-R1] Arbitrary truth closures cannot receive a stable ad hoc ID or canonical checksum

The plan explicitly supports arbitrary closures for interactive and
compatibility use
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1420-1422`
and `:1592-1596`). For a specification without `recipe.id`,
`materialize.synthetic()` derives `dataset.id` from a specification SHA-256
(`:1306-1315`). The checksum payload also includes `truth.spec`, but the
canonicalization rules prohibit closures (`:1693-1719`).

No specification-hash payload is defined separately, and no unhashable or
non-checksummed dataset state exists. This leaves several contradictions:

- hashing the closure/environment is not canonical;
- omitting the closure can give two different functions the same derived
  dataset ID when their labels and parameters match;
- retaining the closure violates the checksum payload contract; and
- a caller-supplied `recipe.id` avoids the ad hoc ID hash but does not make
  `synthetic.dataset.checksum()` valid.

Choose and document one policy. Sound options include:

- remove arbitrary closures from materializable public specifications;
- require a caller-supplied versioned function ID plus a registered evaluator;
- require an explicit caller-supplied source/fingerprint and define its
  canonicalization; or
- mark such datasets non-replayable and non-checksummable with a distinct ID
  and provenance contract.

Until this is resolved, the ordinary public API cannot implement every
specification its constructor accepts.

#### [Major ART-R2] Random-frame cross-platform comparison and “same-platform” scope remain incomplete

The report allows predictor matrices from random-frame datasets to differ by an
ambient orthogonal transformation and compares their pairwise squared-distance
matrices (`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1741-1745`).
The materializer also records the realized frame, which belongs in
`parameters`, and `parameters` is in the checksum/comparison payload
(`:588-593`, `:1697-1703`). The comparator does not say whether the recorded
frame itself is compared elementwise, transformed by the same alignment, or
excluded from scientific parity. Elementwise comparison can reject two
scientifically equivalent frames even when the predictor invariant passes.

In addition, the exact hash's stated scope records R major/minor version, RNG
kind/version, registry version, and BLAS/LAPACK fingerprint (`:1727-1731`) but
does not explicitly record OS/architecture or the math-library/build context.
Truth paths using `sin`, `cos`, `dnorm`, `qnorm`, and related functions can
differ in their final bits without involving BLAS.

Specify the realized-frame comparison invariant and expand the exact-hash
environment fingerprint. The handoff correctly admits that empirical
cross-platform calibration remains future work; the plan should not call this
portion fully resolved until these choices are explicit.

### 6. Estimator and implementation correctness

#### [Blocker IMP-R1] The worked SSRHE specifications name an undefined legacy seed policy

Both worked examples set:

```r
compatibility = list(
  recipe = "ssrhe.<family>.v1",
  seed.policy = "legacy.ssrhe"
)
```

at `synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:899-903`
and `:940-947`. The RNG contract says the materializer dispatches on
`compatibility$seed.policy`, but its implemented names are
`ssrhe.1d.v1`, `ssrhe.flat.v1`, and `ssrhe.quadform.v1`
(`:1619-1633`). It defines no `"legacy.ssrhe"` policy.

As written, the report's canonical flat and quadform examples cannot select
their declared RNG implementation. Change their `seed.policy` values to
`"ssrhe.flat.v1"` and `"ssrhe.quadform.v1"`, respectively, or define and
document how `"legacy.ssrhe"` dispatches through `compatibility$recipe`.

#### [Major IMP-R2] Geometry/sampling support compatibility is not validated

Sphere cap, helix, and torus-patch constructors contain footprint or parameter
ranges (`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1325-1351`),
while the sampling components independently contain disk, interval, or
rectangle bounds (`:1373-1391`). The report says `synthetic.spec()` rejects
missing or incorrectly classed components and registry resolution checks
component dimensions (`:1293-1297`, `:1686-1689`), but it does not define
support compatibility.

For example, a sphere-cap geometry can declare one footprint while its disk
sampler draws outside that footprint—or even outside the sphere radius. Helix
and torus recipes can similarly carry geometry ranges that disagree with their
sampling support. This creates two sources of truth for the parameter domain.

Define whether geometry ranges are:

- hard admissible domains within which sampling support must lie;
- defaults copied into a sampler; or
- redundant fields that should be removed from one component.

Then require `synthetic.spec()` to validate the corresponding cross-component
invariant.

#### [Nonblocking IMP-R3] `uniform.box` documents an option absent from its signature

The signature permits
`order = c("draw", "ascending.first.coordinate")`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1373-1375`),
but the prose says `"ascending"` sorts a one-dimensional uniform-box sample
(`:1522-1526`). Add `"ascending"` to the signature or limit that sentence to
the interval and truncated-normal samplers.

### 7. Rendering fidelity

The prior rendering blocker is resolved.

An independent rerender produced identical bytes after replacing only the two
dynamic Eastern build timestamps. The normalized SHA-256 for both original and
re-rendered HTML was:

```text
a6f869502d171505ee27fdda3b2901b6862756aa0eb4c95a901d41fdb1772df4
```

After waiting for MathJax processing to finish, Chromium at 1,440, 768, and
390 pixels reported document `scrollWidth == clientWidth`, four ownership
cards, no ownership code block, and three reachable headers in every table.
The absorption ledger visibly contains all three columns at desktop and mobile
widths.

A transient horizontal overflow exists while MathJax is in its processing
phase, but it disappears when rendering completes and is not treated as a
report defect.

## End-to-end reproduction and falsification evidence

Raw values independently reproduced from the current SSRHE helper:

```text
V2: n=10, shape S01, replicate=1
seed     = 274202
x[1]     = 0.2909025801394634
truth[1] = 1.0000000000000000
y[1]     = 0.8710440501086796

V3: n=10, shape S01, replicate=1
seed       = 274302
x[1]       = 0.0426684002205729
truth[1]   = 0.0000009467618894
y[1]       = 0.0105267645408216
left count = 4

flat: d=3, n=180, replicate=1
seed     = 287801
U[1,1]   = 0.1601577638648450
truth[1] = 0.9610513819826418
y[1]     = 0.9762859385318423
```

Falsification checks included:

- substituting distributionally valid rejection and exponential-difference
  samplers under the revised signatures, which failed exact parity;
- exercising G4 at `n = 2, fracA = 0.9`, which exposed the absent-B edge;
- checking G7 at `zero.fraction = 0` and `1`, which exposed the need to
  distinguish declared from observed regions;
- tracing arbitrary closures through derived ID and checksum requirements;
- verifying all six recorded repository revisions;
- independently rerendering and byte-normalizing the HTML; and
- checking the completed DOM and ledger screenshots at three viewport widths.

There is no numerical method-ranking headline, quality flag, or dependence-
sensitive interval to stratify or recompute in this design-only phase.

## Artifact validation

The handoff's artifact facts were reproduced:

```text
Rmd
  lines: 2,056
  bytes: 80,608
  SHA-256: 19e04c598523ead9e3a9e61f1cfd2169e904d48394b0a960d8479f1c177ada24

HTML
  lines: 4,017
  bytes: 1,342,153
  SHA-256: a7b0332e0b71f527d098e6ccb8ec2f6ee90b6dca1b2193f0af40e0ac7de7c17a
```

The six recorded source revisions match current `HEAD`. The revised handoff
retains the required facts-and-admissions structure, identifies canonical and
generated files, reports validation not performed, and leaves the verdict to
the auditor.

## Validation commands

Principal commands included:

```sh
shasum -a 256 \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html

wc -l -c \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html

Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file = "reaudit-render.html",
  output_dir = <temporary-directory>,
  quiet = TRUE
)'
```

I also used line-numbered source inspection, targeted R evaluation of the
legacy generator definitions, Beautiful Soup byte/structure comparison, and
Playwright/Google Chrome responsive checks after MathJax completion.

No `make check-fast` or `make check` was run. No package source is changed, and
package checks cannot validate an unimplemented architecture contract.

## Required revision and re-audit

Before acceptance:

1. specify the exact legacy random algorithms and internal draw order required
   by DGP-R1;
2. resolve the arbitrary-closure ID/checksum contradiction in ART-R1;
3. correct the SSRHE seed-policy identifiers in IMP-R1;
4. correct or explicitly adjudicate the major DGP/ART/implementation comments;
5. rerender and update the factual handoff and audit response; and
6. return the revised artifacts for another independent disposition.

