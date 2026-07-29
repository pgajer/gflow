# Synthetic Data API Consolidation Plan Architecture Third Re-audit

Re-audit date: 2026-07-28  
Auditor role: independent package-architecture auditor  
Handoff under audit:
`split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md`  
Response inspected:
`split_audit/synthetic_data_api_plan_architecture_second_reaudit_response_2026-07-28.md`  
Prior re-audit:
`split_audit/synthetic_data_api_plan_architecture_second_reaudit_2026-07-28.md`  
Canonical source:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd`  
Generated report:
`split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html`

## Independence and scope

I did not author the handoff, response, revised canonical report, generated
HTML, or primary generator implementations. I used the handoff as factual
evidence, not as the audit agenda. The audit independently traced the revised
contracts from data generation through identities, checksums, implementation
requirements, source provenance, and rendering.

This is still an architecture-only phase. There are no estimator fits, measured
scores, selection comparisons, or inferential claims. The measurement,
estimation/selection, and inference layers are not applicable rather than
passed.

## Verdict

**Accepted with nonblocking comments.**

All five findings from the second re-audit are resolved:

1. G2 and G4 now have versioned, operation-level legacy algorithms that match
   the current primary source;
2. empty G7 `zero.parts` is explicitly treated as a breaking exclusion with a
   consumer-audit gate;
3. general truth terminology is scope-neutral;
4. every ordinary dataset ID contains the complete canonical specification
   digest; and
5. `rng.policy` is a required, validated top-level dataset field.

No blocker remains in the architecture plan. Implementation acceptance remains
a separate future gate: the proposed API, fixtures, serializer, comparator, and
compatibility wrappers do not yet exist.

Two nonblocking comments should be carried into implementation:

- validate every policy-derived seed before any draw so legacy integer-offset
  overflow cannot cause a partial materialization; and
- refresh the generated HTML/handoff if the accepted artifact is meant to
  describe the repository's new post-handoff `gflow` revision.

## Disposition of the five second re-audit findings

| Finding | Third re-audit disposition | Evidence |
|---|---|---|
| DGP-RR1: exact G2/G4 algorithms | Resolved | Both algorithms are versioned, specified operation by operation, covered by fixture requirements, and independently reproduce the source. |
| DGP-RR2: empty G7 `zero.parts` | Resolved | It is a documented breaking exclusion with a dedicated error contract, consumer audit, and narrowed parity surface. |
| DGP-RR3: scope-neutral truth wording | Resolved | General descriptions now say truth estimand/specification; population wording remains only where scientifically specific. |
| ART-RR1: content-bound ordinary IDs | Resolved | The full canonical specification digest is present in every ordinary ID and is independently validated as a top-level field. |
| IMP-RR1: required `rng.policy` | Resolved | The object, checksum, ID, comparator, validator, and acceptance contracts now agree. |

## Findings by Audit Charter layer

### 1. Data-generating process

#### G2 exact-frame contract — resolved

The report now gives random-frame constructors a `frame.algorithm` field
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:579-608`).
Legacy G2 fixes:

- algorithm ID `legacy.base.qr.gaussian.v1`;
- explicit or `seed + 500000L` frame-seed resolution;
- column-major Gaussian matrix construction;
- base `qr(..., LAPACK = FALSE)`;
- `qr.Q(..., complete = FALSE)`;
- first-column selection without sign or pivot post-processing; and
- recording of the resolved seed, algorithm ID, and realized frame
  (`:630-652`).

This is faithful to `.dgp.orthonormal.frame()` and `dgp.g2()`
(`geosmooth/R/dgp_library.R:177-182` and `:291-314`).

Independent reconstruction returned:

```text
G2 default frame              identical = TRUE
G2 default predictors         identical = TRUE
G2 default resolved seed      identical = TRUE
G2 explicit frame             identical = TRUE
G2 explicit predictors        identical = TRUE
G2 explicit resolved seed     identical = TRUE
```

Phase 1 records the algorithm in frozen fixtures, Phase 2 requires exact frame
and predictor comparison before invariant checks, and acceptance criterion 9
requires the versioned base-QR algorithm (`:2008-2040` and `:2224-2225`).

#### G4 exact traversal contract — resolved

The stratified sampler signature now includes
`legacy.g4.stratum.sequential.v1`
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1491-1505`).
The report fixes:

1. A's latent uniform draw;
2. A's second ambient Gaussian column;
3. A's third ambient Gaussian column;
4. B's first latent uniform draw;
5. B's second latent uniform draw;
6. A-before-B row, region, and mask assembly;
7. draw-free truth evaluation; and
8. the conditional `seed + 900000L` response reset
   (`:1624-1655`).

That sequence matches `dgp.g4()` and `.dgp.add.gaussian.noise()`
(`geosmooth/R/dgp_library.R:195-200` and `:517-547`).

Independent reconstruction returned:

```text
G4 eta=.02, response sd=.1:
  predictors identical = TRUE
  latent     identical = TRUE
  truth      identical = TRUE
  response   identical = TRUE

G4 eta=0, response sd=0:
  predictors identical = TRUE
  latent     identical = TRUE
  truth      identical = TRUE
  response   identical = TRUE
```

The second case falsifies a common shortcut: zero-variance ambient Gaussian
columns still consume their legacy draws, while a zero response standard
deviation does not reset or draw.

#### Empty G7 structural-zero indices — resolved as a breaking exclusion

The report explicitly identifies why the current
`zero.parts = integer(0)` behavior is scientifically misleading: selected rows
receive a `"zero"` label without being placed on a structural-zero face
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1657-1675`).
It requires:

- a dedicated condition and remediation message;
- exclusion from the maintained parity surface;
- a repository-wide pre-release consumer audit; and
- replacement by ordinary Dirichlet sampling or genuine structural-zero
  indices.

Phase 1, Phase 4, and acceptance criterion 12 consistently carry the exclusion
(`:2010-2019`, `:2070-2085`, and `:2230-2233`).

An independent read-only search under `/Users/pgajer/current_projects` found
ordinary nonempty `zero.parts` uses but no empty-subset call. That does not
replace the required release-time consumer record, but it materially reduces
the known migration risk.

#### Truth-scope terminology — resolved

The general DGP definition, ownership card, current-library description, and
component abstraction now use “truth estimand” or “truth specification”
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:342-376`,
`:393-399`, `:493-500`, and `:518-526`). Remaining “population truth” uses
define the population scope or contrast it with finite-design truth.

### 2. Measurement

Not applicable. No measured score or diagnostic is reported.

### 3. Estimation and selection fairness

Not applicable. No estimator comparison or tuning selection is reported.

### 4. Statistical inference

Not applicable. No inferential comparison is reported.

### 5. Artifacts and provenance

#### Content-bound ordinary identifiers — resolved

Every ordinary ID now has the form

```text
<recipe-label-or-adhoc>__spec<full-specification-SHA256>__n<n>__seed<seed>__rng<rng-policy>
```

(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1363-1377`).
The readable recipe label no longer substitutes for content identity.
`specification.sha256` is a required top-level field, is recomputed from
`synthetic-specification-v1`, belongs to the dataset checksum payload, and is
verified against registry recipes (`:688-753`, `:1379-1393`,
`:1868-1875`, and `:1903-1916`).

Two different specifications with the same label therefore differ in their
digest segment and cannot collide under the stated SHA-256 identity assumption.
Frozen short IDs remain a controlled exception because instance
materialization binds both the specification and content checksums.

#### [Nonblocking ART-RRR1] The stored HTML predates an unrelated `gflow` commit

The handoff reports `gflow` revision
`92a61c086f2fa1fa77223edfb02b74a1be3f1a28`, which was the repository revision
at the report's recorded build time. During this audit, `gflow` was at:

```text
5567e11f4904c50fb5829ae04f322a408ce571f3
Add canonical basin summary and provenance API
```

The new commit does not modify the synthetic-data source files inspected by
this report. The other five source repositories remain at the revisions stated
in the report.

An independent rerender had the same line count, byte count, headings, tables,
and preformatted-block count. After normalizing the two build timestamps and
the dynamic eight-character `gflow` revision, the stored and rerendered HTML
both had SHA-256:

```text
43b5753d934e74d68e742141e0bca9d36a27e14e14ce54d86eec3ce89e42dcea
```

Thus the only non-timestamp difference is expected dynamic provenance, not
source/render drift. If the stored HTML is retained as a timestamped report of
the audited revision, no correction is needed. If it is intended to be the
current rolling report, rerender it and refresh the handoff at the new `gflow`
revision. No further architecture re-audit is required solely for that refresh.

### 6. Estimator and implementation correctness

#### Required RNG policy — resolved

`rng.policy` now appears in the required dataset fields and has explicit
validation for presence, known policy, ID/instance agreement, and compatible
stream structure
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:686-753`).
It also appears consistently in the ordinary ID, instance registry, checksum
payload, comparator, and acceptance criterion 25
(`:1363-1367`, `:1877-1888`, `:1903-1916`, `:1956-1958`, and
`:2253-2254`).

#### [Nonblocking IMP-RRR1] Legacy derived-seed bounds should be explicit

The public materializer currently requires only a scalar integer seed
(`synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd:1336-1354`).
Legacy policies derive frame and response seeds by adding `500000L` and
`900000L` (`:630-652`, `:1647-1650`, and `:1831-1845`).

R's maximum integer is itself a valid `set.seed()` input, but integer addition
with either offset produces `NA_integer_`. The current source then fails:

```text
set.seed(.Machine$integer.max)  = accepted
dgp.g2(seed = .Machine$integer.max)
  error: supplied seed is not a valid integer
dgp.g4(seed = .Machine$integer.max, sigma = .1)
  error: supplied seed is not a valid integer
```

Before implementing the materializer, define policy-specific seed validation:
every derived seed that the selected path will use must remain a nonmissing,
representable `set.seed()` value. Validate all such seeds before the first draw
so a bad offset cannot leave a partially consumed legacy sequence. This does
not affect the frozen recipes or the architecture verdict, but it closes a
public-contract edge case.

### 7. Rendering fidelity

No rendering defect was found.

The canonical source rerendered successfully. Both stored and independently
rendered HTML files had:

```text
lines:              4,316
bytes:          1,359,803
h1 headings:           20
h2 headings:           20
h3 headings:           10
report tables:           3
preformatted blocks:    25
```

The retained 390-pixel screenshot shows all three absorption-ledger columns
with readable wrapping. The handoff records completed-MathJax DOM checks at
1,440, 768, and 390 pixels with `scrollWidth == clientWidth`, four ownership
cards, and all table headers reachable. The deterministic rerender and
screenshot corroborate that evidence.

## Reproduction and falsification evidence

The audit attempted to break the revised plan by:

- independently reconstructing G2's frame and predictors for both default and
  explicit frame seeds;
- independently reconstructing G4's predictors, masked latent coordinates,
  truth, and response under noisy and zero-variance branches;
- searching all current projects for the excluded empty G7 call rather than
  assuming the consumer risk was absent;
- tracing ordinary recipe-label reuse through the full specification-digest ID
  contract;
- tracing `rng.policy` through object, ID, registry, checksum, comparator, and
  acceptance contracts;
- probing maximum-integer legacy seeds to expose offset overflow; and
- independently rerendering and normalizing only genuinely dynamic provenance.

There is no numerical method-ranking headline, quality flag, pooled effect,
selection protocol, or dependence-sensitive interval to recompute in this
design-only phase.

## Artifact validation

Artifacts supplied by the report author:

```text
handoff
  lines: 247
  bytes: 10,087
  SHA-256: 2463e68940a1a944b6e6dfb852053862bda056a1563926963eeeedd346d3faa3

second re-audit response
  lines: 110
  bytes: 5,037
  SHA-256: 1d91db02fdc7d65132b01dc7d7e49e20532e212703590c07260b4fd50790901f

R Markdown
  lines: 2,323
  bytes: 95,983
  SHA-256: 0dd07d34ca515f703ece9799c39207920162072933f3eb20c495791210d10c6a

HTML
  lines: 4,316
  bytes: 1,359,803
  SHA-256: 74ba1ea7eba3667c825199b577666e2f640efb44ec7597b2064d457a029aa192
```

The response and handoff follow the worker–auditor protocol: they provide
dispositions, facts, commands, and limitations without defining the audit scope
or asserting a verdict.

## Validation commands

Principal commands included:

```sh
Rscript -e 'source(
  "/Users/pgajer/current_projects/geosmooth/R/dgp_library.R"
); ...'

rg -n \
  --glob "*.{R,Rmd,qmd,py,md}" \
  "zero[._]parts" \
  /Users/pgajer/current_projects

Rscript -e 'rmarkdown::render(
  "split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd",
  output_file = "third-reaudit-render.html",
  output_dir = <temporary-directory>,
  quiet = TRUE
)'

shasum -a 256 \
  split_audit/synthetic_data_api_plan_auditor_handoff_2026-07-28.md \
  split_audit/synthetic_data_api_plan_architecture_second_reaudit_response_2026-07-28.md \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.Rmd \
  split_audit/synthetic_data_ownership_and_geosmooth_expansion_2026-07-28.html

git diff \
  92a61c086f2fa1fa77223edfb02b74a1be3f1a28..\
5567e11f4904c50fb5829ae04f322a408ce571f3 \
  -- R/synthetic_data_utils.R R/random_sampling.R \
     R/generate_knn_example_data.R R/graph_generators.R
```

I also used line-numbered source inspection, structural HTML inspection,
normalized HTML hash comparison, repository revision checks, and inspection of
the retained mobile ledger screenshot.

No `make check-fast` or `make check` was run. No package source, export,
documentation, or test was changed by this audit, and package checks cannot
validate an unimplemented architecture plan.

## Carry-forward comments

No architecture re-audit is required for the following nonblocking work:

1. specify and test policy-specific safe derived-seed validation before
   implementing the materializer;
2. perform and record the G7 empty-subset consumer audit at the breaking-release
   gate; and
3. decide whether the stored HTML remains a historical report of revision
   `92a61c08` or should be rerendered as a rolling report at the current
   `gflow` revision.

Future implementation, migration, and package-release phases require their own
tests and independent acceptance under the worker–auditor workflow.
