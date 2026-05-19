# HTML Report Style Guide

This guide records the recommended structure for research HTML reports in the
`dev/` workspace. The goal is to make reports read like compact methods/results
documents rather than execution notebooks.

For build metadata conventions shared by HTML/PDF report projects, especially
LaTeX build timestamps, see
`/Users/pgajer/.codex/notes/report_build_conventions.md`.

## Core Principle

The main report should tell the scientific story. Code, large tables, and raw
diagnostics belong in an appendix unless they are the object being discussed.

Recommended main-text rhythm:

1. State the question.
2. Define the diagnostic or estimator, with formulas where useful.
3. Show one figure.
4. Interpret what the figure says.
5. Move to the next diagnostic.

## Recommended Document Structure

Use this order for most experimental reports:

1. **Purpose**
   - Explain the scientific or algorithmic question.
   - Define the objects being tested.
   - State the main questions explicitly.

2. **Methods/Diagnostic Sections**
   - One section per figure or tightly related figure group.
   - Each section should have:
     - a short motivation paragraph;
     - any necessary formula;
     - the figure;
     - an interpretation paragraph.

3. **Results Summary And Discussion**
   - Answer each stated question directly.
   - Separate positive results from failure modes.
   - Explain what the results imply for the next implementation step.
   - Avoid checklist language; write this like a methods-paper discussion.

4. **What We Learned**
   - Include a concise synthesis section for exploratory or calibration reports.
   - Ground the claims in the rendered tables/figures, not in developer notes.
   - State which hypotheses survived, which failed, and which remain ambiguous.
   - When possible, include a small computed summary table before the prose.

5. **Appendix**
   - Run metadata.
   - Summary tables.
   - Selected raw diagnostic tables.
   - Reproducibility commands.
   - Code fragments needed to reproduce or audit the report.

## Figure Standards

Figures should be readable without zooming.

- Use wide figures for multi-panel diagnostics.
- Avoid facet grids with too many rows and columns in the main text.
- Prefer wrapped labels over rotated labels when possible.
- If labels are still cramped, flip the axes.
- Use log scales when distributions or errors differ by orders of magnitude.
- Put legends at the bottom unless that makes the figure too tall.
- Make strip labels short and interpretable.
- Ensure y-axis labels, facet strips, legends, and long x-axis labels are not cut off.

Recommended R chunk defaults:

```r
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10.5,
  fig.height = 6.2,
  fig.align = "center",
  out.width = "100%",
  dpi = 140,
  dev = "ragg_png"
)
```

Increase `fig.height` for dense facet grids.

## Figure Section Template

Use this pattern for each figure in the main body:

```markdown
## Diagnostic Name

This diagnostic measures ... The quantity plotted is

\[
Q = ...
\]

```{r diagnostic-name, fig.height=7}
# hidden plotting code
```

The figure shows ... This means ... The relevant failure mode is ...
```

## Tables And Code

Do not show code chunks in the main body unless the report is specifically a
software tutorial. In analysis reports:

- set `echo = FALSE` globally;
- keep all substantial table output in the appendix;
- show only small, curated tables in the main body if they are required for
  interpretation;
- prefer CSV/RDS artifacts for large outputs and link to them from the appendix.

## Discussion Style

The final discussion should not read like development notes. It should:

- answer the original questions directly;
- explain whether the evidence is positive, negative, or mixed;
- distinguish mathematical failure from implementation or sampling artifacts;
- identify the next methodological step;
- describe residual uncertainty plainly.

Avoid phrases like "things to look for" in the final section. Prefer direct
statements such as "The experiment does not support..." or "The residual
diagnostics suggest..."

## Visual QA Checklist

Before treating a report as ready:

1. Render the final self-contained HTML.
2. Render a temporary non-self-contained HTML so the figure PNGs can be inspected.
3. Check every figure for:
   - cut-off axis titles;
   - illegible facet labels;
   - overlapping legends;
   - too-small text;
   - panels that are too compressed to read.
4. Confirm code and large tables start only in the appendix.
5. Confirm the discussion answers the questions posed in the purpose section.

For local visual inspection, a temporary non-self-contained render is useful:

```r
rmarkdown::render(
  "path/to/report.Rmd",
  output_dir = tempfile(),
  output_format = rmarkdown::html_document(self_contained = FALSE)
)
```
