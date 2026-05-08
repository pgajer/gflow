Please work in the gflow repository:

```text
/Users/pgajer/current_projects/gflow
```

Your task is to create a LaTeX report on geodesic-isometry diagnostics for data-derived graphs on quadratic surfaces.

Before implementing anything, read this brief:

```text
/Users/pgajer/current_projects/gflow/dev/geodesic-distance-estimation/notes/geodesic_isometry_diagnostics_report_brief.md
```

Then inspect the key source assets referenced in the brief, especially:

```text
/Users/pgajer/current_projects/geodesic_MDS/manuscript/geodesic_mds_v3.tex
/Users/pgajer/current_projects/geodesic_MDS/gmds_diagnostic_report.tex
/Users/pgajer/current_projects/gflow/R/quadform_geodesics.R
/Users/pgajer/current_projects/gflow/R/graph_geodesic_distances.R
/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/phate_knn_graph_constructions.tex
```

The report should be a self-contained mathematical and practical specification of diagnostics used to evaluate whether a data-derived weighted graph \(G(X)\) recovers the surface geodesic geometry of a sampled quadratic surface \(\Gamma\). It should explain the background, connect the diagnostics to Geodesic MDS, define each diagnostic precisely with LaTeX formulas, and include examples in the style of the existing graph-construction report.

Main diagnostics to define:

1. Scale-calibrated relative geodesic stress comparing \(D_G\) and \(D_\Gamma\).
2. Signed residuals, signed bias, and shortcut fraction.
3. Relative residual tail quantiles.
4. Distance-band / scale-regime diagnostics.
5. Path-level diagnostics comparing graph shortest paths to reference surface geodesic paths.
6. Normal-displacement visualization/loss for 2D quadratic surfaces in \(\mathbb{R}^3\).

Please first produce an action plan with phases and wait for review before editing files. In that plan, make explicit whether you intend to:

- create only the core mathematical report first, or also generate examples;
- implement example-generation scripts now or leave them as a later phase;
- include normal-displacement visualization as formulas only or as an implemented example;
- reuse existing benchmark outputs or generate small new example data.

Preferred output location:

```text
/Users/pgajer/current_projects/gflow/dev/data-geodesic-reconstruction/geodesic-isometry-diagnostics-report/
```

Expected deliverables after the plan is approved:

```text
geodesic_isometry_diagnostics.tex
geodesic_isometry_diagnostics.bib
build/geodesic_isometry_diagnostics.pdf
figures/
scripts/
```

Use the existing LaTeX report style where appropriate. The report should read like a research-methods document, not like a checklist. All formulas should be LaTeX formatted, all terms should be defined, and each diagnostic should include an interpretation paragraph explaining what high/low/positive/negative values mean.

Important constraints:

- Do not edit generated PDFs or LaTeX build intermediates directly.
- Keep source files and figure-generation scripts as the source of truth.
- If you modify any gflow package code, follow package hygiene and run relevant tests; however, the preferred first pass should be a report-only task unless implementation is clearly necessary.
- Do not commit or push unless explicitly asked.

