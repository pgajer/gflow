# Agent Prompt: Expand The Data Graph Constructions Report

You are working in the `gflow` R package repository:

`/Users/pgajer/current_projects/gflow`

Your task is to expand and update the LaTeX report currently located at:

- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/phate_knn_graph_constructions.tex`
- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/build/phate_knn_graph_constructions.pdf`

Before editing the report, read this project brief carefully:

- `/Users/pgajer/current_projects/gflow/dev/phate-knn-graph-constructions/data_graph_constructions_report_brief.md`

Use the brief as the source of context, assets, and task priorities. The report
began as a PHATE/kNN graph construction note, but it should now become a
self-contained report on data-to-graph constructions for geodesic distance
approximation.

Before making report edits, generate a phased action plan for accomplishing the
main objectives. The plan should separate:

1. report reframing and rename decisions;
2. introduction and problem formulation;
3. graph construction method sections;
4. pruning and MST repair sections;
5. new figures and examples;
6. historical perspective and references;
7. build, QA, and final review.

For each phase, list:

- source files to inspect or edit;
- expected output artifacts;
- validation or build steps;
- assumptions or decisions that need confirmation.

After sharing the plan, implement the phases in order unless the user redirects.

Main objectives:

1. Rename/reframe the report around data graph constructions for geodesic
   distance approximation. PHATE should remain in the report, but as an
   adaptive affinity/diffusion construction whose support motivates comparison
   with sKNN and adaptive-radius graphs.
2. Replace the date stamp under the title with a `Build timestamp` line,
   consistent with other LaTeX reports in Codex projects.
3. Add a strong introduction explaining why recovering geodesic distance
   structure matters for highly structured data with nonlinear geometry,
   topology, folds, branches, and near self-approaches.
4. Formulate the data geodesic geometric reconstruction problem precisely.
5. Add detailed algorithmic sections for fixed-radius and adaptive-radius graph
   constructions.
6. Add detailed sections for all pruning methods used in `gflow`, especially
   local geodesic pruning, with formulas and explanatory figures.
7. Add a detailed section for MST/component-MST connectivity repair, with
   formulas and explanatory figures.
8. Add additional 2D examples such as figure eight, crossing or nearly crossing
   line segments, tree/branching data, nonuniform noisy circle, and quadratic
   graph surfaces.
9. Add references and a historical perspective for kNN, mutual kNN, radius
   graphs, adaptive local scaling, shared nearest-neighbor graphs, PHATE,
   geometric pruning/spanners, and MST repair.
10. Add tables or diagrams that clarify terminology, edge-support rules,
    affinity-versus-length conventions, and nesting/non-nesting relationships.

Important implementation/reporting notes:

- Do not hand-edit generated PDFs. Edit the TeX source and figure-generation
  scripts/assets.
- Preserve mathematical precision. All formulas should be LaTeX formatted.
- Verify figure source paths and rebuild the report after changes.
- The report should read as a human-facing research note, not as a log of
  implementation details.
- `dev/*` is ignored by git in this repository, so if any new dev-report files
  should be version controlled, they must be explicitly force-added after user
  approval.
