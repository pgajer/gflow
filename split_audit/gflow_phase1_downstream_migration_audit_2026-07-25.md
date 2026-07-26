# gflow Phase 1 Downstream Migration Audit

## Scope

This audit covers qualified downstream calls to the 19 generic graph
compatibility functions removed from `gflow`. Every replacement is a
namespace-only change from `gflow::<name>` to `dgraphs::<name>`; function names
and arguments are unchanged.

## Migrated tracked repositories

All matching calls in clean tracked files were migrated, committed, and pushed:

| Repository | Commit | Branch |
|---|---|---|
| `gflow` | `8e861103` | `main` |
| `gflowx` | `1e070f9` | `main` |
| `gflowui` | `5e57741` | `main` |
| `CT_clearance` | `f6e488f` | `main` |
| `geodesic_data_geometry` | `d0418eb` | `main` |
| `geodesicMDS` | `d025c3a` | `main` |
| `virgo2/gene_coabundance_modules` | `7d80811` | `main` |
| `graph_modularity` | `b85ed50` | `bfs-collision-harness` |
| `ZB` | `cda073d` | `main` |

The repository-wide sweep found no remaining qualified calls in clean tracked
files under `/Users/pgajer/current_projects`.

## References intentionally not edited

The following files were not changed because they are untracked, ignored,
inside a non-Git project, or part of an active secondary worktree. Overwriting
them would mix the migration with unrelated work that cannot be committed
safely.

### Non-Git projects

- `AGP/R/11_asv_hv_k_gcv_sweep.R`
- `AGP/R/17_evenness_endpoints_k7.R`
- `gflow_examples/_coordination/20260322_generalization_coordinator/graph_api_compliance_amendment_2026-03-23.md`
- `gflow_examples/_coordination/20260322_generalization_coordinator/next_round_shared_worker_note_2026-03-22.md`
- `gflow_examples/_coordination/20260322_generalization_coordinator/next_round_worker_prompt_nhanes_2026-03-22.md`
- `gflow_examples/_coordination/20260322_generalization_coordinator/next_round_worker_prompt_nyc_taxi_2026-03-22.md`
- `gflow_examples/_coordination/20260322_generalization_coordinator/next_round_worker_prompt_retina_2026-03-22.md`
- `gflow_examples/_coordination/20260322_generalization_coordinator/shared_contract_note_2026-03-22.md`
- `gflow_examples/nhanes_metabolic_states/scripts/nhanes_workflow.R`
- `gflow_examples/retina_cell_cycle/scripts/archive/retina_mutual_knn_side_project_2026-03-23.R`
- `gflow_examples/retina_cell_cycle/scripts/retina_workflow.R`

### Untracked, ignored, or already modified files

- `ZB/l_iners_transition_architecture/analysis/asv/run_hmp_u01_iknn_3x3_2026-03-21.R`
- `ZB/l_iners_transition_architecture/analysis/asv/run_valencia_13k_iknn_3x3_2026-03-21.R`
- `ZB/l_iners_transition_architecture/analysis/vaginal_microbiome/iknn_preprocess_utils.R`
- `ZB/l_iners_transition_architecture/docs/iknn_3x3_pipeline_contract_2026-04-09.md`
- `cell_cycle/R/run_hvg_graph_layout_sweep.R`
- `cell_cycle/R/run_notch_cellcycle_condexp_gcv_k_sweep.R`
- `cell_cycle/R/run_targeted_gene_set_graph_layout_sweep.R`
- `cell_cycle/share/yang_k_selection_handoff/run_gene_set_condexp_gcv_k_sweep.R`
- `grip/dev/quadforms/quadform_gflow_graph_constructor_notes.md`
- `gut_microbiome/scripts/run_phase1_iknn_1x1_absorb_depth3.R`
- `symptoms/R/11_asv_hv_k_gcv_sweep.R`
- `symptoms/R/17_fit_vag_odor_gcv_over_asv_graphs.R`
- `symptoms/R/19_generate_vag_odor_cont_html_and_driver.R`
- `symptoms/R/20_evenness_endpoints_all_asv_k05.R`
- `symptoms/R/29_asv_exact_iknn_utils.R`
- `symptoms/results/asv_exact_iknn_preprocess_2026-03-19/scripts/29_asv_exact_iknn_utils.R`

### Active secondary worktree

- `gflow-w2/vignettes/retina_single_cell_workflow_vignette.Rmd`

## Required action before rerunning an unresolved script

Replace each removed `gflow::` graph call with the identically named
`dgraphs::` call. Do not restore a compatibility wrapper in `gflow`.
