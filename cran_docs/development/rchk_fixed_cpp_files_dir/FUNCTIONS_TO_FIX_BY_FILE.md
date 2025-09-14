# Functions Requiring rchk Corrections - Sorted by File

## Summary
- **Total Functions to Fix**: 32 functions
- **Total Files**: 29 files  
- **Already Fixed**: 4 functions (marked with ✅)

## Files Sorted by Number of Functions

### 1. `set_wgraph.cpp` (3 functions)
- [ ] S_find_graph_paths_within_radius
- [ ] S_remove_redundant_edges
- [ ] S_compute_edge_weight_rel_deviations

### 2. `pruning_long_edges.cpp` (2 functions)
- [ ] S_wgraph_prune_long_edges
- [ ] S_find_shortest_alt_path

### 3. `pgmalo.cpp` (2 functions)
- [ ] S_pgmalo
- [ ] S_upgmalo

### Files with 1 Function Each (26 files)

#### ✅ Already Fixed (4 files)
- `agemalo_r.cpp`: ✅ S_agemalo (FIXED)
- `amagelo_r.cpp`: ✅ S_amagelo (FIXED - included in agemalo_r.cpp output)
- `ray_agemalo_r.cpp`: ✅ S_ray_agemalo (FIXED)
- `graph_kernel_smoother_r.cpp`: ✅ S_graph_kernel_smoother (FIXED)

#### Graph Smoothing Functions (8 files)
- [ ] `graph_bw_adaptive_spectral_smoother_r.cpp`: S_graph_bw_adaptive_spectral_smoother
- [ ] `graph_spectral_lowess_r.cpp`: S_graph_spectral_lowess
- [ ] `graph_spectral_lowess_mat_r.cpp`: S_graph_spectral_lowess_mat
- [ ] `graph_spectral_ma_lowess_r.cpp`: S_graph_spectral_ma_lowess
- [ ] `graph_diffusion_smoother.cpp`: S_graph_diffusion_smoother
- [ ] `graph_deg0_lowess_cv.cpp`: S_graph_deg0_lowess_cv
- [ ] `graph_deg0_lowess_buffer_cv_r.cpp`: S_graph_deg0_lowess_buffer_cv
- [ ] `graph_deg0_lowess_cv_mat_r.cpp`: S_graph_deg0_lowess_cv_mat

#### Graph Construction/Utilities (2 files)
- [ ] `graph_gradient_flow.cpp`: S_construct_graph_gradient_flow
- [ ] `graph_utils.cpp`: S_join_graphs

#### Graph Structure Functions (3 files)
- [ ] `hHN_graphs.cpp`: S_create_hHN_graph
- [ ] `mst_completion_graphs_r.cpp`: S_create_mst_completion_graph
- [ ] `nerve_cx_r.cpp`: S_create_nerve_complex

#### Statistical/Analysis Functions (4 files)
- [ ] `density.cpp`: S_estimate_local_density_over_grid
- [ ] `geodesic_stats_r.cpp`: S_compute_geodesic_stats
- [ ] `gflow_basins_r.cpp`: S_create_basin_cx
- [ ] `centered_paths.cpp`: S_ugg_get_path_data

#### Other Functions (4 files)
- [ ] `parameterize_circular_graph_r.cpp`: S_parameterize_circular_graph
- [ ] `pgmalog.cpp`: S_upgmalog
- [ ] `nada_graph_spectral_lowess_r.cpp`: S_nada_graph_spectral_lowess
- [ ] `iknn_graphs.cpp`: create_iknn_graph

## Functions by Category

### Graph Smoothing (8 functions) - High Priority
1. S_graph_bw_adaptive_spectral_smoother
2. S_graph_spectral_lowess
3. S_graph_spectral_lowess_mat
4. S_graph_spectral_ma_lowess
5. S_graph_diffusion_smoother
6. S_graph_deg0_lowess_cv
7. S_graph_deg0_lowess_buffer_cv
8. S_graph_deg0_lowess_cv_mat

### Graph Construction/Utilities (6 functions)
1. S_construct_graph_gradient_flow
2. S_join_graphs
3. S_find_graph_paths_within_radius
4. S_find_shortest_alt_path
5. S_remove_redundant_edges
6. S_compute_edge_weight_rel_deviations

### Graph Structures (5 functions)
1. S_wgraph_prune_long_edges
2. S_create_hHN_graph
3. S_create_mst_completion_graph
4. S_create_nerve_complex
5. create_iknn_graph

### Statistical/Analysis (4 functions)
1. S_estimate_local_density_over_grid
2. S_compute_geodesic_stats
3. S_create_basin_cx
4. S_ugg_get_path_data

### PGMALO Family (3 functions)
1. S_pgmalo
2. S_upgmalo
3. S_upgmalog

### Miscellaneous (2 functions)
1. S_parameterize_circular_graph
2. S_nada_graph_spectral_lowess

## Progress Tracking
- ✅ Fixed: 4 functions (12.5%)
- ❌ Remaining: 28 functions (87.5%)

## Next Steps Priority Order

### High Priority (Multiple issues or critical errors)
1. `set_wgraph.cpp` (3 functions with multiple issues)
2. `pruning_long_edges.cpp` (2 functions with unprotected allocations)
3. `graph_diffusion_smoother.cpp` (lambda capture issues)
4. `graph_gradient_flow.cpp` (unprotected allocations)

### Medium Priority (Single issue per function)
- All graph smoothing functions
- PGMALO family functions
- Graph structure functions

### Lower Priority (Simpler fixes)
- Statistical/analysis functions
- Miscellaneous functions