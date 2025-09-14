# rchk PROTECT/UNPROTECT Corrections Summary

## Status Report

### ✅ Completed Corrections (4 functions)

The following functions have been fully corrected and are available in the `rchk_fixed_cpp_files_dir`:

1. **S_agemalo** (`agemalo_r.cpp`) - Fixed negative depth and stack imbalance
2. **S_amagelo** (`agemalo_r.cpp`) - Fixed unsupported UNPROTECT(variable)
3. **S_ray_agemalo** (`ray_agemalo_r.cpp`) - Fixed negative depth and protection counting
4. **S_graph_kernel_smoother** (`graph_smoothing_functions.cpp`) - Fixed variable UNPROTECT

### ❌ Functions Still Requiring Correction (28+ functions)

#### Graph Smoothing Functions (8)
- S_graph_bw_adaptive_spectral_smoother (`graph_bw_adaptive_spectral_smoother_r.cpp:162`)
- S_graph_spectral_lowess (`graph_spectral_lowess_r.cpp:187`)
- S_graph_spectral_lowess_mat (`graph_spectral_lowess_mat_r.cpp:226, 228`)
- S_graph_spectral_ma_lowess (`graph_spectral_ma_lowess_r.cpp:158, 160`)
- S_graph_diffusion_smoother (`graph_diffusion_smoother.cpp:3291`)
- S_graph_deg0_lowess_cv (`graph_deg0_lowess_cv.cpp:559`)
- S_graph_deg0_lowess_buffer_cv (`graph_deg0_lowess_buffer_cv_r.cpp:209`)
- S_graph_deg0_lowess_cv_mat (`graph_deg0_lowess_cv_mat_r.cpp:198`)

#### Graph Construction/Utilities (6)
- S_construct_graph_gradient_flow (`graph_gradient_flow.cpp:1518-1534, 1683`)
- S_join_graphs (`graph_utils.cpp:275, 278`)
- S_find_graph_paths_within_radius (`set_wgraph.cpp:490-558`)
- S_find_shortest_alt_path (`pruning_long_edges.cpp:794, 797`)
- S_remove_redundant_edges (`set_wgraph.cpp:1248`)
- S_compute_edge_weight_rel_deviations (`set_wgraph.cpp:1031`)

#### Weighted Graph/MST/Nerve (4)
- S_wgraph_prune_long_edges (`pruning_long_edges.cpp:1307-1337`)
- S_create_hHN_graph (`hHN_graphs.cpp:387-389`)
- S_create_mst_completion_graph (`mst_completion_graphs_r.cpp:216`)
- S_create_nerve_complex (`nerve_cx_r.cpp:58, 65`)

#### Density/Stats (2)
- S_estimate_local_density_over_grid (`density.cpp:463, 528`)
- S_compute_geodesic_stats (`geodesic_stats_r.cpp:268`)

#### Basins/Paths (2)
- S_create_basin_cx (`gflow_basins_r.cpp:692`)
- S_ugg_get_path_data (`centered_paths.cpp:2700`)

#### PGMALO Family (3)
- S_pgmalo (`pgmalo.cpp:1406`)
- S_upgmalo (`pgmalo.cpp:791`)
- S_upgmalog (`pgmalog.cpp:1269, 1272`)

#### Other Functions (3)
- S_parameterize_circular_graph (`parameterize_circular_graph_r.cpp:89`)
- S_nada_graph_spectral_lowess (`nada_graph_spectral_lowess_r.cpp:214`)
- create_iknn_graph (`iknn_graphs.cpp:652, 661`)

## Correction Patterns to Apply

### Pattern 1: Coercion Block with PROTECT_WITH_INDEX
```cpp
// WRONG:
int tprot = 0;
SEXP sx = s_x;
if (TYPEOF(sx) != REALSXP) { sx = PROTECT(Rf_coerceVector(sx, REALSXP)); ++tprot; }
// ... use sx ...
if (tprot) UNPROTECT(tprot);  // ❌ Variable UNPROTECT

// CORRECT:
{
    SEXP sx = s_x;
    PROTECT_INDEX px;
    PROTECT_WITH_INDEX(sx, &px);
    if (TYPEOF(sx) != REALSXP) REPROTECT(sx = Rf_coerceVector(sx, REALSXP), px);
    // Copy to std::vector
    const R_xlen_t nx = XLENGTH(sx);
    x.assign(REAL(sx), REAL(sx) + static_cast<size_t>(nx));
    UNPROTECT(1);  // ✅ Fixed constant
}
```

### Pattern 2: Container-First Result Assembly
```cpp
// WRONG:
int nprot = 0;
SEXP result = PROTECT(Rf_allocVector(VECSXP, n)); ++nprot;
// ... add elements with ++nprot ...
UNPROTECT(nprot);  // ❌ Variable UNPROTECT

// CORRECT:
SEXP result = PROTECT(Rf_allocVector(VECSXP, n));
// Add each element with local PROTECT/UNPROTECT(1)
{
    SEXP elem = PROTECT(Rf_allocVector(REALSXP, m));
    // ... fill elem ...
    SET_VECTOR_ELT(result, i, elem);
    UNPROTECT(1);  // ✅ Local unprotect
}
// Set names
SEXP names = PROTECT(Rf_allocVector(STRSXP, n));
// ... set names ...
Rf_setAttrib(result, R_NamesSymbol, names);
UNPROTECT(2);  // ✅ Fixed: result + names
```

### Pattern 3: Remove Lambda Captures
```cpp
// WRONG:
auto create_vector = [&nprot](const std::vector<double>& vec) {
    SEXP r_vec = PROTECT(Rf_allocVector(REALSXP, vec.size())); 
    nprot++;  // ❌ Modifying captured variable
    // ...
};

// CORRECT: Don't use lambdas, or use them without protection counting
// Instead, inline the code or create separate functions
```

### Pattern 4: Matrix Allocation
```cpp
// CORRECT:
if (!matrix_data.empty()) {
    const R_xlen_t nrow = static_cast<R_xlen_t>(matrix_data.size());
    const R_xlen_t ncol = static_cast<R_xlen_t>(matrix_data[0].size());
    SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    double* p = REAL(mat);
    // Fill in column-major order for R
    for (R_xlen_t i = 0; i < nrow; ++i) {
        for (R_xlen_t j = 0; j < ncol; ++j) {
            p[i + j * nrow] = matrix_data[i][j];
        }
    }
    SET_VECTOR_ELT(result, idx, mat);
    UNPROTECT(1);
} else {
    SET_VECTOR_ELT(result, idx, R_NilValue);
}
```

## Common Issues and Fixes

| Issue | Fix |
|-------|-----|
| `UNPROTECT(variable)` | Use literal constants only |
| `if (tprot) UNPROTECT(tprot)` | Use PROTECT_WITH_INDEX/REPROTECT pattern |
| Lambda captures `&nprot` | Remove lambdas or inline the code |
| Multiple protection counters | Use single pattern with local PROTECT/UNPROTECT(1) |
| Unprotected during allocation | PROTECT immediately after allocation |
| Negative depth | Balance all PROTECT/UNPROTECT pairs |

## Next Steps

1. Apply these patterns to each remaining function
2. Test each correction with rchk
3. Verify no memory leaks with valgrind
4. Run package tests to ensure functionality unchanged

## Files Organization

Place corrected functions in separate files by category:
- `graph_smoothing_functions.cpp` - All smoothing functions
- `graph_construction_functions.cpp` - Graph construction/utilities
- `graph_pruning_functions.cpp` - Weighted graph/MST/nerve functions
- `stats_functions.cpp` - Density/stats/basins functions
- `pgmalo_functions.cpp` - PGMALO family functions

## Testing Command

After corrections, test with:
```r
rhub::rhub_check(platforms = c("rchk"))
```