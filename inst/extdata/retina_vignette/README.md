# Retinal vignette assets

`manifest.tsv` records each table or figure, its byte count, MD5 checksum and
original generation time. `source_path` uses a named source root rather than a
machine-specific location: `cell_cycle/` identifies the cell-cycle analysis
inputs, and `retina_example/` identifies the retinal example inputs. Paths
within each root, including run identifiers, are preserved. Multiple inputs
are separated by `; `.

To rebuild from the upstream results, set `GFLOW_RETINA_CELL_CYCLE_DIR` and
`GFLOW_RETINA_EXAMPLE_DIR` to those source directories, then run
`Rscript inst/scripts/build_retina_vignette_assets.R` from the package root.
`GFLOW_RETINA_RUN_ID` selects the example run. The upstream inputs are not
bundled; the tables and figures used by the vignette are included here.
