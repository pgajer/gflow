# Example Batch Backlog

This backlog ranks the next likely-safe example-modernization batches for `gflow`.
It is meant to support small, validated commits that can be checkpointed with
[`tools/run_example_batches.sh`](/Users/pgajer/.codex/worktrees/8174/gflow/tools/run_example_batches.sh)
after the source edits are prepared.

## Selection rules

- Prefer exported `print()`, `summary()`, and `plot()` methods over heavy constructors.
- Prefer files whose methods operate on simple list- or data-frame-backed objects.
- Prefer batches that can be demonstrated with tiny in-memory objects.
- Defer files that still need real model fitting, compiled backends, or complex object invariants.

## Ranked backlog

### 1. `R/fassoc0_test_paired.R`

- Status: ready
- Why first:
  - Strong match for the mock-object pattern already used successfully in [`/Users/pgajer/.codex/worktrees/8174/gflow/R/lslope_test_paired.R`](/Users/pgajer/.codex/worktrees/8174/gflow/R/lslope_test_paired.R)
  - `print.assoc0()`, `summary.assoc0()`, `print.summary.assoc0()`, and `plot.assoc0()` all depend on a manageable list shape
  - High documentation value because this is an exported paired-test result class
- Likely batch scope:
  - `print.assoc0()`
  - `summary.assoc0()`
  - `print.summary.assoc0()`
  - `plot.assoc0()`
- Expected risk: low to medium
- Main watchouts:
  - `plot.assoc0(type = "Exy")` needs small but internally consistent `x`, `y`, `xgrid`, `signal.Eyg`, `null.Eyg`, and `Ey`
  - summary/print paths also need null/signal/diff vectors with sane lengths

### 2. `R/fassoc1_test_paired.R`

- Status: ready
- Why second:
  - Same overall object pattern as `assoc0`, so we can likely reuse the same style of mock object
  - Exported methods have clear user-facing value
  - Good candidate immediately after `assoc0` because the setup is similar
- Likely batch scope:
  - `print.assoc1()`
  - `summary.assoc1()`
  - `print.summary.assoc1()`
  - `plot.assoc1()`
- Expected risk: low to medium
- Main watchouts:
  - `plot.assoc1(type = "dExy")` needs matrix shapes that behave correctly under `diff()`
  - examples should avoid introducing too many fields if a smaller `type = "volcano"`-style equivalent is not available

### 3. `R/fassoc_test.R`

- Status: ready after `assoc0`/`assoc1`
- Why third:
  - Generic methods are simple wrappers, but their examples should follow whatever mock objects we settle on for `assoc0` and `assoc1`
  - Better to modernize this file after the concrete result-class files so the examples stay aligned
- Likely batch scope:
  - `print.fassoc()`
  - `summary.fassoc()`
- Expected risk: low
- Main watchouts:
  - Avoid duplicated example blocks that drift from the concrete class examples
  - Likely best to demonstrate one `assoc0` and one `assoc1` object briefly

### 4. `R/graphics.R`

- Status: investigate
- Why later:
  - `plot.ggraph()` is exported and valuable, but plotting objects here may have more hidden structure than the smaller S3 summary/print methods
  - The file already contains several examples, so the incremental value may be lower than the paired-test files
- Likely batch scope:
  - `plot.ggraph()` only
- Expected risk: medium
- Main watchouts:
  - layout reuse, optional 3D behavior, and igraph/rgl interactions
  - higher chance of example fragility than the result-object methods above

### 5. `R/ulogit.R`

- Status: mostly complete
- Why low priority:
  - Print examples already exist and are lightweight
  - Remaining work here is mainly editorial rather than structural
- Expected risk: low
- Main watchouts:
  - Low payoff compared with the remaining paired-test method files

## Deferred / lower-value targets

### `R/path_graphs.R`

- Status: mostly complete
- Reason deferred:
  - The major missing method examples have already been added

### `R/mknn_graphs.R`

- Status: mostly complete
- Reason deferred:
  - Collection summary/print examples already exist; single-graph print example was just added

### `R/iknn_graph.R`

- Status: mostly complete
- Reason deferred:
  - `summary.IkNN()` already has a direct lightweight example

## Suggested next three commits

1. `Add assoc0 method examples`
2. `Add assoc1 method examples`
3. `Add fassoc generic method examples`
