find.gflow.package.root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate package root containing DESCRIPTION.")
    }
    path <- parent
  }
}

write.tsv.gflow <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

build.retina.pseudotime.long <- function(df) {
  methods <- c("geodesic", "diffusion", "potential")
  out <- lapply(methods, function(method) {
    data.frame(
      hvg_label = as.character(df$hvg.label),
      k = as.integer(df$k),
      status = as.character(df$status),
      n = as.integer(df$n),
      n_components = as.integer(df$n.components),
      lcc_frac = as.numeric(df$lcc.frac),
      method = method,
      pearson_r = as.numeric(df[[paste0(method, "_pearson_r")]]),
      spearman_rho = as.numeric(df[[paste0(method, "_spearman_rho")]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

build.retina.asset.manifest <- function(asset.rows, manifest.path) {
  manifest <- do.call(rbind, lapply(asset.rows, function(x) {
    data.frame(
      asset_id = x$asset_id,
      asset_type = x$asset_type,
      relative_path = x$relative_path,
      source_path = x$source_path,
      source_kind = x$source_kind,
      size_bytes = as.numeric(file.info(x$output_path)$size),
      md5 = unname(tools::md5sum(x$output_path)),
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  }))
  write.tsv.gflow(manifest, manifest.path)
  manifest
}

build_retina_vignette_assets <- function(
  package_root = NULL,
  cell_cycle_dir = Sys.getenv(
    "GFLOW_RETINA_CELL_CYCLE_DIR",
    unset = "/Users/pgajer/current_projects/cell_cycle"
  ),
  retina_example_dir = Sys.getenv(
    "GFLOW_RETINA_EXAMPLE_DIR",
    unset = "/Users/pgajer/current_projects/gflow_examples/retina_cell_cycle"
  ),
  retina_run_id = Sys.getenv(
    "GFLOW_RETINA_RUN_ID",
    unset = "retina_intersection_api_rerun_20260323a__smoke"
  ),
  overwrite = TRUE,
  verbose = TRUE
) {
  if (is.null(package_root) || !nzchar(package_root)) {
    package_root <- find.gflow.package.root()
  } else {
    package_root <- normalizePath(package_root, winslash = "/", mustWork = TRUE)
  }

  cell_cycle_dir <- normalizePath(path.expand(cell_cycle_dir), winslash = "/", mustWork = TRUE)
  retina_example_dir <- normalizePath(path.expand(retina_example_dir), winslash = "/", mustWork = TRUE)
  retina_run_root <- normalizePath(
    file.path(retina_example_dir, "results", retina_run_id),
    winslash = "/",
    mustWork = TRUE
  )

  out.root <- file.path(package_root, "inst", "extdata", "retina_vignette")
  tables.dir <- file.path(out.root, "tables")
  figures.dir <- file.path(out.root, "figures")
  dir.create(tables.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(verbose)) {
    message("Writing retina vignette assets into: ", out.root)
  }

  asset.rows <- list()
  add.asset <- function(asset_id, asset_type, relative_path, source_path, source_kind) {
    asset.rows[[length(asset.rows) + 1L]] <<- list(
      asset_id = asset_id,
      asset_type = asset_type,
      relative_path = relative_path,
      source_path = source_path,
      source_kind = source_kind,
      output_path = file.path(out.root, relative_path)
    )
  }

  copy.or.write <- function(df = NULL, source_path = NULL, out.path, writer = c("tsv", "copy")) {
    writer <- match.arg(writer)
    if (!isTRUE(overwrite) && file.exists(out.path)) {
      return(invisible(NULL))
    }
    if (writer == "copy") {
      file.copy(source_path, out.path, overwrite = TRUE)
    } else {
      write.tsv.gflow(df, out.path)
    }
  }

  src.preproc <- file.path(retina_run_root, "02_prepare_features", "preprocessing_views.tsv")
  preproc <- utils::read.delim(src.preproc, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  out.preproc <- file.path(tables.dir, "preprocessing_views.tsv")
  copy.or.write(df = preproc, out.path = out.preproc, writer = "tsv")
  add.asset("T01", "table", "tables/preprocessing_views.tsv", src.preproc, "retina_example")

  src.pseudotime.full <- file.path(cell_cycle_dir, "data", "pseudotime_method_sweep_summary.csv")
  pseudo.full <- utils::read.csv(src.pseudotime.full, stringsAsFactors = FALSE, check.names = FALSE)
  pseudo.full.long <- build.retina.pseudotime.long(pseudo.full)
  out.pseudotime.full <- file.path(tables.dir, "pseudotime_sweep_full.tsv")
  copy.or.write(df = pseudo.full.long, out.path = out.pseudotime.full, writer = "tsv")
  add.asset("T02", "table", "tables/pseudotime_sweep_full.tsv", src.pseudotime.full, "cell_cycle")

  src.pseudotime.ext <- file.path(cell_cycle_dir, "data", "pseudotime_method_sweep_3k_k09_35_summary.csv")
  pseudo.ext <- utils::read.csv(src.pseudotime.ext, stringsAsFactors = FALSE, check.names = FALSE)
  pseudo.ext.long <- build.retina.pseudotime.long(pseudo.ext)
  out.pseudotime.ext <- file.path(tables.dir, "pseudotime_sweep_3k_extension.tsv")
  copy.or.write(df = pseudo.ext.long, out.path = out.pseudotime.ext, writer = "tsv")
  add.asset("T03", "table", "tables/pseudotime_sweep_3k_extension.tsv", src.pseudotime.ext, "cell_cycle")

  src.targeted.curve <- file.path(
    cell_cycle_dir,
    "data",
    "notch_cellcycle_condexp_gcv_sweep_3k_k03_30",
    "plots",
    "gcv_norm_summary_by_k.csv"
  )
  targeted.curve <- utils::read.csv(src.targeted.curve, stringsAsFactors = FALSE, check.names = FALSE)
  out.targeted.curve <- file.path(tables.dir, "targeted_k_curve.tsv")
  copy.or.write(df = targeted.curve, out.path = out.targeted.curve, writer = "tsv")
  add.asset("T04", "table", "tables/targeted_k_curve.tsv", src.targeted.curve, "cell_cycle")

  src.targeted.selection <- file.path(
    cell_cycle_dir,
    "data",
    "notch_cellcycle_condexp_gcv_sweep_3k_k03_30",
    "tables",
    "k_selection_from_distribution.csv"
  )
  targeted.selection <- utils::read.csv(src.targeted.selection, stringsAsFactors = FALSE, check.names = FALSE)
  names(targeted.selection)[names(targeted.selection) == "k_selected"] <- "selected_k"
  out.targeted.selection <- file.path(tables.dir, "targeted_k_selection.tsv")
  copy.or.write(df = targeted.selection, out.path = out.targeted.selection, writer = "tsv")
  add.asset("T05", "table", "tables/targeted_k_selection.tsv", src.targeted.selection, "cell_cycle")

  selected.files <- c(
    file.path(retina_run_root, "03_fit_base_graph", "lognorm_hvg__pca", "selected_k.tsv"),
    file.path(retina_run_root, "03_fit_base_graph", "residual_hvg__pca", "selected_k.tsv")
  )
  branch.selected <- do.call(rbind, lapply(selected.files, function(path) {
    utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }))
  out.branch.selected <- file.path(tables.dir, "branch_selected_k.tsv")
  copy.or.write(df = branch.selected, out.path = out.branch.selected, writer = "tsv")
  add.asset("T06", "table", "tables/branch_selected_k.tsv", paste(selected.files, collapse = "; "), "retina_example")

  validation.files <- c(
    file.path(retina_run_root, "03_fit_base_graph", "lognorm_hvg__pca", "validation_diagnostics.tsv"),
    file.path(retina_run_root, "03_fit_base_graph", "residual_hvg__pca", "validation_diagnostics.tsv")
  )
  branch.validation <- do.call(rbind, lapply(validation.files, function(path) {
    utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }))
  out.branch.validation <- file.path(tables.dir, "branch_validation_diagnostics.tsv")
  copy.or.write(df = branch.validation, out.path = out.branch.validation, writer = "tsv")
  add.asset("T07", "table", "tables/branch_validation_diagnostics.tsv", paste(validation.files, collapse = "; "), "retina_example")

  src.outcomes <- file.path(retina_run_root, "04_refit_outcomes", "outcome_surface_summary.tsv")
  outcomes <- utils::read.delim(src.outcomes, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  out.outcomes <- file.path(tables.dir, "outcome_surface_summary.tsv")
  copy.or.write(df = outcomes, out.path = out.outcomes, writer = "tsv")
  add.asset("T08", "table", "tables/outcome_surface_summary.tsv", src.outcomes, "retina_example")

  src.modules <- file.path(
    cell_cycle_dir,
    "data",
    "notch_vs_cellcycle_condexp_winner_3k_compare",
    "biologist_question_scan",
    "tables",
    "gene_clusters_by_k10_profile.csv"
  )
  modules <- utils::read.csv(src.modules, stringsAsFactors = FALSE, check.names = FALSE)
  out.modules <- file.path(tables.dir, "gene_module_membership.tsv")
  copy.or.write(df = modules, out.path = out.modules, writer = "tsv")
  add.asset("T09", "table", "tables/gene_module_membership.tsv", src.modules, "cell_cycle")

  src.bifur <- file.path(
    cell_cycle_dir,
    "data",
    "notch_vs_cellcycle_condexp_winner_3k_compare",
    "biologist_followup_scan",
    "tables",
    "bifurcation_pair_summary_for_biology_followup.csv"
  )
  bifur <- utils::read.csv(src.bifur, stringsAsFactors = FALSE, check.names = FALSE)
  out.bifur <- file.path(tables.dir, "module_bifurcation_pairs.tsv")
  copy.or.write(df = bifur, out.path = out.bifur, writer = "tsv")
  add.asset("T10", "table", "tables/module_bifurcation_pairs.tsv", src.bifur, "cell_cycle")

  src.phase.scan <- file.path(
    cell_cycle_dir,
    "data",
    "notch_vs_cellcycle_condexp_winner_3k_compare",
    "biologist_followup_scan",
    "tables",
    "phase_stratified_bifurcation_scan.csv"
  )
  phase.scan <- utils::read.csv(src.phase.scan, stringsAsFactors = FALSE, check.names = FALSE)
  out.phase.scan <- file.path(tables.dir, "phase_bifurcation_scan.tsv")
  copy.or.write(df = phase.scan, out.path = out.phase.scan, writer = "tsv")
  add.asset("T11", "table", "tables/phase_bifurcation_scan.tsv", src.phase.scan, "cell_cycle")

  src.handoff <- file.path(
    cell_cycle_dir,
    "data",
    "notch_vs_cellcycle_condexp_winner_3k_compare",
    "biologist_followup_scan",
    "tables",
    "phase_age_handoff_validation_top_pairs.csv"
  )
  handoff <- utils::read.csv(src.handoff, stringsAsFactors = FALSE, check.names = FALSE)
  out.handoff <- file.path(tables.dir, "phase_age_handoff_validation.tsv")
  copy.or.write(df = handoff, out.path = out.handoff, writer = "tsv")
  add.asset("T12", "table", "tables/phase_age_handoff_validation.tsv", src.handoff, "cell_cycle")

  src.module.figure <- file.path(
    cell_cycle_dir,
    "data",
    "notch_vs_cellcycle_condexp_winner_3k_compare",
    "biologist_question_scan",
    "figures",
    "cluster_mean_trajectories_k10_vs_k20.png"
  )
  out.module.figure <- file.path(figures.dir, "cluster_mean_trajectories_k10_vs_k20.png")
  copy.or.write(source_path = src.module.figure, out.path = out.module.figure, writer = "copy")
  add.asset("F01", "figure", "figures/cluster_mean_trajectories_k10_vs_k20.png", src.module.figure, "cell_cycle")

  manifest.path <- file.path(out.root, "manifest.tsv")
  build.retina.asset.manifest(asset.rows = asset.rows, manifest.path = manifest.path)
  if (isTRUE(verbose)) {
    message("Wrote manifest: ", manifest.path)
  }

  invisible(
    list(
      out_root = out.root,
      tables_dir = tables.dir,
      figures_dir = figures.dir,
      manifest = manifest.path
    )
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  build_retina_vignette_assets()
}
