#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

parse.cli.args <- function(args) {
  out <- list()
  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--[^=]+=", a)) next
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=(.*)$", "\\1", a)
    out[[key]] <- val
  }
  out
}

find.repo.root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Could not locate repository root (missing DESCRIPTION).", call. = FALSE)
    cur <- parent
  }
}

find.latest.run.dir <- function(repo.root) {
  base <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase3_highdim")
  if (!dir.exists(base)) stop("Run base directory does not exist: ", base, call. = FALSE)
  runs <- list.dirs(base, full.names = TRUE, recursive = FALSE)
  if (length(runs) == 0L) stop("No run directories found under: ", base, call. = FALSE)
  info <- file.info(runs)
  runs[which.max(info$mtime)]
}

timestamp.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ensure.dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

safe.snapshot.read <- function(src, dest, log.msg, max.tries = 6L, sleep.sec = 1) {
  if (!file.exists(src)) {
    log.msg("Source file missing (skipped):", src)
    return(data.table())
  }
  for (attempt in seq_len(max.tries)) {
    ok.copy <- file.copy(src, dest, overwrite = TRUE)
    if (!ok.copy) {
      log.msg(sprintf("Copy attempt %d failed for %s", attempt, basename(src)))
      Sys.sleep(sleep.sec)
      next
    }

    dt <- tryCatch(
      fread(dest, showProgress = FALSE, na.strings = c("NA", "")),
      error = function(e) e
    )
    if (!inherits(dt, "error")) return(dt)

    log.msg(sprintf("Read attempt %d failed for %s | %s", attempt, basename(src), conditionMessage(dt)))
    Sys.sleep(sleep.sec)
  }
  stop("Could not snapshot/read file after retries: ", src, call. = FALSE)
}

save.plot <- function(plot.obj, path, width = 12, height = 6.5, dpi = 160) {
  ggsave(filename = path, plot = plot.obj, width = width, height = height, dpi = dpi)
}

format.run.id <- function(run.dir) basename(normalizePath(run.dir, winslash = "/", mustWork = FALSE))

main <- function() {
  args <- parse.cli.args(commandArgs(trailingOnly = TRUE))
  repo.root <- find.repo.root()
  setwd(repo.root)

  run.dir <- args$run_dir %||% find.latest.run.dir(repo.root)
  run.dir <- normalizePath(run.dir, winslash = "/", mustWork = TRUE)
  run.id <- format.run.id(run.dir)
  snapshot.tag <- args$snapshot_tag %||% format(Sys.time(), "%Y%m%d_%H%M%S")

  analysis.dir <- file.path(run.dir, "analysis")
  fig.dir <- file.path(analysis.dir, "figures")
  tab.dir <- file.path(analysis.dir, "tables")
  data.dir <- file.path(analysis.dir, "data")
  log.dir <- file.path(analysis.dir, "logs")
  snap.dir <- file.path(data.dir, "snapshots", snapshot.tag)

  ensure.dir(analysis.dir)
  ensure.dir(fig.dir)
  ensure.dir(tab.dir)
  ensure.dir(data.dir)
  ensure.dir(log.dir)
  ensure.dir(snap.dir)

  log.path <- file.path(log.dir, sprintf("analysis_%s.log", snapshot.tag))
  log.con <- file(log.path, open = "a")
  on.exit(close(log.con), add = TRUE)

  log.msg <- function(...) {
    txt <- paste(...)
    line <- sprintf("[%s] %s", timestamp.now(), txt)
    cat(line, "\n")
    writeLines(line, con = log.con, sep = "\n")
    flush(log.con)
    flush.console()
  }

  log.msg("Starting phase-3 analysis snapshot.")
  log.msg("Run dir:", run.dir)
  log.msg("Snapshot tag:", snapshot.tag)

  cfg.path <- file.path(run.dir, "config", "experiment_config.rds")
  cfg <- if (file.exists(cfg.path)) readRDS(cfg.path) else NULL

  expected.total.tasks <- NA_real_
  expected.rows.per.replicate <- NA_real_
  if (!is.null(cfg)) {
    sg <- if (!is.null(cfg$scenario_grid)) as.data.frame(cfg$scenario_grid) else data.frame()
    if (nrow(sg) > 0L) {
      expected.total.tasks <- sum(as.integer(sg$n_replicates)) * cfg$n_repeats * cfg$n_folds *
        length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
    }
    expected.rows.per.replicate <- cfg$n_repeats * cfg$n_folds *
      length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
  }

  table.files <- c(
    graph_cv_fold_scores = "graph_cv_fold_scores.csv",
    graph_cv_candidate_summary = "graph_cv_candidate_summary.csv",
    graph_selection_by_replicate = "graph_selection_by_replicate.csv",
    classical_selection_by_replicate = "classical_selection_by_replicate.csv",
    classical_cv_candidate_summary = "classical_cv_candidate_summary.csv",
    oracle_metrics_by_replicate = "oracle_metrics_by_replicate.csv",
    runtime_by_replicate = "runtime_by_replicate.csv"
  )

  snap.tables <- list()
  for (nm in names(table.files)) {
    src <- file.path(run.dir, "tables", table.files[[nm]])
    dest <- file.path(snap.dir, table.files[[nm]])
    snap.tables[[nm]] <- safe.snapshot.read(src, dest, log.msg = log.msg)
    log.msg(sprintf("Snapshot loaded | %s | rows=%d", nm, nrow(snap.tables[[nm]])))
  }

  cv <- snap.tables$graph_cv_fold_scores
  cand.g <- snap.tables$graph_cv_candidate_summary
  sel.g <- snap.tables$graph_selection_by_replicate
  sel.c <- snap.tables$classical_selection_by_replicate
  cand.c <- snap.tables$classical_cv_candidate_summary
  oracle <- snap.tables$oracle_metrics_by_replicate
  runtime <- snap.tables$runtime_by_replicate

  if (nrow(cv) == 0L) stop("graph_cv_fold_scores is empty; nothing to analyze.", call. = FALSE)

  num.cols.cv <- c("p", "sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta",
                   "scaled_mse", "raw_mse", "n_test_total", "n_test_modeled", "n_train_total",
                   "n_train_modeled", "modeled_fraction")
  for (cc in intersect(num.cols.cv, names(cv))) set(cv, j = cc, value = suppressWarnings(as.numeric(cv[[cc]])))
  for (cc in intersect(c("embedding_type", "status", "error"), names(cv))) set(cv, j = cc, value = as.character(cv[[cc]]))

  num.cols.sel.g <- c("p", "sigma", "replicate", "k", "n_eigenpairs", "eta", "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")
  for (cc in intersect(num.cols.sel.g, names(sel.g))) set(sel.g, j = cc, value = suppressWarnings(as.numeric(sel.g[[cc]])))
  for (cc in intersect(c("embedding_type", "selected"), names(sel.g))) set(sel.g, j = cc, value = as.character(sel.g[[cc]]))

  num.cols.sel.c <- c("p", "sigma", "replicate", "param_value", "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")
  for (cc in intersect(num.cols.sel.c, names(sel.c))) set(sel.c, j = cc, value = suppressWarnings(as.numeric(sel.c[[cc]])))
  for (cc in intersect(c("embedding_type", "method", "selected", "param_name"), names(sel.c))) set(sel.c, j = cc, value = as.character(sel.c[[cc]]))

  num.cols.oracle <- c("p", "sigma", "replicate", "n_rows", "n_cols", "modeled_fraction",
                       "rmse_all_obs", "rmse_all_hat", "rmse_latent_obs", "rmse_latent_hat", "var_ratio_mean")
  for (cc in intersect(num.cols.oracle, names(oracle))) set(oracle, j = cc, value = suppressWarnings(as.numeric(oracle[[cc]])))
  for (cc in intersect(c("embedding_type", "method", "notes"), names(oracle))) set(oracle, j = cc, value = as.character(oracle[[cc]]))

  num.cols.rt <- c("p", "sigma", "replicate", "elapsed_sec")
  for (cc in intersect(num.cols.rt, names(runtime))) set(runtime, j = cc, value = suppressWarnings(as.numeric(runtime[[cc]])))
  for (cc in intersect(c("embedding_type", "method_group"), names(runtime))) set(runtime, j = cc, value = as.character(runtime[[cc]]))

  log.msg("1) Completion/integrity analysis")
  key.cols <- c("embedding_type", "p", "sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta")
  unique.keys <- if (all(key.cols %in% names(cv))) uniqueN(cv, by = key.cols) else NA_integer_
  observed.rows <- nrow(cv)
  completion.frac <- if (is.finite(expected.total.tasks) && expected.total.tasks > 0) observed.rows / expected.total.tasks else NA_real_

  completion.summary <- data.table(
    run_id = run.id,
    snapshot_tag = snapshot.tag,
    expected_total_rows = expected.total.tasks,
    observed_rows = observed.rows,
    unique_task_keys = unique.keys,
    duplicate_rows = observed.rows - unique.keys,
    completion_fraction = completion.frac
  )
  fwrite(completion.summary, file.path(tab.dir, "completion_summary.csv"))

  status.count <- cv[, .N, by = .(status)][order(-N)]
  fwrite(status.count, file.path(tab.dir, "status_counts.csv"))

  completion.by.group <- cv[, .(
    rows = .N,
    ok_rows = sum(status == "ok", na.rm = TRUE),
    non_ok_rows = sum(status != "ok", na.rm = TRUE)
  ), by = .(embedding_type, p, sigma, replicate)][order(embedding_type, p, sigma, replicate)]

  if (is.finite(expected.rows.per.replicate) && expected.rows.per.replicate > 0) {
    completion.by.group[, completion_fraction := rows / expected.rows.per.replicate]
  } else {
    completion.by.group[, completion_fraction := NA_real_]
  }
  fwrite(completion.by.group, file.path(tab.dir, "completion_by_group_replicate.csv"))

  counts.panel <- data.table(
    panel = "Task counts",
    label = c("observed_rows", "expected_rows"),
    value = c(observed.rows, expected.total.tasks)
  )[is.finite(value)]

  status.panel <- copy(status.count)
  status.panel[, `:=`(panel = "Status counts", label = status, value = N)]
  status.panel <- status.panel[, .(panel, label, value)]

  completion.plot.df <- rbindlist(list(counts.panel, status.panel), fill = TRUE)
  p01 <- ggplot(completion.plot.df, aes(x = label, y = value, fill = label)) +
    geom_col(width = 0.75) +
    facet_wrap(~panel, scales = "free_y", nrow = 1) +
    scale_y_continuous(labels = comma) +
    labs(
      title = sprintf("Phase-3 Run Completion Snapshot (%s)", run.id),
      subtitle = sprintf("Snapshot %s | completion %.2f%%", snapshot.tag, 100 * completion.frac),
      x = NULL, y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save.plot(p01, file.path(fig.dir, "01_completion_status.png"), width = 12, height = 5.5)

  log.msg("2) Graph failure topology")
  fail.by <- cv[, .(
    n_total = .N,
    n_fail = sum(status != "ok", na.rm = TRUE),
    n_graph_disconnected = sum(status == "graph_disconnected", na.rm = TRUE),
    n_fit_error = sum(status == "fit_error", na.rm = TRUE),
    n_refit_error = sum(status == "refit_error", na.rm = TRUE)
  ), by = .(embedding_type, p, sigma, k, n_eigenpairs)][order(embedding_type, p, sigma, n_eigenpairs, k)]
  fail.by[, fail_rate := n_fail / pmax(n_total, 1)]
  fwrite(fail.by, file.path(tab.dir, "graph_failure_by_group_k_eigenpairs.csv"))

  p02 <- ggplot(fail.by, aes(x = k, y = 100 * fail_rate, color = factor(sigma), group = factor(sigma))) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.1) +
    facet_grid(embedding_type + p ~ n_eigenpairs) +
    labs(
      title = "Graph Failure Rate vs k",
      subtitle = "Failure = status != ok",
      x = "k",
      y = "Failure rate (%)",
      color = "sigma"
    ) +
    theme_bw(base_size = 8)
  save.plot(p02, file.path(fig.dir, "02_graph_failrate_vs_k.png"), width = 14, height = 9)

  log.msg("3) Graph selection distribution")
  sel.g.ok <- sel.g[selected %in% c("min", "one_se") & is.finite(k)]
  if (nrow(sel.g.ok) > 0L) {
    freq.k <- sel.g.ok[, .N, by = .(embedding_type, p, sigma, selected, k)][order(embedding_type, p, sigma, selected, k)]
    fwrite(freq.k, file.path(tab.dir, "graph_selection_freq_k_by_group.csv"))

    p03 <- ggplot(freq.k, aes(x = factor(k), y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_grid(embedding_type + p ~ sigma, scales = "free_y") +
      labs(title = "Graph Selected k Distribution", x = "k", y = "Count", fill = "rule") +
      theme_bw(base_size = 8)
    save.plot(p03, file.path(fig.dir, "03_graph_selected_k_distribution.png"), width = 14, height = 8)
  }

  log.msg("4) Oracle method comparison")
  method.order <- c(
    "observed_no_smoothing",
    "classical_knn_mean_cv",
    "classical_gaussian_kernel_cv",
    "cv_one_se",
    "cv_min"
  )
  oracle.sub <- oracle[method %in% method.order]
  if (nrow(oracle.sub) > 0L) {
    oracle.summary <- oracle.sub[, .(
      n = .N,
      rmse_all_hat_mean = mean(rmse_all_hat, na.rm = TRUE),
      rmse_all_hat_sd = sd(rmse_all_hat, na.rm = TRUE),
      rmse_latent_hat_mean = mean(rmse_latent_hat, na.rm = TRUE),
      rmse_latent_hat_sd = sd(rmse_latent_hat, na.rm = TRUE),
      var_ratio_mean_avg = mean(var_ratio_mean, na.rm = TRUE)
    ), by = .(embedding_type, p, sigma, method)]
    fwrite(oracle.summary, file.path(tab.dir, "oracle_summary_by_group_method.csv"))

    plot.all <- copy(oracle.sub)
    plot.all[, method := factor(method, levels = method.order)]

    p04 <- ggplot(plot.all, aes(x = method, y = rmse_all_hat, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_grid(embedding_type + p ~ sigma, scales = "free_y") +
      labs(title = "Oracle RMSE (All Dimensions) by Method", x = NULL, y = "rmse_all_hat") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p04, file.path(fig.dir, "04_oracle_rmse_all_by_method.png"), width = 14, height = 9)

    p05 <- ggplot(plot.all, aes(x = method, y = rmse_latent_hat, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_grid(embedding_type + p ~ sigma, scales = "free_y") +
      labs(title = "Oracle RMSE (Latent 2D Signal) by Method", x = NULL, y = "rmse_latent_hat") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p05, file.path(fig.dir, "05_oracle_rmse_latent_by_method.png"), width = 14, height = 9)
  }

  log.msg("5) Delta and win-rate vs no smoothing")
  delta.tbl <- data.table()
  if (nrow(oracle.sub) > 0L) {
    baseline <- oracle.sub[method == "observed_no_smoothing", .(
      embedding_type, p, sigma, replicate,
      baseline_rmse_all_hat = rmse_all_hat,
      baseline_rmse_latent_hat = rmse_latent_hat
    )]

    delta.tbl <- merge(
      oracle.sub[method != "observed_no_smoothing"],
      baseline,
      by = c("embedding_type", "p", "sigma", "replicate"),
      all.x = TRUE
    )

    delta.tbl[, `:=`(
      delta_rmse_all_hat = rmse_all_hat - baseline_rmse_all_hat,
      delta_rmse_latent_hat = rmse_latent_hat - baseline_rmse_latent_hat,
      improves_all_vs_observed = rmse_all_hat < baseline_rmse_all_hat,
      improves_latent_vs_observed = rmse_latent_hat < baseline_rmse_latent_hat
    )]
    fwrite(delta.tbl, file.path(tab.dir, "delta_vs_observed_by_replicate.csv"))

    win.tbl <- delta.tbl[, .(
      n = .N,
      win_rate_all = mean(improves_all_vs_observed, na.rm = TRUE),
      win_rate_latent = mean(improves_latent_vs_observed, na.rm = TRUE),
      mean_delta_latent = mean(delta_rmse_latent_hat, na.rm = TRUE)
    ), by = .(embedding_type, p, sigma, method)][order(embedding_type, p, sigma, method)]
    fwrite(win.tbl, file.path(tab.dir, "win_rate_vs_observed_by_group_method.csv"))

    p06 <- ggplot(win.tbl, aes(x = sigma, y = 100 * win_rate_latent, color = method)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.8) +
      facet_grid(embedding_type ~ p) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(
        title = "How Often Method Beats No Smoothing (Latent RMSE)",
        subtitle = "Share of replicates with lower latent oracle RMSE than observed baseline",
        x = "sigma",
        y = "Win rate (%)",
        color = "method"
      ) +
      theme_bw(base_size = 9)
    save.plot(p06, file.path(fig.dir, "06_win_rate_vs_observed_latent.png"), width = 14, height = 8)

    rank.tbl <- win.tbl[order(mean_delta_latent), .(
      method = method[1L],
      best_mean_delta_latent = mean_delta_latent[1L],
      best_win_rate_latent = win_rate_latent[1L]
    ), by = .(embedding_type, p, sigma)]
    fwrite(rank.tbl, file.path(tab.dir, "best_method_by_group_sigma.csv"))

    heat.tbl <- win.tbl[, .(
      mean_delta_latent = mean(mean_delta_latent, na.rm = TRUE)
    ), by = .(embedding_type, p, sigma, method)]
    p07 <- ggplot(heat.tbl, aes(x = factor(sigma), y = method, fill = mean_delta_latent)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.3f", mean_delta_latent)), size = 2.5) +
      facet_grid(embedding_type ~ p) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(
        title = "Mean Delta (Latent RMSE) vs No Smoothing",
        subtitle = "Negative values indicate improvement",
        x = "sigma",
        y = "method",
        fill = "mean delta"
      ) +
      theme_bw(base_size = 8)
    save.plot(p07, file.path(fig.dir, "07_mean_delta_latent_heatmap.png"), width = 14, height = 8)
  }

  log.msg("6) Runtime summary")
  if (nrow(runtime) > 0L) {
    rt.summary <- runtime[, .(
      n = .N,
      elapsed_sec_mean = mean(elapsed_sec, na.rm = TRUE),
      elapsed_sec_median = median(elapsed_sec, na.rm = TRUE),
      elapsed_sec_sd = sd(elapsed_sec, na.rm = TRUE)
    ), by = .(embedding_type, p, sigma, method_group)][order(embedding_type, p, sigma, method_group)]
    fwrite(rt.summary, file.path(tab.dir, "runtime_summary_by_group_method_group.csv"))

    p08 <- ggplot(runtime, aes(x = method_group, y = elapsed_sec, fill = method_group)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_grid(embedding_type + p ~ sigma, scales = "free_y") +
      labs(
        title = "Runtime by Method Group",
        x = NULL,
        y = "elapsed seconds"
      ) +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
    save.plot(p08, file.path(fig.dir, "08_runtime_by_method_group.png"), width = 14, height = 8)
  }

  log.msg("7) CV vs oracle agreement (graph methods)")
  if (nrow(sel.g.ok) > 0L && nrow(oracle.sub) > 0L) {
    sel.map <- sel.g.ok[, .(
      embedding_type, p, sigma, replicate, selected,
      cv_scaled_mse = scaled_mse_mean,
      cv_raw_mse = raw_mse_mean
    )]
    sel.map[, method := ifelse(selected == "min", "cv_min", "cv_one_se")]

    cv.or <- merge(
      sel.map,
      oracle.sub[, .(embedding_type, p, sigma, replicate, method, rmse_latent_hat, rmse_all_hat)],
      by = c("embedding_type", "p", "sigma", "replicate", "method"),
      all.x = TRUE
    )
    fwrite(cv.or, file.path(tab.dir, "cv_vs_oracle_joined_graph.csv"))

    corr.tbl <- cv.or[is.finite(cv_scaled_mse) & is.finite(rmse_latent_hat), .(
      n = .N,
      spearman_rho = if (.N > 1L) cor(cv_scaled_mse, rmse_latent_hat, method = "spearman") else NA_real_,
      pearson_r = if (.N > 1L) cor(cv_scaled_mse, rmse_latent_hat, method = "pearson") else NA_real_
    ), by = .(embedding_type, p, sigma, selected)]
    fwrite(corr.tbl, file.path(tab.dir, "cv_vs_oracle_correlation_graph.csv"))

    p09 <- ggplot(cv.or, aes(x = cv_scaled_mse, y = rmse_latent_hat, color = factor(sigma), shape = selected)) +
      geom_point(alpha = 0.75) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.55) +
      facet_grid(selected + embedding_type ~ p, scales = "free") +
      labs(
        title = "Graph CV Loss vs Oracle Latent Error",
        subtitle = "Proxy-risk validity check in high dimensions",
        x = "CV scaled MSE",
        y = "oracle rmse_latent_hat",
        color = "sigma",
        shape = "rule"
      ) +
      theme_bw(base_size = 8)
    save.plot(p09, file.path(fig.dir, "09_graph_cv_vs_oracle_latent.png"), width = 14, height = 9)
  }

  manifest.paths <- list.files(analysis.dir, recursive = TRUE, full.names = TRUE)
  manifest <- data.table(path = manifest.paths)
  manifest[, rel_path := sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", analysis.dir), "/?"), "", path)]
  manifest[, size_bytes := file.info(path)$size]
  manifest <- manifest[order(rel_path)]
  fwrite(manifest, file.path(analysis.dir, "analysis_manifest.csv"))

  summary.txt <- c(
    sprintf("run_id: %s", run.id),
    sprintf("run_dir: %s", run.dir),
    sprintf("snapshot_tag: %s", snapshot.tag),
    sprintf("expected_total_rows: %s", ifelse(is.finite(expected.total.tasks), format(expected.total.tasks, scientific = FALSE), "NA")),
    sprintf("observed_rows: %s", format(observed.rows, scientific = FALSE)),
    sprintf("completion_fraction: %.6f", completion.frac),
    sprintf("analysis_dir: %s", analysis.dir)
  )
  writeLines(summary.txt, file.path(analysis.dir, "analysis_summary.txt"))

  log.msg("Analysis complete.")
  log.msg("Analysis dir:", analysis.dir)
}

main()
