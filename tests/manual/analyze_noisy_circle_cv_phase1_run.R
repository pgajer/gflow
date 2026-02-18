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
    if (identical(parent, cur)) {
      stop("Could not locate repository root (missing DESCRIPTION).", call. = FALSE)
    }
    cur <- parent
  }
}

find.latest.run.dir <- function(repo.root) {
  base <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase1")
  if (!dir.exists(base)) stop("Run base directory does not exist: ", base, call. = FALSE)
  runs <- list.dirs(base, full.names = TRUE, recursive = FALSE)
  if (length(runs) == 0L) stop("No run directories found under: ", base, call. = FALSE)
  info <- file.info(runs)
  runs[which.max(info$mtime)]
}

timestamp.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ensure.dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

extract.n.components <- function(x) {
  if (length(x) == 0L) return(integer(0))
  m <- regexec("has ([0-9]+) connected components", x, perl = TRUE)
  parts <- regmatches(x, m)
  out <- rep(NA_integer_, length(x))
  ok <- lengths(parts) >= 2L
  if (any(ok)) {
    out[ok] <- as.integer(vapply(parts[ok], function(z) z[2L], FUN.VALUE = character(1)))
  }
  out
}

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

    log.msg(sprintf("Read attempt %d failed for %s | %s",
                    attempt, basename(src), conditionMessage(dt)))
    Sys.sleep(sleep.sec)
  }
  stop("Could not snapshot/read file after retries: ", src, call. = FALSE)
}

save.plot <- function(plot.obj, path, width = 11, height = 7, dpi = 160) {
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

  log.msg("Starting analysis snapshot.")
  log.msg("Run dir:", run.dir)
  log.msg("Snapshot tag:", snapshot.tag)

  cfg.path <- file.path(run.dir, "config", "experiment_config.rds")
  cfg <- if (file.exists(cfg.path)) readRDS(cfg.path) else NULL
  expected.total.tasks <- NA_real_
  expected.rows.per.replicate <- NA_real_
  if (!is.null(cfg)) {
    expected.total.tasks <- length(cfg$sigma_grid) * cfg$n_replicates * cfg$n_repeats * cfg$n_folds *
      length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
    expected.rows.per.replicate <- cfg$n_repeats * cfg$n_folds *
      length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
  }

  table.files <- c(
    cv_fold_scores = "cv_fold_scores.csv",
    cv_candidate_summary = "cv_candidate_summary.csv",
    selection_by_replicate = "selection_by_replicate.csv",
    oracle_metrics_by_replicate = "oracle_metrics_by_replicate.csv"
  )

  snap.tables <- list()
  for (nm in names(table.files)) {
    src <- file.path(run.dir, "tables", table.files[[nm]])
    dest <- file.path(snap.dir, table.files[[nm]])
    snap.tables[[nm]] <- safe.snapshot.read(src, dest, log.msg = log.msg)
    log.msg(sprintf("Snapshot loaded | %s | rows=%d", nm, nrow(snap.tables[[nm]])))
  }

  cv <- snap.tables$cv_fold_scores
  cand <- snap.tables$cv_candidate_summary
  sel <- snap.tables$selection_by_replicate
  oracle <- snap.tables$oracle_metrics_by_replicate

  if (nrow(cv) == 0L) stop("cv_fold_scores is empty; nothing to analyze.", call. = FALSE)

  num.cols.cv <- c("sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta", "scaled_mse", "raw_mse")
  for (cc in intersect(num.cols.cv, names(cv))) set(cv, j = cc, value = suppressWarnings(as.numeric(cv[[cc]])))
  if ("status" %in% names(cv)) set(cv, j = "status", value = as.character(cv$status))
  if ("error" %in% names(cv)) set(cv, j = "error", value = as.character(cv$error))

  num.cols.cand <- c("sigma", "replicate", "k", "n_eigenpairs", "eta", "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")
  for (cc in intersect(num.cols.cand, names(cand))) set(cand, j = cc, value = suppressWarnings(as.numeric(cand[[cc]])))

  num.cols.sel <- c("sigma", "replicate", "k", "n_eigenpairs", "eta", "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")
  for (cc in intersect(num.cols.sel, names(sel))) set(sel, j = cc, value = suppressWarnings(as.numeric(sel[[cc]])))
  if ("selected" %in% names(sel)) set(sel, j = "selected", value = as.character(sel$selected))

  num.cols.oracle <- c("sigma", "replicate", "n_rows", "rmse_xy_obs", "rmse_xy_hat", "rmse_rad_obs", "rmse_rad_hat", "var_ratio_mean")
  for (cc in intersect(num.cols.oracle, names(oracle))) set(oracle, j = cc, value = suppressWarnings(as.numeric(oracle[[cc]])))
  if ("method" %in% names(oracle)) set(oracle, j = "method", value = as.character(oracle$method))

  log.msg("1) Completion/integrity analysis")
  key.cols <- c("sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta")
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

  status.count <- if ("status" %in% names(cv)) cv[, .N, by = .(status)][order(-N)] else data.table(status = "unknown", N = observed.rows)
  fwrite(status.count, file.path(tab.dir, "status_counts.csv"))

  if (all(c("sigma", "replicate", "status") %in% names(cv))) {
    progress.by.rep <- cv[, .(
      rows = .N,
      ok_rows = sum(status == "ok", na.rm = TRUE),
      non_ok_rows = sum(status != "ok", na.rm = TRUE)
    ), by = .(sigma, replicate)][order(sigma, replicate)]
    if (is.finite(expected.rows.per.replicate) && expected.rows.per.replicate > 0) {
      progress.by.rep[, completion_fraction := rows / expected.rows.per.replicate]
    } else {
      progress.by.rep[, completion_fraction := NA_real_]
    }
    fwrite(progress.by.rep, file.path(tab.dir, "completion_by_sigma_replicate.csv"))
  }

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
      title = sprintf("Run Completion Snapshot (%s)", run.id),
      subtitle = sprintf("Snapshot %s | completion %.2f%%", snapshot.tag, 100 * completion.frac),
      x = NULL, y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save.plot(p01, file.path(fig.dir, "01_completion_status.png"), width = 12, height = 5.5)

  log.msg("2) Failure topology analysis")
  if (all(c("sigma", "k", "n_eigenpairs", "status") %in% names(cv))) {
    fail.by <- cv[, .(
      n_total = .N,
      n_fail = sum(status != "ok", na.rm = TRUE),
      n_fit_error = sum(status == "fit_error", na.rm = TRUE),
      n_refit_error = sum(status == "refit_error", na.rm = TRUE)
    ), by = .(sigma, k, n_eigenpairs)][order(sigma, n_eigenpairs, k)]
    fail.by[, fail_rate := n_fail / pmax(n_total, 1)]
    fwrite(fail.by, file.path(tab.dir, "failure_by_k_sigma_eigenpairs.csv"))

    p02 <- ggplot(fail.by, aes(x = k, y = 100 * fail_rate, color = factor(sigma), group = factor(sigma))) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.2) +
      facet_wrap(~n_eigenpairs, nrow = 1) +
      labs(
        title = "Failure Rate vs k by Sigma",
        subtitle = "Failure = status != ok",
        x = "k",
        y = "Failure rate (%)",
        color = "sigma"
      ) +
      theme_bw(base_size = 11)
    save.plot(p02, file.path(fig.dir, "02_failrate_vs_k_by_sigma.png"), width = 12, height = 4.5)
  }

  if (all(c("status", "error", "sigma", "k", "n_eigenpairs") %in% names(cv))) {
    fit.err <- cv[status == "fit_error" & !is.na(error) & nzchar(error), .N, by = .(sigma, k, n_eigenpairs, error)]
    fit.err[, n_components := extract.n.components(error)]
    comp.summary <- fit.err[!is.na(n_components), .(N = sum(N)), by = .(sigma, k, n_eigenpairs, n_components)][order(sigma, k, n_components)]
    fwrite(comp.summary, file.path(tab.dir, "failure_connected_components_summary.csv"))
  }

  log.msg("3) CV risk landscape")
  if (all(c("sigma", "k", "n_eigenpairs", "eta", "scaled_mse_mean", "raw_mse_mean") %in% names(cand))) {
    risk <- cand[is.finite(scaled_mse_mean), .(
      scaled_mse_mean_avg = mean(scaled_mse_mean, na.rm = TRUE),
      scaled_mse_mean_median = median(scaled_mse_mean, na.rm = TRUE),
      raw_mse_mean_avg = mean(raw_mse_mean, na.rm = TRUE),
      n_replicates = .N
    ), by = .(sigma, n_eigenpairs, k, eta)][order(sigma, n_eigenpairs, k, eta)]
    risk[, log10_eta := log10(eta)]
    fwrite(risk, file.path(tab.dir, "risk_landscape_mean.csv"))

    p03 <- ggplot(risk, aes(x = k, y = log10_eta, fill = scaled_mse_mean_avg)) +
      geom_tile() +
      facet_grid(n_eigenpairs ~ sigma) +
      scale_fill_viridis_c(option = "C", na.value = "grey90") +
      labs(
        title = "CV Risk Landscape",
        subtitle = "Mean scaled CV MSE across available replicates",
        x = "k",
        y = "log10(eta)",
        fill = "mean scaled MSE"
      ) +
      theme_bw(base_size = 10)
    save.plot(p03, file.path(fig.dir, "03_cv_surface_heatmaps.png"), width = 12, height = 8)
  }

  log.msg("4) Selection distributions")
  sel.ok <- if (all(c("selected", "sigma", "k", "n_eigenpairs", "eta") %in% names(sel))) {
    sel[selected %in% c("min", "one_se") & is.finite(k) & is.finite(n_eigenpairs) & is.finite(eta)]
  } else {
    data.table()
  }

  if (nrow(sel.ok) > 0L) {
    freq.k <- sel.ok[, .N, by = .(sigma, selected, k)][order(sigma, selected, k)]
    freq.eig <- sel.ok[, .N, by = .(sigma, selected, n_eigenpairs)][order(sigma, selected, n_eigenpairs)]
    freq.eta <- sel.ok[, .N, by = .(sigma, selected, eta)][order(sigma, selected, eta)]
    fwrite(freq.k, file.path(tab.dir, "selection_freq_k.csv"))
    fwrite(freq.eig, file.path(tab.dir, "selection_freq_n_eigenpairs.csv"))
    fwrite(freq.eta, file.path(tab.dir, "selection_freq_eta.csv"))

    p04 <- ggplot(freq.k, aes(x = factor(k), y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Selected k Distribution", x = "k", y = "Count", fill = "rule") +
      theme_bw(base_size = 10)
    save.plot(p04, file.path(fig.dir, "04_selected_k_distribution.png"), width = 12, height = 6)

    p05 <- ggplot(freq.eig, aes(x = factor(n_eigenpairs), y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Selected n.eigenpairs Distribution", x = "n.eigenpairs", y = "Count", fill = "rule") +
      theme_bw(base_size = 10)
    save.plot(p05, file.path(fig.dir, "05_selected_eigenpairs_distribution.png"), width = 12, height = 6)

    eta.levels <- sort(unique(freq.eta$eta))
    freq.eta[, eta_label := factor(sprintf("%.4g", eta), levels = sprintf("%.4g", eta.levels))]
    p06 <- ggplot(freq.eta, aes(x = eta_label, y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Selected eta Distribution", x = "eta", y = "Count", fill = "rule") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save.plot(p06, file.path(fig.dir, "06_selected_eta_distribution.png"), width = 12, height = 6)
  }

  log.msg("5) One-SE vs min tradeoff")
  if (nrow(sel.ok) > 0L) {
    min.sel <- sel.ok[selected == "min", .(
      sigma, replicate,
      k_min = k, n_eigenpairs_min = n_eigenpairs, eta_min = eta,
      scaled_mse_mean_min = scaled_mse_mean
    )]
    one.sel <- sel.ok[selected == "one_se", .(
      sigma, replicate,
      k_one_se = k, n_eigenpairs_one_se = n_eigenpairs, eta_one_se = eta,
      scaled_mse_mean_one_se = scaled_mse_mean
    )]
    tradeoff <- merge(min.sel, one.sel, by = c("sigma", "replicate"), all = FALSE)
    if (nrow(tradeoff) > 0L) {
      tradeoff[, `:=`(
        delta_scaled_mse = scaled_mse_mean_one_se - scaled_mse_mean_min,
        delta_k = k_one_se - k_min,
        delta_n_eigenpairs = n_eigenpairs_one_se - n_eigenpairs_min,
        delta_log10_eta = log10(pmax(eta_one_se, 1e-15) / pmax(eta_min, 1e-15))
      )]
      fwrite(tradeoff, file.path(tab.dir, "one_se_vs_min_deltas.csv"))

      p07 <- ggplot(tradeoff, aes(x = delta_n_eigenpairs, y = delta_scaled_mse, color = factor(sigma))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(size = abs(delta_k) + 1), alpha = 0.75) +
        facet_wrap(~sigma, scales = "free") +
        labs(
          title = "One-SE vs Min Tradeoff",
          subtitle = "Point size = |delta k| + 1",
          x = "delta n.eigenpairs (one_se - min)",
          y = "delta scaled CV MSE (one_se - min)",
          color = "sigma",
          size = "|delta k|+1"
        ) +
        theme_bw(base_size = 10)
      save.plot(p07, file.path(fig.dir, "07_one_se_vs_min_tradeoff.png"), width = 12, height = 6)
    }
  }

  log.msg("6) Oracle method comparison")
  method.order <- c("observed_no_smoothing", "data_smoother_default", "cv_one_se", "cv_min")
  oracle.sub <- if ("method" %in% names(oracle)) oracle[method %in% method.order] else data.table()
  if (nrow(oracle.sub) > 0L) {
    oracle.sub[, method := factor(method, levels = method.order)]
    oracle.summary <- oracle.sub[, .(
      n = .N,
      rmse_xy_hat_mean = mean(rmse_xy_hat, na.rm = TRUE),
      rmse_xy_hat_sd = sd(rmse_xy_hat, na.rm = TRUE),
      rmse_rad_hat_mean = mean(rmse_rad_hat, na.rm = TRUE),
      rmse_rad_hat_sd = sd(rmse_rad_hat, na.rm = TRUE),
      var_ratio_mean_avg = mean(var_ratio_mean, na.rm = TRUE),
      var_ratio_mean_sd = sd(var_ratio_mean, na.rm = TRUE)
    ), by = .(sigma, method)]
    fwrite(oracle.summary, file.path(tab.dir, "oracle_summary_by_sigma_method.csv"))

    p08 <- ggplot(oracle.sub, aes(x = method, y = rmse_xy_hat, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Oracle RMSE(x,y) by Method and Sigma", x = NULL, y = "rmse_xy_hat") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p08, file.path(fig.dir, "08_oracle_rmse_xy_by_method_sigma.png"), width = 12, height = 6)

    p09 <- ggplot(oracle.sub, aes(x = method, y = rmse_rad_hat, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Oracle Radial RMSE by Method and Sigma", x = NULL, y = "rmse_rad_hat") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p09, file.path(fig.dir, "09_oracle_rmse_rad_by_method_sigma.png"), width = 12, height = 6)

    p10 <- ggplot(oracle.sub, aes(x = method, y = var_ratio_mean, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_wrap(~sigma, scales = "free_y") +
      labs(title = "Variance Ratio by Method and Sigma", x = NULL, y = "var_ratio_mean") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p10, file.path(fig.dir, "10_var_ratio_by_method_sigma.png"), width = 12, height = 6)
  }

  log.msg("7) Delta vs observed baseline")
  delta.tbl <- data.table()
  if (nrow(oracle.sub) > 0L && all(c("sigma", "replicate", "method", "rmse_xy_hat", "rmse_rad_hat", "var_ratio_mean") %in% names(oracle.sub))) {
    baseline <- oracle.sub[method == "observed_no_smoothing", .(
      sigma, replicate,
      baseline_rmse_xy_hat = rmse_xy_hat,
      baseline_rmse_rad_hat = rmse_rad_hat,
      baseline_var_ratio_mean = var_ratio_mean
    )]
    delta.tbl <- merge(
      oracle.sub[method != "observed_no_smoothing"],
      baseline,
      by = c("sigma", "replicate"),
      all.x = TRUE
    )
    delta.tbl[, `:=`(
      delta_rmse_xy_hat = rmse_xy_hat - baseline_rmse_xy_hat,
      delta_rmse_rad_hat = rmse_rad_hat - baseline_rmse_rad_hat,
      delta_var_ratio_mean = var_ratio_mean - baseline_var_ratio_mean
    )]
    fwrite(delta.tbl, file.path(tab.dir, "delta_vs_observed_by_replicate.csv"))

    p11 <- ggplot(delta.tbl, aes(x = method, y = delta_rmse_xy_hat, fill = method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_wrap(~sigma, scales = "free_y") +
      labs(
        title = "Delta RMSE(x,y) vs Observed Baseline",
        subtitle = "Negative is better than observed/no-smoothing",
        x = NULL,
        y = "delta rmse_xy_hat"
      ) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p11, file.path(fig.dir, "11_delta_vs_observed.png"), width = 12, height = 6)
  }

  log.msg("8) CV vs oracle agreement")
  if (nrow(sel.ok) > 0L && nrow(oracle.sub) > 0L) {
    sel.map <- sel.ok[, .(
      sigma, replicate, selected,
      cv_scaled_mse = scaled_mse_mean,
      cv_raw_mse = raw_mse_mean
    )]
    sel.map[, method := ifelse(selected == "min", "cv_min", "cv_one_se")]
    cv.or <- merge(
      sel.map,
      oracle.sub[, .(sigma, replicate, method, rmse_xy_hat, rmse_rad_hat, var_ratio_mean)],
      by = c("sigma", "replicate", "method"),
      all.x = TRUE
    )
    fwrite(cv.or, file.path(tab.dir, "cv_vs_oracle_joined.csv"))

    corr.tbl <- cv.or[is.finite(cv_scaled_mse) & is.finite(rmse_xy_hat), .(
      n = .N,
      spearman_rho = if (.N > 1L) cor(cv_scaled_mse, rmse_xy_hat, method = "spearman") else NA_real_,
      pearson_r = if (.N > 1L) cor(cv_scaled_mse, rmse_xy_hat, method = "pearson") else NA_real_
    ), by = .(sigma, selected)]
    fwrite(corr.tbl, file.path(tab.dir, "cv_vs_oracle_correlation.csv"))

    p12 <- ggplot(cv.or, aes(x = cv_scaled_mse, y = rmse_xy_hat, color = factor(sigma), shape = selected)) +
      geom_point(alpha = 0.75) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      facet_wrap(~selected, scales = "free") +
      labs(
        title = "CV Loss vs Oracle Error",
        subtitle = "Association check for proxy-risk validity",
        x = "CV scaled MSE (selected candidate)",
        y = "oracle rmse_xy_hat",
        color = "sigma",
        shape = "rule"
      ) +
      theme_bw(base_size = 10)
    save.plot(p12, file.path(fig.dir, "12_cv_vs_oracle_scatter.png"), width = 12, height = 6)
  }

  log.msg("9) Oversmoothing frontier")
  if (nrow(delta.tbl) > 0L) {
    frontier <- delta.tbl[method %in% c("cv_min", "cv_one_se", "data_smoother_default"), .(
      sigma, replicate, method,
      var_ratio_mean,
      improvement_rmse_xy = baseline_rmse_xy_hat - rmse_xy_hat,
      improvement_rmse_rad = baseline_rmse_rad_hat - rmse_rad_hat
    )]
    fwrite(frontier, file.path(tab.dir, "oversmoothing_frontier_points.csv"))

    p13 <- ggplot(frontier, aes(x = var_ratio_mean, y = improvement_rmse_xy, color = method)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.7) +
      facet_wrap(~sigma, scales = "free") +
      labs(
        title = "Oversmoothing Frontier",
        subtitle = "Higher improvement is better; low variance ratio can indicate aggressive smoothing",
        x = "var_ratio_mean",
        y = "improvement in rmse_xy (vs observed baseline)",
        color = "method"
      ) +
      theme_bw(base_size = 10)
    save.plot(p13, file.path(fig.dir, "13_oversmoothing_frontier.png"), width = 12, height = 6)
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
