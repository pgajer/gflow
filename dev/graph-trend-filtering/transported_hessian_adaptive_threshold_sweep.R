source(file.path("/Users/pgajer/current_projects/gflow",
                 "dev/graph-trend-filtering/transported_hessian_backend_validation.R"))

transported.hessian.adaptive.threshold.policies <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  policies <- list(
    list(
      name = "baseline_none",
      label = "No filter",
      rule.family = "baseline",
      args = list(),
      description = "All candidate soft-transport rows are retained."
    ),
    list(
      name = "baseline_fixed_geometric",
      label = "Fixed angle/length",
      rule.family = "baseline",
      args = list(max.match.angle = pi / 3,
                  max.length.relative.error = 0.50),
      description = "Fixed geometric gate: best angle at most pi/3 and relative length error at most 0.50."
    )
  )

  if (identical(mode, "smoke")) {
    quantile.grid <- expand.grid(
      match.score.quantile = c(0.25, 0.50),
      match.margin.quantile = c(0.50, 0.75),
      max.effective.match.fraction = c(NA_real_, 0.50)
    )
    robust.grid <- expand.grid(
      min.best.score.z = c(0, 0.5, 1.0),
      min.margin.z = c(0, 0.5),
      max.effective.match.fraction = c(NA_real_, 0.50)
    )
  } else {
    quantile.grid <- expand.grid(
      match.score.quantile = c(0.25, 0.50),
      match.margin.quantile = c(0.25, 0.50, 0.75, 1.00),
      max.effective.match.fraction = c(NA_real_, 0.50)
    )
    robust.grid <- expand.grid(
      min.best.score.z = c(0, 0.5, 1.0, 1.5),
      min.margin.z = c(0, 0.5, 1.0),
      max.effective.match.fraction = c(NA_real_, 0.50)
    )
  }

  for (idx in seq_len(nrow(quantile.grid))) {
    row <- quantile.grid[idx, ]
    suffix <- if (is.na(row$max.effective.match.fraction)) {
      "noef"
    } else {
      sprintf("ef%.2f", row$max.effective.match.fraction)
    }
    args <- list(
      match.threshold.rule = "local.quantile",
      match.score.quantile = row$match.score.quantile,
      match.margin.quantile = row$match.margin.quantile
    )
    if (!is.na(row$max.effective.match.fraction)) {
      args$max.effective.match.fraction <- row$max.effective.match.fraction
    }
    policies[[length(policies) + 1L]] <- list(
      name = sprintf("quantile_qs%.2f_qm%.2f_%s",
                     row$match.score.quantile,
                     row$match.margin.quantile,
                     suffix),
      label = sprintf("Quantile qs=%.2f qm=%.2f %s",
                      row$match.score.quantile,
                      row$match.margin.quantile,
                      suffix),
      rule.family = "local_quantile",
      args = args,
      description = "Adaptive local quantile gate on best score and best-vs-second margin."
    )
  }

  for (idx in seq_len(nrow(robust.grid))) {
    row <- robust.grid[idx, ]
    suffix <- if (is.na(row$max.effective.match.fraction)) {
      "noef"
    } else {
      sprintf("ef%.2f", row$max.effective.match.fraction)
    }
    args <- list(
      match.threshold.rule = "local.robust.z",
      min.best.score.z = row$min.best.score.z,
      min.margin.z = row$min.margin.z
    )
    if (!is.na(row$max.effective.match.fraction)) {
      args$max.effective.match.fraction <- row$max.effective.match.fraction
    }
    policies[[length(policies) + 1L]] <- list(
      name = sprintf("robust_z_zs%.1f_zm%.1f_%s",
                     row$min.best.score.z,
                     row$min.margin.z,
                     suffix),
      label = sprintf("Robust-z zs=%.1f zm=%.1f %s",
                      row$min.best.score.z,
                      row$min.margin.z,
                      suffix),
      rule.family = "local_robust_z",
      args = args,
      description = "Adaptive local robust-z gate on best score and best-vs-second margin."
    )
  }

  policies
}

transported.hessian.adaptive.threshold.policy.table <- function(policies) {
  do.call(rbind, lapply(policies, function(policy) {
    data.frame(
      transport.filter = policy$name,
      filter.label = policy$label,
      rule.family = policy$rule.family,
      description = policy$description,
      max.match.angle = policy$args$max.match.angle %||% NA_real_,
      max.length.relative.error = policy$args$max.length.relative.error %||% NA_real_,
      match.threshold.rule = policy$args$match.threshold.rule %||% "fixed",
      match.score.quantile = policy$args$match.score.quantile %||% NA_real_,
      match.margin.quantile = policy$args$match.margin.quantile %||% NA_real_,
      min.best.score.z = policy$args$min.best.score.z %||% NA_real_,
      min.margin.z = policy$args$min.margin.z %||% NA_real_,
      max.effective.match.fraction = policy$args$max.effective.match.fraction %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

transported.hessian.drop.counts <- function(run, policy, case, method) {
  dropped <- run$dropped %||% data.frame()
  if (!nrow(dropped)) {
    return(data.frame(
      case.id = case$case.id,
      graph.type = case$graph.type,
      local.embedding.method = method,
      transport.filter = policy$name,
      drop.reason = "none",
      n = 0L,
      stringsAsFactors = FALSE
    ))
  }
  tab <- as.data.frame(table(dropped$drop.reason), stringsAsFactors = FALSE)
  names(tab) <- c("drop.reason", "n")
  tab$case.id <- case$case.id
  tab$graph.type <- case$graph.type
  tab$local.embedding.method <- method
  tab$transport.filter <- policy$name
  tab[, c("case.id", "graph.type", "local.embedding.method",
          "transport.filter", "drop.reason", "n")]
}

run.transported.hessian.adaptive.threshold.sweep <- function(mode = c("smoke", "full"),
                                                             repo.dir = "/Users/pgajer/current_projects/gflow",
                                                             output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
                                                             use.local.grip = TRUE,
                                                             methods = NULL,
                                                             min.row.retention = c(0.25, 0.50, 0.75)) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.load.packages(repo.dir, use.local.grip = use.local.grip)

  if (is.null(methods)) {
    methods <- strsplit(Sys.getenv("GFLOW_ADAPTIVE_THRESHOLD_METHODS",
                                   unset = "cmdscale"), ",", fixed = TRUE)[[1]]
    methods <- trimws(methods)
    methods <- methods[nzchar(methods)]
  }
  if (!length(methods)) {
    stop("At least one local embedding method is required.", call. = FALSE)
  }

  cases <- transported.hessian.backend.validation.cases(mode)
  policies <- transported.hessian.adaptive.threshold.policies(mode)
  policy.table <- transported.hessian.adaptive.threshold.policy.table(policies)

  summaries <- list()
  drop.counts <- list()
  counter <- 0L
  for (case in cases) {
    for (method in methods) {
      for (policy in policies) {
        counter <- counter + 1L
        run <- transported.hessian.run.one(
          case = case,
          method = method,
          filter.name = policy$name,
          filter.args = policy$args
        )
        summaries[[counter]] <- run$summary
        drop.counts[[counter]] <- transported.hessian.drop.counts(run, policy, case, method)
      }
    }
  }

  summary <- do.call(rbind, summaries)
  summary <- merge(summary, policy.table, by = "transport.filter", all.x = TRUE,
                   sort = FALSE)
  summary$degenerate <- with(summary, is.finite(row.retention) &
                               (row.retention <= 0 | n.rows <= 1))
  for (rho in min.row.retention) {
    nm <- sprintf("feasible_retention_%s", gsub("\\.", "_", sprintf("%.2f", rho)))
    summary[[nm]] <- with(summary, status == "ok" &
                            is.finite(row.retention) &
                            row.retention >= rho)
  }
  drop.count <- do.call(rbind, drop.counts)
  drop.count <- merge(drop.count, policy.table[, c("transport.filter", "rule.family")],
                      by = "transport.filter", all.x = TRUE, sort = FALSE)

  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    use.local.grip = use.local.grip,
    grip.available = requireNamespace("grip", quietly = TRUE),
    edge.kk.available = transported.hessian.edge.kk.available(),
    edge.kk.function = transported.hessian.edge.kk.function.name(),
    methods = methods,
    min.row.retention = min.row.retention,
    gflow.repo = repo.dir
  )
  out <- list(
    metadata = metadata,
    policies = policy.table,
    summary = summary,
    drop.count = drop.count
  )

  prefix <- file.path(output.dir, paste0("transported_hessian_adaptive_threshold_sweep_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(policy.table, paste0(prefix, "_policies.csv"), row.names = FALSE)
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(drop.count, paste0(prefix, "_drop_counts.csv"), row.names = FALSE)
  out
}
