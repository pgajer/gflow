source(file.path("/Users/pgajer/current_projects/gflow",
                 "dev/graph-trend-filtering/transported_hessian_backend_validation.R"))

transported.hessian.transport.filter.policies <- function() {
  list(
    list(
      name = "none",
      label = "No filter",
      args = list(),
      description = "All candidate soft-transport rows are retained."
    ),
    list(
      name = "fixed_geometric",
      label = "Fixed angle/length",
      args = list(max.match.angle = pi / 3,
                  max.length.relative.error = 0.50),
      description = "Retain rows only when the best match has angle at most pi/3 and relative length error at most 0.50."
    ),
    list(
      name = "local_quantile",
      label = "Local quantile",
      args = list(match.threshold.rule = "local.quantile",
                  match.score.quantile = 0.25,
                  match.margin.quantile = 0.50),
      description = "Retain rows only when the best score is in the local lower quartile and the best-vs-second margin is at least the local median margin."
    ),
    list(
      name = "local_robust_z",
      label = "Local robust-z",
      args = list(match.threshold.rule = "local.robust.z",
                  min.best.score.z = 1.0,
                  min.margin.z = 0.0),
      description = "Retain rows only when the best match is at least one robust scale better than the local median score and its margin is not below the local median gap."
    ),
    list(
      name = "local_robust_z_entropy",
      label = "Robust-z + entropy",
      args = list(match.threshold.rule = "local.robust.z",
                  min.best.score.z = 1.0,
                  min.margin.z = 0.0,
                  max.effective.match.fraction = 0.35),
      description = "Apply the local robust-z rule and additionally reject diffuse softmax matches by effective-match fraction."
    )
  )
}

run.transported.hessian.transport.filter.validation <- function(mode = c("smoke", "full"),
                                                                repo.dir = "/Users/pgajer/current_projects/gflow",
                                                                output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
                                                                use.local.grip = TRUE) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.load.packages(repo.dir, use.local.grip = use.local.grip)

  methods <- c("cmdscale", "mds.edge.kk", "grip.edge.kk")
  cases <- transported.hessian.backend.validation.cases(mode)
  policies <- transported.hessian.transport.filter.policies()
  runs <- list()
  counter <- 0L
  for (case in cases) {
    for (method in methods) {
      for (policy in policies) {
        counter <- counter + 1L
        runs[[counter]] <- transported.hessian.run.one(
          case = case,
          method = method,
          filter.name = policy$name,
          filter.args = policy$args
        )
      }
    }
  }

  summary <- do.call(rbind, lapply(runs, `[[`, "summary"))
  embedding <- do.call(rbind, lapply(runs, `[[`, "embedding"))
  rows <- do.call(rbind, lapply(runs, `[[`, "rows"))
  dropped <- do.call(rbind, lapply(runs, `[[`, "dropped"))

  policy.table <- do.call(rbind, lapply(policies, function(policy) {
    data.frame(
      transport.filter = policy$name,
      filter.label = policy$label,
      description = policy$description,
      max.match.angle = policy$args$max.match.angle %||% NA_real_,
      max.length.relative.error = policy$args$max.length.relative.error %||% NA_real_,
      min.match.margin = policy$args$min.match.margin %||% NA_real_,
      max.effective.matches = policy$args$max.effective.matches %||% NA_real_,
      match.threshold.rule = policy$args$match.threshold.rule %||% "fixed",
      match.score.quantile = policy$args$match.score.quantile %||% NA_real_,
      match.margin.quantile = policy$args$match.margin.quantile %||% NA_real_,
      min.best.score.z = policy$args$min.best.score.z %||% NA_real_,
      min.margin.z = policy$args$min.margin.z %||% NA_real_,
      max.effective.match.fraction = policy$args$max.effective.match.fraction %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    use.local.grip = use.local.grip,
    grip.available = requireNamespace("grip", quietly = TRUE),
    edge.kk.available = transported.hessian.edge.kk.available(),
    edge.kk.function = transported.hessian.edge.kk.function.name(),
    gflow.repo = repo.dir
  )
  out <- list(
    metadata = metadata,
    policies = policy.table,
    summary = summary,
    embedding = embedding,
    rows = rows,
    dropped = dropped
  )

  prefix <- file.path(output.dir, paste0("transported_hessian_transport_filter_validation_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(policy.table, paste0(prefix, "_policies.csv"), row.names = FALSE)
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(embedding, paste0(prefix, "_embedding.csv"), row.names = FALSE)
  utils::write.csv(rows, paste0(prefix, "_rows.csv"), row.names = FALSE)
  utils::write.csv(dropped, paste0(prefix, "_dropped.csv"), row.names = FALSE)
  out
}
