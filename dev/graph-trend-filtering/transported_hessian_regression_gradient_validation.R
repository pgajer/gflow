source(file.path("/Users/pgajer/current_projects/gflow",
                 "dev/graph-trend-filtering/transported_hessian_backend_validation.R"))

.transported.hessian.bind.rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  all.names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all.names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all.names, drop = FALSE]
  })
  do.call(rbind, rows)
}

.transported.hessian.policy.arg <- function(policy, name, default = NULL) {
  if (name %in% names(policy$args) && !is.null(policy$args[[name]])) {
    policy$args[[name]]
  } else {
    default
  }
}

.transported.hessian.policy.arg.collapse <- function(policy, name) {
  value <- .transported.hessian.policy.arg(policy, name, NULL)
  if (is.null(value)) NA_character_ else paste(value, collapse = ",")
}

transported.hessian.regression.gradient.policies <- function() {
  base.policies <- list(
    list(
      name = "exact_coordinate",
      label = "Exact coordinate",
      family = "exact",
      args = list(transport.rule = "exact.coordinate"),
      description = "Reference dart matching by exact coordinate direction labels."
    ),
    list(
      name = "soft_coordinate",
      label = "Soft coordinate",
      family = "soft",
      args = list(transport.rule = "local.embedding.soft",
                  local.embedding.method = "coordinates",
                  soft.bandwidth = 0.18),
      description = "Soft direction matching using supplied coordinates."
    ),
    list(
      name = "soft_coordinate_fixed_gate",
      label = "Soft + fixed gate",
      family = "soft",
      args = list(transport.rule = "local.embedding.soft",
                  local.embedding.method = "coordinates",
                  soft.bandwidth = 0.18,
                  max.match.angle = pi / 3,
                  max.length.relative.error = 0.50),
      description = "Soft direction matching with the conservative fixed angle/length gate."
    ),
    list(
      name = "regression_gradient_ridge0",
      label = "Regression gradient ridge=0",
      family = "regression_gradient_ambient",
      args = list(transport.rule = "regression.gradient",
                  gradient.coordinate.method = "coordinates",
                  gradient.ridge = 0),
      description = "Local weighted least-squares gradient comparison with no ridge."
    ),
    list(
      name = "regression_gradient_ridge1e-8",
      label = "Regression gradient ridge=1e-8",
      family = "regression_gradient_ambient",
      args = list(transport.rule = "regression.gradient",
                  gradient.coordinate.method = "coordinates",
                  gradient.ridge = 1e-8),
      description = "Local weighted least-squares gradient comparison with tiny ridge."
    )
  )
  local.methods <- c("cmdscale", "mds.edge.kk", "grip.edge.kk")
  local.labels <- c(
    cmdscale = "local cmdscale",
    "mds.edge.kk" = "local MDS+edge-KK",
    "grip.edge.kk" = "local GRIP+edge-KK"
  )
  local.descriptions <- c(
    cmdscale = "graph-derived cmdscale charts",
    "mds.edge.kk" = "local MDS charts refined by edge-KK when available",
    "grip.edge.kk" = "weighted-GRIP charts refined by edge-KK when available"
  )
  local.policies <- list()
  for (method in local.methods) {
    for (hops in 1:5) {
      local.policies[[length(local.policies) + 1L]] <- list(
        name = sprintf("regression_gradient_%s_hop%d",
                       gsub("\\.", "_", method), hops),
        label = sprintf("Regression gradient %s h=%d",
                        local.labels[[method]], hops),
        family = "regression_gradient_local",
        args = list(transport.rule = "regression.gradient",
                    gradient.coordinate.method = "local.embedding",
                    gradient.embedding.method = method,
                    gradient.embedding.dim = NULL,
                    gradient.disk.hops = hops,
                    gradient.max.vertices = 64L,
                    gradient.ridge = 1e-8),
        description = sprintf("Regression-gradient comparison in %s with %d-hop graph disks.",
                              local.descriptions[[method]], hops)
      )
    }
    for (radius.fraction in c(0.05, 0.075, 0.10, 0.15, 0.20)) {
      local.policies[[length(local.policies) + 1L]] <- list(
        name = sprintf("regression_gradient_%s_metric%03d",
                       gsub("\\.", "_", method),
                       as.integer(round(1000 * radius.fraction))),
        label = sprintf("Regression gradient %s metric %.3gD",
                        local.labels[[method]], radius.fraction),
        family = "regression_gradient_local",
        args = list(transport.rule = "regression.gradient",
                    gradient.coordinate.method = "local.embedding",
                    gradient.embedding.method = method,
                    gradient.embedding.dim = NULL,
                    gradient.disk.rule = "metric.diameter.fraction",
                    gradient.disk.radius.fraction = radius.fraction,
                    gradient.max.vertices = 64L,
                    gradient.ridge = 1e-8),
        description = sprintf("Regression-gradient comparison in %s with metric graph-diameter fraction %.3g.",
                              local.descriptions[[method]], radius.fraction)
      )
    }
    for (multiplier in c(1, 1.5, 2, 3, 4)) {
      local.policies[[length(local.policies) + 1L]] <- list(
        name = sprintf("regression_gradient_%s_localscale%03d",
                       gsub("\\.", "_", method),
                       as.integer(round(100 * multiplier))),
        label = sprintf("Regression gradient %s local scale %.3g",
                        local.labels[[method]], multiplier),
        family = "regression_gradient_local",
        args = list(transport.rule = "regression.gradient",
                    gradient.coordinate.method = "local.embedding",
                    gradient.embedding.method = method,
                    gradient.embedding.dim = NULL,
                    gradient.disk.rule = "metric.local.scale",
                    gradient.disk.local.scale.method = "knn.distance",
                    gradient.disk.local.scale.k = 8L,
                    gradient.disk.local.scale.multiplier = multiplier,
                    gradient.disk.min.vertices = 6L,
                    gradient.max.vertices = 64L,
                    gradient.ridge = 1e-8),
        description = sprintf("Regression-gradient comparison in %s with local kNN-distance metric scale multiplier %.3g.",
                              local.descriptions[[method]], multiplier)
      )
    }
  }
  adaptive.policy <- list(list(
    name = "regression_gradient_adaptive_cmdscale_mds_h1_5_metric_local",
    label = "Regression gradient adaptive cmdscale/MDS hop+metric+local",
    family = "regression_gradient_adaptive",
    args = list(transport.rule = "regression.gradient",
                gradient.coordinate.method = "local.embedding",
                gradient.chart.selection = "adaptive",
                gradient.embedding.candidates = c("cmdscale", "mds.edge.kk"),
                gradient.disk.rule.candidates = c("hops", "metric.diameter.fraction",
                                                   "metric.local.scale"),
                gradient.disk.hops.candidates = 1:5,
                gradient.disk.radius.fraction.candidates = c(0.05, 0.075, 0.10, 0.15, 0.20),
                gradient.disk.local.scale.method = "knn.distance",
                gradient.disk.local.scale.k = 8L,
                gradient.disk.local.scale.multiplier.candidates = c(1, 1.5, 2, 3, 4),
                gradient.disk.min.vertices = 6L,
                gradient.embedding.dim = NULL,
                gradient.max.vertices = 64L,
                gradient.ridge = 1e-8),
    description = "Adaptive graph-only selection among cmdscale and MDS+edge-KK charts with hop, metric-diameter-fraction, and local-scale metric disks."
  ))
  c(base.policies, local.policies, adaptive.policy)
}

transported.hessian.regression.gradient.policy.table <- function(policies) {
  do.call(rbind, lapply(policies, function(policy) {
    data.frame(
      transport.policy = policy$name,
      policy.label = policy$label,
      policy.family = policy$family,
      transport.rule = .transported.hessian.policy.arg(policy, "transport.rule"),
      local.embedding.method = .transported.hessian.policy.arg(policy, "local.embedding.method", NA_character_),
      soft.bandwidth = .transported.hessian.policy.arg(policy, "soft.bandwidth", NA_real_),
      max.match.angle = .transported.hessian.policy.arg(policy, "max.match.angle", NA_real_),
      max.length.relative.error = .transported.hessian.policy.arg(policy, "max.length.relative.error", NA_real_),
      gradient.coordinate.method = .transported.hessian.policy.arg(policy, "gradient.coordinate.method", NA_character_),
      gradient.embedding.method = .transported.hessian.policy.arg(policy, "gradient.embedding.method", NA_character_),
      gradient.embedding.dim = .transported.hessian.policy.arg(policy, "gradient.embedding.dim", NA_integer_),
      gradient.disk.rule = .transported.hessian.policy.arg(policy, "gradient.disk.rule", "hops"),
      gradient.disk.radius.fraction = .transported.hessian.policy.arg(policy, "gradient.disk.radius.fraction", NA_real_),
      gradient.disk.local.scale.method = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.method", NA_character_),
      gradient.disk.local.scale.k = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.k", NA_integer_),
      gradient.disk.local.scale.multiplier = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.multiplier", NA_real_),
      gradient.disk.min.vertices = .transported.hessian.policy.arg(policy, "gradient.disk.min.vertices", NA_integer_),
      gradient.disk.hops = .transported.hessian.policy.arg(policy, "gradient.disk.hops", NA_integer_),
      gradient.chart.selection = .transported.hessian.policy.arg(policy, "gradient.chart.selection", "fixed"),
      gradient.embedding.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.embedding.candidates"),
      gradient.disk.rule.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.rule.candidates"),
      gradient.disk.hops.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.hops.candidates"),
      gradient.disk.radius.fraction.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.radius.fraction.candidates"),
      gradient.disk.local.scale.multiplier.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.local.scale.multiplier.candidates"),
      gradient.max.vertices = .transported.hessian.policy.arg(policy, "gradient.max.vertices", NA_integer_),
      gradient.ridge = .transported.hessian.policy.arg(policy, "gradient.ridge", NA_real_),
      description = policy$description,
      stringsAsFactors = FALSE
    )
  }))
}

transported.hessian.regression.gradient.summary <- function(op, case, policy,
                                                            elapsed.sec,
                                                            status = "ok",
                                                            message = NA_character_) {
  row.tab <- op$row.table %||% data.frame()
  unique.rows <- if (nrow(row.tab) && "row" %in% names(row.tab)) {
    row.tab[!duplicated(row.tab$row), , drop = FALSE]
  } else {
    data.frame()
  }
  grad <- op$transport$gradient.diagnostics %||% data.frame()
  embedding <- op$transport$embedding.table %||% data.frame()
  embedding.summary <- if (nrow(embedding) && "chart.selected" %in% names(embedding)) {
    subset(embedding, isTRUE(chart.selected) | chart.selected)
  } else {
    embedding
  }
  residuals <- transported.hessian.residual.summary(op)
  data.frame(
    case.id = case$case.id,
    graph.type = case$graph.type,
    n.vertices = length(case$adj.list),
    expected.dimension = case$expected.dimension,
    transport.policy = policy$name,
    policy.label = policy$label,
    policy.family = policy$family,
    transport.rule = .transported.hessian.policy.arg(policy, "transport.rule"),
    gradient.coordinate.method = .transported.hessian.policy.arg(policy, "gradient.coordinate.method", NA_character_),
    gradient.embedding.method = .transported.hessian.policy.arg(policy, "gradient.embedding.method", NA_character_),
    gradient.disk.rule = .transported.hessian.policy.arg(policy, "gradient.disk.rule", "hops"),
    gradient.disk.radius.fraction = .transported.hessian.policy.arg(policy, "gradient.disk.radius.fraction", NA_real_),
    gradient.disk.local.scale.method = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.method", NA_character_),
    gradient.disk.local.scale.k = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.k", NA_integer_),
    gradient.disk.local.scale.multiplier = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.multiplier", NA_real_),
    gradient.disk.min.vertices = .transported.hessian.policy.arg(policy, "gradient.disk.min.vertices", NA_integer_),
    gradient.disk.hops = .transported.hessian.policy.arg(policy, "gradient.disk.hops", NA_integer_),
    gradient.chart.selection = .transported.hessian.policy.arg(policy, "gradient.chart.selection", "fixed"),
    gradient.embedding.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.embedding.candidates"),
    gradient.disk.rule.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.rule.candidates"),
    gradient.disk.hops.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.hops.candidates"),
    gradient.disk.radius.fraction.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.radius.fraction.candidates"),
    gradient.disk.local.scale.multiplier.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.local.scale.multiplier.candidates"),
    gradient.max.vertices = .transported.hessian.policy.arg(policy, "gradient.max.vertices", NA_integer_),
    embedding.backend = if (nrow(embedding.summary) && "backend.used" %in% names(embedding.summary)) {
      paste(unique(embedding.summary$backend.used), collapse = "; ")
    } else {
      NA_character_
    },
    status = status,
    runtime.sec = elapsed.sec,
    nullity.estimate = op$diagnostics$nullity.estimate,
    n.candidate.rows = op$diagnostics$summary$n.candidate.rows,
    n.rows = op$diagnostics$summary$n.rows,
    n.row.records = op$diagnostics$summary$n.row.records,
    n.dropped = op$diagnostics$summary$n.dropped,
    row.retention = if (op$diagnostics$summary$n.candidate.rows > 0) {
      op$diagnostics$summary$n.rows / op$diagnostics$summary$n.candidate.rows
    } else {
      NA_real_
    },
    mean.transport.entropy = if (nrow(unique.rows) &&
                                 "transport.entropy" %in% names(unique.rows)) {
      mean(unique.rows$transport.entropy, na.rm = TRUE)
    } else {
      NA_real_
    },
    mean.effective.match.fraction = if (nrow(unique.rows) &&
                                         "effective.match.fraction" %in% names(unique.rows)) {
      mean(unique.rows$effective.match.fraction, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.local.rank = if (nrow(grad)) stats::median(grad$rank, na.rm = TRUE) else NA_real_,
    min.local.rank = if (nrow(grad)) min(grad$rank, na.rm = TRUE) else NA_real_,
    rank.deficient.fraction = if (nrow(grad)) {
      mean(grad$rank < case$expected.dimension, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.local.condition = if (nrow(grad)) {
      finite <- grad$condition[is.finite(grad$condition)]
      if (length(finite)) stats::median(finite) else NA_real_
    } else {
      NA_real_
    },
    max.local.condition = if (nrow(grad)) {
      finite <- grad$condition[is.finite(grad$condition)]
      if (length(finite)) max(finite) else NA_real_
    } else {
      NA_real_
    },
    median.local.neighborhood = if (nrow(grad)) {
      stats::median(grad$n.local, na.rm = TRUE)
    } else {
      NA_real_
    },
    n.chart.candidates = nrow(embedding),
    n.selected.charts = if (nrow(embedding.summary)) nrow(embedding.summary) else NA_integer_,
    selected.oracle.fraction = if (nrow(embedding.summary) &&
                                   "chart.oracle" %in% names(embedding.summary)) {
      mean(embedding.summary$chart.oracle, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.disk.size = if (nrow(embedding.summary) &&
                           "disk.size" %in% names(embedding.summary)) {
      stats::median(embedding.summary$disk.size, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.disk.fraction = if (nrow(embedding.summary) &&
                               "disk.fraction" %in% names(embedding.summary)) {
      stats::median(embedding.summary$disk.fraction, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.disk.weighted.radius = if (nrow(embedding.summary) &&
                                      "disk.weighted.radius" %in% names(embedding.summary)) {
      stats::median(embedding.summary$disk.weighted.radius, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.embedding.edge.stress = if (nrow(embedding.summary) &&
                                       "edge.stress" %in% names(embedding.summary)) {
      stats::median(embedding.summary$edge.stress, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.embedding.distance.stress = if (nrow(embedding.summary) &&
                                           "graph.distance.stress" %in% names(embedding.summary)) {
      stats::median(embedding.summary$graph.distance.stress, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.ambient.affine.residual = if (nrow(embedding.summary) &&
                                          "ambient.affine.residual" %in% names(embedding.summary)) {
      stats::median(embedding.summary$ambient.affine.residual, na.rm = TRUE)
    } else {
      NA_real_
    },
    median.chart.linear.gradient.residual = if (nrow(embedding.summary) &&
                                                "ambient.linear.gradient.residual.mean" %in%
                                                names(embedding.summary)) {
      stats::median(embedding.summary$ambient.linear.gradient.residual.mean, na.rm = TRUE)
    } else {
      NA_real_
    },
    max.chart.linear.gradient.residual = if (nrow(embedding.summary) &&
                                             "ambient.linear.gradient.residual.max" %in%
                                             names(embedding.summary)) {
      max(embedding.summary$ambient.linear.gradient.residual.max, na.rm = TRUE)
    } else {
      NA_real_
    },
    message = message,
    residual.constant = residuals$residual.constant,
    residual.linear.mean = residuals$residual.linear.mean,
    residual.quadratic.mean = residuals$residual.quadratic.mean,
    residual.overall = residuals$residual.overall,
    stringsAsFactors = FALSE
  )
}

transported.hessian.regression.gradient.error.summary <- function(case, policy,
                                                                  elapsed.sec,
                                                                  message) {
  data.frame(
    case.id = case$case.id,
    graph.type = case$graph.type,
    n.vertices = length(case$adj.list),
    expected.dimension = case$expected.dimension,
    transport.policy = policy$name,
    policy.label = policy$label,
    policy.family = policy$family,
    transport.rule = .transported.hessian.policy.arg(policy, "transport.rule"),
    gradient.coordinate.method = .transported.hessian.policy.arg(policy, "gradient.coordinate.method", NA_character_),
    gradient.embedding.method = .transported.hessian.policy.arg(policy, "gradient.embedding.method", NA_character_),
    gradient.disk.rule = .transported.hessian.policy.arg(policy, "gradient.disk.rule", "hops"),
    gradient.disk.radius.fraction = .transported.hessian.policy.arg(policy, "gradient.disk.radius.fraction", NA_real_),
    gradient.disk.local.scale.method = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.method", NA_character_),
    gradient.disk.local.scale.k = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.k", NA_integer_),
    gradient.disk.local.scale.multiplier = .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.multiplier", NA_real_),
    gradient.disk.min.vertices = .transported.hessian.policy.arg(policy, "gradient.disk.min.vertices", NA_integer_),
    gradient.disk.hops = .transported.hessian.policy.arg(policy, "gradient.disk.hops", NA_integer_),
    gradient.chart.selection = .transported.hessian.policy.arg(policy, "gradient.chart.selection", "fixed"),
    gradient.embedding.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.embedding.candidates"),
    gradient.disk.rule.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.rule.candidates"),
    gradient.disk.hops.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.hops.candidates"),
    gradient.disk.radius.fraction.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.radius.fraction.candidates"),
    gradient.disk.local.scale.multiplier.candidates = .transported.hessian.policy.arg.collapse(policy, "gradient.disk.local.scale.multiplier.candidates"),
    gradient.max.vertices = .transported.hessian.policy.arg(policy, "gradient.max.vertices", NA_integer_),
    embedding.backend = NA_character_,
    status = "error",
    runtime.sec = elapsed.sec,
    nullity.estimate = NA_integer_,
    n.candidate.rows = NA_integer_,
    n.rows = NA_integer_,
    n.row.records = NA_integer_,
    n.dropped = NA_integer_,
    row.retention = NA_real_,
    mean.transport.entropy = NA_real_,
    mean.effective.match.fraction = NA_real_,
    median.local.rank = NA_real_,
    min.local.rank = NA_real_,
    rank.deficient.fraction = NA_real_,
    median.local.condition = NA_real_,
    max.local.condition = NA_real_,
    median.local.neighborhood = NA_real_,
    n.chart.candidates = NA_integer_,
    n.selected.charts = NA_integer_,
    selected.oracle.fraction = NA_real_,
    median.disk.size = NA_real_,
    median.disk.fraction = NA_real_,
    median.disk.weighted.radius = NA_real_,
    median.embedding.edge.stress = NA_real_,
    median.embedding.distance.stress = NA_real_,
    median.ambient.affine.residual = NA_real_,
    median.chart.linear.gradient.residual = NA_real_,
    max.chart.linear.gradient.residual = NA_real_,
    message = message,
    residual.constant = NA_real_,
    residual.linear.mean = NA_real_,
    residual.quadratic.mean = NA_real_,
    residual.overall = NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.run.regression.gradient.one <- function(case, policy) {
  probes <- transported.hessian.polynomial.probes(case$coordinates)
  op.args <- c(
    list(
      adj.list = case$adj.list,
      weight.list = case$weight.list,
      coordinates = case$coordinates,
      polynomial.probes = probes
    ),
    policy$args
  )
  start <- proc.time()[["elapsed"]]
  fit <- try(do.call(transported.graph.hessian.operator, op.args), silent = TRUE)
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(fit, "try-error")) {
    return(list(
      summary = transported.hessian.regression.gradient.error.summary(
        case, policy, elapsed, as.character(fit)
      ),
      gradient = data.frame(),
      embedding = data.frame()
    ))
  }
  summary <- transported.hessian.regression.gradient.summary(
    op = fit,
    case = case,
    policy = policy,
    elapsed.sec = elapsed
  )
  grad <- fit$transport$gradient.diagnostics %||% data.frame()
  if (nrow(grad)) {
    grad$case.id <- case$case.id
    grad$graph.type <- case$graph.type
    grad$transport.policy <- policy$name
    grad$policy.label <- policy$label
    grad$policy.family <- policy$family
    if (!"gradient.disk.rule" %in% names(grad)) {
      grad$gradient.disk.rule <- .transported.hessian.policy.arg(policy, "gradient.disk.rule", "hops")
    }
    if (!"gradient.disk.hops" %in% names(grad)) {
      grad$gradient.disk.hops <- .transported.hessian.policy.arg(policy, "gradient.disk.hops", NA_integer_)
    }
    if (!"gradient.disk.radius.fraction" %in% names(grad)) {
      grad$gradient.disk.radius.fraction <- .transported.hessian.policy.arg(policy, "gradient.disk.radius.fraction", NA_real_)
    }
    if (!"gradient.disk.local.scale.multiplier" %in% names(grad)) {
      grad$gradient.disk.local.scale.multiplier <- .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.multiplier", NA_real_)
    }
    if (!"gradient.disk.local.scale" %in% names(grad)) {
      grad$gradient.disk.local.scale <- NA_real_
    }
  }
  embedding <- fit$transport$embedding.table %||% data.frame()
  if (nrow(embedding)) {
    embedding$case.id <- case$case.id
    embedding$graph.type <- case$graph.type
    embedding$transport.policy <- policy$name
    embedding$policy.label <- policy$label
    embedding$policy.family <- policy$family
    if (!"gradient.disk.rule" %in% names(embedding)) {
      embedding$gradient.disk.rule <- .transported.hessian.policy.arg(policy, "gradient.disk.rule", "hops")
    }
    if (!"gradient.disk.hops" %in% names(embedding)) {
      embedding$gradient.disk.hops <- .transported.hessian.policy.arg(policy, "gradient.disk.hops", NA_integer_)
    }
    if (!"gradient.disk.radius.fraction" %in% names(embedding)) {
      embedding$gradient.disk.radius.fraction <- .transported.hessian.policy.arg(policy, "gradient.disk.radius.fraction", NA_real_)
    }
    if (!"gradient.disk.local.scale.multiplier" %in% names(embedding)) {
      embedding$gradient.disk.local.scale.multiplier <- .transported.hessian.policy.arg(policy, "gradient.disk.local.scale.multiplier", NA_real_)
    }
    if (!"gradient.disk.local.scale" %in% names(embedding)) {
      embedding$gradient.disk.local.scale <- NA_real_
    }
  }
  list(summary = summary, gradient = grad, embedding = embedding)
}

run.transported.hessian.regression.gradient.validation <- function(mode = c("smoke", "full"),
                                                                   repo.dir = "/Users/pgajer/current_projects/gflow",
                                                                   output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
                                                                   use.local.grip = TRUE) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.load.packages(repo.dir, use.local.grip = use.local.grip)

  cases <- transported.hessian.backend.validation.cases(mode)
  policies <- transported.hessian.regression.gradient.policies()
  policy.table <- transported.hessian.regression.gradient.policy.table(policies)
  runs <- list()
  counter <- 0L
  for (case in cases) {
    for (policy in policies) {
      counter <- counter + 1L
      runs[[counter]] <- transported.hessian.run.regression.gradient.one(case, policy)
    }
  }

  summary <- .transported.hessian.bind.rows(lapply(runs, `[[`, "summary"))
  gradient <- .transported.hessian.bind.rows(lapply(runs, `[[`, "gradient"))
  embedding <- .transported.hessian.bind.rows(lapply(runs, `[[`, "embedding"))
  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    use.local.grip = use.local.grip,
    gflow.repo = repo.dir,
    n.policies = length(policies)
  )
  out <- list(
    metadata = metadata,
    policies = policy.table,
    summary = summary,
    gradient = gradient,
    embedding = embedding
  )

  prefix <- file.path(output.dir, paste0("transported_hessian_regression_gradient_validation_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(policy.table, paste0(prefix, "_policies.csv"), row.names = FALSE)
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(gradient, paste0(prefix, "_gradient.csv"), row.names = FALSE)
  utils::write.csv(embedding, paste0(prefix, "_embedding.csv"), row.names = FALSE)
  out
}
