repo.dir.default <- "/Users/pgajer/current_projects/gflow"

`%||%` <- function(x, y) if (is.null(x)) y else x

source(file.path(repo.dir.default,
                 "dev/graph-trend-filtering/transported_hessian_quadric_coordinate_validation.R"))

transported.hessian.quadric.polynomial.probes <- function(U,
                                                          fractional.alpha = c(1.25,
                                                                               1.50,
                                                                               1.75,
                                                                               2.25,
                                                                               2.50,
                                                                               3.50)) {
  U <- as.matrix(U)
  n <- nrow(U)
  out <- list(constant = rep(1, n))
  for (j in seq_len(ncol(U))) {
    out[[paste0("affine_u", j)]] <- U[, j]
  }
  for (i in seq_len(ncol(U))) {
    for (j in i:ncol(U)) {
      out[[paste0("quadratic_u", i, "u", j)]] <- U[, i] * U[, j]
    }
  }
  for (i in seq_len(ncol(U))) {
    for (j in i:ncol(U)) {
      for (k in j:ncol(U)) {
        out[[paste0("cubic_u", i, "u", j, "u", k)]] <- U[, i] * U[, j] * U[, k]
      }
    }
  }
  centered <- sweep(U, 2L, colMeans(U), "-")
  radius <- sqrt(rowSums(centered^2))
  radius <- radius / max(radius, sqrt(.Machine$double.eps))
  for (alpha in fractional.alpha) {
    out[[sprintf("fractional_radial_%.2f", alpha)]] <- radius^alpha
  }
  as.data.frame(out, check.names = FALSE)
}

transported.hessian.quadric.multimethod.policies <- function() {
  list(
    list(
      method.id = "regression_gradient_intrinsic",
      method.label = "Regression gradient: intrinsic oracle",
      method.family = "regression.gradient",
      args = function(case) list(
        transport.rule = "regression.gradient",
        coordinates = case$U,
        gradient.coordinate.method = "coordinates",
        gradient.ridge = 1e-8
      )
    ),
    list(
      method.id = "regression_gradient_ambient",
      method.label = "Regression gradient: data coordinates",
      method.family = "regression.gradient",
      args = function(case) list(
        transport.rule = "regression.gradient",
        coordinates = case$X,
        gradient.coordinate.method = "coordinates",
        gradient.ridge = 1e-8
      )
    ),
    list(
      method.id = "regression_gradient_pca",
      method.label = "Regression gradient: PCA chart",
      method.family = "regression.gradient",
      args = function(case) list(
        transport.rule = "regression.gradient",
        coordinates = case$pca,
        gradient.coordinate.method = "coordinates",
        gradient.ridge = 1e-8
      )
    ),
    list(
      method.id = "local_embedding_soft_coordinates",
      method.label = "Soft direction matching",
      method.family = "local.embedding.soft",
      args = function(case) list(
        transport.rule = "local.embedding.soft",
        local.embedding.method = "coordinates",
        coordinates = case$X,
        soft.bandwidth = 0.25
      )
    ),
    list(
      method.id = "local_embedding_soft_gated",
      method.label = "Soft direction matching: gated",
      method.family = "local.embedding.soft",
      args = function(case) list(
        transport.rule = "local.embedding.soft",
        local.embedding.method = "coordinates",
        coordinates = case$X,
        soft.bandwidth = 0.25,
        max.match.angle = pi / 3,
        max.length.relative.error = 0.75
      )
    ),
    list(
      method.id = "edge_angle_hard",
      method.label = "Edge-angle hard matching",
      method.family = "edge.angle",
      args = function(case) list(
        transport.rule = "edge.angle.hard",
        coordinates = case$X,
        edge.angle.max.angle.difference = pi / 4,
        edge.angle.max.length.relative.error = 0.75
      )
    ),
    list(
      method.id = "edge_angle_soft",
      method.label = "Edge-angle soft matching",
      method.family = "edge.angle",
      args = function(case) list(
        transport.rule = "edge.angle.soft",
        coordinates = case$X,
        edge.angle.bandwidth = 0.35
      )
    )
  )
}

transported.hessian.quadric.multimethod.cases <- function(mode = c("smoke",
                                                                    "full")) {
  mode <- match.arg(mode)
  if (identical(mode, "smoke")) {
    return(list(
      transported.hessian.quadric.case(2L, "elliptic", 0.75, n = 36L,
                                       seed = 11L, k = 6L),
      transported.hessian.quadric.case(3L, "saddle", 1.50, n = 44L,
                                       seed = 22L, k = 9L)
    ))
  }
  cases <- list()
  idx <- 0L
  for (dim in c(2L, 3L)) {
    n <- if (dim == 2L) 44L else 54L
    k <- if (dim == 2L) 7L else 10L
    for (surface in c("elliptic", "saddle")) {
      for (curvature in c(0.25, 1.50)) {
        idx <- idx + 1L
        cases[[idx]] <- transported.hessian.quadric.case(
          dim = dim,
          surface = surface,
          curvature = curvature,
          n = n,
          seed = 501L + 17L * idx,
          k = k
        )
      }
    }
  }
  cases
}

transported.hessian.quadric.probe.summary <- function(op) {
  tab <- op$diagnostics$polynomial.residuals$per.column
  residual <- stats::setNames(tab$residual, tab$probe)
  summarize <- function(prefix) {
    idx <- grep(paste0("^", prefix), names(residual), value = TRUE)
    if (!length(idx)) {
      return(c(mean = NA_real_, min = NA_real_, max = NA_real_))
    }
    c(mean = mean(residual[idx], na.rm = TRUE),
      min = min(residual[idx], na.rm = TRUE),
      max = max(residual[idx], na.rm = TRUE))
  }
  affine <- summarize("affine_")
  quadratic <- summarize("quadratic_")
  cubic <- summarize("cubic_")
  fractional <- grep("^fractional_", names(residual), value = TRUE)
  fractional.alpha <- as.numeric(sub("^fractional_radial_", "", fractional))
  best.frac <- if (length(fractional)) {
    idx <- which.min(residual[fractional])
    fractional[idx]
  } else {
    NA_character_
  }
  data.frame(
    residual.constant = unname(residual[["constant"]] %||% NA_real_),
    residual.affine.mean = affine[["mean"]],
    residual.affine.max = affine[["max"]],
    residual.quadratic.mean = quadratic[["mean"]],
    residual.quadratic.min = quadratic[["min"]],
    residual.cubic.mean = cubic[["mean"]],
    residual.cubic.min = cubic[["min"]],
    residual.fractional.mean = if (length(fractional)) {
      mean(residual[fractional], na.rm = TRUE)
    } else {
      NA_real_
    },
    residual.fractional.min = if (length(fractional)) {
      min(residual[fractional], na.rm = TRUE)
    } else {
      NA_real_
    },
    best.fractional.probe = best.frac,
    best.fractional.alpha = if (length(fractional)) {
      fractional.alpha[which.min(residual[fractional])]
    } else {
      NA_real_
    },
    residual.overall = op$diagnostics$polynomial.residuals$overall %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.method.fit.summary <- function(op,
                                                           case,
                                                           policy,
                                                           graph.weight,
                                                           elapsed.sec) {
  summ <- op$diagnostics$summary
  grad <- op$transport$gradient.diagnostics %||% data.frame()
  probe.summary <- transported.hessian.quadric.probe.summary(op)
  data.frame(
    case.id = case$case.id,
    intrinsic.dim = case$intrinsic.dim,
    surface = case$surface,
    curvature = case$curvature,
    n = case$n,
    seed = case$seed,
    k = case$k,
    graph.weight = graph.weight,
    method.id = policy$method.id,
    method.label = policy$method.label,
    method.family = policy$method.family,
    status = "ok",
    runtime.sec = elapsed.sec,
    n.edges = case$graph.metadata$n.edges,
    n.mst.edges.added = case$graph.metadata$n_mst_edges_added %||% case$graph.metadata$n.mst.edges.added,
    n.rows = summ$n.rows,
    n.row.records = summ$n.row.records,
    n.dropped = summ$n.dropped,
    row.retention = if (summ$n.candidate.rows > 0) {
      summ$n.rows / summ$n.candidate.rows
    } else {
      NA_real_
    },
    mean.transport.entropy = summ$mean.transport.entropy,
    mean.match.score = summ$mean.match.score,
    nullity.estimate = op$diagnostics$nullity.estimate,
    median.local.condition = if (nrow(grad) && "condition" %in% names(grad)) {
      finite <- grad$condition[is.finite(grad$condition)]
      if (length(finite)) stats::median(finite) else NA_real_
    } else {
      NA_real_
    },
    message = NA_character_,
    probe.summary,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.method.error.summary <- function(case,
                                                             policy,
                                                             graph.weight,
                                                             elapsed.sec,
                                                             message) {
  data.frame(
    case.id = case$case.id,
    intrinsic.dim = case$intrinsic.dim,
    surface = case$surface,
    curvature = case$curvature,
    n = case$n,
    seed = case$seed,
    k = case$k,
    graph.weight = graph.weight,
    method.id = policy$method.id,
    method.label = policy$method.label,
    method.family = policy$method.family,
    status = "error",
    runtime.sec = elapsed.sec,
    n.edges = case$graph.metadata$n.edges,
    n.mst.edges.added = case$graph.metadata$n.mst.edges.added,
    n.rows = NA_integer_,
    n.row.records = NA_integer_,
    n.dropped = NA_integer_,
    row.retention = NA_real_,
    mean.transport.entropy = NA_real_,
    mean.match.score = NA_real_,
    nullity.estimate = NA_integer_,
    median.local.condition = NA_real_,
    message = message,
    residual.constant = NA_real_,
    residual.affine.mean = NA_real_,
    residual.affine.max = NA_real_,
    residual.quadratic.mean = NA_real_,
    residual.quadratic.min = NA_real_,
    residual.cubic.mean = NA_real_,
    residual.cubic.min = NA_real_,
    residual.fractional.mean = NA_real_,
    residual.fractional.min = NA_real_,
    best.fractional.probe = NA_character_,
    best.fractional.alpha = NA_real_,
    residual.overall = NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.probe.table <- function(op, case, policy, graph.weight) {
  tab <- op$diagnostics$polynomial.residuals$per.column
  if (!nrow(tab)) return(data.frame())
  tab$case.id <- case$case.id
  tab$intrinsic.dim <- case$intrinsic.dim
  tab$surface <- case$surface
  tab$curvature <- case$curvature
  tab$graph.weight <- graph.weight
  tab$method.id <- policy$method.id
  tab$method.label <- policy$method.label
  tab$method.family <- policy$method.family
  tab$probe.family <- sub("_.*$", "", tab$probe)
  tab
}

transported.hessian.quadric.method.run.one <- function(case,
                                                       policy,
                                                       graph.weight = "ambient") {
  probes <- as.matrix(transported.hessian.quadric.polynomial.probes(case$U))
  weight.list <- transported.hessian.quadric.weight.list(case, graph.weight)
  args <- policy$args(case)
  start <- proc.time()[["elapsed"]]
  fit <- try(
    do.call(
      transported.graph.hessian.operator,
      c(list(adj.list = case$adj.list,
             weight.list = weight.list,
             polynomial.probes = probes),
        args)
    ),
    silent = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(fit, "try-error")) {
    return(list(
      summary = transported.hessian.quadric.method.error.summary(
        case, policy, graph.weight, elapsed, as.character(fit)
      ),
      probes = data.frame()
    ))
  }
  list(
    summary = transported.hessian.quadric.method.fit.summary(
      fit, case, policy, graph.weight, elapsed
    ),
    probes = transported.hessian.quadric.probe.table(fit, case, policy, graph.weight)
  )
}

run.transported.hessian.quadric.multimethod.validation <- function(
    mode = c("smoke", "full"),
    repo.dir = repo.dir.default,
    output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
    graph.weights = "ambient") {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.quadric.load.packages(repo.dir)
  cases <- transported.hessian.quadric.multimethod.cases(mode)
  policies <- transported.hessian.quadric.multimethod.policies()
  runs <- list()
  counter <- 0L
  for (case in cases) {
    for (graph.weight in graph.weights) {
      for (policy in policies) {
        counter <- counter + 1L
        message(sprintf("[%s %d/%d] %s | %s | %s",
                        mode,
                        counter,
                        length(cases) * length(graph.weights) * length(policies),
                        case$case.id,
                        graph.weight,
                        policy$method.id))
        runs[[counter]] <- transported.hessian.quadric.method.run.one(
          case = case,
          policy = policy,
          graph.weight = graph.weight
        )
      }
    }
  }
  summary <- .transported.hessian.bind.rows(lapply(runs, `[[`, "summary"))
  probes <- .transported.hessian.bind.rows(lapply(runs, `[[`, "probes"))
  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    gflow.repo = repo.dir,
    n.cases = length(cases),
    n.runs = nrow(summary),
    graph.weights = graph.weights
  )
  out <- list(metadata = metadata, summary = summary, probes = probes)
  prefix <- file.path(output.dir,
                      paste0("transported_hessian_quadric_multimethod_validation_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(probes, paste0(prefix, "_probes.csv"), row.names = FALSE)
  out
}
