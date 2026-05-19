repo.dir.default <- "/Users/pgajer/current_projects/gflow"

`%||%` <- function(x, y) if (is.null(x)) y else x

transported.hessian.quadric.load.packages <- function(repo.dir = repo.dir.default) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required for this validation.", call. = FALSE)
  }
  devtools::load_all(repo.dir, quiet = TRUE)
}

transported.hessian.quadric.sample.parameters <- function(dim,
                                                          n,
                                                          domain.radius,
                                                          domain.shape,
                                                          sample.method,
                                                          seed) {
  set.seed(seed)
  if (dim == 2L) {
    sampler <- get(".quadform.sample.parameter", envir = asNamespace("gflow"))
    return(sampler(n, domain.radius, sample.method))
  }
  sampler <- get(".quadform.sample.parameter.3d", envir = asNamespace("gflow"))
  sampler(n, domain.radius, domain.shape)
}

transported.hessian.quadric.intrinsic.probes <- function(U) {
  U <- as.matrix(U)
  n <- nrow(U)
  out <- list(constant = rep(1, n))
  for (j in seq_len(ncol(U))) {
    out[[paste0("u", j)]] <- U[, j]
  }
  for (i in seq_len(ncol(U))) {
    for (j in i:ncol(U)) {
      out[[paste0("u", i, "u", j)]] <- U[, i] * U[, j]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

transported.hessian.quadric.residual.summary <- function(op) {
  tab <- op$diagnostics$polynomial.residuals$per.column
  named <- stats::setNames(tab$residual, tab$probe)
  linear <- grep("^u[0-9]+$", names(named), value = TRUE)
  quadratic <- grep("^u[0-9]+u[0-9]+$", names(named), value = TRUE)
  data.frame(
    residual.constant = unname(named[["constant"]] %||% NA_real_),
    residual.linear.mean = if (length(linear)) mean(named[linear], na.rm = TRUE) else NA_real_,
    residual.linear.max = if (length(linear)) max(named[linear], na.rm = TRUE) else NA_real_,
    residual.quadratic.mean = if (length(quadratic)) mean(named[quadratic], na.rm = TRUE) else NA_real_,
    residual.quadratic.min = if (length(quadratic)) min(named[quadratic], na.rm = TRUE) else NA_real_,
    residual.overall = op$diagnostics$polynomial.residuals$overall %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.graph.surface.weights <- function(graph,
                                                              U,
                                                              index.k,
                                                              coefficients) {
  weight.list <- vector("list", length(graph$adj_list))
  for (i in seq_along(graph$adj_list)) {
    nbrs <- graph$adj_list[[i]]
    if (!length(nbrs)) {
      weight.list[[i]] <- numeric()
      next
    }
    weight.list[[i]] <- quadform.edge.lengths(
      U[rep.int(i, length(nbrs)), , drop = FALSE],
      U[nbrs, , drop = FALSE],
      index.k = index.k,
      coefficients = coefficients
    )
  }
  weight.list
}

transported.hessian.quadric.case <- function(dim,
                                             surface,
                                             curvature,
                                             n,
                                             seed,
                                             k,
                                             sample.method = NULL,
                                             domain.shape = NULL,
                                             domain.radius = 1) {
  surface <- match.arg(surface, c("elliptic", "saddle"))
  if (dim == 2L) {
    sample.method <- sample.method %||% "uniform.parameter.square"
    domain.shape <- if (grepl("\\.disk$", sample.method)) "disk" else "square"
  } else if (dim == 3L) {
    domain.shape <- domain.shape %||% "cube"
    sample.method <- domain.shape
  } else {
    stop("Only intrinsic dimensions 2 and 3 are supported.", call. = FALSE)
  }
  index.k <- if (identical(surface, "elliptic")) dim else max(1L, dim - 1L)
  coefficients <- rep(as.double(curvature), dim)
  U <- transported.hessian.quadric.sample.parameters(
    dim = dim,
    n = n,
    domain.radius = domain.radius,
    domain.shape = domain.shape,
    sample.method = sample.method,
    seed = seed
  )
  colnames(U) <- paste0("u", seq_len(dim))
  X <- quadform.embed(U, index.k = index.k, coefficients = coefficients)
  graph <- create.sknn.graph(X, k = k, connect.components = TRUE)
  surface.weight.list <- transported.hessian.quadric.graph.surface.weights(
    graph = graph,
    U = U,
    index.k = index.k,
    coefficients = coefficients
  )
  pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  pca.coords <- pca$x[, seq_len(dim), drop = FALSE]
  colnames(pca.coords) <- paste0("pc", seq_len(dim))
  list(
    case.id = sprintf("quadric_%dd_%s_c%.2f_n%d_seed%d",
                      dim, surface, curvature, n, seed),
    intrinsic.dim = dim,
    surface = surface,
    curvature = curvature,
    n = n,
    seed = seed,
    k = k,
    index.k = index.k,
    coefficients = coefficients,
    sample.method = sample.method,
    domain.shape = domain.shape,
    U = U,
    X = X,
    pca = pca.coords,
    adj.list = graph$adj_list,
    weight.list.ambient = graph$weight_list,
    weight.list.surface = surface.weight.list,
    graph.metadata = list(
      n.components.before = graph$n_components_before,
      n.components.after = graph$n_components_after,
      n.mst.edges.added = graph$n_mst_edges_added,
      n.edges = sum(lengths(graph$adj_list)) / 2
    )
  )
}

transported.hessian.quadric.validation.cases <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  if (identical(mode, "smoke")) {
    return(list(
      transported.hessian.quadric.case(2L, "elliptic", 0.75, n = 40L,
                                       seed = 11L, k = 6L),
      transported.hessian.quadric.case(2L, "saddle", 1.50, n = 40L,
                                       seed = 12L, k = 6L),
      transported.hessian.quadric.case(3L, "elliptic", 0.75, n = 50L,
                                       seed = 21L, k = 9L),
      transported.hessian.quadric.case(3L, "saddle", 1.50, n = 50L,
                                       seed = 22L, k = 9L)
    ))
  }

  cases <- list()
  idx <- 0L
  for (dim in c(2L, 3L)) {
    n <- if (dim == 2L) 70L else 90L
    k <- if (dim == 2L) 8L else 12L
    for (surface in c("elliptic", "saddle")) {
      for (curvature in c(0.25, 0.75, 1.50, 3.00)) {
        for (seed in c(101L, 202L)) {
          idx <- idx + 1L
          cases[[idx]] <- transported.hessian.quadric.case(
            dim = dim,
            surface = surface,
            curvature = curvature,
            n = n,
            seed = seed + 10L * dim + if (identical(surface, "saddle")) 1000L else 0L,
            k = k
          )
        }
      }
    }
  }
  cases
}

transported.hessian.quadric.coordinate.policies <- function() {
  list(
    list(
      coordinate.source = "intrinsic",
      label = "Intrinsic chart coordinates",
      description = "Oracle coordinates equal to the intrinsic quadric parameters."
    ),
    list(
      coordinate.source = "ambient",
      label = "Ambient data coordinates",
      description = "Observed data matrix coordinates on the quadratic graph hypersurface."
    ),
    list(
      coordinate.source = "pca_intrinsic_dim",
      label = "PCA to intrinsic dimension",
      description = "PCA projection of ambient data to the known intrinsic dimension."
    )
  )
}

transported.hessian.quadric.coordinates <- function(case, coordinate.source) {
  switch(
    coordinate.source,
    intrinsic = case$U,
    ambient = case$X,
    pca_intrinsic_dim = case$pca,
    stop("Unknown coordinate source: ", coordinate.source, call. = FALSE)
  )
}

transported.hessian.quadric.weight.list <- function(case, graph.weight) {
  switch(
    graph.weight,
    ambient = case$weight.list.ambient,
    surface.segment = case$weight.list.surface,
    stop("Unknown graph weight rule: ", graph.weight, call. = FALSE)
  )
}

transported.hessian.quadric.run.one <- function(case,
                                                policy,
                                                graph.weight,
                                                ridge = 1e-8) {
  probes <- as.matrix(transported.hessian.quadric.intrinsic.probes(case$U))
  coords <- transported.hessian.quadric.coordinates(case, policy$coordinate.source)
  weight.list <- transported.hessian.quadric.weight.list(case, graph.weight)
  start <- proc.time()[["elapsed"]]
  fit <- try(
    transported.graph.hessian.operator(
      adj.list = case$adj.list,
      weight.list = weight.list,
      transport.rule = "regression.gradient",
      coordinates = coords,
      polynomial.probes = probes,
      gradient.coordinate.method = "coordinates",
      gradient.ridge = ridge
    ),
    silent = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(fit, "try-error")) {
    return(list(
      summary = transported.hessian.quadric.error.summary(case, policy,
                                                          graph.weight,
                                                          elapsed,
                                                          as.character(fit)),
      gradient = data.frame()
    ))
  }
  list(
    summary = transported.hessian.quadric.fit.summary(fit, case, policy,
                                                      graph.weight, elapsed),
    gradient = transported.hessian.quadric.gradient.table(fit, case, policy,
                                                          graph.weight)
  )
}

transported.hessian.quadric.fit.summary <- function(op,
                                                    case,
                                                    policy,
                                                    graph.weight,
                                                    elapsed.sec) {
  grad <- op$transport$gradient.diagnostics %||% data.frame()
  res <- transported.hessian.quadric.residual.summary(op)
  data.frame(
    case.id = case$case.id,
    intrinsic.dim = case$intrinsic.dim,
    surface = case$surface,
    curvature = case$curvature,
    n = case$n,
    seed = case$seed,
    k = case$k,
    index.k = case$index.k,
    graph.weight = graph.weight,
    coordinate.source = policy$coordinate.source,
    coordinate.label = policy$label,
    coordinate.dim = ncol(transported.hessian.quadric.coordinates(case, policy$coordinate.source)),
    status = "ok",
    runtime.sec = elapsed.sec,
    n.edges = case$graph.metadata$n.edges,
    n.mst.edges.added = case$graph.metadata$n.mst.edges.added,
    n.rows = op$diagnostics$summary$n.rows,
    nullity.estimate = op$diagnostics$nullity.estimate,
    median.local.rank = if (nrow(grad)) stats::median(grad$rank, na.rm = TRUE) else NA_real_,
    min.local.rank = if (nrow(grad)) min(grad$rank, na.rm = TRUE) else NA_real_,
    rank.deficient.fraction = if (nrow(grad)) mean(grad$rank < ncol(transported.hessian.quadric.coordinates(case, policy$coordinate.source)), na.rm = TRUE) else NA_real_,
    median.local.condition = if (nrow(grad)) {
      finite <- grad$condition[is.finite(grad$condition)]
      if (length(finite)) stats::median(finite) else NA_real_
    } else {
      NA_real_
    },
    message = NA_character_,
    res,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.error.summary <- function(case,
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
    index.k = case$index.k,
    graph.weight = graph.weight,
    coordinate.source = policy$coordinate.source,
    coordinate.label = policy$label,
    coordinate.dim = NA_integer_,
    status = "error",
    runtime.sec = elapsed.sec,
    n.edges = case$graph.metadata$n.edges,
    n.mst.edges.added = case$graph.metadata$n.mst.edges.added,
    n.rows = NA_integer_,
    nullity.estimate = NA_integer_,
    median.local.rank = NA_real_,
    min.local.rank = NA_real_,
    rank.deficient.fraction = NA_real_,
    median.local.condition = NA_real_,
    message = message,
    residual.constant = NA_real_,
    residual.linear.mean = NA_real_,
    residual.linear.max = NA_real_,
    residual.quadratic.mean = NA_real_,
    residual.quadratic.min = NA_real_,
    residual.overall = NA_real_,
    stringsAsFactors = FALSE
  )
}

transported.hessian.quadric.gradient.table <- function(op,
                                                       case,
                                                       policy,
                                                       graph.weight) {
  grad <- op$transport$gradient.diagnostics %||% data.frame()
  if (!nrow(grad)) return(data.frame())
  grad$case.id <- case$case.id
  grad$intrinsic.dim <- case$intrinsic.dim
  grad$surface <- case$surface
  grad$curvature <- case$curvature
  grad$graph.weight <- graph.weight
  grad$coordinate.source <- policy$coordinate.source
  grad
}

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

run.transported.hessian.quadric.coordinate.validation <- function(mode = c("smoke", "full"),
                                                                  repo.dir = repo.dir.default,
                                                                  output.dir = file.path(repo.dir, "dev/graph-trend-filtering/reports"),
                                                                  graph.weights = c("ambient", "surface.segment")) {
  mode <- match.arg(mode)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  transported.hessian.quadric.load.packages(repo.dir)
  cases <- transported.hessian.quadric.validation.cases(mode)
  policies <- transported.hessian.quadric.coordinate.policies()
  runs <- list()
  counter <- 0L
  for (case in cases) {
    for (graph.weight in graph.weights) {
      for (policy in policies) {
        counter <- counter + 1L
        runs[[counter]] <- transported.hessian.quadric.run.one(
          case = case,
          policy = policy,
          graph.weight = graph.weight
        )
      }
    }
  }
  summary <- .transported.hessian.bind.rows(lapply(runs, `[[`, "summary"))
  gradient <- .transported.hessian.bind.rows(lapply(runs, `[[`, "gradient"))
  metadata <- list(
    mode = mode,
    generated.at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    gflow.repo = repo.dir,
    n.cases = length(cases),
    n.runs = nrow(summary),
    graph.weights = graph.weights
  )
  out <- list(
    metadata = metadata,
    summary = summary,
    gradient = gradient
  )
  prefix <- file.path(output.dir, paste0("transported_hessian_quadric_coordinate_validation_", mode))
  saveRDS(out, paste0(prefix, ".rds"))
  utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(gradient, paste0(prefix, "_gradient.csv"), row.names = FALSE)
  out
}
