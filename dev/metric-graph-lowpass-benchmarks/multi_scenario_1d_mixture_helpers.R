normalize_to_unit_max <- function(y) {
  if (max(y) <= 0) y else y / max(y)
}

mixture_signal <- function(x, means, sds, weights) {
  y <- numeric(length(x))
  for (j in seq_along(means)) {
    y <- y + weights[j] * stats::dnorm(x, mean = means[j], sd = sds[j])
  }
  normalize_to_unit_max(y)
}

shape_catalog <- function() {
  data.frame(
    shape.id = sprintf("S%02d", 1:16),
    shape.name = c(
      "two separated narrow equal", "two separated unequal",
      "two close overlapping equal", "two close overlapping unequal",
      "narrow left plus broad right", "broad left plus narrow right",
      "left boundary plus interior", "interior plus right boundary",
      "deep valley separated", "nearly unimodal overlap",
      "three separated peaks", "dominant center side bumps",
      "two close plus distant", "narrow center spike",
      "left-heavy triple", "boundary interior boundary"
    ),
    means = I(list(
      c(0.28, 0.72), c(0.25, 0.74), c(0.43, 0.56), c(0.40, 0.54),
      c(0.30, 0.69), c(0.31, 0.70), c(0.08, 0.48), c(0.52, 0.93),
      c(0.22, 0.78), c(0.45, 0.56), c(0.20, 0.50, 0.80),
      c(0.23, 0.52, 0.79), c(0.35, 0.46, 0.82),
      c(0.24, 0.50, 0.76), c(0.16, 0.38, 0.70),
      c(0.08, 0.50, 0.92)
    )),
    sds = I(list(
      c(0.045, 0.045), c(0.050, 0.070), c(0.080, 0.080),
      c(0.060, 0.095), c(0.035, 0.140), c(0.140, 0.035),
      c(0.035, 0.090), c(0.090, 0.035), c(0.040, 0.040),
      c(0.115, 0.115), c(0.040, 0.050, 0.040),
      c(0.060, 0.045, 0.060), c(0.045, 0.045, 0.060),
      c(0.090, 0.025, 0.090), c(0.040, 0.060, 0.080),
      c(0.035, 0.075, 0.035)
    )),
    weights = I(list(
      c(1.00, 1.00), c(1.00, 0.55), c(1.00, 1.00), c(1.00, 0.45),
      c(1.00, 0.85), c(0.85, 1.00), c(0.85, 1.00), c(1.00, 0.85),
      c(1.00, 1.00), c(1.00, 0.90), c(0.85, 1.00, 0.75),
      c(0.45, 1.00, 0.45), c(0.90, 0.75, 1.00),
      c(0.55, 1.00, 0.55), c(1.00, 0.75, 0.45),
      c(0.70, 1.00, 0.70)
    )),
    stringsAsFactors = FALSE
  )
}

variant_catalog <- function() {
  data.frame(
    variant.id = c("V1", "V2", "V3"),
    variant.name = c("uniform gaussian", "center-biased heteroskedastic",
                     "gap-biased laplace outliers"),
    stringsAsFactors = FALSE
  )
}

sample_x <- function(n, variant.id) {
  if (identical(variant.id, "V1")) return(sort(stats::runif(n)))
  if (identical(variant.id, "V2")) {
    p0 <- stats::pnorm(0, mean = 0.50, sd = 0.20)
    p1 <- stats::pnorm(1, mean = 0.50, sd = 0.20)
    return(sort(stats::qnorm(stats::runif(n, p0, p1), mean = 0.50, sd = 0.20)))
  }
  n.left <- floor(0.42 * n)
  sort(c(stats::runif(n.left, 0.02, 0.42), stats::runif(n - n.left, 0.58, 0.98)))
}

rlaplace <- function(n, scale) {
  u <- stats::runif(n) - 0.5
  -scale * sign(u) * log1p(-2 * abs(u))
}

add_noise <- function(x, truth, variant.id) {
  n <- length(x)
  if (identical(variant.id, "V1")) return(truth + stats::rnorm(n, sd = 0.08))
  if (identical(variant.id, "V2")) return(truth + stats::rnorm(n, sd = 0.045 + 0.14 * truth))
  y <- truth + rlaplace(n, scale = 0.055)
  idx <- sample.int(n, max(1L, round(0.04 * n)))
  y[idx] <- y[idx] + stats::rnorm(length(idx), mean = 0, sd = 0.45)
  y
}

make_dataset <- function(shape, variant, replicate, n = 250L, base.seed = 73001L) {
  seed <- base.seed + as.integer(sub("^S", "", shape$shape.id)) * 1000L +
    as.integer(sub("^V", "", variant$variant.id)) * 100L + replicate
  set.seed(seed)
  x <- sample_x(n, variant$variant.id)
  truth <- mixture_signal(x, shape$means[[1]], shape$sds[[1]], shape$weights[[1]])
  list(seed = seed, x = x, X = matrix(x, ncol = 1), truth = truth,
       y = add_noise(x, truth, variant$variant.id))
}

safe_median <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

conductance_label <- function(rule, alpha, sigma.rule, sigma.quantile, local.k) {
  if (identical(rule, "inverse.length.power")) return(sprintf("%s alpha=%s", rule, signif(alpha, 3)))
  if (identical(rule, "exp.length") || identical(rule, "exp.length.squared")) {
    if (identical(sigma.rule, "edge.quantile")) return(sprintf("%s sigma=edge.q%s", rule, signif(sigma.quantile, 3)))
    return(sprintf("%s sigma=%s", rule, sigma.rule))
  }
  if (identical(rule, "self.tuned.gaussian")) return(sprintf("%s local.k=%d", rule, as.integer(local.k)))
  rule
}

metric_grid <- function() {
  out <- rbind(
    data.frame(conductance.rule = "inverse.length.power",
               alpha = c(0.5, 1, 1.5, 2, 3), sigma.rule = NA_character_,
               sigma.quantile = NA_real_, local.k = NA_integer_),
    data.frame(conductance.rule = c("exp.length", "exp.length", "exp.length.squared", "exp.length.squared"),
               alpha = NA_real_, sigma.rule = c("median", "edge.quantile", "median", "edge.quantile"),
               sigma.quantile = c(NA_real_, 0.75, NA_real_, 0.75), local.k = NA_integer_),
    data.frame(conductance.rule = "self.tuned.gaussian", alpha = NA_real_,
               sigma.rule = NA_character_, sigma.quantile = NA_real_,
               local.k = c(3L, 5L, 10L))
  )
  out$laplacian.type <- "unnormalized"
  out$filter.type <- "heat_kernel"
  out$n.eigenpairs <- 60L
  out$n.candidates <- 30L
  out$config.id <- sprintf("M%02d", seq_len(nrow(out)))
  out$config.label <- mapply(conductance_label, out$conductance.rule, out$alpha,
                             out$sigma.rule, out$sigma.quantile, out$local.k,
                             USE.NAMES = FALSE)
  out
}

make_scenario_grid <- function(replicates = 1:5,
                               shape.ids = NULL,
                               variant.ids = NULL) {
  shapes <- shape_catalog()
  variants <- variant_catalog()
  if (!is.null(shape.ids)) {
    shapes <- shapes[shapes$shape.id %in% shape.ids, , drop = FALSE]
  }
  if (!is.null(variant.ids)) {
    variants <- variants[variants$variant.id %in% variant.ids, , drop = FALSE]
  }
  scenario.grid <- merge(
    merge(shapes[, c("shape.id", "shape.name")], variants, by = NULL),
    data.frame(replicate = replicates),
    by = NULL
  )
  scenario.grid$dataset.id <- sprintf("%s_%s_R%02d", scenario.grid$shape.id,
                                      scenario.grid$variant.id,
                                      scenario.grid$replicate)
  scenario.grid
}

make_path_graph <- function(x) {
  n <- length(x)
  if (n < 2L) stop("Path graph requires at least two vertices")
  ord <- order(x)
  adj.list <- vector("list", n)
  weight.list <- vector("list", n)
  for (ii in seq_len(n - 1L)) {
    u <- ord[ii]
    v <- ord[ii + 1L]
    len <- abs(x[v] - x[u])
    adj.list[[u]] <- c(adj.list[[u]], v)
    weight.list[[u]] <- c(weight.list[[u]], len)
    adj.list[[v]] <- c(adj.list[[v]], u)
    weight.list[[v]] <- c(weight.list[[v]], len)
  }
  list(
    adj_list = adj.list,
    weight_list = weight.list,
    graph_type = "path",
    n_components_before = 1L,
    n_components_after = 1L,
    n_mst_edges_added = 0L
  )
}

make_benchmark_graph <- function(X, x, graph.type = c("path", "sknn.mst"),
                                 graph.k = 5L) {
  graph.type <- match.arg(graph.type)
  if (identical(graph.type, "path")) {
    return(make_path_graph(x))
  }
  create.sknn.graph(X, k = graph.k, connect.components = TRUE,
                    connect.method = "component.mst")
}

graph_edge_count <- function(graph) {
  sum(lengths(graph$adj_list)) / 2
}

conductance_adjacency <- function(graph,
                                  conductance.rule = "exp.length",
                                  conductance.sigma.rule = "edge.quantile",
                                  conductance.sigma.quantile = 0.75) {
  op <- metric.graph.lowpass.operator(
    adj.list = graph$adj_list,
    weight.list = graph$weight_list,
    conductance.rule = conductance.rule,
    conductance.sigma.rule = conductance.sigma.rule,
    conductance.sigma.quantile = conductance.sigma.quantile,
    laplacian.type = "unnormalized",
    return.sparse = FALSE,
    verbose = FALSE
  )
  Matrix::sparseMatrix(
    i = c(op$edge.table$u, op$edge.table$v),
    j = c(op$edge.table$v, op$edge.table$u),
    x = rep(op$edge.table$conductance, 2L),
    dims = c(length(graph$adj_list), length(graph$adj_list))
  )
}

score_prediction <- function(method, yhat, elapsed, y, truth,
                             eta = NA_real_, gcv = NA_real_,
                             effective.df = NA_real_,
                             eigen.backend = NA_character_,
                             conductance.median = NA_real_,
                             message = "") {
  yhat <- unname(as.numeric(yhat))
  if (length(yhat) != length(y) || any(!is.finite(yhat))) {
    return(score_nonfinite(method, elapsed))
  }
  data.frame(
    method = method,
    status = "ok",
    rmse.truth = sqrt(mean((yhat - truth)^2)),
    rmse.observed = sqrt(mean((yhat - y)^2)),
    mae.truth = mean(abs(yhat - truth)),
    runtime.sec = unname(elapsed),
    eta = unname(eta),
    gcv = unname(gcv),
    effective.df = unname(effective.df),
    eigen.backend = unname(eigen.backend),
    conductance.median = unname(conductance.median),
    message = message,
    stringsAsFactors = FALSE
  )
}

score_nonfinite <- function(method, elapsed = NA_real_) {
  data.frame(
    method = method,
    status = "score_error",
    rmse.truth = NA_real_,
    rmse.observed = NA_real_,
    mae.truth = NA_real_,
    runtime.sec = elapsed,
    eta = NA_real_,
    gcv = NA_real_,
    effective.df = NA_real_,
    eigen.backend = NA_character_,
    conductance.median = NA_real_,
    message = "Fit returned, but fitted values were missing or non-finite.",
    stringsAsFactors = FALSE
  )
}

score_failure <- function(method, err, elapsed = NA_real_) {
  data.frame(
    method = method,
    status = "fit_error",
    rmse.truth = NA_real_,
    rmse.observed = NA_real_,
    mae.truth = NA_real_,
    runtime.sec = elapsed,
    eta = NA_real_,
    gcv = NA_real_,
    effective.df = NA_real_,
    eigen.backend = NA_character_,
    conductance.median = NA_real_,
    message = if (inherits(err, "condition")) conditionMessage(err) else as.character(err),
    stringsAsFactors = FALSE
  )
}

time_expression <- function(expr) {
  start <- proc.time()
  out <- tryCatch(force(expr), error = function(e) e)
  list(result = out, elapsed = unname((proc.time() - start)[["elapsed"]]))
}

fit_smooth_spline_gcv <- function(x, y) {
  time_expression({
    fit <- stats::smooth.spline(x = x, y = y, cv = FALSE)
    list(fitted.values = stats::predict(fit, x = x)$y,
         gcv = fit$cv.crit,
         effective.df = fit$df)
  })
}

fit_mgcv_gam_reml <- function(x, y) {
  time_expression({
    if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is not installed")
    dat <- data.frame(x = x, y = y)
    fit <- mgcv::gam(y ~ s(x, bs = "cs", k = min(40L, length(x) - 1L)),
                     data = dat, method = "REML")
    list(fitted.values = stats::predict(fit, newdata = dat),
         gcv = fit$gcv.ubre,
         effective.df = sum(fit$edf))
  })
}

fit_loess_gcv <- function(x, y,
                          spans = c(0.15, 0.25, 0.35, 0.50, 0.70)) {
  time_expression({
    dat <- data.frame(x = x, y = y)
    fits <- lapply(spans, function(span) {
      fit <- stats::loess(
        y ~ x,
        data = dat,
        span = span,
        degree = 2L,
        family = "gaussian",
        control = stats::loess.control(surface = "direct", trace.hat = "exact")
      )
      pred <- stats::predict(fit, newdata = dat)
      tr <- if (!is.null(fit$trace.hat)) fit$trace.hat else NA_real_
      gcv <- if (is.finite(tr) && tr < length(y)) {
        mean((y - pred)^2, na.rm = TRUE) / (1 - tr / length(y))^2
      } else {
        Inf
      }
      list(fit = fit, pred = pred, span = span, gcv = gcv, effective.df = tr)
    })
    gcvs <- vapply(fits, `[[`, numeric(1), "gcv")
    best <- fits[[which.min(gcvs)]]
    list(fitted.values = best$pred,
         eta = best$span,
         gcv = best$gcv,
         effective.df = best$effective.df)
  })
}

fit_np_npreg_cv <- function(x, y) {
  time_expression({
    if (!requireNamespace("np", quietly = TRUE)) stop("Package 'np' is not installed")
    dat <- data.frame(x = x)
    bw <- np::npregbw(xdat = dat, ydat = y, regtype = "ll", bwmethod = "cv.ls")
    fit <- np::npreg(bws = bw)
    list(fitted.values = as.numeric(stats::predict(fit, exdat = dat)),
         eta = as.numeric(bw$bw),
         gcv = bw$fval)
  })
}

fit_locfit_fixed <- function(x, y) {
  time_expression({
    if (!requireNamespace("locfit", quietly = TRUE)) stop("Package 'locfit' is not installed")
    suppressPackageStartupMessages(require("locfit", quietly = TRUE, character.only = TRUE))
    fit <- locfit::locfit.raw(x = matrix(x, ncol = 1), y = y,
                              alpha = 0.35, deg = 2L,
                              family = "gaussian")
    list(fitted.values = as.numeric(stats::fitted(fit)),
         eta = 0.35)
  })
}

fit_genlasso_trendfilter_cv <- function(x, y) {
  time_expression({
    if (!requireNamespace("genlasso", quietly = TRUE)) stop("Package 'genlasso' is not installed")
    fit <- genlasso::trendfilter(y, pos = x, ord = 2L)
    cv <- genlasso::cv.trendfilter(fit, k = 5L, mode = "lambda", verbose = FALSE)
    cf <- stats::coef(fit, lambda = cv$lambda.min)
    list(fitted.values = as.numeric(cf$beta),
         eta = cv$lambda.min,
         gcv = min(cv$err, na.rm = TRUE),
         effective.df = as.numeric(cf$df))
  })
}

fit_gsd_gsmoothing <- function(graph, y) {
  time_expression({
    if (!requireNamespace("GSD", quietly = TRUE)) stop("Package 'GSD' is not installed")
    A <- as.matrix(conductance_adjacency(graph))
    list(fitted.values = as.numeric(GSD::gsmoothing(A, y)))
  })
}

fit_gasper_sgwt_sure <- function(graph, y) {
  time_expression({
    if (!requireNamespace("gasper", quietly = TRUE)) stop("Package 'gasper' is not installed")
    A <- as.matrix(conductance_adjacency(graph))
    L <- gasper::laplacian_mat(A, type = "unnormalized")
    eig <- gasper::eigensort(L)
    tf <- gasper::tight_frame(eig$evalues, eig$evectors)
    wcn <- gasper::analysis(y, tf)
    diagWWt <- colSums(t(tf)^2)
    sigma2 <- gasper::GVN(y, A, L)
    sigma <- sqrt(max(sigma2, .Machine$double.eps))
    thresh <- sort(unique(abs(wcn)))
    opt <- gasper::SUREthresh(wcn, thresh, diagWWt, beta = 2,
                              sigma = sigma, keepwc = TRUE)
    wc <- opt$wc[, opt$min[["xminSURE"]]]
    list(fitted.values = as.numeric(gasper::synthesis(wc, tf)),
         eta = opt$thr[opt$min[["xminSURE"]]],
         gcv = min(opt$res$SURE, na.rm = TRUE))
  })
}

score_prediction_result <- function(method, fit.out, y, truth) {
  if (inherits(fit.out$result, "error")) {
    return(score_failure(method, fit.out$result, fit.out$elapsed))
  }
  score_prediction(
    method = method,
    yhat = fit.out$result$fitted.values,
    elapsed = fit.out$elapsed,
    y = y,
    truth = truth,
    eta = if (!is.null(fit.out$result$eta)) fit.out$result$eta else NA_real_,
    gcv = if (!is.null(fit.out$result$gcv)) fit.out$result$gcv else NA_real_,
    effective.df = if (!is.null(fit.out$result$effective.df)) fit.out$result$effective.df else NA_real_
  )
}

fit_baseline_methods <- function(x, y, graph) {
  list(
    smooth.spline.gcv = fit_smooth_spline_gcv(x, y),
    mgcv.gam.reml = fit_mgcv_gam_reml(x, y),
    loess.gcv = fit_loess_gcv(x, y),
    np.npreg.cv = fit_np_npreg_cv(x, y),
    genlasso.trendfilter.cv = fit_genlasso_trendfilter_cv(x, y),
    GSD.gsmoothing = fit_gsd_gsmoothing(graph, y),
    gasper.sgwt.sure = fit_gasper_sgwt_sure(graph, y)
  )
}

method_family <- function(method) {
  if (identical(method, "fit.metric.graph.lowpass")) return("metric low-pass")
  if (identical(method, "fit.rdgraph.regression")) return("rdgraph")
  if (method %in% c("GSD.gsmoothing", "gasper.sgwt.sure")) return("graph denoising")
  "classical 1D"
}

result_row <- function(scenario.row, graph, graph.type, graph.k,
                       graph.elapsed, config.id, config.label, score) {
  cbind(
    scenario.row,
    data.frame(
      graph.type = graph.type,
      graph.k = if (identical(graph.type, "path")) NA_integer_ else graph.k,
      n.components.before = graph$n_components_before,
      n.components.after = graph$n_components_after,
      n.mst.edges.added = graph$n_mst_edges_added,
      n.edges = graph_edge_count(graph),
      graph.runtime.sec = graph.elapsed,
      config.id = config.id,
      config.label = config.label,
      method.family = method_family(score$method[1]),
      score,
      stringsAsFactors = FALSE
    )
  )
}

run_multi_method_benchmark <- function(repo.dir = "/Users/pgajer/current_projects/gflow",
                                       n = 250L,
                                       graph.type = c("path", "sknn.mst"),
                                       graph.k = 5L,
                                       replicates = 1:5,
                                       shape.ids = NULL,
                                       variant.ids = NULL,
                                       include.baselines = TRUE,
                                       progress = TRUE) {
  graph.type <- match.arg(graph.type)
  shapes <- shape_catalog()
  variants <- variant_catalog()
  scenario.grid <- make_scenario_grid(
    replicates = replicates,
    shape.ids = shape.ids,
    variant.ids = variant.ids
  )
  metric.configs <- metric_grid()
  n.methods <- nrow(metric.configs) + 1L + if (isTRUE(include.baselines)) 7L else 0L
  results <- vector("list", nrow(scenario.grid) * n.methods)
  idx <- 0L

  for (s in seq_len(nrow(scenario.grid))) {
    if (isTRUE(progress)) {
      message(sprintf("[%d/%d] %s", s, nrow(scenario.grid), scenario.grid$dataset.id[s]))
    }
    shape <- shapes[match(scenario.grid$shape.id[s], shapes$shape.id), ]
    variant <- variants[match(scenario.grid$variant.id[s], variants$variant.id), ]
    dat <- make_dataset(shape, variant, scenario.grid$replicate[s], n = n)

    graph.out <- time_expression(
      make_benchmark_graph(dat$X, dat$x, graph.type = graph.type, graph.k = graph.k)
    )
    if (inherits(graph.out$result, "error")) {
      fake.graph <- list(
        adj_list = vector("list", n),
        weight_list = vector("list", n),
        n_components_before = NA_integer_,
        n_components_after = NA_integer_,
        n_mst_edges_added = NA_integer_
      )
      idx <- idx + 1L
      results[[idx]] <- result_row(
        scenario.grid[s, ], fake.graph, graph.type, graph.k, graph.out$elapsed,
        "graph", "graph construction",
        score_failure("graph construction", graph.out$result, graph.out$elapsed)
      )
      next
    }
    graph <- graph.out$result

    rd.fit <- time_expression(
      fit.rdgraph.regression(
        X = dat$X,
        y = dat$y,
        k = if (identical(graph.type, "path")) 2L else graph.k,
        adj.list = graph$adj_list,
        weight.list = graph$weight_list,
        verbose.level = 0L
      )
    )
    idx <- idx + 1L
    results[[idx]] <- result_row(
      scenario.grid[s, ], graph, graph.type, graph.k, graph.out$elapsed,
      "rdgraph", paste0("fit.rdgraph.regression.", graph.type),
      score_fit("fit.rdgraph.regression", rd.fit$result, rd.fit$elapsed,
                dat$y, dat$truth)
    )

    for (m in seq_len(nrow(metric.configs))) {
      spec <- metric.configs[m, ]
      metric.out <- fit_metric_one(spec, graph, dat$y, n = n)
      idx <- idx + 1L
      results[[idx]] <- result_row(
        scenario.grid[s, ], graph, graph.type, graph.k, graph.out$elapsed,
        spec$config.id, spec$config.label,
        score_fit("fit.metric.graph.lowpass", metric.out$result,
                  metric.out$elapsed, dat$y, dat$truth)
      )
    }

    if (isTRUE(include.baselines)) {
      baseline.fits <- fit_baseline_methods(dat$x, dat$y, graph)
      for (nm in names(baseline.fits)) {
        idx <- idx + 1L
        results[[idx]] <- result_row(
          scenario.grid[s, ], graph, graph.type, graph.k, graph.out$elapsed,
          nm, nm,
          score_prediction_result(nm, baseline.fits[[nm]], dat$y, dat$truth)
        )
      }
    }
  }

  out <- do.call(rbind, results[seq_len(idx)])
  rownames(out) <- NULL
  out
}

fit_metric_one <- function(spec, graph, y, n) {
  args <- list(
    adj.list = graph$adj_list, weight.list = graph$weight_list, y = y,
    conductance.rule = spec$conductance.rule, laplacian.type = "unnormalized",
    filter.type = spec$filter.type, n.eigenpairs = min(spec$n.eigenpairs, n - 1L),
    n.candidates = spec$n.candidates, eigen.solver = "auto",
    dense.eigen.threshold = 200L, dense.fallback = "auto",
    dense.fallback.threshold = 5000L, verbose = FALSE
  )
  if (identical(spec$conductance.rule, "inverse.length.power")) args$conductance.alpha <- spec$alpha
  if (identical(spec$conductance.rule, "exp.length") || identical(spec$conductance.rule, "exp.length.squared")) {
    args$conductance.sigma.rule <- spec$sigma.rule
    if (!is.na(spec$sigma.quantile)) args$conductance.sigma.quantile <- spec$sigma.quantile
  }
  if (identical(spec$conductance.rule, "self.tuned.gaussian")) args$conductance.local.k <- spec$local.k
  start <- proc.time()
  out <- tryCatch(do.call(fit.metric.graph.lowpass, args), error = function(e) e)
  list(result = out, elapsed = unname((proc.time() - start)[["elapsed"]]))
}

score_fit <- function(method, fit, elapsed, y, truth) {
  if (inherits(fit, "error")) {
    return(data.frame(method = method, status = "fit_error", rmse.truth = NA_real_,
                      rmse.observed = NA_real_, mae.truth = NA_real_,
                      runtime.sec = elapsed, eta = NA_real_, gcv = NA_real_,
                      effective.df = NA_real_, eigen.backend = NA_character_,
                      conductance.median = NA_real_, message = conditionMessage(fit),
                      stringsAsFactors = FALSE))
  }
  yhat <- as.numeric(fit$fitted.values)
  if (length(yhat) != length(y) || any(!is.finite(yhat))) {
    return(score_nonfinite(method, elapsed))
  }
  if (inherits(fit, "metric.graph.lowpass.fit")) {
    conductance <- fit$operator$edge.table$conductance
    eta <- fit$gcv$eta.optimal
    gcv <- fit$gcv$gcv.optimal
    edf <- fit$gcv$effective.df
    backend <- fit$spectral$backend
  } else {
    conductance <- NA_real_
    eta <- tail(fit$gcv$eta.optimal, 1)
    gcv <- tail(fit$gcv$gcv.optimal, 1)
    edf <- NA_real_
    backend <- if (!is.null(fit$spectral$backend)) fit$spectral$backend else NA_character_
  }
  data.frame(method = method, status = "ok",
             rmse.truth = sqrt(mean((yhat - truth)^2)),
             rmse.observed = sqrt(mean((yhat - y)^2)),
             mae.truth = mean(abs(yhat - truth)), runtime.sec = elapsed,
             eta = eta, gcv = gcv, effective.df = edf,
             eigen.backend = backend, conductance.median = safe_median(conductance),
             message = "", stringsAsFactors = FALSE)
}

run_multi_scenario_benchmark <- function(repo.dir = "/Users/pgajer/current_projects/gflow",
                                         n = 250L,
                                         graph.k = 5L,
                                         replicates = 1:5,
                                         progress = TRUE) {
  shapes <- shape_catalog()
  variants <- variant_catalog()
  scenario.grid <- make_scenario_grid(replicates)
  metric.configs <- metric_grid()
  results <- vector("list", nrow(scenario.grid) * (nrow(metric.configs) + 1L))
  idx <- 0L

  for (s in seq_len(nrow(scenario.grid))) {
    if (isTRUE(progress)) {
      message(sprintf("[%d/%d] %s", s, nrow(scenario.grid), scenario.grid$dataset.id[s]))
    }
    shape <- shapes[match(scenario.grid$shape.id[s], shapes$shape.id), ]
    variant <- variants[match(scenario.grid$variant.id[s], variants$variant.id), ]
    dat <- make_dataset(shape, variant, scenario.grid$replicate[s], n = n)

    graph.start <- proc.time()
    graph <- tryCatch(
      create.sknn.graph(dat$X, k = graph.k, connect.components = TRUE,
                        connect.method = "component.mst"),
      error = function(e) e
    )
    graph.elapsed <- unname((proc.time() - graph.start)[["elapsed"]])

    if (inherits(graph, "error")) {
      idx <- idx + 1L
      results[[idx]] <- cbind(
        scenario.grid[s, ],
        data.frame(graph.k = graph.k, n.components.before = NA_integer_,
                   n.components.after = NA_integer_, n.mst.edges.added = NA_integer_,
                   graph.runtime.sec = graph.elapsed, config.id = "rdgraph",
                   config.label = "fit.rdgraph.regression",
                   score_fit("fit.rdgraph.regression", graph, NA_real_, dat$y, dat$truth),
                   stringsAsFactors = FALSE)
      )
      next
    }

    rd.start <- proc.time()
    rd.fit <- tryCatch(
      fit.rdgraph.regression(X = dat$X, y = dat$y, k = graph.k,
                             adj.list = graph$adj_list, weight.list = graph$weight_list,
                             verbose.level = 0L),
      error = function(e) e
    )
    idx <- idx + 1L
    results[[idx]] <- cbind(
      scenario.grid[s, ],
      data.frame(graph.k = graph.k, n.components.before = graph$n_components_before,
                 n.components.after = graph$n_components_after,
                 n.mst.edges.added = graph$n_mst_edges_added,
                 graph.runtime.sec = graph.elapsed, config.id = "rdgraph",
                 config.label = "fit.rdgraph.regression",
                 score_fit("fit.rdgraph.regression", rd.fit,
                           unname((proc.time() - rd.start)[["elapsed"]]), dat$y, dat$truth),
                 stringsAsFactors = FALSE)
    )

    for (m in seq_len(nrow(metric.configs))) {
      spec <- metric.configs[m, ]
      metric.out <- fit_metric_one(spec, graph, dat$y, n = n)
      idx <- idx + 1L
      results[[idx]] <- cbind(
        scenario.grid[s, ],
        data.frame(graph.k = graph.k, n.components.before = graph$n_components_before,
                   n.components.after = graph$n_components_after,
                   n.mst.edges.added = graph$n_mst_edges_added,
                   graph.runtime.sec = graph.elapsed, config.id = spec$config.id,
                   config.label = spec$config.label,
                   score_fit("fit.metric.graph.lowpass", metric.out$result,
                             metric.out$elapsed, dat$y, dat$truth),
                   stringsAsFactors = FALSE)
      )
    }
  }
  out <- do.call(rbind, results[seq_len(idx)])
  rownames(out) <- NULL
  out
}
