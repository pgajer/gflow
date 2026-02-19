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
  base <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase2")
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

    log.msg(sprintf("Read attempt %d failed for %s | %s",
                    attempt, basename(src), conditionMessage(dt)))
    Sys.sleep(sleep.sec)
  }
  stop("Could not snapshot/read file after retries: ", src, call. = FALSE)
}

save.plot <- function(plot.obj, path, width = 12, height = 6.5, dpi = 160) {
  ggsave(filename = path, plot = plot.obj, width = width, height = height, dpi = dpi)
}

format.run.id <- function(run.dir) basename(normalizePath(run.dir, winslash = "/", mustWork = FALSE))

weighted.mean.safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted.se.safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  sw <- sum(w)
  if (sw <= 0) return(NA_real_)
  mu <- sum(w * x) / sw
  var.w <- sum(w * (x - mu)^2) / sw
  ess <- (sw^2) / sum(w^2)
  if (!is.finite(ess) || ess <= 1) return(0.0)
  sqrt(var.w / ess)
}

safe.col.scales <- function(X.train) {
  med <- apply(X.train, 2L, stats::mad, constant = 1.4826)
  sdv <- apply(X.train, 2L, stats::sd)
  out <- med
  bad <- !is.finite(out) | out <= sqrt(.Machine$double.eps)
  out[bad] <- sdv[bad]
  bad <- !is.finite(out) | out <= sqrt(.Machine$double.eps)
  out[bad] <- 1.0
  as.double(out)
}

rlaplace.local <- function(n, location = 0, scale = 1) {
  u <- stats::runif(n, min = -0.5, max = 0.5)
  location - scale * sign(u) * log1p(-2 * abs(u))
}

draw.radial.noise <- function(n, noise, noise.type) {
  if (!is.finite(noise) || noise <= 0) return(rep(0.0, n))
  if (identical(noise.type, "laplace")) {
    rlaplace.local(n, location = 0, scale = noise / sqrt(2))
  } else {
    stats::rnorm(n, mean = 0, sd = noise)
  }
}

generate.circle.random <- function(n, radius, noise, noise.type = "normal", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  angles <- sort(stats::runif(n, min = 0, max = 2 * pi))
  eps <- draw.radial.noise(n, noise = noise, noise.type = noise.type)
  data.frame(
    x = (radius + eps) * cos(angles),
    y = (radius + eps) * sin(angles),
    x.noise.free = radius * cos(angles),
    y.noise.free = radius * sin(angles),
    angles = angles,
    is_outlier = FALSE,
    stringsAsFactors = FALSE
  )
}

generate.circle.nonuniform <- function(n, radius, noise, noise.type = "normal", seed = NULL,
                                       arc1.weight = 0.45, arc2.weight = 0.35,
                                       arc.halfwidth = pi / 6) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::runif(n)
  angles <- numeric(n)

  idx.arc1 <- u <= arc1.weight
  idx.arc2 <- u > arc1.weight & u <= (arc1.weight + arc2.weight)
  idx.bg <- !(idx.arc1 | idx.arc2)

  angles[idx.arc1] <- stats::runif(sum(idx.arc1), -arc.halfwidth, arc.halfwidth)
  angles[idx.arc2] <- pi + stats::runif(sum(idx.arc2), -arc.halfwidth, arc.halfwidth)
  angles[idx.bg] <- stats::runif(sum(idx.bg), 0, 2 * pi)
  angles <- (angles %% (2 * pi))

  eps <- draw.radial.noise(n, noise = noise, noise.type = noise.type)
  radial <- pmax(radius + eps, 1e-6)

  data.frame(
    x = radial * cos(angles),
    y = radial * sin(angles),
    x.noise.free = radius * cos(angles),
    y.noise.free = radius * sin(angles),
    angles = angles,
    is_outlier = FALSE,
    stringsAsFactors = FALSE
  )
}

generate.circle.outliers <- function(n, radius, noise, noise.type = "normal", seed = NULL,
                                     outlier.frac = 0.05, outlier.radial.shift = 0.8,
                                     outlier.radial.sd = 0.25) {
  if (!is.null(seed)) set.seed(seed)
  angles <- stats::runif(n, 0, 2 * pi)
  eps <- draw.radial.noise(n, noise = noise, noise.type = noise.type)
  radial <- pmax(radius + eps, 1e-6)

  n.out <- max(1L, as.integer(round(outlier.frac * n)))
  out.idx <- sort(sample.int(n, size = n.out, replace = FALSE))
  bump <- outlier.radial.shift + abs(stats::rnorm(n.out, mean = 0, sd = outlier.radial.sd))
  signs <- sample(c(-1, 1), size = n.out, replace = TRUE)
  radial[out.idx] <- pmax(radial[out.idx] + signs * bump, 1e-6)

  is.outlier <- rep(FALSE, n)
  is.outlier[out.idx] <- TRUE

  data.frame(
    x = radial * cos(angles),
    y = radial * sin(angles),
    x.noise.free = radius * cos(angles),
    y.noise.free = radius * sin(angles),
    angles = angles,
    is_outlier = is.outlier,
    stringsAsFactors = FALSE
  )
}

generate.dgp.data <- function(dgp.id, cfg, sigma, seed) {
  if (identical(dgp.id, "circle_random")) {
    return(generate.circle.random(
      n = cfg$n_points,
      radius = cfg$radius,
      noise = sigma,
      noise.type = cfg$noise_type,
      seed = seed
    ))
  }
  if (identical(dgp.id, "circle_nonuniform")) {
    return(generate.circle.nonuniform(
      n = cfg$n_points,
      radius = cfg$radius,
      noise = sigma,
      noise.type = cfg$noise_type,
      seed = seed,
      arc1.weight = cfg$nonuniform_arc1_weight,
      arc2.weight = cfg$nonuniform_arc2_weight,
      arc.halfwidth = cfg$nonuniform_arc_halfwidth
    ))
  }
  if (identical(dgp.id, "circle_outliers")) {
    return(generate.circle.outliers(
      n = cfg$n_points,
      radius = cfg$radius,
      noise = sigma,
      noise.type = cfg$noise_type,
      seed = seed,
      outlier.frac = cfg$outlier_fraction,
      outlier.radial.shift = cfg$outlier_radial_shift,
      outlier.radial.sd = cfg$outlier_radial_sd
    ))
  }
  stop("Unsupported DGP id: ", dgp.id, call. = FALSE)
}

build.scenario.grid <- function(cfg) {
  sg <- expand.grid(
    dgp_id = cfg$dgp_grid,
    sigma = as.double(cfg$sigma_grid),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  sg$n_replicates <- ifelse(
    abs(sg$sigma - cfg$edge_sigma) <= 1e-12,
    cfg$n_replicates_edge,
    cfg$n_replicates_main
  )
  sg <- sg[order(match(sg$dgp_id, cfg$dgp_grid), sg$sigma), ]
  rownames(sg) <- NULL
  sg
}

build.cv.splits <- function(X, rep.seed, n.repeats, n.folds) {
  n <- nrow(X)
  split.count <- as.integer(n.repeats * n.folds)
  splits <- vector("list", split.count)
  ptr <- 0L
  for (repeat.id in seq_len(n.repeats)) {
    set.seed(as.integer(rep.seed + repeat.id * 1000L))
    fold.id <- sample(rep(seq_len(n.folds), length.out = n))
    for (fold.idx in seq_len(n.folds)) {
      ptr <- ptr + 1L
      test.idx <- which(fold.id == fold.idx)
      train.idx <- which(fold.id != fold.idx)
      splits[[ptr]] <- list(
        repeat_id = as.integer(repeat.id),
        fold_id = as.integer(fold.idx),
        train.idx = as.integer(train.idx),
        test.idx = as.integer(test.idx),
        col.scale = safe.col.scales(X[train.idx, , drop = FALSE])
      )
    }
  }
  splits
}

sqdist2 <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * (A %*% t(B))
  pmax(d2, 0)
}

predict.knn.mean <- function(X.train, ord.idx, k) {
  X.train <- as.matrix(X.train)
  k.use <- min(max(1L, as.integer(round(k))), nrow(X.train))
  if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)
  n.test <- nrow(ord.idx)
  pred <- matrix(NA_real_, nrow = n.test, ncol = ncol(X.train))
  for (i in seq_len(n.test)) {
    idx <- ord.idx[i, seq_len(k.use)]
    pred[i, ] <- colMeans(X.train[idx, , drop = FALSE])
  }
  pred
}

predict.gaussian.kernel <- function(X.train, d2, bandwidth) {
  X.train <- as.matrix(X.train)
  h <- max(as.numeric(bandwidth), 1e-6)
  W <- exp(-d2 / (2 * h * h))
  sw <- rowSums(W)
  bad <- !is.finite(sw) | sw <= 1e-12
  sw[bad] <- 1.0
  pred <- sweep(W %*% X.train, 1L, sw, "/")
  if (any(bad)) {
    nn.idx <- max.col(-d2, ties.method = "first")
    pred[bad, ] <- X.train[nn.idx[bad], , drop = FALSE]
  }
  pred
}

predict.circular.kernel <- function(X.train, dtheta, angle.bandwidth) {
  X.train <- as.matrix(X.train)
  h <- max(as.numeric(angle.bandwidth), 1e-6)
  W <- exp(-(dtheta^2) / (2 * h * h))
  sw <- rowSums(W)
  bad <- !is.finite(sw) | sw <= 1e-12
  sw[bad] <- 1.0
  pred <- sweep(W %*% X.train, 1L, sw, "/")
  if (any(bad)) {
    nn.idx <- max.col(-dtheta, ties.method = "first")
    pred[bad, ] <- X.train[nn.idx[bad], , drop = FALSE]
  }
  pred
}

evaluate.classical.cv <- function(X, splits, classical.grid) {
  X <- as.matrix(X)
  method.spec <- list(
    classical_knn_mean_cv = list(param_name = "k", param_grid = as.double(classical.grid$knn_k_grid)),
    classical_gaussian_kernel_cv = list(param_name = "bw_mult", param_grid = as.double(classical.grid$gaussian_bw_mult_grid)),
    classical_circular_kernel_cv = list(param_name = "angle_bw", param_grid = as.double(classical.grid$circular_bw_grid))
  )

  candidate.out <- list()
  selection.out <- list()
  cand.ptr <- 0L
  sel.ptr <- 0L

  for (method.id in names(method.spec)) {
    grid <- method.spec[[method.id]]$param_grid
    param.name <- method.spec[[method.id]]$param_name
    if (length(grid) == 0L) next

    scaled.mat <- matrix(NA_real_, nrow = length(splits), ncol = length(grid))
    raw.mat <- matrix(NA_real_, nrow = length(splits), ncol = length(grid))
    w.vec <- numeric(length(splits))

    for (sidx in seq_along(splits)) {
      split <- splits[[sidx]]
      X.train <- X[split$train.idx, , drop = FALSE]
      X.test <- X[split$test.idx, , drop = FALSE]
      w.vec[sidx] <- nrow(X.test)

      d2 <- sqdist2(X.test, X.train)
      ord.idx <- t(apply(d2, 1L, order))
      if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)

      d2.train <- sqdist2(X.train, X.train)
      diag(d2.train) <- Inf
      nn.dist <- sqrt(apply(d2.train, 1L, min))
      base.nn <- stats::median(nn.dist[is.finite(nn.dist)], na.rm = TRUE)
      if (!is.finite(base.nn) || base.nn <= 0) base.nn <- 0.25

      theta.train <- atan2(X.train[, 2L], X.train[, 1L])
      theta.test <- atan2(X.test[, 2L], X.test[, 1L])
      dtheta <- abs(outer(theta.test, theta.train, "-"))
      dtheta <- pmin(dtheta, 2 * pi - dtheta)

      for (pidx in seq_along(grid)) {
        param.val <- grid[pidx]
        pred <- switch(
          method.id,
          classical_knn_mean_cv = predict.knn.mean(X.train = X.train, ord.idx = ord.idx, k = param.val),
          classical_gaussian_kernel_cv = predict.gaussian.kernel(X.train = X.train, d2 = d2, bandwidth = param.val * base.nn),
          classical_circular_kernel_cv = predict.circular.kernel(X.train = X.train, dtheta = dtheta, angle.bandwidth = param.val),
          stop("Unsupported method: ", method.id, call. = FALSE)
        )

        resid <- X.test - pred
        resid.scaled <- sweep(resid, 2L, split$col.scale, "/")
        scaled.mat[sidx, pidx] <- mean(resid.scaled^2)
        raw.mat[sidx, pidx] <- mean(resid^2)
      }
    }

    cand <- data.table(
      method = method.id,
      param_name = param.name,
      param_value = grid,
      scaled_mse_mean = vapply(seq_along(grid), function(j) weighted.mean.safe(scaled.mat[, j], w.vec), numeric(1)),
      scaled_mse_se = vapply(seq_along(grid), function(j) weighted.se.safe(scaled.mat[, j], w.vec), numeric(1)),
      raw_mse_mean = vapply(seq_along(grid), function(j) weighted.mean.safe(raw.mat[, j], w.vec), numeric(1)),
      n_eval_rows = vapply(seq_along(grid), function(j) sum(is.finite(scaled.mat[, j])), integer(1))
    )
    cand <- cand[is.finite(scaled_mse_mean)][order(scaled_mse_mean, raw_mse_mean)]
    if (nrow(cand) == 0L) next

    cand.ptr <- cand.ptr + 1L
    candidate.out[[cand.ptr]] <- cand

    best <- cand[1L]
    sel.ptr <- sel.ptr + 1L
    selection.out[[sel.ptr]] <- data.table(
      method = method.id,
      selected = "min",
      param_name = best$param_name,
      param_value = best$param_value,
      scaled_mse_mean = best$scaled_mse_mean,
      scaled_mse_se = best$scaled_mse_se,
      raw_mse_mean = best$raw_mse_mean
    )
  }

  list(
    candidates = if (cand.ptr > 0L) rbindlist(candidate.out[seq_len(cand.ptr)], fill = TRUE) else data.table(),
    selection = if (sel.ptr > 0L) rbindlist(selection.out[seq_len(sel.ptr)], fill = TRUE) else data.table()
  )
}

predict.classical.full <- function(X, method.id, param.value) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2L) return(X)

  if (identical(method.id, "classical_knn_mean_cv")) {
    d2 <- sqdist2(X, X)
    diag(d2) <- Inf
    ord.idx <- t(apply(d2, 1L, order))
    if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)
    return(predict.knn.mean(X.train = X, ord.idx = ord.idx, k = param.value))
  }

  if (identical(method.id, "classical_gaussian_kernel_cv")) {
    d2 <- sqdist2(X, X)
    diag(d2) <- Inf
    nn.dist <- sqrt(apply(d2, 1L, min))
    base.nn <- stats::median(nn.dist[is.finite(nn.dist)], na.rm = TRUE)
    if (!is.finite(base.nn) || base.nn <= 0) base.nn <- 0.25
    h <- max(as.numeric(param.value) * base.nn, 1e-6)
    W <- exp(-d2 / (2 * h * h))
    diag(W) <- 0
    sw <- rowSums(W)
    bad <- !is.finite(sw) | sw <= 1e-12
    sw[bad] <- 1.0
    pred <- sweep(W %*% X, 1L, sw, "/")
    if (any(bad)) {
      nn.idx <- max.col(-d2, ties.method = "first")
      pred[bad, ] <- X[nn.idx[bad], , drop = FALSE]
    }
    return(pred)
  }

  if (identical(method.id, "classical_circular_kernel_cv")) {
    theta <- atan2(X[, 2L], X[, 1L])
    dtheta <- abs(outer(theta, theta, "-"))
    dtheta <- pmin(dtheta, 2 * pi - dtheta)
    h <- max(as.numeric(param.value), 1e-6)
    W <- exp(-(dtheta^2) / (2 * h * h))
    diag(W) <- 0
    sw <- rowSums(W)
    bad <- !is.finite(sw) | sw <= 1e-12
    sw[bad] <- 1.0
    pred <- sweep(W %*% X, 1L, sw, "/")
    if (any(bad)) {
      nn.idx <- max.col(-dtheta, ties.method = "first")
      pred[bad, ] <- X[nn.idx[bad], , drop = FALSE]
    }
    return(pred)
  }

  stop("Unsupported method id: ", method.id, call. = FALSE)
}

compute.oracle.metrics <- function(X.obs, X.hat, X.true, method, dgp.id, sigma, replicate,
                                   notes = NA_character_, modeled_fraction = NA_real_) {
  X.obs <- as.matrix(X.obs)
  X.hat <- as.matrix(X.hat)
  X.true <- as.matrix(X.true)
  rad.obs <- sqrt(rowSums(X.obs^2))
  rad.hat <- sqrt(rowSums(X.hat^2))
  rad.true <- sqrt(rowSums(X.true^2))
  data.table(
    dgp_id = dgp.id,
    sigma = as.double(sigma),
    replicate = as.integer(replicate),
    method = method,
    n_rows = nrow(X.hat),
    modeled_fraction = modeled_fraction,
    rmse_xy_obs = sqrt(mean((X.obs - X.true)^2)),
    rmse_xy_hat = sqrt(mean((X.hat - X.true)^2)),
    rmse_rad_obs = sqrt(mean((rad.obs - rad.true)^2)),
    rmse_rad_hat = sqrt(mean((rad.hat - rad.true)^2)),
    var_ratio_mean = mean(apply(X.hat, 2L, stats::var) / pmax(apply(X.obs, 2L, stats::var), 1e-12)),
    notes = notes
  )
}

run.classical.baselines <- function(case.tbl, cfg, classical.grid, log.msg) {
  empty.candidates <- data.table(
    dgp_id = character(),
    sigma = numeric(),
    replicate = integer(),
    method = character(),
    param_name = character(),
    param_value = numeric(),
    scaled_mse_mean = numeric(),
    scaled_mse_se = numeric(),
    raw_mse_mean = numeric(),
    n_eval_rows = integer()
  )
  empty.selection <- data.table(
    dgp_id = character(),
    sigma = numeric(),
    replicate = integer(),
    method = character(),
    selected = character(),
    param_name = character(),
    param_value = numeric(),
    scaled_mse_mean = numeric(),
    scaled_mse_se = numeric(),
    raw_mse_mean = numeric()
  )
  empty.oracle <- data.table(
    dgp_id = character(),
    sigma = numeric(),
    replicate = integer(),
    method = character(),
    n_rows = integer(),
    modeled_fraction = numeric(),
    rmse_xy_obs = numeric(),
    rmse_xy_hat = numeric(),
    rmse_rad_obs = numeric(),
    rmse_rad_hat = numeric(),
    var_ratio_mean = numeric(),
    notes = character()
  )

  if (is.null(cfg) || nrow(case.tbl) == 0L) {
    return(list(candidates = empty.candidates, selection = empty.selection, oracle = empty.oracle))
  }

  scenario.grid <- if (!is.null(cfg$scenario_grid)) as.data.table(cfg$scenario_grid) else as.data.table(build.scenario.grid(cfg))
  scenario.grid[, `:=`(
    dgp_id = as.character(dgp_id),
    sigma = as.double(sigma),
    scenario_idx = .I
  )]

  case.tbl <- unique(as.data.table(case.tbl)[, .(dgp_id = as.character(dgp_id), sigma = as.double(sigma), replicate = as.integer(replicate))])
  setorder(case.tbl, dgp_id, sigma, replicate)

  cand.list <- list()
  sel.list <- list()
  oracle.list <- list()
  cand.ptr <- 0L
  sel.ptr <- 0L
  oracle.ptr <- 0L

  total <- nrow(case.tbl)
  for (i in seq_len(total)) {
    dgp.id <- case.tbl$dgp_id[i]
    sigma <- case.tbl$sigma[i]
    sigma.val <- as.double(sigma)
    rep.idx <- case.tbl$replicate[i]

    scenario.idx <- scenario.grid[dgp_id == dgp.id & abs(sigma - sigma.val) <= 1e-12, scenario_idx][1L]
    if (!is.finite(scenario.idx)) {
      log.msg(sprintf("  Classical baseline skipped | missing scenario index for dgp=%s sigma=%.4f rep=%d",
                      dgp.id, sigma, rep.idx))
      next
    }

    rep.seed <- as.integer(cfg$base_seed + as.integer(scenario.idx) * 100000L + rep.idx)
    X.df <- tryCatch(
      generate.dgp.data(dgp.id = dgp.id, cfg = cfg, sigma = sigma, seed = rep.seed),
      error = function(e) e
    )
    if (inherits(X.df, "error")) {
      log.msg(sprintf("  Classical baseline skipped | dgp=%s sigma=%.4f rep=%d | generation error: %s",
                      dgp.id, sigma, rep.idx, conditionMessage(X.df)))
      next
    }

    X <- as.matrix(X.df[, c("x", "y"), drop = FALSE])
    X.true <- as.matrix(X.df[, c("x.noise.free", "y.noise.free"), drop = FALSE])
    splits <- build.cv.splits(X = X, rep.seed = rep.seed, n.repeats = cfg$n_repeats, n.folds = cfg$n_folds)
    cv.res <- evaluate.classical.cv(X = X, splits = splits, classical.grid = classical.grid)

    if (nrow(cv.res$candidates) > 0L) {
      cv.res$candidates[, `:=`(dgp_id = dgp.id, sigma = sigma, replicate = rep.idx)]
      setcolorder(cv.res$candidates, c("dgp_id", "sigma", "replicate", "method", "param_name", "param_value",
                                       "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean", "n_eval_rows"))
      cand.ptr <- cand.ptr + 1L
      cand.list[[cand.ptr]] <- cv.res$candidates
    }

    if (nrow(cv.res$selection) > 0L) {
      cv.res$selection[, `:=`(dgp_id = dgp.id, sigma = sigma, replicate = rep.idx)]
      setcolorder(cv.res$selection, c("dgp_id", "sigma", "replicate", "method", "selected",
                                      "param_name", "param_value", "scaled_mse_mean",
                                      "scaled_mse_se", "raw_mse_mean"))
      sel.ptr <- sel.ptr + 1L
      sel.list[[sel.ptr]] <- cv.res$selection

      for (j in seq_len(nrow(cv.res$selection))) {
        method.id <- as.character(cv.res$selection$method[j])
        param.name <- as.character(cv.res$selection$param_name[j])
        param.value <- as.double(cv.res$selection$param_value[j])

        X.hat <- tryCatch(
          predict.classical.full(X = X, method.id = method.id, param.value = param.value),
          error = function(e) e
        )
        if (inherits(X.hat, "error")) {
          log.msg(sprintf("  Classical oracle prediction failed | dgp=%s sigma=%.4f rep=%d method=%s | %s",
                          dgp.id, sigma, rep.idx, method.id, conditionMessage(X.hat)))
          next
        }

        oracle.ptr <- oracle.ptr + 1L
        oracle.list[[oracle.ptr]] <- compute.oracle.metrics(
          X.obs = X,
          X.hat = X.hat,
          X.true = X.true,
          method = method.id,
          dgp.id = dgp.id,
          sigma = sigma,
          replicate = rep.idx,
          notes = sprintf("%s=%.6f", param.name, param.value),
          modeled_fraction = 1.0
        )
      }
    }

    if ((i %% 16L) == 0L || i == total) {
      log.msg(sprintf("  Classical baseline progress %d/%d (%.1f%%)",
                      i, total, 100 * i / total))
    }
  }

  list(
    candidates = if (cand.ptr > 0L) rbindlist(cand.list[seq_len(cand.ptr)], fill = TRUE) else empty.candidates,
    selection = if (sel.ptr > 0L) rbindlist(sel.list[seq_len(sel.ptr)], fill = TRUE) else empty.selection,
    oracle = if (oracle.ptr > 0L) rbindlist(oracle.list[seq_len(oracle.ptr)], fill = TRUE) else empty.oracle
  )
}

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

  log.msg("Starting phase-2 analysis snapshot.")
  log.msg("Run dir:", run.dir)
  log.msg("Snapshot tag:", snapshot.tag)

  cfg.path <- file.path(run.dir, "config", "experiment_config.rds")
  cfg <- if (file.exists(cfg.path)) readRDS(cfg.path) else NULL

  expected.total.tasks <- NA_real_
  expected.rows.per.replicate <- NA_real_
  if (!is.null(cfg)) {
    if (!is.null(cfg$scenario_grid)) {
      sg <- as.data.frame(cfg$scenario_grid)
      expected.total.tasks <- sum(as.integer(sg$n_replicates)) * cfg$n_repeats * cfg$n_folds *
        length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
    } else {
      expected.total.tasks <- length(cfg$dgp_grid) * length(cfg$sigma_grid) * cfg$n_replicates_main *
        cfg$n_repeats * cfg$n_folds * length(cfg$k_grid) * length(cfg$n_eigenpairs_grid) * length(cfg$eta_grid)
    }
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

  num.cols.cv <- c("sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta",
                   "scaled_mse", "raw_mse", "n_components", "lcc_size", "lcc_frac",
                   "n_rows_model", "n_test_total", "n_test_modeled", "n_train_total",
                   "n_train_modeled", "modeled_fraction")
  for (cc in intersect(num.cols.cv, names(cv))) set(cv, j = cc, value = suppressWarnings(as.numeric(cv[[cc]])))
  for (cc in intersect(c("status", "error", "dgp_id", "graph_policy"), names(cv))) set(cv, j = cc, value = as.character(cv[[cc]]))

  num.cols.cand <- c("sigma", "replicate", "k", "n_eigenpairs", "eta", "scaled_mse_mean",
                     "scaled_mse_se", "raw_mse_mean", "modeled_fraction_mean", "n_eval_rows")
  for (cc in intersect(num.cols.cand, names(cand))) set(cand, j = cc, value = suppressWarnings(as.numeric(cand[[cc]])))
  for (cc in intersect(c("dgp_id"), names(cand))) set(cand, j = cc, value = as.character(cand[[cc]]))

  num.cols.sel <- c("sigma", "replicate", "k", "n_eigenpairs", "eta", "scaled_mse_mean",
                    "scaled_mse_se", "raw_mse_mean", "modeled_fraction_mean")
  for (cc in intersect(num.cols.sel, names(sel))) set(sel, j = cc, value = suppressWarnings(as.numeric(sel[[cc]])))
  for (cc in intersect(c("dgp_id", "selected"), names(sel))) set(sel, j = cc, value = as.character(sel[[cc]]))

  num.cols.oracle <- c("sigma", "replicate", "n_rows", "modeled_fraction", "rmse_xy_obs",
                       "rmse_xy_hat", "rmse_rad_obs", "rmse_rad_hat", "var_ratio_mean")
  for (cc in intersect(num.cols.oracle, names(oracle))) set(oracle, j = cc, value = suppressWarnings(as.numeric(oracle[[cc]])))
  for (cc in intersect(c("dgp_id", "method", "notes"), names(oracle))) set(oracle, j = cc, value = as.character(oracle[[cc]]))

  if (!("dgp_id" %in% names(cv))) cv[, dgp_id := "circle_random"]
  if (!("dgp_id" %in% names(cand))) cand[, dgp_id := "circle_random"]
  if (!("dgp_id" %in% names(sel))) sel[, dgp_id := "circle_random"]
  if (!("dgp_id" %in% names(oracle))) oracle[, dgp_id := "circle_random"]

  log.msg("0) Classical baseline reconstruction")
  classical.grid <- list(
    knn_k_grid = c(5L, 9L, 15L, 25L, 35L),
    gaussian_bw_mult_grid = c(0.40, 0.70, 1.00, 1.50, 2.20),
    circular_bw_grid = c(0.10, 0.20, 0.35, 0.55, 0.85)
  )
  classical.grid.tbl <- rbindlist(list(
    data.table(method = "classical_knn_mean_cv", param_name = "k", param_value = as.double(classical.grid$knn_k_grid)),
    data.table(method = "classical_gaussian_kernel_cv", param_name = "bw_mult", param_value = as.double(classical.grid$gaussian_bw_mult_grid)),
    data.table(method = "classical_circular_kernel_cv", param_name = "angle_bw", param_value = as.double(classical.grid$circular_bw_grid))
  ), fill = TRUE)
  fwrite(classical.grid.tbl, file.path(tab.dir, "classical_baseline_grid.csv"))

  classical.case.tbl <- unique(cv[, .(dgp_id, sigma, replicate)])
  classical.res <- run.classical.baselines(
    case.tbl = classical.case.tbl,
    cfg = cfg,
    classical.grid = classical.grid,
    log.msg = log.msg
  )
  classical.candidates <- classical.res$candidates
  classical.selection <- classical.res$selection
  classical.oracle <- classical.res$oracle

  fwrite(classical.candidates, file.path(tab.dir, "classical_cv_candidate_summary_by_replicate.csv"))
  fwrite(classical.selection, file.path(tab.dir, "classical_cv_selection_by_replicate.csv"))
  fwrite(classical.oracle, file.path(tab.dir, "classical_oracle_metrics_by_replicate.csv"))

  if (nrow(classical.oracle) > 0L) {
    log.msg(sprintf("  Classical oracle rows appended: %d", nrow(classical.oracle)))
  } else {
    log.msg("  Classical oracle rows appended: 0")
  }

  oracle <- rbindlist(list(oracle, classical.oracle), fill = TRUE, use.names = TRUE)
  fwrite(oracle, file.path(tab.dir, "oracle_metrics_extended_by_replicate.csv"))

  log.msg("1) Completion/integrity analysis")
  key.cols <- c("dgp_id", "sigma", "replicate", "repeat_id", "fold_id", "k", "n_eigenpairs", "eta")
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
  ), by = .(dgp_id, sigma, replicate)][order(dgp_id, sigma, replicate)]
  if (is.finite(expected.rows.per.replicate) && expected.rows.per.replicate > 0) {
    completion.by.group[, completion_fraction := rows / expected.rows.per.replicate]
  } else {
    completion.by.group[, completion_fraction := NA_real_]
  }
  fwrite(completion.by.group, file.path(tab.dir, "completion_by_dgp_sigma_replicate.csv"))

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
      title = sprintf("Phase-2 Run Completion Snapshot (%s)", run.id),
      subtitle = sprintf("Snapshot %s | completion %.2f%%", snapshot.tag, 100 * completion.frac),
      x = NULL, y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save.plot(p01, file.path(fig.dir, "01_completion_status.png"), width = 12, height = 5.5)

  log.msg("2) Failure topology by DGP")
  fail.by <- cv[, .(
    n_total = .N,
    n_fail = sum(status != "ok", na.rm = TRUE),
    n_graph_disconnected = sum(status == "graph_disconnected", na.rm = TRUE),
    n_small_lcc = sum(status == "graph_disconnected_small_lcc", na.rm = TRUE),
    n_fit_error = sum(status %in% c("fit_error", "fit_error_lcc"), na.rm = TRUE),
    n_refit_error = sum(status == "refit_error", na.rm = TRUE)
  ), by = .(dgp_id, sigma, k, n_eigenpairs)][order(dgp_id, sigma, n_eigenpairs, k)]
  fail.by[, fail_rate := n_fail / pmax(n_total, 1)]
  fwrite(fail.by, file.path(tab.dir, "failure_by_dgp_sigma_k_eigenpairs.csv"))

  p02 <- ggplot(fail.by, aes(x = k, y = 100 * fail_rate, color = factor(sigma), group = factor(sigma))) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.1) +
    facet_grid(dgp_id ~ n_eigenpairs) +
    labs(
      title = "Failure Rate vs k (by DGP and Sigma)",
      subtitle = "Failure = status != ok",
      x = "k",
      y = "Failure rate (%)",
      color = "sigma"
    ) +
    theme_bw(base_size = 10)
  save.plot(p02, file.path(fig.dir, "02_failrate_vs_k_by_dgp_sigma.png"), width = 13, height = 7.5)

  log.msg("3) Largest-component policy utilization")
  ok.rows <- cv[status == "ok"]
  lcc.util <- ok.rows[, .(
    n = .N,
    avg_modeled_fraction = mean(modeled_fraction, na.rm = TRUE),
    sd_modeled_fraction = sd(modeled_fraction, na.rm = TRUE),
    lcc_policy_share = mean(graph_policy == "largest_component", na.rm = TRUE),
    avg_lcc_frac = mean(lcc_frac, na.rm = TRUE)
  ), by = .(dgp_id, sigma, k)][order(dgp_id, sigma, k)]
  fwrite(lcc.util, file.path(tab.dir, "lcc_policy_utilization_by_dgp_sigma_k.csv"))

  p03 <- ggplot(lcc.util, aes(x = k, y = factor(sigma), fill = avg_modeled_fraction)) +
    geom_tile() +
    facet_wrap(~dgp_id, ncol = 1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), na.value = "grey90") +
    labs(
      title = "Average Modeled Fraction in Held-out Rows",
      subtitle = "1.0 = full smoothing on test fold; below 1 indicates LCC-only smoothing",
      x = "k",
      y = "sigma",
      fill = "modeled\nfraction"
    ) +
    theme_bw(base_size = 10)
  save.plot(p03, file.path(fig.dir, "03_modeled_fraction_heatmap.png"), width = 11, height = 7)

  log.msg("4) CV risk landscape")
  risk <- cand[is.finite(scaled_mse_mean), .(
    scaled_mse_mean_avg = mean(scaled_mse_mean, na.rm = TRUE),
    scaled_mse_mean_median = median(scaled_mse_mean, na.rm = TRUE),
    raw_mse_mean_avg = mean(raw_mse_mean, na.rm = TRUE),
    modeled_fraction_mean_avg = mean(modeled_fraction_mean, na.rm = TRUE),
    n_replicates = .N
  ), by = .(dgp_id, sigma, n_eigenpairs, k, eta)][order(dgp_id, sigma, n_eigenpairs, k, eta)]
  risk[, log10_eta := log10(eta)]
  fwrite(risk, file.path(tab.dir, "risk_landscape_mean_by_dgp_sigma.csv"))

  p04 <- ggplot(risk, aes(x = k, y = log10_eta, fill = scaled_mse_mean_avg)) +
    geom_tile() +
    facet_grid(dgp_id + n_eigenpairs ~ sigma) +
    scale_fill_viridis_c(option = "C", na.value = "grey90") +
    labs(
      title = "CV Risk Landscape (Phase-2)",
      subtitle = "Mean scaled CV MSE across available replicates",
      x = "k",
      y = "log10(eta)",
      fill = "mean\nscaled MSE"
    ) +
    theme_bw(base_size = 8)
  save.plot(p04, file.path(fig.dir, "04_cv_surface_heatmaps_by_dgp.png"), width = 15, height = 10)

  log.msg("5) Selection distributions")
  sel.ok <- sel[selected %in% c("min", "one_se") & is.finite(k) & is.finite(n_eigenpairs) & is.finite(eta)]
  if (nrow(sel.ok) > 0L) {
    freq.k <- sel.ok[, .N, by = .(dgp_id, sigma, selected, k)][order(dgp_id, sigma, selected, k)]
    freq.eig <- sel.ok[, .N, by = .(dgp_id, sigma, selected, n_eigenpairs)][order(dgp_id, sigma, selected, n_eigenpairs)]
    freq.eta <- sel.ok[, .N, by = .(dgp_id, sigma, selected, eta)][order(dgp_id, sigma, selected, eta)]
    fwrite(freq.k, file.path(tab.dir, "selection_freq_k_by_dgp_sigma.csv"))
    fwrite(freq.eig, file.path(tab.dir, "selection_freq_n_eigenpairs_by_dgp_sigma.csv"))
    fwrite(freq.eta, file.path(tab.dir, "selection_freq_eta_by_dgp_sigma.csv"))

    p05 <- ggplot(freq.k, aes(x = factor(k), y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_grid(dgp_id ~ sigma, scales = "free_y") +
      labs(title = "Selected k Distribution", x = "k", y = "Count", fill = "rule") +
      theme_bw(base_size = 8)
    save.plot(p05, file.path(fig.dir, "05_selected_k_distribution_by_dgp_sigma.png"), width = 14, height = 8)

    eta.levels <- sort(unique(freq.eta$eta))
    freq.eta[, eta_label := factor(sprintf("%.4g", eta), levels = sprintf("%.4g", eta.levels))]
    p06 <- ggplot(freq.eta, aes(x = eta_label, y = N, fill = selected)) +
      geom_col(position = "dodge") +
      facet_grid(dgp_id ~ sigma, scales = "free_y") +
      labs(title = "Selected eta Distribution", x = "eta", y = "Count", fill = "rule") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save.plot(p06, file.path(fig.dir, "06_selected_eta_distribution_by_dgp_sigma.png"), width = 14, height = 8)
  }

  log.msg("6) One-SE vs min tradeoff")
  if (nrow(sel.ok) > 0L) {
    min.sel <- sel.ok[selected == "min", .(
      dgp_id, sigma, replicate,
      k_min = k, n_eigenpairs_min = n_eigenpairs, eta_min = eta,
      scaled_mse_mean_min = scaled_mse_mean
    )]
    one.sel <- sel.ok[selected == "one_se", .(
      dgp_id, sigma, replicate,
      k_one_se = k, n_eigenpairs_one_se = n_eigenpairs, eta_one_se = eta,
      scaled_mse_mean_one_se = scaled_mse_mean
    )]
    tradeoff <- merge(min.sel, one.sel, by = c("dgp_id", "sigma", "replicate"), all = FALSE)
    if (nrow(tradeoff) > 0L) {
      tradeoff[, `:=`(
        delta_scaled_mse = scaled_mse_mean_one_se - scaled_mse_mean_min,
        delta_k = k_one_se - k_min,
        delta_n_eigenpairs = n_eigenpairs_one_se - n_eigenpairs_min,
        delta_log10_eta = log10(pmax(eta_one_se, 1e-15) / pmax(eta_min, 1e-15))
      )]
      fwrite(tradeoff, file.path(tab.dir, "one_se_vs_min_deltas_by_dgp_sigma.csv"))

      p07 <- ggplot(tradeoff, aes(x = delta_n_eigenpairs, y = delta_scaled_mse, color = factor(sigma))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(size = abs(delta_k) + 1), alpha = 0.75) +
        facet_wrap(~dgp_id, scales = "free") +
        labs(
          title = "One-SE vs Min Tradeoff by DGP",
          subtitle = "Point size = |delta k| + 1",
          x = "delta n.eigenpairs (one_se - min)",
          y = "delta scaled CV MSE (one_se - min)",
          color = "sigma",
          size = "|delta k|+1"
        ) +
        theme_bw(base_size = 10)
      save.plot(p07, file.path(fig.dir, "07_one_se_vs_min_tradeoff_by_dgp.png"), width = 12, height = 6)
    }
  }

  log.msg("7) Oracle method comparison")
  method.order <- c(
    "observed_no_smoothing",
    "data_smoother_default",
    "classical_knn_mean_cv",
    "classical_gaussian_kernel_cv",
    "classical_circular_kernel_cv",
    "cv_one_se",
    "cv_min"
  )
  oracle.sub <- oracle[method %in% method.order]
  if (nrow(oracle.sub) > 0L) {
    oracle.summary <- oracle.sub[, .(
      n = .N,
      rmse_xy_hat_mean = mean(rmse_xy_hat, na.rm = TRUE),
      rmse_xy_hat_sd = sd(rmse_xy_hat, na.rm = TRUE),
      rmse_rad_hat_mean = mean(rmse_rad_hat, na.rm = TRUE),
      rmse_rad_hat_sd = sd(rmse_rad_hat, na.rm = TRUE),
      var_ratio_mean_avg = mean(var_ratio_mean, na.rm = TRUE),
      var_ratio_mean_sd = sd(var_ratio_mean, na.rm = TRUE),
      modeled_fraction_mean = mean(modeled_fraction, na.rm = TRUE)
    ), by = .(dgp_id, sigma, method)]
    fwrite(oracle.summary, file.path(tab.dir, "oracle_summary_by_dgp_sigma_method.csv"))

    oracle.plot <- copy(oracle.sub)
    oracle.plot[, method := factor(method, levels = method.order)]

    p08 <- ggplot(oracle.plot, aes(x = method, y = rmse_xy_hat, fill = method)) +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_grid(dgp_id ~ sigma, scales = "free_y") +
      labs(title = "Oracle RMSE(x,y) by Method, DGP, and Sigma", x = NULL, y = "rmse_xy_hat") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p08, file.path(fig.dir, "08_oracle_rmse_xy_by_method_dgp_sigma.png"), width = 14, height = 8)

    rank.tbl <- copy(oracle.summary)
    rank.tbl[, rank_rmse_xy := frank(rmse_xy_hat_mean, ties.method = "min"), by = .(dgp_id, sigma)]
    rank.tbl[, method := factor(method, levels = method.order)]
    fwrite(rank.tbl, file.path(tab.dir, "method_rank_by_dgp_sigma.csv"))

    p13 <- ggplot(rank.tbl, aes(x = factor(sigma), y = method, fill = rank_rmse_xy)) +
      geom_tile(color = "white") +
      geom_text(aes(label = rank_rmse_xy), size = 2.7) +
      facet_wrap(~dgp_id, ncol = 1) +
      scale_fill_viridis_c(option = "C", direction = -1) +
      labs(
        title = "Method Rank by Mean Oracle RMSE(x,y)",
        subtitle = "Rank 1 = lowest mean oracle RMSE in a given (DGP, sigma) regime",
        x = "sigma",
        y = "method",
        fill = "rank"
      ) +
      theme_bw(base_size = 9)
    save.plot(p13, file.path(fig.dir, "13_method_rank_heatmap_by_dgp_sigma.png"), width = 11, height = 9)
  }

  log.msg("8) Delta vs observed baseline")
  delta.tbl <- data.table()
  if (nrow(oracle.sub) > 0L) {
    baseline <- oracle.sub[method == "observed_no_smoothing", .(
      dgp_id, sigma, replicate,
      baseline_rmse_xy_hat = rmse_xy_hat,
      baseline_rmse_rad_hat = rmse_rad_hat,
      baseline_var_ratio_mean = var_ratio_mean
    )]
    delta.tbl <- merge(
      oracle.sub[method != "observed_no_smoothing"],
      baseline,
      by = c("dgp_id", "sigma", "replicate"),
      all.x = TRUE
    )
    delta.tbl[, `:=`(
      delta_rmse_xy_hat = rmse_xy_hat - baseline_rmse_xy_hat,
      delta_rmse_rad_hat = rmse_rad_hat - baseline_rmse_rad_hat,
      delta_var_ratio_mean = var_ratio_mean - baseline_var_ratio_mean
    )]
    delta.tbl[, improves_vs_observed := delta_rmse_xy_hat < 0]
    delta.tbl[, method := factor(as.character(method), levels = method.order)]
    fwrite(delta.tbl, file.path(tab.dir, "delta_vs_observed_by_dgp_sigma_replicate.csv"))

    p09 <- ggplot(delta.tbl, aes(x = method, y = delta_rmse_xy_hat, fill = method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_boxplot(outlier.alpha = 0.25) +
      facet_grid(dgp_id ~ sigma, scales = "free_y") +
      labs(
        title = "Delta RMSE(x,y) vs Observed Baseline",
        subtitle = "Negative is better than observed/no-smoothing",
        x = NULL,
        y = "delta rmse_xy_hat"
      ) +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save.plot(p09, file.path(fig.dir, "09_delta_vs_observed_by_dgp_sigma.png"), width = 14, height = 8)

    win.tbl <- delta.tbl[, .(
      n = .N,
      win_rate = mean(improves_vs_observed, na.rm = TRUE),
      median_delta_rmse_xy = median(delta_rmse_xy_hat, na.rm = TRUE)
    ), by = .(dgp_id, sigma, method)][order(dgp_id, sigma, method)]
    fwrite(win.tbl, file.path(tab.dir, "smooth_vs_nosmooth_win_rate_by_dgp_sigma_method.csv"))

    delta.mean.tbl <- delta.tbl[, .(
      n = .N,
      mean_delta_rmse_xy = mean(delta_rmse_xy_hat, na.rm = TRUE),
      median_delta_rmse_xy = median(delta_rmse_xy_hat, na.rm = TRUE),
      win_rate = mean(improves_vs_observed, na.rm = TRUE)
    ), by = .(dgp_id, sigma, method)][order(dgp_id, sigma, method)]
    fwrite(delta.mean.tbl, file.path(tab.dir, "mean_delta_vs_observed_by_dgp_sigma_method.csv"))

    p10 <- ggplot(win.tbl, aes(x = sigma, y = 100 * win_rate, color = method)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.8) +
      facet_wrap(~dgp_id, ncol = 1) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(
        title = "How Often Smoothing Beats No Smoothing",
        subtitle = "Share of replicates with lower oracle RMSE(x,y) than observed baseline",
        x = "sigma",
        y = "Win rate (%)",
        color = "method"
      ) +
      theme_bw(base_size = 10)
    save.plot(p10, file.path(fig.dir, "10_smoothing_win_rate_vs_observed.png"), width = 11, height = 7)

    p14 <- ggplot(delta.mean.tbl, aes(x = factor(sigma), y = method, fill = mean_delta_rmse_xy)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.3f", mean_delta_rmse_xy)), size = 2.6) +
      facet_wrap(~dgp_id, ncol = 1) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(
        title = "Mean Delta RMSE(x,y) vs Observed Baseline",
        subtitle = "Negative values indicate improvement over no-smoothing",
        x = "sigma",
        y = "method",
        fill = "mean delta"
      ) +
      theme_bw(base_size = 9)
    save.plot(p14, file.path(fig.dir, "14_mean_delta_heatmap_by_method_dgp_sigma.png"), width = 11, height = 9)
  }

  log.msg("9) CV vs oracle agreement")
  if (nrow(sel.ok) > 0L && nrow(oracle.sub) > 0L) {
    sel.map <- sel.ok[, .(
      dgp_id, sigma, replicate, selected,
      cv_scaled_mse = scaled_mse_mean,
      cv_raw_mse = raw_mse_mean
    )]
    sel.map[, method := ifelse(selected == "min", "cv_min", "cv_one_se")]
    cv.or <- merge(
      sel.map,
      oracle.sub[, .(dgp_id, sigma, replicate, method, rmse_xy_hat, rmse_rad_hat, var_ratio_mean)],
      by = c("dgp_id", "sigma", "replicate", "method"),
      all.x = TRUE
    )
    fwrite(cv.or, file.path(tab.dir, "cv_vs_oracle_joined_by_dgp_sigma.csv"))

    corr.tbl <- cv.or[is.finite(cv_scaled_mse) & is.finite(rmse_xy_hat), .(
      n = .N,
      spearman_rho = if (.N > 1L) cor(cv_scaled_mse, rmse_xy_hat, method = "spearman") else NA_real_,
      pearson_r = if (.N > 1L) cor(cv_scaled_mse, rmse_xy_hat, method = "pearson") else NA_real_
    ), by = .(dgp_id, sigma, selected)]
    fwrite(corr.tbl, file.path(tab.dir, "cv_vs_oracle_correlation_by_dgp_sigma.csv"))

    p11 <- ggplot(cv.or, aes(x = cv_scaled_mse, y = rmse_xy_hat, color = factor(sigma), shape = selected)) +
      geom_point(alpha = 0.75) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.55) +
      facet_grid(selected ~ dgp_id, scales = "free") +
      labs(
        title = "CV Loss vs Oracle Error",
        subtitle = "Association check for proxy-risk validity",
        x = "CV scaled MSE (selected candidate)",
        y = "oracle rmse_xy_hat",
        color = "sigma",
        shape = "rule"
      ) +
      theme_bw(base_size = 9)
    save.plot(p11, file.path(fig.dir, "11_cv_vs_oracle_scatter_by_dgp.png"), width = 13, height = 8)
  }

  log.msg("10) Oversmoothing frontier")
  if (nrow(delta.tbl) > 0L) {
    frontier <- delta.tbl[, .(
      dgp_id, sigma, replicate, method,
      var_ratio_mean,
      improvement_rmse_xy = baseline_rmse_xy_hat - rmse_xy_hat,
      improvement_rmse_rad = baseline_rmse_rad_hat - rmse_rad_hat
    )]
    fwrite(frontier, file.path(tab.dir, "oversmoothing_frontier_points_by_dgp_sigma.csv"))

    p12 <- ggplot(frontier, aes(x = var_ratio_mean, y = improvement_rmse_xy, color = method)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.7) +
      facet_grid(dgp_id ~ sigma, scales = "free") +
      labs(
        title = "Oversmoothing Frontier by DGP and Sigma",
        subtitle = "Higher improvement is better; low variance ratio indicates stronger smoothing",
        x = "var_ratio_mean",
        y = "improvement in rmse_xy (vs observed baseline)",
        color = "method"
      ) +
      theme_bw(base_size = 8)
    save.plot(p12, file.path(fig.dir, "12_oversmoothing_frontier_by_dgp_sigma.png"), width = 14, height = 8)
  }

  log.msg("11) Classical parameter-selection stability")
  if (nrow(classical.selection) > 0L) {
    classical.freq <- classical.selection[
      is.finite(param_value),
      .N,
      by = .(dgp_id, sigma, method, param_name, param_value)
    ][order(method, dgp_id, sigma, param_value)]
    classical.freq[, param_label := sprintf("%.4g", param_value)]
    fwrite(classical.freq, file.path(tab.dir, "classical_selection_freq_by_dgp_sigma_method.csv"))

    p15 <- ggplot(classical.freq, aes(x = param_label, y = N, fill = factor(sigma))) +
      geom_col(position = "dodge") +
      facet_grid(method ~ dgp_id, scales = "free_x") +
      labs(
        title = "Classical Baseline Selected Parameter Frequencies",
        subtitle = "Counts of CV-selected parameter values across replicates",
        x = "selected parameter value",
        y = "count",
        fill = "sigma"
      ) +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save.plot(p15, file.path(fig.dir, "15_classical_selected_parameter_frequency.png"), width = 14, height = 8)
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
