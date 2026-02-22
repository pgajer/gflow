.gflow.prepare.spline.data <- function(x, y, w = NULL) {
    keep <- is.finite(x) & is.finite(y)
    if (!is.null(w)) {
        keep <- keep & is.finite(w) & (w >= 0)
    }

    x <- as.numeric(x[keep])
    y <- as.numeric(y[keep])
    if (!is.null(w)) w <- as.numeric(w[keep])

    if (length(x) == 0L) {
        return(list(x = numeric(0), y = numeric(0), w = if (is.null(w)) NULL else numeric(0)))
    }

    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (!is.null(w)) w <- w[ord]

    if (anyDuplicated(x)) {
        groups <- split(seq_along(x), x)
        x.u <- as.numeric(names(groups))
        if (is.null(w)) {
            y.u <- vapply(groups, function(idx) mean(y[idx]), numeric(1))
            w.u <- NULL
        } else {
            y.u <- vapply(groups, function(idx) stats::weighted.mean(y[idx], w[idx]), numeric(1))
            w.u <- vapply(groups, function(idx) sum(w[idx]), numeric(1))
        }
        x <- x.u
        y <- y.u
        w <- w.u
    }

    list(x = x, y = y, w = w)
}

.gflow.safe.weighted.mean <- function(x, w = NULL) {
    if (is.null(w)) {
        return(mean(x, na.rm = TRUE))
    }
    w <- as.numeric(w)
    if (length(w) != length(x) || sum(w, na.rm = TRUE) <= 0) {
        return(mean(x, na.rm = TRUE))
    }
    stats::weighted.mean(x, w = w, na.rm = TRUE)
}

.gflow.select.spline.spar <- function(x,
                                      y,
                                      w = NULL,
                                      spar.grid = NULL,
                                      cv.folds = 5L,
                                      cv.repeats = 3L,
                                      one.se = TRUE,
                                      df.max = NULL,
                                      seed = NULL) {
    dat <- .gflow.prepare.spline.data(x = x, y = y, w = w)
    x <- dat$x
    y <- dat$y
    w <- dat$w

    n <- length(x)
    if (n < 6L || length(unique(x)) < 6L) {
        return(list(
            spar = NA_real_,
            cv.table = data.frame(
                spar = numeric(0),
                mean.mse = numeric(0),
                se.mse = numeric(0),
                n.scores = integer(0)
            ),
            method = "cv_unavailable_small_n"
        ))
    }

    if (is.null(spar.grid)) {
        spar.grid <- seq(0.35, 1.25, by = 0.05)
    }
    spar.grid <- sort(unique(as.numeric(spar.grid[is.finite(spar.grid)])))
    if (length(spar.grid) == 0L) {
        return(list(
            spar = NA_real_,
            cv.table = data.frame(
                spar = numeric(0),
                mean.mse = numeric(0),
                se.mse = numeric(0),
                n.scores = integer(0)
            ),
            method = "cv_unavailable_empty_grid"
        ))
    }

    cv.folds <- max(2L, as.integer(cv.folds))
    cv.repeats <- max(1L, as.integer(cv.repeats))
    k <- min(cv.folds, max(2L, n - 1L))
    n.scores <- k * cv.repeats

    if (!is.null(seed) && is.finite(seed)) {
        set.seed(as.integer(seed))
    }

    score.mat <- matrix(NA_real_, nrow = length(spar.grid), ncol = n.scores)
    col.idx <- 1L

    for (r in seq_len(cv.repeats)) {
        fold.id <- sample(rep(seq_len(k), length.out = n), size = n, replace = FALSE)

        for (fold in seq_len(k)) {
            train.idx <- which(fold.id != fold)
            test.idx <- which(fold.id == fold)

            if (length(test.idx) == 0L || length(unique(x[train.idx])) < 4L) {
                col.idx <- col.idx + 1L
                next
            }

            for (j in seq_along(spar.grid)) {
                fit.args <- list(
                    x = x[train.idx],
                    y = y[train.idx],
                    spar = spar.grid[j]
                )
                if (!is.null(w)) {
                    fit.args$w <- w[train.idx]
                }

                fit <- tryCatch(
                    do.call(stats::smooth.spline, fit.args),
                    error = function(e) NULL
                )
                if (is.null(fit)) next

                if (!is.null(df.max) && is.finite(df.max) && fit$df > df.max) {
                    capped.args <- list(x = x[train.idx], y = y[train.idx], df = as.numeric(df.max))
                    if (!is.null(w)) capped.args$w <- w[train.idx]
                    fit <- tryCatch(
                        do.call(stats::smooth.spline, capped.args),
                        error = function(e) fit
                    )
                }

                pred <- tryCatch(
                    as.numeric(stats::predict(fit, x = x[test.idx])$y),
                    error = function(e) rep(NA_real_, length(test.idx))
                )

                err <- (y[test.idx] - pred)^2
                score.mat[j, col.idx] <- .gflow.safe.weighted.mean(err, if (is.null(w)) NULL else w[test.idx])
            }
            col.idx <- col.idx + 1L
        }
    }

    mean.mse <- rowMeans(score.mat, na.rm = TRUE)
    n.valid <- rowSums(is.finite(score.mat))
    sd.mse <- apply(score.mat, 1, stats::sd, na.rm = TRUE)
    se.mse <- sd.mse / sqrt(pmax(1L, n.valid))

    valid <- is.finite(mean.mse) & n.valid > 0L
    if (!any(valid)) {
        return(list(
            spar = NA_real_,
            cv.table = data.frame(
                spar = spar.grid,
                mean.mse = mean.mse,
                se.mse = se.mse,
                n.scores = as.integer(n.valid)
            ),
            method = "cv_failed_all_candidates"
        ))
    }

    best.idx <- which(valid)[which.min(mean.mse[valid])]
    chosen.idx <- best.idx

    if (isTRUE(one.se) && is.finite(se.mse[best.idx])) {
        cutoff <- mean.mse[best.idx] + se.mse[best.idx]
        one.se.candidates <- which(valid & mean.mse <= cutoff)
        if (length(one.se.candidates) > 0L) {
            chosen.idx <- max(one.se.candidates)
        }
    }

    list(
        spar = as.numeric(spar.grid[chosen.idx]),
        cv.table = data.frame(
            spar = spar.grid,
            mean.mse = mean.mse,
            se.mse = se.mse,
            n.scores = as.integer(n.valid)
        ),
        method = if (chosen.idx == best.idx) "cv_min" else "cv_1se"
    )
}

#' Robust Smoothing Spline Wrapper
#'
#' User-facing wrapper around [stats::smooth.spline()] with optional robust
#' repeated K-fold cross-validation over a `spar` grid.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param w Optional non-negative numeric weights.
#' @param df Optional fixed degrees of freedom passed to
#'   [stats::smooth.spline()].
#' @param spar Optional fixed smoothing parameter passed to
#'   [stats::smooth.spline()].
#' @param use.gcv Logical; if `TRUE`, use `smooth.spline` automatic GCV when
#'   neither `df` nor `spar` is fixed.
#' @param df.max Optional upper bound on effective degrees of freedom.
#' @param cv.folds Number of CV folds when `use.gcv = FALSE`.
#' @param cv.repeats Number of repeated CV splits when `use.gcv = FALSE`.
#' @param spar.grid Optional candidate `spar` grid for robust CV selection.
#' @param cv.one.se Logical; if `TRUE`, use a 1-SE rule after repeated CV.
#' @param cv.seed Optional seed for CV reproducibility.
#' @param reuse.spar Optional fixed `spar` value to reuse (e.g. across BB
#'   replicates).
#'
#' @return A `smooth.spline` fit object or `NULL` if fitting is not possible.
#'   Selection metadata is attached as `fit$gflow.selection`.
#'
#' @export
gflow.smooth.spline <- function(x,
                                y,
                                w = NULL,
                                df = NULL,
                                spar = NULL,
                                use.gcv = TRUE,
                                df.max = NULL,
                                cv.folds = 5L,
                                cv.repeats = 3L,
                                spar.grid = NULL,
                                cv.one.se = TRUE,
                                cv.seed = NULL,
                                reuse.spar = NULL) {
    dat <- .gflow.prepare.spline.data(x = x, y = y, w = w)
    x <- dat$x
    y <- dat$y
    w <- dat$w

    if (length(x) < 4L || length(unique(x)) < 4L) {
        return(NULL)
    }

    selection.method <- "gcv_auto"
    cv.table <- NULL

    fit.args <- list(x = x, y = y)
    if (!is.null(w)) fit.args$w <- w

    if (!is.null(reuse.spar) && is.finite(reuse.spar)) {
        fit.args$spar <- as.numeric(reuse.spar)
        selection.method <- "reused_spar"
    } else if (!is.null(df) && is.finite(df)) {
        fit.args$df <- as.numeric(df)
        selection.method <- "fixed_df"
    } else if (!is.null(spar) && is.finite(spar)) {
        fit.args$spar <- as.numeric(spar)
        selection.method <- "fixed_spar"
    } else if (!isTRUE(use.gcv)) {
        sel <- .gflow.select.spline.spar(
            x = x,
            y = y,
            w = w,
            spar.grid = spar.grid,
            cv.folds = cv.folds,
            cv.repeats = cv.repeats,
            one.se = cv.one.se,
            df.max = df.max,
            seed = cv.seed
        )
        cv.table <- sel$cv.table
        if (is.finite(sel$spar)) {
            fit.args$spar <- sel$spar
            selection.method <- sel$method
        } else {
            selection.method <- "cv_fallback_gcv"
        }
    }

    fit <- tryCatch(
        do.call(stats::smooth.spline, fit.args),
        error = function(e) NULL
    )

    if (is.null(fit) && !isTRUE(use.gcv)) {
        fallback.args <- list(x = x, y = y)
        if (!is.null(w)) fallback.args$w <- w
        fit <- tryCatch(
            do.call(stats::smooth.spline, fallback.args),
            error = function(e) NULL
        )
        selection.method <- "hard_fallback_gcv"
    }
    if (is.null(fit)) {
        return(NULL)
    }

    if (!is.null(df.max) && is.finite(df.max) && df.max > 1 && is.finite(fit$df) && fit$df > df.max) {
        capped.args <- list(x = x, y = y, df = as.numeric(df.max))
        if (!is.null(w)) capped.args$w <- w
        fit.capped <- tryCatch(
            do.call(stats::smooth.spline, capped.args),
            error = function(e) NULL
        )
        if (!is.null(fit.capped)) {
            fit <- fit.capped
            selection.method <- paste0(selection.method, "+df_capped")
        }
    }

    fit$gflow.selection <- list(
        method = selection.method,
        use.gcv = isTRUE(use.gcv),
        selected.spar = if (!is.null(fit$spar) && is.finite(fit$spar)) as.numeric(fit$spar)[1] else NA_real_,
        selected.df = if (!is.null(fit$df) && is.finite(fit$df)) as.numeric(fit$df)[1] else NA_real_,
        cv.table = cv.table
    )

    fit
}

# Backward-compatible internal alias.
.gflow.smooth.spline.wrapper <- gflow.smooth.spline

.gflow.safe.spline.predict <- function(x,
                                       y,
                                       xout = NULL,
                                       w = NULL,
                                       spar = NULL,
                                       df = NULL,
                                       use.gcv = TRUE,
                                       df.max = NULL,
                                       cv.folds = 5L,
                                       cv.repeats = 3L,
                                       spar.grid = NULL,
                                       cv.one.se = TRUE,
                                       cv.seed = NULL,
                                       reuse.spar = NULL) {
    dat <- .gflow.prepare.spline.data(x = x, y = y, w = w)
    x <- dat$x
    y <- dat$y
    w <- dat$w

    if (length(x) == 0L) {
        return(list(
            x = x,
            y = y,
            fit = NULL,
            yhat.in = numeric(0),
            yhat.out = numeric(0),
            method = "empty"
        ))
    }

    if (is.null(xout)) {
        xout <- x
    } else {
        xout <- as.numeric(xout)
    }

    if (length(x) == 1L) {
        yhat.out <- rep(y[1], length(xout))
        yhat.in <- y
        return(list(
            x = x,
            y = y,
            fit = NULL,
            yhat.in = yhat.in,
            yhat.out = yhat.out,
            method = "constant"
        ))
    }

    if (length(x) < 4L || length(unique(x)) < 4L) {
        yhat.out <- stats::approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
        yhat.in <- stats::approx(x = x, y = y, xout = x, rule = 2, ties = "ordered")$y
        return(list(
            x = x,
            y = y,
            fit = NULL,
            yhat.in = yhat.in,
            yhat.out = yhat.out,
            method = "approx"
        ))
    }

    fit <- gflow.smooth.spline(
        x = x,
        y = y,
        w = w,
        df = df,
        spar = spar,
        use.gcv = use.gcv,
        df.max = df.max,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        spar.grid = spar.grid,
        cv.one.se = cv.one.se,
        cv.seed = cv.seed,
        reuse.spar = reuse.spar
    )

    if (is.null(fit)) {
        yhat.out <- stats::approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
        yhat.in <- stats::approx(x = x, y = y, xout = x, rule = 2, ties = "ordered")$y
        return(list(
            x = x,
            y = y,
            fit = NULL,
            yhat.in = yhat.in,
            yhat.out = yhat.out,
            method = "approx"
        ))
    }

    yhat.out <- as.numeric(stats::predict(fit, x = xout)$y)
    yhat.in <- as.numeric(stats::predict(fit, x = x)$y)

    list(
        x = x,
        y = y,
        fit = fit,
        yhat.in = yhat.in,
        yhat.out = yhat.out,
        method = if (isTRUE(use.gcv)) "smooth.spline.gcv" else "smooth.spline.cv",
        selected.spar = if (!is.null(fit$gflow.selection$selected.spar)) fit$gflow.selection$selected.spar else NA_real_,
        selected.df = if (!is.null(fit$gflow.selection$selected.df)) fit$gflow.selection$selected.df else NA_real_
    )
}

.gflow.fit.curve.with.ci <- function(x,
                                     y,
                                     grid.size = 200L,
                                     spar = NULL,
                                     df = NULL,
                                     use.gcv = TRUE,
                                     df.max = NULL,
                                     cv.folds = 5L,
                                     cv.repeats = 3L,
                                     spar.grid = NULL,
                                     cv.one.se = TRUE,
                                     cv.seed = NULL,
                                     reuse.spar = NULL,
                                     clip = NULL) {
    x <- as.numeric(x)
    y <- as.numeric(y)

    grid.size <- as.integer(grid.size)
    if (!is.finite(grid.size) || grid.size < 10L) {
        grid.size <- max(10L, length(unique(x)))
    }

    xr <- range(x[is.finite(x)], na.rm = TRUE)
    if (!all(is.finite(xr)) || diff(xr) <= 0) {
        xgrid <- rep(if (all(is.finite(x))) x[1] else 0, grid.size)
    } else {
        xgrid <- seq(xr[1], xr[2], length.out = grid.size)
    }

    fit <- .gflow.safe.spline.predict(
        x = x,
        y = y,
        xout = xgrid,
        spar = spar,
        df = df,
        use.gcv = use.gcv,
        df.max = df.max,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        spar.grid = spar.grid,
        cv.one.se = cv.one.se,
        cv.seed = cv.seed,
        reuse.spar = reuse.spar
    )
    sigma <- stats::sd(fit$y - fit$yhat.in, na.rm = TRUE)
    if (!is.finite(sigma)) sigma <- 0

    lower <- fit$yhat.out - 1.96 * sigma
    upper <- fit$yhat.out + 1.96 * sigma

    if (!is.null(clip) && length(clip) == 2L && all(is.finite(clip))) {
        lower <- pmax(lower, clip[1])
        upper <- pmin(upper, clip[2])
    }

    list(
        xgrid = xgrid,
        gpredictions = fit$yhat.out,
        gpredictions.CrI = rbind(lower, upper),
        xg = xgrid,
        Eyg = fit$yhat.out,
        fit = fit$fit,
        fit.method = fit$method,
        selected.spar = fit$selected.spar,
        selected.df = fit$selected.df
    )
}

.gflow.fit.gcv.smoother <- function(x, y, params = list()) {
    grid.size <- params$grid.size
    if (is.null(grid.size)) {
        grid.size <- max(200L, 5L * length(unique(x)))
    }

    spar <- params$spar
    df <- params$df
    use.gcv <- if (is.null(params$use.gcv)) TRUE else isTRUE(params$use.gcv)
    df.max <- params$df.max
    cv.folds <- if (is.null(params$cv.folds)) 5L else as.integer(params$cv.folds)
    cv.repeats <- if (is.null(params$cv.repeats)) 3L else as.integer(params$cv.repeats)
    spar.grid <- params$spar.grid
    cv.one.se <- if (is.null(params$cv.one.se)) TRUE else isTRUE(params$cv.one.se)
    cv.seed <- params$cv.seed

    out <- .gflow.fit.curve.with.ci(
        x = x,
        y = y,
        grid.size = grid.size,
        spar = spar,
        df = df,
        use.gcv = use.gcv,
        df.max = df.max,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        spar.grid = spar.grid,
        cv.one.se = cv.one.se,
        cv.seed = cv.seed,
        clip = NULL
    )

    class(out) <- c("gflow_spline_fit", "list")
    out
}

.gflow.with.external.BB.spline <- function(x,
                                           y,
                                           lambda,
                                           bw = NULL,
                                           grid.size = 400L,
                                           degree = 1L,
                                           min.K = 5L,
                                           use.gcv = TRUE,
                                           df.max = NULL,
                                           cv.folds = 5L,
                                           cv.repeats = 3L,
                                           spar.grid = NULL,
                                           cv.one.se = TRUE,
                                           cv.seed = NULL) {
    x <- as.numeric(x)
    y <- as.numeric(y)

    if (!is.matrix(lambda)) {
        lambda <- matrix(lambda, ncol = 1L)
    }
    storage.mode(lambda) <- "double"

    if (nrow(lambda) != length(x)) {
        stop("'lambda' must have nrow(lambda) == length(x)")
    }

    grid.size <- as.integer(grid.size)
    if (!is.finite(grid.size) || grid.size < 10L) grid.size <- 400L

    x.range <- range(x[is.finite(x)], na.rm = TRUE)
    xgrid <- seq(x.range[1], x.range[2], length.out = grid.size)

    fixed.spar <- NULL
    if (!is.null(bw) && is.numeric(bw) && length(bw) == 1L && is.finite(bw) && bw > 0 && bw < 1) {
        fixed.spar <- as.numeric(bw)
    }

    # Select smoothing strength once and reuse across BB replicates.
    base.fit <- .gflow.safe.spline.predict(
        x = x,
        y = y,
        xout = xgrid,
        spar = fixed.spar,
        use.gcv = use.gcv,
        df.max = df.max,
        cv.folds = cv.folds,
        cv.repeats = cv.repeats,
        spar.grid = spar.grid,
        cv.one.se = cv.one.se,
        cv.seed = cv.seed
    )
    base.spar <- base.fit$selected.spar
    if (!is.finite(base.spar)) base.spar <- fixed.spar

    n.BB <- ncol(lambda)
    BB.gpredictions <- matrix(NA_real_, nrow = grid.size, ncol = n.BB)
    for (b in seq_len(n.BB)) {
        pred.b <- .gflow.safe.spline.predict(
            x = x,
            y = y,
            xout = xgrid,
            w = lambda[, b],
            spar = fixed.spar,
            use.gcv = use.gcv,
            df.max = df.max,
            cv.folds = cv.folds,
            cv.repeats = cv.repeats,
            spar.grid = spar.grid,
            cv.one.se = cv.one.se,
            cv.seed = cv.seed,
            reuse.spar = base.spar
        )
        BB.gpredictions[, b] <- pred.b$yhat.out
    }

    list(
        xgrid = xgrid,
        gpredictions = base.fit$yhat.out,
        BB.gpredictions = BB.gpredictions,
        opt.bw = if (is.null(base.spar) || !is.finite(base.spar)) NA_real_ else base.spar,
        xg = xgrid,
        Eyg = base.fit$yhat.out,
        degree = degree,
        min.K = min.K,
        fit.method = base.fit$method
    )
}

.gflow.fit.iknn.trend <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]

    if (length(x) < 4L || length(unique(x)) < 4L) {
        return(NULL)
    }

    xout <- sort(unique(x))
    pred <- .gflow.safe.spline.predict(x = x, y = y, xout = xout)
    yhat <- pred$yhat.out

    breakpoint <- NA_real_
    if (length(xout) >= 3L) {
        dx <- diff(xout)
        dx[dx <= 0] <- 1
        d1 <- diff(yhat) / dx
        if (length(d1) >= 2L) {
            x.mid <- (xout[-1] + xout[-length(xout)]) / 2
            dmid <- diff(x.mid)
            dmid[dmid <= 0] <- 1
            d2 <- diff(d1) / dmid
            if (length(d2) > 0L && any(is.finite(d2))) {
                idx <- which.max(abs(d2))[1]
                breakpoint <- xout[idx + 1L]
            }
        }
    }

    structure(
        list(
            x = xout,
            y = yhat,
            breakpoint = breakpoint,
            method = pred$method
        ),
        class = c("gflow_trend_fit", "list")
    )
}
