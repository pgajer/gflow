#' Plot trajectory-wise hurdle association between a transformed response and y.hat
#'
#' @description
#' Visualizes the association between a transformed response \code{z.tr} and a
#' trajectory coordinate \code{y.hat} across one or more trajectories (paths),
#' with special handling of zeros / non-detections in \code{z} via a detection
#' indicator \code{det = I(z > eps)}.
#'
#' @details
#' The function first constructs long-form per-sample data using
#' \code{prep.trajectory.hurdle.association.data()} (one row per retained sample
#' within each selected trajectory). It then provides two display modes:
#'
#' \enumerate{
#' \item \code{display = "grid"}: small multiples (one panel per trajectory),
#' showing carrier and noncarrier points (optional) and optionally adding OLS
#' and Huber M-estimator lines among carriers, plus per-panel Spearman summary.
#' \item \code{display = "overlay"}: a single panel overlaying trajectories as
#' faint polylines/points, optionally adding a robust median curve and an
#' additional fitted summary curve (piecewise-linear or smooth spline) to
#' communicate overall profile shape.
#' }
#'
#' In overlay mode, the summary fit is, by default, applied to the median curve
#' (recommended for clean communication). Alternatively, it may be applied to
#' pooled points across all selected trajectories.
#'
#' @param x Numeric vector (length \code{n.samples}) of response values.
#' @param y.hat Numeric vector (length \code{n.samples}) giving a trajectory
#'   coordinate or score (x-axis).
#' @param trajectories List of integer vectors; each element contains sample
#'   indices defining a trajectory.
#' @param trj.id Integer vector of trajectory identifiers (1-based) to plot.
#' @param eps Numeric scalar; detection threshold. Values with \code{z <= eps}
#'   are treated as non-detections. Default 0.
#' @param z.transform Character; one of \code{"log"}, \code{"logit"},
#'   \code{"arcsin-sqrt"}.
#' @param z.pseudo Numeric pseudo-count used by transforms requiring stabilization
#'   near 0 and/or 1 (notably \code{"log"} and \code{"logit"}). Default \code{1e-6}.
#' @param include.noncarriers Logical; if \code{TRUE}, include noncarriers
#'   (det==0) in the plotting data. If \code{FALSE}, include carriers only.
#' @param add.ols Logical; in grid mode, add OLS line among carriers. Default TRUE.
#' @param add.robust Logical; in grid mode, add Huber M-estimator line (MASS::rlm)
#'   among carriers if MASS is available. Default TRUE.
#' @param add.spearman Logical; in grid mode, annotate each panel with Spearman
#'   rho/p among carriers. Default TRUE.
#' @param display Character; \code{"grid"} or \code{"overlay"}. Default \code{"grid"}.
#'
#' @param overlay.type Character; overlay drawing type: \code{"b"}, \code{"l"},
#'   or \code{"p"}. Default \code{"b"}.
#' @param overlay.draw Character; \code{"carriers"} or \code{"all"}. Default
#'   \code{"carriers"}.
#' @param overlay.col Color for overlay trajectories. Default \code{"grey40"}.
#' @param overlay.lwd Numeric; overlay line width. Default 1.
#' @param overlay.alpha Numeric in (0,1]; overlay alpha. Default 0.5.
#'
#' @param median.curve Logical; in overlay mode, add median curve. Default FALSE.
#' @param median.n.grid Integer; grid size for median curve interpolation. Default 101.
#' @param median.col Color for median curve. Default \code{"red"}.
#' @param median.lwd Numeric; median curve line width. Default 2.
#' @param median.lty Integer; median curve line type. Default 1.
#'
#' @param summary.fit Character; \code{"none"}, \code{"piecewise"}, or \code{"spline"}.
#'   Default \code{"none"}. Overlay mode only.
#' @param summary.fit.to Character; \code{"median"} or \code{"pooled"}.
#'   Default \code{"median"}. Overlay mode only.
#' @param summary.fit.col Color for summary fit curve. Default \code{"blue"}.
#' @param summary.fit.lwd Numeric; summary fit line width. Default 3.
#' @param summary.fit.lty Integer; summary fit line type. Default 3.
#'
#' @param piecewise.n.breaks Integer; number of breakpoints for piecewise fit.
#'   Default 1 (two segments). Overlay mode only.
#' @param piecewise.breaks Optional numeric vector of explicit breakpoints in
#'   \code{y.hat} units. If provided, overrides automatic selection.
#' @param piecewise.min.span.frac Numeric in (0,0.5); excludes candidate
#'   breakpoints too close to range ends (fraction of x-range). Default 0.1.
#' @param piecewise.grid.n Integer; number of candidate breakpoints for
#'   one-break grid search. Default 51.
#'
#' @param spline.df Optional numeric degrees of freedom for \code{smooth.spline}.
#' @param spline.spar Optional numeric smoothing parameter for \code{smooth.spline}.
#'
#' @param pch.carrier Integer; carrier point symbol. Default 19.
#' @param pch.noncarrier Integer; noncarrier point symbol. Default 1.
#' @param cex Numeric point size. Default 0.8.
#'
#' @param main Character main title. Default NULL.
#' @param xlab Character x-axis label. Default \code{"y.hat"}.
#' @param ylab Character y-axis label. Default NULL.
#'
#' @param legend.pos Legend position (passed to \code{legend()}). Default \code{"topleft"}.
#' @param show.legend Logical; show legend. Default TRUE.
#' @param legend.cex Numeric; legend text size. Default 0.8.
#' @param legend.inset Numeric; legend inset parameter. Default 0.05.
#'
#' @param ols.col Color for OLS line. Default \code{"black"}.
#' @param ols.lwd Numeric; OLS line width. Default 1.
#' @param ols.lty Integer; OLS line type. Default 1.
#' @param robust.col Color for robust line. Default \code{"black"}.
#' @param robust.lwd Numeric; robust line width. Default 1.
#' @param robust.lty Integer; robust line type. Default 2.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the long-form plotting \code{data.frame}. When overlay mode
#' computes a median curve and/or summary fit, the returned object has an
#' attribute \code{"summary.fit"} containing fit details.
#'
#' @export
plot.trajectory.hurdle.association <- function(x,
                                               y.hat,
                                               trajectories,
                                               trj.id,
                                               eps = 0,
                                               z.transform = c("log", "logit", "arcsin-sqrt"),
                                               z.pseudo = 1e-6,
                                               include.noncarriers = TRUE,
                                               add.ols = TRUE,
                                               add.robust = TRUE,
                                               add.spearman = TRUE,
                                               display = c("grid", "overlay"),
                                               overlay.type = c("b", "l", "p"),
                                               overlay.draw = c("carriers", "all"),
                                               overlay.col = "grey40",
                                               overlay.lwd = 1,
                                               overlay.alpha = 0.5,
                                               median.curve = FALSE,
                                               median.n.grid = 101L,
                                               median.col = "red",
                                               median.lwd = 2,
                                               median.lty = 1,
                                               summary.fit = c("none", "piecewise", "spline"),
                                               summary.fit.to = c("median", "pooled"),
                                               summary.fit.col = "blue",
                                               summary.fit.lwd = 3,
                                               summary.fit.lty = 3,
                                               piecewise.n.breaks = 1L,
                                               piecewise.breaks = NULL,
                                               piecewise.min.span.frac = 0.1,
                                               piecewise.grid.n = 51L,
                                               spline.df = NULL,
                                               spline.spar = NULL,
                                               pch.carrier = 19,
                                               pch.noncarrier = 1,
                                               cex = 0.8,
                                               main = NULL,
                                               xlab = "y.hat",
                                               ylab = NULL,
                                               legend.pos = "topleft",
                                               show.legend = TRUE,
                                               legend.cex = 0.8,
                                               legend.inset = 0.05,
                                               ols.col = "black",
                                               ols.lwd = 1,
                                               ols.lty = 1,
                                               robust.col = "black",
                                               robust.lwd = 1,
                                               robust.lty = 2,
                                               ...) {
  z <- x

  z.transform <- match.arg(z.transform)
  display <- match.arg(display)
  overlay.type <- match.arg(overlay.type)
  overlay.draw <- match.arg(overlay.draw)
  summary.fit <- match.arg(summary.fit)
  summary.fit.to <- match.arg(summary.fit.to)

  stopifnot(is.numeric(z))
  stopifnot(is.numeric(y.hat))
  stopifnot(length(z) == length(y.hat))
  stopifnot(is.list(trajectories))

  if (!exists("prep.trajectory.hurdle.association.data", mode = "function")) {
    stop("prep.trajectory.hurdle.association.data() not found. Please add it to the package namespace.")
  }

  dat <- prep.trajectory.hurdle.association.data(
    z = z,
    y.hat = y.hat,
    trajectories = trajectories,
    trj.id = trj.id,
    eps = eps,
    z.transform = z.transform,
    z.pseudo = z.pseudo,
    include.noncarriers = include.noncarriers
  )

  if (nrow(dat) == 0L) stop("No data available for requested trj.id.")

  if (is.null(ylab)) ylab <- paste0("z.tr (", z.transform, ")")
  if (is.null(main)) main <- paste0("Trajectory association: ", z.transform)

  trjs <- sort(unique(dat$trajectory))
  has.mass <- requireNamespace("MASS", quietly = TRUE)

  ## ------------------------------------------------------------------------
  ## Helpers (overlay mode): median curve + summary fits
  ## ------------------------------------------------------------------------
  compute.median.curve <- function(dat, trjs, use = c("carriers", "all"), n.grid = 101L) {
    use <- match.arg(use)

    dt.use <- dat
    if (use == "carriers") dt.use <- dt.use[dt.use$det == 1L, , drop = FALSE]
    if (nrow(dt.use) == 0L) return(NULL)

    x.rng <- range(dt.use$y.hat, na.rm = TRUE)
    if (!all(is.finite(x.rng)) || x.rng[2] <= x.rng[1]) return(NULL)

    x.grid <- seq(x.rng[1], x.rng[2], length.out = as.integer(n.grid))
    z.mat <- matrix(NA_real_, nrow = length(x.grid), ncol = length(trjs))

    for (j in seq_along(trjs)) {
      t <- trjs[j]
      dt.t <- dt.use[dt.use$trajectory == t, , drop = FALSE]
      if (nrow(dt.t) < 2L) next

      o <- order(dt.t$y.hat)
      x <- dt.t$y.hat[o]
      zt <- dt.t$z.tr[o]

      ok <- is.finite(x) & is.finite(zt)
      x <- x[ok]
      zt <- zt[ok]
      if (length(x) < 2L) next
      if (length(unique(x)) < 2L) next

      ap <- stats::approx(x = x, y = zt, xout = x.grid, rule = 1, ties = "ordered")
      z.mat[, j] <- ap$y
    }

    med <- apply(z.mat, 1, function(v) stats::median(v, na.rm = TRUE))
    n.ok <- apply(z.mat, 1, function(v) sum(is.finite(v)))

    keep <- is.finite(med) & (n.ok >= 2L)
    if (!any(keep)) return(NULL)

    list(x = x.grid[keep], y = med[keep], n.ok = n.ok[keep])
  }

    fit.spline <- function(x, y, xout,
                           df = NULL,
                           spar = NULL) {

        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]; y <- y[ok]
        if (length(x) < 4L) return(NULL)
        if (length(unique(x)) < 3L) return(NULL)

        fit <- gflow.smooth.spline(
            x = x,
            y = y,
            df = df,
            spar = spar,
            use.gcv = is.null(df) && is.null(spar)
        )
        if (is.null(fit)) return(NULL)

        pr <- stats::predict(fit, x = xout)

        list(
            fit = fit,
            x = pr$x,
            y = pr$y,
            spline.df = as.numeric(fit$df)[1],
            spline.spar = as.numeric(fit$spar)[1]
        )
    }

    fit.piecewise <- function(x, y, xout,
                              n.breaks = 1L,
                              breaks = NULL,
                              min.span.frac = 0.1,
                              grid.n = 51L) {

        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]; y <- y[ok]
        if (length(x) < 6L) return(NULL)
        if (length(unique(x)) < 4L) return(NULL)

        n.breaks <- as.integer(n.breaks)
        if (is.na(n.breaks) || n.breaks < 1L) n.breaks <- 1L

        x.rng <- range(x, na.rm = TRUE)
        if (!all(is.finite(x.rng)) || x.rng[2] <= x.rng[1]) return(NULL)

        ## If explicit breaks provided, use them (sorted, unique)
        if (!is.null(breaks)) {
            br <- sort(unique(as.numeric(breaks)))
            br <- br[is.finite(br)]
            br <- br[br > x.rng[1] & br < x.rng[2]]
            if (length(br) == 0L) return(NULL)
            n.breaks <- length(br)
            breaks.use <- br

        } else {
            ## Automatic selection:
            if (n.breaks == 1L) {
                ## Grid search for one breakpoint minimizing SSE
                span <- (x.rng[2] - x.rng[1])
                lo <- x.rng[1] + min.span.frac * span
                hi <- x.rng[2] - min.span.frac * span
                if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(NULL)

                cand <- seq(lo, hi, length.out = as.integer(grid.n))
                best.sse <- Inf
                best.b <- NA_real_

                for (b in cand) {
                    i1 <- x <= b
                    i2 <- x > b
                    if (sum(i1) < 3L || sum(i2) < 3L) next
                    f1 <- stats::lm(y[i1] ~ x[i1])
                    f2 <- stats::lm(y[i2] ~ x[i2])
                    rss <- sum(stats::residuals(f1)^2) + sum(stats::residuals(f2)^2)
                    if (is.finite(rss) && rss < best.sse) {
                        best.sse <- rss
                        best.b <- b
                    }
                }

                if (!is.finite(best.b)) return(NULL)
                breaks.use <- best.b

            } else {
                ## For >1 breakpoints, use quantile-based breaks (stable, dependency-free)
                probs <- seq(0, 1, length.out = n.breaks + 2L)
                br <- as.numeric(stats::quantile(x, probs = probs[-c(1, length(probs))], na.rm = TRUE))
                br <- sort(unique(br))
                br <- br[br > x.rng[1] & br < x.rng[2]]
                if (length(br) < 1L) return(NULL)
                breaks.use <- br
                n.breaks <- length(breaks.use)
            }
        }

        ## Fit separate OLS in each segment and predict on xout
        breaks.use <- sort(as.numeric(breaks.use))
        seg.cuts <- c(-Inf, breaks.use, Inf)

        yhat <- rep(NA_real_, length(xout))
        models <- vector("list", length(seg.cuts) - 1L)

        for (s in seq_len(length(seg.cuts) - 1L)) {
            lo <- seg.cuts[s]
            hi <- seg.cuts[s + 1L]
            idx <- (x > lo) & (x <= hi)
            if (sum(idx) < 3L) next

            models[[s]] <- stats::lm(y[idx] ~ x[idx])

            idx.out <- (xout > lo) & (xout <= hi)
            if (any(idx.out)) {
                yhat[idx.out] <- stats::predict(models[[s]], newdata = data.frame(x = xout[idx.out]))
            }
        }

        keep <- is.finite(yhat)
        if (!any(keep)) return(NULL)

        list(breaks = breaks.use, models = models, x = xout[keep], y = yhat[keep])
    }

    ## ------------------------------------------------------------------------
    ## Overlay mode
    ## ------------------------------------------------------------------------
    if (display == "overlay") {

        dt.draw <- dat
        if (overlay.draw == "carriers") dt.draw <- dt.draw[dt.draw$det == 1L, , drop = FALSE]
        if (nrow(dt.draw) == 0L) stop("No points to draw under overlay.draw setting.")

        xlim <- range(dt.draw$y.hat, na.rm = TRUE)
        ylim <- range(dt.draw$z.tr, na.rm = TRUE)

        plot(NA, NA,
             xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab,
             main = main)

        col.use <- grDevices::adjustcolor(overlay.col, alpha.f = overlay.alpha)

        ## Draw each trajectory as a sorted curve in y.hat
        for (t in trjs) {
            dt.t <- dt.draw[dt.draw$trajectory == t, , drop = FALSE]
            if (nrow(dt.t) == 0L) next
            o <- order(dt.t$y.hat)
            x <- dt.t$y.hat[o]
            zt <- dt.t$z.tr[o]

            if (overlay.type %in% c("l", "b")) {
                lines(x, zt, col = col.use, lwd = overlay.lwd)
            }
            if (overlay.type %in% c("p", "b")) {
                points(x, zt, pch = pch.carrier, cex = cex, col = col.use)
            }
        }

        ## Median curve (compute if requested OR needed for summary fit to median)
        summary.attr <- list()
        med <- NULL
        need.med <- isTRUE(median.curve) || (summary.fit != "none" && summary.fit.to == "median")

        if (isTRUE(need.med)) {
            med <- compute.median.curve(dat = dat, trjs = trjs, use = overlay.draw, n.grid = median.n.grid)
            if (!is.null(med) && isTRUE(median.curve)) {
                lines(med$x, med$y, col = median.col, lwd = median.lwd, lty = median.lty)
            }
            summary.attr$median.curve <- med
        }

        ## Summary fit (overlay only)
        fit.obj <- NULL
        if (summary.fit != "none") {

            if (summary.fit.to == "median") {
                if (is.null(med)) {
                    stop("summary.fit.to = 'median' requested but median curve could not be computed.")
                }
                x.fit <- med$x
                y.fit <- med$y
            } else {
                ## pooled
                x.fit <- dt.draw$y.hat
                y.fit <- dt.draw$z.tr
            }

            xout <- seq(min(x.fit, na.rm = TRUE), max(x.fit, na.rm = TRUE), length.out = 200L)

            if (summary.fit == "spline") {
                fit.obj <- fit.spline(x = x.fit, y = y.fit, xout = xout, df = spline.df, spar = spline.spar)
                if (!is.null(fit.obj)) {
                    lines(fit.obj$x, fit.obj$y, col = summary.fit.col, lwd = summary.fit.lwd, lty = summary.fit.lty)
                }
            } else if (summary.fit == "piecewise") {
                fit.obj <- fit.piecewise(x = x.fit, y = y.fit, xout = xout,
                                         n.breaks = piecewise.n.breaks,
                                         breaks = piecewise.breaks,
                                         min.span.frac = piecewise.min.span.frac,
                                         grid.n = piecewise.grid.n)
                if (!is.null(fit.obj)) {
                    lines(fit.obj$x, fit.obj$y, col = summary.fit.col, lwd = summary.fit.lwd, lty = summary.fit.lty)
                }
            }
        }

        summary.attr$summary.fit <- list(
            method = summary.fit,
            fit.to = summary.fit.to,
            fit.obj = fit.obj,
            spline.df = if (!is.null(fit.obj$spline.df)) fit.obj$spline.df else NULL,
            spline.spar = if (!is.null(fit.obj$spline.spar)) fit.obj$spline.spar else NULL
        )

        ## Overlay legend (single global legend)
        if (isTRUE(show.legend)) {
            leg <- character(0)
            col <- character(0)
            lty <- integer(0)
            lwd <- numeric(0)

            leg <- c(leg, paste0("Trajectories (n=", length(trjs), ")"))
            col <- c(col, col.use)
            lty <- c(lty, 1L)
            lwd <- c(lwd, overlay.lwd)

            if (isTRUE(median.curve) && !is.null(med)) {
                leg <- c(leg, "Median curve")
                col <- c(col, median.col)
                lty <- c(lty, median.lty)
                lwd <- c(lwd, median.lwd)
            }

            if (summary.fit != "none" && !is.null(fit.obj)) {
                leg <- c(leg, paste0("Summary fit: ", summary.fit, " (", summary.fit.to, ")"))
                col <- c(col, summary.fit.col)
                lty <- c(lty, summary.fit.lty)
                lwd <- c(lwd, summary.fit.lwd)
            }

            legend(legend.pos,
                   legend = leg,
                   col = col,
                   lty = lty,
                   lwd = lwd,
                   bty = "n",
                   inset = legend.inset,
                   cex = legend.cex)
        }

        attr(dat, "summary.fit") <- summary.attr
        return(invisible(dat))
    }

    ## ------------------------------------------------------------------------
    ## Grid mode (small multiples)
    ## ------------------------------------------------------------------------
    n.trj <- length(trjs)
    ncol <- if (n.trj <= 2L) n.trj else ceiling(sqrt(n.trj))
    nrow <- ceiling(n.trj / ncol)

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)
    par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))

    for (t in trjs) {

        dt <- dat[dat$trajectory == t, , drop = FALSE]
        is.carrier <- dt$det == 1L
        dt.car <- dt[is.carrier, , drop = FALSE]
        dt.ncar <- dt[!is.carrier, , drop = FALSE]

        plot(dt$y.hat, dt$z.tr,
             type = "n",
             main = paste0(main, " (trj ", t, ")"),
             xlab = xlab,
             ylab = ylab)

        if (nrow(dt.ncar) > 0L) {
            points(dt.ncar$y.hat, dt.ncar$z.tr, pch = pch.noncarrier, cex = cex)
        }
        if (nrow(dt.car) > 0L) {
            points(dt.car$y.hat, dt.car$z.tr, pch = pch.carrier, cex = cex)
        }

        if (nrow(dt.car) >= 3L && length(unique(dt.car$y.hat)) > 1L) {

            if (isTRUE(add.ols)) {
                fit.lm <- stats::lm(z.tr ~ y.hat, data = dt.car)
                xs <- seq(min(dt.car$y.hat, na.rm = TRUE),
                          max(dt.car$y.hat, na.rm = TRUE),
                          length.out = 100)
                ys <- stats::predict(fit.lm, newdata = data.frame(y.hat = xs))
                lines(xs, ys, lty = ols.lty, lwd = ols.lwd, col = ols.col)
            }

            if (isTRUE(add.robust) && has.mass) {
                fit.rlm <- suppressWarnings(MASS::rlm(z.tr ~ y.hat, data = dt.car, maxit = 100))
                xs <- seq(min(dt.car$y.hat, na.rm = TRUE),
                          max(dt.car$y.hat, na.rm = TRUE),
                          length.out = 100)
                ys <- as.numeric(stats::predict(fit.rlm, newdata = data.frame(y.hat = xs)))
                lines(xs, ys, lty = robust.lty, lwd = robust.lwd, col = robust.col)
            }

            if (isTRUE(add.spearman) && isTRUE(show.legend)) {
                sp <- suppressWarnings(stats::cor.test(dt.car$y.hat, dt.car$z.tr,
                                                       method = "spearman", exact = FALSE))
                rho <- as.numeric(sp$estimate)[1]
                pval <- as.numeric(sp$p.value)[1]

                legend.txt <- c(
                    paste0("k=", nrow(dt.car)),
                    paste0("Spearman rho=", sprintf("%.3f", rho)),
                    paste0("p=", signif(pval, 3))
                )
                if (isTRUE(add.ols)) legend.txt <- c(legend.txt, "OLS")
                if (isTRUE(add.robust) && has.mass) legend.txt <- c(legend.txt, "Huber rlm")

                legend(legend.pos,
                       legend = legend.txt,
                       bty = "n",
                       inset = legend.inset,
                       cex = legend.cex)
            }

        } else if (isTRUE(show.legend)) {
            legend(legend.pos,
                   legend = c(paste0("k=", nrow(dt.car)), "Too few carriers for fit"),
                   inset = legend.inset,
                   bty = "n", cex = legend.cex)
        }
    }

    invisible(dat)
}



#' Prepare trajectory-wise hurdle association plot data
#'
#' @description
#' Prepares long-form, per-sample data for visualizing the association between a
#' transformed response \code{z} and a trajectory coordinate \code{y.hat} within
#' selected trajectories.
#'
#' @details
#' For each selected trajectory, this function:
#' \enumerate{
#'   \item extracts the sample indices for that trajectory,
#'   \item computes a detection indicator \code{det = I(z > eps)},
#'   \item optionally restricts to carriers (\code{det == 1}) when
#'         \code{include.noncarriers = FALSE},
#'   \item applies \code{z.transform} to the retained \code{z} values to produce
#'         \code{z.tr}.
#' }
#'
#' The returned data are suitable for overlay or faceted visualization, and for
#' fitting trajectory-specific models of \code{z.tr ~ y.hat} among carriers.
#'
#' @param z Numeric vector of response values (length \code{n.samples}).
#' @param y.hat Numeric vector (length \code{n.samples}) giving the trajectory
#'   coordinate or score (e.g., a monotone risk score) used on the x-axis.
#' @param trajectories List of integer vectors; each element contains sample
#'   indices defining a trajectory.
#' @param trj.id Integer vector of trajectory identifiers (1-based) to include.
#'   If \code{NULL}, all trajectories are used.
#' @param eps Numeric scalar; detection threshold. Values with \code{z <= eps}
#'   are treated as non-detections. Default is 0.
#' @param z.transform Character string specifying the transform applied to \code{z}
#'   on retained points. One of \code{"log"}, \code{"logit"}, \code{"arcsin-sqrt"}.
#' @param z.pseudo Numeric pseudo-count used by transforms requiring stabilization
#'   near 0 and/or 1 (notably \code{"log"} and \code{"logit"}). Default is \code{1e-6}.
#' @param include.noncarriers Logical; if \code{TRUE}, include both detected and
#'   non-detected points in the output. If \code{FALSE}, include detected points
#'   only (\code{det == 1}). Default is \code{TRUE}.
#'
#' @return
#' A \code{data.frame} with one row per retained sample within each selected
#' trajectory and the following columns:
#' \describe{
#'   \item{trajectory}{Integer trajectory identifier (1-based).}
#'   \item{sample.index}{Integer sample index in the original vectors.}
#'   \item{y.hat}{Numeric trajectory coordinate.}
#'   \item{z}{Numeric original response.}
#'   \item{det}{Integer detection indicator (\code{1} if \code{z > eps}, else \code{0}).}
#'   \item{z.tr}{Numeric transformed response.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Example uses toy data; replace with your objects in practice
#' z <- c(0, 0.1, 0.2, 0, 0.3)
#' y.hat <- seq_along(z)
#' trajectories <- list(c(1L, 2L, 3L), c(3L, 4L, 5L))
#' prep.trajectory.hurdle.association.data(
#'   z = z,
#'   y.hat = y.hat,
#'   trajectories = trajectories,
#'   trj.id = 1L,
#'   eps = 0,
#'   z.transform = "log",
#'   include.noncarriers = FALSE
#' )
#' }
#'
#' @keywords internal
prep.trajectory.hurdle.association.data.v1 <- function(z,
                                                    y.hat,
                                                    trajectories,
                                                    trj.id = NULL,
                                                    eps = 0,
                                                    z.transform = c("log", "logit", "arcsin-sqrt"),
                                                    z.pseudo = 1e-6,
                                                    include.noncarriers = TRUE) {

    z.transform <- match.arg(z.transform)

    stopifnot(is.numeric(z))
    stopifnot(is.numeric(y.hat))
    stopifnot(length(z) == length(y.hat))
    stopifnot(is.list(trajectories))

    if (is.null(trj.id)) trj.id <- seq_along(trajectories)
    trj.id <- as.integer(trj.id)
    trj.id <- trj.id[is.finite(trj.id)]
    trj.id <- trj.id[trj.id >= 1L & trj.id <= length(trajectories)]
    trj.id <- unique(trj.id)

    ## Transform helper (vectorized)
    transform.z <- function(z.vec) {
        if (z.transform == "log") {
            return(log(z.vec + z.pseudo))
        }
        if (z.transform == "logit") {
            z.clamp <- pmin(pmax(z.vec, z.pseudo), 1 - z.pseudo)
            return(stats::qlogis(z.clamp))
        }
        if (z.transform == "arcsin-sqrt") {
            if (any(z.vec < 0, na.rm = TRUE) || any(z.vec > 1, na.rm = TRUE)) {
                stop("z.transform = 'arcsin-sqrt' requires z in [0,1].")
            }
            z.clamp <- pmin(pmax(z.vec, 0), 1)
            return(asin(sqrt(z.clamp)))
        }
        stop("Unknown z.transform")
    }

    out <- vector("list", length(trj.id))

    for (j in seq_along(trj.id)) {
        t <- trj.id[j]

        idx <- trajectories[[t]]
        idx <- idx[!is.na(idx)]
        idx <- idx[idx >= 1L & idx <= length(z)]
        if (length(idx) == 0L) next

        z.t <- z[idx]
        yh.t <- y.hat[idx]
        det <- as.integer(z.t > eps)

        keep <- rep(TRUE, length(idx))
        if (!isTRUE(include.noncarriers)) keep <- det == 1L

        idx.keep <- idx[keep]
        z.keep <- z.t[keep]
        yh.keep <- yh.t[keep]
        det.keep <- det[keep]

        z.tr <- rep(NA_real_, length(z.keep))
        if (length(z.keep) > 0L) {
            ## Apply transform to kept points (carriers-only if requested)
            z.tr <- transform.z(z.keep)
        }

        out[[j]] <- data.frame(
            trajectory = as.integer(t),
            sample.index = as.integer(idx.keep),
            y.hat = as.numeric(yh.keep),
            z = as.numeric(z.keep),
            det = as.integer(det.keep),
            z.tr = as.numeric(z.tr)
        )
    }

    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0L) {
        return(data.frame(trajectory = integer(),
                          sample.index = integer(),
                          y.hat = numeric(),
                          z = numeric(),
                          det = integer(),
                          z.tr = numeric()))
    }

    do.call(rbind, out)
}


## -------------------------------------------------------------------------
## prep.trajectory.hurdle.association.data()
##   Now supports trajectories as either:
##     - integer vectors of vertex indices, OR
##     - trajectory objects with $vertices and optional $y.hat.modified
## -------------------------------------------------------------------------

#' Prepare trajectory-wise hurdle association plot data
#'
#' @description
#' Prepares long-form, per-sample data for visualizing the association between a
#' transformed response \code{z} and a trajectory coordinate \code{y.hat} within
#' selected trajectories.
#'
#' @details
#' For each selected trajectory, this function:
#' \enumerate{
#'   \item extracts the sample indices for that trajectory,
#'   \item computes a detection indicator \code{det = I(z > eps)},
#'   \item optionally restricts to carriers (\code{det == 1}) when
#'         \code{include.noncarriers = FALSE},
#'   \item applies \code{z.transform} to the retained \code{z} values to produce
#'         \code{z.tr}.
#' }
#'
#' The returned data are suitable for overlay or faceted visualization, and for
#' fitting trajectory-specific models of \code{z.tr ~ y.hat} among carriers.
#'
#' @param z Numeric vector of response values (length \code{n.samples}).
#' @param y.hat Numeric vector (length \code{n.samples}) giving the trajectory
#'   coordinate or score (e.g., a monotone risk score) used on the x-axis.
#' @param trajectories List of integer vectors; each element contains sample
#'   indices defining a trajectory.
#' @param trj.id Integer vector of trajectory identifiers (1-based) to include.
#'   If \code{NULL}, all trajectories are used.
#' @param eps Numeric scalar; detection threshold. Values with \code{z <= eps}
#'   are treated as non-detections. Default is 0.
#' @param z.transform Character string specifying the transform applied to \code{z}
#'   on retained points. One of \code{"log"}, \code{"logit"}, \code{"arcsin-sqrt"}.
#' @param z.pseudo Numeric pseudo-count used by transforms requiring stabilization
#'   near 0 and/or 1 (notably \code{"log"} and \code{"logit"}). Default is \code{1e-6}.
#' @param include.noncarriers Logical; if \code{TRUE}, include both detected and
#'   non-detected points in the output. If \code{FALSE}, include detected points
#'   only (\code{det == 1}). Default is \code{TRUE}.
#'
#' @return
#' A \code{data.frame} with one row per retained sample within each selected
#' trajectory and the following columns:
#' \describe{
#'   \item{trajectory}{Integer trajectory identifier (1-based).}
#'   \item{sample.index}{Integer sample index in the original vectors.}
#'   \item{y.hat}{Numeric trajectory coordinate.}
#'   \item{z}{Numeric original response.}
#'   \item{det}{Integer detection indicator (\code{1} if \code{z > eps}, else \code{0}).}
#'   \item{z.tr}{Numeric transformed response.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Example uses toy data; replace with your objects in practice
#' z <- c(0, 0.1, 0.2, 0, 0.3)
#' y.hat <- seq_along(z)
#' trajectories <- list(c(1L, 2L, 3L), c(3L, 4L, 5L))
#' prep.trajectory.hurdle.association.data(
#'   z = z,
#'   y.hat = y.hat,
#'   trajectories = trajectories,
#'   trj.id = 1L,
#'   eps = 0,
#'   z.transform = "log",
#'   include.noncarriers = FALSE
#' )
#' }
#'
#' @keywords internal
prep.trajectory.hurdle.association.data <- function(z,
                                                    y.hat,
                                                    trajectories,
                                                    trj.id = NULL,
                                                    eps = 0,
                                                    z.transform = c("log", "logit", "arcsin-sqrt"),
                                                    z.pseudo = 1e-6,
                                                    include.noncarriers = TRUE) {

    z.transform <- match.arg(z.transform)

    stopifnot(is.numeric(z))
    stopifnot(is.numeric(y.hat))
    stopifnot(length(z) == length(y.hat))
    stopifnot(is.list(trajectories))

    if (is.null(trj.id)) trj.id <- seq_along(trajectories)
    trj.id <- as.integer(trj.id)
    trj.id <- trj.id[is.finite(trj.id)]
    trj.id <- trj.id[trj.id >= 1L & trj.id <= length(trajectories)]
    trj.id <- unique(trj.id)

    ## Transform helper (vectorized)
    transform.z <- function(z.vec) {
        if (z.transform == "log") {
            return(log(z.vec + z.pseudo))
        }
        if (z.transform == "logit") {
            z.clamp <- pmin(pmax(z.vec, z.pseudo), 1 - z.pseudo)
            return(stats::qlogis(z.clamp))
        }
        if (z.transform == "arcsin-sqrt") {
            if (any(z.vec < 0, na.rm = TRUE) || any(z.vec > 1, na.rm = TRUE)) {
                stop("z.transform = 'arcsin-sqrt' requires z in [0,1].")
            }
            z.clamp <- pmin(pmax(z.vec, 0), 1)
            return(asin(sqrt(z.clamp)))
        }
        stop("Unknown z.transform")
    }

    ## Extract trajectory vertices and a y.hat vector aligned to them
    ## If trajectory provides y.hat.modified, prefer it.
    get.traj.idx.and.yhat <- function(tr, y.hat.global) {

        ## Case 1: trajectory is a plain integer vector
        if (is.atomic(tr) && is.integer(tr) || is.numeric(tr)) {
            idx <- as.integer(tr)
            yh <- as.numeric(y.hat.global[idx])
            return(list(idx = idx, yh = yh))
        }

        ## Case 2: trajectory is a list object with $vertices
        if (is.list(tr) && !is.null(tr$vertices)) {
            idx <- as.integer(tr$vertices)

            if (!is.null(tr$y.hat.modified)) {
                yh.mod <- as.numeric(tr$y.hat.modified)
                if (length(yh.mod) == length(idx) && all(is.finite(yh.mod) | is.na(yh.mod))) {
                    return(list(idx = idx, yh = yh.mod))
                }
            }

            ## Fallback: global y.hat on vertices
            yh <- as.numeric(y.hat.global[idx])
            return(list(idx = idx, yh = yh))
        }

        stop("Each element of trajectories must be an integer vector or a trajectory list with $vertices.")
    }

    out <- vector("list", length(trj.id))

    for (j in seq_along(trj.id)) {
        t <- trj.id[j]

        tr <- trajectories[[t]]
        tmp <- get.traj.idx.and.yhat(tr, y.hat)
        idx <- tmp$idx
        yh.t <- tmp$yh

        idx <- idx[!is.na(idx)]
        idx <- idx[idx >= 1L & idx <= length(z)]
        if (length(idx) == 0L) next

        ## IMPORTANT: idx filtering must be applied to yh.t as well (same positions)
        ## Keep only those positions that survived idx bounds filtering.
        ## We do it by rebuilding a mask on the original trajectory positions.
        ## Safer: recompute a keep mask using original positions.
        ## Here: since we re-filtered idx directly, align by taking y.hat on filtered idx,
        ## unless we are using a modified y.hat vector.
        ##
        ## If we used modified y.hat (yh.t came from tr$y.hat.modified), it must be
        ## filtered consistently with idx. We enforce this by reconstructing using match:
        if (is.list(tr) && !is.null(tr$vertices) && !is.null(tr$y.hat.modified) &&
            length(tr$y.hat.modified) == length(tr$vertices)) {

            pos <- match(idx, as.integer(tr$vertices))
            if (any(is.na(pos))) {
                ## Fallback: safest
                yh.t <- as.numeric(y.hat[idx])
            } else {
                yh.t <- as.numeric(tr$y.hat.modified[pos])
            }
        } else {
            yh.t <- as.numeric(y.hat[idx])
        }

        z.t <- z[idx]
        det <- as.integer(z.t > eps)

        keep <- rep(TRUE, length(idx))
        if (!isTRUE(include.noncarriers)) keep <- det == 1L

        idx.keep <- idx[keep]
        z.keep <- z.t[keep]
        yh.keep <- yh.t[keep]
        det.keep <- det[keep]

        z.tr <- rep(NA_real_, length(z.keep))
        if (length(z.keep) > 0L) {
            z.tr <- transform.z(z.keep)
        }

        out[[j]] <- data.frame(
            trajectory = as.integer(t),
            sample.index = as.integer(idx.keep),
            y.hat = as.numeric(yh.keep),
            z = as.numeric(z.keep),
            det = as.integer(det.keep),
            z.tr = as.numeric(z.tr)
        )
    }

    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0L) {
        return(data.frame(trajectory = integer(),
                          sample.index = integer(),
                          y.hat = numeric(),
                          z = numeric(),
                          det = integer(),
                          z.tr = numeric()))
    }

    do.call(rbind, out)
}
