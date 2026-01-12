#' Plot transformed abundance vs y.hat for selected trajectories
#'
#' @details
#' Visualizes z.tr vs y.hat for one or more trajectories, optionally including:
#'   - carrier vs noncarrier points
#'   - OLS and/or Huber M-estimator regression lines (carriers only)
#'   - Spearman rho annotation (carriers only)
#'   - an overlay mode (single panel) for many trajectories
#'   - an optional median association curve in overlay mode
#'
#' @param z Numeric vector (length n.samples) of abundances (relative or absolute).
#' @param y.hat Numeric vector (length n.samples), monotone score along trajectory.
#' @param trajectories List of integer vectors, sample indices per trajectory.
#' @param trj.id Integer vector of trajectory ids to plot.
#' @param eps Numeric scalar; detection threshold. Default 0.
#' @param z.transform Character; "log", "logit", or "arcsin-sqrt".
#' @param z.pseudo Numeric pseudo-count for transforms requiring it. Default 1e-6.
#' @param include.noncarriers Logical; include noncarriers in plot data. Default TRUE.
#' @param add.ols Logical; add OLS line (grid mode). Default TRUE.
#' @param add.robust Logical; add Huber rlm line (grid mode). Default TRUE.
#' @param add.spearman Logical; add Spearman annotation (grid mode). Default TRUE.
#' @param display Character; "grid" or "overlay". Default "overlay".
#' @param overlay.type Character; "b", "l", or "p". Default "b".
#' @param overlay.draw Character; "carriers" or "all". Default "carriers".
#' @param overlay.col Color for overlay trajectories. Default "grey40".
#' @param overlay.lwd Numeric line width for overlay trajectories. Default 1.
#' @param overlay.alpha Numeric in (0,1]; alpha for overlay trajectories. Default 0.5.
#' @param median.curve Logical; add median curve in overlay mode. Default FALSE.
#' @param median.n.grid Integer; grid size for median curve. Default 101.
#' @param median.col Color for median curve. Default "black".
#' @param median.lwd Numeric; median curve line width. Default 2.
#' @param median.lty Integer; median curve line type. Default 1.
#' @param pch.carrier Integer point character. Default 19.
#' @param pch.noncarrier Integer point character. Default 1.
#' @param cex Numeric point size. Default 0.8.
#' @param main Character; main title. Default NULL.
#' @param xlab Character; x label. Default "y.hat".
#' @param ylab Character; y label. Default NULL.
#' @param legend.pos Character; legend position in grid mode. Default "topleft".
#' @param show.legend Logical; show legend. Default TRUE.
#' @param legend.cex Numeric; legend text size. Default 0.8.
#' @param legend.inset Numeric; legend inset parameter. Default 0.05.
#' @param ols.col Color for OLS line. Default "black".
#' @param ols.lwd Numeric; OLS line width. Default 1.
#' @param ols.lty Integer; OLS line type. Default 1.
#' @param robust.col Color for robust line. Default "black".
#' @param robust.lwd Numeric; robust line width. Default 1.
#' @param robust.lty Integer; robust line type. Default 2.
#'
#' @return Invisibly returns the long-form plot data.frame.
#'
#' @export
plot.trajectory.hurdle.association <- function(z,
                                               y.hat,
                                               trajectories,
                                               trj.id,
                                               eps = 0,
                                               z.transform = c("log", "logit", "arcsin-sqrt"),
                                               z.pseudo = 1e-6,
                                               include.noncarriers = FALSE,
                                               add.ols = FALSE,
                                               add.robust = FALSE,
                                               add.spearman = FALSE,
                                               display = c("overlay", "grid"),
                                               overlay.type = c("b", "l", "p"),
                                               overlay.draw = c("carriers", "all"),
                                               overlay.col = "grey40",
                                               overlay.lwd = 1,
                                               overlay.alpha = 0.5,
                                               median.curve = FALSE,
                                               median.n.grid = 101L,
                                               median.col = "black",
                                               median.lwd = 2,
                                               median.lty = 1,
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
                                               robust.lty = 2) {

  z.transform <- match.arg(z.transform)
  display <- match.arg(display)
  overlay.type <- match.arg(overlay.type)
  overlay.draw <- match.arg(overlay.draw)

  dat <- prep.trajectory.hurdle.association.data(
      z = z, y.hat = y.hat, trajectories = trajectories,
      trj.id = trj.id, eps = eps,
      z.transform = z.transform, z.pseudo = z.pseudo,
      include.noncarriers = include.noncarriers
  )

    if (nrow(dat) == 0L) stop("No data available for requested trj.id.")

    if (is.null(ylab)) ylab <- paste0("z.tr (", z.transform, ")")
    if (is.null(main)) main <- paste0("Trajectory association: ", z.transform)

    has.mass <- requireNamespace("MASS", quietly = TRUE)

    ## ----------------------------------------------------------------------
    ## Helper: compute median curve in overlay mode via interpolation to grid
    ## ----------------------------------------------------------------------
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

            ## Remove non-finite
            ok <- is.finite(x) & is.finite(zt)
            x <- x[ok]
            zt <- zt[ok]
            if (length(x) < 2L) next
            if (length(unique(x)) < 2L) next

            ## Interpolate without extrapolation
            ap <- stats::approx(x = x, y = zt, xout = x.grid, rule = 1, ties = "ordered")
            z.mat[, j] <- ap$y
        }

        med <- apply(z.mat, 1, function(v) stats::median(v, na.rm = TRUE))
        n.ok <- apply(z.mat, 1, function(v) sum(is.finite(v)))

        keep <- is.finite(med) & (n.ok >= 2L)
        if (!any(keep)) return(NULL)

        list(x = x.grid[keep], y = med[keep])
    }

    ## ----------------------------------------------------------------------
    ## Display modes
    ## ----------------------------------------------------------------------
    trjs <- sort(unique(dat$trajectory))

    if (display == "overlay") {

        ## Decide what to draw
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

        ## Optional median curve
        if (isTRUE(median.curve)) {
            med <- compute.median.curve(dat = dat, trjs = trjs, use = overlay.draw, n.grid = median.n.grid)
            if (!is.null(med)) {
                lines(med$x, med$y, col = median.col, lwd = median.lwd, lty = median.lty)
            }
        }

        if (isTRUE(show.legend)) {

            leg <- c(paste0("Trajectories (n=", length(trjs), ")"))
            col <- c(grDevices::adjustcolor(overlay.col, alpha.f = overlay.alpha))
            lty <- c(1)
            lwd <- c(overlay.lwd)

            if (isTRUE(median.curve)) {
                leg <- c(leg, "Median curve")
                col <- c(col, median.col)
                lty <- c(lty, median.lty)
                lwd <- c(lwd, median.lwd)
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

        return(invisible(dat))
    }

    ## ----------------------------------------------------------------------
    ## Grid mode (current behavior, with line-style controls)
    ## ----------------------------------------------------------------------
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
