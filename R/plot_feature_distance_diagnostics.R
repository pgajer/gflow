#' Logit Transform with Pseudocount Stabilization
#'
#' @description
#' Applies a numerically-stable logit transform to values in \code{[0,1]} by first
#' adding a pseudocount and rescaling away from the boundaries. This is useful for
#' transforming relative abundances or probabilities prior to linear modeling.
#'
#' Specifically, the transform is:
#' \deqn{x' = \frac{x + \epsilon}{1 + 2\epsilon}, \quad \mathrm{logit}(x') = \log\left(\frac{x'}{1-x'}\right),}
#' where \code{epsilon} is \code{pseudo.count}.
#'
#' @param x Numeric vector. Values typically in \code{[0,1]}.
#' @param pseudo.count Positive numeric scalar. Pseudocount controlling the boundary shrinkage.
#'
#' @return Numeric vector of transformed values.
#'
#' @examples
#' \dontrun{
#' y <- transform.logit(c(0, 0.1, 0.5, 0.9, 1), pseudo.count = 1e-6)
#' }
#'
#' @keywords internal
#' @noRd
transform.logit <- function(x, pseudo.count = 1e-6) {
    ## Clamp away from 0/1 with a pseudocount
    x2 <- (x + pseudo.count) / (1 + 2 * pseudo.count)
    stats::qlogis(x2)
}

#' Compute CLR Values Within a Vertex Set
#'
#' @description
#' Computes the centered log-ratio (CLR) transform within a specified vertex set.
#' For each row (vertex/sample), CLR is computed as
#' \deqn{\mathrm{clr}(x_j) = \log(x_j + \epsilon) - \frac{1}{p}\sum_{k=1}^p \log(x_k + \epsilon),}
#' where \eqn{p} is the number of features (columns) and \code{epsilon} is a pseudocount.
#'
#' This routine is intended as a compositional sensitivity step for analyses that
#' operate on relative abundance matrices.
#'
#' @param X.sub Numeric matrix. Rows correspond to vertices/samples; columns to features.
#' @param pseudo.count Positive numeric scalar. Pseudocount added before log transform.
#'
#' @return Numeric matrix of CLR values with the same dimensions as \code{X.sub}.
#'
#' @examples
#' \dontrun{
#' clr.mat <- compute.clr.disk(phi.zmb[M1.disk.res$vertices, , drop = FALSE], pseudo.count = 1e-6)
#' }
#'
#' @keywords internal
#' @noRd
compute.clr.disk <- function(X.disk, pseudo.count = 1e-6) {
    ## CLR(x_i) = log(x_i + pc) - mean_j log(x_j + pc), per row i
    log.X <- log(X.disk + pseudo.count)
    row.center <- rowMeans(log.X)
    sweep(log.X, 1, row.center, FUN = "-")
}

#' Diagnostic Plot: Abundance vs Distance (Optionally Rank-Transformed)
#'
#' @description
#' Produces a diagnostic plot of (transformed) feature abundance versus distance within a
#' specified set of vertices. The plot is computed over carriers (\code{X > eps}) and can
#' optionally display a hex-binned view when many points are present. A smooth spline trend
#' is overlaid, and a small number of high-leverage / influential points (by hat-values and
#' Cook's distance from a simple linear model) are highlighted and labeled by vertex ID.
#'
#' The distance axis can be either raw geodesic distance or rank-transformed distance
#' (\code{distance.transform="rank"}), which is often more robust to leverage induced by a
#' small number of vertices at very small distances.
#'
#' @param feature Character scalar. Feature (column) name in \code{X} to plot.
#' @param X Numeric matrix (or coercible data.frame) with vertices in rows and features in columns.
#' @param vertices Integer vector of vertex indices (1-based) to include in the plot.
#' @param dists Numeric vector of distances for \code{vertices}. If named, names are treated as
#'   vertex IDs and used to align distances; otherwise \code{dists} must be aligned to \code{vertices}.
#' @param eps Non-negative numeric scalar. Carrier threshold; points are treated as carriers if
#'   \code{X[vertex, feature] > eps}. Default is 0.
#' @param pseudo.count Positive numeric scalar. Pseudocount used for the logit transform when
#'   \code{use.clr = FALSE}. Default is \code{1e-6}.
#' @param use.clr Logical. If TRUE, plot CLR values provided via \code{clr.mat} instead of logit
#'   transformed abundances from \code{X}. Default is FALSE.
#' @param clr.mat Optional numeric matrix of CLR values aligned with \code{vertices} (rows correspond
#'   to \code{vertices} after alignment and filtering). Required if \code{use.clr = TRUE}.
#' @param main.prefix Character scalar. Prefix appended to the plot title. Default is \code{""}.
#' @param hex Logical. If TRUE and the \pkg{hexbin} package is available, draw a hex-binned plot.
#'   Otherwise, draw a standard scatter plot. Default is TRUE.
#' @param hex.bins Integer. Number of bins in the x-direction for hex-binning (\pkg{hexbin}).
#'   Default is 25.
#' @param spline.df Optional numeric. Degrees of freedom passed to \code{\link[stats]{smooth.spline}}.
#'   If NULL (default), \code{smooth.spline} selects smoothing internally.
#' @param label.n Integer. Number of high-leverage/influential points to label (from hat-values
#'   and Cook's distance). Default is 5.
#' @param distance.transform Character string. One of \code{"raw"} or \code{"rank"}.
#'   If \code{"rank"}, ranks are computed within the plotted vertex set and scaled to \code{(0, 1)}.
#' @param xlab Optional character scalar. X-axis label. If NULL (default), a label is chosen based
#'   on \code{distance.transform}.
#'
#' @return
#' Invisibly returns a list with components \code{fit}, \code{hat}, \code{cooks}, and \code{marked}
#' (vertex IDs of labeled points). Intended for interactive diagnostics.
#'
#' @examples
#' \dontrun{
#' ## Raw distance
#' plot.abundance.distance.diagnostics("Prevotella_amnii", phi.zmb,
#'                                    vertices = M1.disk.res$vertices,
#'                                    dists = M1.disk.res$dists)
#'
#' ## Rank distance
#' plot.abundance.distance.diagnostics("Prevotella_amnii", phi.zmb,
#'                                    vertices = M1.disk.res$vertices,
#'                                    dists = M1.disk.res$dists,
#'                                    distance.transform = "rank")
#' }
#'
#' @export
plot.abundance.distance.diagnostics <- function(feature,
                                               X,
                                               vertices,
                                               dists,
                                               eps = 0,
                                               pseudo.count = 1e-6,
                                               use.clr = FALSE,
                                               clr.mat = NULL,
                                               main.prefix = "",
                                               hex = TRUE,
                                               hex.bins = 25L,
                                               spline.df = NULL,
                                               label.n = 5L,
                                               distance.transform = c("raw", "rank"),
                                               xlab = NULL) {

    distance.transform <- match.arg(distance.transform)

    ## ---- align distance vector to vertices ----
    if (!is.null(names(dists))) {
        d <- as.numeric(dists[as.character(vertices)])
    } else {
        if (length(dists) != length(vertices)) stop("If dists is unnamed, it must be aligned with vertices (same length).")
        d <- as.numeric(dists)
    }

    ok <- is.finite(d)
    vertices <- as.integer(vertices[ok])
    d <- d[ok]

    if (length(vertices) < 2L) {
        plot.new()
        title(main = paste0(main.prefix, feature, " (too few vertices)"))
        return(invisible(NULL))
    }

    ## ---- distance transform (rank within displayed set) ----
    if (distance.transform == "rank") {
        d <- rank(d, ties.method = "average") / (length(d) + 1)
        if (is.null(xlab)) xlab <- "rank(distance) within sector"
    } else {
        if (is.null(xlab)) xlab <- "geodesic distance to center"
    }

    x <- X[vertices, feature]
    carrier <- x > eps

    if (sum(carrier) < 5L) {
        plot.new()
        title(main = paste0(main.prefix, feature, " (too few carriers)"))
        return(invisible(NULL))
    }

    d.use <- d[carrier]

    transform.logit <- function(xx) {
        xx2 <- (xx + pseudo.count) / (1 + 2 * pseudo.count)
        stats::qlogis(xx2)
    }

    if (isTRUE(use.clr)) {
        if (is.null(clr.mat)) stop("clr.mat must be provided when use.clr=TRUE.")
        y.use <- clr.mat[carrier, feature]
        y.lab <- "CLR abundance"
    } else {
        y.use <- transform.logit(x[carrier])
        y.lab <- "logit abundance"
    }

    ## ---- leverage diagnostics (LM on the displayed predictor) ----
    fit <- stats::lm(y.use ~ d.use)
    h <- stats::hatvalues(fit)
    ck <- stats::cooks.distance(fit)

    idx.h <- order(h, decreasing = TRUE)[seq_len(min(label.n, length(h)))]
    idx.ck <- order(ck, decreasing = TRUE)[seq_len(min(label.n, length(ck)))]
    idx.mark <- sort(unique(c(idx.h, idx.ck)))

    ## ---- plot: hex or scatter ----
    if (isTRUE(hex) && requireNamespace("hexbin", quietly = TRUE)) {
        hb <- hexbin::hexbin(d.use, y.use, xbins = as.integer(hex.bins))
        plot(hb, main = paste0(main.prefix, feature), xlab = xlab, ylab = y.lab)
    } else {
        graphics::plot(d.use, y.use,
                       pch = 16, cex = 0.7,
                       xlab = xlab, ylab = y.lab,
                       main = paste0(main.prefix, feature))
    }

    ## ---- spline trend ----
    o <- order(d.use)
    dd <- d.use[o]
    yy <- y.use[o]
    if (is.null(spline.df)) {
        sp <- stats::smooth.spline(x = dd, y = yy)
    } else {
        sp <- stats::smooth.spline(x = dd, y = yy, df = spline.df)
    }
    graphics::lines(sp, lwd = 2)

    ## ---- mark leverage/outliers ----
    graphics::points(d.use[idx.mark], y.use[idx.mark], pch = 1, cex = 1.2, lwd = 2)

    v.car <- vertices[carrier]
    labs <- as.character(v.car[idx.mark])
    graphics::text(d.use[idx.mark], y.use[idx.mark], labels = labs, pos = 4, cex = 0.7)

    graphics::mtext(sprintf("n.carriers=%d; marked=%d (hat/Cook)", sum(carrier), length(idx.mark)),
                    side = 3, line = 0.2, cex = 0.8)

    invisible(list(fit = fit, hat = h, cooks = ck, marked = v.car[idx.mark]))
}

#' Diagnostic Plot: Presence vs Distance (Optionally Rank-Transformed)
#'
#' @description
#' Produces a diagnostic plot of presence/absence (carrier indicator) versus distance within
#' a specified set of vertices. Presence is defined as \code{X[vertex, feature] > eps}.
#' The plot includes jittered binary points and a binned prevalence curve computed over
#' equal-count quantile bins of the distance axis.
#'
#' The distance axis can be either raw geodesic distance or rank-transformed distance
#' (\code{distance.transform="rank"}), which is often more robust to leverage induced by a
#' small number of vertices at very small distances.
#'
#' @param feature Character scalar. Feature (column) name in \code{X} to plot.
#' @param X Numeric matrix (or coercible data.frame) with vertices in rows and features in columns.
#' @param vertices Integer vector of vertex indices (1-based) to include in the plot.
#' @param dists Numeric vector of distances for \code{vertices}. If named, names are treated as
#'   vertex IDs and used to align distances; otherwise \code{dists} must be aligned to \code{vertices}.
#' @param eps Non-negative numeric scalar. Carrier threshold; points are treated as carriers if
#'   \code{X[vertex, feature] > eps}. Default is 0.
#' @param n.bins Integer. Number of equal-count quantile bins used to compute the binned prevalence
#'   curve. Default is 10.
#' @param main.prefix Character scalar. Prefix appended to the plot title. Default is \code{""}.
#' @param distance.transform Character string. One of \code{"raw"} or \code{"rank"}.
#'   If \code{"rank"}, ranks are computed within the plotted vertex set and scaled to \code{(0, 1)}.
#' @param xlab Optional character scalar. X-axis label. If NULL (default), a label is chosen based
#'   on \code{distance.transform}.
#'
#' @return
#' Invisibly returns a list with components \code{mid} (bin midpoints) and \code{prev}
#' (binned prevalence estimates).
#'
#' @examples
#' \dontrun{
#' plot.presence.distance.diagnostics("Dialister_sp001553355", phi.zmb,
#'                                   vertices = M1.disk.res$vertices,
#'                                   dists = M1.disk.res$dists,
#'                                   n.bins = 12)
#' }
#'
#' @export
plot.presence.distance.diagnostics <- function(feature,
                                              X,
                                              vertices,
                                              dists,
                                              eps = 0,
                                              n.bins = 10L,
                                              main.prefix = "",
                                              distance.transform = c("raw", "rank"),
                                              xlab = NULL) {

    distance.transform <- match.arg(distance.transform)

    ## ---- align distance vector to vertices ----
    if (!is.null(names(dists))) {
        d <- as.numeric(dists[as.character(vertices)])
    } else {
        if (length(dists) != length(vertices)) stop("If dists is unnamed, it must be aligned with vertices (same length).")
        d <- as.numeric(dists)
    }

    ok <- is.finite(d)
    vertices <- as.integer(vertices[ok])
    d <- d[ok]

    if (length(vertices) < 2L) {
        plot.new()
        title(main = paste0(main.prefix, feature, " (too few vertices)"))
        return(invisible(NULL))
    }

    ## ---- distance transform (rank within displayed set) ----
    if (distance.transform == "rank") {
        d <- rank(d, ties.method = "average") / (length(d) + 1)
        if (is.null(xlab)) xlab <- "rank(distance) within sector"
    } else {
        if (is.null(xlab)) xlab <- "geodesic distance to center"
    }

    x <- X[vertices, feature]
    carrier <- (x > eps)
    y <- as.integer(carrier)

    ## jittered points
    yj <- y + stats::runif(length(y), min = -0.05, max = 0.05)

    graphics::plot(d, yj,
                   pch = 16, cex = 0.7,
                   xlab = xlab,
                   ylab = "carrier indicator (jittered)",
                   main = paste0(main.prefix, feature))

    graphics::abline(h = 0, lty = 3)
    graphics::abline(h = 1, lty = 3)

    ## binned prevalence curve (equal-count bins by quantiles)
    probs <- seq(0, 1, length.out = as.integer(n.bins) + 1L)
    brks <- unique(stats::quantile(d, probs = probs, na.rm = TRUE))
    if (length(brks) < 3L) return(invisible(NULL))

    bin <- cut(d, breaks = brks, include.lowest = TRUE)
    prev <- tapply(y, bin, mean)
    mid <- tapply(d, bin, function(z) mean(range(z)))

    graphics::lines(as.numeric(mid), as.numeric(prev), lwd = 2)
    graphics::points(as.numeric(mid), as.numeric(prev), pch = 1, cex = 1.2, lwd = 2)

    graphics::mtext(sprintf("n=%d; prevalence=%.3f", length(y), mean(y)), side = 3, line = 0.2, cex = 0.8)

    invisible(list(mid = mid, prev = prev))
}

#' Diagnostic Plot for Quantile-Bin Distance Analysis
#'
#' @description
#' Visualizes the results of a distance quantile-bin analysis (equal-count bins) as either:
#' \itemize{
#'   \item \code{"abundance"}: bin midpoint vs bin abundance summary (typically median), optionally with SE bars;
#'   \item \code{"presence"}: bin midpoint vs binned prevalence, optionally with SE bars.
#' }
#' If \code{bin.res$tests} includes a trend test (\code{spearman} or \code{wls}), the plot
#' annotates the corresponding test statistic and p-value (and BH-adjusted p-value when present).
#'
#' @param bin.res A list produced by a quantile-bin analysis routine, containing at least a
#'   data.frame component \code{bin.res$bins}. For abundance plots, \code{bins} must contain
#'   \code{d.mid} and \code{abundance}; for presence plots it must contain \code{d.mid.all} and
#'   \code{presence.p}. Optional standard error columns (\code{abundance.se}, \code{presence.se})
#'   are used when present.
#' @param what Character string. One of \code{"abundance"} or \code{"presence"}.
#' @param main Optional character scalar. Plot title.
#' @param xlab Character scalar. X-axis label. Default is \code{"distance (bin midpoint)"}.
#' @param ylab Optional character scalar. Y-axis label. If NULL, a default is used based on \code{what}.
#' @param show.se Logical. If TRUE (default) and corresponding SE columns are present, draw SE bars.
#' @param add.line Logical. If TRUE (default), connect binned points with a line.
#' @param pch Numeric. Plotting character for binned points. Default is 16.
#' @param cex Numeric. Point expansion factor. Default is 1.0.
#'
#' @return
#' Invisibly returns TRUE. Called for its plotting side effects.
#'
#' @examples
#' \dontrun{
#' ## Suppose br is an entry from spt.qbin.res$bin.results[["Prevotella_amnii||1376"]]
#' plot.quantile.bin.diagnostics(br$clr, what = "abundance",
#'                              main = "Prevotella_amnii (CLR) binned trend")
#' plot.quantile.bin.diagnostics(br$clr, what = "presence",
#'                              main = "Prevotella_amnii prevalence vs distance bins")
#' }
#'
#' @export
plot.quantile.bin.diagnostics <- function(bin.res,
                                         what = c("abundance", "presence"),
                                         main = NULL,
                                         xlab = "distance (bin midpoint)",
                                         ylab = NULL,
                                         show.se = TRUE,
                                         add.line = TRUE,
                                         pch = 16,
                                         cex = 1.0) {

    what <- match.arg(what)

    if (is.null(bin.res) || is.null(bin.res$bins) || !is.data.frame(bin.res$bins)) {
        stop("bin.res must be a list with a data.frame component bin.res$bins (as returned by quantile-bin analysis).")
    }

    df <- bin.res$bins

    if (what == "abundance") {

        if (!all(c("d.mid", "abundance") %in% colnames(df))) {
            stop("bin.res$bins must contain columns d.mid and abundance for what='abundance'.")
        }

        if (is.null(ylab)) ylab <- "abundance (bin median)"

        graphics::plot(df$d.mid, df$abundance,
                       xlab = xlab, ylab = ylab,
                       main = main,
                       pch = pch, cex = cex)

        if (isTRUE(show.se) && "abundance.se" %in% colnames(df)) {
            se <- df$abundance.se
            ok <- is.finite(se)
            if (any(ok)) {
                graphics::segments(df$d.mid[ok], df$abundance[ok] - se[ok],
                                   df$d.mid[ok], df$abundance[ok] + se[ok])
            }
        }

        if (isTRUE(add.line)) {
            o <- order(df$d.mid)
            graphics::lines(df$d.mid[o], df$abundance[o], lwd = 2)
        }

        ## optional p-value annotation if present
        if (!is.null(bin.res$tests) && !is.null(bin.res$tests$abundance)) {
            tt <- bin.res$tests$abundance
            txt <- NULL
            if (!is.null(tt$method) && tt$method == "spearman" && !is.null(tt$estimate) && !is.null(tt$p.value)) {
                txt <- sprintf("Spearman rho=%.3f, p=%.3g", tt$estimate, tt$p.value)
            }
            if (!is.null(tt$method) && tt$method == "wls" && !is.null(tt$slope) && !is.null(tt$p.value)) {
                txt <- sprintf("WLS slope=%.3f, p=%.3g", tt$slope, tt$p.value)
            }
            if (!is.null(tt$p.adj) && is.finite(tt$p.adj)) {
                txt <- paste0(txt, sprintf(", BH=%.3g", tt$p.adj))
            }
            if (!is.null(txt)) graphics::mtext(txt, side = 3, line = 0.2, cex = 0.8)
        }

        return(invisible(TRUE))
    }

    ## presence
    if (!all(c("d.mid.all", "presence.p") %in% colnames(df))) {
        stop("bin.res$bins must contain columns d.mid.all and presence.p for what='presence'.")
    }

    if (is.null(ylab)) ylab <- "prevalence"

    graphics::plot(df$d.mid.all, df$presence.p,
                   xlab = xlab, ylab = ylab,
                   main = main,
                   pch = pch, cex = cex,
                   ylim = c(0, 1))

    if (isTRUE(show.se) && "presence.se" %in% colnames(df)) {
        se <- df$presence.se
        ok <- is.finite(se)
        if (any(ok)) {
            graphics::segments(df$d.mid.all[ok], df$presence.p[ok] - se[ok],
                               df$d.mid.all[ok], df$presence.p[ok] + se[ok])
        }
    }

    if (isTRUE(add.line)) {
        o <- order(df$d.mid.all)
        graphics::lines(df$d.mid.all[o], df$presence.p[o], lwd = 2)
    }

    if (!is.null(bin.res$tests) && !is.null(bin.res$tests$presence)) {
        tt <- bin.res$tests$presence
        txt <- NULL
        if (!is.null(tt$method) && tt$method == "spearman" && !is.null(tt$estimate) && !is.null(tt$p.value)) {
            txt <- sprintf("Spearman rho=%.3f, p=%.3g", tt$estimate, tt$p.value)
        }
        if (!is.null(tt$method) && tt$method == "wls" && !is.null(tt$slope) && !is.null(tt$p.value)) {
            txt <- sprintf("WLS slope=%.3f, p=%.3g", tt$slope, tt$p.value)
        }
        if (!is.null(tt$p.adj) && is.finite(tt$p.adj)) {
            txt <- paste0(txt, sprintf(", BH=%.3g", tt$p.adj))
        }
        if (!is.null(txt)) graphics::mtext(txt, side = 3, line = 0.2, cex = 0.8)
    }

    invisible(TRUE)
}
