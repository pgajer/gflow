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


#' Create a PDF of Distance Diagnostics for Multiple Features
#'
#' @description
#' Generates a multi-page PDF containing diagnostic plots for a set of features:
#' (i) abundance vs distance (logit transform) among carriers, (ii) optional CLR
#' abundance vs distance among carriers, and (iii) presence/absence vs distance
#' with a binned prevalence curve. Optionally, if quantile-bin results are supplied,
#' adds bin-summary diagnostics via \code{\link{plot.quantile.bin.diagnostics}}.
#'
#' @param pdf.file Character scalar. Path to the output PDF file.
#' @param features Character vector of feature names (columns of \code{X}).
#' @param X Numeric matrix (or coercible data.frame) with vertices in rows and features in columns.
#' @param vertices Integer vector of vertex indices (1-based) defining the plotted set (e.g., disk).
#' @param dists Numeric vector of distances for \code{vertices} (named by vertex id or aligned).
#' @param eps Non-negative numeric scalar. Carrier threshold; \code{X > eps}. Default is 0.
#' @param pseudo.count Positive numeric scalar. Pseudocount for logit transform and CLR. Default is \code{1e-6}.
#' @param do.clr Logical. If TRUE (default), compute CLR within \code{vertices} and plot CLR diagnostics.
#' @param spline.df Optional numeric. Degrees of freedom for \code{\link[stats]{smooth.spline}}.
#'   If NULL (default), smoothing is chosen internally.
#' @param distance.transform Character string. One of \code{"raw"} or \code{"rank"}.
#'   Passed to \code{plot.abundance.distance.diagnostics} and \code{plot.presence.distance.diagnostics}.
#' @param xlab Optional character scalar. If provided, used as x-axis label for all plots; if NULL,
#'   plot functions choose a default based on \code{distance.transform}.
#' @param bin.results Optional named list of quantile-bin results. If provided, the function will attempt
#'   to add bin diagnostics for each feature using \code{plot.quantile.bin.diagnostics}. This can be either:
#'   \itemize{
#'     \item a list keyed by feature name with entries that contain \code{$clr} and/or \code{$logit}, or
#'     \item a list keyed by \code{"feature||sector"} as in the SPT quantile-bin output (in which case,
#'       pass \code{bin.key.fun} to specify how to find the key).
#'   }
#' @param bin.key.fun Optional function. If \code{bin.results} is keyed by something other than \code{feature},
#'   provide a function \code{function(feature) key} returning the lookup key in \code{bin.results}.
#'   Default is NULL (use \code{feature} directly).
#' @param include.bin.plots Logical. If TRUE and \code{bin.results} is provided, include bin diagnostics.
#'   Default is TRUE.
#'
#' @return Invisibly returns TRUE.
#'
#' @export
make.disk.diagnostics.pdf <- function(pdf.file,
                                      features,
                                      X,
                                      vertices,
                                      dists,
                                      eps = 0,
                                      pseudo.count = 1e-6,
                                      do.clr = TRUE,
                                      spline.df = NULL,
                                      distance.transform = c("raw", "rank"),
                                      xlab = NULL,
                                      bin.results = NULL,
                                      bin.key.fun = NULL,
                                      include.bin.plots = TRUE) {

    distance.transform <- match.arg(distance.transform)

    if (!is.character(pdf.file) || length(pdf.file) != 1L || is.na(pdf.file)) {
        stop("pdf.file must be a single character path.")
    }
    if (length(features) == 0L) stop("No features supplied.")
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix (or coercible data.frame).")
    if (is.null(colnames(X))) stop("X must have column names (feature names).")

    vertices <- as.integer(vertices)
    vertices <- vertices[!is.na(vertices)]
    vertices <- vertices[vertices >= 1L & vertices <= nrow(X)]
    if (length(vertices) == 0L) stop("vertices is empty after filtering to valid indices.")

    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) stop("eps must be a single non-negative numeric.")
    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric.")
    }

    ## Precompute CLR matrix over vertices if requested
    clr.mat <- NULL
    if (isTRUE(do.clr)) {
        X.sub <- X[vertices, , drop = FALSE]
        clr.mat <- compute.clr.disk(X.sub, pseudo.count = pseudo.count)
        ## clr.mat rows are aligned to 'vertices' order used here
    }

    ## Key resolver for bin.results
    if (!is.null(bin.results) && isTRUE(include.bin.plots)) {
        if (!is.list(bin.results)) stop("bin.results must be a list when provided.")
        if (is.null(bin.key.fun)) {
            bin.key.fun <- function(feature) feature
        }
        if (!is.function(bin.key.fun)) stop("bin.key.fun must be a function(feature) -> key.")
    }

    grDevices::pdf(pdf.file, width = 8.5, height = 6.5)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (f in features) {

        if (!f %in% colnames(X)) {
            ## Skip missing features (defensive)
            next
        }

        ## 1) logit-abundance diagnostics among carriers
        plot.abundance.distance.diagnostics(
            feature = f,
            X = X,
            vertices = vertices,
            dists = dists,
            eps = eps,
            pseudo.count = pseudo.count,
            use.clr = FALSE,
            main.prefix = "logit: ",
            spline.df = spline.df,
            distance.transform = distance.transform,
            xlab = xlab
        )

        ## 2) CLR-abundance diagnostics among carriers (sensitivity)
        if (isTRUE(do.clr)) {
            plot.abundance.distance.diagnostics(
                feature = f,
                X = X,
                vertices = vertices,
                dists = dists,
                eps = eps,
                pseudo.count = pseudo.count,
                use.clr = TRUE,
                clr.mat = clr.mat,
                main.prefix = "CLR: ",
                spline.df = spline.df,
                distance.transform = distance.transform,
                xlab = xlab
            )
        }

        ## 3) presence/absence diagnostic
        plot.presence.distance.diagnostics(
            feature = f,
            X = X,
            vertices = vertices,
            dists = dists,
            eps = eps,
            main.prefix = "presence: ",
            distance.transform = distance.transform,
            xlab = xlab
        )

        ## 4) optional quantile-bin diagnostics
        if (!is.null(bin.results) && isTRUE(include.bin.plots)) {

            key <- bin.key.fun(f)
            br <- bin.results[[key]]

            ## Allow br to be either a bin-res directly or a list with $clr/$logit
            if (!is.null(br)) {

                ## CLR bin plots if present
                if (is.list(br) && !is.null(br$clr)) {
                    plot.quantile.bin.diagnostics(
                        bin.res = br$clr,
                        what = "abundance",
                        main = paste0("qbin CLR abundance: ", f),
                        xlab = "distance (bin midpoint)",
                        ylab = "CLR (bin median)"
                    )
                    ## presence bins (if available)
                    if (!is.null(br$clr$bins) &&
                        all(c("d.mid.all", "presence.p") %in% names(br$clr$bins))) {
                        plot.quantile.bin.diagnostics(
                            bin.res = br$clr,
                            what = "presence",
                            main = paste0("qbin presence: ", f),
                            xlab = "distance (bin midpoint)",
                            ylab = "prevalence"
                        )
                    }
                }

                ## logit bin plots if present
                if (is.list(br) && !is.null(br$logit)) {
                    plot.quantile.bin.diagnostics(
                        bin.res = br$logit,
                        what = "abundance",
                        main = paste0("qbin logit abundance: ", f),
                        xlab = "distance (bin midpoint)",
                        ylab = "logit (bin median)"
                    )
                }

                ## If br itself is a bin-res (rare), plot abundance as default
                if (!is.list(br$clr) && !is.list(br$logit) && !is.null(br$bins)) {
                    plot.quantile.bin.diagnostics(
                        bin.res = br,
                        what = "abundance",
                        main = paste0("qbin abundance: ", f)
                    )
                }
            }
        }
    }

    invisible(TRUE)
}

#' Create a PDF of SPT Sector Diagnostics for Feature–Distance Associations
#'
#' @description
#' Generates a multi-page PDF of diagnostic plots for directional (SPT sector) analyses.
#' Panels are selected from \code{spt.res$sector.results} and plotted using the vertices
#' and distances stored in \code{spt.res$sector.map}. For each selected (feature, sector)
#' pair, the function can produce:
#' \itemize{
#'   \item abundance vs distance diagnostics on the logit scale (among carriers),
#'   \item optional CLR abundance vs distance diagnostics (among carriers),
#'   \item presence/absence vs distance diagnostics with a binned prevalence curve,
#'   \item optional quantile-bin summary diagnostics if \code{spt.res$bin.results} is present.
#' }
#'
#' The distance axis can be either raw geodesic distance or rank-transformed distance
#' (\code{distance.transform="rank"}), which is often more robust when a small number of
#' vertices at very small distances induce leverage.
#'
#' @param pdf.file Character scalar. Path to the output PDF file.
#' @param spt.res List. Result returned by \code{test.disk.feature.distance.association.spt()},
#'   containing at least \code{$sector.results} and \code{$sector.map}. If quantile-bin outputs
#'   are present, \code{$bin.results} may also be used.
#' @param X Numeric matrix (or coercible data.frame) with vertices in rows and features in columns.
#' @param eps Non-negative numeric scalar. Carrier threshold; points are treated as carriers if
#'   \code{X[vertex, feature] > eps}. Default is 0.
#' @param pseudo.count Positive numeric scalar. Pseudocount used for the logit transform and CLR.
#'   Default is \code{1e-6}.
#' @param do.clr Logical. If TRUE (default), compute CLR within each plotted sector and include CLR
#'   abundance diagnostics.
#' @param spline.df Optional numeric. Degrees of freedom passed to \code{\link[stats]{smooth.spline}}.
#'   If NULL (default), \code{smooth.spline} selects smoothing internally.
#' @param distance.transform Character string. One of \code{"raw"} or \code{"rank"}.
#'   Passed to \code{plot.abundance.distance.diagnostics()} and \code{plot.presence.distance.diagnostics()}.
#' @param xlab Optional character scalar. X-axis label used for raw/rank diagnostic plots. If NULL
#'   (default), a label is chosen based on \code{distance.transform} and the sector id.
#'
#' @param alpha Numeric scalar in \code{(0,1)}. Threshold used to select (feature, sector) pairs
#'   from \code{spt.res$sector.results}. Default is 0.05. If no rows satisfy the selection rule,
#'   the function falls back to the top-ranked rows by available p-values.
#' @param select.by Character string controlling how panels are selected. One of:
#'   \code{"any"}, \code{"clr.lm.p.adj"}, \code{"clr.spearman.p.adj"}, \code{"lm.p.adj"},
#'   \code{"spearman.p.adj"}, \code{"presence.p.adj"}.
#'   If \code{"any"} (default), a row is selected if any available adjusted p-value family is
#'   \code{<= alpha}.
#' @param order.by Character string controlling the ordering of panels in the PDF. One of
#'   \code{"sector"} (default), \code{"feature"}, or \code{"p"}.
#' @param max.panels Integer or \code{Inf}. Maximum number of (feature, sector) panels to plot.
#'   Default is \code{Inf}.
#'
#' @param include.bin.plots Logical. If TRUE (default) and \code{spt.res$bin.results} is present,
#'   include quantile-bin diagnostic plots via \code{plot.quantile.bin.diagnostics()}.
#' @param include.bin.logit Logical. If TRUE (default), include logit-bin abundance plots when available.
#' @param include.bin.clr Logical. If TRUE (default), include CLR-bin abundance and presence plots when available.
#'
#' @param include.logit Logical. If TRUE (default), include logit-scale abundance diagnostics.
#' @param include.presence Logical. If TRUE (default), include presence/absence diagnostics.
#'
#' @return Invisibly returns TRUE. Called for its side effects (PDF creation).
#'
#' @examples
#' \dontrun{
#' ## After running an SPT analysis:
#' ## spt.res <- test.disk.feature.distance.association.spt(...)
#'
#' make.spt.diagnostics.pdf(
#'   pdf.file = "M1_spt_rank_diagnostics.pdf",
#'   spt.res = spt.res,
#'   X = phi.zmb,
#'   eps = 0,
#'   distance.transform = "rank",
#'   alpha = 0.05,
#'   select.by = "any",
#'   order.by = "sector"
#' )
#' }
#'
#' @seealso
#' \code{\link{test.disk.feature.distance.association.spt}},
#' \code{\link{plot.abundance.distance.diagnostics}},
#' \code{\link{plot.presence.distance.diagnostics}},
#' \code{\link{plot.quantile.bin.diagnostics}}
#'
#' @export
make.spt.diagnostics.pdf <- function(pdf.file,
                                     spt.res,
                                     X,
                                     eps = 0,
                                     pseudo.count = 1e-6,
                                     do.clr = TRUE,
                                     spline.df = NULL,
                                     distance.transform = c("raw", "rank"),
                                     xlab = NULL,
                                     ## Selection of (feature, sector) panels to plot
                                     alpha = 0.05,
                                     select.by = c("any", "clr.lm.p.adj", "clr.spearman.p.adj",
                                                   "lm.p.adj", "spearman.p.adj", "presence.p.adj"),
                                     order.by = c("sector", "feature", "p"),
                                     max.panels = Inf,
                                     ## Include quantile-bin diagnostics if present
                                     include.bin.plots = TRUE,
                                     include.bin.logit = TRUE,
                                     include.bin.clr = TRUE,
                                     ## Plot toggles
                                     include.logit = TRUE,
                                     include.presence = TRUE) {

    distance.transform <- match.arg(distance.transform)
    select.by <- match.arg(select.by)
    order.by <- match.arg(order.by)

    if (!is.character(pdf.file) || length(pdf.file) != 1L || is.na(pdf.file)) {
        stop("pdf.file must be a single character path.")
    }
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix (or coercible data.frame).")
    if (is.null(colnames(X))) stop("X must have column names (feature names).")

    if (!is.list(spt.res) || is.null(spt.res$sector.results) || is.null(spt.res$sector.map)) {
        stop("spt.res must be an SPT result list with components $sector.results and $sector.map.")
    }

    sector.df <- spt.res$sector.results
    sector.map <- spt.res$sector.map
    bin.results <- spt.res$bin.results

    if (!is.data.frame(sector.df) || nrow(sector.df) == 0L) stop("spt.res$sector.results is empty.")
    if (!is.data.frame(sector.map) || nrow(sector.map) == 0L) stop("spt.res$sector.map is empty.")
    if (!all(c("vertex", "dist", "sector") %in% names(sector.map))) {
        stop("spt.res$sector.map must contain columns: vertex, dist, sector.")
    }

    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) stop("eps must be a single non-negative numeric.")
    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric.")
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
        stop("alpha must be a single numeric in (0, 1).")
    }
    if (!is.numeric(max.panels) || length(max.panels) != 1L || is.na(max.panels) || max.panels < 1) {
        stop("max.panels must be a single integer >= 1 (or Inf).")
    }

    ## ---- helpers ----
    compute.clr.disk <- function(X.sub, pseudo.count = 1e-6) {
        log.X <- log(X.sub + pseudo.count)
        row.center <- rowMeans(log.X)
        sweep(log.X, 1, row.center, FUN = "-")
    }

    get.sector.vertices.dists <- function(sector.map, sector) {
        idx <- which(sector.map$sector == sector)
        v <- as.integer(sector.map$vertex[idx])
        d <- as.numeric(sector.map$dist[idx])
        ok <- is.finite(d) & !is.na(v)
        v <- v[ok]
        d <- d[ok]
        names(d) <- as.character(v)
        list(vertices = v, dists = d)
    }

    ## ---- choose panels to plot ----
    ## We use adjusted p-values when available; otherwise fall back to raw p-values.
    has.col <- function(nm) nm %in% names(sector.df)

    get.p <- function(nm) {
        if (has.col(nm)) return(sector.df[[nm]])
        rep(NA_real_, nrow(sector.df))
    }

    p.pres <- get.p("presence.p.adj")
    p.spear <- if (has.col("spearman.p.adj")) get.p("spearman.p.adj") else get.p("spearman.p")
    p.lm <- if (has.col("lm.p.adj")) get.p("lm.p.adj") else get.p("lm.p")
    p.clr.spear <- if (has.col("clr.spearman.p.adj")) get.p("clr.spearman.p.adj") else get.p("clr.spearman.p")
    p.clr.lm <- if (has.col("clr.lm.p.adj")) get.p("clr.lm.p.adj") else get.p("clr.lm.p")

    keep <- rep(FALSE, nrow(sector.df))

    if (select.by == "any") {
        ## Keep if any available p-value family is <= alpha
        keep <- keep | (is.finite(p.pres) & p.pres <= alpha)
        keep <- keep | (is.finite(p.spear) & p.spear <= alpha)
        keep <- keep | (is.finite(p.lm) & p.lm <= alpha)
        keep <- keep | (is.finite(p.clr.spear) & p.clr.spear <= alpha)
        keep <- keep | (is.finite(p.clr.lm) & p.clr.lm <= alpha)
    } else {
        p.sel <- get.p(select.by)
        keep <- is.finite(p.sel) & (p.sel <= alpha)
    }

    sub.df <- sector.df[keep, , drop = FALSE]

    ## If nothing passes alpha, fall back to the top-ranked rows by clr.lm.p.adj (if present)
    if (nrow(sub.df) == 0L) {
        if (has.col("clr.lm.p.adj") && any(is.finite(sector.df$clr.lm.p.adj))) {
            o <- order(sector.df$clr.lm.p.adj)
        } else if (has.col("clr.lm.p") && any(is.finite(sector.df$clr.lm.p))) {
            o <- order(sector.df$clr.lm.p)
        } else if (has.col("lm.p.adj") && any(is.finite(sector.df$lm.p.adj))) {
            o <- order(sector.df$lm.p.adj)
        } else {
            o <- seq_len(nrow(sector.df))
        }
        sub.df <- sector.df[o, , drop = FALSE]
    }

    ## Ordering
    if (order.by == "sector" && has.col("sector")) {
        ## Within sector order by clr.lm.p.adj if present, else by lm.p.adj / lm.p
        if (has.col("clr.lm.p.adj")) {
            sub.df <- sub.df[order(sub.df$sector, sub.df$clr.lm.p.adj), , drop = FALSE]
        } else if (has.col("lm.p.adj")) {
            sub.df <- sub.df[order(sub.df$sector, sub.df$lm.p.adj), , drop = FALSE]
        } else if (has.col("lm.p")) {
            sub.df <- sub.df[order(sub.df$sector, sub.df$lm.p), , drop = FALSE]
        } else {
            sub.df <- sub.df[order(sub.df$sector), , drop = FALSE]
        }
    } else if (order.by == "feature" && has.col("feature")) {
        if (has.col("clr.lm.p.adj")) {
            sub.df <- sub.df[order(sub.df$feature, sub.df$clr.lm.p.adj), , drop = FALSE]
        } else {
            sub.df <- sub.df[order(sub.df$feature), , drop = FALSE]
        }
    } else if (order.by == "p") {
        if (has.col("clr.lm.p.adj")) {
            sub.df <- sub.df[order(sub.df$clr.lm.p.adj), , drop = FALSE]
        } else if (has.col("lm.p.adj")) {
            sub.df <- sub.df[order(sub.df$lm.p.adj), , drop = FALSE]
        } else if (has.col("lm.p")) {
            sub.df <- sub.df[order(sub.df$lm.p), , drop = FALSE]
        }
    }

    ## Cap number of panels
    if (is.finite(max.panels) && nrow(sub.df) > max.panels) {
        sub.df <- sub.df[seq_len(as.integer(max.panels)), , drop = FALSE]
    }

    ## ---- produce PDF ----
    grDevices::pdf(pdf.file, width = 8.5, height = 6.5)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (i in seq_len(nrow(sub.df))) {

        f <- as.character(sub.df$feature[i])
        s <- as.integer(sub.df$sector[i])

        if (!f %in% colnames(X)) next

        sd <- get.sector.vertices.dists(sector.map, s)
        v.s <- sd$vertices
        d.s <- sd$dists

        if (length(v.s) < 2L) next

        ## Compute sector CLR if requested
        clr.s <- NULL
        if (isTRUE(do.clr)) {
            clr.s <- compute.clr.disk(X[v.s, , drop = FALSE], pseudo.count = pseudo.count)
            ## clr.s rows aligned to v.s order
        }

        ## Choose an informative x-label if not provided
        xlab.use <- xlab
        if (is.null(xlab.use)) {
            if (distance.transform == "rank") {
                xlab.use <- paste0("rank(distance) within sector ", s)
            } else {
                xlab.use <- paste0("geodesic distance within sector ", s)
            }
        }

        ## --- (1) logit abundance scatter/hex diagnostics ---
        if (isTRUE(include.logit)) {
            plot.abundance.distance.diagnostics(
                feature = f,
                X = X,
                vertices = v.s,
                dists = d.s,
                eps = eps,
                pseudo.count = pseudo.count,
                use.clr = FALSE,
                main.prefix = paste0("SPT sector ", s, " | logit: "),
                spline.df = spline.df,
                distance.transform = distance.transform,
                xlab = xlab.use
            )
        }

        ## --- (2) CLR abundance diagnostics ---
        if (isTRUE(do.clr) && !is.null(clr.s)) {
            plot.abundance.distance.diagnostics(
                feature = f,
                X = X,
                vertices = v.s,
                dists = d.s,
                eps = eps,
                pseudo.count = pseudo.count,
                use.clr = TRUE,
                clr.mat = clr.s,
                main.prefix = paste0("SPT sector ", s, " | CLR: "),
                spline.df = spline.df,
                distance.transform = distance.transform,
                xlab = xlab.use
            )
        }

        ## --- (3) presence diagnostics ---
        if (isTRUE(include.presence)) {
            plot.presence.distance.diagnostics(
                feature = f,
                X = X,
                vertices = v.s,
                dists = d.s,
                eps = eps,
                n.bins = 10L,
                main.prefix = paste0("SPT sector ", s, " | presence: "),
                distance.transform = distance.transform,
                xlab = xlab.use
            )
        }

        ## --- (4) quantile-bin diagnostics if available ---
        if (isTRUE(include.bin.plots) && !is.null(bin.results) && length(bin.results) > 0L) {

            key <- paste0(f, "||", s)
            br <- bin.results[[key]]
            if (!is.null(br) && is.list(br)) {

                if (isTRUE(include.bin.clr) && !is.null(br$clr)) {

                    plot.quantile.bin.diagnostics(
                        bin.res = br$clr,
                        what = "abundance",
                        main = paste0("SPT sector ", s, " | qbin CLR abundance: ", f),
                        xlab = "distance (bin midpoint)",
                        ylab = "CLR (bin median)"
                    )

                    if (!is.null(br$clr$bins) &&
                        is.data.frame(br$clr$bins) &&
                        all(c("d.mid.all", "presence.p") %in% names(br$clr$bins))) {

                        plot.quantile.bin.diagnostics(
                            bin.res = br$clr,
                            what = "presence",
                            main = paste0("SPT sector ", s, " | qbin presence: ", f),
                            xlab = "distance (bin midpoint)",
                            ylab = "prevalence"
                        )
                    }
                }

                if (isTRUE(include.bin.logit) && !is.null(br$logit)) {
                    plot.quantile.bin.diagnostics(
                        bin.res = br$logit,
                        what = "abundance",
                        main = paste0("SPT sector ", s, " | qbin logit abundance: ", f),
                        xlab = "distance (bin midpoint)",
                        ylab = "logit (bin median)"
                    )
                }
            }
        }
    }

    invisible(TRUE)
}
