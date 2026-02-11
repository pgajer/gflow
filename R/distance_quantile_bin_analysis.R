#' Distance-Quantile Bin Analysis for Feature–Distance Signal
#'
#' @description
#' Perform a distance-quantile (equal-count) bin analysis of the relationship between
#' a feature's abundance (among carriers) and geodesic distance, and optionally
#' carriage (presence/absence) versus distance.
#'
#' This analysis is useful when models suffer from leverage at small distances:
#' instead of modeling abundance as a smooth function of raw distance, distance
#' is discretized into quantile bins and summaries/trends are computed per bin.
#'
#' The function supports:
#' \itemize{
#'   \item Abundance among carriers: bin-wise mean/median with standard error, plus a simple trend test.
#'   \item Presence/absence: bin-wise prevalence with standard error, plus a simple trend test.
#' }
#'
#' For abundance, the trend test is computed on bin midpoints using either:
#' \itemize{
#'   \item Spearman correlation between bin midpoint and bin summary (default), or
#'   \item Weighted linear regression on bin midpoint using bin counts as weights.
#' }
#'
#' @param x Numeric vector of feature values (e.g., relative abundances) aligned to \code{d}.
#' @param d Numeric vector of distances aligned to \code{x}.
#' @param eps Non-negative threshold defining a carrier: \code{x > eps}. Default 0.
#' @param n.bins Integer number of quantile bins. Default 10.
#' @param min.per.bin Minimum number of observations required in a bin to keep it. Default 5.
#' @param carriers.only Logical; if TRUE (default), abundance summaries are computed only among carriers.
#'   If FALSE, summaries are computed on all points (including zeros), which may be preferable for CLR.
#' @param abundance.summary Character string: \code{"median"} (default) or \code{"mean"}.
#' @param trend.test Character string: \code{"spearman"} (default) or \code{"wls"}.
#' @param with.presence Logical; if TRUE (default), also compute bin-wise prevalence and its trend test.
#' @param na.rm Logical; remove missing/non-finite values in \code{x} or \code{d}. Default TRUE.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{bins}: a data.frame of bin-wise summaries (distance range, midpoint, counts,
#'     abundance summary, abundance SE; and if requested, prevalence and prevalence SE).
#'   \item \code{tests}: a list of simple trend tests for abundance and (optionally) presence.
#'   \item \code{bin.assignment}: integer vector of bin labels (aligned to the filtered input).
#'   \item \code{keep}: logical vector indicating which original indices were used after filtering.
#' }
#'
#' @examples
#' \dontrun{
#' ## Example for one phylotype using raw relative abundance
#' x <- phi.zmb[M1.disk.res$vertices, "Prevotella_amnii"]
#' d <- as.numeric(M1.disk.res$dists[as.character(M1.disk.res$vertices)])
#'
#' res <- distance.quantile.bin.analysis(x, d, eps = 0, n.bins = 10)
#' res$bins
#' res$tests
#'
#' ## For CLR values (already transformed), set carriers.only = FALSE
#' x.clr <- clr.mat[, "Prevotella_amnii"]
#' res.clr <- distance.quantile.bin.analysis(x.clr, d, carriers.only = FALSE)
#' }
#'
#' @export
distance.quantile.bin.analysis <- function(x,
                                          d,
                                          eps = 0,
                                          n.bins = 10L,
                                          min.per.bin = 5L,
                                          carriers.only = TRUE,
                                          abundance.summary = c("median", "mean"),
                                          trend.test = c("spearman", "wls"),
                                          with.presence = TRUE,
                                          na.rm = TRUE) {

    abundance.summary <- match.arg(abundance.summary)
    trend.test <- match.arg(trend.test)

    if (!is.numeric(x) || !is.numeric(d)) {
        stop("x and d must be numeric vectors.")
    }
    if (length(x) != length(d)) {
        stop("x and d must have the same length.")
    }
    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) {
        stop("eps must be a single non-negative numeric value.")
    }
    if (!is.numeric(n.bins) || length(n.bins) != 1L || is.na(n.bins) || n.bins < 2L) {
        stop("n.bins must be an integer >= 2.")
    }
    n.bins <- as.integer(n.bins)

    if (!is.numeric(min.per.bin) || length(min.per.bin) != 1L || is.na(min.per.bin) || min.per.bin < 1L) {
        stop("min.per.bin must be an integer >= 1.")
    }
    min.per.bin <- as.integer(min.per.bin)

    ## ---- filter to finite ----
    if (isTRUE(na.rm)) {
        keep <- is.finite(x) & is.finite(d)
        x <- x[keep]
        d <- d[keep]
    } else {
        keep <- rep(TRUE, length(x))
        if (any(!is.finite(x)) || any(!is.finite(d))) {
            stop("Non-finite values detected in x or d; set na.rm=TRUE to drop them.")
        }
    }

    if (length(x) < 2L) {
        stop("Too few observations after filtering.")
    }

    ## ---- define carriers ----
    carrier <- (x > eps)

    ## ---- choose values used for abundance summaries ----
    if (isTRUE(carriers.only)) {
        idx.ab <- carrier
    } else {
        idx.ab <- rep(TRUE, length(x))
    }

    if (sum(idx.ab) < 2L) {
        stop("Too few observations for abundance summaries (after carriers.only filtering).")
    }

    d.ab <- d[idx.ab]
    x.ab <- x[idx.ab]

    ## ---- create equal-count (quantile) bins on distance (using abundance subset) ----
    probs <- seq(0, 1, length.out = n.bins + 1L)
    brks <- unique(stats::quantile(d.ab, probs = probs, na.rm = TRUE, type = 7))

    ## If too many ties in d, quantiles can collapse
    if (length(brks) < 3L) {
        stop("Distance quantile breaks collapsed (too many ties in d). Consider fewer bins or jitter d slightly.")
    }

    bin.ab <- cut(d.ab, breaks = brks, include.lowest = TRUE, right = TRUE)
    bin.levels <- levels(bin.ab)

    ## ---- per-bin summaries for abundance ----
    sum.ab <- function(z) {
        if (abundance.summary == "median") stats::median(z) else mean(z)
    }

    ## midpoint as mean of range in that bin
    mid.ab <- tapply(d.ab, bin.ab, function(z) mean(range(z)))
    n.ab <- tapply(x.ab, bin.ab, function(z) sum(is.finite(z)))
    val.ab <- tapply(x.ab, bin.ab, function(z) sum.ab(z))
    se.ab <- tapply(x.ab, bin.ab, function(z) {
        z <- z[is.finite(z)]
        if (length(z) <= 1L) return(NA_real_)
        stats::sd(z) / sqrt(length(z))
    })

    ## ---- presence/absence summaries (bins defined on all distances for comparability) ----
    ## Use same breaks brks but apply to all points (not just abundance subset)
    bin.all <- cut(d, breaks = brks, include.lowest = TRUE, right = TRUE)

    n.all <- tapply(d, bin.all, length)
    prev <- tapply(carrier, bin.all, mean)
    prev.se <- tapply(carrier, bin.all, function(z) {
        ## binomial SE
        p <- mean(z)
        n <- length(z)
        if (n <= 1L) return(NA_real_)
        sqrt(p * (1 - p) / n)
    })
    mid.all <- tapply(d, bin.all, function(z) mean(range(z)))

    ## ---- assemble bins table ----
    bins.df <- data.frame(
        bin = bin.levels,
        d.min = as.numeric(tapply(d.ab, bin.ab, min)),
        d.max = as.numeric(tapply(d.ab, bin.ab, max)),
        d.mid = as.numeric(mid.ab),
        n.abundance = as.integer(n.ab),
        abundance = as.numeric(val.ab),
        abundance.se = as.numeric(se.ab),
        stringsAsFactors = FALSE
    )

    if (isTRUE(with.presence)) {
        ## Align presence bins to the same order as bins.df$bin
        bins.df$n.total <- as.integer(n.all[bins.df$bin])
        bins.df$presence.p <- as.numeric(prev[bins.df$bin])
        bins.df$presence.se <- as.numeric(prev.se[bins.df$bin])
        bins.df$d.mid.all <- as.numeric(mid.all[bins.df$bin])
    }

    ## ---- drop bins with too few points ----
    keep.bin <- bins.df$n.abundance >= min.per.bin
    if (isTRUE(with.presence)) {
        keep.bin <- keep.bin & (bins.df$n.total >= min.per.bin)
    }
    bins.df <- bins.df[keep.bin, , drop = FALSE]

    if (nrow(bins.df) < 3L) {
        stop("Too few bins remain after min.per.bin filtering; reduce n.bins or min.per.bin.")
    }

    ## ---- trend tests ----
    tests <- list()

    ## abundance trend (bin midpoints vs abundance summary)
    if (trend.test == "spearman") {
        ct <- suppressWarnings(stats::cor.test(bins.df$d.mid, bins.df$abundance, method = "spearman", exact = FALSE))
        tests$abundance <- list(
            method = "spearman (bin mid vs bin summary)",
            estimate = unname(ct$estimate),
            p.value = ct$p.value
        )
    } else {
        ## WLS on bin midpoints; weights = number of observations in bin
        fit <- stats::lm(abundance ~ d.mid, data = bins.df, weights = n.abundance)
        sm <- summary(fit)$coefficients["d.mid", ]
        tests$abundance <- list(
            method = "weighted least squares (bin mid vs bin summary)",
            slope = unname(sm["Estimate"]),
            p.value = unname(sm["Pr(>|t|)"])
        )
    }

    ## presence trend (bin midpoints vs prevalence)
    if (isTRUE(with.presence)) {
        if (trend.test == "spearman") {
            ct <- suppressWarnings(stats::cor.test(bins.df$d.mid.all, bins.df$presence.p, method = "spearman", exact = FALSE))
            tests$presence <- list(
                method = "spearman (bin mid vs bin prevalence)",
                estimate = unname(ct$estimate),
                p.value = ct$p.value
            )
        } else {
            fit <- stats::lm(presence.p ~ d.mid.all, data = bins.df, weights = n.total)
            sm <- summary(fit)$coefficients["d.mid.all", ]
            tests$presence <- list(
                method = "weighted least squares (bin mid vs bin prevalence)",
                slope = unname(sm["Estimate"]),
                p.value = unname(sm["Pr(>|t|)"])
            )
        }
    }

    ## ---- return ----
    ## bin.assignment returned for the abundance subset (aligned to filtered x/d subset)
    out <- list(
        bins = bins.df,
        tests = tests,
        bin.assignment = bin.ab,
        keep = keep
    )

    out
}

#' @export
plot.distance.quantile.bins <- function(x,
                                       what = c("abundance", "presence"),
                                       main = NULL,
                                       ...) {
    res <- x

    what <- match.arg(what)
    df <- res$bins

    if (what == "abundance") {
        graphics::plot(df$d.mid, df$abundance,
                       xlab = "distance (bin midpoint)",
                       ylab = "abundance (bin summary)",
                       main = if (is.null(main)) "Abundance vs distance bins" else main,
                       pch = 16)
        graphics::segments(df$d.mid, df$abundance - df$abundance.se,
                           df$d.mid, df$abundance + df$abundance.se)
    } else {
        if (!("presence.p" %in% names(df))) stop("Presence summaries not present in res$bins.")
        graphics::plot(df$d.mid.all, df$presence.p,
                       xlab = "distance (bin midpoint)",
                       ylab = "prevalence",
                       main = if (is.null(main)) "Presence vs distance bins" else main,
                       pch = 16, ylim = c(0, 1))
        graphics::segments(df$d.mid.all, df$presence.p - df$presence.se,
                           df$d.mid.all, df$presence.p + df$presence.se)
    }

    invisible(TRUE)
}
