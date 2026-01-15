#' Test Feature–Distance Associations Within a Geodesic Disk
#'
#' @description
#' For features with prevalence at least \code{p.min} inside a supplied vertex set,
#' test association between (transformed) relative abundance and geodesic distance.
#' Performs a two-part assessment:
#' (1) presence/absence vs distance (logistic regression),
#' (2) abundance among carriers vs distance (Spearman, LM, robust LM).
#'
#' Adds a compositional sensitivity analysis using CLR (centered log-ratio)
#' computed within the disk with a pseudocount.
#'
#' @param X Numeric matrix: rows are vertices/samples, columns are features.
#' @param vertices Integer vector of vertex indices (1-based) defining the disk/subset.
#' @param dists Numeric vector of distances for vertices in \code{vertices}; either
#'   named by vertex index (character) or aligned with \code{vertices}.
#' @param eps Non-negative threshold defining a carrier: \code{X > eps}. Default 0.
#' @param p.min Minimum prevalence within \code{vertices} to include a feature. Default 0.5.
#' @param transform One of \code{"log"}, \code{"logit"}, \code{"arcsin_sqrt"} for the
#'   primary (non-CLR) carrier-only analysis.
#' @param pseudo.count Small positive value to stabilize transforms at 0/1 boundaries,
#'   and used for CLR. Default 1e-6.
#' @param min.carriers Minimum number of carriers required to run carrier-only tests.
#'   Default 20.
#' @param do.clr Logical; if TRUE (default), run CLR sensitivity analysis.
#' @param clr.over Character; one of \code{"carriers"} (default) or \code{"all"}.
#'   Controls whether CLR association is computed over carriers only (matching the
#'   primary analysis) or over all disk vertices.
#' @param do.bootstrap.rlm Logical; if TRUE, bootstrap RLM slope p-values. Default FALSE.
#' @param n.boot Number of bootstrap resamples for RLM slope p-values. Default 2000.
#'
#' @return A data.frame with one row per tested feature, including prevalence, sample sizes,
#'   and p-values (raw and BH-adjusted) for presence and carrier-only tests, plus CLR
#'   sensitivity columns when enabled.
#'
#' @export
test.disk.feature.distance.association <- function(X,
                                                   vertices,
                                                   dists,
                                                   eps = 0,
                                                   p.min = 0.5,
                                                   transform = c("log", "logit", "arcsin_sqrt"),
                                                   pseudo.count = 1e-6,
                                                   min.carriers = 20L,
                                                   do.clr = TRUE,
                                                   clr.over = c("carriers", "all"),
                                                   do.bootstrap.rlm = FALSE,
                                                   n.boot = 2000L) {

    transform <- match.arg(transform)
    clr.over <- match.arg(clr.over)

    ## ---- validate ----
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix (or coercible data.frame).")
    if (is.null(colnames(X))) stop("X must have column names.")
    if (!is.numeric(vertices) || length(vertices) == 0L) stop("vertices must be a non-empty numeric/integer vector.")

    vertices <- as.integer(vertices)
    vertices <- vertices[!is.na(vertices)]
    vertices <- vertices[vertices >= 1L & vertices <= nrow(X)]
    if (length(vertices) == 0L) stop("No valid vertices after filtering to [1, nrow(X)].")

    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) stop("eps must be a single non-negative numeric.")
    if (!is.numeric(p.min) || length(p.min) != 1L || is.na(p.min) || p.min < 0 || p.min > 1) stop("p.min must be in [0, 1].")
    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric value.")
    }

    ## ---- align distances to vertices ----
    d <- NULL
    if (!is.null(names(dists))) {
        ## treat names as vertex ids
        d <- as.numeric(dists[as.character(vertices)])
    } else {
        ## assume aligned with vertices order
        if (length(dists) != length(vertices)) stop("If dists is unnamed, it must be aligned with vertices (same length).")
        d <- as.numeric(dists)
    }
    ok <- is.finite(d)
    vertices <- vertices[ok]
    d <- d[ok]
    if (length(vertices) < 2L) stop("Too few vertices with finite distances.")

    ## ---- prevalence filter within disk ----
    prev.df <- vertices.feature.carriers(vertices = vertices, X = X, eps = eps, match.by = "index")
    keep <- prev.df$p >= p.min
    feats <- rownames(prev.df)[keep]
    if (length(feats) == 0L) {
        out <- data.frame()
        return(out)
    }

    ## ---- transformation helper (primary analysis) ----
    transform.vec <- function(x) {
        if (transform == "log") {
            return(log(x + pseudo.count))
        }
        if (transform == "logit") {
            ## clamp away from 0/1 using pseudo.count
            x2 <- (x + pseudo.count) / (1 + 2 * pseudo.count)
            return(qlogis(x2))
        }
        ## arcsin-sqrt (expects proportions in [0,1])
        x2 <- pmin(pmax(x, 0), 1)
        return(asin(sqrt(x2)))
    }

    ## ---- CLR matrix for disk (sensitivity) ----
    clr.mat <- NULL
    if (isTRUE(do.clr)) {
        ## Compute CLR within the disk over all features in X
        ## clr(x_i) = log(x_i + pc) - mean_j log(x_j + pc), per sample (row)
        X.disk <- X[vertices, , drop = FALSE]
        log.X <- log(X.disk + pseudo.count)
        row.center <- rowMeans(log.X)
        clr.mat <- sweep(log.X, 1, row.center, FUN = "-")
        ## rownames(clr.mat) correspond to disk vertices order
    }

    ## ---- per-feature tests ----
    res.list <- vector("list", length(feats))
    names(res.list) <- feats

    for (j in seq_along(feats)) {
        f <- feats[j]
        x <- X[vertices, f]
        carrier <- x > eps

        n.all <- length(x)
        n.car <- sum(carrier)

        ## presence model: carrier ~ distance
        pres.p <- NA_real_
        pres.slope <- NA_real_
        if (n.car > 0L && n.car < n.all) {
            fit.pres <- suppressWarnings(stats::glm(carrier ~ d, family = stats::binomial()))
            sm <- summary(fit.pres)
            pres.slope <- unname(stats::coef(fit.pres)["d"])
            pres.p <- sm$coefficients["d", "Pr(>|z|)"]
        }

        ## carrier-only tests (primary transform)
        rho <- NA_real_
        rho.p <- NA_real_
        lm.slope <- NA_real_
        lm.p <- NA_real_
        rlm.slope <- NA_real_
        rlm.p <- NA_real_

        if (n.car >= as.integer(min.carriers)) {
            y <- transform.vec(x[carrier])
            dd <- d[carrier]

            ## Spearman
            sp <- suppressWarnings(stats::cor.test(y, dd, method = "spearman", exact = FALSE))
            rho <- unname(sp$estimate)
            rho.p <- sp$p.value

            ## LM
            fit.lm <- stats::lm(y ~ dd)
            lm.slope <- unname(stats::coef(fit.lm)["dd"])
            lm.p <- summary(fit.lm)$coefficients["dd", "Pr(>|t|)"]

            ## RLM
            if (requireNamespace("MASS", quietly = TRUE)) {
                fit.rlm <- MASS::rlm(y ~ dd)
                rlm.slope <- unname(stats::coef(fit.rlm)["dd"])

                if (isTRUE(do.bootstrap.rlm)) {
                    set.seed(1)
                    b.slopes <- numeric(n.boot)
                    n.y <- length(y)
                    for (b in seq_len(n.boot)) {
                        ii <- sample.int(n.y, size = n.y, replace = TRUE)
                        fb <- MASS::rlm(y[ii] ~ dd[ii])
                        b.slopes[b] <- unname(stats::coef(fb)[2])
                    }
                    rlm.p <- 2 * min(mean(b.slopes <= 0), mean(b.slopes >= 0))
                }
            }
        }

        ## ---- CLR sensitivity tests ----
        clr.rho <- NA_real_
        clr.rho.p <- NA_real_
        clr.lm.slope <- NA_real_
        clr.lm.p <- NA_real_
        clr.rlm.slope <- NA_real_
        clr.rlm.p <- NA_real_

        if (isTRUE(do.clr)) {
            if (clr.over == "all") {
                idx.use <- rep(TRUE, length(vertices))
            } else {
                idx.use <- carrier
            }

            n.use <- sum(idx.use)
            if (n.use >= 2L) {
                y.clr <- clr.mat[idx.use, f]
                d.clr <- d[idx.use]

                ## Spearman on CLR
                spc <- suppressWarnings(stats::cor.test(y.clr, d.clr, method = "spearman", exact = FALSE))
                clr.rho <- unname(spc$estimate)
                clr.rho.p <- spc$p.value

                ## LM on CLR
                fit.clr.lm <- stats::lm(y.clr ~ d.clr)
                clr.lm.slope <- unname(stats::coef(fit.clr.lm)["d.clr"])
                clr.lm.p <- summary(fit.clr.lm)$coefficients["d.clr", "Pr(>|t|)"]

                ## RLM on CLR
                if (requireNamespace("MASS", quietly = TRUE)) {
                    fit.clr.rlm <- MASS::rlm(y.clr ~ d.clr)
                    clr.rlm.slope <- unname(stats::coef(fit.clr.rlm)["d.clr"])

                    if (isTRUE(do.bootstrap.rlm)) {
                        set.seed(1)
                        b.slopes <- numeric(n.boot)
                        n.y <- length(y.clr)
                        for (b in seq_len(n.boot)) {
                            ii <- sample.int(n.y, size = n.y, replace = TRUE)
                            fb <- MASS::rlm(y.clr[ii] ~ d.clr[ii])
                            b.slopes[b] <- unname(stats::coef(fb)[2])
                        }
                        clr.rlm.p <- 2 * min(mean(b.slopes <= 0), mean(b.slopes >= 0))
                    }
                }
            }
        }

        res.list[[j]] <- data.frame(
            feature = f,
            p.disk = prev.df[f, "p"],
            n.disk = n.all,
            n.carriers = n.car,
            presence.slope = pres.slope,
            presence.p = pres.p,
            spearman.rho = rho,
            spearman.p = rho.p,
            lm.slope = lm.slope,
            lm.p = lm.p,
            rlm.slope = rlm.slope,
            rlm.p = rlm.p,
            clr.over = if (isTRUE(do.clr)) clr.over else NA_character_,
            clr.spearman.rho = clr.rho,
            clr.spearman.p = clr.rho.p,
            clr.lm.slope = clr.lm.slope,
            clr.lm.p = clr.lm.p,
            clr.rlm.slope = clr.rlm.slope,
            clr.rlm.p = clr.rlm.p,
            stringsAsFactors = FALSE
        )
    }

    out <- do.call(rbind, res.list)
    rownames(out) <- out$feature

    ## ---- FDR adjust within families of tests ----
    out$presence.p.adj <- stats::p.adjust(out$presence.p, method = "BH")
    out$spearman.p.adj <- stats::p.adjust(out$spearman.p, method = "BH")
    out$lm.p.adj <- stats::p.adjust(out$lm.p, method = "BH")
    if (!all(is.na(out$rlm.p))) {
        out$rlm.p.adj <- stats::p.adjust(out$rlm.p, method = "BH")
    }

    if (isTRUE(do.clr)) {
        out$clr.spearman.p.adj <- stats::p.adjust(out$clr.spearman.p, method = "BH")
        out$clr.lm.p.adj <- stats::p.adjust(out$clr.lm.p, method = "BH")
        if (!all(is.na(out$clr.rlm.p))) {
            out$clr.rlm.p.adj <- stats::p.adjust(out$clr.rlm.p, method = "BH")
        }
    }

    out
}
