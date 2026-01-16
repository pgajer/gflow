#' Test Feature–Distance Associations Within a Geodesic Disk
#'
#' @description
#' For features with prevalence at least \code{p.min} inside a supplied vertex set,
#' test association between (transformed) relative abundance and geodesic distance.
#'
#' Supports three distance modes:
#' \itemize{
#'   \item \code{"raw"}: use raw distance values.
#'   \item \code{"rank"}: use scaled rank distance (robust to spacing/leverage in predictor).
#'   \item \code{"quantile.bin"}: perform distance-quantile bin analysis (equal-count bins)
#'         and compute bin-level trend tests on bin summaries.
#' }
#'
#' For \code{"raw"} and \code{"rank"}, results are vertex-level (carriers-only for abundance).
#' For \code{"quantile.bin"}, results include additional bin-level summaries in
#' \code{$bin.results}.
#'
#' @param X Numeric matrix: rows are vertices/samples, columns are features.
#' @param vertices Integer vector of vertex indices (1-based) defining the disk/subset.
#' @param dists Numeric vector of distances for vertices in \code{vertices}; either
#'   named by vertex index (character) or aligned with \code{vertices}.
#' @param eps Non-negative threshold defining a carrier: \code{X > eps}. Default 0.
#' @param p.min Minimum prevalence within \code{vertices} to include a feature. Default 0.5.
#' @param transform One of \code{"log"}, \code{"logit"}, \code{"arcsin_sqrt"} for the
#'   primary (non-CLR) carrier-only analysis.
#' @param pseudo.count Small positive value used to stabilize transforms and for CLR. Default 1e-6.
#' @param min.carriers Minimum number of carriers required to run carrier-only tests. Default 20.
#' @param do.clr Logical; if TRUE (default), run CLR sensitivity analysis.
#' @param clr.over One of \code{"carriers"} (default) or \code{"all"}.
#' @param distance.transform One of \code{"raw"} (default), \code{"rank"}, \code{"quantile.bin"}.
#' @param n.bins Integer; number of distance quantile bins for \code{"quantile.bin"}. Default 10.
#' @param min.per.bin Integer; minimum count required per bin to keep it. Default 5.
#' @param bin.trend.test One of \code{"spearman"} (default) or \code{"wls"} for bin-level trend tests.
#' @param with.presence.bins Logical; if TRUE (default), compute binned prevalence trend as well.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{vertex.results}: data.frame of per-feature vertex-level tests (always returned; empty for quantile.bin if \code{min.carriers} not met, but typically present).
#'   \item \code{bin.results}: (only for \code{distance.transform="quantile.bin"}) a named list of per-feature bin tables and bin-level tests.
#'   \item \code{meta}: list with settings used (distance.transform, n.bins, etc.).
#' }
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
                                                  distance.transform = c("raw", "rank", "quantile.bin"),
                                                  n.bins = 10L,
                                                  min.per.bin = 5L,
                                                  bin.trend.test = c("spearman", "wls"),
                                                  with.presence.bins = TRUE) {

    transform <- match.arg(transform)
    clr.over <- match.arg(clr.over)
    distance.transform <- match.arg(distance.transform)
    bin.trend.test <- match.arg(bin.trend.test)

    ## ---- validate X ----
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix (or coercible data.frame).")
    if (is.null(colnames(X))) stop("X must have column names.")

    ## ---- validate vertices ----
    if (!is.numeric(vertices) || length(vertices) == 0L) stop("vertices must be a non-empty numeric/integer vector.")
    vertices <- as.integer(vertices)
    vertices <- vertices[!is.na(vertices)]
    vertices <- vertices[vertices >= 1L & vertices <= nrow(X)]
    if (length(vertices) == 0L) stop("No valid vertices after filtering to [1, nrow(X)].")

    ## ---- validate scalars ----
    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) stop("eps must be a single non-negative numeric.")
    if (!is.numeric(p.min) || length(p.min) != 1L || is.na(p.min) || p.min < 0 || p.min > 1) stop("p.min must be in [0, 1].")
    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric value.")
    }
    n.bins <- as.integer(n.bins)
    min.per.bin <- as.integer(min.per.bin)
    if (distance.transform == "quantile.bin") {
        if (n.bins < 2L) stop("n.bins must be >= 2 for quantile.bin.")
        if (min.per.bin < 1L) stop("min.per.bin must be >= 1.")
    }

    ## ---- align distances to vertices ----
    if (!is.null(names(dists))) {
        d <- as.numeric(dists[as.character(vertices)])
    } else {
        if (length(dists) != length(vertices)) stop("If dists is unnamed, it must be aligned with vertices (same length).")
        d <- as.numeric(dists)
    }
    ok <- is.finite(d)
    vertices <- vertices[ok]
    d <- d[ok]
    if (length(vertices) < 2L) stop("Too few vertices with finite distances.")

    ## ---- distance transforms ----
    distance.rank <- function(dd) {
        r <- rank(dd, ties.method = "average")
        r / (length(r) + 1)
    }
    d.used <- d
    if (distance.transform == "rank") {
        d.used <- distance.rank(d)
    }

    ## ---- prevalence filter within disk ----
    prev.df <- vertices.feature.carriers(vertices = vertices, X = X, eps = eps, match.by = "index")
    keep <- prev.df$p >= p.min
    feats <- rownames(prev.df)[keep]
    if (length(feats) == 0L) {
        return(list(vertex.results = data.frame(),
                    bin.results = NULL,
                    meta = list(distance.transform = distance.transform)))
    }

    ## ---- transformation helper (primary analysis) ----
    transform.vec <- function(x) {
        if (transform == "log") {
            return(log(x + pseudo.count))
        }
        if (transform == "logit") {
            x2 <- (x + pseudo.count) / (1 + 2 * pseudo.count)
            return(qlogis(x2))
        }
        x2 <- pmin(pmax(x, 0), 1)
        asin(sqrt(x2))
    }

    ## ---- CLR matrix for disk (sensitivity) ----
    clr.mat <- NULL
    if (isTRUE(do.clr)) {
        X.disk <- X[vertices, , drop = FALSE]
        log.X <- log(X.disk + pseudo.count)
        row.center <- rowMeans(log.X)
        clr.mat <- sweep(log.X, 1, row.center, FUN = "-")
    }

    ## ---- bin-analysis helper (local, avoids dependency on external function) ----
    quantile.bin.analysis <- function(x, dd, carrier) {

        ## Bin boundaries based on the sample used for abundance summaries
        if (isTRUE(clr.over == "all")) {
            idx.ab <- rep(TRUE, length(dd))
        } else {
            idx.ab <- carrier
        }
        if (sum(idx.ab) < 2L) return(NULL)

        probs <- seq(0, 1, length.out = n.bins + 1L)
        brks <- unique(stats::quantile(dd[idx.ab], probs = probs, na.rm = TRUE, type = 7))
        if (length(brks) < 3L) return(NULL)

        bin.ab <- cut(dd[idx.ab], breaks = brks, include.lowest = TRUE)
        x.ab <- x[idx.ab]
        d.ab <- dd[idx.ab]

        mid <- tapply(d.ab, bin.ab, function(z) mean(range(z)))
        n.ab <- tapply(x.ab, bin.ab, length)
        med <- tapply(x.ab, bin.ab, stats::median)
        se <- tapply(x.ab, bin.ab, function(z) {
            if (length(z) <= 1L) return(NA_real_)
            stats::sd(z) / sqrt(length(z))
        })

        bins.df <- data.frame(
            bin = names(mid),
            d.mid = as.numeric(mid),
            n.abundance = as.integer(n.ab),
            abundance = as.numeric(med),
            abundance.se = as.numeric(se),
            stringsAsFactors = FALSE
        )

        ## Presence bins on all points using same breaks
        if (isTRUE(with.presence.bins)) {
            bin.all <- cut(dd, breaks = brks, include.lowest = TRUE)
            n.all <- tapply(dd, bin.all, length)
            prev <- tapply(carrier, bin.all, mean)
            mid.all <- tapply(dd, bin.all, function(z) mean(range(z)))
            bins.df$n.total <- as.integer(n.all[bins.df$bin])
            bins.df$presence.p <- as.numeric(prev[bins.df$bin])
            bins.df$d.mid.all <- as.numeric(mid.all[bins.df$bin])
        }

        ## Drop sparse bins
        keep.bin <- bins.df$n.abundance >= min.per.bin
        if (isTRUE(with.presence.bins)) keep.bin <- keep.bin & (bins.df$n.total >= min.per.bin)
        bins.df <- bins.df[keep.bin, , drop = FALSE]
        if (nrow(bins.df) < 3L) return(NULL)

        ## Trend tests
        tests <- list()

        if (bin.trend.test == "spearman") {
            ct <- suppressWarnings(stats::cor.test(bins.df$d.mid, bins.df$abundance, method = "spearman", exact = FALSE))
            tests$abundance <- list(method = "spearman", estimate = unname(ct$estimate), p.value = ct$p.value)
        } else {
            fit <- stats::lm(abundance ~ d.mid, data = bins.df, weights = n.abundance)
            sm <- summary(fit)$coefficients["d.mid", ]
            tests$abundance <- list(method = "wls", slope = unname(sm["Estimate"]), p.value = unname(sm["Pr(>|t|)"]))
        }

        if (isTRUE(with.presence.bins)) {
            if (bin.trend.test == "spearman") {
                ct <- suppressWarnings(stats::cor.test(bins.df$d.mid.all, bins.df$presence.p, method = "spearman", exact = FALSE))
                tests$presence <- list(method = "spearman", estimate = unname(ct$estimate), p.value = ct$p.value)
            } else {
                fit <- stats::lm(presence.p ~ d.mid.all, data = bins.df, weights = n.total)
                sm <- summary(fit)$coefficients["d.mid.all", ]
                tests$presence <- list(method = "wls", slope = unname(sm["Estimate"]), p.value = unname(sm["Pr(>|t|)"]))
            }
        }

        list(bins = bins.df, tests = tests, breaks = brks)
    }

    ## ---- per-feature loop ----
    res.list <- vector("list", length(feats))
    names(res.list) <- feats

    bin.results <- NULL
    if (distance.transform == "quantile.bin") {
        bin.results <- vector("list", length(feats))
        names(bin.results) <- feats
    }

    for (j in seq_along(feats)) {
        f <- feats[j]

        x0 <- X[vertices, f]
        carrier <- x0 > eps

        n.all <- length(x0)
        n.car <- sum(carrier)

        ## presence model (raw or rank only; for quantile.bin we also provide bin-level presence tests)
        pres.p <- NA_real_
        pres.slope <- NA_real_
        if (distance.transform != "quantile.bin") {
            if (n.car > 0L && n.car < n.all) {
                fit.pres <- suppressWarnings(stats::glm(carrier ~ d.used, family = stats::binomial()))
                sm <- summary(fit.pres)
                pres.slope <- unname(stats::coef(fit.pres)["d.used"])
                pres.p <- sm$coefficients["d.used", "Pr(>|z|)"]
            }
        } else {
            ## still compute vertex-level logistic if you want (optional); keep it consistent:
            if (n.car > 0L && n.car < n.all) {
                fit.pres <- suppressWarnings(stats::glm(carrier ~ d, family = stats::binomial()))
                sm <- summary(fit.pres)
                pres.slope <- unname(stats::coef(fit.pres)["d"])
                pres.p <- sm$coefficients["d", "Pr(>|z|)"]
            }
        }

        ## carrier-only tests (primary transform; raw/rank use d.used; quantile.bin uses d)
        rho <- NA_real_
        rho.p <- NA_real_
        lm.slope <- NA_real_
        lm.p <- NA_real_

        if (n.car >= as.integer(min.carriers)) {
            y <- transform.vec(x0[carrier])
            dd <- if (distance.transform == "quantile.bin") d[carrier] else d.used[carrier]

            sp <- suppressWarnings(stats::cor.test(y, dd, method = "spearman", exact = FALSE))
            rho <- unname(sp$estimate)
            rho.p <- sp$p.value

            fit.lm <- stats::lm(y ~ dd)
            lm.slope <- unname(stats::coef(fit.lm)["dd"])
            lm.p <- summary(fit.lm)$coefficients["dd", "Pr(>|t|)"]
        }

        ## CLR sensitivity
        clr.rho <- NA_real_
        clr.rho.p <- NA_real_
        clr.lm.slope <- NA_real_
        clr.lm.p <- NA_real_

        if (isTRUE(do.clr)) {
            idx.use <- if (clr.over == "all") rep(TRUE, length(vertices)) else carrier
            if (sum(idx.use) >= 2L) {
                y.clr <- clr.mat[idx.use, f]
                d.clr <- if (distance.transform == "quantile.bin") d[idx.use] else d.used[idx.use]

                spc <- suppressWarnings(stats::cor.test(y.clr, d.clr, method = "spearman", exact = FALSE))
                clr.rho <- unname(spc$estimate)
                clr.rho.p <- spc$p.value

                fit.clr <- stats::lm(y.clr ~ d.clr)
                clr.lm.slope <- unname(stats::coef(fit.clr)["d.clr"])
                clr.lm.p <- summary(fit.clr)$coefficients["d.clr", "Pr(>|t|)"]
            }
        }

        ## Bin-mode: compute bin summaries/tests for logit and CLR
        if (distance.transform == "quantile.bin") {
            bin.entry <- list()

            if (n.car >= 2L) {
                ## logit among carriers (by construction); if you prefer, you can pass y including zeros
                y.logit <- transform.vec(x0)
                bin.entry$logit <- quantile.bin.analysis(
                    x = if (clr.over == "all") y.logit else y.logit, ## y.logit already defined for all; carrier controls subset
                    dd = d,
                    carrier = carrier
                )
            }

            if (isTRUE(do.clr)) {
                y.clr.all <- clr.mat[, f]
                bin.entry$clr <- quantile.bin.analysis(
                    x = y.clr.all,
                    dd = d,
                    carrier = carrier
                )
            }

            bin.results[[f]] <- bin.entry
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
            clr.over = if (isTRUE(do.clr)) clr.over else NA_character_,
            clr.spearman.rho = clr.rho,
            clr.spearman.p = clr.rho.p,
            clr.lm.slope = clr.lm.slope,
            clr.lm.p = clr.lm.p,
            stringsAsFactors = FALSE
        )
    }

    vertex.df <- do.call(rbind, res.list)
    rownames(vertex.df) <- vertex.df$feature

    ## ---- FDR adjustments (vertex-level) ----
    vertex.df$presence.p.adj <- stats::p.adjust(vertex.df$presence.p, method = "BH")
    vertex.df$spearman.p.adj <- stats::p.adjust(vertex.df$spearman.p, method = "BH")
    vertex.df$lm.p.adj <- stats::p.adjust(vertex.df$lm.p, method = "BH")
    if (isTRUE(do.clr)) {
        vertex.df$clr.spearman.p.adj <- stats::p.adjust(vertex.df$clr.spearman.p, method = "BH")
        vertex.df$clr.lm.p.adj <- stats::p.adjust(vertex.df$clr.lm.p, method = "BH")
    }

    ## ---- FDR adjustments for bin-level tests (if present) ----
    if (distance.transform == "quantile.bin") {
        ## Extract p-values across features for each family (logit abundance, clr abundance, presence)
        extract.p <- function(which.slot, which.test) {
            v <- rep(NA_real_, length(feats))
            names(v) <- feats
            for (f in feats) {
                obj <- bin.results[[f]][[which.slot]]
                if (!is.null(obj) && !is.null(obj$tests[[which.test]])) {
                    v[f] <- obj$tests[[which.test]]$p.value
                }
            }
            v
        }

        p.logit.ab <- extract.p("logit", "abundance")
        p.clr.ab <- if (isTRUE(do.clr)) extract.p("clr", "abundance") else NULL
        p.pres <- if (isTRUE(with.presence.bins)) extract.p("clr", "presence") else NULL

        padj.logit.ab <- stats::p.adjust(p.logit.ab, method = "BH")
        if (!is.null(p.clr.ab)) padj.clr.ab <- stats::p.adjust(p.clr.ab, method = "BH") else padj.clr.ab <- NULL
        if (!is.null(p.pres)) padj.pres <- stats::p.adjust(p.pres, method = "BH") else padj.pres <- NULL

        ## Store adjusted p-values back into bin.results
        for (f in feats) {
            if (!is.null(bin.results[[f]]$logit) && !is.null(bin.results[[f]]$logit$tests$abundance)) {
                bin.results[[f]]$logit$tests$abundance$p.adj <- padj.logit.ab[f]
            }
            if (isTRUE(do.clr) && !is.null(bin.results[[f]]$clr) && !is.null(bin.results[[f]]$clr$tests$abundance)) {
                bin.results[[f]]$clr$tests$abundance$p.adj <- padj.clr.ab[f]
            }
            if (isTRUE(with.presence.bins) && !is.null(bin.results[[f]]$clr) && !is.null(bin.results[[f]]$clr$tests$presence)) {
                bin.results[[f]]$clr$tests$presence$p.adj <- padj.pres[f]
            }
        }
    }

    list(
        vertex.results = vertex.df,
        bin.results = bin.results,
        meta = list(
            distance.transform = distance.transform,
            n.bins = if (distance.transform == "quantile.bin") n.bins else NA_integer_,
            min.per.bin = if (distance.transform == "quantile.bin") min.per.bin else NA_integer_,
            bin.trend.test = if (distance.transform == "quantile.bin") bin.trend.test else NA_character_
        )
    )
}
