#' Two-channel (detection + intensity) preprocessing for Metabolon-style metabolomics matrices
#'
#' @description
#' Preprocesses a metabolite matrix using a two-channel representation:
#'   (i) an intensity channel (continuous) for sufficiently observed metabolites, and
#'   (ii) a detection channel (binary) for moderately sparse metabolites (and optionally very sparse).
#'
#' This avoids the common pitfall where replacing NA with 0 and then z-scoring causes
#' missingness to dominate Euclidean geometry.
#'
#' @details
#' The pipeline:
#' \enumerate{
#'   \item Compute per-feature missingness and robust spread; drop near-constants.
#'   \item Split features into tiers by missingness:
#'     \itemize{
#'       \item Tier A: dense (default p.miss <= 0.20) -> intensity only
#'       \item Tier B: moderately sparse (0.20 < p.miss <= 0.80) -> detection + intensity
#'       \item Tier C: very sparse (p.miss > 0.80 or too few observations) -> dropped by default
#'     }
#'   \item For intensity columns: impute missing values with a low value consistent with censoring,
#'         then apply scaled asinh, winsorization, and robust z-scoring.
#'   \item For detection columns: create 0/1 indicator and z-scale to unit variance,
#'         optionally down-weight via det.weight.
#'   \item Produce diagnostics (tables + plots + optional PCA).
#' }
#'
#' @param S Numeric matrix/data.frame: samples in rows, metabolites in columns.
#' @param p.miss.A Numeric in (0,1): Tier A threshold for missingness.
#' @param p.miss.B Numeric in (0,1): Tier B threshold for missingness (Tier C above this).
#' @param min.obs.int Integer: minimum non-missing observations required to keep an intensity feature.
#' @param det.rate.bounds Numeric length-2: allowed detection rate range for keeping detection indicators.
#' @param keep.tierC.det Logical: if TRUE, keep Tier C as detection-only (subject to det.rate.bounds).
#' @param impute.method Character: "half.min.pos", "q01", or "q05" (per-feature low-value imputation).
#' @param a Scale for asinh: "median.pos" (per-feature) or a single positive numeric scalar.
#' @param winsor.probs Numeric length-2: winsorization probabilities.
#' @param det.weight Numeric >= 0: multiplier applied to standardized detection columns.
#' @param do.pca Logical: compute PCA on the final matrix for diagnostics.
#' @param pca.max.rank Integer: maximum number of PCs to compute (prcomp rank cap).
#' @param meta Optional data.frame with sample metadata aligned to rows of S.
#' @param meta.color Optional character: a column name in meta to color the PCA scatter.
#' @param make.plots Logical: if TRUE, produce diagnostic plots.
#' @param plot.file Optional file path to write a multi-page PDF of diagnostic plots.
#' @param verbose Logical: print a short preprocessing summary.
#'
#' @return A list with elements:
#' \itemize{
#'   \item X: processed matrix (intensity + detection columns)
#'   \item feature.table: per-feature diagnostics and tier assignment
#'   \item name.map: mapping from raw feature names to safe names used internally
#'   \item pca: prcomp object (or NULL)
#'   \item diagnostics: list of summary stats and preprocessing parameters
#' }
#'
#' @examples
#' ## res <- preprocess.metabolon.twoplane(mb_ZB$S, plot.file = "metabolon_preproc_diagnostics.pdf")
#' ## X <- res$X
#'
#' @export
preprocess.metabolon.twoplane <- function(S,
                                         p.miss.A = 0.20,
                                         p.miss.B = 0.80,
                                         min.obs.int = 30L,
                                         det.rate.bounds = c(0.05, 0.95),
                                         keep.tierC.det = FALSE,
                                         impute.method = c("half.min.pos", "q01", "q05"),
                                         a = "median.pos",
                                         winsor.probs = c(0.01, 0.99),
                                         det.weight = 1.0,
                                         do.pca = TRUE,
                                         pca.max.rank = 50L,
                                         meta = NULL,
                                         meta.color = NULL,
                                         make.plots = TRUE,
                                         plot.file = NULL,
                                         verbose = TRUE) {

    ## ---- helpers (self-contained, CRAN-safe) ----
    winsorize.vec <- function(x, probs = c(0.01, 0.99)) {
        q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
        x[x < q[1]] <- q[1]
        x[x > q[2]] <- q[2]
        x
    }

    robust.zscore.vec <- function(x, eps = 1e-8) {
        m <- stats::median(x, na.rm = TRUE)
        s <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
        if (!is.finite(s) || s < eps) s <- 1.0
        (x - m) / s
    }

    sanitize.feature.names <- function(nm) {
        ## make syntactically valid + unique, while keeping a mapping to raw names
        if (is.null(nm)) nm <- paste0("V", seq_len(ncol(S)))
        make.names(nm, unique = TRUE)
    }

    impute.low.value <- function(v, method = "half.min.pos") {
        v <- v[is.finite(v)]
        if (length(v) == 0L) return(0)
        vpos <- v[v > 0]
        if (method == "half.min.pos") {
            if (length(vpos) > 0L) return(0.5 * min(vpos))
            return(stats::quantile(v, probs = 0.01, names = FALSE, na.rm = TRUE, type = 7))
        }
        if (method == "q01") return(stats::quantile(v, probs = 0.01, names = FALSE, na.rm = TRUE, type = 7))
        if (method == "q05") return(stats::quantile(v, probs = 0.05, names = FALSE, na.rm = TRUE, type = 7))
        stop("Unknown impute.method.")
    }

    preprocess.matrix.asinh.internal <- function(X,
                                                 a = "median.pos",
                                                 winsor.probs = c(0.01, 0.99)) {
        ## X must be numeric matrix, ideally finite (imputed upstream)
        X <- as.matrix(X)
        storage.mode(X) <- "double"

        ## replace any remaining non-finite with 0 (should be rare)
        X[!is.finite(X)] <- 0

        ## choose scale per feature
        aa <- rep(1.0, ncol(X))
        if (!is.null(colnames(X))) names(aa) <- colnames(X)

        if (is.character(a)) {
            if (length(a) != 1L || is.na(a) || a != "median.pos") stop("`a` must be \"median.pos\" or a single numeric > 0.")
            for (j in seq_len(ncol(X))) {
                v <- X[, j]
                vpos <- v[v > 0 & is.finite(v)]
                aa[j] <- if (length(vpos) > 0L) stats::median(vpos) else 1.0
                if (!is.finite(aa[j]) || aa[j] <= 0) aa[j] <- 1.0
            }
        } else {
            if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 0) stop("`a` numeric value must be a single finite number > 0.")
            aa[] <- a
        }

        aa[!is.finite(aa) | aa <= 0] <- 1.0

        for (j in seq_len(ncol(X))) X[, j] <- asinh(X[, j] / aa[j])
        for (j in seq_len(ncol(X))) X[, j] <- winsorize.vec(X[, j], probs = winsor.probs)
        for (j in seq_len(ncol(X))) X[, j] <- robust.zscore.vec(X[, j])

        attr(X, "asinh.scale") <- aa
        X
    }

    zscale.binary <- function(x, eps = 1e-8) {
        ## center and scale to unit variance using Bernoulli variance p(1-p)
        p <- mean(x)
        s <- sqrt(p * (1 - p))
        if (!is.finite(s) || s < eps) s <- 1.0
        (x - p) / s
    }

    ## ---- validate inputs ----
    if (missing(S) || is.null(S)) stop("`S` must be provided.")
    S0 <- tryCatch(as.matrix(S), error = function(e) NULL)
    if (is.null(S0)) stop("`S` must be coercible to a matrix via as.matrix().")
    suppressWarnings(storage.mode(S0) <- "double")
    if (!is.numeric(S0)) stop("`S` must be numeric (or coercible to numeric).")
    if (length(dim(S0)) != 2L) stop("`S` must be a 2D matrix.")
    if (nrow(S0) < 2L) stop("`S` must have at least 2 rows.")
    if (ncol(S0) < 1L) stop("`S` must have at least 1 column.")

    if (!is.numeric(p.miss.A) || length(p.miss.A) != 1L || p.miss.A <= 0 || p.miss.A >= 1) stop("`p.miss.A` must be in (0,1).")
    if (!is.numeric(p.miss.B) || length(p.miss.B) != 1L || p.miss.B <= 0 || p.miss.B >= 1) stop("`p.miss.B` must be in (0,1).")
    if (p.miss.A >= p.miss.B) stop("Require p.miss.A < p.miss.B.")
    if (!is.numeric(det.rate.bounds) || length(det.rate.bounds) != 2L) stop("`det.rate.bounds` must be numeric length 2.")
    if (det.rate.bounds[1] < 0 || det.rate.bounds[2] > 1 || det.rate.bounds[1] >= det.rate.bounds[2]) {
        stop("`det.rate.bounds` must satisfy 0 <= lo < hi <= 1.")
    }
    if (!is.numeric(winsor.probs) || length(winsor.probs) != 2L) stop("`winsor.probs` must be numeric length 2.")
    if (winsor.probs[1] < 0 || winsor.probs[2] > 1 || winsor.probs[1] >= winsor.probs[2]) stop("`winsor.probs` must satisfy 0 <= p1 < p2 <= 1.")
    if (!is.numeric(det.weight) || length(det.weight) != 1L || !is.finite(det.weight) || det.weight < 0) stop("`det.weight` must be a single finite number >= 0.")
    if (!is.numeric(min.obs.int) || length(min.obs.int) != 1L || min.obs.int < 2) stop("`min.obs.int` must be an integer >= 2.")

    impute.method <- match.arg(impute.method)

    ## ---- name handling ----
    col.names.raw <- colnames(S0)
    col.names.safe <- sanitize.feature.names(col.names.raw)
    colnames(S0) <- col.names.safe

    name.map <- data.frame(
        feature.name.raw = if (is.null(col.names.raw)) col.names.safe else col.names.raw,
        feature.name = col.names.safe,
        stringsAsFactors = FALSE
    )

    ## ---- compute missingness and basic per-feature stats ----
    is.obs <- is.finite(S0)
    p.miss <- 1.0 - colMeans(is.obs)
    n.obs <- colSums(is.obs)
    det.rate <- colMeans(is.obs)

    ## robust spread among observed values (raw scale)
    mad.raw <- rep(NA_real_, ncol(S0))
    iqr.raw <- rep(NA_real_, ncol(S0))
    for (j in seq_len(ncol(S0))) {
        v <- S0[is.obs[, j], j]
        if (length(v) >= 2L) {
            mad.raw[j] <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
            iqr.raw[j] <- stats::IQR(v, na.rm = TRUE, type = 7)
        }
    }

    ## near-constant flag (raw)
    is.constant <- (!is.finite(mad.raw) | mad.raw <= 0) & (!is.finite(iqr.raw) | iqr.raw <= 0)

    ## ---- tiering ----
    tier <- rep("C", ncol(S0))
    tier[p.miss <= p.miss.A] <- "A"
    tier[p.miss > p.miss.A & p.miss <= p.miss.B] <- "B"

    ## enforce minimal observations for intensity
    tier[n.obs < min.obs.int] <- "C"

    ## drop constants from intensity eligibility
    tier[is.constant] <- "C"

    ## sets
    idx.A <- which(tier == "A")
    idx.B <- which(tier == "B")
    idx.C <- which(tier == "C")

    ## ---- choose which detection indicators to keep ----
    keep.det.B <- idx.B[det.rate[idx.B] >= det.rate.bounds[1] & det.rate[idx.B] <= det.rate.bounds[2]]
    keep.det.C <- integer(0)
    if (isTRUE(keep.tierC.det)) {
        keep.det.C <- idx.C[det.rate[idx.C] >= det.rate.bounds[1] & det.rate[idx.C] <= det.rate.bounds[2]]
    }
    keep.det <- sort(unique(c(keep.det.B, keep.det.C)))

    ## ---- choose which intensity features to keep ----
    keep.int <- sort(unique(c(idx.A, idx.B)))
    keep.int <- keep.int[keep.int %in% which(n.obs >= min.obs.int)]
    keep.int <- keep.int[!is.constant[keep.int]]

    ## ---- build intensity matrix (with low-value imputation) ----
    X.int.raw <- S0[, keep.int, drop = FALSE]
    impute.val <- rep(NA_real_, ncol(S0)); names(impute.val) <- colnames(S0)

    for (j in seq_len(ncol(X.int.raw))) {
        jj <- keep.int[j]
        v <- S0[, jj]
        v.obs <- v[is.finite(v)]
        fill <- impute.low.value(v.obs, method = impute.method)
        if (!is.finite(fill)) fill <- 0
        ## enforce a small non-negative fill in typical Metabolon-positive matrices
        if (fill < 0) fill <- stats::quantile(v.obs, probs = 0.01, names = FALSE, na.rm = TRUE, type = 7)
        if (!is.finite(fill)) fill <- 0
        impute.val[jj] <- fill
        v[!is.finite(v)] <- fill
        X.int.raw[, j] <- v
    }

    ## rename intensity columns
    colnames(X.int.raw) <- paste0(colnames(X.int.raw), ".int")

    ## preprocess intensity with scaled asinh + winsor + robust z
    X.int <- preprocess.matrix.asinh.internal(X.int.raw, a = a, winsor.probs = winsor.probs)
    aa <- attr(X.int, "asinh.scale")

    ## ---- build detection matrix ----
    X.det <- NULL
    if (length(keep.det) > 0L) {
        X.det.raw <- matrix(0.0, nrow = nrow(S0), ncol = length(keep.det))
        colnames(X.det.raw) <- paste0(colnames(S0)[keep.det], ".det")
        for (k in seq_along(keep.det)) {
            j <- keep.det[k]
            d <- as.numeric(is.finite(S0[, j]))
            d <- zscale.binary(d)
            X.det.raw[, k] <- det.weight * d
        }
        X.det <- X.det.raw
    }

    ## ---- combine channels ----
    X <- if (is.null(X.det)) X.int else cbind(X.int, X.det)
    rownames(X) <- rownames(S0)

    ## ---- optional PCA for diagnostics ----
    pca <- NULL
    pca.info <- NULL
    if (isTRUE(do.pca)) {
        ## cap rank to avoid needless compute
        rcap <- min(as.integer(pca.max.rank), nrow(X) - 1L, ncol(X))
        if (rcap >= 2L) {
            ## prcomp expects no NA
            pca <- stats::prcomp(X, center = FALSE, scale. = FALSE, rank. = rcap)
            ve <- (pca$sdev^2) / sum(pca$sdev^2)
            pca.info <- list(var.explained = ve, var.explained.cum = cumsum(ve))
        }
    }

    ## ---- feature table ----
    feature.table <- data.frame(
        feature.name = colnames(S0),
        feature.name.raw = name.map$feature.name.raw[match(colnames(S0), name.map$feature.name)],
        p.miss = p.miss,
        n.obs = n.obs,
        det.rate = det.rate,
        mad.raw = mad.raw,
        iqr.raw = iqr.raw,
        tier = tier,
        keep.int = FALSE,
        keep.det = FALSE,
        impute.value = impute.val,
        stringsAsFactors = FALSE
    )
    feature.table$keep.int[keep.int] <- TRUE
    feature.table$keep.det[keep.det] <- TRUE

    ## add asinh scale for intensity features
    feature.table$asinh.scale <- NA_real_
    if (!is.null(aa) && length(keep.int) == length(aa)) {
        ## aa is named by intensity column names; map back
        aa.base <- aa
        names(aa.base) <- sub("\\.int$", "", names(aa.base))
        feature.table$asinh.scale[match(names(aa.base), feature.table$feature.name)] <- as.numeric(aa.base)
    }

    ## ---- diagnostics summary ----
    diagnostics <- list(
        n.samples = nrow(S0),
        n.features.total = ncol(S0),
        n.tier.A = length(idx.A),
        n.tier.B = length(idx.B),
        n.tier.C = length(idx.C),
        n.intensity = ncol(X.int),
        n.detection = if (is.null(X.det)) 0L else ncol(X.det),
        n.final = ncol(X),
        impute.method = impute.method,
        a = a,
        winsor.probs = winsor.probs,
        p.miss.A = p.miss.A,
        p.miss.B = p.miss.B,
        det.rate.bounds = det.rate.bounds,
        det.weight = det.weight,
        constant.features = sum(is.constant, na.rm = TRUE)
    )

    if (isTRUE(verbose)) {
        msg <- paste0(
            "Two-channel preprocessing summary:\n",
            "  samples: ", diagnostics$n.samples, "\n",
            "  features total: ", diagnostics$n.features.total, "\n",
            "  tiers: A=", diagnostics$n.tier.A, ", B=", diagnostics$n.tier.B, ", C=", diagnostics$n.tier.C, "\n",
            "  kept: intensity=", diagnostics$n.intensity, ", detection=", diagnostics$n.detection, ", final=", diagnostics$n.final, "\n"
        )
        cat(msg)
    }

    ## ---- diagnostic plots ----
    if (isTRUE(make.plots)) {
        if (!is.null(plot.file)) grDevices::pdf(plot.file, width = 10, height = 7)

        ## (1) Missingness histogram with tier cutoffs
        graphics::hist(
            p.miss,
            breaks = 30,
            main = "Per-metabolite missingness",
            xlab = "p.miss",
            col = "gray",
            border = "white"
        )
        graphics::abline(v = p.miss.A, lwd = 2, lty = 2)
        graphics::abline(v = p.miss.B, lwd = 2, lty = 2)
        graphics::mtext(paste0("p.miss.A=", p.miss.A, "   p.miss.B=", p.miss.B), side = 3, line = 0.2)

        ## (2) Detection rate vs missingness (raw)
        graphics::plot(
            det.rate, p.miss,
            pch = 16, cex = 0.6,
            xlab = "Detection rate (1 - p.miss)",
            ylab = "p.miss",
            main = "Detection rate vs missingness"
        )
        graphics::abline(h = p.miss.A, lty = 2)
        graphics::abline(h = p.miss.B, lty = 2)
        graphics::abline(v = det.rate.bounds[1], lty = 3)
        graphics::abline(v = det.rate.bounds[2], lty = 3)

        ## (3) Asinh scales (intensity features)
        if (sum(feature.table$keep.int) > 0L) {
            aa.plot <- feature.table$asinh.scale[feature.table$keep.int]
            aa.plot <- aa.plot[is.finite(aa.plot) & aa.plot > 0]
            if (length(aa.plot) > 0L) {
                graphics::hist(
                    log10(aa.plot),
                    breaks = 30,
                    main = "log10(asinh scale) for intensity features",
                    xlab = "log10(scale)",
                    col = "gray",
                    border = "white"
                )
            } else {
                graphics::plot.new()
                graphics::title("No finite asinh scales to plot.")
            }
        }

        ## (4) Sample row-norms after preprocessing
        row.norm <- sqrt(rowSums(X^2))
        graphics::hist(
            row.norm,
            breaks = 30,
            main = "Row norms after preprocessing",
            xlab = "||x_i||_2",
            col = "gray",
            border = "white"
        )

        ## (5) PCA scree + PC1/PC2 scatter (optional)
        if (!is.null(pca)) {
            ve <- pca.info$var.explained
            k <- seq_along(ve)
            graphics::plot(k, ve, type = "b", pch = 16,
                           xlab = "PC", ylab = "Variance explained",
                           main = "PCA scree (on processed matrix)")

            scores <- pca$x
            x1 <- scores[, 1]
            x2 <- scores[, 2]
            main.txt <- "PCA: PC1 vs PC2"

            if (!is.null(meta) && !is.null(meta.color) && meta.color %in% colnames(meta)) {
                g <- meta[[meta.color]]
                g <- as.factor(g)
                col.vec <- as.integer(g)
                graphics::plot(x1, x2, pch = 16, cex = 0.7, col = col.vec,
                               xlab = "PC1", ylab = "PC2", main = paste0(main.txt, " (colored by ", meta.color, ")"))
                graphics::legend("topright", legend = levels(g), col = seq_along(levels(g)), pch = 16, cex = 0.7)
            } else {
                graphics::plot(x1, x2, pch = 16, cex = 0.7,
                               xlab = "PC1", ylab = "PC2", main = main.txt)
            }
        }

        if (!is.null(plot.file)) grDevices::dev.off()
    }

    ## ---- return ----
    list(
        X = X,
        feature.table = feature.table[order(feature.table$p.miss), ],
        name.map = name.map,
        pca = pca,
        diagnostics = diagnostics
    )
}


#' Diagnose geometry of a two-channel (intensity + detection) matrix
#'
#' @description
#' Computes diagnostics to check whether detection indicators dominate geometry:
#'   - per-sample L2 norm contribution from .int vs .det blocks
#'   - per-pair squared-distance contribution from .int vs .det blocks (random pairs)
#'   - correlation of PC1 with detection-rate proxies (optional)
#'
#' @param X Numeric matrix produced by preprocess.metabolon.twoplane().
#' @param n.pairs Integer number of random sample pairs for distance decomposition.
#' @param seed Integer seed for reproducibility.
#' @return A list with summary stats and per-sample contributions.
#' @export
diagnose.twoplane.geometry <- function(X, n.pairs = 5000L, seed = 1L) {
    if (missing(X) || is.null(X)) stop("`X` must be provided.")
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    if (nrow(X) < 3L) stop("`X` must have >= 3 rows.")
    if (ncol(X) < 2L) stop("`X` must have >= 2 cols.")

    cn <- colnames(X)
    if (is.null(cn)) stop("`X` must have colnames ending with .int / .det.")
    idx.int <- grepl("\\.int$", cn)
    idx.det <- grepl("\\.det$", cn)
    if (!any(idx.int)) stop("No .int columns found.")
    if (!any(idx.det)) warning("No .det columns found.")

    X.int <- X[, idx.int, drop = FALSE]
    X.det <- if (any(idx.det)) X[, idx.det, drop = FALSE] else NULL

    ## per-sample L2 norms
    row.norm.int <- sqrt(rowSums(X.int^2))
    row.norm.det <- if (!is.null(X.det)) sqrt(rowSums(X.det^2)) else rep(0, nrow(X))
    row.norm.all <- sqrt(rowSums(X^2))

    ## fraction of norm (note: norm fractions, not variance explained)
    frac.det.norm <- ifelse(row.norm.all > 0, row.norm.det / row.norm.all, NA_real_)
    frac.int.norm <- ifelse(row.norm.all > 0, row.norm.int / row.norm.all, NA_real_)

    ## pairwise squared-distance decomposition (random pairs)
    set.seed(seed)
    n <- nrow(X)
    n.pairs <- as.integer(n.pairs)
    i1 <- sample.int(n, size = n.pairs, replace = TRUE)
    i2 <- sample.int(n, size = n.pairs, replace = TRUE)

    d2.int <- rowSums((X.int[i1, , drop = FALSE] - X.int[i2, , drop = FALSE])^2)
    d2.det <- if (!is.null(X.det)) rowSums((X.det[i1, , drop = FALSE] - X.det[i2, , drop = FALSE])^2) else rep(0, n.pairs)
    d2.all <- d2.int + d2.det
    frac.det.d2 <- ifelse(d2.all > 0, d2.det / d2.all, NA_real_)

    ## summarize
    out <- list(
        n.samples = n,
        p.int = ncol(X.int),
        p.det = if (is.null(X.det)) 0L else ncol(X.det),
        frac.det.norm.summary = stats::quantile(frac.det.norm, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE),
        frac.det.d2.summary = stats::quantile(frac.det.d2, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE),
        per.sample = data.frame(
            row.norm.int = row.norm.int,
            row.norm.det = row.norm.det,
            row.norm.all = row.norm.all,
            frac.int.norm = frac.int.norm,
            frac.det.norm = frac.det.norm
        ),
        per.pair = data.frame(
            i1 = i1, i2 = i2,
            d2.int = d2.int,
            d2.det = d2.det,
            frac.det.d2 = frac.det.d2
        )
    )

    out
}

#' Compute Jaccard distance between two unique integer edge-code vectors
#'
#' @param e1 Integer vector of unique edge codes.
#' @param e2 Integer vector of unique edge codes.
#' @return Numeric scalar in [0,1].
#' @export
jaccard.edge.distance <- function(e1, e2) {
    e1 <- as.integer(e1); e2 <- as.integer(e2)
    if (length(e1) == 0L && length(e2) == 0L) return(0)
    if (length(e1) == 0L || length(e2) == 0L) return(1)
    inter <- sum(!is.na(match(e1, e2)))
    uni <- length(e1) + length(e2) - inter
    if (uni <= 0L) return(0)
    1 - (inter / uni)
}

#' Build an undirected edge list (and codes) from a kNN index matrix
#'
#' @description
#' Converts directed kNN relationships into undirected graphs:
#'   - mode="sym": edge if i->j OR j->i (symmetrized kNN)
#'   - mode="mutual": edge if i->j AND j->i (mutual kNN)
#'
#' @param nn.index Integer matrix n x k.max returned by FNN::get.knn()$nn.index
#' @param k Integer <= ncol(nn.index)
#' @param mode Character: "sym" or "mutual"
#' @return List(edge.mat = two-col integer matrix with i<j, edge.code = integer codes)
#' @export
edges.from.nnindex <- function(nn.index, k, mode = c("sym", "mutual")) {
    mode <- match.arg(mode)
    nn.index <- as.matrix(nn.index)
    n <- nrow(nn.index)
    k <- as.integer(k)
    if (k < 1L || k > ncol(nn.index)) stop("Invalid k for nn.index.")
    if (n < 2L) stop("Need n>=2.")

    ## directed adjacency
    A <- matrix(FALSE, nrow = n, ncol = n)
    i.rep <- rep.int(seq_len(n), times = k)
    j.vec <- as.integer(nn.index[, seq_len(k), drop = FALSE])
    A[cbind(i.rep, j.vec)] <- TRUE
    diag(A) <- FALSE

    ## undirected
    U <- if (mode == "sym") (A | t(A)) else (A & t(A))
    diag(U) <- FALSE

    ## extract i<j edges
    ij <- which(U & upper.tri(U), arr.ind = TRUE)
    if (nrow(ij) == 0L) {
        return(list(edge.mat = matrix(integer(0), ncol = 2), edge.code = integer(0)))
    }
    edge.mat <- cbind(as.integer(ij[, 1]), as.integer(ij[, 2]))

    ## encode edge (i,j) with i<j as unique code
    edge.code <- as.integer((edge.mat[, 1] - 1L) * n + edge.mat[, 2])
    edge.code <- sort(unique(edge.code))

    list(edge.mat = edge.mat, edge.code = edge.code)
}

#' Compute kNN graph sequence diagnostics and edit-distance curve
#'
#' @param X Numeric matrix (rows are samples).
#' @param k.grid Integer vector (ideally consecutive) of k values to evaluate.
#' @param mode "sym" or "mutual"
#' @param min.lcc.frac Minimum largest-connected-component fraction defining the "connected tail".
#' @param eps Relative tolerance for selecting smallest k within (1+eps)*min edit distance.
#' @return List with a data.frame curve and chosen k.
#' @export
knn.edit.distance.curve <- function(X,
                                   k.grid = 5:120,
                                   mode = c("sym", "mutual"),
                                   min.lcc.frac = 1.0,
                                   eps = 0.05) {
    mode <- match.arg(mode)
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    n <- nrow(X)
    if (n < 10L) stop("Too few samples for kNN sweep.")

    k.grid <- sort(unique(as.integer(k.grid)))
    k.grid <- k.grid[k.grid >= 1L & k.grid <= (n - 1L)]
    if (length(k.grid) < 3L) stop("k.grid too short after filtering.")

    k.max <- max(k.grid)
    nn <- FNN::get.knn(X, k = k.max)$nn.index

    ## compute per-k stats + cache edge codes for edit distance
    n.edges <- integer(length(k.grid))
    n.comp <- integer(length(k.grid))
    lcc.frac <- numeric(length(k.grid))
    edge.code.list <- vector("list", length(k.grid))

    for (ii in seq_along(k.grid)) {
        k <- k.grid[ii]
        eg <- edges.from.nnindex(nn, k = k, mode = mode)
        edge.code.list[[ii]] <- eg$edge.code
        n.edges[ii] <- length(eg$edge.code)

        ## IMPORTANT: include isolates by constructing an empty n-vertex graph
        g <- igraph::make_empty_graph(n = n, directed = FALSE)
        if (nrow(eg$edge.mat) > 0L) {
            ## add_edges expects a vector: v1, v2, v3, v4, ...
            g <- igraph::add_edges(g, as.vector(t(eg$edge.mat)))
        }

        comp <- igraph::components(g)
        n.comp[ii] <- comp$no
        lcc.frac[ii] <- max(comp$csize) / n
    }

    ## edit distances between consecutive k values in k.grid
    edit.dist <- rep(NA_real_, length(k.grid))
    for (ii in seq_len(length(k.grid) - 1L)) {
        edit.dist[ii] <- jaccard.edge.distance(edge.code.list[[ii]], edge.code.list[[ii + 1L]])
    }

    curve <- data.frame(
        k = k.grid,
        n.edges = n.edges,
        n.components = n.comp,
        lcc.frac = lcc.frac,
        edit.dist.to.next = edit.dist
    )

    ## connected tail start
    ok.tail <- which(curve$lcc.frac >= min.lcc.frac)
    k.cc <- if (length(ok.tail) > 0L) curve$k[min(ok.tail)] else NA_integer_

    ## select k: smallest within tolerance of min edit distance in tail (exclude last row where edit is NA)
    k.sel <- NA_integer_
    if (is.finite(k.cc)) {
        tail.idx <- which(curve$k >= k.cc & is.finite(curve$edit.dist.to.next))
        if (length(tail.idx) > 0L) {
            m <- min(curve$edit.dist.to.next[tail.idx])
            thr <- (1 + eps) * m
            k.sel <- curve$k[min(tail.idx[curve$edit.dist.to.next[tail.idx] <= thr])]
        }
    }

    list(curve = curve, k.cc = k.cc, k.selected = k.sel, mode = mode,
         min.lcc.frac = min.lcc.frac, eps = eps)
}


#' Plot kNN edit-distance diagnostics
#'
#' @param res Output of knn.edit.distance.curve().
#' @export
plot.knn.edit.curve <- function(res) {
    curve <- res$curve
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    graphics::par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

    ## components / LCC
    graphics::plot(curve$k, curve$lcc.frac, type = "l",
                   xlab = "k", ylab = "Largest CC fraction",
                   main = paste0("Connectivity (mode=", res$mode, ")"))
    if (is.finite(res$k.cc)) graphics::abline(v = res$k.cc, lty = 2)
    if (is.finite(res$k.selected)) graphics::abline(v = res$k.selected, lty = 3)

    ## edges
    graphics::plot(curve$k, curve$n.edges, type = "l",
                   xlab = "k", ylab = "# edges",
                   main = "Edge count")

    ## edit distance
    graphics::plot(curve$k, curve$edit.dist.to.next, type = "l",
                   xlab = "k", ylab = "Jaccard edit dist to next k",
                   main = "Edit distance between consecutive graphs")
    if (is.finite(res$k.cc)) graphics::abline(v = res$k.cc, lty = 2)
    if (is.finite(res$k.selected)) graphics::abline(v = res$k.selected, lty = 3)
    graphics::mtext(paste0("k.cc=", res$k.cc, "   k.selected=", res$k.selected), side = 3, line = 0.2)
}


#' Winsorize a numeric vector
#'
#' @param x Numeric vector.
#' @param probs Numeric length-2 with lower/upper probabilities.
#' @return Winsorized numeric vector.
#' @export
winsorize.vec <- function(x, probs = c(0.01, 0.99)) {
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    x
}

#' Robust z-score a numeric vector using median and MAD
#'
#' @param x Numeric vector.
#' @param eps Small positive scalar for numerical stability.
#' @return Robustly standardized numeric vector.
#' @export
robust.zscore.vec <- function(x, eps = 1e-8) {
    m <- stats::median(x, na.rm = TRUE)
    s <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s < eps) s <- 1.0
    (x - m) / s
}

#' Replace missing (and non-positive) entries by 2/3 of min positive per metabolite
#'
#' @description
#' Implements the McMillan/Reid-style replacement used before log10 transforms:
#' zeros were replaced by 2/3 of the minimum detected value per metabolite.
#' Here we apply the same idea to NA (and optionally non-positive values).
#' :contentReference[oaicite:4]{index=4}
#'
#' @param S Numeric matrix (samples x metabolites).
#' @param frac Fraction of minimum positive used as replacement (default 2/3).
#' @param replace.nonpos Logical; if TRUE, values <= 0 are replaced too.
#' @return List with S.filled and per-feature delta.
#' @export
impute.two.thirds.min.pos <- function(S, frac = 2/3, replace.nonpos = TRUE) {
    S <- as.matrix(S)
    storage.mode(S) <- "double"

    p <- ncol(S)
    delta <- rep(NA_real_, p)
    names(delta) <- colnames(S)

    for (j in seq_len(p)) {
        v <- S[, j]
        ok <- is.finite(v)
        v.ok <- v[ok]

        v.pos <- v.ok[v.ok > 0]
        if (length(v.pos) > 0L) {
            delta[j] <- frac * min(v.pos)
        } else if (length(v.ok) > 0L) {
            ## fallback: small quantile of finite values
            delta[j] <- stats::quantile(v.ok, probs = 0.01, names = FALSE, na.rm = TRUE, type = 7)
        } else {
            delta[j] <- 1e-6
        }

        if (!is.finite(delta[j]) || delta[j] <= 0) delta[j] <- 1e-6

        bad <- !is.finite(v)
        if (isTRUE(replace.nonpos)) bad <- bad | (v <= 0)
        v[bad] <- delta[j]
        S[, j] <- v
    }

    list(S.filled = S, delta = delta)
}

#' Compute CLR-style log-ratio features (graph-friendly ratio embedding)
#'
#' @description
#' Returns centered log-ratio (CLR) features: log(x) minus per-sample mean/median log(x).
#' Euclidean distance on CLR features corresponds to Aitchison distance for compositions. :contentReference[oaicite:5]{index=5}
#'
#' @param S Numeric matrix (samples x metabolites).
#' @param log.base Numeric; use 10 for log10 (matches Reid LC-MS processing). :contentReference[oaicite:6]{index=6}
#' @param center Character: "mean" or "median" for per-sample centering.
#' @param winsor.probs Optional winsorization probs on log scale (set NULL to disable).
#' @param zscore Logical; if TRUE, robust z-score each metabolite after CLR.
#' @return Numeric matrix of CLR log-ratio features.
#' @export
logratio.clr.matrix <- function(S,
                               log.base = 10,
                               center = c("mean", "median"),
                               winsor.probs = c(0.01, 0.99),
                               zscore = TRUE) {
    center <- match.arg(center)

    imp <- impute.two.thirds.min.pos(S, frac = 2/3, replace.nonpos = TRUE)
    S.filled <- imp$S.filled

    L <- if (log.base == 10) log10(S.filled) else log(S.filled, base = log.base)

    if (!is.null(winsor.probs)) {
        for (j in seq_len(ncol(L))) L[, j] <- winsorize.vec(L[, j], probs = winsor.probs)
    }

    c0 <- if (center == "mean") rowMeans(L) else apply(L, 1, stats::median)
    X.clr <- L - c0

    if (isTRUE(zscore)) {
        for (j in seq_len(ncol(X.clr))) X.clr[, j] <- robust.zscore.vec(X.clr[, j])
    }

    X.clr
}

#' Select stable reference metabolites for ALR-style ratios
#'
#' @description
#' Picks reference metabolites with low missingness and low robust variance in log space.
#'
#' @param S Numeric matrix (samples x metabolites).
#' @param p.miss.max Maximum allowed missingness for candidates.
#' @param n.refs Number of references to return.
#' @param log.base Numeric; log base (10 recommended for Reid-style).
#' @return Character vector of selected reference metabolite names.
#' @export
select.stable.references <- function(S,
                                    p.miss.max = 0.05,
                                    n.refs = 10L,
                                    log.base = 10) {
    S <- as.matrix(S)
    storage.mode(S) <- "double"

    p.miss <- colMeans(!is.finite(S))
    keep <- which(p.miss <= p.miss.max)
    if (length(keep) < 2L) stop("Too few candidates under p.miss.max; relax threshold.")

    imp <- impute.two.thirds.min.pos(S[, keep, drop = FALSE], frac = 2/3, replace.nonpos = TRUE)
    S.filled <- imp$S.filled
    L <- if (log.base == 10) log10(S.filled) else log(S.filled, base = log.base)

    mad.log <- apply(L, 2, stats::mad, constant = 1.4826)
    ord <- order(mad.log, decreasing = FALSE, na.last = NA)

    refs <- colnames(S.filled)[ord][seq_len(min(as.integer(n.refs), length(ord)))]
    refs
}

#' Compute ALR-style log ratios to reference metabolites
#'
#' @param S Numeric matrix (samples x metabolites).
#' @param refs Character vector of reference metabolite names (must be colnames of S).
#' @param log.base Numeric; log base (10 recommended).
#' @param zscore Logical; if TRUE, robust z-score each ratio feature.
#' @return Numeric matrix of log ratios log(x/ref).
#' @export
logratio.to.refs.matrix <- function(S, refs, log.base = 10, zscore = TRUE) {
    S <- as.matrix(S)
    storage.mode(S) <- "double"
    if (is.null(colnames(S))) stop("`S` must have column names.")
    if (length(refs) < 1L) stop("`refs` must be non-empty.")
    if (any(!refs %in% colnames(S))) stop("Some `refs` not found in colnames(S).")

    imp <- impute.two.thirds.min.pos(S, frac = 2/3, replace.nonpos = TRUE)
    S.filled <- imp$S.filled
    L <- if (log.base == 10) log10(S.filled) else log(S.filled, base = log.base)

    X <- matrix(0.0, nrow = nrow(L), ncol = ncol(L) * length(refs))
    cn <- character(ncol(X))

    k <- 0L
    for (r in refs) {
        for (j in seq_len(ncol(L))) {
            k <- k + 1L
            X[, k] <- L[, j] - L[, r]
            cn[k] <- paste0(colnames(L)[j], ".over.", r)
        }
    }
    colnames(X) <- cn
    rownames(X) <- rownames(L)

    if (isTRUE(zscore)) {
        for (j in seq_len(ncol(X))) X[, j] <- robust.zscore.vec(X[, j])
    }

    X
}

#' CST mixing statistics on a graph (homophily, assortativity, conductance) with permutation null
#'
#' @description
#' Quantifies how much a categorical labeling (e.g., CST) mixes over a graph.
#' Provides:
#'   - weighted edge homophily (same-label edge fraction),
#'   - Newman nominal assortativity coefficient r (categorical),
#'   - per-label conductance (cut / min(vol, vol_complement)),
#'   - permutation z-scores and p-values (optional stratified permutation).
#'
#' @param g An igraph graph object (undirected recommended).
#' @param labels Vector of CST labels aligned to vertices (length vcount(g)).
#'   Missing labels allowed (NA); edges involving NA endpoints are ignored.
#' @param edge.weights Optional numeric vector of edge weights (length ecount(g)).
#'   If NULL, unweighted (w=1) is used.
#' @param n.perm Integer number of permutations for null distribution (0 disables).
#' @param perm.blocks Optional factor/character vector (length vcount(g)) for within-block permutations
#'   (e.g., subject id, batch). If NULL, global permutation is used.
#' @param seed Integer RNG seed.
#'
#' @return List with observed metrics, per-label metrics, and permutation summaries.
#' @export
cst.graph.mixing.stats <- function(g,
                                  labels,
                                  edge.weights = NULL,
                                  n.perm = 200L,
                                  perm.blocks = NULL,
                                  seed = 1L) {
    if (!inherits(g, "igraph")) stop("`g` must be an igraph object.")
    n <- igraph::vcount(g)
    m <- igraph::ecount(g)
    if (length(labels) != n) stop("`labels` must have length vcount(g).")

    labels <- as.character(labels)
    ok.v <- !is.na(labels)

    ## edge list
    ends <- igraph::ends(g, igraph::E(g), names = FALSE)
    u <- ends[, 1]
    v <- ends[, 2]

    lab.u <- labels[u]
    lab.v <- labels[v]
    ok.e <- !is.na(lab.u) & !is.na(lab.v)

    if (is.null(edge.weights)) {
        w <- rep(1.0, m)
    } else {
        if (length(edge.weights) != m) stop("`edge.weights` length must equal ecount(g).")
        w <- as.double(edge.weights)
        if (any(!is.finite(w)) || any(w < 0)) stop("`edge.weights` must be finite and non-negative.")
    }

    ## weighted homophily on labeled edges
    w.ok <- w[ok.e]
    same <- (lab.u[ok.e] == lab.v[ok.e])
    homophily <- if (sum(w.ok) > 0) sum(w.ok[same]) / sum(w.ok) else NA_real_

    ## assortativity (nominal) on induced labeled subgraph
    ## igraph assortativity_nominal is based on Newman mixing matrix (categorical). :contentReference[oaicite:3]{index=3}
    r.assort <- NA_real_
    if (sum(ok.v) >= 3L && sum(ok.e) >= 1L) {
        g.sub <- igraph::induced_subgraph(g, vids = which(ok.v))
        types <- as.integer(as.factor(labels[ok.v]))
        r.assort <- igraph::assortativity_nominal(g.sub, types = types, directed = FALSE)
    }

    ## conductance per label (on labeled induced subgraph)
    cond.by.label <- NULL
    cond.summary <- NULL
    if (sum(ok.v) >= 3L) {
        g.sub <- igraph::induced_subgraph(g, vids = which(ok.v))
        labs.sub <- labels[ok.v]
        lev <- sort(unique(labs.sub))

        ## degrees/strengths
        if (is.null(edge.weights)) {
            deg <- igraph::degree(g.sub)
        } else {
            ## weights on induced subgraph edges
            ## pull weights for edges in g.sub by mapping original edge ids
            ## easiest: recompute unweighted conductance by default when weighted graphs are not explicit
            deg <- igraph::degree(g.sub)
        }

        ## edge endpoints in subgraph
        ends2 <- igraph::ends(g.sub, igraph::E(g.sub), names = FALSE)
        u2 <- ends2[, 1]
        v2 <- ends2[, 2]
        lab.u2 <- labs.sub[u2]
        lab.v2 <- labs.sub[v2]

        ## compute cut_l = #edges from label l to not-l
        ## compute vol_l = sum(deg in label l)
        vol.total <- sum(deg)
        cut.l <- setNames(rep(0.0, length(lev)), lev)
        vol.l <- setNames(rep(0.0, length(lev)), lev)

        for (l in lev) {
            in.l <- (labs.sub == l)
            vol.l[l] <- sum(deg[in.l])
        }
        for (e in seq_along(u2)) {
            a <- lab.u2[e]; b <- lab.v2[e]
            if (!is.na(a) && !is.na(b) && a != b) {
                cut.l[a] <- cut.l[a] + 1
                cut.l[b] <- cut.l[b] + 1
            }
        }

        cond <- rep(NA_real_, length(lev)); names(cond) <- lev
        for (l in lev) {
            vol.a <- vol.l[l]
            vol.b <- vol.total - vol.a
            denom <- min(vol.a, vol.b)
            cond[l] <- if (denom > 0) cut.l[l] / denom else NA_real_
        }

        cond.by.label <- data.frame(
            cst = lev,
            vol = as.numeric(vol.l[lev]),
            cut = as.numeric(cut.l[lev]),
            conductance = as.numeric(cond[lev]),
            stringsAsFactors = FALSE
        )

        ## summarize conductance (weighted by vol is often sensible)
        w.vol <- cond.by.label$vol
        ok.c <- is.finite(cond.by.label$conductance) & w.vol > 0
        cond.summary <- list(
            conductance.median = stats::median(cond.by.label$conductance[ok.c], na.rm = TRUE),
            conductance.vol.weighted.mean = sum(cond.by.label$conductance[ok.c] * w.vol[ok.c]) / sum(w.vol[ok.c])
        )
    }

    ## permutation nulls (assortativity + homophily)
    perm.res <- NULL
    n.perm <- as.integer(n.perm)
    if (n.perm > 0L && sum(ok.v) >= 10L && sum(ok.e) >= 10L) {
        permute.labels <- function(lbl, blocks = NULL) {
            lbl2 <- lbl
            idx <- which(!is.na(lbl2))
            if (is.null(blocks)) {
                lbl2[idx] <- sample(lbl2[idx], replace = FALSE)
            } else {
                if (length(blocks) != length(lbl2)) stop("perm.blocks must have length vcount(g).")
                b <- blocks
                for (bb in unique(b[idx])) {
                    ii <- idx[b[idx] == bb]
                    if (length(ii) >= 2L) lbl2[ii] <- sample(lbl2[ii], replace = FALSE)
                }
            }
            lbl2
        }

        set.seed(seed)
        r.null <- rep(NA_real_, n.perm)
        h.null <- rep(NA_real_, n.perm)

        for (b in seq_len(n.perm)) {
            lbl.p <- permute.labels(labels, blocks = perm.blocks)

            ## homophily
            lab.u.p <- lbl.p[u]
            lab.v.p <- lbl.p[v]
            ok.e.p <- !is.na(lab.u.p) & !is.na(lab.v.p)
            w.ok.p <- w[ok.e.p]
            same.p <- (lab.u.p[ok.e.p] == lab.v.p[ok.e.p])
            h.null[b] <- if (sum(w.ok.p) > 0) sum(w.ok.p[same.p]) / sum(w.ok.p) else NA_real_

            ## assortativity
            ok.v.p <- !is.na(lbl.p)
            g.sub.p <- igraph::induced_subgraph(g, vids = which(ok.v.p))
            types.p <- as.integer(as.factor(lbl.p[ok.v.p]))
            r.null[b] <- igraph::assortativity_nominal(g.sub.p, types = types.p, directed = FALSE)
        }

        ## z-scores and (one-sided) p-values
        z.from.null <- function(x, x.null) {
            mu <- mean(x.null, na.rm = TRUE)
            sd0 <- stats::sd(x.null, na.rm = TRUE)
            if (!is.finite(sd0) || sd0 <= 0) return(list(z = NA_real_, mu = mu, sd = sd0))
            list(z = (x - mu) / sd0, mu = mu, sd = sd0)
        }
        p.upper <- function(x, x.null) {
            x.null <- x.null[is.finite(x.null)]
            if (length(x.null) < 10L) return(NA_real_)
            (1 + sum(x.null >= x)) / (1 + length(x.null))
        }

        perm.res <- list(
            n.perm = n.perm,
            homophily.null = h.null,
            assortativity.null = r.null,
            homophily.z = z.from.null(homophily, h.null),
            assortativity.z = z.from.null(r.assort, r.null),
            homophily.p = p.upper(homophily, h.null),
            assortativity.p = p.upper(r.assort, r.null)
        )
    }

    list(
        n.vertices = n,
        n.edges = m,
        homophily = homophily,
        assortativity = r.assort,
        conductance.by.label = cond.by.label,
        conductance.summary = cond.summary,
        permutation = perm.res
    )
}
