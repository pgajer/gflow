#' Test Feature–Distance Associations Within a Geodesic Disk Using SPT Sectors
#'
#' @description
#' SPT-sector (direction-aware) counterpart of \code{test.disk.feature.distance.association()}.
#' Given a disk (set of vertices) around \code{center.vertex}, this function builds a
#' shortest-path tree (SPT) rooted at \code{center.vertex} on the full graph and assigns each
#' disk vertex to a \emph{sector} defined as the first hop from \code{center.vertex} along
#' the chosen shortest path in the SPT.
#'
#' Associations between feature signal and distance are then computed \emph{within each sector},
#' providing a graph-intrinsic directional (anisotropic) sensitivity analysis.
#'
#' Supports \code{distance.transform} modes:
#' \itemize{
#'   \item \code{"raw"}: use raw distances.
#'   \item \code{"rank"}: use scaled rank distances within each sector.
#'   \item \code{"quantile.bin"}: perform distance-quantile bin analysis within each sector.
#' }
#'
#' @param X Numeric matrix: rows are vertices/samples, columns are features.
#' @param vertices Integer vector of disk vertex indices (1-based).
#' @param dists Optional numeric vector of distances for \code{vertices}. If named by vertex ID
#'   it will be aligned; otherwise must be same length as \code{vertices}. If NULL, distances are
#'   computed via Dijkstra.
#' @param adj.list Graph adjacency list; a list of integer vectors (neighbors) of length \code{nV}.
#' @param weight.list Graph weight list; a list of numeric vectors (edge weights) matching \code{adj.list}.
#'   If NULL, all edge weights are treated as 1.
#' @param center.vertex Integer; center vertex index (1-based), typically the extremum vertex (e.g., M1).
#' @param eps Non-negative threshold defining a carrier: \code{X > eps}. Default 0.
#' @param p.min Minimum prevalence within the disk to include a feature. Default 0.5.
#' @param logit.model One of \code{"glm"}, \code{"firth"}, \code{"bayesglm"} and \code{"all"} for modeling presence/absence binary response.
#' @param transform One of \code{"log"}, \code{"logit"}, \code{"arcsin_sqrt"} for the primary analysis.
#' @param pseudo.count Small positive value used for transforms and CLR. Default 1e-6.
#' @param min.carriers Minimum number of carriers required for carrier-only abundance tests. Default 20.
#' @param min.sector.size Minimum number of disk vertices in a sector for that sector to be tested. Default 20.
#' @param do.clr Logical; if TRUE (default) run CLR sensitivity analysis.
#' @param clr.over One of \code{"carriers"} (default) or \code{"all"} controlling CLR evaluation subset.
#' @param distance.transform One of \code{"raw"} (default), \code{"rank"}, \code{"quantile.bin"}.
#' @param n.bins Integer number of distance quantile bins (for \code{"quantile.bin"}). Default 10.
#' @param min.per.bin Integer minimum per-bin count (for \code{"quantile.bin"}). Default 5.
#' @param bin.trend.test One of \code{"spearman"} (default) or \code{"wls"} for bin-level trend tests.
#' @param with.presence.bins Logical; if TRUE (default), compute binned prevalence trend as well.
#' @param stop.when.targets.reached Logical; if TRUE (default), Dijkstra stops early once all disk vertices
#'   are finalized. If FALSE, computes full SPT.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{sector.results}: data.frame with one row per (feature, sector).
#'   \item \code{bin.results}: optional nested list of bin summaries/tests (only when \code{distance.transform="quantile.bin"}).
#'   \item \code{sector.map}: data.frame mapping disk vertices to sector ids and distances.
#'   \item \code{meta}: list of settings used.
#' }
#'
test.disk.feature.distance.association.spt <- function(X,
                                                       vertices,
                                                       dists = NULL,
                                                       adj.list,
                                                       weight.list = NULL,
                                                       center.vertex,
                                                       eps = 0,
                                                       p.min = 0.5,
                                                       logit.model = c("glm","firth","bayesglm","all"),
                                                       transform = c("log", "logit", "arcsin_sqrt"),
                                                       pseudo.count = 1e-6,
                                                       min.carriers = 20L,
                                                       min.sector.size = 20L,
                                                       do.clr = TRUE,
                                                       clr.over = c("carriers", "all"),
                                                       distance.transform = c("raw", "rank", "quantile.bin"),
                                                       n.bins = 10L,
                                                       min.per.bin = 5L,
                                                       bin.trend.test = c("spearman", "wls"),
                                                       with.presence.bins = TRUE,
                                                       stop.when.targets.reached = TRUE) {

    ## ----------------------------
    ## Input validation
    ## ----------------------------
    transform <- match.arg(transform)
    clr.over <- match.arg(clr.over)
    distance.transform <- match.arg(distance.transform)
    bin.trend.test <- match.arg(bin.trend.test)
    logit.model <- match.arg(logit.model)

    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix (or coercible data.frame).")
    if (is.null(colnames(X))) stop("X must have column names.")
    nV <- nrow(X)

    if (!is.list(adj.list) || length(adj.list) != nV) {
        stop("adj.list must be a list of length nrow(X).")
    }
    if (is.null(weight.list)) {
        weight.list <- vector("list", length(adj.list))
        for (i in seq_len(nV)) {
            weight.list[[i]] <- rep(1, length(adj.list[[i]]))
        }
    }
    if (!is.list(weight.list) || length(weight.list) != nV) {
        stop("weight.list must be NULL or a list of length nrow(X).")
    }

    center.vertex <- as.integer(center.vertex)
    if (length(center.vertex) != 1L || is.na(center.vertex) || center.vertex < 1L || center.vertex > nV) {
        stop("center.vertex must be a single integer in [1, nrow(X)].")
    }

    vertices <- as.integer(vertices)
    vertices <- vertices[!is.na(vertices)]
    vertices <- vertices[vertices >= 1L & vertices <= nV]
    vertices <- unique(vertices)
    if (length(vertices) < 2L) stop("vertices must contain at least 2 valid indices.")
    if (!(center.vertex %in% vertices)) {
        ## Defensive: include the center in the disk set
        vertices <- sort(unique(c(vertices, center.vertex)))
    }

    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) stop("eps must be a single non-negative numeric.")
    if (!is.numeric(p.min) || length(p.min) != 1L || is.na(p.min) || p.min < 0 || p.min > 1) stop("p.min must be in [0, 1].")
    if (!is.numeric(pseudo.count) || length(pseudo.count) != 1L || is.na(pseudo.count) || pseudo.count <= 0) {
        stop("pseudo.count must be a single positive numeric.")
    }
    min.carriers <- as.integer(min.carriers)
    min.sector.size <- as.integer(min.sector.size)

    n.bins <- as.integer(n.bins)
    min.per.bin <- as.integer(min.per.bin)
    if (distance.transform == "quantile.bin") {
        if (n.bins < 2L) stop("n.bins must be >= 2 for quantile.bin.")
        if (min.per.bin < 1L) stop("min.per.bin must be >= 1 for quantile.bin.")
    }

    ## ----------------------------
    ## Helpers: transforms + CLR
    ## ----------------------------
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

    compute.clr.disk <- function(X.sub) {
        log.X <- log(X.sub + pseudo.count)
        row.center <- rowMeans(log.X)
        sweep(log.X, 1, row.center, FUN = "-")
    }

    distance.rank <- function(dd) {
        r <- rank(dd, ties.method = "average")
        r / (length(r) + 1)
    }

    ## ----------------------------
    ## Dijkstra with parent pointers (SPT)
    ## ----------------------------
    dijkstra.spt <- function(adj.list, weight.list, source, targets = NULL, stop.early = TRUE) {
        n <- length(adj.list)

        dist <- rep(Inf, n)
        parent <- rep(NA_integer_, n)
        visited <- rep(FALSE, n)

        ## ---- binary heap (min-heap on dist) ----
        heap.v <- integer(n)
        heap.k <- rep(Inf, n)
        heap.pos <- rep(0L, n)
        heap.size <- 0L

        heap.swap <- function(i, j) {
            vi <- heap.v[i]; vj <- heap.v[j]
            ki <- heap.k[i]; kj <- heap.k[j]
            heap.v[i] <<- vj; heap.v[j] <<- vi
            heap.k[i] <<- kj; heap.k[j] <<- ki
            heap.pos[heap.v[i]] <<- i
            heap.pos[heap.v[j]] <<- j
        }

        heap.sift.up <- function(i) {
            while (i > 1L) {
                p <- i %/% 2L
                if (heap.k[p] <= heap.k[i]) break
                heap.swap(p, i)
                i <- p
            }
        }

        heap.sift.down <- function(i) {
            while (TRUE) {
                l <- 2L * i
                r <- l + 1L
                if (l > heap.size) break
                m <- l
                if (r <= heap.size && heap.k[r] < heap.k[l]) m <- r
                if (heap.k[i] <= heap.k[m]) break
                heap.swap(i, m)
                i <- m
            }
        }

        heap.push.or.decrease <- function(v, key) {
            p <- heap.pos[v]
            if (p == 0L) {
                heap.size <<- heap.size + 1L
                heap.v[heap.size] <<- v
                heap.k[heap.size] <<- key
                heap.pos[v] <<- heap.size
                heap.sift.up(heap.size)
            } else {
                if (key < heap.k[p]) {
                    heap.k[p] <<- key
                    heap.sift.up(p)
                }
            }
        }

        heap.pop.min <- function() {
            if (heap.size == 0L) return(c(NA_integer_, NA_real_))
            v.min <- heap.v[1]
            k.min <- heap.k[1]
            heap.pos[v.min] <<- 0L
            if (heap.size == 1L) {
                heap.size <<- 0L
                return(c(v.min, k.min))
            }
            heap.v[1] <<- heap.v[heap.size]
            heap.k[1] <<- heap.k[heap.size]
            heap.pos[heap.v[1]] <<- 1L
            heap.size <<- heap.size - 1L
            heap.sift.down(1L)
            c(v.min, k.min)
        }

        dist[source] <- 0
        parent[source] <- source
        heap.push.or.decrease(source, 0)

        targets.set <- NULL
        targets.remaining <- 0L
        if (!is.null(targets)) {
            targets.set <- rep(FALSE, n)
            targets.set[targets] <- TRUE
            targets.remaining <- sum(targets.set)
        }

        while (heap.size > 0L) {
            top <- heap.pop.min()
            v <- as.integer(top[1])
            if (is.na(v)) break
            if (visited[v]) next
            visited[v] <- TRUE

            if (!is.null(targets.set) && targets.set[v]) {
                targets.remaining <- targets.remaining - 1L
                if (isTRUE(stop.early) && targets.remaining <= 0L) break
            }

            nb <- adj.list[[v]]
            w <- weight.list[[v]]
            if (length(nb) != length(w)) stop("Mismatch adj.list and weight.list at vertex ", v, ".")

            for (i in seq_along(nb)) {
                u <- as.integer(nb[i])
                if (u < 1L || u > n) next
                if (visited[u]) next
                alt <- dist[v] + w[i]
                if (alt < dist[u]) {
                    dist[u] <- alt
                    parent[u] <- v
                    heap.push.or.decrease(u, alt)
                }
            }
        }

        list(dist = dist, parent = parent, visited = visited)
    }

    ## ----------------------------
    ## Obtain distances + parent pointers
    ## ----------------------------
    dist.full <- NULL
    parent.full <- NULL

    if (is.null(dists)) {
        dj <- dijkstra.spt(adj.list, weight.list, source = center.vertex,
                           targets = vertices, stop.early = stop.when.targets.reached)
        dist.full <- dj$dist
        parent.full <- dj$parent
    } else {
        ## Use provided distances for the disk; still compute parents for sector assignment
        dj <- dijkstra.spt(adj.list, weight.list, source = center.vertex,
                           targets = vertices, stop.early = stop.when.targets.reached)
        dist.full <- dj$dist
        parent.full <- dj$parent

        ## Replace dist on disk with supplied values (aligned)
        if (!is.null(names(dists))) {
            d.sub <- as.numeric(dists[as.character(vertices)])
        } else {
            if (length(dists) != length(vertices)) stop("If dists is unnamed, it must align with vertices (same length).")
            d.sub <- as.numeric(dists)
        }
        ok <- is.finite(d.sub)
        if (!all(ok)) {
            ## Keep only finite-distance vertices
            vertices <- vertices[ok]
            d.sub <- d.sub[ok]
        }
        dist.full[vertices] <- d.sub
    }

    d.disk <- dist.full[vertices]
    ok <- is.finite(d.disk)
    vertices <- vertices[ok]
    d.disk <- d.disk[ok]
    if (length(vertices) < 2L) stop("Too few disk vertices have finite distances from center.vertex.")

    ## ----------------------------
    ## Sector assignment: first hop on SPT path
    ## ----------------------------
    sector.of.vertex <- rep(NA_integer_, nV)
    sector.of.vertex[center.vertex] <- center.vertex

    get.first.hop <- function(v) {
        if (v == center.vertex) return(center.vertex)
        u <- v
        p <- parent.full[u]
        if (is.na(p)) return(NA_integer_)
        ## Walk up until parent is center or we fail
        while (!is.na(p) && p != center.vertex && p != u) {
            u <- p
            p <- parent.full[u]
        }
        if (!is.na(p) && p == center.vertex) return(u)
        NA_integer_
    }

    for (v in vertices) {
        sector.of.vertex[v] <- get.first.hop(v)
    }

    sector.disk <- sector.of.vertex[vertices]
    ## Drop vertices without a valid sector (disconnected / parent missing)
    ok2 <- !is.na(sector.disk)
    vertices <- vertices[ok2]
    d.disk <- d.disk[ok2]
    sector.disk <- sector.disk[ok2]

    ## Sector sizes in disk
    sector.tab <- sort(table(sector.disk), decreasing = TRUE)
    keep.sectors <- as.integer(names(sector.tab)[sector.tab >= min.sector.size])

    if (length(keep.sectors) == 0L) {
        stop("No sectors meet min.sector.size = ", min.sector.size, " within the supplied disk.")
    }

    ## Keep only vertices in retained sectors
    keep.v <- sector.disk %in% keep.sectors
    vertices <- vertices[keep.v]
    d.disk <- d.disk[keep.v]
    sector.disk <- sector.disk[keep.v]

    ## Map for return
    sector.map <- data.frame(
        vertex = as.integer(vertices),
        dist = as.numeric(d.disk),
        sector = as.integer(sector.disk),
        stringsAsFactors = FALSE
    )

    ## ----------------------------
    ## Prevalence filter over the (retained) disk vertices
    ## ----------------------------
    prev.df <- vertices.feature.carriers(vertices = vertices, X = X, eps = eps, match.by = "index")
    keep.feature <- prev.df$p >= p.min
    feats <- rownames(prev.df)[keep.feature]
    if (length(feats) == 0L) {
        return(list(sector.results = data.frame(),
                    bin.results = NULL,
                    sector.map = sector.map,
                    meta = list(distance.transform = distance.transform)))
    }

    ## CLR (computed once for the retained disk)
    clr.mat <- NULL
    if (isTRUE(do.clr)) {
        clr.mat <- compute.clr.disk(X[vertices, , drop = FALSE])
        ## rows correspond to 'vertices' order in sector.map
    }

    ## ----------------------------
    ## Quantile-bin helper (within sector)
    ## ----------------------------
    quantile.bin.analysis <- function(x, dd, carrier) {

        idx.ab <- if (clr.over == "all") rep(TRUE, length(dd)) else carrier
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

        if (isTRUE(with.presence.bins)) {
            bin.all <- cut(dd, breaks = brks, include.lowest = TRUE)
            n.all <- tapply(dd, bin.all, length)
            prev <- tapply(carrier, bin.all, mean)
            mid.all <- tapply(dd, bin.all, function(z) mean(range(z)))
            bins.df$n.total <- as.integer(n.all[bins.df$bin])
            bins.df$presence.p <- as.numeric(prev[bins.df$bin])
            bins.df$d.mid.all <- as.numeric(mid.all[bins.df$bin])
        }

        keep.bin <- bins.df$n.abundance >= min.per.bin
        if (isTRUE(with.presence.bins)) keep.bin <- keep.bin & (bins.df$n.total >= min.per.bin)
        bins.df <- bins.df[keep.bin, , drop = FALSE]
        if (nrow(bins.df) < 3L) return(NULL)

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
                fit <- stats::lm(presence.p ~ d.mid.all, data = bins.df, weights = bins.df$n.total)
                sm <- summary(fit)$coefficients["d.mid.all", ]
                tests$presence <- list(method = "wls", slope = unname(sm["Estimate"]), p.value = unname(sm["Pr(>|t|)"]))
            }
        }

        list(bins = bins.df, tests = tests, breaks = brks)
    }

    fit.presence.logit.models <- function(carrier, d, logit.model = c("glm", "firth", "bayesglm", "all")) {

        logit.model <- match.arg(logit.model)

        ## Convenience: return NA if degenerate
        if (!(any(carrier == 0L) && any(carrier == 1L))) {
            out <- list()
            if (logit.model == "all") {
                out$glm <- c(slope = NA_real_, p = NA_real_)
                out$firth <- c(slope = NA_real_, p = NA_real_)
                out$bayesglm <- c(slope = NA_real_, p = NA_real_)
            } else {
                out[[logit.model]] <- c(slope = NA_real_, p = NA_real_)
            }
            return(out)
        }

        fit.glm <- function() {
            tryCatch({
                fit <- suppressWarnings(stats::glm(carrier ~ d, family = stats::binomial()))
                sm <- summary(fit)
                slope <- unname(stats::coef(fit)["d"])
                p <- unname(sm$coefficients["d", "Pr(>|z|)"])
                c(slope = slope, p = p)
            }, error = function(e) c(slope = NA_real_, p = NA_real_))
        }

        ## Firth / separation-robust logistic:
        ## Prefer brglm2; fall back to logistf; then fall back to glm
        fit.firth <- function() {

            ## brglm2 path
            if (requireNamespace("brglm2", quietly = TRUE)) {

                out <- tryCatch({

                    ## design matrix for carrier ~ d (intercept + d)
                    Xmat <- stats::model.matrix(~ d)

                    fit <- brglm2::brglmFit(
                                       x = Xmat,
                                       y = carrier,
                                       family = stats::binomial("logit"),
                                       type = "AS_mean"
                                   )

                    ## coefficient for d is the second column of Xmat
                    coef.vec <- stats::coef(fit)
                    slope <- unname(coef.vec[2])

                    ## Wald SE via vcov if available
                    V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
                    if (is.null(V) || any(!is.finite(diag(V)))) {
                        return(c(slope = slope, p = NA_real_))
                    }

                    se <- sqrt(unname(diag(V)[2]))
                    if (!is.finite(se) || se <= 0) {
                        return(c(slope = slope, p = NA_real_))
                    }

                    z <- slope / se
                    p <- 2 * stats::pnorm(-abs(z))

                    c(slope = slope, p = p)

                }, error = function(e) NULL)

                if (!is.null(out)) return(out)
            }

            ## logistf fallback (Firth penalized likelihood)
            if (requireNamespace("logistf", quietly = TRUE)) {
                out <- tryCatch({
                    fit <- logistf::logistf(carrier ~ d)
                    slope <- unname(stats::coef(fit)["d"])
                    ## logistf stores p-values in fit$prob (named)
                    p <- unname(fit$prob["d"])
                    c(slope = slope, p = p)
                }, error = function(e) NULL)

                if (!is.null(out)) return(out)
            }

            ## Final fallback: standard glm (better than all-NA)
            fit.glm()
        }

        ## Weakly-informative-prior logistic (regularized)
        fit.bayesglm <- function() {
            if (!requireNamespace("arm", quietly = TRUE)) {
                return(c(slope = NA_real_, p = NA_real_))
            }
            tryCatch({
                fit <- suppressWarnings(
                    arm::bayesglm(carrier ~ d, family = stats::binomial(),
                                  prior.mean = 0, prior.scale = 2.5, prior.df = Inf)
                )
                cf <- summary(fit)$coefficients
                slope <- unname(cf["d", "Estimate"])
                se <- unname(cf["d", "Std. Error"])
                z <- slope / se
                p <- 2 * stats::pnorm(-abs(z))
                c(slope = slope, p = p)
            }, error = function(e) c(slope = NA_real_, p = NA_real_))
        }

        if (logit.model == "glm") return(list(glm = fit.glm()))
        if (logit.model == "firth") return(list(firth = fit.firth()))
        if (logit.model == "bayesglm") return(list(bayesglm = fit.bayesglm()))

        ## all
        list(
            glm = fit.glm(),
            firth = fit.firth(),
            bayesglm = fit.bayesglm()
        )
    }

    ## ----------------------------
    ## Main loop over (feature, sector)
    ## ----------------------------
    res.rows <- list()
    bin.results <- NULL
    if (distance.transform == "quantile.bin") {
        bin.results <- list()
    }

    ## Precompute per-sector indices (in sector.map order)
    sector.levels <- sort(unique(sector.map$sector))
    sector.idx.list <- lapply(sector.levels, function(s) which(sector.map$sector == s))
    names(sector.idx.list) <- as.character(sector.levels)

    for (f in feats) {
        ## feature vector on retained disk (aligned with sector.map rows)
        x.disk <- X[vertices, f]
        carrier.disk <- x.disk > eps

        ## transformed (non-CLR) abundance on retained disk
        y.disk <- transform.vec(x.disk)

        ## CLR vector on retained disk
        y.clr.disk <- NULL
        if (isTRUE(do.clr)) {
            y.clr.disk <- clr.mat[, f]
        }

        for (s in sector.levels) {

            idx <- sector.idx.list[[as.character(s)]]
            if (length(idx) < min.sector.size) next

            d.s <- sector.map$dist[idx]
            if (distance.transform == "rank") {
                d.s.used <- distance.rank(d.s)
            } else {
                d.s.used <- d.s
            }

            carrier.s <- carrier.disk[idx]
            n.all <- length(idx)
            n.car <- sum(carrier.s)

            ## Presence model(s)
            ## Primary outputs always exist
            pres.slope <- NA_real_
            pres.p <- NA_real_

            ## Per-model outputs ONLY when logit.model == "all"
            if (logit.model == "all") {
                pres.slope.glm <- NA_real_
                pres.p.glm <- NA_real_
                pres.slope.firth <- NA_real_
                pres.p.firth <- NA_real_
                pres.slope.bayesglm <- NA_real_
                pres.p.bayesglm <- NA_real_
            }

            ## IMPORTANT: d.s.used is the predictor used for the presence model in this sector
            ## carrier.s is the 0/1 response in this sector

            pres.fit <- fit.presence.logit.models(
                carrier = as.integer(carrier.s),
                d = d.s.used,
                logit.model = logit.model
            )

            if (logit.model == "all") {

                pres.slope.glm <- as.numeric(pres.fit$glm["slope"])
                pres.p.glm <- as.numeric(pres.fit$glm["p"])

                pres.slope.firth <- as.numeric(pres.fit$firth["slope"])
                pres.p.firth <- as.numeric(pres.fit$firth["p"])

                pres.slope.bayesglm <- as.numeric(pres.fit$bayesglm["slope"])
                pres.p.bayesglm <- as.numeric(pres.fit$bayesglm["p"])

                ## Choose a primary presence p-value/slope for backwards compatibility:
                ## prefer firth if available (handles separation), else glm, else bayesglm.
                if (is.finite(pres.p.firth)) {
                    pres.slope <- pres.slope.firth
                    pres.p <- pres.p.firth
                } else if (is.finite(pres.p.glm)) {
                    pres.slope <- pres.slope.glm
                    pres.p <- pres.p.glm
                } else {
                    pres.slope <- pres.slope.bayesglm
                    pres.p <- pres.p.bayesglm
                }

            } else if (logit.model == "glm") {
                pres.slope <- as.numeric(pres.fit$glm["slope"])
                pres.p <- as.numeric(pres.fit$glm["p"])
            } else if (logit.model == "firth") {
                pres.slope <- as.numeric(pres.fit$firth["slope"])
                pres.p <- as.numeric(pres.fit$firth["p"])
            } else if (logit.model == "bayesglm") {
                pres.slope <- as.numeric(pres.fit$bayesglm["slope"])
                pres.p <- as.numeric(pres.fit$bayesglm["p"])
            }

            ## Abundance among carriers
            rho <- NA_real_
            rho.p <- NA_real_
            lm.slope <- NA_real_
            lm.p <- NA_real_

            if (n.car >= min.carriers) {
                yy <- y.disk[idx][carrier.s]
                dd <- d.s.used[carrier.s]

                sp <- suppressWarnings(stats::cor.test(yy, dd, method = "spearman", exact = FALSE))
                rho <- unname(sp$estimate)
                rho.p <- sp$p.value

                fit.lm <- stats::lm(yy ~ dd)
                lm.slope <- unname(stats::coef(fit.lm)["dd"])
                lm.p <- summary(fit.lm)$coefficients["dd", "Pr(>|t|)"]
            }

            ## CLR sensitivity
            clr.rho <- NA_real_
            clr.rho.p <- NA_real_
            clr.lm.slope <- NA_real_
            clr.lm.p <- NA_real_

            if (isTRUE(do.clr)) {
                idx.use <- if (clr.over == "all") rep(TRUE, length(idx)) else carrier.s
                if (sum(idx.use) >= 2L) {
                    yy.clr <- y.clr.disk[idx][idx.use]
                    dd.clr <- d.s.used[idx.use]

                    spc <- suppressWarnings(stats::cor.test(yy.clr, dd.clr, method = "spearman", exact = FALSE))
                    clr.rho <- unname(spc$estimate)
                    clr.rho.p <- spc$p.value

                    fit.clr <- stats::lm(yy.clr ~ dd.clr)
                    clr.lm.slope <- unname(stats::coef(fit.clr)["dd.clr"])
                    clr.lm.p <- summary(fit.clr)$coefficients["dd.clr", "Pr(>|t|)"]
                }
            }

            ## Quantile-bin mode: store bin summaries + bin-level tests
            if (distance.transform == "quantile.bin") {
                key <- paste0(f, "||", s)
                bin.entry <- list()

                ## Non-CLR bins (use transformed y.disk; carriers define subset)
                bin.entry$logit <- quantile.bin.analysis(
                    x = y.disk[idx],
                    dd = d.s,
                    carrier = carrier.s
                )

                if (isTRUE(do.clr)) {
                    bin.entry$clr <- quantile.bin.analysis(
                        x = y.clr.disk[idx],
                        dd = d.s,
                        carrier = carrier.s
                    )
                }

                bin.results[[key]] <- bin.entry
            }

            ## Collect row
            row <- data.frame(
                feature = f,
                sector = as.integer(s),
                n.sector = as.integer(n.all),
                n.carriers = as.integer(n.car),
                p.disk = as.numeric(prev.df[f, "p"]),
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

            if (logit.model == "all") {
                row$presence.slope.glm <- pres.slope.glm
                row$presence.p.glm <- pres.p.glm
                row$presence.slope.firth <- pres.slope.firth
                row$presence.p.firth <- pres.p.firth
                row$presence.slope.bayesglm <- pres.slope.bayesglm
                row$presence.p.bayesglm <- pres.p.bayesglm
            }

            res.rows[[length(res.rows) + 1L]] <- row
        }
    }

    sector.df <- do.call(rbind, res.rows)
    if (is.null(sector.df) || nrow(sector.df) == 0L) {
        return(list(sector.results = data.frame(),
                    bin.results = bin.results,
                    sector.map = sector.map,
                    meta = list(distance.transform = distance.transform)))
    }

    ## ----------------------------
    ## Multiple-testing correction across all (feature, sector)
    ## ----------------------------
    sector.df$presence.p.adj <- stats::p.adjust(sector.df$presence.p, method = "BH")
    sector.df$spearman.p.adj <- stats::p.adjust(sector.df$spearman.p, method = "BH")
    sector.df$lm.p.adj <- stats::p.adjust(sector.df$lm.p, method = "BH")

    if (isTRUE(do.clr)) {
        sector.df$clr.spearman.p.adj <- stats::p.adjust(sector.df$clr.spearman.p, method = "BH")
        sector.df$clr.lm.p.adj <- stats::p.adjust(sector.df$clr.lm.p, method = "BH")
    }

    if ("presence.p.glm" %in% names(sector.df)) {
        sector.df$presence.p.glm.adj <- stats::p.adjust(sector.df$presence.p.glm, method = "BH")
    }
    if ("presence.p.firth" %in% names(sector.df)) {
        sector.df$presence.p.firth.adj <- stats::p.adjust(sector.df$presence.p.firth, method = "BH")
    }
    if ("presence.p.bayesglm" %in% names(sector.df)) {
        sector.df$presence.p.bayesglm.adj <- stats::p.adjust(sector.df$presence.p.bayesglm, method = "BH")
    }

    ## Drop columns that are all NA (but keep essential keys)
    keep.col <- vapply(sector.df, function(z) any(!is.na(z)), logical(1))

    ## Always keep these even if all NA (defensive; usually they are not all NA)
    keep.always <- names(sector.df) %in% c("feature", "sector", "n.sector", "n.carriers", "p.disk")
    keep.col <- keep.col | keep.always

    sector.df <- sector.df[, keep.col, drop = FALSE]
    rownames(sector.df) <- NULL

    ## Bin-level p.adjust for quantile.bin (stored inside bin.results)
    if (distance.transform == "quantile.bin" && length(bin.results) > 0L) {

        ## Extract p-values across keys
        keys <- names(bin.results)

        get.p <- function(slot, test.name) {
            p <- rep(NA_real_, length(keys))
            for (i in seq_along(keys)) {
                obj <- bin.results[[keys[i]]][[slot]]
                if (!is.null(obj) && !is.null(obj$tests[[test.name]])) {
                    p[i] <- obj$tests[[test.name]]$p.value
                }
            }
            p
        }

        p.logit.ab <- get.p("logit", "abundance")
        p.clr.ab <- if (isTRUE(do.clr)) get.p("clr", "abundance") else NULL

        padj.logit.ab <- stats::p.adjust(p.logit.ab, method = "BH")
        padj.clr.ab <- if (!is.null(p.clr.ab)) stats::p.adjust(p.clr.ab, method = "BH") else NULL

        ## Store back
        for (i in seq_along(keys)) {
            k <- keys[i]
            if (!is.null(bin.results[[k]]$logit) && !is.null(bin.results[[k]]$logit$tests$abundance)) {
                bin.results[[k]]$logit$tests$abundance$p.adj <- padj.logit.ab[i]
            }
            if (isTRUE(do.clr) && !is.null(bin.results[[k]]$clr) && !is.null(bin.results[[k]]$clr$tests$abundance)) {
                bin.results[[k]]$clr$tests$abundance$p.adj <- padj.clr.ab[i]
            }
        }
    }

    ## ----------------------------
    ## Return
    ## ----------------------------
    list(
        sector.results = sector.df,
        bin.results = bin.results,
        sector.map = sector.map,
        meta = list(
            center.vertex = center.vertex,
            distance.transform = distance.transform,
            min.sector.size = min.sector.size,
            n.bins = if (distance.transform == "quantile.bin") n.bins else NA_integer_,
            min.per.bin = if (distance.transform == "quantile.bin") min.per.bin else NA_integer_,
            bin.trend.test = if (distance.transform == "quantile.bin") bin.trend.test else NA_character_,
            clr.over = if (isTRUE(do.clr)) clr.over else NA_character_
        )
    )
}
