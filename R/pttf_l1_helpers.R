# Private L1 helper subset retained for fit.pttf.trend.filtering().
# Public SSRHE smoother APIs moved to geosmooth.

.ssrhe.operator.timing.row <- function(phase, elapsed.sec) {
    data.frame(phase = phase, elapsed.sec = as.numeric(elapsed.sec), stringsAsFactors = FALSE)
}

.validate.ssrhe.X <- function(X) {
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    storage.mode(X) <- "double"
    if (!is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L) {
        stop("X must be a numeric matrix with at least two rows and one column.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("X must contain only finite values.", call. = FALSE)
    }
    X
}

.validate.ssrhe.positive.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a positive integer scalar.", name), call. = FALSE)
    }
    as.integer(round(x))
}

.validate.ssrhe.derivative.order <- function(x) {
    x <- .validate.ssrhe.positive.integer(x, "derivative.order")
    if (!x %in% c(2L, 3L)) {
        stop("derivative.order must be either 2L or 3L.", call. = FALSE)
    }
    x
}

.validate.ssrhe.numeric.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) {
        stop(sprintf("%s must be a finite numeric scalar.", name), call. = FALSE)
    }
    as.double(x)
}

.validate.ssrhe.nn.index <- function(nn.index, n, k) {
    if (is.null(nn.index)) {
        return(NULL)
    }
    if (!is.matrix(nn.index)) {
        nn.index <- as.matrix(nn.index)
    }
    storage.mode(nn.index) <- "integer"
    if (!identical(dim(nn.index), c(n, k))) {
        stop("nn.index must be an nrow(X) by k integer matrix.", call. = FALSE)
    }
    if (any(is.na(nn.index)) || any(nn.index < 1L) || any(nn.index > n)) {
        stop("nn.index must contain integer indices in 1:nrow(X).", call. = FALSE)
    }
    has.center <- vapply(seq_len(n), function(i) any(nn.index[i, ] == i), logical(1))
    if (!all(has.center)) {
        stop("Each row of nn.index must contain its center vertex.", call. = FALSE)
    }
    nn.index
}

.validate.ssrhe.support.index <- function(support.index, n) {
    if (is.null(support.index)) {
        return(NULL)
    }
    if (is.matrix(support.index)) {
        support.index <- lapply(seq_len(nrow(support.index)), function(i) {
            support.index[i, !is.na(support.index[i, ])]
        })
    }
    if (!is.list(support.index) || length(support.index) != n) {
        stop("support.index must be a list with length nrow(X).", call. = FALSE)
    }
    out <- vector("list", n)
    for (i in seq_len(n)) {
        ids <- as.integer(support.index[[i]])
        if (length(ids) < 2L || any(is.na(ids)) || any(ids < 1L) || any(ids > n)) {
            stop("Each support.index element must contain at least two valid vertex indices.",
                call. = FALSE)
        }
        ids <- unique(ids)
        if (!i %in% ids) {
            stop("Each support.index element must contain its center vertex.", call. = FALSE)
        }
        out[[i]] <- ids
    }
    out
}

.ssrhe.design.ncol <- function(m, derivative.order = 2L) {
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    q2 <- m * (m + 1L)/2L
    if (derivative.order == 2L) {
        q2 + m
    }
    else {
        q3 <- m * (m + 1L) * (m + 2L)/6L
        q3 + q2 + m
    }
}

.ssrhe.default.min.support <- function(tangent.dim, tangent.dim.rule, ambient.dim, support.buffer,
    n, derivative.order = 2L) {
    support.buffer <- .validate.ssrhe.nonnegative.integer(support.buffer, "support.buffer")
    if (!is.null(tangent.dim)) {
        m <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
    }
    else if (identical(tangent.dim.rule, "eigen.cumulative")) {
        m <- ambient.dim
    }
    else {
        stop("tangent.dim is required to choose the default min.support.", call. = FALSE)
    }
    min(n, .ssrhe.design.ncol(m, derivative.order = derivative.order) + support.buffer)
}

.validate.ssrhe.support.grid <- function(support.grid, n, tangent.dim, derivative.order, support.buffer,
    max.candidates) {
    if (is.null(support.grid)) {
        return(ssrhe.support.grid(n = n, tangent.dim = tangent.dim, derivative.order = derivative.order,
            support.buffer = support.buffer, max.candidates = max.candidates))
    }
    if (!is.data.frame(support.grid)) {
        stop("support.grid must be a data frame.", call. = FALSE)
    }
    required <- c("adaptive.k.scale", "min.support")
    if (!all(required %in% names(support.grid))) {
        stop("support.grid must contain adaptive.k.scale and min.support columns.", call. = FALSE)
    }
    grid <- support.grid[, unique(c(required, intersect("max.support", names(support.grid)))),
        drop = FALSE]
    if (!"max.support" %in% names(grid)) {
        grid$max.support <- NA_integer_
    }
    grid$adaptive.k.scale <- as.integer(round(grid$adaptive.k.scale))
    grid$min.support <- as.integer(round(grid$min.support))
    grid$max.support <- ifelse(is.na(grid$max.support), NA_integer_, as.integer(round(grid$max.support)))
    bad <- is.na(grid$adaptive.k.scale) | grid$adaptive.k.scale < 1L | grid$adaptive.k.scale >=
        n | is.na(grid$min.support) | grid$min.support < 2L | grid$min.support > n | (!is.na(grid$max.support) &
        (grid$max.support < grid$min.support | grid$max.support > n))
    if (any(bad)) {
        stop("support.grid contains invalid adaptive.k.scale/min.support/max.support values.",
            call. = FALSE)
    }
    grid <- unique(grid)
    if (nrow(grid) > max.candidates) {
        grid <- grid[seq_len(max.candidates), , drop = FALSE]
    }
    rownames(grid) <- NULL
    grid
}

.ssrhe.support.diagnostics <- function(operator) {
    support.size <- operator$neighborhoods$support.size
    topup <- operator$neighborhoods$n.topup
    truncated <- operator$neighborhoods$n.truncated
    data.frame(support.size.min = min(support.size), support.size.median = stats::median(support.size),
        support.size.max = max(support.size), support.topup.median = stats::median(topup, na.rm = TRUE),
        support.topup.max = max(topup, na.rm = TRUE), support.truncated.max = max(truncated,
            na.rm = TRUE), operator.rows = operator$A.triplet$dim[[1L]], operator.cols = operator$A.triplet$dim[[2L]],
        stringsAsFactors = FALSE)
}

.ssrhe.local.geometry.diagnostics <- function(X, support.index, diagnostics) {
    n <- nrow(X)
    if (!is.list(support.index) || length(support.index) != n) {
        stop("Internal error: invalid SSRHE support.index.", call. = FALSE)
    }
    out <- diagnostics
    support.size <- out$k
    if (is.null(support.size)) {
        support.size <- lengths(support.index)
    }
    out$support.size <- support.size
    out$design.rank.deficiency <- pmax(0L, out$design.ncol - out$design.rank)
    out$pca.variance.explained <- out$local.variance.ratio
    out$pca.discarded.variance.ratio <- pmax(0, 1 - out$local.variance.ratio)
    chart.distortion <- rep(NA_real_, n)
    boundary.asymmetry <- rep(NA_real_, n)
    curvature.bias.proxy <- out$pca.discarded.variance.ratio
    for (center in seq_len(n)) {
        ids <- as.integer(support.index[[center]])
        if (!length(ids) || !center %in% ids)
            next
        local <- X[ids, , drop = FALSE]
        centered <- sweep(local, 2L, colMeans(local), "-")
        svd.local <- tryCatch(svd(centered), error = function(e) NULL)
        if (is.null(svd.local) || !length(svd.local$d))
            next
        m <- as.integer(out$tangent.dim[center])
        m <- max(1L, min(m, ncol(svd.local$v), length(svd.local$d)))
        coords <- centered %*% svd.local$v[, seq_len(m), drop = FALSE]
        base.local <- match(center, ids)
        coords <- sweep(coords, 2L, coords[base.local, ], "-")
        radius <- sqrt(max(rowSums(coords^2), na.rm = TRUE))
        centroid.norm <- sqrt(sum(colMeans(coords)^2))
        boundary.asymmetry[center] <- if (is.finite(radius) && radius > 0) {
            centroid.norm/radius
        }
        else {
            0
        }
        if (length(ids) >= 3L) {
            d.ambient <- stats::dist(local)
            d.chart <- stats::dist(coords)
            d.ambient <- as.numeric(d.ambient)
            d.chart <- as.numeric(d.chart)
            ok <- is.finite(d.ambient) & is.finite(d.chart) & d.ambient > 0
            denom <- mean(d.ambient[ok])
            chart.distortion[center] <- if (any(ok) && is.finite(denom) && denom > 0) {
                sqrt(mean((d.chart[ok] - d.ambient[ok])^2))/denom
            }
            else {
                0
            }
        }
        else {
            chart.distortion[center] <- 0
        }
    }
    out$chart.distortion <- chart.distortion
    out$boundary.asymmetry <- boundary.asymmetry
    out$curvature.bias.proxy <- curvature.bias.proxy
    out
}

.validate.ssrhe.nonnegative.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a nonnegative integer scalar.", name), call. = FALSE)
    }
    as.integer(round(x))
}

.prepare.ssrhe.neighborhood <- function(X, k, nn.index, neighborhood.type, support.index, adaptive.k.scale,
    radius.rule, radius.factor, min.support, max.support, support.buffer, support.topup, tangent.dim,
    tangent.dim.rule, derivative.order = 2L, return.timing = FALSE) {
    return.timing <- isTRUE(return.timing)
    timing.start <- proc.time()[["elapsed"]]
    timing.phase.start <- timing.start
    timing.rows <- list()
    add.timing <- function(subphase) {
        if (return.timing) {
            timing.rows[[length(timing.rows) + 1L]] <<- .ssrhe.operator.timing.row(phase = paste0("neighborhood.",
                subphase), elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start)
            timing.phase.start <<- proc.time()[["elapsed"]]
        }
        invisible(NULL)
    }
    finalize.timing <- function() {
        if (!return.timing) {
            return(NULL)
        }
        timing <- do.call(rbind, timing.rows)
        if (is.null(timing)) {
            timing <- data.frame(phase = character(), elapsed.sec = numeric(), stringsAsFactors = FALSE)
        }
        total.elapsed <- proc.time()[["elapsed"]] - timing.start
        timing$total.elapsed.sec <- total.elapsed
        timing$fraction.of.total <- if (total.elapsed > 0) {
            timing$elapsed.sec/total.elapsed
        }
        else {
            0
        }
        rownames(timing) <- NULL
        timing
    }
    n <- nrow(X)
    if (!is.null(support.index) && identical(neighborhood.type, "knn") && is.null(nn.index) &&
        is.null(k)) {
        neighborhood.type <- "supplied"
    }
    if (identical(neighborhood.type, "knn")) {
        if (is.null(k)) {
            if (!is.null(nn.index)) {
                k <- ncol(as.matrix(nn.index))
            }
            else {
                stop("k is required when neighborhood.type = 'knn'.", call. = FALSE)
            }
        }
        k <- .validate.ssrhe.positive.integer(k, "k")
        if (k < 2L || k > n) {
            stop("k must be between 2 and nrow(X).", call. = FALSE)
        }
        nn.index <- .validate.ssrhe.nn.index(nn.index, n, k)
        support.size <- rep(k, n)
        add.timing("validation")
        return(list(type = "knn", k = k, nn.index = nn.index, support.index = NULL, support.size = support.size,
            adaptive.k.scale = NA_integer_, radius.rule = NA_character_, radius.factor = NA_real_,
            metadata = list(type = "knn", support.size = support.size, n.topup = rep(0L, n),
                n.truncated = rep(0L, n), timing = finalize.timing())))
    }
    if (identical(neighborhood.type, "supplied")) {
        support.index <- .validate.ssrhe.support.index(support.index, n)
        support.size <- lengths(support.index)
        add.timing("validation")
        return(list(type = "supplied", k = 0L, nn.index = NULL, support.index = support.index,
            support.size = support.size, adaptive.k.scale = NA_integer_, radius.rule = NA_character_,
            radius.factor = NA_real_, metadata = list(type = "supplied", support.size = support.size,
                n.topup = rep(0L, n), n.truncated = rep(0L, n), timing = finalize.timing())))
    }
    if (!identical(neighborhood.type, "adaptive.radius")) {
        stop("Unknown SSRHE neighborhood type.", call. = FALSE)
    }
    if (is.null(adaptive.k.scale)) {
        if (is.null(k)) {
            stop("adaptive.k.scale is required when neighborhood.type = 'adaptive.radius'.",
                call. = FALSE)
        }
        adaptive.k.scale <- k
    }
    adaptive.k.scale <- .validate.ssrhe.positive.integer(adaptive.k.scale, "adaptive.k.scale")
    if (adaptive.k.scale >= n) {
        stop("adaptive.k.scale must be smaller than nrow(X).", call. = FALSE)
    }
    radius.factor <- .validate.ssrhe.numeric.scalar(radius.factor, "radius.factor")
    if (radius.factor <= 0) {
        stop("radius.factor must be positive.", call. = FALSE)
    }
    if (is.null(min.support)) {
        min.support <- .ssrhe.default.min.support(tangent.dim = tangent.dim, tangent.dim.rule = tangent.dim.rule,
            ambient.dim = ncol(X), support.buffer = support.buffer, n = n, derivative.order = derivative.order)
    }
    else {
        min.support <- .validate.ssrhe.positive.integer(min.support, "min.support")
    }
    if (min.support < 2L || min.support > n) {
        stop("min.support must be between 2 and nrow(X).", call. = FALSE)
    }
    if (!is.null(max.support)) {
        max.support <- .validate.ssrhe.positive.integer(max.support, "max.support")
        if (max.support < min.support || max.support > n) {
            stop("max.support must be between min.support and nrow(X).", call. = FALSE)
        }
    }
    add.timing("validation")
    graph <- create.rknn.graph(X = X, type = "adaptive.radius", k.scale = adaptive.k.scale,
        radius.factor = radius.factor, radius.rule = radius.rule, prune.method = "none", connect.components = FALSE,
        return.timing = return.timing, graph.detail = "minimal")
    add.timing("create.graph")
    support.index <- lapply(seq_len(n), function(i) {
        unique(c(i, graph$adj_list[[i]]))
    })
    add.timing("initial.supports")
    n.topup <- integer(n)
    if (identical(support.topup, "nearest")) {
        undersized <- which(lengths(support.index) < min.support)
        for (i in undersized) {
            ids <- support.index[[i]]
            need <- min.support - length(ids)
            add <- .ssrhe.nearest.outside.support(X = X, center = i, exclude = ids, n.add = need)
            ids <- unique(c(ids, add))
            n.topup[[i]] <- max(0L, length(ids) - length(support.index[[i]]))
            support.index[[i]] <- ids
        }
    }
    add.timing("topup")
    n.truncated <- integer(n)
    for (i in seq_len(n)) {
        ids <- support.index[[i]]
        if (!is.null(max.support) && length(ids) > max.support) {
            non.center <- setdiff(ids, i)
            d <- rowSums((t(t(X[non.center, , drop = FALSE]) - X[i, ]))^2)
            keep <- non.center[order(d, non.center)][seq_len(max.support - 1L)]
            n.truncated[[i]] <- length(ids) - max.support
            ids <- c(i, keep)
        }
        else if (ids[[1L]] != i) {
            ids <- c(i, setdiff(ids, i))
        }
        support.index[[i]] <- as.integer(ids)
    }
    add.timing("truncate.reorder")
    support.index <- .validate.ssrhe.support.index(support.index, n)
    support.size <- lengths(support.index)
    add.timing("final.validation")
    list(type = "adaptive.radius", k = 0L, nn.index = NULL, support.index = support.index,
        support.size = support.size, adaptive.k.scale = adaptive.k.scale, radius.rule = radius.rule,
        radius.factor = radius.factor, metadata = list(type = "adaptive.radius", graph = graph,
            adaptive.k.scale = adaptive.k.scale, radius.rule = radius.rule, radius.factor = radius.factor,
            min.support = min.support, max.support = max.support, support.topup = support.topup,
            support.size = support.size, n.topup = n.topup, n.truncated = n.truncated, sigma = graph$sigma,
            graph.timing = graph$timing, timing = finalize.timing()))
}

.ssrhe.nearest.outside.support <- function(X, center, exclude, n.add) {
    if (n.add <= 0L) {
        return(integer())
    }
    n <- nrow(X)
    d <- rowSums((t(t(X) - X[center, ]))^2)
    d[unique(c(center, exclude))] <- Inf
    candidates <- order(d, seq_len(n))
    candidates <- candidates[is.finite(d[candidates])]
    as.integer(candidates[seq_len(min(n.add, length(candidates)))])
}

.ssrhe.triplet.to.sparse <- function(triplet) {
    Matrix::sparseMatrix(i = triplet$i, j = triplet$j, x = triplet$x, dims = triplet$dim)
}

.fit.ssrhe.hessian.gcv.from.operator <- function(operator, y, lambda1.grid, lambda2.grid, weights,
    ridge, trace.method = c("exact", "hutchinson"), trace.n.probes = 50L, trace.seed = NULL,
    verbose = FALSE) {
    n <- ncol(operator$B)
    trace.method <- match.arg(trace.method)
    trace.n.probes <- .validate.ssrhe.positive.integer(trace.n.probes, "trace.n.probes")
    trace.seed <- .validate.ssrhe.optional.seed(trace.seed, "trace.seed")
    if (length(y) != n)
        stop("Internal error: y length mismatch.", call. = FALSE)
    if (length(weights) != n) {
        stop("Internal error: weights length mismatch.", call. = FALSE)
    }
    if (any(!is.finite(y))) {
        stop("GCV requires a finite response vector.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop("GCV requires finite strictly positive weights.", call. = FALSE)
    }
    grid <- expand.grid(lambda1 = lambda1.grid, lambda2 = lambda2.grid, KEEP.OUT.ATTRS = FALSE)
    rows <- vector("list", nrow(grid))
    fits <- vector("list", nrow(grid))
    for (g in seq_len(nrow(grid))) {
        fit <- tryCatch(.fit.ssrhe.hessian.from.operator(operator = operator, y = y, lambda1 = grid$lambda1[g],
            lambda2 = grid$lambda2[g], weights = weights, ridge = ridge, verbose = verbose),
            error = function(e) e)
        if (inherits(fit, "error")) {
            rows[[g]] <- data.frame(lambda1 = grid$lambda1[g], lambda2 = grid$lambda2[g], rss = NA_real_,
                trace.S = NA_real_, trace.se = NA_real_, trace.method = trace.method, trace.n.probes = if (identical(trace.method,
                  "hutchinson")) {
                  trace.n.probes
                }
                else {
                  NA_integer_
                }, edf = NA_real_, gcv = Inf, status = "error", message = conditionMessage(fit),
                stringsAsFactors = FALSE)
            next
        }
        system <- .ssrhe.hessian.system.matrix(B = operator$B, BS = operator$BS, weights = weights,
            lambda1 = grid$lambda1[g], lambda2 = grid$lambda2[g], ridge = ridge)
        grid.seed <- .ssrhe.trace.grid.seed(trace.seed, g)
        trace.info <- tryCatch(.ssrhe.hessian.smoother.trace(system = system, weights = weights,
            method = trace.method, n.probes = trace.n.probes, seed = grid.seed, verbose = verbose),
            error = function(e) e)
        if (inherits(trace.info, "error")) {
            rows[[g]] <- data.frame(lambda1 = grid$lambda1[g], lambda2 = grid$lambda2[g], rss = NA_real_,
                trace.S = NA_real_, trace.se = NA_real_, trace.method = trace.method, trace.n.probes = if (identical(trace.method,
                  "hutchinson")) {
                  trace.n.probes
                }
                else {
                  NA_integer_
                }, edf = NA_real_, gcv = Inf, status = "error", message = conditionMessage(trace.info),
                stringsAsFactors = FALSE)
            next
        }
        trace.S <- trace.info$trace
        rss <- sum(weights * (y - fit$fitted.values)^2)
        denom <- (1 - trace.S/n)^2
        gcv <- if (denom > 0)
            (rss/n)/denom
        else Inf
        rows[[g]] <- data.frame(lambda1 = grid$lambda1[g], lambda2 = grid$lambda2[g], rss = rss,
            trace.S = trace.S, trace.se = trace.info$trace.se, trace.method = trace.info$method,
            trace.n.probes = trace.info$n.probes, edf = trace.S, gcv = gcv, status = "ok",
            message = "", stringsAsFactors = FALSE)
        fits[[g]] <- fit
    }
    gcv.table <- do.call(rbind, rows)
    if (!any(is.finite(gcv.table$gcv))) {
        stop("All SSRHE GCV fits failed.", call. = FALSE)
    }
    selected.idx <- which.min(gcv.table$gcv)
    final <- fits[[selected.idx]]
    final$gcv.table <- gcv.table
    final$selection <- list(rule = "gcv", selected.index = selected.idx, lambda1 = gcv.table$lambda1[selected.idx],
        lambda2 = gcv.table$lambda2[selected.idx], gcv = gcv.table$gcv[selected.idx], rss = gcv.table$rss[selected.idx],
        trace.S = gcv.table$trace.S[selected.idx], trace.se = gcv.table$trace.se[selected.idx],
        trace.method = gcv.table$trace.method[selected.idx], trace.n.probes = gcv.table$trace.n.probes[selected.idx],
        n = n)
    final
}

.fit.ssrhe.hessian.from.operator <- function(operator, y, lambda1, lambda2, weights, ridge,
    verbose = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE regression solves.", call. = FALSE)
    }
    if (is.null(operator$B)) {
        stop("operator must contain B. Rebuild with return.B = TRUE.", call. = FALSE)
    }
    n <- ncol(operator$B)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    Y.raw <- y.info$Y
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    Y.clean <- Y.raw
    Y.clean[!is.finite(Y.clean)] <- 0
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    if (lambda2 > 0 && is.null(operator$BS)) {
        stop("lambda2 > 0 requires an operator with BS. Rebuild with stabilizer = TRUE.", call. = FALSE)
    }
    n.responses <- ncol(Y.clean)
    Y.hat <- matrix(NA_real_, nrow = n, ncol = n.responses)
    residuals <- matrix(NA_real_, nrow = n, ncol = n.responses)
    data.loss <- hessian.energy <- stabilizer.energy <- objective <- numeric(n.responses)
    solve.method <- character(n.responses)
    common.weights <- n.responses == 1L || all(vapply(seq_len(n.responses), function(j) {
        identical(W[, j], W[, 1L])
    }, logical(1)))
    if (common.weights) {
        solved <- .solve.ssrhe.hessian.system(B = operator$B, BS = operator$BS, weights = W[,
            1L], Y = Y.clean, lambda1 = lambda1, lambda2 = lambda2, ridge = ridge, verbose = verbose)
        Y.hat[, ] <- solved$fitted
        solve.method[] <- solved$method
    }
    else {
        for (j in seq_len(n.responses)) {
            solved <- .solve.ssrhe.hessian.system(B = operator$B, BS = operator$BS, weights = W[,
                j], Y = Y.clean[, j, drop = FALSE], lambda1 = lambda1, lambda2 = lambda2, ridge = ridge,
                verbose = verbose)
            Y.hat[, j] <- solved$fitted[, 1L]
            solve.method[j] <- solved$method
        }
    }
    for (j in seq_len(n.responses)) {
        observed <- y.info$observed[, j] & W[, j] > 0
        residuals[observed, j] <- Y.raw[observed, j] - Y.hat[observed, j]
        residuals[!observed, j] <- NA_real_
        data.loss[j] <- 0.5 * sum(W[, j] * (Y.clean[, j] - Y.hat[, j])^2)
        hessian.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$B %*% Y.hat[,
            j]))
        if (!is.null(operator$BS)) {
            stabilizer.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$BS %*%
                Y.hat[, j]))
        }
        else {
            stabilizer.energy[j] <- 0
        }
        objective[j] <- data.loss[j] + 0.5 * lambda1 * hessian.energy[j] + 0.5 * lambda2 *
            stabilizer.energy[j]
    }
    if (!is.null(y.info$col.names)) {
        colnames(Y.hat) <- y.info$col.names
        colnames(residuals) <- y.info$col.names
        names(data.loss) <- y.info$col.names
        names(hessian.energy) <- y.info$col.names
        names(stabilizer.energy) <- y.info$col.names
        names(objective) <- y.info$col.names
    }
    single <- n.responses == 1L
    list(fitted.values = if (single) as.vector(Y.hat) else Y.hat, residuals = if (single) as.vector(residuals) else residuals,
        y = if (single) as.vector(Y.raw) else Y.raw, weights = if (single) as.vector(W) else W,
        lambda = list(lambda1 = lambda1, lambda2 = lambda2, ridge = ridge), objective = if (single) objective[1L] else objective,
        energies = list(data.loss = if (single) data.loss[1L] else data.loss, hessian = if (single) hessian.energy[1L] else hessian.energy,
            stabilizer = if (single) stabilizer.energy[1L] else stabilizer.energy), operator = operator,
        n.responses = n.responses, solver = list(method = if (length(unique(solve.method)) ==
            1L) solve.method[1L] else solve.method), observed = if (single) as.vector(y.info$observed &
            W > 0) else y.info$observed & W > 0)
}

.solve.ssrhe.hessian.system <- function(B, BS, weights, Y, lambda1, lambda2, ridge, verbose = FALSE) {
    n <- nrow(B)
    if (length(weights) != n)
        stop("Internal error: weights length mismatch.", call. = FALSE)
    system <- .ssrhe.hessian.system.matrix(B = B, BS = BS, weights = weights, lambda1 = lambda1,
        lambda2 = lambda2, ridge = ridge)
    rhs <- weights * Y
    fit <- tryCatch({
        as.matrix(Matrix::solve(system, rhs))
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse solve failed; retrying with dense solve().")
        }
        tryCatch(solve(as.matrix(system), as.matrix(rhs)), error = function(e2) {
            stop("SSRHE linear solve failed: ", conditionMessage(e2), call. = FALSE)
        })
    })
    list(fitted = fit, method = "linear.solve")
}

.ssrhe.hessian.system.matrix <- function(B, BS, weights, lambda1, lambda2, ridge) {
    n <- nrow(B)
    system <- lambda1 * B
    if (lambda2 > 0) {
        if (is.null(BS)) {
            stop("lambda2 > 0 requires BS.", call. = FALSE)
        }
        system <- system + lambda2 * BS
    }
    if (ridge > 0) {
        system <- system + Matrix::Diagonal(n, ridge)
    }
    if (any(weights > 0)) {
        system <- system + Matrix::Diagonal(n, weights)
    }
    system
}

.ssrhe.hessian.smoother.trace <- function(system, weights, method = c("exact", "hutchinson"),
    n.probes = 50L, seed = NULL, verbose = FALSE) {
    method <- match.arg(method)
    if (identical(method, "exact")) {
        trace <- .ssrhe.hessian.smoother.trace.exact(system, weights, verbose = verbose)
        return(list(trace = trace, trace.se = NA_real_, method = "exact", n.probes = NA_integer_))
    }
    .ssrhe.hessian.smoother.trace.hutchinson(system = system, weights = weights, n.probes = n.probes,
        seed = seed, verbose = verbose)
}

.ssrhe.hessian.smoother.trace.exact <- function(system, weights, verbose = FALSE) {
    rhs <- Matrix::Diagonal(length(weights), x = weights)
    S <- tryCatch({
        Matrix::solve(system, rhs)
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse smoother-trace solve failed; retrying with dense solve().")
        }
        tryCatch(solve(as.matrix(system), as.matrix(rhs)), error = function(e2) {
            stop("SSRHE smoother-trace solve failed: ", conditionMessage(e2), call. = FALSE)
        })
    })
    if (inherits(S, "Matrix")) {
        sum(as.numeric(Matrix::diag(S)))
    }
    else {
        sum(diag(S))
    }
}

.ssrhe.hessian.smoother.trace.hutchinson <- function(system, weights, n.probes, seed = NULL,
    verbose = FALSE) {
    n.probes <- .validate.ssrhe.positive.integer(n.probes, "n.probes")
    seed <- .validate.ssrhe.optional.seed(seed, "seed")
    n <- length(weights)
    probes <- .ssrhe.with.optional.seed(seed, {
        matrix(sample(c(-1, 1), n * n.probes, replace = TRUE), nrow = n, ncol = n.probes)
    })
    rhs <- weights * probes
    solved <- tryCatch({
        as.matrix(Matrix::solve(system, rhs))
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse Hutchinson trace solve failed; retrying with dense solve().")
        }
        tryCatch(solve(as.matrix(system), as.matrix(rhs)), error = function(e2) {
            stop("SSRHE Hutchinson trace solve failed: ", conditionMessage(e2), call. = FALSE)
        })
    })
    trace.samples <- colSums(probes * solved)
    list(trace = mean(trace.samples), trace.se = if (n.probes > 1L) {
        stats::sd(trace.samples)/sqrt(n.probes)
    } else {
        NA_real_
    }, method = "hutchinson", n.probes = n.probes)
}

.ssrhe.with.optional.seed <- function(seed, expr) {
    if (is.null(seed)) {
        return(force(expr))
    }
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had.seed) {
        old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if (had.seed) {
            assign(".Random.seed", old.seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    force(expr)
}

.ssrhe.trace.grid.seed <- function(seed, grid.index) {
    if (is.null(seed))
        return(NULL)
    as.integer(((seed + grid.index - 1L)%%2147483646L) + 1L)
}

.prepare.ssrhe.response.matrix <- function(y, n, name) {
    is.matrix.input <- is.matrix(y) || inherits(y, "Matrix") || (is.data.frame(y) && ncol(y) >
        1L)
    if (is.matrix.input) {
        Y <- if (inherits(y, "Matrix"))
            as.matrix(y)
        else as.matrix(y)
        if (nrow(Y) != n)
            stop(sprintf("nrow(%s) must be %d.", name, n), call. = FALSE)
        col.names <- colnames(Y)
    }
    else {
        if (length(y) != n)
            stop(sprintf("%s must have length %d.", name, n), call. = FALSE)
        Y <- matrix(y, ncol = 1L)
        col.names <- NULL
    }
    storage.mode(Y) <- "double"
    if (any(is.infinite(Y))) {
        stop(sprintf("%s cannot contain infinite values.", name), call. = FALSE)
    }
    observed <- is.finite(Y)
    list(Y = Y, observed = observed, col.names = col.names)
}

.prepare.ssrhe.weight.matrix <- function(weights, y.info, n) {
    p <- ncol(y.info$Y)
    if (is.null(weights)) {
        W <- matrix(1, nrow = n, ncol = p)
    }
    else if (is.matrix(weights) || inherits(weights, "Matrix") || is.data.frame(weights)) {
        W <- if (inherits(weights, "Matrix"))
            as.matrix(weights)
        else as.matrix(weights)
        if (!identical(dim(W), dim(y.info$Y))) {
            stop("weights matrix must have the same dimensions as y.", call. = FALSE)
        }
    }
    else {
        if (length(weights) != n)
            stop("weights vector must have length nrow(X).", call. = FALSE)
        W <- matrix(as.double(weights), nrow = n, ncol = p)
    }
    storage.mode(W) <- "double"
    if (any(!is.finite(W)) || any(W < 0)) {
        stop("weights must be finite and nonnegative.", call. = FALSE)
    }
    W[!y.info$observed] <- 0
    W
}

.validate.ssrhe.nonnegative.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
        stop(sprintf("%s must be a finite nonnegative numeric scalar.", name), call. = FALSE)
    }
    as.double(x)
}

.validate.ssrhe.optional.seed <- function(x, name) {
    if (is.null(x))
        return(NULL)
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be NULL or a finite nonnegative integer scalar.", name), call. = FALSE)
    }
    as.integer(round(x))
}

.validate.ssrhe.lambda.grid <- function(x, name) {
    if (!is.numeric(x) || !length(x) || any(is.na(x)) || any(!is.finite(x)) || any(x < 0)) {
        stop(sprintf("%s must be a nonempty finite nonnegative numeric vector.", name), call. = FALSE)
    }
    sort(unique(as.double(x)))
}

.prepare.ssrhe.cv.folds <- function(observed, nfolds, fold.id) {
    n <- length(observed)
    if (is.null(fold.id)) {
        nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
        observed.idx <- which(observed)
        if (length(observed.idx) < 2L) {
            stop("At least two observed positive-weight labels are required for CV.", call. = FALSE)
        }
        nfolds <- min(nfolds, length(observed.idx))
        out <- integer(n)
        out[observed.idx] <- rep(seq_len(nfolds), length.out = length(observed.idx))
    }
    else {
        if (length(fold.id) != n) {
            stop("fold.id must have length nrow(X).", call. = FALSE)
        }
        out <- as.integer(fold.id)
        out[is.na(out) | out < 1L | !observed] <- 0L
    }
    folds <- sort(unique(out[out > 0L]))
    if (length(folds) < 2L) {
        stop("At least two nonempty validation folds are required.", call. = FALSE)
    }
    for (ff in folds) {
        n.validation <- sum(out == ff)
        n.training <- sum(observed & out != ff)
        if (n.validation < 1L || n.training < 1L) {
            stop("Each validation fold must leave at least one training label.", call. = FALSE)
        }
    }
    out
}

.fit.ssrhe.hessian.l1.from.operator <- function(operator, y, weights, lambda.grid, lambda.selection,
    n.lambda, fold.id, loss, selection, solver.args, verbose = FALSE) {
    if (is.null(operator$A)) {
        stop("operator must contain A. Rebuild with return.A = TRUE.", call. = FALSE)
    }
    n <- ncol(operator$A)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("SSRHE Hessian L1 fitting currently supports one response vector.", call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    y.raw <- as.vector(y.info$Y)
    weights <- as.vector(W)
    y.clean <- y.raw
    y.clean[!is.finite(y.clean)] <- 0
    D.info <- .ssrhe.hessian.l1.penalty.matrix(operator$A, row.scaling = solver.args$row.scaling)
    diagnostics <- .ssrhe.hessian.l1.diagnostics(A = operator$A, D.info = D.info)
    backend <- solver.args$solver
    path.info <- NULL
    path <- NULL
    if (!identical(backend, "admm")) {
        path.info <- tryCatch(.fit.ssrhe.hessian.l1.path(y = y.clean, weights = weights, D = D.info$D,
            solver.args = solver.args, use.svd = D.info$svd), error = function(e) {
            if (identical(backend, "auto")) {
                structure(list(error = e), class = "ssrhe.l1.path.error")
            }
            else {
                stop(e)
            }
        })
        if (!inherits(path.info, "ssrhe.l1.path.error")) {
            path <- path.info$path
        }
        if (is.null(lambda.grid)) {
            if (inherits(path.info, "ssrhe.l1.path.error")) {
                stop("lambda.grid is required when solver = 'auto' falls back before a genlasso path is available.",
                  call. = FALSE)
            }
            lambda.grid <- .ssrhe.hessian.l1.default.lambda.grid(path, n.lambda)
        }
    }
    else if (is.null(lambda.grid)) {
        stop("lambda.grid is required when solver = 'admm'.", call. = FALSE)
    }
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.", call. = FALSE)
    }
    cv <- NULL
    if (identical(lambda.selection, "cv")) {
        if (is.null(fold.id)) {
            fold.id <- .prepare.ssrhe.cv.folds(observed = as.vector(y.info$observed & weights >
                0), nfolds = solver.args$nfolds, fold.id = NULL)
        }
        cv <- .ssrhe.hessian.l1.cv(y = y.clean, weights = weights, D = D.info$D, lambda.grid = lambda.grid,
            fold.id = fold.id, loss = loss, selection = selection, solver.args = solver.args,
            use.svd = D.info$svd)
        selected.idx <- cv$selected.idx
    }
    else {
        selected.idx <- 1L
    }
    fit.source <- backend
    beta.grid <- NULL
    if (!identical(backend, "admm") && !inherits(path.info, "ssrhe.l1.path.error")) {
        beta.grid <- .ssrhe.hessian.l1.coef.matrix(path, lambda.grid, n)
    }
    if (identical(backend, "admm") || inherits(path.info, "ssrhe.l1.path.error") || any(!is.finite(beta.grid))) {
        if (identical(backend, "genlasso")) {
            stop("SSRHE Hessian L1 path produced non-finite fitted values.", call. = FALSE)
        }
        beta.grid <- .ssrhe.hessian.l1.admm.grid(y = y.clean, weights = weights, D = D.info$D.sparse,
            lambda.grid = lambda.grid, solver.args = solver.args)
        fit.source <- if (identical(backend, "auto"))
            "admm_fallback"
        else "admm"
    }
    fitted <- beta.grid[, selected.idx]
    observed <- as.vector(y.info$observed & weights > 0)
    residuals <- rep(NA_real_, n)
    residuals[observed] <- y.raw[observed] - fitted[observed]
    data.loss <- 0.5 * sum(weights * (y.clean - fitted)^2)
    hessian.l1 <- sum(abs(as.vector(operator$A %*% fitted)))
    lambda <- lambda.grid[selected.idx]
    list(fitted.values = as.vector(fitted), residuals = residuals, y = y.raw, weights = weights,
        lambda = lambda, lambda.grid = lambda.grid, lambda.selection = lambda.selection, objective = data.loss +
            lambda * hessian.l1, energies = list(data.loss = data.loss, hessian.l1 = hessian.l1),
        operator = operator, beta.grid = beta.grid, path = path, cv = cv, fold.id = fold.id,
        observed = observed, n.responses = 1L, diagnostics = diagnostics, solver = list(backend = fit.source,
            requested = backend, representation = D.info$representation, svd = D.info$svd,
            row.scaling = D.info$row.scaling, n.observed = if (is.null(path.info) || inherits(path.info,
                "ssrhe.l1.path.error")) {
                sum(observed)
            } else {
                path.info$n.observed
            }, admm = if (grepl("admm", fit.source, fixed = TRUE)) {
                attr(beta.grid, "admm")
            } else {
                NULL
            }), selection = list(rule = selection, loss = loss, selected.index = selected.idx,
            lambda = lambda))
}

.ssrhe.hessian.l1.penalty.matrix <- function(A, row.scaling = c("none", "l2")) {
    row.scaling <- match.arg(row.scaling)
    A.sparse <- methods::as(A, "dgCMatrix")
    row.scale <- rep(1, nrow(A.sparse))
    if (identical(row.scaling, "l2")) {
        row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
        nz <- is.finite(row.norm) & row.norm > 0
        row.scale[nz] <- 1/row.norm[nz]
        A.sparse <- Matrix::Diagonal(x = row.scale) %*% A.sparse
        A.sparse <- methods::as(A.sparse, "dgCMatrix")
    }
    use.svd <- nrow(A) >= ncol(A)
    list(D = if (use.svd) as.matrix(A.sparse) else A.sparse, D.sparse = A.sparse, svd = use.svd,
        representation = if (use.svd) "dense" else "sparse", row.scaling = row.scaling, row.scale = row.scale)
}

.ssrhe.hessian.l1.diagnostics <- function(A, D.info) {
    A.sparse <- methods::as(A, "dgCMatrix")
    scaled <- D.info$D.sparse
    row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
    scaled.row.norm <- sqrt(Matrix::rowSums(scaled^2))
    list(nrow = nrow(A.sparse), ncol = ncol(A.sparse), nnzero = Matrix::nnzero(A.sparse), density = Matrix::nnzero(A.sparse)/(nrow(A.sparse) *
        ncol(A.sparse)), row.norm = list(min = suppressWarnings(min(row.norm, na.rm = TRUE)),
        median = stats::median(row.norm), max = suppressWarnings(max(row.norm, na.rm = TRUE)),
        zero = sum(!is.finite(row.norm) | row.norm == 0)), scaled.row.norm = list(min = suppressWarnings(min(scaled.row.norm,
        na.rm = TRUE)), median = stats::median(scaled.row.norm), max = suppressWarnings(max(scaled.row.norm,
        na.rm = TRUE)), zero = sum(!is.finite(scaled.row.norm) | scaled.row.norm == 0)), row.scaling = D.info$row.scaling,
        representation = D.info$representation, svd = D.info$svd)
}

.fit.ssrhe.hessian.l1.path <- function(y, weights, D, solver.args, use.svd = FALSE) {
    observed <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(observed)) {
        stop("At least one observed positive-weight response is required.", call. = FALSE)
    }
    if (all(observed) && all(abs(weights - 1) < sqrt(.Machine$double.eps))) {
        path <- suppressWarnings(genlasso::genlasso(y = y, D = D, approx = solver.args$approx,
            maxsteps = solver.args$maxsteps, minlam = solver.args$minlam, rtol = solver.args$rtol,
            btol = solver.args$btol, eps = solver.args$eps, verbose = solver.args$verbose,
            svd = use.svd))
    }
    else {
        obs <- which(observed)
        sw <- sqrt(weights[obs])
        X.design <- as.matrix(Matrix::sparseMatrix(i = seq_along(obs), j = obs, x = sw, dims = c(length(obs),
            length(y)), giveCsparse = TRUE))
        path <- suppressWarnings(genlasso::genlasso(y = sw * y[obs], X = X.design, D = D, approx = solver.args$approx,
            maxsteps = solver.args$maxsteps, minlam = solver.args$minlam, rtol = solver.args$rtol,
            btol = solver.args$btol, eps = solver.args$eps, verbose = solver.args$verbose,
            svd = use.svd))
    }
    list(path = path, n.observed = sum(observed))
}

.ssrhe.hessian.l1.default.lambda.grid <- function(path, n.lambda) {
    lambda.path <- path$lambda
    lambda.path <- lambda.path[is.finite(lambda.path) & lambda.path >= 0]
    lambda.max <- suppressWarnings(max(lambda.path, na.rm = TRUE))
    if (!is.finite(lambda.max) || lambda.max <= 0) {
        lambda.max <- 1
    }
    n.lambda <- max(2L, as.integer(n.lambda))
    lambda.positive <- lambda.path[lambda.path > 0]
    complete.path <- isTRUE(path$completepath)
    lambda.min <- if (complete.path || !length(lambda.positive)) {
        lambda.max * 1e-04
    }
    else {
        max(min(lambda.positive) * (1 + 1e-08), lambda.max * 1e-04)
    }
    grid <- exp(seq(log(lambda.max), log(lambda.min), length.out = n.lambda - 1L))
    if (complete.path)
        grid <- c(grid, 0)
    unique(grid)
}

.ssrhe.hessian.l1.coef.matrix <- function(path, lambda.grid, n) {
    beta <- tryCatch(as.matrix(stats::coef(path, lambda = lambda.grid)$beta), error = function(e) NULL)
    if (is.null(beta) || !identical(dim(beta), c(as.integer(n), length(lambda.grid)))) {
        beta <- vapply(lambda.grid, function(lambda) {
            tryCatch(as.vector(stats::coef(path, lambda = lambda)$beta), error = function(e) rep(NA_real_,
                n))
        }, numeric(n))
    }
    beta
}

.ssrhe.soft.threshold <- function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
}

.ssrhe.hessian.l1.admm <- function(y, weights, D, lambda, solver.args) {
    D <- methods::as(D, "dgCMatrix")
    n <- length(y)
    m <- nrow(D)
    weights <- as.numeric(weights)
    y <- as.numeric(y)
    rho <- solver.args$admm.rho
    AtA <- Matrix::crossprod(D)
    system <- AtA * rho + Matrix::Diagonal(n, x = weights)
    ridge <- max(1e-10, sqrt(.Machine$double.eps) * max(1, mean(Matrix::diag(system))))
    system <- system + Matrix::Diagonal(n, x = rep(ridge, n))
    factor <- Matrix::Cholesky(system, LDL = FALSE, Imult = 0)
    wy <- weights * y
    beta <- rep(0, n)
    z <- rep(0, m)
    u <- rep(0, m)
    converged <- FALSE
    iter <- 0L
    primal <- dual <- Inf
    for (iter in seq_len(solver.args$admm.maxiter)) {
        z.old <- z
        rhs <- wy + rho * as.vector(Matrix::crossprod(D, z - u))
        beta <- as.vector(Matrix::solve(factor, rhs))
        Dbeta <- as.vector(D %*% beta)
        z <- .ssrhe.soft.threshold(Dbeta + u, lambda/rho)
        u <- u + Dbeta - z
        primal <- sqrt(sum((Dbeta - z)^2))
        dual <- rho * sqrt(sum(as.vector(Matrix::crossprod(D, z - z.old))^2))
        eps.pri <- sqrt(m) * solver.args$admm.abstol + solver.args$admm.reltol * max(sqrt(sum(Dbeta^2)),
            sqrt(sum(z^2)))
        eps.dual <- sqrt(n) * solver.args$admm.abstol + solver.args$admm.reltol * sqrt(sum(as.vector(rho *
            Matrix::crossprod(D, u))^2))
        if (primal <= eps.pri && dual <= eps.dual) {
            converged <- TRUE
            break
        }
    }
    attr(beta, "admm") <- list(iterations = iter, converged = converged, primal.residual = primal,
        dual.residual = dual, rho = rho)
    beta
}

.ssrhe.hessian.l1.admm.grid <- function(y, weights, D, lambda.grid, solver.args) {
    admm.info <- vector("list", length(lambda.grid))
    beta <- vapply(seq_along(lambda.grid), function(ii) {
        fit <- .ssrhe.hessian.l1.admm(y = y, weights = weights, D = D, lambda = lambda.grid[ii],
            solver.args = solver.args)
        admm.info[[ii]] <<- attr(fit, "admm")
        as.vector(fit)
    }, numeric(length(y)))
    attr(beta, "admm") <- admm.info
    beta
}

.ssrhe.hessian.l1.cv <- function(y, weights, D, lambda.grid, fold.id, loss, selection, solver.args,
    use.svd = FALSE) {
    n <- length(y)
    folds <- sort(unique(fold.id[fold.id > 0L]))
    cv.errors <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda.grid))
    rownames(cv.errors) <- paste0("Fold", folds)
    colnames(cv.errors) <- format(signif(lambda.grid, 6), scientific = TRUE)
    for (ii in seq_along(folds)) {
        fold <- folds[ii]
        test <- which(fold.id == fold)
        train.weights <- weights
        train.weights[test] <- 0
        if (identical(solver.args$solver, "admm")) {
            beta <- tryCatch(.ssrhe.hessian.l1.admm.grid(y = y, weights = train.weights, D = D,
                lambda.grid = lambda.grid, solver.args = solver.args), error = function(e) NULL)
        }
        else {
            path.info <- tryCatch(.fit.ssrhe.hessian.l1.path(y = y, weights = train.weights,
                D = D, solver.args = solver.args, use.svd = use.svd), error = function(e) NULL)
            if (is.null(path.info)) {
                beta <- NULL
            }
            else {
                beta <- .ssrhe.hessian.l1.coef.matrix(path.info$path, lambda.grid, n)
            }
            if (identical(solver.args$solver, "auto") && (is.null(beta) || any(!is.finite(beta)))) {
                beta <- tryCatch(.ssrhe.hessian.l1.admm.grid(y = y, weights = train.weights,
                  D = methods::as(D, "dgCMatrix"), lambda.grid = lambda.grid, solver.args = solver.args),
                  error = function(e) NULL)
            }
        }
        if (is.null(beta)) {
            next
        }
        beta.test <- beta[test, , drop = FALSE]
        finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
        target <- matrix(y[test], nrow = length(test), ncol = length(lambda.grid))
        w.test <- matrix(weights[test], nrow = length(test), ncol = length(lambda.grid))
        if (identical(loss, "mse")) {
            fold.error <- colSums(w.test * (target - beta.test)^2)/colSums(w.test)
        }
        else {
            fold.error <- colSums(w.test * abs(target - beta.test))/colSums(w.test)
        }
        fold.error[!finite.pred] <- NA_real_
        cv.errors[ii, ] <- fold.error
    }
    mean.error <- colMeans(cv.errors, na.rm = TRUE)
    mean.error[!is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L)
            return(Inf)
        stats::sd(x)/sqrt(length(x))
    })
    best.idx <- which.min(mean.error)
    if (!is.finite(mean.error[best.idx])) {
        stop("All SSRHE Hessian L1 CV fits failed for the supplied lambda grid.", call. = FALSE)
    }
    if (identical(selection, "one.se")) {
        threshold <- mean.error[best.idx] + se[best.idx]
        eligible <- which(mean.error <= threshold)
        selected.idx <- eligible[which.max(lambda.grid[eligible])]
    }
    else {
        selected.idx <- best.idx
    }
    list(fold.id = fold.id, lambda.grid = lambda.grid, fold.errors = cv.errors, mean.error = mean.error,
        se = se, best.idx = best.idx, selected.idx = selected.idx, lambda.min = lambda.grid[best.idx],
        lambda = lambda.grid[selected.idx], error.min = mean.error[best.idx], error.selected = mean.error[selected.idx],
        selection = selection, loss = loss)
}

.validate.ssrhe.hessian.l1.lambda.grid <- function(lambda.grid, lambda.selection) {
    if (is.null(lambda.grid)) {
        if (identical(lambda.selection, "fixed")) {
            stop("lambda.grid is required when lambda.selection = 'fixed'.", call. = FALSE)
        }
        return(NULL)
    }
    lambda.grid <- .validate.ssrhe.lambda.grid(lambda.grid, "lambda.grid")
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.", call. = FALSE)
    }
    lambda.grid
}

.validate.ssrhe.hessian.l1.solver.args <- function(n.lambda, nfolds, fold.id, observed, maxsteps,
    minlam, approx, rtol, btol, eps, solver = c("genlasso", "admm", "auto"), row.scaling = c("none",
        "l2"), admm.rho = 1, admm.maxiter = 2000L, admm.abstol = 1e-04, admm.reltol = 0.001,
    verbose) {
    n.lambda <- .validate.ssrhe.positive.integer(n.lambda, "n.lambda")
    nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
    maxsteps <- .validate.ssrhe.positive.integer(maxsteps, "maxsteps")
    solver <- match.arg(solver)
    row.scaling <- match.arg(row.scaling)
    minlam <- .validate.ssrhe.nonnegative.scalar(minlam, "minlam")
    rtol <- .validate.ssrhe.nonnegative.scalar(rtol, "rtol")
    btol <- .validate.ssrhe.nonnegative.scalar(btol, "btol")
    eps <- .validate.ssrhe.nonnegative.scalar(eps, "eps")
    admm.rho <- .validate.ssrhe.numeric.scalar(admm.rho, "admm.rho")
    if (!is.finite(admm.rho) || admm.rho <= 0) {
        stop("admm.rho must be a finite positive numeric scalar.", call. = FALSE)
    }
    admm.maxiter <- .validate.ssrhe.positive.integer(admm.maxiter, "admm.maxiter")
    admm.abstol <- .validate.ssrhe.nonnegative.scalar(admm.abstol, "admm.abstol")
    admm.reltol <- .validate.ssrhe.nonnegative.scalar(admm.reltol, "admm.reltol")
    if (!is.null(fold.id)) {
        fold.id <- .prepare.ssrhe.cv.folds(observed, nfolds, fold.id)
    }
    list(n.lambda = n.lambda, nfolds = nfolds, fold.id = fold.id, maxsteps = maxsteps, minlam = minlam,
        approx = isTRUE(approx), rtol = rtol, btol = btol, eps = eps, solver = solver, row.scaling = row.scaling,
        admm.rho = admm.rho, admm.maxiter = admm.maxiter, admm.abstol = admm.abstol, admm.reltol = admm.reltol,
        verbose = isTRUE(verbose))
}
