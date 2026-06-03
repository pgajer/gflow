# Shared local PCA chart-dimension utilities.

.local.pca.auto.chart.dim <- function(X, support.size = NULL,
                                      min.support = NULL,
                                      degree = 2L,
                                      variance.threshold = 0.95,
                                      eigengap.threshold = 4,
                                      max.anchors = 60L,
                                      distance.matrix = NULL,
                                      support.index = NULL,
                                      support.metric = "coordinates",
                                      return.diagnostics = FALSE) {
    result <- .local.pca.auto.chart.dim.result(
        X = X,
        support.size = support.size,
        min.support = min.support,
        degree = degree,
        variance.threshold = variance.threshold,
        eigengap.threshold = eigengap.threshold,
        max.anchors = max.anchors,
        distance.matrix = distance.matrix,
        support.index = support.index,
        support.metric = support.metric
    )
    if (isTRUE(return.diagnostics)) return(result)
    result$chart.dim
}

.local.pca.auto.chart.dim.with.metric <- function(
    X, support.size = NULL,
    min.support = NULL,
    degree = 2L,
    variance.threshold = 0.95,
    eigengap.threshold = 4,
    max.anchors = 60L,
    operator.distance.matrix = NULL,
    operator.support.metric = "coordinates",
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator")) {

    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    operator.support.metric <- as.character(operator.support.metric[[1L]])
    if (identical(operator.support.metric, "graph.geodesic") &&
        is.null(operator.distance.matrix)) {
        stop("Operator-support auto chart dimension requires operator distances.",
             call. = FALSE)
    }

    coordinate <- NULL
    operator <- NULL
    need.coordinate <- auto.chart.support.metric %in% c("coordinates", "both") ||
        identical(auto.chart.selection.metric, "coordinates") ||
        identical(operator.support.metric, "coordinates")
    need.operator <- auto.chart.support.metric %in% c("operator", "both") ||
        identical(auto.chart.selection.metric, "operator")

    if (need.coordinate) {
        coordinate <- .local.pca.auto.chart.dim(
            X = X,
            support.size = support.size,
            min.support = min.support,
            degree = degree,
            variance.threshold = variance.threshold,
            eigengap.threshold = eigengap.threshold,
            max.anchors = max.anchors,
            distance.matrix = NULL,
            support.metric = "coordinates",
            return.diagnostics = TRUE
        )
    }
    if (need.operator) {
        operator <- if (identical(operator.support.metric, "coordinates")) {
            if (is.null(coordinate)) {
                .local.pca.auto.chart.dim(
                    X = X,
                    support.size = support.size,
                    min.support = min.support,
                    degree = degree,
                    variance.threshold = variance.threshold,
                    eigengap.threshold = eigengap.threshold,
                    max.anchors = max.anchors,
                    distance.matrix = NULL,
                    support.metric = "coordinates",
                    return.diagnostics = TRUE
                )
            } else {
                coordinate
            }
        } else {
            .local.pca.auto.chart.dim(
                X = X,
                support.size = support.size,
                min.support = min.support,
                degree = degree,
                variance.threshold = variance.threshold,
                eigengap.threshold = eigengap.threshold,
                max.anchors = max.anchors,
                distance.matrix = operator.distance.matrix,
                support.metric = operator.support.metric,
                return.diagnostics = TRUE
            )
        }
    }

    selected <- if (identical(auto.chart.selection.metric, "operator")) {
        operator
    } else {
        coordinate
    }
    if (is.null(selected)) {
        stop("Internal error: no auto chart-dimension diagnostic was selected.",
             call. = FALSE)
    }
    selected$auto.chart.support.metric <- auto.chart.support.metric
    selected$auto.chart.selection.metric <- auto.chart.selection.metric
    selected$operator.support.metric <- operator.support.metric
    selected$summary$auto.chart.support.metric <- auto.chart.support.metric
    selected$summary$auto.chart.selection.metric <- auto.chart.selection.metric
    selected$summary$operator.support.metric <- operator.support.metric
    selected$all.diagnostics <- list(
        coordinates = coordinate,
        operator = operator
    )
    selected
}

.local.pca.auto.chart.dim.result <- function(X, support.size = NULL,
                                             min.support = NULL,
                                             degree = 2L,
                                             variance.threshold = 0.95,
                                             eigengap.threshold = 4,
                                             max.anchors = 60L,
                                             distance.matrix = NULL,
                                             support.index = NULL,
                                             support.metric = "coordinates") {
    X <- as.matrix(X)
    n <- nrow(X)
    m <- ncol(X)
    if (is.null(support.metric) || !length(support.metric) ||
        all(is.na(support.metric))) {
        support.metric <- "coordinates"
    }
    support.metric <- as.character(support.metric[[1L]])
    empty.diag <- data.frame(
        anchor.index = integer(),
        support.size = integer(),
        support.metric = character(),
        variance.dim = integer(),
        eigengap.dim = integer(),
        selected.local.dim = integer(),
        support.cap = integer(),
        max.eigengap.ratio = numeric(),
        participation.rank = numeric(),
        status = character(),
        stringsAsFactors = FALSE
    )
    make.result <- function(chart.dim, diagnostics = empty.diag,
                            fallback = FALSE) {
        chart.dim <- as.integer(max(1L, min(m, chart.dim)))
        capped <- diagnostics$selected.local.dim >= diagnostics$support.cap
        if (!length(capped)) capped <- logical(0)
        list(
            chart.dim = chart.dim,
            diagnostics = diagnostics,
            summary = list(
                support.metric = support.metric,
                variance.threshold = variance.threshold,
                eigengap.threshold = eigengap.threshold,
                max.anchors = as.integer(max.anchors),
                n.anchors = nrow(diagnostics),
                median.local.dim = if (nrow(diagnostics)) {
                    stats::median(diagnostics$selected.local.dim, na.rm = TRUE)
                } else {
                    NA_real_
                },
                iqr.local.dim = if (nrow(diagnostics)) {
                    stats::IQR(diagnostics$selected.local.dim, na.rm = TRUE)
                } else {
                    NA_real_
                },
                min.local.dim = if (nrow(diagnostics)) {
                    min(diagnostics$selected.local.dim, na.rm = TRUE)
                } else {
                    NA_integer_
                },
                max.local.dim = if (nrow(diagnostics)) {
                    max(diagnostics$selected.local.dim, na.rm = TRUE)
                } else {
                    NA_integer_
                },
                n.capped.by.support = sum(capped, na.rm = TRUE),
                fallback.used = isTRUE(fallback)
            )
        )
    }
    if (m <= 1L || n <= 2L) return(make.result(1L, fallback = TRUE))
    k <- if (!is.null(support.size) && is.finite(support.size)) {
        as.integer(support.size)
    } else if (!is.null(min.support) && is.finite(min.support)) {
        as.integer(min.support + 1L)
    } else {
        min(n, max(12L, min(30L, n)))
    }
    k <- max(3L, min(n, k))
    if (!is.null(distance.matrix)) {
        distance.matrix <- as.matrix(distance.matrix)
        if (!identical(dim(distance.matrix), c(n, n))) {
            stop("'distance.matrix' must be an nrow(X) by nrow(X) matrix.",
                 call. = FALSE)
        }
    }
    if (!is.null(support.index)) {
        if (!is.list(support.index) || length(support.index) != n) {
            stop("'support.index' must be a list with one element per row of X.",
                 call. = FALSE)
        }
        support.index <- lapply(seq_len(n), function(i) {
            ids <- as.integer(support.index[[i]])
            ids <- unique(ids[is.finite(ids) & ids >= 1L & ids <= n])
            if (!i %in% ids) ids <- c(i, ids)
            ids
        })
    }
    anchor.index <- if (n <= max.anchors) {
        seq_len(n)
    } else {
        unique(as.integer(round(seq(1, n, length.out = max.anchors))))
    }
    D <- if (is.null(distance.matrix) && is.null(support.index)) {
        as.matrix(stats::dist(X))
    } else {
        distance.matrix
    }
    rows <- vector("list", length(anchor.index))
    for (aa in seq_along(anchor.index)) {
        ii <- anchor.index[[aa]]
        cand <- if (!is.null(support.index)) {
            support.index[[ii]]
        } else {
            ord <- order(D[ii, ], na.last = NA)
            ord[seq_len(min(k, length(ord)))]
        }
        cand <- unique(as.integer(cand))
        cand <- cand[is.finite(cand) & cand >= 1L & cand <= n]
        if (!length(cand)) {
            rows[[aa]] <- .local.pca.auto.chart.dim.row(
                anchor = ii, support.size = 0L, support.metric = support.metric,
                status = "empty_support"
            )
            next
        }
        centered <- sweep(X[cand, , drop = FALSE], 2L,
                          X[ii, , drop = TRUE], "-")
        sv <- tryCatch(svd(centered, nu = 0L, nv = 0L)$d,
                       error = function(e) numeric(0))
        sv <- sv[is.finite(sv) & sv > 0]
        if (!length(sv)) {
            rows[[aa]] <- .local.pca.auto.chart.dim.row(
                anchor = ii, support.size = length(cand),
                support.metric = support.metric, status = "zero_spectrum"
            )
            next
        }
        energy <- sv^2
        total <- sum(energy)
        if (!is.finite(total) || total <= 0) {
            rows[[aa]] <- .local.pca.auto.chart.dim.row(
                anchor = ii, support.size = length(cand),
                support.metric = support.metric, status = "bad_spectrum"
            )
            next
        }
        support.cap <- .local.pca.max.chart.dim.for.support(
            n.support = max(1L, length(cand) - 1L),
            degree = degree,
            ambient.dim = m
        )
        max.dim <- min(m, support.cap, length(sv), max(1L, length(cand) - 1L))
        var.dim <- which(cumsum(energy) / total >= variance.threshold)[1L]
        var.dim <- min(max.dim, max(1L, var.dim))
        gap.dim <- NA_integer_
        max.gap <- NA_real_
        if (max.dim >= 2L) {
            ratios <- sv[seq_len(max.dim - 1L)] /
                pmax(sv[seq_len(max.dim - 1L) + 1L], .Machine$double.eps)
            best <- which.max(ratios)
            max.gap <- ratios[[best]]
            if (length(best) && is.finite(ratios[[best]]) &&
                ratios[[best]] >= eigengap.threshold) {
                gap.dim <- as.integer(best)
            }
        }
        selected <- if (is.na(gap.dim)) var.dim else min(var.dim, gap.dim)
        participation <- sum(energy)^2 / sum(energy^2)
        rows[[aa]] <- .local.pca.auto.chart.dim.row(
            anchor = ii,
            support.size = length(cand),
            support.metric = support.metric,
            variance.dim = var.dim,
            eigengap.dim = gap.dim,
            selected.local.dim = selected,
            support.cap = support.cap,
            max.eigengap.ratio = max.gap,
            participation.rank = participation,
            status = "ok"
        )
    }
    diagnostics <- do.call(rbind, rows)
    if (is.null(diagnostics)) diagnostics <- empty.diag
    rownames(diagnostics) <- NULL
    ok <- diagnostics$status == "ok" &
        is.finite(diagnostics$selected.local.dim)
    fallback.cap <- .local.pca.max.chart.dim.for.support(
        n.support = max(1L, k - 1L),
        degree = degree,
        ambient.dim = m
    )
    if (!any(ok)) {
        return(make.result(as.integer(max(1L, min(fallback.cap, degree))),
                           diagnostics, fallback = TRUE))
    }
    dim <- as.integer(max(
        1L,
        min(m, fallback.cap,
            stats::median(diagnostics$selected.local.dim[ok], na.rm = TRUE))
    ))
    make.result(dim, diagnostics, fallback = FALSE)
}

.local.pca.auto.chart.dim.row <- function(anchor, support.size,
                                          support.metric,
                                          variance.dim = NA_integer_,
                                          eigengap.dim = NA_integer_,
                                          selected.local.dim = NA_integer_,
                                          support.cap = NA_integer_,
                                          max.eigengap.ratio = NA_real_,
                                          participation.rank = NA_real_,
                                          status) {
    data.frame(
        anchor.index = as.integer(anchor),
        support.size = as.integer(support.size),
        support.metric = as.character(support.metric),
        variance.dim = as.integer(variance.dim),
        eigengap.dim = as.integer(eigengap.dim),
        selected.local.dim = as.integer(selected.local.dim),
        support.cap = as.integer(support.cap),
        max.eigengap.ratio = as.numeric(max.eigengap.ratio),
        participation.rank = as.numeric(participation.rank),
        status = as.character(status),
        stringsAsFactors = FALSE
    )
}

.local.pca.max.chart.dim.for.support <- function(n.support, degree,
                                                 ambient.dim) {
    n.support <- as.integer(max(1L, n.support))
    degree <- as.integer(max(0L, degree))
    ambient.dim <- as.integer(max(1L, ambient.dim))
    ok <- vapply(seq_len(ambient.dim), function(d) {
        choose(d + degree, degree) <= n.support
    }, logical(1))
    if (!any(ok)) return(1L)
    as.integer(max(which(ok)))
}
