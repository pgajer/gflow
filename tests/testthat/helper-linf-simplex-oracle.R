.linf.simplex.active.faces.oracle <- function(x, tol = sqrt(.Machine$double.eps)) {
    which(abs(x - 1) <= tol)
}

.linf.simplex.pair.dist.oracle <- function(x, y, tol = sqrt(.Machine$double.eps)) {
    ax <- .linf.simplex.active.faces.oracle(x, tol)
    ay <- .linf.simplex.active.faces.oracle(y, tol)

    best <- Inf
    for (face.x in ax) {
        for (face.y in ay) {
            if (identical(face.x, face.y)) {
                keep <- seq_along(x) != face.x
                dist2 <- sum((x[keep] - y[keep])^2)
            } else {
                keep <- seq_along(x) != face.x & seq_along(x) != face.y
                crossing <- 2 - x[face.y] - y[face.x]
                dist2 <- crossing^2 + sum((x[keep] - y[keep])^2)
            }
            best <- min(best, dist2)
        }
    }

    sqrt(best)
}

.linf.simplex.knn.oracle <- function(X, k, tol = sqrt(.Machine$double.eps)) {
    stopifnot(is.matrix(X), is.numeric(X), k >= 1L, k <= nrow(X))
    n <- nrow(X)
    indices <- matrix(NA_integer_, nrow = n, ncol = k)
    distances <- matrix(NA_real_, nrow = n, ncol = k)

    for (i in seq_len(n)) {
        d <- vapply(seq_len(n), function(j) {
            .linf.simplex.pair.dist.oracle(X[i, ], X[j, ], tol)
        }, numeric(1))
        ord <- order(d, seq_len(n) - 1L)
        keep <- ord[seq_len(k)]
        indices[i, ] <- keep - 1L
        distances[i, ] <- d[keep]
    }

    list(indices = indices, distances = distances)
}

.linf.simplex.backend.knn <- function(X, k, tol = sqrt(.Machine$double.eps)) {
    .Call("S_linf_simplex_knn",
          X,
          as.integer(k),
          as.double(tol),
          PACKAGE = "gflow")
}

.random.linf.simplex.matrix <- function(n, p, ridge.every = 0L, corner.every = 0L, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    X <- matrix(stats::runif(n * p), nrow = n, ncol = p)
    X <- X / apply(X, 1L, max)

    if (ridge.every > 0L) {
        ridge.rows <- seq.int(ridge.every, n, by = ridge.every)
        for (i in ridge.rows) {
            faces <- sample.int(p, 2L)
            X[i, faces] <- 1
        }
    }

    if (corner.every > 0L) {
        corner.rows <- seq.int(corner.every, n, by = corner.every)
        for (i in corner.rows) {
            n.faces <- min(p, 3L)
            faces <- sample.int(p, n.faces)
            X[i, faces] <- 1
        }
    }

    storage.mode(X) <- "double"
    X
}

.compare.linf.simplex.backend.to.oracle <- function(X,
                                                    k,
                                                    tol = sqrt(.Machine$double.eps),
                                                    distance.tol = 1e-10,
                                                    allow.tie.index.drift = TRUE) {
    oracle <- .linf.simplex.knn.oracle(X, k, tol)
    backend <- .linf.simplex.backend.knn(X, k, tol)

    distance.diff <- abs(backend$distances - oracle$distances)
    max.distance.diff <- max(distance.diff)

    index.identical <- identical(unname(backend$indices), unname(oracle$indices))
    distance.ok <- isTRUE(all(distance.diff <= distance.tol))

    if (allow.tie.index.drift) {
        membership.ok <- TRUE
        for (i in seq_len(nrow(X))) {
            kth <- oracle$distances[i, k]
            backend.d <- vapply(backend$indices[i, ] + 1L, function(j) {
                .linf.simplex.pair.dist.oracle(X[i, ], X[j, ], tol)
            }, numeric(1))
            if (any(backend.d > kth + distance.tol)) {
                membership.ok <- FALSE
                break
            }
        }
    } else {
        membership.ok <- index.identical
    }

    list(
        oracle = oracle,
        backend = backend,
        max.distance.diff = max.distance.diff,
        distance.ok = distance.ok,
        index.identical = index.identical,
        membership.ok = membership.ok,
        passed = distance.ok && membership.ok
    )
}
