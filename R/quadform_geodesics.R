.validate.quadform.index <- function(p, index.k) {
    if (!is.numeric(index.k) || length(index.k) != 1L || !is.finite(index.k) ||
        index.k != floor(index.k) || index.k < 0L || index.k > p) {
        stop("'index.k' must be an integer between 0 and ncol(X).", call. = FALSE)
    }
    as.integer(index.k)
}

.validate.quadform.data.matrix <- function(X) {
    if (!(is.matrix(X) || is.data.frame(X))) {
        stop("'X' must be a matrix or data frame.", call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X)) {
        stop("'X' must contain numeric data.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("'X' cannot contain NA, NaN, or Inf values.", call. = FALSE)
    }
    if (nrow(X) < 1L || ncol(X) < 1L) {
        stop("'X' must contain at least one observation and one column.", call. = FALSE)
    }
    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }
    X
}

.quadform.signs <- function(p, index.k) {
    index.k <- .validate.quadform.index(p, index.k)
    c(rep.int(1, index.k), rep.int(-1, p - index.k))
}

.quadform.value <- function(X, index.k) {
    signs <- .quadform.signs(ncol(X), index.k)
    as.numeric((X^2) %*% signs)
}

#' Embed a Quadratic Hypersurface
#'
#' @description
#' Embeds parameter points into the graph of the quadratic form
#' \deqn{q(x)=\sum_{i=1}^k x_i^2-\sum_{i=k+1}^n x_i^2.}
#'
#' @param X Numeric matrix or data frame with parameter points in rows.
#' @param index.k Integer between 0 and `ncol(X)`. Number of positive-square
#'   terms in the quadratic form.
#'
#' @return A numeric matrix with `ncol(X) + 1` columns. The first `ncol(X)`
#'   columns are `X`; the last column is `q(X)`.
#'
#' @examples
#' X <- rbind(c(1, 2), c(2, 1))
#' quadform.embed(X, index.k = 1)
#'
#' @export
quadform.embed <- function(X, index.k) {
    X <- .validate.quadform.data.matrix(X)
    q <- .quadform.value(X, index.k)
    cbind(X, q = q)
}

#' Compute Gradients of a Quadratic Form
#'
#' @description
#' Computes gradients of
#' \deqn{q(x)=\sum_{i=1}^k x_i^2-\sum_{i=k+1}^n x_i^2}
#' at parameter points.
#'
#' @inheritParams quadform.embed
#'
#' @return A numeric matrix whose rows are \eqn{\nabla q(x)}.
#'
#' @examples
#' quadform.gradient(matrix(c(1, 2), nrow = 1), index.k = 1)
#'
#' @export
quadform.gradient <- function(X, index.k) {
    X <- .validate.quadform.data.matrix(X)
    signs <- .quadform.signs(ncol(X), index.k)
    t(t(2 * X) * signs)
}

#' Compute the Induced Metric of a Quadratic Hypersurface
#'
#' @description
#' For the graph embedding \eqn{F(x)=(x,q(x))}, the induced Riemannian metric is
#' \deqn{G(x)=I+\nabla q(x)\nabla q(x)^T.}
#'
#' @inheritParams quadform.embed
#'
#' @return If `X` has one row, a numeric matrix. Otherwise, an array with
#'   dimensions `c(nrow(X), ncol(X), ncol(X))`.
#'
#' @examples
#' quadform.metric(matrix(c(1, 2), nrow = 1), index.k = 1)
#'
#' @export
quadform.metric <- function(X, index.k) {
    X <- .validate.quadform.data.matrix(X)
    grad <- quadform.gradient(X, index.k)
    p <- ncol(X)
    out <- array(0, dim = c(nrow(X), p, p))
    eye <- diag(p)
    for (i in seq_len(nrow(X))) {
        g <- grad[i, ]
        out[i, , ] <- eye + tcrossprod(g)
    }
    if (nrow(X) == 1L) {
        return(out[1L, , ])
    }
    out
}

#' Compute Exact Segment Lengths on a Quadratic Hypersurface
#'
#' @description
#' Computes the length of the embedded straight parameter segment from `u` to
#' `v` on the quadratic graph hypersurface. For
#' \eqn{x(t)=u+t(v-u)}, the length is evaluated in closed form as
#' \deqn{\int_0^1 \sqrt{\|v-u\|_2^2 + (a+bt)^2}\,dt,}
#' where \eqn{a = 2u^T S(v-u)}, \eqn{b = 2(v-u)^T S(v-u)}, and
#' \eqn{S = diag(1,\ldots,1,-1,\ldots,-1)} with `index.k` positive entries.
#'
#' @param u Numeric vector, one endpoint in parameter coordinates.
#' @param v Numeric vector, the other endpoint in parameter coordinates.
#' @inheritParams quadform.embed
#' @param tol Positive numeric tolerance for treating the affine term slope as
#'   zero.
#'
#' @return Numeric scalar segment length.
#'
#' @examples
#' quadform.edge.length(c(0, 0), c(1, 0), index.k = 2)
#'
#' @export
quadform.edge.length <- function(u, v, index.k, tol = sqrt(.Machine$double.eps)) {
    if (!is.numeric(u) || !is.numeric(v) || length(u) != length(v) ||
        any(!is.finite(u)) || any(!is.finite(v))) {
        stop("'u' and 'v' must be finite numeric vectors of equal length.", call. = FALSE)
    }
    p <- length(u)
    signs <- .quadform.signs(p, index.k)
    h <- v - u
    A <- sum(h^2)
    if (A == 0) {
        return(0)
    }
    a <- 2 * sum(u * signs * h)
    b <- 2 * sum(h * signs * h)
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("'tol' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (abs(b) <= tol) {
        return(sqrt(A + a^2))
    }
    sqrt.A <- sqrt(A)
    antiderivative <- function(y) {
        0.5 * (y * sqrt(A + y^2) + A * asinh(y / sqrt.A))
    }
    as.numeric((antiderivative(a + b) - antiderivative(a)) / b)
}

.quadform.reference.vertices.2d <- function(X, domain.radius, grid.size) {
    coords <- seq(-domain.radius, domain.radius, length.out = grid.size)
    grid <- expand.grid(x = coords, y = coords)
    inside <- rowSums(as.matrix(grid)^2) <= domain.radius^2 * (1 + 1e-12)
    grid <- as.matrix(grid[inside, , drop = FALSE])
    storage.mode(grid) <- "double"
    list(vertices = rbind(X, grid), grid = grid, coords = coords, inside = inside)
}

.quadform.reference.grid.edges.2d <- function(offset, inside, grid.size) {
    id <- matrix(NA_integer_, nrow = grid.size, ncol = grid.size)
    id[inside] <- offset + seq_len(sum(inside))
    offsets <- rbind(c(1L, 0L), c(0L, 1L), c(1L, 1L), c(1L, -1L))
    rows <- list()
    cursor <- 0L
    for (i in seq_len(grid.size)) {
        for (j in seq_len(grid.size)) {
            from <- id[i, j]
            if (is.na(from)) {
                next
            }
            for (r in seq_len(nrow(offsets))) {
                ii <- i + offsets[r, 1L]
                jj <- j + offsets[r, 2L]
                if (ii < 1L || ii > grid.size || jj < 1L || jj > grid.size) {
                    next
                }
                to <- id[ii, jj]
                if (!is.na(to)) {
                    cursor <- cursor + 1L
                    rows[[cursor]] <- c(from, to)
                }
            }
        }
    }
    if (!cursor) {
        return(matrix(integer(), ncol = 2L))
    }
    do.call(rbind, rows)
}

.quadform.reference.sample.edges <- function(X, grid, sample.connection.k) {
    n <- nrow(X)
    ng <- nrow(grid)
    rows <- vector("list", n * sample.connection.k)
    cursor <- 0L
    for (i in seq_len(n)) {
        d <- rowSums((t(t(grid) - X[i, ]))^2)
        nbrs <- order(d, seq_len(ng))[seq_len(sample.connection.k)]
        for (j in nbrs) {
            cursor <- cursor + 1L
            rows[[cursor]] <- c(i, n + j)
        }
    }
    do.call(rbind, rows[seq_len(cursor)])
}

.with.quadform.seed <- function(seed, expr) {
    if (is.null(seed)) {
        return(force(expr))
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed)) {
        stop("'seed' must be NULL or a finite integer scalar.", call. = FALSE)
    }
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old.seed <- if (had.seed) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had.seed) {
            assign(".Random.seed", old.seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(as.integer(seed))
    force(expr)
}

.quadform.sample.uniform.parameter.disk <- function(n, domain.radius) {
    theta <- stats::runif(n, min = 0, max = 2 * pi)
    radius <- domain.radius * sqrt(stats::runif(n))
    X <- cbind(
        x1 = radius * cos(theta),
        x2 = radius * sin(theta)
    )
    storage.mode(X) <- "double"
    X
}

#' Estimate Reference Geodesic Distances on a 2D Quadratic Surface
#'
#' @description
#' Computes a numerical reference geodesic distance matrix for points sampled
#' from a two-dimensional quadratic graph surface. The method builds a dense
#' reference graph on a regular parameter-disk grid, includes the sample points
#' as vertices, weights all reference edges by exact quadratic-surface segment
#' length, and returns shortest-path distances between the sample vertices.
#'
#' @param X Numeric matrix or data frame with two parameter-coordinate columns.
#' @param index.k Integer between 0 and 2. Number of positive-square terms.
#' @param domain.radius Positive numeric scalar or `NULL`. Radius of the
#'   parameter disk covered by the reference grid. If `NULL`, uses the maximum
#'   sample radius.
#' @param grid.size Odd or even integer at least 5. Number of grid coordinates
#'   along each axis before clipping to the disk.
#' @param sample.connection.k Positive integer. Number of nearest reference-grid
#'   vertices connected to each sample point.
#'
#' @return A list with `distances`, the sample-by-sample reference geodesic
#'   distance matrix, plus reference graph diagnostics.
#'
#' @examples
#' X <- rbind(c(0, 0), c(0.5, 0), c(0, 0.5))
#' ref <- quadform.reference.geodesics(X, index.k = 2, grid.size = 9)
#' ref$distances
#'
#' @export
quadform.reference.geodesics <- function(X,
                                         index.k,
                                         domain.radius = NULL,
                                         grid.size = 51,
                                         sample.connection.k = 8) {
    X <- .validate.numeric.data.matrix(X)
    if (ncol(X) != 2L) {
        stop("'quadform.reference.geodesics' currently supports only 2D parameter coordinates.",
             call. = FALSE)
    }
    .validate.quadform.index(2L, index.k)
    if (is.null(domain.radius)) {
        domain.radius <- max(sqrt(rowSums(X^2)))
        if (domain.radius == 0) {
            domain.radius <- 1
        }
    }
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (any(sqrt(rowSums(X^2)) > domain.radius * (1 + 1e-10))) {
        stop("All rows of 'X' must lie inside the parameter disk.", call. = FALSE)
    }
    if (!is.numeric(grid.size) || length(grid.size) != 1L ||
        !is.finite(grid.size) || grid.size != floor(grid.size) || grid.size < 5L) {
        stop("'grid.size' must be an integer at least 5.", call. = FALSE)
    }
    grid.size <- as.integer(grid.size)
    ref <- .quadform.reference.vertices.2d(X, domain.radius, grid.size)
    if (!is.numeric(sample.connection.k) || length(sample.connection.k) != 1L ||
        !is.finite(sample.connection.k) || sample.connection.k != floor(sample.connection.k) ||
        sample.connection.k < 1L) {
        stop("'sample.connection.k' must be a positive integer.", call. = FALSE)
    }
    sample.connection.k <- min(sample.connection.k, nrow(ref$grid))
    sample.connection.k <- as.integer(sample.connection.k)
    n <- nrow(X)
    grid.edges <- .quadform.reference.grid.edges.2d(n, ref$inside, grid.size)
    sample.edges <- .quadform.reference.sample.edges(X, ref$grid, sample.connection.k)
    edge.matrix <- unique(rbind(grid.edges, sample.edges))
    edge.weight <- apply(edge.matrix, 1L, function(e) {
        quadform.edge.length(ref$vertices[e[[1L]], ], ref$vertices[e[[2L]], ], index.k)
    })
    graph <- igraph::graph_from_edgelist(edge.matrix, directed = FALSE)
    distances <- igraph::distances(graph, v = seq_len(n), to = seq_len(n),
                                   weights = edge.weight)
    list(
        distances = as.matrix(distances),
        vertices = ref$vertices,
        edge_matrix = edge.matrix,
        edge_weight = as.numeric(edge.weight),
        n_sample_vertices = n,
        n_reference_vertices = nrow(ref$vertices),
        n_edges = nrow(edge.matrix),
        domain_radius = domain.radius,
        grid_size = grid.size,
        sample_connection_k = sample.connection.k,
        index_k = as.integer(index.k)
    )
}

#' Estimate 2D Quadratic-Surface Geodesics with the C++ Grid Engine
#'
#' @description
#' Computes numerical reference geodesic distances for points on a
#' two-dimensional quadratic graph surface using the package's C++ grid
#' shortest-path engine. A disk-clipped regular parameter grid is built once,
#' grid edges are weighted by exact quadratic-surface segment length, and each
#' sample point is connected to nearby grid vertices before running Dijkstra
#' shortest paths.
#'
#' @inheritParams quadform.reference.geodesics
#' @param oracle Character scalar. `"none"` returns only the continuum grid
#'   estimate. `"sample.path"` also returns a finite-sample oracle distance
#'   matrix guided by the grid geodesic path.
#' @param oracle.tube.radius Positive numeric scalar or `NULL`. Ambient
#'   embedded-space tube radius used to select sample points close to each
#'   grid-geodesic path. If `NULL`, uses the median `oracle.tube.k`-nearest
#'   neighbor distance among embedded sample points.
#' @param oracle.tube.k Positive integer. Local scale used when
#'   `oracle.tube.radius = NULL`.
#' @param return.oracle.paths Logical scalar. If `TRUE`, return the selected
#'   sample-index sequence for each ordered source-target pair.
#'
#' @return A list with `distances`, the sample-by-sample continuum grid
#'   geodesic distance matrix, plus grid diagnostics. If
#'   `oracle = "sample.path"`, the list also contains `oracle_distances`,
#'   `oracle_n_points`, `oracle_status`, `oracle_tube_radius`, and
#'   `oracle_method`.
#'
#' @examples
#' X <- rbind(c(0, 0), c(0.5, 0), c(0, 0.5))
#' ref <- quadform.grid.geodesic.distances(X, index.k = 2, grid.size = 21)
#' ref$distances
#'
#' @export
quadform.grid.geodesic.distances <- function(X,
                                             index.k,
                                             domain.radius = NULL,
                                             grid.size = 51,
                                             sample.connection.k = 8,
                                             oracle = c("none", "sample.path"),
                                             oracle.tube.radius = NULL,
                                             oracle.tube.k = 8,
                                             return.oracle.paths = FALSE) {
    X <- .validate.numeric.data.matrix(X)
    if (ncol(X) != 2L) {
        stop("'quadform.grid.geodesic.distances' currently supports only 2D parameter coordinates.",
             call. = FALSE)
    }
    .validate.quadform.index(2L, index.k)
    if (is.null(domain.radius)) {
        domain.radius <- max(sqrt(rowSums(X^2)))
        if (domain.radius == 0) {
            domain.radius <- 1
        }
    }
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (any(sqrt(rowSums(X^2)) > domain.radius * (1 + 1e-10))) {
        stop("All rows of 'X' must lie inside the parameter disk.", call. = FALSE)
    }
    if (!is.numeric(grid.size) || length(grid.size) != 1L ||
        !is.finite(grid.size) || grid.size != floor(grid.size) || grid.size < 5L) {
        stop("'grid.size' must be an integer at least 5.", call. = FALSE)
    }
    if (!is.numeric(sample.connection.k) || length(sample.connection.k) != 1L ||
        !is.finite(sample.connection.k) || sample.connection.k != floor(sample.connection.k) ||
        sample.connection.k < 1L) {
        stop("'sample.connection.k' must be a positive integer.", call. = FALSE)
    }
    oracle <- match.arg(oracle)
    if (!is.numeric(oracle.tube.k) || length(oracle.tube.k) != 1L ||
        !is.finite(oracle.tube.k) || oracle.tube.k != floor(oracle.tube.k) ||
        oracle.tube.k < 1L) {
        stop("'oracle.tube.k' must be a positive integer.", call. = FALSE)
    }
    if (!is.logical(return.oracle.paths) || length(return.oracle.paths) != 1L ||
        is.na(return.oracle.paths)) {
        stop("'return.oracle.paths' must be TRUE or FALSE.", call. = FALSE)
    }
    with.oracle <- identical(oracle, "sample.path")
    if (with.oracle && is.null(oracle.tube.radius)) {
        oracle.tube.radius <- .quadform.default.oracle.tube.radius(
            X, index.k, as.integer(oracle.tube.k)
        )
    }
    if (is.null(oracle.tube.radius)) {
        oracle.tube.radius <- NA_real_
    } else if (!is.numeric(oracle.tube.radius) || length(oracle.tube.radius) != 1L ||
        !is.finite(oracle.tube.radius) || oracle.tube.radius <= 0) {
        stop("'oracle.tube.radius' must be NULL or a positive finite numeric scalar.",
             call. = FALSE)
    }
    out <- rcpp_quadform_grid_geodesic_distances(
        X,
        as.integer(index.k),
        as.numeric(domain.radius),
        as.integer(grid.size),
        as.integer(sample.connection.k),
        with.oracle,
        as.numeric(oracle.tube.radius),
        as.integer(oracle.tube.k),
        isTRUE(return.oracle.paths)
    )
    out$distances <- as.matrix(out$distances)
    if (with.oracle) {
        out$oracle_distances <- as.matrix(out$oracle_distances)
        out$oracle_n_points <- as.matrix(out$oracle_n_points)
        status.labels <- c("same", "ok", "direct", "no_path")
        status <- matrix(
            status.labels[pmax(0L, out$oracle_status_code) + 1L],
            nrow = nrow(out$oracle_status_code),
            ncol = ncol(out$oracle_status_code)
        )
        out$oracle_status <- status
    }
    out
}

.quadform.default.oracle.tube.radius <- function(X, index.k, oracle.tube.k) {
    n <- nrow(X)
    if (n < 2L) {
        stop("At least two sample points are required for sample-path oracle distances.",
             call. = FALSE)
    }
    k <- min(as.integer(oracle.tube.k), n - 1L)
    X.embed <- quadform.embed(X, index.k = index.k)
    kth <- numeric(n)
    for (i in seq_len(n)) {
        d <- sqrt(rowSums((t(t(X.embed) - X.embed[i, ]))^2))
        d[[i]] <- Inf
        kth[[i]] <- sort(d, partial = k)[[k]]
    }
    as.numeric(stats::median(kth))
}

.quadform.reference.grid.2d <- function(domain.radius, grid.size) {
    coords <- seq(-domain.radius, domain.radius, length.out = grid.size)
    grid <- expand.grid(x1 = coords, x2 = coords)
    X <- as.matrix(grid)
    inside <- rowSums(X^2) <= domain.radius^2 * (1 + 1e-12)
    X <- X[inside, , drop = FALSE]
    storage.mode(X) <- "double"
    X
}

.sample.quadform.reference.pairs <- function(domain.radius,
                                             reference.grid.size,
                                             n.sources,
                                             targets.per.source,
                                             seed) {
    grid <- .quadform.reference.grid.2d(domain.radius, reference.grid.size)
    n.grid <- nrow(grid)
    pair.index <- .with.quadform.seed(seed, {
        sources <- sample.int(n.grid, n.sources, replace = n.sources > n.grid)
        targets <- matrix(sample.int(n.grid, n.sources * targets.per.source,
                                     replace = TRUE),
                          nrow = n.sources, ncol = targets.per.source)
        for (i in seq_len(n.sources)) {
            same <- targets[i, ] == sources[[i]]
            while (any(same)) {
                targets[i, same] <- sample.int(n.grid, sum(same), replace = TRUE)
                same <- targets[i, ] == sources[[i]]
            }
        }
        cbind(
            source = rep(sources, each = targets.per.source),
            target = as.vector(t(targets))
        )
    })
    pair.points <- cbind(
        grid[pair.index[, "source"], , drop = FALSE],
        grid[pair.index[, "target"], , drop = FALSE]
    )
    colnames(pair.points) <- c("x1_source", "x2_source", "x1_target", "x2_target")
    list(grid = grid, pair_index = pair.index, pair_points = pair.points)
}

.summarize.quadform.grid.error <- function(grid.size, result, reference.distances,
                                           runtime.sec) {
    d <- as.numeric(result$distances)
    ref <- as.numeric(reference.distances)
    keep <- is.finite(d) & is.finite(ref) & ref > 0
    rel <- abs(d[keep] - ref[keep]) / ref[keep]
    abs.err <- abs(d[keep] - ref[keep])
    data.frame(
        grid_size = as.integer(grid.size),
        n_vertices = as.integer(result$n_vertices),
        n_edges = as.integer(result$n_edges),
        n_pairs = length(rel),
        frob_rel_error = sqrt(sum((d[keep] - ref[keep])^2) / sum(ref[keep]^2)),
        rms_rel_error = sqrt(mean(rel^2)),
        mean_rel_error = mean(rel),
        median_rel_error = unname(stats::median(rel)),
        q90_rel_error = unname(stats::quantile(rel, 0.90)),
        q95_rel_error = unname(stats::quantile(rel, 0.95)),
        q99_rel_error = unname(stats::quantile(rel, 0.99)),
        max_rel_error = max(rel),
        median_abs_error = unname(stats::median(abs.err)),
        q95_abs_error = unname(stats::quantile(abs.err, 0.95)),
        max_abs_error = max(abs.err),
        mean_endpoint_snap = mean(c(result$source_snap_distance,
                                    result$target_snap_distance)),
        q95_endpoint_snap = unname(stats::quantile(c(result$source_snap_distance,
                                                     result$target_snap_distance), 0.95)),
        runtime_sec = as.numeric(runtime.sec),
        stringsAsFactors = FALSE
    )
}

#' Calibrate Grid Error for 2D Quadratic-Surface Geodesics
#'
#' @description
#' Estimates discretization error for the C++ grid geodesic engine on a
#' two-dimensional quadratic graph surface. The function samples many
#' source-target point pairs from a high-resolution reference grid, computes
#' their graph-geodesic distances on a set of candidate grids, and compares
#' them to distances on a much finer reference grid.
#'
#' @param index.k Integer between 0 and 2. Number of positive-square terms in
#'   the quadratic form.
#' @param domain.radius Positive numeric scalar. Radius of the parameter disk.
#' @param grid.sizes Integer vector of candidate grid sizes.
#' @param reference.grid.size Integer scalar. High-resolution grid size used as
#'   the numerical reference.
#' @param n.sources Positive integer. Number of source points sampled from the
#'   reference grid.
#' @param targets.per.source Positive integer. Number of target points sampled
#'   per source. The total number of pairwise distances is
#'   `n.sources * targets.per.source`.
#' @param seed `NULL` or finite integer scalar. If supplied, pair sampling is
#'   reproducible and the caller's random-number state is restored.
#'
#' @return A list of class `"quadform_grid_geodesic_calibration"` with
#'   `summary`, `pair_results`, sampled `pair_points`, and `metadata`.
#'
#' @examples
#' cal <- quadform.grid.geodesic.calibration(
#'   index.k = 2,
#'   grid.sizes = c(11, 21),
#'   reference.grid.size = 31,
#'   n.sources = 3,
#'   targets.per.source = 2,
#'   seed = 1
#' )
#' cal$summary
#'
#' @export
quadform.grid.geodesic.calibration <- function(index.k,
                                               domain.radius = 1,
                                               grid.sizes = c(101, 201, 251, 501),
                                               reference.grid.size = 1001,
                                               n.sources = 50,
                                               targets.per.source = 20,
                                               seed = NULL) {
    .validate.quadform.index(2L, index.k)
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    validate.grid.size <- function(x, name) {
        if (!is.numeric(x) || any(!is.finite(x)) || any(x != floor(x)) ||
            any(x < 5L)) {
            stop("'", name, "' must contain integer values at least 5.", call. = FALSE)
        }
        as.integer(x)
    }
    grid.sizes <- unique(validate.grid.size(grid.sizes, "grid.sizes"))
    reference.grid.size <- validate.grid.size(reference.grid.size, "reference.grid.size")
    if (length(reference.grid.size) != 1L) {
        stop("'reference.grid.size' must be a single integer at least 5.",
             call. = FALSE)
    }
    if (!reference.grid.size %in% grid.sizes) {
        grid.sizes <- sort(c(grid.sizes, reference.grid.size))
    } else {
        grid.sizes <- sort(grid.sizes)
    }
    if (!is.numeric(n.sources) || length(n.sources) != 1L ||
        !is.finite(n.sources) || n.sources != floor(n.sources) || n.sources < 1L) {
        stop("'n.sources' must be a positive integer.", call. = FALSE)
    }
    if (!is.numeric(targets.per.source) || length(targets.per.source) != 1L ||
        !is.finite(targets.per.source) ||
        targets.per.source != floor(targets.per.source) ||
        targets.per.source < 1L) {
        stop("'targets.per.source' must be a positive integer.", call. = FALSE)
    }
    pairs <- .sample.quadform.reference.pairs(
        domain.radius = domain.radius,
        reference.grid.size = reference.grid.size,
        n.sources = as.integer(n.sources),
        targets.per.source = as.integer(targets.per.source),
        seed = seed
    )

    pair.results <- list()
    runtimes <- numeric(length(grid.sizes))
    names(runtimes) <- as.character(grid.sizes)
    for (grid.size in grid.sizes) {
        elapsed <- system.time({
            pair.results[[as.character(grid.size)]] <-
                rcpp_quadform_grid_pair_distances(
                    as.integer(index.k),
                    as.numeric(domain.radius),
                    as.integer(grid.size),
                    pairs$pair_points
                )
        })
        runtimes[[as.character(grid.size)]] <- unname(elapsed[["elapsed"]])
    }
    reference.result <- pair.results[[as.character(reference.grid.size)]]
    reference.distances <- reference.result$distances
    summary <- do.call(rbind, lapply(grid.sizes, function(grid.size) {
        .summarize.quadform.grid.error(
            grid.size = grid.size,
            result = pair.results[[as.character(grid.size)]],
            reference.distances = reference.distances,
            runtime.sec = runtimes[[as.character(grid.size)]]
        )
    }))
    rownames(summary) <- NULL

    out <- list(
        summary = summary,
        pair_results = pair.results,
        pair_points = pairs$pair_points,
        pair_index = pairs$pair_index,
        metadata = list(
            index_k = as.integer(index.k),
            domain_radius = as.numeric(domain.radius),
            grid_sizes = as.integer(grid.sizes),
            reference_grid_size = as.integer(reference.grid.size),
            n_sources = as.integer(n.sources),
            targets_per_source = as.integer(targets.per.source),
            n_pairs = nrow(pairs$pair_points),
            seed = seed
        )
    )
    class(out) <- c("quadform_grid_geodesic_calibration", "list")
    out
}

#' Sample a 2D Quadratic Surface Dataset with Reference Geodesics
#'
#' @description
#' Samples parameter points from the disk
#' \eqn{\{x \in \mathbb{R}^2 : \|x\|_2 \le r\}} and embeds them into the
#' quadratic graph surface
#' \deqn{q(x)=\sum_{i=1}^k x_i^2-\sum_{i=k+1}^2 x_i^2.}
#' The currently supported sampling method is explicitly named
#' `"uniform.parameter.disk"`: points are uniform in the parameter disk, not
#' with respect to the induced surface-area measure. The function also builds
#' the regular parameter-disk reference grid used by
#' [quadform.reference.geodesics()] and returns sample-by-sample reference
#' geodesic distances.
#'
#' @param n Positive integer. Number of sample points.
#' @param index.k Integer between 0 and 2. Number of positive-square terms in
#'   the quadratic form.
#' @param domain.radius Positive numeric scalar. Radius of the parameter disk.
#' @param sample.method Character scalar. Currently only
#'   `"uniform.parameter.disk"` is supported.
#' @param grid.size Integer at least 5. Number of reference-grid coordinates
#'   along each axis before clipping to the parameter disk.
#' @param sample.connection.k Positive integer. Number of nearest reference-grid
#'   vertices connected to each sample point when computing reference geodesics.
#' @param seed `NULL` or finite integer scalar. If supplied, sampling is
#'   reproducible and the caller's random-number state is restored.
#'
#' @return A list of class `"quadform_sample_dataset"` with entries:
#' \describe{
#'   \item{X_param}{Sample points in two-dimensional parameter coordinates.}
#'   \item{X_embed}{Sample points embedded as \code{(x1, x2, q)} in
#'     \eqn{\mathbb{R}^3}.}
#'   \item{q}{Quadratic-form values at the sample points.}
#'   \item{D_geodesic}{Sample-by-sample reference geodesic distance matrix.}
#'   \item{distances}{Alias of \code{D_geodesic}.}
#'   \item{reference}{Reference grid and graph payload used to compute
#'     \code{D_geodesic}.}
#'   \item{metadata}{Dataset parameters, including the explicit sampling
#'     method.}
#' }
#'
#' @examples
#' ds <- quadform.sample.dataset(
#'   n = 10,
#'   index.k = 1,
#'   grid.size = 9,
#'   seed = 1
#' )
#' dim(ds$X_embed)
#' dim(ds$D_geodesic)
#'
#' @export
quadform.sample.dataset <- function(n,
                                    index.k,
                                    domain.radius = 1,
                                    sample.method = c("uniform.parameter.disk"),
                                    grid.size = 51,
                                    sample.connection.k = 8,
                                    seed = NULL) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n != floor(n) || n < 1L) {
        stop("'n' must be a positive integer.", call. = FALSE)
    }
    n <- as.integer(n)
    .validate.quadform.index(2L, index.k)
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    sample.method <- match.arg(sample.method)

    X.param <- .with.quadform.seed(
        seed,
        .quadform.sample.uniform.parameter.disk(n, domain.radius)
    )
    X.embed <- quadform.embed(X.param, index.k = index.k)
    ref <- quadform.reference.geodesics(
        X.param,
        index.k = index.k,
        domain.radius = domain.radius,
        grid.size = grid.size,
        sample.connection.k = sample.connection.k
    )
    grid.idx <- seq.int(n + 1L, nrow(ref$vertices))
    grid.param <- ref$vertices[grid.idx, , drop = FALSE]
    vertices.embed <- quadform.embed(ref$vertices, index.k = index.k)

    out <- list(
        X_param = X.param,
        X_embed = X.embed,
        q = as.numeric(X.embed[, "q"]),
        D_geodesic = ref$distances,
        distances = ref$distances,
        reference = list(
            grid_param = grid.param,
            grid_embed = quadform.embed(grid.param, index.k = index.k),
            vertices_param = ref$vertices,
            vertices_embed = vertices.embed,
            sample_vertex = seq_len(n),
            grid_vertex = grid.idx,
            edge_matrix = ref$edge_matrix,
            edge_weight = ref$edge_weight,
            n_reference_vertices = ref$n_reference_vertices,
            n_edges = ref$n_edges,
            grid_size = ref$grid_size,
            sample_connection_k = ref$sample_connection_k
        ),
        metadata = list(
            dim = 2L,
            index_k = as.integer(index.k),
            domain_radius = as.numeric(domain.radius),
            sample_method = sample.method,
            grid_size = as.integer(ref$grid_size),
            sample_connection_k = as.integer(ref$sample_connection_k),
            seed = seed
        )
    )
    class(out) <- c("quadform_sample_dataset", "list")
    out
}
