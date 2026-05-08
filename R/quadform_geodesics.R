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

.validate.quadform.coefficients <- function(coefficients, p) {
    if (is.null(coefficients)) {
        return(rep.int(1, p))
    }
    if (!is.numeric(coefficients) || length(coefficients) != p ||
        any(!is.finite(coefficients)) || any(coefficients <= 0)) {
        stop("'coefficients' must be NULL or a positive finite numeric vector with length ncol(X).",
             call. = FALSE)
    }
    as.numeric(coefficients)
}

.validate.quadform.domain.shape <- function(domain.shape) {
    match.arg(domain.shape, c("disk", "square"))
}

.quadform.default.domain.radius <- function(X, domain.shape) {
    if (identical(domain.shape, "disk")) {
        radius <- max(sqrt(rowSums(X^2)))
    } else {
        radius <- max(abs(X))
    }
    if (radius == 0) {
        radius <- 1
    }
    radius
}

.validate.quadform.reference.domain.shape <- function(domain.shape, p) {
    if (p == 2L) {
        return(.validate.quadform.domain.shape(domain.shape))
    }
    if (p == 3L) {
        return(match.arg(domain.shape, c("ball", "cube")))
    }
    stop("Reference-domain helpers currently support only p = 2 or p = 3.",
         call. = FALSE)
}

.quadform.points.inside.domain <- function(X, domain.radius, domain.shape,
                                           tol = 1e-10) {
    if (identical(domain.shape, "disk")) {
        sqrt(rowSums(X^2)) <= domain.radius * (1 + tol)
    } else {
        apply(abs(X), 1L, max) <= domain.radius * (1 + tol)
    }
}

.quadform.points.inside.domain.nd <- function(X, domain.radius, domain.shape,
                                              tol = 1e-10) {
    if (domain.shape %in% c("disk", "ball")) {
        sqrt(rowSums(X^2)) <= domain.radius * (1 + tol)
    } else {
        apply(abs(X), 1L, max) <= domain.radius * (1 + tol)
    }
}

.quadform.domain.label <- function(domain.shape) {
    if (identical(domain.shape, "disk")) "parameter disk" else "parameter square"
}

.quadform.value <- function(X, index.k, coefficients = NULL) {
    signs <- .quadform.signs(ncol(X), index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, ncol(X))
    as.numeric((X^2) %*% (signs * coefficients))
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
#' @param coefficients `NULL` or positive numeric vector with length
#'   `ncol(X)`. If supplied, the quadratic graph is
#'   \eqn{q(x)=\sum_i s_i c_i x_i^2}, where `s_i` is determined by
#'   `index.k`.
#'
#' @return A numeric matrix with `ncol(X) + 1` columns. The first `ncol(X)`
#'   columns are `X`; the last column is `q(X)`.
#'
#' @examples
#' X <- rbind(c(1, 2), c(2, 1))
#' quadform.embed(X, index.k = 1)
#'
#' @export
quadform.embed <- function(X, index.k, coefficients = NULL) {
    X <- .validate.quadform.data.matrix(X)
    q <- .quadform.value(X, index.k, coefficients = coefficients)
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
quadform.gradient <- function(X, index.k, coefficients = NULL) {
    X <- .validate.quadform.data.matrix(X)
    signs <- .quadform.signs(ncol(X), index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, ncol(X))
    t(t(2 * X) * signs * coefficients)
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
quadform.metric <- function(X, index.k, coefficients = NULL) {
    X <- .validate.quadform.data.matrix(X)
    grad <- quadform.gradient(X, index.k, coefficients = coefficients)
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
quadform.edge.length <- function(u, v, index.k, coefficients = NULL,
                                 tol = sqrt(.Machine$double.eps)) {
    if (!is.numeric(u) || !is.numeric(v) || length(u) != length(v) ||
        any(!is.finite(u)) || any(!is.finite(v))) {
        stop("'u' and 'v' must be finite numeric vectors of equal length.", call. = FALSE)
    }
    p <- length(u)
    signs <- .quadform.signs(p, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, p)
    h <- v - u
    A <- sum(h^2)
    if (A == 0) {
        return(0)
    }
    a <- 2 * sum(u * signs * coefficients * h)
    b <- 2 * sum(h * signs * coefficients * h)
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

#' Compute Vectorized Exact Segment Lengths on a Quadratic Hypersurface
#'
#' @description
#' Vectorized C++ implementation of `quadform.edge.length()` for matching rows
#' of endpoint matrices. This is mainly intended for building reference graphs
#' whose edge weights are exact lengths of embedded straight parameter
#' segments.
#'
#' @param U,V Numeric matrices or data frames with matching dimensions. Row
#'   `i` of `U` and row `i` of `V` define the two endpoints of segment `i`.
#' @inheritParams quadform.embed
#'
#' @return Numeric vector of segment lengths.
#'
#' @examples
#' U <- rbind(c(0, 0), c(0, 0))
#' V <- rbind(c(1, 0), c(0, 1))
#' quadform.edge.lengths(U, V, index.k = 2)
#'
#' @export
quadform.edge.lengths <- function(U, V, index.k, coefficients = NULL) {
    U <- .validate.quadform.data.matrix(U)
    V <- .validate.quadform.data.matrix(V)
    if (!identical(dim(U), dim(V))) {
        stop("'U' and 'V' must have the same dimensions.", call. = FALSE)
    }
    p <- ncol(U)
    index.k <- .validate.quadform.index(p, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, p)
    rcpp_quadform_edge_lengths(U, V, index.k, coefficients)
}

.quadform.reference.vertices.2d <- function(X, domain.radius, grid.size,
                                            domain.shape = c("disk", "square")) {
    domain.shape <- .validate.quadform.domain.shape(domain.shape)
    coords <- seq(-domain.radius, domain.radius, length.out = grid.size)
    grid <- expand.grid(x = coords, y = coords)
    if (identical(domain.shape, "disk")) {
        inside <- rowSums(as.matrix(grid)^2) <= domain.radius^2 * (1 + 1e-12)
    } else {
        inside <- rep(TRUE, nrow(grid))
    }
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

.quadform.sample.radial.parameter.disk <- function(n, domain.radius, bias) {
    theta <- stats::runif(n, min = 0, max = 2 * pi)
    u <- stats::runif(n)
    if (identical(bias, "center")) {
        radius <- domain.radius * u
    } else if (identical(bias, "boundary")) {
        radius <- domain.radius * sqrt(sqrt(u))
    } else {
        stop("Unknown disk radial bias.", call. = FALSE)
    }
    X <- cbind(
        x1 = radius * cos(theta),
        x2 = radius * sin(theta)
    )
    storage.mode(X) <- "double"
    X
}

.quadform.sample.uniform.parameter.square <- function(n, domain.radius) {
    X <- cbind(
        x1 = stats::runif(n, min = -domain.radius, max = domain.radius),
        x2 = stats::runif(n, min = -domain.radius, max = domain.radius)
    )
    storage.mode(X) <- "double"
    X
}

.quadform.sample.radial.parameter.square <- function(n, domain.radius, bias) {
    theta <- stats::runif(n, min = 0, max = 2 * pi)
    u <- stats::runif(n)
    r.max <- 1 / pmax(abs(cos(theta)), abs(sin(theta)))
    if (identical(bias, "center")) {
        radius <- domain.radius * r.max * u
    } else if (identical(bias, "boundary")) {
        radius <- domain.radius * r.max * sqrt(sqrt(u))
    } else {
        stop("Unknown square radial bias.", call. = FALSE)
    }
    X <- cbind(
        x1 = radius * cos(theta),
        x2 = radius * sin(theta)
    )
    storage.mode(X) <- "double"
    X
}

.quadform.sample.parameter <- function(n, domain.radius, sample.method) {
    switch(
        sample.method,
        uniform.parameter.disk = .quadform.sample.uniform.parameter.disk(n, domain.radius),
        radial.center.parameter.disk = .quadform.sample.radial.parameter.disk(n, domain.radius, "center"),
        radial.boundary.parameter.disk = .quadform.sample.radial.parameter.disk(n, domain.radius, "boundary"),
        uniform.parameter.square = .quadform.sample.uniform.parameter.square(n, domain.radius),
        radial.center.parameter.square = .quadform.sample.radial.parameter.square(n, domain.radius, "center"),
        radial.boundary.parameter.square = .quadform.sample.radial.parameter.square(n, domain.radius, "boundary"),
        stop("Unsupported sample method: ", sample.method, call. = FALSE)
    )
}

#' Estimate Reference Geodesic Distances on a 2D Quadratic Surface
#'
#' @description
#' Computes a numerical reference geodesic distance matrix for points sampled
#' from a two-dimensional quadratic graph surface. The method builds a dense
#' reference graph on a regular parameter-domain grid, includes the sample
#' points as vertices, weights all reference edges by exact quadratic-surface
#' segment length, and returns shortest-path distances between the sample
#' vertices.
#'
#' @param X Numeric matrix or data frame with two parameter-coordinate columns.
#' @param index.k Integer between 0 and 2. Number of positive-square terms.
#' @param domain.radius Positive numeric scalar or `NULL`. Radius of the
#'   parameter disk or half-side length of the parameter square covered by the
#'   reference grid. If `NULL`, uses the smallest radius/half-side covering the
#'   sample.
#' @param domain.shape Character scalar, either `"disk"` or `"square"`.
#' @param grid.size Odd or even integer at least 5. Number of grid coordinates
#'   along each axis before optional clipping to the disk.
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
                                         coefficients = NULL,
                                         domain.radius = NULL,
                                         domain.shape = c("disk", "square"),
                                         grid.size = 51,
                                         sample.connection.k = 8) {
    X <- .validate.numeric.data.matrix(X)
    if (ncol(X) != 2L) {
        stop("'quadform.reference.geodesics' currently supports only 2D parameter coordinates.",
             call. = FALSE)
    }
    .validate.quadform.index(2L, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, 2L)
    domain.shape <- .validate.quadform.domain.shape(domain.shape)
    if (is.null(domain.radius)) {
        domain.radius <- .quadform.default.domain.radius(X, domain.shape)
    }
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (any(!.quadform.points.inside.domain(X, domain.radius, domain.shape))) {
        stop("All rows of 'X' must lie inside the ",
             .quadform.domain.label(domain.shape), ".", call. = FALSE)
    }
    if (!is.numeric(grid.size) || length(grid.size) != 1L ||
        !is.finite(grid.size) || grid.size != floor(grid.size) || grid.size < 5L) {
        stop("'grid.size' must be an integer at least 5.", call. = FALSE)
    }
    grid.size <- as.integer(grid.size)
    ref <- .quadform.reference.vertices.2d(X, domain.radius, grid.size,
                                           domain.shape = domain.shape)
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
        quadform.edge.length(ref$vertices[e[[1L]], ], ref$vertices[e[[2L]], ],
                             index.k, coefficients = coefficients)
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
        domain_shape = domain.shape,
        grid_size = grid.size,
        sample_connection_k = sample.connection.k,
        index_k = as.integer(index.k),
        coefficients = coefficients
    )
}

#' Estimate 2D Quadratic-Surface Geodesics with the C++ Grid Engine
#'
#' @description
#' Computes numerical reference geodesic distances for points on a
#' two-dimensional quadratic graph surface using the package's C++ grid
#' shortest-path engine. A regular parameter grid is built once, optionally
#' clipped to the disk domain, grid edges are weighted by exact
#' quadratic-surface segment length, and each sample point is connected to
#' nearby grid vertices before running Dijkstra shortest paths.
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
                                             coefficients = NULL,
                                             domain.radius = NULL,
                                             domain.shape = c("disk", "square"),
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
    coefficients <- .validate.quadform.coefficients(coefficients, 2L)
    domain.shape <- .validate.quadform.domain.shape(domain.shape)
    if (is.null(domain.radius)) {
        domain.radius <- .quadform.default.domain.radius(X, domain.shape)
    }
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    if (any(!.quadform.points.inside.domain(X, domain.radius, domain.shape))) {
        stop("All rows of 'X' must lie inside the ",
             .quadform.domain.label(domain.shape), ".", call. = FALSE)
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
            X, index.k, as.integer(oracle.tube.k),
            coefficients = coefficients
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
        as.numeric(coefficients),
        domain.shape,
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

.quadform.default.oracle.tube.radius <- function(X, index.k, oracle.tube.k,
                                                coefficients = NULL) {
    n <- nrow(X)
    if (n < 2L) {
        stop("At least two sample points are required for sample-path oracle distances.",
             call. = FALSE)
    }
    k <- min(as.integer(oracle.tube.k), n - 1L)
    X.embed <- quadform.embed(X, index.k = index.k, coefficients = coefficients)
    kth <- numeric(n)
    for (i in seq_len(n)) {
        d <- sqrt(rowSums((t(t(X.embed) - X.embed[i, ]))^2))
        d[[i]] <- Inf
        kth[[i]] <- sort(d, partial = k)[[k]]
    }
    as.numeric(stats::median(kth))
}

.quadform.reference.grid.2d <- function(domain.radius, grid.size,
                                        domain.shape = c("disk", "square")) {
    domain.shape <- .validate.quadform.domain.shape(domain.shape)
    coords <- seq(-domain.radius, domain.radius, length.out = grid.size)
    grid <- expand.grid(x1 = coords, x2 = coords)
    X <- as.matrix(grid)
    if (identical(domain.shape, "disk")) {
        inside <- rowSums(X^2) <= domain.radius^2 * (1 + 1e-12)
    } else {
        inside <- rep(TRUE, nrow(X))
    }
    X <- X[inside, , drop = FALSE]
    storage.mode(X) <- "double"
    X
}

.sample.quadform.reference.pairs <- function(domain.radius,
                                             domain.shape,
                                             reference.grid.size,
                                             n.sources,
                                             targets.per.source,
                                             seed) {
    grid <- .quadform.reference.grid.2d(domain.radius, reference.grid.size,
                                        domain.shape = domain.shape)
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
#' @param domain.shape Character scalar, either `"disk"` or `"square"`.
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
                                               coefficients = NULL,
                                               domain.radius = 1,
                                               domain.shape = c("disk", "square"),
                                               grid.sizes = c(101, 201, 251, 501),
                                               reference.grid.size = 1001,
                                               n.sources = 50,
                                               targets.per.source = 20,
                                               seed = NULL) {
    .validate.quadform.index(2L, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, 2L)
    domain.shape <- .validate.quadform.domain.shape(domain.shape)
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
        domain.shape = domain.shape,
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
                    as.numeric(coefficients),
                    domain.shape,
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
            coefficients = coefficients,
            domain_radius = as.numeric(domain.radius),
            domain_shape = domain.shape,
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

.quadform.sample.parameter.3d <- function(n, domain.radius, domain.shape) {
    if (identical(domain.shape, "ball")) {
        dirs <- matrix(stats::rnorm(n * 3L), ncol = 3L)
        norms <- sqrt(rowSums(dirs^2))
        while (any(norms == 0)) {
            zero <- which(norms == 0)
            dirs[zero, ] <- matrix(stats::rnorm(length(zero) * 3L), ncol = 3L)
            norms <- sqrt(rowSums(dirs^2))
        }
        dirs <- dirs / norms
        radii <- domain.radius * stats::runif(n)^(1 / 3)
        return(dirs * radii)
    }
    matrix(stats::runif(n * 3L, -domain.radius, domain.radius), ncol = 3L)
}

.quadform.sample.boundary.3d <- function(n, domain.radius, domain.shape) {
    if (n <= 0L) {
        return(matrix(numeric(), ncol = 3L))
    }
    if (identical(domain.shape, "ball")) {
        dirs <- matrix(stats::rnorm(n * 3L), ncol = 3L)
        norms <- sqrt(rowSums(dirs^2))
        while (any(norms == 0)) {
            zero <- which(norms == 0)
            dirs[zero, ] <- matrix(stats::rnorm(length(zero) * 3L), ncol = 3L)
            norms <- sqrt(rowSums(dirs^2))
        }
        return((dirs / norms) * domain.radius)
    }

    out <- matrix(stats::runif(n * 3L, -domain.radius, domain.radius), ncol = 3L)
    face <- sample.int(6L, n, replace = TRUE)
    axis <- ((face - 1L) %% 3L) + 1L
    sign <- ifelse(face <= 3L, -1, 1)
    out[cbind(seq_len(n), axis)] <- sign * domain.radius
    if (n >= 8L) {
        corners <- as.matrix(expand.grid(
            x1 = c(-domain.radius, domain.radius),
            x2 = c(-domain.radius, domain.radius),
            x3 = c(-domain.radius, domain.radius)
        ))
        out[seq_len(8L), ] <- corners
    }
    out
}

.quadform.domain.volume.3d <- function(domain.radius, domain.shape) {
    if (identical(domain.shape, "ball")) {
        4 * pi * domain.radius^3 / 3
    } else {
        (2 * domain.radius)^3
    }
}

.quadform.cell.key <- function(x, epsilon) {
    paste(floor(x / epsilon), collapse = ",")
}

.quadform.greedy.epsilon.net <- function(candidates, epsilon, forced.n = 0L) {
    if (!is.numeric(epsilon) || length(epsilon) != 1L ||
        !is.finite(epsilon) || epsilon <= 0) {
        stop("'epsilon' must be a positive finite numeric scalar.", call. = FALSE)
    }
    candidates <- .validate.quadform.data.matrix(candidates)
    p <- ncol(candidates)
    forced.n <- as.integer(forced.n)
    if (forced.n < 0L || forced.n > nrow(candidates)) {
        stop("'forced.n' must be between 0 and nrow(candidates).", call. = FALSE)
    }
    neighbor.offsets <- as.matrix(expand.grid(rep(list(-1L:1L), p)))
    selected <- integer()
    grid <- new.env(parent = emptyenv(), hash = TRUE)
    eps2 <- epsilon^2

    insert <- function(i) {
        key <- .quadform.cell.key(candidates[i, ], epsilon)
        old <- grid[[key]]
        grid[[key]] <- c(old, i)
        selected <<- c(selected, i)
    }
    has.close <- function(i) {
        cell <- floor(candidates[i, ] / epsilon)
        for (r in seq_len(nrow(neighbor.offsets))) {
            key <- paste(cell + neighbor.offsets[r, ], collapse = ",")
            idx <- grid[[key]]
            if (length(idx)) {
                d2 <- rowSums((candidates[idx, , drop = FALSE] -
                               matrix(candidates[i, ], nrow = length(idx),
                                      ncol = p, byrow = TRUE))^2)
                if (any(d2 < eps2 * (1 - 1e-12))) {
                    return(TRUE)
                }
            }
        }
        FALSE
    }

    if (forced.n > 0L) {
        for (i in seq_len(forced.n)) {
            insert(i)
        }
    }
    if (forced.n < nrow(candidates)) {
        for (i in seq.int(forced.n + 1L, nrow(candidates))) {
            if (!has.close(i)) {
                insert(i)
            }
        }
    }
    candidates[selected, , drop = FALSE]
}

.quadform.epsilon.net.3d <- function(X,
                                     domain.radius,
                                     domain.shape,
                                     n.ref,
                                     seed = NULL,
                                     candidate.multiplier = 6,
                                     boundary.fraction = 0.2,
                                     epsilon = NULL) {
    X <- .validate.quadform.data.matrix(X)
    if (ncol(X) != 3L) {
        stop("'X' must have three parameter-coordinate columns.", call. = FALSE)
    }
    if (!is.numeric(n.ref) || length(n.ref) != 1L || !is.finite(n.ref) ||
        n.ref != floor(n.ref) || n.ref < nrow(X)) {
        stop("'n.ref' must be an integer at least nrow(X).", call. = FALSE)
    }
    if (!is.numeric(candidate.multiplier) || length(candidate.multiplier) != 1L ||
        !is.finite(candidate.multiplier) || candidate.multiplier < 1) {
        stop("'candidate.multiplier' must be a finite scalar at least 1.",
             call. = FALSE)
    }
    if (!is.numeric(boundary.fraction) || length(boundary.fraction) != 1L ||
        !is.finite(boundary.fraction) || boundary.fraction < 0 ||
        boundary.fraction > 0.9) {
        stop("'boundary.fraction' must be a finite scalar in [0, 0.9].",
             call. = FALSE)
    }
    if (!all(.quadform.points.inside.domain.nd(X, domain.radius, domain.shape))) {
        stop("All rows of 'X' must lie inside the parameter domain.",
             call. = FALSE)
    }
    if (is.null(epsilon)) {
        volume <- .quadform.domain.volume.3d(domain.radius, domain.shape)
        epsilon <- 0.62 * (volume / n.ref)^(1 / 3)
    } else if (!is.numeric(epsilon) || length(epsilon) != 1L ||
               !is.finite(epsilon) || epsilon <= 0) {
        stop("'epsilon' must be NULL or a positive finite numeric scalar.",
             call. = FALSE)
    }
    n.boundary <- as.integer(ceiling(n.ref * boundary.fraction))
    n.interior <- as.integer(ceiling(n.ref * candidate.multiplier))
    candidates <- .with.quadform.seed(seed, {
        rbind(
            X,
            .quadform.sample.boundary.3d(n.boundary, domain.radius, domain.shape),
            .quadform.sample.parameter.3d(n.interior, domain.radius, domain.shape)
        )
    })
    vertices <- .quadform.greedy.epsilon.net(candidates, epsilon,
                                             forced.n = nrow(X))
    rownames(vertices) <- NULL
    list(
        vertices = vertices,
        sample_vertex = seq_len(nrow(X)),
        epsilon = epsilon,
        n_target = as.integer(n.ref),
        n_candidates = nrow(candidates),
        n_boundary_candidates = n.boundary,
        n_vertices = nrow(vertices),
        domain_shape = domain.shape,
        domain_radius = domain.radius,
        seed = seed
    )
}

.quadform.delaunay.edges.3d <- function(vertices) {
    if (!requireNamespace("geometry", quietly = TRUE)) {
        stop("The optional 'geometry' package is required for 3D Delaunay reference geodesics.",
             call. = FALSE)
    }
    tess <- geometry::delaunayn(vertices, options = "Qt Qbb Qc")
    if (is.null(dim(tess))) {
        tess <- matrix(tess, nrow = 1L)
    }
    if (ncol(tess) != 4L) {
        stop("Delaunay tessellation did not return tetrahedra.", call. = FALSE)
    }
    edges <- rbind(
        tess[, c(1L, 2L), drop = FALSE],
        tess[, c(1L, 3L), drop = FALSE],
        tess[, c(1L, 4L), drop = FALSE],
        tess[, c(2L, 3L), drop = FALSE],
        tess[, c(2L, 4L), drop = FALSE],
        tess[, c(3L, 4L), drop = FALSE]
    )
    edges <- t(apply(edges, 1L, sort))
    edges <- unique(edges)
    storage.mode(edges) <- "integer"
    edges
}

.quadform.delaunay.edges.3d.cpp <- function(vertices,
                                            qhull.options = "Qt Qbb Qc") {
    vertices <- .validate.quadform.data.matrix(vertices)
    if (ncol(vertices) != 3L) {
        stop("'vertices' must have three columns.", call. = FALSE)
    }
    out <- rcpp_quadform_delaunay_edges_3d(vertices, qhull.options)
    edge.matrix <- out$edge_matrix
    storage.mode(edge.matrix) <- "integer"
    edge.matrix
}

.quadform.edge.list.to.adj <- function(n.vertices, edge.matrix, edge.weight) {
    adj <- vector("list", n.vertices)
    weight <- vector("list", n.vertices)
    for (i in seq_len(n.vertices)) {
        adj[[i]] <- integer()
        weight[[i]] <- numeric()
    }
    if (nrow(edge.matrix)) {
        for (r in seq_len(nrow(edge.matrix))) {
            a <- edge.matrix[r, 1L]
            b <- edge.matrix[r, 2L]
            w <- edge.weight[r]
            adj[[a]] <- c(adj[[a]], b)
            weight[[a]] <- c(weight[[a]], w)
            adj[[b]] <- c(adj[[b]], a)
            weight[[b]] <- c(weight[[b]], w)
        }
    }
    list(adj_list = adj, weight_list = weight)
}

.quadform.adj.components <- function(adj.list) {
    n <- length(adj.list)
    comp <- rep.int(0L, n)
    current <- 0L
    for (i in seq_len(n)) {
        if (comp[i] != 0L) {
            next
        }
        current <- current + 1L
        queue <- i
        comp[i] <- current
        head <- 1L
        while (head <= length(queue)) {
            v <- queue[head]
            head <- head + 1L
            nbr <- adj.list[[v]]
            unseen <- nbr[comp[nbr] == 0L]
            if (length(unseen)) {
                comp[unseen] <- current
                queue <- c(queue, unseen)
            }
        }
    }
    comp
}

.quadform.filter.delaunay.edges <- function(vertices,
                                            edge.matrix,
                                            edge.length.factor,
                                            epsilon) {
    if (!nrow(edge.matrix)) {
        attempts <- data.frame(
            filter_factor = edge.length.factor,
            filter_factor_label = ifelse(is.infinite(edge.length.factor), "Inf",
                                         as.character(edge.length.factor)),
            n_edges = 0L,
            retained_edge_fraction = NA_real_,
            n_components = nrow(vertices),
            connected = nrow(vertices) <= 1L,
            stringsAsFactors = FALSE
        )
        return(list(edge_matrix = edge.matrix,
                    filter_factor_used = edge.length.factor,
                    n_components = nrow(vertices),
                    filter_attempts = attempts,
                    filter_relaxation_happened = FALSE,
                    relaxed_to_inf = FALSE,
                    retained_edge_fraction = NA_real_))
    }
    parameter.length <- sqrt(rowSums(
        (vertices[edge.matrix[, 1L], , drop = FALSE] -
         vertices[edge.matrix[, 2L], , drop = FALSE])^2
    ))
    factors <- unique(c(edge.length.factor,
                        edge.length.factor * 1.5,
                        edge.length.factor * 2,
                        Inf))
    best <- NULL
    attempts <- vector("list", length(factors))
    for (factor in factors) {
        keep <- if (is.infinite(factor)) {
            rep(TRUE, nrow(edge.matrix))
        } else {
            parameter.length <= factor * epsilon
        }
        candidate.edges <- edge.matrix[keep, , drop = FALSE]
        graph <- .quadform.edge.list.to.adj(
            nrow(vertices), candidate.edges, rep(1, nrow(candidate.edges))
        )
        comp <- .quadform.adj.components(graph$adj_list)
        attempts[[which(factors == factor)[1L]]] <- data.frame(
            filter_factor = factor,
            filter_factor_label = ifelse(is.infinite(factor), "Inf",
                                         as.character(factor)),
            n_edges = nrow(candidate.edges),
            retained_edge_fraction = nrow(candidate.edges) / nrow(edge.matrix),
            n_components = max(comp),
            connected = max(comp) == 1L,
            stringsAsFactors = FALSE
        )
        best <- list(edge_matrix = candidate.edges,
                     filter_factor_used = factor,
                     n_components = max(comp))
        if (best$n_components == 1L) {
            break
        }
    }
    attempts <- do.call(rbind, Filter(Negate(is.null), attempts))
    relaxed <- !(identical(is.infinite(edge.length.factor),
                           is.infinite(best$filter_factor_used)) &&
                 (is.infinite(edge.length.factor) ||
                  isTRUE(all.equal(edge.length.factor, best$filter_factor_used))))
    best$filter_attempts <- attempts
    best$filter_relaxation_happened <- relaxed
    best$relaxed_to_inf <- isTRUE(is.finite(edge.length.factor) &&
                                  is.infinite(best$filter_factor_used))
    best$retained_edge_fraction <- nrow(best$edge_matrix) / nrow(edge.matrix)
    best
}

#' Compute 3D Quadratic-Hypersurface Reference Geodesics by Delaunay Graphs
#'
#' @description
#' Builds an experimental three-dimensional parameter-domain reference graph for
#' a quadratic graph hypersurface. The reference vertices are an approximate
#' epsilon-net that forcibly includes the supplied sample points. The one-
#' skeleton of the 3D Delaunay tessellation is weighted by exact quadratic-
#' hypersurface segment lengths, and shortest-path distances between sample
#' vertices are returned.
#'
#' @details
#' This function is intended as the first high-dimensional reference-geodesic
#' oracle for the data-geodesic reconstruction experiments. It currently uses
#' `geometry::delaunayn()` as a Qhull-backed Delaunay prototype. The exact edge
#' lengths are computed by a C++ vectorization of `quadform.edge.length()`. Long
#' Delaunay edges can be removed with `edge.length.factor`; if that disconnects
#' the reference graph, the factor is progressively relaxed, falling back to the
#' unfiltered Delaunay one-skeleton.
#'
#' @param X Numeric matrix or data frame with three-dimensional parameter
#'   points in rows.
#' @inheritParams quadform.embed
#' @param domain.radius Positive numeric scalar. Radius of the parameter ball or
#'   half-side length of the parameter cube. If `NULL`, it is inferred from `X`.
#' @param domain.shape Parameter-domain shape, either `"ball"` or `"cube"`.
#' @param n.ref Target number of reference vertices, including sample vertices.
#'   The epsilon-net is approximate, so the actual number is reported.
#' @param seed `NULL` or finite integer scalar for reproducible candidate
#'   sampling.
#' @param candidate.multiplier Number of interior candidate points per target
#'   reference vertex before greedy epsilon-net filtering.
#' @param boundary.fraction Approximate fraction of target reference vertices
#'   used as boundary candidates before filtering.
#' @param epsilon `NULL` or positive scalar. If supplied, this separation radius
#'   is used directly instead of estimating one from `n.ref`.
#' @param edge.length.factor Positive scalar or `Inf` used to remove Delaunay
#'   edges whose parameter-space length exceeds `edge.length.factor * epsilon`.
#'   The default, `4`, is the current provisional reference-oracle setting from
#'   the 3D Delaunay stress tests. Use `Inf` to disable long-edge filtering.
#'   The filter is relaxed automatically if needed for connectivity.
#'
#' @return A list of class `"quadform_delaunay_geodesics"` containing the sample
#'   distance matrix, reference vertices, Delaunay edge set, edge weights, and
#'   diagnostics.
#'
#' @examples
#' if (requireNamespace("geometry", quietly = TRUE)) {
#'   X <- matrix(c(0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5), ncol = 3,
#'               byrow = TRUE)
#'   ref <- quadform.delaunay.geodesic.distances(X, index.k = 3, n.ref = 80,
#'                                               domain.shape = "cube", seed = 1)
#'   dim(ref$distances)
#' }
#'
#' @export
quadform.delaunay.geodesic.distances <- function(X,
                                                 index.k,
                                                 coefficients = NULL,
                                                 domain.radius = NULL,
                                                 domain.shape = c("ball", "cube"),
                                                 n.ref = 5000,
                                                 seed = NULL,
                                                 candidate.multiplier = 6,
                                                 boundary.fraction = 0.2,
                                                 epsilon = NULL,
                                                 edge.length.factor = 4) {
    X <- .validate.quadform.data.matrix(X)
    if (ncol(X) != 3L) {
        stop("'X' must have three parameter-coordinate columns.", call. = FALSE)
    }
    index.k <- .validate.quadform.index(3L, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, 3L)
    domain.shape <- match.arg(domain.shape)
    if (is.null(domain.radius)) {
        domain.radius <- if (identical(domain.shape, "ball")) {
            max(sqrt(rowSums(X^2)))
        } else {
            max(abs(X))
        }
        if (domain.radius == 0) {
            domain.radius <- 1
        }
    }
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.",
             call. = FALSE)
    }
    if (!is.numeric(edge.length.factor) || length(edge.length.factor) != 1L ||
        is.na(edge.length.factor) || edge.length.factor <= 0) {
        stop("'edge.length.factor' must be a positive numeric scalar or Inf.",
             call. = FALSE)
    }

    net <- .quadform.epsilon.net.3d(
        X = X,
        domain.radius = domain.radius,
        domain.shape = domain.shape,
        n.ref = n.ref,
        seed = seed,
        candidate.multiplier = candidate.multiplier,
        boundary.fraction = boundary.fraction,
        epsilon = epsilon
    )
    all.edges <- .quadform.delaunay.edges.3d(net$vertices)
    filtered <- .quadform.filter.delaunay.edges(
        net$vertices, all.edges, edge.length.factor, net$epsilon
    )
    edge.matrix <- filtered$edge_matrix
    edge.weight <- quadform.edge.lengths(
        net$vertices[edge.matrix[, 1L], , drop = FALSE],
        net$vertices[edge.matrix[, 2L], , drop = FALSE],
        index.k = index.k,
        coefficients = coefficients
    )
    graph <- .quadform.edge.list.to.adj(nrow(net$vertices), edge.matrix,
                                        edge.weight)
    D <- shortest.path(graph$adj_list, graph$weight_list, net$sample_vertex)

    out <- list(
        distances = D,
        vertices_param = net$vertices,
        vertices_embed = quadform.embed(net$vertices, index.k = index.k,
                                        coefficients = coefficients),
        sample_vertex = net$sample_vertex,
        edge_matrix = edge.matrix,
        edge_weight = edge.weight,
        adj_list = graph$adj_list,
        weight_list = graph$weight_list,
        n_sample_vertices = nrow(X),
        n_reference_vertices = nrow(net$vertices),
        n_edges = nrow(edge.matrix),
        n_delaunay_edges = nrow(all.edges),
        n_edges_unfiltered = nrow(all.edges),
        retained_edge_fraction = filtered$retained_edge_fraction,
        epsilon = net$epsilon,
        n_target = net$n_target,
        n_candidates = net$n_candidates,
        n_boundary_candidates = net$n_boundary_candidates,
        filter_factor_requested = edge.length.factor,
        filter_factor_used = filtered$filter_factor_used,
        filter_attempts = filtered$filter_attempts,
        filter_relaxation_happened = filtered$filter_relaxation_happened,
        relaxed_to_inf = filtered$relaxed_to_inf,
        n_components = filtered$n_components,
        index_k = as.integer(index.k),
        coefficients = coefficients,
        domain_shape = domain.shape,
        domain_radius = as.numeric(domain.radius),
        seed = seed
    )
    class(out) <- c("quadform_delaunay_geodesics", "list")
    out
}

#' Sample a 2D Quadratic Surface Dataset with Reference Geodesics
#'
#' @description
#' Samples parameter points from either the disk
#' \eqn{\{x \in \mathbb{R}^2 : \|x\|_2 \le r\}} or the square
#' \eqn{[-r,r]^2} and embeds them into the quadratic graph surface
#' \deqn{q(x)=\sum_{i=1}^k x_i^2-\sum_{i=k+1}^2 x_i^2.}
#' The sampling methods are explicitly named by domain and density pattern:
#' `"uniform.parameter.disk"`, `"radial.center.parameter.disk"`,
#' `"radial.boundary.parameter.disk"`, `"uniform.parameter.square"`,
#' `"radial.center.parameter.square"`, and
#' `"radial.boundary.parameter.square"`. Points are sampled in parameter
#' coordinates, not with respect to the induced surface-area measure. The
#' function also builds the regular parameter-domain reference grid used by
#' `quadform.reference.geodesics()` and returns sample-by-sample reference
#' geodesic distances.
#'
#' @param n Positive integer. Number of sample points.
#' @param index.k Integer between 0 and 2. Number of positive-square terms in
#'   the quadratic form.
#' @param domain.radius Positive numeric scalar. Radius of the parameter disk
#'   or half-side length of the parameter square.
#' @param sample.method Character scalar naming the parameter-domain sampling
#'   method.
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
                                    coefficients = NULL,
                                    domain.radius = 1,
                                    sample.method = c(
                                        "uniform.parameter.disk",
                                        "radial.center.parameter.disk",
                                        "radial.boundary.parameter.disk",
                                        "uniform.parameter.square",
                                        "radial.center.parameter.square",
                                        "radial.boundary.parameter.square"
                                    ),
                                    grid.size = 51,
                                    sample.connection.k = 8,
                                    seed = NULL) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n != floor(n) || n < 1L) {
        stop("'n' must be a positive integer.", call. = FALSE)
    }
    n <- as.integer(n)
    .validate.quadform.index(2L, index.k)
    coefficients <- .validate.quadform.coefficients(coefficients, 2L)
    if (!is.numeric(domain.radius) || length(domain.radius) != 1L ||
        !is.finite(domain.radius) || domain.radius <= 0) {
        stop("'domain.radius' must be a positive finite numeric scalar.", call. = FALSE)
    }
    sample.method <- match.arg(sample.method)
    domain.shape <- if (grepl("\\.square$", sample.method)) "square" else "disk"

    X.param <- .with.quadform.seed(
        seed,
        .quadform.sample.parameter(n, domain.radius, sample.method)
    )
    X.embed <- quadform.embed(X.param, index.k = index.k,
                              coefficients = coefficients)
    ref <- quadform.reference.geodesics(
        X.param,
        index.k = index.k,
        coefficients = coefficients,
        domain.radius = domain.radius,
        domain.shape = domain.shape,
        grid.size = grid.size,
        sample.connection.k = sample.connection.k
    )
    grid.idx <- seq.int(n + 1L, nrow(ref$vertices))
    grid.param <- ref$vertices[grid.idx, , drop = FALSE]
    vertices.embed <- quadform.embed(ref$vertices, index.k = index.k,
                                     coefficients = coefficients)

    out <- list(
        X_param = X.param,
        X_embed = X.embed,
        q = as.numeric(X.embed[, "q"]),
        D_geodesic = ref$distances,
        distances = ref$distances,
        reference = list(
            grid_param = grid.param,
            grid_embed = quadform.embed(grid.param, index.k = index.k,
                                        coefficients = coefficients),
            vertices_param = ref$vertices,
            vertices_embed = vertices.embed,
            sample_vertex = seq_len(n),
            grid_vertex = grid.idx,
            edge_matrix = ref$edge_matrix,
            edge_weight = ref$edge_weight,
            n_reference_vertices = ref$n_reference_vertices,
            n_edges = ref$n_edges,
            domain_shape = ref$domain_shape,
            grid_size = ref$grid_size,
            sample_connection_k = ref$sample_connection_k
        ),
        metadata = list(
            dim = 2L,
            index_k = as.integer(index.k),
            coefficients = coefficients,
            domain_radius = as.numeric(domain.radius),
            domain_shape = domain.shape,
            sample_method = sample.method,
            grid_size = as.integer(ref$grid_size),
            sample_connection_k = as.integer(ref$sample_connection_k),
            seed = seed
        )
    )
    class(out) <- c("quadform_sample_dataset", "list")
    out
}
