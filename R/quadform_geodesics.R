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
