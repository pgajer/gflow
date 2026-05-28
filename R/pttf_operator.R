#' Construct PTTF Transported-Difference Operators
#'
#' Assembles experimental parallel-transport trend-filtering (PTTF)
#' transported-difference operators from a \code{\link{pttf.geometry}} object.
#' This is an operator-construction function only: it does not fit a response,
#' tune a penalty parameter, or run an \eqn{\ell_1}/\eqn{\ell_2} smoother.
#'
#' @param geometry A \code{"pttf_geometry"} object from \code{\link{pttf.geometry}}.
#' @param derivative.order Integer derivative order, one of \code{1L},
#'   \code{2L}, or \code{3L}.
#' @param edge.status.policy Edge transport policy. \code{"ok.only"} uses only
#'   Phase 1 transports with status \code{"ok"}. \code{"frame.fallback"} also
#'   allows direct-frame fallback transports for non-OK shared-support statuses.
#' @param regression.weight.rule Local regression edge-weight rule.
#' @param edge.length.epsilon Positive numeric edge-length floor.
#' @param tensor.scaling Symmetric tensor coordinate scaling. \code{"hs"} uses
#'   square-root multiplicity scaling so coordinate Euclidean norms match
#'   Hilbert--Schmidt norms. \code{"raw"} stores unscaled unique symmetric
#'   components.
#' @param row.mass.rule Final row-mass multiplier. \code{"node.mass"} multiplies
#'   rows for vertex \eqn{i} by \eqn{\sqrt{\mu_i}} from the Phase 1 density
#'   metadata. \code{"none"} leaves rows unweighted.
#' @param row.normalize Optional final row normalization.
#' @param min.operator.rank.tol Relative singular-value tolerance for local
#'   derivative regression rank.
#' @param max.operator.condition Maximum allowed weighted local-design
#'   condition number.
#' @param return.intermediate Logical. If \code{TRUE}, include full logical
#'   intermediate derivative matrices.
#' @param return.A Logical. If \code{TRUE}, include the compact sparse row
#'   operator.
#' @param return.B Logical. If \code{TRUE}, include \eqn{A^\top A}.
#' @param diagnostics Logical. If \code{TRUE}, include assembly diagnostics.
#'
#' @details
#' Phase 2 keeps a full logical vertex/component layout internally. For order
#' \eqn{r}, the logical field has one block of
#' \eqn{q_r=\binom{m+r-1}{r}} components at each vertex. The returned
#' \code{A} is compact: it includes rows only for accepted final-order
#' vertex/component blocks, and \code{row.table$full.row} records the original
#' full logical row.
#'
#' For derivative level \eqn{r}, the local regression uses blocks
#' \eqn{\Phi_r(u_{ij})\in\mathbb R^{q_{r-1}\times q_r}},
#' lifted weights
#' \eqn{W_i^{(r)}=\operatorname{diag}(w_{ij})\otimes I_{q_{r-1}}},
#' and transported predecessor differences
#' \eqn{K_{r-1}(O_{ij})z_j-z_i}. Unavailable predecessor blocks are never filled
#' with zero; affected neighbors are dropped and recorded.
#'
#' @return A list of class \code{"pttf_operator"} containing sparse operators,
#'   row/vertex provenance, tensor-basis metadata, and diagnostics.
#' @export
pttf.operator <- function(
    geometry,
    derivative.order = 3L,
    edge.status.policy = c("ok.only", "frame.fallback"),
    regression.weight.rule = c("inverse.length.squared", "inverse.length", "unit"),
    edge.length.epsilon = 1e-8,
    tensor.scaling = c("hs", "raw"),
    row.mass.rule = c("node.mass", "none"),
    row.normalize = c("none", "l2"),
    min.operator.rank.tol = 1e-10,
    max.operator.condition = 1e8,
    return.intermediate = TRUE,
    return.A = TRUE,
    return.B = TRUE,
    diagnostics = TRUE) {

    if (!inherits(geometry, "pttf_geometry")) {
        stop("'geometry' must be a pttf_geometry object from pttf.geometry().",
             call. = FALSE)
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for pttf.operator().", call. = FALSE)
    }

    derivative.order <- .pttf.operator.match.order(derivative.order)
    edge.status.policy <- match.arg(edge.status.policy)
    regression.weight.rule <- match.arg(regression.weight.rule)
    tensor.scaling <- match.arg(tensor.scaling)
    row.mass.rule <- match.arg(row.mass.rule)
    row.normalize <- match.arg(row.normalize)
    .pttf.operator.require.positive(edge.length.epsilon, "edge.length.epsilon")
    .pttf.operator.require.positive(min.operator.rank.tol, "min.operator.rank.tol")
    .pttf.operator.require.positive(max.operator.condition, "max.operator.condition")

    n <- nrow(geometry$frames$centers)
    m <- geometry$frames$tangent.dim
    bases <- lapply(0:3, function(order) {
        .pttf.symmetric.tensor.basis(m, order, tensor.scaling)
    })
    names(bases) <- paste0("order", 0:3)

    oriented <- .pttf.operator.oriented.edges(
        geometry = geometry,
        edge.status.policy = edge.status.policy,
        regression.weight.rule = regression.weight.rule,
        edge.length.epsilon = edge.length.epsilon
    )
    state <- .pttf.operator.initial.state(n)
    level.results <- vector("list", derivative.order)
    vertex.tables <- list()

    for (level in seq_len(derivative.order)) {
        assembled <- .pttf.operator.assemble.level(
            geometry = geometry,
            oriented = oriented,
            bases = bases,
            level = level,
            prev.full = state$prev.full,
            available.prev = state$available.prev,
            min.operator.rank.tol = min.operator.rank.tol,
            max.operator.condition = max.operator.condition
        )
        level.results[[level]] <- assembled
        vertex.tables[[level]] <- assembled$vertex.table
        state$prev.full <- assembled$D.full
        state$available.prev <- assembled$available
    }

    final <- level.results[[derivative.order]]
    compact <- .pttf.operator.compact.final(
        D.full = final$D.full,
        available = final$available,
        bases = bases,
        derivative.order = derivative.order,
        geometry = geometry,
        row.mass.rule = row.mass.rule,
        row.normalize = row.normalize
    )
    A <- compact$A
    B <- if (isTRUE(return.B)) Matrix::crossprod(A) else NULL

    out <- list(
        A = if (isTRUE(return.A)) A else NULL,
        B = B,
        A.triplet = if (isTRUE(return.A)) .pttf.operator.triplet(A) else NULL,
        row.table = compact$row.table,
        vertex.table = .pttf.operator.vertex.table(vertex.tables, edge.status.policy),
        tensor.basis = bases,
        geometry.summary = .pttf.operator.geometry.summary(geometry),
        intermediate = if (isTRUE(return.intermediate)) {
            .pttf.operator.intermediate(level.results, derivative.order)
        } else NULL,
        diagnostics = if (isTRUE(diagnostics)) {
            .pttf.operator.diagnostics(
                level.results = level.results,
                final = final,
                compact = compact,
                oriented = oriented,
                derivative.order = derivative.order,
                n = n,
                bases = bases
            )
        } else NULL,
        parameters = list(
            derivative.order = derivative.order,
            edge.status.policy = edge.status.policy,
            regression.weight.rule = regression.weight.rule,
            edge.length.epsilon = edge.length.epsilon,
            tensor.scaling = tensor.scaling,
            row.mass.rule = row.mass.rule,
            row.normalize = row.normalize,
            min.operator.rank.tol = min.operator.rank.tol,
            max.operator.condition = max.operator.condition
        ),
        geometry = geometry
    )
    class(out) <- c("pttf_operator", "list")
    out
}

.pttf.operator.match.order <- function(x) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != as.integer(x) ||
        !(as.integer(x) %in% 1:3)) {
        stop("'derivative.order' must be one of 1L, 2L, or 3L.", call. = FALSE)
    }
    as.integer(x)
}

.pttf.operator.vertex.table <- function(vertex.tables, edge.status.policy) {
    out <- do.call(rbind, vertex.tables)
    out$edge.policy <- edge.status.policy
    row.names(out) <- NULL
    out
}

.pttf.operator.require.positive <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(sprintf("'%s' must be a positive numeric scalar.", name), call. = FALSE)
    }
    invisible(x)
}

.pttf.operator.initial.state <- function(n) {
    list(
        prev.full = Matrix::Diagonal(n),
        available.prev = rep(TRUE, n)
    )
}

.pttf.operator.q <- function(m, order) {
    if (order == 0L) return(1L)
    as.integer(choose(m + order - 1L, order))
}

.pttf.multi.indices <- function(m, order) {
    if (order == 0L) {
        return(matrix(integer(), nrow = 1L, ncol = 0L))
    }
    out <- utils::combn(seq_len(m + order - 1L), order)
    out <- t(out - matrix(seq_len(order) - 1L, nrow = order, ncol = ncol(out)))
    storage.mode(out) <- "integer"
    out
}

.pttf.index.multiplicity <- function(idx, m) {
    if (!length(idx)) return(1L)
    tab <- tabulate(idx, nbins = m)
    as.integer(factorial(length(idx)) / prod(factorial(tab)))
}

.pttf.symmetric.tensor.basis <- function(m, order, scaling = c("hs", "raw")) {
    scaling <- match.arg(scaling)
    multi <- .pttf.multi.indices(m, order)
    mult <- apply(multi, 1L, .pttf.index.multiplicity, m = m)
    scale <- if (identical(scaling, "hs")) sqrt(mult) else rep(1, length(mult))
    list(
        m = m,
        order = order,
        q = nrow(multi),
        multi.index = multi,
        multiplicity = mult,
        scale = scale,
        scaling = scaling
    )
}

.pttf.tensor.coord.from.full <- function(full, basis) {
    if (basis$order == 0L) return(as.numeric(full))
    out <- numeric(basis$q)
    for (b in seq_len(basis$q)) {
        idx <- as.list(basis$multi.index[b, ])
        out[b] <- basis$scale[b] * do.call(`[`, c(list(full), idx, list(drop = TRUE)))
    }
    out
}

.pttf.tensor.full.from.coord <- function(coord, basis) {
    if (basis$order == 0L) return(as.numeric(coord)[1L])
    dims <- rep(basis$m, basis$order)
    full <- array(0, dim = dims)
    for (b in seq_len(basis$q)) {
        idx <- basis$multi.index[b, ]
        raw <- coord[b] / basis$scale[b]
        perms <- unique(.pttf.permutations(idx))
        for (r in seq_len(nrow(perms))) {
            full <- do.call("[<-", c(list(full), as.list(perms[r, ]), list(raw)))
        }
    }
    full
}

.pttf.permutations <- function(x) {
    if (length(x) <= 1L) return(matrix(x, nrow = 1L))
    out <- do.call(rbind, lapply(seq_along(x), function(i) {
        cbind(x[i], .pttf.permutations(x[-i]))
    }))
    unique(out)
}

.pttf.tensor.contract.last <- function(full, u) {
    order <- length(dim(full))
    if (order == 0L) stop("Cannot contract an order-0 tensor.", call. = FALSE)
    if (order == 1L) return(sum(full * u))
    if (order == 2L) return(as.vector(full %*% u))
    if (order == 3L) {
        m <- dim(full)[1L]
        out <- matrix(0, nrow = m, ncol = m)
        for (a in seq_len(m)) {
            for (b in seq_len(m)) {
                out[a, b] <- sum(full[a, b, ] * u)
            }
        }
        return(out)
    }
    stop("Tensor contraction currently supports orders up to 3.", call. = FALSE)
}

.pttf.symmetric.tensor.direction.design <- function(u, input.order, output.order,
                                                    scaling = c("hs", "raw")) {
    scaling <- match.arg(scaling)
    r <- as.integer(output.order)
    if (input.order != r - 1L) {
        stop("'input.order' must equal output.order - 1.", call. = FALSE)
    }
    m <- length(u)
    in.basis <- .pttf.symmetric.tensor.basis(m, input.order, scaling)
    out.basis <- .pttf.symmetric.tensor.basis(m, output.order, scaling)
    Phi <- matrix(0, nrow = in.basis$q, ncol = out.basis$q)
    for (b in seq_len(out.basis$q)) {
        coord <- numeric(out.basis$q)
        coord[b] <- 1
        full <- .pttf.tensor.full.from.coord(coord, out.basis)
        contracted <- .pttf.tensor.contract.last(full, u)
        Phi[, b] <- .pttf.tensor.coord.from.full(contracted, in.basis)
    }
    Phi
}

.pttf.symmetric.tensor.transport <- function(O, order, scaling = c("hs", "raw")) {
    scaling <- match.arg(scaling)
    order <- as.integer(order)
    m <- nrow(O)
    basis <- .pttf.symmetric.tensor.basis(m, order, scaling)
    if (order == 0L) return(matrix(1, nrow = 1L, ncol = 1L))
    K <- matrix(0, nrow = basis$q, ncol = basis$q)
    for (b in seq_len(basis$q)) {
        coord <- numeric(basis$q)
        coord[b] <- 1
        full <- .pttf.tensor.full.from.coord(coord, basis)
        transported <- .pttf.tensor.transport.full(full, O, order)
        K[, b] <- .pttf.tensor.coord.from.full(transported, basis)
    }
    K
}

.pttf.tensor.transport.full <- function(full, O, order) {
    if (order == 1L) return(as.vector(O %*% full))
    if (order == 2L) return(O %*% full %*% t(O))
    if (order == 3L) {
        m <- nrow(O)
        out <- array(0, dim = c(m, m, m))
        for (a in seq_len(m)) {
            for (b in seq_len(m)) {
                for (c in seq_len(m)) {
                    total <- 0
                    for (p in seq_len(m)) {
                        for (q in seq_len(m)) {
                            for (r in seq_len(m)) {
                                total <- total + O[a, p] * O[b, q] *
                                    O[c, r] * full[p, q, r]
                            }
                        }
                    }
                    out[a, b, c] <- total
                }
            }
        }
        return(out)
    }
    stop("Tensor transport currently supports orders 1, 2, and 3.", call. = FALSE)
}

.pttf.operator.oriented.edges <- function(geometry, edge.status.policy,
                                          regression.weight.rule,
                                          edge.length.epsilon) {
    edges <- geometry$edges
    rows <- vector("list", 2L * nrow(edges))
    cursor <- 0L
    for (e in seq_len(nrow(edges))) {
        for (pair in list(c(edges$from[e], edges$to[e]),
                          c(edges$to[e], edges$from[e]))) {
            status <- edges$status[e]
            use <- identical(status, "ok")
            source <- if (use) "procrustes" else NA_character_
            if (!use && identical(edge.status.policy, "frame.fallback") &&
                status %in% c("too_few_shared_points", "shared_rank_deficient",
                              "shared_ill_conditioned")) {
                use <- TRUE
                source <- "frame.fallback"
            }
            ell <- edges$length[e]
            w <- switch(regression.weight.rule,
                        inverse.length.squared = (ell + edge.length.epsilon)^(-2),
                        inverse.length = (ell + edge.length.epsilon)^(-1),
                        unit = 1)
            cursor <- cursor + 1L
            rows[[cursor]] <- data.frame(
                from = pair[1L],
                to = pair[2L],
                edge.index = e,
                length = ell,
                regression.weight = w,
                use = use,
                transport.source = source,
                transport.status = status,
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows)) do.call(rbind, rows) else data.frame()
}

.pttf.operator.neighbors <- function(oriented, i) {
    rows <- oriented[oriented$from == i & oriented$use, , drop = FALSE]
    rows[order(rows$to), , drop = FALSE]
}

.pttf.operator.assemble.level <- function(geometry, oriented, bases, level,
                                          prev.full, available.prev,
                                          min.operator.rank.tol,
                                          max.operator.condition) {
    n <- nrow(geometry$frames$centers)
    q.prev <- bases[[level]]$q
    q.cur <- bases[[level + 1L]]$q
    full.rows <- n * q.cur
    D.full <- Matrix::Matrix(0, nrow = full.rows, ncol = n, sparse = TRUE)
    available <- rep(FALSE, n)
    vertex.rows <- vector("list", n)

    for (i in seq_len(n)) {
        assembled <- .pttf.operator.assemble.vertex(
            geometry = geometry,
            oriented = oriented,
            bases = bases,
            level = level,
            vertex = i,
            prev.full = prev.full,
            available.prev = available.prev,
            min.operator.rank.tol = min.operator.rank.tol,
            max.operator.condition = max.operator.condition
        )
        vertex.rows[[i]] <- assembled$vertex.row
        if (identical(assembled$status, "ok")) {
            available[i] <- TRUE
            rows <- ((i - 1L) * q.cur + 1L):(i * q.cur)
            D.full[rows, ] <- assembled$block
        }
    }

    list(
        D.full = D.full,
        available = available,
        vertex.table = do.call(rbind, vertex.rows),
        derivative.order = level
    )
}

.pttf.operator.assemble.vertex <- function(geometry, oriented, bases, level,
                                           vertex, prev.full, available.prev,
                                           min.operator.rank.tol,
                                           max.operator.condition) {
    n <- nrow(geometry$frames$centers)
    q.prev <- bases[[level]]$q
    q.cur <- bases[[level + 1L]]$q
    frame.status <- geometry$frames$support.diagnostics$status[vertex]
    candidates <- .pttf.operator.neighbors(oriented, vertex)
    dropped.predecessor <- 0L
    center.available <- TRUE
    status <- "ok"

    if (!identical(frame.status, "ok")) {
        status <- "geometry_support_not_ok"
    }
    if (level > 1L && !isTRUE(available.prev[vertex])) {
        center.available <- FALSE
        status <- "predecessor_block_unavailable"
    }
    if (!nrow(candidates) && identical(status, "ok")) {
        status <- "no_valid_transport"
    }

    if (identical(status, "ok") && level > 1L) {
        keep <- available.prev[candidates$to]
        dropped.predecessor <- sum(!keep)
        candidates <- candidates[keep, , drop = FALSE]
        if (!nrow(candidates)) {
            status <- "predecessor_block_unavailable"
        }
    }

    if (!identical(status, "ok")) {
        return(.pttf.operator.empty.vertex(level, vertex, nrow(candidates),
                                           dropped.predecessor, center.available,
                                           status, q.cur, n))
    }

    design <- .pttf.operator.local.design(
        geometry = geometry,
        candidates = candidates,
        vertex = vertex,
        bases = bases,
        level = level,
        n = n
    )
    sv <- svd(design$weighted.Phi, nu = 0L, nv = 0L)$d
    rank <- if (length(sv)) sum(sv > min.operator.rank.tol * max(sv)) else 0L
    condition <- if (rank > 0L && length(sv) >= q.cur) {
        max(sv) / max(sv[q.cur], .Machine$double.eps)
    } else Inf

    if (rank < q.cur) {
        status <- "operator_rank_deficient"
    } else if (condition > max.operator.condition) {
        status <- "operator_ill_conditioned"
    }
    if (!identical(status, "ok")) {
        return(.pttf.operator.empty.vertex(level, vertex, nrow(candidates),
                                           dropped.predecessor, center.available,
                                           status, q.cur, n, rank, condition))
    }

    pinv <- .pttf.operator.pinv(design$weighted.Phi, min.operator.rank.tol)
    G <- pinv %*% design$weighted.Delta
    block <- Matrix::Matrix(G %*% as.matrix(prev.full), sparse = TRUE)
    vertex.row <- .pttf.operator.vertex.row(
        level = level,
        vertex = vertex,
        candidate.degree = nrow(.pttf.operator.neighbors(oriented, vertex)),
        usable.degree = nrow(candidates),
        design.rank = rank,
        design.condition = condition,
        dropped.edge.count = sum(!oriented$use & oriented$from == vertex),
        n.predecessor.neighbors.dropped = dropped.predecessor,
        predecessor.center.available = center.available,
        status = "ok"
    )
    list(block = block, status = "ok", vertex.row = vertex.row)
}

.pttf.operator.empty.vertex <- function(level, vertex, usable.degree,
                                        dropped.predecessor, center.available,
                                        status, q.cur, n, rank = 0L,
                                        condition = Inf) {
    list(
        block = Matrix::Matrix(0, nrow = q.cur, ncol = n, sparse = TRUE),
        status = status,
        vertex.row = .pttf.operator.vertex.row(
            level = level,
            vertex = vertex,
            candidate.degree = usable.degree,
            usable.degree = 0L,
            design.rank = rank,
            design.condition = condition,
            dropped.edge.count = NA_integer_,
            n.predecessor.neighbors.dropped = dropped.predecessor,
            predecessor.center.available = center.available,
            status = status
        )
    )
}

.pttf.operator.vertex.row <- function(level, vertex, candidate.degree,
                                      usable.degree, design.rank,
                                      design.condition, dropped.edge.count,
                                      n.predecessor.neighbors.dropped,
                                      predecessor.center.available,
                                      status) {
    data.frame(
        vertex = vertex,
        derivative.order = level,
        candidate.degree = candidate.degree,
        usable.degree = usable.degree,
        design.rank = design.rank,
        design.condition = design.condition,
        dropped.edge.count = dropped.edge.count,
        n.predecessor.neighbors.dropped = n.predecessor.neighbors.dropped,
        predecessor.center.available = predecessor.center.available,
        edge.policy = NA_character_,
        status = status,
        stringsAsFactors = FALSE
    )
}

.pttf.operator.local.design <- function(geometry, candidates, vertex, bases,
                                        level, n) {
    q.prev <- bases[[level]]$q
    q.cur <- bases[[level + 1L]]$q
    n.blocks <- nrow(candidates)
    Phi <- matrix(0, nrow = n.blocks * q.prev, ncol = q.cur)
    Delta <- matrix(0, nrow = n.blocks * q.prev, ncol = n * q.prev)
    weights <- numeric(n.blocks * q.prev)

    for (k in seq_len(n.blocks)) {
        j <- candidates$to[k]
        key <- paste(vertex, j, sep = "->")
        u <- geometry$transport$uij[[key]]
        O <- geometry$transport$Oij[[key]]
        K <- .pttf.symmetric.tensor.transport(O, level - 1L,
                                               bases[[level]]$scaling)
        phi <- .pttf.symmetric.tensor.direction.design(
            u = u,
            input.order = level - 1L,
            output.order = level,
            scaling = bases[[level]]$scaling
        )
        rr <- ((k - 1L) * q.prev + 1L):(k * q.prev)
        Phi[rr, ] <- phi
        weights[rr] <- candidates$regression.weight[k]
        center.cols <- ((vertex - 1L) * q.prev + 1L):(vertex * q.prev)
        neighbor.cols <- ((j - 1L) * q.prev + 1L):(j * q.prev)
        Delta[rr, center.cols] <- Delta[rr, center.cols] - diag(1, q.prev)
        Delta[rr, neighbor.cols] <- Delta[rr, neighbor.cols] + K
    }
    sqrtw <- sqrt(weights)
    list(
        Phi = Phi,
        Delta = Delta,
        weighted.Phi = Phi * sqrtw,
        weighted.Delta = Delta * sqrtw,
        weights = weights
    )
}

.pttf.operator.pinv <- function(X, tol) {
    sv <- svd(X)
    if (!length(sv$d)) return(matrix(0, ncol(X), nrow(X)))
    keep <- sv$d > tol * max(sv$d)
    if (!any(keep)) return(matrix(0, ncol(X), nrow(X)))
    sv$v[, keep, drop = FALSE] %*%
        diag(1 / sv$d[keep], nrow = sum(keep)) %*%
        t(sv$u[, keep, drop = FALSE])
}

.pttf.operator.compact.final <- function(D.full, available, bases,
                                         derivative.order, geometry,
                                         row.mass.rule, row.normalize) {
    basis <- bases[[derivative.order + 1L]]
    q <- basis$q
    n <- length(available)
    rows <- unlist(lapply(which(available), function(i) {
        ((i - 1L) * q + 1L):(i * q)
    }), use.names = FALSE)
    if (!length(rows)) {
        A <- Matrix::Matrix(0, nrow = 0L, ncol = ncol(D.full), sparse = TRUE)
        return(list(A = A, row.table = data.frame()))
    }
    A <- D.full[rows, , drop = FALSE]
    row.table <- .pttf.operator.row.table(rows, basis, derivative.order, geometry)
    mass <- .pttf.operator.row.mass(row.table$vertex, geometry, row.mass.rule)
    norms.before <- sqrt(Matrix::rowSums(A^2))
    A <- Matrix::Diagonal(x = mass) %*% A
    if (identical(row.normalize, "l2") && nrow(A)) {
        norms <- sqrt(Matrix::rowSums(A^2))
        scale <- ifelse(norms > 0, 1 / norms, 1)
        A <- Matrix::Diagonal(x = scale) %*% A
    }
    norms.after <- sqrt(Matrix::rowSums(A^2))
    row.table$row.mass <- mass
    row.table$row.norm.before.normalization <- as.numeric(norms.before)
    row.table$row.norm.after.normalization <- as.numeric(norms.after)
    row.table$compact.row <- seq_len(nrow(row.table))
    list(A = A, row.table = row.table)
}

.pttf.operator.row.table <- function(rows, basis, derivative.order, geometry) {
    q <- basis$q
    vertex <- ((rows - 1L) %/% q) + 1L
    component <- ((rows - 1L) %% q) + 1L
    multi <- apply(basis$multi.index[component, , drop = FALSE], 1L,
                   function(x) paste(x, collapse = ","))
    data.frame(
        row = seq_along(rows),
        vertex = vertex,
        derivative.order = derivative.order,
        tensor.component = component,
        component = component,
        tensor.multi.index = multi,
        tensor.multiplicity = basis$multiplicity[component],
        tensor.scaling = basis$scaling,
        row.mass = NA_real_,
        row.norm.before.normalization = NA_real_,
        row.norm.after.normalization = NA_real_,
        status = "ok",
        full.row = rows,
        compact.row = seq_along(rows),
        stringsAsFactors = FALSE
    )
}

.pttf.operator.row.mass <- function(vertex, geometry, row.mass.rule) {
    if (identical(row.mass.rule, "none")) return(rep(1, length(vertex)))
    sqrt(pmax(geometry$density$node.mass[vertex], 0))
}

.pttf.operator.triplet <- function(A) {
    s <- methods::as(A, "TsparseMatrix")
    list(i = s@i + 1L, j = s@j + 1L, x = s@x, dim = dim(A))
}

.pttf.operator.geometry.summary <- function(geometry) {
    list(
        n.vertices = nrow(geometry$frames$centers),
        tangent.dim = geometry$frames$tangent.dim,
        n.edges = nrow(geometry$edges),
        geometry.class = class(geometry)[1L]
    )
}

.pttf.operator.intermediate <- function(level.results, derivative.order) {
    out <- vector("list", derivative.order)
    for (level in seq_len(derivative.order)) {
        out[[level]] <- list(
            D.full = level.results[[level]]$D.full,
            available = level.results[[level]]$available
        )
    }
    names(out) <- paste0("D", seq_len(derivative.order), ".full")
    out
}

.pttf.operator.diagnostics <- function(level.results, final, compact, oriented,
                                       derivative.order, n, bases) {
    availability <- lapply(level.results, `[[`, "available")
    names(availability) <- paste0("available.level", seq_along(availability))
    list(
        operator.dim = dim(compact$A),
        full.logical.dim = c(n * bases[[derivative.order + 1L]]$q, n),
        n.rows.by.level = vapply(level.results, function(x) {
            sum(x$available) * bases[[x$derivative.order + 1L]]$q
        }, numeric(1L)),
        availability = availability,
        dropped.edges.by.status = as.data.frame(table(
            oriented$transport.status[!oriented$use]
        )),
        row.norm.summary = summary(as.numeric(sqrt(Matrix::rowSums(compact$A^2)))),
        tensor.q = vapply(bases, `[[`, integer(1L), "q")
    )
}

#' @export
print.pttf_operator <- function(x, ...) {
    cat("PTTF transported-difference operator\n")
    cat("  derivative order:", x$parameters$derivative.order, "\n")
    if (!is.null(x$A)) {
        cat("  A dim:", paste(dim(x$A), collapse = " x "), "\n")
    } else {
        cat("  A dim: not returned\n")
    }
    invisible(x)
}

#' Filter Rows Of A PTTF Operator
#'
#' Creates a new \code{"pttf_operator"} object with a subset of compact
#' operator rows. The sparse triplet payload is rebuilt from the filtered
#' matrix so triplet row indices always refer to the filtered fit rows, while
#' original compact and full logical row identities remain in \code{row.table}.
#'
#' @param operator A \code{"pttf_operator"} object from \code{\link{pttf.operator}}.
#' @param rows Integer, logical, or character row selector. Character selectors
#'   are matched against \code{row.table$status}.
#' @param reason Optional character reason stored in \code{row.filter}.
#' @param preserve.original Logical. If \code{TRUE}, preserve original row
#'   identities in \code{row.table$compact.row} and \code{row.table$full.row}.
#'
#' @return A filtered \code{"pttf_operator"} object.
#' @export
pttf.operator.filter.rows <- function(operator,
                                      rows,
                                      reason = NULL,
                                      preserve.original = TRUE) {
    if (!inherits(operator, "pttf_operator")) {
        stop("'operator' must be a pttf_operator object.", call. = FALSE)
    }
    if (is.null(operator$A)) {
        stop("'operator' must contain A.", call. = FALSE)
    }
    if (is.null(operator$row.table)) {
        stop("'operator' must contain row.table.", call. = FALSE)
    }
    A <- methods::as(operator$A, "dgCMatrix")
    n.rows <- nrow(A)
    if (is.character(rows)) {
        if (!"status" %in% names(operator$row.table)) {
            stop("Character row selection requires row.table$status.",
                 call. = FALSE)
        }
        keep <- which(operator$row.table$status %in% rows)
    } else if (is.logical(rows)) {
        if (length(rows) != n.rows) {
            stop("Logical 'rows' must have length nrow(operator$A).",
                 call. = FALSE)
        }
        keep <- which(rows)
    } else if (is.numeric(rows) || is.integer(rows)) {
        keep <- as.integer(rows)
        if (anyNA(keep) || any(keep < 1L) || any(keep > n.rows)) {
            stop("Integer 'rows' must be valid row indices for operator$A.",
                 call. = FALSE)
        }
    } else {
        stop("'rows' must be an integer, logical, or character selector.",
             call. = FALSE)
    }
    keep <- unique(keep)

    out <- operator
    out$A <- A[keep, , drop = FALSE]
    out$B <- Matrix::crossprod(out$A)
    out$A.triplet <- .pttf.operator.triplet(out$A)
    row.table <- operator$row.table[keep, , drop = FALSE]
    row.names(row.table) <- NULL
    if (!isTRUE(preserve.original)) {
        row.table$compact.row <- seq_len(nrow(row.table))
    }
    row.table$fit.row <- seq_len(nrow(row.table))
    out$row.table <- row.table
    out$row.filter <- list(
        original.nrow = n.rows,
        retained.nrow = nrow(out$A),
        dropped.nrow = n.rows - nrow(out$A),
        retained.rows = keep,
        reason = reason %||% NA_character_,
        preserve.original = isTRUE(preserve.original)
    )
    if (!is.null(out$diagnostics)) {
        out$diagnostics$operator.dim <- dim(out$A)
        out$diagnostics$row.norm.summary <- summary(as.numeric(
            sqrt(Matrix::rowSums(out$A^2))
        ))
    }
    class(out) <- class(operator)
    out
}
