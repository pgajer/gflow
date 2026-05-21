#' Construct a Graph Trend-Filtering Operator
#'
#' Builds the weighted graph difference operators used by graph trend
#' filtering. The default \code{operator.family = "graph.laplacian.recursive"}
#' uses the phase-2 recursive graph operators
#' \deqn{k=0:\quad \Delta_w^{(1)} = D_w,}
#' \deqn{k=1:\quad \Delta_w^{(2)} = L_w = D_w^\top D_w,}
#' and
#' \deqn{k=2:\quad \Delta_w^{(3)} = D_w L_w.}
#' The experimental \code{operator.family = "path.divided.difference"} builds
#' canonical path rows and computes position-aware finite differences along
#' graph paths using \code{weight.list} as metric edge lengths.
#'
#' @param adj.list List of integer neighbor vectors using 1-based vertex
#'   indices. The graph must be undirected: every edge must appear in both
#'   endpoint adjacency lists.
#' @param weight.list Optional list of positive edge weights parallel to
#'   \code{adj.list}. For \code{operator.family = "graph.laplacian.recursive"},
#'   weights are graph trend-filtering weights or conductances. For
#'   \code{operator.family = "path.divided.difference"}, weights are metric edge
#'   lengths used to define path coordinates. If \code{NULL}, all values are one.
#' @param order Integer trend-filtering order. Supported values are
#'   \code{0L}, \code{1L}, and \code{2L}.
#' @param weight.rule Character scalar. \code{"conductance"} uses
#'   \code{weight.list} directly as \eqn{\omega_e}; \code{"sqrt.conductance"}
#'   uses \eqn{\sqrt{w_e}}; \code{"unit"} ignores \code{weight.list}. This
#'   applies to \code{operator.family = "graph.laplacian.recursive"}.
#' @param operator.family Character scalar selecting the penalty operator family.
#'   The default \code{"graph.laplacian.recursive"} preserves the existing
#'   graph trend-filtering semantics. \code{"path.divided.difference"} uses
#'   local pathwise finite differences and is designed to reproduce classical
#'   one-dimensional trend filtering on path graphs.
#' @param path.family Character scalar used for
#'   \code{operator.family = "path.divided.difference"}. \code{"branch.continuation"}
#'   keeps paths whose internal vertices do not pass through branch points.
#'   \code{"all.simple"} keeps all simple paths of the required length.
#' @param path.weighting Character scalar used for
#'   \code{operator.family = "path.divided.difference"}. \code{"unit"} assigns
#'   every canonical path row weight one. \code{"anchor.mean"} divides rows by
#'   the number of retained canonical paths sharing the same canonical start
#'   vertex.
#' @param return.sparse Logical. If \code{TRUE}, attach \pkg{Matrix} sparse
#'   incidence, Laplacian, and penalty matrices.
#'
#' @return A list of class \code{"graph.trend.filtering.operator"} containing
#'   the canonical graph, edge table, trend-filtering weights, incidence,
#'   weighted Laplacian, selected penalty operator, and nullity estimate.
#' @export
graph.trend.filtering.operator <- function(adj.list,
                                           weight.list = NULL,
                                           order = 0L,
                                           weight.rule = c("conductance",
                                                           "sqrt.conductance",
                                                           "unit"),
                                           operator.family = c("graph.laplacian.recursive",
                                                               "path.divided.difference"),
                                           path.family = c("branch.continuation",
                                                           "all.simple"),
                                           path.weighting = c("unit",
                                                              "anchor.mean"),
                                           return.sparse = TRUE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for graph.trend.filtering.operator().",
             call. = FALSE)
    }
    args <- .validate.graph.trend.filtering.operator.args(
        order = order,
        weight.rule = weight.rule,
        operator.family = operator.family,
        path.family = path.family,
        path.weighting = path.weighting
    )
    graph <- .validate.graph.trend.filtering.graph(adj.list, weight.list)
    edge.table <- .graph.trend.filtering.edge.table(graph, args$weight.rule)
    if (nrow(edge.table) < 1L) {
        stop("The graph must contain at least one undirected edge.", call. = FALSE)
    }
    if (identical(args$operator.family, "path.divided.difference")) {
        matrices <- .graph.trend.filtering.path.matrices(
            graph = graph,
            edge.table = edge.table,
            order = args$order,
            path.family = args$path.family,
            path.weighting = args$path.weighting
        )
    } else {
        matrices <- .graph.trend.filtering.matrices(
            edge.table = edge.table,
            n.vertices = graph$n.vertices,
            order = args$order
        )
    }

    out <- list(
        graph = list(
            n.vertices = graph$n.vertices,
            n.edges = nrow(edge.table),
            adj.list = graph$adj.list,
            weight.list = graph$weight.list
        ),
        order = args$order,
        weight.rule = args$weight.rule,
        operator.family = args$operator.family,
        path.family = args$path.family,
        path.weighting = args$path.weighting,
        edge.table = edge.table,
        incidence = .graph.trend.filtering.matrix.payload(matrices$incidence, return.sparse),
        laplacian = .graph.trend.filtering.matrix.payload(matrices$laplacian, return.sparse),
        penalty = .graph.trend.filtering.matrix.payload(matrices$penalty, return.sparse),
        nullity.estimate = .estimate.graph.trend.filtering.nullity(matrices$penalty),
        path.table = matrices$path.table
    )
    out$penalty$order <- args$order
    if (identical(args$operator.family, "path.divided.difference")) {
        out$penalty$kind <- "path_divided_difference"
        out$penalty$difference.order <- args$order + 1L
    } else {
        out$penalty$kind <- switch(
            as.character(args$order),
            "0" = "incidence",
            "1" = "laplacian",
            "2" = "incidence_laplacian"
        )
    }
    class(out) <- c("graph.trend.filtering.operator", "list")
    out
}

#' Fit Graph Trend Filtering
#'
#' Fits graph trend filtering on a supplied undirected graph. Supported
#' phase-2 orders are \code{0L}, \code{1L}, and \code{2L}:
#' \deqn{
#' \widehat\beta_\lambda =
#' \arg\min_{\beta\in\mathbb R^n}
#' \left\{
#' \frac{1}{2}\|y-\beta\|_2^2
#' + \lambda\|\Delta_w^{(k+1)}\beta\|_1
#' \right\}.
#' }
#'
#' @inheritParams graph.trend.filtering.operator
#' @param y Numeric response vector of length \code{length(adj.list)}.
#' @param lambda.grid Optional non-negative lambda grid. For
#'   \code{lambda.selection = "fixed"}, this must contain exactly one value. For
#'   \code{lambda.selection = "cv"}, \code{NULL} builds a default grid from the
#'   fitted full-data solution path.
#' @param lambda.selection \code{"cv"} for vertex K-fold cross-validation or
#'   \code{"fixed"} for a supplied fixed lambda.
#' @param n.lambda Number of default lambda candidates when
#'   \code{lambda.grid = NULL} and \code{lambda.selection = "cv"}.
#' @param nfolds Number of CV folds.
#' @param foldid Optional integer fold assignments of length \code{length(y)}.
#' @param maxsteps Maximum number of path steps passed to \pkg{genlasso}.
#' @param minlam Minimum lambda passed to \pkg{genlasso}.
#' @param approx Logical; passed to \pkg{genlasso}.
#' @param rtol,btol,eps Numerical controls passed to \pkg{genlasso}.
#' @param verbose Logical. If \code{TRUE}, pass verbose output to
#'   \pkg{genlasso}.
#'
#' @details
#' Graph trend filtering is intended as an \eqn{\ell_1}-adaptive comparator to
#' the \eqn{\ell_2} graph smoothers in \code{\link{fit.metric.graph.lowpass}}
#' and \code{\link{fit.rdgraph.regression}}. The low-pass smoother penalizes
#' broad quadratic graph roughness through a term like
#' \deqn{\eta f^\top L f,}
#' which shrinks high-frequency variation everywhere. Graph trend filtering
#' instead penalizes absolute graph differences. For \code{order = 0L}, this is
#' \deqn{\lambda \sum_e \omega_e |\beta_j-\beta_i|,}
#' so many edge differences can become exactly zero while selected edges carry
#' jumps. In this sense \code{order = 0L} is locally adaptive: fitted values can
#' be piecewise constant over graph regions. Higher orders use
#' \deqn{\Delta_w^{(1)} = D_w,\qquad
#'       \Delta_w^{(2)} = L_w,\qquad
#'       \Delta_w^{(3)} = D_w L_w,}
#' with \eqn{L_w = D_w^\top D_w}.
#'
#' The supplied \code{weight.list}, when used, is interpreted as a graph
#' trend-filtering edge weight or conductance for
#' \code{operator.family = "graph.laplacian.recursive"}. For
#' \code{operator.family = "path.divided.difference"}, \code{weight.list} is
#' interpreted as metric edge length and is used to form local path coordinates.
#' Use \code{weight.rule = "sqrt.conductance"} when the recursive graph
#' square-root convention is desired, and \code{weight.rule = "unit"} for an
#' unweighted recursive graph fused-lasso penalty.
#'
#' The solver backend is currently \pkg{genlasso}. Cross-validation refits the
#' generalized-lasso problem on vertex-held-out training sets using a selection
#' matrix \eqn{X}; this is useful for small and moderate diagnostic problems but
#' is not yet a scalable production graph-trend-filtering solver.
#'
#' @return A list of class \code{"graph.trend.filtering.fit"} containing fitted
#'   values, residuals, selected lambda, the operator, path metadata, and CV
#'   diagnostics when requested.
#' @references
#' Wang, Y.-X., Sharpnack, J., Smola, A., and Tibshirani, R. J. (2016).
#' Trend Filtering on Graphs. \emph{Journal of Machine Learning Research}.
#'
#' @export
fit.graph.trend.filtering <- function(adj.list,
                                      weight.list = NULL,
                                      y,
                                      order = 0L,
                                      lambda.grid = NULL,
                                      lambda.selection = c("cv", "fixed"),
                                      weight.rule = c("conductance",
                                                      "sqrt.conductance",
                                                      "unit"),
                                      operator.family = c("graph.laplacian.recursive",
                                                          "path.divided.difference"),
                                      path.family = c("branch.continuation",
                                                      "all.simple"),
                                      path.weighting = c("unit",
                                                         "anchor.mean"),
                                      n.lambda = 40L,
                                      nfolds = 5L,
                                      foldid = NULL,
                                      maxsteps = 2000L,
                                      minlam = 0,
                                      approx = FALSE,
                                      rtol = 1e-7,
                                      btol = 1e-7,
                                      eps = 1e-4,
                                      verbose = FALSE) {
    if (!requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for fit.graph.trend.filtering().",
             call. = FALSE)
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for fit.graph.trend.filtering().",
             call. = FALSE)
    }

    lambda.selection <- match.arg(lambda.selection)
    operator <- graph.trend.filtering.operator(
        adj.list = adj.list,
        weight.list = weight.list,
        order = order,
        weight.rule = weight.rule,
        operator.family = operator.family,
        path.family = path.family,
        path.weighting = path.weighting,
        return.sparse = TRUE
    )
    n <- operator$graph$n.vertices
    y <- .validate.metric.graph.lowpass.response(y, n, "y")
    solver.args <- .validate.graph.trend.filtering.solver.args(
        n.lambda = n.lambda,
        nfolds = nfolds,
        foldid = foldid,
        n = n,
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        verbose = verbose
    )
    lambda.grid <- .validate.graph.trend.filtering.lambda.grid(
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection
    )
    if (lambda.selection == "fixed" && is.null(lambda.grid)) {
        stop("lambda.grid must be supplied when lambda.selection = 'fixed'.",
             call. = FALSE)
    }

    penalty.solver <- .graph.trend.filtering.solver.penalty(operator)

    path <- genlasso::genlasso(
        y = y,
        D = penalty.solver$D,
        approx = solver.args$approx,
        maxsteps = solver.args$maxsteps,
        minlam = solver.args$minlam,
        rtol = solver.args$rtol,
        btol = solver.args$btol,
        eps = solver.args$eps,
        verbose = solver.args$verbose,
        svd = penalty.solver$svd
    )

    if (is.null(lambda.grid)) {
        lambda.grid <- .graph.trend.filtering.default.lambda.grid(path, solver.args$n.lambda)
    }
    if (lambda.selection == "fixed" && length(lambda.grid) != 1L) {
        stop("lambda.grid must contain exactly one value when lambda.selection = 'fixed'.",
             call. = FALSE)
    }

    coef.full <- tryCatch(
        stats::coef(path, lambda = lambda.grid),
        error = function(e) NULL
    )
    beta.grid <- if (is.null(coef.full)) NULL else as.matrix(coef.full$beta)
    if (is.null(beta.grid) || ncol(beta.grid) != length(lambda.grid)) {
        beta.grid <- vapply(lambda.grid, function(lambda) {
            tryCatch(
                as.vector(stats::coef(path, lambda = lambda)$beta),
                error = function(e) rep(NA_real_, length(y))
            )
        }, numeric(length(y)))
    }
    colnames(beta.grid) <- format(signif(lambda.grid, 6), scientific = TRUE)
    finite.beta <- colSums(is.finite(beta.grid)) == nrow(beta.grid)

    cv <- NULL
    if (lambda.selection == "cv") {
        foldid <- solver.args$foldid
        cv <- .graph.trend.filtering.cv(
            y = y,
            D = penalty.solver$D,
            lambda.grid = lambda.grid,
            foldid = foldid,
            solver.args = solver.args,
            use.svd = penalty.solver$svd
        )
        cv$mean.error[!finite.beta] <- Inf
        best.idx <- cv$best.idx
        if (!is.finite(cv$mean.error[best.idx])) {
            best.idx <- which.min(cv$mean.error)
        }
        if (!is.finite(cv$mean.error[best.idx])) {
            stop("No graph trend-filtering lambda produced finite CV and fitted values.",
                 call. = FALSE)
        }
        cv$best.idx <- best.idx
        cv$lambda.min <- lambda.grid[best.idx]
        cv$error.min <- cv$mean.error[best.idx]
    } else {
        best.idx <- 1L
        if (!finite.beta[best.idx]) {
            stop("The requested graph trend-filtering lambda produced non-finite fitted values.",
                 call. = FALSE)
        }
    }

    fitted.values <- as.vector(beta.grid[, best.idx])
    if (any(!is.finite(fitted.values))) {
        stop("Graph trend filtering produced non-finite fitted values.",
             call. = FALSE)
    }
    residuals <- y - fitted.values
    result <- list(
        fitted.values = fitted.values,
        residuals = residuals,
        y = y,
        lambda = lambda.grid[best.idx],
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        order = operator$order,
        weight.rule = operator$weight.rule,
        operator.family = operator$operator.family,
        path.family = operator$path.family,
        path.weighting = operator$path.weighting,
        operator = operator,
        graph = operator$graph,
        edge.table = operator$edge.table,
        path.table = operator$path.table,
        penalty = operator$penalty,
        nullity.estimate = operator$nullity.estimate,
        beta.grid = beta.grid,
        path = list(
            lambda = path$lambda,
            df = path$df,
            completepath = path$completepath,
            backend = "genlasso",
            solver.penalty = penalty.solver$representation,
            svd = penalty.solver$svd
        ),
        cv = cv,
        parameters = c(
            list(lambda.selection = lambda.selection),
            solver.args[setdiff(names(solver.args), "foldid")]
        )
    )
    attr(result, "call") <- match.call()
    class(result) <- c("graph.trend.filtering.fit", "list")
    result
}

.validate.graph.trend.filtering.graph <- function(adj.list, weight.list) {
    if (is.null(weight.list)) {
        if (!is.list(adj.list)) stop("adj.list must be a list.")
        weight.list <- lapply(adj.list, function(nbrs) rep(1, length(nbrs)))
    }
    graph <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    graph$n.vertices <- length(graph$adj.list)
    graph
}

.validate.graph.trend.filtering.operator.args <- function(order,
                                                         weight.rule = c("conductance",
                                                                         "sqrt.conductance",
                                                                         "unit"),
                                                         operator.family = c("graph.laplacian.recursive",
                                                                             "path.divided.difference"),
                                                         path.family = c("branch.continuation",
                                                                         "all.simple"),
                                                         path.weighting = c("unit",
                                                                            "anchor.mean")) {
    if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
        !is.finite(order) || order != floor(order) || order < 0) {
        stop("order must be a non-negative integer scalar.", call. = FALSE)
    }
    order <- as.integer(order)
    if (!(order %in% 0:2)) {
        stop("Phase-2 graph trend filtering supports only order = 0L, 1L, or 2L.",
             call. = FALSE)
    }
    weight.rule <- match.arg(weight.rule)
    operator.family <- match.arg(operator.family)
    path.family <- match.arg(path.family)
    path.weighting <- match.arg(path.weighting)
    list(
        order = order,
        weight.rule = weight.rule,
        operator.family = operator.family,
        path.family = path.family,
        path.weighting = path.weighting
    )
}

.validate.graph.trend.filtering.solver.args <- function(n.lambda, nfolds, foldid,
                                                       n, maxsteps, minlam,
                                                       approx, rtol, btol,
                                                       eps, verbose) {
    n.lambda <- .validate.positive.integer.scalar(n.lambda, "n.lambda")
    nfolds <- .validate.positive.integer.scalar(nfolds, "nfolds")
    maxsteps <- .validate.positive.integer.scalar(maxsteps, "maxsteps")
    if (nfolds < 2L) stop("nfolds must be at least 2.", call. = FALSE)
    if (nfolds > n) nfolds <- n
    if (!is.null(foldid)) {
        if (!is.numeric(foldid) && !is.integer(foldid)) stop("foldid must be an integer vector.", call. = FALSE)
        if (length(foldid) != n || anyNA(foldid)) stop("foldid must have length length(y) with no NA.", call. = FALSE)
        foldid <- as.integer(foldid)
        if (any(foldid < 1L)) stop("foldid values must be positive integers.", call. = FALSE)
        fold.levels <- sort(unique(foldid))
        foldid <- match(foldid, fold.levels)
        nfolds <- length(fold.levels)
        if (nfolds < 2L) stop("foldid must contain at least two folds.", call. = FALSE)
    } else {
        foldid <- sample(rep(seq_len(nfolds), length.out = n))
    }
    if (!is.numeric(minlam) || length(minlam) != 1L ||
        !is.finite(minlam) || minlam < 0) {
        stop("minlam must be a finite non-negative numeric scalar.", call. = FALSE)
    }
    for (nm in c("rtol", "btol", "eps")) {
        val <- get(nm)
        if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
            stop(sprintf("%s must be a finite positive numeric scalar.", nm), call. = FALSE)
        }
    }
    for (nm in c("approx", "verbose")) {
        val <- get(nm)
        if (!is.logical(val) || length(val) != 1L || is.na(val)) {
            stop(sprintf("%s must be TRUE or FALSE.", nm), call. = FALSE)
        }
    }
    list(
        n.lambda = n.lambda,
        nfolds = nfolds,
        foldid = foldid,
        maxsteps = maxsteps,
        minlam = as.double(minlam),
        approx = as.logical(approx),
        rtol = as.double(rtol),
        btol = as.double(btol),
        eps = as.double(eps),
        verbose = as.logical(verbose)
    )
}

.validate.graph.trend.filtering.lambda.grid <- function(lambda.grid, lambda.selection) {
    if (is.null(lambda.grid)) return(NULL)
    if (!is.numeric(lambda.grid) || length(lambda.grid) < 1L ||
        any(!is.finite(lambda.grid)) || any(lambda.grid < 0)) {
        stop("lambda.grid must be NULL or a non-negative finite numeric vector.",
             call. = FALSE)
    }
    lambda.grid <- unique(as.double(lambda.grid))
    sort(lambda.grid, decreasing = TRUE)
}

.graph.trend.filtering.edge.table <- function(graph, weight.rule) {
    from <- integer()
    to <- integer()
    base.weight <- numeric()
    for (i in seq_len(graph$n.vertices)) {
        nbrs <- graph$adj.list[[i]]
        wts <- graph$weight.list[[i]]
        keep <- nbrs > i
        if (any(keep)) {
            from <- c(from, rep.int(i, sum(keep)))
            to <- c(to, nbrs[keep])
            base.weight <- c(base.weight, wts[keep])
        }
    }
    trend.weight <- switch(
        weight.rule,
        conductance = base.weight,
        sqrt.conductance = sqrt(base.weight),
        unit = rep(1, length(base.weight))
    )
    data.frame(
        edge = seq_along(from),
        from = from,
        to = to,
        base.weight = base.weight,
        trend.weight = trend.weight
    )
}

.graph.trend.filtering.incidence.triplets <- function(edge.table, n.vertices) {
    n.edges <- nrow(edge.table)
    i <- rep(seq_len(n.edges), each = 2L)
    j <- as.integer(rbind(edge.table$from, edge.table$to))
    x <- as.double(rbind(-edge.table$trend.weight, edge.table$trend.weight))
    list(
        i = i,
        j = j,
        x = x,
        dim = c(n.edges, n.vertices)
    )
}

.graph.trend.filtering.matrices <- function(edge.table, n.vertices, order) {
    inc <- .graph.trend.filtering.incidence.triplets(edge.table, n.vertices)
    D <- Matrix::sparseMatrix(
        i = inc$i,
        j = inc$j,
        x = inc$x,
        dims = inc$dim,
        giveCsparse = TRUE
    )
    L <- Matrix::t(D) %*% D
    penalty <- switch(
        as.character(order),
        "0" = D,
        "1" = L,
        "2" = D %*% L
    )
    list(incidence = D, laplacian = L, penalty = penalty)
}

.graph.trend.filtering.path.matrices <- function(graph, edge.table, order,
                                                 path.family, path.weighting) {
    base <- .graph.trend.filtering.matrices(
        edge.table = edge.table,
        n.vertices = graph$n.vertices,
        order = 0L
    )
    difference.order <- order + 1L
    path.info <- .graph.trend.filtering.enumerate.paths(
        graph = graph,
        difference.order = difference.order,
        path.family = path.family
    )
    if (length(path.info$paths) < 1L) {
        stop(sprintf(
            "No length-%d paths are available for path.divided.difference.",
            difference.order
        ), call. = FALSE)
    }
    path.weights <- .graph.trend.filtering.path.row.weights(
        paths = path.info$paths,
        path.weighting = path.weighting
    )
    mat <- .graph.trend.filtering.path.penalty.matrix(
        paths = path.info$paths,
        graph = graph,
        row.weights = path.weights
    )
    path.table <- .graph.trend.filtering.path.table(
        paths = path.info$paths,
        graph = graph,
        row.weights = path.weights,
        path.family = path.family,
        path.weighting = path.weighting,
        duplicates.removed = path.info$duplicates.removed
    )
    list(
        incidence = base$incidence,
        laplacian = base$laplacian,
        penalty = mat,
        path.table = path.table
    )
}

.graph.trend.filtering.enumerate.paths <- function(graph, difference.order,
                                                   path.family) {
    if (difference.order < 1L) stop("difference.order must be positive.", call. = FALSE)
    degrees <- lengths(graph$adj.list)
    directed <- list()
    for (start in seq_len(graph$n.vertices)) {
        directed <- c(
            directed,
            .graph.trend.filtering.extend.paths(
                graph = graph,
                path = start,
                remaining = difference.order
            )
        )
    }
    keep <- vapply(directed, function(path) {
        if (!identical(path.family, "branch.continuation")) return(TRUE)
        if (length(path) <= 2L) return(TRUE)
        internal <- path[seq.int(2L, length(path) - 1L)]
        all(degrees[internal] <= 2L)
    }, logical(1))
    directed <- directed[keep]

    canonical <- list()
    seen <- character()
    duplicates.removed <- 0L
    for (path in directed) {
        can <- .graph.trend.filtering.canonical.path(path)
        key <- paste(can, collapse = ":")
        if (key %in% seen) {
            duplicates.removed <- duplicates.removed + 1L
        } else {
            seen <- c(seen, key)
            canonical[[length(canonical) + 1L]] <- can
        }
    }
    canonical <- canonical[order(vapply(canonical, paste, character(1), collapse = ":"))]
    list(paths = canonical, duplicates.removed = duplicates.removed)
}

.graph.trend.filtering.extend.paths <- function(graph, path, remaining) {
    if (remaining == 0L) return(list(as.integer(path)))
    current <- path[length(path)]
    nbrs <- graph$adj.list[[current]]
    nbrs <- nbrs[!(nbrs %in% path)]
    if (length(nbrs) < 1L) return(list())
    out <- list()
    for (nbr in nbrs) {
        out <- c(
            out,
            .graph.trend.filtering.extend.paths(
                graph = graph,
                path = c(path, nbr),
                remaining = remaining - 1L
            )
        )
    }
    out
}

.graph.trend.filtering.canonical.path <- function(path) {
    rev.path <- rev(path)
    path.key <- paste(sprintf("%010d", path), collapse = ":")
    rev.key <- paste(sprintf("%010d", rev.path), collapse = ":")
    if (rev.key < path.key) rev.path else path
}

.graph.trend.filtering.path.row.weights <- function(paths, path.weighting) {
    if (identical(path.weighting, "unit")) {
        return(rep(1, length(paths)))
    }
    anchors <- vapply(paths, `[`, integer(1), 1L)
    counts <- table(anchors)
    as.double(1 / counts[as.character(anchors)])
}

.graph.trend.filtering.path.penalty.matrix <- function(paths, graph, row.weights) {
    ii <- integer()
    jj <- integer()
    xx <- numeric()
    for (row in seq_along(paths)) {
        path <- paths[[row]]
        positions <- .graph.trend.filtering.path.positions(path, graph)
        coefs <- .graph.trend.filtering.path.coefficients(positions) * row.weights[row]
        keep <- abs(coefs) > 0
        ii <- c(ii, rep.int(row, sum(keep)))
        jj <- c(jj, path[keep])
        xx <- c(xx, coefs[keep])
    }
    Matrix::sparseMatrix(
        i = ii,
        j = jj,
        x = xx,
        dims = c(length(paths), graph$n.vertices),
        giveCsparse = TRUE
    )
}

.graph.trend.filtering.path.positions <- function(path, graph) {
    if (length(path) < 2L) return(0)
    lengths <- numeric(length(path) - 1L)
    for (idx in seq_along(lengths)) {
        from <- path[idx]
        to <- path[idx + 1L]
        nbrs <- graph$adj.list[[from]]
        pos <- match(to, nbrs)
        if (is.na(pos)) stop("Internal error: path contains a missing graph edge.", call. = FALSE)
        lengths[idx] <- graph$weight.list[[from]][pos]
    }
    c(0, cumsum(lengths))
}

.graph.trend.filtering.path.coefficients <- function(positions) {
    n <- length(positions)
    if (n == 2L) return(c(-1, 1))
    rows <- matrix(0, nrow = n - 1L, ncol = n)
    for (idx in seq_len(n - 1L)) {
        rows[idx, idx] <- -1
        rows[idx, idx + 1L] <- 1
    }
    for (level in seq_len(n - 2L)) {
        n.old <- nrow(rows)
        weights <- level / (positions[(level + 1L):n] - positions[seq_len(n - level)])
        new.rows <- matrix(0, nrow = n.old - 1L, ncol = n)
        for (idx in seq_len(n.old - 1L)) {
            new.rows[idx, ] <- -weights[idx] * rows[idx, ] +
                weights[idx + 1L] * rows[idx + 1L, ]
        }
        rows <- new.rows
    }
    as.double(rows[1L, ])
}

.graph.trend.filtering.path.table <- function(paths, graph, row.weights,
                                              path.family, path.weighting,
                                              duplicates.removed) {
    vertices <- vapply(paths, paste, character(1), collapse = "-")
    path.length <- vapply(paths, function(path) {
        positions <- .graph.trend.filtering.path.positions(path, graph)
        positions[length(positions)]
    }, numeric(1))
    data.frame(
        row = seq_along(paths),
        start = vapply(paths, `[`, integer(1), 1L),
        end = vapply(paths, function(path) path[length(path)], integer(1)),
        vertices = vertices,
        metric.length = path.length,
        row.weight = row.weights,
        path.family = path.family,
        path.weighting = path.weighting,
        duplicates.removed = c(duplicates.removed, rep(NA_integer_, length(paths) - 1L))
    )
}

.graph.trend.filtering.matrix.payload <- function(mat, return.sparse) {
    trip <- methods::as(mat, "TsparseMatrix")
    payload <- list(
        i = as.integer(trip@i + 1L),
        j = as.integer(trip@j + 1L),
        x = as.double(trip@x),
        dim = dim(mat)
    )
    if (isTRUE(return.sparse)) payload$matrix <- mat
    payload
}

.estimate.graph.trend.filtering.nullity <- function(mat) {
    if (ncol(mat) > 200L) return(NA_integer_)
    if (nrow(mat) == 0L) return(as.integer(ncol(mat)))
    rank <- Matrix::rankMatrix(as.matrix(mat))[1]
    as.integer(ncol(mat) - rank)
}

.graph.trend.filtering.solver.penalty <- function(operator) {
    if (identical(operator$operator.family, "graph.laplacian.recursive") &&
        (identical(operator$order, 1L) ||
         nrow(operator$penalty$matrix) > ncol(operator$penalty$matrix))) {
        return(list(
            D = as.matrix(operator$penalty$matrix),
            svd = TRUE,
            representation = "dense"
        ))
    }
    list(
        D = operator$penalty$matrix,
        svd = FALSE,
        representation = "sparse"
    )
}

.graph.trend.filtering.default.lambda.grid <- function(path, n.lambda) {
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
        lambda.max * 1e-4
    } else {
        max(min(lambda.positive) * (1 + 1e-8), lambda.max * 1e-4)
    }
    grid <- exp(seq(log(lambda.max), log(lambda.min), length.out = n.lambda - 1L))
    if (complete.path) grid <- c(grid, 0)
    unique(grid)
}

.graph.trend.filtering.cv <- function(y, D, lambda.grid, foldid, solver.args,
                                      use.svd = FALSE) {
    n <- length(y)
    folds <- sort(unique(foldid))
    cv.errors <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda.grid))
    rownames(cv.errors) <- paste0("Fold", folds)
    colnames(cv.errors) <- format(signif(lambda.grid, 6), scientific = TRUE)

    for (ii in seq_along(folds)) {
        fold <- folds[ii]
        test <- which(foldid == fold)
        train <- which(foldid != fold)
        X <- as.matrix(Matrix::sparseMatrix(
            i = seq_along(train),
            j = train,
            x = 1,
            dims = c(length(train), n),
            giveCsparse = TRUE
        ))
        path <- suppressWarnings(genlasso::genlasso(
            y = y[train],
            X = X,
            D = D,
            approx = solver.args$approx,
            maxsteps = solver.args$maxsteps,
            minlam = solver.args$minlam,
            rtol = solver.args$rtol,
            btol = solver.args$btol,
            eps = solver.args$eps,
            verbose = solver.args$verbose,
            svd = use.svd
        ))
        beta <- tryCatch(
            as.matrix(stats::coef(path, lambda = lambda.grid)$beta),
            error = function(e) NULL
        )
        if (is.null(beta) || ncol(beta) != length(lambda.grid)) {
            beta <- vapply(lambda.grid, function(lambda) {
                tryCatch(
                    as.vector(stats::coef(path, lambda = lambda)$beta),
                    error = function(e) rep(NA_real_, n)
                )
            }, numeric(n))
        }
        beta.test <- beta[test, , drop = FALSE]
        finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
        fold.error <- colMeans((matrix(y[test], nrow = length(test), ncol = length(lambda.grid)) -
                                beta.test)^2)
        fold.error[!finite.pred] <- NA_real_
        cv.errors[ii, ] <- fold.error
    }

    mean.error <- colMeans(cv.errors)
    mean.error[!is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, 2, function(x) {
        if (any(!is.finite(x)) || length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    best.idx <- which.min(mean.error)
    if (!is.finite(mean.error[best.idx])) {
        stop("All graph trend-filtering CV fits failed for the supplied lambda grid.",
             call. = FALSE)
    }
    list(
        foldid = foldid,
        lambda.grid = lambda.grid,
        fold.errors = cv.errors,
        mean.error = mean.error,
        se = se,
        best.idx = best.idx,
        lambda.min = lambda.grid[best.idx],
        error.min = mean.error[best.idx]
    )
}

#' @export
print.graph.trend.filtering.fit <- function(x, ...) {
    cat("\nGraph Trend-Filtering Fit\n")
    cat("=========================\n\n")
    cat(sprintf("Vertices: %d; edges: %d\n", x$graph$n.vertices, x$graph$n.edges))
    cat(sprintf("Order: %d; operator family: %s\n", x$order, x$operator.family))
    cat(sprintf("Weight rule: %s\n", x$weight.rule))
    if (identical(x$operator.family, "path.divided.difference")) {
        cat(sprintf("Path family: %s; path weighting: %s\n",
                    x$path.family, x$path.weighting))
    }
    cat(sprintf("Lambda selection: %s; lambda: %.4e\n", x$lambda.selection, x$lambda))
    if (!is.null(x$cv)) {
        cat(sprintf("CV error: %.4e\n", x$cv$error.min))
    }
    invisible(x)
}
