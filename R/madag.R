#' Construct Monotonic Ascent DAG (MADAG) from a Source Vertex
#'
#' @description
#' Constructs a Monotonic Ascent DAG rooted at a specified local minimum vertex.
#' The MADAG captures all vertices reachable via monotonically ascending paths,
#' enabling identification of reachable maxima (cells) and enumeration of
#' trajectories within each cell.
#'
#' @details
#' The MADAG is a directed acyclic graph where edges connect vertices u to v
#' when (u,v) is an edge in the original graph and \code{y[v] > y[u]}. The DAG is
#' rooted at the specified source vertex (typically a local minimum) and
#' contains all vertices reachable via monotonically ascending paths.
#'
#' For each reachable vertex, the MADAG tracks all monotonic predecessors
#' (enabling trajectory enumeration) and successors. Vertices with no
#' successors in the MADAG are local maxima (sinks), defining the cells
#' of the Morse-Smale decomposition.
#'
#' Unlike smooth manifold settings where trajectory uniqueness is guaranteed
#' by ODE theory, discrete graphs can have multiple trajectories passing
#' through the same vertex, leading to multi-valued cell memberships.
#'
#' @param adj.list List of integer vectors. Each element \code{adj.list[[i]]}
#'   contains the 1-based indices of vertices adjacent to vertex \code{i}.
#' @param weight.list List of numeric vectors. Each element
#'   \code{weight.list[[i]]} contains the edge weights (distances) corresponding
#'   to \code{adj.list[[i]]}.
#' @param y Numeric vector of function values at each vertex.
#' @param source.vertex Integer. The 1-based index of the source vertex
#'   (typically a local minimum) from which to construct the MADAG.
#' @param max.trajectories.per.cell Integer. Maximum number of trajectories
#'   to enumerate per cell. If exceeded, trajectories are not enumerated
#'   (use sampling instead). Default is 10000.
#' @param compute.path.counts Logical. Whether to compute path counts for
#'   trajectory counting without full enumeration. Default is \code{TRUE}.
#' @param enumerate.trajectories Logical. Whether to enumerate trajectories
#'   during construction (up to \code{max.trajectories.per.cell}). Default
#'   is \code{TRUE}.
#' @param edge.length.quantile.thld Numeric in \eqn{(0,1]}. Quantile threshold
#'   for edge length filtering. Edges longer than this quantile are excluded.
#'   Default is 1.0 (no filtering).
#' @param min.cell.support Integer. Minimum cell support size to retain.
#'   Default is 1.
#' @param verbose Logical. Whether to print progress messages. Default is
#'   \code{FALSE}.
#'
#' @return A list of class \code{"madag"} containing:
#'   \item{source.vertex}{The 1-based index of the source vertex.}
#'   \item{source.value}{The y value at the source vertex.}
#'   \item{reachable.vertices}{Integer vector of 1-based indices of all
#'     vertices reachable from source via ascending paths.}
#'   \item{reachable.maxima}{Integer vector of 1-based indices of local
#'     maxima reachable from source (sinks of the MADAG).}
#'   \item{predecessors}{Named list where \code{predecessors[[as.character(v)]]}
#'     contains the 1-based indices of monotonic predecessors of vertex v.}
#'   \item{successors}{Named list where \code{successors[[as.character(v)]]}
#'     contains the 1-based indices of monotonic successors of vertex v.}
#'   \item{topological.order}{Integer vector of vertices in topological order
#'     (source first, sinks last).}
#'   \item{path.count.from.source}{Named integer vector where
#'     \code{path.count.from.source[[as.character(v)]]} is the number of
#'     distinct paths from source to v.}
#'   \item{path.count.to.sinks}{Named integer vector where
#'     \code{path.count.to.sinks[[as.character(v)]]} is the number of
#'     distinct paths from v to any sink.}
#'   \item{cells}{List of cell structures, one for each reachable maximum.
#'     Each cell contains: \code{min.vertex}, \code{max.vertex},
#'     \code{min.value}, \code{max.value}, \code{support},
#'     \code{n.trajectories}, \code{explicitly.enumerated},
#'     \code{bottlenecks}, \code{n.clusters}.}
#'   \item{n.vertices}{Number of vertices in the MADAG.}
#'   \item{n.cells}{Number of cells (reachable maxima).}
#'
#' @seealso \code{\link{enumerate.cell.trajectories}},
#'   \code{\link{sample.cell.trajectories}},
#'   \code{\link{identify.bottlenecks}}
#'
#' @examples
#' \dontrun{
#' ## Create a simple graph
#' adj.list <- list(
#'   c(2L, 3L),
#'   c(1L, 3L, 4L),
#'   c(1L, 2L, 4L),
#'   c(2L, 3L)
#' )
#' weight.list <- list(
#'   c(1.0, 1.5),
#'   c(1.0, 1.2, 1.0),
#'   c(1.5, 1.2, 1.0),
#'   c(1.0, 1.0)
#' )
#' y <- c(0.5, 1.2, 0.8, 2.0)
#'
#' ## Construct MADAG from vertex 1 (local minimum)
#' madag <- construct.madag(adj.list, weight.list, y, source.vertex = 1L)
#'
#' ## Print summary
#' print(madag)
#' }
#'
#' @export
construct.madag <- function(
    adj.list,
    weight.list,
    y,
    source.vertex,
    max.trajectories.per.cell = 10000L,
    compute.path.counts = TRUE,
    enumerate.trajectories = TRUE,
    edge.length.quantile.thld = 1.0,
    min.cell.support = 1L,
    verbose = FALSE
) {
    ## -------------------------------------------------------------------------
    ## Input validation
    ## -------------------------------------------------------------------------

    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)
    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }
    if (length(y) != n.vertices) {
        stop("y must have the same length as adj.list")
    }

    if (!is.numeric(y)) {
        stop("y must be a numeric vector")
    }

    if (!is.numeric(source.vertex) || length(source.vertex) != 1) {
        stop("source.vertex must be a single integer")
    }
    source.vertex <- as.integer(source.vertex)
    if (source.vertex < 1L || source.vertex > n.vertices) {
        stop("source.vertex must be between 1 and ", n.vertices)
    }

    ## -------------------------------------------------------------------------
    ## Convert adjacency list to 0-based for C++
    ## -------------------------------------------------------------------------

    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            return(integer(0))
        }
        as.integer(x - 1L)
    })

    ## -------------------------------------------------------------------------
    ## Build parameter list
    ## -------------------------------------------------------------------------

    params <- list(
        max.trajectories.per.cell = as.integer(max.trajectories.per.cell),
        compute.path.counts = as.logical(compute.path.counts),
        enumerate.trajectories = as.logical(enumerate.trajectories),
        edge.length.quantile.thld = as.double(edge.length.quantile.thld),
        min.cell.support = as.integer(min.cell.support)
    )

    ## -------------------------------------------------------------------------
    ## Call C++ implementation
    ## -------------------------------------------------------------------------

    result <- .Call(
        S_construct_madag,
        adj.list.0based,
        weight.list,
        as.double(y),
        as.integer(source.vertex),
        params,
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## -------------------------------------------------------------------------
    ## Post-process result
    ## -------------------------------------------------------------------------

    ## Convert names to dot.snake convention
    names(result) <- gsub("_", ".", names(result))

    ## Convert cell names
    if (length(result$cells) > 0) {
        result$cells <- lapply(result$cells, function(cell) {
            names(cell) <- gsub("_", ".", names(cell))
            cell
        })
    }

    ## Add class
    class(result) <- c("madag", "list")

    return(result)
}


#' Enumerate Trajectories in a Cell
#'
#' @description
#' Enumerates all monotonically ascending trajectories from the MADAG source
#' vertex to the specified maximum vertex.
#'
#' @param madag A \code{madag} object from \code{\link{construct.madag}}.
#' @param y Numeric vector of function values at each vertex.
#' @param max.vertex Integer. The 1-based index of the maximum vertex defining
#'   the cell.
#' @param max.trajectories Integer. Maximum number of trajectories to enumerate.
#'   Use 0 for unlimited. Default is 0.
#'
#' @return A list of trajectory structures, each containing:
#'   \item{vertices}{Integer vector of 1-based vertex indices along the
#'     trajectory.}
#'   \item{source.min}{The source minimum vertex (1-based).}
#'   \item{sink.max}{The sink maximum vertex (1-based).}
#'   \item{total.ascent}{Total y change along the trajectory.}
#'   \item{cluster.id}{Cluster assignment (-1 if unassigned).}
#'
#' @seealso \code{\link{construct.madag}}, \code{\link{sample.cell.trajectories}}
#'
#' @export
enumerate.cell.trajectories <- function(
    madag,
    y,
    max.vertex,
    max.trajectories = 0L
) {
    if (!inherits(madag, "madag")) {
        stop("madag must be a 'madag' object")
    }

    max.vertex <- as.integer(max.vertex)
    if (max.vertex < 1L || max.vertex > length(y)) {
        stop("max.vertex out of range")
    }

    result <- .Call(
        S_enumerate_cell_trajectories,
        madag,
        as.double(y),
        max.vertex,
        as.integer(max.trajectories),
        PACKAGE = "gflow"
    )

    ## Convert names
    result <- lapply(result, function(traj) {
        names(traj) <- gsub("_", ".", names(traj))
        class(traj) <- c("madag.trajectory", "list")
        traj
    })

    class(result) <- c("madag.trajectory.list", "list")
    return(result)
}


#' Sample Trajectories from a Cell
#'
#' @description
#' Samples trajectories from a cell by random walks from source to sink.
#' Useful when the number of trajectories is too large for full enumeration.
#'
#' @param madag A \code{madag} object from \code{\link{construct.madag}}.
#' @param y Numeric vector of function values at each vertex.
#' @param max.vertex Integer. The 1-based index of the maximum vertex defining
#'   the cell.
#' @param n.samples Integer. Number of trajectories to sample.
#' @param seed Integer. Random seed (0 for random). Default is 0.
#'
#' @return A list of trajectory structures (same format as
#'   \code{\link{enumerate.cell.trajectories}}).
#'
#' @seealso \code{\link{construct.madag}}, \code{\link{enumerate.cell.trajectories}}
#'
#' @export
sample.cell.trajectories <- function(
    madag,
    y,
    max.vertex,
    n.samples,
    seed = 0L
) {
    if (!inherits(madag, "madag")) {
        stop("madag must be a 'madag' object")
    }

    max.vertex <- as.integer(max.vertex)
    n.samples <- as.integer(n.samples)
    seed <- as.integer(seed)

    result <- .Call(
        S_sample_cell_trajectories,
        madag,
        as.double(y),
        max.vertex,
        n.samples,
        seed,
        PACKAGE = "gflow"
    )

    ## Convert names
    result <- lapply(result, function(traj) {
        names(traj) <- gsub("_", ".", names(traj))
        class(traj) <- c("madag.trajectory", "list")
        traj
    })

    class(result) <- c("madag.trajectory.list", "list")
    return(result)
}


#' Compute Trajectory Similarity Matrix
#'
#' @description
#' Computes pairwise similarity between trajectories based on vertex overlap.
#'
#' @param trajectories A list of trajectories from
#'   \code{\link{enumerate.cell.trajectories}} or
#'   \code{\link{sample.cell.trajectories}}.
#' @param similarity.type Character. Type of similarity measure.
#'   Currently only "jaccard" is supported. Default is "jaccard".
#'
#' @return A symmetric numeric matrix of pairwise similarities in \eqn{[0, 1]}.
#'
#' @seealso \code{\link{enumerate.cell.trajectories}}
#'
#' @export
trajectory.similarity.matrix <- function(
    trajectories,
    similarity.type = "jaccard"
) {
    if (!is.list(trajectories)) {
        stop("trajectories must be a list")
    }

    similarity.type <- match.arg(similarity.type, choices = c("jaccard"))

    result <- .Call(
        S_trajectory_similarity_matrix,
        trajectories,
        as.character(similarity.type),
        PACKAGE = "gflow"
    )

    return(result)
}


#' Identify Bottleneck Vertices in a Cell
#'
#' @description
#' Identifies vertices through which a large fraction of trajectories must pass.
#' These represent obligate intermediate states in the transition from minimum
#' to maximum.
#'
#' @param x A \code{madag} object from \code{\link{construct.madag}}.
#' @param max.vertex Integer. The 1-based index of the maximum vertex defining
#'   the cell.
#' @param min.fraction Numeric. Minimum fraction of trajectories passing through
#'   a vertex for it to be considered a bottleneck. Default is 0.5.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Integer vector of 1-based indices of bottleneck vertices.
#'
#' @seealso \code{\link{construct.madag}}
#'
#' @export
identify.bottlenecks <- function(
    x,
    max.vertex,
    min.fraction = 0.5,
    ...
) {
    madag <- x
    if (!inherits(madag, "madag")) {
        stop("madag must be a 'madag' object")
    }

    max.vertex <- as.integer(max.vertex)
    min.fraction <- as.double(min.fraction)

    result <- .Call(
        S_identify_bottlenecks,
        madag,
        max.vertex,
        min.fraction,
        PACKAGE = "gflow"
    )

    return(result)
}


#' Print Method for madag Objects
#'
#' @param x A \code{madag} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.madag <- function(x, ...) {
    cat("Monotonic Ascent DAG (MADAG)\n")
    cat("============================\n")
    cat(sprintf("Source vertex: %d (y = %.4f)\n",
                x$source.vertex, x$source.value))
    cat(sprintf("Reachable vertices: %d\n", x$n.vertices))
    cat(sprintf("Reachable maxima (cells): %d\n", x$n.cells))

    if (length(x$cells) > 0) {
        cat("\nCells:\n")
        for (i in seq_along(x$cells)) {
            cell <- x$cells[[i]]
            cat(sprintf("  (m%d, M%d): support=%d, trajectories=%d%s\n",
                        cell$min.vertex,
                        cell$max.vertex,
                        length(cell$support),
                        cell$n.trajectories,
                        if (cell$explicitly.enumerated) " (enumerated)" else ""))
        }
    }

    invisible(x)
}


#' Summary Method for madag Objects
#'
#' @param object A \code{madag} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.madag <- function(object, ...) {
    cat("MADAG Summary\n")
    cat("=============\n\n")

    cat("Source:\n")
    cat(sprintf("  Vertex: %d\n", object$source.vertex))
    cat(sprintf("  Value: %.4f\n", object$source.value))

    cat("\nReachability:\n")
    cat(sprintf("  Vertices: %d\n", object$n.vertices))
    cat(sprintf("  Maxima: %d\n", object$n.cells))

    if (length(object$cells) > 0) {
        cat("\nCell Statistics:\n")

        ## Create summary data frame
        cell.df <- data.frame(
            cell = sapply(object$cells, function(c)
                sprintf("(m%d,M%d)", c$min.vertex, c$max.vertex)),
            support = sapply(object$cells, function(c) length(c$support)),
            trajectories = sapply(object$cells, function(c) c$n.trajectories),
            enumerated = sapply(object$cells, function(c) c$explicitly.enumerated),
            y.change = sapply(object$cells, function(c) c$max.value - c$min.value)
        )

        print(cell.df, row.names = FALSE)

        cat(sprintf("\nTotal trajectories across all cells: %d\n",
                    sum(cell.df$trajectories)))
    }

    invisible(object)
}


#' Get Cell from MADAG by Maximum Vertex
#'
#' @param madag A \code{madag} object.
#' @param max.vertex Integer. The 1-based index of the maximum vertex.
#'
#' @return The cell structure, or NULL if not found.
#'
#' @export
get.cell <- function(madag, max.vertex) {
    if (!inherits(madag, "madag")) {
        stop("madag must be a 'madag' object")
    }

    for (cell in madag$cells) {
        if (cell$max.vertex == max.vertex) {
            return(cell)
        }
    }

    return(NULL)
}


#' Create Cell Summary Data Frame
#'
#' @param madag A \code{madag} object.
#'
#' @return A data frame with one row per cell.
#'
#' @export
cell.summary.df <- function(madag) {
    if (!inherits(madag, "madag")) {
        stop("madag must be a 'madag' object")
    }

    if (length(madag$cells) == 0) {
        return(data.frame())
    }

    data.frame(
        min.vertex = sapply(madag$cells, `[[`, "min.vertex"),
        max.vertex = sapply(madag$cells, `[[`, "max.vertex"),
        min.value = sapply(madag$cells, `[[`, "min.value"),
        max.value = sapply(madag$cells, `[[`, "max.value"),
        support.size = sapply(madag$cells, function(c) length(c$support)),
        n.trajectories = sapply(madag$cells, `[[`, "n.trajectories"),
        enumerated = sapply(madag$cells, `[[`, "explicitly.enumerated"),
        n.bottlenecks = sapply(madag$cells, function(c) length(c$bottlenecks)),
        stringsAsFactors = FALSE
    )
}
