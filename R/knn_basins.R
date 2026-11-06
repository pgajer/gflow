#' Identify Extrema and Basins Directly from a Data Matrix
#'
#' Find maxima and minima of a supplied outcome on exact, union-symmetrized
#' nearest-neighbor neighborhoods of the rows of `X`. Assign ascent/descent
#' basins, merge them by topological persistence, then independently filter
#' by the number of assigned observations.
#'
#' @param X Finite numeric matrix with at least two rows and one column.
#'   All columns are used as supplied: no PCA, scaling, or normalization.
#' @param y Finite numeric vector with `length(y) == nrow(X)`. May contain
#'   outcomes smoothed elsewhere; this function does not smooth them.
#' @param k Number of other nearest neighbors, an integer in
#'   `1:(nrow(X) - 1)`. Required explicitly because it sets the spatial scale.
#' @param metric Either `"L2"` (Euclidean) or `"Linfty"` (maximum absolute
#'   coordinate difference). This is ordinary Linfty distance, without a
#'   simplex constraint or row normalization.
#' @param direction `"both"`, `"max"`, or `"min"`.
#' @param min.persistence Nonnegative threshold in units of `y`. Basins with
#'   persistence strictly below it merge into their surviving ancestors.
#'   `Inf` merges all nonessential basins. Default zero performs no merges.
#' @param min.basin.size Positive integer. After persistence merging, basins
#'   with fewer assigned rows are filtered out; their final labels are `NA`.
#'   This filter does not reassign rows or alter the persistence hierarchy.
#'
#' @details
#' Each row selects exactly `k` other rows, with distance ties resolved by row
#' index. An undirected edge is present if either endpoint selects the other.
#' Consequently, the symmetrized degree can exceed `k`. The same adjacency
#' defines extrema, flow, and persistence. Duplicate rows of `X` are allowed
#' and remain separate observations, connected by zero-length edges when
#' selected. Ties make results dependent on row order; no random jitter is used.
#'
#' Connected rows with exactly equal `y` form a plateau. A maximum plateau
#' has only lower-valued neighbors outside it; minima reverse the inequality.
#' A component with constant `y` is one essential plateau in each direction.
#' The smallest row index represents a plateau. Each nonterminal plateau
#' chooses its shortest boundary edge to an improving plateau, breaking ties
#' by target row and then source row. Following these choices gives disjoint
#' initial basins, with one label per observation and direction.
#'
#' Persistence uses the exact superlevel (maxima) or sublevel (minima) merge
#' tree. At a merge the higher maximum or lower minimum survives; equal
#' extrema use the smaller representative row as elder. Persistence is peak
#' minus merge level for maxima, or merge level minus trough for minima.
#' One essential extremum per connected component never dies: it has
#' `merge.level = NA` and `persistence = Inf`. This differs from the finite
#' component-range convention of [create.basin.complex()]. Essential basins
#' survive persistence merging but are still subject to the size filter.
#'
#' Size means assigned observation count, not hierarchical superlevel support
#' or geometric volume. A large `k` does not guarantee a minimum basin size.
#' Merged basins are unions of initial flow basins; trajectories are not
#' recomputed. Such unions can be disconnected when an intervening basin
#' survives the persistence threshold. `parent.row` records the topological
#' elder, which need not be an adjacent basin in the initial flow partition.
#' Persistence is a feature-height measure, not a significance test. Feature
#' scaling and a scientifically meaningful metric are the caller's choices.
#'
#' Exact neighbor search costs quadratic time in the number of rows (and
#' linear time in the number of columns per distance). It stores neighbor
#' lists rather than a dense pairwise distance matrix. The merge-tree backend
#' additionally constructs hierarchical supports internally. This interface
#' is intended for data sets of a few thousand observations.
#'
#' @return An object of class `gflow_knn_basins` with:
#' \describe{
#'   \item{max, min}{One result per requested direction; `NULL` otherwise.
#'     Each contains `labels` (final representative-row labels, with `NA`
#'     for size-filtered observations), `initial.labels` (before merging),
#'     `persistence.labels` (after merging, before size filtering), `basins`
#'     (a named list of retained row-index sets), and `extrema`.}
#'   \item{extrema tables}{Each direction's table includes every original
#'     extremum: `row`, `value`, `component`, `plateau.size`, `plateau.rows`,
#'     `local.gap` (height difference to the most competitive neighbor outside
#'     the plateau, or `NA` without a boundary), `merge.level`, `persistence`,
#'     `parent.row` (original elder at death, `NA` for essential extrema),
#'     `essential`, `initial.size`, `persistence.root` (surviving representative
#'     after persistence merging), `merged.size` (assigned count before the
#'     size filter; zero for merged extrema), `retained`, and `status`
#'     (`"retained"`, `"merged"`, or `"filtered_size"`).}
#'   \item{neighborhood}{Directed `knn.index` and `knn.distance` matrices,
#'     `radius` (kth-neighbor distances), union-symmetrized `adj.list` and
#'     parallel `edge.length.list`, and `component` and `n.components`.}
#'   \item{y, row.names, parameters, call}{Supplied field and row names,
#'     resolved settings, and the call. The input matrix is not duplicated
#'     in the result.}
#' }
#'
#' @examples
#' X <- matrix(seq_len(9), ncol = 1)
#' y <- c(0, 4, 3, 1, 2, 3, 2, 1, 0)
#' fit <- knn.basins(X, y, k = 2, metric = "L2",
#'                   min.persistence = 2, min.basin.size = 2)
#' fit
#' fit$max$extrema
#' fit$max$labels
#' fit$min$basins
#' # Minima of y are maxima of -y, with the same geometry.
#' inverted <- knn.basins(X, -y, k = 2, direction = "max")
#' original <- knn.basins(X, y, k = 2, direction = "min")
#' stopifnot(identical(inverted$max$labels, original$min$labels))
#'
#' @seealso [create.basin.complex()], [detect.adaptive.extrema()]
#' @export
knn.basins <- function(X, y, k, metric = c("L2", "Linfty"),
                       direction = c("both", "max", "min"),
                       min.persistence = 0, min.basin.size = 1L) {
    call <- match.call()
    if (!is.matrix(X) || !is.numeric(X) || is.complex(X) ||
        nrow(X) < 2L || ncol(X) < 1L || any(!is.finite(X))) {
        stop("'X' must be a finite numeric matrix with at least two rows and one column.")
    }
    if (!is.null(dim(y))) stop("'y' must be a numeric vector, not an array.")
    y <- .basin.assert.vector(y, "y", nrow(X))
    if (!is.finite(max(y) - min(y))) {
        stop("The range of 'y' overflows; rescale y explicitly.")
    }
    k <- .basin.assert.integer(k, "k", 1L, nrow(X) - 1L)
    metric <- .basin.assert.choice(metric, c("L2", "Linfty"), "metric")
    direction <- .basin.assert.choice(direction, c("both", "max", "min"), "direction")
    min.persistence <- .basin.assert.number(
        min.persistence, "min.persistence", lower = 0, finite = FALSE
    )
    min.basin.size <- .basin.assert.integer(
        min.basin.size, "min.basin.size", 1L, .Machine$integer.max
    )
    graph <- rcpp_knn_basins_neighbors(X, k, metric == "Linfty")
    graph$radius <- graph$knn.distance[, k]
    plateaus <- .build.plateau.graph(graph$adj.list, y, graph$component)
    requested <- if (direction == "both") c("max", "min") else direction
    result <- list(max = NULL, min = NULL)
    for (orientation in requested) {
        result[[orientation]] <- .knn.basins.direction(
            graph, plateaus, y, orientation, min.persistence, min.basin.size
        )
    }
    result$neighborhood <- graph
    result$y <- y
    result["row.names"] <- list(rownames(X))
    result$parameters <- list(
        k = k, metric = metric, direction = direction,
        min.persistence = min.persistence, min.basin.size = min.basin.size,
        symmetrization = "union", flow = "nearest_improving_plateau",
        scaling = "none", projection = "none", smoothing = "none",
        plateau.tolerance = 0, distance.ties = "row_index"
    )
    result$call <- call
    class(result) <- "gflow_knn_basins"
    result
}

# Choose a single improving edge per plateau, then follow the resulting DAG
# in field order. Work on flat edge vectors, avoiding a scan of all edges
# separately for every plateau or storing every observation's trajectory.
.knn.basins.initial.labels <- function(graph, plateaus, direction) {
    source <- rep.int(seq_along(graph$adj.list), lengths(graph$adj.list))
    target <- unlist(graph$adj.list, use.names = FALSE)
    distance <- unlist(graph$edge.length.list, use.names = FALSE)
    from <- plateaus$plateau.of.vertex[source]
    to <- plateaus$plateau.of.vertex[target]
    level <- if (direction == "max") plateaus$level else -plateaus$level
    improving <- which(level[to] > level[from])
    ordered <- improving[order(from[improving], distance[improving],
                               target[improving], source[improving])]
    selected <- ordered[!duplicated(from[ordered])]
    next.plateau <- rep(NA_integer_, length(level))
    next.plateau[from[selected]] <- to[selected]
    root <- plateaus$representative
    for (i in order(-level, plateaus$representative)) {
        if (!is.na(next.plateau[[i]])) root[[i]] <- root[[next.plateau[[i]]]]
    }
    root[plateaus$plateau.of.vertex]
}

.knn.basins.direction <- function(graph, plateaus, y, direction,
                                  min.persistence, min.basin.size) {
    tree <- .build.level.set.tree(plateaus, y, direction, NULL)
    basins <- tree$basins
    rows <- vapply(basins, `[[`, integer(1), "extremum.vertex")
    ids <- vapply(basins, `[[`, character(1), "basin.id")
    parent <- match(vapply(basins, `[[`, character(1), "parent.basin.id"), ids)
    essential <- is.na(parent)
    persistence <- vapply(basins, `[[`, numeric(1), "persistence")
    merge.level <- vapply(basins, `[[`, numeric(1), "death.level")
    # Essential extrema have no death; do not use the canonical backend's
    # finite component-range endpoint as evidence for a merge.
    persistence[essential] <- Inf
    merge.level[essential] <- NA_real_
    initial <- .knn.basins.initial.labels(graph, plateaus, direction)
    initial.index <- match(initial, rows)
    if (anyNA(initial.index)) stop("Internal error: flow and merge-tree extrema disagree.")
    next.basin <- rep(NA_integer_, length(rows))
    merge <- !essential & persistence < min.persistence
    next.basin[merge] <- parent[merge]
    survivor <- .trajectory.follow.roots(next.basin, length(rows) + 1L)
    merged <- rows[survivor[initial.index]]
    merged.size <- tabulate(match(merged, rows), nbins = length(rows))
    retained <- !merge & merged.size >= min.basin.size
    labels <- merged
    labels[!retained[match(merged, rows)]] <- NA_integer_
    gap <- vapply(basins, function(basin) {
        p <- plateaus$plateau.of.vertex[[basin$extremum.vertex]]
        neighbors <- plateaus$adj.list[[p]]
        if (!length(neighbors)) return(NA_real_)
        if (direction == "max") {
            plateaus$level[[p]] - max(plateaus$level[neighbors])
        } else {
            min(plateaus$level[neighbors]) - plateaus$level[[p]]
        }
    }, numeric(1))
    extrema <- data.frame(
        row = rows, value = y[rows], component = graph$component[rows],
        plateau.size = lengths(lapply(basins, `[[`, "plateau.vertices")),
        local.gap = gap, merge.level = merge.level, persistence = persistence,
        parent.row = rows[parent], essential = essential,
        initial.size = tabulate(initial.index, nbins = length(rows)),
        persistence.root = rows[survivor], merged.size = merged.size,
        retained = retained,
        status = ifelse(merge, "merged", ifelse(retained, "retained", "filtered_size"))
    )
    extrema$plateau.rows <- I(lapply(basins, `[[`, "plateau.vertices"))
    list(
        labels = labels, initial.labels = initial, persistence.labels = merged,
        basins = split(seq_along(y), factor(labels, levels = rows[retained])),
        extrema = extrema
    )
}

.validate.knn.basins <- function(x) {
    .validate.s3.contract(
        x, "gflow_knn_basins", storage = "list",
        fields = c("max", "min", "neighborhood", "y", "row.names", "parameters", "call")
    )
    if (!is.list(x$parameters) ||
        !all(c("k", "metric", "direction") %in% names(x$parameters))) {
        stop("Invalid 'gflow_knn_basins' object; missing neighborhood parameters.")
    }
    for (direction in c("max", "min")) {
        result <- x[[direction]]
        if (is.null(result)) next
        if (!is.list(result) || length(result$labels) != length(x$y) ||
            !is.data.frame(result$extrema) ||
            !all(c("status", "retained") %in% names(result$extrema))) {
            stop("Invalid 'gflow_knn_basins' object; malformed ", direction, " result.")
        }
    }
    invisible(x)
}

#' @rdname knn.basins
#' @param x A `gflow_knn_basins` result.
#' @param ... Unused.
#' @export
print.gflow_knn_basins <- function(x, ...) {
    .validate.knn.basins(x)
    cat(sprintf("kNN basins: %d rows, k = %d, %s, union-symmetrized neighborhoods\n",
                length(x$y), x$parameters$k, x$parameters$metric))
    for (direction in c("max", "min")) {
        result <- x[[direction]]
        if (is.null(result)) next
        cat(sprintf("  %s: %d original extrema, %d persistence survivors, %d retained; %d unassigned rows\n",
                    direction, nrow(result$extrema),
                    sum(result$extrema$status != "merged"),
                    sum(result$extrema$retained), sum(is.na(result$labels))))
    }
    invisible(x)
}
