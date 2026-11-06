#' Identify Extrema and Basins Directly from a Data Matrix
#'
#' Find maxima and minima of a supplied outcome on exact nearest-neighbor
#' neighborhoods of the rows of `X`, with optional union symmetrization.
#' Assign shortest-improving-edge ascent/descent basins, merge them by
#' persistence or directed escape prominence, then independently filter
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
#' @param symmetrize Logical, default `TRUE`, preserving union-symmetrized
#'   neighborhoods and connected-plateau extrema with level-set persistence.
#'   With `FALSE`, each row is tested against exactly its own `k` nearest
#'   other rows. Basin flow and escape paths follow only outgoing selections;
#'   equal-valued rows are not contracted. See Details for the different
#'   prominence and tie conventions.
#'
#' @details
#' Each row selects exactly `k` other rows, with distance ties resolved by row
#' index. With `symmetrize = TRUE`, an undirected edge is present if either
#' endpoint selects the other; consequently, degree can exceed `k`. With
#' `FALSE`, every outgoing neighborhood contains exactly `k` rows. The same
#' selected adjacency defines extrema, flow, and the mode's prominence.
#' Duplicate rows of `X` are allowed
#' and remain separate observations, connected by zero-length edges when
#' selected. Ties make results dependent on row order; no random jitter is used.
#'
#' In the default symmetric mode, connected rows with exactly equal `y`
#' form a plateau. A maximum plateau
#' has only lower-valued neighbors outside it; minima reverse the inequality.
#' A component with constant `y` is one essential plateau in each direction.
#' The smallest row index represents a plateau. Each nonterminal plateau
#' chooses its shortest boundary edge to an improving plateau, breaking ties
#' by target row and then source row. Following these choices gives disjoint
#' initial basins, with one label per observation and direction.
#'
#' In symmetric mode, persistence uses the exact superlevel (maxima) or sublevel (minima) merge
#' tree. At a merge the higher maximum or lower minimum survives; equal
#' extrema use the smaller representative row as elder. Persistence is peak
#' minus merge level for maxima, or merge level minus trough for minima.
#' One essential extremum per connected component never dies: it has
#' `merge.level = NA` and `persistence = Inf`. This differs from the finite
#' component-range convention of [create.basin.complex()]. Essential basins
#' survive persistence merging but are still subject to the size filter.
#'
#' With `symmetrize = FALSE`, a maximum is a row with no strictly higher
#' selected neighbor; a minimum has no strictly lower selected neighbor.
#' These are weak inequalities: equal-valued neighbors do not disqualify a
#' row, and a flat region can contain many terminal extrema. The query row
#' is implicit in its comparison. Each nonterminal row chooses its shortest
#' outgoing edge to a strictly improving neighbor, breaking distance ties
#' by target row. This is the `CLOSEST` rule used by graph-based trajectory
#' flow, without an edge-length cutoff. Following these successors partitions
#' all rows into initial basins. Equal-valued edges are not flow steps.
#' Choosing the shortest improving edge prevents a larger outcome gain from
#' overriding a nearer improving neighbor, but does not impose a maximum
#' step length or minimize the total path length to a terminal extremum.
#'
#' Directed mode reports **directed escape prominence**, not ordinary
#' connected-component topological persistence. For a maximum `m`, consider
#' all outgoing directed paths from `m` to any strictly higher terminal
#' maximum. `merge.level` is the largest achievable minimum outcome along
#' such a path, and `persistence = y[m] - merge.level`. For minima, minimize
#' the largest outcome along paths to strictly lower terminal minima and
#' use `persistence = merge.level - y[m]`. Escape paths may go both uphill
#' and downhill, including through equal values; they are distinct from
#' the strictly improving basin trajectories. If several target extrema
#' attain the same escape level, `parent.row` selects the most extreme,
#' then the smallest row index. Parents are always strictly more extreme,
#' so these links form an acyclic forest. Equal-height extrema never serve
#' as each other's escape targets.
#'
#' If no strictly more extreme terminal row is reachable, directed prominence
#' is `Inf`, `merge.level` is `NA`, and `essential` is `TRUE`. Several such
#' extrema can occur in one weakly connected component, even at different
#' heights. An entirely flat field has one essential extremum per row in
#' each direction. A weak extremum that can reach a better extremum through
#' a level path can have zero prominence. Threshold equality is retained
#' in both modes. Merging follows parent links of extrema strictly below
#' `min.persistence`; the same independent size filter then applies.
#'
#' Size means assigned observation count, not hierarchical superlevel support
#' or geometric volume. A large `k` does not guarantee a minimum basin size.
#' Merged basins are unions of initial flow basins; trajectories are not
#' recomputed. Such unions can be disconnected when an intervening basin
#' survives the persistence threshold. `parent.row` records the topological
#' elder (symmetric mode) or escape target (directed mode), which need not
#' be an adjacent basin in the initial flow partition.
#' Persistence is a feature-height measure, not a significance test. Feature
#' scaling and a scientifically meaningful metric are the caller's choices.
#'
#' Exact neighbor search costs quadratic time in the number of rows (and
#' linear time in the number of columns per distance). It stores neighbor
#' lists rather than a dense pairwise distance matrix. The merge-tree backend
#' additionally constructs hierarchical supports internally. This interface
#' is intended for data sets of a few thousand observations. Directed escape
#' uses a widest-path search per non-global terminal extremum, with worst-case
#' time `O(m * n * k * log(n * k))` for `m` such extrema and `O(n * k)`
#' workspace; many terminal extrema can therefore be expensive.
#'
#' @return An object of class `gflow_knn_basins` with:
#' \describe{
#'   \item{max, min}{One result per requested direction; `NULL` otherwise.
#'     Each contains `labels` (final representative-row labels, with `NA`
#'     for size-filtered observations), `initial.labels` (before merging),
#'     `persistence.labels` (after merging, before size filtering), `basins`
#'     (a named list of retained row-index sets), and `extrema`. Directed mode
#'     additionally returns `next.vertex`, the initial shortest-improving-edge
#'     successor of each row (`NA` at terminal extrema), unchanged by merging.}
#'   \item{extrema tables}{Each direction's table includes every original
#'     extremum: `row`, `value`, `component`, `plateau.size`, `plateau.rows`,
#'     `local.gap` (height difference to the most competitive neighbor outside
#'     the plateau, or `NA` without a boundary), `merge.level`, `persistence`,
#'     `parent.row` (original elder at death, `NA` for essential extrema),
#'     `essential`, `initial.size`, `persistence.root` (surviving representative
#'     after persistence merging), `merged.size` (assigned count before the
#'     size filter; zero for merged extrema), `retained`, and `status`
#'     (`"retained"`, `"merged"`, or `"filtered_size"`). In directed mode,
#'     `plateau.size` is always one and `plateau.rows` contains the terminal
#'     row alone; `local.gap` uses only that row's selected neighbors, and
#'     `persistence` and `merge.level` describe directed escape as above.}
#'   \item{neighborhood}{Directed `knn.index` and `knn.distance` matrices,
#'     `radius` (kth-neighbor distances), selected `adj.list` and parallel
#'     `edge.length.list`, and `component` and `n.components`. Adjacency is
#'     union-symmetrized by default and outgoing-only when `symmetrize = FALSE`.
#'     Component labels always describe the underlying undirected connectivity
#'     (weak components in directed mode), not directed reachability.}
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
#' # A point selected by a higher row can still be an outgoing-kNN maximum.
#' directed <- knn.basins(matrix(c(0, 1, 3.5)), c(0, 5, 10), k = 1,
#'                        symmetrize = FALSE)
#' directed$max$extrema
#' directed$max$next.vertex
#'
#' @seealso [create.basin.complex()], [detect.adaptive.extrema()]
#' @export
knn.basins <- function(X, y, k, metric = c("L2", "Linfty"),
                       direction = c("both", "max", "min"),
                       min.persistence = 0, min.basin.size = 1L,
                       symmetrize = TRUE) {
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
    symmetrize <- .basin.assert.logical(symmetrize, "symmetrize")
    graph <- rcpp_knn_basins_neighbors(X, k, metric == "Linfty")
    graph$radius <- graph$knn.distance[, k]
    if (symmetrize) {
        plateaus <- .build.plateau.graph(graph$adj.list, y, graph$component)
    } else {
        # Keep the weak-component labels, but use ONLY selected outgoing
        # neighbors for extrema, flow, and directed escape paths.
        for (i in seq_len(nrow(X))) {
            ordered <- order(graph$knn.index[i, ])
            graph$adj.list[[i]] <- graph$knn.index[i, ordered]
            graph$edge.length.list[[i]] <- graph$knn.distance[i, ordered]
        }
        # Singleton records reuse the same shortest-improving-edge flow rule
        # without pooling equal-valued rows' neighborhoods.
        plateaus <- list(plateau.of.vertex = seq_along(y),
                         representative = seq_along(y), level = y)
    }
    requested <- if (direction == "both") c("max", "min") else direction
    result <- list(max = NULL, min = NULL)
    for (orientation in requested) {
        result[[orientation]] <- if (symmetrize) {
            .knn.basins.direction(graph, plateaus, y, orientation,
                                  min.persistence, min.basin.size)
        } else {
            .knn.basins.directed.direction(graph, plateaus, y, orientation,
                                           min.persistence, min.basin.size)
        }
    }
    result$neighborhood <- graph
    result$y <- y
    result["row.names"] <- list(rownames(X))
    result$parameters <- list(
        k = k, metric = metric, direction = direction,
        min.persistence = min.persistence, min.basin.size = min.basin.size,
        symmetrization = if (symmetrize) "union" else "none",
        flow = if (symmetrize) "nearest_improving_plateau" else "nearest_improving_vertex",
        persistence = if (symmetrize) "level_set" else "directed_escape",
        plateau.policy = if (symmetrize) "connected_exact" else "per_row",
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
.knn.basins.initial.flow <- function(graph, plateaus, direction) {
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
    list(labels = root[plateaus$plateau.of.vertex], next.plateau = next.plateau)
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
    initial <- .knn.basins.initial.flow(graph, plateaus, direction)$labels
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
    .knn.basins.assemble(graph, y, rows, parent, persistence, merge.level,
                         initial, lapply(basins, `[[`, "plateau.vertices"), gap,
                         min.persistence, min.basin.size)
}

.knn.basins.directed.direction <- function(graph, plateaus, y, direction,
                                           min.persistence, min.basin.size) {
    field <- if (direction == "max") y else -y
    flow <- .knn.basins.initial.flow(graph, plateaus, direction)
    rows <- which(is.na(flow$next.plateau))
    rows <- rows[order(-field[rows], rows)]
    escape <- rcpp_knn_basins_escape(graph$adj.list, field, rows)
    essential <- is.na(escape$parent)
    persistence <- field[rows] - escape$merge.level
    persistence[essential] <- Inf
    merge.level <- if (direction == "max") escape$merge.level else -escape$merge.level
    gap <- vapply(rows, function(i) field[i] - max(field[graph$adj.list[[i]]]), numeric(1))
    result <- .knn.basins.assemble(
        graph, y, rows, escape$parent, persistence, merge.level,
        flow$labels, lapply(rows, as.integer), gap, min.persistence, min.basin.size
    )
    # Mirrors the graph-based trajectory representation. A missing successor
    # means the per-row extremum test passed, including weak extrema on ties.
    result$next.vertex <- flow$next.plateau
    result
}

.knn.basins.assemble <- function(graph, y, rows, parent, persistence, merge.level,
                                 initial, plateau.rows, gap,
                                 min.persistence, min.basin.size) {
    essential <- is.na(parent)
    initial.index <- match(initial, rows)
    if (anyNA(initial.index)) stop("Internal error: flow and reported extrema disagree.")
    next.basin <- rep(NA_integer_, length(rows))
    merge <- !essential & persistence < min.persistence
    next.basin[merge] <- parent[merge]
    survivor <- .trajectory.follow.roots(next.basin, length(rows) + 1L)
    merged <- rows[survivor[initial.index]]
    merged.size <- tabulate(match(merged, rows), nbins = length(rows))
    retained <- !merge & merged.size >= min.basin.size
    labels <- merged
    labels[!retained[match(merged, rows)]] <- NA_integer_
    extrema <- data.frame(
        row = rows, value = y[rows], component = graph$component[rows],
        plateau.size = lengths(plateau.rows),
        local.gap = gap, merge.level = merge.level, persistence = persistence,
        parent.row = rows[parent], essential = essential,
        initial.size = tabulate(initial.index, nbins = length(rows)),
        persistence.root = rows[survivor], merged.size = merged.size,
        retained = retained,
        status = ifelse(merge, "merged", ifelse(retained, "retained", "filtered_size"))
    )
    extrema$plateau.rows <- I(plateau.rows)
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
    neighborhood <- if (identical(x$parameters$symmetrization, "none")) {
        "directed neighborhoods, directed escape prominence"
    } else "union-symmetrized neighborhoods"
    cat(sprintf("kNN basins: %d rows, k = %d, %s, %s\n",
                length(x$y), x$parameters$k, x$parameters$metric, neighborhood))
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
