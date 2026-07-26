# Exact connected-plateau support for canonical CLOSEST trajectory flow.

.trajectory.plateau.boundary.edges <- function(plateaus,
                                               adj.list,
                                               edge.length.list) {
    rows <- vector("list", sum(lengths(adj.list)))
    cursor <- 0L
    for (source in seq_along(adj.list)) {
        source.plateau <- plateaus$plateau.of.vertex[[source]]
        for (edge.index in seq_along(adj.list[[source]])) {
            target <- adj.list[[source]][[edge.index]]
            target.plateau <- plateaus$plateau.of.vertex[[target]]
            if (source.plateau == target.plateau) {
                next
            }
            cursor <- cursor + 1L
            rows[[cursor]] <- list(
                source.plateau = as.integer(source.plateau),
                target.plateau = as.integer(target.plateau),
                source.vertex = as.integer(source),
                target.vertex = as.integer(target),
                edge.length = as.numeric(
                    edge.length.list[[source]][[edge.index]]
                )
            )
        }
    }
    if (cursor == 0L) {
        return(data.frame(
            source.plateau = integer(),
            target.plateau = integer(),
            source.vertex = integer(),
            target.vertex = integer(),
            edge.length = numeric()
        ))
    }
    rows <- rows[seq_len(cursor)]
    data.frame(
        source.plateau = as.integer(vapply(
            rows, `[[`, integer(1), "source.plateau"
        )),
        target.plateau = as.integer(vapply(
            rows, `[[`, integer(1), "target.plateau"
        )),
        source.vertex = as.integer(vapply(
            rows, `[[`, integer(1), "source.vertex"
        )),
        target.vertex = as.integer(vapply(
            rows, `[[`, integer(1), "target.vertex"
        )),
        edge.length = as.numeric(vapply(
            rows, `[[`, numeric(1), "edge.length"
        ))
    )
}

.trajectory.plateau.internal.next <- function(adj.list,
                                               vertices,
                                               anchor) {
    n <- length(adj.list)
    in.plateau <- rep(FALSE, n)
    in.plateau[vertices] <- TRUE
    distance <- rep(Inf, n)
    distance[[anchor]] <- 0
    queue <- as.integer(anchor)
    cursor <- 1L
    while (cursor <= length(queue)) {
        vertex <- queue[[cursor]]
        cursor <- cursor + 1L
        neighbors <- sort(adj.list[[vertex]][
            in.plateau[adj.list[[vertex]]] &
                is.infinite(distance[adj.list[[vertex]]])
        ])
        if (length(neighbors) == 0L) {
            next
        }
        distance[neighbors] <- distance[[vertex]] + 1
        queue <- c(queue, neighbors)
    }
    if (any(is.infinite(distance[vertices]))) {
        .stop.basin.complex(
            "A contracted plateau is not connected.",
            "field",
            class = "gflow_basin_backend_error"
        )
    }
    next.vertex <- rep(NA_integer_, n)
    for (vertex in setdiff(vertices, anchor)) {
        candidates <- adj.list[[vertex]][
            in.plateau[adj.list[[vertex]]] &
                distance[adj.list[[vertex]]] == distance[[vertex]] - 1
        ]
        next.vertex[[vertex]] <- min(candidates)
    }
    next.vertex
}

.trajectory.plateau.outgoing.edge <- function(plateau.index,
                                              direction,
                                              plateaus,
                                              boundary.edges,
                                              edge.length.threshold,
                                              long.edge.fallback) {
    candidates <- boundary.edges[
        boundary.edges$source.plateau == plateau.index,
        ,
        drop = FALSE
    ]
    if (nrow(candidates) > 0L) {
        target.level <- plateaus$level[candidates$target.plateau]
        improving <- if (direction == "max") {
            target.level > plateaus$level[[plateau.index]]
        } else {
            target.level < plateaus$level[[plateau.index]]
        }
        candidates <- candidates[improving, , drop = FALSE]
    }
    if (nrow(candidates) == 0L) {
        return(list(
            edge = NULL,
            used.long.edge.fallback = FALSE,
            blocked.by.long.edge = FALSE
        ))
    }
    candidates <- candidates[order(
        candidates$edge.length,
        candidates$target.vertex,
        candidates$source.vertex
    ), , drop = FALSE]
    admissible <- candidates$edge.length <= edge.length.threshold
    if (any(admissible)) {
        return(list(
            edge = candidates[which(admissible)[[1L]], , drop = FALSE],
            used.long.edge.fallback = FALSE,
            blocked.by.long.edge = FALSE
        ))
    }
    if (long.edge.fallback == "forbid") {
        return(list(
            edge = NULL,
            used.long.edge.fallback = FALSE,
            blocked.by.long.edge = TRUE
        ))
    }
    list(
        edge = candidates[1L, , drop = FALSE],
        used.long.edge.fallback = TRUE,
        blocked.by.long.edge = FALSE
    )
}

.trajectory.follow.roots <- function(next.vertex, max.length) {
    n <- length(next.vertex)
    roots <- rep(NA_integer_, n)
    for (start in seq_len(n)) {
        current <- start
        path <- integer()
        for (step in seq_len(max.length)) {
            if (!is.na(roots[[current]])) {
                root <- roots[[current]]
                break
            }
            if (current %in% path) {
                .stop.basin.complex(
                    "Expanded plateau trajectory contains a directed cycle.",
                    "field",
                    class = "gflow_basin_backend_error"
                )
            }
            path <- c(path, current)
            target <- next.vertex[[current]]
            if (is.na(target)) {
                root <- current
                break
            }
            current <- target
            if (step == max.length) {
                .stop.basin.complex(
                    "Expanded plateau trajectory exceeded the maximum length.",
                    "method.params$max.trajectory.length",
                    class = "gflow_basin_backend_error"
                )
            }
        }
        roots[path] <- root
    }
    roots
}

.trajectory.paths.from.next <- function(next.vertex, roots, direction) {
    lapply(seq_along(next.vertex), function(start) {
        vertices <- as.integer(start)
        current <- start
        while (!is.na(next.vertex[[current]])) {
            current <- next.vertex[[current]]
            vertices <- c(vertices, current)
        }
        list(
            direction = direction,
            vertices = vertices,
            start.vertex = as.integer(start),
            end.vertex = as.integer(roots[[start]])
        )
    })
}

.trajectory.plateau.direction <- function(direction,
                                          plateaus,
                                          boundary.edges,
                                          adj.list,
                                          edge.length.threshold,
                                          long.edge.fallback,
                                          max.trajectory.length) {
    n <- length(plateaus$plateau.of.vertex)
    next.vertex <- rep(NA_integer_, n)
    used.long.edge.fallback <- rep(FALSE, n)
    blocked.by.long.edge <- rep(FALSE, n)
    outgoing.edges <- vector("list", length(plateaus$id))

    for (plateau.index in seq_along(plateaus$id)) {
        selection <- .trajectory.plateau.outgoing.edge(
            plateau.index,
            direction,
            plateaus,
            boundary.edges,
            edge.length.threshold,
            long.edge.fallback
        )
        edge <- selection$edge
        anchor <- if (is.null(edge)) {
            plateaus$representative[[plateau.index]]
        } else {
            edge$source.vertex[[1L]]
        }
        internal.next <- .trajectory.plateau.internal.next(
            adj.list,
            plateaus$vertices[[plateau.index]],
            anchor
        )
        plateau.vertices <- plateaus$vertices[[plateau.index]]
        next.vertex[plateau.vertices] <- internal.next[plateau.vertices]
        if (!is.null(edge)) {
            next.vertex[[anchor]] <- edge$target.vertex[[1L]]
            used.long.edge.fallback[[anchor]] <-
                selection$used.long.edge.fallback
            outgoing.edges[[plateau.index]] <- edge
        } else if (selection$blocked.by.long.edge) {
            blocked.by.long.edge[[anchor]] <- TRUE
        }
    }

    roots <- .trajectory.follow.roots(
        next.vertex,
        max.trajectory.length
    )
    root.vertices <- sort(unique(roots))
    basins <- lapply(root.vertices, function(root) {
        list(
            basin.id = .basin.direction.id(direction, root),
            extremum.id = .basin.extremum.id(direction, root),
            direction = direction,
            extremum.vertex = as.integer(root),
            extremum.value = as.numeric(
                plateaus$level[plateaus$plateau.of.vertex[[root]]]
            ),
            raw.vertices = as.integer(which(roots == root))
        )
    })
    assignment.index <- match(roots, root.vertices)
    names(outgoing.edges) <- plateaus$id

    list(
        basins = basins,
        assignment.index = as.integer(assignment.index),
        next.vertex = as.integer(next.vertex),
        root.vertex = as.integer(roots),
        used.long.edge.fallback = used.long.edge.fallback,
        blocked.by.long.edge = blocked.by.long.edge,
        outgoing.edges = outgoing.edges,
        trajectories = .trajectory.paths.from.next(
            next.vertex,
            roots,
            direction
        )
    )
}

.trajectory.plateau.long.edge.telemetry <- function(results,
                                                    requested,
                                                    policy) {
    n <- length(results[[1L]]$next.vertex)
    empty <- rep(FALSE, n)
    max.result <- if ("max" %in% requested) {
        results[["max"]]
    } else {
        NULL
    }
    min.result <- if ("min" %in% requested) {
        results[["min"]]
    } else {
        NULL
    }
    next.up.used <- if (is.null(max.result)) {
        empty
    } else {
        max.result$used.long.edge.fallback
    }
    next.down.used <- if (is.null(min.result)) {
        empty
    } else {
        min.result$used.long.edge.fallback
    }
    next.up.blocked <- if (is.null(max.result)) {
        empty
    } else {
        max.result$blocked.by.long.edge
    }
    next.down.blocked <- if (is.null(min.result)) {
        empty
    } else {
        min.result$blocked.by.long.edge
    }
    n.used <- sum(next.up.used) + sum(next.down.used)
    list(
        policy = policy,
        attention.required =
            policy == "allow_and_flag" && n.used > 0L,
        next.up.used = next.up.used,
        next.down.used = next.down.used,
        next.up.blocked = next.up.blocked,
        next.down.blocked = next.down.blocked,
        n.next.up = as.integer(sum(next.up.used)),
        n.next.down = as.integer(sum(next.down.used)),
        n.next.up.blocked = as.integer(sum(next.up.blocked)),
        n.next.down.blocked = as.integer(sum(next.down.blocked)),
        trajectory.step.counts = integer(),
        max.basin.vertex.counts = integer(),
        min.basin.vertex.counts = integer()
    )
}

.create.connected.plateau.trajectory.complex <- function(direction,
                                                         graph.input,
                                                         field,
                                                         parameters,
                                                         backend.result,
                                                         warnings) {
    params <- parameters$method.params
    plateaus <- .build.plateau.graph(
        graph.input$adj.list,
        field$input.values,
        graph.input$validation$component
    )
    plateaus$input.level <- plateaus$level
    plateaus$level <- field$construction.values[
        plateaus$representative
    ]
    boundary.edges <- .trajectory.plateau.boundary.edges(
        plateaus,
        graph.input$adj.list,
        graph.input$edge.length.list
    )
    requested <- .basin.requested.directions(direction)
    results <- lapply(requested, function(current) {
        .trajectory.plateau.direction(
            current,
            plateaus,
            boundary.edges,
            graph.input$adj.list,
            backend.result$edge.length.thld,
            params$long.edge.fallback,
            params$max.trajectory.length
        )
    })
    names(results) <- requested

    basins.by.direction <- lapply(results, `[[`, "basins")
    basins <- unlist(basins.by.direction, recursive = FALSE)
    memberships <- lapply(requested, function(current) {
        .basin.membership.from.supports(
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current
        )
    })
    membership <- do.call(rbind, memberships)
    row.names(membership) <- seq_len(nrow(membership))
    assignments <- lapply(requested, function(current) {
        .basin.assignment.from.membership(
            membership,
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current,
            params$primary.assignment.policy,
            results[[current]]$assignment.index,
            results[[current]]$next.vertex
        )
    })
    assignment <- do.call(rbind, assignments)
    row.names(assignment) <- seq_len(nrow(assignment))

    long.edge.telemetry <- .trajectory.plateau.long.edge.telemetry(
        results,
        requested,
        params$long.edge.fallback
    )
    forest <- list(
        modulation = params$modulation,
        plateau.policy = params$plateau.policy,
        plateau = list(
            id = plateaus$id,
            representative.vertex = plateaus$representative,
            vertices = plateaus$vertices,
            field.value = plateaus$level,
            plateau.of.vertex = plateaus$plateau.of.vertex,
            outgoing.edges = lapply(results, `[[`, "outgoing.edges")
        ),
        long.edge.fallback = long.edge.telemetry,
        next.vertex = lapply(results, `[[`, "next.vertex"),
        root.vertex = lapply(results, `[[`, "root.vertex"),
        trajectories = lapply(results, `[[`, "trajectories"),
        trajectory.status = "complete",
        tie.breaking = params$tie.breaking,
        tie.seed = params$tie.seed
    )

    backend.result$next.up <- if ("max" %in% requested) {
        results[["max"]]$next.vertex
    } else {
        rep(NA_integer_, length(graph.input$adj.list))
    }
    backend.result$next.down <- if ("min" %in% requested) {
        results[["min"]]$next.vertex
    } else {
        rep(NA_integer_, length(graph.input$adj.list))
    }
    backend.result$long.edge.fallback <- long.edge.telemetry
    backend.record <- list(
        name = "compute.gfc.trajectory+connected_exact_plateau_contraction",
        object = backend.result
    )
    object <- .basin.success.object(
        method = "trajectory_flow",
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        basins = basins,
        membership = membership,
        assignment = assignment,
        backend = backend.record,
        warnings = warnings,
        trajectory.forest = forest
    )

    for (row in seq_len(nrow(object$extrema))) {
        representative <- object$extrema$representative.vertex[[row]]
        plateau.index <- plateaus$plateau.of.vertex[[representative]]
        object$extrema$plateau.id[[row]] <- plateaus$id[[plateau.index]]
        object$extrema$plateau.vertices[[row]] <-
            plateaus$vertices[[plateau.index]]
        object$extrema$n.plateau.vertices[[row]] <-
            length(plateaus$vertices[[plateau.index]])
    }
    object$provenance$construction$plateau.policy <- params$plateau.policy
    object$provenance$construction$n.plateaus <- length(plateaus$id)
    object$provenance$construction$n.contracted.plateaus <- sum(
        lengths(plateaus$vertices) > 1L
    )
    object$provenance$construction$n.contracted.vertices <- sum(
        pmax(lengths(plateaus$vertices) - 1L, 0L)
    )
    .validate.basin.complex(object)
    object
}
