# Exact graph level-set merge trees for canonical basin complexes.

.build.plateau.graph <- function(adj.list, field, graph.component) {
    n <- length(adj.list)
    plateau.of.vertex <- integer(n)
    plateau.vertices <- list()
    plateau.level <- numeric()
    plateau.component <- integer()
    n.plateaus <- 0L

    for (start in seq_len(n)) {
        if (plateau.of.vertex[[start]] != 0L) {
            next
        }
        n.plateaus <- n.plateaus + 1L
        queue <- start
        plateau.of.vertex[[start]] <- n.plateaus
        cursor <- 1L
        while (cursor <= length(queue)) {
            vertex <- queue[[cursor]]
            cursor <- cursor + 1L
            same.level <- adj.list[[vertex]][
                field[adj.list[[vertex]]] == field[[vertex]]
            ]
            unseen <- same.level[plateau.of.vertex[same.level] == 0L]
            if (length(unseen) > 0L) {
                unseen <- sort(unique(unseen))
                plateau.of.vertex[unseen] <- n.plateaus
                queue <- c(queue, unseen)
            }
        }
        vertices <- sort(unique(as.integer(queue)))
        plateau.vertices[[n.plateaus]] <- vertices
        plateau.level[[n.plateaus]] <- field[[start]]
        plateau.component[[n.plateaus]] <- graph.component[[start]]
    }

    plateau.adj <- vector("list", n.plateaus)
    for (vertex in seq_len(n)) {
        from <- plateau.of.vertex[[vertex]]
        to <- plateau.of.vertex[adj.list[[vertex]]]
        to <- sort(unique(to[to != from]))
        plateau.adj[[from]] <- sort(unique(c(plateau.adj[[from]], to)))
    }
    representatives <- as.integer(vapply(
        plateau.vertices,
        min,
        integer(1)
    ))
    ids <- sprintf("plateau_%08d", representatives)
    list(
        id = ids,
        vertices = plateau.vertices,
        representative = representatives,
        level = as.numeric(plateau.level),
        component = as.integer(plateau.component),
        adj.list = plateau.adj,
        plateau.of.vertex = plateau.of.vertex
    )
}

.merge.tree.elder.order <- function(indices, basins, direction) {
    births <- vapply(
        indices,
        function(index) basins[[index]]$birth.level,
        numeric(1)
    )
    representatives <- vapply(
        indices,
        function(index) basins[[index]]$extremum.vertex,
        integer(1)
    )
    if (direction == "max") {
        order(-births, representatives)
    } else {
        order(births, representatives)
    }
}

.build.level.set.tree <- function(plateaus,
                                  field,
                                  direction,
                                  normalized.mass) {
    n.plateaus <- length(plateaus$id)
    process.order <- if (direction == "max") {
        order(-plateaus$level, plateaus$representative)
    } else {
        order(plateaus$level, plateaus$representative)
    }
    parent <- seq_len(n.plateaus)
    active <- rep(FALSE, n.plateaus)
    current.basin <- rep(NA_integer_, n.plateaus)
    basins <- list()
    events <- list()
    assignment.index <- rep(NA_integer_, length(field))

    find.root <- function(node) {
        root <- node
        while (parent[[root]] != root) {
            root <- parent[[root]]
        }
        while (parent[[node]] != node) {
            next.node <- parent[[node]]
            parent[[node]] <<- root
            node <- next.node
        }
        root
    }

    for (plateau.index in process.order) {
        neighbor.roots <- unique(vapply(
            plateaus$adj.list[[plateau.index]][
                active[plateaus$adj.list[[plateau.index]]]
            ],
            find.root,
            integer(1)
        ))

        if (length(neighbor.roots) == 0L) {
            basin.index <- length(basins) + 1L
            representative <- plateaus$representative[[plateau.index]]
            basins[[basin.index]] <- list(
                basin.id = .basin.direction.id(direction, representative),
                extremum.id =
                    .basin.extremum.id(direction, representative),
                direction = direction,
                extremum.vertex = representative,
                extremum.value = plateaus$level[[plateau.index]],
                plateau.id = plateaus$id[[plateau.index]],
                plateau.vertices =
                    plateaus$vertices[[plateau.index]],
                component = plateaus$component[[plateau.index]],
                birth.level = plateaus$level[[plateau.index]],
                death.level = NA_real_,
                persistence = NA_real_,
                parent.basin.id = NA_character_,
                raw.vertices = integer()
            )
            survivor <- basin.index
        } else {
            candidates <- unique(current.basin[neighbor.roots])
            elder.order <- .merge.tree.elder.order(
                candidates,
                basins,
                direction
            )
            survivor <- candidates[[elder.order[[1L]]]]
            losers <- candidates[elder.order[-1L]]
            for (loser in losers) {
                birth <- basins[[loser]]$birth.level
                death <- plateaus$level[[plateau.index]]
                persistence <- if (direction == "max") {
                    birth - death
                } else {
                    death - birth
                }
                basins[[loser]]$death.level <- death
                basins[[loser]]$persistence <- persistence
                basins[[loser]]$parent.basin.id <-
                    basins[[survivor]]$basin.id
                basins[[survivor]]$raw.vertices <- sort(unique(c(
                    basins[[survivor]]$raw.vertices,
                    basins[[loser]]$raw.vertices
                )))
                events[[length(events) + 1L]] <- list(
                    event.id = sprintf(
                        "merge_%s_%08d",
                        direction,
                        length(events) + 1L
                    ),
                    direction = direction,
                    losing.basin.id = basins[[loser]]$basin.id,
                    surviving.basin.id = basins[[survivor]]$basin.id,
                    merge.plateau.id = plateaus$id[[plateau.index]],
                    merge.vertices =
                        plateaus$vertices[[plateau.index]],
                    merge.level = death,
                    birth.level = birth,
                    death.level = death,
                    persistence = persistence,
                    event.status = if (
                        birth == basins[[survivor]]$birth.level
                    ) {
                        "elder_tie"
                    } else {
                        "ok"
                    }
                )
            }
        }

        active[[plateau.index]] <- TRUE
        parent[[plateau.index]] <- plateau.index
        if (length(neighbor.roots) > 0L) {
            for (root in neighbor.roots) {
                parent[[root]] <- plateau.index
            }
        }
        current.basin[[plateau.index]] <- survivor
        plateau.vertices <- plateaus$vertices[[plateau.index]]
        basins[[survivor]]$raw.vertices <- sort(unique(c(
            basins[[survivor]]$raw.vertices,
            plateau.vertices
        )))
        assignment.index[plateau.vertices] <- survivor
    }

    for (basin.index in seq_along(basins)) {
        if (!is.na(basins[[basin.index]]$death.level)) {
            next
        }
        component.vertices <- which(
            plateaus$component[plateaus$plateau.of.vertex] ==
                basins[[basin.index]]$component
        )
        death <- if (direction == "max") {
            min(field[component.vertices])
        } else {
            max(field[component.vertices])
        }
        basins[[basin.index]]$death.level <- death
        basins[[basin.index]]$persistence <- if (direction == "max") {
            basins[[basin.index]]$birth.level - death
        } else {
            death - basins[[basin.index]]$birth.level
        }
    }

    membership <- .basin.membership.from.supports(
        basins,
        length(field),
        direction
    )
    assignment.rows <- lapply(seq_along(assignment.index), function(vertex) {
        basin.index <- assignment.index[[vertex]]
        list(
            vertex = as.integer(vertex),
            direction = direction,
            basin.id = basins[[basin.index]]$basin.id,
            assignment.weight = 1,
            assignment.status = "assigned",
            assignment.policy = "elder_at_merge",
            root.vertex = basins[[basin.index]]$extremum.vertex,
            next.vertex = NA_integer_
        )
    })
    assignment <- .basin.bind.assignment.rows(assignment.rows)
    membership <- .basin.mark.primary.membership(membership, assignment)

    list(
        basins = basins,
        membership = membership,
        assignment = assignment,
        events = events,
        plateaus = plateaus,
        normalized.mass = normalized.mass
    )
}

.bind.merge.tree.extrema <- function(basins) {
    if (length(basins) == 0L) {
        return(.empty.extrema.table())
    }
    structure(
        list(
            extremum.id = vapply(
                basins, `[[`, character(1), "extremum.id"
            ),
            type = vapply(basins, `[[`, character(1), "direction"),
            representative.vertex = as.integer(vapply(
                basins, `[[`, integer(1), "extremum.vertex"
            )),
            extremum.value = as.numeric(vapply(
                basins, `[[`, numeric(1), "extremum.value"
            )),
            plateau.id = vapply(basins, `[[`, character(1), "plateau.id"),
            plateau.vertices = I(lapply(
                basins, `[[`, "plateau.vertices"
            )),
            n.plateau.vertices = as.integer(vapply(
                basins,
                function(basin) length(basin$plateau.vertices),
                integer(1)
            )),
            is.retained = rep(TRUE, length(basins)),
            retention.status = rep("retained", length(basins))
        ),
        class = "data.frame",
        row.names = seq_along(basins)
    )
}

.bind.merge.tree.basins <- function(basins,
                                    membership,
                                    assignment,
                                    normalized.mass) {
    table <- .basin.bind.basin.table(
        basins,
        membership,
        assignment,
        "superlevel_merge_tree",
        normalized.mass
    )
    if (nrow(table) == 0L) {
        return(table)
    }
    table$birth.level <- as.numeric(vapply(
        basins, `[[`, numeric(1), "birth.level"
    ))
    table$death.level <- as.numeric(vapply(
        basins, `[[`, numeric(1), "death.level"
    ))
    table$persistence <- as.numeric(vapply(
        basins, `[[`, numeric(1), "persistence"
    ))
    table$parent.basin.id <- vapply(
        basins, `[[`, character(1), "parent.basin.id"
    )
    table
}

.bind.merge.events <- function(events) {
    if (length(events) == 0L) {
        return(.empty.merge.table())
    }
    structure(
        list(
            event.id = vapply(events, `[[`, character(1), "event.id"),
            direction = vapply(events, `[[`, character(1), "direction"),
            losing.basin.id = vapply(
                events, `[[`, character(1), "losing.basin.id"
            ),
            surviving.basin.id = vapply(
                events, `[[`, character(1), "surviving.basin.id"
            ),
            merge.plateau.id = vapply(
                events, `[[`, character(1), "merge.plateau.id"
            ),
            merge.vertices = I(lapply(events, `[[`, "merge.vertices")),
            merge.level = as.numeric(vapply(
                events, `[[`, numeric(1), "merge.level"
            )),
            birth.level = as.numeric(vapply(
                events, `[[`, numeric(1), "birth.level"
            )),
            death.level = as.numeric(vapply(
                events, `[[`, numeric(1), "death.level"
            )),
            persistence = as.numeric(vapply(
                events, `[[`, numeric(1), "persistence"
            )),
            event.status = vapply(
                events, `[[`, character(1), "event.status"
            )
        ),
        class = "data.frame",
        row.names = seq_along(events)
    )
}

.validate.merge.tree.complex <- function(object) {
    table <- object$basin.table
    merge <- object$merge.table
    if (is.null(merge)) {
        .stop.basin.complex(
            "A successful merge-tree complex requires a merge table.",
            "object$merge.table",
            class = "gflow_basin_schema_error"
        )
    }
    if (any(!is.finite(table$persistence)) ||
        any(table$persistence < 0) ||
        any(!is.finite(table$birth.level)) ||
        any(!is.finite(table$death.level))) {
        .stop.basin.complex(
            "Merge-tree birth, death, and persistence values must be finite with nonnegative persistence.",
            "object$basin.table",
            class = "gflow_basin_schema_error"
        )
    }
    roots <- is.na(table$parent.basin.id)
    expected.roots <- object$graph.input$validation$n.components *
        length(.basin.requested.directions(object$direction))
    if (sum(roots) != expected.roots) {
        .stop.basin.complex(
            "A merge tree must have one root per graph component and requested direction.",
            "object$basin.table$parent.basin.id",
            class = "gflow_basin_schema_error"
        )
    }
    nonroot.ids <- table$basin.id[!roots]
    if (anyDuplicated(merge$losing.basin.id) ||
        !setequal(nonroot.ids, merge$losing.basin.id)) {
        .stop.basin.complex(
            "Every nonroot basin must have exactly one merge event.",
            "object$merge.table$losing.basin.id",
            class = "gflow_basin_schema_error"
        )
    }
    parent.index <- match(table$parent.basin.id[!roots], table$basin.id)
    if (anyNA(parent.index) ||
        any(table$type[parent.index] != table$type[!roots])) {
        .stop.basin.complex(
            "Every nonroot basin parent must exist in the same direction.",
            "object$basin.table$parent.basin.id",
            class = "gflow_basin_schema_error"
        )
    }
    if (any(object$extrema$n.plateau.vertices !=
        lengths(object$extrema$plateau.vertices)) ||
        any(!vapply(
            seq_len(nrow(object$extrema)),
            function(index) {
                object$extrema$representative.vertex[[index]] %in%
                    object$extrema$plateau.vertices[[index]]
            },
            logical(1)
        ))) {
        .stop.basin.complex(
            "Extremum plateau metadata is inconsistent.",
            "object$extrema",
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}

.create.superlevel.merge.tree.complex <- function(direction,
                                                  graph.input,
                                                  field,
                                                  parameters) {
    plateaus <- .build.plateau.graph(
        graph.input$adj.list,
        field$construction.values,
        graph.input$validation$component
    )
    requested <- .basin.requested.directions(direction)
    trees <- lapply(requested, function(current) {
        .build.level.set.tree(
            plateaus,
            field$construction.values,
            current,
            field$vertex.mass.normalized
        )
    })
    names(trees) <- requested
    basins <- unlist(lapply(trees, `[[`, "basins"), recursive = FALSE)
    names(basins) <- NULL
    membership <- do.call(rbind, lapply(trees, `[[`, "membership"))
    assignment <- do.call(rbind, lapply(trees, `[[`, "assignment"))
    events <- unlist(lapply(trees, `[[`, "events"), recursive = FALSE)
    names(events) <- NULL
    row.names(membership) <- seq_len(nrow(membership))
    row.names(assignment) <- seq_len(nrow(assignment))

    .new.basin.complex(
        method = "superlevel_merge_tree",
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = list(
            package.version = .basin.package.version(),
            method.backend = ".build.level.set.tree",
            construction = list(
                status = "complete",
                completed = TRUE,
                membership.weight.policy = "equal_share_hierarchical_support",
                assignment.policy = "elder_at_merge",
                plateau.policy = "connected_exact_equal_field"
            ),
            refinement.stages = .empty.refinement.stages(),
            random.seed = NULL
        ),
        status = "ok",
        extrema = .bind.merge.tree.extrema(basins),
        basin.table = .bind.merge.tree.basins(
            basins,
            membership,
            assignment,
            field$vertex.mass.normalized
        ),
        membership = membership,
        assignment = assignment,
        merge.table = .bind.merge.events(events),
        diagnostics = .empty.diagnostics.table(),
        raw = list(plateau.graph = plateaus)
    )
}
