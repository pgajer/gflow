# Canonical post-construction basin refinement.

.basin.refinement.applicability <- function() {
    list(
        relative.value = c(
            "trajectory_flow",
            "geodesic_reachability",
            "rtcb"
        ),
        maxima.clustering = c(
            "trajectory_flow",
            "geodesic_reachability",
            "rtcb"
        ),
        minima.clustering = c(
            "trajectory_flow",
            "geodesic_reachability",
            "rtcb"
        ),
        geometric.filter = c(
            "trajectory_flow",
            "geodesic_reachability",
            "rtcb"
        ),
        support.filter = c(
            "trajectory_flow",
            "superlevel_merge_tree",
            "geodesic_reachability",
            "rtcb"
        ),
        expansion = c("geodesic_reachability", "rtcb")
    )
}

.basin.retained.counts <- function(object) {
    table <- object$basin.table
    c(
        max = sum(table$retained & table$type == "max"),
        min = sum(table$retained & table$type == "min")
    )
}

.basin.retained.membership.count <- function(object) {
    retained <- object$basin.table$basin.id[object$basin.table$retained]
    sum(object$membership$basin.id %in% retained)
}

.basin.refinement.snapshot <- function(object) {
    table <- object$basin.table
    table[
        ,
        c(
            "basin.id",
            "type",
            "retained",
            "retained.support.size",
            "retained.support.mass",
            "primary.support.size",
            "primary.support.mass",
            "retention.status"
        ),
        drop = FALSE
    ]
}

.basin.refinement.record <- function(stage,
                                     status,
                                     parameters,
                                     before,
                                     after,
                                     removed = character(),
                                     warnings = character(),
                                     diagnostics = list()) {
    before.counts <- .basin.retained.counts(before)
    after.counts <- .basin.retained.counts(after)
    structure(
        list(
            stage = stage,
            status = status,
            parameters = I(list(parameters)),
            input.basin.count.max = as.integer(before.counts[["max"]]),
            input.basin.count.min = as.integer(before.counts[["min"]]),
            output.basin.count.max = as.integer(after.counts[["max"]]),
            output.basin.count.min = as.integer(after.counts[["min"]]),
            input.membership.count =
                as.integer(.basin.retained.membership.count(before)),
            output.membership.count =
                as.integer(.basin.retained.membership.count(after)),
            retained.basin.ids = I(list(
                after$basin.table$basin.id[after$basin.table$retained]
            )),
            removed.basin.ids = I(list(as.character(removed))),
            warnings = I(list(as.character(warnings))),
            diagnostics = I(list(diagnostics)),
            summary.snapshot = I(list(.basin.refinement.snapshot(after)))
        ),
        class = "data.frame",
        row.names = 1L
    )
}

.basin.set.retained <- function(object,
                                retain,
                                stage,
                                merged.into = NULL) {
    table <- object$basin.table
    currently.retained <- table$retained
    remove <- currently.retained & !retain
    table$retained[remove] <- FALSE
    table$retention.status[remove] <- if (is.null(merged.into)) {
        paste0("filtered_by_", stage)
    } else {
        paste0("merged_into_", merged.into[remove])
    }
    table$retained.support.vertices[remove] <-
        rep(list(integer()), sum(remove))
    table$retained.support.size[remove] <- 0L
    table$retained.support.mass[remove] <- if (
        is.null(object$field$vertex.mass.normalized)
    ) {
        NA_real_
    } else {
        0
    }
    object$basin.table <- table

    removed.ids <- table$basin.id[remove]
    if (length(removed.ids) > 0L) {
        assignment.remove <- object$assignment$basin.id %in% removed.ids &
            object$assignment$assignment.status == "assigned"
        object$assignment$assignment.status[assignment.remove] <- "filtered"
        object$assignment$assignment.weight[assignment.remove] <- 0
        object$assignment$root.vertex[assignment.remove] <- NA_integer_
        object$assignment$next.vertex[assignment.remove] <- NA_integer_
    }
    object
}

.basin.refresh.extrema <- function(object) {
    retained <- object$basin.table$retained[
        match(object$extrema$extremum.id, object$basin.table$extremum.id)
    ]
    status <- object$basin.table$retention.status[
        match(object$extrema$extremum.id, object$basin.table$extremum.id)
    ]
    mapped <- !is.na(retained)
    object$extrema$is.retained[mapped] <- retained[mapped]
    object$extrema$retention.status[mapped] <- status[mapped]
    object
}

.basin.refresh.primary.support <- function(object) {
    table <- object$basin.table
    mass <- object$field$vertex.mass.normalized
    for (index in seq_len(nrow(table))) {
        assigned <- object$assignment$vertex[
            object$assignment$basin.id == table$basin.id[[index]] &
                object$assignment$assignment.status == "assigned"
        ]
        assigned <- sort(unique(as.integer(assigned)))
        table$primary.support.vertices[[index]] <- assigned
        table$primary.support.size[[index]] <- length(assigned)
        direction.has.assignment <- any(
            object$assignment$direction == table$type[[index]] &
                object$assignment$assignment.status == "assigned"
        )
        table$primary.support.mass[[index]] <- if (
            is.null(mass) || !direction.has.assignment
        ) {
            NA_real_
        } else {
            sum(mass[assigned])
        }
        table$assignment.status[[index]] <- if (!direction.has.assignment) {
            "not_assigned"
        } else if (length(assigned) == 0L) {
            "empty"
        } else {
            "assigned"
        }
    }
    object$basin.table <- table
    object$membership <- .basin.mark.primary.membership(
        object$membership,
        object$assignment
    )
    .basin.refresh.extrema(object)
}

.refine.relative.value <- function(object, params) {
    table <- object$basin.table
    mean.field <- mean(object$field$construction.values)
    relative <- table$extremum.value / mean.field
    retain <- table$retained
    max.rows <- table$type == "max" & table$retained
    min.rows <- table$type == "min" & table$retained
    retain[max.rows] <- !is.na(relative[max.rows]) &
        relative[max.rows] >= params$min.relative.value.max
    retain[min.rows] <- !is.na(relative[min.rows]) &
        relative[min.rows] <= params$max.relative.value.min
    .basin.set.retained(object, retain, "relative.value")
}

.basin.overlap.distance <- function(left, right) {
    denominator <- min(length(left), length(right))
    if (denominator == 0L) {
        return(1)
    }
    1 - length(intersect(left, right)) / denominator
}

.basin.overlap.components <- function(supports, threshold) {
    n <- length(supports)
    if (n == 0L) {
        return(integer())
    }
    adjacency <- vector("list", n)
    if (n > 1L) {
        for (left in seq_len(n - 1L)) {
            for (right in (left + 1L):n) {
                distance <- .basin.overlap.distance(
                    supports[[left]],
                    supports[[right]]
                )
                if (distance < threshold) {
                    adjacency[[left]] <- c(adjacency[[left]], right)
                    adjacency[[right]] <- c(adjacency[[right]], left)
                }
            }
        }
    }
    component <- integer(n)
    n.components <- 0L
    for (start in seq_len(n)) {
        if (component[[start]] != 0L) {
            next
        }
        n.components <- n.components + 1L
        queue <- start
        component[[start]] <- n.components
        cursor <- 1L
        while (cursor <= length(queue)) {
            vertex <- queue[[cursor]]
            cursor <- cursor + 1L
            unseen <- adjacency[[vertex]][
                component[adjacency[[vertex]]] == 0L
            ]
            if (length(unseen) > 0L) {
                unseen <- unique(unseen)
                component[unseen] <- n.components
                queue <- c(queue, unseen)
            }
        }
    }
    component
}

.refine.cluster.direction <- function(object, params, direction) {
    table <- object$basin.table
    active <- which(table$retained & table$type == direction)
    if (length(active) < 2L) {
        return(object)
    }
    components <- .basin.overlap.components(
        table$retained.support.vertices[active],
        params$overlap.threshold
    )
    merged.into <- rep(NA_character_, nrow(table))
    retain <- table$retained
    for (component in unique(components)) {
        members <- active[components == component]
        if (length(members) < 2L) {
            next
        }
        representative.order <- if (direction == "max") {
            order(
                -table$extremum.value[members],
                table$extremum.vertex[members]
            )
        } else {
            order(
                table$extremum.value[members],
                table$extremum.vertex[members]
            )
        }
        representative <- members[[representative.order[[1L]]]]
        losers <- setdiff(members, representative)
        merged.vertices <- sort(unique(unlist(
            table$retained.support.vertices[members],
            use.names = FALSE
        )))
        table$retained.support.vertices[[representative]] <- merged.vertices
        table$retained.support.size[[representative]] <-
            length(merged.vertices)
        table$retained.support.mass[[representative]] <-
            .basin.support.mass(
                merged.vertices,
                object$field$vertex.mass.normalized
            )
        retain[losers] <- FALSE
        merged.into[losers] <- table$basin.id[[representative]]

        remap <- object$assignment$basin.id %in% table$basin.id[losers] &
            object$assignment$assignment.status == "assigned"
        object$assignment$basin.id[remap] <-
            table$basin.id[[representative]]
        object$assignment$root.vertex[remap] <-
            table$extremum.vertex[[representative]]
    }
    object$basin.table <- table
    .basin.set.retained(
        object,
        retain,
        paste0(direction, ".clustering"),
        merged.into
    )
}

.basin.vertex.geometry <- function(object, hop.k) {
    adj <- object$graph.input$adj.list
    edge <- object$graph.input$edge.length.list
    n <- object$n.vertices
    mean.neighbor <- vapply(seq_len(n), function(vertex) {
        if (length(edge[[vertex]]) == 0L) {
            Inf
        } else {
            mean(edge[[vertex]])
        }
    }, numeric(1))
    mean.hop <- vapply(seq_len(n), function(vertex) {
        compute_mean_hopk_distance(vertex, adj, edge, hop.k)
    }, numeric(1))
    degree <- lengths(adj)
    list(
        mean.neighbor.percentile = vapply(
            mean.neighbor,
            function(value) sum(mean.neighbor <= value) / n,
            numeric(1)
        ),
        mean.hop.percentile = vapply(
            mean.hop,
            function(value) sum(mean.hop <= value) / n,
            numeric(1)
        ),
        degree.percentile = vapply(
            degree,
            function(value) sum(degree >= value) / n,
            numeric(1)
        )
    )
}

.refine.geometric <- function(object, params) {
    geometry <- .basin.vertex.geometry(object, params$hop.k)
    table <- object$basin.table
    vertices <- table$extremum.vertex
    retain <- table$retained
    active <- which(retain)
    retain[active] <-
        geometry$mean.neighbor.percentile[vertices[active]] <
            params$mean.neighbor.distance.percentile.max &
        geometry$mean.hop.percentile[vertices[active]] <
            params$mean.hop.distance.percentile.max &
        geometry$degree.percentile[vertices[active]] <=
            params$degree.percentile.max
    .basin.set.retained(object, retain, "geometric.filter")
}

.basin.trajectory.counts <- function(object) {
    counts <- setNames(integer(nrow(object$basin.table)), object$basin.table$basin.id)
    if (object$method != "trajectory_flow" ||
        is.null(object$raw$legacy.object$trajectories)) {
        return(counts)
    }
    legacy <- object$raw$legacy.object
    if (!is.null(object$trajectory.forest$plateau)) {
        for (trajectory in legacy$trajectories) {
            endpoints <- c(
                max = trajectory$end.vertex,
                min = trajectory$start.vertex
            )
            for (direction in names(endpoints)) {
                row <- which(
                    object$assignment$direction == direction &
                        object$assignment$vertex == endpoints[[direction]] &
                        object$assignment$assignment.weight > 0
                )
                if (length(row) == 1L) {
                    basin.id <- object$assignment$basin.id[[row]]
                    counts[[basin.id]] <- counts[[basin.id]] + 1L
                }
            }
        }
        return(counts)
    }
    max.ids <- object$basin.table$basin.id[
        object$basin.table$type == "max"
    ]
    min.ids <- object$basin.table$basin.id[
        object$basin.table$type == "min"
    ]
    for (trajectory in legacy$trajectories) {
        max.index <- trajectory$end.basin.idx
        min.index <- trajectory$start.basin.idx
        if (!is.null(max.index) && !is.na(max.index) &&
            max.index >= 1L && max.index <= length(max.ids)) {
            counts[[max.ids[[max.index]]]] <-
                counts[[max.ids[[max.index]]]] + 1L
        }
        if (!is.null(min.index) && !is.na(min.index) &&
            min.index >= 1L && min.index <= length(min.ids)) {
            counts[[min.ids[[min.index]]]] <-
                counts[[min.ids[[min.index]]]] + 1L
        }
    }
    counts
}

.refine.support <- function(object, params) {
    table <- object$basin.table
    retain <- table$retained
    active <- which(retain)
    trajectory.counts <- .basin.trajectory.counts(object)
    mass.ok <- if (params$min.basin.mass == 0) {
        rep(TRUE, length(active))
    } else {
        !is.na(table$retained.support.mass[active]) &
            table$retained.support.mass[active] >= params$min.basin.mass
    }
    retain[active] <-
        table$retained.support.size[active] >= params$min.basin.size &
        trajectory.counts[table$basin.id[active]] >=
            params$min.trajectory.count &
        mass.ok
    .basin.set.retained(object, retain, "support.filter")
}

.refine.expansion <- function(object, params) {
    table <- object$basin.table
    for (direction in .basin.requested.directions(object$direction)) {
        active <- which(table$retained & table$type == direction)
        if (length(active) == 0L) {
            next
        }
        supports <- table$retained.support.vertices[active]
        names(supports) <- table$basin.id[active]
        expanded <- expand.basins.to.cover(
            basins.vertices.list = supports,
            adj.list = object$graph.input$adj.list,
            weight.list = object$graph.input$edge.length.list,
            n.vertices = object$n.vertices
        )
        for (index in active) {
            vertices <- sort(unique(as.integer(
                expanded[[table$basin.id[[index]]]]
            )))
            table$retained.support.vertices[[index]] <- vertices
            table$retained.support.size[[index]] <- length(vertices)
            table$retained.support.mass[[index]] <- .basin.support.mass(
                vertices,
                object$field$vertex.mass.normalized
            )
        }
    }
    object$basin.table <- table
    object
}

.apply.basin.refinements <- function(object) {
    if (object$status != "ok") {
        return(object)
    }
    simplify <- object$parameters$simplify.params
    applicability <- .basin.refinement.applicability()
    records <- list()
    stage.functions <- list(
        relative.value = .refine.relative.value,
        maxima.clustering = function(object, params) {
            .refine.cluster.direction(object, params, "max")
        },
        minima.clustering = function(object, params) {
            .refine.cluster.direction(object, params, "min")
        },
        geometric.filter = .refine.geometric,
        support.filter = .refine.support,
        expansion = .refine.expansion
    )

    for (stage in names(stage.functions)) {
        before <- object
        params <- simplify[[stage]]
        if (!object$method %in% applicability[[stage]]) {
            status <- "not_applicable"
            removed <- character()
        } else if (!params$enabled) {
            status <- "disabled"
            removed <- character()
        } else {
            before.ids <- object$basin.table$basin.id[
                object$basin.table$retained
            ]
            object <- stage.functions[[stage]](object, params)
            object <- .basin.refresh.primary.support(object)
            after.ids <- object$basin.table$basin.id[
                object$basin.table$retained
            ]
            removed <- setdiff(before.ids, after.ids)
            status <- "completed"
        }
        records[[length(records) + 1L]] <- .basin.refinement.record(
            stage,
            status,
            params,
            before,
            object,
            removed
        )
    }
    object$provenance$refinement.stages <- do.call(rbind, records)
    row.names(object$provenance$refinement.stages) <-
        seq_len(nrow(object$provenance$refinement.stages))
    completed <- vapply(
        records,
        function(record) identical(record$status[[1L]], "completed"),
        logical(1)
    )
    object$provenance$allocation <- list(
        raw.current = !any(completed),
        reason = if (any(completed)) {
            "refinement_completed_without_raw_allocation_recompute"
        } else {
            "raw_membership_allocation_current"
        }
    )
    external <- .basin.attach.external.vertex.ids(
        extrema = object$extrema,
        basin.table = object$basin.table,
        membership = object$membership,
        assignment = object$assignment,
        merge.table = object$merge.table,
        vertex.id = object$graph.input$vertex.id
    )
    object$extrema <- external$extrema
    object$basin.table <- external$basin.table
    object$membership <- external$membership
    object$assignment <- external$assignment
    object["merge.table"] <- list(external$merge.table)
    .validate.basin.complex(object)
    object
}
