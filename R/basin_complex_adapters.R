# Canonical adapters for legacy basin-construction backends.

.basin.requested.directions <- function(direction) {
    if (identical(direction, "both")) c("max", "min") else direction
}

.basin.capture.backend <- function(expr) {
    warnings <- character()
    value <- tryCatch(
        withCallingHandlers(
            expr,
            warning = function(condition) {
                warnings <<- c(warnings, conditionMessage(condition))
                invokeRestart("muffleWarning")
            }
        ),
        error = identity
    )
    list(
        value = value,
        warnings = unique(warnings),
        error = if (inherits(value, "error")) value else NULL
    )
}

.basin.with.local.seed <- function(seed, expr) {
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had.seed) {
        old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if (had.seed) {
            assign(".Random.seed", old.seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    force(expr)
}

.basin.support.mass <- function(vertices, normalized.mass) {
    if (is.null(normalized.mass)) {
        return(NA_real_)
    }
    sum(normalized.mass[vertices])
}

.basin.direction.id <- function(direction, vertex) {
    sprintf("basin_%s_v%08d", direction, as.integer(vertex))
}

.basin.extremum.id <- function(direction, vertex) {
    sprintf("extremum_%s_v%08d", direction, as.integer(vertex))
}

.basin.membership.from.index <- function(membership.by.vertex,
                                         basins,
                                         direction,
                                         source.stage = "construction") {
    rows <- vector("list", length(membership.by.vertex))
    cursor <- 0L
    for (vertex in seq_along(membership.by.vertex)) {
        indices <- unique(as.integer(membership.by.vertex[[vertex]]))
        indices <- indices[
            !is.na(indices) & indices >= 1L & indices <= length(basins)
        ]
        if (length(indices) == 0L) {
            next
        }
        weight <- 1 / length(indices)
        for (index in indices) {
            cursor <- cursor + 1L
            rows[[cursor]] <- list(
                vertex = as.integer(vertex),
                direction = direction,
                basin.id = basins[[index]]$basin.id,
                membership.weight = weight,
                membership.status = "raw",
                source.stage = source.stage,
                is.primary = FALSE
            )
        }
    }
    .basin.bind.membership.rows(rows[seq_len(cursor)])
}

.basin.membership.from.supports <- function(basins, n.vertices, direction) {
    membership.by.vertex <- vector("list", n.vertices)
    for (index in seq_along(basins)) {
        vertices <- basins[[index]]$raw.vertices
        for (vertex in vertices) {
            membership.by.vertex[[vertex]] <- c(
                membership.by.vertex[[vertex]],
                index
            )
        }
    }
    .basin.membership.from.index(
        membership.by.vertex,
        basins,
        direction
    )
}

.basin.bind.membership.rows <- function(rows) {
    if (length(rows) == 0L) {
        return(.empty.membership.table())
    }
    structure(
        list(
            vertex = as.integer(vapply(rows, `[[`, integer(1), "vertex")),
            direction = vapply(rows, `[[`, character(1), "direction"),
            basin.id = vapply(rows, `[[`, character(1), "basin.id"),
            membership.weight = as.numeric(vapply(
                rows, `[[`, numeric(1), "membership.weight"
            )),
            membership.status = vapply(
                rows, `[[`, character(1), "membership.status"
            ),
            source.stage = vapply(rows, `[[`, character(1), "source.stage"),
            is.primary = as.logical(vapply(
                rows, `[[`, logical(1), "is.primary"
            ))
        ),
        class = "data.frame",
        row.names = seq_along(rows)
    )
}

.basin.choose.membership <- function(candidates, basins, direction) {
    if (length(candidates) == 0L) {
        return(NA_integer_)
    }
    values <- vapply(
        candidates,
        function(index) basins[[index]]$extremum.value,
        numeric(1)
    )
    vertices <- vapply(
        candidates,
        function(index) basins[[index]]$extremum.vertex,
        integer(1)
    )
    order.index <- if (direction == "max") {
        order(-values, vertices)
    } else {
        order(values, vertices)
    }
    candidates[[order.index[[1L]]]]
}

.basin.assignment.from.membership <- function(membership,
                                               basins,
                                               n.vertices,
                                               direction,
                                               policy,
                                               backend.assignment = NULL,
                                               next.vertex = NULL) {
    if (is.null(next.vertex)) {
        next.vertex <- rep(NA_integer_, n.vertices)
    }
    rows <- vector("list", n.vertices)
    for (vertex in seq_len(n.vertices)) {
        member.rows <- which(
            membership$vertex == vertex &
                membership$direction == direction
        )
        candidate.ids <- membership$basin.id[member.rows]
        candidates <- match(
            candidate.ids,
            vapply(basins, `[[`, character(1), "basin.id")
        )
        candidates <- candidates[!is.na(candidates)]

        if (policy == "none") {
            selected <- NA_integer_
            status <- "not_applicable"
            assignment.weight <- NA_real_
        } else if (policy == "backend_primary") {
            selected <- if (length(backend.assignment) >= vertex) {
                as.integer(backend.assignment[[vertex]])
            } else {
                NA_integer_
            }
            if (is.na(selected) || selected < 1L ||
                selected > length(basins)) {
                selected <- NA_integer_
                status <- "unassigned"
                assignment.weight <- 0
            } else {
                status <- "assigned"
                assignment.weight <- 1
            }
        } else {
            selected <- .basin.choose.membership(
                candidates,
                basins,
                direction
            )
            if (is.na(selected)) {
                status <- "unassigned"
                assignment.weight <- 0
            } else {
                status <- "assigned"
                assignment.weight <- 1
            }
        }

        rows[[vertex]] <- list(
            vertex = as.integer(vertex),
            direction = direction,
            basin.id = if (is.na(selected)) {
                NA_character_
            } else {
                basins[[selected]]$basin.id
            },
            assignment.weight = assignment.weight,
            assignment.status = status,
            assignment.policy = policy,
            root.vertex = if (is.na(selected)) {
                NA_integer_
            } else {
                basins[[selected]]$extremum.vertex
            },
            next.vertex = as.integer(next.vertex[[vertex]])
        )
    }
    .basin.bind.assignment.rows(rows)
}

.basin.bind.assignment.rows <- function(rows) {
    if (length(rows) == 0L) {
        return(.empty.assignment.table())
    }
    structure(
        list(
            vertex = as.integer(vapply(rows, `[[`, integer(1), "vertex")),
            direction = vapply(rows, `[[`, character(1), "direction"),
            basin.id = vapply(rows, `[[`, character(1), "basin.id"),
            assignment.weight = as.numeric(vapply(
                rows, `[[`, numeric(1), "assignment.weight"
            )),
            assignment.status = vapply(
                rows, `[[`, character(1), "assignment.status"
            ),
            assignment.policy = vapply(
                rows, `[[`, character(1), "assignment.policy"
            ),
            root.vertex = as.integer(vapply(
                rows, `[[`, integer(1), "root.vertex"
            )),
            next.vertex = as.integer(vapply(
                rows, `[[`, integer(1), "next.vertex"
            ))
        ),
        class = "data.frame",
        row.names = seq_along(rows)
    )
}

.basin.mark.primary.membership <- function(membership, assignment) {
    if (nrow(membership) == 0L) {
        return(membership)
    }
    assignment.key <- paste(
        assignment$vertex,
        assignment$direction,
        assignment$basin.id,
        sep = "\r"
    )
    membership.key <- paste(
        membership$vertex,
        membership$direction,
        membership$basin.id,
        sep = "\r"
    )
    membership$is.primary <- membership.key %in% assignment.key[
        assignment$assignment.status == "assigned"
    ]
    membership
}

.basin.bind.extrema <- function(basins) {
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
            plateau.id = rep(NA_character_, length(basins)),
            plateau.vertices = I(lapply(
                basins,
                function(basin) as.integer(basin$extremum.vertex)
            )),
            n.plateau.vertices = rep.int(1L, length(basins)),
            is.retained = rep(TRUE, length(basins)),
            retention.status = rep("retained", length(basins))
        ),
        class = "data.frame",
        row.names = seq_along(basins)
    )
}

.basin.bind.basin.table <- function(basins,
                                    membership,
                                    assignment,
                                    method,
                                    normalized.mass) {
    if (length(basins) == 0L) {
        return(.empty.basin.table())
    }
    rows <- lapply(basins, function(basin) {
        raw.rows <- membership$basin.id == basin$basin.id
        raw.vertices <- sort(unique(membership$vertex[raw.rows]))
        primary.vertices <- sort(unique(assignment$vertex[
            assignment$basin.id == basin$basin.id &
                assignment$assignment.status == "assigned"
        ]))
        raw.allocated.mass <- if (is.null(normalized.mass)) {
            NA_real_
        } else {
            sum(
                normalized.mass[membership$vertex[raw.rows]] *
                    membership$membership.weight[raw.rows]
            )
        }
        list(
            basin.id = basin$basin.id,
            extremum.id = basin$extremum.id,
            type = basin$direction,
            method = method,
            retained = TRUE,
            extremum.vertex = basin$extremum.vertex,
            extremum.value = basin$extremum.value,
            birth.level = basin$extremum.value,
            death.level = NA_real_,
            persistence = NA_real_,
            parent.basin.id = NA_character_,
            raw.support.vertices = raw.vertices,
            raw.support.size = length(raw.vertices),
            raw.support.mass = .basin.support.mass(
                raw.vertices,
                normalized.mass
            ),
            retained.support.vertices = raw.vertices,
            retained.support.size = length(raw.vertices),
            retained.support.mass = .basin.support.mass(
                raw.vertices,
                normalized.mass
            ),
            primary.support.vertices = primary.vertices,
            primary.support.size = length(primary.vertices),
            primary.support.mass = if (length(primary.vertices) == 0L) {
                NA_real_
            } else {
                .basin.support.mass(primary.vertices, normalized.mass)
            },
            raw.allocated.mass = raw.allocated.mass,
            assignment.status = if (length(primary.vertices) == 0L) {
                "not_assigned"
            } else {
                "assigned"
            },
            retention.status = "retained"
        )
    })
    structure(
        list(
            basin.id = vapply(rows, `[[`, character(1), "basin.id"),
            extremum.id = vapply(rows, `[[`, character(1), "extremum.id"),
            type = vapply(rows, `[[`, character(1), "type"),
            method = vapply(rows, `[[`, character(1), "method"),
            retained = as.logical(vapply(
                rows, `[[`, logical(1), "retained"
            )),
            extremum.vertex = as.integer(vapply(
                rows, `[[`, integer(1), "extremum.vertex"
            )),
            extremum.value = as.numeric(vapply(
                rows, `[[`, numeric(1), "extremum.value"
            )),
            birth.level = as.numeric(vapply(
                rows, `[[`, numeric(1), "birth.level"
            )),
            death.level = as.numeric(vapply(
                rows, `[[`, numeric(1), "death.level"
            )),
            persistence = as.numeric(vapply(
                rows, `[[`, numeric(1), "persistence"
            )),
            parent.basin.id = vapply(
                rows, `[[`, character(1), "parent.basin.id"
            ),
            raw.support.vertices = I(lapply(
                rows, `[[`, "raw.support.vertices"
            )),
            raw.support.size = as.integer(vapply(
                rows, `[[`, integer(1), "raw.support.size"
            )),
            raw.support.mass = as.numeric(vapply(
                rows, `[[`, numeric(1), "raw.support.mass"
            )),
            retained.support.vertices = I(lapply(
                rows, `[[`, "retained.support.vertices"
            )),
            retained.support.size = as.integer(vapply(
                rows, `[[`, integer(1), "retained.support.size"
            )),
            retained.support.mass = as.numeric(vapply(
                rows, `[[`, numeric(1), "retained.support.mass"
            )),
            primary.support.vertices = I(lapply(
                rows, `[[`, "primary.support.vertices"
            )),
            primary.support.size = as.integer(vapply(
                rows, `[[`, integer(1), "primary.support.size"
            )),
            primary.support.mass = as.numeric(vapply(
                rows, `[[`, numeric(1), "primary.support.mass"
            )),
            raw.allocated.mass = as.numeric(vapply(
                rows, `[[`, numeric(1), "raw.allocated.mass"
            )),
            assignment.status = vapply(
                rows, `[[`, character(1), "assignment.status"
            ),
            retention.status = vapply(
                rows, `[[`, character(1), "retention.status"
            )
        ),
        class = "data.frame",
        row.names = seq_along(rows)
    )
}

.basin.legacy.trajectory.basins <- function(result, direction) {
    source <- if (direction == "max") {
        result$max.basins.all
    } else {
        result$min.basins.all
    }
    lapply(source, function(basin) {
        vertex <- as.integer(basin$extremum.vertex)
        list(
            basin.id = .basin.direction.id(direction, vertex),
            extremum.id = .basin.extremum.id(direction, vertex),
            direction = direction,
            extremum.vertex = vertex,
            extremum.value = as.numeric(basin$extremum.value),
            raw.vertices = sort(unique(as.integer(basin$vertices)))
        )
    })
}

.basin.legacy.geodesic.basins <- function(result, direction) {
    source <- if (direction == "max") result$lmax_basins else result$lmin_basins
    lapply(source, function(basin) {
        vertex <- as.integer(basin$vertex)
        list(
            basin.id = .basin.direction.id(direction, vertex),
            extremum.id = .basin.extremum.id(direction, vertex),
            direction = direction,
            extremum.vertex = vertex,
            extremum.value = as.numeric(basin$value),
            raw.vertices = sort(unique(as.integer(basin$basin_df[, 1L])))
        )
    })
}

.basin.success.object <- function(method,
                                  direction,
                                  graph.input,
                                  field,
                                  parameters,
                                  basins,
                                  membership,
                                  assignment,
                                  backend,
                                  warnings,
                                  trajectory.forest = NULL) {
    membership <- .basin.mark.primary.membership(membership, assignment)
    normalized.mass <- field$vertex.mass.normalized
    .new.basin.complex(
        method = method,
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = list(
            package.version = .basin.package.version(),
            method.backend = backend$name,
            construction = list(
                status = "complete",
                completed = TRUE,
                membership.weight.policy = "equal_share_set_membership",
                assignment.policy =
                    parameters$method.params$primary.assignment.policy
            ),
            refinement.stages = .empty.refinement.stages(),
            random.seed = .basin.or(
                parameters$method.params$tie.seed,
                NULL
            )
        ),
        status = "ok",
        extrema = .basin.bind.extrema(basins),
        basin.table = .basin.bind.basin.table(
            basins,
            membership,
            assignment,
            method,
            normalized.mass
        ),
        membership = membership,
        assignment = assignment,
        trajectory.forest = trajectory.forest,
        diagnostics = .empty.diagnostics.table(),
        warnings = warnings,
        raw = list(legacy.object = backend$object)
    )
}

.create.geodesic.reachability.complex <- function(direction,
                                                  graph.input,
                                                  field,
                                                  parameters) {
    params <- parameters$method.params
    capture <- .basin.capture.backend(
        .compute.basins.of.attraction.backend(
            graph.input$adj.list,
            graph.input$edge.length.list,
            field$construction.values,
            edge.length.quantile.thld = params$edge.length.quantile.thld,
            with.trajectories = params$store.trajectories
        )
    )
    if (!is.null(capture$error)) {
        return(.new.failed.basin.complex(
            method = "geodesic_reachability",
            direction = direction,
            graph.input = graph.input,
            field = field,
            parameters = parameters,
            condition.class = "gflow_basin_backend_error",
            message = conditionMessage(capture$error)
        ))
    }

    requested <- .basin.requested.directions(direction)
    basins.by.direction <- lapply(
        requested,
        function(current) {
            .basin.legacy.geodesic.basins(capture$value, current)
        }
    )
    names(basins.by.direction) <- requested
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
            params$primary.assignment.policy
        )
    })
    assignment <- do.call(rbind, assignments)
    row.names(assignment) <- seq_len(nrow(assignment))

    backend.record <- list(
        name = "compute.basins.of.attraction",
        object = capture$value
    )
    .basin.success.object(
        method = "geodesic_reachability",
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        basins = basins,
        membership = membership,
        assignment = assignment,
        backend = backend.record,
        warnings = capture$warnings
    )
}

.create.trajectory.flow.complex <- function(direction,
                                            graph.input,
                                            field,
                                            parameters,
                                            verbose) {
    params <- parameters$method.params
    backend.call <- function() {
        .compute.gfc.trajectory.backend(
            adj.list = graph.input$adj.list,
            weight.list = graph.input$edge.length.list,
            y = field$input.values,
            density = field$vertex.density,
            modulation = params$modulation,
            edge.length.quantile.thld = params$edge.length.quantile.thld,
            long.edge.fallback = params$long.edge.fallback,
            apply.relvalue.filter = FALSE,
            apply.maxima.clustering = FALSE,
            apply.minima.clustering = FALSE,
            apply.geometric.filter = FALSE,
            min.basin.size = 0L,
            min.n.trajectories = 0L,
            store.trajectories = params$store.trajectories,
            max.trajectory.length = params$max.trajectory.length,
            symmetric.seeding = params$symmetric.seeding,
            tie.breaking = params$tie.breaking,
            verbose = verbose
        )
    }
    capture <- if (params$tie.breaking) {
        .basin.with.local.seed(
            params$tie.seed,
            .basin.capture.backend(backend.call())
        )
    } else {
        .basin.capture.backend(backend.call())
    }
    if (!is.null(capture$error)) {
        return(.new.failed.basin.complex(
            method = "trajectory_flow",
            direction = direction,
            graph.input = graph.input,
            field = field,
            parameters = parameters,
            condition.class = "gflow_basin_backend_error",
            message = conditionMessage(capture$error)
        ))
    }

    result <- capture$value
    field$construction.values <- as.numeric(result$y)
    perturbation <- field$construction.values - field$input.values
    field$tie.policy$applied <- any(perturbation != 0)
    field$tie.policy$perturbation.scale <- if (length(perturbation) == 0L) {
        0
    } else {
        max(abs(perturbation))
    }
    if (params$modulation == "CLOSEST" &&
        params$plateau.policy == "connected_exact") {
        plateau.graph <- .build.plateau.graph(
            graph.input$adj.list,
            field$input.values,
            graph.input$validation$component
        )
        if (any(lengths(plateau.graph$vertices) > 1L)) {
            return(.create.connected.plateau.trajectory.complex(
                direction = direction,
                graph.input = graph.input,
                field = field,
                parameters = parameters,
                backend.result = result,
                warnings = capture$warnings
            ))
        }
    }

    requested <- .basin.requested.directions(direction)
    basins.by.direction <- lapply(
        requested,
        function(current) {
            .basin.legacy.trajectory.basins(result, current)
        }
    )
    names(basins.by.direction) <- requested
    basins <- unlist(basins.by.direction, recursive = FALSE)
    memberships <- lapply(requested, function(current) {
        source <- if (current == "max") {
            result$max.membership.all
        } else {
            result$min.membership.all
        }
        .basin.membership.from.index(
            source,
            basins.by.direction[[current]],
            current
        )
    })
    membership <- do.call(rbind, memberships)
    row.names(membership) <- seq_len(nrow(membership))
    assignments <- lapply(requested, function(current) {
        backend.assignment <- if (current == "max") {
            result$max.assignment
        } else {
            result$min.assignment
        }
        next.vertex <- if (current == "max") {
            result$next.up
        } else {
            result$next.down
        }
        .basin.assignment.from.membership(
            membership,
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current,
            params$primary.assignment.policy,
            backend.assignment,
            next.vertex
        )
    })
    assignment <- do.call(rbind, assignments)
    row.names(assignment) <- seq_len(nrow(assignment))

    forest <- list(
        modulation = params$modulation,
        plateau.policy = params$plateau.policy,
        long.edge.fallback = result$long.edge.fallback,
        next.vertex = lapply(requested, function(current) {
            assignment$next.vertex[assignment$direction == current]
        }),
        root.vertex = lapply(requested, function(current) {
            assignment$root.vertex[assignment$direction == current]
        }),
        trajectories = result$trajectories,
        trajectory.status = "complete",
        tie.breaking = params$tie.breaking,
        tie.seed = params$tie.seed
    )
    names(forest$next.vertex) <- requested
    names(forest$root.vertex) <- requested

    backend.record <- list(
        name = "compute.gfc.trajectory",
        object = result
    )
    .basin.success.object(
        method = "trajectory_flow",
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        basins = basins,
        membership = membership,
        assignment = assignment,
        backend = backend.record,
        warnings = capture$warnings,
        trajectory.forest = forest
    )
}

.create.rtcb.complex <- function(direction,
                                 graph.input,
                                 field,
                                 parameters) {
    params <- parameters$method.params
    capture <- .basin.capture.backend(
        compute.basins.of.attraction.rtcb(
            adj.list = graph.input$adj.list,
            weight.list = graph.input$edge.length.list,
            y = field$construction.values,
            edge.length.quantile.thld =
                params$edge.length.quantile.thld,
            with.trajectories = params$store.trajectories,
            n.min = params$n.min,
            m.min = params$m.min,
            q.min = params$q.min,
            run.max = params$run.max,
            tau0 = params$tau0,
            kappa = params$kappa,
            k.max = params$k.max,
            h.max = params$h.max,
            d.path.max = params$d.path.max,
            eta.step = params$eta.step,
            epsilon.d = params$epsilon.d,
            eps.num = params$eps.num,
            sink.prune = params$sink.prune,
            sink.prune.max.iter = params$sink.prune.max.iter,
            max.paths.per.sink = params$max.paths.per.sink
        )
    )
    if (!is.null(capture$error)) {
        return(.new.failed.basin.complex(
            method = "rtcb",
            direction = direction,
            graph.input = graph.input,
            field = field,
            parameters = parameters,
            condition.class = "gflow_basin_backend_error",
            message = conditionMessage(capture$error)
        ))
    }

    requested <- .basin.requested.directions(direction)
    basins.by.direction <- lapply(requested, function(current) {
        .basin.legacy.geodesic.basins(capture$value, current)
    })
    names(basins.by.direction) <- requested
    basins <- unlist(basins.by.direction, recursive = FALSE)
    names(basins) <- NULL
    membership <- do.call(rbind, lapply(requested, function(current) {
        .basin.membership.from.supports(
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current
        )
    }))
    row.names(membership) <- seq_len(nrow(membership))
    assignment <- do.call(rbind, lapply(requested, function(current) {
        .basin.assignment.from.membership(
            membership,
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current,
            "none"
        )
    }))
    row.names(assignment) <- seq_len(nrow(assignment))

    .new.basin.complex(
        method = "rtcb",
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = list(
            package.version = .basin.package.version(),
            method.backend = "compute.basins.of.attraction.rtcb",
            construction = list(
                status = "complete",
                completed = TRUE,
                membership.weight.policy = "equal_share_set_membership",
                assignment.policy = "none",
                requested.assignment.policy =
                    params$primary.assignment.policy,
                assignment.backend.available = FALSE,
                backend.parameters = capture$value$rtcb.params
            ),
            refinement.stages = .empty.refinement.stages(),
            random.seed = NULL
        ),
        status = "ok",
        extrema = .basin.bind.extrema(basins),
        basin.table = .basin.bind.basin.table(
            basins,
            membership,
            assignment,
            "rtcb",
            field$vertex.mass.normalized
        ),
        membership = membership,
        assignment = assignment,
        diagnostics = .empty.diagnostics.table(),
        warnings = capture$warnings,
        raw = list(legacy.object = capture$value)
    )
}

.basin.legacy.cell.basins <- function(object, direction) {
    source <- if (direction == "min") {
        object$basins$ascending
    } else {
        object$basins$descending
    }
    lapply(seq_along(source), function(index) {
        basin <- source[[index]]
        vertex <- as.integer(basin$vertex)
        list(
            basin.id = .basin.direction.id(direction, vertex),
            extremum.id = .basin.extremum.id(direction, vertex),
            direction = direction,
            extremum.vertex = vertex,
            extremum.value = as.numeric(basin$value),
            raw.vertices = sort(unique(as.integer(basin$vertices))),
            legacy.label = names(source)[[index]]
        )
    })
}

.basin.cell.summary <- function(cells) {
    rows <- list()
    for (family in names(cells)) {
        for (cell.id in names(cells[[family]])) {
            cell <- cells[[family]][[cell.id]]
            rows[[length(rows) + 1L]] <- list(
                family = family,
                cell.id = cell.id,
                n.vertices = length(cell$vertices),
                vertices = sort(unique(as.integer(cell$vertices)))
            )
        }
    }
    if (length(rows) == 0L) {
        return(structure(
            list(
                family = character(),
                cell.id = character(),
                n.vertices = integer(),
                vertices = I(list())
            ),
            class = "data.frame",
            row.names = integer()
        ))
    }
    structure(
        list(
            family = vapply(rows, `[[`, character(1), "family"),
            cell.id = vapply(rows, `[[`, character(1), "cell.id"),
            n.vertices = as.integer(vapply(
                rows, `[[`, integer(1), "n.vertices"
            )),
            vertices = I(lapply(rows, `[[`, "vertices"))
        ),
        class = "data.frame",
        row.names = seq_along(rows)
    )
}

.bind.cell.extrema <- function(object, basins) {
    extrema <- object$local_extrema
    if (nrow(extrema) == 0L) {
        return(.empty.extrema.table())
    }
    direction <- ifelse(extrema$is_maximum == 1, "max", "min")
    retained.vertices <- vapply(
        basins,
        `[[`,
        integer(1),
        "extremum.vertex"
    )
    retained <- extrema$vertex_index %in% retained.vertices
    mapping <- object$cluster_mappings
    merged.label <- mapping$merged_label[
        match(extrema$label, mapping$original_label)
    ]
    structure(
        list(
            extremum.id = mapply(
                .basin.extremum.id,
                direction,
                extrema$vertex_index,
                USE.NAMES = FALSE
            ),
            type = direction,
            representative.vertex = as.integer(extrema$vertex_index),
            extremum.value = as.numeric(extrema$fn_value),
            plateau.id = rep(NA_character_, nrow(extrema)),
            plateau.vertices = I(lapply(
                extrema$vertex_index,
                as.integer
            )),
            n.plateau.vertices = rep.int(1L, nrow(extrema)),
            is.retained = retained,
            retention.status = ifelse(
                retained,
                "retained",
                paste0("merged_into_", merged.label)
            )
        ),
        class = "data.frame",
        row.names = seq_len(nrow(extrema))
    )
}

.create.overlap.cell.complex <- function(graph.input,
                                         field,
                                         parameters) {
    params <- parameters$method.params
    capture <- .basin.capture.backend(
        .create.basin.cx.backend(
            adj.list = graph.input$adj.list,
            weight.list = graph.input$edge.length.list,
            y = field$construction.values,
            basin.merge.overlap.thld =
                params$basin.merge.overlap.thld,
            min.asc.desc.cell.size.thld =
                params$min.asc.desc.cell.size.thld,
            min.asc.asc.cell.size.thld =
                params$min.asc.asc.cell.size.thld,
            min.desc.desc.cell.size.thld =
                params$min.desc.desc.cell.size.thld,
            graph.params = params$cell.graph.params
        )
    )
    if (!is.null(capture$error)) {
        return(.new.failed.basin.complex(
            method = "overlap_cell_complex",
            direction = "both",
            graph.input = graph.input,
            field = field,
            parameters = parameters,
            condition.class = "gflow_basin_backend_error",
            message = conditionMessage(capture$error)
        ))
    }

    legacy <- capture$value
    basins.by.direction <- list(
        max = .basin.legacy.cell.basins(legacy, "max"),
        min = .basin.legacy.cell.basins(legacy, "min")
    )
    basins <- c(basins.by.direction$max, basins.by.direction$min)
    membership <- do.call(rbind, lapply(c("max", "min"), function(current) {
        .basin.membership.from.supports(
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current
        )
    }))
    row.names(membership) <- seq_len(nrow(membership))
    assignment <- do.call(rbind, lapply(c("max", "min"), function(current) {
        .basin.assignment.from.membership(
            membership,
            basins.by.direction[[current]],
            length(graph.input$adj.list),
            current,
            "none"
        )
    }))
    row.names(assignment) <- seq_len(nrow(assignment))

    initial <- legacy$initial_basin_cx
    cell.complex <- list(
        cells = legacy$cells,
        basin.intersections = list(
            initial.min.intersection =
                initial$lmin_intersection_dists,
            initial.max.intersection =
                initial$lmax_intersection_dists,
            initial.min.overlap = initial$lmin_overlap_dists,
            initial.max.overlap = initial$lmax_overlap_dists,
            merged = legacy$graph$intersection_matrix
        ),
        complex.graph = legacy$graph,
        cell.summary = .basin.cell.summary(legacy$cells),
        cluster.assignments = legacy$cluster_assignments,
        cluster.mappings = legacy$cluster_mappings,
        simplified.field = as.numeric(legacy$harmonic_predictions)
    )

    .new.basin.complex(
        method = "overlap_cell_complex",
        direction = "both",
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = list(
            package.version = .basin.package.version(),
            method.backend = "create.basin.cx",
            construction = list(
                status = "complete",
                completed = TRUE,
                membership.weight.policy = "equal_share_set_membership",
                assignment.policy = "none",
                legacy.class = class(legacy)
            ),
            refinement.stages = .empty.refinement.stages(),
            random.seed = NULL
        ),
        status = "ok",
        extrema = .bind.cell.extrema(legacy, basins),
        basin.table = .basin.bind.basin.table(
            basins,
            membership,
            assignment,
            "overlap_cell_complex",
            field$vertex.mass.normalized
        ),
        membership = membership,
        assignment = assignment,
        cell.complex = cell.complex,
        diagnostics = .empty.diagnostics.table(),
        warnings = capture$warnings,
        raw = list(
            legacy.object = legacy,
            initial.basin.complex = initial,
            merged.basins = legacy$basins
        )
    )
}

.validate.rtcb.complex <- function(object) {
    params <- object$parameters$method.params
    observed <- object$provenance$construction$backend.parameters
    expected <- list(
        n_min = params$n.min,
        m_min = params$m.min,
        q_min = if (is.null(params$q.min)) -1 else params$q.min,
        run_max = params$run.max,
        tau0 = params$tau0,
        kappa = params$kappa,
        k_max = params$k.max,
        h_max = params$h.max,
        d_path_max = params$d.path.max,
        eta_step = params$eta.step,
        epsilon_d = params$epsilon.d,
        eps_num = params$eps.num,
        sink_prune = params$sink.prune,
        sink_prune_max_iter = params$sink.prune.max.iter,
        max_paths_per_sink = params$max.paths.per.sink
    )
    if (!isTRUE(all.equal(observed, expected, check.attributes = TRUE))) {
        .stop.basin.complex(
            "Stored RTCB backend parameters do not match resolved canonical parameters.",
            "object$provenance$construction$backend.parameters",
            class = "gflow_basin_schema_error"
        )
    }
    if (any(object$assignment$assignment.status != "not_applicable") ||
        any(!is.na(object$assignment$assignment.weight))) {
        .stop.basin.complex(
            "RTCB assignments must be not-applicable while the backend supplies no primary assignment.",
            "object$assignment",
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}

.validate.overlap.cell.complex <- function(object) {
    required <- c(
        "cells",
        "basin.intersections",
        "complex.graph",
        "cell.summary",
        "cluster.assignments",
        "cluster.mappings",
        "simplified.field"
    )
    cell <- object$cell.complex
    if (!is.list(cell) || !identical(names(cell), required)) {
        .stop.basin.complex(
            "The overlap cell complex is missing required canonical fields.",
            "object$cell.complex",
            class = "gflow_basin_schema_error"
        )
    }
    families <- c(
        "asc_desc_cells",
        "asc_asc_cells",
        "desc_desc_cells"
    )
    if (!is.list(cell$cells) ||
        !all(families %in% names(cell$cells)) ||
        !is.data.frame(cell$cell.summary) ||
        !is.data.frame(cell$cluster.mappings) ||
        length(cell$simplified.field) != object$n.vertices ||
        any(!is.finite(cell$simplified.field)) ||
        is.null(cell$complex.graph)) {
        .stop.basin.complex(
            "The overlap cell-complex payload is inconsistent.",
            "object$cell.complex",
            class = "gflow_basin_schema_error"
        )
    }
    if (!inherits(object$raw$legacy.object, "basin_cx") ||
        inherits(object, "basin_cx") ||
        is.null(object$raw$initial.basin.complex) ||
        is.null(object$raw$merged.basins)) {
        .stop.basin.complex(
            "Legacy cell-complex provenance is missing or canonical classes were contaminated.",
            "object$raw",
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}
