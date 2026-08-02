# Public merge-tree objects, exact cuts, and crossing-free visualizations.

.new.basin.merge.tree <- function(object) {
    .validate.basin.complex(object)
    if (object$method != "superlevel_merge_tree") {
        .stop.basin.complex(
            "A basin merge tree requires method = 'superlevel_merge_tree'.",
            "object"
        )
    }
    tree <- list(
        schema.version = 1L,
        method = object$method,
        direction = object$direction,
        status = object$status,
        n.vertices = object$n.vertices,
        graph.input = object$graph.input,
        field = object$field,
        parameters = object$parameters,
        provenance = object$provenance,
        extrema = object$extrema,
        basin.table = object$basin.table,
        merge.table = object$merge.table,
        membership = object$membership,
        assignment = object$assignment,
        plateau.graph = object$raw$plateau.graph,
        diagnostics = object$diagnostics,
        warnings = object$warnings
    )
    class(tree) <- c("basin.merge.tree", "gflow_basin_merge_tree", "list")
    .validate.basin.merge.tree(tree)
    tree
}

.validate.basin.merge.tree.relations <- function(object,
                                                 direction,
                                                 branches,
                                                 events) {
    branch.columns <- c(
        "basin.id", "type", "extremum.vertex", "birth.level",
        "death.level", "persistence", "parent.basin.id"
    )
    event.columns <- c(
        "event.id", "direction", "losing.basin.id",
        "surviving.basin.id", "merge.level", "birth.level",
        "death.level", "persistence"
    )
    if (!all(branch.columns %in% names(branches)) ||
        !all(event.columns %in% names(events))) {
        .stop.basin.complex(
            "The merge-tree relational fields are incomplete.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    branch.values <- branches[
        ,
        c("birth.level", "death.level", "persistence"),
        drop = FALSE
    ]
    if (any(!is.finite(as.matrix(branch.values))) ||
        any(branches$persistence < 0)) {
        .stop.basin.complex(
            paste(
                "Merge-tree branch birth, death, and persistence values",
                "must be finite with nonnegative persistence."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    event.values <- events[
        ,
        c("merge.level", "birth.level", "death.level", "persistence"),
        drop = FALSE
    ]
    if (nrow(events) > 0L &&
        (any(!is.finite(as.matrix(event.values))) ||
            any(events$persistence < 0))) {
        .stop.basin.complex(
            paste(
                "Merge-tree event levels and persistence values must be",
                "finite with nonnegative persistence."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    expected.persistence <- if (direction == "max") {
        branches$birth.level - branches$death.level
    } else {
        branches$death.level - branches$birth.level
    }
    if (any(branches$persistence != expected.persistence)) {
        .stop.basin.complex(
            paste(
                "Merge-tree branch persistence does not equal the exact",
                "direction-specific birth-to-death difference."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    graph.component <- object$graph.input$validation$component
    field <- object$field$construction.values
    extremum.vertex <- branches$extremum.vertex
    if (length(graph.component) != object$n.vertices ||
        length(field) != object$n.vertices ||
        any(!is.finite(field)) ||
        !is.numeric(extremum.vertex) ||
        anyNA(extremum.vertex) ||
        any(!is.finite(extremum.vertex)) ||
        any(extremum.vertex != floor(extremum.vertex)) ||
        any(extremum.vertex < 1L | extremum.vertex > object$n.vertices)) {
        .stop.basin.complex(
            paste(
                "The merge-tree field, graph components, and extremum",
                "vertices are inconsistent."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    branch.component <- as.integer(graph.component[extremum.vertex])
    if (anyNA(branch.component)) {
        .stop.basin.complex(
            "The merge-tree branch component mapping is incomplete.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    roots <- is.na(branches$parent.basin.id)
    for (component in sort(unique(branch.component))) {
        component.rows <- which(branch.component == component)
        component.roots <- component.rows[roots[component.rows]]
        if (length(component.roots) != 1L) {
            .stop.basin.complex(
                paste(
                    "Each merge-tree graph component must contain exactly",
                    "one canonical root."
                ),
                "object",
                class = "gflow_basin_schema_error"
            )
        }
        component.vertices <- which(graph.component == component)
        expected.root.death <- if (direction == "max") {
            min(field[component.vertices])
        } else {
            max(field[component.vertices])
        }
        if (branches$death.level[[component.roots]] != expected.root.death) {
            .stop.basin.complex(
                paste(
                    "The merge-tree component-root death level does not",
                    "equal the finite component field boundary."
                ),
                "object",
                class = "gflow_basin_schema_error"
            )
        }
    }

    nonroot.rows <- which(!roots)
    if (length(nonroot.rows) == 0L) {
        return(invisible(TRUE))
    }
    event.index <- match(
        branches$basin.id[nonroot.rows],
        events$losing.basin.id
    )
    if (anyNA(event.index) ||
        anyDuplicated(events$losing.basin.id) ||
        nrow(events) != length(nonroot.rows)) {
        .stop.basin.complex(
            "Each nonroot merge-tree branch must match exactly one event.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    parent.index <- match(
        branches$parent.basin.id[nonroot.rows],
        branches$basin.id
    )
    if (anyNA(parent.index) ||
        any(branch.component[parent.index] !=
            branch.component[nonroot.rows])) {
        .stop.basin.complex(
            paste(
                "Each nonroot merge-tree parent must be a branch in the",
                "same direction and graph component."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    matched.events <- events[event.index, , drop = FALSE]
    matched.branches <- branches[nonroot.rows, , drop = FALSE]
    if (any(matched.events$direction != direction) ||
        any(matched.events$surviving.basin.id !=
            matched.branches$parent.basin.id)) {
        .stop.basin.complex(
            paste(
                "A merge-tree event survivor disagrees with its losing",
                "branch parent."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    relation.mismatch <-
        matched.events$birth.level != matched.branches$birth.level |
        matched.events$death.level != matched.branches$death.level |
        matched.events$merge.level != matched.branches$death.level |
        matched.events$persistence != matched.branches$persistence
    if (any(relation.mismatch)) {
        .stop.basin.complex(
            paste(
                "A merge-tree event's birth, death, merge, or persistence",
                "value disagrees with its losing branch."
            ),
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}

.validate.basin.merge.tree <- function(object) {
    required <- c(
        "schema.version", "method", "direction", "status", "n.vertices",
        "graph.input", "field", "parameters", "provenance", "extrema",
        "basin.table", "merge.table", "membership", "assignment",
        "plateau.graph", "diagnostics", "warnings"
    )
    expected.class <- c(
        "basin.merge.tree", "gflow_basin_merge_tree", "list"
    )
    if (!is.list(object) || !identical(class(object), expected.class) ||
        !all(required %in% names(object))) {
        .stop.basin.complex(
            "'object' is not a complete canonical basin merge tree.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    if (!identical(object$schema.version, 1L) ||
        object$method != "superlevel_merge_tree" ||
        !object$direction %in% .basin.complex.directions ||
        !object$status %in% c("ok", "partial", "failed")) {
        .stop.basin.complex(
            "The basin merge-tree identity fields are invalid.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    if (length(object$graph.input$adj.list) != object$n.vertices ||
        length(object$field$construction.values) != object$n.vertices) {
        .stop.basin.complex(
            "The basin merge-tree graph and field sizes are inconsistent.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    if (!is.data.frame(object$basin.table) ||
        !is.data.frame(object$merge.table) ||
        !is.data.frame(object$membership) ||
        !is.data.frame(object$assignment)) {
        .stop.basin.complex(
            "The basin merge-tree canonical tables are invalid.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    if (object$status == "ok") {
        directions <- .basin.requested.directions(object$direction)
        if (any(!object$basin.table$type %in% directions) ||
            any(!object$merge.table$direction %in% directions)) {
            .stop.basin.complex(
                "The merge-tree tables contain an unrequested direction.",
                "object",
                class = "gflow_basin_schema_error"
            )
        }
        for (direction in directions) {
            branches <- object$basin.table[
                object$basin.table$type == direction, , drop = FALSE
            ]
            events <- object$merge.table[
                object$merge.table$direction == direction, , drop = FALSE
            ]
            roots <- is.na(branches$parent.basin.id)
            if (anyDuplicated(branches$basin.id) ||
                anyDuplicated(events$losing.basin.id) ||
                !setequal(
                    branches$basin.id[!roots],
                    events$losing.basin.id
                )) {
                .stop.basin.complex(
                    "The merge-tree branch hierarchy is inconsistent.",
                    "object",
                    class = "gflow_basin_schema_error"
                )
            }
            .validate.basin.merge.tree.relations(
                object, direction, branches, events
            )
        }
    }
    invisible(TRUE)
}

.as.basin.merge.tree <- function(object) {
    if (inherits(object, "basin.merge.tree")) {
        .validate.basin.merge.tree(object)
        return(object)
    }
    if (inherits(object, "basin_complex")) {
        return(get.basin.merge.tree(object, required = TRUE))
    }
    .stop.basin.complex(
        "'object' must be a basin.merge.tree or basin_complex.",
        "object"
    )
}

.basin.merge.tree.assert.ready <- function(object) {
    .validate.basin.merge.tree(object)
    if (object$status != "ok") {
        .stop.basin.complex(
            "Merge-tree operations require status = 'ok'.",
            "object"
        )
    }
    invisible(TRUE)
}

.basin.merge.tree.direction <- function(object, direction) {
    direction <- .basin.assert.choice(
        direction, c("max", "min"), "direction"
    )
    requested <- .basin.requested.directions(object$direction)
    if (!direction %in% requested) {
        .stop.basin.complex(
            sprintf(
                "Direction '%s' was not constructed in this merge tree.",
                direction
            ),
            "direction"
        )
    }
    direction
}

.basin.merge.tree.branch.table <- function(object, direction) {
    branches <- object$basin.table[
        object$basin.table$type == direction, , drop = FALSE
    ]
    branches$graph.component <- as.integer(
        object$graph.input$validation$component[branches$extremum.vertex]
    )
    branches
}

.basin.merge.tree.component <- function(branches, component) {
    available <- sort(unique(branches$graph.component))
    if (is.null(component)) {
        if (length(available) != 1L) {
            .stop.basin.complex(
                paste(
                    "The selected direction is a merge forest.",
                    "Supply one graph component for dendrogram conversion",
                    "or plotting."
                ),
                "component",
                details = list(available.components = available)
            )
        }
        return(available[[1L]])
    }
    component <- .basin.assert.integer(
        component, "component", lower = 1
    )
    if (!component %in% available) {
        .stop.basin.complex(
            sprintf(
                "Graph component %d is not present; available components: %s.",
                component,
                paste(available, collapse = ", ")
            ),
            "component"
        )
    }
    component
}

.basin.merge.tree.labels <- function(branches,
                                     label,
                                     labels = NULL) {
    choices <- c("basin.id", "extremum.id", "extremum.vertex")
    label <- .basin.assert.choice(label, choices, "label")
    if (is.null(labels)) {
        out <- branches[[label]]
        if (label == "extremum.vertex") {
            out <- paste0("v", out)
        }
    } else {
        if (!is.atomic(labels) || is.object(labels) || anyNA(labels)) {
            .stop.basin.complex(
                "'labels' must be a nonmissing atomic vector.",
                "labels"
            )
        }
        if (!is.null(names(labels)) && all(nzchar(names(labels)))) {
            index <- match(branches$basin.id, names(labels))
            if (anyNA(index)) {
                .stop.basin.complex(
                    "Named 'labels' must cover every displayed basin id.",
                    "labels"
                )
            }
            out <- labels[index]
        } else {
            if (length(labels) != nrow(branches)) {
                .stop.basin.complex(
                    "Unnamed 'labels' must have one value per displayed basin.",
                    "labels"
                )
            }
            out <- labels
        }
    }
    out <- as.character(out)
    if (anyNA(out) || any(!nzchar(out)) || anyDuplicated(out)) {
        .stop.basin.complex(
            "Resolved merge-tree labels must be unique and nonempty.",
            "labels"
        )
    }
    out
}

.basin.merge.tree.leaf.order <- function(merge, reference) {
    if (reference < 0L) {
        return(-reference)
    }
    c(
        .basin.merge.tree.leaf.order(merge, merge[reference, 1L]),
        .basin.merge.tree.leaf.order(merge, merge[reference, 2L])
    )
}

.basin.merge.tree.branch.depth <- function(branches) {
    parent <- setNames(branches$parent.basin.id, branches$basin.id)
    vapply(branches$basin.id, function(basin.id) {
        depth <- 0L
        current <- basin.id
        visited <- character()
        repeat {
            next.id <- parent[[current]]
            if (is.na(next.id)) {
                break
            }
            if (next.id %in% visited || !next.id %in% names(parent)) {
                .stop.basin.complex(
                    "The merge-tree hierarchy contains a cycle or missing parent.",
                    "object",
                    class = "gflow_basin_schema_error"
                )
            }
            visited <- c(visited, current)
            current <- next.id
            depth <- depth + 1L
        }
        depth
    }, integer(1))
}

.basin.merge.tree.requested.ids <- function(basin.ids) {
    if (!is.character(basin.ids) ||
        length(basin.ids) < 1L ||
        anyNA(basin.ids) ||
        any(!nzchar(basin.ids)) ||
        anyDuplicated(basin.ids)) {
        .stop.basin.complex(
            paste(
                "'basin.ids' must be a nonempty character vector of unique,",
                "nonmissing canonical basin ids."
            ),
            "basin.ids"
        )
    }
    basin.ids
}

.basin.merge.tree.ancestor.closure <- function(branches, basin.ids) {
    parent <- setNames(branches$parent.basin.id, branches$basin.id)
    closure <- basin.ids
    for (basin.id in basin.ids) {
        current <- basin.id
        visited <- character()
        repeat {
            parent.id <- parent[[current]]
            if (is.null(parent.id)) {
                .stop.basin.complex(
                    "The merge-tree hierarchy contains a missing parent.",
                    "object",
                    class = "gflow_basin_schema_error"
                )
            }
            if (is.na(parent.id)) {
                break
            }
            if (parent.id %in% visited) {
                .stop.basin.complex(
                    "The merge-tree hierarchy contains a cycle.",
                    "object",
                    class = "gflow_basin_schema_error"
                )
            }
            visited <- c(visited, current)
            closure <- c(closure, parent.id)
            current <- parent.id
        }
    }
    branches$basin.id[branches$basin.id %in% unique(closure)]
}

.basin.merge.tree.layout.selected <- function(tree,
                                              direction,
                                              component,
                                              branches,
                                              label,
                                              labels = NULL) {
    resolved.labels <- .basin.merge.tree.labels(
        branches, label, labels
    )
    n.leaf <- nrow(branches)
    leaf.reference <- setNames(-seq_len(n.leaf), branches$basin.id)
    roots <- which(is.na(branches$parent.basin.id))
    if (length(roots) != 1L) {
        .stop.basin.complex(
            "A plotted graph component must contain exactly one tree root.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    events <- tree$merge.table[
        tree$merge.table$direction == direction &
            tree$merge.table$losing.basin.id %in% branches$basin.id,
        ,
        drop = FALSE
    ]

    if (n.leaf == 1L) {
        return(list(
            tree = tree,
            direction = direction,
            component = component,
            branches = branches,
            labels = resolved.labels,
            events = events,
            merge = matrix(integer(), nrow = 0L, ncol = 2L),
            merge.level = numeric(),
            transformed.height = numeric(),
            order = 1L,
            hclust = NULL,
            root.leaf = roots[[1L]]
        ))
    }

    child.index <- match(events$losing.basin.id, branches$basin.id)
    transformed <- if (direction == "max") {
        max(branches$birth.level) - events$merge.level
    } else {
        events$merge.level - min(branches$birth.level)
    }
    branch.depth <- .basin.merge.tree.branch.depth(branches)
    event.order <- order(
        transformed,
        -branch.depth[child.index],
        branches$extremum.vertex[child.index],
        events$losing.basin.id
    )
    events <- events[event.order, , drop = FALSE]
    transformed <- transformed[event.order]
    if (nrow(events) != n.leaf - 1L ||
        any(diff(transformed) < -sqrt(.Machine$double.eps))) {
        .stop.basin.complex(
            "Merge events cannot be represented as a monotone dendrogram.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }

    merge <- matrix(NA_integer_, nrow = nrow(events), ncol = 2L)
    cluster.reference <- leaf.reference
    for (index in seq_len(nrow(events))) {
        child <- events$losing.basin.id[[index]]
        parent <- events$surviving.basin.id[[index]]
        if (!child %in% names(cluster.reference) ||
            !parent %in% names(cluster.reference)) {
            .stop.basin.complex(
                "A merge event references a basin outside its component.",
                "object",
                class = "gflow_basin_schema_error"
            )
        }
        merge[index, ] <- c(
            cluster.reference[[parent]],
            cluster.reference[[child]]
        )
        cluster.reference[[parent]] <- index
    }
    root.id <- branches$basin.id[[roots]]
    if (cluster.reference[[root.id]] != nrow(merge)) {
        .stop.basin.complex(
            "The merge-event hierarchy does not terminate at its root.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    leaf.order <- .basin.merge.tree.leaf.order(merge, nrow(merge))
    if (!identical(sort(leaf.order), seq_len(n.leaf))) {
        .stop.basin.complex(
            "The dendrogram leaf order is not a permutation.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    hclust <- structure(
        list(
            merge = merge,
            height = as.numeric(transformed),
            order = as.integer(leaf.order),
            labels = resolved.labels,
            method = "exact graph level-set filtration",
            call = match.call(),
            dist.method = "not applicable"
        ),
        class = "hclust"
    )
    list(
        tree = tree,
        direction = direction,
        component = component,
        branches = branches,
        labels = resolved.labels,
        events = events,
        merge = merge,
        merge.level = as.numeric(events$merge.level),
        transformed.height = as.numeric(transformed),
        order = as.integer(leaf.order),
        hclust = hclust,
        root.leaf = roots[[1L]]
    )
}

.basin.merge.tree.layout.coordinates <- function(layout) {
    coordinates <- .basin.merge.tree.plot.coordinates(layout)
    branch.coordinates <- data.frame(
        basin.id = layout$branches$basin.id,
        x = as.numeric(coordinates$leaf.x),
        birth.level = as.numeric(layout$branches$birth.level),
        death.level = as.numeric(layout$branches$death.level),
        stringsAsFactors = FALSE
    )
    if (nrow(layout$events) == 0L) {
        event.coordinates <- data.frame(
            event.id = character(),
            losing.basin.id = character(),
            surviving.basin.id = character(),
            losing.x = numeric(),
            surviving.x = numeric(),
            merge.level = numeric(),
            stringsAsFactors = FALSE
        )
    } else {
        child.index <- match(
            layout$events$losing.basin.id,
            layout$branches$basin.id
        )
        parent.index <- match(
            layout$events$surviving.basin.id,
            layout$branches$basin.id
        )
        event.coordinates <- data.frame(
            event.id = layout$events$event.id,
            losing.basin.id = layout$events$losing.basin.id,
            surviving.basin.id = layout$events$surviving.basin.id,
            losing.x = as.numeric(coordinates$leaf.x[child.index]),
            surviving.x = as.numeric(coordinates$leaf.x[parent.index]),
            merge.level = as.numeric(layout$events$merge.level),
            stringsAsFactors = FALSE
        )
    }
    coordinates$branches <- branch.coordinates
    coordinates$events <- event.coordinates
    coordinates
}

#' Extract a Validated Basin Merge-Tree Layout
#'
#' Returns the canonical branches, merge events, crossing-free leaf order, and
#' plotting coordinates for one direction and graph component without
#' drawing. A restricted set of canonical basin ids must contain its component
#' root and every ancestor unless `close.ancestors = TRUE`, in which case the
#' missing ancestors are added and reported. Scalar-field birth, death, merge,
#' and persistence values are preserved exactly from the canonical merge tree.
#'
#' This accessor is the common computational boundary for static and
#' interactive merge-tree consumers. Restricting the display never constructs
#' an induced graph or recomputes the canonical basin complex.
#'
#' @param x A `basin.merge.tree`, or a compatible `basin_complex`.
#' @param direction Tree orientation, `"max"` or `"min"`.
#' @param component Graph-component number. It may be omitted when the selected
#'   direction has one component, or when `basin.ids` identify exactly one
#'   component.
#' @param basin.ids Optional nonempty character vector of canonical basin ids.
#'   `NULL` selects the complete component.
#' @param close.ancestors Add missing canonical ancestors of `basin.ids`.
#' @param label Canonical basin-table field used for branch labels.
#' @param labels Optional custom labels. A named vector is matched by basin id;
#'   an unnamed vector follows deterministic selected-branch order.
#'
#' @return A `basin.merge.tree.layout` list containing the requested and
#'   closure-added ids, selected canonical branch and event tables, component
#'   root, restricted crossing-free leaf order, dendrogram representation,
#'   branch and event coordinates, and validation status.
#' @export
get.basin.merge.tree.layout <- function(
    x,
    direction = "max",
    component = NULL,
    basin.ids = NULL,
    close.ancestors = FALSE,
    label = c("basin.id", "extremum.id", "extremum.vertex"),
    labels = NULL
) {
    tree <- .as.basin.merge.tree(x)
    .basin.merge.tree.assert.ready(tree)
    direction <- .basin.merge.tree.direction(tree, direction)
    .basin.assert.logical(close.ancestors, "close.ancestors")
    label <- match.arg(label)

    direction.branches <- .basin.merge.tree.branch.table(tree, direction)
    if (is.null(basin.ids)) {
        component <- .basin.merge.tree.component(
            direction.branches, component
        )
    } else {
        basin.ids <- .basin.merge.tree.requested.ids(basin.ids)
        all.ids <- tree$basin.table$basin.id
        unknown <- setdiff(basin.ids, all.ids)
        if (length(unknown) > 0L) {
            .stop.basin.complex(
                sprintf(
                    "Unknown canonical basin id%s: %s.",
                    if (length(unknown) == 1L) "" else "s",
                    paste(unknown, collapse = ", ")
                ),
                "basin.ids",
                details = list(unknown.ids = unknown)
            )
        }
        wrong.direction <- setdiff(
            basin.ids, direction.branches$basin.id
        )
        if (length(wrong.direction) > 0L) {
            .stop.basin.complex(
                "The requested basin ids contain a different direction.",
                "basin.ids",
                details = list(direction.mismatch.ids = wrong.direction)
            )
        }
        requested.index <- match(
            basin.ids, direction.branches$basin.id
        )
        requested.components <- unique(
            direction.branches$graph.component[requested.index]
        )
        if (length(requested.components) != 1L) {
            .stop.basin.complex(
                "The requested basin ids span multiple graph components.",
                "basin.ids",
                details = list(
                    requested.components = sort(requested.components)
                )
            )
        }
        if (is.null(component)) {
            component <- requested.components[[1L]]
        } else {
            component <- .basin.merge.tree.component(
                direction.branches, component
            )
            if (!identical(component, requested.components[[1L]])) {
                .stop.basin.complex(
                    paste(
                        "The requested basin ids do not belong to the",
                        "selected graph component."
                    ),
                    "basin.ids",
                    details = list(
                        component = component,
                        requested.components = requested.components
                    )
                )
            }
        }
    }

    component.branches <- direction.branches[
        direction.branches$graph.component == component,
        ,
        drop = FALSE
    ]
    component.branches <- component.branches[
        order(
            component.branches$extremum.vertex,
            component.branches$basin.id
        ),
        ,
        drop = FALSE
    ]
    row.names(component.branches) <- seq_len(nrow(component.branches))

    requested.ids <- if (is.null(basin.ids)) {
        component.branches$basin.id
    } else {
        basin.ids
    }
    closed.ids <- .basin.merge.tree.ancestor.closure(
        component.branches, requested.ids
    )
    closure.added.ids <- setdiff(closed.ids, requested.ids)
    if (!close.ancestors && length(closure.added.ids) > 0L) {
        .stop.basin.complex(
            paste(
                "The requested basin ids are not ancestor-closed and do not",
                "contain the component root."
            ),
            "basin.ids",
            details = list(missing.ancestor.ids = closure.added.ids)
        )
    }
    selected.ids <- if (close.ancestors) closed.ids else {
        component.branches$basin.id[
            component.branches$basin.id %in% requested.ids
        ]
    }
    branches <- component.branches[
        component.branches$basin.id %in% selected.ids,
        ,
        drop = FALSE
    ]
    row.names(branches) <- seq_len(nrow(branches))
    layout <- .basin.merge.tree.layout.selected(
        tree, direction, component, branches, label, labels
    )

    if (length(selected.ids) < nrow(component.branches)) {
        complete.layout <- .basin.merge.tree.layout.selected(
            tree,
            direction,
            component,
            component.branches,
            "basin.id"
        )
        complete.order <- component.branches$basin.id[
            complete.layout$order
        ]
        expected.order <- complete.order[
            complete.order %in% layout$branches$basin.id
        ]
        actual.order <- layout$branches$basin.id[layout$order]
        if (!identical(actual.order, expected.order)) {
            .stop.basin.complex(
                paste(
                    "The restricted merge-tree leaf order does not preserve",
                    "the canonical complete-tree order."
                ),
                "object",
                class = "gflow_basin_schema_error"
            )
        }
    }

    layout$requested.basin.ids <- requested.ids
    layout$closure.added.ids <- component.branches$basin.id[
        component.branches$basin.id %in% closure.added.ids
    ]
    layout$basin.ids <- layout$branches$basin.id
    layout$component.root.basin.id <-
        layout$branches$basin.id[[layout$root.leaf]]
    layout$leaf.order <- layout$branches$basin.id[layout$order]
    layout$coordinates <- .basin.merge.tree.layout.coordinates(layout)
    layout$validation.status <- "ok"
    class(layout) <- c("basin.merge.tree.layout", "list")
    layout
}

.basin.merge.tree.layout <- function(object,
                                     direction,
                                     component,
                                     label,
                                     labels = NULL) {
    get.basin.merge.tree.layout(
        object,
        direction = direction,
        component = component,
        label = label,
        labels = labels
    )
}

.basin.merge.tree.single.dendrogram <- function(label) {
    structure(
        1L,
        members = 1L,
        height = 0,
        label = label,
        leaf = TRUE,
        class = "dendrogram"
    )
}

#' Convert a Basin Merge Tree to a Dendrogram
#'
#' Converts one connected component and one orientation of a canonical basin
#' merge tree to a crossing-free `dendrogram`. The conversion preserves the
#' exact branch hierarchy and deterministic leaf order. Its nonnegative
#' `hclust` heights are display coordinates measured from the most extreme
#' birth level; exact scalar-field merge levels remain attached as attributes.
#' Consequently, the returned dendrogram must not be interpreted as a tree
#' inferred from pairwise distances.
#'
#' @param object A `basin.merge.tree`, or a compatible `basin_complex`.
#' @param direction Tree orientation, `"max"` or `"min"`.
#' @param component Graph-component number. It may be omitted when the selected
#'   direction has exactly one graph component.
#' @param label Canonical basin-table field used for leaf labels.
#' @param labels Optional custom labels. A named vector is matched by basin id;
#'   an unnamed vector follows deterministic basin-table order.
#' @param ... Reserved for S3 compatibility.
#'
#' @return A `dendrogram`. Exact branch metadata are available in the
#'   `basin.merge.tree.branches`, `basin.merge.tree.events`, `direction`, and
#'   `graph.component` attributes.
#' @method as.dendrogram basin.merge.tree
#' @export
#' @rawNamespace export(as.dendrogram.basin.merge.tree)
as.dendrogram.basin.merge.tree <- function(
    object,
    direction = "max",
    component = NULL,
    label = c("basin.id", "extremum.id", "extremum.vertex"),
    labels = NULL,
    ...
) {
    label <- match.arg(label)
    layout <- .basin.merge.tree.layout(
        object, direction, component, label, labels
    )
    dendrogram <- if (nrow(layout$merge) == 0L) {
        .basin.merge.tree.single.dendrogram(layout$labels[[1L]])
    } else {
        stats::as.dendrogram(layout$hclust)
    }
    attr(dendrogram, "basin.merge.tree.branches") <- layout$branches
    events <- layout$events
    events$transformed.height <- layout$transformed.height
    attr(dendrogram, "basin.merge.tree.events") <- events
    attr(dendrogram, "direction") <- layout$direction
    attr(dendrogram, "graph.component") <- layout$component
    dendrogram
}

.basin.merge.tree.induced.components <- function(adj.list,
                                                 active,
                                                 graph.component = NULL) {
    n <- length(adj.list)
    component <- integer(n)
    active.vertices <- which(active)
    count <- 0L
    for (start in active.vertices) {
        if (component[[start]] != 0L) {
            next
        }
        count <- count + 1L
        queue <- start
        component[[start]] <- count
        cursor <- 1L
        while (cursor <= length(queue)) {
            vertex <- queue[[cursor]]
            cursor <- cursor + 1L
            neighbors <- adj.list[[vertex]]
            neighbors <- neighbors[
                active[neighbors] & component[neighbors] == 0L
            ]
            if (!is.null(graph.component)) {
                neighbors <- neighbors[
                    graph.component[neighbors] ==
                        graph.component[[start]]
                ]
            }
            if (length(neighbors) > 0L) {
                neighbors <- sort(unique(neighbors))
                component[neighbors] <- count
                queue <- c(queue, neighbors)
            }
        }
    }
    component
}

#' Cut a Basin Merge Tree at a Scalar-Field Height
#'
#' Computes connected components of the exact graph superlevel set
#' `field >= height` for maximum trees, or sublevel set `field <= height` for
#' minimum trees. Equality is included, so a saddle plateau has already joined
#' its incident components at its exact scalar-field level. Every returned
#' component is labeled by the elder basin that survives at the cut.
#'
#' @param x A `basin.merge.tree`, or a compatible `basin_complex`.
#' @param height Finite scalar-field cut height.
#' @param direction Tree orientation, `"max"` or `"min"`.
#' @param component Optional graph-component number. When omitted, all graph
#'   components are cut.
#' @param ... Reserved for S3 compatibility.
#'
#' @return A `basin.merge.tree.cut` list with component and vertex-membership
#'   tables, the number of active vertices and components, and cut metadata.
#' @method cut basin.merge.tree
#' @export
#' @rawNamespace export(cut.basin.merge.tree)
cut.basin.merge.tree <- function(
    x,
    height,
    direction = "max",
    component = NULL,
    ...
) {
    tree <- .as.basin.merge.tree(x)
    .basin.merge.tree.assert.ready(tree)
    direction <- .basin.merge.tree.direction(tree, direction)
    height <- .basin.assert.number(height, "height")
    graph.component <- as.integer(
        tree$graph.input$validation$component
    )
    available <- sort(unique(graph.component))
    if (!is.null(component)) {
        component <- .basin.assert.integer(
            component, "component", lower = 1
        )
        if (!component %in% available) {
            .stop.basin.complex(
                sprintf("Graph component %d is not present.", component),
                "component"
            )
        }
    }

    field <- tree$field$construction.values
    active <- if (direction == "max") {
        field >= height
    } else {
        field <= height
    }
    if (!is.null(component)) {
        active <- active & graph.component == component
    }
    component.index <- .basin.merge.tree.induced.components(
        tree$graph.input$adj.list,
        active,
        graph.component
    )
    active.components <- sort(unique(component.index[component.index > 0L]))
    branches <- .basin.merge.tree.branch.table(tree, direction)
    mass <- tree$field$vertex.mass.normalized
    vertex.id <- tree$graph.input$vertex.id

    rows <- lapply(seq_along(active.components), function(index) {
        vertices <- which(component.index == active.components[[index]])
        candidates <- branches[
            branches$extremum.vertex %in% vertices, , drop = FALSE
        ]
        if (nrow(candidates) == 0L) {
            .stop.basin.complex(
                "A cut component contains no active extremum.",
                "object",
                class = "gflow_basin_schema_error"
            )
        }
        elder.order <- if (direction == "max") {
            order(-candidates$birth.level, candidates$extremum.vertex)
        } else {
            order(candidates$birth.level, candidates$extremum.vertex)
        }
        elder <- candidates[elder.order[[1L]], , drop = FALSE]
        list(
            component.id = sprintf(
                "component_%s_%04d", direction, index
            ),
            graph.component = graph.component[[vertices[[1L]]]],
            basin.id = elder$basin.id[[1L]],
            extremum.id = elder$extremum.id[[1L]],
            extremum.vertex = elder$extremum.vertex[[1L]],
            extremum.value = elder$extremum.value[[1L]],
            support.size = length(vertices),
            support.mass = if (is.null(mass)) {
                NA_real_
            } else {
                sum(mass[vertices])
            },
            vertices = as.integer(vertices),
            vertex.ids = vertex.id[vertices]
        )
    })

    components <- if (length(rows) == 0L) {
        structure(
            list(
                component.id = character(),
                graph.component = integer(),
                basin.id = character(),
                extremum.id = character(),
                extremum.vertex = integer(),
                extremum.value = numeric(),
                support.size = integer(),
                support.mass = numeric(),
                vertices = I(list()),
                vertex.ids = I(list())
            ),
            class = "data.frame",
            row.names = integer()
        )
    } else {
        structure(
            list(
                component.id = vapply(
                    rows, `[[`, character(1), "component.id"
                ),
                graph.component = as.integer(vapply(
                    rows, `[[`, integer(1), "graph.component"
                )),
                basin.id = vapply(
                    rows, `[[`, character(1), "basin.id"
                ),
                extremum.id = vapply(
                    rows, `[[`, character(1), "extremum.id"
                ),
                extremum.vertex = as.integer(vapply(
                    rows, `[[`, integer(1), "extremum.vertex"
                )),
                extremum.value = as.numeric(vapply(
                    rows, `[[`, numeric(1), "extremum.value"
                )),
                support.size = as.integer(vapply(
                    rows, `[[`, integer(1), "support.size"
                )),
                support.mass = as.numeric(vapply(
                    rows, `[[`, numeric(1), "support.mass"
                )),
                vertices = I(lapply(rows, `[[`, "vertices")),
                vertex.ids = I(lapply(rows, `[[`, "vertex.ids"))
            ),
            class = "data.frame",
            row.names = seq_along(rows)
        )
    }
    membership <- if (nrow(components) == 0L) {
        data.frame(
            vertex = integer(),
            vertex.id = character(),
            component.id = character(),
            basin.id = character(),
            stringsAsFactors = FALSE
        )
    } else {
        do.call(rbind, lapply(seq_len(nrow(components)), function(index) {
            vertices <- components$vertices[[index]]
            data.frame(
                vertex = vertices,
                vertex.id = vertex.id[vertices],
                component.id = components$component.id[[index]],
                basin.id = components$basin.id[[index]],
                stringsAsFactors = FALSE
            )
        }))
    }
    row.names(membership) <- seq_len(nrow(membership))
    result <- list(
        direction = direction,
        height = height,
        relation = if (direction == "max") ">=" else "<=",
        graph.component = component,
        n.active.vertices = sum(active),
        n.components = nrow(components),
        components = components,
        membership = membership
    )
    class(result) <- c("basin.merge.tree.cut", "list")
    result
}

.basin.merge.tree.colors <- function(branches, branch.col) {
    palette <- c(
        "#A4243B", "#2F6690", "#3D8B6D", "#D9A441",
        "#E56B2F", "#7A5195", "#4C78A8", "#59A14F",
        "#B07AA1", "#F28E2B", "#76B7B2", "#E15759"
    )
    if (is.null(branch.col)) {
        return(rep(palette, length.out = nrow(branches)))
    }
    if (!is.character(branch.col) || anyNA(branch.col)) {
        .stop.basin.complex(
            "'branch.col' must contain nonmissing colors.",
            "branch.col"
        )
    }
    if (!is.null(names(branch.col)) && all(nzchar(names(branch.col)))) {
        index <- match(branches$basin.id, names(branch.col))
        if (anyNA(index)) {
            .stop.basin.complex(
                "Named 'branch.col' must cover every displayed basin id.",
                "branch.col"
            )
        }
        return(branch.col[index])
    }
    if (length(branch.col) == 1L) {
        return(rep(branch.col, nrow(branches)))
    }
    if (length(branch.col) != nrow(branches)) {
        .stop.basin.complex(
            "Unnamed 'branch.col' must have one color per displayed basin.",
            "branch.col"
        )
    }
    branch.col
}

.basin.merge.tree.format.annotation <- function(x, digits = 3L) {
    ifelse(
        is.na(x),
        "NA",
        formatC(x, digits = digits, format = "fg", flag = "#")
    )
}

.basin.merge.tree.plot.coordinates <- function(layout) {
    n.leaf <- nrow(layout$branches)
    leaf.x <- numeric(n.leaf)
    leaf.x[layout$order] <- seq_len(n.leaf)
    if (nrow(layout$merge) == 0L) {
        return(list(
            leaf.x = leaf.x,
            node.x = numeric(),
            node.leaves = list()
        ))
    }
    node.x <- numeric(nrow(layout$merge))
    node.leaves <- vector("list", nrow(layout$merge))
    reference.leaves <- function(reference) {
        if (reference < 0L) {
            return(-reference)
        }
        node.leaves[[reference]]
    }
    for (index in seq_len(nrow(layout$merge))) {
        node.leaves[[index]] <- c(
            reference.leaves(layout$merge[index, 1L]),
            reference.leaves(layout$merge[index, 2L])
        )
        parent.index <- match(
            layout$events$surviving.basin.id[[index]],
            layout$branches$basin.id
        )
        node.x[[index]] <- leaf.x[[parent.index]]
    }
    list(
        leaf.x = leaf.x,
        node.x = node.x,
        node.leaves = node.leaves
    )
}

.basin.merge.tree.plot.top <- function(layout,
                                       coordinates,
                                       colors,
                                       mass.measure,
                                       support.measure,
                                       show.mass,
                                       show.support,
                                       show.leaf.labels,
                                       height,
                                       main,
                                       field.label,
                                       annotation.cex) {
    branches <- layout$branches
    merge <- layout$merge
    merge.level <- layout$merge.level
    leaf.x <- coordinates$leaf.x
    node.x <- coordinates$node.x
    n.leaf <- nrow(branches)
    value.range <- range(c(branches$birth.level, branches$death.level))
    span <- diff(value.range)
    if (!is.finite(span) || span == 0) {
        span <- max(1, abs(value.range[[1L]]))
    }
    graphics::plot(
        c(0.45, n.leaf + 0.55),
        value.range + c(-0.03, 0.03) * span,
        type = "n",
        xaxt = "n",
        xlab = "",
        ylab = field.label,
        main = main
    )
    graphics::grid(
        nx = NA, ny = NULL, col = "#DCE2E3", lty = 3
    )

    if (nrow(merge) > 0L) {
        graphics::segments(
            leaf.x, branches$birth.level,
            leaf.x, branches$death.level,
            col = colors, lwd = 2.4
        )
        child.index <- match(
            layout$events$losing.basin.id, branches$basin.id
        )
        parent.index <- match(
            layout$events$surviving.basin.id, branches$basin.id
        )
        for (index in seq_len(nrow(merge))) {
            graphics::segments(
                leaf.x[[child.index[[index]]]], merge.level[[index]],
                leaf.x[[parent.index[[index]]]], merge.level[[index]],
                col = colors[[child.index[[index]]]], lwd = 1.4
            )
            graphics::points(
                leaf.x[[parent.index[[index]]]], merge.level[[index]],
                pch = 21, bg = "white", col = "#25343B",
                cex = 0.72, lwd = 1
            )
        }
        root.leaf <- layout$root.leaf
        graphics::points(
            leaf.x[[root.leaf]], branches$death.level[[root.leaf]],
            pch = 22, bg = "white", col = "#25343B", cex = 0.72
        )
    } else {
        graphics::segments(
            leaf.x[[1L]], branches$birth.level[[1L]],
            leaf.x[[1L]], branches$death.level[[1L]],
            col = colors[[1L]], lwd = 2.4
        )
        graphics::points(
            leaf.x[[1L]], branches$death.level[[1L]],
            pch = 22, bg = "white", col = "#25343B", cex = 0.72
        )
    }
    graphics::points(
        leaf.x, branches$birth.level,
        pch = 21, bg = colors, col = colors, cex = 0.82
    )
    if (show.leaf.labels) {
        direction.sign <- if (layout$direction == "max") 1 else -1
        graphics::text(
            leaf.x,
            branches$birth.level + direction.sign * 0.015 * span,
            labels = layout$labels,
            pos = if (layout$direction == "max") 3 else 1,
            cex = 0.67,
            col = colors,
            xpd = NA
        )
    }
    if (!is.null(height)) {
        graphics::abline(
            h = height, col = "#B32239", lty = 2, lwd = 1.3
        )
    }

    usr <- graphics::par("usr")
    annotation.row <- 0L
    if (show.support) {
        annotation.row <- annotation.row + 1L
        y <- usr[[4L]] + span * (0.045 + 0.055 * annotation.row)
        graphics::text(
            leaf.x, y, branches[[support.measure]],
            cex = annotation.cex, xpd = NA
        )
        graphics::text(
            usr[[1L]], y, "Support", pos = 2,
            font = 2, cex = annotation.cex, xpd = NA
        )
    }
    if (show.mass) {
        annotation.row <- annotation.row + 1L
        y <- usr[[4L]] + span * (0.045 + 0.055 * annotation.row)
        graphics::text(
            leaf.x, y,
            .basin.merge.tree.format.annotation(
                branches[[mass.measure]]
            ),
            cex = annotation.cex, xpd = NA
        )
        graphics::text(
            usr[[1L]], y, "Mass", pos = 2,
            font = 2, cex = annotation.cex, xpd = NA
        )
    }
}

.basin.merge.tree.plot.barcode <- function(layout,
                                           colors,
                                           show.guides,
                                           show.birth.labels,
                                           show.parent.labels,
                                           height,
                                           main,
                                           field.label) {
    branches <- layout$branches
    display.index <- rev(layout$order)
    barcode <- branches[display.index, , drop = FALSE]
    barcode.labels <- layout$labels[display.index]
    barcode.colors <- colors[display.index]
    y <- seq_len(nrow(barcode))
    value.range <- range(c(barcode$birth.level, barcode$death.level))
    span <- diff(value.range)
    if (!is.finite(span) || span == 0) {
        span <- max(1, abs(value.range[[1L]]))
    }
    x.padding <- c(
        if (show.parent.labels) 0.35 else 0.04,
        if (show.birth.labels) 0.35 else 0.04
    ) * span
    graphics::plot(
        value.range + c(-x.padding[[1L]], x.padding[[2L]]),
        range(y) + c(-0.65, 0.65),
        type = "n",
        yaxt = "n",
        xlab = field.label,
        ylab = "",
        main = main
    )
    if (show.guides) {
        graphics::abline(
            h = y, col = "#CCD2D3", lty = 3, lwd = 0.8
        )
    }
    graphics::segments(
        barcode$death.level, y, barcode$birth.level, y,
        col = barcode.colors, lwd = 2.6
    )
    root <- is.na(barcode$parent.basin.id)
    graphics::points(
        barcode$death.level[!root], y[!root],
        pch = 21, bg = "white", col = barcode.colors[!root],
        cex = 0.7
    )
    graphics::points(
        barcode$death.level[root], y[root],
        pch = 22, bg = "white", col = barcode.colors[root],
        cex = 0.7
    )
    graphics::points(
        barcode$birth.level, y,
        pch = 21, bg = barcode.colors, col = barcode.colors,
        cex = 0.7
    )
    graphics::axis(
        2, at = y, labels = barcode.labels, las = 1, cex.axis = 0.66
    )
    label.by.id <- setNames(layout$labels, branches$basin.id)
    if (show.birth.labels) {
        graphics::text(
            barcode$birth.level, y, barcode.labels,
            pos = if (layout$direction == "max") 4 else 2,
            cex = 0.62, col = barcode.colors, xpd = NA
        )
    }
    if (show.parent.labels) {
        parent.labels <- rep("root", nrow(barcode))
        parent.labels[!root] <- paste(
            "dies into",
            unname(label.by.id[barcode$parent.basin.id[!root]])
        )
        graphics::text(
            barcode$death.level, y, parent.labels,
            pos = if (layout$direction == "max") 2 else 4,
            cex = 0.58, col = "#59676D", xpd = NA
        )
    }
    if (!is.null(height)) {
        graphics::abline(
            v = height, col = "#B32239", lty = 2, lwd = 1.3
        )
    }
    graphics::legend(
        "bottomright",
        legend = c(
            "Saddle death", "Extremum birth", "Connected root"
        ),
        pch = c(21, 21, 22),
        pt.bg = c("white", "#25343B", "white"),
        col = "#25343B",
        bty = "n",
        cex = 0.66
    )
}

#' Plot a Crossing-Free Basin Merge Tree and Persistence Barcode
#'
#' Draws an exact graph level-set merge tree, its extremum-to-saddle persistence
#' barcode, or both. Tree topology and scalar-field heights come from the
#' canonical merge tree. An `hclust`-compatible representation determines only
#' a deterministic crossing-free leaf order. The rendered tree is directed:
#' every elder basin remains a continuous vertical trunk, and each dying
#' branch terminates horizontally on the trunk that survives its saddle.
#'
#' The optional top-margin rows report one mass and support quantity for each
#' branch. Defaults use the non-overlapping primary assignment. Raw support
#' measures are hierarchical and can overlap across ancestor and descendant
#' branches. In the barcode, dotted guides, labels beside extremum births, and
#' labels identifying the elder basin into which each branch dies can be
#' controlled independently.
#'
#' @param x A `basin.merge.tree`, or a compatible `basin_complex`.
#' @param direction Tree orientation, `"max"` or `"min"`.
#' @param component Graph-component number. It may be omitted when the selected
#'   direction has exactly one graph component.
#' @param type Panels to draw.
#' @param label Canonical basin-table field used for branch labels.
#' @param labels Optional custom labels. A named vector is matched by basin id;
#'   an unnamed vector follows deterministic basin-table order.
#' @param mass.measure Numeric canonical basin-table column shown as mass.
#' @param support.measure Numeric canonical basin-table column shown as support.
#' @param show.mass,show.support Show the aligned top-margin annotation rows.
#' @param show.leaf.labels Label extrema in the tree panel.
#' @param show.barcode.guides Draw a gray dotted guide for every barcode row.
#' @param show.barcode.birth.labels Add basin labels beside filled birth
#'   endpoints.
#' @param show.barcode.parent.labels Add `"dies into <parent>"` beside each
#'   open death endpoint, where `<parent>` is the elder basin that survives the
#'   merge. Root endpoints are labeled `"root"`.
#' @param height Optional exact scalar-field cut shown in both panels.
#' @param branch.col Optional branch colors. A named vector is matched by basin
#'   id; an unnamed vector must have one color per displayed basin.
#' @param main.tree,main.barcode Panel titles.
#' @param field.label Scalar-field axis label.
#' @param annotation.cex Character expansion for mass and support annotations.
#' @param basin.ids Optional canonical basin ids to display. The selection must
#'   be ancestor-closed unless `close.ancestors = TRUE`.
#' @param close.ancestors Add and report missing canonical ancestors of
#'   `basin.ids`.
#' @param ... Reserved for S3 compatibility.
#'
#' @return Invisibly, a list containing the tree, selected branch table,
#'   crossing-free layout, and plotting coordinates.
#' @method plot basin.merge.tree
#' @export
#' @rawNamespace export(plot.basin.merge.tree)
plot.basin.merge.tree <- function(
    x,
    direction = "max",
    component = NULL,
    type = c("tree_and_barcode", "tree", "barcode"),
    label = c("basin.id", "extremum.id", "extremum.vertex"),
    labels = NULL,
    mass.measure = "primary.support.mass",
    support.measure = "primary.support.size",
    show.mass = TRUE,
    show.support = TRUE,
    show.leaf.labels = TRUE,
    show.barcode.guides = TRUE,
    show.barcode.birth.labels = FALSE,
    show.barcode.parent.labels = FALSE,
    height = NULL,
    branch.col = NULL,
    main.tree = "Crossing-free elder-rule merge tree",
    main.barcode = "Extremum-to-saddle persistence barcode",
    field.label = "Scalar-field value",
    annotation.cex = 0.6,
    basin.ids = NULL,
    close.ancestors = FALSE,
    ...
) {
    type <- match.arg(type)
    label <- match.arg(label)
    for (name in c(
        "show.mass", "show.support", "show.leaf.labels",
        "show.barcode.guides", "show.barcode.birth.labels",
        "show.barcode.parent.labels"
    )) {
        .basin.assert.logical(get(name), name)
    }
    if (!is.null(height)) {
        height <- .basin.assert.number(height, "height")
    }
    annotation.cex <- .basin.assert.number(
        annotation.cex, "annotation.cex", lower = 0, lower.open = TRUE
    )
    layout <- get.basin.merge.tree.layout(
        x,
        direction = direction,
        component = component,
        basin.ids = basin.ids,
        close.ancestors = close.ancestors,
        label = label,
        labels = labels
    )
    branches <- layout$branches
    if (!mass.measure %in% names(branches) ||
        !is.numeric(branches[[mass.measure]])) {
        .stop.basin.complex(
            "'mass.measure' must name a numeric basin-table column.",
            "mass.measure"
        )
    }
    if (!support.measure %in% names(branches) ||
        !is.numeric(branches[[support.measure]])) {
        .stop.basin.complex(
            "'support.measure' must name a numeric basin-table column.",
            "support.measure"
        )
    }
    colors <- .basin.merge.tree.colors(branches, branch.col)
    coordinates <- layout$coordinates
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)

    if (type == "tree_and_barcode") {
        graphics::layout(
            matrix(1:2, nrow = 2L), heights = c(1.35, 1)
        )
    }
    if (type %in% c("tree_and_barcode", "tree")) {
        graphics::par(
            family = "sans", las = 1,
            mar = c(2.1, 5.1, if (show.mass || show.support) 6.2 else 3.5, 1.3)
        )
        .basin.merge.tree.plot.top(
            layout = layout,
            coordinates = coordinates,
            colors = colors,
            mass.measure = mass.measure,
            support.measure = support.measure,
            show.mass = show.mass,
            show.support = show.support,
            show.leaf.labels = show.leaf.labels,
            height = height,
            main = main.tree,
            field.label = field.label,
            annotation.cex = annotation.cex
        )
    }
    if (type %in% c("tree_and_barcode", "barcode")) {
        graphics::par(
            family = "sans", las = 1, mar = c(4.2, 8.2, 3.5, 1.3)
        )
        .basin.merge.tree.plot.barcode(
            layout = layout,
            colors = colors,
            show.guides = show.barcode.guides,
            show.birth.labels = show.barcode.birth.labels,
            show.parent.labels = show.barcode.parent.labels,
            height = height,
            main = main.barcode,
            field.label = field.label
        )
    }
    invisible(list(
        tree = layout$tree,
        direction = layout$direction,
        graph.component = layout$component,
        branches = branches,
        layout = layout,
        coordinates = coordinates
    ))
}
