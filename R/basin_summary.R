# Ranked canonical basin-complex summaries.

.basin.summary.rank.measures <- c(
    "primary.support.mass",
    "raw.allocated.mass",
    "retained.support.mass",
    "raw.support.mass",
    "primary.support.size",
    "retained.support.size",
    "raw.support.size"
)

.basin.summary.rank.definitions <- c(
    primary.support.mass =
        "Normalized mass of vertices uniquely assigned to the basin.",
    raw.allocated.mass =
        "Conserved membership-weighted raw mass; available only while raw allocation is current.",
    retained.support.mass =
        "Normalized coverage mass of retained support; overlapping basins may double-count mass.",
    raw.support.mass =
        "Normalized coverage mass of raw support; overlapping basins may double-count mass.",
    primary.support.size =
        "Number of vertices uniquely assigned to the basin.",
    retained.support.size =
        "Number of vertices in retained support after refinement.",
    raw.support.size =
        "Number of vertices in raw support before refinement."
)

.basin.summary.top.k <- function(x, argument) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
        x < 0 || (!is.infinite(x) && (!is.finite(x) || x != floor(x)))) {
        .stop.basin.complex(
            sprintf(
                "'%s' must be a nonnegative whole number or Inf.",
                argument
            ),
            argument
        )
    }
    if (is.infinite(x)) Inf else as.integer(x)
}

.basin.summary.measure.availability <- function(rows,
                                                measure,
                                                raw.allocation.current) {
    if (nrow(rows) == 0L) {
        return(list(available = FALSE, reason = "empty_direction"))
    }
    if (!measure %in% names(rows)) {
        return(list(available = FALSE, reason = "missing_column"))
    }
    if (identical(measure, "raw.allocated.mass") &&
        !isTRUE(raw.allocation.current)) {
        return(list(available = FALSE, reason = "stale_raw_allocation"))
    }
    value <- suppressWarnings(as.numeric(rows[[measure]]))
    if (length(value) != nrow(rows) || any(!is.finite(value))) {
        return(list(available = FALSE, reason = "nonfinite_or_partial"))
    }
    if (any(value < 0)) {
        return(list(available = FALSE, reason = "negative_value"))
    }
    if (!any(value > 0)) {
        return(list(available = FALSE, reason = "all_zero"))
    }
    list(available = TRUE, reason = "usable")
}

.basin.summary.column.definitions <- function(table) {
    definitions <- c(
        basin.id = "Stable canonical basin identifier.",
        extremum.id = "Stable canonical identifier of the representative extremum.",
        type = "Flow direction: max for ascending flow or min for descending flow.",
        method = "Canonical basin-construction method.",
        retained = "Whether the basin remains after canonical refinement.",
        extremum.vertex = "Internal integer index of the representative extremum.",
        extremum.vertex.id = "External ID of the representative extremum.",
        extremum.value = "Raw scalar-field value at the representative extremum.",
        birth.level = "Field level at which the basin is born.",
        death.level = "Field level at which the basin merges or dies.",
        persistence = "Birth-to-death field-level persistence.",
        parent.basin.id = "Canonical parent basin identifier when applicable.",
        raw.support.size =
            unname(.basin.summary.rank.definitions[["raw.support.size"]]),
        raw.support.mass =
            unname(.basin.summary.rank.definitions[["raw.support.mass"]]),
        retained.support.size =
            unname(.basin.summary.rank.definitions[["retained.support.size"]]),
        retained.support.mass =
            unname(.basin.summary.rank.definitions[["retained.support.mass"]]),
        primary.support.size =
            unname(.basin.summary.rank.definitions[["primary.support.size"]]),
        primary.support.mass =
            unname(.basin.summary.rank.definitions[["primary.support.mass"]]),
        raw.allocated.mass =
            unname(.basin.summary.rank.definitions[["raw.allocated.mass"]]),
        assignment.status = "Canonical primary-assignment availability status.",
        retention.status = "Canonical reason the basin remains or was removed.",
        rank = "Direction-local rank under the resolved ranking measure.",
        rank.measure = "Direction-specific measure used to compute rank."
    )
    scalar <- names(table)[!vapply(table, is.list, logical(1))]
    missing <- setdiff(scalar, names(definitions))
    if (length(missing) > 0L) {
        fallback <- gsub(".", " ", missing, fixed = TRUE)
        fallback <- paste0(
            toupper(substr(fallback, 1L, 1L)),
            substr(fallback, 2L, nchar(fallback)),
            "."
        )
        names(fallback) <- missing
        definitions <- c(definitions, fallback)
    }
    fields <- scalar
    data.frame(
        field = fields,
        label = gsub(".", " ", fields, fixed = TRUE),
        definition = unname(definitions[fields]),
        unit = ifelse(
            grepl("\\.mass$", fields),
            "normalized mass",
            ifelse(grepl("\\.size$|^rank$", fields), "count", "field-specific")
        ),
        availability = ifelse(
            fields %in% .basin.summary.rank.measures,
            "See measure.availability by direction.",
            "Canonical column."
        ),
        compact = fields %in% c(
            "basin.id", "type", "rank", "extremum.vertex.id",
            "extremum.value", "primary.support.size",
            "primary.support.mass", "raw.allocated.mass"
        ),
        stringsAsFactors = FALSE
    )
}

.basin.summary.direction <- function(table,
                                     direction,
                                     rank.by,
                                     top.k,
                                     raw.allocation.current) {
    rows <- table[table$type == direction, , drop = FALSE]
    availability <- lapply(
        .basin.summary.rank.measures,
        function(measure) .basin.summary.measure.availability(
            rows,
            measure,
            raw.allocation.current
        )
    )
    names(availability) <- .basin.summary.rank.measures
    availability.table <- data.frame(
        measure = .basin.summary.rank.measures,
        available = vapply(availability, `[[`, logical(1), "available"),
        reason = vapply(availability, `[[`, character(1), "reason"),
        stringsAsFactors = FALSE
    )
    if (nrow(rows) == 0L) {
        rows$rank <- integer()
        rows$rank.measure <- character()
        return(list(
            table = rows,
            resolved = NA_character_,
            status = "empty",
            availability = availability.table
        ))
    }
    resolved <- if (identical(rank.by, "auto")) {
        usable <- availability.table$measure[availability.table$available]
        if (length(usable) == 0L) {
            .stop.basin.complex(
                sprintf(
                    "No usable ranking measure is available for direction '%s'.",
                    direction
                ),
                "rank.by",
                class = "gflow_basin_ranking_error",
                details = list(
                    direction = direction,
                    availability = availability.table
                )
            )
        }
        usable[[1L]]
    } else {
        selected <- availability[[rank.by]]
        if (!isTRUE(selected$available)) {
            .stop.basin.complex(
                sprintf(
                    "Ranking measure '%s' is unavailable for direction '%s': %s.",
                    rank.by,
                    direction,
                    selected$reason
                ),
                "rank.by",
                class = "gflow_basin_ranking_error",
                details = list(
                    direction = direction,
                    measure = rank.by,
                    reason = selected$reason
                )
            )
        }
        rank.by
    }
    value <- as.numeric(rows[[resolved]])
    order.index <- order(
        -value,
        as.character(rows$basin.id),
        method = "radix"
    )
    rows <- rows[order.index, , drop = FALSE]
    rows$rank <- seq_len(nrow(rows))
    rows$rank.measure <- rep.int(resolved, nrow(rows))
    if (!is.infinite(top.k)) {
        keep <- seq_len(min(as.integer(top.k), nrow(rows)))
        if (top.k == 0L) {
            keep <- integer()
        }
        rows <- rows[keep, , drop = FALSE]
    }
    list(
        table = rows,
        resolved = resolved,
        status = "ranked",
        availability = availability.table
    )
}

.basin.summary.ranked <- function(object,
                                  rank.by,
                                  top.k.max,
                                  top.k.min,
                                  include.unretained,
                                  include.vertex.lists) {
    object <- .basin.accessor.object(object)
    rank.by <- .basin.assert.choice(
        rank.by,
        c("auto", .basin.summary.rank.measures),
        "rank.by"
    )
    top.k.max <- .basin.summary.top.k(top.k.max, "top.k.max")
    top.k.min <- .basin.summary.top.k(top.k.min, "top.k.min")
    include.unretained <- .basin.assert.logical(
        include.unretained,
        "include.unretained"
    )
    include.vertex.lists <- .basin.assert.logical(
        include.vertex.lists,
        "include.vertex.lists"
    )

    table <- object$basin.table
    if (!include.unretained) {
        table <- table[table$retained, , drop = FALSE]
    }
    raw.current <- isTRUE(object$provenance$allocation$raw.current)
    requested <- .basin.requested.directions(object$direction)
    results <- list()
    for (direction in c("max", "min")) {
        if (direction %in% requested) {
            results[[direction]] <- .basin.summary.direction(
                table,
                direction,
                rank.by,
                if (direction == "max") top.k.max else top.k.min,
                raw.current
            )
        } else {
            empty <- table[FALSE, , drop = FALSE]
            empty$rank <- integer()
            empty$rank.measure <- character()
            results[[direction]] <- list(
                table = empty,
                resolved = NA_character_,
                status = "not_requested",
                availability = data.frame(
                    measure = .basin.summary.rank.measures,
                    available = FALSE,
                    reason = "direction_not_requested",
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    combined <- rbind(results$max$table, results$min$table)
    row.names(combined) <- seq_len(nrow(combined))
    if (!include.vertex.lists) {
        list.columns <- names(combined)[vapply(combined, is.list, logical(1))]
        combined <- combined[, setdiff(names(combined), list.columns), drop = FALSE]
        for (direction in names(results)) {
            results[[direction]]$table <- results[[direction]]$table[
                ,
                setdiff(names(results[[direction]]$table), list.columns),
                drop = FALSE
            ]
        }
    }
    resolved <- c(
        max = results$max$resolved,
        min = results$min$resolved
    )
    status <- c(max = results$max$status, min = results$min$status)
    measure.definition <- vapply(resolved, function(measure) {
        if (is.na(measure)) NA_character_ else
            unname(.basin.summary.rank.definitions[[measure]])
    }, character(1))
    structure(
        list(
            method = object$method,
            direction = object$direction,
            status = object$status,
            n.vertices = object$n.vertices,
            n.components = object$graph.input$validation$n.components,
            n.basins = nrow(object$basin.table),
            n.memberships = nrow(object$membership),
            n.assignments = nrow(object$assignment),
            n.diagnostics = nrow(object$diagnostics),
            has.vertex.mass =
                !is.null(object$field$vertex.mass.normalized),
            rank.requested = rank.by,
            rank.resolved = resolved,
            rank.status = status,
            rank.measure.definition = measure.definition,
            measure.availability = list(
                max = results$max$availability,
                min = results$min$availability
            ),
            raw.allocation.current = raw.current,
            mass.available =
                !is.null(object$field$vertex.mass.normalized),
            mass.provenance = object$provenance$mass,
            build.identity = object$provenance$build.identity,
            top.k.max = top.k.max,
            top.k.min = top.k.min,
            basin.table = combined,
            maxima = results$max$table,
            minima = results$min$table,
            column.definitions =
                .basin.summary.column.definitions(combined)
        ),
        class = c("summary.basin_complex", "list")
    )
}
