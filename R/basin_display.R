# Presentation helpers consume canonical data without changing construction.
.basin.coverage <- function(object) {
    rows <- lapply(.basin.requested.directions(object$direction), function(direction) {
        table <- object$basin.table
        retained <- table[table$type == direction & table$retained, , drop = FALSE]
        membership <- object$membership
        assignment <- object$assignment
        data.frame(
            direction = direction,
            retained.basins = nrow(retained),
            raw.vertices = length(unique(membership$vertex[membership$direction == direction])),
            retained.vertices = length(unique(unlist(retained$retained.support.vertices))),
            assigned.vertices = length(unique(assignment$vertex[
                assignment$direction == direction & !is.na(assignment$basin.id)])),
            total.vertices = object$n.vertices,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

.print.basin.summary <- function(x, n, digits) {
    n <- .basin.summary.top.k(n, "n")
    digits <- .basin.assert.number(digits, "digits", lower = 1, upper = 22)
    cat("Canonical Basin Complex Summary\n")
    cat("  Method: ", x$method, "\n  Direction: ", x$direction,
        "\n  Status: ", x$status, "\n", sep = "")
    cat("  Vertices/components: ", x$n.vertices, "/", x$n.components, "\n", sep = "")
    cat("  Basin rows: ", x$n.basins, "; membership rows: ", x$n.memberships,
        "; assignment rows: ", x$n.assignments, "\n", sep = "")
    if (!is.null(x$diagnostics) && nrow(x$diagnostics)) {
        cat("  Diagnostics: ", paste(unique(x$diagnostics$message), collapse = "; "), "\n", sep = "")
    }
    if (length(x$warnings)) cat("  Warnings: ", paste(x$warnings, collapse = "; "), "\n", sep = "")
    if (!is.null(x$coverage)) {
        cat("\nCoverage (unique vertices; assignment rows may have no label):\n")
        print(x$coverage, row.names = FALSE)
    }
    for (direction in c("max", "min")) {
        table <- if (direction == "max") x$maxima else x$minima
        measure <- x$rank.resolved[[direction]]
        if (!nrow(table) || is.na(measure) || n == 0) next
        shown <- utils::head(table, n = min(n, nrow(table)))
        columns <- unique(c("rank", "extremum.vertex.id", measure, "persistence", "retained"))
        columns <- intersect(columns, names(shown))
        cat("\n", if (direction == "max") "Maxima" else "Minima", ": ",
            nrow(shown), " of ", nrow(table), " returned basins; ranked by ", measure, "\n", sep = "")
        definition <- x$rank.measure.definition[[direction]]
        if (!is.na(definition)) cat("  ", definition, "\n", sep = "")
        print(shown[, columns, drop = FALSE], row.names = FALSE, digits = as.integer(digits))
    }
    cat("\nPersistence is in field units; mass measures are normalized support fractions.\n")
    invisible(x)
}

.plot.basin.view <- function(x, view, direction, coordinates, stage,
                              label.vertices, xlab, ylab, main, ...) {
    x <- .basin.accessor.object(x)
    label.vertices <- .basin.assert.logical(label.vertices, "label.vertices")
    suffix <- if (x$status == "ok") "" else paste0(" [construction ", x$status, "]")
    diagnostic <- if (nrow(x$diagnostics)) paste(unique(x$diagnostics$message), collapse = "; ") else "See object diagnostics."
    if (x$status == "failed" && view != "field") {
        stop("Cannot plot a basin result: construction failed. ", diagnostic,
             " Use view = 'field' to inspect the input.", call. = FALSE)
    }
    if (x$status != "ok") warning("Construction ", x$status, ": ", diagnostic, call. = FALSE)
    if (view == "field") {
        graphics::plot(seq_len(x$n.vertices), x$field$input.values,
                       xlab = if (is.null(xlab)) "Vertex index (graph order)" else xlab,
                       ylab = if (is.null(ylab)) "Scalar field" else ylab,
                       main = paste0(if (is.null(main)) "Input scalar field" else main, suffix), ...)
        return(invisible(x))
    }
    if (is.null(direction)) {
        if (x$direction == "both") stop("Choose direction = 'max' or 'min' for this view.", call. = FALSE)
        direction <- x$direction
    }
    direction <- .basin.assert.choice(direction, c("max", "min"), "direction")
    if (!direction %in% .basin.requested.directions(x$direction)) {
        stop("The requested direction was not constructed.", call. = FALSE)
    }
    if (view == "merge_tree") {
        args <- list(...)
        if (is.null(args$main.tree)) args$main.tree <- if (is.null(main)) "Basin merge tree" else main
        args$main.tree <- paste0(args$main.tree, suffix)
        if (is.null(args$main.barcode)) args$main.barcode <- "Persistence in field units"
        args$main.barcode <- paste0(args$main.barcode, suffix)
        do.call(plot, c(list(x = get.basin.merge.tree(x, required = TRUE), direction = direction), args))
        return(invisible(x))
    }
    if (!is.matrix(coordinates) || !is.numeric(coordinates) ||
        !identical(dim(coordinates), c(x$n.vertices, 2L)) || any(!is.finite(coordinates))) {
        stop("coordinates must be a finite numeric n-by-2 matrix in graph vertex order.", call. = FALSE)
    }
    if (!is.null(rownames(coordinates)) &&
        !identical(rownames(coordinates), x$graph.input$vertex.id)) {
        stop("coordinates row names must match external vertex IDs in graph order.", call. = FALSE)
    }
    groups <- .basin.plot.groups(x, view, direction, stage)
    colors <- grDevices::hcl.colors(max(3L, length(groups$legend)), "Dark 3")[seq_along(groups$legend)]
    colors[1L] <- "grey85"
    dots <- list(...)
    defaults <- list(x = coordinates[, 1L], y = coordinates[, 2L], type = "n", asp = 1,
                     ylim = grDevices::extendrange(coordinates[, 2L], f = c(0.1, 0.35)),
                     xlab = if (!is.null(xlab)) xlab else if (!is.null(colnames(coordinates))) colnames(coordinates)[1L] else "Drawing coordinate 1",
                     ylab = if (!is.null(ylab)) ylab else if (!is.null(colnames(coordinates))) colnames(coordinates)[2L] else "Drawing coordinate 2",
                     main = paste0(if (is.null(main)) groups$title else main, suffix))
    defaults[names(dots)] <- dots
    do.call(graphics::plot, defaults)
    for (i in seq_along(x$graph.input$adj.list)) {
        neighbors <- x$graph.input$adj.list[[i]]
        neighbors <- neighbors[neighbors > i]
        if (length(neighbors)) graphics::segments(coordinates[i, 1L], coordinates[i, 2L],
                                                coordinates[neighbors, 1L], coordinates[neighbors, 2L], col = "grey75")
    }
    graphics::points(coordinates, pch = ifelse(groups$index == 1L, 4, 21),
                     bg = colors[groups$index], cex = if (is.null(dots$cex)) 1.2 else dots$cex)
    if (label.vertices) graphics::text(coordinates, labels = x$graph.input$vertex.id, pos = 3, cex = 0.7)
    graphics::legend("topright", legend = groups$legend, pt.bg = colors,
                     pch = c(4, rep(21, length(colors) - 1L)), bty = "n", cex = 0.7)
    invisible(x)
}

.basin.plot.groups <- function(x, view, direction, stage) {
    if (view == "assignment") {
        assignment <- x$assignment[x$assignment$direction == direction, , drop = FALSE]
        if (!nrow(assignment) || all(assignment$assignment.status == "not_applicable")) {
            stop("Primary assignment is unavailable for this method; use view = 'overlap' to inspect support.", call. = FALSE)
        }
        ids <- unique(assignment$basin.id[!is.na(assignment$basin.id)])
        index <- rep(1L, x$n.vertices)
        index[assignment$vertex] <- match(assignment$basin.id, ids, nomatch = 0L) + 1L
        peaks <- x$basin.table$extremum.vertex.id[match(ids, x$basin.table$basin.id)]
        return(list(index = index, legend = c("No assignment", paste(if (direction == "max") "Peak" else "Minimum", peaks)),
                    title = paste("Primary basin assignments:", direction)))
    }
    if (stage == "raw") {
        membership <- x$membership[x$membership$direction == direction, , drop = FALSE]
        pairs <- unique(membership[, c("vertex", "basin.id"), drop = FALSE])
        count <- tabulate(pairs$vertex, nbins = x$n.vertices)
    } else {
        table <- x$basin.table[x$basin.table$type == direction & x$basin.table$retained, , drop = FALSE]
        count <- tabulate(as.integer(unlist(lapply(table$retained.support.vertices, unique))), nbins = x$n.vertices)
    }
    levels <- sort(unique(c(0L, count)))
    list(index = match(count, levels), legend = paste(levels, ifelse(levels == 1L, "basin", "basins")),
         title = paste(if (stage == "raw") "Raw" else "Retained", "support overlap:", direction))
}
