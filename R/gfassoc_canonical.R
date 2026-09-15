## Canonical trajectory supports are packed for the existing association kernel.
## This does not create or impersonate an archived basin object.
.gfassoc.pack <- function(x, support.stage) {
    if (!identical(x$status, "ok") || !identical(x$method, "trajectory_flow") ||
        !identical(x$direction, "both")) {
        stop("Canonical association requires a successful trajectory_flow complex with direction = 'both'.", call. = FALSE)
    }
    if (is.null(support.stage)) {
        stop("Choose support.stage = 'raw' or 'retained' explicitly for canonical basins.", call. = FALSE)
    }
    support.stage <- match.arg(support.stage, c("raw", "retained"))
    refinements <- x$parameters$simplify.params
    if (support.stage == "retained" && any(vapply(
        refinements[setdiff(names(refinements), "support.filter")],
        function(p) isTRUE(p$enabled), logical(1)))) {
        stop("Retained association supports only support.filter refinement; use raw supports or reconstruct without other refinements.", call. = FALSE)
    }
    tab <- get.basin.table(x)
    if (support.stage == "retained") tab <- tab[tab$retained, , drop = FALSE]
    values <- x$field$construction.values
    n <- length(x$graph.input$vertex.id)
    if (!is.numeric(values) || length(values) != n || any(!is.finite(values)))
        stop("Canonical construction values must be finite and match the graph.", call. = FALSE)
    supports <- tab[[paste0(support.stage, ".support.vertices")]]
    packed <- lapply(seq_len(nrow(tab)), function(i) {
        v <- supports[[i]]
        e <- tab$extremum.vertex[i]
        if (!is.numeric(v) || anyNA(v) || any(v != floor(v) | v < 1 | v > n) ||
            anyDuplicated(v) || length(e) != 1L || is.na(e) || e < 1 || e > n ||
            !identical(as.numeric(tab$extremum.value[i]), as.numeric(values[e])))
            stop("Invalid canonical support or extremum metadata.", call. = FALSE)
        if ((tab$type[i] == "max" && any(values[v] > values[e])) ||
            (tab$type[i] == "min" && any(values[v] < values[e])))
            stop("Selected support is not bounded by its recorded extremum.", call. = FALSE)
        list(vertex = as.integer(e), value = values[e], hop_idx = 0L,
             basin_df = cbind(as.integer(v), integer(length(v))))
    })
    ids <- list(max = tab$basin.id[tab$type == "max"], min = tab$basin.id[tab$type == "min"])
    list(lmax_basins = packed[tab$type == "max"], lmin_basins = packed[tab$type == "min"],
         n_vertices = n, metadata = list(vertex.id = x$graph.input$vertex.id,
             graph = x$graph.input[c("adj.list", "edge.length.list")],
             construction.values = as.numeric(values), basin.ids = ids,
             support.stage = support.stage, method = x$method))
}

.gfassoc.check.field <- function(y, metadata) {
    if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y)) ||
        !identical(as.numeric(y), metadata$construction.values) ||
        (!is.null(names(y)) && !identical(names(y), metadata$vertex.id)))
        stop("Field values and any vertex names must match canonical construction values in vertex order.", call. = FALSE)
}

.gfassoc.check.pair <- function(a, b) {
    if (!identical(a$vertex.id, b$vertex.id) || !identical(a$graph, b$graph))
        stop("Canonical inputs must use the same graph, lengths, vertex IDs and vertex order; align and reconstruct them first.", call. = FALSE)
    if (!identical(a$support.stage, b$support.stage))
        stop("Canonical memberships must use the same support.stage.", call. = FALSE)
}

.gfassoc.mass <- function(mass, n) {
    if (is.null(mass)) return(rep(1 / n, n))
    if (!is.numeric(mass) || !is.null(dim(mass)) || length(mass) != n ||
        any(!is.finite(mass)) || any(mass < 0) || !is.finite(sum(mass)) || sum(mass) <= 0)
        stop("vertex.mass must contain one finite nonnegative weight per vertex with a positive finite total.", call. = FALSE)
    mass / sum(mass)
}
