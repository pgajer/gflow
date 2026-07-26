# Archived basin-object conversion.

.basin.conversion.graph <- function(adj.list, edge.length.list) {
    if (is.null(adj.list) || is.null(edge.length.list)) {
        .stop.basin.complex(
            paste(
                "Archived basin conversion requires 'adj.list' and",
                "'edge.length.list' because legacy objects do not retain a",
                "complete canonical graph input."
            ),
            "adj.list",
            class = "gflow_basin_conversion_error"
        )
    }
    list(adj.list = adj.list, edge.length.list = edge.length.list)
}

.basin.record.conversion <- function(object, archived, source.class) {
    object$raw$archived.object <- archived
    object$provenance$conversion <- list(
        source.class = source.class,
        mode = "reconstructed_from_archived_field_and_supplied_graph",
        archived.object.preserved = TRUE
    )
    .validate.basin.complex(object)
    object
}

#' @method as.basin.complex gfc.flow
#' @export
as.basin.complex.gfc.flow <- function(
    object,
    adj.list = NULL,
    edge.length.list = NULL,
    vertex.mass = NULL,
    vertex.density = NULL,
    direction = "both",
    method.params = list(),
    simplify.params = list(),
    verbose = FALSE,
    ...
) {
    graph <- .basin.conversion.graph(adj.list, edge.length.list)
    inferred <- list(
        modulation = .basin.or(object$modulation, "CLOSEST"),
        store.trajectories = !is.null(object$trajectories),
        tie.breaking = FALSE,
        primary.assignment.policy = "backend_primary"
    )
    method.params <- utils::modifyList(inferred, method.params)
    converted <- create.basin.complex(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        field = object$y,
        method = "trajectory_flow",
        direction = direction,
        vertex.mass = vertex.mass,
        vertex.density = vertex.density,
        method.params = method.params,
        simplify.params = simplify.params,
        verbose = verbose
    )
    .basin.record.conversion(converted, object, "gfc.flow")
}

#' @method as.basin.complex basins_of_attraction
#' @export
as.basin.complex.basins_of_attraction <- function(
    object,
    adj.list = NULL,
    edge.length.list = NULL,
    field = NULL,
    vertex.mass = NULL,
    direction = "both",
    method.params = list(),
    simplify.params = list(),
    verbose = FALSE,
    ...
) {
    graph <- .basin.conversion.graph(adj.list, edge.length.list)
    field <- .basin.or(object$y, field)
    if (is.null(field)) {
        .stop.basin.complex(
            paste(
                "This archived basins_of_attraction object does not retain",
                "its construction field; supply 'field' explicitly."
            ),
            "field",
            class = "gflow_basin_conversion_error"
        )
    }
    converted <- create.basin.complex(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        field = field,
        method = "geodesic_reachability",
        direction = direction,
        vertex.mass = vertex.mass,
        method.params = method.params,
        simplify.params = simplify.params,
        verbose = verbose
    )
    .basin.record.conversion(
        converted,
        object,
        "basins_of_attraction"
    )
}

#' @method as.basin.complex basin_cx
#' @export
as.basin.complex.basin_cx <- function(
    object,
    adj.list = NULL,
    edge.length.list = NULL,
    vertex.mass = NULL,
    method.params = list(),
    verbose = FALSE,
    ...
) {
    graph <- .basin.conversion.graph(adj.list, edge.length.list)
    converted <- create.basin.complex(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        field = object$original_y,
        method = "overlap_cell_complex",
        direction = "both",
        vertex.mass = vertex.mass,
        method.params = method.params,
        verbose = verbose
    )
    .basin.record.conversion(converted, object, "basin_cx")
}
