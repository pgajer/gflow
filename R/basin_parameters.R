## Help tables use the constructor's resolved defaults, not a second copy.
.basin.parameter.rd <- function() {
    table <- function(values) {
        items <- vapply(names(values), function(name) {
            value <- paste(utils::capture.output(dput(values[[name]])), collapse = " ")
            paste0("\\item{\\code{", name, "}}{\\code{", value, "}}")
        }, character(1))
        paste0("\\describe{\n", paste(items, collapse = "\n"), "\n}")
    }
    methods <- c("trajectory_flow", "superlevel_merge_tree", "geodesic_reachability", "rtcb", "overlap_cell_complex")
    sections <- vapply(methods, function(method) paste0(
        "\\subsection{", method, "}{", table(.basin.resolve.method.params(method, list(), 100L, NULL)), "}"), character(1))
    stages <- .basin.simplify.defaults()
    paste(c("\\section{Graph defaults}{", table(.basin.graph.defaults()), "}",
            "\\section{Resolved method defaults}{", sections, "}",
            "\\section{Refinement defaults}{", vapply(names(stages), function(stage)
                paste0("\\subsection{", stage, "}{", table(stages[[stage]]), "}"), character(1)), "}"), collapse = "\n")
}

#' Basin Construction Parameters
#'
#' Reference for the named parameter lists accepted by [create.basin.complex()].
#' Tables below are generated from the same defaults and resolvers as the
#' constructor. They show the default for a graph with 100 vertices. The RTCB
#' `n.min` default is `max(20, ceiling(sqrt(n.vertices)))`; inspect
#' `object$parameters` for the exact settings used in your analysis.
#' Unknown names, invalid values, and method-inapplicable settings are errors.
#'
#' @section Graph lengths and validation:
#' `edge.length.symmetry.tolerance` is the nonnegative tolerance for matching
#' lengths on the two orientations of an undirected edge. Neighbor lists and
#' lengths are parallel; zero lengths are allowed. The inherited `weight.list`
#' name in [lcor()] and [lslope()] also supplies edge lengths, not vertex masses
#' or generic affinity weights. Derivative local correlation skips lengths at
#' or below 1e-10. Vertex mass weights support, whereas vertex density can change
#' the trajectory direction. Neither is inferred from the other.
#'
#' @section Trajectory flow:
#' `modulation` selects CLOSEST, NONE, DENSITY, EDGELEN or DENSITY_EDGELEN.
#' CLOSEST uses exact connected plateaus and the nearest improving boundary
#' edge; other modes require `plateau.policy = "none"`. Density modes require
#' `vertex.density`. `edge.length.quantile.thld` in (0,1] controls the preferred
#' local edge-length threshold. `long.edge.fallback` is `"allow_and_flag"`,
#' `"allow"`, or `"forbid"` when only longer improving edges are available.
#' `symmetric.seeding` requests the backend's symmetric seed rule;
#' `max.trajectory.length` bounds stored path length and `store.trajectories`
#' controls path storage. `tie.breaking = TRUE` requires an explicit integer
#' `tie.seed` and perturbs construction values; FALSE requires a NULL seed.
#' `primary.assignment.policy` is `"backend_primary"`,
#' `"largest_membership_then_extremum"`, or `"none"`. It chooses a single
#' convenience label without replacing overlapping membership.
#'
#' @section Merge trees:
#' `plateau.tolerance = 0` is currently required: plateaus use exact equality.
#' The only elder rule is `"birth_then_representative_vertex"`: higher birth
#' wins, with representative vertex breaking equal-birth ties. Assignment is
#' `"elder_at_merge"`; root death uses `"component_opposite_extreme"`.
#' Persistence is in field units and is not a significance score.
#'
#' @section Geodesic and relaxed trajectory basins:
#' Geodesic reachability uses the edge-length quantile in (0,1], optional stored
#' trajectories, and `"none"` or `"largest_membership_then_extremum"` assignment.
#' RTCB exposes stopping and path-search controls: `n.min` limits the initial Dijkstra candidate region
#' (clamped to graph size); `m.min` is the minimum monotonicity in (-1,1); optional `q.min`
#' is the required favorable-change fraction in (0,1] and takes precedence over `m.min`;
#' `run.max` limits consecutive adverse steps.
#' `tau0`, `kappa`, `k.max`, `h.max`, `d.path.max`, and `eta.step` control the
#' path search: local tolerance is tau0*(1-progress)^kappa, k.max bounds
#' labels per vertex, h.max bounds hops, d.path.max bounds metric path length,
#' and eta.step is the field-change deadband for classifying steps. `epsilon.d` and `eps.num` are numerical
#' tolerances. `sink.prune`, `sink.prune.max.iter`, and `max.paths.per.sink`
#' control sink pruning and retained path counts. See
#' the details below for the path-budget interpretation.
#' With cumulative favorable field change G and adverse change L, the
#' path budget is lambda*G - L + tau >= -eps.num. Here lambda is
#' (1-m.min)/(1+m.min), or (1-q.min)/q.min when q.min is supplied.
#' These controls are method-specific, not interchangeable measures of scale.
#'
#' @section Overlap cells:
#' This family requires both directions. `basin.merge.overlap.thld` in `[0,1]`
#' controls basin merging. The three `min.*.cell.size.thld` integers control
#' minimum ascending/descending, ascending/ascending, and descending/descending
#' cell sizes. `cell.graph.params` is forwarded to the cell graph backend.
#' Primary assignment is `"none"`; overlapping cells are the authoritative result.
#'
#' @section Refinement applicability and order:
#' Every stage is disabled by default. Relative-value filtering, maxima and
#' minima clustering, and geometric filtering apply to trajectory, geodesic and
#' RTCB methods. Support filtering also applies to merge trees. Expansion applies
#' only to geodesic and RTCB methods. No refinement stage applies to overlap cells.
#' Nondefault settings in inapplicable stages are rejected even if disabled.
#'
#' Stages run in the order above, followed by support filtering and expansion.
#' Relative thresholds compare extrema with the mean construction field;
#' clustering uses the overlap threshold in `[0,1]`. Geometric filters use the
#' listed percentiles and positive integer hop radius. Support filters require
#' minimum vertex count, trajectory count and normalized support mass; mass
#' thresholds need supplied mass. Expansion assigns uncovered vertices to the
#' nearest retained basin. Raw supports remain available independently of these
#' retained results. See the canonical workflow for worked refinement examples.
#'
#' @evalRd .basin.parameter.rd()
#' @seealso [create.basin.complex()], [get.basin.table()], [gfcor()],
#'   [gflow-package]. Run `vignette("function-guide", package = "gflow")` for
#'   the task map and `vignette("basin_complex_workflow_vignette", package = "gflow")`
#'   for construction and refinement examples.
#' @name basin-parameters
NULL
