#' Detect probabilistic extrema in an rdgraph fitted surface
#'
#' Retained in `gflow` because this is an extrema computation API, even though
#' its input is the archived `riem.dcx` fit class now produced by `gflowx`.
#' The function filters candidate vertices by extremum probability and
#' effective degree, then measures the hop radius over which the probabilistic
#' extremum condition persists.
#'
#' @param dcx A `riem.dcx` object, normally returned by
#'   `gflowx::fit.rdgraph.regression()`.
#' @param extremp_quantile Quantile in `[0, 1]` used to select initial maximum
#'   and minimum candidates.
#' @param hop_extremp_threshold Probability threshold in `(0, 1]` used when
#'   expanding hop neighborhoods.
#' @param eff_degree_quantile Quantile in `[0, 1]` below which weakly connected
#'   candidate vertices are discarded.
#' @param max_hop Positive integer maximum hop distance to explore.
#' @param compute_hop_idx Whether to add the classical hop index computed by
#'   [compute.extrema.hop.nbhds()].
#'
#' @return A data frame with one row per selected extremum and columns
#'   `vertex`, `value`, `rel_value`, `type`, `extremp`, `eff_degree`,
#'   `hop_extremp_radius`, `hop_idx`, and `label`.
#'
#' @seealso [maxp()], [minp()], [compute.extrema.hop.nbhds()]
#' @export
compute.pextrema.nbhds <- function(dcx,
                                   extremp_quantile = 0.95,
                                   hop_extremp_threshold = 0.90,
                                   eff_degree_quantile = 0.10,
                                   max_hop = 20L,
                                   compute_hop_idx = FALSE) {
    if (!inherits(dcx, "riem.dcx")) {
        stop(
            "dcx must be an object of class 'riem.dcx' ",
            "(output from gflowx::fit.rdgraph.regression())"
        )
    }

    required_fields <- c("fitted.values", "graph")
    missing_fields <- setdiff(required_fields, names(dcx$optimal.fit))
    if (length(missing_fields) > 0L) {
        stop(
            "dcx$optimal.fit is missing required fields: ",
            paste(missing_fields, collapse = ", ")
        )
    }

    required_graph_fields <- c(
        "adj.list",
        "vertex.densities",
        "edge.densities"
    )
    graph <- dcx$optimal.fit$graph
    missing_graph_fields <- setdiff(required_graph_fields, names(graph))
    if (length(missing_graph_fields) > 0L) {
        stop(
            "dcx$optimal.fit$graph is missing required fields: ",
            paste(missing_graph_fields, collapse = ", ")
        )
    }

    yhat <- dcx$optimal.fit$fitted.values
    expected_y <- mean(dcx$optimal.fit$y)
    vertex_density <- graph$vertex.densities
    adjacency <- graph$adj.list
    n <- length(adjacency)

    if (length(yhat) != n || length(vertex_density) != n) {
        stop("fitted values and vertex densities must match graph size")
    }
    if (extremp_quantile < 0 || extremp_quantile > 1) {
        stop("extremp_quantile must be in [0, 1]")
    }
    if (hop_extremp_threshold <= 0 || hop_extremp_threshold > 1) {
        stop("hop_extremp_threshold must be in (0, 1]")
    }
    if (eff_degree_quantile < 0 || eff_degree_quantile > 1) {
        stop("eff_degree_quantile must be in [0, 1]")
    }
    if (length(max_hop) != 1L || is.na(max_hop) || max_hop <= 0) {
        stop("max_hop must be a positive integer")
    }

    vertex_maxp <- numeric(n)
    vertex_minp <- numeric(n)
    for (i in seq_len(n)) {
        vertex_maxp[i] <- .maxp_no_check(
            i,
            adjacency[[i]],
            yhat,
            vertex_density
        )
        vertex_minp[i] <- .minp_no_check(
            i,
            adjacency[[i]],
            yhat,
            vertex_density
        )
    }

    effective_degree <- numeric(n)
    edge_index <- 1L
    for (i in seq_len(n)) {
        degree <- length(adjacency[[i]])
        if (degree > 0L) {
            indices <- edge_index:(edge_index + degree - 1L)
            effective_degree[i] <- sum(graph$edge.densities[indices])
            edge_index <- edge_index + degree
        }
    }

    maxp_threshold <- quantile(
        vertex_maxp,
        probs = extremp_quantile,
        na.rm = TRUE
    )
    minp_threshold <- quantile(
        vertex_minp,
        probs = extremp_quantile,
        na.rm = TRUE
    )
    degree_threshold <- quantile(
        effective_degree,
        probs = eff_degree_quantile,
        na.rm = TRUE
    )

    max_candidates <- which(
        vertex_maxp >= maxp_threshold &
            effective_degree >= degree_threshold
    )
    min_candidates <- which(
        vertex_minp >= minp_threshold &
            effective_degree >= degree_threshold
    )

    empty_result <- function() {
        data.frame(
            vertex = integer(),
            value = numeric(),
            rel_value = numeric(),
            type = character(),
            extremp = numeric(),
            eff_degree = numeric(),
            hop_extremp_radius = numeric(),
            hop_idx = numeric(),
            stringsAsFactors = FALSE
        )
    }
    results <- empty_result()

    add_candidates <- function(candidates, type) {
        if (length(candidates) == 0L) {
            return(empty_result())
        }

        radii <- .Call(
            "S_compute_hop_extremp_radii_batch",
            adjacency,
            graph$edge.densities,
            vertex_density,
            as.integer(candidates),
            as.numeric(yhat),
            as.numeric(hop_extremp_threshold),
            identical(type, "max"),
            as.integer(max_hop),
            PACKAGE = "gflow"
        )
        radii[radii == -1L] <- Inf
        scores <- if (identical(type, "max")) {
            vertex_maxp[candidates]
        } else {
            vertex_minp[candidates]
        }

        data.frame(
            vertex = candidates,
            value = yhat[candidates],
            rel_value = yhat[candidates] / expected_y,
            type = type,
            extremp = scores,
            eff_degree = effective_degree[candidates],
            hop_extremp_radius = radii,
            hop_idx = rep(NA_real_, length(candidates)),
            stringsAsFactors = FALSE
        )
    }

    results <- rbind(
        results,
        add_candidates(max_candidates, "max"),
        add_candidates(min_candidates, "min")
    )
    if (nrow(results) == 0L) {
        results$label <- character()
        return(results)
    }

    maxima <- results[results$type == "max", , drop = FALSE]
    minima <- results[results$type == "min", , drop = FALSE]
    if (nrow(maxima) > 0L) {
        maxima <- maxima[order(-maxima$value), , drop = FALSE]
        maxima$label <- paste0("M", seq_len(nrow(maxima)))
    }
    if (nrow(minima) > 0L) {
        minima <- minima[order(minima$value), , drop = FALSE]
        minima$label <- paste0("m", seq_len(nrow(minima)))
    }

    results <- rbind(maxima, minima)
    rownames(results) <- NULL

    if (compute_hop_idx) {
        classical <- compute.extrema.hop.nbhds(
            adjacency,
            graph$edge.length.list,
            yhat
        )
        common <- intersect(
            results$vertex,
            classical$extrema_df$vertex
        )
        for (vertex in common) {
            results$hop_idx[results$vertex == vertex] <-
                classical$extrema_df$hop_idx[
                    classical$extrema_df$vertex == vertex
                ]
        }
    }

    results
}
