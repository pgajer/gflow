#' Select Graph Endpoints by Core-Eccentricity Geometry
#'
#' Builds an endpoint set from graph geometry in three steps:
#' 1) estimate vertex eccentricity, 2) define a graph core by low eccentricity
#' quantile, 3) pick endpoint candidates as local maxima of distance-to-core.
#'
#' @param adj.list Graph adjacency list (1-based vertex indices).
#' @param weight.list Edge-length list aligned with `adj.list`.
#' @param core.quantile Numeric in (0,1). Vertices with eccentricity in this
#'   lower quantile define the graph core.
#' @param endpoint.quantile Numeric in \eqn{[0,1]}. Endpoints are selected among
#'   local maxima of distance-to-core that are at or above this quantile.
#' @param use.approx.eccentricity Logical. If `TRUE`, use landmark-based
#'   approximation of eccentricity for scalability.
#' @param n.landmarks Integer >= 1. Number of landmarks when approximation is used.
#' @param max.endpoints Optional positive integer cap on returned endpoints.
#' @param seed Integer seed controlling landmark initialization.
#' @param verbose Logical; print progress from C++ backend.
#'
#' @return A list with components:
#' \describe{
#'   \item{endpoints}{Selected endpoint vertices (1-based).}
#'   \item{core.vertices}{Core vertices (1-based).}
#'   \item{eccentricity}{Per-vertex eccentricity estimate.}
#'   \item{distance.to.core}{Per-vertex geodesic distance to core.}
#'   \item{summary}{Per-vertex diagnostics data frame.}
#' }
#'
#' @export
geodesic.core.endpoints <- function(adj.list,
                                    weight.list,
                                    core.quantile = 0.10,
                                    endpoint.quantile = 0.90,
                                    use.approx.eccentricity = TRUE,
                                    n.landmarks = 64L,
                                    max.endpoints = NULL,
                                    seed = 1L,
                                    verbose = FALSE) {
    if (!is.list(adj.list)) stop("'adj.list' must be a list.")
    if (!is.list(weight.list)) stop("'weight.list' must be a list.")
    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length.")
    }
    if (!is.numeric(core.quantile) || length(core.quantile) != 1L ||
        !is.finite(core.quantile) || core.quantile <= 0 || core.quantile >= 1) {
        stop("'core.quantile' must be a finite scalar in (0, 1).")
    }
    if (!is.numeric(endpoint.quantile) || length(endpoint.quantile) != 1L ||
        !is.finite(endpoint.quantile) || endpoint.quantile < 0 || endpoint.quantile > 1) {
        stop("'endpoint.quantile' must be a finite scalar in [0, 1].")
    }
    if (!is.logical(use.approx.eccentricity) || length(use.approx.eccentricity) != 1L) {
        stop("'use.approx.eccentricity' must be a scalar logical.")
    }
    if (!is.numeric(n.landmarks) || length(n.landmarks) != 1L || !is.finite(n.landmarks) || n.landmarks < 1) {
        stop("'n.landmarks' must be a finite scalar >= 1.")
    }
    if (!is.null(max.endpoints)) {
        if (!is.numeric(max.endpoints) || length(max.endpoints) != 1L ||
            !is.finite(max.endpoints) || max.endpoints < 1) {
            stop("'max.endpoints' must be NULL or a finite scalar >= 1.")
        }
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
        stop("'seed' must be a finite scalar.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("'verbose' must be a scalar logical.")
    }

    adj.list.0 <- lapply(adj.list, function(x) as.integer(x - 1L))
    max.endpoints.int <- if (is.null(max.endpoints)) 0L else as.integer(max.endpoints)

    res <- .Call(
        "S_geodesic_core_endpoints",
        adj.list.0,
        weight.list,
        as.double(core.quantile),
        as.double(endpoint.quantile),
        as.logical(use.approx.eccentricity),
        as.integer(n.landmarks),
        as.integer(max.endpoints.int),
        as.integer(seed),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    if (!is.null(res$endpoints)) {
        res$endpoints <- as.integer(res$endpoints) + 1L
    }
    if (!is.null(res$core_vertices)) {
        res$core_vertices <- as.integer(res$core_vertices) + 1L
    }
    if (!is.null(res$landmarks)) {
        res$landmarks <- as.integer(res$landmarks) + 1L
    }
    if (!is.null(res$summary) && is.data.frame(res$summary) && "vertex" %in% names(res$summary)) {
        res$summary$vertex <- as.integer(res$summary$vertex) + 1L
    }

    names(res)[names(res) == "core_vertices"] <- "core.vertices"
    names(res)[names(res) == "distance_to_core"] <- "distance.to.core"
    names(res)[names(res) == "is_core"] <- "is.core"
    names(res)[names(res) == "is_endpoint"] <- "is.endpoint"
    names(res)[names(res) == "is_local_max"] <- "is.local.max"
    names(res)[names(res) == "endpoint_rank"] <- "endpoint.rank"
    names(res)[names(res) == "core_threshold"] <- "core.threshold"
    names(res)[names(res) == "endpoint_threshold"] <- "endpoint.threshold"
    names(res)[names(res) == "used_approx_eccentricity"] <- "used.approx.eccentricity"
    names(res)[names(res) == "n_landmarks_used"] <- "n.landmarks.used"

    class(res) <- c("geodesic_core_endpoints", class(res))
    res
}
