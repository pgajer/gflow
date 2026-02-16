#' Compute Gradient Flow Basins with Optional Modulation
#'
#' @description
#' Computes basins of attraction for all local extrema of a function defined
#' on a weighted graph. Supports modulated gradient flow that accounts for
#' data density and edge length distributions to produce more robust basins.
#'
#' @param adj.list A list where each element i contains the indices of vertices
#'   adjacent to vertex i (1-based indexing).
#' @param weight.list A list where each element i contains the edge weights
#'   (lengths) from vertex i to its neighbors.
#' @param y Numeric vector of function values at each vertex.
#' @param modulation Character string specifying gradient flow modulation:
#'   \describe{
#'     \item{"none"}{Standard gradient flow (steepest ascent/descent)}
#'     \item{"density"}{Density-modulated: \eqn{\rho(u) \cdot \Delta \hat{y}}, favors high-density regions}
#'     \item{"edgelen"}{Edge-length-modulated: \eqn{dl([v,u]) \cdot \Delta \hat{y}}, penalizes atypical edges}
#'     \item{"density_edgelen"}{Combined: \eqn{\rho(u) \cdot dl([v,u]) \cdot \Delta \hat{y}}}
#'   }
#' @param density Optional numeric vector of density values at each vertex.
#'   If NULL and modulation includes "density", computed from nearest neighbor
#'   distances as \eqn{\rho(v) = 1/d_1(v)}.
#' @param edgelen.bandwidth Bandwidth for edge length KDE. If negative (default),
#'   uses Silverman's rule of thumb.
#' @param verbose Logical; if TRUE, print progress information.
#'
#' @return An object of class \code{"gfc_basins"}, a named list where each
#'   element corresponds to a local extremum and contains:
#'   \describe{
#'     \item{extremum.vertex}{Vertex index of the local extremum (1-based)}
#'     \item{value}{Function value at the extremum}
#'     \item{is.maximum}{Logical; TRUE for maxima, FALSE for minima}
#'     \item{vertices}{Integer vector of vertices in the basin (1-based)}
#'     \item{boundary}{Integer vector of boundary vertices (1-based)}
#'   }
#'   Names are formatted as "min_V" or "max_V" where V is the vertex index.
#'
#' @details
#' The gradient flow at vertex v follows the edge \eqn{[v,u]} that maximizes the
#' modulated gradient score. For ascending flow (to maxima), this is:
#'
#' \itemize{
#'   \item none: \eqn{\max(\hat{y}(u) - \hat{y}(v))}
#'   \item density: \eqn{\max(\rho(u) · (\hat{y}(u) - \hat{y}(v)))}
#'   \item edgelen: \eqn{\max(dl([v,u]) · (\hat{y}(u) - \hat{y}(v)))}
#'   \item density_edgelen: \eqn{\max(\rho(u) · dl([v,u]) · (\hat{y}(u) - \hat{y}(v)))}
#' }
#'
#' The density modulation favors flow through high-density regions, avoiding
#' spurious paths through sparse areas. The edge length modulation penalizes
#' very short and very long edges, favoring edges near the mode of the length
#' distribution. This addresses the "long edge basin jumping" problem.
#'
#' Boundary vertices are neighbors of basin vertices that do not belong to
#' the basin. These are essential for Dirichlet harmonic extension problems.
#'
#' @examples
#' \dontrun{
#' # Compute standard basins
#' basins <- compute.gfc.basins(adj.list, weight.list, fitted.values,
#'                               modulation = "none")
#'
#' # Compute density-modulated basins
#' basins.density <- compute.gfc.basins(adj.list, weight.list, fitted.values,
#'                                       modulation = "density")
#'
#' # Access a specific basin
#' m1.basin <- basins[["max_1147"]]
#' print(m1.basin$vertices)
#' print(m1.basin$boundary)
#' }
#'
#' @export
compute.gfc.basins <- function(adj.list,
                               weight.list,
                               y,
                               modulation = c("none", "density", "edgelen", "density_edgelen"),
                               density = NULL,
                               edgelen.bandwidth = -1.0,
                               verbose = TRUE) {

    ## Input validation
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != n.vertices) {
        stop("y must have the same length as adj.list")
    }

    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    modulation <- match.arg(modulation)

    ## Validate density if provided
    if (!is.null(density)) {
        if (!is.numeric(density) || length(density) != n.vertices) {
            stop("density must be a numeric vector of length n.vertices")
        }
    }

    ## Convert adjacency list to 0-based for C++
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            integer(0)
        } else {
            as.integer(x - 1)
        }
    })

    ## Prepare density argument
    if (is.null(density)) {
        density.arg <- numeric(0)
    } else {
        density.arg <- as.numeric(density)
    }

    ## Call C++
    result <- .Call(
        S_compute_gfc_basins,
        adj.list.0based,
        weight.list,
        as.numeric(y),
        modulation,
        density.arg,
        as.numeric(edgelen.bandwidth),
        as.logical(verbose),
        PACKAGE = "gflow"
    )

    ## Add metadata
    attr(result, "n.vertices") <- n.vertices
    attr(result, "modulation") <- modulation

    return(result)
}

#' Print Method for gfc_basins Objects
#'
#' @param x A gfc_basins object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gfc_basins <- function(x, ...) {
    cat("Gradient Flow Basins\n")
    cat("====================\n")
    cat("Modulation:", attr(x, "modulation"), "\n")
    cat("Total vertices:", attr(x, "n.vertices"), "\n")

    n.min <- sum(sapply(x, function(b) !b$is.maximum))
    n.max <- sum(sapply(x, function(b) b$is.maximum))

    cat(sprintf("Basins: %d minima, %d maxima\n\n", n.min, n.max))

    ## Summary table
    df <- data.frame(
        name = names(x),
        vertex = sapply(x, function(b) b$extremum.vertex),
        type = sapply(x, function(b) ifelse(b$is.maximum, "max", "min")),
        value = sapply(x, function(b) round(b$value, 4)),
        n.vertices = sapply(x, function(b) length(b$vertices)),
        n.boundary = sapply(x, function(b) length(b$boundary))
    )

    ## Sort by type then value
    df <- df[order(df$type, -df$value), ]
    rownames(df) <- NULL

    print(df, row.names = FALSE)

    invisible(df)
}

#' Extract Basin by Extremum Vertex or Name
#'
#' @param basins A gfc_basins object
#' @param query Either a vertex index (integer) or basin name (character)
#'
#' @return The basin structure, or NULL if not found
#'
get.basin.of.gfc.basins <- function(basins, query) {
    if (is.character(query)) {
        return(basins[[query]])
    } else if (is.numeric(query)) {
        vertex <- as.integer(query)
        ## Search for basin with this extremum vertex
        for (basin in basins) {
            if (basin$extremum.vertex == vertex) {
                return(basin)
            }
        }
        return(NULL)
    } else {
        stop("query must be a character name or integer vertex index")
    }
}
