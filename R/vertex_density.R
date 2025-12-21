#' Compute Raw Vertex Density from Nearest Neighbor Distances
#'
#' @description
#' Computes the raw density estimate \eqn{\rho(v) = d_1(v)^{-1}} at each vertex, where
#' \eqn{d_1(v)} is the distance to the nearest neighbor of v (minimum edge weight).
#' This raw estimate is typically noisy and should be smoothed using
#' \code{refit.rdgraph.regression()} to obtain the conditional expectation.
#'
#' @param adj.list Graph adjacency list (1-based indexing).
#' @param weight.list Edge weight list containing edge lengths.
#' @param normalize Logical; if TRUE (default), normalize values to \eqn{[0,1]}.
#'   If FALSE, return raw \eqn{d_1^{-1}} values.
#'
#' @return Numeric vector of density values, one per vertex.
#'
#' @details
#' The nearest neighbor distance \eqn{d_1(v)} is the minimum edge weight among all
#' edges incident to vertex v. In regions of high data density, vertices are
#' close together so \eqn{d_1} is small and \eqn{\rho = d_1^{-1}} is large. In sparse regions,
#' \eqn{d_1} is large and \eqn{\rho} is small.
#'
#' The raw values are noisy estimates. For gradient flow modulation, it is
#' recommended to compute the conditional expectation using the regression
#' framework:
#'
#' \preformatted{
#' raw.density <- compute.vertex.density(adj.list, weight.list, normalize = FALSE)
#' smoothed.density <- refit.rdgraph.regression(fit, raw.density)$fitted.values
#' }
#'
#' @examples
#' \dontrun{
#' # Compute raw density
#' raw.rho <- compute.vertex.density(adj.list, weight.list, normalize = FALSE)
#'
#' # Smooth using conditional expectation
#' density.fit <- refit.rdgraph.regression(fit, raw.rho)
#' smoothed.rho <- density.fit$fitted.values
#'
#' # Use in basin computation
#' basins <- compute.gfc.basins(
#'     adj.list, weight.list, fitted.values,
#'     modulation = "density",
#'     density = smoothed.rho
#' )
#' }
#'
#' @seealso \code{\link{compute.gfc.basins}}, \code{\link{refit.rdgraph.regression}}
#'
#' @export
compute.vertex.density <- function(adj.list,
                                    weight.list,
                                    normalize = TRUE) {

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

    ## Convert adjacency list to 0-based for C++
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            integer(0)
        } else {
            as.integer(x - 1)
        }
    })

    ## Call C++
    result <- .Call(
        S_compute_vertex_density,
        adj.list.0based,
        weight.list,
        as.logical(normalize),
        PACKAGE = "gflow"
    )

    return(result)
}


#' Compute Smoothed Vertex Density for Gradient Flow Modulation
#'
#' @description
#' Convenience function that computes the conditional expectation of \eqn{d_1^{-1}}
#' using the graph regression framework. This provides a smoothed density
#' estimate suitable for gradient flow modulation.
#'
#' @param fit An rdgraph regression fit object (from \code{fit.rdgraph.regression()}).
#' @param adj.list Graph adjacency list (1-based). If NULL, extracted from fit.
#' @param weight.list Edge weight list. If NULL, extracted from fit.
#'
#' @return Numeric vector of smoothed density values (conditional expectation
#'   of d_1^{-1}), normalized to \eqn{[0,1]}.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Computes raw density \eqn{\rho = d_1^{-1}} at each vertex
#'   \item Refits the regression to obtain \eqn{E[\rho | G]}, the conditional expectation
#'   \item Normalizes to \eqn{[0,1]} for use in gradient flow modulation
#' }
#'
#' @examples
#' \dontrun{
#' # Compute smoothed density from existing regression fit
#' smoothed.rho <- compute.smoothed.density(fit)
#'
#' # Use in density-modulated basin computation
#' basins <- compute.gfc.basins(
#'     adj.list, weight.list, fitted.values,
#'     modulation = "density",
#'     density = smoothed.rho
#' )
#' }
#'
#' @seealso \code{\link{compute.vertex.density}}, \code{\link{compute.gfc.basins}}
#'
#' @export
compute.smoothed.density <- function(fit,
                                      adj.list = NULL,
                                      weight.list = NULL) {

    ## Extract graph from fit if not provided
    if (is.null(adj.list)) {
        if (!"adj.list" %in% names(fit)) {
            stop("adj.list must be provided or present in fit object")
        }
        adj.list <- fit$adj.list
    }

    if (is.null(weight.list)) {
        if (!"weight.list" %in% names(fit)) {
            stop("weight.list must be provided or present in fit object")
        }
        weight.list <- fit$weight.list
    }

    ## Step 1: Compute raw density (unnormalized)
    raw.rho <- compute.vertex.density(adj.list, weight.list, normalize = FALSE)

    ## Step 2: Refit regression to get conditional expectation
    density.fit <- refit.rdgraph.regression(fit, raw.rho)

    ## Step 3: Extract fitted values and normalize to [0,1]
    smoothed.rho <- density.fit$fitted.values

    ## Normalize
    rho.min <- min(smoothed.rho)
    rho.max <- max(smoothed.rho)

    if (rho.max > rho.min) {
        smoothed.rho <- (smoothed.rho - rho.min) / (rho.max - rho.min)
    } else {
        smoothed.rho <- rep(1.0, length(smoothed.rho))
    }

    return(smoothed.rho)
}
