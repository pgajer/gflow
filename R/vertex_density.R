#' Compute Raw Vertex Density from Nearest Neighbor Distances
#'
#' @description
#' Computes the raw density estimate \eqn{\rho(v) = d_1(v)^{-1}} at each vertex, where
#' \eqn{d_1(v)} is the distance to the nearest neighbor of v (minimum edge weight).
#' This is an unsmoothed graph-density diagnostic. Conditional-expectation
#' smoothing based on the retired rdgraph estimator is available from
#' `gflowx::compute.smoothed.density()`.
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
#' The raw values are noisy estimates. Callers that require smoothing should
#' estimate it outside `gflow` and pass the resulting density vector to the
#' basin/flow operation.
#'
#' @examples
#' adj.list <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
#' weight.list <- list(c(1, 2), c(1, 1.5), c(2, 1.5))
#' compute.vertex.density(adj.list, weight.list, normalize = FALSE)
#'
#' @seealso \code{\link{compute.gfc.basins}}
#'
#' @noRd
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
