#' Extract Cell Trajectories
#'
#' Generic for extracting gradient-flow trajectories for a specified
#' minimum/maximum cell from supported objects.
#'
#' @param x An object containing trajectory information.
#' @param min.id Minimum extremum identifier.
#' @param max.id Maximum extremum identifier.
#' @param ... Additional arguments passed to methods.
#'
#' @keywords internal
cell.trajectories <- function(x, min.id, max.id, ...) {
    UseMethod("cell.trajectories")
}
