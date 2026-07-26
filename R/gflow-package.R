#' gflow: Basin-Complex Exploration and Flow-Aware Association
#'
#' `gflow` constructs canonical basin and gradient-flow complexes and provides
#' stable methods for exploring their basins, cells, trajectories, and
#' associations. Generic graph construction and algorithms are supplied by
#' `dgraphs`; retired smoothing and conditional-expectation estimators are not
#' part of this package.
#'
#' @section Basin and flow objects:
#' \itemize{
#'   \item \code{\link{create.basin.complex}} constructs a canonical
#'     \code{basin_complex}.
#'   \item \code{\link{as.basin.complex}} converts supported archived objects.
#'   \item \code{\link{get.basin.table}}, \code{\link{get.basin.membership}},
#'     and \code{\link{get.basin.assignment}} expose stable tables.
#'   \item \code{\link{get.basin.merge.tree}},
#'     \code{\link{get.basin.trajectory.forest}}, and
#'     \code{\link{get.basin.cells}} expose method-specific structure.
#' }
#'
#' @section Exploration and analysis:
#' \itemize{
#'   \item \code{\link{construct.gflow.graph}} and \code{\link{construct.madag}}
#'     summarize flow structure.
#'   \item \code{\link{compute.harmonic.extension}} extends trajectory
#'     coordinates into their neighborhoods.
#'   \item \code{\link{compute.gfc.modulation}} and
#'     \code{\link{extremality.summary}} support complex interpretation.
#' }
#'
#' @section Local and flow-aware association:
#' \itemize{
#'   \item \code{\link{lcor}} and \code{\link{lslope}} provide graph-local
#'     association and directed response.
#'   \item \code{\link{gfcor}} provides basin- and flow-aware association.
#'   \item \code{\link{gfassoc.membership}},
#'     \code{\link{gfassoc.polarity}}, \code{\link{gfassoc.overlap}}, and
#'     \code{\link{gfassoc.deviation}} expose flow-aware components.
#' }
#'
#' @section Package boundaries:
#' Build graphs with `dgraphs` and pass adjacency and weight lists to `gflow`.
#' Use `gflowx` only when reproducing archived graph-regression or smoothing
#' analyses. PHATE, diffusion/potential pseudotime, generic graph utilities,
#' interactive selection widgets, and domain-specific pipelines are not part
#' of the supported `gflow` surface.
#'
#' @author Pawel Gajer \email{pgajer@@gmail.com}
#'
#' @references
#' Gerber, S., Rubel, O., Bremer, P. T., Pascucci, V., & Whitaker, R. T. (2013).
#' Morse--Smale regression. Journal of Computational and Graphical Statistics,
#' 22(1), 193-214.
#'
#' Chen, Y. C., Genovese, C. R., & Wasserman, L. (2017). Statistical inference
#' using the Morse-Smale complex.
#'
#' @import stats
#' @import graphics
#' @import grDevices
"_PACKAGE"
