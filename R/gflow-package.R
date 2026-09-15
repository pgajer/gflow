#' gflow: Basin-Complex Exploration and Flow-Aware Association
#'
#' `gflow` constructs canonical basin and gradient-flow complexes and provides
#' stable methods for exploring their basins, cells, trajectories, and
#' associations. Generic graph construction and algorithms are supplied by
#' `dgraphs`; retired smoothing and conditional-expectation estimators are not
#' part of this package.
#'
#' @section Start here:
#' Start with [create.basin.complex()] and inspect the result with
#' [summary.basin_complex()] and the `get.basin.*()` accessors. The small
#' example below needs only the installed package: its two peaks merge at
#' field height 1. Persistence measures lifetime across field thresholds,
#' not statistical significance. Use `plot(bc, view = "merge_tree")` to explore
#' the hierarchy, or change the field and reconstruct it.
#'
#' The documented source-build installation generates help and all four
#' vignettes. A direct Git installation that skips documentation generation
#' may lack them; follow the repository's installation instructions to build
#' a complete source archive. Package loading does not open guides automatically.
#'
#' @section User guides:
#' Run `vignette("function-guide", package = "gflow")` for the task map and
#' complete API catalog, and `vignette("example-graphs-and-fields", package =
#' "gflow")` for reproducible graph/field recipes. The canonical basin and
#' noisy-circle workflow vignettes provide longer analyses.
#' See \code{\link{gflow-migration}} for installed migration guidance.
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
#'   \item \code{\link{gfcor}} provides flow-aware association on canonical
#'     trajectory complexes with both directions and explicit raw/retained
#'     supports, or compatible archived \code{basins_of_attraction} inputs.
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
#' @examples
#' adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
#' lengths <- lapply(adjacency, function(v) rep(1, length(v)))
#' bc <- create.basin.complex(adjacency, lengths, c(0, 3, 1, 2, 0),
#'                            method = "superlevel_merge_tree", direction = "max")
#' get.basin.table(bc)[, c("extremum.vertex", "persistence", "raw.support.size")]
#' # Peaks at vertices 2 and 4 have persistence 3 and 1 field units.
#' # The root's support contains all five vertices; supports can overlap.
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
