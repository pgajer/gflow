#' Migrating to the Basin and Flow API
#'
#' @description
#' For new basin analyses use \code{\link{create.basin.complex}} with explicit
#' adjacency, aligned edge lengths, and a scalar field. The function guide
#' vignette contains a complete task map and status-labeled export catalog.
#'
#' @section Retired constructors:
#' `compute.gfc()`, `compute.basins.of.attraction()`, `compute.gfc.trajectory()`,
#' `compute.gfc.flow()`, and `create.basin.cx()` raise migration errors.
#' Choose canonical `geodesic_reachability`, `trajectory_flow`, `rtcb`,
#' `superlevel_merge_tree`, or `overlap_cell_complex` as appropriate, and
#' translate refinement controls into `simplify.params`.
#' There is no `merge_tree` alias. Arguments are `edge.length.list` and `field`,
#' not `weight.list` and `f`.
#' \code{\link{as.basin.complex}} reconstructs compatible archived results
#' from their original graph and an available field; it is not a guarantee of
#' identical historical output.
#'
#' @section Package ownership:
#' Generic graph construction and algorithms belong to `dgraphs`; use its
#' current public help and object schemas. Examples include
#' `dgraphs::create.path.graph()` and `dgraphs::detect.graph.endpoints()`.
#' Adaptive extrema use \code{\link{detect.adaptive.extrema}} in gflow;
#' `dgraphs::detect.local.extrema()` is a distinct fixed-radius method.
#' Both packages share the `dgraphs::vertices()` generic.
#' Layout packages such as `grip` own their layouts. Generic 3D rendering uses
#' `ivue::plot3D.plain()`, `ivue::plot3D.cont()`, `ivue::plot3D.groups()`,
#' `ivue::layer3D.path()`, and `ivue::layer3D.edges()`.
#'
#' @section Removed estimators:
#' `gflowx` preserves archived response-smoothing and conditional-expectation
#' estimators. Current `geosmooth` is a separate package, not an automatic
#' replacement for all retired estimators. The one-dimensional `magelo`,
#' `mabilo`, `fit.pwlm` and related families were extracted to `malo`; no gflow
#' forwarders remain. Check the destination package's public help.
#' PHATE, diffusion/potential pseudotime, quadratic-form geodesics and
#' random-walk smoothing experiments are removed; no supported successor is
#' claimed here. Generic weighted-p-value helpers and `fassoc*` tests have no
#' public core adapter.
#'
#' @section Association:
#' Use \code{\link{lcor}}, \code{\link{lslope}}, and
#' \code{\link{lslope.neighborhood}} for graph-local estimates.
#' \code{\link{lcor.with.posterior}} summarizes supplied field draws and fits
#' no posterior model. \code{\link{gfcor}} and
#' \code{\link{gfassoc.membership}} accept successful canonical trajectory
#' complexes with both directions and an explicit support stage, as well as
#' archived `basins_of_attraction` inputs. Other canonical method families are
#' rejected; these descriptive summaries do not replace permutation tests.
#'
#' @seealso [create.basin.complex()], [as.basin.complex()], [gflow-package]
#' @name gflow-migration
NULL
