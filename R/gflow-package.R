#' gflow: Geometric Data Analysis Through Gradient Flow
#'
#' @details
#' `gflow` constructs basin and gradient-flow complexes and provides methods for
#' exploring their cells, extrema, trajectories, and flow-aware associations.
#' Generic graph construction is delegated to `dgraphs`.
#'
#' @section Key Features:
#'
#' \strong{1. Geometric Data Representation}
#' \itemize{
#'   \item Pruned intersection k-nearest neighbor (ikNN) graphs that capture data geometry
#'   \item Simplicial complex construction for modeling higher-order relationships
#' }
#'
#' \strong{2. Gradient Flow Domain Decomposition}
#' \itemize{
#'   \item Automatic identification of critical points (local minima, maxima, and saddles)
#'   \item Morse-Smale complex construction for function analysis
#'   \item Natural partitioning into regions of monotonic behavior
#' }
#'
#' \strong{3. Local and Flow-Aware Analysis}
#' \itemize{
#'   \item Basin-, cell-, and trajectory-level summaries
#'   \item Local and flow-aware association analysis
#'   \item Visualization and diagnostics for complex structure
#' }
#'
#' @section Main Function Categories:
#'
#' \strong{Graph and Simplicial Complex Construction:}
#' \itemize{
#'   \item \code{dgraphs::create.iknn.graphs()} - Builds intersection k-nearest neighbor graphs
#'   \item \code{dgraphs::create.single.iknn.graph()} - Creates mutual k-nearest neighbor graphs
#'   \item \code{dgraphs::create.cmst.graph()} - Constructs a minimal spanning tree completion graph
#'   \item \code{\link{create.nerve.complex}} - Creates a nerve complex for k-nearest neighbor covering
#' }
#'
#' \strong{Methods maintained outside gflow:}
#' \itemize{
#'   \item `gflowx::fit.rdgraph.regression()` and related retired graph-response
#'         smoothers
#'   \item \code{malo::amagelo()} - Adaptive 1D smoother with extrema diagnostics
#'   \item \code{malo::magelo()} - Disk-neighborhood 1D model averaging smoother
#'   \item \code{malo::mabilo()} - Symmetric k-hop 1D model averaging smoother
#'   \item \code{malo::magelog()} - 1D local logistic smoother
#'   \item `geosmooth` for local-polynomial, trend-filtering, harmonic, and
#'         graph-spectral smoothing workflows
#' }
#'
#' \strong{Morse-Smale Complex Analysis:}
#' \itemize{
#'   \item \code{\link{create.basin.cx}} - Constructs basin complex from gradient flow
#'   \item \code{\link{compute.graph.gradient.flow}} - Computes gradient flow trajectories
#'   \item \code{\link{find.critical.points}} - Identifies local minima, maxima, and saddles
#' }
#'
#' \strong{Statistical Inference:}
#' \itemize{
#'   \item \code{\link{lcor}} - Symmetric graph-local association
#'   \item \code{\link{lslope}} - Directed graph-local response
#'   \item \code{\link{gfcor}} - Gradient-flow-aware association
#' }
#'
#' @author Pawel Gajer \email{pgajer@@gmail.com}
#'
#' @references
#' Gajer, P. and Ravel, J. (2025). The Geometry of Machine Learning Models.
#' \emph{arXiv preprint} arXiv:2501.01234.
#' \url{https://arxiv.org/abs/2501.01234}
#'
#' Gerber, S., Rubel, O., Bremer, P. T., Pascucci, V., & Whitaker, R. T. (2013).
#' Morse--smale regression. Journal of Computational and Graphical Statistics,
#' 22(1), 193-214.
#'
#' Chen, Y. C., Genovese, C. R., & Wasserman, L. (2017). Statistical inference
#' using the Morse-Smale complex.
#'
#' Gerber, S., & Potter, K. (2012). Data analysis with the morse-smale complex:
#' The msr package for r. Journal of Statistical Software, 50, 1-22.
#'
#' Farrelly, C. M. (2017). Extensions of morse-smale regression with application
#' to actuarial science. arXiv preprint arXiv:1708.05712.
#'
#' Gerber, S., Bremer, P. T., Pascucci, V., & Whitaker, R. (2010). Visual
#' exploration of high dimensional scalar functions. IEEE transactions on
#' visualization and computer graphics, 16(6), 1271-1280.
#'
#' @import stats
#' @import graphics
#' @import grDevices
"_PACKAGE"
