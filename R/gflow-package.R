#' gflow: Geometric Data Analysis Through Gradient Flow
#'
#' @description
#' The gflow package implements geometric methods for analyzing high-dimensional
#' data that naturally inhabits lower-dimensional manifolds. It provides tools
#' for discovering intrinsic geometric structures, building adaptive
#' nearest-neighbor graphs, performing graph-based smoothing, and partitioning
#' data domains through Morse-Smale complexes.
#'
#' @details
#' Modern datasets, particularly in biological sciences, often contain thousands
#' of features with complex non-linear relationships and multi-way interactions
#' that cannot be captured by examining only pairwise associations. Traditional
#' regression models struggle with this complexity, while machine learning
#' approaches, though powerful for prediction, sacrifice interpretability.
#'
#' The gflow package addresses this challenge by exploiting a fundamental
#' observation: despite existing in high-dimensional spaces, real-world data
#' typically lives on much lower-dimensional geometric structures. Like a
#' twisted ribbon in 3D space that is fundamentally 2-dimensional, biological
#' systems are highly constrained, and their data traces out specific geometric
#' shapes within the ambient high-dimensional space.
#'
#' @section Key Features:
#'
#' \strong{1. Geometric Model Construction}
#' \itemize{
#'   \item Adaptive k-nearest neighbor (ikNN) graphs with intelligent edge pruning
#'   \item Intersection-weighted graphs capturing local data geometry
#'   \item Simplicial complex construction for higher-order relationships
#' }
#'
#' \strong{2. Conditional Expectation Estimation on Geometric Structures}
#' \itemize{
#'   \item Graph-based smoothing using harmonic functions
#'   \item Adaptive local polynomial regression (MAGELO family)
#'   \item Spectral methods and diffusion-based smoothing
#' }
#'
#' \strong{3. Morse-Smale Complex Analysis}
#' \itemize{
#'   \item Automatic identification of critical points (minima, maxima, saddles)
#'   \item Basin construction through gradient flow
#'   \item Domain partitioning based on function behavior
#' }
#'
#' \strong{4. Statistical Analysis Within Geometric Regions}
#' \itemize{
#'   \item Within-basin regression and classification
#'   \item Feature importance assessment in local regions
#'   \item Interpretable decomposition of complex models
#' }
#'
#' @section Main Function Categories:
#'
#' \strong{Graph and Simplicial Complex Construction:}
#' \itemize{
#'   \item \code{\link{create.iknn.graphs}} - Builds intersection k-nearest neighbor graphs
#'   \item \code{\link{create.single.iknn.graph}} - Creates mutual k-nearest neighbor graphs
#'   \item \code{\link{create.cmst.graph}} - Constructs a minimal spanning tree completion graph
#'   \item \code{\link{create.nerve.complex}} - Creates a nerve complex associated with a k-nearest neighbor covering of a dataset
#' }
#'
#' \strong{Smoothing and Estimation:}
#' 
#' \emph{Adaptive Local Polynomial Methods:}
#' \itemize{
#'   \item \code{\link{amagelo}} - Adaptive Multi-resolution Averaging with Geometric Local 
#'         weighting. Performs 1D nonparametric smoothing using model-averaged local linear 
#'         regression with automatic bandwidth selection and local extrema identification.
#'   \item \code{\link{magelo}} - Model Averaged Geometric Local regression. Uses disk-shaped 
#'         neighborhoods and averages local polynomial models fitted at uniformly spaced grid 
#'         points, providing smoother transitions in regions with varying data density.
#'   \item \code{\link{agemalo}} - Adaptive Geodesic Model Averaged Local. Creates hierarchies 
#'         of local geodesics through maximal packing vertices and fits weighted linear models 
#'         along geodesic paths with model averaging.
#'   \item \code{\link{pgmalo}} - Piecewise Geometric Model Averaged Local. Combines piecewise 
#'         segmentation with local polynomial regression for handling data with structural breaks.
#' }
#' 
#' \emph{Graph-based LOWESS Methods:}
#' \itemize{
#'   \item \code{\link{deg0.lowess.graph.smoothing}} - Degree-0 LOWESS on graphs. Performs 
#'         locally weighted averaging using graph distances with adaptive bandwidth selection.
#'   \item \code{\link{graph.spectral.lowess}} - Spectral LOWESS. Uses spectral embedding to 
#'         transform graph distances into Euclidean space for local linear regression with 
#'         leave-one-out cross-validation.
#'   \item \code{\link{spectral.lowess.graph.smoothing}} - Alternative spectral LOWESS 
#'         implementation with Cleveland's robustness iterations and multiple kernel options.
#' }
#' 
#' \emph{Kernel and Diffusion Methods:}
#' \itemize{
#'   \item \code{\link{graph.kernel.smoother}} - Graph kernel smoothing with spatially-aware 
#'         cross-validation using buffer zones to prevent spatial autocorrelation bias.
#'   \item \code{\link{graph.diffusion.smoother}} - Iterative diffusion smoothing. Denoises 
#'         graph signals and interpolates missing values while preserving graph structure.
#'   \item \code{\link{harmonic.smoother}} - Harmonic function smoothing. Solves discrete 
#'         Laplace equation by iteratively updating interior vertices as weighted averages 
#'         while keeping boundary values fixed.
#'   \item \code{\link{meanshift.data.smoother}} - Mean shift algorithm on graphs. Non-parametric
#'         technique for finding density maxima with adaptive step sizes and momentum updates.
#' }
#' 
#' \emph{Specialized Methods:}
#' \itemize{
#'   \item \code{\link{amagelogit}} - Model-averaged logistic regression for binary outcomes 
#'         using local polynomial fitting with automatic bandwidth selection.
#'   \item \code{\link{uggmalo}} - Uniform Grid Geometric Model Averaged Local. Combines graph 
#'         structure preservation with local linear modeling on uniform grids.
#'   \item \code{\link{adaptive.uggmalo}} - Adaptive version of UGGMALO with cross-validated 
#'         bandwidth selection and optional bootstrap confidence intervals.
#'   \item \code{\link{graph.spectral.filter}} - Spectral filtering on graphs using eigenvalue 
#'         decomposition for frequency-based smoothing.
#'   \item \code{\link{nerve.cx.spectral.filter}} - Spectral filtering on nerve complexes for 
#'         higher-order topological smoothing.
#'   \item \code{\link{mabilo}} - Model Averaged Bilateral Local. Extends LOWESS with model 
#'         averaging and Bayesian bootstrap for uncertainty quantification.
#' }
#'
#' \strong{Morse-Smale Complex:}
#' \itemize{
#'   \item \code{\link{create.basin.cx}} - Construct basin complex
#'   \item \code{\link{compute.graph.gradient.flow}} - Compute gradient flow trajectories
#'   \item \code{\link{find.critical.points}} - Identify critical points
#' }
#'
#' \strong{Statistical Analysis:}
#' \itemize{
#'   \item \code{\link{fassoc.test}} - Functional association testing
#'   \item \code{\link{compute.bayesian.effects}} - Bayesian effect estimation
#'   \item \code{\link{wasserstein.distance}} - Geometric distance metrics
#' }
#'
#' @section Theoretical Foundation:
#' The theoretical foundation is described in:
#' Gajer, P. and Ravel, J. (2025). "The Geometry of Machine Learning Models".
#' arXiv preprint arXiv:2501.01234. Available at \url{https://arxiv.org/abs/2501.01234}
#'
#' @section Getting Started:
#' For a comprehensive introduction, see the package vignettes:
#' \itemize{
#'   \item \code{vignette("gflow-overview", package = "gflow")} - Package overview
#'   \item \code{vignette("iknn-graphs", package = "gflow")} - Graph construction tutorial
#'   \item \code{vignette("morse-smale", package = "gflow")} - Morse-Smale analysis
#' }
#'
#' @author
#' Pawel Gajer \email{pgajer@@gmail.com}
#'
#' Maintainer: Pawel Gajer \email{pgajer@@gmail.com}
#'
#' @references
#' Gajer, P. and Ravel, J. (2025). The Geometry of Machine Learning Models.
#' \emph{arXiv preprint} arXiv:2501.01234.
#'
#' @keywords package
#' @docType package
#' @name gflow-package
#' @aliases gflow
NULL
