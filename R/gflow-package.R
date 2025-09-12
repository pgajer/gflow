#' gflow: Geometric Data Analysis Through Gradient Flow
#'
#' @details
#' Modern datasets, particularly in biological sciences, often contain
#' thousands of features with complex non-linear relationships and
#' multi-way interactions that cannot be captured by examining only
#' pairwise associations. Traditional regression models struggle with
#' this complexity, while machine learning approaches, though powerful
#' for prediction, sacrifice interpretability.
#'
#' The gflow package addresses this challenge by exploiting a fundamental
#' observation: despite existing in high-dimensional spaces, real-world
#' data typically lives on much lower-dimensional geometric structures.
#' Like a twisted ribbon in 3D space that is fundamentally 2-dimensional,
#' most real-world systems generate highly constrained data that traces
#' out specific geometric shapes within the ambient space. The geometry
#' of these shapes implicitly encodes complex associations between features.
#'
#' Rather than working in a fixed coordinate system, gflow adopts a
#' coordinate-free framework that models data as a weighted graph (or more
#' generally, a Riemannian simplicial complex) capturing the underlying
#' geometric structure. This transforms the inference problem from \eqn{E[Y|X]}
#' to \eqn{E[Y|G(X)]}, where \eqn{G(X)} represents the geometric object
#' constructed from your data matrix \eqn{X}. This reformulation provides a more
#' flexible representation that adapts to the inherent complexity of the data
#' without imposing predetermined functional forms.
#'
#' A key strength of gflow is its ability to provide interpretable results from
#' complex data. The gradient-flow decomposition identifies natural regions
#' where relationships between predictors and outcomes are locally simple and
#' monotonic. This allows one to "peer inside" black-box models and understand
#' how predictions change across different parts of the data space. By
#' respecting the intrinsic geometry of the data, gflow bridges the gap between
#' the flexibility of modern machine learning and the interpretability of
#' classical statistical methods.
#'
#' @section Key Features:
#'
#' \strong{1. Geometric Data Representation}
#' \itemize{
#'   \item Pruned intersection k-nearest neighbor (ikNN) graphs that capture data geometry
#'   \item Simplicial complex construction for modeling higher-order relationships
#' }
#'
#' \strong{2. Robust Signal Recovery and Smoothing}
#' \itemize{
#'   \item Kernel graph Laplacian smoothing for noise reduction
#'   \item Spectral filtering and diffusion-based approaches
#' }
#'
#' \strong{3. Gradient Flow Domain Decomposition}
#' \itemize{
#'   \item Automatic identification of critical points (local minima, maxima, and saddles)
#'   \item Morse-Smale complex construction for function analysis
#'   \item Natural partitioning into regions of monotonic behavior
#' }
#'
#' \strong{4. Interpretable Statistical Inference}
#' \itemize{
#'   \item Within-region regression and classification
#'   \item Bootstrap-Wasserstein testing for non-linear associations
#'   \item Feature importance assessment across different data regions
#'   \item Visualization tools for understanding model behavior
#' }
#'
#' @section Main Function Categories:
#'
#' \strong{Graph and Simplicial Complex Construction:}
#' \itemize{
#'   \item \code{\link{create.iknn.graphs}} - Builds intersection k-nearest neighbor graphs
#'   \item \code{\link{create.single.iknn.graph}} - Creates mutual k-nearest neighbor graphs
#'   \item \code{\link{create.cmst.graph}} - Constructs a minimal spanning tree completion graph
#'   \item \code{\link{create.nerve.complex}} - Creates a nerve complex for k-nearest neighbor covering
#' }
#'
#' \strong{Conditional Expectation Estimation Methods:}
#'
#' \emph{Model-Averaged Local Regression (1D):}
#' \itemize{
#'   \item \code{\link{amagelo}} - Adaptive MAGELO with automatic bandwidth selection, extrema
#'         detection, and robust fitting for continuous responses
#'   \item \code{\link{magelo}} - Disk-based neighborhoods with grid-centered models for smooth
#'         transitions in varying density data
#'   \item \code{\link{mabilo}} - Symmetric k-hop neighborhoods with model averaging and
#'         Bayesian bootstrap for uncertainty quantification
#'   \item \code{\link{amagelogit}} - Grid-based logistic regression with model averaging for
#'         binary outcomes
#'   \item \code{\link{maelog}} - Data-centered logistic regression with adaptive disk
#'         neighborhoods and k-NN fallback
#' }
#'
#' \emph{Graph-Based Local Regression Methods:}
#' \itemize{
#'   \item \code{\link{deg0.lowess.graph.smoothing}} - Locally weighted averaging using graph
#'         distances with adaptive bandwidth selection
#'   \item \code{\link{graph.spectral.lowess}} - Spectral embedding transforms graph distances
#'         to Euclidean space for local linear regression
#'   \item \code{\link{spectral.lowess.graph.smoothing}} - Includes Cleveland's robustness
#'         iterations and multiple kernel options
#' }
#'
#' \emph{Diffusion and Kernel Methods:}
#' \itemize{
#'   \item \code{\link{graph.kernel.smoother}} - Spatially-aware cross-validation with buffer
#'         zones to prevent autocorrelation bias
#'   \item \code{\link{graph.diffusion.smoother}} - Iterative diffusion for denoising and
#'         interpolating missing values
#'   \item \code{\link{harmonic.smoother}} - Solves discrete Laplace equation for smooth
#'         interpolation with fixed boundary values
#'   \item \code{\link{graph.spectral.filter}} - Frequency-based smoothing using Laplacian
#'         eigendecomposition
#' }
#'
#' \emph{Path and Geodesic Methods:}
#' \itemize{
#'   \item \code{\link{pgmalo}} - Piecewise segmentation with automatic change point detection
#'         for data with structural breaks
#'   \item \code{\link{agemalo}} - Adaptive geodesic regression using hierarchies of local
#'         geodesics through maximal packing vertices
#'   \item \code{\link{adaptive.uggmalo}} - Uniform grid approach with path-based local models
#'         and cross-validated bandwidth selection
#' }
#'
#' \emph{Specialized Methods:}
#' \itemize{
#'   \item \code{\link{magelog}}, \code{\link{mabilog}}, \code{\link{uggmalog}} - Log-space
#'         variants for positive-valued responses
#'   \item \code{\link{ulogit}} - Uniform grid logistic regression for binary outcomes on graphs
#'   \item \code{\link{nerve.cx.spectral.filter}} - Spectral filtering over simplicial complexes
#'   \item \code{\link{meanshift.data.smoother}} - Data denoising using mean-shift algorithm with adaptive step sizes
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
#'   \item \code{\link{fassoc.test}} - Tests for non-linear functional associations
#'   \item \code{\link{compute.bayesian.effects}} - Bayesian effect size estimation
#'   \item \code{\link{wasserstein.distance}} - Computes optimal transport distances
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
