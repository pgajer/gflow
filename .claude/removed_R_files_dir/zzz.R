# Package imports to resolve R CMD check NOTEs

#' @importFrom grDevices colorRampPalette contourLines heat.colors rgb
#' @importFrom graphics box contour hist image layout persp points rect segments strwidth text
#' @importFrom methods hasArg
#' @importFrom stats approxfun density deriv dist dnorm line lm loess mad median na.omit optimize pnorm quantile rbeta setNames splinefun weighted.mean
#' @importFrom utils head
NULL

# Global variables used in the package
# This section declares global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  # Variables used in various functions
  "x",
  "y", 
  "n",
  "L",
  "R",
  
  # External functions from packages
  "Matrix",
  "clusterExport",
  "detectCores",
  "get.knnx",
  "makeCluster",
  "parLapply",
  "persp3d",
  "stopCluster",
  "wasserstein1d",
  
  # Function names that are dynamically called or missing
  "C_wasserstein_distance_1D",
  "IW.kNN.graph",
  "compute.degree.js.divergence",
  "compute.edit_distancesy",
  "compute_euclidean_distances",
  "compute_shortest_paths",
  "create.ED.grid.boxes",
  "discretize",
  "evaluate.function.on.grid",
  "find.critical.points",
  "mabilo",
  "magelo",
  "mstree",
  "mutinformation",
  "plot.function.contours",
  "plot.graph.3d",
  "plot.grid.critical.points",
  "prof.fn",
  "rllm.1D",
  "rm.duplicate.grid.elements"
))