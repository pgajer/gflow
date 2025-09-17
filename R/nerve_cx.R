#' Creates a Nerve Complex Associated With kNN Covering of a Dataset
#'
#' This function creates a nerve complex from a set of points in \eqn{R^n} based on their
#' k-nearest neighbor covering.
#'
#' @param X Matrix of point coordinates, where rows are points and columns are dimensions
#' @param k Number of nearest neighbors to use for the covering
#' @param max.dim Maximum dimension of simplices to compute
#'
#' @return A nerve complex object
#' @export
create.nerve.complex <- function(X, k, max.dim = 2) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error"))
            stop("X must be a matrix or coercible to a numeric matrix.")
    }
    if (!is.numeric(X)) stop("X must be numeric.")
    if (any(!is.finite(X))) stop("X cannot contain NA/NaN/Inf.")
    storage.mode(X) <- "double"

    if (k < 2) {
        stop("k must be at least 2")
    }

    if (max.dim < 1) {
        stop("max.dim must be at least 1")
    }

    result <- .Call(S_create_nerve_complex,
                    X,
                    as.integer(k + 1), # Note that ANN library is configured so that it includes the query point within the set of kNN's so we need to increase k by 1 to really get kNN's
                    as.integer(max.dim))

    class(result) <- "nerve_complex"

    return(result)
}

#' Set Function Values on Nerve Complex Vertices
#'
#' This function sets function values at the vertices of the nerve complex.
#'
#' @param complex A nerve complex object
#' @param values Vector of function values (one per vertex)
#'
#' @return The updated nerve complex object (invisibly)
#' @export
set.complex.function.values <- function(complex, values) {
  # Check inputs
  if (!inherits(complex, "nerve_complex")) {
    stop("complex must be a nerve_complex object")
  }

  if (length(values) != complex$n_vertices) {
    stop("Number of values must match number of vertices")
  }

  # Call C++ function
    .Call(S_set_function_values,
          complex$complex_ptr,
          as.numeric(values))

  return(invisible(complex))
}

#' Set Weight Scheme for Nerve Complex
#'
#' This function sets the weighting scheme for simplices in the nerve complex.
#'
#' @param complex A nerve complex object
#' @param weight.type Type of weight scheme to use
#' @param params Parameters for the weight scheme (if needed)
#'
#' @details Available weight types:
#' \itemize{
#'   \item "uniform": All simplices have weight 1.0
#'   \item "inverse_distance": Weights are inversely proportional to average squared distance
#'   \item "gaussian": Weights use Gaussian kernel based on distances (param: sigma)
#'   \item "volume": Weights based on approximate simplex volume (param: alpha)
#'   \item "gradient": Weights based on function gradient over simplex (param: gamma)
#' }
#'
#' @return The updated nerve complex object (invisibly)
#' @export
set.complex.weight.scheme <- function(complex, weight.type, params = numeric(0)) {
  # Check inputs
  if (!inherits(complex, "nerve_complex")) {
    stop("complex must be a nerve_complex object")
  }

  valid_types <- c("uniform", "inverse_distance", "gaussian", "volume", "gradient")
  if (!weight.type %in% valid_types) {
    stop(paste0("weight.type must be one of: ", paste(valid_types, collapse = ", ")))
  }

  # Check params based on weight type
  if (weight.type == "gaussian" && length(params) < 1) {
    stop("Gaussian weight scheme requires sigma parameter")
  }

  if (weight.type == "volume" && length(params) < 1) {
    stop("Volume weight scheme requires alpha parameter")
  }

  if (weight.type == "gradient" && length(params) < 1) {
    stop("Gradient weight scheme requires gamma parameter")
  }

  # Call C++ function
    .Call(S_set_weight_scheme,
          complex$complex_ptr,
          weight.type,
          as.numeric(params))

  return(invisible(complex))
}

#' Solve Full Laplacian Problem on Nerve Complex
#'
#' This function solves the full Laplacian problem on the nerve complex to
#' extend the function values.
#'
#' @param complex A nerve complex object
#' @param lambda Regularization parameter
#' @param dim.weights Weights for each dimension's contribution
#'
#' @return Vector of extended function values
#' @export
complex.laplacian.solve <- function(complex, lambda = 1.0,
                                   dim.weights = rep(1, complex$max_dimension + 1)) {
  # Check inputs
  if (!inherits(complex, "nerve_complex")) {
    stop("complex must be a nerve_complex object")
  }

  if (length(dim.weights) < complex$max_dimension + 1) {
    stop("Not enough dimension weights provided")
  }

  # Call C++ function
    result <- .Call(S_solve_full_laplacian,
                    complex$complex_ptr,
                    as.numeric(lambda),
                    as.numeric(dim.weights))

  return(result)
}

#' Get Simplex Counts for Nerve Complex
#'
#' This function returns the number of simplices in each dimension.
#'
#' @param complex A nerve complex object
#'
#' @return Vector of simplex counts
#' @export
get.simplex.counts <- function(complex) {
                                        # Check inputs
    if (!inherits(complex, "nerve_complex")) {
        stop("complex must be a nerve_complex object")
    }

    ## Call C++ function
    .Call(S_get_simplex_counts,
          complex$complex_ptr)
}

#' Extract 1-Skeleton Graph from Nerve Complex
#'
#' This function extracts the 1-skeleton of the nerve complex as a graph.
#'
#' @param complex A nerve complex object
#'
#' @return A graph representation of the 1-skeleton
#' @export
extract.skeleton.graph <- function(complex) {
                                        # Check inputs
    if (!inherits(complex, "nerve_complex")) {
        stop("complex must be a nerve_complex object")
    }

    ## Call C++ function
    result <- .Call(S_extract_skeleton_graph,
                    complex$complex_ptr)

    return(result)
}

#' Print Method for Nerve Complex Objects
#'
#' @param x A nerve complex object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the object
#' @export
print.nerve_complex <- function(x, ...) {
  cat("Nerve Complex:\n")
  cat("  Number of vertices:", x$n_vertices, "\n")
  cat("  Maximum dimension:", x$max_dimension, "\n")
  cat("  Simplex counts by dimension:\n")

  for (d in 0:x$max_dimension) {
    cat(sprintf("    %d-simplices: %d\n", d, x$simplex_counts[d+1]))
  }

  return(invisible(x))
}

#' Compare Graph vs Simplicial Complex Regression
#'
#' Compares regression performance using graph Laplacian vs full simplicial complex Laplacian.
#'
#' @param n.points Number of points to generate
#' @param k Number of nearest neighbors
#' @param max.dim Maximum simplex dimension
#' @param lambda Regularization parameter
#' @param n.test Number of test points
#' @param weight.type Weight scheme to use
#' @param weight.params Parameters for weight scheme
#' @param dim.weights Weights for each dimension's contribution
#'
#' @return A list containing comparison results
#' @export
graph.vs.complex.regression.compare <- function(n.points = 200,
                                                k = 5,
                                                max.dim = 2,
                                                lambda = 1.0,
                                                n.test = 50,
                                                weight.type = "gaussian",
                                                weight.params = c(0.2),
                                                dim.weights = c(1, 0.5, 0.25)) {
  # Generate random points in [0,1]²
  set.seed(123)
  train.points <- matrix(runif(n.points * 2), ncol = 2)
  test_points <- matrix(runif(n.test * 2), ncol = 2)

  # Create target function (mixture of Gaussians)
  target_fn <- get.gaussian.mixture(n.components = 5, sd.knot = 0.1)

  # Evaluate function at points
  train_values <- target_fn(train.points)
  test_true_values <- target_fn(test_points)

  # Create nerve complex
  complex <- create.nerve.complex(train.points, k, max.dim)
  complex <- set.complex.function.values(complex, train_values)
  complex <- set.complex.weight.scheme(complex, weight.type, weight.params)

  # Solve using only graph Laplacian (1-skeleton)
  graph_only_weights <- c(1, rep(0, max.dim))
  graph_pred <- complex.laplacian.solve(complex, lambda, graph_only_weights)

  # Solve using full complex
  complex_pred <- complex.laplacian.solve(complex, lambda, dim.weights)

  # Perform nearest neighbor interpolation to test points
  nn_interp <- function(pred_values, train_pts, test_pts) {
    result <- numeric(nrow(test_pts))

    for (i in 1:nrow(test_pts)) {
      # Find nearest neighbor
      dists <- apply(train_pts, 1, function(p) sum((p - test_pts[i,])^2))
      nn_idx <- which.min(dists)
      result[i] <- pred_values[nn_idx]
    }

    return(result)
  }

  graph_test_pred <- nn_interp(graph_pred, train.points, test_points)
  complex_test_pred <- nn_interp(complex_pred, train.points, test_points)

  # Calculate MSE
  graph_mse <- mean((graph_test_pred - test_true_values)^2)
  complex_mse <- mean((complex_test_pred - test_true_values)^2)

  # Return results
  return(list(
    graph_mse = graph_mse,
    complex_mse = complex_mse,
    improvement_pct = 100 * (graph_mse - complex_mse) / graph_mse,
    train.points = train.points,
    test_points = test_points,
    train_values = train_values,
    test_true_values = test_true_values,
    graph_pred = graph_pred,
    complex_pred = complex_pred,
    graph_test_pred = graph_test_pred,
    complex_test_pred = complex_test_pred,
    complex_info = complex
  ))
}
