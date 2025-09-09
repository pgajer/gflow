##
## Collection of functions for estimating the Mean Absolute Bias (MAB) for
## different types of regression models on different types of multi-dimensional
## synthetic data
##

#' Create k-fold cross-validation folds
#'
#' This function creates indices for k-fold cross-validation, matching the behavior
#' of \code{\link[caret]{createFolds}} but without the dependency.
#'
#' @param y A vector of outcomes
#' @param k The number of folds (default: 10)
#' @param list Logical. If TRUE, returns a list of fold indices. If FALSE, returns
#'   a vector of fold assignments (default: TRUE)
#' @param returnTrain Logical. If TRUE, returns training set indices for each fold.
#'   If FALSE, returns test set indices (default: FALSE)
#'
#' @return If list = TRUE, a list with k elements where each element contains the
#'   indices for that fold. If list = FALSE, a vector of integers indicating fold
#'   membership for each observation.
#'
#' @examples
#' # Create 5 folds for a vector of 100 observations
#' y <- rnorm(100)
#' folds <- create.folds(y, k = 5, list = TRUE, returnTrain = TRUE)
#'
#' @noRd
create.folds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {
  # Handle edge cases
  if (length(y) < k) {
    stop("k must be less than or equal to the number of observations")
  }

  n <- length(y)

  # For classification, try to maintain class balance
  if (is.factor(y) || length(unique(y)) < min(30, n/2)) {
    # Stratified sampling for classification
    folds <- vector("list", k)

    # Get unique classes and their indices
    classes <- split(seq_along(y), y)

    # Assign each class's observations to folds
    for (cls in classes) {
      if (length(cls) >= k) {
        # Randomly shuffle indices within class
        cls_shuffled <- sample(cls)
        # Assign to folds in round-robin fashion
        fold_assignments <- rep(1:k, length.out = length(cls_shuffled))
        for (i in 1:k) {
          folds[[i]] <- c(folds[[i]], cls_shuffled[fold_assignments == i])
        }
      } else {
        # If fewer observations than folds, assign randomly
        fold_assignments <- sample(1:k, length(cls), replace = FALSE)
        for (i in seq_along(cls)) {
          folds[[fold_assignments[i]]] <- c(folds[[fold_assignments[i]]], cls[i])
        }
      }
    }
  } else {
    # Regular k-fold for regression
    indices <- sample(seq_len(n))
    folds <- split(indices, cut(seq_along(indices), k, labels = FALSE))
  }

  # Convert to returnTrain format if requested
  if (returnTrain) {
    folds <- lapply(seq_along(folds), function(i) {
      setdiff(seq_len(n), folds[[i]])
    })
  }

  # Return as vector if list = FALSE
  if (!list) {
    fold_vec <- integer(n)
    for (i in seq_along(folds)) {
      fold_vec[folds[[i]]] <- i
    }
    return(fold_vec)
  }

  return(folds)
}

#' Estimate Mean Absolute Bias (MAB) of MAGELO Model
#'
#' Estimates the Mean Absolute Bias (MAB) of a MAGELO model.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values.
#' @param xt A numeric vector of true predictor values for evaluation.
#' @param yt A numeric vector of true response values for evaluation.
#' @param deg The degree of the MAGELO polynomial model (default = 2).
#' @param y.binary Logical indicating if y is a binary variable (default = FALSE).
#' @param n.cores Number of cores to use for parallel computation (default = 10).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predicted values at xt}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the degree used}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- seq(0, 10, length.out = 50)
#' yt <- sin(xt)
#'
#' # Fit magelo model and compute MAB
#' result <- get.magelo.MAB(x, y, xt, yt, deg = 2)
#' print(result$MAB)
#' }
#'
#' @export
get.magelo.MAB <- function(x, y, xt, yt, deg = 2, y.binary = FALSE, n.cores = 10) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    if (!is.numeric(deg) || length(deg) != 1 || deg < 1 || deg != round(deg)) {
        stop("'deg' must be a positive integer")
    }

    if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1 || n.cores != round(n.cores)) {
        stop("'n.cores' must be a positive integer")
    }

    # Fit MAGELO model
    magelo.m <- magelo(x, y, degree = deg, y.binary = y.binary, n.BB = 0,
                       min.bw.f = 0.01, n.bws = 100, n.cores = n.cores)

    # Generate predictions
    magelo.pred <- predict(magelo.m, newdata = xt)$Ey

    # Calculate residuals and metrics
    residuals <- abs(magelo.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((magelo.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = magelo.pred,
        residuals = residuals,
        parameters = list(deg = deg)
    )
}

#' Estimate Mean Absolute Bias (MAB) of LOESS Model
#'
#' Estimates the Mean Absolute Bias (MAB) of a LOESS (Locally Estimated Scatterplot Smoothing)
#' model with the span parameter chosen based on cross-validated MAE estimates.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values.
#' @param xt A numeric vector of true predictor values for evaluation.
#' @param yt A numeric vector of true response values for evaluation.
#' @param folds A list of indices for cross-validation folds
#' @param deg The degree of the local polynomial used (default = 2).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predicted values at xt}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the optimal span}
#' }
#'
#' @importFrom stats loess predict
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- seq(0, 10, length.out = 50)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Fit LOESS model and compute MAB
#' result <- get.loess.MAB(x, y, xt, yt, folds, deg = 2)
#' print(result$MAB)
#' }
#'
#' @export
get.loess.MAB <- function(x, y, xt, yt, folds, deg = 2) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.numeric(deg) || length(deg) != 1 || deg < 1 || deg > 2 || deg != round(deg)) {
        stop("'deg' must be 1 or 2")
    }

    # Internal function to compute MAE for a given span
    loess.MAE.fn <- function(x, y, span, folds) {
        perform.loess.cv <- function(training.set, test.set, span) {
            tryCatch({
                loess.fit <- loess(y ~ x, data = training.set, span = span, degree = deg)
                predictions <- predict(loess.fit, newdata = test.set)
                mae <- mean(abs(predictions - test.set$y), na.rm = TRUE)
                return(mae)
            }, error = function(e) {
                return(NA)
            })
        }

        mae.values <- numeric(length(folds))
        for(i in seq_along(folds)) {
            train.indices <- folds[[i]]
            test.indices <- setdiff(seq_along(x), train.indices)

            if (length(test.indices) == 0) {
                stop("Cross-validation fold resulted in empty test set")
            }

            training.set <- data.frame(x = x[train.indices], y = y[train.indices])
            test.set <- data.frame(x = x[test.indices], y = y[test.indices])
            mae.values[i] <- perform.loess.cv(training.set, test.set, span)
        }

        mean.mae <- mean(mae.values, na.rm = TRUE)
        return(mean.mae)
    }

    # Define span values to test
    spans <- seq(0.1, 1, by = 0.01)

    # Find optimal span value
    loess.MAE <- numeric(length(spans))
    for (i in seq_along(spans)) {
        span <- spans[i]
        loess.MAE[i] <- loess.MAE.fn(x, y, span, folds)
    }

    # Remove NA values if any
    valid.indices <- !is.na(loess.MAE)
    if (sum(valid.indices) == 0) {
        stop("All span values resulted in errors during cross-validation")
    }

    spans <- spans[valid.indices]
    loess.MAE <- loess.MAE[valid.indices]

    # Select optimal span
    opt.span <- spans[which.min(loess.MAE)]

    # Fit the final model using the optimal span
    loess.m <- loess(y ~ x, span = opt.span, degree = deg)

    # Predict using the loess model
    loess.pred <- predict(loess.m, newdata = data.frame(x = xt))

    # Handle potential NA predictions
    if (any(is.na(loess.pred))) {
        warning("Some predictions are NA. These will be excluded from MAB calculation.")
    }

    # Compute metrics
    residuals <- abs(loess.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((loess.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = loess.pred,
        residuals = residuals,
        parameters = list(opt.span = opt.span)
    )
}

#' Estimate Mean Absolute Bias (MAB) of Local Polynomial Regression
#'
#' Fits a local polynomial regression model to the given dataset and evaluates its
#' performance using cross-validation. The mean absolute error (MAE) is used as the
#' performance metric, and the function selects the optimal bandwidth for the local
#' polynomial regression.
#'
#' @param x A numeric vector representing the predictor variable for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A numeric vector representing the predictor variable for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param folds A list of indices for cross-validation folds
#' @param bandwidths A numeric vector of bandwidth values to be tested.
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the local polynomial regression model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the optimal bandwidth}
#' }
#'
#' @details The function performs cross-validation to find the optimal bandwidth that
#' minimizes the MAE. It uses the \code{locpoly} function from the KernSmooth package
#' for fitting local polynomial models. The optimal bandwidth is then used to fit the
#' final model and make predictions on the test set.
#'
#' @importFrom stats approx
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- seq(0, 10, length.out = 50)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Define bandwidth values
#' bandwidths <- seq(0.1, 1, by = 0.1)
#'
#' # Run the function
#' result <- get.locpoly.MAB(x, y, xt, yt, folds, bandwidths)
#' print(result$MAB)
#' }
#'
#' @export
get.locpoly.MAB <- function(x, y, xt, yt, folds, bandwidths) {

    if (!requireNamespace("KernSmooth", quietly = TRUE)) {
        stop("This function requires the suggested package 'KernSmooth'. ",
             "Install it with install.packages('KernSmooth').",
             call. = FALSE)
    }

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.numeric(bandwidths) || length(bandwidths) == 0) {
        stop("'bandwidths' must be a non-empty numeric vector")
    }

    if (any(bandwidths <= 0)) {
        stop("All bandwidth values must be positive")
    }

    # Internal function to compute MAE for a given bandwidth
    locpoly.MAE.fn <- function(x, y, bandwidth, folds) {
        perform.locpoly.cv <- function(training.set, test.set, bandwidth) {
            tryCatch({
                locpoly.fit <- KernSmooth::locpoly(x = training.set$x, y = training.set$y,
                                                   bandwidth = bandwidth)

                # Predict using interpolation
                predictions <- approx(locpoly.fit$x, locpoly.fit$y,
                                    xout = test.set$x, rule = 2)$y

                mae <- mean(abs(predictions - test.set$y), na.rm = TRUE)
                return(mae)
            }, error = function(e) {
                return(NA)
            })
        }

        mae.values <- numeric(length(folds))
        for (i in seq_along(folds)) {
            train.indices <- folds[[i]]
            test.indices <- setdiff(seq_along(x), train.indices)

            if (length(test.indices) == 0) {
                stop("Cross-validation fold resulted in empty test set")
            }

            training.set <- data.frame(x = x[train.indices], y = y[train.indices])
            test.set <- data.frame(x = x[test.indices], y = y[test.indices])
            mae.values[i] <- perform.locpoly.cv(training.set, test.set, bandwidth)
        }

        mean.MAE <- mean(mae.values, na.rm = TRUE)
        return(mean.MAE)
    }

    # Find optimal bandwidth
    mae.scores <- setNames(numeric(length(bandwidths)), as.character(bandwidths))

    for (i in seq_along(bandwidths)) {
        bandwidth <- bandwidths[i]
        mae.scores[i] <- locpoly.MAE.fn(x, y, bandwidth, folds)
    }

    # Remove NA values if any
    valid.indices <- !is.na(mae.scores)
    if (sum(valid.indices) == 0) {
        stop("All bandwidth values resulted in errors during cross-validation")
    }

    mae.scores <- mae.scores[valid.indices]

    # Select optimal bandwidth
    opt.bandwidth <- as.numeric(names(mae.scores)[which.min(mae.scores)])

    # Fit the final model using the optimal bandwidth
    final.locpoly <- KernSmooth::locpoly(x = x, y = y, bandwidth = opt.bandwidth)

    # Predict using the locpoly model
    locpoly.pred <- approx(final.locpoly$x, final.locpoly$y, xout = xt, rule = 2)$y

    # Handle potential NA predictions
    if (any(is.na(locpoly.pred))) {
        warning("Some predictions are NA. These will be excluded from MAB calculation.")
    }

    # Compute metrics
    residuals <- abs(locpoly.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((locpoly.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = locpoly.pred,
        residuals = residuals,
        parameters = list(opt.bandwidth = opt.bandwidth)
    )
}

#' Estimate Mean Absolute Bias (MAB) of Gaussian Process Regression Model
#'
#' Fits a Gaussian Process Regression (GPR) model to a given dataset and evaluates
#' its performance. The Mean Absolute Bias (MAB) is calculated using the test set.
#'
#' @param x A matrix or data frame representing the predictor variables for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A matrix or data frame representing the predictor variables for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param folds A list of indices for cross-validation folds (currently not used in implementation).
#' @param kernel A kernel function to be used in GPR (default = 'rbfdot').
#'        Options include 'rbfdot', 'polydot', 'laplacedot', etc. from kernlab package.
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the GPR model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the kernel used}
#' }
#'
#' @details The function uses Gaussian Process Regression to model the relationship
#' between the predictors and the response. Users can specify different kernel
#' functions to explore various types of non-linear relationships in the data.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' library(kernlab)
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- matrix(rnorm(50), ncol = 1)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Run the function
#' result <- get.gausspr.MAB(x, y, xt, yt, folds)
#' print(result$MAB)
#' }
#'
#' @export
get.gausspr.MAB <- function(x, y, xt, yt, folds, kernel = 'rbfdot') {

    # Input validation
    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    # Convert to matrix if needed
    if (is.vector(x)) {
        x <- as.matrix(x)
    } else if (!is.matrix(x) && !is.data.frame(x)) {
        stop("'x' must be a vector, matrix, or data frame")
    }

    if (is.vector(xt)) {
        xt <- as.matrix(xt)
    } else if (!is.matrix(xt) && !is.data.frame(xt)) {
        stop("'xt' must be a vector, matrix, or data frame")
    }

    # Convert data frames to matrices
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    if (is.data.frame(xt)) {
        xt <- as.matrix(xt)
    }

    if (nrow(x) != length(y)) {
        stop("Number of rows in 'x' must equal length of 'y'")
    }

    if (nrow(xt) != length(yt)) {
        stop("Number of rows in 'xt' must equal length of 'yt'")
    }

    if (ncol(x) != ncol(xt)) {
        stop("'x' and 'xt' must have the same number of columns")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.character(kernel) || length(kernel) != 1) {
        stop("'kernel' must be a single character string")
    }

    # Valid kernels in kernlab
    valid.kernels <- c('rbfdot', 'polydot', 'laplacedot', 'tanhdot',
                      'vanilladot', 'besseldot', 'anovadot', 'splinedot')

    if (!(kernel %in% valid.kernels)) {
        warning(paste("'kernel' may not be valid. Common options are:",
                     paste(valid.kernels, collapse = ", ")))
    }

    # Fit the GPR model on the entire dataset
    tryCatch({
        gausspr.m <- kernlab::gausspr(x = x, y = y, kernel = kernel)
    }, error = function(e) {
        stop(paste("Error fitting Gaussian Process model:", e$message))
    })

    # Generate predictions
    tryCatch({
        gausspr.pred <- kernlab::predict(gausspr.m, newdata = xt)
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Convert predictions to numeric vector if needed
    if (is.matrix(gausspr.pred) && ncol(gausspr.pred) == 1) {
        gausspr.pred <- as.vector(gausspr.pred)
    }

    # Compute metrics
    residuals <- abs(gausspr.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((gausspr.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = gausspr.pred,
        residuals = residuals,
        parameters = list(kernel = kernel)
    )
}

#' Estimate Mean Absolute Bias (MAB) of Kernel Ridge Regression Model
#'
#' Fits a Kernel Ridge Regression (KRR) model to a given dataset and evaluates
#' its performance using cross-validation. The Mean Absolute Bias (MAB) is calculated
#' using the test set with optimized hyperparameters.
#'
#' @param x A matrix or data frame representing the predictor variables for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A matrix or data frame representing the predictor variables for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param kernel A kernel function to be used in KRR (default = 'rbfdot').
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the KRR model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing optimal sigma and lambda values}
#' }
#'
#' @details The function uses Kernel Ridge Regression to model the relationship
#' between the predictors and the response. It performs cross-validation to
#' find optimal hyperparameters (sigma and lambda) and computes the MAB using
#' the true x,y values. Users can specify different kernel functions to explore
#' various types of non-linear relationships in the data.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- matrix(rnorm(50), ncol = 1)
#' yt <- sin(xt)
#'
#' # Run the function
#' result <- get.krr.MAB(x, y, xt, yt)
#' print(result$MAB)
#' }
#'
#' @export
get.krr.MAB <- function(x, y, xt, yt, kernel = 'rbfdot') {

    if (!requireNamespace("CVST", quietly = TRUE)) {
        stop("This function requires the suggested package 'CVST'. ",
             "Install it with install.packages('CVST').",
             call. = FALSE)
    }

    # Input validation
    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    # Convert to matrix if needed
    if (is.vector(x)) {
        x <- as.matrix(x)
    } else if (!is.matrix(x) && !is.data.frame(x)) {
        stop("'x' must be a vector, matrix, or data frame")
    }

    if (is.vector(xt)) {
        xt <- as.matrix(xt)
    } else if (!is.matrix(xt) && !is.data.frame(xt)) {
        stop("'xt' must be a vector, matrix, or data frame")
    }

    # Convert data frames to matrices
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    if (is.data.frame(xt)) {
        xt <- as.matrix(xt)
    }

    if (nrow(x) != length(y)) {
        stop("Number of rows in 'x' must equal length of 'y'")
    }

    if (nrow(xt) != length(yt)) {
        stop("Number of rows in 'xt' must equal length of 'yt'")
    }

    if (ncol(x) != ncol(xt)) {
        stop("'x' and 'xt' must have the same number of columns")
    }

    if (!is.character(kernel) || length(kernel) != 1) {
        stop("'kernel' must be a single character string")
    }

    # Valid kernels for KRR
    valid.kernels <- c('rbfdot', 'polydot', 'laplacedot', 'tanhdot', 'vanilladot')

    if (!(kernel %in% valid.kernels)) {
        warning(paste("'kernel' may not be valid. Common options are:",
                     paste(valid.kernels, collapse = ", ")))
    }

    # Construct KRR learner
    krr <- CVST::constructKRRLearner()

    # Define hyperparameter grid
    lambdas <- 10^(-8:0)
    sigmas <- 10^((1:9)/3)

    # Construct parameters
    tryCatch({
        params <- CVST::constructParams(kernel = kernel, sigma = sigmas, lambda = lambdas)
    }, error = function(e) {
        stop(paste("Error constructing parameters:", e$message))
    })

    # Construct data
    dat <- CVST::constructData(x, y)

    # Perform cross-validation to find optimal parameters
    tryCatch({
        opt <- CVST::CV(dat, krr, params, fold = 10, verbose = FALSE)
    }, error = function(e) {
        stop(paste("Error during cross-validation:", e$message))
    })

    # Extract optimal parameters
    opt.sigma <- opt[[1]]$sigma
    opt.lambda <- opt[[1]]$lambda

    # Set final parameters
    param <- list(kernel = kernel, sigma = opt.sigma, lambda = opt.lambda)

    # Train final model with optimal parameters
    tryCatch({
        krr.model <- krr$learn(dat, param)
    }, error = function(e) {
        stop(paste("Error training KRR model:", e$message))
    })

    # Predict on test data
    dat.tst <- CVST::constructData(xt, rep(0, nrow(xt)))  # Dummy y values for prediction

    tryCatch({
        krr.pred <- krr$predict(krr.model, dat.tst)[, 1]
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Compute metrics
    residuals <- abs(krr.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((krr.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = krr.pred,
        residuals = residuals,
        parameters = list(
            kernel = kernel,
            sigma = opt.sigma,
            lambda = opt.lambda
        )
    )
}

#' Estimate Mean Absolute Bias (MAB) of Support Vector Machine Model
#'
#' Fits a Support Vector Machine (SVM) model to a given dataset and evaluates
#' its performance using cross-validation. The Mean Absolute Bias (MAB) is used
#' as the performance metric. The function selects the optimal hyperparameters
#' for the SVM model based on cross-validated MAE.
#'
#' @param x A matrix, data frame, or vector representing the predictor variables for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A matrix, data frame, or vector representing the predictor variables for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param folds A list of indices for cross-validation folds
#' @param cost A numeric vector of candidate cost values to evaluate (default = 10^seq(-1, 1, length.out = 3)).
#' @param gamma A numeric vector of candidate gamma values for the radial basis kernel
#'        (default = 10^seq(-1, 1, length.out = 3)).
#' @param y.binary Logical indicating if y is binary (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the SVM model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{model}{The fitted SVM model object}
#'   \item{parameters}{List containing optimal cost and gamma values}
#' }
#'
#' @details The function iterates over a range of candidate cost and gamma values, evaluates
#' each model using cross-validation, and selects the hyperparameters that minimize the MAE.
#' Finally, the function fits an SVM model with the optimal hyperparameters to the entire
#' dataset and evaluates its performance on the test set.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 1)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- matrix(rnorm(50), ncol = 1)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Define hyperparameter grid
#' cost <- 10^seq(-1, 1, length.out = 3)
#' gamma <- 10^seq(-1, 1, length.out = 3)
#'
#' # Run the function
#' result <- get.svm.MAB(x, y, xt, yt, folds, cost, gamma)
#' print(result$MAB)
#' }
#'
#' @export
get.svm.MAB <- function(x, y, xt, yt, folds,
                        cost = 10^seq(-1, 1, length.out = 3),
                        gamma = 10^seq(-1, 1, length.out = 3),
                        y.binary = FALSE) {

    # Input validation
    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    # Handle vector inputs - convert to data frame with column name
    if (is.vector(x)) {
        x <- as.data.frame(x)
        colnames(x) <- "x"
    } else if (is.matrix(x)) {
        x <- as.data.frame(x)
    } else if (!is.data.frame(x)) {
        stop("'x' must be a vector, matrix, or data frame")
    }

    if (is.vector(xt)) {
        xt <- as.data.frame(xt)
        colnames(xt) <- "x"
    } else if (is.matrix(xt)) {
        xt <- as.data.frame(xt)
    } else if (!is.data.frame(xt)) {
        stop("'xt' must be a vector, matrix, or data frame")
    }

    if (nrow(x) != length(y)) {
        stop("Number of rows in 'x' must equal length of 'y'")
    }

    if (nrow(xt) != length(yt)) {
        stop("Number of rows in 'xt' must equal length of 'yt'")
    }

    if (ncol(x) != ncol(xt)) {
        stop("'x' and 'xt' must have the same number of columns")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.numeric(cost) || length(cost) == 0 || any(cost <= 0)) {
        stop("'cost' must be a non-empty numeric vector with positive values")
    }

    if (!is.numeric(gamma) || length(gamma) == 0 || any(gamma <= 0)) {
        stop("'gamma' must be a non-empty numeric vector with positive values")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Function to perform cross-validation
    svm.MAB.cv <- function(x, y, fold, cost, gamma) {
        train.indices <- unlist(folds[-fold])
        test.indices <- folds[[fold]]

        if (length(test.indices) == 0) {
            return(NA)
        }

        tryCatch({
            svm.model <- e1071::svm(x = x[train.indices, , drop = FALSE],
                           y = y[train.indices],
                           cost = cost,
                           gamma = gamma,
                           probability = y.binary,
                           kernel = "radial")

            predictions <- predict(svm.model, newdata = x[test.indices, , drop = FALSE])
            mae <- mean(abs(predictions - y[test.indices]), na.rm = TRUE)
            return(mae)
        }, error = function(e) {
            return(NA)
        })
    }

    # Grid search over cost and gamma
    best.MAB <- Inf
    best.params <- list(cost = NA, gamma = NA)

    for (c in cost) {
        for (g in gamma) {
            cv.MAEs <- sapply(seq_along(folds), function(fold)
                svm.MAB.cv(x, y, fold, c, g))

            # Remove NA values
            cv.MAEs <- cv.MAEs[!is.na(cv.MAEs)]

            if (length(cv.MAEs) > 0) {
                avg.MAB <- mean(cv.MAEs)

                if (avg.MAB < best.MAB) {
                    best.MAB <- avg.MAB
                    best.params <- list(cost = c, gamma = g)
                }
            }
        }
    }

    # Check if we found valid parameters
    if (is.na(best.params$cost) || is.na(best.params$gamma)) {
        stop("Failed to find optimal parameters during cross-validation")
    }

    # Fit the final model using the optimal parameters
    tryCatch({
        best.model <- e1071::svm(x = x, y = y,
                         cost = best.params$cost,
                         gamma = best.params$gamma,
                         probability = y.binary,
                         kernel = "radial")
    }, error = function(e) {
        stop(paste("Error fitting final SVM model:", e$message))
    })

    # Predict on the test set
    tryCatch({
        predictions <- predict(best.model, newdata = xt)
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Compute metrics
    residuals <- abs(predictions - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((predictions - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = predictions,
        residuals = residuals,
        model = best.model,
        parameters = list(
            cost = best.params$cost,
            gamma = best.params$gamma
        )
    )
}

#' Estimate Mean Absolute Bias (MAB) of GAM Spline Model
#'
#' Fits a Generalized Additive Model (GAM) with spline smoothing to a given dataset
#' and evaluates its performance. For binary outcomes, uses logistic GAM.
#'
#' @param x A numeric vector representing the predictor variable for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A numeric vector representing the predictor variable for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param folds A list of indices for cross-validation folds (currently not used in implementation).
#' @param y.binary Logical indicating if y is a binary variable (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the GAM model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{Empty list (for consistency with other functions)}
#' }
#'
#' @details This function fits a GAM model using spline smoothing. For binary outcomes,
#' it uses a binomial family with logit link and returns predicted probabilities.
#' For continuous outcomes, it uses a Gaussian family.
#'
#' @importFrom stats binomial
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- seq(0, 10, length.out = 50)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Run the function
#' result <- get.gam.spline.MAB(x, y, xt, yt, folds)
#' print(result$MAB)
#'
#' # Example with binary data
#' y.binary <- rbinom(100, 1, plogis(sin(x)))
#' yt.binary <- plogis(sin(xt))
#' result.binary <- get.gam.spline.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
#' print(result.binary$MAB)
#' }
#'
#' @export
get.gam.spline.MAB <- function(x, y, xt, yt, folds, y.binary = FALSE) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check for required package
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required. Please install it.")
    }

    # Create data frame for model fitting
    train.data <- data.frame(x = x, y = y)
    test.data <- data.frame(x = xt)

    # Fit GAM model
    if (y.binary) {
        # Check if y contains only 0s and 1s
        if (!all(y %in% c(0, 1))) {
            stop("For binary outcome, 'y' must contain only 0s and 1s")
        }

        # Fit GAM model with binomial family for binary outcomes
        tryCatch({
            gam.fit <- mgcv::gam(y ~ mgcv::s(x), family = binomial(), data = train.data)
        }, error = function(e) {
            stop(paste("Error fitting binary GAM model:", e$message))
        })

        # Predict probabilities
        tryCatch({
            spline.pred <- predict(gam.fit, newdata = test.data, type = "response")
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })

    } else {
        # Fit GAM model for continuous outcomes
        tryCatch({
            gam.fit <- mgcv::gam(y ~ mgcv::s(x), data = train.data)
        }, error = function(e) {
            stop(paste("Error fitting GAM model:", e$message))
        })

        # Predict values
        tryCatch({
            spline.pred <- predict(gam.fit, newdata = test.data, type = "response")
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })
    }

    # Ensure predictions are numeric vector
    if (!is.numeric(spline.pred)) {
        spline.pred <- as.numeric(spline.pred)
    }

    # Compute metrics
    residuals <- abs(spline.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((spline.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = spline.pred,
        residuals = residuals,
        parameters = list()  # Empty for consistency
    )
}

#' Estimate Mean Absolute Bias (MAB) of Smooth Spline Model
#'
#' Fits a smooth spline model to a given dataset and evaluates its performance using
#' cross-validation. The mean absolute error (MAE) is used as the performance metric.
#' The function selects the optimal degrees of freedom for the spline model based on MAE.
#'
#' @param x A numeric vector representing the predictor variable for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A numeric vector representing the predictor variable for the test set.
#' @param yt A numeric vector representing the response variable for the test set.
#' @param folds A list of indices for cross-validation folds
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the smooth spline model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the optimal degrees of freedom}
#' }
#'
#' @details The function internally performs cross-validation to select the optimal
#' degrees of freedom (between 3 and 5) that minimizes the MAE. It then fits a
#' smooth spline model with the optimal degrees of freedom to the entire dataset
#' and evaluates its performance on the test set.
#'
#' @importFrom stats smooth.spline predict
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' x <- runif(100, 0, 10)
#' y <- sin(x) + rnorm(100, sd = 0.1)
#' xt <- seq(0, 10, length.out = 50)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Run the function
#' result <- get.smooth.spline.MAB(x, y, xt, yt, folds)
#' print(result$MAB)
#' }
#'
#' @export
get.smooth.spline.MAB <- function(x, y, xt, yt, folds) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    # Check for minimum sample size
    if (length(x) < 4) {
        stop("Need at least 4 observations to fit a smooth spline")
    }

    # Range of degrees of freedom to try for model selection
    dfs <- seq(3, min(5, length(unique(x)) - 1), by = 1)

    if (length(dfs) == 0) {
        stop("Not enough unique x values to fit smooth spline with df >= 3")
    }

    # Internal function to compute MAE for given degrees of freedom
    spline.MAE.fn <- function(x, y, df, folds) {
        perform.spline.cv <- function(training.set, test.set, df) {
            tryCatch({
                # Check if we have enough data points
                if (nrow(training.set) < df + 1) {
                    return(NA)
                }

                spline.fit <- smooth.spline(x = training.set$x,
                                          y = training.set$y,
                                          df = df)
                predictions <- predict(spline.fit, x = test.set$x)$y
                mae <- mean(abs(predictions - test.set$y), na.rm = TRUE)
                return(mae)
            }, error = function(e) {
                return(NA)
            })
        }

        mae.values <- numeric(length(folds))
        for(i in seq_along(folds)) {
            train.indices <- folds[[i]]
            test.indices <- setdiff(seq_along(x), train.indices)

            if (length(test.indices) == 0) {
                mae.values[i] <- NA
                next
            }

            training.set <- data.frame(x = x[train.indices], y = y[train.indices])
            test.set <- data.frame(x = x[test.indices], y = y[test.indices])
            mae.values[i] <- perform.spline.cv(training.set, test.set, df)
        }

        # Return mean MAE, excluding NA values
        mean.MAE <- mean(mae.values, na.rm = TRUE)
        return(mean.MAE)
    }

    # Find optimal degrees of freedom
    spline.MAE <- setNames(numeric(length(dfs)), as.character(dfs))

    for (i in seq_along(dfs)) {
        df <- dfs[i]
        spline.MAE[i] <- spline.MAE.fn(x, y, df, folds)
    }

    # Remove NA values if any
    valid.indices <- !is.na(spline.MAE)
    if (sum(valid.indices) == 0) {
        stop("All degrees of freedom resulted in errors during cross-validation")
    }

    spline.MAE <- spline.MAE[valid.indices]

    # Find the degrees of freedom with the lowest MAE
    opt.df <- as.numeric(names(spline.MAE)[which.min(spline.MAE)])

    # Fit the final model using the optimal degrees of freedom
    tryCatch({
        spline.m <- smooth.spline(x = x, y = y, df = opt.df)
    }, error = function(e) {
        stop(paste("Error fitting final smooth spline model:", e$message))
    })

    # Predict using the smooth spline model
    tryCatch({
        spline.pred <- predict(spline.m, x = xt)$y
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Compute metrics
    residuals <- abs(spline.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((spline.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = spline.pred,
        residuals = residuals,
        parameters = list(opt.df = opt.df)
    )
}

#' Estimate Mean Absolute Bias (MAB) of Random Forest Model
#'
#' Fits a Random Forest model to a given dataset and evaluates its performance
#' by computing the Mean Absolute Bias on test data. Supports both regression
#' and classification tasks with optional optimization of the number of trees.
#'
#' @param x A vector, matrix, or data frame of predictor values for training.
#' @param y A vector of response values (numeric for regression, factor for classification).
#' @param xt A vector, matrix, or data frame of predictor values for testing.
#' @param yt A vector of true response values for testing. For classification,
#'   these should be probabilities of the positive class.
#' @param ntree Number of trees to grow in the random forest (default = 500).
#' @param optimize.ntree Logical. If TRUE, uses OOB error to find optimal number of trees
#'   between 100 and max.trees (default = FALSE).
#' @param max.trees Maximum number of trees to consider when optimize.ntree = TRUE
#'   (default = 1000).
#' @param plot.oob Logical. If TRUE and optimize.ntree = TRUE, plots the OOB error curve
#'   (default = FALSE).
#' @param ... Additional arguments passed to randomForest().
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the random forest model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{model}{The fitted random forest model object}
#'   \item{parameters}{List containing ntree value used and optimal.ntree if optimization was performed}
#' }
#'
#' @details For classification tasks (when y is a factor), the function returns
#' predicted probabilities for the second class level. For regression tasks,
#' it returns the predicted values directly.
#'
#' When optimize.ntree = TRUE, the function fits a random forest with max.trees
#' and examines the OOB error curve to find where the error stabilizes. It selects
#' the smallest number of trees where the OOB error is within 1% of the minimum.
#'
#' @importFrom graphics plot abline legend
#' @importFrom stats predict
#'
#' @examples
#' \dontrun{
#' # Regression example
#' set.seed(123)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- x[,1] + x[,2]^2 + rnorm(100, sd = 0.5)
#' xt <- matrix(rnorm(100), ncol = 2)
#' yt <- xt[,1] + xt[,2]^2
#'
#' # Basic usage
#' result <- get.random.forest.MAB(x, y, xt, yt, ntree = 500)
#' print(result$MAB)
#'
#' # With optimization
#' result.opt <- get.random.forest.MAB(x, y, xt, yt, optimize.ntree = TRUE)
#' print(result.opt$parameters$optimal.ntree)
#'
#' # Classification example
#' y.class <- factor(ifelse(y > median(y), "high", "low"))
#' yt.prob <- pnorm(yt, mean = mean(yt), sd = sd(yt))
#' result.class <- get.random.forest.MAB(x, y.class, xt, yt.prob)
#' print(result.class$MAB)
#' }
#'
#' @export
get.random.forest.MAB <- function(x, y, xt, yt, ntree = 500,
                                  optimize.ntree = FALSE, max.trees = 1000,
                                  plot.oob = FALSE, ...) {

    # Input validation
    if (!is.numeric(ntree) || length(ntree) != 1 || ntree < 1 || ntree != round(ntree)) {
        stop("'ntree' must be a positive integer")
    }

    if (!is.logical(optimize.ntree) || length(optimize.ntree) != 1) {
        stop("'optimize.ntree' must be TRUE or FALSE")
    }

    if (!is.numeric(max.trees) || length(max.trees) != 1 || max.trees < 1 || max.trees != round(max.trees)) {
        stop("'max.trees' must be a positive integer")
    }

    if (!is.logical(plot.oob) || length(plot.oob) != 1) {
        stop("'plot.oob' must be TRUE or FALSE")
    }

    if (!is.numeric(yt)) {
        stop("'yt' must be numeric")
    }

    # Handle different input types for x
    if (is.vector(x)) {
        x.df <- as.data.frame(x)
        colnames(x.df) <- "x"
    } else if (is.matrix(x)) {
        x.df <- as.data.frame(x)
        if (is.null(colnames(x))) {
            colnames(x.df) <- paste0("x", seq_len(ncol(x)))
        }
    } else if (is.data.frame(x)) {
        x.df <- x
    } else {
        stop("'x' must be a vector, matrix, or data frame")
    }

    # Handle different input types for xt
    if (is.vector(xt)) {
        xt.df <- as.data.frame(xt)
        colnames(xt.df) <- "x"
    } else if (is.matrix(xt)) {
        xt.df <- as.data.frame(xt)
        if (is.null(colnames(xt))) {
            colnames(xt.df) <- paste0("x", seq_len(ncol(xt)))
        }
    } else if (is.data.frame(xt)) {
        xt.df <- xt
    } else {
        stop("'xt' must be a vector, matrix, or data frame")
    }

    # Check dimensions
    if (nrow(x.df) != length(y)) {
        stop("Number of rows in 'x' must equal length of 'y'")
    }

    if (nrow(xt.df) != length(yt)) {
        stop("Number of rows in 'xt' must equal length of 'yt'")
    }

    if (ncol(x.df) != ncol(xt.df)) {
        stop("'x' and 'xt' must have the same number of columns")
    }

    # Ensure column names match
    if (!all(colnames(x.df) == colnames(xt.df))) {
        colnames(xt.df) <- colnames(x.df)
    }

    # Check y type and yt compatibility
    if (is.factor(y)) {
        if (any(yt < 0 | yt > 1)) {
            warning("'yt' contains values outside [0,1] for classification task")
        }
    } else if (!is.numeric(y)) {
        stop("'y' must be numeric for regression or factor for classification")
    }

    # Check for required package
    if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required. Please install it.")
    }

    # Initialize parameters list
    params <- list(ntree = ntree)

    # Optimize ntree if requested
    if (optimize.ntree) {
        # Fit model with max trees to get OOB error curve
        tryCatch({
            rf.full <- randomForest::randomForest(x.df, y, ntree = max.trees, keep.forest = TRUE, ...)
        }, error = function(e) {
            stop(paste("Error fitting random forest for optimization:", e$message))
        })

        # Extract OOB errors
        if (is.factor(y)) {
            # For classification, use overall OOB error rate
            oob.errors <- rf.full$err.rate[, "OOB"]
        } else {
            # For regression, use MSE
            oob.errors <- rf.full$mse
        }

        # Find where OOB error stabilizes (within 1% of minimum)
        min.error <- min(oob.errors, na.rm = TRUE)
        threshold <- min.error * 1.01

        # Find first ntree where error is below threshold
        optimal.ntree <- which(oob.errors <= threshold)[1]

        # Ensure we have at least 100 trees
        optimal.ntree <- max(100, optimal.ntree)

        if (plot.oob) {
            plot(1:max.trees, oob.errors, type = 'l',
                 xlab = "Number of Trees",
                 ylab = ifelse(is.factor(y), "OOB Error Rate", "OOB MSE"),
                 main = "OOB Error vs Number of Trees")
            abline(v = optimal.ntree, col = "red", lty = 2)
            abline(h = threshold, col = "blue", lty = 2)
            legend("topright",
                   legend = c(paste("Optimal:", optimal.ntree),
                             paste("Threshold:", round(threshold, 4))),
                   col = c("red", "blue"), lty = 2)
        }

        # Use the existing model if optimal = max, otherwise refit
        if (optimal.ntree < max.trees) {
            tryCatch({
                rf.model <- randomForest::randomForest(x.df, y, ntree = optimal.ntree, ...)
            }, error = function(e) {
                stop(paste("Error fitting final random forest model:", e$message))
            })
        } else {
            rf.model <- rf.full
        }

        # Update parameters
        params$optimal.ntree <- optimal.ntree
        params$oob.error <- oob.errors[optimal.ntree]

    } else {
        # Fit with specified ntree
        tryCatch({
            rf.model <- randomForest::randomForest(x.df, y, ntree = ntree, ...)
        }, error = function(e) {
            stop(paste("Error fitting random forest model:", e$message))
        })
    }

    # Generate predictions
    if (is.factor(y)) {
        # For classification, get probabilities
        tryCatch({
            rf.pred <- predict(rf.model, xt.df, type = "prob")[, 2]
        }, error = function(e) {
            # If binary classification with only one class in predictions
            if (ncol(predict(rf.model, xt.df, type = "prob")) == 1) {
                warning("Only one class predicted. Using class probabilities.")
                rf.pred <- predict(rf.model, xt.df, type = "prob")[, 1]
            } else {
                stop(paste("Error generating predictions:", e$message))
            }
        })
    } else {
        # For regression, get predicted values
        tryCatch({
            rf.pred <- predict(rf.model, xt.df)
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })
    }

    # Ensure predictions are numeric
    rf.pred <- as.numeric(rf.pred)

    # Compute metrics
    residuals <- abs(rf.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((rf.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = rf.pred,
        residuals = residuals,
        model = rf.model,
        parameters = params
    )
}


#' Estimate Mean Absolute Bias (MAB) of Spline Model
#'
#' Fits a spline model to a given dataset and evaluates its performance using
#' cross-validation. The function supports both continuous and binary outcomes,
#' using smooth splines for continuous data and B-splines with GLM for binary data.
#' The optimal degrees of freedom is selected based on minimizing the Mean Absolute Error.
#'
#' @param x A numeric vector representing the predictor variable for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A numeric vector representing the predictor variable for the test set.
#' @param yt A numeric vector representing the true response values for the test set.
#' @param folds A list of indices for cross-validation folds.
#' @param y.binary Logical indicating if y is a binary variable (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the spline model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing the optimal degrees of freedom}
#' }
#'
#' @details
#' For continuous outcomes (\code{y.binary = FALSE}), the function uses \code{smooth.spline}
#' to fit smoothing splines with different degrees of freedom (3, 4, 5) and selects the
#' optimal value via cross-validation.
#'
#' For binary outcomes (\code{y.binary = TRUE}), the function uses generalized linear models
#' with B-spline basis functions via \code{glm} and \code{bs}. The degrees of freedom
#' parameter controls the number of basis functions.
#'
#' @importFrom splines bs
#' @importFrom stats glm binomial predict smooth.spline
#'
#' @examples
#' \dontrun{
#' # Example with continuous outcome
#' set.seed(123)
#' x <- seq(0, 2*pi, length.out = 100)
#' y <- sin(x) + rnorm(100, sd = 0.2)
#' xt <- seq(0, 2*pi, length.out = 50)
#' yt <- sin(xt)
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 5, list = TRUE, returnTrain = TRUE)
#'
#' # Fit spline model
#' result <- get.spline.MAB(x, y, xt, yt, folds, y.binary = FALSE)
#' print(result$MAB)
#'
#' # Example with binary outcome
#' y.binary <- rbinom(100, 1, plogis(sin(x)))
#' yt.binary <- plogis(sin(xt))
#' result.binary <- get.spline.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
#' print(result.binary$MAB)
#' }
#'
#' @export
get.spline.MAB <- function(x, y, xt, yt, folds, y.binary = FALSE) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check for minimum sample size
    if (length(x) < 4) {
        stop("Need at least 4 observations to fit a spline")
    }

    # Check for required packages
    if (y.binary && !requireNamespace("splines", quietly = TRUE)) {
        stop("Package 'splines' is required for binary outcomes. Please install it.")
    }

    # Range of degrees of freedom to try for model selection
    dfs <- seq(3, min(5, length(unique(x)) - 1), by = 1)

    if (length(dfs) == 0) {
        stop("Not enough unique x values to fit spline with df >= 3")
    }

    if (y.binary) {
        # Check if y contains only 0s and 1s
        if (!all(y %in% c(0, 1))) {
            stop("For binary outcome, 'y' must contain only 0s and 1s")
        }

        # Internal function to compute MAE for binary outcomes
        spline.MAE.fn <- function(x, y, df, folds) {
            perform.spline.cv <- function(training.set, test.set, df) {
                tryCatch({
                    # Check if we have enough data points
                    if (nrow(training.set) < df + 1) {
                        return(NA)
                    }

                    spline.fit <- glm(y ~ bs(x, df = df,
                                           Boundary.knots = range(c(training.set$x, test.set$x))),
                                     family = binomial(), data = training.set)
                    spline.pred <- predict(spline.fit, newdata = test.set, type = "response")
                    mae <- mean(abs(spline.pred - test.set$y), na.rm = TRUE)
                    return(mae)
                }, error = function(e) {
                    return(NA)
                })
            }

            mae.values <- numeric(length(folds))
            for(i in seq_along(folds)) {
                train.indices <- folds[[i]]
                test.indices <- setdiff(seq_along(x), train.indices)

                if (length(test.indices) == 0) {
                    mae.values[i] <- NA
                    next
                }

                training.set <- data.frame(x = x[train.indices], y = y[train.indices])
                test.set <- data.frame(x = x[test.indices], y = y[test.indices])
                mae.values[i] <- perform.spline.cv(training.set, test.set, df)
            }

            mean.MAE <- mean(mae.values, na.rm = TRUE)
            return(mean.MAE)
        }

        # Find optimal degrees of freedom
        spline.MAE <- setNames(numeric(length(dfs)), as.character(dfs))
        for (i in seq_along(dfs)) {
            df <- dfs[i]
            spline.MAE[i] <- spline.MAE.fn(x, y, df, folds)
        }

        # Remove NA values if any
        valid.indices <- !is.na(spline.MAE)
        if (sum(valid.indices) == 0) {
            stop("All degrees of freedom resulted in errors during cross-validation")
        }

        spline.MAE <- spline.MAE[valid.indices]

        # Find the degrees of freedom with the lowest MAE
        opt.df <- as.numeric(names(spline.MAE)[which.min(spline.MAE)])

        # Fit the final model using the optimal degrees of freedom
        tryCatch({
            spline.fit <- glm(y ~ bs(x, df = opt.df, Boundary.knots = range(x)),
                             family = binomial(), data = data.frame(x = x, y = y))
        }, error = function(e) {
            stop(paste("Error fitting final spline GLM model:", e$message))
        })

        # Predict probabilities
        tryCatch({
            spline.pred <- predict(spline.fit, newdata = data.frame(x = xt), type = "response")
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })

    } else {
        # Internal function to compute MAE for continuous outcomes
        spline.MAE.fn <- function(x, y, df, folds) {
            perform.spline.cv <- function(training.set, test.set, df) {
                tryCatch({
                    # Check if we have enough data points
                    if (nrow(training.set) < df + 1) {
                        return(NA)
                    }

                    spline.fit <- smooth.spline(x = training.set$x, y = training.set$y, df = df)
                    predictions <- predict(spline.fit, x = test.set$x)$y
                    mae <- mean(abs(predictions - test.set$y), na.rm = TRUE)
                    return(mae)
                }, error = function(e) {
                    return(NA)
                })
            }

            mae.values <- numeric(length(folds))
            for(i in seq_along(folds)) {
                train.indices <- folds[[i]]
                test.indices <- setdiff(seq_along(x), train.indices)

                if (length(test.indices) == 0) {
                    mae.values[i] <- NA
                    next
                }

                training.set <- data.frame(x = x[train.indices], y = y[train.indices])
                test.set <- data.frame(x = x[test.indices], y = y[test.indices])
                mae.values[i] <- perform.spline.cv(training.set, test.set, df)
            }

            mean.MAE <- mean(mae.values, na.rm = TRUE)
            return(mean.MAE)
        }

        # Find optimal degrees of freedom
        spline.MAE <- setNames(numeric(length(dfs)), as.character(dfs))
        for (i in seq_along(dfs)) {
            df <- dfs[i]
            spline.MAE[i] <- spline.MAE.fn(x, y, df, folds)
        }

        # Remove NA values if any
        valid.indices <- !is.na(spline.MAE)
        if (sum(valid.indices) == 0) {
            stop("All degrees of freedom resulted in errors during cross-validation")
        }

        spline.MAE <- spline.MAE[valid.indices]

        # Find the degrees of freedom with the lowest MAE
        opt.df <- as.numeric(names(spline.MAE)[which.min(spline.MAE)])

        # Fit the final model using the optimal degrees of freedom
        tryCatch({
            spline.m <- smooth.spline(x = x, y = y, df = opt.df)
        }, error = function(e) {
            stop(paste("Error fitting final smooth spline model:", e$message))
        })

        # Predict using the smooth spline model
        tryCatch({
            spline.pred <- predict(spline.m, x = xt)$y
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })
    }

    # Ensure predictions are numeric
    if (!is.numeric(spline.pred)) {
        spline.pred <- as.numeric(spline.pred)
    }

    # Handle potential NA predictions
    if (any(is.na(spline.pred))) {
        warning("Some predictions are NA. These will be excluded from MAB calculation.")
    }

    # Compute metrics
    residuals <- abs(spline.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((spline.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = spline.pred,
        residuals = residuals,
        parameters = list(opt.df = opt.df)
    )
}


#' Estimate Mean Absolute Bias (MAB) of Regularized Polynomial Regression Model
#'
#' Fits a regularized polynomial regression model using elastic net regularization
#' via the glmnet package. The function creates polynomial features (x, x², x³) and
#' uses cross-validation to select the optimal regularization parameter lambda.
#'
#' @param x A numeric vector representing the predictor variable for the training set.
#' @param y A numeric vector representing the response variable for the training set.
#' @param xt A numeric vector representing the predictor variable for the test set.
#' @param yt A numeric vector representing the true response values for the test set.
#' @param folds A list of indices for cross-validation folds (currently not used).
#' @param lambda.sequence A numeric vector of lambda values to test (default = NULL, uses glmnet default).
#' @param alpha The elastic net mixing parameter (default = 1 for lasso, 0 for ridge).
#' @param y.binary Logical indicating if y is binary (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the test set}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the glmnet model on the test set}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing optimal lambda value}
#' }
#'
#' @details This function creates polynomial features (x, x², x³) to capture non-linear
#' relationships. It uses cv.glmnet for automatic cross-validation to select the optimal
#' lambda parameter. The alpha parameter controls the type of regularization:
#' \itemize{
#'   \item alpha = 1: Lasso regression (L1 penalty)
#'   \item alpha = 0: Ridge regression (L2 penalty)
#'   \item 0 < alpha < 1: Elastic net (combination of L1 and L2)
#' }
#'
#' For binary outcomes, it uses binomial family with logit link.
#' For continuous outcomes, it uses gaussian family.
#'
#' @note Consider renaming this function to get.glmnet.poly.MAB to better reflect
#' that it uses polynomial features rather than splines.
#'
#' @importFrom stats predict
#'
#' @examples
#' \dontrun{
#' # Example with continuous outcome
#' set.seed(123)
#' x <- seq(-2, 2, length.out = 100)
#' y <- x^2 + rnorm(100, sd = 0.5)
#' xt <- seq(-2, 2, length.out = 50)
#' yt <- xt^2
#'
#' # Create cross-validation folds (not used internally)
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Fit model
#' result <- get.glmnet.poly.MAB(x, y, xt, yt, folds, alpha = 1, y.binary = FALSE)
#' print(result$MAB)
#'
#' # Example with binary outcome
#' y.binary <- rbinom(100, 1, plogis(x^2))
#' yt.binary <- plogis(xt^2)
#' result.binary <- get.glmnet.poly.MAB(x, y.binary, xt, yt.binary, folds, y.binary = TRUE)
#' print(result.binary$MAB)
#' }
#'
#' @export
get.glmnet.poly.MAB <- function(x, y, xt, yt, folds,
                                lambda.sequence = NULL,
                                alpha = 1,
                                y.binary = TRUE) {

    # Input validation
    if (!is.numeric(x) || !is.numeric(y) || !is.numeric(xt) || !is.numeric(yt)) {
        stop("All inputs 'x', 'y', 'xt', and 'yt' must be numeric")
    }

    if (length(x) != length(y)) {
        stop("'x' and 'y' must have the same length")
    }

    if (length(xt) != length(yt)) {
        stop("'xt' and 'yt' must have the same length")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.null(lambda.sequence) && (!is.numeric(lambda.sequence) || any(lambda.sequence < 0))) {
        stop("'lambda.sequence' must be NULL or a numeric vector of non-negative values")
    }

    if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
        stop("'alpha' must be a single numeric value between 0 and 1")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check for required package
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required. Please install it.")
    }

    # Check sample size
    if (length(x) < 10) {
        stop("Need at least 10 observations for reliable cross-validation")
    }

    # Ensure x is a matrix
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (!is.matrix(xt)) {
        xt <- as.matrix(xt)
    }

    # If lambda.sequence is not provided, use default
    if (is.null(lambda.sequence)) {
        lambda.sequence <- 10^seq(2, -2, length = 100)
    }

    # Create polynomial features
    x.mat <- cbind(x, x^2, x^3)
    xt.mat <- cbind(xt, xt^2, xt^3)

    # Add column names to avoid warnings
    colnames(x.mat) <- c("x", "x2", "x3")
    colnames(xt.mat) <- c("x", "x2", "x3")

    # Set family based on y.binary
    family <- ifelse(y.binary, "binomial", "gaussian")

    if (y.binary) {
        # Check if y contains only 0s and 1s
        if (!all(y %in% c(0, 1))) {
            stop("For binary outcome, 'y' must contain only 0s and 1s")
        }
    }

    # Cross-validation using glmnet
    tryCatch({
        cvfit <- glmnet::cv.glmnet(x.mat, y, alpha = alpha, lambda = lambda.sequence,
                          family = family, type.measure = "mae")
    }, error = function(e) {
        stop(paste("Error in cross-validation:", e$message))
    })

    # Best lambda
    best.lambda <- cvfit$lambda.min

    # Fit the final model using the optimal lambda
    tryCatch({
        final.model <- glmnet::glmnet(x.mat, y, alpha = alpha, lambda = best.lambda,
                             family = family)
    }, error = function(e) {
        stop(paste("Error fitting final model:", e$message))
    })

    # Predict on the test set
    tryCatch({
        poly.pred <- predict(final.model, newx = xt.mat, s = best.lambda,
                            type = "response")
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Convert predictions to vector
    poly.pred <- as.vector(poly.pred)

    # Handle potential NA predictions
    if (any(is.na(poly.pred))) {
        warning("Some predictions are NA. These will be excluded from MAB calculation.")
    }

    # Compute metrics
    residuals <- abs(poly.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((poly.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = poly.pred,
        residuals = residuals,
        parameters = list(opt.lambda = best.lambda)
    )
}

#' Estimate Mean Absolute Bias for Multiple Non-linear Regression Models
#'
#' Compares the performance of different non-linear regression models by estimating
#' their Mean Absolute Bias (MAB) on synthetic data with added noise.
#'
#' @param xt A numeric vector of true predictor values. Its length must equal the number of rows of df.
#' @param df A data frame or matrix with rows corresponding to predictor values and columns
#'        corresponding to different functional relationships between predictor and response variables.
#' @param error A character string indicating the error distribution (currently only "norm" supported).
#' @param sd The standard deviation of the normal error distribution (default = 0.5).
#' @param n.subsamples The number of elements to randomly sample from the predictor and response
#'        variables. Must be less than or equal to the number of rows of df (default = 100).
#' @param x.min The minimum of the range of x values (default = 0).
#' @param x.max The maximum of the range of x values (default = 10).
#' @param n.bw The number of bandwidth values in the bandwidth grid (default = 100).
#' @param min.bw.p The minimal proportion of the x range to use as bandwidth (default = 0.01).
#' @param max.bw.p The maximal proportion of the x range to use as bandwidth (default = 0.9).
#' @param n.cores The number of cores to use for parallel computation (default = 10).
#' @param verbose Logical indicating whether to print progress messages (default = TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB.df}{A data frame of MAB estimates with rows corresponding to columns of df
#'                 and columns corresponding to different models}
#'   \item{RMSE.df}{A data frame of RMSE estimates with the same structure as MAB.df}
#'   \item{model.run.times.df}{A data frame of run times for each model (if n.cores = 1)}
#' }
#'
#' @details This function tests the following models:
#' \itemize{
#'   \item rf: Random Forest
#'   \item magelo1: Robust Local Linear Model (degree 1)
#'   \item magelo2: Robust Local Linear Model (degree 2)
#'   \item smooth.spline: Smooth spline regression
#'   \item loess: Local polynomial regression
#'   \item locpoly: Local polynomial regression (KernSmooth)
#'   \item svm: Support Vector Machine
#'   \item krr: Kernel Ridge Regression
#' }
#'
#' @importFrom stats loess smooth.spline rnorm
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 200
#' xt <- seq(0, 10, length.out = n)
#' df <- data.frame(
#'   linear = xt,
#'   quadratic = xt^2 / 10,
#'   sine = sin(xt),
#'   exponential = exp(xt/5) - 1
#' )
#'
#' # Compare models
#' results <- get.1d.models.MABs(xt, df, sd = 0.5, n.subsamples = 100, n.cores = 1)
#' print(results$MAB.df)
#' }
#'
#' @export
get.1d.models.MABs <- function(xt, df, error = "norm", sd = 0.5, n.subsamples = 100,
                               x.min = 0, x.max = 10, n.bw = 100,
                               min.bw.p = 0.01, max.bw.p = 0.9,
                               n.cores = 10, verbose = TRUE) {

    ## Input validation
    if (!requireNamespace("KernSmooth", quietly = TRUE)) {
        stop("This function requires the suggested package 'KernSmooth'. ",
             "Install it with install.packages('KernSmooth').",
             call. = FALSE)
    }

    if (!is.numeric(xt)) {
        stop("'xt' must be numeric")
    }

    if (!is.data.frame(df) && !is.matrix(df)) {
        stop("'df' must be a data frame or matrix")
    }

    if (length(xt) != nrow(df)) {
        stop("Length of 'xt' must equal the number of rows in 'df'")
    }

    if (error != "norm") {
        stop("Currently only 'norm' error distribution is supported")
    }

    if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
        stop("'sd' must be a positive numeric value")
    }

    if (!is.numeric(n.subsamples) || length(n.subsamples) != 1 ||
        n.subsamples < 1 || n.subsamples != round(n.subsamples)) {
        stop("'n.subsamples' must be a positive integer")
    }

    if (n.subsamples > nrow(df)) {
        stop("'n.subsamples' must be less than or equal to nrow(df)")
    }

    if (!is.numeric(x.min) || !is.numeric(x.max) || x.min >= x.max) {
        stop("'x.min' must be less than 'x.max'")
    }

    if (!is.numeric(n.bw) || n.bw < 1 || n.bw != round(n.bw)) {
        stop("'n.bw' must be a positive integer")
    }

    if (!is.numeric(min.bw.p) || !is.numeric(max.bw.p) ||
        min.bw.p <= 0 || max.bw.p >= 1 || min.bw.p >= max.bw.p) {
        stop("'min.bw.p' and 'max.bw.p' must be between 0 and 1, with min.bw.p < max.bw.p")
    }

    if (!is.numeric(n.cores) || n.cores < 1 || n.cores != round(n.cores)) {
        stop("'n.cores' must be a positive integer")
    }

    if (!is.logical(verbose)) {
        stop("'verbose' must be logical")
    }

    # Setup
    n.datasets <- ncol(df)
    n.samples <- nrow(df)

    if (verbose) {
        routine.ptm <- proc.time()
    }

    models <- c("rf", "magelo1", "magelo2", "smooth.spline", "loess", "locpoly", "svm", "krr")
    n.models <- length(models)

    # locpoly parameters
    min.bw <- min.bw.p * (x.max - x.min)
    max.bw <- max.bw.p * (x.max - x.min)
    bws <- seq(min.bw, max.bw, length.out = n.bw)

    # SVM parameters
    cost <- 10^seq(-1, 1, length.out = 3)
    gamma <- 10^seq(-1, 1, length.out = 3)

    # Initialize results matrices
    MAB.df <- matrix(nrow = n.datasets, ncol = n.models)
    colnames(MAB.df) <- models
    RMSE.df <- MAB.df
    model.run.times.df <- MAB.df

    if (n.cores == 1) {
        # Sequential processing
        for (dataset.i in seq_len(n.datasets)) {
            if (verbose) {
                cat("\rProcessing", dataset.i, "/", n.datasets, "dataset")
            }

            yt <- df[, dataset.i]
            ii <- sample(n.samples, size = n.subsamples)
            xp <- xt[ii]
            eps <- rnorm(n.subsamples, sd = sd)
            yp <- yt[ii] + eps

            # Creating folds
            folds <- create.folds(yp, k = 10, list = TRUE, returnTrain = TRUE)

            # Random Forest
            rf.run.time <- system.time({
                rf.res <- get.random.forest.MAB(xp, yp, xt, yt, ntree = 1000)
                MAB.df[dataset.i, "rf"] <- rf.res$MAB
                RMSE.df[dataset.i, "rf"] <- rf.res$RMSE
            })
            model.run.times.df[dataset.i, "rf"] <- rf.run.time["elapsed"]

            # MAGELO degree 1
            magelo1.run.time <- system.time({
                magelo.deg1.res <- get.magelo.MAB(xp, yp, xt, yt, deg = 1, n.cores = 1)
                MAB.df[dataset.i, "magelo1"] <- magelo.deg1.res$MAB
                RMSE.df[dataset.i, "magelo1"] <- magelo.deg1.res$RMSE
            })
            model.run.times.df[dataset.i, "magelo1"] <- magelo1.run.time["elapsed"]

            # MAGELO degree 2
            magelo2.run.time <- system.time({
                magelo.deg2.res <- get.magelo.MAB(xp, yp, xt, yt, deg = 2, n.cores = 1)
                MAB.df[dataset.i, "magelo2"] <- magelo.deg2.res$MAB
                RMSE.df[dataset.i, "magelo2"] <- magelo.deg2.res$RMSE
            })
            model.run.times.df[dataset.i, "magelo2"] <- magelo2.run.time["elapsed"]

            # Smooth spline
            spline.run.time <- system.time({
                spline.res <- get.smooth.spline.MAB(xp, yp, xt, yt, folds)
                MAB.df[dataset.i, "smooth.spline"] <- spline.res$MAB
                RMSE.df[dataset.i, "smooth.spline"] <- spline.res$RMSE
            })
            model.run.times.df[dataset.i, "smooth.spline"] <- spline.run.time["elapsed"]

            # LOESS
            loess.run.time <- system.time({
                loess.res <- get.loess.MAB(xp, yp, xt, yt, folds)
                MAB.df[dataset.i, "loess"] <- loess.res$MAB
                RMSE.df[dataset.i, "loess"] <- loess.res$RMSE
            })
            model.run.times.df[dataset.i, "loess"] <- loess.run.time["elapsed"]

            # locpoly
            locpoly.run.time <- system.time({
                locpoly.res <- get.locpoly.MAB(xp, yp, xt, yt, folds, bws)
                MAB.df[dataset.i, "locpoly"] <- locpoly.res$MAB
                RMSE.df[dataset.i, "locpoly"] <- locpoly.res$RMSE
            })
            model.run.times.df[dataset.i, "locpoly"] <- locpoly.run.time["elapsed"]

            # SVM
            svm.run.time <- system.time({
                svm.res <- get.svm.MAB(xp, yp, xt, yt, folds, cost, gamma)
                MAB.df[dataset.i, "svm"] <- svm.res$MAB
                RMSE.df[dataset.i, "svm"] <- svm.res$RMSE
            })
            model.run.times.df[dataset.i, "svm"] <- svm.run.time["elapsed"]

            # KRR
            krr.run.time <- system.time({
                krr.res <- get.krr.MAB(xp, yp, xt, yt)
                MAB.df[dataset.i, "krr"] <- krr.res$MAB
                RMSE.df[dataset.i, "krr"] <- krr.res$RMSE
            })
            model.run.times.df[dataset.i, "krr"] <- krr.run.time["elapsed"]
        }

    } else {
        # Parallel processing - model.run.times.df will be NA
        model.run.times.df <- NA

        # Process each model type in parallel across datasets
        for (model.i in seq_len(n.models)) {
            model.name <- models[model.i]

            if (model.name %in% c("magelo1", "magelo2")) {
                # MAGELO models need special handling due to internal parallelization
                deg <- ifelse(model.name == "magelo1", 1, 2)

                if (verbose) {
                    cat("\nProcessing", model.name, "models...")
                    ptm <- proc.time()
                }

                for (dataset.i in seq_len(n.datasets)) {
                    if (verbose) {
                        cat("\r  Dataset", dataset.i, "/", n.datasets)
                    }

                    yt <- df[, dataset.i]
                    ii <- sample(n.samples, size = n.subsamples)
                    xp <- xt[ii]
                    eps <- rnorm(n.subsamples, sd = sd)
                    yp <- yt[ii] + eps

                    magelo.res <- get.magelo.MAB(xp, yp, xt, yt, deg = deg, n.cores = n.cores)
                    MAB.df[dataset.i, model.name] <- magelo.res$MAB
                    RMSE.df[dataset.i, model.name] <- magelo.res$RMSE
                }

                if (verbose) {
                    cat("\r", model.name, "models completed")
                    elapsed.time(ptm, " ")
                }

            } else {
                # Other models can be parallelized across datasets
                if (verbose) {
                    cat("\nProcessing", model.name, "models...")
                    ptm <- proc.time()
                }

                registerDoParallel(n.cores)

                results <- foreach(dataset.i = seq_len(n.datasets),
                                 .combine = rbind,
                                 .packages = c("randomForest", "KernSmooth", "e1071",
                                             "CVST", "mgcv")) %dopar% {

                    yt <- df[, dataset.i]
                    ii <- sample(n.samples, size = n.subsamples)
                    xp <- xt[ii]
                    eps <- rnorm(n.subsamples, sd = sd)
                    yp <- yt[ii] + eps

                    folds <- create.folds(yp, k = 10, list = TRUE, returnTrain = TRUE)

                    if (model.name == "rf") {
                        res <- get.random.forest.MAB(xp, yp, xt, yt, ntree = 1000)
                    } else if (model.name == "smooth.spline") {
                        res <- get.smooth.spline.MAB(xp, yp, xt, yt, folds)
                    } else if (model.name == "loess") {
                        res <- get.loess.MAB(xp, yp, xt, yt, folds)
                    } else if (model.name == "locpoly") {
                        res <- get.locpoly.MAB(xp, yp, xt, yt, folds, bws)
                    } else if (model.name == "svm") {
                        res <- get.svm.MAB(xp, yp, xt, yt, folds, cost, gamma)
                    } else if (model.name == "krr") {
                        res <- get.krr.MAB(xp, yp, xt, yt)
                    }

                    c(MAB = res$MAB, RMSE = res$RMSE)
                }

                stopImplicitCluster()

                # Store results
                MAB.df[, model.name] <- results[, "MAB"]
                RMSE.df[, model.name] <- results[, "RMSE"]

                if (verbose) {
                    elapsed.time(ptm, "")
                }
            }
        }
    }

    if (verbose) {
        cat("\n")
        elapsed.time(routine.ptm, "Total elapsed time: ", with.brackets = FALSE)
    }

    # Return results
    list(
        MAB.df = as.data.frame(MAB.df),
        RMSE.df = as.data.frame(RMSE.df),
        model.run.times.df = if (n.cores == 1) as.data.frame(model.run.times.df) else NULL
    )
}

#' Estimate Mean Absolute Bias for Multiple Non-linear Regression Models with Binary Outcome
#'
#' Compares the performance of different non-linear regression models for binary outcomes
#' by estimating their Mean Absolute Bias (MAB) on synthetic data.
#'
#' @param xt A numeric vector of predictor values. Its length must equal the number of rows of df.
#' @param df A data frame or matrix of smooth function values (columns) over xt values.
#'        These will be transformed via min-max normalization and inverse logit function.
#' @param y.min The target minimum for min-max transformation (default = -3).
#' @param y.max The target maximum for min-max transformation (default = 3).
#' @param n.subsamples The number of elements to randomly sample from the predictor and response
#'        variables. Must be less than or equal to the number of rows of df (default = 100).
#' @param verbose Logical indicating whether to print progress messages (default = TRUE).
#' @param n.cores The number of cores to use for parallel computation (default = 10).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB.df}{A data frame of MAB estimates with rows corresponding to columns of df
#'                 and columns corresponding to different models}
#'   \item{RMSE.df}{A data frame of RMSE estimates with the same structure as MAB.df}
#'   \item{model.run.times.df}{A data frame of run times for each model (if n.cores = 1)}
#' }
#'
#' @details This function tests the following models for binary outcomes:
#' \itemize{
#'   \item rf: Random Forest
#'   \item magelo1: Robust Local Linear Model (degree 1)
#'   \item magelo2: Robust Local Linear Model (degree 2)
#'   \item spline: GAM spline regression
#'   \item glmnet.poly: Regularized polynomial regression
#'   \item loess: Local polynomial regression
#' }
#'
#' The function transforms the smooth functions in df to probabilities using
#' min-max normalization followed by the inverse logit transformation, then
#' generates binary outcomes by sampling from these probabilities.
#'
#' @importFrom stats loess smooth.spline rnorm rbinom plogis
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 200
#' xt <- seq(0, 10, length.out = n)
#' df <- data.frame(
#'   linear = xt,
#'   quadratic = xt^2 / 10,
#'   sine = sin(xt) * 5,
#'   exponential = exp(xt/5) - 1
#' )
#'
#' # Compare models for binary outcomes
#' results <- get.1d.binary.models.MABs(xt, df, n.subsamples = 100, n.cores = 1)
#' print(results$MAB.df)
#' }
#'
#' @export
get.1d.binary.models.MABs <- function(xt, df, y.min = -3, y.max = 3,
                                      n.subsamples = 100, verbose = TRUE,
                                      n.cores = 10) {

    # Input validation
    if (!is.numeric(xt)) {
        stop("'xt' must be numeric")
    }

    if (!is.data.frame(df) && !is.matrix(df)) {
        stop("'df' must be a data frame or matrix")
    }

    n.samples <- length(xt)
    n.datasets <- ncol(df)

    if (n.samples != nrow(df)) {
        stop("Length of 'xt' must equal the number of rows in 'df'")
    }

    if (!is.numeric(n.subsamples) || length(n.subsamples) != 1 ||
        n.subsamples < 1 || n.subsamples != round(n.subsamples)) {
        stop("'n.subsamples' must be a positive integer")
    }

    if (n.subsamples > n.samples) {
        stop("'n.subsamples' must be less than or equal to nrow(df)")
    }

    if (!is.numeric(y.min) || !is.numeric(y.max) || y.min >= y.max) {
        stop("'y.min' must be less than 'y.max'")
    }

    if (!is.logical(verbose)) {
        stop("'verbose' must be logical")
    }

    if (!is.numeric(n.cores) || n.cores < 1 || n.cores != round(n.cores)) {
        stop("'n.cores' must be a positive integer")
    }

    # Helper function for normalization and inverse logit
    normalize.and.inv.logit <- function(y, y.min, y.max) {
        # Min-max normalization to [y.min, y.max]
        y.norm <- (y - min(y)) / (max(y) - min(y)) * (y.max - y.min) + y.min
        # Inverse logit transformation
        plogis(y.norm)
    }

    if (verbose) {
        routine.ptm <- proc.time()
    }

    models <- c("rf", "magelo1", "magelo2", "spline", "glmnet.poly", "loess")
    n.models <- length(models)

    # Initialize results matrices
    MAB.df <- matrix(nrow = n.datasets, ncol = n.models)
    colnames(MAB.df) <- models
    RMSE.df <- MAB.df
    model.run.times.df <- MAB.df

    if (n.cores == 1) {
        # Sequential processing
        for (dataset.i in seq_len(n.datasets)) {
            if (verbose) {
                cat("\rProcessing", dataset.i, "/", n.datasets, "dataset")
            }

            yt <- df[, dataset.i]
            dy <- normalize.and.inv.logit(yt, y.min, y.max)  # True probability of success

            # Sample binary outcomes from probabilities
            by <- numeric(n.samples)
            for (j in seq_len(n.samples)) {
                by[j] <- rbinom(1, size = 1, prob = dy[j])
            }

            ii <- sample(n.samples, size = n.subsamples)
            xp <- xt[ii]
            yp <- by[ii]

            # Creating folds
            folds <- create.folds(yp, k = 10, list = TRUE, returnTrain = TRUE)

            # Random Forest
            rf.run.time <- system.time({
                yf <- as.factor(yp)
                rf.res <- get.random.forest.MAB(xp, yf, xt, dy, ntree = 1000)
                MAB.df[dataset.i, "rf"] <- rf.res$MAB
                RMSE.df[dataset.i, "rf"] <- rf.res$RMSE
            })
            model.run.times.df[dataset.i, "rf"] <- rf.run.time["elapsed"]

            # MAGELO degree 1
            magelo1.run.time <- system.time({
                magelo.deg1.res <- get.magelo.MAB(xp, yp, xt, dy, deg = 1, y.binary = TRUE, n.cores = 1)
                MAB.df[dataset.i, "magelo1"] <- magelo.deg1.res$MAB
                RMSE.df[dataset.i, "magelo1"] <- magelo.deg1.res$RMSE
            })
            model.run.times.df[dataset.i, "magelo1"] <- magelo1.run.time["elapsed"]

            # MAGELO degree 2
            magelo2.run.time <- system.time({
                magelo.deg2.res <- get.magelo.MAB(xp, yp, xt, dy, deg = 2, y.binary = TRUE, n.cores = 1)
                MAB.df[dataset.i, "magelo2"] <- magelo.deg2.res$MAB
                RMSE.df[dataset.i, "magelo2"] <- magelo.deg2.res$RMSE
            })
            model.run.times.df[dataset.i, "magelo2"] <- magelo2.run.time["elapsed"]

            # GAM spline
            spline.run.time <- system.time({
                spline.res <- get.gam.spline.MAB(xp, yp, xt, dy, folds, y.binary = TRUE)
                MAB.df[dataset.i, "spline"] <- spline.res$MAB
                RMSE.df[dataset.i, "spline"] <- spline.res$RMSE
            })
            model.run.times.df[dataset.i, "spline"] <- spline.run.time["elapsed"]

            # GLMnet spline
            glmnet.poly.run.time <- system.time({
                glmnet.poly.res <- get.glmnet.poly.MAB(xp, yp, xt, dy, folds, y.binary = TRUE)
                MAB.df[dataset.i, "glmnet.poly"] <- glmnet.poly.res$MAB
                RMSE.df[dataset.i, "glmnet.poly"] <- glmnet.poly.res$RMSE
            })
            model.run.times.df[dataset.i, "glmnet.poly"] <- glmnet.poly.run.time["elapsed"]

            # LOESS
            loess.run.time <- system.time({
                loess.res <- get.loess.MAB(xp, yp, xt, dy, folds)
                MAB.df[dataset.i, "loess"] <- loess.res$MAB
                RMSE.df[dataset.i, "loess"] <- loess.res$RMSE
            })
            model.run.times.df[dataset.i, "loess"] <- loess.run.time["elapsed"]
        }

    } else {
        # Parallel processing - model.run.times.df will be NA
        model.run.times.df <- NA

        # Process each model type in parallel across datasets
        for (model.i in seq_len(n.models)) {
            model.name <- models[model.i]

            if (model.name %in% c("magelo1", "magelo2")) {
                # MAGELO models need special handling due to internal parallelization
                deg <- ifelse(model.name == "magelo1", 1, 2)

                if (verbose) {
                    cat("\nProcessing", model.name, "models...")
                    ptm <- proc.time()
                }

                for (dataset.i in seq_len(n.datasets)) {
                    if (verbose) {
                        cat("\r  Dataset", dataset.i, "/", n.datasets)
                    }

                    yt <- df[, dataset.i]
                    dy <- normalize.and.inv.logit(yt, y.min, y.max)

                    by <- numeric(n.samples)
                    for (j in seq_len(n.samples)) {
                        by[j] <- rbinom(1, size = 1, prob = dy[j])
                    }

                    ii <- sample(n.samples, size = n.subsamples)
                    xp <- xt[ii]
                    yp <- by[ii]

                    magelo.res <- get.magelo.MAB(xp, yp, xt, dy, deg = deg, y.binary = TRUE, n.cores = n.cores)
                    MAB.df[dataset.i, model.name] <- magelo.res$MAB
                    RMSE.df[dataset.i, model.name] <- magelo.res$RMSE
                }

                if (verbose) {
                    cat("\r", model.name, "models completed")
                    elapsed.time(ptm, " ")
                }

            } else {
                # Other models can be parallelized across datasets
                if (verbose) {
                    cat("\nProcessing", model.name, "models...")
                    ptm <- proc.time()
                }

                registerDoParallel(n.cores)

                results <- foreach(dataset.i = seq_len(n.datasets),
                                 .combine = rbind,
                                 .packages = c("randomForest", "mgcv", "glmnet", "stats")) %dopar% {

                    yt <- df[, dataset.i]
                    dy <- normalize.and.inv.logit(yt, y.min, y.max)

                    by <- numeric(n.samples)
                    for (j in seq_len(n.samples)) {
                        by[j] <- rbinom(1, size = 1, prob = dy[j])
                    }

                    ii <- sample(n.samples, size = n.subsamples)
                    xp <- xt[ii]
                    yp <- by[ii]

                    folds <- create.folds(yp, k = 10, list = TRUE, returnTrain = TRUE)

                    if (model.name == "rf") {
                        yf <- as.factor(yp)
                        res <- get.random.forest.MAB(xp, yf, xt, dy, ntree = 1000)
                    } else if (model.name == "spline") {
                        res <- get.gam.spline.MAB(xp, yp, xt, dy, folds, y.binary = TRUE)
                    } else if (model.name == "glmnet.poly") {
                        res <- get.glmnet.poly.MAB(xp, yp, xt, dy, folds, y.binary = TRUE)
                    } else if (model.name == "loess") {
                        res <- get.loess.MAB(xp, yp, xt, dy, folds)
                    }

                    c(MAB = res$MAB, RMSE = res$RMSE)
                }

                stopImplicitCluster()

                # Store results
                MAB.df[, model.name] <- results[, "MAB"]
                RMSE.df[, model.name] <- results[, "RMSE"]

                if (verbose) {
                    elapsed.time(ptm, "")
                }
            }
        }
    }

    if (verbose) {
        cat("\n")
        elapsed.time(routine.ptm, "Total elapsed time: ", with.brackets = FALSE)
    }

    # Return results
    list(
        MAB.df = as.data.frame(MAB.df),
        RMSE.df = as.data.frame(RMSE.df),
        model.run.times.df = if (n.cores == 1) as.data.frame(model.run.times.df) else NULL
    )
}

##
## dim > 1
##

#' Estimate Mean Absolute Bias (MAB) of Random Forest Model for Multi-dimensional Data
#'
#' Fits a Random Forest model to multi-dimensional data and evaluates its performance
#' by computing the Mean Absolute Bias on the same data (in-sample evaluation).
#'
#' @param X A matrix or data frame of predictor variables (dimension > 1).
#' @param y A vector of response values (numeric for regression, factor for classification).
#' @param yt A vector of true response values for evaluation.
#' @param ntree Number of trees to grow in the random forest (default = 500).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the dataset}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the random forest model}
#'   \item{residuals}{Absolute residuals}
#'   \item{model}{The fitted random forest model object}
#'   \item{parameters}{List containing ntree value used}
#' }
#'
#' @details For classification tasks (when y is a factor), the function returns
#' predicted probabilities for the second class level. For regression tasks,
#' it returns the predicted values directly. Note that this function evaluates
#' the model on the training data (in-sample).
#'
#' @examples
#' \dontrun{
#' # Regression example
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
#' yt <- X[,1] + X[,2]^2  # True values
#'
#' result <- get.random.forest.MAB.xD(X, y, yt, ntree = 500)
#' print(result$MAB)
#'
#' # Classification example
#' y.class <- factor(ifelse(y > median(y), "high", "low"))
#' yt.prob <- plogis(yt - median(yt))  # True probabilities
#' result.class <- get.random.forest.MAB.xD(X, y.class, yt.prob, ntree = 500)
#' print(result.class$MAB)
#' }
#'
#' @export
get.random.forest.MAB.xD <- function(X, y, yt, ntree = 500) {

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("'X' must be a matrix or data frame")
    }

    if (is.data.frame(X)) {
        # Check if all columns are numeric
        if (!all(sapply(X, is.numeric))) {
            stop("All columns of 'X' must be numeric")
        }
    }

    if (!is.numeric(ntree) || length(ntree) != 1 || ntree < 1 || ntree != round(ntree)) {
        stop("'ntree' must be a positive integer")
    }

    if (nrow(X) != length(y)) {
        stop("Number of rows in 'X' must equal length of 'y'")
    }

    if (nrow(X) != length(yt)) {
        stop("Number of rows in 'X' must equal length of 'yt'")
    }

    # Check y type and yt compatibility
    if (is.factor(y)) {
        if (!is.numeric(yt)) {
            stop("When 'y' is a factor, 'yt' must be numeric (probabilities)")
        }
        if (any(yt < 0 | yt > 1)) {
            warning("'yt' contains values outside [0,1] for classification task")
        }
    } else if (!is.numeric(y)) {
        stop("'y' must be numeric for regression or factor for classification")
    }

    if (!is.numeric(yt)) {
        stop("'yt' must be numeric")
    }

    # Check for required package
    if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required. Please install it.")
    }

    # Fit random forest model
    tryCatch({
        rf.model <- randomForest::randomForest(X, y, ntree = ntree)
    }, error = function(e) {
        stop(paste("Error fitting random forest model:", e$message))
    })

    # Generate predictions (in-sample)
    if (is.factor(y)) {
        # For classification, get probabilities
        tryCatch({
            rf.pred <- predict(rf.model, type = "prob")[, 2]
        }, error = function(e) {
            # If binary classification with only one class in predictions
            if (ncol(predict(rf.model, type = "prob")) == 1) {
                warning("Only one class predicted. Using class probabilities.")
                rf.pred <- predict(rf.model, type = "prob")[, 1]
            } else {
                stop(paste("Error generating predictions:", e$message))
            }
        })
    } else {
        # For regression, get predicted values
        tryCatch({
            rf.pred <- predict(rf.model)
        }, error = function(e) {
            stop(paste("Error generating predictions:", e$message))
        })
    }

    # Ensure predictions are numeric
    rf.pred <- as.numeric(rf.pred)

    # Compute metrics
    residuals <- abs(rf.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((rf.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = rf.pred,
        residuals = residuals,
        model = rf.model,
        parameters = list(ntree = ntree)
    )
}

#' Estimate Mean Absolute Bias (MAB) of Nonparametric Kernel Regression Model
#'
#' Fits a nonparametric kernel regression model to multi-dimensional data and evaluates
#' its performance using in-sample predictions.
#'
#' @param X A matrix or data frame of predictor variables (dimension > 1).
#' @param y A numeric vector of response values.
#' @param yt A numeric vector of true response values for evaluation.
#' @param y.binary Logical indicating if y is a binary variable (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the dataset}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the nonparametric model}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{Empty list (for consistency with other functions)}
#' }
#'
#' @details This function fits a nonparametric kernel regression model using the np package.
#' For continuous outcomes, it uses local linear regression with bandwidth selection via
#' cross-validated AIC. For binary outcomes, it uses conditional mode estimation.
#' Note that bandwidth selection can be computationally intensive for large datasets.
#'
#' @importFrom stats as.formula
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' library(np)
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2)
#' colnames(X) <- c("x1", "x2")
#' y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
#' yt <- X[,1] + X[,2]^2  # True values
#'
#' # Fit nonparametric model and compute MAB
#' result <- get.np.MAB.xD(X, y, yt)
#' print(result$MAB)
#'
#' # Example with binary data
#' y.binary <- rbinom(n, 1, plogis(y))
#' yt.binary <- plogis(yt)
#' result.binary <- get.np.MAB.xD(X, y.binary, yt.binary, y.binary = TRUE)
#' print(result.binary$MAB)
#' }
#'
#' @export
get.np.MAB.xD <- function(X, y, yt, y.binary = FALSE) {

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("'X' must be a matrix or data frame")
    }

    # Convert matrix to data frame with column names if needed
    if (is.matrix(X)) {
        X <- as.data.frame(X)
        if (is.null(colnames(X))) {
            colnames(X) <- paste0("x", seq_len(ncol(X)))
        }
    }

    # Check if all columns are numeric
    if (!all(sapply(X, is.numeric))) {
        stop("All columns of 'X' must be numeric")
    }

    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    if (nrow(X) != length(y)) {
        stop("Number of rows in 'X' must equal length of 'y'")
    }

    if (nrow(X) != length(yt)) {
        stop("Number of rows in 'X' must equal length of 'yt'")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check dimensions
    if (ncol(X) < 2) {
        stop("'X' must have at least 2 columns for multi-dimensional analysis")
    }

    # Check for required package
    if (!requireNamespace("np", quietly = TRUE)) {
        stop("Package 'np' is required. Please install it.")
    }

    # Create formula
    predictors <- colnames(X)
    formula.str <- paste0("y ~ ", paste(predictors, collapse = " + "))
    np.formula <- as.formula(formula.str)

    if (y.binary) {
        # Check if y contains only 0s and 1s
        if (!all(y %in% c(0, 1))) {
            stop("For binary outcome, 'y' must contain only 0s and 1s")
        }

        # Create data frame with factor outcome for binary case
        data.df <- data.frame(X, y = factor(y))

        # Fit nonparametric conditional mode model for binary outcomes
        tryCatch({
            np.bin.m <- np::npconmode(np.formula, data = data.df)
            np.pred <- np.bin.m$condens
        }, error = function(e) {
            stop(paste("Error fitting binary nonparametric model:", e$message))
        })

    } else {
        # Create data frame for continuous case
        data.df <- data.frame(X, y = y)

        # Fit nonparametric regression model
        tryCatch({
            # Bandwidth selection
            bw.all <- np::npregbw(formula = np.formula,
                            regtype = "ll",
                            bwmethod = "cv.aic",
                            data = data.df)
        }, error = function(e) {
            stop(paste("Error in bandwidth selection:", e$message))
        })

        # Fit the model with selected bandwidths
        tryCatch({
            np.model <- np::npreg(bws = bw.all)
            np.pred <- predict(np.model)
        }, error = function(e) {
            stop(paste("Error fitting nonparametric model:", e$message))
        })
    }

    # Ensure predictions are numeric
    if (!is.numeric(np.pred)) {
        np.pred <- as.numeric(np.pred)
    }

    # Compute metrics
    residuals <- abs(np.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((np.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = np.pred,
        residuals = residuals,
        parameters = list()  # Empty for consistency
    )
}

#' Estimate Mean Absolute Bias (MAB) of GAM Spline Model for Multi-dimensional Data
#'
#' Fits a Generalized Additive Model (GAM) with spline smoothing to multi-dimensional
#' data and evaluates its performance using in-sample predictions.
#'
#' @param X A matrix or data frame of predictor variables (dimension > 1).
#' @param y A numeric vector of response values.
#' @param yt A numeric vector of true response values for evaluation.
#' @param folds A list of indices for cross-validation folds (currently not used in implementation).
#' @param y.binary Logical indicating if y is a binary variable (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the dataset}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the GAM model}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{Empty list (for consistency with other functions)}
#' }
#'
#' @details This function fits a GAM model with smooth terms for each predictor.
#' For binary outcomes, it uses a binomial family with logit link and returns
#' predicted probabilities. For continuous outcomes, it uses a Gaussian family.
#' Note that this creates a tensor product smooth of all predictors.
#'
#' @importFrom stats binomial as.formula
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2)
#' colnames(X) <- c("x1", "x2")
#' y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
#' yt <- X[,1] + X[,2]^2  # True values
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Fit GAM model and compute MAB
#' result <- get.spline.MAB.xD(X, y, yt, folds)
#' print(result$MAB)
#'
#' # Example with binary data
#' y.binary <- rbinom(n, 1, plogis(y))
#' yt.binary <- plogis(yt)
#' result.binary <- get.spline.MAB.xD(X, y.binary, yt.binary, folds, y.binary = TRUE)
#' print(result.binary$MAB)
#' }
#'
#' @export
get.spline.MAB.xD <- function(X, y, yt, folds, y.binary = FALSE) {

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("'X' must be a matrix or data frame")
    }

    # Convert matrix to data frame with column names if needed
    if (is.matrix(X)) {
        X <- as.data.frame(X)
        if (is.null(colnames(X))) {
            colnames(X) <- paste0("x", seq_len(ncol(X)))
        }
    }

    # Check if all columns are numeric
    if (!all(sapply(X, is.numeric))) {
        stop("All columns of 'X' must be numeric")
    }

    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    if (nrow(X) != length(y)) {
        stop("Number of rows in 'X' must equal length of 'y'")
    }

    if (nrow(X) != length(yt)) {
        stop("Number of rows in 'X' must equal length of 'yt'")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check dimensions
    if (ncol(X) < 2) {
        stop("'X' must have at least 2 columns for multi-dimensional analysis")
    }

    # Check for required package
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required. Please install it.")
    }

    # Create data frame for model fitting
    data.df <- cbind(y = y, X)

    # Create formula for GAM
    predictors <- colnames(X)
    # Create tensor product smooth of all predictors
    formula.str <- paste0("y ~ mgcv::s(", paste(predictors, collapse = ", "), ")")
    gam.formula <- as.formula(formula.str)

    # Fit GAM model
    if (y.binary) {
        # Check if y contains only 0s and 1s
        if (!all(y %in% c(0, 1))) {
            stop("For binary outcome, 'y' must contain only 0s and 1s")
        }

        # Fit GAM model with binomial family for binary outcomes
        tryCatch({
            gam.fit <- mgcv::gam(gam.formula, data = data.df, family = binomial())
        }, error = function(e) {
            stop(paste("Error fitting binary GAM model:", e$message))
        })

    } else {
        # Fit GAM model for continuous outcomes
        tryCatch({
            gam.fit <- mgcv::gam(gam.formula, data = data.df)
        }, error = function(e) {
            stop(paste("Error fitting GAM model:", e$message))
        })
    }

    # Generate predictions (in-sample)
    tryCatch({
        if (y.binary) {
            spline.pred <- predict(gam.fit, type = "response")
        } else {
            spline.pred <- predict(gam.fit, type = "response")
        }
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Ensure predictions are numeric vector
    if (!is.numeric(spline.pred)) {
        spline.pred <- as.numeric(spline.pred)
    }

    # Compute metrics
    residuals <- abs(spline.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((spline.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = spline.pred,
        residuals = residuals,
        parameters = list()  # Empty for consistency
    )
}

#' Estimate Mean Absolute Bias (MAB) of SVM Model for Multi-dimensional Data
#'
#' Fits a Support Vector Machine (SVM) model to multi-dimensional data and evaluates
#' its performance using cross-validation. The optimal hyperparameters are selected
#' based on cross-validated MAE.
#'
#' @param X A matrix or data frame of predictor variables (dimension > 1).
#' @param y A numeric vector of response values.
#' @param yt A numeric vector of true response values for evaluation.
#' @param folds A list of indices for cross-validation folds.
#' @param cost A numeric vector of candidate cost values to evaluate.
#' @param gamma A numeric vector of candidate gamma values for the radial basis kernel.
#' @param y.binary Logical indicating if y is binary (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the dataset}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the SVM model}
#'   \item{residuals}{Absolute residuals}
#'   \item{model}{The fitted SVM model object}
#'   \item{parameters}{List containing optimal cost and gamma values}
#' }
#'
#' @details The function performs grid search over cost and gamma parameters using
#' cross-validation to find the optimal values. It then fits a final SVM model
#' with these parameters and evaluates on the training data (in-sample).
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
#' yt <- X[,1] + X[,2]^2  # True values
#'
#' # Create cross-validation folds
#' folds <- create.folds(y, k = 10, list = TRUE, returnTrain = TRUE)
#'
#' # Define hyperparameter grid
#' cost <- 10^seq(-1, 1, length.out = 3)
#' gamma <- 10^seq(-1, 1, length.out = 3)
#'
#' # Fit SVM model and compute MAB
#' result <- get.svm.MAB.xD(X, y, yt, folds, cost, gamma)
#' print(result$MAB)
#' }
#'
#' @export
get.svm.MAB.xD <- function(X, y, yt, folds, cost, gamma, y.binary = FALSE) {

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("'X' must be a matrix or data frame")
    }

    if (is.data.frame(X)) {
        # Check if all columns are numeric
        if (!all(sapply(X, is.numeric))) {
            stop("All columns of 'X' must be numeric")
        }
    }

    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    if (nrow(X) != length(y)) {
        stop("Number of rows in 'X' must equal length of 'y'")
    }

    if (nrow(X) != length(yt)) {
        stop("Number of rows in 'X' must equal length of 'yt'")
    }

    if (!is.list(folds)) {
        stop("'folds' must be a list of fold indices")
    }

    if (!is.numeric(cost) || length(cost) == 0 || any(cost <= 0)) {
        stop("'cost' must be a non-empty numeric vector with positive values")
    }

    if (!is.numeric(gamma) || length(gamma) == 0 || any(gamma <= 0)) {
        stop("'gamma' must be a non-empty numeric vector with positive values")
    }

    if (!is.logical(y.binary)) {
        stop("'y.binary' must be logical (TRUE/FALSE)")
    }

    # Check dimensions
    if (ncol(X) < 2) {
        stop("'X' must have at least 2 columns for multi-dimensional analysis")
    }

    # Check for required package
    if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' is required. Please install it.")
    }

    # Function to perform cross-validation
    svm.MAB.cv <- function(X, y, fold, cost, gamma) {
        train.indices <- unlist(folds[-fold])
        test.indices <- folds[[fold]]

        if (length(test.indices) == 0) {
            return(NA)
        }

        tryCatch({
            svm.model <- e1071::svm(x = X[train.indices, , drop = FALSE],
                           y = y[train.indices],
                           cost = cost,
                           gamma = gamma,
                           probability = y.binary,
                           kernel = "radial")

            svm.pred <- predict(svm.model, newdata = X[test.indices, , drop = FALSE])
            mae <- mean(abs(svm.pred - y[test.indices]), na.rm = TRUE)
            return(mae)
        }, error = function(e) {
            return(NA)
        })
    }

    # Grid search over cost and gamma
    best.MAB <- Inf
    best.params <- list(cost = NA, gamma = NA)

    for (c in cost) {
        for (g in gamma) {
            cv.MAEs <- sapply(seq_along(folds), function(fold)
                svm.MAB.cv(X, y, fold, c, g))

            # Remove NA values
            cv.MAEs <- cv.MAEs[!is.na(cv.MAEs)]

            if (length(cv.MAEs) > 0) {
                avg.MAB <- mean(cv.MAEs)

                if (avg.MAB < best.MAB) {
                    best.MAB <- avg.MAB
                    best.params <- list(cost = c, gamma = g)
                }
            }
        }
    }

    # Check if we found valid parameters
    if (is.na(best.params$cost) || is.na(best.params$gamma)) {
        stop("Failed to find optimal parameters during cross-validation")
    }

    # Fit the final model using the optimal parameters
    tryCatch({
        best.model <- e1071::svm(x = X, y = y,
                         cost = best.params$cost,
                         gamma = best.params$gamma,
                         probability = y.binary,
                         kernel = "radial")
    }, error = function(e) {
        stop(paste("Error fitting final SVM model:", e$message))
    })

    # Predict on the data (in-sample)
    tryCatch({
        svm.pred <- predict(best.model)
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Ensure predictions are numeric
    if (!is.numeric(svm.pred)) {
        svm.pred <- as.numeric(svm.pred)
    }

    # Compute metrics
    residuals <- abs(svm.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((svm.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = svm.pred,
        residuals = residuals,
        model = best.model,
        parameters = list(
            cost = best.params$cost,
            gamma = best.params$gamma
        )
    )
}


#' Estimate Mean Absolute Bias (MAB) of Kernel Ridge Regression for Multi-dimensional Data
#'
#' Fits a Kernel Ridge Regression (KRR) model to multi-dimensional data and evaluates
#' its performance using cross-validation to select optimal hyperparameters.
#'
#' @param X A matrix or data frame of predictor variables (dimension > 1).
#' @param y A numeric vector of response values.
#' @param yt A numeric vector of true response values for evaluation.
#' @param kernel A kernel function to be used in KRR (default = 'rbfdot').
#' @param lambdas A numeric vector of regularization parameters to test
#'        (default = 10^(-8:0)).
#' @param sigmas A numeric vector of kernel width parameters to test
#'        (default = 10^((1:9)/3)).
#'
#' @return A list containing:
#' \describe{
#'   \item{MAB}{Mean Absolute Bias calculated on the dataset}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{predictions}{Predictions made by the KRR model}
#'   \item{residuals}{Absolute residuals}
#'   \item{parameters}{List containing optimal sigma and lambda values}
#' }
#'
#' @details The function uses Kernel Ridge Regression to model the relationship
#' between the predictors and the response. It performs 10-fold cross-validation
#' to find optimal hyperparameters (sigma and lambda) and computes the MAB using
#' in-sample predictions.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' y <- X[,1] + X[,2]^2 + rnorm(n, sd = 0.1)
#' yt <- X[,1] + X[,2]^2  # True values
#'
#' # Fit KRR model and compute MAB
#' result <- get.krr.MAB.xD(X, y, yt)
#' print(result$MAB)
#'
#' # With custom parameters
#' result2 <- get.krr.MAB.xD(X, y, yt, kernel = 'polydot',
#'                          lambdas = 10^seq(-6, -2, length = 5))
#' print(result2$MAB)
#' }
#'
#' @export
get.krr.MAB.xD <- function(X, y, yt, kernel = 'rbfdot',
                           lambdas = 10^(-8:0),
                           sigmas = 10^((1:9)/3)) {

    if (!requireNamespace("CVST", quietly = TRUE)) {
        stop("This function requires the suggested package 'CVST'. ",
             "Install it with install.packages('CVST').",
             call. = FALSE)
    }

    # Input validation
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("'X' must be a matrix or data frame")
    }

    # Convert data frame to matrix if needed
    if (is.data.frame(X)) {
        # Check if all columns are numeric
        if (!all(sapply(X, is.numeric))) {
            stop("All columns of 'X' must be numeric")
        }
        X <- as.matrix(X)
    }

    if (!is.numeric(y) || !is.numeric(yt)) {
        stop("'y' and 'yt' must be numeric")
    }

    if (nrow(X) != length(y)) {
        stop("Number of rows in 'X' must equal length of 'y'")
    }

    if (nrow(X) != length(yt)) {
        stop("Number of rows in 'X' must equal length of 'yt'")
    }

    if (!is.character(kernel) || length(kernel) != 1) {
        stop("'kernel' must be a single character string")
    }

    # Valid kernels for KRR
    valid.kernels <- c('rbfdot', 'polydot', 'laplacedot', 'tanhdot', 'vanilladot')

    if (!(kernel %in% valid.kernels)) {
        warning(paste("'kernel' may not be valid. Common options are:",
                     paste(valid.kernels, collapse = ", ")))
    }

    if (!is.numeric(lambdas) || length(lambdas) == 0 || any(lambdas < 0)) {
        stop("'lambdas' must be a non-empty numeric vector of non-negative values")
    }

    if (!is.numeric(sigmas) || length(sigmas) == 0 || any(sigmas <= 0)) {
        stop("'sigmas' must be a non-empty numeric vector of positive values")
    }

    # Check dimensions
    if (ncol(X) < 2) {
        stop("'X' must have at least 2 columns for multi-dimensional analysis")
    }

    # Check for required package
    if (!requireNamespace("CVST", quietly = TRUE)) {
        stop("Package 'CVST' is required. Please install it.")
    }

    # Construct KRR learner
    krr <- CVST::constructKRRLearner()

    # Construct parameters
    tryCatch({
        params <- CVST::constructParams(kernel = kernel, sigma = sigmas, lambda = lambdas)
    }, error = function(e) {
        stop(paste("Error constructing parameters:", e$message))
    })

    # Construct data
    dat <- CVST::constructData(X, y)

    # Perform cross-validation to find optimal parameters
    tryCatch({
        opt <- CVST::CV(dat, krr, params, fold = 10, verbose = FALSE)
    }, error = function(e) {
        stop(paste("Error during cross-validation:", e$message))
    })

    # Extract optimal parameters
    opt.sigma <- opt[[1]]$sigma
    opt.lambda <- opt[[1]]$lambda

    # Set final parameters
    param <- list(kernel = kernel, sigma = opt.sigma, lambda = opt.lambda)

    # Train final model with optimal parameters
    tryCatch({
        krr.model <- krr$learn(dat, param)
    }, error = function(e) {
        stop(paste("Error training KRR model:", e$message))
    })

    # Predict on training data (in-sample)
    dat.train <- CVST::constructData(X, rep(0, nrow(X)))  # Dummy y values for prediction

    tryCatch({
        krr.pred <- krr$predict(krr.model, dat.train)[, 1]
    }, error = function(e) {
        stop(paste("Error generating predictions:", e$message))
    })

    # Ensure predictions are numeric
    if (!is.numeric(krr.pred)) {
        krr.pred <- as.numeric(krr.pred)
    }

    # Compute metrics
    residuals <- abs(krr.pred - yt)
    MAB <- mean(residuals, na.rm = TRUE)
    RMSE <- sqrt(mean((krr.pred - yt)^2, na.rm = TRUE))

    # Return standardized output
    list(
        MAB = MAB,
        RMSE = RMSE,
        predictions = krr.pred,
        residuals = residuals,
        parameters = list(
            kernel = kernel,
            sigma = opt.sigma,
            lambda = opt.lambda
        )
    )
}
