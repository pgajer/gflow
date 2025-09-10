# Global variable declarations for R CMD check
if(getRversion() >= "2.15.1") utils::globalVariables(c("idx"))

#' Model Averaged Grid-based Epsilon LOwess (MAGELO)
#'
#' @description
#' A non-parametric regression method that combines local polynomial regression with model
#' averaging using disk-shaped neighborhoods. Unlike traditional LOWESS which uses k-nearest
#' neighbors, MAGELO employs fixed-radius neighborhoods (disks) and averages local polynomial
#' models fitted at uniformly spaced grid points. This approach provides smoother transitions
#' between local models and more stable estimates in regions with varying data density.
#'
#' @param x Numeric vector of predictor values
#' @param y Numeric vector of response values, must be same length as x
#' @param grid.size Integer specifying the number of grid points for model fitting.
#'        Default is 400, which provides good balance between accuracy and computation speed
#' @param degree Integer specifying polynomial degree (0, 1, or 2):
#'        - 0: weighted mean regression
#'        - 1: local linear regression (y ~ x)
#'        - 2: local quadratic regression (y ~ poly(x, 2))
#' @param min.K Integer specifying minimum number of points required in each local neighborhood
#' @param f Numeric between 0 and 1, specifying the proportion of x range to use for each local window.
#'        If NULL, optimal value is determined by cross-validation
#' @param bw Numeric specifying bandwidth (radius of disk neighborhood).
#'        If NULL, optimal value is determined by cross-validation
#' @param min.bw.f Numeric specifying minimum bandwidth as fraction of x range.
#'        min_bandwidth = min.bw.f * range(x). Default: 0.025
#' @param method Character string specifying bandwidth optimization method:
#'        - "LOOCV": Leave-one-out cross-validation (for n < 1000)
#'        - "CV": K-fold cross-validation (for n >= 1000)
#' @param n.bws Integer specifying number of bandwidths to try during optimization
#' @param n.cv.folds Integer specifying number of folds for cross-validation
#' @param n.cv.reps Integer specifying number of cross-validation repetitions
#' @param nn.kernel Character string specifying kernel function:
#'        "epanechnikov", "triangular", "tr.exponential", or "normal"
#' @param n.BB Integer specifying number of Bayesian bootstrap iterations
#' @param get.predictions.CrI Logical; if TRUE, compute credible intervals for fitted values at x
#' @param get.gpredictions.CrI Logical; if TRUE, compute credible intervals for fitted values at grid points
#' @param get.BB.predictions Logical; if TRUE, return matrix of bootstrap estimates at x
#' @param get.BB.gpredictions Logical; if TRUE, return matrix of bootstrap estimates at grid points
#' @param level Numeric between 0 and 1 specifying credible interval level
#' @param n.C.itr Integer specifying number of Cleveland's iterative re-weighting steps
#' @param C Numeric specifying scaling factor for residuals in robust fitting
#' @param stop.C.itr.on.min Logical; if TRUE, stop iterations when improvement plateaus
#' @param y.binary Logical; if TRUE, use binary loss function for optimization
#' @param cv.nNN Integer specifying number of nearest neighbors for grid interpolation
#' @param n.perms Integer specifying number of y permutations for p-value calculation
#' @param y.true Optional numeric vector of true response values (for validation)
#' @param use.binloss Logical; if TRUE, use binary loss function for optimization
#' @param n.cores Integer specifying number of CPU cores for parallel processing
#' @param verbose Logical; if TRUE, print progress information
#'
#' @section Uncertainty Estimation:
#' The function provides several options for uncertainty quantification through Bayesian bootstrap.
#' Use \code{n.BB} to specify the number of bootstrap iterations, and control which intervals
#' to compute with \code{get.predictions.CrI} and \code{get.gpredictions.CrI}. The credible
#' interval level can be adjusted with the \code{level} parameter.
#'
#' @section Robustness Parameters:
#' MAGELO includes options for robust fitting using Cleveland's iterative re-weighting procedure.
#' The \code{n.C.itr} parameter controls the number of iterations, while \code{C} sets the
#' scaling factor for residuals. Use \code{stop.C.itr.on.min} to enable early stopping when
#' improvements plateau.
#'
#' @section Permutation Testing: Statistical significance can be assessed
#'     through permutation testing. Set \code{n.perms} for standard permutation
#'     tests. If true response values are available, provide them via
#'     \code{y.true} for validation purposes.
#'
#' @return A list with class "magelo" containing:
#' \itemize{
#'   \item Fitted values at grid points (gpredictions) and input locations (predictions)
#'   \item Credible intervals if requested
#'   \item Bootstrap estimates if requested
#'   \item Optimal bandwidth and optimization results
#'   \item Model coefficients
#'   \item Input parameters
#' }
#'
#' @examples
#' x <- seq(0, 10, length.out = 100)
#' y <- sin(x) + rnorm(100, 0, 0.1)
#' fit <- magelo(x, y, degree = 1)
#' plot(x, y)
#' lines(fit$xgrid, fit$gpredictions, col = "red")
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel stopImplicitCluster registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats approx weighted.mean lm predict quantile sd var
#' @importFrom graphics lines par plot points polygon
#' @export
magelo <- function(x,
                   y,
                   y.true = NULL,
                   grid.size = 400,
                   degree = 1,
                   min.K = 5,
                   f = NULL,
                   bw = NULL,
                   min.bw.f = 0.025,
                   method = ifelse(length(x) < 1000, "LOOCV", "CV"),
                   n.bws = 100,
                   n.BB = 1000,
                   get.predictions.CrI = TRUE,
                   get.gpredictions.CrI = TRUE,
                   get.BB.predictions = FALSE,
                   get.BB.gpredictions = FALSE,
                   level = 0.95,
                   n.C.itr = 0,
                   C = 6,
                   stop.C.itr.on.min = TRUE,
                   n.cv.folds = 10,
                   n.cv.reps = 20,
                   nn.kernel = "epanechnikov",
                   y.binary = FALSE,
                   cv.nNN = 3,
                   n.perms  =  0,
                   n.cores = 1,
                   use.binloss = FALSE,
                   verbose = FALSE) {
    stopifnot(is.numeric(x))
    nx <- length(x)

    stopifnot(is.numeric(y))
    stopifnot(length(y) == nx)

    stopifnot(min.K < nx)

    if ( !is.null(f) )
    {
        stopifnot(is.numeric(f))
        stopifnot(f>0 && f <=1)
    }

    stopifnot(as.integer(degree)==degree)
    stopifnot(degree==0 || degree==1 || degree==2)

    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    nn.kernel <- match.arg(nn.kernel, kernels)
    ikernel <- pmatch(nn.kernel, kernels)

    if ( y.binary && n.C.itr != 0 ) {
        warning("In y binary mode, n.C.itr needs to be 0. Setting it to 0.")
        n.C.itr <- 0
    }

    ## Sorting x and y
    o <- order(x)
    x <- x[o]
    y <- y[o]

    ng <- grid.size

    ## Defining uniform grid over the range of x values
    x.range <- range(x)
    x.width <- diff(x.range)
    xgrid <- seq(x.range[1], x.range[2], length = ng)
    x.widthg <- xgrid[2] - xgrid[1]

    ## Defining bandwidth, bw, as the half of the window size defined as f*x.width
    if ( is.null(bw) && !is.null(f) ) {
        window.size <- f * x.width
        bw <- window.size / 2
    }

    ## Identifying x NN's from xgrid
    ## They are used for training of linear models and prediction of predictions over x
    nn <- get.knnx(x, xgrid, k = nx)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    errors <- NULL
    min.error <- NULL
    min.bw <- min.bw.f * x.width
    max.bw <- 0.9 * x.width
    opt.bw.i <- NULL
    opt.bw <- NULL
    log.bws <- seq(log(min.bw), log(max.bw), length=n.bws)

    # Define opt.bw.fn functions first, before using them
    if ( is.null(bw) && degree == 0 && !use.binloss) {

        K <- ncol(nn.i)

        if ( verbose ) {
            cat("Optimizing bw using cv_deg0_mae() ")
            ptm <- proc.time()
        }

        # Define the optimization function for degree 0, MAE
        opt.bw.fn.deg0.mae <- function(bw) {
            mae <- 0
            out <- .C("C_cv_deg0_mae",
                     as.integer(n.cv.folds),
                     as.integer(n.cv.reps),
                     as.integer(cv.nNN),
                     as.integer(y.binary),
                     as.integer(t(nn.i-1)),
                     as.double(t(nn.d)),
                     as.double(y),
                     as.integer(K),
                     as.integer(ng),
                     as.integer(nx),
                     as.double(bw),
                     as.integer(min.K),
                     as.integer(ikernel),
                     mae = as.double(mae))

            if ( verbose ) cat(".")

            out$mae
        }

        loop.ptm <- proc.time()
        if (n.cores > 1) {
            doParallel::registerDoParallel(n.cores)
        } else {
            foreach::registerDoSEQ()
        }
        errors <- foreach::foreach (idx = seq(n.bws),
                           .combine = c,
                           .export = c(
                               "nn.i","nn.d","y","n.cv.folds","n.cv.reps","cv.nNN","y.binary",
                               "ng","nx","min.K","ikernel","log.bws","opt.bw.fn.deg0.mae"
                           )
                           ) %dopar% {
            if ( verbose ) { cat("\ri=", idx) }
            opt.bw.fn.deg0.mae(log.bws[idx])
        }
        doParallel::stopImplicitCluster()
        if ( verbose ) {
            cat("\r\n")
            elapsed.time(loop.ptm)
        }

        opt.bw.i <- which.min(errors)
        opt.bw <- exp(log.bws[opt.bw.i])
        min.error <- errors[opt.bw.i]

        if ( verbose ) {
            cat("opt.bw.i: ", opt.bw.i, "\n")
            cat("opt.bw: ", opt.bw, "\n")
            elapsed.time(ptm)
        }

    } else if ( is.null(bw) && degree == 0 && use.binloss) {

        K <- ncol(nn.i)

        if ( verbose ) {
            cat("Optimizing bw using cv_deg0_binloss() ")
            ptm <- proc.time()
        }

        # Define the optimization function for degree 0, binary loss
        opt.bw.fn.deg0.binloss <- function(bw) {
            bl <- 0
            out <- .C("C_cv_deg0_binloss",
                     as.integer(n.cv.folds),
                     as.integer(n.cv.reps),
                     as.integer(cv.nNN),
                     as.integer(y.binary),
                     as.integer(t(nn.i-1)),
                     as.double(t(nn.d)),
                     as.double(y),
                     as.integer(K),
                     as.integer(ng),
                     as.integer(nx),
                     as.double(bw),
                     as.integer(min.K),
                     as.integer(ikernel),
                     bl = as.double(bl))

            if ( verbose ) cat(".")

            out$bl
        }

        loop.ptm <- proc.time()
        if (n.cores > 1) {
            doParallel::registerDoParallel(n.cores)
        } else {
            foreach::registerDoSEQ()
        }
        errors <- foreach::foreach (idx = seq(n.bws),
                                    .combine = c,
                                    .export = c(
                                        "nn.i","nn.d","y","n.cv.folds","n.cv.reps","cv.nNN","y.binary",
                                        "ng","nx","min.K","ikernel","log.bws","opt.bw.fn.deg0.binloss"
                                    )
                                    ) %dopar% {
            if ( verbose ) { cat("\ri=", idx) }
            opt.bw.fn.deg0.binloss(log.bws[idx])
        }
        doParallel::stopImplicitCluster()
        if ( verbose ) {
            cat("\r\n")
            elapsed.time(loop.ptm)
        }

        opt.bw.i <- which.min(errors)
        opt.bw <- exp(log.bws[opt.bw.i])
        min.error <- errors[opt.bw.i]

        if ( verbose ) {
            cat("opt.bw.i: ", opt.bw.i, "\n")
            cat("opt.bw: ", opt.bw, "\n")
            elapsed.time(ptm)
        }

    } else if ( is.null(bw) && degree > 0 ) {

        if ( method == "LOOCV" ) {
            if ( verbose ) {
                cat("Optimizing bw using LOOCV\n")
                ptm <- proc.time()
            }

            # Define the optimization function for LOOCV
            opt.bw.fn.loocv <- function(logbw)
            {
                bw <- exp(logbw)
                max.models <- 10*ceiling( 2*bw / x.widthg )
                nn.r <- get.bws(nn.d, min.K, bw)
                rw <- row.weighting(nn.d, nn.r, nn.kernel)
                nn.w <- rw$nn.w
                max.K <- rw$max.K # max.K[i] is the index of the first 0 weight

                ## Eliminating 0 columns from nn.i, nn.d and nn.w
                cs <- colSums(nn.w)
                idx <- cs > 0
                nn.i <- nn.i[,idx]
                nn.w <- nn.w[,idx]

                ## making weights to sum to 1 row-wise
                nn.w <- nn.w / rowSums(nn.w)

                ## x and y over x NN's of xgrid
                nn.x <- row.eval(nn.i, x)
                nn.y <- row.eval(nn.i, y)

                predictions <- numeric(nx)

                out <- .C("C_loo_llm_1D",
                         as.double(x),
                         as.integer(nx),
                         as.integer(t(nn.i-1)),
                         as.double(t(nn.w)),
                         as.double(t(nn.x)),
                         as.double(t(nn.y)),
                         as.integer(ncol(nn.i)),
                         as.integer(ng),
                         as.integer(degree),
                         as.integer(max.models),
                         predictions = as.double(predictions))

                ae <- abs(y - out$predictions)

                if ( verbose ) cat(".")

                mean(ae)
            }

            loop.ptm <- proc.time()
            if (n.cores > 1) {
                doParallel::registerDoParallel(n.cores)
            } else {
                foreach::registerDoSEQ()
            }
            errors <- foreach::foreach (idx = seq(n.bws),
                                        .combine = c,
                                        .export = c(
                                            "nn.i","nn.d","y","n.cv.folds","n.cv.reps","cv.nNN","y.binary",
                                            "ng","nx","min.K","ikernel","log.bws","opt.bw.fn.loocv"
                                        )
                                        ) %dopar% {
                if ( verbose ) { cat("\ri=", idx) }
                opt.bw.fn.loocv(log.bws[idx])
            }
            doParallel::stopImplicitCluster()
            if ( verbose ) {
                cat("\r\n")
                elapsed.time(loop.ptm)
            }

            opt.bw.i <- which.min(errors)
            opt.bw <- exp(log.bws[opt.bw.i])
            min.error <- errors[opt.bw.i]

            if ( verbose ) {
                cat("opt.bw.i: ", opt.bw.i, "\n")
                cat("opt.bw: ", opt.bw, "\n")
                elapsed.time(ptm)
            }

        } else if ( method == "CV" ) {

            if ( verbose ) {
                cat("Optimizing bw using CV\n")
                ptm <- proc.time()
            }

            # Define the optimization function for CV
            opt.bw.fn.cv <- function(logbw) {

                bw <- exp(logbw)

                max.models <- 10*ceiling( 2*bw / x.widthg )
                nn.r <- get.bws(nn.d, min.K, bw)
                rw <- row.weighting(nn.d, nn.r, nn.kernel)
                nn.w <- rw$nn.w
                max.K <- rw$max.K # max.K[i] is the index of the first 0 weight

                ## Eliminating 0 columns from nn.i, nn.d and nn.w
                cs <- colSums(nn.w)
                idx <- cs > 0
                nn.i <- nn.i[,idx]
                nn.w <- nn.w[,idx]

                ## making weights to sum to 1 row-wise
                nn.w <- nn.w / rowSums(nn.w)

                ## x and y over x NN's of xgrid
                nn.x <- row.eval(nn.i, x)
                nn.y <- row.eval(nn.i, y)

                mae = 0
                K <- ncol(nn.i)

                out <- .C("C_cv_mae_1D",
                         as.integer(n.cv.folds),
                         as.integer(n.cv.reps),
                         as.double(y),
                         as.integer(t(nn.i-1)),
                         as.double(t(nn.d)),
                         as.double(t(nn.x)),
                         as.double(t(nn.y)),
                         as.integer(y.binary),
                         as.integer(ng),
                         as.integer(K),
                         as.integer(nx),
                         as.integer(degree),
                         as.double(bw),
                         as.integer(min.K),
                         as.integer(ikernel),
                         mae = as.double(mae))

                if ( verbose ) cat(".")

                out$mae
            }

            loop.ptm <- proc.time()
            if (n.cores > 1) {
                doParallel::registerDoParallel(n.cores)
            } else {
                foreach::registerDoSEQ()
            }
            errors <- foreach::foreach (idx = seq(n.bws),
                                        .combine = c,
                                        .export = c(
                                            "nn.i","nn.d","y","n.cv.folds","n.cv.reps","cv.nNN","y.binary",
                                            "ng","nx","min.K","ikernel","log.bws","opt.bw.fn.cv"
                                        )
                                        ) %dopar% {
                if ( verbose ) { cat("\ri=", idx) }
                opt.bw.fn.cv(log.bws[idx])
            }
            doParallel::stopImplicitCluster()
            if ( verbose ) {
                cat("\r\n")
                elapsed.time(loop.ptm)
            }

            opt.bw.i <- which.min(errors)
            opt.bw <- exp(log.bws[opt.bw.i])

            if ( verbose ) {
                cat("opt.bw: ", opt.bw, "\n")
                cat("opt.bw.i: ", opt.bw.i, "\n")
                elapsed.time(ptm)
            }

        } else {
            stop(paste("Unknown method:",method))
        }
    } ## END OF if ( is.null(bw) && degree == 0 && !y.binary)  ... else ...

    if ( verbose ) {
        cat("Running magelo.fit() ... ")
        ptm <- proc.time()
    }

    if (!is.null(bw)) {
        opt.bw <- bw
    }

    r <- magelo.fit(opt.bw, degree, x, y, xgrid, nn.i, nn.d, min.K = min.K, n.C.itr = n.C.itr,
                    C = C, stop.C.itr.on.min = stop.C.itr.on.min,
                    n.BB = n.BB,
                    get.BB.predictions = get.BB.predictions,
                    get.BB.gpredictions = get.BB.gpredictions,
                    get.gpredictions.CrI = get.gpredictions.CrI,
                    get.predictions.CrI = get.predictions.CrI,
                    level = level,
                    nn.kernel = nn.kernel,
                    y.binary = y.binary,
                    n.perms  =  n.perms)

    if ( verbose ) elapsed.time(ptm)

    params <- list(x = x,
                 y = y,
                 y.true = y.true,
                 grid.size = grid.size,
                 degree = degree,
                 max.K = r$max.K,
                 min.K = min.K,
                 nn.kernel = nn.kernel,
                 y.binary = y.binary,
                 f = f,
                 opt.bw = opt.bw,
                 min.bw.f = min.bw.f,
                 method = method,
                 n.bws = n.bws,
                 n.BB = n.BB,
                 get.predictions.CrI = get.predictions.CrI,
                 get.gpredictions.CrI = get.gpredictions.CrI,
                 get.BB.predictions = get.BB.predictions,
                 get.BB.gpredictions = get.BB.gpredictions,
                 level = level,
                 n.C.itr = n.C.itr,
                 C = C,
                 stop.C.itr.on.min = stop.C.itr.on.min,
                 n.cv.folds = n.cv.folds,
                 n.cv.reps = n.cv.reps,
                 cv.nNN = cv.nNN,
                 n.perms  =  n.perms,
                 n.cores = n.cores)

    output <- list(gpredictions = r$gpredictions,
                   xgrid = xgrid,
                   predictions = r$predictions,
                   n.BB = n.BB,
                   BB.predictions = r$BB.predictions,
                   BB.gpredictions = r$BB.gpredictions,
                   BB.dgpredictions = r$BB.dgpredictions,
                   gpredictions.CrI = r$gpredictions.CrI,
                   gpredictions.CrI.smooth = r$gpredictions.CrI.smooth,
                   predictions.CrI = r$predictions.CrI,
                   min.bw = min.bw,
                   max.bw = max.bw,
                   log.bws = log.bws,
                   errors = errors,
                   min.error = min.error,
                   opt.bw = opt.bw,
                   opt.bw.i = opt.bw.i,
                   beta = r$beta,
                   params = params)

    class(output) <- "magelo"

    return(output)
}

#' Plot method for magelo objects
#'
#' @param x An output from magelo()
#' @param type Type of plot: "fit", "diagnostic", "residuals" or "residuals.hist"
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param with.y.true Logical, whether to plot true values if available
#' @param with.pts Logical, whether to plot data points
#' @param with.CrI Logical, whether to plot credible intervals
#' @param true.col Color for true values
#' @param ma.col Color for magelo fit
#' @param pts.col Color for data points
#' @param with.predictions.pts Logical, whether to show prediction points
#' @param predictions.pts.col Color for prediction points
#' @param predictions.pts.pch Point character for predictions
#' @param CrI.as.polygon Logical, whether to show CrI as polygon
#' @param CrI.polygon.col Color for CrI polygon
#' @param CrI.line.col Color for CrI lines
#' @param CrI.line.lty Line type for CrI lines
#' @param ylim Y-axis limits
#' @param legend.cex Legend text size
#' @param bw Bandwidth value to use (analogous to k in malowess)
#' @param with.legend Logical, whether to show legend
#' @param legend.position Character string indicating legend position (default: "topright").
#'        One of "bottomright", "bottom", "bottomleft", "left", "topleft",
#'        "top", "topright", "right", or "center"
#' @param legend.inset Numeric value for inset distance from the margins as a fraction
#'        of the plot region (default: 0.05)
#' @param ... Additional parameters passed to plot methods
#'
#' @export
plot.magelo <- function(x, type = "fit",
                       title = "", xlab = "", ylab = "",
                       with.y.true = TRUE,
                       with.pts = FALSE,
                       with.CrI = TRUE,
                       true.col = "red",
                       ma.col = "blue",
                       pts.col = "gray60",
                       with.predictions.pts = FALSE,
                       predictions.pts.col = "blue",
                       predictions.pts.pch = 20,
                       CrI.as.polygon = TRUE,
                       CrI.polygon.col = "gray85",
                       CrI.line.col = "gray10",
                       CrI.line.lty = 2,
                       ylim = NULL,
                       legend.cex = 0.8,
                       bw = NULL,
                       with.legend = TRUE,
                       legend.position = "topright",
                       legend.inset = 0.05,
                       ...) {

    if (!inherits(x, "magelo")) {
        stop("Input must be a 'magelo' object")
    }

    # Define valid plot types
    valid_types <- c("fit", "diagnostic")

    # Find partial matches for the provided type
    type_match <- grep(paste0("^", type), valid_types, value = TRUE)

    if (length(type_match) == 0) {
        stop("Invalid plot type. Use 'fit' or 'diagnostic'")
    } else if (length(type_match) > 1) {
        stop(sprintf("Ambiguous type. '%s' matches multiple types: %s",
                    type, paste(type_match, collapse = ", ")))
    }

    type <- type_match

    # Validate bandwidth if provided
    if (!is.null(bw)) {
        if (!is.numeric(bw) || length(bw) != 1) {
            stop("bw must be a single numeric value")
        }
        if (bw < x$min.bw || bw > x$max.bw) {
            stop(sprintf("bw must be between %f and %f", x$min.bw, x$max.bw))
        }
    }

    switch(type,
           "fit" = {
               magelo.plot.fit(x, title, xlab, ylab, with.y.true,
                             with.pts, with.CrI, true.col, ma.col, pts.col,
                             with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                             CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                             CrI.line.lty, ylim, legend.cex, bw = bw,
                             with.legend = with.legend,
                             legend.position = legend.position,
                             legend.inset = legend.inset,
                             ...)
           },
           "diagnostic" = {
               magelo.plot.diagnostic(x, title, xlab, ylab,
                                    ma.col, true.col, legend.cex, ...)
           }
    )
    invisible(NULL)
}

# Helper function for fit plots
magelo.plot.fit <- function(x, title, xlab, ylab, with.y.true,
                           with.pts, with.CrI, true.col, ma.col, pts.col,
                           with.predictions.pts, predictions.pts.col, predictions.pts.pch,
                           CrI.as.polygon, CrI.polygon.col, CrI.line.col,
                           CrI.line.lty, ylim, legend.cex, bw = NULL,
                           with.legend = TRUE,
                           legend.position = "topright",
                           legend.inset = 0.05,
                           ...) {

    # Determine which predictions to use
    predictions <- x$gpredictions
    pred_label <- "Predictions"

    # Calculate y-limits if not provided
    if (is.null(ylim)) {
        ylim_data <- c(x$predictions, predictions)
        if (with.CrI && !is.null(x$gpredictions.CrI)) {
            ylim_data <- c(ylim_data, x$gpredictions.CrI[1,], x$gpredictions.CrI[2,])
        }
        ylim <- range(ylim_data, na.rm = TRUE)
    }

    # Initialize plot
    plot(x$xgrid, x$gpredictions, type = "n",
         las = 1, ylim = ylim, xlab = xlab, ylab = ylab,
         main = title, ...)

    # Add credible intervals if available and requested
    has_cri <- with.CrI && (!is.null(x$gpredictions.CrI) || !is.null(x$gpredictions.CrI.smooth))
    if (has_cri) {
        cri_data <- if(!is.null(x$gpredictions.CrI.smooth)) x$gpredictions.CrI.smooth else x$gpredictions.CrI
        if (CrI.as.polygon) {
            polygon(c(x$xgrid, rev(x$xgrid)),
                   c(cri_data[2,], rev(cri_data[1,])),
                   col = CrI.polygon.col, border = NA)
        } else {
            lines(x$xgrid, cri_data[1,], col = CrI.line.col, lty = CrI.line.lty)
            lines(x$xgrid, cri_data[2,], col = CrI.line.col, lty = CrI.line.lty)
        }
    }

    # Add predictions
    lines(x$xgrid, predictions, col = ma.col, ...)

    # Add original data points if requested
    if (with.pts) {
        points(x$params$x, x$params$y, col = pts.col)
    }

    # Add prediction points if requested
    if (with.predictions.pts) {
        points(x$xgrid, predictions,
              col = predictions.pts.col,
              pch = predictions.pts.pch)
    }

    ## Add true values if available and requested
    if (with.y.true && !is.null(x$pars$y.true) && length(x$pars$y.true) == length(x$pars$x)) {
        lines(x$pars$x, x$pars$y.true, col = true.col)
    }

    # Add legend if requested
    if (with.legend) {
        legend_items <- c()
        legend_cols <- c()
        legend_ltys <- c()
        legend_pchs <- c()

        if (with.pts) {
            legend_items <- c(legend_items, "Data Points")
            legend_cols <- c(legend_cols, pts.col)
            legend_ltys <- c(legend_ltys, NA)
            legend_pchs <- c(legend_pchs, 1)
        }

        legend_items <- c(legend_items, pred_label)
        legend_cols <- c(legend_cols, ma.col)
        legend_ltys <- c(legend_ltys, 1)
        legend_pchs <- c(legend_pchs, NA)

        if (has_cri) {
            legend_items <- c(legend_items, "95% Credible Interval")
            legend_cols <- c(legend_cols, if(CrI.as.polygon) CrI.polygon.col else CrI.line.col)
            legend_ltys <- c(legend_ltys, if(CrI.as.polygon) 1 else CrI.line.lty)
            legend_pchs <- c(legend_pchs, NA)
        }

        legend(legend.position, legend = legend_items,
               col = legend_cols, lty = legend_ltys, pch = legend_pchs,
               bg = "white", inset = legend.inset, cex = legend.cex)
    }
}

# Helper function for diagnostic plots
magelo.plot.diagnostic <- function(x, title, xlab, ylab, ma.col, true.col, legend.cex, ...) {
    # Plot errors vs bandwidth
    plot(x$log.bws, x$errors, type = 'l',
         xlab = if(xlab == "") "log(bandwidth)" else xlab,
         ylab = if(ylab == "") "Error" else ylab,
         main = if(title == "") "Error Diagnostic Plot" else title,
         col = ma.col,
         pch = 1, las = 1,
         ...)

    # Add optimal bandwidth value
    abline(v = log(x$opt.bw), col = ma.col, lty = 2)
    mtext(sprintf("Optimal bw = %.3f", x$opt.bw),
          side = 3, line = 0.25, at = log(x$opt.bw), col = ma.col)
}


#' A helper 1D local linear model routine for smoothing predictions.CI's
#'
#' @param x           A numeric vector of a predictor variable.
#' @param y           A numeric vector of an outcome variable.
#' @param grid.size   A number of grid points; was grid.size = 10*length(x), but the
#'                      results don't seem to be different from 400 which is much faster.
#' @param degree      A degree of the polynomial of x in the linear regression; 0
#'                      means weighted mean, 1 is regular linear model lm(y ~ x), and deg = d is
#'                      lm(y ~ poly(x, d)). The only allowed values are 1 and 2.
#'
#' @param f           The proportion of the range of x that is used within moving window
#'                      to train the model. It can have NULL value in which case the optimal
#'                      value of f will be found using minimum median absolute error optimization algorithm.
#' @param bw          A bandwidth parameter.
#' @param min.K       The minimal number of x NN's that must be present in each window.
#'
#' @param n.cv.folds  The number of cross-validation folds. Used only when f = NULL. Default value: 10.
#' @param n.cv.reps   The number of repetitions of cross-validation. Used only when f = NULL. Default value: 5.

#' @param nn.kernel   A kernel.
#' @param n.BB    The number of Bayesian bootstrap iterations for estimates of CI's of beta's.
#'
#' @param n.C.itr     The number of Cleveland's absolute residue based reweighting iterations for a robust estimates of mean y values.
#' @param C           A scaling of |res| parameter changing |res| to |res|/C  before applying ae.kernel to |res|'s.
#'
#' @param stop.C.itr.on.min A logical variable, if TRUE, the Cleveland's iterative reweighting stops when the maximum of the absolute values of
#'                          differences of the old and new predictions estimates are reach a local minimum.
#'

#' @return list of parameters and residues of all linear models
#'
rllmf.1D <- function(x, y, grid.size = 400, degree = 2, f = 0.3, bw = NULL,
                    n.C.itr = 0, C = 6, stop.C.itr.on.min = TRUE,
                    n.cv.folds = 10, n.cv.reps = 1,
                    n.BB = 0, min.K = 5, nn.kernel = "epanechnikov")
{
    stopifnot(is.numeric(x))
    nx <- length(x)

    stopifnot(is.numeric(y))
    stopifnot(length(y)==nx)

    stopifnot(min.K < nx)

    if ( !is.null(f) )
    {
        stopifnot(is.numeric(f))
        stopifnot(f>0 && f <=1)
    }

    stopifnot(as.integer(degree)==degree)
    stopifnot(degree==0 || degree==1 || degree==2)

    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    nn.kernel <- match.arg(nn.kernel, kernels)
    ikernel <- pmatch(nn.kernel, kernels)

    ## Sorting x and y
    o <- order(x)
    x <- x[o]
    y <- y[o]

    ng <- grid.size

    ## Defining uniform grid over the range of x values
    x.range <- range(x)
    x.width <- diff(x.range)
    xgrid <- seq(x.range[1], x.range[2], length = ng)
    x.widthg <- xgrid[2] - xgrid[1]

    ## Defining bandwidth, bw, as the half of the window size defined as f*x.width
    if ( is.null(bw) )
    {
        window.size <- f * x.width
        bw <- window.size / 2
    }

    ## Identifying x NN's from xgrid
    ##
    ## They are used for training of linear models and prediction of predictions over x
    nn <- get.knnx(x, xgrid, k = nx)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist
    init.nn.i <- nn.i
    init.nn.d <- nn.d

    r <- magelo.fit(bw, degree, x, y, xgrid, nn.i, nn.d, min.K = min.K, n.C.itr = n.C.itr,
                    C = C, stop.C.itr.on.min = stop.C.itr.on.min, n.BB = n.BB,
                    nn.kernel = nn.kernel)

    list(x = x,
         y = y,
         xgrid = xgrid,
         gpredictions = r$gpredictions,
         predictions = r$predictions)
}

#' Core Fitting Function for MAGELO (Model Averaged Grid-based Epsilon LOwess)
#'
#' @description
#' Internal fitting function for MAGELO that performs local polynomial regression using
#' pre-computed nearest neighbor information. Supports weighted mean (degree 0), local linear
#' (degree 1), and local quadratic (degree 2) regression with optional robust fitting via
#' iterative reweighting and uncertainty estimation via Bayesian bootstrap.
#'
#' @param bw Numeric specifying bandwidth (radius of disk neighborhood)
#' @param degree Integer (0, 1, or 2) specifying polynomial degree:
#'        - 0: weighted mean regression
#'        - 1: local linear regression (y ~ x)
#'        - 2: local quadratic regression (y ~ poly(x, 2))
#' @param x Numeric vector of predictor values
#' @param y Numeric vector of response values, must be same length as x
#' @param xgrid Numeric vector of uniformly spaced grid points spanning range(x)
#' @param nn.i Integer matrix (n_grid x K) of indices for x neighbors of each grid point
#' @param nn.d Numeric matrix (n_grid x K) of distances to x neighbors of each grid point
#' @param min.K Integer specifying minimum number of points required in each local neighborhood
#' @param n.C.itr Integer specifying number of Cleveland's iterative reweighting steps
#' @param C Numeric scaling factor for residuals in robust fitting.
#'        Higher values reduce influence of outliers
#' @param stop.C.itr.on.min Logical; if TRUE, stop iterations when improvements in
#'        fitted values become negligible
#' @param n.BB Integer specifying number of Bayesian bootstrap iterations
#' @param get.BB.gpredictions Logical; if TRUE, return matrix of bootstrap estimates at grid points
#' @param get.BB.predictions Logical; if TRUE, return matrix of bootstrap estimates at x
#' @param get.gpredictions.CrI Logical; if TRUE, compute credible intervals at grid points
#' @param get.predictions.CrI Logical; if TRUE, compute credible intervals at x
#' @param level Numeric between 0 and 1 specifying credible interval level
#' @param nn.kernel Character string specifying kernel function:
#'        "epanechnikov", "triangular", "tr.exponential", or "normal"
#' @param y.binary Logical; if TRUE, treat y as binary outcome
#' @param n.perms Integer specifying number of y permutations for p-value calculation
#'
#' @section Robust Fitting Parameters:
#' The robust fitting procedure implements Cleveland's iterative reweighting algorithm to
#' reduce the influence of outliers. The \code{n.C.itr} parameter controls the maximum number
#' of reweighting iterations, while \code{C} determines how aggressively outliers are
#' downweighted. Setting \code{stop.C.itr.on.min = TRUE} enables early stopping when the
#' algorithm converges, saving computation time without sacrificing accuracy.
#'
#' @section Uncertainty Estimation:
#' Uncertainty quantification is performed using Bayesian bootstrap, which generates posterior
#' distributions for fitted values by resampling with exponential weights. The \code{n.BB}
#' parameter controls the number of bootstrap samples. Credible intervals can be computed at
#' both the original data points and grid points. The bootstrap estimates themselves can also
#' be returned for custom uncertainty analyses.
#'
#' @return A list containing:
#' \itemize{
#'   \item x: Input predictor values
#'   \item y: Input response values
#'   \item xgrid: Grid points
#'   \item max.K: Maximum number of neighbors used for each point
#'   \item gpredictions: Fitted values at grid points
#'   \item predictions: Fitted values at x
#'   \item BB.predictions: Bootstrap estimates at x (if requested)
#'   \item BB.gpredictions: Bootstrap estimates at grid points (if requested)
#'   \item BB.dgpredictions: Bootstrap derivatives at grid points (if requested)
#'   \item gpredictions.CrI: Credible intervals at grid points (if requested)
#'   \item predictions.CrI: Credible intervals at x (if requested)
#'   \item beta: Model coefficients at each grid point
#' }
#'
#' @note
#' This is an internal function called by magelo() after bandwidth optimization.
#' It should not typically be called directly unless you have pre-computed nearest
#' neighbor information and know the appropriate bandwidth.
magelo.fit <- function(bw,
                       degree,
                       x,
                       y,
                       xgrid,
                       nn.i,
                       nn.d,
                       min.K = 15,
                       n.C.itr = 0,
                       C = 6,
                       stop.C.itr.on.min = TRUE,
                       n.BB = 0,
                       get.BB.gpredictions = FALSE,
                       get.BB.predictions = FALSE,
                       get.gpredictions.CrI = FALSE,
                       get.predictions.CrI = TRUE,
                       level = 0.95,
                       nn.kernel = "epanechnikov",
                       y.binary = FALSE,
                       n.perms  =  0) {
    nx <- length(x)
    ng <- nrow(nn.i)

    nn.r <- get.bws(nn.d, min.K, bw)

    rw <- row.weighting(nn.d, nn.r, nn.kernel)
    nn.w <- rw$nn.w
    max.K <- rw$max.K # max.K[i] is the index of the last non-zero weight

    ## Eliminating 0 columns from nn.i, nn.d and nn.w
    cs <- colSums(nn.w)
    idx <- cs > 0
    nn.i <- nn.i[,idx]
    nn.d <- nn.d[,idx]
    nn.w <- nn.w[,idx]

    ## x and y over x NN's of xgrid
    nn.x <- row.eval(nn.i, x)
    nn.y <- row.eval(nn.i, y)

    ## Local linear models and corresonding wmean gpredictions estimates
    beta <- NA
    BB.gpredictions <- NULL
    BB.dgpredictions <- NULL
    BB.predictions <- NULL
    gpredictions.CrI <- NULL
    predictions.CrI <- NULL
    gpredictions.CrI.smooth <- NULL

    if ( degree==0 ) # weighed mean regression
    {
        ## making weights to sum to 1 row-wise
        nn.w <- nn.w / rowSums(nn.w)

        gpredictions <- row.weighted.mean(nn.y, nn.w, max.K)
        predictions <- approx(xgrid, gpredictions, xout = x)$y

        if ( n.C.itr > 0 ) {
            predictions <- Clevelands.1D.deg0.loop(n.C.itr, x, y, xgrid, predictions, nn.i, nn.w, nn.y, max.K, C, stop.C.itr.on.min)
        }

        if ( n.BB > 0 && get.predictions.CrI ) {
            gpredictions.CrI <- row.weighted.mean.BB.qCrI(y.binary, nn.y, nn.w, max.K, n.BB, 1 - level)

            ## smoothing gpredictions.CrI
            rr <- rllmf.1D(xgrid, gpredictions.CrI[1, ], degree = 0, bw = bw/3)
            gpredictions.CrI.smooth.u <- rr$gpredictions

            rr <- rllmf.1D(xgrid, gpredictions.CrI[2, ], degree = 0, bw = bw/3)
            gpredictions.CrI.smooth.d <- rr$gpredictions

            gpredictions.CrI.smooth <- rbind(gpredictions.CrI.smooth.u, gpredictions.CrI.smooth.d)

            predictions.CrI.u <- approx(xgrid, gpredictions.CrI.smooth[1, ], xout = x)$y
            predictions.CrI.l <- approx(xgrid, gpredictions.CrI.smooth[2, ], xout = x)$y

            predictions.CrI <- rbind(predictions.CrI.u, predictions.CrI.l)

        } else if ( n.BB > 0 && !get.predictions.CrI ) {

            gpredictions.mad.CI <- row.weighted.mean.BB.CI.v2(gpredictions, nn.y, nn.w, max.K, n.BB = n.BB)
            gpredictions.CrI <- rbind(gpredictions + gpredictions.mad.CI, gpredictions - gpredictions.mad.CI)

            ## smoothing gpredictions.CrI
            rr <- rllmf.1D(xgrid, gpredictions.mad.CI, degree = 0, bw = bw)
            gpredictions.mad.CI.smooth <- rr$gpredictions
            gpredictions.CrI.smooth <- rbind(gpredictions + gpredictions.mad.CI.smooth, gpredictions - gpredictions.mad.CI.smooth)

            predictions.mad.CI <- approx(xgrid, gpredictions.mad.CI.smooth, xout = x)$y
            predictions.CrI <- rbind(predictions + predictions.mad.CI, predictions - predictions.mad.CI)
        }

    } else { # degrees 2 and 3

        beta <- llm.1D.beta(nn.x, nn.y, nn.w, max.K, degree)

        ## NOTE: normalization of nn.w is done within llm_1D_beta() function
        ## fp <- llm.1D.fit.and.predict(nn.i, nn.w, nn.x, nn.y, max.K, degree, nx)
        ## beta <- fp$beta
        ## predictions <- fp$predictions

        ## packaging variables of the model for prediction of predictions over xgrid
        r <- list(beta = beta,
                 y.binary = y.binary,
                 xgrid = xgrid,
                 ng = ng,
                 min.K = min.K,
                 bw = bw,
                 nn.r = nn.r,
                 max.K = max.K,
                 degree = degree,
                 nn.kernel = nn.kernel)

        class(r) <- "magelo"

        ##r.pred <- pred.llm.1D(r, new.x = xgrid)
        r.pred <- predict(r, newdata = xgrid)

        gpredictions <- r.pred$predictions
        grid.nn.i <- r.pred$nn.i
        grid.nn.w <- r.pred$nn.w
        grid.nn.x <- r.pred$nn.x
        n.grid <- r.pred$n
        grid.max.K <- r.pred$max.K

        predictions <- approx(xgrid, gpredictions, xout = x)$y

        if ( n.C.itr > 0 )
        {
            predictions <- Clevelands.1D.deg12.loop(n.C.itr, degree, nx, y, predictions, nn.i, nn.w, nn.x, nn.y, max.K, C, stop.C.itr.on.min)
        }

        if ( n.BB > 0 )
        {

            if ( get.BB.gpredictions || get.BB.predictions )
            {
                BB.res <- get.BB.gpredictions(n.BB,
                                    nn.i,
                                    nn.x,
                                    nn.y,
                                    y.binary,
                                    nn.w,
                                    nx,
                                    max.K,
                                    degree,
                                    grid.nn.i,
                                    grid.nn.x,
                                    grid.nn.w,
                                    grid.max.K)

                BB.gpredictions <- BB.res$bb.gpredictions
                BB.dgpredictions <- BB.res$bb.dgpredictions

                if ( get.BB.predictions ) {
                    BB.predictions <- apply(BB.gpredictions, 2, function(z) approx(xgrid, z, xout = x)$y )
                }

                ## R implementation
                ## loop over n.BB iterations
                ## at each iterations
                ## 1) select random sample, lambda, from a simple of dim nx-1
                ## 2) change nn.w by lambda
                ## 3) fit local linear models at the grid points using the modified nn.w's
                ## 4) use the resulting model coefficients to predict bb.gpredictions (BB version of gpredictions)
            }

            if ( get.gpredictions.CrI || get.predictions.CrI )
            {
                gpredictions.CrI <- get.gpredictions.CrI(n.BB,
                                      nn.i,
                                      nn.x,
                                      nn.y,
                                      nn.w,
                                      nx,
                                      max.K,
                                      degree,
                                      grid.nn.i,
                                      grid.nn.x,
                                      grid.nn.w,
                                      grid.max.K,
                                      y.binary,
                                      1 - level)

                if ( get.predictions.CrI )
                {
                    predictions.CrI <- matrix(nrow = 2, ncol = nx)
                    predictions.CrI[1,] <- approx(xgrid, gpredictions.CrI[1,], xout = x)$y
                    predictions.CrI[2,] <- approx(xgrid, gpredictions.CrI[2,], xout = x)$y
                }
            }
        }

        if ( n.perms > 0 ) {

            beta.perms <- llm.1D.beta.perms(nn.x, nn.i, y, nn.w, max.K, degree, n.perms)
        }
    }

    list(x = x,
         y = y,
         xgrid = xgrid,
         max.K = max.K,
         gpredictions = gpredictions,
         predictions = predictions,
         BB.predictions = BB.predictions,
         BB.gpredictions = BB.gpredictions,
         BB.dgpredictions = BB.dgpredictions,
         gpredictions.CrI = gpredictions.CrI,
         gpredictions.CrI.smooth = gpredictions.CrI.smooth,
         predictions.CrI = predictions.CrI,
         nn.kernel = nn.kernel,
         min.K = min.K,
         bw = bw,
         beta = beta)
}

#' Eestimates the coefficients of linear models in the vicinity of each grid point of radius nn.r.
#'
#' The max.K parameter is used to avoid looping over indices where weights are 0.
#'
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param nn.w      A matrix of weights over K nearest neighbors of each element of the grid.
#' @param max.K      A vector of indices indicating the range where weights are not 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#'
#' @return list with components: beta, gpredictions.
#'
#' NOTE: Weights must sum up to 1!
#'
llm.1D.beta <- function(nn.x, nn.y, nn.w, max.K, degree) {

    ng <- nrow(nn.x) # number of grid points = ncol(t(nn's))
    K <- ncol(nn.x)  # number of NN's = nrow(t(nn.'s))

    ncoef <- degree + 1

    beta <- matrix(0, nrow = ng, ncol = ncoef)

    out <- .C("C_llm_1D_beta",
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(degree),
             beta = as.double(t(beta)))

    beta <- matrix(out$beta, nrow = ng, ncol = ncoef, byrow = TRUE)

    beta
}

#' Estimates the coefficients of linear models in the vicinity of each grid point of radius nn.r on the permuted version
#'
#' The max.K parameter is used to avoid looping over indices where weights are 0.
#'
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.i      A matrix of indices of the K nearest neighbors of each element of the grid.
#' @param y         A vector of response values.
#' @param nn.w      A matrix of weights over K nearest neighbors of each element of the grid.
#' @param max.K     A vector of indices indicating the range where weights are not 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param n.perms   Number of permutations to perform. Must be a positive integer.
#'
#' @return A numeric vector of beta coefficients from all permutations.
#'
#' NOTE: Weights must sum up to 1!
#'
llm.1D.beta.perms <- function(nn.x, nn.i, y, nn.w, max.K, degree, n.perms) {
    if ( n.perms <= 0 ) {
        stop("n.perms has to be positive integer; n.perms:", n.perms)
    }
    ng <- nrow(nn.x) # number of grid points = ncol(t(nn's))
    K <- ncol(nn.x)  # number of NN's = nrow(t(nn.'s))
    ncoef <- degree + 1
    beta <- numeric(ng * ncoef * n.perms)
    out <- .C("C_llm_1D_beta_perms",
             as.double(t(nn.x)),
             as.integer(t(nn.i)),
             as.double(y),
             as.integer(length(y)),
             as.double(t(nn.w)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(degree),
             as.integer(n.perms),
             beta = as.double(t(beta)))
    out$beta
}

#' Compute predictions from local linear model coefficients
#'
#' @param beta       Coefficients of the local linear model.
#' @param nn.i       A matrix of indices of x NN's of xgrid.
#' @param nn.d       A matrix of NN distances.
#' @param nn.x       The values of predictor variable over NN's of xgrid.
#' @param nx         The length of x.
#' @param max.K      A vector of indices indicating the range where weights are not 0.
#' @param nn.kernel  A kernel used for generating weights.
#'
#' @return A numeric vector of predictions
#'
llm.predict.1D <- function(beta, nn.i, nn.d, nn.x, nx, max.K, nn.kernel = "epanechnikov") {
    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    ikernel <- pmatch(nn.kernel, kernels)

    ng <- nrow(beta) # number of grid points = ncol(t(nn's))
    K <- ncol(nn.x)  # number of NN's = nrow(t(nn.'s))
    degree <- ncol(beta) - 1

    predictions <- numeric(nx)

    out <- .C("C_predict_1D",
             as.double(t(beta)),
             as.integer(t(nn.i-1)),
             as.double(t(nn.d)),
             as.double(t(nn.x)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(degree),
             as.integer(ikernel),
             as.integer(nx),
             predictions = as.double(predictions))

    out$predictions
}

#' Predicts the mean values of linear models, beta, over grid points of x
#'
#' It is a variation on predict.1D() with nn.w passed to it instead of nn.d
#' and nn.kernel for calculation of nn.w inside the corresponding C routine.
#'
#' @param beta       Coefficients of the local linear model.
#' @param nn.i       A matrix of indices of x NN's of xgrid.
#' @param nn.w       A matrix of NN weights.
#' @param nn.x       The values of predictor variable over NN's of xgrid.
#' @param nx         The length of x.
#' @param max.K      A vector of __indices__ indicating the range where weights are not 0.
#' @param y.binary   Set to TRUE if y is a binary variable.
wpredict.1D <- function(beta, nn.i, nn.w, nn.x, nx, max.K, y.binary = FALSE)
{
    ng <- nrow(beta) # number of grid points = ncol(t(nn's))
    K <- ncol(nn.x)  # number of NN's = nrow(t(nn.'s))

    degree <- ncol(beta) - 1

    predictions <- numeric(nx)

    out <- .C("C_wpredict_1D",
             as.double(t(beta)),
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.integer(max.K-1),
             as.integer(K),
             as.integer(ng),
             as.integer(degree),
             as.integer(nx),
             as.integer(y.binary),
             predictions = as.double(predictions))

    out$predictions
}

#' Predicting values from a magelo model
#'
#' @param object   an object of class "magelo", typically the result of a call to \code{\link{magelo}}
#' @param newdata  a numeric vector of values at which predictions are required
#' @param ...      additional arguments (currently unused)
#' @return numeric vector of predicted values
#' @method predict magelo
#' @examples
#' \dontrun{
#' model <- magelo(x, y)
#' predictions <- predict(model, newdata = c(1, 2, 3))
#' }
#' @export
predict.magelo <- function(object, newdata, ...)
{
    if ( !is.vector(newdata) )
    {
        stop("newdata has to be a vector")
    }

    y.binary <- object$y.binary
    n.newdata <- length(newdata)
    nn.i <- NULL
    nn.w <- NULL
    nn.x <- NULL
    max.K <- NULL

    if ( object$degree > 0 )
    {
        ## newdata NN's from xgrid
        nn <- get.knnx(newdata, object$xgrid, k = n.newdata)
        nn.i  <- nn$nn.index
        nn.d  <- nn$nn.dist

        rw <- row.weighting(nn.d, object$nn.r, object$nn.kernel)
        nn.w <- rw$nn.w
        max.K <- rw$max.K

        ## Eliminating 0 columns from nn.i, nn.d, nn.w
        if ( ncol(nn.w) > 1 ) {
            cs <- colSums(nn.w)
            idx <- cs > 0
            nn.i <- nn.i[,idx]
            nn.d <- nn.d[,idx]
            nn.w <- nn.w[,idx]
        }

        if ( is.null(ncol(nn.i)) ) {
            nn.i <- cbind(nn.i)
            nn.d <- cbind(nn.d)
            nn.w <- cbind(nn.w)
        }

        if ( ncol(nn.w) > 1 ) {
            nn.w <- row.TS.norm( nn.w )
        }

        ## newdata NN's of xgrid
        nn.x <- row.eval(nn.i, newdata)

        if ( is.null(ncol(nn.x)) )
        {
            nn.x <- cbind(nn.x)
        }

        predictions.newdata <- wpredict.1D(object$beta, nn.i, nn.w, nn.x, n.newdata, max.K, y.binary)

    } else {
        predictions.newdata <- approx(object$xgrid, object$gpredictions, xout = newdata, rule = 2)$y
    }

    if ( y.binary ) {
        predictions.newdata <- ifelse(predictions.newdata<0, 0, predictions.newdata)
        predictions.newdata <- ifelse(predictions.newdata>1, 1, predictions.newdata)
    }

    list(predictions = predictions.newdata,
         nn.i = nn.i,
         nn.w = nn.w,
         nn.x = nn.x,
         n = n.newdata,
         max.K = max.K)
}

#' Leave-One-Out (LOO) cross-validation of 1D local linear models
#'
#' @param x           A numeric vector of a predictor variable.
#' @param y           A numeric vector of an outcome variable.
#' @param grid.size   A number of grid points; was grid.size = 10*length(x), but the
#'                      results don't seem to be different from 400 which is much faster.
#' @param degree      A degree of the polynomial of x in the linear regression; 0
#'                      means weighted mean, 1 is regular linear model lm(y ~ x), and deg = d is
#'                      lm(y ~ poly(x, d)). The only allowed values are 1 and 2.
#' @param f           The proportion of the range of x that is used within moving window to train the model.
#' @param bw          A bandwidth parameter.
#' @param min.K       The minimal number of x NN's that must be present in each window.
#' @param nn.kernel  A kernel.
#'
#' @return list of parameters and residues of all linear models
#'
#' @export
loo.llm.1D <- function(x, y, grid.size = 400, degree = 2, f = 0.2, bw = NULL, min.K = 5, nn.kernel = "epanechnikov")
{
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))

    nx <- length(x)
    stopifnot(length(y) == nx)

    stopifnot(min.K < nx)

    stopifnot(is.numeric(f))
    stopifnot(f>0 && f <=1)

    stopifnot(as.integer(degree)==degree)
    stopifnot(degree==0 || degree==1 || degree==2)

    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    nn.kernel <- match.arg(nn.kernel, kernels)

    ## Sorting x and y
    o <- order(x)
    x <- x[o]
    y <- y[o]

    ## Defining uniform grid over the range of x values
    ng <- grid.size
    x.range <- range(x)
    dx <- diff(x.range)
    xgrid <- seq(x.range[1], x.range[2], length = ng)
    dxgrid <- xgrid[2] - xgrid[1]

    ## Defining bandwidth, bw, as the half of the window size defined as f*dx
    if ( is.null(bw) )
    {
        window.size <- f * dx
        bw <- window.size / 2
    }

    ##
    ## Identifying x NN's from xgrid
    ##
    ## They are used for training of linear models and prediction of predictions over x
    ##
    nn <- get.knnx(x, xgrid, k = nx-1)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    ## Making sure bw is not too small. That is resulting in window
    ## sizes that have too few points to train linear modela.
    nn.r <- get.bws(nn.d, min.K, bw)
    rw <- row.weighting(nn.d, nn.r, nn.kernel)
    nn.w <- rw$nn.w
    max.K <- rw$max.K

    ## Eliminating 0 columns from nn.i, nn.d and nn.w
    cs <- colSums(nn.w)
    idx <- cs > 0
    nn.i <- nn.i[,idx]
    nn.d <- nn.d[,idx]
    nn.w <- nn.w[,idx]

    ## making weights to sum to 1 row-wise
    nn.w <- nn.w / rowSums(nn.w)

    ## x and y over x NN's of xgrid
    nn.x <- row.eval(nn.i, x)
    nn.y <- row.eval(nn.i, y)

    ## An estimate of the maximum of the number of models associated with a
    ## single value of x.
    max.models <- ceiling( 2*bw / dxgrid )

    K <- ncol(nn.i) # number of NN's = nrow(t(nn.'s))

    predictions <- numeric(nx)

    out <- .C("C_loo_llm_1D",
             as.double(x),
             as.integer(nx),
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(K),
             as.integer(ng),
             as.integer(degree),
             as.integer(max.models),
             predictions = as.double(predictions))

    out$predictions
}

#' Degree 0 case of Leave-One-Out (LOO) cross-validation of 1D local linear models
#'
#' @param nx         The number of elements of x.
#' @param nn.i       A matrix of indices of x NN's of xgrid.
#' @param nn.w       A matrix of NN weights.
#' @param nn.y       The values of y over nn.i.
#'
#' @return predictions
#'
deg0.loo.llm.1D <- function(nx, nn.i, nn.w, nn.y)
{
    predictions <- numeric(nx)

    out <- .C("C_deg0_loo_llm_1D",
              as.integer(nx),
              as.integer(t(nn.i-1)),
              as.double(t(nn.w)),
              as.double(t(nn.y)),
              as.integer(ncol(nn.i)),
              as.integer(nrow(nn.i)),
              predictions = as.double(predictions))

    out$predictions
}

#' Fits 1D rllm model and generates predictions estimates.
#'
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param y.binary  Set to TRUE if y's values are within the interval \code{[0,1]}.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param nx        The number of elements of x.
#'
#' @return list with components: beta, predictions.
#'
#' NOTE: Weights must sum up to 1!
#'
llm.1D.fit.and.predict <- function(nn.i, nn.w, nn.x, nn.y, y.binary, max.K, degree, nx)
{
    ## stopifnot(is.finite(nn.i))
    ## stopifnot(is.finite(nn.y))
    ## stopifnot(is.finite(nn.w))

    ng <- nrow(nn.i) # number of grid points
    ncoef <- degree + 1

    beta <- matrix(0, nrow = ng, ncol = ncoef)
    predictions <- numeric(nx)

    out <- .C("C_llm_1D_fit_and_predict",
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(y.binary),
             as.integer(max.K-1),
             as.integer(ncol(nn.i)),
             as.integer(ng),
             as.integer(nx),
             as.integer(degree),
             predictions = as.double(predictions),
             beta = as.double(t(beta)))

    beta <- matrix(out$beta, nrow = ng, ncol = ncoef, byrow = TRUE)

    list(beta = beta,
         predictions = out$predictions)
}

#' Creates BB CI's of a 1D rllm model's predictions estimates.
#'
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param y.binary  Set to TRUE if y values are in the interval \code{[0,1]}.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param nx        The number of elements of x.
#' @param predictions        An predictions estimates.
#' @param n.BB      The number of BB iterations.
#'
#' @return predictions.CrI of length nx.
#'
llm.1D.fit.and.predict.BB.CrI <- function(nn.i, nn.w, nn.x, nn.y, y.binary, max.K, degree, nx, predictions, n.BB)
{
    ## stopifnot(is.finite(nn.i))
    ## stopifnot(is.finite(nn.y))
    ## stopifnot(is.finite(nn.w))

    ng <- nrow(nn.i) # number of grid points

    ncoef <- degree + 1

    predictions.CrI <- numeric(nx)

    out <- .C("C_llm_1D_fit_and_predict_BB_CrI",
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(y.binary),
             as.integer(max.K-1),
             as.integer(ncol(nn.i)),
             as.integer(ng),
             as.integer(nx),
             as.integer(degree),
             as.integer(n.BB),
             as.double(predictions),
             predictions.CrI = as.double(predictions.CrI))

    return( out$predictions.CrI )
}

#' Creates BB CI's of a 1D rllm model's predictions estimates using global reweighting of the elements of x.
#'
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param nx        The number of elements of x.
#' @param predictions        An predictions estimates.
#' @param y.binary  Set to TRUE if y is a binary variable.
#' @param n.BB      The number of BB iterations.
#'
#' @return predictions.CrI of length nx.
#'
llm.1D.fit.and.predict.global.BB.CrI <- function(nn.i, nn.w, nn.x, nn.y, max.K, degree, nx, predictions, y.binary, n.BB)
{
    ## stopifnot(is.finite(nn.i))
    ## stopifnot(is.finite(nn.y))
    ## stopifnot(is.finite(nn.w))

    ng <- nrow(nn.i) # number of grid points

    ncoef <- degree + 1

    predictions.CrI <- numeric(nx)

    out <- .C("C_llm_1D_fit_and_predict_global_BB_CrI",
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(y.binary),
             as.integer(max.K-1),
             as.integer(ncol(nn.i)),
             as.integer(ng),
             as.integer(nx),
             as.integer(degree),
             as.integer(n.BB),
             as.double(predictions),
             predictions.CrI = as.double(predictions.CrI))

    return( out$predictions.CrI )
}

#' Creates Bayesian bootstrap (BB) estimates of predictions using global reweighting of the elements of x.
#'
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param nx        The number of elements of x.
#' @param n.BB      The number of BB iterations.
#'
#' @return BB predictions matrix of dim nx-by-n.BB.
#'
llm.1D.fit.and.predict.global.BB <- function(nn.i, nn.w, nn.x, nn.y, max.K, degree, nx, n.BB)
{
    ng <- nrow(nn.i) # number of grid points

    ncoef <- degree + 1

    bbpredictions <- numeric(nx * n.BB)

    out <- .C("C_llm_1D_fit_and_predict_global_BB",
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(max.K-1),
             as.integer(ncol(nn.i)),
             as.integer(ng),
             as.integer(nx),
             as.integer(degree),
             as.integer(n.BB),
             bbpredictions = as.double(bbpredictions))

    bbpredictions <- matrix(out$bbpredictions, nrow = nx, ncol = n.BB, byrow = FALSE)

    return( bbpredictions )
}

#' Creates BB quantile-based estimates of predictions CI's.
#'
#' @param y.binary  A logical variable. If TRUE, predicted predictions will be trimmed to the closed interval \code{[0, 1]}.
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param nn.y      A matrix of y values over K nearest neighbors of each element of the grid.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#' @param nx        The number of elements of x.
#' @param n.BB      The number of BB iterations.
#' @param alpha     The confidence level.
#'
#' @return predictions.qCI matrix with two rows (upper and lower CI) and nx columns.
#'
llm.1D.fit.and.predict.global.BB.qCrI <- function(y.binary, nn.i, nn.w, nn.x, nn.y, max.K, degree, nx, n.BB, alpha = 0.05)
{
    stopifnot(is.logical(y.binary))

    ng <- nrow(nn.i) # number of grid points

    ncoef <- degree + 1

    predictions.qCI <- numeric(nx * 2)

    out <- .C("C_llm_1D_fit_and_predict_global_BB_qCrI",
             as.integer(y.binary),
             as.integer(t(nn.i-1)),
             as.double(t(nn.w)),
             as.double(t(nn.x)),
             as.double(t(nn.y)),
             as.integer(max.K-1),
             as.integer(ncol(nn.i)),
             as.integer(ng),
             as.integer(nx),
             as.integer(degree),
             as.integer(n.BB),
             as.double(alpha),
             predictions.qCI = as.double(predictions.qCI))

    predictions.qCI <- matrix(out$predictions.qCI, nrow = 2, ncol = nx, byrow = FALSE)

    return(predictions.qCI)
}

#' Cleveland's iterative reweighting for 1D degree 0 models
#'
#' Internal function that performs Cleveland's iterative reweighting on 1D degree 0 models
#' to obtain robust estimates by downweighting outliers.
#'
#' @param n.C.itr The number of Cleveland's weights updating iterations.
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of response values.
#' @param xgrid A uniform grid within the range of x.
#' @param predictions Current estimate of predictions.
#' @param nn.i A matrix of indices of x nearest neighbors of xgrid.
#' @param nn.w A matrix of weights of x nearest neighbors of xgrid.
#' @param nn.y A matrix of y values over nn.i.
#' @param max.K A vector of indices indicating the range where weights are not 0.
#' @param C Cleveland's weights updating scaling factor. Default is 6.
#' @param stop.C.itr.on.min Logical; if TRUE, stops iteration when improvement plateaus.
#'
#' @return Updated predictions vector
#'
#' @details
#' By default, this function attempts to use Gaussian mixture modeling via
#' \pkg{mclust} (function \code{Mclust}) to identify the main component of the
#' residual distribution and estimate a robust standard deviation. If the
#' \pkg{mclust} package is not available, it falls back to using the median
#' absolute deviation (\code{\link[stats]{mad}}). In both cases, observations
#' are then iteratively reweighted based on their residuals.
#'
#' @keywords internal
#' @noRd
Clevelands.1D.deg0.loop <- function(n.C.itr, x, y, xgrid, predictions, nn.i, nn.w,
                                   nn.y, max.K, C = 6, stop.C.itr.on.min = TRUE) {
    ae.kernel <- "normal"

    res <- y - predictions

    ## Robust standard deviation estimate
    if (requireNamespace("mclust", quietly = TRUE)) {
        ## Use Gaussian mixture model to estimate robust standard deviation
        m <- mclust::Mclust(res, modelNames = "V", verbose = FALSE)
        i.max <- which.max(m$parameters$pro)
        sigma <- sqrt(m$parameters$variance$sigmasq[i.max])
    } else {
        warning("Package 'mclust' not installed; falling back to MAD-based estimate.")
        sigma <- mad(res, constant = 1)  # constant=1 so it's comparable to sd()
    }

    itr <- 1
    old.Dpredictions <- 1000
    predictions.old <- predictions

    while (itr < n.C.itr) {
        # Calculate weights based on absolute errors
        ae <- abs(res) / (C * sigma)
        ae.weights <- kernel.fn(ae, ae.kernel)
        ae.weights <- ae.weights / sum(ae.weights)

        # Update neighbor weights
        for (j in seq(nrow(nn.w))) {
            nn.w[j,] <- nn.w[j,] * ae.weights[nn.i[j,]]
            s <- sum(nn.w[j,])
            if (s > 0) {
                nn.w[j,] <- nn.w[j,] / s
            }
        }

        # Compute new predictions
        gpredictions <- row.weighted.mean(nn.y, nn.w, max.K)
        predictions <- approx(xgrid, gpredictions, xout = x)$y

        # Check for convergence
        Dpredictions <- max(abs(predictions - predictions.old))
        if (!is.finite(Dpredictions)) {
            stop("Non-finite change in predictions")
        }

        # Stop if improvement has plateaued
        if (stop.C.itr.on.min && Dpredictions > old.Dpredictions) {
            predictions <- predictions.old
            break
        }

        # Update for next iteration
        res <- y - predictions
        predictions.old <- predictions
        old.Dpredictions <- Dpredictions
        itr <- itr + 1
    }

    predictions
}

#' Cleveland's iterative reweighting for 1D degree 1 and 2 models
#'
#' Internal function that performs Cleveland's iterative reweighting on 1D degree 1 and 2
#' models to obtain robust estimates by downweighting outliers.
#'
#' @param n.C.itr The number of Cleveland's weights updating iterations.
#' @param degree The degree of the polynomial model (1 or 2).
#' @param nx The number of elements of x and y.
#' @param y A numeric vector of response values.
#' @param predictions Current estimate of predictions.
#' @param nn.i A matrix of indices of x nearest neighbors of xgrid.
#' @param nn.w A matrix of weights of x nearest neighbors of xgrid.
#' @param nn.x A matrix of x values over nn.i.
#' @param nn.y A matrix of y values over nn.i.
#' @param max.K A vector of indices indicating the range where weights are not 0.
#' @param C Cleveland's weights updating scaling factor. Default is 6.
#' @param stop.C.itr.on.min Logical; if TRUE, stops iteration when improvement plateaus.
#'
#' @return Updated predictions vector
#'
#' @details
#' By default, this function attempts to use Gaussian mixture modeling via
#' \pkg{mclust} (function \code{Mclust}) to identify the main component of the
#' residual distribution and estimate a robust standard deviation. If the
#' \pkg{mclust} package is not available, it falls back to using the median
#' absolute deviation (\code{\link[stats]{mad}}). In both cases, observations
#' are then iteratively reweighted based on their residuals.
#'
#' @keywords internal
#' @noRd
Clevelands.1D.deg12.loop <- function(n.C.itr, degree, nx, y, predictions, nn.i, nn.w,
                                    nn.x, nn.y, max.K, C = 6, stop.C.itr.on.min = TRUE) {
    ae.kernel <- "normal"

    res <- y - predictions

    ## Robust standard deviation estimate
    if (requireNamespace("mclust", quietly = TRUE)) {
        ## Use Gaussian mixture model to estimate robust standard deviation
        m <- mclust::Mclust(res, modelNames = "V", verbose = FALSE)
        i.max <- which.max(m$parameters$pro)
        sigma <- sqrt(m$parameters$variance$sigmasq[i.max])
    } else {
        warning("Package 'mclust' not installed; falling back to MAD-based estimate.")
        sigma <- mad(res, constant = 1)  # constant=1 so it's comparable to sd()
    }

    # Initialize for iterations
    itr <- 1
    old.Dpredictions <- 1000
    predictions.old <- predictions

    # Store original weights for final prediction
    orig.nn.w <- nn.w

    # Initial beta calculation (note: beta.old was used before being defined)
    beta <- llm.1D.beta(nn.x, nn.y, nn.w, max.K, degree)
    beta.old <- beta

    while (itr < n.C.itr) {
        # Calculate weights based on absolute errors
        ae <- abs(res) / (C * sigma)
        ae.weights <- kernel.fn(ae, ae.kernel)
        ae.weights <- ae.weights / sum(ae.weights)

        # Update neighbor weights
        for (j in seq(nrow(nn.w))) {
            nn.w[j,] <- nn.w[j,] * ae.weights[nn.i[j,]]
            s <- sum(nn.w[j,])
            if (s > 0) {
                nn.w[j,] <- nn.w[j,] / s
            }
        }

        # Refit local linear models with updated weights
        beta <- llm.1D.beta(nn.x, nn.y, nn.w, max.K, degree)

        # Compute predictions using original weights
        predictions <- wpredict.1D(beta, nn.i, orig.nn.w, nn.x, nx, max.K)

        # Check for convergence
        Dpredictions <- max(abs(predictions - predictions.old))

        # Stop if improvement has plateaued
        if (stop.C.itr.on.min && Dpredictions > old.Dpredictions) {
            predictions <- predictions.old
            beta <- beta.old
            break
        }

        # Update for next iteration
        res <- y - predictions
        predictions.old <- predictions
        beta.old <- beta
        old.Dpredictions <- Dpredictions
        itr <- itr + 1
    }

    predictions
}

#' 1D local weighed Pearson correlation model
#'
#' @param x           A numeric vector of a predictor variable.
#' @param y1          A numeric vector of the first outcome variable.
#' @param y2          A numeric vector of the second outcome variable.
#' @param grid.size   A number of grid points; was grid.size = 10*length(x), but the
#'                      results don't seem to be different from 400 which is much faster.
#'
#' @param f           The proportion of the range of x that is used within a moving window
#'                      to train the model. If NULL, the optimal value of f will be found using
#'                      minimum median absolute error optimization algorithm.
#'
#' @param bw          A bandwidth parameter.
#' @param smooth      Set to TRUE (default) to smooth the estimates of local weighted correlations.
#' @param min.K       The minimal number of x NN's that must be present in each window.
#'
#' @param n.cv.folds  The number of cross-validation folds. Used only when f = NULL. Default value: 10.
#' @param n.cv.reps   The number of repetitions of cross-validation. Used only when f = NULL. Default value: 5.
#'
#' @param nn.kernel   A kernel.
#'
#' @param n.BB        The number of Bayesian bootstrap (BB) iterations for estimates of CI's of beta's.
#' @param get.predictions.CrI    A logical parameter. If TRUE, BB with quantile determined by the value of 'level' will be used to determine the upper and lower limits of CI's.
#' @param level       A confidence level.
#' @param n.C.itr     The number of Cleveland's absolute residue based reweighting iterations for a robust estimates of mean y values.
#' @param C           A scaling of |res| parameter changing |res| to |res|/C  before applying ae.kernel to |res|'s.
#'
#' @param stop.C.itr.on.min A logical variable, if TRUE, the Cleveland's iterative reweighting stops when the maximum of the absolute values of
#'                          differences of the old and new predictions estimates are reach a local minimum.
#' @param y1.binary   Set to TRUE if y1 is a binary variable.
#' @param cv.nNN      The number of nearest neighbors in interpolate_gpredictions() used to find predictions given gpredictions in the cv_deg0_ routines.
#' @param verbose     Prints info about what is being done.
#'
#' @return A list of input parameters as well as coefficients and residues of all linear models
#'
lcor.1D <- function(x, y1, y2, grid.size = 400, f = NULL, bw = NULL, smooth = TRUE,
                   n.BB = 0, get.predictions.CrI = TRUE, level = 0.95, n.C.itr = 100,
                   C = 6, stop.C.itr.on.min = TRUE, n.cv.folds = 10, n.cv.reps = 5,
                   min.K = 5, nn.kernel = "epanechnikov", y1.binary = FALSE, cv.nNN = 3, verbose = FALSE)
{
    stopifnot(is.numeric(x))
    nx <- length(x)

    stopifnot(is.numeric(y1))
    stopifnot(length(y1)==nx)

    stopifnot(is.numeric(y2))
    stopifnot(length(y2)==nx)

    stopifnot(min.K < nx)

    if ( !is.null(f) )
    {
        stopifnot(is.numeric(f))
        stopifnot(f>0 && f <=1)
    }

    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    nn.kernel <- match.arg(nn.kernel, kernels)
    ikernel <- pmatch(nn.kernel, kernels)

    ## Sorting x and y
    o <- order(x)
    x <- x[o]
    y1 <- y1[o]
    y2 <- y2[o]

    ng <- grid.size

    ## Defining uniform grid over the range of x values
    x.range <- range(x)
    dx <- diff(x.range)
    xgrid <- seq(x.range[1], x.range[2], length = ng)
    dxgrid <- xgrid[2] - xgrid[1]

    ## Defining bandwidth, bw, as the half of the window size defined as f*dx
    if ( is.null(bw) && !is.null(f) )
    {
        window.size <- f * dx
        bw <- window.size / 2
    }

    ##
    ## Identifying x NN's from xgrid
    ##
    ## They are used for training of linear models and prediction of predictions over x
    ##
    nn <- get.knnx(x, xgrid, k = nx)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    nn.r <- get.bws(nn.d, min.K, bw)

    rw <- row.weighting(nn.d, nn.r, nn.kernel)
    nn.w <- rw$nn.w
    max.K <- rw$max.K # max.K[i] is the index of the last non-zero weight

    ## Eliminating 0 columns from nn.i, nn.d and nn.w
    cs <- colSums(nn.w)
    idx <- cs > 0
    nn.i <- nn.i[,idx]
    nn.d <- nn.d[,idx]
    nn.w <- nn.w[,idx]

    ##
    ## y1 and y2 over x NN's of xgrid
    ##
    nn.y1 <- row.eval(nn.i, y1)
    nn.y2 <- row.eval(nn.i, y2)

    lwcor.grid <- numeric(ng)
    for ( i in seq(ng) ) {
        lwcor.grid[i] <- pearson.wcor(nn.y1[i,], nn.y2[i,], nn.w[i,])
    }

    lwcor <- approx(xgrid, lwcor.grid, xout = x)$y

    smooth.lwcor.grid <- NULL
    smooth.lwcor <- NULL
    if ( smooth )
    {
        rr <- rllmf.1D(xgrid, lwcor.grid, degree = 1, bw = bw/2)
        smooth.lwcor.grid <- rr$gpredictions

        smooth.lwcor.grid <- ifelse(smooth.lwcor.grid< -1, -1, smooth.lwcor.grid)
        smooth.lwcor.grid <- ifelse(smooth.lwcor.grid>1, 1, smooth.lwcor.grid)

        smooth.lwcor <- approx(xgrid, smooth.lwcor.grid, xout = x)$y
    }

    lwcor.CI <- NULL
    smooth.lwcor.CI <- NULL
    smooth.lwcor.grid.CI <- NULL
    if ( n.BB > 0 )
    {
        if ( get.predictions.CrI )
        {
            lwcor.CI <- pearson.wcor.BB.qCrI(nn.y1, nn.y2, nn.i, nn.w, nx, n.BB, 1 - level)

            if ( smooth )
            {
                rr <- rllmf.1D(xgrid, lwcor.CI[1,], degree = 1, bw = bw/2)
                smooth.lwcor.grid.CI.lw <- rr$gpredictions

                smooth.lwcor.grid.CI.lw <- ifelse(smooth.lwcor.grid.CI.lw < -1, -1, smooth.lwcor.grid.CI.lw)
                smooth.lwcor.grid.CI.lw <- ifelse(smooth.lwcor.grid.CI.lw > 1, 1, smooth.lwcor.grid.CI.lw)

                smooth.lwcor.CI.lw <- approx(xgrid, smooth.lwcor.grid.CI.lw, xout = x)$y

                rr <- rllmf.1D(xgrid, lwcor.CI[2,], degree = 1, bw = bw/2)
                smooth.lwcor.grid.CI.up <- rr$gpredictions

                smooth.lwcor.grid.CI.up <- ifelse(smooth.lwcor.grid.CI.up < -1, -1, smooth.lwcor.grid.CI.up)
                smooth.lwcor.grid.CI.up <- ifelse(smooth.lwcor.grid.CI.up > 1, 1, smooth.lwcor.grid.CI.up)

                smooth.lwcor.CI.up <- approx(xgrid, smooth.lwcor.grid.CI.up, xout = x)$y

                smooth.lwcor.grid.CI <- rbind(smooth.lwcor.grid.CI.lw, smooth.lwcor.grid.CI.up)
                smooth.lwcor.CI <- rbind(smooth.lwcor.CI.lw, smooth.lwcor.CI.up)
            }

        } else {
            warning("Only get.predictions.CrI implemented")
        }
    }

    list(x = x,
         y1 = y1,
         y2 = y2,
         xgrid = xgrid,
         lwcor.grid = lwcor.grid,
         lwcor = lwcor,
         smooth.lwcor.grid = smooth.lwcor.grid,
         smooth.lwcor = smooth.lwcor,
         lwcor.CI = lwcor.CI,
         smooth.lwcor.CI = smooth.lwcor.CI,
         smooth.lwcor.grid.CI = smooth.lwcor.grid.CI)
}

#' 1D local linear model of two outcomes
#'
#' @param x           A numeric vector of a predictor variable.
#' @param y1          A numeric vector of the first outcome variable.
#' @param y2          A numeric vector of the second outcome variable.
#' @param grid.size   A number of grid points; was grid.size = 10*length(x), but the
#'                      results don't seem to be different from 400 which is much faster.
#' @param degree      A degree of the polynomial of x in the linear regression; 0
#'                      means weighted mean, 1 is regular linear model lm(y ~ x), and deg = d is
#'                      lm(y ~ poly(x, d)). The only allowed values are 1 and 2.
#'
#' @param f           The proportion of the range of x that is used within a moving window
#'                      to train the model. If NULL, the optimal value of f will be found using
#'                      minimum median absolute error optimization algorithm.
#'
#' @param bw          A bandwidth parameter.
#' @param min.K       The minimal number of x NN's that must be present in each window.
#'
#' @param n.cv.folds  The number of cross-validation folds. Used only when f = NULL. Default value: 10.
#' @param n.cv.reps   The number of repetitions of cross-validation. Used only when f = NULL. Default value: 5.
#'
#' @param nn.kernel   A kernel.
#'
#' @param n.BB        The number of Bayesian bootstrap (BB) iterations for estimates of CI's of beta's.
#' @param get.predictions.CrI    A logical parameter. If TRUE, BB with quantile determined by the value of 'level' will be used to determine the upper and lower limits of CI's.
#' @param level       A confidence level.
#' @param n.C.itr     The number of Cleveland's absolute residue based reweighting iterations for a robust estimates of mean y values.
#' @param C           A scaling of |res| parameter changing |res| to |res|/C  before applying ae.kernel to |res|'s.
#'
#' @param stop.C.itr.on.min A logical variable, if TRUE, the Cleveland's iterative reweighting stops when the maximum of the absolute values of
#'                          differences of the old and new predictions estimates are reach a local minimum.
#' @param y1.binary   Set to TRUE if y1 is binary.
#' @param cv.nNN      The number of nearest neighbors in interpolate_gpredictions() used to find predictions given gpredictions in the cv_deg0_ routines.
#' @param verbose     Prints info about what is being done.
#'
#' @return A list of input parameters as well as coefficients and residues of all linear models
#'
#' @export
rllm.2os.1D <- function(x, y1, y2, grid.size = 400, degree = 2, f = NULL, bw = NULL,
                       n.BB = 1000, get.predictions.CrI = TRUE, level = 0.95, n.C.itr = 100,
                       C = 6, stop.C.itr.on.min = TRUE, n.cv.folds = 10, n.cv.reps = 5,
                       min.K = 5, nn.kernel = "epanechnikov", y1.binary = FALSE, cv.nNN = 3, verbose = FALSE)
{
    stopifnot(is.numeric(x))
    nx <- length(x)

    stopifnot(is.numeric(y1))
    stopifnot(length(y1)==nx)

    stopifnot(is.numeric(y2))
    stopifnot(length(y2)==nx)

    stopifnot(min.K < nx)

    if ( !is.null(f) )
    {
        stopifnot(is.numeric(f))
        stopifnot(f>0 && f <=1)
    }

    stopifnot(as.integer(degree)==degree)
    stopifnot(degree==1 || degree==2)

    kernels <- c("epanechnikov", "triangular", "tr.exponential","normal")
    nn.kernel <- match.arg(nn.kernel, kernels)
    ikernel <- pmatch(nn.kernel, kernels)

    ## Sorting x and y
    o <- order(x)
    x <- x[o]
    y <- y[o]

    ng <- grid.size

    ## Defining uniform grid over the range of x values
    x.range <- range(x)
    dx <- diff(x.range)
    xgrid <- seq(x.range[1], x.range[2], length = ng)
    dxgrid <- xgrid[2] - xgrid[1]

    ## Defining bandwidth, bw, as the half of the window size defined as f*dx
    if ( is.null(bw) && !is.null(f) )
    {
        window.size <- f * dx
        bw <- window.size / 2
    }

    ##
    ## Identifying x NN's from xgrid
    ##
    ## They are used for training of linear models and prediction of predictions over x
    ##
    nn <- get.knnx(x, xgrid, k = nx)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    nn.r <- get.bws(nn.d, min.K, bw)

    rw <- row.weighting(nn.d, nn.r, nn.kernel)
    nn.w <- rw$nn.w
    max.K <- rw$max.K # max.K[i] is the index of the last non-zero weight

    ## Eliminating 0 columns from nn.i, nn.d and nn.w
    cs <- colSums(nn.w)
    idx <- cs > 0
    nn.i <- nn.i[,idx]
    nn.d <- nn.d[,idx]
    nn.w <- nn.w[,idx]

    ## y1 and y2 over x NN's of xgrid
    nn.y1 <- row.eval(nn.i, y1)
    nn.y2 <- row.eval(nn.i, y2)

    fp <- llm.1D.fit.and.predict(nn.i, nn.w, nn.y2, nn.y1, max.K, degree, nx)
    beta <- fp$beta
    predictions1 <- fp$predictions

    ## packaging variables of the model for prediction of predictions over xgrid
    r <- list(beta = beta,
             y.binary = y1.binary,
             xgrid = xgrid,
             ng = ng,
             min.K = min.K,
             bw = bw,
             nn.r = nn.r,
             max.K = max.K,
             degree = degree,
             nn.kernel = nn.kernel)

    class(r) <- "magelo"

    r.pred <- predict.magelo(r, newdata = xgrid)
    predictions1g <- r.pred$predictions

    predictions1.CI <- NULL
    if ( n.BB > 0 )
    {
        if ( get.predictions.CrI )
        {
            predictions1.CI <- llm.1D.fit.and.predict.global.BB.qCrI(y1.binary, nn.i, nn.w, nn.y2, nn.y1, max.K, degree, nx, n.BB, 1 - level)
        }
    }

    list(x = x,
         y1 = y1,
         y2 = y2,
         predictions1 = predictions1,
         predictions1g = predictions1g,
         predictions1.CI = predictions1.CI)
}

#' Fits 1D rllm model and generates predictions estimates for each column of a matrix Y.
#'
#' @param Y         A matrix with the number of rows that is the same as the length of x that was used to construct nn.* matrices. The main application of this routine is for the case when the columns of Y are permutations some y.
#' @param y.binary  Set to TRUE, if the values of all columns of Y are within the interval \code{[0,1]}.
#' @param nn.i      A matrix of the indices of K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.w      A matrix of weights.
#' @param nn.x      A matrix of x values over K nearest neighbors of each element of the grid.
#' @param max.K     An array of indices indicating the range where weights are not 0. Indices < max.K at i have weights > 0.
#' @param degree    A degree of the polynomial of x in the linear regression. The only allowed values are 1 and 2.
#'
#' @return  EY.grid
#'
mllm.1D.fit.and.predict <- function(Y, y.binary, nn.i, nn.w, nn.x, max.K, degree)
{
    stopifnot(is.matrix(Y))

    nrY <- nrow(Y)
    ncY <- ncol(Y)

    Tnn.i <- t(nn.i-1)
    Tnn.w <- t(nn.w)
    Tnn.x <- t(nn.x)
    nrTnn <- nrow(Tnn.i)
    ncTnn <- ncol(Tnn.i) # number of grid points

    ncoef <- degree + 1

    EY.grid <- numeric(ncTnn * ncY)

    out <- .C("C_mllm_1D_fit_and_predict",
             as.double(Y),
             as.integer(nrY),
             as.integer(ncY),
             as.integer(y.binary),

             as.integer(Tnn.i),
             as.double(Tnn.w),
             as.double(Tnn.x),
             as.integer(nrTnn),
             as.integer(ncTnn),

             as.integer(max.K-1),
             as.integer(degree),
             EY.grid = as.double(EY.grid))

    matrix(out$EY.grid, nrow = ncTnn, ncol = ncY, byrow = FALSE)
}

#' A local linear 1D model with y being a matrix. It expects bw value.
#'
#' @param x              A numeric vector of a predictor variable.
#' @param Y              A matrix such that nrow(Y) = length(x) representation many instances of the outcome variable.
#' @param bw             A bandwidth parameter.
#' @param y.binary       A logical variable. If TRUE, bw optimization is going to use a binary loss function mean(y(1-p) + (1-y)p).
#' @param with.BB        A logical parameter. Set to TRUE if Bayesian bootstraps of the column's mean function are to be returned.
# @param degree         A degree of the polynomial of x in the linear regression; 0
#' @param grid.size      A number of grid points; was grid.size = 10*length(x), but the results don't seem to be different from 400 which is much faster.
#' @param min.K          The minimal number of x NN's that must be present in each window.
#' @param nn.kernel      The name of a kernel that will be applied to NN distances of local linear regression models.
# @param n.cores        The number of cores to use.
#'
#' @return EY.grid
#'
mllm.1D <- function(x,
                   Y,
                   bw,
                   y.binary = FALSE,
                   ##degree = 0,
                   with.BB  =  FALSE,
                   grid.size = 400,
                   min.K = 15,
                   nn.kernel = "epanechnikov") # n.cores = 10)
{
    ae.kernel = "normal" ## The name (a character string) of a kernel that will be applied to absolute residues of the llm.1D().

    stopifnot(is.numeric(x))
    nx <- length(x)

    stopifnot(min.K < nx)

    stopifnot(is.matrix(Y))
    stopifnot(is.numeric(Y))
    stopifnot(is.finite(Y)) # no non-finite values allowed

    nrY <- nrow(Y)
    ncY <- ncol(Y)

    stopifnot(nx == nrY)

    ## defining a uniform grid over x range
    x.range <- range(x[is.finite(x)])
    xgrid <- seq(x.range[1], x.range[2], length = grid.size)
    ng <- grid.size

    ## kNN's
    nn <- get.knnx(x, xgrid, k = nx)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    nn.r <- get.bws(nn.d, min.K, bw)

    rw <- row.weighting(nn.d, nn.r, nn.kernel)
    nn.w <- rw$nn.w
    max.K <- rw$max.K # max.K[i] is the index of the last non-zero weight

    ## Eliminating 0 columns from nn.i, nn.d and nn.w
    cs <- colSums(nn.w)
    idx <- cs > 0
    nn.i <- nn.i[,idx]
    nn.w <- nn.w[,idx]

    ## Local linear models and corresonding wmean gpredictions estimates
    ## for now only 0 degree models are implemented
    nn.w <- nn.w / rowSums(nn.w)

    EY.grid <- matrix.weighted.means(Y, nn.i, nn.w, max.K)

    if ( with.BB ) {
        for ( j in seq(ncY) ) {
            nn.y <- row.eval(nn.i, Y[,j])
            EY.grid[,j] <- row.weighted.mean.BB(nn.y, nn.w, max.K, n.BB = 1)[,1]
        }
    }

    EY.grid
}

#' Generates Bayesian bootstrap estimates of gpredictions
#'
#' @param n.BB       The number of Bayesian bootstrap iterations
#' @param nn.i       A matrix of NN indices.
#' @param nn.x       A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.y       A matrix of y values over K nearest neighbors of each element of the grid.
#' @param y.binary   Set to TRUE if y values are within the interval \code{[0,1]}.
#' @param nn.w       A matrix of NN weights.
#' @param nx         The number of elements of x.
#' @param max.K      A vector of __indices__ indicating the range where weights are not 0.
#' @param degree     The degree of the local models.
#' @param grid.nn.i  A matrix of grid associated NN indices.
#' @param grid.nn.w  A matrix of grid associated NN weights.
#' @param grid.nn.x  A matrix of grid associated x values over NNs.
#' @param grid.max.K A vector of grid associated max.K values.
get.BB.gpredictions <- function(n.BB,
                      nn.i,
                      nn.x,
                      nn.y,
                      y.binary,
                      nn.w,
                      nx,
                      max.K,
                      degree,
                      grid.nn.i,
                      grid.nn.x,
                      grid.nn.w,
                      grid.max.K) {
    stopifnot( degree == 1 || degree == 2 )

    Tnn.i <- t(nn.i - 1)
    Tnn.x <- t(nn.x)
    Tnn.y <- t(nn.y)
    Tnn.w <- t(nn.w)
    nrTnn <- nrow(Tnn.i)
    ncTnn <- ncol(Tnn.i) # number of grid points

    Tgrid.nn.i <- t(grid.nn.i - 1)
    Tgrid.nn.x <- t(grid.nn.x)
    Tgrid.nn.w <- t(grid.nn.w)
    nrTgrid.nn <- nrow(Tgrid.nn.i)
    ncTgrid.nn <- ncol(Tgrid.nn.i)

    bb.gpredictions  <- numeric(ncTnn * n.BB)
    bb.dgpredictions <- numeric(ncTnn * n.BB)

    out <- .C("C_get_BB_Eyg",
             as.integer(n.BB),
             as.integer(Tnn.i),
             as.double(Tnn.x),
             as.double(Tnn.y),
             as.integer(y.binary),
             as.double(Tnn.w),
             as.integer(nx),
             as.integer(nrTnn),
             as.integer(ncTnn),
             as.integer(max.K-1),
             as.integer(degree),
             as.integer(Tgrid.nn.i),
             as.double(Tgrid.nn.x),
             as.double(Tgrid.nn.w),
             as.integer(nrTgrid.nn),
             as.integer(ncTgrid.nn),
             as.integer(grid.max.K-1),
             bb.dgpredictions = as.double(bb.dgpredictions),
             bb.gpredictions = as.double(bb.gpredictions))

    bb.gpredictions  <- matrix(out$bb.gpredictions, nrow = ncTnn, ncol = n.BB, byrow = FALSE)
    bb.dgpredictions <- matrix(out$bb.dgpredictions, nrow = ncTnn, ncol = n.BB, byrow = FALSE)

    list(bb.gpredictions = bb.gpredictions,
         bb.dgpredictions = bb.dgpredictions)
}

#' Estimates gpredictions Means Matrix Over a Uniform Grid for Different bandwidths.
#'
#' Given a vector, bws, of bandwidths this routine estimates a matrix of gpredictions's
#' where each column corresponds to the gpredictions estimate for a different bandwidth.
#' This is an R interface to a C routine, C_get_gpredictionss(), where the estimates of
#' gpredictions's are done. In practice, it is not intended to be used directly from R,
#' but rather to test the correctness of C_get_gpredictionss(). This routine is used in
#' different bandwidth optimization processes based on gpredictions characteristics (like
#' total absolute curvature or the number of inflection points).
#'
#' @param bws        A vector of bandwidths.
#' @param nn.i       A matrix of indices for K nearest neighbors.
#' @param nn.d       A matrix of distances to the nearest neighbors.
#' @param nn.x       A matrix of x values over K nearest neighbors of each grid point, where K is determined in the parent routine.
#' @param nn.y       A matrix of y values over K nearest neighbors of each grid point.
#' @param y.binary   Set to TRUE if y values are in the interval \code{[0,1]}.
#' @param degree     The degree of the local models.
#' @param min.K      The minimal number of elements in each set of nearest neighbors.
#' @param xgrid         A vector of grid points.
#'
#' @return A matrix of the estimates of the means, gpredictions, of y over a uniform grid where the i-th column consists of gpredictions estimates using the i-th value of the 'bws' vector.
#'
#' @examples
#' \dontrun{
#' res <- get.gpredictionss(bws, nn.i, nn.d, nn.x, nn.y, y.binary, degree, min.K, xgrid)
#' str(res)
#' }
get.gpredictionss <- function(bws,
                    nn.i,
                    nn.d,
                    nn.x,
                    nn.y,
                    y.binary,
                    degree,
                    min.K,
                    xgrid) {
    n.bws <- length(bws)

    Tnn.i <- t(nn.i - 1)
    Tnn.d <- t(nn.d)
    Tnn.x <- t(nn.x)
    Tnn.y <- t(nn.y)
    nrTnn <- nrow(Tnn.i)
    ncTnn <- ncol(Tnn.i) # number of grid points

    n.grid <- length(xgrid)

    nn <- get.knnx(xgrid, xgrid, k = n.grid)
    grid.nn.i  <- nn$nn.index
    grid.nn.d  <- nn$nn.dist

    grid.nn.x <- row.eval(nn.i, xgrid)

    Tgrid.nn.i <- t(grid.nn.i - 1)
    Tgrid.nn.d <- t(grid.nn.d)
    Tgrid.nn.x <- t(grid.nn.x)
    nrTgrid.nn <- nrow(Tgrid.nn.i)
    ncTgrid.nn <- ncol(Tgrid.nn.i)

    gpredictionss  <- numeric(ncTnn * n.bws)

    out <- .C("C_get_gpredictionss",
             as.double(bws),
             as.integer(n.bws),
             as.integer(Tnn.i),
             as.double(Tnn.d),
             as.double(Tnn.x),
             as.double(Tnn.y),
             as.integer(y.binary),
             as.integer(nrTnn),
             as.integer(ncTnn),
             as.integer(degree),
             as.integer(min.K),
             as.integer(Tgrid.nn.i),
             as.double(Tgrid.nn.d),
             as.double(Tgrid.nn.x),
             as.integer(nrTgrid.nn),
             as.integer(ncTgrid.nn),
             gpredictionss = as.double(gpredictionss))

    matrix(out$gpredictionss, nrow = ncTnn, ncol = n.bws, byrow = FALSE)
}

#' Generates Bayesian bootstrap credible intervals of gpredictions.
#'
#' @param n.BB       The number of Bayesian bootstrap iterations
#' @param nn.i       A matrix of indices of X of the nearest neighbors of each point of X.
#' @param nn.x       A matrix of x values over K nearest neighbors of each element of the grid, where K is determined in the parent routine.
#' @param nn.y       A matrix of y values over K nearest neighbors of each element of the grid.
#' @param nn.w       A matrix of NN weights.
#' @param nx         The number of elements of x.
#' @param max.K      A vector of __indices__ indicating the range where weights are not 0.
#' @param degree     The degree of the local models.
#' @param grid.nn.i  A matrix of grid associated NN indices.
#' @param grid.nn.w  A matrix of grid associated NN weights.
#' @param grid.nn.x  A matrix of grid associated x values over NNs.
#' @param grid.max.K A vector of grid associated max.K values.
#' @param y.binary   A logical variable. If TRUE, bw optimization is going to use a binary loss function mean(y(1-p) + (1-y)p).
#' @param alpha      The confidence level.
#'
#' @return A matrix of Bayesian bootstrap credible intervals of gpredictions.
#'
#' @examples
#' # TBD
get.gpredictions.CrI <- function(n.BB,
                       nn.i,
                       nn.x,
                       nn.y,
                       nn.w,
                       nx,
                       max.K,
                       degree,
                       grid.nn.i,
                       grid.nn.x,
                       grid.nn.w,
                       grid.max.K,
                       y.binary = FALSE,
                       alpha = 0.05) {
    Tnn.i <- t(nn.i - 1)
    Tnn.x <- t(nn.x)
    Tnn.y <- t(nn.y)
    Tnn.w <- t(nn.w)
    nrTnn <- nrow(Tnn.i)
    ncTnn <- ncol(Tnn.i) # number of grid points

    Tgrid.nn.i <- t(grid.nn.i - 1)
    Tgrid.nn.x <- t(grid.nn.x)
    Tgrid.nn.w <- t(grid.nn.w)
    nrTgrid.nn <- nrow(Tgrid.nn.i)
    ncTgrid.nn <- ncol(Tgrid.nn.i)

    gpredictions.CrI <- numeric(ncTnn * 2)

    out <- .C("C_get_Eyg_CrI",
             as.integer(y.binary),
             as.integer(n.BB),
             as.integer(Tnn.i),
             as.double(Tnn.x),
             as.double(Tnn.y),
             as.double(Tnn.w),
             as.integer(nx),
             as.integer(nrTnn),
             as.integer(ncTnn),
             as.integer(max.K-1),
             as.integer(degree),
             as.integer(Tgrid.nn.i),
             as.double(Tgrid.nn.x),
             as.double(Tgrid.nn.w),
             as.integer(nrTgrid.nn),
             as.integer(ncTgrid.nn),
             as.integer(grid.max.K-1),
             as.double(alpha),
             gpredictions.CrI = as.double(gpredictions.CrI))

    matrix(out$gpredictions.CrI, nrow = 2, ncol = ncTnn, byrow = FALSE)
}


#' Row-wise evaluates x at nn.i.
#'
#' @param nn.i      An array of indices of K nearest neighbors of the i-th element of x.
#' @param x         An array of nx elements.
#'
#' @return A matrix obtained from x by evaluating it at the indices of nn.i.
row.eval <- function(nn.i, x) {
    nx <- length(x)

    Tnn.i <- t(nn.i-1)
    nr.Tnn.i <- nrow(Tnn.i)
    nc.Tnn.i <- ncol(Tnn.i)

    nn.x <- numeric(nr.Tnn.i * nc.Tnn.i)

    out <- .C("C_columnwise_eval",
             as.integer(Tnn.i),
             as.integer(nr.Tnn.i),
             as.integer(nc.Tnn.i),
             as.double(x),
             nn.x=as.double(nn.x))

    nn.x <- matrix(out$nn.x, nrow=nrow(nn.i), byrow=TRUE)

    return(nn.x)
}

#' 1D MAGELO (Model Averaged Grid-based Epsilon LOwess) with fixed boundary condition.
#' Given a vector 'y', this routine finds 'predictions' such that the first element of
#' 'predictions' equals the first element of 'y', and the last element of 'predictions'
#' equals the last element of 'y'.
#'
#' @param x           A numeric vector of a predictor variable.
#' @param y           A numeric vector of an outcome variable.
#' @param grid.size   The number of grid points; default 400 which provides good balance
#'                    between accuracy and computation speed.
#' @param degree      A degree of the polynomial of x in the linear regression; 0
#'                    means weighted mean, 1 is regular linear model lm(y ~ x), and deg = d is
#'                    lm(y ~ poly(x, d)). The only allowed values are 0, 1 and 2.
#' @param min.K       The minimal number of x points that must be present in each disk neighborhood.
#' @param f           The proportion of the range of x that is used within a moving window
#'                    to train the model. If NULL, the optimal value of f will be found using
#'                    cross-validation optimization algorithm.
#' @param bw          A bandwidth parameter (radius of disk neighborhood).
#' @param min.bw.f    The min.bw factor, such that, min.bw = min.bw.f * dx, where dx <- diff(x.range).
#'                    The default value is 0.025.
#' @param method      A method of estimating the optimal value of the bandwidth.
#'                    Possible choices are "LOOCV" (for small datasets) and "CV".
#' @param n.bws       The number of bandwidths in the optimization process. Default: 100.
#' @param n.cv.folds  The number of cross-validation folds. Used only when f = NULL. Default value: 10.
#' @param n.cv.reps   The number of repetitions of cross-validation. Used only when f = NULL. Default value: 20.
#'
#' @param nn.kernel   A kernel for weighting neighbors. Options: "epanechnikov", "triangular",
#'                    "tr.exponential", "normal".
#'
#' @param n.BB        The number of Bayesian bootstrap (BB) iterations for estimates of CI's.
#' @param get.predictions.CrI  A logical parameter. If TRUE, BB with quantile determined by the
#'                            value of 'level' will be used to determine the upper and lower limits of CI's.
#' @param get.gpredictions.CrI Set to TRUE if gpredictions.CI's need to be estimated.
#' @param level       A confidence level.
#' @param n.C.itr     The number of Cleveland's absolute residue based reweighting iterations
#'                    for robust estimates of mean y values.
#' @param C           A scaling of |res| parameter changing |res| to |res|/C before applying
#'                    ae.kernel to |res|'s.
#'
#' @param stop.C.itr.on.min A logical variable, if TRUE, the Cleveland's iterative reweighting
#'                          stops when the maximum of the absolute values of differences of
#'                          the old and new predictions estimates reach a local minimum.
#' @param y.binary    A logical variable. If TRUE, bw optimization is going to use a binary
#'                    loss function mean(y(1-p) + (1-y)p).
#' @param cv.nNN      The number of nearest neighbors in interpolation used to find predictions
#'                    given gpredictions in the cv routines.
#' @param get.BB.predictions   Set to TRUE if a matrix of Bayesian bootstrap estimates of
#'                            predictions is to be returned.
#' @param get.BB.gpredictions  Set to TRUE if a matrix of Bayesian bootstrap estimates of
#'                            gpredictions is to be returned.
#' @param fb.C        The number of bw's away from the boundary points of the domain of x
#'                    that the adjusting of gpredictions will take place, so that gpredictions
#'                    satisfies the boundary condition.
#' @param n.perms     Number of y permutations for p-value calculation. Default: 0.
#' @param n.cores     Number of CPU cores for parallel processing. Default: 1.
#' @param use.binloss Logical; if TRUE, use binary loss function for optimization.
#' @param verbose     Prints info about what is being done.
#'
#' @return A list of class "magelo" containing input parameters as well as fitted values,
#'         credible intervals, bootstrap estimates, and model coefficients
#'
#' @export
fb.magelo <- function(x,
                      y,
                      grid.size = 400,
                      degree = 1,
                      min.K = 5,
                      f = NULL,
                      bw = NULL,
                      min.bw.f = 0.025,
                      method = ifelse(length(x) < 1000, "LOOCV", "CV"),
                      n.bws = 100,
                      n.BB = 1000,
                      get.predictions.CrI = TRUE,
                      get.gpredictions.CrI = TRUE,
                      get.BB.predictions = FALSE,
                      get.BB.gpredictions = FALSE,
                      level = 0.95,
                      n.C.itr = 0,
                      C = 6,
                      stop.C.itr.on.min = TRUE,
                      n.cv.folds = 10,
                      n.cv.reps = 20,
                      nn.kernel = "epanechnikov",
                      y.binary = FALSE,
                      cv.nNN = 3,
                      n.perms = 0,
                      n.cores = 1,
                      use.binloss = FALSE,
                      verbose = FALSE,
                      fb.C = 1)
{
    ## Call magelo with the provided parameters
    r <- magelo(x = x,
                y = y,
                y.true = NULL,
                grid.size = grid.size,
                degree = degree,
                min.K = min.K,
                f = f,
                bw = bw,
                min.bw.f = min.bw.f,
                method = method,
                n.bws = n.bws,
                n.BB = n.BB,
                get.predictions.CrI = get.predictions.CrI,
                get.gpredictions.CrI = get.gpredictions.CrI,
                get.BB.predictions = get.BB.predictions,
                get.BB.gpredictions = get.BB.gpredictions,
                level = level,
                n.C.itr = n.C.itr,
                C = C,
                stop.C.itr.on.min = stop.C.itr.on.min,
                n.cv.folds = n.cv.folds,
                n.cv.reps = n.cv.reps,
                nn.kernel = nn.kernel,
                y.binary = y.binary,
                cv.nNN = cv.nNN,
                n.perms = n.perms,
                n.cores = n.cores,
                use.binloss = use.binloss,
                verbose = verbose)

    ## lambdaL is defined over [x0, x0 + C*bw]
    lambdaL <- function(x, x0, bw, C = 1) {
        -1/(C*bw) * (x - x0) + 1
    }

    ## lambdaR is defined over [x1 - C*bw, x1]
    lambdaR <- function(x, x1, bw, C = 1) {
        1/(C*bw) * (x - x1 + C*bw)
    }

    ##
    ## left end modification
    ##
    idx <- r$xgrid < r$xgrid[1] + fb.C*r$opt.bw
    xgridL <- r$xgrid[idx]

    fbgpredictions <- r$gpredictions
    for ( i in seq(xgridL) )
    {
        t <- xgridL[i]
        lambda <- lambdaL(t, xgridL[1], r$opt.bw, fb.C)
        lambda <- lambda^2
        fbgpredictions[i] <- lambda * y[1] + (1 - lambda) * r$gpredictions[i]
    }

    ##
    ## right end modification
    ##
    n <- length(r$xgrid)
    idx <- r$xgrid > r$xgrid[n] - fb.C*r$opt.bw
    xgridR <- r$xgrid[idx]

    i0 <- which(idx)[1] - 1
    t1 <- xgridR[length(xgridR)]
    for ( i in seq(xgridR) )
    {
        t <- xgridR[i]
        lambda <- lambdaR(t, t1, r$opt.bw, fb.C)
        lambda <- lambda^2
        fbgpredictions[i0 + i] <- lambda * y[length(y)] + (1 - lambda) * r$gpredictions[i0 + i]
    }

    predictions <- approx(r$xgrid, fbgpredictions, xout = x)$y

    output <- list(gpredictions = fbgpredictions,
                  predictions = predictions,
                  n.BB = n.BB,
                  BB.predictions = r$BB.predictions,
                  BB.gpredictions = r$BB.gpredictions,
                  BB.dgpredictions = r$BB.dgpredictions,
                  gpredictions.CrI = r$gpredictions.CrI,
                  gpredictions.CrI.smooth = r$gpredictions.CrI.smooth,
                  predictions.CrI = r$predictions.CrI,
                  min.error = r$min.error,
                  beta = r$beta,
                  min.bw = r$min.bw,
                  max.bw = r$max.bw,
                  opt.bw = r$opt.bw,
                  opt.bw.i = r$opt.bw.i,
                  log.bws = r$log.bws,
                  errors = r$errors,
                  x = r$params$x,
                  y = r$params$y,
                  grid.size = r$params$grid.size,
                  degree = r$params$degree,
                  nn.kernel = r$params$nn.kernel,
                  min.K = r$params$min.K,
                  max.K = r$params$max.K,
                  y.binary = y.binary,
                  xgrid = r$xgrid,
                  params = r$params)

    class(output) <- "magelo"

    return(output)
}
