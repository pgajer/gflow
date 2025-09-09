#' Fit a Piecewise Linear Model
#'
#' @description
#' Fits a piecewise linear model to data using the segmented package. Can fit models
#' with single or multiple breakpoints. If the segmented regression fails, falls back
#' to simple linear regression.
#'
#' @param x Numeric vector containing the predictor variable values
#' @param y Numeric vector containing the response variable values
#' @param n_breakpoints Integer specifying the number of breakpoints to estimate (default = 1)
#' @param breakpoint_bounds List of vectors specifying the bounds for each breakpoint (optional)
#'
#' @return An object of class "pwlm" containing:
#'   \item{model}{Either a segmented model object or a linear model object if segmented regression failed}
#'   \item{breakpoints}{Numeric vector of estimated breakpoints. NA if using fallback linear model}
#'   \item{x}{The original x values}
#'   \item{y}{The original y values}
#'   \item{type}{Character string indicating "segmented" or "linear"}
#'   \item{n_breakpoints}{Number of breakpoints specified}
#'
#' @examples
#' # Generate sample data
#' x <- 1:20
#' y <- x + 2*x^2 + rnorm(20, 0, 10)
#'
#' # Single breakpoint
#' fit1 <- fit.pwlm(x, y)
#'
#' # Two breakpoints
#' fit2 <- fit.pwlm(x, y, n_breakpoints = 2)
#' @export
fit.pwlm <- function(x,
                     y,
                     n_breakpoints = 1,
                     breakpoint_bounds = NULL) {
    ## Input validation
    if (!is.numeric(n_breakpoints) || n_breakpoints < 1) {
        stop("n_breakpoints must be a positive integer")
    }

    if (!is.null(breakpoint_bounds)) {
        if (!is.list(breakpoint_bounds) || length(breakpoint_bounds) != n_breakpoints) {
            stop("breakpoint_bounds must be a list with length equal to n_breakpoints")
        }
        # Validate each bound
        for (bound in breakpoint_bounds) {
            if (length(bound) != 2 || !all(bound >= min(x)) || !all(bound <= max(x))) {
                stop("Each breakpoint bound must be a vector of length 2 within the range of x")
            }
        }
    }

    if (requireNamespace("segmented", quietly = TRUE)) {
        init.lm <- lm(y ~ x)

        # Initialize breakpoints
        if (is.null(breakpoint_bounds)) {
            # Equally spaced initial guesses
            breaks <- seq(min(x), max(x), length.out = n_breakpoints + 2)[2:(n_breakpoints + 1)]
            psi <- list(x = breaks)
        } else {
            # Use midpoints of specified bounds as initial guesses
            breaks <- sapply(breakpoint_bounds, function(bound) mean(bound))
            psi <- list(x = breaks)
        }

        # Try fitting the segmented model
        pwlm <- try(segmented::segmented(init.lm, seg.Z = ~x, psi = psi,
                                       control = list(K = n_breakpoints)),
                    silent = TRUE)

        if (!inherits(pwlm, "try-error")) {
            result <- list(
                model = pwlm,
                breakpoints = pwlm$psi[, 2],
                x = x,
                y = y,
                type = "segmented",
                n_breakpoints = n_breakpoints
            )
            class(result) <- "pwlm"
            return(result)
        }
    }

    # Fallback if segmented regression fails
    result <- list(
        model = lm(y ~ x),
        breakpoints = rep(NA, n_breakpoints),
        x = x,
        y = y,
        type = "linear",
        n_breakpoints = n_breakpoints
    )
    class(result) <- "pwlm"
    return(result)
}

#' Plot a Piecewise Linear Model
#'
#' @description
#' Creates a visualization of a piecewise linear model fit, showing the data points,
#' fitted lines, and breakpoints. This function handles both segmented models with
#' breakpoints and simple linear models.
#'
#' @param x An object of class "pwlm", typically the output of fit.pwlm().
#' @param main Character string for plot title. Default is "Piecewise Linear Regression".
#' @param xlab Character string for x-axis label. Default is "X".
#' @param ylab Character string for y-axis label. Default is "Y".
#' @param point_color Color for data points. Default is "black".
#' @param line_color Color for fitted lines. Default is "blue".
#' @param breakpoint_color Color for breakpoint vertical lines. Default is "red".
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' x <- 1:20
#' y <- c(1:10, 20:11) + rnorm(20, 0, 2)
#'
#' # Fit piecewise linear model with one breakpoint
#' pwlm_fit <- fit.pwlm(x, y)
#'
#' # Plot the model
#' \dontrun{
#'   plot(pwlm_fit)
#'
#'   # Customize plot appearance
#'   plot(pwlm_fit, main = "My Custom PWLM Plot",
#'        point_color = "blue", line_color = "red")
#' }
#'
#' @export
#' @importFrom graphics lines abline legend
#' @importFrom stats predict
plot.pwlm <- function(x,
                      main = "Piecewise Linear Regression",
                      xlab = "X",
                      ylab = "Y",
                      point_color = "black",
                      line_color = "blue",
                      breakpoint_color = "red",
                     ...) {
    ## Create the base plot with data points
    plot(x$x, x$y,
         main = main,
         xlab = xlab,
         ylab = ylab,
         pch = 16,
         col = point_color,
         ...)

    # If we have a segmented model
    if (x$type == "segmented") {
        # Get predicted values for smooth line
        new_x <- seq(min(x$x), max(x$x), length.out = 200)
        pred_y <- predict(x$model, newdata = data.frame(x = new_x))

        # Add the fitted lines
        lines(new_x, pred_y, col = line_color, lwd = 2)

        # Add vertical lines at breakpoints
        for (bp in x$breakpoints) {
            abline(v = bp, col = breakpoint_color, lty = 2, lwd = 1)
        }

        # Add legend
        legend("topleft",
               legend = c("Data", "Fitted Line", "Breakpoints"),
               col = c(point_color, line_color, breakpoint_color),
               pch = c(16, NA, NA),
               lty = c(NA, 1, 2),
               lwd = c(NA, 2, 1))
    } else {
        # If no breakpoint, just plot simple linear regression
        abline(x$model, col = line_color, lwd = 2)

        # Add legend for simple case
        legend("topleft",
               legend = c("Data", "Fitted Line"),
               col = c(point_color, line_color),
               pch = c(16, NA),
               lty = c(NA, 1),
               lwd = c(NA, 2))
    }
}

#' Fit a Piecewise Linear Model with Optimal Number of Breakpoints
#'
#' @description
#' Fits a piecewise linear model to data, automatically determining the optimal
#' number of breakpoints using model selection criteria.
#'
#' @param x Numeric vector containing the predictor variable values
#' @param y Numeric vector containing the response variable values
#' @param max_breakpoints Maximum number of breakpoints to consider (default = 5)
#' @param method Character string specifying the selection method ("aic", "bic", or "davies")
#' @param alpha Significance level for Davies' test (default = 0.05)
#' @param plot_selection Logical indicating whether to plot selection criteria (default = FALSE)
#'
#' @return An object of class "pwlm" with additional component 'selection_results'
#'
#' @export
fit.pwlm.optimal <- function(x,
                             y,
                             max_breakpoints = 5,
                             method = c("aic", "bic", "davies"),
                             alpha = 0.05,
                             plot_selection = FALSE) {
    method <- match.arg(method)

    # Find optimal number of breakpoints
    optimal <- find_optimal_breakpoints(x, y, max_breakpoints, method, alpha)

    # Fit model with optimal number of breakpoints
    model <- fit.pwlm(x, y, n_breakpoints = optimal$optimal_breakpoints)

    # Add selection results to model object
    model$selection_results <- optimal$results

    # Plot selection criteria if requested
    if (plot_selection) {
        par(mfrow = c(2, 1))

        # Plot AIC/BIC
        plot(optimal$results$n_breakpoints, optimal$results$aic,
             type = "b", col = "blue", ylab = "Criterion Value",
             xlab = "Number of Breakpoints", main = "Model Selection Criteria")
        lines(optimal$results$n_breakpoints, optimal$results$bic,
              type = "b", col = "red")
        legend("topright", c("AIC", "BIC"), col = c("blue", "red"), lty = 1)

        # Plot Davies' test p-values
        plot(optimal$results$n_breakpoints[-1], optimal$results$davies_pvalue[-1],
             type = "b", col = "purple", ylab = "P-value",
             xlab = "Number of Breakpoints", main = "Davies' Test P-values")
        abline(h = alpha, lty = 2, col = "gray")

        par(mfrow = c(1, 1))
    }

    return(model)
}


#' Print Method for PWLM Objects
#'
#' @param x An object of class "pwlm"
#' @param ... Additional arguments passed to print
#'
#' @export
print.pwlm <- function(x, ...) {
    cat("Piecewise Linear Model\n")
    cat("---------------------\n")
    cat("Model type:", x$type, "\n")
    if (!is.na(x$breakpoint)) {
        cat("Breakpoint at x =", round(x$breakpoint, 4), "\n")
    }
    cat("\nModel Summary:\n")
    print(summary(x$model))
}

#' Summary Method for PWLM Objects
#'
#' @param object An object of class "pwlm"
#' @param ... Additional arguments passed to summary
#'
#' @export
summary.pwlm <- function(object, ...) {
    result <- list(
        type = object$type,
        breakpoint = object$breakpoint,
        model_summary = summary(object$model)
    )
    class(result) <- "summary.pwlm"
    return(result)
}

#' Print Method for PWLM Summary Objects
#'
#' @param x An object of class "summary.pwlm"
#' @param ... Additional arguments passed to print
#'
#' @export
print.summary.pwlm <- function(x, ...) {
    cat("Piecewise Linear Model Summary\n")
    cat("-----------------------------\n")
    cat("Model type:", x$type, "\n")
    if (!is.na(x$breakpoint)) {
        cat("Breakpoint at x =", round(x$breakpoint, 4), "\n")
    }
    cat("\nDetailed Model Summary:\n")
    print(x$model_summary)
}

#' Find Optimal Number of Breakpoints
#'
#' @param x Numeric vector containing the predictor variable values
#' @param y Numeric vector containing the response variable values
#' @param max_breakpoints Maximum number of breakpoints to consider (default = 5)
#' @param method Character string specifying the selection method ("aic", "bic", or "davies")
#' @param alpha Significance level for Davies' test (default = 0.05)
#'
#' @return A list containing the optimal number of breakpoints and model comparison results
#'
#' @importFrom stats AIC BIC
find_optimal_breakpoints <- function(x,
                                     y,
                                     max_breakpoints = 5,
                                     method = c("aic", "bic", "davies"),
                                     alpha = 0.05) {
    method <- match.arg(method)

    # Initialize results storage
    results <- data.frame(
        n_breakpoints = 0:max_breakpoints,
        aic = NA_real_,
        bic = NA_real_,
        davies_pvalue = NA_real_
    )

    # Fit linear model (0 breakpoints)
    lm_fit <- lm(y ~ x)
    results$aic[1] <- AIC(lm_fit)
    results$bic[1] <- BIC(lm_fit)

    # Fit models with increasing number of breakpoints
    for (i in 1:max_breakpoints) {
        model <- try(fit.pwlm(x, y, n_breakpoints = i), silent = TRUE)

        if (!inherits(model, "try-error") && model$type == "segmented") {
            results$aic[i + 1] <- AIC(model$model)
            results$bic[i + 1] <- BIC(model$model)

            # Davies' test comparing to model with one fewer breakpoint
            if (i == 1) {
                davies_test <- try(segmented::davies.test(lm_fit, ~x, k = 1), silent = TRUE)
            } else {
                prev_model <- fit.pwlm(x, y, n_breakpoints = i - 1)
                davies_test <- try(segmented::davies.test(prev_model$model, ~x, k = 1), silent = TRUE)
            }

            if (!inherits(davies_test, "try-error")) {
                results$davies_pvalue[i + 1] <- davies_test$p.value
            }
        }
    }

    # Determine optimal number based on selected method
    optimal <- switch(method,
        "aic" = {
            which.min(results$aic) - 1
        },
        "bic" = {
            which.min(results$bic) - 1
        },
        "davies" = {
            # Find last significant breakpoint
            max(0, max(which(results$davies_pvalue <= alpha)) - 1)
        }
    )

    return(list(
        optimal_breakpoints = optimal,
        results = results
    ))
}
