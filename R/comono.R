## Co-monotonicity Coefficient Computation
##
## Functions for computing co-monotonicity coefficients between two functions
## defined on graph vertices. These measures quantify the extent to which two
## functions vary together across graph edges.

#' Compute Co-monotonicity Coefficient
#'
#' @description
#' Computes the vertex-wise co-monotonicity coefficient between two functions
#' defined on the vertices of a weighted graph. The co-monotonicity coefficient
#' measures the extent to which two functions vary together across graph edges,
#' providing a graph-based analogue of correlation that respects local geometric
#' structure.
#'
#' @details
#' Consider two real-valued functions \eqn{y} and \eqn{z} defined on the vertices
#' of a weighted graph \eqn{G = (V, E)}. For an edge \eqn{e = [v,u]} connecting
#' vertices \eqn{v} and \eqn{u}, we define the edge difference operator
#' \eqn{\Delta_e f = f(u) - f(v)}.
#'
#' The co-monotonicity coefficient at vertex \eqn{v} is defined as:
#' \deqn{\text{comono}(y,z;w)(v) = \frac{\sum_{u \in N(v)} w_e \Delta_e y \Delta_e z}{\sum_{u \in N(v)} w_e |\Delta_e y \Delta_e z|}}
#'
#' where \eqn{N(v)} denotes the neighborhood of vertex \eqn{v} and \eqn{w_e \geq 0}
#' are edge weights. The coefficient takes values in \eqn{[-1, 1]}, where:
#' \itemize{
#'   \item +1 indicates perfect positive co-monotonicity (both functions always
#'         increase or decrease together)
#'   \item -1 indicates perfect negative co-monotonicity (functions always change
#'         in opposite directions)
#'   \item 0 indicates no systematic relationship
#' }
#'
#' @section Weighting Schemes:
#'
#' Three weighting schemes are supported:
#'
#' \strong{Unit weights} (\code{type = "unit"}): Uses \eqn{w_e = 1} for all edges.
#' This treats all edges equally regardless of their length. Appropriate when edge
#' lengths are roughly comparable or when counting directional agreements without
#' geometric normalization.
#'
#' \strong{Derivative weights} (\code{type = "derivative"}): Uses
#' \eqn{w_e = 1/(\Delta_e)^2} where \eqn{\Delta_e} is the edge length. This
#' normalizes by edge length squared, making the measure analogous to comparing
#' derivatives rather than absolute changes. Natural when functions represent
#' continuous quantities sampled at irregular spatial positions.
#'
#' \strong{Sign-based} (\code{type = "sign"}): Uses only the sign of
#' \eqn{\Delta_e y \cdot \Delta_e z}, computing:
#' \deqn{\text{comono}_{\pm}(y,z)(v) = \frac{\sum_{u \in N(v)} \text{sign}(\Delta_e y \Delta_e z)}{|N(v)|}}
#' This robust measure is insensitive to outliers and magnitude of change.
#'
#' @section Interpretation:
#'
#' The vertex-wise coefficients provide a spatially-resolved view of co-monotonicity,
#' revealing where in the graph the two functions are most strongly associated.
#' The summary statistics characterize the global pattern of association.
#'
#' For example, when comparing observed data \eqn{y} with model predictions
#' \eqn{\hat{y}}, vertices with low co-monotonicity indicate regions where the
#' model may poorly represent local directional behavior.
#'
#' @param adj.list A list of integer vectors representing the graph's adjacency
#'   structure. Element \code{i} contains the indices of vertices adjacent to
#'   vertex \code{i}. Indices should be 1-based (R convention).
#' @param weight.list A list of numeric vectors containing edge weights.
#'   Element \code{i} contains weights corresponding to the edges in
#'   \code{adj.list[[i]]}. Must have the same structure as \code{adj.list}.
#' @param y Numeric vector of function values at each vertex. Length must equal
#'   the number of vertices (i.e., \code{length(adj.list)}).
#' @param z Numeric vector of function values at each vertex. Length must equal
#'   the number of vertices.
#' @param type Character string specifying the weighting scheme. One of:
#'   \code{"unit"} (default), \code{"derivative"}, or \code{"sign"}.
#'
#' @return A list with class \code{"comono_result"} containing:
#'   \item{vertex.coefficients}{Numeric vector of co-monotonicity coefficients at
#'     each vertex}
#'   \item{mean.coefficient}{Mean co-monotonicity across all vertices}
#'   \item{median.coefficient}{Median co-monotonicity across all vertices}
#'   \item{n.positive}{Number of vertices with positive co-monotonicity}
#'   \item{n.negative}{Number of vertices with negative co-monotonicity}
#'   \item{n.zero}{Number of vertices with zero co-monotonicity}
#'   \item{type}{Weighting scheme used}
#'   \item{n.vertices}{Total number of vertices in the graph}
#'
#' @examples
#' \dontrun{
#' ## Create example graph
#' library(gflow)
#' data("sPTB_microbiome")
#'
#' ## Build k-NN graph
#' graph <- build.knn.graph(X, k = 10)
#' adj.list <- graph$adj.list
#' weight.list <- graph$weight.list
#'
#' ## Fit regression models with different smoothing parameters
#' y.hat.1 <- gflow.regression(adj.list, weight.list, y, lambda = 0.1)
#' y.hat.2 <- gflow.regression(adj.list, weight.list, y, lambda = 0.5)
#'
#' ## Compare co-monotonicity with observations
#' comono.1 <- comono(adj.list, weight.list, y, y.hat.1, type = "unit")
#' comono.2 <- comono(adj.list, weight.list, y, y.hat.2, type = "unit")
#'
#' ## Model 1 has higher co-monotonicity
#' comono.1$mean.coefficient  # 0.87
#' comono.2$mean.coefficient  # 0.72
#'
#' ## Identify vertices with poor directional fit
#' poor.fit <- which(comono.1$vertex.coefficients < 0.5)
#'
#' ## Visualize spatial distribution of co-monotonicity
#' plot(graph, vertex.color = comono.1$vertex.coefficients,
#'      vertex.size = 3, main = "Co-monotonicity: y vs ŷ")
#' }
#'
#' @references
#' Gajer, P. (2025). Morse-Smale regression: A geometric approach to
#' statistical analysis on graphs. \emph{In preparation}.
#'
#' @seealso
#' \code{\link{comono.global}} for scalar summary,
#' \code{\link{plot.comono_result}} for visualization
#'
#' @export
comono <- function(adj.list, weight.list, y, z, type = c("unit", "derivative", "sign")) {
    ## Validate type argument
    type <- match.arg(type)
    
    ## Input validation
    if (!is.list(adj.list) || !is.list(weight.list)) {
        stop("adj.list and weight.list must be lists")
    }
    
    if (length(adj.list) != length(weight.list)) {
        stop("adj.list and weight.list must have the same length")
    }
    
    n.vertices <- length(adj.list)
    
    if (length(y) != n.vertices) {
        stop("Length of y must match the number of vertices")
    }
    
    if (length(z) != n.vertices) {
        stop("Length of z must match the number of vertices")
    }
    
    ## Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))
    
    ## Call C++ function
    result <- .Call(S_comono,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    as.numeric(z),
                    type,
                    PACKAGE = "gflow")
    
    ## Add metadata
    result$type <- type
    result$n.vertices <- n.vertices
    
    ## Set class for S3 method dispatch
    class(result) <- c("comono_result", "list")
    
    return(result)
}


#' Compute Global Co-monotonicity Coefficient
#'
#' @description
#' Computes a single scalar co-monotonicity coefficient summarizing the
#' association between two functions across the entire graph. This is defined
#' as the mean of the vertex-wise coefficients.
#'
#' @details
#' This function provides a convenient summary when a single global measure of
#' co-monotonicity is desired. For spatially-resolved analysis revealing local
#' patterns of association, use \code{\link{comono}} which returns vertex-wise
#' coefficients.
#'
#' @inheritParams comono
#'
#' @return A scalar numeric value in \eqn{[-1, 1]} representing the global
#'   co-monotonicity coefficient
#'
#' @examples
#' \dontrun{
#' ## Global co-monotonicity between observations and predictions
#' global.comono <- comono.global(adj.list, weight.list, y, y.hat, type = "unit")
#'
#' ## Compare different models
#' comono.model.1 <- comono.global(adj.list, weight.list, y, y.hat.1)
#' comono.model.2 <- comono.global(adj.list, weight.list, y, y.hat.2)
#'
#' if (comono.model.1 > comono.model.2) {
#'     cat("Model 1 has better directional agreement\n")
#' }
#' }
#'
#' @seealso \code{\link{comono}} for vertex-wise coefficients
#'
#' @export
comono.global <- function(adj.list, weight.list, y, z, type = c("unit", "derivative", "sign")) {
    result <- comono(adj.list, weight.list, y, z, type)
    return(result$mean.coefficient)
}


#' Print Method for Co-monotonicity Results
#'
#' @param x An object of class \code{"comono_result"}
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @method print comono_result
#' @export
print.comono_result <- function(x, digits = 4, ...) {
    cat("Co-monotonicity Coefficient Results\n")
    cat("====================================\n\n")
    
    cat("Weighting scheme:", x$type, "\n")
    cat("Number of vertices:", x$n.vertices, "\n\n")
    
    cat("Summary Statistics:\n")
    cat(sprintf("  Mean:   %.*f\n", digits, x$mean.coefficient))
    cat(sprintf("  Median: %.*f\n", digits, x$median.coefficient))
    cat("\n")
    
    cat("Directional Agreement:\n")
    total <- x$n.positive + x$n.negative + x$n.zero
    cat(sprintf("  Positive: %d (%.1f%%)\n",
                x$n.positive, 100 * x$n.positive / total))
    cat(sprintf("  Negative: %d (%.1f%%)\n",
                x$n.negative, 100 * x$n.negative / total))
    cat(sprintf("  Zero:     %d (%.1f%%)\n",
                x$n.zero, 100 * x$n.zero / total))
    
    invisible(x)
}


#' Summary Method for Co-monotonicity Results
#'
#' @param object An object of class \code{"comono_result"}
#' @param ... Additional arguments (ignored)
#'
#' @return A summary data frame
#'
#' @method summary comono_result
#' @export
summary.comono_result <- function(object, ...) {
    coeffs <- object$vertex.coefficients
    
    summary.df <- data.frame(
        Statistic = c("Minimum", "1st Quartile", "Median", "Mean",
                      "3rd Quartile", "Maximum", "Std. Dev."),
        Value = c(
            min(coeffs),
            quantile(coeffs, 0.25),
            object$median.coefficient,
            object$mean.coefficient,
            quantile(coeffs, 0.75),
            max(coeffs),
            sd(coeffs)
        )
    )
    
    cat("\nCo-monotonicity Coefficient Summary\n")
    cat("===================================\n")
    cat("Weighting scheme:", object$type, "\n")
    cat("Number of vertices:", object$n.vertices, "\n\n")
    
    print(summary.df, row.names = FALSE, digits = 4)
    
    cat("\nDirectional Composition:\n")
    total <- object$n.positive + object$n.negative + object$n.zero
    cat(sprintf("  Positive: %d vertices (%.1f%%)\n",
                object$n.positive, 100 * object$n.positive / total))
    cat(sprintf("  Negative: %d vertices (%.1f%%)\n",
                object$n.negative, 100 * object$n.negative / total))
    cat(sprintf("  Zero:     %d vertices (%.1f%%)\n",
                object$n.zero, 100 * object$n.zero / total))
    
    invisible(summary.df)
}


#' Compute Co-monotonicity Coefficients Between Vector and Matrix Columns
#'
#' @description
#' Efficiently computes vertex-wise co-monotonicity coefficients between a vector
#' y and each column of a matrix Z. This function is optimized for computing
#' multiple co-monotonicity measures simultaneously by reusing graph structure
#' information across all columns.
#'
#' @details
#' This function extends the scalar \code{\link{comono}} to the multivariate case,
#' computing co-monotonicity between y and each column z_j of Z. The implementation
#' is optimized to avoid redundant computation: graph-dependent quantities (edge
#' differences for y, edge weights) are computed once and reused for all columns.
#'
#' The vertex-wise co-monotonicity at vertex v for column j is:
#' \deqn{\text{comono}(y,z_j)(v) = \frac{\sum_{u \in N(v)} w_e \Delta_e y \Delta_e z_j}{\sum_{u \in N(v)} w_e |\Delta_e y \Delta_e z_j|}}
#'
#' This formulation is identical to the scalar case but applied independently to
#' each column of Z.
#'
#' @section Use Cases:
#'
#' \strong{Feature Selection}: When Z contains multiple candidate predictors,
#' columns with high mean absolute co-monotonicity with y indicate features that
#' strongly capture y's directional behavior. This provides a graph-based alternative
#' to correlation-based feature selection that respects local neighborhood structure.
#'
#' \strong{Model Comparison}: When Z contains predictions from multiple models,
#' comparing mean co-monotonicity values reveals which model best captures the
#' directional patterns in y. Higher co-monotonicity indicates better alignment
#' with local trends.
#'
#' \strong{Multivariate Regression Diagnostics}: When Z contains multiple covariates,
#' examining co-monotonicity between residuals and each covariate helps identify
#' nonlinear relationships or model misspecification patterns.
#'
#' \strong{Time Series Analysis}: When columns of Z represent temporal snapshots,
#' co-monotonicity reveals which time points exhibit strongest concordance with y.
#'
#' @section Performance:
#'
#' For q columns and a graph with n vertices and m edges, the function performs:
#' \itemize{
#'   \item One-time setup: O(m) to pre-compute y-edge differences and weights
#'   \item Per-column processing: O(m) to compute z-edge differences and products
#'   \item Total: O(m + q·m) compared to O(q·m) for q separate calls to \code{comono}
#' }
#'
#' The speedup is modest (approximately 2x for large q) but memory access patterns
#' are more cache-friendly, and the single function call simplifies user code.
#'
#' @inheritParams comono
#' @param Z Numeric matrix where each column represents a function on vertices.
#'   Number of rows must equal the number of vertices (\code{length(adj.list)}).
#'
#' @return A list with class \code{"comono_matrix_result"} containing:
#'   \item{column.coefficients}{List of numeric vectors, one per column of Z.
#'     Each vector contains vertex-wise coefficients.}
#'   \item{mean.coefficients}{Numeric vector of length \code{ncol(Z)} containing
#'     mean co-monotonicity for each column.}
#'   \item{median.coefficients}{Numeric vector of length \code{ncol(Z)} containing
#'     median co-monotonicity for each column.}
#'   \item{n.positive}{Integer vector of length \code{ncol(Z)} containing count
#'     of vertices with positive co-monotonicity for each column.}
#'   \item{n.negative}{Integer vector of length \code{ncol(Z)} containing count
#'     of vertices with negative co-monotonicity for each column.}
#'   \item{n.zero}{Integer vector of length \code{ncol(Z)} containing count
#'     of vertices with zero co-monotonicity for each column.}
#'   \item{type}{Weighting scheme used.}
#'   \item{n.vertices}{Total number of vertices.}
#'   \item{n.columns}{Number of columns in Z.}
#'
#' @examples
#' \dontrun{
#' library(gflow)
#'
#' ## Feature selection example
#' ## Z contains multiple ASV abundances as candidate predictors
#' data("sPTB_microbiome")
#' y <- sPTB_data$outcome
#' Z <- sPTB_data$composition  # n x p matrix of abundances
#'
#' ## Build graph
#' graph <- build.knn.graph(X, k = 15)
#'
#' ## Compute co-monotonicity with each feature
#' result <- comono.matrix(graph$adj.list, graph$weight.list, y, Z,
#'                         type = "derivative")
#'
#' ## Identify most concordant features
#' feature.importance <- abs(result$mean.coefficients)
#' top.features <- order(feature.importance, decreasing = TRUE)[1:10]
#'
#' cat("Top 10 features by co-monotonicity:\n")
#' print(data.frame(
#'   Feature = top.features,
#'   CoMonotonicity = result$mean.coefficients[top.features]
#' ))
#'
#' ## Visualize spatial pattern for top feature
#' plot(graph, vertex.color = result$column.coefficients[[top.features[1]]],
#'      main = paste("Co-monotonicity: y vs Feature", top.features[1]))
#'
#'
#' ## Model comparison example
#' ## Compare multiple regression models
#' fit1 <- gflow.regression(graph$adj.list, graph$weight.list, y, lambda = 0.01)
#' fit2 <- gflow.regression(graph$adj.list, graph$weight.list, y, lambda = 0.1)
#' fit3 <- gflow.regression(graph$adj.list, graph$weight.list, y, lambda = 1.0)
#'
#' predictions <- cbind(fit1$fitted, fit2$fitted, fit3$fitted)
#' colnames(predictions) <- c("lambda=0.01", "lambda=0.1", "lambda=1.0")
#'
#' result <- comono.matrix(graph$adj.list, graph$weight.list, y, predictions)
#'
#' ## Compare models
#' model.comparison <- data.frame(
#'   Model = colnames(predictions),
#'   MeanComono = result$mean.coefficients,
#'   MedianComono = result$median.coefficients,
#'   PropPositive = result$n.positive / result$n.vertices
#' )
#' print(model.comparison)
#'
#' ## Best model has highest mean co-monotonicity
#' best.model <- which.max(result$mean.coefficients)
#' cat("Best model:", colnames(predictions)[best.model], "\n")
#'
#'
#' ## Diagnostic example: check residual structure
#' fit <- gflow.regression(graph$adj.list, graph$weight.list, y, lambda = 0.1)
#' residuals <- y - fit$fitted
#'
#' ## Check co-monotonicity between residuals and covariates
#' result <- comono.matrix(graph$adj.list, graph$weight.list,
#'                         residuals, Z[, 1:5])  # First 5 covariates
#'
#' ## Significant co-monotonicity suggests model misspecification
#' for (j in 1:5) {
#'   if (abs(result$mean.coefficients[j]) > 0.2) {
#'     warning(sprintf("Residuals show co-monotonicity with covariate %d", j))
#'   }
#' }
#' }
#'
#' @references
#' Gajer, P. (2025). Morse-Smale regression: A geometric approach to
#' statistical analysis on graphs. \emph{In preparation}.
#'
#' @seealso
#' \code{\link{comono}} for scalar co-monotonicity,
#' \code{\link{comono.global}} for single summary value,
#' \code{\link{print.comono_matrix_result}} for printing results,
#' \code{\link{summary.comono_matrix_result}} for detailed summaries
#'
#' @export
comono.matrix <- function(adj.list, weight.list, y, Z,
                          type = c("unit", "derivative", "sign")) {
    type <- match.arg(type)

    ## Validation
    if (!is.matrix(Z)) {
        stop("Z must be a matrix")
    }

    n.vertices <- length(adj.list)

    if (length(y) != n.vertices) {
        stop("Length of y must match number of vertices")
    }

    if (nrow(Z) != n.vertices) {
        stop("Number of rows in Z must match number of vertices")
    }

    ## Convert to 0-based indexing
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    ## Call C++
    result <- .Call(S_comono_matrix,
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    Z,
                    type,
                    PACKAGE = "gflow")

    ## Add metadata
    result$type <- type
    result$n.vertices <- n.vertices
    result$n.columns <- ncol(Z)

    ## Set class
    class(result) <- c("comono_matrix_result", "list")

    return(result)
}

#' @export
print.comono_matrix_result <- function(x, digits = 4, ...) {
    cat("Co-monotonicity Matrix Results\n")
    cat("==============================\n\n")

    cat("Weighting scheme:", x$type, "\n")
    cat("Number of vertices:", x$n.vertices, "\n")
    cat("Number of columns:", x$n.columns, "\n\n")

    cat("Mean coefficients by column:\n")
    print(round(x$mean.coefficients, digits))

    invisible(x)
}
