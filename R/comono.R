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


#' Fast permutation test using C++ batch computation
comono.perm.test.fast <- function(adj.list, weight.list, y, z,
                                   n.perm = 1000, type = "derivative") {
    ## Observed
    obs <- comono(adj.list, weight.list, y, z, type = type)

    ## Create permutation matrix
    Z.perm <- replicate(n.perm, sample(z))

    ## Call C++ batch function (MUCH faster than loop)
    null.matrix <- .Call(S_comono_batch,
                         adj.list.0based,
                         weight.list,
                         as.numeric(y),
                         Z.perm,  # Matrix of permuted z vectors
                         type,
                         PACKAGE = "gflow")

    ## Compute p-values in R
    p.values <- sapply(1:length(y), function(v) {
        mean(abs(null.matrix[v, ]) >= abs(obs$vertex.coefficients[v]))
    })

    q.values <- p.adjust(p.values, method = "BH")

    list(observed = obs$vertex.coefficients,
         p.values = p.values,
         q.values = q.values,
         null.distribution = null.matrix)
}


#' Compute Correlation-Type Co-Monotonicity Coefficients
#'
#' Computes vertex-level co-monotonicity using Pearson-style normalization,
#' which measures correlation between directional changes. This variant is
#' robust to sparse signal patterns where one function varies on only a subset
#' of edges.
#'
#' @param adj.list List of integer vectors (1-based). Each element contains
#'   neighbor indices for a vertex.
#' @param weight.list List of numeric vectors. Each element contains edge
#'   weights/lengths corresponding to neighbors in adj.list.
#' @param y Numeric vector of response values (length = number of vertices).
#' @param z Numeric vector of feature values (length = number of vertices).
#' @param type Character string specifying weighting scheme:
#'   \itemize{
#'     \item \code{"unit"}: Equal weight (w = 1) for all edges
#'     \item \code{"derivative"}: Geometric weight (w = 1/length^2) for edges
#'   }
#'   Default is "derivative".
#'
#' @return Named list with components:
#'   \item{vertex.coefficients}{Numeric vector of correlation-type
#'     co-monotonicity coefficients, one per vertex}
#'   \item{mean.coefficient}{Mean of vertex coefficients}
#'   \item{median.coefficient}{Median of vertex coefficients}
#'   \item{counts}{Named integer vector with n.positive, n.negative, n.zero}
#'
#' @details
#' The correlation-type co-monotonicity coefficient at vertex v is:
#' \deqn{cm_{cor}(y,z)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z}
#'   {\sqrt{\sum w_e (\Delta_e y)^2} \sqrt{\sum w_e (\Delta_e z)^2}}}
#'
#' This normalization separates the variation in y and z, making the measure
#' naturally robust to sparse signals. If z varies on only k out of n edges,
#' the coefficient scales approximately by sqrt(k/n), automatically down-weighting
#' associations with sparse support.
#'
#' For smooth functions on manifolds with derivative weighting, this measure
#' converges to cos(theta), where theta is the angle between gradient vectors.
#'
#' @section When to Use:
#' \itemize{
#'   \item Default choice for most applications
#'   \item Continuous functions with real variation
#'   \item When gradient correlation is the scientific question
#'   \item Transformed microbiome data (e.g., log-transformed abundances)
#' }
#'
#' @section Comparison to Other Variants:
#' \itemize{
#'   \item vs. \code{comono()}: cm_cor has smaller magnitude when signal is
#'     sparse; more robust to heterogeneous feature presence
#'   \item vs. \code{comono.proportion()}: cm_cor is magnitude-weighted;
#'     proportion is count-based
#' }
#'
#' @examples
#' \dontrun{
#' ## Create simple graph
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
#'
#' ## Response and feature values
#' y <- c(0.0, 1.0, 0.5)
#' z <- c(0.0, 0.8, 0.6)
#'
#' ## Compute correlation-type coefficients
#' result <- comono.cor(adj.list, weight.list, y, z, type = "derivative")
#'
#' print(result$vertex.coefficients)
#' print(result$mean.coefficient)
#' }
#'
#' @seealso \code{\link{comono.proportion}} for threshold-based variant,
#'   \code{\link{comono}} for absolute-value normalization
#'
#' @export
comono.cor <- function(adj.list, weight.list, y, z, type = "derivative") {

  ## Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }
  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }
  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have same length")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.numeric(z)) {
    stop("z must be numeric")
  }
  if (length(y) != length(adj.list)) {
    stop("Length of y must equal length of adj.list")
  }
  if (length(z) != length(adj.list)) {
    stop("Length of z must equal length of adj.list")
  }
  if (!type %in% c("unit", "derivative")) {
    stop("type must be 'unit' or 'derivative'")
  }

  ## Convert to 0-based indexing for C++
  adj.list.0 <- lapply(adj.list, function(x) {
    if (length(x) == 0) return(integer(0))
    as.integer(x - 1)
  })

  ## Call C++ function
  result <- .Call(
    S_comono_cor,
    adj.list.0,
    weight.list,
    y,
    z,
    type,
    PACKAGE = "gflow"
  )

  ## Add S3 class
  class(result) <- c("comono_result", "list")

  return(result)
}

#' Compute Thresholded Proportion-of-Agreements Co-Monotonicity
#'
#' Computes vertex-level co-monotonicity by counting directional agreements
#' among edges where both functions exhibit meaningful variation above specified
#' thresholds. Divides by total number of edges, so sparse signal naturally
#' produces values near zero.
#'
#' @param adj.list List of integer vectors (1-based). Each element contains
#'   neighbor indices for a vertex.
#' @param weight.list List of numeric vectors. Each element contains edge
#'   weights/lengths corresponding to neighbors in adj.list.
#' @param y Numeric vector of response values (length = number of vertices).
#' @param z Numeric vector of feature values (length = number of vertices).
#' @param tau.y Numeric scalar. Threshold for |Delta_y|. Edges with |Delta_y| <= tau.y
#'   contribute zero. Default 0.0 (no filtering).
#' @param tau.z Numeric scalar. Threshold for |Delta_z|. Edges with |Delta_z| <= tau.z
#'   contribute zero. Default 0.0 (no filtering).
#'
#' @return Named list with components:
#'   \item{vertex.coefficients}{Numeric vector of proportion-based
#'     co-monotonicity coefficients, one per vertex}
#'   \item{mean.coefficient}{Mean of vertex coefficients}
#'   \item{median.coefficient}{Median of vertex coefficients}
#'   \item{counts}{Named integer vector with n.positive, n.negative, n.zero}
#'
#' @details
#' Each edge receives score +1 (agreement), -1 (disagreement), or 0 (insufficient
#' signal). The coefficient at vertex v is the average score across ALL edges:
#' \deqn{cm_{prop}(y,z)(v) = \frac{n_{agree} - n_{disagree}}{|N(v)|}}
#'
#' where n_agree and n_disagree count edges passing both thresholds with
#' concordant/discordant directions, and |N(v)| is the total number of edges
#' (not just those passing thresholds).
#'
#' @section Threshold Selection:
#' For response y:
#' \itemize{
#'   \item Absolute: \code{tau.y = 0.01 * diff(range(y))}
#'   \item Relative: \code{tau.y = 0.05 * sd(y)}
#'   \item Data-driven: \code{tau.y = quantile(abs(diff.y), 0.10)}
#' }
#'
#' For feature z (adaptive per-feature recommended):
#' \itemize{
#'   \item \code{tau.z = quantile(abs(delta.z), 0.25)} filters bottom 25\% as noise
#' }
#'
#' @section When to Use:
#' \itemize{
#'   \item Sparse/binary features (e.g., phylotype presence/absence)
#'   \item When signal filtering is scientifically justified
#'   \item Noisy data where subthreshold changes should be ignored
#'   \item When "prevalence of agreement" is the scientific question
#' }
#'
#' @section Comparison to Other Variants:
#' \itemize{
#'   \item vs. \code{comono.cor()}: proportion is count-based, cor is
#'     magnitude-weighted
#'   \item vs. \code{comono()}: proportion divides by all edges, comono divides
#'     by valid edges only
#' }
#'
#' @examples
#' \dontrun{
#' ## Sparse feature example
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight.list <- list(c(1.0, 1.0), c(1.0, 1.0), c(1.0, 1.0))
#'
#' y <- c(0.0, 1.0, 0.5)
#' z <- c(0.0, 0.001, 0.0005)  ## Sparse signal
#'
#' ## Set thresholds
#' tau.y <- 0.05 * sd(y)
#' tau.z <- quantile(abs(z), 0.25)
#'
#' ## Compute proportion-based coefficients
#' result <- comono.proportion(adj.list, weight.list, y, z, tau.y, tau.z)
#'
#' ## Sparse signal produces small coefficients
#' print(result$vertex.coefficients)
#' }
#'
#' @seealso \code{\link{comono.cor}} for correlation-type variant,
#'   \code{\link{comono}} for absolute-value normalization
#'
#' @export
comono.proportion <- function(adj.list, weight.list, y, z,
                               tau.y = 0.0, tau.z = 0.0) {

  ## Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }
  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }
  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have same length")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.numeric(z)) {
    stop("z must be numeric")
  }
  if (length(y) != length(adj.list)) {
    stop("Length of y must equal length of adj.list")
  }
  if (length(z) != length(adj.list)) {
    stop("Length of z must equal length of adj.list")
  }
  if (!is.numeric(tau.y) || length(tau.y) != 1 || tau.y < 0) {
    stop("tau.y must be a non-negative numeric scalar")
  }
  if (!is.numeric(tau.z) || length(tau.z) != 1 || tau.z < 0) {
    stop("tau.z must be a non-negative numeric scalar")
  }

  ## Convert to 0-based indexing for C++
  adj.list.0 <- lapply(adj.list, function(x) {
    if (length(x) == 0) return(integer(0))
    as.integer(x - 1)
  })

  ## Call C++ function
  result <- .Call(
    S_comono_proportion,
    adj.list.0,
    weight.list,
    y,
    z,
    tau.y,
    tau.z,
    PACKAGE = "gflow"
  )

  ## Add S3 class
  class(result) <- c("comono_result", "list")

  return(result)
}

#' Helper Function: Compute Adaptive Thresholds for Phylotype Co-Monotonicity
#'
#' Computes per-phylotype adaptive thresholds based on quantiles of edge
#' differences. Useful for filtering noise in sparse biological data.
#'
#' @param X Numeric matrix of phylotype abundances (samples x phylotypes).
#' @param adj.list List of neighbor indices (1-based).
#' @param quantile.cutoff Quantile for threshold (default 0.25 = filter bottom 25\%).
#'
#' @return Numeric vector of thresholds, one per phylotype (column of X).
#'
#' @export
compute.adaptive.thresholds <- function(X, adj.list, quantile.cutoff = 0.25) {

  n.phylotypes <- ncol(X)
  thresholds <- numeric(n.phylotypes)

  for (j in 1:n.phylotypes) {
    z <- X[, j]

    ## Collect all |delta_z| values for this phylotype
    all.deltas <- c()
    for (i in 1:length(adj.list)) {
      neighbors <- adj.list[[i]]
      if (length(neighbors) > 0) {
        deltas <- abs(z[neighbors] - z[i])
        all.deltas <- c(all.deltas, deltas)
      }
    }

    ## Set threshold at specified quantile
    if (length(all.deltas) > 0) {
      thresholds[j] <- quantile(all.deltas, quantile.cutoff)
    } else {
      thresholds[j] <- 0.0
    }
  }

  return(thresholds)
}

#' Compute Correlation-Type Co-Monotonicity Matrix
#'
#' Efficiently computes correlation-type co-monotonicity between a response
#' vector and multiple feature vectors (columns of a matrix). More efficient
#' than repeated calls to \code{comono.cor} because edge differences for y
#' are computed only once.
#'
#' @param adj.list List of integer vectors (1-based). Each element contains
#'   neighbor indices for a vertex.
#' @param weight.list List of numeric vectors. Each element contains edge
#'   weights/lengths corresponding to neighbors in adj.list.
#' @param y Numeric vector of response values (length = number of vertices).
#' @param Z Numeric matrix of feature values (vertices x features). Each column
#'   represents one feature.
#' @param type Character string specifying weighting scheme:
#'   \itemize{
#'     \item \code{"unit"}: Equal weight (w = 1) for all edges
#'     \item \code{"derivative"}: Geometric weight (w = 1/length^2) for edges
#'   }
#'   Default is "derivative".
#'
#' @return Named list with components:
#'   \item{column.coefficients}{Numeric matrix (n x m) of vertex-level
#'     coefficients, where n = number of vertices and m = number of features}
#'   \item{column.means}{Numeric vector of length m with mean coefficient per feature}
#'   \item{column.medians}{Numeric vector of length m with median coefficient per feature}
#'   \item{column.counts}{Numeric matrix (m x 3) with counts of positive, negative,
#'     and zero coefficients for each feature}
#'
#' @details
#' For each feature (column j of Z), computes the correlation-type co-monotonicity
#' coefficient at each vertex v:
#' \deqn{cm_{cor}(y,z_j)(v) = \frac{\sum w_e \Delta_e y \cdot \Delta_e z_j}
#'   {\sqrt{\sum w_e (\Delta_e y)^2} \sqrt{\sum w_e (\Delta_e z_j)^2}}}
#'
#' The resulting matrix has dimensions n x m, where entry (v,j) is the
#' co-monotonicity between y and feature j at vertex v.
#'
#' @section Computational Efficiency:
#' Computing co-monotonicity for m features requires:
#' \itemize{
#'   \item Separate calls: O(m * |E|) where |E| is total number of edges
#'   \item Matrix version: O(|E| + m * |E|) with shared preprocessing
#' }
#' The matrix version is approximately 2x faster for large m.
#'
#' @examples
#' \dontrun{
#' ## Microbiome data example
#' adj.list <- ... ## k-NN graph
#' weight.list <- ... ## edge lengths
#' y <- sPTB.response  ## spontaneous preterm birth prevalence
#' Z <- phylotype.abundances  ## n x p matrix of phylotype abundances
#'
#' ## Compute correlation-type co-monotonicity for all phylotypes
#' result <- comono.cor.matrix(adj.list, weight.list, y, Z, type = "derivative")
#'
#' ## Extract coefficients for phylotype j
#' phylo.j.coeffs <- result$column.coefficients[, j]
#'
#' ## Find phylotypes with strong positive association (mean > 0.5)
#' strong.positive <- which(result$column.means > 0.5)
#' }
#'
#' @seealso \code{\link{comono.cor}} for single-feature version,
#'   \code{\link{comono.proportion.matrix}} for threshold-based variant
#'
#' @export
comono.cor.matrix <- function(adj.list, weight.list, y, Z, type = "derivative") {

  ## Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }
  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }
  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have same length")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("Z must be a numeric matrix")
  }
  if (length(y) != nrow(Z)) {
    stop("Length of y must equal number of rows in Z")
  }
  if (length(y) != length(adj.list)) {
    stop("Length of y must equal length of adj.list")
  }
  if (!type %in% c("unit", "derivative")) {
    stop("type must be 'unit' or 'derivative'")
  }

  ## Convert to 0-based indexing for C++
  adj.list.0 <- lapply(adj.list, function(x) {
    if (length(x) == 0) return(integer(0))
    as.integer(x - 1)
  })

  ## Call C++ function
  result <- .Call(
    S_comono_cor_matrix,
    adj.list.0,
    weight.list,
    y,
    Z,
    type,
    PACKAGE = "gflow"
  )

  ## Preserve Z column names if present
  if (!is.null(colnames(Z))) {
    colnames(result$column.coefficients) <- colnames(Z)
    names(result$column.means) <- colnames(Z)
    names(result$column.medians) <- colnames(Z)
    rownames(result$column.counts) <- colnames(Z)
  }

  ## Add S3 class
  class(result) <- c("comono_matrix_result", "list")

  return(result)
}

#' Compute Thresholded Proportion Co-Monotonicity Matrix
#'
#' Efficiently computes threshold-filtered proportion-based co-monotonicity
#' between a response vector and multiple feature vectors. Accepts per-feature
#' adaptive thresholds for optimal noise filtering.
#'
#' @param adj.list List of integer vectors (1-based). Each element contains
#'   neighbor indices for a vertex.
#' @param weight.list List of numeric vectors. Each element contains edge
#'   weights/lengths corresponding to neighbors in adj.list.
#' @param y Numeric vector of response values (length = number of vertices).
#' @param Z Numeric matrix of feature values (vertices x features). Each column
#'   represents one feature.
#' @param tau.y Numeric scalar. Threshold for |Delta_y|. Edges with |Delta_y| <= tau.y
#'   contribute zero. Default 0.0 (no filtering).
#' @param tau.z Numeric vector of length m (number of features) OR scalar.
#'   Thresholds for |Delta_z| per feature. If scalar, same threshold is used
#'   for all features. Default 0.0 (no filtering).
#'
#' @return Named list with components:
#'   \item{column.coefficients}{Numeric matrix (n x m) of vertex-level
#'     coefficients}
#'   \item{column.means}{Numeric vector of length m with mean coefficient per feature}
#'   \item{column.medians}{Numeric vector of length m with median coefficient per feature}
#'   \item{column.counts}{Numeric matrix (m x 3) with counts of positive, negative,
#'     and zero coefficients for each feature}
#'
#' @details
#' For each feature (column j of Z), computes the proportion-based co-monotonicity
#' at each vertex v by counting directional agreements among edges where both
#' |Delta_y| > tau.y and |Delta_z_j| > tau.z[j].
#'
#' Per-feature adaptive thresholds (vector tau.z) are recommended for microbiome
#' data where different phylotypes have different abundance scales and noise levels.
#'
#' @section Threshold Selection:
#' For optimal results with microbiome data:
#' \enumerate{
#'   \item Set tau.y based on response variation:
#'     \code{tau.y <- 0.05 * sd(y)}
#'   \item Compute adaptive per-feature thresholds:
#'     \code{tau.z <- compute.adaptive.thresholds(Z, adj.list, quantile.cutoff = 0.25)}
#'   \item Pass vector tau.z to this function
#' }
#'
#' @examples
#' \dontrun{
#' ## Sparse microbiome data example
#' adj.list <- ... ## k-NN graph
#' weight.list <- ... ## edge lengths
#' y <- sPTB.response
#' Z <- phylotype.abundances  ## Many near-zero values
#'
#' ## Set thresholds
#' tau.y <- 0.05 * sd(y)
#' tau.z <- compute.adaptive.thresholds(Z, adj.list, quantile.cutoff = 0.25)
#'
#' ## Compute proportion-based coefficients
#' result <- comono.proportion.matrix(adj.list, weight.list, y, Z,
#'                                    tau.y, tau.z)
#'
#' ## Features with sparse signal will have small mean coefficients
#' hist(result$column.means)
#'
#' ## Identify features with prevalent directional concordance (mean > 0.2)
#' concordant.features <- which(result$column.means > 0.2)
#' }
#'
#' @seealso \code{\link{comono.proportion}} for single-feature version,
#'   \code{\link{comono.cor.matrix}} for correlation-type variant,
#'   \code{\link{compute.adaptive.thresholds}} for threshold selection
#'
#' @export
comono.proportion.matrix <- function(adj.list, weight.list, y, Z,
                                     tau.y = 0.0, tau.z = 0.0) {

  ## Input validation
  if (!is.list(adj.list)) {
    stop("adj.list must be a list")
  }
  if (!is.list(weight.list)) {
    stop("weight.list must be a list")
  }
  if (length(adj.list) != length(weight.list)) {
    stop("adj.list and weight.list must have same length")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("Z must be a numeric matrix")
  }
  if (length(y) != nrow(Z)) {
    stop("Length of y must equal number of rows in Z")
  }
  if (length(y) != length(adj.list)) {
    stop("Length of y must equal length of adj.list")
  }
  if (!is.numeric(tau.y) || length(tau.y) != 1 || tau.y < 0) {
    stop("tau.y must be a non-negative numeric scalar")
  }
  if (!is.numeric(tau.z) || any(tau.z < 0)) {
    stop("tau.z must be non-negative numeric (scalar or vector)")
  }
  if (length(tau.z) != 1 && length(tau.z) != ncol(Z)) {
    stop("tau.z must be scalar or vector of length ncol(Z)")
  }

  ## Convert to 0-based indexing for C++
  adj.list.0 <- lapply(adj.list, function(x) {
    if (length(x) == 0) return(integer(0))
    as.integer(x - 1)
  })

  ## Call C++ function
  result <- .Call(
    S_comono_proportion_matrix,
    adj.list.0,
    weight.list,
    y,
    Z,
    tau.y,
    tau.z,
    PACKAGE = "gflow"
  )

  ## Preserve Z column names if present
  if (!is.null(colnames(Z))) {
    colnames(result$column.coefficients) <- colnames(Z)
    names(result$column.means) <- colnames(Z)
    names(result$column.medians) <- colnames(Z)
    rownames(result$column.counts) <- colnames(Z)
  }

  ## Add S3 class
  class(result) <- c("comono_matrix_result", "list")

  return(result)
}


#' Find Vertex-Phylotype Pairs with Discrepant Unit vs Derivative Correlation Coefficients
#'
#' Identifies vertex-phylotype pairs where correlation-type co-monotonicity
#' differs substantially between unit and derivative weighting. Used for
#' debugging and understanding when geometric normalization matters.
#'
#' @param cm.unit Result from comono.cor.matrix() with type="unit"
#' @param cm.deriv Result from comono.cor.matrix() with type="derivative"
#' @param partition Optional vector of cell assignments for vertices (for filtering)
#' @param cell.id Optional cell ID to restrict analysis (requires partition)
#' @param n.pairs Number of discrepant pairs to extract (default 10)
#' @param diff.threshold Minimum absolute difference to consider (default 0.75)
#'
#' @return List with components:
#'   \item{discrepant.pairs}{Data frame with vertex.idx, phylotype.idx, phylotype.name,
#'     unit.value, deriv.value, difference, case.type}
#'   \item{vertex.idx}{Vector of vertex indices (0-based for C++)}
#'   \item{phylotype.idx}{Vector of phylotype indices (0-based for C++)}
#'   \item{csv.path}{Path to saved CSV file for C++ debugging}
#'
#' @details
#' Saves a CSV file to /tmp/comono_debugging_dir/comono_cor_debug_pairs.csv
#' that will be loaded by the C++ debugging code when DEBUG_COMONO_COR is enabled.
#'
#' The case.type categorizes discrepancies:
#' \itemize{
#'   \item "opposite_sign": Unit and derivative have opposite signs
#'   \item "large_positive": Both positive, large difference
#'   \item "large_negative": Both negative, large difference
#'   \item "unit_larger": Unit coefficient has larger magnitude
#'   \item "deriv_larger": Derivative coefficient has larger magnitude
#' }
#'
#' @examples
#' \dontrun{
#' ## Compute both types
#' cm.unit <- comono.cor.matrix(adj.list, weight.list, y, Z, type = "unit")
#' cm.deriv <- comono.cor.matrix(adj.list, weight.list, y, Z, type = "derivative")
#'
#' ## Find discrepant pairs
#' discrepant <- find.discrepant.pairs.cor(
#'   cm.unit = cm.unit,
#'   cm.deriv = cm.deriv,
#'   n.pairs = 10,
#'   diff.threshold = 0.75
#' )
#'
#' ## Now recompile with DEBUG_COMONO_COR = 1 and run specific phylotypes
#' phylo.idx <- discrepant$phylotype.idx[1]
#' z.values <- Z[, phylo.idx]
#'
#' result.unit <- comono.cor(adj.list, weight.list, y, z.values, type = "unit")
#' result.deriv <- comono.cor(adj.list, weight.list, y, z.values, type = "derivative")
#' }
#'
#' @export
find.discrepant.pairs.cor <- function(cm.unit,
                                      cm.deriv,
                                      partition = NULL,
                                      cell.id = NULL,
                                      n.pairs = 10,
                                      diff.threshold = 0.75) {

  ## Input validation
  if (!inherits(cm.unit, "comono_matrix_result")) {
    stop("cm.unit must be result from comono.cor.matrix()")
  }
  if (!inherits(cm.deriv, "comono_matrix_result")) {
    stop("cm.deriv must be result from comono.cor.matrix()")
  }

  ## Extract coefficient matrices
  unit.mat <- cm.unit$column.coefficients
  deriv.mat <- cm.deriv$column.coefficients

  if (!identical(dim(unit.mat), dim(deriv.mat))) {
    stop("Unit and derivative matrices must have same dimensions")
  }

  n.vertices <- nrow(unit.mat)
  n.phylotypes <- ncol(unit.mat)

  ## Filter by cell if requested
  vertex.mask <- rep(TRUE, n.vertices)
  if (!is.null(partition) && !is.null(cell.id)) {
    if (length(partition) != n.vertices) {
      stop("partition length must equal number of vertices")
    }
    vertex.mask <- (partition == cell.id)
    cat(sprintf("Filtering to cell %d: %d vertices\n", cell.id, sum(vertex.mask)))
  }

  ## Compute differences
  diff.mat <- abs(unit.mat - deriv.mat)

  ## Find pairs exceeding threshold within filtered vertices
  candidate.pairs <- data.frame(
    vertex.idx = integer(),
    phylotype.idx = integer(),
    unit.value = numeric(),
    deriv.value = numeric(),
    difference = numeric()
  )

  for (j in 1:n.phylotypes) {
    for (v in which(vertex.mask)) {
      if (diff.mat[v, j] >= diff.threshold) {
        candidate.pairs <- rbind(candidate.pairs, data.frame(
          vertex.idx = v,
          phylotype.idx = j,
          unit.value = unit.mat[v, j],
          deriv.value = deriv.mat[v, j],
          difference = diff.mat[v, j]
        ))
      }
    }
  }

  if (nrow(candidate.pairs) == 0) {
    cat(sprintf("No pairs found with difference >= %.2f\n", diff.threshold))
    return(NULL)
  }

  ## Sort by difference (descending)
  candidate.pairs <- candidate.pairs[order(-candidate.pairs$difference), ]

  ## Take top n.pairs
  n.extract <- min(n.pairs, nrow(candidate.pairs))
  top.pairs <- candidate.pairs[1:n.extract, ]

  ## Add phylotype names if available
  if (!is.null(colnames(unit.mat))) {
    top.pairs$phylotype.name <- colnames(unit.mat)[top.pairs$phylotype.idx]
  } else {
    top.pairs$phylotype.name <- paste0("Phylotype_", top.pairs$phylotype.idx)
  }

  ## Categorize discrepancy type
  top.pairs$case.type <- sapply(1:nrow(top.pairs), function(i) {
    u <- top.pairs$unit.value[i]
    d <- top.pairs$deriv.value[i]
    diff <- top.pairs$difference[i]

    if (sign(u) != sign(d)) {
      return("opposite_sign")
    } else if (u > 0 && d > 0) {
      if (abs(u) > abs(d)) {
        return("unit_larger_positive")
      } else {
        return("deriv_larger_positive")
      }
    } else if (u < 0 && d < 0) {
      if (abs(u) > abs(d)) {
        return("unit_larger_negative")
      } else {
        return("deriv_larger_negative")
      }
    } else {
      return("near_zero")
    }
  })

  ## Create output directory
  debug.dir <- "/tmp/comono_debugging_dir"
  if (!dir.exists(debug.dir)) {
    dir.create(debug.dir, recursive = TRUE)
  }

  ## Save CSV for C++ debugging
  csv.path <- file.path(debug.dir, "comono_cor_debug_pairs.csv")

  ## Prepare output (convert to 0-based for C++)
  output.df <- data.frame(
    vertex_idx = top.pairs$vertex.idx - 1,  # 0-based
    phylotype_idx = top.pairs$phylotype.idx - 1,  # 0-based
    phylotype_name = top.pairs$phylotype.name,
    unit_value = top.pairs$unit.value,
    deriv_value = top.pairs$deriv.value,
    difference = top.pairs$difference,
    case = top.pairs$case.type
  )

  write.csv(output.df, csv.path, row.names = FALSE, quote = TRUE)

  ## Print summary
  cat("\n========================================\n")
  cat(sprintf("Found %d discrepant pairs (threshold >= %.2f)\n",
              nrow(candidate.pairs), diff.threshold))
  cat(sprintf("Saved top %d pairs to:\n%s\n", n.extract, csv.path))
  cat("========================================\n\n")

  cat("Summary by case type:\n")
  print(table(output.df$case))
  cat("\n")

  cat("Top 5 pairs:\n")
  print(output.df[1:min(5, nrow(output.df)),
                  c("vertex_idx", "phylotype_name", "unit_value",
                    "deriv_value", "difference", "case")])
  cat("\n")

  ## Return result
  result <- list(
    discrepant.pairs = output.df,
    vertex.idx = output.df$vertex_idx,
    phylotype.idx = output.df$phylotype_idx,
    csv.path = csv.path
  )

  class(result) <- "discrepant_pairs_cor"
  return(result)
}

#' Compute Scale-Consistent Co-Monotonicity Coefficients
#'
#' Computes correlation-type co-monotonicity followed by geometric smoothing
#' to remove scale artifacts from edge length heterogeneity. This ensures
#' associations are measured at the same geometric scale as the response.
#'
#' @param adj.list Graph adjacency list (1-based indexing)
#' @param weight.list Edge weight/length list
#' @param y Response values (should be pre-smoothed)
#' @param Z Feature matrix (vertices x features)
#' @param type Weighting type: "unit" or "derivative"
#' @param smooth.fit Fitted rdgraph regression object from response smoothing
#'   (output of rdcx.regression or similar). If NULL, no smoothing applied.
#'
#' @return List with:
#'   - raw.coefficients: n x m matrix of unsmoothed coefficients
#'   - smoothed.coefficients: n x m matrix after geometric filtering (if smooth.fit provided)
#'   - column.means/medians/counts: Summary statistics on smoothed coefficients
#'
#' @details
#' This function implements the recommended workflow for scale-consistent
#' association analysis: compute co-monotonicity at the vertex level, then
#' apply the same geometric smoothing used for the response variable. This
#' removes micro-scale artifacts while preserving biologically meaningful
#' large-scale association patterns.
#'
#' @examples
#' \dontrun{
#' ## First smooth response
#' y.fit <- rdcx.regression(adj.list, weight.list, y.raw, lambda = 0.1)
#' y.smooth <- y.fit$fitted.values
#'
#' ## Compute scale-consistent associations
#' cm.result <- comono.cor.smooth(
#'   adj.list, weight.list,
#'   y = y.smooth,
#'   Z = Z.raw,
#'   type = "derivative",
#'   smooth.fit = y.fit
#' )
#'
#' ## Smoothed coefficients ready for downstream analysis
#' associations <- cm.result$smoothed.coefficients
#' }
#'
#' @export
comono.cor.smooth <- function(adj.list, weight.list, y, Z,
                               type = "derivative",
                               smooth.fit = NULL) {

  ## Compute raw co-monotonicity
  cm.raw <- comono.cor.matrix(adj.list, weight.list, y, Z, type)

  result <- list(
    raw.coefficients = cm.raw$column.coefficients,
    raw.means = cm.raw$column.means,
    raw.medians = cm.raw$column.medians,
    raw.counts = cm.raw$column.counts
  )

  ## Apply smoothing if requested
  if (!is.null(smooth.fit)) {
    cat("Smoothing co-monotonicity coefficients using provided filter...\n")

    ## Smooth each phylotype's coefficient map
    cm.smooth.fit <- refit.rdgraph.regression(
      smooth.fit,
      cm.raw$column.coefficients
    )

    result$smoothed.coefficients <- cm.smooth.fit$fitted.values

    ## Recompute summary statistics on smoothed coefficients
    n.vertices <- nrow(result$smoothed.coefficients)
    n.features <- ncol(result$smoothed.coefficients)

    result$smoothed.means <- colMeans(result$smoothed.coefficients)
    result$smoothed.medians <- apply(result$smoothed.coefficients, 2, median)

    ## Count positive/negative/zero (using same threshold as raw)
    epsilon <- 1e-10
    result$smoothed.counts <- matrix(0, nrow = n.features, ncol = 3)
    colnames(result$smoothed.counts) <- c("n.positive", "n.negative", "n.zero")

    for (j in 1:n.features) {
      coeffs <- result$smoothed.coefficients[, j]
      result$smoothed.counts[j, 1] <- sum(coeffs > epsilon)
      result$smoothed.counts[j, 2] <- sum(coeffs < -epsilon)
      result$smoothed.counts[j, 3] <- sum(abs(coeffs) <= epsilon)
    }

    if (!is.null(colnames(Z))) {
      rownames(result$smoothed.counts) <- colnames(Z)
    }

    cat(sprintf("Smoothing complete: %d vertices x %d features\n",
                n.vertices, n.features))
  } else {
    cat("No smoothing requested - returning raw coefficients only\n")
  }

  class(result) <- c("comono_smooth_result", "list")
  return(result)
}
