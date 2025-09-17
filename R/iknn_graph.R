#' Create a Single Intersection-weighted k-Nearest Neighbors Graph
#'
#' @description
#' Creates and prunes an intersection-weighted k-nearest neighbors (IWD-kNN) graph from
#' input data. The graph is constructed based on feature space distances and intersection
#' sizes between neighbor sets.
#'
#' @param X A numeric matrix where rows represent data points and columns represent features.
#'        Data frames will be coerced to matrices.
#' @param k Integer scalar. The number of nearest neighbors to consider for each point.
#'        Must be positive and less than the number of data points.
#' @param pruning.thld numeric, controlling intensity of the geometric edge pruning.
#'         Edge weight relative deviation is computed as
#'         rel_dev = (w_ij + w_jk) / w_ik - 1.0;
#'         Geometric pruning is performed on all edges with rel_dev < pruning_thld
#' @param compute.full Logical scalar. If TRUE, computes additional graph components and metrics.
#'        If FALSE, computes only essential components. Default value: FALSE.
#' @param pca.dim Maximum number of principal components to use if dimensionality reduction
#'        is applied (default: 100). Set to NULL to skip dimensionality reduction.
#' @param variance.explained Percentage of variance to be explained by the principal components
#'        (default: 0.99). If this threshold can be met with fewer components than pca.dim,
#'        the smaller number will be used. Set to NULL to use exactly pca.dim components.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#'
#' @return An object of class "IkNN" (inheriting from "list") containing:
#' \describe{
#'   \item{pruned_adj_list}{Adjacency lists after pruning (1-based indices)}
#'   \item{pruned_weight_list}{Distances for edges in pruned graph}
#'   \item{n_edges}{Total number of edges in original graph}
#'   \item{n_edges_in_pruned_graph}{Number of edges after pruning}
#'   \item{n_removed_edges}{Number of edges removed by pruning}
#'   \item{edge_reduction_ratio}{Proportion of edges removed}
#'   \item{call}{The matched function call}
#'   \item{k}{Number of nearest neighbors used}
#'   \item{n}{Number of data points}
#'   \item{d}{Number of features}
#' }
#'
#' If compute.full = TRUE, additional components include:
#' \describe{
#'   \item{adj_list}{Original adjacency lists (1-based indices)}
#'   \item{isize_list}{Intersection sizes for original edges}
#'   \item{weight_list}{Distances for original edges}
#'   \item{conn_comps}{Connected components identification}
#'   \item{connected_components}{Alternative format of connected components}
#' }
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 2)
#' result <- create.single.iknn.graph(X, k = 3)
#' summary(result)
#' @export
create.single.iknn.graph <- function(X,
                                     k,
                                     pruning.thld = 0.1,
                                     compute.full = FALSE,
                                     pca.dim = 100,
                                     variance.explained = 0.99,
                                     verbose = TRUE) {
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    n <- nrow(X)
    if (n < 2) {
        stop("X must contain at least 2 data points")
    }
    storage.mode(X) <- "double"

    if (!is.numeric(k) || length(k) != 1 || k != round(k) || k < 1 || k >= n) {
        stop("k must be a positive integer less than the number of data points")
    }

    if (!is.numeric(pruning.thld) || length(pruning.thld) != 1 || pruning.thld <= 0 || pruning.thld >= 0.2) {
        stop("pruning.thld must be a positive numeric value less then 0.2")
    }

    if (!is.logical(compute.full) || length(compute.full) != 1) {
        stop("compute.full must be a single logical value")
    }

    ## Check pca.dim
    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim)) {
            stop("pca.dim must be a positive integer or NULL")
        }
    }

    ## Check variance.explained
    if (!is.null(variance.explained)) {
        if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
            variance.explained <= 0 || variance.explained > 1) {
            stop("variance.explained must be a numeric value between 0 and 1, or NULL")
        }
    }

    ## PCA dimensionality reduction if needed
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose) {
            message("High-dimensional data detected. Performing dimensionality reduction.")
        }

        ## Store original dimensions for reporting
        original_dim <- ncol(X)

        ## If variance.explained is specified, analyze variance
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(X,
                                                   variance.threshold = variance.explained,
                                                   max.components = pca.dim)

            ## Number of components to use (based on variance or pca.dim)
            n_components <- pca_analysis$n.components

            if (verbose) {
                message(sprintf("Using %d principal components (explains %.2f%% of variance)",
                                n_components, pca_analysis$variance.explained * 100))
            }

            ## Project data onto selected components
            X <- pca.project(X, pca_analysis$pca.result, n_components)

            ## Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance,
                pca_projected_X = X
            )
        } else {
            ## Use fixed number of components (pca.dim)
            if (verbose) {
                message(sprintf("Projecting data onto first %d principal components", pca.dim))
            }

            ## Perform PCA and projection
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)

            ## Calculate variance explained by pca.dim components
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)

            ## Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained,
                pca_projected_X = X
            )
        }
    }

    result <- .Call(S_create_single_iknn_graph,
                    X,
                    as.integer(k + 1),
                    as.double(pruning.thld),
                    as.logical(compute.full),
                    PACKAGE = "gflow")

    attr(result, "k") <- k
    attr(result, "pruning_threshold") <- pruning.thld
    attr(result, "call") <- match.call()

                                        # Add PCA-related attributes if PCA was performed
    if (!is.null(pca_info)) {
        attr(result, "pca") <- pca_info
    }

    class(result) <- c("IkNN", "list")

    return(result)
}

#' Summarize IkNN Graph Object
#'
#' @description
#' Prints a summary of an IkNN graph object, including the number of vertices,
#' edges, and pruning statistics.
#'
#' @param object An object of class "IkNN", typically output from create.single.iknn.graph()
#' @param ... Additional arguments passed to summary().
#'
#' @return Invisibly returns NULL while printing summary information to the console
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' result <- create.single.iknn.graph(X, k = 3)
#' summary(result)
#'
#' @method summary IkNN
#' @export
summary.IkNN <- function(object, ...) {
    cat("Graph Summary:\n")
    cat("Number of vertices:", length(object$pruned_adj_list), "\n")
    cat("Number of edges:", object$n_edges, "\n")                  # Total number of edges in the original graph
    cat("Number of edges after pruning:", object$n_pruned_edges, "\n")    # Number of edges after pruning
    cat("Number of removed edges:", object$n_removed_edges, "\n")  # Number of edges removed during pruning
    cat("Proportion of edges removed:", object$edge_reduction_ratio,"\n") # Proportion of edges removed (n_removed_edges / n_edges)
}


