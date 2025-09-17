#' Create intersection k-nearest neighbor graphs with dual pruning
#'
#' @description
#' For each \eqn{k \in [k_{\mathrm{min}},\,k_{\mathrm{max}}]}, builds an
#' intersection-weighted k-NN graph and applies two pruning schemes:
#' (1) geometric (path-to-edge ratio) and (2) intersection-size.
#' Optionally performs PCA before graph construction.
#'
#' @param X A numeric matrix (or object coercible to a numeric matrix) with rows
#'     = observations and columns = features.
#' @param kmin Integer \eqn{\ge 1}, the minimum k.
#' @param kmax Integer \eqn{> k_{\mathrm{min}}}, the maximum k.
#' @param max.path.edge.ratio.deviation.thld Numeric in \eqn{[0, 0.2)}.
#'     Geometric pruning removes an edge \eqn{(i,j)} when there exists an
#'     alternative path between \eqn{i} and \eqn{j} whose path/edge length ratio
#'     minus 1.0 is \emph{less than} this threshold. This is a deviation
#'     threshold \eqn{\delta} in \eqn{[0, 0.2)}. Internally we compare the
#'     path-to-edge ratio R to \eqn{1 + \delta}.
#' @param path.edge.ratio.percentile Numeric in \eqn{[0,1]}. Only edges with
#'     length above this percentile are considered for geometric pruning.
#' @param compute.full Logical. If `TRUE`, return the pruned graphs; if `FALSE`,
#'     return only edge statistics.
#' @param pca.dim Positive integer or `NULL`. If not `NULL` and `ncol(X) >
#'     pca.dim`, PCA is used to reduce to at most `pca.dim` components.
#' @param variance.explained Numeric in \eqn{(0,1]} or `NULL`. If not `NULL`,
#'     choose the smallest number of PCs whose cumulative variance explained
#'     exceeds this threshold, capped by `pca.dim`.
#' @param n.cores Integer or `NULL`. Number of CPU cores. `NULL` uses the
#'     maximum available (OpenMP build only).
#' @param verbose Logical; print progress and timing.
#'
#' @return A list of class `"iknn_graphs"` with entries:
#' \describe{
#'   \item{k_statistics}{Matrix of per-\eqn{k} edge counts and reductions.
#'     (If the C++ side supplies column names, they’re preserved.
#'     Otherwise we add names consistent with what the C++ returns.)}
#'   \item{geom_pruned_graphs}{If `compute.full=TRUE`, list of geometrically
#'     pruned graphs (adjacency + weights); otherwise `NULL`.}
#'   \item{isize_pruned_graphs}{If `compute.full=TRUE`, list of intersection-size
#'     pruned graphs; otherwise `NULL`.}
#'   \item{edge_pruning_stats}{List (per \eqn{k}) of matrices with edge-level
#'     statistics (lengths, path/edge ratios, etc.).}
#' }
#'
#' @details
#' Geometric pruning uses the deviation threshold
#' `max.path.edge.ratio.deviation.thld` and the filtering percentile
#' `path.edge.ratio.percentile`. Intersection-size pruning currently uses
#' a maximum alternative path length of 2.
#'
#' @examples
#' # Generate sample data
#' X <- matrix(rnorm(100 * 5), 100, 5)
#'
#' # Basic usage
#' res1 <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10, n.cores = 1,
#'   compute.full = FALSE
#' )
#'
#' # With custom pruning parameters
#' res2 <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   max.path.edge.ratio.deviation.thld = 0.1,
#'   path.edge.ratio.percentile = 0.5,
#'   compute.full = TRUE,
#'   n.cores = 1,
#'   verbose = TRUE
#' )
#'
#' # View statistics for each k
#' print(res2$k_statistics)
#'
#' @export
create.iknn.graphs <- function(X,
                               kmin,
                               kmax,
                               max.path.edge.ratio.deviation.thld = 0.1,
                               path.edge.ratio.percentile = 0.5,
                               compute.full = TRUE,
                               pca.dim = 100,
                               variance.explained = 0.99,
                               n.cores = NULL,
                               verbose = FALSE) {

    ## Coerce & basic checks
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error"))
            stop("X must be a matrix or coercible to a numeric matrix.")
    }
    if (!is.numeric(X)) stop("X must be numeric.")
    if (any(!is.finite(X))) stop("X cannot contain NA/NaN/Inf.")

    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    n <- nrow(X)
    if (n < 2) stop("X must have at least 2 rows (observations).")

    if (!is.numeric(kmin) || length(kmin) != 1 || kmin < 1 || kmin != floor(kmin))
        stop("kmin must be a positive integer.")
    if (!is.numeric(kmax) || length(kmax) != 1 || kmax < kmin || kmax != floor(kmax))
        stop("kmax must be an integer not smaller than kmin.")
    if (n <= kmax)
        stop("Number of observations (nrow(X)) must be greater than kmax.")

    if (!is.numeric(max.path.edge.ratio.deviation.thld) || length(max.path.edge.ratio.deviation.thld) != 1)
        stop("max.path.edge.ratio.deviation.thld must be numeric.")
    if (max.path.edge.ratio.deviation.thld < 0 || max.path.edge.ratio.deviation.thld >= 0.2)
        stop("max.path.edge.ratio.deviation.thld must be in [0, 0.2).")

    if (!is.numeric(path.edge.ratio.percentile) || length(path.edge.ratio.percentile) != 1 ||
        path.edge.ratio.percentile < 0 || path.edge.ratio.percentile > 1)
        stop("path.edge.ratio.percentile must be in [0, 1].")

    if (!is.logical(compute.full) || length(compute.full) != 1)
        stop("compute.full must be TRUE/FALSE.")
    if (!is.logical(verbose) || length(verbose) != 1)
        stop("verbose must be TRUE/FALSE.")

    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim))
            stop("pca.dim must be a positive integer or NULL.")
    }
    if (!is.null(variance.explained)) {
        if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
            variance.explained <= 0 || variance.explained > 1)
            stop("variance.explained must be in (0, 1], or NULL.")
    }

    ## PCA (optional)
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose) message("High-dimensional data detected. Performing PCA.")
        original_dim <- ncol(X)
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(
                X, variance.threshold = variance.explained, max.components = pca.dim
            )
            n_components <- pca_analysis$n.components
            if (verbose) {
                message(sprintf("Using %d PCs (explains %.2f%% variance)",
                                n_components, 100 * pca_analysis$variance.explained))
            }
            X <- pca.project(X, pca_analysis$pca.result, n_components)
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance
            )
        } else {
            if (verbose) message(sprintf("Projecting to first %d PCs", pca.dim))
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained
            )
        }
    }

    ## Note on k: ANN returns self in its kNN sets, this is why k+1 is passed here (for kmin and kmax).
    ## We need to ensure the C++ labels/columns reflect the *original* k.
    result <- .Call(S_create_iknn_graphs,
                    X,
                    as.integer(kmin + 1L),
                    as.integer(kmax + 1L),
                    as.double(max.path.edge.ratio.deviation.thld + 1.0),
                    as.double(path.edge.ratio.percentile),
                    as.logical(compute.full),
                    if (is.null(n.cores)) NULL else as.integer(n.cores),
                    as.logical(verbose),
                    PACKAGE = "gflow")

    ## Normalize optional matrices (placeholders may be empty)
    if (!is.null(result$birth_death_matrix)) {
        if (!is.matrix(result$birth_death_matrix) || nrow(result$birth_death_matrix) == 0) {
            result$birth_death_matrix <- matrix(
                numeric(0), nrow = 0, ncol = 4,
                dimnames = list(NULL, c("start","end","birth_time","death_time"))
            )
        } else if (is.null(colnames(result$birth_death_matrix))) {
            colnames(result$birth_death_matrix) <- c("start","end","birth_time","death_time")
        }
    }
    if (!is.null(result$double_birth_death_matrix)) {
        if (!is.matrix(result$double_birth_death_matrix) || nrow(result$double_birth_death_matrix) == 0) {
            result$double_birth_death_matrix <- matrix(
                numeric(0), nrow = 0, ncol = 4,
                dimnames = list(NULL, c("start","end","birth_time","death_time"))
            )
        }
    }

    ## Add column names to k_statistics only if missing, and match column count
    if (!is.null(result$k_statistics) && is.matrix(result$k_statistics) &&
        is.null(colnames(result$k_statistics))) {
        nc <- ncol(result$k_statistics)
                                        # Common layouts: with k (8 cols) or without k (7 cols)
        if (nc == 8L) {
            colnames(result$k_statistics) <- c("k",
                                               "n_edges",
                                               "n_edges_in_geom_pruned_graph",
                                               "n_geom_removed_edges",
                                               "geom_edge_reduction_ratio",
                                               "n_edges_in_isize_pruned_graph",
                                               "n_isize_removed_edges",
                                               "isize_edge_reduction_ratio"
                                               )
        } else if (nc == 7L) {
            colnames(result$k_statistics) <- c(
                "n_edges",
                "n_edges_in_geom_pruned_graph",
                "n_geom_removed_edges",
                "geom_edge_reduction_ratio",
                "n_edges_in_isize_pruned_graph",
                "n_isize_removed_edges",
                "isize_edge_reduction_ratio"
            )
        }
    }

    attr(result, "kmin") <- kmin
    attr(result, "kmax") <- kmax
    attr(result, "max_path_edge_ratio_deviation_thld") <- max.path.edge.ratio.deviation.thld
    attr(result, "path_edge_ratio_percentile") <- path.edge.ratio.percentile
    if (!is.null(pca_info)) attr(result, "pca") <- pca_info
    class(result) <- "iknn_graphs"
    result
}


#' Summarize an iknn_graphs Object
#'
#' @description
#' Provides a detailed summary of an iknn_graphs object created by the create.iknn.graphs() function.
#' The summary includes statistics for each intersection kNN graph in the sequence, displaying information
#' about the connectivity and structure of the graphs for different k values.
#'
#' @param object An object of class 'iknn_graphs', typically the output of create.iknn.graphs().
#' @param use_isize_pruned Logical. If TRUE, computes and displays statistics for the intersection-size
#'        pruned graphs (isize_pruned_graphs). If FALSE (default), computes statistics for the
#'        geometrically pruned graphs (geom_pruned_graphs).
#' @param ... Additional arguments passed to or from other methods (not currently used).
#'
#' @return Invisibly returns a data frame containing statistics for each graph. The data frame has the following columns:
#'   \item{idx}{The index of the given k value}
#'   \item{k}{The k value for the intersection kNN graph}
#'   \item{n_ccomp}{Number of connected components of the graph}
#'   \item{edges}{Number of edges in the graph}
#'   \item{mean_degree}{Average number of connections per vertex}
#'   \item{min_degree}{Minimum vertex degree (least connected vertex)}
#'   \item{max_degree}{Maximum vertex degree (most connected vertex)}
#'   \item{sparsity}{Graph sparsity, calculated as 1 - density. It measures how many (proportion) potential connections are missing.}
#'
#' @details
#' The summary function extracts and presents key statistics about the structure of each
#' intersection kNN graph in the iknn_graphs object. All graphs share the same number of vertices
#' (equal to the number of rows in the input data matrix), but they differ in the number of edges
#' due to varying k values and the pruning methods applied during graph creation.
#'
#' The function displays general information about the graph sequence, including the number of vertices,
#' the range of k values, and the pruning parameters used. It then presents a tabular summary of statistics
#' for each individual graph, showing how the graph structure changes as k increases.
#'
#' For each intersection kNN graph, the following metrics are calculated:
#' - Number of edges: Total number of connections in the graph
#' - Mean degree: Average number of connections per vertex
#' - Min/Max degree: Range of vertex connectivity
#' - Density: Proportion of potential connections that are actually present
#' - Sparsity: 1 - density, measuring the proportion of missing connections
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' x <- matrix(rnorm(1000), ncol = 5)
#'
#' # Generate intersection kNN graphs
#' iknn.res <- create.iknn.graphs(x, kmin = 3, kmax = 10, n.cores = 1)
#'
#' # Summarize the geometrically pruned graphs
#' summary(iknn.res)
#'
#' # Summarize the intersection-size pruned graphs
#' summary(iknn.res, use_isize_pruned = TRUE)
#'
#' @seealso \code{\link{create.iknn.graphs}} for creating intersection kNN graphs.
#'
#' @export
summary.iknn_graphs <- function(object,
                                use_isize_pruned = FALSE,
                                ...) {

    ## Check if the object is of the correct class
    if (!inherits(object, "iknn_graphs")) {
        stop("Object must be of class 'iknn_graphs'")
    }

    ## Determine which graphs to use
    if (use_isize_pruned) {
        graphs_to_use <- object$isize_pruned_graphs
        graph_type <- "intersection-size pruned"
        if (is.null(graphs_to_use)) {
            stop("Intersection-size pruned graphs not available. Set compute.full = TRUE when creating the graphs.")
        }
    } else {
        graphs_to_use <- object$geom_pruned_graphs
        graph_type <- "geometrically pruned"
        if (is.null(graphs_to_use)) {
            stop("Geometrically pruned graphs not available. Set compute.full = TRUE when creating the graphs.")
        }
    }

    ## Extract relevant information
    kmin <- attr(object, "kmin")
    kmax <- attr(object, "kmax")
    max_path_edge_ratio_deviation_thld <- attr(object, "max_path_edge_ratio_deviation_thld")
    path_edge_ratio_percentile <- attr(object, "path_edge_ratio_percentile")

    ## Get number of vertices (same for all graphs)
    n_vertices <- length(graphs_to_use[[1]]$adj_list)

    ## Initialize table of statistics
    k_values <- kmin:kmax
    n_graphs <- length(k_values)
    stats_table <- data.frame(
        idx = seq_along(k_values),
        k = k_values,
        n_ccomp = numeric(n_graphs),
        edges = numeric(n_graphs),
        mean_degree = numeric(n_graphs),
        min_degree = numeric(n_graphs),
        max_degree = numeric(n_graphs),
        sparsity = numeric(n_graphs),
        stringsAsFactors = FALSE
    )

    ## Calculate statistics for each graph
    for (i in 1:n_graphs) {
        graph <- graphs_to_use[[i]]
        adj_list <- graph$adj_list
        weight_list <- graph$weight_list

        ## Calculate number of edges (sum of adjacency list lengths divided by 2 because each edge is counted twice)
        edge_count <- sum(sapply(adj_list, length)) / 2

        ## Calculate degrees for each vertex
        degrees <- sapply(adj_list, length)

        ## Calculate mean degree
        mean_deg <- mean(degrees)

        ## Calculate min and max degree
        min_deg <- min(degrees)
        max_deg <- max(degrees)

        ## Calculate graph density (ratio of actual edges to potential edges)
        density <- edge_count / (n_vertices * (n_vertices - 1) / 2)

        ## Sparsity measures how many potential connections are missing
        sparsity <- 1 - density

        ## Number of connected components
        n.ccomp <- length(table(graph.connected.components(adj_list)))

        ## Store statistics
        stats_table$edges[i] <- edge_count
        stats_table$mean_degree[i] <- mean_deg
        stats_table$min_degree[i] <- min_deg
        stats_table$max_degree[i] <- max_deg
        stats_table$sparsity[i] <- sparsity
        stats_table$n_ccomp[i] <- n.ccomp
    }

    ## Print summary
    cat("Summary of iknn_graphs object\n")
    cat("----------------------------\n")
    cat("Number of vertices:", n_vertices, "\n")
    cat("k range:", kmin, "to", kmax, "\n")
    cat("Max path-edge ratio deviation threshold:", max_path_edge_ratio_deviation_thld, "\n")
    cat("Path-edge ratio percentile:", path_edge_ratio_percentile, "\n")
    cat("Graph type:", graph_type, "\n\n")

    ## Round numeric columns for cleaner display
    stats_table$mean_degree <- round(stats_table$mean_degree, 2)
    stats_table$sparsity <- round(stats_table$sparsity, 5)

    ## Print table
    print(stats_table, row.names = FALSE)

    ## Return the statistics table invisibly
    invisible(stats_table)
}


# Updated helper function to compute stability metrics
compute.stability.metrics <- function(graphs, k.values, graph_type = "geom") {

    # Initialize results
    n.edges <- numeric(length(k.values))
    n.edges.in.pruned.graph <- numeric(length(k.values))
    edge.reduction.ratio <- numeric(length(k.values))
    js.div <- numeric(length(k.values) - 1)

    # Get the appropriate column name based on graph type
    if (graph_type == "geom") {
        edge_col <- "n_edges_in_geom_pruned_graph"
        ratio_col <- "geom_edge_reduction_ratio"
    } else {
        edge_col <- "n_edges_in_isize_pruned_graph"
        ratio_col <- "isize_edge_reduction_ratio"
    }

    # Compute basic metrics
    for (i in seq_along(k.values)) {
        graph <- graphs[[i]]

        # Count edges in the graph
        if (!is.null(graph$adj_list)) {
            n.edges.in.pruned.graph[i] <- sum(sapply(graph$adj_list, length)) / 2
        }

        # Compute JS divergence for consecutive graphs
        if (i < length(k.values)) {
            next_graph <- graphs[[i + 1]]
            js.div[i] <- compute.degrees.js.divergence(graph, next_graph)
        }
    }

    ## Compute edit distances
    edit.distances <- compute.edit.distances(graphs, k.values)
    edit.distances.lmin <- internal.find.local.minima(edit.distances, k.values)

    ## Fit piecewise linear models
    edit.distances.model <- fit.pwlm(k.values[-length(k.values)], edit.distances)
    edge.model <- fit.pwlm(k.values, n.edges.in.pruned.graph)
    js.model <- fit.pwlm(k.values[-length(k.values)], js.div)

    ## Find local minima
    edit.distances.lmin <- internal.find.local.minima(edit.distances, k.values)
    edge.lmin <- internal.find.local.minima(n.edges.in.pruned.graph, k.values)
    js.lmin <- internal.find.local.minima(js.div, k.values)

    list(
        edit.distances = edit.distances,
        edit.distances.lmin = edit.distances.lmin,
        edit.distances.pwlm = edit.distances.model$model,
        edit.distances.breakpoint = edit.distances.model$breakpoint,
        n.edges = n.edges,
        n.edges.in.pruned.graph = n.edges.in.pruned.graph,
        edge.reduction.ratio = edge.reduction.ratio,
        n.edges.lmin = edge.lmin,
        edge.pwlm = edge.model$model,
        edge.breakpoint = edge.model$breakpoint,
        js.div = js.div,
        js.lmin = js.lmin,
        js.div.pwlm = js.model$model,
        js.div.breakpoint = js.model$breakpoint
    )
}


# Updated helper function to compute JS divergence between degree distributions
compute.degrees.js.divergence <- function(g1, g2) {
    # Handle both old and new graph structures
    if (!is.null(g1$pruned_adj_list)) {
        g1.degrees <- lengths(g1$pruned_adj_list)
    } else if (!is.null(g1$adj_list)) {
        g1.degrees <- lengths(g1$adj_list)
    } else {
        stop("No adjacency list found in g1")
    }

    if (!is.null(g2$pruned_adj_list)) {
        g2.degrees <- lengths(g2$pruned_adj_list)
    } else if (!is.null(g2$adj_list)) {
        g2.degrees <- lengths(g2$adj_list)
    } else {
        stop("No adjacency list found in g2")
    }

    g1.rel.degrees <- g1.degrees / max(g1.degrees)
    g2.rel.degrees <- g2.degrees / max(g2.degrees)

    jensen.shannon.divergence(g1.rel.degrees, g2.rel.degrees)
}


# Helper function to compute edit distances between graphs
compute.edit.distances <- function(graphs, k.values) {
    n <- length(graphs)
    if (n < 2) return(numeric(0))

    edit.distances <- numeric(n - 1)

    for (i in 1:(n-1)) {
        g1 <- graphs[[i]]
        g2 <- graphs[[i+1]]

        # Get adjacency lists
        adj1 <- g1$adj_list
        adj2 <- g2$adj_list

        # Count edge differences
        edges1 <- set()
        edges2 <- set()

        for (v in seq_along(adj1)) {
            for (neighbor in adj1[[v]]) {
                if (v < neighbor) {  # Avoid counting edges twice
                    edges1 <- edges1 + set(paste(v, neighbor, sep="-"))
                }
            }
        }

        for (v in seq_along(adj2)) {
            for (neighbor in adj2[[v]]) {
                if (v < neighbor) {  # Avoid counting edges twice
                    edges2 <- edges2 + set(paste(v, neighbor, sep="-"))
                }
            }
        }

        # Edit distance is the size of symmetric difference
        edit.distances[i] <- length(setdiff(edges1, edges2)) + length(setdiff(edges2, edges1))
    }

    return(edit.distances)
}


# Define set operations if not available
if (!exists("set")) {
    set <- function(...) {
        unique(c(...))
    }
}


#' Find Optimal k Parameter Using Edge Persistence Analysis
#'
#' @description
#' Analyzes the stability of graph structures across different k values by examining
#' the persistence patterns of edges. The function combines multiple stability metrics
#' to identify the k value that produces the most stable graph structure.
#'
#' @param birth.death.matrix A numeric matrix with columns "birth_time" and "death_time"
#'        containing the birth and death times of edges. Each row represents one edge.
#'        For edges that persist through kmax, death_time equals kmax + 1.
#' @param kmin The minimum k value to consider
#' @param kmax The maximum k value to consider
#' @param matrix_type Character string describing the type of matrix (default: "geom").
#'        Used for informative warning messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{stability.scores}{Numeric vector of stability scores for each k value}
#'   \item{k.values}{Vector of k values analyzed (kmin:kmax)}
#'   \item{opt.k}{The k value with highest stability score}
#' }
#'
#' @details
#' The function calculates stability scores for each k value using three components:
#' 1. Average persistence (lifetime) of edges present at k
#' 2. Proportion of persistent edges (those that continue to kmax+1)
#' 3. Edge stability measure that penalizes recently born or soon-to-die edges
#'
#' The final stability score for each k is the product of these three components.
#' Higher scores indicate more stable graph structures.
#'
#' @examples
#' \dontrun{
#' # Assuming we have birth-death data from create.iknn.graphs()
#' result <- create.iknn.graphs(X, kmin = 3, kmax = 10)
#'
#' # Find optimal k
#' stability <- find.optimal.k(result$birth_death_matrix, kmin = 3, kmax = 10)
#'
#' # Plot stability scores
#' plot(stability$k.values, stability$stability.scores, type = "l",
#'      xlab = "k", ylab = "Stability Score")
#' abline(v = stability$opt.k, col = "red", lty = 2)
#' }
#'
#' @seealso
#' \code{\link{create.iknn.graphs}} for generating the birth-death matrix
#'
#' @export
find.optimal.k <- function(birth.death.matrix, kmin, kmax, matrix_type = "geom") {
    # Handle empty birth/death matrix
    if (is.null(birth.death.matrix) || nrow(birth.death.matrix) == 0) {
        warning(paste("Empty", matrix_type, "birth/death matrix. Returning middle k value."))
        return(list(
            stability.scores = rep(0, kmax - kmin + 1),
            k.values = kmin:kmax,
            opt.k = floor((kmin + kmax) / 2)
        ))
    }

    # Calculate persistence of each edge
    persistence <- birth.death.matrix[,"death_time"] - birth.death.matrix[,"birth_time"]

    # For each k, calculate stability metrics
    stability.scores <- numeric(kmax - kmin + 1)

    for(k in kmin:kmax) {
        # Find edges that exist at k (birth_time <= k < death_time)
        edges.at.k <- birth.death.matrix[,"birth_time"] <= k &
                     birth.death.matrix[,"death_time"] > k

        if (sum(edges.at.k) > 0) {
            # Calculate stability score combining multiple factors:
            # 1. Average persistence of edges present at k
            avg.persistence <- mean(persistence[edges.at.k])

            # 2. Proportion of persistent edges (those that continue to kmax+1)
            persistent.ratio <- mean(birth.death.matrix[edges.at.k, "death_time"] == (kmax + 1))

            # 3. Edge stability measure that penalizes recently born or soon-to-die edges
            edge.stability <- mean(pmin(k - birth.death.matrix[edges.at.k, "birth_time"],
                                        birth.death.matrix[edges.at.k, "death_time"] - k))

            stability.scores[k - kmin + 1] <- avg.persistence * persistent.ratio * edge.stability
        }
    }

    # Return k with highest stability score
    opt.k <- kmin - 1 + which.max(stability.scores)

    list(stability.scores = stability.scores,
         k.values = kmin:kmax,
         opt.k = opt.k)
}

#' Plot Diagnostics for Intersection k-NN Graph Analysis
#'
#' @description
#' Creates diagnostic plots for analyzing the properties of intersection k-NN graphs
#' across different k values. Can display edit distances, edge counts, and
#' Jensen-Shannon divergence metrics.
#'
#' @param x A list containing analysis results for intersection k-NN graphs, typically
#'        with components like 'k.values', 'edit.distances', 'n.edges.in.pruned.graph',
#'        and 'js.div'.
#' @param type Character string specifying the type of plot. Either "diag" for diagnostic
#'        plots (default) or "graph" for network visualization.
#' @param diags Character vector specifying which diagnostics to show. Options include
#'        "edist" (edit distances), "edge" (edge counts), and "deg" (degree distribution).
#'        Default is c("edist", "edge", "deg").
#' @param with.pwlm Logical. If TRUE, overlays piecewise linear model fits on the plots.
#'        Default is TRUE.
#' @param with.lmin Logical. If TRUE, shows vertical lines at local minimum points.
#'        Default is FALSE.
#' @param breakpoint.col Color for breakpoint vertical lines. Default is "blue".
#' @param lmin.col Color for local minima vertical lines. Default is "gray".
#' @param k Optional numeric value to highlight a specific k value on the plots.
#' @param mar Numeric vector of plot margins (bottom, left, top, right). Default is
#'        c(2.5, 2.5, 0.5, 0.5).
#' @param mgp Numeric vector for axis title, labels, and line positions.
#'        Default is c(2.5, 0.5, 0).
#' @param tcl Numeric value for tick mark length. Default is -0.3.
#' @param xline Numeric value for position of x-axis label. Default is 2.4.
#' @param yline Numeric value for position of y-axis label. Default is 3.15.
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating plots.
#'
#' @examples
#' \dontrun{
#' # Create sample data for IkNNgraphs analysis
#' set.seed(123)
#' res <- list(
#'   k.values = 3:10,
#'   edit.distances = c(25, 18, 14, 12, 10, 9, 8, 7),
#'   n.edges.in.pruned.graph = c(50, 80, 100, 115, 125, 132, 138, 142),
#'   js.div = c(0.5, 0.3, 0.2, 0.15, 0.12, 0.1, 0.09, 0.08),
#'   edit.distances.breakpoint = 5,
#'   n.edges.in.pruned.graph.breakpoint = 6,
#'   js.div.breakpoint = 5.5
#' )
#'
#' # Plot diagnostics with default settings
#'   plot.IkNNgraphs(res)
#'
#'   # Plot only edit distances and edge counts
#'   plot.IkNNgraphs(res, diags = c("edist", "edge"))
#' }
#'
#' @importFrom graphics par mtext abline
#'
#' @export
plot.IkNNgraphs <- function(x,
                            type = "diag",
                            diags = c("edist","edge","deg"),
                            with.pwlm = TRUE,
                            with.lmin = FALSE,
                            breakpoint.col = "blue",
                            lmin.col = "gray",
                            k = NA,
                            mar = c(2.5, 2.5, 0.5, 0.5),
                            mgp = c(2.5,0.5,0),
                            tcl = -0.3,
                            xline = 2.4,
                            yline = 3.15,
                            ...) {

    types <- c("diag")
    type <- match.arg(type, types)

    old.par <- par(no.readonly = TRUE)  # Save old par settings
    on.exit(par(old.par), add = TRUE)   # Restore on exit

    switch(type,
           "diag" = {

               if (!"k.values" %in% names(x)) {
                   stop("k.values not in x")
               }

               if (setequal(diags, c("edist","edge","deg"))) {
                   par(mfrow = c(1,3), mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = x$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(x$k.values, x$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = x$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

                   plot(x$k.values, x$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = x$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = x$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist","edge"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = x$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(x$k.values, x$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = x$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist","deg"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = x$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(x$k.values, x$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = x$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = x$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edge","deg"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = x$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

                   plot(x$k.values, x$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = x$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = x$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = x$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edge"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = x$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = x$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("deg"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(x$k.values, x$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(x$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = x$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = x$js.div.lmin, lty = 2, col = lmin.col)
                   }
               }
           })
}
