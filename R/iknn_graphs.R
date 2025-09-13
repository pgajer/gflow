#' Create Intersection k-Nearest Neighbor Graphs with Dual Pruning Methods
#'
#' @description
#' Computes a sequence of intersection-weighted k-nearest neighbor graphs for k in \eqn{[\text{kmin, \text{kmax}]}}
#' with two different edge pruning methods: geometric pruning and intersection-size pruning.
#' The function can optionally perform dimensionality reduction via PCA before graph construction.
#'
#' @param X numeric matrix where rows represent observations and columns represent features.
#'        Cannot be a data frame.
#' @param kmin integer, minimum number of nearest neighbors (>= 1)
#' @param kmax integer, maximum number of nearest neighbors (> kmin)
#' @param max.path.edge.ratio.deviation.thld numeric, threshold for geometric pruning based on path-to-edge ratio.
#'         If > 0, removes edges where the ratio of alternative path length to direct edge length
#'         minus 1.0 is less than this value. If <= 0, geometric pruning uses a default method.
#'         Must be in the interval [0, 0.2).
#' @param path.edge.ratio.percentile numeric in \eqn{[0,1]}, percentile threshold for edge lengths
#'         considered in geometric pruning. Only edges with length greater than this
#'         percentile are evaluated for path-ratio pruning.
#' @param compute.full logical, if TRUE returns all pruned graphs, if FALSE returns only
#'        edge statistics
#' @param pca.dim Maximum number of principal components to use if dimensionality reduction
#'        is applied (default: 100). Set to NULL to skip dimensionality reduction.
#' @param variance.explained Percentage of variance to be explained by the principal components
#'        (default: 0.99). If this threshold can be met with fewer components than pca.dim,
#'        the smaller number will be used. Set to NULL to use exactly pca.dim components.
#' @param verbose Logical. If TRUE, print progress messages and timing information.
#'        Default is FALSE.
#'
#' @return A list of class "iknn_graphs" containing:
#' \describe{
#'   \item{k_statistics}{Matrix with columns: k, number of edges in original graph,
#'         number of edges after geometric pruning, number of removed edges,
#'         edge reduction ratio, number of edges after intersection-size pruning,
#'         additional edges removed in intersection-size pruning,
#'         intersection-size edge reduction ratio}
#'   \item{geom_pruned_graphs}{If compute_full=TRUE, list of geometrically pruned graphs for each k.
#'         Each graph contains adjacency lists and edge weights.
#'         If compute_full=FALSE, NULL}
#'   \item{isize_pruned_graphs}{If compute_full=TRUE, list of intersection-size pruned graphs for each k.
#'         Each graph contains adjacency lists and edge weights.
#'         If compute_full=FALSE, NULL}
#'   \item{edge_pruning_stats}{List of matrices, one for each k value, containing edge pruning statistics
#'         including edge lengths and path-to-edge length ratios}
#' }
#'
#' @details
#' The function applies two different pruning methods to construct efficient graph representations:
#'
#' 1. Geometric pruning: Based on path-to-edge ratios
#'    - Removes edges where alternative paths exist with similar or better geometric properties
#'    - Controlled by max.path.edge.ratio.deviation.thld and path.edge.ratio.percentile parameters
#'
#' 2. Intersection-size pruning: Based on intersection sizes of k-NN sets
#'    - Removes edges where alternative paths exist with all edges having larger intersection sizes
#'    - Uses a fixed maximum alternative path length of 2
#'
#' Note: The current implementation does not compute edge birth/death times. The birth_death_matrix
#' and double_birth_death_matrix fields may be present but will be empty.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(100 * 5), 100, 5)
#'
#' # Basic usage
#' result <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   compute.full = FALSE
#' )
#'
#' # With custom pruning parameters
#' result <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   max.path.edge.ratio.deviation.thld = 0.1,
#'   path.edge.ratio.percentile = 0.5,
#'   compute.full = TRUE,
#'   verbose = TRUE
#' )
#'
#' # View statistics for each k
#' print(result$k_statistics)
#' }
#'
#' @export
create.iknn.graphs <- function(X,
                               kmin,
                               kmax,
                               ## pruning parameters
                               max.path.edge.ratio.deviation.thld = 0.1,
                               path.edge.ratio.percentile = 0.5,
                               ## other
                               compute.full = TRUE,
                               pca.dim = 100,
                               variance.explained = 0.99,
                               verbose = FALSE) {
    ## Input validation
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

    ## Check kmin
    if (!is.numeric(kmin) || length(kmin) != 1 || kmin < 1 || kmin != floor(kmin)) {
        stop("kmin must be a positive integer")
    }
    ## Check kmax
    if (!is.numeric(kmax) || length(kmax) != 1 || kmax < kmin || kmax != floor(kmax)) {
        stop("kmax must be an integer not smaller than kmin")
    }
    ## Check max.path.edge.ratio.deviation.thld
    if (!is.numeric(max.path.edge.ratio.deviation.thld) || length(max.path.edge.ratio.deviation.thld) != 1) {
        stop("max.path.edge.ratio.deviation.thld must be a numeric value")
    }

    if (max.path.edge.ratio.deviation.thld < 0 || max.path.edge.ratio.deviation.thld >= 0.2) {
        stop("max.path.edge.ratio.deviation.thld must be in the interval [0, 0.2)")
    }

    ## Check path.edge.ratio.percentile
    if (!is.numeric(path.edge.ratio.percentile) || length(path.edge.ratio.percentile) != 1 ||
        path.edge.ratio.percentile < 0.0 || path.edge.ratio.percentile > 1.0) {
        stop("path.edge.ratio.percentile must be a numeric value between 0.0 and 1.0")
    }
    ## Check compute.full
    if (!is.logical(compute.full) || length(compute.full) != 1) {
        stop("compute.full must be a logical value (TRUE/FALSE)")
    }
    ## Check verbose
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a logical value (TRUE/FALSE)")
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
    ## Check if number of observations is sufficient
    if (nrow(X) <= kmax) {
        stop("Number of observations must be greater than kmax")
    }

    ## PCA dimensionality reduction if needed
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose) {
            message("High-dimensional data detected. Performing dimensionality reduction.")
        }

        # Store original dimensions for reporting
        original_dim <- ncol(X)

        # If variance.explained is specified, analyze variance
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(X,
                                                   variance.threshold = variance.explained,
                                                   max.components = pca.dim)

            # Number of components to use (based on variance or pca.dim)
            n_components <- pca_analysis$n.components

            if (verbose) {
                message(sprintf("Using %d principal components (explains %.2f%% of variance)",
                                n_components, pca_analysis$variance.explained * 100))
            }

            # Project data onto selected components
            X <- pca.project(X, pca_analysis$pca.result, n_components)

            # Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance,
                pca_projected_X = X
            )
        } else {
            # Use fixed number of components (pca.dim)
            if (verbose) {
                message(sprintf("Projecting data onto first %d principal components", pca.dim))
            }

            # Perform PCA and projection
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)

            # Calculate variance explained by pca.dim components
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)

            # Store PCA information
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained,
                pca_projected_X = X
            )
        }
    }

    ## Call the C++ function with updated signature
    ## Note: The C++ function expects k values incremented by 1 to account for
    ## ANN library including the reference vertex in its kNN set
    result <- .Call(S_create_iknn_graphs,
                    X,
                    as.integer(kmin + 1),
                    as.integer(kmax + 1),
                    as.double(max.path.edge.ratio.deviation.thld + 1.0),
                    as.double(path.edge.ratio.percentile),
                    as.logical(compute.full),
                    as.logical(verbose),
                    PACKAGE = "gflow")

    ## The C++ function may return placeholder birth_death matrices
    ## We'll keep them for backward compatibility but note they may be empty
    if (!is.null(result$birth_death_matrix)) {
        if (!is.matrix(result$birth_death_matrix) || nrow(result$birth_death_matrix) == 0) {
            ## Create an empty matrix with correct structure
            result$birth_death_matrix <- matrix(numeric(0),
                                                nrow = 0,
                                                ncol = 4,
                                                dimnames = list(NULL,
                                                                c("start", "end",
                                                                  "birth_time", "death_time")))
        } else if (is.null(colnames(result$birth_death_matrix))) {
            ## Set column names if they weren't set in C++
            colnames(result$birth_death_matrix) <- c("start", "end",
                                                     "birth_time", "death_time")
        }
    }

    ## Similar handling for double_birth_death_matrix if it exists
    if (!is.null(result$double_birth_death_matrix)) {
        if (!is.matrix(result$double_birth_death_matrix) || nrow(result$double_birth_death_matrix) == 0) {
            result$double_birth_death_matrix <- matrix(numeric(0),
                                                       nrow = 0,
                                                       ncol = 4,
                                                       dimnames = list(NULL,
                                                                       c("start", "end",
                                                                         "birth_time", "death_time")))
        }
    }

    ## Add names to k_statistics columns if they weren't set in C++
    if (!is.null(result$k_statistics) && is.null(colnames(result$k_statistics))) {
        colnames(result$k_statistics) <- c("k",
                                           "n_edges",
                                           "n_edges_in_geom_pruned_graph",
                                           "n_geom_removed_edges",
                                           "geom_edge_reduction_ratio",
                                           "n_edges_in_isize_pruned_graph",
                                           "n_isize_removed_edges",
                                           "isize_edge_reduction_ratio")
    }

    ## Parameter attributes
    attr(result, "kmin") <- kmin
    attr(result, "kmax") <- kmax
    attr(result, "max_path_edge_ratio_deviation_thld") <- max.path.edge.ratio.deviation.thld
    attr(result, "path_edge_ratio_percentile") <- path.edge.ratio.percentile

    # Add PCA-related attributes if PCA was performed
    if (!is.null(pca_info)) {
        attr(result, "pca") <- pca_info
    }

    class(result) <- "iknn_graphs"

    return(result)
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
#' iknn.res <- create.iknn.graphs(x, kmin = 3, kmax = 10)
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
