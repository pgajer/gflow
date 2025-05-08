#' @useDynLib gflow, .registration = TRUE
NULL

#' Create Intersection k-Nearest Neighbor Graphs with Multi-Stage Pruning
#'
#' @description
#' Computes a sequence of intersection-weighted k-nearest neighbor graphs for k in [kmin, kmax]
#' with multi-stage edge pruning to produce efficient graph representations. The function
#' tracks birth and death times of edges across different k values and can optionally
#' perform dimensionality reduction via PCA before graph construction.
#'
#' @param X numeric matrix where rows represent observations and columns represent features.
#'        Cannot be a data frame.
#' @param kmin integer, minimum number of nearest neighbors (>= 1)
#' @param kmax integer, maximum number of nearest neighbors (> kmin)
#' @param pruning.thld numeric in (0,0.9), controlling intensity of the first-stage geometric edge pruning.
#'         Edge weight relative deviation is computed as rel_dev = (w_ij + w_jk) / w_ik - 1.0.
#'         Edges with rel_dev < pruning.thld are removed.
#' @param outlier.long.edge.thld numeric in [0,1], threshold percentile for second-stage pruning of long edges.
#'         Default is 0.1, which removes the top 10% longest edges while preserving connectivity.
#' @param max.path.edge.ratio.deviation.thld numeric, threshold for third-stage pruning based on path-to-edge ratio.
#'         If > 0, removes edges where the ratio of alternative path length to direct edge length
#'         is less than or equal to this value. If <= 0, third-stage pruning is skipped.
#' @param path.edge.ratio.percentile numeric in [0,1], percentile threshold for edge lengths
#'         considered in third-stage pruning. Only edges with length greater than this
#'         percentile are evaluated for path-ratio pruning.
#' @param compute.full logical, if TRUE returns all pruned graphs, if FALSE returns only
#'        edge statistics and birth/death times
#' @param pca.dim Maximum number of principal components to use if dimensionality reduction
#'        is applied (default: 100). Set to NULL to skip dimensionality reduction.
#' @param variance.explained Percentage of variance to be explained by the principal components
#'        (default: 0.99). If this threshold can be met with fewer components than pca.dim,
#'        the smaller number will be used. Set to NULL to use exactly pca.dim components.
#' @param verbose Logical. If TRUE, print progress messages and timing information.
#'        Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{birth_death_matrix}{Matrix with columns: start vertex (1-based), end vertex (1-based),
#'         birth time, and death time. For each edge, birth time is k when edge first appears,
#'         death time is (last k where edge exists) + 1}
#'   \item{k_statistics}{Matrix with columns: k, number of edges in original graph,
#'         number of edges after first-stage pruning, number of removed edges,
#'         edge reduction ratio, number of edges after second-stage pruning,
#'         additional edges removed in second-stage, second-stage edge reduction ratio}
#'   \item{pruned_graphs}{If compute_full=TRUE, list of first-stage pruned graphs for each k.
#'         Each graph contains adjacency lists and distances.
#'         If compute_full=FALSE, NULL}
#'   \item{double_birth_death_matrix}{Similar to birth_death_matrix but for the multi-stage
#'         pruned graphs}
#'   \item{double_pruned_graphs}{If compute_full=TRUE, list of multi-stage pruned graphs for each k.
#'         If compute_full=FALSE, NULL}
#'   \item{edge_pruning_stats}{List of matrices, one for each k value, containing edge pruning statistics
#'         including edge lengths and path-to-edge length ratios}
#' }
#'
#' @details
#' The function applies a multi-stage pruning process to construct efficient graph representations:
#'
#' 1. First stage: Geometric pruning based on relative deviation threshold
#'    - Removes edges where indirect paths through common neighbors provide a good alternative
#'    - Controlled by pruning.thld parameter
#'
#' 2. Second stage: Long-edge pruning
#'    - Removes edges longer than the specified percentile threshold if they can be replaced
#'      by alternate paths
#'    - Controlled by outlier.long.edge.thld parameter
#'
#' 3. Third stage (optional): Path-ratio pruning
#'    - Removes edges where the ratio of alternative path length to direct edge length
#'      is below the specified threshold
#'    - Controlled by max.path.edge.ratio.thld and path.edge.ratio.percentile parameters
#'
#' For each edge, birth and death times are tracked:
#' * Birth time b(e) is the smallest k where edge e appears in the pruned graph
#' * Death time d(e) is one plus the largest k where edge e appears in the pruned graph
#' * If edge persists through kmax, its death time will be kmax + 1
#' * The interval [b(e), d(e)) contains exactly the k values where edge e exists
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
#' # With additional pruning options
#' result <- create.iknn.graphs(
#'   X, kmin = 3, kmax = 10,
#'   pruning.thld = 0.1,
#'   outlier.long.edge.thld = 0.2,
#'   max.path.edge.ratio.deviation.thld = 0.2,
#'   path.edge.ratio.percentile = 0.5,
#'   compute.full = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Access birth-death statistics
#' head(result$birth.death.matrix)
#'
#' # View statistics for each k
#' print(result$k.statistics)
#' }
#'
#' @export
create.iknn.graphs <- function(X,
                               kmin,
                               kmax,
                               ## pruning parameters
                               pruning.thld = 0.1,
                               outlier.long.edge.thld = 0.1,
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
    ## Check pruning.thld
    if (!is.numeric(pruning.thld) || length(pruning.thld) != 1 || pruning.thld <= 0 || pruning.thld >= 0.9) {
        stop("pruning.thld must be a positive numeric value less than 0.9")
    }
    ## Check outlier.long.edge.thld
    if (!is.numeric(outlier.long.edge.thld) || length(outlier.long.edge.thld) != 1 ||
        outlier.long.edge.thld < 0.0 || outlier.long.edge.thld > 1.0) {
        stop("outlier.long.edge.thld must be a numeric value between 0.0 and 1.0")
    }
    ## Check max.path.edge.ratio.deviation.thld
    if (!is.numeric(max.path.edge.ratio.deviation.thld) || length(max.path.edge.ratio.deviation.thld) != 1) {
        stop("max.path.edge.ratio.thld must be a numeric value")
    }

    if (max.path.edge.ratio.deviation.thld < 0 || max.path.edge.ratio.deviation.thld >= 1) {
        stop("max.path.edge.ratio.deviation.thld must be in the interval [0,1).")
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

    ## Call the C++ function
    result <- .Call("S_create_iknn_graphs",
                    X,
                    as.integer(kmin + 1),
                    as.integer(kmax + 1),
                    as.double(pruning.thld),
                    as.double(outlier.long.edge.thld),
                    as.double(max.path.edge.ratio.deviation.thld + 1.0),
                    as.double(path.edge.ratio.percentile),
                    as.logical(compute.full),
                    as.logical(verbose),
                    PACKAGE = "gflow")

    ## Check if birth_death_matrix is actually a matrix
    if (!is.null(result$birth_death_matrix)) {
        if (!is.matrix(result$birth_death_matrix)) {
            ## If no edges were found, create an empty matrix with correct structure
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

    ## Add names to k_statistics columns if they weren't set in C++
    if (!is.null(result$k_statistics) && is.null(colnames(result$k_statistics))) {
        colnames(result$k_statistics) <- c("k", "n_edges", "n_pruned_edges",
                                           "n_removed_edges", "edge_reduction_ratio",
                                           "n_edges_in_double_pruned_graph",
                                           "n_removed_edges_in_double_pruning",
                                           "double_edge_reduction_ratio")
    }

    ## Parameter attributes
    attr(result, "kmin") <- kmin
    attr(result, "kmax") <- kmax
    attr(result, "pruning_threshold") <- pruning.thld
    attr(result, "long_edge_threshold") <- outlier.long.edge.thld
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
#' @param use_double_pruned Logical. If TRUE, computes and displays statistics for the multi-stage
#'        pruned graphs (double_pruned_graphs). If FALSE (default), computes statistics for the
#'        first-stage pruned graphs (pruned_graphs). The multi-stage pruned graphs incorporate
#'        additional long-edge and path-ratio pruning beyond the basic geometric pruning.
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
#' due to varying k values and the pruning threshold applied during graph creation.
#'
#' The function displays general information about the graph sequence, including the number of vertices,
#' the range of k values, and the pruning threshold used. It then presents a tabular summary of statistics
#' for each individual graph, showing how the graph structure changes as k increases.
#'
#' For each intersection kNN graph, the following metrics are calculated:
#' - Number of edges: Total number of connections in the graph
#' - Mean degree: Average number of connections per vertex
#' - Min/Max degree: Range of vertex connectivity
#' - Density: Proportion of potential connections that are actually present
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' x <- matrix(rnorm(1000), ncol = 5)
#'
#' # Generate intersection kNN graphs
#' iknn.res <- create.iknn.graphs(x, kmin = 3, kmax = 10, pruning.thld = 0.7)
#'
#' # Summarize the graphs
#' summary(iknn.res)
#'
#' # Store the summary statistics for further analysis
#' graph_stats <- summary(iknn.res)
#'
#' @seealso \code{\link{create.iknn.graphs}} for creating intersection kNN graphs.
#'
#' @export
summary.iknn_graphs <- function(object,
                                use_double_pruned = TRUE,
                                ...) {

    ## Check if the object is of the correct class
    if (!inherits(object, "iknn_graphs")) {
        stop("Object must be of class 'iknn_graphs'")
    }

    ## Determine which graphs to use
    if (use_double_pruned) {
        graphs_to_use <- object$double_pruned_graphs
        graph_type <- "double_pruned"
        if (is.null(graphs_to_use)) {
            stop("Double pruned graphs not available. Set compute.full = TRUE when creating the graphs.")
        }
    } else {
        graphs_to_use <- object$pruned_graphs
        graph_type <- "pruned"
        if (is.null(graphs_to_use)) {
            stop("Pruned graphs not available. Set compute.full = TRUE when creating the graphs.")
        }
    }

    ## Extract relevant information
    kmin <- attr(object, "kmin")
    kmax <- attr(object, "kmax")
    pruning_threshold <- attr(object, "pruning_threshold")

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
        ## density = numeric(n_graphs),
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

        ## Sparsity measures how many potential connections are missing, rather than how many are present. It ranges from 0 to 1, but with the opposite interpretation to density:
        ## - A sparsity of 0 means the graph is complete (no missing edges)
        ## - A sparsity of 1 means the graph has no edges
        ## - Values in between indicate the proportion of potential edges that are absent
        sparsity <- 1 - density

        ## Number of connected components
        n.ccomp <- length(table(graph.connected.components(adj_list)))

        ## Store statistics
        stats_table$edges[i] <- edge_count
        stats_table$mean_degree[i] <- mean_deg
        stats_table$min_degree[i] <- min_deg
        stats_table$max_degree[i] <- max_deg
        ## stats_table$density[i] <- density
        stats_table$sparsity[i] <- sparsity
        stats_table$n_ccomp[i] <- n.ccomp
    }

    ## Print summary
    cat("Summary of iknn_graphs object\n")
    cat("----------------------------\n")
    cat("Number of vertices:", n_vertices, "\n")
    cat("k range:", kmin, "to", kmax, "\n")
    cat("Pruning threshold:", pruning_threshold, "\n")
    cat("Graph type:", if(use_double_pruned) "Multi-stage pruned graphs" else "First-stage pruned graphs", "\n\n")

    ## Round numeric columns for cleaner display
    stats_table$mean_degree <- round(stats_table$mean_degree, 2)
    ## stats_table$density <- round(stats_table$density, 5)
    stats_table$sparsity <- round(stats_table$sparsity, 5)

    ## Print table
    print(stats_table, row.names = FALSE)

    ## Return the statistics table invisibly
    invisible(stats_table)
}



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

    result <- .Call("S_create_single_iknn_graph",
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
#' @param x An object of class "IkNN", typically output from create.iknn.graph()
#'
#' @return Invisibly returns NULL while printing summary information to the console
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' result <- create.iknn.graph(X, k = 3)
#' summary(result)
#'
#' @method summary IkNN
#' @export
summary.IkNN <- function(x, ...) {
    cat("Graph Summary:\n")
    cat("Number of vertices:", length(x$pruned_adj_list), "\n")
    cat("Number of edges:", x$n_edges, "\n")                  # Total number of edges in the original graph
    cat("Number of edges after pruning:", x$n_pruned_edges, "\n")    # Number of edges after pruning
    cat("Number of removed edges:", x$n_removed_edges, "\n")  # Number of edges removed during pruning
    cat("Proportion of edges removed:", x$edge_reduction_ratio,"\n") # Proportion of edges removed (n_removed_edges / n_edges)
}



#' Verify Graph Pruning Implementation
#'
#' @description
#' Verifies the equivalence between old and new implementations of graph pruning
#' algorithms by comparing their outputs and identifying any discrepancies.
#'
#' @param X A numeric matrix where rows represent data points and columns represent features.
#'          Data frames will be coerced to matrices.
#' @param k Integer scalar. The number of nearest neighbors to consider.
#' @param max.alt.path.length Integer. Maximum allowed length for alternative paths when
#'        determining if an edge can be pruned.
#' @return
#' A list with the following components:
#' \describe{
#'   \item{identical}{Logical. TRUE if no discrepancies were found between implementations.}
#'   \item{total_discrepancies}{Integer. The total number of vertices with discrepancies.}
#'   \item{discrepancies}{List of length n (number of vertices). NULL for vertices without
#'        discrepancies. For vertices with discrepancies, contains a list with:
#'        \itemize{
#'          \item vertex: Integer index of the vertex
#'          \item missing: Matrix of edges (vertex, weight) present in old but not new implementation
#'          \item extra: Matrix of edges (vertex, weight) present in new but not old implementation
#'        }}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(100 * 5), 100, 5)
#'
#' # Verify pruning with k=5 neighbors and max path length of 3
#' result <- verify.pruning(X, k = 5, max.alt.path.length = 3)
#'
#' # Check if implementations match
#' if (result$identical) {
#'   message("Implementations produce identical results")
#' } else {
#'   message(sprintf("Found %d vertices with discrepancies",
#'                  result$total_discrepancies))
#' }
#' }
#'
#' @seealso
#' Related functions for graph construction and manipulation
#'
#' @export
verify.pruning <- function(X, k, max.alt.path.length) {

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

    if (!is.numeric(k) || length(k) != 1 || k != round(k) || k < 1) {
        stop("k must be a positive integer scalar")
    }

    if (!is.numeric(max.alt.path.length) ||
        length(max.alt.path.length) != 1 ||
        max.alt.path.length != round(max.alt.path.length) ||
        max.alt.path.length < 1) {
        stop("max.alt.path.length must be a positive integer scalar")
    }

    .Call("S_verify_pruning",
          X,
          as.integer(k + 1),
          as.integer(max.alt.path.length),
          PACKAGE = "gflow")
}



# Helper function to compute stability metrics
compute.stability.metrics <- function(graphs, k.values) {

    # Initialize results
    n.edges <- numeric(length(k.values))
    n.edges.in.pruned.graph <- numeric(length(k.values))
    edge.reduction.ratio <- numeric(length(k.values))
    js.div <- numeric(length(k.values) - 1)

    # Compute basic metrics
    for (i in seq_along(k.values)) {
        graph <- graphs[[i]]
        n.edges[i] <- graph$n_edges
        n.edges.in.pruned.graph[i] <- graph$n_edges_in_pruned_graph
        edge.reduction.ratio[i] <- graph$edge_reduction_ratio

        # Compute JS divergence for consecutive graphs
        if (i < length(k.values)) {
            next_graph <- graphs[[i + 1]]
            js.div[i] <- compute.degree.js.divergence(graph, next_graph)
        }
    }

    ## Compute edit distances
    edit.distances <- compute.edit_distancesy(graphs, k.values)
    edit.distances.lmin <- find.local.minima(edit.distances, k.values)

    ## Fit piecewise linear models
    edit.distances.model <- fit.pwlm(k.values[-length(k.values)], edit.distances)
    edge.model <- fit.pwlm(k.values, n.edges.in.pruned.graph)
    js.model <- fit.pwlm(k.values[-length(k.values)], js.div)

    ## Find local minima
    edit.distances.lmin <- find.local.minima(edit.distances, k.values)
    edge.lmin <- find.local.minima(n.edges.in.pruned.graph, k.values)
    js.lmin <- find.local.minima(js.div, k.values)

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

# Helper function to compute JS divergence between degree distributions
compute.degree.rel.freq.js.divergence <- function(g1, g2) {
    # Compute degrees
    g1.degrees <- lengths(g1$pruned_adj_list)
    g2.degrees <- lengths(g2$pruned_adj_list)

    # Compute degree frequencies
    g1.deg.freqs <- table(factor(g1.degrees, levels = 0:max(g1.degrees, g2.degrees)))
    g2.deg.freqs <- table(factor(g2.degrees, levels = 0:max(g1.degrees, g2.degrees)))

    # Convert to proportions
    g1.deg.props <- as.vector(g1.deg.freqs / sum(g1.deg.freqs))
    g2.deg.props <- as.vector(g2.deg.freqs / sum(g2.deg.freqs))

    jensen.shannon.divergence(g1.deg.props, g2.deg.props)
}

compute.degrees.js.divergence <- function(g1, g2) {

    g1.degrees <- lengths(g1$pruned_adj_list)
    g2.degrees <- lengths(g2$pruned_adj_list)

    g1.rel.degrees <- g1.degrees / max(g1.degrees)
    g2.rel.degrees <- g2.degrees / max(g2.degrees)

    jensen.shannon.divergence(g1.rel.degrees, g2.rel.degrees)
}

#' Fit a Piecewise Linear Model
#'
#' @description
#' Fits a piecewise linear model to data using the segmented package. If the segmented
#' regression fails, falls back to simple linear regression.
#'
#' @param x Numeric vector containing the predictor variable values
#' @param y Numeric vector containing the response variable values
#'
#' @return An object of class "pwlm" containing:
#'   \item{model}{Either a segmented model object or a linear model object if segmented regression failed}
#'   \item{breakpoint}{Numeric value indicating the estimated breakpoint. NA if using fallback linear model}
#'   \item{x}{The original x values}
#'   \item{y}{The original y values}
#'   \item{type}{Character string indicating "segmented" or "linear"}
#'
#' @export
fit.pwlm1 <- function(x, y) {
    if (requireNamespace("segmented", quietly = TRUE)) {
        init.lm <- lm(y ~ x)
        pwlm <- try(segmented::segmented(init.lm, seg.Z = ~x, psi = list(x = mean(x))),
                    silent = TRUE)

        if (!inherits(pwlm, "try-error")) {
            result <- list(
                model = pwlm,
                breakpoint = pwlm$psi[, 2],
                x = x,
                y = y,
                type = "segmented"
            )
            class(result) <- "pwlm"
            return(result)
        }
    }

    # Fallback if segmented regression fails
    result <- list(
        model = lm(y ~ x),
        breakpoint = NA,
        x = x,
        y = y,
        type = "linear"
    )
    class(result) <- "pwlm1"
    return(result)
}

#' Plot Method for PWLM Objects
#'
#' @param x An object of class "pwlm"
#' @param main Character string for the plot title (default: "Piecewise Linear Regression")
#' @param xlab Character string for x-axis label (default: "X")
#' @param ylab Character string for y-axis label (default: "Y")
#' @param point_color Color for data points (default: "black")
#' @param line_color Color for fitted lines (default: "blue")
#' @param breakpoint_color Color for breakpoint vertical line (default: "red")
#' @param ... Additional arguments passed to plot
#'
#' @export
plot.pwlm1 <- function(x,
                       main = "Piecewise Linear Regression",
                       xlab = "X", ylab = "Y",
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

    ## If we have a segmented model
    if (x$type == "segmented") {
        # Get predicted values for smooth line
        new_x <- seq(min(x$x), max(x$x), length.out = 200)
        pred_y <- predict(x$model, newdata = data.frame(x = new_x))

        # Add the fitted lines
        lines(new_x, pred_y, col = line_color, lwd = 2)

        # Add vertical line at breakpoint
        abline(v = x$breakpoint, col = breakpoint_color,
               lty = 2, lwd = 1)

        # Add legend
        legend("topleft",
               legend = c("Data", "Fitted Line", "Breakpoint"),
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


# Helper function to find local minima
find.local.minima <- function(x, k.values) {
    if (length(x) < 3) return(numeric(0))

    is.min <- c(FALSE, diff(diff(x) > 0) == -1, FALSE)
    k.values[is.min]
}

# Helper function to wrap single result
wrap.single.result <- function(res) {
    list(
        graphs = list(res),
        edit.distances = numeric(0),
        edit.distances.lmin = numeric(0),
        edit.distances.pwlm = NULL,
        edit.distances.breakpoint = NULL,
        n.edges = res$n_edges,
        n.edges.in.pruned.graph = res$n_edges_in_pruned_graph,
        edge.reduction.ratio = res$edge_reduction_ratio,
        edge.pwlm = NULL,
        edge.breakpoint = NULL,
        js.div = numeric(0),
        js.div.pwlm = NULL,
        js.div.breakpoint = NULL
    )
}

#' Plot Diagnostics for Intersection k-NN Graph Analysis
#'
#' @description
#' Creates diagnostic plots for analyzing the properties of intersection k-NN graphs
#' across different k values. Can display edit distances, edge counts, and
#' Jensen-Shannon divergence metrics.
#'
#' @param res A list containing analysis results for intersection k-NN graphs, typically
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
#' \dontrun{
#'   plot.IkNNgraphs(res)
#'
#'   # Plot only edit distances and edge counts
#'   plot.IkNNgraphs(res, diags = c("edist", "edge"))
#' }
#'
#' @export
#' @importFrom graphics par mtext abline
#' @importFrom rgl open3d
plot.IkNNgraphs <- function(res,
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
    types <- c("diag", "graph")
    type <- match.arg(type, types)

    old.par <- par(no.readonly = TRUE)  # Save old par settings
    on.exit(par(old.par), add = TRUE)   # Restore on exit

    switch(type,
           "diag" = {

               if (!"k.values" %in% names(res)) {
                   stop("k.values not in res")
               }

               if (setequal(diags, c("edist","edge","deg"))) {
                   par(mfrow = c(1,3), mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = res$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(res$k.values, res$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = res$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

                   plot(res$k.values, res$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = res$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = res$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist","edge"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = res$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(res$k.values, res$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = res$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist","deg"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = res$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

                   plot(res$k.values, res$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = res$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = res$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edge","deg"))) {
                   par(mfrow = c(1,2), mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = res$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

                   plot(res$k.values, res$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = res$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = res$js.div.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edist"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Edit Distance", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$edit.distances.pwlm, add = TRUE, col = "red")
                       abline(v = res$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$edit.distances.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("edge"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Number of Edges in Pruned Graphs", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
                       abline(v = res$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
                   }
                   if (with.lmin) {
                       abline(v = res$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
                   }

               } else if (setequal(diags, c("deg"))) {
                   par(mar = mar, mgp = mgp, tcl = tcl)

                   plot(res$k.values, res$js.div, las = 1, type = "b", xlab = "", ylab = "")
                   mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
                   mtext("Jensen-Shannon Divergence", side = 2, line = yline, outer = FALSE)
                   if (with.pwlm) {
                       plot(res$js.div.pwlm, add = TRUE, col = "red")
                       abline(v = res$js.div.breakpoint, lty = 2, col = breakpoint.col)  # Add a vertical line at the breakpoint
                   }
                   if (with.lmin) {
                       abline(v = res$js.div.lmin, lty = 2, col = lmin.col)
                   }
               }
           },
           "graph" = {
               open3d()

           })
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
find.optimal.k <- function(birth.death.matrix, kmin, kmax) {
    ## Calculate persistence of each edge
    persistence <- birth.death.matrix[,"death_time"] - birth.death.matrix[,"birth_time"]

    ## For each k, calculate stability metrics
    stability.scores <- numeric(kmax - kmin + 1)
    for(k in kmin:kmax) {
        ## Find edges that exist at k (birth_time ≤ k < death_time)
        edges.at.k <- birth.death.matrix[,"birth_time"] <= k &
            birth.death.matrix[,"death_time"] > k

        ## Calculate stability score combining multiple factors:
        ## 1. Average persistence of edges present at k
        ## 2. Proportion of persistent edges (those that continue to kmax+1)
        ## 3. Penalty for recently born or soon-to-die edges
        avg.persistence <- mean(persistence[edges.at.k])
        persistent.ratio <- mean(birth.death.matrix[edges.at.k, "death_time"] == (kmax + 1))
        edge.stability <- mean(pmin(k - birth.death.matrix[edges.at.k, "birth_time"],
                                    birth.death.matrix[edges.at.k, "death_time"] - k))

        stability.scores[k - kmin + 1] <- avg.persistence * persistent.ratio * edge.stability
    }

    ##  Return k with highest stability score
    opt.k <- kmin - 1 + which.max(stability.scores)

    list(stability.scores = stability.scores,
         k.values = kmin:kmax,
         opt.k = opt.k)
}

#' Find Optimal k Using Local Stability Analysis
#'
#' @description
#' Analyzes the stability of graph structures by examining edge persistence within
#' local windows around each k value. This approach focuses on identifying k values
#' where the graph structure remains consistent across neighboring k values.
#'
#' @param birth.death.matrix A numeric matrix with columns "birth_time" and "death_time"
#'        containing the birth and death times of edges
#' @param kmin The minimum k value to consider
#' @param kmax The maximum k value to consider
#' @param window.size Integer specifying the size of the window around each k value
#'        Default is 2
#'
#' @return A list containing:
#' \describe{
#'   \item{stability.scores}{Numeric vector of stability scores for each k value}
#'   \item{opt.k}{The k value with highest stability score}
#' }
#'
#' @details
#' For each k value, the function:
#' 1. Defines a window of size [k - window.size, k + window.size]
#' 2. Counts edges that exist throughout this window
#' 3. Normalizes the count by total number of edges
#'
#' Higher scores indicate k values where edge structure remains more stable
#' across neighboring k values.
#'
#' @export
find.optimal.k.local <- function(birth.death.matrix,
                                 kmin,
                                 kmax,
                                 window.size = 2) {
    ## For each k, calculate how many edges remain stable in a window around k
    stability.scores <- numeric(kmax - kmin + 1)

    for(k in kmin:kmax) {
        ## Define window boundaries
        k.low <- max(k - window.size, kmin)
        k.high <- min(k + window.size, kmax)

        ## Find edges that exist throughout the window
        stable.edges <- birth.death.matrix[,"birth_time"] <= k.low &
            birth.death.matrix[,"death_time"] > k.high

        ## Calculate stability score
        stability.scores[k - kmin + 1] <- sum(stable.edges) / nrow(birth.death.matrix)
    }

    opt.k <- kmin - 1 + which.max(stability.scores)

    list(stability.scores = stability.scores,
         opt.k = opt.k)
}

#' Analyze Graph Stability Using Barcode Density
#'
#' @description
#' Creates and analyzes a barcode representation of edge persistence patterns
#' to identify stable graph structures. Uses density patterns and their changes
#' to measure stability.
#'
#' @param birth.death.matrix A numeric matrix with columns "birth_time" and "death_time"
#'        containing the birth and death times of edges
#' @param kmin The minimum k value to consider
#' @param kmax The maximum k value to consider
#'
#' @return A list containing:
#' \describe{
#'   \item{bars}{Binary matrix representing edge existence across k values}
#'   \item{density.scores}{Number of edges present at each k value}
#'   \item{stability.scores}{Stability scores based on density changes}
#'   \item{opt.k}{The k value with highest stability score}
#' }
#'
#' @details
#' The function:
#' 1. Creates a binary matrix where rows represent k values and columns represent edges
#' 2. Calculates edge density at each k value
#' 3. Computes stability scores based on density changes between consecutive k values
#' 4. Higher scores indicate k values where edge density remains more consistent
#'
#' @export
analyze.barcode.density <- function(birth.death.matrix,
                                    kmin,
                                    kmax) {
    ## Create "barcode" representation
    bars <- matrix(0, nrow = kmax - kmin + 1, ncol = nrow(birth.death.matrix))

    ## Fill in bars for each edge
    for(i in 1:nrow(birth.death.matrix)) {
        birth <- birth.death.matrix[i, "birth_time"]
        death <- birth.death.matrix[i, "death_time"]
        ## if(birth <= kmax && death > kmin) {
        start.idx <- max(1, birth - kmin + 1)
        end.idx <- min(kmax - kmin + 1, death - kmin)
        bars[start.idx:end.idx, i] <- 1
        ##}
    }

    ## Analyze density and stability of bars
    density.scores <- rowSums(bars)
    stability.scores <- numeric(length(density.scores))

    ## Calculate stability using density changes
    for(i in 2:(length(density.scores)-1)) {
        stability.scores[i] <- density.scores[i] /
            (1 + abs(density.scores[i+1] - density.scores[i]) +
             abs(density.scores[i] - density.scores[i-1]))
    }

    opt.k <- kmin - 1 + which.max(stability.scores)

    list(bars = bars,
         density.scores = density.scores,
         stability.scores = stability.scores,
         opt.k = opt.k)
}

#' Find Optimal k Using Combined Methods
#'
#' @description
#' Combines multiple stability analysis approaches to identify the optimal k value
#' by evaluating candidates suggested by different methods.
#'
#' @param birth.death.matrix A numeric matrix with columns "birth_time" and "death_time"
#'        containing the birth and death times of edges
#' @param kmin The minimum k value to consider
#' @param kmax The maximum k value to consider
#'
#' @return The optimal k value selected from candidates suggested by different methods
#'
#' @details
#' The function:
#' 1. Gets candidate k values from three different methods:
#'    - Edge persistence analysis
#'    - Local stability analysis
#'    - Barcode density analysis
#' 2. Evaluates each candidate using multiple stability metrics:
#'    - Average edge persistence
#'    - Future stability (remaining lifetime of edges)
#'    - Past stability (time edges have existed)
#'    - Edge count
#' 3. Selects the candidate with best overall performance across metrics
#'
#' @seealso
#' \code{\link{find.optimal.k}},
#' \code{\link{find.optimal.k.local}},
#' \code{\link{analyze.barcode.density}}
#'
#' @export
find.optimal.k.combined <- function(birth.death.matrix,
                                    kmin,
                                    kmax) {
    ## Get recommendations from different methods
    k1 <- find.optimal.k(birth.death.matrix, kmin, kmax)$opt.k
    k2 <- find.optimal.k.local(birth.death.matrix, kmin, kmax)$opt.k
    k3 <- analyze.barcode.density(birth.death.matrix, kmin, kmax)$opt.k

    ## Calculate additional metrics for each suggested k
    evaluate.k <- function(k) {
        edges <- birth.death.matrix[,"birth_time"] <= k &
            birth.death.matrix[,"death_time"] > k

        ## Calculate various stability metrics
        persistence <- mean(birth.death.matrix[edges, "death_time"] -
                            birth.death.matrix[edges, "birth_time"])
        future.stability <- mean(birth.death.matrix[edges, "death_time"] - k)
        past.stability <- mean(k - birth.death.matrix[edges, "birth_time"])
        edge.count <- sum(edges)

        return(c(persistence, future.stability, past.stability, edge.count))
    }

    ## Evaluate each candidate k
    scores <- rbind(
        evaluate.k(k1),
        evaluate.k(k2),
        evaluate.k(k3)
    )

    ## Choose the k with best overall metrics
    best.idx <- which.max(rowMeans(scale(scores)))
    return(c(k1, k2, k3)[best.idx])
}

#' Analyze Graph Stability Across k Values
#'
#' @description
#' Evaluates multiple stability metrics for each k value in the range [kmin, kmax]
#' to assess the stability characteristics of the graph structure at each k, focusing
#' on the temporal persistence patterns of edges.
#'
#' @param birth.death.matrix A numeric matrix with columns "birth_time" and "death_time"
#'        containing the birth and death times of edges
#' @param kmin The minimum k value to consider
#' @param kmax The maximum k value to consider
#'
#' @return A list containing:
#' \describe{
#'   \item{stability.matrix}{Matrix where rows correspond to k values and columns to different
#'         stability metrics:
#'         - persistence: death_time - birth_time (total lifespan of edges)
#'         - future_stability: death_time - k (remaining lifetime after k)
#'         - past_stability: k - birth_time (time existed before k)
#'         - edge_count: number of edges present at k}
#'   \item{stability.scores}{Vector of stability scores combining temporal metrics.
#'         For each k, score = persistence * min(future_stability, past_stability),
#'         emphasizing k values where edges show balanced stability in both directions}
#'   \item{k.values}{Vector of k values analyzed (kmin:kmax)}
#'   \item{opt.k}{The k value with highest stability score}
#' }
#'
#' @details
#' The function assesses graph stability by analyzing three temporal aspects of edges:
#' persistence (total lifespan), future stability (remaining lifetime), and past
#' stability (time existed). The stability score for each k combines these metrics
#' to identify values where edges demonstrate balanced stability - they have existed
#' for a substantial time (high past stability) and continue to exist (high future
#' stability), with the overall effect weighted by persistence to favor long-lived
#' edges. This approach is inspired by concepts from persistent homology and
#' topological data analysis, where the stability of features across different
#' scales helps identify meaningful structures.
#'
#' @examples
#' \dontrun{
#' result <- create.iknn.graphs.with.bd.stats(X, kmin = 3, kmax = 10)
#' stability <- analyze.graph.stability(result$birth_death_matrix, kmin = 3, kmax = 10)
#'
#' # Plot stability metrics
#' par(mfrow = c(2,2))
#' for(metric in c("persistence", "future_stability", "past_stability")) {
#'   plot(stability$k.values, stability$stability.matrix[,metric], type = "l",
#'        xlab = "k", ylab = metric)
#'   abline(v = stability$opt.k, col = "red", lty = 2)
#' }
#' plot(stability$k.values, stability$stability.scores, type = "l",
#'      xlab = "k", ylab = "Combined Stability Score")
#' abline(v = stability$opt.k, col = "red", lty = 2)
#' }
#'
#' @seealso
#' \code{\link{create.iknn.graphs.with.bd.stats}},
#' \code{\link{find.optimal.k}},
#' \code{\link{find.optimal.k.local}},
#' \code{\link{analyze.barcode.density}}
#'
#' @export
analyze.edge.birth.death.graph.stability <- function(birth.death.matrix,
                                                     kmin,
                                                     kmax) {
    # Initialize matrix to store metrics for each k
    k.range <- kmin:kmax
    n.k <- length(k.range)
    stability.matrix <- matrix(0, nrow = n.k, ncol = 4,
                             dimnames = list(NULL,
                                           c("persistence", "future_stability",
                                             "past_stability", "edge_count")))

    # Evaluate metrics for each k
    for(i in seq_len(n.k)) {
        k <- k.range[i]

        # Find edges present at this k
        edges.at.k <- birth.death.matrix[,"birth_time"] <= k &
                     birth.death.matrix[,"death_time"] > k

        if(sum(edges.at.k) > 0) {
            # Calculate stability metrics
            stability.matrix[i, "persistence"] <-
                mean(birth.death.matrix[edges.at.k, "death_time"] -
                     birth.death.matrix[edges.at.k, "birth_time"])

            stability.matrix[i, "future_stability"] <-
                mean(birth.death.matrix[edges.at.k, "death_time"] - k)

            stability.matrix[i, "past_stability"] <-
                mean(k - birth.death.matrix[edges.at.k, "birth_time"])

            stability.matrix[i, "edge_count"] <- sum(edges.at.k)
        }
    }

    stability.scores <- apply(stability.matrix, 1, function(row) {
        persistence <- row["persistence"]
        future <- row["future_stability"]
        past <- row["past_stability"]
        ## We want high values for both future and past stability
        temporal_balance <- min(future, past)
        ## We also want these stable edges to have high persistence
        score <- persistence * temporal_balance
        return(score)
    })

    # Find optimal k
    opt.k <- k.range[which.max(stability.scores)]

    # Return results
    list(stability.matrix = stability.matrix,
         stability.scores = stability.scores,
         k.values = k.range,
         opt.k = opt.k)
}




# Function to calculate rolling standard deviation
roll_sd <- function(x, window) {
    n <- length(x)
    if (n < window) return(rep(NA, n))
    sd_vals <- numeric(n - window + 1)
    for(i in 1:(n - window + 1)) {
        sd_vals[i] <- sd(x[i:(i + window - 1)])
    }
    return(sd_vals)
}

# Function to detect stabilization using moving window statistics
detect.stabilization <- function(k.values,
                                 n.edges.in.pruned.graph,
                                 window_size = 5,
                                 std_threshold = 0.1,
                                 consec_periods = 5) {

    # Calculate moving statistics
    n <- length(k.values)
    mov_mean <- zoo::rollmean(n.edges.in.pruned.graph, window_size)
    mov_sd <- roll_sd(n.edges.in.pruned.graph, window_size)

    mov_stats <- data.frame(
        k = k.values[window_size:n],
        mean = mov_mean,
        sd = mov_sd
    )

    # Calculate the range of values for relative threshold
    total_range <- diff(range(n.edges.in.pruned.graph))
    sd_threshold <- total_range * std_threshold

    # Find where standard deviation stays below threshold
    stable_periods <- which(mov_stats$sd < sd_threshold)

    # Find consecutive periods
    if(length(stable_periods) > 1) {
        runs <- rle(diff(stable_periods) == 1)
        long_runs <- which(runs$lengths >= (consec_periods - 1))

        if(length(long_runs) > 0) {
            # Get the start of the first long run
            start_pos <- sum(runs$lengths[1:long_runs[1]]) - runs$lengths[long_runs[1]] + 1
            k_stable <- mov_stats$k[stable_periods[start_pos]]
            return(list(
                k = k_stable,
                method = "moving_window",
                stats = mov_stats
            ))
        }
    }
    return(NULL)
}


## Function to detect change points using the changepoint package
detect.changepoints <- function(k.values,
                                n.edges.in.pruned.graph) {
    # Detect multiple change points using PELT algorithm
    cp <- changepoint::cpt.mean(n.edges.in.pruned.graph,
                               method = "PELT",
                               penalty = "MBIC")

    # Get change point locations
    cp_locations <- changepoint::cpts(cp)

    # Find the last major change point
    if(length(cp_locations) > 0) {
        k_stable <- k.values[cp_locations[length(cp_locations)]]
        return(list(
            k = k_stable,
            method = "changepoint",
            changepoints = k.values[cp_locations]
        ))
    } else {
        return(NULL)
    }
}

## Function to detect plateau using derivatives
detect.plateau <- function(k.values,
                           n.edges.in.pruned.graph,
                           smoothing_span = 0.25,
                           deriv_threshold = 0.1) {

    # Smooth the data using LOESS
    smooth_fit <- loess(n.edges.in.pruned.graph ~ k.values, span = smoothing_span)
    smooth_values <- predict(smooth_fit)

    # Calculate approximate derivatives
    derivatives <- diff(smooth_values) / diff(k.values)

    # Find where absolute derivative is consistently small
    stable_regions <- which(abs(derivatives) < deriv_threshold)

    if(length(stable_regions) > 0) {
        # Get the first point where derivative becomes small
        k_stable <- k.values[stable_regions[1] + 1]  # +1 because diff reduces length by 1
        return(list(
            k = k_stable,
            method = "plateau",
            derivatives = derivatives
        ))
    } else {
        return(NULL)
    }
}

## Main function to analyze stability using all methods
analyze.n.edges.in.pruned.graph.stability <- function(k.values,
                                                      n.edges.in.pruned.graph) {
    # Store results from all methods
    results <- list()

    # 1. Moving window analysis
    results$moving_window <- detect.stabilization(k.values, n.edges.in.pruned.graph)

    # 2. Change point detection
    results$changepoint <- detect.changepoints(k.values, n.edges.in.pruned.graph)

    # 3. Plateau detection
    results$plateau <- detect.plateau(k.values, n.edges.in.pruned.graph)

    # Create summary of results
    k_values <- sapply(results, function(x) if(!is.null(x)) x$k else NA)
    methods <- names(results)

    summary_df <- data.frame(
        method = methods,
        k_stable = k_values,
        stringsAsFactors = FALSE
    )

    # Add to results
    results$summary <- summary_df

    return(results)
}

#' Plot Graph Stabilization Analysis
#'
#' @description
#' Creates a visualization of graph stabilization analysis, showing the number of edges in a pruned graph
#' across different k values and highlighting detected stabilization points using various methods.
#'
#' @param k.values Numeric vector containing the k values (number of neighbors) used in the analysis.
#' @param n.edges.in.pruned.graph Numeric vector containing the number of edges in the pruned graph
#'        for each corresponding k value.
#' @param results List containing the results of stability analysis methods. Should include components
#'        for each method ('moving_window', 'changepoint', 'plateau') with detected stabilization
#'        points.
#' @param point_color Character string specifying the color for data points. Default is "black".
#' @param method_colors Named character vector specifying colors for each method.
#'        Default is c("moving_window" = "blue", "changepoint" = "red", "plateau" = "green").
#' @param main Character string for plot title. Default is "Graph Stabilization Analysis".
#' @param xlab Character string for x-axis label. Default is "k values".
#' @param ylab Character string for y-axis label. Default is "Number of Edges in Pruned Graph".
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns NULL. The function is called for its side effect of creating a plot.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' k.values <- 3:15
#' n.edges <- c(30, 45, 65, 95, 120, 140, 155, 168, 175, 180, 183, 184, 185)
#'
#' # Create mock results
#' results <- list(
#'   moving_window = list(k = 8),
#'   changepoint = list(k = 10),
#'   plateau = list(k = 11)
#' )
#'
#' # Create the plot
#' plot.stability.analysis(k.values, n.edges, results)
#'
#' @export
#' @importFrom graphics abline legend
#' @importFrom stats median
plot.stability.analysis <- function(k.values,
                                    n.edges.in.pruned.graph,
                                    results) {
    # Create base plot
    plot(k.values, n.edges.in.pruned.graph, type = "p",
         xlab = "k values", ylab = "Number of Edges in Pruned Graph",
         main = "Graph Stabilization Analysis")

    # Add median line
    abline(h = median(n.edges.in.pruned.graph), col = "gray", lty = 2)

    # Add detected stabilization points
    methods_colors <- c("moving_window" = "blue", "changepoint" = "red", "plateau" = "green")

    for(method in names(methods_colors)) {
        if(!is.null(results[[method]])) {
            abline(v = results[[method]]$k,
                   col = methods_colors[method],
                   lty = 2)
        }
    }

    # Add legend
    legend("topright",
           legend = c("Data", "Median", names(methods_colors)),
           col = c("black", "gray", unlist(methods_colors)),
           lty = c(NA, 2, 2, 2),
           pch = c(1, NA, NA, NA))
}

#' Plot a Matrix of Binary Time Series as Line Segments
#'
#' @description
#' Creates a visualization of a binary (0/1) matrix where each column represents a time series
#' and consecutive 1s are shown as line segments. This is particularly useful for visualizing
#' temporal events or periods of activity across multiple samples.
#'
#' @param bars A matrix of 0s and 1s where rows represent time points and columns represent samples
#' @param kmin Numeric value indicating the starting time point (used for x-axis scaling)
#' @param time.range Optional numeric vector of time points. If NULL, defaults to kmin + 0:(nrow(bars) - 1)
#' @param col Character string specifying the color of the line segments. Default is "gray"
#' @param xlab Character string for x-axis label. Default is "Time"
#' @param ylab Character string for y-axis label. Default is "Sample"
#' @param main Character string for plot title. Default is "Bar Plot Matrix"
#'
#' @return A plot is produced on the current graphics device
#'
#' @details
#' The function visualizes each column of the input matrix as a horizontal series of line segments,
#' where consecutive 1s in the matrix are represented by continuous lines. The y-axis represents
#' different samples (columns), while the x-axis represents time points (rows).
#'
#' @examples
#' # Create a sample matrix
#' bars <- matrix(0, nrow = 10, ncol = 5)
#' bars[2:4, 1] <- 1
#' bars[6:8, 2] <- 1
#'
#' # Basic plot
#' plot.bars.matrix(bars, kmin = 1)
#'
#' # Customized plot
#' plot.bars.matrix(bars, kmin = 1,
#'                 col = "blue",
#'                 xlab = "Time Point",
#'                 ylab = "Sample ID",
#'                 main = "Sample Timeline")
#'
plot.bars.matrix <- function(bars,
                             kmin,
                             time.range = NULL,
                             col = "gray",
                             xlab = "Time",
                             ylab = "Sample",
                             main = "Bar Plot Matrix") {
    # Get dimensions
    n.times <- nrow(bars)
    n.samples <- ncol(bars)

    # Set up time range if not provided
    if (is.null(time.range)) {
        time.range <- kmin + 0:(n.times - 1)
    }

    # Create empty plot
    plot(NA, NA,
         xlim = range(time.range),
         ylim = c(0, n.samples),
         xlab = xlab,
         ylab = ylab,
         main = main)

    # Add gray segments for each column
    for (i in 1:n.samples) {
        # Find runs of 1s
        ones <- which(bars[, i] == 1)
        if (length(ones) > 0) {
            # Get start and end points of runs
            runs <- rle(diff(ones) == 1)
            starts <- c(ones[1])
            if (length(runs$lengths) > 0) {
                starts <- c(starts, ones[cumsum(runs$lengths) + 1])
            }
            ends <- c()
            for (j in 1:length(starts)) {
                run.end <- starts[j]
                while (run.end + 1 <= length(ones) && ones[run.end + 1] == ones[run.end] + 1) {
                    run.end <- run.end + 1
                }
                ends <- c(ends, ones[run.end])
            }

            # Plot segments
            segments(time.range[starts], rep(i, length(starts)),
                    time.range[ends], rep(i, length(ends)),
                    col = col, lwd = 1)
        }
    }
}



## ---------------------------------------------------------------------------
##
## Implementation Neighborhood Preservation Goodness of Fit Metric
##
## ---------------------------------------------------------------------------

# Helper function to compute pairwise Euclidean distances
compute.euclidean.distances <- function(X) {
  n <- nrow(X)
  dists <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        dist <- sqrt(sum((X[i,] - X[j,])^2))
        dists[i,j] <- dist
        dists[j,i] <- dist
      }
    }
  }
  return(dists)
}

# Helper function to compute shortest path distances using Dijkstra's algorithm
compute.shortest.paths <- function(adj_list, weight_list) {
  n <- length(adj_list)
  dists <- matrix(Inf, n, n)
  diag(dists) <- 0

  # Initialize direct connections
  for (i in 1:n) {
    neighbors <- adj_list[[i]]
    weights <- weight_list[[i]]
    dists[i, neighbors] <- weights
  }

  # Floyd-Warshall algorithm
  for (k in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (dists[i,j] > dists[i,k] + dists[k,j]) {
          dists[i,j] <- dists[i,k] + dists[k,j]
        }
      }
    }
  }

  return(dists)
}

# Helper function to get points within epsilon distance
get_epsilon_neighbors <- function(dists, i, epsilon) {
  which(dists[i,] <= epsilon)
}

# Main function to compute neighborhood preservation metric
compute.np <- function(X_true, adj_list, weight_list, epsilon) {
  n <- nrow(X_true)

  # Compute true distances
  true_dists <- compute_euclidean_distances(X_true)

  # Compute graph distances
  graph_dists <- compute_shortest_paths(adj_list, weight_list)

  # Compute neighborhood preservation score
  np_score <- 0
  for (i in 1:n) {
    # Get epsilon neighborhoods
    true_neighbors <- get_epsilon_neighbors(true_dists, i, epsilon)
    graph_neighbors <- get_epsilon_neighbors(graph_dists, i, epsilon)

    ## # Remove the point itself from both sets
    ## true_neighbors <- setdiff(true_neighbors, i)
    ## graph_neighbors <- setdiff(graph_neighbors, i)

    ## # Take only k nearest neighbors in graph space if we have more than k
    ## if (length(graph_neighbors) > k) {
    ##   graph_dists_i <- graph_dists[i, graph_neighbors]
    ##   graph_neighbors <- graph_neighbors[order(graph_dists_i)[1:k]]
    ## }

    # Compute intersection size
    intersection_size <- length(intersect(true_neighbors, graph_neighbors))

    # Add to score
    np_score <- np_score + intersection_size
  }

  return(np_score / n)
}

# library(FNN)
# Helper function to get points within epsilon distance using FNN
## get_epsilon_neighbors_fnn <- function(X, point, epsilon) {
##   # Get all neighbors within epsilon radius
##   nn <- get.knnx(X, matrix(point, nrow=1), k=nrow(X), radius=epsilon)
##   # Return indices of points within epsilon (excluding the query point itself)
##   neighbors <- which(nn$nn.dist[1,] <= epsilon)
##   return(neighbors)
## }

compute.np.from.dists <- function(true_dists, graph_dists, epsilon) {y
    n <- nrow(true_dists)
    np_score <- 0
    for (i in 1:n) {
        ## Get epsilon neighborhoods
        ## true_neighbors <- get_epsilon_neighbors_fnn(X_true, X_true[i,], epsilon)
        true_neighbors <- get_epsilon_neighbors(true_dists, i, epsilon)
        graph_neighbors <- get_epsilon_neighbors(graph_dists, i, epsilon)

        ## Compute intersection size
        intersection_size <- length(intersect(true_neighbors, graph_neighbors))

        ## Add to score
        np_score <- np_score + intersection_size
    }

    return(np_score / n)
}
