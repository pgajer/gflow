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
#'
#' @param kmin Integer \eqn{\ge 1}, the minimum k.
#'
#' @param kmax Integer \eqn{> k_{\mathrm{min}}}, the maximum k.
#'
#' @param max.path.edge.ratio.deviation.thld Numeric in \eqn{[0, 0.2)}.
#'     Geometric pruning removes an edge \eqn{(i,j)} when there exists an
#'     alternative path between \eqn{i} and \eqn{j} whose path/edge length ratio
#'     minus 1.0 is \emph{less than} this threshold. This is a deviation
#'     threshold \eqn{\delta} in \eqn{[0, 0.2)}. Internally we compare the
#'     path-to-edge ratio R to \eqn{1 + \delta}.
#'
#' @param path.edge.ratio.percentile Numeric in \eqn{[0,1]}. Only edges with
#'     length above this percentile are considered for geometric pruning.
#'
#' @param threshold.percentile Numeric in \eqn{[0, 0.5]}. Percentile threshold for
#'     quantile-based edge length pruning. Default is 0 (no quantile pruning).
#'     When greater than 0, edges in the top (1 - threshold.percentile) quantile
#'     by length are removed if their removal preserves graph connectivity.
#'     For example, threshold.percentile = 0.9 removes edges in the top 10\% by length.
#'     This pruning is applied after geometric pruning and targets unusually long edges
#'     based on absolute edge lengths rather than path-to-edge ratios.
#'
#' @param compute.full Logical. If `TRUE`, return the pruned graphs; if `FALSE`,
#'     return only edge statistics.
#'
#' @param pca.dim Positive integer or `NULL`. If not `NULL` and `ncol(X) >
#'     pca.dim`, PCA is used to reduce to at most `pca.dim` components.
#'
#' @param variance.explained Numeric in \eqn{(0,1]} or `NULL`. If not `NULL`,
#'     choose the smallest number of PCs whose cumulative variance explained
#'     exceeds this threshold, capped by `pca.dim`.
#'
#' @param n.cores Integer or `NULL`. Number of CPU cores. `NULL` uses the
#'     maximum available (OpenMP build only).
#'
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
                               threshold.percentile = 0,
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

    if (!is.numeric(threshold.percentile) || length(threshold.percentile) != 1)
        stop("threshold.percentile must be numeric.")
    if (threshold.percentile < 0 || threshold.percentile > 0.5)
        stop("threshold.percentile must be in [0, 0.5].")

    if (!is.logical(compute.full) || length(compute.full) != 1) {
        stop("compute.full must be a single logical value")
    }

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
                    as.double(threshold.percentile),
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
        ## Common layouts: with k (8 cols) or without k (7 cols)
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
#' @param use.isize.pruned Logical. If TRUE, computes and displays statistics for the intersection-size
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
#' summary(iknn.res, use.isize.pruned = TRUE)
#'
#' @seealso \code{\link{create.iknn.graphs}} for creating intersection kNN graphs.
#'
#' @export
summary.iknn_graphs <- function(object,
                                use.isize.pruned = FALSE,
                                ...) {

    ## Check if the object is of the correct class
    if (!inherits(object, "iknn_graphs")) {
        stop("Object must be of class 'iknn_graphs'")
    }

    ## Determine which graphs to use
    if (use.isize.pruned) {
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


#' Compute Stability Metrics Across a Sequence of IkNN Graphs
#'
#' @description
#' Computes stability diagnostics across the pruned IkNN graphs produced by
#' \code{\link{create.iknn.graphs}}. Metrics include:
#' \itemize{
#'   \item Edit distances between consecutive graphs (edge-set symmetric difference)
#'   \item Jensen-Shannon divergence between degree profiles of consecutive graphs
#'   \item Edge counts in pruned graphs across k
#'   \item Piecewise linear fits and breakpoints for each diagnostic curve
#'   \item Local minima locations for each diagnostic curve
#' }
#'
#' The function derives \code{k.values} from the \code{"kmin"} and \code{"kmax"}
#' attributes attached to the input \code{iknn_graphs} object.
#'
#' @param graphs An object of class \code{"iknn_graphs"} returned by
#'   \code{\link{create.iknn.graphs}}.
#' @param graph.type Character string, either \code{"geom"} or \code{"isize"},
#'   selecting which pruned graph sequence is analyzed.
#'
#' @return An object of class \code{"iknn_stability_metrics"} (a list) with fields:
#' \describe{
#'   \item{k.values}{Integer vector of k values analyzed.}
#'   \item{graph.type}{Which pruned graph sequence was used (\code{"geom"} or \code{"isize"}).}
#'   \item{edit.distances}{Numeric vector of length \code{length(k.values)-1}.}
#'   \item{js.div}{Numeric vector of length \code{length(k.values)-1}.}
#'   \item{n.edges}{Numeric vector of length \code{length(k.values)} (if available).}
#'   \item{n.edges.in.pruned.graph}{Numeric vector of length \code{length(k.values)}.}
#'   \item{edge.reduction.ratio}{Numeric vector of length \code{length(k.values)} (if available).}
#'   \item{*.pwlm}{Piecewise linear model objects for each curve (if \code{fit.pwlm} exists).}
#'   \item{*.breakpoint}{Estimated breakpoint for each curve (if \code{fit.pwlm} exists).}
#'   \item{*.lmin}{Local minima k values for each curve (if \code{internal.find.local.minima} exists).}
#' }
#'
#' @examples
#' \dontrun{
#' graphs <- create.iknn.graphs(X, kmin = 6, kmax = 15, n.cores = 4)
#' stab <- compute.stability.metrics(graphs, graph.type = "geom")
#' plot(stab)
#' ok <- find.optimal.k(stab)
#' ok$opt.k
#' }
#'
#' @export
compute.stability.metrics <- function(graphs, graph.type = c("geom", "isize")) {

    graph.type <- match.arg(graph.type)

    if (!inherits(graphs, "iknn_graphs")) {
        stop("graphs must be an object of class 'iknn_graphs' returned by create.iknn.graphs().")
    }

    kmin <- attr(graphs, "kmin")
    kmax <- attr(graphs, "kmax")

    if (!is.numeric(kmin) || !is.numeric(kmax) || length(kmin) != 1 || length(kmax) != 1) {
        stop("graphs must have numeric scalar attributes 'kmin' and 'kmax'.")
    }

    k.values <- as.integer(kmin:kmax)

    graphs.list <- if (graph.type == "geom") graphs$geom_pruned_graphs else graphs$isize_pruned_graphs

    if (is.null(graphs.list)) {
        stop("Requested pruned graphs are not available. Recompute create.iknn.graphs(..., compute.full = TRUE).")
    }

    if (length(graphs.list) != length(k.values)) {
        stop("Length mismatch: pruned graph list length does not match kmin:kmax.")
    }

    ## ------------------------------------------------------------------------
    ## Prefer k_statistics for edge counts if present
    ## ------------------------------------------------------------------------

    k.stats <- graphs$k_statistics
    have.k.stats <- !is.null(k.stats) && is.matrix(k.stats) && nrow(k.stats) >= length(k.values)

    n.edges <- rep(NA_real_, length(k.values))
    n.edges.in.pruned.graph <- rep(NA_real_, length(k.values))
    edge.reduction.ratio <- rep(NA_real_, length(k.values))

    if (have.k.stats) {

        ## Determine row mapping: either explicit 'k' column, or assume order
        if (!is.null(colnames(k.stats)) && ("k" %in% colnames(k.stats))) {
            row.idx <- match(k.values, as.integer(k.stats[, "k"]))
        } else {
            ## Assume the matrix rows already correspond to kmin:kmax
            row.idx <- seq_along(k.values)
        }

        ## Choose columns based on graph.type
        if (!is.null(colnames(k.stats))) {
            if (graph.type == "geom") {
                edge.col <- "n_edges_in_geom_pruned_graph"
                ratio.col <- "geom_edge_reduction_ratio"
            } else {
                edge.col <- "n_edges_in_isize_pruned_graph"
                ratio.col <- "isize_edge_reduction_ratio"
            }

            if ("n_edges" %in% colnames(k.stats)) {
                n.edges <- as.numeric(k.stats[row.idx, "n_edges"])
            }
            if (edge.col %in% colnames(k.stats)) {
                n.edges.in.pruned.graph <- as.numeric(k.stats[row.idx, edge.col])
            }
            if (ratio.col %in% colnames(k.stats)) {
                edge.reduction.ratio <- as.numeric(k.stats[row.idx, ratio.col])
            }
        }
    }

    ## ------------------------------------------------------------------------
    ## Fallback: count pruned edges directly if missing
    ## ------------------------------------------------------------------------

    if (any(!is.finite(n.edges.in.pruned.graph))) {
        for (i in seq_along(k.values)) {
            g <- graphs.list[[i]]
            if (!is.null(g$adj_list)) {
                n.edges.in.pruned.graph[i] <- sum(vapply(g$adj_list, length, integer(1))) / 2
            }
        }
    }

    if (any(!is.finite(edge.reduction.ratio)) && any(is.finite(n.edges)) && all(is.finite(n.edges.in.pruned.graph))) {
        edge.reduction.ratio <- (n.edges - n.edges.in.pruned.graph) / n.edges
    }

    ## ------------------------------------------------------------------------
    ## Consecutive-graph metrics: JS divergence and edit distance
    ## ------------------------------------------------------------------------

    js.div <- numeric(max(0L, length(k.values) - 1L))
    if (length(k.values) >= 2L) {
        for (i in seq_len(length(k.values) - 1L)) {
            js.div[i] <- compute.degrees.js.divergence(graphs.list[[i]], graphs.list[[i + 1L]])
        }
    }

    edit.distances <- compute.edit.distances(graphs.list)

    ## ------------------------------------------------------------------------
    ## Local minima + piecewise linear models (if helpers exist)
    ## ------------------------------------------------------------------------

    have.lmin <- exists("internal.find.local.minima", mode = "function")
    have.pwlm <- exists("fit.pwlm", mode = "function")

    edit.distances.lmin <- integer(0)
    edge.lmin <- integer(0)
    js.lmin <- integer(0)

    edit.distances.model <- NULL
    edge.model <- NULL
    js.model <- NULL

    if (have.lmin && length(k.values) >= 2L) {
        edit.distances.lmin <- internal.find.local.minima(edit.distances, k.values[-length(k.values)])
        js.div.lmin <- internal.find.local.minima(js.div, k.values[-length(k.values)])
    }
    if (have.lmin) {
        edge.lmin <- internal.find.local.minima(n.edges.in.pruned.graph, k.values)
    }

    if (have.pwlm && length(k.values) >= 2L) {
        edit.distances.model <- fit.pwlm(k.values[-length(k.values)], edit.distances)
        js.model <- fit.pwlm(k.values[-length(k.values)], js.div)
    }
    if (have.pwlm) {
        edge.model <- fit.pwlm(k.values, n.edges.in.pruned.graph)
    }

    result <- list(
        k.values = k.values,
        k.tr = k.values[-length(k.values)],
        graph.type = graph.type,
        n.edges = n.edges,
        n.edges.in.pruned.graph = n.edges.in.pruned.graph,
        edge.reduction.ratio = edge.reduction.ratio,
        edit.distances = edit.distances,
        js.div = js.div,
        edit.distances.lmin = edit.distances.lmin,
        n.edges.lmin = edge.lmin,
        js.div.lmin = js.div.lmin,
        edit.distances.pwlm = if (!is.null(edit.distances.model)) edit.distances.model$model else NULL,
        edit.distances.breakpoint = if (!is.null(edit.distances.model)) edit.distances.model$breakpoint else NA_real_,
        n.edges.in.pruned.graph.pwlm = if (!is.null(edge.model)) edge.model$model else NULL,
        n.edges.in.pruned.graph.breakpoint = if (!is.null(edge.model)) edge.model$breakpoint else NA_real_,
        js.div.pwlm = if (!is.null(js.model)) js.model$model else NULL,
        js.div.breakpoint = if (!is.null(js.model)) js.model$breakpoint else NA_real_
    )

    class(result) <- c("iknn_stability_metrics", "list")
    return(result)
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



## ============================================================================
## Helper: edit distance between consecutive graphs
## ============================================================================

##' @keywords internal
compute.edit.distances <- function(graphs.list) {

    n.graphs <- length(graphs.list)
    if (n.graphs < 2L) {
        return(numeric(0))
    }

    edit.distances <- numeric(n.graphs - 1L)

    edge.keys <- function(adj.list) {
        ## Build undirected edge keys "i-j" with i < j, using 1-based vertex ids
        keys <- character(0)

        for (i in seq_along(adj.list)) {
            nbrs <- adj.list[[i]]
            if (length(nbrs) == 0) next
            j <- nbrs[nbrs > i]
            if (length(j) > 0) {
                keys <- c(keys, paste(i, j, sep = "-"))
            }
        }

        unique(keys)
    }

    for (i in seq_len(n.graphs - 1L)) {
        g1 <- graphs.list[[i]]
        g2 <- graphs.list[[i + 1L]]

        if (is.null(g1$adj_list) || is.null(g2$adj_list)) {
            stop("Each graph must contain an 'adj_list'.")
        }

        e1 <- edge.keys(g1$adj_list)
        e2 <- edge.keys(g2$adj_list)

        ## Symmetric difference size
        edit.distances[i] <- length(setdiff(e1, e2)) + length(setdiff(e2, e1))
    }

    edit.distances
}


# Define set operations if not available
if (!exists("set")) {
    set <- function(...) {
        unique(c(...))
    }
}

#' Find Optimal k
#'
#' @description
#' Two modes:
#' \itemize{
#'   \item If \code{x} is a birth-death matrix (legacy), uses edge persistence analysis.
#'   \item If \code{x} is an \code{"iknn_stability_metrics"} object, selects k by
#'     combining edit-distance, JS-divergence, and edge-count stability.
#' }
#'
#' @param x Either a birth-death matrix (legacy interface) or an object returned by
#'   \code{\link{compute.stability.metrics}}.
#' @param ... Additional arguments (see details).
#'
#' @return A list with \code{opt.k} and diagnostic vectors. The exact fields depend on mode.
#'
#' @export
find.optimal.k <- function(x, ...) {

    if (inherits(x, "iknn_stability_metrics")) {
        return(find.optimal.k.from.stability(x, ...))
    }

    ## Legacy behavior: treat x as birth.death.matrix
    args <- list(...)
    kmin <- args$kmin
    kmax <- args$kmax
    matrix.type <- args$matrix_type %||% "geom"

    if (is.null(kmin) || is.null(kmax)) {
        stop("For birth-death input, you must supply kmin and kmax (e.g., find.optimal.k(bd, kmin=..., kmax=...)).")
    }

    find.optimal.k.from.birth.death(x, kmin = kmin, kmax = kmax, matrix_type = matrix.type)
}

## ---------------------------------------------------------------------------
## Stability-based optimal k
## ---------------------------------------------------------------------------

##' @keywords internal
find.optimal.k.from.stability <- function(x,
                                         weights = c(edist = 1, js = 1, edges = 1),
                                         k.range = NULL) {

    k.values <- x$k.values
    n <- length(k.values)

    if (n < 2L) {
        return(list(
            k.values = k.values,
            stability.scores = numeric(0),
            opt.k = if (n == 1L) k.values[1] else NA_integer_
        ))
    }

    ## Comparable ks correspond to transitions k -> k+1
    k.comp <- k.values[-n]
    ed <- x$edit.distances
    js <- x$js.div
    ne <- x$n.edges.in.pruned.graph[-n]

    ## Optional restriction of candidate range
    if (!is.null(k.range)) {
        keep <- (k.comp >= k.range[1]) & (k.comp <= k.range[2])
        k.comp <- k.comp[keep]
        ed <- ed[keep]
        js <- js[keep]
        ne <- ne[keep]
    }

    scale01 <- function(v) {
        if (length(v) == 0) return(v)
        r <- range(v, finite = TRUE)
        if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) {
            return(rep(0.5, length(v)))
        }
        (v - r[1]) / (r[2] - r[1])
    }

    ed.bad <- scale01(ed)   ## higher is worse
    js.bad <- scale01(js)   ## higher is worse
    ne.good <- scale01(ne)  ## higher is better

    ## Combine into a single stability score (higher is better)
    w <- weights
    score <- (1 - ed.bad)^w["edist"] * (1 - js.bad)^w["js"] * (ne.good)^w["edges"]

    opt.k <- k.comp[which.max(score)]

    list(
        k.values = k.comp,
        stability.scores = score,
        opt.k = as.integer(opt.k),
        components = list(
            edit.distances = ed,
            js.div = js,
            n.edges.in.pruned.graph = ne
        ),
        weights = w
    )
}

## ---------------------------------------------------------------------------
## Legacy birth-death optimal k (your existing implementation, moved verbatim)
## ---------------------------------------------------------------------------

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
#' @keywords internal
find.optimal.k.from.birth.death <- function(birth.death.matrix, kmin, kmax, matrix_type = "geom") {

    if (is.null(birth.death.matrix) || nrow(birth.death.matrix) == 0) {
        warning(paste("Empty", matrix_type, "birth/death matrix. Returning middle k value."))
        return(list(
            stability.scores = rep(0, kmax - kmin + 1),
            k.values = kmin:kmax,
            opt.k = floor((kmin + kmax) / 2)
        ))
    }

    persistence <- birth.death.matrix[, "death_time"] - birth.death.matrix[, "birth_time"]
    stability.scores <- numeric(kmax - kmin + 1)

    for (k in kmin:kmax) {

        edges.at.k <- birth.death.matrix[, "birth_time"] <= k &
            birth.death.matrix[, "death_time"] > k

        if (sum(edges.at.k) > 0) {

            avg.persistence <- mean(persistence[edges.at.k])
            persistent.ratio <- mean(birth.death.matrix[edges.at.k, "death_time"] == (kmax + 1))

            edge.stability <- mean(pmin(
                k - birth.death.matrix[edges.at.k, "birth_time"],
                birth.death.matrix[edges.at.k, "death_time"] - k
            ))

            stability.scores[k - kmin + 1] <- avg.persistence * persistent.ratio * edge.stability
        }
    }

    opt.k <- kmin - 1 + which.max(stability.scores)

    list(
        stability.scores = stability.scores,
        k.values = kmin:kmax,
        opt.k = opt.k
    )
}

#' Plot Method for IkNN Stability Metrics
#'
#' @param x An object returned by \code{\link{compute.stability.metrics}}.
#' @param ... Passed to \code{\link{plot.IkNNgraphs}}.
#'
#' @export
plot.iknn_stability_metrics <- function(x, ...) {
    plot.IkNNgraphs(x, ...)
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
    type <- match.arg(type, choices = types)

    if (type != "diag") {
        stop("Only type = 'diag' is currently supported.")
    }

    if (!"k.values" %in% names(x)) {
        stop("k.values not in x")
    }

    k.edge <- x$k.values
    if (length(k.edge) < 1L) {
        stop("k.values must have positive length.")
    }

    ## Transition k values correspond to metrics between consecutive graphs (k -> k+1)
    k.tr <- if (length(k.edge) >= 2L) k.edge[-length(k.edge)] else integer(0)

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par), add = TRUE)

    if (setequal(diags, c("edist","edge","deg"))) {

        par(mfrow = c(1,3), mar = mar, mgp = mgp, tcl = tcl)

        ## ---- Edit distance (transition metric) ----
        if (!"edit.distances" %in% names(x)) stop("edit.distances not in x")
        if (length(x$edit.distances) != length(k.tr)) {
            stop("Length mismatch: edit.distances must have length length(k.values)-1.")
        }

        plot(k.tr, x$edit.distances, las = 1, type = "b", xlab = "", ylab = "")
        mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
        mtext("Edit Distance", side = 2, line = yline, outer = FALSE)

        if (with.pwlm && "edit.distances.pwlm" %in% names(x)) {
            plot(x$edit.distances.pwlm, add = TRUE, col = "red")
            if ("edit.distances.breakpoint" %in% names(x)) {
                abline(v = x$edit.distances.breakpoint, lty = 2, col = breakpoint.col)
            }
        }
        if (with.lmin && "edit.distances.lmin" %in% names(x)) {
            abline(v = x$edit.distances.lmin, lty = 2, col = lmin.col)
        }

        ## ---- Edge count (per-graph metric) ----
        if (!"n.edges.in.pruned.graph" %in% names(x)) stop("n.edges.in.pruned.graph not in x")
        if (length(x$n.edges.in.pruned.graph) != length(k.edge)) {
            stop("Length mismatch: n.edges.in.pruned.graph must have length length(k.values).")
        }

        plot(k.edge, x$n.edges.in.pruned.graph, las = 1, type = "b", xlab = "", ylab = "")
        mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
        mtext("Num. Edges in Pruned Graph", side = 2, line = yline, outer = FALSE)

        if (with.pwlm && "n.edges.in.pruned.graph.pwlm" %in% names(x)) {
            plot(x$n.edges.in.pruned.graph.pwlm, add = TRUE, col = "red")
            if ("n.edges.in.pruned.graph.breakpoint" %in% names(x)) {
                abline(v = x$n.edges.in.pruned.graph.breakpoint, lty = 2, col = breakpoint.col)
            }
        }
        if (with.lmin && "n.edges.in.pruned.graph.lmin" %in% names(x)) {
            abline(v = x$n.edges.in.pruned.graph.lmin, lty = 2, col = lmin.col)
        }

        ## ---- JS divergence of degree profiles (transition metric) ----
        if (!"js.div" %in% names(x)) stop("js.div not in x")
        if (length(x$js.div) != length(k.tr)) {
            stop("Length mismatch: js.div must have length length(k.values)-1.")
        }

        plot(k.tr, x$js.div, las = 1, type = "b", xlab = "", ylab = "")
        mtext("Number of Nearest Neighbors (k)", side = 1, line = xline, outer = FALSE)
        mtext("JS Divergence (Degrees)", side = 2, line = yline, outer = FALSE)

        if (with.pwlm && "js.div.pwlm" %in% names(x)) {
            plot(x$js.div.pwlm, add = TRUE, col = "red")
            if ("js.div.breakpoint" %in% names(x)) {
                abline(v = x$js.div.breakpoint, lty = 2, col = breakpoint.col)
            }
        }
        if (with.lmin && "js.div.lmin" %in% names(x)) {
            abline(v = x$js.div.lmin, lty = 2, col = lmin.col)
        }

    } else {
        stop("Currently supported diags combination is exactly c('edist','edge','deg').")
    }

    invisible(TRUE)
}

#' Compute edit distances between consecutive graphs
#'
#' @description
#' Computes a simple graph edit distance between consecutive graphs in a sequence.
#' For each adjacent pair \eqn{(G_i, G_{i+1})}, the distance is defined as the size
#' of the symmetric difference between their undirected edge sets.
#'
#' Each graph is expected to provide an adjacency list \code{adj_list} using
#' \strong{1-based} vertex indices.
#'
#' @param graphs A list of graph objects. Each element must contain an
#'   \code{adj_list} component, which is a list of integer neighbor vectors.
#'
#' @return A numeric vector of length \code{length(graphs) - 1L}. Entry \code{i}
#'   gives the edit distance between \code{graphs[[i]]} and \code{graphs[[i + 1L]]}.
#'   If \code{length(graphs) < 2L}, returns \code{numeric(0)}.
#'
#' @keywords internal
#' @noRd
internal.compute.edit.distances <- function(graphs) {
    ## Graph edit distance between consecutive graphs = symmetric difference of edge sets
    ## graphs: list of graph objects with adj_list (1-based)
    n <- length(graphs)
    if (n < 2) return(numeric(0))

    edge.set <- function(adj.list) {
        ## Build undirected edge keys "i-j" with i<j
        keys <- character(0)
        for (i in seq_along(adj.list)) {
            nbrs <- adj.list[[i]]
            if (length(nbrs) == 0) next
            j <- nbrs[nbrs > i]
            if (length(j) > 0) {
                keys <- c(keys, paste(i, j, sep = "-"))
            }
        }
        unique(keys)
    }

    out <- numeric(n - 1)
    for (i in seq_len(n - 1)) {
        e1 <- edge.set(graphs[[i]]$adj_list)
        e2 <- edge.set(graphs[[i + 1]]$adj_list)
        out[i] <- length(setdiff(e1, e2)) + length(setdiff(e2, e1))
    }
    out
}

#' Trim rows of \code{X} to the main connected component
#'
#' @description
#' Given a data matrix \code{X} whose rows correspond to vertices in a graph,
#' trims \code{X} to the largest connected component of the graph defined by an
#' adjacency list.
#'
#' Connected components are computed by \code{graph.connected.components()}, which
#' must return an integer component label for each vertex. Component labels need
#' not be contiguous.
#'
#' @param X A numeric matrix. Rows correspond to vertices.
#' @param adj.list A list representing the graph adjacency list using \strong{1-based}
#'   vertex indices.
#' @param verbose Logical. If \code{TRUE}, prints the number of vertices before and
#'   after trimming.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{X}: The trimmed matrix \code{X[in.main, , drop = FALSE]}.
#'   \item \code{kept}: A logical vector of length \code{nrow(X)} indicating which
#'   rows/vertices were retained.
#' }
#'
#' @keywords internal
#' @noRd
trim.X.to.main.cc <- function(X, adj.list, verbose = FALSE) {
    cc <- graph.connected.components(adj.list)
    cc.tbl <- table(cc)
    main.cc <- as.integer(names(sort(cc.tbl, decreasing = TRUE)[1]))
    in.main <- (cc == main.cc)

    if (verbose) {
        cat("Trimming to main connected component:\n")
        cat("  vertices before:", nrow(X), "\n")
        cat("  vertices kept  :", sum(in.main), "\n")
    }

    list(
        X = X[in.main, , drop = FALSE],
        kept = in.main
    )
}

#' Pick the smallest k within eps of the global optimum (max or min), optionally requiring a local extremum
#'
#' @description
#' Selects a k value from a candidate grid based on a scalar metric evaluated at each k.
#' For \code{direction="max"}, chooses the smallest k whose metric is within
#' \code{(1 - eps)} of the global maximum. For \code{direction="min"}, chooses the
#' smallest k whose metric is within \code{(1 + eps)} of the global minimum.
#'
#' Optionally restricts to local extrema (maxima/minima) using a sliding window.
#' Handles non-finite prefixes/suffixes safely (no \code{max(which(...))} pitfalls).
#'
#' @param metric Numeric vector of metric values (one per candidate k, aligned with \code{k.values}).
#' @param k.values Integer/numeric vector of candidate k values. If NULL, uses \code{seq_along(metric)}.
#' @param eps Non-negative numeric tolerance. For max: keep metric >= (1-eps)*max.
#'   For min: keep metric <= (1+eps)*min.
#' @param direction Character: "max" or "min".
#' @param idx.ok Optional integer indices specifying eligible positions in \code{metric}.
#' @param k.min Optional lower bound on k (inclusive).
#' @param k.max Optional upper bound on k (inclusive).
#' @param require.local.extremum Logical; if TRUE, restrict candidates to local extrema.
#' @param window Integer >= 1, neighborhood half-width for local extremum detection.
#' @param return.details Logical; if TRUE, returns a list with details and diagnostics.
#'
#' @return If \code{return.details=FALSE}, returns the selected k (scalar) or NA.
#'   If \code{return.details=TRUE}, returns a list with fields:
#'   \code{k.opt}, \code{idx.opt}, \code{threshold}, \code{idx.candidates}, \code{idx.local}.
#'
#' @export
pick.k.within.eps.global.max <- function(metric,
                                         k.values = NULL,
                                         eps = 0.05,
                                         direction = c("max", "min"),
                                         idx.ok = NULL,
                                         k.min = -Inf,
                                         k.max = Inf,
                                         require.local.extremum = FALSE,
                                         window = 1L,
                                         return.details = FALSE) {
    ## ---- validate ----
    if (missing(metric) || is.null(metric)) stop("`metric` must be provided.")
    metric <- as.double(metric)
    n <- length(metric)
    if (n < 1L) stop("`metric` must have length >= 1.")
    if (is.null(k.values)) {
        k.values <- seq_len(n)
    } else {
        if (length(k.values) != n) stop("`k.values` must have the same length as `metric`.")
    }
    k.values <- as.double(k.values)

    direction <- match.arg(direction)
    if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps < 0) {
        stop("`eps` must be a single finite number >= 0.")
    }
    window <- as.integer(window)
    if (!is.finite(window) || window < 1L) stop("`window` must be an integer >= 1.")

    ## ---- eligibility mask ----
    finite <- is.finite(metric) & is.finite(k.values)
    keep <- finite & (k.values >= k.min) & (k.values <= k.max)
    if (!is.null(idx.ok)) {
        idx.ok <- as.integer(idx.ok)
        idx.ok <- idx.ok[idx.ok >= 1L & idx.ok <= n]
        keep2 <- rep(FALSE, n)
        keep2[idx.ok] <- TRUE
        keep <- keep & keep2
    }

    idx.keep <- which(keep)
    if (length(idx.keep) == 0L) {
        if (isTRUE(return.details)) {
            return(list(k.opt = NA_real_, idx.opt = NA_integer_, threshold = NA_real_,
                        idx.candidates = integer(0), idx.local = integer(0)))
        }
        return(NA_real_)
    }

    ## ---- compute global optimum on eligible set ----
    m.keep <- metric[idx.keep]
    if (direction == "max") {
        m.opt <- max(m.keep)
        thr <- (1 - eps) * m.opt
        idx.cand <- idx.keep[metric[idx.keep] >= thr]
    } else {
        m.opt <- min(m.keep)
        thr <- (1 + eps) * m.opt
        idx.cand <- idx.keep[metric[idx.keep] <= thr]
    }

    if (length(idx.cand) == 0L) {
        if (isTRUE(return.details)) {
            return(list(k.opt = NA_real_, idx.opt = NA_integer_, threshold = thr,
                        idx.candidates = integer(0), idx.local = integer(0)))
        }
        return(NA_real_)
    }

    ## ---- optionally restrict to local extrema ----
    idx.local <- integer(0)
    if (isTRUE(require.local.extremum)) {
        ## local extremum test uses the full metric vector but respects finiteness
        is.local <- rep(FALSE, n)
        for (ii in idx.cand) {
            lo <- max(1L, ii - window)
            hi <- min(n, ii + window)
            nb <- metric[lo:hi]
            nb.ok <- is.finite(nb)
            if (!any(nb.ok)) next
            nb <- nb[nb.ok]

            if (direction == "max") {
                ## plateau-safe: point is local max if it attains neighborhood max
                is.local[ii] <- isTRUE(all(metric[ii] >= nb))
            } else {
                is.local[ii] <- isTRUE(all(metric[ii] <= nb))
            }
        }
        idx.local <- idx.cand[is.local[idx.cand]]

        ## if no local extrema among candidates, fall back to candidates
        if (length(idx.local) == 0L) idx.local <- idx.cand
    } else {
        idx.local <- idx.cand
    }

    ## ---- choose smallest k among remaining ----
    ## (ties broken by smallest k; if multiple entries share same k, earliest index)
    k.sub <- k.values[idx.local]
    k.min.val <- min(k.sub, na.rm = TRUE)
    idx.opt <- idx.local[which(k.sub == k.min.val)[1L]]
    k.opt <- k.values[idx.opt]

    if (isTRUE(return.details)) {
        return(list(k.opt = k.opt, idx.opt = idx.opt, threshold = thr,
                    idx.candidates = idx.cand, idx.local = idx.local))
    }
    k.opt
}

#' Build ikNN graphs over a k range and select k by edit-distance and/or CST mixing
#'
#' @description
#' Builds a sequence of ikNN-derived graphs over \code{kmin:kmax} using \code{create.iknn.graphs()},
#' optionally trims to the largest connected component if the sequence is disconnected, and selects
#' an "optimal" k using:
#' \itemize{
#'   \item structural stability (minimum edit distance between consecutive graphs) and/or
#'   \item external label coherence (e.g., CST mixing metrics; homophily/assortativity effect or z-score)
#' }
#'
#' This function fixes common pitfalls:
#' \enumerate{
#'   \item Trimming uses the graph with the largest LCC in the range (not necessarily kmin).
#'   \item After trimming, graphs are rebuilt with the same PCA settings (\code{pca.dim}, \code{variance.explained}).
#'   \item Mixing selection defaults to metrics that can peak (e.g., effect/z/adjusted), rather than raw homophily.
#'   \item Edge list conversion is duplicate-safe; optional \code{igraph::simplify()} is applied.
#'   \item For edge lengths, weights can be converted to affinities using a fixed sigma across k.
#' }
#'
#' @param X Numeric matrix (samples in rows, features in columns).
#' @param kmin Integer scalar, minimum k.
#' @param kmax Integer scalar, maximum k.
#' @param method Character: "edit", "mixing", "both", or "none".
#' @param pca.dim Optional integer for PCA dimension inside \code{create.iknn.graphs()}.
#' @param variance.explained Optional numeric in (0,1] for PCA variance inside \code{create.iknn.graphs()}.
#' @param trim.disconnected Logical; if TRUE and the sequence fails connectivity requirements,
#'   trims to largest connected component of the best (largest-LCC) graph and rebuilds graphs.
#'
#' @param edit.min.lcc.frac Numeric in (0,1], default 1.0 (fully connected tail).
#' @param edit.eps Numeric >=0, relative tolerance for edit-distance selection (minimization).
#'
#' @param labels Optional categorical labels (e.g., CST), either length nrow(X) or named by rownames(X).
#'   Required for method "mixing" or "both".
#' @param perm.blocks Optional vector/factor for block-wise label permutation (same length/names as labels).
#' @param mixing.metric Character, one of:
#'   \itemize{
#'     \item "homophily" "homophily.z" "homophily.effect" "homophily.adjusted"
#'     \item "assortativity" "assortativity.z" "assortativity.effect"
#'     \item "conductance.median" "conductance.wmean" (minimization metrics)
#'   }
#' @param mixing.min.lcc.frac Numeric in (0,1], minimum LCC fraction required for mixing evaluation.
#' @param mixing.eps Numeric >=0 tolerance used with \code{pick.k.within.eps.global.max()}.
#' @param mixing.require.local.extremum Logical; if TRUE, prefers local extrema for mixing metric.
#' @param mixing.window Integer neighborhood for local extremum detection.
#' @param n.perm Integer number of permutations for mixing nulls.
#'
#' @param use.edge.weights Logical; if TRUE, uses graph edge weights when available.
#' @param weights.are.edge.lengths Logical; if TRUE, interpret weights as lengths and convert to affinities.
#' @param affinity.method "exp" or "inv". For exp: exp(-(d/sigma)^2). For inv: 1/(d+eps).
#' @param affinity.sigma Optional positive scalar; if NULL, estimated once from a reference k graph.
#' @param affinity.sigma.from Character: "k.cc.mixing" or "k.trim" or "k.max.lcc" to choose reference graph.
#' @param affinity.eps Small positive scalar for inverse affinity.
#' @param simplify.multiple Logical; if TRUE, simplifies multi-edges/loops and merges weights.
#' @param seed Numeric; seed of a random number generator.
#' @param n.cores Integer number of cores passed to \code{create.iknn.graphs()}.
#' @param verbose Logical.
#' @param ... Additional arguments forwarded to \code{create.iknn.graphs()} ONLY.
#'
#' @return Object of class "build_iknn_graphs_and_selectk" with fields:
#'   \itemize{
#'     \item X.graphs
#'     \item k.values
#'     \item connectivity (data.frame)
#'     \item edit (data.frame or NULL)
#'     \item mixing (data.frame or NULL)
#'     \item k.cc.edit, k.opt.edit
#'     \item k.cc.mixing, k.opt.mixing
#'     \item trim (list)
#'     \item params (list)
#'   }
#'
#' @export
build.iknn.graphs.and.selectk <- function(X,
                                         kmin,
                                         kmax,
                                         method = c("both", "edit", "mixing", "none"),
                                         pca.dim = 100,
                                         variance.explained = 0.99,
                                         trim.disconnected = TRUE,
                                         edit.min.lcc.frac = 1.0,
                                         edit.eps = 0.05,
                                         labels = NULL,
                                         perm.blocks = NULL,
                                         mixing.metric = c("homophily.effect",
                                                          "homophily.z",
                                                          "homophily.adjusted",
                                                          "assortativity.effect",
                                                          "assortativity.z",
                                                          "homophily",
                                                          "assortativity",
                                                          "conductance.median",
                                                          "conductance.wmean"),
                                         mixing.min.lcc.frac = 0.98,
                                         mixing.eps = 0.05,
                                         mixing.require.local.extremum = TRUE,
                                         mixing.window = 1L,
                                         n.perm = 200L,
                                         use.edge.weights = TRUE,
                                         weights.are.edge.lengths = FALSE,
                                         affinity.method = c("exp", "inv"),
                                         affinity.sigma = NULL,
                                         affinity.sigma.from = c("k.cc.mixing", "k.max.lcc", "k.trim"),
                                         affinity.eps = 1e-8,
                                         simplify.multiple = TRUE,
                                         seed = 1L,
                                         n.cores = 1L,
                                         verbose = TRUE,
                                         ...) {

    ## ---- validate X ----
    if (missing(X) || is.null(X)) stop("`X` must be provided.")
    X <- tryCatch(as.matrix(X), error = function(e) NULL)
    if (is.null(X)) stop("`X` must be coercible to a matrix via as.matrix().")
    suppressWarnings(storage.mode(X) <- "double")
    if (!is.numeric(X)) stop("`X` must be numeric (or coercible to numeric).")
    if (nrow(X) < 5L) stop("`X` must have at least 5 rows.")
    if (ncol(X) < 1L) stop("`X` must have at least 1 column.")

    sample.ids <- rownames(X)
    if (is.null(sample.ids)) sample.ids <- as.character(seq_len(nrow(X)))

    kmin <- as.integer(kmin)
    kmax <- as.integer(kmax)
    if (!is.finite(kmin) || !is.finite(kmax) || kmin < 1L) stop("`kmin` must be an integer >= 1.")
    if (kmax < kmin) stop("`kmax` must be >= kmin.")
    if (kmax >= nrow(X)) {
        kmax <- nrow(X) - 1L
        if (isTRUE(verbose)) {
            cat("NOTE: reducing kmax to nrow(X)-1 =", kmax, "\n")
        }
    }

    method <- match.arg(method)
    mixing.metric <- match.arg(mixing.metric)
    affinity.method <- match.arg(affinity.method)
    affinity.sigma.from <- match.arg(affinity.sigma.from)

    n.cores <- as.integer(n.cores)
    if (!is.finite(n.cores) || n.cores < 1L) stop("`n.cores` must be an integer >= 1.")

    if (!is.numeric(edit.min.lcc.frac) || length(edit.min.lcc.frac) != 1L ||
        !is.finite(edit.min.lcc.frac) || edit.min.lcc.frac <= 0 || edit.min.lcc.frac > 1) {
        stop("`edit.min.lcc.frac` must be in (0,1].")
    }
    if (!is.numeric(mixing.min.lcc.frac) || length(mixing.min.lcc.frac) != 1L ||
        !is.finite(mixing.min.lcc.frac) || mixing.min.lcc.frac <= 0 || mixing.min.lcc.frac > 1) {
        stop("`mixing.min.lcc.frac` must be in (0,1].")
    }

    ## ---- helper: align named labels to sample.ids robustly ----
    align.labels.to.sample.ids <- function(labels, sample.ids, min.match = 10L) {
        if (is.null(names(labels))) stop("`labels` must be named (or length nrow(X) without names).")
        min.match <- as.integer(min.match)
        if (!is.finite(min.match) || min.match < 0L) stop("`min.match` must be an integer >= 0.")

        ## strict match
        lab0 <- labels[sample.ids]
        n.match0 <- sum(!is.na(lab0))

        if (n.match0 >= min.match) {
            return(list(labels.aligned = as.character(lab0), used.make.names = FALSE))
        }

        ## fallback: make.names normalization
        nm.sample <- make.names(sample.ids, unique = FALSE)
        nm.labels <- make.names(names(labels), unique = FALSE)
        labels.mn <- labels
        names(labels.mn) <- nm.labels
        lab1 <- labels.mn[nm.sample]
        n.match1 <- sum(!is.na(lab1))

        if (n.match1 >= min.match) {
            warning("Label alignment succeeded only after make.names() normalization. ",
                    "Consider normalizing rownames(X) and names(labels) consistently upstream.")
            return(list(labels.aligned = as.character(lab1), used.make.names = TRUE))
        }

        stop(
            "Label alignment failed: too few matches between rownames(X) and names(labels).\n",
            "  matches strict: ", n.match0, " / ", length(sample.ids), "\n",
            "  matches make.names: ", n.match1, " / ", length(sample.ids), "\n",
            "  head(sample.ids): ", paste(utils::head(sample.ids, 5), collapse = ", "), "\n",
            "  head(names(labels)): ", paste(utils::head(names(labels), 5), collapse = ", "), "\n",
            "Fix: ensure rownames(X) equals names(labels) (or pass labels as an unnamed vector in row order)."
        )
    }

    need.mixing <- method %in% c("mixing", "both")
    labels.aligned <- NULL
    blocks.aligned <- NULL

    if (need.mixing) {
        if (is.null(labels)) stop("`labels` must be provided when method includes 'mixing'.")

        if (!is.null(names(labels))) {
            ## if X has no rownames but labels are named and length matches, adopt label names
            if (is.null(rownames(X)) && length(labels) == nrow(X) && length(unique(names(labels))) == nrow(X)) {
                rownames(X) <- names(labels)
                sample.ids <- rownames(X)
            }

            al <- align.labels.to.sample.ids(labels, sample.ids, min.match = 10L)
            labels.aligned <- al$labels.aligned
        } else {
            if (length(labels) != nrow(X)) stop("`labels` must have length nrow(X) or be named by rownames(X).")
            labels.aligned <- as.character(labels)
            names(labels.aligned) <- sample.ids
        }

        ## optional: warn if many NA
        n.lab <- sum(!is.na(labels.aligned))
        if (n.lab < 10L) stop("Too few non-NA labels after alignment (n=", n.lab, ").")

        if (!is.null(perm.blocks)) {
            if (!is.null(names(perm.blocks))) {
                blocks.aligned <- perm.blocks[sample.ids]
            } else {
                if (length(perm.blocks) != nrow(X)) stop("`perm.blocks` must have length nrow(X) or be named.")
                blocks.aligned <- perm.blocks
                names(blocks.aligned) <- sample.ids
            }
        }
    }

    ## ---- helper: build igraph + edge codes from (adj.list, weight.list) ----
    adjlist.to.edge.mat <- function(adj.list, weight.list = NULL, n) {
        has.w <- !is.null(weight.list)
        e1 <- integer(0); e2 <- integer(0); ew <- numeric(0)

        for (i in seq_len(n)) {
            nb <- adj.list[[i]]
            if (length(nb) == 0L) next
            nb <- as.integer(nb)

            if (!has.w) {
                jj <- nb[nb > i]
                if (length(jj) > 0L) {
                    e1 <- c(e1, rep.int(i, length(jj)))
                    e2 <- c(e2, jj)
                    ew <- c(ew, rep.int(1.0, length(jj)))
                }
            } else {
                wv <- as.double(weight.list[[i]])
                if (length(wv) != length(nb)) stop("weight.list[[i]] length must match adj.list[[i]].")
                keep <- (nb > i)
                if (any(keep)) {
                    e1 <- c(e1, rep.int(i, sum(keep)))
                    e2 <- c(e2, nb[keep])
                    ew <- c(ew, wv[keep])
                }
            }
        }

        if (length(e1) == 0L) {
            return(list(edge.mat = matrix(integer(0), ncol = 2), weights = numeric(0)))
        }

        edge.mat <- cbind(e1, e2)
        code <- (edge.mat[, 1] - 1L) * n + edge.mat[, 2]

        ## merge duplicates defensively
        if (length(code) != length(unique(code))) {
            ## combine weights depending on semantics
            comb <- if (isTRUE(weights.are.edge.lengths)) min else max
            w.by.code <- tapply(ew, code, comb)
            code.u <- as.integer(names(w.by.code))
            ## decode back to (i,j)
            i.u <- (code.u - 1L) %/% n + 1L
            j.u <- (code.u - 1L) %% n + 1L
            edge.mat <- cbind(as.integer(i.u), as.integer(j.u))
            ew <- as.double(w.by.code)
        }

        list(edge.mat = edge.mat, weights = ew)
    }

    edge.codes.from.graph <- function(g.obj, n) {
        el <- adjlist.to.edge.mat(g.obj$adj_list, g.obj$weight_list, n = n)
        if (nrow(el$edge.mat) == 0L) return(integer(0))
        code <- (el$edge.mat[, 1] - 1L) * n + el$edge.mat[, 2]
        sort(unique(as.integer(code)))
    }

    jaccard.distance.codes <- function(a, b) {
        a <- as.integer(a); b <- as.integer(b)
        if (length(a) == 0L && length(b) == 0L) return(0)
        if (length(a) == 0L || length(b) == 0L) return(1)
        inter <- sum(!is.na(match(a, b)))
        uni <- length(a) + length(b) - inter
        if (uni <= 0L) return(0)
        1 - (inter / uni)
    }

    make.igraph.from.graph <- function(g.obj, n, w.use) {
        el <- adjlist.to.edge.mat(g.obj$adj_list, g.obj$weight_list, n = n)
        g <- igraph::make_empty_graph(n = n, directed = FALSE)
        if (nrow(el$edge.mat) > 0L) {
            g <- igraph::add_edges(g, as.vector(t(el$edge.mat)))
            if (!is.null(w.use) && length(w.use) == nrow(el$edge.mat)) {
                igraph::E(g)$weight <- as.double(w.use)
            } else {
                ## keep raw weights if present
                if (length(el$weights) == nrow(el$edge.mat)) igraph::E(g)$weight <- as.double(el$weights)
            }
        }
        g
    }

    estimate.sigma.from.lengths <- function(d) {
        d <- as.double(d)
        d <- d[is.finite(d) & d > 0]
        if (length(d) == 0L) return(1.0)
        stats::median(d)
    }

    lengths.to.affinity <- function(d, sigma, method = "exp") {
        d <- as.double(d)
        if (!is.finite(sigma) || sigma <= 0) sigma <- estimate.sigma.from.lengths(d)
        if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
        if (method == "exp") {
            w <- exp(- (d / sigma)^2)
        } else {
            w <- 1 / (d + affinity.eps)
        }
        w[!is.finite(w)] <- 0
        w[w < 0] <- 0
        w
    }

    ## ---- call mixing stats robustly across possible function signatures ----
    call.mixing.stats <- function(igraph.obj, labels.vec, blocks.vec = NULL, w.vec = NULL, seed = 1L) {
        fml <- tryCatch(names(formals(cst.graph.mixing.stats)), error = function(e) character(0))

        args <- list(igraph.obj = igraph.obj, labels = labels.vec, n.perm = as.integer(n.perm), seed = as.integer(seed))
        if ("perm.blocks" %in% fml) args$perm.blocks <- blocks.vec

        ## Prefer passing edge.weights if supported; else rely on E(g)$weight
        if (!is.null(w.vec)) {
            if ("edge.weights" %in% fml) {
                args$edge.weights <- w.vec
            } else {
                igraph::E(igraph.obj)$weight <- w.vec
            }
        }

        do.call(cst.graph.mixing.stats, args)
    }

    extract.metric <- function(ms, metric.name, igraph.obj = NULL, labels.vec = NULL) {
        ## helper to safely extract from old/new return structures
        get1 <- function(x, path) {
            cur <- x
            for (nm in path) {
                if (is.null(cur) || is.null(cur[[nm]])) return(NULL)
                cur <- cur[[nm]]
            }
            cur
        }

        if (metric.name == "homophily") return(ms$homophily)

        if (metric.name == "homophily.z") {
            z <- get1(ms, c("permutation", "homophily", "z"))
            if (is.null(z)) z <- get1(ms, c("permutation", "homophily.z", "z"))
            if (is.null(z)) return(NA_real_)
            return(as.double(z))
        }

        if (metric.name == "homophily.effect") {
            eff <- get1(ms, c("permutation", "homophily", "effect"))
            if (!is.null(eff)) return(eff)
            mu <- get1(ms, c("permutation", "homophily.z", "mu"))
            if (!is.null(mu) && is.finite(ms$homophily)) return(ms$homophily - mu)
            return(NA_real_)
        }

        if (metric.name == "homophily.adjusted") {
            adj <- ms$homophily.adjusted
            if (!is.null(adj)) return(adj)

            ## fallback: compute baseline using strength-weighted label proportions if graph provided
            if (!is.null(igraph.obj) && !is.null(labels.vec)) {
                w <- if ("weight" %in% igraph::edge_attr_names(igraph.obj)) igraph::E(igraph.obj)$weight else NULL
                s <- igraph::strength(igraph.obj, weights = w)
                ok <- !is.na(labels.vec)
                if (sum(s[ok]) > 0) {
                    p <- tapply(s[ok], labels.vec[ok], sum)
                    p <- p / sum(s[ok])
                    h0 <- sum(as.numeric(p)^2)
                    if (is.finite(h0) && h0 < 1 && is.finite(ms$homophily)) return((ms$homophily - h0) / (1 - h0))
                }
            }
            return(NA_real_)
        }

        if (metric.name == "assortativity") return(ms$assortativity)

        if (metric.name == "assortativity.z") {
            z <- get1(ms, c("permutation", "assortativity", "z"))
            if (is.null(z)) z <- get1(ms, c("permutation", "assortativity.z", "z"))
            if (is.null(z)) return(NA_real_)
            return(as.double(z))
        }

        if (metric.name == "assortativity.effect") {
            eff <- get1(ms, c("permutation", "assortativity", "effect"))
            if (!is.null(eff)) return(eff)
            mu <- get1(ms, c("permutation", "assortativity.z", "mu"))
            if (!is.null(mu) && is.finite(ms$assortativity)) return(ms$assortativity - mu)
            return(NA_real_)
        }

        if (metric.name == "conductance.median") {
            v <- get1(ms, c("conductance.summary", "conductance.median"))
            if (is.null(v)) return(NA_real_)
            return(as.double(v))
        }
        if (metric.name == "conductance.wmean") {
            v <- get1(ms, c("conductance.summary", "conductance.vol.weighted.mean"))
            if (is.null(v)) return(NA_real_)
            return(as.double(v))
        }

        NA_real_
    }

    metric.direction.default <- function(metric.name) {
        if (metric.name %in% c("conductance.median", "conductance.wmean")) "min" else "max"
    }

    ## ---- build graphs (initial) ----
    if (!exists("create.iknn.graphs", mode = "function")) {
        stop("create.iknn.graphs() not found in the environment/namespace.")
    }

    X.graphs <- create.iknn.graphs(
        X,
        kmin = kmin,
        kmax = kmax,
        pca.dim = pca.dim,
        variance.explained = variance.explained,
        compute.full = TRUE,
        n.cores = n.cores,
        verbose = verbose,
        ...
    )

    ## ---- locate k.values and graph list ----
    k.values <- NULL
    if (!is.null(X.graphs$k_statistics)) k.values <- X.graphs$k_statistics[,"k"]
    if (!is.null(X.graphs$k.values)) k.values <- X.graphs$k.values
    if (is.null(k.values) && !is.null(X.graphs$k)) k.values <- X.graphs$k
    if (is.null(k.values)) k.values <- seq.int(kmin, kmax)

    g.list <- NULL
    if (!is.null(X.graphs$geom_pruned_graphs)) g.list <- X.graphs$geom_pruned_graphs
    if (is.null(g.list)) stop("X.graphs$geom_pruned_graphs not found; cannot proceed.")

    if (length(g.list) != length(k.values)) {
        stop("Length mismatch: geom_pruned_graphs and k.values.")
    }

    ## ---- connectivity diagnostics ----
    compute.connectivity <- function(g.list, k.values, n) {
        n.comp <- integer(length(k.values))
        lcc.size <- integer(length(k.values))
        lcc.frac <- numeric(length(k.values))
        n.edges <- integer(length(k.values))

        for (i in seq_along(k.values)) {
            g.obj <- g.list[[i]]
            el <- adjlist.to.edge.mat(g.obj$adj_list, g.obj$weight_list, n = n)
            n.edges[i] <- nrow(el$edge.mat)

            gi <- igraph::make_empty_graph(n = n, directed = FALSE)
            if (nrow(el$edge.mat) > 0L) gi <- igraph::add_edges(gi, as.vector(t(el$edge.mat)))

            comp <- igraph::components(gi)
            n.comp[i] <- comp$no
            lcc.size[i] <- max(comp$csize)
            lcc.frac[i] <- lcc.size[i] / n
        }

        data.frame(
            k = as.integer(k.values),
            n.edges = n.edges,
            n.components = n.comp,
            lcc.size = lcc.size,
            lcc.frac = lcc.frac
        )
    }

    conn <- compute.connectivity(g.list, k.values, n = nrow(X))

    ## ---- connected-tail helper ----
    find.k.cc <- function(conn.df, min.lcc.frac) {
        bad <- which(conn.df$lcc.frac < min.lcc.frac)
        if (length(bad) == 0L) return(conn.df$k[1L])
        last.bad <- max(bad)
        if (last.bad >= nrow(conn.df)) return(NA_integer_)
        k.cc <- conn.df$k[last.bad + 1L]
        ## verify tail
        ok.tail <- which(conn.df$k >= k.cc)
        if (!all(conn.df$lcc.frac[ok.tail] >= min.lcc.frac)) return(NA_integer_)
        k.cc
    }

    k.cc.edit <- find.k.cc(conn, min.lcc.frac = edit.min.lcc.frac)

    ## ---- trimming if needed ----
    trim.info <- list(trimmed = FALSE, keep.idx = seq_len(nrow(X)), dropped.idx = integer(0), k.trim = NA_integer_)
    if (isTRUE(trim.disconnected) && (is.na(k.cc.edit) && method %in% c("edit", "both"))) {
        ## choose k with max LCC (tie -> smallest k)
        best.idx <- which(conn$lcc.size == max(conn$lcc.size))
        best.idx <- best.idx[which.min(conn$k[best.idx])]
        k.trim <- conn$k[best.idx]
        trim.info$k.trim <- k.trim

        ## compute LCC vertices at k.trim
        g.trim <- g.list[[best.idx]]
        el <- adjlist.to.edge.mat(g.trim$adj_list, g.trim$weight_list, n = nrow(X))
        gi <- igraph::make_empty_graph(n = nrow(X), directed = FALSE)
        if (nrow(el$edge.mat) > 0L) gi <- igraph::add_edges(gi, as.vector(t(el$edge.mat)))
        comp <- igraph::components(gi)
        keep <- which(comp$membership == which.max(comp$csize))

        if (length(keep) < 5L) stop("Trimming would leave too few vertices; aborting.")

        if (isTRUE(verbose)) {
            cat("Trimming to largest CC at k =", k.trim, " (n=", length(keep), " of ", nrow(X), ")\n", sep = "")
        }

        X <- X[keep, , drop = FALSE]
        sample.ids <- sample.ids[keep]

        if (need.mixing) {
            labels.aligned <- labels.aligned[sample.ids]
            if (!is.null(blocks.aligned)) blocks.aligned <- blocks.aligned[sample.ids]
        }

        trim.info$trimmed <- TRUE
        trim.info$keep.idx <- keep
        trim.info$dropped.idx <- setdiff(seq_len(nrow(gi)), keep)

        ## rebuild graphs with same PCA settings
        if (kmax >= nrow(X)) kmax <- nrow(X) - 1L
        X.graphs <- create.iknn.graphs(
            X,
            kmin = kmin,
            kmax = kmax,
            pca.dim = pca.dim,
            variance.explained = variance.explained,
            compute.full = TRUE,
            n.cores = n.cores,
            verbose = verbose,
            ...
        )

        ## refresh extracted fields
        if (!is.null(X.graphs$k.values)) k.values <- X.graphs$k.values else if (!is.null(X.graphs$k)) k.values <- X.graphs$k else k.values <- seq.int(kmin, kmax)
        if (!is.null(X.graphs$geom_pruned_graphs)) g.list <- X.graphs$geom_pruned_graphs else stop("Rebuilt X.graphs missing geom_pruned_graphs.")

        conn <- compute.connectivity(g.list, k.values, n = nrow(X))
        k.cc.edit <- find.k.cc(conn, min.lcc.frac = edit.min.lcc.frac)
    }

    ## k.cc for mixing eligibility (can differ from edit)
    k.cc.mixing <- find.k.cc(conn, min.lcc.frac = mixing.min.lcc.frac)

    ## ---- edit-distance curve + selection ----
    edit.df <- NULL
    k.opt.edit <- NA_integer_

    pick.k.within.eps.global.min.internal <- function(metric, k.values, eps = 0.05, idx.ok = NULL) {
        pick.k.within.eps.global.max(metric = metric, k.values = k.values,
                                     eps = eps, direction = "min", idx.ok = idx.ok,
                                     require.local.extremum = FALSE, window = 1L,
                                     return.details = FALSE)
    }

    if (method %in% c("edit", "both")) {
        if (is.na(k.cc.edit)) stop("No connected tail found for edit selection; consider trim.disconnected=TRUE or relax edit.min.lcc.frac.")

        ## compute edge codes per k
        edge.codes <- vector("list", length(k.values))
        for (i in seq_along(k.values)) {
            edge.codes[[i]] <- edge.codes.from.graph(g.list[[i]], n = nrow(X))
        }

        edit.dist <- rep(NA_real_, length(k.values))
        for (i in seq_len(length(k.values) - 1L)) {
            edit.dist[i] <- jaccard.distance.codes(edge.codes[[i]], edge.codes[[i + 1L]])
        }

        edit.df <- data.frame(
            k = as.integer(k.values),
            edit.dist.to.next = edit.dist
        )

        idx.ok <- which(edit.df$k >= k.cc.edit & is.finite(edit.df$edit.dist.to.next))
        if (length(idx.ok) > 0L) {
            k.opt.edit <- as.integer(pick.k.within.eps.global.min.internal(edit.df$edit.dist.to.next, edit.df$k, eps = edit.eps, idx.ok = idx.ok))
        }
    }

    ## ---- mixing curve + selection ----
    mixing.df <- NULL
    k.opt.mixing <- NA_integer_
    sigma.used <- affinity.sigma

    if (need.mixing) {
        ## decide direction for selection
        dir0 <- metric.direction.default(mixing.metric)

        ## eligible indices for mixing evaluation
        idx.mix <- which(conn$lcc.frac >= mixing.min.lcc.frac)
        if (!is.na(k.cc.mixing)) idx.mix <- idx.mix[conn$k[idx.mix] >= k.cc.mixing]
        if (length(idx.mix) == 0L) stop("No k values satisfy mixing connectivity constraint; relax mixing.min.lcc.frac or trim.disconnected.")

        ## estimate sigma once if needed and weights are lengths
        if (isTRUE(use.edge.weights) && isTRUE(weights.are.edge.lengths) && is.null(affinity.sigma)) {
            idx.ref <- idx.mix[1L]
            if (affinity.sigma.from == "k.max.lcc") {
                idx.ref <- which(conn$lcc.size == max(conn$lcc.size))[1L]
            } else if (affinity.sigma.from == "k.trim" && isTRUE(trim.info$trimmed)) {
                idx.ref <- which(conn$k == trim.info$k.trim)
                if (length(idx.ref) == 0L) idx.ref <- idx.mix[1L]
            } else if (affinity.sigma.from == "k.cc.mixing") {
                idx.ref <- idx.mix[1L]
            }

            ## pull raw lengths from reference graph edgelist
            g.ref <- g.list[[idx.ref]]
            el.ref <- adjlist.to.edge.mat(g.ref$adj_list, g.ref$weight_list, n = nrow(X))
            sigma.used <- estimate.sigma.from.lengths(el.ref$weights)
            if (isTRUE(verbose)) cat("Estimated affinity.sigma =", signif(sigma.used, 5), "from k =", conn$k[idx.ref], "\n")
        }

        ## compute metrics across k
        k.out <- conn$k[idx.mix]
        val <- rep(NA_real_, length(idx.mix))
        val.z <- rep(NA_real_, length(idx.mix))
        val.effect <- rep(NA_real_, length(idx.mix))
        val.adj <- rep(NA_real_, length(idx.mix))
        assort <- rep(NA_real_, length(idx.mix))
        cond.med <- rep(NA_real_, length(idx.mix))
        cond.wm <- rep(NA_real_, length(idx.mix))

        for (jj in seq_along(idx.mix)) {
            ii <- idx.mix[jj]

            k0 <- conn$k[ii]
            g.obj <- g.list[[ii]]
            n0 <- nrow(X)

            ## build igraph and prepare weights (affinity if needed)
            el <- adjlist.to.edge.mat(g.obj$adj_list, g.obj$weight_list, n = n0)
            if (nrow(el$edge.mat) == 0L) next

            w.use <- NULL
            if (isTRUE(use.edge.weights)) {
                if (isTRUE(weights.are.edge.lengths)) {
                    w.use <- lengths.to.affinity(el$weights, sigma = sigma.used, method = affinity.method)
                } else {
                    w.use <- el$weights
                    w.use[!is.finite(w.use)] <- 0
                    w.use[w.use < 0] <- 0
                }
            }

            ig <- igraph::make_empty_graph(n = n0, directed = FALSE)
            ig <- igraph::add_edges(ig, as.vector(t(el$edge.mat)))
            igraph::E(ig)$weight <- el$weights

            if (!is.null(w.use) && length(w.use) == nrow(el$edge.mat)) igraph::E(ig)$weight <- w.use

            if (isTRUE(simplify.multiple)) {
                comb.fun <- if (isTRUE(weights.are.edge.lengths)) "max" else "max"
                ## after length->affinity conversion, use max for duplicates
                ig <- igraph::simplify(ig, remove.multiple = TRUE, remove.loops = TRUE,
                                       edge.attr.comb = list(weight = comb.fun, "ignore"))
            }

            ms <- call.mixing.stats(ig, labels.vec = labels.aligned, blocks.vec = blocks.aligned, w.vec = if (!("edge.weights" %in% names(formals(cst.graph.mixing.stats)))) NULL else igraph::E(ig)$weight, seed = seed)

            k.out <- c(k.out, k0)

            val[jj] <- extract.metric(ms, mixing.metric, igraph.obj = ig, labels.vec = labels.aligned)
            val.z[jj] <- extract.metric(ms, "homophily.z", igraph.obj = ig, labels.vec = labels.aligned)
            val.effect[jj] <- extract.metric(ms, "homophily.effect", igraph.obj = ig, labels.vec = labels.aligned)
            val.adj[jj] <- extract.metric(ms, "homophily.adjusted", igraph.obj = ig, labels.vec = labels.aligned)
            assort[jj] <- extract.metric(ms, "assortativity", igraph.obj = ig, labels.vec = labels.aligned)
            cond.med[jj] <- extract.metric(ms, "conductance.median", igraph.obj = ig, labels.vec = labels.aligned)
            cond.wm[jj] <- extract.metric(ms, "conductance.wmean", igraph.obj = ig, labels.vec = labels.aligned)
        }

        mixing.df <- data.frame(
            k = as.integer(k.out),
            metric = as.double(val),
            homophily.z = as.double(val.z),
            homophily.effect = as.double(val.effect),
            homophily.adjusted = as.double(val.adj),
            assortativity = as.double(assort),
            conductance.median = as.double(cond.med),
            conductance.wmean = as.double(cond.wm)
        )
        mixing.df <- mixing.df[order(mixing.df$k), , drop = FALSE]

        ## select k on the chosen mixing.metric curve
        metric.vec <- mixing.df$metric
        k.vec <- mixing.df$k
        idx.ok <- which(is.finite(metric.vec))

        if (length(idx.ok) > 0L) {
            k.opt.mixing <- pick.k.within.eps.global.max(
                metric = metric.vec,
                k.values = k.vec,
                eps = mixing.eps,
                direction = dir0,
                idx.ok = idx.ok,
                require.local.extremum = isTRUE(mixing.require.local.extremum),
                window = mixing.window,
                return.details = FALSE
            )
            k.opt.mixing <- as.integer(k.opt.mixing)
        }
    }

    ## ---- assemble return ----
    out <- list(
        X.graphs = X.graphs,
        k.values = as.integer(k.values),
        connectivity = conn,
        edit = edit.df,
        mixing = mixing.df,
        k.cc.edit = as.integer(k.cc.edit),
        k.opt.edit = as.integer(k.opt.edit),
        k.cc.mixing = as.integer(k.cc.mixing),
        k.opt.mixing = as.integer(k.opt.mixing),
        trim = trim.info,
        params = list(
            method = method,
            kmin = kmin, kmax = kmax,
            pca.dim = pca.dim,
            variance.explained = variance.explained,
            n.cores = n.cores,
            edit.min.lcc.frac = edit.min.lcc.frac,
            edit.eps = edit.eps,
            mixing.metric = mixing.metric,
            mixing.min.lcc.frac = mixing.min.lcc.frac,
            mixing.eps = mixing.eps,
            mixing.require.local.extremum = mixing.require.local.extremum,
            mixing.window = mixing.window,
            n.perm = n.perm,
            use.edge.weights = use.edge.weights,
            weights.are.edge.lengths = weights.are.edge.lengths,
            affinity.method = affinity.method,
            affinity.sigma = sigma.used,
            affinity.eps = affinity.eps,
            simplify.multiple = simplify.multiple
        )
    )
    class(out) <- "build_iknn_graphs_and_selectk"
    out
}

#' Print method for build_iknn_graphs_and_selectk
#'
#' @param x Object from build.iknn.graphs.and.selectk().
#' @param ... Unused.
#' @export
print.build_iknn_graphs_and_selectk <- function(x, ...) {
    if (!inherits(x, "build_iknn_graphs_and_selectk")) stop("x must be class 'build_iknn_graphs_and_selectk'.")
    cat("build.iknn.graphs.and.selectk result\n")
    cat("  k range: ", min(x$k.values), " .. ", max(x$k.values), "\n", sep = "")
    cat("  trimmed: ", isTRUE(x$trim$trimmed), "\n", sep = "")
    if (isTRUE(x$trim$trimmed)) cat("  trim k:  ", x$trim$k.trim, "\n", sep = "")
    cat("  k.cc.edit:    ", x$k.cc.edit, "\n", sep = "")
    cat("  k.opt.edit:   ", x$k.opt.edit, "\n", sep = "")
    cat("  k.cc.mixing:  ", x$k.cc.mixing, "\n", sep = "")
    cat("  k.opt.mixing: ", x$k.opt.mixing, "\n", sep = "")
    invisible(x)
}

#' Plot method for build_iknn_graphs_and_selectk
#'
#' @description
#' Produces diagnostic plots without forwarding \code{...} to multiple panels.
#' Customize panels via \code{connect.args}, \code{edit.args}, \code{mixing.args}.
#'
#' @param x Object from build.iknn.graphs.and.selectk().
#' @param which Character vector selecting panels among: "connect", "edit", "mixing".
#' @param connect.args Named list of arguments forwarded to the connectivity plot only.
#' @param edit.args Named list of arguments forwarded to the edit plot only.
#' @param mixing.args Named list of arguments forwarded to the mixing plot only.
#' @param par.args Named list of arguments forwarded to \code{par()} (e.g., mfrow, mar).
#' @export
plot.build_iknn_graphs_and_selectk <- function(x,
                                              which = c("connect", "edit", "mixing"),
                                              connect.args = list(),
                                              edit.args = list(),
                                              mixing.args = list(),
                                              par.args = list()) {
    if (!inherits(x, "build_iknn_graphs_and_selectk")) stop("x must be class 'build_iknn_graphs_and_selectk'.")
    which <- unique(which)

    ## default layout: stacked panels
    np <- length(which)
    if (np < 1L) return(invisible(NULL))

    par.default <- list(mfrow = c(np, 1), mar = c(3.2, 3.2, 1.5, 0.8), mgp = c(2.0, 0.6, 0), tcl = -0.3)
    par.use <- utils::modifyList(par.default, par.args)

    op <- do.call(graphics::par, par.use)
    on.exit(graphics::par(op), add = TRUE)

    ## ---- connectivity panel ----
    if ("connect" %in% which) {
        df <- x$connectivity
        args <- list(df$k, df$lcc.frac,
                     type = "l", las = 1,
                     xlab = "k", ylab = "LCC fraction",
                     main = "Connectivity (LCC fraction)")
        args <- utils::modifyList(args, connect.args)
        do.call(graphics::plot, args)

        if (is.finite(x$k.cc.edit)) graphics::abline(v = x$k.cc.edit, lty = 2)
        if (is.finite(x$k.cc.mixing)) graphics::abline(v = x$k.cc.mixing, lty = 3)
    }

    ## ---- edit panel ----
    if ("edit" %in% which) {
        if (is.null(x$edit)) {
            graphics::plot.new()
            graphics::title("Edit curve (not computed)")
        } else {
            df <- x$edit
            args <- list(df$k, df$edit.dist.to.next,
                         type = "l", las = 1,
                         xlab = "k", ylab = "edit dist to next",
                         main = "Edit distance between consecutive graphs")
            args <- utils::modifyList(args, edit.args)
            do.call(graphics::plot, args)

            if (is.finite(x$k.cc.edit)) graphics::abline(v = x$k.cc.edit, lty = 2)
            if (is.finite(x$k.opt.edit)) graphics::abline(v = x$k.opt.edit, lty = 3)
        }
    }

    ## ---- mixing panel ----
    if ("mixing" %in% which) {
        if (is.null(x$mixing)) {
            graphics::plot.new()
            graphics::title("Mixing curve (not computed)")
        } else {
            df <- x$mixing
            args <- list(df$k, df$metric,
                         type = "l", las = 1,
                         xlab = "k", ylab = "mixing metric",
                         main = paste0("Mixing metric: ", x$params$mixing.metric))
            args <- utils::modifyList(args, mixing.args)
            do.call(graphics::plot, args)

            if (is.finite(x$k.cc.mixing)) graphics::abline(v = x$k.cc.mixing, lty = 2)
            if (is.finite(x$k.opt.mixing)) graphics::abline(v = x$k.opt.mixing, lty = 3)
        }
    }

    invisible(x)
}
