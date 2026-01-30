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

#' Select \code{k} within a relative tolerance of the global minimum edit distance
#'
#' @description
#' Selects the smallest \code{k} in a connected-tail regime whose edit distance is
#' within \code{(1 + eps)} times the minimum edit distance over that regime.
#'
#' The inputs \code{edit.distances} and \code{k.for.edit} are assumed to represent
#' edit distances between consecutive graphs \eqn{(G_k, G_{k+1})}, with
#' \code{k.for.edit[i]} giving the left endpoint \eqn{k} associated with
#' \code{edit.distances[i]}.
#'
#' @param edit.distances Numeric vector of edit distances between consecutive graphs.
#' @param k.for.edit Integer vector of the same length as \code{edit.distances},
#'   giving the \eqn{k} values associated with each edit distance.
#' @param k.cc Integer scalar. Start of the terminal connected tail; only values with
#'   \code{k.for.edit >= k.cc} are eligible.
#' @param eps Nonnegative numeric scalar. Relative tolerance; candidates satisfy
#'   \code{d <= (1 + eps) * d.min} where \code{d.min} is the minimum over the eligible
#'   regime.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{k.pick}: The selected \code{k} value.
#'   \item \code{d.pick}: The edit distance at \code{k.pick}.
#'   \item \code{d.min}: The minimum edit distance over eligible \code{k}.
#'   \item \code{d.thld}: The threshold \code{(1 + eps) * d.min}.
#'   \item \code{eps}: The input tolerance \code{eps}.
#'   \item \code{candidates}: A data frame of eligible candidate \code{k} values and
#'   their edit distances.
#' }
#'
#' @keywords internal
#' @noRd
pick.k.within.eps.global.min <- function(edit.distances, k.for.edit, k.cc, eps = 0.05) {

    if (!is.numeric(k.cc) || length(k.cc) != 1 || is.na(k.cc)) {
        stop("k.cc must be a single non-NA numeric")
    }
    if (k.cc != as.integer(k.cc)) stop("k.cc must be integer-valued.")
    k.cc <- as.integer(k.cc)

    if (!is.numeric(edit.distances) || length(edit.distances) < 1) {
        stop("edit.distances must be a non-empty numeric vector")
    }

    if (!is.numeric(k.for.edit) || length(k.for.edit) != length(edit.distances)) {
        stop("k.for.edit must be numeric/integer and same length as edit.distances")
    }
    if (!all(k.for.edit == as.integer(k.for.edit))) {
        stop("k.for.edit must be integer-valued")
    }
    k.for.edit <- as.integer(k.for.edit)
    if (anyNA(k.for.edit)) stop("k.for.edit contains NA after coercion to integer.")
    if (is.unsorted(k.for.edit, strictly = TRUE)) stop("k.for.edit must be strictly increasing.")

    if (!is.numeric(eps) || length(eps) != 1 || eps < 0) {
        stop("eps must be a single nonnegative numeric value")
    }

    keep <- (k.for.edit >= k.cc)
    if (!any(keep)) {
        stop("No k values satisfy k >= k.cc")
    }

    d.keep <- edit.distances[keep]
    k.keep <- k.for.edit[keep]

    ## Must check NA/finite before min/thresholding
    if (all(is.na(d.keep))) {
        stop("All edit distances are NA in the k >= k.cc regime.")
    }
    if (!any(is.finite(d.keep))) {
        stop("No finite edit distances in the k >= k.cc regime.")
    }

    d.min <- min(d.keep, na.rm = TRUE)
    d.thld <- (1 + eps) * d.min

    cand.idx <- which(is.finite(d.keep) & d.keep <= d.thld)
    if (length(cand.idx) == 0) {
        stop("No candidates found. Check eps / inputs.")
    }

    k.pick <- as.integer(min(k.keep[cand.idx]))
    d.pick <- d.keep[match(k.pick, k.keep)]
    if (length(d.pick) != 1 || !is.finite(d.pick)) {
        stop("Internal error computing d.pick")
    }

    list(
        k.pick = k.pick,
        d.pick = as.numeric(d.pick),
        d.min = as.numeric(d.min),
        d.thld = as.numeric(d.thld),
        eps = eps,
        candidates = data.frame(
            k = k.keep[cand.idx],
            edit.distance = d.keep[cand.idx]
        )
    )
}

#' Build ikNN graph sequence and select \code{k} by edit-distance stability
#'
#' @description
#' Builds an intersection k-nearest-neighbor (ikNN) graph sequence for
#' \code{k = kmin:kmax} using \code{create.iknn.graphs()}, then selects an
#' optimal \code{k} based on stability of graph structure across the terminal
#' connected tail.
#'
#' The procedure:
#' \enumerate{
#'   \item Construct a sequence of geometrically pruned ikNN graphs for \code{kmin:kmax}.
#'   \item If all graphs in the range are disconnected, trim \code{X} to the largest
#'         connected component of the graph at \code{k = kmin} and rebuild the graph sequence.
#'   \item Define \code{k.cc} as the smallest \code{k} such that all graphs for \code{k' >= k}
#'         are connected (terminal connected tail). Stop if no such tail exists in the range.
#'   \item Compute edit distances between consecutive graphs in the tail regime and select
#'         \code{k.opt} as the smallest \code{k >= k.cc} within a relative tolerance of the
#'         minimum edit distance (via \code{pick.k.within.eps.global.min()}).
#' }
#'
#' @param X Numeric matrix with rows corresponding to vertices/observations.
#' @param kmin Integer \eqn{\ge 1}. Minimum k for the ikNN graph sequence.
#' @param kmax Integer \eqn{\ge kmin}. Maximum k for the ikNN graph sequence.
#'   Requires \code{nrow(X) > kmax}. If trimming is applied, this requirement must also
#'   hold after trimming.
#' @param pruning A named list of pruning parameters passed to \code{create.iknn.graphs()}:
#'   \itemize{
#'     \item \code{max.path.edge.ratio.deviation.thld}
#'     \item \code{path.edge.ratio.percentile}
#'     \item \code{threshold.percentile}
#'   }
#' @param n.cores Integer scalar or \code{NULL}. Number of cores passed to
#'   \code{create.iknn.graphs()}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{k.opt}: Selected k value.
#'   \item \code{X.graphs}: The object returned by \code{create.iknn.graphs()}.
#'   \item \code{X.graphs.stats}: Summary statistics from \code{summary(X.graphs)}.
#'   \item \code{edit.distances}: Numeric vector of edit distances between consecutive
#'   graphs in the regime used for selection.
#'   \item \code{k.for.edit}: Integer vector of k values corresponding to \code{edit.distances}.
#'   \item \code{pick.k.within.eps.global.min.res}: The full result returned by
#'   \code{pick.k.within.eps.global.min()}.
#'   \item \code{graph.opt}: The selected graph object at \code{k = k.opt} from the
#'   geometrically pruned graph sequence.
#' }
#'
#' @export
build.iknn.graphs.and.selectk <- function(X,
                                         kmin,
                                         kmax,
                                         pruning = list(
                                             max.path.edge.ratio.deviation.thld = 0.1,
                                             path.edge.ratio.percentile = 0.5,
                                             threshold.percentile = 0
                                         ),
                                         n.cores = NULL,
                                         verbose = TRUE) {

    stopifnot(is.matrix(X))
    stopifnot(nrow(X) > kmax)
    stopifnot(kmin >= 1, kmax >= kmin)

    stopifnot(kmin == as.integer(kmin), kmax == as.integer(kmax))
    kmin <- as.integer(kmin); kmax <- as.integer(kmax)

    k.values <- kmin:kmax

    X.graphs <- create.iknn.graphs(
        X,
        kmin = kmin,
        kmax = kmax,
        max.path.edge.ratio.deviation.thld = pruning$max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = pruning$path.edge.ratio.percentile,
        threshold.percentile = pruning$threshold.percentile,
        compute.full = TRUE,
        pca.dim = NULL,
        variance.explained = NULL,
        n.cores = n.cores,
        verbose = verbose
    )

    X.graphs.stats <- summary(X.graphs)

    any.connected <- any(X.graphs.stats$n_ccomp == 1)

    if (!any.connected) {
        if (verbose) {
            cat("\nAll graphs in k range have >1 connected component.\n")
            cat("Applying outlier trimming (largest CC) using graph at k = kmin.\n")
        }

        g0 <- X.graphs$geom_pruned_graphs[[1]]

        if (is.null(g0$adj_list)) stop("geom_pruned_graphs[[1]] is missing adj_list.")

        trim.res <- trim.X.to.main.cc(
            X = X,
            adj.list = g0$adj_list,
            verbose = verbose
        )
        X <- trim.res$X

        if (nrow(X) <= kmax) {
            stop("After trimming to main CC, nrow(X) <= kmax. Reduce kmax or use a less aggressive trimming/pruning.")
        }

        if (verbose) cat("\nRebuilding ikNN graph sequence (after trimming)\n")

        X.graphs <- create.iknn.graphs(
            X,
            kmin = kmin,
            kmax = kmax,
            max.path.edge.ratio.deviation.thld = pruning$max.path.edge.ratio.deviation.thld,
            path.edge.ratio.percentile = pruning$path.edge.ratio.percentile,
            threshold.percentile = pruning$threshold.percentile,
            compute.full = TRUE,
            n.cores = n.cores,
            verbose = verbose
        )
        X.graphs.stats <- summary(X.graphs)
    }

    graphs <- X.graphs$geom_pruned_graphs
    if (is.null(graphs)) stop("create.iknn.graphs() did not return geom_pruned_graphs (compute.full=TRUE required).")

    ## Define k.cc as the start of the terminal connected tail
    ## i.e., the smallest k such that all graphs for k' >= k are connected
    n.ccomp.vec <- X.graphs.stats$n_ccomp

    if (length(n.ccomp.vec) != length(k.values)) {
        stop("Length mismatch: X.graphs.stats$n_ccomp must match kmin:kmax")
    }

    disc.idx <- which(n.ccomp.vec > 1)

    if (length(disc.idx) == 0) {
        ## All graphs are connected across the full k range
        k.cc <- as.integer(kmin)
    } else if (max(disc.idx) == length(k.values)) {
        ## Last k is disconnected -> no terminal connected tail in this range
        k.cc <- NA_integer_
    } else {
        ## First k after the last disconnected k
        k.cc <- as.integer(k.values[max(disc.idx) + 1L])
    }

    ## Enforce that the tail is indeed connected
    if (!is.na(k.cc)) {
        tail.idx <- which(k.values >= k.cc)
        if (any(n.ccomp.vec[tail.idx] > 1)) {
            stop("Internal error: computed k.cc does not define a connected tail")
        }
    } else {
        stop("No terminal connected tail in this k range. Increase kmax and/or relax pruning.")
    }

    if (verbose && k.cc > kmin)  {
        cat("\nFound connected graphs starting at k =", k.cc, "\n")
        cat("Restricting stability analysis to k in [", k.cc, ",", kmax, "]\n", sep = "")
    }

    keep.idx <- NULL
    if (!is.na(k.cc) && k.cc > kmin) {
        keep.idx <- which(k.values >= k.cc)
        graphs.use <- X.graphs$geom_pruned_graphs[keep.idx]
        k.values.use <- k.values[keep.idx]
    } else {
        graphs.use <- X.graphs$geom_pruned_graphs
        k.values.use <- k.values
    }

    if (length(graphs.use) < 2L) {
        stop("Need at least two graphs to compute edit distances; increase k range or adjust k.cc.")
    }

    if (verbose) {
        cat("\nComputing edit distance between consecutive graphs\n")
    }

    edit.distances <- internal.compute.edit.distances(graphs.use)

    if (length(edit.distances) != (length(graphs.use) - 1L)) {
        stop("internal.compute.edit.distances() returned unexpected length.")
    }

    ## k values corresponding to edit.distances (pairs G_k vs G_{k+1})
    k.for.edit <- k.values.use[-length(k.values.use)]

    res.10 <- pick.k.within.eps.global.min(edit.distances, k.for.edit, k.cc = k.cc, eps = 0.10)

    if (verbose) {
        cat("\nSelected k.opt =", res.10$k.pick,
            " (eps=0.10, d.pick=", res.10$d.pick,
            ", d.min=", res.10$d.min, ")\n", sep = "")
    }

    k.opt <- res.10$k.pick

    ## extract chosen graph
    idx <- which(k.values == k.opt)
    if (length(idx) != 1L) stop("Internal error: k.opt not found uniquely in kmin:kmax.")
    g.opt <- graphs[[idx]]

    if (verbose) {
        cat(sprintf("DONE: k.opt=%d\n", k.opt))
    }

    out <- list(
        k.opt = k.opt,
        X.graphs = X.graphs,
        X.graphs.stats = X.graphs.stats,
        edit.distances = edit.distances,
        k.for.edit = k.for.edit,
        pick.k.within.eps.global.min.res = res.10,
        graph.opt = g.opt
    )

    class(out) <- "build_iknn_graphs_and_selectk"

    out
}

#' Plot method for \code{build_iknn_graphs_and_selectk} objects
#'
#' @description
#' Produces a diagnostic stability plot for objects returned by
#' \code{build.iknn.graphs.and.selectk()} with class
#' \code{"build_iknn_graphs_and_selectk"}.
#'
#' The plot shows the edit distance between consecutive graphs in the
#' selected regime (typically the terminal connected tail), as a function of
#' \code{k} (the left endpoint of the consecutive pair \eqn{(G_k, G_{k+1})}).
#' Vertical reference lines indicate \code{k.cc} (start of the connected tail)
#' and \code{k.opt} (selected k).
#'
#' @param x An object of class \code{"build_iknn_graphs_and_selectk"}.
#'   Must contain components \code{k.for.edit} and \code{edit.distances}. If present,
#'   \code{k.cc} and \code{k.opt} are used for reference lines.
#' @param ... Additional arguments passed to \code{plot()}.
#' @param type Plot type passed to \code{plot()}. Default \code{"b"}.
#' @param pch Plotting character passed to \code{plot()}. Default \code{16}.
#' @param las Axis label style passed to \code{plot()}. Default \code{1}.
#' @param xlab X-axis label. Default \code{"k (for consecutive pair)"}.
#' @param ylab Y-axis label. Default \code{"Edit distance between consecutive graphs"}.
#' @param main Plot title. Default \code{"ikNN graph stability: edit distance curve"}.
#' @param add.grid Logical. If \code{TRUE}, adds a \code{grid()}.
#' @param lty.k.cc Line type for \code{k.cc} reference line. Default \code{2}.
#' @param lty.k.opt Line type for \code{k.opt} reference line. Default \code{3}.
#' @param legend.pos Legend position passed to \code{legend()}. Default \code{"topright"}.
#' @param legend.bty Legend box type passed to \code{legend()}. Default \code{"n"}.
#' @param show.threshold A horizontal line corresponding to \eqn{(1+\varepsilon)d_{\min}} value.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
plot.build_iknn_graphs_and_selectk <- function(
    x,
    ...,
    type = "b",
    pch = 16,
    las = 1,
    xlab = "k (for consecutive pair)",
    ylab = "Edit distance between consecutive graphs",
    main = "ikNN graph stability: edit distance curve",
    add.grid = TRUE,
    lty.k.cc = 2,
    lty.k.opt = 3,
    legend.pos = "topright",
    legend.bty = "n",
    show.threshold = FALSE
) {

    if (is.null(x) || !inherits(x, "build_iknn_graphs_and_selectk")) {
        stop("x must be an object of class 'build_iknn_graphs_and_selectk'.")
    }

    if (is.null(x$k.for.edit) || is.null(x$edit.distances)) {
        stop("x must contain components 'k.for.edit' and 'edit.distances'.")
    }

    k.for.edit <- x$k.for.edit
    edit.distances <- x$edit.distances

    if (!is.numeric(k.for.edit) || !is.numeric(edit.distances)) {
        stop("'k.for.edit' and 'edit.distances' must be numeric.")
    }

    if (length(k.for.edit) != length(edit.distances)) {
        stop("Length mismatch: length(k.for.edit) must equal length(edit.distances).")
    }

    ## Basic plot
    graphics::plot(
        k.for.edit,
        edit.distances,
        type = type,
        pch = pch,
        las = las,
        xlab = xlab,
        ylab = ylab,
        main = main,
        ...
    )

    if (isTRUE(add.grid)) {
        graphics::grid()
    }

    ## Optional: horizontal threshold line at (1+eps)*d.min from stored picker result
    if (isTRUE(show.threshold)) {
        if (!is.null(x$pick.k.within.eps.global.min.res) &&
            is.list(x$pick.k.within.eps.global.min.res) &&
            !is.null(x$pick.k.within.eps.global.min.res$d.thld) &&
            length(x$pick.k.within.eps.global.min.res$d.thld) == 1L &&
            is.finite(x$pick.k.within.eps.global.min.res$d.thld)) {

            d.thld <- x$pick.k.within.eps.global.min.res$d.thld
            graphics::abline(h = d.thld, lty = 4)

        } else {
            warning("show.threshold=TRUE, but x$pick.k.within.eps.global.min.res$d.thld is missing or invalid.")
        }
    }

    ## Optional reference lines
    has.k.cc <- !is.null(x$k.cc) && is.finite(x$k.cc) && length(x$k.cc) == 1
    has.k.opt <- !is.null(x$k.opt) && is.finite(x$k.opt) && length(x$k.opt) == 1

    if (has.k.cc) {
        graphics::abline(v = x$k.cc, lty = lty.k.cc)
    }
    if (has.k.opt) {
        graphics::abline(v = x$k.opt, lty = lty.k.opt)
    }

    ## Legend
    legend.items <- character(0)
    if (has.k.cc) legend.items <- c(legend.items, paste0("k.cc = ", x$k.cc))
    if (has.k.opt) legend.items <- c(legend.items, paste0("k.opt = ", x$k.opt))

    if (isTRUE(show.threshold) &&
        !is.null(x$pick.k.within.eps.global.min.res$d.thld) &&
        length(x$pick.k.within.eps.global.min.res$d.thld) == 1L &&
        is.finite(x$pick.k.within.eps.global.min.res$d.thld)) {
        legend.items <- c(legend.items, paste0("d.thld = ", signif(x$pick.k.within.eps.global.min.res$d.thld, 6)))
    }

    if (length(legend.items) > 0) {
        graphics::legend(
            legend.pos,
            legend = legend.items,
            bty = legend.bty
        )
    }

    invisible(x)
}

#' Print method for \code{build_iknn_graphs_and_selectk} objects
#'
#' @description
#' Prints a compact textual summary for objects returned by
#' \code{build.iknn.graphs.and.selectk()} with class
#' \code{"build_iknn_graphs_and_selectk"}. This method does not produce any plots.
#'
#' @param x An object of class \code{"build_iknn_graphs_and_selectk"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.build_iknn_graphs_and_selectk <- function(x, ...) {

    if (is.null(x) || !inherits(x, "build_iknn_graphs_and_selectk")) {
        stop("x must be an object of class 'build_iknn_graphs_and_selectk'.")
    }

    ## Safe getters
    .get1 <- function(obj, nm) {
        v <- obj[[nm]]
        if (is.null(v)) return(NULL)
        if (length(v) == 0) return(NULL)
        v
    }

    k.opt <- .get1(x, "k.opt")
    k.cc <- .get1(x, "k.cc")
    k.for.edit <- .get1(x, "k.for.edit")
    edit.distances <- .get1(x, "edit.distances")

    ## Header
    cat("build_iknn_graphs_and_selectk\n")

    ## k summary
    if (!is.null(k.cc)) cat("  k.cc :", k.cc, "\n")
    if (!is.null(k.opt)) cat("  k.opt:", k.opt, "\n")

    ## Edit distance summary
    if (!is.null(k.for.edit) && !is.null(edit.distances) &&
        is.numeric(k.for.edit) && is.numeric(edit.distances) &&
        length(k.for.edit) == length(edit.distances) && length(edit.distances) > 0) {

        cat("  edit distances (n =", length(edit.distances), ")\n")
        cat("    k range:", min(k.for.edit), "to", max(k.for.edit), "\n")

        fin <- is.finite(edit.distances)
        if (any(fin)) {
            d.min <- min(edit.distances[fin])
            d.med <- stats::median(edit.distances[fin])
            d.max <- max(edit.distances[fin])

            cat("    min/median/max:", signif(d.min, 6), "/",
                signif(d.med, 6), "/", signif(d.max, 6), "\n", sep = "")
        } else {
            cat("    all distances are non-finite\n")
        }
    } else {
        cat("  edit distances: <not available>\n")
    }

    ## Picker summary if present
    pick.res <- .get1(x, "pick.k.within.eps.global.min.res")
    if (!is.null(pick.res) && is.list(pick.res)) {
        eps <- pick.res$eps
        d.thld <- pick.res$d.thld
        d.min <- pick.res$d.min
        if (!is.null(eps) || !is.null(d.min) || !is.null(d.thld)) {
            cat("  selection rule\n")
            if (!is.null(eps)) cat("    eps  :", eps, "\n")
            if (!is.null(d.min)) cat("    d.min:", signif(d.min, 6), "\n")
            if (!is.null(d.thld)) cat("    d.thld:", signif(d.thld, 6), "\n")
        }
    }

    ## Graph sequence summary if present
    stats <- .get1(x, "X.graphs.stats")
    if (!is.null(stats) && is.data.frame(stats)) {
        cat("  graph sequence stats: ", nrow(stats), " graphs\n", sep = "")
        if (all(c("k", "n_ccomp", "edges") %in% names(stats))) {
            ## Show first/last rows compactly
            cat("    k:", stats$k[1], "to", stats$k[nrow(stats)], "\n")
            cat("    n_ccomp (min/max):", min(stats$n_ccomp), "/", max(stats$n_ccomp), "\n")
        }
    }

    invisible(x)
}
