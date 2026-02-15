#' Adaptive Graph Geodesic Model-Averaged Local Linear Regression
#'
#' @description Performs geodesic model-averaged local linear regression on
#'     graph-structured data using an adaptive uniform grid approach.
#'
#' @details
#' The AGEMALO algorithm proceeds through several steps:
#' 1. Creates a maximal packing of the input graph
#' 2. For each packing vertex it
#' - creates a hierarchy of local geodesics passing through that vertex
#' - computes minimal and maximal bandwidths and bandwidth candidates
#' - fits local weighted linear models along geodesic paths
#' 3. Combines predictions using model averaging
#'
#' The function supports optional bootstrap confidence intervals and permutation testing
#' for statistical inference.
#'
#' @param adj.list List of integer vectors. Each vector contains indices of vertices
#'        adjacent to the corresponding vertex. Indices should be 1-based.
#' @param weight.list List of numeric vectors. Each vector contains weights of edges
#'        corresponding to adjacencies in adj.list.
#' @param y Numeric vector of response values at each vertex.
#' @param min.path.size Integer. Minimum number of vertices required in valid paths.
#' @param n.packing.vertices Integer. Number of vertices to use in the maximal packing.
#'        Defaults to the number of vertices (length of y).
#' @param max.packing.iterations Integer. Maximum number of iterations for the packing algorithm.
#' @param packing.precision Numeric. Precision threshold for packing convergence.
#' @param n.bws Integer. Number of candidate bandwidths to evaluate.
#' @param log.grid Logical. Whether to use logarithmic spacing for bandwidth candidates.
#' @param min.bw.factor Numeric. Factor multiplied by graph diameter for minimum bandwidth.
#' @param max.bw.factor Numeric. Factor multiplied by graph diameter for maximum bandwidth.
#' @param dist.normalization.factor Numeric. Factor for normalizing graph distances.
#' @param kernel.type Integer. Type of kernel function (1-7, with 7 as default).
#' @param model.tolerance Numeric. Convergence tolerance for model fitting.
#' @param model.blending.coef Numeric. Blending coefficient for model averaging.
#' @param n.bb Integer. Number of bootstrap iterations (0 for no bootstrap).
#' @param cri.probability Numeric. Confidence level for bootstrap intervals (0-1).
#' @param n.perms Integer. Number of permutation test iterations (0 for no testing).
#' @param verbose Logical. Whether to print progress information.
#'
#' @return A list containing:
#'   \item{graph.diameter}{Numeric. Computed diameter of input graph.}
#'   \item{grid.opt.bw}{Numeric vector. Optimal bandwidth for each grid vertex.}
#'   \item{predictions}{Numeric vector. Model-averaged predictions for original vertices.}
#'   \item{grid.predictions}{Numeric vector. Model-averaged predictions for grid vertices.}
#'   \item{bb.predictions}{Matrix. Bootstrap predictions (if n.bb > 0).}
#'   \item{cri.lower}{Numeric vector. Lower confidence bounds (if n.bb > 0).}
#'   \item{cri.upper}{Numeric vector. Upper confidence bounds (if n.bb > 0).}
#'   \item{null.predictions}{Matrix. Permutation test predictions (if n.perms > 0).}
#'   \item{p.values}{Numeric vector. Vertex-wise p-values (if n.perms > 0).}
#'   \item{effect.sizes}{Numeric vector. Effect sizes (if n.perms > 0).}
#'   \item{significant.vertices}{Logical vector. Significance indicators (if n.perms > 0).}
#'
#' @examples
#' \dontrun{
#' # Create a simple graph with 3 vertices
#' adj.list <- list(c(2), c(1, 3), c(2))
#' weight.list <- list(c(1), c(1, 1), c(1))
#' y <- c(1, 2, 3)
#'
#' # Run basic analysis
#' result <- agemalo(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   min.path.size = 2,
#'   n.packing.vertices = 5,
#'   n.bws = 10
#' )
#'
#' # Run with bootstrap confidence intervals
#' result.boot <- agemalo(
#'   adj.list = adj.list,
#'   weight.list = weight.list,
#'   y = y,
#'   min.path.size = 2,
#'   n.packing.vertices = 5,
#'   n.bws = 10,
#'   n.bb = 100
#' )
#' }
#'
#' @export
agemalo <- function(adj.list,
                    weight.list,
                    y,
                    min.path.size = 6,
                    n.packing.vertices = length(y),
                    max.packing.iterations = 20,
                    packing.precision = 0.0001,
                    n.bws = 50,
                    log.grid = TRUE,
                    min.bw.factor = 0.025,
                    max.bw.factor = 0.99,
                    dist.normalization.factor = 1.1,
                    kernel.type = 7L,
                    model.tolerance = 1e-6,
                    model.blending.coef = 0.1,
                    n.bb = 0L,
                    cri.probability = 0.95,
                    n.perms = 0L,
                    verbose = FALSE) {

    # Input validation
    if (!is.list(adj.list) || !is.list(weight.list))
        stop("adj.list and weight.list must be lists")

    if (length(adj.list) != length(weight.list))
        stop("adj.list and weight.list must have the same length")

    if (length(y) != length(adj.list))
        stop("y must have the same length as adj.list")

    if (!is.numeric(y))
        stop("y must be numeric")

    if (min.path.size < 2)
        stop("min.path.size must be at least 2")

    if (n.packing.vertices < 1)
        stop("n.packing.vertices must be positive")

    if (!is.numeric(max.packing.iterations) || length(max.packing.iterations) != 1 ||
        max.packing.iterations != floor(max.packing.iterations) ||
        max.packing.iterations < 1) {
        stop("max.packing.iterations must be an integer greater than 0")
    }

    if (n.bws < 1)
        stop("n.bws must be positive")

    if (!is.numeric(min.bw.factor))
        stop("min.bw.factor must be numeri")

    if (!is.numeric(max.bw.factor) || max.bw.factor <= 0)
        stop("max.bw.factor must be positive")

    if (min.bw.factor >= max.bw.factor)
        stop("max.bw.factor must be greater than min.bw.factor")

    if (packing.precision <= 0 || packing.precision > 0.1)
        stop("packing.precision must be positive not greater than 0.1")

    if (dist.normalization.factor < 1.1)
        stop("dist.normalization.factor must be greater than or equal to 1.1")

    if (!is.numeric(kernel.type) || kernel.type < 0 || kernel.type > 7)
        stop("kernel.type must be an integer between 0 and 7")

    if (model.tolerance <= 0)
        stop("model.tolerance must be positive")

    if (n.bb < 0)
        stop("n.bb must be non-negative")

    if (cri.probability <= 0 || cri.probability >= 1)
        stop("cri.probability must be between 0 and 1")

    if (n.perms < 0)
        stop("n.perms must be non-negative")

    if (model.blending.coef < 0 || model.blending.coef > 1) {
        stop("model.blending.coef must be between 0 and 1")
    }

    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    result <- .Call("S_agemalo",
                   adj.list.0based,
                   weight.list,
                   as.double(y),
                   as.integer(min.path.size),
                   as.integer(n.packing.vertices),
                   as.integer(max.packing.iterations),
                   as.double(packing.precision),
                   as.integer(n.bws),
                   as.logical(log.grid),
                   as.double(min.bw.factor),
                   as.double(max.bw.factor),
                   as.double(dist.normalization.factor),
                   as.integer(kernel.type),
                   as.double(model.tolerance),
                   as.double(model.blending.coef),
                   as.integer(n.bb),
                   as.double(cri.probability),
                   as.integer(n.perms),
                   as.logical(verbose))

    result
}

#' Test Bidirectional Dijkstra and Composite Path Validation
#'
#' @description
#' Validates the correctness of bidirectional Dijkstra algorithm and
#' composite path geodesic verification using igraph.
#'
#' @param grid.vertex Numeric index of the grid vertex to validate
#' @param debug.dir Character string specifying the directory containing debugging data
#' @param full.report Logical; if TRUE, returns detailed validation information
#' @param graph List; An output of graph_from_data_frame()
#'
#' @return A list containing test results and validation summaries
#'
#' @importFrom igraph graph_from_data_frame distances shortest_paths
#' @importFrom utils read.csv
#'
#' @export
test.bidirectional.dijkstra <- function(grid.vertex, debug.dir, full.report = FALSE, graph = NULL) {

    if (is.null(graph)) {
        ## Load the full graph
        graph.edges <- read.csv(file.path(debug.dir, "full_graph.csv"), header = TRUE)
        graph <- graph_from_data_frame(graph.edges, directed = FALSE)
    }

    ## Base path for files
    base_path <- file.path(debug.dir, paste0("grid.vertex_", grid.vertex))

    ## 1. Test bidirectional Dijkstra algorithm
    bidirectional_results <- list(passed = TRUE)
    if (file.exists(paste0(base_path, "_bidirectional.csv"))) {
        bidi_data <- read.csv(paste0(base_path, "_bidirectional.csv"))

        ## Calculate ground truth distances using igraph
        bidi_data$igraph_length <- 0
        for (i in 1:nrow(bidi_data)) {
            source <- bidi_data$source[i]
            target <- bidi_data$target[i]
            bidi_data$igraph_length[i] <- distances(graph, v = source, to = target)[1,1]
        }

        ## Compare results
        bidi_data$difference <- abs(bidi_data$bidirectional_length - bidi_data$igraph_length)
        bidi_data$match <- bidi_data$difference < 1e-10

        bidirectional_results$total_tests <- nrow(bidi_data)
        bidirectional_results$passed_tests <- sum(bidi_data$match)
        bidirectional_results$passed <- all(bidi_data$match)

        if (full.report) {
            bidirectional_results$details <- bidi_data
            if (!bidirectional_results$passed) {
                bidirectional_results$failures <- bidi_data[!bidi_data$match,]
            }
        }
    }

    ## 2. Test composite path validation
    composite_results <- list(passed = TRUE)
    if (file.exists(paste0(base_path, "_composite.csv"))) {
        comp_data <- read.csv(paste0(base_path, "_composite.csv"))

        ## For each composite path, check if it's truly a geodesic
        comp_data$igraph_length <- 0
        for (i in 1:nrow(comp_data)) {
            start_i <- comp_data$start_i[i]
            end_j <- comp_data$start_j[i]  ## The composite path goes from start_i to start_j

            ## Get the shortest path length between endpoints
            comp_data$igraph_length[i] <- distances(graph, v = start_i, to = end_j)[1,1]
        }

        ## Compare with composite path length
        comp_data$is_true_geodesic <- abs(comp_data$composite_length - comp_data$igraph_length) < 1e-10
        comp_data$validation_correct <- (comp_data$is_geodesic == comp_data$is_true_geodesic)

        composite_results$total_tests <- nrow(comp_data)
        composite_results$passed_tests <- sum(comp_data$validation_correct)
        composite_results$passed <- all(comp_data$validation_correct)

        if (full.report) {
            composite_results$details <- comp_data
            if (!composite_results$passed) {
                composite_results$false_positives <- comp_data[comp_data$is_geodesic & !comp_data$is_true_geodesic,]
                composite_results$false_negatives <- comp_data[!comp_data$is_geodesic & comp_data$is_true_geodesic,]
            }
        }
    }

    ## Return combined results
    return(list(
        bidirectional_dijkstra = bidirectional_results,
        composite_path_validation = composite_results,
        passed = bidirectional_results$passed && composite_results$passed
    ))
}
