#' Create Graph Flow Complex with Harmonic Extension
#'
#' @description
#' Computes a graph flow complex by identifying local extrema and their neighborhoods,
#' then applies harmonic extension to smooth spurious extrema while preserving
#' significant topological features. The algorithm uses hop index thresholding to
#' identify and smooth spurious extrema through various smoothing methods.
#'
#' @param adj.list A list where each element i contains the indices of vertices
#'   adjacent to vertex i. Indices should be 1-based (R convention). Must have
#'   the same length as \code{weight.list} and \code{y}.
#' @param weight.list A list where each element i contains the weights of edges
#'   from vertex i to its neighbors. Must have the same structure and length as
#'   \code{adj.list}.
#' @param y A numeric vector of function values defined on the graph vertices.
#'   Must have the same length as \code{adj.list}.
#' @param hop.idx.thld Numeric scalar specifying the hop index threshold. Extrema
#'   with hop index less than or equal to this value are considered spurious and
#'   will be smoothed. Default is 5. Must be non-negative.
#' @param smoother.type Integer specifying the smoothing algorithm:
#'   \itemize{
#'     \item 0: Weighted Mean (default) - Simple weighted averaging
#'     \item 1: Harmonic Iterative - Iterative harmonic smoothing
#'     \item 2: Harmonic Eigen - Eigenfunction-based harmonic smoothing
#'     \item 3: Hybrid Biharmonic-Harmonic - Combined biharmonic and harmonic
#'     \item 4: Boundary Smoothed Harmonic - Harmonic with boundary smoothing
#'   }
#' @param max.outer.iterations Maximum number of outer iterations for the smoothing
#'   process. Default is 5.
#' @param max.inner.iterations Maximum number of iterations for each individual
#'   smoothing operation. Default is 100.
#' @param smoothing.tolerance Numeric scalar for convergence tolerance in smoothing
#'   algorithms. Default is 1e-6.
#' @param sigma Numeric scalar controlling the width of the smoothing kernel.
#'   Default is 1.0. Larger values produce more smoothing.
#' @param process.in.order Logical; if TRUE (default), processes extrema in
#'   ascending order of hop index, smoothing the most spurious extrema first.
#' @param verbose Logical; if TRUE (default), prints progress information during
#'   the computation.
#' @param detailed.recording Logical; if TRUE, records detailed information about
#'   each smoothing step for later visualization. Default is FALSE. Note that this
#'   increases memory usage.
#'
#' @return An object of class \code{"gflow_cx"} containing:
#'   \describe{
#'     \item{harmonic_predictions}{Numeric vector of smoothed function values at
#'       each vertex after applying harmonic extension.}
#'     \item{lmin_hop_nbhds}{List of hop neighborhoods for local minima. Each
#'       element contains vertex index, hop index, and neighborhood information.}
#'     \item{lmax_hop_nbhds}{List of hop neighborhoods for local maxima. Each
#'       element contains vertex index, hop index, and neighborhood information.}
#'     \item{extrema_df}{Data frame summarizing all extrema with columns: vertex,
#'       hop_idx, is_max, label, and fn_value.}
#'     \item{extrema_df2}{Alternative data frame format with spurious extrema
#'       information included.}
#'     \item{smoothing_history}{(If detailed.recording=TRUE) List of records
#'       detailing each smoothing step for visualization.}
#'   }
#'   The object also has attributes:
#'   \itemize{
#'     \item \code{smoother.type}: The integer smoother type used
#'     \item \code{smoother.name}: Human-readable name of the smoother
#'     \item \code{hop.idx.threshold}: The hop index threshold used
#'   }
#'
#' @details
#' The algorithm proceeds in several steps:
#' \enumerate{
#'   \item Identifies local minima and maxima in the function values on the graph
#'   \item Computes hop neighborhoods around each extremum
#'   \item Calculates hop indices measuring the significance of each extremum
#'   \item Identifies spurious extrema (those with hop index <= threshold)
#'   \item Applies harmonic extension to smooth out spurious extrema
#'   \item Returns the smoothed function values and extrema information
#' }
#'
#' The hop index measures how many graph hops are needed before the function
#' value exceeds the extremum value, providing a scale-free measure of extremum
#' significance. Spurious extrema typically have small hop indices and represent
#' noise or insignificant features.
#'
#' @examples
#' \dontrun{
#' # Create a simple triangle graph
#' adj.list <- list(c(2,3), c(1,3), c(1,2))
#' weight.list <- list(c(1,1), c(1,1), c(1,1))
#'
#' # Function with a spurious maximum at vertex 2
#' y <- c(0.5, 1.0, 0.7)
#'
#' # Smooth out spurious extrema
#' result <- create.gflow.cx(
#'   adj.list, weight.list, y,
#'   hop.idx.thld = 1,
#'   smoother.type = 0,  # Weighted Mean
#'   verbose = FALSE
#' )
#'
#' # Check smoothed values
#' print(result$harmonic_predictions)
#'
#' # More complex example with detailed recording
#' # Create a larger graph (grid)
#' n <- 5
#' adj.list <- vector("list", n*n)
#' weight.list <- vector("list", n*n)
#'
#' for(i in 1:n) {
#'   for(j in 1:n) {
#'     idx <- (i-1)*n + j
#'     neighbors <- c()
#'     weights <- c()
#'
#'     # Add edges to grid neighbors
#'     if(i > 1) { neighbors <- c(neighbors, (i-2)*n + j); weights <- c(weights, 1) }
#'     if(i < n) { neighbors <- c(neighbors, i*n + j); weights <- c(weights, 1) }
#'     if(j > 1) { neighbors <- c(neighbors, (i-1)*n + j-1); weights <- c(weights, 1) }
#'     if(j < n) { neighbors <- c(neighbors, (i-1)*n + j+1); weights <- c(weights, 1) }
#'
#'     adj.list[[idx]] <- neighbors
#'     weight.list[[idx]] <- weights
#'   }
#' }
#'
#' # Create function with multiple extrema
#' y <- sin(seq(0, 2*pi, length.out = n*n)) + 0.1*rnorm(n*n)
#'
#' # Apply smoothing with detailed recording
#' result <- create.gflow.cx(
#'   adj.list, weight.list, y,
#'   hop.idx.thld = 2,
#'   smoother.type = 2,  # Harmonic Eigen
#'   detailed.recording = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Examine extrema
#' print(result$extrema_df)
#' }
#'
#' @seealso
#' \code{\link{create.hop.nbhd.extrema.df}} for extracting extrema information,
#'
#' @export
create.gflow.cx <- function(adj.list,
                            weight.list,
                            y,
                            hop.idx.thld = 5,
                            smoother.type = 0,
                            max.outer.iterations = 5,
                            max.inner.iterations = 100,
                            smoothing.tolerance = 1e-6,
                            sigma = 1.0,
                            process.in.order = TRUE,
                            verbose = TRUE,
                            detailed.recording = FALSE) {

    ## Input validation
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }

    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }

    if (length(y) != n.vertices) {
        stop("y must have the same length as adj.list")
    }

    ## Check that each adjacency/weight pair has matching lengths
    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("adj.list[[%d]] and weight.list[[%d]] must have the same length", i, i))
        }

        ## Check for valid indices
        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n.vertices)) {
                stop(sprintf("Invalid vertex indices in adj.list[[%d]]. Indices must be between 1 and %d",
                            i, n.vertices))
            }
        }
    }

    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    if (any(is.na(y))) {
        stop("y cannot contain NA values")
    }

    if (!is.numeric(hop.idx.thld) || length(hop.idx.thld) != 1 || hop.idx.thld < 0) {
        stop("hop.idx.thld must be a single non-negative numeric value")
    }

    if (!smoother.type %in% 0:4) {
        stop("smoother.type must be an integer between 0 and 4")
    }

    if (!is.numeric(max.outer.iterations) || length(max.outer.iterations) != 1 ||
        max.outer.iterations < 1) {
        stop("max.outer.iterations must be a positive integer")
    }

    if (!is.numeric(max.inner.iterations) || length(max.inner.iterations) != 1 ||
        max.inner.iterations < 1) {
        stop("max.inner.iterations must be a positive integer")
    }

    if (!is.numeric(smoothing.tolerance) || length(smoothing.tolerance) != 1 ||
        smoothing.tolerance <= 0) {
        stop("smoothing.tolerance must be a positive numeric value")
    }

    if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
        stop("sigma must be a positive numeric value")
    }

    if (!is.logical(process.in.order) || length(process.in.order) != 1) {
        stop("process.in.order must be a single logical value")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a single logical value")
    }

    if (!is.logical(detailed.recording) || length(detailed.recording) != 1) {
        stop("detailed.recording must be a single logical value")
    }

    ## Create a mapping between smoother.type integers and descriptive names
    smoother.names <- c(
        "Weighted Mean",
        "Harmonic Iterative",
        "Harmonic Eigen",
        "Hybrid Biharmonic-Harmonic",
        "Boundary Smoothed Harmonic"
    )

    ## Print information about the chosen smoother
    if (verbose) {
        message("Using smoother: ", smoother.names[smoother.type + 1])
        message(sprintf("Processing %d vertices with hop index threshold %.1f",
                       n.vertices, hop.idx.thld))
    }

    ## Convert to 0-based indices for C++
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) == 0) {
            integer(0)
        } else {
            as.integer(x - 1)
        }
    })

    result <- .Call(
        "S_create_gflow_cx",
        adj.list.0based,
        weight.list,
        as.numeric(y),
        as.integer(hop.idx.thld),
        as.integer(smoother.type),
        as.integer(max.outer.iterations),
        as.integer(max.inner.iterations),
        as.numeric(smoothing.tolerance),
        as.numeric(sigma),
        as.logical(process.in.order),
        as.logical(verbose),
        as.logical(detailed.recording)
    )

    ## Create extrema data frames
    result$extrema_df <- create.hop.nbhd.extrema.df(result)

    result$extrema_df2 <- create.hop.nbhd.extrema.df(
        result,
        include_spurious = TRUE,
        threshold = hop.idx.thld,
        sort_by_value = FALSE,
        vertex_column_name = "vertex"
    )

    ## Add attributes
    attr(result, "smoother.type") <- smoother.type
    attr(result, "smoother.name") <- smoother.names[smoother.type + 1]
    attr(result, "hop.idx.threshold") <- hop.idx.thld

    ## Add class for custom print/plot methods
    class(result) <- c("gflow_cx", class(result))

    return(result)
}

#' Print Method for gflow_cx Objects
#'
#' @description
#' Provides a concise summary of a graph flow complex object, including information
#' about the smoother used, extrema counts, and spurious extrema identification.
#'
#' @param x A \code{gflow_cx} object returned by \code{\link{create.gflow.cx}}
#' @param ... Additional arguments passed to internal print functions (currently unused)
#'
#' @return Invisibly returns the input object \code{x}
#'
#' @examples
#' \dontrun{
#' # Create simple example
#' adj.list <- list(c(2,3), c(1,3), c(1,2))
#' weight.list <- list(c(1,1), c(1,1), c(1,1))
#' y <- c(0.5, 1.0, 0.7)
#'
#' result <- create.gflow.cx(adj.list, weight.list, y, verbose = FALSE)
#' print(result)
#' }
#' @export
print.gflow_cx <- function(x, ...) {
    cat("Graph Flow Complex with Harmonic Extension\n")
    cat("------------------------------------------\n")
    cat("Smoother:", attr(x, "smoother.name"), "\n")
    cat("Hop index threshold:", attr(x, "hop.idx.threshold"), "\n")

    n_lmin <- length(x$lmin_hop_nbhds)
    n_lmax <- length(x$lmax_hop_nbhds)

    threshold <- attr(x, "hop.idx.threshold")
    if (is.null(threshold)) threshold <- 2

    n_spurious_lmin <- sum(vapply(x$lmin_hop_nbhds,
                                  function(nbhd) nbhd$hop_idx <= threshold,
                                  logical(1)))
    n_spurious_lmax <- sum(vapply(x$lmax_hop_nbhds,
                                  function(nbhd) nbhd$hop_idx <= threshold,
                                  logical(1)))

    cat(sprintf("Local minima: %d (%d spurious)\n", n_lmin, n_spurious_lmin))
    cat(sprintf("Local maxima: %d (%d spurious)\n", n_lmax, n_spurious_lmax))

    if ("smoothing_history" %in% names(x)) {
        cat("Detailed smoothing history:", length(x$smoothing_history), "steps recorded\n")
    }

    cat("Predictions:", length(x$harmonic_predictions), "values\n")

    invisible(x)
}

#' Create Data Frame of Extrema Information from Hop Neighborhoods
#'
#' @description
#' Extracts and formats extrema information from hop neighborhood results into a
#'     convenient data frame format. This function works with results from
#'     \code{\link{create.gflow.cx}}.
#'
#' @param result An object containing hop neighborhood information, typically
#'     from \code{create.gflow.cx()}. Must contain components named
#'     \code{lmin_hop_nbhds} and \code{lmax_hop_nbhds}.
#' @param include_spurious Logical; whether to include spurious extrema in the
#'     output. Default is TRUE. Spurious extrema are those with hop index less
#'     than or equal to the threshold.
#' @param threshold Numeric scalar specifying the threshold for identifying
#'     spurious extrema. Extrema with hop_idx <= threshold are marked as
#'     spurious. If NULL (default), attempts to extract the threshold from the
#'     result object or uses 2.
#' @param sort_by_value Logical; if TRUE, orders extrema by their function
#'     values (ascending for minima, descending for maxima). If FALSE (default),
#'     preserves the original order.
#' @param vertex_column_name Character string specifying the name for the vertex
#'     column in the output data frame. Default is "vertex".
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{vertex}{Integer vertex indices (1-based)}
#'     \item{hop_idx}{Numeric hop index values for each extremum}
#'     \item{is_max}{Integer indicator: 0 for minima, 1 for maxima}
#'     \item{label}{Character labels for extrema (e.g., "m1", "m2" for minima,
#'       "M1", "M2" for maxima)}
#'     \item{value}{Numeric function values at extrema (NA if not available)}
#'     \item{spurious}{Logical indicating whether the extremum is spurious}
#'   }
#'
#'   The vertex column may be renamed based on \code{vertex_column_name}.
#'
#' @details
#' Labels are assigned to extrema based on their function values:
#' \itemize{
#'   \item Minima are labeled "m1", "m2", etc., in ascending order of function value
#'   \item Maxima are labeled "M1", "M2", etc., in descending order of function value
#' }
#'
#' This labeling scheme ensures that "m1" is the global minimum among detected
#' minima, and "M1" is the global maximum among detected maxima.
#'
#' @examples
#' \dontrun{
#' # Create example graph and compute hop neighborhoods
#' adj.list <- list(c(2,3), c(1,3,4), c(1,2,4), c(2,3))
#' weight.list <- list(c(1,1), c(1,1,1), c(1,1,1), c(1,1))
#' y <- c(0, 1, 0.5, 2)  # Minima at vertices 1,3 and maxima at vertices 2,4
#'
#' # Using with create.gflow.cx
#' result <- create.gflow.cx(adj.list, weight.list, y, verbose = FALSE)
#' extrema_df <- create.hop.nbhd.extrema.df(result)
#' print(extrema_df)
#'
#' # Exclude spurious extrema
#' significant_extrema <- create.hop.nbhd.extrema.df(
#'   result,
#'   include_spurious = FALSE,
#'   threshold = 2
#' )
#' print(significant_extrema)
#' }
#' @export
create.hop.nbhd.extrema.df <- function(result,
                                       include_spurious = TRUE,
                                       threshold = NULL,
                                       sort_by_value = FALSE,
                                       vertex_column_name = "vertex") {

    ## Input validation
    if (!is.list(result)) {
        stop("result must be a list")
    }

    if (!all(c("lmin_hop_nbhds", "lmax_hop_nbhds") %in% names(result))) {
        stop("Input must have 'lmin_hop_nbhds' and 'lmax_hop_nbhds' components")
    }

    if (!is.logical(include_spurious) || length(include_spurious) != 1) {
        stop("include_spurious must be a single logical value")
    }

    if (!is.logical(sort_by_value) || length(sort_by_value) != 1) {
        stop("sort_by_value must be a single logical value")
    }

    if (!is.character(vertex_column_name) || length(vertex_column_name) != 1) {
        stop("vertex_column_name must be a single character string")
    }

    ## Determine if we're dealing with gflow_cx object
    is_gflow <- inherits(result, "gflow_cx")

    ## Determine threshold for spurious extrema
    if (is.null(threshold)) {
        if (is_gflow && !is.null(attr(result, "hop.idx.threshold"))) {
            threshold <- attr(result, "hop.idx.threshold")
        } else {
            threshold <- 2  # Default if not found
        }
    }

    if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
        stop("threshold must be a single non-negative numeric value")
    }

    ## Extract data from minima
    lmin_data <- lapply(result$lmin_hop_nbhds, function(x) {
        list(
            vertex = x$vertex,
            hop_idx = x$hop_idx,
            value = if ("value" %in% names(x)) x$value else NA_real_
        )
    })

    ## Extract data from maxima
    lmax_data <- lapply(result$lmax_hop_nbhds, function(x) {
        list(
            vertex = x$vertex,
            hop_idx = x$hop_idx,
            value = if ("value" %in% names(x)) x$value else NA_real_
        )
    })

    n_min <- length(lmin_data)
    n_max <- length(lmax_data)

    ## Handle empty case
    if (n_min + n_max == 0) {
        df <- data.frame(
            vertex = integer(),
            hop_idx = numeric(),
            is_max = integer(),
            label = character(),
            value = numeric(),
            spurious = logical(),
            stringsAsFactors = FALSE
        )
        names(df)[1] <- vertex_column_name
        return(df)
    }

    ## Extract vectors
    if (n_min > 0) {
        min_vertices <- vapply(lmin_data, `[[`, integer(1), "vertex")
        min_hop_indices <- vapply(lmin_data, `[[`, numeric(1), "hop_idx")
        min_values <- vapply(lmin_data, `[[`, numeric(1), "value")
    } else {
        min_vertices <- integer()
        min_hop_indices <- numeric()
        min_values <- numeric()
    }

    if (n_max > 0) {
        max_vertices <- vapply(lmax_data, `[[`, integer(1), "vertex")
        max_hop_indices <- vapply(lmax_data, `[[`, numeric(1), "hop_idx")
        max_values <- vapply(lmax_data, `[[`, numeric(1), "value")
    } else {
        max_vertices <- integer()
        max_hop_indices <- numeric()
        max_values <- numeric()
    }

    ## Combine data
    vertices <- c(min_vertices, max_vertices)
    hop_indices <- c(min_hop_indices, max_hop_indices)
    is_max <- c(rep(0L, n_min), rep(1L, n_max))
    values <- c(min_values, max_values)

    ## Determine spurious extrema
    spurious <- hop_indices <= threshold

    ## Filter spurious extrema if requested
    if (!include_spurious) {
        keep <- !spurious
        vertices <- vertices[keep]
        hop_indices <- hop_indices[keep]
        is_max <- is_max[keep]
        values <- values[keep]
        spurious <- spurious[keep]

        # Update counts
        n_min <- sum(is_max == 0)
        n_max <- sum(is_max == 1)
    }

    ## Create labels
    labels <- character(length(vertices))

    if (n_min > 0) {
        min_idx <- which(is_max == 0)
        if (all(is.na(values[min_idx]))) {
            # No values available, use original order
            labels[min_idx] <- paste0("m", seq_len(n_min))
        } else {
            # Sort by value (ascending for minima)
            ord <- order(values[min_idx], na.last = TRUE)
            labels[min_idx[ord]] <- paste0("m", seq_len(n_min))
        }
    }

    if (n_max > 0) {
        max_idx <- which(is_max == 1)
        if (all(is.na(values[max_idx]))) {
            # No values available, use original order
            labels[max_idx] <- paste0("M", seq_len(n_max))
        } else {
            # Sort by value (descending for maxima)
            ord <- order(values[max_idx], decreasing = TRUE, na.last = TRUE)
            labels[max_idx[ord]] <- paste0("M", seq_len(n_max))
        }
    }

    ## Create data frame
    df <- data.frame(
        vertex = vertices,
        hop_idx = hop_indices,
        is_max = is_max,
        label = labels,
        value = values,
        spurious = spurious,
        stringsAsFactors = FALSE
    )

    ## Rename vertex column if requested
    if (vertex_column_name != "vertex") {
        names(df)[1] <- vertex_column_name
    }

    ## Sort if requested
    if (sort_by_value && !all(is.na(values))) {
        # Sort minima ascending, maxima descending
        min_order <- if (n_min > 0) order(df$value[df$is_max == 0], na.last = TRUE) else integer()
        max_order <- if (n_max > 0) order(df$value[df$is_max == 1], decreasing = TRUE, na.last = TRUE) else integer()

        min_rows <- which(df$is_max == 0)[min_order]
        max_rows <- which(df$is_max == 1)[max_order]

        df <- df[c(min_rows, max_rows), ]
        rownames(df) <- NULL
    }

    return(df)
}

#' Visualize Smoothing Steps from Gradient Flow Complex
#'
#' Creates an animated 3D visualization showing how spurious extrema are smoothed
#' during the gradient flow complex construction process.
#'
#' @param gflow_result Result from \code{create.gflow.cx(..., detailed.recording = TRUE)}.
#'        Must contain \code{$smoothing_history} (list of steps).
#' @param plot_res A list containing at least \code{$layout}: an \eqn{n\times 2} numeric matrix
#'        of 2D coordinates for vertices. If \code{$graph} (igraph) or adjacency data
#'        are present, edges will also be rendered when available.
#' @param step_by_step Logical; if TRUE, waits for user input between steps (interactive only).
#'        In non-interactive/headless sessions this is ignored and snapshot mode is used.
#' @param animation_delay Numeric; delay in seconds between frames (used when not waiting).
#' @param out.dir Optional output directory for snapshot PNGs (default: \code{tempdir()}).
#' @param filename.prefix Character; prefix for snapshot filenames (default: \code{"smoothing_step"}).
#'
#' @return NULL invisibly.
#'
#' @details
#' The function opens a single rgl window and reuses it for all frames.
#' In headless/CI environments, an off-screen rgl device is used automatically and
#' images are written to \code{out.dir}.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   # result <- create.gflow.cx(..., detailed.recording = TRUE)
#'   # plot_res <- list(layout = <n x 2 matrix>, graph = <igraph optional>)
#'   # visualize.smoothing.steps(result, plot_res)
#' }
#' }
#' @export
visualize.smoothing.steps <- function(gflow_result,
                                      plot_res,
                                      step_by_step = TRUE,
                                      animation_delay = 0.5,
                                      out.dir = NULL,
                                      filename.prefix = "smoothing_step") {
  ## ---- Preconditions --------------------------------------------------------
  if (!"smoothing_history" %in% names(gflow_result)) {
    stop("No detailed smoothing history available. ",
         "Run create.gflow.cx(..., detailed.recording = TRUE).", call. = FALSE)
  }
  steps <- gflow_result$smoothing_history
  n_steps <- length(steps)

  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
  }

  if (!is.list(plot_res) || is.null(plot_res$layout)) {
    stop("'plot_res' must be a list with component 'layout' (n x 2 numeric matrix).", call. = FALSE)
  }
  layout_2d <- as.matrix(plot_res$layout)
  if (!is.numeric(layout_2d) || ncol(layout_2d) != 2L) {
    stop("'plot_res$layout' must be numeric with exactly 2 columns.", call. = FALSE)
  }

  if (!is.logical(step_by_step) || length(step_by_step) != 1L) {
    stop("'step_by_step' must be a single logical.", call. = FALSE)
  }
  if (!is.numeric(animation_delay) || length(animation_delay) != 1L || animation_delay < 0) {
    stop("'animation_delay' must be a non-negative numeric.", call. = FALSE)
  }

  ## ---- Edge extraction (optional) ------------------------------------------
  # Try to get edges for rendering (best-effort; OK if absent)
  get_edges <- function() {
    if (!is.null(plot_res$graph)) {
      e <- try(igraph::as_edgelist(plot_res$graph, names = FALSE), silent = TRUE)
      if (!inherits(e, "try-error") && is.matrix(e) && ncol(e) == 2L) return(e)
    }
    if (!is.null(gflow_result$adj.list)) {
      el <- lapply(seq_along(gflow_result$adj.list), function(i) {
        v <- gflow_result$adj.list[[i]]
        if (!length(v)) return(NULL)
        cbind(i, as.integer(v))
      })
      e <- do.call(rbind, el)
      if (!is.null(e)) {
        e <- e[e[, 1] < e[, 2], , drop = FALSE] # undirected unique
        return(e)
      }
    }
    NULL
  }
  edges <- get_edges()

  ## ---- Device mode: real vs off-screen -------------------------------------
  # Use a real window only if interactive & display available; otherwise off-screen.
  use_null <- (!interactive()) ||
              identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
              (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
  old_opt <- options(rgl.useNULL = use_null)
  on.exit(options(old_opt), add = TRUE)

  # In non-interactive/headless, don't block: force snapshot mode
  if (isTRUE(step_by_step) && (use_null || !interactive())) {
    message("Non-interactive/headless context detected; switching to snapshot mode.")
    step_by_step <- FALSE
  }

  # Snapshot directory (dont pollute working dir by default)
  if (is.null(out.dir)) out.dir <- getOption("gflow.snapshot.dir", default = tempdir())
  if (!dir.exists(out.dir)) {
    ok <- dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop("Unable to create output directory: ", out.dir, call. = FALSE)
  }

  ## ---- Open a single rgl device --------------------------------------------
  rgl::open3d()
  on.exit(try(rgl::rgl.close(), silent = TRUE), add = TRUE)

  ## ---- 3D drawing helper ----------------------------------------------------
  draw_frame <- function(values, title_text, highlight_vertex = NULL) {
    # validate values length
    if (length(values) != nrow(layout_2d)) {
      stop("Length of 'values' does not match number of vertices in layout.", call. = FALSE)
    }

    layout_3d <- cbind(layout_2d, as.numeric(values))

    rgl::clear3d()

    # base plane points
    rgl::points3d(cbind(layout_2d, 0), size = 4, col = "gray70", alpha = 0.35)

    # base plane edges (optional)
    if (!is.null(edges) && nrow(edges) > 0) {
      base_seg <- matrix(NA_real_, nrow = nrow(edges) * 2L, ncol = 3L)
      for (i in seq_len(nrow(edges))) {
        v1 <- edges[i, 1]; v2 <- edges[i, 2]
        base_seg[(i - 1L) * 2L + 1L, ] <- c(layout_2d[v1, ], 0)
        base_seg[(i - 1L) * 2L + 2L, ] <- c(layout_2d[v2, ], 0)
      }
      rgl::segments3d(base_seg, col = "gray60", alpha = 0.35, lwd = 1)
    }

    # 3D vertices (heights)
    # color by rank for a quick gradient
    pal <- grDevices::hcl.colors(nrow(layout_3d), palette = "Spectral")
    cols <- pal[rank(values, ties.method = "average")]
    rgl::spheres3d(layout_3d, radius = 0.02, col = cols, alpha = 0.9)

    # 3D edges (optional)
    if (!is.null(edges) && nrow(edges) > 0) {
      for (i in seq_len(nrow(edges))) {
        v1 <- edges[i, 1]; v2 <- edges[i, 2]
        rgl::lines3d(rbind(layout_3d[v1, ], layout_3d[v2, ]), col = "gray50", lwd = 1)
      }
    }

    # Highlight current vertex (optional)
    if (!is.null(highlight_vertex) && is.finite(highlight_vertex)) {
      hv <- as.integer(highlight_vertex)
      if (hv >= 1L && hv <= nrow(layout_3d)) {
        rgl::spheres3d(layout_3d[hv, , drop = FALSE], radius = 0.05,
                       col = "yellow", alpha = 0.8)
      }
    }

    # axes + title
    rgl::axes3d()
    rgl::title3d(xlab = "X", ylab = "Y", zlab = "Value", main = title_text)
  }

  ## ---- Run animation --------------------------------------------------------
  if (n_steps == 0) {
    message("No smoothing steps recorded (no spurious extrema found).")
    return(invisible(NULL))
  }

  message(sprintf("Found %d smoothing step%s.", n_steps, if (n_steps == 1) "" else "s"))

  # Initial frame
  initial_values <- steps[[1]]$before
  draw_frame(initial_values, "Initial State")

  if (isTRUE(step_by_step)) {
    if (interactive()) {
      invisible(readline("Press [Enter] to start smoothing process..."))
      for (i in seq_len(n_steps)) {
        st <- steps[[i]]
        title_before <- sprintf("Step %d/%d: Before smoothing %s at vertex %d (hop_idx=%.2f)",
                                i, n_steps, ifelse(isTRUE(st$is_minimum), "minimum", "maximum"),
                                st$vertex, as.numeric(st$hop_idx))
        draw_frame(st$before, title_before, st$vertex)
        if (animation_delay > 0) Sys.sleep(animation_delay)

        delta <- abs(st$after[st$vertex] - st$before[st$vertex])
        title_after <- sprintf("Step %d/%d: After smoothing (\u0394=%.4f)", i, n_steps, delta)
        draw_frame(st$after, title_after, st$vertex)

        if (i < n_steps) invisible(readline("Press [Enter] to continue..."))
      }
    }
  } else {
    # Snapshot mode
    # helper for filenames
    fpath <- function(tag) file.path(out.dir, sprintf("%s_%s.png", filename.prefix, tag))

    message("Saving snapshots to: ", normalizePath(out.dir, winslash = "/"))
    # initial
    draw_frame(initial_values, "Initial State")
    rgl::snapshot3d(filename = fpath("000_initial"))

    for (i in seq_len(n_steps)) {
      st <- steps[[i]]

      title_before <- sprintf("Step %d/%d: Before smoothing %s at vertex %d",
                              i, n_steps, ifelse(isTRUE(st$is_minimum), "minimum", "maximum"), st$vertex)
      draw_frame(st$before, title_before, st$vertex)
      rgl::snapshot3d(filename = fpath(sprintf("%03d_before", i)))
      if (animation_delay > 0) Sys.sleep(animation_delay)

      title_after <- sprintf("Step %d/%d: After smoothing", i, n_steps)
      draw_frame(st$after, title_after, st$vertex)
      rgl::snapshot3d(filename = fpath(sprintf("%03d_after", i)))
      if (animation_delay > 0) Sys.sleep(animation_delay)
    }

    message(sprintf("Saved %d snapshot%s.",
                    2 * n_steps + 1L, if (2 * n_steps + 1L == 1L) "" else "s"))
  }

  # Final frame
  final_values <- steps[[n_steps]]$after
  draw_frame(final_values, "Final State (All Spurious Extrema Removed)")
  if (isTRUE(step_by_step) && interactive()) {
    invisible(readline("Press [Enter] to close..."))
  }

  invisible(NULL)
}
