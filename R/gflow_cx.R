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
        S_create_gflow_cx,
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
#' Extracts and formats local extrema information from hop-neighborhood results
#' into a data frame.
#'
#' @param result A list-like object containing hop neighborhood information.
#'   Must contain components \code{lmin_hop_nbhds} and \code{lmax_hop_nbhds}.
#'   If \code{inherits(result, "gflow_cx")}, the spurious threshold may be read
#'   from an attribute (see \code{threshold}).
#' @param include_spurious Logical scalar. If FALSE, remove extrema with
#'   \code{hop_idx <= threshold}.
#' @param threshold Numeric scalar or NULL. If NULL and \code{result} is a
#'   \code{gflow_cx}, attempts to read \code{attr(result, "hop.idx.threshold")}
#'   and then \code{attr(result, "hop_idx_threshold")} for backward compatibility.
#'   If still NULL, defaults to 2.
#' @param sort_by_value Logical scalar. If TRUE, labels (and optionally rows)
#'   are ordered by function value: minima ascending, maxima descending.
#' @param vertex_column_name Character scalar. Name of the vertex column.
#' @param label_style Character scalar. Either \code{"mM"} (default; labels m1/M1)
#'   or \code{"minmax"} (labels min1/max1).
#' @param reorder Logical scalar. If TRUE, reorder rows so minima (ascending value)
#'   come first, followed by maxima (descending value). Defaults to \code{sort_by_value}.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item vertex (or renamed): integer vertex index
#'     \item hop_idx: numeric hop index
#'     \item is_max: integer 0 (min) / 1 (max)
#'     \item label: character extremum label
#'     \item value: numeric extremum value (NA if unavailable)
#'     \item spurious: logical flag (\code{hop_idx <= threshold})
#'   }
#'
#' @export
create.hop.nbhd.extrema.df <- function(result,
                                       include_spurious = TRUE,
                                       threshold = NULL,
                                       sort_by_value = FALSE,
                                       vertex_column_name = "vertex",
                                       label_style = c("mM", "minmax"),
                                       reorder = sort_by_value) {

    ## ------------------------------------------------------------------------
    ## Input validation
    ## ------------------------------------------------------------------------

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
    if (!is.logical(reorder) || length(reorder) != 1) {
        stop("reorder must be a single logical value")
    }
    if (!is.character(vertex_column_name) || length(vertex_column_name) != 1) {
        stop("vertex_column_name must be a single character string")
    }

    label_style <- match.arg(label_style)

    ## ------------------------------------------------------------------------
    ## Threshold resolution (support both attribute names)
    ## ------------------------------------------------------------------------

    is.gflow <- inherits(result, "gflow_cx")

    if (is.null(threshold)) {
        if (is.gflow) {
            threshold <- attr(result, "hop.idx.threshold")
            if (is.null(threshold)) {
                threshold <- attr(result, "hop_idx_threshold")
            }
        }
        if (is.null(threshold)) {
            threshold <- 2
        }
    }

    if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold < 0) {
        stop("threshold must be a single non-negative numeric value")
    }

    ## ------------------------------------------------------------------------
    ## Extract extrema vectors
    ## ------------------------------------------------------------------------

    extract.extrema <- function(x.list) {
        n <- length(x.list)
        if (n == 0) {
            return(list(
                vertex = integer(0),
                hop_idx = numeric(0),
                value = numeric(0)
            ))
        }

        vertex <- vapply(x.list, function(x) as.integer(x$vertex), integer(1))
        hop.idx <- vapply(x.list, function(x) as.double(x$hop_idx), numeric(1))
        value <- vapply(
            x.list,
            function(x) {
                if ("value" %in% names(x)) as.double(x$value) else NA_real_
            },
            numeric(1)
        )

        list(vertex = vertex, hop_idx = hop.idx, value = value)
    }

    min.data <- extract.extrema(result$lmin_hop_nbhds)
    max.data <- extract.extrema(result$lmax_hop_nbhds)

    n.min <- length(min.data$vertex)
    n.max <- length(max.data$vertex)

    ## Empty case
    if (n.min + n.max == 0) {
        df <- data.frame(
            vertex = integer(0),
            hop_idx = numeric(0),
            is_max = integer(0),
            label = character(0),
            value = numeric(0),
            spurious = logical(0),
            stringsAsFactors = FALSE
        )
        names(df)[1] <- vertex_column_name
        return(df)
    }

    vertex <- c(min.data$vertex, max.data$vertex)
    hop.idx <- c(min.data$hop_idx, max.data$hop_idx)
    is.max <- c(rep(0L, n.min), rep(1L, n.max))
    value <- c(min.data$value, max.data$value)

    spurious <- hop.idx <= threshold

    if (!include_spurious) {
        keep <- !spurious
        vertex <- vertex[keep]
        hop.idx <- hop.idx[keep]
        is.max <- is.max[keep]
        value <- value[keep]
        spurious <- spurious[keep]
        n.min <- sum(is.max == 0L)
        n.max <- sum(is.max == 1L)
    }

    ## ------------------------------------------------------------------------
    ## Labeling
    ## ------------------------------------------------------------------------

    min.prefix <- if (label_style == "mM") "m" else "min"
    max.prefix <- if (label_style == "mM") "M" else "max"

    labels <- character(length(vertex))

    if (n.min > 0) {
        idx <- which(is.max == 0L)
        if (sort_by_value && !all(is.na(value[idx]))) {
            ord <- order(value[idx], na.last = TRUE)
            labels[idx[ord]] <- paste0(min.prefix, seq_len(n.min))
        } else {
            labels[idx] <- paste0(min.prefix, seq_len(n.min))
        }
    }

    if (n.max > 0) {
        idx <- which(is.max == 1L)
        if (sort_by_value && !all(is.na(value[idx]))) {
            ord <- order(value[idx], decreasing = TRUE, na.last = TRUE)
            labels[idx[ord]] <- paste0(max.prefix, seq_len(n.max))
        } else {
            labels[idx] <- paste0(max.prefix, seq_len(n.max))
        }
    }

    df <- data.frame(
        vertex = vertex,
        hop_idx = hop.idx,
        is_max = is.max,
        label = labels,
        value = value,
        spurious = spurious,
        stringsAsFactors = FALSE
    )

    ## Rename vertex column
    if (vertex_column_name != "vertex") {
        names(df)[1] <- vertex_column_name
    }

    ## Optional reordering
    if (reorder && sort_by_value && !all(is.na(df$value))) {
        min.rows <- which(df$is_max == 0L)
        max.rows <- which(df$is_max == 1L)

        min.ord <- if (length(min.rows) > 0) order(df$value[min.rows], na.last = TRUE) else integer(0)
        max.ord <- if (length(max.rows) > 0) order(df$value[max.rows], decreasing = TRUE, na.last = TRUE) else integer(0)

        df <- df[c(min.rows[min.ord], max.rows[max.ord]), , drop = FALSE]
        rownames(df) <- NULL
    }

    ## Ensure vertex column is first (already is), but keep explicit ordering stable
    col.order <- c(vertex_column_name, "hop_idx", "is_max", "label", "value", "spurious")
    df <- df[, col.order, drop = FALSE]

    df
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
    ## Try to get edges for rendering (best-effort; OK if absent)
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
    ## Use a real window only if interactive & display available; otherwise off-screen.
    use_null <- (!interactive()) ||
        identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
        (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
    old_opt <- options(rgl.useNULL = use_null)
    on.exit(options(old_opt), add = TRUE)

    ## In non-interactive/headless, don't block: force snapshot mode
    if (isTRUE(step_by_step) && (use_null || !interactive())) {
        message("Non-interactive/headless context detected; switching to snapshot mode.")
        step_by_step <- FALSE
    }

    ## Snapshot directory (dont pollute working dir by default)
    if (is.null(out.dir)) out.dir <- getOption("gflow.snapshot.dir", default = tempdir())
    if (!dir.exists(out.dir)) {
        ok <- dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok) stop("Unable to create output directory: ", out.dir, call. = FALSE)
    }

    ## ---- Open a single rgl device --------------------------------------------
    rgl::open3d()
    if (use_null) {
        ## Only close the device automatically if using null device
        on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
    }
    rgl::clear3d()

    ## ---- 3D drawing helper ----------------------------------------------------
    draw_frame <- function(values, title_text, highlight_vertex = NULL) {
        ## validate values length
        if (length(values) != nrow(layout_2d)) {
            stop("Length of 'values' does not match number of vertices in layout.", call. = FALSE)
        }

        layout_3d <- cbind(layout_2d, as.numeric(values))

        rgl::clear3d()

        ## base plane points
        rgl::points3d(cbind(layout_2d, 0), size = 4, col = "gray70", alpha = 0.35)

        ## base plane edges (optional)
        if (!is.null(edges) && nrow(edges) > 0) {
            base_seg <- matrix(NA_real_, nrow = nrow(edges) * 2L, ncol = 3L)
            for (i in seq_len(nrow(edges))) {
                v1 <- edges[i, 1]; v2 <- edges[i, 2]
                base_seg[(i - 1L) * 2L + 1L, ] <- c(layout_2d[v1, ], 0)
                base_seg[(i - 1L) * 2L + 2L, ] <- c(layout_2d[v2, ], 0)
            }
            rgl::segments3d(base_seg, col = "gray60", alpha = 0.35, lwd = 1)
        }

        ## 3D vertices (heights)
        ## color by rank for a quick gradient
        pal <- grDevices::hcl.colors(nrow(layout_3d), palette = "Spectral")
        cols <- pal[rank(values, ties.method = "average")]
        rgl::spheres3d(layout_3d, radius = 0.02, col = cols, alpha = 0.9)

        ## 3D edges (optional)
        if (!is.null(edges) && nrow(edges) > 0) {
            for (i in seq_len(nrow(edges))) {
                v1 <- edges[i, 1]; v2 <- edges[i, 2]
                rgl::lines3d(rbind(layout_3d[v1, ], layout_3d[v2, ]), col = "gray50", lwd = 1)
            }
        }

        ## Highlight current vertex (optional)
        if (!is.null(highlight_vertex) && is.finite(highlight_vertex)) {
            hv <- as.integer(highlight_vertex)
            if (hv >= 1L && hv <= nrow(layout_3d)) {
                rgl::spheres3d(layout_3d[hv, , drop = FALSE], radius = 0.05,
                               col = "yellow", alpha = 0.8)
            }
        }

        ## axes + title
        rgl::axes3d()
        rgl::title3d(xlab = "X", ylab = "Y", zlab = "Value", main = title_text)
    }

    ## ---- Run animation --------------------------------------------------------
    if (n_steps == 0) {
        message("No smoothing steps recorded (no spurious extrema found).")
        return(invisible(NULL))
    }

    message(sprintf("Found %d smoothing step%s.", n_steps, if (n_steps == 1) "" else "s"))

    ## Initial frame
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
        ## Snapshot mode
        ## helper for filenames
        fpath <- function(tag) file.path(out.dir, sprintf("%s_%s.png", filename.prefix, tag))

        message("Saving snapshots to: ", normalizePath(out.dir, winslash = "/"))
        ## initial
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

    ## Final frame
    final_values <- steps[[n_steps]]$after
    draw_frame(final_values, "Final State (All Spurious Extrema Removed)")
    if (isTRUE(step_by_step) && interactive()) {
        invisible(readline("Press [Enter] to close..."))
    }

    invisible(NULL)
}

#' Classify Spurious Extrema Using Multiple Criteria
#'
#' @description
#' Determines which extrema should be considered spurious using various criteria.
#' Provides more sophisticated classification than a simple hop_idx threshold.
#'
#' @param hop_indices Numeric vector of hop indices for all extrema
#' @param values Optional numeric vector of function values at extrema
#' @param method Character string specifying classification method:
#'   - "threshold": hop_idx <= threshold (default, backward compatible)
#'   - "quantile": hop_idx <= quantile(hop_idx, probs = quantile_level)
#'   - "outlier": statistical outlier detection based on hop_idx distribution
#'   - "combined": uses both threshold and relative criteria
#' @param threshold Numeric threshold for "threshold" and "combined" methods (default 2)
#' @param quantile_level Quantile level for "quantile" method (default 0.25)
#' @param verbose Logical; print diagnostic information (default FALSE)
#'
#' @return Logical vector indicating which extrema are spurious
#'
#' @examples
#' \dontrun{
#' # Simple threshold (backward compatible)
#' spurious <- classify.spurious.extrema(hop_indices, method = "threshold", threshold = 2)
#'
#' # Quantile-based (adaptive to data)
#' spurious <- classify.spurious.extrema(hop_indices, method = "quantile", quantile_level = 0.25)
#'
#' # Combined approach
#' spurious <- classify.spurious.extrema(hop_indices, values, method = "combined")
#' }
#'
#' @export
classify.spurious.extrema <- function(hop_indices,
                                     values = NULL,
                                     method = c("threshold", "quantile", "outlier", "combined"),
                                     threshold = 2,
                                     quantile_level = 0.25,
                                     verbose = FALSE) {

  method <- match.arg(method)
  n <- length(hop_indices)

  # Handle edge cases
  if (n == 0) return(logical(0))
  if (n == 1) return(FALSE)  # Single extremum is not spurious

  # Remove infinite values (global extrema) for analysis
  finite_indices <- hop_indices[is.finite(hop_indices)]

  if (length(finite_indices) == 0) {
    # All extrema are global - none spurious
    return(rep(FALSE, n))
  }

  if (verbose) {
    cat("Spurious extrema classification:\n")
    cat("  Method:", method, "\n")
    cat("  Total extrema:", n, "\n")
    cat("  Global extrema (hop_idx=Inf):", sum(is.infinite(hop_indices)), "\n")
    cat("  Finite hop_idx range: [", min(finite_indices), ",", max(finite_indices), "]\n")
    cat("  Median hop_idx:", median(finite_indices), "\n")
  }

  spurious <- switch(method,

    threshold = {
      # Simple threshold approach (original)
      result <- hop_indices <= threshold
      if (verbose) {
        cat("  Threshold:", threshold, "\n")
        cat("  Spurious count:", sum(result), "\n")
      }
      result
    },

    quantile = {
      # Quantile-based approach (adaptive to data distribution)
      cutoff <- quantile(finite_indices, probs = quantile_level)
      result <- is.finite(hop_indices) & hop_indices <= cutoff
      if (verbose) {
        cat("  Quantile level:", quantile_level, "\n")
        cat("  Cutoff value:", cutoff, "\n")
        cat("  Spurious count:", sum(result), "\n")
      }
      result
    },

    outlier = {
      # Statistical outlier detection using IQR method
      # Extrema with unusually low hop_idx are spurious
      q1 <- quantile(finite_indices, 0.25)
      q3 <- quantile(finite_indices, 0.75)
      iqr <- q3 - q1
      lower_bound <- q1 - 1.5 * iqr

      # In this context, we're looking for LOW outliers
      result <- is.finite(hop_indices) & hop_indices < lower_bound

      if (verbose) {
        cat("  Q1:", q1, ", Q3:", q3, ", IQR:", iqr, "\n")
        cat("  Lower bound:", lower_bound, "\n")
        cat("  Spurious count:", sum(result), "\n")
      }
      result
    },

    combined = {
      # Combined approach: must satisfy BOTH criteria
      # 1. hop_idx <= threshold (absolute criterion)
      # 2. hop_idx in lower quantile (relative criterion)

      cutoff <- quantile(finite_indices, probs = quantile_level)
      abs_criterion <- hop_indices <= threshold
      rel_criterion <- is.finite(hop_indices) & hop_indices <= cutoff

      result <- abs_criterion & rel_criterion

      if (verbose) {
        cat("  Absolute threshold:", threshold, "\n")
        cat("  Relative cutoff (", quantile_level, " quantile):", cutoff, "\n")
        cat("  Meeting absolute criterion:", sum(abs_criterion), "\n")
        cat("  Meeting relative criterion:", sum(rel_criterion), "\n")
        cat("  Meeting both (spurious):", sum(result), "\n")
      }
      result
    }
  )

  # Never mark global extrema as spurious
  spurious[is.infinite(hop_indices)] <- FALSE

  return(spurious)
}


#' Enhanced Version of create.hop.nbhd.extrema.df with Improved Spurious Classification
#'
#' @description
#' Creates a data frame containing information about local extrema with multiple
#' options for classifying spurious extrema.
#'
#' @param result An object from compute.extrema.hop.nbhds() or create.gflow.cx()
#' @param include_spurious Logical; whether to include spurious extrema (default TRUE)
#' @param spurious_method Method for determining spurious extrema (see classify.spurious.extrema)
#' @param threshold Threshold for hop_idx (default 2, or from result object)
#' @param quantile_level Quantile level for quantile-based spurious detection (default 0.25)
#' @param sort_by_value Logical; if TRUE, orders extrema by function value (default FALSE)
#' @param vertex_column_name Character; name of the vertex column (default "vertex")
#' @param verbose Logical; print diagnostic information (default FALSE)
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item vertex: Vertex index
#'     \item hop_idx: Hop index of the extremum
#'     \item is_max: 0 for minima, 1 for maxima
#'     \item label: Label for the extremum
#'     \item spurious: Logical; whether the extremum is considered spurious
#'     \item value: Function value at the extremum (if available)
#'   }
#'
#' @examples
#' \dontrun{
#' # Original behavior (backward compatible)
#' extrema_df <- create.hop.nbhd.extrema.df2(result)
#'
#' # Adaptive quantile-based spurious detection
#' extrema_df <- create.hop.nbhd.extrema.df2(result,
#'                                            spurious_method = "quantile",
#'                                            quantile_level = 0.33)
#'
#' # Combined approach with diagnostics
#' extrema_df <- create.hop.nbhd.extrema.df2(result,
#'                                            spurious_method = "combined",
#'                                            verbose = TRUE)
#' }
#'
#' @export
create.hop.nbhd.extrema.df2 <- function(result,
                                       include_spurious = TRUE,
                                       spurious_method = "threshold",
                                       threshold = NULL,
                                       quantile_level = 0.25,
                                       sort_by_value = FALSE,
                                       vertex_column_name = "vertex",
                                       verbose = FALSE) {

  # Check if the input has required components
  if (!all(c("lmin_hop_nbhds", "lmax_hop_nbhds") %in% names(result))) {
    stop("Input must have 'lmin_hop_nbhds' and 'lmax_hop_nbhds' components")
  }

  # Determine if we're dealing with gflow_cx object
  is_gflow <- inherits(result, "gflow_cx")

  # Determine threshold
  if (is.null(threshold)) {
    if (is_gflow) {
      threshold <- attr(result, "hop_idx_threshold")
    }
    if (is.null(threshold)) {
      threshold <- 2
    }
  }

  # Extract data from minima
  lmin_data <- lapply(result$lmin_hop_nbhds, function(x) {
    list(
      vertex = x$vertex,
      hop_idx = x$hop_idx,
      value = if ("value" %in% names(x)) x$value else NA
    )
  })

  # Extract data from maxima
  lmax_data <- lapply(result$lmax_hop_nbhds, function(x) {
    list(
      vertex = x$vertex,
      hop_idx = x$hop_idx,
      value = if ("value" %in% names(x)) x$value else NA
    )
  })

  n_min <- length(lmin_data)
  n_max <- length(lmax_data)

  # Handle empty case
  if (n_min + n_max == 0) {
    return(data.frame(
      vertex = integer(),
      hop_idx = numeric(),
      is_max = integer(),
      label = character(),
      spurious = logical(),
      value = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  # Extract data
  vertices <- c(sapply(lmin_data, `[[`, "vertex"), sapply(lmax_data, `[[`, "vertex"))
  hop_indices <- c(sapply(lmin_data, `[[`, "hop_idx"), sapply(lmax_data, `[[`, "hop_idx"))
  is_max <- c(rep(0, n_min), rep(1, n_max))
  values <- c(sapply(lmin_data, `[[`, "value"), sapply(lmax_data, `[[`, "value"))
  has_values <- !all(is.na(values))

  # Classify spurious extrema using enhanced method
  spurious <- classify.spurious.extrema(
    hop_indices = hop_indices,
    values = if (has_values) values else NULL,
    method = spurious_method,
    threshold = threshold,
    quantile_level = quantile_level,
    verbose = verbose
  )

  # Filter spurious extrema if requested
  if (!include_spurious) {
    keep <- !spurious
    vertices <- vertices[keep]
    hop_indices <- hop_indices[keep]
    is_max <- is_max[keep]
    values <- values[keep]
    spurious <- spurious[keep]

    # Adjust counts after filtering
    n_min <- sum(is_max == 0)
    n_max <- sum(is_max == 1)
  }

  # Create labels
  if (sort_by_value && has_values) {
    labels <- character(length(vertices))

    # Label maxima in descending order
    if (n_max > 0) {
      max_indices <- which(is_max == 1)
      max_values <- values[max_indices]
      ord <- order(max_values, decreasing = TRUE)
      labels[max_indices[ord]] <- paste0("M", seq_len(n_max))
    }

    # Label minima in ascending order
    if (n_min > 0) {
      min_indices <- which(is_max == 0)
      min_values <- values[min_indices]
      ord <- order(min_values, decreasing = FALSE)
      labels[min_indices[ord]] <- paste0("m", seq_len(n_min))
    }
  } else {
    # Sequential labeling
    labels <- c(
      paste0("m", seq_len(n_min)),
      paste0("M", seq_len(n_max))
    )
  }

  # Create data frame
  df <- data.frame(
    hop_idx = hop_indices,
    is_max = is_max,
    label = labels,
    spurious = spurious,
    stringsAsFactors = FALSE
  )

  if (has_values) {
    df$value <- values
  }

  # Add vertex column with requested name
  df[[vertex_column_name]] <- vertices

  # Reorder columns
  col_order <- c(vertex_column_name, "hop_idx", "is_max", "label", "spurious")
  if (has_values) {
    col_order <- c(col_order, "value")
  }

  df <- df[, col_order]

  return(df)
}


#' Analyze Extrema Distribution and Suggest Spurious Classification Parameters
#'
#' @description
#' Provides diagnostic information about the distribution of hop indices
#' to help choose appropriate parameters for spurious extrema classification.
#'
#' @param result An object from compute.extrema.hop.nbhds() or create.gflow.cx()
#' @param plot Logical; if TRUE, create diagnostic plots (default TRUE)
#'
#' @return List with diagnostic information:
#'   \itemize{
#'     \item summary: Summary statistics of hop_idx distribution
#'     \item suggested_threshold: Suggested threshold based on data
#'     \item recommended_method: Recommended classification method
#'   }
#'
#' @examples
#' \dontrun{
#' diagnostics <- analyze.extrema.distribution(result)
#' print(diagnostics$summary)
#' }
#'
#' @export
analyze.extrema.distribution <- function(result, plot = TRUE) {

  # Extract hop indices
  lmin_hop <- sapply(result$lmin_hop_nbhds, function(x) x$hop_idx)
  lmax_hop <- sapply(result$lmax_hop_nbhds, function(x) x$hop_idx)
  all_hop <- c(lmin_hop, lmax_hop)

  # Remove infinite values for analysis
  finite_hop <- all_hop[is.finite(all_hop)]

  if (length(finite_hop) == 0) {
    return(list(
      summary = "All extrema are global (hop_idx = Inf)",
      suggested_threshold = Inf,
      recommended_method = "threshold"
    ))
  }

  # Compute statistics
  hop_summary <- summary(finite_hop)
  q25 <- quantile(finite_hop, 0.25)
  q75 <- quantile(finite_hop, 0.75)
  iqr <- q75 - q25

  # Suggest threshold and method
  if (iqr < 1) {
    recommended_method <- "quantile"
    suggested_threshold <- q25
    rationale <- "Low IQR suggests using quantile-based method"
  } else if (max(finite_hop) > 10) {
    recommended_method <- "combined"
    suggested_threshold <- min(3, q25)
    rationale <- "Wide range suggests combined absolute/relative criteria"
  } else {
    recommended_method <- "threshold"
    suggested_threshold <- max(1, floor(q25))
    rationale <- "Moderate distribution suitable for simple threshold"
  }

  # Create diagnostic plot
  if (plot) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # Histogram
    hist(finite_hop, breaks = max(10, length(unique(finite_hop))),
         main = "Distribution of hop_idx",
         xlab = "hop_idx", col = "lightblue", las = 1)
    abline(v = suggested_threshold, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Suggested:", suggested_threshold),
           col = "red", lty = 2, lwd = 2)

    # Boxplot
    boxplot(finite_hop, horizontal = TRUE, col = "lightgreen",
            main = "Boxplot of hop_idx", xlab = "hop_idx", las = 1)
    abline(v = suggested_threshold, col = "red", lwd = 2, lty = 2)
  }

  # Print summary
  cat("\n=== Extrema Distribution Analysis ===\n\n")
  cat("Total extrema:", length(all_hop), "\n")
  cat("  Minima:", length(lmin_hop), "\n")
  cat("  Maxima:", length(lmax_hop), "\n")
  cat("  Global extrema (hop_idx=Inf):", sum(is.infinite(all_hop)), "\n\n")

  cat("Finite hop_idx distribution:\n")
  print(hop_summary)
  cat("\nIQR:", iqr, "\n\n")

  cat("RECOMMENDATIONS:\n")
  cat("  Method:", recommended_method, "\n")
  cat("  Threshold:", suggested_threshold, "\n")
  cat("  Rationale:", rationale, "\n\n")

  invisible(list(
    summary = hop_summary,
    iqr = iqr,
    n_global = sum(is.infinite(all_hop)),
    suggested_threshold = suggested_threshold,
    recommended_method = recommended_method,
    rationale = rationale
  ))
}

### 000000000
#' Interactive Threshold Selection Based on Expected Number of Extrema
#'
#' @description
#' Helps determine appropriate hop_idx threshold by showing how many extrema
#' would be retained at different threshold levels. Particularly useful when
#' you have domain knowledge about the expected number of true extrema.
#'
#' @param extrema_df Data frame from create.hop.nbhd.extrema.df with hop_idx column
#' @param expected_n_extrema Expected number of true extrema (NULL for no target)
#' @param expected_range Optional vector c(min, max) for expected range
#' @param plot Logical; create diagnostic plots (default TRUE)
#' @param exclude_global Logical; exclude global extrema (hop_idx = Inf) from counts
#'
#' @return List with suggested thresholds and diagnostic information
#'
#' @examples
#' \dontrun{
#' # You expect 5-7 true extrema total
#' suggestions <- suggest.threshold(
#'   lextr.res$extrema_df,
#'   expected_range = c(5, 7),
#'   plot = TRUE
#' )
#' }
#'
#' @export
suggest.threshold <- function(extrema_df,
                             expected_n_extrema = NULL,
                             expected_range = NULL,
                             plot = TRUE,
                             exclude_global = TRUE) {

  # Extract hop indices
  hop_idx <- extrema_df$hop_idx
  is_max <- extrema_df$is_max

  # Separate finite and infinite
  finite_hop <- hop_idx[is.finite(hop_idx)]
  n_global <- sum(is.infinite(hop_idx))

  if (length(finite_hop) == 0) {
    stop("All extrema are global (hop_idx = Inf)")
  }

  # Get unique hop_idx values
  unique_hop <- sort(unique(finite_hop))

  # Compute retention counts at each threshold
  thresholds <- unique_hop
  n_retained <- sapply(thresholds, function(t) {
    base_count <- sum(hop_idx > t)
    if (!exclude_global) base_count <- base_count + n_global
    base_count
  })

  # Compute retention by type
  n_retained_min <- sapply(thresholds, function(t) {
    sum(hop_idx > t & is_max == 0)
  })
  n_retained_max <- sapply(thresholds, function(t) {
    sum(hop_idx > t & is_max == 1)
  })

  # Create results table
  retention_df <- data.frame(
    threshold = thresholds,
    n_retained = n_retained,
    n_minima = n_retained_min,
    n_maxima = n_retained_max,
    pct_retained = round(100 * n_retained / length(hop_idx), 1)
  )

  # Find suggested thresholds based on expected counts
  suggestions <- list()

  if (!is.null(expected_n_extrema)) {
    # Find threshold that gives closest to expected count
    idx <- which.min(abs(n_retained - expected_n_extrema))
    suggestions$threshold_for_expected <- thresholds[idx]
    suggestions$n_at_threshold <- n_retained[idx]
  }

  if (!is.null(expected_range)) {
    # Find thresholds bracketing the expected range
    in_range_idx <- which(n_retained >= expected_range[1] &
                          n_retained <= expected_range[2])
    if (length(in_range_idx) > 0) {
      suggestions$threshold_range <- range(thresholds[in_range_idx])
      suggestions$recommended_threshold <- max(thresholds[in_range_idx])
      suggestions$n_at_recommended <- n_retained[max(in_range_idx)]
    } else {
      # Find closest
      below_range <- which(n_retained < expected_range[1])
      above_range <- which(n_retained > expected_range[2])

      if (length(below_range) > 0) {
        suggestions$threshold_below_range <- thresholds[max(below_range)]
        suggestions$n_below <- n_retained[max(below_range)]
      }
      if (length(above_range) > 0) {
        suggestions$threshold_above_range <- thresholds[min(above_range)]
        suggestions$n_above <- n_retained[min(above_range)]
      }
    }
  }

  # Statistical suggestions
  q75 <- quantile(finite_hop, 0.75)
  q90 <- quantile(finite_hop, 0.90)

  suggestions$q75_threshold <- q75
  suggestions$q90_threshold <- q90
  suggestions$n_above_q75 <- sum(hop_idx > q75)
  suggestions$n_above_q90 <- sum(hop_idx > q90)

  # Print summary
  cat("\n=== THRESHOLD SELECTION ANALYSIS ===\n\n")
  cat("Total extrema:", length(hop_idx), "\n")
  cat("  Minima:", sum(is_max == 0), "\n")
  cat("  Maxima:", sum(is_max == 1), "\n")
  if (n_global > 0) {
    cat("  Global extrema (hop_idx=Inf):", n_global, "\n")
  }
  cat("\n")

  cat("Hop index distribution:\n")
  cat("  Range: [", min(finite_hop), ",", max(finite_hop), "]\n")
  cat("  Median:", median(finite_hop), "\n")
  cat("  Q75:", q75, "\n")
  cat("  Q90:", q90, "\n\n")

  cat("RETENTION TABLE (extrema with hop_idx > threshold):\n")
  print(retention_df, row.names = FALSE)
  cat("\n")

  if (!is.null(expected_range) && !is.null(suggestions$recommended_threshold)) {
    cat("RECOMMENDATION:\n")
    cat("  To retain", expected_range[1], "-", expected_range[2], "extrema:\n")
    cat("  Use threshold =", suggestions$recommended_threshold, "\n")
    cat("  This retains:", suggestions$n_at_recommended, "extrema\n")
    cat("    (", n_retained_min[thresholds == suggestions$recommended_threshold],
        "minima,",
        n_retained_max[thresholds == suggestions$recommended_threshold],
        "maxima )\n\n")
  } else if (!is.null(expected_range)) {
    cat("NOTE: No threshold gives exactly", expected_range[1], "-",
        expected_range[2], "extrema\n")
    if (!is.null(suggestions$threshold_below_range)) {
      cat("  Threshold", suggestions$threshold_below_range,
          "gives", suggestions$n_below, "extrema (below range)\n")
    }
    if (!is.null(suggestions$threshold_above_range)) {
      cat("  Threshold", suggestions$threshold_above_range,
          "gives", suggestions$n_above, "extrema (above range)\n")
    }
    cat("\n")
  }

  # Create plots
  if (plot) {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

    # Plot 1: Histogram of hop_idx
    hist(finite_hop, breaks = max(20, length(unique(finite_hop))),
         main = "Distribution of hop_idx",
         xlab = "hop_idx", col = "lightblue", las = 1)
    if (!is.null(suggestions$recommended_threshold)) {
      abline(v = suggestions$recommended_threshold, col = "red", lwd = 2, lty = 2)
      text(suggestions$recommended_threshold, par("usr")[4] * 0.9,
           paste("threshold =", suggestions$recommended_threshold),
           pos = 4, col = "red", cex = 0.9)
    }

    # Plot 2: Retention curve
    plot(thresholds, n_retained, type = "b", pch = 19,
         main = "Extrema Retained vs Threshold",
         xlab = "hop_idx threshold", ylab = "N extrema retained",
         las = 1, col = "blue")
    grid()
    if (!is.null(expected_range)) {
      abline(h = expected_range, col = "green", lty = 2, lwd = 2)
      if (!is.null(suggestions$recommended_threshold)) {
        abline(v = suggestions$recommended_threshold, col = "red", lwd = 2, lty = 2)
        points(suggestions$recommended_threshold, suggestions$n_at_recommended,
               pch = 19, col = "red", cex = 1.5)
      }
    }

    # Plot 3: Minima vs Maxima retention
    plot(thresholds, n_retained_min, type = "b", pch = 19,
         main = "Retention by Type",
         xlab = "hop_idx threshold", ylab = "N extrema retained",
         las = 1, col = "blue", ylim = range(c(n_retained_min, n_retained_max)))
    lines(thresholds, n_retained_max, type = "b", pch = 17, col = "red")
    legend("topright", legend = c("Minima", "Maxima"),
           col = c("blue", "red"), pch = c(19, 17), lty = 1)
    grid()
    if (!is.null(suggestions$recommended_threshold)) {
      abline(v = suggestions$recommended_threshold, col = "darkgreen", lwd = 2, lty = 2)
    }

    # Plot 4: Cumulative distribution
    ecdf_hop <- ecdf(finite_hop)
    plot(ecdf_hop, main = "Cumulative Distribution of hop_idx",
         xlab = "hop_idx", ylab = "Cumulative proportion",
         las = 1, col = "blue", lwd = 2)
    grid()
    abline(h = c(0.75, 0.90), col = "gray", lty = 2)
    abline(v = c(q75, q90), col = "gray", lty = 2)
    if (!is.null(suggestions$recommended_threshold)) {
      abline(v = suggestions$recommended_threshold, col = "red", lwd = 2, lty = 2)
    }
  }

  # Return results
  result <- list(
    retention_table = retention_df,
    suggestions = suggestions,
    hop_idx_summary = summary(finite_hop),
    n_global = n_global
  )

  invisible(result)
}


#' Compare Function Values of Retained Extrema at Different Thresholds
#'
#' @description
#' Shows the actual extrema (with their function values) that would be retained
#' at different threshold choices. Helps assess whether retained extrema make
#' biological/domain sense.
#'
#' @param extrema_df Data frame with vertex, hop_idx, value, label columns
#' @param thresholds Vector of thresholds to compare
#' @param sort_by Character; how to sort results ("hop_idx", "value", or "vertex")
#'
#' @return List of data frames, one per threshold
#'
#' @examples
#' \dontrun{
#' # Compare what you'd retain at thresholds 2, 4, 6
#' retained <- compare.retained.extrema(
#'   lextr.res$extrema_df,
#'   thresholds = c(2, 4, 6),
#'   sort_by = "value"
#' )
#' }
#'
#' @export compare.retained.extrema
compare.retained.extrema <- function(extrema_df,
                                    thresholds = c(2, 3, 4, 5),
                                    sort_by = "hop_idx") {

  results <- list()

  for (t in thresholds) {
    retained <- extrema_df[extrema_df$hop_idx > t, ]

    # Sort according to preference
    if (sort_by == "hop_idx") {
      retained <- retained[order(-retained$hop_idx), ]
    } else if (sort_by == "value" && "value" %in% names(retained)) {
      retained <- retained[order(retained$value), ]
    } else if (sort_by == "vertex") {
      retained <- retained[order(retained$vertex), ]
    }

    cat("\n========================================\n")
    cat("THRESHOLD =", t, " (retain", nrow(retained), "extrema)\n")
    cat("========================================\n")

    if (nrow(retained) > 0) {
      # Print summary by type
      n_min <- sum(retained$is_max == 0)
      n_max <- sum(retained$is_max == 1)
      cat("  Minima:", n_min, "\n")
      cat("  Maxima:", n_max, "\n\n")

      # Print the extrema
      print(retained, row.names = FALSE)
    } else {
      cat("  (no extrema retained)\n")
    }

    results[[paste0("threshold_", t)]] <- retained
  }

  invisible(results)
}


#' Visualize Extrema on hop_idx vs value plot
#'
#' @description
#' Creates a scatter plot showing hop_idx vs function value for all extrema,
#' colored by type and marked by spurious classification. Helps visualize
#' the relationship between stability (hop_idx) and extremal values.
#'
#' @param extrema_df Data frame with hop_idx, value, is_max columns
#' @param threshold Optional threshold to draw as reference line
#' @param main Plot title
#'
#' @export
plot_extrema_stability <- function(extrema_df,
                                   threshold = NULL,
                                   main = "Extrema Stability vs Function Value") {

  if (!"value" %in% names(extrema_df)) {
    stop("extrema_df must have 'value' column")
  }

  # Separate minima and maxima
  minima <- extrema_df[extrema_df$is_max == 0, ]
  maxima <- extrema_df[extrema_df$is_max == 1, ]

  # Replace Inf with a large value for plotting
  finite_max <- max(extrema_df$hop_idx[is.finite(extrema_df$hop_idx)])
  plot_hop <- extrema_df$hop_idx
  plot_hop[is.infinite(plot_hop)] <- finite_max * 1.2

  minima_plot_hop <- plot_hop[extrema_df$is_max == 0]
  maxima_plot_hop <- plot_hop[extrema_df$is_max == 1]

  # Create plot
  plot(minima_plot_hop, minima$value,
       pch = 19, col = "blue", cex = 1.2,
       xlab = "hop_idx (Inf shown at right edge)",
       ylab = "Function value",
       main = main,
       las = 1,
       xlim = range(plot_hop),
       ylim = range(extrema_df$value))

  points(maxima_plot_hop, maxima$value,
         pch = 17, col = "red", cex = 1.2)

  grid()

  # Add threshold line if specified
  if (!is.null(threshold)) {
    abline(v = threshold, col = "darkgreen", lwd = 2, lty = 2)
    text(threshold, par("usr")[4], paste("threshold =", threshold),
         pos = 4, col = "darkgreen", cex = 0.9, offset = 0.5)
  }

  # Add legend
  legend("topleft",
         legend = c("Minima", "Maxima"),
         pch = c(19, 17),
         col = c("blue", "red"),
         pt.cex = 1.2,
         bg = "white")

  # Mark infinite hop_idx
  if (any(is.infinite(extrema_df$hop_idx))) {
    text(par("usr")[2], par("usr")[3],
         "← Global extrema",
         pos = 2, cex = 0.8, col = "gray40")
  }
}

