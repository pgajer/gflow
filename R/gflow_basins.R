#' Find Gradient-Flow Basins on a Weighted Graph
#'
#' Identifies local minima and maxima of a scalar field on a graph
#' and computes their basins via geodesic flow analysis.
#'
#' @param adj.list A list of integer vectors. Each vector contains indices of vertices
#'   adjacent to the corresponding vertex. Indices must be 1-based.
#' @param weight.list A list of numeric vectors. Each vector contains weights of edges
#'   corresponding to adjacencies in \code{adj.list}. Must have the same structure
#'   as \code{adj.list}.
#' @param y A numeric vector of length \eqn{n}, the response values at each graph vertex.
#' @param min.basin.size An integer scalar, the minimum number of vertices a basin must
#'   contain to be considered significant. Default is 1.
#' @param min.path.size An integer scalar, the minimum number of vertices along a
#'   shortest path for which monotonicity is tested. Default is 2.
#' @param q.edge.thld A numeric scalar in \eqn{[0,1]}, the quantile threshold for maximum
#'   allowed edge weight along geodesics. Default is 0.5.
#'
#' @return An object of class \code{"gflow_basins"}, which is a list with the following components:
#' \describe{
#'   \item{basins}{A list with two elements:
#'     \describe{
#'       \item{ascending}{A list of local-minimum basins}
#'       \item{descending}{A list of local-maximum basins}
#'     }
#'     Each basin is a list with components:
#'     \describe{
#'       \item{vertex}{Integer, the 1-based vertex index of the extremum}
#'       \item{value}{Numeric, the function value at the extremum}
#'       \item{basin}{A matrix with columns "vertex" and "distance"}
#'       \item{label}{Character, a unique label for the extremum}
#'     }
#'   }
#'   \item{local_extrema}{A data frame summarizing all detected extrema with columns:
#'     \describe{
#'       \item{vertex_index}{Integer, 1-based vertex index}
#'       \item{is_maximum}{Integer, 1 if maximum, 0 if minimum}
#'       \item{label}{Character, unique label (e.g., "M1", "m1")}
#'       \item{fn_value}{Numeric, function value at extremum}
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{summary.gflow_basins}}, \code{\link{find.local.extrema}}
#'
#' @examples
#' # Create a simple graph
#' adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
#' weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
#' y <- c(1, 0, 2, 1.5)
#'
#' \donttest{
#' # Find basins (requires compiled C++ backend)
#' result <- find.gflow.basins(adj.list, weight.list, y,
#'                             min.basin.size = 1,
#'                             min.path.size = 2,
#'                             q.edge.thld = 0.5)
#' summary(result)
#' }
#'
#' @export
find.gflow.basins <- function(adj.list,
                              weight.list,
                              y,
                              min.basin.size = 1L,
                              min.path.size = 2L,
                              q.edge.thld = 0.5) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    n <- length(adj.list)

    if (length(weight.list) != n) {
        stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
    }

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    if (length(y) != n) {
        stop("'y' must have the same length as 'adj.list'", call. = FALSE)
    }

    # Validate list contents
    for (i in seq_len(n)) {
        if (!is.numeric(adj.list[[i]]) && !is.integer(adj.list[[i]])) {
            stop(sprintf("'adj.list[[%d]]' must contain numeric or integer values", i),
                 call. = FALSE)
        }
        if (!is.numeric(weight.list[[i]])) {
            stop(sprintf("'weight.list[[%d]]' must contain numeric values", i),
                 call. = FALSE)
        }
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("'adj.list[[%d]]' and 'weight.list[[%d]]' must have the same length", i, i),
                 call. = FALSE)
        }
        # Check for valid vertex indices
        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n)) {
                stop(sprintf("'adj.list[[%d]]' contains invalid vertex indices (must be between 1 and %d)", i, n),
                     call. = FALSE)
            }
        }
    }

    # Validate scalar parameters
    min.basin.size <- as.integer(min.basin.size)
    if (length(min.basin.size) != 1L || is.na(min.basin.size) || min.basin.size < 1L) {
        stop("'min.basin.size' must be a single positive integer", call. = FALSE)
    }

    min.path.size <- as.integer(min.path.size)
    if (length(min.path.size) != 1L || is.na(min.path.size) || min.path.size < 1L) {
        stop("'min.path.size' must be a single positive integer", call. = FALSE)
    }

    if (!is.numeric(q.edge.thld) || length(q.edge.thld) != 1L ||
        is.na(q.edge.thld) || q.edge.thld < 0 || q.edge.thld > 1) {
        stop("'q.edge.thld' must be a single numeric value in [0, 1]", call. = FALSE)
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) > 0) as.integer(x - 1L) else integer(0)
    })

    # Call C++ backend
    result <- .Call("S_find_gflow_basins",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    min.basin.size,
                    min.path.size,
                    as.numeric(q.edge.thld))

    # Process results
    result <- process_basin_results(result)

    class(result) <- "gflow_basins"
    return(result)
}


#' Process Basin Results from C++ Backend
#'
#' Internal function to process raw basin results from C++ and create
#' the final R object structure.
#'
#' @param result Raw result from C++ backend
#' @return Processed result list
#' @keywords internal
process_basin_results <- function(result) {
    # Extract minima and maxima information
    mins <- lapply(result$lmin_basins, function(b) {
        list(vertex = b$vertex, value = b$value)
    })
    maxs <- lapply(result$lmax_basins, function(b) {
        list(vertex = b$vertex, value = b$value)
    })

    n_min <- length(mins)
    n_max <- length(maxs)

    # Create extrema data frame
    if (n_min + n_max > 0) {
        e_idx <- c(
            if (n_min > 0) vapply(mins, `[[`, numeric(1), "vertex") else numeric(0),
            if (n_max > 0) vapply(maxs, `[[`, numeric(1), "vertex") else numeric(0)
        )

        is_max <- c(rep(FALSE, n_min), rep(TRUE, n_max))

        fn_val <- c(
            if (n_min > 0) vapply(mins, `[[`, numeric(1), "value") else numeric(0),
            if (n_max > 0) vapply(maxs, `[[`, numeric(1), "value") else numeric(0)
        )

        # Create labels
        labels <- character(length(e_idx))

        # Label maxima (M1, M2, ...) in decreasing order of function value
        if (n_max > 0) {
            max_idx <- (n_min + 1):(n_min + n_max)
            ord <- order(fn_val[max_idx], decreasing = TRUE)
            labels[max_idx[ord]] <- paste0("M", seq_len(n_max))
        }

        # Label minima (m1, m2, ...) in increasing order of function value
        if (n_min > 0) {
            min_idx <- 1:n_min
            ord <- order(fn_val[min_idx], decreasing = FALSE)
            labels[min_idx[ord]] <- paste0("m", seq_len(n_min))
        }

        local_extrema <- data.frame(
            vertex_index = as.integer(e_idx),
            is_maximum = as.integer(is_max),
            label = labels,
            fn_value = fn_val,
            stringsAsFactors = FALSE
        )
    } else {
        local_extrema <- data.frame(
            vertex_index = integer(),
            is_maximum = integer(),
            label = character(),
            fn_value = numeric(),
            stringsAsFactors = FALSE
        )
    }

    # Attach labels to basin lists
    if (n_min > 0) {
        for (i in seq_len(n_min)) {
            result$lmin_basins[[i]]$label <- local_extrema$label[
                local_extrema$vertex_index == result$lmin_basins[[i]]$vertex
            ][1]
        }
    }

    if (n_max > 0) {
        for (i in seq_len(n_max)) {
            result$lmax_basins[[i]]$label <- local_extrema$label[
                local_extrema$vertex_index == result$lmax_basins[[i]]$vertex
            ][1]
        }
    }

    # Reorganize results
    result$basins <- list(
        ascending = result$lmin_basins,
        descending = result$lmax_basins
    )
    result$local_extrema <- local_extrema
    result$lmin_basins <- NULL
    result$lmax_basins <- NULL

    return(result)
}


#' Summarize Gradient-Flow Basin Results
#'
#' @param object An object of class \code{"gflow_basins"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return An object of class \code{"summary.gflow_basins"} containing:
#' \describe{
#'   \item{n_extrema}{Total number of extrema detected}
#'   \item{n_minima}{Number of local minima}
#'   \item{n_maxima}{Number of local maxima}
#'   \item{n_ascending_basins}{Number of ascending basins}
#'   \item{n_descending_basins}{Number of descending basins}
#'   \item{ascending_basin_sizes}{Summary statistics for ascending basin sizes}
#'   \item{descending_basin_sizes}{Summary statistics for descending basin sizes}
#' }
#'
#' @export
#' @method summary gflow_basins
summary.gflow_basins <- function(object, ...) {
    if (!inherits(object, "gflow_basins")) {
        stop("'object' must be of class 'gflow_basins'", call. = FALSE)
    }

    extrema <- object$local_extrema
    n_extrema <- nrow(extrema)
    n_maxima <- sum(extrema$is_maximum == 1L)
    n_minima <- n_extrema - n_maxima

    # Get basin sizes
    asc_sizes <- if (length(object$basins$ascending) > 0) {
        vapply(object$basins$ascending, function(b) {
            if (!is.null(b$basin) && is.matrix(b$basin)) {
                nrow(b$basin)
            } else {
                0L
            }
        }, integer(1))
    } else {
        integer(0)
    }

    desc_sizes <- if (length(object$basins$descending) > 0) {
        vapply(object$basins$descending, function(b) {
            if (!is.null(b$basin) && is.matrix(b$basin)) {
                nrow(b$basin)
            } else {
                0L
            }
        }, integer(1))
    } else {
        integer(0)
    }

    structure(
        list(
            n_extrema = n_extrema,
            n_minima = n_minima,
            n_maxima = n_maxima,
            n_ascending_basins = length(object$basins$ascending),
            n_descending_basins = length(object$basins$descending),
            ascending_basin_sizes = if (length(asc_sizes) > 0) summary(asc_sizes) else NULL,
            descending_basin_sizes = if (length(desc_sizes) > 0) summary(desc_sizes) else NULL
        ),
        class = "summary.gflow_basins"
    )
}


#' Print Summary of Gradient-Flow Basin Results
#'
#' @param x An object of class \code{"summary.gflow_basins"}
#' @param ... Additional arguments (currently ignored)
#'
#' @export
#' @method print summary.gflow_basins
print.summary.gflow_basins <- function(x, ...) {
    cat("Gradient-Flow Basins Summary\n")
    cat(sprintf("Total extrema: %d (minima: %d, maxima: %d)\n",
                x$n_extrema, x$n_minima, x$n_maxima))

    if (x$n_ascending_basins > 0 && !is.null(x$ascending_basin_sizes)) {
        cat("\nAscending basins (", x$n_ascending_basins, "):\n", sep = "")
        cat("  Size summary: ")
        cat(sprintf("Min=%d, 1st Qu=%d, Median=%d, Mean=%.1f, 3rd Qu=%d, Max=%d\n",
                    x$ascending_basin_sizes["Min."],
                    x$ascending_basin_sizes["1st Qu."],
                    x$ascending_basin_sizes["Median"],
                    x$ascending_basin_sizes["Mean"],
                    x$ascending_basin_sizes["3rd Qu."],
                    x$ascending_basin_sizes["Max."]))
    } else {
        cat("\nNo ascending basins found.\n")
    }

    if (x$n_descending_basins > 0 && !is.null(x$descending_basin_sizes)) {
        cat("\nDescending basins (", x$n_descending_basins, "):\n", sep = "")
        cat("  Size summary: ")
        cat(sprintf("Min=%d, 1st Qu=%d, Median=%d, Mean=%.1f, 3rd Qu=%d, Max=%d\n",
                    x$descending_basin_sizes["Min."],
                    x$descending_basin_sizes["1st Qu."],
                    x$descending_basin_sizes["Median"],
                    x$descending_basin_sizes["Mean"],
                    x$descending_basin_sizes["3rd Qu."],
                    x$descending_basin_sizes["Max."]))
    } else {
        cat("\nNo descending basins found.\n")
    }

    invisible(x)
}


#' Detect Local Extrema on a Graph with Labeled Basins
#'
#' Identifies local minima and maxima of a scalar function defined
#' on the vertices of a graph, using a basin-based definition of extrema.
#'
#' This function uses a simpler algorithm than \code{\link{find.gflow.basins}},
#' focusing on local extrema detection without the geodesic flow analysis.
#' Each extremum is associated with a basin of attraction obtained via a
#' gradient flow process.
#'
#' @param adj.list A list of integer vectors, where each element represents the
#'   neighbors of a vertex (1-based indices). Must match the length of \code{y}.
#' @param weight.list A list of numeric vectors representing edge weights
#'   corresponding to the adjacency list structure.
#' @param y A numeric vector of function values defined on the vertices of the graph.
#' @param min.basin.size An integer specifying the minimum number of vertices
#'   required for a valid basin to be considered an extremum. Default is 1.
#'
#' @return An object of class \code{"local_extrema_basins"}, which is a list with
#'   the following components:
#' \describe{
#'   \item{basins}{A list with elements \code{ascending} (minima) and
#'     \code{descending} (maxima), each containing a list of basins with
#'     vertex index, function value, and basin information.}
#'   \item{local_extrema}{A data frame summarizing all detected extrema with columns:
#'     \code{vertex_index}, \code{is_maximum}, \code{label}, and \code{fn_value}.}
#' }
#'
#' @seealso \code{\link{find.gflow.basins}} for more comprehensive basin analysis
#'
#' @examples
#' # Create a simple graph
#' adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
#' weight.list <- list(c(1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1))
#' y <- c(1, 0, 2, 1.5)
#'
#' \donttest{
#' # Find local extrema (requires compiled C++ backend)
#' result <- find.local.extrema(adj.list, weight.list, y, min.basin.size = 1)
#' print(result$local_extrema)
#' }
#'
#' @export
find.local.extrema <- function(adj.list,
                               weight.list,
                               y,
                               min.basin.size = 1L) {

    # Input validation (similar to find.gflow.basins but simpler)
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list", call. = FALSE)
    }

    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list", call. = FALSE)
    }

    n <- length(adj.list)

    if (length(weight.list) != n) {
        stop("'adj.list' and 'weight.list' must have the same length", call. = FALSE)
    }

    if (!is.numeric(y) || !is.vector(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }

    if (length(y) != n) {
        stop("'y' must have the same length as 'adj.list'", call. = FALSE)
    }

    # Validate list contents
    for (i in seq_len(n)) {
        if (!is.numeric(adj.list[[i]]) && !is.integer(adj.list[[i]])) {
            stop(sprintf("'adj.list[[%d]]' must contain numeric or integer values", i),
                 call. = FALSE)
        }
        if (!is.numeric(weight.list[[i]])) {
            stop(sprintf("'weight.list[[%d]]' must contain numeric values", i),
                 call. = FALSE)
        }
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("'adj.list[[%d]]' and 'weight.list[[%d]]' must have the same length", i, i),
                 call. = FALSE)
        }
        # Check for valid vertex indices
        if (length(adj.list[[i]]) > 0) {
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n)) {
                stop(sprintf("'adj.list[[%d]]' contains invalid vertex indices (must be between 1 and %d)", i, n),
                     call. = FALSE)
            }
        }
    }

    # Validate min.basin.size
    min.basin.size <- as.integer(min.basin.size)
    if (length(min.basin.size) != 1L || is.na(min.basin.size) || min.basin.size < 1L) {
        stop("'min.basin.size' must be a single positive integer", call. = FALSE)
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) {
        if (length(x) > 0) as.integer(x - 1L) else integer(0)
    })

    # Call C++ backend
    result <- .Call("S_find_local_extrema",
                    adj.list.0based,
                    weight.list,
                    as.numeric(y),
                    min.basin.size)

    # Process results (reuse the same processing function)
    result <- process_basin_results(result)

    class(result) <- "local_extrema_basins"
    return(result)
}
