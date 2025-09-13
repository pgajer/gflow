#' Detect Local Extrema in a Graph
#'
#' @description
#' Identifies local maxima or minima in a graph based on vertex function values.
#' A vertex is considered a local extremum if it has the highest (for maxima) or
#' lowest (for minima) function value within a neighborhood of specified radius,
#' and the neighborhood contains at least a minimum number of vertices.
#'
#' @param adj.list A list where each element contains integer indices of vertices
#'   adjacent to the corresponding vertex. Must have length equal to the number
#'   of vertices in the graph.
#' @param weight.list A list where each element contains numeric weights of edges
#'   from the corresponding vertex. \code{weight.list[[i]][j]} is the weight of
#'   the edge from vertex \code{i} to vertex \code{adj.list[[i]][j]}.
#' @param y A numeric vector of function values at each vertex. Must have the
#'   same length as \code{adj.list}.
#' @param max.radius Positive numeric value specifying the maximum radius for
#'   neighborhood search.
#' @param min.neighborhood.size Positive integer specifying the minimum number
#'   of vertices required in a neighborhood for a vertex to be considered an extremum.
#' @param detect.maxima Logical; if \code{TRUE} (default), detect local maxima;
#'   if \code{FALSE}, detect local minima.
#' @param custom.prefix Character string to use as prefix for extrema labels.
#'   If \code{NULL} (default), uses "M" for maxima and "m" for minima.
#'
#' @return An object of class \code{"local_extrema"}, which is a list containing:
#'   \describe{
#'     \item{vertices}{Integer vector of vertex indices identified as extrema}
#'     \item{values}{Numeric vector of function values at the extrema}
#'     \item{radii}{Numeric vector of neighborhood radii where extremum property holds}
#'     \item{neighborhood_sizes}{Integer vector of the number of vertices in each extremum's neighborhood}
#'     \item{is_maxima}{Logical vector indicating whether each extremum is a maximum (\code{TRUE}) or minimum (\code{FALSE})}
#'     \item{type}{Character vector with values "Maximum" or "Minimum" for each extremum}
#'     \item{labels}{Character vector of labels for each extremum (e.g., "M1", "M2" for maxima)}
#'   }
#'
#' @details
#' The algorithm uses a graph-based approach to identify local extrema by examining
#' neighborhoods defined by graph distance. For each vertex, it searches within
#' increasing radii up to \code{max.radius} to find a neighborhood where the vertex
#' has the extreme value among all vertices in that neighborhood.
#'
#' The implementation uses a C++ backend for computational efficiency, particularly
#' beneficial for large graphs.
#'
#' @examples
#' # Create a simple chain graph
#' adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
#' weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
#' y <- c(1, 3, 2, 5, 1)  # Function values with peaks at vertices 2 and 4
#'
#' # Detect maxima
#' maxima <- detect.local.extrema(adj.list, weight.list, y,
#'                                max.radius = 2,
#'                                min.neighborhood.size = 2)
#' print(maxima$vertices)  # Should identify vertices 2 and 4
#'
#' # Detect minima
#' minima <- detect.local.extrema(adj.list, weight.list, y,
#'                                max.radius = 2,
#'                                min.neighborhood.size = 2,
#'                                detect.maxima = FALSE)
#' print(minima$vertices)  # Should identify vertices 1, 3, and 5
#'
#' @seealso
#' \code{\link{summary.local_extrema}} for summarizing results,
#' \code{\link{plot.local_extrema}} for visualization
#'
#' @export
detect.local.extrema <- function(adj.list,
                                 weight.list,
                                 y,
                                 max.radius,
                                 min.neighborhood.size,
                                 detect.maxima = TRUE,
                                 custom.prefix = NULL) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }
    if (!is.list(weight.list)) {
        stop("'weight.list' must be a list")
    }
    if (length(adj.list) != length(weight.list)) {
        stop("'adj.list' and 'weight.list' must have the same length")
    }
    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop("'y' must be a numeric vector with length equal to the number of vertices")
    }
    if (!is.numeric(max.radius) || length(max.radius) != 1 || max.radius <= 0) {
        stop("'max.radius' must be a positive numeric value")
    }
    if (!is.numeric(min.neighborhood.size) || length(min.neighborhood.size) != 1 ||
        min.neighborhood.size < 1 || min.neighborhood.size != floor(min.neighborhood.size)) {
        stop("'min.neighborhood.size' must be a positive integer")
    }
    if (!is.logical(detect.maxima) || length(detect.maxima) != 1) {
        stop("'detect.maxima' must be a single logical value")
    }
    if (!is.null(custom.prefix) && (!is.character(custom.prefix) || length(custom.prefix) != 1)) {
        stop("'custom.prefix' must be NULL or a single character string")
    }

    # Convert to 0-based indexing for C++
    adj.list.0based <- lapply(adj.list, function(x) as.integer(x - 1))

    # Call the C++ function
    result <- .Call(S_detect_local_extrema,
                    adj.list.0based,
                    weight.list,
                    y,
                    as.double(max.radius),
                    as.integer(min.neighborhood.size),
                    as.logical(detect.maxima))
    
    # Add extrema type label
    result$type <- ifelse(result$is_maxima, "Maximum", "Minimum")

    # Add labels based on function values
    if (length(result$vertices) > 0) {
        # Determine prefix for labels
        if (!is.null(custom.prefix)) {
            prefix <- custom.prefix
        } else {
            prefix <- ifelse(detect.maxima, "M", "m")
        }

        # Order based on function values
        if (detect.maxima) {
            value_order <- order(result$values, decreasing = TRUE)
        } else {
            value_order <- order(result$values, decreasing = FALSE)
        }

        # Assign labels in order
        result$labels <- character(length(result$vertices))
        result$labels[value_order] <- paste0(prefix, seq_along(value_order))
    } else {
        result$labels <- character(0)
    }

    class(result) <- "local_extrema"
    return(result)
}

#' Summarize Local Extrema Detection Results
#'
#' @description
#' Provides a comprehensive summary of local extrema detection results, including
#' statistics on the number of extrema found, their function values, neighborhood
#' sizes, and radii.
#'
#' @param object An object of class \code{"local_extrema"}, as returned by
#'   \code{\link{detect.local.extrema}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{"summary.local_extrema"}, which is a list containing:
#'   \describe{
#'     \item{n_extrema}{Integer; total number of extrema found}
#'     \item{extrema_type}{Character; type of extrema ("Maximum" or "Minimum")}
#'     \item{fn_values_summary}{Named numeric vector; summary statistics of function values at extrema}
#'     \item{neighborhood_sizes_summary}{Named numeric vector; summary statistics of neighborhood sizes}
#'     \item{radius_summary}{Named numeric vector; summary statistics of neighborhood radii}
#'     \item{extrema_details}{Data frame with detailed information about each extremum}
#'   }
#'
#' @examples
#' # Create example data
#' adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
#' weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
#' y <- c(1, 3, 2, 5, 1)
#'
#' # Detect and summarize maxima
#' maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)
#' summary(maxima)
#'
#' @method summary local_extrema
#' @export
summary.local_extrema <- function(object, ...) {
    if (!inherits(object, "local_extrema")) {
        stop("Object must be of class 'local_extrema'")
    }

    result <- list()
    n_extrema <- length(object$vertices)
    result$n_extrema <- n_extrema

    # Determine type of extrema
    if (n_extrema > 0) {
        result$extrema_type <- unique(object$type)[1]
    } else {
        result$extrema_type <- ifelse(all(object$is_maxima), "Maximum", "Minimum")
    }

    # Handle case with no extrema
    if (n_extrema == 0) {
        result$fn_values_summary <- NA
        result$neighborhood_sizes_summary <- NA
        result$radius_summary <- NA
        result$extrema_details <- data.frame()
        class(result) <- "summary.local_extrema"
        return(result)
    }

    # Calculate summary statistics
    result$fn_values_summary <- summary(object$values)
    result$neighborhood_sizes_summary <- summary(object$neighborhood_sizes)
    result$radius_summary <- summary(object$radii)

    # Create detailed data frame
    extrema_df <- data.frame(
        label = object$labels,
        vertex_index = object$vertices,
        fn_value = object$values,
        radius = object$radii,
        neighborhood_size = object$neighborhood_sizes,
        stringsAsFactors = FALSE
    )

    # Order by function value
    if (result$extrema_type == "Maximum") {
        extrema_df <- extrema_df[order(extrema_df$fn_value, decreasing = TRUE), ]
    } else {
        extrema_df <- extrema_df[order(extrema_df$fn_value), ]
    }

    result$extrema_details <- extrema_df
    class(result) <- "summary.local_extrema"
    return(result)
}

#' Print Summary of Local Extrema Detection Results
#'
#' @description
#' Prints a formatted summary of local extrema detection results.
#'
#' @param x An object of class \code{"summary.local_extrema"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.local_extrema
#' @export
print.summary.local_extrema <- function(x, ...) {
    cat("Local Extrema Detection Summary\n")
    cat("==============================\n\n")

    cat("Extrema type:", x$extrema_type, "\n")
    cat("Number of extrema found:", x$n_extrema, "\n\n")

    if (x$n_extrema == 0) {
        cat("No extrema were found with the specified parameters.\n")
        return(invisible(x))
    }

    cat("Function values statistics:\n")
    print(x$fn_values_summary)
    cat("\n")

    cat("Neighborhood sizes statistics:\n")
    print(x$neighborhood_sizes_summary)
    cat("\n")

    cat("Neighborhood radii statistics:\n")
    print(x$radius_summary)
    cat("\n")

    cat("Top extrema by function value:\n")
    print(head(x$extrema_details, 10))

    if (nrow(x$extrema_details) > 10) {
        cat("\n(Showing first 10 out of", nrow(x$extrema_details), "extrema)\n")
    }

    invisible(x)
}

#' Extract Vertices from a Local Extrema Object
#'
#' @description
#' Generic function to extract vertices from various objects.
#'
#' @param object An object from which to extract vertices.
#' @param ... Additional arguments passed to methods.
#'
#' @return The extracted vertices (format depends on the method).
#'
#' @export
vertices <- function(object, ...) {
    UseMethod("vertices")
}

#' Extract Vertices of a Specific Local Extremum
#'
#' @description
#' Extracts the vertices belonging to a specific local extremum identified by its label.
#' Returns all vertices in the neighborhood of the specified extremum.
#'
#' @param object An object of class \code{"local_extrema"}.
#' @param label Character string specifying the label of the extremum (e.g., "M1", "m2").
#' @param include.center Logical; if \code{TRUE} (default), include the center vertex.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector of vertex indices in the extremum's neighborhood.
#'
#' @note This method requires that the C++ backend provides neighborhood vertex
#'   information in the \code{neighborhood_vertices} component of the result.
#'
#' @examples
#' # Create example data
#' adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
#' weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
#' y <- c(1, 3, 2, 5, 1)
#'
#' # Detect maxima
#' maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)
#'
#' # Extract vertices for the first maximum (if it exists)
#' if (length(maxima$vertices) > 0) {
#'   v <- vertices(maxima, maxima$labels[1])
#' }
#'
#' @method vertices local_extrema
#' @export
vertices.local_extrema <- function(object, label, include.center = TRUE, ...) {
    if (!inherits(object, "local_extrema")) {
        stop("Object must be of class 'local_extrema'")
    }

    if (missing(label)) {
        stop("'label' argument is required")
    }

    if (!is.character(label) || length(label) != 1) {
        stop("'label' must be a single character string")
    }

    # Find the index of the extremum with the specified label
    label_index <- which(object$labels == label)

    if (length(label_index) == 0) {
        stop(paste("Label '", label, "' not found in the local extrema object", sep = ""))
    }

    # Check if neighborhood_vertices component exists
    if (is.null(object$neighborhood_vertices)) {
        stop("Neighborhood vertices information not available in this object")
    }

    # Extract vertices for the specified extremum
    vertices <- object$neighborhood_vertices[[label_index]]

    # Remove center vertex if requested
    if (!include.center) {
        center_vertex <- object$vertices[label_index]
        vertices <- vertices[vertices != center_vertex]
    }

    return(vertices)
}

#' Plot Local Extrema Detection Results
#'
#' @description
#' Creates visualizations showing the distribution of function values and
#' neighborhood characteristics of detected extrema.
#'
#' @param x An object of class \code{"local_extrema"}.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @details
#' This function creates a two-panel plot:
#' \itemize{
#'   \item Left panel: Function values at extrema, ordered by rank
#'   \item Right panel: Neighborhood radii for each extremum
#' }
#' Point sizes in both panels are proportional to neighborhood sizes.
#'
#' @examples
#' # Create example data
#' adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
#' weight.list <- list(c(1), c(1,1), c(1,1), c(1,1), c(1))
#' y <- c(1, 3, 2, 5, 1)
#'
#' # Detect and plot maxima
#' maxima <- detect.local.extrema(adj.list, weight.list, y, 2, 2)
#' plot(maxima)
#'
#' @method plot local_extrema
#' @export
plot.local_extrema <- function(x, ...) {
    if (!inherits(x, "local_extrema")) {
        stop("Object must be of class 'local_extrema'")
    }

    if (length(x$vertices) == 0) {
        message("No extrema to plot")
        return(invisible(x))
    }

    # Determine if we're dealing with maxima or minima
    detect.maxima <- all(x$is_maxima)

    # Sort extrema by function value
    idx <- order(x$values, decreasing = detect.maxima)
    sorted_values <- x$values[idx]
    sorted_radii <- x$radii[idx]
    sorted_sizes <- x$neighborhood_sizes[idx]
    sorted_labels <- x$labels[idx]

    # Save current par settings
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # Create two-panel plot
    par(mfrow = c(1, 2))

    # Plot 1: Function values
    plot(seq_along(sorted_values), sorted_values,
         main = paste("Local", ifelse(detect.maxima, "Maxima", "Minima"), "Values"),
         xlab = "Rank",
         ylab = "Function Value",
         pch = 19,
         col = "blue",
         cex = sqrt(sorted_sizes) / max(sqrt(sorted_sizes)) * 2,
         ...)

    # Add labels to top points
    n_labels <- min(5, length(sorted_values))
    if (n_labels > 0) {
        text(1:n_labels, sorted_values[1:n_labels],
             labels = sorted_labels[1:n_labels],
             pos = 4,
             offset = 0.5,
             cex = 0.8)
    }

    # Add legend
    legend("topright",
           legend = c("Point size proportional to", "neighborhood size"),
           bty = "n",
           cex = 0.9)

    # Plot 2: Neighborhood radii
    plot(seq_along(sorted_radii), sorted_radii,
         main = "Neighborhood Radii",
         xlab = "Rank",
         ylab = "Radius",
         pch = 19,
         col = "darkgreen",
         cex = sqrt(sorted_sizes) / max(sqrt(sorted_sizes)) * 2,
         ...)

    # Add labels to top points
    if (n_labels > 0) {
        text(1:n_labels, sorted_radii[1:n_labels],
             labels = sorted_labels[1:n_labels],
             pos = 4,
             offset = 0.5,
             cex = 0.8)
    }

    invisible(x)
}

#' Compute Persistent Homology for Graphs
#'
#' @description
#' Computes 0-dimensional persistent homology for a graph based on vertex function
#' values. This tracks the birth and death of connected components during a
#' filtration process.
#'
#' @param adj.list A list where each element contains the indices of adjacent vertices.
#' @param y A numeric vector of function values at vertices.
#' @param n.top.basins Integer; number of top persistence basins to return (default: 4).
#' @param alpha Numeric between 0 and 1; weight parameter for weighted persistence
#'   (default: 0.5). Only used when \code{weighted = TRUE}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{persistence}{Data frame with persistence information for each component}
#'     \item{basins}{List of vertex sets for the top persistence basins}
#'   }
#'
#' @details
#' For unweighted graphs, the algorithm processes vertices in decreasing order of
#' function values, tracking when connected components merge. For weighted graphs,
#' the death values incorporate both vertex values and edge weights.
#'
#' The persistence of a component quantifies its significance, with higher values
#' indicating more prominent topological features.
#'
#' @examples
#' # Simple chain graph
#' adj.list <- list(c(2), c(1,3), c(2,4), c(3,5), c(4))
#' y <- c(10, 7, 9, 2, 5)
#'
#' # Compute persistence
#' result <- compute.persistence(adj.list, y)
#' print(result$persistence)
#'
#' @export
compute.persistence <- function(adj.list, y, n.top.basins = 4, alpha = 0.5) {
    if (!is.list(adj.list)) {
        stop("'adj.list' must be a list")
    }
    if (!is.numeric(y) || length(y) != length(adj.list)) {
        stop("'y' must be a numeric vector with length equal to the number of vertices")
    }
    if (!is.numeric(n.top.basins) || length(n.top.basins) != 1 || n.top.basins < 1) {
        stop("'n.top.basins' must be a positive integer")
    }
    if ((!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1)) {
        stop("'alpha' must be a single numeric value between 0 and 1")
    }

    unweighted.compute.persistence(adj.list, y, as.integer(n.top.basins))
}

# Internal function for unweighted persistence (not exported)
unweighted.compute.persistence <- function(adj.list, y, n.top.basins = 4) {
    n.vertices <- length(y)

    # Sort vertices by function value (decreasing)
    sorted.vertices <- order(y, decreasing = TRUE)

    # Initialize union-find structure
    parent <- seq_len(n.vertices)

    # Initialize component tracking
    components <- data.frame(
        creator = integer(0),
        birth.value = numeric(0),
        death.value = numeric(0),
        merged.at = integer(0),
        persistence = numeric(0),
        stringsAsFactors = FALSE
    )
    components$basin <- vector("list", nrow(components))

    # Track current vertex being processed
    current.vertex <- 0

    # Mapping from union-find representative to row index
    comp.row <- list()

    # Find root with path compression
    find.root <- function(x) {
        if (parent[x] != x) {
            parent[x] <<- find.root(parent[x])
        }
        return(parent[x])
    }

    # Union function
    union <- function(x, other) {
        root.x <- find.root(x)
        root.other <- find.root(other)
        if (root.x == root.other) return(NULL)

        # Determine which component dies
        is.x.less <- y[root.x] < y[root.other]
        dying.root <- if (is.x.less) root.x else root.other
        surviving.root <- if (is.x.less) root.other else root.x

        key.dying <- as.character(dying.root)
        key.surviving <- as.character(surviving.root)

        if (!(key.dying %in% names(comp.row))) {
            warning("Dying component not found in mapping.")
            return(NULL)
        }
        dying.row <- comp.row[[key.dying]]

        # Record merge event
        components$death.value[dying.row] <<- y[current.vertex]
        components$merged.at[dying.row] <<- current.vertex
        components$persistence[dying.row] <<- components$birth.value[dying.row] - y[current.vertex]

        # Update parent pointer
        parent[dying.root] <<- surviving.root

        # Update surviving component's basin
        if (key.surviving %in% names(comp.row)) {
            surviving.row <- comp.row[[key.surviving]]
            new.basin <- unique(c(components$basin[[surviving.row]], components$basin[[dying.row]]))
            components$basin[[surviving.row]] <<- new.basin
        }

        # Remove dying component from mapping
        comp.row[[key.dying]] <<- NULL
        return(list(dying = dying.root, surviving = surviving.root))
    }

    # Process vertices
    for (i in 1:n.vertices) {
        current.vertex <- sorted.vertices[i]

        # Get processed neighbors
        neighbors <- adj.list[[current.vertex]]
        processed.neighbors <- integer(0)
        for (neighbor in neighbors) {
            if (neighbor <= n.vertices) {
                neighbor.pos <- which(sorted.vertices == neighbor)
                if (length(neighbor.pos) > 0 && neighbor.pos < i) {
                    processed.neighbors <- c(processed.neighbors, neighbor)
                }
            }
        }

        if (length(processed.neighbors) == 0) {
            # New local maximum
            new.row <- nrow(components) + 1
            new.entry <- data.frame(
                creator = current.vertex,
                birth.value = y[current.vertex],
                death.value = NA,
                merged.at = NA,
                persistence = NA,
                stringsAsFactors = FALSE
            )
            new.entry$basin <- list(current.vertex)
            components <- rbind(components, new.entry)
            comp.row[[as.character(current.vertex)]] <- new.row
        } else {
            # Connect to existing components
            neighbor.components <- sapply(processed.neighbors, find.root)
            unique.components <- unique(neighbor.components)

            if (length(unique.components) == 1) {
                # Attach to single component
                rep <- unique.components[1]
                parent[current.vertex] <- rep
                key <- as.character(rep)
                if (!is.null(comp.row[[key]])) {
                    row.idx <- comp.row[[key]]
                    components$basin[[row.idx]] <- unique(c(components$basin[[row.idx]], current.vertex))
                }
            } else {
                # Merge multiple components
                sorted.indices <- order(y[unique.components], decreasing = TRUE)
                sorted.components <- unique.components[sorted.indices]
                surviving.component <- sorted.components[1]
                parent[current.vertex] <- surviving.component
                key <- as.character(surviving.component)
                if (!is.null(comp.row[[key]])) {
                    row.idx <- comp.row[[key]]
                    components$basin[[row.idx]] <- unique(c(components$basin[[row.idx]], current.vertex))
                }
                # Merge other components
                for (comp in sorted.components[-1]) {
                    union(comp, surviving.component)
                }
            }
        }
    }

    # Handle global maximum
    global.max.idx <- which.max(y)
    global.max.root <- find.root(global.max.idx)
    global.key <- as.character(global.max.root)
    if (global.key %in% names(comp.row)) {
        global.row <- comp.row[[global.key]]
        components$death.value[global.row] <- min(y)
        components$persistence[global.row] <- components$birth.value[global.row] - min(y)
    }

    # Sort by persistence and assign ranks
    if (nrow(components) > 0) {
        sorted.components <- components[order(-components$persistence), ]
        sorted.components$rank <- 1:nrow(sorted.components)
        sorted.components <- sorted.components[, c("rank", "creator", "birth.value",
                                                   "death.value", "merged.at",
                                                   "persistence", "basin")]
        components <- sorted.components
    } else {
        components <- data.frame(
            rank = integer(0),
            creator = integer(0),
            birth.value = numeric(0),
            death.value = numeric(0),
            merged.at = integer(0),
            persistence = numeric(0),
            basin = I(list())
        )
    }

    # Extract top basins
    top.basins <- list()
    if (nrow(components) > 0) {
        for (i in 1:min(n.top.basins, nrow(components))) {
            top.basins[[i]] <- components$basin[[i]]
        }
    }

    return(list(persistence = components, basins = top.basins))
}
