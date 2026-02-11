#' Compute Boundary Vertices of a Subset in a Graph
#'
#' @description
#' Computes the boundary of a subset U of vertices in a graph. A vertex v is in the boundary
#' of U if and only if v is in U and its expanded neighborhood (itself and its adjacent vertices)
#' intersects both U and the complement of U.
#'
#' @param U Integer vector of vertex indices representing a subset of the graph vertices.
#'   Must contain valid indices between 1 and the number of vertices in the graph.
#' @param adj.list List where each element contains the adjacency list for a vertex.
#'   Element \code{adj.list[[i]]} contains indices of vertices adjacent to vertex i.
#'   Self-loops are ignored.
#'
#' @return Integer vector containing the indices of vertices in the boundary of U,
#'   sorted in ascending order. Returns an empty vector if U has no boundary vertices.
#'
#' @details
#' The boundary of a subset U consists of those vertices in U that are adjacent to
#' at least one vertex outside of U. Mathematically, v is in the boundary of U if
#' and only if v is in U and there exists a vertex w not in U such that w is
#' adjacent to v.
#'
#' This definition is consistent with the topological notion of boundary adapted
#' to the discrete setting of graphs.
#'
#' @examples
#' # Create a simple graph adjacency list
#' adj.list <- list(c(2, 3), c(1, 3, 4), c(1, 2), c(2, 5), c(4))
#'
#' # Define a subset of vertices
#' U <- c(1, 2, 3)
#'
#' # Compute the boundary
#' boundary <- set.boundary(U, adj.list)
#' print(boundary)  # Returns c(2) since only vertex 2 is adjacent to vertex 4 (outside U)
#'
#' @seealso \code{\link{lmax.basins}}, \code{\link{lmin.basins}}
#' @export
set.boundary <- function(U, adj.list) {
    ## Input validation
    if (!is.list(adj.list) || length(adj.list) == 0) {
        stop("adj.list must be a non-empty list")
    }

    n_vertices <- length(adj.list)

    if (!is.numeric(U) || length(U) == 0) {
        stop("U must be a non-empty numeric vector")
    }

    ## Convert to integer and check validity
    U <- as.integer(U)
    if (!all(U %in% seq_len(n_vertices))) {
        stop("All elements of U must be valid vertex indices between 1 and ", n_vertices)
    }

    ## Remove duplicates
    U <- unique(U)

    ## Initialize boundary set
    boundary <- integer(0)

    ## Compute complement of U (V - U)
    V_minus_U <- setdiff(seq_len(n_vertices), U)

    ## If U is the entire vertex set, it has no boundary
    if (length(V_minus_U) == 0) {
        return(boundary)
    }

    ## For each vertex in U
    for (v in U) {
        ## Get neighbors, excluding self-loops
        neighbors <- setdiff(adj.list[[v]], v)

        ## Check if v has any neighbors in V-U
        if (any(neighbors %in% V_minus_U)) {
            boundary <- c(boundary, v)
        }
    }

    return(sort(unique(boundary)))
}

#' Identify Basins of Attraction of Local Maxima on a Weighted Graph
#'
#' @description
#' Computes the basins of attraction for local maxima on a weighted graph using
#' a gradient flow algorithm. Each basin consists of vertices that would flow
#' towards a particular local maximum following the steepest ascent path.
#'
#' @param object List where \code{object[[i]]} contains indices of vertices
#'   adjacent to vertex i. Must be a valid adjacency list representation.
#' @param weight.list List where \code{weight.list[[i]]} contains positive weights
#'   of edges from vertex i to corresponding vertices in \code{object[[i]]}.
#'   If \code{NULL}, uniform weights of 1 are used.
#' @param y Numeric vector of values at each vertex. Must have the same length
#'   as \code{object}. Can contain \code{NA} values which are treated as
#'   negative infinity.
#' @param lmax.list List where each element contains:
#'  \describe{
#'     \item{\code{lmax}: }{Index of the local maximum vertex}
#'     \item{\code{vertices}: }{Set of vertices forming the initial basin
#'       (usually the local maximum and its neighborhood)}
#'     \item{\code{label}: }{Character label for the basin}
#'   }
#' @param verbose Logical; if \code{TRUE}, prints progress messages during
#'   computation. Default is \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named list of basins of attraction, where names correspond to labels
#'   in \code{lmax.list}. Each element is an integer vector of vertex indices
#'   belonging to that basin, sorted in ascending order.
#'
#' @details
#' The algorithm expands basins iteratively from local maxima by examining boundary
#' vertices and adding neighboring vertices that satisfy a threshold criterion.
#' The threshold is computed as a weighted average of values at vertices already
#' in the basin.
#'
#' Edge weights influence both the threshold calculation and the decision to
#' include new vertices. Higher weights indicate stronger connections and give
#' more influence in the weighted averaging process.
#'
#' The algorithm terminates when no new vertices can be added to any basin,
#' ensuring convergence in finite graphs.
#'
#' @examples
#' \dontrun{
#' # Create a weighted graph
#' adj.list <- list(
#'   c(2, 3),      # vertex 1 neighbors
#'   c(1, 3, 4),   # vertex 2 neighbors
#'   c(1, 2, 5),   # vertex 3 neighbors
#'   c(2, 5),      # vertex 4 neighbors
#'   c(3, 4)       # vertex 5 neighbors
#' )
#'
#' # Edge weights (optional)
#' weight.list <- list(
#'   c(1.0, 0.5),     # weights for edges from vertex 1
#'   c(1.0, 0.8, 1.2), # weights for edges from vertex 2
#'   c(0.5, 0.8, 1.0), # weights for edges from vertex 3
#'   c(1.2, 1.5),      # weights for edges from vertex 4
#'   c(1.0, 1.5)       # weights for edges from vertex 5
#' )
#'
#' # Values at vertices
#' y <- c(0.5, 0.8, 0.3, 1.2, 1.0)
#'
#' # Define local maxima
#' lmax.list <- list(
#'   list(lmax = 4, vertices = c(4), label = "max1"),
#'   list(lmax = 5, vertices = c(5), label = "max2")
#' )
#'
#' # Compute basins
#' basins <- lmax.basins(adj.list, weight.list, y, lmax.list, verbose = TRUE)
#' }
#'
#' @seealso \code{\link{set.boundary}}, \code{\link{lmin.basins}},
#'   \code{\link{compute.graph.gradient.flow}}
#' @export
lmax.basins <- function(object,
                        weight.list = NULL,
                        y,
                        lmax.list,
                        verbose = FALSE,
                        ...) {
    adj.list <- object
    ## Input validation
    if (!is.list(adj.list) || length(adj.list) == 0) {
        stop("adj.list must be a non-empty list")
    }

    n_vertices <- length(adj.list)

    if (!is.numeric(y) || length(y) != n_vertices) {
        stop("y must be a numeric vector with length equal to adj.list")
    }

    if (!is.list(lmax.list) || length(lmax.list) == 0) {
        stop("lmax.list must be a non-empty list")
    }

    ## Handle NA values in y
    if (any(is.na(y))) {
        if (verbose) {
            message("Found ", sum(is.na(y)), " NA values in y, treating as -Inf")
        }
        y[is.na(y)] <- -Inf
    }

    ## Process weight list
    if (is.null(weight.list)) {
        weight.list <- lapply(adj.list, function(neighbors) rep(1, length(neighbors)))
    } else {
        if (!is.list(weight.list) || length(weight.list) != n_vertices) {
            stop("weight.list must be a list with length equal to adj.list")
        }

        ## Validate weight list structure
        for (i in seq_len(n_vertices)) {
            if (length(adj.list[[i]]) != length(weight.list[[i]])) {
                stop("Length mismatch at vertex ", i, ": adj.list has ",
                     length(adj.list[[i]]), " neighbors but weight.list has ",
                     length(weight.list[[i]]), " weights")
            }

            ## Check for positive weights
            if (any(weight.list[[i]] <= 0)) {
                stop("All weights must be positive. Found non-positive weight at vertex ", i)
            }
        }
    }

    ## Validate lmax.list structure
    required_fields <- c("lmax", "vertices", "label")
    for (i in seq_along(lmax.list)) {
        if (!all(required_fields %in% names(lmax.list[[i]]))) {
            stop("Element ", i, " of lmax.list must contain fields: ",
                 paste(required_fields, collapse = ", "))
        }

        if (!is.character(lmax.list[[i]]$label)) {
            stop("Label in element ", i, " of lmax.list must be a character string")
        }

        ## Check that vertices are valid
        vertices <- lmax.list[[i]]$vertices
        if (!all(vertices %in% seq_len(n_vertices))) {
            stop("Invalid vertex indices in element ", i, " of lmax.list")
        }
    }

    ## Initialize basins list
    basin.list <- list()

    ## Track iteration count for convergence monitoring
    max_iterations <- n_vertices  # Safety limit

    ## Process each local maximum
    for (i in seq_along(lmax.list)) {
        lmax_info <- lmax.list[[i]]
        basin_id <- lmax_info$label

        ## Initial basin
        basin <- unique(as.integer(lmax_info$vertices))
        basin_size <- length(basin)

        if (verbose) {
            message("Processing basin '", basin_id, "' starting with ",
                    basin_size, " vertices")
        }

        iteration <- 0

        ## Iterative expansion
        repeat {
            iteration <- iteration + 1

            ## Safety check
            if (iteration > max_iterations) {
                warning("Basin '", basin_id, "' reached maximum iterations (",
                        max_iterations, "). Stopping expansion.")
                break
            }

            ## Get boundary vertices
            boundary_vertices <- set.boundary(basin, adj.list)

            if (length(boundary_vertices) == 0) {
                if (verbose) {
                    message("  Basin '", basin_id, "' has no boundary vertices")
                }
                break
            }

            ## Collect candidates for addition
            candidates <- integer(0)

            ## Process each boundary vertex
            for (b in boundary_vertices) {
                ## Get neighbors and weights, excluding self-loops
                neighbors <- setdiff(adj.list[[b]], b)
                if (length(neighbors) == 0) next

                ## Get corresponding weights
                edge_indices <- which(adj.list[[b]] %in% neighbors)
                neighbor_weights <- weight.list[[b]][edge_indices]

                ## Split into neighbors inside and outside basin
                in_basin <- neighbors %in% basin
                neighbors_in <- neighbors[in_basin]
                neighbors_out <- neighbors[!in_basin]

                if (length(neighbors_in) == 0 || length(neighbors_out) == 0) {
                    next
                }

                ## Weights for neighbors in basin
                weights_in <- neighbor_weights[in_basin]

                ## Calculate weighted threshold
                threshold <- sum(y[neighbors_in] * weights_in) / sum(weights_in)

                ## Check which outside neighbors to add
                ## We add vertices with values below the threshold
                for (j in which(!in_basin)) {
                    if (y[neighbors[j]] < threshold) {
                        candidates <- c(candidates, neighbors[j])
                    }
                }
            }

            ## Add unique candidates to basin
            candidates <- unique(candidates)

            if (length(candidates) == 0) {
                if (verbose) {
                    message("  No new candidates found for basin '", basin_id, "'")
                }
                break
            }

            ## Update basin
            basin <- sort(unique(c(basin, candidates)))
            new_size <- length(basin)

            if (verbose) {
                message("  Iteration ", iteration, ": Added ",
                        length(candidates), " vertices (total: ", new_size, ")")
            }

            ## Check for convergence
            if (new_size == basin_size) {
                break
            }

            basin_size <- new_size
        }

        ## Store final basin
        basin.list[[basin_id]] <- basin

        if (verbose) {
            message("  Basin '", basin_id, "' finalized with ",
                    length(basin), " vertices")
        }
    }

    return(basin.list)
}

#' Identify Basins of Attraction of Local Minima on a Weighted Graph
#'
#' @description
#' Computes the basins of attraction for local minima on a weighted graph using
#' a gradient flow algorithm. Each basin consists of vertices that would flow
#' towards a particular local minimum following the steepest descent path.
#'
#' @param object List where \code{object[[i]]} contains indices of vertices
#'   adjacent to vertex i. Must be a valid adjacency list representation.
#' @param weight.list List where \code{weight.list[[i]]} contains positive weights
#'   of edges from vertex i to corresponding vertices in \code{object[[i]]}.
#'   If \code{NULL}, uniform weights of 1 are used.
#' @param y Numeric vector of values at each vertex. Must have the same length
#'   as \code{object}. Can contain \code{NA} values which are treated as
#'   positive infinity.
#' @param lmin.list List where each element contains:
#'   \describe{
#'     \item{\code{lmin}: }{Index of the local minimum vertex}
#'     \item{\code{vertices}: }{Set of vertices forming the initial basin
#'       (usually the local minimum and its neighborhood)}
#'     \item{\code{label}: }{Character label for the basin}
#'   }
#' @param verbose Logical; if \code{TRUE}, prints progress messages during
#'   computation. Default is \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named list of basins of attraction, where names correspond to labels
#'   in \code{lmin.list}. Each element is an integer vector of vertex indices
#'   belonging to that basin, sorted in ascending order.
#'
#' @details
#' The algorithm expands basins iteratively from local minima by examining boundary
#' vertices and adding neighboring vertices that satisfy a threshold criterion.
#' The threshold is computed as a weighted average of values at vertices already
#' in the basin.
#'
#' Edge weights influence both the threshold calculation and the decision to
#' include new vertices. Higher weights indicate stronger connections and give
#' more influence in the weighted averaging process.
#'
#' For local minima, vertices are added if their values are above the threshold,
#' representing flow towards lower values.
#'
#' @examples
#' \dontrun{
#' # Create a weighted graph
#' adj.list <- list(
#'   c(2, 3),      # vertex 1 neighbors
#'   c(1, 3, 4),   # vertex 2 neighbors
#'   c(1, 2, 5),   # vertex 3 neighbors
#'   c(2, 5),      # vertex 4 neighbors
#'   c(3, 4)       # vertex 5 neighbors
#' )
#'
#' # Values at vertices
#' y <- c(0.5, 0.8, 0.3, 1.2, 0.2)
#'
#' # Define local minima
#' lmin.list <- list(
#'   list(lmin = 5, vertices = c(5), label = "min1"),
#'   list(lmin = 3, vertices = c(3), label = "min2")
#' )
#'
#' # Compute basins
#' basins <- lmin.basins(adj.list, weight.list = NULL, y, lmin.list)
#' }
#'
#' @seealso \code{\link{set.boundary}}, \code{\link{lmax.basins}},
#'   \code{\link{compute.graph.gradient.flow}}
#' @export
lmin.basins <- function(object, weight.list = NULL, y, lmin.list, verbose = FALSE, ...) {
    adj.list <- object
    ## Input validation
    if (!is.list(adj.list) || length(adj.list) == 0) {
        stop("adj.list must be a non-empty list")
    }

    n_vertices <- length(adj.list)

    if (!is.numeric(y) || length(y) != n_vertices) {
        stop("y must be a numeric vector with length equal to adj.list")
    }

    if (!is.list(lmin.list) || length(lmin.list) == 0) {
        stop("lmin.list must be a non-empty list")
    }

    ## Handle NA values in y
    if (any(is.na(y))) {
        if (verbose) {
            message("Found ", sum(is.na(y)), " NA values in y, treating as +Inf")
        }
        y[is.na(y)] <- Inf
    }

    ## Process weight list
    if (is.null(weight.list)) {
        weight.list <- lapply(adj.list, function(neighbors) rep(1, length(neighbors)))
    } else {
        if (!is.list(weight.list) || length(weight.list) != n_vertices) {
            stop("weight.list must be a list with length equal to adj.list")
        }

        ## Validate weight list structure
        for (i in seq_len(n_vertices)) {
            if (length(adj.list[[i]]) != length(weight.list[[i]])) {
                stop("Length mismatch at vertex ", i, ": adj.list has ",
                     length(adj.list[[i]]), " neighbors but weight.list has ",
                     length(weight.list[[i]]), " weights")
            }

            ## Check for positive weights
            if (any(weight.list[[i]] <= 0)) {
                stop("All weights must be positive. Found non-positive weight at vertex ", i)
            }
        }
    }

    ## Validate lmin.list structure
    required_fields <- c("lmin", "vertices", "label")
    for (i in seq_along(lmin.list)) {
        if (!all(required_fields %in% names(lmin.list[[i]]))) {
            stop("Element ", i, " of lmin.list must contain fields: ",
                 paste(required_fields, collapse = ", "))
        }

        if (!is.character(lmin.list[[i]]$label)) {
            stop("Label in element ", i, " of lmin.list must be a character string")
        }

        ## Check that vertices are valid
        vertices <- lmin.list[[i]]$vertices
        if (!all(vertices %in% seq_len(n_vertices))) {
            stop("Invalid vertex indices in element ", i, " of lmin.list")
        }
    }

    ## Initialize basins list
    basin.list <- list()

    ## Track iteration count for convergence monitoring
    max_iterations <- n_vertices  # Safety limit

    ## Process each local minimum
    for (i in seq_along(lmin.list)) {
        lmin_info <- lmin.list[[i]]
        basin_id <- lmin_info$label

        ## Initial basin
        basin <- unique(as.integer(lmin_info$vertices))
        basin_size <- length(basin)

        if (verbose) {
            message("Processing basin '", basin_id, "' starting with ",
                    basin_size, " vertices")
        }

        iteration <- 0

        ## Iterative expansion
        repeat {
            iteration <- iteration + 1

            ## Safety check
            if (iteration > max_iterations) {
                warning("Basin '", basin_id, "' reached maximum iterations (",
                        max_iterations, "). Stopping expansion.")
                break
            }

            ## Get boundary vertices
            boundary_vertices <- set.boundary(basin, adj.list)

            if (length(boundary_vertices) == 0) {
                if (verbose) {
                    message("  Basin '", basin_id, "' has no boundary vertices")
                }
                break
            }

            ## Collect candidates for addition
            candidates <- integer(0)

            ## Process each boundary vertex
            for (b in boundary_vertices) {
                ## Get neighbors and weights, excluding self-loops
                neighbors <- setdiff(adj.list[[b]], b)
                if (length(neighbors) == 0) next

                ## Get corresponding weights
                edge_indices <- which(adj.list[[b]] %in% neighbors)
                neighbor_weights <- weight.list[[b]][edge_indices]

                ## Split into neighbors inside and outside basin
                in_basin <- neighbors %in% basin
                neighbors_in <- neighbors[in_basin]
                neighbors_out <- neighbors[!in_basin]

                if (length(neighbors_in) == 0 || length(neighbors_out) == 0) {
                    next
                }

                ## Weights for neighbors in basin
                weights_in <- neighbor_weights[in_basin]

                ## Calculate weighted threshold
                threshold <- sum(y[neighbors_in] * weights_in) / sum(weights_in)

                ## Check which outside neighbors to add
                ## For minima, we add vertices with values above the threshold
                for (j in which(!in_basin)) {
                    if (y[neighbors[j]] > threshold) {
                        candidates <- c(candidates, neighbors[j])
                    }
                }
            }

            ## Add unique candidates to basin
            candidates <- unique(candidates)

            if (length(candidates) == 0) {
                if (verbose) {
                    message("  No new candidates found for basin '", basin_id, "'")
                }
                break
            }

            ## Update basin
            basin <- sort(unique(c(basin, candidates)))
            new_size <- length(basin)

            if (verbose) {
                message("  Iteration ", iteration, ": Added ",
                        length(candidates), " vertices (total: ", new_size, ")")
            }

            ## Check for convergence
            if (new_size == basin_size) {
                break
            }

            basin_size <- new_size
        }

        ## Store final basin
        basin.list[[basin_id]] <- basin

        if (verbose) {
            message("  Basin '", basin_id, "' finalized with ",
                    length(basin), " vertices")
        }
    }

    return(basin.list)
}

#' Compute Gradient Flow Information for Vertices in a Graph
#'
#' @description
#' Analyzes the gradient flow structure of a weighted graph by computing basins
#' of attraction for local maxima and minima, determining flow directions at
#' each vertex, and identifying ambiguous vertices that belong to multiple basins.
#'
#' @param adj.list List where \code{adj.list[[i]]} contains indices of vertices
#'   adjacent to vertex i. Must be a valid adjacency list representation.
#' @param weight.list List where \code{weight.list[[i]]} contains positive weights
#'   of edges from vertex i to corresponding vertices in \code{adj.list[[i]]}.
#'   If \code{NULL}, uniform weights of 1 are used.
#' @param y Numeric vector of values at each vertex. Must have the same length
#'   as \code{adj.list}.
#' @param lmax.list List of local maxima information. See \code{\link{lmax.basins}}
#'   for format details.
#' @param lmin.list List of local minima information. See \code{\link{lmin.basins}}
#'   for format details.
#' @param verbose Logical; if \code{TRUE}, prints progress messages during
#'   computation. Default is \code{FALSE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{flow_directions}{List where each element describes the gradient flow at a vertex.}
#'   \item{lmax_basin_assignment}{Integer vector indicating which maximum basin
#'     each vertex belongs to, with \code{NA} for unassigned vertices}
#'   \item{lmin_basin_assignment}{Integer vector indicating which minimum basin
#'     each vertex belongs to, with \code{NA} for unassigned vertices}
#'   \item{ambiguous_vertices}{Integer vector of vertices that belong to multiple
#'     basins of the same type}
#'   \item{basin_stats}{List containing statistics about the basins:
#'       - n_max_basins: Number of maximum basins
#'       - n_min_basins: Number of minimum basins
#'       - coverage: Proportion of vertices assigned to at least one basin
#'       - max_basin_sizes: Named vector of maximum basin sizes
#'       - min_basin_sizes: Named vector of minimum basin sizes
#'   }
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Computes basins of attraction for all local maxima using
#'     \code{\link{lmax.basins}}
#'   \item Computes basins of attraction for all local minima using
#'     \code{\link{lmin.basins}}
#'   \item Assigns each vertex to appropriate basins, identifying conflicts
#'   \item Determines gradient flow directions based on neighboring values
#'   \item For unassigned vertices, attempts to infer basin membership from
#'     neighbors
#' }
#'
#' Gradient flow directions indicate the steepest ascent and descent paths
#' from each vertex. When edge weights are provided, they influence the
#' determination of these paths by scaling the effective gradients.
#'
#' @examples
#' \dontrun{
#' # Create example graph
#' adj.list <- list(c(2), c(1,3), c(2,4,5), c(3), c(3))
#' y <- c(0.2, 0.5, 0.8, 1.0, 0.1)
#'
#' # Define extrema
#' lmax.list <- list(list(lmax=4, vertices=c(4), label="peak"))
#' lmin.list <- list(list(lmin=5, vertices=c(5), label="valley"))
#'
#' # Compute gradient flow
#' flow_info <- compute.graph.gradient.flow(adj.list, NULL, y,
#'                                          lmax.list, lmin.list)
#'
#' # Examine results
#' print(flow_info$basin_stats)
#' }
#'
#' @references
#' Morse, M. (1934). The Calculus of Variations in the Large.
#' American Mathematical Society.
#'
#' Forman, R. (1998). Morse theory for cell complexes.
#' Advances in Mathematics, 134(1), 90-145.
#'
#' @seealso \code{\link{lmax.basins}}, \code{\link{lmin.basins}},
#'   \code{\link{set.boundary}}
#' @export
compute.graph.gradient.flow <- function(adj.list,
                                        weight.list = NULL,
                                        y,
                                        lmax.list,
                                        lmin.list,
                                        verbose = FALSE) {
    ## Input validation
    if (!is.list(adj.list) || length(adj.list) == 0) {
        stop("adj.list must be a non-empty list")
    }

    n_vertices <- length(adj.list)

    if (!is.numeric(y) || length(y) != n_vertices) {
        stop("y must be a numeric vector with length equal to adj.list")
    }

    ## Compute basins of attraction
    if (verbose) {
        message("Computing basins of attraction for local maxima...")
    }
    max_basins <- lmax.basins(adj.list, weight.list, y, lmax.list, verbose)

    if (verbose) {
        message("Computing basins of attraction for local minima...")
    }
    min_basins <- lmin.basins(adj.list, weight.list, y, lmin.list, verbose)

    ## Initialize result structures
    lmax_basin_assignment <- rep(NA_character_, n_vertices)
    lmin_basin_assignment <- rep(NA_character_, n_vertices)
    flow_directions <- vector("list", n_vertices)
    ambiguous_vertices <- integer(0)

    ## Assign vertices to maximum basins
    for (basin_id in names(max_basins)) {
        vertices <- max_basins[[basin_id]]
        for (v in vertices) {
            if (!is.na(lmax_basin_assignment[v])) {
                ambiguous_vertices <- c(ambiguous_vertices, v)
                if (verbose) {
                    message("Vertex ", v, " assigned to multiple max basins: ",
                            lmax_basin_assignment[v], " and ", basin_id)
                }
            }
            lmax_basin_assignment[v] <- basin_id
        }
    }

    ## Assign vertices to minimum basins
    for (basin_id in names(min_basins)) {
        vertices <- min_basins[[basin_id]]
        for (v in vertices) {
            if (!is.na(lmin_basin_assignment[v])) {
                ambiguous_vertices <- c(ambiguous_vertices, v)
                if (verbose) {
                    message("Vertex ", v, " assigned to multiple min basins: ",
                            lmin_basin_assignment[v], " and ", basin_id)
                }
            }
            lmin_basin_assignment[v] <- basin_id
        }
    }

    ## Remove duplicates from ambiguous vertices
    ambiguous_vertices <- unique(ambiguous_vertices)

    ## Compute gradient flow directions
    if (verbose) {
        message("Computing gradient flow directions...")
    }

    for (v in seq_len(n_vertices)) {
        ## Get neighbors, excluding self-loops
        neighbors <- setdiff(adj.list[[v]], v)

        if (length(neighbors) == 0) {
            flow_directions[[v]] <- list(type = "isolated")
            next
        }

        ## Get neighbor values
        neighbor_values <- y[neighbors]

        ## Check if vertex is an extremum
        if (all(y[v] >= neighbor_values)) {
            ## Local maximum
            flow_directions[[v]] <- list(
                type = "source",
                basin_id = lmax_basin_assignment[v]
            )
        } else if (all(y[v] <= neighbor_values)) {
            ## Local minimum
            flow_directions[[v]] <- list(
                type = "sink",
                basin_id = lmin_basin_assignment[v]
            )
        } else {
            ## Regular flow point
            ## Find steepest ascent and descent
            if (is.null(weight.list)) {
                ## Unweighted case
                ascent_idx <- which.max(neighbor_values - y[v])
                descent_idx <- which.max(y[v] - neighbor_values)
            } else {
                ## Weighted case
                edge_indices <- which(adj.list[[v]] %in% neighbors)
                weights <- weight.list[[v]][edge_indices]

                ## Compute effective gradients
                ascent_gradients <- (neighbor_values - y[v]) * weights
                descent_gradients <- (y[v] - neighbor_values) * weights

                ascent_idx <- which.max(ascent_gradients)
                descent_idx <- which.max(descent_gradients)
            }

            flow_directions[[v]] <- list(
                type = "regular",
                ascent_to = neighbors[ascent_idx],
                ascent_basin = lmax_basin_assignment[v],
                descent_to = neighbors[descent_idx],
                descent_basin = lmin_basin_assignment[v]
            )
        }
    }

    ## Handle uncovered vertices
    uncovered <- which(is.na(lmax_basin_assignment) & is.na(lmin_basin_assignment))

    if (length(uncovered) > 0) {
        if (verbose) {
            message("Found ", length(uncovered),
                    " vertices not covered by any basin. Attempting to infer...")
        }

        for (v in uncovered) {
            neighbors <- setdiff(adj.list[[v]], v)

            if (length(neighbors) > 0) {
                ## Get basin assignments of neighbors
                neighbor_max_basins <- lmax_basin_assignment[neighbors]
                neighbor_min_basins <- lmin_basin_assignment[neighbors]

                ## Remove NAs
                neighbor_max_basins <- neighbor_max_basins[!is.na(neighbor_max_basins)]
                neighbor_min_basins <- neighbor_min_basins[!is.na(neighbor_min_basins)]

                ## Use most common basin among neighbors
                if (length(neighbor_max_basins) > 0) {
                    basin_counts <- table(neighbor_max_basins)
                    lmax_basin_assignment[v] <- names(basin_counts)[which.max(basin_counts)]
                }

                if (length(neighbor_min_basins) > 0) {
                    basin_counts <- table(neighbor_min_basins)
                    lmin_basin_assignment[v] <- names(basin_counts)[which.max(basin_counts)]
                }

                ## Update flow direction
                if (!is.na(lmax_basin_assignment[v]) || !is.na(lmin_basin_assignment[v])) {
                    flow_directions[[v]]$type <- "inferred"
                    flow_directions[[v]]$ascent_basin <- lmax_basin_assignment[v]
                    flow_directions[[v]]$descent_basin <- lmin_basin_assignment[v]
                }
            }
        }
    }

    ## Compute basin statistics
    max_basin_sizes <- sapply(max_basins, length)
    min_basin_sizes <- sapply(min_basins, length)

    n_assigned <- sum(!is.na(lmax_basin_assignment) | !is.na(lmin_basin_assignment))
    coverage <- n_assigned / n_vertices

    basin_stats <- list(
        n_max_basins = length(max_basins),
        n_min_basins = length(min_basins),
        coverage = coverage,
        max_basin_sizes = max_basin_sizes,
        min_basin_sizes = min_basin_sizes
    )

    ## Return results
    return(list(
        flow_directions = flow_directions,
        lmax_basin_assignment = lmax_basin_assignment,
        lmin_basin_assignment = lmin_basin_assignment,
        ambiguous_vertices = sort(ambiguous_vertices),
        basin_stats = basin_stats
    ))
}
