##
##  Morse-Smale complexe utility functions
##

#' Constructs the Morse-Smale complex of a function defined on graph vertices
#'
#' @details This function computes gradient flow trajectories and Morse-Smale complex
#' components for a function Ey defined on vertices of a graph G = pIG_k^h(X).
#'
#' Key concepts and definitions:
#' 1. Gradient Flow Trajectories:
#'    - For each vertex v, computes ascending and descending trajectories
#'    - Ascending trajectory follows steepest ascent path to local maximum
#'    - Descending trajectory follows steepest descent path to local minimum
#'    - Trajectories are constrained to valid paths in core graph pIG_k(X)
#'
#' 2. Morse-Smale Pro-cells:
#'    - A pro-cell is defined by a pair (M,m) where M is a local maximum and m is a local minimum
#'    - Contains all vertices that lie on gradient trajectories starting at m and ending at M
#'    - Formally: pro-cell(M,m) = \{v \eqn{\in} V | there exists a trajectory through v from m to M\}
#'
#' 3. Morse-Smale Cells:
#'    - Connected components of pro-cells after removing their defining extrema
#'    - For pro-cell(M,m), compute components of subgraph induced by:
#'      vertices in pro-cell(M,m) \\ \{M,m\}
#'    - Unlike classical Morse theory, cells may not be disjoint
#'
#' Implementation Overview:
#' 1. Trajectory Computation:
#'    a. For each vertex v:
#'       - Compute ascending trajectory following steepest ascent
#'       - Compute descending trajectory following steepest descent
#'       - Store endpoints as local extrema
#'       - Update connectivity maps between extrema
#'
#' 2. Pro-cell Construction:
#'    - During trajectory computation, for each vertex v:
#'      - If trajectory connects local minimum m to maximum M
#'      - Add all trajectory vertices to pro-cell(M,m)
#'
#' 3. Cell Decomposition:
#'    - For each pro-cell(M,m):
#'      a. Remove M and m from vertex set
#'      b. Compute connected components of remaining subgraph
#'      c. Each component is a Morse-Smale cell
#'
#' @param graph A list representing the graph adjacency list. Each element of
#'              the list should be an integer vector specifying the neighboring vertex
#'              indices for a given vertex. The vertex indices should be 1-based.
#'
#' @param hop.list Optional list of hop constraints for gradient flow. If NULL, no
#'                 constraints are applied. If provided, should contain hop information
#'                 for each vertex.
#'
#' @param core.graph A list representing the core graph adjacency list used to constrain
#'                   gradient flow paths. Each element should be an integer vector of
#'                   1-based neighbor indices.
#'
#' @param Ey A numeric vector of function values at vertices.
#'
#' @return An R list of class "graphMScx" with the following components:
#'   \itemize{
#'     \item trajectories: A list where element i contains the gradient flow trajectory
#'           for vertex i (both descending and ascending paths)
#'     \item lmax_to_lmin: A list where element i contains indices of local minima
#'           connected to the i-th local maximum by gradient trajectories
#'     \item lmin_to_lmax: A list where element i contains indices of local maxima
#'           connected to the i-th local minimum by gradient trajectories
#'     \item local_maxima: Integer vector containing indices of all local maxima
#'     \item local_minima: Integer vector containing indices of all local minima
#'     \item procell_keys: List of integer pairs \code{(max_idx, min_idx)} identifying each pro-cell
#'     \item procells: List where element i contains all vertices in the i-th pro-cell
#'     \item cells: List where \code{cells[[i]]} contains the connected components of
#'           pro-cell i after removing its defining extrema
#'   }
#'
#' @details This function computes the Morse-Smale complex components:
#'   \itemize{
#'     \item Gradient flow trajectories following steepest ascent/descent paths
#'     \item Local extrema (maxima and minima) of the function
#'     \item Morse-Smale pro-cells: sets of vertices whose trajectories flow from
#'           a specific minimum to a specific maximum
#'     \item Morse-Smale cells: connected components of pro-cells with extrema removed
#'   }
#'
#'   The function converts indices between R's 1-based and C++'s 0-based systems.
#'   All returned indices are 1-based following R convention.
#'
#' @note Unlike classical Morse-Smale complexes on manifolds, the cells of the
#'       graph Morse-Smale complex may not be disjoint.
#'
graph.MS.cx <- function(graph, hop.list = NULL, core.graph, Ey) {
    if (!is.numeric(Ey)) {
        stop("Ey has to be a numeric vector")
    }

    if (length(graph) != length(Ey)) {
        cat("length(graph):", length(graph),"\n")
        cat("length(Ey):", length(Ey),"\n")
        stop("The length of Ey has to be the same as the number of components of the graph's adjacency matrix.")
    }

    if (length(core.graph) != length(Ey)) {
        cat("length(core.graph):", length(core.graph), "\n")
        cat("length(Ey):", length(Ey), "\n")
        stop("The length of Ey has to be the same as the number of components of the core.graph's adjacency matrix.")
    }

    ## Convert adjacency lists to 0-based indexing
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))
    core.graph.0based <- lapply(core.graph, function(x) as.integer(x - 1))

    ## Call C++ function
    results <- NULL
    if (!is.null(hop.list)) {
        hop.list <- lapply(hop.list, function(x) as.integer(x))
        results <- .Call("S_graph_MS_cx_using_short_h_hops",
                         graph.0based,
                         hop.list,
                         core.graph.0based,
                         Ey)
    } else {
        results <- .Call("S_graph_MS_cx_with_path_search",
                         graph.0based,
                         core.graph.0based,
                         Ey)
    }

    ## Convert all indices back to 1-based
    results$trajectories <- lapply(results$trajectories, function(x) as.integer(x + 1))
    results$local_maxima <- as.integer(results$local_maxima + 1)
    results$local_minima <- as.integer(results$local_minima + 1)

    ## Convert indices in connectivity maps
    results$lmax_to_lmin <- lapply(results$lmax_to_lmin, function(x) as.integer(x + 1))
    results$lmin_to_lmax <- lapply(results$lmin_to_lmax, function(x) as.integer(x + 1))

    ## Convert procell keys (pairs of indices)
    results$procell_keys <- lapply(results$procell_keys, function(x) as.integer(x + 1))

    ## Convert vertices in procells
    results$procells <- lapply(results$procells, function(x) as.integer(x + 1))

    ## Convert vertices in cells (list of lists of components)
    results$cells <- lapply(results$cells, function(x) {
        lapply(x, function(comp) as.integer(comp + 1))
    })

    class(results) <- "graphMScx"
    results
}


## #' Estimates graph gradient flow trajectories of a response variable for each vertex of a graph.
## #'
## #' This function computes the gradient flow trajectories of a response variable over a graph,
## #' identifying the ascending and descending paths from each vertex to its local minimum and maximum.
## #'
## #' @param graph A list representing the graph adjacency list. Each element of the list should be
## #'              an integer vector specifying the neighboring vertex indices for a given vertex.
## #'              The vertex indices should be 1-based.
## #' @param Ey A numeric vector containing estimates of the conditional expectation of a response
## #'           variable over the graph vertices.
## #'
## #' @return An R list with the following components:
## #'         - "lext": An integer matrix with dimensions (n_vertices, 2), where the first column
## #'                   contains the local minimum indices and the second column contains the local
## #'                   maximum indices for each vertex.
## #'         - "trajectories": A list of integer vectors representing the gradient flow trajectories
## #'                           for each vertex of the graph. Each trajectory includes both the
## #'                           descending and ascending paths from the vertex to its local minimum
## #'                           and maximum, respectively.
## #'
## #' @examples
## #' # Create a sample graph adjacency list
## #' graph <- list(c(2, 3), c(1, 3), c(1, 2))
## #'
## #' # Create a sample response variable estimate
## #' Ey <- c(0.5, 0.2, 0.8)
## #'
## #' # Compute the gradient flow trajectories
## #' result <- graph.gradient.flow.trajectories(graph, Ey)
## #'
## #' # Access the local extrema matrix
## #' lext <- result$lext
## #'
## #' # Access the gradient flow trajectories
## #' trajectories <- result$trajectories
## #'
## #' @export
## graph.gradient.flow.trajectories <- function(graph, Ey) {
##     if (!is.numeric(Ey)) {
##         stop("Ey has to be a numeric vector")
##     }

##     if (length(graph) != length(Ey)) {
##         cat("length(graph):", length(graph), "\n")
##         cat("length(Ey):", length(Ey), "\n")
##         stop("The length of Ey has to be the same as the number of components of the graph's adjacency matrix.")
##     }

##     ## Converting each component of graph adjacency list to an integer vector
##     graph.0based <- lapply(graph, function(x) as.integer(x - 1))

##     res <- .Call("S_graph_gradient_flow_trajectories", graph.0based, Ey)

##     res$trajectories <- lapply(res$trajectories, function(x) as.integer(x + 1))

##     res
## }


#' Compute Constrained Gradient Flow Trajectories
#'
#' @description
#' Computes gradient flow trajectories on a graph while preventing basin-jumping between
#' different regions of attraction by enforcing monotonicity along paths in a core graph.
#'
#' @param graph List of integer vectors representing the h-th power graph pIG_k^h(X)
#' @param core.graph List of integer vectors representing the core graph pIG_k(X)
#' @param Ey Numeric vector of function values at vertices
#'
#' @details
#' The function ensures that trajectories follow paths that are monotonic in the core graph,
#' preventing undesirable jumps between basins of attraction of different local maxima.
#' Indices in graph and core.graph are converted from 1-based (R) to 0-based (C++) indexing.
#'
#' @return List containing:
#' \itemize{
#'   \item lext - Matrix with two columns containing local minima and maxima for each vertex
#'   \item trajectories - List of integer vectors containing complete trajectories
#' }
#'
#' @examples
#' \dontrun{
#' trajectories <- graph.constrained.gradient.flow.trajectories(graph, core.graph, Ey)
#' local.extrema <- trajectories$lext
#' paths <- trajectories$trajectories
#' }
#'
#' @export
graph.constrained.gradient.flow.trajectories <- function(graph, core.graph, Ey) {
    if (!is.numeric(Ey)) {
        stop("Ey has to be a numeric vector")
    }

    if (length(graph) != length(Ey)) {
        cat("length(graph):", length(graph), "\n")
        cat("length(Ey):", length(Ey), "\n")
        stop("The length of Ey has to be the same as the number of components of the graph's adjacency matrix.")
    }

    if (length(core.graph) != length(Ey)) {
        cat("length(core.graph):", length(core.graph), "\n")
        cat("length(Ey):", length(Ey), "\n")
        stop("The length of Ey has to be the same as the number of components of the core.graph's adjacency matrix.")
    }

    ## Converting each component of graph and core.graph adjacency lists to an integer vector
    graph.0based <- lapply(graph, function(x) as.integer(x - 1))
    core.graph.0based <- lapply(core.graph, function(x) as.integer(x - 1))

    res <- .Call("S_graph_constrained_gradient_flow_trajectories",
                 graph.0based,
                 core.graph.0based,
                 Ey)

    ## res$trajectories <- lapply(res$trajectories, function(x) as.integer(x + 1))

    res
}



#' Finds the Longest Stretch of Consecutive TRUE Values in a Logical Vector
#'
#' This function takes a logical vector as input and returns the starting index
#' and length of the longest sub-interval where the values are TRUE.
#'
#' @param l A logical vector.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{start_index}{The index at which the longest stretch of consecutive TRUE values starts.}
#'   \item{length}{The length of the longest stretch of consecutive TRUE values.}
#' }
#'
#' @note If there are multiple disjoined stretch of TRUE values of the same length, then the first one is returned.
#'
#' @examples
#' l <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' result <- find.longest.true.stretch(l)
#' print(result)
#' # $start_index
#' # [1] 4
#' # $length
#' # [1] 3
#'
#' ## Two longest TRUE stretches
#' l <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' find.longest.true.stretch(l)
#' # $start.index
#' # [1] 1
#' # $length
#' # [1] 3
#'
#' ## No TRUE values
#' l <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
#' find.longest.true.stretch(l)
#' # $start.index
#' # [1] NA
#' # $length
#' # [1] 0
#'
#' @export
find.longest.true.stretch <- function(l) {
    ## Ensure the input is a logical vector
    if (!is.logical(l)) {
        stop("Input must be a logical vector.")
    }

    ## Initialize variables
    max.length <- 0
    max.start.index <- NA
    current.length <- 0
    current.start.index <- NA

    for (i in seq_along(l)) {
        if (l[i]) {
            ## If current element is TRUE, increment current length
            if (current.length == 0) {
                current.start.index <- i
            }
            current.length <- current.length + 1
        } else {
            ## If current element is FALSE, check if current length is the longest
            if (current.length > max.length) {
                max.length <- current.length
                max.start.index <- current.start.index
            }
            ## Reset current length
            current.length <- 0
            current.start.index <- NA
        }
    }

    ## Final check at the end of the vector
    if (current.length > max.length) {
        max.length <- current.length
        max.start.index <- current.start.index
    }

    ## Return the result as a list
    list(start.index = max.start.index, length = max.length)
}

#' Morse-Smale Graph Persistence Analysis
#'
#' This function analyzes the persistence of Morse-Smale graphs over a range of nearest neighbor values (k) for a given data matrix X and response variable y.
#' It computes the nerve graphs and their associated Morse-Smale complex graphs for each k value and determines the isomorphism and deviation from isomorphism
#' between consecutive Morse-Smale graphs.
#'
#' @param X            A data matrix or data frame.
#' @param y            A numeric vector representing the response variable over X.
#' @param Ks           A vector of positive integer values specifying the range of k values.
#' @param ref.MS.graph A reference Morse-Smale graph. Usually, the Morse-Smale
#'     graph of the y over the full graph or of full dataset before subsampling.
#'     If ref.MS.graph is not NULL, is.isomorphic.with.ref.MS.graph and
#'     dev.from.isomorphism.with.MS.graph are computed.
#' @param allow.disconnected.graph Logical. If TRUE (default), the algorithm proceeds even if none of the nerve graphs associated with (X, y) are connected.
#' @param verbose If set to TRUE, it prints progress messages.
#'
#' @return A list with four components:
#'   - nerve.graph: A list of nerve graphs computed for each k value.
#'   - MS.graph: A list of Morse-Smale complex graphs computed for each k value.
#'   - are.consecutive.graphs.isomorphic: A vector indicating whether consecutive Morse-Smale graphs are isomorphic.
#'   - is.isomorphic.with.ref.MS.graph: A vector indicating whether the computed Morse-Smale graphs are isomorphic with the reference Morse-Smale graph.
#'   - first.index.of.longest.stretch.of.isomorphic.MS.graphs: The first index of the longest stretch of isomorphic MS graphs.
#'   - length.of.longest.stretch.of.isomorphic.MS.graphs: The length of the longest stretch of isomorphic MS graphs.
#'
morse.smale.graph.persistence <- function(X,
                                          y,
                                          Ks,
                                          ref.MS.graph = NULL,
                                          allow.disconnected.graph = TRUE,
                                          verbose = FALSE) {

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

    if (!is.numeric(y)) {
        stop("y has to be a numeric vector")
    }

    if (!is.vector(y)) {
        stop("y has to be a vector")
    }

    if (length(y) != nrow(X)) {
        stop("The number of rows of X has to be the same as the length of y")
    }

    if (!is.vector(Ks)) {
        stop("Ks has to be a vector")
    }

    if (!all(Ks == floor(Ks))) {
        stop("Ks have to be a vector of integers")
    }

    if (min(Ks) < 2) {
        stop("The smallest value of Ks has to be at least 2")
    }

    if ( verbose ) {
        routine.ptm <- proc.time()
        ptm <- proc.time()
        cat("Creating nerve graphs ... ")
    }

    ## Computing nerve graphs of X for all k within Ks
    nerve.graph <- list()
    n.cc <- numeric(max(Ks))

    for (k in Ks) {
        if ( verbose ) {
            cat("\rCreating nerve graphs:",k)
        }
        g <- IW.kNN.graph(X, k)
        nerve.graph[[k]] <- g$adj_list
        n.cc[k] <- length(table(g$conn_comps))
    }
    if ( verbose ) {
        cat("\rCreating nerve graphs ")
        elapsed.time(ptm)
    }

    ## Looking for the smallest k within Ks for which the nerve graph of X given
    ## k is connected (if such k exists within Ks).
    k0 <- min(which(n.cc == 1))

    if (is.infinite(k0)) {
        cat("WARNING: None of the nerve graphs computed for k's in Ks is connected!\n")
        if (!allow.disconnected.graph) {
            ret <- list(
                nerve.graph = nerve.graph,
                MS.graph = NULL,
                are.consecutive.graphs.isomorphic = NULL,
                dev.from.isomorphism = NULL
            )
            return(ret)
        }
        k0 <- min(Ks)
    }

    if ( verbose ) {
        cat("Computing the Morse-Smale complex graph ... ")
    }

    ## Computing the Morse-Smale complex graph for every nerve graph for k's in k0:max(Ks) range
    ## Computing the value of
    ## are.consecutive.graphs.isomorphic[k]
    ## dev.from.isomorphism[k]
    ## for k > k0
    MS.graph <- list()
    MS.igraph <- list()
    are.consecutive.graphs.isomorphic <- logical(max(Ks))
    is.isomorphic.with.ref.MS.graph <- logical(max(Ks))
    ## dev.from.isomorphism <- numeric(max(Ks))
    ## dev.from.isomorphism.with.MS.graph <- numeric(max(Ks))

    for (k in k0:max(Ks)) {
        if ( verbose ) {
            cat("\rComputing the Morse-Smale complex graph:",k)
        }

        MS.graph[[k]] <- graph.MS.cx(nerve.graph[[k]], y)$MS_graph

        ## Creating an igraph from the adjacency list MS.graph[[k]]
        g.m <- convert.adjacency.to.edge.matrix(MS.graph[[k]])
        MS.igraph[[k]] <- igraph::graph_from_edgelist(g.m, directed = FALSE)

        ## Check if MS.igraph[[k]] is an igraph object
        if (!inherits(MS.igraph[[k]], "igraph")) {
            stop(paste("MS.igraph at k =", k, "is not an igraph object."))
        }

        if (k > k0) {
            are.consecutive.graphs.isomorphic[k] <- igraph::graph.isomorphic(MS.igraph[[k - 1]], MS.igraph[[k]])
            ##dev.from.isomorphism[k] <- hungarian.frobenius.graph.matching(MS.graph[[k - 1]], MS.graph[[k]])
        }
    }
    if ( verbose ) {
        cat("\rComputing the Morse-Smale complex graph ")
        elapsed.time(ptm)
    }

    if (!is.null(ref.MS.graph)) {

        if ( verbose ) {
            cat("\rTesting isomorphism with the reference graph ... ")
        }

        g.m <- convert.adjacency.to.edge.matrix(ref.MS.graph)
        ref.MS.igraph <- igraph::graph_from_edgelist(g.m, directed = FALSE)

        for (k in k0:max(Ks)) {
            if ( verbose ) {
                cat("\rk:",k)
            }

            is.isomorphic.with.ref.MS.graph[k] <- igraph::graph.isomorphic(MS.igraph[[k]], ref.MS.igraph)
            ##dev.from.isomorphism.with.MS.graph[k] <- hungarian.frobenius.graph.matching(MS.graph[[k]], ref.MS.graph)
        }
        if ( verbose ) {
            cat("\rTesting isomorphism with the reference graph ... ")
            elapsed.time(ptm)
        }
    }

    if ( verbose ) {
        txt <- sprintf("Total elapsed time")
        elapsed.time(routine.ptm, txt, with.brackets = FALSE)
    }

    longest.true.stretch <- find.longest.true.stretch(are.consecutive.graphs.isomorphic)

    list(
        nerve.graph = nerve.graph,
        MS.graph = MS.graph,
        are.consecutive.graphs.isomorphic = are.consecutive.graphs.isomorphic,
        first.index.of.longest.stretch.of.isomorphic.MS.graphs = longest.true.stretch$start.index - 1, ## subtracting 1 as TRUE at the start.index position of are.consecutive.graphs.isomorphic means that the start.index - 1 graph is isomorphic to the start.index graph.
        length.of.longest.stretch.of.isomorphic.MS.graphs = longest.true.stretch$length + 1,
        is.isomorphic.with.ref.MS.graph = is.isomorphic.with.ref.MS.graph
    )
}


#' Create Labels for Local Maxima and Minima in Morse-Smale Complex
#'
#' @description
#' Creates unique labels for local maxima and minima in a Morse-Smale complex by
#' concatenating two-letter shortcuts of the most abundant taxa. The function ensures
#' uniqueness of labels by iteratively including additional taxa if needed.
#'
#' @param MS.res A list containing Morse-Smale complex results with component MS_cx
#' @param state.space Matrix of ASV counts/abundances with samples as columns and ASVs as rows
#' @param taxonomy Taxonomy information for ASVs
#' @param freq.thld Threshold for frequency filtering (default: 100)
#' @param min.relAb.thld Minimum relative abundance threshold for including taxa in labels (default: 0.05)
#'
#' @return A list with four components:
#'   \item{lmax.labels}{Named character vector of local maxima labels}
#'   \item{lmin.labels}{Named character vector of local minima labels}
#'   \item{lmax.profiles}{Named list of relative abundance profiles for local maxima, where names are point indices and each profile is a matrix with taxa in rows and two columns: species names and their relative abundances}
#'   \item{lmin.profiles}{Named list of relative abundance profiles for local minima, where names are point indices and each profile is a matrix with taxa in rows and two columns: species names and their relative abundances}
#'
#' @details
#' The function processes both local maxima and minima, creating unique labels by:
#' 1. Filtering points by frequency threshold
#' 2. For each point, identifying taxa above the relative abundance threshold
#' 3. Creating initial labels using two-letter shortcuts
#' 4. Ensuring uniqueness by incorporating additional taxa if needed
#'
#' @examples
#' \dontrun{
#' labels <- create.lmin.lmax.labels(MS.res, state.space, taxonomy,
#'                                   freq.thld = 100, min.relAb.thld = 0.05)
#' print(labels$lmax.labels)
#' print(labels$lmin.labels)
#' }
#' @export
create.lmin.lmax.labels <- function(MS.res, state.space, taxonomy, freq.thld = 100, min.relAb.thld = 0.05) {
    # Helper function to create two-letter shortcuts from species names
    create.shortcut <- function(sp.name) {
        parts <- strsplit(sp.name, "_")[[1]]
        if (length(parts) >= 2) {
            paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))
        } else {
            paste0(substr(parts[1], 1, 1), "")
        }
    }

    # Helper function to create labels and store profiles
    create.labels <- function(points, type = "max") {
        # Get frequencies and filter
        freq <- sort(table(MS.res$MS_cx[, paste0("local_", type)]), decreasing = TRUE)
        freq <- freq[freq > freq.thld]

        # Initialize labels vector and profiles list
        labels <- character(length(freq))
        names(labels) <- names(freq)
        profiles <- list()

        # First pass: create initial labels and store profiles
        for (i in seq_along(freq)) {
            point.i <- as.integer(names(freq)[i])
            id <- rownames(state.space)[point.i]

            # Get profile for this point
            profile <- prof.fn(id, state.space, taxonomy)
            # Store profile
            profiles[[as.character(point.i)]] <- profile

            abundances <- as.numeric(profile[, 2])
            above.threshold <- which(abundances >= min.relAb.thld)

            if (length(above.threshold) > 0) {
                taxa.names <- profile[above.threshold, 1]
                shortcuts <- sapply(taxa.names, create.shortcut)
                labels[i] <- paste(shortcuts, collapse = "")
            } else {
                labels[i] <- create.shortcut(profile[1, 1])
            }
        }

        # Second pass: ensure uniqueness
        while (any(duplicated(labels))) {
            dups <- labels[duplicated(labels)]
            for (dup.label in unique(dups)) {
                dup.indices <- which(labels == dup.label)
                for (idx in dup.indices) {
                    point.i <- as.integer(names(freq)[idx])
                    profile <- profiles[[as.character(point.i)]]

                    current.taxa <- ceiling(nchar(labels[idx])/2)
                    if (current.taxa < nrow(profile)) {
                        new.taxon <- create.shortcut(profile[current.taxa + 1, 1])
                        labels[idx] <- paste0(labels[idx], new.taxon)
                    }
                }
            }
        }

        return(list(
            labels = labels,
            profiles = profiles
        ))
    }

    # Create labels and store profiles for both maxima and minima
    lmax.results <- create.labels(MS.res$MS_cx[, "local_max"], "max")
    lmin.results <- create.labels(MS.res$MS_cx[, "local_min"], "min")

    return(list(
        lmax.labels = lmax.results$labels,
        lmin.labels = lmin.results$labels,
        lmax.profiles = lmax.results$profiles,
        lmin.profiles = lmin.results$profiles
    ))
}


#' Create Label Tables and Indicator Vectors for Morse-Smale Complex Critical Points
#'
#' @description
#' Creates label tables and indicator vectors for local minima and maxima in a Morse-Smale
#' complex. The label tables map point IDs to their corresponding labels, while indicator
#' vectors mark the presence/absence of critical points in the state space.
#'
#' @param lmin.labels Named character vector of labels for local minima
#' @param lmax.labels Named character vector of labels for local maxima
#' @param state.space State space matrix where rownames correspond to point IDs
#'
#' @return A list with four components:
#'   \item{lmin.lab.tbl}{Named character vector mapping point IDs to local minima labels}
#'   \item{lmax.lab.tbl}{Named character vector mapping point IDs to local maxima labels}
#'   \item{lmin.ind}{Numeric vector indicating local minima (1) and non-minima (0)}
#'   \item{lmax.ind}{Numeric vector indicating local maxima (1) and non-maxima (0)}
#'
#' @details
#' The function processes both local minima and maxima labels to create:
#' 1. Label tables that map point IDs to their corresponding labels
#' 2. Binary indicator vectors marking the presence (1) or absence (0) of critical points
#' All vectors maintain the same length as the number of rows in state.space and use consistent naming.
#'
#' @examples
#' \dontrun{
#' # Given lmin.labels, lmax.labels and state space matrix state.space
#' result <- create.lmin.lmax.label.indicators(lmin.labels, lmax.labels, state.space)
#' print(head(result$lmin.lab.tbl))
#' print(sum(result$lmin.ind)) # Number of local minima
#' }
#' @export
create.lmin.lmax.label.indicators <- function(lmin.labels, lmax.labels, state.space) {
    # Helper function to process one type of critical points
    create.single.indicators <- function(labels) {
        # Convert label names to integers
        point.indices <- as.integer(names(labels))

        # Get corresponding row names from state.space
        point.ids <- rownames(state.space)[point.indices]

        # Create label table
        lab.tbl <- c()
        lab.tbl[point.ids] <- labels[as.character(point.indices)]

        # Create indicator vector
        ind <- numeric(nrow(state.space))
        names(ind) <- rownames(state.space)
        ind[point.ids] <- 1

        return(list(lab.tbl = lab.tbl, ind = ind))
    }

    # Process local minima
    lmin.result <- create.single.indicators(lmin.labels)

    # Process local maxima
    lmax.result <- create.single.indicators(lmax.labels)

    # Return combined results
    return(list(
        lmin.lab.tbl = lmin.result$lab.tbl,
        lmax.lab.tbl = lmax.result$lab.tbl,
        lmin.ind = lmin.result$ind,
        lmax.ind = lmax.result$ind
    ))
}

#' Create Frequency Tables and LaTeX Captions for Morse-Smale Complex Critical Points
#'
#' @description
#' Generates detailed frequency tables and corresponding LaTeX captions for local maxima
#' and minima in a Morse-Smale complex. The tables include basin of attraction sizes,
#' relative conditional expectation statistics, and relationships between critical points.
#'
#' @param MS.res List containing Morse-Smale complex results with component MS_cx
#' @param lmin.labels Named character vector of labels for local minima
#' @param lmax.labels Named character vector of labels for local maxima
#' @param rel.condEy Numeric vector of relative conditional expectations
#' @param freq.thld Frequency threshold for filtering (default: 100)
#' @param outcome.name Character string describing the outcome variable (default: "outcome")
#'
#' @return A list with four components:
#'   \item{lmax.freq.df}{Data frame containing local maxima statistics}
#'   \item{lmin.freq.df}{Data frame containing local minima statistics}
#'   \item{lmax.caption}{LaTeX caption for local maxima table}
#'   \item{lmin.caption}{LaTeX caption for local minima table}
#'
#' @details
#' For each critical point type, the function computes:
#' 1. Basin of attraction sizes
#' 2. Maximum and minimum relative conditional expectations within each basin
#' 3. Range of relative conditional expectations (Delta)
#' 4. Number of complementary critical points (minima for maxima and vice versa)
#' @export
create.lmin.lmax.frequency.tables <- function(MS.res,
                                              lmin.labels,
                                              lmax.labels,
                                              rel.condEy,
                                              freq.thld = 100,
                                              outcome.name = "outcome") {

    # Helper function to create frequency table for one type of critical point
    create.single.freq.table <- function(point.type, labels, other.type) {
        # Get frequencies
        freq <- sort(table(MS.res$MS_cx[, paste0("local_", point.type)]),
                    decreasing = TRUE)
        freq <- freq[freq > freq.thld]

        # Initialize frequency dataframe
        freq.df <- data.frame(BoA = as.integer(freq))
        rownames(freq.df) <- labels[names(freq)]

        # Get point map
        point.map <- MS.res$MS_cx[, paste0("local_", point.type)]

        # Calculate max and min relative conditional expectations
        max.rel <- min.rel <- numeric(length(freq))
        names(max.rel) <- names(min.rel) <- names(freq)

        for (i in seq_along(freq)) {
            point.id <- names(freq)[i]
            BoA <- which(as.integer(point.id) == point.map)
            max.rel[point.id] <- max(rel.condEy[BoA])
            min.rel[point.id] <- min(rel.condEy[BoA])
        }

        # Calculate number of complementary critical points
        n.other <- numeric(length(freq))
        names(n.other) <- names(freq)

        for (i in seq_along(freq)) {
            point.id <- names(freq)[i]
            idx <- MS.res$MS_cx[, paste0("local_", point.type)] == as.integer(point.id)
            f <- table(MS.res$MS_cx[idx, paste0("local_", other.type)])
            n.other[point.id] <- length(f)
        }

        # Combine all statistics
        freq.df <- cbind(
            freq.df,
            format(round(max.rel, 2), nsmall = 2),
            format(round(min.rel, 2), nsmall = 2),
            format(round(max.rel - min.rel, 2), nsmall = 2),
            n.other
        )

        colnames(freq.df) <- c("|BoA|",
                              "max(rel cond)",
                              "min(rel cond)",
                              "Delta(rel cond)",
                              paste0("n(", other.type, ")"))

        return(freq.df)
    }

    # Create frequency tables
    lmax.freq.df <- create.single.freq.table("max", lmax.labels, "min")
    lmin.freq.df <- create.single.freq.table("min", lmin.labels, "max")

    ## Creating captions
    # Create captions with no line breaks or extra spaces
    lmax.caption <- sprintf(
        "Key metrics for each local maximum of the %s conditional expectation: size of the basin of attraction \\(|\\text{BoA}|\\) (first column), maximum relative conditional expectation within the basin (second column), minimum relative conditional expectation within the basin (third column), the difference (\\(\\Delta\\)) between the maximum and minimum relative conditional expectations (fourth column), and the number of local minima within the basin (fifth column).",
        outcome.name
    )

    lmin.caption <- sprintf(
        "Key metrics for each local minimum of the %s conditional expectation: size of the basin of attraction \\(|\\text{BoA}|\\) (first column), maximum relative conditional expectation within the basin (second column), minimum relative conditional expectation within the basin (third column), the difference (\\(\\Delta\\)) between the maximum and minimum relative conditional expectations (fourth column), and the number of local maxima within the basin (fifth column).",
        outcome.name
    )

    return(list(
        lmax.freq.df = lmax.freq.df,
        lmin.freq.df = lmin.freq.df,
        lmax.caption = lmax.caption,
        lmin.caption = lmin.caption
    ))
}

#' Create Contingency Table of Local Maxima and Minima Relationships
#'
#' @description
#' Generates a contingency table showing the relationships between local maxima and
#' minima in a Morse-Smale complex, ordered by marginal frequencies. Also provides
#' a LaTeX caption describing the table contents.
#'
#' @param MS.res List containing Morse-Smale complex results with component MS_cx
#' @param lmin.labels Named character vector of labels for local minima
#' @param lmax.labels Named character vector of labels for local maxima
#'
#' @return A list with two components:
#'   \item{contingency_table}{Table showing frequencies of local maxima-minima combinations}
#'   \item{caption}{LaTeX caption describing the contingency table}
#'
#' @details
#' The function:
#' 1. Maps critical points to their labels
#' 2. Creates a contingency table of label combinations
#' 3. Orders rows and columns by marginal frequencies
#' 4. Generates an appropriate LaTeX caption
#'
#' @examples
#' \dontrun{
#' result <- create.lmin.lmax.contingency.table(MS.res, lmin.labels, lmax.labels)
#' print(result$contingency_table)
#' cat(result$caption)
#' }
#' @export
create.lmin.lmax.contingency.table <- function(MS.res, lmin.labels, lmax.labels) {
    # Get point mappings
    lmax.map <- MS.res$MS_cx[, "local_max"]
    lmin.map <- MS.res$MS_cx[, "local_min"]

    # Map points to labels
    lmax.labs <- lmax.labels[as.character(lmax.map)]
    lmin.labs <- lmin.labels[as.character(lmin.map)]

    # Create frequency tables
    lmax.labs.freq <- sort(table(lmax.labs), decreasing = TRUE)
    lmin.labs.freq <- sort(table(lmin.labs), decreasing = TRUE)

    # Create contingency table
    lmaxlmin.freq <- table(lmax.labs, lmin.labs)

    # Order by marginal frequencies
    row.order <- order(rowSums(lmaxlmin.freq), decreasing = TRUE)
    col.order <- order(colSums(lmaxlmin.freq), decreasing = TRUE)

    # Reorder table
    contingency.table <- lmaxlmin.freq[row.order, col.order]

    # Create caption
    caption <- paste(
        "Contingency table showing the relationships between local maxima (rows) and local minima (columns) in the Morse-Smale complex.",
        "Each entry represents the number of points in the state space that flow from the corresponding local minimum to the corresponding local maximum.",
        "Rows and columns are ordered by their marginal frequencies (decreasing order)."
    )

    return(list(
        contingency.table = contingency.table,
        caption = caption
    ))
}


#' Process Complete Morse-Smale Complex Analysis Pipeline
#'
#' @description
#' Executes a complete analysis pipeline for Morse-Smale complex data by sequentially
#' processing critical points, creating labels, generating indicators, computing
#' frequency tables, and forming contingency relationships.
#'
#' @param MS.res List containing Morse-Smale complex results with component MS_cx
#' @param state.space Matrix representing the state space (e.g., ASV abundance matrix)
#' @param taxonomy Taxonomy information for state space features
#' @param rel.condE Numeric vector of relative conditional expectations
#' @param freq.thld Frequency threshold for filtering (default: 100)
#' @param min.relAb.thld Minimum relative abundance threshold for labels (default: 0.05)
#' @param outcome.name Character string describing the outcome variable (default: "outcome")
#'
#' @return A list with four components:
#'   \item{labels}{Results from create.lmin.lmax.labels including lmin.labels and lmax.labels}
#'   \item{indicators}{Results from create.lmin.lmax.label.indicators including indicator vectors and tables}
#'   \item{frequencies}{Results from create.lmin.lmax.frequency.tables including frequency tables and captions}
#'   \item{contingency}{Results from create.lmin.lmax.contingency.table including contingency table and caption}
#'
#' @details
#' The function executes the following steps:
#' 1. Creates labels for local minima and maxima based on taxonomic composition
#' 2. Generates indicator vectors and label tables for critical points
#' 3. Computes frequency tables with statistics for critical points
#' 4. Creates contingency table showing relationships between local minima and maxima
#'
#' @examples
#' \dontrun{
#' results <- morse.smale.complex.lmin.lmax.summary(
#'     MS.res = smoothed.condE.MS.res,
#'     state.space = asv.matrix,
#'     taxonomy = asv.taxonomy,
#'     rel.condE = rel.smoothed.condE,
#'     outcome.name = "smoothed conditional expectation"
#' )
#' }
#' @export
morse.smale.complex.lmin.lmax.summary <- function(MS.res,
                                                  state.space,
                                                  taxonomy,
                                                  rel.condE,
                                                  freq.thld = 100,
                                                  min.relAb.thld = 0.05,
                                                  outcome.name = "outcome") {

    # Step 1: Create labels for critical points
    label.results <- create.lmin.lmax.labels(
        MS.res = MS.res,
        state.space = state.space,
        taxonomy = taxonomy,
        freq.thld = freq.thld,
        min.relAb.thld = min.relAb.thld
    )

    # Step 2: Generate indicator vectors and label tables
    indicator.results <- create.lmin.lmax.label.indicators(
        lmin.labels = label.results$lmin.labels,
        lmax.labels = label.results$lmax.labels,
        state.space = state.space
    )

    # Step 3: Compute frequency tables and statistics
    frequency.results <- create.lmin.lmax.frequency.tables(
        MS.res = MS.res,
        lmin.labels = label.results$lmin.labels,
        lmax.labels = label.results$lmax.labels,
        rel.condEy = rel.condE,
        freq.thld = freq.thld,
        outcome.name = outcome.name
    )

    # Step 4: Create contingency table
    contingency.results <- create.lmin.lmax.contingency.table(
        MS.res = MS.res,
        lmin.labels = label.results$lmin.labels,
        lmax.labels = label.results$lmax.labels
    )

    # Return all results in a structured list
    return(list(
        labels = label.results,
        indicators = indicator.results,
        frequencies = frequency.results,
        contingency = contingency.results
    ))
}
