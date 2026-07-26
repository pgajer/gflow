##
## Graph utilities
##

#' Create Pruned Intersection Size List
#'
#' This function generates a list of intersection sizes for the edges of a pruned adjacency list.
#' It takes the original adjacency list, intersection size list, and pruned adjacency list as inputs,
#' and returns a new list containing intersection sizes corresponding to the pruned graph structure.
#'
#' @param adj.list A list where each element represents a node and contains indices of its neighbors
#'   in the original graph.
#' @param isize.list A list of the same length as adj.list, where each element contains intersection
#'   sizes corresponding to the edges in adj.list.
#' @param pruned.adj.list A list representing the pruned graph, where each element contains indices
#'   of neighbors after removing redundant edges.
#'
#' @return A list of the same length as pruned.adj.list, where each element contains intersection
#'   sizes corresponding to the edges in the pruned graph.
#'
#' @details
#' The function iterates through each node in the pruned adjacency list. For each node, it identifies
#' the neighbors in the pruned graph and retrieves their corresponding intersection sizes from the
#' original isize.list. The resulting pruned.isize.list maintains the structure of pruned.adj.list
#' but contains intersection sizes instead of neighbor indices.
#'
#' @note
#' - The function assumes that all nodes present in pruned.adj.list are also present in adj.list
#'   and isize.list.
#' - If a neighbor in pruned.adj.list is not found in the original adj.list (which should not happen
#'   under normal circumstances), the corresponding intersection size will be set to NA.
#'
#' @examples
#' adj.list <- list(c(2,3,4), c(1,3), c(1,2,4), c(1,3))
#' isize.list <- list(c(2,1,3), c(2,1), c(1,1,2), c(3,2))
#' pruned.adj.list <- list(c(2,4), c(1), c(4), c(1,3))
#' pruned.isize.list <- create.pruned.isize.list(adj.list, isize.list, pruned.adj.list)
#' print(pruned.isize.list)
#'
#' @export
create.pruned.isize.list <- function(adj.list, isize.list, pruned.adj.list) {
    ## Initialize an empty list to store the pruned intersection sizes
    pruned.isize.list <- vector("list", length(pruned.adj.list))

    ## Iterate through each node in the pruned adjacency list
    for (i in seq_along(pruned.adj.list)) {
        ## Get the neighbors of the current node in the pruned graph
        pruned.neighbors <- pruned.adj.list[[i]]

        ## Initialize a vector to store intersection sizes for the current node
        pruned.isizes <- numeric(length(pruned.neighbors))

        ## Iterate through each neighbor in the pruned graph
        for (j in seq_along(pruned.neighbors)) {
            ## Get the index of the neighbor in the original adjacency list
            neighbor.index <- which(adj.list[[i]] == pruned.neighbors[j])

            ## If the neighbor is found, get its intersection size from isize.list
            if (length(neighbor.index) > 0) {
                pruned.isizes[j] <- isize.list[[i]][neighbor.index]
            } else {
                ## If not found (shouldn't happen), set to NA or 0
                pruned.isizes[j] <- NA
            }
        }

        ## Store the pruned intersection sizes for the current node
        pruned.isize.list[[i]] <- pruned.isizes
    }

    return(pruned.isize.list)
}



#' Converts a graph adjacency matrix to an adjacency list
#'
#' This function converts a given graph adjacency matrix into an adjacency list,
#' where each entry in the list corresponds to a node and contains the indices
#' of nodes it is connected to.
#'
#' @param M A matrix representing the adjacency matrix of a graph.
#' @return A list representing the adjacency list of the graph.
convert.adjacency.matrix.to.adjacency.list <- function(M) {

    if (!is.matrix(M)) {
        stop("M has to be a matrix.")
    }

    L <- list()
    for (i in seq(nrow(M))) {
        L[[i]] <- which(M[i,] != 0)
    }

    L
}

#' Converts a graph weighted adjacency matrix to an adjacency list and weights list
#'
#' This function converts a given weighted graph adjacency matrix into two lists:
#' 1. An adjacency list where each entry corresponds to a node and contains the indices
#'    of nodes it is connected to.
#' 2. A weights list where each entry corresponds to a node and contains the weights
#'    of the edges connecting it to other nodes.
#'
#' @param M A matrix representing the weighted adjacency matrix of a graph.
#' @return A list with two elements: `adjacency.list` and `weights.list`.
#'         `adjacency.list` is a list where each element is a vector of connected node indices.
#'         `weights.list` is a list where each element is a vector of weights corresponding to the edges in the adjacency list.
convert.weighted.adjacency.matrix.to.adjacency.list <- function(M) {

    if (!is.matrix(M)) {
        stop("M has to be a matrix.")
    }

    adjacency.list <- list()
    weights.list <- list()
    for (i in seq(nrow(M))) {
        connections <- which(M[i,] != 0)
        adjacency.list[[i]] <- connections
        weights.list[[i]] <- M[i, connections]
    }

    list(adjacency.list = adjacency.list, weights.list = weights.list)
}

#' Convert an Adjacency List to an Edgelist
#'
#' This function takes an adjacency list representation of a graph and converts it
#' into an edgelist format suitable for use with plotting functions like `igraph`.
#'
#' @param adj.list An adjacency list where each element is a vector of indices
#'   representing the neighbors of a node.
#' @return A list where each element is a vector of length two, representing
#'   an edge in the graph (source node, target node).
convert.adjacency.to.edgelist <- function(adj.list) {
    edges <- list()
    ## Initialize empty list to store edges
    j <- 1
    for (i in seq_along(adj.list)) {
        source.node <- i
        for (target.node in adj.list[[i]]) {
            edges[[j]] <- c(source.node, target.node)
            j <- j + 1
        }
    }
    return(edges)
}


#' Converts an Edge Label List to an Edge Vector
#'
#' This function takes a list of edge labels and turns it into a vector.
#'
#' @param edge.label.list A list of edge labels.
#' @param rm.duplicates Set to TRUE to allow for duplicate edges.
#' @return A vector of edge labels.
convert.edge.label.list.to.edge.label.vector <- function(edge.label.list, rm.duplicates = TRUE) {
    edge.label <- c()
    for (i in seq_along(edge.label.list)) {
        source.node <- i
        for (label in edge.label.list[[i]]) {
            edge.label <- c(edge.label, label)
        }
    }

    if (rm.duplicates) {
        edge.label <- unique(edge.label)
    }

    edge.label
}



## -------------------------------------------------------------------------------------------
##
## Distance between graphs
##
## -------------------------------------------------------------------------------------------

#' Hungarian-Frobenius Graph Matching Algorithm
#'
#' This function implements the Hungarian-Frobenius Graph Matching Algorithm to measure the similarity
#' between two graphs. It uses the Hungarian algorithm to find the optimal permutation of vertices that
#' minimizes the cost of matching the vertices between the graphs, and then calculates the Frobenius norm
#' distance between the distance matrices of the graphs after applying the optimal permutation.
#'
#' @param graph1 A list representing the adjacency list of the first graph.
#' @param graph2 A list representing the adjacency list of the second graph.
#' @return The similarity score between the two graphs, ranging from 0 to 1, where 0 indicates perfect similarity.
hungarian.frobenius.graph.matching <- function(graph1, graph2) {
  g1.m <- dgraphs::convert.adjacency.to.edge.matrix(graph1)$edge.matrix
  g2.m <- dgraphs::convert.adjacency.to.edge.matrix(graph2)$edge.matrix

  g1 <- igraph::graph_from_edgelist(g1.m, directed = FALSE)
  g2 <- igraph::graph_from_edgelist(g2.m, directed = FALSE)

  igraph::V(g1)$label <- as.character(1:igraph::vcount(g1))
  igraph::V(g2)$label <- as.character(1:igraph::vcount(g2))

  ## Calculate the minimal path distance matrices
  D1 <- igraph::distances(g1)
  D2 <- igraph::distances(g2)

  if (igraph::vcount(g1) == igraph::vcount(g2)) {
    ## Compute the cost matrix
    C <- matrix(0, nrow = igraph::vcount(g1), ncol = igraph::vcount(g2))
    for (i in 1:igraph::vcount(g1)) {
      for (j in 1:igraph::vcount(g2)) {
        C[i, j] <- sum((D1[i, ] - D2[j, ])^2)
      }
    }

    ## Find the optimal permutation using the Hungarian algorithm
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' is required for this function. Please install it with install.packages('clue')")
    }
    mapping <- clue::solve_LSAP(C)

    ## Permute the vertices of graph 2 according to the optimal mapping
    perm_g2 <- igraph::permute(g2, mapping)

    ## Recompute the distance matrix for the permuted graph 2
    perm_D2 <- igraph::distances(perm_g2)  # Corrected: Changed shortest.paths to distances

    ## Calculate the Frobenius distance between the original and permuted distance matrices
    deviation <- norm(D1 - perm_D2, type = "F")
  } else {
    ## Determine the maximum number of vertices between the two graphs
    max_vertices <- max(igraph::vcount(g1), igraph::vcount(g2))

    ## Pad the distance matrices with dummy vertices
    padded_D1 <- matrix(max_vertices, nrow = max_vertices, ncol = max_vertices)
    padded_D1[1:igraph::vcount(g1), 1:igraph::vcount(g1)] <- D1

    padded_D2 <- matrix(max_vertices, nrow = max_vertices, ncol = max_vertices)
    padded_D2[1:igraph::vcount(g2), 1:igraph::vcount(g2)] <- D2

    ## Compute the cost matrix
    C <- matrix(0, nrow = max_vertices, ncol = max_vertices)
    for (i in 1:max_vertices) {
      for (j in 1:max_vertices) {
        C[i, j] <- sum((padded_D1[i, ] - padded_D2[j, ])^2)
      }
    }

    ## Find the optimal permutation using the Hungarian algorithm
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' is required for this function. Please install it with install.packages('clue')")
    }
    mapping <- clue::solve_LSAP(C, maximum = FALSE)

      S <- intersect(mapping[1:igraph::vcount(g2)], seq(igraph::vcount(g2)))
      valid_mapping <- c(S, setdiff(seq(igraph::vcount(g2)), S))

    ## Ensure valid_mapping is within the correct range
    if (length(valid_mapping) != igraph::vcount(g2)) {
      stop("Invalid mapping length. length(valid_mapping):", length(valid_mapping), " vcount(g2): ", igraph::vcount(g2))
    }

    ## Permute the vertices of graph 2 according to the optimal mapping
    perm_g2 <- igraph::permute(g2, valid_mapping)

    ## Recompute the distance matrix for the permuted graph 2
    perm_D2 <- igraph::distances(perm_g2)  # Corrected: Changed shortest.paths to distances

    ## Calculate the Frobenius distance between the original and permuted distance matrices
    deviation <- norm(padded_D1 - padded_D2, type = "F")
  }

  ## Normalize the deviation
  n <- vcount(g1)
  max_distance <- n - 1
  max_deviation <- sqrt(n * max_distance^2)
  normalized_deviation <- deviation / max_deviation

  return(normalized_deviation)
}




#' Calculate Degree Distribution Properties for Random Points on a Sphere
#'
#' @description
#' Simulates points uniformly on a sphere and computes the degree distribution
#' properties of their k-nearest neighbor graph. For each simulation, points are
#' generated, a k-NN graph is constructed, and the proportion of vertices with
#' each degree is calculated. The function returns mean proportions and confidence
#' intervals across all simulations.
#'
#' @param n.pts numeric; Number of points to generate on the sphere for each simulation.
#' @param n.sims numeric; Number of simulations to run.
#' @param k numeric; Number of nearest neighbors to use in constructing the graph.
#' @param dim numeric; Dimension of the sphere (e.g., 2 for circle, 3 for sphere).
#'
#' @return A list containing:
#' \describe{
#'   \item{mean.props}{Vector of mean proportions for each degree}
#'   \item{ci.lower}{Vector of lower 95% confidence interval bounds}
#'   \item{ci.upper}{Vector of upper 95% confidence interval bounds}
#'   \item{degrees}{Vector of degree values corresponding to the proportions}
#' }
#'
#' @details
#' The function generates uniform random points on a sphere using \code{runif.sphere}
#' and constructs k-nearest neighbor graphs using \code{dgraphs::create.single.iknn.graph()}. For each
#' simulation, it computes the proportion of vertices with each degree. The final
#' results include means and 95% confidence intervals for these proportions across
#' all simulations.
#'
#' The confidence intervals are computed using the normal approximation:
#' mean +/- 1.96 * (standard deviation / sqrt(n.sims))
#'
#' @examples
#' \dontrun{
#' # Calculate degree distribution properties for 1000 points on a circle
#' circle_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
#'
#' # Calculate for points on a sphere
#' sphere_props <- get.sphere.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 3)
#' }
#' @export
get.sphere.degree.props <- function(n.pts = 1000, n.sims = 100, k = 10, dim = 2) {

    if(!all(sapply(list(n.pts, n.sims, k, dim), is.numeric))) {
        stop("All parameters must be numeric")
    }
    if(!all(sapply(list(n.pts, n.sims, k, dim), function(x) x > 0))) {
        stop("All parameters must be positive")
    }
    if(dim < 1) stop("dimension must be at least 1")
    if(k >= n.pts) stop("k must be less than n.pts")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        X <- runif.sphere(n.pts, dim = dim)
        graph <- dgraphs::create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)  # Added to show which degrees correspond to the proportions
}

#' Generate Degree Distribution Properties in Tubular Neighborhood of Unit Circle
#'
#' @description
#' Generates random samples in a tubular neighborhood of the unit circle and computes
#' degree distribution properties of the resulting k-nearest neighbor graph. The function
#' performs multiple simulations to estimate mean proportions and confidence intervals
#' for each degree.
#'
#' @param n.pts Positive integer. Number of points to generate in each simulation.
#' @param n.sims Positive integer. Number of simulations to run.
#' @param k Positive integer. Number of nearest neighbors for graph construction. Must be less than n.pts.
#' @param noise Non-negative numeric. Standard deviation of the noise added to the radius.
#' @param noise.type Character string. Type of noise distribution to use ("laplace" or "normal").
#'
#' @return A list containing:
#' \itemize{
#'   \item mean.props: Vector of mean proportions for each degree
#'   \item ci.lower: Vector of lower 95% confidence interval bounds
#'   \item ci.upper: Vector of upper 95% confidence interval bounds
#'   \item degrees: Vector of degrees corresponding to the proportions
#' }
#'
#' @details
#' The function uses `generate.circle.data()` to create points and `dgraphs::create.single.iknn.graph()`
#' to construct the k-nearest neighbor graph. It computes degree distributions for each
#' simulation and aggregates the results to estimate population parameters.
#'
#' @examples
#' res <- get.TN.S1.degree.props(
#'   n.pts = 100,
#'   n.sims = 10,
#'   k = 5,
#'   noise = 0.05,
#'   noise.type = "normal"
#' )
#' str(res)
#' @importFrom stats sd
#' @export
get.TN.S1.degree.props <- function(n.pts = 100,
                                   n.sims = 10,
                                   k = 10,
                                   noise = 0.1,
                                   noise.type = "laplace") {

    noise.type <- match.arg(noise.type, c("laplace", "normal"))

    if (!all(sapply(list(n.pts, n.sims, k, noise), is.numeric))) {
        stop("The first four parameters must be numeric")
    }
    if (!all(sapply(list(n.pts, n.sims, k), function(x) x > 0))) {
        stop("The first three parameters must be positive")
    }
    if (k >= n.pts) stop("k must be less than n.pts")
    if (noise < 0) stop("noise has to be non-negative")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        cX <- generate.circle.data(n.pts, radius = 1, noise = noise, type = "random", noise.type = noise.type)
        X <- cX[,1:2]
        graph <- dgraphs::create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)
}

#' Calculate Degree Distribution Properties for Random Points on a Torus
#'
#' @description
#' Simulates points uniformly on a torus and computes the degree distribution
#' properties of their k-nearest neighbor graph. For each simulation, points are
#' generated, a k-NN graph is constructed, and the proportion of vertices with
#' each degree is calculated. The function returns mean proportions and confidence
#' intervals across all simulations.
#'
#' @param n.pts numeric; Number of points to generate on the torus for each simulation.
#' @param n.sims numeric; Number of simulations to run.
#' @param k numeric; Number of nearest neighbors to use in constructing the graph.
#' @param dim numeric; Dimension of the torus (e.g., 1 for circle, 2 for 2-torus).
#'
#' @return A list containing:
#' \describe{
#'   \item{mean.props}{Vector of mean proportions for each degree}
#'   \item{ci.lower}{Vector of lower 95% confidence interval bounds}
#'   \item{ci.upper}{Vector of upper 95% confidence interval bounds}
#'   \item{degrees}{Vector of degree values corresponding to the proportions}
#' }
#'
#' @details
#' The function generates uniform random points on a torus using \code{\link{runif.torus}}
#' and constructs k-nearest neighbor graphs using \code{dgraphs::create.single.iknn.graph()}. For each
#' simulation, it computes the proportion of vertices with each degree. The final
#' results include means and 95% confidence intervals for these proportions across
#' all simulations.
#'
#' The confidence intervals are computed using the normal approximation:
#' \eqn{\text{mean} \pm 1.96 \times (\text{standard deviation} / \sqrt{n.sims})}
#'
#' @examples
#' \dontrun{
#' # Calculate degree distribution properties for 1000 points on a circle
#' circle_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 1)
#'
#' # Calculate for points on a 2-torus
#' torus_props <- get.torus.degree.props(n.pts = 1000, n.sims = 100, k = 10, dim = 2)
#' }
#' @export
get.torus.degree.props <- function(n.pts = 1000, n.sims = 100, k = 10, dim = 1) {
    if(!all(sapply(list(n.pts, n.sims, k, dim), is.numeric))) {
        stop("All parameters must be numeric")
    }
    if(!all(sapply(list(n.pts, n.sims, k, dim), function(x) x > 0))) {
        stop("All parameters must be positive")
    }
    if(round(dim) != dim) {
        stop("dimension must be a positive integer")
    }
    if(k >= n.pts) stop("k must be less than n.pts")

    # Initialize list to store degree tables
    props <- list()

    # First pass to collect degree information and find max degree
    max.degree <- 0
    for(i in 1:n.sims) {
        X <- runif.torus(n.pts, dim = dim)
        graph <- dgraphs::create.single.iknn.graph(X, k = k, compute.full = FALSE, verbose = FALSE)
        degrees <- sapply(graph$pruned_adj_list, FUN = length)
        freq.table <- table(degrees)
        props[[i]] <- freq.table / sum(freq.table)
        max.degree <- max(max.degree, max(as.numeric(names(freq.table))))
    }

    # Create and fill matrix
    deg.matrix <- matrix(0, nrow = n.sims, ncol = max.degree)
    for(i in 1:n.sims) {
        i.props <- props[[i]]
        degrees <- as.numeric(names(i.props))
        for(d in degrees) {
            deg.matrix[i, d] <- i.props[as.character(d)]
        }
    }

    # Calculate statistics
    mean.props <- colMeans(deg.matrix, na.rm = TRUE)
    sd.props <- apply(deg.matrix, 2, sd, na.rm = TRUE)
    ci.lower <- mean.props - 1.96 * sd.props/sqrt(n.sims)
    ci.upper <- mean.props + 1.96 * sd.props/sqrt(n.sims)

    list(mean.props = mean.props,
         ci.lower = ci.lower,
         ci.upper = ci.upper,
         degrees = 1:max.degree)
}

#' Create Labels for Vertices in State Space
#'
#' @description
#' Creates unique labels for specified vertices in state space by concatenating
#' two-letter shortcuts of the most abundant taxa. The function ensures uniqueness
#' of labels by iteratively including additional taxa if needed.
#'
#' @param vertices Vector of indices corresponding to rows in state.space
#' @param state.space Matrix of ASV counts/abundances. When taxonomy is NULL,
#'        column names should contain species names (e.g., "Genus_species")
#' @param taxonomy Taxonomy information for ASVs (default: NULL). If NULL,
#'        species names are taken directly from state.space column names
#' @param min.relAb.thld Minimum relative abundance threshold for including taxa in labels (default: 0.05)
#' @param profile.length Integer specifying the number of top species to keep in each profile (default: 5)
#'
#' @return A list with two components:
#'   \item{labels}{Named character vector of vertex labels}
#'   \item{profiles}{Named list of relative abundance profiles for vertices, where names are point indices and each profile is a matrix containing the top profile.length species}
#'
#' @details
#' The function processes the specified vertices, creating unique labels by:
#' 1. For each vertex, identifying taxa above the relative abundance threshold
#' 2. Creating initial labels using two-letter shortcuts
#' 3. Ensuring uniqueness by incorporating additional taxa if needed
#'
#' When taxonomy is NULL, the function uses column names of state.space directly as species names.
#' The profile.length parameter controls how many species are kept in each profile, keeping the most abundant ones.
#'
#' @examples
#' \dontrun{
#' # Keep only top 5 species in profiles
#' result <- create.vertex.labels(c(1,3,5), state.space, taxonomy,
#'                               min.relAb.thld = 0.05, profile.length = 5)
#' }
#' @export
create.vertex.labels <- function(vertices, state.space, taxonomy = NULL, min.relAb.thld = 0.05, profile.length = 5) {
    # Helper function to create two-letter shortcuts from species names
    create.shortcut <- function(sp.name) {
        parts <- strsplit(sp.name, "_")[[1]]
        if (length(parts) >= 2) {
            paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))
        } else {
            paste0(substr(parts[1], 1, 1), "")
        }
    }

    # Helper function to get profile, handling both taxonomy and direct column names cases
    get.profile <- function(id, state.space, taxonomy) {
        if (is.null(taxonomy)) {
            # Use state.space directly
            abundances <- state.space[id, ]
            profile <- cbind(
                species = colnames(state.space),
                abundance = as.numeric(abundances)
            )
            # Sort by abundance in decreasing order
            profile <- profile[order(as.numeric(profile[, 2]), decreasing = TRUE), , drop = FALSE]
        } else {
            # Use prof.fn when taxonomy is provided
            profile <- prof.fn(id, state.space, taxonomy)
        }

        # Limit profile length
        if (nrow(profile) > profile.length) {
            profile <- profile[1:profile.length, , drop = FALSE]
        }

        return(profile)
    }

    # Initialize labels vector and profiles list
    labels <- character(length(vertices))
    names(labels) <- vertices
    profiles <- list()

    # First pass: create initial labels and store profiles
    for (i in seq_along(vertices)) {
        point.i <- vertices[i]
        id <- if (is.null(taxonomy)) point.i else rownames(state.space)[point.i]

        # Get profile for this point
        profile <- get.profile(id, state.space, taxonomy)
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
                point.i <- vertices[idx]
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

#' Create Label Tables and Indicator Vectors for Vertices
#'
#' @description
#' Creates a label table and indicator vector for specified vertices in the state space.
#' The label table maps point IDs to their corresponding labels, while the indicator
#' vector marks the presence/absence of vertices in the state space.
#'
#' @param vertex.labels Named character vector of labels for vertices
#' @param state.space State space matrix where rownames correspond to point IDs
#'
#' @return A list with two components:
#'   \item{lab.tbl}{Named character vector mapping point IDs to vertex labels}
#'   \item{ind}{Numeric vector indicating vertices (1) and non-vertices (0)}
#'
#' @details
#' The function processes vertex labels to create:
#' 1. A label table that maps point IDs to their corresponding labels
#' 2. A binary indicator vector marking the presence (1) or absence (0) of vertices
#' The vector maintains the same length as the number of rows in state.space and uses consistent naming.
#'
#' @examples
#' vertex.labels <- c("1" = "A", "3" = "C")
#' state.space <- matrix(1:12, nrow = 4, dimnames = list(paste0("p", 1:4), NULL))
#' result <- create.vertex.label.indicators(vertex.labels, state.space)
#' print(head(result$lab.tbl))
#' print(sum(result$ind)) # Number of vertices
#' @export
create.vertex.label.indicators <- function(vertex.labels, state.space) {
    # Convert label names to integers
    point.indices <- as.integer(names(vertex.labels))

    # Get corresponding row names from state.space
    point.ids <- rownames(state.space)[point.indices]

    # Create label table
    lab.tbl <- c()
    lab.tbl[point.ids] <- vertex.labels[as.character(point.indices)]

    # Create indicator vector
    ind <- numeric(nrow(state.space))
    names(ind) <- rownames(state.space)
    ind[point.ids] <- 1

    # Return results
    return(list(
        lab.tbl = lab.tbl,
        ind = ind
    ))
}




#' Generate a Gaussian Mixture Function on a Graph
#'
#' @description
#' Creates a function on a graph that is a mixture of Gaussian-like components centered
#' at specified vertices. Each component follows a decay based on the shortest path
#' distance from its center.
#'
#' @param adj.list List of adjacency lists, where \code{adj.list[[i]]} contains indices of
#'        vertices adjacent to vertex \code{i}.
#' @param weight.list List of edge weights, where \code{weight.list[[i]][j]} is the weight
#'        of the edge from vertex \code{i} to \code{adj.list[[i]][j]}.
#' @param centers Vector of vertex indices to use as centers for the Gaussian components
#' @param amplitudes Vector of amplitudes for each Gaussian component
#' @param sigmas Vector of sigma values (spread) for each Gaussian component
#' @param normalize Logical; if TRUE, normalize the resulting function to \code{[0,1]}.
#'
#' @return A numeric vector of function values at each vertex of the graph
#'
#' @examples
#' adj.list <- list(c(2L), c(1L, 3L), c(2L, 4L), c(3L))
#' weight.list <- list(1, c(1, 1), c(1, 1), 1)
#' centers <- c(1, 4)
#' amplitudes <- c(1.0, 0.5)
#' sigmas <- c(1.0, 1.5)
#' generate.graph.gaussian.mixture(adj.list, weight.list, centers, amplitudes, sigmas)
#'
#' @export
generate.graph.gaussian.mixture <- function(
  adj.list,
  weight.list,
  centers,
  amplitudes = rep(1.0, length(centers)),
  sigmas = rep(2.0, length(centers)),
  normalize = TRUE
) {
  # Verify inputs
  if (length(centers) != length(amplitudes) || length(centers) != length(sigmas)) {
    stop("centers, amplitudes, and sigmas must have the same length")
  }

  n_vertices <- length(adj.list)
  n_components <- length(centers)

  # Initialize result vector
  result <- numeric(n_vertices)

  # Compute distance matrices from each center
  distance_matrices <- list()
  for (center in centers) {
    # Initialize distances for Dijkstra's algorithm
    distances <- rep(Inf, n_vertices)
    distances[center] <- 0

    # Set of vertices whose shortest distance is determined
    visited <- logical(n_vertices)

    # Priority queue (implemented as a while loop with min selection)
    while (!all(visited)) {
      # Find unvisited vertex with minimum distance
      current <- which.min(ifelse(visited, Inf, distances))

      # If all remaining vertices are unreachable, break
      if (is.infinite(distances[current])) {
        break
      }

      # Mark as visited
      visited[current] <- TRUE

      # Update distances to neighbors
      for (i in seq_along(adj.list[[current]])) {
        neighbor <- adj.list[[current]][i]
        weight <- weight.list[[current]][i]

        # Calculate potential new distance
        new_dist <- distances[current] + weight

        # Update if better path found
        if (new_dist < distances[neighbor]) {
          distances[neighbor] <- new_dist
        }
      }
    }

    # Store distance matrix
    distance_matrices[[length(distance_matrices) + 1]] <- distances
  }

  # Compute mixture of Gaussians
  for (i in 1:n_components) {
    # Compute Gaussian function using distances
    gauss_values <- amplitudes[i] * exp(-(distance_matrices[[i]]^2) / (2 * sigmas[i]^2))

    # Add to mixture
    result <- result + gauss_values
  }

  # Normalize if requested
  if (normalize && max(result) > 0) {
    result <- result / max(result)
  }

  return(result)
}

#' Visualize a Function on a Grid Graph
#'
#' Creates multiple visualizations of a function defined on a grid graph including
#' heatmap with contours, 3D perspective plot, and optionally an interactive 3D plot.
#'
#' @param grid.size Integer; size of the square grid (grid.size x grid.size).
#' @param z Numeric vector; function values at each vertex in row-wise order (length = grid.size^2).
#' @param centers Optional integer vector; vertex indices to highlight (e.g., peaks).
#' @param title Character string; title for the plots (default: "Function on Grid Graph").
#'
#' @return Invisibly returns the input \code{z}.
#'
#' @details
#' Left panel: heatmap with contour lines (and optional centers).
#' Right panel: 3D perspective plot (base graphics).
#' If \pkg{rgl} is available, a separate off-screen rgl window is opened to render a 3D point view.
#'
#' @examples
#' \dontrun{
#' grid.size <- 20
#' x <- rep(1:grid.size, grid.size) / grid.size
#' y <- rep(1:grid.size, each = grid.size) / grid.size
#' z <- exp(-10*((x-0.3)^2 + (y-0.3)^2)) + 0.5*exp(-8*((x-0.7)^2 + (y-0.6)^2))
#' centers <- which(z > 0.9)
#' visualize.grid.function(grid.size, z, centers)
#'
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'   old <- options(rgl.useNULL = TRUE); on.exit(options(old), add = TRUE)
#'   visualize.grid.function(grid.size, z, centers)
#' }
#' }
#' @export
visualize.grid.function <- function(grid.size, z, centers = NULL, title = "Function on Grid Graph") {
    ## ---- Validation -----------------------------------------------------------
    if (length(grid.size) != 1L || !is.numeric(grid.size) || grid.size < 2 || !is.finite(grid.size)) {
        stop("'grid.size' must be a single finite numeric value >= 2", call. = FALSE)
    }
    grid.size <- as.integer(grid.size)
    n_vertices <- grid.size^2L

    if (!is.numeric(z) || length(z) != n_vertices) {
        stop(sprintf("'z' must be a numeric vector of length grid.size^2 (= %d)", n_vertices), call. = FALSE)
    }
    if (any(!is.finite(z))) {
        stop("'z' contains non-finite values; please remove or impute NA/Inf/NaN", call. = FALSE)
    }

    if (!is.null(centers)) {
        centers <- as.integer(centers)
        centers <- centers[is.finite(centers)]
        centers <- centers[centers >= 1L & centers <= n_vertices]
        if (!length(centers)) centers <- NULL
    }

    ## ---- Coordinates & helpers -----------------------------------------------
    x_coords <- rep(seq_len(grid.size), grid.size) / grid.size
    y_coords <- rep(seq_len(grid.size), each = grid.size) / grid.size

    vertex_to_coords <- function(vertices) {
        xv <- (vertices - 1L) %% grid.size + 1L
        yv <- ceiling(vertices / grid.size)
        list(x = xv / grid.size, y = yv / grid.size)
    }

    center_coords <- if (!is.null(centers)) vertex_to_coords(centers) else NULL

    ## ---- Set layout and par safely -------------------------------------------
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)

    graphics::layout(matrix(c(1, 2), nrow = 1L, ncol = 2L))
    on.exit(try(graphics::layout(1), silent = TRUE), add = TRUE)

    ## ---- Plot 1: Heatmap + contours ------------------------------------------
    graphics::par(mar = c(4, 4, 2, 1))
    z_matrix <- matrix(z, nrow = grid.size, ncol = grid.size, byrow = FALSE)

    pal2d <- grDevices::hcl.colors(100, palette = "YlOrRd")
    graphics::image(
                  x = seq_len(grid.size) / grid.size,
                  y = seq_len(grid.size) / grid.size,
                  z = z_matrix,
                  col = pal2d,
                  main = title,
                  xlab = "X", ylab = "Y"
              )

    graphics::contour(
                  x = seq_len(grid.size) / grid.size,
                  y = seq_len(grid.size) / grid.size,
                  z = z_matrix,
                  add = TRUE,
                  col = "black"
              )

    if (!is.null(centers)) {
        graphics::points(center_coords$x, center_coords$y, pch = 19, col = "blue", cex = 1.5)
        graphics::text(center_coords$x, center_coords$y,
                       labels = paste("Center", seq_along(centers)),
                       pos = 3, offset = 0.7, cex = 0.8)
    }

    ## ---- Plot 2: Base-graphics 3D perspective -------------------------------
    graphics::par(mar = c(4, 4, 2, 1))
    graphics::persp(
                  x = seq_len(grid.size) / grid.size,
                  y = seq_len(grid.size) / grid.size,
                  z = z_matrix,
                  theta = 30, phi = 30,
                  expand = 0.7,
                  col = "lightblue",
                  shade = 0.5,
                  main = "3D Perspective",
                  xlab = "X", ylab = "Y", zlab = "Value"
              )

    ## ---- Optional rgl 3D view (headless-safe) --------------------------------
    if (requireNamespace("rgl", quietly = TRUE)) {

        use_null <- (!interactive()) ||
            identical(Sys.getenv("RGL_USE_NULL"), "TRUE") ||
            (Sys.getenv("DISPLAY") == "" && .Platform$OS.type != "windows")
        old_opt <- options(rgl.useNULL = use_null)
        on.exit(options(old_opt), add = TRUE)

        rgl::open3d()
        if (use_null) {
            ## Only close the device automatically if using null device
            on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
        }
        rgl::clear3d()

        pal3d <- grDevices::hcl.colors(length(z), palette = "Spectral")
        rgl::plot3d(
                 x_coords, y_coords, z,
                 col = pal3d[rank(z, ties.method = "average")],
                 size = 3,
                 xlab = "X", ylab = "Y", zlab = "Z",
                 type = "p"
             )
        rgl::title3d(main = title)

        if (!is.null(centers)) {
            for (i in seq_along(centers)) {
                rgl::spheres3d(
                         center_coords$x[i], center_coords$y[i], z[centers[i]],
                         radius = 0.02,
                         color = "blue"
                     )
            }
        }

        rgl::axes3d()
    }

    invisible(z)
}

#' Assess Fidelity of Graph-Based Geodesic Distances to Euclidean Geometry
#'
#' This function compares local neighborhood structures defined by graph-based distances
#' to those defined in the original Euclidean space. It evaluates how well a given graph
#' preserves the local geometry of a dataset by computing two metrics for each data point:
#'
#' \itemize{
#'   \item \strong{Jaccard Index:} Measures the similarity between the sets of neighbors
#'         within a radius \eqn{\tau} in both the Euclidean and graph-based spaces.
#'   \item \strong{Mean Absolute Deviation (MAD):} Measures the distortion in local distances
#'         over the intersection of the two neighborhood sets.
#' }
#'
#' The graph is assumed to be provided as an adjacency list and a corresponding weight list.
#' Any duplicate edges (i.e., undirected edges appearing in both directions) are automatically deduplicated.
#'
#' @param X A numeric matrix of shape \code{[n, d]}, where each row represents a data point in \code{d}-dimensional space.
#' @param adj.list A list of integer vectors of length \code{n}, giving the graph adjacency list. Each entry contains the 1-based indices of neighbors for that vertex.
#' @param weight.list A list of numeric vectors of same length as \code{adj.list}, where each element contains edge weights for the corresponding neighbors.
#' @param taus A numeric vector of radius values \eqn{\tau} used to define local neighborhoods around each point.
#' @param max.k Integer; the number of nearest neighbors to compute in the Euclidean space (used to approximate local neighborhoods).
#'
#' @return A named list of length equal to \code{length(taus)}. Each element is a list with:
#' \describe{
#'   \item{\code{mean.jaccard}}{Mean Jaccard index over all data points for a given \eqn{\tau}.}
#'   \item{\code{mean.mad}}{Mean absolute deviation in distances between Euclidean and graph metrics over intersecting neighbors.}
#'   \item{\code{jaccard.vals}}{A numeric vector of Jaccard index values per vertex.}
#'   \item{\code{mad.vals}}{A numeric vector of MAD values per vertex.}
#' }
#'
#' @examples
#' X <- matrix(
#'   c(0, 0,
#'     1, 0,
#'     0, 1,
#'     1, 1),
#'   ncol = 2,
#'   byrow = TRUE
#' )
#' adj.list <- list(c(2L, 3L), c(1L, 4L), c(1L, 4L), c(2L, 3L))
#' weight.list <- list(c(1, 1), c(1, 1), c(1, 1), c(1, 1))
#' taus <- c(1.1, 1.5)
#' fidelity <- compute.local.distance.fidelity(X, adj.list, weight.list, taus, max.k = 2)
#' names(fidelity)
#'
#' @importFrom FNN get.knn
#' @importFrom igraph graph_from_adj_list as_edgelist graph_from_data_frame E distances
#'
#' @export
compute.local.distance.fidelity <- function(X, adj.list, weight.list, taus, max.k = 50) {
    n <- nrow(X)

    ## 1. Euclidean distances
    nn.X <- FNN::get.knn(X, k = max.k)
    D.X <- as.matrix(dist(X))  ## For true Euclidean distances

    ## 2. Graph distances
    G <- igraph::graph_from_adj_list(adj.list, mode="all")
    ## Get edge list and remove duplicates
    edgelist <- igraph::as_edgelist(G)
    edge.df <- data.frame(from = pmin(edgelist[,1], edgelist[,2]),
                          to   = pmax(edgelist[,1], edgelist[,2]),
                          weight = unlist(weight.list))
    ## Remove duplicates
    edge.df <- edge.df[!duplicated(edge.df), ]
    ## Rebuild the graph and assign deduplicated weights
    G <- igraph::graph_from_data_frame(edge.df, directed = FALSE)
    igraph::E(G)$weight <- edge.df$weight

    D.G <- igraph::distances(G, v=1:n, to=1:n, weights=E(G)$weight)

    results <- list()
    for (tau.idx in seq(taus)) {
        tau <- taus[tau.idx]
        mad_vals <- numeric(n)
        jaccard_vals <- numeric(n)

        for (i in 1:n) {
            A <- which(D.X[i, ] < tau)
            B <- which(D.G[i, ] < tau)

            I <- intersect(A, B)
            U <- union(A, B)

            if (length(U) > 0) {
                jaccard_vals[i] <- length(I) / length(U)
            }

            if (length(I) > 0) {
                d_x <- D.X[i, I]
                d_g <- D.G[i, I]
                mad_vals[i] <- mean(abs(d_x - d_g))
            }
        }

        results[[tau.idx]] <- list(
            mean.jaccard = mean(jaccard_vals, na.rm = TRUE),
            mean.mad = mean(mad_vals, na.rm = TRUE),
            jaccard.vals = jaccard_vals,
            mad.vals = mad_vals
        )
    }

    return(results)
}

compute.kernel.graph.laplacian.eigenfunctions.I.minus.L.powered <- function(x.min,
                                                                            x.max,
                                                                            grid.size,
                                                                            x.grid = NULL,
                                                                            random = TRUE,
                                                                            n.evects,
                                                                            tau,
                                                                            power = 1) {
    if (is.null(x.grid)) {
        if (random) {
            x.grid <- sort(runif(grid.size, min = x.min, max = x.max))
        } else {
            x.grid <- seq(x.min, x.max, length.out = grid.size) ## Uniform grid
        }
    }

    ## Compute pairwise squared distances
    distsq.mat <- outer(x.grid, x.grid, function(a, b) (a - b)^2)

    ## Construct kernel-based affinity matrix W
    W <- exp(-distsq.mat / tau^2)
    diag(W) <- 0  # No self-loops

    ## Degree matrix D
    D <- diag(rowSums(W))

    ## Unnormalized Graph Laplacian
    L <- D - tau * W

    ## Construct (I - L)
    I_minus_L <- diag(grid.size) - L

    ## Raise (I - L) to the given power
    L_mod <- I_minus_L
    if (power > 1) {
        for (i in 2:power) {
            L_mod <- L_mod %*% I_minus_L
        }
    }

    ## Eigen decomposition
    eig <- eigen(L_mod, symmetric = TRUE)

    list(
        x.grid = x.grid,
        eigenvectors = eig$vectors[, 1:n.evects, drop = FALSE],
        eigenvalues = eig$values[1:n.evects],
        L_mod = L_mod,
        W = W,
        D = D,
        power = power
    )
}

## compare.kernel.graph.laplacian.with.continuum.laplacian <- function(result, x.min, x.max, n.evects = NULL) {
##     x.grid <- result$x.grid

##     if (is.null(n.evects)) {
##         n.evects <- ncol(result$eigenvectors)
##     }

##     if (n.evects == 1) {
##         op <- par(mfrow=c(1,1))
##     } else {
##         op <- par(mfrow = c(ceiling(n.evects / 2), 2), mar = c(2, 2, 2, 1))
##     }

##     for (j in 1:n.evects) {
##         ## Classical continuum Laplacian eigenfunctions (sine waves)
##         phi.j <- sin(j * pi * (x.grid - x.min) / (x.max - x.min))

##         ## Normalize both to have unit norm (L2)
##         phi.j <- phi.j / sqrt(sum(phi.j^2))
##         evect.j <- result$eigenvectors[, j]
##         evect.j <- evect.j / sqrt(sum(evect.j^2))

##         plot(x.grid, evect.j, type = "l", col = "blue", lwd = 2, las = 1,
##              ylab = "", xlab = "", main = paste("Eigenfunction", j))
##         lines(x.grid, phi.j, col = "red", lty = 2, lwd = 2)
##         ## legend("topright", legend = c(paste0(\"(I-L)^\", result$power, \" eigen\"), \"Classical (sine)\"),
##         ##        col = c(\"blue\", \"red\"), lty = c(1, 2), cex = 0.8, inset = 0.025)
##     }

##     par(op)
## }

#' Create Threshold Distance Graph from Distance Matrix
#'
#' This function constructs a graph where vertices are connected by an edge
#' if and only if their distance is less than a specified threshold.
#' The weight of each edge is the corresponding distance.
#'
#' @param dist.matrix A symmetric distance matrix where rows and columns correspond to vertices
#' @param threshold A numeric threshold value; vertices i and j are connected if dist(i,j) < threshold
#' @param include.names Logical; whether to include vertex names in the output (default: TRUE)
#'
#' @return A list with two components:
#'   \item{adj_list}{A list of integer vectors. Each vector contains the indices of vertices
#'         adjacent to the corresponding vertex.}
#'   \item{weight_list}{A list of numeric vectors. Each vector contains weights of edges
#'         corresponding to adjacencies in adj_list.}
#'
#' Compute local cluster evenness with full cluster support
#'
#' Computes the evenness (normalized entropy) of cluster label distribution
#' in the expanded neighborhood of each vertex. The entropy is computed using
#' the full set of cluster labels, including those not present in the neighborhood
#' (i.e., with zero frequency), to ensure comparability across vertices.
#'
#' @param adj.list List of adjacency vectors (1-based indices) for each vertex.
#' @param cltr A vector of cluster labels (factor, character, or numeric), one per vertex.
#'
#' @return Numeric vector of evenness values for each vertex.
graph.cltr.evenness <- function(adj.list, cltr) {
    n.vertices <- length(adj.list)
    uq.clusters <- unique(cltr)

    cltr.evenness <- numeric(n.vertices)
    for (i in seq_len(n.vertices)) {
        nbrs <- c(i, adj.list[[i]])
                                        # Include all clusters, even if count is 0
        cltr.labels <- factor(cltr[nbrs], levels = uq.clusters)
        nbrs.cltr.freq <- table(cltr.labels)
        cltr.evenness[i] <- evenness(as.numeric(nbrs.cltr.freq))
    }

    return(cltr.evenness)
}


#' Select a geodesic corridor (thickened path) in a weighted graph
#'
#' @description
#' Given a vertex path (sequence of vertex ids), selects vertices that lie on or
#' near the geodesic segment between consecutive path vertices, using the metric
#' condition d(u,x) + d(x,v) <= d(u,v) * (1 + rel.tol) + abs.tol.
#'
#' @param graph An igraph object (undirected or directed).
#' @param path Integer vector of vertex ids defining the polyline on the graph.
#' @param weights Numeric vector of edge weights interpreted as lengths/distances
#'   (same order as \code{E(graph)}), or \code{NULL} for unweighted.
#' @param rel.tol Non-negative numeric scalar. Relative corridor thickness.
#'   Typical values: 0.02 to 0.15.
#' @param abs.tol Non-negative numeric scalar. Absolute corridor thickness in the
#'   same units as \code{weights}.
#' @param include.path Logical; if TRUE, always include the original \code{path}.
#'
#' @return Integer vector of selected vertex ids.
select.path.corridor <- function(graph,
                                 path,
                                 weights = NULL,
                                 rel.tol = 0.05,
                                 abs.tol = 0,
                                 include.path = TRUE) {
  ## Basic checks
  stopifnot(igraph::is_igraph(graph))
  stopifnot(is.numeric(path), length(path) >= 2)
  stopifnot(rel.tol >= 0, abs.tol >= 0)

  n.v <- igraph::vcount(graph)
  stopifnot(all(path >= 1L), all(path <= n.v))

  path.u <- unique(as.integer(path))

  ## Distances from each path vertex to all vertices (matrix: |path.u| x n.v)
  d.to.all <- igraph::distances(graph,
                                v = path.u,
                                to = igraph::V(graph),
                                weights = weights)

  ## Distances among path vertices (matrix: |path.u| x |path.u|)
  d.between <- igraph::distances(graph,
                                 v = path.u,
                                 to = path.u,
                                 weights = weights)

  keep <- rep(FALSE, n.v)

  for (i in seq_len(length(path) - 1L)) {
    u <- as.integer(path[i])
    v <- as.integer(path[i + 1L])

    iu <- match(u, path.u)
    iv <- match(v, path.u)

    d.u <- as.numeric(d.to.all[iu, ])
    d.v <- as.numeric(d.to.all[iv, ])
    d.uv <- as.numeric(d.between[iu, iv])

    ## Corridor condition
    d.sum <- d.u + d.v
    ok <- is.finite(d.sum) & is.finite(d.uv) &
      (d.sum <= d.uv * (1 + rel.tol) + abs.tol)

    keep[ok] <- TRUE
  }

  if (isTRUE(include.path)) keep[as.integer(path)] <- TRUE
  which(keep)
}

#' Select a graph-metric tube around a path
#'
#' @description
#' Selects vertices within \code{max.dist} (graph shortest-path distance) of any
#' vertex in \code{path}.
#'
#' @param graph An igraph object.
#' @param path Integer vector of vertex ids.
#' @param weights Edge weights as in \code{igraph::distances}.
#' @param max.dist Non-negative numeric scalar; tube radius in graph-distance units.
#'
#' @return Integer vector of selected vertex ids.
select.path.neighborhood <- function(graph,
                                    path,
                                    weights = NULL,
                                    max.dist) {
  stopifnot(igraph::is_igraph(graph))
  stopifnot(is.numeric(path), length(path) >= 1)
  stopifnot(is.numeric(max.dist), length(max.dist) == 1, max.dist >= 0)

  path.u <- unique(as.integer(path))
  d.mat <- igraph::distances(graph,
                             v = path.u,
                             to = igraph::V(graph),
                             weights = weights)
  d.min <- apply(d.mat, 2, min, na.rm = TRUE)
  which(is.finite(d.min) & d.min <= max.dist)
}

#' Distance from points to a 3D polyline (minimum point-to-segment distance)
#'
#' @param x Numeric matrix n x 3 of point coordinates.
#' @param poly Numeric matrix m x 3 of polyline vertices (in order).
#'
#' @return Numeric vector length n with distances to the polyline.
dist.to.polyline.3d <- function(x, poly) {
  stopifnot(is.matrix(x), ncol(x) == 3)
  stopifnot(is.matrix(poly), ncol(poly) == 3, nrow(poly) >= 2)

  n <- nrow(x)
  d.min <- rep(Inf, n)

  for (i in seq_len(nrow(poly) - 1L)) {
    a <- poly[i, ]
    b <- poly[i + 1L, ]
    ab <- b - a
    ab2 <- sum(ab * ab)

    ## Degenerate segment guard
    if (ab2 <= 0) next

    ap <- sweep(x, 2, a, "-")
    t <- drop(ap %*% ab) / ab2
    t <- pmax(0, pmin(1, t))

    proj <- sweep(matrix(t, n, 3) * rep(ab, each = n), 2, a, "+")
    d <- sqrt(rowSums((x - proj) * (x - proj)))

    d.min <- pmin(d.min, d)
  }

  d.min
}
