## -------------------------------------------------------------------------------------------
##
## Graph construction functions
##
## -------------------------------------------------------------------------------------------

#' Constructs the 1-skeleton (nerve graph) of the nerve simplicial complex associate with a covering of some finite set
#'
#' This function constructs the 1-skeleton of the nerve of a given covering of
#' a set. The covering is represented as a list of components, each being an
#' integer index of the elements of the set. The function returns an adjacency
#' list representing the graph. Parallel processing can be utilized by
#' specifying the number of cores.
#'
#' @param covering.list A list where each component is an integer vector
#'     representing the indices of elements in a subset of the set being
#'     covered. The list represents a covering of the set.
#' @param n.cores An integer specifying the number of cores for parallel
#'     processing. If NULL, it defaults to `detectCores() - 1`. If set to 1,
#'     parallel processing is not used. Based on limited experimental data,
#'     the parallel processing becomes a viable option when the length of the
#'     list is more than 100.
#' @return A list containing:
#'   \item{adjacency.list}{The adjacency list representation of the nerve graph}
#'   \item{weights.list}{The weights (intersection sizes) for each edge}
#'   \item{adjacency.matrix}{The sparse adjacency matrix with intersection sizes as entries}
#' @details
#' The nerve of a covering is a simplicial complex where each set in the covering
#' corresponds to a vertex, and vertices are connected by an edge if and only if
#' the corresponding sets have a non-empty intersection. This function constructs
#' the 1-skeleton (graph) of this nerve complex, with edge weights equal to the
#' size of the intersection between sets.
#' @examples
#' covering.list <- list(c(1, 2), c(2, 3), c(1, 3, 4))
#' nerve.graph(covering.list, n.cores = 2)
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
#' @export
nerve.graph <- function(covering.list, n.cores = 1) {
  n <- length(covering.list)
  if (n < 2L) {
    return(list(
      adjacency.list = vector("list", n),
      weights.list   = vector("list", n),
      adjacency.matrix = if (requireNamespace("Matrix", quietly = TRUE)) {
        Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
      } else {
        matrix(0, nrow = n, ncol = n)
      }
    ))
  }

  if (is.null(n.cores)) {
    n.cores <- max(1L, parallel::detectCores() - 1L)
  }

  # Build a dense adjacency first (fast/simple), convert to sparse at the end if Matrix is present
  adj <- matrix(0, nrow = n, ncol = n)

  if (n.cores == 1L) {
    for (i in seq_len(n - 1L)) {
      ci <- covering.list[[i]]
      for (j in (i + 1L):n) {
        w <- length(intersect(ci, covering.list[[j]]))
        if (w) adj[i, j] <- adj[j, i] <- w
      }
    }
  } else {
    cl <- parallel::makeCluster(n.cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(cl, c("covering.list"), envir = environment())
    results <- parallel::parLapply(cl, 1:(n - 1L), function(i) {
      ci <- covering.list[[i]]
      v <- integer(n)
      for (j in (i + 1L):n) {
        v[j] <- length(intersect(ci, covering.list[[j]]))
      }
      v
    })
    for (i in seq_len(n - 1L)) {
      idx <- which(results[[i]] > 0L)
      if (length(idx)) {
        adj[i, idx] <- results[[i]][idx]
        adj[idx, i] <- results[[i]][idx]
      }
    }
  }

  # Convert to sparse if Matrix is available
  A <- if (requireNamespace("Matrix", quietly = TRUE)) {
    Matrix::Matrix(adj, sparse = TRUE)
  } else {
    adj
  }

  # Build adjacency list + weights from dense adj (works for both)
  L <- convert.weighted.adjacency.matrix.to.adjacency.list(adj)

  list(adjacency.list = L$adjacency.list,
       weights.list   = L$weights.list,
       adjacency.matrix = A)
}

#' Create a Complete Graph
#'
#' @description
#' Generates a complete graph where every node is connected to every other node.
#'
#' @param n An integer specifying the number of nodes in the graph.
#'
#' @return A list representing the adjacency list of the complete graph.
#'         Each element of the list corresponds to a node and contains
#'         a vector of all other nodes it's connected to.
#'
#' @details
#' In a complete graph of n nodes, each node has n-1 connections (to all other nodes).
#' The resulting graph is undirected.
#'
#' @examples
#' complete_graph <- create.complete.graph(5)
#'
#' @export
create.complete.graph <- function(n) {
  graph <- vector("list", n)
  for (i in 1:n) {
    graph[[i]] <- setdiff(1:n, i)
  }
  graph
}

#' Create an Empty Graph
#'
#' @description
#' Generates an empty graph with a specified number of nodes but no edges.
#'
#' @param n An integer specifying the number of nodes in the graph.
#'
#' @return A list representing the adjacency list of the empty graph.
#'         Each element of the list is an empty vector, indicating no connections.
#'
#' @details
#' An empty graph contains nodes but no edges connecting these nodes.
#' This function is useful for initializing graphs or testing edge cases in graph algorithms.
#'
#' @examples
#' empty_graph <- create.empty.graph(5)
#'
#' @export
create.empty.graph <- function(n) {
  replicate(n, integer(0), simplify = FALSE)
}

#' Create a Bipartite Graph
#'
#' @description
#' Generates a bipartite graph with two disjoint sets of nodes.
#'
#' @param n1 An integer specifying the number of nodes in the first set.
#' @param n2 An integer specifying the number of nodes in the second set.
#'
#' @return A list representing the adjacency list of the bipartite graph.
#'         The first n1 elements correspond to nodes in the first set,
#'         and the next n2 elements correspond to nodes in the second set.
#'
#' @details
#' In a bipartite graph, nodes are divided into two disjoint sets,
#' and each edge connects a node in the first set to a node in the second set.
#' There are no edges between nodes in the same set.
#'
#' @examples
#' bipartite_graph <- create.bipartite.graph(3, 4)
#'
#' @export
create.bipartite.graph <- function(n1, n2) {
  graph <- vector("list", n1 + n2)
  for (i in 1:n1) {
    graph[[i]] <- (n1+1):(n1+n2)
  }
  for (i in (n1+1):(n1+n2)) {
    graph[[i]] <- 1:n1
  }
  graph
}


#' Creates a bi-k-NN chain graph
#'
#' This function constructs a bi-directional k-nearest neighbor chain graph.
#' Each vertex is connected to up to k neighbors on both sides, forming a
#' chain-like structure. The graph can be based on sorted x-coordinates if provided.
#'
#' @param n.vertices An integer. The number of vertices in the graph.
#' @param k An integer. The number of neighbors to connect on each side of a vertex.
#' @param x An optional numeric vector of length n.vertices. If provided, vertices
#'   will be ordered based on these values.
#' @param y An optional numeric vector of length n.vertices. If provided along with x,
#'   it will be sorted according to x.
#'
#' @return A list containing:
#'   \item{adj.list}{A list representing the graph structure. Each element is a vector
#'     of indices of neighboring vertices.}
#'   \item{edge.lengths}{A list of edge lengths corresponding to the graph structure.
#'     If x is provided, these are based on differences in x values; otherwise, they
#'     are set to 1.}
#'   \item{x.sorted}{The sorted x vector if x was provided, otherwise NULL.}
#'   \item{y.sorted}{The y vector sorted according to x if both x and y were provided,
#'     otherwise NULL.}
#'
#' @examples
#' # Create a 5-chain graph with 10 vertices
#' graph <- create.bi.kNN.chain.graph(10, 5)
#'
#' # Create a graph based on x-coordinates
#' x <- runif(10)
#' y <- sin(x)
#' graph.with.coords <- create.bi.kNN.chain.graph(10, 2, x, y)
#'
#' @export
create.bi.kNN.chain.graph <- function(n.vertices = 5, k = 1, x = NULL, y = NULL) {

    if (!is.null(x)) {
        n.vertices <- length(x)
    }

    if (n.vertices < 2) {
        stop("A chain graph must have at least two vertices.")
    }


    adj.list <- vector("list", n.vertices)
    edge.lengths <- vector("list", n.vertices)
    x.sorted <- NULL
    y.sorted <- NULL

    if (!is.null(x)) {
        o <- order(x)
        x.sorted <- x[o]
        if (!is.null(y)) {
            if (length(y) != n.vertices) {
                stop("Length of y must equal n.vertices.")
            }
            y.sorted <- y[o]
        }
    }

    for (i in seq_len(n.vertices)) {
        neighbors <- max(1, i - k):min(n.vertices, i + k)
        neighbors <- setdiff(neighbors, i)
        adj.list[[i]] <- neighbors

        if (!is.null(x.sorted)) {
            edge.lengths[[i]] <- as.numeric(abs(x.sorted[i] - x.sorted[neighbors]))
        } else {
            edge.lengths[[i]] <- rep(1.0, length(neighbors))
        }
    }

    list(adj.list = adj.list,
         edge.lengths = edge.lengths,
         x.sorted = x.sorted,
         y.sorted = y.sorted)
}

#' Create a Chain Graph Structure
#'
#' Generates a chain graph where each vertex is connected to its immediate neighbors,
#' forming a linear chain structure. The graph can be constructed based on either
#' the number of vertices or a set of ordered x-coordinates.
#'
#' @param n.vertices Integer specifying the number of vertices. Required if x is NULL.
#' @param x Optional numeric vector. If provided, vertices will be ordered based on
#'          these values, and n.vertices will be set to length(x).
#' @param y Optional numeric vector corresponding to x values. Must be the same
#'          length as x if provided.
#'
#' @return A list containing:
#' \describe{
#'   \item{adj.list}{List of adjacent vertices for each vertex}
#'   \item{edge.lengths}{List of edge lengths between adjacent vertices (1 if x
#'                       not provided, otherwise based on x-coordinate differences)}
#'   \item{x.sorted}{Sorted x values if x was provided, otherwise NULL}
#'   \item{y.sorted}{y values sorted according to x order if both provided,
#'                   otherwise NULL}
#' }
#'
#' @examples
#' # Create a simple 4-vertex chain
#' chain1 <- create.chain.graph(4)
#' @export
create.chain.graph <- function(n.vertices = NULL, x = NULL, y = NULL) {
    # Handle x input and validation
    if (!is.null(x)) {
        n.vertices <- length(x)
        o <- order(x)
        x.sorted <- x[o]
        y.sorted <- if (!is.null(y)) {
            if (length(y) != n.vertices) stop("Length of y must equal length of x.")
            y[o]
        } else NULL
    } else {
        if (is.null(n.vertices)) stop("Either n.vertices or x must be provided.")
        x.sorted <- NULL
        y.sorted <- NULL
    }

    # Validate n.vertices
    if (n.vertices < 2) stop("A chain has to have at least two vertices.")

    # Initialize lists for all vertices
    adj.list <- vector("list", n.vertices)
    edge.lengths <- vector("list", n.vertices)

    # Set first and last vertices
    adj.list[[1]] <- 2
    adj.list[[n.vertices]] <- n.vertices - 1

    # Set middle vertices if they exist
    if (n.vertices > 2) {
        middle_vertices <- 2:(n.vertices - 1)
        adj.list[middle_vertices] <- lapply(middle_vertices, function(i) c(i - 1, i + 1))
    }

    # Calculate edge lengths
    if (!is.null(x.sorted)) {
        # Edge lengths based on x values
        edge_calc <- function(i) abs(x.sorted[i] - x.sorted[adj.list[[i]]])
        edge.lengths <- lapply(seq_len(n.vertices), edge_calc)
    } else {
        # Unit edge lengths
        edge.lengths <- lapply(adj.list, function(neighbors) rep(1, length(neighbors)))
    }

    list(
        adj.list = adj.list,
        edge.lengths = edge.lengths,
        x.sorted = x.sorted,
        y.sorted = y.sorted
    )
}

#' Creates a chain graph with n vertices
#'
#' @param n The number of vertices in the graph
#' @param offset An offset in indexing the vertices of the graph.
#'
#' @return A chain graph.
#' @export
create.chain.graph.with.offset <- function(n, offset = 0) {

    if(n < 2) {
        stop("A chain has to have at least two vertices.")
    }

    graph <- list()
    for (i in seq(n - 1)) {
        graph[[i + offset]] <- c(i + 1 + offset)
    }
    graph[[n + offset]] <- c()

    convert.to.undirected(graph)
}


#' Creates a circular graph with n vertices
#'
#' @param n The number of vertices in the graph
#'
#' @return A circular graph.
#' @export
create.circular.graph <- function(n) {

    graph <- list()
    for (i in seq(n - 1)) {
        graph[[i]] <- c(i + 1)
    }
    graph[[n]] <- c(1)

    convert.to.undirected(graph)
}


#' Join two graphs at specified vertices
#'
#' This function joins two graphs by identifying a vertex from each graph. The resulting graph
#' combines the vertices and edges from both input graphs, with the specified vertices being
#' treated as a single vertex in the joined graph.
#'
#' @param graph1 A list representing the adjacency list of the first graph.
#' @param graph2 A list representing the adjacency list of the second graph.
#' @param i1 An integer specifying the index of the vertex in the first graph to be joined.
#' @param i2 An integer specifying the index of the vertex in the second graph to be joined.
#' @return A list representing the adjacency list of the joined graph.
#'
#' @examples
#' graph1 <- list(c(2, 3), c(1), c(1))
#' graph2 <- list(c(2), c(1, 3), c(2))
#' joined_graph <- join.graphs(graph1, graph2, 2, 1)
#' print(joined_graph)
#'
#' @export
join.graphs <- function(graph1, graph2, i1, i2) {

    if (!is.list(graph1) || !is.list(graph2)) {
        stop("graph1 and graph2 must be lists representing adjacency lists.")
    }

    if (!is.numeric(i1) || !is.numeric(i2) || length(i1) != 1 || length(i2) != 1) {
        stop("i1 and i2 must be single integer values.")
    }

    if (i1 < 1 || i1 > length(graph1) || i2 < 1 || i2 > length(graph2)) {
        stop("Invalid vertex indices. i1 and i2 must be within the range of vertices in graph1 and graph2, respectively.")
    }

    ## Converting each component of graph1 and graph2 to an integer vector
    graph1 <- lapply(graph1, function(x) as.integer(x - 1))
    graph2 <- lapply(graph2, function(x) as.integer(x - 1))

    joined.graph <- .Call(S_join_graphs,
                    graph1,
                    graph2,
                    as.integer(i1 - 1),
                    as.integer(i2 - 1))

    joined.graph <- lapply(joined.graph, function(x) as.integer(x + 1))

    return(joined.graph)
}


#' Generate a Circle Graph
#'
#' Creates a circle graph where each vertex is connected to the vertices to its left and right.
#' The edge weight is the angular distance (shortest angle) between connected vertices.
#'
#' @param n A positive integer specifying the number of vertices in the graph.
#' @param type Character string specifying the type of angle distribution. Either "uniform" for
#'             equally spaced angles or "random" for random angles. Default is "random".
#' @param seed Optional seed for the random number generator. Default is NULL.
#'
#' @return A list with two components:
#'   \item{adj.list}{Adjacency list of the circle graph. Each element \code{adj.list[[i]]} contains
#'                   the vertices adjacent to vertex i.}
#'   \item{weight.list}{Weight list corresponding to the adjacency list. Each element
#'                      \code{weight.list[[i]]} contains the weights of the edges connecting
#'                      vertex i to the vertices in \code{adj.list[[i]]}.}
#'
#' @examples
#' # Generate a circle graph with 5 vertices and uniform angles
#' graph <- generate.circle.graph(5, type = "uniform")
#'
#' # Generate a circle graph with 10 vertices and random angles
#' graph <- generate.circle.graph(10, type = "random", seed = 123)
#'
#' @export
generate.circle.graph <- function(n,
                                  type = "random",
                                  seed = NULL) {
    ## Parameter checks
    if (!is.numeric(n) || n <= 0 || n != round(n))
        stop("n must be a positive integer.")
    type <- match.arg(type, c("uniform", "random"))
    ## Set seed if provided
    if (!is.null(seed)) set.seed(seed)
    ## Generate angles
    if (type == "uniform") {
        angles <- seq(0, 2 * pi, length.out = n + 1)[-1]
    } else {
        angles <- sort(stats::runif(n, min = 0, max = 2 * pi))
    }

    ## Construct adjacency list for a circle graph
    adj.list <- vector("list", n)
    for (i in 1:n) {
        adj.list[[i]] <- c(ifelse(i == 1, n, i - 1), ifelse(i == n, 1, i + 1))
    }
    names(adj.list) <- 1:n

    ## Construct weight list corresponding to adjacency list
    weight.list <- vector("list", n)
    for (i in 1:n) {
        next_vertex <- adj.list[[i]][2]  ## Second element is the next vertex
        prev_vertex <- adj.list[[i]][1]  ## First element is the previous vertex

        ## Calculate angle differences (always take the smaller angle)
        angle_to_next <- abs(angles[i] - angles[next_vertex])
        if (angle_to_next > pi) angle_to_next <- 2 * pi - angle_to_next

        angle_to_prev <- abs(angles[i] - angles[prev_vertex])
        if (angle_to_prev > pi) angle_to_prev <- 2 * pi - angle_to_prev

        weight.list[[i]] <- c(angle_to_prev, angle_to_next)
        names(weight.list[[i]]) <- as.character(adj.list[[i]])
    }
    names(weight.list) <- 1:n

    return(list(adj.list = adj.list, weight.list = weight.list))
}




#' Create a star graph by joining chain graphs
#'
#' This function creates a star graph by joining chain graphs of specified sizes to a central
#' vertex. The resulting star graph consists of a central vertex connected to multiple chains
#' of vertices.
#'
#' @param sizes A vector of positive integers specifying the sizes of the chain graphs to be
#'              joined to the central vertex. Each size represents the number of vertices in
#'              a chain graph.
#' @return A list representing the adjacency list of the created star graph.
#'
#' @examples
#' star_graph <- create.star.graph(c(3, 4, 2))
#' print(star_graph)
#'
#' @export
create.star.graph <- function(sizes) {

    if (!is.numeric(sizes) || any(sizes <= 0)) {
        stop("sizes must be a vector of positive integers.")
    }

    ## sizes <- as.integer(sizes)
    ## result <- .Call(S_create_star_graph, sizes)

    n.chains <- length(sizes)
    if (n.chains < 2) {
        stop("Size vector has to be of length at least three.")
    }

    star.graph <- create.chain.graph(sizes[1] + 1)$adj.list

    for (i in 2:n.chains) {
        chain.graph <- create.chain.graph(sizes[i] + 1)$adj.list
        star.graph <- join.graphs(star.graph, chain.graph, 1, 1)
    }

    return(star.graph)
}

#' Generate Synthetic Dataset on a Star-Shaped Geometric Space
#'
#' @description
#' Creates a synthetic dataset consisting of three main components:
#' 1. A geometric realization of a star graph embedded in a d-dimensional Euclidean space
#' 2. A smooth function defined on the vertices of this geometric realization
#' 3. Noisy observations obtained by adding random noise to the smooth function values
#'
#' The star graph consists of n.arms line segments (arms) meeting at a central point,
#' with each arm having potentially different lengths. Points are distributed along
#' these arms to create the graph vertices.
#'
#' @param n.points Integer. Total number of points to generate in the star graph, including
#'        the center point and arm endpoints. Must be greater than n.arms * min.n.pts.within.arm
#' @param n.arms Integer. Number of arms in the star graph. Must be at least 2
#' @param min.n.pts.within.arm Integer. Minimum number of interior points to place within each arm
#' @param min.arm.length Numeric. Minimum length for graph arms. Must be positive
#' @param max.arm.length Numeric. Maximum length for graph arms. Must be greater than min.arm.length
#' @param arm.lengths Numeric vector or NULL. If provided, must be of length n.arms with values
#'        between min.arm.length and max.arm.length. If NULL, lengths are randomly sampled
#'        from a uniform distribution
#' @param dim Integer. Dimension of the Euclidean space (>= 2) in which the star graph is embedded
#' @param fn Character. Type of smooth function to use:
#'        "exp" - Gaussian function exp(-((x - center)^2)/(2 * scale^2))
#'        "poly" - Cubic polynomial peaking at specified center
#'        "sin" - Damped sinusoidal function
#'        "wave" - Mexican hat wavelet
#' @param fn.center Controls the center point of the smooth function. Can be specified in three ways:
#'        1. NULL: Centers the function at the origin (graph center)
#'        2. Single number: Can be either:
#'           - An integer: Uses the coordinates of the vertex at that index
#'           - A non-integer: Places the center along the first arm at a fractional
#'             distance (e.g., 1.3 places it 30% along the first arm)
#'        3. Numeric vector of length dim: Uses these exact coordinates as the center
#'        Default is NULL, centering the function at the origin.
#' @param fn.scale Numeric. Scale parameter controlling the spread of the function.
#'        For "exp": standard deviation of the Gaussian
#'        For "poly": controls the rate of polynomial decay
#'        For "sin": controls the frequency and damping rate
#'        For "wave": controls the width of the wavelet
#' @param noise Character. Type of noise to add:
#'        "norm" - Normal distribution
#'        "laplace" - Laplace distribution
#'        "t" - Student's t-distribution with 3 degrees of freedom
#'        "unif" - Uniform distribution
#' @param noise.sd Numeric. Scale parameter for the noise distribution
#' @param rand.directions Logical. If TRUE, generates random directions for the arms.
#'        If FALSE and dim=2, uses evenly spaced points on a circle. Default is FALSE.
#' @param eps Numeric. Minimum relative distance between points within arms (0 < eps < 1)
#'
#' @return List containing:
#'   \item{points}{Matrix of vertex coordinates, where rows correspond to vertices}
#'   \item{adj.list}{List of adjacency lists defining the graph structure}
#'   \item{edge.lengths}{List of edge lengths corresponding to adj.list}
#'   \item{y.smooth}{Vector of smooth function values at the vertices}
#'   \item{y.noisy}{Vector of noisy observations (y.smooth + noise)}
#'   \item{arm.lengths}{Vector of actual arm lengths used}
#'   \item{center}{Coordinates of the star center (origin)}
#'   \item{endpoints}{Matrix of arm endpoints coordinates}
#'   \item{directions}{Matrix of unit vectors corresponding to the directions of the arms}
#'   \item{partition}{An integer vector showing the number of elements in each arm}
#'   \item{fn.center}{An index of vector of the center of the y.smooth}
#'
#' @examples
#' # Generate a 2D star graph with 5 arms and exponential function centered at origin
#' result <- generate.star.dataset(
#'   n.points = 20,
#'   n.arms = 5,
#'   min.n.pts.within.arm = 3,
#'   min.arm.length = 1,
#'   max.arm.length = 3,
#'   dim = 2,
#'   fn = "exp",
#'   fn.center = c(0, 0),
#'   fn.scale = 1,
#'   noise = "norm",
#'   noise.sd = 0.1
#' )
#'
#' @importFrom stats runif rnorm rexp rt dunif
#' @export
generate.star.dataset <- function(n.points,
                                  n.arms = 3,
                                  min.n.pts.within.arm = 3,
                                  min.arm.length = 0.5,
                                  max.arm.length = 2,
                                  arm.lengths = NULL,
                                  dim = 2,
                                  fn = "exp",
                                  fn.center = NULL,
                                  fn.scale = 1,
                                  noise = "norm",
                                  noise.sd = 0.1,
                                  rand.directions = FALSE,
                                  eps = 0.05) {

    ## Input validation
    if (dim < 2) stop("Dimension must be at least 2")
    if (n.points < n.arms * min.n.pts.within.arm) stop("n.points must be at least n.arms * min.n.pts.within.arm")
    if (min.arm.length <= 0) stop("min.arm.length must be positive")
    if (max.arm.length <= min.arm.length) stop("max.arm.length must be greater than min.arm.length")

    ## Generate or validate arm lengths
    arm.lengths.were.not.specified <- TRUE
    specified.arm.lengths <- arm.lengths
    if (is.null(arm.lengths)) {
        arm.lengths <- runif(n.arms, min = min.arm.length, max = max.arm.length)
    } else {
        if (length(arm.lengths) != n.arms) stop("arm.lengths must have length n.arms")
        arm.lengths.were.not.specified <- FALSE
    }

    ## Generate random directions on (dim-1)-sphere
    generate.sphere.points <- function(n, d) {
        points <- matrix(rnorm(n * d), n, d)
        points <- points / sqrt(rowSums(points^2))
        return(points)
    }

    ## Generate arm endpoints
    center <- rep(0, dim)

    if (rand.directions) {
        directions <- generate.sphere.points(n.arms, dim)
    } else if (!rand.directions && dim == 2) {
        directions <- as.matrix(generate.circle(n.arms, radius = 1))
    }

    endpoints <- matrix(nrow = n.arms, ncol = dim)
    for (j in seq(n.arms)) {
        endpoints[j,] <- directions[j,] * arm.lengths[j]
    }

    partition <- generate.partition(n.points, n.arms, min.n.pts.within.arm)

    ## Processing the first arm
    i <- 1
    lambda <- sort(runif(partition[i] - 2, min = eps, max = 1 - eps))

    ## distances to center
    dist.to.center <- c(arm.lengths[i] * lambda, arm.lengths[i])

    ## Creating adjacency and edge lengths lists
    graph.points <- matrix(nrow = n.points, ncol = dim)
    graph.points[1,] <- center

    adj.list <- list()
    edge.lengths <- list()
    adj.list[[1]] <- c(2) # 'center' neighbor
    edge.lengths[[1]] <- dist.to.center[1]
    for (j in seq(lambda)) {
        vertex.idx <- j + 1
        graph.points[vertex.idx,] <- directions[i,] * arm.lengths[i] * lambda[j]
        adj.list[[vertex.idx]] <- c(vertex.idx - 1, vertex.idx + 1)
        if (j == 1) {
            edge.lengths[[vertex.idx]] <- c(dist.to.center[j], dist.to.center[j + 1] - dist.to.center[j])
        } else {
            edge.lengths[[vertex.idx]] <- c(dist.to.center[j] - dist.to.center[j - 1], dist.to.center[j + 1] - dist.to.center[j])
        }
    }
    j <- length(lambda) + 1 # the end point of that arm
    vertex.idx <- j + 1
    graph.points[vertex.idx,] <- directions[i,] * arm.lengths[i]
    adj.list[[vertex.idx]] <- c(vertex.idx - 1)
    edge.lengths[[vertex.idx]] <- c(dist.to.center[j] - dist.to.center[j - 1])

    ## processing other arms
    for (i in 2:n.arms) {

        m <- length(adj.list) # the number of vertices constructed so far
        lambda <- sort(runif(partition[i] - 1, min = eps, max = 1 - eps))

        ## distances to center
        dist.to.center <- c(arm.lengths[i] * lambda, arm.lengths[i])

        adj.list[[1]] <- c(adj.list[[1]], m + 1) # 'center' neighbor
        edge.lengths[[1]] <- c(edge.lengths[[1]], dist.to.center[1])
        for (j in seq(lambda)) {
            vertex.idx <- j + m
            graph.points[vertex.idx,] <- directions[i,] * arm.lengths[i] * lambda[j]
            if (j == 1) {
                adj.list[[vertex.idx]] <- c(1, vertex.idx + 1)
                edge.lengths[[vertex.idx]] <- c(dist.to.center[j], dist.to.center[j + 1] - dist.to.center[j])
            } else {
                adj.list[[vertex.idx]] <- c(vertex.idx - 1, vertex.idx + 1)
                edge.lengths[[vertex.idx]] <- c(dist.to.center[j] - dist.to.center[j - 1], dist.to.center[j + 1] - dist.to.center[j])
            }
        }

        j <- length(lambda) + 1 # the end point of that arm
        vertex.idx <- j + m
        graph.points[vertex.idx,] <- directions[i,] * arm.lengths[i]
        adj.list[[vertex.idx]] <- c(vertex.idx - 1)
        edge.lengths[[vertex.idx]] <- c(dist.to.center[j] - dist.to.center[j - 1])
    }

    ## Determine function center
    ## Determine function center
    mu <- NA
    if (is.null(fn.center)) {
        mu <- rep(0, dim)  # Use graph center as default
    } else if (length(fn.center) == 1) {
        if (fn.center == 0) {
            mu <- rep(0, dim)  # Use graph center
        } else if (fn.center %% 1 == 0 && fn.center > 0) {
            mu <- graph.points[fn.center,]  # Use specified vertex
        } else if (fn.center %% 1 != 0) {
            lambda = fn.center %% 1
            mu <- lambda * endpoints[1,]
        }
    } else if (length(fn.center) == dim) {
        mu <- fn.center  # Use specified coordinates
    } else {
        stop("fn.center must be NULL, a vertex index, or a coordinate vector of length dim")
    }

    ## Generate smooth function values
    generate.smooth.values <- function(points, fn.type) {
        max.radius <- max(arm.lengths)

        ## Calculate distances for all function types
        dists.to.mu <- sqrt(rowSums(sweep(points, 2, mu, "-")^2))

        switch(fn.type,
               "exp" = {
                   exp(-((dists.to.mu)^2)/(2 * fn.scale^2))
               },
               "poly" = {
                   ## x is already the distance from the center point (dists.to.mu)
                   ## We use negative quadratic to create a peak (rather than a valley)
                   ## and scale it to have maximum value 1 at the center
                   -((dists.to.mu)/(max.radius * fn.scale))^2 + 1
               },
               "sin" = {
                   2 * sin(2 * pi * dists.to.mu/(max.radius * fn.scale)) *
                       exp(-dists.to.mu/(max.radius * fn.scale))
               },
               "wave" = {
                   r <- dists.to.mu/(max.radius * fn.scale)
                   (1 - r^2) * exp(-r^2/2)
               },
               stop("Unknown function type")
        )
    }

    ## Generate noise
    generate.noise <- function(n, noise.type, sd) {
        switch(noise.type,
               "norm" = rnorm(n, 0, sd),
               "laplace" = {
                   u <- runif(n, -0.5, 0.5)
                   -sd * sign(u) * log(1 - 2*abs(u))
               },
               "t" = rt(n, df=3) * sd,
               "unif" = runif(n, -sd, sd),
               stop("Unknown noise type")
               )
    }

    if (any(is.na(graph.points))) {
        print("graph.points")
        print(graph.points)
        stop("Some graph points are not initialized properly")
    }

    ## Generate final values
    y.smooth <- generate.smooth.values(graph.points, fn)
    noise.vals <- generate.noise(nrow(graph.points), noise, noise.sd)
    y.noisy <- y.smooth + noise.vals

    ## Return results
    result <- list(
        points = graph.points,
        adj.list = adj.list,
        edge.lengths = edge.lengths,
        y.smooth = y.smooth,
        y.noisy = y.noisy,
        arm.lengths = arm.lengths,
        center = center,
        endpoints = endpoints,
        directions = directions,
        partition = partition,
        fn.center = fn.center,
        specified.arm.lengths = specified.arm.lengths,
        mu = mu
    )

    class(result) <- "star_object"

    return(result)
}


#' Visualize 2D Star Graph with Smooth Function Values
#'
#' @description
#' Creates a visualization of a 2D star graph where vertices are colored
#' according to the values of the smooth function. The visualization includes
#' a reference circle, edges connecting vertices, and optional path highlighting.
#' Vertices can be labeled and the function center (mu) can be displayed.
#'
#' @param x A star_object (output from generate.star.dataset function)
#' @param y Optional numeric vector. If not NULL, it has have the same length as graph.data$y.smooth.
#' @param point.size Numeric. Size of vertex points (default = 1.5)
#' @param edge.col Character. Color of edges when not using color.edges (default = "gray70")
#' @param edge.lwd Numeric. Line width for edges (default = 1)
#' @param title Character. Plot title (default = "")
#' @param color.palette Character vector. Colors for the gradient (default = c("blue", "yellow", "red"))
#' @param color.edges Logical. If TRUE, edges are colored based on average function
#'        values of connected vertices (default = FALSE)
#' @param circle.color Character. Color of the reference circle (default = "gray")
#' @param pt.lab.adj Numeric vector of length 2. Text adjustment for vertex labels (default = c(1.0, 1.0))
#' @param pt.lab.cex Numeric. Character expansion factor for vertex labels (default = 0.75)
#' @param grid.pt.lab.adj Numeric vector of length 2. Text adjustment for grid point labels (default = c(1.0, 1.0))
#' @param grid.pt.lab.cex Numeric. Character expansion factor for grid point labels (default = 0.75)
#' @param grid.pt.color Character. Color for grid points (default = "cyan")
#' @param grid.pt.text.color Character. Color for grid point text labels (default = "cyan")
#' @param grid.pt.cex.factor Numeric. Size factor for grid points (default = 1.5)
#' @param show.vertex.labs Logical. Whether to show vertex labels (default = FALSE)
#' @param zero.based Logical. If TRUE, print vertex labels in 0-based format; otherwise 1-based (default = FALSE)
#' @param show.arm.labs Logical. Whether to show arm labels (default = TRUE)
#' @param arm.cex Numeric. Character expansion factor for arm labels (default = 2)
#' @param arms.adj Numeric vector of length 2. Text adjustment for arm labels (default = c(-0.5, 1))
#' @param arm.label.offset Numeric. It controls label distance from arm end.
#' @param show.mu Logical. Whether to show the function center point (default = TRUE)
#' @param mu.cex Numeric. Character expansion factor for function center point (default = 3)
#' @param epsilon Numeric. Margin factor for plot boundaries (default = 0.1)
#' @param legend.inset Numeric. Inset factor for legend position (default = 0.05)
#' @param legend.bty Character. Box type for legend ("o" for box, "n" for no box) (default = "n")
#' @param legend.title Character. Title for the color scale legend (default = "Function Value")
#' @param axes Logical. Whether to show plot axes (default = FALSE)
#' @param xlab Character. Label for x-axis (default = "")
#' @param ylab Character. Label for y-axis (default = "")
#' @param gpd.obj List. Optional path data for highlighting vertices. Each element should be
#'        a list containing $vertices (vector of vertex indices) and optionally $ref_vertex
#'        (index of reference vertex) (default = NULL)
#' @param gpd.path.index Integer. Index of the path to highlight from gpd.obj (default = NULL)
#' @param ugg.obj Object. Uniform grid graph object for additional visualization (default = NULL)
#' @param ... Additional graphical parameters (currently unused)
#'
#' @details
#' The function creates a visualization where vertices are arranged according to the star
#' graph structure. Vertices are colored based on their function values using the specified
#' color palette. A reference circle is drawn to show the maximum arm length. When gpd.obj
#' is provided, vertices in paths are highlighted with red circles and the reference vertex
#' (if specified) is highlighted with a larger blue circle.
#'
#' @return
#' NULL (generates a plot as a side effect)
#'
#' @examples
#' # Generate basic star graph data
#' data <- generate.star.dataset(
#'   n.points = 20,
#'   n.arms = 5,
#'   min.arm.length = 1,
#'   max.arm.length = 3,
#'   dim = 2,
#'   fn = "exp",
#'   noise = "norm",
#'   noise.sd = 0.1
#' )
#'
#' # Basic visualization
#' plot(data)
#'
#' # Visualization with custom styling and path highlighting
#' path_data <- list(
#'   list(vertices = c(1, 2, 3), ref_vertex = 1)
#' )
#' plot(data,
#'   point.size = 2,
#'   color.edges = TRUE,
#'   show.vertex.labs = TRUE,
#'   gpd.obj = path_data
#' )
#'
#' @importFrom graphics box text points lines
#' @importFrom grDevices colorRampPalette
#' @export
plot.star_object <- function(x,
                             y = NULL,
                             point.size = 1.5,
                             edge.col = "gray70",
                             edge.lwd = 1,
                             title = "",
                             color.palette = c("blue", "yellow", "red"),
                             color.edges = FALSE,
                             circle.color = "gray",
                             pt.lab.adj = c(1.0, 1.0),
                             pt.lab.cex = 0.75,
                             grid.pt.lab.adj = c(1.0, 1.0),
                             grid.pt.lab.cex = 0.75,
                             grid.pt.color = "cyan",
                             grid.pt.text.color = "cyan",
                             grid.pt.cex.factor = 1.5,
                             show.vertex.labs = FALSE,
                             zero.based = FALSE, # print vertex labels in 0-based or 1-based format
                             show.arm.labs = TRUE,
                             arm.cex = 2,
                             arms.adj = c(-0.5,1),
                             arm.label.offset = 0,
                             show.mu = TRUE,
                             mu.cex = 3,
                             epsilon = 0.1,
                             legend.inset = 0.05,
                             legend.bty = "n",
                             legend.title = "Function Value",
                             axes = FALSE,
                             xlab = "",
                             ylab = "",
                             gpd.obj = NULL,
                             gpd.path.index = NULL,
                             ugg.obj = NULL,
                             ...) {
    ## For S3 compatibility, rename x to graph.data internally
    graph.data <- x
    
    ## Extract required data
    points <- graph.data$points
    if (ncol(points) != 2) {
        stop("This visualization function only works for 2D star graphs (dim = 2)")
    }

    adj.list <- graph.data$adj.list

    if (is.null(y)) {
        y.smooth <- graph.data$y.smooth
    } else if(length(y) == length(graph.data$y.smooth)) {
        y.smooth <- y
    }

    ## Set up the plot area
    n <- 100
    radius <- 1
    if (!is.null(graph.data$specified.arm.lengths)) {
        radius <- max(graph.data$specified.arm.lengths)
    }
    circle.n100.df <- generate.circle(n, radius = radius)

    xlim <- range(c(as.numeric(points), as.numeric(as.matrix(circle.n100.df))))
    plot.margin <- max(c(abs(xlim[1]), xlim[2])) * (1 + epsilon)
    xlim <- c(-plot.margin, plot.margin)
    ylim <- xlim

    plot(circle.n100.df, type = 'l', xlab = xlab, ylab = ylab,
         las = 1, xlim = xlim, ylim = ylim,
         col = circle.color, axes = axes,
         main = title)

    if (!axes) {
        box()
    }

    if (show.arm.labs) {
        ## Calculate label positions with offset from the endpoints
        label_positions <- graph.data$endpoints + graph.data$directions * arm.label.offset

        text(label_positions, labels = seq(nrow(graph.data$directions)),
             cex = arm.cex, adj = arms.adj)

        ## text(radius * graph.data$directions, labels = seq(nrow(graph.data$directions)),
        ##      cex = arm.cex, adj = directions.adj)
    }

    if (show.mu) {
        points(graph.data$mu[1], graph.data$mu[2], cex = mu.cex)
    }

    # Create color scale
    n.colors <- 100
    color.ramp <- colorRampPalette(color.palette)
    colors <- color.ramp(n.colors)

    # Scale function values to [0,1] for color mapping
    y.scaled <- (y.smooth - min(y.smooth)) / (max(y.smooth) - min(y.smooth))
    point.colors <- colors[1 + floor(y.scaled * (n.colors - 1))]

    # Draw edges first (so they appear behind vertices)
    for (i in seq_along(adj.list)) {
        if (length(adj.list[[i]]) > 0) {
            for (j in adj.list[[i]]) {
                if (color.edges) {
                    # For colored edges, use average color of connected vertices
                    edge.color <- colors[1 + floor(mean(y.scaled[c(i,j)]) * (n.colors - 1))]
                } else {
                    edge.color <- edge.col
                }

                segments(points[i,1], points[i,2],
                        points[j,1], points[j,2],
                        col = edge.color, lwd = edge.lwd)
            }
        }
    }

    ## Draw vertices
    points(points, pch = 19, cex = point.size, col = point.colors)

    if (show.vertex.labs && !zero.based) {
        text(points, labels = seq(nrow(points)), cex = pt.lab.cex, adj = pt.lab.adj)
    } else if (show.vertex.labs && zero.based) {
        text(points, labels = seq(nrow(points)) - 1, cex = pt.lab.cex, adj = pt.lab.adj)
    }

    # Highlight paths from gpd.obj if provided
    if (!is.null(gpd.obj) && length(gpd.obj) > 0) {

        if (!is.null(gpd.path.index) && gpd.path.index > 0 && gpd.path.index <= length(gpd.obj) && gpd.path.index %% 1 == 0) {
            path <- gpd.obj[[gpd.path.index]]
            points(points[path$vertices, , drop = FALSE],
                   pch = 1,  # Circle
                   col = "red",
                   cex = 2 * point.size)
        } else {
            for (path in gpd.obj) {
                if (!is.null(path$vertices)) {
                    ## Draw red circles around vertices in the path
                    points(points[path$vertices, , drop = FALSE],
                           pch = 1,  # Circle
                           col = "red",
                           cex = 2 * point.size)
                }
            }
        }

        # Highlight reference vertex from first path with a larger blue circle
        path <- gpd.obj[[1]]
        if (!is.null(path$ref_vertex)) {
            points(points[path$ref_vertex, , drop = FALSE],
                   pch = 1,  # Circle
                   col = "blue",
                   cex = 3 * point.size,
                   lwd = 2)
        }
    }

    ## Helper function to find the second original vertex
    find_second_original_vertex <- function(start_vertex, neighbor, n.orig.vertices, adj_list, visited = NULL) {
        ## Initialize visited if this is the first call
        if (is.null(visited)) {
            visited <- c(start_vertex)
        }

        ## Get neighbors of current vertex
        curr_neighbors <- adj_list[[neighbor]]

        ## Check each neighbor
        for (next_neighbor in curr_neighbors) {
            ## Skip if we've already visited this vertex
            if (next_neighbor %in% visited) {
                next
            }

            ## If this neighbor is an original vertex, we found it
            if (next_neighbor <= n.orig.vertices) {
                return(list(
                    vertex = next_neighbor,
                    path = c(neighbor, next_neighbor)
                ))
            }

            ## Otherwise, recursively search through this neighbor's connections
            visited <- c(visited, next_neighbor)
            result <- find_second_original_vertex(start_vertex, next_neighbor,
                                                  n.orig.vertices, adj_list, visited)
            if (!is.null(result)) {
                return(list(
                    vertex = result$vertex,
                    path = c(neighbor, result$path)
                ))
            }
        }

        ## If we get here, no original vertex was found
        return(NULL)
    }

    find_two_original_vertex_neighbors <- function(start_vertex, n.orig.vertices, adj_list, weight_list) {
        ## Initialize queue for BFS and tracking structures
        queue <- list(list(
            vertex = start_vertex,
            path = numeric(0),
            total_weight = 0
        ))
        visited <- c(start_vertex)
        found_paths <- list()

        ## Continue BFS until we find two original vertices or exhaust all paths
        while (length(queue) > 0) {
            ## Get current vertex info from queue
            current <- queue[[1]]
            queue <- queue[-1]

            ## Get neighbors of current vertex
            neighbors <- adj_list[[current$vertex]]

            for (neighbor in neighbors) {
                ## Skip if we've visited this vertex
                if (neighbor %in% visited) {
                    next
                }

                ## Calculate weight to this neighbor
                weight_idx <- which(adj_list[[current$vertex]] == neighbor)
                new_weight <- current$total_weight + weight_list[[current$vertex]][weight_idx]

                ## Create new path including this neighbor
                new_path <- c(current$path, neighbor)

                if (neighbor <= n.orig.vertices) {
                    ## Found an original vertex - add to found_paths
                    found_paths[[length(found_paths) + 1]] <- list(
                        vertex = neighbor,
                        path = new_path,
                        total_weight = new_weight
                    )

                    ## If we found two original vertices, we can process them
                    if (length(found_paths) >= 2) {
                        ## Sort paths by total weight to get closest original vertices
                        found_paths <- found_paths[order(sapply(found_paths, function(x) x$total_weight))]

                        ## Return the two closest original vertices and their paths
                        return(list(
                            vertex1 = found_paths[[1]]$vertex,
                            path1 = found_paths[[1]]$path,
                            weight1 = found_paths[[1]]$total_weight,
                            vertex2 = found_paths[[2]]$vertex,
                            path2 = found_paths[[2]]$path,
                            weight2 = found_paths[[2]]$total_weight
                        ))
                    }
                } else {
                    ## Add unvisited grid vertex to queue
                    queue[[length(queue) + 1]] <- list(
                        vertex = neighbor,
                        path = new_path,
                        total_weight = new_weight
                    )
                    visited <- c(visited, neighbor)
                }
            }
        }

        ## Return NULL if we couldn't find two original vertices
        return(NULL)
    }

    ## Modified grid vertex plotting code
    if (!is.null(ugg.obj)) {

        n.orig.vertices <- length(adj.list)

        ## Get grid vertices from ugg.obj
        if ("ugg_grid_vertices" %in% names(ugg.obj)) {
            grid_vertices <- ugg.obj$ugg_grid_vertices
        } else {
            grid_vertices <- ugg.obj$grid_vertices
        }

        if ("ugg_adj_list" %in% names(ugg.obj)) {
            ugg.adj.list <- ugg.obj$ugg_adj_list
            ugg.weight.list <- ugg.obj$ugg_weight_list
        } else {
            ugg.adj.list <- ugg.obj$adj_list
            ugg.weight.list <- ugg.obj$weight_list
        }

        ## plotting grid vertices
        for (grid.vertex.idx in grid_vertices) {
            if (grid.vertex.idx <= n.orig.vertices) {
                points(points[grid.vertex.idx, ,drop = FALSE],
                       pch = 1,  # Circle
                       col = grid.pt.color,
                       cex = grid.pt.cex.factor * point.size,
                       lwd = 2)
                if (zero.based) {
                    text(points[grid.vertex.idx, ,drop = FALSE], labels = grid.vertex.idx - 1, cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                } else {
                    text(points[grid.vertex.idx, ,drop = FALSE], labels = grid.vertex.idx    , cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                }

            } else {
                neighbors <- ugg.adj.list[[grid.vertex.idx]]
                weights <- ugg.weight.list[[grid.vertex.idx]]
                ii <- which(neighbors <= n.orig.vertices)

                ## cat("\nProcessing vertex: ",grid.vertex.idx - 1,"\n")
                ## cat("neighbors: ")
                ## cat(neighbors - 1)
                ## cat("\n")
                ## cat("weights: ")
                ## cat(weights)
                ## cat("\n")
                ## cat("ii: ")
                ## cat(ii - 1)
                ## cat("\n")

                if (length(ii) == 2) {
                    ## Case 1: Both neighbors are original vertices
                    nhbr1 <- neighbors[ii[1]]
                    nhbr2 <- neighbors[ii[2]]
                    w1 <- weights[ii[1]]
                    w2 <- weights[ii[2]]
                    lambda <- w1 / (w1 + w2)
                    pt <- points[nhbr1, ,drop = FALSE] * (1 - lambda) + points[nhbr2, ,drop = FALSE] * lambda

                    points(pt,
                           pch = 1,  ## Circle
                           col = grid.pt.color,
                           cex = grid.pt.cex.factor * point.size,
                           lwd = 2)

                    if (zero.based) {
                        text(pt, labels = grid.vertex.idx - 1, cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                    } else {
                        text(pt, labels = grid.vertex.idx    , cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                    }
                } else if (length(ii) == 1) {
                    ## Case 2: Only one neighbor is an original vertex
                    orig_vertex <- neighbors[ii]
                    grid_neighbor <- neighbors[which(neighbors > n.orig.vertices)]
                    orig_weight <- weights[ii]

                    ## Find path to second original vertex
                    result <- find_second_original_vertex(grid.vertex.idx, grid_neighbor,
                                                          n.orig.vertices, ugg.adj.list)


                    ## cat("orig_vertex: ", orig_vertex - 1, "\n")
                    ## cat("grid_neighbor: ", grid_neighbor - 1, "\n")
                    ## cat("orig_weight: ", orig_weight, "\n")
                    ## cat("is.null(result): ", is.null(result), "\n")

                    if (!is.null(result)) {
                        second_orig_vertex <- result$vertex
                        path <- result$path

                        ## cat("second_orig_vertex: ", second_orig_vertex - 1, "\n")
                        ## cat("path: ")
                        ## cat(path - 1)
                        ## cat("\n")

                        ## Calculate total weight to second original vertex
                        total_weight <- orig_weight
                        prev_vertex <- grid.vertex.idx
                        for (i in 1:(length(path))) {
                            curr_vertex <- path[i]
                            weight_idx <- which(ugg.adj.list[[curr_vertex]] == prev_vertex)
                            total_weight <- total_weight + ugg.weight.list[[curr_vertex]][weight_idx]
                            prev_vertex <- path[i]
                        }

                        ## Calculate interpolation parameter
                        lambda <- orig_weight / total_weight

                        ## Calculate grid point position
                        pt <- points[orig_vertex,,drop = FALSE] * (1 - lambda) +
                            points[second_orig_vertex,,drop = FALSE] * lambda

                        points(pt,
                               pch = 1,  ## Circle
                               col = grid.pt.color,
                               cex = grid.pt.cex.factor * point.size,
                               lwd = 2)

                        if (zero.based) {
                            text(pt, labels = grid.vertex.idx - 1, cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                        } else {
                            text(pt, labels = grid.vertex.idx    , cex = grid.pt.lab.cex, adj = grid.pt.lab.adj, col = grid.pt.text.color)
                        }
                    } else {
                        cat("Could not find path to second original vertex for grid.vertex.idx:",
                            grid.vertex.idx, "\n")
                    }
                } else { ## all neighbors of grid.vertex.idx are grid vertices
                    result <- find_two_original_vertex_neighbors(grid.vertex.idx, n.orig.vertices, ugg.adj.list, ugg.weight.list)

                    if (!is.null(result)) {
                        ## Calculate interpolation parameter based on weights to both original vertices
                        total_weight <- result$weight1 + result$weight2
                        lambda <- result$weight1 / total_weight

                        ## Calculate grid point position
                        pt <- points[result$vertex1, ,drop = FALSE] * (1 - lambda) +
                            points[result$vertex2, ,drop = FALSE] * lambda

                        points(pt,
                               pch = 1,  #### Circle
                               col = grid.pt.color,
                               cex = grid.pt.cex.factor * point.size,
                               lwd = 2)

                        if (zero.based) {
                            text(pt, labels = grid.vertex.idx - 1,
                                 cex = grid.pt.lab.cex,
                                 adj = grid.pt.lab.adj,
                                 col = grid.pt.text.color)
                        } else {
                            text(pt, labels = grid.vertex.idx,
                                 cex = grid.pt.lab.cex,
                                 adj = grid.pt.lab.adj,
                                 col = grid.pt.text.color)
                        }

                        ## Add debug output
                        ## cat("Found two original vertices:\n")
                        ## cat("First vertex:", result$vertex1 - 1,
                        ##     "with total weight:", result$weight1, "\n")
                        ## cat("Path 1:", paste(result$path1 - 1, collapse = " -> "), "\n")
                        ## cat("Second vertex:", result$vertex2 - 1,
                        ##     "with total weight:", result$weight2, "\n")
                        ## cat("Path 2:", paste(result$path2 - 1, collapse = " -> "), "\n")
                    } else {
                        cat("Could not find two original vertices for grid.vertex.idx:",
                            grid.vertex.idx - 1, "\n")
                    }
                }
            }
        }
    }

    ## Add color legend
    legend.colors <- color.ramp(5)
    legend.labels <- round(seq(min(y.smooth), max(y.smooth), length.out = 5), 2)
    legend("topright",
           legend = legend.labels,
           pch = 19,
           col = legend.colors,
           title = legend.title,
           bty = legend.bty,
           inset = legend.inset)
}
