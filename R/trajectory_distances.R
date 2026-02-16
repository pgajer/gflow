.make.dist.from.lower <- function(d.vec, n, labels = NULL, method = "DTW") {
  ## d.vec: numeric vector of length n*(n-1)/2 storing lower triangle distances
  ## n: number of objects
  ## returns: an object of class 'dist'

  n <- as.integer(n)
  if (n < 2L) stop("n must be >= 2.")
  if (!is.numeric(d.vec)) stop("d.vec must be numeric.")

  expected.len <- n * (n - 1L) / 2L
  if (length(d.vec) != expected.len) {
    stop("length(d.vec) must equal n*(n-1)/2: expected ", expected.len,
         ", got ", length(d.vec), ".")
  }
  if (anyNA(d.vec) || any(!is.finite(d.vec)) || any(d.vec < 0)) {
    stop("d.vec must be finite, nonnegative, and without NA.")
  }

  if (is.null(labels)) {
    labels <- as.character(seq_len(n))
  } else {
    if (length(labels) != n) stop("labels must have length n.")
    labels <- as.character(labels)
  }

  structure(d.vec,
            Size = n,
            Labels = labels,
            Diag = FALSE,
            Upper = FALSE,
            method = method,
            class = "dist")
}

#' DTW distance between two trajectories represented as feature sequences
#'
#' @param x.mat Numeric matrix (len_x x p) for trajectory x.
#' @param y.mat Numeric matrix (len_y x p) for trajectory y.
#' @param window.type DTW global constraint; see dtw::dtw() (e.g., "sakoechiba", "none").
#' @param window.size Optional window size (passed via ... to dtw::dtw()).
#' @param step.pattern DTW step pattern; see dtw::stepPattern().
#' @param normalize Logical; return normalized distance if available.
#' @param window.frac Fraction of trajectory length used to derive the default
#'   Sakoe-Chiba window width when \code{window.size} is NULL.
#' @return Numeric scalar distance.
#'
traj.dtw.distance <- function(x.mat,
                              y.mat,
                              window.type = "sakoechiba",
                              window.size = NULL,
                              step.pattern = dtw::symmetric2,
                              normalize = TRUE,
                              window.frac = 0.15) {

    ## local cost matrix (len_x x len_y)
    local.cost <- proxy::dist(x.mat, y.mat, method = "Euclidean")
    local.cost <- as.matrix(local.cost)

    ## Handle Sakoe-Chiba default window size robustly
    if (identical(window.type, "sakoechiba") && is.null(window.size)) {
        nx <- nrow(x.mat)
        ny <- nrow(y.mat)
        window.size <- max(abs(nx - ny), ceiling(window.frac * max(nx, ny)))
    }

    align <- dtw::dtw(local.cost,
                      step.pattern = step.pattern,
                      window.type = window.type,
                      window.size = window.size,
                      distance.only = TRUE)

    if (normalize && !is.null(align$normalizedDistance)) {
        return(align$normalizedDistance)
    }
    align$distance
}

#' Compute pairwise distances between trajectories
#'
#' @description
#' Computes a pairwise trajectory distance matrix (returned as a \code{dist} object)
#' for a list of trajectories. Distances are computed using one of:
#' \itemize{
#'   \item \code{"frechet"}: discrete Fr\'echet distance (order-respecting),
#'   \item \code{"lcss"}: LCSS-based distance (order-respecting, tolerance-based),
#'   \item \code{"hausdorff"}: Hausdorff distance (order-agnostic).
#' }
#' Distances are evaluated on either:
#' \itemize{
#'   \item a graph metric (shortest-path distances) provided via \code{pre}, or
#'   \item an embedding/feature matrix \code{E} defined on vertices (Euclidean metric).
#' }
#'
#' @param traj.list List of integer vectors; each element is a trajectory as vertex indices (1-based).
#' @param method Character string specifying the distance measure. One of
#'   \code{"lcss"}, \code{"frechet"}, \code{"hausdorff"}.
#' @param pre Optional list returned by \code{\link{graph.metric.precompute}}. If provided,
#'   the underlying point-to-point distance is the graph shortest-path metric.
#' @param E Optional numeric matrix of vertex features/coordinates. If provided, the underlying
#'   point-to-point distance is Euclidean in feature space. Must have \code{nrow(E)} equal
#'   to the number of vertices.
#' @param lcss.eps Optional nonnegative numeric scalar. LCSS match threshold (in the units
#'   of the underlying metric). If \code{NULL} and \code{method="lcss"} with \code{pre} provided,
#'   \code{lcss.eps} is selected automatically using \code{lcss.eps.auto.matchrate.graph()}.
#' @param lcss.match.rate Numeric scalar in \eqn{(0,1)} used for automatic LCSS threshold selection
#'   when \code{lcss.eps} is \code{NULL} and \code{pre} is provided.
#' @param lcss.delta Optional temporal constraint for LCSS. If finite, matches are only allowed
#'   when \code{abs(i - j) <= lcss.delta}. Use \code{Inf} to disable.
#' @param lcss.normalize Normalization of LCSS length used to form a distance in \eqn{[0,1]}.
#'   One of \code{"min"} or \code{"mean"}.
#' @param n.cores Integer scalar giving number of workers. For CRAN safety, parallel execution
#'   uses a PSOCK cluster via \code{parallel::makeCluster()} when \code{n.cores > 1}.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return A list with fields:
#' \itemize{
#'   \item \code{dist}: a \code{dist} object of pairwise trajectory distances.
#'   \item \code{params}: list of parameters used (including chosen \code{lcss.eps} if applicable).
#' }
#'
#' @examples
#' \dontrun{
#' ## Graph-metric example
#' pre <- graph.metric.precompute(g, traj.list)
#' res <- trajectory.dist(traj.list, method = "frechet", pre = pre)
#' hc <- hclust(res$dist, method = "ward.D2")
#' plot(hc)
#' }
#'
#' @export
trajectory.dist <- function(traj.list,
                            method = c("lcss", "frechet", "hausdorff"),
                            pre = NULL,
                            E = NULL,
                            lcss.eps = NULL,
                            lcss.match.rate = 0.10,
                            lcss.delta = Inf,
                            lcss.normalize = c("min", "mean"),
                            n.cores = 1L,
                            verbose = TRUE) {

  method <- match.arg(method)
  lcss.normalize <- match.arg(lcss.normalize)

  if (!is.list(traj.list) || length(traj.list) < 2L) {
    stop("traj.list must be a list of at least 2 trajectories.")
  }

  ## Exactly one of pre or E must be supplied
  if (is.null(pre) && is.null(E)) {
    stop("You must provide either 'pre' (graph metric) or 'E' (vertex feature matrix).")
  }
  if (!is.null(pre) && !is.null(E)) {
    stop("Provide only one of 'pre' or 'E', not both.")
  }

  ## Choose LCSS eps automatically if requested and feasible
  if (method == "lcss" && is.null(lcss.eps)) {
    if (!is.null(pre)) {
      lcss.eps <- lcss.eps.auto.matchrate.graph(pre,
                                                match.rate = lcss.match.rate,
                                                verbose = verbose)
    } else {
      stop("For method='lcss' with E provided, please supply lcss.eps explicitly (auto-selection currently implemented for graph metric via pre).")
    }
  }

  n <- length(traj.list)
  pairs <- utils::combn(n, 2L)
  npairs <- ncol(pairs)

  if (verbose) {
    message("Computing ", method, " distances for ", n, " trajectories (", npairs, " pairs).")
  }

  ## Worker that computes one pair distance
  pair.fun <- function(k) {
    i <- pairs[1L, k]
    j <- pairs[2L, k]

    if (!is.null(pre)) {
      C <- traj.local.cost.graph(traj.list[[i]], traj.list[[j]], pre)
    } else {
      ## Euclidean cost in feature space E
      Xi <- traj.points(traj.list[[i]], E)
      Xj <- traj.points(traj.list[[j]], E)
      C <- point.dist.matrix(Xi, Xj)
    }

    out <- switch(method,
      frechet = frechet.distance.discrete.cost(C),
      hausdorff = hausdorff.distance.cost(C),
      lcss = {
        L <- lcss.length.cost(C, eps = lcss.eps, delta = lcss.delta)
        if (lcss.normalize == "min") {
          s <- L / min(nrow(C), ncol(C))
        } else {
          s <- L / ((nrow(C) + ncol(C)) / 2)
        }
        1 - s
      }
    )

    out <- as.numeric(out)[1L]
    if (!is.finite(out) || out < 0) NA_real_ else out
  }

  ## Compute lower-triangle vector
  if (as.integer(n.cores) <= 1L) {
    d.vec <- vapply(seq_len(npairs), pair.fun, numeric(1))
  } else {
    n.cores <- as.integer(n.cores)
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    ## Export required objects/functions
    parallel::clusterExport(
      cl,
      varlist = c("traj.list", "pairs", "pre", "E", "method",
                  "lcss.eps", "lcss.delta", "lcss.normalize",
                  "traj.local.cost.graph", "traj.points", "point.dist.matrix",
                  "frechet.distance.discrete.cost", "lcss.length.cost",
                  "hausdorff.distance.cost", "pair.fun"),
      envir = environment()
    )

    d.list <- parallel::parLapply(cl, seq_len(npairs), pair.fun)
    d.vec <- unlist(d.list, use.names = FALSE)
  }

  expected.len <- n * (n - 1L) / 2L
  if (length(d.vec) != expected.len) {
    stop("Internal error: wrong d.vec length: expected ", expected.len, ", got ", length(d.vec), ".")
  }
  if (!is.numeric(d.vec) || anyNA(d.vec) || any(!is.finite(d.vec))) {
    stop("Distance computation produced non-finite values; inspect inputs and parameters.")
  }

  D <- .make.dist.from.lower(d.vec, n, method = method)

  if (verbose) {
    cat("typeof(D) =", typeof(D), "\n")
    cat("len(D)    =", length(D), "\n")
    cat("Size(D)   =", attr(D, "Size"), "\n")
  }

  list(dist = D,
       params = list(method = method,
                     metric = if (!is.null(pre)) "graph" else "features",
                     lcss.eps = lcss.eps,
                     lcss.match.rate = lcss.match.rate,
                     lcss.delta = lcss.delta,
                     lcss.normalize = lcss.normalize,
                     n.cores = as.integer(n.cores)))
}

#' Precompute graph shortest-path distances on trajectory vertices
#'
#' @description
#' Given an \code{igraph} object and a list of trajectories represented as vertex
#' index sequences, this function constructs the set of unique vertices appearing
#' in the trajectories and computes the all-pairs shortest-path distance matrix on
#' this vertex set. The resulting distance matrix can then be used to define a
#' graph-metric local cost matrix for trajectory-to-trajectory distances such as
#' DTW, discrete Fr\'echet, LCSS, and Hausdorff.
#'
#' @param g An \code{igraph} graph object.
#' @param traj.list List of integer vectors. Each element is a trajectory as a
#'   sequence of vertex indices (1-based, consistent with \code{igraph} vertex IDs).
#' @param weights Optional edge weights passed to \code{igraph::distances()}.
#'   May be \code{NULL} (unweighted), a numeric vector of edge weights, or the name
#'   of an edge attribute containing weights.
#' @param mode Character string specifying the direction of shortest paths, passed
#'   to \code{igraph::distances()}. One of \code{"all"}, \code{"out"}, \code{"in"}.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{vertex.set}: Integer vector of unique trajectory vertices (sorted).
#'   \item \code{index.map}: Integer vector mapping global vertex IDs to indices in
#'         \code{vertex.set} (with \code{NA} for vertices not in \code{vertex.set}).
#'   \item \code{D}: Numeric matrix of shortest-path distances between vertices in
#'         \code{vertex.set}.
#' }
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' traj.list <- list(c(1,2,3), c(1,10,9,8))
#' pre <- graph.metric.precompute(g, traj.list)
#' str(pre)
#' }
#'
#' @export
graph.metric.precompute <- function(g,
                                    traj.list,
                                    weights = NULL,
                                    mode = c("all", "out", "in"),
                                    verbose = TRUE) {
  mode <- match.arg(mode)

  v.set <- sort(unique(as.integer(unlist(traj.list, use.names = FALSE))))
  if (length(v.set) < 2L) stop("Too few vertices in traj.list.")

  if (verbose) message("Precomputing shortest-path distances for ", length(v.set), " vertices.")

  ## igraph::distances supports restricting to vertex subsets
  D <- igraph::distances(g,
                         v = v.set,
                         to = v.set,
                         weights = weights,
                         mode = mode)

  if (any(!is.finite(D))) {
    stop("Non-finite shortest-path distances detected (graph disconnected or invalid weights).")
  }

  index.map <- match(seq_len(max(v.set)), v.set)
  ## index.map[v] gives row/col index in D for vertex v, or NA if v not in v.set

  list(vertex.set = v.set, index.map = index.map, D = D)
}

#' Construct a graph-metric local cost matrix for a pair of trajectories
#'
#' @description
#' Builds the local cost matrix \eqn{C} for two trajectories \eqn{\gamma_i} and
#' \eqn{\gamma_j} using a precomputed graph shortest-path distance matrix.
#' If trajectories are \eqn{(v_1,\dots,v_n)} and \eqn{(w_1,\dots,w_m)}, the output
#' matrix is \eqn{C_{ab} = d_G(v_a, w_b)} where \eqn{d_G} is the graph shortest-path
#' metric.
#'
#' @param traj.i Integer vector of vertex indices (1-based) for the first trajectory.
#' @param traj.j Integer vector of vertex indices (1-based) for the second trajectory.
#' @param pre A list as returned by \code{\link{graph.metric.precompute}}.
#'
#' @return Numeric matrix of dimension \code{length(traj.i) x length(traj.j)}
#'   containing graph-metric distances between trajectory vertices.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' traj.list <- list(c(1,2,3), c(1,10,9,8))
#' pre <- graph.metric.precompute(g, traj.list)
#' C <- traj.local.cost.graph(traj.list[[1]], traj.list[[2]], pre)
#' dim(C)
#' }
#'
#' @export
traj.local.cost.graph <- function(traj.i, traj.j, pre) {
  ## pre: output of graph.metric.precompute()
  vi <- as.integer(traj.i)
  vj <- as.integer(traj.j)

  ii <- pre$index.map[vi]
  jj <- pre$index.map[vj]

  if (anyNA(ii) || anyNA(jj)) stop("Trajectory contains vertices not included in precomputed set.")

  pre$D[ii, jj, drop = FALSE]
}

#' Discrete Fr\'echet distance from a local cost matrix
#'
#' @description
#' Computes the discrete Fr\'echet distance between two discrete curves given their
#' local cost matrix \eqn{C}. The discrete Fr\'echet distance is an order-respecting
#' curve similarity measure that can be computed by dynamic programming.
#'
#' @details
#' Let \eqn{C} be an \eqn{n \times m} matrix with \eqn{C_{ij}} representing the
#' distance between the \eqn{i}-th point of the first curve and the \eqn{j}-th point
#' of the second curve under a chosen metric. The discrete Fr\'echet distance is the
#' value of a dynamic programming recurrence that aggregates costs using \code{max}
#' over matched pairs and \code{min} over predecessor states.
#'
#' @param C Numeric matrix of nonnegative local costs (distances) with dimensions
#'   \code{n x m}.
#'
#' @return Numeric scalar giving the discrete Fr\'echet distance.
#'
#' @examples
#' C <- matrix(c(0, 1, 2,
#'               1, 0, 1), nrow = 2, byrow = TRUE)
#' frechet.distance.discrete.cost(C)
#'
#' @export
frechet.distance.discrete.cost <- function(C) {

  n <- nrow(C)
  m <- ncol(C)
  ca <- matrix(NA_real_, n, m)

  ca[1, 1] <- C[1, 1]
  if (n >= 2L) for (i in 2:n) ca[i, 1] <- max(ca[i - 1L, 1L], C[i, 1L])
  if (m >= 2L) for (j in 2:m) ca[1, j] <- max(ca[1L, j - 1L], C[1L, j])

  if (n >= 2L && m >= 2L) {
    for (i in 2:n) {
      for (j in 2:m) {
        ca[i, j] <- max(C[i, j],
                        min(ca[i - 1L, j],
                            ca[i - 1L, j - 1L],
                            ca[i, j - 1L]))
      }
    }
  }
  ca[n, m]
}

#' Graph-metric discrete Fr\'echet distance between two trajectories
#'
#' @description
#' Computes the discrete Fr\'echet distance between two trajectories when the
#' underlying data graph is treated as a metric space. The point-to-point distance
#' between vertices is the shortest-path distance on the graph, precomputed by
#' \code{\link{graph.metric.precompute}}.
#'
#' @details
#' Let \eqn{\gamma_i = (v_1,\dots,v_n)} and \eqn{\gamma_j = (w_1,\dots,w_m)} be
#' vertex sequences. This function constructs the local cost matrix
#' \eqn{C_{ab} = d_G(v_a,w_b)} where \eqn{d_G} is the graph shortest-path metric,
#' then evaluates the discrete Fr\'echet distance via dynamic programming using
#' \code{\link{frechet.distance.discrete.cost}}.
#'
#' @param traj.i Integer vector of vertex indices (1-based) for the first trajectory.
#' @param traj.j Integer vector of vertex indices (1-based) for the second trajectory.
#' @param pre A list as returned by \code{\link{graph.metric.precompute}}.
#'
#' @return Numeric scalar giving the discrete Fr\'echet distance under the graph metric.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' traj.list <- list(c(1,2,3), c(1,10,9,8))
#' pre <- graph.metric.precompute(g, traj.list)
#' d <- traj.frechet.distance.graph(traj.list[[1]], traj.list[[2]], pre)
#' d
#' }
#'
#' @export
traj.frechet.distance.graph <- function(traj.i, traj.j, pre) {
  C <- traj.local.cost.graph(traj.i, traj.j, pre)
  frechet.distance.discrete.cost(C)
}

#' LCSS length from a local cost matrix
#'
#' @description
#' Computes the length of the Longest Common Subsequence (LCSS) alignment between
#' two trajectories given a local cost matrix \eqn{C}. A match occurs when
#' \eqn{C_{ij} \le \epsilon}. Optionally, matches can be restricted to a temporal
#' band \eqn{|i-j| \le \delta}.
#'
#' @details
#' LCSS is an order-respecting similarity measure that is robust to noise and
#' local mismatches. This function returns the raw LCSS length (an integer),
#' which can be converted to a normalized similarity or distance by the caller.
#'
#' @param C Numeric matrix of local costs (distances) with dimensions \code{n x m}.
#' @param eps Nonnegative numeric scalar. A pair \code{(i,j)} is considered a match
#'   if \code{C[i,j] <= eps}.
#' @param delta Optional temporal constraint. If finite, matches are only allowed
#'   when \code{abs(i - j) <= delta}. Use \code{Inf} to disable.
#'
#' @return Integer scalar giving the LCSS length.
#'
#' @examples
#' C <- matrix(c(0.0, 0.3, 0.9,
#'               0.2, 0.1, 0.4), nrow = 2, byrow = TRUE)
#' lcss.length.cost(C, eps = 0.25, delta = Inf)
#'
#' @export
lcss.length.cost <- function(C, eps, delta = Inf) {
  ## C: n x m local cost matrix
  ## eps: match threshold, match if C[i,j] <= eps
  ## delta: optional temporal band, match allowed only if |i-j| <= delta
  n <- nrow(C); m <- ncol(C)
  L <- matrix(0L, n + 1L, m + 1L)

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      if (is.finite(delta) && abs(i - j) > as.integer(delta)) {
        L[i + 1L, j + 1L] <- max(L[i, j + 1L], L[i + 1L, j])
      } else if (C[i, j] <= eps) {
        L[i + 1L, j + 1L] <- L[i, j] + 1L
      } else {
        L[i + 1L, j + 1L] <- max(L[i, j + 1L], L[i + 1L, j])
      }
    }
  }
  L[n + 1L, m + 1L]
}

#' Graph-metric LCSS distance between two trajectories
#'
#' @description
#' Computes an LCSS-based distance between two trajectories on a graph treated as a
#' metric space. The local cost matrix is constructed from precomputed graph
#' shortest-path distances, and LCSS matches are defined by a graph-distance
#' tolerance \code{eps}.
#'
#' @details
#' The function first forms the local cost matrix \eqn{C_{ab} = d_G(v_a, w_b)} where
#' \eqn{d_G} is the shortest-path distance on the graph and \eqn{(v_a)} and \eqn{(w_b)}
#' are trajectory vertex sequences. The LCSS length \eqn{L} is computed with match
#' threshold \code{eps} and optional temporal band \code{delta}. A normalized
#' similarity \eqn{s \in [0,1]} is then formed and returned as a distance
#' \eqn{1 - s}.
#'
#' @param traj.i Integer vector of vertex indices (1-based) for the first trajectory.
#' @param traj.j Integer vector of vertex indices (1-based) for the second trajectory.
#' @param pre A list as returned by \code{\link{graph.metric.precompute}}.
#' @param eps Nonnegative numeric scalar giving the graph-distance tolerance used
#'   to declare an LCSS match.
#' @param delta Optional temporal constraint. If finite, matches are only allowed
#'   when \code{abs(i - j) <= delta}. Use \code{Inf} to disable.
#' @param normalize Character string specifying normalization of LCSS length:
#'   \code{"min"} divides by \code{min(n,m)}; \code{"mean"} divides by \code{(n+m)/2}.
#'
#' @return Numeric scalar in \eqn{[0,1]} giving the LCSS distance.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' traj.list <- list(c(1,2,3), c(1,10,9,8))
#' pre <- graph.metric.precompute(g, traj.list)
#' d <- traj.lcss.distance.graph(traj.list[[1]], traj.list[[2]], pre,
#'                               eps = 2, delta = Inf, normalize = "min")
#' d
#' }
#'
#' @export
traj.lcss.distance.graph <- function(traj.i, traj.j, pre, eps, delta = Inf,
                                     normalize = c("min", "mean")) {
  normalize <- match.arg(normalize)
  C <- traj.local.cost.graph(traj.i, traj.j, pre)
  L <- lcss.length.cost(C, eps = eps, delta = delta)

  if (normalize == "min") {
    s <- L / min(nrow(C), ncol(C))
  } else {
    s <- L / ((nrow(C) + ncol(C)) / 2)
  }
  1 - s
}

#' Choose LCSS match tolerance on a graph metric by distance quantile
#'
#' @description
#' Selects the LCSS spatial match tolerance \code{eps} for trajectories on a graph
#' metric space as a quantile of shortest-path distances among vertices appearing
#' in the trajectories. This supports automatic, data-adaptive selection of the
#' LCSS matching threshold.
#'
#' @details
#' The LCSS matching rule declares a match between trajectory positions when the
#' graph distance \eqn{d_G(u,v)} is at most \code{eps}. This function chooses
#' \code{eps} as the \code{q}-quantile of sampled off-diagonal entries of the
#' shortest-path distance matrix \code{pre$D}. Distances are taken among vertices
#' that appear in the trajectory list used to build \code{pre} (i.e., the local
#' vertex set for a cell), which yields a cell-local and geometry-aware threshold.
#'
#' @param pre A list as returned by \code{\link{graph.metric.precompute}} containing
#'   a finite distance matrix \code{pre$D}.
#' @param q Numeric scalar in \eqn{(0,1)} giving the desired distance quantile used
#'   to define \code{eps}. Typical values are in \code{0.05}--\code{0.20}.
#' @param sample.pairs Integer scalar giving the number of off-diagonal distances
#'   to sample from \code{pre$D} for estimating the quantile. Sampling is used to
#'   reduce memory/time when the vertex set is large.
#' @param seed Optional integer seed for reproducible sampling.
#' @param verbose Logical; if \code{TRUE}, print messages describing the selection.
#'
#' @return Numeric scalar \code{eps} (nonnegative) to be used as the LCSS match threshold.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(50)
#' traj.list <- list(1:10, c(1,50:42))
#' pre <- graph.metric.precompute(g, traj.list)
#' eps <- lcss.eps.auto.graph(pre, q = 0.10, seed = 1)
#' eps
#' }
#'
#' @export
lcss.eps.auto.graph <- function(pre,
                                q = 0.10,
                                sample.pairs = 200000L,
                                seed = NULL,
                                verbose = TRUE) {

  if (is.null(pre$D) || !is.matrix(pre$D) || !is.numeric(pre$D)) {
    stop("pre must contain a numeric distance matrix pre$D.")
  }
  q <- as.numeric(q)
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0,1).")

  D <- pre$D
  n <- nrow(D)
  if (n < 2L) stop("pre$D must be at least 2 x 2.")

  ## Collect off-diagonal distances; sample if large
  if (!is.null(seed)) set.seed(as.integer(seed))

  ## Number of off-diagonal entries
  n.off <- n * (n - 1L)

  if (verbose) {
    message("Selecting LCSS eps from graph distances: n = ", n,
            ", off-diagonal entries = ", n.off,
            ", sampling up to ", sample.pairs, ".")
  }

  if (n.off <= sample.pairs) {
    vals <- D[upper.tri(D) | lower.tri(D)]
  } else {
    ## Sample indices uniformly from off-diagonal positions
    ## Sample row/col pairs and reject diagonal
    vals <- numeric(0)
    need <- as.integer(sample.pairs)
    while (length(vals) < need) {
      m <- need - length(vals)
      ii <- sample.int(n, size = m, replace = TRUE)
      jj <- sample.int(n, size = m, replace = TRUE)
      keep <- ii != jj
      if (any(keep)) {
        vals <- c(vals, D[cbind(ii[keep], jj[keep])])
      }
    }
    vals <- vals[seq_len(need)]
  }

  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) stop("No finite distances available to select eps.")

  eps <- as.numeric(stats::quantile(vals, probs = q, names = FALSE, type = 7))
  if (verbose) message("Selected eps = ", format(eps, digits = 6), " (q = ", q, ").")

  eps
}

#' Choose LCSS match tolerance on a graph metric by target match rate
#'
#' @description
#' Selects the LCSS spatial match tolerance \code{eps} for trajectories on a graph
#' metric space by targeting a desired match rate between random vertex pairs.
#' Specifically, \code{eps} is chosen as the \code{match.rate}-quantile of sampled
#' off-diagonal shortest-path distances among vertices appearing in the input
#' trajectories.
#'
#' @details
#' This function is a convenience wrapper around \code{\link{lcss.eps.auto.graph}}
#' that interprets the tuning parameter as a target match probability:
#' \eqn{\Pr(d_G(u,v) \le \epsilon) \approx \text{match.rate}} for random trajectory
#' vertices \eqn{u,v}. The distances \eqn{d_G} are taken from the precomputed matrix
#' \code{pre$D} returned by \code{\link{graph.metric.precompute}}.
#'
#' @param pre A list as returned by \code{\link{graph.metric.precompute}} containing
#'   a finite distance matrix \code{pre$D}.
#' @param match.rate Numeric scalar in \eqn{(0,1)} giving the target fraction of
#'   random vertex pairs that should be considered a match.
#' @param sample.pairs Integer scalar giving the number of off-diagonal distances
#'   to sample from \code{pre$D} for estimating the quantile. Sampling is used to
#'   reduce memory/time when the vertex set is large.
#' @param seed Optional integer seed for reproducible sampling.
#' @param verbose Logical; if \code{TRUE}, print messages describing the selection.
#'
#' @return Numeric scalar \code{eps} (nonnegative) to be used as the LCSS match threshold.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(50)
#' traj.list <- list(1:10, c(1,50:42))
#' pre <- graph.metric.precompute(g, traj.list)
#' eps <- lcss.eps.auto.matchrate.graph(pre, match.rate = 0.10, seed = 1)
#' eps
#' }
#'
#' @export
lcss.eps.auto.matchrate.graph <- function(pre,
                                          match.rate = 0.10,
                                          sample.pairs = 200000L,
                                          seed = NULL,
                                          verbose = TRUE) {
  ## Choose eps so that approximately match.rate of random vertex pairs satisfy d_G(u,v) <= eps
  match.rate <- as.numeric(match.rate)
  if (!is.finite(match.rate) || match.rate <= 0 || match.rate >= 1) {
    stop("match.rate must be in (0,1).")
  }
  lcss.eps.auto.graph(pre, q = match.rate, sample.pairs = sample.pairs,
                      seed = seed, verbose = verbose)
}

#' Hausdorff distance from a local cost matrix
#'
#' @description
#' Computes the (symmetric) Hausdorff distance between two discrete point sets
#' given their local cost matrix \eqn{C}. Here \eqn{C_{ij}} represents the distance
#' between the \eqn{i}-th point of the first set and the \eqn{j}-th point of the
#' second set under a chosen metric.
#'
#' @details
#' The directed Hausdorff distance from the first set to the second is
#' \deqn{h(X,Y) = \max_i \min_j C_{ij}.}
#' The (symmetric) Hausdorff distance is
#' \deqn{H(X,Y) = \max\{h(X,Y), h(Y,X)\}.}
#' This function implements the discrete version using \code{max} of row-wise minima
#' and column-wise minima of \code{C}.
#'
#' @param C Numeric matrix of nonnegative local costs (distances) with dimensions
#'   \code{n x m}.
#'
#' @return Numeric scalar giving the Hausdorff distance.
#'
#' @examples
#' C <- matrix(c(0, 1, 2,
#'               1, 0, 1), nrow = 2, byrow = TRUE)
#' hausdorff.distance.cost(C)
#'
#' @export
hausdorff.distance.cost <- function(C) {
    ## C: n x m local cost matrix
    hxy <- max(apply(C, 1L, min))
    hyx <- max(apply(C, 2L, min))
    max(hxy, hyx)
}

#' Extract trajectory point/feature matrix from a vertex-defined embedding
#'
#' @description
#' Given a trajectory represented as a sequence of vertex indices and a vertex
#' feature/coordinate matrix \code{E}, this function returns the corresponding
#' sequence of feature vectors \eqn{(E_{v_1}, \ldots, E_{v_k})} as a matrix.
#'
#' @param traj Integer vector of vertex indices (1-based).
#' @param E Numeric matrix of vertex features/coordinates with vertices in rows.
#'   Must satisfy \code{nrow(E) >= max(traj)}.
#'
#' @return Numeric matrix with \code{length(traj)} rows and \code{ncol(E)} columns.
#'
#' @examples
#' E <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
#' traj <- c(1, 3, 5, 2)
#' X <- traj.points(traj, E)
#' dim(X)
#'
#' @export
traj.points <- function(traj, E) {
  if (!is.matrix(E) || !is.numeric(E)) {
    stop("E must be a numeric matrix.")
  }
  v <- as.integer(traj)
  if (length(v) < 1L) {
    stop("traj must contain at least one vertex index.")
  }
  if (anyNA(v)) {
    stop("traj contains NA vertex indices.")
  }
  if (any(v < 1L) || any(v > nrow(E))) {
    stop("Trajectory contains vertex indices outside 1..nrow(E).")
  }
  E[v, , drop = FALSE]
}

#' Compute an Euclidean point-to-point distance matrix between two sequences
#'
#' @description
#' Computes the pairwise Euclidean distance matrix between the rows of \code{X}
#' and the rows of \code{Y}. This is commonly used as a local cost matrix for
#' trajectory alignment/distance methods (e.g., DTW, discrete Fr\'echet, LCSS, or
#' Hausdorff) when trajectories are embedded in a Euclidean feature space.
#'
#' @param X Numeric matrix of dimension \code{n x p}.
#' @param Y Numeric matrix of dimension \code{m x p}.
#'
#' @return Numeric matrix \code{C} of dimension \code{n x m} with entries
#'   \eqn{C_{ij} = \|X_i - Y_j\|_2}.
#'
#' @examples
#' X <- matrix(rnorm(5 * 2), nrow = 5)
#' Y <- matrix(rnorm(3 * 2), nrow = 3)
#' C <- point.dist.matrix(X, Y)
#' dim(C)
#'
#' @export
point.dist.matrix <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric.")
  }
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.")
  }
  if (nrow(X) < 1L || nrow(Y) < 1L) {
    stop("X and Y must each have at least one row.")
  }

  ## Squared Euclidean distances:
  ## ||x - y||^2 = ||x||^2 + ||y||^2 - 2 x^T y
  x2 <- rowSums(X * X)
  y2 <- rowSums(Y * Y)
  D2 <- outer(x2, y2, "+") - 2 * (X %*% t(Y))

  ## Numerical safety: clamp tiny negatives to 0 before sqrt
  sqrt(pmax(D2, 0))
}
