#' Compute Graph Spectrum
#'
#' Computes the first \code{nev} eigenvalues and eigenvectors of the
#' Laplacian matrix of a given graph. Supports both R-based and C-based
#' implementations for performance optimization.
#'
#' @param graph A list of integer vectors representing the adjacency list of
#'   the graph. Each element of the list corresponds to a vertex, and contains
#'   the indices of vertices adjacent to it.
#' @param nev An integer specifying the number of eigenvalues and eigenvectors
#'   to compute. If \code{NULL} (default), computes all eigenvalues except
#'   the trivial zero eigenvalue.
#' @param use.R Logical; if \code{TRUE}, uses the R-based implementation
#'   (slower but more portable). If \code{FALSE} (default), uses the C-based
#'   implementation for better performance.
#' @param return.Laplacian Logical; if \code{TRUE}, includes the Laplacian
#'   matrix in the return value. Default is \code{FALSE}.
#' @param return.dense Logical; if \code{TRUE} and \code{return.Laplacian} is
#'   \code{TRUE}, returns the Laplacian as a dense matrix. If \code{FALSE}
#'   (default), returns it as a sparse matrix. Only applicable when
#'   \code{use.R = FALSE}.
#'
#' @return A list containing:
#'   \item{evalues}{A numeric vector of eigenvalues in descending order.}
#'   \item{evectors}{A matrix of eigenvectors, where each column corresponds
#'     to an eigenvalue.}
#'   \item{laplacian}{(Optional) The Laplacian matrix of the graph, included
#'     only if \code{return.Laplacian = TRUE}.}
#'
#' @details
#' The Laplacian matrix L of a graph is defined as L = D - A, where D is the
#' degree matrix and A is the adjacency matrix. The eigenvalues and eigenvectors
#' of the Laplacian matrix provide important spectral properties of the graph.
#'
#' The function requires the \pkg{igraph} package when \code{use.R = TRUE}, and
#' the \pkg{Matrix} package when returning sparse Laplacian matrices.
#'
#' When \code{return.Laplacian = TRUE} on the C path and \pkg{Matrix} is not
#' installed, the function automatically returns a **dense** Laplacian (sets
#' \code{return.dense = TRUE}) to avoid a hard dependency on Matrix.
#'
#' @examples
#' # Create a simple graph adjacency list
#' graph <- list(c(2, 3), c(1, 3, 4), c(1, 2, 4), c(2, 3))
#'
#' # Compute spectrum using R implementation
#' spec_r <- graph.spectrum(graph, nev = 3, use.R = TRUE)
#'
#' # Compute spectrum using C implementation
#' spec_c <- graph.spectrum(graph, nev = 3, use.R = FALSE)
#'
#' # Get spectrum with Laplacian matrix
#' spec_lap <- graph.spectrum(graph, return.Laplacian = TRUE)
#'
#' @seealso \code{\link{graph.spectral.embedding}}, \code{\link{graph.low.pass.filter}}
#' @export
graph.spectrum <- function(graph,
                           nev = NULL,
                           use.R = FALSE,
                           return.Laplacian = FALSE,
                           return.dense = FALSE) {
  # Basic checks
  if (!is.list(graph)) stop("'graph' must be a list of integer vectors")
  n <- length(graph)
  if (n == 0L) stop("'graph' must contain at least one vertex")

  if (!is.null(nev)) {
    if (!is.numeric(nev) || length(nev) != 1L || nev < 1) stop("'nev' must be a positive integer")
    nev <- as.integer(nev)
    if (nev >= n) {
      warning("'nev' >= number of vertices; setting to n - 1")
      nev <- n - 1L
    }
  }

  stopifnot(is.logical(use.R), length(use.R) == 1L)
  stopifnot(is.logical(return.Laplacian), length(return.Laplacian) == 1L)
  stopifnot(is.logical(return.dense), length(return.dense) == 1L)

  if (use.R) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' is required when use.R = TRUE. Install it with install.packages('igraph').", call. = FALSE)
    }
    g.m <- convert.adjacency.list.to.adjacency.matrix(graph)
    g <- igraph::graph_from_adjacency_matrix(g.m, mode = "undirected")
    L  <- igraph::laplacian_matrix(g, normalized = FALSE)
    ed <- eigen(L)
    res <- list(evalues = ed$values, evectors = ed$vectors)
    if (return.Laplacian) res$laplacian <- L
    return(res)
  }

  # C implementation
  if (is.null(nev)) nev <- n - 1L

  # If user asked for sparse Laplacian but Matrix is missing, auto-switch to dense
  want_sparse <- return.Laplacian && !return.dense
  has_Matrix  <- requireNamespace("Matrix", quietly = TRUE)
  if (want_sparse && !has_Matrix) {
    warning("Matrix not installed: returning a dense Laplacian instead of sparse.")
    return.dense <- TRUE
  }

  # 0-based adjacency for C
  graph.0 <- lapply(graph, function(x) if (length(x)) as.integer(x - 1L) else integer(0))

  if (return.Laplacian) {
      ans <- .Call("S_graph_spectrum_plus",
                   graph.0,
                   as.integer(nev),
                   as.logical(return.dense))
    if (isTRUE(return.dense)) {
      ans$laplacian <- ans$dense_laplacian
      ans$dense_laplacian <- NULL
    } else {
      # Build sparse Laplacian using Matrix (guaranteed present here)
      ans$laplacian <- Matrix::sparseMatrix(
        i = ans$laplacian[[2]] + 1L,
        j = ans$laplacian[[3]] + 1L,
        x = ans$laplacian[[4]],
        dims = ans$laplacian[[1]],
        giveCsparse = TRUE
      )
    }
    return(list(evalues = ans$evalues, evectors = ans$evectors, laplacian = ans$laplacian))
  } else {
      ans <- .Call("S_graph_spectrum",
                   graph.0,
                   as.integer(nev))
    return(list(evalues = ans$evalues, evectors = ans$evectors))
  }
}

#' Compute Low-Pass Filter of a Function Over a Graph
#'
#' Computes the low-pass filter of a function over a graph using
#' the Graph Fourier Transform (GFT). The filter is applied by summing
#' the contributions of the eigenvectors starting from a specified index,
#' effectively filtering out high-frequency components.
#'
#' @param init.ev An integer specifying the index (1-based) of the first
#'   eigenvalue to include in the low-pass filter. Must be between 1 and
#'   the number of eigenvectors.
#' @param evectors A numeric matrix of eigenvectors of the graph Laplacian.
#'   Each column corresponds to an eigenvector, ordered by their associated
#'   eigenvalues.
#' @param y.gft A numeric matrix representing the Graph Fourier Transform (GFT)
#'   of the function over the graph. The first column should contain
#'   the GFT coefficients corresponding to each eigenvector.
#'
#' @return A numeric vector of length equal to the number of vertices,
#'   representing the filtered function values over the graph.
#'
#' @details
#' The Graph Fourier Transform decomposes a signal defined on the vertices
#' of a graph into its frequency components using the eigenvectors of the
#' graph Laplacian. Low-pass filtering retains only the low-frequency
#' components (associated with smaller eigenvalues), which typically
#' represent smooth variations over the graph structure.
#'
#' @examples
#' # Create example eigenvectors (2 vertices, 2 eigenvectors)
#' evectors <- matrix(c(1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)),
#'                   nrow = 2, ncol = 2)
#'
#' # Example GFT coefficients
#' y.gft <- matrix(c(3, 1), nrow = 2, ncol = 1)
#'
#' # Apply low-pass filter starting from the first eigenvector
#' filtered <- graph.low.pass.filter(1, evectors, y.gft)
#'
#' # Apply low-pass filter using only the second eigenvector
#' filtered_high <- graph.low.pass.filter(2, evectors, y.gft)
#'
#' @seealso \code{\link{graph.spectrum}}
#' @export
graph.low.pass.filter <- function(init.ev, evectors, y.gft) {

    # Input validation
    if (!is.numeric(init.ev) || length(init.ev) != 1 || init.ev < 1) {
        stop("'init.ev' must be a positive integer")
    }
    init.ev <- as.integer(init.ev)

    if (!is.matrix(evectors) && !is.numeric(evectors)) {
        stop("'evectors' must be a numeric matrix")
    }

    if (!is.matrix(y.gft) && !is.numeric(y.gft)) {
        stop("'y.gft' must be a numeric matrix")
    }

    # Ensure matrices
    if (!is.matrix(evectors)) {
        evectors <- as.matrix(evectors)
    }

    if (!is.matrix(y.gft)) {
        y.gft <- as.matrix(y.gft)
    }

    nev <- ncol(evectors)

    if (init.ev > nev) {
        stop("'init.ev' cannot exceed the number of eigenvectors")
    }

    if (nrow(y.gft) < nev) {
        stop("'y.gft' must have at least as many rows as there are eigenvectors")
    }

    if (ncol(y.gft) < 1) {
        stop("'y.gft' must have at least one column")
    }

    # Apply low-pass filter
    low.pass.y <- numeric(nrow(evectors))

    for (k in init.ev:nev) {
        low.pass.y <- low.pass.y + y.gft[k, 1] * as.numeric(evectors[, k])
    }

    return(low.pass.y)
}

#' Generate Spectral Embedding of a Graph
#'
#' Generates a spectral embedding of a graph into \eqn{R^d} using the
#' eigenvectors corresponding to the smallest non-zero eigenvalues of
#' the graph's Laplacian matrix. This embedding preserves the graph's
#' structural properties in a low-dimensional space.
#'
#' @param evectors A numeric matrix of eigenvectors of the graph Laplacian.
#'   Each column corresponds to an eigenvector, ordered by their corresponding
#'   eigenvalues in descending order (largest to smallest).
#' @param dim An integer specifying the dimension of the embedding space
#'   (\eqn{R^{dim}}). Must be between 1 and \code{ncol(evectors) - 1}.
#' @param evalues An optional numeric vector of eigenvalues corresponding to
#'   the eigenvectors, also in descending order. If provided, the eigenvectors
#'   will be scaled by dividing by the square root of their corresponding
#'   eigenvalues, which can improve the embedding quality.
#'
#' @return A numeric matrix of dimension \code{nrow(evectors)} by \code{dim},
#'   where each row represents the spectral embedding coordinates of a vertex
#'   in \eqn{R^{dim}}.
#'
#' @details
#' Spectral embedding is a dimensionality reduction technique that uses the
#' eigenvectors of the graph Laplacian to embed vertices in a low-dimensional
#' space. The embedding preserves the connectivity structure of the graph,
#' placing connected vertices close to each other in the embedded space.
#'
#' The function excludes the trivial constant eigenvector (associated with
#' eigenvalue 0) and uses the next \code{dim} eigenvectors corresponding to
#' the smallest non-zero eigenvalues.
#'
#' When eigenvalues are provided, the eigenvectors are normalized by dividing
#' by the square root of their eigenvalues, which is a common practice in
#' spectral clustering and can lead to better separation of clusters.
#'
#' @examples
#' # Create example eigenvectors for 5 vertices
#' set.seed(123)
#' evectors <- matrix(rnorm(5 * 5), nrow = 5, ncol = 5)
#'
#' # Apply Gram-Schmidt to ensure orthogonality
#' evectors <- qr.Q(qr(evectors))
#'
#' # Generate 2D embedding without eigenvalue scaling
#' embedding_2d <- graph.spectral.embedding(evectors, dim = 2)
#'
#' # Generate 3D embedding with eigenvalue scaling
#' evalues <- c(5, 3, 2, 0.5, 0)  # Example eigenvalues in descending order
#' embedding_3d_scaled <- graph.spectral.embedding(evectors, dim = 3, evalues)
#'
#' @references
#' von Luxburg, U. (2007). A tutorial on spectral clustering.
#' \emph{Statistics and Computing}, 17(4), 395-416.
#'
#' @seealso \code{\link{graph.spectrum}}
#' @export
graph.spectral.embedding <- function(evectors, dim, evalues = NULL) {

    # Input validation
    if (!is.matrix(evectors) && !is.numeric(evectors)) {
        stop("'evectors' must be a numeric matrix")
    }

    if (!is.matrix(evectors)) {
        evectors <- as.matrix(evectors)
    }

    if (!is.numeric(dim) || length(dim) != 1 || dim < 1) {
        stop("'dim' must be a positive integer")
    }
    dim <- as.integer(dim)

    if (dim >= ncol(evectors)) {
        stop("'dim' must be less than the number of eigenvectors")
    }

    # Validate evalues if provided
    if (!is.null(evalues)) {
        if (!is.numeric(evalues)) {
            stop("'evalues' must be a numeric vector")
        }

        if (length(evalues) != ncol(evectors)) {
            stop("The length of 'evalues' must match the number of columns in 'evectors'")
        }

        # Check for non-positive eigenvalues in the selected range
        selected_indices <- (ncol(evectors) - dim):(ncol(evectors) - 1)
        selected_evalues <- evalues[selected_indices]

        if (any(selected_evalues <= 0)) {
            warning("Some selected eigenvalues are non-positive. ",
                   "This may lead to numerical issues.")
        }
    }

    # Select the dim eigenvectors, excluding the last one (constant eigenvector)
    # The eigenvectors are ordered by descending eigenvalues, so we want the
    # ones just before the last (which corresponds to eigenvalue 0)
    embedding <- evectors[, (ncol(evectors) - dim):(ncol(evectors) - 1), drop = FALSE]

    # Scale by eigenvalues if provided
    if (!is.null(evalues)) {
        # Normalize each eigenvector by the square root of its eigenvalue
        for (i in 1:dim) {
            eval_idx <- ncol(evectors) - dim - 1 + i
            if (evalues[eval_idx] > 0) {
                embedding[, i] <- embedding[, i] / sqrt(evalues[eval_idx])
            }
        }
    }

    # Ensure column names are meaningful
    colnames(embedding) <- paste0("Dim", 1:dim)

    return(embedding)
}

