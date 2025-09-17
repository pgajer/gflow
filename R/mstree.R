#' C implementation of a Minimal Spanning Tree on a matrix X with > 1 column
#'
#' @param X A numeric matrix or data frame with more than one column.
#' @param i An index of the row of X from which the algorithm starts. Default is 1.
#' @return A list containing two elements:
#'   \item{edges}{A matrix of edge indices, where each row represents an edge in the MST.}
#'   \item{edge.lens}{A vector of edge lengths corresponding to the edges in the MST.}
#' @importFrom FNN get.knn
#' @export
mstree <- function(X, i = 1) {

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

    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    if (!is.numeric(i) || length(i) != 1 || i != round(i) || i < 1 || i > nrow(X)) {
        stop("i must be a positive integer not greater than the number of points in X")
    }

    n <- nrow(X)
    nn <- get.knn(X, k = n-1)
    nn.i <- nn$nn.index
    nn.d <- nn$nn.dist

    ldist <- max(nn.d[,n - 1]) + 1 ## a distance larger than all distances
    edges <- numeric((n - 1) * 2)
    edge.lens <- numeric(n - 1)

    out <- .C(C_mstree,
             as.integer(i - 1),
             as.integer(t(nn.i - 1)),
             as.double(t(nn.d)),
             as.double(ldist),
             as.integer(n),
             edges=as.integer(edges),
             edge.lens=as.double(edge.lens))

    edges <- matrix(out$edges, ncol=2, byrow=TRUE) + 1

    list(edges=edges,
         edge.lens=out$edge.lens)
}
