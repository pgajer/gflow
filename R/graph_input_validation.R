.validate.numeric.data.matrix <- function(X) {
    if (!(is.matrix(X) || is.data.frame(X))) {
        stop("'X' must be a matrix or data frame.", call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X)) {
        stop("'X' must contain numeric data.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("'X' cannot contain NA, NaN, or Inf values.", call. = FALSE)
    }
    if (nrow(X) < 2L) {
        stop("'X' must contain at least two observations.", call. = FALSE)
    }
    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }
    X
}
