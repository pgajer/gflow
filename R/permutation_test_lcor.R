#' Permutation Test for Local Correlation Screening
#'
#' Performs feature-wise permutation testing for local correlation statistics
#' computed by \code{\link{lcor}}. By default, vertex labels of feature columns
#' are permuted (\code{permute = "z"}), which is appropriate for screening
#' associations between a fixed response field and candidate features.
#'
#' @param adj.list Adjacency list (1-based indices).
#' @param weight.list Edge-length (weight) list matching \code{adj.list}.
#' @param y Numeric response vector of length \code{n}.
#' @param z Numeric feature vector, matrix, or data frame with \code{n} rows.
#' @param type Local-correlation weighting mode passed to \code{\link{lcor}}.
#' @param y.diff.type Edge difference type for \code{y}; passed to \code{\link{lcor}}.
#' @param z.diff.type Edge difference type for \code{z}; passed to \code{\link{lcor}}.
#' @param statistic Summary statistic computed per feature from vertex-wise
#'   local-correlation coefficients.
#'   \describe{
#'     \item{"mean.abs"}{Mean absolute local correlation (default).}
#'     \item{"mean"}{Signed mean local correlation.}
#'     \item{"max.abs"}{Maximum absolute local correlation.}
#'   }
#' @param permute Permutation target:
#'   \describe{
#'     \item{"z"}{Permute feature rows (default).}
#'     \item{"y"}{Permute response values.}
#'   }
#' @param n.perm Integer \eqn{\ge 1}; number of permutations.
#' @param seed Optional integer RNG seed. If \code{NULL}, uses current RNG state.
#' @param return.perm.stats Logical; if \code{TRUE}, include permutation
#'   statistics matrix in the result.
#' @param verbose Logical; if \code{TRUE}, prints progress every 25 permutations.
#'
#' @return A list with class \code{"lcor_permutation_test"} containing:
#'   \describe{
#'     \item{table}{Data frame with feature, observed statistic, p-value, q-value.}
#'     \item{stat.obs}{Named numeric vector of observed statistics.}
#'     \item{p.value}{Named numeric vector of permutation p-values.}
#'     \item{q.value}{Named numeric vector of BH-adjusted p-values.}
#'     \item{statistic}{Statistic definition used.}
#'     \item{permute}{Permutation target used.}
#'     \item{n.perm}{Number of permutations.}
#'     \item{stat.perm}{Permutation statistics matrix (optional).}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- fit.rdgraph.regression(X, y, k = 10)
#' y.hat <- fit$fitted.values
#'
#' res <- permutation.test.lcor(
#'   adj.list = fit$graph$adj.list,
#'   weight.list = fit$graph$edge.length.list,
#'   y = y.hat,
#'   z = Z,
#'   type = "derivative",
#'   statistic = "mean.abs",
#'   permute = "z",
#'   n.perm = 200,
#'   seed = 1
#' )
#'
#' head(res$table[order(res$table$q.value), ])
#' }
#'
#' @export
permutation.test.lcor <- function(adj.list,
                                  weight.list,
                                  y,
                                  z,
                                  type = c("derivative", "unit", "sign"),
                                  y.diff.type = c("difference", "logratio"),
                                  z.diff.type = c("difference", "logratio"),
                                  statistic = c("mean.abs", "mean", "max.abs"),
                                  permute = c("z", "y"),
                                  n.perm = 200L,
                                  seed = 1L,
                                  return.perm.stats = FALSE,
                                  verbose = FALSE) {
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)
    statistic <- match.arg(statistic)
    permute <- match.arg(permute)

    if (!is.numeric(y)) {
        stop("y must be numeric.")
    }
    y <- as.double(y)
    n <- length(y)
    if (n < 2L) {
        stop("length(y) must be >= 2.")
    }

    if (!is.numeric(n.perm) || length(n.perm) != 1L || n.perm < 1L || n.perm != floor(n.perm)) {
        stop("n.perm must be a positive integer.")
    }
    n.perm <- as.integer(n.perm)

    if (!is.logical(return.perm.stats) || length(return.perm.stats) != 1L) {
        stop("return.perm.stats must be TRUE/FALSE.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("verbose must be TRUE/FALSE.")
    }

    if (is.null(dim(z))) {
        z <- matrix(as.double(z), ncol = 1L)
    } else if (is.data.frame(z) || is.matrix(z)) {
        z <- as.matrix(z)
        storage.mode(z) <- "double"
    } else {
        stop("z must be a numeric vector, matrix, or data.frame.")
    }

    if (nrow(z) != n) {
        stop("nrow(z) must equal length(y).")
    }

    if (is.null(colnames(z))) {
        colnames(z) <- paste0("feature", seq_len(ncol(z)))
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || seed != floor(seed)) {
            stop("seed must be an integer or NULL.")
        }
        set.seed(as.integer(seed))
    }

    lcor.to.matrix <- function(x, p) {
        if (is.list(x) && !is.null(x$column.coefficients)) {
            out <- as.matrix(x$column.coefficients)
        } else if (is.null(dim(x))) {
            out <- matrix(as.double(x), ncol = 1L)
        } else {
            out <- as.matrix(x)
        }
        if (ncol(out) != p) {
            stop("Unexpected lcor output shape during permutation testing.")
        }
        out
    }

    summarize.stat <- function(mat) {
        if (identical(statistic, "mean.abs")) {
            return(colMeans(abs(mat), na.rm = TRUE))
        }
        if (identical(statistic, "mean")) {
            return(colMeans(mat, na.rm = TRUE))
        }
        apply(abs(mat), 2L, max, na.rm = TRUE)
    }

    lcor.obs <- lcor(
        adj.list = adj.list,
        weight.list = weight.list,
        y = y,
        z = z,
        type = type,
        y.diff.type = y.diff.type,
        z.diff.type = z.diff.type
    )
    lcor.obs <- lcor.to.matrix(lcor.obs, p = ncol(z))
    stat.obs <- summarize.stat(lcor.obs)
    names(stat.obs) <- colnames(z)

    stat.perm <- matrix(NA_real_, nrow = n.perm, ncol = ncol(z))
    colnames(stat.perm) <- colnames(z)

    if (verbose) {
        message(sprintf("Running %d permutations (%s permutation)...", n.perm, permute))
    }

    for (b in seq_len(n.perm)) {
        if (identical(permute, "z")) {
            z.b <- z[sample.int(n), , drop = FALSE]
            y.b <- y
        } else {
            z.b <- z
            y.b <- y[sample.int(n)]
        }

        lcor.b <- lcor(
            adj.list = adj.list,
            weight.list = weight.list,
            y = y.b,
            z = z.b,
            type = type,
            y.diff.type = y.diff.type,
            z.diff.type = z.diff.type
        )
        lcor.b <- lcor.to.matrix(lcor.b, p = ncol(z))
        stat.perm[b, ] <- summarize.stat(lcor.b)

        if (verbose && (b %% 25L == 0L || b == n.perm)) {
            message(sprintf("  completed %d/%d permutations", b, n.perm))
        }
    }

    p.value <- (1 + colSums(t(t(stat.perm) >= stat.obs))) / (n.perm + 1L)
    q.value <- stats::p.adjust(p.value, method = "BH")

    out <- list(
        table = data.frame(
            feature = colnames(z),
            stat.obs = as.double(stat.obs),
            p.value = as.double(p.value),
            q.value = as.double(q.value),
            row.names = NULL
        ),
        stat.obs = stat.obs,
        p.value = p.value,
        q.value = q.value,
        statistic = statistic,
        permute = permute,
        n.perm = n.perm,
        call = match.call()
    )

    if (isTRUE(return.perm.stats)) {
        out$stat.perm <- stat.perm
    }

    class(out) <- c("lcor_permutation_test", "list")
    out
}
