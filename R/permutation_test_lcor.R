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
#' @param seed Optional integer RNG seed. A supplied seed is local: the caller's
#'   random-number state is restored, even on error. `NULL` advances the current stream.
#' @param return.perm.stats Logical; if \code{TRUE}, include permutation
#'   statistics matrix in the result.
#' @param verbose Logical; if \code{TRUE}, prints progress every 25 permutations.
#'
#' @param epsilon,winsorize.quantile,hop.radius Settings passed unchanged to
#'   [lcor()] for both observed and permuted statistics.
#' @param strata Optional vector of group labels, one per vertex, without missing
#'   values. Permutations stay within groups; singleton groups remain fixed.
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
#' @details
#' Permuted rows must be exchangeable under the null conditional on the fixed
#' graph, edge lengths, and unpermuted field. A common permutation is applied
#' across feature columns. Local correlations and their feature summaries are
#' recomputed, but graphs, fitted fields, and neighborhood selection are not.
#' With `strata`, rows must be exchangeable within the supplied groups under
#' that same conditional null. Grouping alone does not make spatial dependence
#' exchangeable, and is not a replacement for an appropriate null model. The signed `mean` statistic uses an upper-tail comparison.
#' P-values use the add-one Monte Carlo correction. BH-adjusted `q.value`
#' refers to feature tests, not vertex tests, and its false discovery rate
#' interpretation requires valid p-values and suitable dependence conditions.
#'
#' @examples
#' # Independent null features on a fixed path: rows are exchangeable.
#' local({
#'   had <- exists(".Random.seed", envir = .GlobalEnv)
#'   if (had) old <- .Random.seed
#'   on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else
#'     rm(".Random.seed", envir = .GlobalEnv))
#'   set.seed(42)
#'   a <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
#'   w <- lapply(a, function(v) rep(1, length(v)))
#'   z <- matrix(rnorm(15), 5, 3, dimnames = list(NULL, c("A", "B", "C")))
#'   result <- permutation.test.lcor(a, w, c(0, 1, 2, 1, 0), z,
#'                                    n.perm = 19, seed = 7)
#'   print(result)
#' })
#' # Nineteen draws give resolution 1/20; this tiny run illustrates the null,
#' # not a well-powered study. BH values concern features, not vertices.
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
                                  verbose = FALSE,
                                  epsilon = 0,
                                  winsorize.quantile = 0,
                                  hop.radius = 1L,
                                  strata = NULL) {
    type <- match.arg(type)
    y.diff.type <- match.arg(y.diff.type)
    z.diff.type <- match.arg(z.diff.type)
    statistic <- match.arg(statistic)
    permute <- match.arg(permute)

    if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y))) {
        stop("y must be a finite numeric vector.")
    }
    y <- as.double(y)
    n <- length(y)
    if (n < 2L) {
        stop("length(y) must be >= 2.")
    }

    if (!is.numeric(n.perm) || length(n.perm) != 1L || !is.finite(n.perm) || n.perm > .Machine$integer.max || n.perm < 1L || n.perm != floor(n.perm)) {
        stop("n.perm must be a positive integer.")
    }
    n.perm <- as.integer(n.perm)

    if (!is.logical(return.perm.stats) || length(return.perm.stats) != 1L || is.na(return.perm.stats)) {
        stop("return.perm.stats must be TRUE/FALSE.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE/FALSE.")
    }

    if (!is.numeric(z) && !(is.data.frame(z) && all(vapply(z, is.numeric, logical(1))))) {
        stop("z must contain only numeric values.")
    }
    if (is.null(dim(z))) {
        z <- matrix(as.double(z), ncol = 1L)
    } else if (is.data.frame(z) || is.matrix(z)) {
        z <- as.matrix(z)
        storage.mode(z) <- "double"
    } else {
        stop("z must be a numeric vector, matrix, or data.frame.")
    }

    if (any(!is.finite(z)) || ncol(z) < 1L) stop("z must have finite values and at least one feature.")
    if (nrow(z) != n) {
        stop("nrow(z) must equal length(y).")
    }

    if (is.null(colnames(z))) {
        colnames(z) <- paste0("feature", seq_len(ncol(z)))
    }

    if (!is.null(strata) && (length(strata) != n || !is.atomic(strata) ||
                             !is.null(dim(strata)) || anyNA(strata))) {
        stop("strata must have one nonmissing group label per vertex.")
    }
    groups <- if (is.null(strata)) list(seq_len(n)) else split(seq_len(n), as.character(strata))
    permutation <- function() {
        if (is.null(strata)) return(sample.int(n))
        index <- seq_len(n)
        for (g in groups) index[g] <- g[sample.int(length(g))]
        index
    }
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || abs(seed) > .Machine$integer.max || seed != floor(seed)) {
            stop("seed must be an integer or NULL.")
        }
        had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (had.seed) assign(".Random.seed", old.seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        }, add = TRUE)
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
        z.diff.type = z.diff.type,
        epsilon = epsilon, winsorize.quantile = winsorize.quantile,
        hop.radius = hop.radius
    )
    lcor.obs <- lcor.to.matrix(lcor.obs, p = ncol(z))
    stat.obs <- summarize.stat(lcor.obs)
    names(stat.obs) <- colnames(z)

    stat.perm <- if (return.perm.stats) matrix(NA_real_, nrow = n.perm, ncol = ncol(z),
                                             dimnames = list(NULL, colnames(z))) else NULL
    exceed <- numeric(ncol(z))

    if (verbose) {
        message(sprintf("Running %d permutations (%s permutation)...", n.perm, permute))
    }

    for (b in seq_len(n.perm)) {
        if (identical(permute, "z")) {
            z.b <- z[permutation(), , drop = FALSE]
            y.b <- y
        } else {
            z.b <- z
            y.b <- y[permutation()]
        }

        lcor.b <- lcor(
            adj.list = adj.list,
            weight.list = weight.list,
            y = y.b,
            z = z.b,
            type = type,
            y.diff.type = y.diff.type,
            z.diff.type = z.diff.type,
        epsilon = epsilon, winsorize.quantile = winsorize.quantile,
        hop.radius = hop.radius
        )
        lcor.b <- lcor.to.matrix(lcor.b, p = ncol(z))
        simulated <- summarize.stat(lcor.b)
        exceed <- exceed + (simulated >= stat.obs)
        if (return.perm.stats) stat.perm[b, ] <- simulated

        if (verbose && (b %% 25L == 0L || b == n.perm)) {
            message(sprintf("  completed %d/%d permutations", b, n.perm))
        }
    }

    p.value <- setNames((1 + exceed) / (n.perm + 1), colnames(z))
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
        tail = "upper",
        settings = list(type = type, y.diff.type = y.diff.type, z.diff.type = z.diff.type,
                        epsilon = epsilon, winsorize.quantile = winsorize.quantile,
                        hop.radius = hop.radius, seed = seed),
        strata = strata,
        group.sizes = lengths(groups),
        call = match.call()
    )

    if (isTRUE(return.perm.stats)) {
        out$stat.perm <- stat.perm
    }

    class(out) <- c("lcor_permutation_test", "list")
    out
}


#' Print a Local-Correlation Permutation Test
#'
#' @param x A result from [permutation.test.lcor()].
#' @param ... Unused.
#' @param n Maximum number of feature rows to print.
#' @return `x`, invisibly. Complete results remain in `x$table`.
#' @export
print.lcor_permutation_test <- function(x, ..., n = 5L) {
    .validate.lcor.permutation.test(x)
    n <- .basin.summary.top.k(n, "n")
    cat("Local-correlation permutation test\n")
    cat("  Statistic: ", x$statistic, "; upper tail; target: ", x$permute, "\n", sep = "")
    cat("  Fixed: graph, edge lengths, unpermuted field, statistic settings.\n")
    cat("  Null: exchangeable ", x$permute, " rows",
        if (is.null(x$strata)) ".\n" else " within supplied groups.\n", sep = "")
    cat("  Monte Carlo permutations: ", x$n.perm, "; minimum p-value: ",
        format(1 / (x$n.perm + 1)), "; tested features: ", nrow(x$table), "\n", sep = "")
    if (!is.null(x$settings)) cat("  Weighting: ", x$settings$type, "; hop radius: ",
                                x$settings$hop.radius, "\n", sep = "")
    tab <- x$table[order(x$table$q.value, x$table$p.value), , drop = FALSE]
    print(utils::head(tab, n), row.names = FALSE)
    cat("BH q-values apply to features and require valid null p-values.\n")
    invisible(x)
}

.validate.lcor.permutation.test <- function(x) {
    .validate.s3.contract(x, "lcor_permutation_test",
        fields = c("table", "stat.obs", "p.value", "q.value", "statistic", "permute", "n.perm"),
        storage = "list")
    if (!is.data.frame(x$table) ||
        !all(c("feature", "stat.obs", "p.value", "q.value") %in% names(x$table)))
        stop("Invalid lcor permutation result table.")
    invisible(x)
}
