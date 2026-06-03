#' Kernel Local Polynomial Regression with CV-Selected Neighborhoods
#'
#' Fits a kernel local polynomial smoother and selects its support size,
#' polynomial degree, and kernel by cross-validation.  By default, the smoother
#' works in the observed ambient coordinates: each prediction point uses its
#' nearest training points in Euclidean distance, centers the support at the
#' prediction point, fits a weighted local polynomial, and uses the fitted
#' intercept as the prediction.
#'
#' The optional \code{coordinate.method = "local.pca"} mode keeps the same
#' support and kernel weighting rule, but builds the local polynomial in a local
#' PCA chart centered at each prediction point.  With \code{chart.dim = "auto"},
#' the chart dimension is estimated from observed \code{X} only, using the same
#' shared local-PCA dimension helper used by LPL-TF and S-LPL-TF.
#'
#' @param X Numeric design/coordinate matrix with one observation per row.
#' @param y Numeric response vector with length \code{nrow(X)}.
#' @param foldid Optional positive integer vector assigning rows to CV folds.
#' @param support.grid Integer candidate neighborhood sizes.
#' @param degree.grid Integer polynomial degrees. Currently degrees 0, 1, and 2
#'   are supported.
#' @param kernel.grid Candidate kernels. Supported kernels are
#'   \code{"gaussian"}, \code{"tricube"}, \code{"epanechnikov"}, and
#'   \code{"triangular"}.
#' @param cv.folds Number of folds used when \code{foldid} is not supplied.
#' @param cv.seed Random seed used to generate folds when \code{foldid} is not
#'   supplied.
#' @param X.eval Optional matrix of prediction locations. Defaults to \code{X}.
#' @param coordinate.method Local coordinate system. \code{"coordinates"} uses
#'   centered ambient coordinates. \code{"local.pca"} uses a local PCA chart.
#' @param chart.dim Chart dimension for \code{coordinate.method = "local.pca"}.
#'   If \code{NULL}, defaults to \code{ncol(X)}. The special value
#'   \code{"auto"} estimates chart dimension from observed \code{X} only.
#' @param auto.chart.support.metric Support system used when
#'   \code{chart.dim = "auto"}. Included for consistency with LPL-TF and
#'   S-LPL-TF; because this smoother uses coordinate supports,
#'   \code{"operator"} is equivalent to \code{"coordinates"}.
#' @param auto.chart.selection.metric Which auto chart-dimension diagnostic to
#'   select when both diagnostics are requested.
#'
#' @return A list of class \code{"kernel.local.polynomial.cv"} with fitted
#'   values, selected parameters, and a candidate CV table.
#' @export
kernel.local.polynomial.cv <- function(
    X, y, foldid = NULL,
    support.grid = c(10L, 15L, 20L),
    degree.grid = 0:2,
    kernel.grid = c("gaussian", "tricube"),
    cv.folds = 5L,
    cv.seed = 1L,
    X.eval = NULL,
    coordinate.method = c("coordinates", "local.pca"),
    chart.dim = NULL,
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator")) {

    X <- as.matrix(X)
    y <- as.numeric(y)
    if (!is.numeric(X) || !length(X) || any(!is.finite(X))) {
        stop("'X' must be a finite numeric matrix.", call. = FALSE)
    }
    if (length(y) != nrow(X) || any(!is.finite(y))) {
        stop("'y' must be a finite numeric vector with length nrow(X).",
             call. = FALSE)
    }
    X.eval <- if (is.null(X.eval)) X else as.matrix(X.eval)
    if (ncol(X.eval) != ncol(X) || any(!is.finite(X.eval))) {
        stop("'X.eval' must be a finite matrix with ncol(X.eval) = ncol(X).",
             call. = FALSE)
    }
    coordinate.method <- match.arg(coordinate.method)
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    support.grid <- .klp.clean.support.grid(support.grid, nrow(X))
    degree.grid <- .klp.clean.degree.grid(degree.grid)
    kernel.grid <- .klp.clean.kernel.grid(kernel.grid)
    foldid <- .klp.prepare.foldid(nrow(X), foldid, cv.folds, cv.seed)

    cand <- expand.grid(
        support.size = support.grid,
        degree = degree.grid,
        kernel = kernel.grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    cv.result <- .klp.cv.table(
        X = X,
        y = y,
        foldid = foldid,
        cand = cand,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
    cv.table <- cv.result$cv.table
    best.idx <- order(
        cv.table$cv.rmse.observed,
        cv.table$support.size,
        cv.table$degree,
        cv.table$kernel
    )[[1L]]
    selected <- cv.table[best.idx, , drop = FALSE]
    selected.dim <- .klp.resolve.chart.dim(
        X = X,
        support.size = selected$support.size[[1L]],
        degree = selected$degree[[1L]],
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
    fitted <- .klp.predict.local.polynomial(
        X.train = X,
        y.train = y,
        X.eval = X.eval,
        support.size = selected$support.size[[1L]],
        degree = selected$degree[[1L]],
        kernel = selected$kernel[[1L]],
        coordinate.method = coordinate.method,
        chart.dim = selected.dim$chart.dim
    )
    out <- list(
        method.id = "kernel_local_polynomial_cv",
        method.family = "kernel_local_polynomial",
        X = X,
        y = y,
        X.eval = X.eval,
        fitted.values = fitted,
        selected = selected,
        cv.table = cv.table,
        foldid = foldid,
        coordinate.method = coordinate.method,
        requested.chart.dim = chart.dim,
        chart.dim = selected.dim$chart.dim,
        auto.chart.dim = identical(chart.dim, "auto"),
        auto.chart.dim.diagnostics = selected.dim$diagnostics,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        call = match.call()
    )
    class(out) <- c("kernel.local.polynomial.cv", "list")
    out
}

#' @export
predict.kernel.local.polynomial.cv <- function(object, newdata = NULL, ...) {
    dots <- list(...)
    if (length(dots)) {
        stop("Unused arguments: ", paste(names(dots), collapse = ", "),
             call. = FALSE)
    }
    X.eval <- if (is.null(newdata)) object$X.eval else as.matrix(newdata)
    .klp.predict.local.polynomial(
        X.train = object$X,
        y.train = object$y,
        X.eval = X.eval,
        support.size = object$selected$support.size[[1L]],
        degree = object$selected$degree[[1L]],
        kernel = object$selected$kernel[[1L]],
        coordinate.method = object$coordinate.method,
        chart.dim = object$chart.dim
    )
}

#' @export
print.kernel.local.polynomial.cv <- function(x, ...) {
    cat("Kernel local polynomial CV fit\n")
    cat("  observations:", nrow(x$X), "\n")
    cat("  coordinate method:", x$coordinate.method, "\n")
    cat("  selected support.size:", x$selected$support.size[[1L]], "\n")
    cat("  selected degree:", x$selected$degree[[1L]], "\n")
    cat("  selected kernel:", x$selected$kernel[[1L]], "\n")
    cat("  selected CV RMSE:",
        signif(x$selected$cv.rmse.observed[[1L]], 5), "\n")
    invisible(x)
}

.klp.rmse <- function(x, y) {
    sqrt(mean((as.numeric(x) - as.numeric(y))^2, na.rm = TRUE))
}

.klp.cv.table <- function(X, y, foldid, cand, coordinate.method, chart.dim,
                          auto.chart.support.metric,
                          auto.chart.selection.metric) {
    cand$chart.dim <- NA_integer_
    dim.lookup <- list()
    combos <- unique(cand[, c("support.size", "degree"), drop = FALSE])
    for (ii in seq_len(nrow(combos))) {
        info <- .klp.resolve.chart.dim(
            X = X,
            support.size = combos$support.size[[ii]],
            degree = combos$degree[[ii]],
            coordinate.method = coordinate.method,
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        key <- paste(combos$support.size[[ii]], combos$degree[[ii]], sep = "_")
        dim.lookup[[key]] <- info$chart.dim
    }
    for (rr in seq_len(nrow(cand))) {
        key <- paste(cand$support.size[[rr]], cand$degree[[rr]], sep = "_")
        cand$chart.dim[[rr]] <- dim.lookup[[key]]
    }
    pred <- matrix(NA_real_, nrow = length(y), ncol = nrow(cand))
    support.sizes <- sort(unique(cand$support.size))
    max.support.size <- max(support.sizes)
    for (fold in sort(unique(foldid))) {
        test <- which(foldid == fold)
        train <- which(foldid != fold)
        X.train <- X[train, , drop = FALSE]
        y.train <- y[train]
        fold.max.support <- min(max.support.size, length(train))
        for (ii in seq_along(test)) {
            target <- test[[ii]]
            center <- X[target, , drop = TRUE]
            ordered <- .klp.local.order(
                X.train = X.train,
                center = center,
                support.size = fold.max.support
            )
            for (support.size in support.sizes) {
                support.rows <- which(cand$support.size == support.size)
                max.chart.dim <- max(cand$chart.dim[support.rows],
                                     na.rm = TRUE)
                local <- .klp.local.neighborhood.from.order(
                    X.train = X.train,
                    y.train = y.train,
                    center = center,
                    ordered = ordered,
                    support.size = support.size,
                    coordinate.method = coordinate.method,
                    chart.dim = max.chart.dim
                )
                kernel.weights <- lapply(
                    unique(cand$kernel[support.rows]),
                    function(kernel) .klp.kernel.weights(local$distances, kernel)
                )
                names(kernel.weights) <- unique(cand$kernel[support.rows])
                design.cache <- new.env(parent = emptyenv())
                for (rr in support.rows) {
                    w <- kernel.weights[[cand$kernel[[rr]]]]
                    pred[target, rr] <- .klp.fit.intercept.lazy(
                        z = local$z,
                        y = local$y,
                        weights = w,
                        degree = cand$degree[[rr]],
                        chart.dim = cand$chart.dim[[rr]],
                        design.cache = design.cache
                    )
                }
            }
        }
    }
    cv.table <- cand
    cv.table$cv.rmse.observed <- vapply(
        seq_len(ncol(pred)),
        function(j) .klp.rmse(pred[, j], y),
        numeric(1L)
    )
    list(cv.table = cv.table, predictions = pred)
}

.klp.clean.support.grid <- function(support.grid, n) {
    out <- sort(unique(as.integer(support.grid)))
    out <- out[is.finite(out) & out >= 2L & out <= n]
    if (!length(out)) {
        stop("'support.grid' has no valid support sizes.", call. = FALSE)
    }
    out
}

.klp.clean.degree.grid <- function(degree.grid) {
    out <- sort(unique(as.integer(degree.grid)))
    out <- out[is.finite(out) & out %in% 0:2]
    if (!length(out)) {
        stop("'degree.grid' must contain at least one of 0, 1, or 2.",
             call. = FALSE)
    }
    out
}

.klp.clean.kernel.grid <- function(kernel.grid) {
    allowed <- c("gaussian", "tricube", "epanechnikov", "triangular")
    out <- unique(as.character(kernel.grid))
    out <- out[nzchar(out)]
    if (!length(out) || any(!out %in% allowed)) {
        stop("'kernel.grid' contains unsupported kernels.", call. = FALSE)
    }
    out
}

.klp.prepare.foldid <- function(n, foldid, cv.folds, cv.seed) {
    if (!is.null(foldid)) {
        if (!is.numeric(foldid) || length(foldid) != n ||
            any(is.na(foldid)) || any(foldid != as.integer(foldid)) ||
            any(foldid < 1L)) {
            stop("'foldid' must be a positive integer vector of length nrow(X).",
                 call. = FALSE)
        }
        return(as.integer(foldid))
    }
    cv.folds <- as.integer(cv.folds)
    if (!is.finite(cv.folds) || cv.folds < 2L || cv.folds > n) {
        stop("'cv.folds' must be an integer between 2 and nrow(X).",
             call. = FALSE)
    }
    set.seed(cv.seed)
    sample(rep(seq_len(cv.folds), length.out = n))
}

.klp.resolve.chart.dim <- function(X, support.size, degree, coordinate.method,
                                   chart.dim, auto.chart.support.metric,
                                   auto.chart.selection.metric) {
    if (identical(coordinate.method, "coordinates")) {
        if (!is.null(chart.dim) &&
            !(length(chart.dim) == 1L && is.numeric(chart.dim) &&
              as.integer(chart.dim) == ncol(X))) {
            stop("'chart.dim' must be NULL when coordinate.method = 'coordinates'.",
                 call. = FALSE)
        }
        return(list(chart.dim = ncol(X), diagnostics = NULL))
    }
    if (is.null(chart.dim)) {
        return(list(chart.dim = ncol(X), diagnostics = NULL))
    }
    if (identical(chart.dim, "auto")) {
        diagnostics <- .local.pca.auto.chart.dim.with.metric(
            X = X,
            support.size = support.size,
            degree = degree,
            operator.support.metric = "coordinates",
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        return(list(chart.dim = diagnostics$chart.dim,
                    diagnostics = diagnostics))
    }
    dim <- as.integer(chart.dim)
    if (!is.finite(dim) || dim < 1L || dim > ncol(X)) {
        stop("'chart.dim' must be between 1 and ncol(X), or 'auto'.",
             call. = FALSE)
    }
    list(chart.dim = dim, diagnostics = NULL)
}

.klp.predict.local.polynomial <- function(X.train, y.train, X.eval,
                                          support.size, degree, kernel,
                                          coordinate.method, chart.dim) {
    X.train <- as.matrix(X.train)
    X.eval <- as.matrix(X.eval)
    y.train <- as.numeric(y.train)
    support.size <- min(as.integer(support.size), nrow(X.train))
    out <- rep(NA_real_, nrow(X.eval))
    for (i in seq_len(nrow(X.eval))) {
        center <- X.eval[i, , drop = TRUE]
        d <- sqrt(rowSums((X.train -
            matrix(center, nrow(X.train), ncol(X.train), byrow = TRUE))^2))
        idx <- order(d, seq_along(d))[seq_len(support.size)]
        local.d <- d[idx]
        weights <- .klp.kernel.weights(local.d, kernel)
        if (!any(weights > 0)) weights[] <- 1
        z <- .klp.local.coordinates(
            X.support = X.train[idx, , drop = FALSE],
            center = center,
            coordinate.method = coordinate.method,
            chart.dim = chart.dim
        )
        design <- .malps.design.matrix(z, degree)
        ok <- is.finite(y.train[idx]) & is.finite(weights) & weights > 0
        if (sum(ok) < ncol(design)) {
            out[[i]] <- stats::weighted.mean(y.train[idx], weights,
                                             na.rm = TRUE)
            next
        }
        fit <- tryCatch(
            stats::lm.wfit(design[ok, , drop = FALSE],
                           y.train[idx][ok], weights[ok]),
            error = function(e) NULL
        )
        out[[i]] <- if (is.null(fit) || !length(fit$coefficients) ||
                         !is.finite(fit$coefficients[[1L]])) {
            stats::weighted.mean(y.train[idx], weights, na.rm = TRUE)
        } else {
            fit$coefficients[[1L]]
        }
    }
    out
}

.klp.local.neighborhood <- function(X.train, y.train, center, support.size,
                                    coordinate.method, chart.dim) {
    ordered <- .klp.local.order(
        X.train = X.train,
        center = center,
        support.size = support.size
    )
    .klp.local.neighborhood.from.order(
        X.train = X.train,
        y.train = y.train,
        center = center,
        ordered = ordered,
        support.size = support.size,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim
    )
}

.klp.local.order <- function(X.train, center, support.size) {
    d <- sqrt(rowSums((X.train -
        matrix(center, nrow(X.train), ncol(X.train), byrow = TRUE))^2))
    idx <- order(d, seq_along(d))[seq_len(min(as.integer(support.size),
                                             nrow(X.train)))]
    list(index = idx, distances = d[idx])
}

.klp.local.neighborhood.from.order <- function(X.train, y.train, center,
                                               ordered, support.size,
                                               coordinate.method, chart.dim) {
    support.size <- min(as.integer(support.size), length(ordered$index))
    idx <- ordered$index[seq_len(support.size)]
    distances <- ordered$distances[seq_len(support.size)]
    z <- .klp.local.coordinates(
        X.support = X.train[idx, , drop = FALSE],
        center = center,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim
    )
    list(
        index = idx,
        distances = distances,
        y = y.train[idx],
        z = z
    )
}

.klp.fit.intercept <- function(z, y, weights, degree) {
    .klp.fit.intercept.lazy(
        z = z,
        y = y,
        weights = weights,
        degree = degree,
        chart.dim = ncol(z),
        design.cache = new.env(parent = emptyenv())
    )
}

.klp.fit.intercept.lazy <- function(z, y, weights, degree, chart.dim,
                                    design.cache) {
    ok <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(weights > 0)) {
        weights[] <- 1
        ok <- is.finite(y) & is.finite(weights) & weights > 0
    }
    n.design <- .klp.design.ncol(degree, chart.dim)
    if (sum(ok) < n.design) {
        return(stats::weighted.mean(y, weights, na.rm = TRUE))
    }
    design <- .klp.get.local.design(z, degree, chart.dim, design.cache)
    .klp.fit.intercept.design(design, y, weights)
}

.klp.fit.intercept.design <- function(design, y, weights) {
    ok <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(weights > 0)) {
        weights[] <- 1
        ok <- is.finite(y) & is.finite(weights) & weights > 0
    }
    if (sum(ok) < ncol(design)) {
        return(stats::weighted.mean(y, weights, na.rm = TRUE))
    }
    fit <- tryCatch(
        stats::lm.wfit(design[ok, , drop = FALSE], y[ok], weights[ok]),
        error = function(e) NULL
    )
    if (is.null(fit) || !length(fit$coefficients) ||
        !is.finite(fit$coefficients[[1L]])) {
        stats::weighted.mean(y, weights, na.rm = TRUE)
    } else {
        fit$coefficients[[1L]]
    }
}

.klp.design.ncol <- function(degree, chart.dim) {
    degree <- as.integer(degree)
    chart.dim <- as.integer(chart.dim)
    if (degree == 0L) return(1L)
    if (degree == 1L) return(1L + chart.dim)
    if (degree == 2L) return(1L + chart.dim + chart.dim * (chart.dim + 1L) / 2L)
    stop("Unsupported local polynomial degree: ", degree, call. = FALSE)
}

.klp.design.cache.key <- function(degree, chart.dim) {
    paste(as.integer(degree), as.integer(chart.dim), sep = "_")
}

.klp.get.local.design <- function(z, degree, chart.dim, design.cache) {
    key <- .klp.design.cache.key(degree, chart.dim)
    if (!exists(key, envir = design.cache, inherits = FALSE)) {
        design <- .malps.design.matrix(
            z[, seq_len(chart.dim), drop = FALSE],
            degree
        )
        assign(key, design, envir = design.cache)
    }
    get(key, envir = design.cache, inherits = FALSE)
}

.klp.local.design.cache <- function(z, cand, rows) {
    combos <- unique(cand[rows, c("degree", "chart.dim"), drop = FALSE])
    out <- new.env(parent = emptyenv())
    for (ii in seq_len(nrow(combos))) {
        dim <- combos$chart.dim[[ii]]
        degree <- combos$degree[[ii]]
        .klp.get.local.design(z, degree, dim, out)
    }
    out
}

.klp.local.coordinates <- function(X.support, center, coordinate.method,
                                   chart.dim) {
    centered <- sweep(X.support, 2L, center, "-")
    if (identical(coordinate.method, "coordinates")) return(centered)
    sv <- svd(centered, nu = 0L, nv = chart.dim)
    basis <- sv$v[, seq_len(chart.dim), drop = FALSE]
    centered %*% basis
}

.klp.kernel.weights <- function(distances, kernel) {
    if (!length(distances)) return(numeric(0))
    h <- max(distances[is.finite(distances)], 0)
    if (!is.finite(h) || h <= 0) h <- 1
    u <- as.numeric(distances) / (h + sqrt(.Machine$double.eps))
    w <- switch(
        kernel,
        gaussian = exp(-0.5 * u^2),
        tricube = ifelse(u < 1, (1 - u^3)^3, 0),
        epanechnikov = pmax(0, 1 - u^2),
        triangular = pmax(0, 1 - u)
    )
    w[!is.finite(w)] <- 0
    as.numeric(w)
}
