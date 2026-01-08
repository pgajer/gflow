#' Cluster cell trajectories using DTW distance on vertex features and HDBSCAN
#'
#' @description
#' For a fixed gradient-flow cell (m, M), each trajectory induces a sequence of feature vectors
#' taken from a vertex-feature matrix Z. Pairwise trajectory distances are computed using
#' Dynamic Time Warping (DTW) and clustered via HDBSCAN.
#'
#' @param traj.list List of integer vectors; each element is a trajectory as vertex indices (1-based).
#' @param Z Numeric matrix of vertex features with nrow(Z) == number of vertices.
#' @param feature.cols Optional integer vector selecting columns of Z to use.
#' @param scale.features Logical; scale each feature column before extracting sequences.
#' @param dtw.window.type DTW global constraint (passed to dtw::dtw()).
#' @param dtw.window.size Optional window size parameter for dtw::dtw().
#' @param dtw.step.pattern Step pattern object for dtw::dtw() (e.g., dtw::symmetric2).
#' @param minPts HDBSCAN minPts parameter.
#' @param n.cores Number of parallel workers for pairwise DTW computation (uses parallel).
#' @param verbose Logical.
#' @return A list with fields:
#'   \item{dist}{A dist object of DTW distances.}
#'   \item{hdbscan}{The dbscan::hdbscan() result.}
#'   \item{cluster}{Integer vector of cluster labels (0 = noise).}
#'
#' @export
cluster.cell.trajectories.dtw <- function(traj.list,
                                          Z,
                                          feature.cols = NULL,
                                          scale.features = TRUE,
                                          dtw.window.type = "sakoechiba",
                                          dtw.window.size = NULL,
                                          dtw.step.pattern = dtw::symmetric2,
                                          minPts = 5L,
                                          n.cores = 1L,
                                          verbose = TRUE) {

    if (!is.list(traj.list) || length(traj.list) < 2L) {
        stop("traj.list must be a list of at least 2 trajectories.")
    }
    if (!is.matrix(Z) || !is.numeric(Z)) {
        stop("Z must be a numeric matrix.")
    }
    if (!is.null(feature.cols)) {
        Z.use <- Z[, feature.cols, drop = FALSE]
    } else {
        Z.use <- Z
    }

    if (scale.features) {
        ## scale columns globally (not per-trajectory)
        Z.use <- scale(Z.use)
    }

    ## build per-trajectory feature matrices
    seq.list <- lapply(traj.list, function(v) {
        v <- as.integer(v)
        if (any(v < 1L) || any(v > nrow(Z.use))) {
            stop("Trajectory contains vertex indices outside 1..nrow(Z).")
        }
        Z.use[v, , drop = FALSE]
    })

    n <- length(seq.list)
    d.vec <- numeric(n * (n - 1L) / 2L)
    idx <- 1L

    ## pair indexing helper
    pairs <- utils::combn(n, 2L)

    if (verbose) {
        message("Computing DTW distances for ", n, " trajectories (", ncol(pairs), " pairs).")
    }

    if (n.cores <= 1L) {
        for (k in seq_len(ncol(pairs))) {
            i <- pairs[1L, k]
            j <- pairs[2L, k]
            d.vec[idx] <- traj.dtw.distance(seq.list[[i]], seq.list[[j]],
                                            window.type = dtw.window.type,
                                            window.size = dtw.window.size,
                                            step.pattern = dtw.step.pattern,
                                            normalize = TRUE)
            idx <- idx + 1L
        }
    } else {
        ## CRAN-safe: use parallel::mclapply only on Unix; otherwise fall back to PSOCK.
        worker.fun <- function(k) {
            i <- pairs[1L, k]
            j <- pairs[2L, k]

            out <- tryCatch(
                traj.dtw.distance(seq.list[[i]], seq.list[[j]],
                                  window.type = dtw.window.type,
                                  window.size = dtw.window.size,
                                  step.pattern = dtw.step.pattern,
                                  normalize = TRUE),
                error = function(e) NA_real_
            )

            ## Defensive: guarantee scalar numeric
            out <- as.numeric(out)[1L]
            if (!is.finite(out)) NA_real_ else out
        }

        ## if (.Platform$OS.type == "unix") {
        ##     d.list <- parallel::mclapply(seq_len(ncol(pairs)), worker.fun,
        ##                                  mc.cores = n.cores)
        ## } else {
        ##     cl <- parallel::makeCluster(n.cores)
        ##     on.exit(parallel::stopCluster(cl), add = TRUE)
        ##     d.list <- parallel::parLapply(cl, seq_len(ncol(pairs)), worker.fun)
        ## }

        if (n.cores <= 1L) {
            d.vec <- vapply(seq_len(ncol(pairs)), worker.fun, numeric(1))
        } else {
            cl <- parallel::makeCluster(n.cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
            d.list <- parallel::parLapply(cl, seq_len(ncol(pairs)), worker.fun)
            d.vec <- unlist(d.list, use.names = FALSE)
        }

        is.try.error <- vapply(d.list, inherits, logical(1), what = "try-error")
        sum(is.try.error)
        if (any(is.try.error)) {
            cat("First DTW error:\n")
            print(d.list[[which(is.try.error)[1L]]])
        }

        d.vec <- unlist(d.list, use.names = FALSE)

        expected.len <- n * (n - 1L) / 2L
        if (length(d.vec) != expected.len) {
            stop("Internal error: wrong d.vec length: expected ", expected.len,
                 ", got ", length(d.vec), ".")
        }

        if (verbose) {
            cat("typeof(d.vec) = ", typeof(d.vec), "\n")
            cat("anyNA(d.vec) = ", anyNA(d.vec), "\n")
            cat("any(!is.finite(d.vec)) = ", any(!is.finite(d.vec)), "\n")
            if (is.numeric(d.vec)) {
                cat("range(d.vec) = ", paste(range(d.vec, na.rm = TRUE), collapse = " .. "), "\n")
            }
        }

        if (!is.numeric(d.vec)) {
            stop("Internal error: DTW worker returned non-numeric results; d.vec is not numeric.")
        }
        if (anyNA(d.vec) || any(!is.finite(d.vec))) {
            stop("DTW produced NA/Inf distances; inspect DTW parameters or enable error diagnostics.")
        }
    }

    D <- .make.dist.from.lower(d.vec, n, method = "DTW")

    stopifnot(inherits(D, "dist"))

    if (verbose) {
        cat("typeof(D) =", typeof(D), "\n")          ## should be "double"
        cat("len(D)    =", length(D), "\n")          ## should be 7875
        cat("Size(D)   =", attr(D, "Size"), "\n")    ## should be 126

        v <- as.vector(D)
        cat("typeof(v) =", typeof(v), "\n")          ## should be "double"
        cat("finite    =", all(is.finite(v)), "\n")
        cat("nonneg    =", all(v >= 0), "\n")
        cat("range     =", paste(range(v), collapse = " .. "), "\n")
    }

    ## HDBSCAN supports dist objects directly
    hdb <- dbscan::hdbscan(D, minPts = as.integer(minPts))

    list(dist = D, hdbscan = hdb, cluster = hdb$cluster)
}
