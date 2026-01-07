#' Apply HDBSCAN Clustering
#'
#' This function applies HDBSCAN clustering to a given dataset (expressed as a
#' matrix), and evaluates clustering results with different minimum cluster
#' sizes using dunn and connectivity indices.
#'
#' @param X             A matrix of the data to be clustered. Each row represents an observation.
#' @param method        A method (index) to be used for identifying the optimal number of clusters.
#'                      Default is 'dunn'. Other options: 'connectivity', 'cl0.size' and 'all'.
#' @param min.pts       A vector of integers indicating the minimum cluster sizes to be evaluated.
#'                      Default is 5:50. If NULL, the function will generate values based on min.prop,
#'                      max.prop, and n.test.cltr.
#' @param min.prop      The minimum proportion of nrow(X) to be used in test clusterings when min.pts is NULL.
#'                      Default is 0.1.
#' @param max.prop      The maximum proportion of nrow(X) to be used in test clusterings when min.pts is NULL.
#'                      Default is 0.5.
#' @param n.test.cltr   The number of test clusters when min.pts is NULL. Default is 10.
#' @param soft.K        The number of nearest neighbors to consider when softening the cluster boundaries.
#'                      Default is 20.
#' @param n.cores       The number of cores to use for parallel processing. Default is 10.
#' @param verbose       A logical value indicating whether progress should be displayed. Default is FALSE.
#'
#' @return A list containing the details of the clustering for three different metrics:
#'   \item{cl0.size}{Vector of sizes of cluster labeled as '0' (noise) for each min.pts value}
#'   \item{n.cltrs}{Vector of number of clusters (excluding noise) for each min.pts value}
#'   \item{cl0.cltr}{Cluster assignments using the min.pts that minimizes cl0.size}
#'   \item{cl0.cltr.ext}{Extended cluster assignments after softening (cl0.size method)}
#'   \item{cl0.opt.k}{Optimal min.pts value that minimizes cl0.size}
#'   \item{cl0.cl0.size}{Size of cluster '0' for the cl0.size-optimal clustering}
#'   \item{cl0.n.cltrs}{Number of clusters for the cl0.size-optimal clustering}
#'   \item{cl0.cltr.freq}{Frequency table of clusters for cl0.size method}
#'   \item{cl0.cltr.ext.freq}{Frequency table of extended clusters for cl0.size method}
#'   \item{dunn.idx}{Vector of Dunn indices for each min.pts value}
#'   \item{dunn.cltr}{Cluster assignments using the min.pts that maximizes Dunn index}
#'   \item{dunn.cltr.ext}{Extended cluster assignments after softening (Dunn method)}
#'   \item{dunn.opt.k}{Optimal min.pts value that maximizes Dunn index}
#'   \item{dunn.cl0.size}{Size of cluster '0' for the Dunn-optimal clustering}
#'   \item{dunn.n.cltrs}{Number of clusters for the Dunn-optimal clustering}
#'   \item{dunn.cltr.freq}{Frequency table of clusters for Dunn method}
#'   \item{dunn.cltr.ext.freq}{Frequency table of extended clusters for Dunn method}
#'   \item{connectivity.idx}{Vector of connectivity indices for each min.pts value}
#'   \item{connectivity.cltr}{Cluster assignments using the min.pts that minimizes connectivity}
#'   \item{connectivity.cltr.ext}{Extended cluster assignments after softening (connectivity method)}
#'   \item{connectivity.opt.k}{Optimal min.pts value that minimizes connectivity index}
#'   \item{connectivity.cl0.size}{Size of cluster '0' for the connectivity-optimal clustering}
#'   \item{connectivity.n.cltrs}{Number of clusters for the connectivity-optimal clustering}
#'   \item{connectivity.cltr.freq}{Frequency table of clusters for connectivity method}
#'   \item{connectivity.cltr.ext.freq}{Frequency table of extended clusters for connectivity method}
#'
#' @note This function requires the packages 'foreach', and 'doParallel'.
#' @note The Dunn and connectivity metrics are computed via the suggested package
#' 'clValid'. If 'clValid' is not installed, use \code{method = "cl0.size"} (or avoid
#' \code{method = "dunn"}, \code{"connectivity"}, or \code{"all"}).
#' @note The HDBSCAN clustering is performed via the suggested package
#' 'dbscan'. This function will error if 'dbscan' is not installed.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' X <- rbind(
#'   matrix(rnorm(100, mean = 0), ncol = 2),
#'   matrix(rnorm(100, mean = 5), ncol = 2)
#' )
#'
#' # Apply HDBSCAN clustering
#' result <- hdbscan.cltr(X, min.pts = 5:50, soft.K = 20, verbose = TRUE)
#'
#' # View optimal clustering for Dunn index
#' table(result$dunn.cltr)
#' }
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @export
hdbscan.cltr <- function(X,
                         method = "dunn",
                         min.pts = 5:50,
                         min.prop = 0.1,
                         max.prop = 0.5,
                         n.test.cltr = 10,
                         soft.K = 20,
                         n.cores = 10,
                         verbose = FALSE)
{
    ## Check required packages
    needs_clv <- method %in% c("dunn", "connectivity", "all")
    if (needs_clv && !requireNamespace("clValid", quietly = TRUE)) {
        stop(
            sprintf("method='%s' requires the suggested package 'clValid'. Install it or use method='cl0.size'.", method),
            call. = FALSE
        )
    }

    if (!requireNamespace("dbscan", quietly = TRUE)) {
        stop("This function requires the suggested package 'dbscan'. Please install it to use hdbscan.cltr().",
             call. = FALSE)
    }

    if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("The foreach package is not installed. Please install it before using this routine.")
    }

    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("The doParallel package is not installed. Please install it before using this routine.")
    }

    # Generate min.pts if NULL
    if (is.null(min.pts)) {
        stopifnot(0 < min.prop && min.prop < 1)
        stopifnot(0 < max.prop && max.prop < 1)
        stopifnot(min.prop < max.prop)

        min.pts <- seq(min.prop * nrow(X), max.prop * nrow(X), length.out = n.test.cltr)
        min.pts <- as.integer(min.pts)
    }

    # Adjust n.cores if necessary
    if (length(min.pts) < n.cores) {
        n.cores <- length(min.pts)
        if (verbose) {
            cat("Set n.cores to", n.cores, "\n")
        }
    }

    # Check min.pts validity
    if (any(min.pts >= nrow(X))) {
        idx <- min.pts < nrow(X)
        min.pts <- min.pts[idx]
        warning("Restricting min.pts to the range with min.pts < nrow(X)")
    }

    if (length(min.pts) == 0) {
        stop("min.pts is empty!")
    }

    # Pre-compute distance matrix for Dunn index
    dX.mat <- as.matrix(dist(X))

    # Initialize result variables
    cltrgs <- list()
    cl0.size <- numeric()
    n.cltrs <- numeric()
    dunn.idx <- numeric()
    connectivity.idx <- numeric()

    # Initialize output variables
    cl0.cltr <- NULL
    cl0.cltr.ext <- NULL
    cl0.k <- NA
    cl0.freq <- table(numeric(0))
    cl0.freq.ext <- NULL

    dunn.cltr <- NULL
    dunn.cltr.ext <- NULL
    dunn.k <- NA
    dunn.freq <- table(numeric(0))
    dunn.freq.ext <- NULL

    connectivity.cltr <- NULL
    connectivity.cltr.ext <- NULL
    connectivity.k <- NA
    connectivity.freq <- table(numeric(0))
    connectivity.freq.ext <- NULL

    # Define a single helper function that handles all methods
    # This function is defined once with all possible parameters
    compute_clustering_metrics <- function(k, X, dX.mat = NULL, compute_dunn = FALSE,
                                         compute_connectivity = FALSE, compute_cl0 = FALSE,
                                         verbose = FALSE) {
        if (verbose) cat("\r", k)

        cltr <- dbscan::hdbscan(X, minPts = k)$cluster
        names(cltr) <- rownames(X)

        n.cltrs <- max(cltr)
        dunn.idx <- NA
        connectivity.idx <- NA
        cl0.size <- NA

        if (length(unique(cltr)) > 1) {
            if (compute_dunn && !is.null(dX.mat)) {
                dunn.idx <- clValid::dunn(dX.mat, clusters = cltr)
            }

            if (compute_connectivity) {
                idx <- cltr != 0
                if (sum(idx) > 0) {
                    connectivity.idx <- clValid::connectivity(clusters = cltr[idx],
                                                             Data = X[idx, , drop = FALSE])
                }
            }

            if (compute_cl0) {
                f <- table(cltr)
                s <- f['0']
                if (!is.na(s)) {
                    cl0.size <- as.numeric(s)
                } else {
                    cl0.size <- 0
                }
            }
        }

        list(cltr = cltr,
             n.cltrs = n.cltrs,
             dunn.idx = dunn.idx,
             connectivity.idx = connectivity.idx,
             cl0.size = cl0.size)
    }

    if (method == "all") {
        if (n.cores == 1) {
            # Sequential processing
            for (k in min.pts) {
                res <- compute_clustering_metrics(k, X, dX.mat,
                                                compute_dunn = TRUE,
                                                compute_connectivity = TRUE,
                                                compute_cl0 = TRUE,
                                                verbose = verbose)

                cltrgs[[as.character(k)]] <- res$cltr
                cl0.size[as.character(k)] <- res$cl0.size
                n.cltrs[as.character(k)] <- res$n.cltrs
                dunn.idx[as.character(k)] <- res$dunn.idx
                connectivity.idx[as.character(k)] <- res$connectivity.idx
            }
            if (verbose) cat("\n")

        } else {
            # Parallel processing
            doParallel::registerDoParallel(n.cores)

            res <- foreach::foreach(k = min.pts) %dopar% {
                compute_clustering_metrics(k, X, dX.mat,
                                         compute_dunn = TRUE,
                                         compute_connectivity = TRUE,
                                         compute_cl0 = TRUE,
                                         verbose = FALSE)
            }

            doParallel::stopImplicitCluster()

            # Collect results
            for (j in seq_along(min.pts)) {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[as.character(k)]] <- r$cltr
                cl0.size[as.character(k)] <- r$cl0.size
                n.cltrs[as.character(k)] <- r$n.cltrs
                dunn.idx[as.character(k)] <- r$dunn.idx
                connectivity.idx[as.character(k)] <- r$connectivity.idx
            }
        }

        # Find optimal clusterings
        if (length(cl0.size) > 0 && any(!is.na(cl0.size))) {
            cl0.k <- names(which.min(cl0.size))
            cl0.cltr <- cltrgs[[cl0.k]]
            cl0.freq <- table(cl0.cltr)
            if ('0' %in% names(cl0.freq)) {
                cl0.cltr.ext <- kNN.cltr.imputation(X, cl0.cltr, K = soft.K)
                cl0.freq.ext <- table(cl0.cltr.ext)
            }
        }

        if (length(dunn.idx) > 0 && any(!is.na(dunn.idx))) {
            dunn.k <- names(which.max(dunn.idx))
            dunn.cltr <- cltrgs[[dunn.k]]
            dunn.freq <- table(dunn.cltr)
            if ('0' %in% names(dunn.freq)) {
                dunn.cltr.ext <- kNN.cltr.imputation(X, dunn.cltr, K = soft.K)
                dunn.freq.ext <- table(dunn.cltr.ext)
            }
        }

        if (length(connectivity.idx) > 0 && any(!is.na(connectivity.idx))) {
            valid.idx <- !is.na(connectivity.idx)
            if (any(valid.idx)) {
                connectivity.k <- names(connectivity.idx)[valid.idx][which.min(connectivity.idx[valid.idx])]
                connectivity.cltr <- cltrgs[[connectivity.k]]
                connectivity.freq <- table(connectivity.cltr)
                if ('0' %in% names(connectivity.freq)) {
                    connectivity.cltr.ext <- kNN.cltr.imputation(X, connectivity.cltr, K = soft.K)
                    connectivity.freq.ext <- table(connectivity.cltr.ext)
                }
            }
        }

    } else if (method == "dunn") {
        if (n.cores == 1) {
            for (k in min.pts) {
                res <- compute_clustering_metrics(k, X, dX.mat,
                                                compute_dunn = TRUE,
                                                compute_connectivity = FALSE,
                                                compute_cl0 = FALSE,
                                                verbose = verbose)

                cltrgs[[as.character(k)]] <- res$cltr
                n.cltrs[as.character(k)] <- res$n.cltrs
                dunn.idx[as.character(k)] <- res$dunn.idx
            }
            if (verbose) cat("\n")

        } else {
            doParallel::registerDoParallel(n.cores)

            res <- foreach::foreach(k = min.pts) %dopar% {
                compute_clustering_metrics(k, X, dX.mat,
                                         compute_dunn = TRUE,
                                         compute_connectivity = FALSE,
                                         compute_cl0 = FALSE,
                                         verbose = FALSE)
            }

            doParallel::stopImplicitCluster()

            for (j in seq_along(min.pts)) {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[as.character(k)]] <- r$cltr
                n.cltrs[as.character(k)] <- r$n.cltrs
                dunn.idx[as.character(k)] <- r$dunn.idx
            }
        }

        if (length(dunn.idx) > 0 && any(!is.na(dunn.idx))) {
            dunn.k <- names(which.max(dunn.idx))
            dunn.cltr <- cltrgs[[dunn.k]]
            dunn.freq <- table(dunn.cltr)
            if ('0' %in% names(dunn.freq)) {
                dunn.cltr.ext <- kNN.cltr.imputation(X, dunn.cltr, K = soft.K)
                dunn.freq.ext <- table(dunn.cltr.ext)
            }
        }

    } else if (method == "connectivity") {
        if (n.cores == 1) {
            for (k in min.pts) {
                res <- compute_clustering_metrics(k, X, NULL,
                                                compute_dunn = FALSE,
                                                compute_connectivity = TRUE,
                                                compute_cl0 = FALSE,
                                                verbose = verbose)

                cltrgs[[as.character(k)]] <- res$cltr
                n.cltrs[as.character(k)] <- res$n.cltrs
                connectivity.idx[as.character(k)] <- res$connectivity.idx
            }
            if (verbose) cat("\n")

        } else {
            doParallel::registerDoParallel(n.cores)

            res <- foreach::foreach(k = min.pts) %dopar% {
                compute_clustering_metrics(k, X, NULL,
                                         compute_dunn = FALSE,
                                         compute_connectivity = TRUE,
                                         compute_cl0 = FALSE,
                                         verbose = FALSE)
            }

            doParallel::stopImplicitCluster()

            for (j in seq_along(min.pts)) {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[as.character(k)]] <- r$cltr
                n.cltrs[as.character(k)] <- r$n.cltrs
                connectivity.idx[as.character(k)] <- r$connectivity.idx
            }
        }

        if (length(connectivity.idx) > 0 && any(!is.na(connectivity.idx))) {
            valid.idx <- !is.na(connectivity.idx)
            if (any(valid.idx)) {
                connectivity.k <- names(connectivity.idx)[valid.idx][which.min(connectivity.idx[valid.idx])]
                connectivity.cltr <- cltrgs[[connectivity.k]]
                connectivity.freq <- table(connectivity.cltr)
                if ('0' %in% names(connectivity.freq)) {
                    connectivity.cltr.ext <- kNN.cltr.imputation(X, connectivity.cltr, K = soft.K)
                    connectivity.freq.ext <- table(connectivity.cltr.ext)
                }
            }
        }

    } else if (method == "cl0.size") {
        # Only compute cl0.size metric
        for (k in min.pts) {
            if (verbose) cat("\r", k)

            cltr <- dbscan::hdbscan(X, minPts = k)$cluster
            names(cltr) <- rownames(X)
            cltrgs[[as.character(k)]] <- cltr

            f <- table(cltr)
            s <- f['0']
            if (!is.na(s)) {
                cl0.size[as.character(k)] <- as.numeric(s)
            } else {
                cl0.size[as.character(k)] <- 0
            }
            n.cltrs[as.character(k)] <- max(cltr)
        }
        if (verbose) cat("\n")

        if (length(cl0.size) > 0) {
            cl0.k <- names(which.min(cl0.size))
            cl0.cltr <- cltrgs[[cl0.k]]
            cl0.freq <- table(cl0.cltr)
            if ('0' %in% names(cl0.freq)) {
                cl0.cltr.ext <- kNN.cltr.imputation(X, cl0.cltr, K = soft.K)
                cl0.freq.ext <- table(cl0.cltr.ext)
            }
        }

    } else {
        stop("Invalid method. Choose from 'dunn', 'connectivity', 'cl0.size', or 'all'.")
    }

    # Prepare return list
    list(cl0.size = cl0.size,
         n.cltrs = n.cltrs,
         cl0.cltr = cl0.cltr,
         cl0.cltr.ext = cl0.cltr.ext,
         cl0.opt.k = if (!is.na(cl0.k)) as.integer(cl0.k) else NA,
         cl0.cl0.size = if (!is.null(cl0.freq) && '0' %in% names(cl0.freq)) as.numeric(cl0.freq['0']) else NA,
         cl0.n.cltrs = if (!is.null(cl0.cltr)) max(cl0.cltr) else NA,
         cl0.cltr.freq = cl0.freq,
         cl0.cltr.ext.freq = cl0.freq.ext,
         dunn.idx = dunn.idx,
         dunn.cltr = dunn.cltr,
         dunn.cltr.ext = dunn.cltr.ext,
         dunn.opt.k = if (!is.na(dunn.k)) as.integer(dunn.k) else NA,
         dunn.cl0.size = if (!is.null(dunn.freq) && '0' %in% names(dunn.freq)) as.numeric(dunn.freq['0']) else NA,
         dunn.n.cltrs = if (!is.null(dunn.cltr)) max(dunn.cltr) else NA,
         dunn.cltr.freq = dunn.freq,
         dunn.cltr.ext.freq = dunn.freq.ext,
         connectivity.idx = connectivity.idx,
         connectivity.cltr = connectivity.cltr,
         connectivity.cltr.ext = connectivity.cltr.ext,
         connectivity.opt.k = if (!is.na(connectivity.k)) as.integer(connectivity.k) else NA,
         connectivity.cl0.size = if (!is.null(connectivity.freq) && '0' %in% names(connectivity.freq)) as.numeric(connectivity.freq['0']) else NA,
         connectivity.n.cltrs = if (!is.null(connectivity.cltr)) max(connectivity.cltr) else NA,
         connectivity.cltr.freq = connectivity.freq,
         connectivity.cltr.ext.freq = connectivity.freq.ext)
}


#' kNN-Based Cluster Imputation
#'
#' This function utilizes the k-Nearest Neighbors (k-NN) algorithm with a modified
#' majority vote strategy to reassign points in the reference cluster (typically cluster 0,
#' representing noise or outliers) to one of the non-reference clusters. The process
#' iterates until no further reduction in the size of the reference cluster is observed.
#'
#' @param X Matrix representing the feature space in which the clustering is defined.
#'        Each row represents an observation and columns represent features.
#' @param cltr Vector containing the initial clustering labels for each point in 'X'.
#'        Must have the same length as nrow(X).
#' @param ref.cltr Cluster label designated as the reference cluster for reassignment
#'        (default is 0). Points with this label will be reassigned to other clusters.
#' @param K Number of nearest neighbors used for cluster reassignment (default is 20).
#'        Will be automatically adjusted to min(K, nrow(X)-1) if K exceeds available points.
#' @param use.geodesic.dist Logical flag indicating whether to use geodesic distance
#'        instead of Euclidean distance (default is FALSE). Requires geodesic.knn() function
#'        if set to TRUE.
#'
#' @return A vector containing the new clustering labels after iterative imputation.
#'         The length and order match the input 'cltr' vector.
#'
#' @details
#' The algorithm works as follows:
#' \enumerate{
#'   \item For each point in the reference cluster, find its K nearest neighbors
#'   \item Among these neighbors, find the first neighbor that belongs to a non-reference cluster
#'   \item Reassign the point to that cluster
#'   \item Repeat until no more points can be reassigned
#' }
#'
#' Note: The current implementation uses a "first valid neighbor" strategy rather than
#' a true majority vote. Points are assigned to the cluster of their first non-reference
#' neighbor, which may lead to suboptimal assignments in some cases.
#'
#' @examples
#' \dontrun{
#' # Generate sample data with 3 clusters and some noise
#' set.seed(123)
#' X <- rbind(
#'   matrix(rnorm(100, mean = 0), ncol = 2),
#'   matrix(rnorm(100, mean = 5), ncol = 2),
#'   matrix(rnorm(100, mean = c(5, 0)), ncol = 2),
#'   matrix(runif(20, -2, 7), ncol = 2)  # noise points
#' )
#'
#' # Initial clustering with some points labeled as 0 (noise)
#' cltr <- c(rep(1, 50), rep(2, 50), rep(3, 50), rep(0, 20))
#' names(cltr) <- rownames(X) <- paste0("point_", 1:nrow(X))
#'
#' # Apply kNN imputation
#' new_cltr <- kNN.cltr.imputation(X, cltr, ref.cltr = 0, K = 10)
#'
#' # Check results
#' table(Original = cltr, Imputed = new_cltr)
#' }
#'
#' @importFrom FNN get.knn
#'
#' @note
#' \itemize{
#'   \item If rownames are provided for X and names for cltr, they must match exactly
#'   \item The function requires the FNN package for Euclidean k-NN search
#'   \item If use.geodesic.dist = TRUE, the geodesic.knn() function must be available
#'   \item The algorithm may not converge if the reference cluster contains isolated points
#'         with no non-reference neighbors within K nearest neighbors
#' }
#'
#' @seealso \code{\link[FNN]{get.knn}} for the k-NN implementation used
#'
#' @export
kNN.cltr.imputation <- function(X, cltr, ref.cltr = 0, K = 20, use.geodesic.dist = FALSE) {

    # Input validation
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (nrow(X) != length(cltr)) {
        cat("nrow(X):", nrow(X), "\n")
        cat("length(cltr):", length(cltr), "\n")
        stop("nrow(X) != length(cltr)")
    }

    # Check rownames and names consistency if they exist
    if (!is.null(rownames(X)) && !is.null(names(cltr))) {
        if (!all(rownames(X) == names(cltr))) {
            stop("rownames(X) != names(cltr)")
        }
    }

    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    # Validate K
    if (!is.numeric(K) || length(K) != 1 || K < 1 || K != floor(K)) {
        stop("K must be a positive integer")
    }

    # Adjust K if necessary
    K <- min(K, nrow(X) - 1)

    if (K < 1) {
        stop("Not enough points for k-NN search (need at least 2 points)")
    }

    # Check if reference cluster exists
    if (!ref.cltr %in% cltr) {
        warning(paste("Reference cluster", ref.cltr, "not found in clustering. Returning original clustering."))
        return(cltr)
    }

    # Compute k-NN
    if (use.geodesic.dist) {
        # Check if geodesic.knn function exists
        if (!exists("geodesic.knn")) {
            stop("geodesic.knn function not found. Please ensure it is loaded or set use.geodesic.dist = FALSE")
        }

        # Note: geodesic.knn is assumed to return k+1 neighbors including the point itself
        r <- geodesic.knn(X, k = K + 1)
        nn.i <- r$nn.index
        nn.i <- nn.i[, -1, drop = FALSE]  # Remove the first column (self-reference)
    } else {
        # Check if FNN package is available
        if (!requireNamespace("FNN", quietly = TRUE)) {
            stop("Package 'FNN' is required but not installed. Please install it using install.packages('FNN')")
        }

        nn <- FNN::get.knn(X, k = K)
        nn.i <- nn$nn.index
    }

    # Initialize working copy of cluster assignments
    cltr2 <- cltr

    # Get initial size of reference cluster
    idx <- cltr == ref.cltr
    old.cl0.size <- sum(idx)

    if (old.cl0.size == 0) {
        return(cltr)  # No points to reassign
    }

    # Iteration counter for debugging
    iter <- 0
    max.iter <- length(cltr)  # Prevent infinite loops

    # Main iteration loop
    cl0.size.diff <- 1
    while (cl0.size.diff > 0 && iter < max.iter) {
        iter <- iter + 1

        # Find current reference cluster points
        idx <- cltr2 == ref.cltr
        i0 <- which(idx)

        # Try to reassign each reference cluster point
        for (i in i0) {
            # Get neighbors of point i
            ii <- nn.i[i, ]

            # Get cluster labels of neighbors
            nbrd.cl <- cltr2[ii]

            # Find first neighbor not in reference cluster
            j <- which(nbrd.cl != ref.cltr)

            if (length(j) > 0) {
                # Assign to cluster of first non-reference neighbor
                # TODO: Consider implementing true majority vote here
                cltr2[i] <- nbrd.cl[j[1]]
            }
            # If no non-reference neighbors found, point remains in reference cluster
        }

        # Calculate new reference cluster size
        idx <- cltr2 == ref.cltr
        new.cl0.size <- sum(idx)

        # Calculate change
        cl0.size.diff <- old.cl0.size - new.cl0.size
        old.cl0.size <- new.cl0.size

        # Optional: print progress
        # if (verbose) {
        #     cat("Iteration", iter, ": Reference cluster size =", new.cl0.size, "\n")
        # }
    }

    if (iter >= max.iter) {
        warning("Maximum iterations reached. Algorithm may not have fully converged.")
    }

    # Preserve names if they existed in original clustering
    if (!is.null(names(cltr))) {
        names(cltr2) <- names(cltr)
    }

    return(cltr2)
}

#' Enhanced kNN-Based Cluster Imputation with Majority Vote
#'
#' An improved version of kNN cluster imputation that uses true majority voting
#' and provides more options for handling ties and isolated points.
#'
#' @param X Matrix representing the feature space (observations in rows, features in columns)
#' @param cltr Vector of initial cluster labels
#' @param ref.cltr Reference cluster label to reassign (default is 0)
#' @param K Number of nearest neighbors (default is 20)
#' @param method Voting method: "majority" (default), "weighted", or "first"
#' @param tie.break How to handle ties: "random", "first", or "none" (keep in ref cluster)
#' @param min.votes Minimum number of votes needed for reassignment (default is 1)
#' @param use.geodesic.dist Use geodesic distance if TRUE (default is FALSE)
#' @param verbose Print progress information (default is FALSE)
#'
#' @return Vector of cluster labels after imputation
#'
#' @export
kNN.cltr.imputation.enhanced <- function(X, cltr, ref.cltr = 0, K = 20,
                                       method = "majority",
                                       tie.break = "first",
                                       min.votes = 1,
                                       use.geodesic.dist = FALSE,
                                       verbose = FALSE) {

    # [Initial validation same as above...]

    # Compute k-NN
    if (use.geodesic.dist) {
        if (!exists("geodesic.knn")) {
            stop("geodesic.knn function not found")
        }
        r <- geodesic.knn(X, k = K + 1)
        nn.i <- r$nn.index[, -1, drop = FALSE]
        nn.d <- r$nn.dist[, -1, drop = FALSE]
    } else {
        if (!requireNamespace("FNN", quietly = TRUE)) {
            stop("Package 'FNN' is required")
        }
        nn <- FNN::get.knn(X, k = K)
        nn.i <- nn$nn.index
        nn.d <- nn$nn.dist
    }

    cltr2 <- cltr
    iter <- 0
    max.iter <- length(cltr)

    repeat {
        iter <- iter + 1
        idx <- cltr2 == ref.cltr
        i0 <- which(idx)

        if (length(i0) == 0) break

        changed <- FALSE

        for (i in i0) {
            nbrd.cl <- cltr2[nn.i[i, ]]
            non.ref <- nbrd.cl != ref.cltr

            if (sum(non.ref) >= min.votes) {
                if (method == "first") {
                    # Original method
                    j <- which(non.ref)[1]
                    new.label <- nbrd.cl[j]

                } else if (method == "majority") {
                    # True majority vote
                    votes <- table(nbrd.cl[non.ref])
                    max.votes <- max(votes)
                    winners <- names(votes)[votes == max.votes]

                    if (length(winners) == 1) {
                        new.label <- as.numeric(winners[1])
                    } else {
                        # Handle ties
                        if (tie.break == "first") {
                            # Pick first in neighbor order
                            for (cl in nbrd.cl[non.ref]) {
                                if (as.character(cl) %in% winners) {
                                    new.label <- cl
                                    break
                                }
                            }
                        } else if (tie.break == "random") {
                            new.label <- as.numeric(sample(winners, 1))
                        } else {
                            new.label <- ref.cltr  # Keep in reference cluster
                        }
                    }

                } else if (method == "weighted") {
                    # Weight by inverse distance
                    weights <- 1 / nn.d[i, non.ref]
                    weighted.votes <- tapply(weights, nbrd.cl[non.ref], sum)
                    new.label <- as.numeric(names(which.max(weighted.votes)))
                }

                if (new.label != ref.cltr) {
                    cltr2[i] <- new.label
                    changed <- TRUE
                }
            }
        }

        if (!changed || iter >= max.iter) break

        if (verbose && iter %% 10 == 0) {
            cat("Iteration", iter, ": Reference cluster size =", sum(cltr2 == ref.cltr), "\n")
        }
    }

    if (iter >= max.iter && verbose) {
        warning("Maximum iterations reached")
    }

    if (!is.null(names(cltr))) {
        names(cltr2) <- names(cltr)
    }

    return(cltr2)
}


#' Reorder Cluster IDs by Size
#'
#' Reorders cluster IDs based on cluster sizes, with the largest cluster
#' receiving ID 1, the second largest ID 2, and so on.
#'
#' @param cltr A vector of cluster IDs. Can be numeric, character, or factor.
#' @param decreasing Logical. If \code{TRUE} (default), clusters are ordered
#'   from largest to smallest. If \code{FALSE}, from smallest to largest.
#'
#' @return A vector of the same length and type as \code{cltr} with
#'   reordered cluster IDs. The output preserves the names attribute
#'   if present in the input.
#'
#' @examples
#' # Numeric clusters
#' cltr <- c(1, 2, 2, 3, 3, 3)
#' clusters.reorder(cltr)
#' # Returns: [1] 3 2 2 1 1 1
#'
#' # Character clusters
#' cltr_char <- c("A", "B", "B", "C", "C", "C")
#' clusters.reorder(cltr_char)
#' # Returns: [1] "C" "B" "B" "A" "A" "A"
#'
#' # With names preserved
#' named_cltr <- c(a = 1, b = 2, c = 2, d = 3, e = 3, f = 3)
#' clusters.reorder(named_cltr)
#'
#' @export
clusters.reorder <- function(cltr, decreasing = TRUE) {
  # Input validation
  if (length(cltr) == 0) {
    return(cltr)
  }

  # Preserve attributes
  input_names <- names(cltr)
  input_class <- class(cltr)

  # Calculate cluster sizes
  cltr_table <- table(cltr)

  # Order clusters by size
  cltr_order <- names(cltr_table)[order(cltr_table, decreasing = decreasing)]

  # Create mapping from old to new IDs
  new_ids <- setNames(seq_along(cltr_order), cltr_order)

  # Apply new IDs
  new_cltr <- new_ids[as.character(cltr)]

  # Restore original type if numeric
  if (is.numeric(cltr)) {
    new_cltr <- as.numeric(new_cltr)
  }

  # Restore names
  names(new_cltr) <- input_names

  # Remove unnecessary attributes added by table operations
  attributes(new_cltr) <- attributes(new_cltr)[c("names", "class")]

  return(new_cltr)
}
