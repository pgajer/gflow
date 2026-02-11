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

#' Graph kNN Imputation for Categorical Vertex Labels
#'
#' Imputes missing (NA) vertex labels using k-nearest neighbors defined by
#' shortest-path distances in a weighted, undirected graph.
#'
#' @param adj.list List of integer vectors. Each vector contains indices of vertices
#'        adjacent to the corresponding vertex. Indices should be 1-based.
#' @param weight.list List of numeric vectors. Each vector contains weights (edge lengths)
#'        corresponding to adjacencies in adj.list.
#' @param y A categorical variable defined on the vertices of the undirected graph
#'        given by (adj.list, weight.list). Missing values are coded as NA.
#' @param k Number of nearest labeled neighbors used for imputation (default is 10).
#'
#' @return A vector of labels with NA values imputed where possible. The output
#'         preserves the input type and names.
#'
#' @details
#' For each vertex with \code{y[i] = NA}, the function computes shortest-path distances
#' to all vertices, selects the \code{k} nearest vertices with non-NA labels (excluding
#' the vertex itself), and assigns the majority label among those neighbors. Ties are
#' broken by the nearest-neighbor order (smallest distance first). Imputation uses the
#' original non-NA labels only; newly imputed values are not reused. If fewer than \code{k}
#' labeled neighbors are reachable, all available labeled neighbors are used. If none are
#' reachable, the value remains NA.
#'
#' @examples
#' adj.list <- list(c(2, 3), c(1, 3), c(1, 2))
#' weight.list <- list(c(1, 1), c(1, 1), c(1, 1))
#' y <- c("A", NA, "B")
#' graph.cltr.imputation(adj.list, weight.list, y, k = 2)
#'
#' @importFrom igraph graph_from_edgelist distances E make_empty_graph
#'
#' @export
graph.cltr.imputation <- function(adj.list, weight.list, y, k = 10) {

    # Input validation
    if (!is.list(adj.list)) {
        stop("adj.list must be a list")
    }
    if (!is.list(weight.list)) {
        stop("weight.list must be a list")
    }

    n.vertices <- length(adj.list)

    if (length(weight.list) != n.vertices) {
        stop("adj.list and weight.list must have the same length")
    }
    if (length(y) != n.vertices) {
        stop("y must have the same length as adj.list")
    }

    if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
        stop("k must be a positive integer")
    }

    for (i in seq_len(n.vertices)) {
        if (length(adj.list[[i]]) != length(weight.list[[i]])) {
            stop(sprintf("Mismatch in lengths for vertex %d: adj.list has %d entries but weight.list has %d",
                         i, length(adj.list[[i]]), length(weight.list[[i]])))
        }

        if (length(adj.list[[i]]) > 0) {
            if (!is.numeric(adj.list[[i]])) {
                stop(sprintf("adj.list[[%d]] must be a numeric vector", i))
            }
            if (any(is.na(adj.list[[i]]))) {
                stop(sprintf("adj.list[[%d]] contains NA values", i))
            }
            if (any(adj.list[[i]] != floor(adj.list[[i]]))) {
                stop(sprintf("adj.list[[%d]] must contain only integer values", i))
            }
            if (any(adj.list[[i]] < 1) || any(adj.list[[i]] > n.vertices)) {
                stop(sprintf("adj.list[[%d]] contains invalid vertex indices", i))
            }
        }

        if (length(weight.list[[i]]) > 0) {
            if (!is.numeric(weight.list[[i]])) {
                stop(sprintf("weight.list[[%d]] must be a numeric vector", i))
            }
            if (any(is.na(weight.list[[i]]))) {
                stop(sprintf("weight.list[[%d]] contains NA values", i))
            }
            if (any(weight.list[[i]] <= 0)) {
                stop(sprintf("All weights in weight.list[[%d]] must be positive", i))
            }
        }
    }

    # Quick return if there is nothing to impute
    na.idx <- which(is.na(y))
    if (length(na.idx) == 0) {
        return(y)
    }

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required but not installed. Please install it using install.packages('igraph')")
    }

    # Build igraph object
    graph.obj <- convert.adjacency.to.edge.matrix(adj.list, weight.list)
    if (nrow(graph.obj$edge.matrix) == 0) {
        g <- igraph::make_empty_graph(n = n.vertices, directed = FALSE)
    } else {
        g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
        igraph::E(g)$weight <- graph.obj$weights
    }

    y2 <- y
    y.obs <- y

    # If all labels are missing, nothing can be imputed
    if (all(is.na(y.obs))) {
        warning("All y values are NA; nothing to impute.")
        return(y)
    }

    # Compute distances from NA vertices to all vertices
    dist.mat <- igraph::distances(
        g,
        v = na.idx,
        to = seq_len(n.vertices),
        weights = igraph::E(g)$weight
    )

    for (row.idx in seq_along(na.idx)) {
        i <- na.idx[row.idx]
        d <- dist.mat[row.idx, ]

        # Exclude the vertex itself
        d[i] <- Inf

        ord <- order(d, na.last = NA)
        if (length(ord) == 0) {
            next
        }

        # Keep only finite distances
        ord <- ord[is.finite(d[ord])]
        if (length(ord) == 0) {
            next
        }

        # Keep only neighbors with non-NA labels from the observed data
        ord <- ord[!is.na(y.obs[ord])]
        if (length(ord) == 0) {
            next
        }

        nbr.idx <- ord[seq_len(min(k, length(ord)))]
        nbr.vals <- y.obs[nbr.idx]
        nbr.vals.chr <- as.character(nbr.vals)

        votes <- table(nbr.vals.chr)
        max.votes <- max(votes)
        winners <- names(votes)[votes == max.votes]

        # Tie-break by nearest-neighbor order
        pick.idx <- which(nbr.vals.chr %in% winners)[1]
        y2[i] <- nbr.vals[pick.idx]
    }

    # Preserve names if they existed
    if (!is.null(names(y))) {
        names(y2) <- names(y)
    }

    return(y2)
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


#' Select Number of Clusters for a Hierarchical Clustering
#'
#' Evaluate candidate cuts (numbers of clusters) for a given \code{hclust} object
#' and select the best cut using internal cluster validity criteria.
#'
#' This routine is designed to support (i) a focused "core" panel of indices that
#' are commonly used for k-selection, and (ii) an optional extended diagnostics
#' panel obtained from \code{fpc::cluster.stats()}.
#'
#' @param hc A \code{hclust} object, typically from \code{hclust()}.
#' @param X Optional numeric matrix (rows = observations). Required for computing
#'   \code{connectivity} via \code{clValid::connectivity()} and for computing
#'   \code{d} if \code{d} is not provided.
#' @param d Optional \code{dist} object. Used for Dunn/silhouette and for
#'   \code{fpc::cluster.stats()}.
#' @param method Criterion used to select the "best" number of clusters. One of
#'   \code{"dunn"}, \code{"silhouette"}, \code{"ch"}, \code{"wb.ratio"},
#'   \code{"connectivity"}, or \code{"all"}.
#' @param k Integer vector of candidate numbers of clusters. If \code{NULL},
#'   defaults to \code{2:min(50, n-1)} where \code{n} is the number of observations.
#' @param neighb.size Integer kNN size used for \code{clValid::connectivity()}.
#'   Default is 10.
#' @param stats.level Controls how many indices are returned in \code{scores}:
#'   \code{"core"}, \code{"extended"}, or \code{"all"}.
#'   \code{"core"} returns a focused set of k-selection indices: \code{"dunn"},
#'   \code{"avg.silwidth"}, \code{"ch"}, and \code{"wb.ratio"}.
#'   \code{"extended"} returns \code{"core"} plus additional diagnostics:
#'   \code{"pearsongamma"}, \code{"entropy"}, \code{"widestgap"}, and \code{"sindex"}.
#'   \code{"all"} attempts to return all scalar indices produced by
#'   \code{fpc::cluster.stats()} (subject to \code{include.slow}); note that this
#'   can substantially increase computation time for large \code{n}, and some
#'   returned components are descriptive diagnostics rather than k-selection criteria.
#'   When \code{stats.level="all"}, indices known to be slow (e.g., \code{"g2"},
#'   \code{"g3"}) are omitted unless \code{include.slow=TRUE}. Default is \code{"core"}.
#' @param with.fpc.stats Logical; if \code{TRUE} and \code{fpc} is available,
#'   compute indices via \code{fpc::cluster.stats()}. Default is \code{TRUE}.
#' @param with.connectivity Logical; if \code{TRUE}, compute connectivity when
#'   possible (requires \code{X} and \code{clValid}). Default is \code{TRUE}.
#' @param include.slow Logical; if \code{TRUE} and \code{stats.level="all"},
#'   attempt to include indices that can be slow for large \code{n} (e.g., G2/G3).
#'   Default is \code{FALSE}.
#' @param n.cores Integer; if > 1 and \code{foreach/doParallel} are available,
#'   attempt parallel evaluation across \code{k}. Default is 1.
#' @param verbose Logical; print progress. Default is \code{FALSE}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{scores}: data.frame with one row per candidate \code{k} and
#'         columns for each computed index.
#'   \item \code{opt.k}: named integer vector of selected \code{k} values (per metric).
#'   \item \code{best}: list of best cluster assignments per metric (per \code{opt.k}).
#'   \item \code{clusters}: list of cluster assignments for each \code{k}.
#'   \item \code{directions}: named character vector indicating optimization
#'         direction (maximize/minimize) for indices used in automated selection.
#' }
#'
#' @references
#' \url{https://cran.r-project.org/package=fpc}
#' \url{https://cran.r-project.org/package=cluster}
#' \url{https://cran.r-project.org/package=clValid}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- rbind(matrix(rnorm(200), ncol = 2),
#'            matrix(rnorm(200, mean = 3), ncol = 2))
#' d <- dist(X)
#' hc <- hclust(d, method = "ward.D2")
#'
#' res <- hclust.select.k(hc, X = X, d = d,
#'                        method = "all",
#'                        stats.level = "extended",
#'                        with.connectivity = TRUE,
#'                        verbose = TRUE)
#'
#' res$opt.k
#' head(res$scores)
#' table(res$best$silhouette)
#' }
#'
#' @export
hclust.select.k <- function(hc,
                            X = NULL,
                            d = NULL,
                            method = c("dunn", "silhouette", "ch", "wb.ratio",
                                       "connectivity", "all"),
                            k = NULL,
                            neighb.size = 10L,
                            stats.level = c("core", "extended", "all"),
                            with.fpc.stats = TRUE,
                            with.connectivity = TRUE,
                            include.slow = FALSE,
                            n.cores = 1L,
                            verbose = FALSE) {

    ## ---------------------------------------------------------------------
    ## Helpers
    ## ---------------------------------------------------------------------

    null.coalesce <- function(a, b) {
        if (!is.null(a)) a else b
    }

    is.scalar.int <- function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x)
    }

    take.numeric <- function(x) {
        if (is.null(x)) return(NA_real_)
        if (length(x) == 0L) return(NA_real_)
        as.numeric(x[1L])
    }

    ## Known optimization directions for automated selection
    ## (Only include indices with stable interpretation across use cases.)
    directions <- c(
        dunn = "maximize",
        silhouette = "maximize",
        ch = "maximize",
        wb.ratio = "minimize",
        connectivity = "minimize"
    )

    ## Map selection metric name to score column name
    method <- match.arg(method)
    stats.level <- match.arg(stats.level)

    call <- match.call()

    args <- list(
        method = method,
        k = k,
        neighb.size = neighb.size,
        stats.level = stats.level,
        with.fpc.stats = with.fpc.stats,
        with.connectivity = with.connectivity,
        include.slow = include.slow,
        n.cores = n.cores,
        verbose = verbose
    )

    ## Tracking what data inputs were actually provided
    args$has.X <- !is.null(X)
    args$has.d <- !is.null(d)

    ## ---------------------------------------------------------------------
    ## Input validation
    ## ---------------------------------------------------------------------

    if (!inherits(hc, "hclust")) {
        stop("hc must be a 'hclust' object.")
    }

    if (!is.scalar.int(neighb.size) || neighb.size < 1L) {
        stop("neighb.size must be a positive integer.")
    }
    neighb.size <- as.integer(neighb.size)

    if (!is.scalar.int(n.cores) || n.cores < 1L) {
        stop("n.cores must be a positive integer.")
    }
    n.cores <- as.integer(n.cores)

    ## Infer n
    n <- length(hc$order)
    if (!is.null(X)) {
        if (!is.matrix(X)) X <- as.matrix(X)
        if (!is.numeric(X)) stop("X must be numeric.")
        if (nrow(X) != n) {
            stop("nrow(X) does not match the number of observations in hc.")
        }
    }

    if (is.null(k)) {
        k.max <- min(50L, n - 1L)
        if (k.max < 2L) stop("Need at least 2 observations.")
        k <- 2L:k.max
    } else {
        k <- as.integer(k)
        k <- k[is.finite(k)]
        k <- sort(unique(k))
        k <- k[k >= 2L & k <= (n - 1L)]
        if (length(k) == 0L) stop("No valid k values after filtering (need 2 <= k <= n-1).")
    }

    ## Distances
    if (!is.null(d) && !inherits(d, "dist")) {
        stop("d must be a 'dist' object.")
    }
    if (is.null(d)) {
        if (is.null(X)) {
            stop("Provide either d (dist) or X (to compute dist(X)).")
        }
        d <- stats::dist(X)
    }

    ## Suggested-package availability
    has.fpc <- requireNamespace("fpc", quietly = TRUE)
    has.cluster <- requireNamespace("cluster", quietly = TRUE)
    has.clValid <- requireNamespace("clValid", quietly = TRUE)
    has.foreach <- requireNamespace("foreach", quietly = TRUE)
    has.doParallel <- requireNamespace("doParallel", quietly = TRUE)

    if (with.fpc.stats && !has.fpc) {
        with.fpc.stats <- FALSE
        if (verbose) cat("Note: package 'fpc' not available; skipping fpc::cluster.stats().\n")
    }
    if (with.connectivity && !has.clValid) {
        with.connectivity <- FALSE
        if (verbose) cat("Note: package 'clValid' not available; skipping connectivity.\n")
    }
    if (with.connectivity && is.null(X)) {
        with.connectivity <- FALSE
        if (verbose) cat("Note: X not provided; skipping connectivity.\n")
    }

    ## If silhouette is requested but cluster not installed, we can still use
    ## fpc::cluster.stats() when available (it provides avg.silwidth).
    if ((method %in% c("silhouette", "all")) && !has.cluster && !with.fpc.stats) {
        stop("Selecting by silhouette requires either package 'cluster' or package 'fpc'.")
    }

    ## ---------------------------------------------------------------------
    ## Cluster assignments for each k (computed once)
    ## ---------------------------------------------------------------------

    clusters <- vector("list", length(k))
    names(clusters) <- as.character(k)

    for (i in seq_along(k)) {
        clusters[[i]] <- stats::cutree(hc, k = k[i])
    }

    ## ---------------------------------------------------------------------
    ## Decide which fpc stats to keep
    ## ---------------------------------------------------------------------

    core.keep <- c("dunn", "avg.silwidth", "ch", "wb.ratio")

    ## Extended: add a small number of additional, often-useful diagnostics
    extended.keep <- unique(c(
        core.keep,
        "pearsongamma",
        "entropy",
        "widestgap",
        "sindex"
    ))

    ## All: try to include everything that cluster.stats returns, but optionally
    ## filter out known slow ones unless include.slow = TRUE.
    ## Note: We do not attempt to compute indices requiring alt.clustering.
    slow.names <- c("g2", "g3")

    keep.names <- switch(stats.level,
        core = core.keep,
        extended = extended.keep,
        all = NULL
    )

    ## ---------------------------------------------------------------------
    ## Worker to compute metrics for one k
    ## ---------------------------------------------------------------------

    compute.one <- function(ki, cl) {

        out <- list()
        out$k <- ki

        ## fpc stats (distance-based)
        fpc.stats <- NULL
        if (with.fpc.stats) {
            fpc.stats <- tryCatch(
                fpc::cluster.stats(d = d, clustering = cl),
                error = function(e) NULL
            )
        }

        ## Dunn
        ## Prefer fpc when available; otherwise try clValid::dunn on distance matrix.
        out$dunn <- NA_real_
        if (!is.null(fpc.stats) && "dunn" %in% names(fpc.stats)) {
            out$dunn <- take.numeric(fpc.stats$dunn)
        } else if (has.clValid) {
            out$dunn <- tryCatch(
                clValid::dunn(as.matrix(d), clusters = cl),
                error = function(e) NA_real_
            )
        }

        ## Silhouette (mean)
        ## Prefer cluster::silhouette if available; otherwise fall back to fpc avg.silwidth.
        out$silhouette <- NA_real_
        if (has.cluster) {
            out$silhouette <- tryCatch({
                si <- cluster::silhouette(cl, d)
                mean(si[, "sil_width"])
            }, error = function(e) NA_real_)
        } else if (!is.null(fpc.stats) && "avg.silwidth" %in% names(fpc.stats)) {
            out$silhouette <- take.numeric(fpc.stats$avg.silwidth)
        }

        ## Calinski-Harabasz
        out$ch <- NA_real_
        if (!is.null(fpc.stats) && "ch" %in% names(fpc.stats)) {
            out$ch <- take.numeric(fpc.stats$ch)
        }

        ## W/B ratio
        out$wb.ratio <- NA_real_
        if (!is.null(fpc.stats) && "wb.ratio" %in% names(fpc.stats)) {
            out$wb.ratio <- take.numeric(fpc.stats$wb.ratio)
        }

        ## Connectivity (requires X)
        out$connectivity <- NA_real_
        if (with.connectivity) {
            out$connectivity <- tryCatch(
                clValid::connectivity(clusters = cl, Data = X, neighbSize = neighb.size),
                error = function(e) NA_real_
            )
        }

        ## Optional extra fpc stats into out$fpc (for table assembly)
        out$fpc <- NULL
        if (!is.null(fpc.stats)) {
            ## Keep either specified subset or everything (minus slow unless requested).
            if (!is.null(keep.names)) {
                keep <- intersect(keep.names, names(fpc.stats))
            } else {
                keep <- names(fpc.stats)
                if (!include.slow) keep <- setdiff(keep, slow.names)
            }

            ## Convert to a named numeric vector when possible
            v <- suppressWarnings(as.numeric(fpc.stats[keep]))
            names(v) <- keep

            ## Some components are not numeric scalars; coerce failures to NA
            bad <- !is.finite(v)
            v[bad] <- NA_real_

            out$fpc <- v
        }

        out
    }

    ## ---------------------------------------------------------------------
    ## Evaluate across k (optionally parallel)
    ## ---------------------------------------------------------------------

    can.parallel <- (n.cores > 1L && has.foreach && has.doParallel)
    if (can.parallel) {
        doParallel::registerDoParallel(n.cores)
        on.exit(doParallel::stopImplicitCluster(), add = TRUE)

        res <- foreach::foreach(i = seq_along(k), .inorder = TRUE) %dopar% {
            compute.one(k[i], clusters[[i]])
        }
    } else {
        res <- vector("list", length(k))
        for (i in seq_along(k)) {
            if (verbose) cat(sprintf("\rEvaluating k = %d", k[i]))
            res[[i]] <- compute.one(k[i], clusters[[i]])
        }
        if (verbose) cat("\n")
    }

    ## ---------------------------------------------------------------------
    ## Assemble score table
    ## ---------------------------------------------------------------------

    base.scores <- data.frame(
        k = vapply(res, function(r) r$k, integer(1L)),
        dunn = vapply(res, function(r) r$dunn, numeric(1L)),
        silhouette = vapply(res, function(r) r$silhouette, numeric(1L)),
        ch = vapply(res, function(r) r$ch, numeric(1L)),
        wb.ratio = vapply(res, function(r) r$wb.ratio, numeric(1L)),
        connectivity = vapply(res, function(r) r$connectivity, numeric(1L)),
        stringsAsFactors = FALSE
    )

    ## Add fpc-derived columns according to stats.level
    if (stats.level != "core" && with.fpc.stats) {
        ## Collect union of fpc names across k
        fpc.names <- unique(unlist(lapply(res, function(r) names(r$fpc)), use.names = FALSE))
        fpc.names <- setdiff(fpc.names, c("dunn", "avg.silwidth", "ch", "wb.ratio"))

        if (length(fpc.names) > 0L) {
            fpc.mat <- matrix(NA_real_, nrow = nrow(base.scores), ncol = length(fpc.names))
            colnames(fpc.mat) <- fpc.names

            for (i in seq_along(res)) {
                v <- res[[i]]$fpc
                if (!is.null(v)) {
                    common <- intersect(names(v), fpc.names)
                    if (length(common) > 0L) {
                        fpc.mat[i, common] <- v[common]
                    }
                }
            }

            scores <- cbind(base.scores, as.data.frame(fpc.mat, stringsAsFactors = FALSE))
        } else {
            scores <- base.scores
        }
    } else {
        scores <- base.scores
    }

    ## ---------------------------------------------------------------------
    ## Select optima and prepare "best" clusterings
    ## ---------------------------------------------------------------------

    pick.best.k <- function(score.vec, direction) {
        ok <- is.finite(score.vec)
        if (!any(ok)) return(NA_integer_)
        if (direction == "maximize") {
            as.integer(scores$k[which.max(ifelse(ok, score.vec, -Inf))])
        } else if (direction == "minimize") {
            as.integer(scores$k[which.min(ifelse(ok, score.vec, Inf))])
        } else {
            NA_integer_
        }
    }

    opt.k <- integer(0)
    best <- list()

    add.best <- function(metric.name, score.col) {
        dir <- directions[[metric.name]]
        kk <- pick.best.k(scores[[score.col]], dir)
        if (!is.na(kk)) {
            opt.k[[metric.name]] <<- kk
            best[[metric.name]] <<- clusters[[as.character(kk)]]
        }
    }

    if (method == "all") {
        add.best("dunn", "dunn")
        add.best("silhouette", "silhouette")
        add.best("ch", "ch")
        add.best("wb.ratio", "wb.ratio")
        if (with.connectivity) add.best("connectivity", "connectivity")
    } else {
        if (method == "dunn") add.best("dunn", "dunn")
        if (method == "silhouette") add.best("silhouette", "silhouette")
        if (method == "ch") add.best("ch", "ch")
        if (method == "wb.ratio") add.best("wb.ratio", "wb.ratio")
        if (method == "connectivity") {
            if (!with.connectivity) {
                stop("connectivity selection requested, but connectivity could not be computed (need X and clValid).")
            }
            add.best("connectivity", "connectivity")
        }
    }

    ## ---------------------------------------------------------------------
    ## Return
    ## ---------------------------------------------------------------------

    out <- list(
        scores = scores,
        opt.k = opt.k,
        best = best,
        clusters = clusters,
        directions = directions,
        call = call,
        args = args
    )

    class(out) <- c("hclust_select_k", "list")
    out
}

#' Print Method for hclust_select_k Objects
#'
#' Prints a concise report for objects returned by \code{hclust.select.k()},
#' including available metrics, selected \code{k} values, and a brief
#' disagreement summary across metrics (when applicable).
#'
#' @param x An object of class \code{"hclust_select_k"}, typically returned by
#'   \code{hclust.select.k()}.
#' @param ... Additional arguments (ignored).
#' @param top.n Integer; number of top-scoring \code{k} values to display per
#'   selection-grade metric (when available). Default is 5.
#' @param show.call Logical; if TRUE, prints function call.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.hclust_select_k <- function(x, ..., top.n = 0L, show.call = FALSE) {

    if (is.null(x$scores) || !is.data.frame(x$scores)) {
        stop("Invalid object: missing scores.")
    }

    scores <- x$scores
    k <- scores$k

    cat("hclust.select.k result\n")
    cat("======================\n")

    if (!is.null(x$call) && show.call) {
        cat("Call:\n")
        cat("  ")
        dput(x$call)
    }

    cat(sprintf("Candidate k: %d to %d (n = %d)\n",
                min(k), max(k), length(k)))

    ## Report which metrics are present (non-NA)
    metric.cols <- setdiff(names(scores), "k")
    avail <- vapply(metric.cols, function(nm) any(is.finite(scores[[nm]])), logical(1L))
    avail.metrics <- metric.cols[avail]

    cat("Available metrics:\n")
    if (length(avail.metrics) == 0L) {
        cat("  (none)\n")
    } else {
        cat("  ", paste(avail.metrics, collapse = ", "), "\n", sep = "")
    }

    ## opt.k summary
    if (!is.null(x$opt.k) && length(x$opt.k) > 0L) {
        cat("Selected k (opt.k):\n")
        for (nm in names(x$opt.k)) {
            cat(sprintf("  %-13s %d\n", nm, x$opt.k[[nm]]))
        }

        ## Disagreement summary across metrics
        opt.k.vec <- as.integer(unlist(x$opt.k, use.names = FALSE))
        opt.k.names <- names(x$opt.k)
        ok <- is.finite(opt.k.vec)

        if (sum(ok) >= 2L) {
            opt.k.ok <- opt.k.vec[ok]
            rng <- range(opt.k.ok)
            cat("Disagreement summary (across selected metrics):\n")
            cat(sprintf("  range: %d to %d (spread = %d)\n",
                        rng[1L], rng[2L], rng[2L] - rng[1L]))
            cat(sprintf("  mean: %.3f; sd: %.3f\n",
                        mean(opt.k.ok), stats::sd(opt.k.ok)))

            ## Optional: show which metric picked min/max k
            min.idx <- which.min(opt.k.ok)
            max.idx <- which.max(opt.k.ok)
            ## Map back to names (handle ties by first)
            opt.k.ok.names <- opt.k.names[ok]
            cat(sprintf("  min k chosen by: %s\n", opt.k.ok.names[min.idx]))
            cat(sprintf("  max k chosen by: %s\n", opt.k.ok.names[max.idx]))
        } else {
            cat("Disagreement summary: not available (fewer than 2 selected metrics).\n")
        }

    } else {
        cat("Selected k (opt.k): (none)\n")
    }

    ## Show top rows per selection-grade metric if present
    top.n <- as.integer(top.n)
    if (is.finite(top.n) && top.n > 0L) {

        ## Directions are stored in object; fall back to defaults
        directions <- x$directions
        if (is.null(directions)) {
            directions <- c(dunn = "maximize",
                            silhouette = "maximize",
                            ch = "maximize",
                            wb.ratio = "minimize",
                            connectivity = "minimize")
        }

        show.metrics <- intersect(names(directions), names(scores))
        show.metrics <- show.metrics[vapply(show.metrics, function(nm) any(is.finite(scores[[nm]])), logical(1L))]

        if (length(show.metrics) > 0L) {
            cat("\nTop k by metric:\n")
            for (m in show.metrics) {
                v <- scores[[m]]
                ok <- is.finite(v)
                if (!any(ok)) next

                ord <- if (directions[[m]] == "maximize") {
                    order(v, decreasing = TRUE, na.last = NA)
                } else {
                    order(v, decreasing = FALSE, na.last = NA)
                }
                ord <- head(ord, top.n)

                cat(sprintf("  %s (%s):\n", m, directions[[m]]))
                for (ii in ord) {
                    cat(sprintf("    k=%d  %s=%.6g\n", scores$k[ii], m, scores[[m]][ii]))
                }
            }
        }
    }

    invisible(x)
}


#' Plot Method for hclust_select_k Objects
#'
#' Produces base R diagnostic plots of cluster validity indices (scores) versus
#' candidate numbers of clusters \code{k}.
#'
#' @param x An object of class \code{"hclust_select_k"}, typically returned by
#'   \code{hclust.select.k()}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @param metric Character; which metric to plot. One of
#'   \code{"dunn"}, \code{"silhouette"}, \code{"ch"}, \code{"wb.ratio"},
#'   \code{"connectivity"}, or \code{"all"}. Default is \code{"dunn"}.
#' @param show.opt Logical; if TRUE, draw a vertical line at the selected
#'   \code{k} (if available) for the plotted metric(s). Default is TRUE.
#' @param pch Plotting character used for points. Default is 16.
#' @param lty Line type used for connecting lines. Default is 1.
#' @param main Optional plot title. If NULL, a default title is used.
#' @param xlab X-axis label. Default is \code{"k"}.
#' @param ylab Y-axis label. If NULL, defaults to the metric name.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
plot.hclust_select_k <- function(x,
                                 ...,
                                 metric = c("dunn", "silhouette", "ch", "wb.ratio", "connectivity", "all"),
                                 show.opt = TRUE,
                                 pch = 16,
                                 lty = 1,
                                 main = NULL,
                                 xlab = "k",
                                 ylab = NULL) {

    if (is.null(x$scores) || !is.data.frame(x$scores)) {
        stop("Invalid object: missing scores.")
    }

    scores <- x$scores
    metric <- match.arg(metric)

    directions <- x$directions
    if (is.null(directions)) {
        directions <- c(dunn = "maximize",
                        silhouette = "maximize",
                        ch = "maximize",
                        wb.ratio = "minimize",
                        connectivity = "minimize")
    }

    plot.one <- function(m) {

        if (!m %in% names(scores)) {
            stop(sprintf("Metric '%s' is not present in scores.", m))
        }

        v <- scores[[m]]
        ok <- is.finite(v)

        if (!any(ok)) {
            plot.new()
            title(main = sprintf("%s (not available)", m))
            return(invisible(NULL))
        }

        yy <- v
        yy.lab <- ylab %||% m

        mm.main <- if (is.null(main)) {
            sprintf("%s (%s)", m, directions[[m]] %||% "n/a")
        } else {
            main
        }

        plot(scores$k[ok], yy[ok],
             type = "b",
             pch = pch,
             lty = lty,
             xlab = xlab,
             ylab = yy.lab,
             main = mm.main,
             ...)

        if (isTRUE(show.opt) && !is.null(x$opt.k) && m %in% names(x$opt.k)) {
            kk <- x$opt.k[[m]]
            abline(v = kk, lty = 2)
        }

        invisible(NULL)
    }

    if (metric == "all") {

        ## Plot only selection-grade metrics that exist and have finite values
        plot.metrics <- intersect(names(directions), names(scores))
        plot.metrics <- plot.metrics[vapply(plot.metrics, function(nm) any(is.finite(scores[[nm]])), logical(1L))]

        if (length(plot.metrics) == 0L) {
            stop("No plottable metrics found in scores.")
        }

        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par), add = TRUE)

        n.panels <- length(plot.metrics)
        if (n.panels <= 4L) {
            par(mfrow = c(2L, 2L))
        } else {
            par(mfrow = c(3L, 2L))
        }

        for (m in plot.metrics) plot.one(m)

    } else {
        plot.one(metric)
    }

    invisible(x)
}

#' Describe an Object
#'
#' S3 generic that returns a compact, class-specific description intended to
#' assist interpretation of an object and its key diagnostics.
#'
#' @param x An R object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A class-specific object (often a \code{data.frame}) describing \code{x}.
#'
#' @export
describe <- function(x, ...) {
    UseMethod("describe")
}

#' Describe a hclust_select_k Object
#'
#' Returns a compact reference describing cluster validity metrics available in
#' a \code{"hclust_select_k"} object, including optimization direction and brief
#' interpretation notes.
#'
#' @param x An object of class \code{"hclust_select_k"}, typically returned by
#'   \code{hclust.select.k()}.
#' @param ... Additional arguments (ignored).
#' @param include.unavailable Logical; if TRUE, include all known metrics in the
#'   output even if not computed for \code{x}. Default is FALSE (show only
#'   computed metrics).
#' @param metric Optional character; if provided, restrict output to a single
#'   metric (e.g., \code{"silhouette"}). If \code{NULL}, include all selected metrics.
#' @param format Character; one of \code{"table"} or \code{"cards"}. \code{"table"}
#'   returns a data.frame (default). \code{"cards"} prints one metric at a time
#'   and invisibly returns the same data.frame.
#'
#' @return A \code{data.frame} with columns \code{metric}, \code{direction},
#'   \code{what.it.measures}, \code{pros}, \code{cons}. If \code{format="cards"},
#'   the table is returned invisibly.
#'
#' @export
describe.hclust_select_k <- function(x,
                                     ...,
                                     include.unavailable = FALSE,
                                     metric = NULL,
                                     format = c("cards", "table")) {

    if (is.null(x$scores) || !is.data.frame(x$scores)) {
        stop("Invalid object: missing scores.")
    }

    format <- match.arg(format)

    tab <- data.frame(
        metric = c("silhouette", "ch", "dunn", "wb.ratio", "connectivity"),
        direction = c("maximize", "maximize", "maximize", "minimize", "minimize"),
        what.it.measures = c(
            "Mean silhouette width: separation vs within-cluster cohesion (distance-based).",
            "Calinski-Harabasz: between-cluster dispersion relative to within-cluster dispersion.",
            "Dunn index: min inter-cluster separation / max intra-cluster diameter.",
            "Within/between dispersion ratio (smaller indicates better separation relative to spread).",
            "Neighbor consistency: penalizes nearby points assigned to different clusters (kNN-based)."
        ),
        pros = c(
            "Often a good general-purpose internal check; penalizes over-splitting when clusters become too similar.",
            "Fast and often stable; often selects a smaller number of clusters when the data contain a few well-separated groups.",
            "Rewards compact, well-separated clusters; can detect clear separation.",
            "Useful as a compactness/separation diagnostic; sometimes aligns with 'spherical cluster' intuition.",
            "Captures local structure; useful for manifold-like data where local neighborhoods matter."
        ),
        cons = c(
            "Can favor convex / evenly sized clusters; may be conservative under chaining or strong density gradients.",
            "Can favor small k in some geometries; sensitive to distance scaling and to elongated clusters.",
            "Sensitive to outliers/singletons because it uses extreme distances; can push to large k to shrink diameters.",
            "Can monotonically improve with larger k in some cases (risk of over-splitting); interpret with other metrics.",
            "Depends on neighb.size and requires X; can disagree with distance-based compactness indices."
        ),
        stringsAsFactors = FALSE
    )

    ## Filter to computed metrics unless include.unavailable=TRUE
    if (!isTRUE(include.unavailable)) {
        scores <- x$scores
        keep <- tab$metric %in% names(scores)
        keep <- keep & vapply(tab$metric, function(m) any(is.finite(scores[[m]])), logical(1L))
        tab <- tab[keep, , drop = FALSE]
    }

    ## Restrict to one metric if requested
    if (!is.null(metric)) {
        metric <- as.character(metric)[1L]
        if (!metric %in% tab$metric) {
            stop(sprintf("Unknown or unavailable metric '%s'. Available: %s",
                         metric, paste(tab$metric, collapse = ", ")))
        }
        tab <- tab[tab$metric == metric, , drop = FALSE]
    }

    if (format == "table") {
        return(tab)
    }

    ## format == "cards"
    for (i in seq_len(nrow(tab))) {
        cat(sprintf("%s (%s)\n", tab$metric[i], tab$direction[i]))
        cat(strrep("-", nchar(tab$metric[i]) + nchar(tab$direction[i]) + 3L), "\n", sep = "")
        cat("Measures: ", tab$what.it.measures[i], "\n", sep = "")
        cat("Pros:     ", tab$pros[i], "\n", sep = "")
        cat("Cons:     ", tab$cons[i], "\n", sep = "")
        if (i < nrow(tab)) cat("\n")
    }

    ## -------------------------------------------------------------------------
    ## Interpretation note: metric disagreement patterns (quantile-based grouping)
    ## -------------------------------------------------------------------------

    if (!is.null(x$opt.k) && length(x$opt.k) > 1L) {

        opt <- x$opt.k

        ## Heuristic grouping: "small-k leaning" vs "large-k leaning" metrics
        ## (based on typical behavior in practice, not a theorem)
        small.k.metrics <- intersect(c("silhouette", "ch"), names(opt))
        large.k.metrics <- intersect(c("dunn", "wb.ratio", "connectivity"), names(opt))

        small.k <- opt[small.k.metrics]
        large.k <- opt[large.k.metrics]

        if (length(small.k) > 0L && length(large.k) > 0L) {

            ## Trigger note when the groups are clearly separated
            if (max(small.k) <= min(large.k)) {

                fmt <- function(v) {
                    paste(sprintf("%s=%d", names(v), as.integer(v)), collapse = ", ")
                }

                cat("\nInterpretation note\n")
                cat("-------------------\n")
                cat("Some metrics favor small k: ", fmt(small.k), "\n", sep = "")
                cat("Other metrics favor large k: ", fmt(large.k), "\n", sep = "")
                cat("This pattern can indicate chaining/gradual structure, elongated clusters,\n")
                cat("outlier sensitivity, or a distance scaling effect. Consider inspecting the\n")
                cat("dendrogram and cluster size distribution at candidate k values.\n")
            }
        }
    }


    if (!is.null(x$opt.k) && length(x$opt.k) > 1L && !is.null(x$scores)) {

        scores <- x$scores
        k.min <- min(scores$k)
        k.max <- max(scores$k)

        opt <- x$opt.k
        opt.k <- as.integer(unlist(opt, use.names = TRUE))
        opt.k <- opt.k[is.finite(opt.k)]

        if (length(opt.k) >= 2L) {

            ## Define small/large groups by quartiles of selected k values
            q1 <- as.numeric(stats::quantile(opt.k, probs = 0.25, type = 7, names = FALSE))
            q3 <- as.numeric(stats::quantile(opt.k, probs = 0.75, type = 7, names = FALSE))

            small.idx <- opt.k <= q1
            large.idx <- opt.k >= q3

            small <- opt.k[small.idx]
            large <- opt.k[large.idx]

            ## Only print if there is a meaningful separation and both groups exist
            if (length(small) > 0L && length(large) > 0L) {

                ## Explicit separation check (avoid printing for mild spread)
                if (max(small) < min(large)) {

                    fmt <- function(v) {
                        v <- v[order(v)]
                        paste(sprintf("%s=%d", names(v), as.integer(v)), collapse = ", ")
                    }

                    cat("\nInterpretation note\n")
                    cat("-------------------\n")
                    cat(sprintf("Selected k values span %d to %d (k.min=%d, k.max=%d).\n",
                                min(opt.k), max(opt.k), k.min, k.max))
                    cat(sprintf("Small-k group: k <= %d: %s\n",
                                as.integer(q1), fmt(small)))
                    cat(sprintf("Large-k group: k >= %d: %s\n",
                                as.integer(q3), fmt(large)))

                    cat("This split often indicates either (i) gradual/chaining structure in the dendrogram,\n")
                    cat("(ii) elongated clusters being split into smaller pieces, (iii) outlier sensitivity,\n")
                    cat("or (iv) distance scaling effects.\n")

                    ## Boundary-driven recommendation
                    at.max <- names(opt.k)[opt.k == k.max]
                    at.min <- names(opt.k)[opt.k == k.min]

                    if (length(at.max) > 0L) {
                        cat(sprintf("Boundary flag: %s selected k.max=%d.\n",
                                    paste(at.max, collapse = ", "), k.max))
                        cat("Recommendation: treat this as a potential over-splitting signature.\n")
                        cat("Consider (a) reducing k.max (e.g., set k to 2:20), (b) inspecting cluster sizes\n")
                        cat("for singletons/tiny clusters, and (c) checking whether the corresponding index\n")
                        cat("improves almost monotonically with k.\n")
                    }

                    if (length(at.min) > 0L) {
                        cat(sprintf("Boundary flag: %s selected k.min=%d.\n",
                                    paste(at.min, collapse = ", "), k.min))
                        cat("Recommendation: verify whether the data may genuinely have very few clusters,\n")
                        cat("or whether the distance definition is compressing structure (e.g., too aggressive\n")
                        cat("normalization / scaling).\n")
                    }
                }
            }
        }
    }

    invisible(tab)
}

#' Summary Method for hclust_select_k Objects
#'
#' Returns a ranked shortlist of candidate \code{k} values using rank aggregation
#' across computed selection-grade metrics. Optionally returns a table for all
#' candidate \code{k} values.
#'
#' Metrics that are not computed (i.e., have no finite values in
#' \code{object$scores}) are automatically excluded from aggregation and omitted
#' from the returned columns.
#'
#' @param object An object of class \code{"hclust_select_k"}, typically returned
#'   by \code{hclust.select.k()}.
#' @param ... Additional arguments (ignored).
#' @param metrics Character vector of metric column names to include in the rank
#'   aggregation. If \code{NULL}, uses selection-grade metrics found in
#'   \code{object$directions} that are present and computed in \code{object$scores}.
#' @param top.n Integer; number of \code{k} values to return in the shortlist
#'   (before adding forced \code{opt.k} rows). Default is 10.
#' @param min.metrics Integer; require at least this many metrics to contribute
#'   (finite ranks) for a given \code{k} to be eligible. Default is 1.
#' @param ties.method Character; passed to \code{rank()} (default \code{"average"}).
#' @param include.opt.k Logical; if TRUE, force inclusion of \code{k} values in
#'   \code{object$opt.k} in the returned shortlist (when present). Default is TRUE.
#' @param return Character; one of \code{"shortlist"}, \code{"all"}, or \code{"both"}.
#'   Default is \code{"shortlist"}.
#'
#' @return If \code{return="shortlist"}, a data.frame shortlist. If
#'   \code{return="all"}, a data.frame with one row per candidate \code{k}. If
#'   \code{return="both"}, a list with components \code{shortlist} and \code{all}.
#'
#' @export
summary.hclust_select_k <- function(object,
                                    ...,
                                    metrics = NULL,
                                    top.n = 10L,
                                    min.metrics = 1L,
                                    ties.method = "average",
                                    include.opt.k = TRUE,
                                    return = c("shortlist", "all", "both")) {

    if (is.null(object$scores) || !is.data.frame(object$scores)) {
        stop("Invalid object: missing scores.")
    }

    scores <- object$scores

    directions <- object$directions
    if (is.null(directions)) {
        directions <- c(dunn = "maximize",
                        silhouette = "maximize",
                        ch = "maximize",
                        wb.ratio = "minimize",
                        connectivity = "minimize")
    }

    return <- match.arg(return)

    ## ---------------------------------------------------------------------
    ## Decide which metrics to aggregate, skipping any not computed
    ## ---------------------------------------------------------------------

    if (is.null(metrics)) {
        metrics <- intersect(names(directions), names(scores))
    } else {
        metrics <- intersect(as.character(metrics), names(scores))
        metrics <- intersect(metrics, names(directions))
    }

    metrics <- metrics[vapply(metrics, function(m) any(is.finite(scores[[m]])), logical(1L))]

    if (length(metrics) == 0L) {
        stop("No computed metrics available for aggregation.")
    }

    top.n <- as.integer(top.n)
    if (!is.finite(top.n) || top.n < 1L) top.n <- 10L

    min.metrics <- as.integer(min.metrics)
    if (!is.finite(min.metrics) || min.metrics < 1L) min.metrics <- 1L

    ## ---------------------------------------------------------------------
    ## Compute per-metric ranks (lower is better)
    ## ---------------------------------------------------------------------

    rank.mat <- matrix(NA_real_, nrow = nrow(scores), ncol = length(metrics))
    colnames(rank.mat) <- paste0("rank.", metrics)

    for (j in seq_along(metrics)) {
        m <- metrics[j]
        v <- scores[[m]]
        ok <- is.finite(v)
        if (!any(ok)) next

        r <- rep(NA_real_, length(v))

        if (directions[[m]] == "maximize") {
            r[ok] <- rank(-v[ok], ties.method = ties.method)
        } else if (directions[[m]] == "minimize") {
            r[ok] <- rank(v[ok], ties.method = ties.method)
        } else {
            next
        }

        rank.mat[, j] <- r
    }

    n.metrics.used <- rowSums(is.finite(rank.mat))
    agg.rank <- rowMeans(rank.mat, na.rm = TRUE)

    out <- data.frame(
        k = scores$k,
        agg.rank = agg.rank,
        n.metrics = n.metrics.used,
        stringsAsFactors = FALSE
    )

    out <- cbind(out, as.data.frame(rank.mat, stringsAsFactors = FALSE))
    out <- cbind(out, scores[, metrics, drop = FALSE])

    ## Mark k values that match any per-metric optimum
    out$is.opt <- FALSE
    if (isTRUE(include.opt.k) && !is.null(object$opt.k) && length(object$opt.k) > 0L) {
        opt.k.vec <- as.integer(unlist(object$opt.k, use.names = FALSE))
        opt.k.vec <- opt.k.vec[is.finite(opt.k.vec)]
        if (length(opt.k.vec) > 0L) {
            out$is.opt <- out$k %in% unique(opt.k.vec)
        }
    }

    ## Filter to eligible rows (for ranking output)
    keep <- is.finite(out$agg.rank) & (out$n.metrics >= min.metrics)
    out.eligible <- out[keep, , drop = FALSE]

    ## Full table request
    if (return == "all") {
        ## Return all k rows (still omitting non-computed metric columns by design)
        ## If you want to keep ineligible rows, return 'out' instead of out.eligible.
        return(out.eligible)
    }

    ## Shortlist
    if (nrow(out.eligible) == 0L) {
        if (return == "both") {
            return(list(shortlist = out.eligible, all = out.eligible))
        }
        return(out.eligible)
    }

    out.eligible <- out.eligible[order(out.eligible$agg.rank, -out.eligible$n.metrics, out.eligible$k), , drop = FALSE]
    shortlist <- head(out.eligible, top.n)

    if (isTRUE(include.opt.k) && any(out.eligible$is.opt)) {
        opt.rows <- out.eligible[out.eligible$is.opt, , drop = FALSE]
        shortlist <- rbind(shortlist, opt.rows)
        shortlist <- shortlist[!duplicated(shortlist$k), , drop = FALSE]
        shortlist <- shortlist[order(shortlist$agg.rank, -shortlist$n.metrics, shortlist$k), , drop = FALSE]
    }

    if (return == "both") {
        list(shortlist = shortlist, all = out.eligible)
    } else {
        shortlist
    }
}


#' Louvain clustering wrapper with optional seeding, weights, and repeated runs
#'
#' A convenience wrapper around \code{igraph::cluster_louvain()} that supports:
#' \itemize{
#' \item deterministic RNG control via \code{seed}
#' \item optional edge weights (explicit vector, auto-detect, or forced unweighted)
#' \item repeated runs (\code{n.itrs}) returning a membership matrix
#' }
#'
#' @param graph An \code{igraph} graph object.
#' @param weights Optional edge weights. If \code{NULL} (default), uses the
#'   \code{weight} edge attribute if present; otherwise clustering is unweighted.
#'   If \code{FALSE}, forces unweighted clustering even if \code{weight} exists.
#'   If numeric, must be a vector of length \code{igraph::ecount(graph)}.
#' @param seed Optional integer seed. If not \code{NULL}, the RNG is set for each
#'   iteration as \code{seed + itr - 1}, yielding reproducible but not identical
#'   results across iterations.
#' @param n.itrs Positive integer number of clustering runs. If \code{n.itrs > 1},
#'   returns a membership matrix with one column per run.
#'
#' @return If \code{n.itrs == 1}, an integer membership vector of length
#'   \code{igraph::vcount(graph)}. If \code{n.itrs > 1}, an integer matrix of
#'   dimension \code{vcount(graph)} by \code{n.itrs}; each column is a membership
#'   vector.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_gnp(50, 0.08)
#' E(g)$weight <- runif(ecount(g))
#'
#' ## Uses E(g)$weight by default
#' mem <- cluster.graph.louvain(g)
#'
#' ## Force unweighted
#' mem0 <- cluster.graph.louvain(g, weights = FALSE)
#'
#' ## Repeat runs, reproducible
#' mem.mat <- cluster.graph.louvain(g, seed = 1, n.itrs = 5)
#' }
#'
#' @importFrom igraph cluster_louvain ecount vcount edge_attr_names membership
#' @export
cluster.graph.louvain <- function(graph, weights = NULL, seed = NULL, n.itrs = 1) {
    if (!inherits(graph, "igraph")) stop("'graph' must be an igraph object.")
    n.itrs <- as.integer(n.itrs)
    if (!is.finite(n.itrs) || n.itrs < 1L) stop("'n.itrs' must be a positive integer.")

    ## Resolve weights
    w.use <- NULL
    if (isFALSE(weights)) {
        w.use <- NULL
    } else if (is.null(weights)) {
        if ("weight" %in% igraph::edge_attr_names(graph)) {
            w.use <- igraph::E(graph)$weight
        } else {
            w.use <- NULL
        }
    } else {
        if (!is.numeric(weights)) stop("'weights' must be NULL, FALSE, or a numeric vector.")
        if (length(weights) != igraph::ecount(graph)) stop("'weights' must have length igraph::ecount(graph).")
        w.use <- as.numeric(weights)
    }

    ## CRAN-safe: save/restore .Random.seed if we touch it
    set.seed.local <- function(seed.value) {
        had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        old.seed <- if (had.seed) get(".Random.seed", envir = .GlobalEnv) else NULL
        on.exit({
            if (had.seed) {
                assign(".Random.seed", old.seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed.value)
        invisible(NULL)
    }

    n <- igraph::vcount(graph)

    if (n.itrs == 1L) {
        if (!is.null(seed)) set.seed.local(as.integer(seed))
        cl <- igraph::cluster_louvain(graph, weights = w.use)
        return(as.integer(igraph::membership(cl)))
    }

    mem.mat <- matrix(NA_integer_, nrow = n, ncol = n.itrs)
    for (itr in seq_len(n.itrs)) {
        if (!is.null(seed)) set.seed.local(as.integer(seed) + itr - 1L)
        cl <- igraph::cluster_louvain(graph, weights = w.use)
        mem.mat[, itr] <- as.integer(igraph::membership(cl))
    }
    colnames(mem.mat) <- paste0("itr.", seq_len(n.itrs))

    mem.mat
}

#' Compute Congruence Between Graph Clusters and Labels
#'
#' @description
#' Aligns graph-derived cluster labels to external Labels (for example Community State Type (CST))
#' annotations using \code{sample.ids}, then computes agreement via the Adjusted
#' Rand Index (ARI) and returns a CST-by-cluster contingency table.
#'
#' @details
#' The function subsets \code{cst.labels} to \code{sample.ids} (by name matching),
#' drops samples with missing CST labels, and requires a minimum number of labeled
#' samples (default: 10) to compute ARI. Cluster labels are treated as a
#' partition of \code{sample.ids}; CST labels are coerced to a factor partition.
#'
#' @param sample.ids Character vector of sample identifiers defining the ordering
#'   of \code{cluster.labels}. Typically \code{rownames(X)} used for graph construction.
#' @param cluster.labels Vector of cluster/community assignments aligned to
#'   \code{sample.ids}. Must have the same length as \code{sample.ids}. Can be
#'   numeric/integer or character/factor.
#' @param labels Named vector of external labels, with \code{names(labels)}
#'   giving sample identifiers. Values may be character, factor, or numeric.
#'
#' @return A list with components:
#' \describe{
#'   \item{n}{Number of samples retained after alignment and removal of missing labels.}
#'   \item{ari}{Adjusted Rand Index between \code{cluster.labels} and labels (over retained samples).}
#'   \item{cst.table}{Contingency table of labels by cluster labels (over retained samples).}
#' }
#'
#' @seealso \code{\link[mclust]{adjustedRandIndex}}
#'
#' @examples
#' ## Two clusters vs two CSTs (toy)
#' sample.ids <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10")
#' cluster.labels <- c(1,1,1,1,1, 2,2,2,2,2)
#' labels <- c(s1="I", s2="I", s3="I", s4="I", s5="I",
#'                 s6="III", s7="III", s8="III", s9="III", s10="III")
#' congruence.with.labels(sample.ids, cluster.labels, labels)
#'
#' @export
congruence.with.labels <- function(sample.ids,
                                cluster.labels,
                                labels) {
    ## ---- validate inputs (CRAN-safe defensive checks) ----
    if (missing(sample.ids) || is.null(sample.ids)) {
        stop("`sample.ids` must be a non-null vector of sample identifiers.")
    }
    if (!is.atomic(sample.ids)) stop("`sample.ids` must be an atomic vector.")
    sample.ids <- as.character(sample.ids)
    if (length(sample.ids) < 1L) stop("`sample.ids` must have positive length.")
    if (anyNA(sample.ids)) stop("`sample.ids` must not contain NA.")
    if (any(sample.ids == "")) stop("`sample.ids` must not contain empty strings.")
    if (anyDuplicated(sample.ids)) stop("`sample.ids` must not contain duplicates.")

    if (missing(cluster.labels) || is.null(cluster.labels)) {
        stop("`cluster.labels` must be provided and non-null.")
    }
    if (!is.atomic(cluster.labels)) stop("`cluster.labels` must be an atomic vector.")
    if (length(cluster.labels) != length(sample.ids)) {
        stop("`cluster.labels` must have the same length as `sample.ids`.")
    }
    if (anyNA(cluster.labels)) stop("`cluster.labels` must not contain NA.")

    if (missing(labels) || is.null(labels)) {
        stop("`labels` must be provided and non-null.")
    }
    if (!is.atomic(labels)) stop("`labels` must be an atomic vector.")
    nms <- names(labels)
    if (is.null(nms) || length(nms) != length(labels)) {
        stop("`labels` must be a named vector with names as sample identifiers.")
    }
    if (anyNA(nms) || any(nms == "")) {
        stop("`labels` must have non-missing, non-empty names.")
    }

    ## ---- align CST labels to sample.ids ----
    cst <- labels[sample.ids]
    ok <- !is.na(cst)

    if (sum(ok) < 10L) {
        stop("Too few samples with CST labels after alignment (need at least 10).")
    }

    ## ---- compute ARI ----
    ## coerce to partitions
    cl.ok <- cluster.labels[ok]
    cst.ok <- cst[ok]

    ari <- mclust::adjustedRandIndex(as.integer(as.factor(cl.ok)),
                                     as.integer(as.factor(cst.ok)))

    list(
        n = sum(ok),
        ari = ari,
        cst.table = table(cst.ok, cl.ok)
    )
}
