#' kNN-Based Cluster Imputation
#'
#' This function utilizes the k-Nearest Neighbors (k-NN) algorithm and majority vote to
#' reassign points in the zero cluster to one of the non-zero clusters.
#' The process iterates until no further reduction in the size of the zero cluster is observed.
#'
#' @param X Matrix representing the state space in which the clustering is defined.
#' @param cltr Vector containing the initial clustering labels for each point in 'X'.
#' @param ref.cltr Cluster label designated as the 'zero' cluster for reassignment (default is 0).
#' @param K Number of nearest neighbors used for majority vote-based cluster reassignment (default is 20).
#' @param use.geodesic.dist Logical flag indicating whether to use geodesic distance instead of Euclidean distance (default is FALSE).
#'
#' @return A vector containing the new clustering labels after iterative imputation.
#'
#' @examples
#' # Add illustrative examples here to demonstrate the function's usage
#'
#' @export
## was kNN.soft.cltr
kNN.cltr.imputation <- function(X, cltr, ref.cltr = 0, K = 20, use.geodesic.dist = FALSE) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if ( nrow(X) != length(cltr) ) {
        cat("nrow(X): ", nrow(X), "\n")
        cat("length(cltr): ", length(cltr), "\n")
        stop("nrow(X) != length(cltr)")
    }

    if ( !all(rownames(X) == names(cltr)) ) {
        stop("rownames(X) != names(cltr)")
    }

    K <- min(c(K, nrow(X)-1))

    if ( use.geodesic.dist ) {

        r <- geodesic.knn(X, k = K + 1)
        nn.i <- r$nn.index
        nn.i <- nn.i[,-1] # removing the first column that is the index of the given point

    } else {

        nn <- get.knn(X, k = K)
        nn.i <- nn$nn.index
    }

    cltr2 <- cltr
    idx <- cltr == ref.cltr
    i0 <- which(idx)

    ## size of cluster 0 berfore cltr imputation
    old.cl0.size <- sum(idx)
    cl0.size.diff <- 1

    while ( cl0.size.diff > 0 )
    {
        idx <- cltr2 == ref.cltr
        i0 <- which(idx)
        for ( i in i0 )
        {
            ii <- nn.i[i,]
            nbrd.cl <- cltr2[ii]
            j <- which(nbrd.cl>0)
            if ( length(j)>0 )
            {
                cltr2[i] <- nbrd.cl[j[1]]
            }
        }
        idx <- cltr2==ref.cltr
        new.cl0.size <- sum(idx)
        cl0.size.diff <- old.cl0.size - new.cl0.size
        old.cl0.size <- new.cl0.size
        ##print(new.cl0.size)
    }

    cltr2
}

#' NN-Based Cluster Imputation
#'
#' This function assigns a point not in any cluster to the closest to it cluster.
#'
#' @param X Matrix representing the state space in which the clustering is defined.
#' @param cltr     A vector of clustering labels.
#' @param ref.cltr A cluster label designated as the 'zero' cluster for reassignment. Default is 0.
#'
#' @return A vector containing the new clustering labels after iterative imputation.
#'
#' @examples
#' # Add illustrative examples here to demonstrate the function's usage
#'
#' @export
NN.cltr.imputation <- function(X, cltr, ref.cltr = 0) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if ( nrow(X) != length(cltr) ) {
        cat("nrow(X): ", nrow(X), "\n")
        cat("length(cltr): ", length(cltr), "\n")
        stop("nrow(X) != length(cltr)")
    }

    if ( !all(rownames(X) == names(cltr)) ) {
        stop("rownames(X) != names(cltr)")
    }

    ## Points in the reference cluster
    idx <- cltr == ref.cltr
    i0 <- which(idx)
    Z <- X[i0,]

    ## X - Z
    XmZ <- X[-i0,]

    ## Restricting cltr to XmZ
    XmZ.cltr <- cltr[-i0]

    ## NN's of Z in X
    nn <- get.knnx(XmZ, Z, k = 1)
    nn.i <- nn$nn.index[,1]

    ## Assigning cluster ID of the nearest neighbor NN(x) \in XmZ of x \in Z
    cltr0 <- XmZ.cltr[nn.i]

    cltr[idx] <- cltr0

    cltr
}


#' Reorder cluster IDs
#'
#' This function reorders the cluster IDs of a given vector in decreasing order of cluster sizes.
#' The largest cluster is assigned the ID 1, the second largest gets the ID 2, and so forth.
#'
#' @param cltr A numeric vector where each element represents the cluster ID of a data point.
#'
#' @param decreasing A logical parameter. Set to TRUE if clusters are to be sorted in the decreasing order of their sizes.
#'
#' @return A numeric vector of the same length as the input. The cluster IDs are replaced such that
#' the largest cluster has an ID of 1, the second largest has an ID of 2, and so on.
#'
#' @examples
#' cltr <- c(1,2,2,3,3,3)
#' clusters.reorder(cltr)
#'
clusters.reorder <- function(cltr, decreasing = TRUE)
{
    ## check if labels are numeric
    is.cltr.numeric <- all(is.numeric(cltr))

    ## save current cluster names
    nm <- names(cltr)

    ## Determine the size of each cluster
    cltr.table <- table(cltr)

    ## Order cluster sizes in decreasing order
    cltr.order <- names(cltr.table)[order(cltr.table, decreasing = decreasing)]

    ## Create a named vector that maps the old IDs to the new IDs
    new.ids <- setNames(names(cltr.table), cltr.order)

    ## Relabel the clusters using the 'new.ids' vector
    new.cltr <- new.ids[as.character(cltr)]

    if ( is.cltr.numeric ) {
        new.cltr <- as.numeric(new.cltr)
    }

    names(new.cltr) <- nm

    return(new.cltr)
}


#' Apply HDBSCAN Clustering
#'
#' This function applies HDBSCAN clustering to a given dataset (expressed as a
#' matrix), and evaluates clustering results with different minimum cluster
#' sizes using dunn and connectivity indices.
#'
#' @param X             A matrix of the data to be clustered.
#' @param method        A method (index) to be used for identifying the optimal number of clusters. Default is 'dunn'. Other options: 'connectivity', 'cl0.size' and 'all'.
#' @param min.pts       A vector of integers indicating the minimum cluster sizes to be evaluated. Default is 5:50.
#' @param min.prop      The minimum proportion of nrow(X) to be used in test clusterings.
#' @param max.prop      The maxinum proportion of nrow(X) to be used in test clusterings.
#' @param n.test.cltr   The number of test clusters.
#' @param soft.K        The number of nearest neighbors to consider when softening the cluster boundaries. Default is 20.
#' @param n.cores       The number of cores to use.
#' @param verbose       A logical value indicating whether progress should be displayed. Default is FALSE.
#'
#' @return A list containing the details of the clustering for three different metrics:
#' cl0.size (size of the cluster labeled as '0'), dunn.idx (Dunn index, a measure of cluster quality),
#' and connectivity.idx (connectivity index). For each metric, the optimal cluster size,
#' the extended cluster (after softening), and the frequencies of the clusters are returned.
#'
#' @examples
#' # Assuming 'X' is your data matrix
#' # X <- matrix(...)
#' # Then you can call the function as:
#' # hdbscan.cltr(X, min.pts = 5:50, soft.K = 20, verbose = FALSE)
#'
#' @importFrom dbscan hdbscan
#' @importFrom clValid dunn connectivity
#'
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
    if ( !requireNamespace("clValid", quietly = TRUE) ) {
        stop("The clValid package is not installed. Please install it before using this routine.")
    }

    if ( !requireNamespace("dbscan", quietly = TRUE) ) {
        stop("The dbscan package is not installed. Please install it before using this routine.")
    }

    if ( is.null(min.pts) ) {

        stopifnot( 0 < min.prop && min.prop < 1 )
        stopifnot( 0 < max.prop && max.prop < 1 )

        min.pts <- seq(min.prop * nrow(X), max.prop * nrow(X), length.out = n.test.cltr)
        min.pts <- as.integer(min.pts)
    }

    ## cat("length(min.pts): ", length(min.pts), "\n")

    if ( length(min.pts) < n.cores ) {
        n.cores <- length(min.pts)
        ## cat("Set n.cores to ", n.cores, "\n")
    }

    if ( any(min.pts >= nrow(X)) ) {
        idx <- min.pts < nrow(X)
        min.pts <- min.pts[idx]
        cat("WARNING: Restricting min.pts to the range with min.pts < nrow(X)\n")
    }

    if ( length(min.pts) == 0 ) {
        stop("min.pts is empty!")
    }

    dX.mat <- as.matrix(dist(X))

    cl0.size <- NA
    n.cltrs <- NA
    cl0.cltr <- NA
    cl0.cltr.ext <- NA
    cl0.k <- NA
    cl0.freq <- c()
    cl0.freq['0'] <- NA
    cl0.freq.ext <- NA
    ##
    dunn.idx <- NA
    dunn.cltr <- NA
    dunn.cltr.ext <- NA
    dunn.k <- NA
    dunn.freq <- c()
    dunn.freq['0'] <- NA
    dunn.freq.ext <- NA
    ##
    connectivity.idx <- NA
    connectivity.cltr <- NA
    connectivity.cltr.ext <- NA
    connectivity.k <- NA
    connectivity.freq <- c()
    connectivity.freq['0'] <- NA
    connectivity.freq.ext <- NA

    cltrgs <- list()
    cl0.size <- c()
    n.cltrs <- c()
    dunn.idx <- c()
    connectivity.idx <- c()

    if ( method == "all" )
    {
        if ( n.cores == 1 )
        {
            for ( k in min.pts )
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)
                cltrgs[[k]] <- cltr
                if ( length(unique(cltr)) > 1 )
                {
                    dunn.idx[k] <- clValid::dunn(dX.mat, clusters=cltr)
                    idx <- cltr != 0
                    if ( sum(idx) != 0 )
                    {
                        connectivity.idx[k] <- clValid::connectivity(clusters=cltr[idx], Data=X[idx,])
                    } else {
                        connectivity.idx[k] <- NA
                    }
                    f <- table(cltr)
                    s <- f['0'][[1]]
                    if ( !is.na(s) ) {
                        cl0.size[k] <- s
                    } else {
                        cl0.size[k] <- 0
                    }
                    n.cltrs[k] <- max(cltr)
                }
            }
        } else {

            idx.fn <- function(k, X, verbose = TRUE)
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)

                n.cltrs <- max(cltr)

                dunn.idx <- 0
                connectivity.idx <- 0
                cl0.size <- 0

                if ( length(unique(cltr)) > 1 )
                {
                    dunn.idx <- clValid::dunn(dX.mat, clusters=cltr)
                    idx <- cltr != 0
                    if ( sum(idx) != 0 )
                    {
                        connectivity.idx <- clValid::connectivity(clusters=cltr[idx], Data=X[idx,])
                    }
                    f <- table(cltr)
                    s <- f['0'][[1]]
                    if ( !is.na(s) ) {
                        cl0.size <- s
                    }
                }

                list(cltr=cltr,
                     n.cltrs=n.cltrs,
                     dunn.idx=dunn.idx,
                     connectivity.idx=connectivity.idx,
                     cl0.size=cl0.size)
            }


            registerDoParallel(n.cores)
            ##ptm <- proc.time()
            res <- foreach ( k = min.pts ) %dopar% {
                idx.fn(k, X, verbose)
            }
            ##elapsed.time(ptm)
            stopImplicitCluster()

            for ( j in seq(min.pts) )
            {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[k]] <- r$cltr
                cl0.size[k] <- r$cl0.size
                n.cltrs[k] <- r$n.cltrs
                dunn.idx[k] <- r$dunn.idx
                connectivity.idx[k] <- r$connectivity.idx
            }
        }

        cl0.k <- which.min(cl0.size)
        cl0.cltr <- cltrgs[[cl0.k]]
        cl0.freq <- table(cl0.cltr)
        cl0.cltr.ext <- NULL
        cl0.freq.ext <- NULL
        if ( '0' %in% names(cl0.freq) )
        {
            cl0.cltr.ext <- kNN.cltr.imputation(X, cl0.cltr, K=soft.K)
            cl0.freq.ext <- table(cl0.cltr.ext)
        }

        dunn.k <- which.max(dunn.idx)
        dunn.cltr <- cltrgs[[dunn.k]]
        dunn.freq <- table(dunn.cltr)
        dunn.cltr.ext <- NULL
        dunn.freq.ext <- NULL
        if ( '0' %in% names(dunn.freq) )
        {
            dunn.cltr.ext <- kNN.cltr.imputation(X, dunn.cltr, K=soft.K)
            dunn.freq.ext <- table(dunn.cltr.ext)
        }

        connectivity.k <- which.min(connectivity.idx)

        connectivity.cltr <- NULL
        connectivity.freq <- NULL
        connectivity.cltr.ext <- NULL
        connectivity.freq.ext <- NULL

        if ( length(connectivity.k) > 0 )
        {
            connectivity.cltr <- cltrgs[[connectivity.k]]
            connectivity.freq <- table(connectivity.cltr)
            if ( '0' %in% names(connectivity.freq) )
            {
                connectivity.cltr.ext <- kNN.cltr.imputation(X, connectivity.cltr, K=soft.K)
                connectivity.freq.ext <- table(connectivity.cltr.ext)
            }
        }

    } else if ( method == "dunn" ) {

        if ( n.cores == 1 )
        {
            for ( k in min.pts )
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)
                cltrgs[[k]] <- cltr
                if ( length(unique(cltr)) > 1 ) {
                    dunn.idx[k] <- clValid::dunn(dX.mat, clusters=cltr)
                    n.cltrs[k] <- max(cltr)
                }
            }
        } else {

            idx.fn <- function(k, X, verbose = TRUE)
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)

                n.cltrs <- max(cltr)

                dunn.idx <- 0
                if ( length(unique(cltr)) > 1 ) {
                    dunn.idx <- clValid::dunn(dX.mat, clusters=cltr)
                }

                list(cltr=cltr,
                     n.cltrs=n.cltrs,
                     dunn.idx=dunn.idx,
                     connectivity.idx=NA,
                     cl0.size=NA)
            }


            registerDoParallel(n.cores)
            ##ptm <- proc.time()
            res <- foreach ( k = min.pts ) %dopar% {
                idx.fn(k, X, verbose)
            }
            ##elapsed.time(ptm)
            stopImplicitCluster()

            for ( j in seq(min.pts) )
            {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[k]] <- r$cltr
                n.cltrs[k] <- r$n.cltrs
                dunn.idx[k] <- r$dunn.idx
            }
        }

        dunn.k <- which.max(dunn.idx)
        dunn.cltr <- cltrgs[[dunn.k]]
        dunn.freq <- table(dunn.cltr)
        dunn.cltr.ext <- NULL
        dunn.freq.ext <- NULL
        if ( '0' %in% names(dunn.freq) )
        {
            dunn.cltr.ext <- kNN.cltr.imputation(X, dunn.cltr, K=soft.K)
            dunn.freq.ext <- table(dunn.cltr.ext)
        }

    } else if ( method == "connectivity" ) {

        if ( n.cores == 1 )
        {
            for ( k in min.pts )
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)
                cltrgs[[k]] <- cltr
                n.cltrs[k] <- max(cltr)
                if ( length(unique(cltr)) > 1 )
                {
                    idx <- cltr != 0
                    if ( sum(idx) != 0 )
                    {
                        connectivity.idx[k] <- clValid::connectivity(clusters=cltr[idx], Data=X[idx,])
                    } else {
                        connectivity.idx[k] <- NA
                    }
                }
            }
        } else {

            idx.fn <- function(k, X, verbose = TRUE)
            {
                if ( verbose )
                    cat("\r",k)
                cltr <- dbscan::hdbscan(X, minPts=k)$cluster
                names(cltr) <- rownames(X)

                n.cltrs <- max(cltr)

                dunn.idx <- 0
                connectivity.idx <- 0
                cl0.size <- 0

                if ( length(unique(cltr)) > 1 )
                {
                    idx <- cltr != 0
                    if ( sum(idx) != 0 )
                    {
                        connectivity.idx <- clValid::connectivity(clusters=cltr[idx], Data=X[idx,])
                    }
                }

                list(cltr=cltr,
                     n.cltrs=n.cltrs,
                     dunn.idx=NA,
                     connectivity.idx=connectivity.idx,
                     cl0.size=NA)
            }

            registerDoParallel(n.cores)
            ptm <- proc.time()
            res <- foreach ( k = min.pts ) %dopar% {
                idx.fn(k, X, verbose)
            }
            elapsed.time(ptm)
            stopImplicitCluster()

            for ( j in seq(min.pts) )
            {
                k <- min.pts[j]
                r <- res[[j]]
                cltrgs[[k]] <- r$cltr
                n.cltrs[k] <- r$n.cltrs
                connectivity.idx[k] <- r$connectivity.idx
            }
        }

        connectivity.k <- which.min(connectivity.idx)
        connectivity.cltr <- NULL
        connectivity.freq <- NULL
        connectivity.cltr.ext <- NULL
        connectivity.freq.ext <- NULL

        if ( length(connectivity.k) > 0 )
        {
            connectivity.cltr <- cltrgs[[connectivity.k]]
            connectivity.freq <- table(connectivity.cltr)
            if ( '0' %in% names(connectivity.freq) )
            {
                connectivity.cltr.ext <- kNN.cltr.imputation(X, connectivity.cltr, K=soft.K)
                connectivity.freq.ext <- table(connectivity.cltr.ext)
            }
        }
    }

    list(cl0.size=cl0.size,
         n.cltrs=n.cltrs,
         cl0.cltr=cl0.cltr,
         cl0.cltr.ext=cl0.cltr.ext,
         cl0.opt.k=cl0.k,
         cl0.cl0.size=cl0.freq['0'][[1]],
         cl0.n.cltrs=length(cl0.freq),
         cl0.cltr.freq=cl0.freq,
         cl0.cltr.ext.freq=cl0.freq.ext,
         ##
         dunn.idx=dunn.idx,
         dunn.cltr=dunn.cltr,
         dunn.cltr.ext=dunn.cltr.ext,
         dunn.opt.k=dunn.k,
         dunn.cl0.size=dunn.freq['0'][[1]],
         dunn.n.cltrs=max(dunn.cltr),
         dunn.cltr.freq=dunn.freq,
         dunn.cltr.ext.freq=dunn.freq.ext,
         ##
         connectivity.idx=connectivity.idx,
         connectivity.cltr=connectivity.cltr,
         connectivity.cltr.ext=connectivity.cltr.ext,
         connectivity.opt.k=connectivity.k,
         connectivity.cl0.size=connectivity.freq['0'][[1]],
         connectivity.n.cltrs=max(connectivity.cltr),
         connectivity.cltr.freq=connectivity.freq,
         connectivity.cltr.ext.freq=connectivity.freq.ext)
}



#' Performs Inductive Clustering on a List of Clusterings
#'
#' This function takes as input a list of clusterings, where each item is a
#' sub-clustering of the largest cluster in the previous item, and outputs a
#' new clustering of the entire data set constructed by successively subdividing
#' the largest cluster at each step.
#'
#' @param cltr.list A list of integer vectors representing the clusterings.
#' The i-th element of the k-th vector indicates the cluster of the i-th data
#' point in the k-th clustering.
#'
#' @return An integer vector representing the final clustering.
#'
#' @examples
#' # Suppose cltr.list is a list of clusterings
#' # cltr.list <- list(c(1,1,2,2), c(1,2), c(1,2,2), ...)
#' # final_clustering <- inductive.clustering(cltr.list)
#'
inductive.clustering <- function(cltr.list)
{

    if (!is.list(cltr.list)) {
        stop("Input must be a list")
    }

    for ( i in seq(cltr.list) ) {

        item <- cltr.list[[i]]

        if (!is.vector(item) || !all(as.character(as.integer(item)) == as.character(item)) ) {

            if ( !is.vector(item) ) {
                print(paste("Cluster ",i," is not a vector!"))
            }

            if ( !all(as.character(as.integer(item)) == as.character(item)) ) {
                print(paste("Cluster ",i," is not integer valued"))
            }

            file <- tempfile(fileext = ".rda")
            save(cltr.list, file=file)
            print(paste("cltr.list saved to ",file))


            stop("All items in the list must be integer vectors")
        }
    }

    main.cltr <- cltr.list[[1]]

    for ( k in 2:length(cltr.list) ) {

        sub.cltr <- cltr.list[[k]]
        sub.cltr.inds <- which(main.cltr == 1)

        if ( length(sub.cltr.inds) != length(sub.cltr) ) {

            cat("k=",k,"\n")
            cat("length(sub.cltr.inds)=", length(sub.cltr.inds), "\n")
            cat("length(sub.cltr)=", length(sub.cltr), "\n")

            file <- tempfile(fileext = ".rda")
            save(cltr.list, file=file)
            print(paste("cltr.list saved to ",file))

            stop("Sub-clusterings size cannot be different from the size of cluster 1 of main clustering.")
        }

        ## replace cluster 1 with sub.cltr
        main.cltr[sub.cltr.inds] <- sub.cltr

        ## change the indices of the remaining clusters
        other.inds <- setdiff(seq(main.cltr), sub.cltr.inds)
        zero.inds <- which(main.cltr == 0)
        if ( length(zero.inds) > 0 ) {
            other.inds <- setdiff(other.inds, zero.inds)
        }
        main.cltr[other.inds] <- main.cltr[other.inds] + max(sub.cltr) - 1
    }

    return(main.cltr)
}

## #' Splits recursively the largest cluster into smaller ones.
## #'
## #' The algorithm splits the largest cluster of Y (if its size is greater than
## #' cltr.size.thld) by selecting the samples of the cluster in X and then
## #' generating the low dimensional embedding of the resulting space and then
## #' clustering it.
## #'
## #' @param X              A source data table from which Y was constructed.
## #' @param Y              A low dimensional representation of X.
## #' @param cltr           A clustering of Y.
## #' @param cltr.size.thld The cluster size threshold. The clustering stops when the largest cluster size is smaller than cltr.size.thld or the resulting space has only one cluster. Default is 100.
## #' @param min.pts        A vector of minimal cluster size values for hdbscan.cltr(). Default is 5*(1:10).
## #' @param min.prop       The proportion of the number of samples in the current maximal cluster that will be taken as the min.pts in HDBSCAN's k parameter. More precisely, k = min( c(cltr.size.thld, min.prop * nrow(max.cltr.Y)) ). Default is 0.1.
## #' @param max.prop       The maxinum proportion of nrow(X) to be used in test clusterings.
## #' @param n.test.cltr    The number of test clusters.
## #' @param distance       The distance in the PaCMAP algorithm. Default is 'angular'.
## #' @param dim            The dimenion into which \code{X\[max cltr samples,\]} is embedded. Default is 3.
## #' @param n.cores        The number of cores to use in hdbscan.cltr(). Default is 10.
## #' @param plot.it        Logical. If TRUE a plot showing the clustering of the maximal sub-cluster is shown. Default is TRUE.
## #' @param verbose        A logical value indicating whether detailed output should be printed during the function's execution. If TRUE, the function will print out additional details about the calculations it's performing. Default is TRUE.
## maxcltr.split <- function(X,
##                           Y,
##                           cltr,
##                           cltr.size.thld = 100,
##                           min.pts = 5*(1:10),
##                           min.prop = 0.1,
##                           max.prop = 0.5,
##                           n.test.cltr = 10,
##                           distance = "angular",
##                           dim = 3,
##                           n.cores = 10,
##                           plot.it = TRUE,
##                           verbose = TRUE)
## {
##     stopifnot(rownames(X) == names(cltr))
##     stopifnot(rownames(Y) == names(cltr))
##     stopifnot(as.character(as.integer(cltr)) == as.character(cltr))

##     if ( verbose ) {
##         fn.ptm <- proc.time()
##     }

##     ## Since clusters are integers, let's explicitly make them integers
##     cn <- names(cltr)
##     cltr <- as.integer(cltr)
##     names(cltr) <- cn

##     ## Computing sizes of different clusters and sorting them by size
##     cltr.freq <- sort(table(cltr[cltr != 0]), decreasing = TRUE)
##     n.cltrs <- length(cltr.freq)

##     ## Selecting the maximal cluster sample IDs
##     max.cltr.i <- as.integer(names(cltr.freq)[which.max(cltr.freq)])
##     idx <- cltr == max.cltr.i
##     max.cltr.ids <- names(cltr)[idx]
##     max.cltr.size <- length(max.cltr.ids)

##     cltr.list <- list()
##     cltr.list[[1]] <- cltr

##     Y.cltr.list <- list()
##     Y.cltr.list[[1]] <- cltr

##     Y.list <- list()
##     Y.list[[1]] <- Y

##     k <- 2
##     init.reticulate <- TRUE
##     while ( n.cltrs > 1 && max.cltr.size > cltr.size.thld ) {

##         if ( verbose ) {
##             print(paste("Splitting level ", k-1))
##         }

##         ## Generating PaCMAP embedding of X[cl.ids,] into R^dim
##         max.cltr.Y <- Rpacmap(X[max.cltr.ids,], ndim = 3, init.reticulate = init.reticulate, distance = distance)
##         init.reticulate <- FALSE

##         if ( verbose ) {
##             cat("Generating clustering of the gene embedding ... \n")
##             ptm <- proc.time()
##         }

##         if ( !is.null(min.prop) && !is.null(max.prop) ) {

##             r <- hdbscan.cltr(max.cltr.Y,
##                              min.prop = min.prop,
##                              max.prop = max.prop,
##                              n.test.cltr = n.test.cltr)
##             cltr <- r$dunn.cltr

##         } else if ( !is.null(min.prop) ) {

##             stopifnot( 0 < min.prop && min.prop < 1 )

##             L <- as.integer(min( c(cltr.size.thld, min.prop * nrow(max.cltr.Y)) ))
##             cltr <- dbscan::hdbscan(max.cltr.Y, minPts = L)$cluster
##             names(cltr) <- rownames(max.cltr.Y)

##         } else {

##             max.cltr.min.pts <- min.pts
##             idx <- max.cltr.min.pts < max.cltr.size
##             max.cltr.min.pts <- max.cltr.min.pts[idx]
##             max.soft.K <- min(c(20, length(max.cltr.ids))) - 1

##             r <- hdbscan.cltr(max.cltr.Y, min.pts = max.cltr.min.pts, soft.K = max.soft.K, n.cores = n.cores, verbose = FALSE)
##             cltr <- r$dunn.cltr
##         }

##         cltr[cltr != 0] <- clusters.reorder(cltr[cltr != 0])

##         if ( verbose ) {
##             elapsed.time(ptm)
##             print(table(cltr))
##         }

##         ## Computing sizes of different clusters and sorting them by size
##         cltr.freq <- sort(table(cltr[cltr != 0]), decreasing = TRUE)
##         n.cltrs <- length(cltr.freq)

##         if ( 0 %in% cltr && n.cltrs > 1 ) {

##             if ( verbose ) {
##                 cat("Reducing the size of 0 cluster with soft clustering ... \n")
##                 ptm <- proc.time()
##             }

##             cltr <- kNN.cltr.imputation(max.cltr.Y, cltr)

##             if ( verbose ) {
##                 elapsed.time(ptm)
##                 print(table(cltr))
##             }
##         }

##         cltr <- as.integer(cltr)
##         names(cltr) <- max.cltr.ids

##         if ( plot.it ) {
##             plot3d.cltrs(max.cltr.Y, cltr, show.cltr.labels = TRUE, sort.legend.labs.by.freq = TRUE)
##         }

##         Y.list[[k]] <- max.cltr.Y
##         Y.cltr.list[[k]] <- cltr
##         cltr.list[[k]] <- inductive.clustering(Y.cltr.list)
##         k <- k + 1

##         ## Selecting the maximal cluster sample IDs
##         max.cltr.i <- as.integer(names(cltr.freq)[which.max(cltr.freq)])
##         idx <- cltr == max.cltr.i
##         max.cltr.ids <- names(cltr)[idx]
##         max.cltr.size <- length(max.cltr.ids)
##     }

##     new.cltr <- inductive.clustering(Y.cltr.list)

##     if ( 0 %in% new.cltr ) {

##         if ( verbose ) {
##             cat("Reducing the size of 0 cluster in new.cltr with soft clustering ... \n")
##             ptm <- proc.time()
##         }

##         new.cltr <- kNN.cltr.imputation(Y, new.cltr)

##         if ( verbose ) {
##             elapsed.time(ptm)
##             print(table(cltr))
##         }
##     }

##     if ( verbose ) {
##         txt <- "maxcltr.split() total elapsed time:"
##         elapsed.time(fn.ptm, txt, with.brackets=FALSE)
##     }

##     list(Y.list=Y.list,
##          Y.cltr.list=Y.cltr.list,
##          cltr.list=cltr.list,
##          new.cltr=new.cltr)
## }


#' Splits recursively the cluster of the given clustering into smaller clusters.
#'
#' The algorithm iteratively applies maxcltr.split() function to the given
#' clustering so that all clusters of size greater than cltr.size.thld are split
#' (if possible) into clusters of smaller size.
#'
#' @param X              A source data table from which Y was constructed.
#' @param Y              A low dimensional representation of X.
#' @param cltr           A clustering of Y.
#' @param cltr.size.thld The cluster size threshold. The clustering stops when the largest cluster size is smaller than cltr.size.thld or the resulting space has only one cluster. Default is 100.
#' @param min.pts        A vector of minimal cluster size values for hdbscan.cltr(). Default is 5*(1:10).
#' @param min.prop       The proportion of the number of samples in the current maximal cluster that will be taken as the min.pts in HDBSCAN's k parameter. More precisely, k = min( c(cltr.size.thld, min.prop * nrow(max.cltr.Y)) ). Default is 0.1.
#' @param max.prop       The maxinum proportion of nrow(X) to be used in test clusterings.
#' @param n.test.cltr    The number of test clusters.
#' @param distance       The distance in the PaCMAP algorithm. Default is 'angular'.
#' @param dim            The dimenion into which \code{X\[max cltr samples,\]} is embedded. Default is 3.
#' @param n.cores        The number of cores to use in hdbscan.cltr(). Default is 10.
#' @param plot.it        Logical. If TRUE a plot showing the clustering of the maximal sub-cluster is shown. Default is TRUE.
#' @param verbose        A logical value indicating whether detailed output should be printed during the function's execution. If TRUE, the function will print out additional details about the calculations it's performing. Default is TRUE.
cltrs.split <- function(X,
                       Y,
                       cltr,
                       cltr.size.thld = 100,
                       min.pts = 5*(1:10),
                       min.prop = 0.1,
                       max.prop = 0.5,
                       n.test.cltr = 10,
                       distance = "angular",
                       dim = 3,
                       n.cores = 10,
                       plot.it = TRUE,
                       verbose = TRUE)
{
    if ( verbose ) {
        fn.ptm <- proc.time()
    }

    cltr[cltr != 0] <- clusters.reorder(cltr[cltr != 0])
    cltr.freq <- table(cltr[cltr != 0])
    max.cltr.size <- cltr.freq["1"][[1]]

    while ( max.cltr.size > cltr.size.thld ) {

        r <- maxcltr.split(X,
                           Y,
                           cltr,
                           cltr.size.thld,
                           min.pts,
                           min.prop,
                           max.prop,
                           n.test.cltr,
                           distance,
                           dim,
                           n.cores,
                           plot.it,
                           verbose)

        cltr <- r$new.cltr
        cltr[cltr != 0] <- clusters.reorder(cltr[cltr != 0])
        cltr.freq <- table(cltr[cltr != 0])
        max.cltr.size <- cltr.freq["1"][[1]]
    }

    if ( verbose ) {
        txt <- "cltrs.split() total elapsed time"
        elapsed.time(fn.ptm, txt, with.brackets=FALSE)
    }

    return(cltr)
}


#' Generates Synthetic Clustering
#'
#' This function generates a synthetic clustering by randomly assigning each point to a cluster.
#' The clusters are then labeled with consecutive integers, with the largest cluster having an index of 1.
#'
#' @param size Integer. The total number of points in the clustering.
#' @param num_clusters Integer. The total number of clusters.
#'
#' @return An integer vector of length 'size' where each element represents the cluster assignment of the corresponding point.
#' The cluster labels are integers from 1 to 'num_clusters' with the largest cluster labeled as 1.
#'
#' @examples
#' generate.cluster(100, 5)
#'
#' @export
generate.cluster <- function(size, num_clusters) {
    ## Randomly assign each point to a cluster
    cltr <- sample(1:num_clusters, size, replace = TRUE)
    ## make sure clusters are labeled with consecutive integers
    cltr.freq <- table(cltr)
    if ( length(cltr.freq) != max(cltr) ) {
        cltr <- relabel.int.cltrg(cltr)
    }

    clusters.reorder(cltr) # now the largest cluster has index 1 and cluster indices are ordered in decreasing order of their size
}

#' Relabel Integer Clustering
#'
#' This function takes an integer-valued clustering vector and relabels the clusters
#' with consecutive integers starting from 1. The function first checks that the
#' input clustering is integer-valued and then performs the relabeling.
#'
#' @param cltr A named vector representing the clustering, where the names of the elements
#'        correspond to the objects being clustered, and the values (integers) represent
#'        the cluster assignments.
#'
#' @return A named integer vector with the same names as the input, but with clusters
#'         relabeled with consecutive integers.
#'
#' @examples
#' # Create a clustering with non-consecutive integer labels
#' cltr <- c(A = 3, B = 3, C = 5)
#' # Relabel the clustering with consecutive integers
#' relabeled_cltr <- relabel.int.cltrg(cltr)
#'
#' @note The function assumes that the input clustering is integer-valued, and it will
#'       stop with an error message if this condition is not met.
#'
#' @export
relabel.int.cltrg <- function(cltr) {

    ## checking if the clustering is integer-valued
    stopifnot(as.character(as.integer(cltr)) == as.character(cltr))

    cn <- names(cltr)
    cltr <- as.integer(cltr)
    names(cltr) <- cn

    cltr.freq <- table(cltr)

    relabel.tbl <- seq(cltr.freq)
    names(relabel.tbl) <- names(cltr.freq)
    cltr <- as.integer(relabel.tbl[as.character(cltr)])

    cltr
}

#' Plot Taxons
#'
#' This function plots the taxons given the index 'i' from the taxon's file. It also has an option to show or hide the cluster labels.
#'
#' @param i Integer. The index of the taxon in the taxon's file.
#' @param taxons.g10 Vector. A vector of taxon names.
#' @param show.cltr.labels Logical. Whether to show the cluster labels in the plot. Defaults to TRUE.
#'
#' @return Returns an invisible list containing the taxon, the file path of the clustering, the coordinates of the points, and the cluster assignments.
#' This function is mainly used for its side effect of creating a plot.
#'
#' @examples
#' tx.plot(1, taxons.g10)
#'
#' @export
tx.plot <- function(i, taxons.g10, show.cltr.labels = TRUE)
{
    if ( is.null(taxons.g10) ) stop("taxons.g10 cannot be NULL")

    taxon <- taxons.g10[i]

    data.dir = "~/projects/rllm-projects/projects/ZB/data/virgo2_sp_tbls/"
    pacmap.file <- paste0(data.dir,taxon,"_pacmap.csv")
    X <- read.csv(pacmap.file, row.names=1)

    cltr.ext.file <- paste0(data.dir,taxon,"_pacmap_cltr_ext.csv")
    cltr <- read.csv(cltr.ext.file, row.names=1)[,1]

    plot3d.cltrs(X, cltr, show.cltr.labels = show.cltr.labels, sort.legend.labs.by.freq = TRUE, legend.title = taxon)

    names(cltr) <- rownames(X)

    invisible(list(taxon=taxon,
                   cltr.ext.file=cltr.ext.file,
                   X=X,
                   cltr=cltr))
}


#' Function to perform a BFS and split clusters based on log likelihood criteria
#'
#' Function to perform a Breadth-First Search (BFS) and split clusters based on log likelihood criteria.
#' The function iteratively traverses a dendrogram, grouping leaves into clusters where the maximum difference
#' in log likelihood values within the cluster is no more than a specified threshold (epsilon).
#'
#' @param tree A hierarchical clustering dendrogram (typically created using 'hclust') to be divided into clusters.
#' @param logL A numerical vector of log likelihood values corresponding to the leaves in the dendrogram. The order
#'             of the values in logL should match the original order of observations in the input to 'hclust'.
#' @param epsilon A numerical threshold for the maximum allowed difference in log likelihood values within a cluster.
#'                If the difference between the maximum and minimum log likelihood values within a potential cluster
#'                exceeds this threshold, the subtree is further divided.
#'
#' @return A list of clusters, where each cluster is represented by a vector of indices corresponding to the leaves
#'         in the original dendrogram. The clusters satisfy the condition that the maximum difference in log likelihood
#'         values within each cluster is no more than the specified epsilon.
#'
#' @examples
#' data <- matrix(rnorm(100), ncol = 2)
#' logL <- rnorm(50)
#' hc <- hclust(dist(data))
#' epsilon <- 1.0
#' clusters <- loglikelihood.cut.tree(hc, logL, epsilon)
#'
loglikelihood.cut.tree <- function(tree, logL, epsilon) {

    tree <- as.dendrogram(tree)

    clusters <- list()
    queue <- list(tree)

    ## Function to check if the log-likelihood difference within a cluster exceeds epsilon
    check.cluster.log.likelihood.range <- function(cluster) {
        max.diff <- max(logL[cluster]) - min(logL[cluster])
        max.diff <= epsilon
    }

    while (length(queue) > 0) {

        current_tree <- queue[[1]]
        queue <- queue[-1]

        ## Get the leaves of the current subtree
        leaves <- order.dendrogram(current_tree)

        ## Diagnostic print statements
        ## print("Leaves:")
        ## print(leaves)
        ## print("logL values:")
        ## print(logL[leaves])

        ## Check log likelihood condition
        if ( check.cluster.log.likelihood.range(leaves) ) {
            ## Add the leaves as a finalized cluster
            clusters <- append(clusters, list(leaves))
        } else {
            ## Get the children of the current subtree
            children <- list(current_tree[[1]], current_tree[[2]])
            ## Add the children to the queue if they are dendrogram objects
            queue <- append(queue, lapply(children, function(x) if(class(x) == "dendrogram") x))
        }
    }

    return(clusters)
}


#' Function to create a vector of cluster indices based on cluster sizes
#'
#' @param clusters A list of clusters, where each cluster is represented by a vector of indices
#' @return A vector of cluster indices, indexed in the decreasing order of their size
#'
cluster.indices.ordered.by.size <- function(clusters) {
  ## Calculate the size of each cluster
  cluster.sizes <- sapply(clusters, length)

  ## Order the clusters by size in decreasing order
  ordered.cluster.indices <- order(cluster.sizes, decreasing = TRUE)

  ## Create a result vector to store the final cluster indices
  result <- numeric(length(clusters))

  ## Assign the ordered indices to the corresponding clusters
  for (i in seq(ordered.cluster.indices)) {
    cluster.index <- ordered.cluster.indices[i]
    result[clusters[[cluster.index]]] <- i
  }

  return(result)
}

#' Create Species-Level Cluster of co-Abundant Genes (CAGs)
#'
#' This function generates species-level CAGs using a defined taxon. It identifies
#' non-outlier and outlier data points using a specified method and forms CAGs accordingly.
#'
#' @param i An integer, indicating the index or identifier of the taxon.
#' @param taxon An optional parameter to specify a taxon. Default is NULL.
#' @param min.pts       The minimal number of points in the hdbscan() algorithm.
#' @param cltr.from.scratch A logical value indicating whether to cluster r$X from scratch.
#' @param n.cores  The number of cores to use.
#' @param verbose       A logical value indicating whether detailed output should be printed during the function's execution. If TRUE, the function will print out additional details about the calculations it's performing. Default is TRUE.
#'
#' @return A list with two components: the output of show.tx(i) and the given taxon newly formed CAGs.
#'
#' @examples
#' # Assuming 'i' is the identifier of your taxon
#' # i <- ...
#' # Then you can call the function as:
#' # create.sp.CAGs(i, taxon = NULL)
#'
create.sp.CAGs <- function(i,
                          taxon=NULL,
                          min.pts=5*(1:10),
                          cltr.from.scratch = TRUE,
                          n.cores = 10,
                          verbose = TRUE)
{
    if ( verbose ) {
        ptm <- proc.time()
    }

    r <- show.tx(i, taxon, cltr.from.scratch = cltr.from.scratch,
                reorder.cltrs = TRUE,
                show.plot = FALSE,
                n.cores = n.cores,
                min.pts = min.pts,
                verbose = verbose)

    cltr <- r$cltr.ext
    names(cltr) <- rownames(r$X)

    ## Identifying outlier genes
    K <- min(100, nrow(r$X)) - 1
    rr <- rm.SS.outliers(r$X, p=0.98, K=K, dist.factor=75)#, method = "dist.factor")
    outlier.ids <- setdiff(rownames(r$X), rownames(rr$S.q))

    if ( length(outlier.ids) > 1 ) {

        if ( 0 %in% cltr ) {
            zero.ids <- names(cltr)[cltr == 0]
            outlier.ids <- c(outlier.ids, zero.ids)
            outlier.ids <- unique(outlier.ids)
        }

        ## Clustering of outlier.ids
        max.min.pts <- min(c(50, length(outlier.ids))) - 1
        max.soft.K <- min(c(20, length(outlier.ids))) - 1
        r3 <- hdbscan.cltr(r$X[outlier.ids,], min.pts=2:max.min.pts, soft.K=max.soft.K, n.cores=n.cores, verbose=FALSE, method="dunn")
        outlier.cltr <- clusters.reorder(r3$dunn.cltr)

        sp.cags.ext <- paste0(r$taxon, "__CAG_E", outlier.cltr)
        names(sp.cags.ext) <- outlier.ids
        ##head(sp.cags.ext)

        non.outlier.ids <- setdiff(rownames(r$X), outlier.ids)
        sp.cags <- paste0(r$taxon, "__CAG_", cltr[non.outlier.ids])
        names(sp.cags) <- non.outlier.ids

        sp.cags <- c(sp.cags, sp.cags.ext)
        sp.cags <- sp.cags[rownames(r$X)]

    } else if ( length(outlier.ids) == 1 ) {

        sp.cags.ext <- paste0(r$taxon, "__CAG_E1")
        names(sp.cags.ext) <- outlier.ids

        non.outlier.ids <- setdiff(rownames(r$X), outlier.ids)
        sp.cags <- paste0(r$taxon, "__CAG_", cltr[non.outlier.ids])
        names(sp.cags) <- non.outlier.ids

        sp.cags <- c(sp.cags, sp.cags.ext)
        sp.cags <- sp.cags[rownames(r$X)]

    } else {
        sp.cags <- paste0(r$taxon, "__CAG_", cltr)
        names(sp.cags) <- names(cltr)
    }

    if ( verbose ) {
        elapsed.time(ptm)
    }

    list(r=r,
         sp.cags=sp.cags)
}


#' Identifies Neighboring Clusters for Each Cluster in a State Space Clusering.
#'
#' This function determines the neighboring clusters for each cluster in a given state space \( X \).
#' A cluster \( cl_j \) is defined as a neighbor of \( cl_i \) if there are at least \( \code{min.hits} \)
#' points in \( cl_i \) such that one of its \( K \) nearest neighbors belongs to \( cl_j \).
#'
#' @param cltr     A named vector representing the clustering in the state space \( X \). The names should match the row names of \( X \).
#' @param X        A matrix representing the state space, with each row corresponding to a point in the space.
#' @param K        An integer specifying the number of nearest neighbors to consider when identifying neighboring clusters. Must be positive.
#' @param min.hits An integer specifying the minimum number of points in \( cl_i \) that must have a nearest neighbor in \( cl_j \) for \( cl_j \) to be considered a neighbor. Must be non-negative.
#'
#' @return         Returns a list of frequency tables, each corresponding to a cluster \( cl_i \), indicating the number of elements in \( cl_i \) whose \( K \) nearest neighbors intersect with a neighbor cluster \( cl_j \).
#'
#' @examples
#' # Add illustrative examples to demonstrate the functionality of the function.
#'
#' @export
neighbor.cltrs <- function(cltr, X, K = 10, min.hits = 5) {

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if ( !all(names(cltr) == rownames(X)) ) {
        stop(paste("Mismatch between names(cltr) and rownames(X):",
                   length(setdiff(names(cltr), rownames(X))), "differences found."))
    }

    if ( as.integer(K) != K ) {
        stop("K has to be an integer")
    }

    if ( K < 1 ) {
        stop("K has to be at least 1")
    }

    if (as.integer(min.hits) != min.hits || min.hits < 0) {
        stop("min.hits has to be a non-negative integer")
    }

    nn <- get.knn(X, k = K)
    nn.i <- nn$nn.index
    rownames(nn.i) <- rownames(X)

    cl.ids <- names(table(cltr))
    cl.ids <- setdiff(cl.ids, 0)

    cltr.nbr <- list()
    for ( cl.id in cl.ids ) {
        idx <- cltr == cl.id
        ids <- names(cltr)[idx]
        cl.nn.i <- nn.i[ids,]
        ii <- as.vector(cl.nn.i)
        ii <- unique(ii)
        cltr.freq <- table(cltr[ii])
        cltr.freq <- cltr.freq[setdiff(names(cltr.freq[cltr.freq >= min.hits]), cl.id)]
        cltr.nbr[[cl.id]] <- cltr.freq
        ##cltr.nbr[[cl.id]] <- setdiff(names(cltr.freq[cltr.freq >= min.hits]), cl.id)
    }

    return(cltr.nbr)
}

## #' Re-embeds each BoA cluster and its neighbors into 11d space
## #'
## #' This routine is used to check for # disconnected components of a single cluster.
## #'
## #' @param id The ID of cltr.nbrs list identifying cluster that will be
## #'     re-embedded into dim-dimensional space together with its neighbors.
## #' @param X The original abundance matrix before low dimensional embedding.
## #' @param cltr A clustering of a low dimensional model of X.
## #' @param cltr.nbrs A list of cluster neighbors. An output from neighbor.cltrs().
## #' @param dim The dimension of the Euclidean space into which the the samples of
## #'     the given cluster (plus its neighbors) will be embedded. Default is 3.
## #' @param init.reticulate A logical value indicating whether reticulate pythyon
## #'     environment should be initiated. Default value: FALSE.
## #' @param exclude.cltrs A vector of cluster IDs that are to be excluded from the re-embedding.
## #'
## reembed.cltr <- function(id, X, cltr, cltr.nbrs, dim = 3, init.reticulate = FALSE, exclude.cltrs = NULL ) {

##     if (!is.matrix(X)) {
##         X <- try(as.matrix(X), silent = TRUE)
##         if (inherits(X, "try-error")) {
##             stop("X must be a matrix or coercible to a matrix")
##         }
##     }

##     if (!is.numeric(X)) {
##         stop("X must contain numeric values")
##     }

##     if (any(is.na(X)) || any(is.infinite(X))) {
##         stop("X cannot contain NA, NaN, or Inf values")
##     }

##     if ( nrow(X) != length(cltr) ) {
##         cat("nrow(X): ", nrow(X), "\n")
##         cat("length(cltr): ", length(cltr), "\n")
##         stop("The number of rows of X has to be the same as the length of the cluster vector")
##     }

##     if ( !all(rownames(X) == names(cltr)) ) {
##         stop("rownames(X) != names(clrt)")
##     }

##     n.cltrs <- length(table(cltr[cltr != 0]))

##     if ( length(cltr.nbrs) != n.cltrs ) {
##         cat("length(cltr.nbrs): ", length(cltr.nbrs), "\n")
##         cat("Number of non-zero clusters: ", n.cltrs, "\n")
##         stop("length(cltr.nbrs) != n.cltrs")
##     }

##     if (  !is.character(id) ) {
##         stop("id has to be a character string")
##     }

##     if ( !(id %in% names(cltr.nbrs)) ) {
##         cat("id: ", id, "\n")
##         stop("id is not a name of a component of the cltr.nbrs list.")
##     }

##     ## Identifying sample IDs of the specified cluster and its neighbors
##     if ( !is.null(cltr.nbrs) ) {
##         cl.nbrs <- names(cltr.nbrs[[id]])
##         cl.nbrs <- setdiff(cl.nbrs, 0)
##         cl.nbrs <- c(id, cl.nbrs)
##     } else {
##         cl.nbrs <- c(id)
##     }

##     if ( !is.null(exclude.cltrs) ) {
##         cl.nbrs <- setdiff(cl.nbrs, exclude.cltrs)
##     }

##     ## Collecting ids of all the clusters that will be re-embedded
##     ids <- c()
##     for ( cl in cl.nbrs ) {
##         idx <- cltr == cl
##         ids <- c(ids, names(cltr)[idx])
##     }

##     clX <- Rpacmap(X[ids,], ndim = dim, init.reticulate = init.reticulate)

##     return(clX)
## }


#' Generates a synthetic clustering in R^2 using Gaussian distributions
#'
#' @param n.clusters The number of clusters to generate
#' @param pts.per.cluster The number of points to generate per cluster; if NULL, number is chosen randomly for each cluster
#' @param min.pts The minimum number of points per cluster when pts.per.cluster is NULL
#' @param max.pts The maximum number of points per cluster when pts.per.cluster is NULL
#'
#' @return Returns a data frame with columns 'x', 'y' for coordinates and 'cluster' for the cluster label.
synthetic.Gaussian.2D.clustering <- function(n.clusters = 5, pts.per.cluster = NULL, min.pts = 50, max.pts = 150) {

  # Validate inputs
  if (!is.numeric(n.clusters) || n.clusters <= 0) stop("Invalid number of clusters.")
  if (!is.null(pts.per.cluster) && (pts.per.cluster <= 0)) stop("Invalid number of points per cluster.")
  if (is.null(pts.per.cluster) && (min.pts <= 0 || max.pts <= 0 || min.pts > max.pts)) stop("Invalid min.pts or max.pts.")

  # Initialize an empty data frame
  synthetic_data <- data.frame(x = numeric(0), y = numeric(0), cluster = character(0))

  for (i in 1:n.clusters) {
    # Decide the number of points for this cluster
    n_pts <- if (is.null(pts.per.cluster)) {
      sample(min.pts:max.pts, 1)
    } else {
      pts.per.cluster
    }

    # Generate random mean and covariance matrix for Gaussian distribution
    mean_vec <- runif(2, -10, 10)
    cov_matrix <- matrix(runif(4, 0.5, 2), nrow = 2)
    cov_matrix <- (cov_matrix + t(cov_matrix)) / 2  ## Make it symmetric

    # Make the matrix positive-definite if needed
    eigen_values <- eigen(cov_matrix)$values
    if (any(eigen_values <= 0)) {
      cov_matrix <- cov_matrix + (-min(eigen_values) + 0.1) * diag(2)
    }

    # Sample points from this Gaussian distribution
    cluster_data <- MASS::mvrnorm(n = n_pts, mu = mean_vec, Sigma = cov_matrix)
    cluster_df <- data.frame(cluster_data, cluster = as.character(i))
    colnames(cluster_df) <- c("x", "y", "cluster")

    # Append to the existing data
    synthetic_data <- rbind(synthetic_data, cluster_df)
  }

  return(synthetic_data)
}

## Generate synthetic clustering in R^2 using Gaussian distributions with controlled orientation
##
## @param n.clusters The number of clusters to generate
## @param pts.per.cluster The number of points to generate per cluster; if NULL, number is chosen randomly for each cluster
## @param min.pts The minimum number of points per cluster when pts.per.cluster is NULL
## @param max.pts The maximum number of points per cluster when pts.per.cluster is NULL
##
## @return Returns a data frame with columns 'x', 'y' for coordinates and 'cluster' for the cluster label.
controlled.synthetic.Gaussian.2D.clustering <- function(n.clusters = 5, pts.per.cluster = NULL, min.pts = 50, max.pts = 150) {

  synthetic_data <- data.frame(x = numeric(0), y = numeric(0), cluster = character(0))

  for (i in 1:n.clusters) {

    n_pts <- if (is.null(pts.per.cluster)) {
      sample(min.pts:max.pts, 1)
    } else {
      pts.per.cluster
    }

    mean_vec <- runif(2, -10, 10)

    ## Create controlled eigenvalues and eigenvectors
    eigen_values <- runif(2, 0.5, 2)
    angle <- runif(1, 0, 2 * pi)
    eigenvectors <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)

    ## Create the covariance matrix
    cov_matrix <- eigenvectors %*% diag(eigen_values) %*% t(eigenvectors)

    cluster_data <- MASS::mvrnorm(n = n_pts, mu = mean_vec, Sigma = cov_matrix)
    cluster_df <- data.frame(cluster_data, cluster = as.character(i))
    colnames(cluster_df) <- c("x", "y", "cluster")

    synthetic_data <- rbind(synthetic_data, cluster_df)
  }

  return(synthetic_data)
}

## #' Performs an embedding and clustering of a set of points.
## #'
## #' If a factor F is not NULL, then the function attempts to find a clustering of
## #' X (or rather its low dimensional model) that is concordant with the levels of
## #' F, where the conconrdance suppose to have the following properties:
## #'
## #' 1) A single level of F can be split into several clusters.
## #'
## #' 2) A cluster should have high level of F-levels' purity. That is, it is OK
## #' for a cluster to have vast majority of the elements to come from one level of
## #' F. It would be nice for the purity to be normalized to have values between 0
## #' and 1 with 1 being the highest purity. Of course, a cluster of size 1 would
## #' always have perfect purity. This is why I need the following extra condition.
## #'
## #' 3) The algorithm needs to find the smallest number of clusters with the
## #' highest possible purity.
## #'
## #' When F is not NULL, the algorithm performs clustering of low dim model of X
## #' for different values of min.pts and then choses the one with the maximal
## #' value of the following objective function:
## #'
## #' \deqn{F(m) = \frac{1}{C} \sum_{i=1}^{C} Purity(i) * Size(i)^\lambda}
## #'
## #' where \eqn{Purity(i) = 1 - Ev(i)} with Ev(i) defined as the levels of F
## #' evenness within the i-th cluster, where the evenness is defined as the
## #' entropy of the within cluster levels proportions divided by the log of the
## #' number of levels within the cluster. In this way the purity is always
## #' between 0 and 1 with the value of 1 corresponding to the case of a single
## #' level within the given cluster.
## #'
## #' @param X             A set of points (in the form of a matrix) for which the embedding and clustdering will be peformed.
## #' @param dim           The dimension of the low dimensional model of X. The algorithm uses PaCMAP algorithm.
## #' @param min.pts       A vector of integers indicating the minimum cluster sizes to be evaluated. Default is NULL which means min/max prop parameters will be used to determine the values of min.pts.
## #' @param min.prop      The minimum proportion of nrow(X) to be used in test clusterings.
## #' @param max.prop      The maxinum proportion of nrow(X) to be used in test clusterings.
## #' @param n.test.cltr   The number of test clusters.
## #' @param soft.K        The number of nearest neighbors to consider when softening the cluster boundaries. Default is 20.
## #' @param F             A factor with which we want to the clustering to be as concordant as possible.
## #' @param lambda        A scaling factor of the objective function. Large values of lambda promote smaller number of clusters. Default is 1.
## #' @param n.cores       The number of cores to use.
## #' @param distance      A string indicating the distance to use in the dimensionality reduction algorithm. Default is 'angular'.
## #' @param verbose       A logical value indicating whether progress should be displayed. Default is FALSE.
## #'
## #' @return A list with the following components:
## #' \describe{
## #'   \item{cltr}{A vector representing the cluster assignment of each point in the low-dimensional space.}
## #'   \item{cltr.list}{A list of ids within each cluster.}
## #'   \item{cltrgs}{A list containing cluster assignments obtained for different values of `min.pts` if the factor `F` is not NULL. Otherwise, this will be NULL.}
## #'   \item{obj.values}{A numeric vector containing the values of the objective function for each clustering tested when `F` is not NULL. Otherwise, this will be NULL.}
## #'   \item{opt.i}{The index corresponding to the optimal clustering according to the objective function when `F` is not NULL. Otherwise, this will be NULL.}
## #' }
## embedd.and.cltr <- function(X,
##                            dim = 3,
##                            min.prop = 0.1,
##                            max.prop = 0.5,
##                            n.test.cltr = 10,
##                            min.pts = NULL,
##                            soft.K = 20,
##                            F = NULL,
##                            lambda = 1,
##                            n.cores = 10,
##                            distance = "angular",
##                            verbose = FALSE)
## {
##     if ( !requireNamespace("clValid", quietly = TRUE) ) {
##         stop("The clValid package is not installed. Please install it before using this routine.")
##     }

##     if ( !requireNamespace("dbscan", quietly = TRUE) ) {
##         stop("The dbscan package is not installed. Please install it before using this routine.")
##     }

##     if ( lambda <= 0 ) {
##         stop("lambda needs to be a positive number.")
##     }

##     if ( is.null(min.pts) ) {

##         stopifnot( is.numeric(min.prop) && 0 < min.prop && min.prop < 1 )
##         stopifnot( is.numeric(max.prop) && 0 < max.prop && max.prop < 1 )
##         stopifnot( as.integer(n.test.cltr) == n.test.cltr && n.test.cltr > 3 )

##         min.pts <- seq(min.prop * nrow(X), max.prop * nrow(X), length.out = n.test.cltr)
##         min.pts <- as.integer(min.pts)

##         if ( verbose ) {
##             cat("min.pts was set to: ")
##             cat(min.pts)
##             cat("\n")
##         }
##     }

##     if ( length(min.pts) < n.cores ) {
##         n.cores <- length(min.pts)
##         if ( verbose ) {
##             cat("WARNING:  n.cores was set to ", n.cores, "\n")
##         }
##     }

##     if ( verbose ) {
##         fn.ptm <- proc.time()
##         cat("Generating a low dimensional embedding ... \n")
##         ptm <- proc.time()
##     }

##     X.dim <- Rpacmap(X, ndim = dim, distance = distance, verbose = FALSE)

##     if ( verbose ) {
##         elapsed.time(ptm)
##         cat("Clustering the low dimensional embedding ... \n")
##         ptm <- proc.time()
##     }

##     cltrg.res <- optimize.min.pts(X.dim, min.pts, soft.K, F, lambda, n.cores)

##     if ( verbose ) {
##         elapsed.time(ptm)
##         txt <- "Total elapsed time:"
##         elapsed.time(fn.ptm, txt, with.brackets = FALSE)
##     }

##     res <- list(X.dim = X.dim,
##                cltr = cltrg.res$cltr,
##                cltr.list = cltrg.res$cltr.list,
##                cltrgs = cltrg.res$cltrg.res,
##                obj.values = cltrg.res$obj.values,
##                opt.i = cltrg.res$opt.i)

##     class(res) <- "embedAndCltr"

##     return(res)
## }

#' Optimizes min.pts based on the objective function
#'
#' @param X.dim         A set of points (in the form of a matrix) for which the clustering will be peformed.
#' @param min.pts       A vector of integers indicating the minimum cluster sizes to be evaluated. Default is NULL which means min/max prop parameters will be used to determine the values of min.pts.
#' @param soft.K        The number of nearest neighbors to consider when softening the cluster boundaries. Default is 20.
#' @param F             A factor with which we want to the clustering to be as concordant as possible.
#' @param lambda        A scaling factor of the objective function. Large values of lambda promote smaller number of clusters. Default is 1.
#' @param n.cores       The number of cores to use.
#' @param verbose       A logical value indicating whether progress should be displayed. Default is FALSE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{cltr}{A vector representing the cluster assignment of each point in the low-dimensional space.}
#'   \item{cltr.list}{A list of ids within each cluster.}
#'   \item{cltrg.res}{A list containing cluster assignments obtained for different values of `min.pts` if the factor `F` is not NULL. Otherwise, this will be NULL.}
#'   \item{obj.values}{A numeric vector containing the values of the objective function for each clustering tested when `F` is not NULL. Otherwise, this will be NULL.}
#'   \item{opt.i}{The index corresponding to the optimal clustering according to the objective function when `F` is not NULL. Otherwise, this will be NULL.}
#' }'
optimize.min.pts <- function(X.dim, min.pts, soft.K = 20, F = NULL, lambda = 1, n.cores = 10) {

    cltrg.res <- NULL
    obj.values <- NULL
    opt.i <- NULL

    if ( is.null(F) ) {

        cltr.res <- hdbscan.cltr(X.dim,
                                method = "dunn",
                                min.pts = min.pts,
                                soft.K = soft.K,
                                n.cores = n.cores,
                                verbose = FALSE)

        if ( !is.null(cltr.res$dunn.cltr.ext) ) {
            cltr <- cltr.res$dunn.cltr.ext
        } else {
            cltr <- cltr.res$dunn.cltr
        }

    } else {

        if ( !is.factor(F) ) {
            F <- as.factor(F)
            cat("WARNING: F was turned into a factor!")
        }

        if ( length(levels(F)) > nrow(X) / 2 ) {
            cat("n(levels(F)):", length(levels(F)), "\n")
            cat("nrow(X):", nrow(X), "\n")
            cat("WARNING: The number of levels of F is greater than nrow(X) / 2.")
        }

        ## Clustering function
        cltr.fn <- function(k, X.dim) {

            cltr <- dbscan::hdbscan(X.dim, minPts = k)$cluster
            names(cltr) <- rownames(X.dim)

            return(cltr)
        }

        registerDoParallel(n.cores)
        ptm <- proc.time()
        cltrg.res <- foreach ( k = min.pts ) %dopar% {
            cltr.fn(k, X.dim)
        }
        elapsed.time(ptm)
        stopImplicitCluster()

        ## Defining an objective function
        ## \[ F(m) = \frac{1}{C} \sum_{i=1}^{C}  Purity(i) \cdot Size(i)^\lambda \]
        obj.fn <- function(cltrs, F, lambda) {

            uq.cltrs <- names(table(cltrs[cltrs != 0]))
            n.cltrs <- length(uq.cltrs)
            weighted.purity <- numeric(n.cltrs)

            for ( i in seq(n.cltrs) ) {
                cltr <- uq.cltrs[i]
                idx <- cltrs == cltr
                ##
                cltr.F.freq <- table(F[idx])        ## frequency of F levels within the given cluster
                p <- cltr.F.freq / sum(cltr.F.freq) ## proportion of F levels within the given cluster
                ##
                cltr.size <- sum(idx)
                ##
                weighted.purity[i] <- purity(p) * cltr.size^lambda
            }

            sum(weighted.purity) / n.cltrs
        }

        obj.values <- numeric(length(cltrg.res))
        for ( i in seq(length(cltrg.res)) ) {
            cltrs <- cltrg.res[[i]]
            cltrs.freq <- table(cltrs)
            if ( '0' %in% names(cltrs.freq) ) {
                cltrs <- kNN.cltr.imputation(X.dim, cltrs, K = soft.K)
            }
            obj.values[i] <- obj.fn(cltrs, F, lambda)
        }

        opt.i <- which.max(obj.values)
        cltr <- cltrg.res[[opt.i]]
        cltr[cltr != 0] <- clusters.reorder(cltr[cltr != 0])
    }

    ## Creating a list of ids within each cluster
    uq.cltrs <- names(table(cltr))
    n.cltrs <- length(uq.cltrs)
    cltr.list <- list()
    for ( i in seq(n.cltrs) ) {
        idx <- cltr == uq.cltrs[i]
        cltr.list[[as.character(uq.cltrs[i])]] <- names(cltr)[idx]
    }

    list(cltr = cltr,
         cltr.list = cltr.list,
         cltrg.res = cltrg.res,
         obj.values = obj.values,
         opt.i = opt.i)
}


## #' Generates a list of deeper and deeper clusterings of a state space.
## #'
## #' The function first clusters a low dimensional model, X.dim, of X. This is
## #' level 1 clustering. At subsequent level, the function embeds and clusters X
## #' restricted to samples of each cluster from the previous level. The
## #' sub-clustering stops either if the size of the given cluster goes below the
## #' threshold of 'cltr.size.thld' or the reult of clustering at some level is
## #' only a single non-zero cluster.
## #'
## #' @param X              A source data table from which Y was constructed.
## #' @param X.dim          A low dimensional model of X.
## #' @param n.levels       The number of levels of clustering of the points of X.
## #' @param cltr.size.thld The cluster size threshold. The clustering stops when the largest cluster size is smaller than cltr.size.thld or the resulting space has only one cluster. Default is 100.
## #' @param min.pts        A vector of minimal cluster size values for hdbscan.cltr(). Default is 5*(1:10).
## #' @param min.prop       The proportion of the number of samples in the current maximal cluster that will be taken as the min.pts in HDBSCAN's k parameter. More precisely, k = min( c(cltr.size.thld, min.prop * nrow(max.cltr.X.dim)) ). Default is 0.1.
## #' @param max.prop       The maxinum proportion of nrow(X) to be used in test clusterings.
## #' @param n.test.cltr    The number of test clusters.
## #' @param soft.K         The number of nearest neighbors to consider when softening the cluster boundaries. Default is 20.
## #' @param distance       The distance in the PaCMAP algorithm. Default is 'angular'.
## #' @param dim            The dimenion into which \code{X\[max cltr samples,\]} is embedded. Default is 3.
## #' @param n.cores        The number of cores to use in hdbscan.cltr(). Default is 10.
## #' @param verbose        A logical value indicating whether detailed output should be printed during the function's execution. If TRUE, the function will print out additional details about the calculations it's performing. Default is TRUE.
## #'
## #' ToDo: add a possibility to call local intrinsic dim estimator and determine dimension of the embedding that way (possibly having the dimension to be different for different clusters).
## ehclust <- function(X,
##                    X.dim = NULL,
##                    n.levels = 3,
##                    cltr.size.thld = 50,
##                    min.pts = 5*(1:10),
##                    min.prop = 0.1,
##                    max.prop = 0.5,
##                    n.test.cltr = 10,
##                    soft.K = 20,
##                    distance = "angular",
##                    dim = 3,
##                    n.cores = 10,
##                    verbose = TRUE) {

##     if (!is.matrix(X)) {
##         X <- try(as.matrix(X), silent = TRUE)
##         if (inherits(X, "try-error")) {
##             stop("X must be a matrix or coercible to a matrix")
##         }
##     }

##     if (!is.numeric(X)) {
##         stop("X must contain numeric values")
##     }

##     if (any(is.na(X)) || any(is.infinite(X))) {
##         stop("X cannot contain NA, NaN, or Inf values")
##     }

##     if ( !is.null(X.dim) && !(is.matrix(X.dim) || is.data.frame(X.dim)) ) {
##         stop("X.dim has to be a matrix or data frame")
##     }

##     if ( !is.null(X.dim) && nrow(X) != nrow(X.dim) ) {
##         stop(paste("The number of rows in X and X.dim must be the same. Currently, nrow(X):", nrow(X), "and nrow(X.dim):", nrow(X.dim), "."))
##     }

##     if ( !requireNamespace("clValid", quietly = TRUE) ) {
##         stop("The clValid package is not installed. Please install it before using this routine.")
##     }

##     if ( !requireNamespace("dbscan", quietly = TRUE) ) {
##         stop("The dbscan package is not installed. Please install it before using this routine.")
##     }

##     if ( is.null(min.pts) ) {

##         stopifnot( is.numeric(min.prop) && 0 < min.prop && min.prop < 1 )
##         stopifnot( is.numeric(max.prop) && 0 < max.prop && max.prop < 1 )
##         stopifnot( as.integer(n.test.cltr) == n.test.cltr && n.test.cltr > 3 )

##         min.pts <- seq(min.prop * nrow(X), max.prop * nrow(X), length.out = n.test.cltr)
##         min.pts <- as.integer(min.pts)

##         if ( verbose ) {
##             cat("min.pts was set to: ")
##             cat(min.pts)
##             cat("\n")
##         }
##     }

##     if ( is.null(rownames(X)) ) {
##         rownames(X) <- seq(nrow(X))
##     }

##     ##
##     ## Performing dimensionality reduction if X.dim is NULL
##     ##
##     if ( is.null(X.dim) ) {
##         X.dim <- Rpacmap(X, ndim = dim, distance = distance)
##     }

##     ##
##     ## Clustering of X.dim
##     ##
##     if ( verbose ) {
##         fn.ptm <- proc.time()
##         cat("Clustering of X.dim the low dimensional embedding ... ")
##         ptm <- proc.time()
##     }

##     cltr.res <- hdbscan.cltr(X.dim,
##                              method = "dunn",
##                              min.pts = min.pts,
##                              soft.K = soft.K,
##                              n.cores = n.cores,
##                              verbose = FALSE)

##     if ( !is.null(cltr.res$dunn.cltr.ext) ) {
##         cltr <- cltr.res$dunn.cltr.ext
##     } else {
##         cltr <- cltr.res$dunn.cltr
##     }

##     ##names(cltr) <- rownames(X)
##     cltr[cltr != 0] <- clusters.reorder(cltr[cltr != 0])

##     ## Turning vector, cltr, of cluster IDs over indices of X into a list with
##     ## components: cluster IDs and values lists of row of X ids
##     cltr.list <- list()
##     for ( id in names(table(cltr)) ) {
##         ids <- names(cltr[cltr == id])
##         to.be.cltred <- length(ids) > cltr.size.thld
##         cltr.list[[id]] <- list("ids" = ids,
##                                 "to.be.cltred" = to.be.cltred) # this flag is used to stop
##                                         # 'embed-and-cluster'
##                                         # process on the given
##                                         # cluster either because the
##                                         # size of cluster is less
##                                         # than the threshold or
##                                         # 'embed-and-cluster'
##                                         # produces only one non-zero
##                                         # cluster
##     }

##     if ( verbose ) elapsed.time(ptm)

##     ## A list of embeddings
##     X.dim.level <- list()
##     X.dim.level[[1]] <- list() # at ech level, it is a list of embeddings corresponding to clusters of X
##     X.dim.level[[1]][[1]] <- X.dim

##     ## clustering at level 1
##     cltrg.level <- list()
##     cltrg.level[[1]] <- cltr

##     ## cltr list is level 1 clustering
##     cltr.list.level <- list()
##     cltr.list.level[[1]] <- cltr.list

##     ##
##     Ecltr.list <- list()
##     Ecltr.list[[1]] <- NULL

##     if ( n.levels > 1 ) {
##         for ( level in 2:n.levels ) {

##             if ( verbose ) {
##                 cat("Creating level",level,"clusters\n")
##             }

##             X.dim.level[[level]] <- list()
##             ##
##             cltr.list.level[[level]] <- list()
##             level.cltr.list <- cltr.list.level[[level - 1]]

##             if ( "0" %in% names(level.cltr.list) ) {
##                 level.cltr <- level.cltr.list[["0"]]
##                 cltr.list <- list()
##                 cltr.list[["0"]] <- level.cltr
##                 cltr.list.level[[level]] <- c(cltr.list.level[[level]], cltr.list)
##             }

##             ##
##             Ecltr.list[[level]] <- list()

##             for ( id in setdiff(names(level.cltr.list), "0") ) { # sub-clustering all clusters, except cluster "0"

##                 if ( verbose ) {
##                     cat("Clustering", id, "... ")
##                     ptm <- proc.time()
##                 }

##                 level.cltr <- level.cltr.list[[id]]
##                 ids <- level.cltr$ids ## sample ids of the given cluster

##                 if ( level.cltr$to.be.cltred ) {

##                     cltr.res <- embedd.and.cltr(X[ids,],
##                                                 dim = dim,
##                                                 min.pts = min.pts,
##                                                 soft.K = soft.K,
##                                                 distance = distance,
##                                                 n.cores = n.cores,
##                                                 verbose = FALSE)

##                     X.dim.level[[level]][[id]] <- cltr.res$X.dim
##                     cltr.list <- cltr.res$cltr.list

##                     Ecltr.list[[level]][[id]] <- cltr.res

##                     ## Turning cltr.list into a list of lists, each list carrying
##                     ## cluster IDs and also a component 'to.be.cltred' indicating if
##                     ## this cluster needs to be further subclustered
##                     for ( name in names(cltr.list) ) {
##                         ids <- cltr.list[[name]]
##                         cltr.list[[name]] <- NULL
##                         new.name <- paste0(id,":",name)
##                         to.be.cltred <- TRUE
##                         if ( length(ids) < cltr.size.thld ||
##                              (length(cltr.list) == 1 && names(cltr.list) == c("1")) ||
##                              (length(cltr.list) == 2 && all(names(cltr.list) == c("0","1"))) ) {
##                             to.be.cltred <- FALSE
##                         }
##                         cltr.list[[new.name]] <- list("ids" = ids,
##                                                       "to.be.cltred" = to.be.cltred) # set to FALSE if names(cltr.list) are {"1"} or {"0", "1"} or length(ids) < thld
##                     }

##                 } else {

##                     X.dim.level[[level]][[id]] <- X.dim.level[[level - 1]][[id]]
##                     cltr.list <- list()
##                     cltr.list[[id]] <- level.cltr
##                 }

##                 cltr.list.level[[level]] <- c(cltr.list.level[[level]], cltr.list)

##                 if ( verbose ) elapsed.time(ptm)

##             } # END OF for ( id in setdiff(names(level.cltr.list), "0") )

##             ## Turning cltr.list.level[[level]] into a clustering vector
##             cltr <- c()
##             cltr.level <- cltr.list.level[[level]]
##             for ( name in names(cltr.level) ) {
##                 cl <- cltr.level[[name]]
##                 cltr[cl$ids] <- name
##             }

##             cltrg.level[[level]] <- cltr

##         } # END OF for ( level in 2:n.levels ) {
##     }

##     if ( verbose ) {
##         elapsed.time(ptm)
##         txt <- "Total elapsed time:"
##         elapsed.time(fn.ptm, txt, with.brackets = FALSE)
##     }

##     res <- list(X.dim.level = X.dim.level,
##                 cltrg.level = cltrg.level,
##                 cltr.list.level = cltr.list.level,
##                 Ecltr.list = Ecltr.list,
##                 ## selected input parameters of the function
##                 cltr.size.thld = cltr.size.thld,
##                 min.pts = min.pts,
##                 soft.K = soft.K,
##                 distance = distance,
##                 dim = dim)

##     class(res) <- "ehclust"

##     return(res)
## }

#' Plots in 3d the results of ehclust
#'
#' @param obj A result of ehclust().
#' @param level A level to be plotted.
#' @param show.cltr.labels Set to TRUE to show cluster labels. Default is TRUE.
#' @param ... A placeholder for parameters of plot3d.cltrs().
ehclust.plot <- function(obj, level = 1, show.cltr.labels = TRUE, ...) {

    if ( as.integer(level) != level ) {
        stop("level has to be an integer")
    }

    if ( level < 1 ) {
        stop("level has to be an integer greater than 0")
    }

    if ( level > length(obj$cltr.list.level) ) {
        cat("level:", level, "\n")
        cat("length(obj$cltr.list.level):", length(obj$cltr.list.level), "\n")
        stop("level has to be an integer no greater than obj$cltr.list.level")
    }

    X.dim <- obj$X.dim.level[[1]][[1]]
    cltr <- obj$cltrg.level[[level]]
    cltr <- cltr[rownames(X.dim)]
    plot3d.cltrs(X.dim, cltr, show.cltr.labels = show.cltr.labels, ...)
}

#' Taxonomic purity of a gene clustering
#'
#' @param cltr   A clustering.
#' @param n.min  The minimal number genes with the given species assignment for the species to be a part of the species-clltr-frequency table.
gene.cltr.tx.purity <- function(cltr, n.min = 20) {

    if ( !exists("gene.tx.map", envir = environment()) ) {
        file <- "~/projects/rllm-projects/projects/ZB/data/virgo2/VIRGO2_taxon_map.rda"
        ## save(gene.tx.map, file = file)
        load(file)
    }

    ## Extracting gene taxonomic annotation for clustered genes
    gene.tx <- gene.tx.map[names(cltr)]
    has.tx <- !is.na(gene.tx)
    gene.tx <- gene.tx[has.tx]

    ## Restricting to genes with species taxonomy
    gene.with.tx <- intersect(names(gene.tx), names(cltr))
    cltr <- cltr[gene.with.tx]
    gene.tx <- gene.tx[gene.with.tx]

    ## species-clustering frequency table
    tx.cltr.freq <- table(gene.tx, cltr)

    ## Restricting to species with at least n.min genes
    spp.n.det <- rowSums(tx.cltr.freq)
    idx <- spp.n.det >= n.min
    tx.cltr.freq <- tx.cltr.freq[idx,]
    core.tx.cltr.freq <- tx.cltr.freq

    ## Purity of species assignment
    i <- grep("MultiGenera", rownames(tx.cltr.freq))
    tx.purity <- apply(tx.cltr.freq[-i,], 2, purity)

    ## Adding purity as the last column
    ## tx.cltr.freq <- rbind(tx.cltr.freq, signif(tx.purity, digits = 2))
    ## tx.cltr.freq <- apply(tx.cltr.freq, 2, as.character)
    ## rownames(tx.cltr.freq)[nrow(tx.cltr.freq)] <- "Taxon Purity"
    ## rownames(tx.cltr.freq)[1:(nrow(tx.cltr.freq)-1)] <- rownames(core.tx.cltr.freq)
    tx.cltr.freq <- add.row.to.matrix(tx.cltr.freq, signif(tx.purity, digits = 2))
    rownames(tx.cltr.freq)[nrow(tx.cltr.freq)] <- "Taxon Purity"

    list(tx.cltr.freq = tx.cltr.freq,
         core.tx.cltr.freq = core.tx.cltr.freq,
         tx.purity = tx.purity,
         gene.tx = gene.tx)
}

#' Species Cluster Frequency Analysis
#'
#' @description
#' Generic function for generating species-cluster frequency tables from clustering results.
#'
#' @param object A clustering object (e.g., from ehclust).
#' @param ... Additional arguments passed to methods.
#' @export
spp.cltr.freq <- function(object, ...) {
  UseMethod("spp.cltr.freq")
}

#' Species Cluster Frequency Table for ehclust Objects
#'
#' @description
#' Returns species-cluster frequency table for the clustering at given level listing
#' only species detected in at least n.min genes.
#'
#' @param object A result of ehclust().
#' @param level The level of clustering for which the species-cluster-frequency table
#'        will be created. Default is 1.
#' @param base The base of the logarithm in the entropy formula. Default is 50.
#' @param n.min The minimal number of genes with the given species assignment for
#'        the species to be a part of the species-cluster-frequency table. Default is 20.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing:
#'   \item{tx.freq}{Frequency table of taxonomic assignments}
#'   \item{spp.cltr.freq}{Species-cluster frequency table}
#'   \item{tx.purity}{Purity scores for each cluster}
#'   \item{median.tx.purity}{Median purity across all clusters}
#'   \item{spp.cltr.freq.with.purity}{Combined table with purity scores}
#'   \item{gene.tx}{Gene to taxonomy mapping}
#'
#' @export
#' @rdname spp.cltr.freq
spp.cltr.freq.ehclust <- function(object, level = 1, base = 50, n.min = 20, ...) {
    if (!inherits(object, "ehclust")) {
        stop("This function operates on ehclust objects only!")
    }

    # Check if gene.tx.map exists in the current environment
    if (!exists("gene.tx.map", envir = parent.frame())) {
        # Try to find it in the package environment or load from file
        if (exists("gene.tx.map", envir = .GlobalEnv)) {
            gene.tx.map <- get("gene.tx.map", envir = .GlobalEnv)
        } else {
            file <- "~/projects/rllm-projects/projects/ZB/data/virgo2/VIRGO2_taxon_map.rda"
            if (file.exists(file)) {
                load(file, envir = environment())
            } else {
                stop("gene.tx.map not found and data file does not exist: ", file)
            }
        }
    }

    ## Extracting a clustering at the given level
    if (length(object$cltrg.level) < level) {
        stop("Requested level ", level, " not available. Maximum level is ",
             length(object$cltrg.level))
    }
    cltr <- object$cltrg.level[[level]]

    ## Extracting gene taxonomic annotation for clustered genes
    gene.tx <- gene.tx.map[names(cltr)]
    has.tx <- !is.na(gene.tx)
    gene.tx <- gene.tx[has.tx]

    if (length(gene.tx) == 0) {
        warning("No taxonomic annotations found for clustered genes")
        return(list(tx.freq = NULL,
                   spp.cltr.freq = NULL,
                   tx.purity = NULL,
                   median.tx.purity = NA,
                   spp.cltr.freq.with.purity = NULL,
                   gene.tx = NULL))
    }

    tx.freq <- sort(table(gene.tx), decreasing = TRUE)
    idx <- tx.freq >= n.min
    tx.freq <- tx.freq[idx]
    tx.freq.df <- cbind(tx.freq)
    tx.freq.df <- as.data.frame(tx.freq.df)
    colnames(tx.freq.df)[1] <- "Count"

    ## Restricting to genes with species taxonomy
    gene.with.tx <- intersect(names(gene.tx), names(cltr))
    cltr <- cltr[gene.with.tx]
    gene.tx <- gene.tx[gene.with.tx]

    ## species-clustering frequency table
    tx.cltr.freq <- table(gene.tx, cltr)

    ## Restricting to species with at least n.min genes
    spp.n.det <- rowSums(tx.cltr.freq)
    idx <- spp.n.det >= n.min
    tx.cltr.freq <- tx.cltr.freq[idx,]

    ## Purity of species assignment
    i <- grep("MultiGenera", rownames(tx.cltr.freq))
    if (length(i) > 0) {
        tx.purity <- apply(tx.cltr.freq[-i,], 2, purity, base = base)
    } else {
        tx.purity <- apply(tx.cltr.freq, 2, purity, base = base)
    }

    ## Adding purity as the last column
    spp.cltr.freq.with.purity <- rbind(tx.cltr.freq, signif(tx.purity, digits = 2))
    rownames(spp.cltr.freq.with.purity)[nrow(spp.cltr.freq.with.purity)] <- "Taxon Purity"

    list(tx.freq = tx.freq,
         spp.cltr.freq = tx.cltr.freq,
         tx.purity = tx.purity,
         median.tx.purity = median(tx.purity),
         spp.cltr.freq.with.purity = spp.cltr.freq.with.purity,
         gene.tx = gene.tx)
}

#' Performes EH-clustering of genes based on the similarity of their cross-distance basin matrices.
#'
#' First, we construct the distance matrix between cross.dist's of the ids, where the distance between two matrices is
#' the Euclidean distance between the corresponding vectors. Then multi-dimensional scaling is used to embedd the genes
#' in 100 dimensional space using the given distance matrix. Finally, a 3D embedding is created using PaCMAP and HDBSCAN
#' is used to cluster it. The process is iterated n.levels time in effect constructing an EH-clustering.
#'
#' @param ids Gene IDs.
#' @param cross.dist A list whose components are numerical matrices (no missing values allowed) of the same shape. That
#'     is, the same number of rows and columns.
#' @param n.levels The number of levels of EH-clustering.
#' @param with.cmdscaling Set to TRUE to use multi-dimensional scaling approach.
#'
#' @return The function returns an object of class ehclust with the EH-clustering of the genes.
ehclust.cross.dist <- function(ids, cross.dist, n.levels = 2, with.cmdscaling = FALSE) {

    n.ids <- length(ids)

    if ( with.cmdscaling ) {
        ## Constructing the distance matrix between cross.dist's of the ids.
        ## The distance between two matrices is the Euclidean distance between the corresponding vectors.
        dist.mat <- matrix(0, nrow = n.ids, ncol = n.ids)
        for ( i in 1:(n.ids - 1) ) {
            for ( j in (i + 1):n.ids ) {
                M1 <- cross.dist[[ids[i]]]
                M2 <- cross.dist[[ids[j]]]
                dist.mat[i,j] <- L2.norm(as.numeric(M1) - as.numeric(M2))
                dist.mat[j,i] <- dist.mat[i,j]
            }
        }
        ## Creating a 100d embedding using Multi-dimensional scaling
        X <- cmdscale(as.dist(dist.mat), eig = FALSE, k = 100)

    } else {

        M <- cross.dist[[ids[1]]]
        X <- matrix(0, nrow = n.ids, ncol = nrow(M) * ncol(M))
        for ( i in seq(n.ids) ) {
            M <- cross.dist[[ids[i]]]
            x <- as.numeric(M)
            min.log.x <- min(log(x[x > 0]))
            x[x > 0] <- log(x[x > 0])
            x[x == 0] <- min.log.x
            X[i,] <- x
        }
    }

    ## Creating a ehclust of X
    res <- ehclust(X, n.levels = n.levels)

    return(res)
}



#' Changes cluster assignment based on the majoriy vote over kNN's of each point.
#'
#' This function changes cluster assignment based on the majoriy vote over kNN's
#' of each point. For each point a purity index, pIdx(x), of the clustering at x
#' is computes as the proportion of elements of the given set, N(x,k), of kNN of
#' x that were assigned to the same cluster as x. All points are sorted in the
#' ascending order of the purity index of the clustering. For each point the
#' clustering assignment is reassessed based on the majority vote over its
#' kNN's.
#'
#' @param X        A matrix representing the state space, with each row corresponding to a point in the space.
#' @param cltr     A named vector representing the clustering in the state space \( X \). The names should match the row names of \( X \).
#' @param K        An integer specifying the number of nearest neighbors to consider when identifying neighboring clusters. Must be positive.
#'
#' @return         Returns a new clustering of X.
#'
knn.recltrg <- function(X, cltr, K = 10) {

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

    if (!is.vector(cltr)) {
        stop("'cltr' should be a vector.")
    }

    if ( nrow(X) != length(cltr) ) {
        cat("nrow(X):", nrow(X),"\n")
        cat("length(cltr):", length(cltr),"\n")
        stop("The length of cltr has to be the same as the number of rows of X")
    }

    ## Ensure 'k' is an integer
    k <- as.integer(K)
    ## if ( as.integer(K) != K ) {
    ##   stop("K has to be an integer")
    ## }

    if ( !all(names(cltr) == rownames(X)) ) {
        stop(paste("Mismatch between names(cltr) and rownames(X):",
                   length(setdiff(names(cltr), rownames(X))), "differences found."))
    }

    if ( K < 2 ) {
        stop("K has to be at least 2")
    }


    nn <- get.knn(X, k = k)
    nn.i <- nn$nn.index
    ##nn.d <- nn$nn.dist
    nn.i <- cbind(seq(nrow(X)), nn.i)
    ##nn.d <- cbind(rep(0,nrow(X)), nn.d)
    ##rownames(nn.i) <- rownames(X)

    nrX <- nrow(X)

    ## cluster purity
    k1 <- k + 1
    purity.idx <- numeric(nrX)
    for ( i in seq(nrX) ) {
        knn.cltr <- cltr[nn.i[i,]]
        ref.cltr <- cltr[i]
        knn.cltr <- ifelse(knn.cltr == ref.cltr, 1, 0)
        knn.cltr.freq <- table(knn.cltr)
        purity.idx[i] <- knn.cltr.freq["1"][[1]] / k1
    }

    o <- order(purity.idx)
    ##purity.idx <- purity.idx[o]
    cltr <- cltr[o,]
    X <- X[o,]

    nn <- get.knn(X, k = k)
    nn.i <- nn$nn.index
    nn.i <- cbind(seq(nrow(X)), nn.i)

    for ( i in seq(nrX) ) {
        knn.cltr <- cltr[nn.i[i,]]
        knn.cltr.freq <- sort(table(knn.cltr), decreasing = TRUE)
    }

    list(purity.idx = purity.idx)
}


#' Calculate Within-Cluster Means and Relative Means
#'
#' @description This function computes the mean values of a numeric variable within
#' each cluster defined by a grouping variable, and returns the results as a
#' data frame with absolute and relative means.
#'
#' @param y A numeric vector containing the values to be averaged.
#' @param cltr A character vector of the same length as \code{y}, containing
#' cluster identifiers.
#'
#' @return A data frame with three columns:
#' \itemize{
#'   \item \code{Cluster} - The unique cluster identifiers from the input \code{cltr}
#'   \item \code{Mean} - The mean value of \code{y} within each cluster
#'   \item \code{Relative_Mean} - The ratio of each cluster's mean to the overall mean of \code{y}
#' }
#'
#' @examples
#' # Simple example
#' y.values <- c(10, 15, 12, 8, 22, 14, 9, 18)
#' clusters <- c("A", "B", "A", "B", "C", "A", "C", "B")
#' cluster.means(y.values, clusters)
#'
#' @details The function performs validation to ensure that \code{y} is numeric and
#' that both input vectors have the same length. It uses the \code{aggregate()}
#' function internally to calculate cluster means. The Relative.Mean column shows
#' how each cluster's mean compares to the overall mean of all values.
#'
#' @export
cluster.means <- function(y, cltr) {
    ## Input validation
    if(!is.numeric(y)) {
        stop("y must be a numeric vector")
    }
    if(length(y) != length(cltr)) {
        stop("y and cltr must have the same length")
    }

    ## Calculate overall mean
    overall_mean <- mean(y, na.rm = TRUE)

    ## Create a data frame with the input vectors
    data <- data.frame(cluster = cltr, value = y)

    ## Compute mean for each cluster
    result <- aggregate(value ~ cluster, data = data, FUN = mean)

    ## Add relative mean column
    result$Relative_Mean <- result$value / overall_mean

    ## Rename columns for clarity
    colnames(result) <- c("Cluster", "Mean", "Relative_Mean")

    return(result)
}
