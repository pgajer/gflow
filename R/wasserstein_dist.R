##
## Implementations of Wasserstein distance and Wasserstein distance test for two samples.
##

#' Calculate Wasserstein Distance Between 1D Samples
#'
#' This function computes the Wasserstein distance (also known as Earth Mover's Distance)
#' between two samples from one-dimensional distributions.
#'
#' @param x A numeric vector representing a sample from the first distribution.
#' @param y A numeric vector representing a sample from the second distribution.
#'
#' @return A single numeric value representing the Wasserstein distance between the two samples.
#'
#' @details
#' The Wasserstein distance is calculated by sorting both input samples and computing
#' the average absolute difference between the sorted values. This implementation uses
#' a C function for efficient computation.
#'
#' Both input vectors must have the same length, contain only numeric values, and all
#' elements must be finite.
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000, mean = 1)
#' dist <- wasserstein.distance.1D(x, y)
#' print(dist)
#'
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Wasserstein_metric} for more information on the Wasserstein distance.
#'
#' @export
wasserstein.distance.1D <- function(x, y)
{
    # Check if x and y are numeric
    if (!is.numeric(x)) stop("x must be a numeric vector")
    if (!is.numeric(y)) stop("y must be a numeric vector")

    # Check for finite values
    if (!all(is.finite(x))) stop("All elements in x must be finite")
    if (!all(is.finite(y))) stop("All elements in y must be finite")

    # Check for equal length
    n <- length(x)
    if (length(y) != n) stop("x and y must have the same length")

    W <- 0
    out <- .C(C_wasserstein_distance_1D,
             as.double(x),
             as.double(y),
             as.integer(n),
             W=as.double(W))
    out$W
}

#' Performs Permutation Test for the Wasserstein Distance Between Two 1D Samples
#'
#' Performs a permutation test to assess the significance of the Wasserstein distance
#' between two 1D samples. The test determines if the observed Wasserstein distance
#' is significantly different from the distances obtained by permuting the samples,
#' under the null hypothesis that the samples come from the same distribution.
#' The function supports parallel computation for faster processing.
#'
#' @param x A numeric vector representing the first 1D sample.
#' @param y A numeric vector representing the second 1D sample.
#' @param n.perms An integer specifying the number of permutations to use for the
#'   permutation test. Default value is 10,000.
#' @param n.cores An integer specifying the number of cores to use for parallel
#'   computation. If `n.cores` is greater than 1, the null Wasserstein distances are
#'   computed in parallel. If `n.cores` is 1 (default), the computation is done serially.
#' @return A list containing:
#'   \itemize{
#'     \item \code{wasserstein1d}: The Wasserstein distance between `x` and `y`.
#'     \item \code{null.wasserstein1d}: A numeric vector of Wasserstein distances between
#'       permuted `x` and `y`, obtained from `n.perms` permutations.
#'     \item \code{p.value}: The p-value of the permutation test, representing the proportion
#'       of permuted Wasserstein distances that are greater than or equal to the observed distance.
#'   }
#' @details The function computes the actual Wasserstein distance between the two provided
#'   samples, then performs `n.perms` permutations of the concatenated samples and computes
#'   the Wasserstein distance for each permutation. The p-value is computed as the proportion
#'   of permuted distances that are at least as large as the actual distance.
#' @seealso \code{\link[transport]{wasserstein1d}} for computing the Wasserstein distance between two 1D samples.
#' @importFrom transport wasserstein1d
#' @examples
#' ## Two samples from the same distribution
#' x <- rnorm(100); y <- rnorm(100)
#' res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 1)
#' str(res)
#'
#' ## Two samples from different distributions
#' x <- rnorm(100)
#' y <- rnorm(100, mean = 2)
#' ## Perform the permutation test with 1000 permutations, serially
#' res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 1)
#' str(res)
#' ## Perform the permutation test with 1000 permutations, using 2 cores
#' res <- wasserstein1d.test(x, y, n.perms = 1000, n.cores = 2)
#' str(res)
#' @export
wasserstein1d.test <- function(x, y, n.perms = 10000, n.cores = 7)
{
  ## --- helpers -------------------------------------------------------------
  choose_cores <- function(requested) {
    if (!is.numeric(requested) || length(requested) != 1L ||
        is.na(requested) || !is.finite(requested)) {
      requested <- 1L
    }
    requested <- max(1L, as.integer(requested))

    ## Respect CRAN check limiter
    check_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(check_limit)) {
      ## CRAN asks you to limit to at most 2 in checks
      requested <- min(requested, 2L)
    }

    ## Respect getOption("mc.cores") if set; otherwise detect
    opt_mc <- getOption("mc.cores")
    if (!is.null(opt_mc) && is.numeric(opt_mc) && length(opt_mc) == 1L &&
        is.finite(opt_mc) && !is.na(opt_mc)) {
      requested <- min(requested, max(1L, as.integer(opt_mc)))
    }

    dc <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) 1L)
    if (!is.numeric(dc) || length(dc) != 1L || is.na(dc) || !is.finite(dc) || dc < 1L) {
      dc <- 1L
    }
    requested <- min(requested, max(1L, dc))

    max(1L, as.integer(requested))
  }

  use_fork <- function() {
    ## mclapply (fork) is only on Unix-like and is disabled on macOS + RStudio with certain settings;
    ## a safe heuristic is .Platform$OS.type == "unix"
    .Platform$OS.type == "unix"
  }

  ## --- compute observed stat ----------------------------------------------
  actual.wasserstein <- transport::wasserstein1d(x, y)

  ## --- permutation worker --------------------------------------------------
  combined.sample <- c(x, y)
  n.x <- length(x); n.y <- length(y)
  permute.and.compute <- function(i) {
    permuted.sample <- sample(combined.sample)
    permuted.x <- permuted.sample[1:n.x]
    permuted.y <- permuted.sample[(n.x + 1L):(n.x + n.y)]
    transport::wasserstein1d(permuted.x, permuted.y)
  }

  ## --- run permutations (parallel if allowed) -----------------------------
  n.cores <- choose_cores(n.cores)
  if (n.cores > 1L) {
    if (use_fork()) {
      null.wasserstein <- unlist(parallel::mclapply(seq_len(n.perms), permute.and.compute,
                                                   mc.cores = n.cores), use.names = FALSE)
    } else {
      cl <- parallel::makePSOCKcluster(n.cores)  # .check_ncores() will allow up to 2 during checks
      on.exit(parallel::stopCluster(cl), add = TRUE)
      null.wasserstein <- unlist(parallel::parLapply(cl, seq_len(n.perms), permute.and.compute),
                                 use.names = FALSE)
    }
  } else {
    null.wasserstein <- numeric(n.perms)
    for (i in seq_len(n.perms)) null.wasserstein[i] <- permute.and.compute(i)
  }

  ## --- p-value -------------------------------------------------------------
  p.value <- mean(null.wasserstein >= actual.wasserstein)

  list(
    wasserstein1d = actual.wasserstein,
    p.value = p.value,
    null.wasserstein1d = null.wasserstein,
    n.cores.used = n.cores
  )
}

#' A function to compute the Wasserstein distance between two samples from one-dimensional distributions using nearest neighbors strategy
#'
#' The samples are assumed to be of the same size.
#'
#' @param x    A sample from a one-dimensional distribution.
#' @param y    A sample from a one-dimensional distribution.
## NN.wasserstein.distance.1D <- function(x, y)
greedy.NN.transport.distance.1D <- function(x, y)
{
    n <- length(x)
    stopifnot(length(y)==n)
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(all(is.finite(x)))
    stopifnot(all(is.finite(y)))

    W <- 0
    i <- 1
    while ( length(x) > 1 )
    {
        nn <- get.knnx(cbind(y), cbind(x), k=1)
        d <- nn$nn.dist[,1]
        ##
        ## cat("----",i,"\n")
        ## cat("x:", x)
        ## cat("\ny:", y)
        ## cat("\nd:", d)
        ##
        nn.i <- nn$nn.index[,1]
        i.min <- which.min(d)
        j.min <- nn.i[i.min]
        ##print(d[i.min])
        W <- W + d[i.min]
        x <- x[-i.min]
        y <- y[-j.min]
    }

    W <- W + abs(x[1] - y[1])

    W / n
}

#' Local Wasserstein distance to the uniform and delta at 1 of the inner products of the NN unit vectors
#'
#' For each point x of a state space, S, the unit vectors starting at x and
#' ending at x's K-NN's are defined. The routine estimates the entropy of the
#' inner products <n_i, n_K> of all these unit vectors with the last one.
#'
#' @param S                 A state space.
#' @param K                 The number of nearest neighbors to use for the identification of the end points
#' @param use.transport     Set to TRUE, to use wasserstein1d() from transport package.
#' @param n.breaks          The number of bin for the estimate of distribution of cosine of the angle values.
#' @param ref.i             The index of the reference point whose cosine angles will be used as the reference for the Wd1 estimates.
#' @param use.geodesic.knn  Set to TRUE, if geodesic.knn() is to be used instead of get.knn().
#' @param nn.i              A matrix of indices of the nearest neighbors.
#' @param nn.d              A matrix of distances to the nearest neighbors.
#' @param verbose           Set to TRUE for messages indicating which part of the routine is currently run.
#'
#' @return A list with cosine of the angle values for each state of the state space.
##@param n.u               The number of uniformly distributed points over the interval from -1 to 1.
wasserstein.ipNNuv <- function(S,
                              K=25,
                              use.transport=TRUE,
                              n.breaks=50,
                              ref.i=2813,
                              use.geodesic.knn=TRUE,
                              nn.i=NULL,
                              nn.d=NULL,
                              verbose=TRUE)
{
    ## get K-NN's of each point
    if ( is.null(nn.i) || is.null(nn.d) )
    {
        if ( use.geodesic.knn )
        {
            if ( verbose ) {
                cat("Calculating geodesic kNN's ... ")
                ptm <- proc.time()
            }

            stop("geodesic.knn function is not available in this version. Please set use.geodesic.knn = FALSE")
            # nn <- geodesic.knn(S, k=K)
            nn.i <- nn$nn.index[,-1]
            nn.d <- nn$nn.dist[,-1]

            if ( verbose ) {
                elapsed.time(ptm)
            }
        } else {

            if ( verbose ) {
                cat("Calculating ambient space kNN's ... ")
                ptm <- proc.time()
            }

            nn <- get.knn(S, k=K)
            nn.i <- nn$nn.index
            nn.d <- nn$nn.dist

            if ( verbose ) {
                elapsed.time(ptm)
            }
        }
    }

    ## For each point compute the cosine of the angle between the last nearest
    ## neighbor and all other nearest neighbors

    if ( verbose ) {
        cat("Calculating local entropy of the inner products of NN's unit vectors over [-1,1] ... ")
        ptm <- proc.time()
    }

    ## determining the reference point coss
    if ( use.transport ) {
        U <- seq(-1,1, length.out = 1000)
        delta1 <- rep(1, 1000)
    } else {

        i <- ref.i
        x <- as.numeric(S[i,])
        x.nn <- S[nn.i[i,],]
        n <- nrow(x.nn)
        L.nn <- numeric(n)
        for ( j in seq(n) )
        {
            nj <- x.nn[j,] - x
            L.nn[j] <- sqrt(sum(nj^2))
        }
        idx <- L.nn>0
        x.nn <- x.nn[idx,]
        L.nn <- L.nn[idx]
        m <- length(L.nn)
        nK <- x.nn[m,] - x
        nK <- nK / sqrt(sum(nK^2))
        jj <- seq(m-1)
        coss <- numeric(m-1)
        for ( j in jj )
        {
            nj <- x.nn[j,] - x
            nj <- nj / sqrt(sum(nj^2))
            coss[j] <- sum( nK * nj ) # inner product b/w nK and nj = cosine of the angle between them
        }

        coss[!is.finite(coss)] <- 1

        n <- ncol(nn.i)-1
        if ( length(coss) < n ) {
            m <- length(coss) + 1
            coss[m:n] <- 1
        }

        delta1 <- coss
    }

    n.samples <- nrow(nn.i)
    Wd1 <- rep(NA, n.samples)
    Wu <- rep(NA, n.samples)
    H <- rep(NA, n.samples)
    cos.df <- matrix(NA, nrow=n.samples, ncol=n)
    for ( i in seq(n.samples) )
    {
        cat("\r",i)
        x <- as.numeric(S[i,])
        x.nn <- S[nn.i[i,],]
        ##rgl::spheres3d(x.nn, col='blue', radius=0.0002)
        n <- nrow(x.nn)
        L.nn <- numeric(n)
        for ( j in seq(n) )
        {
            nj <- x.nn[j,] - x
            L.nn[j] <- sqrt(sum(nj^2))
        }
        idx <- L.nn>0
        x.nn <- x.nn[idx,]
        L.nn <- L.nn[idx]
        m <- length(L.nn)
        nK <- x.nn[m,] - x
        nK <- nK / sqrt(sum(nK^2))
        jj <- seq(m-1)
        coss <- numeric(m-1)
        for ( j in jj )
        {
            nj <- x.nn[j,] - x
            nj <- nj / sqrt(sum(nj^2))
            coss[j] <- sum( nK * nj ) # inner product b/w nK and nj = cosine of the angle between them
        }
        cos.df[i,jj] <- coss
        ##
        d <- nn.d[i, jj]
        if ( length(unique(d)) > 10 && length(d[d>0]) > 10 ) {

            coss <- coss[is.finite(coss)]
            if ( use.transport ) {
                Wd1[i] <- transport::wasserstein1d(coss, delta1)
                Wu[i] <- transport::wasserstein1d(coss, U)
            } else {
                n.coss <- length(coss)
                d <- 0
                out <- .C(C_wasserstein_distance_1D,
                         as.double(coss),
                         as.double(delta1),
                         as.integer(n.coss),
                         d=as.double(d))
                Wd1[i] <- out$d
                ##
                u <- seq(-1,1, length.out = n.coss)
                out <- .C(C_wasserstein_distance_1D,
                         as.double(coss),
                         as.double(u),
                         as.integer(n.coss),
                         d=as.double(d))
                Wu[i] <- out$d
            }
            ## computing entropy
            h <- graphics::hist(coss, breaks=seq(-1,1, length.out = n.breaks), plot=FALSE)
            p <- h$density
            p <- p[p>0]
            p <- p / sum(p)
            H[i] <- -sum(p * log(p))

        }
    }

    if ( verbose ) {
        elapsed.time(ptm)
    }

    list(Wd1=Wd1,
         Wu=Wu,
         H=H,
         cos.df=cos.df,
         delta1=delta1,
         nn.i=nn.i,
         nn.d=nn.d)
}
