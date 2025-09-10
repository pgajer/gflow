##
## Library of grid forming routines
##

#' Creates a 2D grid of points around X.
#'
#' If X is NULL, x1.range and x2.range will be used to create a rectangular grid of n*n points.
#'
#' @param n    The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'               on each axis that are the basis of the grid.
#'
#' @param X    A set of points around which the grid is created if not NULL.
#'
#' @param f    A fraction of x1 and x2 range that the grid is extended to.
#'
#' @param eSDf An edge subdivision factor that is a scaling factor, such that edges of mstree(X) are
#'                subdivided if their length is greater than eSDf * mode(edge.len).
#'
#' @param gRf  A grid radius factor that is a scaling factor, such that points of the rectangular grid
#'                that are further than gRf*mode(edge.len) away from the closest
#'                subdivision of the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#'
#' @return     A data frame, grid, with n*n rows and 2 columns of the grid points in
#'             2D. If X is not NULL, grid contains only points nor more than eps away
#'             from the mstree(X).
#'
#' @importFrom FNN get.knnx
#' @export
create.2D.grid <- function(n, X, f=0.2, eSDf=1.5, gRf=4) {

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

    stopifnot(is.numeric(n))
    stopifnot(length(n)==1)
    stopifnot(ncol(X)==2)

    stopifnot(is.numeric(f))
    stopifnot(f>=0 && f<0.3)

    stopifnot(is.numeric(eSDf))
    stopifnot(is.numeric(gRf))

    ##
    ## Building a rectangular grid
    ##
    x1.range <- range(X[,1])
    x2.range <- range(X[,2])

    dx1 <- diff(x1.range)*f/2
    dx2 <- diff(x2.range)*f/2

    x1g <- seq(from=x1.range[1]-dx1, to=x1.range[2]+dx1, length.out=n)
    x2g <- seq(from=x2.range[1]-dx2, to=x2.range[2]+dx2, length.out=n)

    X.grid <- matrix(nrow=n*n, ncol=2)
    colnames(X.grid) <- c("x1","x2")
    X.grid[,1] <- rep(x1g, n)

    z <- c()
    for ( i in seq(x2g) )
        z <- c(z, rep(x2g[i], n))

    X.grid[,2] <- z

    ##
    ## Removing elments of X.grid far away from X
    ##
    r <- mstree(X)
    edges <- r$edges
    edge.lens <- r$edge.lens

    mode.edge.len <- mode.1D(edge.lens)

    e.pts <- c()
    edge.2eps <- eSDf * mode.edge.len
    for ( j in seq(edge.lens) )
    {
        e.d <- edge.lens[j]
        if ( e.d > edge.2eps )
        {
            e.i <- edges[j,]
            e.start <- X[e.i[1],]
            e.end <- X[e.i[2],]
            ##
            t <- seq(0, e.d, by=mode.edge.len)
            t <- t[-1]
            for ( i in seq(t) )
            {
                p <- t[i] / e.d
                e.pts <- rbind(e.pts, p*e.start + (1-p)*e.end)
            }
        }
    }

    mst.grid <- rbind(X, e.pts)
    ## nrow(mst.grid)

    nn <- get.knnx(mst.grid, X.grid, k=1)
    nn.d <- nn$nn.dist

    idx <- nn.d[,1] < gRf * mode.edge.len
    ##sum(idx)
    X.grid <- X.grid[idx,]

    list(X.grid=X.grid,
         mst.grid=mst.grid,
         T.edges=edges,
         dx1=dx1,
         dx2=dx2)
}

#' Creates a equi-number (within each axis) 2D grid using a C routine
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#'
create.ENPs.grid.2D <- function(n, x1.range, x2.range, f=0.2)
{
    X.grid <- numeric(2*n*n)

    out <- .C("C_create_ENPs_grid_2D",
             as.integer(n),
             as.double(x1.range[1]),
             as.double(x1.range[2]),
             as.double(x2.range[1]),
             as.double(x2.range[2]),
             as.double(f),
             grid=as.double(X.grid))


    matrix(out$grid, nrow=n^2, ncol=2, byrow=TRUE)
}

#' Creates a equi-numer (within each axis) 3D grid using a C routine
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param x3.range   A range of x3 values - the third coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#'
create.ENPs.grid.3D <- function(n, x1.range, x2.range, x3.range, f=0.2)
{
    X.grid <- numeric(3*n^3)

    out <- .C("C_create_ENPs_grid_3D",
             as.integer(n),
             as.double(x1.range[1]),
             as.double(x1.range[2]),
             as.double(x2.range[1]),
             as.double(x2.range[2]),
             as.double(x3.range[1]),
             as.double(x3.range[2]),
             as.double(f),
             grid=as.double(X.grid))


    matrix(out$grid, nrow=n^3, ncol=3, byrow=TRUE)
}


#' Creates a equi-distant (all edges are of equal length) 2D grid using a C routine
#'
#' @param dx         The length of each edge.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#' @param verbose    If TRUE, some debugging messages will be printed.
create.ED.grid.2D <- function(dx, x1.range, x2.range, f=0.2, verbose=FALSE)
{
    Dx1 <- diff(x1.range)
    Dx2 <- diff(x2.range)

    x1L <- x1.range[1]
    x1R <- x1.range[2]
    x2L <- x2.range[1]
    x2R <- x2.range[2]

    if ( f > 0 )
    {
        f <- f/2

        dx1 <- Dx1 * f
        dx2 <- Dx2 * f

        x1L <- x1L - dx1
        x1R <- x1R + dx1

        x2L <- x2L - dx2
        x2R <- x2R + dx2

        Dx1 = x1R - x1L
        Dx2 = x2R - x2L
    }

    n1 <- ceiling( Dx1 / dx )
    n2 <- ceiling( Dx2 / dx )

    if ( verbose ) {
        print(sprintf("n1=%d\tn2=%d\n", n1, n2))
    }

    x1R <- x1L + n1 * dx
    x2R <- x2L + n2 * dx

    Dx1 <- x1R - x1L
    Dx2 <- x2R - x2L

    if ( verbose ) {
        print(sprintf("Updated range widths: Dx1=%.3f\tDx2=%.3f\n", Dx1, Dx2))
    }

    n1 <- n1 + 1
    n2 <- n2 + 1

    if ( verbose ) {
        print(sprintf("Updated ns: n1=%d\tn2=%d\n", n1, n2))
    }

    X.grid <- numeric(2*n1*n2)

    out <- .C("C_create_ED_grid_2D",
             as.double(dx),
             as.double(x1L),
             as.integer(n1),
             as.double(x2L),
             as.integer(n2),
             grid=as.double(X.grid))

    matrix(out$grid, nrow=n1*n2, ncol=2, byrow=TRUE)
}

#' Creates a equi-distant (all edges are of equal length) 3D grid using a C routine
#'
#' @param dx         The length of the grid edge.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param x3.range   A range of x3 values - the third coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#' @param verbose    If TRUE, some debugging messages will be printed.
create.ED.grid.3D <- function(dx, x1.range, x2.range, x3.range, f=0.1, verbose=FALSE)
{
    Dx1 <- diff(x1.range)
    Dx2 <- diff(x2.range)
    Dx3 <- diff(x3.range)

    x1L <- x1.range[1]
    x1R <- x1.range[2]

    x2L <- x2.range[1]
    x2R <- x2.range[2]

    x3L <- x3.range[1]
    x3R <- x3.range[2]

    if ( f > 0 )
    {
        f <- f/2

        dx1 <- Dx1 * f
        dx2 <- Dx2 * f
        dx3 <- Dx3 * f

        x1L <- x1L - dx1
        x1R <- x1R + dx1

        x2L <- x2L - dx2
        x2R <- x2R + dx2

        x3L <- x3L - dx3
        x3R <- x3R + dx3

        Dx1 = x1R - x1L
        Dx2 = x2R - x2L
        Dx3 = x3R - x3L
    }

    n1 <- ceiling( Dx1 / dx )
    n2 <- ceiling( Dx2 / dx )
    n3 <- ceiling( Dx3 / dx )

    if ( verbose )
        print(sprintf("n1=%d\tn2=%d\tn3=%d\n", n1, n2, n3))

    x1R <- x1L + n1 * dx
    x2R <- x2L + n2 * dx
    x3R <- x3L + n3 * dx

    Dx1 <- x1R - x1L
    Dx2 <- x2R - x2L
    Dx3 <- x3R - x3L

    if ( verbose )
        print(sprintf("Updated range widths: Dx1=%.3f\tDx2=%.3f\n", Dx1, Dx2))

    n1 <- n1 + 1
    n2 <- n2 + 1
    n3 <- n3 + 1

    if ( verbose )
        print(sprintf("Updated ns: n1=%d\tn2=%d\n", n1, n2))

    X.grid <- numeric(3*n1*n2*n3)

    out <- .C("C_create_ED_grid_3D",
             as.double(dx),
             as.double(x1L),
             as.integer(n1),
             as.double(x2L),
             as.integer(n2),
             as.double(x3L),
             as.integer(n3),
             grid=as.double(X.grid))

    matrix(out$grid, nrow=n1*n2*n3, ncol=3, byrow=TRUE)
}


#' Creates a rectangular 2D grid
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#' @param type       A type of distribution the points are sampled from.
#'
create.2D.rect.grid <- function(n, x1.range, x2.range, type=c("unif","runif","norm"), f=0.2)
{
    types <- c("unif","runif","norm")
    type <- match.arg(type, types)

    stopifnot(is.numeric(n))
    stopifnot(length(n)==1)

    stopifnot(is.numeric(f))
    stopifnot(f>=0 && f<0.3)

    ##
    ## Building a rectangular grid
    ##
    dx1 <- diff(x1.range)*f/2
    dx2 <- diff(x2.range)*f/2

    X.grid <- matrix(nrow=n*n, ncol=2)
    colnames(X.grid) <- c("x1","x2")

    if ( type == "unif" ){
        x1g <- seq(from=x1.range[1]-dx1, to=x1.range[2]+dx1, length.out=n)
        x2g <- seq(from=x2.range[1]-dx2, to=x2.range[2]+dx2, length.out=n)

        X.grid[,1] <- rep(x1g, n)
        z <- c()
        for ( i in seq(x2g) )
            z <- c(z, rep(x2g[i], n))
        X.grid[,2] <- z

    } else if ( type == "runif" ){
        ## x1g <- runif(n, min=x1.range[1]-dx1, max=x1.range[2]+dx1)
        ## x2g <- runif(n, min=x2.range[1]-dx2, max=x2.range[2]+dx2)
        for ( i in seq(nrow(X.grid)) )
        {
            X.grid[i,] <- c(runif(1, min=x1.range[1]-dx1, max=x1.range[2]+dx1),
                           runif(1, min=x2.range[1]-dx2, max=x2.range[2]+dx2))
        }

    } else if ( type == "rnorm" ){
        ## x1g <- rnorm(n, mean=mean(x1.range), sd=2)
        ## x2g <- rnorm(n, mean=mean(x2.range), sd=2)
        for ( i in seq(nrow(X.grid)) )
        {
            X.grid[i,] <- c(rnorm(1, mean=mean(x1.range), sd=2),
                           rnorm(1, mean=mean(x2.range), sd=2))
        }
    }

    X.grid
}


#' Creates a 3D grid of points around X
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param X          A set of points around which the grid is created if not NULL.
#' @param f          A fraction of x1, x2 and x3 range that the grid is extended to.
#'
#' @param min.gSf A scaling factor, indicating the minimal size of
#'     the grid around X, with respect to the size of X. That is the minimal
#'     size of the grid needs to be at least min.gSf*|X|, where |X|
#'     is the size of X. For example, with min.gSf set to 1.5, the
#'     minimal size of the grid will be, if possible, at least 1.5*|X|. For this
#'     to be possible the initial rectangular grid size has to be at least of
#'     that minimal size, but in practice much bigger. It overrides
#'     gRf if the size of the grid given by the value of
#'     gRf is not
#'
#' @param eSDf  An edge subdivision factor that is a scaling factor, such that edges of mstree(X) are
#'                     subdivided is their length is greater than eSDf * mode(edge.len).
#'
#' @param gRf A scaling factor, such that points of the rectangular grid
#'                    that are further than gRf*mode.edge.len away from the closest
#'                    subdivision of the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#'
#' @return           A data frame, grid, with n*n*n rows and 3 columns of the grid points in
#'                    3D. If X is not NULL, grid contains only points nor more than eps away
#'                    from the mstree(X).
#'
create.3D.grid <- function(n, X, f=0.2, min.gSf=1.5, eSDf=1.5, gRf=2) {

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

    stopifnot(is.numeric(n))
    stopifnot(length(n)==1)
    stopifnot(ncol(X)==3)

    stopifnot(is.numeric(f))
    stopifnot(f>=0 && f<0.3)

    stopifnot(is.numeric(eSDf))
    stopifnot(is.numeric(gRf))

    x1.range <- range(X[,1])
    x2.range <- range(X[,2])
    x3.range <- range(X[,3])

    ## dx1 <- diff(x1.range)*f/2
    ## dx2 <- diff(x2.range)*f/2
    dx3 <- diff(x3.range)*f/2

    ## x1g <- seq(from=x1.range[1]-dx1, to=x1.range[2]+dx1, length.out=n)
    ## x2g <- seq(from=x2.range[1]-dx2, to=x2.range[2]+dx2, length.out=n)
    x3g <- seq(from=x3.range[1]-dx3, to=x3.range[2]+dx3, length.out=n)

    X.grid <- matrix(nrow=n^3, ncol=3)
    X.grid <- as.data.frame(X.grid)
    colnames(X.grid) <- c("x1","x2","x3")

    g2 <- create.2D.rect.grid(n, x1.range=x1.range, x2.range=x2.range, f=f)

    X.grid <- c()
    for ( i in seq(x3g) )
        X.grid <- rbind(X.grid, cbind(g2,x3g[i]))

    r <- mstree(X)
    edges <- r$edges
    edge.lens <- r$edge.lens

    mode.edge.len <- mode.1D(edge.lens)

    e.pts <- c()
    edge.2eps <- eSDf * mode.edge.len
    for ( j in seq(edge.lens) )
    {
        e.d <- edge.lens[j]
        if ( e.d > edge.2eps )
        {
            e.i <- edges[j,]
            e.start <- X[e.i[1],]
            e.end <- X[e.i[2],]
            ##
            t <- seq(0, e.d, by=mode.edge.len)
            t <- t[-1]
            for ( i in seq(t) )
            {
                p <- t[i] / e.d
                e.pts <- rbind(e.pts, p*e.start + (1-p)*e.end)
            }
        }
    }

    mst.grid <- rbind(X, e.pts)
    ##nrow(mst.grid)

    nn <- get.knnx(mst.grid, X.grid, k=1)
    nn.d <- nn$nn.dist

    min.grid.size <- min.gSf * nrow(X)

    idx <- nn.d[,1] < gRf * mode.edge.len

    if ( sum(idx) >= min.grid.size )
    {
        X.grid <- X.grid[idx,]

    } else {

        max.d <- max(nn.d[,1])
        max.gRf <- max.d / mode.edge.len # this value of
                                                       # gRf
                                                       # accepts all elements of
                                                       # X.grid in the final grid.
        n <- 1000
        GRFs <- seq(gRf, max.gRf, length=n)

        i <- 2
        grf <- GRFs[i]
        idx <- (nn.d[,1] < grf * mode.edge.len)

        while ( i < n && sum(idx) < min.grid.size ) {
            i <- i + 1
            grf <- GRFs[i]
            idx <- (nn.d[,1] < grf * mode.edge.len)
        }

        X.grid <- X.grid[idx,]

        gRf <- grf
    }

    list(X.grid=X.grid,
         mst.grid=mst.grid,
         T.edges=edges,
         gRf=gRf)
}


#' Creates a tubular neighborhood 3D grid of X.
#'
#' This is a complete rewrite of the predecessor routine create.3D.grid.v2()
#' with gSf replaced by the length of the grid edge (all edges in these grids
#' have the same length).
#'
#' @param X         A set of points around which the grid is created if not NULL.
#' @param mst.grid     A matrix of 3D points that are subdivisions of edges of the mstree(X).
#' @param mode.edge.len The mode of the edge lengths.#'
#'
#' @param dx       The length of the grid edge.
#' @param gRf      A scaling factor, such that points of the rectangular grid
#'                  that are further than gRf*mode.edge.len away from the closest
#'                  subdivision of the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#'
#' @param f        A fraction of x1, x2 and x3 range that the grid is extended to.
#'
#' @param verbose  Set to TRUE to see info on what is being done.
#'
#' @return         A data frame, grid, with n*n*n rows and 3 columns of the grid points in
#'                    3D. If X is not NULL, grid contains only points nor more than eps away
#'                    from the mstree(X).
#'
#' @importFrom FNN get.knnx
#' @export
create.3D.TN.grid <- function(X, mst.grid, dx, mode.edge.len, gRf, f=0.05, verbose=FALSE)
{
    grid <- create.ED.grid.3D(dx, x1.range=range(X[,1]), x2.range=range(X[,2]), x3.range=range(X[,3]), f)

    nn <- get.knnx(mst.grid, grid, k=1)
    d <- nn$nn.dist[,1]

    gRf.thld <- gRf * mode.edge.len
    idx <- d < gRf.thld

    X.grid <- grid[idx,]

    X.grid
}

#' Creates a tubular neighborhood 3D grid of X returning a list with a grid
#' around X and the distance to the boundary of the grid.
#'
#' @param X         A set of points around which the grid is created if not NULL.
#' @param mst.grid     A matrix of 3D points that are subdivisions of edges of the mstree(X).
#' @param mode.edge.len The mode of the edge lengths.#'
#'
#' @param dx       The length of the grid edge.
#' @param gRf      A scaling factor, such that points of the rectangular grid
#'                  that are further than gRf*mode.edge.len away from the closest
#'                  subdivision of the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#'
#' @param f        A fraction of x1, x2 and x3 range that the grid is extended to.
#'
#' @param verbose  Set to TRUE to see info on what is being done.
#'
#' @return         A data frame, grid, with n*n*n rows and 3 columns of the grid points in
#'                    3D. If X is not NULL, grid contains only points nor more than eps away
#'                    from the mstree(X).
#'
#' @importFrom FNN get.knnx
#' @export
create.3D.TN.grid.plus <- function(X, mst.grid, dx, mode.edge.len, gRf, f=0.05, verbose=FALSE)
{
    grid <- create.ED.grid.3D(dx, x1.range=range(X[,1]), x2.range=range(X[,2]), x3.range=range(X[,3]), f)

    nn <- get.knnx(mst.grid, grid, k=1)
    d <- nn$nn.dist[,1]

    gRf.thld <- gRf * mode.edge.len
    idx <- d < gRf.thld

    X.grid <- grid[idx,]
    d.grid <- 1 - d[idx]/gRf.thld

    list(X.grid=X.grid,
         d.grid=d.grid)
}


#' Creates a grid around a state space.
#'
#' This function generates a uniform grid in a tubular neighborhood of a state space.
#'
#' @param X            A state space for which the grid is to be created.
#' @param gSf          A grid size scaling factor.
#' @param gRf          A grid radius scaling factor. Points on the rectangular grid that
#'                     are farther than gRf*mode.edge.len from the closest subdivision of
#'                     the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#' @param min.K        The minimum number of nearest neighbors ('x' NNs) that must be present in each window.
#' @param med.dK.divf  A division factor for the lower bound of bandwidth (bw) and dx. It is equal to
#'                     med.dK/med.dK.divf, where \code{med.dK=median(dK)} and \code{dK=nn.dist[,min.K]}.
#' @param max.dx.C     A division factor for estimating max.dx, defined as the length of the diagonal in
#'                     the X enclosing box divided by max.dx.C.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{X.grid}: A uniform grid in a tubular neighborhood of X.
#' \item \code{d.grid}: The distance to the boundary of the grid.
#' \item \code{mst.grid}: A uniform grid over the edges of the minimum spanning tree of X.
#' \item \code{opt.dx}: The optimal value of the dx parameter for given gSf and gRf values.
#' \item \code{mode.edge.len}: The mode of the lengths of the edges of the minimum spanning tree of X.
#' }
#'
#' @examples
#' \dontrun{
#' # Let X be a low-dimensional model of a state space.
#' res <- create.X.grid(X, gSf=5, gRf=5, min.K=10, med.dK.divf=5, max.dx.C=1)
#' str(res)
#' }
create.X.grid <- function(X, gSf, gRf, min.K=10, med.dK.divf=5, max.dx.C=1)
{
    nrX <- nrow(X)

    rr <- mstree.grid(X)
    mst.grid <- rr$mst.grid
    mode.edge.len <- rr$mode.edge.len

    log.dx.fn <- function(logdx, gSf=2, gRf=1)
    {
        X.grid <- create.3D.TN.grid(X, mst.grid, exp(logdx), mode.edge.len, gRf)
        if ( is.null(dim(X.grid)) ) {
            X.grid <- rbind(X.grid)
        }
        d <- abs(nrow(X.grid)/nrX - gSf)
        ##cat(sprintf("dx=%0.3g gRf=%.2g X.grid=%d deviation=%.3g\n", dx, gRf, nrX.grid, d))
        d
    }

    nn <- get.knn(X, k=min.K)
    dK <- nn$nn.dist[,min.K]
    med.dK <- median(dK)
    ## dx range
    min.dx <- med.dK / med.dK.divf

    r <- numeric(3)
    for ( i in 1:3 ) {
        r[i] <- diff(range(X[,i]))
    }
    max.dx <- sqrt(sum(r^2))
    max.dx <- max.dx / max.dx.C

    log.min.dx <- log(min.dx)
    log.max.dx <- log(max.dx)

    ##
    ## Checking that nrow(X.grid|min.dx) > nrX and nrow(X.grid|max.dx) < nrX
    ##

    ## min.dx
    min.X.grid <- create.3D.TN.grid(X, mst.grid, min.dx, mode.edge.len, gRf=3)
    if ( is.null(ncol(min.X.grid)) ) {
        min.X.grid <- rbind(min.X.grid)
    }
    min.dev <- abs(nrow(min.X.grid)/nrX - 1)

    ## max.dx
    max.X.grid <- create.3D.TN.grid(X, mst.grid, max.dx, mode.edge.len, gRf=3)
    if ( is.null(ncol(max.X.grid)) ) {
        max.X.grid <- rbind(max.X.grid)
    }
    max.dev <- abs(nrow(max.X.grid)/nrX - 1)

    dev.thld <- 0.005
    opt.dx <- NULL

    if ( min.dev < dev.thld )
    {
        X.grid <- min.X.grid
        opt.dx <- min.dx

    } else if ( max.dev < dev.thld ) {

        X.grid <- max.X.grid
        opt.dx <- max.dx

    } else {

        ## In case min.dx is not small enough
        while ( nrow(min.X.grid) < nrX )  {
            min.dx <- 0.1 * min.dx
            ##min.X.grid <- create.3D.TN.grid(X, mst.grid, min.dx, mode.edge.len, gRf=3)
            grid.obj <- create.3D.TN.grid.plus(X, mst.grid, min.dx, mode.edge.len, gRf=3)
            min.X.grid <- grid.obj$X.grid
            d.grid <- grid.obj$d.grid
        }

        ## In case max.dx is not large enough
        while ( nrow(max.X.grid) > nrX )  {
            max.dx <- 2 * max.dx
            ##max.X.grid <- create.3D.TN.grid(X, mst.grid, max.dx, mode.edge.len, gRf=3)
            grid.obj <- create.3D.TN.grid.plus(X, mst.grid, max.dx, mode.edge.len, gRf=3)
            max.X.grid <- grid.obj$X.grid
            d.grid <- grid.obj$d.grid
        }

        ## Finging dx for which gSf=1 and gRf=3
        r <- optimize(log.dx.fn, interval=c(log.min.dx, log.max.dx), gSf=gSf, gRf=gRf)
        opt.dx <- exp(r$minimum)

        ##X.grid <- create.3D.TN.grid(X, mst.grid, opt.dx, mode.edge.len, gRf)
        grid.obj <- create.3D.TN.grid.plus(X, mst.grid, opt.dx, mode.edge.len, gRf)
        X.grid <- grid.obj$X.grid
        d.grid <- grid.obj$d.grid
    }

    grid.obj <- list()
    grid.obj$X.grid <- X.grid
    grid.obj$d.grid <- d.grid
    grid.obj$mst.grid <- mst.grid
    grid.obj$opt.dx <- opt.dx
    grid.obj$mode.edge.len <- mode.edge.len

    grid.obj
}


#' Creates a 3D grid of points around X.
#'
#' In this version if the min.gSf condition is not satisfied, n is
#' increases until it satisfies that condition. Thus, gRf is never modified.
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param X          A set of points around which the grid is created if not NULL.
#' @param f          A fraction of x1, x2 and x3 range that the grid is extended to.
#'
#' @param min.gSf A scaling factor, indicating the minimal size of
#'     the grid around X, with respect to the size of X. That is the minimal
#'     size of the grid needs to be at least min.gSf*|X|, where |X|
#'     is the size of X. For example, with min.gSf set to 1.5, the
#'     minimal size of the grid will be, if possible, at least 1.5*|X|. For this
#'     to be possible the initial rectangular grid size has to be at least of
#'     that minimal size, but in practice much bigger. It overrides
#'     gRf if the size of the grid given by the value of
#'     gRf is not
#'
#' @param min.min.gSf The minimal allowable value of min.gSf.
#'
#' @param eSDf A scaling factor, such that edges of mstree(X) are
#'                     subdivided is their length is greater than eSDf * mode(edge.len).
#'
#' @param gRf A scaling factor, such that points of the rectangular grid
#'                    that are further than gRf*mode.edge.len away from the closest
#'                    subdivision of the minimal spanning tree of X, mstree(X), are eliminated from the grid.
#'
#' @param min.gRf The minimal allowable value of gRf.
#'
#' @param max.n       The maximal value of n that can be reached as n is adjusted to
#'                    satisfy the size and radius factor conditions.
#'
#' @param mst.grid      A matrix of 3D points that are subdivisions of edges of the mstree(X).
#' @param mode.edge.len The mode of the edge lengths.
#' @param verbose    Set to TRUE to see info on what is being done.
#'
#' @return           A data frame, grid, with n*n*n rows and 3 columns of the grid points in
#'                    3D. If X is not NULL, grid contains only points nor more than eps away
#'                    from the mstree(X).
#'
create.3D.grid.v2 <- function(n, X, f=0.2, mst.grid=NULL, mode.edge.len=NULL,
                             min.gSf=1.5,
                             min.min.gSf=1,
                             gRf=2,
                             min.gRf=1,
                             eSDf=1.5,
                             max.n=120, verbose=FALSE) {
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

    stopifnot(is.finite(n))
    ## stopifnot(is.numeric(n))
    ## stopifnot(length(n)==1)

    ## stopifnot(!is.null(X))
    ## stopifnot(is.numeric(X))
    ## stopifnot(is.matrix(X))
    ## stopifnot(ncol(X)==3)

    ## stopifnot(is.numeric(f))
    ## stopifnot(f>=0 && f<0.3)

    ## stopifnot(is.numeric(eSDf))
    ## stopifnot(is.numeric(gRf))

    if ( gRf < min.gRf ){
        warning("gRf < min.gRf: gRf set to 0.1")
        gRf <- min.gRf
    }

    if ( min.gSf < min.min.gSf ){
        warning("min.gSf < min.min.gSf; min.gSf set to 1.")
        min.gSf <- min.min.gSf
    }

    ## browser()

    min.grid.size <- min.gSf * nrow(X)

    rr <- get.3D.rect.grid.NN(n, X, f, eSDf, mst.grid, mode.edge.len)
    idx <- rr$d < gRf * rr$mode.edge.len
    proposed.grid.size <- sum(idx)

    if ( proposed.grid.size < min.grid.size )
    {
        if ( verbose )
            print(sprintf("n=%d  sum(idx)=%d  min.grid.size=%d; sum(idx) <  min.grid.size; min.gSf=%f  nrow(X)=%d",
                          n, proposed.grid.size, as.integer(min.grid.size), min.gSf, nrow(X)))

        if ( proposed.grid.size > 0 ){
            good.grid.size.factor <- as.integer((min.grid.size / proposed.grid.size)^(1/3))
            n <- n * good.grid.size.factor
        } else {
            ##n <- n + 1
            ##n <- n * as.integer( min.grid.size^(1/3) )
            n <- as.integer( round( nrow(X) / 10 ) )
        }

        rr <- get.3D.rect.grid.NN(n, X, f, eSDf, rr$mst.grid, rr$mode.edge.len)
        idx <- rr$d < gRf * rr$mode.edge.len
        proposed.grid.size <- sum(idx)

        if ( verbose )
            print(sprintf("n=%d  sum(idx)=%d", n, proposed.grid.size))

        if ( proposed.grid.size < min.grid.size )
        {
            if ( verbose )
                print("In proposed.grid.size < min.grid.size")

            while ( n < max.n && (proposed.grid.size < min.grid.size)) {
                n <- n + 1
                rr <- get.3D.rect.grid.NN(n, X, f, eSDf, rr$mst.grid, rr$mode.edge.len)
                idx <- rr$d < gRf * rr$mode.edge.len
                proposed.grid.size <- sum(idx)

                if ( verbose )
                    print(sprintf("n=%d  sum(idx)=%d", n, proposed.grid.size))
            }

        } else { # finding the smallest n that satisfies the condition proposed.grid.size >= min.grid.size

            if ( verbose )
                print("In proposed.grid.size >= min.grid.size")

            while ( n > 5 && (proposed.grid.size > min.grid.size)) {
                n <- n - 1
                rr <- get.3D.rect.grid.NN(n, X, f, eSDf, rr$mst.grid, rr$mode.edge.len)
                idx <- rr$d < gRf * rr$mode.edge.len
                proposed.grid.size <- sum(idx)

                if ( verbose )
                    print(sprintf("n=%d  sum(idx)=%d", n, proposed.grid.size))
            }

            n <- n + 1
            rr <- get.3D.rect.grid.NN(n, X, f, eSDf, rr$mst.grid, rr$mode.edge.len)
            idx <- rr$d < gRf * rr$mode.edge.len

            if ( verbose )
                print(sprintf("n=%d  sum(idx)=%d", n, sum(idx)))
        }
    }

    X.grid <- rr$X.grid[idx,]

    list(X.grid=X.grid,
         mst.grid=rr$mst.grid,
         mode.edge.len=rr$mode.edge.len,
         ##T.edges=rr$edges,
         n=n,
         gRf=gRf)
}


#' Creates a rectangular 3D grid and returns a vector of NN's
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param X          A set of points around which the grid is created if not NULL.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#'
#' @param eSDf A scaling factor, such that edges of mstree(X) are
#'                     subdivided is their length is greater than eSDf * mode(edge.len).
#' @param mst.grid      A matrix of 3D points that are subdivisions of edges of the mstree(X).
#' @param mode.edge.len The mode of the edge lengths.
#'
get.3D.rect.grid.NN <- function(n, X, f=0.2, eSDf=1.5, mst.grid=NULL, mode.edge.len=NULL)
{
    stopifnot(is.finite(n))

    x1.range <- range(X[,1])
    x2.range <- range(X[,2])
    x3.range <- range(X[,3])

    ## dx1 <- diff(x1.range)*f/2
    ## dx2 <- diff(x2.range)*f/2
    dx3 <- diff(x3.range)*f/2

    ##browser()

    ## x1g <- seq(from=x1.range[1]-dx1, to=x1.range[2]+dx1, length.out=n)
    ## x2g <- seq(from=x2.range[1]-dx2, to=x2.range[2]+dx2, length.out=n)
    x3g <- seq(from=x3.range[1]-dx3, to=x3.range[2]+dx3, length.out=n)

    X.grid <- matrix(nrow=n^3, ncol=3)
    X.grid <- as.data.frame(X.grid)
    colnames(X.grid) <- c("x1","x2","x3")

    g2 <- create.2D.rect.grid(n, x1.range=x1.range, x2.range=x2.range, f=f)

    ## X.grid <- c()
    ## for ( i in seq(x3g) )
    ##     X.grid <- rbind(X.grid, cbind(g2,x3g[i]))
    X.grid <- matrix(nrow=(nrow(g2)*length(x3g)), ncol=3)
    nr.g2 <- nrow(g2)
    for ( i in seq(x3g) )
    {
        ## X.grid <- rbind(X.grid, cbind(g2,x3g[i]))
        ii <- ((i - 1)*nr.g2 + 1):(i*nr.g2)
        X.grid[ii,] <- cbind(g2,x3g[i])
    }

    if ( is.null(mst.grid) || is.null(mode.edge.len) )
    {
        rr <- mstree.grid(X, eSDf=eSDf)
        mst.grid <- rr$mst.grid
        mode.edge.len <- rr$mode.edge.len
    }

    ## cat("Running get.knnx(mst.grid, X.grid, k=1) ... ")
    ## ptm <- proc.time()
    nn <- get.knnx(mst.grid, X.grid, k=1)
    ## elapsed.time(ptm)

    list(d=nn$nn.dist[, 1],
         X.grid=X.grid,
         mode.edge.len=mode.edge.len,
         mst.grid=mst.grid)
    ## edges=edges)
}

#' Creates a quasi-uniform grid over the edges of the minimal spanning tree of a state space
#'
#' For each edge of the the minimal spanning tree of X whose length is greater
#' than eSDf * m a uniform grid with the distance between consecutive points
#' equal to m, where m = mode(edge lengths).
#'
#' @param X A set of points for which minimal spanning tree edge subdivision is
#'     to be created.
#' @param eSDf An edge subdivision factor, such that a uniform grid is created over the
#'     given edge of mstree(X) if its length is greater than eSDf *
#'     mode(edge.len).
mstree.grid <- function(X, eSDf = 1.5) {

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

    r <- mstree(X)
    edges <- r$edges
    edge.lens <- r$edge.lens

    ## the mode of the edge lengths
    n.edge.thld <- 50
    if ( length(edge.lens) > n.edge.thld ) {
        mode.edge.len <- mode.1D(edge.lens)
    } else {
        mode.edge.len <- median(edge.lens)
    }

    if ( is.na(mode.edge.len) ) {
        stop("mode.edge.len is NULL")
    }

    grid.size <- floor(sum(edge.lens) / mode.edge.len)
    grid <- matrix(nrow=grid.size, ncol = ncol(X))
    edge.len.thld <- eSDf * mode.edge.len

    grid.index <- 1
    for ( j in seq(edge.lens) ) {

        e.d <- edge.lens[j]

        if ( e.d > edge.len.thld ) {

            e.i <- edges[j,]
            e.start <- X[e.i[1],]
            e.end <- X[e.i[2],]

            t <- seq(0, e.d, by = mode.edge.len)
            t <- t[-1]

            for ( i in seq(t) ) {
                p <- t[i] / e.d
                grid[grid.index,] <- p*e.start + (1-p)*e.end
                grid.index <- grid.index + 1
            }
        }
    }

    ## removing NA rows/pts
    idx <- !is.na(grid[,1])
    grid <- grid[idx,]

    mst.grid <- rbind(X, grid)

    list(mode.edge.len = mode.edge.len,
         mst.grid = mst.grid,
         grid = grid,
         edges = edges,
         edge.lens = edge.lens,
         edge.len.thld = edge.len.thld)
}

#' Creates a rectangular 3D grid
#'
#' @param n          The number of uniformly spaced points, seq(min(xi), max(xi), length=n),
#'                    on each axis that are the basis of the grid.
#' @param x1.range   A range of x1 values - the first coordinate.
#' @param x2.range   A range of x2 values - the second coordinate.
#' @param x3.range   A range of x2 values - the second coordinate.
#' @param f          A fraction of x1 and x2 range that the grid is extended to.
#' @param type       A type of distribution the points are sampled from.
#'
create.3D.rect.grid <- function(n, x1.range, x2.range, x3.range, type=c("unif","runif","norm"), f=0.2)
{
    types <- c("unif","runif","norm")
    type <- match.arg(type, types)

    stopifnot(is.numeric(n))
    stopifnot(length(n)==1)

    stopifnot(is.numeric(f))
    stopifnot(f>=0 && f<0.3)

    ##
    ## Building a rectangular grid
    ##
    dx1 <- diff(x1.range)*f/2
    dx2 <- diff(x2.range)*f/2
    dx3 <- diff(x3.range)*f/2

    X.grid <- matrix(nrow=n^3, ncol=3)
    X.grid <- as.data.frame(X.grid)
    colnames(X.grid) <- c("x1", "x2", "x3")

    if ( type == "unif" ){

        x3g <- seq(from=x3.range[1]-dx3, to=x3.range[2]+dx3, length.out=n)
        g2 <- create.2D.rect.grid(n, x1.range=x1.range, x2.range=x2.range, f=f)

        X.grid <- c()
        for ( i in seq(x3g) )
            X.grid <- rbind(X.grid, cbind(g2,x3g[i]))

    } else if ( type == "runif" ){

        for ( i in seq(nrow(X.grid)) )
        {
            X.grid[i,] <- c(runif(1, min=x1.range[1]-dx1, max=x1.range[2]+dx1),
                           runif(1, min=x2.range[1]-dx2, max=x2.range[2]+dx2),
                           runif(1, min=x3.range[1]-dx3, max=x3.range[2]+dx3))
        }

    } else if ( type == "rnorm" ){

        for ( i in seq(nrow(X.grid)) )
        {
            X.grid[i,] <- c(rnorm(1, mean=mean(x1.range), sd=2),
                           rnorm(1, mean=mean(x2.range), sd=2),
                           rnorm(1, mean=mean(x3.range), sd=2))
        }
    }

    X.grid
}

#' Given gSf and gRf, dx are found so that corresponding grid(X) has give gSf.
#' The resulting data is used to build models allowing to construct grid(X) with
#' the prespecified values of gRf and gSf.
#'
#' @param X           A state space matrix.
#' @param min.K       The number of nearest neighbors used to estimate median min.K.
#' @param med.dK.divf A med.dK division factor for the lower bound of bw and dx to be med.dK/med.dK.divf, where \code{med.dK=median(dK)} with \code{dK=nn.dist[,min.K]}.
#' @param max.dx.C    A division factor for max.dx estimate such that max.dx = (length of the diagonal in the X encolosing box) / max.dx.C.
#' @param max.dx      The maximal value of dx.
#' @param min.gRf     The minimum of gRf.
#' @param max.gRf     The maximum of gRf.
#' @param n.gRfs      The number of values of gRf between min.gRf and max.gRf.
#' @param min.gSf     The minimum of gSf.
#' @param max.gSf     The maximum of gSf.
#' @param n.gSfs      The number of values of gSf between min.gSf and max.gSf.
#' @param n.cores     The number of cores in foreach loop.
#' @param verbose     A logical flag. Set it to TRUE to report timing of the foreach loop.
dx.vs.gRf.mod.gSf <- function(X,
                             min.K=15,
                             med.dK.divf=10,
                             max.dx.C=5,
                             max.dx=NULL,
                             min.gRf=3, max.gRf=20, n.gRfs=20,
                             min.gSf=1, max.gSf=5,  n.gSfs=10,
                             n.cores=10,
                             verbose=TRUE)
{
    if ( verbose ) {
        cat("\nEntering dx.vs.gRf.mod.gSf()\n")
        routine.ptm <- proc.time()
    }

    nrX <- nrow(X)

    rr <- mstree.grid(X)
    mst.grid <- rr$mst.grid
    mode.edge.len <- rr$mode.edge.len

    dx.fn <- function(dx, gSf=2, gRf=1)
    {
        X.grid <- create.3D.TN.grid(X, mst.grid, dx, mode.edge.len, gRf)
        if ( is.null(dim(X.grid)) ) {
            X.grid <- rbind(X.grid)
        }
        d <- abs(nrow(X.grid)/nrX - gSf)
        ##cat(sprintf("dx=%0.3g gRf=%.2g X.grid=%d deviation=%.3g\n", dx, gRf, nrX.grid, d))
        d
    }

    nn <- get.knn(X, k=min.K)
    dK <- nn$nn.dist[,min.K]
    med.dK <- median(dK)
    ## dx range
    min.dx <- med.dK / med.dK.divf

    if ( is.null(max.dx) ) {

        ##max.dx <- max.dx.C * min.dx
        r <- numeric(3)
        for ( i in 1:3 ) {
            r[i] <- diff(range(X[,i]))
        }

        max.dx <- sqrt(sum(r^2))
        max.dx <- max.dx / max.dx.C
    }

    gRfs <- seq(min.gRf, max.gRf, length=n.gRfs)
    gSfs <- seq(min.gSf, max.gSf, length=n.gSfs)

    if ( verbose ) {
        cat("Given gSf and gRf, dx's are found so that corresponding grid(X) has give gSf ... ")
    }
    ptm <- proc.time()
    registerDoParallel(n.cores)
    rres <- list()
    for ( i in seq(gSfs) )
    {
        gSf <- gSfs[i]
        res <- foreach ( i=seq(gRfs) ) %dopar% {
            r <- optimize(dx.fn, interval=c(min.dx, max.dx), gSf=gSf, gRf=gRfs[i])
            r$minimum
        }
        rres[[i]] <- res
    }
    stopImplicitCluster()
    if ( verbose ) {
        elapsed.time(ptm)
    }

    robust.lines <- matrix(nrow=length(rres), ncol=2)
    for ( i in seq(rres) )
    {
        dxs <- unlist(rres[[i]])
        ## if ( i == 1 ) {
        ##     idx <- log(gRfs) < 2.5
        ##     l <- line(log(gRfs[idx]), log(dxs[idx]))
        ## } else {
        l <- line(log(gRfs), log(dxs))
        ##}
        robust.lines[i,] <- l$coefficients
    }


    ##
    ## creating functions based on r1 and r2 models
    ##
    ## r1 <- mabilo(seq(gSfs), robust.lines[,1], k.min=5)
    ## r2 <- mabilo(seq(gSfs), robust.lines[,2], k.min=5)
    ## log.yint.vs.log.gSf.fn <- approxfun(r1$xg, r1$Eyg, rule=2)
    ## log.slope.vs.log.gSf.fn <- approxfun(r2$xg, r2$Eyg, rule=2)

    if ( verbose ) {
        txt <- sprintf("Total elapsed time in dx.vs.gRf.mod.gSf()")
        elapsed.time(routine.ptm, txt, with.brackets=FALSE)
    }

    list(min.K=min.K,
         med.dK.divf=med.dK.divf,
         max.dx.C=max.dx.C,
         min.gRf=min.gRf, max.gRf=max.gRf, n.gRfs=n.gRfs,
         min.gSf=min.gSf, max.gSf=max.gSf, n.gSfs=n.gSfs,
         ## estimates
         mst.grid=mst.grid,
         mode.edge.len=mode.edge.len,
         med.dK=med.dK,
         min.dx=min.dx,
         max.dx=max.dx,
         gRfs=gRfs,
         gSfs=gSfs,
         rres=rres,
         robust.lines=robust.lines)
         ## r1=r1,
         ## r2=r2,
         ## log.yint.vs.log.gSf.fn=log.yint.vs.log.gSf.fn,
         ## log.slope.vs.log.gSf.fn=log.slope.vs.log.gSf.fn)
}


#' Create an Equidistant Grid in x-Dimensional Bounding Box (C Interface)
#'
#' This function generates an equidistant grid within a specified x-dimensional bounding box.
#' The grid consists of points that are uniformly spaced, and each edge of the grid has the same length.
#'
#' @param edge.length      A positive numeric value representing the common length of each edge in the grid.
#' @param lower.bounds     A numeric vector containing the left (lower) end of the range for each dimension.
#' @param upper.bounds     A numeric vector containing the right (upper) end of the range for each dimension.
#'
#' @return A matrix representing the equidistant grid, where each row corresponds to a point in the grid,
#'         and each column corresponds to a dimension.
#'
#' @examples
#' \dontrun{
#' edge.length <- 1.0
#' lower.bounds <- c(0, 0, 0)
#' upper.bounds <- c(2, 3, 4)
#' grid <- create.ED.grid.xD(edge.length, lower.bounds, upper.bounds, round.up = TRUE)
#' print(grid)
#' }
## It provides an option to either round up or down the number of grid elements along each axis.
## @param round.up         Set to TRUE to compute the number of grid elements along the axis using
##                         \code{ceiling((upper_bounds - lower_bounds) / edge_length)}. If FALSE,
##                         \code{floor((upper_bounds - lower_bounds) / edge_length)} will be used.
##                         This determines whether the number of grid elements is rounded up or down.
create.ED.grid.xD <- function(edge.length, lower.bounds, upper.bounds ) { # , round.up = TRUE

    if (length(lower.bounds) != length(upper.bounds)) {
        stop("lower.bounds and upper.bounds must have the same length.")
    }

    if ( edge.length <= 0 ) {
        stop("Edge length must be greater than zero.")
    }

    if ( !all(lower.bounds < upper.bounds) ) {
        stop("Right end of range must be greater than left end for each dimension.")
    }

    dim <- length(upper.bounds)

    ##if ( round.up ) {
    size.d <- ceiling((upper.bounds - lower.bounds) / edge.length)
    ## } else {
    ##     size.d <- floor((upper.bounds - lower.bounds) / edge.length)
    ## }
    ## size.d[size.d == 0] <- 1

    n.grid.nodes <- prod(size.d)

    grid <- matrix(0, ncol = dim, nrow = n.grid.nodes)

    out <- .C("C_create_ED_grid_xD",
             as.double(edge.length),
             as.double(lower.bounds),
             as.integer(dim),
             as.integer(size.d),
             as.integer(n.grid.nodes),
             grid = as.double(grid))

    return(matrix(out$grid, ncol = dim, byrow = TRUE))
}

#' Creates a Uniform Grid in a Tubular Neighborhood of a State Space
#'
#' This function generates a uniform grid in a tubular neighborhood of a given state space X, considering specified parameters for grid size, radius scaling, and other factors.
#'
#' @param X          A numeric matrix representing the state space for which the grid is to be created.
#' @param gSf        A numeric value representing the grid size scaling factor, such that the grid size \eqn{n(X.grid) = gSf * n(X)}, where \eqn{n(X)} is the number of rows of X.
#' @param gRf        A numeric value representing the grid radius scaling factor. Points on the rectangular grid that are farther than \eqn{gRf * mode.edge.len} from the closest subdivision of the minimal spanning tree of X, \eqn{mstree(X)}, are eliminated from the grid.
#' @param p.exp      A numeric value representing the expansion factor of the bounding box. Must be between 0 and 1. Default is 0.05.
#' @param wC         A numeric value representing a multiplication factor controlling the precision of the estimation of the volume of X. Must be a positive number. Default is 0.1.
#' @param wF         A numeric value representing a weighting factor. Default is 1.25.
#' @param X.vol.frac The volume fraction of X in its bounding box.
#' @param max.tmp.grid.size Maximum size of temporary grid. Default is 10^7.
#' @param min.tree.obj A minimal spanning tree object. This parameter is useful when running this function with different values of gSf and gRf parameters and so only the first call may have min.tree.obj set to NULL and all remaining may use a minimal tree object from the first call.
#' @param n.threads  The number of threads to use in the corresponding C function call.
#' @param verbose    A logical value. Set to TRUE to see information on what is being done. Default is FALSE.
#'
#' @return A list with class \code{"gridX"} and the following components:
#' \itemize{
#'   \item \code{X.grid}: A numeric matrix representing a uniform grid in a tubular neighborhood of X.
#'   \item \code{d.grid}: A numeric vector representing the distance to the boundary of the grid.
#'   \item \code{mst.grid}: A numeric matrix representing a uniform grid over the edges of the minimum spanning tree of X.
#'   \item \code{mode.edge.len}: A numeric value representing the mode of the lengths of the edges of the minimum spanning tree of X.
#'   \item \code{L}: A numeric vector representing the left end of the range for each dimension.
#'   \item \code{R}: A numeric vector representing the right end of the range for each dimension.
#'   \item \code{dim}: An integer representing the number of dimensions.
#'   \item \code{X}: A numeric matrix representing the original state space.
#'   \item \code{X.vol.frac}: A numeric value representing the fraction of the volume of X in the bounding box.
#'   \item \code{X.vol}: A numeric value representing the estimated volume of X.
#'   \item \code{gRf}: A numeric value representing the grid radius scaling factor.
#'   \item \code{gSf}: A numeric value representing the grid size scaling factor.
#'   \item \code{p.exp}: A numeric value representing the expansion factor of the bounding box.
#'   \item \code{wC}: A numeric value representing the multiplication factor controlling the precision of the volume estimation.
#' }
#'
#' @examples
#' \dontrun{
#' # Let X be a low-dimensional model of a state space.
#' res <- create.X.grid.xD(X, gSf=5, gRf=5, p.exp=0.05, wC=0.1)
#' str(res)
#' }
create.X.grid.xD <- function(X,
                            gSf,
                            gRf,
                            p.exp = 0.05,
                            wC = 0.1,
                            wF = 1.25,
                            X.vol.frac = NULL,
                            max.tmp.grid.size = 10^7,
                            min.tree.obj = NULL,
                            n.threads = 1,
                            verbose = FALSE)
{
    if (p.exp < 0 || p.exp > 1) {
        stop("The expansion factor 'p.exp' must be between 0 and 1.")
    }

    if (wC <= 0) {
        stop("The multiplication factor 'wC' must be a positive number.")
    }

    nrX <- nrow(X)

    if ( verbose ) {
        fn.ptm <- proc.time()
        cat("n(X): ", nrX, "\n")
    }

    ## Target grid size
    target.n.grid <- as.integer(gSf * nrX)

    ## Extreme vertices of the bounding box
    L <- apply(X, 2, min)
    R <- apply(X, 2, max)
    dim <- length(L)

    ## Expanding the bounding box (L, R)
    r <- expand.box(L, R, p.exp)
    L <- r$L
    R <- r$R

    if ( !is.null(min.tree.obj) && !is.null(min.tree.obj$mst.grid) && ncol(min.tree.obj$mst.grid) == ncol(X)
         && !is.null(min.tree.obj$mode.edge.len) && is.numeric(min.tree.obj$mode.edge.len) ) {

        mst.grid <- min.tree.obj$mst.grid
        mode.edge.len <- min.tree.obj$mode.edge.len

    } else {

        ## Mimal spanning tree
        if ( verbose ) {
            ptm <- proc.time()
            cat("Creating a minimal spanning tree of X ... ")
        }

        min.tree.obj <- mstree.grid(X)
        mst.grid <- min.tree.obj$mst.grid
        mode.edge.len <- min.tree.obj$mode.edge.len

        if ( verbose ) {
            elapsed.time(ptm)
            cat("mode.edge.len: ", mode.edge.len, "\n")
        }
    }

    gRf.thld <- gRf * mode.edge.len

    if ( is.null(X.vol.frac) ) {

        ## Estimating the volume of X in the bounding box
        w <- wC * gRf * mode.edge.len
        if ( verbose )
            cat("w: ", w, "\n")

        size.d <- as.integer((R - L) / w) + 1
        n.tmp.grid <- prod(size.d)
        if ( verbose )
            cat("n.tmp.grid: ", n.tmp.grid, "\n")

        while ( n.tmp.grid > max.tmp.grid.size ) {
            w <- wF * w
            if ( verbose )
                cat("w: ", w, "\n")
            size.d <- as.integer((R - L) / w) + 1
            n.tmp.grid <- prod(size.d)
            if ( verbose )
                cat("n.tmp.grid: ", n.tmp.grid, "\n")
        }

        if ( verbose ) {
            ptm <- proc.time()
            cat("Creating tmp.grid ... ")
        }

        tmp.grid <- create.ED.grid.xD(w, L, R)

        if ( verbose )
            elapsed.time(ptm)


        size.d <- as.integer((R - L) / w) + 2
        n.tmp.grid <- prod(size.d)
        if ( verbose )
            cat("n.tmp.grid: ", n.tmp.grid, "\n")

        ## Identifying points of tmp.grid that are gRf * mode.edge.len way from the min spanning tree
        nn <- get.knnx(mst.grid, tmp.grid, k=1)
        d <- nn$nn.dist[,1]

        idx <- d < gRf.thld
        n.tmp.grid.around.X <- sum(idx)
        X.vol.frac <- n.tmp.grid.around.X / n.tmp.grid
        if ( verbose )
            cat("X.vol.frac: ", X.vol.frac, "\n")
    }

    ## Compute the volume of the bounding box
    bounding.box.vol <- prod(R-L)
    X.vol <- X.vol.frac * bounding.box.vol
    if ( verbose )
        cat("X.vol: ", X.vol, "\n")

    ## Estimating the Desired Grid Density
    density = target.n.grid / X.vol

    ## Estimating the Trarget Grid Width w
    p <- 1/dim
    w <- (1/density)^p
    if ( verbose )
        cat("w: ", w, "\n")

    size.d <- as.integer((R - L) / w) + 1
    est.n.grid <- prod(size.d)
    if ( verbose )
        cat("est.n.grid: ", est.n.grid, "\n")

    while ( est.n.grid > max.tmp.grid.size ) {
        w <- wF * w
        if ( verbose )
            cat("w: ", w, "\n")
        size.d <- as.integer((R - L) / w) + 1
        est.n.grid <- prod(size.d)
        if ( verbose )
            cat("est.n.grid: ", est.n.grid, "\n")
    }

    X.grid <- create.ED.grid.xD(w, L, R)
    if ( verbose )
        cat("Before filtering n(X.grid): ", nrow(X.grid), "\n")

    nn <- get.knnx(mst.grid, X.grid, k=1)
    d <- nn$nn.dist[,1]
    idx <- d < gRf.thld
    d.grid <- 1 - d[idx] / gRf.thld
    X.grid <- X.grid[idx,]
    if ( verbose ) {
        cat("After filtering n(X.grid): ", nrow(X.grid), "\n")
        cat("Target n(X.grid): ", target.n.grid, "\n")
    }

    grid.obj <- list()
    grid.obj$X.grid <- X.grid
    grid.obj$d.grid <- d.grid
    grid.obj$mst.grid <- mst.grid
    grid.obj$mode.edge.len <- mode.edge.len
    ##grid.obj$tmp.grid <- tmp.grid
    grid.obj$L <- L
    grid.obj$R <- R
    grid.obj$dim <- dim
    grid.obj$X <- X
    grid.obj$X.vol.frac <- X.vol.frac
    grid.obj$X.vol <- X.vol
    grid.obj$gRf <- gRf
    grid.obj$gSf <- gSf
    grid.obj$p.exp <- p.exp
    grid.obj$wC <- wC

    class(grid.obj) <- "gridX"

    if ( verbose )
        elapsed.time(fn.ptm, message="Total elapsed time")

    grid.obj
}

#' Create an Adaptive Tiled Grid Representation of X
#'
#' This function creates a grid representation of the input data X using an adaptive
#' tiling approach. It's particularly effective for large or complex datasets.
#'
#' @param X A numeric matrix where each row represents a point in the space.
#' @param gSf Grid size factor. Determines the target number of grid points relative to the number of input points.
#' @param gRf Grid radius factor. Controls the proximity of grid points to the minimal spanning tree of X.
#' @param n Number of boxes for the initial tiling. Default is 5.
#' @param n.itrs Number of iterations for the box tiling process. Default is 1.
#' @param p.exp Expansion factor for the bounding box. Must be between 0 and 1. Default is 0.075.
#' @param K Number of nearest neighbors to consider in the box tiling process. Default is 10.
#' @param wC Multiplication factor for temporary grid creation. Must be positive. Default is 0.1.
#' @param wF Factor for increasing grid width if necessary. Default is 1.25.
#' @param X.vol.frac Pre-computed volume fraction of X, if available. Default is NULL.
#' @param vol.box.grid.size Target size for the grid used in volume estimation. Default is 1000.
#' @param max.box.grid.size Maximum allowed size for a single box grid. Default is 10^7.
#' @param min.tree.obj Pre-computed minimal spanning tree object, if available. Default is NULL.
#' @param verbose Logical; if TRUE, print progress information. Default is FALSE.
#'
#' @return A list of class "gridX" containing:
#'   \item{X.grid}{The generated grid points}
#'   \item{d.grid}{Distance values for the grid points}
#'   \item{mst.grid}{Minimal spanning tree of X}
#'   \item{mode.edge.len}{Mode of the edge lengths in the minimal spanning tree}
#'   \item{L, R}{Lower and upper bounds of the bounding box}
#'   \item{dim}{Dimension of the space}
#'   \item{X}{Original input data}
#'   \item{X.vol.frac}{Estimated volume fraction of X}
#'   \item{X.vol}{Estimated volume of X}
#'   \item{gRf, gSf, p.exp, wC}{Input parameters}
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(1000), ncol = 2)
#' grid_obj <- create.adaptive.tiled.X.grid.xD(X, gSf = 2, gRf = 1.5, verbose = TRUE)
#' plot(grid_obj$X.grid, col = "red", pch = 20)
#' points(X, col = "blue", pch = 1)
#' }
#'
#' @seealso \code{\link{create.X.grid.xD}} for the original, non-tiled version
#'
#' @export
create.adaptive.tiled.X.grid.xD <- function(X,
                                            gSf,
                                            gRf,
                                            n = 5,
                                            n.itrs = 1,
                                            p.exp = 0.075,
                                            K = 10,
                                            wC = 0.1,
                                            wF = 1.25,
                                            X.vol.frac = NULL,
                                            vol.box.grid.size = 10^3,
                                            max.box.grid.size = 10^7,
                                            min.tree.obj = NULL,
                                            verbose = FALSE)
{
    if (p.exp < 0 || p.exp > 1) {
        stop("The expansion factor 'p.exp' must be between 0 and 1.")
    }

    if (wC <= 0) {
        stop("The multiplication factor 'wC' must be a positive number.")
    }

    nrX <- nrow(X)
    dim <- ncol(X)

    if ( verbose ) {
        fn.ptm <- proc.time()
        cat("n(X): ", nrX, "\n")
    }

    if ( !is.null(min.tree.obj) && !is.null(min.tree.obj$mst.grid) && ncol(min.tree.obj$mst.grid) == ncol(X)
         && !is.null(min.tree.obj$mode.edge.len) && is.numeric(min.tree.obj$mode.edge.len) ) {

        mst.grid <- min.tree.obj$mst.grid
        mode.edge.len <- min.tree.obj$mode.edge.len

    } else {

        ## Mimal spanning tree
        if ( verbose ) {
            ptm <- proc.time()
            cat("Creating a minimal spanning tree of X ... ")
        }

        min.tree.obj <- mstree.grid(X)
        mst.grid <- min.tree.obj$mst.grid
        mode.edge.len <- min.tree.obj$mode.edge.len

        if ( verbose ) {
            elapsed.time(ptm)
            cat("mode.edge.len: ", mode.edge.len, "\n")
        }
    }

    gRf.thld <- gRf * mode.edge.len

    ## Target grid size
    target.n.grid <- as.integer(gSf * nrX)

    ## Calculate bounding box of X
    L <- apply(X, 2, min)
    R <- apply(X, 2, max)

    ## Expand the bounding box by p.exp
    box.range <- R - L
    L <- L - p.exp * box.range
    R <- R + p.exp * box.range

    ## Creating a box coverage of the bounding box of X of approximately size
    ## 'n.boxes' and then selects from that coverage only boxes contating
    ## elements of X.

    boxes.with.elts.of.X <- box.tiling(X, n, n.itrs, p.exp , K, verbose)

    if ( verbose ) {
        if ( exists("ptm") ) elapsed.time(ptm)
        cat("n(boxes.with.elts.of.X): ", length(boxes.with.elts.of.X), "\n")
    }

    ##
    ## Estimating the volume of X in the boxes containing elements of X
    ##if ( is.null(X.vol.frac) ) {
    ##

    ## Volume of a single box
    box <- boxes.with.elts.of.X[[1]]
    box.vol <- prod(box$R - box$L)
    boxes.with.elts.of.X.vol <- length(boxes.with.elts.of.X) * box.vol

    if ( verbose ) {
        if ( exists("ptm") ) elapsed.time(ptm)
        cat("Vol(boxes.with.elts.of.X): ", boxes.with.elts.of.X.vol, "\n")
    }

    ## Within single box grid width, w, for the purpose of volume of X estimate
    w.vol <- (box.vol / vol.box.grid.size)^(1/dim)
    if ( verbose ) {
        cat("w.vol: ", w.vol, "\n")
        ptm <- proc.time()
        cat("Estimating X volume ration within the boxes containing elements of X ... ")
    }

    total.n.grid.elements <- 0
    X.n.grid.elements <- 0 # elements of the grid close to X
    for ( i in seq(length(boxes.with.elts.of.X)) ) {
        box <- boxes.with.elts.of.X[[i]]
        ids <- box$ids

        grid.within.box <- create.ED.grid.xD(w.vol, box$L, box$R)
        total.n.grid.elements <- total.n.grid.elements + nrow(grid.within.box)

        if ( 0 ) {

            cat("i: ", i, "\n")
            cat("nrow(grid.within.box): ", nrow(grid.within.box), "\n")
            cat("ncol(grid.within.box): ", ncol(grid.within.box), "\n")
            cat("ncol(X): ", ncol(X), "\n")
            cat("box\n")
            print(box)
            cat("ids: ", ids, "\n")
        }

        gX <- X[ids,]
        if ( length(ids) == 1 ) {
            gX <- rbind(gX)
        }

        nn <- get.knnx(gX, grid.within.box, k = 1)
        d <- nn$nn.dist[,1]
        idx <- d < gRf.thld
        X.n.grid.elements <- X.n.grid.elements + sum(idx)
    }

    X.vol.frac <- X.n.grid.elements / total.n.grid.elements
    X.vol <- X.vol.frac * boxes.with.elts.of.X.vol

    if ( verbose ) {
        elapsed.time(ptm)
        cat("X.vol.frac: ", X.vol.frac, "\n")
        cat("X.vol: ", X.vol, "\n")
    }

    ## Estimating the Desired Grid Density
    density = target.n.grid / X.vol

    ## Estimating the Trarget Grid Width w
    w <- (1/density)^(1/dim)
    if ( verbose )
        cat("Grid width for the target X.grid: ", w, "\n")

    est.n.grid <- box.vol / w^dim
    if ( verbose )
        cat("est.n.grid: ", est.n.grid, "\n")

    while ( est.n.grid > max.box.grid.size ) {
        w <- wF * w
        if ( verbose )
            cat("w: ", w, "\n")
        est.n.grid <- box.vol / w^dim
        if ( verbose )
            cat("est.n.grid: ", est.n.grid, "\n")
    }

    X.grid <- hgrid(X, w, gRf.thld, p.exp, K, n)

    if ( verbose ) {
        cat("After filtering n(X.grid): ", nrow(X.grid), "\n")
        cat("Target n(X.grid): ", target.n.grid, "\n")
    }

    nn <- get.knnx(mst.grid, X.grid, k=1)
    d <- nn$nn.dist[,1]
    idx <- d < gRf.thld
    d.grid <- 1 - d[idx] / gRf.thld

    # Update X.grid to only include points within threshold
    X.grid <- X.grid[idx, , drop = FALSE]

    grid.obj <- list()
    grid.obj$X.grid <- X.grid
    grid.obj$d.grid <- d.grid
    grid.obj$mst.grid <- mst.grid
    grid.obj$mode.edge.len <- mode.edge.len
    grid.obj$L <- L
    grid.obj$R <- R
    grid.obj$dim <- dim
    grid.obj$X <- X
    grid.obj$X.vol.frac <- X.vol.frac
    grid.obj$X.vol <- X.vol
    grid.obj$gRf <- gRf
    grid.obj$gSf <- gSf
    grid.obj$p.exp <- p.exp
    grid.obj$wC <- wC

    class(grid.obj) <- "gridX"

    if ( verbose )
        elapsed.time(fn.ptm, message="Total elapsed time")

    grid.obj
}

#' Plots the grid of a set around which the grid was created
#'
#' @param x An output from v1.create.X.grid.xD() or create.X.grid.xD().
#' @param with.bounding.box.pts Set to TRUE to show the bounding box.
#' @param ... Additional arguments passed to the plotting functions.
#'
plot.gridX <- function(x, with.bounding.box.pts = TRUE, ...)
{
    grid.obj <- x  # for backward compatibility
    L <- grid.obj$L
    R <- grid.obj$R
    X <- grid.obj$X
    X.grid <- grid.obj$X.grid
    mst.grid <- grid.obj$mst.grid

    if ( grid.obj$dim == 2 ) {

        xrg <- range(c(L[1], R[1], X[,1]))
        yrg <- range(c(L[2], R[2], X[,2]))
        plot(X[,1], X[,2], pch=20, las=1, xlim=xrg, ylim=yrg)
        ## for ( i in seq(nrow(T.edges)) ) {
        ##     s <- T.edges[i,1]
        ##     e <- T.edges[i,2]
        ##     segments(X[s,1], X[s,2], X[e,1], X[e,2], col="red")
        ## }
        points(mst.grid, pch=".")
        points(X.grid, pch=".")

        if ( with.bounding.box.pts ) {
            points(L[1], L[2], col="red", pch=19)
            points(R[1], R[2], col="red", pch=19)
        }

    } else if ( grid.obj$dim == 3 ) {

    } else {

    }
}


#' Construct a Hierarchical Uniform Grid Around a State Space
#'
#' This function constructs a uniform grid G(X, w, epsilon) of width
#' w in the k-dimensional Euclidean space R^k. The grid consists of
#' points that are not more than epsilon away from the closest point of
#' the given state space X.
#'
#' The algorithm first constructs a box tiling of X. For each box, elements of X
#' that are within that box are identified. The order in which the boxes are
#' tested for the presence and identity of elements of X is done starting from
#' the box with the highest density of X elements that is estimated using the
#' distance, dK, to the K-th nearest neighbor. If a box has elements of X in it,
#' it is added to a list of boxes-with-elements-of-X and the corresponding
#' elements of X are deleted from dK. The process is iterated until dK is empty.
#'
#' This process can be iterated, by performing a second round of rough
#' subdivisions of the boxes in the list of boxes-with-elements-of-X.
#'
#' When the process is finished a uniform grid of width w is created within each
#' box from the list of boxes-with-elements-of-X.
#'
#'
#' @param X         A numeric matrix representing the state space for which the grid is to be created.
#' @param w         The desired edge length of the grid. Points in the grid are spaced at intervals of w.
#' @param epsilon   The maximum allowed distance from the closest point of X for points in the grid.
#' @param p.exp     A numeric value representing the expansion factor of the bounding box. Must be between 0 and 1. Default is 0.05.
#' @param K         The number of nearest neighbors to estimate the density of X at the given point.
#' @param n.segments.per.axis   An integer number specifying the number of sub-intervals (of equal
#'                  length within each axis, but of potentially of different lengths between
#'                  different axes) of each axis that will be used to construct box tiling of
#'                  the given box.
#' @param n.itrs    The number of iterations of box sub-divisions. Default 1.
#' @param verbose   Set to TRUE for progress messages.
#'
#' @return A numeric matrix containing the coordinates of the points in the hierarchical uniform grid.
#'         Each row represents a point in the grid, and the columns correspond to the dimensions of the space.
#'
#' @examples
#' \dontrun{
#' # Generate a synthetic 2D state space
#' X <- matrix(runif(200), ncol = 2)
#' # Construct the hierarchical grid with edge length 0.1 and epsilon 0.2
#' grid <- hgrid(X, w = 0.1, epsilon = 0.2)
#' plot(X, col = "red")
#' points(grid, col = "blue", pch = 3)
#' }
#' @export
hgrid <- function(X,
                 w,
                 epsilon,
                 p.exp = 0.05,
                 K = 10,
                 n.segments.per.axis = 10,
                 n.itrs = 1,
                 verbose = TRUE)
{
    if (p.exp < 0 || p.exp > 1) {
        stop("The expansion factor 'p.exp' must be between 0 and 1.")
    }

    stopifnot(as.integer(n.itrs) == n.itrs)
    stopifnot(n.itrs > 0)

    ## Creating a box tiling of the bounding box of X by 'n^dim' boxes and then
    ## selects from that tiling only boxes contating elements of X.
    boxes.with.elts.of.X <- box.tiling(X, n.segments.per.axis, n.itrs, p.exp , K, verbose)

    ## Create grids over selected sub-boxes and filter points by distance to X
    final.grid <- create.final.grid(boxes.with.elts.of.X, X, w, epsilon)

    return(final.grid)
}

#' Expands a bounding box given by the opposite vertices  (L, R)
#'
#' @param L A numeric vector representing the left vertex that is the point of the box with the smallest coordinates among all other points of the box.
#' @param R A numeric vector representing the right vertex that is the point of the box with the largest coordinates among all other points of the box.
#' @param p.exp An expansion factor for the bounding box. Each edge of the box is scaled by the factor (1 + p.exp).
#' @param margin A numeric value specifying a fixed margin to add to the box. If NULL (default), margin is calculated based on p.exp.
#'
expand.box <- function(L, R, p.exp, margin = NULL) {

    if ( is.null(margin) ) {

        M <- (L + R) / 2 # the center of the bounding box
        L <- (1 + p.exp) * (L - M) + M
        R <- (1 + p.exp) * (R - M) + M

    } else {

        L <- L - margin
        R <- R + margin
    }

    list(L = L, R = R)
}

#' Creates a box tiling of the bounding box of some set.
#'
#' This function creates a box tiling of the bounding box of some set by
#' selecting only boxes of the box tiling containing elements of that set. Box
#' tiling is an arrangement of boxes without gaps or overlaps, meaning that two
#' boxes can intersect only along their faces.
#'
#' @param X A matrix of points in some Euclidean space.
#' @param n.segments.per.axis An integer number specifying the number of
#'     segments (sub-intervals) of each axis that will be used to construct box
#'     tiling of the given box. The segments are of equal length within each
#'     axis, but potentially of different lengths between different axes)
#' @param eps A distance threshold (non-negative real number) such that each box
#'     is expanded by eps and then tested for the presence of the elements of X.
#' @param n.itrs The number of iterations of box sub-divisions. Default 1.
#' @param p.exp An expansion factor for the bounding box. Each edge of the box
#'     is scaled by the factor (1 + p.exp), where p.exp is a non-negative real
#'     number.
#' @param plot.it Set to TRUE to see in 2D the progress of the tiling construct.
#' @param verbose Set to TRUE for progress messages.
#'
#' @return A list with two components:
#'     \itemize{
#'       \item \code{X.boxes}: A list of boxes covering X, where each box is represented by its left and right vertices.
#'       \item \code{X.vol.ratio}: A numerical vector containing estimates of X's volume ratio over the iterations within this function.
#'     }
#'
#' @examples
#' \dontrun{
#'   X <- matrix(runif(100), ncol = 2)
#'   boxes <- box.tiling(X, n = 50, eps = 0.05, p.exp = 0.05)
#' }
#' @importFrom graphics points rect
#' @export
box.tiling <- function(X, n.segments.per.axis, eps, n.itrs = 1, p.exp = 0.05, plot.it = FALSE, verbose = FALSE) {

    ## Finding the bounding box of X
    L <- apply(X, 2, min)
    R <- apply(X, 2, max)
    dim <- length(L)

    ## Expanding the bounding box (L, R)
    r <- expand.box(L, R, p.exp)
    L <- r$L
    R <- r$R

    X.vol <- prod(R - L)

    if ( verbose ) {
        ptm <- proc.time()
        cat("Creating box tiling of X ... ")
    }

    bbox <- list() # bounding box
    bbox$L <- L
    bbox$R <- R

    ##boxes <- box.tiling.of.box(n.segments.per.axis, bbox)

    ## Creating a covering of the bounding box with hypercubes with the edge
    ## length equal to the half of the shortest edge of the bounding box
    bbox.edge.lengths <- R - L
    bbox.min.edge.length <- min(bbox.edge.lengths)
    w <- bbox.min.edge.length / 2
    boxes <- bbox.hcube.tiling(w, bbox$L, bbox$R)

    n.single.division.boxes <- n.segments.per.axis^dim

    if ( verbose ) {
        elapsed.time(ptm)
        n.boxes <- length(boxes)
        cat("itr = 1: n(boxes): ", n.boxes, "\n")
        cat("Scanning ",n.boxes," boxes for the elements of X ... ")
    }

    find.X.boxes <- function(boxes, eps, X) {

        X.boxes <- list()
        for ( i in seq(length(boxes)) ) {

            box <- boxes[[i]]

            within.box.ids <- find.points.within.box(X, box, eps)

            if ( length(within.box.ids) > 0 ) {
                box$ids <- within.box.ids
                boxes[[i]] <- box
                X.boxes[[length(X.boxes) + 1]] <- box
            }
        }

        return(X.boxes)
    }

    X.boxes <- find.X.boxes(boxes, eps, X)

    if ( plot.it ) {
        ##plot.box.tiling(X.boxes, L, R, with.box.ids = TRUE)
        plot.box.tiling(X.boxes, with.box.ids = TRUE)
        points(X, pch=20, cex=0.5)
    }


    box.tiling.volume <- function(boxes) {

        ## this can be simplified as all boxes of the tiling are of the same size
        ## for ( box in boxes ) {
        ##     vol <- vol + prod(box$R - box$L)
        ## }

        box <- boxes[[1]]
        box.vol <- prod(box$R - box$L)
        vol <- box.vol * length(boxes)

        return(vol)
    }

    X.vol.ratio <- c()
    X.vol.ratio[1] <- box.tiling.volume(X.boxes) / X.vol

    n.X.boxes <- c()
    n.X.boxes[1] <- length(X.boxes)

    X.boxes.list <- list()
    X.boxes.list[[1]] <- X.boxes

    if ( verbose ) {
        elapsed.time(ptm)
        cat("itr = 1: n(X.boxes): ", length(X.boxes), " out of ",length(boxes),"\n")
        cat("itr = 1: X.vol.ratio: ", X.vol.ratio, "\n")
        ## take one box and compute the ratio min(edge) / eps, max(edge) / eps
    }

    itr <- 2
    while ( n.itrs > 1 ) {

        if ( verbose ) {
            k <- n.single.division.boxes * length(X.boxes)
            k <- format(k, big.mark = ",")
            cat("itr: ",itr," - Scanning ",k," boxes for the elements of X ... ")
        }

        sub.X.boxes <- list()

        for ( i in seq(length(X.boxes)) ) {

            box <- X.boxes[[i]]
            ids <- box$ids # getting ids of the points of X that are in the given box

            ## subdividing the box into sub-boxes
            sub.boxes <- box.tiling.of.box(n.segments.per.axis, box)

            if ( plot.it ) {
                for ( j in seq(sub.boxes) ) {
                    box <- sub.boxes[[j]]
                    rect(box$L[1], box$L[2], box$R[1], box$R[2], border = "gray")
                }
            }

            ## looking for sub-boxes of 'box' that (after eps-expansion) contain elements of X
            X.sub.boxes <- find.X.boxes(sub.boxes, eps, X)

            if ( plot.it ) {
                for ( j in seq(X.sub.boxes) ) {
                    box <- X.sub.boxes[[j]]
                    rect(box$L[1], box$L[2], box$R[1], box$R[2], border = "blue")
                }
            }

            sub.X.boxes <- c(sub.X.boxes, X.sub.boxes)
        }

        X.boxes <- sub.X.boxes
        X.boxes.list[[itr]] <- X.boxes

        X.vol.ratio[itr] <- box.tiling.volume(X.boxes) / X.vol
        n.X.boxes[itr] <- length(X.boxes)

        if ( verbose ) {
            elapsed.time(ptm)
            k <- n.segments.per.axis^{dim * itr}
            k <- format(k, big.mark = ",")
            cat("itr = ",itr,": n(X.boxes): ", length(X.boxes), " out of ",k,"\n")
            cat("itr = ",itr,": X.vol.ratio: ", X.vol.ratio, "\n")
            ## take one box and compute the ratio min(edge) / eps, max(edge) / eps
        }

        itr <- itr + 1
        n.itrs <- n.itrs - 1
    }

    list(X.boxes = X.boxes,
         X.boxes.list = X.boxes.list,
         X.vol.ratio = X.vol.ratio,
         n.X.boxes = n.X.boxes,
         n.segments.per.axis = n.segments.per.axis,
         eps = eps,
         n.itrs = n.itrs,
         p.exp = p.exp,
         bbox = bbox)
}


#' Extracts Numerical x and y Values from a Character Vector.
#'
#' This function takes a character vector of the form \code{c("x1,y1", "x2,y2", ..., "xN,yN")}
#' and extracts two numerical vectors, \code{x} and \code{y}, containing the corresponding
#' x and y values.
#'
#' @param s A character vector containing strings of the form "xi,yi". Each string must
#'   contain exactly one comma, and the values before and after the comma must be numerical.
#'
#' @return A list containing two numerical vectors:
#' \itemize{
#'   \item \code{x}: A vector containing the x values (c(x1, x2, ..., xN)).
#'   \item \code{y}: A vector containing the y values (c(y1, y2, ..., yN)).
#' }
#'
#' @examples
#' \dontrun{
#' s <- c("1,2", "3,4", "5,6")
#' result <- extract_xy(s)
#' x <- result$x
#' y <- result$y
#' }
#' @export
extract.xy <- function(s) {
    ## Split each string in the vector by the comma
    split.values <- strsplit(s, split = ",")

    ## Extract x and y values
    x <- sapply(split.values, function(val) as.numeric(val[1]))
    y <- sapply(split.values, function(val) as.numeric(val[2]))

    ## Return as a list
    list(x = x, y = y)
}

#' Creates box tiling of a box.
#'
#' This function constructs a box tiling of a box specified by vectors 'L' and 'R'. The boxes forming the tiling are
#' created by subdividing each axis into the same number of equal length intervals. A box is (also called hyperrectangle
#' or orthotope) is defined as the Cartesian product of k intervals, where k is the dimension of the orthotope. A box is
#' defined by two opposite vertices L and R.
#'
#' @param n An integer number specifying the number of sub-intervals (of equal
#'     length within each axis, but of potentially of different lengths between
#'     different axes) of each axis that will be used to construct box tiling of
#'     the given box.
#' @param box         A box with L, R, i components.
#'
#' @return A list of sub-boxes, where each sub-box is represented as a list with components 'L' and 'R', containing
#' the left and right vertices of the corresponding sub-box.
#'
#' @examples
#' \dontrun{
#' w <- 1
#' box <- list()
#' box$L <- c(0, 0)
#' box$R <- c(2, 3)
#' box$i <- 1
#' sub_boxes <- box.tiling.of.box(w, box)
#' }
#'
#' @seealso \code{\link{extract.xy}} for extracting x and y values from a character vector.
#' @export
## @param edge.length The length of the edge of a hypercube that will be use for the tiling of the box. The tiling will be most likely going ouside of the right faces of the box.
## create.ED.grid.boxes <- function(w, L, R)
box.tiling.of.box <- function(n, box) { #, edge.length = NULL) {

    L <- box$L
    R <- box$R

    ## Ensure that L and R have the same length as dims
    if (length(box$L) != length(R) ) {
        stop("L and R must have the same length.")
    }

    if ( as.integer(n) != n ) {
        stop("n has to be an integer.")
    }

    if ( n <= 1 ) {
        stop("n has to be greater than 1.")
    }

    dim <- length(L)

    ## Create corresponding intervals
    BBox.Is <- list()
    for ( d in seq(dim) ) {

        ## if ( !is.null(edge.length) && is.numeric(edge.length) && edge.length > 0 ) {
        ##     n <- ceiling( (R[d] - L[d]) / edge.length )
        ##     x <- L[d] + edge.length * 0:n
        ## } else {

        ## }
        x <- seq(L[d], R[d], length.out = n + 1)
        edge.Is <- c()
        for ( i in seq(n) ) {
            edge.Is[i] <- paste0(x[i],",",x[i+1])
        }
        BBox.Is[[d]] <- edge.Is
    }

    BBox.edge.seqs <- expand.grid(BBox.Is)
    BBox.edge.seqs <- apply(BBox.edge.seqs, 2, FUN=as.character)

    BBox.boxes <- list()
    for ( i in seq(nrow(BBox.edge.seqs)) ) {

        s <- BBox.edge.seqs[i,]
        names(s) <- NULL
        r <- extract.xy(s)
        BBox.boxes[[i]] <- list(L = r$x, R = r$y, i = ifelse(!is.null(box$i), paste0(box$i,":",i), i), ids = box$ids)
    }

    return(BBox.boxes)
}

#' Creates a hypercube tiling of a bounding box
#'
#' This function constructs a hypercube tiling of a bounding box using hypercubes of width 'w'.
#'
#' @param w A numeric value representing the width of each sub-box, i.e., the edge length.
#' @param L A numeric vector representing the left vertices of the bounding box.
#' @param R A numeric vector representing the right vertices of the bounding box.
#'
#' @return A list of sub-boxes, where each sub-box is represented as a list with components 'L' and 'R', containing
#' the left and right vertices of the corresponding sub-box.
#'
#' @details
#' The function first ensures that 'L' and 'R' have the same length and then calculates the number of grid elements along
#' each coordinate. It then creates corresponding intervals and constructs the sub-boxes within the bounding box.
#'
#' @examples
#' \dontrun{
#' w <- 1
#' L <- c(0, 0)
#' R <- c(2, 3)
#' sub_boxes <- bbox.hcube.tiling(w, L, R)
#' }
#'
#' @seealso \code{\link{extract.xy}} for extracting x and y values from a character vector.
#' @export
## create.ED.grid.boxes <- function(w, L, R) {
bbox.hcube.tiling <- function(w, L, R) {

    ## Ensure that L and R have the same lengths
    if (length(L) != length(R) ) {
        stop("L and R must have the same length")
    }

    dim <- length(L)

    ## Calculate the number of intervals of length w that will be covering each coordinate
    n.intervals <- ceiling((R - L) / w)

    ## Create corresponding intervals
    BBox.Is <- list()
    for ( d in seq(dim) ) {

        x <- L[d] + w * 0:n.intervals[d]
        edge.Is <- c()
        for ( i in seq(n.intervals[d]) ) {
            edge.Is[i] <- paste0(x[i],",",x[i+1])
        }
        BBox.Is[[d]] <- edge.Is
    }

    BBox.edge.seqs <- expand.grid(BBox.Is)
    BBox.edge.seqs <- apply(BBox.edge.seqs, 2, FUN=as.character)

    BBox.boxes <- list()
    for ( i in seq(nrow(BBox.edge.seqs)) ) {

        s <- BBox.edge.seqs[i,]
        names(s) <- NULL
        r <- extract.xy(s)
        BBox.boxes[[i]] <- list(L=r$x, R=r$y)
    }

    return(BBox.boxes)
}

## Creates a hypercube tiling of a bounding box. An R interfe to a C function.
##
## This function constructs a hypercube tiling of a bounding box using hypercubes of width 'w'.
##
## @param w A numeric value representing the width of each sub-box, i.e., the edge length.
## @param L A numeric vector representing the left vertices of the bounding box.
## @param R A numeric vector representing the right vertices of the bounding box.
##
## @return A list of sub-boxes, where each sub-box is represented as a list with components 'L' and 'R', containing
## the left and right vertices of the corresponding sub-box.
##
## @details This function creates a tiling of a bounding box defined by the vertices L and R using hypercubes of width
##     w. The result is a list of hypercubes represented as an m-by-2*n matrix. Each row of the matrix corresponds to a
##     sub-hypercube, where the first half of the row corresponds to the 'L' values and the second half to the 'R'
##     values.
##
##   For example, given the inputs:
##   w = 1.0
##   L = [0, 0]
##   R = [2, 3]
##
##   The output matrix would be:
##   | 0 | 0 | 1 | 1 |
##   | 0 | 1 | 1 | 2 |
##   | 0 | 2 | 1 | 3 |
##   | 1 | 0 | 2 | 1 |
##   | 1 | 1 | 2 | 2 |
##   | 1 | 2 | 2 | 3 |
##
## @examples
## \dontrun{
## w <- 1
## L <- c(0, 0)
## R <- c(2, 3)
## sub_boxes <- C.bbox.hcube.tiling(w, L, R)
## }
##
## @seealso \code{\link{extract.xy}} for extracting x and y values from a character vector.
## C.bbox.hcube.tiling <- function(w, L, R) {

##     ## Ensure that L and R have the same lengths
##     if (length(L) != length(R) ) {
##         stop("L and R must have the same length")
##     }

##     dim <- length(L)

##     ## The number of intervals into which each coordinate is split into
##     n.intervals <- ceiling((R - L) / w)

##     n.hcubes <- prod(n.intervals)

##     hcubes <- numeric(n.hcubes * 2 * dim)

##     res <- .C("C_bbox_hcube_tiling",
##              as.double(w),
##              as.double(L),
##              as.double(R),
##              as.integer(dim),
##              as.integer(n.hcubes),
##              hcubes = hcubes)

##     ## Turning res$hcubes matrix into a list
##     hcubes.list <- vector("list", n.hcubes)
##     for (i in 1:n.hcubes) {
##         L.component <- res$hcubes[i, 1:dim]
##         R.component <- res$hcubes[i, (dim+1):(2*dim)]
##         hcubes.list[[i]] <- list(L = L.component, R = R.component)
##     }

##     return(hcubes.list)
## }

#' Identify Points Within a Box
#'
#' This function identifies the points within a given box specified by the lower and upper bounds.
#'
#' @param X A matrix or data frame containing the points in the space.
#' @param L A numeric vector specifying the lower bounds of the box.
#' @param R A numeric vector specifying the upper bounds of the box.
#'
#' @return A matrix or data frame containing only the points that fall within the specified box.
#'
#' @examples
#' \dontrun{
#' X <- matrix(runif(100), ncol = 2)
#' L <- c(0.2, 0.2)
#' R <- c(0.5, 0.5)
#' points_within_box <- identify_points_within_box(X, L, R)
#' }
#'
#' @export
identify_points_within_box <- function(X, L, R) {

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

    if (length(L) != ncol(X) || length(R) != ncol(X)) {
        stop("Dimensions of L and R must match the number of columns in X.")
    }

    within_box <- apply(X, 1, function(row) all(L < row & row <= R))

    return(names(within_box)[within_box])
}

#' Find the points of a set X within a box.
#'
#' This function finds ids of the points of X that are found within the given box.
#'
#' @param box A matrix or data frame containing the points in the space.
#' @param X   A numeric matrix corresponding to a set of points in some Euclidean space.
#' @param eps A distance threshold (non-negative real number) such that each box is expanded by eps and then tested for the presence of the elements of X.
#'
#' @return A vector of row names of X corresponding to the points fund within the given box.
#'
#' @examples
#' \dontrun{
#' X <- matrix(runif(100), ncol = 2)
#' box <- list()
#' box$L <- c(0.2, 0.2)
#' box$R <- c(0.5, 0.5)
#' points_within_box <- find.points.within.box(X, box, eps = 0)
#' }
#'
#' @export
find.points.within.box <- function(X, box, eps) {

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

    if (length(box$L) != ncol(X) || length(box$R) != ncol(X)) {
        stop("Dimensions of L and R must match the number of columns in X.")
    }

    if ( eps > 0 ) {
        r <- expand.box(box$L, box$R, margin = eps)
        L <- r$L
        R <- r$R
    } else {
        L <- box$L
        R <- box$R
    }

    ##if ( is.null(box$ids) || ( !is.null(box$ids) && length(box$ids) == 0 ) ) {
    if ( 1 ) {

        within.box <- apply(X, 1, function(row) all(L < row & row <= R))

    } else {

        box.X <- X[box$ids,]
        if ( length(box$ids) == 1 ) {
            box.X <- rbind(box.X)
        }

        within.box <- apply(box.X, 1, function(row) all(L < row & row <= R))
    }


    return(names(within.box)[within.box])
}


#' Finds the box within a box list covering some data set that contains a given point.
#'
#' This function finds the sub-box within a given set of boxes that contains a specified point \code{x}.
#'
#' @param boxes A list of boxes, where each box is a list containing vectors \code{L} and \code{R} representing the left and right boundaries of the box.
#' @param x A vector representing the point for which the containing sub-box is to be found.
#'
#' @return A list representing the sub-box containing \code{x}, or \code{NULL} if \code{x} is not found within any sub-box.
#' @examples
#' \dontrun{
#' boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
#' x <- c(0.2, 0.3)
#' containing_box <- find.box.containing.x(boxes, x)
#' }
#' @export
find.box.containing.x <- function(boxes, x) {

    for ( i in seq(length(boxes)) ) {
        box <- boxes[[i]]
        within.box <- all(mapply(function(Li, Ri, xi) Li < xi & xi <= Ri, box$L, box$R, x))
        if (within.box) {
            ## box$i <- i
            return(box)
        }
    }

    # Returns NULL if x is not found within any sub-box
    return(NULL)
}

#' Finds boxes within a box list covering some data set that (after expansion by a factor p.exp) contain a given point.
#'
#' @param boxes A list of boxes, where each box is a list containing vectors \code{L} and \code{R} representing the left and right boundaries of the box.
#' @param x     A vector representing the point for which the containing sub-box is to be found.
#' @param p.exp An expansion factor for the bounding box. Each edge of the box is scaled by the factor (1 + p.exp). Default is 0.1.
#' @param margin An explicitly expansion amount of the box in along each axis.
#'
#' @return A list representing the sub-box containing \code{x}, or \code{NULL} if \code{x} is not found within any sub-box.
#' @examples
#' \dontrun{
#' boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
#' x <- c(0.2, 0.3)
#' containing_box <- find.box.containing.x(boxes, x)
#' }
#' @export
find.boxes.containing.x <- function(boxes, x, p.exp = 0.1, margin = NULL) {

    if ( !is.null(margin) && length(margin) != length(x) ) {
        stop("margin has to have the same length as x.")
    }

    x.boxes <- list()
    for ( i in seq(length(boxes)) ) {

        box <- boxes[[i]]
        r <- expand.box(box$L, box$R, p.exp, margin)
        L <- r$L
        R <- r$R

        within.box <- all(mapply(function(Li, Ri, xi) Li < xi & xi <= Ri, L, R, x))

        if (within.box) {
            ##box$i <- i
            x.boxes[[length(x.boxes) + 1]] <- box
        }
    }

    return(x.boxes)
}


#' Creates a Final Grid within Selected Boxes
#'
#' This function creates a uniform grid within a set of selected sub-boxes, then filters the grid to remove points that are farther than a given distance \code{epsilon} from the original dataset \code{X}, and removes duplicated grid elements.
#'
#' @param boxes A list of selected sub-boxes, each represented as a list containing vectors \code{L} and \code{R} for left and right boundaries.
#' @param X A matrix representing the original dataset for which the grid is to be created.
#' @param w A numeric value representing the width of the grid.
#' @param epsilon A numeric value representing the maximum distance from \code{X} that a grid point can be to be included in the final grid.
#'
#' @return A matrix representing the final filtered grid.
#' @examples
#' \dontrun{
#' boxes <- create.ED.boxes(0.5, c(0,0), c(1,1))
#' X <- matrix(runif(20), ncol = 2)
#' w <- 0.1
#' epsilon <- 0.5
#' final_grid <- create.final.grid(boxes, X, w, epsilon)
#' }
#' @export
create.final.grid <- function(boxes, X, w, epsilon)
{
    stopifnot(is.list(boxes))
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(w))
    stopifnot(length(w) == 1)
    stopifnot(is.numeric(epsilon))
    stopifnot(length(epsilon) == 1)

    final.grid <- c() ## Final grid
    ## Loop through the selected sub-boxes
    for (box in boxes) {
        ## Create a uniform grid within the sub-box using the provided width 'w'
        grid.within.box <- create.ED.grid.xD(w, box$L, box$R)
        ## Add the filtered grid to the final grid
        final.grid <- rbind(final.grid, as.matrix(grid.within.box))
    }
    ## Check the distance to the closest point of 'X' and filter by 'epsilon'
    nn <- get.knnx(X, final.grid, k = 1)
    d <- nn$nn.dist[,1]
    idx <- d < epsilon
    filtered.grid <- final.grid[idx,]

    ## By the way grid within a box is created, in the case of points on the
    ## boundary of two boxes, there are grid elements from adjacent boxes
    ## created that are closer to each other than the width of the grid. The
    ## following function removes one of the duplicated grid elements for each
    ## duplicate pair.
    idices.to.rm <- rm.duplicate.grid.elements(filtered.grid, w, w.factor = 0.9)
    if ( length(idices.to.rm) > 0 ) {
        filtered.grid <- filtered.grid[setdiff(seq(nrow(filtered.grid)), idices.to.rm), ]
    }

    return(filtered.grid)
}



#' Plots boxes of a box tiling of a 2D state space.
#'
#' This function creates a plot showing the edges of boxes of the given box
#' tiling and their indices in the center of each box. It is designed to work
#' with the output of the \code{box.tiling} function.
#'
#' @param x              A list of sub-boxes, as returned by \code{create.ED.boxes}.
#' @param L              The lower left vertex of the bounding box  whose grid is plotted. Default is NULL.
#' @param R              The upper right vertex of the bounding box  whose grid is plotted. Default is NULL.
#' @param with.box.ids   Set to TRUE if box IDs are to be displayed.
#' @param BBox.col       A color to be used for the edges of bounding box whose grid is plotted. Default is "red".
#' @param col            A color to be used for the edges of sub-boxes. Default is "gray".
#' @param xlab           Label for x-axis. Default is "".
#' @param ylab           Label for y-axis. Default is "".
#' @param sleep          The number of seconds to sleep before drawing the next box.
#' @param ...            Additional graphical parameters to be passed to the plotting functions.
#'
#' @examples
#' \dontrun{
#' w <- 1
#' L <- c(0, 0)
#' R <- c(2, 3)
#' grid_boxes <- create.ED.boxes(w, L, R)
#' plot.boxes(grid_boxes, L, R)
#' }
#'
#' @export
## was plot.boxes
plot.box.tiling <- function(x,
                           L = NULL,
                           R = NULL,
                           with.box.ids = TRUE,
                           BBox.col = "red",
                           col = "gray",
                           xlab = "",
                           ylab = "",
                           sleep = 0,
                           ...)
{
    boxes <- x  # for backward compatibility
    xmax <- numeric(length(boxes))
    ymax <- xmax
    xmin <- xmax
    ymin <- xmax

    for (i in seq(boxes) ) {
        xmax[i] <- boxes[[i]]$R[1]
        ymax[i] <- boxes[[i]]$R[2]
        xmin[i] <- boxes[[i]]$L[1]
        ymin[i] <- boxes[[i]]$L[2]
    }
    xmax <- max(xmax)
    ymax <- max(ymax)
    xmin <- min(xmin)
    ymin <- min(ymin)

    xlim <- c(xmin, xmax)
    ylim <- c(ymin, ymax)

    plot(0, 0, type = "n", las = 1, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)

    if ( !is.null(L) && !is.null(R) && length(L) == length(R) ) {
        rect(L[1], L[2], R[1], R[2], border = BBox.col)
    }

    for (i in seq(boxes) ) {

        box <- boxes[[i]]

        if ( sleep ) {
            Sys.sleep(sleep)
        }

        rect(box$L[1], box$L[2], box$R[1], box$R[2], border = col)

        if ( with.box.ids ) {
            if ( !is.null(box[["i"]]) ) {
                graphics::text(mean(c(box$L[1], box$R[1])), mean(c(box$L[2], box$R[2])), labels = box[["i"]], col = col)
            } else {
                graphics::text(mean(c(box$L[1], box$R[1])), mean(c(box$L[2], box$R[2])), labels = i, col = col)
            }
        }
    }
}
