edge.cosine.reference <- function(a, w, y, z, type = "derivative", log = FALSE, epsilon = 0.2) {
    if (log) { y <- base::log(y + epsilon); z <- base::log(z + epsilon) }
    vapply(seq_along(a), function(v) {
        dy <- y[a[[v]]] - y[v]; dz <- z[a[[v]]] - z[v]
        weights <- if (type == "derivative") ifelse(w[[v]] > 1e-10, 1 / w[[v]]^2, 0) else rep(1, length(dy))
        ny <- sqrt(sum(weights * dy^2)); nz <- sqrt(sum(weights * dz^2))
        if (ny <= 1e-10 || nz <= 1e-10) 0 else sum(weights * dy * dz) / (ny * nz)
    }, numeric(1))
}

test_that("edge cosine agrees with independent formulas and an anisotropic star", {
    a <- list(c(2L,3L),1L,1L,integer())
    w <- list(c(1,sqrt(2)),1,sqrt(2),numeric())
    y <- c(0,1,1,5); z <- c(0,0,1,2)
    for (type in c("derivative","unit","sign")) {
        expect_equal(as.numeric(lcor(a,w,y,z,type=type)), edge.cosine.reference(a,w,y,z,type), tolerance=1e-12)
        expect_equal(as.numeric(lcor(a,w,3*y+7,2*z-1,type=type)), edge.cosine.reference(a,w,y,z,type), tolerance=1e-12)
        expect_equal(as.numeric(lcor(a,w,y,-z,type=type)), -edge.cosine.reference(a,w,y,z,type), tolerance=1e-12)
    }
    expect_equal(as.numeric(lcor(a,w,y,z))[1], 1/sqrt(3), tolerance=1e-12)
    expect_equal(as.numeric(lcor(a,lapply(w,`*`,0.01),y,z))[1], 1/sqrt(3), tolerance=1e-12)
    expect_equal(as.numeric(lcor(a,w,y,rep(2,4))), rep(0,4))
    expect_equal(as.numeric(lcor(a,w,y,y)), c(1,1,1,0))
})

test_that("matrix paths agree for specified difference and pseudocount contracts", {
    a <- list(c(2L,3L),c(1L,3L),c(1L,2L),integer())
    w <- list(c(1,2),c(1,3),c(2,3),numeric())
    y <- c(0,1,4,2); Z <- cbind(a=c(1,1,2,0),b=c(2,0,5,1))
    for (diff in c("difference","logratio")) {
        expected <- vapply(seq_len(ncol(Z)),function(j) edge.cosine.reference(a,w,y,Z[,j],log=diff=="logratio"), numeric(4))
        vm <- lcor(a,w,y,Z,y.diff.type=diff,z.diff.type=diff,epsilon=.2)
        expect_equal(unname(as.matrix(unclass(vm))), expected,tolerance=1e-12,ignore_attr=TRUE)
        for(j in 1:2) expect_equal(as.numeric(lcor(a,w,y,Z[,j],y.diff.type=diff,z.diff.type=diff,epsilon=.2)),expected[,j],tolerance=1e-12)
    }
    expect_error(lcor(a,w,-y,Z,y.diff.type="logratio"),"nonnegative")
    expect_error(lcor(a,w,y,Z,epsilon=NA),"epsilon")
    expect_error(lcor(a,w,y,Z,winsorize.quantile=.5),"winsorize")
})

test_that("directed slopes recover a known affine response", {
    a <- list(2L,c(1L,3L),2L); w <- list(1,c(1,1),1); y <- c(1,2,4)
    ans <- lslope(a,w,y,3*y+2,y.diff.type="difference",z.diff.type="difference")
    expect_equal(as.numeric(ans)[1:2],c(3,3),tolerance=1e-12)
})

test_that("selected pairs match complete arrays and preserve requested order", {
    a <- list(2L, c(1L,3L), c(2L,4L), 3L)
    w <- lapply(a, function(v) rep(1, length(v)))
    y <- cbind(a = c(0,1,3,2), b = c(3,2,1,0))
    z <- cbind(c = c(1,4,2,3), d = c(4,3,2,1), e = c(2,1,4,3))
    pairs <- rbind(c(2,3), c(1,2), c(2,3))
    for (difference in c("difference", "logratio")) {
        full <- lcor(a,w,y,z,y.diff.type=difference,z.diff.type=difference,epsilon=.2)
        selected <- lcor(a,w,y,z,y.diff.type=difference,z.diff.type=difference,epsilon=.2,pairs=pairs)
        for (k in 1:3) expect_equal(unname(selected[,k]), unname(full[,pairs[k,1],pairs[k,2]]))
        blocks <- lapply(1:3, function(k) lcor(a,w,y,z,y.diff.type=difference,z.diff.type=difference,epsilon=.2,pairs=pairs[k,,drop=FALSE]))
        expect_equal(unname(do.call(cbind, blocks)), unname(as.matrix(selected)), ignore_attr=TRUE)
        expect_identical(attr(selected,"pairs"),pairs)
        expect_identical(colnames(selected),c("b:e","a:d","b:e"))
    }
    expect_error(lcor(a,w,y,z,pairs=matrix(c(0,1),1)),"valid y/z")
    expect_error(lcor(a,w,y,z[,1],pairs=matrix(c(1,1),1)),"two matrix")
})

test_that("isotropic stars recover gradient cosine and shrinking smooth-field limits", {
    a <- list(2:5,1L,1L,1L,1L)
    xy <- rbind(c(0,0),c(1,0),c(-1,0),c(0,1),c(0,-1))
    errors <- vapply(c(.2,.02,.002), function(h) {
        X <- h*xy; w <- lapply(a,function(v) rep(h,length(v)))
        y <- X[,1]+X[,1]^2
        z <- X[,1]+X[,2]+2*X[,2]^2
        abs(as.numeric(lcor(a,w,y,z))[1]-1/sqrt(2))
    },numeric(1))
    expect_true(all(diff(errors)<0))
    expect_lt(tail(errors,1),1e-5)
    w <- lapply(a,function(v) rep(1,length(v)))
    expect_equal(as.numeric(lcor(a,w,xy[,1],xy[,1]+xy[,2]))[1],1/sqrt(2),tolerance=1e-12)
    w[[1]][1] <- 0; w[[2]] <- 0
    expect_equal(as.numeric(lcor(a,w,xy[,1],xy[,1]+xy[,2])),
                 edge.cosine.reference(a,w,xy[,1],xy[,1]+xy[,2]),tolerance=1e-12)
})
