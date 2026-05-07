#!/usr/bin/env Rscript

script.path <- if (length(commandArgs(trailingOnly = FALSE))) {
    marker <- "--file="
    file.arg <- commandArgs(trailingOnly = FALSE)
    file.arg <- file.arg[startsWith(file.arg, marker)]
    if (length(file.arg)) sub(marker, "", file.arg[[1]]) else getwd()
} else {
    getwd()
}
report.dir <- if (file.exists(script.path)) dirname(normalizePath(script.path)) else getwd()
figure.dir <- file.path(report.dir, "figures")
dir.create(figure.dir, showWarnings = FALSE, recursive = TRUE)

pairwise.dist <- function(X) {
    n <- nrow(X)
    D <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            d <- sqrt(sum((X[i, , drop = TRUE] - X[j, , drop = TRUE])^2))
            D[i, j] <- d
            D[j, i] <- d
        }
    }
    D
}

knn.index <- function(D, k) {
    n <- nrow(D)
    out <- matrix(NA_integer_, n, k)
    for (i in seq_len(n)) {
        d <- D[i, ]
        d[[i]] <- Inf
        out[i, ] <- order(d, seq_len(n))[seq_len(k)]
    }
    out
}

support.sknn <- function(D, k) {
    idx <- knn.index(D, k)
    n <- nrow(D)
    A <- matrix(FALSE, n, n)
    for (i in seq_len(n)) A[i, idx[i, ]] <- TRUE
    A <- A | t(A)
    diag(A) <- FALSE
    A
}

support.mknn <- function(D, k) {
    idx <- knn.index(D, k)
    n <- nrow(D)
    A <- matrix(FALSE, n, n)
    for (i in seq_len(n)) A[i, idx[i, ]] <- TRUE
    A <- A & t(A)
    diag(A) <- FALSE
    A
}

support.iknn <- function(D, k) {
    idx <- knn.index(D, k)
    n <- nrow(D)
    closed <- vector("list", n)
    for (i in seq_len(n)) closed[[i]] <- sort(unique(c(i, idx[i, ])))
    A <- matrix(FALSE, n, n)
    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            A[i, j] <- length(intersect(closed[[i]], closed[[j]])) > 0L
            A[j, i] <- A[i, j]
        }
    }
    A
}

support.phate <- function(D, k, a = 40, tau = 1e-4) {
    idx <- knn.index(D, k)
    sigma <- vapply(seq_len(nrow(D)), function(i) D[i, idx[i, k]], numeric(1))
    r.tau <- (-log(tau))^(1 / a)
    S <- outer(sigma, sigma, pmax)
    A <- D <= r.tau * S
    diag(A) <- FALSE
    A
}

support.fixed <- function(D, radius) {
    A <- D <= radius
    diag(A) <- FALSE
    A
}

support.adaptive <- function(D, k.scale, radius.factor = 1, radius.rule = "max") {
    idx <- knn.index(D, k.scale)
    sigma <- vapply(seq_len(nrow(D)), function(i) D[i, idx[i, k.scale]], numeric(1))
    f <- if (identical(radius.rule, "min")) pmin else pmax
    S <- outer(sigma, sigma, f)
    A <- D <= radius.factor * S
    diag(A) <- FALSE
    A
}

plot.edges <- function(X, A, edge.col = "#6b7280", edge.lwd = 0.5, edge.lty = 1) {
    ij <- which(upper.tri(A) & A, arr.ind = TRUE)
    if (!nrow(ij)) return(invisible())
    for (r in seq_len(nrow(ij))) {
        i <- ij[r, 1]
        j <- ij[r, 2]
        lines(X[c(i, j), 1], X[c(i, j), 2], col = edge.col, lwd = edge.lwd, lty = edge.lty)
    }
}

plot.graph <- function(X, A, main, point.col = "#111827", edge.col = "#6b7280") {
    plot(X[, 1], X[, 2], type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)
    plot.edges(X, A, edge.col = edge.col, edge.lwd = 0.6)
    points(X[, 1], X[, 2], pch = 19, cex = 0.65, col = point.col)
    box(col = "#d1d5db")
}

connected.components <- function(A) {
    n <- nrow(A)
    comp <- rep(NA_integer_, n)
    comp.id <- 0L
    for (start in seq_len(n)) {
        if (!is.na(comp[[start]])) next
        comp.id <- comp.id + 1L
        comp[[start]] <- comp.id
        stack <- start
        while (length(stack)) {
            u <- stack[[length(stack)]]
            stack <- stack[-length(stack)]
            nbrs <- which(A[u, ])
            for (v in nbrs) {
                if (is.na(comp[[v]])) {
                    comp[[v]] <- comp.id
                    stack <- c(stack, v)
                }
            }
        }
    }
    list(id = comp, n = comp.id)
}

component.mst.edges <- function(X, A) {
    comps <- connected.components(A)
    if (comps$n <= 1L) return(matrix(integer(), ncol = 2))
    D <- pairwise.dist(X)
    candidates <- list()
    cursor <- 0L
    for (a in seq_len(comps$n - 1L)) {
        for (b in (a + 1L):comps$n) {
            ia <- which(comps$id == a)
            ib <- which(comps$id == b)
            sub <- D[ia, ib, drop = FALSE]
            pos <- arrayInd(which.min(sub), dim(sub))
            cursor <- cursor + 1L
            candidates[[cursor]] <- data.frame(
                component.from = a,
                component.to = b,
                from = ia[pos[1]],
                to = ib[pos[2]],
                weight = sub[pos[1], pos[2]]
            )
        }
    }
    C <- do.call(rbind, candidates)
    C <- C[order(C$weight, C$component.from, C$component.to), , drop = FALSE]
    parent <- seq_len(comps$n)
    find <- function(x) {
        while (parent[[x]] != x) {
            parent[[x]] <<- parent[[parent[[x]]]]
            x <- parent[[x]]
        }
        x
    }
    keep <- list()
    for (r in seq_len(nrow(C))) {
        ra <- find(C$component.from[[r]])
        rb <- find(C$component.to[[r]])
        if (ra != rb) {
            parent[[rb]] <- ra
            keep[[length(keep) + 1L]] <- as.integer(c(C$from[[r]], C$to[[r]]))
            if (length(keep) == comps$n - 1L) break
        }
    }
    do.call(rbind, keep)
}

fermat.distance <- function(X, k.base = 8, alpha = 2) {
    D <- pairwise.dist(X)
    A <- support.sknn(D, k.base)
    W <- matrix(Inf, nrow(D), nrow(D))
    diag(W) <- 0
    W[A] <- D[A]^alpha
    for (mid in seq_len(nrow(W))) {
        W <- pmin(W, outer(W[, mid], W[mid, ], "+"))
    }
    W
}

write.roots <- function(n, k, path) {
    theta <- 2 * pi * seq(0, n - 1L) / n
    X <- cbind(cos(theta), sin(theta))
    D <- pairwise.dist(X)
    panels <- list(
        PHATE = support.phate(D, k),
        sKNN = support.sknn(D, k),
        mKNN = support.mknn(D, k),
        iKNN = support.iknn(D, k)
    )
    pdf(path, width = 9.5, height = 2.6)
    par(mfrow = c(1, 4), mar = c(0.8, 0.8, 2.0, 0.8))
    for (nm in names(panels)) plot.graph(X, panels[[nm]], nm)
    dev.off()
}

write.noisy.circle <- function(n, k, seed, path, diagnostics.path = NULL) {
    set.seed(seed)
    theta <- sort(runif(n, 0, 2 * pi))
    radius <- 1 + 0.08 * rnorm(n)
    X <- cbind(radius * cos(theta), radius * sin(theta))
    Df <- fermat.distance(X, k.base = min(8, n - 1L), alpha = 2)
    panels <- list(
        PHATE = support.phate(Df, k),
        sKNN = support.sknn(Df, k),
        mKNN = support.mknn(Df, k),
        iKNN = support.iknn(Df, k)
    )
    pdf(path, width = 9.5, height = 2.6)
    par(mfrow = c(1, 4), mar = c(0.8, 0.8, 2.0, 0.8))
    for (nm in names(panels)) plot.graph(X, panels[[nm]], nm)
    dev.off()

    if (!is.null(diagnostics.path)) {
        phate.extra <- panels$PHATE & !panels$sKNN
        sknn.missing <- panels$sKNN & !panels$PHATE
        bridges <- component.mst.edges(X, panels$sKNN)
        pdf(diagnostics.path, width = 8.0, height = 2.7)
        par(mfrow = c(1, 3), mar = c(0.8, 0.8, 2.3, 0.8))
        plot.graph(X, panels$sKNN, "PHATE edges beyond sKNN", edge.col = "#cbd5e1")
        plot.edges(X, phate.extra, edge.col = "#dc2626", edge.lwd = 1.0)
        plot.graph(X, panels$sKNN, "sKNN absent from PHATE", edge.col = "#cbd5e1")
        plot.edges(X, sknn.missing, edge.col = "#7c3aed", edge.lwd = 1.0)
        plot.graph(X, panels$sKNN, "component-MST repair", edge.col = "#9ca3af")
        if (length(bridges)) {
            for (r in seq_len(nrow(bridges))) {
                lines(X[bridges[r, ], 1], X[bridges[r, ], 2], col = "#2563eb", lwd = 1.3)
            }
        }
        dev.off()
    }
}

write.construction.rules <- function(path) {
    set.seed(4)
    theta <- sort(c(runif(45, 0, 1.25 * pi), runif(25, 1.25 * pi, 2 * pi)))
    X <- cbind((1 + 0.035 * rnorm(length(theta))) * cos(theta),
               (1 + 0.035 * rnorm(length(theta))) * sin(theta))
    D <- pairwise.dist(X)
    eps <- median(apply(D + diag(Inf, nrow(D)), 1, function(z) sort(z)[5]))
    panels <- list(
        "fixed radius" = support.fixed(D, eps),
        "adaptive max radius" = support.adaptive(D, 5, radius.rule = "max"),
        "adaptive min radius" = support.adaptive(D, 5, radius.rule = "min"),
        "sKNN" = support.sknn(D, 5),
        "mKNN" = support.mknn(D, 5),
        "PHATE support" = support.phate(D, 5)
    )
    pdf(path, width = 8.2, height = 5.2)
    par(mfrow = c(2, 3), mar = c(0.7, 0.7, 2.0, 0.7))
    for (nm in names(panels)) plot.graph(X, panels[[nm]], nm)
    dev.off()
}

write.local.pruning <- function(path) {
    t <- seq(0.12, 2.82, length.out = 13)
    X <- cbind(cos(t), sin(t))
    A <- matrix(FALSE, nrow(X), nrow(X))
    for (i in seq_len(nrow(X) - 1L)) A[i, i + 1L] <- A[i + 1L, i] <- TRUE
    A[1, nrow(X)] <- A[nrow(X), 1] <- TRUE
    local <- seq_len(nrow(X))
    pdf(path, width = 8.5, height = 2.8)
    par(mfrow = c(1, 3), mar = c(0.8, 0.8, 2.2, 0.8))
    plot.graph(X, A, "candidate graph", edge.col = "#9ca3af")
    lines(X[c(1, nrow(X)), 1], X[c(1, nrow(X)), 2], col = "#dc2626", lwd = 2.0)
    plot(X[, 1], X[, 2], type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
         main = "local alternative path")
    points(X[local, 1], X[local, 2], pch = 21, bg = "#dbeafe", col = "#2563eb", cex = 1.1)
    for (i in seq_len(nrow(X) - 1L)) lines(X[c(i, i + 1L), 1], X[c(i, i + 1L), 2],
                                           col = "#16a34a", lwd = 1.8)
    lines(X[c(1, nrow(X)), 1], X[c(1, nrow(X)), 2], col = "#dc2626", lwd = 2.0, lty = 2)
    box(col = "#d1d5db")
    Apruned <- A
    Apruned[1, nrow(X)] <- Apruned[nrow(X), 1] <- FALSE
    plot.graph(X, Apruned, "after pruning", edge.col = "#16a34a")
    dev.off()
}

write.component.mst <- function(path) {
    set.seed(11)
    centers <- rbind(c(-1.5, 0.6), c(1.5, 0.7), c(0, -1.2))
    X <- do.call(rbind, lapply(seq_len(3), function(i) {
        cbind(rnorm(8, centers[i, 1], 0.18), rnorm(8, centers[i, 2], 0.15))
    }))
    D <- pairwise.dist(X)
    A <- support.fixed(D, 0.43)
    comps <- connected.components(A)
    bridges <- component.mst.edges(X, A)
    comp.cols <- c("#2563eb", "#16a34a", "#dc2626")
    pdf(path, width = 8.5, height = 2.8)
    par(mfrow = c(1, 3), mar = c(0.8, 0.8, 2.2, 0.8))
    for (panel in c("disconnected graph", "component candidates", "MST-repaired graph")) {
        plot(X[, 1], X[, 2], type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "", main = panel)
        plot.edges(X, A, edge.col = "#9ca3af")
        if (panel != "disconnected graph") {
            for (a in 1:2) for (b in (a + 1):3) {
                ia <- which(comps$id == a)
                ib <- which(comps$id == b)
                sub <- D[ia, ib, drop = FALSE]
                pos <- arrayInd(which.min(sub), dim(sub))
                lines(X[c(ia[pos[1]], ib[pos[2]]), 1], X[c(ia[pos[1]], ib[pos[2]]), 2],
                      col = "#64748b", lty = 2, lwd = 1.0)
            }
        }
        if (panel == "MST-repaired graph") {
            for (r in seq_len(nrow(bridges))) {
                lines(X[bridges[r, ], 1], X[bridges[r, ], 2], col = "#111827", lwd = 2.0)
            }
        }
        points(X[, 1], X[, 2], pch = 19, col = comp.cols[comps$id], cex = 0.75)
        box(col = "#d1d5db")
    }
    dev.off()
}

write.failure.examples <- function(path) {
    t <- seq(0, 2 * pi, length.out = 90)
    fig8 <- cbind(sin(t), sin(t) * cos(t))
    cross <- rbind(cbind(seq(-1, 1, length.out = 35), 0),
                   cbind(0, seq(-1, 1, length.out = 35)))
    near <- rbind(cbind(seq(-1, 1, length.out = 35), 0.08),
                  cbind(0, seq(-1, 1, length.out = 35)))
    tree <- rbind(cbind(seq(0, 1, length.out = 28), 0),
                  cbind(seq(1, 1.8, length.out = 24), seq(0, 0.9, length.out = 24)),
                  cbind(seq(1, 1.8, length.out = 24), seq(0, -0.9, length.out = 24)))
    theta <- c(seq(0, 1.3 * pi, length.out = 65), seq(1.3 * pi, 2 * pi, length.out = 25))
    nonuniform <- cbind((1 + 0.05 * sin(5 * theta)) * cos(theta), (1 + 0.05 * cos(3 * theta)) * sin(theta))
    u <- seq(-0.9, 0.9, length.out = 18)
    grid <- expand.grid(u = u, v = u)
    quad <- as.matrix(grid[grid$u^2 + grid$v^2 <= 0.9^2, ])
    examples <- list(
        "figure eight" = fig8,
        "crossing segments" = cross,
        "near crossing" = near,
        "branching tree" = tree,
        "nonuniform circle" = nonuniform,
        "quadratic parameter disk" = quad
    )
    pdf(path, width = 8.4, height = 5.4)
    par(mfrow = c(2, 3), mar = c(0.8, 0.8, 2.0, 0.8))
    for (nm in names(examples)) {
        X <- examples[[nm]]
        D <- pairwise.dist(X)
        A <- support.sknn(D, min(5, nrow(X) - 1L))
        plot.graph(X, A, nm, edge.col = "#9ca3af")
    }
    dev.off()
}

write.quadform.examples <- function(path) {
    u <- seq(-1, 1, length.out = 28)
    grid <- expand.grid(u = u, v = u)
    grid <- grid[grid$u^2 + grid$v^2 <= 1, ]
    make.z <- function(kind) if (kind == "paraboloid") grid$u^2 + grid$v^2 else grid$u^2 - grid$v^2
    pdf(path, width = 7.6, height = 4.0)
    par(mfrow = c(1, 2), mar = c(1.0, 1.0, 2.3, 1.0))
    for (kind in c("paraboloid", "saddle")) {
        z <- make.z(kind)
        plot(grid$u, grid$v, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
             main = if (kind == "paraboloid") expression(y == x[0]^2 + x[1]^2) else expression(y == x[0]^2 - x[1]^2))
        col <- grDevices::hcl.colors(100, "Viridis")[as.integer(cut(z, 100))]
        D <- pairwise.dist(as.matrix(grid))
        A <- support.sknn(D, 5)
        plot.edges(as.matrix(grid), A, edge.col = "#cbd5e1", edge.lwd = 0.35)
        points(grid$u, grid$v, pch = 19, cex = 0.55, col = col)
        box(col = "#d1d5db")
    }
    dev.off()
}

paths <- c(
    file.path(figure.dir, "roots_n3_k2.pdf"),
    file.path(figure.dir, "roots_n4_k2.pdf"),
    file.path(figure.dir, "roots_n5_k2.pdf"),
    file.path(figure.dir, "roots_n100_k10.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n50_k8.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n50_k8_diagnostics.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n100_k8.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n100_k8_diagnostics.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n150_k8.pdf"),
    file.path(figure.dir, "noisy_circle_fermat_n150_k8_diagnostics.pdf"),
    file.path(figure.dir, "construction_rules_nonuniform_circle.pdf"),
    file.path(figure.dir, "local_geodesic_pruning_diagram.pdf"),
    file.path(figure.dir, "component_mst_repair_diagram.pdf"),
    file.path(figure.dir, "failure_mode_examples.pdf"),
    file.path(figure.dir, "quadform_surface_examples.pdf")
)

write.roots(3, 2, paths[[1]])
write.roots(4, 2, paths[[2]])
write.roots(5, 2, paths[[3]])
write.roots(100, 10, paths[[4]])
write.noisy.circle(50, 8, 10, paths[[5]], paths[[6]])
write.noisy.circle(100, 8, 11, paths[[7]], paths[[8]])
write.noisy.circle(150, 8, 12, paths[[9]], paths[[10]])
write.construction.rules(paths[[11]])
write.local.pruning(paths[[12]])
write.component.mst(paths[[13]])
write.failure.examples(paths[[14]])
write.quadform.examples(paths[[15]])

log.path <- file.path(report.dir, "figure_generation.log")
writeLines(c("Wrote figures:", normalizePath(paths)), con = log.path)
message("Wrote ", length(paths), " figures")
