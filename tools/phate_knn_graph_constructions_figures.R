#!/usr/bin/env Rscript

## Figures for the PHATE/sKNN/mKNN/ikNN graph-construction note.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required.")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required.")
}

out.dir <- file.path("dev", "phate-knn-graph-constructions", "figures")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

decay <- 40
threshold <- 1e-4

pairwise.dist <- function(X) {
  as.matrix(stats::dist(X))
}

nn.sets <- function(D, k) {
  n <- nrow(D)
  lapply(seq_len(n), function(i) {
    ord <- order(D[i, ], seq_len(n))
    ord <- ord[ord != i]
    ord[seq_len(min(k, n - 1L))]
  })
}

closed.nn.sets <- function(D, k) {
  nn <- nn.sets(D, k)
  Map(c, seq_along(nn), nn)
}

edge.df.from.matrix <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(data.frame(i = integer(), j = integer()))
  }
  data.frame(i = pmin(edges[, 1], edges[, 2]), j = pmax(edges[, 1], edges[, 2]))
}

edges.mknn <- function(D, k) {
  n <- nrow(D)
  nn <- nn.sets(D, k)
  out <- list()
  m <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (j %in% nn[[i]] && i %in% nn[[j]]) {
        m <- m + 1L
        out[[m]] <- c(i, j)
      }
    }
  }
  if (m == 0L) matrix(integer(), ncol = 2L) else do.call(rbind, out)
}

edges.sknn <- function(D, k) {
  n <- nrow(D)
  nn <- nn.sets(D, k)
  out <- list()
  m <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (j %in% nn[[i]] || i %in% nn[[j]]) {
        m <- m + 1L
        out[[m]] <- c(i, j)
      }
    }
  }
  if (m == 0L) matrix(integer(), ncol = 2L) else do.call(rbind, out)
}

edges.iknn <- function(D, k) {
  n <- nrow(D)
  nn <- closed.nn.sets(D, k)
  out <- list()
  m <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (length(intersect(nn[[i]], nn[[j]])) > 0L) {
        m <- m + 1L
        out[[m]] <- c(i, j)
      }
    }
  }
  if (m == 0L) matrix(integer(), ncol = 2L) else do.call(rbind, out)
}

phate.kernel <- function(D, k, decay = 40, threshold = 1e-4) {
  n <- nrow(D)
  sigma <- numeric(n)
  for (i in seq_len(n)) {
    row <- D[i, ]
    row[i] <- Inf
    sigma[i] <- sort(row, partial = k)[k]
  }
  A <- matrix(0, n, n)
  diag(A) <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      val <- exp(-((D[i, j] / sigma[i])^decay))
      if (is.finite(val) && val >= threshold) A[i, j] <- val
    }
  }
  K <- 0.5 * (A + t(A))
  list(A = A, K = K, sigma = sigma)
}

edges.phate <- function(D, k, decay = 40, threshold = 1e-4) {
  K <- phate.kernel(D, k, decay, threshold)$K
  idx <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) matrix(integer(), ncol = 2L) else idx
}

graph.constructions <- function(X, k, D = pairwise.dist(X)) {
  list(
    PHATE = edges.phate(D, k, decay = decay, threshold = threshold),
    sKNN = edges.sknn(D, k),
    mKNN = edges.mknn(D, k),
    ikNN = edges.iknn(D, k)
  )
}

roots.unity <- function(n) {
  theta <- 2 * pi * (seq_len(n) - 1L) / n
  data.frame(
    id = seq_len(n),
    x = cos(theta),
    y = sin(theta),
    theta = theta
  )
}

edge.plot.data <- function(points, edges) {
  if (nrow(edges) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
  }
  data.frame(
    x = points$x[edges[, 1]],
    y = points$y[edges[, 1]],
    xend = points$x[edges[, 2]],
    yend = points$y[edges[, 2]]
  )
}

plot.graph <- function(points, edges, title, subtitle = NULL) {
  ep <- edge.plot.data(points, edges)
  ggplot2::ggplot(points, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(
      data = ep,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.35,
      alpha = 0.55,
      color = "#1f4e79"
    ) +
    ggplot2::geom_point(size = 1.7, color = "#b2182b") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8)
    )
}

edge.keys <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(character(0))
  }
  edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  sort(paste(edges[, 1], edges[, 2], sep = "-"))
}

edges.from.keys <- function(keys) {
  if (length(keys) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  vals <- do.call(rbind, strsplit(keys, "-", fixed = TRUE))
  matrix(as.integer(vals), ncol = 2L)
}

component.mst.edges <- function(D, edges) {
  n <- nrow(D)
  if (nrow(edges) == 0L) {
    comp <- seq_len(n)
    n.comp <- n
  } else {
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    if (igraph::vcount(g) < n) {
      g <- igraph::add_vertices(g, n - igraph::vcount(g))
    }
    comp <- igraph::components(g)$membership
    n.comp <- max(comp)
  }
  if (n.comp <= 1L) {
    return(matrix(integer(), ncol = 2L))
  }

  candidates <- data.frame()
  for (a in seq_len(n.comp - 1L)) {
    ia <- which(comp == a)
    for (b in (a + 1L):n.comp) {
      ib <- which(comp == b)
      grid <- expand.grid(i = ia, j = ib)
      grid$d <- D[cbind(grid$i, grid$j)]
      grid <- grid[order(grid$d, grid$i, grid$j), ]
      candidates <- rbind(candidates, grid[1L, c("i", "j", "d")])
    }
  }
  cg <- igraph::graph_from_data_frame(
    data.frame(from = comp[candidates$i], to = comp[candidates$j], weight = candidates$d),
    directed = FALSE,
    vertices = data.frame(name = seq_len(n.comp))
  )
  mst <- igraph::mst(cg, weights = igraph::E(cg)$weight)
  if (igraph::ecount(mst) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  ends <- igraph::ends(mst, igraph::E(mst), names = TRUE)
  out <- matrix(integer(), nrow = nrow(ends), ncol = 2L)
  for (r in seq_len(nrow(ends))) {
    a <- as.integer(ends[r, 1])
    b <- as.integer(ends[r, 2])
    idx <- which((comp[candidates$i] == a & comp[candidates$j] == b) |
                   (comp[candidates$i] == b & comp[candidates$j] == a))
    best <- candidates[idx[which.min(candidates$d[idx])], ]
    out[r, ] <- as.integer(c(best$i, best$j))
  }
  cbind(pmin(out[, 1], out[, 2]), pmax(out[, 1], out[, 2]))
}

plot.edge.diff <- function(points,
                           base.edges,
                           highlight.edges,
                           title,
                           subtitle = NULL,
                           highlight.color = "#d73027",
                           base.color = "#c9c9c9",
                           highlight.label = "highlight") {
  base.ep <- edge.plot.data(points, base.edges)
  high.ep <- edge.plot.data(points, highlight.edges)
  ggplot2::ggplot(points, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(
      data = base.ep,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.45,
      color = base.color
    ) +
    ggplot2::geom_segment(
      data = high.ep,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.5,
      alpha = 0.9,
      color = highlight.color
    ) +
    ggplot2::geom_point(size = 1.5, color = "#333333") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8)
    )
}

save.root.example <- function(n, k) {
  pts <- roots.unity(n)
  X <- as.matrix(pts[, c("x", "y")])
  graphs <- graph.constructions(X, k)
  panels <- lapply(names(graphs), function(name) {
    plot.graph(pts, graphs[[name]], name, sprintf("%d edges", nrow(graphs[[name]])))
  })
  p <- (panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) +
    patchwork::plot_annotation(
      title = sprintf("%d roots of unity, k = %d", n, k),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  path <- file.path(out.dir, sprintf("roots_n%d_k%d.pdf", n, k))
  ggplot2::ggsave(path, p, width = 9.4, height = 2.7, device = grDevices::cairo_pdf)
  invisible(path)
}

save.root.large <- function(n = 100L, k = 10L) {
  pts <- roots.unity(n)
  X <- as.matrix(pts[, c("x", "y")])
  graphs <- graph.constructions(X, k)
  panels <- lapply(names(graphs), function(name) {
    plot.graph(pts, graphs[[name]], name, sprintf("%d edges", nrow(graphs[[name]])))
  })
  p <- (panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) +
    patchwork::plot_annotation(
      title = sprintf("%d roots of unity, k = %d", n, k),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  path <- file.path(out.dir, sprintf("roots_n%d_k%d.pdf", n, k))
  ggplot2::ggsave(path, p, width = 9.4, height = 2.7, device = grDevices::cairo_pdf)
  invisible(path)
}

fermat.distance <- function(X, base.k = 12L, alpha = 2) {
  D <- pairwise.dist(X)
  n <- nrow(D)
  nn <- nn.sets(D, min(base.k, n - 1L))
  edges <- do.call(rbind, lapply(seq_len(n), function(i) {
    cbind(i, nn[[i]])
  }))
  edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  edges <- unique(edges)
  weights <- D[cbind(edges[, 1], edges[, 2])]^alpha
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  if (igraph::vcount(g) < n) {
    g <- igraph::add_vertices(g, n - igraph::vcount(g))
  }
  igraph::E(g)$weight <- weights
  as.matrix(igraph::distances(g, weights = igraph::E(g)$weight))
}

noisy.circle <- function(n, seed) {
  set.seed(seed)
  theta <- sort(stats::runif(n, 0, 2 * pi))
  radius <- 1 + stats::rnorm(n, sd = 0.075)
  data.frame(
    id = seq_len(n),
    x = radius * cos(theta) + stats::rnorm(n, sd = 0.025),
    y = radius * sin(theta) + stats::rnorm(n, sd = 0.025)
  )
}

save.noisy.fermat <- function(n, k = 8L, base.k = 12L, alpha = 2) {
  pts <- noisy.circle(n, seed = 71000L + n)
  X <- as.matrix(pts[, c("x", "y")])
  Df <- fermat.distance(X, base.k = base.k, alpha = alpha)
  graphs <- graph.constructions(X, k, D = Df)
  panels <- lapply(names(graphs), function(name) {
    plot.graph(pts, graphs[[name]], name, sprintf("%d edges", nrow(graphs[[name]])))
  })
  p <- (panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) +
    patchwork::plot_annotation(
      title = sprintf("Noisy circle, n = %d, Fermat distance alpha = %g, graph k = %d", n, alpha, k),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  path <- file.path(out.dir, sprintf("noisy_circle_fermat_n%d_k%d.pdf", n, k))
  ggplot2::ggsave(path, p, width = 9.4, height = 2.7, device = grDevices::cairo_pdf)

  phate.keys <- edge.keys(graphs$PHATE)
  sknn.keys <- edge.keys(graphs$sKNN)
  phate.only <- edges.from.keys(setdiff(phate.keys, sknn.keys))
  sknn.only <- edges.from.keys(setdiff(sknn.keys, phate.keys))
  mst.edges <- component.mst.edges(Df, graphs$sKNN)
  sknn.mst <- unique(rbind(graphs$sKNN, mst.edges))
  sknn.mst <- sknn.mst[order(sknn.mst[, 1], sknn.mst[, 2]), , drop = FALSE]
  mst.subtitle <- if (nrow(mst.edges) == 0L) {
    "sKNN connected; no MST bridges"
  } else {
    sprintf("%d MST bridge edges highlighted", nrow(mst.edges))
  }
  diff.panels <- list(
    plot.edge.diff(
      pts,
      graphs$sKNN,
      phate.only,
      "PHATE-only",
      sprintf("%d edges in PHATE not in sKNN", nrow(phate.only)),
      highlight.color = "#d73027"
    ),
    plot.edge.diff(
      pts,
      graphs$PHATE,
      sknn.only,
      "sKNN-only",
      sprintf("%d edges in sKNN not in PHATE", nrow(sknn.only)),
      highlight.color = "#1a9850"
    ),
    plot.edge.diff(
      pts,
      sknn.mst,
      mst.edges,
      "sKNN + component-MST",
      mst.subtitle,
      highlight.color = "#fdae61"
    )
  )
  dp <- (diff.panels[[1]] | diff.panels[[2]] | diff.panels[[3]]) +
    patchwork::plot_annotation(
      title = sprintf("PHATE versus sKNN diagnostics, n = %d, graph k = %d", n, k),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  diag.path <- file.path(out.dir, sprintf("noisy_circle_fermat_n%d_k%d_diagnostics.pdf", n, k))
  ggplot2::ggsave(diag.path, dp, width = 9.4, height = 2.7, device = grDevices::cairo_pdf)
  invisible(c(path, diag.path))
}

paths <- c(
  save.root.example(3L, 2L),
  save.root.example(4L, 2L),
  save.root.example(5L, 2L),
  save.root.large(100L, 10L),
  save.noisy.fermat(50L),
  save.noisy.fermat(100L),
  save.noisy.fermat(150L)
)

message("Wrote figures:")
message(paste(normalizePath(paths), collapse = "\n"))
