##plot3D.graph <- function(x,
                         z = NULL;
                         layout = "kk";
                         conn.points = TRUE;
                         use.spheres = TRUE;
                         graph.alpha = 0.7;
                         z.point.size = 0.5;
                         z.color = NULL;
                         z.alpha = 1;
                         edge.color = "gray70";
                         edge.width = 1;
                         base.plane = TRUE;
                         base.vertex.size = 0.5;
                         z.scale = 1;
                         show.axes = TRUE;
                         vertical.lines = TRUE;
                         vertical.line.style = "dashed";
                         dash.length = 0.05;
                         gap.length = 0.05;
                         vertical.line.color = "darkgray";
                         vertical.line.alpha = 0.5;
                         vertical.line.width = 0.5;
                         basins.df = NULL;
                         evertex.sphere.radius = 0.2;
                         evertex.min.color = "blue";
                         evertex.max.color = "red";
                         evertex.cex = 1;
                         evertex.adj = c(0.5, 0.5);
                         evertex.label.offset = 0.3;
##                         ...) {

x <- plot.res
                                        # rgl gating
if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("This function requires the optional package 'rgl'. ",
         "Install with install.packages('rgl').", call. = FALSE)
}

                                        # Headless/CI-safe: render to the null device; harmless on desktops
old <- options(rgl.useNULL = TRUE)
on.exit(options(old), add = TRUE)

                                        # ---- extract/construct graph and 2D layout ----
if (inherits(x, "graph.3d")) {
    g <- x$graph
    layout.2d <- x$layout
    if (is.null(z)) z <- x$z
    vertex.color <- if (!is.null(x$vertex.color)) x$vertex.color else rep("gray", igraph::vcount(g))

} else if (inherits(x, "ggraph")) {
    graph.adj.list     <- x$adj.list
    graph.edge.lengths <- x$weight.list

    graph.obj <- convert.adjacency.to.edge.matrix(graph.adj.list, graph.edge.lengths)
    g <- igraph::graph_from_edgelist(graph.obj$edge.matrix, directed = FALSE)
    if (!is.null(graph.edge.lengths)) igraph::E(g)$weight <- graph.obj$weights

    if (is.character(layout)) {
        layout_func <- switch(layout,
                              "nicely"    = igraph::layout_nicely,
                              "randomly"  = igraph::layout_randomly,
                              "in_circle" = igraph::layout_in_circle,
                              "on_sphere" = igraph::layout_on_sphere,
                              "fr"        = igraph::layout_with_fr,
                              "kk"        = igraph::layout_with_kk,
                              "as_tree"   = igraph::layout_as_tree,
                              "lgl"       = igraph::layout_with_lgl,
                              "graphopt"  = igraph::layout_with_graphopt,
                              "drl"       = igraph::layout_with_drl,
                              "dh"        = igraph::layout_with_dh,
                              "mds"       = igraph::layout_with_mds,
                              stop("Invalid layout specified")
                              )
        layout.2d <- layout_func(g, dim = 2)
    } else if (is.matrix(layout) && ncol(layout) == 2) {
        layout.2d <- layout
    } else {
        stop("'layout' must be a known layout name or an n x 2 numeric matrix")
    }

    vertex.color <- rep("gray", igraph::vcount(g))

} else if (is.list(x) && !is.null(x$graph) && !is.null(x$layout)) {
    g <- x$graph
    layout.2d <- x$layout
    vertex.color <- if (!is.null(x$vertex.color)) x$vertex.color else rep("gray", igraph::vcount(g))

} else {
    stop("`x` must be a graph.3d object, a ggraph object, or a list with components 'graph' and 'layout'")
}

                                        # ---- validate z and layout ----
if (is.null(z)) {
    stop("`z` values must be provided (argument `z` or inside a graph.3d object)")
}
if (!is.numeric(z)) stop("`z` must be a numeric vector", call. = FALSE)
layout.2d <- as.matrix(layout.2d)
if (!is.numeric(layout.2d)) stop("`layout` must be numeric coordinates", call. = FALSE)
if (ncol(layout.2d) != 2L) stop("`layout` must have exactly 2 columns", call. = FALSE)
if (length(z) != nrow(layout.2d)) {
    stop(sprintf("Length of `z` (%d) must match number of vertices (%d).",
                 length(z), nrow(layout.2d)), call. = FALSE)
}
if (!is.numeric(graph.alpha) || graph.alpha < 0 || graph.alpha > 1) {
    stop("`graph.alpha` must be in [0, 1]", call. = FALSE)
}
if (!is.numeric(z.alpha) || z.alpha < 0 || z.alpha > 1) {
    stop("`z.alpha` must be in [0, 1]", call. = FALSE)
}

                                        # ---- scale & colors ----
z <- as.numeric(z) * z.scale
layout.3d <- cbind(layout.2d, z)

if (is.null(z.color)) {
    if (diff(range(z)) > 0) {
        z.norm <- (z - min(z)) / diff(range(z))
        pal <- grDevices::hcl.colors(100, palette = "Spectral")
        idx <- pmin(100L, pmax(1L, 1L + floor(z.norm * 99)))
        z.color <- pal[idx]
    } else {
        z.color <- rep("red", length(z))
    }
} else if (length(z.color) == 1L) {
    z.color <- rep(z.color, length(z))
} else if (length(z.color) != length(z)) {
    stop("`z.color` must be a single color or have the same length as `z`", call. = FALSE)
}

                                        # ---- open/clear device ----
rgl::open3d()
on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
rgl::clear3d()

                                        # ---- base plane (z=0) ----
if (isTRUE(base.plane)) {
    base.layout <- cbind(layout.2d, 0)

    if (isTRUE(use.spheres)) {
        rgl::spheres3d(
                 base.layout,
                 radius = base.vertex.size * 0.02,
                 col = vertex.color,
                 alpha = graph.alpha
             )
    } else {
        rgl::points3d(
                 base.layout,
                 size = base.vertex.size * 10,
                 col = vertex.color,
                 alpha = graph.alpha
             )
    }

    edges <- igraph::as_edgelist(g, names = FALSE)
    if (NROW(edges) > 0) {
        edge.coords <- matrix(NA_real_, nrow = NROW(edges) * 2L, ncol = 3L)
        for (i in seq_len(NROW(edges))) {
            v1 <- edges[i, 1]; v2 <- edges[i, 2]
            edge.coords[(i - 1L) * 2L + 1L, ] <- base.layout[v1, ]
            edge.coords[(i - 1L) * 2L + 2L, ] <- base.layout[v2, ]
        }
        rgl::segments3d(edge.coords, col = "gray50", alpha = graph.alpha, lwd = edge.width)
    }
}

## ---- z points ----
if (isTRUE(use.spheres)) {
    rgl::spheres3d(
             layout.3d,
             radius = z.point.size * 0.02,
             col = z.color,
             alpha = z.alpha
         )
} else {
    rgl::points3d(
             layout.3d,
             size = z.point.size * 10,
             col = z.color,
             alpha = z.alpha
         )
}

rgl::spheres3d(layout.3d, radius = z.point.size, col = z.color, alpha = z.alpha)

                                        # ---- 3D edges between z points ----
if (isTRUE(conn.points)) {
    edges <- igraph::as_edgelist(g, names = FALSE)
    if (NROW(edges) > 0) {
        for (i in seq_len(NROW(edges))) {
            v1 <- edges[i, 1]; v2 <- edges[i, 2]
            rgl::lines3d(rbind(layout.3d[v1, ], layout.3d[v2, ]),
                         col = edge.color, lwd = edge.width)
        }
    }
}

                                        # ---- vertical lines from base to z points ----
if (isTRUE(vertical.lines) && isTRUE(base.plane)) {
    if (identical(vertical.line.style, "dashed")) {
        for (i in seq_len(nrow(layout.2d))) {
            z_height <- z[i]
            if (z_height <= 0) next
            seg_len <- dash.length + gap.length
            n_segments <- floor(z_height / seg_len)
            if (n_segments < 1L) n_segments <- 1L
            for (j in 0:(n_segments - 1L)) {
                start_z <- j * seg_len
                end_z   <- min(start_z + dash.length, z_height)
                rgl::lines3d(
                         rbind(c(layout.2d[i, ], start_z),
                               c(layout.2d[i, ], end_z)),
                         col   = vertical.line.color,
                         alpha = vertical.line.alpha,
                         lwd   = vertical.line.width
                     )
            }
        }
    } else {
        base.layout <- cbind(layout.2d, 0)
        for (i in seq_len(nrow(layout.2d))) {
            rgl::lines3d(
                     rbind(base.layout[i, ], layout.3d[i, ]),
                     col   = vertical.line.color,
                     alpha = vertical.line.alpha,
                     lwd   = vertical.line.width
                 )
        }
    }
}

                                        # ---- extrema (basins.df) ----
if (!is.null(basins.df)) {
    if (!all(c("evertex", "is_max") %in% names(basins.df))) {
        stop("`basins.df` must contain columns 'evertex' and 'is_max'", call. = FALSE)
    }
    labels <- if (!is.null(basins.df$label)) basins.df$label else paste0("E", seq_len(nrow(basins.df)))
    for (i in seq_len(nrow(basins.df))) {
        v <- basins.df$evertex[i]
        if (!is.finite(v) || v < 1L || v > nrow(layout.3d)) next
        colv <- if (isTRUE(basins.df$is_max[i] == 1)) evertex.max.color else evertex.min.color
        rgl::spheres3d(layout.3d[v, 1], layout.3d[v, 2], layout.3d[v, 3],
                       radius = evertex.sphere.radius * 0.02, col = colv, alpha = 0.8)
        rgl::text3d(layout.3d[v, 1], layout.3d[v, 2], layout.3d[v, 3] + evertex.label.offset,
                    texts = labels[i], col = colv, cex = evertex.cex, adj = evertex.adj)
    }
}

                                        # ---- axes ----
if (isTRUE(show.axes)) {
    rgl::axes3d()
    rgl::title3d(xlab = "X", ylab = "Y", zlab = "Z")
}

invisible(NULL)
##}
