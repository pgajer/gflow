#!/usr/bin/env Rscript

# Purpose: Demonstrate the macOS Tahoe rgl WebGL workaround in a simple 3D scene.
# Usage:   Rscript tests/manual/rgl_webgl_workaround_demo_driver.R [output_html]
# Output:  A standalone HTML file with an interactive WebGL 3D plot.

args <- commandArgs(trailingOnly = TRUE)
output_html <- if (length(args) >= 1L) args[[1L]] else "rgl_webgl_workaround_demo.html"

if (!nzchar(system.file(package = "rgl"))) {
  stop("Package 'rgl' is required. Install it first (install.packages('rgl')).", call. = FALSE)
}
if (!nzchar(system.file(package = "htmlwidgets"))) {
  stop("Package 'htmlwidgets' is required. Install it first (install.packages('htmlwidgets')).", call. = FALSE)
}

# Set workaround options BEFORE loading rgl.
options(
  rgl.useNULL = TRUE,
  rgl.printRglwidget = TRUE
)

suppressPackageStartupMessages(library(rgl))

set.seed(42)

# Build a small smooth surface with overlay points.
x <- seq(-pi, pi, length.out = 60L)
y <- seq(-pi, pi, length.out = 60L)
z <- outer(x, y, function(xi, yi) cos(xi) + sin(yi))

invisible(open3d())
bg3d(color = "white")
material3d(specular = "black")

persp3d(
  x = x, y = y, z = z,
  col = "steelblue", alpha = 0.85,
  xlab = "x", ylab = "y", zlab = "z",
  main = "rgl WebGL workaround demo"
)

pts_x <- runif(120L, min = -pi, max = pi)
pts_y <- runif(120L, min = -pi, max = pi)
pts_z <- cos(pts_x) + sin(pts_y)
points3d(pts_x, pts_y, pts_z, color = "firebrick", size = 6)
text3d(0, 0, 2.2, texts = "WebGL fallback active", color = "black", cex = 1.2)

widget <- rglwidget(width = 900, height = 700)
htmlwidgets::saveWidget(widget, file = output_html, selfcontained = TRUE)

cat("Saved interactive WebGL demo to:", normalizePath(output_html), "\n")
cat("Workaround options: rgl.useNULL=TRUE, rgl.printRglwidget=TRUE\n")
