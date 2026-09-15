#!/usr/bin/env Rscript
# Reproduce the README/guide illustration from explicit, unit-length path data.
library(gflow)
adjacency <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), c(4L, 6L), c(5L, 7L), 6L)
lengths <- lapply(adjacency, function(v) rep(1, length(v)))
field <- c(0, 3, 3, 1, 2, 4, 0)
bc <- create.basin.complex(adjacency, lengths, field,
                           method = "superlevel_merge_tree", direction = "max")
stopifnot(bc$status == "ok")
dir.create("vignettes/figures", recursive = TRUE, showWarnings = FALSE)
local({
    png("vignettes/figures/basin-showcase.png", width = 1400, height = 650, res = 150)
    on.exit(dev.off())
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), bg = "#fbfcfa", fg = "#20343c", col.axis = "#20343c", col.lab = "#20343c")
    plot(seq_along(field), field, type = "b", pch = 21, bg = "#16828d", col = "#20343c",
         lwd = 2, cex = 1.4, xlab = "Vertex on a unit-length path", ylab = "Scalar field",
         ylim = c(-0.3, 4.8), main = "Two peaks, one connected plateau")
    text(seq_along(field), field, labels = seq_along(field), pos = 3, cex = 0.8)
    abline(h = 1, lty = 3, col = "#687c80")
    plot(get.basin.merge.tree(bc), direction = "max", type = "tree",
         label = "extremum.vertex", show.mass = FALSE, show.support = FALSE,
         main.tree = "The peaks join at field height 1", field.label = "Scalar field")
})
cat("Wrote vignettes/figures/basin-showcase.png\n")
