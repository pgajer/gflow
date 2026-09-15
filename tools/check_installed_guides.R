#!/usr/bin/env Rscript
# Validate the installed distribution, independently of checkout-generated help.
library(gflow)
for (topic in c("gflow", "gflow-package", "gflow-migration", "create.basin.complex", "lcor")) {
    if (length(help(topic, package = "gflow")) != 1L) stop("Missing installed help: ", topic)
}
guides <- c("function-guide", "example-graphs-and-fields",
            "basin_complex_workflow_vignette", "noisy_circle_core_workflow_vignette")
index <- vignette(package = "gflow")$results
stopifnot(all(guides %in% index[, "Item"]))
for (guide in guides) {
    stopifnot(file.exists(system.file("doc", paste0(guide, ".html"), package = "gflow")))
}
local({
    scratch <- tempfile("gflow-guides-")
    dir.create(scratch)
    old <- getwd()
    on.exit({setwd(old); unlink(scratch, recursive = TRUE)}, add = TRUE)
    setwd(scratch)
    grDevices::pdf("guide-example.pdf")
    on.exit(grDevices::dev.off(), add = TRUE)
    source(system.file("doc", "function-guide.R", package = "gflow"),
           local = new.env(parent = globalenv()))
})
cat("Installed help, all four guides, and the introductory workflow passed.\n")
